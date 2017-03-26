#include <psi4/psi4-dec.h>
#include <psi4/libparallel/parallel.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libmints/matrix.h>
#include <psi4/libmints/molecule.h>
#include <psi4/libmints/wavefunction.h>
#include <psi4/libpsio/psio.hpp>
#include <psi4/libdpd/dpd.h>
#include <psi4/libtrans/integraltransform.h>
#include <psi4/libiwl/iwl.hpp>
#include <psi4/libciomr/libciomr.h>

// This allows us to be lazy in getting the spaces in DPD calls
#define ID(x) ints.DPD_ID(x)

namespace psi{ namespace fcidump{

extern "C"
int read_options(std::string name, Options &options)
{
    if (name == "FCIDUMP"|| options.read_globals()) {

        options.add_str("DUMPFILENAME", "FCIDUMP");

    }

    return true;
}

extern "C"
SharedWavefunction fcidump(SharedWavefunction wfn, Options& options)
{

    const std::string filenamefcidump = options.get_str("DUMPFILENAME");

    // Grab the global (default) PSIO object, for file I/O
    std::shared_ptr<PSIO> psio(_default_psio_lib_);
    const int nIrreps  = wfn->nirrep();

    // For now, we'll just transform for closed shells and generate all integrals.  For more elaborate use of the
    // LibTrans object, check out the plugin_mp2 example.
    std::vector<std::shared_ptr<MOSpace> > spaces;
    spaces.push_back(MOSpace::all);
    IntegralTransform ints(wfn, spaces, IntegralTransform::Restricted);
    ints.transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);
    dpd_set_default(ints.get_dpd_id());

    //Readin the MO OEI in moOei & print everything
    int nmo       = wfn->nmo();
    int *orbspi   = init_int_array(nIrreps);
    int *docc     = init_int_array(nIrreps);
    int *socc     = init_int_array(nIrreps);
    for ( int h = 0; h < nIrreps; ++h ){
        orbspi[h] = wfn->nmopi()[h];
        docc[h] = wfn->doccpi()[h];
        socc[h] = wfn->soccpi()[h];
    }

    int nTriMo = nmo * (nmo + 1) / 2;
    double *temp = new double[nTriMo];
    Matrix moOei("MO OEI", nIrreps, orbspi, orbspi);
    IWL::read_one(psio.get(), PSIF_OEI, PSIF_MO_OEI, temp, nTriMo, 0, 0, "outfile");
    moOei.set(temp);

    int numElectrons = 0;
    int numMS2 = 0;
    int targetIrrep = 0;
    for (int h=0; h<nIrreps; h++){
        numElectrons += 2*docc[h] + socc[h];
        numMS2 += socc[h];
        if ( socc[h] % 2 != 0 ){ targetIrrep = targetIrrep ^ h; }
    }
    std::string SymmLabel = Process::environment.molecule()->sym_label();
    int * symm_psi2molpro = new int[nIrreps];
    if ( SymmLabel.compare("c1")==0 ){
        symm_psi2molpro[0] = 1;
    }
    if ( ( SymmLabel.compare("ci")==0 ) || ( SymmLabel.compare("c2")==0 ) || ( SymmLabel.compare("cs")==0 ) ){
        symm_psi2molpro[0] = 1;
        symm_psi2molpro[1] = 2;
    }
    if ( ( SymmLabel.compare("d2")==0 ) ){
        symm_psi2molpro[0] = 1;
        symm_psi2molpro[1] = 4;
        symm_psi2molpro[2] = 3;
        symm_psi2molpro[3] = 2;
    }
    if ( ( SymmLabel.compare("c2v")==0 ) || ( SymmLabel.compare("c2h")==0 ) ){
        symm_psi2molpro[0] = 1;
        symm_psi2molpro[1] = 4;
        symm_psi2molpro[2] = 2;
        symm_psi2molpro[3] = 3;
    }
    if ( ( SymmLabel.compare("d2h")==0 ) ){
        symm_psi2molpro[0] = 1;
        symm_psi2molpro[1] = 4;
        symm_psi2molpro[2] = 6;
        symm_psi2molpro[3] = 7;
        symm_psi2molpro[4] = 8;
        symm_psi2molpro[5] = 5;
        symm_psi2molpro[6] = 3;
        symm_psi2molpro[7] = 2;
    }

    FILE * capturing;
    capturing = fopen( filenamefcidump.c_str(), "w" ); // "w" with fopen means truncate file
    fprintf( capturing, " &FCI NORB= %d,NELEC= %d,MS2= %d,\n", nmo, numElectrons, numMS2 );
    fprintf( capturing, "  ORBSYM=" );
    for (int h=0; h<nIrreps; h++){
        for (int orb=0; orb<orbspi[h]; orb++){
            fprintf( capturing, "%d,", symm_psi2molpro[h] );
        }
    }
    fprintf( capturing, "\n  ISYM=%d,\n /\n", symm_psi2molpro[targetIrrep] );
    delete [] symm_psi2molpro;

    // Now, loop over the DPD buffer, printing the integrals
    dpdbuf4 K;
    psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    // To only process the permutationally unique integrals, change the ID("[A,A]") to ID("[A>=A]+")
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[A,A]"), ID("[A,A]"), ID("[A>=A]+"), ID("[A>=A]+"), 0, "MO Ints (AA|AA)");
    for ( int h = 0; h < nIrreps; ++h ){
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for ( int pq = 0; pq < K.params->rowtot[h]; ++pq ){
            const int p = K.params->roworb[h][pq][0];
            const int q = K.params->roworb[h][pq][1];
            if ( p >= q ){
                for ( int rs = 0; rs < K.params->coltot[h]; ++rs ){
                    const int r = K.params->colorb[h][rs][0];
                    const int s = K.params->colorb[h][rs][1];
                    if ( ( r >= s ) && ( ( p > r ) || ( ( p == r ) && ( q >= s ) ) ) ){
                        fprintf( capturing, " % 23.16E %3d %3d %3d %3d\n", K.matrix[h][pq][rs], p+1, q+1, r+1, s+1 );
                    }
                }
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);
    psio->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    int jump = 0;
    for ( int h = 0; h < nIrreps; ++h ){
       for ( int p = 0; p < moOei.rowspi(h); ++p ){
          for ( int q = 0; q <= p; ++q ){
              fprintf( capturing, " % 23.16E %3d %3d %3d %3d\n", moOei[h][p][q], jump+p+1, jump+q+1, 0, 0 );
          }
       }
       jump += moOei.rowspi(h);
    }

    delete [] temp;

    fprintf( capturing, " % 23.16E %3d %3d %3d %3d", Process::environment.molecule()->nuclear_repulsion_energy(), 0, 0, 0, 0 );
    fclose( capturing );

    outfile->Printf( "Created the file %s", filenamefcidump.c_str() );

    return wfn;
}

}} // End Namespaces
