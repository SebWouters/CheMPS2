/*
   CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
   Copyright (C) 2013-2016 Sebastian Wouters

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along
   with this program; if not, write to the Free Software Foundation, Inc.,
   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#include <iostream>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <string.h>
#include <getopt.h>
#include <sys/stat.h>

#include "Initialize.h"
#include "DMRG.h"
#include "MPIchemps2.h"

using namespace std;

void fetch_ints( const string rawdata, int * result, const int num ){

   int pos  = 0;
   int pos2 = 0;
   for ( int no = 0; no < num; no++ ){
      pos2 = rawdata.find( ",", pos );
      if ( pos2 == string::npos ){ pos2 = rawdata.length(); }
      result[ no ] = atoi(rawdata.substr( pos, pos2-pos ).c_str());
      pos = pos2 + 1;
   }

}

void fetch_doubles( const string rawdata, double * result, const int num ){

   int pos  = 0;
   int pos2 = 0;
   for ( int no = 0; no < num; no++ ){
      pos2 = rawdata.find( ",", pos );
      if ( pos2 == string::npos ){ pos2 = rawdata.length(); }
      result[ no ] = atof(rawdata.substr( pos, pos2-pos ).c_str());
      pos = pos2 + 1;
   }

}

void print_help(){

cout << "\n"
"CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry\n"
"Copyright (C) 2013-2016 Sebastian Wouters\n"
"\n"
"Usage: chemps2 [OPTIONS]\n"
"\n"

/**************************************************
* The following is copied directly from manpage.1 *
**************************************************/

"   SYMMETRY\n"
"       Conventions for the symmetry group and irrep numbers (same as psi4):\n"
"\n"
"                        |  0    1    2    3    4    5    6    7\n"
"               ---------|-----------------------------------------\n"
"               0 : c1   |  A\n"
"               1 : ci   |  Ag   Au\n"
"               2 : c2   |  A    B\n"
"               3 : cs   |  Ap   App\n"
"               4 : d2   |  A    B1   B2   B3\n"
"               5 : c2v  |  A1   A2   B1   B2\n"
"               6 : c2h  |  Ag   Bg   Au   Bu\n"
"               7 : d2h  |  Ag   B1g  B2g  B3g  Au   B1u  B2u  B3u\n"
"\n"
"   ARGUMENTS\n"
"       -f, --fcidump=filename\n"
"              Set the fcidump filename. Note that orbital irreps in this file follow molpro convention!\n"
"\n"
"       -g, --group=int\n"
"              Set the psi4 symmetry group number [0-7] which corresponds to the fcidump file.\n"
"\n"
"       -m, --multiplicity=int\n"
"              Overwrite the spin multiplicity [2S+1] of the fcidump file.\n"
"\n"
"       -n, --nelectrons=int\n"
"              Overwrite the number of electrons of the fcidump file.\n"
"\n"
"       -i, --irrep=int\n"
"              Overwrite the target wavefunction irrep [0-7] of the fcidump file (psi4 convention).\n"
"\n"
"       -D, --sweep_d=int,int,int\n"
"              Set the bond dimensions for the successive sweep instructions (positive integers).\n"
"\n"
"       -E, --sweep_econv=flt,flt,flt\n"
"              Set the energy convergence to stop sweep instructions (positive floats).\n"
"\n"
"       -M, --sweep_maxit=int,int,int\n"
"              Set the maximum number of sweeps for the sweep instructions (positive integers).\n"
"\n"
"       -N, --sweep_noise=flt,flt,flt\n"
"              Set the noise prefactors for the successive sweep instructions (floats).\n"
"\n"
"       -e, --excitation=int\n"
"              Set which excitation should be calculated (positive integer). If not set, the ground state is calculated.\n"
"\n"
"       -o, --twodmfile=filename\n"
"              Set the filename to dump the 2-RDM. If not set, the 2-RDM is not dumped.\n"
"\n"
"       -c, --checkpoint\n"
"              Read and create MPS checkpoints.\n"
"\n"
"       -p, --print_corr\n"
"              Print correlation functions.\n"
"\n"
"       -t, --tmpfolder=path\n"
"              Overwrite the tmp folder for the renormalized operators (default /tmp).\n"
"\n"
"       -r, --reorder=int,int,int\n"
"              Specify an orbital reordering w.r.t. the fcidump file (counting starts at 0).\n"
"\n"
"       -h, --help\n"
"              Display this help.\n"
"\n"
"   EXAMPLE\n"
"        $ cd /tmp\n"
"        $ wget \'https://github.com/SebWouters/CheMPS2/raw/master/tests/matrixelements/H2O.631G.FCIDUMP\'\n"
"        $ ls -al H2O.631G.FCIDUMP\n"
"        $ chemps2 --fcidump=H2O.631G.FCIDUMP \\\n"
"                  --group=5 \\\n"
"                  --sweep_d=200,1000 \\\n"
"                  --sweep_econv=1e-8,1e-8 \\\n"
"                  --sweep_maxit=2,10 \\\n"
"                  --sweep_noise=0.05,0.0 \\\n"
"                  --twodmfile=2dm.out \\\n"
"                  --print_corr \\\n"
"                  --reorder=6,5,4,3,2,1,0,7,8,9,10,11,12\n"
" " << endl;


}

int main(int argc, char **argv){

   #ifdef CHEMPS2_MPI_COMPILATION
      CheMPS2::MPIchemps2::mpi_init();
      const bool output = ( CheMPS2::MPIchemps2::mpi_rank() == MPI_CHEMPS2_MASTER );
   #else
      const bool output = true;
   #endif

   /*****************************
   *  Reading in the arguments  *
   *****************************/

   string fcidump     = "";
   int group          = -1;
   int multiplicity   = -1;
   int nelectrons     = -1;
   int irrep          = -1;
   string sweep_d     = "";
   string sweep_econv = "";
   string sweep_maxit = "";
   string sweep_noise = "";
   int excitation     = 0; // If nothing is passed the ground state is calculated
   string twodmfile   = "";
   bool checkpoint    = false;
   bool print_corr    = false;
   string tmpfolder   = CheMPS2::defaultTMPpath;
   string reorder     = "";

   struct option long_options[] =
   {
      {"fcidump",      required_argument, 0, 'f'},
      {"group",        required_argument, 0, 'g'},
      {"multiplicity", required_argument, 0, 'm'},
      {"nelectrons",   required_argument, 0, 'n'},
      {"irrep",        required_argument, 0, 'i'},
      {"sweep_d",      required_argument, 0, 'D'},
      {"sweep_econv",  required_argument, 0, 'E'},
      {"sweep_maxit",  required_argument, 0, 'M'},
      {"sweep_noise",  required_argument, 0, 'N'},
      {"excitation",   required_argument, 0, 'e'},
      {"twodmfile",    required_argument, 0, 'o'},
      {"checkpoint",   no_argument,       0, 'c'},
      {"print_corr",   no_argument,       0, 'p'},
      {"tmpfolder",    required_argument, 0, 't'},
      {"reorder",      required_argument, 0, 'r'},
      {"help",         no_argument,       0, 'h'},
      {0, 0, 0, 0}
   };

   int option_index = 0;
   int c;
   while((c = getopt_long(argc, argv, "hf:g:m:n:i:D:E:M:N:e:o:cpt:r:", long_options, &option_index)) != -1){
      switch(c){
         case 'h':
         case '?':
            if ( output ){ print_help(); }
            return 0;
            break;
         case 'f':
            fcidump = optarg;
            {
               struct stat file_info;
               const bool file_exists = ( stat( fcidump.c_str(), &file_info ) == 0 );
               if ( file_exists == false ){
                  if ( output ){ cerr << "fcidump file " << fcidump << " does not exist!" << endl; }
                  return -1;
               }
            }
            break;
         case 'g':
            group = atoi(optarg);
            if (( group < 0 ) || ( group > 7 )){
               if ( output ){ cerr << "Invalid group number!" << endl; }
               return -1;
            }
            break;
         case 'm':
            multiplicity = atoi(optarg);
            if ( multiplicity < 1 ){
               if ( output ){ cerr << "Invalid multiplicity!" << endl; }
               return -1;
            }
            break;
         case 'n':
            nelectrons = atoi(optarg);
            if ( nelectrons < 2 ){
               if ( output ){ cerr << "Invalid number of electrons!" << endl; }
               return -1;
            }
            break;
         case 'i':
            irrep = atoi(optarg);
            if (( irrep < 0 ) || ( irrep > 7 )){
               if ( output ){ cerr << "Invalid irrep!" << endl; }
               return -1;
            }
            break;
         case 'D':
            sweep_d = optarg;
            break;
         case 'E':
            sweep_econv = optarg;
            break;
         case 'M':
            sweep_maxit = optarg;
            break;
         case 'N':
            sweep_noise = optarg;
            break;
         case 'e':
            excitation = atoi(optarg);
            if ( excitation < 1 ){
               if ( output ){ cerr << "Invalid excitation number!" << endl; }
               return -1;
            }
            break;
         case 'o':
            twodmfile = optarg;
            break;
         case 'c':
            checkpoint = true;
            break;
         case 'p':
            print_corr = true;
            break;
         case 't':
            tmpfolder = optarg;
            if ( tmpfolder.length()==0 ){
               if ( output ){ cerr << "Invalid tmp path!" << endl; }
               return -1;
            }
            {
               struct stat file_info;
               const bool file_exists = ( stat( tmpfolder.c_str(), &file_info ) == 0 );
               if ( file_exists == false ){
                  if ( output ){ cerr << "tmp folder " << tmpfolder << " does not exist!" << endl; }
                  return -1;
               }
            }
            break;
         case 'r':
            reorder = optarg;
            break;
      }
   }
   
   /*******************************************
   *  Checking argument consistency (part 1)  *
   *******************************************/
   
   if ( fcidump.length()==0 ){
      if ( output ){ cerr << "The fcidump file should be specified!" << endl; }
      return -1;
   }
   if ( group == -1 ){
      if ( output ){ cerr << "The group number should be specified!" << endl; }
      return -1;
   }
   
   /******************************************
   *  Fetching unset arguments from FCIDUMP  *
   ******************************************/
   
   int fcidump_norb  = -1;
   int fcidump_nelec = -1;
   int fcidump_two_s = -1;
   int fcidump_irrep = -1;
   {
      ifstream thefcidump( fcidump.c_str() );
      string line;
      int pos, pos2;
      getline( thefcidump, line ); // &FCI NORB= X,NELEC= Y,MS2= Z,
      pos = line.find( "FCI" );
      if ( pos == string::npos ){
         if ( output ){ cerr << "The file " << fcidump << " is not a fcidump file!" << endl; }
         return -1;
      }
      pos = line.find( "NORB"  ); pos = line.find( "=", pos ); pos2 = line.find( ",", pos );
      fcidump_norb = atoi( line.substr( pos+1, pos2-pos-1 ).c_str() );
      pos = line.find( "NELEC" ); pos = line.find( "=", pos ); pos2 = line.find( ",", pos );
      fcidump_nelec = atoi( line.substr( pos+1, pos2-pos-1 ).c_str() );
      pos = line.find( "MS2"   ); pos = line.find( "=", pos ); pos2 = line.find( ",", pos );
      fcidump_two_s = atoi( line.substr( pos+1, pos2-pos-1 ).c_str() );
      do { getline( thefcidump, line ); } while ( line.find( "ISYM" ) == string::npos );
      pos = line.find( "ISYM"  ); pos = line.find( "=", pos ); pos2 = line.find( ",", pos );
      const int molpro_wfn_irrep = atoi( line.substr( pos+1, pos2-pos-1 ).c_str() );
      thefcidump.close();
      
      CheMPS2::Irreps Symmhelper(group);
      const int nIrreps = Symmhelper.getNumberOfIrreps();
      int * psi2molpro = new int[ nIrreps ];
      Symmhelper.symm_psi2molpro( psi2molpro );
      for ( int irrep = 0; irrep < nIrreps; irrep++ ){
         if ( molpro_wfn_irrep == psi2molpro[ irrep ] ){ fcidump_irrep = irrep; }
      }
      if ( fcidump_irrep == -1 ){
         if ( output ){ cerr << "Could not find the molpro wavefunction symmetry (ISYM) in the fcidump file!" << endl; }
         return -1;
      }
      delete [] psi2molpro;
   }
   if ( multiplicity == -1 ){ multiplicity = fcidump_two_s + 1; }
   if ( nelectrons   == -1 ){   nelectrons = fcidump_nelec;     }
   if ( irrep        == -1 ){        irrep = fcidump_irrep;     }
   
   /*******************************************
   *  Checking argument consistency (part 2)  *
   *******************************************/
   
   if (( sweep_d.length() == 0 ) || ( sweep_econv.length() == 0 ) || ( sweep_maxit.length() == 0 ) || ( sweep_noise.length() == 0 )){
      if ( output ){ cerr << "All sweep instructions should be specified!" << endl; }
      return -1;
   }
   
   const int ni_d     = count(sweep_d.begin(),     sweep_d.end(),     ',') + 1;
   const int ni_econv = count(sweep_econv.begin(), sweep_econv.end(), ',') + 1;
   const int ni_maxit = count(sweep_maxit.begin(), sweep_maxit.end(), ',') + 1;
   const int ni_noise = count(sweep_noise.begin(), sweep_noise.end(), ',') + 1;
   const bool num_eq  = (( ni_d == ni_econv ) && ( ni_d == ni_maxit ) && ( ni_d == ni_noise ));
   
   if ( num_eq == false ){
      if ( output ){ cerr << "The number of instruction lines in sweep_* should be equal!" << endl; }
      return -1;
   }
   
   int    * value_d     = new int[ni_d];       fetch_ints( sweep_d,     value_d,     ni_d );
   double * value_econv = new double[ni_d]; fetch_doubles( sweep_econv, value_econv, ni_d );
   int    * value_maxit = new int[ni_d];       fetch_ints( sweep_maxit, value_maxit, ni_d );
   double * value_noise = new double[ni_d]; fetch_doubles( sweep_noise, value_noise, ni_d );
   
   int * val_reorder = NULL;
   int ni_reo = -1;
   if ( reorder.length() > 0 ){
      ni_reo = count( reorder.begin(), reorder.end(), ',') + 1;
      val_reorder = new int[ ni_reo ];
      fetch_ints( reorder, val_reorder, ni_reo );
      if ( fcidump_norb != ni_reo ){
         if ( output ){ cerr << "The orbital reordering should contain as many elements as there are orbitals in the fcidump!" << endl; }
         return -1;
      }
      int * doublecheck = new int[ ni_reo ];
      for ( int cnt = 0; cnt < ni_reo; cnt++ ){ doublecheck[ cnt ] = 0; }
      for ( int cnt = 0; cnt < ni_reo; cnt++ ){
         if (( val_reorder[ cnt ] >= 0 ) && ( val_reorder[ cnt ] < ni_reo )){ doublecheck[ val_reorder[ cnt ] ] += 1; }
      }
      bool is_ok = true;
      for ( int cnt = 0; cnt < ni_reo; cnt++ ){ if ( doublecheck[ cnt ] != 1 ){ is_ok = false; } }
      delete [] doublecheck;
      if ( is_ok == false ){
         if ( output ){ cerr << "The orbital reordering should have all orbitals from 0 to " << fcidump_norb - 1 << " exactly once!" << endl; }
         return -1;
      }
   }
   
   if ( output ){
      CheMPS2::Irreps Symmhelper(group);
      cout << "\nRunning chemps2 with the following options:\n" << endl;
      cout << "  --fcidump = "      << fcidump      << endl;
      cout << "  --group = "        << Symmhelper.getGroupName() << endl;
      cout << "  --multiplicity = " << multiplicity << endl;
      cout << "  --nelectrons = "   << nelectrons   << endl;
      cout << "  --irrep = "        << Symmhelper.getIrrepName(irrep) << endl;
      cout << "  --sweep_d     = [ "; for (int cnt=0; cnt<ni_d-1; cnt++){ cout << value_d[cnt]     << " ; "; } cout << value_d[ni_d-1]     << " ]" << endl;
      cout << "  --sweep_econv = [ "; for (int cnt=0; cnt<ni_d-1; cnt++){ cout << value_econv[cnt] << " ; "; } cout << value_econv[ni_d-1] << " ]" << endl;
      cout << "  --sweep_maxit = [ "; for (int cnt=0; cnt<ni_d-1; cnt++){ cout << value_maxit[cnt] << " ; "; } cout << value_maxit[ni_d-1] << " ]" << endl;
      cout << "  --sweep_noise = [ "; for (int cnt=0; cnt<ni_d-1; cnt++){ cout << value_noise[cnt] << " ; "; } cout << value_noise[ni_d-1] << " ]" << endl;
      if ( excitation > 0 ){           cout << "  --excitation = "   << excitation   << endl; }
      if ( twodmfile.length() > 0 ){   cout << "  --twodmfile = "    << twodmfile    << endl; }
      if ( checkpoint ){               cout << "  --checkpoint"      << endl; }
      if ( print_corr ){               cout << "  --print_corr"      << endl; }
      cout << "  --tmpfolder = "    << tmpfolder    << endl;
      if ( ni_reo > 0 ){ cout << "  --reorder = [ "; for (int cnt=0; cnt<ni_reo-1; cnt++){ cout << val_reorder[cnt] << " ; "; } cout << val_reorder[ni_reo-1] << " ]" << endl; }
      cout << " " << endl;
   }
   
   /********************************
   *  Running the DMRG calculation *
   ********************************/

   //Initialize a bunch of stuff
   CheMPS2::Initialize::Init();
   CheMPS2::Hamiltonian * Ham = new CheMPS2::Hamiltonian( fcidump, group );
   CheMPS2::Problem * Prob = new CheMPS2::Problem( Ham, multiplicity-1, nelectrons, irrep );
   if ( ni_reo > 0 ){
      Prob->setup_reorder_custom( val_reorder );
      delete [] val_reorder;
   }
   
   //The convergence scheme
   CheMPS2::ConvergenceScheme * OptScheme = new CheMPS2::ConvergenceScheme(ni_d);
   //OptScheme->setInstruction(instruction, DSU(2), Econvergence, maxSweeps, noisePrefactor);
   for ( int instruction = 0; instruction < ni_d; instruction++ ){
      OptScheme->setInstruction( instruction, value_d[ instruction ],
                                          value_econv[ instruction ],
                                          value_maxit[ instruction ],
                                          value_noise[ instruction ] );
   }
   
   delete [] value_d;
   delete [] value_econv;
   delete [] value_maxit;
   delete [] value_noise;
   
   //Run the DMRG calculations
   CheMPS2::DMRG * theDMRG = new CheMPS2::DMRG(Prob, OptScheme, checkpoint, tmpfolder);
   double Energy = 0.0;
   for (int state = 0; state <= excitation; state++){
      if (state > 0){ theDMRG->newExcitation( fabs( Energy ) ); }
      Energy = theDMRG->Solve();
      if ((state == 0) && (excitation >= 1)){ theDMRG->activateExcitations( excitation ); }
   }
   
   //Calculate the 2-RDM and correlation functions
   if (( twodmfile.length() > 0 ) || ( print_corr == true )){
      theDMRG->calc2DMandCorrelations();
      if (( twodmfile.length() > 0 ) && ( output )){ theDMRG->get2DM()->write2DMAfile( twodmfile ); }
      if (( print_corr == true ) && ( output )){ theDMRG->getCorrelations()->Print(); }
   }
   
   //Clean up DMRG
   if (CheMPS2::DMRG_storeRenormOptrOnDisk){ theDMRG->deleteStoredOperators(); }
   delete theDMRG;
   delete OptScheme;
   delete Prob;
   delete Ham;
   
   #ifdef CHEMPS2_MPI_COMPILATION
   CheMPS2::MPIchemps2::mpi_finalize();
   #endif
   
   return 0;

}


