/*
   CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
   Copyright (C) 2013-2018 Sebastian Wouters

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
#include <assert.h>
#include <sstream>

#include "Initialize.h"
#include "CASSCF.h"
#include "Molden.h"
#include "MPIchemps2.h"
#include "EdmistonRuedenberg.h"

using namespace std;

void fetch_ints( const string rawdata, int * result, const int num ){

   int pos  = 0;
   int pos2 = 0;
   for ( int no = 0; no < num; no++ ){
      pos2 = rawdata.find( ",", pos );
      if ( pos2 == string::npos ){ pos2 = rawdata.length(); }
      result[ no ] = atoi( rawdata.substr( pos, pos2-pos ).c_str() );
      pos = pos2 + 1;
   }

}

void fetch_doubles( const string rawdata, double * result, const int num ){

   int pos  = 0;
   int pos2 = 0;
   for ( int no = 0; no < num; no++ ){
      pos2 = rawdata.find( ",", pos );
      if ( pos2 == string::npos ){ pos2 = rawdata.length(); }
      result[ no ] = atof( rawdata.substr( pos, pos2-pos ).c_str() );
      pos = pos2 + 1;
   }

}

bool file_exists( const string filename, const string tag ){

   #ifdef CHEMPS2_MPI_COMPILATION
      const bool am_i_master = ( CheMPS2::MPIchemps2::mpi_rank() == MPI_CHEMPS2_MASTER );
   #else
      const bool am_i_master = true;
   #endif

   struct stat file_info;
   const bool on_disk = (( filename.length() > 0 ) && ( stat( filename.c_str(), &file_info ) == 0 ));
   if (( on_disk == false ) && ( am_i_master )){
      cerr << "Unable to retrieve file " << filename << "!" << endl;
      cerr << "Invalid option for " << tag << "!" << endl;
   }
   return on_disk;

}

bool find_integer( int * result, const string line, const string tag, const bool lower_bound, const int val_lower, const bool upper_bound, const int val_upper ){

   #ifdef CHEMPS2_MPI_COMPILATION
      const bool am_i_master = ( CheMPS2::MPIchemps2::mpi_rank() == MPI_CHEMPS2_MASTER );
   #else
      const bool am_i_master = true;
   #endif

   if ( line.find( tag ) != string::npos ){

      const int pos = line.find( "=" ) + 1;
      result[ 0 ] = atoi( line.substr( pos, line.length() - pos ).c_str() );

      const bool lower_ok = (( lower_bound == false ) || ( result[ 0 ] >= val_lower ));
      const bool upper_ok = (( upper_bound == false ) || ( result[ 0 ] <= val_upper ));

      if (( lower_ok == false ) || ( upper_ok == false )){
         if ( am_i_master ){
            cerr << line << endl;
            cerr << "Invalid option for " << tag << "!" << endl;
         }
         return false;
      }
   }

   return true;

}

bool find_double( double * result, const string line, const string tag, const bool lower_bound, const double val_lower ){

   #ifdef CHEMPS2_MPI_COMPILATION
      const bool am_i_master = ( CheMPS2::MPIchemps2::mpi_rank() == MPI_CHEMPS2_MASTER );
   #else
      const bool am_i_master = true;
   #endif

   if ( line.find( tag ) != string::npos ){

      const int pos = line.find( "=" ) + 1;
      result[ 0 ] = atof( line.substr( pos, line.length() - pos ).c_str() );

      const bool lower_ok = (( lower_bound == false ) || ( result[ 0 ] >= val_lower ));

      if ( lower_ok == false ){
         if ( am_i_master ){
            cerr << line << endl;
            cerr << "Invalid option for " << tag << "!" << endl;
         }
         return false;
      }
   }

   return true;

}

bool find_character( char * result, const string line, const string tag, char * options, const int num_options ){

   #ifdef CHEMPS2_MPI_COMPILATION
      const bool am_i_master = ( CheMPS2::MPIchemps2::mpi_rank() == MPI_CHEMPS2_MASTER );
   #else
      const bool am_i_master = true;
   #endif

   if ( line.find( tag ) != string::npos ){

      const int pos = line.find( "=" ) + 1;
      string temp = line.substr( pos, line.length() - pos );
      temp.erase( std::remove( temp.begin(), temp.end(), ' ' ), temp.end() );
      result[ 0 ] = temp.c_str()[ 0 ];

      bool encountered = false;
      for ( int cnt = 0; cnt < num_options; cnt++ ){
         if ( options[ cnt ] == result[ 0 ] ){ encountered = true; }
      }

      if ( encountered == false ){
         if ( am_i_master ){
            cerr << line << endl;
            cerr << "Invalid option for " << tag << "!" << endl;
         }
         return false;
      }
   }

   return true;

}

bool find_boolean( bool * result, const string line, const string tag ){

   #ifdef CHEMPS2_MPI_COMPILATION
      const bool am_i_master = ( CheMPS2::MPIchemps2::mpi_rank() == MPI_CHEMPS2_MASTER );
   #else
      const bool am_i_master = true;
   #endif

   if ( line.find( tag ) != string::npos ){

      const int pos       = line.find( "=" ) + 1;
      const int pos_true  = line.substr( pos, line.length() - pos ).find( "TRUE" );
      const int pos_false = line.substr( pos, line.length() - pos ).find( "FALSE" );
      result[ 0 ] = ( pos_true != string::npos );

      if (( pos_true == string::npos ) && ( pos_false == string::npos )){
         if ( am_i_master ){
            cerr << line << endl;
            cerr << "Invalid option for " << tag << "!" << endl;
         }
         return false;
      }
   }

   return true;

}

int clean_exit( const int return_code ){

   #ifdef CHEMPS2_MPI_COMPILATION
   CheMPS2::MPIchemps2::mpi_finalize();
   #endif

   return return_code;

}

bool print_molcas_reorder( int * dmrg2ham, const int L, const string filename, const bool read ){

   bool on_disk = false;

   if ( read ){
      struct stat file_info;
      on_disk = (( filename.length() > 0 ) && ( stat( filename.c_str(), &file_info ) == 0 ));
      if ( on_disk ){
         ifstream input( filename.c_str() );
         string line;
         getline( input, line );
         const int num = count( line.begin(), line.end(), ',' ) + 1;
         assert( num == L );
         fetch_ints( line, dmrg2ham, L );
         input.close();
         cout << "Read orbital reordering = [ ";
         for ( int orb = 0; orb < L - 1; orb++ ){ cout << dmrg2ham[ orb ] << ", "; }
         cout << dmrg2ham[ L - 1 ] << " ]." << endl;
      }
   } else { // write
      FILE * capturing;
      capturing = fopen( filename.c_str(), "w" ); // "w" with fopen means truncate file
      for ( int orb = 0; orb < L - 1; orb++ ){
         fprintf( capturing, "%d, ", dmrg2ham[ orb ] );
      }
      fprintf( capturing, "%d \n", dmrg2ham[ L - 1 ] );
      fclose( capturing );
      cout << "Orbital reordering written to " << filename << "." << endl;
      on_disk = true;
   }

   return on_disk;

}

void print_help(){

cout << "\n"
"CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry\n"
"Copyright (C) 2013-2018 Sebastian Wouters\n"
"\n"
"Usage: chemps2 [OPTIONS]\n"
"\n"

/**************************************************
* The following is copied directly from chemps2.1 *
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
"       -f, --file=inputfile\n"
"              Specify the input file.\n"
"\n"
"       -v, --version\n"
"              Print the version of chemps2.\n"
"\n"
"       -h, --help\n"
"              Display this help.\n"
"\n"
"   INPUT FILE\n"
"       FCIDUMP = /path/to/fcidump\n"
"              Note that orbital irreps in the FCIDUMP file follow molpro convention!\n"
"\n"
"       GROUP = int\n"
"              Set the psi4 symmetry group number [0-7] which corresponds to the FCIDUMP file.\n"
"\n"
"       MULTIPLICITY = int\n"
"              Overwrite the spin multiplicity [2S+1] of the FCIDUMP file.\n"
"\n"
"       NELECTRONS = int\n"
"              Overwrite the number of electrons of the FCIDUMP file.\n"
"\n"
"       IRREP = int\n"
"              Overwrite the target wavefunction irrep [0-7] of the FCIDUMP file (psi4 convention).\n"
"\n"
"       EXCITATION = int\n"
"              Set which excitation should be calculated. If zero, the ground state is calculated (default 0).\n"
"\n"
"       SWEEP_STATES = int, int, int\n"
"              Set the number of reduced renormalized basis states for the successive sweep instructions (positive integers).\n"
"\n"
"       SWEEP_ENERGY_CONV = flt, flt, flt\n"
"              Set the energy convergence to stop the successive sweep instructions (positive floats).\n"
"\n"
"       SWEEP_MAX_SWEEPS = int, int, int\n"
"              Set the maximum number of sweeps for the successive sweep instructions (positive integers).\n"
"\n"
"       SWEEP_NOISE_PREFAC = flt, flt, flt\n"
"              Set the noise prefactors for the successive sweep instructions (floats).\n"
"\n"
"       SWEEP_DVDSON_RTOL = flt, flt, flt\n"
"              Set the residual norm tolerance for the Davidson algorithm for the successive sweep instructions (positive floats).\n"
"\n"
"       NOCC = int, int, int, int\n"
"              Set the number of occupied (external core) orbitals per irrep (psi4 irrep ordering).\n"
"\n"
"       NACT = int, int, int, int\n"
"              Set the number of active orbitals per irrep (psi4 irrep ordering).\n"
"\n"
"       NVIR = int, int, int, int\n"
"              Set the number of virtual (secondary) orbitals per irrep (psi4 irrep ordering).\n"
"\n"
"       MOLCAS_2RDM = /path/to/2rdm/output\n"
"              When all orbitals are active orbitals, write out the 2-RDM in HDF5 format when specified (default unspecified).\n"
"\n"
"       MOLCAS_3RDM = /path/to/3rdm/output\n"
"              When all orbitals are active orbitals, write out the 3-RDM in HDF5 format when specified (default unspecified).\n"
"\n"
"       MOLCAS_F4RDM = /path/to/f4rdm/output\n"
"              When all orbitals are active orbitals, write out the 4-RDM contracted with the Fock operator in HDF5 format when specified (default unspecified).\n"
"\n"
"       MOLCAS_FOCK = /path/to/fock/input\n"
"              When all orbitals are active orbitals, read in this file containing the Fock operator (default unspecified).\n"
"\n"
"       MOLCAS_FIEDLER = bool\n"
"              When all orbitals are active orbitals, switch on orbital reordering based on the Fiedler vector of the exchange matrix (TRUE or FALSE; default FALSE).\n"
"\n"
"       MOLCAS_ORDER = int, int, int, int\n"
"              When all orbitals are active orbitals, provide a custom orbital reordering (default unspecified). When specified, this option takes precedence over MOLCAS_FIEDLER.\n"
"\n"
"       MOLCAS_OCC = int, int, int, int\n"
"              When all orbitals are active orbitals, set initial guess to an ROHF determinant (default unspecified). The occupancy integers should be 0, 1 or 2 and the orbital ordering convention is FCIDUMP.\n"
"\n"
"       MOLCAS_MPS = bool\n"
"              When all orbitals are active orbitals, switch on the creation of MPS checkpoints (TRUE or FALSE; default FALSE).\n"
"\n"
"       MOLCAS_STATE_AVG = bool\n"
"              Switch on writing to disk of N-RDMs of intermediate roots (TRUE or FALSE; default FALSE).\n"
"\n"
"       SCF_STATE_AVG = bool\n"
"              Switch on state-averaging (TRUE or FALSE; default FALSE).\n"
"\n"
"       SCF_DIIS_THR = flt\n"
"              Switch on DIIS when the update norm is smaller than the given threshold (default 0.0).\n"
"\n"
"       SCF_GRAD_THR = flt\n"
"              Gradient norm threshold for convergence of the DMRG-SCF orbital rotations (default 1e-6).\n"
"\n"
"       SCF_MAX_ITER = int\n"
"              Specify the maximum number of DMRG-SCF iterations (default 100).\n"
"\n"
"       SCF_ACTIVE_SPACE = char\n"
"              Rotate the active space orbitals: no additional rotations (I), natural orbitals (N), localized and ordered orbitals (L), or ordered orbitals only (F) (default I).\n"
"\n"
"       SCF_MOLDEN = /path/to/molden\n"
"              Rotate the FCIDUMP orbitals to the DMRG-SCF occupied (external core), active, and virtual (secondary) orbitals.\n"
"\n"
"       CASPT2_CALC = bool\n"
"              Switch on the CASPT2 calculation (TRUE or FALSE; default FALSE).\n"
"\n"
"       CASPT2_ORBS = char\n"
"              Perform the DMRG calculation for the 4-RDM in the SCF_ACTIVE_SPACE orbitals (A) or in the pseudocanonical orbitals (P) (default A).\n"
"\n"
"       CASPT2_IPEA = flt\n"
"              Ionization potential - electron affinity shift (default 0.0).\n"
"\n"
"       CASPT2_IMAG = flt\n"
"              Imaginary level shift (default 0.0).\n"
"\n"
"       CASPT2_CHECKPT = bool\n"
"              Create checkpoints to continue the CASPT2 4-RDM calculation over multiple runs (TRUE or FALSE; default FALSE).\n"
"\n"
"       CASPT2_CUMUL = bool\n"
"              Use a cumulant approximation for the CASPT2 4-RDM and overwrite CASPT2_CHECKPT to FALSE (TRUE or FALSE; default FALSE).\n"
"\n"
"       PRINT_CORR = bool\n"
"              Print correlation functions (TRUE or FALSE; default FALSE).\n"
"\n"
"       TMP_FOLDER = /path/to/tmp/folder\n"
"              Overwrite the tmp folder for the renormalized operators. With MPI, separate folders per process can (but do not have to) be used (default /tmp).\n"
"\n"
"   EXAMPLE\n"
"       $ cd /tmp\n"
"       $ wget \'https://github.com/SebWouters/CheMPS2/raw/master/tests/matrixelements/N2.CCPVDZ.FCIDUMP\'\n"
"       $ ls -al N2.CCPVDZ.FCIDUMP\n"
"       $ wget \'https://github.com/SebWouters/CheMPS2/raw/master/tests/test14.input\'\n"
"       $ sed -i \"s/path\\/to/tmp/\" test14.input\n"
"       $ cat test14.input\n"
"       $ chemps2 --file=test14.input\n"
" " << endl;

}

int main( int argc, char ** argv ){

   #ifdef CHEMPS2_MPI_COMPILATION
      CheMPS2::MPIchemps2::mpi_init();
      const bool am_i_master = ( CheMPS2::MPIchemps2::mpi_rank() == MPI_CHEMPS2_MASTER );
   #else
      const bool am_i_master = true;
   #endif

   /************************
   *  Read in the options  *
   *************************/

   string inputfile = "";
   string fcidump   = "";

   int group        = -1;
   int multiplicity = -1;
   int nelectrons   = -1;
   int irrep        = -1;
   int excitation   = 0;

   string sweep_states = "";
   string sweep_econv  = "";
   string sweep_maxit  = "";
   string sweep_noise  = "";
   string sweep_rtol   = "";

   string nocc = "";
   string nact = "";
   string nvir = "";

   string molcas_2rdm      = "";
   string molcas_3rdm      = "";
   string molcas_f4rdm     = "";
   string molcas_fock      = "";
   bool   molcas_fiedler   = false;
   bool   molcas_mps       = false;
   bool   molcas_state_avg = false;
   string molcas_order     = "";
   string molcas_occ       = "";

   bool   scf_state_avg    = false;
   double scf_diis_thr     = 0.0;
   double scf_grad_thr     = 1e-6;
   int    scf_max_iter     = 100;
   char   scf_active_space = 'I';
   string scf_molden       = "";

   bool   caspt2_calc    = false;
   char   caspt2_orbs    = 'A';
   double caspt2_ipea    = 0.0;
   double caspt2_imag    = 0.0;
   bool   caspt2_checkpt = false;
   bool   caspt2_cumul   = false;

   bool   print_corr = false;
   string tmp_folder = "/tmp";

   struct option long_options[] =
   {
      {"file",    required_argument, 0, 'f'},
      {"version", no_argument,       0, 'v'},
      {"help",    no_argument,       0, 'h'},
      {0, 0, 0, 0}
   };

   int option_index = 0;
   int c;
   while (( c = getopt_long( argc, argv, "hvf:", long_options, &option_index )) != -1 ){
      switch( c ){
         case 'h':
         case '?':
            if ( am_i_master ){ print_help(); }
            return clean_exit( 0 );
            break;
         case 'v':
            if ( am_i_master ){ cout << "chemps2 version " << CHEMPS2_VERSION << endl; }
            return clean_exit( 0 );
            break;
         case 'f':
            inputfile = optarg;
            if ( file_exists( inputfile, "--file" ) == false ){ return clean_exit( -1 ); }
            break;
      }
   }

   if ( inputfile.length() == 0 ){
      if ( am_i_master ){ cerr << "The input file should be specified!" << endl; }
      return clean_exit( -1 );
   }

   ifstream input( inputfile.c_str() );
   string line;
   while ( input.eof() == false ){

      getline( input, line );

      if ( line.find( "FCIDUMP" ) != string::npos ){
         const int pos = line.find( "=" ) + 1;
         fcidump = line.substr( pos, line.length() - pos );
         fcidump.erase( remove( fcidump.begin(), fcidump.end(), ' ' ), fcidump.end() );
         if ( file_exists( fcidump, "FCIDUMP" ) == false ){ return clean_exit( -1 ); }
      }

      if ( line.find( "MOLCAS_2RDM" ) != string::npos ){
         const int pos = line.find( "=" ) + 1;
         molcas_2rdm = line.substr( pos, line.length() - pos );
         molcas_2rdm.erase( remove( molcas_2rdm.begin(), molcas_2rdm.end(), ' ' ), molcas_2rdm.end() );
      }

      if ( line.find( "MOLCAS_3RDM" ) != string::npos ){
         const int pos = line.find( "=" ) + 1;
         molcas_3rdm = line.substr( pos, line.length() - pos );
         molcas_3rdm.erase( remove( molcas_3rdm.begin(), molcas_3rdm.end(), ' ' ), molcas_3rdm.end() );
      }

      if ( line.find( "MOLCAS_F4RDM" ) != string::npos ){
         const int pos = line.find( "=" ) + 1;
         molcas_f4rdm = line.substr( pos, line.length() - pos );
         molcas_f4rdm.erase( remove( molcas_f4rdm.begin(), molcas_f4rdm.end(), ' ' ), molcas_f4rdm.end() );
      }

      if ( line.find( "MOLCAS_FOCK" ) != string::npos ){
         const int pos = line.find( "=" ) + 1;
         molcas_fock = line.substr( pos, line.length() - pos );
         molcas_fock.erase( remove( molcas_fock.begin(), molcas_fock.end(), ' ' ), molcas_fock.end() );
         if ( file_exists( molcas_fock, "MOLCAS_FOCK" ) == false ){ return clean_exit( -1 ); }
      }

      if ( line.find( "SCF_MOLDEN" ) != string::npos ){
         const int pos = line.find( "=" ) + 1;
         scf_molden = line.substr( pos, line.length() - pos );
         scf_molden.erase( remove( scf_molden.begin(), scf_molden.end(), ' ' ), scf_molden.end() );
         if ( file_exists( scf_molden, "SCF_MOLDEN" ) == false ){ return clean_exit( -1 ); }
      }

      if ( line.find( "TMP_FOLDER" ) != string::npos ){
         const int pos = line.find( "=" ) + 1;
         tmp_folder = line.substr( pos, line.length() - pos );
         tmp_folder.erase( remove( tmp_folder.begin(), tmp_folder.end(), ' ' ), tmp_folder.end() );
         if ( file_exists( tmp_folder, "TMP_FOLDER" ) == false ){ return clean_exit( -1 ); }
      }

      if ( find_integer( &group,        line, "GROUP",        true, 0, true,   7 ) == false ){ return clean_exit( -1 ); }
      if ( find_integer( &multiplicity, line, "MULTIPLICITY", true, 1, false, -1 ) == false ){ return clean_exit( -1 ); }
      if ( find_integer( &nelectrons,   line, "NELECTRONS",   true, 2, false, -1 ) == false ){ return clean_exit( -1 ); }
      if ( find_integer( &irrep,        line, "IRREP",        true, 0, true,   7 ) == false ){ return clean_exit( -1 ); }
      if ( find_integer( &excitation,   line, "EXCITATION",   true, 0, false, -1 ) == false ){ return clean_exit( -1 ); }
      if ( find_integer( &scf_max_iter, line, "SCF_MAX_ITER", true, 1, false, -1 ) == false ){ return clean_exit( -1 ); }

      if ( find_double( &scf_diis_thr, line, "SCF_DIIS_THR", true, 0.0 ) == false ){ return clean_exit( -1 ); }
      if ( find_double( &scf_grad_thr, line, "SCF_GRAD_THR", true, 0.0 ) == false ){ return clean_exit( -1 ); }
      if ( find_double( &caspt2_ipea,  line, "CASPT2_IPEA",  true, 0.0 ) == false ){ return clean_exit( -1 ); }
      if ( find_double( &caspt2_imag,  line, "CASPT2_IMAG",  true, 0.0 ) == false ){ return clean_exit( -1 ); }

      char options1[] = { 'I', 'N', 'L', 'F' };
      char options2[] = { 'A', 'P' };
      if ( find_character( &scf_active_space, line, "SCF_ACTIVE_SPACE", options1, 4 ) == false ){ return clean_exit( -1 ); }
      if ( find_character( &caspt2_orbs,      line, "CASPT2_ORBS",      options2, 2 ) == false ){ return clean_exit( -1 ); }

      if ( find_boolean( &molcas_fiedler,   line, "MOLCAS_FIEDLER"   ) == false ){ return clean_exit( -1 ); }
      if ( find_boolean( &molcas_mps,       line, "MOLCAS_MPS"       ) == false ){ return clean_exit( -1 ); }
      if ( find_boolean( &molcas_state_avg, line, "MOLCAS_STATE_AVG" ) == false ){ return clean_exit( -1 ); }
      if ( find_boolean( &scf_state_avg,    line, "SCF_STATE_AVG"    ) == false ){ return clean_exit( -1 ); }
      if ( find_boolean( &caspt2_calc,      line, "CASPT2_CALC"      ) == false ){ return clean_exit( -1 ); }
      if ( find_boolean( &caspt2_checkpt,   line, "CASPT2_CHECKPT"   ) == false ){ return clean_exit( -1 ); }
      if ( find_boolean( &caspt2_cumul,     line, "CASPT2_CUMUL"     ) == false ){ return clean_exit( -1 ); }
      if ( find_boolean( &print_corr,       line, "PRINT_CORR"       ) == false ){ return clean_exit( -1 ); }

      if ( line.find( "SWEEP_STATES" ) != string::npos ){
         const int pos = line.find( "=" ) + 1;
         sweep_states = line.substr( pos, line.length() - pos );
      }

      if ( line.find( "SWEEP_ENERGY_CONV" ) != string::npos ){
         const int pos = line.find( "=" ) + 1;
         sweep_econv = line.substr( pos, line.length() - pos );
      }

      if ( line.find( "SWEEP_MAX_SWEEPS" ) != string::npos ){
         const int pos = line.find( "=" ) + 1;
         sweep_maxit = line.substr( pos, line.length() - pos );
      }

      if ( line.find( "SWEEP_NOISE_PREFAC" ) != string::npos ){
         const int pos = line.find( "=" ) + 1;
         sweep_noise = line.substr( pos, line.length() - pos );
      }

      if ( line.find( "SWEEP_DVDSON_RTOL" ) != string::npos ){
         const int pos = line.find( "=" ) + 1;
         sweep_rtol = line.substr( pos, line.length() - pos );
      }

      if ( line.find( "NOCC" ) != string::npos ){
         const int pos = line.find( "=" ) + 1;
         nocc = line.substr( pos, line.length() - pos );
      }

      if ( line.find( "NACT" ) != string::npos ){
         const int pos = line.find( "=" ) + 1;
         nact = line.substr( pos, line.length() - pos );
      }

      if ( line.find( "NVIR" ) != string::npos ){
         const int pos = line.find( "=" ) + 1;
         nvir = line.substr( pos, line.length() - pos );
      }

      if ( line.find( "MOLCAS_ORDER" ) != string::npos ){
         const int pos = line.find( "=" ) + 1;
         molcas_order = line.substr( pos, line.length() - pos );
      }

      if ( line.find( "MOLCAS_OCC" ) != string::npos ){
         const int pos = line.find( "=" ) + 1;
         molcas_occ = line.substr( pos, line.length() - pos );
      }

      if ( line.find( "MOLCAS_REORDER" ) != string::npos ){
         if ( am_i_master ){ cerr << "MOLCAS_REORDER is deprecated. Please use MOLCAS_ORDER and/or MOLCAS_FIEDLER." << endl; }
         return clean_exit( -1 );
      }

   }
   input.close();

  /*******************************
   *  Check the target symmetry  *
   *******************************/

   if ( group == -1 ){
      if ( am_i_master ){ cerr << "GROUP is a mandatory option!" << endl; }
      return clean_exit( -1 );
   }
   CheMPS2::Irreps Symmhelper( group );
   const int num_irreps = Symmhelper.getNumberOfIrreps();

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
         if ( am_i_master ){ cerr << "The file " << fcidump << " is not a fcidump file!" << endl; }
         return clean_exit( -1 );
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

      int * psi2molpro = new int[ num_irreps ];
      Symmhelper.symm_psi2molpro( psi2molpro );
      for ( int cnt = 0; cnt < num_irreps; cnt++ ){
         if ( molpro_wfn_irrep == psi2molpro[ cnt ] ){ fcidump_irrep = cnt; }
      }
      if ( fcidump_irrep == -1 ){
         if ( am_i_master ){ cerr << "Could not find the molpro wavefunction symmetry (ISYM) in the fcidump file!" << endl; }
         return clean_exit( -1 );
      }
      delete [] psi2molpro;
   }
   if ( multiplicity == -1 ){ multiplicity = fcidump_two_s + 1; }
   if ( nelectrons   == -1 ){   nelectrons = fcidump_nelec;     }
   if ( irrep        == -1 ){        irrep = fcidump_irrep;     }

   /*********************************
   *  Check the sweep instructions  *
   **********************************/

   if (( sweep_states.length() == 0 ) || ( sweep_econv.length() == 0 ) || ( sweep_maxit.length() == 0 ) || ( sweep_noise.length() == 0 ) || ( sweep_rtol.length() == 0 )){
      if ( am_i_master ){ cerr << "SWEEP_* are mandatory options!" << endl; }
      return clean_exit( -1 );
   }
   const int ni_d     = count( sweep_states.begin(), sweep_states.end(), ',' ) + 1;
   const int ni_econv = count( sweep_econv.begin(),  sweep_econv.end(),  ',' ) + 1;
   const int ni_maxit = count( sweep_maxit.begin(),  sweep_maxit.end(),  ',' ) + 1;
   const int ni_noise = count( sweep_noise.begin(),  sweep_noise.end(),  ',' ) + 1;
   const int ni_rtol  = count( sweep_rtol.begin(),   sweep_rtol.end(),   ',' ) + 1;
   const bool num_eq  = (( ni_d == ni_econv ) && ( ni_d == ni_maxit ) && ( ni_d == ni_noise ) && ( ni_d == ni_rtol ));

   if ( num_eq == false ){
      if ( am_i_master ){ cerr << "The number of instructions in SWEEP_* should be equal!" << endl; }
      return clean_exit( -1 );
   }

   int    * value_states = new int   [ ni_d ];    fetch_ints( sweep_states, value_states, ni_d );
   double * value_econv  = new double[ ni_d ]; fetch_doubles( sweep_econv,  value_econv,  ni_d );
   int    * value_maxit  = new int   [ ni_d ];    fetch_ints( sweep_maxit,  value_maxit,  ni_d );
   double * value_noise  = new double[ ni_d ]; fetch_doubles( sweep_noise,  value_noise,  ni_d );
   double * value_rtol   = new double[ ni_d ]; fetch_doubles( sweep_rtol,   value_rtol,   ni_d );

   /*****************************************
   *  Check the active space specification  *
   ******************************************/

   if (( nocc.length() == 0 ) || ( nact.length() == 0 ) || ( nvir.length() == 0 )){
      if ( am_i_master ){ cerr << "NOCC, NACT, and NVIR are mandatory options!" << endl; }
      return clean_exit( -1 );
   }

   const int ni_occ  = count( nocc.begin(), nocc.end(), ',' ) + 1;
   const int ni_act  = count( nact.begin(), nact.end(), ',' ) + 1;
   const int ni_vir  = count( nvir.begin(), nvir.end(), ',' ) + 1;
   const bool cas_ok = (( ni_occ == ni_act ) && ( ni_occ == ni_vir ) && ( ni_occ == num_irreps ));

   if ( cas_ok == false ){
      if ( am_i_master ){ cerr << "There should be " << num_irreps << " numbers in NOCC, NACT, and NVIR!" << endl; }
      return clean_exit( -1 );
   }

   int * nocc_parsed = new int[ ni_occ ]; fetch_ints( nocc, nocc_parsed, ni_occ );
   int * nact_parsed = new int[ ni_act ]; fetch_ints( nact, nact_parsed, ni_act );
   int * nvir_parsed = new int[ ni_vir ]; fetch_ints( nvir, nvir_parsed, ni_vir );

   /**************************************
   *  Check consistency MOLCAS_ options  *
   ***************************************/

   bool full_active_space_calculation = true;
   for ( int cnt = 0; cnt < num_irreps; cnt++ ){
      if ( nocc_parsed[ cnt ] != 0 ){ full_active_space_calculation = false; }
      if ( nvir_parsed[ cnt ] != 0 ){ full_active_space_calculation = false; }
   }

   if ( ( molcas_2rdm.length() != 0 ) || ( molcas_3rdm.length() != 0 ) || ( molcas_f4rdm.length() != 0 ) || ( molcas_fock.length() != 0 ) || ( molcas_order.length() != 0 ) ){
      if ( full_active_space_calculation == false ){
         if ( am_i_master ){ cerr << "The options MOLCAS_* can only be specified for full active space calculations (when NOCC = NVIR = 0)!" << endl; }
         return clean_exit( -1 );
      }
   }

   if ( ( molcas_f4rdm.length() != 0 ) && ( molcas_fock.length() == 0 ) ){
      if ( am_i_master ){ cerr << "When MOLCAS_F4RDM should be written, MOLCAS_FOCK should be specified as well!" << endl; }
      return clean_exit( -1 );
   }

   /*********************************
   *  Parse reordering if required  *
   *********************************/

   int * dmrg2ham = NULL;
   if (( full_active_space_calculation == true ) && ( molcas_order.length() > 0 )){
      const int list_length = count( molcas_order.begin(), molcas_order.end(), ',' ) + 1;
      if ( list_length != fcidump_norb ){
         if ( am_i_master ){ cerr << "The number of integers specified in MOLCAS_ORDER should be equal to the number of orbitals in the FCIDUMP file!" << endl; }
         return clean_exit( -1 );
      }
      dmrg2ham = new int[ fcidump_norb ];
      fetch_ints( molcas_order, dmrg2ham, fcidump_norb );
   }

   /**********************************
   *  Parse occupancies if required  *
   **********************************/

   int * occupancies = NULL;
   if (( full_active_space_calculation == true ) && ( molcas_occ.length() > 0 )){
      const int list_length = count( molcas_occ.begin(), molcas_occ.end(), ',' ) + 1;
      if ( list_length != fcidump_norb ){
         if ( am_i_master ){ cerr << "The number of integers specified in MOLCAS_OCC should be equal to the number of orbitals in the FCIDUMP file!" << endl; }
         return clean_exit( -1 );
      }
      occupancies = new int[ fcidump_norb ];
      fetch_ints( molcas_occ, occupancies, fcidump_norb );
      int occ_n_tot = 0;
      int occ_2s_tot = 0;
      for ( int cnt = 0; cnt < fcidump_norb; cnt++ ){
         if (( occupancies[ cnt ] < 0 ) || ( occupancies[ cnt ] > 2 )){
            if ( am_i_master ){ cerr << "The integers specified in MOLCAS_OCC should be 0, 1 or 2!" << endl; }
            return clean_exit( -1 );
         }
         occ_n_tot += occupancies[ cnt ];
         if ( occupancies[ cnt ] == 1 ){ occ_2s_tot += 1; }
      }
      if ( occ_n_tot != nelectrons ){
         if ( am_i_master ){ cerr << "The sum of the integers specified in MOLCAS_OCC should be equal to the number of electrons specified in the FCIDUMP file!" << endl; }
         return clean_exit( -1 );
      }
      if ( ( occ_2s_tot + 1 ) != multiplicity ){
         if ( am_i_master ){ cerr << "The number of singly occupied orbitals in MOLCAS_OCC should be equal to the value 2S specified in the FCIDUMP file!" << endl; }
         return clean_exit( -1 );
      }
   }

   /**********************
   *  Print the options  *
   ***********************/

   if ( am_i_master ){
      cout << "\nRunning chemps2 version " << CHEMPS2_VERSION << " with the following options:\n" << endl;
      cout << "   FCIDUMP            = " << fcidump << endl;
      cout << "   GROUP              = " << Symmhelper.getGroupName() << endl;
      cout << "   MULTIPLICITY       = " << multiplicity << endl;
      cout << "   NELECTRONS         = " << nelectrons << endl;
      cout << "   IRREP              = " << Symmhelper.getIrrepName( irrep ) << endl;
      cout << "   EXCITATION         = " << excitation << endl;
      cout << "   SWEEP_STATES       = [ " << value_states[ 0 ]; for ( int cnt = 1; cnt < ni_d; cnt++ ){ cout << " ; " << value_states[ cnt ]; } cout << " ]" << endl;
      cout << "   SWEEP_ENERGY_CONV  = [ " << value_econv [ 0 ]; for ( int cnt = 1; cnt < ni_d; cnt++ ){ cout << " ; " << value_econv [ cnt ]; } cout << " ]" << endl;
      cout << "   SWEEP_MAX_SWEEPS   = [ " << value_maxit [ 0 ]; for ( int cnt = 1; cnt < ni_d; cnt++ ){ cout << " ; " << value_maxit [ cnt ]; } cout << " ]" << endl;
      cout << "   SWEEP_NOISE_PREFAC = [ " << value_noise [ 0 ]; for ( int cnt = 1; cnt < ni_d; cnt++ ){ cout << " ; " << value_noise [ cnt ]; } cout << " ]" << endl;
      cout << "   SWEEP_DVDSON_RTOL  = [ " << value_rtol  [ 0 ]; for ( int cnt = 1; cnt < ni_d; cnt++ ){ cout << " ; " << value_rtol  [ cnt ]; } cout << " ]" << endl;
      cout << "   NOCC               = [ " << nocc_parsed[ 0 ]; for ( int cnt = 1; cnt < num_irreps; cnt++ ){ cout << " ; " << nocc_parsed[ cnt ]; } cout << " ]" << endl;
      cout << "   NACT               = [ " << nact_parsed[ 0 ]; for ( int cnt = 1; cnt < num_irreps; cnt++ ){ cout << " ; " << nact_parsed[ cnt ]; } cout << " ]" << endl;
      cout << "   NVIR               = [ " << nvir_parsed[ 0 ]; for ( int cnt = 1; cnt < num_irreps; cnt++ ){ cout << " ; " << nvir_parsed[ cnt ]; } cout << " ]" << endl;
   if ( full_active_space_calculation ){
      cout << "   MOLCAS_2RDM        = " << molcas_2rdm << endl;
      cout << "   MOLCAS_3RDM        = " << molcas_3rdm << endl;
      cout << "   MOLCAS_F4RDM       = " << molcas_f4rdm << endl;
      cout << "   MOLCAS_FOCK        = " << molcas_fock << endl;
   if ( molcas_order.length() > 0 ){
      cout << "   MOLCAS_ORDER       = [ " << dmrg2ham[ 0 ]; for ( int cnt = 1; cnt < fcidump_norb; cnt++ ){ cout << " ; " << dmrg2ham[ cnt ]; } cout << " ]" << endl;
   } else {
      cout << "   MOLCAS_FIEDLER     = " << (( molcas_fiedler ) ? "TRUE" : "FALSE" ) << endl;
   }
   if ( molcas_occ.length() > 0 ){
      cout << "   MOLCAS_OCC (HAM)   = [ " << occupancies[ 0 ]; for ( int cnt = 1; cnt < fcidump_norb; cnt++ ){ cout << " ; " << occupancies[ cnt ]; } cout << " ]" << endl;
   }
      cout << "   MOLCAS_MPS         = " << (( molcas_mps ) ? "TRUE" : "FALSE" ) << endl;
      cout << "   MOLCAS_STATE_AVG   = " << (( molcas_state_avg ) ? "TRUE" : "FALSE" ) << endl;
   } else {
      cout << "   SCF_STATE_AVG      = " << (( scf_state_avg ) ? "TRUE" : "FALSE" ) << endl;
      cout << "   SCF_DIIS_THR       = " << scf_diis_thr << endl;
      cout << "   SCF_GRAD_THR       = " << scf_grad_thr << endl;
      cout << "   SCF_MAX_ITER       = " << scf_max_iter << endl;
      cout << "   SCF_ACTIVE_SPACE   = " << (( scf_active_space == 'I' ) ? "I : no additional rotations" :
                                            (( scf_active_space == 'N' ) ? "N : natural orbitals" :
                                            (( scf_active_space == 'L' ) ? "L : localized and ordered orbitals" : "F : ordered orbitals only" ))) << endl;
      cout << "   SCF_MOLDEN         = " << scf_molden << endl;
      cout << "   CASPT2_CALC        = " << (( caspt2_calc ) ? "TRUE" : "FALSE" ) << endl;
      cout << "   CASPT2_ORBS        = " << (( caspt2_orbs == 'A' ) ? "A : as specified in SCF_ACTIVE_SPACE" : "P : pseudocanonical orbitals" ) << endl;
      cout << "   CASPT2_IPEA        = " << caspt2_ipea << endl;
      cout << "   CASPT2_IMAG        = " << caspt2_imag << endl;
      cout << "   CASPT2_CHECKPT     = " << (( caspt2_checkpt ) ? "TRUE" : "FALSE" ) << endl;
      cout << "   CASPT2_CUMUL       = " << (( caspt2_cumul   ) ? "TRUE" : "FALSE" ) << endl;
   }
      cout << "   PRINT_CORR         = " << (( print_corr     ) ? "TRUE" : "FALSE" ) << endl;
      cout << "   TMP_FOLDER         = " << tmp_folder << endl;
      cout << " " << endl;
   }

   /********************************
   *  Running the DMRG calculation *
   ********************************/

   CheMPS2::Initialize::Init();
   CheMPS2::Hamiltonian * ham = new CheMPS2::Hamiltonian( fcidump, group );
   CheMPS2::ConvergenceScheme * opt_scheme = new CheMPS2::ConvergenceScheme( ni_d );
   for ( int count = 0; count < ni_d; count++ ){
      opt_scheme->set_instruction( count, value_states[ count ],
                                          value_econv [ count ],
                                          value_maxit [ count ],
                                          value_noise [ count ],
                                          value_rtol  [ count ] );
   }
   delete [] value_states;
   delete [] value_econv;
   delete [] value_maxit;
   delete [] value_noise;
   delete [] value_rtol;

   if ( full_active_space_calculation ){

      CheMPS2::Problem * prob = new CheMPS2::Problem( ham, multiplicity - 1, nelectrons, irrep );

      // Reorder the orbitals if desired
      if (( group == 7 ) && ( molcas_fiedler == false ) && ( molcas_order.length() == 0 )){ prob->SetupReorderD2h(); }
      if (( molcas_fiedler ) && ( molcas_order.length() == 0 )){
         dmrg2ham = new int[ ham->getL() ];
         if ( am_i_master ){
            const bool read_success = (( molcas_mps ) ? print_molcas_reorder( dmrg2ham, ham->getL(), "molcas_fiedler.txt", true ) : false );
            if ( read_success == false ){
               CheMPS2::EdmistonRuedenberg * fiedler = new CheMPS2::EdmistonRuedenberg( ham->getVmat(), group );
               fiedler->FiedlerGlobal( dmrg2ham );
               delete fiedler;
               if ( molcas_mps ){ print_molcas_reorder( dmrg2ham, ham->getL(), "molcas_fiedler.txt", false ); }
            }
         }
         #ifdef CHEMPS2_MPI_COMPILATION
         CheMPS2::MPIchemps2::broadcast_array_int( dmrg2ham, ham->getL(), MPI_CHEMPS2_MASTER );
         #endif
         prob->setup_reorder_custom( dmrg2ham );
         delete [] dmrg2ham;
      }
      if ( molcas_order.length() > 0 ){
         assert( fcidump_norb == ham->getL() );
         prob->setup_reorder_custom( dmrg2ham );
         delete [] dmrg2ham;
      }

      CheMPS2::DMRG * dmrgsolver = new CheMPS2::DMRG( prob, opt_scheme, molcas_mps, tmp_folder, occupancies );
      if ( molcas_occ.length() > 0 ){ delete [] occupancies; }

      // Solve for the correct root
      double DMRG_ENERGY;
      for ( int state = 0; state < ( excitation + 1 ); state++ ){
         if ( state > 0 ){ dmrgsolver->newExcitation( fabs( DMRG_ENERGY ) ); }
         DMRG_ENERGY = dmrgsolver->Solve();
         if (( state == 0 ) && ( excitation > 0 )){ dmrgsolver->activateExcitations( excitation ); }

         // Only if state specific or last state N-RDMs should be calculated
         if (( molcas_state_avg == true ) || ( state == excitation )){

            const bool calc_3rdm = (( molcas_3rdm.length() != 0 ) || ( molcas_f4rdm.length() != 0 ));
            const bool calc_2rdm = (( print_corr == true ) || ( molcas_2rdm.length() != 0 ));
            if (( calc_2rdm ) || ( calc_3rdm )){

               dmrgsolver->calc_rdms_and_correlations( calc_3rdm, false );
               std::stringstream result_filename;

               // 2-RDM
               if ( molcas_2rdm.length() != 0 ){
                  result_filename.str("");
                  result_filename << molcas_2rdm << ".r" << state;
                  dmrgsolver->get2DM()->save_HAM( result_filename.str() );
               }

               // 3-RDM
               if ( molcas_3rdm.length() != 0 ){
                  result_filename.str("");
                  result_filename << molcas_3rdm << ".r" << state;
                  dmrgsolver->get3DM()->save_HAM( result_filename.str() );
               }

               // F . 4-RDM
               if ( molcas_f4rdm.length() != 0 ){
                  const int LAS      = ham->getL();
                  const int LAS_pow6 = LAS * LAS * LAS * LAS * LAS * LAS;
                  double * fockmx = new double[ LAS * LAS ];
                  double * work   = new double[ LAS_pow6  ];
                  double * result = new double[ LAS_pow6  ];
                  for ( int cnt = 0; cnt < LAS_pow6; cnt++ ){ result[ cnt ] = 0.0; }
                  ham->readfock( molcas_fock, fockmx, true );
                  CheMPS2::CASSCF::fock_dot_4rdm( fockmx, dmrgsolver, ham, 0, 0, work, result, false, false );

                  result_filename.str("");
                  result_filename << molcas_f4rdm << ".r" << state;
                  CheMPS2::ThreeDM::save_HAM_generic( result_filename.str(), LAS, "F.4-RDM", result );
                  delete [] fockmx;
                  delete [] work;
                  delete [] result;
               }
               if (( print_corr ) && ( state == excitation )){ dmrgsolver->getCorrelations()->Print(); }
            }
         }
      }

      // Clean up
      if ( CheMPS2::DMRG_storeRenormOptrOnDisk ){ dmrgsolver->deleteStoredOperators(); }
      delete dmrgsolver;
      delete prob;

   } else {

      CheMPS2::CASSCF koekoek( ham, NULL, NULL, nocc_parsed, nact_parsed, nvir_parsed, tmp_folder );

      const int root_num = excitation + 1;
      CheMPS2::DMRGSCFoptions * scf_options = new CheMPS2::DMRGSCFoptions();
      scf_options->setDoDIIS( true );
      scf_options->setDIISGradientBranch( scf_diis_thr );
      scf_options->setStoreDIIS( true );
      scf_options->setMaxIterations( scf_max_iter );
      scf_options->setGradientThreshold( scf_grad_thr );
      scf_options->setStoreUnitary( true );
      scf_options->setStateAveraging( scf_state_avg );
      if ( scf_active_space == 'I' ){ scf_options->setWhichActiveSpace( 0 ); }
      if ( scf_active_space == 'N' ){ scf_options->setWhichActiveSpace( 1 ); }
      if ( scf_active_space == 'L' ){ scf_options->setWhichActiveSpace( 2 ); }
      if ( scf_active_space == 'F' ){ scf_options->setWhichActiveSpace( 3 ); }
      scf_options->setDumpCorrelations( print_corr );
      scf_options->setStartLocRandom( true );

      const double E_CASSCF = koekoek.solve( nelectrons, multiplicity - 1, irrep, opt_scheme, root_num, scf_options );
      double E_CASPT2 = 0.0;
      if ( caspt2_calc ){
         E_CASPT2 = koekoek.caspt2( nelectrons, multiplicity - 1, irrep, opt_scheme, root_num, scf_options, caspt2_ipea, caspt2_imag, ( caspt2_orbs == 'P' ), caspt2_checkpt, caspt2_cumul );
         if ( am_i_master ){
            cout << "E_CASSCF + E_CASPT2 = E_0 + E_1 + E_2 = " << E_CASSCF + E_CASPT2 << endl;
         }
      }

      if ( ( am_i_master ) && ( scf_molden.length() > 0 ) ){
         int * norb_ham = new int[ Symmhelper.getNumberOfIrreps() ];
         for ( int irrep = 0; irrep < Symmhelper.getNumberOfIrreps(); irrep++ ){ norb_ham[ irrep ] = 0; }
         for ( int orb = 0; orb < ham->getL(); orb++ ){ norb_ham[ ham->getOrbitalIrrep( orb ) ] += 1; }
         CheMPS2::Molden * molden_rotator = new CheMPS2::Molden( ham->getL(), Symmhelper.getGroupNumber(), norb_ham );
         molden_rotator->read_molden( scf_molden );
         molden_rotator->read_unitary( scf_options->getUnitaryStorageName() );
         molden_rotator->print( scf_molden, scf_molden + ".rotated" );
         delete [] norb_ham;
         delete molden_rotator;
      }

      // Clean up
      if (( scf_options->getStoreDIIS() ) && ( am_i_master )){ koekoek.deleteStoredDIIS( scf_options->getDIISStorageName() ); }
      delete scf_options;

   }

   delete [] nocc_parsed;
   delete [] nact_parsed;
   delete [] nvir_parsed;
   delete opt_scheme;
   delete ham;

   return clean_exit( 0 );

}



