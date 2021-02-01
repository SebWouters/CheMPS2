/*
   CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
   Copyright (C) 2013-2021 Sebastian Wouters

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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sys/stat.h>

#include <hdf5.h>

#include "Irreps.h"
#include "TwoIndex.h"
#include "FourIndex.h"
#include "Hamiltonian.h"

using std::cout;
using std::endl;
using std::string;
using std::ifstream;

CheMPS2::Hamiltonian::Hamiltonian(const int Norbitals, const int nGroup, const int * OrbIrreps){

   L = Norbitals;
   assert( nGroup>=0 );
   assert( nGroup<=7 );
   SymmInfo.setGroup(nGroup);
   
   orb2irrep = new int[L];
   orb2indexSy = new int[L];
   int nIrreps = SymmInfo.getNumberOfIrreps();
   irrep2num_orb = new int[nIrreps];
   for (int cnt=0; cnt<nIrreps; cnt++) irrep2num_orb[cnt] = 0;
   for (int cnt=0; cnt<L; cnt++){
      assert( OrbIrreps[cnt]>=0 );
      assert( OrbIrreps[cnt]<nIrreps );
      orb2irrep[cnt] = OrbIrreps[cnt];
      orb2indexSy[cnt] = irrep2num_orb[orb2irrep[cnt]];
      irrep2num_orb[orb2irrep[cnt]]++;
   }
   
   Econst = 0.0;
   Tmat = new TwoIndex(SymmInfo.getGroupNumber(),irrep2num_orb);
   Vmat = new FourIndex(SymmInfo.getGroupNumber(),irrep2num_orb);

}

CheMPS2::Hamiltonian::Hamiltonian( const string filename, const int psi4groupnumber ){

    SymmInfo.setGroup( psi4groupnumber );
    CreateAndFillFromFCIDUMP( filename );

}

CheMPS2::Hamiltonian::Hamiltonian( const bool fileh5, const string main_file, const string file_tmat, const string file_vmat ){

   if ( fileh5 ){
      CreateAndFillFromH5( main_file, file_tmat, file_vmat );
   } else {
      cout << "CheMPS2::Hamiltonian::Hamiltonian( false, const string , const string , const string ) was deprecated." << endl;
      assert( fileh5 == true );
   }

}

CheMPS2::Hamiltonian::~Hamiltonian(){
   
   delete [] orb2irrep;
   delete [] orb2indexSy;
   delete [] irrep2num_orb;
   delete Tmat;
   delete Vmat;
   
}

int CheMPS2::Hamiltonian::getL() const{ return L; }

int CheMPS2::Hamiltonian::getNGroup() const{ return SymmInfo.getGroupNumber(); }

int CheMPS2::Hamiltonian::getOrbitalIrrep(const int nOrb) const{ return orb2irrep[nOrb]; }

void CheMPS2::Hamiltonian::setEconst(const double val){ Econst = val; }

double CheMPS2::Hamiltonian::getEconst() const{ return Econst; }

void CheMPS2::Hamiltonian::setTmat(const int index1, const int index2, const double val){

   assert( orb2irrep[index1]==orb2irrep[index2] );
   Tmat->set(orb2irrep[index1], orb2indexSy[index1], orb2indexSy[index2], val);

}

double CheMPS2::Hamiltonian::getTmat(const int index1, const int index2) const{

   if (orb2irrep[index1]==orb2irrep[index2]){
      return Tmat->get(orb2irrep[index1], orb2indexSy[index1], orb2indexSy[index2]);
   }

   return 0.0;
   
}

const CheMPS2::TwoIndex * CheMPS2::Hamiltonian::getTmat(){ return Tmat; }

void CheMPS2::Hamiltonian::setVmat(const int index1, const int index2, const int index3, const int index4, const double val){

   assert( Irreps::directProd(orb2irrep[index1],orb2irrep[index2]) == Irreps::directProd(orb2irrep[index3],orb2irrep[index4]) );
   Vmat->set(orb2irrep[index1], orb2irrep[index2], orb2irrep[index3], orb2irrep[index4], orb2indexSy[index1], orb2indexSy[index2], orb2indexSy[index3], orb2indexSy[index4], val);

}

void CheMPS2::Hamiltonian::addToVmat(const int index1, const int index2, const int index3, const int index4, const double val){

   assert( Irreps::directProd(orb2irrep[index1],orb2irrep[index2]) == Irreps::directProd(orb2irrep[index3],orb2irrep[index4]) );
   Vmat->add(orb2irrep[index1], orb2irrep[index2], orb2irrep[index3], orb2irrep[index4], orb2indexSy[index1], orb2indexSy[index2], orb2indexSy[index3], orb2indexSy[index4], val);

}

double CheMPS2::Hamiltonian::getVmat(const int index1, const int index2, const int index3, const int index4) const{

   if ( Irreps::directProd(orb2irrep[index1],orb2irrep[index2]) == Irreps::directProd(orb2irrep[index3],orb2irrep[index4]) ){
      return Vmat->get(orb2irrep[index1], orb2irrep[index2], orb2irrep[index3], orb2irrep[index4], orb2indexSy[index1], orb2indexSy[index2], orb2indexSy[index3], orb2indexSy[index4]);
   }

   return 0.0;
   
}

CheMPS2::FourIndex * CheMPS2::Hamiltonian::getVmat(){ return Vmat; }

void CheMPS2::Hamiltonian::save(const string file_parent, const string file_tmat, const string file_vmat) const{

   Tmat->save(file_tmat);
   Vmat->save(file_vmat);

   //The hdf5 file
   hid_t file_id = H5Fcreate(file_parent.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      
      //The data
      hid_t group_id = H5Gcreate(file_id, "/Data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
       
         //The chain length
         hsize_t dimarray       = 1;
         hid_t dataspace_id     = H5Screate_simple(1, &dimarray, NULL);
         hid_t dataset_id       = H5Dcreate(group_id, "L", H5T_STD_I32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
         H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &L);
         
         //The group number
         hsize_t dimarray2      = 1;
         hid_t dataspace_id2    = H5Screate_simple(1, &dimarray2, NULL);
         hid_t dataset_id2      = H5Dcreate(group_id, "nGroup", H5T_STD_I32LE, dataspace_id2, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
         int nGroup = SymmInfo.getGroupNumber();
         H5Dwrite(dataset_id2, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nGroup);
         
         //orb2irrep
         hsize_t dimarray3      = L;
         hid_t dataspace_id3    = H5Screate_simple(1, &dimarray3, NULL);
         hid_t dataset_id3      = H5Dcreate(group_id, "orb2irrep", H5T_STD_I32LE, dataspace_id3, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
         H5Dwrite(dataset_id3, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, orb2irrep);
         
         //Econst
         hsize_t dimarray4      = 1;
         hid_t dataspace_id4    = H5Screate_simple(1, &dimarray4, NULL);
         hid_t dataset_id4      = H5Dcreate(group_id, "Econst", H5T_IEEE_F64LE, dataspace_id4, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
         H5Dwrite(dataset_id4, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Econst);
    
         H5Dclose(dataset_id);
         H5Sclose(dataspace_id);
         H5Dclose(dataset_id2);
         H5Sclose(dataspace_id2);
         H5Dclose(dataset_id3);
         H5Sclose(dataspace_id3);
         H5Dclose(dataset_id4);
         H5Sclose(dataspace_id4);

      H5Gclose(group_id);
      
   H5Fclose(file_id);

}

void CheMPS2::Hamiltonian::read(const string file_parent, const string file_tmat, const string file_vmat){

   Tmat->read(file_tmat);
   Vmat->read(file_vmat);

   //The hdf5 file
   hid_t file_id = H5Fopen(file_parent.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      
      //The data
      hid_t group_id = H5Gopen(file_id, "/Data",H5P_DEFAULT);
       
         //The chain length
         hid_t dataset_id1 = H5Dopen(group_id, "L", H5P_DEFAULT);
         int Lagain;
         H5Dread(dataset_id1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Lagain);
         assert( L==Lagain );
         
         //The group number
         hid_t dataset_id2 = H5Dopen(group_id, "nGroup", H5P_DEFAULT);
         int nGroup;
         H5Dread(dataset_id2, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nGroup);
         assert( nGroup==SymmInfo.getGroupNumber() );
         
         //orb2irrep
         hid_t dataset_id3 = H5Dopen(group_id, "orb2irrep", H5P_DEFAULT);
         int * orb2irrepAgain = new int[L];
         H5Dread(dataset_id3, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, orb2irrepAgain);
         for (int cnt=0; cnt<L; cnt++){
            assert( orb2irrep[cnt]==orb2irrepAgain[cnt] );
         }
         delete [] orb2irrepAgain;
         
         //Econst
         hid_t dataset_id4 = H5Dopen(group_id, "Econst", H5P_DEFAULT);
         H5Dread(dataset_id4, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Econst);
    
         H5Dclose(dataset_id1);
         H5Dclose(dataset_id2);
         H5Dclose(dataset_id3);
         H5Dclose(dataset_id4);

      H5Gclose(group_id);
      
   H5Fclose(file_id);
   
   if (CheMPS2::HAMILTONIAN_debugPrint) debugcheck();

}

void CheMPS2::Hamiltonian::CreateAndFillFromH5(const string file_parent, const string file_tmat, const string file_vmat){

   //The hdf5 file
   hid_t file_id = H5Fopen(file_parent.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      
      //The data
      hid_t group_id = H5Gopen(file_id, "/Data",H5P_DEFAULT);
       
         //The number of orbitals: Set L directly!
         hid_t dataset_id1 = H5Dopen(group_id, "L", H5P_DEFAULT);
         H5Dread(dataset_id1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &L);
         
         //The group number: SymmInfo.setGroup()
         hid_t dataset_id2 = H5Dopen(group_id, "nGroup", H5P_DEFAULT);
         int nGroup_LOADH5;
         H5Dread(dataset_id2, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nGroup_LOADH5);
         SymmInfo.setGroup( nGroup_LOADH5 );
         
         //The irrep of each orbital: Create and fill orb2irrep directly!
         hid_t dataset_id3 = H5Dopen(group_id, "orb2irrep", H5P_DEFAULT);
         orb2irrep = new int[ L ];
         H5Dread(dataset_id3, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, orb2irrep);
    
         H5Dclose(dataset_id1);
         H5Dclose(dataset_id2);
         H5Dclose(dataset_id3);

      H5Gclose(group_id);
      
   H5Fclose(file_id);
   
   orb2indexSy = new int[ L ];
   int nIrreps = SymmInfo.getNumberOfIrreps();
   irrep2num_orb = new int[ nIrreps ];
   for (int irrep = 0; irrep < nIrreps; irrep++){ irrep2num_orb[ irrep ] = 0; }
   for (int orb = 0; orb < L; orb++){
      orb2indexSy[ orb ] = irrep2num_orb[ orb2irrep[ orb ] ];
      irrep2num_orb[ orb2irrep[ orb ] ]++;
   }
   
   Tmat = new TwoIndex(  SymmInfo.getGroupNumber(), irrep2num_orb );
   Vmat = new FourIndex( SymmInfo.getGroupNumber(), irrep2num_orb );

   read(file_parent, file_tmat, file_vmat);

}

void CheMPS2::Hamiltonian::CreateAndFillFromFCIDUMP( const string fcidumpfile ){

    struct stat file_info;
    const bool on_disk = (( fcidumpfile.length() > 0 ) && ( stat( fcidumpfile.c_str(), &file_info ) == 0 ));
    if ( on_disk == false ){
       cout << "CheMPS2::Hamiltonian : Unable to find FCIDUMP file " << fcidumpfile << "!" << endl;
    }
    assert( on_disk );

    cout << "CheMPS2::Hamiltonian : Reading FCIDUMP file " << fcidumpfile << endl;

    const int nIrreps = SymmInfo.getNumberOfIrreps();
    int * psi2molpro = new int[ nIrreps ];
    SymmInfo.symm_psi2molpro( psi2molpro );

    ifstream thefcidump( fcidumpfile.c_str() );
    string line, part;
    int pos, pos2;

    // Get the number of orbitals
    getline( thefcidump, line ); // &FCI NORB= X,NELEC= Y,MS2= Z,   
    pos  = line.find( "NORB" );
    pos  = line.find( "=", pos ); //1
    pos2 = line.find( ",", pos ); //4
    part = line.substr( pos+1, pos2-pos-1 );
    L = atoi( part.c_str() );
    if ( CheMPS2::HAMILTONIAN_debugPrint ){ cout << "The number of orbitals <<" << part << ">> or " << L << "." << endl; }

    // Get the orbital irreps in psi4 convention (XOR, see Irreps.h).
    orb2irrep = new int[ L ];
    getline( thefcidump, line ); //  ORBSYM=A,B,C,D,
    getline( thefcidump, part );
    while ( part.find( "ISYM" ) == string::npos ){
       pos = line.find("\n");
       if ( pos != string::npos ){ line.erase( pos ); }
       pos = part.find(" ");
       if ( pos != string::npos ){ part.erase( pos, 1 ); }
       line.append( part );
       getline( thefcidump, part );
    }
    pos = line.find( "ORBSYM" );
    pos = line.find( "=", pos ); //1
    for ( int orb = 0; orb < L; orb++ ){
        pos2 = line.find( ",", pos+1 ); //3
        part = line.substr( pos+1, pos2-pos-1 );
        const int molproirrep = atoi( part.c_str() );
        if ( CheMPS2::HAMILTONIAN_debugPrint ){ cout << "This molpro irrep <<" << part << ">> or " << molproirrep << "." << endl; }
        orb2irrep[ orb ] = -1;
        for ( int irrep = 0; irrep < nIrreps; irrep++ ){
            if ( molproirrep == psi2molpro[ irrep ] ){
                orb2irrep[ orb ] = irrep;
            }
        }
        assert( orb2irrep[ orb ] != -1 );
        pos = pos2;
    }

    getline( thefcidump, line ); // /
    assert( line.size() < 16 );

    orb2indexSy = new int[ L ];
    irrep2num_orb = new int[ nIrreps ];
    for ( int cnt = 0; cnt < nIrreps; cnt++){ irrep2num_orb[cnt] = 0; }
    for ( int cnt = 0; cnt < L; cnt++){
        orb2indexSy[cnt] = irrep2num_orb[orb2irrep[cnt]];
        irrep2num_orb[orb2irrep[cnt]]++;
    }
    Tmat = new TwoIndex(  SymmInfo.getGroupNumber(), irrep2num_orb ); // Constructor ends with Clear(); call
    Vmat = new FourIndex( SymmInfo.getGroupNumber(), irrep2num_orb ); // Constructor ends with Clear(); call

    // Read the Hamiltonian in
    bool stop = false;
    while ( stop == false ){
    
        getline( thefcidump, line ); // value i1 i2 i3 i4
        pos  = line.find( " " );
        pos2 = line.find( "." );
        pos2 = line.find( " ", pos2 );
        part = line.substr( pos, pos2-pos );
        if ( CheMPS2::HAMILTONIAN_debugPrint ){ cout << "Next line: <<" << part << ">> <<"; }
        const double value = atof( part.c_str() );
        pos  = pos2;
        while ( line.substr( pos, 1 ).compare(" ")==0 ){ pos++; }
        pos2 = line.find( " ", pos );
        part = line.substr( pos, pos2-pos );
        if ( CheMPS2::HAMILTONIAN_debugPrint ){ cout << part << ">> <<"; }
        const int index1 = atoi( part.c_str() );
        pos  = pos2;
        while ( line.substr( pos, 1 ).compare(" ")==0 ){ pos++; }
        pos2 = line.find( " ", pos );
        part = line.substr( pos, pos2-pos );
        if ( CheMPS2::HAMILTONIAN_debugPrint ){ cout << part << ">> <<"; }
        const int index2 = atoi( part.c_str() );
        pos  = pos2;
        while ( line.substr( pos, 1 ).compare(" ")==0 ){ pos++; }
        pos2 = line.find( " ", pos );
        part = line.substr( pos, pos2-pos );
        if ( CheMPS2::HAMILTONIAN_debugPrint ){ cout << part << ">> <<"; }
        const int index3 = atoi( part.c_str() );
        pos  = pos2;
        while ( line.substr( pos, 1 ).compare(" ")==0 ){ pos++; }
        part = line.substr( pos, line.size()-pos );
        if ( CheMPS2::HAMILTONIAN_debugPrint ){ cout << part << ">>." << endl; }
        const int index4 = atoi( part.c_str() );
        if ( CheMPS2::HAMILTONIAN_debugPrint ){
           cout << "Same line: " << value << " " << index1 << " " << index2 << " " << index3 << " " << index4 << endl;
        }
        
        if ( index4 != 0 ){
            setVmat( index1-1, index3-1, index2-1, index4-1, value ); // From chemists to physicist notation!
        } else {
            if ( index2 != 0 ){ setTmat( index1-1, index2-1, value ); }
            else {
                Econst = value;
                stop = true;
            }
        }
    
    }
    
    if ( CheMPS2::HAMILTONIAN_debugPrint ){ debugcheck(); }

    delete [] psi2molpro;
    thefcidump.close();

    cout << "CheMPS2::Hamiltonian : Finished reading FCIDUMP file " << fcidumpfile << endl;

}

void CheMPS2::Hamiltonian::readfock( const string fockfile, double * fockmx, const bool printinfo ) const{

/************************
 *   FOCK file format   *
************************************
 &FOCK NACT= X,
  ORBSYM=X,X,X,X,
 /
  1.1234567890123456E-03   I   J
  1.1234567890123456E-03   I   J
************************************
     Remark: ORBSYM are the MOLPRO irreps !!!
*/

    struct stat file_info;
    const bool on_disk = (( fockfile.length() > 0 ) && ( stat( fockfile.c_str(), &file_info ) == 0 ));
    if ( on_disk == false ){
       cout << "Unable to find the FOCK file " << fockfile << "!" << endl;
    }
    assert( on_disk );

    cout << "Reading FOCK file " << fockfile << endl;

    const int nIrreps = SymmInfo.getNumberOfIrreps();
    int * psi2molpro = new int[ nIrreps ];
    SymmInfo.symm_psi2molpro( psi2molpro );

    ifstream thedump( fockfile.c_str() );
    string line, part;
    int pos, pos2;

    // Double check the number of orbitals
    getline( thedump, line ); // &FOCK NACT= X,
    pos  = line.find( "NACT" );
    pos  = line.find( "=", pos ); //1
    pos2 = line.find( ",", pos ); //4
    part = line.substr( pos+1, pos2-pos-1 );
    const int LAS = atoi( part.c_str() );
    if ( LAS != getL() ){
       cout << "The number of orbitals in the FOCK file and Hamiltonian object does not match!" << endl;
       assert( LAS == getL() );
    }

    // Double check the orbital irreps
    getline( thedump, line ); //  ORBSYM=A,B,C,D,
    getline( thedump, part );
    while ( part.find( "/" ) == string::npos ){
        pos = line.find( "\n" );
        if ( pos != string::npos ){ line.erase( pos ); }
        pos = part.find( " " );
        if ( pos != string::npos ){ part.erase( pos, 1 ); }
        line.append( part );
        getline( thedump, part );
    }
    pos = line.find( "ORBSYM" );
    pos = line.find( "=", pos ); //1
    for ( int orb = 0; orb < LAS; orb++ ){
        pos2 = line.find( ",", pos+1 ); //3
        part = line.substr( pos+1, pos2-pos-1 );
        const int molproirrep = atoi( part.c_str() );
        if ( molproirrep != psi2molpro[ getOrbitalIrrep( orb ) ] ){
            cout << "The irrep of orbital " << orb << " in the FOCK file and Hamiltonian object does not match!" << endl;
            assert( molproirrep == psi2molpro[ getOrbitalIrrep( orb ) ] );
        }
        pos = pos2;
    }

    // Read the FOCK matrix in
    for ( int cnt = 0; cnt < LAS * LAS; cnt++ ){ fockmx[ cnt ] = 0.0; }
    while ( getline( thedump, line ) ){ // value i1 i2

        if ( line.length() > 2 ){

            pos  = line.find( " " );
            pos2 = line.find( "." );
            pos2 = line.find( " ", pos2 );
            part = line.substr( pos, pos2-pos );
            const double value = atof( part.c_str() );
            pos  = pos2;
            while ( line.substr( pos, 1 ).compare(" ") == 0 ){ pos++; }
            pos2 = line.find( " ", pos );
            part = line.substr( pos, pos2-pos );
            const int index1 = atoi( part.c_str() );
            pos  = pos2;
            while ( line.substr( pos, 1 ).compare(" ") == 0 ){ pos++; }
            pos2 = line.find( " ", pos );
            part = line.substr( pos, pos2-pos );
            const int index2 = atoi( part.c_str() );

            if ( printinfo ){
                cout << "        Processed FOCK( " << index1 << ", " << index2 << " ) = " << value << endl;
            }

            const int orb1 = index1 - 1;
            const int orb2 = index2 - 1;

            if ( getOrbitalIrrep( orb1 ) != getOrbitalIrrep( orb2 ) ){
                cout << "In the FOCK file a specific value is given for orbitals " << index1 << " and " << index2 << " which have different irreps. This is not allowed!" << endl;
                assert( getOrbitalIrrep( orb1 ) == getOrbitalIrrep( orb2 ) );
            }

            fockmx[ orb1 + LAS * orb2 ] = value;
            fockmx[ orb2 + LAS * orb1 ] = value;

        }
    }

    delete [] psi2molpro;
    thedump.close();

    cout << "Finished reading FOCK file " << fockfile << endl;

}

void CheMPS2::Hamiltonian::writeFCIDUMP( const string fcidumpfile, const int Nelec, const int TwoS, const int TargetIrrep ) const{

   int * psi2molpro = new int[ SymmInfo.getNumberOfIrreps() ];
   SymmInfo.symm_psi2molpro( psi2molpro );
   
   FILE * capturing;
   capturing = fopen( fcidumpfile.c_str(), "w" ); // "w" with fopen means truncate file
   fprintf( capturing, " &FCI NORB= %d,NELEC= %d,MS2= %d,\n", getL(), Nelec, TwoS );
   fprintf( capturing, "  ORBSYM=" );
   for (int orb=0; orb<getL(); orb++){
      fprintf( capturing, "%d,", psi2molpro[getOrbitalIrrep(orb)] );
   }
   fprintf( capturing, "\n  ISYM=%d,\n /\n", psi2molpro[TargetIrrep] );
   delete [] psi2molpro;
   
   for (int p=0; p<getL(); p++){
      for (int q=0; q<=p; q++){ // p>=q
         const int irrep_pq = Irreps::directProd( getOrbitalIrrep(p), getOrbitalIrrep(q) );
         for (int r=0; r<=p; r++){ // p>=r
            for (int s=0; s<=r; s++){ // r>=s
               const int irrep_rs = Irreps::directProd( getOrbitalIrrep(r), getOrbitalIrrep(s) );
               if ( irrep_pq == irrep_rs ){
                  if ( ( p > r ) || ( ( p == r ) && ( q >= s ) ) ){
                     fprintf( capturing, " % 23.16E %3d %3d %3d %3d\n", getVmat(p,r,q,s), p+1, q+1, r+1, s+1 );
                  }
               }
            }
         }
      }
   }

   for (int p=0; p<getL(); p++){
      for (int q=0; q<=p; q++){ // p>=q
         if ( getOrbitalIrrep(p) == getOrbitalIrrep(q) ){
            fprintf( capturing, " % 23.16E %3d %3d %3d %3d\n", getTmat(p,q), p+1, q+1, 0, 0 );
         }
      }
   }
   
   fprintf( capturing, " % 23.16E %3d %3d %3d %3d", getEconst(), 0, 0, 0, 0 );
   fclose( capturing );
   cout << "Created the FCIDUMP file " << fcidumpfile << "." << endl;

}

void CheMPS2::Hamiltonian::debugcheck() const{

   cout << "Econst = " << Econst << endl;
   
   double test = 0.0;
   double test2 = 0.0;
   double test3 = 0.0;
   for (int i=0; i<L; i++){
      test3 += getTmat(i,i);
      for (int j=0; j<L; j++){
         test += getTmat(i,j);
         if (i<=j) test2 += getTmat(i,j);
      }
   }
   cout << "1-electron integrals: Trace                  : " << test3 << endl;
   cout << "1-electron integrals: Sum over all elements  : " << test << endl;
   cout << "1-electron integrals: Sum over Tij with i<=j : " << test2 << endl;
      
   test = 0.0;
   test2 = 0.0;
   test3 = 0.0;
   for (int i=0; i<L; i++){
      test3 += getVmat(i,i,i,i);
      for (int j=0; j<L; j++){
         for (int k=0; k<L; k++){
            for (int l=0; l<L; l++){
               test += getVmat(i,j,k,l);
               if ((i<=j) && (j<=k) && (k<=l)) test2 += getVmat(i,j,k,l);
            }
         }
      }
   }
   cout << "2-electron integrals: Trace                          : " << test3 << endl;
   cout << "2-electron integrals: Sum over all elements          : " << test << endl;
   cout << "2-electron integrals: Sum over Vijkl with i<=j<=k<=l : " << test2 << endl;

}


