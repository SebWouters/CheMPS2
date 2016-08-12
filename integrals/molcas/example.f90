
PROGRAM main

  USE HDF5
  USE ISO_C_BINDING

  IMPLICIT NONE

  CHARACTER(LEN=14), PARAMETER :: file_2rdm  = "molcas_2rdm.h5"
  CHARACTER(LEN=14), PARAMETER :: file_3rdm  = "molcas_3rdm.h5"
  CHARACTER(LEN=15), PARAMETER :: file_f4rdm = "molcas_f4rdm.h5"

  INTEGER( HID_T )   :: file_h5, group_h5, space_h5, dset_h5 ! Handles
  INTEGER( HSIZE_T ) :: dimension
  INTEGER            :: hdferr
  TYPE( C_PTR )      :: f_ptr

  INTEGER, PARAMETER :: L = 10 ! Number of orbitals
  INTEGER, PARAMETER :: Lpow4  = L * L * L * L
  INTEGER, PARAMETER :: Lpow6  = Lpow4  * L * L

  REAL*8, DIMENSION( 1 : Lpow4 ), TARGET :: two_rdm
  REAL*8, DIMENSION( 1 : Lpow6 ), TARGET :: three_rdm
  REAL*8, DIMENSION( 1 : Lpow6 ), TARGET :: f_dot_4rdm

  CALL h5open_f( hdferr )

  ! Process the 2-RDM
  CALL h5fopen_f( file_2rdm, H5F_ACC_RDONLY_F, file_h5, hdferr )
  CALL h5gopen_f( file_h5, "2-RDM", group_h5, hdferr )
  CALL h5dopen_f( group_h5, "elements", dset_h5, hdferr )
  CALL h5dget_space_f( dset_h5, space_h5, hdferr )
  f_ptr = C_LOC( two_rdm( 1 ) )
  CALL h5dread_f( dset_h5, H5T_NATIVE_DOUBLE, f_ptr, hdferr )
  CALL h5dclose_f( dset_h5 , hdferr )
  CALL h5sclose_f( space_h5, hdferr )
  CALL h5gclose_f( group_h5, hdferr )
  CALL h5fclose_f( file_h5,  hdferr )

  WRITE( *, * ) "2-RDM"
  WRITE( *, * ) "   two_rdm( 1 + i + L * ( j + L * ( k + L * l ) ) )"
  WRITE( *, * ) "     = sum_{sigma,tau} < a^+_{i,sigma} a^+_{j,tau} a_{l,tau} a_{k,sigma} >"
  WRITE( *, * ) "     = 2-RDM[ i, j, k, l ]"
  WRITE( *, * ) "   2-RDM[ 0, 0, 0, 0 ] = ", two_rdm( 1 )
  WRITE( *, * ) "   2-RDM[ 9, 9, 9, 9 ] = ", two_rdm( Lpow4 )

  ! Process the 3-RDM
  CALL h5fopen_f( file_3rdm, H5F_ACC_RDONLY_F, file_h5, hdferr )
  CALL h5gopen_f( file_h5, "3-RDM", group_h5, hdferr )
  CALL h5dopen_f( group_h5, "elements", dset_h5, hdferr )
  CALL h5dget_space_f( dset_h5, space_h5, hdferr )
  f_ptr = C_LOC( three_rdm( 1 ) )
  CALL h5dread_f( dset_h5, H5T_NATIVE_DOUBLE, f_ptr, hdferr )
  CALL h5dclose_f( dset_h5 , hdferr )
  CALL h5sclose_f( space_h5, hdferr )
  CALL h5gclose_f( group_h5, hdferr )
  CALL h5fclose_f( file_h5,  hdferr )

  WRITE( *, * ) "3-RDM"
  WRITE( *, * ) "   three_rdm( 1 + i + L * ( j + L * ( k + L * ( l + L * ( m + L * n ) ) ) ) )"
  WRITE( *, * ) "     = sum_{sigma,tau,z} < a^+_{i,sigma} a^+_{j,tau} a^+_{k,z} a_{n,z} a_{m,tau} a_{l,sigma} >"
  WRITE( *, * ) "     = 3-RDM[ i, j, k, l, m, n ]"
  WRITE( *, * ) "   3-RDM[ 1, 0, 0, 1, 0, 0 ] = ", three_rdm( 1001 + 1 )
  WRITE( *, * ) "   3-RDM[ 9, 9, 8, 9, 9, 8 ] = ", three_rdm( 998998 + 1 )

  ! Process the Fock oeprator contracted with the 4-RDM
  CALL h5fopen_f( file_f4rdm, H5F_ACC_RDONLY_F, file_h5, hdferr )
  CALL h5gopen_f( file_h5, "F.4-RDM", group_h5, hdferr )
  CALL h5dopen_f( group_h5, "elements", dset_h5, hdferr )
  CALL h5dget_space_f( dset_h5, space_h5, hdferr )
  f_ptr = C_LOC( f_dot_4rdm( 1 ) )
  CALL h5dread_f( dset_h5, H5T_NATIVE_DOUBLE, f_ptr, hdferr )
  CALL h5dclose_f( dset_h5 , hdferr )
  CALL h5sclose_f( space_h5, hdferr )
  CALL h5gclose_f( group_h5, hdferr )
  CALL h5fclose_f( file_h5,  hdferr )

  WRITE( *, * ) "F.4-RDM"
  WRITE( *, * ) "   f_dot_4rdm( 1 + i + L * ( j + L * ( k + L * ( l + L * ( m + L * n ) ) ) ) )"
  WRITE( *, * ) "     = sum_{p,q} sum_{s,tau,z} < a^+_{i,s} a^+_{j,tau} a^+_{k,z} E_{p,q} a_{n,z} a_{m,tau} a_{l,s} > * Fock[p,q]"
  WRITE( *, * ) "     = F.4-RDM[ i, j, k, l, m, n ]"
  WRITE( *, * ) "   F.4-RDM[ 1, 0, 0, 1, 0, 0 ] = ", f_dot_4rdm( 1001 + 1 )
  WRITE( *, * ) "   F.4-RDM[ 9, 9, 8, 9, 9, 8 ] = ", f_dot_4rdm( 998998 + 1 )

END PROGRAM main

