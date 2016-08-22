!
!  MOLCAS conventions (counting starts from 1)
!
!  'c1 '  ag
!  'ci '  ag  au
!  'c2 '  a   b
!  'cs '  a'  a"
!  'd2 '  a   b2  b1  b3
!  'c2v'  a1  b1  a2  b2
!  'c2h'  ag  bg  au  bu
!  'd2h'  ag  b2g b1g b3g au  b2u b1u b3u
!
!  MOLPRO conventions (counting starts from 1)
!
!  'c1 '  ag
!  'ci '  ag  au
!  'c2 '  a   b
!  'cs '  a'  a"
!  'd2 '  a   b3  b2  b1
!  'c2v'  a1  b1  b2  a2
!  'c2h'  ag  au  bu  bg
!  'd2h'  ag  b3u b2u b1g b1u b2g b3g au
!
!  PSI4 conventions (counting starts from 0)
!
!  'c1 '  ag
!  'ci '  ag  au
!  'c2 '  a   b
!  'cs '  a'  a"
!  'd2 '  a   b1  b2  b3
!  'c2v'  a1  a2  b1  b2
!  'c2h'  ag  bg  au  bu
!  'd2h'  ag  b1g b2g b3g au  b1u b2u b3u
!

subroutine group_psi4number( groupname, psi4number )

  implicit none
  character(3), intent(in)  :: groupname
  integer,      intent(out) :: psi4number

  psi4number = -1 ! Error output

  if ( groupname .eq. 'c1 ' ) psi4number = 0
  if ( groupname .eq. 'ci ' ) psi4number = 1
  if ( groupname .eq. 'c2 ' ) psi4number = 2
  if ( groupname .eq. 'cs ' ) psi4number = 3
  if ( groupname .eq. 'd2 ' ) psi4number = 4
  if ( groupname .eq. 'c2v' ) psi4number = 5
  if ( groupname .eq. 'c2h' ) psi4number = 6
  if ( groupname .eq. 'd2h' ) psi4number = 7

end subroutine

subroutine molpro2psi( groupname, conversion )

  implicit none
  character(3), intent(in)  :: groupname
  integer,      intent(out) :: conversion(8)
  integer,      parameter   :: X = -1

  conversion(1:8) = (/ X, X, X, X, X, X, X, X /) ! Error output

  if ( groupname .eq. 'c1 ' ) conversion(1:8) = (/ 1, X, X, X, X, X, X, X /)
  if ( groupname .eq. 'ci ' ) conversion(1:8) = (/ 1, 2, X, X, X, X, X, X /)
  if ( groupname .eq. 'c2 ' ) conversion(1:8) = (/ 1, 2, X, X, X, X, X, X /)
  if ( groupname .eq. 'cs ' ) conversion(1:8) = (/ 1, 2, X, X, X, X, X, X /)
  if ( groupname .eq. 'd2 ' ) conversion(1:8) = (/ 1, 4, 3, 2, X, X, X, X /)
  if ( groupname .eq. 'c2v' ) conversion(1:8) = (/ 1, 3, 4, 2, X, X, X, X /)
  if ( groupname .eq. 'c2h' ) conversion(1:8) = (/ 1, 3, 4, 2, X, X, X, X /)
  if ( groupname .eq. 'd2h' ) conversion(1:8) = (/ 1, 8, 7, 2, 6, 3, 4, 5 /)

end subroutine

subroutine psi2molpro( groupname, conversion )

  implicit none
  character(3), intent(in)  :: groupname
  integer,      intent(out) :: conversion(8)
  integer,      parameter   :: X = -1

  conversion(1:8) = (/ X, X, X, X, X, X, X, X /) ! Error output

  if ( groupname .eq. 'c1 ' ) conversion(1:8) = (/ 1, X, X, X, X, X, X, X /)
  if ( groupname .eq. 'ci ' ) conversion(1:8) = (/ 1, 2, X, X, X, X, X, X /)
  if ( groupname .eq. 'c2 ' ) conversion(1:8) = (/ 1, 2, X, X, X, X, X, X /)
  if ( groupname .eq. 'cs ' ) conversion(1:8) = (/ 1, 2, X, X, X, X, X, X /)
  if ( groupname .eq. 'd2 ' ) conversion(1:8) = (/ 1, 4, 3, 2, X, X, X, X /)
  if ( groupname .eq. 'c2v' ) conversion(1:8) = (/ 1, 4, 2, 3, X, X, X, X /)
  if ( groupname .eq. 'c2h' ) conversion(1:8) = (/ 1, 4, 2, 3, X, X, X, X /)
  if ( groupname .eq. 'd2h' ) conversion(1:8) = (/ 1, 4, 6, 7, 8, 5, 3, 2 /)

end subroutine

subroutine molcas2molpro( groupname, conversion )

  implicit none
  character(3), intent(in)  :: groupname
  integer,      intent(out) :: conversion(8)
  integer,      parameter   :: X = -1

  conversion(1:8) = (/ X, X, X, X, X, X, X, X /) ! Error output

  if ( groupname .eq. 'c1 ' ) conversion(1:8) = (/ 1, X, X, X, X, X, X, X /)
  if ( groupname .eq. 'ci ' ) conversion(1:8) = (/ 1, 2, X, X, X, X, X, X /)
  if ( groupname .eq. 'c2 ' ) conversion(1:8) = (/ 1, 2, X, X, X, X, X, X /)
  if ( groupname .eq. 'cs ' ) conversion(1:8) = (/ 1, 2, X, X, X, X, X, X /)
  if ( groupname .eq. 'd2 ' ) conversion(1:8) = (/ 1, 3, 4, 2, X, X, X, X /)
  if ( groupname .eq. 'c2v' ) conversion(1:8) = (/ 1, 2, 4, 3, X, X, X, X /)
  if ( groupname .eq. 'c2h' ) conversion(1:8) = (/ 1, 4, 2, 3, X, X, X, X /)
  if ( groupname .eq. 'd2h' ) conversion(1:8) = (/ 1, 6, 4, 7, 8, 3, 5, 2 /)

end subroutine

subroutine molpro2molcas( groupname, conversion )

  implicit none
  character(3), intent(in)  :: groupname
  integer,      intent(out) :: conversion(8)
  integer,      parameter   :: X = -1

  conversion(1:8) = (/ X, X, X, X, X, X, X, X /) ! Error output

  if ( groupname .eq. 'c1 ' ) conversion(1:8) = (/ 1, X, X, X, X, X, X, X /)
  if ( groupname .eq. 'ci ' ) conversion(1:8) = (/ 1, 2, X, X, X, X, X, X /)
  if ( groupname .eq. 'c2 ' ) conversion(1:8) = (/ 1, 2, X, X, X, X, X, X /)
  if ( groupname .eq. 'cs ' ) conversion(1:8) = (/ 1, 2, X, X, X, X, X, X /)
  if ( groupname .eq. 'd2 ' ) conversion(1:8) = (/ 1, 4, 2, 3, X, X, X, X /)
  if ( groupname .eq. 'c2v' ) conversion(1:8) = (/ 1, 2, 4, 3, X, X, X, X /)
  if ( groupname .eq. 'c2h' ) conversion(1:8) = (/ 1, 3, 4, 2, X, X, X, X /)
  if ( groupname .eq. 'd2h' ) conversion(1:8) = (/ 1, 8, 6, 3, 7, 2, 4, 5 /)

end subroutine

subroutine molcas2psi( groupname, conversion )

  implicit none
  character(3), intent(in)  :: groupname
  integer,      intent(out) :: conversion(8)
  integer,      parameter   :: X = -1

  conversion(1:8) = (/ X, X, X, X, X, X, X, X /) ! Error output

  if ( groupname .eq. 'c1 ' ) conversion(1:8) = (/ 1, X, X, X, X, X, X, X /)
  if ( groupname .eq. 'ci ' ) conversion(1:8) = (/ 1, 2, X, X, X, X, X, X /)
  if ( groupname .eq. 'c2 ' ) conversion(1:8) = (/ 1, 2, X, X, X, X, X, X /)
  if ( groupname .eq. 'cs ' ) conversion(1:8) = (/ 1, 2, X, X, X, X, X, X /)
  if ( groupname .eq. 'd2 ' ) conversion(1:8) = (/ 1, 3, 2, 4, X, X, X, X /)
  if ( groupname .eq. 'c2v' ) conversion(1:8) = (/ 1, 3, 2, 4, X, X, X, X /)
  if ( groupname .eq. 'c2h' ) conversion(1:8) = (/ 1, 2, 3, 4, X, X, X, X /)
  if ( groupname .eq. 'd2h' ) conversion(1:8) = (/ 1, 3, 2, 4, 5, 7, 6, 8 /)

end subroutine

subroutine psi2molcas( groupname, conversion )

  implicit none
  character(3), intent(in)  :: groupname
  integer,      intent(out) :: conversion(8)
  integer,      parameter   :: X = -1

  conversion(1:8) = (/ X, X, X, X, X, X, X, X /) ! Error output

  if ( groupname .eq. 'c1 ' ) conversion(1:8) = (/ 1, X, X, X, X, X, X, X /)
  if ( groupname .eq. 'ci ' ) conversion(1:8) = (/ 1, 2, X, X, X, X, X, X /)
  if ( groupname .eq. 'c2 ' ) conversion(1:8) = (/ 1, 2, X, X, X, X, X, X /)
  if ( groupname .eq. 'cs ' ) conversion(1:8) = (/ 1, 2, X, X, X, X, X, X /)
  if ( groupname .eq. 'd2 ' ) conversion(1:8) = (/ 1, 3, 2, 4, X, X, X, X /)
  if ( groupname .eq. 'c2v' ) conversion(1:8) = (/ 1, 3, 2, 4, X, X, X, X /)
  if ( groupname .eq. 'c2h' ) conversion(1:8) = (/ 1, 2, 3, 4, X, X, X, X /)
  if ( groupname .eq. 'd2h' ) conversion(1:8) = (/ 1, 3, 2, 4, 5, 7, 6, 8 /)

end subroutine

program main

  implicit none

  integer :: conv_psi_cas(8)
  integer :: conv_cas_psi(8)
  integer :: conv_psi_pro(8)
  integer :: conv_pro_psi(8)
  integer :: conv_pro_cas(8)
  integer :: conv_cas_pro(8)
  integer :: group, num_irreps, irrep
  character(3) :: groupname

  do group=0,7

    if ( group == 0 ) groupname = 'c1 '
    if ( group == 1 ) groupname = 'ci '
    if ( group == 2 ) groupname = 'c2 '
    if ( group == 3 ) groupname = 'cs '
    if ( group == 4 ) groupname = 'd2 '
    if ( group == 5 ) groupname = 'c2v'
    if ( group == 6 ) groupname = 'c2h'
    if ( group == 7 ) groupname = 'd2h'

    if ( group == 0 ) num_irreps = 1
    if ( group == 1 ) num_irreps = 2
    if ( group == 2 ) num_irreps = 2
    if ( group == 3 ) num_irreps = 2
    if ( group == 4 ) num_irreps = 4
    if ( group == 5 ) num_irreps = 4
    if ( group == 6 ) num_irreps = 4
    if ( group == 7 ) num_irreps = 8

    call psi2molcas(    groupname, conv_psi_cas )
    call molcas2psi(    groupname, conv_cas_psi )
    call psi2molpro(    groupname, conv_psi_pro )
    call molpro2psi(    groupname, conv_pro_psi )
    call molpro2molcas( groupname, conv_pro_cas )
    call molcas2molpro( groupname, conv_cas_pro )

    write( 6, '(A8, A3)' ) "Group = ", groupname
    do irrep=1,num_irreps
      write( 6, * ) conv_cas_psi( conv_psi_cas( irrep ) )
      write( 6, * ) conv_pro_psi( conv_psi_pro( irrep ) )
      write( 6, * ) conv_pro_cas( conv_cas_pro( irrep ) )
      write( 6, * ) conv_cas_psi( conv_pro_cas( conv_psi_pro( irrep ) ) )
      write( 6, * ) conv_pro_psi( conv_cas_pro( conv_psi_cas( irrep ) ) )
      write( 6, * ) "*************************************"
    end do
  end do

end program main


