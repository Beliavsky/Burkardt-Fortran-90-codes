program main

!*****************************************************************************80
!
!! fem3d_pack_test() tests fem3d_pack().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 March 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM3D_PACK_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the FEM3D_PACK library.'

  call basis_mn_tet4_test ( )
  call basis_mn_tet10_test ( )
  call basis_brick8_test ( )
  call basis_brick20_test ( )
  call basis_brick27_test ( )
  call test03 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM3D_PACK_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests PHYSICAL_TO_REFERENCE_TET4 and REFERENCE_TO_PHYSICAL_TET4.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 10

  integer j
  real ( kind = rk ) phy(3,n)
  real ( kind = rk ) ref(3,n)
  real ( kind = rk ) ref2(3,n)
  real ( kind = rk ), dimension(3,4) :: t = reshape ( (/ &
    1.0D+00,  2.0D+00,  3.0D+00, &
    4.0D+00,  1.0D+00,  2.0D+00, &
    2.0D+00,  4.0D+00,  4.0D+00, &
    3.0D+00,  2.0D+00,  5.0D+00 /), (/ 3, 4 /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  For an order 4 tetrahedron,'
  write ( *, '(a)' ) '  PHYSICAL_TO_REFERENCE_TET4 maps a physical point to'
  write ( *, '(a)' ) '    a reference point.'
  write ( *, '(a)' ) '  REFERENCE_TO_PHYSICAL_TET4 maps a reference point to'
  write ( *, '(a)' ) '    a physical point.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    ( R       S       T ) ==> ( X       Y       Z ) ' // &
    '==> ( R2      S2      T2 )'
  write ( *, '(a)' ) ' '

  call reference_tet4_uniform ( n, ref )

  call reference_to_physical_tet4 ( t, n, ref, phy )
  call physical_to_reference_tet4 ( t, n, phy, ref2 )

  do j = 1, n

    write ( *, '(2x,3f8.4,2x,3f8.4,2x,3f8.4)' ) &
      ref(1:3,j), phy(1:3,j), ref2(1:3,j)
  end do

  return
end
