program main

!*****************************************************************************80
!
!! MONOMIAL_VALUE_TEST() tests MONOMIAL_VALUE().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 August 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MONOMIAL_VALUE_TEST():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test MONOMIAL_VALUE().'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MONOMIAL_VALUE_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests MONOMIAL_VALUE on sets of data in various dimensions.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 August 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer, allocatable :: e(:)
  integer e_max
  integer e_min
  integer i
  integer j
  real ( kind = rk ), allocatable :: v(:)
  real ( kind = rk ), allocatable :: x(:,:)
  real ( kind = rk ) x_max
  real ( kind = rk ) x_min

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  monomial_value() evaluates monomials in'
  write ( *, '(a)' ) '  dimensions 1 through 3.' 

  e_min = -3
  e_max = 6
  n = 5
  x_min = -2.0D+00
  x_max = +10.0D+00

  do m = 1, 3

    write ( *, '(a)' ) ''
    write ( *, '(a,i1)' ) '  Spatial dimension M = ', m

    allocate ( e(1:m) )
    allocate ( x(1:m,1:n) )
    allocate ( v(1:n) )

    call i4vec_uniform_ab ( m, e_min, e_max, e )
    call i4vec_transpose_print ( m, e, '  Exponents:' )
    call r8mat_uniform_ab ( m, n, x_min, x_max, x )
!
!  To make checking easier, make the X values integers.
!
    call r8mat_nint ( m, n, x )
    call monomial_value ( m, n, e, x, v )

    write ( *, '(a)' ) ''
    write ( *, '(a)', advance = 'no' ) '   V(X)         '
    do i = 1, m
      write ( *, '(a,i1,a)', advance = 'no' ) '      X(', i, ')'
    end do
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) ''
    do j = 1, n
      write ( *, '(g14.6,2x,3f10.4)' ) v(j), x(1:m,j)
    end do

    deallocate ( e )
    deallocate ( x )
    deallocate ( v )

  end do

  return
end
