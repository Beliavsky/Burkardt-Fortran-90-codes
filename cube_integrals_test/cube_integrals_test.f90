program main

!*****************************************************************************80
!
!! cube_integrals_test() tests cube_integrals().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 January 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CUBE_INTEGRALS_TEST():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test CUBE_INTEGRALS().'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CUBE_INTEGRALS_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 compares exact and estimated integrals over the unit cube in 3D.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 January 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 3
  integer, parameter :: n = 4192

  real ( kind = rk ) cube01_volume
  integer e(m)
  real ( kind = rk ) error
  real ( kind = rk ) exact
  real ( kind = rk ) result
  integer test
  integer, parameter :: test_num = 20
  real ( kind = rk ) value(n)
  real ( kind = rk ) x(m,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Compare exact and estimated integrals'
  write ( *, '(a)' ) '  over the interior of the unit cube in 3D.'
!
!  Get sample points.
!
  call cube01_sample ( n, x )

  write ( *, '(a)' ) ''
  write ( *, '(a,i6)' ) '  Number of sample points used is ', n
!
!  Randomly choose exponents.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Ex  Ey  Ez     MC-Estimate           Exact      Error'
  write ( *, '(a)' ) ''

  do test = 1, test_num

    call i4vec_uniform_ab ( m, 0, 7, e )

    call monomial_value ( m, n, e, x, value )

    result = cube01_volume ( ) * sum ( value(1:n) ) &
      / real ( n, kind = rk )
    call cube01_monomial_integral ( e, exact )
    error = abs ( result - exact )

    write ( *, '(3(2x,i2),2x,g14.6,2x,g14.6,2x,e10.2)' ) &
      e(1:m), result, exact, error

  end do

  return
end
