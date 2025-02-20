program main

!*****************************************************************************80
!
!! circle_integrals_test() tests circle_integrals().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 January 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CIRCLE_INTEGRALS_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the CIRCLE_INTEGRALS library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CIRCLE_INTEGRALS_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 uses CIRCLE01_SAMPLE to compare exact and estimated monomial integrals.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 January 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 2
  integer, parameter :: n = 4192

  real ( kind = rk ) circle01_length
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
  write ( *, '(a)' ) '  Use CIRCLE01_SAMPLE to compare exact and'
  write ( *, '(a)' ) '  estimated integrals along the circumference '
  write ( *, '(a)' ) '  of the unit circle in 2D.'
!
!  Get sample points.
!
  call circle01_sample ( n, x )

  write ( *, '(a)' ) ''
  write ( *, '(a,i6)' ) '  Number of sample points used is ', n
!
!  Randomly choose X, Y exponents.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  If any exponent is odd, the integral is zero.'
  write ( *, '(a)' ) '  We restrict this test to randomly chosen even exponents.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Ex  Ey     MC-Estimate           Exact      Error'
  write ( *, '(a)' ) ''

  do test = 1, test_num

    call i4vec_uniform_ab ( m, 0, 5, e )

    e(1:m) = e(1:m) * 2

    call monomial_value ( m, n, e, x, value )

    result = circle01_length ( ) * sum ( value(1:n) ) &
      / real ( n, kind = rk )
    call circle01_monomial_integral ( e, exact )
    error = abs ( result - exact )

    write ( *, '(2x,i2,2x,i2,2x,g14.6,2x,g14.6,2x,e10.2)' ) &
      e(1:m), result, exact, error

  end do

  return
end
 
