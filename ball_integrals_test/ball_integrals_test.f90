program main

!*****************************************************************************80
!
!! ball_integrals_test() tests ball_integrals().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ball_integrals_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test ball_integrals().'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ball_integrals_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! test01() uses ball01_sample() to compare estimated and exact integrals.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 3
  integer, parameter :: n = 4192

  real ( kind = rk ) ball01_volume
  integer e(m)
  real ( kind = rk ) error
  real ( kind = rk ) exact
  real ( kind = rk ) result
  integer test
  integer, parameter :: test_num = 20
  real ( kind = rk ) value(n)
  real ( kind = rk ) x(m,n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'test01():'
  write ( *, '(a)' ) '  Use ball01_sample() to estimate monomial integrals'
  write ( *, '(a)' ) '  over the interior of the unit ball in 3D.'
!
!  Get sample points.
!
  call ball01_sample ( n, x )

  write ( *, '(a)' ) ''
  write ( *, '(a,i6)' ) '  Number of sample points used is ', n
!
!  Randomly choose X,Y exponents between 0 and 8.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  If any exponent is odd, the integral is zero.'
  write ( *, '(a)' ) '  We will restrict this test to randomly chosen even exponents.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Ex  Ey  Ez     MC-Estimate           Exact      Error'
  write ( *, '(a)' ) ''

  do test = 1, test_num

    call i4vec_uniform_ab ( m, 0, 4, e )

    e(1:m) = e(1:m) * 2

    call monomial_value ( m, n, e, x, value )

    result = ball01_volume ( ) * sum ( value(1:n) ) &
      / real ( n, kind = rk )
    call ball01_monomial_integral ( e, exact )
    error = abs ( result - exact )

    write ( *, '(2x,i2,2x,i2,2x,i2,2x,g14.6,2x,g14.6,2x,e10.2)' ) &
      e(1:m), result, exact, error

  end do

  return
end
