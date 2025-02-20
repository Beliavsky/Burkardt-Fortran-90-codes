program main

!*****************************************************************************80
!
!! sphere_integrals_test() tests sphere_integrals().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 June 2024
!
!  Author:
!
!    John Burkardt
!
  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'sphere_integrals_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test SPHERE_INTEGRALS().'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'SPHERE_INTEGRALS_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 uses SPHERE01_SAMPLE to estimate monomial integrands.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 January 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 3
  integer, parameter :: n = 4192

  integer e(m)
  real ( kind = rk ) error
  real ( kind = rk ) exact
  real ( kind = rk ) result
  real ( kind = rk ) sphere01_area
  integer test
  integer, parameter :: test_num = 20
  real ( kind = rk ) value(n)
  real ( kind = rk ) x(m,n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Estimate monomial integrals using Monte Carlo'
  write ( *, '(a)' ) '  over the surface of the unit sphere in 3D.'
!
!  Get sample points.
!
  call sphere01_sample ( n, x )

  write ( *, '(a)' ) ''
  write ( *, '(a,i6)' ) '  Number of sample points used is ', n
!
!  Randomly choose X,Y,Z exponents between (0,0,0) and (9,9,9).
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

    result = sphere01_area ( ) * sum ( value(1:n) ) &
      / real ( n, kind = rk )
    call sphere01_monomial_integral ( e, exact )
    error = abs ( result - exact )

    write ( *, '(2x,i2,2x,i2,2x,i2,2x,g14.6,2x,g14.6,2x,e10.2)' ) &
      e(1:3), result, exact, error

  end do

  return
end

