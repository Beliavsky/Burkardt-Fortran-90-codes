program main

!*****************************************************************************80
!
!! HYPERSPHERE_INTEGRALS_TEST() tests HYPERSPHERE_INTEGRALS().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 January 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HYPERSPHERE_INTEGRALS_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the HYPERSPHERE_INTEGRALS library.'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HYPERSPHERE_INTEGRALS_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 uses HYPERSPHERE01_SAMPLE to estimate monomial integrands in 3D.
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
  real ( kind = rk ) hypersphere01_area
  real ( kind = rk ) result
  integer test
  integer, parameter :: test_num = 20
  real ( kind = rk ) value(n)
  real ( kind = rk ) x(m,n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Estimate monomial integrals using Monte Carlo'
  write ( *, '(a)' ) '  over the surface of the unit hypersphere in 3D.'
!
!  Get sample points.
!
  call hypersphere01_sample ( m, n, x )

  write ( *, '(a)' ) ''
  write ( *, '(a,i6)' ) '  Number of sample points used is ', n
!
!  Randomly choose X,Y,Z exponents between (0,0,0) and (8,8,8).
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

    result = hypersphere01_area ( m ) * sum ( value(1:n) ) &
      / real ( n, kind = rk )
    call hypersphere01_monomial_integral ( m, e, exact )
    error = abs ( result - exact )

    write ( *, '(2x,i2,2x,i2,2x,i2,2x,g14.6,2x,g14.6,2x,e10.2)' ) &
      e(1:m), result, exact, error

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 uses HYPERSPHERE01_SAMPLE to estimate monomial integrands in 6D.
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

  integer, parameter :: m = 6
  integer, parameter :: n = 4192

  integer e(m)
  real ( kind = rk ) error
  real ( kind = rk ) exact
  real ( kind = rk ) hypersphere01_area
  real ( kind = rk ) result
  integer test
  integer, parameter :: test_num = 20
  real ( kind = rk ) value(n)
  real ( kind = rk ) x(m,n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Estimate monomial integrals using Monte Carlo'
  write ( *, '(a)' ) '  over the surface of the unit hypersphere in 3D.'
!
!  Get sample points.
!
  call hypersphere01_sample ( m, n, x )

  write ( *, '(a)' ) ''
  write ( *, '(a,i6)' ) '  Number of sample points used is ', n
!
!  Randomly choose X,Y,Z exponents between (0,0,0) and (6,6,6).
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  If any exponent is odd, the integral is zero.'
  write ( *, '(a)' ) '  We will restrict this test to randomly chosen even exponents.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  E1  E2  E3  E4  E5  E6     MC-Estimate           Exact      Error'
  write ( *, '(a)' ) ''

  do test = 1, test_num

    call i4vec_uniform_ab ( m, 0, 3, e )

    e(1:m) = e(1:m) * 2

    call monomial_value ( m, n, e, x, value )

    result = hypersphere01_area ( m ) * sum ( value(1:n) ) &
      / real ( n, kind = rk )
    call hypersphere01_monomial_integral ( m, e, exact )
    error = abs ( result - exact )

    write ( *, '(6(2x,i2),2x,g14.6,2x,g14.6,2x,e10.2)' ) &
      e(1:m), result, exact, error

  end do

  return
end

