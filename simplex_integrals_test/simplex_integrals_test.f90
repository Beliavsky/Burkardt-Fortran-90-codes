program main

!*****************************************************************************80
!
!! simplex_integrals_test() tests simplex_integrals().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 January 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'simplex_integrals_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test SIMPLEX_INTEGRALS().'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'simplex_integrals_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! test01() compares exact and estimated integrals in 3D.
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 3
  integer, parameter :: n = 4192

  integer e(m)
  real ( kind = rk ) error
  real ( kind = rk ) exact
  real ( kind = rk ) result
  real ( kind = rk ) simplex01_volume
  integer test
  integer, parameter :: test_num = 20
  real ( kind = rk ) value(n)
  real ( kind = rk ), allocatable :: x(:,:)

  allocate ( x(1:m,1:n) )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Estimate monomial integrals using Monte Carlo'
  write ( *, '(a)' ) '  over the interior of the unit simplex in M dimensions.'
!
!  Get sample points.
!
  call simplex01_sample ( m, n, x )

  write ( *, '(a)' ) ''
  write ( *, '(a,i6)' ) '  Number of sample points used is ', n
!
!  Randomly choose exponents.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  We randomly choose the exponents.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Ex  Ey  Ez     MC-Estimate      Exact           Error'
  write ( *, '(a)' ) ''

  do test = 1, test_num

    call i4vec_uniform_ab ( m, 0, 4, e )

    call monomial_value ( m, n, e, x, value )

    result = simplex01_volume ( m ) * sum ( value(1:n) ) &
      / real ( n, kind = rk )
    call simplex01_monomial_integral ( m, e, exact )
    error = abs ( result - exact )

    write ( *, '(2x,i2,2x,i2,2x,i2,2x,g14.6,2x,g14.6,2x,e10.2)' ) &
      e(1:m), result, exact, error

  end do

  deallocate ( x )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! test02() compares exact and estimated integrals in 6D.
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 6
  integer, parameter :: n = 4192

  integer e(m)
  real ( kind = rk ) error
  real ( kind = rk ) exact
  real ( kind = rk ) result
  real ( kind = rk ) simplex01_volume
  integer test
  integer, parameter :: test_num = 20
  real ( kind = rk ) value(n)
  real ( kind = rk ), allocatable :: x(:,:)

  allocate ( x(1:m,1:n) )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'test02()'
  write ( *, '(a)' ) '  Estimate monomial integrals using Monte Carlo'
  write ( *, '(a)' ) '  over the interior of the unit simplex in M dimensions.'
!
!  Get sample points.
!
  call simplex01_sample ( m, n, x )

  write ( *, '(a)' ) ''
  write ( *, '(a,i6)' ) '  Number of sample points used is ', n
!
!  Randomly choose exponents.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  We randomly choose the exponents.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  E1  E2  E3  E4  E5  E6     MC-Estimate      Exact           Error'
  write ( *, '(a)' ) ''

  do test = 1, test_num

    call i4vec_uniform_ab ( m, 0, 3, e )

    call monomial_value ( m, n, e, x, value )

    result = simplex01_volume ( m ) * sum ( value(1:n) ) &
      / real ( n, kind = rk )
    call simplex01_monomial_integral ( m, e, exact )
    error = abs ( result - exact )

    write ( *, '(6(2x,i2),2x,g14.6,2x,g14.6,2x,e10.2)' ) &
      e(1:m), result, exact, error

  end do

  deallocate ( x )

  return
end
