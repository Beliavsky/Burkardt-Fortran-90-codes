program main

!*****************************************************************************80
!
!! MAIN is the main program for TETRAHEDRON_INTEGRALS_TEST.
!
!  Discussion:
!
!    TETRAHEDRON_INTEGRALS_TEST tests the TETRAHEDRON_INTEGRALS library.
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

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TETRAHEDRON_INTEGRALS_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TETRAHEDRON_INTEGRALS library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TETRAHEDRON_INTEGRALS_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 uses TETRAHEDRON01_SAMPLE to compare exact and estimated integrals.
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

  integer, parameter :: m = 3
  integer, parameter :: n = 4192

  integer e(m)
  real ( kind = rk ) error
  real ( kind = rk ) exact
  integer i
  integer j
  integer k
  real ( kind = rk ) result
  integer seed
  integer, parameter :: test_num = 20
  real ( kind = rk ) tetrahedron01_volume
  real ( kind = rk ) value(n)
  real ( kind = rk ) x(m,n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Estimate monomial integrals using Monte Carlo'
  write ( *, '(a)' ) '  over the interior of the unit tetrahedron in 3D.'
!
!  Get sample points.
!
  seed = 123456789
  call tetrahedron01_sample ( n, seed, x )

  write ( *, '(a)' ) ''
  write ( *, '(a,i6)' ) '  Number of sample points used is ', n
!
!  Run through the exponents.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Ex  Ey  Ez     MC-Estimate      Exact           Error'
  write ( *, '(a)' ) ''

  do i = 0, 3
    e(1) = i
    do j = 0, 3
      e(2) = j
      do k = 0, 3
        e(3) = k

        call monomial_value ( m, n, e, x, value )

        result = tetrahedron01_volume ( ) * sum ( value(1:n) ) &
          / real ( n, kind = rk )
        call tetrahedron01_monomial_integral ( e, exact )
        error = abs ( result - exact )

        write ( *, '(2x,i2,2x,i2,2x,i2,2x,g14.6,2x,g14.6,2x,e10.2)' ) &
          e(1:m), result, exact, error

      end do
    end do
  end do

  return
end
