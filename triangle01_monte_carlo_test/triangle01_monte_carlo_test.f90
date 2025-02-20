program main

!*****************************************************************************80
!
!! MAIN is the main program for TRIANGLE01_MONTE_CARLO_TEST.
!
!  Discussion:
!
!    TRIANGLE01_MONTE_CARLO_TEST tests the TRIANGLE01_MONTE_CARLO library.
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
  write ( *, '(a)' ) 'TRIANGLE01_MONTE_CARLO_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TRIANGLE01_MONTE_CARLO library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIANGLE01_MONTE_CARLO_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 uses TRIANGLE01_SAMPLE with an increasing number of points.
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

  integer e(m)
  integer :: e_test(m,7) = reshape ( (/ &
    0, 0, &
    1, 0, &
    0, 1, &
    2, 0, &
    1, 1, &
    0, 2, &
    3, 0 /), (/ m, 7 /) )
  integer j
  integer n
  real ( kind = rk ) result(7)
  integer seed
  real ( kind = rk ) triangle01_area
  real ( kind = rk ), allocatable :: value(:)
  real ( kind = rk ), allocatable :: x(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Use TRIANGLE01_SAMPLE for a Monte Carlo estimate of an'
  write ( *, '(a)' ) '  integral over the interior of the unit triangle in 2D.'

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N        1               X               Y ' // &
    '             X^2               XY             Y^2             X^3'
  write ( *, '(a)' ) ' '

  n = 1

  do while ( n <= 65536 )

    allocate ( value(1:n) )
    allocate ( x(1:m,1:n) )

    call triangle01_sample ( n, seed, x )

    do j = 1, 7

      e(1:m) = e_test(1:m,j)

      call monomial_value ( m, n, e, x, value )

      result(j) = triangle01_area ( ) * sum ( value(1:n) ) / real ( n, kind = rk )

    end do

    write ( *, '(2x,i8,7(2x,g14.6))' ) n, result(1:7)

    deallocate ( value )
    deallocate ( x )

    n = 2 * n

  end do

  write ( *, '(a)' ) ' '

  do j = 1, 7

    e(1:m) = e_test(1:m,j)

    call triangle01_monomial_integral ( e, result(j) )

  end do

  write ( *, '(2x,a8,7(2x,g14.6))' ) '   Exact', result(1:7)

  return
end
