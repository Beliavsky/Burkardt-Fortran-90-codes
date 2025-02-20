program main

!*****************************************************************************80
!
!! WEDGE_MONTE_CARLO_TEST() tests WEDGE_MONTE_CARLO().
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
  write ( *, '(a)' ) 'WEDGE_MONTE_CARLO_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the WEDGE_MONTE_CARLO library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WEDGE_MONTE_CARLO_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 uses WEDGE01_SAMPLE with an increasing number of points.
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

  integer e(m)
  integer :: e_test(m,8) = reshape ( (/ &
    0, 0, 0, &
    1, 0, 0, &
    0, 1, 0, &
    0, 0, 1, &
    2, 0, 0, &
    1, 1, 0, &
    0, 0, 2, &
    3, 0, 0 /), (/ m, 8 /) )
  integer j
  integer n
  real ( kind = rk ) result(8)
  real ( kind = rk ), allocatable :: value(:)
  real ( kind = rk ), allocatable :: x(:,:)
  real ( kind = rk ) wedge01_volume

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Use WEDGE01_SAMPLE for a Monte Carlo estimate of an'
  write ( *, '(a)' ) '  integral over the interior of the unit wedge in 3D.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N        1               X               Y ' // &
    '              Z                X^2            XY              Z^2    ' // &
    '        X^3'
  write ( *, '(a)' ) ' '

  n = 1

  do while ( n <= 65536 )

    allocate ( value(1:n) )
    allocate ( x(1:m,1:n) )

    call wedge01_sample ( n, x )

    do j = 1, 8

      e(1:m) = e_test(1:m,j)

      call monomial_value ( m, n, e, x, value )

      result(j) = wedge01_volume ( ) * sum ( value(1:n) ) / real ( n, kind = rk )

    end do

    write ( *, '(2x,i8,8(2x,g14.6))' ) n, result(1:8)

    deallocate ( value )
    deallocate ( x )

    n = 2 * n

  end do

  write ( *, '(a)' ) ' '

  do j = 1, 8

    e(1:m) = e_test(1:m,j)

    call wedge01_integral ( e, result(j) )

  end do

  write ( *, '(2x,a8,8(2x,g14.6))' ) '   Exact', result(1:8)

  return
end
