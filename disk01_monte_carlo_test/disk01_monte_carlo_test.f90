program main

!*****************************************************************************80
!
!! disk01_monte_carlo_test() tests disk01_monte_carlo().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 January 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DISK01_MONTE_CARLO_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the DISK01_MONTE_CARLO library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DISK01_MONTE_CARLO_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 uses DISK01_SAMPLE with an increasing number of points.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 January 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) disk01_area
  integer e(2)
  integer :: e_test(2,7) = reshape ( (/ &
    0, 0, &
    2, 0, &
    0, 2, &
    4, 0, &
    2, 2, &
    0, 4, &
    6, 0 /), (/ 2, 7 /) )
  integer j
  integer n
  real ( kind = rk ) result(7)
  real ( kind = rk ), allocatable :: value(:)
  real ( kind = rk ), allocatable :: x(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Use DISK01_SAMPLE to estimate integrals in the unit disk.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N        1              X^2             Y^2' // &
    '             X^4             X^2Y^2           Y^4             X^6'
  write ( *, '(a)' ) ' '

  n = 1

  do while ( n <= 65536 )

    allocate ( value(1:n) )
    allocate ( x(1:2,1:n) )

    call disk01_sample ( n, x )

    do j = 1, 7

      e(1:2) = e_test(1:2,j)

      call monomial_value ( 2, n, e, x, value )

      result(j) = disk01_area ( ) * sum ( value(1:n) ) / real ( n, kind = rk )

    end do

    write ( *, '(2x,i8,7(2x,g14.6))' ) n, result(1:7)

    deallocate ( value )
    deallocate ( x )

    n = 2 * n

  end do

  write ( *, '(a)' ) ' '

  do j = 1, 7

    e(1:2) = e_test(1:2,j)

    call disk01_monomial_integral ( e, result(j) )

  end do

  write ( *, '(2x,a8,7(2x,g14.6))' ) '   Exact', result(1:7)

  return
end
