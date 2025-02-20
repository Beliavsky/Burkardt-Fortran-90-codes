program main

!*****************************************************************************80
!
!! ball_monte_carlo_test() tests ball_monte_carlo().
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
  write ( *, '(a)' ) 'ball_monte_carlo_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test ball_monte_carlo().'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ball_monte_carlo_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! test01() uses ball01_sample() with an increasing number of points.
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

  real ( kind = rk ) ball01_volume
  integer e(3)
  integer :: e_test(3,7) = reshape ( (/ &
    0, 0, 0, &
    2, 0, 0, &
    0, 2, 0, &
    0, 0, 2, &
    4, 0, 0, &
    2, 2, 0, &
    0, 0, 4 /), (/ 3, 7 /) )
  integer j
  integer n
  real ( kind = rk ) result(7)
  real ( kind = rk ), allocatable :: value(:)
  real ( kind = rk ), allocatable :: x(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test01():'
  write ( *, '(a)' ) '  Use ball01_sample() to estimate integrals over the interior'
  write ( *, '(a)' ) '  of the unit ball using the Monte Carlo method.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N        1              X^2             Y^2' // &
    '             Z^2             X^4           X^2Y^2           Z^4'
  write ( *, '(a)' ) ' '

  n = 1

  do while ( n <= 65536 )

    allocate ( value(1:n) )
    allocate ( x(1:3,1:n) )

    call ball01_sample ( n, x )

    do j = 1, 7

      e(1:3) = e_test(1:3,j)

      call monomial_value ( 3, n, e, x, value )

      result(j) = ball01_volume ( ) * sum ( value(1:n) ) &
        / real ( n, kind = rk )

    end do

    write ( *, '(2x,i8,7(2x,g14.6))' ) n, result(1:7)

    deallocate ( value )
    deallocate ( x )

    n = 2 * n

  end do

  write ( *, '(a)' ) ' '

  do j = 1, 7

    e(1:3) = e_test(1:3,j)

    call ball01_monomial_integral ( e, result(j) )

  end do

  write ( *, '(2x,a8,7(2x,g14.6))' ) '   Exact', result(1:7)

  return
end
