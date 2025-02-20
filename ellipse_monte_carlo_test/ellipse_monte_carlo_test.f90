program main

!*****************************************************************************80
!
!! ELLIPSE_MONTE_CARLO_TEST() tests ELLIPSE_MONTE_CARLO().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 April 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ELLIPSE_MONTE_CARLO_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the ELLIPSE_MONTE_CARLO library.'

  call ellipse_area1_test ( )
  call ellipse_area2_test ( )
  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ELLIPSE_MONTE_CARLO_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine ellipse_area1_test ( )

!*****************************************************************************80
!
!! ELLIPSE_AREA1_TEST tests ELLIPSE_AREA1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 November 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(2,2)
  real ( kind = rk ) area
  real ( kind = rk ) ellipse_area1
  real ( kind = rk ) r

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'ELLIPSE_AREA1_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  ELLIPSE_AREA1 computes the area of an ellipse.'

  r = 10.0D+00

  a = reshape ( (/ 5.0D+00, 1.0D+00, 1.0D+00, 2.0D+00 /), (/ 2, 2 /) )

  area = ellipse_area1 ( a, r )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  R = ', r
  call r8mat_print ( 2, 2, a, '  Matrix A in ellipse definition x*A*x=r^2' )
  write ( *, '(a,g14.6)' ) '  Area = ', area
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'ELLIPSE_AREA1_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'

  return
end
subroutine ellipse_area2_test ( )

!*****************************************************************************80
!
!! ELLIPSE_AREA2_TEST tests ELLIPSE_AREA2.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 November 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) area
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) d
  real ( kind = rk ) ellipse_area2
 
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'ELLIPSE_AREA2_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  ELLIPSE_AREA2 computes the area of an ellipse.'

  a = 5.0D+00
  b = 2.0D+00
  c = 2.0D+00
  d = 10.0D+00

  area = ellipse_area2 ( a, b, c, d )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6,a,g14.6,a,g14.6,a,g14.6)' ) &
    '  Ellipse: ', a, ' * x^2 + ', b, ' * xy + ', c, ' * y^2 = ', d
  write ( *, '(a,g14.6)' ) '  Area = ', area
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'ELLIPSE_AREA2_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'

  return
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 uses ELLIPSE01_SAMPLE with an increasing number of points.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 April 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 2

  real ( kind = rk ), dimension ( m, m ) :: a = reshape ( (/ &
    9.0, 1.0, &
    1.0, 4.0 /), (/ m, m /) )
  real ( kind = rk ) area
  integer e(m)
  integer :: e_test(m,7) = reshape ( (/ &
    0, 0, &
    1, 0, &
    0, 1, &
    2, 0, &
    1, 1, &
    0, 2, &
    3, 0 /), (/ m, 7 /) )
  real ( kind = rk ) ellipse_area1
  integer j
  integer n
  real ( kind = rk ), parameter :: r = 2.0D+00
  real ( kind = rk ) result(7)
  real ( kind = rk ), allocatable :: value(:)
  real ( kind = rk ), allocatable :: x(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Use ELLIPSE01_SAMPLE to estimate integrals'
  write ( *, '(a)' ) '  in the ellipse x'' * A * x <= r^2.'

  area = ellipse_area1 ( a, r )
  write ( *, '(a,g14.6)' ) '  Ellipse area = ', area
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N        1              X               Y  ' // &
    '             X^2               XY             Y^2             X^3'
  write ( *, '(a)' ) ' '

  n = 1

  do while ( n <= 65536 )

    allocate ( value(1:n) )
    allocate ( x(1:m,1:n) )

    call ellipse_sample ( n, a, r, x )

    do j = 1, 7

      e(1:m) = e_test(1:m,j)

      call monomial_value ( m, n, e, x, value )

      result(j) = area * sum ( value(1:n) ) / real ( n, kind = rk )

    end do

    write ( *, '(2x,i8,7(2x,g14.6))' ) n, result(1:7)

    deallocate ( value )
    deallocate ( x )

    n = 2 * n

  end do

  return
end
