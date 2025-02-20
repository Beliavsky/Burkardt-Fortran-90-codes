program main

!*****************************************************************************80
!
!! polygon_monte_carlo_test() tests polygon_monte_carlo().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 January 2025
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nv1 = 4

  real ( kind = rk ), dimension ( 2, nv1 ) :: v1 = reshape ( (/ &
    -1.0, -1.0, &
     1.0, -1.0, &
     1.0,  1.0, &
    -1.0,  1.0 /), (/ 2, nv1 /) )

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'polygon_monte_carlo_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test polygon_monte_carlo().'

  call test01 ( nv1, v1 )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'polygon_monte_carlo_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( nv, v )

!*****************************************************************************80
!
!! TEST01 estimates integrals over a polygon in 2D.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 May 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nv

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
  real ( kind = rk ) polygon_area
  real ( kind = rk ) result(7)
  real ( kind = rk ) v(2,nv)
  real ( kind = rk ), allocatable :: value(:)
  real ( kind = rk ), allocatable :: x(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Use POLYGON_SAMPLE to estimate integrals '
  write ( *, '(a)' ) '  over the interior of a polygon in 2D.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '         N' // &
    '        1' // &
    '              X^2 ' // &
    '             Y^2' // &
    '             X^4' // &
    '           X^2Y^2' // &
    '             Y^4' // &
    '           X^6'
  write ( *, '(a)' ) ' '

  n = 1

  do while ( n <= 65536 )

    allocate ( value(1:n) )
    allocate ( x(1:2,1:n) )

    call polygon_sample ( nv, v, n, x )

    do j = 1, 7

      e(1:2) = e_test(1:2,j)

      call monomial_value ( 2, n, e, x, value )

      result(j) = polygon_area ( nv, v ) * sum ( value(1:n) ) &
        / real ( n, kind = rk )

    end do

    write ( *, '(2x,i8,7(2x,g14.6))' ) n, result(1:7)

    deallocate ( value )
    deallocate ( x )

    n = 2 * n

  end do

  write ( *, '(a)' ) ' '

  do j = 1, 7

    e(1:2) = e_test(1:2,j)

    call polygon_monomial_integral ( nv, v, e, result(j) )

  end do

  write ( *, '(2x,a8,7(2x,g14.6))' ) '   Exact', result(1:7)

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! timestamp() prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = 8 ) ampm
  integer d
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  integer values(8)
  integer y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end

