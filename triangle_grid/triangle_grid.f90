subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine triangle_grid ( n, t, tg )  

!*****************************************************************************80
!
!! TRIANGLE_GRID computes points on a triangular grid.
!
!  Discussion:
!
!    The grid is defined by specifying the coordinates of an enclosing
!    triangle T, and the number of subintervals each side of the triangle
!    should be divided into.
!
!    Choosing N = 10, for instance, breaks each side into 10 subintervals,
!    and produces a grid of ((10+1)*(10+2))/2 = 66 points.
!
!              X
!             9 X
!            8 9 X
!           7 8 9 X
!          6 7 8 9 X
!         5 6 7 8 9 X
!        4 5 6 7 8 9 X
!       3 4 5 6 7 8 9 X
!      2 3 4 5 6 7 8 9 X
!     1 2 3 4 5 6 7 8 9 X
!    0 1 2 3 4 5 6 7 8 9 X
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of subintervals.
!
!    Input, real ( kind = rk ) T(2,3), the coordinates of the points
!    defining the triangle.
!
!    Output, real ( kind = rk ) TG(2,((N+1)*(N+2))/2), the coordinates
!    of the points in the triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  real ( kind = rk ) ir
  integer j
  real ( kind = rk ) jr
  integer k
  real ( kind = rk ) kr
  real ( kind = rk ) nr
  integer p
  real ( kind = rk ) t(2,3)
  real ( kind = rk ) tg(2,((n+1)*(n+2))/2)

  p = 0
  nr = real ( n, kind = rk )

  do i = 0, n
    ir = real ( i, kind = rk )
    do j = 0, n - i
      jr = real ( j, kind = rk )
      k = n - i - j
      kr = real ( k, kind = rk )
      p = p + 1
      tg(1:2,p) = ( ir * t(1:2,1) + jr * t(1:2,2) + kr * t(1:2,3) ) / nr
    end do
  end do

  return
end
