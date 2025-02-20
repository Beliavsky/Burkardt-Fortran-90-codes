subroutine line_grid ( n, a, b, c, x )

!*****************************************************************************80
!
!! LINE_GRID: grid points over the interior of a line segment in 1D.
!
!  Discussion:
!
!    In 1D, a grid is to be created, using N points.
!
!    Over the interval [A,B], we have 5 choices for grid centering:
!      1: 0,   1/3, 2/3, 1
!      2: 1/5, 2/5, 3/5, 4/5
!      3: 0,   1/4, 2/4, 3/4
!      4: 1/4, 2/4, 3/4, 1
!      5: 1/8, 3/8, 5/8, 7/8
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    31 August 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = rk ) A, B, the endpoints.
!
!    Input, integer ( kind = 4 ) C, the grid centering.
!    1 <= C <= 5.
!
!    Output, real ( kind = rk ) X(N), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ( kind = 4 ) n

  real ( kind = rk ) a
  real ( kind = rk ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) j
  real ( kind = rk ) x(n)

  do j = 1, n

    if ( c == 1 ) then

      if ( n == 1 ) then
        x(j) = 0.5D+00 * ( a + b )
      else
        x(j) = (   real ( n - j,     kind = rk ) * a   &
                 + real (     j - 1, kind = rk ) * b ) & 
                 / real ( n     - 1, kind = rk )
      end if
    else if ( c == 2 ) then
      x(j) = (   real ( n - j + 1, kind = rk ) * a   &
               + real (     j,     kind = rk ) * b ) & 
               / real ( n     + 1, kind = rk )
    else if ( c == 3 ) then
      x(j) = (   real ( n - j + 1, kind = rk ) * a   &
               + real (     j - 1, kind = rk ) * b ) & 
               / real ( n,         kind = rk )
    else if ( c == 4 ) then
      x(j) = (   real ( n - j, kind = rk ) * a   &
               + real (     j, kind = rk ) * b ) & 
               / real ( n,     kind = rk )
    else if ( c == 5 ) then
      x(j) = (   real ( 2 * n - 2 * j + 1, kind = rk ) * a   &
               + real (         2 * j - 1, kind = rk ) * b ) & 
               / real ( 2 * n,             kind = rk )
    end if

  end do

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = rk ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ( kind = 4 ) n

  real ( kind = rk ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
  end do

  return
end
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
!  Parameters:
!
!    None
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

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
