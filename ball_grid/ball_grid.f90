subroutine ball_grid ( n, r, c, ng, bg )

!*****************************************************************************80
!
!! ball_grid() computes grid points inside a ball.
!
!  Discussion:
!
!    The grid is defined by specifying the radius and center of the ball,
!    and the number of subintervals N into which the horizontal radius
!    should be divided.  Thus, a value of N = 2 will result in 5 points
!    along that horizontal line.
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
!  Input:
!
!    integer N, the number of subintervals.
!
!    real ( kind = rk ) R, the radius of the ball.
!
!    real ( kind = rk ) C(3), the coordinates of the center of the ball.
!
!    integer NG, the number of grid points, as determined by
!    BALL_GRID_COUNT.
!
!  Output:
!
!    real ( kind = rk ) BG(3,NG), the grid points inside the ball.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ng

  real ( kind = rk ) bg(3,ng)
  real ( kind = rk ) c(3)
  integer i
  integer j
  integer k
  integer n
  integer p
  real ( kind = rk ) r
  real ( kind = rk ) x
  real ( kind = rk ) y
  real ( kind = rk ) z

  p = 0

  do i = 0, n

    x = c(1) + r * real ( 2 * i, kind = rk ) / real ( 2 * n + 1, kind = rk )
    
    do j = 0, n

      y = c(2) + r * real ( 2 * j, kind = rk ) / real ( 2 * n + 1, kind = rk )

      do k = 0, n

        z = c(3) + r * real ( 2 * k, kind = rk ) / real ( 2 * n + 1, kind = rk )

        if ( r * r < ( x - c(1) )**2 &
                   + ( y - c(2) )**2 &
                   + ( z - c(3) )**2 ) then
          exit
        end if

        p = p + 1
        bg(1,p) = x
        bg(2,p) = y
        bg(3,p) = z

        if ( 0 < i ) then
          p = p + 1
          bg(1,p) = 2.0D+00 * c(1) - x
          bg(2,p) = y
          bg(3,p) = z
        end if

        if ( 0 < j ) then
          p = p + 1
          bg(1,p) = x
          bg(2,p) = 2.0D+00 * c(2) - y
          bg(3,p) = z
        end if

        if ( 0 < k ) then
          p = p + 1
          bg(1,p) = x
          bg(2,p) = y
          bg(3,p) = 2.0D+00 * c(3) - z
        end if

        if ( 0 < i .and. 0 < j ) then
          p = p + 1
          bg(1,p) = 2.0D+00 * c(1) - x
          bg(2,p) = 2.0D+00 * c(2) - y
          bg(3,p) = z
        end if

        if ( 0 < i .and. 0 < k ) then
          p = p + 1
          bg(1,p) = 2.0D+00 * c(1) - x
          bg(2,p) = y
          bg(3,p) = 2.0D+00 * c(3) - z
        end if

        if ( 0 < j .and. 0 < k ) then
          p = p + 1
          bg(1,p) = x
          bg(2,p) = 2.0D+00 * c(2) - y
          bg(3,p) = 2.0D+00 * c(3) - z
        end if

        if ( 0 < i .and. 0 < j .and. 0 < k ) then
          p = p + 1
          bg(1,p) = 2.0D+00 * c(1) - x
          bg(2,p) = 2.0D+00 * c(2) - y
          bg(3,p) = 2.0D+00 * c(3) - z
        end if

      end do
    end do
  end do

  return
end
subroutine ball_grid_count ( n, r, c, ng )

!*****************************************************************************80
!
!! ball_grid_count() counts grid points inside a ball.
!
!  Discussion:
!
!    The grid is defined by specifying the radius and center of the ball,
!    and the number of subintervals N into which the horizontal radius
!    should be divided.  Thus, a value of N = 2 will result in 5 points
!    along that horizontal line.
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
!  Input:
!
!    integer N, the number of subintervals.
!
!    real ( kind = rk ) R, the radius of the ball.
!
!    real ( kind = rk ) C(3), the coordinates of the center of the ball.
!
!  Output:
!
!    integer NG, the number of grid points inside the ball.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) c(3)
  integer i
  integer j
  integer k
  integer n
  integer ng
  real ( kind = rk ) r
  real ( kind = rk ) x
  real ( kind = rk ) y
  real ( kind = rk ) z

  ng = 0

  do i = 0, n

    x = c(1) + r * real ( 2 * i, kind = rk ) / real ( 2 * n + 1, kind = rk )
    
    do j = 0, n

      y = c(2) + r * real ( 2 * j, kind = rk ) / real ( 2 * n + 1, kind = rk )

      do k = 0, n

        z = c(3) + r * real ( 2 * k, kind = rk ) / real ( 2 * n + 1, kind = rk )

        if ( r * r < ( x - c(1) )**2 &
                   + ( y - c(2) )**2 &
                   + ( z - c(3) )**2 ) then
          exit
        end if

        ng = ng + 1

        if ( 0 < i ) then
          ng = ng + 1
        end if

        if ( 0 < j ) then
          ng = ng + 1
        end if

        if ( 0 < k ) then
          ng = ng + 1
        end if

        if ( 0 < i .and. 0 < j ) then
          ng = ng + 1
        end if

        if ( 0 < i .and. 0 < k ) then
          ng = ng + 1
        end if

        if ( 0 < j .and. 0 < k ) then
          ng = ng + 1
        end if

        if ( 0 < i .and. 0 < j .and. 0 < k ) then
          ng = ng + 1
        end if

      end do
    end do
  end do

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! get_unit() returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
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
!  Output:
!
!    integer IUNIT, the free unit number.
!
  implicit none

  integer i
  integer ios
  integer iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
subroutine r83vec_print_part ( n, a, max_print, title )

!*****************************************************************************80
!
!! r83vec_print_part() prints "part" of an R83VEC.
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vector, is no more than MAX_PRINT, then
!    the entire vector is printed, one entry per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
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
!  Input:
!
!    integer N, the number of entries of the vector.
!
!    real ( kind = rk ) A(3,N), the vector to be printed.
!
!    integer MAX_PRINT, the maximum number of lines
!    to print.
!
!    character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(3,n)
  integer i
  integer max_print
  character ( len = * ) title

  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  if ( n <= max_print ) then

    do i = 1, n
      write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6,2x,g14.6)' ) i, ':', a(1:3,i)
    end do

  else if ( 3 <= max_print ) then

    do i = 1, max_print - 2
      write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6,2x,g14.6)' ) i, ':', a(1:3,i)
    end do
    write ( *, '(a)' ) &
      '  ........  ..............  ..............  ..............'
    i = n
    write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6,2x,g14.6)' ) i, ':', a(1:3,i)

  else

    do i = 1, max_print - 1
      write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6,2x,g14.6)' ) i, ':', a(1:3,i)
    end do
    i = max_print
    write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6,2x,g14.6,2x,a)' ) i, ':', a(1:3,i), &
      '...more entries...'

  end if

  return
end
subroutine r8mat_write ( output_filename, m, n, table )

!*****************************************************************************80
!
!! r8mat_write() writes an R8MAT file.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
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
!  Input:
!
!    character ( len = * ) OUTPUT_FILENAME, the output file name.
!
!    integer M, the spatial dimension.
!
!    integer N, the number of points.
!
!    real ( kind = rk ) TABLE(M,N), the data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer j
  character ( len = * ) output_filename
  integer output_status
  integer output_unit
  character ( len = 30 ) string
  real ( kind = rk ) table(m,n)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_filename ) // '" on unit ', output_unit
    output_unit = -1
    stop
  end if
!
!  Create a format string.
!
!  For less precision in the output file, try:
!
!                                            '(', m, 'g', 14, '.', 6, ')'
!
  if ( 0 < m .and. 0 < n ) then

    write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 24, '.', 16, ')'
!
!  Write the data.
!
    do j = 1, n
      write ( output_unit, string ) table(1:m,j)
    end do

  end if
!
!  Close the file.
!
  close ( unit = output_unit )

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
!    01 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

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
 
