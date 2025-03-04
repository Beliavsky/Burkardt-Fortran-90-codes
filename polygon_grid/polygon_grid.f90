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
!    26 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IUNIT, the free unit number.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer ios
  integer iunit
  logical ( kind = 4 ) lopen

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
subroutine i4_fake_use ( n )

!*****************************************************************************80
!
!! i4_fake_use() pretends to use a variable.
!
!  Discussion:
!
!    Some compilers will issue a warning if a variable is unused.
!    Sometimes there's a good reason to include a variable in a program,
!    but not to use it.  Calling this function with that variable as
!    the argument will shut the compiler up.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 April 2020
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the variable to be "used".
!
  implicit none

  integer n

  if ( n /= n ) then
    write ( *, '(a)' ) '  i4_fake_use: variable is NAN.'
  end if

  return
end
subroutine polygon_grid_count ( n, nv, ng )

!*****************************************************************************80
!
!! POLYGON_GRID_COUNT counts the grid points inside a polygon.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 May 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of subintervals on a side.
!
!    Input, integer NV, the number of vertices.
!    3 <= NV.
!
!    Output, integer NG, the number of grid points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer ng
  integer nv

  ng = 1 + nv * ( n * ( n + 1 ) ) / 2

  return
end
subroutine polygon_grid_display ( n, nv, v, ng, xg, prefix )

!*****************************************************************************80
!
!! POLYGON_GRID_DISPLAY displays grid points inside a polygon.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 May 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of subintervals.
!
!    Input, integer NV, the number of vertices in the polygon.
!
!    Input, real ( kind = rk ) V(2,NV), the coordinates of the vertices.
!
!    Input, integer NG, the number of grid points.
!
!    Input, real ( kind = rk ) XG(2,NG), the grid points.
!
!    Input, character ( len = * ) PREFIX, a string used to name the files.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ng
  integer nv

  integer command_unit
  character ( len = 255 ) command_filename
  integer grid_unit
  character ( len = 255 ) grid_filename
  integer j
  integer n
  character ( len = 255 ) plot_filename
  character ( len = * ) prefix
  real ( kind = rk ) v(2,nv)
  real ( kind = rk ) vc(2)
  integer vertex_unit
  character ( len = 255 ) vertex_filename
  real ( kind = rk ) xg(2,ng)

  call i4_fake_use ( n )
!
!  Write the vertex file.
!
  vc(1) = sum ( v(1,1:nv) ) / real ( nv, kind = rk )
  vc(2) = sum ( v(2,1:nv) ) / real ( nv, kind = rk )

  call get_unit ( vertex_unit )
  vertex_filename = trim ( prefix ) // '_vertex.txt'
  open ( unit = vertex_unit, file = vertex_filename, status = 'replace' )
  do j = 1, nv
    write ( vertex_unit, '(2x,g14.6,2x,g14.6)' ) v(1,j), v(2,j)
  end do
  write ( vertex_unit, '(2x,g14.6,2x,g14.6)' ) v(1,1), v(2,1)
  do j = 1, nv
    write ( vertex_unit, '(a)' ) ''
    write ( vertex_unit, '(2x,g14.6,2x,g14.6)' ) v(1,j), v(2,j)
    write ( vertex_unit, '(2x,g14.6,2x,g14.6)' ) vc(1), vc(2)
  end do
  close ( unit = vertex_unit )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Created vertex file "' // trim ( vertex_filename ) // '".'
!
!  Write the gridpoint file.
!
  call get_unit ( grid_unit )
  grid_filename = trim ( prefix ) // '_grid.txt'
  open ( unit = grid_unit, file = grid_filename, status = 'replace' )
  do j = 1, ng
    write ( grid_unit, '(2x,g14.6,2x,g14.6)' ) xg(1,j), xg(2,j)
  end do
  close ( unit = grid_unit )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Created grid file "' // trim ( grid_filename ) // '".'
!
!  Write the command file.
!
  plot_filename = trim ( prefix ) // '.png'
  call get_unit ( command_unit )
  command_filename = trim ( prefix ) // '_commands.txt'
  open ( unit = command_unit, file = command_filename, status = 'replace' )
  write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) 'set output "' // trim ( plot_filename ) // '"'
  write ( command_unit, '(a)' ) 'set xlabel "<--- X --->"'
  write ( command_unit, '(a)' ) 'set ylabel "<--- Y --->"'
  write ( command_unit, '(a)' ) 'set title "' // trim ( prefix ) // '"'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set key off'
  write ( command_unit, '(a)' ) 'set size ratio -1'
  write ( command_unit, '(a)' ) 'set style data lines'

  write ( command_unit, '(a)' ) 'plot "' // trim ( grid_filename ) // &
    '" using 1:2 with points lt 3 pt 3,\'
  write ( command_unit, '(a)' ) '    "' // trim ( vertex_filename ) // &
    '" using 1:2 lw 3 linecolor rgb "black"'
  write ( command_unit, '(a)' ) 'quit'
  close ( unit = command_unit )

  write ( *, '(a)' ) &
    '  Created command file "' // trim ( command_filename ) // '".'

  return
end
subroutine polygon_grid_points ( n, nv, v, ng, xg )

!*****************************************************************************80
!
!! POLYGON_GRID_POINTS computes points on a polygonal grid.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 May 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of subintervals.
!
!    Input, integer NV, the number of vertices in the polygon.
!
!    Input, real ( kind = rk ) V(2,NV), the coordinates of the vertices.
!
!    Input, integer NG, the number of grid points.
!
!    Output, real ( kind = rk ) XG(2,NG), the coordinates of the grid points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ng
  integer nv

  integer i
  integer j
  integer k
  integer l
  integer lp1
  integer n
  integer p
  real ( kind = rk ) v(2,nv)
  real ( kind = rk ) vc(2)
  real ( kind = rk ) xg(2,ng)
!
!  Determine the centroid, and use it as the first grid point.
!
  p = 1
  vc(1) = sum ( v(1,1:nv) ) / real ( nv, kind = rk )
  vc(2) = sum ( v(2,1:nv) ) / real ( nv, kind = rk )
  xg(1:2,p) = vc(1:2)
!
!  Consider each triangle formed by two consecutive vertices and the centroid,
!  but skip the first line of points.
!
  do l = 1, nv
    lp1 = mod ( l, nv ) + 1
    do i = 1, n
      do j = 0, n - i
        k = n - i - j
        p = p + 1
        xg(1:2,p) = ( real ( i, kind = rk ) * v(1:2,l)   &
                    + real ( j, kind = rk ) * v(1:2,lp1) & 
                    + real ( k, kind = rk ) * vc(1:2) )  &
                    / real ( n, kind = rk )
      end do
    end do
  end do

  return
end
subroutine r8mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = rk ) A(M,N), an M by N matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  character ( len = * ) title

  call r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = rk ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ILO, JLO, the first row and column to print.
!
!    Input, integer IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: incx = 5
  integer m
  integer n

  real ( kind = rk ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer i
  integer i2
  integer i2hi
  integer i2lo
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)'
    return
  end if

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(i8,6x)' ) i
    end do

    write ( *, '(''  Row   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc
        i = i2lo - 1 + i2
        write ( ctemp(i2), '(g14.6)' ) a(i,j)
      end do

      write ( *, '(i5,a,5a14)' ) j, ':', ( ctemp(i), i = 1, inc )

    end do

  end do

  return
end
subroutine r8mat_write ( output_filename, m, n, table )

!*****************************************************************************80
!
!! R8MAT_WRITE writes an R8MAT file.
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
!    31 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the output file name.
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) TABLE(M,N), the data.
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
    stop 1
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
