subroutine pentomino_matrix ( name, p_m, p_n, p )

!*****************************************************************************80
!
!! PENTOMINO_MATRIX returns a 0/1 matrix defining a particular pentomino.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 April 2018
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( length = 1 ) NAME, is 'f', 'i', 'l', 'n', 'p', 't', 
!    'u', 'v', 'w', 'x', 'y' or 'z'.
!
!    Output, integer ( kind = 4 ) P_M, P_N, the number of rows and columns of 
!    the representation.
!
!    Output, integer ( kind = 4 ) P(P_Mp_N), a P_MxP_N matrix of 0's and 1's 
!    that indicates the shape of the pentomino.
!
  implicit none

  character ( len = 1 ) name
  integer ( kind = 4 ) p(5*5)
  integer ( kind = 4 ) p_m
  integer ( kind = 4 ) p_n

  integer ( kind = 4 ), save :: f_mat(3*3) = (/ &
     0, 1, 0, &
     1, 1, 1, &
     1, 0, 0 /)
  integer ( kind = 4 ), save :: i_mat(5*1) = (/ &
    1, 1, 1, 1, 1 /)
  integer ( kind = 4 ), save :: l_mat(4*2) = (/ &
    1, 1, 1, 1, &
    0, 0, 0, 1 /)
  integer ( kind = 4 ), save :: n_mat(2*4) = (/ &
    1, 0, &
    1, 1, &
    0, 1, &
    0, 1 /)
  integer ( kind = 4 ), save :: p_mat(3*2) = (/ &
    1, 1, 1, &
    1, 1, 0 /)
  integer ( kind = 4 ), save :: t_mat(3*3) = (/ &
    1, 0, 0, &
    1, 1, 1, &
    1, 0, 0 /)
  integer ( kind = 4 ), save :: u_mat(2*3) = (/ &
    1, 1, &
    0, 1, &
    1, 1 /)
  integer ( kind = 4 ), save :: v_mat(3*3) = (/ &
    1, 1, 1, &
    0, 0, 1, &
    0, 0, 1 /)
  integer ( kind = 4 ), save :: w_mat(3*3) = (/ &
    1, 1, 0, &
    0, 1, 1, &
    0, 0, 1 /)
  integer ( kind = 4 ), save :: x_mat(3*3) = (/ &
    0, 1, 0, &
    1, 1, 1, &
    0, 1, 0 /)
  integer ( kind = 4 ), save :: y_mat(2*4) = (/ &
    0, 1, &
    0, 1, &
    1, 1, &
    0, 1 /)
  integer ( kind = 4 ), save :: z_mat(3*3) = (/ &
    1, 0, 0, &
    1, 1, 1, &
    0, 0, 1 /)

  p(1:5*5) = 0

  if ( name == 'f' .or. name == 'F' ) then
    p_m = 3
    p_n = 3
    p(1:p_m*p_n) = f_mat(1:p_m*p_n)
  else if ( name == 'i' .or. name == 'I' ) then
    p_m = 5
    p_n = 1
    p(1:p_m*p_n) = i_mat(1:p_m*p_n)
  else if ( name == 'l' .or. name == 'L' ) then
    p_m = 4
    p_n = 2
    p(1:p_m*p_n) = l_mat(1:p_m*p_n)
  else if ( name == 'n' .or. name == 'N' ) then
    p_m = 2
    p_n = 4
    p(1:p_m*p_n) = n_mat(1:p_m*p_n)
  else if ( name == 'p' .or. name == 'P' ) then
    p_m = 3
    p_n = 2
    p(1:p_m*p_n) = p_mat(1:p_m*p_n)
  else if ( name == 't' .or. name == 'T' ) then
    p_m = 3
    p_n = 3
    p(1:p_m*p_n) = t_mat(1:p_m*p_n)
  else if ( name == 'u' .or. name == 'U' ) then
    p_m = 2
    p_n = 3
    p(1:p_m*p_n) = u_mat(1:p_m*p_n)
  else if ( name == 'v' .or. name == 'V' ) then
    p_m = 3
    p_n = 3
    p(1:p_m*p_n) = v_mat(1:p_m*p_n)
  else if ( name == 'w' .or. name == 'W' ) then
    p_m = 3
    p_n = 3
    p(1:p_m*p_n) = w_mat(1:p_m*p_n)
  else if ( name == 'x' .or. name == 'X' ) then
    p_m = 3
    p_n = 3
    p(1:p_m*p_n) = x_mat(1:p_m*p_n)
  else if ( name == 'y' .or. name == 'Y' ) then
    p_m = 2
    p_n = 4
    p(1:p_m*p_n) = y_mat(1:p_m*p_n)
  else if ( name == 'z' .or. name == 'Z' ) then
    p_m = 3
    p_n = 3
    p(1:p_m*p_n) = z_mat(1:p_m*p_n)
  else
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'PENTOMINO_MATRIX - Fatal error!'
    write ( *, '(a)' ) '  Illegal name = "' // name // '"'
    write ( *, '(a)' ) '  Legal names: f, i, l, n, p, t, u, v, w, x, y, z.'
    stop 1
  end if

  return
end
subroutine pentomino_plot ( p_m, p_n, p, label )

!*****************************************************************************80
!
!! PENTOMINO_PLOT plots a particular pentomino in a 5x5 grid.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 April 2018
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P_M, P_N, the number of rows and columns of the
!    representation.
!
!    Input, integer ( kind = 4 ) P(P_M*P_N), a matrix of 0's and 1's.
!    1 <= P_M, P_N <= 5.  There should be exactly 5 values of one.
!
!    Input, string LABEL, a title for the plot.
!
  implicit none

  integer ( kind = 4 ) p_m
  integer ( kind = 4 ) p_n

  character ( len = 16 ) color
  integer ( kind = 4 ) color_index(5,5)
  character ( len = 80 ) command_filename
  integer ( kind = 4 ) command_unit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_reverse
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  character ( len = * ) label
  integer ( kind = 4 ) :: m = 5
  integer ( kind = 4 ) :: n = 5
  integer ( kind = 4 ) p(p_m*p_n)
  character ( len = 80 ) plot_filename

  command_filename = trim ( label ) // '_commands.txt'
  plot_filename = trim ( label ) // '.png'
!
!  Initially, the grid is entirely white (color 0)
!
  color_index(1:m,1:n) = 0
!
!  Place the pentomino on the grid, so that it is "snug" in the upper left corner.
!
  do j = 1, p_n
    do i = 1, p_m
      color_index(i,j) = p(i+(j-1)*p_m)
    end do
  end do
!
!  Create the command file.
!
  call get_unit ( command_unit )
  open ( unit = command_unit, file = command_filename, &
    status = 'replace' )
  write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < ' // &
    trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) 'set output "' // &
    trim ( plot_filename ) // '"'
  write ( command_unit, '(a)' ) 'set title "' // label // '"'
!
!  Get a plot of TRUE SQUARES.
!
  write ( command_unit, '(a)' ) 'set xrange [ 0 : 5 ]'
  write ( command_unit, '(a)' ) 'set yrange [ 0 : 5 ]'
  write ( command_unit, '(a)' ) 'set size square'
  write ( command_unit, '(a)' ) 'unset border'
  write ( command_unit, '(a)' ) 'unset tics'
  write ( command_unit, '(a)' ) 'set nokey'

  k = 0
  do i = 1, m
    do j = 1, n
      k = k + 1

      if ( color_index(i,j) == 0 ) then
        color = "white"
      else if ( color_index(i,j) == 1 ) then
        color = "black"
      end if

      i_reverse = m + 1 - i

      write ( command_unit, '(a,i4,a,i2,a,i2,a,i2,a,i2,a)' ) &
        'set object ', k, ' rect from ', j-1, ',', i_reverse-1, &
        ' to ', j, ',', i_reverse, ' back'
      write ( command_unit, '(a,i4,a)' ) &
        'set object ', k, ' rect fc rgb "' // trim ( color ) // &
        '" fillstyle solid 1.0'

    end do
  end do
!
!  If you don't have some bogus PLOT command here, all the previous work
!  results in no plot all.  Way to go, gnuplot!
!  Here, we plot the function y = -1, which is out of range and won't show up.
!
  write ( command_unit, '(a)' ) 'plot -1 with lines'
  close ( unit = command_unit )

  write ( *, '(a)' ) &
    '  Created command file "' // trim ( command_filename ) // '".'

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
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
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
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
