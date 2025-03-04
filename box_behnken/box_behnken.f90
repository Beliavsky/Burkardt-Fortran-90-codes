subroutine box_behnken ( dim_num, x_num, range, x )

!*****************************************************************************80
!
!! box_behnken() returns a Box-Behnken design for the given number of factors.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    George Box, Donald Behnken,
!    Some new three level designs for the study of quantitative variables,
!    Technometrics,
!    Volume 2, pages 455-475, 1960.
!
!  Input:
!
!    integer DIM_NUM, the spatial dimension.
!
!    integer X_NUM, the number of elements of the design.
!    X_NUM should be equal to DIM_NUM * 2**(DIM_NUM-1) + 1.
!
!    real ( kind = rk ) RANGE(DIM_NUM,2), the minimum and maximum
!    value for each component.
!
!  Output:
!
!    real ( kind = rk ) X(DIM_NUM,X_NUM), the elements of the design.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num
  integer x_num

  integer i
  integer i2
  integer j
  integer last_low
  real ( kind = rk ) range(dim_num,2)
  real ( kind = rk ) x(dim_num,x_num)
!
!  Ensure that the range is legal.
!
  if ( any ( range(1:dim_num,2) <= range(1:dim_num,1) ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BOX_BEHNKEN - Fatal error!'
    write ( *, '(a)' ) '  For some index I,'
    write ( *, '(a)' ) '  RANGE(I,2) <= RANGE(I,1).'
    stop
  end if
!
!  The first point is the center.
!
  j = 1

  x(1:dim_num,j) = ( range(1:dim_num,1) + range(1:dim_num,2) ) / 2.0D+00
!
!  For subsequent elements, one entry is fixed at the middle of the range.
!  The others are set to either extreme.
!
  do i = 1, dim_num

    j = j + 1

    x(1:dim_num,j) = range(1:dim_num,1)
    x(i,j) = ( range(i,1) + range(i,2) ) / 2.0D+00
!
!  The next element is made by finding the last low value, making it
!  high, and all subsequent high values low.
!
    do

      last_low = -1

      do i2 = 1, dim_num
        if ( x(i2,j) == range(i2,1) ) then
          last_low = i2
        end if
      end do

      if ( last_low == -1 ) then
        exit
      end if

      j = j + 1
      x(1:dim_num,j) = x(1:dim_num,j-1)
      x(last_low,j) = range(last_low,2)

      do i2 = last_low + 1, dim_num
        if ( x(i2,j) == range(i2,2) ) then
          x(i2,j) = range(i2,1)
        end if
      end do

    end do

  end do

  return
end
subroutine box_behnken_size ( dim_num, x_num )

!*****************************************************************************80
!
!! box_behnken_size() returns the size of a Box-Behnken design.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    George Box, Donald Behnken,
!    Some new three level designs for the study of quantitative variables,
!    Technometrics,
!    Volume 2, pages 455-475, 1960.
!
!  Input:
!
!    integer DIM_NUM, the spatial dimension.
!
!  Output:
!
!    integer X_NUM, the number of elements of the design.
!    X_NUM will be equal to DIM_NUM * 2^(DIM_NUM-1) + 1.
!
  implicit none

  integer dim_num
  integer x_num

  if ( 1 <= dim_num ) then
    x_num = 1 + dim_num * 2**( dim_num - 1 )
  else
    x_num = -1
  end if

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! get_unit() returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
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
subroutine r8mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! r8mat_transpose_print() prints an R8MAT, transposed.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, N, the number of rows and columns.
!
!    real ( kind = rk ) A(M,N), an M by N matrix to be printed.
!
!    character ( len = * ) TITLE, a title.
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
!! r8mat_transpose_print_some() prints some of an R8MAT, transposed.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, N, the number of rows and columns.
!
!    real ( kind = rk ) A(M,N), an M by N matrix to be printed.
!
!    integer ILO, JLO, the first row and column to print.
!
!    integer IHI, JHI, the last row and column to print.
!
!    character ( len = * ) TITLE, a title.
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

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(i7,7x)') i
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

      write ( *, '(i5,1x,5a14)' ) j, ( ctemp(i), i = 1, inc )

    end do

  end do

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
!    03 September 2021
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
!    03 September 2021
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
 
