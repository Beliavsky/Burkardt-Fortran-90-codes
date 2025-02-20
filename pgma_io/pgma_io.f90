subroutine ch_cap ( c )

!*****************************************************************************80
!
!! ch_cap() capitalizes a single character.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, the character to capitalize.
!
  implicit none

  character c
  integer itemp

  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

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
!    18 September 2005
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
subroutine getint ( done, ierror, inunit, ival, string )

!*****************************************************************************80
!
!! GETINT reads an integer from a file.
!
!  Discussion:
!
!    The file, or at least the part read by GETINT, is assumed to
!    contain nothing but integers.  These integers may be separated
!    by spaces, or appear on separate lines.  Comments, which begin
!    with "#" and extend to the end of the line, may appear anywhere.
!
!    Each time GETINT is called, it tries to read the next integer
!    it can find.  It remembers where it was in the current line
!    of text.
!
!    The user should open a text file on FORTRAN unit INUNIT,
!    set STRING = ' ' and DONE = TRUE.  The GETINT routine will take
!    care of reading in a new STRING as necessary, and extracting
!    as many integers as possible from the line of text before 
!    reading in the next line.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, logical DONE.
!    On input, if this is the first call, or the user has changed
!    STRING, then set DONE = TRUE.
!    On output, if there is no more data to be read from STRING,
!    then DONE is TRUE.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred while trying to read the integer.
!
!    Input, integer INUNIT, the FORTRAN unit from which to read.
!
!    Output, integer IVAL, the integer that was read.
!
!    Input/output, character ( len = * ) STRING, the text of the most recently 
!    read line of the file.
!
  implicit none

  logical done
  integer i
  integer ierror
  integer inunit
  integer ios
  integer ival
  integer last
  character ( len = * ) string
  character ( len = 80 ) word

  do

    call word_next_rd ( string, word, done )

    if ( .not. done ) then
      exit
    end if

    read ( inunit, '(a)', iostat = ios ) string

    if ( ios /= 0 ) then
      ierror = 1
      return
    end if

    i = index ( string, '#' )
    if ( i /= 0 ) then
      string(i:) = ' '
    end if

  end do

  call s_to_i4 ( word, ival, ierror, last )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GETINT - Fatal error!'
    write ( *, '(a)' ) '  Error trying to convert string to integer.'
    stop
  end if

  return
end
subroutine pgma_check_data ( row_num, col_num, g_max, g, ierror )

!*****************************************************************************80
!
!! PGMA_CHECK_DATA checks ASCII PGM data.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 December 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ROW_NUM,  COL_NUM, the number of rows 
!    and columns of data.
!
!    Input, integer G_MAX, the maximum gray value.
!
!    Input, integer G(ROW_NUM,COL_NUM), contains the gray data.
!
!    Output, integer IERROR, error flag.
!    0, no error detected.
!    1, the data is illegal.
!
  implicit none

  integer col_num
  integer row_num

  integer g(row_num,col_num)
  integer ierror
  integer g_max

  ierror = 0

  if ( minval ( g(1:row_num,1:col_num) ) < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PGMA_CHECK_DATA - Fatal error!'
    write ( *, '(a)' ) '  At least one gray value is below 0.'
    ierror = 1
    stop
  end if

  if ( g_max < maxval ( g(1:row_num,1:col_num) ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PGMA_CHECK_DATA - Fatal error!'
    write ( *, '(a)' ) '  At least one gray value exceeds G_MAX.'
    write ( *, '(a,i12)' ) '  G_MAX = ', g_max
    ierror = 1
    stop
  end if

  return
end
subroutine pgma_example ( row_num, col_num, g )

!*****************************************************************************80
!
!! PGMA_EXAMPLE sets up sample ASCII PGM data.
!
!  Discussion:
!
!    The data is based on three periods of a sine curve.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 December 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ROW_NUM, COL_NUM, the number of rows and 
!    columns of data.  A reasonable value is 200 for ROW_NUM and 600 
!    for COL_NUM.
!
!    Output, integer G(ROW_NUM,COL_NUM), the gray scale data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer row_num
  integer col_num

  integer g(row_num,col_num)
  integer i
  integer j
  integer, parameter :: periods = 3
  real ( kind = rk ), parameter :: pi = 3.14159265D+00
  real ( kind = rk ) x
  real ( kind = rk ) y

  do i = 1, row_num
    y = 2.0D+00 * real ( i - 1, kind = rk ) &
      / real ( row_num - 1, kind = rk ) - 1.0D+00
    do j = 1, col_num
      x = 2.0D+00 * pi * real ( periods * ( j - 1 ), kind = rk ) &
        / real ( col_num - 1, kind = rk )
      g(i,j) = int ( 20.0D+00 * ( sin ( x ) - y + 2.0D+00 ) )
    end do
  end do

  return
end
subroutine pgma_read_data ( file_in_unit, row_num, col_num, g )

!*****************************************************************************80
!
!! PGMA_READ_DATA reads the data in an ASCII PGM file.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer FILE_IN_UNIT, the unit number of the file.
!
!    Input, integer ROW_NUM, COL_NUM, the number of rows and 
!    columns of data.
!
!    Output, integer G(ROW_NUM,COL_NUM), the gray data.
!
  implicit none

  integer col_num
  integer row_num

  logical done
  integer file_in_unit
  integer g(row_num,col_num)
  integer i
  integer ierror
  integer j
  character ( len = 80 ) string

  ierror = 0
  done = .true.
  string = ' '

  do i = 1, row_num
    do j = 1, col_num

      call getint ( done, ierror, file_in_unit, g(i,j), string )

      if ( ierror /= 0 ) then
        close ( unit = file_in_unit )
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PGMA_READ_DATA - Fatal error!'
        write ( *, '(a)' ) '  Problem reading G data.'
        stop
      end if

    end do
  end do

  return
end
subroutine pgma_read_header ( file_in_unit, row_num, col_num, g_max )

!*****************************************************************************80
!
!! PGMA_READ_HEADER reads the header of an ASCII PGM file.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer FILE_IN_UNIT, the unit number of the file.
!
!    Output, integer ROW_NUM, COL_NUM, the number of rows and 
!    columns of data.
!
!    Output, integer G_MAX, the maximum gray value.
!
  implicit none

  logical done
  integer file_in_unit
  integer ierror
  integer ios
  character ( len = 2 )  magic
  integer g_max
  integer col_num
  integer row_num
  logical s_eqi
  character ( len = 80 ) string
!
!  Read the first line of data, which must begin with the magic number.
!
  read ( file_in_unit, '(a)', iostat = ios ) magic

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PGMA_READ_HEADER - Fatal error!'
    write ( *, '(a)' ) '  End or error while reading file.'
    ierror = 2
    stop
  end if

  if ( .not. s_eqi ( magic, 'P2' ) ) then
    ierror = 3
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PGMA_READ_HEADER - Fatal error.'
    write ( *, '(a)' ) '  First two bytes are not magic number "P2".'
    write ( *, '(a)' ) '  First two bytes are: "' // magic // '".'
    stop
  end if
!
!  Now search for COL_NUM, ROW_NUM, and G_MAX.
!
  done = .true.
  string = ' '

  call getint ( done, ierror, file_in_unit, col_num, string )

  if ( ierror /= 0 ) then
    close ( unit = file_in_unit )
    ierror = 4
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PGMA_READ_HEADER - Fatal error!'
    write ( *, '(a)' ) '  Problem reading COL_NUM.'
    stop
  end if

  call getint ( done, ierror, file_in_unit, row_num, string )

  if ( ierror /= 0 ) then
    ierror = 4
    close ( unit = file_in_unit )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PGMA_READ_HEADER - Fatal error!'
    write ( *, '(a)' ) '  Problem reading ROW_NUM.'
    stop
  end if

  call getint ( done, ierror, file_in_unit, g_max, string )

  if ( ierror /= 0 ) then
    ierror = 4
    close ( unit = file_in_unit )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PGMA_READ_HEADER - Fatal error!'
    write ( *, '(a)' ) '  Problem reading G_MAX.'
    stop
  end if

  return
end
subroutine pgma_read_test ( file_in_name, ierror )

!*****************************************************************************80
!
!! PGMA_READ_TEST tests an ASCII PGM file.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_IN_NAME, the name of the file 
!    containing the ASCII PGM data.
!
!    Output, integer IERROR, an error flag which is nonzero if
!    there was an error.
!
  implicit none

  character ( len = * ) file_in_name
  integer file_in_unit
  integer, allocatable, dimension ( :, : ) :: g
  integer ierror
  integer ios
  integer g_max
  integer col_num
  integer row_num

  call get_unit ( file_in_unit )

  open ( unit = file_in_unit, file = file_in_name, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PGMA_READ_TEST - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file.'
    stop
  end if
!
!  Read the header.
!
  call pgma_read_header ( file_in_unit, row_num, col_num, g_max )
!
!  Allocate the data.
!
  allocate ( g(row_num,col_num) )
!
!  Read the data.
!
  call pgma_read_data ( file_in_unit, row_num, col_num, g )

  close ( unit = file_in_unit )
!
!  Check the data.
!
  call pgma_check_data ( row_num, col_num, g_max, g, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PGMA_READ_TEST - Warning!'
    write ( *, '(a)' ) '  PGMA_CHECK_DATA did not approve the data.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PGMA_READ_TEST:'
    write ( *, '(a)' ) '  PGMA_CHECK_DATA has approved the data from the file.'
  end if

  deallocate ( g )

  return
end
subroutine pgma_write ( file_out_name, row_num, col_num, g, ierror )

!*****************************************************************************80
!
!! PGMA_WRITE writes an ASCII PGM file.
!
!  Example:
!
!    P2
!    # feep.pgm
!    24 7
!    15
!    0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
!    0  3  3  3  3  0  0  7  7  7  7  0  0 11 11 11 11  0  0 15 15 15 15  0
!    0  3  0  0  0  0  0  7  0  0  0  0  0 11  0  0  0  0  0 15  0  0 15  0
!    0  3  3  3  0  0  0  7  7  7  0  0  0 11 11 11  0  0  0 15 15 15 15  0
!    0  3  0  0  0  0  0  7  0  0  0  0  0 11  0  0  0  0  0 15  0  0  0  0
!    0  3  0  0  0  0  0  7  7  7  7  0  0 11 11 11 11  0  0 15  0  0  0  0
!    0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
!    
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_OUT_NAME, the name of the file 
!    to which the data should be written.
!
!    Input, integer ROW_NUM, COL_NUM, the number of rows and 
!    columns of data.
!
!    Input, integer G(ROW_NUM,COL_NUM), the gray value of each 
!    pixel.  These should be nonnegative.
!
!    Output, integer IERROR, an error flag.
!    0, no error.
!    1, the data was illegal.
!    2, the file could not be opened.
!
  implicit none

  integer col_num
  integer row_num

  logical, parameter :: debug = .false.
  character ( len = * )  file_out_name
  integer file_out_unit
  integer g(row_num,col_num)
  integer ierror
  integer ios
  integer g_max

  ierror = 0
!
!  Compute the maximum color value.
!
  g_max = maxval ( g(1:row_num,1:col_num) )
!
!  Open the file.
!
  call get_unit ( file_out_unit )

  open ( unit = file_out_unit, file = file_out_name, status = 'replace', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PGMA_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file.'
    ierror = 2
    stop
  end if
!
!  Write the header.
!
  call pgma_write_header ( file_out_name, file_out_unit, row_num, col_num, &
    g_max, ierror )
!
!  Write the data.
!
  call pgma_write_data ( file_out_unit, row_num, col_num, g, ierror )
!
!  Close the file.
!
  close ( unit = file_out_unit )
!
!  Report
!
  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PGMA_WRITE - Note:'
    write ( *, '(a)' ) '  The data was checked and written.'
    write ( *, '(a,i8)' ) '  Number of data rows ROW_NUM =    ', row_num
    write ( *, '(a,i8)' ) '  Number of data columns COL_NUM = ', col_num
    write ( *, '(a,i8)' ) '  Maximum gray value G_MAX =       ', g_max
  end if

  return
end
subroutine pgma_write_data ( file_out_unit, row_num, col_num, g, ierror )

!*****************************************************************************80
!
!! PGMA_WRITE_DATA writes the data of an ASCII PGM file.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    28 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer FILE_OUT_UNIT, the output file unit number.
!
!    Input, integer ROW_NUM, COL_NUM, the number of rows and 
!    columns of data.
!
!    Input, integer G(ROW_NUM,COL_NUM), the gray value of each 
!    pixel.  These should be nonnegative.
!
!    Output, integer IERROR, an error flag.
!    0, no error.
!    1, the data was illegal.
!    2, the file could not be opened.
!
  implicit none

  integer col_num
  integer row_num

  integer file_out_unit
  integer g(row_num,col_num)
  integer i
  integer ierror
  integer jhi
  integer jlo

  ierror = 0

  do i = 1, row_num
    do jlo = 1, col_num, 12
      jhi = min ( jlo + 11, col_num )
      write ( file_out_unit, '(12i5)' ) g(i,jlo:jhi)
    end do
  end do

  return
end
subroutine pgma_write_header ( file_out_name, file_out_unit, row_num, col_num, &
  g_max, ierror )

!*****************************************************************************80
!
!! PGMA_WRITE_HEADER writes the header of an ASCII PGM file.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_OUT_NAME, the name of the output file.
!
!    Input, integer FILE_OUT_UNIT, the output file unit number.
!
!    Input, integer ROW_NUM, COL_NUM, the number of rows and 
!    columns of data.
!
!    Input, integer G_MAX, the maximum gray value.
!
!    Output, integer IERROR, an error flag.
!    0, no error.
!    1, the data was illegal.
!    2, the file could not be opened.
!
  implicit none

  character ( len = * ) file_out_name
  integer file_out_unit
  integer ierror
  character ( len = 2 ) :: magic = 'P2'
  integer g_max
  integer col_num
  integer row_num

  ierror = 0
!
!  Write the header.
!
  write ( file_out_unit, '(a2)' ) magic
  write ( file_out_unit, '(a)' ) '# ' // trim ( file_out_name ) &
    // ' created by PGMA_IO::PGMA_WRITE.'
  write ( file_out_unit, '(i8,2x,i8)' ) col_num, row_num
  write ( file_out_unit, '(i8)' ) g_max

  return
end
subroutine pgma_write_test ( file_out_name )

!*****************************************************************************80
!
!! PGMA_WRITE_TEST tests the ASCII PGM write routines.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_OUT_NAME, the name of the file 
!    to contain the ASCII PGM data.
!
  implicit none

  character ( len = * ) file_out_name
  integer, allocatable, dimension ( :, : ) :: g
  integer ierror
  integer col_num
  integer row_num

  row_num = 300
  col_num = 300
!
!  Allocate memory.
!
  allocate ( g(row_num,col_num) )
!
!  Set the data.
!
  call pgma_example ( row_num, col_num, g )
!
!  Write the data to the file.
!
  call pgma_write ( file_out_name, row_num, col_num, g, ierror )

  deallocate ( g );

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PGMA_WRITE_TEST - Fatal error!'
    write ( *, '(a)' ) '  PGMA_WRITE failed.'
    stop
  end if

  return
end
function s_eqi ( strng1, strng2 )

!*****************************************************************************80
!
!! S_EQI is a case insensitive comparison of two strings for equality.
!
!  Example:
!
!    S_EQI ( 'Anjana', 'ANJANA' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRNG1, STRNG2, the strings to compare.
!
!    Output, logical S_EQI, the result of the comparison.
!
  implicit none

  integer i
  integer len1
  integer len2
  integer lenc
  logical s_eqi
  character s1
  character s2
  character ( len = * ) strng1
  character ( len = * ) strng2

  len1 = len ( strng1 )
  len2 = len ( strng2 )
  lenc = min ( len1, len2 )

  s_eqi = .false.

  do i = 1, lenc

    s1 = strng1(i:i)
    s2 = strng2(i:i)
    call ch_cap ( s1 )
    call ch_cap ( s2 )

    if ( s1 /= s2 ) then
      return
    end if

  end do

  do i = lenc + 1, len1
    if ( strng1(i:i) /= ' ' ) then
      return
    end if
  end do

  do i = lenc + 1, len2
    if ( strng2(i:i) /= ' ' ) then
      return
    end if
  end do

  s_eqi = .true.

  return
end
subroutine s_to_i4 ( s, ival, ierror, last )

!*****************************************************************************80
!
!! S_TO_I4 reads an I4 from a string.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string to be examined.
!
!    Output, integer IVAL, the integer value read from the string.
!    If blank, then IVAL will be returned 0.
!
!    Output, integer IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer LAST, the last character that was
!    part of the representation of IVAL.
!
  implicit none

  character c
  integer i
  integer ierror
  integer isgn
  integer istate
  integer ival
  integer last
  integer lens
  character ( len = * ) s

  ierror = 0
  istate = 0

  isgn = 1
  ival = 0

  lens = len ( s )

  i = 0

  do

    i = i + 1

    c = s(i:i)

    if ( istate == 0 ) then

      if ( c == ' ' ) then

      else if ( c == '-' ) then
        istate = 1
        isgn = -1
      else if ( c == '+' ) then
        istate = 1
        isgn = + 1
      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        exit
      end if

    else if ( istate == 1 ) then

      if ( c == ' ' ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        exit
      end if

    else if ( istate == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        ival = 10 * ival + ichar ( c ) - ichar ( '0' )
      else
        istate = 3
      end if

    end if
!
!  Continue or exit?
!
    if ( istate == 3 ) then
      ival = isgn * ival
      last = i - 1
      exit
    else if ( lens <= i ) then
      if ( istate == 2 ) then
        ival = isgn * ival
        last = lens
      else
        ierror = 1
        last = 0
      end if
      exit
    end if

  end do

  return
end
subroutine word_next_rd ( line, word, done )

!*****************************************************************************80
!
!! WORD_NEXT_RD "reads" words from a string, one at a time.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 June 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) LINE, a string, presumably containing
!    words separated by spaces.
!
!    Output, character ( len = * ) WORD.
!    If DONE is FALSE,
!      WORD contains the "next" word read from LINE.
!    Else
!      WORD is blank.
!
!    Input/output, logical DONE.
!    On input, on the first call, or with a fresh value of LINE,
!      set DONE to TRUE.
!    Else
!      leave it at the output value of the previous call.
!    On output, if a new nonblank word was extracted from LINE
!      DONE is FALSE
!    ELSE
!      DONE is TRUE.
!    If DONE is TRUE, then you need to provide a new LINE of data.
!
!  Local Parameters:
!
!    NEXT is the next location in LINE that should be searched.
!
  implicit none

  logical done
  integer ilo
  integer lenl
  character ( len = * ) line
  integer, save :: next = 1
  character ( len = 1 ), parameter :: TAB = char(9)
  character ( len = * ) word

  lenl = len_trim ( line )

  if ( done ) then
    next = 1
    done = .false.
  end if
!
!  Beginning at index NEXT, search LINE for the next nonblank.
!
  ilo = next

  do
!
!  ...LINE(NEXT:LENL) is blank.  Return with WORD=' ', and DONE=TRUE.
!
    if ( lenl < ilo ) then
      word = ' '
      done = .true.
      next = lenl + 1
      return
    end if
!
!  ...If the current character is blank, skip to the next one.
!
    if ( line(ilo:ilo) /= ' ' .and. line(ilo:ilo) /= TAB ) then
      exit
    end if

    ilo = ilo + 1

  end do
!
!  To get here, ILO must be the index of the nonblank starting
!  character of the next word.
!
!  Now search for the LAST nonblank character.
!
  next = ilo + 1

  do

    if ( lenl < next ) then
      word = line(ilo:next-1)
      return
    end if

    if ( line(next:next) == ' ' .or. line(next:next) == TAB ) then
      exit
    end if

    next = next + 1

  end do

  word = line(ilo:next-1)

  return
end
