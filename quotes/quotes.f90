program main

!*****************************************************************************80
!
!! quotes() extracts a random quote from a file.
!
!  Discussion:
!
!    This program is descended from a program written by David Moses.
!
!    The quote file contains a series of quotes, separated by a blank line.
!    When printing out a quote from the file, any line that ends with a
!    quotation mark is followed by an extra blank line.
!
!    Note that my version of Linux includes a useless and undocumented
!    command "quote" which quotes a string.  Hence this program is
!    renamed "quotes".
!
!  Usage:
!
!    quotes
!
!      extracts a quote from the default quote file.
!
!    quotes FILE
!
!      extracts a quote from the file FILE.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  logical, parameter :: DEBUG = .false.
  integer i4_uniform_ab
  integer iarg
  integer iargc
  integer ierror
  character ( len = 255 ) line
  integer line_index
  integer line_num
  character ( len = 255 ) list_file_name
  integer numarg
  character ( len = 255 ) quote_file_name
  character ( len = 255 ) quote_file_title
  integer quote_index
  integer quote_num

  if ( DEBUG ) then
    call timestamp ( )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'quotes():'
    write ( *, '(a)' ) '  Fortran90 version'
  end if

  write ( *, '(a)' ) ' '

  call timestamp ( )
!
!  Count the command line arguments.
!
  numarg = iargc ( )
!
!  If a quote file was not specified, open QUOTE_FILES.TXT and
!  pick a file at random.
!
!  Linux laptop: /home/john/public_html/f_src/quotes/quote_files.txt
!
  if ( numarg == 0 ) then

    list_file_name = '/home/john/public_html/f_src/quotes/quote_files.txt'

    call file_line_count ( list_file_name, line_num )

    if ( line_num < 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'quotes(): Fatal error.'
      write ( *, '(a)' ) '  The list of quote files is missing!'
      stop 1
    end if

    if ( DEBUG ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'quotes(): Debug'
      write ( *, '(a,i8)' ) &
        '  The number of quote file listing entries is ', line_num
    end if

    line_index = i4_uniform_ab ( 1, line_num )

    if ( DEBUG ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'quotes(): Debug'
      write ( *, '(a,i8)' ) '  Choosing file number ', line_index
    end if

    call file_line_get ( list_file_name, line_index, line )

    call word_next2 ( line, quote_file_name, quote_file_title )

    if ( DEBUG ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'quotes(): Debug'
      write ( *, '(a)' ) trim ( quote_file_title )
    end if

  else

    iarg = 1
    call getarg ( iarg, quote_file_name )

    quote_file_title = 'One of my favorite quotations.'

  end if
!
!  Compute QUOTE_NUM, the number of quotes in the file.
!
  call file_para_count ( quote_file_name, quote_num )

  if ( quote_num <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'quotes(): Fatal error!'
    write ( *, '(a)' ) '  file_para_count() could not count the quotes.'
    stop 1
  end if
!
!  Set QUOTE_INDEX, the number of the quote to be displayed.
!
  quote_index = i4_uniform_ab ( 1, quote_num )
 
  if ( quote_index <= 0 ) then
    quote_index = 1
  end if
 
  if ( quote_num < quote_index ) then
    quote_index = quote_num
  end if

  if ( DEBUG ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'quotes(): Debug'
    write ( *, '(a,i8)' ) '  Extract quote number ', quote_index
  end if
!
!  Read from the quote file until you reach the desired quote.
!  QUOTE_INDEX tells us which quote we are currently reading.
!  Quotes are presumed to be separated by a single blank line.
!
  call quote_file_print ( quote_file_name, quote_file_title, &
    quote_index, quote_num, ierror )
!
!  Terminate.
!
  if ( DEBUG ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'quotes():'
    write ( *, '(a)' ) '  Normal end of execution.'
    write ( *, '(a)' ) ' '
    call timestamp ( )
  end if

  stop
end
subroutine file_line_count ( file_name, line_num )

!*****************************************************************************80
!
!! file_line_count() counts the number of lines in a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.  
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 August 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    character ( len = * ) FILE_NAME, the name of the file.
!
!  Output:
!
!    integer LINE_NUM, the number of lines found in the file.
!
  implicit none

  character ( len = * ) file_name
  integer ios
  integer iunit
  integer line_num

  line_num = 0
!
!  Open the file.
!
  call get_unit ( iunit )

  if ( iunit == 0 ) then
    line_num = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'file_line_count(): Fatal error!'
    write ( *, '(a)' ) '  Could not find a free FORTRAN unit.'
    return
  end if

  open ( unit = iunit, file = file_name, status = 'old', form = 'formatted', &
    access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    line_num = - 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_LINE_COUNT - Fatal error!'
    write ( *, '(a)' ) &
      '  Could not open the file "' // trim ( file_name ) // '".'
    return
  end if
!
!  Count the lines.
!
  do

    read ( iunit, '(a)', iostat = ios )

    if ( ios /= 0 ) then
      exit
    end if

    line_num = line_num + 1

  end do

  close ( unit = iunit )

  return
end
subroutine file_line_get ( file_name, line_index, line )

!*****************************************************************************80
!
!! file_line_get() gets a particular line of a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.  
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 August 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    character ( len = * ) FILE_NAME, the name of the file.
!
!    integer LINE_INDEX, the line number to be read.
!
!  Output:
!
!    character ( len = * ) LINE, the text of the line.
!
  implicit none

  character ( len = * ) file_name
  integer i
  integer ios
  integer iunit
  character ( len = * ) line
  integer line_index
!
!  Open the file.
!
  call get_unit ( iunit )

  if ( iunit == 0 ) then
    line = ' '
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_LINE_GET - Fatal error!'
    write ( *, '(a)' ) '  Could not find a free FORTRAN unit.'
    return
  end if

  open ( unit = iunit, file = file_name, status = 'old', form = 'formatted', &
    access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    line = ' '
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_LINE_GET - Fatal error!'
    write ( *, '(a)' ) &
      '  Could not open the the file "' // trim ( file_name ) // '".'
    return
  end if
!
!  Count the lines.
!
  do i = 1, line_index

    read ( iunit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      line = ' '
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FILE_LINE_GET - Fatal error!'
      write ( *, '(a)' ) '  Unexpected end of file.'
      return
    end if

  end do

  close ( unit = iunit )

  return
end
subroutine file_para_count ( file_name, para_num )

!*****************************************************************************80
!
!! file_para_count() counts the number of paragraphs in a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.  A paragraph is
!    a sequence of nonblank lines.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 August 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    character ( len = * ) FILE_NAME, the name of the file.
!
!  Output:
!
!    integer PARA_NUM, the number of paragraphs found in the file.
!
  implicit none

  character ( len = * ) file_name
  integer ios
  integer iunit
  integer lenc
  integer lenc_old
  character ( len = 255 ) line
  integer para_num

  para_num = 0
!
!  Open the file.
!
  call get_unit ( iunit )

  if ( iunit == 0 ) then
    para_num = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'file_para_count(): Fatal error!'
    write ( *, '(a)' ) '  Could not find a free FORTRAN unit.'
    return
  end if

  open ( unit = iunit, file = file_name, status = 'old', form = 'formatted', &
    access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    para_num = - 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'file_para_count(): Fatal error!'
    write ( *, '(a)') &
      '  Could not open the the file "' // trim ( file_name ) // '".'
    return
  end if
!
!  Count the paragraphs.
!
  lenc = 0

  do

    read ( iunit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    lenc_old = lenc
    lenc = len_trim ( line )

    if ( 0 < lenc .and. lenc_old <= 0 ) then
      para_num = para_num + 1
    end if

  end do

  close ( unit = iunit )

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
!    27 August 2021
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
function i4_uniform_ab ( a, b )

!*****************************************************************************80
!
!! i4_uniform_ab() returns a scaled pseudorandom I4.
!
!  Discussion:
!
!    An I4 is an integer value.
!
!    The pseudorandom number will be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 August 2021
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Input:
!
!    integer A, B, the limits of the interval.
!
!  Output:
!
!    integer I4_UNIFORM_AB, a number between A and B.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer a
  integer b
  integer i4_uniform_ab
  real ( kind = rk ) r

  call random_number ( harvest = r )
  i4_uniform_ab = int ( a + ( b + 1 - a ) * r )

  return
end
subroutine quote_file_print ( quote_file_name, quote_file_title, quote_index, &
  quote_num, ierror )

!*****************************************************************************80
!
!! quote_file_print() prints a given quote from a quote file.
!
!  Discussion:
!
!    The index of the desired quote is specified.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 August 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    character ( len = * ) QUOTE_FILE_NAME, the name of the quote file.
!
!    character ( len = * ) QUOTE_FILE_TITLE, the title of the quote file.
!
!    integer QUOTE_INDEX, the index of the quote to 
!    be printed.  The first quote has QUOTE_INDEX = 1, and so on.  If there 
!    are fewer that QUOTE_INDEX quotes in the file, then nothing will 
!    be printed.
!
!    integer QUOTE_NUM, the number of quotes in the file.
!
!  Output:
!
!    integer IERROR, error flag.
!    0, no error was encountered.
!    nonzero, an error occurred.
!
  implicit none

  integer ierror
  integer ios
  integer jquote
  integer lenc
  integer lencold
  character ( len = 255 ) line
  character ( len = * ) quote_file_name
  character ( len = * ) quote_file_title
  integer quote_file_unit
  integer quote_index
  integer quote_num
  character ( len = 6 ) s
!
!  Open the quote file.
!
  call get_unit ( quote_file_unit )

  if ( quote_file_unit == 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'QUOTE_FILE_PRINT - Fatal error!'
    write ( *, '(a)' ) '  Could not find a free FORTRAN unit.'
    return
  end if

  open ( unit = quote_file_unit, file = quote_file_name, status = 'old', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'QUOTE_FILE_PRINT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the quote the file "' &
      // trim ( quote_file_name ) // '".'
    return
  end if

  jquote = 1
  lenc = 0
 
  write ( *, '(a)' ) ' '
 
  do
 
    read ( quote_file_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then

      if ( quote_index /= jquote ) then
        ierror = 2
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'QUOTE_FILE_PRINT - Fatal error!'
        write ( *, '(a)' ) '  Could not find the desired quote.'
      end if
      close ( unit = quote_file_unit )
      return
    end if

    lencold = lenc
    lenc = len_trim ( line )
 
    if ( 0 < lenc .and. lencold <= 0 ) then
      jquote = jquote + 1
    end if

    if ( jquote < quote_index ) then

    else if ( jquote == quote_index ) then

      write ( *, '(2x,a)' ) trim ( line )
      if ( 0 < lenc ) then
        if ( line(lenc:lenc) == '"' ) then
          write ( *, '(a)' ) ' '
        end if
      end if
 
    else if ( quote_index < jquote ) then

      exit

    end if
 
  end do
 
  write ( s, '(i6)' ) quote_num + 1 - jquote
  write ( *, '(2x,a)' ) 'Quotation ' // trim ( adjustl ( s ) ) 
  write ( *, '(2x,a)' ) trim ( quote_file_title )
  write ( *, '(a)' ) ' '

  close ( unit = quote_file_unit )

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
!    27 August 2021
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
subroutine word_next2 ( s, first, last )

!*****************************************************************************80
!
!! word_next2() returns the first word in a string.
!
!  Discussion:
!
!    "Words" are any string of characters, separated by commas or blanks.
!
!    The routine returns:
!    * FIRST, the first string of nonblank, noncomma characters;
!    * LAST, the characters of the string that occur after FIRST and
!      the commas and blanks.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 August 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    character ( len = * ) S, the string of words to be analyzed.
!
!  Output:
!
!    character ( len = * ) FIRST, the next word in the string.
!
!    character ( len = * ) LAST, the remaining string.
!
  implicit none

  character c
  character ( len = * ) first
  integer i
  integer ido
  integer ifirst
  integer ilast
  character ( len = * ) last
  integer lenf
  integer lenl
  integer lens
  character ( len = * ) s

  first = ' '
  last = ' '

  ifirst = 0
  ilast = 0

  lens = len_trim ( s )
  lenf = len ( first )
  lenl = len ( last )

  ido = 0

  do i = 1, lens

    c = s(i:i)

    if ( ido == 0 ) then
      if ( c /= ' ' .and. c /= ',' ) then
        ido = 1
      end if
    end if

    if ( ido == 1 ) then
      if ( c /= ' ' .and. c /= ',' ) then
        ifirst = ifirst + 1
        if ( ifirst <= lenf ) then
          first(ifirst:ifirst) = c
        end if
      else
        ido = 2
      end if
    end if

    if ( ido == 2 ) then
      if ( c /= ' ' .and. c /= ',' ) then
        ido = 3
      end if
    end if

    if ( ido == 3 ) then
      ilast = ilast + 1
      if ( ilast <= lenl ) then
        last(ilast:ilast) = c
      end if
    end if

  end do

  return
end
 
