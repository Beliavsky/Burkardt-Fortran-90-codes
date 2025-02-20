program main

!*****************************************************************************80
!
!! c_comment() copies a file, converting C++ comments to C comments.
!
!  Discussion:
!
!    Occasionally, I need to modify a C or C++ file so that, instead of using
!    the C++ comment style, it uses the C comment style.
!
!    Although C++ comments can start in any column, I only expect to deal with
!    comments that start in column 1.  The program tries to correctly remark
!    the file so that the C++ comments are replaced by equivalent C comments.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 April 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer  arg_num
  integer  iarg
  integer iargc
  integer ierror
  character ( len = 255 ) input_filename
  character ( len = 255 ) output_filename

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'c_comment():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Modify a C or C++ file so that C++ comments are replaced'
  write ( *, '(a)' ) '  by C-style comments.'

  ierror = 0

  arg_num = iargc ( )
!
!  If at least one command line argument, it's the input file name.
!
  if ( 1 <= arg_num ) then

    iarg = 1

    call getarg ( iarg, input_filename )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'c_comment():'
    write ( *, '(a)' ) '  Please enter the name of the input file.'

    read ( *, '(a)' ) input_filename

  end if
!
!  If at least two command line arguments, it's the output file name.
!
  if ( 2 <= arg_num ) then

    iarg = 2

    call getarg ( iarg, output_filename )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'c_comment():'
    write ( *, '(a)' ) '  Please enter the name of the output file.'

    read ( *, '(a)' ) output_filename

  end if

  call handle_file ( input_filename, output_filename )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'c_comment():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
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
!    This code is distributed under the GNU LGPL license. 
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
subroutine handle_file ( input_filename, output_filename )

!*****************************************************************************80
!
!! handle_file() rewrites a file so it uses C-style comments instead of C++.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 June 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the name of the output file.
!
  implicit none

  logical comment1
  logical comment2
  logical comment3
  character ( len = * ) input_filename
  integer input_unit
  integer ios
  character ( len = * ) output_filename
  integer output_unit
  character ( len = 255 ) string1
  character ( len = 255 ) string2
  character ( len = 255 ) string3

  call get_unit ( input_unit )
  open ( unit = input_unit, file = input_filename, status = 'old' )

  call get_unit ( output_unit )
  open ( unit = output_unit, file = output_filename, status = 'replace' )

  string1 = ' '
  comment1 = .false.
  string2 = '/*  Output from c_comment */'
  comment2 = .false.
  string3 = ' '
  comment3 = .false.
  
  do

    string1 = string2
    comment1 = comment2
    string2 = string3
    comment2 = comment3
    read ( input_unit, '(a)', iostat = ios ) string3
!
!  Special actions when end of file is reached.
!
    if ( ios /= 0 ) then
      exit
    end if
!
!  Is this a comment?
!
    if ( string3(1:2) == '//' ) then
      comment3 = .true.
    else
      comment3 = .false.
    end if
!
!  Decide how to handle string 2.
!
    if ( comment2 ) then
      if ( comment1 .and. comment3 ) then
        string2 = string2(3:)
      else if ( comment1 .and. ( .not. comment3 ) ) then
        string2 = trim ( string2(3:) ) // '*/'
      else if ( ( .not. comment1 ) .and. comment3 ) then
        string2 = '/*' // trim ( string2(3:) )
      else
        string2 = '/*' // trim ( string2(3:) ) // '*/'
      end if
    end if

    write ( output_unit, '(a)' ) trim ( string1 )

  end do

  write ( output_unit, '(a)' ) trim ( string1 )
  write ( output_unit, '(a)' ) trim ( string2 )

  close ( unit = input_unit )
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
!    This code is distributed under the GNU LGPL license.
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
