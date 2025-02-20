program main

!*****************************************************************************80
!
!! MAIN is the main program for SPHERE_DELAUNAY.
!
!  Discussion:
!
!    SPHERE_DELAUNAY determines the Delaunay triangulation of points 
!    on the unit sphere in 3D.
!
!    The user supplies an XYZ file containing points on the unit sphere.
!    This routine reads that file and calls the appropriate STRIPACK
!    routines to compute and plot the Delaunay triangulation. 
!
!    These plots are stored as PostScript files.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 December 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer arg_num
  integer iarg
  integer iargc
  integer ierror
  integer ios
  character ( len = 255 ) point_file_name

  ierror = 0

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPHERE_DELAUNAY'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Given a set of XYZ points on the sphere,'
  write ( *, '(a)' ) '  compute the Delaunay triangulation.'
!
!  Get the number of command line arguments.
!
  arg_num = iargc ( )
!
!  If at least one command line argument, it's the input file name.
!
  if ( arg_num < 1 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Enter the input XYZ file name:'
    read ( *, '(a)', iostat = ios ) point_file_name

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPHERE_DELAUNAY - Fatal error!'
      write ( *, '(a)' ) '  Unexpected read error!'
      stop
    end if

  else

    iarg = 1
    call getarg ( iarg, point_file_name )

  end if

  call stripack_interface ( point_file_name )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPHERE_DELAUNAY'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine stripack_interface ( point_file_name )

!*****************************************************************************80
!
!! STRIPACK_INTERFACE calls STRIPACK routines.
!
!  Discussion:
!
!    10 June 2002: Changed this routine so that it can handle any size
!    problem, by making all arrays allocatable.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) POINT_FILE_NAME, the name of the input file.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  character ( len = 255 ) delaunay_plot_file_name
  real ( kind = rk ), allocatable, dimension ( : ) :: ds
  real ( kind = rk ) elat
  real ( kind = rk ) elon
  integer i
  integer ierror
  integer iunit
  integer, allocatable, dimension ( : ) :: iwk
  integer, allocatable, dimension ( :, : ) :: lbtri
  integer, allocatable, dimension ( : ) :: lend
  integer, allocatable, dimension ( : ) :: list
  integer, allocatable, dimension ( : ) :: listc
  integer lnew
  integer, allocatable, dimension ( : ) :: lptr
  integer, allocatable, dimension ( :, : ) :: ltri
  integer n
  real ( kind = rk ) norm
  integer nt
  logical numbr
  real ( kind = rk ), parameter :: pltsiz = 7.5D+00
  character ( len = * ) point_file_name
  real ( kind = rk ), allocatable, dimension ( : ) :: rc
  character ( len = 255 ) trplot_title
  real ( kind = rk ), allocatable, dimension ( : ) :: x
  real ( kind = rk ), allocatable, dimension ( : ) :: xc
  real ( kind = rk ), allocatable, dimension ( : ) :: y
  real ( kind = rk ), allocatable, dimension ( : ) :: yc
  real ( kind = rk ), allocatable, dimension ( : ) :: z
  real ( kind = rk ), allocatable, dimension ( : ) :: zc
!
!  Count the number of lines of (X,Y,Z) data.
!
  call file_row_count ( point_file_name, n )

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRIPACK_INTERFACE - Fatal error!'
    write ( *, '(a)' ) '  The input file has no data.'
    return
  end if
!
!  Allocate everything.
!
  allocate ( ds(1:n) )
  allocate ( iwk(1:2*n) )
  allocate ( lbtri(1:6,1:n) )
  allocate ( lend(1:n) )
  allocate ( list(1:6*(n-2)) )
  allocate ( listc(1:6*(n-2)) )
  allocate ( lptr(1:6*(n-2)) )
  allocate ( ltri(1:9,1:2*(n-2)) )
  allocate ( rc(1:2*(n-2)) )
  allocate ( x(1:n) )
  allocate ( xc(1:2*(n-2)) )
  allocate ( y(1:n) )
  allocate ( yc(1:2*(n-2)) )
  allocate ( z(1:n) )
  allocate ( zc(1:2*(n-2)) )
!
!  Read the (X,Y,Z) data from a file.
!
  call xyz_read ( point_file_name, n, x, y, z, ierror )
!
!  Make sure the data is on the unit sphere.
!
  do i = 1, n
    norm = sqrt ( x(i)**2 + y(i)**2 + z(i)**2 )
    x(i) = x(i) / norm
    y(i) = y(i) / norm
    z(i) = z(i) / norm
  end do
!
!  Create the triangulation.
!
  call trmesh ( n, x, y, z, list, lptr, lend, lnew, iwk, iwk(n+1), ds, ierror )

  if ( ierror == -2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRIPACK_INTERFACE - Fatal error!'
    write ( *, '(a)' ) '  Error in TRMESH.'
    write ( *, '(a)' ) '  The first 3 nodes are collinear.'
    stop
  end if

  if ( 0 < ierror ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRIPACK_INTERFACE - Fatal error!'
    write ( *, '(a)' ) '  Error in TRMESH.'
    write ( *, '(a)' ) '  Duplicate nodes encountered.'
    stop
  end if
!
!  Create a triangle list.
!
  call trlist ( n, list, lptr, lend, 9, nt, ltri, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRIPACK_INTERFACE - Fatal error!'
    write ( *, '(a)' ) '  Error in TRLIST.'
    stop
  end if
!
!  Plot the portion of the triangulation contained 
!  in the hemisphere centered at E = (ELAT,ELON), where ELAT and ELON
!  are taken to be the center of the range of
!  the nodal latitudes and longitudes.
!
  elat = 0.0D+00
  elon = 0.0D+00
  a = 90.0D+00
  numbr = ( n <= 200 )

  trplot_title = '(' // trim ( point_file_name ) // ')'

  delaunay_plot_file_name = point_file_name
  call file_name_ext_swap ( delaunay_plot_file_name, 'eps' )

  call get_unit ( iunit )

  open ( unit = iunit, file = delaunay_plot_file_name )

  call trplot ( iunit, pltsiz, elat, elon, a, n, x, y, z, list, &
    lptr, lend, trplot_title, numbr, ierror )

  close ( unit = iunit )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRIPACK_INTERFACE - Warning!'
    write ( *, '(a,i8)' ) '  TRPLOT returned error code ', ierror
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRPLOT created the triangulation plot file: "' // &
    trim ( delaunay_plot_file_name ) // '".'
!
!  Write the XYZL file that indexes the points that form Delaunay 
!  triangulation lines.
!
  call tr_to_xyzl ( n, x, y, z, list, lptr, lend, point_file_name )
!
!  Free memory.
!
  deallocate ( ds )
  deallocate ( iwk )
  deallocate ( lbtri )
  deallocate ( lend )
  deallocate ( list )
  deallocate ( listc )
  deallocate ( lptr )
  deallocate ( ltri )
  deallocate ( rc )
  deallocate ( x )
  deallocate ( xc )
  deallocate ( y )
  deallocate ( yc )
  deallocate ( z )
  deallocate ( zc )

  return
end
subroutine ch_cap ( c )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
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

  integer, parameter :: rk = kind ( 1.0D+00 )

  character c
  integer itemp

  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end
function ch_eqi ( c1, c2 )

!*****************************************************************************80
!
!! CH_EQI is a case insensitive comparison of two characters for equality.
!
!  Example:
!
!    CH_EQI ( 'A', 'a' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    28 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C1, C2, the characters to compare.
!
!    Output, logical CH_EQI, the result of the comparison.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character c1
  character c1_cap
  character c2
  character c2_cap
  logical ch_eqi

  c1_cap = c1
  c2_cap = c2

  call ch_cap ( c1_cap )
  call ch_cap ( c2_cap )

  if ( c1_cap == c2_cap ) then
    ch_eqi = .true.
  else
    ch_eqi = .false.
  end if

  return
end
subroutine ch_to_digit ( c, digit )

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the numeric value of a base 10 digit.
!
!  Example:
!
!     C   DIGIT
!    ---  -----
!    '0'    0
!    '1'    1
!    ...  ...
!    '9'    9
!    ' '    0
!    'X'   -1
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the decimal digit, '0' through '9' or blank
!    are legal.
!
!    Output, integer DIGIT, the numeric value.  If C was
!    'illegal', then DIGIT is -1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character c
  integer digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

    digit = ichar ( c ) - 48

  else if ( c == ' ' ) then

    digit = 0

  else

    digit = -1

  end if

  return
end
subroutine file_name_ext_get ( file_name, i, j )

!*****************************************************************************80
!
!! FILE_NAME_EXT_GET determines the "extension" of a file name.
!
!  Discussion:
!
!    The "extension" of a filename is the string of characters
!    that appears after the LAST period in the name.  A file
!    with no period, or with a period as the last character
!    in the name, has a "null" extension.
!
!    Blanks are unusual in filenames.  This routine ignores all
!    trailing blanks, but will treat initial or internal blanks
!    as regular characters acceptable in a file name.
!
!  Example:
!
!    FILE_NAME   I  J
!
!    bob.for     4  7
!    N.B.C.D     6  7
!    Naomi.      6  6
!    Arthur     -1 -1 
!    .com        1  1
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, a file name to be examined.
!
!    Output, integer I, J, the indices of the first and last 
!    characters in the file extension.  
!
!    If no period occurs in FILE_NAME, then
!      I = J = -1;
!    Otherwise,
!      I is the position of the LAST period in FILE_NAME, and J is the
!      position of the last nonblank character following the period.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = * ) file_name
  integer i
  integer j
  integer s_index_last_c

  i = s_index_last_c ( file_name, '.' )

  if ( i == -1 ) then

    j = -1

  else

    j = len_trim ( file_name )

  end if

  return
end
subroutine file_name_ext_swap ( file_name, ext )

!*****************************************************************************80
!
!! FILE_NAME_EXT_SWAP replaces the current "extension" of a file name.
!
!  Discussion:
!
!    The "extension" of a filename is the string of characters
!    that appears after the LAST period in the name.  A file
!    with no period, or with a period as the last character
!    in the name, has a "null" extension.
!
!  Example:
!
!          Input           Output
!    ================     =========
!    FILE_NAME    EXT     FILE_NAME
!
!    bob.for      obj     bob.obj
!    bob.bob.bob  txt     bob.bob.txt
!    bob          yak     bob.yak
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) FILE_NAME, a file name.
!    On output, the extension of the file has been changed.
!
!    Input, character ( len = * ) EXT, the extension to be used on the output
!    copy of FILE_NAME, replacing the current extension if any.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = * ) ext
  character ( len = * ) file_name
  integer i
  integer j
  integer len_max
  integer len_name

  len_max = len ( file_name )
  len_name = len_trim ( file_name )

  call file_name_ext_get ( file_name, i, j )

  if ( i == -1 ) then

    if ( len_max < len_name + 1 ) then
      return
    end if

    len_name = len_name + 1
    file_name(len_name:len_name) = '.'
    i = len_name + 1

  else

    i = i + 1
    file_name(i:j) = ' '

  end if

  file_name(i:) = ext

  return
end
subroutine file_row_count ( file_name, nline )

!*****************************************************************************80
!
!! FILE_ROW_COUNT counts the number of rows in a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.
!
!    Blank lines and comment lines, which begin with '#', are not counted.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
!    Output, integer NLINE, the number of lines found in the file.
!    If the file could not be opened, then NLINE is returned as -1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = * ) file_name
  integer ios
  integer iunit
  character ( len = 255 ) line
  integer nline
  logical, parameter :: verbose = .false.

  nline = 0
!
!  Open the file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = file_name, status = 'old', form = 'formatted', &
    access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then

    nline = -1

    if ( verbose ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FILE_ROW_COUNT - Fatal error!'
      write ( *, '(a)' ) '  Could not open the file:'
      write ( *, '(4x,a)' ) '"' // trim ( file_name ) // '".'
    end if

    return

  end if
!
!  Count the lines.
!
  do

    read ( iunit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    if ( len_trim ( line ) == 0 ) then
      cycle
    end if

    if ( line(1:1) == '#' ) then
      cycle
    end if

    nline = nline + 1

  end do

  close ( unit = iunit )

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
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

  integer, parameter :: rk = kind ( 1.0D+00 )

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
function s_index_last_c ( s, c )

!*****************************************************************************80
!
!! S_INDEX_LAST_C finds the LAST occurrence of a given character.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 December 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be searched.
!
!    Input, character C, the character to search for.
!
!    Output, integer S_INDEX_LAST_C, the index in S where C occurs
!    last, or -1 if it does not occur.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character c
  integer i
  character ( len = * ) s
  integer s_len
  integer s_index_last_c

  if ( c == ' ' ) then
    s_len = len ( s )
  else
    s_len = len_trim ( s )
  end if

  do i = s_len, 1, -1

    if ( s(i:i) == c ) then
      s_index_last_c = i
      return
    end if

  end do

  s_index_last_c = -1

  return
end
subroutine s_to_r8 ( s, r, ierror, lchar )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 from a string.
!
!  Discussion:
!
!    This routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the real number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 spaces
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semicolon.
!
!    with most quantities optional.
!
!  Example:
!
!    S                 R
!
!    '1'               1.0
!    '     1   '       1.0
!    '1A'              1.0
!    '12,34,56'        12.0
!    '  34 7'          34.0
!    '-1E2ABCD'        -100.0
!    '-1X2ABCD'        -1.0
!    ' 2E-1'           0.2
!    '23.45'           23.45
!    '-4.2E+2'         -420.0
!    '17d2'            1700.0
!    '-14e-2'         -0.14
!    'e2'              100.0
!    '-12.73e-9.23'   -12.73 * 10.0^(-9.23)
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string containing the
!    data to be read.  Reading will begin at position 1 and
!    terminate at the end of the string, or when no more
!    characters can be read to form a legal real.  Blanks,
!    commas, or other nonnumeric data will, in particular,
!    cause the conversion to halt.
!
!    Output, real ( kind = rk ) R, the real value that was read from the string.
!
!    Output, integer IERROR, error flag.
!
!    0, no errors occurred.
!
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer LCHAR, the number of characters read from
!    the string to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character c
  logical ch_eqi
  integer ierror
  integer ihave
  integer isgn
  integer iterm
  integer jbot
  integer jsgn
  integer jtop
  integer lchar
  integer nchar
  integer ndig
  real ( kind = rk ) r
  real ( kind = rk ) rbot
  real ( kind = rk ) rexp
  real ( kind = rk ) rtop
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

  nchar = len_trim ( s )
  ierror = 0
  r = 0.0D+00
  lchar = - 1
  isgn = 1
  rtop = 0.0D+00
  rbot = 1.0D+00
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

  do

    lchar = lchar + 1
    c = s(lchar+1:lchar+1)
!
!  Blank or TAB character.
!
    if ( c == ' ' .or. c == TAB ) then

      if ( ihave == 2 ) then

      else if ( ihave == 6 .or. ihave == 7 ) then
        iterm = 1
      else if ( ihave > 1 ) then
        ihave = 11
      end if
!
!  Comma.
!
    else if ( c == ',' .or. c == ';' ) then

      if ( ihave /= 1 ) then
        iterm = 1
        ihave = 12
        lchar = lchar + 1
      end if
!
!  Minus sign.
!
    else if ( c == '-' ) then

      if ( ihave == 1 ) then
        ihave = 2
        isgn = - 1
      else if ( ihave == 6 ) then
        ihave = 7
        jsgn = - 1
      else
        iterm = 1
      end if
!
!  Plus sign.
!
    else if ( c == '+' ) then

      if ( ihave == 1 ) then
        ihave = 2
      else if ( ihave == 6 ) then
        ihave = 7
      else
        iterm = 1
      end if
!
!  Decimal point.
!
    else if ( c == '.' ) then

      if ( ihave < 4 ) then
        ihave = 4
      else if ( ihave >= 6 .and. ihave <= 8 ) then
        ihave = 9
      else
        iterm = 1
      end if
!
!  Exponent marker.
!
    else if ( ch_eqi ( c, 'E' ) .or. ch_eqi ( c, 'D' ) ) then

      if ( ihave < 6 ) then
        ihave = 6
      else
        iterm = 1
      end if
!
!  Digit.
!
    else if ( ihave < 11 .and. lge ( c, '0' ) .and. lle ( c, '9' ) ) then

      if ( ihave <= 2 ) then
        ihave = 3
      else if ( ihave == 4 ) then
        ihave = 5
      else if ( ihave == 6 .or. ihave == 7 ) then
        ihave = 8
      else if ( ihave == 9 ) then
        ihave = 10
      end if

      call ch_to_digit ( c, ndig )

      if ( ihave == 3 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = rk )
      else if ( ihave == 5 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = rk )
        rbot = 10.0D+00 * rbot
      else if ( ihave == 8 ) then
        jtop = 10 * jtop + ndig
      else if ( ihave == 10 ) then
        jtop = 10 * jtop + ndig
        jbot = 10 * jbot
      end if
!
!  Anything else is regarded as a terminator.
!
    else
      iterm = 1
    end if
!
!  If we haven't seen a terminator, and we haven't examined the
!  entire string, go get the next character.
!
    if ( iterm == 1 .or. lchar+1 >= nchar ) then
      exit
    end if

  end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LCHAR is equal to NCHAR.
!
  if ( iterm /= 1 .and. lchar+1 == nchar ) then
    lchar = nchar
  end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
  if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then

    ierror = ihave

    return
  end if
!
!  Number seems OK.  Form it.
!
  if ( jtop == 0 ) then
    rexp = 1.0D+00
  else

    if ( jbot == 1 ) then
      rexp = 10.0D+00**( jsgn * jtop )
    else
      rexp = jsgn * jtop
      rexp = rexp / jbot
      rexp = 10.0D+00**rexp
    end if

  end if

  r = isgn * rexp * rtop / rbot

  return
end
subroutine tr_to_xyzl ( n, x, y, z, list, lptr, lend, point_file_name )

!*****************************************************************************80
!
!! TR_TO_XYZL makes an XYZL file of Delaunay triangulation data.
!
!  Discussion:
!
!    The XYZL file, in combination with the original XYZ file, can be used
!    to draw the lines that form the Delaunay triangulation.
!
!    The XYZL file is simply a list of point indices from the XYZ file
!    which are to be connected.  The termination of a line is
!    indicated by an index value of -1.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of nodes.
!
!    Input, real ( kind = rk ) X(N), Y(N), Z(N), the coordinates of the nodes.
!
!    Input, integer LIST(6*(N-2)), LPTR(6*(N-2)), LEND(N),
!    information defining the triangulation, created by TRMESH.
!
!    Input, character ( len = * ) POINT_FILE_NAME, the name of the
!    file containing the point coordinates.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  character ( len = 80 ) delaunay_file_name
  logical, parameter :: header = .false.
  integer iunit
  integer lend(n)
  integer list(6*(n-2))
  integer lp
  integer lpl
  integer lptr(6*(n-2))
  integer node1
  integer node2
  character ( len = *  ) point_file_name
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)
  real ( kind = rk ) z(n)

  delaunay_file_name = point_file_name

  call file_name_ext_swap ( delaunay_file_name, 'xyzl' )

  call get_unit ( iunit )

  open ( unit = iunit, file = delaunay_file_name )
!
!  Write header.
!
  if ( header ) then
    write ( iunit, '(a)' ) '# "' // trim ( delaunay_file_name ) // '".'
    write ( iunit, '(a)' ) '# created by SPHERE_DELAUNAY.F90.'
    write ( iunit, '(a)' ) '#'
    write ( iunit, '(a)' ) &
      '# This file contains the indices of points to be connected'
    write ( iunit, '(a)' ) '# by lines, to form a Delaunay triangulation.'
    write ( iunit, '(a)' ) &
      '# The points are in the file "' // trim ( point_file_name ) // '".'
    write ( iunit, '(a)' ) '#'
    write ( iunit, '(a)' ) &
      '# This file lists the indices of points to be connected'
    write ( iunit, '(a)' ) '# by lines.  A line may include several points.'
    write ( iunit, '(a)' ) '# Each line is terminated by an index of -1.'
    write ( iunit, '(a)' ) '#'
  end if
!
!  List all the line segments that emanate from a point.
!  This will involve listing some segments twice.
!
  do node1 = 1, n

    lpl = lend(node1)
    lp = lpl

    do

      lp = lptr(lp)
      node2 = abs ( list(lp) )

      write ( iunit, '(i8)' ) node1
      write ( iunit, '(i8)' ) node2
      write ( iunit, '(i8)' ) -1

      if ( lp == lpl ) then
        exit
      end if

    end do

  end do

  close ( unit = iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TR_TO_XYZL:'
  write ( *, '(a)' ) '  Wrote the Delaunay XYZL file "' // &
    trim ( delaunay_file_name ) // '".'

  return
end
subroutine word_next_read ( s, word, done )

!*****************************************************************************80
!
!! WORD_NEXT_READ "reads" words from a string, one at a time.
!
!  Discussion:
!
!    The following characters are considered to be a single word,
!    whether surrounded by spaces or not:
!
!      " ( ) { } [ ]
!
!    Also, if there is a trailing comma on the word, it is stripped off.
!    This is to facilitate the reading of lists.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string, presumably containing words
!    separated by spaces.
!
!    Output, character ( len = * ) WORD.
!
!    If DONE is FALSE, then WORD contains the "next" word read.
!    If DONE is TRUE, then WORD is blank, because there was no more to read.
!
!    Input/output, logical DONE.
!
!    On input with a fresh string, set DONE to TRUE.
!
!    On output, the routine sets DONE:
!      FALSE if another word was read,
!      TRUE if no more words could be read.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  logical done
  integer ilo
  integer, save :: lenc = 0
  integer, save :: next = 1
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )
  character ( len = * ) word
!
!  We "remember" LENC and NEXT from the previous call.
!
!  An input value of DONE = TRUE signals a new line of text to examine.
!
  if ( done ) then

    next = 1
    done = .false.
    lenc = len_trim ( s )

    if ( lenc <= 0 ) then
      done = .true.
      word = ' '
      return
    end if

  end if
!
!  Beginning at index NEXT, search the string for the next nonblank,
!  which signals the beginning of a word.
!
  ilo = next
!
!  ...S(NEXT:) is blank.  Return with WORD = ' ' and DONE = TRUE.
!
  do

    if ( ilo > lenc ) then
      word = ' '
      done = .true.
      next = lenc + 1
      return
    end if
!
!  If the current character is blank, skip to the next one.
!
    if ( s(ilo:ilo) /= ' ' .and. s(ilo:ilo) /= TAB ) then
      exit
    end if

    ilo = ilo + 1

  end do
!
!  ILO is the index of the next nonblank character in the string.
!
!  If this initial nonblank is a special character,
!  then that's the whole word as far as we're concerned,
!  so return immediately.
!
  if ( s(ilo:ilo) == '"' .or. &
       s(ilo:ilo) == '(' .or. &
       s(ilo:ilo) == ')' .or. &
       s(ilo:ilo) == '{' .or. &
       s(ilo:ilo) == '}' .or. &
       s(ilo:ilo) == '[' .or. &
       s(ilo:ilo) == ']' ) then

    word = s(ilo:ilo)
    next = ilo + 1
    return

  end if
!
!  Now search for the last contiguous character that is not a
!  blank, TAB, or special character.
!
  next = ilo + 1

  do while ( next <= lenc )

    if ( s(next:next) == ' ' ) then
      exit
    else if ( s(next:next) == TAB ) then
      exit
    else if ( s(next:next) == '"' ) then
      exit
    else if ( s(next:next) == '(' ) then
      exit
    else if ( s(next:next) == ')' ) then
      exit
    else if ( s(next:next) == '{' ) then
      exit
    else if ( s(next:next) == '}' ) then
      exit
    else if ( s(next:next) == '[' ) then
      exit
    else if ( s(next:next) == ']' ) then
      exit
    end if

    next = next + 1

  end do
!
!  Ignore a trailing comma.
!
  if ( s(next-1:next-1) == ',' ) then
    word = s(ilo:next-2)
  else
    word = s(ilo:next-1)
  end if

  return
end
subroutine xyz_read ( point_file_name, n, x, y, z, ierror )

!*****************************************************************************80
!
!! XYZ_READ reads graphics information from an XYZ file.
!
!  Discussion:
!
!    Comment lines begin with '#";
!    The XYZ coordinates of a point are written on a single line.
!
!  Example:
!
!     # cube.xyz
!     #
!     0 0 0
!     0 0 1
!     0 1 0
!     0 1 1
!     1 0 0
!     1 0 1
!     1 1 0
!     1 1 1
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) XYZ_FILE_NAME, the name of the input file.
!
!    Input, integer N, the number of points.
!
!    Output, real ( kind = rk ) X(N), Y(N), Z(N), the point coordinates.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    nonzero, an error occurred.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  logical done
  integer i
  integer ierror
  integer ios
  integer iunit
  integer lchar
  character ( len = 256 ) line
  integer n2
  real ( kind = rk ) temp(3)
  integer text_num
  character ( len = 100 ) word
  real ( kind = rk ) x(n)
  character ( len = * ) point_file_name
  real ( kind = rk ) y(n)
  real ( kind = rk ) z(n)

  n2 = 0
  ierror = 0
  word = ' '
  text_num = 0

  call get_unit ( iunit )

  open ( unit = iunit, file = point_file_name, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'XYZ_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file.'
    ierror = 1
    stop
  end if
!
!  Read a line of text from the file.
!
  do

    read ( iunit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    text_num = text_num + 1
!
!  If this line begins with '#' , then it's a comment.  Read a new line.
!
    if ( line(1:1) == '#' ) then
      cycle
    end if
!
!  If this line is blank, then record that information.
!
    if ( len_trim ( line ) == 0 ) then
      cycle
    end if
!
!  This line records a node's coordinates.
!
    n2 = n2 + 1

    if ( n2 <= n ) then

      done = .true.

      do i = 1, 3

        call word_next_read ( line, word, done )

        call s_to_r8 ( word, temp(i), ierror, lchar )

        if ( ierror /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'XYZ_READ - Fatal error!'
          write ( *, '(a,i8)' ) '  S_TO_R8 returned IERROR = ', ierror
          write ( *, '(a,i8)' ) '  Reading (X,Y,Z) component ', i
          exit
        end if

      end do

      if ( ierror /= 0 ) then
        exit
      end if

      x(n2) = temp(1)
      y(n2) = temp(2)
      z(n2) = temp(3)

    end if

  end do

  close ( unit = iunit )
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) 'XYZ_READ:'
  write ( *, '(a,i8,a)' ) '  Read ', text_num, ' text lines from ' &
    // trim ( point_file_name )
  write ( *, '(a,i8,a)' ) '  Read ', n2, ' sets of (X,Y,Z) coordinates.'

  return
end
