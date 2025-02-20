program main

!*****************************************************************************80
!
!! CLENSHAW_CURTIS_RULE computes a Clenshaw Curtis quadrature rule.
!
!  Discussion:
!
!    The user specifies:
!    * the ORDER (number of points) in the rule
!    * the root name of the output files.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 February 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  integer arg_num
  real ( kind = rk ) b
  character ( len = 255 ) filename
  integer iarg
  integer iargc
  integer ierror
  integer last
  integer order
  real ( kind = rk ) r(2)
  character ( len = 255 ) string
  real ( kind = rk ), allocatable, dimension ( : ) :: w
  real ( kind = rk ), allocatable, dimension ( : ) :: x

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CLENSHAW_CURTIS_RULE'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Compute a Clenshaw Curtis rule for approximating'
  write ( *, '(a)' ) '    Integral ( A <= x <= B ) f(x) dx'
  write ( *, '(a)' ) '  of order ORDER.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The user specifies ORDER, A, B and FILENAME.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  ORDER is the number of points;'
  write ( *, '(a)' ) '  A is the left endpoint;'
  write ( *, '(a)' ) '  B is the right endpoint;'
  write ( *, '(a)' ) '  FILENAME is used to generate 3 files:'
  write ( *, '(a)' ) '    filename_w.txt - the weight file'
  write ( *, '(a)' ) '    filename_x.txt - the abscissa file.'
  write ( *, '(a)' ) '    filename_r.txt - the region file.'
!
!  Get the number of command line arguments.
!
  arg_num = iargc ( )
!
!  Get ORDER.
!
  if ( 1 <= arg_num ) then
    iarg = 1
    call getarg ( iarg, string )
    call s_to_i4 ( string, order, ierror, last )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter the rule order ORDER:'
    read ( *, * ) order
  end if
!
!  Get A.
!
  if ( 2 <= arg_num ) then
    iarg = 2
    call getarg ( iarg, string )
    call s_to_r8 ( string, a, ierror, last )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter the left endpoint A:'
    read ( *, * ) a
  end if
!
!  Get B.
!
  if ( 3 <= arg_num ) then
    iarg = 3
    call getarg ( iarg, string )
    call s_to_r8 ( string, b, ierror, last )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter the right endpoint B:'
    read ( *, * ) b
  end if
!
!  Get FILENAME.
!
  if ( 4 <= arg_num ) then
    iarg = 4
    call getarg ( iarg, filename )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) &
      '  Enter FILENAME, the "root name" of the quadrature files).'
    read ( *, '(a)' ) filename
  end if
!
!  Input summary.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  ORDER = ', order
  write ( *, '(a,g14.6)' ) '  A = ', a
  write ( *, '(a,g14.6)' ) '  B = ', b
  write ( *, '(a)'    ) '  FILENAME = "' // trim ( filename ) // '".'
!
!  Construct the rule.
!
  r(1) = a
  r(2) = b

  allocate ( w(order) )
  allocate ( x(order) )

  call clenshaw_curtis_compute ( order, x, w )
!
!  Rescale the rule.
!
  call rescale ( a, b, order, x, w )
!
!  Write the rule.
!
  call rule_write ( order, x, w, r, filename )
!
!  Free memory.
!
  deallocate ( w )
  deallocate ( x )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CLENSHAW_CURTIS_RULE:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine ch_cap ( ch )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
!
!  Discussion:
!
!    Instead of CHAR and ICHAR, we now use the ACHAR and IACHAR functions, 
!    which guarantee the ASCII collating sequence.
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
!    Input/output, character CH, the character to capitalize.
!
  implicit none

  character ch
  integer itemp

  itemp = iachar ( ch )
 
  if ( 97 <= itemp .and. itemp <= 122 ) then
    ch = achar ( itemp - 32 )
  end if
 
  return
end
function ch_eqi ( c1, c2 )

!*****************************************************************************80
!
!! CH_EQI is a case insensitive comparison of two characters for equality.  
!
!  Discussion:
!
!    CH_EQI ( 'A', 'a' ) is TRUE.
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
subroutine ch_to_digit ( ch, digit )

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the integer value of a base 10 digit.
!
!  Discussion:
!
!    Instead of ICHAR, we now use the IACHAR function, which
!    guarantees the ASCII collating sequence.
!
!  Example:
!
!     CH  DIGIT
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
!    Input, character CH, the decimal digit, '0' through '9' or blank
!    are legal. 
!
!    Output, integer DIGIT, the corresponding integer value.  
!    If CH was 'illegal', then DIGIT is -1.
!
  implicit none

  character ch
  integer digit

  if ( lle ( '0', ch ) .and. lle ( ch, '9' ) ) then
 
    digit = iachar ( ch ) - 48
 
  else if ( ch == ' ' ) then
 
    digit = 0
 
  else

    digit = -1

  end if
 
  return
end
subroutine clenshaw_curtis_compute ( order, x, w )

!*****************************************************************************80
!
!! CLENSHAW_CURTIS_COMPUTE computes a Clenshaw Curtis quadrature rule.
!
!  Discussion:
!
!    Our convention is that the abscissas are numbered from left to right.
!
!    The rule is defined on [-1,1].
!
!    The integral to approximate:
!
!      Integral ( -1 <= X <= 1 ) F(X) dX
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 February 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ORDER, the order of the rule.
!    1 <= ORDER.
!
!    Output, real ( kind = rk ) X(ORDER), the abscissas.
!
!    Output, real ( kind = rk ) W(ORDER), the weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer order

  real ( kind = rk ) b
  integer i
  integer j
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) theta
  real ( kind = rk ) w(order)
  real ( kind = rk ) x(order)

  if ( order < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CLENSHAW_CURTIS_COMPUTE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of ORDER = ', order
    stop
  end if

  if ( order == 1 ) then
    x(1) = 0.0D+00
    w(1) = 2.0D+00
    return
  end if

  do i = 1, order
    x(i) = cos ( real ( order - i, kind = rk ) * r8_pi &
               / real ( order - 1, kind = rk ) )
  end do

  x(1) = -1.0D+00
  if ( mod ( order, 2 ) == 1 ) then
    x((order+1)/2) = 0.0D+00
  end if
  x(order) = +1.0D+00

  do i = 1, order

    theta = real ( i - 1, kind = rk ) * r8_pi &
          / real ( order - 1, kind = rk )

    w(i) = 1.0D+00

    do j = 1, ( order - 1 ) / 2

      if ( 2 * j == ( order - 1 ) ) then
        b = 1.0D+00
      else
        b = 2.0D+00
      end if

      w(i) = w(i) - b * cos ( 2.0D+00 * real ( j, kind = rk ) * theta ) &
           / real ( 4 * j * j - 1, kind = rk )

    end do

  end do

  w(1)         =           w(1)         / real ( order - 1, kind = rk )
  w(2:order-1) = 2.0D+00 * w(2:order-1) / real ( order - 1, kind = rk )
  w(order)     =           w(order)     / real ( order - 1, kind = rk )

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
!    15 January 2008
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
subroutine r8mat_write ( output_filename, m, n, table )

!*****************************************************************************80
!
!! R8MAT_WRITE writes an R8MAT file.
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
!    Input, real ( kind = rk ) TABLE(M,N), the table data.
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
subroutine rescale ( a, b, n, x, w )

!*****************************************************************************80
!
!! RESCALE rescales a quadrature rule from [-1,+1] to [A,B].
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, the endpoints of the new interval.
!
!    Input, integer N, the order.
!
!    Input/output, real ( kind = rk ) X(N), on input, the abscissas for [-1,+1].
!    On output, the abscissas for [A,B].
!
!    Input/output, real ( kind = rk ) W(N), on input, the weights for [-1,+1].
!    On output, the weights for [A,B].
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n)

  x(1:n) = ( ( a + b ) + ( b - a ) * x(1:n) ) / 2.0D+00
  w(1:n) = ( b - a ) * w(1:n) / 2.0D+00

  return
end
subroutine rule_write ( order, x, w, r, filename )

!*****************************************************************************80
!
!! RULE_WRITE writes a quadrature rule to a file.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ORDER, the order of the rule.
!
!    Input, real ( kind = rk ) X(ORDER), the abscissas.
!
!    Input, real ( kind = rk ) W(ORDER), the weights.
!
!    Input, real ( kind = rk ) R(2), defines the region.
!
!    Input, character ( len = * ) FILENAME, specifies the output.
!    'filename_w.txt', 'filename_x.txt', 'filename_r.txt' defining weights,
!    abscissas, and region.
! 
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer order

  character ( len = * ) filename
  character ( len = 255 ) filename_r
  character ( len = 255 ) filename_w
  character ( len = 255 ) filename_x
  real ( kind = rk ) r(2)
  real ( kind = rk ) w(order)
  real ( kind = rk ) x(order)

  filename_w = trim ( filename ) // '_w.txt'
  filename_x = trim ( filename ) // '_x.txt'
  filename_r = trim ( filename ) // '_r.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Creating quadrature files.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  "Root" file name is   "' // trim ( filename ) // '".'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Weight file will be   "' // trim ( filename_w ) // '".'
  write ( *, '(a)' ) '  Abscissa file will be "' // trim ( filename_x ) // '".'
  write ( *, '(a)' ) '  Region file will be   "' // trim ( filename_r ) // '".'
            
  call r8mat_write ( filename_w, 1, order, w )
  call r8mat_write ( filename_x, 1, order, x )
  call r8mat_write ( filename_r, 1, 2,     r )

  return
end
subroutine s_to_i4 ( s, ival, ierror, length )

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
!    15 January 2008
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
!    If the string is blank, then IVAL will be returned 0.
!
!    Output, integer IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer LENGTH, the number of characters of S 
!    used to make IVAL.
!
  implicit none

  character c
  integer i
  integer ierror
  integer isgn
  integer istate
  integer ival
  integer length
  character ( len = * ) s

  ierror = 0
  istate = 0
  isgn = 1
  ival = 0

  do i = 1, len_trim ( s )

    c = s(i:i)
!
!  Haven't read anything.
!
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
        return
      end if
!
!  Have read the sign, expecting digits.
!
    else if ( istate == 1 ) then

      if ( c == ' ' ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read at least one digit, expecting more.
!
    else if ( istate == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        ival = 10 * ival + ichar ( c ) - ichar ( '0' )
      else
        ival = isgn * ival
        length = i - 1
        return
      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( istate == 2 ) then
    ival = isgn * ival
    length = len_trim ( s )
  else
    ierror = 1
    length = 0
  end if

  return
end
subroutine s_to_r8 ( s, dval, ierror, length )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 value from a string.
!
!  Discussion:
!
!    An "R8" value is simply a real number to be stored as a
!    variable of type "real ( kind = rk )".
!
!    The routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 blanks
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semicolon,
!
!    with most quantities optional.
!
!  Example:
!
!    S                 DVAL
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
!    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 January 2009
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
!    Output, real ( kind = rk ) DVAL, the value read from the string.
!
!    Output, integer IERROR, error flag.
!    0, no errors occurred.
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer LENGTH, the number of characters read
!    to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character c
  logical ch_eqi
  real ( kind = rk ) dval
  integer ierror
  integer ihave
  integer isgn
  integer iterm
  integer jbot
  integer jsgn
  integer jtop
  integer length
  integer ndig
  real ( kind = rk ) rbot
  real ( kind = rk ) rexp
  real ( kind = rk ) rtop
  character ( len = * ) s
  integer s_length
  character :: TAB = achar ( 9 )

  s_length = len_trim ( s )

  ierror = 0
  dval = 0.0D+00
  length = -1
  isgn = 1
  rtop = 0
  rbot = 1
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

  do

    length = length + 1

    if ( s_length < length+1 ) then
      exit
    end if

    c = s(length+1:length+1)
!
!  Blank character.
!
    if ( c == ' ' .or. c == TAB ) then

      if ( ihave == 2 ) then

      else if ( ihave == 6 .or. ihave == 7 ) then
        iterm = 1
      else if ( 1 < ihave ) then
        ihave = 11
      end if
!
!  Comma.
!
    else if ( c == ',' .or. c == ';' ) then

      if ( ihave /= 1 ) then
        iterm = 1
        ihave = 12
        length = length + 1
      end if
!
!  Minus sign.
!
    else if ( c == '-' ) then

      if ( ihave == 1 ) then
        ihave = 2
        isgn = -1
      else if ( ihave == 6 ) then
        ihave = 7
        jsgn = -1
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
      else if ( 6 <= ihave .and. ihave <= 8 ) then
        ihave = 9
      else
        iterm = 1
      end if
!
!  Scientific notation exponent marker.
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
    else if (  ihave < 11 .and. lle ( '0', c ) .and. lle ( c, '9' ) ) then

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
    if ( iterm == 1 ) then
      exit
    end if

  end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LENGTH is equal to S_LENGTH.
!
  if ( iterm /= 1 .and. length+1 == s_length ) then
    length = s_length
  end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
  if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then
    ierror = ihave
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'S_TO_R8 - Serious error!'
    write ( *, '(a)' ) '  Illegal or nonnumeric input:'
    write ( *, '(a)' ) '    ' // trim ( s )
    return
  end if
!
!  Number seems OK.  Form it.
!
  if ( jtop == 0 ) then
    rexp = 1.0D+00
  else
    if ( jbot == 1 ) then
      rexp = 10.0D+00 ** ( jsgn * jtop )
    else
      rexp = 10.0D+00 ** ( real ( jsgn * jtop, kind = rk ) &
        / real ( jbot, kind = rk ) )
    end if
  end if

  dval = real ( isgn, kind = rk ) * rexp * rtop / rbot

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 January 2008
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
  integer d
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = 10 ) time
  integer values(8)
  integer y
  character ( len = 5 ) zone

  call date_and_time ( date, time, zone, values )

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

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
