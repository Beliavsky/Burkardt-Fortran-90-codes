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

  logical ch_eqi
  character c1
  character c1_cap
  character c2
  character c2_cap

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
!! CH_TO_DIGIT returns the integer value of a base 10 digit.
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
!    Output, integer DIGIT, the corresponding integer value.
!    If C was 'illegal', then DIGIT is -1.
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
subroutine i4mat_copy ( m, n, a1, a2 )

!*****************************************************************************80
!
!! I4MAT_COPY copies an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer A1(M,N), the matrix to be copied.
!
!    Output, integer A2(M,N), the copied matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer a1(m,n)
  integer a2(m,n)

  a2(1:m,1:n) = a1(1:m,1:n)

  return
end
subroutine i4mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! I4MAT_TRANSPOSE_PRINT prints an I4MAT, transposed.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer A(M,N), an M by N matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer a(m,n)
  character ( len = * ) title

  call i4mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine i4mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! I4MAT_TRANSPOSE_PRINT_SOME prints some of the transpose of an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4's.
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
!    Input, integer A(M,N), an M by N matrix to be printed.
!
!    Input, integer ILO, JLO, the first row and column to print.
!
!    Input, integer IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: incx = 10
  integer m
  integer n

  integer a(m,n)
  character ( len = 8 ) ctemp(incx)
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
      write ( ctemp(i2), '(i8)' ) i
    end do

    write ( *, '(''  Row '',10a8)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc

        i = i2lo - 1 + i2

        write ( ctemp(i2), '(i8)' ) a(i,j)

      end do

      write ( *, '(i5,a,10a8)' ) j, ':', ( ctemp(i), i = 1, inc )

    end do

  end do

  return
end
subroutine r8mat_copy ( m, n, a, b )

!*****************************************************************************80
!
!! R8MAT_COPY copies an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> (I+J*M).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 July 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the order of the matrix.
!
!    Input, real ( kind = rk ) A(M,N), the matrix to be copied.
!
!    Output, real ( kind = rk ) B(M,N), a copy of the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) b(m,n)

  b(1:m,1:n) = a(1:m,1:n)

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
!    28 June 2000
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

  integer, parameter :: rk = kind ( 1.0D+00 )

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
subroutine s_to_i4vec ( s, n, ivec, ierror )

!*****************************************************************************80
!
!! S_TO_I4VEC reads an I4VEC from a string.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be read.
!
!    Input, integer N, the number of values expected.
!
!    Output, integer IVEC(N), the values read from the string.
!
!    Output, integer IERROR, error flag.
!    0, no errors occurred.
!    -K, could not read data for entries -K through N.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer ierror
  integer ilo
  integer ivec(n)
  integer length
  character ( len = * ) s

  i = 0
  ierror = 0
  ilo = 1

  do while ( i < n )

    i = i + 1

    call s_to_i4 ( s(ilo:), ivec(i), ierror, length )

    if ( ierror /= 0 ) then
      ierror = -i
      exit
    end if

    ilo = ilo + length

  end do

  return
end
subroutine s_to_r8 ( s, dval, ierror, length )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 from a string.
!
!  Discussion:
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
!    '-12.73e-9.23'   -12.73 * 10.0^(-9.23)
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 September 2004
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
  integer nchar
  integer ndig
  real ( kind = rk ) rbot
  real ( kind = rk ) rexp
  real ( kind = rk ) rtop
  character ( len = * ) s

  nchar = len_trim ( s )

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

    if ( nchar < length+1 ) then
      exit
    end if

    c = s(length+1:length+1)
!
!  Blank character.
!
    if ( c == ' ' ) then

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
!  entire string, then we're done, and LENGTH is equal to NCHAR.
!
  if ( iterm /= 1 .and. length + 1 == nchar ) then
    length = nchar
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
subroutine s_to_r8vec ( s, n, rvec, ierror )

!*****************************************************************************80
!
!! S_TO_R8VEC reads an R8VEC from a string.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be read.
!
!    Input, integer N, the number of values expected.
!
!    Output, real ( kind = rk ) RVEC(N), the values read from the string.
!
!    Output, integer IERROR, error flag.
!    0, no errors occurred.
!    -K, could not read data for entries -K through N.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer ierror
  integer ilo
  integer lchar
  real ( kind = rk ) rvec(n)
  character ( len = * ) s

  i = 0
  ierror = 0
  ilo = 1

  do while ( i < n )

    i = i + 1

    call s_to_r8 ( s(ilo:), rvec(i), ierror, lchar )

    if ( ierror /= 0 ) then
      ierror = -i
      exit
    end if

    ilo = ilo + lchar

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
subroutine triangle_node_data_example ( node_num, node_dim, node_att_num, &
  node_marker_num, node_coord, node_att, node_marker )

!*****************************************************************************80
!
!! TRIANGLE_NODE_DATA_EXAMPLE returns the node information for the example.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 October 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer NODE_DIM, the spatial dimension.
!
!    Input, integer NODE_ATT_NUM, number of node attributes 
!    listed on each node record.
!
!    Input, integer NODE_MARKER_NUM, 1 if every node record 
!    includes a final boundary marker value.
!
!    Output, real ( kind = rk ) NODE_COORD(NODE_DIM,NODE_NUM), the nodal 
!    coordinates.
!
!    Output, real ( kind = rk ) NODE_ATT(NODE_ATT_NUM,NODE_NUM), the nodal 
!    attributes.
!
!    Output, integer NODE_MARKER(NODE_MARKER_NUM,NODE_NUM), 
!    the node markers.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer node_att_num
  integer node_dim
  integer node_marker_num
  integer node_num

  real ( kind = rk ) node_att(node_att_num,node_num)
  real ( kind = rk ) node_coord(node_dim,node_num)
  real ( kind = rk ), dimension ( 2, 21 ) :: node_coord_save = reshape ( (/ &
    0.0, 0.0, &
    1.0, 0.0, &
    2.0, 0.0, &
    3.0, 0.0, &
    4.0, 0.0, &
    0.0, 1.0, &
    1.0, 1.0, &
    2.0, 1.0, &
    3.0, 1.0, &
    4.0, 1.0, &
    0.0, 2.0, &
    1.0, 2.0, &
    2.0, 2.0, &
    3.0, 2.0, &
    4.0, 2.0, &
    0.0, 3.0, &
    1.0, 3.0, &
    2.0, 3.0, &
    0.0, 4.0, &
    1.0, 4.0, &
    2.0, 4.0 /), (/ 2, 21 /) )
  integer node_marker(node_marker_num,node_num)
  integer, dimension ( 1, 21 ) :: node_marker_save = reshape ( (/ &
    1, 1, 1, 1, 1, 1, 0, 0, 0, 1, &
    1, 0, 0, 1, 1, 1, 0, 1, 1, 1, &
    1 /), (/ 1, 21 /) )

  call r8mat_copy ( node_dim, node_num, node_coord_save, node_coord )

  call i4mat_copy ( node_marker_num, node_num, node_marker_save, &
    node_marker )

  return
end
subroutine triangle_node_data_read ( node_filename, node_num, node_dim, &
  node_att_num, node_marker_num, node_coord, node_att, node_marker )

!*****************************************************************************80
!
!! TRIANGLE_NODE_DATA_READ reads the data from a node file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 October 2014
!
!  Author:
!
!    John Burkardt
!
! Parameters:
!
!   Input, character ( len = * ) NODE_FILENAME, the name of the node file.
!
!   Input, integer NODE_NUM, the number of nodes.
!
!   Input, integer NODE_DIM, the spatial dimension.
!
!   Input, integer NODE_ATT_NUM, number of node attributes 
!   listed on each node record.
!
!   Input, integer NODE_MARKER_NUM, 1 if every node record 
!   includes a final boundary marker value.
!
!   Output, real ( kind = rk ) NODE_COORD(NODE_DIM,NODE_NUM), the nodal 
!   coordinates.
!
!   Output, real ( kind = rk ) NODE_ATT(NODE_ATT_NUM,NODE_NUM), the nodal 
!   attributes.
!
!   Output, integer NODE_MARKER(NODE_MARKER_NUM,NODE_NUM), the 
!   node markers.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer node_att_num
  integer node_dim
  integer node_marker_num
  integer node_num

  integer i
  integer ierror
  integer input
  integer ios
  integer ival
  integer length
  integer node
  real ( kind = rk ) node_att(node_att_num,node_num)
  real ( kind = rk ) node_coord(node_dim,node_num)
  character ( len = * ) node_filename
  integer node_marker(node_marker_num,node_num)
  character ( len = 255 ) text
  real ( kind = rk ) value

  node = 0

  call get_unit ( input )

  open ( unit = input, file = node_filename, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'TRIANGLE_NODE_DATA_READ - Fatal error!'
    write ( *, '(a)' ) '  Unable to open file.'
    stop 1
  end if

  do

    read ( input, '(a)', iostat = ios ) text

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TRIANGLE_NODE_DATA_READ - Fatal error!'
      write ( *, '(a)' ) '  Unexpected end of file while reading.'
      stop 1
    end if

    if ( len_trim ( text ) == 0 ) then
      cycle
    end if

    if ( text(1:1) == '#' ) then
      cycle
    end if
!
!  Ignore the dimension line.
!
    if ( node == 0 ) then

    else

      call s_to_i4 ( text, ival, ierror, length )
      text = text(length+1:)

      do i = 1, node_dim
        call s_to_r8 ( text, value, ierror, length )
        text = text(length+1:)
        node_coord(i,node) = value
      end do

      do i = 1, node_att_num
        call s_to_r8 ( text, value, ierror, length )
        text = text(length+1:)
        node_att(i,node) = value;
      end do

      do i = 1, node_marker_num
        call s_to_i4 ( text, ival, ierror, length )
        text = text(length+1:)
        node_marker(i,node) = ival
      end do

    end if

    node = node + 1

    if ( node_num < node ) then
      exit
    end if

  end do

  close ( unit = input )

  return
end
subroutine triangle_node_size_example ( node_num, node_dim, node_att_num, &
  node_marker_num )

!*****************************************************************************80
!
!! TRIANGLE_NODE_SIZE_EXAMPLE returns the sizes of node information for the example.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 October 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer NODE_NUM, the number of nodes.
!
!    Output, integer NODE_DIM, the spatial dimension.
!
!    Output, integer NODE_ATT_NUM, number of node attributes 
!    listed on each node record.
!
!    Output, integer NODE_MARKER_NUM, 1 if every node record 
!    includes a final boundary marker value.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer node_att_num
  integer node_dim
  integer node_marker_num
  integer node_num

  node_num = 21
  node_dim = 2
  node_att_num = 0
  node_marker_num = 1

  return
end
subroutine triangle_node_size_read ( node_filename, node_num, node_dim, &
  node_att_num, node_marker_num )

!*****************************************************************************80
!
!! NODE_SIZE_READ reads the header information from a node file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 October 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) NODE_FILENAME, the name of the node file.
!
!    Output, integer NODE_NUM, the number of nodes.
!
!    Output, integer NODE_DIM, the spatial dimension.
!
!    Output, integer NODE_ATT_NUM, number of node attributes 
!    listed on each node record.
!
!    Output, integer NODE_MARKER_NUM, 1 if every node record 
!    includes a final boundary marker value.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ierror
  integer input
  integer ios
  integer length
  integer node_att_num
  integer node_dim
  character ( len = * ) node_filename
  integer node_marker_num
  integer node_num
  character ( len = 255 ) text

  call get_unit ( input )

  open ( unit = input, file = node_filename, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'TRIANGLE_NODE_SIZE_READ - Fatal error!'
    write ( *, '(a)' ) '  Unable to open file.'
    stop 1
  end if

  do

    read ( input, '(a)', iostat = ios ) text

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TRIANGLE_NODE_SIZE_READ - Fatal error!'
      write ( *, '(a)' ) '  Unexpected end of file while reading.'
      stop 1
    end if

    if ( len_trim ( text ) == 0 ) then
      cycle
    end if

    if ( text(1:1) == '#' ) then
      cycle
    end if

    call s_to_i4 ( text, node_num, ierror, length )
    text = text(length+1:)
    call s_to_i4 ( text, node_dim, ierror, length )
    text = text(length+1:)
    call s_to_i4 ( text, node_att_num, ierror, length )
    text = text(length+1:)
    call s_to_i4 ( text, node_marker_num, ierror, length )
    text = text(length+1:)

    exit

  end do

  close ( unit = input )

  return
end
subroutine triangle_node_write ( node_filename, node_num, node_dim, &
  node_att_num, node_marker_num, node_coord, node_att, node_marker )

!*****************************************************************************80
!
!! TRIANGLE_NODE_WRITE writes a TRIANGLE ".node" file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 October 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) NODE_FILENAME, the name of the node file.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer NODE_DIM, the spatial dimension.
!
!    Input, integer NODE_ATT_NUM, number of node attributes
!    listed on each node record.
!
!    Input, integer NODE_MARKER_NUM, 1 if every node record 
!    includes a final boundary marker value.
!
!    Output, real ( kind = rk ) NODE_COORD(NODE_DIM*NODE_NUM), the nodal 
!    coordinates.
!
!    Output, real ( kind = rk ) NODE_ATT(NODE_ATT_NUM*NODE_NUM), the nodal
!    attributes.
!
!    Output, integer NODE_MARKER(NODE_MARKER_NUM,NODE_NUM), 
!    the node markers.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer node_att_num
  integer node_dim
  integer node_marker_num
  integer node_num

  integer att
  integer dim
  integer ios
  integer marker
  integer node
  real ( kind = rk ) node_att(node_att_num,node_num)
  real ( kind = rk ) node_coord(node_dim,node_num)
  character ( len = * ) node_filename
  integer node_marker(node_marker_num,node_num)
  integer output

  call get_unit ( output )

  open ( unit = output, file = node_filename, status = 'replace', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'TRIANGLE_NODE_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Unable to open file.'
    stop 1
  end if

  write ( output, '(2x,i4,2x,i4,2x,i4,2x,i4)' ) &
    node_num, node_dim, node_att_num, node_marker_num 

  do node = 1, node_num

    write ( output, '(i4)', advance = 'no' ) node

    do dim = 1, node_dim
      write ( output, '(2x,g14.6)', advance = 'no' ) node_coord(dim,node)
    end do

    do att = 1, node_att_num
      write ( output, '(2x,g14.6)', advance = 'no' ) node_att(att,node)
    end do

    do marker = 1, node_marker_num
      write ( output, '(2x,i4)', advance = 'no' ) node_marker(marker,node)
    end do

    write ( output, '(a)', advance = 'yes' ) 

  end do

  close ( unit = output )

  return
end
subroutine triangle_element_data_example ( element_num, element_order, &
  element_att_num, element_node, element_att )

!*****************************************************************************80
!
!! TRIANGLE_ELEMENT_DATA_EXAMPLE returns the element information for the example.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 October 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer ELEMENT_ORDER, the order of the elements.
!
!    Input, integer ELEMENT_ATT_NUM, the number of element 
!    attributes.
!
!    Output, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM), 
!    the indices of the nodes that make up each element.
!
!    Output, real ( kind = rk ) ELEMENT_ATT(ELEMENT_ATT_NUM,ELEMENT_NUM), 
!    the attributes of each element.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer element_att_num
  integer element_num
  integer element_order

  real ( kind = rk ) element_att(element_att_num,element_num)
  integer element_node(element_order,element_num)
  integer, dimension ( 3, 24 ) :: element_node_save = &
  reshape ( (/ &
    1,  2,  6, &
    7,  6,  2, &
    2,  3,  7, &
    8,  7,  3, &
    3,  4,  8, &
    9,  8,  4, &
    4,  5,  9, &
   10,  9,  5, &
    6,  7, 11, &
   12, 11,  7, &
    7,  8, 12, &
   13, 12,  8, &
    8,  9, 13, &
   14, 13,  9, &
    9, 10, 14, &
   15, 14, 10, &
   11, 12, 16, &
   17, 16, 12, &
   12, 13, 17, &
   18, 17, 13, &
   16, 17, 19, &
   20, 19, 17, &
   17, 18, 20, &
   21, 20, 18 /), (/ 3, 24 /) )

  call i4mat_copy ( element_order, element_num, element_node_save, &
    element_node )

  return
end
subroutine triangle_element_data_read ( element_filename, element_num, &
  element_order, element_att_num, element_node, element_att )

!*****************************************************************************80
!
!! TRIANGLE_ELEMENT_DATA_READ reads the data from an element file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 October 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ELEMENT_FILENAME, the name of the file.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer ELEMENT_ORDER, the order of the elements.
!
!    Input, integer ELEMENT_ATT_NUM, number of element attributes 
!    listed on each node record.
!
!    Output, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM), 
!    the indices of the nodes that make up each element.
!
!    Output, real ( kind = rk ) ELEMENT_ATT(ELEMENT_ATT_NUM,ELEMENT_NUM), the 
!    attributes of each element.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer element_att_num
  integer element_num
  integer element_order

  integer element
  real ( kind = rk ) element_att(element_att_num,element_num)
  character ( len = * ) element_filename
  integer element_node(element_order,element_num)
  integer i
  integer i1
  integer i2
  integer i3
  integer ierror
  integer input
  integer ios
  integer ival
  integer length
  character ( len = 255 ) text
  real ( kind = rk ) value

  element = 0

  call get_unit ( input )

  open ( unit = input, file = element_filename, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'TRIANGLE_ELEMENT_DATA_READ - Fatal error!'
    write ( *, '(a)' ) '  Unable to open file.'
    stop 1
  end if

  do

    read ( input, '(a)', iostat = ios ) text

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TRIANGLE_ELEMENT_DATA_READ - Fatal error!'
      write ( *, '(a)' ) '  Unexpected end of file while reading.'
      stop 1
    end if

    if ( len_trim ( text ) == 0 ) then
      cycle
    end if

    if ( text(1:1) == '#' ) then
      cycle
    end if

    if ( element == 0 ) then

      call s_to_i4 ( text, i1, ierror, length )
      text = text(length+1:)
      call s_to_i4 ( text, i2, ierror, length )
      text = text(length+1:)
      call s_to_i4 ( text, i3, ierror, length )
      text = text(length+1:)

    else

      call s_to_i4 ( text, ival, ierror, length )
      text = text(length+1:)

      do i = 1, element_order
        call s_to_i4 ( text, ival, ierror, length )
        text = text(length+1:)
        element_node(i,element) = ival
      end do

      do i = 1, element_att_num
        call s_to_r8 ( text, value, length, ierror )
        text = text(length+1:)
        element_att(i,element) = value
      end do

    end if

    element = element + 1

    if ( element_num < element ) then
      exit
    end if

  end do

  close ( unit = input )

  return
end
subroutine triangle_element_size_example ( element_num, element_order, &
  element_att_num )

!*****************************************************************************80
!
!! ELEMENT_SIZE_EXAMPLE returns the element size information for the example.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 October 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ELEMENT_NUM, the number of elements.
!
!    Output, integer ELEMENT_ORDER, the order of the elements.
!
!    Output, integer ELEMENT_ATT_NUM, the number of element 
!    attributes.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer element_att_num
  integer element_num
  integer element_order

  element_num = 24
  element_order = 3
  element_att_num = 0

  return
end
subroutine triangle_element_size_read ( element_filename, element_num, &
  element_order, element_att_num )

!*****************************************************************************80
!
!! ELEMENT_SIZE_READ reads the header information from an element file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 October 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ELEMENT_FILENAME, the name of the 
!    element file.
!
!    Output, integer ELEMENT_NUM, the number of elements.
!
!    Output, integer ELEMENT_ORDER, the order of the elements.
!
!    Output, integer ELEMENT_ATT_NUM, the number of 
!    element attributes.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer element_att_num
  character ( len = * ) element_filename
  integer element_num
  integer element_order
  integer ierror
  integer input
  integer ios
  integer length
  character ( len = 255 ) text

  call get_unit ( input )

  open ( unit = input, file = element_filename, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'TRIANGLE_ELEMENT_SIZE_READ - Fatal error!'
    write ( *, '(a)' ) '  Unable to open file.'
    stop 1
  end if

  do

    read ( input, '(a)', iostat = ios ) text

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TRIANGLE_ELEMENT_SIZE_READ - Fatal error!'
      write ( *, '(a)' ) '  Unexpected end of file while reading.'
      stop 1
    end if

    if ( len_trim ( text ) == 0 ) then
      cycle
    end if

    if ( text(1:1) == '#' ) then
      cycle
    end if

    call s_to_i4 ( text, element_num, ierror, length )
    text = text(length+1:)
    call s_to_i4 ( text, element_order, ierror, length )
    text = text(length+1:)
    call s_to_i4 ( text, element_att_num, ierror, length )
    text = text(length+1:)

    exit

  end do

  close ( unit = input )

  return
end
subroutine triangle_element_write ( element_filename, element_num, &
  element_order, element_att_num, element_node, element_att )

!*****************************************************************************80
!
!! TRIANGLE_ELEMENT_WRITE writes a TRIANGLE ".ele" file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 October 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ELEMENT_FILENAME, the name of the file
!    to be written.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer ELEMENT_ORDER, the order of the elements.
!
!    Input, integer ELEMENT_ATT_NUM, the number of element 
!    attributes.
!
!    Input, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM), 
!    the indices of the nodes that make up each element.
!
!    Input, real ( kind = rk ) ELEMENT_ATT(ELEMENT_ATT_NUM,ELEMENT_NUM), 
!    the attributes of each element.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer element_att_num
  integer element_num
  integer element_order

  integer att
  integer element
  real ( kind = rk ) element_att(element_att_num,element_num)
  character ( len = * ) element_filename
  integer element_node(element_order,element_num)
  integer ios
  integer order
  integer output

  call get_unit ( output )

  open ( unit = output, file = element_filename, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'TRIANGLE_ELEMENT_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Unable to open file.'
    stop 1
  end if

  write ( output, '(i4,2x,i4,2x,i4,2x,i4)' ) &
    element_num, element_order, element_att_num

  do element = 1, element_num

    write ( output, '(i4)', advance = 'no' ) element

    do order = 1, element_order
      write ( output, '(2x,i4)', advance = 'no' ) element_node(order,element)
    end do

    do att = 1, element_att_num
      write ( output, '(2x,g14.6)', advance = 'no' ) element_att(att,element)
    end do

    write ( output, '(a)', advance = 'yes' )

  end do

  close ( unit = output )

  return
end
