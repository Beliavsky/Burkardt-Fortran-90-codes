program main

!*****************************************************************************80
!
!! analemma() computes the analemma, equation of time, and declination.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 August 2021
!
!  Author:
!
!    Original C version by Brian Tung.
!    FORTRAN90 version by John Burkardt.
!
!  Local:
!
!    real ( kind = rk ) ECC, the orbital eccentricity.
!
!    real ( kind = rk ) LON, the longitude of the perihelion in radians.
!
!    real ( kind = rk ) OBLIQ, the obliquity in radians.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = 255 ) arg
  character ( len = 80 ) command_filename
  integer command_unit
  character ( len = 80 ) data_filename
  integer data_unit
  real ( kind = rk ) :: days = 365.242D+00
  real ( kind = rk ) dec
  real ( kind = rk ) :: degrees = ( 3.141592653589793D+00 / 180.0D+00 )
  real ( kind = rk ) ecc
  real ( kind = rk ) eot
  real ( kind = rk ) f
  integer i
  real ( kind = rk ) lon
  integer n
  integer numarg
  real ( kind = rk ) obliq
  real ( kind = rk ) :: pi = 3.141592653589793D+00
  real ( kind = rk ) t
  real ( kind = rk ) tau
  real ( kind = rk ) theta
  real ( kind = rk ) x1
  real ( kind = rk ) x2
  real ( kind = rk ) x3
  real ( kind = rk ) y1
  real ( kind = rk ) y2
  real ( kind = rk ) y3
  real ( kind = rk ) z1
  real ( kind = rk ) z2
  real ( kind = rk ) z3

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'analemma():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) &
    '  Compute and plot the analemma, equation of time, and declination.'
  write ( *, '(a)' ) '  This program is based on a C program by Brian Tung.'
!
!  Parse the arguments 
!
  numarg = iargc ( )

  if ( numarg < 1 ) then
    ecc = 0.01671D+00
  else
    i = 1
    call getarg ( i, arg )
    call s_to_r8 ( arg, ecc )
  end if

  if ( numarg < 2 ) then
    lon = 1.347D+00
  else
    i = 2
    call getarg ( i, arg )
    call s_to_r8 ( arg, lon )
    lon = lon * pi / 180.0D+00
  end if

  if ( numarg < 3 ) then
    obliq = 0.4091D+00
  else
    i = 3
    call getarg ( i, arg )
    call s_to_r8 ( arg, obliq )
    obliq = obliq * pi / 180.0D+00
  end if
!
!  Compute the data.
!
  data_filename = 'analemma_data.txt'
  call get_unit ( data_unit )
  open ( unit = data_unit, file = data_filename, status = 'replace' )

  n = 10000

  do i = 0, n - 1

    f = real ( i, kind = rk ) / real ( n, kind = rk )
    tau = 2.0D+00 * pi * f
!
!  Set theta to the current longitude. 
!
    theta = atan2 ( sqrt ( 1.0D+00 - ecc * ecc ) * sin ( tau ), &
      cos ( tau ) - ecc )
! 
!  Rotate clockwise in XY plane by theta, corrected by lon.
!
    x1 = cos ( theta - ( lon - pi / 2.0D+00 ) )
    y1 = sin ( theta - ( lon - pi / 2.0D+00 ) )
    z1 = 0.0D+00
! 
!  Rotate counter-clockwise in XZ plane by obliq.
!
    x2 =   cos ( obliq ) * x1 + sin ( obliq ) * z1
    y2 = y1
    z2 = - sin ( obliq ) * x1 + cos ( obliq ) * z1
! 
!  Set t equal to real time from tau and
!  rotate counter-clockwise by t, corrected by lon 
!
    t = tau - ecc * sin ( tau );
    x3 =   cos ( t - ( lon - pi / 2.0D+00 ) ) * x2 &
      + sin ( t - ( lon - pi / 2.0D+00 ) ) * y2;
    y3 = - sin ( t - ( lon - pi / 2.0D+00 ) ) * x2 &
      + cos ( t - ( lon - pi / 2.0D+00 ) ) * y2;
    z3 = z2;

    eot = - atan2 ( y3, x3 ) * 4.0D+00 / degrees * days / ( days + 1.0D+00 )
    dec = asin ( z3 ) / degrees;
! 
!  Print results in minutes early/late and degrees north/south 
!
    write ( data_unit, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) t &
      / ( 2.0D+00 * pi ), eot, dec

  end do

  close ( unit = data_unit )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Created data file "' // trim ( data_filename ) // '".'
!
!  Create the command file.
!
  command_filename = 'analemma_commands.txt'
  call get_unit ( command_unit )
  open ( unit = command_unit, file = command_filename, status = 'replace' )
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set style data lines'
  write ( command_unit, '(a)' ) 'unset key'
  write ( command_unit, '(a)' ) 'set output "eot.png"'
  write ( command_unit, '(a)' ) 'set xlabel "<---Normalized Date--->"'
  write ( command_unit, '(a)' ) 'set ylabel "<---Minutes Early/Late--->"'
  write ( command_unit, '(a)' ) 'set title "The equation of time"'
  write ( command_unit, '(a)' ) 'plot "' // trim ( data_filename ) &
    // '" using 1:2 with lines'
  write ( command_unit, '(a)' ) 'set output "declination.png"'
  write ( command_unit, '(a)' ) 'set xlabel "<---Normalized Date--->"'
  write ( command_unit, '(a)' ) 'set ylabel "<---Degrees North/South--->"'
  write ( command_unit, '(a)' ) 'set title "Declination"'
  write ( command_unit, '(a)' ) 'plot "' // trim ( data_filename ) &
    // '" using 1:3 with lines'
  write ( command_unit, '(a)' ) 'set output "analemma.png"'
  write ( command_unit, '(a)' ) 'set xlabel "<---Minutes Early/Late--->"'
  write ( command_unit, '(a)' ) 'set ylabel "<---Degrees North/South--->"'
  write ( command_unit, '(a)' ) 'set title "The analemma"'
  write ( command_unit, '(a)' ) 'plot "' // trim ( data_filename ) &
    // '" using 2:3 with lines'
  write ( command_unit, '(a)' ) 'quit'

  close ( unit = command_unit )

  write ( *, '(a)' ) '  Created command file "' &
    // trim ( command_filename ) // '".'
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'analemma:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
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
!    26 August 2021
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
subroutine s_to_r8 ( s, r8 )

!*****************************************************************************80
!
!! s_to_r8() reads a real value from a string.
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
!    S                 R8
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
!    26 August 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    character ( len = * ) S, the string containing the
!    data to be read.  Reading will begin at position 1 and
!    terminate at the end of the string, or when no more
!    characters can be read to form a legal real.  Blanks,
!    commas, or other nonnumeric data will, in particular,
!    cause the conversion to halt.
!
!  Output:
!
!    real ( kind = rk ) R8, the value read from the string.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character c
  integer ierror
  integer ihave
  integer isgn
  integer iterm
  integer jbot
  integer jsgn
  integer jtop
  integer length
  integer ndig
  real ( kind = rk ) r8
  real ( kind = rk ) rbot
  real ( kind = rk ) rexp
  real ( kind = rk ) rtop
  character ( len = * ) s
  integer s_length
  character  :: TAB = achar ( 9 )

  s_length = len_trim ( s )

  ierror = 0
  r8 = 0.0D+00
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

    if ( s_length < length + 1 ) then
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
    else if ( c == 'E' .or. c == 'e' .or. c == 'D' .or. c == 'd' ) then

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

      ndig = iachar ( c ) - 48

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
  if ( iterm /= 1 .and. length + 1 == s_length ) then
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
    stop
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

  r8 = real ( isgn, kind = rk ) * rexp * rtop / rbot

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
!    26 August 2021
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
 
