subroutine r4_machar ( ibeta, it, irnd, ngrd, machep, negep, iexp, &
  minexp, maxexp, eps, epsneg, xmin, xmax )

!*****************************************************************************80
!
!! R4_MACHAR determines single precision machine constants.
!
!  Discussion:
!
!    This routine determines the parameters of the floating-point 
!    arithmetic system specified below.  The determination of the first 
!    three uses an extension of an algorithm due to Malcolm, 
!    incorporating some of the improvements suggested by Gentleman and 
!    Marovich.  
!
!    This routine appeared as ACM algorithm 665.
!
!    An earlier version of this program was published in Cody and Waite.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 June 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody,
!    ACM Algorithm 665, MACHAR, a subroutine to dynamically determine 
!    machine parameters,
!    ACM Transactions on Mathematical Software,
!    Volume 14, Number 4, pages 303-311, 1988.
!
!    William Cody, William Waite,
!    Software Manual for the Elementary Functions,
!    Prentice Hall, 1980.
!
!    Morven Gentleman, Scott Marovich,
!    Communications of the ACM,
!    Volume 17, pages 276-277, 1974.
!
!    Michael Malcolm,
!    Communications of the ACM,
!    Volume 15, pages 949-951, 1972.
!
!  Parameters:
!
!    Output, integer IBETA, the radix for the floating-point 
!    representation.
!
!    Output, integer IT, the number of base IBETA digits 
!    in the floating-point significand.
!
!    Output, integer IRND:
!    0, if floating-point addition chops.
!    1, if floating-point addition rounds, but not in the IEEE style.
!    2, if floating-point addition rounds in the IEEE style.
!    3, if floating-point addition chops, and there is partial underflow.
!    4, if floating-point addition rounds, but not in the IEEE style, and 
!      there is partial underflow.
!    5, if floating-point addition rounds in the IEEE style, and there is 
!      partial underflow.
!
!    Output, integer NGRD, the number of guard digits for 
!    multiplication with truncating arithmetic.  It is
!    0, if floating-point arithmetic rounds, or if it truncates and only 
!      IT base IBETA digits participate in the post-normalization shift of the
!      floating-point significand in multiplication;
!    1, if floating-point arithmetic truncates and more than IT base IBETA
!      digits participate in the post-normalization shift of the floating-point
!      significand in multiplication.
!
!    Output, integer MACHEP, the largest negative integer such that
!      1.0 + real ( IBETA ) ^ MACHEP /= 1.0, 
!    except that MACHEP is bounded below by - ( IT + 3 ).
!
!    Output, integer NEGEPS, the largest negative integer such that
!      1.0 - real ( IBETA ) ^ NEGEPS /= 1.0, 
!    except that NEGEPS is bounded below by - ( IT + 3 ).
!
!    Output, integer IEXP, the number of bits (decimal places 
!    if IBETA = 10) reserved for the representation of the exponent (including 
!    the bias or sign) of a floating-point number.
!
!    Output, integer MINEXP, the largest in magnitude negative 
!    integer such that
!      real ( IBETA ) ^ MINEXP 
!    is positive and normalized.
!
!    Output, integer MAXEXP, the smallest positive power of 
!    BETA that overflows.
!
!    Output, real ( kind = rk ) EPS, the smallest positive floating-point 
!    number such that  
!      1.0 + EPS /= 1.0. 
!    in particular, if either IBETA = 2  or IRND = 0, 
!      EPS = real ( IBETA ) ^ MACHEP.
!    Otherwise,  
!      EPS = ( real ( IBETA ) ^ MACHEP ) / 2.
!
!    Output, real ( kind = rk ) EPSNEG, a small positive floating-point number 
!    such that
!      1.0 - EPSNEG /= 1.0. 
!    In particular, if IBETA = 2 or IRND = 0, 
!      EPSNEG = real ( IBETA ) ^ NEGEPS.
!    Otherwise,  
!      EPSNEG = ( real ( IBETA ) ^ NEGEPS ) / 2.  
!    Because NEGEPS is bounded below by - ( IT + 3 ), EPSNEG might not be the
!    smallest number that can alter 1.0 by subtraction.
!
!    Output, real ( kind = rk ) XMIN, the smallest non-vanishing normalized 
!    floating-point power of the radix:
!      XMIN = real ( IBETA ) ^ MINEXP
!
!    Output, real ( kind = rk ) XMAX, the largest finite floating-point number.
!    In particular,
!      XMAX = ( 1.0 - EPSNEG ) * real ( IBETA ) ^ MAXEXP
!    On some machines, the computed value of XMAX will be only the second, 
!    or perhaps third, largest number, being too small by 1 or 2 units in 
!    the last digit of the significand.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0E+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) beta
  real ( kind = rk ) betah
  real ( kind = rk ) betain
  real ( kind = rk ) eps
  real ( kind = rk ) epsneg
  integer i
  integer ibeta
  integer iexp
  integer irnd
  integer it
  integer itemp
  integer iz
  integer j
  integer k
  integer machep
  integer maxexp
  integer minexp
  integer mx
  integer negep
  integer ngrd
  integer nxres
  real ( kind = rk ) one
  real ( kind = rk ) t
  real ( kind = rk ) temp
  real ( kind = rk ) temp1
  real ( kind = rk ) tempa
  real ( kind = rk ) two
  real ( kind = rk ) xmax
  real ( kind = rk ) xmin
  real ( kind = rk ) y
  real ( kind = rk ) z
  real ( kind = rk ) zero

  one = real ( 1, kind = rk )
  two = one + one
  zero = one - one
!
!  Determine IBETA and BETA ala Malcolm.
!
  a = one

  do
  
    a = a + a
    temp = a + one
    temp1 = temp - a

    if ( temp1 - one /= zero ) then
      exit
    end if

  end do

  b = one

  do

    b = b + b
    temp = a + b
    itemp = int ( temp - a )

    if ( itemp /= 0 ) then
      exit
    end if

  end do

  ibeta = itemp
  beta = real ( ibeta, kind = rk )
!
!  Determine IT and IRND.
!
  it = 0
  b = one

  do

    it = it + 1
    b = b * beta
    temp = b + one
    temp1 = temp - b

    if ( temp1 - one /= zero ) then
      exit
    end if

  end do

  irnd = 0
  betah = beta / two
  temp = a + betah

  if ( temp - a /= zero ) then
    irnd = 1
  end if

  tempa = a + beta
  temp = tempa + betah

  if ( irnd == 0 .and. temp - tempa /= zero ) then
    irnd = 2
  end if
!
!  Determine NEGEP and EPSNEG.
!
  negep = it + 3
  betain = one / beta
  a = one
  do i = 1, negep
    a = a * betain
  end do

  b = a

  do

    temp = one - a

    if ( temp - one /= zero ) then
      exit
    end if

    a = a * beta
    negep = negep - 1

  end do

  negep = -negep
  epsneg = a

  if ( ibeta /= 2 .and. irnd /= 0 ) then

    a = ( a * ( one + a ) ) / two
    temp = one - a

    if ( temp - one /= zero ) then
      epsneg = a
    end if

  end if
!
!  Determine MACHEP and EPS.
!
  machep = -it - 3
  a = b

  do

    temp = one + a

    if ( temp - one /= zero ) then
      exit
    end if

    a = a * beta
    machep = machep + 1

  end do

  eps = a
  temp = tempa + beta * ( one + eps )

  if ( ibeta /= 2 .and. irnd /= 0 ) then

    a = ( a * ( one + a ) ) / two
    temp = one + a

    if ( temp - one /= zero ) then
      eps = a
    end if

  end if
!
!  Determine NGRD.
!
  ngrd = 0
  temp = one + eps

  if ( irnd == 0 .and. temp * one - one /= zero ) then
    ngrd = 1
  end if
!
!  Determine IEXP, MINEXP and XMIN.
!
!  Loop to determine largest I and K = 2^I such that (1/BETA) ^ (2^(I))
!  does not underflow.  Exit from loop is signaled by an underflow.
!
  i = 0
  k = 1
  z = betain
  t = one + eps
  nxres = 0

  do

    y = z
    z = y * y

    a = z * one
    temp = z * t

    if ( a + a == zero .or. y <= abs ( z ) ) then
      exit
    end if

    temp1 = temp * betain

    if ( temp1 * beta == z ) then
      exit
    end if

    i = i + 1
    k = k + k

  end do
!
!  This segment is for nondecimal machines.
!
  if ( ibeta /= 10 ) then

    iexp = i + 1
    mx = k + k
!
!  This segment is for decimal machines only.
!
  else

    iexp = 2
    iz = ibeta

    do

      if ( k < iz ) then
        exit
      end if

      iz = iz * ibeta
      iexp = iexp + 1

    end do

    mx = iz + iz - 1

  end if
!
!  Loop to determine MINEXP, XMIN.
!  Exit from loop is signaled by an underflow.
!
  do

    xmin = y
    y = y * betain

    a = y * one
    temp = y * t

    if ( a + a == zero .or. xmin <= abs ( y ) ) then
      exit
    end if

    k = k + 1
    temp1 = temp * betain

    if ( temp1 * beta == y ) then
      nxres = 3
      xmin = y
      exit
    end if

  end do

  minexp = -k
!
!  Determine MAXEXP and XMAX.
!
  if ( mx <= k + k - 3 .and. ibeta /= 10 ) then
    mx = mx + mx
    iexp = iexp + 1
  end if

  maxexp = mx + minexp
!
!  Adjust IRND to reflect partial underflow.
!
  irnd = irnd + nxres
!
!  Adjust for IEEE-style machines.
!
  if ( irnd == 2 .or. irnd == 5 ) then
    maxexp = maxexp - 2
  end if
!
!  Adjust for non-IEEE machines with partial underflow.
!
  if ( irnd == 3 .or. irnd == 4 ) then
    maxexp = maxexp - it
  end if
!
!  Adjust for machines with implicit leading bit in binary significand, 
!  and machines with radix point at extreme right of significand.
!
  i = maxexp + minexp

  if ( ibeta == 2 .and. i == 0 ) then
    maxexp = maxexp - 1
  end if

  if ( 20 < i ) then
    maxexp = maxexp - 1
  end if

  if ( a /= y ) then
    maxexp = maxexp - 2
  end if

  xmax = one - epsneg

  if ( xmax * one /= xmax ) then
    xmax = one - beta * epsneg
  end if

  xmax = xmax / ( beta * beta * beta * xmin )

  i = maxexp + minexp + 3

  do j = 1, i

    if ( ibeta == 2 ) then
      xmax = xmax + xmax
    else
      xmax = xmax * beta
    end if

  end do

  return
end
subroutine r8_machar ( ibeta, it, irnd, ngrd, machep, negep, iexp, &
  minexp, maxexp, eps, epsneg, xmin, xmax )

!*****************************************************************************80
!
!! R8_MACHAR determines double precision machine constants.
!
!  Discussion:
!
!    This routine determines the parameters of the floating-point 
!    arithmetic system specified below.  The determination of the first 
!    three uses an extension of an algorithm due to Malcolm, 
!    incorporating some of the improvements suggested by Gentleman and 
!    Marovich.  
!
!    This routine appeared as ACM algorithm 665.
!
!    An earlier version of this program was published in Cody and Waite.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 June 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody,
!    ACM Algorithm 665, MACHAR, a subroutine to dynamically determine 
!    machine parameters,
!    ACM Transactions on Mathematical Software,
!    Volume 14, Number 4, pages 303-311, 1988.
!
!    William Cody, William Waite,
!    Software Manual for the Elementary Functions,
!    Prentice Hall, 1980.
!
!    Morven Gentleman, Scott Marovich,
!    Communications of the ACM,
!    Volume 17, pages 276-277, 1974.
!
!    Michael Malcolm,
!    Communications of the ACM,
!    Volume 15, pages 949-951, 1972.
!
!  Parameters:
!
!    Output, integer IBETA, the radix for the floating-point
!    representation.
!
!    Output, integer IT, the number of base IBETA digits in 
!    the floating-point significand.
!
!    Output, integer IRND:
!    0, if floating-point addition chops.
!    1, if floating-point addition rounds, but not in the IEEE style.
!    2, if floating-point addition rounds in the IEEE style.
!    3, if floating-point addition chops, and there is partial underflow.
!    4, if floating-point addition rounds, but not in the IEEE style, and 
!      there is partial underflow.
!    5, if floating-point addition rounds in the IEEE style, and there is 
!      partial underflow.
!
!    Output, integer NGRD, the number of guard digits for 
!    multiplication with truncating arithmetic.  It is
!    0, if floating-point arithmetic rounds, or if it truncates and only 
!      IT base IBETA digits participate in the post-normalization shift of the
!      floating-point significand in multiplication;
!    1, if floating-point arithmetic truncates and more than IT base IBETA
!      digits participate in the post-normalization shift of the floating-point
!      significand in multiplication.
!
!    Output, integer MACHEP, the largest negative integer 
!    such that
!      1.0 < 1.0 + real ( IBETA, kind = rk ) ^ MACHEP, 
!    except that MACHEP is bounded below by - ( IT + 3 ).
!
!    Output, integer NEGEPS, the largest negative integer 
!    such that
!      1.0 - real ( IBETA, kind = rk ) ^ NEGEPS < 1.0, 
!    except that NEGEPS is bounded below by - ( IT + 3 ).
!
!    Output, integer IEXP, the number of bits (decimal places 
!    if IBETA = 10) reserved for the representation of the exponent (including
!    the bias or sign) of a floating-point number.
!
!    Output, integer MINEXP, the largest in magnitude negative
!    integer such that
!      real ( IBETA, kind = rk ) ^ MINEXP 
!    is positive and normalized.
!
!    Output, integer MAXEXP, the smallest positive power of
!    BETA that overflows.
!
!    Output, real ( kind = rk ) EPS, the smallest positive floating-point number
!    such that  
!      1.0 + EPS /= 1.0. 
!    in particular, if either IBETA = 2  or IRND = 0, 
!      EPS = real ( IBETA, kind = rk ) ^ MACHEP.
!    Otherwise,  
!      EPS = ( real ( IBETA, kind = rk ) ^ MACHEP ) / 2.
!
!    Output, real ( kind = rk ) EPSNEG, a small positive floating-point number
!    such that
!      1.0 - EPSNEG < 1.0. 
!    In particular, if IBETA = 2 or IRND = 0, 
!      EPSNEG = real ( IBETA, kind = rk ) ^ NEGEPS.
!    Otherwise,  
!      EPSNEG = ( real ( IBETA, kind = rk ) ^ NEGEPS ) / 2.  
!    Because NEGEPS is bounded below by - ( IT + 3 ), EPSNEG might not be the
!    smallest number that can alter 1.0 by subtraction.
!
!    Output, real ( kind = rk ) XMIN, the smallest non-vanishing normalized
!    floating-point power of the radix:
!      XMIN = real ( IBETA, kind = rk ) ^ MINEXP
!
!    Output, real ( kind = rk ) XMAX, the largest finite floating-point number.
!    In particular,
!      XMAX = ( 1.0 - EPSNEG ) * real ( IBETA, kind = rk ) ^ MAXEXP
!    On some machines, the computed value of XMAX will be only the second, 
!    or perhaps third, largest number, being too small by 1 or 2 units in 
!    the last digit of the significand.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) beta
  real ( kind = rk ) betah
  real ( kind = rk ) betain
  real ( kind = rk ) eps
  real ( kind = rk ) epsneg
  integer i
  integer ibeta
  integer iexp
  integer irnd
  integer it
  integer itemp
  integer iz
  integer j
  integer k
  integer machep
  integer maxexp
  integer minexp
  integer mx
  integer negep
  integer ngrd
  integer nxres
  real ( kind = rk ) one
  real ( kind = rk ) t
  real ( kind = rk ) temp
  real ( kind = rk ) temp1
  real ( kind = rk ) tempa
  real ( kind = rk ) two
  real ( kind = rk ) xmax
  real ( kind = rk ) xmin
  real ( kind = rk ) y
  real ( kind = rk ) z
  real ( kind = rk ) zero

  one = real ( 1, kind = rk )
  two = one + one
  zero = one - one
!
!  Determine IBETA and BETA ala Malcolm.
!
  a = one

  do

    a = a + a
    temp = a + one
    temp1 = temp - a

    if ( temp1 - one /= zero ) then
      exit
    end if

  end do

  b = one

  do

    b = b + b
    temp = a + b
    itemp = int ( temp - a )

    if ( itemp /= 0 ) then
      exit
    end if

  end do

  ibeta = itemp
  beta = real ( ibeta, kind = rk )
!
!  Determine IT and IRND.
!
  it = 0
  b = one

  do

    it = it + 1
    b = b * beta
    temp = b + one
    temp1 = temp - b

    if ( temp1 - one /= zero ) then
      exit
    end if

  end do

  irnd = 0
  betah = beta / two
  temp = a + betah

  if ( temp - a /= zero ) then
    irnd = 1
  end if

  tempa = a + beta
  temp = tempa + betah

  if ( irnd == 0 .and. temp - tempa /= zero ) then
    irnd = 2
  end if
!
!  Determine NEGEP and EPSNEG.
!
  negep = it + 3
  betain = one / beta
  a = one
  do i = 1, negep
    a = a * betain
  end do

  b = a

  do

    temp = one - a

    if ( temp - one /= zero ) then
      exit
    end if

    a = a * beta
    negep = negep - 1

  end do

  negep = -negep
  epsneg = a

  if ( ibeta /= 2 .and. irnd /= 0 ) then

    a = ( a * ( one + a ) ) / two
    temp = one - a

    if ( temp - one /= zero ) then
      epsneg = a
    end if

  end if
!
!  Determine MACHEP and EPS.
!
  machep = -it - 3
  a = b

  do

    temp = one + a

    if ( temp - one /= zero ) then
      exit
    end if

    a = a * beta
    machep = machep + 1

  end do

  eps = a
  temp = tempa + beta * ( one + eps )

  if ( ibeta /= 2 .and. irnd /= 0 ) then

    a = ( a * ( one + a ) ) / two
    temp = one + a

    if ( temp - one /= zero ) then
      eps = a
    end if

  end if
!
!  Determine NGRD.
!
  ngrd = 0
  temp = one + eps

  if ( irnd == 0 .and. temp * one - one /= zero ) then
    ngrd = 1
  end if
!
!  Determine IEXP, MINEXP and XMIN.
!
!  Loop to determine largest I and K = 2^I such that (1/BETA) ^ (2^I)
!  does not underflow.  Exit from loop is signaled by an underflow.
!
  i = 0
  k = 1
  z = betain
  t = one + eps
  nxres = 0

  do

    y = z
    z = y * y

    a = z * one
    temp = z * t

    if ( a + a == zero .or. y <= abs ( z ) ) then
      exit
    end if

    temp1 = temp * betain

    if ( temp1 * beta == z ) then
      exit
    end if

    i = i + 1
    k = k + k

  end do
!
!  This segment is for nondecimal machines.
!
  if ( ibeta /= 10 ) then

    iexp = i + 1
    mx = k + k
!
!  This segment is for decimal machines only.
!
  else

    iexp = 2
    iz = ibeta

    do

      if ( k < iz ) then
        exit
      end if

      iz = iz * ibeta
      iexp = iexp + 1

    end do

    mx = iz + iz - 1

  end if
!
!  Loop to determine MINEXP, XMIN.
!  Exit from loop is signaled by an underflow.
!
  do

    xmin = y
    y = y * betain

    a = y * one
    temp = y * t

    if ( a + a == zero .or. xmin <= abs ( y ) ) then
      exit
    end if

    k = k + 1
    temp1 = temp * betain

    if ( temp1 * beta == y ) then
      nxres = 3
      xmin = y
      exit
    end if

  end do

  minexp = -k
!
!  Determine MAXEXP and XMAX.
!
  if ( mx <= k + k - 3 .and. ibeta /= 10 ) then
    mx = mx + mx
    iexp = iexp + 1
  end if

  maxexp = mx + minexp
!
!  Adjust IRND to reflect partial underflow.
!
  irnd = irnd + nxres
!
!  Adjust for IEEE-style machines.
!
  if ( irnd == 2 .or. irnd == 5 ) then
    maxexp = maxexp - 2
  end if
!
!  Adjust for non-IEEE machines with partial underflow.
!
  if ( irnd == 3 .or. irnd == 4 ) then
    maxexp = maxexp - it
  end if
!
!  Adjust for machines with implicit leading bit in binary significand, 
!  and machines with radix point at extreme right of significand.
!
  i = maxexp + minexp

  if ( ibeta == 2 .and. i == 0 ) then
    maxexp = maxexp - 1
  end if

  if ( 20 < i ) then
    maxexp = maxexp - 1
  end if

  if ( a /= y ) then
    maxexp = maxexp - 2
  end if

  xmax = one - epsneg

  if ( xmax * one /= xmax ) then
    xmax = one - beta * epsneg
  end if

  xmax = xmax / ( beta * beta * beta * xmin )

  i = maxexp + minexp + 3

  do j = 1, i

    if ( ibeta == 2 ) then
      xmax = xmax + xmax
    else
      xmax = xmax * beta
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
!    06 August 2005
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
