function circle01_length ( )

!*****************************************************************************80
!
!! circle01_length(): length of the circumference of the unit circle in 2D.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 January 2014
!
!  Author:
!
!    John Burkardt
!
!  Output:
!
!    real ( kind = rk ) CIRCLE01_LENGTH, the length.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) circle01_length
  real ( kind = rk ), parameter :: r = 1.0D+00
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00

  circle01_length = 2.0 * r8_pi * r

  return
end
subroutine circle01_monomial_integral ( e, integral )

!*****************************************************************************80
!
!! circle01_monomial_integral(): integral on circumference of unit circle in 2D.
!
!  Discussion:
!
!    The integration region is 
!
!      X^2 + Y^2 = 1.
!
!    The monomial is F(X,Y) = X^E(1) * Y^E(2).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 January 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Academic Press, 1984, page 263.
!
!  Input:
!
!    integer E(2), the exponents of X and Y in the 
!    monomial.  Each exponent must be nonnegative.
!
!  Output:
!
!    real ( kind = rk ) INTEGRAL, the integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer e(2)
  integer i
  real ( kind = rk ) integral
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00

  if ( any ( e(1:2) < 0 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CIRCLE01_MONOMIAL_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  All exponents must be nonnegative.'
    stop 1
  end if

  if ( any ( mod ( e(1:2), 2 ) == 1 ) ) then

    integral = 0.0D+00

  else

    integral = 2.0D+00

    do i = 1, 2
      integral = integral * gamma ( 0.5D+00 * real ( e(i) + 1, kind = rk ) )
    end do

    integral = integral &
      / gamma ( 0.5D+00 * ( real ( sum ( e(1:2) + 1 ), kind = rk ) ) )

  end if

  return
end
subroutine circle01_sample ( n, x )

!*****************************************************************************80
!
!! circle01_sample() samples the circumference of the unit circle in 2D.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 January 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Russell Cheng,
!    Random Variate Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998, pages 168.
!
!    Reuven Rubinstein,
!    Monte Carlo Optimization, Simulation, and Sensitivity 
!    of Queueing Networks,
!    Krieger, 1992,
!    ISBN: 0894647644,
!    LC: QA298.R79.
!
!  Input:
!
!    integer N, the number of points.
!
!  Output:
!
!    real ( kind = rk ) X(2,N), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ), dimension ( 2 ) :: c = (/ 0.0D+00, 0.0D+00 /)
  real ( kind = rk ) :: r = 1.0D+00
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) theta(n)
  real ( kind = rk ) x(2,n)

  call random_number ( harvest = theta(1:n) )

  x(1,1:n) = c(1) + r * cos ( 2.0D+00 * r8_pi * theta(1:n) )
  x(2,1:n) = c(2) + r * sin ( 2.0D+00 * r8_pi * theta(1:n) )

  return
end
subroutine i4vec_uniform_ab ( n, a, b, x )

!*****************************************************************************80
!
!! i4vec_uniform_ab() returns a scaled pseudorandom I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    The pseudorandom numbers should be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the dimension of the vector.
!
!    integer A, B, the limits of the interval.
!
!  Output:
!
!    integer X(N), a vector of numbers between A and B.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer a
  integer b
  integer i
  real ( kind = rk ) r
  integer value
  integer x(n)

  do i = 1, n

    call random_number ( harvest = r )
!
!  Scale R to lie between A-0.5 and B+0.5.
!
    r = ( 1.0E+00 - r ) * ( real ( min ( a, b ) ) - 0.5E+00 ) &
      +             r   * ( real ( max ( a, b ) ) + 0.5E+00 )
!
!  Use rounding to convert R to an integer between A and B.
!
    value = nint ( r )

    value = max ( value, min ( a, b ) )
    value = min ( value, max ( a, b ) )

    x(i) = value

  end do

  return
end
subroutine monomial_value ( m, n, e, x, value )

!*****************************************************************************80
!
!! monomial_value() evaluates a monomial.
!
!  Discussion:
!
!    This routine evaluates a monomial of the form
!
!      product ( 1 <= i <= m ) x(i)^e(i)
!
!    where the exponents are nonnegative integers.  Note that
!    if the combination 0^0 is encountered, it should be treated
!    as 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 May 2014
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, the spatial dimension.
!
!    integer N, the number of points at which the
!    monomial is to be evaluated.
!
!    integer E(M), the exponents.
!
!    real ( kind = rk ) X(M,N), the point coordinates.
!
!  Output:
!
!    real ( kind = rk ) VALUE(N), the value of the monomial.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer e(m)
  integer i
  real ( kind = rk ) value(n)
  real ( kind = rk ) x(m,n)

  value(1:n) = 1.0D+00

  do i = 1, m
    if ( 0 /= e(i) ) then
      value(1:n) = value(1:n) * x(i,1:n) ** e(i)
    end if
  end do

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

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
 
