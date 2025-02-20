subroutine hyperball01_monomial_integral ( m, e, integral )

!*****************************************************************************80
!
!! hyperball01_monomial_integral(): monomial integrals in unit hyperball.
!
!  Discussion:
!
!    The integration region is 
!
!      sum ( 1 <= I <= M ) X(I)^2 <= 1.
!
!    The monomial is F(X) = product ( 1 <= I <= M ) X(I)^E(I)
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 January 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Gerald Folland,
!    How to Integrate a Polynomial Over a Sphere,
!    American Mathematical Monthly,
!    Volume 108, Number 5, May 2001, pages 446-448.
!
!  Parameters:
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer E(M), the exponents.  
!    Each exponent must be nonnegative.
!
!    Output, real ( kind = rk ) INTEGRAL, the integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m

  real ( kind = rk ) arg
  integer e(m)
  integer i
  real ( kind = rk ) integral
  real ( kind = rk ), parameter :: r = 1.0D+00
  real ( kind = rk ) s

  if ( any ( e(1:m) < 0 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HYPERBALL01_MONOMIAL_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  All exponents must be nonnegative.'
    stop 1
  end if
!
!  Integrate over the surface.
!
  if ( any ( mod ( e(1:m), 2 ) == 1 ) ) then

    integral = 0.0D+00

  else

    integral = 2.0D+00

    do i = 1, m
      arg = 0.5D+00 * real ( e(i) + 1, kind = rk )
      integral = integral * gamma ( arg )
    end do

    arg = 0.5D+00 * real ( sum ( e(1:m) + 1 ), kind = rk )
    integral = integral / gamma ( arg )

  end if
!
!  The surface integral is now adjusted to give the volume integral.
!
  s = sum ( e(1:m) ) + m

  integral = integral * r ** s / real ( s, kind = rk )

  return
end
subroutine hyperball01_sample ( m, n, x )

!*****************************************************************************80
!
!! HYPERBALL01_SAMPLE uniformly samples points from the unit hyperball.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 January 2014
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
!  Parameters:
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of points.
!
!    Output, real ( kind = rk ) X(M,N), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) exponent
  integer j
  real ( kind = rk ) norm
  real ( kind = rk ) r
  real ( kind = rk ) x(m,n)

  exponent = 1.0D+00 / real ( m, kind = rk )

  do j = 1, n
!
!  Fill a vector with normally distributed values.
!
    call r8vec_normal_01 ( m, x(1:m,j) )
!
!  Compute the length of the vector.
!
    norm = sqrt ( sum ( x(1:m,j)**2 ) )
!
!  Normalize the vector.
!
    x(1:m,j) = x(1:m,j) / norm
!
!  Now compute a value to map the point ON the sphere INTO the sphere.
!
    call random_number ( harvest = r )

    x(1:m,j) = r ** exponent * x(1:m,j)

  end do

  return
end
function hyperball01_volume ( m )

!*****************************************************************************80
!
!! HYPERBALL01_VOLUME returns the volume of the unit hyperball.
!
!  Discussion:
!
!     M  Volume
!
!     1    2
!     2    1        * PI
!     3  ( 4 /   3) * PI
!     4  ( 1 /   2) * PI^2
!     5  ( 8 /  15) * PI^2
!     6  ( 1 /   6) * PI^3
!     7  (16 / 105) * PI^3
!     8  ( 1 /  24) * PI^4
!     9  (32 / 945) * PI^4
!    10  ( 1 / 120) * PI^5
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 January 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the spatial dimension.
!
!    Output, real ( kind = rk ) HYPERBALL01_VOLUME, the volume.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) hyperball01_volume
  integer i
  integer m
  integer m_half
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) volume

  if ( mod ( m, 2 ) == 0 ) then
    m_half = m / 2
    volume = r8_pi ** m_half
    do i = 1, m_half
      volume = volume / real ( i, kind = rk )
    end do
  else
    m_half = ( m - 1 ) / 2
    volume = r8_pi ** m_half * 2.0D+00 ** m
    do i = m_half + 1, 2 * m_half + 1
      volume = volume / real ( i, kind = rk )
    end do
  end if

  hyperball01_volume = volume

  return
end
subroutine monomial_value ( m, n, e, x, value )

!*****************************************************************************80
!
!! MONOMIAL_VALUE evaluates a monomial.
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
!    20 April 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of points at which the
!    monomial is to be evaluated.
!
!    Input, integer E(M), the exponents.
!
!    Input, real ( kind = rk ) X(M,N), the point coordinates.
!
!    Output, real ( kind = rk ) VALUE(N), the value of the monomial.
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
subroutine r8vec_normal_01 ( n, x )

!*****************************************************************************80
!
!! R8VEC_NORMAL_01 returns a unit pseudonormal R8VEC.
!
!  Discussion:
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!    The Box-Muller method is used, which is efficient, but
!    generates an even number of values each time.  On any call
!    to this routine, an even number of new values are generated.
!    Depending on the situation, one value may be left over.
!    In that case, it is saved for the next call.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of values desired.  If N is 
!    negative, then the code will flush its internal memory; in particular,
!    if there is a saved value to be used on the next call, it is
!    instead discarded.  
!
!    Output, real ( kind = rk ) X(N), a sample of the standard normal PDF.
!
!  Local:
!
!    Local, integer MADE, records the number of values that have
!    been computed.  On input with negative N, this value overwrites
!    the return value of N, so the user can get an accounting of
!    how much work has been done.
!
!    Local, real ( kind = rk ) R(N+1), is used to store some uniform
!    random values.  Its dimension is N+1, but really it is only needed
!    to be the smallest even number greater than or equal to N.
!
!    Local, integer SAVED, is 0 or 1 depending on whether there is a
!    single saved value left over from the previous call.
!
!    Local, integer X_LO_INDEX, X_HI_INDEX, records the range of entries of
!    X that we need to compute.  This starts off as 1:N, but is adjusted
!    if we have a saved value that can be immediately stored in X(1),
!    and so on.
!
!    Local, real ( kind = rk ) Y, the value saved from the previous call, if
!    SAVED is 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer m
  integer, save :: made = 0
  real ( kind = rk ) r(n+1)
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  integer, save :: saved = 0
  real ( kind = rk ) x(n)
  integer x_hi_index
  integer x_lo_index
  real ( kind = rk ), save :: y = 0.0D+00
!
!  I'd like to allow the user to reset the internal data.
!  But this won't work properly if we have a saved value Y.
!  I'm making a crock option that allows the user to signal
!  explicitly that any internal memory should be flushed,
!  by passing in a negative value for N.
!
  if ( n < 0 ) then
    n = made
    made = 0
    saved = 0
    y = 0.0D+00
    return
  else if ( n == 0 ) then
    return
  end if
!
!  Record the range of X we need to fill in.
!
  x_lo_index = 1
  x_hi_index = n
!
!  Use up the old value, if we have it.
!
  if ( saved == 1 ) then
    x(1) = y
    saved = 0
    x_lo_index = 2
  end if
!
!  Maybe we don't need any more values.
!
  if ( x_hi_index - x_lo_index + 1 == 0 ) then
!
!  If we need just one new value, do that here to avoid null arrays.
!
  else if ( x_hi_index - x_lo_index + 1 == 1 ) then

    call random_number ( harvest = r(1:2) )

    x(x_hi_index) = &
             sqrt ( - 2.0D+00 * log ( r(1) ) ) * cos ( 2.0D+00 * r8_pi * r(2) )
    y =      sqrt ( - 2.0D+00 * log ( r(1) ) ) * sin ( 2.0D+00 * r8_pi * r(2) )

    saved = 1

    made = made + 2
!
!  If we require an even number of values, that's easy.
!
  else if ( mod ( x_hi_index - x_lo_index + 1, 2 ) == 0 ) then

    m = ( x_hi_index - x_lo_index + 1 ) / 2

    call random_number ( harvest = r(1:2*m) )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * cos ( 2.0D+00 * r8_pi * r(2:2*m:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * sin ( 2.0D+00 * r8_pi * r(2:2*m:2) )

    made = made + x_hi_index - x_lo_index + 1
!
!  If we require an odd number of values, we generate an even number,
!  and handle the last pair specially, storing one in X(N), and
!  saving the other for later.
!
  else

    x_hi_index = x_hi_index - 1

    m = ( x_hi_index - x_lo_index + 1 ) / 2 + 1

    call random_number ( harvest = r(1:2*m) )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * cos ( 2.0D+00 * r8_pi * r(2:2*m-2:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * sin ( 2.0D+00 * r8_pi * r(2:2*m-2:2) )

    x(n) = sqrt ( -2.0D+00 * log ( r(2*m-1) ) ) &
      * cos ( 2.0D+00 * r8_pi * r(2*m) )

    y = sqrt ( -2.0D+00 * log ( r(2*m-1) ) ) &
      * sin ( 2.0D+00 * r8_pi * r(2*m) )

    saved = 1

    made = made + x_hi_index - x_lo_index + 2

  end if

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

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
