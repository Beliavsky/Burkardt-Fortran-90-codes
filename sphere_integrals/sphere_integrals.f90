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
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of the vector.
!
!    Input, integer A, B, the limits of the interval.
!
!    Output, integer X(N), a vector of numbers between A and B.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer a
  integer b
  integer i
  integer, parameter :: i4_huge = 2147483647
  integer k
  real r
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
!    07 May 2014
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
!    An R8VEC is a vector of R8's.
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!    This routine can generate a vector of values on one call.  It
!    has the feature that it should provide the same results
!    in the same order no matter how we break up the task.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of values desired.  
!
!    Output, real ( kind = rk ) X(N), a sample of the standard normal PDF.
!
!  Local:
!
!    real ( kind = rk ) R(N+1), is used to store some uniform
!    random values.  Its dimension is N+1, but really it is only needed
!    to be the smallest even number greater than or equal to N.
!
!    integer X_LO_INDEX, X_HI_INDEX, records the range of entries of
!    X that we need to compute.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer m
  real ( kind = rk ) r(n+1)
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) x(n)
  integer x_hi_index
  integer x_lo_index
!
!  Record the range of X we need to fill in.
!
  x_lo_index = 1
  x_hi_index = n
!
!  Maybe we don't need any more values.
!
  if ( x_hi_index - x_lo_index + 1 == 1 ) then

    call random_number ( harvest = r(1:2) )

    x(x_hi_index) = &
             sqrt ( -2.0D+00 * log ( r(1) ) ) * cos ( 2.0D+00 * r8_pi * r(2) )
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

  end if

  return
end
function sphere01_area ( )

!*****************************************************************************80
!
!! SPHERE01_AREA returns the area of the surface of the unit sphere in 3D.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 January 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) SPHERE01_AREA, the area of the unit sphere.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) sphere01_area
  real ( kind = rk ), parameter :: r = 1.0D+00
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00

  sphere01_area = 4.0 * r8_pi * r ** 2

  return
end
subroutine sphere01_monomial_integral ( e, integral )

!*****************************************************************************80
!
!! SPHERE01_MONOMIAL_INTEGRAL: integrals on surface of the unit sphere in 3D.
!
!  Discussion:
!
!    The integration region is 
!
!      X^2 + Y^2 + Z^2 = 1.
!
!    The monomial is F(X,Y,Z) = X^E(1) * Y^E(2) * Z^E(3).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    24 June 2002
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
!  Parameters:
!
!    Input, integer E(3), the exponents of X, Y and Z in the 
!    monomial.  Each exponent must be nonnegative.
!
!    Output, real ( kind = rk ) INTEGRAL, the integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer e(3)
  integer i
  real ( kind = rk ) integral
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00

  if ( any ( e(1:3) < 0 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPHERE01_MONOMIAL_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  All exponents must be nonnegative.'
    stop 1
  end if

  if ( any ( mod ( e(1:3), 2 ) == 1 ) ) then

    integral = 0.0D+00

  else

    integral = 2.0D+00

    do i = 1, 3
      integral = integral * gamma ( 0.5D+00 * real ( e(i) + 1, kind = rk ) )
    end do

    integral = integral &
      / gamma ( 0.5D+00 * ( real ( sum ( e(1:3) + 1 ), kind = rk ) ) )

  end if

  return
end
subroutine sphere01_sample ( n, x )

!*****************************************************************************80
!
!! SPHERE01_SAMPLE samples points from the surface of the unit sphere in 3D.
!
!  Discussion:
!
!    The sphere has center 0 and radius 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 August 2004
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
!    Input, integer N, the number of points.
!
!    Output, real ( kind = rk ) X(3,N), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer j
  real ( kind = rk ) norm
  real ( kind = rk ) x(3,n)

  do j = 1, n
!
!  Fill a vector with normally distributed values.
!
    call r8vec_normal_01 ( 3, x(1:3,j) )
!
!  Compute the length of the vector.
!
    norm = sqrt ( sum ( x(1:3,j)**2 ) )
!
!  Normalize the vector.
!
    x(1:3,j) = x(1:3,j) / norm

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
