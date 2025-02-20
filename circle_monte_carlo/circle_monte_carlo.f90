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
!    05 September 2021
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

  circle01_length = 2.0D+00 * r8_pi * r

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
!    05 September 2021
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
    write ( *, '(a)' ) 'circle01_monomial_integral(): Fatal error!'
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
subroutine circle01_sample_ergodic ( n, angle, x )

!*****************************************************************************80
!
!! circle01_sample_ergodic() samples the circumference of the unit circle in 2D.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of points.
!
!    real ( kind = rk ) ANGLE, an angle between 0 and 2*PI.
!
!  Output:
!
!    real ( kind = rk ) ANGLE, a new angle between 0 and 2*PI.
!
!    real ( kind = rk ) X(2,N), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) angle
  real ( kind = rk ), dimension ( 2 ) :: c = (/ 0.0D+00, 0.0D+00 /)
  real ( kind = rk ) golden_angle
  real ( kind = rk ) golden_ratio
  integer j
  real ( kind = rk ) :: r = 1.0D+00
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) x(2,n)

  golden_ratio = ( 1.0D+00 + sqrt ( 5.0D+00 ) ) / 2.0D+00

  golden_angle = 2.0D+00 * r8_pi / golden_ratio ** 2

  do j = 1, n
    x(1,j) = c(1) + r * cos ( angle )
    x(2,j) = c(2) + r * sin ( angle )
    angle = mod ( angle + golden_angle, 2.0D+00 * r8_pi )
  end do

  return
end
subroutine circle01_sample_random ( n, x )

!*****************************************************************************80
!
!! circle01_sample_random() samples the circumference of the unit circle in 2D.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 September 2021
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
!    05 September 2021
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
!    05 September 2021
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
 
