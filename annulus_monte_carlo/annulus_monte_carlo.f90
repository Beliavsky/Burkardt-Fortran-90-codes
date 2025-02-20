function annulus_area ( r1, r2 )

!*****************************************************************************80
!
!! annulus_area() computes the area of a circular annulus in 2D.
!
!  Discussion:
!
!    A circular annulus with center (XC,YC), inner radius R1 and
!    outer radius R2, is the set of points (X,Y) so that
!
!      R1^2 <= (X-XC)^2 + (Y-YC)^2 <= R2^2
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
!    real ( kind = rk ) R1, R2, the inner and outer radii of the annulus.
!    0.0 <= R1 <= R2.
!
!  Output:
!
!    real ( kind = rk ) ANNULUS_AREA, the area of the annulus.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) annulus_area
  real ( kind = rk ) r1
  real ( kind = rk ) r2
  real ( kind = rk ), parameter :: r8_pi = 3.1415926535897932384626434D+00
  real ( kind = rk ) value

  if ( r1 < 0.0_rk ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'ANNULUS_AREA - Fatal error!'
    write ( *, '(a)' ) '  Inner radius R1 < 0.0.'
    stop 1
  end if

  if ( r2 < r1 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'ANNULUS_AREA - Fatal error!'
    write ( *, '(a)' ) '  Outer radius R1 < R1 = inner radius.'
    stop 1
  end if

  value = r8_pi * ( r2 + r1 ) * ( r2 - r1 )

  annulus_area = value

  return
end
subroutine annulus_sample ( center, r1, r2, n, p )

!*****************************************************************************80
!
!! annulus_sample() samples a circular annulus.
!
!  Discussion:
!
!    A circular annulus with center PC, inner radius R1 and
!    outer radius R2, is the set of points P so that
!
!      R1^2 <= (P(1)-PC(1))^2 + (P(2)-PC(2))^2 <= R2^2
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
!  Reference:
!
!    Peter Shirley,
!    Nonuniform Random Point Sets Via Warping,
!    Graphics Gems, Volume III,
!    edited by David Kirk,
!    AP Professional, 1992, 
!    ISBN: 0122861663,
!    LC: T385.G6973.
!
!  Input:
!
!    real ( kind = rk ) CENTER(2), the center.
!
!    real ( kind = rk ) R1, R2, the inner and outer radii.
!    0.0 <= R1 <= R2.
!
!    integer N, the number of points to generate.
!
!  Output:
!
!    real ( kind = rk ) P(2,N), sample points.
!
  implicit none

  integer n
  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) center(2)
  real ( kind = rk ) p(2,n)
  real ( kind = rk ) r(n)
  real ( kind = rk ) r1
  real ( kind = rk ) r2
  real ( kind = rk ), parameter :: r8_pi = 3.1415926535897932384626434D+00
  real ( kind = rk ) theta(n)
  real ( kind = rk ) u(n)
  real ( kind = rk ) v(n)

  if ( r1 < 0.0D+00 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'ANNULUS_SAMPLE - Fatal error!'
    write ( *, '(a)' ) '  Inner radius R1 < 0.0.'
    stop 1
  end if

  if ( r2 < r1 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'ANNULUS_SAMPLE - Fatal error!'
    write ( *, '(a)' ) '  Outer radius R1 < R1 = inner radius.'
    stop 1
  end if

  call random_number ( harvest = u )

  theta(1:n) = u(1:n) * 2.0_rk * r8_pi

  call random_number ( harvest = v )

  r(1:n) = sqrt ( ( 1.0_rk - v(1:n) ) * r1 * r1 &
       +                     v(1:n)   * r2 * r2 )

  p(1,1:n) = center(1) + r(1:n) * cos ( theta(1:n) )
  p(2,1:n) = center(2) + r(1:n) * sin ( theta(1:n) )

  return
end
subroutine disk01_monomial_integral ( e, integral )

!*****************************************************************************80
!
!! disk01_monomial_integral() returns monomial integrals in the unit disk in 2D.
!
!  Discussion:
!
!    The integration region is 
!
!      X^2 + Y^2 <= 1.
!
!    The monomial is F(X,Y) = X^E(1) * Y^E(2).
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

  real ( kind = rk ) arg
  integer e(2)
  integer i
  real ( kind = rk ) integral
  real ( kind = rk ), parameter :: r = 1.0D+00
  integer s

  if ( any ( e(1:2) < 0 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DISK01_MONOMIAL_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  All exponents must be nonnegative.'
    stop 1
  end if

  if ( any ( mod ( e(1:2), 2 ) == 1 ) ) then

    integral = 0.0_rk

  else

    integral = 2.0_rk

    do i = 1, 2
      arg = 0.5_rk * real ( e(i) + 1, kind = rk )
      integral = integral * gamma ( arg )
    end do

    arg = 0.5_rk * ( real ( sum ( e(1:2) + 1 ), kind = rk ) )
    integral = integral / gamma ( arg )

  end if
!
!  The surface integral is now adjusted to give the volume integral.
!
  s = sum ( e(1:2) ) + 2

  integral = integral * r ** s / real ( s, kind = rk )

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
!    02 January 2014
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

  integer m
  integer n
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer e(m)
  integer i
  real ( kind = rk ) value(n)
  real ( kind = rk ) x(m,n)

  value(1:n) = 1.0_rk

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
!    02 January 2014
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

