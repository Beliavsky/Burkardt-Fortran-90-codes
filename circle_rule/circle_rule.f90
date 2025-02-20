subroutine circle_rule ( nt, w, t )

!*****************************************************************************80
!
!! circle_rule() computes a quadrature rule for the unit circle.
!
!  Discussion:
!
!    The unit circle is the region:
!
!      x * x + y * y = 1.
!
!    The integral I(f) is then approximated by
!
!      Q(f) = 2 * pi * sum ( 1 <= i <= NT ) W(i) * F ( cos(T(i)), sin(T(i)) ).
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
!    integer NT, the number of angles to use.
!
!  Output:
!
!    real ( kind = rk ) W(NT), the weights for the rule.
!
!    real ( kind = rk ) T(NT), the angles for the rule.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nt

  integer it
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) t(nt)
  real ( kind = rk ) w(nt)

  w(1:nt) = 1.0D+00 / real ( nt, kind = rk )
  do it = 1, nt
    t(it) = 2.0D+00 * r8_pi * real ( it - 1, kind = rk ) / real ( nt, kind = rk )
  end do

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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
