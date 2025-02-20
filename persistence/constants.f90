subroutine constants ( g, c, h )

!*****************************************************************************80
!
!! constants() stores, and returns constants "g", "c" and "h".
!
!  Discussion:
!
!    Calling [g,c,h]=constants() returns the values of g, c, and h.
!
!    Because the values never change, and don't need to be computed,
!    we use assignment statements here, instead of persistent data.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 May 2021
!
!  Author:
!
!    John Burkardt
!
!  Output:
!
!    real ( kind = rk ) G: gravitational constant m^3/s^2/kg
!
!    real ( kind = rk ) C: light speed, m/s.
!
!    real ( kind = rk ) H: Planck's constant, j s;
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) c
  real ( kind = rk ) g
  real ( kind = rk ) h

  g = 6.67384D-11;
  c = 2.99792458D+8;
  h = 6.626070040D-34;

  return
end

