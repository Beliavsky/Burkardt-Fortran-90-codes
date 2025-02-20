subroutine sine_gordon_exact ( a, x, y, u, uxy )

!*****************************************************************************80
!
!! sine_gordon_exact() evaluates an exact solution of the Sine-Gordeon PDE.
!
!  Discussion:
!
!    uxy = sin ( u )
!
!    This is a one-soliton solution.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 May 2024
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Daniel Arrigo,
!    Analytical Techniques for Solving Nonlinear Partial Differential Equations,
!    Morgan and Clayfoot, 2019,
!    ISBN: 978 168 173 5351.
!
!  Input:
!
!    real ( kind = rk ) A: a parameter.
!
!    real ( kind = rk ) X, Y: the X and Y coordinates of a point.
!
!  Output:
!
!    real ( kind = rk ) U, UXY: the values of the solution, and its first 
!    mixed derivative at (X,Y).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) u
  real ( kind = rk ) uxy
  real ( kind = rk ) x
  real ( kind = rk ) y

  u = 4.0 * atan ( exp ( a * x + y / a ) )

  uxy = 4.0 * ( exp ( a * x + y / a ) - exp ( 3 * a * x + 3 * y / a ) ) &
    / ( 1.0 + exp ( 2 * a * x + 2 * y / a ) ) ** 2

  return
end
subroutine sine_gordon_residual ( u, uxy, r )

!*****************************************************************************80
!
!! sine_gordon_residual() evaluates the residual for the Sine-Gordeon PDE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 May 2024
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Daniel Arrigo,
!    Analytical Techniques for Solving Nonlinear Partial Differential Equations,
!    Morgan and Clayfoot, 2019,
!    ISBN: 978 168 173 5351.
!
!  Input:
!
!    real ( kind = rk ) U, UXY: the values of the solution, and its
!    first mixed derivative.
!
!  Output:
!
!    real ( kind = rk ) R, the residual.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r
  real ( kind = rk ) u
  real ( kind = rk ) uxy
  r = uxy - sin ( u );

  return
end

