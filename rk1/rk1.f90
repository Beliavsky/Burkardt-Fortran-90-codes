subroutine rk1 ( dydt, tspan, y0, n, m, t, y )

!*****************************************************************************80
!
!! rk1() solves an ODE using an explicit Runge-Kutta method of order 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 November 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    external dydt: a subroutine that evaluates the right
!    hand side of the ODE.
!
!    real ( kind = rk8 ) tspan(2): contains the initial and final times.
!
!    real ( kind = rk8 ) y0(m): a column vector containing the initial condition.
!
!    integer n: the number of steps to take.
!
!    integer m: the number of variables.
!
!  Output:
!
!    real ( kind = rk8 ) t(n+1), y(n+1,m): the times and solution values.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk8 ) dt
  real ( kind = rk8 ) dy(m)
  external dydt
  integer i
  real ( kind = rk8 ) t(n+1)
  real ( kind = rk8 ) tfirst
  real ( kind = rk8 ) tlast
  real ( kind = rk8 ) tspan(2)
  real ( kind = rk8 ) y(n+1,m)
  real ( kind = rk8 ) y0(m)

  tfirst = tspan(1)
  tlast = tspan(2)
  dt = ( tlast - tfirst ) / real ( n, kind = rk8 )
  t(1) = tspan(1)
  y(1,1:m) = y0(1:m)

  do i = 1, n
    t(i+1) = t(i) + dt
    call dydt ( t(i), y(i,1:m), dy )
    y(i+1,1:m) = y(i,1:m) + dt * dy(1:m)
  end do

  return
end

