subroutine midpoint_fixed ( dydt, tspan, y0, n, m, t, y )

!*****************************************************************************80
!
!! midpoint_fixed() uses a fixed-point midpoint method to solve an ODE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 September 2023
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
!    real ( kind = rk ) tspan(2): contains the initial and final times.
!
!    real ( kind = rk ) y0(m): a column vector containing the initial condition.
!
!    integer n: the number of steps to take.
!
!    integer m: the number of variables.
!
!  Output:
!
!    real ( kind = rk ) t(n+1), y(n+1,m): the times and solution values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) dt
  external dydt
  real ( kind = rk ) f(m)
  integer i
  integer it_max
  integer j
  real ( kind = rk ) t(n+1)
  real ( kind = rk ) theta
  real ( kind = rk ) tm
  real ( kind = rk ) tspan(2)
  real ( kind = rk ) y(n+1,m)
  real ( kind = rk ) y0(m)
  real ( kind = rk ) ym(m)

  dt = ( tspan(2) - tspan(1) ) / n

  it_max = 10
  theta = 0.5D+00

  t(1) = tspan(1)
  y(1,1:m) = y0(1:m)

  do i = 1, n
    tm = t(i) + theta * dt 
    ym(1:m) = y(i,1:m)
    do j = 1, it_max
      call dydt ( tm, ym(1:m), f )
      ym(1:m) = y(i,1:m) + theta * dt * f(1:m)
    end do
    t(i+1) = t(i) + dt
    y(i+1,1:m) = (           1.0D+00 / theta ) * ym(1:m) &
               + ( 1.0D+00 - 1.0D+00 / theta ) * y(i,1:m)
  end do

  return
end

