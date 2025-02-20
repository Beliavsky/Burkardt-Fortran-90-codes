subroutine dirichlet_condition ( node_num, node_xy, time, u )

!*****************************************************************************80
!
!! DIRICHLET_CONDITION sets the value of a Dirichlet boundary condition.
!
!  Discussion:
!
!    The input points (X,Y) are assumed to lie on the boundary of the
!    region.
!
!    This routine is set for the unit square.
!
!    We assume that the equation to be solved is
!
!      dUdT - Laplacian U + K * U = F
!
!    with K = 0, and F = (2*pi*pi-1)*sin(pi*x)*sin(pi*y)*exp(-t).
!
!    The exact solution is:
!
!      U = sin(pi*x) * sin(pi*y) * exp(-t)
!
!  Modified:
!
!    08 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of points.
!
!    Input, real ( kind = rk ) NODE_XY(2,NODE_NUM), 
!    the coordinates of the points.
!
!    Input, real ( kind = rk ) TIME, the current time.
!
!    Output, real ( kind = rk ) U(NODE_NUM), the value of the solution at the
!    the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer node_num

  real ( kind = rk ) node_xy(2,node_num)
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) time
  real ( kind = rk ) u(node_num)

  u(1:node_num) = sin ( pi * node_xy(1,1:node_num) ) &
                * sin ( pi * node_xy(2,1:node_num) ) * exp ( - time )

  return
end
subroutine initial_condition ( node_num, node_xy, time, u )

!*****************************************************************************00
!
!! INITIAL_CONDITION sets the initial condition.
!
!  Discussion:
!
!    The input value TIME is assumed to be the initial time.
!
!    We assume that the equation to be solved is
!
!      dUdT - Laplacian U + K * U = F
!
!    with K = 0, and F = (2*pi*pi-1)*sin(pi*x)*sin(pi*y)*exp(-t).
!
!    The exact solution is:
!
!      U = sin(pi*x) * sin(pi*y) * exp(-t)
!
!  Modified:
!
!    08 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of points.
!
!    Input, real ( kind = rk ) NODE_XY(2,NODE_NUM), 
!    the coordinates of the points.
!
!    Input, real ( kind = rk ) TIME, the current time (assumed to be
!    the initial time).
!
!    Output, real ( kind = rk ) U(NODE_NUM), the value of the solution at the
!    the initial time.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer node_num

  real ( kind = rk ) node_xy(2,node_num)
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) time
  real ( kind = rk ) u(node_num)

  u(1:node_num) = sin ( pi * node_xy(1,1:node_num) ) &
                * sin ( pi * node_xy(2,1:node_num) ) * exp ( - time )

  return
end
subroutine k_coef ( node_num, node_xy, time, k )

!*****************************************************************************80
!
!! K_COEF evaluates the coefficient K(X,Y,T) function.
!
!  Discussion:
!
!    We assume that the equation to be solved is
!
!      dUdT - Laplacian U + K * U = F
!
!    with K = 0, and F = (2*pi*pi-1)*sin(pi*x)*sin(pi*y)*exp(-t).
!
!    The exact solution is:
!
!      U = sin(pi*x) * sin(pi*y) * exp(-t)
!
!  Modified:
!
!    08 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of points.
!
!    Input, real ( kind = rk ) NODE_XY(2,NODE_NUM), 
!    the coordinates of the points.
!
!    Input, real ( kind = rk ) TIME, the current time (assumed to be
!    the initial time).
!
!    Output, real ( kind = rk ) K(NODE_NUM), the value of the coefficient.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer node_num

  real ( kind = rk ) k(node_num)
  real ( kind = rk ) node_xy(2,node_num)
  real ( kind = rk ) time

  call r8_fake_use ( node_xy(1,1) )
  call r8_fake_use ( time )

  k(1:node_num) = 0.0D+00

  return 
end
subroutine rhs ( node_num, node_xy, time, f )

!*****************************************************************************80
!
!! RHS gives the right-hand side of the differential equation.
!
!  Discussion:
!
!    We assume that the equation to be solved is
!
!      dUdT - Laplacian U + K * U = F
!
!    with 
!
!      K = 0, 
!
!    and 
!
!      F = (2*pi*pi-1)*sin(pi*x)*sin(pi*y)*exp(-t).
!
!    The exact solution is:
!
!      U = sin(pi*x) * sin(pi*y) * exp(-t)
!
!  Modified:
!
!    08 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of points.
!
!    Input, real ( kind = rk ) NODE_XY(2,NODE_NUM), 
!    the coordinates of the points.
!
!    Input, real ( kind = rk ) TIME, the current time (assumed to be
!    the initial time).
!
!    Output, real ( kind = rk ) F(NODE_NUM), the value of the right hand side.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer node_num

  real ( kind = rk ) f(node_num)
  real ( kind = rk ) node_xy(2,node_num)
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) time

  f(1:node_num) = ( 2.0D+00 * pi * pi - 1.0D+00 ) &
    * sin ( pi * node_xy(1,1:node_num) ) &
    * sin ( pi * node_xy(2,1:node_num) ) * exp ( - time )

  return 
end
