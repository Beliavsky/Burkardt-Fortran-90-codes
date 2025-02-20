subroutine u_init ( alpha, beta, gamma, delta, m, n, x, y, u0 )

!*****************************************************************************80
!
!! U_INIT supplies the initial value of U at each node.
!
!  Modified:
!
!    01 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) ALPHA, BETA, GAMMA, DELTA, the parameters.
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = rk ) X(M,N), Y(M,N), the X and Y coordinates
!    of the nodes.
!
!    Output, real ( kind = rk ) U0(M,N), the initial value of U at the nodes.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) alpha
  real ( kind = rk ) beta
  real ( kind = rk ) delta
  real ( kind = rk ) gamma
  real ( kind = rk ) u0(m,n)
  real ( kind = rk ) ustar
  real ( kind = rk ) x(m,n)
  real ( kind = rk ) y(m,n)

  call r8_fake_use ( delta )

  ustar = gamma * alpha / ( beta - gamma )

  u0(1:m,1:n) = ustar - 2.0D-07 &
    * ( x(1:m,1:n) - 0.1D+00 * y(1:m,1:n) - 225.0D+00 ) &
    * ( x(1:m,1:n) - 0.1D+00 * y(1:m,1:n) - 675.0D+00 )

  return
end
subroutine v_init ( alpha, beta, gamma, delta, m, n, x, y, v0 )

!*****************************************************************************80
!
!! V_INIT supplies the initial value of U at each node.
!
!  Modified:
!
!    30 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) ALPHA, BETA, GAMMA, DELTA, the parameters.
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = rk ) X(M,N), Y(M,N), the X and Y coordinates
!    of the nodes.
!
!    Output, real ( kind = rk ) V0(M,N), the initial value of U at the nodes.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) alpha
  real ( kind = rk ) beta
  real ( kind = rk ) delta
  real ( kind = rk ) gamma
  real ( kind = rk ) ustar
  real ( kind = rk ) v0(m,n)
  real ( kind = rk ) vstar
  real ( kind = rk ) x(m,n)
  real ( kind = rk ) y(m,n)

  call r8_fake_use ( delta ) 

  ustar = gamma * alpha / ( beta - gamma )
  vstar = ( 1.0D+00 - ustar ) * ( alpha + ustar )

  v0(1:m,1:n) = vstar - 3.0D-05 * ( x(1:m,1:n) - 450.0D+00 ) &
    - 1.2D-04 * ( y(1:m,1:n) - 150.0D+00 )

  return
end
