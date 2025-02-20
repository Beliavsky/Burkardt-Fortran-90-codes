program main

!*****************************************************************************80
!
!! ornstein_uhlenbeck_test() tests ornstein_uhlenbeck().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 January 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)') ' '
  write ( *, '(a)' ) 'ornstein_uhlenbeck_test():'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) '  Test ornstein_uhlenbeck().'

  call ou_euler_test ( )
  call ou_euler_maruyama_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'ornstein_uhlenbeck_test:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine ou_euler_test ( )

!*****************************************************************************80
!
!! OU_EULER_TEST tests OU_EULER.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 January 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) mu
  integer n
  real ( kind = rk ) sigma
  real ( kind = rk ) theta
  real ( kind = rk ) tmax
  real ( kind = rk ) x0

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'OU_EULER_TEST:'
  write ( *, '(a)' ) '  Estimate a solution to the Ornstein-Uhlenbeck equation'
  write ( *, '(a)' ) '  using the Euler method for stochastic differential equations.'
  write ( *, '(a)' ) ''

  theta = 2.0D+00
  write ( *, '(a,g14.6)' ) '  Using decay rate THETA = ', theta
  mu = 1.0D+00
  write ( *, '(a,g14.6)' ) '  Using mean MU = ', mu
  sigma = 0.15D+00
  write ( *, '(a,g14.6)' ) '  Using variance SIGMA = ', sigma
  x0 = 2.0D+00
  write ( *, '(a,g14.6)' ) '  Using initial value X0 = ', x0
  tmax = 3.0D+00
  write ( *, '(a,g14.6)' ) '  Using final time TMAX = ', tmax
  n = 10000
  write ( *, '(a,i8)' ) '  Using number of timesteps N = ', n

  call ou_euler ( theta, mu, sigma, x0, tmax, n )

  return
end
subroutine ou_euler_maruyama_test ( )

!*****************************************************************************80
!
!! OU_EULER_MARUYAMA_TEST tests OU_EULER_MARUYAMA.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 January 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) mu
  integer n
  integer r
  real ( kind = rk ) sigma
  real ( kind = rk ) theta
  real ( kind = rk ) tmax
  real ( kind = rk ) x0

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'OU_EULER_MARUYAMA_TEST:'
  write ( *, '(a)' ) '  Estimate a solution to the Ornstein-Uhlenbeck equation'
  write ( *, '(a)' ) '  using the Euler-Maruyama method for stochastic '
  write ( *, '(a)' ) '  differential equations.'
  write ( *, '(a)' ) ''

  theta = 2.0D+00
  write ( *, '(a,g14.6)' ) '  Using decay rate THETA = ', theta
  mu = 1.0D+00
  write ( *, '(a,g14.6)' ) '  Using mean MU = ', mu
  sigma = 0.15D+00
  write ( *, '(a,g14.6)' ) '  Using variance SIGMA = ', sigma
  x0 = 2.0D+00
  write ( *, '(a,g14.6)' ) '  Using initial value X0 = ', x0
  tmax = 3.0D+00
  write ( *, '(a,g14.6)' ) '  Using final time TMAX = ', tmax
  n = 10000
  write ( *, '(a,i6)' ) '  Using number of large timesteps N = ', n
  r = 16
  write ( *, '(a,i6)' ) '  Using number small time steps per one large time step R = ', r

  call ou_euler_maruyama ( theta, mu, sigma, x0, tmax, n, r )

  return
end

