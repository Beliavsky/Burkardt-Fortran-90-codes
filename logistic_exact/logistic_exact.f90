subroutine logistic_exact ( n, t, y )

!*****************************************************************************80
!
!! logistic_exact() evaluates the exact solution for logistic_ode().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 May 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer n: the number of evaluation points.
!
!    real ( kind = rk8 ) t(n): the evaluation points.
!
!  Output:
!
!    real ( kind = rk8 ) y(n): the function values.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  interface
    subroutine logistic_parameters ( r_in, k_in, t0_in, y0_in, tstop_in, &
      r_out, k_out, t0_out, y0_out, tstop_out )
      integer, parameter :: rk8 = kind ( 1.0D+00 )
      real ( kind = rk8 ), optional :: r_in
      real ( kind = rk8 ), optional :: r_out
      real ( kind = rk8 ), optional :: k_in
      real ( kind = rk8 ), optional :: k_out
      real ( kind = rk8 ), optional :: t0_in
      real ( kind = rk8 ), optional :: t0_out
      real ( kind = rk8 ), optional :: y0_in
      real ( kind = rk8 ), optional :: y0_out
      real ( kind = rk8 ), optional :: tstop_in
      real ( kind = rk8 ), optional :: tstop_out
    end subroutine
  end interface

  integer n

  real ( kind = rk8 ) k
  real ( kind = rk8 ) r
  real ( kind = rk8 ) t(n)
  real ( kind = rk8 ) t0
  real ( kind = rk8 ) y(n)
  real ( kind = rk8 ) y0

  call logistic_parameters ( r_out = r, k_out = k, t0_out = t0, y0_out = y0 )

  y = ( k * y0 * exp ( r * ( t - t0 ) ) ) &
    / ( k + y0 * ( exp ( r * ( t - t0 ) ) - 1.0D+00 ) )

  return
end
subroutine logistic_parameters ( r_in, k_in, t0_in, y0_in, tstop_in, r_out, &
  k_out, t0_out, y0_out, tstop_out )

!*****************************************************************************80
!
!! logistic_parameters() returns parameters for the logistic ODE.
!
!  Discussion:
!
!    If input values are specified, this resets the default parameters.
!    Otherwise, the output will be the current defaults.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 May 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real r_in: the rate of increase.
!
!    real k_in: the carrying capacity.
!
!    real t0_in: the initial time.
!
!    real y0_in: the initial condition at time T0.
!
!    real tstop_in: the final time.
!
!  Output:
!
!    real r_out: the rate of increase..
!
!    real k_out: the carrying capacity..
!
!    real t0_out: the initial time.
!
!    real y0_out: the initial condition at time T0.
!
!    real tstop_out: the final time.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ), save     :: r_default = 1.0D+00
  real ( kind = rk8 ), optional :: r_in
  real ( kind = rk8 ), optional :: r_out
  real ( kind = rk8 ), save     :: k_default = 1.0D+00
  real ( kind = rk8 ), optional :: k_in
  real ( kind = rk8 ), optional :: k_out
  real ( kind = rk8 ), save     :: t0_default = 0.0D+00
  real ( kind = rk8 ), optional :: t0_in
  real ( kind = rk8 ), optional :: t0_out
  real ( kind = rk8 ), save     :: y0_default = 0.5D+00
  real ( kind = rk8 ), optional :: y0_in
  real ( kind = rk8 ), optional :: y0_out
  real ( kind = rk8 ), save     :: tstop_default = 8.0D+00
  real ( kind = rk8 ), optional :: tstop_in
  real ( kind = rk8 ), optional :: tstop_out
!
!  New values, if supplied on input, overwrite the current values.
!
  if ( present ( r_in ) ) then
    r_default = r_in
  end if

  if ( present ( k_in ) ) then
    k_default = k_in
  end if

  if ( present ( t0_in ) ) then
    t0_default = t0_in
  end if

  if ( present ( y0_in ) ) then
    y0_default = y0_in
  end if

  if ( present ( tstop_in ) ) then
    tstop_default = tstop_in
  end if
!
!  The current values are copied to the output.
!
  if ( present ( r_out ) ) then
    r_out = r_default
  end if

  if ( present ( k_out ) ) then
    k_out = k_default
  end if

  if ( present ( t0_out ) ) then
    t0_out = t0_default
  end if

  if ( present ( y0_out ) ) then
    y0_out = y0_default
  end if

  if ( present ( tstop_out ) ) then
    tstop_out = tstop_default
  end if

  return
end

