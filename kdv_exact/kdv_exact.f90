subroutine kdv_exact_rational ( x, t, u, ut, ux, uxx, uxxx )

!*****************************************************************************80
!
!! kdv_exact_rational(): exact solution of KDV PDE using a rational function.
!
!  Discussion:
!
!    This solution u(x,t) satisfies the Korteweg-Devries PDE:
!      u' - 6 u ux + uxxx = 0
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 April 2024
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    John D Cook,
!    Rational solution to the Korteweg-De Vries equation,
!    13 November 2023,
!    https://www.johndcook.com/blog/2023/11/13/rational-kdv/
!
!  Input:
!
!    real ( kind = rk ) x: the position.
!
!    real ( kind = rk ) t: the time.
!
!  Output:
!
!    real ( kind = rk ) u, ut, ux, uxx, uxxx:
!    the values of the solution and its derivatives.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) t
  real ( kind = rk ) u
  real ( kind = rk ) ut
  real ( kind = rk ) ux
  real ( kind = rk ) uxx
  real ( kind = rk ) uxxx
  real ( kind = rk ) x

  u = 6.0 * x * ( x**3 - 24.0 * t ) &
    / ( x**3 + 12.0 * t )**2

  ut = - 288.0 * x * ( x**3 - 6.0 * t ) &
    / ( x**3 + 12.0 * t )**3

  ux = - 12.0 * ( x**6 - 84.0 * t * x**3 + 144.0 * t**2 ) &
    / ( x**3 + 12.0 * t )**3

  uxx =  36.0 * ( x**8 - 192.0 * t * x**5 + 1440.0 * t**2 * x**2 ) &
    / ( x**3 + 12.0 * t )**4

  uxxx = - 144.0 * x * ( x**9 - 360.0 * t * x**6 &
    + 6480.0 * t**2 * x**3 - 8640.0 * t**3 ) &
    / ( x**3 + 12.0 * t )**5

  return
end
subroutine kdv_exact_sech ( x, t, u, ut, ux, uxx, uxxx )

!*****************************************************************************80
!
!! kdv_exact_sech(): exact solution of KDV PDE using sech().
!
!  Discussion:
!
!    This solution u(x,t) satisfies the Korteweg-Devries PDE:
!      u' - 6 u ux + uxxx = 0
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 April 2024
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    John D Cook,
!    Solitons and the KdV equation,
!    03 November 2023,
!    https://www.johndcook.com/blog/2023/11/03/solitons-and-the-kdv-equation/
!
!  Input:
!
!    real ( kind = rk ) x: the position.
!
!    real ( kind = rk ) t: the time.
!
!  Output:
!
!    real ( kind = rk ) u, ut, ux, uxx, uxxx:
!    the values of the solution and its derivatives.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  interface
    subroutine kdv_parameters ( a_in, v_in, t0_in, tstop_in, &
      a_out, v_out, t0_out, tstop_out )
      integer, parameter :: rk = kind ( 1.0D+00 )
      real ( kind = rk ), optional :: a_in
      real ( kind = rk ), optional :: a_out
      real ( kind = rk ), optional :: t0_in
      real ( kind = rk ), optional :: t0_out
      real ( kind = rk ), optional :: tstop_in
      real ( kind = rk ), optional :: tstop_out
      real ( kind = rk ), optional :: v_in
      real ( kind = rk ), optional :: v_out
    end subroutine
  end interface

  real ( kind = rk ) a
  real ( kind = rk ) argument
  real ( kind = rk ) sa
  real ( kind = rk ) t
  real ( kind = rk ) ta
  real ( kind = rk ) u
  real ( kind = rk ) ut
  real ( kind = rk ) ux
  real ( kind = rk ) uxx
  real ( kind = rk ) uxxx
  real ( kind = rk ) v
  real ( kind = rk ) x
!
!  Retrieve parameters a and v.
!
  call kdv_parameters ( a_out = a, v_out = v )

  argument = 0.5 * sqrt ( v ) * ( x - v * t - a )
  sa = 1.0 / cosh ( argument )
  ta = tanh ( argument )

  u =    - 0.5  * v      * sa**2

  ut =   - 0.5  * v**2.5 * sa**2 * ta

  ux =     0.5  * v**1.5 * sa**2 * ta

  uxx =  + 0.25 * v**2   * sa**4 &
         - 0.5  * v**2   * sa**2 * ta**2

  uxxx =        - v**2.5 * sa**4 * ta &
         + 0.5  * v**2.5 * sa**2 * ta**3

  return
end
subroutine kdv_parameters ( a_in, v_in, t0_in, tstop_in, a_out, v_out, &
  t0_out, tstop_out )

!*****************************************************************************80
!
!! kdv_parameters(): parameters for the Korteweg-Devries PDE.
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
!    30 April 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) a_user: the phase
!
!    real ( kind = rk ) v_user: the velocity
!
!    real ( kind = rk ) t0_user: the initial time
!
!    real ( kind = rk ) tstop_user: the final time.
!
!  Output:
!
!    real ( kind = rk ) a_out: the phase.
!
!    real ( kind = rk ) v_out: the velocity
!
!    real ( kind = rk ) t0_out: the initial time
!
!    real ( kind = rk ) tstop_out: the final time.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), save     :: a_default = 0.0D+00
  real ( kind = rk ), optional :: a_in
  real ( kind = rk ), optional :: a_out
  real ( kind = rk ), save     :: t0_default = 0.0D+00
  real ( kind = rk ), optional :: t0_in
  real ( kind = rk ), optional :: t0_out
  real ( kind = rk ), save     :: tstop_default = 10.0D+00
  real ( kind = rk ), optional :: tstop_in
  real ( kind = rk ), optional :: tstop_out
  real ( kind = rk ), save     :: v_default = 1.0D+00
  real ( kind = rk ), optional :: v_in
  real ( kind = rk ), optional :: v_out
!
!  New values, if supplied on input, overwrite the current values.
!
  if ( present ( a_in ) ) then
    a_default = a_in
  end if

  if ( present ( v_in ) ) then
    v_default = v_in
  end if

  if ( present ( t0_in ) ) then
    t0_default = t0_in
  end if

  if ( present ( tstop_in ) ) then
    tstop_default = tstop_in
  end if
!
!  The current values are copied to the output.
!
  if ( present ( a_out ) ) then
    a_out = a_default
  end if

  if ( present ( v_out ) ) then
    v_out = v_default
  end if

  if ( present ( t0_out ) ) then
    t0_out = t0_default
  end if

  if ( present ( tstop_out ) ) then
    tstop_out = tstop_default
  end if
  
  return
end
subroutine kdv_residual ( u, ut, ux, uxxx, r )

!*****************************************************************************80
!
!! kdv_residual(): evaluate the residual of KDV PDE.
!
!  Discussion:
!
!    The residual of the Korteweg-Devries PDE:
!      u' - 6 u ux + uxxx
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 April 2024
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    John D Cook,
!    Solitons and the KdV equation,
!    03 November 2023,
!    https://www.johndcook.com/blog/2023/11/03/solitons-and-the-kdv-equation/
!
!  Input:
!
!    real ( kind = rk ) u, ut, ux, uxxx:
!    the values of the solution and its derivatives.
!
!  Output:
!
!    real ( kind = rk ) r: the residual.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r
  real ( kind = rk ) u
  real ( kind = rk ) ut
  real ( kind = rk ) ux
  real ( kind = rk ) uxxx

  r = ut - 6.0 * u * ux + uxxx

  return
end

