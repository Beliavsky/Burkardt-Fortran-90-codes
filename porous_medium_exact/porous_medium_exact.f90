subroutine porous_medium_exact ( x, t, u, ut, ux, uxx )

!*****************************************************************************80
!
!! porous_medium_exact() evaluates an exact solution of the porous medium equation.
!
!  Discussion:
!
!    du/dt = Del**2 ( u**m )
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 May 2024
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Grigory Barenblatt,
!    On some unsteady fluid and gas motions in a porous medium,
!    Prikladnaya Matematika i Mekhanika (Applied Mathematics and Mechanics, 
!    Volume 16, Number 1, pages 67-78, 1952.
!
!    Rouben Rostamian,
!    Programming Projects in C 
!    for Students of Engineering, Science, and Mathematics,
!    SIAM, 2014,
!    ISBN: 978-1-611973-49-5
!
!  Input:
!
!    real ( kind = rk8 ) x, t: the position and time.
!
!  Output:
!
!    real ( kind = rk8 ) u, ut, ux, uxx: the values of the exact solution, its 
!    time derivative, and its first and second spatial derivatives at (x,t).
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  interface
    subroutine porous_medium_parameters ( c_in, delta_in, m_in, t0_in, &
      tstop_in, c_out, delta_out, m_out, t0_out, tstop_out )
      integer, parameter :: rk8 = kind ( 1.0D+00 )
      real ( kind = rk8 ), optional :: c_in
      real ( kind = rk8 ), optional :: c_out
      real ( kind = rk8 ), optional :: delta_in
      real ( kind = rk8 ), optional :: delta_out
      real ( kind = rk8 ), optional :: m_in
      real ( kind = rk8 ), optional :: m_out
      real ( kind = rk8 ), optional :: t0_in
      real ( kind = rk8 ), optional :: t0_out
      real ( kind = rk8 ), optional :: tstop_in
      real ( kind = rk8 ), optional :: tstop_out
    end subroutine
  end interface

  real ( kind = rk8 ) alpha
  real ( kind = rk8 ) beta
  real ( kind = rk8 ) bot
  real ( kind = rk8 ) c
  real ( kind = rk8 ) delta
  real ( kind = rk8 ) factor
  real ( kind = rk8 ) gamma
  real ( kind = rk8 ) m
  real ( kind = rk8 ) t
  real ( kind = rk8 ) u
  real ( kind = rk8 ) ut
  real ( kind = rk8 ) ux
  real ( kind = rk8 ) uxx
  real ( kind = rk8 ) x

  call porous_medium_parameters ( c_out = c, delta_out = delta, m_out = m )

  alpha =         1.0D+00 / ( m - 1.0D+00 )
  beta =          1.0D+00 / ( m + 1.0D+00 )
  gamma = ( m - 1.0D+00 ) / ( 2.0D+00 * m * ( m + 1.0D+00 ) )

  bot = ( t + delta )**beta
  factor = c - gamma * ( x / bot )**2

  if ( 0.0D+00 < factor ) then

    u = 1.0D+00 / ( t + delta)**beta * factor**alpha

    ut = 2.0D+00 * alpha * beta * gamma * ( t + delta )**(-1.0D+00-3.0D+00*beta) &
      * x**2 * factor**(alpha - 1.0D+00) &
      - beta * ( t + delta )**(-1.0D+00-beta) * factor**alpha

    ux = - 2.0D+00 * alpha * gamma &
      * ( t + delta )**(-3.0D+00 * beta ) * x * factor**(alpha-1.0D+00)

    uxx = 4.0 * ( alpha - 1.0D+00 ) * alpha * gamma**2 &
      * ( t + delta )**(-5.0D+00*beta) * x**2 * factor**(alpha-2.0D+00) &
      - 2.0 * alpha * gamma * ( t + delta )**(-3.0D+00 * beta ) &
      * factor**(alpha-1.0D+00)

  else

    u = 0.0D+00
    ut = 0.0D+00
    ux = 0.0D+00
    uxx = 0.0D+00

  end if

  return
end
subroutine porous_medium_parameters ( c_in, delta_in, m_in, t0_in, tstop_in, &
  c_out, delta_out, m_out, t0_out, tstop_out )

!*****************************************************************************80
!
!! porous_medium_parameters() returns parameters for the porous medium equation.
!
!  Discussion:
!
!    du/dt = Del**2 ( u**m )
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
!    19 May 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk8 ) c_in: a linear factor for the volume of the solution.
!
!    real ( kind = rk8 ) delta_in: an offset to the time t.
!
!    real ( kind = rk8 ) m_in: the power of u in the equation.
!    m must be greater than 1.
!
!    real ( kind = rk8 ) t0_in: the initial time.
!
!    real ( kind = rk8 ) tstop_in: the final time.
!
!  Output:
!
!    real ( kind = rk8 ) c_out: a linear factor for the volume of the solution. 
!
!    real ( kind = rk8 ) delta_out: an offset to the time t.
!
!    real ( kind = rk8 ) m_out: the power of u in the equation.  
!    m must be greater than 1.
!
!    real ( kind = rk8 ) t0_out: the initial time.
!
!    real ( kind = rk8 ) tstop_out: the final time.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ), save     :: c_default = sqrt ( 3.0D+00 ) / 15.0D+00
  real ( kind = rk8 ), optional :: c_in
  real ( kind = rk8 ), optional :: c_out
  real ( kind = rk8 ), save     :: delta_default = 1.0 / 75.0
  real ( kind = rk8 ), optional :: delta_in
  real ( kind = rk8 ), optional :: delta_out
  real ( kind = rk8 ), save     :: m_default = 3.0D+00
  real ( kind = rk8 ), optional :: m_in
  real ( kind = rk8 ), optional :: m_out
  real ( kind = rk8 ), save     :: t0_default = 0.0D+00
  real ( kind = rk8 ), optional :: t0_in
  real ( kind = rk8 ), optional :: t0_out
  real ( kind = rk8 ), save     :: tstop_default = 4.0D+00
  real ( kind = rk8 ), optional :: tstop_in
  real ( kind = rk8 ), optional :: tstop_out
!
!  Update defaults if input was supplied.
!
  if ( present ( c_in ) ) then
    c_default = c_in
  end if

  if ( present ( delta_in ) ) then
    delta_default = delta_in
  end if

  if ( present ( m_in ) ) then
    m_default = m_in
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
  if ( present ( c_out ) ) then
    c_out = c_default
  end if

  if ( present ( delta_out ) ) then
    delta_out = delta_default
  end if

  if ( present ( m_out ) ) then
    m_out = m_default
  end if

  if ( present ( t0_out ) ) then
    t0_out = t0_default
  end if

  if ( present ( tstop_out ) ) then
    tstop_out = tstop_default
  end if

  return
end
subroutine porous_medium_residual ( t, x, r )

!*****************************************************************************80
!
!! porous_medium_residual() computes the residual of the porous medium equation.
!
!  Discussion:
!
!    ut = Del**2 ( u**m )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 May 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk8 ) T, X: the time and position where the solution
!    is evaluated.
!
!    real ( kind = rk8 ) M:
!
!  Output:
!
!    real ( kind = rk8 ) R: the residual at that time and position.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  interface
    subroutine porous_medium_parameters ( c_in, delta_in, m_in, t0_in, &
      tstop_in, c_out, delta_out, m_out, t0_out, tstop_out )
      integer, parameter :: rk8 = kind ( 1.0D+00 )
      real ( kind = rk8 ), optional :: c_in
      real ( kind = rk8 ), optional :: c_out
      real ( kind = rk8 ), optional :: delta_in
      real ( kind = rk8 ), optional :: delta_out
      real ( kind = rk8 ), optional :: m_in
      real ( kind = rk8 ), optional :: m_out
      real ( kind = rk8 ), optional :: t0_in
      real ( kind = rk8 ), optional :: t0_out
      real ( kind = rk8 ), optional :: tstop_in
      real ( kind = rk8 ), optional :: tstop_out
    end subroutine
  end interface

  real ( kind = rk8 ) m
  real ( kind = rk8 ) r
  real ( kind = rk8 ) t
  real ( kind = rk8 ) u
  real ( kind = rk8 ) ut
  real ( kind = rk8 ) ux
  real ( kind = rk8 ) uxx
  real ( kind = rk8 ) x

  call porous_medium_parameters ( m_out = m )

  call porous_medium_exact ( t, x, u, ut, ux, uxx )

  r = ut - m * ( m - 1 ) * u**(m-2) * ux**2 - m * u**(m-1) * uxx

  return
end
