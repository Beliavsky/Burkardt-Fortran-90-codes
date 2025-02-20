program main

!*****************************************************************************80
!
!! midpoint_adaptive_test() tests midpoint_adaptive().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 June 2024
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Catalin Trenchea, John Burkardt,
!    Refactorization of the midpoint rule,
!    Applied Mathematics Letters,
!    Volume 107, September 2020.
!
  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'midpoint_adaptive_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test midpoint_adaptive() on several ODE''s.'

  call lotka_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'midpoint_adaptive_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  return
end
subroutine lotka_test ( )

!*****************************************************************************80
!
!! lotka_test(): midpoint_adaptive() solves the lotka ODE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 July 2024
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    George Lindfield, John Penny,
!    Numerical Methods Using MATLAB,
!    Second Edition,
!    Prentice Hall, 1999,
!    ISBN: 0-13-012641-1,
!    LC: QA297.P45.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: nmax = 400
  integer, parameter :: m = 2

  real ( kind = rk8 ) abstol
  character ( len = 5 ) :: label = 'lotka'
  external lotka_deriv
  integer n_fsolve
  integer n_rejected
  integer nstep
  real ( kind = rk8 ) reltol
  real ( kind = rk8 ) t(0:nmax)
  real ( kind = rk8 ) t0
  real ( kind = rk8 ) tau0
  real ( kind = rk8 ) tmax
  real ( kind = rk8 ) y0(m)
  real ( kind = rk8 ) y(0:nmax,m)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'lotka_test():'
  write ( *, '(a)' ) '  midpoint_adaptive() solves the lotka ODE.'
  write ( *, '(a)' ) '  A pair of ordinary differential equations for a population'
  write ( *, '(a)' ) '  of predators and prey are solved using midpoint().'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  The exact solution shows periodic behavior, with a fixed'
  write ( *, '(a)' ) '  period and amplitude.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Steady state at r = 5000, f = 2000.'

  t0 = 0.0
  tmax = 20.0
  y0 = (/ 5000.0, 100.0 /)
  tau0 = 2.5D-2
  reltol = 1.0D-3
  abstol = 1.0D-3

  call mad_compute ( lotka_deriv, t0, tmax, m, y0, tau0, reltol, abstol, nmax, &
    nstep, n_rejected, n_fsolve, t, y )

  call mad_solution_plot ( nstep, m, t, y, label )

  call mad_phase_plot ( nstep, m, t, y, label )

  call mad_timestep_plot ( nstep, t, label )

  call mad_stats ( nstep, n_rejected, n_fsolve, t )

  return
end
subroutine lotka_deriv ( t, y, dydt )

!*****************************************************************************80
!
!! lotka_deriv() evaluates the right hand side of a Lotka-Volterra ODE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 July 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real T: the time;
!
!    real Y(2): the current solution values.
!
!  Output:
!
!    real DYDT(2): the derivative values.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) dydt(2)
  real ( kind = rk8 ) t
  real ( kind = rk8 ) y(2)
  
  call r8_fake_use ( t )

  dydt(1) =    2.0 * y(1) - 0.001 * y(1) * y(2)
  dydt(2) = - 10.0 * y(2) + 0.002 * y(1) * y(2)

  return
end
subroutine r8_fake_use ( x )

!*****************************************************************************80
!
!! r8_fake_use() pretends to use an R8 variable.
!
!  Discussion:
!
!    Some compilers will issue a warning if a variable is unused.
!    Sometimes there's a good reason to include a variable in a program,
!    but not to use it.  Calling this function with that variable as
!    the argument will shut the compiler up.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 April 2020
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real X, the variable to be "used".
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) x

  if ( x /= x ) then
    write ( *, '(a)' ) '  r8_fake_use(): variable is NAN.'
  end if

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! timestamp() prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 8 ) ampm
  integer d
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  integer values(8)
  integer y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
