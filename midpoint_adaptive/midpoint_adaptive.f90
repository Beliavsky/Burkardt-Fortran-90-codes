subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! get_unit() returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free FORTRAN unit.  The code assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Output:
!
!    integer IUNIT, the free unit number.
!
  implicit none

  integer i
  integer ios
  integer iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
subroutine mad_compute ( f, t0, tmax, m, y0, tau0, reltol, abstol, nmax, &
  nstep, n_rejected, n_fsolve, t, y )

!*****************************************************************************80
!
!! mad_compute() estimates the solution of an ordinary differential equation.
!
!  Discussion:
!
!    mad_compute() uses an adaptive time stepping algorithm for the
!    implicit midpoint method, applied to a system of ordinary differential
!    equations.
!
!    The implicit midpoint method has local truncation error O(h^3).
! 
!    The time step is adaptively controlled with respect to relative and
!    absolute tolerances applied to the estimated local truncation error (LTE).
!    
!    The adaptivity with respect to relative and absolute error tolerances
!    is described on page 168 in the Hairer reference.
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
!    John Burkardt, Catalin Trenchea
!
!  Reference:
!
!    Catalin Trenchea, John Burkardt,
!    Refactorization of the midpoint rule,
!    Applied Mathematics Letters,
!    Volume 107, September 2020.
!
!    Ernst Hairer, Syvert Norsett, Gerhard Wanner,
!    Solving ordinary differential equations, I. Nonstiff problems, 
!    Springer Series in Computational Mathematics, Number 8. 
!    Springer-Verlag, Berlin, 1987.
!
!  Input:
!
!    external f: the function that evaluates the right hand side 
!    of the differential equation, of the form 
!      subroutine f ( t, y, dydt )
!
!    real t0, tmax: the interval of integration.
!
!    integer m: the dimension of a single solution vector.
!
!    real y0(m): the vector of initial conditions at the starting time. 
!
!    real tau0: the initial timestep.
!
!    real reltol, abstol: the tolerances for the local truncation error.
!
!    integer nmax: the maximum number of time steps allowed.
!
!  Output:
!
!    integer nstep: the number of time steps taken.
!
!    integer n_rejected: the number of rejected time steps.
!
!    integer n_fsolve: the number of failed calls to fsolve().
!
!    real t(0:nmax), y(0:nmax,1:m): contains in entries 0:n the 
!    the sequence of times and solution estimates actually computed.
!
!  Local:
!
!    integer count: the number of steps, whether accepted or rejected.
!    It should be the case that n <= count.
!
!    real kappa: a factor for the stepsize update.
!    0 < kappa <= 1.  The value kappa = 0.85 is suggested.
!
!    real theta: the theta-method parameter.
!    theta = 0.0 for the backward Euler method.
!    theta = 0.5 for the midpoint method.
!    theta = 1.0 for the forward Euler method.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer m
  integer nmax

  real ( kind = rk8 ) abstol
  integer count
  logical, parameter :: debug = .true.
  external f
  real ( kind = rk8 ) factor
  logical fsolve_fail
  integer iter
  real ( kind = rk8 ) kappa
  integer n_fsolve
  integer n_rejected
  integer nstep
  real ( kind = rk8 ) reltol
  real ( kind = rk8 ) t(0:nmax)
  real ( kind = rk8 ) t0
  real ( kind = rk8 ) tau(0:nmax)
  real ( kind = rk8 ) tau0
  real ( kind = rk8 ) theta
  real ( kind = rk8 ) tmax
  real ( kind = rk8 ) tn
  real ( kind = rk8 ) tnmax
  real ( kind = rk8 ) tnmin
  real ( kind = rk8 ) y(0:nmax,m)
  real ( kind = rk8 ) y0(m)
  real ( kind = rk8 ) yn(m)

  kappa = 0.85D+00
  theta = 0.5D+00

  n_fsolve = 0
  n_rejected = 0
!
!  On starting this loop with a given value of n,
!  we have t(0:n-1), y(0:n-1,1:m) and tau(0:n-1).
!  We want to take step n, and compute t(n), y(n,1:m), tau(n).
!
  do iter = 0, nmax
!
!  Idiot fact.  If we "do n = 0, nmax", then on typical finish, n is nmax+1!
!
    nstep = iter
!
!  n = 0: set n data
!
    if ( nstep == 0 ) then

      t(nstep) = t0
      y(nstep,1:m) = y0(1:m)
      tau(nstep) = tau0
      count = 0

      if ( debug ) then
        write ( *, '(a,i8)' ) '  nstep = ', nstep
      end if
!
!  Exit if we reached or surpassed final time b.
!
    else if ( tmax <= t(nstep-1) ) then

      exit
!
!  n = 1, 2: y1, y2
!  Steps 1 and 2 use fixed tau.
!
    else if ( nstep <= 2 ) then

      if ( debug ) then
        write ( *, '(a,i8)' ) '  nstep = ', nstep
      end if

      call mad_step ( m, t(nstep-1), y(nstep-1,1:m), tau(nstep-1), f, &
        theta, reltol, t(nstep), y(nstep,1:m), fsolve_fail )

      if ( fsolve_fail ) then
        n_fsolve = n_fsolve + 1
        write ( *, '(a)' ) ''
        write ( *, '(a)' ) 'mad_compute(): Fatal error!'
        write ( *, '(a)' ) '  fsolve() failure on step 1 or 2.'
        stop ( 1 )
      end if

      tau(nstep) = tau0;
      count = count + 1;
!
!  n = 3, ..., nmax: use data from step n-1 to compute data for step n.
!  Tau gets updated.
!  Trying a step to t(n) = t(n-1) + theta * tau(n-1).
!  If step fails, may need to reduce theta and try again.
!
    else

      if ( debug ) then
        write ( *, '(a,i8)' ) '  nstep = ', nstep
      end if

      do while ( .true. )

        if ( debug ) then
          write ( *, '(a,g14.6)' ) '  tau(nstep-1) = ', tau(nstep-1)
        end if

        call mad_step ( m, t(nstep-1), y(nstep-1,1:m), tau(nstep-1), f, &
          theta, reltol, tn, yn, fsolve_fail )

        if ( fsolve_fail ) then
          n_fsolve = n_fsolve + 1
          write ( *, '(a)' ) ''
          write ( *, '(a)' ) 'mad_compute(): Fatal error!'
          write ( *, '(a)' ) '  fsolve() failure.'
          stop ( 1 )
        end if

        t(nstep) = tn
        y(nstep,1:m) = yn(1:m)
        count = count + 1
        if ( debug ) then
          write ( *, '(a,i8)' ) '  count = ', count
        end if
!
!  Estimate local truncation error.
! 
        call mad_lte ( m, y(nstep-3:nstep,1:m), tau(nstep-3:nstep-1), reltol, &
          abstol, tnmin, tnmax )
!
!  If the LTE test was met, we accept T and Y, and set the next time step.
!
        if ( tnmin <= 1.0 ) then
          factor = kappa * ( 1.0 / tnmax ) ** (1.0D+00/3.0D+00)
          factor = min ( factor, 1.5D+00 )
          factor = max ( factor, 0.02D+00 )
          tau(nstep) = factor * tau(nstep-1)
          exit
        else
          n_rejected = n_rejected + 1
        end if

      end do

    end if

  end do

  return
end
subroutine mad_lte ( m, y, tau, reltol, abstol, tnmin, tnmax )

!*****************************************************************************80
!
!! mad_lte() estimates local truncation error for midpoint adaptive method.
!
!  Discussion:
!
!    This function is called by mad_compute() to estimate the local truncation
!    error incurred during a single step of the implicit midpoint method.
!
!    We assume the step has been taken over an interval which we will 
!    represent as going from (t1,y1) to (t2,y2).  We may regard y2 as the
!    estimated solution of the local ODE:
!      dy/dt = f(t,y), y(t1) = y1.
!
!    We suppose that a function y*() is the exact solution of this system
!    and we define the local truncation error for this step as
!      lte = y*(t2) - y2
!    We wish to accept this step as long as lte is small, that is:
!      lte = y*(t2) - y2 < abstol + reltol * max ( y1, y2 )
!
!    Of course, to carry out this check, we need to estimate lte or y*(t2).
!    This is done using the Milne Device, described in the Milne reference.
!    Essentially, lte is estimated using a weighted difference of solution
!    estimates y2 and an estimate formed by an Adams-Bashforth AB2-like 
!    method.  See page 5 of the Burkardt/Trenchea reference for details.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 June 2024
!
!  Author:
!
!    John Burkardt, Catalin Trenchea
!
!  Reference:
!
!    John Burkardt, Catalin Trenchea,
!    Refactorization of the midpoint rule,
!    Applied Mathematics Letters,
!    Volume 107, 106438, 2020.
!
!    William Milne,
!    Numerical Integration of Ordinary Differential Equations,
!    American Mathematical Monthly,
!    Volume 33, number 9, pages 455–460, 1926.
!
!  Input:
!
!    integer m: the spatial dimension of the solution.
!
!    real y(4,m): the last four solution vectors.
!
!    real tau(3): the last three stepsizes.
!
!    real reltol, abstol: local truncation error tolerances.
!
!  Output:
!
!    real tnmin, tnmax: the minimum and maximum entries of the
!    local truncation error estimator.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer m

  real ( kind = rk8 ) abstol
  real ( kind = rk8 ) c1
  real ( kind = rk8 ) c2
  real ( kind = rk8 ) c3
  real ( kind = rk8 ) r
  real ( kind = rk8 ) reltol
  real ( kind = rk8 ) tau(3)
  real ( kind = rk8 ) tn(m)
  real ( kind = rk8 ) tnmax
  real ( kind = rk8 ) tnmin
  real ( kind = rk8 ) uab2(m)
  real ( kind = rk8 ) y(4,m)
!
!  Evaluate the AB2-like solution.
!
  c1 = ( tau(3) + tau(2) ) * ( tau(3) + tau(2) + tau(1) ) &
    / ( tau(2) * ( tau(2) + tau(1) ) )

  c2 = - tau(3) * ( tau(3) + tau(2) + tau(1) ) / ( tau(2) * tau(1) )

  c3 = tau(3) * ( tau(3) + tau(2) ) / ( tau(1) * ( tau(2) + tau(1) ) )

  uab2(1:m) = c1 * y(3,1:m) + c2 * y(2,1:m) + c3 * y(1,1:m)
!
!  Evaluate R, the variable error coefficient in the LTE.
!
  r = 1.0 / 24.0 + 1.0 / 8.0 * ( 1.0 + tau(2) / tau(3) ) &
    * ( 1.0 + 2.0 * tau(2) / tau(3) + tau(1) / tau(3) )
!
!  Use AB2 estimator to approximate the LTE vector.
!
  tn(1:m) = ( y(4,1:m) - uab2(1:m) ) * 24.0 * r / ( 24.0 * r - 1.0 ) &
    / ( abstol + abs ( y(4,1:m) ) * reltol )
   
  tnmax = maxval ( abs ( tn ) )
  tnmin = minval ( abs ( tn ) )

  return
end
subroutine mad_phase_plot ( nstep, m, t, y, label )

!*****************************************************************************80
!
!! mad_phase_plot() writes data and command files for a phase plot.
!
!  Discussion:
!
!    We assume a simple phase plane plot of y(1) versus y(2).
!
!    A better choice would plot any y(i) versus y'(i).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 December 2023
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer nstep: the number of steps.
!
!    integer m: the size of a solution vector.
!
!    real t(0:nstep): the sequence of times.
!
!    real y(0:nstep,1:m): the sequence of solutions.
!
!    character * ( * ) label: a title.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer m
  integer nstep

  character * ( 255 ) command_filename
  integer command_unit
  character * ( 255 ) data_filename
  integer data_unit
  integer i
  integer j
  character * ( * ) label
  character * ( 255 ) png_filename
  real ( kind = rk8 ) t(0:nstep)
  real ( kind = rk8 ) y(0:nstep,1:m)

  data_filename = trim ( label ) //'_phase_data.txt'
  call get_unit ( data_unit )
  open ( unit = data_unit, file = data_filename, status = 'replace' )
  do i = 0, nstep - 1
    write ( data_unit, '(2x,g14.6)', advance = 'no' ) t(i)
    do j = 1, m
      write ( data_unit, '(2x,g14.6)', advance = 'no' ) y(i,j)
    end do
    write ( data_unit, * )
  end do
  close ( unit = data_unit )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Created graphics data file "' &
    // trim ( data_filename ) // '".'

  png_filename = trim ( label ) // '_phase.png'

  command_filename = trim ( label ) // '_phase_commands.txt'
  call get_unit ( command_unit )
  open ( unit = command_unit, file = command_filename, status = 'replace' )
  write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) 'set output "' // trim ( png_filename ) // '"'
  write ( command_unit, '(a)' ) 'set xlabel "Y(1)"'
  write ( command_unit, '(a)' ) 'set ylabel "Y(2)"'
  write ( command_unit, '(a)' ) 'set title "MAD phase plane"'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set style data lines'
  write ( command_unit, '(a,i2,a)' ) 'plot "' // trim ( data_filename ) // &
    '" using 2:3 lw 3 linecolor rgb "blue"'
  write ( command_unit, '(a)' ) 'quit'
  close ( unit = command_unit )
  write ( *, '(a)' ) &
    '  Created command file "' // trim ( command_filename ) // '".'

  return
end
subroutine mad_solution_plot ( nstep, m, t, y, label )

!*****************************************************************************80
!
!! mad_solution_plot() writes data and command files for a solution plot.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 December 2023
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer nstep: the number of steps.
!
!    integer m: the size of a solution vector.
!
!    real t(0:nstep): the sequence of time values.
!
!    real y(0:nstep,1:m): the sequence of solutions.
!
!    character * ( * ) label: a title.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer m
  integer nstep

  character * ( 255 ) command_filename
  integer command_unit
  character * ( 255 ) data_filename
  integer data_unit
  integer i
  integer j
  character * ( * ) label
  character * ( 255 ) png_filename
  real ( kind = rk8 ) t(0:nstep)
  real ( kind = rk8 ) y(0:nstep,1:m)

  data_filename = trim ( label ) //'_solution_data.txt'
  call get_unit ( data_unit )
  open ( unit = data_unit, file = data_filename, status = 'replace' )
  do i = 0, nstep - 1
    write ( data_unit, '(2x,g14.6)', advance = 'no' ) t(i)
    do j = 1, m
      write ( data_unit, '(2x,g14.6)', advance = 'no' ) y(i,j)
    end do
    write ( data_unit, * )
  end do
  close ( unit = data_unit )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Created graphics data file "' &
    // trim ( data_filename ) // '".'

  png_filename = trim ( label ) // '_solution.png'

  command_filename = trim ( label ) // '_solution_commands.txt'
  call get_unit ( command_unit )
  open ( unit = command_unit, file = command_filename, status = 'replace' )
  write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) 'set output "' // trim ( png_filename ) // '"'
  write ( command_unit, '(a)' ) 'set xlabel "Time"'
  write ( command_unit, '(a)' ) 'set ylabel "Y(T)"'
  write ( command_unit, '(a)' ) 'set title "MAD solution Y(T)"'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set style data lines'

  write ( command_unit, '(a)'   ) 'plot \'
  do i = 2, m
    write ( command_unit, '(a,i2,a)' ) &
      '  "' // trim ( data_filename ) // '" using 1:', i, ' lw 3,\'
  end do
  write ( command_unit, '(a,i2,a)' ) &
    '  "' // trim ( data_filename ) // '" using 1:', m + 1, ' lw 3'
  write ( command_unit, '(a)' ) 'quit'
  close ( unit = command_unit )
  write ( *, '(a)' ) &
    '  Created command file "' // trim ( command_filename ) // '".'

  return
end
subroutine mad_stats ( nstep, n_rejected, n_fsolve, t )

!*****************************************************************************80
!
!! mad_stats() reports time and timestep statistics for mad_compute().
!
!  Discussion:
!
!    This function prints some information about the time nodes and time steps.
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
!    John Burkardt, Catalin Trenchea
!
!  Reference:
!
!    Catalin Trenchea, John Burkardt,
!    Refactorization of the midpoint rule,
!    Applied Mathematics Letters,
!    Volume 107, September 2020.
!
!  Input:
!
!    integer nstep: the number of time steps.
!
!    integer n_rejected: the number of rejected time steps.
!
!    integer n_fsolve: the number of failed calls to fsolve().
!
!    real t(0:nstep): the sequence of time values.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer nstep

  integer n_fsolve
  integer n_rejected
  real ( kind = rk8 ) t(0:nstep)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'mad_stats():'
  write ( *, '(a,i6)' ) '  Number of timesteps = ', nstep
  write ( *, '(a,i6)' ) '  Rejected timesteps  = ', n_rejected
  write ( *, '(a,i6)' ) '  fsolve() failures   = ', n_fsolve
  write ( *, '(a,g14.6,a,g14.6,a)' ) &
    '  First time interval is [', t(0), ',', t(1), ']'
  write ( *, '(a,g14.6,a,g14.6,a)' ) &
    '  Final time interval is [', t(nstep-2), ',', t(nstep-1), ']'
  write ( *, '(a,g14.6)' ) '  Minimum timestep is ', minval ( t(1:nstep-1) - t(0:nstep-2) )
  write ( *, '(a,g14.6)' ) '  Maximum timestep is ', maxval ( t(1:nstep-1) - t(0:nstep-2) )

  return
end
subroutine mad_step ( m, to, yo, tauo, f, theta, tol, tn, yn, fsolve_fail )

!*****************************************************************************80
!
!! mad_step() extends the (t,y) solution sequence by one more step.
!
!  Discussion:
!
!    This code is requested to add one more pair of (t,y) data to a string
!    of estimated solutions to the differential equation.
!
!    On input, an initial condition and n solutions have already been 
!    computed and stored in the arrays t and y.  
!
!    A stepsize tau(n) and factor theta are provided.  These define a solution
!    ym of the implicit backward Euler method at tm = t(n) + theta * tau(n). 
!    Then ym is used to construct the implicit midpoint solution y(n+1,1:n) 
!    at time t(n+1) = t(n) + tau(n).
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
!    John Burkardt, Catalin Trenchea
!
!  Reference:
!
!    Catalin Trenchea, John Burkardt,
!    Refactorization of the midpoint rule,
!    Applied Mathematics Letters,
!    Volume 107, September 2020.
!
!  Input:
!
!    integer m: the size of a single solution vector y.
!
!    real to, yo(1:m): the "old" time and solution.
!
!    real tau(0:nmax): the sequence of time steps, indexed from 1 to n.  
!    The entry tau(n) is tentative.
!
!    external f: the function that evaluates the right hand side 
!    of the differential equation, of the form 
!      subroutine f ( t, y, dydt )
!
!    real theta: the theta-method parameter.  Common values include:
!    theta = 0.0 for the backward Euler method.
!    theta = 0.5 for the midpoint method.
!    theta = 1.0 for the forward Euler method.
!
!  Output:
!
!    real tn, yn(1:m): the new time and solution information.
!
!    logical fsolve_fail: true if the call to fsolve() failed.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer m

  external f
  logical fsolve_fail
  integer info
  real ( kind = rk8 ) fm(m)
  real ( kind = rk8 ) tauo
  real ( kind = rk8 ) theta
  real ( kind = rk8 ) tm
  real ( kind = rk8 ) tn
  real ( kind = rk8 ) to
  real ( kind = rk8 ) tol
  real ( kind = rk8 ) ym(m)
  real ( kind = rk8 ) yn(m)
  real ( kind = rk8 ) yo(m)

  fsolve_fail = .false.
!
!  Make a forward solution prediction at to + theta * tau
!
  tm = to + theta * tauo
  call f ( tm, yo(1:m), fm(1:m) )
  ym(1:m) = yo(1:m) + theta * tauo * fm(1:m)
!
!  Starting from the forward prediction, 
!  solve the implicit backward Euler equation for ym.
!
  call fsolve_be ( f, m, to, yo, tm, ym, fm, tol, info )

  if ( info /= 1 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'mad_step(): Fatal error!'
    write ( *, '(a,i2)' ) '  fsolve_be() failed with info = ', info
    fsolve_fail = .true.
    return
  end if
!
!  Now predict the solution at tn.
!
  tn = to + tauo
  yn(1:m) = (       1.0 / theta ) * ym(1:m) &
          + ( 1.0 - 1.0 / theta ) * yo(1:m)

  return
end
subroutine mad_timestep_plot ( nstep, t, label )

!*****************************************************************************80
!
!! mad_timestep_plot() writes data and command files for a timestep plot.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 December 2023
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer nstep: the number of steps.
!
!    real t(0:nstep): the sequence of time values.
!
!    character * ( * ) label: a title.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer nstep

  character * ( 255 ) command_filename
  integer command_unit
  character * ( 255 ) data_filename
  integer data_unit
  integer i
  character * ( * ) label
  character * ( 255 ) png_filename
  real ( kind = rk8 ) t(0:nstep)

  data_filename = trim ( label ) //'_timestep_data.txt'
  call get_unit ( data_unit )
  open ( unit = data_unit, file = data_filename, status = 'replace' )
  do i = 1, nstep - 1
    write ( data_unit, '(2x,i6,2x,g14.6)' ) i, t(i) - t(i-1)
  end do
  close ( unit = data_unit )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Created graphics data file "' &
    // trim ( data_filename ) // '".'

  png_filename = trim ( label ) // '_timestep.png'

  command_filename = trim ( label ) // '_timestep_commands.txt'
  call get_unit ( command_unit )
  open ( unit = command_unit, file = command_filename, status = 'replace' )
  write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) 'set output "' // trim ( png_filename ) // '"'
  write ( command_unit, '(a)' ) 'set xlabel "Index"'
  write ( command_unit, '(a)' ) 'set ylabel "DT"'
  write ( command_unit, '(a)' ) 'set title "MAD timesteps"'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set style data lines'
  write ( command_unit, '(a)' ) 'plot "' // trim ( data_filename ) // &
    '" using 1:2 lw 3 linecolor rgb "blue"'
  write ( command_unit, '(a)' ) 'quit'
  close ( unit = command_unit )
  write ( *, '(a)' ) &
    '  Created command file "' // trim ( command_filename ) // '".'

  return
end

