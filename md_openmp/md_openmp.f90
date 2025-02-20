program main

!*****************************************************************************80
!
!! md_openmp() simulates molecular dynamics, using OpenMP for parallel execution.
!
!  Discussion:
!
!    The velocity Verlet time integration scheme is used. 
!
!    The particles interact with a central pair potential.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 July 2009
!
!  Author:
!
!    Original FORTRAN90 version by Bill Magro.
!    This FORTRAN90 version by John Burkardt.
!
  use omp_lib

  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nd = 3
  integer, parameter :: np = 500

  real ( kind = rk ) acc(nd,np)
  real ( kind = rk ) box(nd)
  real ( kind = rk ), parameter :: dt = 0.01D+00
  real ( kind = rk ) e0
  real ( kind = rk ) force(nd,np)
  real ( kind = rk ) kinetic
  real ( kind = rk ), parameter :: mass = 1.0D+00
  real ( kind = rk ) pos(nd,np)
  real ( kind = rk ) potential
  integer proc_num
  integer step
  integer, parameter :: step_num = 5000
  integer step_print
  integer step_print_index
  integer step_print_num
  integer thread_num
  real ( kind = rk ) vel(nd,np)
  real ( kind = rk ) wtime

  call timestamp ( )

  proc_num = omp_get_num_procs ( )
  thread_num = omp_get_max_threads ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MD_OPENMP'
  write ( *, '(a)' ) '  FORTRAN90/OpenMP version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A molecular dynamics program.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) &
    '  NP, the number of particles in the simulation is ', np
  write ( *, '(a,i8)' ) '  STEP_NUM, the number of time steps, is ', step_num
  write ( *, '(a,g14.6)' ) '  DT, the size of each time step, is ', dt
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of processors available is: ', proc_num
  write ( *, '(a,i8)' ) '  The number of threads available is:    ', thread_num
!
!  Set the dimensions of the box.
!
  box(1:nd) = 10.0D+00
!
!  Set initial positions, velocities, and accelerations.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Initializing positions, velocities, and accelerations.'

  call initialize ( np, nd, box, pos, vel, acc )
!
!  Compute the forces and energies.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computing initial forces and energies.'

  call compute ( np, nd, pos, vel, mass, force, potential, kinetic )
!
!  Save the initial total energy for use in the accuracy check.
!
  e0 = potential + kinetic
!
!  This is the main time stepping loop:
!    Compute forces and energies,
!    Update positions, velocities, accelerations.
!
  step_print = 0
  step_print_index = 0
  step_print_num = 10
  
  step = 0
  write ( *, '(2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    step, potential, kinetic, ( potential + kinetic - e0 ) / e0
  step_print_index = step_print_index + 1
  step_print = ( step_print_index * step_num ) / step_print_num

  wtime = omp_get_wtime ( )

  do step = 1, step_num

    call compute ( np, nd, pos, vel, mass, force, potential, kinetic )

    if ( step == step_print ) then

      write ( *, '(2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
        step, potential, kinetic, ( potential + kinetic - e0 ) / e0

      step_print_index = step_print_index + 1
      step_print = ( step_print_index * step_num ) / step_print_num

    end if

    call update ( np, nd, pos, vel, force, acc, mass, dt )

  end do

  wtime = omp_get_wtime ( ) - wtime
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Elapsed time for main computation:'
  write ( *, '(2x,g14.6,a)' ) wtime, ' seconds'
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MD_OPENMP'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine compute ( np, nd, pos, vel, mass, f, pot, kin )

!*****************************************************************************80
!
!! COMPUTE computes the forces and energies.
!
!  Discussion:
!
!    The computation of forces and energies is fully parallel.
!
!    The potential function V(X) is a harmonic well which smoothly
!    saturates to a maximum value at PI/2:
!
!      v(x) = ( sin ( min ( x, PI2 ) ) )^2
!
!    The derivative of the potential is:
!
!      dv(x) = 2.0D+00 * sin ( min ( x, PI2 ) ) * cos ( min ( x, PI2 ) )
!            = sin ( 2.0 * min ( x, PI2 ) )
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 July 2008
!
!  Author:
!
!    Original FORTRAN90 version by Bill Magro.
!    This FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer NP, the number of particles.
!
!    Input, integer ND, the number of spatial dimensions.
!
!    Input, real ( kind = rk ) POS(ND,NP), the position of each particle.
!
!    Input, real ( kind = rk ) VEL(ND,NP), the velocity of each particle.
!
!    Input, real ( kind = rk ) MASS, the mass of each particle.
!
!    Output, real ( kind = rk ) F(ND,NP), the forces.
!
!    Output, real ( kind = rk ) POT, the total potential energy.
!
!    Output, real ( kind = rk ) KIN, the total kinetic energy.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer np
  integer nd

  real ( kind = rk ) d
  real ( kind = rk ) d2
  real ( kind = rk ) f(nd,np)
  integer i
  integer j
  real ( kind = rk ) kin
  real ( kind = rk ) mass
  real ( kind = rk ), parameter :: PI2 = 3.141592653589793D+00 / 2.0D+00
  real ( kind = rk ) pos(nd,np)
  real ( kind = rk ) pot
  real ( kind = rk ) rij(nd)
  real ( kind = rk ) vel(nd,np)

  pot = 0.0D+00
  kin = 0.0D+00

!$omp parallel &
!$omp shared ( f, nd, np, pos, vel ) &
!$omp private ( d, d2, i, j, rij )

!$omp do reduction ( + : pot, kin )

  do i = 1, np
!
!  Compute the potential energy and forces.
!
    f(1:nd,i) = 0.0D+00

    do j = 1, np

      if ( i /= j ) then

        call dist ( nd, pos(1,i), pos(1,j), rij, d )
!
!  Attribute half of the potential energy to particle J.
!
        d2 = min ( d, PI2 )

        pot = pot + 0.5D+00 * ( sin ( d2 ) )**2

        f(1:nd,i) = f(1:nd,i) - rij(1:nd) * sin ( 2.0D+00 * d2 ) / d

      end if

    end do
!
!  Compute the kinetic energy.
!
    kin = kin + sum ( vel(1:nd,i)**2 )

  end do
!$omp end do

!$omp end parallel

  kin = kin * 0.5D+00 * mass
  
  return
end
subroutine dist ( nd, r1, r2, dr, d )

!*****************************************************************************80
!
!! DIST computes the displacement and distance between two particles.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 March 2002
!
!  Author:
!
!    Original FORTRAN90 version by Bill Magro.
!    This FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ND, the number of spatial dimensions.
!
!    Input, real ( kind = rk ) R1(ND), R2(ND), the positions of the particles.
!
!    Output, real ( kind = rk ) DR(ND), the displacement vector.
!
!    Output, real ( kind = rk ) D, the Euclidean norm of the displacement,
!    in other words, the distance between the two particles.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nd

  real ( kind = rk ) d
  real ( kind = rk ) dr(nd)
  real ( kind = rk ) r1(nd)
  real ( kind = rk ) r2(nd)

  dr(1:nd) = r1(1:nd) - r2(1:nd)

  d = sqrt ( sum ( dr(1:nd)**2 ) )

  return
end
subroutine initialize ( np, nd, box, pos, vel, acc )

!*****************************************************************************80
!
!! INITIALIZE initializes the positions, velocities, and accelerations.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 November 2007
!
!  Author:
!
!    Original FORTRAN90 version by Bill Magro.
!    This FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer NP, the number of particles.
!
!    Input, integer ND, the number of spatial dimensions.
!
!    Input, real ( kind = rk ) BOX(ND), specifies the maximum position
!    of particles in each dimension.
!
!    Output, real ( kind = rk ) POS(ND,NP), the position of each particle.
!
!    Output, real ( kind = rk ) VEL(ND,NP), the velocity of each particle.
!
!    Output, real ( kind = rk ) ACC(ND,NP), the acceleration of each particle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer np
  integer nd

  real ( kind = rk ) acc(nd,np)
  real ( kind = rk ) box(nd)
  integer j
  real ( kind = rk ) pos(nd,np)
  real ( kind = rk ) vel(nd,np)
!
!  Start by setting the positions to random numbers between 0 and 1.
!
  call random_number ( harvest = pos(1:nd,1:np) )
!
!  Use these random values as scale factors to pick random locations
!  inside the box.
!
!$omp parallel &
!$omp shared ( box, nd, np, pos ) &
!$omp private ( j )

!$omp do

  do j = 1, np
    pos(1:nd,j) = box(1:nd) * pos(1:nd,j)
  end do

!$omp end do
!$omp end parallel
!
!  Velocities and accelerations begin at 0.
!
!$omp parallel &
!$omp shared ( acc, nd, np, vel )

!$omp workshare

  vel(1:nd,1:np) = 0.0D+00
  acc(1:nd,1:np) = 0.0D+00

!$omp end workshare

!$omp end parallel

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
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
!    18 May 2013
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
subroutine update ( np, nd, pos, vel, f, acc, mass, dt )

!*****************************************************************************80
!
!! UPDATE updates positions, velocities and accelerations.
!
!  Discussion:
!
!    The time integration is fully parallel.
!
!    A velocity Verlet algorithm is used for the updating.
!
!    x(t+dt) = x(t) + v(t) * dt + 0.5 * a(t) * dt * dt
!    v(t+dt) = v(t) + 0.5 * ( a(t) + a(t+dt) ) * dt
!    a(t+dt) = f(t) / m
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 November 2007
!
!  Author:
!
!    Original FORTRAN90 version by Bill Magro.
!    This FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer NP, the number of particles.
!
!    Input, integer ND, the number of spatial dimensions.
!
!    Input/output, real ( kind = rk ) POS(ND,NP), the position of each particle.
!
!    Input/output, real ( kind = rk ) VEL(ND,NP), the velocity of each particle.
!
!    Input, real ( kind = rk ) F(ND,NP), the force on each particle.
!
!    Input/output, real ( kind = rk ) ACC(ND,NP), the acceleration of each
!    particle.
!
!    Input, real ( kind = rk ) MASS, the mass of each particle.
!
!    Input, real ( kind = rk ) DT, the time step.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer np
  integer nd

  real ( kind = rk ) acc(nd,np)
  real ( kind = rk ) dt
  real ( kind = rk ) f(nd,np)
  integer i
  integer j
  real ( kind = rk ) mass
  real ( kind = rk ) pos(nd,np)
  real ( kind = rk ) rmass
  real ( kind = rk ) vel(nd,np)

  rmass = 1.0D+00 / mass

!$omp parallel &
!$omp shared ( acc, dt, f, nd, np, pos, rmass, vel ) &
!$omp private ( i, j )

!$omp do
  do j = 1, np
    do i = 1, nd
      pos(i,j) = pos(i,j) + vel(i,j) * dt + 0.5D+00 * acc(i,j) * dt * dt
      vel(i,j) = vel(i,j) + 0.5D+00 * dt * ( f(i,j) * rmass + acc(i,j) )
      acc(i,j) = f(i,j) * rmass
    end do
  end do
!$omp end do

!$omp end parallel

  return
end
