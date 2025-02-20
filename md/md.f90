subroutine md ( nd, np, step_num, dt )

!*****************************************************************************80
!
!! md() is a molecular dynamics probram.
!
!  Discussion:
!
!    MD implements a simple molecular dynamics simulation.
!
!    The velocity Verlet time integration scheme is used. 
!
!    The particles interact with a central pair potential.
!
!    Based on a FORTRAN90 program by Bill Magro.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 December 2014
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer nd, the spatial dimension (2 or 3);
!
!    integer np, the number of particles (500, for instance);
!
!    integer step_num, the number of time steps (500, for instance).
!
!    real ( kind = rk ) dt, the time step (0.1 for instance )
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: acc(:,:)
  real ( kind = rk ) dt
  real ( kind = rk ) e0
  real ( kind = rk ), allocatable :: force(:,:)
  real ( kind = rk ) kinetic
  real ( kind = rk ), parameter :: mass = 1.0D+00
  integer nd
  integer np
  real ( kind = rk ), allocatable :: pos(:,:)
  real ( kind = rk ) potential
  real ( kind = rk ) rel
  integer seed
  integer step
  integer step_num
  integer step_print
  integer step_print_index
  integer step_print_num
  real ( kind = rk ), allocatable :: vel(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MD'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  A molecular dynamics program.'
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  ND, the spatial dimension, is ', nd
  write ( *, '(a,i8)' ) &
    '  NP, the number of particles in the simulation is ', np
  write ( *, '(a,i8)' ) '  STEP_NUM, the number of time steps, is ', step_num
  write ( *, '(a,g14.6)' ) '  DT, the size of each time step, is ', dt
!
!  Allocate memory.
!
  allocate ( acc(nd,np) )
  allocate ( force(nd,np) )
  allocate ( pos(nd,np) )
  allocate ( vel(nd,np) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  At each step, we report the potential and kinetic energies.'
  write ( *, '(a)' ) '  The sum of these energies should be a constant.'
  write ( *, '(a)' ) '  As an accuracy check, we also print the relative error'
  write ( *, '(a)' ) '  in the total energy.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '      Step      Potential       Kinetic        (P+K-E0)/E0'
  write ( *, '(a)' ) &
    '                Energy P        Energy K       Relative Energy Error'
  write ( *, '(a)' ) ' '
!
!  This is the main time stepping loop:
!    Initialize or update positions, velocities, accelerations.
!    Compute forces and energies,
!
  step_print = 0
  step_print_index = 0
  step_print_num = 10

  do step = 0, step_num

    if ( step == 0 ) then
      seed = 123456789
      call r8mat_uniform_ab ( nd, np, 0.0D+00, 10.0D+00, seed, pos )
      vel(1:nd,1:np) = 0.0D+00
      acc(1:nd,1:np) = 0.0D+00
    else
      call update ( np, nd, pos, vel, force, acc, mass, dt )
    end if

    call compute ( np, nd, pos, vel, mass, force, potential, kinetic )

    if ( step == 0 ) then
      e0 = potential + kinetic
    end if

    if ( step == step_print ) then
      rel = ( potential + kinetic - e0 ) / e0
      write ( *, '(2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
        step, potential, kinetic, rel 
      step_print_index = step_print_index + 1
      step_print = ( step_print_index * step_num ) / step_print_num
    end if

  end do
!
!  Free memory.
!
  deallocate ( acc )
  deallocate ( force )
  deallocate ( pos )
  deallocate ( vel )

  return
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
!      v(x) = ( sin ( min ( x, PI/2 ) ) )^2
!
!    The derivative of the potential is:
!
!      dv(x) = 2.0 * sin ( min ( x, PI/2 ) ) * cos ( min ( x, PI/2 ) )
!            = sin ( 2.0 * min ( x, PI/2 ) )
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
!    John Burkardt
!
!  Input:
!
!    integer NP, the number of particles.
!
!    integer ND, the number of spatial dimensions.
!
!    real ( kind = rk ) POS(ND,NP), the positions.
!
!    real ( kind = rk ) VEL(ND,NP), the velocities.
!
!    real ( kind = rk ) MASS, the mass.
!
!  Output:
!
!    real ( kind = rk ) F(ND,NP), the forces.
!
!    real ( kind = rk ) POT, the total potential energy.
!
!    real ( kind = rk ) KIN, the total kinetic energy.
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

  do i = 1, np
!
!  Compute the potential energy and forces.
!
    f(1:nd,i) = 0.0D+00

    do j = 1, np

      if ( i /= j ) then

        rij(1:nd) = pos(1:nd,i) - pos(1:nd,j)

        d = sqrt ( sum ( rij(1:nd)**2 ) )
!
!  Truncate the distance.
!
        d2 = min ( d, PI2 )
!
!  Attribute half of the total potential energy to particle J.
!
        pot = pot + 0.5D+00 * sin ( d2 ) * sin ( d2 )
!
!  Add particle J's contribution to the force on particle I.
!
        f(1:nd,i) = f(1:nd,i) - rij(1:nd) * sin ( 2.0D+00 * d2 ) / d

      end if

    end do

  end do
!
!  Compute the total kinetic energy.
!
  kin = 0.5D+00 * mass * sum ( vel(1:nd,1:np)**2 )
  
  return
end
subroutine r8mat_uniform_ab ( m, n, a, b, seed, r )

!*****************************************************************************80
!
!! R8MAT_UNIFORM_AB returns a scaled pseudorandom R8MAT.
!
!  Discussion:
!
!    A <= R(I,J) <= B.
!
!    An R8MAT is an array of R8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Input:
!
!    integer M, N, the number of rows and columns
!    in the array.
!
!    real ( kind = rk ) A, B, the lower and upper limits.
!
!    integer SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!  Output:
!
!    real ( kind = rk ) R(M,N), the array of pseudorandom values.
!
!    integer SEED, a seed for the random number generator.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a
  real ( kind = rk ) b
  integer i
  integer, parameter :: i4_huge = 2147483647
  integer j
  integer k
  integer seed
  real ( kind = rk ) r(m,n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_UNIFORM_AB - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  do j = 1, n

    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge
      end if

      r(i,j) = a + ( b - a ) * real ( seed, kind = rk ) * 4.656612875D-10

    end do
  end do

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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
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
!    John Burkardt
!
!  Input:
!
!    integer NP, the number of particles.
!
!    integer ND, the number of spatial dimensions.
!
!    real ( kind = rk ) POS(ND,NP), the positions.
!
!    real ( kind = rk ) VEL(ND,NP), the velocities.
!
!    real ( kind = rk ) F(ND,NP), the forces.
!
!    real ( kind = rk ) ACC(ND,NP), the accelerations.
!
!    real ( kind = rk ) MASS, the mass of each particle.
!
!    real ( kind = rk ) DT, the time step.
!
!  Output:
!
!    real ( kind = rk ) POS(ND,NP), the positions.
!
!    real ( kind = rk ) VEL(ND,NP), the velocities.
!
!    real ( kind = rk ) ACC(ND,NP), the accelerations.
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
  real ( kind = rk ) vel(nd,np)

  do j = 1, np
    do i = 1, nd
      pos(i,j) = pos(i,j) + vel(i,j) * dt + 0.5D+00 * acc(i,j) * dt * dt
      vel(i,j) = vel(i,j) + 0.5D+00 * dt * ( f(i,j) / mass + acc(i,j) )
      acc(i,j) = f(i,j) / mass
    end do
  end do

  return
end
