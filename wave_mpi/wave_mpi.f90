program main

!*****************************************************************************80
!
!! WAVE_MPI() solves the wave equation in parallel using MPI.
!
!  Discussion:
!
!    Discretize the equation for u(x,t):
!      d^2 u/dt^2 - c^2 * d^2 u/dx^2 = 0  for 0 < x < 1, 0 < t
!    with boundary conditions:
!      u(0,t) = u0(t) = sin ( 2 * pi * ( 0 - c * t ) )
!      u(1,t) = u1(t) = sin ( 2 * pi * ( 1 - c * t ) )
!    and initial conditions:
!         u(x,0) = g(x,t=0) =                sin ( 2 * pi * ( x - c * t ) )
!      dudt(x,0) = h(x,t=0) = - 2 * pi * c * cos ( 2 * pi * ( x - c * t ) ) 
!
!    by:
!
!      alpha = c * dt / dx.
!
!      U(x,t+dt) = 2 U(x,t) - U(x,t-dt) 
!        + alpha^2 ( U(x-dx,t) - 2 U(x,t) + U(x+dx,t) ).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 November 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Geoffrey Fox, Mark Johnson, Gregory Lyzenga, Steve Otto, John Salmon, 
!    David Walker,
!    Solving problems on concurrent processors, 
!    Volume 1: General Techniques and Regular Problems,
!    Prentice Hall, 1988,
!    ISBN: 0-13-8230226,
!    LC: QA76.5.F627.
!
!  Local parameters:
!
!    Local, real ( kind = rk ) DT, the time step.
!
!    Local, integer ID, the MPI process ID.
!
!    Local, integer N_GLOBAL, the total number of points.
!
!    Local, integer N_LOCAL, the number of points visible to 
!    this process.
!
!    Local, integer NSTEPS, the number of time steps.
!
!    Local, integer P, the number of MPI processes.
!
  use mpi

  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), parameter :: dt = 0.00125
  integer error
  integer i_global_hi
  integer i_global_lo
  integer id
  integer, parameter :: n_global = 401
! integer, parameter :: n_global = 101
  integer n_local
  integer, parameter :: nsteps = 4000
! integer, parameter :: nsteps = 2
  integer p
  real ( kind = rk ), allocatable :: u1_local(:)
  real ( kind = rk ) wtime
!
!  Initialize MPI.
!
  call MPI_Init ( error )

  call MPI_Comm_rank ( MPI_COMM_WORLD, id, error )

  call MPI_Comm_size ( MPI_COMM_WORLD, p, error )

  if ( id == 0 ) then
    call timestamp ( )
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'mpi_wave():'
    write ( *, '(a)' ) '  FORTRAN90 version.'
    write ( *, '(a)' ) '  Estimate a solution of the wave equation using MPI.'
    write ( *, '(a)' ) ''
    write ( *, '(a,i2,a)' ) '  Using ', p, ' processes.'
    write ( *, '(a,i6,a)' ) '  Using a total of ', n_global, ' points.'
    write ( *, '(a,i6,a,g14.6)' ) '  Using ', nsteps, ' time steps of size ', dt
    write ( *, '(a,g14.6)' ) '  Computing final solution at time ', dt * nsteps
  end if

  wtime = MPI_Wtime ( )
!
!  Determine N_LOCAL
!
  i_global_lo = (   id       * ( n_global - 1 ) ) / p
  i_global_hi = ( ( id + 1 ) * ( n_global - 1 ) ) / p
  if ( 0 < id ) then
    i_global_lo = i_global_lo - 1
  end if

  n_local = i_global_hi + 1 - i_global_lo
  allocate ( u1_local(1:n_local) )
!
!  Update N_LOCAL values.
!
  call update ( id, p, n_global, n_local, nsteps, dt, u1_local )
!
!  Collect local values into global array.
!
  call collect ( id, p, n_global, n_local, nsteps, dt, u1_local )
!
!  Report elapsed wallclock time.
!
  wtime = MPI_Wtime ( ) - wtime

  if ( id == 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a,g14.6,a)' ) &
      '  Elapsed wallclock time was ', wtime, ' seconds.'
  end if
!
!  Terminate MPI.
!
  call MPI_Finalize ( error )
!
!  Free memory.
!
  deallocate ( u1_local )
!
!  Terminate.
!
  if ( id == 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'wave_mpi():'
    write ( *, '(a)' ) '  Normal end of execution.'
    write ( *, '(a)' ) ''
    call timestamp ( )
  end if

  stop
end
subroutine update ( id, p, n_global, n_local, nsteps, dt, u1_local ) 

!*****************************************************************************80
!
!! UPDATE advances the solution a given number of time steps.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 November 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ID, the identifier of this process.
!
!    Input, integer P, the number of processes.
!
!    Input, integer N_GLOBAL, the total number of points.
!
!    Input, integer N_LOCAL, the number of points visible to 
!    this process.
!
!    Input, integer NSTEPS, the number of time steps.
!
!    Input, real ( kind = rk ) DT, the size of the time step.
!
!    Output, real ( kind = rk ) U1_LOCAL(N_LOCAL), the portion of the solution
!    at the last time, as evaluated by this process.
!
  use mpi

  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n_local

  real ( kind = rk ) alpha
  real ( kind = rk ) c
  real ( kind = rk ) dt
  real ( kind = rk ) dudt
  real ( kind = rk ) dx
  integer error
  real ( kind = rk ) exact
  integer i
  integer i_global
  integer i_global_hi
  integer i_global_lo
  integer i_local
  integer i_local_hi
  integer i_local_lo
  integer id
  integer, parameter :: ltor = 20
  integer n_global
  integer nsteps
  integer p
  integer, parameter :: rtol = 10
  integer status(MPI_STATUS_SIZE)
  real ( kind = rk ) t
  real ( kind = rk ) u0_local(n_local)
  real ( kind = rk ) u1_local(n_local)
  real ( kind = rk ) u2_local(n_local)
  real ( kind = rk ) x
!
!  Determine the value of ALPHA.
!
  c = 1.0D+00
  dx = 1.0D+00 / real ( n_global - 1, kind = rk )
  alpha = c * dt / dx

  if ( 1.0D+00 <= abs ( alpha ) ) then

    if ( id == 0 ) then
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'UPDATE - Warning!'
      write ( *, '(a)' ) '  1 <= |ALPHA| = | C * dT / dX |.'
      write ( *, '(a,g14.6)' ) '  C = ', c
      write ( *, '(a,g14.6)' ) '  dT = ', dt
      write ( *, '(a,g14.6)' ) '  dX = ', dx
      write ( *, '(a,g14.6)' ) '  ALPHA = ', alpha
      write ( *, '(a)' ) '  Computation will not be stable!'
    end if

    call MPI_Finalize ( error )
    stop 1

  end if
!
!  The global array of N_GLOBAL points must be divided up among the processes.
!  Each process stores about 1/P of the total + 2 extra slots.
!
  i_global_lo = (   id       * ( n_global - 1 ) ) / p
  i_global_hi = ( ( id + 1 ) * ( n_global - 1 ) ) / p
  if ( 0 < id ) then
    i_global_lo = i_global_lo - 1
  end if

  i_local_lo = 0
  i_local_hi = i_global_hi - i_global_lo

  t = 0.0D+00
  do i_global = i_global_lo, i_global_hi
    x = real ( i_global, kind = rk ) / real ( n_global - 1, kind = rk )
    i_local = i_global - i_global_lo
    u1_local(i_local+1) = exact ( x, t )
  end do

  do i_local = i_local_lo, i_local_hi
    u0_local(i_local+1) = u1_local(i_local+1)
  end do
!
!  Take NSTEPS time steps.
!
  do i = 1, nsteps

    t = dt * real ( i, kind = rk )
! 
!  For the first time step, we need to use the initial derivative information.
!
    if ( i == 1 ) then

      do i_local = i_local_lo + 1, i_local_hi - 1
        i_global = i_global_lo + i_local
        x = real ( i_global, kind = rk ) / real ( n_global - 1, kind = rk )
        u2_local(i_local+1) = &
          +         0.5D+00 * alpha * alpha   * u1_local(i_local-1+1) &
          + ( 1.0D+00 -       alpha * alpha ) * u1_local(i_local+0+1)   &
          +         0.5D+00 * alpha * alpha   * u1_local(i_local+1+1) &
          +                                dt * dudt ( x, t )
      end do
!
!  After the first time step, we can use the previous two solution estimates.
!
    else

      do i_local = i_local_lo + 1, i_local_hi - 1
        u2_local(i_local+1) = &
          +               alpha * alpha   * u1_local(i_local-1+1) &
          + 2.0 * ( 1.0 - alpha * alpha ) * u1_local(i_local+0+1)   &
          +               alpha * alpha   * u1_local(i_local+1+1) &
          -                                 u0_local(i_local+0+1)
      end do

    end if
!
!  Exchange data with "left-hand" neighbor. 
!
    if ( 0 < id ) then
      call MPI_Send ( u2_local(i_local_lo+2), 1, MPI_DOUBLE_PRECISION, &
        id - 1, rtol, MPI_COMM_WORLD, error )
      call MPI_Recv ( u2_local(i_local_lo+1), 1, MPI_DOUBLE_PRECISION, &
        id - 1, ltor, MPI_COMM_WORLD, status, error )
    else
      x = 0.0D+00
      u2_local(i_local_lo+1) = exact ( x, t )
    end if
!
!  Exchange data with "right-hand" neighbor.
!
    if ( id < p - 1 ) then
      call MPI_Send ( u2_local(i_local_hi), 1, MPI_DOUBLE_PRECISION, id + 1, &
        ltor, MPI_COMM_WORLD, error )
      call MPI_Recv ( u2_local(i_local_hi+1),   1, MPI_DOUBLE_PRECISION, &
        id + 1, rtol, MPI_COMM_WORLD, status, error )
    else
      x = 1.0D+00
      u2_local(i_local_hi+1) = exact ( x, t )
    end if
!
!  Shift data for next time step.
!
    do i_local = i_local_lo, i_local_hi
      u0_local(i_local+1) = u1_local(i_local+1)
      u1_local(i_local+1) = u2_local(i_local+1)
    end do

  end do

  return
end
subroutine collect ( id, p, n_global, n_local, nsteps, dt, u_local ) 

!*****************************************************************************80
!
!! COLLECT has workers send results to the master, which prints them.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 November 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ID, the identifier of this process.
!
!    Input, integer P, the number of processes.
!
!    Input, integer N_GLOBAL, the total number of points.
!
!    Input, integer N_LOCAL, the number of points visible
!    to this process.
!
!    Input, integer NSTEPS, the number of time steps.
!
!    Input, real ( kind = rk ) DT, the size of the time step.
!
!    Input, real ( kind = rk ) U_LOCAL(N_LOCAL), the final solution estimate 
!    computed by this process.
!
  use mpi

  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n_global
  integer n_local

  integer buffer(2)
  integer, parameter :: collect1 = 10
  integer, parameter :: collect2 = 20
  real ( kind = rk ) dt
  integer error
  real ( kind = rk ) exact
  integer i
  integer i_global
  integer i_global_hi
  integer i_global_lo
  integer i_local
  integer i_local_hi
  integer i_local_lo
  integer id
  integer n_local2
  integer nsteps
  integer p
  integer status(MPI_STATUS_SIZE)
  real ( kind = rk ) t
  real ( kind = rk ), allocatable :: u_global(:)
  real ( kind = rk ) u_local(n_local)
  real ( kind = rk ) x

  i_global_lo = (   id       * ( n_global - 1 ) ) / p
  i_global_hi = ( ( id + 1 ) * ( n_global - 1 ) ) / p
  if ( 0 < id ) then
    i_global_lo = i_global_lo - 1
  end if

  i_local_lo = 0
  i_local_hi = i_global_hi - i_global_lo
!
!  Master collects worker results into the U_GLOBAL array.
!
  if ( id == 0 ) then
!
!  Create the global array.
!
    allocate ( u_global(1:n_global) )
!
!  Copy the master's results into the global array.
!
    do i_local = i_local_lo, i_local_hi
      i_global = i_global_lo + i_local - i_local_lo
      u_global(i_global+1) = u_local(i_local+1)
    end do
!
!  Contact each worker.
!
    do i = 1, p - 1
!
!  Message "collect1" contains the global index and number of values.
!
      call MPI_Recv ( buffer, 2, MPI_INTEGER, i, collect1, MPI_COMM_WORLD, &
        status, error )

      i_global_lo = buffer(1)
      n_local2 = buffer(2)

      if ( i_global_lo < 0 ) then
        write ( *, '(a,i6)' ) '  Illegal I_GLOBAL_LO = ', i_global_lo
        call MPI_Finalize ( error )
        stop 1
      else if ( n_global <= i_global_lo + n_local2 - 1 ) then
        write ( *, '(a,i6)' ) '  Illegal I_GLOBAL_LO + N_LOCAL = ', &
          i_global_lo + n_local2
        call MPI_Finalize ( error )
        stop 1
      end if
!
!  Message "collect2" contains the values.
!
      call MPI_Recv ( u_global(i_global_lo+1), n_local2, MPI_DOUBLE_PRECISION, &
        i, collect2, MPI_COMM_WORLD, status, error )

    end do
!
!  Print the results.
!
    t = dt * real ( nsteps, kind = rk )
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) '    I      X     F(X)   Exact'
    write ( *, '(a)' ) ''
    do i_global = 0, n_global - 1
      x = real ( i_global, kind = rk ) / real ( n_global - 1, kind = rk )
      write ( *, '(2x,i3,2x,f6.3,2x,f6.3,2x,f6.3)' ) &
        i_global, x, u_global(i_global+1), exact ( x, t )
    end do

    deallocate ( u_global )
!
!  Workers send results to process 0.
!
  else
!
!  Message "collect1" contains the global index and number of values.
!
    buffer(1) = i_global_lo
    buffer(2) = n_local
    call MPI_Send ( buffer, 2, MPI_INTEGER, 0, collect1, MPI_COMM_WORLD, error )
!
!  Message "collect2" contains the values.
!
    call MPI_Send ( u_local, n_local, MPI_DOUBLE_PRECISION, 0, collect2, &
      MPI_COMM_WORLD, error )

  end if

  return
end
function exact ( x, t )

!*****************************************************************************80
!
!! EXACT evaluates the exact solution
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 November 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, the location.
!
!    Input, real ( kind = rk ) T, the time.
!
!    Output, real ( kind = rk ) EXACT, the value of the exact solution.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), parameter :: c = 1.0D+00
  real ( kind = rk ) exact
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) t
  real ( kind = rk ) x

  exact = sin ( 2.0D+00 * pi * ( x - c * t ) )

  return
end
function dudt ( x, t )

!*****************************************************************************80
!
!! DUDT evaluates the partial derivative dudt.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 November 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, the location.
!
!    Input, real ( kind = rk ) T, the time.
!
!    Output, real ( kind = rk ) DUDT, the value of the time derivative of 
!    the solution.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), parameter :: c = 1.0D+00
  real ( kind = rk ) dudt
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) t
  real ( kind = rk ) x

  dudt = - 2.0D+00 * pi * c * cos ( 2.0D+00 * pi * ( x - c * t ) )

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
