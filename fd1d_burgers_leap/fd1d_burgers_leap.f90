program main

!*****************************************************************************80
!
!! FD1D_BURGERS_LEAP solves the nonviscous Burgers equation using leapfrogging.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) dt
  real ( kind = rk ) dx
  integer ihi
  integer ilo
  integer n
  integer step
  integer step_num
  real ( kind = rk ) t
  real ( kind = rk ) t_init
  real ( kind = rk ) t_last
  real ( kind = rk ), allocatable :: uc(:)
  real ( kind = rk ), allocatable :: un(:)
  real ( kind = rk ), allocatable :: uo(:)
  real ( kind = rk ), allocatable :: x(:)

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FD1D_BURGERS_LEAP:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Solve the non-viscous time-dependent Burgers equation,'
  write ( *, '(a)' ) '  using the leap-frog method.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Equation to be solved:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    du/dt + u * du/dx = 0'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  for x in [ a, b ], for t in [t_init, t_last]'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  with initial conditions:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    u(x,0) = u_init'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  and boundary conditions:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    u(a,t) = 0.5, u(b,t) = -0.5'
!
!  Set and report the problem parameters.
!
  n = 21
  a = -1.0D+00
  b = +1.0D+00
  dx = ( b - a ) / real ( n - 1, kind = rk )
  step_num = 30
  t_init = 0.0D+00
  t_last = 3.0D+00
  dt = ( t_last - t_init ) / real ( step_num, kind = rk )

  write ( *, '(a)' ) ' '
  write ( *, '(2x,g14.6,a,g14.6)' ) a, ' <= X <= ', b
  write ( *, '(a,i8)' ) '  Number of nodes = ', n
  write ( *, '(a,g14.6)' ) '  DX = ', dx
  write ( *, '(a)' ) ' '
  write ( *, '(2x,g14.6,a,g14.6)' ) t_init, ' <= T <= ', t_last
  write ( *, '(a,i8)' ) '  Number of time steps = ', step_num
  write ( *, '(a,g14.6)' ) '  DT = ', dt

  allocate ( uc(1:n) )
  allocate ( un(1:n) )
  allocate ( uo(1:n) )
  allocate ( x(1:n) )

  call r8vec_even ( n, a, b, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X:'
  write ( *, '(a)' ) ' '
  do ilo = 1, n, 5
    ihi = min ( ilo + 4, n )
    write ( *, '(2x,5g14.6)' ) x(ilo:ihi)
  end do
!
!  Set the initial condition,
!  and apply boundary conditions to first and last entries.
!
  step = 0
  t = t_init
  call u_init ( n, x, un )
  un(1) = +0.5D+00
  un(n) = -0.5D+00

  call report ( step, n, t, un )
!
!  Use Euler's method to get the first step.
!
  step = 1
  t = ( real ( step_num - step, kind = rk ) * t_init   &
      + real (            step, kind = rk ) * t_last ) &
      / real ( step_num,        kind = rk )

  uc(1:n) = un(1:n)

  un(2:n-1) = uc(2:n-1) - dt * uc(2:n-1) * ( uc(3:n) - uc(1:n-2) ) &
    / 2.0D+00 / dx

  un(1) = +0.5D+00
  un(n) = -0.5D+00

  call report ( step, n, t, un )
!
!  Subsequent steps use the leapfrog method.
!
  do step = 2, step_num
 
    t = ( real ( step_num - step, kind = rk ) * t_init   &
        + real (            step, kind = rk ) * t_last ) &
        / real ( step_num,        kind = rk )

    uo(1:n) = uc(1:n)
    uc(1:n) = un(1:n)

    un(2:n-1) = uo(2:n-1) - dt * uc(2:n-1) * ( uc(3:n) - uc(1:n-2) ) / dx

    un(1) = +0.5D+00
    un(n) = -0.5D+00

    call report ( step, n, t, un )

  end do
!
!  Free memory.
!
  deallocate ( uc )
  deallocate ( un )
  deallocate ( uo )
  deallocate ( x )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FD1D_BURGERS_LEAP:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine r8vec_even ( n, alo, ahi, a )

!*****************************************************************************80
!
!! R8VEC_EVEN returns an R8VEC of evenly spaced values.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    If N is 1, then the midpoint is returned.
!
!    Otherwise, the two endpoints are returned, and N-2 evenly
!    spaced points between them.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of values.
!
!    Input, real ( kind = rk ) ALO, AHI, the low and high values.
!
!    Output, real ( kind = rk ) A(N), N evenly spaced values.
!    Normally, A(1) = ALO and A(N) = AHI.
!    However, if N = 1, then A(1) = 0.5*(ALO+AHI).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) ahi
  real ( kind = rk ) alo
  integer i

  if ( n == 1 ) then

    a(1) = 0.5D+00 * ( alo + ahi )

  else

    do i = 1, n
      a(i) = ( real ( n - i,     kind = rk ) * alo   &
             + real (     i - 1, kind = rk ) * ahi ) &
             / real ( n     - 1, kind = rk )
    end do

  end if

  return
end
subroutine report ( step, n, t, u )

!*****************************************************************************80
!
!! REPORT prints or plots or saves the data at the current time step.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer STEP, the index of the current step,
!    between 0 and STEP_NUM.
!
!    Input, integer N, the number of nodes.
!
!    Input, real ( kind = rk ) T, the current time.
!
!    Input, real ( kind = rk ) U(N), the initial values U(X,T).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer ihi
  integer ilo
  integer step
  real ( kind = rk ) t
  real ( kind = rk ) u(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  STEP = ', step
  write ( *, '(a,g14.6)' ) '  TIME = ', t
  write ( *, '(a)' ) ' '
  do ilo = 1, n, 5
    ihi = min ( ilo + 4, n )
    write ( *, '(2x,5g14.6)' ) u(ilo:ihi)
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
!  Parameters:
!
!    None
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

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
subroutine u_init ( n, x, u )

!*****************************************************************************80
!
!! U_INIT sets the initial condition for U.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of nodes.
!
!    Input, real ( kind = rk ) X(N), the coordinates of the nodes.
!
!    Output, real ( kind = rk ) U(N), the initial values U(X,T).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) q
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) u(n)
  real ( kind = rk ) ua
  real ( kind = rk ) ub
  real ( kind = rk ) x(n)

  ua = + 0.5D+00
  ub = - 0.5D+00

  q = 2.0D+00 * ( ua - ub ) / pi
  r = ( ua + ub ) / 2.0D+00
!
!  S can be varied.  It is the slope of the initial condition at the midpoint.
!
  s = 1.0D+00

  u(1:n) = - q * atan ( s * ( 2.0D+00 * x(1:n) - x(1) - x(n) ) &
    / ( x(n) - x(1) ) ) + r

  return
end

