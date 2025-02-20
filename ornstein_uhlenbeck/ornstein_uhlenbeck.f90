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
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
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
!  Parameters:
!
!    Output, integer IUNIT, the free unit number.
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
subroutine ou_euler ( theta, mu, sigma, x0, tmax, n )

!*****************************************************************************80
!
!! OU_EULER applies the Euler method to the Ornstein-Uhlenbeck SDE.
!
!  Discussion:
!
!    The stochastic differential equation (SDE) is:
!
!      dx(t) = theta * ( mu - x(t) ) dt + sigma dW,   
!      x(0) = x0.
!
!    The discretized Brownian path uses a constant stepsize.
!
!    For an SDE of the form:
!
!      dx = f(x(t)) dt + g(x(t)) dW(t),
!
!    the Euler method has the form:
!
!      x(j) = x(j-1) + f(x(j-1)) * dt + g(x(j-1)) * dW(j-1)
!
!    Note that if SIGMA is zero, the problem becomes deterministic,
!    with solution:
!
!      x(t) = mu + ( x0 - mu ) * exp ( - theta * t )
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
!  Reference:
!
!    Desmond Higham,
!    An Algorithmic Introduction to Numerical Simulation of
!    Stochastic Differential Equations,
!    SIAM Review,
!    Volume 43, Number 3, September 2001, pages 525-546
!
!  Parameters:
!
!    Input, real ( kind = rk ) THETA, MU, SIGMA, the value of problem parameters.
!
!    Input, real ( kind = rk ) X0, the initial condition.  When studying many
!    realizations of this problem, it is usual for X0 to be chosen
!    from a normal distribution.
!
!    Input, real ( kind = rk ) TMAX, the final time.
!
!    Input, integer N, the number of time steps.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  character ( len = 255 ) command_filename
  integer command_unit
  character ( len = 255 ) data_filename
  integer data_unit
  real ( kind = rk ) dt
  real ( kind = rk ) dw(n)
  integer j
  real ( kind = rk ) mu
  real ( kind = rk ) sigma
  real ( kind = rk ) t(0:n)
  real ( kind = rk ) theta
  real ( kind = rk ) tmax
  real ( kind = rk ) x(0:n)
  real ( kind = rk ) x0

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'OU_EULER:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Use an Euler method to approximate the solution of'
  write ( *, '(a)' ) '  the Ornstein-Uhlenbeck stochastic differential equation:'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '    d x(t) = theta * ( mu - x(t) ) dt + sigma dW'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  with initial condition x(0) = x0.'
!
!  Set the discrete time stepsize.
!
  dt = tmax / real ( n, kind = rk )
!
!  Compute the Brownian increments.
!
  call r8vec_normal_01 ( n, dw )
  dw(1:n) = dw(1:n) * sqrt ( dt )
!
!  Carry out the Euler approximate integration process.
!
  call r8vec_linspace ( n + 1, 0.0D+00, tmax, t )

  x(0) = x0
  do j = 1, n
    x(j) = x(j-1) + dt * theta * ( mu - x(j-1) ) + sigma * dw(j)
  end do
!
!  Create the data file.
!
  call get_unit ( data_unit )
  data_filename = 'ou_euler_data.txt'
  open ( unit = data_unit, file = data_filename, status = 'replace' )
  do j = 0, n
    write ( data_unit, '(2x,g14.6,2x,g14.6)' ) t(j), x(j)
  end do
  close ( unit = data_unit )
  write ( *, '(a)' ) '  Created data file "' // trim ( data_filename ) // '".'
!
!  Create the command file.
!
  call get_unit ( command_unit )
  command_filename = 'ou_euler_commands.txt'
  open ( unit = command_unit, file = command_filename, status = 'replace' )
  write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) &
    'set output "ou_euler.png"'
  write ( command_unit, '(a)' ) 'set xlabel "<--- T --->"'
  write ( command_unit, '(a)' ) 'set ylabel "<--- X(T) --->"'
  write ( command_unit, '(a)' ) &
    'set title "Euler Solution of Ornstein-Uhlenbeck SDE"'
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
subroutine ou_euler_maruyama ( theta, mu, sigma, x0, tmax, n, r )

!*****************************************************************************80
!
!! OU_EULER_MARUYAMA applies Euler-Maruyama to the Ornstein-Uhlenbeck SDE.
!
!  Discussion:
!
!    The stochastic differential equation (SDE) is:
!
!      dx = theta * ( mu - x(t) ) dt + sigma dW,   
!      x(0) = x0,
!
!    The discretized Brownian path uses a constant stepsize.
!
!    A "large" time step DT_LARGE is used for the smooth portion of the
!    equation, while a smaller time step DT_SMALL is used for the
!    discretized Brownian path.  We take R small steps to equal one 
!    large step, so that:
!
!      dt_large = r * dt_small = tmax / n
!
!    For an SDE of the form:
!
!      dx = f(x(t)) dt + g(x(t)) dW(t)
!
!    the Euler-Maruyama method has the form:
!
!      x(j) = x(j-1) + f(X(j-1)) * dt_large + g(X(j-1)) * dW(j-1)
!
!    where dW(j-1) is approximated by the sum of R normal random values
!    multiplied by the square root of DT_SMALL.
!
!    Note that if SIGMA is zero, the problem becomes deterministic,
!    with solution
!
!      x(t) = mu + ( x0 - mu ) * exp ( - theta * t )
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
!  Reference:
!
!    Desmond Higham,
!    An Algorithmic Introduction to Numerical Simulation of
!    Stochastic Differential Equations,
!    SIAM Review,
!    Volume 43, Number 3, September 2001, pages 525-546
!
!  Parameters:
!
!    Input, real ( kind = rk ) THETA, MU, SIGMA, the value of problem parameters.
!
!    Input, real ( kind = rk ) X0, the initial condition.  When studying many
!    realizations of this problem, it is usual for X0 to be chosen
!    from a normal distribution.
!
!    Input, real ( kind = rk ) TMAX, the final time.
!
!    Input, integer N, the number of large scale time steps.
!
!    Input, integer R, the number of small scale time steps per single
!    large scale time step.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer r

  character ( len = 255 ) command_filename
  integer command_unit
  character ( len = 255 ) data_filename
  integer data_unit
  real ( kind = rk ) dt_large
  real ( kind = rk ) dt_small
  real ( kind = rk ) dw(r)
  integer j
  real ( kind = rk ) mu
  real ( kind = rk ) sigma
  real ( kind = rk ) t(0:n)
  real ( kind = rk ) theta
  real ( kind = rk ) tmax
  real ( kind = rk ) x(0:n)
  real ( kind = rk ) x0

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'OU_EULER_MARUYAMA:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Use an Euler-Maruyama method to approximate the solution of'
  write ( *, '(a)' ) '  the Ornstein-Uhlenbeck stochastic differential equation:'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '    d x(t) = theta * ( mu - x(t) ) dt + sigma dW'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  with initial condition x(0) = x0.'
!
!  Set time steps.
!
  dt_large = tmax / real ( n, kind = rk )
  dt_small = tmax / real ( n, kind = rk ) / real ( r, kind = rk )
!
!  Carry out the Euler-Maruyama approximate integration process.
!
  call r8vec_linspace ( n + 1, 0.0D+00, tmax, t )

  x(0) = x0
  do j = 1, n
    call r8vec_normal_01 ( r, dw )
    dw(1:r) = dw(1:r) * sqrt ( dt_small )
    x(j) = x(j-1) + dt_large * theta * ( mu - x(j-1) ) + sigma * sum ( dw(1:r) )
  end do
!
!  Plot the approximate solution.
!
  call get_unit ( data_unit )
  data_filename = 'ou_euler_maruyama_data.txt'
  open ( unit = data_unit, file = data_filename, status = 'replace' )
  do j = 0, n
    write ( data_unit, '(2x,g14.6,2x,g14.6)' ) t(j), x(j)
  end do
  close ( unit = data_unit )
  write ( *, '(a)' ) '  Created data file "' // trim ( data_filename ) // '".'

  call get_unit ( command_unit )
  command_filename = 'ou_euler_maruyama_commands.txt'
  open ( unit = command_unit, file = command_filename, status = 'replace' )
  write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) &
    'set output "ou_euler_maruyama.png"'
  write ( command_unit, '(a)' ) 'set xlabel "<--- T --->"'
  write ( command_unit, '(a)' ) 'set ylabel "<--- X(T) --->"'
  write ( command_unit, '(a)' ) &
    'set title "Euler-Maruyama Solution of Ornstein-Uhlenbeck SDE"'
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
subroutine r8vec_linspace ( n, a, b, x )

!*****************************************************************************80
!
!! R8VEC_LINSPACE creates a vector of linearly spaced values.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    4 points evenly spaced between 0 and 12 will yield 0, 4, 8, 12.
!
!    In other words, the interval is divided into N-1 even subintervals,
!    and the endpoints of intervals are used as the points.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real ( kind = rk ) A_FIRST, A_LAST, the first and last entries.
!
!    Output, real ( kind = rk ) X(N), a vector of linearly spaced data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a
  real ( kind = rk ) b
  integer i
  real ( kind = rk ) x(n)

  if ( n == 1 ) then

    x(1) = ( a + b ) / 2.0D+00

  else

    do i = 1, n
      x(i) = ( real ( n - i,     kind = rk ) * a   &
             + real (     i - 1, kind = rk ) * b ) &
             / real ( n     - 1, kind = rk )
    end do

  end if

  return
end
subroutine r8vec_normal_01 ( n, x )

!*****************************************************************************80
!
!! R8VEC_NORMAL_01 returns a unit pseudonormal R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!    This routine can generate a vector of values on one call.  It
!    has the feature that it should provide the same results
!    in the same order no matter how we break up the task.
!
!    The Box-Muller method is used, which is efficient, but
!    generates an even number of values each time.  On any call
!    to this routine, an even number of new values are generated.
!    Depending on the situation, one value may be left over.
!    In that case, it is saved for the next call.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of values desired.  If N is
!    negative,then the code will flush its internal memory; in particular,
!    if there is a saved value to be used on the next call, it is
!    instead discarded.
!
!    Output, real ( kind = rk ) X(N), a sample of the standard normal PDF.
!
!  Local parameters:
!
!    Local, integer MADE, records the number of values that have
!    been computed.  On input with negative N, this value overwrites
!    the return value of N, so the user can get an accounting of
!    how much work has been done.
!
!    Local, real ( kind = rk ) R(N+1), is used to store some uniform
!    random values.  Its dimension is N+1, but really it is only needed
!    to be the smallest even number greater than or equal to N.
!
!    Local, integer SAVED, is 0 or 1 depending on whether there is a
!    single saved value left over from the previous call.
!
!    Local, integer X_LO_INDEX, X_HI_INDEX, records the range of entries of
!    X that we need to compute.  This starts off as 1:N, but is adjusted
!    if we have a saved value that can be immediately stored in X(1),
!    and so on.
!
!    Local, real ( kind = rk ) Y, the value saved from the previous call, if
!    SAVED is 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer m
  integer, save :: made = 0
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) r(n+1)
  integer, save :: saved = 0
  real ( kind = rk ) x(n)
  integer x_hi_index
  integer x_lo_index
  real ( kind = rk ), save :: y = 0.0D+00
!
!  I'd like to allow the user to reset the internal data.
!  But this won't work properly if we have a saved value Y.
!  I'm making a crock option that allows the user to signal
!  explicitly that any internal memory should be flushed,
!  by passing in a negative value for N.
!
  if ( n < 0 ) then
    n = made
    made = 0
    saved = 0
    y = 0.0D+00
    return
  else if ( n == 0 ) then
    return
  end if
!
!  Record the range of X we need to fill in.
!
  x_lo_index = 1
  x_hi_index = n
!
!  Use up the old value, if we have it.
!
  if ( saved == 1 ) then
    x(1) = y
    saved = 0
    x_lo_index = 2
  end if
!
!  Maybe we don't need any more values.
!
  if ( x_hi_index - x_lo_index + 1 == 0 ) then
!
!  If we need just one new value, do that here to avoid null arrays.
!
  else if ( x_hi_index - x_lo_index + 1 == 1 ) then

    call random_number ( harvest = r(1:2) )

    x(x_hi_index) = &
             sqrt ( -2.0D+00 * log ( r(1) ) ) * cos ( 2.0D+00 * pi * r(2) )
    y =      sqrt ( -2.0D+00 * log ( r(1) ) ) * sin ( 2.0D+00 * pi * r(2) )

    saved = 1

    made = made + 2
!
!  If we require an even number of values, that's easy.
!
  else if ( mod ( x_hi_index - x_lo_index + 1, 2 ) == 0 ) then

    m = ( x_hi_index - x_lo_index + 1 ) / 2

    call random_number ( harvest = r(1:2*m) )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * cos ( 2.0D+00 * pi * r(2:2*m:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * sin ( 2.0D+00 * pi * r(2:2*m:2) )

    made = made + x_hi_index - x_lo_index + 1
!
!  If we require an odd number of values, we generate an even number,
!  and handle the last pair specially, storing one in X(N), and
!  saving the other for later.
!
  else

    x_hi_index = x_hi_index - 1

    m = ( x_hi_index - x_lo_index + 1 ) / 2 + 1

    call random_number ( harvest = r(1:2*m) )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * cos ( 2.0D+00 * pi * r(2:2*m-2:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * sin ( 2.0D+00 * pi * r(2:2*m-2:2) )

    x(n) = sqrt ( -2.0D+00 * log ( r(2*m-1) ) ) &
      * cos ( 2.0D+00 * pi * r(2*m) )

    y = sqrt ( -2.0D+00 * log ( r(2*m-1) ) ) &
      * sin ( 2.0D+00 * pi * r(2*m) )

    saved = 1

    made = made + x_hi_index - x_lo_index + 2

  end if

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
