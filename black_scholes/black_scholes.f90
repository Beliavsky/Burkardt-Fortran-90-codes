subroutine asset_path ( s0, mu, sigma, t1, n, s )

!*****************************************************************************80
!
!! asset_path() simulates the behavior of an asset price over time.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    Original MATLAB version by Desmond Higham.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Desmond Higham,
!    Black-Scholes for Scientific Computing Students,
!    Computing in Science and Engineering,
!    November/December 2004, Volume 6, Number 6, pages 72-79.
!
!  Input:
!
!    real ( kind = rk ) S0, the asset price at time 0.
!
!    real ( kind = rk ) MU, the expected growth rate.
!
!    real ( kind = rk ) SIGMA, the volatility of the asset.
!
!    real ( kind = rk ) T1, the expiry date.
!
!    integer N, the number of steps to take 
!    between 0 and T1.
!
!  Output:
!
!    real ( kind = rk ) S(0:N), the option values from time 0 to T1 
!    in equal steps.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) dt
  integer i
  real ( kind = rk ) mu
  real ( kind = rk ) p
  real ( kind = rk ) r(n)
  real ( kind = rk ) s(0:n)
  real ( kind = rk ) s0
  real ( kind = rk ) sigma
  real ( kind = rk ) t1

  dt = t1 / real ( n, kind = rk )

  call r8vec_normal_01 ( n, r )

  s(0) = s0
  p = s0
  do i = 1, n
    p = p * exp ( ( mu - sigma * sigma ) * dt + sigma * sqrt ( dt ) * r(i) )
    s(i) = p
  end do

  return
end
subroutine binomial ( s0, e, r, sigma, t1, m, c )

!*****************************************************************************80
!
!! binomial() uses the binomial method for a European call.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    Original MATLAB version by Desmond Higham.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Desmond Higham,
!    Black-Scholes for Scientific Computing Students,
!    Computing in Science and Engineering,
!    November/December 2004, Volume 6, Number 6, pages 72-79.
!
!  Input:
!
!    real ( kind = rk ) S0, the asset price at time 0.
!
!    real ( kind = rk ) E, the exercise price.
!
!    real ( kind = rk ) R, the interest rate.
!
!    real ( kind = rk ) SIGMA, the volatility of the asset.
!
!    real ( kind = rk ) T1, the expiry date.
!
!    integer M, the number of steps to take 
!    between 0 and T1.
!
!  Output:
!
!    real ( kind = rk ) C, the option value at time 0.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) c
  real ( kind = rk ) d
  real ( kind = rk ) dt
  real ( kind = rk ) e
  integer i
  integer m
  integer n
  real ( kind = rk ) p
  real ( kind = rk ) r
  real ( kind = rk ) s0
  real ( kind = rk ) sigma
  real ( kind = rk ) t1
  real ( kind = rk ) u
  real ( kind = rk ) w(1:m+1)
!
!  Time stepsize.
!
  dt = t1 / real ( m, kind = rk )

  a = 0.5D+00 * ( exp ( - r * dt ) + exp ( ( r + sigma**2 ) * dt ) )

  d = a - sqrt ( a * a - 1.0D+00 )
  u = a + sqrt ( a * a - 1.0D+00 )

  p = ( exp ( r * dt ) - d ) / ( u - d )

  do i = 1, m + 1
    w(i) = max ( s0 * d**(m+1-i) * u**(i-1) - e, 0.0D+00 )
  end do
!
!  Trace backwards to get the option value at time 0.
!
  do n = m, 1, -1
    do i = 1, n
      w(i) = ( 1.0D+00 - p ) * w(i) + p * w(i+1)
    end do
  end do

  w(1:m+1) = exp ( - r * t1 ) * w(1:m+1)

  c = w(1)

  return
end
subroutine bsf ( s0, t0, e, r, sigma, t1, c )

!*****************************************************************************80
!
!! bsf() evaluates the Black-Scholes formula for a European call.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    Original MATLAB version by Desmond Higham.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Desmond Higham,
!    Black-Scholes for Scientific Computing Students,
!    Computing in Science and Engineering,
!    November/December 2004, Volume 6, Number 6, pages 72-79.
!
!  Input:
!
!    real ( kind = rk ) S0, the asset price at time T0.
!
!    real ( kind = rk ) T0, the time at which the asset price is known.
!
!    real ( kind = rk ) E, the exercise price.
!
!    real ( kind = rk ) R, the interest rate.
!
!    real ( kind = rk ) SIGMA, the volatility of the asset.
!
!    real ( kind = rk ) T1, the expiry date.
!
!  Output:
!
!    real ( kind = rk ) C, the value of the call option.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) c
  real ( kind = rk ) d1
  real ( kind = rk ) d2
  real ( kind = rk ) e
  real ( kind = rk ) n1
  real ( kind = rk ) n2
  real ( kind = rk ) r
  real ( kind = rk ) s0
  real ( kind = rk ) sigma
  real ( kind = rk ) t0
  real ( kind = rk ) t1
  real ( kind = rk ) tau

  tau = t1 - t0

  if ( 0.0D+00 < tau ) then

    d1 = ( log ( s0 / e ) + ( r + 0.5D+00 * sigma * sigma ) * tau ) &
      / ( sigma * sqrt ( tau ) )

    d2 = d1 - sigma * sqrt ( tau )

    n1 = 0.5D+00 * ( 1.0D+00 + erf ( d1 / sqrt ( 2.0D+00 ) ) )
    n2 = 0.5D+00 * ( 1.0D+00 + erf ( d2 / sqrt ( 2.0D+00 ) ) )

    c = s0 * n1 - e * exp ( - r * tau ) * n2

  else

    c = max ( s0 - e, 0.0D+00 )

  end if

  return
end
subroutine forward ( e, r, sigma, t1, nx, nt, smax, u )

!*****************************************************************************80
!
!! forward() uses the forward difference method to value a European call option.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    Original MATLAB version by Desmond Higham.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Desmond Higham,
!    Black-Scholes for Scientific Computing Students,
!    Computing in Science and Engineering,
!    November/December 2004, Volume 6, Number 6, pages 72-79.
!
!  Input:
!
!    real ( kind = rk ) E, the exercise price.
!
!    real ( kind = rk ) R, the interest rate.
!
!    real ( kind = rk ) SIGMA, the volatility of the asset.
!
!    real ( kind = rk ) T1, the expiry date.
!
!    integer NX, the number of "space" steps used to 
!    divide the interval [0,L].
!
!    integer NT, the number of time steps.
!
!    real ( kind = rk ) SMAX, the maximum value of S to consider.
!
!  Output:
!
!    real ( kind = rk ) U(NX-1,NT+1), the value of the European 
!    call option.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nt
  integer nx

  real ( kind = rk ) a(2:nx-1)
  real ( kind = rk ) b(1:nx-1)
  real ( kind = rk ) c(1:nx-2)
  real ( kind = rk ) dt
  real ( kind = rk ) dx
  real ( kind = rk ) e
  integer i
  integer j
  real ( kind = rk ) p
  real ( kind = rk ) r
  real ( kind = rk ) sigma
  real ( kind = rk ) smax
  real ( kind = rk ) t
  real ( kind = rk ) t1
  real ( kind = rk ) u(nx-1,nt+1)
  real ( kind = rk ) u0

  dt = t1 / real ( nt, kind = rk )
  dx = smax / real ( nx, kind = rk )

  do i = 1, nx - 1
    b(i) = 1.0D+00 - r * dt - dt * ( sigma * i )**2
  end do

  do i = 1, nx - 2
    c(i) = 0.5D+00 * dt * ( sigma * i )**2 + 0.5D+00 * dt * r * i
  end do

  do i = 2, nx - 1
    a(i) = 0.5D+00 * dt * ( sigma * i )**2 - 0.5D+00 * dt * r * i
  end do

  u0 = 0.0D+00
  do i = 1, nx - 1
    u0 = u0 + dx
    u(i,1) = max ( u0 - e, 0.0D+00 )
  end do
  
  do j = 1, nt

    t = real ( j - 1, kind = rk ) * t1 / real ( nt, kind = rk )

    p = 0.5D+00 * dt * ( nx - 1 ) * ( sigma * sigma * ( nx - 1 ) + r ) &
      * ( smax - e * exp ( - r * t ) )

    u(1:nx-1,j+1) =                 b(1:nx-1) * u(1:nx-1,j)
    u(1:nx-2,j+1) = u(1:nx-2,j+1) + c(1:nx-2) * u(2:nx-1,j)
    u(2:nx-1,j+1) = u(2:nx-1,j+1) + a(2:nx-1) * u(1:nx-2,j)

    u(nx-1,j+1) = u(nx-1,j+1) + p

  end do

  return
end
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
!    03 September 2021
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
subroutine mc ( s0, e, r, sigma, t1, m, conf )

!*****************************************************************************80
!
!! mc() uses Monte Carlo valuation on a European call.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    Original MATLAB version by Desmond Higham.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Desmond Higham,
!    Black-Scholes for Scientific Computing Students,
!    Computing in Science and Engineering,
!    November/December 2004, Volume 6, Number 6, pages 72-79.
!
!  Input:
!
!    real ( kind = rk ) S0, the asset price at time 0.
!
!    real ( kind = rk ) E, the exercise price.
!
!    real ( kind = rk ) R, the interest rate.
!
!    real ( kind = rk ) SIGMA, the volatility of the asset.
!
!    real ( kind = rk ) T1, the expiry date.
!
!    integer M, the number of simulations.
!
!  Output:
!
!    real ( kind = rk ) CONF(2), the estimated range of the valuation.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m

  real ( kind = rk ) conf(2)
  real ( kind = rk ) e
  real ( kind = rk ) pmean
  real ( kind = rk ) pvals(m)
  real ( kind = rk ) r
  real ( kind = rk ) s0
  real ( kind = rk ) sigma
  real ( kind = rk ) std
  real ( kind = rk ) svals(m)
  real ( kind = rk ) t1
  real ( kind = rk ) u(m)
  real ( kind = rk ) width

  call r8vec_normal_01 ( m, u )

  svals(1:m) = s0 * exp ( ( r - 0.5D+00 * sigma * sigma ) * t1 &
    + sigma * sqrt ( t1 ) * u(1:m) )

  pvals(1:m) = exp ( - r * t1 ) * max ( svals(1:m) - e, 0.0D+00 )

  pmean = sum ( pvals(1:m) ) / real ( m, kind = rk )

  std = sqrt ( sum ( ( pvals(1:m) - pmean )**2 ) / real ( m - 1, kind = rk ) )

  width = 1.96D+00 * std / sqrt ( real ( m, kind = rk ) )

  conf(1) = pmean - width
  conf(2) = pmean + width

  return
end
function r8_normal_01 ( )

!*****************************************************************************80
!
!! r8_normal_01() returns a unit pseudonormal R8.
!
!  Discussion:
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!    Because this routine uses the Box Muller method, it requires pairs
!    of uniform random values to generate a pair of normal random values.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Output:
!
!    real ( kind = rk ) R8_NORMAL_01, a sample of the standard normal PDF.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) r1
  real ( kind = rk ) r2
  real ( kind = rk ) r8_normal_01
  integer, save :: used = 0
  real ( kind = rk ) x
  real ( kind = rk ), save :: y = 0.0D+00
!
!  On odd numbered calls, generate two uniforms, create two normals,
!  return the first normal.
!
  if ( mod ( used, 2 ) == 0 ) then

    call random_number ( harvest = r1 )
    call random_number ( harvest = r2 )

    x = sqrt ( - 2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * pi * r2 )
    y = sqrt ( - 2.0D+00 * log ( r1 ) ) * sin ( 2.0D+00 * pi * r2 )
!
!  On odd calls, return the second normal.
!
  else

    x = y

  end if

  used = used + 1

  r8_normal_01 = x

  return
end
subroutine r8vec_normal_01 ( n,  x )

!*****************************************************************************80
!
!! r8vec_normal_01() returns a unit pseudonormal R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
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
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of values desired.  If N is
!    negative, then the code will flush its internal memory; in particular,
!    if there is a saved value to be used on the next call, it is
!    instead discarded.
!
!  Output:
!
!    real ( kind = rk ) X(N), a sample of the standard normal PDF.
!
!  Local:
!
!    integer MADE, records the number of values that have
!    been computed.  On input with negative N, this value overwrites
!    the return value of N, so the user can get an accounting of
!    how much work has been done.
!
!    real ( kind = rk ) R(N+1), is used to store some uniform
!    random values.  Its dimension is N+1, but really it is only needed
!    to be the smallest even number greater than or equal to N.
!
!    integer SAVED, is 0 or 1 depending on whether there
!    is a single saved value left over from the previous call.
!
!    integer X_LO_INDEX, X_HI_INDEX, records the range
!    of entries of X that we need to compute.  This starts off as 1:N, but
!    is adjusted if we have a saved value that can be immediately stored
!    in X(1), and so on.
!
!    real ( kind = rk ) Y, the value saved from the previous call, if
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
             sqrt ( - 2.0D+00 * log ( r(1) ) ) * cos ( 2.0D+00 * pi * r(2) )
    y =      sqrt ( - 2.0D+00 * log ( r(1) ) ) * sin ( 2.0D+00 * pi * r(2) )

    saved = 1

    made = made + 2
!
!  If we require an even number of values, that's easy.
!
  else if ( mod ( x_hi_index - x_lo_index, 2 ) == 1 ) then

    m = ( x_hi_index - x_lo_index + 1 ) / 2

    call random_number ( harvest = r(1:2*m) )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( - 2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * cos ( 2.0D+00 * pi * r(2:2*m:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( - 2.0D+00 * log ( r(1:2*m-1:2) ) ) &
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
      sqrt ( - 2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * cos ( 2.0D+00 * pi * r(2:2*m-2:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( - 2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * sin ( 2.0D+00 * pi * r(2:2*m-2:2) )

    x(n) = sqrt ( - 2.0D+00 * log ( r(2*m-1) ) ) &
      * cos ( 2.0D+00 * pi * r(2*m) )

    y = sqrt ( - 2.0D+00 * log ( r(2*m-1) ) ) &
      * sin ( 2.0D+00 * pi * r(2*m) )

    saved = 1

    made = made + x_hi_index - x_lo_index + 2

  end if

  return
end
subroutine r8vec_print_part ( n, a, max_print, title )

!*****************************************************************************80
!
!! r8vec_print_part() prints "part" of an R8VEC.
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vector, is no more than MAX_PRINT, then
!    the entire vector is printed, one entry per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of entries of the vector.
!
!    real ( kind = rk ) A(N), the vector to be printed.
!
!    integer MAX_PRINT, the maximum number of lines
!    to print.
!
!    character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  integer i
  integer max_print
  character ( len = * ) title

  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  if ( n <= max_print ) then

    do i = 1, n
      write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)
    end do

  else if ( 3 <= max_print ) then

    do i = 1, max_print - 2
      write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)
    end do
    write ( *, '(a)' ) '  ........  ..............'
    i = n
    write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)

  else

    do i = 1, max_print - 1
      write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)
    end do
    i = max_print
    write ( *, '(2x,i8,a,1x,g14.6,2x,a)' ) i, ':', a(i), '...more entries...'

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
!    03 September 2021
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
 
