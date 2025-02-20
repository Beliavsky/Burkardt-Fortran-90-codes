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
subroutine i4vec_mean ( n, a, mean )

!*****************************************************************************80
!
!! I4VEC_MEAN returns the mean of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, integer A(N), the vector whose mean is desired.
!
!    Output, real ( kind = rk ) MEAN, the mean of the vector entries.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer a(n)
  real ( kind = rk ) mean

  mean = real ( sum ( a(1:n) ), kind = rk ) &
       / real ( n, kind = rk )

  return
end
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 May 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, integer A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer n

  integer a(n)
  integer i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,a,2x,i12)' ) i, ':', a(i)
  end do

  return
end
subroutine i4vec_variance ( n, a, variance )

!*****************************************************************************80
!
!! I4VEC_VARIANCE returns the variance of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, integer A(N), the vector whose variance is desired.
!
!    Output, real ( kind = rk ) VARIANCE, the variance of the vector entries.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer a(n)
  real ( kind = rk ) mean
  real ( kind = rk ) variance

  if ( n < 2 ) then

    variance = 0.0D+00

  else

    mean = real ( sum ( a(1:n) ), kind = rk ) / real ( n, kind = rk )

    variance = sum ( ( real ( a(1:n), kind = rk ) - mean )**2 )

    variance = variance / real ( n - 1, kind = rk )

  end if

  return
end
subroutine poisson_fixed_events ( lambda, event_num, t, w )

!*****************************************************************************80
!
!! POISSON_FIXED_EVENTS waits for a given number of Poisson events.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) LAMBDA, the average number of events per 
!    unit time.
!
!    Input, integer EVENT_NUM, the number of events to wait for.
!
!    Output, real ( kind = rk ) T(0:EVENT_NUM), the time at which a total 
!    of 0, 1, 2, ... and EVENT_NUM events were observed.
!
!    Output, real ( kind = rk ) W(0:EVENT_NUM), the waiting time until the
!    I-th event occurred.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer event_num

  real ( kind = rk ) lambda
  real ( kind = rk ) t(0:event_num)
  real ( kind = rk ) w(0:event_num)
!
!  Poisson waiting times follow an exponential distribution.
!
  w(0) = 0.0D+00
  call random_number ( harvest = w(1:event_num) )
  w(1:event_num) = - log ( w(1:event_num) ) / lambda
!
!  The time til event I is the sum of the waiting times 0 through I.
!
  call r8vec_cum ( event_num + 1, w, t )

  return
end
subroutine poisson_fixed_time ( lambda, time, n )

!*****************************************************************************80
!
!! POISSON_FIXED_TIME counts the Poisson events in a fixed time.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) LAMBDA, the average number of events 
!    per unit time.
!
!    Input, real ( kind = rk ) TIME, the amount of time to observe.
!
!    Output, integer N, the number of Poisson events observed.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) dt
  real ( kind = rk ) lambda
  integer n
  real ( kind = rk ) t
  real ( kind = rk ) time
  real ( kind = rk ) u

  n = 0
  t = 0.0D+00

  do while ( t < time )
    call random_number ( harvest = u )
    dt = - log ( u ) / lambda
    n = n + 1
    t = t + dt
  end do

  return
end
subroutine r8vec_cum ( n, a, a_cum )

!*****************************************************************************80
!
!! R8VEC_CUM computes the cumulutive sums of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    Input:
!
!      A = (/ 1.0, 2.0, 3.0, 4.0 /)
!
!    Output:
!
!      A_CUM = (/ 1.0, 3.0, 6.0, 10.0 /)
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real ( kind = rk ) A(N), the vector to be summed.
!
!    Output, real ( kind = rk ) A_CUM(1:N), the cumulative sums.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) a_cum(n)
  integer i

  a_cum(1) = a(1)

  do i = 2, n
    a_cum(i) = a_cum(i-1) + a(i)
  end do

  return
end
subroutine r8vec_mean ( n, a, mean )

!*****************************************************************************80
!
!! R8VEC_MEAN returns the mean of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real ( kind = rk ) A(N), the vector whose mean is desired.
!
!    Output, real ( kind = rk ) MEAN, the mean of the vector entries.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) mean

  mean = sum ( a(1:n) ) / real ( n, kind = rk )

  return
end
subroutine r8vec_midspace ( n, a, b, x )

!*****************************************************************************80
!
!! R8VEC_MIDSPACE creates a vector of linearly spaced values.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    This function divides the interval [a,b] into n subintervals, and then
!    returns the midpoints of those subintervals.
!
!  Example:
!
!    N = 5, A = 10, B = 20
!    X = [ 11, 13, 15, 17, 19 ]
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 June 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real ( kind = rk ) A, B, the endpoints of the interval.
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

  do i = 1, n
    x(i) = ( real ( 2 * n - 2 * i + 1, kind = rk ) * a &
           + real (         2 * i - 1, kind = rk ) * b ) &
           / real ( 2 * n,             kind = rk )
  end do

  return
end
subroutine r8vec_variance ( n, a, variance )

!*****************************************************************************80
!
!! R8VEC_VARIANCE returns the variance of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The variance of a vector X of length N is defined as
!
!      mean ( X(1:n) ) = sum ( X(1:n) ) / n
!
!      var ( X(1:n) ) = sum ( ( X(1:n) - mean )^2 ) / ( n - 1 )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!    N should be at least 2.
!
!    Input, real ( kind = rk ) A(N), the vector.
!
!    Output, real ( kind = rk ) VARIANCE, the variance of the vector.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) mean
  real ( kind = rk ) variance

  if ( n < 2 ) then

    variance = 0.0D+00

  else

    mean = sum ( a(1:n) ) / real ( n, kind = rk )

    variance = sum ( ( a(1:n) - mean )**2 )

    variance = variance / real ( n - 1, kind = rk )

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

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
