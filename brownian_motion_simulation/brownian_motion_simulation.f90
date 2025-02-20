subroutine brownian_displacement_display ( k, n, m, d, t, dsq, header )

!*****************************************************************************80
!
!! brownian_displacement_display() displays average Brownian motion displacement.
!
!  Discussion:
!
!    Thanks to Feifei Xu for pointing out a missing factor of 2 in the
!    displacement calculation.
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
!    integer K, the number of repetitions.
!
!    integer N, the number of time steps.  
!
!    integer M, the spatial dimension. 
!
!    real ( kind = rk ) D, the diffusion coefficient.
!
!    real ( kind = rk ) T, the total time.
!
!    real ( kind = rk ) DSQ(K,N), the displacements over time for 
!    each repetition.
!
!    character ( len = * ) HEADER, an identifier for the output files.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer k
  integer n

  character ( len = 255 ) command_filename
  integer command_unit
  real ( kind = rk ) d
  character ( len = 255 ) data_filename
  integer data_unit
  real ( kind = rk ) dsq(k,n)
  real ( kind = rk ) dsq_ave
  real ( kind = rk ) dsq_ideal
  character ( len = * ) header
  integer i4_uniform_ab
  integer ii(5)
  integer j
  integer m
  real ( kind = rk ) t
  real ( kind = rk ) ti
!
!  Choose 5 paths at random.
!
  do j = 1, 5
    ii(j) = i4_uniform_ab ( 1, k )
  end do
!
!  Create the data file.
!
  call get_unit ( data_unit )

  data_filename = trim ( header ) // '_data.txt'

  open ( unit = data_unit, file = data_filename, status = 'replace' )

  do j = 1, n
    ti = real ( j - 1, kind = rk ) * t / real ( n - 1, kind = rk )
    dsq_ave = sum ( dsq(1:k,j) ) / real ( k, kind = rk )
    dsq_ideal = 2.0D+00 * real ( m, kind = rk ) * d * ti
    write ( data_unit, '(8(2x,g14.6))' ) &
      ti, dsq(ii(1:5),j), dsq_ave, dsq_ideal
  end do

  close ( unit = data_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  BROWNIAN_DISPLACEMENT data stored in "' &
    // trim ( data_filename ) // '".'
!
!  Create the command file.
!
  call get_unit ( command_unit )

  command_filename = trim ( header ) // '_commands.txt'

  open ( unit = command_unit, file = command_filename, status = 'replace' )

  write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) 'set output "' // trim ( header ) // '.png"'
  write ( command_unit, '(a)' ) 'set xlabel "T"'
  write ( command_unit, '(a)' ) 'set ylabel "D^2"'
  write ( command_unit, '(a)' ) &
    'set title "Squared displacement (Red), Predicted (Black), Samples (Blue)"'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set style data lines'
  write ( command_unit, '(a)' ) 'plot "' // trim ( data_filename ) // &
    '" using 1:2 title "sample 1" linecolor rgb "blue", \'
  write ( command_unit, '(a)' ) '     "' // trim ( data_filename ) // &
    '" using 1:3 title "sample 2" linecolor rgb "blue", \'
  write ( command_unit, '(a)' ) '     "' // trim ( data_filename ) // &
    '" using 1:4 title "sample 3" linecolor rgb "blue", \'
  write ( command_unit, '(a)' ) '     "' // trim ( data_filename ) // &
    '" using 1:5 title "sample 4" linecolor rgb "blue", \'
  write ( command_unit, '(a)' ) '     "' // trim ( data_filename ) // &
    '" using 1:6 title "sample 5" linecolor rgb "blue", \'
  write ( command_unit, '(a)' ) '     "' // trim ( data_filename ) // &
    '" using 1:7 title "Averaged" lw 3 linecolor rgb "red", \'
  write ( command_unit, '(a)' ) '     "' // trim ( data_filename ) // &
    '" using 1:8 title "Ideal" lw 3 linecolor rgb "black"'
  write ( command_unit, '(a)' ) 'quit'

  close ( unit = command_unit )

  write ( *, '(a)' ) '  BROWNIAN_DISPLACEMENT plot commands stored in "' &
    // trim ( command_filename ) // '".'

  return
end
subroutine brownian_displacement_simulation ( k, n, m, d, t, dsq )

!*****************************************************************************80
!
!! brownian_displacement_simulation() simulates Brownian displacement.
!
!  Discussion:
!
!    This function computes the square of the distance of the Brownian
!    particle from the starting point, repeating this calculation 
!    several times.
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
!    integer K, the number of repetitions.
!
!    integer N, the number of time steps to take, plus 1.
!
!    integer M, the spatial dimension. 
!
!    real ( kind = rk ) D, the diffusion coefficient. 
!    Computationally, this is simply a scale factor between time and space.
!
!    real ( kind = rk ) T, the total time.
!
!  Output:
!
!    real ( kind = rk ) DSQ(K,N), the displacements over time for each 
!    repetition.  DSQ(:,1) is 0.0, because we include the displacement at the 
!    initial time. 
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer k
  integer m
  integer n

  real ( kind = rk ) d
  real ( kind = rk ) dsq(k,n)
  integer i
  real ( kind = rk ) t
  real ( kind = rk ) x(m,n)

  do i = 1, k

    call brownian_motion_simulation ( m, n, d, t, x )

    dsq(i,1:n) = sum ( x(1:m,1:n) ** 2, dim = 1 )

  end do

  return
end
subroutine brownian_motion_display ( m, n, x, header )

!*****************************************************************************80
!
!! brownian_motion_display() displays successive Brownian motion positions.
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
!    integer M, the spatial dimension.
!    M should be 1 or 2.
!
!    integer N, the number of time steps. 
!
!    real ( kind = rk ) X(M,N), the particle positions.
!
!    character ( len = * ) HEADER, an identifier for the output files.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  character ( len = 255 ) command_filename
  integer command_unit
  character ( len = 255 ) data_filename
  integer data_unit
  character ( len = * ) header
  integer i
  real ( kind = rk ) t
  real ( kind = rk ) x(m,n)

  if ( m /= 1 .and. m /= 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BROWNIAN_MOTION_DISPLAY - Fatal error!'
    write ( *, '(a)' ) '  This routine can only handle M = 1 or 2.'
    stop
  end if
!
!  Create the data file.
!
  call get_unit ( data_unit )

  data_filename = trim ( header ) // '_data.txt'

  open ( unit = data_unit, file = data_filename, status = 'replace' )

  if ( m == 1 ) then
    do i = 1, n
      t = real ( i - 1, kind = rk ) / real ( n - 1, kind = rk )
      write ( data_unit, '(2x,g14.6,2x,g14.6)' ) t, x(1,i)
    end do
  else if ( m == 2 ) then
    do i = 1, n
      write ( data_unit, '(2x,g14.6,2x,g14.6)' ) x(1,i), x(2,i)
    end do
  end if

  close ( unit = data_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  BROWNIAN_MOTION data stored in "' &
    // trim ( data_filename ) // '".'
!
!  Create the command file.
!
  call get_unit ( command_unit )

  command_filename = trim ( header ) // '_commands.txt'

  open ( unit = command_unit, file = command_filename, status = 'replace' )

  write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) 'set output "' // trim ( header ) // '.png"'
  write ( command_unit, '(a)' ) 'set xlabel "X"'
  write ( command_unit, '(a)' ) 'set ylabel "T"'
  write ( command_unit, '(a)' ) 'set title "Brownian motion in 1D"'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set style data lines'
  write ( command_unit, '(a)' ) 'plot "' // trim ( data_filename ) &
    // '" using 1:2'
  write ( command_unit, '(a)' ) 'quit'

  close ( unit = command_unit )

  write ( *, '(a)' ) '  BROWNIAN_MOTION plot commands stored in "' &
    // trim ( command_filename ) // '".'

  return
end
subroutine brownian_motion_simulation ( m, n, d, t, x )

!*****************************************************************************80
!
!! brownian_motion_simulation() simulates Brownian motion.
!
!  Discussion:
!
!    Thanks to Feifei Xu for pointing out a missing factor of 2 in the
!    stepsize calculation, 08 March 2016.
!
!    Thanks to Joerg Peter Pfannmoeller for pointing out a missing factor
!    of M in the stepsize calculation, 23 April 2018.
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
!    integer M, the spatial dimension.
!
!    integer N, the number of time steps to take, plus 1. 
!
!    real ( kind = rk ) D, the diffusion coefficient.  
!
!    real ( kind = rk ) T, the total time.
!
!  Output:
!
!    real ( kind = rk ) X(M,N), the successive locations of the particle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) d
  real ( kind = rk ) dt
  real ( kind = rk ) dx(m)
  integer j
  real ( kind = rk ) norm_dx
  real ( kind = rk ) r8_normal_01
  real ( kind = rk ) s
  real ( kind = rk ) t
  real ( kind = rk ) x(m,n)
!
!  Set the time step.
!
  dt = t / real ( n - 1, kind = rk )
!
!  Start at the origin.
!
  x(1:m,1) = 0.0D+00
!
!  Take N - 1 steps.
!
  do j = 2, n
!
!  S is the stepsize.
!
    s = sqrt ( 2.0D+00 * real ( m, kind = rk ) * d * dt ) * r8_normal_01 ( )
!
!  Direction DX is random, unit norm.
!
    if ( m == 1 ) then
      dx(1) = s
    else
      call r8vec_normal_01 ( m, dx )
      norm_dx = sqrt ( sum ( dx(1:m) ** 2 ) )
      dx(1:m) = s * dx(1:m) / norm_dx
    end if
!
!  Add the step to the current position.
!
    x(1:m,j) = x(1:m,j-1) + dx(1:m)

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
!  Input:
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
function i4_uniform_ab ( a, b )

!*****************************************************************************80
!
!! i4_uniform_ab() returns a scaled pseudorandom I4.
!
!  Discussion:
!
!    An I4 is an integer value.
!
!    The pseudorandom number will be scaled to be uniformly distributed
!    between A and B.
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
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Input:
!
!    integer A, B, the limits of the interval.
!
!  Output:
!
!    integer I4_UNIFORM_AB, a number between A and B.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer a
  integer b
  integer, parameter :: i4_huge = 2147483647
  integer i4_uniform_ab
  real r
  integer value

  call random_number ( harvest = r )
!
!  Scale R to lie between A-0.5 and B+0.5.
!
  r = ( 1.0E+00 - r ) * ( real ( min ( a, b ) ) - 0.5E+00 ) &
    +             r   * ( real ( max ( a, b ) ) + 0.5E+00 )
!
!  Use rounding to convert R to an integer between A and B.
!
  value = nint ( r )

  value = max ( value, min ( a, b ) )
  value = min ( value, max ( a, b ) )

  i4_uniform_ab = value

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
subroutine r8vec_normal_01 ( n, x )

!*****************************************************************************80
!
!! r8vec_normal_01() returns a unit pseudonormal R8VEC.
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
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of values desired.  If N is
!    negative,then the code will flush its internal memory; in particular,
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
!    integer SAVED, is 0 or 1 depending on whether there is a
!    single saved value left over from the previous call.
!
!    integer X_LO_INDEX, X_HI_INDEX, records the range of entries of
!    X that we need to compute.  This starts off as 1:N, but is adjusted
!    if we have a saved value that can be immediately stored in X(1),
!    and so on.
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
 
