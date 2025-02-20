subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
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
subroutine grid_2d ( x_num, x_lo, x_hi, y_num, y_lo, y_hi, x, y )

!*****************************************************************************80
!
!! GRID_2D returns a regular 2D grid.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 January 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer X_NUM, the number of X values to use.
!
!    Input, real ( kind = rk ) X_LO, X_HI, the range of X values.
!
!    Input, integer Y_NUM, the number of Y values to use.
!
!    Input, real ( kind = rk ) Y_LO, Y_HI, the range of Y values.
!
!    Output, real ( kind = rk ) X(X_NUM,Y_NUM), Y(X_NUM,Y_NUM), 
!    the coordinates of the grid.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer x_num
  integer y_num

  integer i
  integer j
  real ( kind = rk ) x(x_num,y_num)
  real ( kind = rk ) x_hi
  real ( kind = rk ) x_lo
  real ( kind = rk ) xi
  real ( kind = rk ) y(x_num,y_num)
  real ( kind = rk ) y_hi
  real ( kind = rk ) y_lo
  real ( kind = rk ) yj

  if ( x_num == 1 ) then
    x(1:x_num,1:y_num) = ( x_lo + x_hi ) / 2.0D+00
  else
    do i = 1, x_num
      xi = ( real ( x_num - i,     kind = rk ) * x_lo   &
           + real (         i - 1, kind = rk ) * x_hi ) &
           / real ( x_num     - 1, kind = rk )
      x(i,1:y_num) = xi
    end do
  end if

  if ( y_num == 1 ) then
    y(1:x_num,1:y_num) = ( y_lo + y_hi ) / 2.0D+00
  else
    do j = 1, y_num
      yj = ( real ( y_num - j,     kind = rk ) * y_lo   &
           + real (         j - 1, kind = rk ) * y_hi ) &
           / real ( y_num     - 1, kind = rk )
      y(1:x_num,j) = yj
    end do
  end if

  return
end
subroutine r8vec_uniform_ab ( n, a, b, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_AB returns a scaled pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    Each dimension ranges from A to B.
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
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real ( kind = rk ) A, B, the lower and upper limits.
!
!    Input/output, integer SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = rk ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a
  real ( kind = rk ) b
  integer i
  integer, parameter :: i4_huge = 2147483647
  integer k
  integer seed
  real ( kind = rk ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_AB - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    r(i) = a + ( b - a ) * real ( seed, kind = rk ) * 4.656612875D-10

  end do

  return
end
subroutine resid_spiral ( n, x, y, c, pr )

!*****************************************************************************80
!
!! RESID_SPIRAL computes the residual for a spiral velocity vector field.
!
!  Discussion:
!
!    Note that the continuous velocity field (U,V)(X,Y) that is discretely
!    sampled here satisfies the homogeneous continuity equation, that is,
!    it has zero divergence.  In other words:
!
!      dU/dX + dV/dY = 0.
!
!    This is by construction, since we have
!
!      U(X,Y) =  10 * d/dY ( PHI(X) * PHI(Y) )
!      V(X,Y) = -10 * d/dX ( PHI(X) * PHI(Y) )
!
!    which guarantees zero divergence.
!
!    The underlying function PHI is defined by
!
!      PHI(Z) = ( 1 - cos ( C * pi * Z ) ) * ( 1 - Z )^2
!
!    where C is a parameter.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 January 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of evaluation points.
!
!    Input, real ( kind = rk )  X(N), Y(N), the coordinates of the 
!    evaluation points.
!
!    Input, real ( kind = rk ) C, a parameter, typically between 0 and 2 * PI.
!
!    Output, real ( kind = rk ) PR(N), the residual in the continuity equation.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) c
  integer i
  real ( kind = rk ) pr(n)
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) u
  real ( kind = rk ) ux
  real ( kind = rk ) v
  real ( kind = rk ) vy
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)

  do i = 1, n

    u =   10.0D+00 * ( 1.0D+00 - cos ( c * r8_pi * x(i) ) ) &
           * ( 1.0D+00 - x(i) ) ** 2 &
           * ( &
               c * r8_pi * sin ( c * r8_pi * y(i) ) * ( 1.0D+00 - y(i) ) ** 2 &
             - ( 1.0D+00 - cos ( c * r8_pi * y(i) ) ) &
             * 2.0D+00 * ( 1.0D+00 - y(i) ) &
             )

    ux =   10.0D+00 * &
      ( &
        c * r8_pi * sin ( c * r8_pi * x(i) ) * ( 1.0D+00 - x(i) ) ** 2 &
        - ( 1.0D+00 - cos ( c * r8_pi * x(i) ) ) &
        * 2.0D+00 * ( 1.0D+00 - x(i) ) &
      ) &
      * &
      ( &
        c * r8_pi * sin ( c * r8_pi * y(i) ) * ( 1.0D+00 - y(i) ) ** 2 &
        - ( 1.0D+00 - cos ( c * r8_pi * y(i) ) ) &
        * 2.0D+00 * ( 1.0D+00 - y(i) ) &
      )

    v = - 10.0D+00 * ( 1.0D+00 - cos ( c * r8_pi * y(i) ) ) &
      * ( 1.0D+00 - y(i) ) ** 2 &
      * ( &
          c * r8_pi * sin ( c * r8_pi * x(i) ) * ( 1.0D+00 - x(i) ) ** 2 &
        - ( 1.0D+00 - cos ( c * r8_pi * x(i) ) ) &
        * 2.0D+00 * ( 1.0D+00 - x(i) ) &
        )

    vy =  - 10.0D+00 * &
      ( &
        c * r8_pi * sin ( c * r8_pi * x(i) ) * ( 1.0D+00 - x(i) ) ** 2 &
        - ( 1.0D+00 - cos ( c * r8_pi * x(i) ) ) &
        * 2.0D+00 * ( 1.0D+00 - x(i) ) &
      ) &
      * &
      ( &
        c * r8_pi * sin ( c * r8_pi * y(i) ) * ( 1.0D+00 - y(i) ) ** 2 &
        - ( 1.0D+00 - cos ( c * r8_pi * y(i) ) ) &
        * 2.0D+00 * ( 1.0D+00 - y(i) ) &
      )

    pr(i) = ux + vy

  end do

  return
end
subroutine spiral_gnuplot ( header, n, x, y, u, v, s )

!*****************************************************************************80
!
!! SPIRAL_GNUPLOT writes the spiral vector field to files for GNUPLOT.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 January 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) HEADER, a header to be used to name the files.
!
!    Input, integer N, the number of evaluation points.
!
!    Input, real ( kind = rk ) X(N), Y(N), the coordinates of the evaluation points.
!
!    Input, real ( kind = rk ) U(N), V(N), the velocity components.
!
!    Input, real ( kind = rk ) S, a scale factor for the velocity vectors.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  character ( len = 255 ) command_filename
  integer command_unit
  character ( len = 255 ) data_filename
  integer data_unit
  character ( len = * ) header
  integer i
  character ( len = 255 ) plot_filename
  real ( kind = rk ) s
  real ( kind = rk ) u(n)
  real ( kind = rk ) v(n)
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)
!
!  Write the data file.
!
  data_filename = trim ( header ) // '_data.txt'

  call get_unit ( data_unit )

  open ( unit = data_unit, file = data_filename, status = 'replace' )

  do i = 1, n
    write ( data_unit, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      x(i), y(i), s * u(i), s * v(i)
  end do

  close ( unit = data_unit )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Data written to "' // trim ( data_filename ) // '".'
!
!  Write the command file.
!
  command_filename = trim ( header ) // '_commands.txt'
  call get_unit ( command_unit )

  plot_filename = trim ( header ) // '.png'

  open ( unit = command_unit, file = command_filename, status = 'replace' )

  write ( command_unit, '(a)' ) '#  ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) 'set output "' // trim ( plot_filename ) // '"'
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '#  Add titles and labels.'
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set xlabel "<--- X --->"'
  write ( command_unit, '(a)' ) 'set ylabel "<--- Y --->"'
  write ( command_unit, '(a)' ) 'set title "Spiral velocity flow"'
  write ( command_unit, '(a)' ) 'unset key'
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '#  Add grid lines.'
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set size ratio -1'
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '#  Timestamp the plot.'
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set timestamp'
  write ( command_unit, '(a)' ) 'plot "' // trim ( data_filename ) // &
    '" using 1:2:3:4 with vectors \'
  write ( command_unit, '(a)' ) '  head filled lt 2 linecolor rgb "blue"'
  write ( command_unit, '(a)' ) 'quit'

  close ( unit = command_unit )

  write ( *, '(a)' ) '  Commands written to "' // &
    trim ( command_filename ) // '".'

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
subroutine uv_spiral ( n, x, y, c, u, v )

!*****************************************************************************80
!
!! UV_SPIRAL computes a spiral velocity vector field.
!
!  Discussion:
!
!    Note that the continuous velocity field (U,V)(X,Y) that is discretely
!    sampled here satisfies the homogeneous continuity equation, that is,
!    it has zero divergence.  In other words:
!
!      dU/dX + dV/dY = 0.
!
!    This is by construction, since we have
!
!      U(X,Y) =  10 * d/dY ( PHI(X) * PHI(Y) )
!      V(X,Y) = -10 * d/dX ( PHI(X) * PHI(Y) )
!
!    which guarantees zero divergence.
!
!    The underlying function PHI is defined by
!
!      PHI(Z) = ( 1 - cos ( C * pi * Z ) ) * ( 1 - Z )^2
!
!    where C is a parameter.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 January 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of evaluation points.
!
!    Input, real ( kind = rk ) X(N), Y(N), the coordinates of the 
!    evaluation points.
!
!    Input, real ( kind = rk ) C, a parameter, typically between 0 and 2 * PI.
!
!    Output, real ( kind = rk ) U(N), V(N), the velocity components.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) c
  integer i
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) u(n)
  real ( kind = rk ) v(n)
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)

  do i = 1, n

    u(i) =   10.0D+00 * ( 1.0D+00 - cos ( c * r8_pi * x(i) ) ) &
      * ( 1.0D+00 - x(i) ) ** 2 &
      * ( &
          c * r8_pi * sin ( c * r8_pi * y(i) ) * ( 1.0D+00 - y(i) ) ** 2 &
        - ( 1.0D+00 - cos ( c * r8_pi * y(i) ) ) &
          * 2.0D+00 * ( 1.0D+00 - y(i) ) &
        )

    v(i) = - 10.0D+00 * ( 1.0D+00 - cos ( c * r8_pi * y(i) ) ) &
      * ( 1.0D+00 - y(i) ) ** 2 &
      * ( &
          c * r8_pi * sin ( c * r8_pi * x(i) ) * ( 1.0D+00 - x(i) ) ** 2 &
        - ( 1.0D+00 - cos ( c * r8_pi * x(i) ) ) &
        * 2.0D+00 * ( 1.0D+00 - x(i) ) &
        )

  end do

  return
end
