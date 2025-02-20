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
subroutine naca4_cambered ( m, p, t, c, n, xc, xu, yu, xl, yl )

!*****************************************************************************80
!
!! naca4_cambered(): (xu,yu), (xl,yl) for a NACA cambered 4-digit airfoil.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 May 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Eastman Jacobs, Kenneth Ward, Robert Pinkerton,
!    "The characteristics of 78 related airfoil sections from tests in
!    the variable-density wind tunnel",
!    NACA Report 460, 1933.
!
!  Input:
!
!    real ( kind = rk8 ) M, the maximum camber.
!    0.0 < M.
!
!    real ( kind = rk8 ) P, the location of maximum camber.
!    0.0 < P < 1.0
!
!    real ( kind = rk8 ) T, the maximum relative thickness.
!    0.0 < T <= 1.0
!
!    real ( kind = rk8 ) C, the chord length.
!    0.0 < C.
!
!    integer N, the number of sample points.
!
!    real ( kind = rk8 ) XC(N), points along the chord length.  
!    0.0 <= XC(*) <= C.
!
!  Output:
!
!    real ( kind = rk8 ) XU(N), YU(N), XL(N), YL(N), for each value of 
!    XC, measured along the camber line, the corresponding values (XU,YU) 
!    on the upper airfoil surface and (XL,YL) on the lower airfoil surface.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  real ( kind = rk8 ) c
  real ( kind = rk8 ) divisor
  real ( kind = rk8 ) dycdx
  integer i
  real ( kind = rk8 ) m
  real ( kind = rk8 ) p
  real ( kind = rk8 ) t
  real ( kind = rk8 ) theta
  real ( kind = rk8 ) xc(n)
  real ( kind = rk8 ) xl(n)
  real ( kind = rk8 ) xu(n)
  real ( kind = rk8 ) yc
  real ( kind = rk8 ) yl(n)
  real ( kind = rk8 ) yt
  real ( kind = rk8 ) yu(n)

  do i = 1, n

    if ( 0.0D+00 <= xc(i) / c .and. xc(i) / c <= p ) then
      divisor = p ** 2
    else if ( p <= xc(i) / c .and. xc(i) / c <= 1.0D+00 ) then
      divisor = ( 1.0D+00 - p ) ** 2
    else
      divisor = 1.0D+00
    end if

    dycdx = 2.0D+00 * m * ( p - xc(i) / c ) / divisor

    theta = atan ( dycdx )
   
    yt = 5.0D+00 * t * c * ( &
       0.2969D+00 * sqrt ( xc(i) / c ) &
       + (((( &
         - 0.1015D+00 ) * ( xc(i) / c ) &
         + 0.2843D+00 ) * ( xc(i) / c ) &
         - 0.3516D+00 ) * ( xc(i) / c ) &
         - 0.1260D+00 ) * ( xc(i) / c ) )

    if ( 0.0D+00 <= xc(i) / c .and. xc(i) / c <= p ) then
      yc = m * xc(i) * ( 2.0D+00 * p - xc(i) / c ) / p ** 2
    else if ( p <= xc(i) / c .and. xc(i) / c <= 1.0D+00 ) then
      yc = m * ( xc(i) - c ) * ( 2.0D+00 * p - xc(i) / c - 1.0D+00 ) &
        / ( 1.0D+00 - p ) ** 2
    else
      yc = 0.0D+00
    end if

    xu(i) = xc(i) - yt * sin ( theta )
    yu(i) = yc + yt * cos ( theta )
    xl(i) = xc(i) + yt * sin ( theta )
    yl(i) = yc - yt * cos ( theta )

  end do

  return
end
subroutine naca4_symmetric ( t, c, n, x, y )

!*****************************************************************************80
!
!! naca4_symmetric() evaluates y(x) for a NACA symmetric 4-digit airfoil.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 May 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Eastman Jacobs, Kenneth Ward, Robert Pinkerton,
!    "The characteristics of 78 related airfoil sections from tests in
!    the variable-density wind tunnel",
!    NACA Report 460, 1933.
!
!  Input:
!
!    real ( kind = rk8 ) T, the maximum relative thickness.
!
!    real ( kind = rk8 ) C, the chord length.
!
!    integer N, the number of sample points.
!
!    real ( kind = rk8 ) X(N), points along the chord length.  
!    0.0 <= X(*) <= C.
!
!  Output:
!
!    real ( kind = rk8 ) Y(N), for each value of X, the corresponding
!    value of Y so that (X,Y) is on the upper wing surface, and (X,-Y) is on the
!    lower wing surface.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  real ( kind = rk8 ) c
  real ( kind = rk8 ) t
  real ( kind = rk8 ) x(n)
  real ( kind = rk8 ) y(n)

  y(1:n) = 5.0D+00 * t * c * ( &
    0.2969D+00 * sqrt ( x(1:n) / c ) &
    + (((( &
      - 0.1015D+00 ) * ( x(1:n) / c ) &
      + 0.2843D+00 ) * ( x(1:n) / c ) &
      - 0.3516D+00 ) * ( x(1:n) / c ) &
      - 0.1260D+00 ) * ( x(1:n) / c ) )

  return
end
subroutine r8mat_write ( output_filename, m, n, table )

!*****************************************************************************80
!
!! r8mat_write() writes an R8MAT file.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    character ( len = * ) OUTPUT_FILENAME, the output file name.
!
!    integer M, the spatial dimension.
!
!    integer N, the number of points.
!
!    real ( kind = rk8 ) TABLE(M,N), the data.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer m
  integer n

  integer j
  character ( len = * ) output_filename
  integer output_status
  integer output_unit
  character ( len = 30 ) string
  real ( kind = rk8 ) table(m,n)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_filename ) // '" on unit ', output_unit
    output_unit = -1
    stop 1
  end if
!
!  Create a format string.
!
!  For less precision in the output file, try:
!
!                                            '(', m, 'g', 14, '.', 6, ')'
!
  if ( 0 < m .and. 0 < n ) then

    write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 24, '.', 16, ')'
!
!  Write the data.
!
    do j = 1, n
      write ( output_unit, string ) table(1:m,j)
    end do

  end if
!
!  Close the file.
!
  close ( unit = output_unit )

  return
end
subroutine r8vec_linspace ( n, a, b, x )

!*****************************************************************************80
!
!! r8vec_linspace() creates a vector of linearly spaced values.
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
!  Input:
!
!    integer N, the number of entries in the vector.
!
!    real ( kind = rk8 ) A, B, the first and last entries.
!
!  Output:
!
!    real ( kind = rk8 ) X(N), a vector of linearly spaced data.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  integer i
  real ( kind = rk8 ) x(n)

  if ( n == 1 ) then

    x(1) = ( a + b ) / 2.0D+00

  else

    do i = 1, n
      x(i) = ( real ( n - i,     kind = rk8 ) * a   &
             + real (     i - 1, kind = rk8 ) * b ) &
             / real ( n     - 1, kind = rk8 )
    end do

  end if

  return
end
function r8vec_max ( n, a )

!*****************************************************************************80
!
!! r8vec_max() returns the maximum value in an R8VEC.
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
!    30 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of entries in the array.
!
!    real ( kind = rk8 ) A(N), the array.
!
!  Output:
!
!    real ( kind = rk8 ) R8VEC_MAX, the value of the largest entry.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  real ( kind = rk8 ) a(n)
  real ( kind = rk8 ) r8vec_max
  real ( kind = rk8 ) value

  value = maxval ( a(1:n) )

  r8vec_max = value

  return
end
function r8vec_min ( n, a )

!*****************************************************************************80
!
!! r8vec_min() returns the minimum value of an R8VEC.
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
!    17 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of entries in the array.
!
!    real ( kind = rk8 ) A(N), the array.
!
!  Output:
!
!    real ( kind = rk8 ) R8VEC_MIN, the value of the smallest entry.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  real ( kind = rk8 ) a(n)
  real ( kind = rk8 ) r8vec_min
  real ( kind = rk8 ) value

  value = minval ( a(1:n) )

  r8vec_min = value

  return
end

