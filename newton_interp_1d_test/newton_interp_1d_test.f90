program main

!*****************************************************************************80
!
!! newton_interp_1d_test() tests newton_interp_1d().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 July 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer prob
  integer prob_num

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'newton_interp_1d_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test newton_interp_1d().'

  call newton_coef_1d_test ( )

  call newton_value_1d_test ( )

  call p00_prob_num ( prob_num )

  do prob = 1, prob_num
    call newton_interp_1d_test01 ( prob )
  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'newton_interp_1d_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine newton_coef_1d_test ( )

!*****************************************************************************80
!
!! newton_coef_1d_test() tests newton_coef_1d().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 July 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nd = 5

  real ( kind = rk ) cd(nd)
  real ( kind = rk ) xd(nd)
  real ( kind = rk ) yd(nd)

  xd = (/ 0.0, 1.0, 2.0, 3.0, 4.0 /)
  yd = (/ 24.0, 0.0, 0.0, 0.0, 0.0 /)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'newton_coef_1d_test():'
  write ( *, '(a)' ) '  newton_coef_1d() sets the coefficients for a 1D Newton interpolation.'

  call r8vec2_print ( nd, xd, yd, '  Interpolation data:' )

  call newton_coef_1d ( nd, xd, yd, cd )

  call r8vec_print ( nd, cd, '  Newton interpolant coefficients:' )

  return
end
subroutine newton_value_1d_test ( )

!*****************************************************************************80
!
!! newton_value_1d_test() tests newton_value_1d().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 July 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nd = 5
  integer, parameter :: ni = 16

  real ( kind = rk ) cd(nd)
  real ( kind = rk ) x_hi
  real ( kind = rk ) x_lo
  real ( kind = rk ) xd(nd)
  real ( kind = rk ) xi(ni)
  real ( kind = rk ) yi(ni)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'newton_value_1d_test():'
  write ( *, '(a)' ) '  newton_value_1d() evaluates a Newton 1d interpolant.'

  xd = (/ 0.0, 1.0, 2.0, 3.0, 4.0 /)
  cd = (/ 24.0, -24.0, +12.0, -4.0, 1.0 /)
  call r8vec2_print ( nd, xd, cd, '  The Newton interpolant data:' )

  x_lo = 0.0
  x_hi = 5.0
  call r8vec_linspace ( ni, x_lo, x_hi, xi )

  call newton_value_1d ( nd, xd, cd, ni, xi, yi )

  call r8vec2_print ( ni, xi, yi, '  Newton interpolant values:' )

  return
end
subroutine newton_interp_1d_test01 ( prob )

!*****************************************************************************80
!
!! newton_interp_1d_test01() tests newton_interp_1d().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 July 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: cd(:)
  character ( len = 255 ) command_filename
  integer command_unit
  character ( len = 255 ) data_filename
  integer data_unit
  character ( len = 255 ) interp_filename
  real ( kind = rk ) interp_error
  integer interp_unit
  integer j
  real ( kind = rk ) ld
  real ( kind = rk ) li
  integer nd
  integer ni
  character ( len = 255 ) output_filename
  integer prob
  real ( kind = rk ) r8vec_norm_affine
  real ( kind = rk ), allocatable :: xd(:)
  real ( kind = rk ), allocatable :: xi(:)
  real ( kind = rk ) xmax
  real ( kind = rk ) xmin
  real ( kind = rk ), allocatable :: xy(:,:)
  real ( kind = rk ), allocatable :: yd(:)
  real ( kind = rk ), allocatable :: yi(:)
  real ( kind = rk ) ymax
  real ( kind = rk ) ymin

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'newton_interp_1d_test01():'
  write ( *, '(a,i1)' ) '  Interpolate data from TEST_INTERP problem #', prob

  call p00_data_num ( prob, nd )
  write ( *, '(a,i3)' ) '  Number of data points = ', nd

  allocate ( xy(1:2,1:nd) )
  call p00_data ( prob, 2, nd, xy )
  
  allocate ( xd(1:nd) )
  allocate ( yd(1:nd) )
  xd(1:nd) = xy(1,1:nd)
  yd(1:nd) = xy(2,1:nd)
  deallocate ( xy )

  call r8vec2_print ( nd, xd, yd, '  X, Y data:' )
!
!  Get the Newton coefficients.
!
  allocate ( cd(1:nd) )
  call newton_coef_1d ( nd, xd, yd, cd )
!
!  #1:  Does interpolant match function at interpolation points?
!
  ni = nd
  allocate ( xi(1:ni) )
  allocate ( yi(1:ni) )
  xi(1:ni) = xd(1:ni)

  call newton_value_1d ( nd, xd, cd, ni, xi, yi )

  interp_error = r8vec_norm_affine ( ni, yi, yd ) / real ( ni, kind = rk )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  L2 interpolation error averaged per interpolant node = ', interp_error

  deallocate ( xi )
  deallocate ( yi )
!
!  #2: Compare estimated curve length to piecewise linear (minimal) curve length.
!  Assume data is sorted, and normalize X and Y dimensions by (XMAX-XMIN) and
!  (YMAX-YMIN).
!
  xmin = minval ( xd(1:nd) )
  xmax = maxval ( xd(1:nd) )
  ymin = minval ( yd(1:nd) )
  ymax = maxval ( yd(1:nd) )

  ni = 501
  allocate ( xi(1:ni) )
  allocate ( yi(1:ni) )
  call r8vec_linspace ( ni, xmin, xmax, xi )
  call newton_value_1d ( nd, xd, cd, ni, xi, yi )

  ld = sum ( sqrt ( ( ( xd(2:nd) - xd(1:nd-1) ) / ( xmax - xmin ) ) ** 2 &
                  + ( ( yd(2:nd) - yd(1:nd-1) ) / ( ymax - ymin ) ) ** 2 ) )

  li = sum ( sqrt ( ( ( xi(2:ni) - xi(1:ni-1) ) / ( xmax - xmin ) ) ** 2 &
                  + ( ( yi(2:ni) - yi(1:ni-1) ) / ( ymax - ymin ) ) ** 2 ) )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Normalized length of piecewise linear interpolant = ', ld
  write ( *, '(a,g14.6)' ) '  Normalized length of Newton interpolant           = ', li

  deallocate ( xi )
  deallocate ( yi )
!
!  Create data file.
!
  write ( data_filename, '(a,i2.2,a)' ) 'data', prob, '.txt'
! call get_unit ( data_unit )
  data_unit = 99
  open ( unit = data_unit, file = data_filename, status = 'replace' )
  do j = 1, nd
    write ( data_unit, '(2x,g14.6,2x,g14.6)' ) xd(j), yd(j)
  end do
  close ( unit = data_unit )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Created graphics data file "' // trim ( data_filename ) // '".'
!
!  Create interp file.
!
  ni = 501
  allocate ( xi(1:ni) )
  allocate ( yi(1:ni) )
  call r8vec_linspace ( ni, xmin, xmax, xi )
  call newton_value_1d ( nd, xd, cd, ni, xi, yi )

  write ( interp_filename, '(a,i2.2,a)' ) 'interp', prob, '.txt'
! call get_unit ( interp_unit )
  interp_unit = 99
  open ( unit = interp_unit, file = interp_filename, status = 'replace' )
  do j = 1, ni
    write ( interp_unit, '(2x,g14.6,2x,g14.6)' ) xi(j), yi(j)
  end do
  close ( unit = interp_unit )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Created graphics interp file "' // trim ( interp_filename ) // '".'
!
!  Plot the data and the interpolant.
!
  write ( command_filename, '(a,i2.2,a)' ) 'commands', prob, '.txt'
! call get_unit ( command_unit )
  command_unit = 99
  open ( unit = command_unit, file = command_filename, status = 'replace' )

  write ( output_filename, '(a,i2.2,a)' ) 'plot', prob, '.png'

  write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) 'set output "' // trim ( output_filename ) // '"'
  write ( command_unit, '(a)' ) 'set xlabel "<---X--->"'
  write ( command_unit, '(a)' ) 'set ylabel "<---Y--->"'
  write ( command_unit, '(a)' ) &
    'set title "Data versus Newton Polynomial Interpolant"'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set style data lines'
  write ( command_unit, '(a)' ) 'plot "' // trim ( data_filename ) // &
    '" using 1:2 with points pt 7 ps 2 lc rgb "blue",\'
  write ( command_unit, '(a)' ) '     "' // trim ( interp_filename ) // &
    '" using 1:2 lw 3 linecolor rgb "red"'

  close ( unit = command_unit )
  write ( *, '(a)' ) &
    '  Created graphics command file "' // trim ( command_filename ) // '".'

  deallocate ( xd )
  deallocate ( xi )
  deallocate ( yd )
  deallocate ( yi )

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
!    real ( kind = rk ) A, B, the first and last entries.
!
!  Output:
!
!    real ( kind = rk ) X(N), a vector of linearly spaced data.
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
function r8vec_norm_affine ( n, v0, v1 )

!*****************************************************************************80
!
!! r8vec_norm_affine() returns the affine norm of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The affine vector L2 norm is defined as:
!
!      R8VEC_NORM_AFFINE(V0,V1)
!        = sqrt ( sum ( 1 <= I <= N ) ( V1(I) - V0(I) )^2
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the order of the vectors.
!
!    real ( kind = rk ) V0(N), the base vector.
!
!    real ( kind = rk ) V1(N), the vector whose affine norm is desired.
!
!  Output:
!
!    real ( kind = rk ) R8VEC_NORM_AFFINE, the L2 norm of V1-V0.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) r8vec_norm_affine
  real ( kind = rk ) v0(n)
  real ( kind = rk ) v1(n)

  r8vec_norm_affine = sqrt ( sum ( ( v0(1:n) - v1(1:n) ) ** 2 ) )

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! r8vec_print() prints an R8VEC.
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
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of components of the vector.
!
!    real ( kind = rk ) A(N), the vector to be printed.
!
!    character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  integer i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
  end do

  return
end
subroutine r8vec2_print ( n, a1, a2, title )

!*****************************************************************************80
!
!! r8vec2_print() prints an R8VEC2.
!
!  Discussion:
!
!    An R8VEC2 is a dataset consisting of N pairs of R8's, stored
!    as two separate vectors A1 and A2.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of components of the vector.
!
!    real ( kind = rk ) A1(N), A2(N), the vectors to be printed.
!
!    character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a1(n)
  real ( kind = rk ) a2(n)
  integer i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i4,2x,g14.6,2x,g14.6)' ) i, a1(i), a2(i)
  end do

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
!    15 August 2021
!
!  Author:
!
!    John Burkardt
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

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end

