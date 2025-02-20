program main

!*****************************************************************************80
!
!! NEAREST_INTERP_1D_TEST() tests NEAREST_INTERP_1D().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 May 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ni
  integer prob
  integer prob_num

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NEAREST_INTERP_1D_TEST:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the NEAREST_INTERP_1D library.'
  write ( *, '(a)' ) '  The R8LIB library is needed.'
  write ( *, '(a)' ) '  The test needs the TEST_INTERP library.'

  call p00_prob_num ( prob_num )

  ni = 11
  do prob = 1, prob_num
    call nearest_interp_1d_test01 ( prob, ni )
  end do

  do prob = 1, prob_num
    call nearest_interp_1d_test02 ( prob )
  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NEAREST_INTERP_1D_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine nearest_interp_1d_test01 ( prob, ni )

!*****************************************************************************80
!
!! NEAREST_INTERP_1D_TEST01 tests NEAREST_INTERP_1D.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer PROB, the index of the problem.
!
!    Input, integer NI, the number of interpolation points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: d(:,:)
  integer ni
  integer nd
  integer prob
  character ( len = 80 ) title
  real ( kind = rk ), allocatable :: xd(:)
  real ( kind = rk ), allocatable :: xi(:)
  real ( kind = rk ) xd_max
  real ( kind = rk ) xd_min
  real ( kind = rk ), allocatable :: yd(:)
  real ( kind = rk ), allocatable :: yi(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NEAREST_INTERP_1D_TEST01'
  write ( *, '(a,i2)' ) &
    '  Sample the nearest neighbor interpolant for problem # ', prob

  call p00_data_num ( prob, nd )

  allocate ( d(1:2,1:nd) )
  allocate ( xd(1:nd) )
  allocate ( yd(1:nd) )

  call p00_data ( prob, 2, nd, d )

  xd(1:nd) = d(1,1:nd)
  yd(1:nd) = d(2,1:nd)

  xd_min = minval ( xd(1:nd) )
  xd_max = maxval ( xd(1:nd) )

  allocate ( xi(1:ni) )
  allocate ( yi(1:ni) )

  call r8vec_linspace ( ni, xd_min, xd_max, xi )
  call nearest_interp_1d ( nd, xd, yd, ni, xi, yi )

  write ( title, '(a,i2)' ) '  X, Y for problem ', prob

  call r8vec2_print ( ni, xi, yi, title )

  deallocate ( d )
  deallocate ( xd )
  deallocate ( xi )
  deallocate ( yd )
  deallocate ( yi )

  return
end
subroutine nearest_interp_1d_test02 ( prob )

!*****************************************************************************80
!
!! NEAREST_INTERP_1D_TEST02 tests NEAREST_INTERP_1D.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer PROB, the index of the problem.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = 255 ) command_filename
  integer command_unit
  character ( len = 255 ) data_filename
  integer data_unit
  real ( kind = rk ) interp_error
  character ( len = 255 ) interp_filename
  integer interp_unit
  integer j
  integer ni
  integer nd
  character ( len = 255 ) output_filename
  integer prob
  real ( kind = rk ) r8vec_diff_norm
  real ( kind = rk ), allocatable :: xd(:)
  real ( kind = rk ), allocatable :: xi(:)
  real ( kind = rk ) xmax
  real ( kind = rk ) xmin
  real ( kind = rk ), allocatable :: xy(:,:)
  real ( kind = rk ), allocatable :: yd(:)
  real ( kind = rk ), allocatable :: yi(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'NEAREST_INTERP_1D_TEST02:'
  write ( *, '(a,i2)' ) '  Interpolate data from TEST_INTERP problem #', prob

  call p00_data_num ( prob, nd )
  write ( *, '(a,i6)' ) '  Number of data points = ', nd

  allocate ( xy(1:2,1:nd ) )

  call p00_data ( prob, 2, nd, xy )
  
  call r8mat_transpose_print ( 2, nd, xy, '  Data array:' )

  allocate ( xd(1:nd) )
  allocate ( yd(1:nd) )
  xd(1:nd) = xy(1,1:nd)
  yd(1:nd) = xy(2,1:nd)
!
!  #1:  Does interpolant match function at interpolation points?
!
  ni = nd
  allocate ( xi(1:ni) )
  allocate ( yi(1:ni) )
  xi(1:ni) = xd(1:ni)
  call nearest_interp_1d ( nd, xd, yd, ni, xi, yi )

  interp_error = r8vec_diff_norm ( ni, yi, yd ) / real ( ni, kind = rk )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) &
    '  Node-averaged L2 interpolation error = ', interp_error
!
!  Create data file.
!
  write ( data_filename, '(a,i2.2,a)' ) 'data', prob, '.txt'
  call get_unit ( data_unit )
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
  deallocate ( xi )
  deallocate ( yi )

  ni = 501
  call r8vec_min ( nd, xd, xmin )
  call r8vec_max ( nd, xd, xmax )
  allocate ( xi(1:ni) )
  allocate ( yi(1:ni) )
  call r8vec_linspace ( ni, xmin, xmax, xi )
  call nearest_interp_1d ( nd, xd, yd, ni, xi, yi )

  write ( interp_filename, '(a,i2.2,a)' ) 'interp', prob, '.txt'
  call get_unit ( interp_unit )
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
  call get_unit ( command_unit )
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
    'set title "Data versus Nearest Neighbor Interpolant"'
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
  deallocate ( xy )
  deallocate ( yd )
  deallocate ( yi )

  return
end
