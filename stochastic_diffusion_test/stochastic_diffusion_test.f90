program main

!*****************************************************************************80
!
!! stochastic_diffusion_test() tests stochastic_diffusion().
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 July 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'stochastic_diffusion_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  test the stochastic_diffusion().'

  call bnt_contour ( )
  call elman_contour ( )
  call ntw_contour ( )
  call xk_contour ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'STOCHASTIC_DIFFUSION_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine bnt_contour ( )

!*****************************************************************************80
!
!! BNT_CONTOUR displays contour plots of a 2D stochastic diffusivity function.
!
!  Discussion:
!
!    The diffusivity function is compute by DIFFUSIVITY_2D_BNT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 July 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Ivo Babuska, Fabio Nobile, Raul Tempone,
!    A stochastic collocation method for elliptic partial differential equations
!    with random input data,
!    SIAM Journal on Numerical Analysis,
!    Volume 45, Number 3, 2007, pages 1005-1034.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: m = 4
  integer, parameter :: nx = 41
  integer, parameter :: ny = 31

  character ( len = 255 ) command_filename
  integer command_unit
  character ( len = 255 ) data_filename
  integer data_unit
  real ( kind = rk8 ) dc(nx,ny)
  real ( kind = rk8 ) dc0
  integer i
  integer j
  integer n
  real ( kind = rk8 ) omega(m)
  real ( kind = rk8 ) xmat(nx,ny)
  real ( kind = rk8 ) xvec(nx)
  real ( kind = rk8 ) ymat(nx,ny)
  real ( kind = rk8 ) yvec(ny)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'BNT_CONTOUR'
  write ( *, '(a)' ) '  Display contour or surface plots of the stochastic'
  write ( *, '(a)' ) '  diffusivity function defined by DIFFUSIVITY_2D_BNT.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  The first plot uses uniform random values for OMEGA.'
  write ( *, '(a)' ) '  The second uses Gaussian (normal) random values.'
!
!  Set the spatial grid.
!
  call r8vec_linspace ( nx, -1.5D+00, 0.0D+00, xvec )
  call r8vec_linspace ( ny, -0.4D+00, 0.8D+00, yvec )

  call r8vec_mesh_2d ( nx, ny, xvec, yvec, xmat, ymat )
!
!  Sample OMEGA.
!
  call random_number ( harvest = omega(1:m) )
!
!  Compute the diffusivity field.
!
  dc0 = 10.0D+00
  n = nx * ny
  call diffusivity_2d_bnt ( dc0, omega, n, xmat, ymat, dc )
!
!  Create a data file.
!
  data_filename = 'bnt_data.txt'
  call get_unit ( data_unit )
  open ( unit = data_unit, file = data_filename, status = 'replace' )
  do j = 1, ny
    do i = 1, nx
      write ( data_unit, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) xmat(i,j), ymat(i,j), dc(i,j)
    end do
    write ( data_unit, '(a)' ) ''
  end do
  close ( unit = data_unit )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Created graphics data file "' // trim ( data_filename ) // '".'
!
!  Create the command file.
!
  command_filename = 'bnt_commands.txt'
  call get_unit ( command_unit )
  open ( unit = command_unit, file = command_filename, status = 'replace' )

  write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) 'set output "bnt_contour.png"'
  write ( command_unit, '(a)' ) 'set xlabel "<---X--->"'
  write ( command_unit, '(a)' ) 'set ylabel "<---Y--->"'
  write ( command_unit, '(a)' ) 'set zlabel "<---DC(X,Y)--->"'
  write ( command_unit, '(a)' ) &
    'set title "BNT Stochastic diffusivity function"'
  write ( command_unit, '(a)' ) 'set contour'
  write ( command_unit, '(a)' ) 'set timestamp'
  write ( command_unit, '(a)' ) 'set cntrparam levels 10'
  write ( command_unit, '(a)' ) '#set view map'
  write ( command_unit, '(a)' ) 'set view 75, 75'
  write ( command_unit, '(a)' ) 'unset key'
  write ( command_unit, '(a)' ) 'splot "' // trim ( data_filename ) // '"'

  close ( unit = command_unit )

  write ( *, '(a)' ) &
    '  Created graphics command file "' // trim ( command_filename ) // '".'

  return
end
subroutine elman_contour ( )

!*****************************************************************************80
!
!! ELMAN_CONTOUR displays a contour plot of a 2D stochastic diffusivity function.
!
!  Discussion:
!
!    The diffusivity function is compute by DIFFUSIVITY_2D_ELMAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 July 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Howard Elman, Darran Furnaval,
!    Solving the stochastic steady-state diffusion problem using multigrid,
!    IMA Journal on Numerical Analysis,
!    Volume 27, Number 4, 2007, pages 675-688.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: m_1d = 5
  integer, parameter :: nx = 51
  integer, parameter :: ny = 51

  real ( kind = rk8 ) a
  real ( kind = rk8 ) cl
  character ( len = 255 ) command_filename
  integer command_unit
  character ( len = 255 ) data_filename
  integer data_unit
  real ( kind = rk8 ) dc(nx,ny)
  real ( kind = rk8 ) dc0
  integer i
  integer j
  real ( kind = rk8 ) omega(m_1d,m_1d)
  real ( kind = rk8 ) xmat(nx,ny)
  real ( kind = rk8 ) xvec(nx)
  real ( kind = rk8 ) ymat(nx,ny)
  real ( kind = rk8 ) yvec(ny)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'ELMAN_CONTOUR'
  write ( *, '(a)' ) '  Display contour or surface plots of the stochastic'
  write ( *, '(a)' ) '  diffusivity function defined by DIFFUSIVITY_2D_ELMAN.'
!
!  Set the spatial grid.
!
  a = 1.0D+00
  call r8vec_linspace ( nx, -a, +a, xvec )
  call r8vec_linspace ( ny, -a, +a, yvec )

  call r8vec_mesh_2d ( nx, ny, xvec, yvec, xmat, ymat )
!
!  Sample OMEGA.
!
  call r8vec_normal_01 ( m_1d * m_1d, omega )
!
!  Compute the diffusivity field.
!
  cl = 0.1D+00
  dc0 = 10.0D+00
  call diffusivity_2d_elman ( a, cl, dc0, m_1d, omega, nx, nx, xmat, ymat, dc )
!
!  Create a data file.
!
  data_filename = 'elman_data.txt'
  call get_unit ( data_unit )
  open ( unit = data_unit, file = data_filename, status = 'replace' )
  do j = 1, ny
    do i = 1, nx
      write ( data_unit, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) xmat(i,j), ymat(i,j), dc(i,j)
    end do
    write ( data_unit, '(a)' ) ''
  end do
  close ( unit = data_unit )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Created graphics data file "' // trim ( data_filename ) // '".'
!
!  Create the command file.
!
  command_filename = 'elman_commands.txt'
  call get_unit ( command_unit )
  open ( unit = command_unit, file = command_filename, status = 'replace' )

  write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) 'set output "elman_contour.png"'
  write ( command_unit, '(a)' ) 'set xlabel "<---X--->"'
  write ( command_unit, '(a)' ) 'set ylabel "<---Y--->"'
  write ( command_unit, '(a)' ) 'set zlabel "<---DC(X,Y)--->"'
  write ( command_unit, '(a)' ) &
    'set title "Elman Stochastic diffusivity function"'
  write ( command_unit, '(a)' ) 'set contour'
  write ( command_unit, '(a)' ) 'set timestamp'
  write ( command_unit, '(a)' ) 'set cntrparam levels 10'
  write ( command_unit, '(a)' ) '#set view map'
  write ( command_unit, '(a)' ) 'set view 75, 75'
  write ( command_unit, '(a)' ) 'unset key'
  write ( command_unit, '(a)' ) 'splot "' // trim ( data_filename ) // '"'

  close ( unit = command_unit )

  write ( *, '(a)' ) &
    '  Created graphics command file "' // trim ( command_filename ) // '".'

  return
end
subroutine ntw_contour ( )

!*****************************************************************************80
!
!! NTW_CONTOUR displays a contour plot of a 2D stochastic diffusivity function.
!
!  Discussion:
!
!    The diffusivity function is compute by DIFFUSIVITY_2D_NTW.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 July 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Fabio Nobile, Raul Tempone, Clayton Webster,
!    A Sparse Grid Stochastic Collocation Method for Partial Differential
!    Equations with Random Input Data,
!    SIAM Journal on Numerical Analysis,
!    Volume 46, Number 5, 2008, pages 2309-2345.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: m = 21
  integer, parameter :: nx = 101
  integer, parameter :: ny = 101

  real ( kind = rk8 ) cl
  character ( len = 255 ) command_filename
  integer command_unit
  real ( kind = rk8 ) d
  character ( len = 255 ) data_filename
  integer data_unit
  real ( kind = rk8 ), allocatable :: dc(:,:)
  real ( kind = rk8 ) dc0
  integer i
  integer j
  real ( kind = rk8 ) omega(m)
  real ( kind = rk8 ), allocatable :: xmat(:,:)
  real ( kind = rk8 ) xvec(nx)
  real ( kind = rk8 ), allocatable :: ymat(:,:)
  real ( kind = rk8 ) yvec(ny)

  allocate ( dc(1:nx,1:ny) )
  allocate ( xmat(1:nx,1:ny) )
  allocate ( ymat(1:nx,1:ny) )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'NTW_CONTOUR'
  write ( *, '(a)' ) '  Display contour or surface plots of the stochastic'
  write ( *, '(a)' ) '  diffusivity function defined by DIFFUSIVITY_2D_NTW.'
!
!  Set the spatial grid.
!
  d = 1.0D+00
  call r8vec_linspace ( nx, 0.0D+00, d, xvec )
  call r8vec_linspace ( ny, 0.0D+00, d, yvec )

  call r8vec_mesh_2d ( nx, ny, xvec, yvec, xmat, ymat )
!
!  Sample OMEGA.
!  We rescale to  [-sqrt(3),sqrt(3)].
!
  call random_number ( harvest = omega(1:m) )
  omega = ( 1.0D+00 - omega ) * ( - sqrt ( 3.0D+00 ) ) &
        +             omega   *     sqrt ( 3.0D+00 )
!
!  Evaluate the diffusivity field.
!
  cl = 0.1D+00
  dc0 = 0.5D+00
  call diffusivity_2d_ntw ( cl, dc0, m, omega, nx * ny, xmat, ymat, dc )
!
!  Create a data file.
!
  data_filename = 'ntw_data.txt'
  call get_unit ( data_unit )
  open ( unit = data_unit, file = data_filename, status = 'replace' )
  do j = 1, ny
    do i = 1, nx
      write ( data_unit, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) xmat(i,j), ymat(i,j), dc(i,j)
    end do
    write ( data_unit, '(a)' ) ''
  end do
  close ( unit = data_unit )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Created graphics data file "' // trim ( data_filename ) // '".'
!
!  Create the command file.
!
  command_filename = 'ntw_commands.txt'
  call get_unit ( command_unit )
  open ( unit = command_unit, file = command_filename, status = 'replace' )

  write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) 'set output "ntw_contour.png"'
  write ( command_unit, '(a)' ) 'set xlabel "<---X--->"'
  write ( command_unit, '(a)' ) 'set ylabel "<---Y--->"'
  write ( command_unit, '(a)' ) 'set zlabel "<---DC(X,Y)--->"'
  write ( command_unit, '(a)' ) &
    'set title "NTW Stochastic diffusivity function"'
  write ( command_unit, '(a)' ) 'set contour'
  write ( command_unit, '(a)' ) 'set timestamp'
  write ( command_unit, '(a)' ) 'set cntrparam levels 15'
  write ( command_unit, '(a)' ) '#set view map'
  write ( command_unit, '(a)' ) 'set view 65, 65'
  write ( command_unit, '(a)' ) 'set key'
  write ( command_unit, '(a)' ) 'splot "' // trim ( data_filename ) // '"'

  close ( unit = command_unit )

  write ( *, '(a)' ) &
    '  Created graphics command file "' // trim ( command_filename ) // '".'

  deallocate ( dc )
  deallocate ( xmat )
  deallocate ( ymat )

  return
end
subroutine xk_contour ( )

!*****************************************************************************80
!
!! xk_contour() displays contour plots of a 1D stochastic diffusivity function.
!
!  Discussion:
!
!    The diffusivity function is compute by DIFFUSIVITY_1D_XK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 July 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu, George Karniadakis,
!    Modeling uncertainty in steady state diffusion problems via
!    generalized polynomial chaos,
!    Computer Methods in Applied Mechanics and Engineering,
!    Volume 191, 2002, pages 4927-4948.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ), allocatable :: dc(:)
  real ( kind = rk8 ) dc_max
  real ( kind = rk8 ) dc0
  character ( len = 255 ) command_filename
  integer command_unit
  character ( len = 255 ) data_filename
  integer data_unit
  integer j
  integer m
  integer n
  real ( kind = rk8 ) r8vec_max
  real ( kind = rk8 ), allocatable :: omega(:)
  real ( kind = rk8 ), allocatable :: x(:)
  real ( kind = rk8 ) x_min
  real ( kind = rk8 ) x_max

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'XK_CONTOUR'
  write ( *, '(a)' ) '  Plot the stochastic diffusivity function'
  write ( *, '(a)' ) '  defined by DIFFUSIVITY_1D_XK.'
!
!  Set up the spatial grid.
!
  n = 51
  x_min = -1.0D+00
  x_max = +1.0D+00
  allocate ( x(1:n) )
  call r8vec_linspace ( n, x_min, x_max, x )
!
!  Sample the OMEGA values.
!
  m = 5
  allocate ( omega(1:m) )
  call r8vec_normal_01 ( m, omega )
!
!  Compute the diffusivity field.
!
  dc0 = 10.0D+00
  allocate ( dc(1:n) )
  call diffusivity_1d_xk ( dc0, m, omega, n, x, dc )
!
!  Create data file.
!
  data_filename = 'xk_data.txt'
  call get_unit ( data_unit )
  open ( unit = data_unit, file = data_filename, status = 'replace' )
  do j = 1, n
    write ( data_unit, '(2x,g14.6,2x,g14.6)' ) x(j), dc(j)
  end do
  close ( unit = data_unit )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Created graphics data file "' // trim ( data_filename ) // '".'
!
!  Create the command file.
!
  dc_max = r8vec_max ( n, dc )

  command_filename = 'xk_commands.txt'
  call get_unit ( command_unit )
  open ( unit = command_unit, file = command_filename, status = 'replace' )

  write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) 'set output "xk_contour.png"'
  write ( command_unit, '(a)' ) 'set xlabel "<---X--->"'
  write ( command_unit, '(a)' ) 'set ylabel "<---DC(X)--->"'
  write ( command_unit, '(a,g14.6,a)' ) 'set yrange [0.0:', dc_max, ']'
  write ( command_unit, '(a)' ) &
    'set title "XK Stochastic diffusivity function"'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set style data lines'
  write ( command_unit, '(a)' ) 'plot "' // trim ( data_filename ) // &
    '" using 1:2 lw 3 linecolor rgb "red"'

  close ( unit = command_unit )

  write ( *, '(a)' ) &
    '  Created graphics command file "' // trim ( command_filename ) // '".'

  deallocate ( dc )
  deallocate ( omega )
  deallocate ( x )

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

