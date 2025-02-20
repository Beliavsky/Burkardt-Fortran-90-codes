program main

!*****************************************************************************80
!
!! STOCHASTIC_HEAT2D_TEST tests the STOCHASTIC_HEAT2D library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 August 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'stochastic_heat2d_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  test  stochastic_heat2d().'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'stochastic_heat2d_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 plots a sample solution of a 2D stochastic diffusivity equation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  character ( len = 80 ) command_filename
  integer command_unit
  character ( len = 80 ) data_filename
  integer data_unit
  integer i
  integer j
  integer nx
  integer ny
  real ( kind = rk8 ) omega(4)
  real ( kind = rk8 ), external :: test01_f
  real ( kind = rk8 ), allocatable :: umat(:,:)
  real ( kind = rk8 ) u_mean
  real ( kind = rk8 ), allocatable :: xmat(:,:)
  real ( kind = rk8 ) xmax
  real ( kind = rk8 ) xmin
  real ( kind = rk8 ), allocatable :: xvec(:)
  real ( kind = rk8 ), allocatable :: ymat(:,:)
  real ( kind = rk8 ) ymax
  real ( kind = rk8 ) ymin
  real ( kind = rk8 ), allocatable :: yvec(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  Consider the steady heat equation in the unit square,'
  write ( *, '(a)' ) '  with 0 Dirichlet boundary conditions, '
  write ( *, '(a)' ) '  and a heat source term F that is a Gaussian centered at (0.60,0.80).'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Model the diffusivity coefficient as spatially varying,'
  write ( *, '(a)' ) '  with a stochastic dependence on parameters OMEGA(1:4),'
  write ( *, '(a)' ) '  as described in Babuska, Nobile, Tempone (BNT).'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Compute and display the solution U for a given choice'
  write ( *, '(a)' ) '  of the parameters OMEGA.'
!
!  Create the X and Y coordinate vectors.
!
  nx = 21
  xmin = 0.0D+00
  xmax = 1.0D+00
  allocate ( xvec(1:nx) )
  call r8vec_linspace ( nx, xmin, xmax, xvec )

  ny = 21
  ymin = 0.0D+00
  ymax = 1.0D+00
  allocate ( yvec(1:ny) )
  call r8vec_linspace ( ny, ymin, ymax, yvec )
!
!  Create the X and Y coordinate matrices.
!
  allocate ( xmat(1:nx,1:ny) )
  allocate ( ymat(1:nx,1:ny) )
  call r8vec_mesh_2d ( nx, ny, xvec, yvec, xmat, ymat )
!
!  Sample OMEGA:
!
  call r8vec_normal_01 ( 4, omega )
  omega(1:4) = 2.0D+00 * omega(1:4)

  call r8vec_print ( 4, omega, '  Sampled OMEGA values:' )
!
!  Solve the finite difference approximation to the steady 2D heat equation
!  for this set of OMEGA values.
!
  allocate ( umat(1:nx,1:ny) )

  call stochastic_heat2d ( omega, nx, ny, xvec, yvec, test01_f, umat )
!
!  Create a data file.
!
  data_filename = 'solution_data.txt'
  call get_unit ( data_unit )
  open ( unit = data_unit, file = data_filename, status = 'replace' )
  do j = 1, ny
    do i = 1, nx
      write ( data_unit, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) xmat(i,j), ymat(i,j), umat(i,j)
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
  command_filename = 'solution_commands.txt'
  call get_unit ( command_unit )
  open ( unit = command_unit, file = command_filename, status = 'replace' )

  write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) 'set output "solution.png"'
  write ( command_unit, '(a)' ) 'set xlabel "<---X--->"'
  write ( command_unit, '(a)' ) 'set ylabel "<---Y--->"'
  write ( command_unit, '(a)' ) 'set zlabel "<---U(X,Y)--->"'
  write ( command_unit, '(a)' ) &
    'set title "Sample Solution"'
  write ( command_unit, '(a)' ) 'set contour'
  write ( command_unit, '(a)' ) 'set timestamp'
  write ( command_unit, '(a)' ) 'set cntrparam levels 10'
  write ( command_unit, '(a)' ) 'set view 75, 75'
  write ( command_unit, '(a)' ) 'unset key'
  write ( command_unit, '(a)' ) 'splot "' // trim ( data_filename ) // '"'

  close ( unit = command_unit )

  write ( *, '(a)' ) &
    '  Created graphics command file "' // trim ( command_filename ) // '".'
!
!  Report the average value of U.
!
  u_mean = sum ( umat(1:nx,1:ny) ) / real ( nx * ny, kind = rk8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Mean value of U is ', u_mean
!
!  Free memory.
!
  deallocate ( umat )
  deallocate ( xmat )
  deallocate ( xvec )
  deallocate ( ymat )
  deallocate ( yvec )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 looks at mean temperature as a function of OMEGA(1) and OMEGA(2).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  character ( len = 80 ) command_filename
  integer command_unit
  character ( len = 80 ) data_filename
  integer data_unit
  integer i
  integer j
  integer nx
  integer ny
  real ( kind = rk8 ) omega(4)
  real ( kind = rk8 ), allocatable :: omega1_mat(:,:)
  real ( kind = rk8 ) omega1_max
  real ( kind = rk8 ) omega1_min
  integer omega1_num
  real ( kind = rk8 ), allocatable :: omega1_vec(:)
  real ( kind = rk8 ), allocatable :: omega2_mat(:,:)
  real ( kind = rk8 ) omega2_max
  real ( kind = rk8 ) omega2_min
  integer omega2_num
  real ( kind = rk8 ), allocatable :: omega2_vec(:)
  real ( kind = rk8 ), external :: test01_f
  real ( kind = rk8 ), allocatable :: umat(:,:)
  real ( kind = rk8 ), allocatable :: u_mean_mat(:,:)
  real ( kind = rk8 ) u_mean_max
  real ( kind = rk8 ), allocatable :: xmat(:,:)
  real ( kind = rk8 ) xmax
  real ( kind = rk8 ) xmin
  real ( kind = rk8 ), allocatable :: xvec(:)
  real ( kind = rk8 ), allocatable :: ymat(:,:)
  real ( kind = rk8 ) ymax
  real ( kind = rk8 ) ymin
  real ( kind = rk8 ), allocatable :: yvec(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST02:'
  write ( *, '(a)' ) '  Fix OMEGA(3) = 4, OMEGA(4) = 0, and'
  write ( *, '(a)' ) '  examine dependence of average temperature on OMEGA(1) and OMEGA(2)'
  write ( *, '(a)' ) '  over the range [-10,+10].'
!
!  Create the X and Y coordinate vectors.
!
  nx = 21
  xmin = 0.0D+00
  xmax = 1.0D+00
  allocate ( xvec(1:nx) )
  call r8vec_linspace ( nx, xmin, xmax, xvec )

  ny = 21
  ymin = 0.0D+00
  ymax = 1.0D+00
  allocate ( yvec(1:ny) )
  call r8vec_linspace ( ny, ymin, ymax, yvec )
!
!  Create the X and Y coordinate matrices.
!
  allocate ( xmat(1:nx,1:ny) )
  allocate ( ymat(1:nx,1:ny) )
  call r8vec_mesh_2d ( nx, ny, xvec, yvec, xmat, ymat )
!
!  Create OMEGA1 and OMEGA2 vectors.
!
  omega1_num = 21
  omega1_min = -10.0D+00
  omega1_max = +10.0D+00
  allocate ( omega1_vec(1:omega1_num))
  call r8vec_linspace ( omega1_num, omega1_min, omega1_max, omega1_vec )

  omega2_num = 21
  omega2_min = -10.0D+00
  omega2_max = +10.0D+00
  allocate ( omega2_vec(1:omega2_num) )
  call r8vec_linspace ( omega2_num, omega2_min, omega2_max, omega2_vec )
!
!  Create the OMEGA1 and OMEGA2 coordinate matrices.
!
  allocate ( omega1_mat(1:omega1_num,1:omega2_num) )
  allocate ( omega2_mat(1:omega1_num,1:omega2_num) )
  call r8vec_mesh_2d ( omega1_num, omega2_num, omega1_vec, omega2_vec, omega1_mat, omega2_mat )
!
!  Set OMEGA(3) and OMEGA(4).
!
  omega(3) = 4.0D+00
  omega(4) = 0.0D+00

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Omega(3) fixed at ', omega(3)
  write ( *, '(a,g14.6)' ) '  Omega(4) fixed at ', omega(4)
!
!  Solve the finite difference approximation to the steady 2D heat equation,
!  and save the mean value of the solution, which is a slightly biased
!  estimate of the heat integral over the unit square.
!
  allocate ( u_mean_mat(1:omega1_num,1:omega2_num) )
  allocate ( umat(1:nx,1:ny) )

  do j = 1, omega2_num
    omega(2) = omega2_vec(j)
    do i = 1, omega1_num
      omega(1) = omega1_vec(i)
      call stochastic_heat2d ( omega, nx, ny, xvec, yvec, test01_f, umat )
      u_mean_mat(i,j) = sum ( umat(1:nx,1:ny) ) / real ( nx * ny, kind = rk8 )
    end do
  end do
!
!  Create a data file.
!
  data_filename = 'umean_data.txt'
  call get_unit ( data_unit )
  open ( unit = data_unit, file = data_filename, status = 'replace' )
  do j = 1, ny
    do i = 1, nx
      write ( data_unit, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) &
        omega1_mat(i,j), omega2_mat(i,j), u_mean_mat(i,j)
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
  command_filename = 'umean_commands.txt'
  call get_unit ( command_unit )
  open ( unit = command_unit, file = command_filename, status = 'replace' )

  write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) 'set output "umean.png"'
  write ( command_unit, '(a)' ) 'set xlabel "<---OMEGA1--->"'
  write ( command_unit, '(a)' ) 'set ylabel "<---OMEGA2--->"'
  write ( command_unit, '(a)' ) 'set zlabel "<---U_MEAN(OMEGA1,OMEGA2)--->"'
  write ( command_unit, '(a)' ) &
    'set title "Solution Mean as Function of Omega1, Omega2"'
  write ( command_unit, '(a)' ) 'set contour'
  write ( command_unit, '(a)' ) 'set timestamp'
  write ( command_unit, '(a)' ) 'set cntrparam levels 10'
  write ( command_unit, '(a)' ) 'set view 75, 75'
  write ( command_unit, '(a)' ) 'unset key'
  write ( command_unit, '(a)' ) 'splot "' // trim ( data_filename ) // '"'

  close ( unit = command_unit )

  write ( *, '(a)' ) &
    '  Created graphics command file "' // trim ( command_filename ) // '".'
!
!  Print the maximum value of the mean.
!
  u_mean_max = maxval ( u_mean_mat )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  U_Mean_Max = ', u_mean_max
!
!  Free memory.
!
  deallocate ( omega1_mat )
  deallocate ( omega1_vec )
  deallocate ( omega2_mat )
  deallocate ( omega2_vec )
  deallocate ( umat )
  deallocate ( u_mean_mat )
  deallocate ( xmat )
  deallocate ( xvec )
  deallocate ( ymat )
  deallocate ( yvec )

  return
end
subroutine boundary ( nx, ny, x, y, n, a, rhs )

!*****************************************************************************80
!
!! BOUNDARY sets up the matrix and right hand side at boundary nodes.
!
!  Discussion:
!
!    For this simple problem, the boundary conditions specify that the solution
!    is 100 on the left side, and insulated on the right, top and bottom.
!
!    Nodes are assigned a single index K, which increases as:
!
!    (NY-1)*NX+1  (NY-1)*NX+2  ...  NY * NX
!           ....         ....  ...    .....
!           NX+1         NX+2  ...   2 * NX
!              1            2  ...       NX
!
!    The index K of a node on the lower boundary satisfies:
!      1 <= K <= NX
!    The index K of a node on the upper boundary satisfies:
!      (NY-1)*NX+1 <= K <= NY * NX
!    The index K of a node on the left boundary satisfies:
!      mod ( K, NX ) = 1
!    The index K of a node on the right boundary satisfies:
!      mod ( K, NX ) = 0
!
!    If we number rows from bottom I = 1 to top I = NY
!    and columns from left J = 1 to right J = NX, then the relationship
!    between the single index K and the row and column indices I and J is:
!      K = ( I - 1 ) * NX + J
!    and
!      J = 1 + mod ( K - 1, NX )
!      I = 1 + ( K - J ) / NX
!      
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NX, NY, the number of grid points in X and Y.
!
!    Input, real ( kind = rk8 ) X(NX), Y(NY), the coordinates of grid lines.
!
!    Input, integer N, the number of nodes.
!
!    Input/output, real ( kind = rk8 ) A(N,N).  On input, the system matrix, with the 
!    entries for the interior nodes filled in.  On output, the entries for
!    the boundary nodes have been set as well.
!
!    Input, real ( kind = rk8 ) RHS(N), on input, the system right hand side, 
!    with the entries for the interior nodes filled in.  On output, the entries for
!    the boundary nodes have been set as well.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n
  integer nx
  integer ny

  real ( kind = rk8 ) a(n,n)
  integer i
  integer j
  integer kc
  real ( kind = rk8 ) rhs(n)
  real ( kind = rk8 ) x(n)
  real ( kind = rk8 ) y(n)
!
!  Left boundary.
!
  j = 1
  do i = 2, ny - 1
    kc = ( i - 1 ) * nx + j
    a(kc,kc) = a(kc,kc) + 1.0D+00
    rhs(kc) = 0.0D+00
  end do
!
!  Right boundary.
!
  j = nx
  do i = 2, ny - 1
    kc = ( i - 1 ) * nx + j
    a(kc,kc) = a(kc,kc) + 1.0D+00
    rhs(kc) = 0.0D+00
  end do
!
!  Lower boundary.
!
  i = 1
  do j = 1, nx
    kc = ( i - 1 ) * nx + j
    a(kc,kc) = a(kc,kc) + 1.0D+00
    rhs(kc) = 0.0D+00
  end do
!
!  Upper boundary.
!
  i = ny
  do j = 1, nx
    kc = ( i - 1 ) * nx + j
    a(kc,kc) = a(kc,kc) + 1.0D+00
    rhs(kc) = 0.0D+00
  end do

  return
end
function test01_f ( x, y )

!*****************************************************************************80
!
!! TEST01_F evaluates the heat source term.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk8 ) X, Y, the evaluation point.
!
!    Output, real ( kind = rk8 ) TEST01_F, the value of the heat source term at (X,Y).
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) arg
  real ( kind = rk8 ) test01_f
  real ( kind = rk8 ) v
  real ( kind = rk8 ) x
  real ( kind = rk8 ) y

  v = 0.05D+00
  arg = ( ( x - 0.60D+00 )**2 + ( y - 0.80D+00 )**2 ) / v**2
  test01_f = 2000.0D+00 * exp ( - arg )

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
