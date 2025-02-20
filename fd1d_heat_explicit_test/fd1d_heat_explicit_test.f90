program main

!*****************************************************************************80
!
!! fd1d_heat_explicit_test() tests fd1d_heat_explicit().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 January 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'fd1d_heat_explicit_test():'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) '  Test fd1d_heat_explicit().'

  call fd1d_heat_explicit_test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'fd1d_heat_explicit_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine fd1d_heat_explicit_test01 ( )

!*****************************************************************************80
!
!! fd1d_heat_explicit_test01() does a simple test problem
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 January 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  external bc_test01
  real ( kind = rk ) cfl
  real ( kind = rk ) dt
  real ( kind = rk ), allocatable :: h(:)
  real ( kind = rk ), allocatable :: h_new(:)
  real ( kind = rk ), allocatable :: hmat(:,:)
  integer j
  real ( kind = rk ) k
  external rhs_test01
  real ( kind = rk ), allocatable :: t(:)
  real ( kind = rk ) t_max
  real ( kind = rk ) t_min
  integer t_num
  real ( kind = rk ), allocatable :: x(:)
  real ( kind = rk ) x_max
  real ( kind = rk ) x_min
  integer x_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'fd1d_heat_explicit_test01():'
  write ( *, '(a)' ) '  Compute an approximate solution to the time-dependent'
  write ( *, '(a)' ) '  one dimensional heat equation:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    dH/dt - K * d2H/dx2 = f(x,t)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Run a simple test case.'
!
!  Heat coefficient.
!
  k = 0.002D+00
!
!  X_NUM is the number of equally spaced nodes to use between 0 and 1.
!
  x_num = 21
  x_min = 0.0D+00
  x_max = 1.0D+00
  allocate ( x(1:x_num) )
  call r8vec_linspace ( x_num, x_min, x_max, x )
!
!  T_NUM is the number of equally spaced time points between 0 and 10.0.
!
  t_num = 201
  t_min = 0.0D+00
  t_max = 80.0D+00
  dt = ( t_max - t_min ) / real ( t_num - 1, kind = rk )
  allocate ( t(1:t_num) )
  call r8vec_linspace ( t_num, t_min, t_max, t )
!
!  Get the CFL coefficient.
!
  call fd1d_heat_explicit_cfl ( k, t_num, t_min, t_max, x_num, x_min, x_max, cfl )
!
!  Running the code produces an array H of temperatures H(t,x),
!  and vectors x and t.
!
  allocate ( h(1:x_num) )
  allocate ( h_new(1:x_num) )

  call ic_test01 ( x_num, x, t(1), h )
  call bc_test01 ( x_num, x, t(1), h )

  allocate ( hmat(1:x_num,1:t_num) )
  hmat(1:x_num,1) = h(1:x_num)

  do j = 2, t_num
    call fd1d_heat_explicit ( x_num, x, t(j-1), dt, cfl, rhs_test01, &
      bc_test01, h, h_new )
    hmat(1:x_num,j) = h_new(1:x_num)
    h(1:x_num) = h_new(1:x_num)
  end do
!
!  Write the data to files.
!
  call r8mat_write ( 'h_test01.txt', x_num, t_num, hmat )
  call r8vec_write ( 't_test01.txt', t_num, t )
  call r8vec_write ( 'x_test01.txt', x_num, x )

  deallocate ( h )
  deallocate ( h_new )
  deallocate ( hmat )
  deallocate ( t )
  deallocate ( x )

  return
end
subroutine bc_test01 ( x_num, x, t, h )

!*****************************************************************************80
!
!! bc_test01() evaluates the boundary conditions for problem 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 January 2022
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer X_NUM, the number of nodes.
!
!    real ( kind = rk ) X(X_NUM), the node coordinates.
!
!    real ( kind = rk ) T, the current time.
!
!    real ( kind = rk ) H(X_NUM), the current heat values.
!
!  Output:
!
!    real ( kind = rk ) H(X_NUM), the current heat values, after boundary
!    conditions have been imposed.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer x_num

  real ( kind = rk ) h(x_num)
  real ( kind = rk ) t
  real ( kind = rk ) x(x_num)

  call r8_fake_use ( t )
  call r8_fake_use ( x(1) )

  h(1)  = 90.0D+00
  h(x_num) = 70.0D+00

  return
end
subroutine ic_test01 ( x_num, x, t, h  )

!*****************************************************************************80
!
!! ic_test01() evaluates the initial condition for problem 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 January 2022
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer X_NUM, the number of nodes.
!
!    real ( kind = rk ) X(X_NUM), the node coordinates.
!
!    real ( kind = rk ) T, the initial time.
!
!  Output:
!
!    real ( kind = rk ) H(X_NUM), the heat values at the initial time.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer x_num

  real ( kind = rk ) h(x_num)
  real ( kind = rk ) t
  real ( kind = rk ) x(x_num)

  call r8_fake_use ( t )
  call r8_fake_use ( x(1) )

  h(1:x_num) = 50.0D+00

  return
end
subroutine rhs_test01 ( x_num, x, t, value )

!*****************************************************************************80
!
!! rhs_test01() evaluates the right hand side for problem 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 January 2022
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer X_NUM, the number of nodes.
!
!    real ( kind = rk ) X(X_NUM), the node coordinates.
!
!    real ( kind = rk ) T, the current time.
!
!  Output:
!
!    real ( kind = rk ) VALUE(X_NUM), the source term.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer x_num

  real ( kind = rk ) t
  real ( kind = rk ) value(x_num)
  real ( kind = rk ) x(x_num)

  call r8_fake_use ( t )
  call r8_fake_use ( x(1) )

  value(1:x_num) = 0.0D+00

  return
end
subroutine r8_fake_use ( x )

!*****************************************************************************80
!
!! r8_fake_use() pretends to use a variable.
!
!  Discussion:
!
!    Some compilers will issue a warning if a variable is unused.
!    Sometimes there's a good reason to include a variable in a program,
!    but not to use it.  Calling this function with that variable as
!    the argument will shut the compiler up.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 January 2022
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) X, the variable to be "used".
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) x

  if ( x /= x ) then
    write ( *, '(a)' ) '  r8_fake_use: variable is NAN.'
  end if

  return
end
 
