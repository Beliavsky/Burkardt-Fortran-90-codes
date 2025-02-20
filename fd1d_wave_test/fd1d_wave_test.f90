program main

!*****************************************************************************80
!
!! FD1D_WAVE_TEST() tests FD1D_WAVE().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 January 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FD1D_WAVE_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the FD1D_WAVE library.'

  call fd1d_wave_test01 ( )
  call fd1d_wave_test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FD1D_WAVE_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine fd1d_wave_test01 ( )

!*****************************************************************************80
!
!! FD1D_WAVE_TEST_01 tests the FD1D finite difference wave computation.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 January 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: t_num = 41
  integer, parameter :: x_num = 16

  real ( kind = rk ) alpha
  real ( kind = rk ) c
  character ( len = 255 ) data_filename
  integer i
  real ( kind = rk ) t
  real ( kind = rk ) t_delta
  real ( kind = rk ) t_vec(t_num)
  real ( kind = rk ) t1
  real ( kind = rk ) t2
  real ( kind = rk ) u(t_num,x_num)
  external u_x1_01
  external u_x2_01
  external u_t1_01
  external ut_t1_01
  real ( kind = rk ) u1(x_num)
  real ( kind = rk ) u2(x_num)
  real ( kind = rk ) u3(x_num)
  real ( kind = rk ) x_vec(x_num)
  real ( kind = rk ) x1
  real ( kind = rk ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FD1D_WAVE_TEST01'
  write ( *, '(a)' ) '  Try the "shark" wave.'

  x1 = 0.0D+00
  x2 = 1.5D+00
  call r8vec_linspace ( x_num, x1, x2, x_vec )

  t1 = 0.0D+00
  t2 = 4.0D+00
  call r8vec_linspace ( t_num, t1, t2, t_vec )
  t_delta = ( t2 - t1 ) / real ( t_num - 1, kind = rk )

  c = 1.0D+00
  call fd1d_wave_alpha ( x_num, x1, x2, t_num, t1, t2, c, alpha )
!
!  Load the initial condition.
!
  call u_t1_01 ( x_num, x_vec, u1 )
  u(1,1:x_num) = u1(1:x_num)
!
!  Take the first step.
!
  t = t_vec(2)
  call fd1d_wave_start ( x_num, x_vec, t, t_delta, alpha, u_x1_01, u_x2_01, &
    ut_t1_01, u1, u2 )
  u(2,1:x_num) = u2(1:x_num)
!
!  Take all the other steps.
!
  do i = 3, t_num
    t = t_vec(i)
    call fd1d_wave_step ( x_num, t, alpha, u_x1_01, u_x2_01, u1, u2, u3 )
    u(i,1:x_num) = u3(1:x_num)
    u1(1:x_num) = u2(1:x_num)
    u2(1:x_num) = u3(1:x_num)
  end do
!
!  Write the solution to a file.
!
  data_filename = 'test01_data.txt'
  call r8mat_write ( data_filename, t_num, x_num, u )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Plot data written to "' // trim ( data_filename ) // '".'

  return
end
subroutine u_x1_01 ( t, u )

!*****************************************************************************80
!
!! U_X1_01 evaluates U at the boundary X1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) T, the time.
!
!    Output, real ( kind = rk ) U, the value of U(T,X1).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nd = 6

  real ( kind = rk ) t
  real ( kind = rk ), dimension ( nd ) :: td = (/ &
    0.0D+00, 0.10D+00, 0.20D+00, 0.30D+00, 0.40D+00, 0.50D+00 /)
  real ( kind = rk ) u
  real ( kind = rk ), dimension ( nd ) :: ud = (/ &
    0.0D+00, 2.0D+00, 10.0D+00, 8.0D+00, 5.0D+00, 0.0D+00 /)

  call piecewise_linear ( nd, td, ud, 1, t, u )

  return
end
subroutine u_x2_01 ( t, u )

!*****************************************************************************80
!
!! U_X2_01 evaluates U at the boundary X2.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) T, the time.
!
!    Output, real ( kind = rk ) U, the value of U(T,X2).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) t
  real ( kind = rk ) u

  call r8_fake_use ( t )

  u = 0.0D+00

  return
end
subroutine u_t1_01 ( x_num, x_vec, u )

!*****************************************************************************80
!
!! U_T1_01 evaluates U at the initial time T1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer X_NUM, the number of nodes.
!
!    Input, real ( kind = rk ) X_VEC(X_NUM), the coordinates of the nodes.
!
!    Output, real ( kind = rk ) U(X_NUM), the value of U at the initial time. 
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer x_num

  real ( kind = rk ) u(x_num)
  real ( kind = rk ) x_vec(x_num)

  call r8_fake_use ( x_vec(1) )

  u(1:x_num) = 0.0D+00

  return
end
subroutine ut_t1_01 ( x_num, x_vec, ut )

!*****************************************************************************80
!
!! UT_T1_01 evaluates dUdT at the initial time T1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer X_NUM, the number of nodes.
!
!    Input, real ( kind = rk ) X_VEC(X_NUM), the coordinates of the nodes.
!
!    Output, real ( kind = rk ) UT(X_NUM), the value of dUdT at the initial time. 
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer x_num

  real ( kind = rk ) ut(x_num)
  real ( kind = rk ) x_vec(x_num)

  call r8_fake_use ( x_vec(1) )

  ut(1:x_num) = 0.0D+00

  return
end
subroutine fd1d_wave_test02 ( )

!*****************************************************************************80
!
!! FD1D_WAVE_TEST_02 tests the FD1D finite difference wave computation.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 January 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: t_num = 41
  integer, parameter :: x_num = 16

  real ( kind = rk ) alpha
  real ( kind = rk ) c
  character ( len = 255 ) data_filename
  integer i
  real ( kind = rk ) t
  real ( kind = rk ) t_delta
  real ( kind = rk ) t_vec(t_num)
  real ( kind = rk ) t1
  real ( kind = rk ) t2
  real ( kind = rk ) u(t_num,x_num)
  real ( kind = rk ) u1(x_num)
  real ( kind = rk ) u2(x_num)
  real ( kind = rk ) u3(x_num)
  external u_x1_02
  external u_x2_02
  external u_t1_02
  external ut_t1_02
  real ( kind = rk ) x_vec(x_num)
  real ( kind = rk ) x1
  real ( kind = rk ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FD1D_WAVE_TEST02'
  write ( *, '(a)' ) '  Try a sine curve.'

  x1 = 0.0D+00
  x2 = 1.5D+00
  call r8vec_linspace ( x_num, x1, x2, x_vec )

  t1 = 0.0D+00
  t2 = 4.0D+00
  call r8vec_linspace ( t_num, t1, t2, t_vec )
  t_delta = ( t2 - t1 ) / real ( t_num - 1, kind = rk )
!
!  Changing T2 to 4.5 is enough to push the algorithm into instability.
!
!  t2 = 4.5D+00
!
  c = 1.0D+00
  call fd1d_wave_alpha ( x_num, x1, x2, t_num, t1, t2, c, alpha )
!
!  Load the initial condition.
!
  call u_t1_02 ( x_num, x_vec, u1 )
  u(1,1:x_num) = u1(1:x_num)
!
!  Take the first step.
!
  t = t_vec(2)
  call fd1d_wave_start ( x_num, x_vec, t, t_delta, alpha, u_x1_02, u_x2_02, &
    ut_t1_02, u1, u2 )
  u(2,1:x_num) = u2(1:x_num)
!
!  Take all the other steps.
!
  do i = 3, t_num
    t = t_vec(i)
    call fd1d_wave_step ( x_num, t, alpha, u_x1_02, u_x2_02, u1, u2, u3 )
    u(i,1:x_num) = u3(1:x_num)
    u1(1:x_num) = u2(1:x_num)
    u2(1:x_num) = u3(1:x_num)
  end do
!
!  Write the solution to a file.
!
  data_filename = 'test02_data.txt'
  call r8mat_write ( data_filename, t_num, x_num, u )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Plot data written to "' // trim ( data_filename ) // '".'

  return
end
subroutine u_x1_02 ( t, u )

!*****************************************************************************80
!
!! U_X1_02 evaluates U at the boundary X1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) T, the time.
!
!    Output, real ( kind = rk ) U, the value of U(T,X1).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) t
  real ( kind = rk ) u

  call r8_fake_use ( t )

  u = 0.0D+00

  return
end
subroutine u_x2_02 ( t, u )

!*****************************************************************************80
!
!! U_X2_02 evaluates U at the boundary X2.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) T, the time.
!
!    Output, real ( kind = rk ) U, the value of U(T,X2).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) t
  real ( kind = rk ) u

  call r8_fake_use ( t )

  u = 0.0D+00

  return
end
subroutine u_t1_02 ( x_num, x_vec, u )

!*****************************************************************************80
!
!! U_T1_02 evaluates U at the initial time T1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer X_NUM, the number of nodes.
!
!    Input, real ( kind = rk ) X_VEC(X_NUM), the spatial node coordinates.
!
!    Output, real ( kind = rk ) U(X_NUM), the value of U at the initial time,
!    and every node.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer x_num

  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) u(x_num)
  real ( kind = rk ) x_vec(x_num)

  u(1:x_num) = sin ( 2.0D+00 * pi * x_vec(1:x_num) )

  return
end
subroutine ut_t1_02 ( x_num, x_vec, ut )

!*****************************************************************************80
!
!! UT_T1_02 evaluates dUdT at the initial time T1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer X_NUM, the number of spatial intervals.
!
!    Input, real ( kind = rk ) X_VEC(X_NUM), the spatial node coordinates.
!
!    Output, real ( kind = rk ) UT(X_NUM), the value of dUdT at the initial time,
!    and every node.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer x_num

  real ( kind = rk ) ut(x_num)
  real ( kind = rk ) x_vec(x_num)

  call r8_fake_use ( x_vec(1) )

  ut(1:x_num) = 0.0D+00

  return
end
