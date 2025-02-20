program main

!*****************************************************************************80
!
!! r83t_test() tests r83t().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r83t_test():'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test r83t().'

  call r83t_cg_test ( )
  call r83t_dif2_test ( )
  call r83t_gs_sl_test ( )
  call r83t_indicator_test ( )
  call r83t_jac_sl_test ( )
  call r83t_mtv_test ( )
  call r83t_mv_test ( )
  call r83t_print_test ( )
  call r83t_print_some_test ( )
  call r83t_random_test ( )
  call r83t_res_test ( )
  call r83t_to_r8ge_test ( )
  call r83t_zeros_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r83t_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine r83t_cg_test ( )

!*****************************************************************************80
!
!! R83T_CG_TEST tests R83T_CG.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 5

  real ( kind = rk ) a(m,3)
  real ( kind = rk ) b(m)
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83T_CG_TEST'
  write ( *, '(a)' ) '  R83T_CG solves an R83T linear system using'
  write ( *, '(a)' ) '  the conjugate gradient method.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n
!
!  Set the matrix.
!
  call r83t_dif2 ( m, n, a )

  call r83t_print ( m, n, a, '  The R83T matrix A:' )

  call r8vec_indicator1 ( n, x )

  call r83t_mv ( m, n, a, x, b )

  call r8vec_print ( m, b, '  The right hand side B:' )

  x(1:n) = 0.0D+00
  call r83t_cg ( n, a, b, x )

  call r8vec_print ( n, x, '  The solution X:' )

  return
end
subroutine r83t_dif2_test ( )

!*****************************************************************************80
!
!! R83T_DIF2_TEST tests R83T_DIF2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5

  real ( kind = rk ) a(m,3)
  integer n

  n = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83T_DIF2_TEST'
  write ( *, '(a)' ) '  R83T_DIF2 sets up an R83T second difference matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n
!
!  Set the matrix.
!
  call r83t_dif2 ( m, n, a )

  call r83t_print ( m, n, a, '  The R83T second difference matrix:' )

  return
end
subroutine r83t_gs_sl_test ( )

!*****************************************************************************80
!
!! R83T_GS_SL_TEST tests R83T_GS_SL.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 10

  real ( kind = rk ) a(n,3)
  real ( kind = rk ) b(n)
  integer i
  integer :: it_max = 25
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83T_GS_SL_TEST'
  write ( *, '(a)' ) '  R83T_GS_SL solves a linear system using'
  write ( *, '(a)' ) '  Gauss-Seidel iteration, with R83T matrix storage.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =      ', n
  write ( *, '(a,i8)' ) '  Iterations per call = ', it_max
!
!  Set the matrix values.
!
  call r83t_dif2 ( n, n, a )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r83t_mv ( n, n, a, x, b )
!
!  Set the starting solution.
!
  x(1:n) = 0.0D+00
!
!  Solve the linear system.
!
  do i = 1, 3

    call r83t_gs_sl ( n, a, b, x, it_max )

    call r8vec_print ( n, x, '  Current solution estimate:' )

  end do

  return
end
subroutine r83t_indicator_test ( )

!*****************************************************************************80
!
!! R83T_INDICATOR_TEST tests R83T_INDICATOR.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5

  real ( kind = rk ) a(m,3)
  integer n

  n = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83T_INDICATOR_TEST'
  write ( *, '(a)' ) '  R83T_INDICATOR sets up an R83T indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n
!
!  Set the matrix.
!
  call r83t_indicator ( m, n, a )

  call r83t_print ( m, n, a, '  The R83T indicator matrix:' )

  return
end
subroutine r83t_jac_sl_test ( )

!*****************************************************************************80
!
!! R83T_JAC_SL_TEST tests R83T_JAC_SL.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 10

  real ( kind = rk ) a(n,3)
  real ( kind = rk ) b(n)
  integer i
  integer :: it_max = 25
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83T_JAC_SL_TEST'
  write ( *, '(a)' ) '  R83T_JAC_SL solves a linear system using'
  write ( *, '(a)' ) '  Jacobi iteration, with R83T matrix storage.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =      ', n
  write ( *, '(a,i8)' ) '  Iterations per call = ', it_max
!
!  Set the matrix values.
!
  call r83t_dif2 ( n, n, a )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r83t_mv ( n, n, a, x, b )
!
!  Set the starting solution.
!
  x(1:n) = 0.0D+00
!
!  Solve the linear system.
!
  do i = 1, 3

    call r83t_jac_sl ( n, a, b, x, it_max )

    call r8vec_print ( n, x, '  Current solution estimate:' )

  end do

  return
end
subroutine r83t_mtv_test ( )

!*****************************************************************************80
!
!! R83T_MTV_TEST tests R83T_MTV.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 6

  real ( kind = rk ) a(m,3)
  real ( kind = rk ) b(n)
  real ( kind = rk ) x(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83T_MTV_TEST'
  write ( *, '(a)' ) '  R83T_MTV multiplies an R83T matrix transposed times a vector.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n
!
!  Set the matrix.
!
  call r83t_indicator ( m, n, a )

  call r83t_print ( m, n, a, '  The R83T matrix A:' )

  call r8vec_indicator1 ( m, x )

  call r8vec_print ( m, x, '  The vector x:' )

  call r83t_mtv ( m, n, a, x, b )

  call r8vec_print ( n, b, '  The product b = A''*x:' )

  return
end
subroutine r83t_mv_test ( )

!*****************************************************************************80
!
!! R83T_MV_TEST tests R83T_MV.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 6

  real ( kind = rk ) a(m,3)
  real ( kind = rk ) b(m)
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83T_MV_TEST'
  write ( *, '(a)' ) '  R83T_MV multiplies an R83T matrix times a vector.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n
!
!  Set the matrix.
!
  call r83t_indicator ( m, n, a )

  call r83t_print ( m, n, a, '  The R83T matrix A:' )

  call r8vec_indicator1 ( n, x )

  call r8vec_print ( n, x, '  The vector x:' )

  call r83t_mv ( m, n, a, x, b )

  call r8vec_print ( m, b, '  The product b = A*x:' )

  return
end
subroutine r83t_print_test ( )

!*****************************************************************************80
!
!! R83T_PRINT_TEST tests R83T_PRINT.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5

  real ( kind = rk ) a(m,3)
  integer n

  n = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83T_PRINT_TEST'
  write ( *, '(a)' ) '  R83T_PRINT prints an R83T matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n
!
!  Set the matrix.
!
  call r83t_indicator ( m, n, a )

  call r83t_print ( m, n, a, '  The R83T matrix:' )

  return
end
subroutine r83t_print_some_test ( )

!*****************************************************************************80
!
!! R83T_PRINT_SOME_TEST tests R83T_PRINT_SOME.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 9

  real ( kind = rk ) a(m,3)
  integer n

  n = 9

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83T_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R83T_PRINT_SOME prints some of an R83T matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n
!
!  Set the matrix.
!
  call r83t_indicator ( m, n, a )

  call r83t_print_some ( m, n, a, 3, 5, 6, 8, '  Rows 3:6, Cols 5:8:' )

  return
end
subroutine r83t_random_test ( )

!*****************************************************************************80
!
!! R83T_RANDOM_TEST tests R83T_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5

  real ( kind = rk ) a(m,3)
  integer n
  integer seed

  n = 5
  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83T_RANDOM_TEST'
  write ( *, '(a)' ) '  R83T_RANDOM sets up an R83T random matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n
!
!  Set the matrix.
!
  call r83t_random ( m, n, seed, a )

  call r83t_print ( m, n, a, '  The R83T random matrix:' )

  return
end
subroutine r83t_res_test ( )

!*****************************************************************************80
!
!! R83T_RES_TEST tests R83T_RES.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 5

  real ( kind = rk ) a(m,3)
  real ( kind = rk ) b(m)
  real ( kind = rk ) r(m)
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83T_RES_TEST'
  write ( *, '(a)' ) '  R83T_RES evaluates the residual given an approximate'
  write ( *, '(a)' ) '  solution of a linear system A*x=b.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n
!
!  Set the matrix.
!
  call r83t_dif2 ( m, n, a )

  call r83t_print ( m, n, a, '  The R83T matrix A:' )

  call r8vec_indicator1 ( n, x )

  call r83t_mv ( m, n, a, x, b )

  call r8vec_print ( m, b, '  The right hand side B:' )

  x(1:n) = 0.0D+00
  call r83t_cg ( n, a, b, x )

  call r8vec_print ( n, x, '  The solution X:' )

  call r83t_res ( m, n, a, x, b, r )

  call r8vec_print ( m, r, '  The residual b-A*x:' )

  return
end
subroutine r83t_to_r8ge_test ( )

!*****************************************************************************80
!
!! R83T_TO_R8GE_TEST tests R83T_TO_R8GE.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 5

  real ( kind = rk ) a_r83t(m,3)
  real ( kind = rk ) a_r8ge(m,n)
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83T_TO_R8GE_TEST'
  write ( *, '(a)' ) '  R83T_TO_R8GE converts an R83T matrix to R8GE format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n
!
!  Set the matrix.
!
  call r83t_indicator ( m, n, a_r83t )

  call r83t_print ( m, n, a_r83t, '  The R83T indicator matrix:' )

  call r83t_to_r8ge ( m, n, a_r83t, a_r8ge )

  call r8ge_print ( m, n, a_r8ge, '  The R8GE format matrix:' )

  return
end
subroutine r83t_zeros_test ( )

!*****************************************************************************80
!
!! R83T_ZEROS_TEST tests R83T_ZEROS.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5

  real ( kind = rk ) a(m,3)
  integer n

  n = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83T_ZEROS_TEST'
  write ( *, '(a)' ) '  R83T_ZEROS sets up an R83T zero matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n
!
!  Set the matrix.
!
  call r83t_zeros ( m, n, a )

  call r83t_print ( m, n, a, '  The R83T zero matrix:' )

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

