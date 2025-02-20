program main

!*****************************************************************************80
!
!! r83s_test() tests r83s().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 August 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r83s_test():'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test r83s().'

  call r83s_cg_test ( )
  call r83s_dif2_test ( )
  call r83s_gs_sl_test ( )
  call r83s_indicator_test ( )
  call r83s_jac_sl_test ( )
  call r83s_mtv_test ( )
  call r83s_mv_test ( )
  call r83s_print_test ( )
  call r83s_print_some_test ( )
  call r83s_random_test ( )
  call r83s_res_test ( )
  call r83s_to_r8ge_test ( )
  call r83s_zeros_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'r83s_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine r83s_cg_test ( )

!*****************************************************************************80
!
!! R83S_CG_TEST tests R83S_CG.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 July 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(3)
  real ( kind = rk ), allocatable :: b(:)
  real ( kind = rk ) e_norm
  integer n
  real ( kind = rk ), allocatable :: r(:)
  real ( kind = rk ) r_norm
  real ( kind = rk ) r8vec_norm_affine
  real ( kind = rk ) r8vec_norm
  integer seed
  real ( kind = rk ), allocatable :: x1(:)
  real ( kind = rk ), allocatable :: x2(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83S_CG_TEST'
  write ( *, '(a)' ) '  R83S_CG applies the conjugate gradient method'
  write ( *, '(a)' ) '  to solve a linear system with an R83S matrix.'

  seed = 123456789
  n = 10
!
!  Let A be the -1 2 -1 matrix.
!
  call r83s_dif2 ( n, n, a )
!
!  Choose a random solution.
!
  allocate ( x1(1:n) )
  call r8vec_uniform_01 ( n, seed, x1 )
!
!  Compute the corresponding right hand side.
!
  allocate ( b(1:n) )
  call r83s_mv ( n, n, a, x1, b )
!
!  Call the CG routine.
!
  allocate ( x2(1:n) )
  x2(1:n) = 1.0D+00
  call r83s_cg ( n, a, b, x2 )
!
!  Compute the residual.
!
  allocate ( r(1:n) )
  call r83s_res ( n, n, a, x2, b, r )
  r_norm = r8vec_norm ( n, r )
!
!  Compute the error.
!
  e_norm = r8vec_norm_affine ( n, x1, x2 )
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  Number of variables N = ', n
  write ( *, '(a,g14.6)' ) '  Norm of residual ||Ax-b|| = ', r_norm
  write ( *, '(a,g14.6)' ) '  Norm of error ||x1-x2|| = ', e_norm
!
!  Free memory.
!
  deallocate ( b )
  deallocate ( r )
  deallocate ( x1 )
  deallocate ( x2 )

  return
end
subroutine r83s_dif2_test ( )

!*****************************************************************************80
!
!! R83S_DIF2_TEST tests R83S_DIF2.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(3)
  integer i
  integer m
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83S_DIF2_TEST'
  write ( *, '(a)' ) '  R83S_DIF2 sets an R83S matrix to the second difference.'
  write ( *, '(a)' ) '  We check three cases, M<N, M=N, M>N.'

  do i = 1, 3

    if ( i == 1 ) then
      m = 3
      n = 5
    else if ( i == 2 ) then
      m = 5
      n = 5
    else if ( i == 3 ) then
      m = 5
      n = 3
    end if

    call r83s_dif2 ( m, n, a )

    call r83s_print ( m, n, a, '  Second difference in R83S format:' )

  end do

  return
end
subroutine r83s_gs_sl_test ( )

!*****************************************************************************80
!
!! R83S_GS_SL_TEST tests R83S_GS_SL.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 10

  real ( kind = rk ) a(3)
  real ( kind = rk ) b(n)
  real ( kind = rk ) diff
  integer i
  integer it
  integer :: it_max = 25
  real ( kind = rk ) :: tol = 0.000001D+00
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83S_GS_SL_TEST'
  write ( *, '(a)' ) '  R83S_GS_SL solves a linear system using'
  write ( *, '(a)' ) '  Gauss-Seidel iteration, with R83S matrix storage.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =      ', n
  write ( *, '(a,i8)' ) '  Iterations per call = ', it_max
!
!  Set the matrix values.
!
  call r83s_dif2 ( n, n, a )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r83s_mv ( n, n, a, x, b )
!
!  Set the starting solution.
!
  x(1:n) = 0.0D+00
!
!  Solve the linear system.
!
  do i = 1, 3

    call r83s_gs_sl ( n, a, b, x, tol, it_max, it, diff )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Number of iterations taken = ', it
    write ( *, '(a,g14.6)' ) '  Maximum solution change on last step = ', diff

    call r8vec_print ( n, x, '  Current solution estimate:' )

  end do

  return
end
subroutine r83s_indicator_test ( )

!*****************************************************************************80
!
!! R83S_INDICATOR_TEST tests R83S_INDICATOR.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(3)
  integer i
  integer m
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83S_INDICATOR_TEST'
  write ( *, '(a)' ) '  R83S_INDICATOR sets an R83S matrix to an indicator matrix.'
  write ( *, '(a)' ) '  We check three cases, M<N, M=N, M>N.'

  do i = 1, 3

    if ( i == 1 ) then
      m = 3
      n = 5
    else if ( i == 2 ) then
      m = 5
      n = 5
    else if ( i == 3 ) then
      m = 5
      n = 3
    end if

    call r83s_indicator ( m, n, a )

    call r83s_print ( m, n, a, '  R83S indicator matrix:' )

  end do

  return
end
subroutine r83s_jac_sl_test ( )

!*****************************************************************************80
!
!! R83S_JAC_SL_TEST tests R83S_JAC_SL.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 10

  real ( kind = rk ) a(3)
  real ( kind = rk ) b(n)
  real ( kind = rk ) diff
  integer i
  integer it
  integer :: it_max = 25
  real ( kind = rk ) :: tol = 0.000001D+00
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83S_JAC_SL_TEST'
  write ( *, '(a)' ) '  R83S_JAC_SL solves a linear system using'
  write ( *, '(a)' ) '  Jacobi iteration, with R83S matrix storage.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =      ', n
  write ( *, '(a,i8)' ) '  Iterations per call = ', it_max
!
!  Set the matrix values.
!
  call r83s_dif2 ( n, n, a )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r83s_mv ( n, n, a, x, b )
!
!  Set the starting solution.
!
  x(1:n) = 0.0D+00
!
!  Solve the linear system.
!
  do i = 1, 3

    call r83s_jac_sl ( n, a, b, x, tol, it_max, it, diff )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Number of iterations taken = ', it
    write ( *, '(a,g14.6)' ) '  Maximum solution change on last step = ', diff

    call r8vec_print ( n, x, '  Current solution estimate:' )

  end do

  return
end
subroutine r83s_mtv_test ( )

!*****************************************************************************80
!
!! R83S_MTV_TEST tests R83S_MTV.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a_83s(3)
  real ( kind = rk ), allocatable :: a_ge(:,:)
  real ( kind = rk ), allocatable :: ax_83s(:)
  real ( kind = rk ), allocatable :: ax_ge(:)
  integer i
  integer m
  integer n
  integer seed
  real ( kind = rk ), allocatable :: x(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83S_MTV_TEST'
  write ( *, '(a)' ) '  R83S_MTV computes b=A''*x, where A is an R83S matrix.'
  write ( *, '(a)' ) '  We check three cases, M<N, M=N, M>N.'

  do i = 1, 3

    if ( i == 1 ) then
      m = 3
      n = 5
    else if ( i == 2 ) then
      m = 5
      n = 5
    else if ( i == 3 ) then
      m = 5
      n = 3
    end if

    allocate ( a_ge(m,n) )
    allocate ( ax_83s(n) )
    allocate ( ax_ge(n) )
    allocate ( x(m) )

    seed = 123456789
    call r83s_random ( m, n, seed, a_83s )
    call r8vec_indicator1 ( m, x )
    call r83s_mtv ( m, n, a_83s, x, ax_83s )
    call r83s_to_r8ge ( m, n, a_83s, a_ge )
    call r8ge_mtv ( m, n, a_ge, x, ax_ge )
    call r8vec2_print ( n, ax_83s, ax_ge, '  Product comparison:' )

    deallocate ( a_ge )
    deallocate ( ax_83s )
    deallocate ( ax_ge )
    deallocate ( x )

  end do

  return
end
subroutine r83s_mv_test ( )

!*****************************************************************************80
!
!! R83S_MV_TEST tests R83S_MV.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a_83s(3)
  real ( kind = rk ), allocatable :: a_ge(:,:)
  real ( kind = rk ), allocatable :: ax_83s(:)
  real ( kind = rk ), allocatable :: ax_ge(:)
  integer i
  integer m
  integer n
  integer seed
  real ( kind = rk ), allocatable :: x(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83S_MV_TEST'
  write ( *, '(a)' ) '  R83S_MV computes b=A*x, where A is an R83S matrix.'
  write ( *, '(a)' ) '  We check three cases, M<N, M=N, M>N.'

  do i = 1, 3

    if ( i == 1 ) then
      m = 3
      n = 5
    else if ( i == 2 ) then
      m = 5
      n = 5
    else if ( i == 3 ) then
      m = 5
      n = 3
    end if

    allocate ( a_ge(m,n) )
    allocate ( ax_83s(m) )
    allocate ( ax_ge(m) )
    allocate ( x(n) )

    seed = 123456789
    call r83s_random ( m, n, seed, a_83s )
    call r8vec_indicator1 ( n, x )
    call r83s_mv ( m, n, a_83s, x, ax_83s )
    call r83s_to_r8ge ( m, n, a_83s, a_ge )
    call r8ge_mv ( m, n, a_ge, x, ax_ge )
    call r8vec2_print ( m, ax_83s, ax_ge, '  Product comparison:' )

    deallocate ( a_ge )
    deallocate ( ax_83s )
    deallocate ( ax_ge )
    deallocate ( x )

  end do

  return
end
subroutine r83s_print_test ( )

!*****************************************************************************80
!
!! R83S_PRINT_TEST tests R83S_PRINT.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(3)
  integer m
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83S_PRINT_TEST'
  write ( *, '(a)' ) '  R83S_PRINT prints an R83S matrix.'

  m = 5
  n = 4

  call r83s_indicator ( m, n, a )

  call r83s_print ( m, n, a, '  R83S indicator matrix:' )

  return
end
subroutine r83s_print_some_test ( )

!*****************************************************************************80
!
!! R83S_PRINT_SOME_TEST tests R83S_PRINT_SOME.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 5

  real ( kind = rk ) a(3)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83S_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R83S_PRINT_SOME prints some of an R83S matrix.'
!
!  Set the matrix.
!
  call r83s_indicator ( m, n, a )

  call r83s_print_some ( m, n, a, 2, 2, 5, 4, '  Rows 2-5, Cols 2-4:' )

  return
end
subroutine r83s_random_test ( )

!*****************************************************************************80
!
!! R83S_RANDOM_TEST tests R83S_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(3)
  integer m
  integer n
  integer seed

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83S_RANDOM_TEST'
  write ( *, '(a)' ) '  R83S_RANDOM randomizes an R83S matrix.'

  m = 5
  n = 4
  seed = 123456789

  call r83s_random ( m, n, seed, a )

  call r83s_print ( m, n, a, '  R83S matrix:' )

  return
end
subroutine r83s_res_test ( )

!*****************************************************************************80
!
!! R83S_RES_TEST tests R83S_RES.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(3)
  real ( kind = rk ), allocatable :: b(:)
  integer i
  integer m
  integer n
  real ( kind = rk ), allocatable :: r(:)
  integer seed
  real ( kind = rk ), allocatable :: x(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83S_RES_TEST'
  write ( *, '(a)' ) '  R83S_RES computes b-A*x, where A is an R83S matrix.'
  write ( *, '(a)' ) '  We check three cases, M<N, M=N, M>N.'

  do i = 1, 3

    if ( i == 1 ) then
      m = 3
      n = 5
    else if ( i == 2 ) then
      m = 5
      n = 5
    else if ( i == 3 ) then
      m = 5
      n = 3
    end if

    allocate ( b(1:m) )
    allocate ( r(1:m) )
    allocate ( x(1:n) )

    seed = 123456789
    call r83s_random ( m, n, seed, a )
    call r8vec_indicator1 ( n, x )
    call r83s_mv ( m, n, a, x, b )
    call r83s_res ( m, n, a, x, b, r )
    call r8vec_print ( m, r, '  Residual A*x-b:' )

    deallocate ( b )
    deallocate ( r )
    deallocate ( x )

  end do

  return
end
subroutine r83s_to_r8ge_test ( )

!*****************************************************************************80
!
!! R83S_TO_R8GE_TEST tests R83S_TO_R8GE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 4

  real ( kind = rk ) a(3)
  real ( kind = rk ) a_ge(m,n)
  integer seed

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83S_TO_R8GE_TEST'
  write ( *, '(a)' ) '  R83S_TO_R8GE converts an R83S matrix to R8GE format.'

  seed = 123456789

  call r83s_random ( m, n, seed, a )

  call r83s_print ( m, n, a, '  R83S matrix:' )

  call r83s_to_r8ge ( m, n, a, a_ge )

  call r8ge_print ( m, n, a_ge, '  R8GE matrix:' )

  return
end
subroutine r83s_zeros_test ( )

!*****************************************************************************80
!
!! R83S_ZEROS_TEST tests R83S_ZEROS.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(3)
  integer m
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83S_ZEROS_TEST'
  write ( *, '(a)' ) '  R83S_ZEROS zeros an R83S matrix.'

  m = 5
  n = 4

  call r83s_zeros ( m, n, a )

  call r83s_print ( m, n, a, '  R83S matrix:' )

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

