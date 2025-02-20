program main

!*****************************************************************************80
!
!! r83_np_test() tests r83_np().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    25 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r83_np_test():'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test r83_np().'

  call r83_np_det_test ( )
  call r83_np_fa_test ( )
  call r83_np_fs_test ( )
  call r83_np_fss_test ( )
  call r83_np_ml_test ( )
  call r83_np_sl_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r83_np_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine r83_np_det_test ( )

!*****************************************************************************80
!
!! R83_NP_DET_TEST tests R83_NP_DET.
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

  real ( kind = rk ) a(3,n)
  real ( kind = rk ) det
  integer info

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83_NP_DET_TEST'
  write ( *, '(a)' ) '  R83_NP_DET computes the determinant of an R83 matrix'
  write ( *, '(a)' ) '  that was factored by R83_NP_FA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r83_dif2 ( n, n, a )
!
!  Factor the matrix.
!
  call r83_np_fa ( n, a, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R83_NP_DET - Warning!'
    write ( *, '(a)' ) '  R83_NP_FA returns INFO = ', info
  end if

  call r83_print ( n, n, a, '  The factored R83 matrix:' )
!
!  Compute the determinant.
!
  call r83_np_det ( n, a, det )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  R83_NP_DET computes determinant = ', det
  write ( *, '(a,g14.6)' ) '  Exact determinant = ', real ( n + 1, kind = rk )

  return
end
subroutine r83_np_fa_test ( )

!*****************************************************************************80
!
!! R83_NP_FA_TEST tests R83_NP_FA.
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

  real ( kind = rk ) a(3,n)
  real ( kind = rk ) b(n)
  integer info
  integer job
  integer :: seed = 123456789
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83_NP_FA_TEST'
  write ( *, '(a)' ) '  R83_NP_FA factors a tridiagonal matrix with no pivoting,'
  write ( *, '(a)' ) '  after which, R83_NP_SL can solve linear systems.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r83_random ( n, n, seed, a )

  call r83_print ( n, n, a, '  The tridiagonal matrix:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r83_mv ( n, n, a, x, b )
  x(1:n) = 0.0D+00
!
!  Factor the matrix.
!
  call r83_np_fa ( n, a, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R83_NP_FA_TEST - Fatal error!'
    write ( *, '(a)' ) '  The test matrix is singular.'
    return
  end if
!
!  Solve the linear system.
!
  job = 0
  call r83_np_sl ( n, a, b, job )
 
  call r8vec_print ( n, b, '  Solution:' )
!
!  Set the desired solution
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side, using the factored matrix.
!
  job = 1
  call r83_np_ml ( n, a, x, b, job )
!
!  Solve the linear system.
!
  job = 1
  call r83_np_sl ( n, a, b, job )
 
  call r8vec_print ( n, b, '  Solution to transposed system:' )
 
  return
end
subroutine r83_np_fs_test ( )

!*****************************************************************************80
!
!! R83_NP_FS_TEST tests R83_NP_FS.
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

  real ( kind = rk ) a(3,n)
  real ( kind = rk ) b(n)
  integer :: seed = 123456789
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83_NP_FS_TEST'
  write ( *, '(a)' ) '  R83_NP_FS factors and solves a tridiagonal'
  write ( *, '(a)' ) '  linear system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix elements.
!
  call r83_random ( n, n, seed, a )
  call r83_print ( n, n, a, 'What is this?' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute b = A * x.
!
  call r83_mv ( n, n, a, x, b )
!
!  Wipe out the solution.
!
  x(1:n) = 0.0D+00
!
!  Solve the system.
!
  call r83_np_fs ( n, a, b, x )

  call r8vec_print ( n, x, '  Solution:' )

  return
end
subroutine r83_np_fss_test ( )

!*****************************************************************************80
!
!! R83_NP_FSS_TEST tests R83_NP_FSS.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 10
  integer, parameter :: nb = 2

  real ( kind = rk ) a(3,n)
  real ( kind = rk ) b(n,nb)
  real ( kind = rk ) b1(n)
  real ( kind = rk ) b2(n)
  integer :: seed = 123456789
  real ( kind = rk ) x(n,nb)
  real ( kind = rk ) x1(n)
  real ( kind = rk ) x2(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83_NP_FSS_TEST'
  write ( *, '(a)' ) '  R83_NP_FSS factors a tridiagonal linear system without'
  write ( *, '(a)' ) '  pivoting, and solves multiple linear systems.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix elements.
!
  call r83_random ( n, n, seed, a )
!
!  Set the desired solutions.
!
  call r8vec_indicator0 ( n, x1 )
  call r8vec_indicator1 ( n, x2 )
!
!  Compute corresponding right hand sides.
!
  call r83_mv ( n, n, a, x1, b1 )
  call r83_mv ( n, n, a, x2, b2 )
!
!  Merge right hand sides into one array.
!
  b(1:n,1) = b1(1:n)
  b(1:n,2) = b2(1:n)
!
!  Solve the systems.
!
  call r83_np_fss ( n, a, nb, b, x )

  call r8ge_print ( n, nb, x, '  Solutions:' )

  return
end
subroutine r83_np_ml_test ( )

!*****************************************************************************80
!
!! R83_NP_ML_TEST tests R83_NP_ML.
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

  real ( kind = rk ) a(3,n)
  real ( kind = rk ) b(n)
  real ( kind = rk ) b2(n)
  integer info
  integer job
  integer :: seed = 123456789
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83_NP_ML_TEST'
  write ( *, '(a)' ) '  R83_NP_ML computes A*x or A''*x'
  write ( *, '(a)' ) '  where A has been factored by R83_FA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  do job = 0, 1
!
!  Set the matrix.
!
    call r83_random ( n, n, seed, a )
!
!  Set the desired solution.
!
    call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call r83_mv ( n, n, a, x, b )
    else
      call r83_mtv ( n, n, a, x, b )
    end if
!
!  Factor the matrix.
!
    call r83_np_fa ( n, a, info )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST08 - Fatal error!'
      write ( *, '(a)' ) '  R83_NP_FA declares the matrix is singular!'
      write ( *, '(a,i8)' ) '  The value of INFO is ', info
      return
    end if
!
!  Now multiply factored matrix times solution to get right hand side again.
!
    call r83_np_ml ( n, a, x, b2, job )

    if ( job == 0 ) then
      call r8vec2_print_some ( n, b, b2, 10, '  A*x and PLU*x:' )
    else
      call r8vec2_print_some ( n, b, b2, 10, '  A''*x and (PLU)''*x' )
    end if

  end do

  return
end
subroutine r83_np_sl_test ( )

!*****************************************************************************80
!
!! R83_NP_SL_TEST tests R83_NP_SL.
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

  real ( kind = rk ) a(3,n)
  real ( kind = rk ) b(n)
  integer info
  integer job
  integer :: seed = 123456789
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83_NP_SL_TEST'
  write ( *, '(a)' ) '  R83_NP_SL solves a linear system that has been'
  write ( *, '(a)' ) '  factored by R83_NP_FA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r83_random ( n, n, seed, a )

  call r83_print ( n, n, a, '  The tridiagonal matrix:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r83_mv ( n, n, a, x, b )
  x(1:n) = 0.0D+00
!
!  Factor the matrix.
!
  call r83_np_fa ( n, a, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R83_NP_SL_TEST - Fatal error!'
    write ( *, '(a)' ) '  The test matrix is singular.'
    return
  end if
!
!  Solve the linear system.
!
  job = 0
  call r83_np_sl ( n, a, b, job )
 
  call r8vec_print ( n, b, '  Solution:' )
!
!  Set the desired solution
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side, using the factored matrix.
!
  job = 1
  call r83_np_ml ( n, a, x, b, job )
!
!  Solve the linear system.
!
  job = 1
  call r83_np_sl ( n, a, b, job )
 
  call r8vec_print ( n, b, '  Solution to transposed system:' )
 
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

