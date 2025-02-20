program main

!*****************************************************************************80
!
!! r83p_test() tests r83p().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r83p_test():'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test r83p().'

  call r83p_det_test ( )
  call r83p_fa_test ( )
  call r83p_indicator_test ( )
  call r83p_ml_test ( )
  call r83p_mtv_test ( )
  call r83p_mv_test ( )
  call r83p_print_test ( )
  call r83p_print_some_test ( )
  call r83p_random_test ( )
  call r83p_sl_test ( )
  call r83p_to_r8ge_test ( )
  call r83p_zeros_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r83p_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine r83p_det_test ( )

!*****************************************************************************80
!
!! R83P_DET_TEST tests R83P_DET.
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

  integer, parameter :: n = 12

  real ( kind = rk ) a(3,n)
  real ( kind = rk ) b(n,n)
  real ( kind = rk ) det
  integer info
  integer pivot(n)
  integer :: seed = 123456789
  real ( kind = rk ) work2(n-1)
  real ( kind = rk ) work3(n-1)
  real ( kind = rk ) work4

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83P_DET_TEST'
  write ( *, '(a)' ) '  R83P_DET, determinant of a tridiagonal'
  write ( *, '(a)' ) '  periodic matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r83p_random ( n, seed, a )

  call r83p_print ( n, a, '  The periodic tridiagonal matrix:' )
!
!  Copy the matrix into a general array.
!
  call r83p_to_r8ge ( n, a, b )
!
!  Factor the matrix.
!
  call r83p_fa ( n, a, info, work2, work3, work4 )
!
!  Compute the determinant.
!
  call r83p_det ( n, a, work4, det )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  R83P_DET computes the determinant = ', det
!
!  Factor the general matrix.
!
  call r8ge_fa ( n, b, pivot, info )
!
!  Compute the determinant.
!
  call r8ge_det ( n, b, pivot, det )

  write ( *, '(a,g14.6)' ) '  R8GE_DET computes the determinant = ', det

  return
end
subroutine r83p_fa_test ( )

!*****************************************************************************80
!
!! R83P_FA_TEST tests R83P_FA.
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
  real ( kind = rk ) work2(n-1)
  real ( kind = rk ) work3(n-1)
  real ( kind = rk ) work4
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83P_FA_TEST'
  write ( *, '(a)' ) '  R83P_FA factors a tridiagonal periodic system'
  write ( *, '(a)' ) '  which can then be solved by R83P_SL.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r83p_random ( n, seed, a )
!
!  Factor the matrix.
!
  call r83p_fa ( n, a, info, work2, work3, work4 )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Fatal error!'
    write ( *, '(a,i8)' ) '  R83P_FA returns INFO = ', info
    return
  end if

  do job = 0, 1
!
!  Set the desired solution.
!
    call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
    call r83p_ml ( n, a, x, b, job )
!
!  Solve the linear system.
!
    call r83p_sl ( n, a, b, x, job, work2, work3, work4 )

    if ( job == 0 ) then
      call r8vec_print ( n, x, '  Solution:' )
    else
      call r8vec_print ( n, x, '  Solution to transposed system:' )
    end if

  end do
 
  return
end
subroutine r83p_indicator_test ( )

!*****************************************************************************80
!
!! R83P_INDICATOR_TEST tests R83P_INDICATOR.
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

  integer, parameter :: n = 5

  real ( kind = rk ) a(3,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83P_INDICATOR_TEST'
  write ( *, '(a)' ) '  R83P_INDICATOR sets up an R83P indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r83p_indicator ( n, a )

  call r83p_print ( n, a, '  The R83P indicator matrix:' )

  return
end
subroutine r83p_ml_test ( )

!*****************************************************************************80
!
!! R83P_ML_TEST tests R83P_ML.
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
  real ( kind = rk ) work2(n-1)
  real ( kind = rk ) work3(n-1)
  real ( kind = rk ) work4
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83P_ML_TEST'
  write ( *, '(a)' ) '  R83P_ML computes A*x or A''*X'
  write ( *, '(a)' ) '  where A has been factored by R83P_FA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  do job = 0, 1
!
!  Set the matrix.
!
    call r83p_random ( n, seed, a )
!
!  Set the desired solution.
!
    call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call r83p_mv ( n, a, x, b )
    else
      call r83p_mtv ( n, a, x, b )
    end if
!
!  Factor the matrix.
!
    call r83p_fa ( n, a, info, work2, work3, work4 )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST11 - Fatal error!'
      write ( *, '(a)' ) '  R83P_FA declares the matrix is singular!'
      write ( *, '(a,i8)' ) '  The value of INFO is ', info
      return
    end if
!
!  Now multiply factored matrix times solution to get right hand side again.
!
    call r83p_ml ( n, a, x, b2, job )

    if ( job == 0 ) then
      call r8vec2_print_some ( n, b, b2, 10, '  A*x and PLU*x' )
    else
      call r8vec2_print_some ( n, b, b2, 10, '  A''*x and (PLU)''*x' )
    end if

  end do

  return
end
subroutine r83p_mtv_test ( )

!*****************************************************************************80
!
!! R83P_MTV_TEST tests R83P_MTV.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a(3,n)
  real ( kind = rk ) b(n)
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83P_MTV_TEST'
  write ( *, '(a)' ) '  R83P_MTV computes A''*x = b for an R83P matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r83p_indicator ( n, a )

  call r83p_print ( n, a, '  The R83P indicator matrix:' )

  call r8vec_indicator1 ( n, x )

  call r8vec_print ( n, x, '  Vector x:' )

  call r83p_mtv ( n, a, x, b )

  call r8vec_print ( n, b, '  b=A''*x:' )
  
  return
end
subroutine r83p_mv_test ( )

!*****************************************************************************80
!
!! R83P_MV_TEST tests R83P_MV.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a(3,n)
  real ( kind = rk ) b(n)
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83P_MV_TEST'
  write ( *, '(a)' ) '  R83P_MV computes A*x = b for an R83P matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r83p_indicator ( n, a )

  call r83p_print ( n, a, '  The R83P indicator matrix:' )

  call r8vec_indicator1 ( n, x )

  call r8vec_print ( n, x, '  Vector x:' )

  call r83p_mv ( n, a, x, b )

  call r8vec_print ( n, b, '  b=A*x:' )
  
  return
end
subroutine r83p_print_test ( )

!*****************************************************************************80
!
!! R83P_PRINT_TEST tests R83P_PRINT.
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

  integer, parameter :: n = 5

  real ( kind = rk ) a(3,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83P_PRINT_TEST'
  write ( *, '(a)' ) '  R83P_PRINT prints an R83P matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r83p_indicator ( n, a )

  call r83p_print ( n, a, '  An R83P matrix:' )

  return
end
subroutine r83p_print_some_test ( )

!*****************************************************************************80
!
!! R83P_PRINT_SOME_TEST tests R83P_PRINT_SOME.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 10

  real ( kind = rk ) a(3,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83P_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R83P_PRINT_SOME prints some of an R83P matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r83p_indicator ( n, a )

  call r83p_print_some ( n, a, 1, 1, n, 2, '  Rows 1:N, Cols 1:2:' )

  return
end
subroutine r83p_random_test ( )

!*****************************************************************************80
!
!! R83P_RANDOM_TEST tests R83P_RANDOM.
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

  integer, parameter :: n = 5

  real ( kind = rk ) a(3,n)
  integer seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83P_RANDOM_TEST'
  write ( *, '(a)' ) '  R83P_RANDOM returns a random R83P matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
  seed = 123456789
!
!  Set the matrix.
!
  call r83p_random ( n, seed, a )

  call r83p_print ( n, a, '  A random R83P matrix:' )

  return
end
subroutine r83p_sl_test ( )

!*****************************************************************************80
!
!! R83P_SL_TEST tests R83P_SL.
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
  real ( kind = rk ) work2(n-1)
  real ( kind = rk ) work3(n-1)
  real ( kind = rk ) work4
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83P_SL_TEST'
  write ( *, '(a)' ) '  R83P_SL solves a tridiagonal periodic system'
  write ( *, '(a)' ) '  after it has been factored bu R83P_FA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r83p_random ( n, seed, a )
!
!  Factor the matrix.
!
  call r83p_fa ( n, a, info, work2, work3, work4 )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Fatal error!'
    write ( *, '(a,i8)' ) '  R83P_FA returns INFO = ', info
    return
  end if

  do job = 0, 1
!
!  Set the desired solution.
!
    call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
    call r83p_ml ( n, a, x, b, job )
!
!  Solve the linear system.
!
    call r83p_sl ( n, a, b, x, job, work2, work3, work4 )

    if ( job == 0 ) then
      call r8vec_print ( n, x, '  Solution:' )
    else
      call r8vec_print ( n, x, '  Solution to transposed system:' )
    end if

  end do
 
  return
end
subroutine r83p_to_r8ge_test ( )

!*****************************************************************************80
!
!! R83P_TO_R8GE_TEST tests R83P_TO_R8GE.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a_r83p(3,n)
  real ( kind = rk ) a_r8ge(n,n)
  integer seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83P_TO_R8GE_TEST'
  write ( *, '(a)' ) '  R83P_TO_R8GE converts an R83P matrix to R8GE format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
  seed = 123456789
!
!  Set the matrix.
!
  call r83p_random ( n, seed, a_r83p )

  call r83p_print ( n, a_r83p, '  The R83P matrix:' )

  call r83p_to_r8ge ( n, a_r83p, a_r8ge )

  call r8ge_print ( n, n, a_r8ge, '  The R8GE matrix:' )

  return
end
subroutine r83p_zeros_test ( )

!*****************************************************************************80
!
!! R83P_ZEROS_TEST tests R83P_ZEROS.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a(3,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R83P_ZEROS_TEST'
  write ( *, '(a)' ) '  R83P_ZEROS returns a zero R83P matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r83p_zeros ( n, a )

  call r83p_print ( n, a, '  A zero R83P matrix:' )

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

