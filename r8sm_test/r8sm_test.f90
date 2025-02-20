program main

!*****************************************************************************80
!
!! r8sm_test() tests r8sm().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 August 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8sm_test():'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test rsm().'

  call r8sm_indicator_test ( )
  call r8sm_ml_test ( )
  call r8sm_mtv_test ( )
  call r8sm_mv_test ( )
  call r8sm_print_test ( )
  call r8sm_print_some_test ( )
  call r8sm_random_test ( )
  call r8sm_sl_test ( )
  call r8sm_slt_test ( )
  call r8sm_to_r8ge_test ( )
  call r8sm_zeros_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8sm_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine r8sm_indicator_test ( )

!*****************************************************************************80
!
!! R8SM_INDICATOR_TEST tests R8SM_INDICATOR.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 4

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) u(m)
  real ( kind = rk ) v(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SM_INDICATOR_TEST'
  write ( *, '(a)' ) '  R8SM_INDICATOR sets up an R8SM indicator matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n

  call r8sm_indicator ( m, n, a, u, v )

  call r8sm_print ( m, n, a, u, v, '  The R8SM indicator matrix:' )

  return
end
subroutine r8sm_ml_test ( )

!*****************************************************************************80
!
!! R8SM_ML_TEST tests R8SM_ML.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 7
  integer, parameter :: n = m

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) b(n)
  real ( kind = rk ) b2(n)
  integer info
  integer job
  integer pivot(n)
  integer :: seed = 123456789
  real ( kind = rk ) u(m)
  real ( kind = rk ) v(n)
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SM_ML_TEST'
  write ( *, '(a)' ) '  R8SM_ML computes A*x or A''*X'
  write ( *, '(a)' ) '  where A is a Sherman Morrison matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n

  do job = 0, 1
!
!  Set the matrix.
!
    call r8sm_random ( m, n, seed, a, u, v )

    call r8sm_print ( m, n, a, u, v, '  The Sherman Morrison matrix:' )
!
!  Set the desired solution.
!
    call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call r8sm_mv ( m, n, a, u, v, x, b )
    else
      call r8sm_mtv ( m, n, a, u, v, x, b )
    end if
!
!  Factor the matrix.
!
    call r8ge_fa ( n, a, pivot, info )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8SM_ML_TEST - Fatal error!'
      write ( *, '(a)' ) '  R8GE_FA declares the matrix is singular!'
      write ( *, '(a,i8)' ) '  The value of INFO is ', info
      return
    end if
!
!  Now multiply factored matrix times solution to get right hand side again.
!
    call r8sm_ml ( n, a, u, v, pivot, x, b2, job )

    if ( job == 0 ) then
      call r8vec2_print_some ( n, b, b2, 10, '  A*x and PLU*x' )
    else
      call r8vec2_print_some ( n, b, b2, 10, '  A''*x and (PLU)''*x' )
    end if

  end do

  return
end
subroutine r8sm_mtv_test ( )

!*****************************************************************************80
!
!! R8SM_MTV_TEST tests R8SM_MTV.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 4

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) b(n)
  real ( kind = rk ) u(m)
  real ( kind = rk ) v(n)
  real ( kind = rk ) x(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SM_MTV_TEST'
  write ( *, '(a)' ) '  R8SM_MTV computes A''*x=b, where A is an R8SM matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n

  call r8sm_indicator ( m, n, a, u, v )

  call r8sm_print ( m, n, a, u, v, '  The R8SM matrix:' )

  call r8vec_indicator1 ( m, x )

  call r8vec_print ( m, x, '  The vector x:' )

  call r8sm_mtv ( m, n, a, u, v, x, b )

  call r8vec_print ( n, b, '  The product b=A''*x:' )

  return
end
subroutine r8sm_mv_test ( )

!*****************************************************************************80
!
!! R8SM_MV_TEST tests R8SM_MV.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 4

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) b(m)
  real ( kind = rk ) u(m)
  real ( kind = rk ) v(n)
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SM_MV_TEST'
  write ( *, '(a)' ) '  R8SM_MV computes A*x=b, where A is an R8SM matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n

  call r8sm_indicator ( m, n, a, u, v )

  call r8sm_print ( m, n, a, u, v, '  The R8SM matrix:' )

  call r8vec_indicator1 ( n, x )

  call r8vec_print ( n, x, '  The vector x:' )

  call r8sm_mv ( m, n, a, u, v, x, b )

  call r8vec_print ( m, b, '  The product b=A*x:' )

  return
end
subroutine r8sm_print_test ( )

!*****************************************************************************80
!
!! R8SM_PRINT_TEST tests R8SM_PRINT.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 4

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) u(m)
  real ( kind = rk ) v(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SM_PRINT_TEST'
  write ( *, '(a)' ) '  R8SM_PRINT prints an R8SM matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n

  call r8sm_indicator ( m, n, a, u, v )

  call r8sm_print ( m, n, a, u, v, '  The R8SM matrix:' )

  return
end
subroutine r8sm_print_some_test ( )

!*****************************************************************************80
!
!! R8SM_PRINT_SOME_TEST tests R8SM_PRINT_SOME.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 9
  integer, parameter :: n = 9

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) u(m)
  real ( kind = rk ) v(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SM_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R8SM_PRINT_SOME prints some of an R8SM matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n

  call r8sm_indicator ( m, n, a, u, v )

  call r8sm_print_some ( m, n, a, u, v, 2, 3, 5, 7, '  Rows 2-5, Cols 3-7:' )

  return
end
subroutine r8sm_random_test ( )

!*****************************************************************************80
!
!! R8SM_RANDOM_TEST tests R8SM_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 4

  real ( kind = rk ) a(m,n)
  integer seed
  real ( kind = rk ) u(m)
  real ( kind = rk ) v(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SM_RANDOM_TEST'
  write ( *, '(a)' ) '  R8SM_RANDOM sets up a random R8SM matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n

  seed = 123456789
  call r8sm_random ( m, n, seed, a, u, v )

  call r8sm_print ( m, n, a, u, v, '  The R8SM matrix:' )

  return
end
subroutine r8sm_sl_test ( )

!*****************************************************************************80
!
!! R8SM_SL_TEST tests R8SM_SL.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = m

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) b(n)
  integer ierror
  integer info
  integer pivot(n)
  integer :: seed = 123456789
  real ( kind = rk ) u(m)
  real ( kind = rk ) v(n)
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SM_SL_TEST'
  write ( *, '(a)' ) '  R8SM_SL solves B*x=b, where B is an R8SM matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n
!
!  Set the matrix.
!
  call r8sm_random ( m, n, seed, a, u, v )

  call r8sm_print ( m, n, a, u, v, '  The Sherman-Morrison matrix A:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r8sm_mv ( m, n, a, u, v, x, b )

  call r8vec_print ( n, b, '  The right hand side vector B:' )
!
!  Factor the matrix.
!
  call r8ge_fa ( n, a, pivot, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8SM_SL_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8GE_FA declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if
!
!  Solve the linear system.
!
  call r8sm_sl ( n, a, u, v, b, ierror, pivot )
 
  call r8vec_print ( n, b, '  Solution to A * X = B:' )

  return
end
subroutine r8sm_slt_test ( )

!*****************************************************************************80
!
!! R8SM_SLT_TEST tests R8SM_SLT.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = m

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) b(n)
  integer ierror
  integer info
  integer pivot(n)
  integer :: seed = 123456789
  real ( kind = rk ) u(m)
  real ( kind = rk ) v(n)
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SM_SLT_TEST'
  write ( *, '(a)' ) '  R8SM_SLT solves B''*x=b, where B is an R8SM matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n
!
!  Set the matrix.
!
  call r8sm_random ( m, n, seed, a, u, v )

  call r8sm_print ( m, n, a, u, v, '  The Sherman-Morrison matrix A:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r8sm_mtv ( m, n, a, u, v, x, b )

  call r8vec_print ( n, b, '  The right hand side vector B:' )
!
!  Factor the matrix.
!
  call r8ge_fa ( n, a, pivot, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8SM_SLT_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8GE_FA declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if
!
!  Solve the linear system.
!
  call r8sm_slt ( n, a, u, v, b, ierror, pivot )

  call r8vec_print ( n, b, '  Solution to A'' * X = B:' )

  return
end
subroutine r8sm_to_r8ge_test ( )

!*****************************************************************************80
!
!! R8SM_TO_R8GE_TEST tests R8SM_TO_R8GE.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 4

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) a_r8ge(m,n)
  real ( kind = rk ) u(m)
  real ( kind = rk ) v(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SM_TO_R8GE_TEST'
  write ( *, '(a)' ) '  R8SM_TO_R8GE converts an R8SM matrix to R8GE format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n

  call r8sm_indicator ( m, n, a, u, v )

  call r8sm_print ( m, n, a, u, v, '  The R8SM matrix:' )

  call r8sm_to_r8ge ( m, n, a, u, v, a_r8ge )

  call r8ge_print ( m, n, a_r8ge, '  The R8GE matrix:' )

  return
end
subroutine r8sm_zeros_test ( )

!*****************************************************************************80
!
!! R8SM_ZEROS_TEST tests R8SM_ZEROS.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 4

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) u(m)
  real ( kind = rk ) v(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SM_ZEROS_TEST'
  write ( *, '(a)' ) '  R8SM_ZEROS sets up a zero R8SM matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n

  call r8sm_zeros ( m, n, a, u, v )

  call r8sm_print ( m, n, a, u, v, '  The R8SM zero matrix:' )

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

