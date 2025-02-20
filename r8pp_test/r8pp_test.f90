program main

!*****************************************************************************80
!
!! r8pp_test() tests r8pp().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 August 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8pp_test():'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test r8pp().'

  call r8pp_det_test ( )
  call r8pp_dif2_test ( )
  call r8pp_fa_test ( )
  call r8pp_indicator_test ( )
  call r8pp_mv_test ( )
  call r8pp_print_test ( )
  call r8pp_print_some_test ( )
  call r8pp_random_test ( )
  call r8pp_sl_test ( )
  call r8pp_to_r8ge_test ( )
  call r8pp_zeros_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8pp_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine r8pp_det_test ( )

!*****************************************************************************80
!
!! R8PP_DET_TEST tests R8PP_DET.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    13 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a((n*(n+1))/2)
  real ( kind = rk ) det
  integer info

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PP_DET_TEST'
  write ( *, '(a)' ) '  R8PP_DET computes the determinant of an R8PP matrix'
  write ( *, '(a)' ) '  factored by R8PP_FA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8pp_dif2 ( n, a )

  call r8pp_print ( n, a, '  The R8PP matrix:' )
!
!  Factor the matrix.
!
  call r8pp_fa ( n, a, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'R8PP_DET_TEST - Warning!'
    write ( *, '(a)' ) '  R8PP_FA failed to factor the matrix.'
    return
  end if
!
!  Compute the determinant.
!
  call r8pp_det ( n, a, det )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Determinant = ', det
  write ( *, '(a,g14.6)' ) '  Exact determinant = ', real ( n + 1, kind = rk )

  return
end
subroutine r8pp_dif2_test ( )

!*****************************************************************************80
!
!! R8PP_DIF2_TEST tests R8PP_DIF2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    13 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a((n*(n+1))/2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PP_DIF2_TEST'
  write ( *, '(a)' ) '  R8PP_DIF2 sets up an R8PP second difference matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8pp_dif2 ( n, a )

  call r8pp_print ( n, a, '  The R8PP second difference matrix:' )
 
  return
end
subroutine r8pp_fa_test ( )

!*****************************************************************************80
!
!! R8PP_FA_TEST tests R8PP_FA.
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

  integer, parameter :: n = 5

  real ( kind = rk ) a((n*(n+1))/2)
  real ( kind = rk ) b(n)
  integer info
  integer :: seed = 123456789
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PP_FA_TEST'
  write ( *, '(a)' ) '  R8PP_FA factors an R8PP system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8pp_random ( n, seed, a )

  call r8pp_print ( n, a, '  The R8PP matrix:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )

  call r8vec_print ( n, x, '  The desired solution:' )
!
!  Compute the corresponding right hand side.
!
  call r8pp_mv ( n, a, x, b )

  call r8vec_print ( n, b, '  The right hand side:' )
!
!  Factor the matrix.
!
  call r8pp_fa ( n, a, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8PP_FA_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8PP_FA declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The R8PP matrix has been factored.'
  end if
!
!  Solve the linear system.
!
  call r8pp_sl ( n, a, b )
 
  call r8vec_print ( n, b, '  Solution:' )

  return
end
subroutine r8pp_indicator_test ( )

!*****************************************************************************80
!
!! R8PP_INDICATOR_TEST tests R8PP_INDICATOR.
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

  integer, parameter :: n = 5

  real ( kind = rk ) a((n*(n+1))/2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PP_INDICATOR_TEST'
  write ( *, '(a)' ) '  R8PP_INDICATOR sets up an R8PP indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8pp_indicator ( n, a )

  call r8pp_print ( n, a, '  The R8PP indicator matrix:' )
 
  return
end
subroutine r8pp_mv_test ( )

!*****************************************************************************80
!
!! R8PP_MV_TEST tests R8PP_MV.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a((n*(n+1))/2)
  real ( kind = rk ) b(n)
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PP_MV_TEST'
  write ( *, '(a)' ) '  R8PP_MV computes b=A*x, where A is an R8PP matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8pp_indicator ( n, a )
  call r8pp_print ( n, a, '  The R8PP indicator matrix:' )

  call r8vec_indicator1 ( n, x )
  call r8vec_print ( n, x, '  Vector x:' )

  call r8pp_mv ( n, a, x, b )
  call r8vec_print ( n, b, '  Result b=A*x:' )
 
  return
end
subroutine r8pp_print_test ( )

!*****************************************************************************80
!
!! R8PP_PRINT_TEST tests R8PP_PRINT.
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

  integer, parameter :: n = 5

  real ( kind = rk ) a((n*(n+1))/2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PP_PRINT_TEST'
  write ( *, '(a)' ) '  R8PP_PRINT prints an R8PP matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8pp_indicator ( n, a )

  call r8pp_print ( n, a, '  The R8PP matrix:' )
 
  return
end
subroutine r8pp_print_some_test ( )

!*****************************************************************************80
!
!! R8PP_PRINT_SOME_TEST tests R8PP_PRINT_SOME.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    13 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 10

  real ( kind = rk ) a((n*(n+1))/2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PP_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R8PP_PRINT_SOME prints some of an R8PP matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8pp_indicator ( n, a )

  call r8pp_print_some ( n, a, 2, 3, 6, 5, '  Rows 2-6, Cols 3-5:' )
 
  return
end
subroutine r8pp_random_test ( )

!*****************************************************************************80
!
!! R8PP_RANDOM_TEST tests R8PP_RANDOM.
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

  integer, parameter :: n = 5

  real ( kind = rk ) a((n*(n+1))/2)
  integer :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PP_RANDOM_TEST'
  write ( *, '(a)' ) '  R8PP_RANDOM, compute a random positive definite'
  write ( *, '(a)' ) '  symmetric packed matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8pp_random ( n, seed, a )

  call r8pp_print ( n, a, '  The random R8PP matrix:' )
 
  return
end
subroutine r8pp_sl_test ( )

!*****************************************************************************80
!
!! R8PP_SL_TEST tests R8PP_SL.
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

  integer, parameter :: n = 5

  real ( kind = rk ) a((n*(n+1))/2)
  real ( kind = rk ) b(n)
  integer info
  integer :: seed = 123456789
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PP_SL_TEST'
  write ( *, '(a)' ) '  R8PP_SL solves a linear system factored by R8PP_FA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8pp_random ( n, seed, a )

  call r8pp_print ( n, a, '  The R8PP matrix:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )

  call r8vec_print ( n, x, '  The desired solution:' )
!
!  Compute the corresponding right hand side.
!
  call r8pp_mv ( n, a, x, b )

  call r8vec_print ( n, b, '  The right hand side:' )
!
!  Factor the matrix.
!
  call r8pp_fa ( n, a, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8PP_SL_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8PP_FA declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The R8PP matrix has been factored.'
  end if
!
!  Solve the linear system.
!
  call r8pp_sl ( n, a, b )
 
  call r8vec_print ( n, b, '  Solution:' )

  return
end
subroutine r8pp_to_r8ge_test ( )

!*****************************************************************************80
!
!! R8PP_TO_R8GE_TEST tests R8PP_TO_R8GE.
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

  integer, parameter :: n = 5

  real ( kind = rk ) a((n*(n+1))/2)
  real ( kind = rk ) a_r8ge(n,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PP_TO_R8GE_TEST'
  write ( *, '(a)' ) '  R8PP_TO_R8GE converts an R8PP matrix to R8GE format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8pp_indicator ( n, a )

  call r8pp_print ( n, a, '  The R8PP matrix:' )

  call r8pp_to_r8ge ( n, a, a_r8ge )

  call r8ge_print ( n, n, a_r8ge, '  The R8GE matrix:' )

  return
end
subroutine r8pp_zeros_test ( )

!*****************************************************************************80
!
!! R8PP_ZEROS_TEST tests R8PP_ZEROS.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a((n*(n+1))/2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PP_ZEROS_TEST'
  write ( *, '(a)' ) '  R8PP_ZEROS sets a zero R8PP matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8pp_zeros ( n, a )

  call r8pp_print ( n, a, '  The R8PP matrix:' )
 
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

