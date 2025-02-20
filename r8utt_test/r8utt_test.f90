program main

!*****************************************************************************80
!
!! r8utt_test() tests r8utt().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 August 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'r8utt_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test r8utt().'

  call r8utt_det_test ( )
  call r8utt_indicator_test ( )
  call r8utt_inverse_test ( )
  call r8utt_mm_test ( )
  call r8utt_mtm_test ( )
  call r8utt_mtv_test ( )
  call r8utt_mv_test ( )
  call r8utt_print_test ( )
  call r8utt_print_some_test ( )
  call r8utt_random_test ( )
  call r8utt_sl_test ( )
  call r8utt_slt_test ( )
  call r8utt_to_r8ge_test ( )
  call r8utt_zeros_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'r8utt_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine r8utt_det_test ( )

!*****************************************************************************80
!
!! R8UTT_DET_TEST tests R8UTT_DET.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 November 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a(n)
  real ( kind = rk ) det
  integer seed

  seed = 123456789

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8UTT_DET_TEST'
  write ( *, '(a)' ) '  R8UTT_DET computes the determinant of an R8UTT matrix.'
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Matrix order N = ', n

  call  r8utt_random ( n, seed, a )

  call r8utt_print ( n, a, '  The matrix:' )
!
!  Compute the determinant.
!
  call r8utt_det ( n, a, det )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  The determinant = ', det

  return
end
subroutine r8utt_indicator_test ( )

!*****************************************************************************80
!
!! R8UTT_INDICATOR_TEST tests R8UTT_INDICATOR.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 November 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8UTT_INDICATOR_TEST'
  write ( *, '(a)' ) '  R8UTT_INDICATOR sets up an indicator matrix in R8UTT format'
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Matrix order N = ', n

  call r8utt_indicator ( n, a )

  call r8utt_print ( n, a, '  The indicator matrix:' )

  return
end
subroutine r8utt_inverse_test ( )

!*****************************************************************************80
!
!! R8UTT_INVERSE_TEST tests R8UTT_INVERSE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 November 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a(n)
  real ( kind = rk ) b(n)
  real ( kind = rk ) c(n)
  integer seed

  seed = 123456789

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8UTT_INVERSE_TEST'
  write ( *, '(a)' ) '  R8UTT_INVERSE computes the inverse of an R8UTT matrix.'

  call r8utt_random ( n, seed, a )

  call r8utt_print ( n, a, '  The matrix A:' )
!
!  Compute the inverse matrix B.
!
  call r8utt_inverse ( n, a, b )

  call r8utt_print ( n, b, '  The inverse matrix B:' )
!
!  Check
!
  call r8utt_mm ( n, a, b, c )

  call r8utt_print ( n, c, '  The product A * B:' )

  return
end
subroutine r8utt_mm_test ( )

!*****************************************************************************80
!
!! R8UTT_MM_TEST tests R8UTT_MM.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 November 2015
!
!  Author:
!
!    John Burkardt
!
  integer, parameter :: n = 5

  real ( kind = rk ) a(n)
  real ( kind = rk ) a_ge(n,n)
  real ( kind = rk ) b(n)
  real ( kind = rk ) b_ge(n,n)
  real ( kind = rk ) c(n)
  real ( kind = rk ) c_ge(n,n)
  integer seed
 
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8UTT_MM_TEST'
  write ( *, '(a)' ) '  R8UTT_MM computes C = A * B for R8UTT matrices.'
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Matrix order N = ', n

  seed = 123456789
  call r8utt_random ( n, seed, a )
  call r8utt_print ( n, a, '  Factor A:' )
  call r8utt_random ( n, seed, b )
  call r8utt_print ( n, b, '  Factor B:' )
  call r8utt_mm ( n, a, b, c )
  call r8utt_print ( n, c, '  The product C = A * B' )

  call r8utt_to_r8ge ( n, a, a_ge )
  call r8utt_to_r8ge ( n, b, b_ge )
  c_ge = matmul ( a_ge(1:n,1:n), b_ge(1:n,1:n) )
  call r8ge_print ( n, n, c_ge, '  The R8GE product C:' )

  return
end
subroutine r8utt_mtm_test ( )

!*****************************************************************************80
!
!! R8UTT_MTM_TEST tests R8UTT_MTM.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 November 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a(n)
  real ( kind = rk ) a_ge(n,n)
  real ( kind = rk ) b(n)
  real ( kind = rk ) b_ge(n,n)
  real ( kind = rk ) c_ge(n,n)
  integer seed

  seed = 123456789
 
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8UTT_MTM_TEST'
  write ( *, '(a)' ) '  R8UTT_MTM computes C = A'' * B for R8UTT matrices.'
 
  call r8utt_random ( n, seed, a )
  call r8utt_print ( n, a, '  The matrix A:' )
  call r8utt_random ( n, seed, b )
  call r8utt_print ( n, b, '  The matrix B:' )

  call r8utt_mtm ( n, a, b, c_ge )
  call r8ge_print ( n, n, c_ge, '  The product C = A'' * B:' )
!
!  Compare with an R8GE calculation.
!
  call r8utt_to_r8ge ( n, a, a_ge )
  call r8utt_to_r8ge ( n, b, b_ge )
  c_ge = matmul ( transpose ( a_ge ), b_ge )
  call r8ge_print ( n, n, c_ge, '  The R8GE product C = A'' * B:' )

  return
end
subroutine r8utt_mtv_test ( )

!*****************************************************************************80
!
!! R8UTT_MTV_TEST tests R8UTT_MTV.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 November 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a(n)
  real ( kind = rk ) b(n)
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8UTT_MTV_TEST'
  write ( *, '(a)' ) '  R8UTT_MTV computes a matrix product b=A''*x for an R8UTT matrix.'
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Matrix order N = ', n

  call r8utt_indicator ( n, a )
  call r8utt_print ( n, a, '  The matrix A:' )

  call r8vec_indicator1 ( n, x )
  call r8vec_print ( n, x, '  The vector X:' )

  call r8utt_mtv ( n, a, x, b )
  call r8vec_print ( n, b, '  The vector b=A''*x:' )

  return
end
subroutine r8utt_mv_test ( )

!*****************************************************************************80
!
!! R8UTT_MV_TEST tests R8UTT_MV
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 November 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a(n)
  real ( kind = rk ) b(n)
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8UTT_MV_TEST'
  write ( *, '(a)' ) '  R8UTT_MV computes a product b=A*x for an R8UTT matrix.'
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Matrix order N = ', n

  call r8utt_indicator ( n, a )
  call r8utt_print ( n, a, '  The R8UTT matrix A:' )

  call r8vec_indicator1 ( n, x )
  call r8vec_print ( n, x, '  Vector x:' )

  call r8utt_mv ( n, a, x, b )
  call r8vec_print ( n, b, '  Vector b = A*x:' )

  return
end
subroutine r8utt_print_some_test ( )

!*****************************************************************************80
!
!! R8UTT_PRINT_SOME_TEST tests R8UTT_PRINT_SOME.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 November 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 6

  real ( kind = rk ) a(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8UTT_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R8UTT_PRINT_SOME prints some of an R8UTT matrix.'
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Matrix order N = ', n

  call r8utt_indicator ( n, a )

  call r8utt_print_some ( n, a, 1, 4, 3, 6, '  Some of the matrix:' )

  return
end
subroutine r8utt_print_test ( )

!*****************************************************************************80
!
!! R8UTT_PRINT_TEST tests R8UTT_PRINT.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 November 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8UTT_PRINT_TEST'
  write ( *, '(a)' ) '  R8UTT_PRINT prints an R8UTT matrix.'
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Matrix order N = ', n

  call r8utt_indicator ( n, a )

  call r8utt_print ( n, a, '  The matrix:' )

  return
end
subroutine r8utt_random_test ( )

!*****************************************************************************80
!
!! R8UTT_RANDOM_TEST tests R8UTT_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 November 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a(n)
  integer seed

  seed = 123456789

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8UTT_RANDOM_TEST'
  write ( *, '(a)' ) '  R8UTT_RANDOM randomizes an R8UTT matrix.'
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Matrix order N = ', n

  call r8utt_random ( n, seed, a )

  call r8utt_print ( n, a, '  Matrix A:' )

  return
end
subroutine r8utt_sl_test ( )

!*****************************************************************************80
!
!! R8UTT_SL_TEST tests R8UTT_SL.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 November 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a(n)
  real ( kind = rk ) b(n)
  integer seed
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8UTT_SL_TEST'
  write ( *, '(a)' ) '  R8UTT_SL solves a linear system A*x=b with R8UTT matrix'
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Matrix order N = ', n

  seed = 123456789
  call r8utt_random ( n, seed, a )

  call r8utt_print ( n, a, '  Matrix A:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r8utt_mv ( n, a, x, b )
  call r8vec_print ( n, b, '  Right hand side b:' )
!
!  Solve the linear system.
!
  call r8utt_sl ( n, a, b, x )

  call r8vec_print ( n, x, '  Solution x:' )

  return
end
subroutine r8utt_slt_test ( )

!*****************************************************************************80
!
!! R8UTT_SLT_TEST tests R8UTT_SLT.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 November 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a(n)
  real ( kind = rk ) b(n)
  integer seed
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8UTT_SLT_TEST'
  write ( *, '(a)' ) '  R8UTT_SLT solves a linear system A''x=b with R8UTT matrix'
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Matrix order N = ', n

  seed = 123456789
  call r8utt_random ( n, seed, a )

  call r8utt_print ( n, a, '  Matrix A:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r8utt_mtv ( n, a, x, b )

  call r8vec_print ( n, b, '  Right hand side b:' )
!
!  Solve the linear system.
!
  call r8utt_slt ( n, a, b, x )

  call r8vec_print ( n, x, '  Solution x:' )

  return
end
subroutine r8utt_to_r8ge_test ( )

!*****************************************************************************80
!
!! R8UTT_TO_R8GE_TEST tests R8UTT_TO_R8GE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 November 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a(n)
  real ( kind = rk ) a_ge(n,n)
  integer seed

  seed = 123456789

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8UTT_TO_R8GE_TEST'
  write ( *, '(a)' ) '  R8UTT_TO_R8GE converts an R8UTT matrix to R8GE format.'

  call r8utt_random ( n, seed, a )

  call r8utt_print ( n, a, '  The random R8UTT matrix:' )
 
  call r8utt_to_r8ge ( n, a, a_ge )

  call r8ge_print ( n, n, a_ge, '  The R8GE matrix:' )

  return
end
subroutine r8utt_zeros_test ( )

!*****************************************************************************80
!
!! R8UTT_ZEROS_TEST tests R8UTT_ZEROS.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 November 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8UTT_ZEROS_TEST'
  write ( *, '(a)' ) '  R8UTT_ZEROS zeros out space for an R8UTT matrix.'
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Matrix order N = ', n

  call r8utt_zeros ( n, a )

  call r8utt_print ( n, a, '  The matrix:' )

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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
