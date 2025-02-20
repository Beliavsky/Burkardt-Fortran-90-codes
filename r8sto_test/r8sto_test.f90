program main

!*****************************************************************************80
!
!! r8sto_test() tests r8sto().
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
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8sto_test():'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test r8sto().'

  call r8sto_dif2_test ( )
  call r8sto_indicator_test ( )
  call r8sto_inverse_test ( )
  call r8sto_mv_test ( )
  call r8sto_print_test ( )
  call r8sto_print_some_test ( )
  call r8sto_random_test ( )
  call r8sto_sl_test ( )
  call r8sto_to_r8ge_test ( )
  call r8sto_yw_sl_test ( )
  call r8sto_zeros_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8sto_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine r8sto_dif2_test ( )

!*****************************************************************************80
!
!! R8STO_DIF2_TEST tests R8STO_DIF2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 4

  real ( kind = rk ), dimension ( n ) :: a

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8STO_DIF2_TEST'
  write ( *, '(a)' ) '  R8STO_DIF2 sets the second difference as an R8STO matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8sto_dif2 ( n, a )

  call r8sto_print ( n, a, '  The R8STO matrix:' )

  return
end
subroutine r8sto_indicator_test ( )

!*****************************************************************************80
!
!! R8STO_INDICATOR_TEST tests R8STO_INDICATOR.
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

  integer, parameter :: n = 4

  real ( kind = rk ), dimension ( n ) :: a

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8STO_INDICATOR_TEST'
  write ( *, '(a)' ) '  R8STO_INDICATOR sets up an R8STO indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8sto_indicator ( n, a )

  call r8sto_print ( n, a, '  The R8STO indicator matrix:' )

  return
end
subroutine r8sto_inverse_test ( )

!*****************************************************************************80
!
!! R8STO_INVERSE_TEST tests R8STO_INVERSE.
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

  integer, parameter :: n = 3

  real ( kind = rk ), dimension ( n ) :: a = (/ 4.0D+00, 2.0D+00, 0.8D+00 /)
  real ( kind = rk ) a2(n,n)
  real ( kind = rk ) b(n,n)
  real ( kind = rk ) c(n,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8STO_INVERSE_TEST'
  write ( *, '(a)' ) '  R8STO_INVERSE computes the inverse of a positive'
  write ( *, '(a)' ) '  definite symmetric Toeplitz matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8sto_print ( n, a, '  The symmetric Toeplitz matrix A:' )

  call r8sto_inverse ( n, a, b )

  call r8ge_print ( n, n, b, '  The inverse matrix B:' )

  call r8sto_to_r8ge ( n, a, a2 )

  c(1:n,1:n) = matmul ( a2(1:n,1:n), b(1:n,1:n) )

  call r8ge_print ( n, n, c, '  The product C = A * B:' )

  return
end
subroutine r8sto_mv_test ( )

!*****************************************************************************80
!
!! R8STO_MV_TEST tests R8STO_MV.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 September 2015
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

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8STO_MV_TEST'
  write ( *, '(a)' ) '  R8STO_MV_TEST computes b=A*x, where A is an R8STO matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  seed = 123456789
  call r8sto_random ( n, seed, a )

  call r8sto_print ( n, a, '  The symmetric Toeplitz matrix:' )

  call r8vec_indicator1 ( n, x )
  call r8vec_print ( n, x, '  x:' )

  call r8sto_mv ( n, a, x, b )
  call r8vec_print ( n, b, '  b=A*x:' )

  return
end
subroutine r8sto_print_test ( )

!*****************************************************************************80
!
!! R8STO_PRINT_TEST tests R8STO_PRINT.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8STO_PRINT_TEST'
  write ( *, '(a)' ) '  R8STO_PRINT prints an R8STO matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8sto_indicator ( n, a )

  call r8sto_print ( n, a, '  The R8STO indicator matrix:' )

  return
end
subroutine r8sto_print_some_test ( )

!*****************************************************************************80
!
!! R8STO_PRINT_SOME_TEST tests R8STO_PRINT_SOME.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8STO_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R8STO_PRINT_SOME prints some of an R8STO matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8sto_indicator ( n, a )

  call r8sto_print_some ( n, a, 2, 1, 5, 3, '  Rows 2:5, Cols 1:3' )

  return
end
subroutine r8sto_random_test ( )

!*****************************************************************************80
!
!! R8STO_RANDOM_TEST tests R8STO_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 September 2015
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

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8STO_RANDOM_TEST'
  write ( *, '(a)' ) '  R8STO_RANDOM_TEST randomizes an R8STO matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  seed = 123456789
  call r8sto_random ( n, seed, a )

  call r8sto_print ( n, a, '  The symmetric Toeplitz matrix:' )

  return
end
subroutine r8sto_sl_test ( )

!*****************************************************************************80
!
!! R8STO_SL_TEST tests R8STO_SL.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 3

  real ( kind = rk ), dimension ( 0:n-1 ) :: a = (/ 1.0D+00, 0.5D+00, 0.2D+00 /)
  real ( kind = rk ), dimension ( n ) :: b = (/ 4.0D+00, -1.0D+00, 3.0D+00 /)
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8STO_SL_TEST'
  write ( *, '(a)' ) '  R8STO_SL solves a positive definite symmetric '
  write ( *, '(a)' ) '  Toeplitz system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8sto_print ( n, a, '  The symmetric Toeplitz matrix A:' )

  call r8vec_print ( n, b, '  The right hand side vector B:' )

  call r8sto_sl ( n, a, b, x )

  call r8vec_print ( n, x, '  The solution X:' )

  call r8sto_mv ( n, a, x, b )

  call r8vec_print ( n, b, '  The product vector  B = A * X:' )

  return
end
subroutine r8sto_to_r8ge_test ( )

!*****************************************************************************80
!
!! R8STO_TO_R8GE_TEST tests R8STO_TO_R8GE.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a_r8ge(n,n)
  real ( kind = rk ) a_r8sto(n)
  integer :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8STO_TO_R8GE_TEST'
  write ( *, '(a)' ) '  R8STO_TO_R8GE converts a matrix from R8STO to R8GE format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8sto_random ( n, seed, a_r8sto )

  call r8sto_print ( n, a_r8sto, '  The matrix in R8STO format:' )

  call r8sto_to_r8ge ( n, a_r8sto, a_r8ge )

  call r8ge_print ( n, n, a_r8ge, '  The matrix in R8GE format:' )
 
  return
end
subroutine r8sto_yw_sl_test ( )

!*****************************************************************************80
!
!! R8STO_YW_SL_TEST tests R8STO_YW_SL.
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

  integer, parameter :: n = 3

  real ( kind = rk ) a(n)
  real ( kind = rk ) b(n)
  real ( kind = rk ), dimension ( 0:n ) :: r = (/ &
    1.0D+00, 0.5D+00, 0.2D+00, 0.1D+00 /)
  real ( kind = rk ) x(n)

  a(1:n) = r(0:n-1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8STO_YW_SL_TEST'
  write ( *, '(a)' ) '  R8STO_YW_SL solves the Yule-Walker equations for a'
  write ( *, '(a)' ) '  symmetric Toeplitz system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8sto_print ( n, a, '  The symmetric Toeplitz matrix:' )

  b(1:n) = -r(1:n)
  call r8vec_print ( n, b, '  The right hand side, B:' )

  b(1:n) = -b(1:n)
  call r8sto_yw_sl ( n, b, x )

  call r8vec_print ( n, x, '  The computed solution, X:' )

  call r8sto_mv ( n, a, x, b )

  call r8vec_print ( n, b, '  The product A*X:' )

  return
end
subroutine r8sto_zeros_test ( )

!*****************************************************************************80
!
!! R8STO_ZEROS_TEST tests R8STO_ZEROS.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 4

  real ( kind = rk ), dimension ( n ) :: a

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8STO_ZEROS_TEST'
  write ( *, '(a)' ) '  R8STO_ZEROS zeros an R8STO matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8sto_zeros ( n, a )

  call r8sto_print ( n, a, '  The R8STO matrix:' )

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

