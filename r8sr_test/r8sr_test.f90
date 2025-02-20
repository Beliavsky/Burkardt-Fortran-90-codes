program main

!*****************************************************************************80
!
!! r8sr_test() tests r8sr().
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
  write ( *, '(a)' ) 'r8sr_test():'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test r8sr().'

  call r8sr_dif2_test ( )
  call r8sr_indicator_test ( )
  call r8sr_mtv_test ( )
  call r8sr_mv_test ( )
  call r8sr_print_test ( )
  call r8sr_print_some_test ( )
  call r8sr_random_test ( )
  call r8sr_to_r8ge_test ( )
  call r8sr_zeros_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8sr_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine r8sr_dif2_test ( )

!*****************************************************************************80
!
!! R8SR_DIF2_TEST tests R8SR_DIF2.
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
  integer, parameter :: nz = 8

  integer col(nz)
  real ( kind = rk ) diag(n)
  real ( kind = rk ) off(nz)
  integer row(n+1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SR_DIF2_TEST'
  write ( *, '(a)' ) '  R8SR_DIF2 sets up an R8SR second difference matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8sr_dif2 ( n, nz, row, col, diag, off )

  call r8sr_print ( n, nz, row, col, diag, off, '  The R8SR indicator matrix:' )

  return
end
subroutine r8sr_indicator_test ( )

!*****************************************************************************80
!
!! R8SR_INDICATOR_TEST tests R8SR_INDICATOR.
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
  integer, parameter :: nz = 7

  integer, dimension ( nz ) :: col = (/ 2, 5, 5, 1, 1, 2, 3 /)
  real ( kind = rk ) diag(n)
  real ( kind = rk ) off(nz)
  integer, dimension (n+1) :: row = (/ 1, 3, 4, 5, 6, 8 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SR_INDICATOR_TEST'
  write ( *, '(a)' ) '  R8SR_INDICATOR sets up an R8SR indicator matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8sr_indicator ( n, nz, row, col, diag, off )

  call r8sr_print ( n, nz, row, col, diag, off, '  The R8SR indicator matrix:' )

  return
end
subroutine r8sr_mtv_test ( )

!*****************************************************************************80
!
!! R8SR_MTV_TEST tests R8SR_MTV.
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
  integer, parameter :: nz = 7

  real ( kind = rk ) b(n)
  integer, dimension ( nz ) :: col = (/ 2, 5, 5, 1, 1, 2, 3 /)
  real ( kind = rk ) diag(n)
  real ( kind = rk ) off(nz)
  integer, dimension (n+1) :: row = (/ 1, 3, 4, 5, 6, 8 /)
  integer :: seed = 123456789
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SR_MTV_TEST'
  write ( *, '(a)' ) '  R8SR_MTV multiplies a vector by an R8SR matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8sr_random ( n, nz, row, col, diag, off, seed )
!
!  Print the matrix.
!
  call r8sr_print ( n, nz, row, col, diag, off, '  The R8SR matrix:' )

  x(1) = 1.0D+00
  x(2:n-1) = 0.0D+00
  x(n) = -1.0D+00

  call r8vec_print ( n, x, '  The vector x:' )

  call r8sr_mtv ( n, nz, row, col, diag, off, x, b )

  call r8vec_print ( n, b, '  The product A'' * x:' )

  return
end
subroutine r8sr_mv_test ( )

!*****************************************************************************80
!
!! R8SR_MV_TEST tests R8SR_MV.
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
  integer, parameter :: nz = 7

  real ( kind = rk ) b(n)
  integer, dimension ( nz ) :: col = (/ 2, 5, 5, 1, 1, 2, 3 /)
  real ( kind = rk ) diag(n)
  real ( kind = rk ) off(nz)
  integer, dimension (n+1) :: row = (/ 1, 3, 4, 5, 6, 8 /)
  integer :: seed = 123456789
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SR_MV_TEST'
  write ( *, '(a)' ) '  R8SR_MV multiplies an R8SR matrix by a vector;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8sr_random ( n, nz, row, col, diag, off, seed )
!
!  Print the matrix.
!
  call r8sr_print ( n, nz, row, col, diag, off, '  The R8SR matrix:' )

  x(1) = 1.0D+00
  x(2:n-1) = 0.0D+00
  x(n) = -1.0D+00

  call r8vec_print ( n, x, '  The vector x:' )

  call r8sr_mv ( n, nz, row, col, diag, off, x, b )

  call r8vec_print ( n, b, '  The product A * x:' )

  return
end
subroutine r8sr_print_test ( )

!*****************************************************************************80
!
!! R8SR_PRINT_TEST tests R8SR_PRINT.
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
  integer, parameter :: nz = 7

  integer, dimension ( nz ) :: col = (/ 2, 5, 5, 1, 1, 2, 3 /)
  real ( kind = rk ) diag(n)
  real ( kind = rk ) off(nz)
  integer, dimension (n+1) :: row = (/ 1, 3, 4, 5, 6, 8 /)
  integer :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SR_PRINT_TEST'
  write ( *, '(a)' ) '  R8SR_PRINT prints an R8SR matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8sr_random ( n, nz, row, col, diag, off, seed )
!
!  Print the matrix.
!
  call r8sr_print ( n, nz, row, col, diag, off, '  The R8SR matrix:' )

  return
end
subroutine r8sr_print_some_test ( )

!*****************************************************************************80
!
!! R8SR_PRINT_SOME_TEST tests R8SR_PRINT_SOME.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5
  integer, parameter :: nz = 7

  integer, dimension ( nz ) :: col = (/ 2, 5, 5, 1, 1, 2, 3 /)
  real ( kind = rk ) diag(n)
  real ( kind = rk ) off(nz)
  integer, dimension (n+1) :: row = (/ 1, 3, 4, 5, 6, 8 /)
  integer :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SR_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R8SR_PRINT_SOME prints some of an R8SR matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8sr_random ( n, nz, row, col, diag, off, seed )
!
!  Print the matrix.
!
  call r8sr_print_some ( n, nz, row, col, diag, off, &
    1, 5, n, 5, '  Rows 1:N, column 5' )

  return
end
subroutine r8sr_random_test ( )

!*****************************************************************************80
!
!! R8SR_RANDOM_TEST tests R8SR_RANDOM.
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
  integer, parameter :: nz = 7

  integer, dimension ( nz ) :: col = (/ 2, 5, 5, 1, 1, 2, 3 /)
  real ( kind = rk ) diag(n)
  real ( kind = rk ) off(nz)
  integer, dimension (n+1) :: row = (/ 1, 3, 4, 5, 6, 8 /)
  integer :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SR_RANDOM_TEST'
  write ( *, '(a)' ) '  R8SR_RANDOM randomizes an R8SR matrix'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8sr_random ( n, nz, row, col, diag, off, seed )
!
!  Print the matrix.
!
  call r8sr_print ( n, nz, row, col, diag, off, '  The R8SR matrix:' )

  return
end
subroutine r8sr_to_r8ge_test ( )

!*****************************************************************************80
!
!! R8SR_TO_R8GE_TEST tests R8SR_TO_R8GE.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5
  integer, parameter :: nz = 7

  real ( kind = rk ) a_r8ge(n,n)
  integer, dimension ( nz ) :: col = (/ 2, 5, 5, 1, 1, 2, 3 /)
  real ( kind = rk ) diag(n)
  real ( kind = rk ) off(nz)
  integer, dimension (n+1) :: row = (/ 1, 3, 4, 5, 6, 8 /)
  integer :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SR_TO_R8GE_TEST'
  write ( *, '(a)' ) '  R8SR_TO_R8GE converts a matrix from R8SR to R8GE format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8sr_random ( n, nz, row, col, diag, off, seed )
!
!  Print the matrix.
!
  call r8sr_print ( n, nz, row, col, diag, off, '  The R8SR matrix:' )
!
!  Convert the matrix.
!
  call r8sr_to_r8ge ( n, nz, row, col, diag, off, a_r8ge )
!
!  Print the R8GE matrix.
!
  call r8ge_print ( n, n, a_r8ge, '  The R8GE matrix:' )

  return
end
subroutine r8sr_zeros_test ( )

!*****************************************************************************80
!
!! R8SR_ZEROS_TEST tests R8SR_ZEROS.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5
  integer, parameter :: nz = 7

  integer, dimension ( nz ) :: col = (/ 2, 5, 5, 1, 1, 2, 3 /)
  real ( kind = rk ) diag(n)
  real ( kind = rk ) off(nz)
  integer, dimension (n+1) :: row = (/ 1, 3, 4, 5, 6, 8 /)
  integer :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SR_ZEROS_TEST'
  write ( *, '(a)' ) '  R8SR_ZEROS zeros an R8SR matrix'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8sr_zeros ( n, nz, row, col, diag, off, seed )
!
!  Print the matrix.
!
  call r8sr_print ( n, nz, row, col, diag, off, '  The R8SR matrix:' )

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

