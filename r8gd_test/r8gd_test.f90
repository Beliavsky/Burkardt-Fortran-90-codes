program main

!*****************************************************************************80
!
!! r8gd_test() tests r8gd().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 August 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8gd_test():'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test r8gd().'

  call r8gd_dif2_test ( )
  call r8gd_indicator_test ( )
  call r8gd_mtv_test ( )
  call r8gd_mv_test ( )
  call r8gd_print_test ( )
  call r8gd_print_some_test ( )
  call r8gd_random_test ( )
  call r8gd_to_r8ge_test ( )
  call r8gd_zeros_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8gd_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine r8gd_dif2_test ( )

!*****************************************************************************80
!
!! R8GD_DIF2_TEST tests R8GD_DIF2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5
  integer, parameter :: ndiag = 3

  real ( kind = rk ) a(n,ndiag)
  integer offset(ndiag)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GD_DIF2_TEST'
  write ( *, '(a)' ) '  R8GD_DIF2 sets up an R8GD second difference matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N            = ', n
  write ( *, '(a,i8)' ) '  Number of diagonals NDIAG = ', ndiag

  call r8gd_dif2 ( n, ndiag, offset, a )

  call r8gd_print ( n, ndiag, offset, a, '  The R8GS second difference matrix:' )

  return
end
subroutine r8gd_indicator_test ( )

!*****************************************************************************80
!
!! R8GD_INDICATOR_TEST tests R8GD_INDICATOR.
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

  integer, parameter :: n = 10
  integer, parameter :: ndiag = 4

  real ( kind = rk ) a(n,ndiag)
  integer offset(ndiag)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GD_INDICATOR_TEST'
  write ( *, '(a)' ) '  R8GD_INDICATOR sets up an indicator matrix'
  write ( *, '(a)' ) '  for a general diagonal matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N            = ', n
  write ( *, '(a,i8)' ) '  Number of diagonals NDIAG = ', ndiag
!
!  Set the matrix.
!
  offset(1:4) = (/ -2, 0, 1, n - 1 /)

  call i4vec_print ( ndiag, offset, '  The offset vector:' )

  call r8gd_indicator ( n, ndiag, offset, a )

  call r8gd_print ( n, ndiag, offset, a, '  The general diagonal matrix:' )

  return
end
subroutine r8gd_mtv_test ( )

!*****************************************************************************80
!
!! R8GD_MTV_TEST tests R8GD_MTV.
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

  integer, parameter :: n = 10
  integer, parameter :: ndiag = 4

  real ( kind = rk ) a(n,ndiag)
  real ( kind = rk ) b(n)
  integer offset(ndiag)
  integer :: seed = 123456789
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GD_MTV_TEST'
  write ( *, '(a)' ) '  R8GD_MTV computes A'' * x for a general diagonal matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N            = ', n
  write ( *, '(a,i8)' ) '  Number of diagonals NDIAG = ', ndiag
!
!  Set the matrix.
!
  offset(1:4) = (/ -2, 0, 1, n - 1 /)

  call r8gd_random ( n, ndiag, offset, seed, a )

  call r8gd_print ( n, ndiag, offset, a, '  The general diagonal matrix:' )

  call r8vec_indicator1 ( n, x )
  call r8vec_print ( n, x, '  x:' )

  call r8gd_mtv ( n, ndiag, offset, a, x, b )
  call r8vec_print ( n, b, '  b=A''*x:' )

  return
end
subroutine r8gd_mv_test ( )

!*****************************************************************************80
!
!! R8GD_MV_TEST tests R8GD_MV.
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

  integer, parameter :: n = 10
  integer, parameter :: ndiag = 4

  real ( kind = rk ) a(n,ndiag)
  real ( kind = rk ) b(n)
  integer offset(ndiag)
  integer :: seed = 123456789
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GD_MV_TEST'
  write ( *, '(a)' ) '  R8GD_MV computes A * x for a general diagonal matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N            = ', n
  write ( *, '(a,i8)' ) '  Number of diagonals NDIAG = ', ndiag
!
!  Set the matrix.
!
  offset(1:4) = (/ -2, 0, 1, n - 1 /)

  call r8gd_random ( n, ndiag, offset, seed, a )

  call r8gd_print ( n, ndiag, offset, a, '  The general diagonal matrix:' )

  call r8vec_indicator1 ( n, x )
  call r8vec_print ( n, x, '  x:' )

  call r8gd_mv ( n, ndiag, offset, a, x, b )
  call r8vec_print ( n, b, '  b=A*x:' )

  return
end
subroutine r8gd_print_test ( )

!*****************************************************************************80
!
!! R8GD_PRINT_TEST tests R8GD_PRINT.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 10
  integer, parameter :: ndiag = 4

  real ( kind = rk ) a(n,ndiag)
  integer offset(ndiag)
  integer :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GD_PRINT_TEST'
  write ( *, '(a)' ) '  R8GD_PRINT prints a general diagonal matrix'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N            = ', n
  write ( *, '(a,i8)' ) '  Number of diagonals NDIAG = ', ndiag
!
!  Set the matrix.
!
  offset(1:4) = (/ -2, 0, 1, n - 1 /)

  call i4vec_print ( ndiag, offset, '  The offset vector:' )

  call r8gd_random ( n, ndiag, offset, seed, a )

  call r8gd_print ( n, ndiag, offset, a, '  The general diagonal matrix:' )

  return
end
subroutine r8gd_print_some_test ( )

!*****************************************************************************80
!
!! R8GD_PRINT_SOME_TEST tests R8GD_PRINT_SOME.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 10
  integer, parameter :: ndiag = 4

  real ( kind = rk ) a(n,ndiag)
  integer offset(ndiag)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GD_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R8GD_PRINT_SOME prints some of an R8GD matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N            = ', n
  write ( *, '(a,i8)' ) '  Number of diagonals NDIAG = ', ndiag
!
!  Set the matrix.
!
  offset(1:4) = (/ -2, 0, 1, n - 1 /)

  call i4vec_print ( ndiag, offset, '  The offset vector:' )

  call r8gd_indicator ( n, ndiag, offset, a )

  call r8gd_print_some ( n, ndiag, offset, a, 3, 3, 6, 6, '  Rows 3-6, Cols 3-6:' )

  return
end
subroutine r8gd_random_test ( )

!*****************************************************************************80
!
!! R8GD_RANDOM_TEST tests R8GD_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 10
  integer, parameter :: ndiag = 4

  real ( kind = rk ) a(n,ndiag)
  integer offset(ndiag)
  integer :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GD_RANDOM_TEST'
  write ( *, '(a)' ) '  R8GD_RANDOM randomly generates a general diagonal matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N            = ', n
  write ( *, '(a,i8)' ) '  Number of diagonals NDIAG = ', ndiag
!
!  Set the matrix.
!
  offset(1:4) = (/ -2, 0, 1, n - 1 /)

  call i4vec_print ( ndiag, offset, '  The offset vector:' )

  call r8gd_random ( n, ndiag, offset, seed, a )

  call r8gd_print ( n, ndiag, offset, a, '  The general diagonal matrix:' )

  return
end
subroutine r8gd_to_r8ge_test ( )

!*****************************************************************************80
!
!! R8GD_TO_R8GE_TEST tests R8GD_TO_R8GE.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 August 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 10
  integer, parameter :: ndiag = 4

  real ( kind = rk ) a(n,ndiag)
  real ( kind = rk ) a_r8ge(n,n)
  integer offset(ndiag)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GD_TO_R8GE_TEST'
  write ( *, '(a)' ) '  R8GD_TO_R8GE converts an R8GD matrix to R8GD format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N            = ', n
  write ( *, '(a,i8)' ) '  Number of diagonals NDIAG = ', ndiag
!
!  Set the matrix.
!
  offset(1:4) = (/ -2, 0, 1, n - 1 /)

  call i4vec_print ( ndiag, offset, '  The offset vector:' )

  call r8gd_indicator ( n, ndiag, offset, a )

  call r8gd_print ( n, ndiag, offset, a, '  The R8GD matrix:' )

  call r8gd_to_r8ge ( n, ndiag, offset, a, a_r8ge )

  call r8ge_print ( n, n, a_r8ge, '  The R8GE matrix:' )

  return
end
subroutine r8gd_zeros_test ( )

!*****************************************************************************80
!
!! R8GD_ZEROS_TEST tests R8GD_ZEROS.
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
  integer, parameter :: ndiag = 4

  real ( kind = rk ) a(n,ndiag)
  integer offset(ndiag)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GD_ZEROS_TEST'
  write ( *, '(a)' ) '  R8GD_ZEROS returns a zero R8GD matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N            = ', n
  write ( *, '(a,i8)' ) '  Number of diagonals NDIAG = ', ndiag
!
!  Set the matrix.
!
  offset(1:4) = (/ -2, 0, 1, n - 1 /)

  call i4vec_print ( ndiag, offset, '  The offset vector:' )

  call r8gd_zeros ( n, ndiag, offset, a )

  call r8gd_print ( n, ndiag, offset, a, '  The zero R8GD matrix:' )

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

