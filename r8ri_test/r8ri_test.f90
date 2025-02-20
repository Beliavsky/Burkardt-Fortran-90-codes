program main

!*****************************************************************************80
!
!! r8ri_test() tests r8ri().
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
  write ( *, '(a)' ) 'r8ri_test():'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test r8ri().'

  call r8ri_dif2_test ( )
  call r8ri_indicator_test ( )
  call r8ri_mtv_test ( )
  call r8ri_mv_test ( )
  call r8ri_print_test ( )
  call r8ri_print_some_test ( )
  call r8ri_random_test ( )
  call r8ri_to_r8ge_test ( )
  call r8ri_zeros_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8ri_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine r8ri_dif2_test ( )

!*****************************************************************************80
!
!! R8RI_DIF2_TEST tests R8RI_DIF2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5
  integer, parameter :: nz = 3 * n - 1

  real ( kind = rk ) a(nz)
  integer ija(nz)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8RI_DIF2_TEST'
  write ( *, '(a)' ) '  R8RI_DIF2 sets up an R8RI indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
  write ( *, '(a,i8)' ) '  Storage NZ =     ', nz

  call r8ri_dif2 ( n, nz, ija, a )

  call r8ri_print ( n, nz, ija, a, '  The R8RI second difference matrix:' )

  return
end
subroutine r8ri_indicator_test ( )

!*****************************************************************************80
!
!! R8RI_INDICATOR_TEST tests R8RI_INDICATOR.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5
  integer, parameter :: nz = 11

  real ( kind = rk ) :: a(nz)
  integer :: ija(nz) = (/ &
     7,  8,  8, 10, 11, &
    12,  3,  2,  4,  5, &
     4 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8RI_INDICATOR_TEST'
  write ( *, '(a)' ) '  R8RI_INDICATOR returns an R8RI indicator matrix'
  write ( *, '(a)' ) '  for a given sparsity pattern.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
  write ( *, '(a,i8)' ) '  Storage NZ =     ', nz

  call r8ri_indicator ( n, nz, ija, a )

  call r8ri_print ( n, nz, ija, a, '  The R8RI matrix:' )

  return
end
subroutine r8ri_mtv_test ( )

!*****************************************************************************80
!
!! R8RI_MTV_TEST tests R8RI_MTV.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5
  integer, parameter :: nz = 11

  real ( kind = rk ) :: a(nz) = (/ &
    3.0, 4.0, 5.0, 0.0, 8.0, &
    0.0, 1.0, 7.0, 9.0, 2.0, &
    6.0 /)
  real ( kind = rk ) b(n)
  integer :: ija(nz) = (/ &
     7,  8,  8, 10, 11, &
    12,  3,  2,  4,  5, &
     4 /)
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8RI_MTV_TEST'
  write ( *, '(a)' ) '  R8RI_MTV computes b=A''*x, where A is an R8RI matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
  write ( *, '(a,i8)' ) '  Storage NZ =     ', nz

  call r8ri_print ( n, nz, ija, a, '  The R8RI matrix:' )

  call r8vec_indicator1 ( n, x )

  call r8vec_print ( n, x, '  The vector x:' )

  call r8ri_mtv ( n, nz, ija, a, x, b )

  call r8vec_print ( n, b, '  The product b=A''*x' )

  return
end
subroutine r8ri_mv_test ( )

!*****************************************************************************80
!
!! R8RI_MV_TEST tests R8RI_MV.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5
  integer, parameter :: nz = 11

  real ( kind = rk ) :: a(nz) = (/ &
    3.0, 4.0, 5.0, 0.0, 8.0, &
    0.0, 1.0, 7.0, 9.0, 2.0, &
    6.0 /)
  real ( kind = rk ) b(n)
  integer :: ija(nz) = (/ &
     7,  8,  8, 10, 11, &
    12,  3,  2,  4,  5, &
     4 /)
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8RI_MV_TEST'
  write ( *, '(a)' ) '  R8RI_MV computes b=A*x, where A is an R8RI matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
  write ( *, '(a,i8)' ) '  Storage NZ =     ', nz

  call r8ri_print ( n, nz, ija, a, '  The R8RI matrix:' )

  call r8vec_indicator1 ( n, x )

  call r8vec_print ( n, x, '  The vector x:' )

  call r8ri_mv ( n, nz, ija, a, x, b )

  call r8vec_print ( n, b, '  The product b=A*x' )

  return
end
subroutine r8ri_print_test ( )

!*****************************************************************************80
!
!! R8RI_PRINT_TEST tests R8RI_PRINT.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5
  integer, parameter :: nz = 11

  real ( kind = rk ) :: a(nz) = (/ &
    3.0, 4.0, 5.0, 0.0, 8.0, &
    0.0, 1.0, 7.0, 9.0, 2.0, &
    6.0 /)
  integer :: ija(nz) = (/ &
     7,  8,  8, 10, 11, &
    12,  3,  2,  4,  5, &
     4 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8RI_PRINT_TEST'
  write ( *, '(a)' ) '  R8RI_PRINT prints an R8RI matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
  write ( *, '(a,i8)' ) '  Storage NZ =     ', nz

  call r8ri_print ( n, nz, ija, a, '  The R8RI matrix:' )

  return
end
subroutine r8ri_print_some_test ( )

!*****************************************************************************80
!
!! R8RI_PRINT_SOME_TEST tests R8RI_PRINT_SOME.
!
!  Discussion:
!
!    The matrix is related to the discrete 2D Laplacian on a 3x3 grid.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 9
  integer, parameter :: nz = 34

  real ( kind = rk ) :: a(nz) = (/ &
    4.0,  4.0,  4.0,  4.0,  4.0, &
    4.0,  4.0,  4.0,  4.0,  0.0, &
   -1.0, -1.0, -1.0, -1.0, -1.0, &
   -1.0, -1.0, -1.0, -1.0, -1.0, &
   -1.0, -1.0, -1.0, -1.0, -1.0, &
   -1.0, -1.0, -1.0, -1.0, -1.0, &
   -1.0, -1.0, -1.0, -1.0 /)

  integer :: ija(nz) = (/ &
     11, 13, 16, 18, 21, &
     25, 28, 30, 33, 35, &
      2,  4,  1,  3,  5, &
      2,  6,  1,  5,  7, &
      2,  4,  6,  8,  3, &
      5,  9,  4,  8,  5, &
      7,  9,  6,  8 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8RI_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R8RI_PRINT_SOME prints some of an R8RI matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
  write ( *, '(a,i8)' ) '  Storage NZ =     ', nz

  call r8ri_print_some ( n, nz, ija, a, 1, 4, 9, 6, '  Rows 1-9, Cols 4-6:' )

  return
end
subroutine r8ri_random_test ( )

!*****************************************************************************80
!
!! R8RI_RANDOM_TEST tests R8RI_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5
  integer, parameter :: nz = 11

  real ( kind = rk ) :: a(nz)
  integer :: ija(nz) = (/ &
     7,  8,  8, 10, 11, &
    12,  3,  2,  4,  5, &
     4 /)
  integer seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8RI_RANDOM_TEST'
  write ( *, '(a)' ) '  R8RI_RANDOM randomizes an R8RI matrix'
  write ( *, '(a)' ) '  for a given sparsity pattern.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
  write ( *, '(a,i8)' ) '  Storage NZ =     ', nz

  seed = 123456789
  call r8ri_random ( n, nz, ija, seed, a )

  call r8ri_print ( n, nz, ija, a, '  The R8RI matrix:' )

  return
end
subroutine r8ri_to_r8ge_test ( )

!*****************************************************************************80
!
!! R8RI_TO_R8GE_TEST tests R8RI_TO_R8GE.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5
  integer, parameter :: nz = 11

  real ( kind = rk ) :: a(nz)
  real ( kind = rk ) a_r8ge(n,n)
  integer :: ija(nz) = (/ &
     7,  8,  8, 10, 11, &
    12,  3,  2,  4,  5, &
     4 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8RI_TO_R8GE_TEST'
  write ( *, '(a)' ) '  R8RI_TO_R8GE converts an R8RI matrix to R8GE format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
  write ( *, '(a,i8)' ) '  Storage NZ =     ', nz

  call r8ri_indicator ( n, nz, ija, a )

  call r8ri_print ( n, nz, ija, a, '  The R8RI matrix:' )

  call r8ri_to_r8ge ( n, nz, ija, a, a_r8ge )

  call r8ge_print ( n, n, a_r8ge, '  The R8GE matrix:' )

  return
end
subroutine r8ri_zeros_test ( )

!*****************************************************************************80
!
!! R8RI_ZEROS_TEST tests R8RI_ZEROS.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5
  integer, parameter :: nz = 11

  real ( kind = rk ) :: a(nz)
  integer :: ija(nz) = (/ &
     7,  8,  8, 10, 11, &
    12,  3,  2,  4,  5, &
     4 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8RI_ZEROS_TEST'
  write ( *, '(a)' ) '  R8RI_ZEROS zeros an R8RI indicator matrix'
  write ( *, '(a)' ) '  for a given sparsity pattern.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
  write ( *, '(a,i8)' ) '  Storage NZ =     ', nz

  call r8ri_zeros ( n, nz, ija, a )

  call r8ri_print ( n, nz, ija, a, '  The zero R8RI matrix:' )

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

