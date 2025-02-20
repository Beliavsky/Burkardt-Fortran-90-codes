program main

!*****************************************************************************80
!
!! r8ncf_test() tests r8ncf().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    24 August 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8ncf_test():'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test r8ncf().'

  call r8ncf_dif2_test ( )
  call r8ncf_indicator_test ( )
  call r8ncf_mtv_test ( )
  call r8ncf_mv_test ( )
  call r8ncf_print_test ( )
  call r8ncf_print_some_test ( )
  call r8ncf_random_test ( )
  call r8ncf_to_r8ge_test ( )
  call r8ncf_zeros_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8ncf_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine r8ncf_dif2_test ( )

!*****************************************************************************80
!
!! R8NCF_DIF2_TEST tests R8NCF_DIF2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 July 2016
!
!  Author:
!
!   John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 7

  real ( kind = rk ), allocatable :: a(:)
  integer nz_num
  integer, allocatable :: rowcol(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8NCF_DIF2_TEST'
  write ( *, '(a)' ) '  R8NCF_DIF2 sets up an R8NCF second difference matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n

  call r8ncf_dif2_nz_num ( m, n, nz_num )

  write ( *, '(a,i8)' ) '  Matrix nonzeros NZ_NUM = ', nz_num

  allocate ( rowcol(1:2,1:nz_num) )

  call r8ncf_dif2_rowcol ( m, n, nz_num, rowcol )

  allocate ( a(1:nz_num) )

  call r8ncf_dif2 ( m, n, nz_num, rowcol, a )

  call r8ncf_print ( m, n, nz_num, rowcol, a, '  The R8NCF second difference matrix:' )

  deallocate ( a )
  deallocate ( rowcol )

  return
end
subroutine r8ncf_indicator_test ( )

!*****************************************************************************80
!
!! R8NCF_INDICATOR_TEST tests R8NCF_INDICATOR.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 July 2007
!
!  Author:
!
!   John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 7
  integer, parameter :: nz_num = 15

  real ( kind = rk ), dimension ( nz_num ) :: a
  integer, dimension ( 2, nz_num ) :: rowcol = reshape ( &
    (/ &
    1, 1, &
    2, 2, &
    3, 3, &
    4, 4, &
    5, 5, &
    2, 1, &
    5, 1, &
    1, 2, &
    5, 2, &
    1, 4, &
    2, 4, &
    3, 4, &
    4, 5, &
    4, 6, &
    1, 7 /), (/ 2, nz_num /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8NCF_INDICATOR_TEST'
  write ( *, '(a)' ) '  R8NCF_INDICATOR sets up a R8NCF indicator matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n
  write ( *, '(a,i8)' ) '  Matrix nonzeros =  ', nz_num

  call r8ncf_indicator ( m, n, nz_num, rowcol, a )

  call r8ncf_print ( m, n, nz_num, rowcol, a, '  The R8NCF indicator matrix:' )

  return
end
subroutine r8ncf_mtv_test ( )

!*****************************************************************************80
!
!! R8NCF_MTV_TEST tests R8NCF_MTV.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: a(:)
  real ( kind = rk ), allocatable :: b(:)
  integer m
  integer n
  integer nz_num
  integer, allocatable :: rowcol(:,:)
  real ( kind = rk ), allocatable :: x(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8NCF_MTV_TEST'
  write ( *, '(a)' ) '  R8NCF_MTV computes b=A''*x, where A is an R8NCF matrix.'

  m = 5
  n = 7

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order M =         ', m
  write ( *, '(a,i8)' ) '  Matrix order N =         ', n

  call r8ncf_dif2_nz_num ( m, n, nz_num )

  write ( *, '(a,i8)' ) '  Matrix nonzeros =  ', nz_num

  allocate ( rowcol(1:2,1:nz_num) )

  call r8ncf_dif2_rowcol ( m, n, nz_num, rowcol )

  allocate ( a(1:nz_num) )

  call r8ncf_dif2 ( m, n, nz_num, rowcol, a )

  allocate ( x(1:m) )
  call r8vec_indicator1 ( m, x )
  call r8vec_print ( m, x, '  x:' )

  allocate ( b(1:n) )
  call r8ncf_mtv ( m, n, nz_num, rowcol, a, x, b )

  call r8vec_print ( n, b, '  b=A''*x:' )

  deallocate ( a )
  deallocate ( b )
  deallocate ( rowcol )
  deallocate ( x )

  return
end
subroutine r8ncf_mv_test ( )

!*****************************************************************************80
!
!! R8NCF_MV_TEST tests R8NCF_MV.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: a(:)
  real ( kind = rk ), allocatable :: b(:)
  integer m
  integer n
  integer nz_num
  integer, allocatable :: rowcol(:,:)
  real ( kind = rk ), allocatable :: x(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8NCF_MV_TEST'
  write ( *, '(a)' ) '  R8NCF_MV computes b=A*x, where A is an R8NCF matrix.'

  m = 5
  n = 4

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order M =         ', m
  write ( *, '(a,i8)' ) '  Matrix order N =         ', n

  call r8ncf_dif2_nz_num ( m, n, nz_num )

  write ( *, '(a,i8)' ) '  Matrix nonzeros =  ', nz_num

  allocate ( rowcol(1:2,1:nz_num) )

  call r8ncf_dif2_rowcol ( m, n, nz_num, rowcol )

  allocate ( a(1:nz_num) )

  call r8ncf_dif2 ( m, n, nz_num, rowcol, a )

  allocate ( x(1:n) )
  call r8vec_indicator1 ( n, x )
  call r8vec_print ( n, x, '  x:' )

  allocate ( b(1:m) )
  call r8ncf_mv ( m, n, nz_num, rowcol, a, x, b )

  call r8vec_print ( m, b, '  b=A*x:' )

  deallocate ( a )
  deallocate ( b )
  deallocate ( rowcol )
  deallocate ( x )

  return
end
subroutine r8ncf_print_test ( )

!*****************************************************************************80
!
!! R8NCF_PRINT_TEST tests R8NCF_PRINT.
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
!   John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 7
  integer, parameter :: nz_num = 15

  real ( kind = rk ), dimension ( nz_num ) :: a
  integer, dimension ( 2, nz_num ) :: rowcol = reshape ( &
    (/ &
    1, 1, &
    2, 2, &
    3, 3, &
    4, 4, &
    5, 5, &
    2, 1, &
    5, 1, &
    1, 2, &
    5, 2, &
    1, 4, &
    2, 4, &
    3, 4, &
    4, 5, &
    4, 6, &
    1, 7 /), (/ 2, nz_num /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8NCF_PRINT_TEST'
  write ( *, '(a)' ) '  R8NCF_PRINT prints an R8NCF matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n
  write ( *, '(a,i8)' ) '  Matrix nonzeros =  ', nz_num

  call r8ncf_indicator ( m, n, nz_num, rowcol, a )

  call r8ncf_print ( m, n, nz_num, rowcol, a, '  The R8NCF matrix:' )

  return
end
subroutine r8ncf_print_some_test ( )

!*****************************************************************************80
!
!! R8NCF_PRINT_SOME_TEST tests R8NCF_PRINT_SOME.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 July 2016
!
!  Author:
!
!   John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 7
  integer, parameter :: nz_num = 15

  real ( kind = rk ), dimension ( nz_num ) :: a
  integer, dimension ( 2, nz_num ) :: rowcol = reshape ( &
    (/ &
    1, 1, &
    2, 2, &
    3, 3, &
    4, 4, &
    5, 5, &
    2, 1, &
    5, 1, &
    1, 2, &
    5, 2, &
    1, 4, &
    2, 4, &
    3, 4, &
    4, 5, &
    4, 6, &
    1, 7 /), (/ 2, nz_num /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8NCF_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R8NCF_PRINT_SOME prints some of an R8NCF matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n
  write ( *, '(a,i8)' ) '  Matrix nonzeros =  ', nz_num

  call r8ncf_indicator ( m, n, nz_num, rowcol, a )

  call r8ncf_print_some ( m, n, nz_num, rowcol, a, 1, 2, 5, 4, &
    '  Rows 1-5, Cols 2-4:' )

  return
end
subroutine r8ncf_random_test ( )

!*****************************************************************************80
!
!! R8NCF_RANDOM_TEST tests R8NCF_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 July 2016
!
!  Author:
!
!   John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 7
  integer, parameter :: nz_num = 15

  real ( kind = rk ), dimension ( nz_num ) :: a
  integer, dimension ( 2, nz_num ) :: rowcol = reshape ( &
    (/ &
    1, 1, &
    2, 2, &
    3, 3, &
    4, 4, &
    5, 5, &
    2, 1, &
    5, 1, &
    1, 2, &
    5, 2, &
    1, 4, &
    2, 4, &
    3, 4, &
    4, 5, &
    4, 6, &
    1, 7 /), (/ 2, nz_num /) )
  integer seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8NCF_RANDOM_TEST'
  write ( *, '(a)' ) '  R8NCF_RANDOM randomizes an R8NCF matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n
  write ( *, '(a,i8)' ) '  Matrix nonzeros =  ', nz_num

  seed = 123456789
  call r8ncf_random ( m, n, nz_num, rowcol, seed, a )

  call r8ncf_print ( m, n, nz_num, rowcol, a, '  The R8NCF matrix:' )

  return
end
subroutine r8ncf_to_r8ge_test ( )

!*****************************************************************************80
!
!! R8NCF_TO_R8GE_TEST tests R8NCF_TO_R8GE.
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
!   John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 7
  integer, parameter :: nz_num = 15

  real ( kind = rk ), dimension ( nz_num ) :: a
  real ( kind = rk ) a_r8ge(m,n)
  integer, dimension ( 2, nz_num ) :: rowcol = reshape ( &
    (/ &
    1, 1, &
    2, 2, &
    3, 3, &
    4, 4, &
    5, 5, &
    2, 1, &
    5, 1, &
    1, 2, &
    5, 2, &
    1, 4, &
    2, 4, &
    3, 4, &
    4, 5, &
    4, 6, &
    1, 7 /), (/ 2, nz_num /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8NCF_TO_R8GE_TEST'
  write ( *, '(a)' ) '  R8NCF_TO_R8GE converts an R8NCF matrix to R8GE format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n
  write ( *, '(a,i8)' ) '  Matrix nonzeros =  ', nz_num

  call r8ncf_indicator ( m, n, nz_num, rowcol, a )

  call r8ncf_print ( m, n, nz_num, rowcol, a, '  The R8NCF matrix:' )

  call r8ncf_to_r8ge ( m, n, nz_num, rowcol, a, a_r8ge )

  call r8ge_print ( m, n, a_r8ge, '  The R8GE matrix:' )

  return
end
subroutine r8ncf_zeros_test ( )

!*****************************************************************************80
!
!! R8NCF_ZEROS_TEST tests R8NCF_ZEROS.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 July 2016
!
!  Author:
!
!   John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 7
  integer, parameter :: nz_num = 15

  real ( kind = rk ), dimension ( nz_num ) :: a
  integer, dimension ( 2, nz_num ) :: rowcol = reshape ( &
    (/ &
    1, 1, &
    2, 2, &
    3, 3, &
    4, 4, &
    5, 5, &
    2, 1, &
    5, 1, &
    1, 2, &
    5, 2, &
    1, 4, &
    2, 4, &
    3, 4, &
    4, 5, &
    4, 6, &
    1, 7 /), (/ 2, nz_num /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8NCF_ZEROS_TEST'
  write ( *, '(a)' ) '  R8NCF_ZEROS zeros an R8NCF matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n
  write ( *, '(a,i8)' ) '  Matrix nonzeros =  ', nz_num

  call r8ncf_zeros ( m, n, nz_num, rowcol, a )

  call r8ncf_print ( m, n, nz_num, rowcol, a, '  The R8NCF zero matrix:' )

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


