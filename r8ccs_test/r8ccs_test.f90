program main

!*****************************************************************************80
!
!! r8ccs_test() tests r8ccs().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8ccs_test():'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test r8ccs().'

  call r8ccs_dif2_test ( )
  call r8ccs_get_test ( )
  call r8ccs_ijk_test ( )
  call r8ccs_inc_test ( )
  call r8ccs_indicator_test ( )
  call r8ccs_kij_test ( )
  call r8ccs_mtv_test ( )
  call r8ccs_mv_test ( )
  call r8ccs_print_test ( )
  call r8ccs_print_some_test ( )
  call r8ccs_random_test ( )
  call r8ccs_read_test ( )
  call r8ccs_set_test ( )
  call r8ccs_to_r8ge_test ( )
  call r8ccs_write_test ( )
  call r8ccs_zeros_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8ccs_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine r8ccs_dif2_test ( )

!*****************************************************************************80
!
!! r8ccs_DIF2_TEST tests r8ccs_DIF2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: a(:)
  integer, allocatable :: colptr(:)
  integer m
  integer n
  integer nz_num
  integer, allocatable :: rowind(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8ccs_DIF2_TEST'
  write ( *, '(a)' ) '  r8ccs_DIF2 sets the second difference as an r8ccs matrix;'

  m = 5
  n = 5
  if ( m == n ) then
    nz_num = 3 * m - 2
  else
    nz_num = 3 * min ( m, n ) - 1
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M     = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N  = ', n
  write ( *, '(a,i8)' ) '  Nonzeros NZ_NUM   = ', nz_num

  allocate ( a(1:nz_num) )
  allocate ( rowind(1:nz_num) )
  allocate ( colptr(1:n+1) )

  call r8ccs_dif2 ( m, n, nz_num, colptr, rowind, a )

  call r8ccs_print ( m, n, nz_num, colptr, rowind, a, &
    '  The r8ccs matrix:' )

  return
end
subroutine r8ccs_get_test ( )

!*****************************************************************************80
!
!! r8ccs_GET_TEST tests r8ccs_GET.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 5
  integer, parameter :: nz_num = 12

  real ( kind = rk ) a(nz_num)
  integer, dimension (n+1) :: colptr = (/ 1, 4, 6, 8, 10, 13 /)
  integer i
  integer i4_uniform_ab
  integer j
  integer k
  integer, dimension ( nz_num ) :: rowind = (/ &
    1, 2, 4, 1, 2, 3, 5, 4, 5, 1, 2, 5 /)
  integer :: seed = 123456789
  integer test
  real ( kind = rk ) value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8ccs_GET_TEST'
  write ( *, '(a)' ) '  r8ccs_GET gets an entry of a matrix in the r8ccs format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M     = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N  = ', n
  write ( *, '(a,i8)' ) '  Nonzeros NZ_NUM   = ', nz_num

  call i4vec_print ( n + 1, colptr, '  The COLPTR vector:' )

  call i4vec_print ( nz_num, rowind, '  The ROWIND vector:' )

  a(1:nz_num) = 0.0D+00
!
!  Initialize the matrix to random values.
!
  call r8ccs_random ( m, n, nz_num, colptr, rowind, a, seed )

  call r8ccs_print ( m, n, nz_num, colptr, rowind, a, &
    '  The  r8ccs matrix:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  r8ccs_GET retrieves 10 entries.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I         J         K      VALUE'
  write ( *, '(a)' ) ' '

  do test = 1, 10
    k = i4_uniform_ab ( 1, nz_num, seed )
    call r8ccs_kij ( m, n, nz_num, colptr, rowind, k, i, j )
    call r8ccs_get ( m, n, nz_num, colptr, rowind, a, i, j, value )
    write ( *, '(2x,i8,2x,i8,2x,i8,2x,g14.6)' ) i, j, k, value
  end do

  return
end
subroutine r8ccs_ijk_test ( )

!*****************************************************************************80
!
!! r8ccs_IJK_TEST tests r8ccs_IJK.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    29 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 5
  integer, parameter :: nz_num = 12

  real ( kind = rk ) a(nz_num)
  integer, dimension (n+1) :: colptr = (/ 1, 4, 6, 8, 10, 13 /)
  integer i
  integer i4_uniform_ab
  integer j
  integer k
  integer, dimension ( nz_num ) :: rowind = (/ &
    1, 2, 4, 1, 2, 3, 5, 4, 5, 1, 2, 5 /)
  integer :: seed = 123456789
  integer test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8ccs_IJK_TEST'
  write ( *, '(a)' ) '  r8ccs_IJK gets K from (I,J)'
  write ( *, '(a)' ) '  for a matrix in the r8ccs format,'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M     = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N  = ', n
  write ( *, '(a,i8)' ) '  Nonzeros NZ_NUM   = ', nz_num

  call i4vec_print ( n + 1, colptr, '  The COLPTR vector:' )

  call i4vec_print ( nz_num, rowind, '  The ROWIND vector:' )

  a(1:nz_num) = 0.0D+00
!
!  Initialize the matrix to random values.
!
  call r8ccs_random ( m, n, nz_num, colptr, rowind, a, seed )

  call r8ccs_print ( m, n, nz_num, colptr, rowind, a, &
    '  The initial r8ccs matrix:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  r8ccs_IJK locates some (I,J) entries.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I         J         K'
  write ( *, '(a)' ) ' '

  do test = 1, 20
    i = i4_uniform_ab ( 1, m, seed )
    j = i4_uniform_ab ( 1, n, seed )
    call r8ccs_ijk ( m, n, nz_num, colptr, rowind, i, j, k )
    write ( *, '(2x,i8,2x,i8,2x,i8)' ) i, j, k
  end do

  return
end
subroutine r8ccs_inc_test ( )

!*****************************************************************************80
!
!! r8ccs_INC_TEST tests r8ccs_INC.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 5
  integer, parameter :: nz_num = 12

  real ( kind = rk ) a(nz_num)
  integer, dimension (n+1) :: colptr = (/ 1, 4, 6, 8, 10, 13 /)
  integer i
  integer i4_uniform_ab
  integer j
  integer k
  integer, dimension ( nz_num ) :: rowind = (/ &
    1, 2, 4, 1, 2, 3, 5, 4, 5, 1, 2, 5 /)
  integer :: seed = 123456789
  integer test
  integer, parameter :: test_num = 20
  real ( kind = rk ) value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8ccs_INC_TEST'
  write ( *, '(a)' ) '  r8ccs_INC increments entries in an r8ccs matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M     = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N  = ', n
  write ( *, '(a,i8)' ) '  Nonzeros NZ_NUM   = ', nz_num

  call i4vec_print ( n + 1, colptr, '  The COLPTR vector:' )

  call i4vec_print ( nz_num, rowind, '  The ROWIND vector:' )

  a(1:nz_num) = 0.0D+00
!
!  Initialize the matrix to random values.
!
  call r8ccs_random ( m, n, nz_num, colptr, rowind, a, seed )

  call r8ccs_print ( m, n, nz_num, colptr, rowind, a, &
    '  The initial r8ccs matrix:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  r8ccs_INC increments 10 entries at random.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I         J         K       NEW_VALUE'
  write ( *, '(a)' ) ' '

  do test = 1, 10
    k = i4_uniform_ab ( 1, nz_num, seed )
    call r8ccs_kij ( m, n, nz_num, colptr, rowind, k, i, j )
    value = 20.0D+00 + real ( test, kind = rk )
    call r8ccs_inc ( m, n, nz_num, colptr, rowind, a, i, j, value )
    call r8ccs_get ( m, n, nz_num, colptr, rowind, a, i, j, value )
    write ( *, '(2x,i8,2x,i8,2x,i8,2x,g14.6)' ) i, j, k, value
  end do

  call r8ccs_print ( m, n, nz_num, colptr, rowind, a, &
    '  The final r8ccs matrix:' )

  return
end
subroutine r8ccs_indicator_test ( )

!*****************************************************************************80
!
!! r8ccs_INDICATOR_TEST tests r8ccs_INDICATOR.
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
  integer, parameter :: n = 5
  integer, parameter :: nz_num = 12

  real ( kind = rk ) a(nz_num)
  integer, dimension (n+1) :: colptr = (/ 1, 4, 6, 8, 10, 13 /)
  integer, dimension ( nz_num ) :: rowind = (/ &
    1, 2, 4, 1, 2, 3, 5, 4, 5, 1, 2, 5 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8ccs_INDICATOR_TEST'
  write ( *, '(a)' ) '  r8ccs_INDICATOR sets an indicator r8ccs matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M     = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N  = ', n
  write ( *, '(a,i8)' ) '  Nonzeros NZ_NUM   = ', nz_num

  call r8ccs_indicator ( m, n, nz_num, colptr, rowind, a )

  call r8ccs_print ( m, n, nz_num, colptr, rowind, a, &
    '  The r8ccs indicator matrix:' )

  return
end
subroutine r8ccs_kij_test ( )

!*****************************************************************************80
!
!! r8ccs_KIJ_TEST tests r8ccs_KIJ.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 5
  integer, parameter :: nz_num = 12

  real ( kind = rk ) a(nz_num)
  integer, dimension (n+1) :: colptr = (/ 1, 4, 6, 8, 10, 13 /)
  integer i
  integer i4_uniform_ab
  integer j
  integer k
  integer, dimension ( nz_num ) :: rowind = (/ &
    1, 2, 4, 1, 2, 3, 5, 4, 5, 1, 2, 5 /)
  integer :: seed = 123456789
  integer test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8ccs_KIJ_TEST'
  write ( *, '(a)' ) '  r8ccs_KIJ gets (I,J) from K'
  write ( *, '(a)' ) '  for a matrix in the r8ccs format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M     = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N  = ', n
  write ( *, '(a,i8)' ) '  Nonzeros NZ_NUM   = ', nz_num

  call i4vec_print ( n + 1, colptr, '  The COLPTR vector:' )

  call i4vec_print ( nz_num, rowind, '  The ROWIND vector:' )

  a(1:nz_num) = 0.0D+00
!
!  Initialize the matrix to random values.
!
  call r8ccs_random ( m, n, nz_num, colptr, rowind, a, seed )

  call r8ccs_print ( m, n, nz_num, colptr, rowind, a, &
    '  The initial r8ccs matrix:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  r8ccs_KIJ locates some K entries.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         K         I         J'
  write ( *, '(a)' ) ' '

  do test = 1, 20  
    k = i4_uniform_ab ( 1, nz_num, seed )
    call r8ccs_kij ( m, n, nz_num, colptr, rowind, k, i, j )
    write ( *, '(2x,i8,2x,i8,2x,i8)' ) k, i, j
  end do

  return
end
subroutine r8ccs_mtv_test ( )

!*****************************************************************************80
!
!! r8ccs_MTV_TEST tests r8ccs_MTV.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    24 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 5
  integer, parameter :: nz_num = 12

  real ( kind = rk ) a(nz_num)
  real ( kind = rk ) b(n)
  integer, dimension (n+1) :: colptr = (/ 1, 4, 6, 8, 10, 13 /)
  integer, dimension ( nz_num ) :: rowind = (/ &
    1, 2, 4, 1, 2, 3, 5, 4, 5, 1, 2, 5 /)
  integer :: seed = 123456789
  real ( kind = rk ) x(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8ccs_MTV_TEST'
  write ( *, '(a)' ) '  r8ccs_MTV compute b=A''*x, where A is an r8ccs matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M     = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N  = ', n
  write ( *, '(a,i8)' ) '  Nonzeros NZ_NUM   = ', nz_num
!
!  Set the matrix.
!
  call r8ccs_random ( m, n, nz_num, colptr, rowind, a, seed )
!
!  Compute BN = XM * A(MxN)
!
  x(1) = 2.0D+00
  x(2:m-1) = 0.0D+00
  x(m) = -3.0D+00

  call r8vec_print ( m, x, '  x:' )

  call r8ccs_mtv ( m, n, nz_num, colptr, rowind, a, x, b )

  call r8vec_print ( n, b, '  b=A''*x:' )

  return
end
subroutine r8ccs_mv_test ( )

!*****************************************************************************80
!
!! r8ccs_MV_TEST tests r8ccs_MV.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    24 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 5
  integer, parameter :: nz_num = 12

  real ( kind = rk ) a(nz_num)
  real ( kind = rk ) b(m)
  integer, dimension (n+1) :: colptr = (/ 1, 4, 6, 8, 10, 13 /)
  integer, dimension ( nz_num ) :: rowind = (/ &
    1, 2, 4, 1, 2, 3, 5, 4, 5, 1, 2, 5 /)
  integer :: seed = 123456789
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8ccs_MV_TEST'
  write ( *, '(a)' ) '  r8ccs_MV computes b=A*x, where A is an r8ccs matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M     = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N  = ', n
  write ( *, '(a,i8)' ) '  Nonzeros NZ_NUM   = ', nz_num
!
!  Set the matrix.
!
  call r8ccs_random ( m, n, nz_num, colptr, rowind, a, seed )
!
!  Compute B = A(MxN) * X
!
  x(1) = 1.0D+00
  x(2:n-1) = 0.0D+00
  x(n) = -1.0D+00

  call r8vec_print ( n, x, '  x:' )

  call r8ccs_mv ( m, n, nz_num, colptr, rowind, a, x, b )

  call r8vec_print ( m, b, '  b=A*x:' )

  return
end
subroutine r8ccs_print_test ( )

!*****************************************************************************80
!
!! r8ccs_PRINT_TEST tests r8ccs_PRINT.
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
  integer, parameter :: n = 5
  integer, parameter :: nz_num = 12

  real ( kind = rk ) a(nz_num)
  integer, dimension (n+1) :: colptr = (/ 1, 4, 6, 8, 10, 13 /)
  integer, dimension ( nz_num ) :: rowind = (/ &
    1, 2, 4, 1, 2, 3, 5, 4, 5, 1, 2, 5 /)
  integer :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8ccs_PRINT_TEST'
  write ( *, '(a)' ) '  r8ccs_PRINT prints an r8ccs matrix.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M     = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N  = ', n
  write ( *, '(a,i8)' ) '  Nonzeros NZ_NUM   = ', nz_num
!
!  Set the matrix.
!
  call r8ccs_random ( m, n, nz_num, colptr, rowind, a, seed )
!
!  Print the matrix.
!
  call r8ccs_print ( m, n, nz_num, colptr, rowind, a, '  The r8ccs matrix:' )

  return
end
subroutine r8ccs_print_some_test ( )

!*****************************************************************************80
!
!! r8ccs_PRINT_SOME_TEST tests r8ccs_PRINT_SOME.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 10
  integer, parameter :: n = 10
  integer, parameter :: nz_num = 28

  real ( kind = rk ) a(nz_num)
  integer, dimension (n+1) :: colptr = (/ &
    1, 3, 6, 9, 12, 15, 18, 21, 24, 27, 29 /)
  integer, dimension ( nz_num ) :: rowind = (/ &
    1,  2,  &
    1,  2,  3, &
    2,  3,  4, &
    3,  4,  5, &
    4,  5,  6, &
    5,  6,  7, &
    6,  7,  8, &
    7,  8,  9, &
    8,  9, 10, &
    9, 10 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8ccs_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  r8ccs_PRINT_SOME prints some of an r8ccs matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M     = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N  = ', n
  write ( *, '(a,i8)' ) '  Nonzeros NZ_NUM   = ', nz_num

  call r8ccs_indicator ( m, n, nz_num, colptr, rowind, a )

  call r8ccs_print_some ( m, n, nz_num, colptr, rowind, a, &
    2, 5, 6, 8, '  Rows 2-6, Cols 5-8:' )

  return
end
subroutine r8ccs_random_test ( )

!*****************************************************************************80
!
!! r8ccs_RANDOM_TEST tests r8ccs_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    29 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 5
  integer, parameter :: nz_num = 12

  real ( kind = rk ) a(nz_num)
  integer, dimension (n+1) :: colptr = (/ 1, 4, 6, 8, 10, 13 /)
  integer, dimension ( nz_num ) :: rowind = (/ &
    1, 2, 4, 1, 2, 3, 5, 4, 5, 1, 2, 5 /)
  integer seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8ccs_RANDOM_TEST'
  write ( *, '(a)' ) '  r8ccs_RANDOM randomizes an r8ccs matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M     = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N  = ', n
  write ( *, '(a,i8)' ) '  Nonzeros NZ_NUM   = ', nz_num

  seed = 123456789
  call r8ccs_random ( m, n, nz_num, colptr, rowind, a, seed )

  call r8ccs_print ( m, n, nz_num, colptr, rowind, a, &
    '  The r8ccs matrix:' )

  return
end
subroutine r8ccs_read_test ( )

!*****************************************************************************80
!
!! r8ccs_READ_TEST tests r8ccs_READ.
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

  real ( kind = rk ), allocatable, dimension ( : ) :: a
  character ( len = 127 ) :: a_file = 'r8ccs_a.txt'
  integer base
  integer, allocatable, dimension ( : ) :: col
  character ( len = 127 ) :: col_file = 'r8ccs_col.txt'
  integer m
  integer n
  integer nz_num
  integer, allocatable, dimension ( : ) :: row
  character ( len = 127 ) :: row_file = 'r8ccs_row.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8ccs_READ_TEST'
  write ( *, '(a)' ) '  r8ccs_READ reads an r8ccs matrix from 3 files.'

  call r8ccs_read_size ( col_file, row_file, m, n, nz_num, base )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M     = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N  = ', n
  write ( *, '(a,i8)' ) '  Nonzeros NZ_NUM   = ', nz_num
  write ( *, '(a,i8)' ) '  Index base (0/1)  = ', base

  allocate ( a(1:nz_num) )
  allocate ( col(1:n+1) )
  allocate ( row(1:nz_num) )

  call r8ccs_read ( col_file, row_file, a_file, m, n, nz_num, col, row, a )

  call i4vec_print ( n + 1, col, '  The COL vector:' )

  call i4vec_print ( nz_num, row, '  The ROW vector:' )

  call r8ccs_print ( m, n, nz_num, col, row, a, '  The r8ccs matrix:' )

  deallocate ( a )
  deallocate ( col )
  deallocate ( row )

  return
end
subroutine r8ccs_set_test ( )

!*****************************************************************************80
!
!! r8ccs_SET_TEST tests r8ccs_SET.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 5
  integer, parameter :: nz_num = 12

  real ( kind = rk ) a(nz_num)
  integer, dimension (n+1) :: colptr = (/ 1, 4, 6, 8, 10, 13 /)
  integer i
  integer i4_uniform_ab
  integer j
  integer k
  integer, dimension ( nz_num ) :: rowind = (/ &
    1, 2, 4, 1, 2, 3, 5, 4, 5, 1, 2, 5 /)
  integer :: seed = 123456789
  integer test
  integer, parameter :: test_num = 20
  real ( kind = rk ) value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8ccs_SET_TEST'
  write ( *, '(a)' ) '  r8ccs_SET sets entries in an r8ccs matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M     = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N  = ', n
  write ( *, '(a,i8)' ) '  Nonzeros NZ_NUM   = ', nz_num

  call i4vec_print ( n + 1, colptr, '  The COLPTR vector:' )

  call i4vec_print ( nz_num, rowind, '  The ROWIND vector:' )

  a(1:nz_num) = 0.0D+00
!
!  Initialize the matrix to random values.
!
  call r8ccs_random ( m, n, nz_num, colptr, rowind, a, seed )

  call r8ccs_print ( m, n, nz_num, colptr, rowind, a, &
    '  The initial r8ccs matrix:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  r8ccs_SET sets 10 entries at random.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I         J         K      NEW_VALUE'
  write ( *, '(a)' ) ' '

  do test = 1, 10
    k = i4_uniform_ab ( 1, nz_num, seed )
    call r8ccs_kij ( m, n, nz_num, colptr, rowind, k, i, j )
    value = 100.0D+00 + real ( test, kind = rk )
    call r8ccs_set ( m, n, nz_num, colptr, rowind, a, i, j, value )
    write ( *, '(2x,i8,2x,i8,2x,i8,2x,g14.6)' ) i, j, k, value
  end do

  return
end
subroutine r8ccs_to_r8ge_test ( )

!*****************************************************************************80
!
!! r8ccs_TO_R8GE_TEST tests r8ccs_TO_R8GE.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 5
  integer, parameter :: nz_num = 12

  real ( kind = rk ) a_r8ccs(nz_num)
  real ( kind = rk ) a_r8ge(m,n)
  integer, dimension (n+1) :: colptr = (/ 1, 4, 6, 8, 10, 13 /)
  integer, dimension ( nz_num ) :: rowind = (/ &
    1, 2, 4, 1, 2, 3, 5, 4, 5, 1, 2, 5 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8ccs_TO_R8GE_TEST'
  write ( *, '(a)' ) '  r8ccs_TO_R8GE converts a matrix from r8ccs to R8GE format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M     = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N  = ', n
  write ( *, '(a,i8)' ) '  Nonzeros NZ_NUM   = ', nz_num

  call r8ccs_indicator ( m, n, nz_num, colptr, rowind, a_r8ccs )

  call r8ccs_print ( m, n, nz_num, colptr, rowind, a_r8ccs, &
    '  The r8ccs matrix:' )

  call r8ccs_to_r8ge ( m, n, nz_num, colptr, rowind, a_r8ccs, a_r8ge )

  call r8ge_print ( m, n, a_r8ge, '  The R8GE matrix:' )

  return
end
subroutine r8ccs_write_test ( )

!*****************************************************************************80
!
!! r8ccs_WRITE_TEST tests r8ccs_WRITE.
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

  integer, parameter :: m = 5
  integer, parameter :: n = 5
  integer, parameter :: nz_num = 12

  real ( kind = rk ) a(nz_num)
  character ( len = 127 ) :: a_file = 'r8ccs_a.txt'
  integer, dimension (n+1) :: col = (/ 1, 4, 6, 8, 10, 13 /)
  character ( len = 127 ) :: col_file = 'r8ccs_col.txt'
  integer, dimension ( nz_num ) :: row = (/ &
    1, 2, 4, 1, 2, 3, 5, 4, 5, 1, 2, 5 /)
  character ( len = 127 ) :: row_file = 'r8ccs_row.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8ccs_WRITE_TEST'
  write ( *, '(a)' ) '  r8ccs_WRITE writes an r8ccs matrix to 3 files.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M     = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N  = ', n
  write ( *, '(a,i8)' ) '  Nonzeros NZ_NUM   = ', nz_num

  call i4vec_print ( n + 1, col, '  The COL vector:' )

  call i4vec_print ( nz_num, row, '  The ROW vector:' )

  call r8ccs_indicator ( m, n, nz_num, col, row, a )

  call r8ccs_print ( m, n, nz_num, col, row, a, '  The r8ccs matrix:' )

  call r8ccs_write ( col_file, row_file, a_file, m, n, nz_num, col, row, a )

  return
end
subroutine r8ccs_zeros_test ( )

!*****************************************************************************80
!
!! r8ccs_ZEROS_TEST tests r8ccs_ZEROS.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    29 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 5
  integer, parameter :: nz_num = 12

  real ( kind = rk ) a(nz_num)
  integer, dimension (n+1) :: colptr = (/ 1, 4, 6, 8, 10, 13 /)
  integer, dimension ( nz_num ) :: rowind = (/ &
    1, 2, 4, 1, 2, 3, 5, 4, 5, 1, 2, 5 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8ccs_ZEROS_TEST'
  write ( *, '(a)' ) '  r8ccs_ZEROS zeros an r8ccs matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M     = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N  = ', n
  write ( *, '(a,i8)' ) '  Nonzeros NZ_NUM   = ', nz_num

  call r8ccs_zeros ( m, n, nz_num, colptr, rowind, a )

  call r8ccs_print ( m, n, nz_num, colptr, rowind, a, &
    '  The r8ccs matrix:' )

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
