program main

!*****************************************************************************80
!
!! ccs_to_st_test() tests ccs_to_st().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 July 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ccs_to_st_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test ccs_to_st().'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ccs_to_st_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! test01() tests ccs_to_st() using a 1-based matrix.
!
!  Discussion:
!
!    This test uses a trivial matrix whose full representation is:
!
!          2  3  0  0  0
!          3  0  4  0  6
!      A = 0 -1 -3  2  0
!          0  0  1  0  0
!          0  4  2  0  1
!
!    The 1-based CCS representation is
!
!      #  ICC  CCC  ACC
!     --  ---  ---  ---
!      1    1    1    2
!      2    2         3
!
!      3    1    3    3
!      4    3        -1
!      5    5         4
!
!      6    2    6    4
!      7    3        -3
!      8    4         1
!      9    5         2
!
!     10    3   10    2
!
!     11    2   11    6
!     12    5         1
!
!     13    *   13
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 July 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 5
  integer, parameter :: ncc = 12

  real ( kind = rk ), dimension ( ncc ) :: acc = (/ &
    2.0,  3.0, &
    3.0, -1.0,  4.0, &
    4.0, -3.0,  1.0, 2.0, &
    2.0, &
    6.0, 1.0 /)
  real ( kind = rk ) ast(ncc)
  integer, dimension ( n + 1 ) :: ccc = (/ &
    1, 3, 6, 10, 11, 13 /)
  integer, dimension ( ncc ) :: icc = (/ &
    1, 2, &
    1, 3, 5, &
    2, 3, 4, 5, &
    3, &
    2, 5 /)
  integer ist(ncc)
  integer jst(ncc)
  integer nst

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Convert a 1-based CCS matrix to ST format.'
!
!  Print the CCS matrix.
!
  call ccs_print ( m, n, ncc, icc, ccc, acc, '  The CCS matrix:' )
!
!  Convert it.
!
  call ccs_to_st ( m, n, ncc, icc, ccc, acc, nst, ist, jst, ast )
!
!  Print the ST matrix.
!
  call st_print ( m, n, nst, ist, jst, ast, '  The ST matrix:' )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! test02() tests ccs_to_st() using a 0-based matrix.
!
!  Discussion:
!
!    This test uses a trivial matrix whose full representation is:
!
!          2  3  0  0  0
!          3  0  4  0  6
!      A = 0 -1 -3  2  0
!          0  0  1  0  0
!          0  4  2  0  1
!
!    The 0-based CCS representation is
!
!      #  ICC  CCC  ACC
!     --  ---  ---  ---
!      0    0    0    2
!      1    1         3
!
!      2    0    2    3
!      3    2        -1
!      4    4         4
!
!      5    1    5    4
!      6    2        -3
!      7    3         1
!      8    4         2
!
!      9    2    9    2
!
!     10    1   10    6
!     11    4         1
!
!     12    *   12
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 July 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 5
  integer, parameter :: ncc = 12

  real ( kind = rk ), dimension ( ncc ) :: acc = (/ &
    2.0,  3.0, &
    3.0, -1.0,  4.0, &
    4.0, -3.0,  1.0, 2.0, &
    2.0, &
    6.0, 1.0 /)
  real ( kind = rk ) ast(ncc)
  integer, dimension ( n + 1 ) :: ccc = (/ &
    0, 2, 5, 9, 10, 12 /)
  integer, dimension ( ncc ) :: icc = (/ &
    0, 1, &
    0, 2, 4, &
    1, 2, 3, 4, &
    2, &
    1, 4 /)
  integer ist(ncc)
  integer jst(ncc)
  integer nst

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Convert a 0-based CCS matrix to ST format.'
!
!  Print the CCS matrix.
!
  call ccs_print ( m, n, ncc, icc, ccc, acc, '  The CCS matrix:' )
!
!  Convert it.
!
  call ccs_to_st ( m, n, ncc, icc, ccc, acc, nst, ist, jst, ast )
!
!  Print the ST matrix.
!
  call st_print ( m, n, nst, ist, jst, ast, '  The ST matrix:' )

  return
end
