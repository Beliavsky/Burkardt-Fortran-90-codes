program main

!*****************************************************************************80
!
!! ccs_io_test() tests ccs_io().
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

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'ccs_io_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test ccs_io().'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'ccs_io_test()'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! test01() tests ccs_write() using a tiny matrix.
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
  integer, dimension ( n + 1 ) :: ccc = (/ &
    1, 3, 6, 10, 11, 13 /)
  integer, dimension ( ncc ) :: icc = (/ &
    1, 2, &
    1, 3, 5, &
    2, 3, 4, 5, &
    3, &
    2, 5 /)
  character ( len = 255 ) :: prefix = 'simple'

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Write a sparse matrix in CCS format to 3 files.'
!
!  Full storage statistics
!
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Full rows    M = ', m
  write ( *, '(a,i4)' ) '  Full columns N = ', n
  write ( *, '(a,i4)' ) '  Full storage   = ', m * n
!
!  Print the CCS matrix.
!
  call ccs_print ( m, n, ncc, icc, ccc, acc, '  The matrix in 1-based CCS format:' )
!
!  Write the matrix to 3 files.
!
  call ccs_write ( prefix, ncc, n, icc, ccc, acc )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! test02() tests ccs_header_read() and ccs_data_read().
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

  real ( kind = rk ), allocatable :: acc(:)
  integer, allocatable :: ccc(:)
  integer, allocatable :: icc(:)
  integer m
  integer n
  integer ncc
  character ( len = 255 ) :: prefix = 'simple'

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Read a sparse matrix in CCS format from 3 files.'
!
!  Read the header.
!
  call ccs_header_read ( prefix, ncc, n )
!
!  Allocate space.
!
  allocate ( acc(1:ncc) )
  allocate ( ccc(1:n+1) )
  allocate ( icc(1:ncc) )
!
!  Read the matrix data.
!
  call ccs_data_read ( prefix, ncc, n, icc, ccc, acc )
!
!  Print the CCS matrix.
!
  m = n
  call ccs_print ( m, n, ncc, icc, ccc, acc, '  The matrix in 1-based CCS format:' )
!
!  Free memory.
!
  deallocate ( acc )
  deallocate ( ccc )
  deallocate ( icc )

  return
end
