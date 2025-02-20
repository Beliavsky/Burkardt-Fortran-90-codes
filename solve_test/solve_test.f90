program main

!*****************************************************************************80
!
!! solve_test() tests solve().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    29 May 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'solve_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test solve().'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'solve_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! test01() demonstrates how a 3x3 linear system can be set up and solved.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 3

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) b(n)
  integer i
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'test01():'
  write ( *, '(a)' ) '  r8mat_fs() solves a linear system with Gauss elimination.'
!
!  Set the array.
!
  a(1,1) = 1.0
  a(1,2) = 2.0
  a(1,3) = 3.0

  a(2,1) = 4.0
  a(2,2) = 5.0
  a(2,3) = 6.0

  a(3,1) = 7.0
  a(3,2) = 8.0
  a(3,3) = 0.0
!
!  Set the right hand side.
!
  b(1) = 14.0
  b(2) = 32.0
  b(3) = 23.0
!
!  Request the solution of A*x=b.
!
  call r8mat_fs ( n, a, b, x )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Computed solution of linear system:'
  write ( *, '(a)' ) ''
  do i = 1, n
    write ( *, '(g14.6)' ) x(i)
  end do

  return
end


