program main

!*****************************************************************************80
!
!! test_lls_test() tests test_lls().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 October 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_LLS_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TEST_LLS library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_LLS_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 summarizes the test data.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 April 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: a(:,:)
  real ( kind = rk ), allocatable :: b(:)
  real ( kind = rk ) b_norm
  integer m
  integer n
  integer prob
  integer prob_num
  real ( kind = rk ), allocatable :: r(:)
  real ( kind = rk ) r_norm
  real ( kind = rk ) r8vec_norm
  real ( kind = rk ), allocatable :: x(:)
  real ( kind = rk ) x_norm

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Get each least squares test and compute the maximum residual.'
  write ( *, '(a)' ) '  The L2 norm of the residual MUST be no greater than'
  write ( *, '(a)' ) '  the L2 norm of the right hand side, else 0 is a better solution.'

  call p00_prob_num ( prob_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  Number of problems = ', prob_num
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Index     M     N     ||B||         ||X||         ||R||'
  write ( *, '(a)' ) ' '

  do prob = 1, prob_num

    call p00_m ( prob, m )
    call p00_n ( prob, n )

    allocate ( a(1:m,1:n) )
    allocate ( b(1:m) )
    allocate ( x(1:n) )
    allocate ( r(1:m) )

    call p00_a ( prob, m, n, a )
    call p00_b ( prob, m, b )
    call p00_x ( prob, n, x )

    r(1:m) = matmul ( a(1:m,1:n), x(1:n) ) - b(1:m)

    b_norm = r8vec_norm ( m, b )
    x_norm = r8vec_norm ( n, x )
    r_norm = r8vec_norm ( m, r )

    write ( *, '(2x,i5,2x,i4,2x,i4,2x,g12.4,2x,g12.4,2x,g12.4)' ) prob, m, n, b_norm, x_norm, r_norm

    deallocate ( a )
    deallocate ( b )
    deallocate ( r )
    deallocate ( x )

  end do

  return
end
