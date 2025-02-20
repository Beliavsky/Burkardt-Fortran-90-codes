program main

!*****************************************************************************80
!
!! jacobi_test() tests jacobi().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'jacobi_test():'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) '  Test jacobi().'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'jacobi_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!  Purpose:
!
!    jacobi_test01() tests jacobi1().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0d+00 )

  real( kind = rk ), allocatable :: a(:,:)
  real( kind = rk ), allocatable :: b(:)
  real( kind = rk ) d_max
  real( kind = rk ) e_max
  integer i
  integer it
  integer it_max
  integer n
  real( kind = rk ), allocatable :: r(:)
  real( kind = rk ) r_max
  real( kind = rk ) t
  real( kind = rk ) tol
  real( kind = rk ), allocatable :: x(:)
  real( kind = rk ), allocatable :: x_exact(:)
  real( kind = rk ), allocatable :: x_new(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'jacobi_test01():'
  write ( *, '(a)' ) '  Try the Jacobi iteration on the second difference matrix'

  it_max = 2000
  n = 33
  tol = 1.0D-05

  allocate ( a(1:n,1:n) )
  allocate ( b(1:n) )
  allocate ( r(1:n) )
  allocate ( x(1:n) )
  allocate ( x_exact(1:n) )
  allocate ( x_new(1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  Matrix order N = ', n
  write ( *, '(a,i6)' ) '  Maximum number of iterations = ', it_max
!
!  Set the matrix A.
!
  call dif2 ( n, n, a )
!
!  Determine the right hand side vector B.
!
  do i = 1, n
    t = real ( i - 1, kind = rk ) / real ( n - 1, kind = rk )
    x_exact(i) = exp ( t ) * ( t - 1 ) * t
  end do

  b = matmul ( a(1:n,1:n), x_exact(1:n) )
!
!  Set the initial estimate for the solution.
!
  it = 0

  x(1:n) = 0.0D+00

  r = matmul ( a(1:n,1:n), x(1:n) )
  r(1:n) = b(1:n) - r(1:n)
  r_max = maxval ( r )

  e_max = maxval ( abs ( x - x_exact ) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       I    Resid           X-Error         X-Change'
  write ( *, '(a)' ) ' '

  write ( *, '(2x,i6,2x,g14.6,2x,g14.6)' ) it, r_max, e_max
!
!  Carry out the iteration.
!
  do it = 1, it_max

    call jacobi1 ( n, a, b, x, x_new )

    r = matmul ( a(1:n,1:n), x(1:n) )
    r(1:n) = b(1:n) - r(1:n)
    r_max = maxval ( abs ( r ) )
!
!  Compute the average point motion.
!
    d_max = maxval ( abs ( x - x_new ) )
!
!  Compute the average point motion.
!
    e_max = maxval ( abs ( x - x_exact ) )
!
!  Update the solution
!
    x(1:n) = x_new(1:n)

    write ( *, '(2x,i6,2x,g14.6,2x,g14.6,2x,g14.6)' ) it, r_max, e_max, d_max

    if ( r_max <= tol ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Convergence criterion satisifed on step ', it
      exit
    end if

  end do

  call r8vec_print ( n, x, '  Estimated solution:' )

  deallocate ( a )
  deallocate ( b )
  deallocate ( r )
  deallocate ( x )
  deallocate ( x_exact )
  deallocate ( x_new )

  return
end
