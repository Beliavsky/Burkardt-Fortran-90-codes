subroutine gauss_seidel ( n, r, u, dif_l1 )

!*****************************************************************************80
!                                                    
!! gauss_seidel() carries out one step of a Gauss-Seidel iteration.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Hager,
!    Applied Numerical Linear Algebra,
!    Prentice-Hall, 1988,
!    ISBN13: 978-0130412942,
!    LC: QA184.H33.
!
!  Input:
!
!    integer N, the number of unknowns.
!
!    real ( kind = rk8 ) R(N), the right hand side.
!
!    real ( kind = rk8 ) U(N), the estimated solution.
!
!  Output:
!
!    real ( kind = rk8 ) U(N), the estimated solution.
!
!    real ( kind = rk8 ) DIF_L1, the L1 norm of the difference between the
!    input and output solution estimates.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  real ( kind = rk8 ) dif_l1
  integer i
  real ( kind = rk8 ) r(n)
  real ( kind = rk8 ) u(n)
  real ( kind = rk8 ) u_old

  dif_l1 = 0.0D+00

  do i = 2, n - 1
    u_old = u(i)
    u(i) = 0.5D+00 * ( u(i-1) + u(i+1) + r(i) )
    dif_l1 = dif_l1 + abs ( u(i) - u_old )
  end do

  return
end
subroutine poisson_1d ( n, a, b, ua, ub, force, it_num, u )

!*****************************************************************************80
!                                                    
!! poisson_1d() solves a 1D PDE, using the Gauss-Seidel method.
!
!  Discussion:
!
!    This routine solves a 1D boundary value problem of the form
!
!      - U''(X) = F(X) for A < X < B,
!
!    with boundary conditions U(A) = UA, U(B) = UB.
!
!    The Gauss-Seidel method is used. 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 October 2024
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Hager,
!    Applied Numerical Linear Algebra,
!    Prentice-Hall, 1988,
!    ISBN13: 978-0130412942,
!    LC: QA184.H33.
!
!  Input:
!
!    integer N, the number of intervals.
!
!    real ( kind = rk8 ) A, B, the endpoints.
!
!    real ( kind = rk8 ) UA, UB, the boundary values at the endpoints.
!
!    external real ( kind = rk8 ) FORCE, the name of the function 
!    which evaluates the right hand side.
!
!  Output:
!
!    integer IT_NUM, the number of iterations.
!
!    real ( kind = rk8 ) U(N+1), the computed solution.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  real ( kind = rk8 ) d1
  real ( kind = rk8 ), external :: force
  real ( kind = rk8 ) h
  integer i
  integer it_num
  real ( kind = rk8 ) r(n+1)
  real ( kind = rk8 ) tol
  real ( kind = rk8 ) u(n+1)
  real ( kind = rk8 ) ua
  real ( kind = rk8 ) ub
  real ( kind = rk8 ) x(n+1)
!
!  Initialization.
!
  tol = 0.0001D+00
!
!  Set the nodes.
!
  do i = 1, n + 1
    x(i) = ( ( n + 1 - i ) * a + ( i - 1 ) * b ) / n
  end do
!
!  Set the right hand side.
!
  r(1) = ua
  h = ( b - a ) / real ( n, kind = rk8 )
  do i = 2, n
    r(i) = h**2 * force ( x(i) )
  end do
  r(n+1) = ub
!
!  Initialize the solution.
!
  u(1) = ua
  u(2:n) = 0.0D+00
  u(n+1) = ub
!
!  Gauss-Seidel iteration.
!
  it_num = 0

  do

    it_num = it_num + 1

    call gauss_seidel ( n + 1, r, u, d1 )

    if ( d1 <= tol ) then
      exit
    end if

  end do

  return
end

