subroutine ctof ( nc, uc, nf, uf )

!*****************************************************************************80
!                                                    
!! ctof() transfers data from a coarse to a finer grid.
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
!    integer NC, the number of coarse nodes.
!
!    real ( kind = rk8 ) UC(NC), the coarse correction data.
!
!    integer NF, the number of fine nodes.
!
!    real ( kind = rk8 ) UF(NF); the fine grid data.
!
!  Output:
!
!    real ( kind = rk8 ) UF(NF): the updated fine grid data.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer nc
  integer nf

  integer ic
  integer if
  real ( kind = rk8 ) uc(nc)
  real ( kind = rk8 ) uf(nf)

  do ic = 1, nc
    if = 2 * ic - 1
    uf(if) = uf(if) + uc(ic)
  end do

  do ic = 1, nc - 1
    if = 2 * ic
    uf(if) = uf(if) + 0.5D+00 * ( uc(ic) + uc(ic+1) )
  end do

  return
end
subroutine ftoc ( nf, uf, rf, nc, uc, rc )

!*****************************************************************************80
!                                                    
!! ftoc() transfers data from a fine grid to a coarser grid.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 December 2011
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
!    integer NF, the number of fine nodes.
!
!    real ( kind = rk8 ) UF(NF), the fine data.
!
!    real ( kind = rk8 ) RF(NF), the right hand side for the fine grid.
!
!    integer NC, the number of coarse nodes.
!
!  Output:
!
!    real ( kind = rk8 ) UC(NC), the coarse grid data, set to zero.
!
!    real ( kind = rk8 ) RC(NC), the right hand side for the coarse grid.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer nc
  integer nf

  integer ic
  integer if
  real ( kind = rk8 ) rc(nc)
  real ( kind = rk8 ) rf(nf)
  real ( kind = rk8 ) uc(nc)
  real ( kind = rk8 ) uf(nf)

  uc(1) = 0.0D+00
  rc(1) = 0.0D+00
  do ic = 2, nc - 1
    if = 2 * ic - 1
    uc(ic) = 0.0D+00
    rc(ic) = 4.0D+00 * ( rf(if) + uf(if-1) - 2.0D+00 * uf(if) + uf(if+1) )
  end do
  uc(nc) = 0.0D+00
  rc(nc) = 0.0D+00

  return
end
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
function i4_log_2 ( i )

!*****************************************************************************80
!
!! i4_log_2() returns the integer part of the logarithm base 2 of an I4.
!
!  Discussion:
!
!    For positive I4_LOG_2(I), it should be true that
!      2^I4_LOG_2(X) <= |I| < 2^(I4_LOG_2(I)+1).
!    The special case of I4_LOG_2(0) returns -HUGE().
!
!    An I4 is an integer value.
!
!  Example:
!
!     I  I4_LOG_2
!
!     0  -1
!     1,  0
!     2,  1
!     3,  1
!     4,  2
!     5,  2
!     6,  2
!     7,  2
!     8,  3
!     9,  3
!    10,  3
!   127,  6
!   128,  7
!   129,  7
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 January 2003
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer I, the number whose logarithm base 2 is desired.
!
!  Output:
!
!    integer I4_LOG_2, the integer part of the
!    logarithm base 2 of the absolute value of I.
!
  implicit none

  integer i
  integer i_abs
  integer i4_log_2
  integer, parameter :: i4_huge = 2147483647

  if ( i == 0 ) then

    i4_log_2 = - i4_huge

  else

    i4_log_2 = 0

    i_abs = abs ( i )

    do while ( 2 <= i_abs )
      i_abs = i_abs / 2
      i4_log_2 = i4_log_2 + 1
    end do

  end if

  return
end
subroutine poisson_1d_multigrid ( n, a, b, ua, ub, force, it_num, u )

!*****************************************************************************80
!                                                    
!! poisson_1d_multigrid() solves a 1D PDE using the multigrid method.
!
!  Discussion:
!
!    This routine solves a 1D boundary value problem of the form
!
!      - U''(X) = F(X) for A < X < B,
!
!    with boundary conditions U(A) = UA, U(B) = UB.
!
!    The multigrid method is used. 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 July 2014
!
!  Author:
!
!    Original FORTRAN77 version by William Hager.
!    This version by John Burkardt.
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
!    N must be a power of 2.
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
  real ( kind = rk8 ) d0
  real ( kind = rk8 ) d1
  real ( kind = rk8 ), external :: force
  real ( kind = rk8 ) h
  integer i
  integer i4_log_2
  integer it
  integer it_num
  integer j
  integer k
  integer l
  integer ll
  integer m
  integer nl
  real ( kind = rk8 ), allocatable :: r(:)
  real ( kind = rk8 ) tol
  real ( kind = rk8 ) u(n+1)
  real ( kind = rk8 ) ua
  real ( kind = rk8 ) ub
  real ( kind = rk8 ) utol
  real ( kind = rk8 ), allocatable :: uu(:)
  real ( kind = rk8 ) x(n+1)
!
!  Determine if we have enough storage.
!
  k = i4_log_2 ( n )

  if ( n /= 2**k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'poisson_1d_multigrid(): Fatal error!'
    write ( *, '(a)' ) '  N is not a power of 2.'
    stop 1
  end if

  nl = n + n + k - 2

  allocate ( uu(1:nl) )
  allocate ( r(1:nl) )
!
!  Initialization.
!
  it = 4
  it_num = 0
  tol = 0.0001D+00
  utol = 0.7D+00
  m = n
!
!  Set the nodes.
!
  call r8vec_linspace ( n + 1, a, b, x )
! 
!  Set the right hand side.
! 
  r(1) = ua
  h = ( b - a ) / real ( n, kind = rk8 )
  do i = 2, n
    r(i) = h**2 * force ( x(i) )
  end do
  r(n+1) = ub

  uu(1:nl) = 0.0D+00
!
!  L points to first entry of solution
!  LL points to penultimate entry.
!
  l = 1
  ll = n
! 
!  Gauss-Seidel iteration
!
  d1 = 0.0D+00
  j = 0

  do

    d0 = d1
    j = j + 1
    call gauss_seidel ( n + 1, r(l), uu(l), d1 )
    it_num = it_num + 1
!
!  Do at least 4 iterations at each level.
!
    if ( j < it ) then

      cycle
!
!  Enough iterations, satisfactory decrease, on finest grid, exit.
!
    else if ( d1 < tol .and. n == m ) then

      exit
!
!  Enough iterations, satisfactory convergence, go finer.
!
    else if ( d1 < tol ) then

      call ctof ( n + 1, uu(l), n + n + 1, uu(l-1-n-n) )

      n = n + n
      ll = l - 2
      l = l - 1 - n
      j = 0
!
!  Enough iterations, slow convergence, 2 < N, go coarser.
!
    else if ( utol * d0 <= d1 .and. 2 < n ) then

      call ftoc ( n + 1, uu(l), r(l), (n/2)+1, uu(l+n+1), r(l+n+1) )

      n = n / 2
      l = ll + 2
      ll = ll + n + 1
      j = 0

    end if

  end do

  u(1:n+1) = uu(1:n+1)

  deallocate ( r )
  deallocate ( uu )

  return
end
subroutine r8vec_linspace ( n, a, b, x )

!*****************************************************************************80
!
!! r8vec_linspace() creates a vector of linearly spaced values.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    4 points evenly spaced between 0 and 12 will yield 0, 4, 8, 12.
!
!    In other words, the interval is divided into N-1 even subintervals,
!    and the endpoints of intervals are used as the points.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of entries in the vector.
!
!    real ( kind = rk8 ) A, B, the first and last entries.
!
!  Output:
!
!    real ( kind = rk8 ) X(N), a vector of linearly spaced data.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  integer i
  real ( kind = rk8 ) x(n)

  if ( n == 1 ) then

    x(1) = ( a + b ) / 2.0D+00

  else

    do i = 1, n
      x(i) = ( real ( n - i,     kind = rk8 ) * a   &
             + real (     i - 1, kind = rk8 ) * b ) &
             / real ( n     - 1, kind = rk8 )
    end do

  end if

  return
end

