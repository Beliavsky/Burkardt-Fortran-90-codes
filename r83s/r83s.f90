function i4_log_10 ( i )

!*****************************************************************************80
!
!! i4_log_10() returns the integer part of the logarithm base 10 of an I4.
!
!  Example:
!
!        I  I4_LOG_10
!    -----  --------
!        0    0
!        1    0
!        2    0
!        9    0
!       10    1
!       11    1
!       99    1
!      100    2
!      101    2
!      999    2
!     1000    3
!     1001    3
!     9999    3
!    10000    4
!
!  Discussion:
!
!    I4_LOG_10 ( I ) + 1 is the number of decimal digits in I.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    08 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, the number whose logarithm base 10 
!    is desired.
!
!    Output, integer I4_LOG_10, the integer part of the 
!    logarithm base 10 of the absolute value of X.
!
  implicit none

  integer i
  integer i_abs
  integer i4_log_10
  integer ten_pow

  if ( i == 0 ) then

    i4_log_10 = 0

  else

    i4_log_10 = 0
    ten_pow = 10

    i_abs = abs ( i )

    do while ( ten_pow <= i_abs )
      i4_log_10 = i4_log_10 + 1
      ten_pow = ten_pow * 10
    end do

  end if

  return
end
function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = rk ) value.
!
!    For now, the input quantity SEED is an integer variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer SEED, the "seed" value, 
!    which should NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = rk ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer k
  real ( kind = rk ) r8_uniform_01
  integer seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = rk ) * 4.656612875D-10

  return
end
subroutine r83s_cg ( n, a, b, x )

!*****************************************************************************80
!
!! R83S_CG uses the conjugate gradient method on an R83S system.
!
!  Discussion:
!
!    The R83S storage format is used for a tridiagonal scalar matrix.
!    The vector A(3) contains the subdiagonal, diagonal, and superdiagonal
!    values that occur on every row.
!    RGE A(I,J) = R83S A(I-J+2).
!
!    The matrix A must be a positive definite symmetric band matrix.
!
!    The method is designed to reach the solution after N computational
!    steps.  However, roundoff may introduce unacceptably large errors for
!    some problems.  In such a case, calling the routine again, using
!    the computed solution as the new starting estimate, should improve
!    the results.
!
!  Example:
!
!    Here is how an R83S matrix of order 5, stored as (A1,A2,A3), would
!    be interpreted:
!
!      A2  A1   0   0   0
!      A3  A2  A1   0   0
!       0  A3  A2  A1   0 
!       0   0  A3  A2  A1
!       0   0   0  A3  A2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Frank Beckman,
!    The Solution of Linear Equations by the Conjugate Gradient Method,
!    in Mathematical Methods for Digital Computers,
!    edited by John Ralston, Herbert Wilf,
!    Wiley, 1967,
!    ISBN: 0471706892,
!    LC: QA76.5.R3.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real ( kind = rk ) A(3), the matrix.
!
!    Input, real ( kind = rk ) B(N), the right hand side vector.
!
!    Input/output, real ( kind = rk ) X(N).
!    On input, an estimate for the solution, which may be 0.
!    On output, the approximate solution vector.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(3)
  real ( kind = rk ) alpha
  real ( kind = rk ) ap(n)
  real ( kind = rk ) b(n)
  real ( kind = rk ) beta
  integer it
  real ( kind = rk ) p(n)
  real ( kind = rk ) pap
  real ( kind = rk ) pr
  real ( kind = rk ) r(n)
  real ( kind = rk ) rap
  real ( kind = rk ) x(n)
!
!  Initialize
!    AP = A * x,
!    R  = b - A * x,
!    P  = b - A * x.
!
  call r83s_mv ( n, n, a, x, ap )

  r(1:n) = b(1:n) - ap(1:n)
  p(1:n) = b(1:n) - ap(1:n)
!
!  Do the N steps of the conjugate gradient method.
!
  do it = 1, n
!
!  Compute the matrix*vector product AP=A*P.
!
    call r83s_mv ( n, n, a, p, ap )
!
!  Compute the dot products
!    PAP = P*AP,
!    PR  = P*R
!  Set
!    ALPHA = PR / PAP.
!
    pap = dot_product ( p, ap )
    pr = dot_product ( p, r )

    if ( pap == 0.0D+00 ) then
      return
    end if

    alpha = pr / pap
!
!  Set
!    X = X + ALPHA * P
!    R = R - ALPHA * AP.
!
    x(1:n) = x(1:n) + alpha * p(1:n)
    r(1:n) = r(1:n) - alpha * ap(1:n)
!
!  Compute the vector dot product
!    RAP = R*AP
!  Set
!    BETA = - RAP / PAP.
!
    rap = dot_product ( r, ap )

    beta = - rap / pap
!
!  Update the perturbation vector
!    P = R + BETA * P.
!
    p(1:n) = r(1:n) + beta * p(1:n)

  end do

  return
end
subroutine r83s_dif2 ( m, n, a )

!*****************************************************************************80
!
!! R83S_DIF2 returns the DIF2 matrix in R83S format.
!
!  Discussion:
!
!    The R83S storage format is used for a tridiagonal scalar matrix.
!    The vector A(3) contains the subdiagonal, diagonal, and superdiagonal
!    values that occur on every row.
!    RGE A(I,J) = R83S A(I-J+2).
!
!  Example:
!
!    Here is how an R83S matrix of order 5, stored as (A1,A2,A3), would
!    be interpreted:
!
!      A2  A1   0   0   0
!      A3  A2  A1   0   0
!       0  A3  A2  A1   0 
!       0   0  A3  A2  A1
!       0   0   0  A3  A2
!
!  Properties:
!
!    A is banded, with bandwidth 3.
!    A is tridiagonal.
!    Because A is tridiagonal, it has property A (bipartite).
!    A is a special case of the TRIS or tridiagonal scalar matrix.
!    A is integral, therefore det ( A ) is integral, and 
!    det ( A ) * inverse ( A ) is integral.
!    A is Toeplitz: constant along diagonals.
!    A is symmetric: A' = A.
!    Because A is symmetric, it is normal.
!    Because A is normal, it is diagonalizable.
!    A is persymmetric: A(I,J) = A(N+1-J,N+1-I).
!    A is positive definite.
!    A is an M matrix.
!    A is weakly diagonally dominant, but not strictly diagonally dominant.
!    A has an LU factorization A = L * U, without pivoting.
!      The matrix L is lower bidiagonal with subdiagonal elements:
!        L(I+1,I) = -I/(I+1)
!      The matrix U is upper bidiagonal, with diagonal elements
!        U(I,I) = (I+1)/I
!      and superdiagonal elements which are all -1.
!    A has a Cholesky factorization A = L * L', with L lower bidiagonal.
!      L(I,I) =    sqrt ( (I+1) / I )
!      L(I,I-1) = -sqrt ( (I-1) / I )
!    The eigenvalues are
!      LAMBDA(I) = 2 + 2 * COS(I*PI/(N+1))
!                = 4 SIN^2(I*PI/(2*N+2))
!    The corresponding eigenvector X(I) has entries
!       X(I)(J) = sqrt(2/(N+1)) * sin ( I*J*PI/(N+1) ).
!    Simple linear systems:
!      x = (1,1,1,...,1,1),   A*x=(1,0,0,...,0,1)
!      x = (1,2,3,...,n-1,n), A*x=(0,0,0,...,0,n+1)
!    det ( A ) = N + 1.
!    The value of the determinant can be seen by induction,
!    and expanding the determinant across the first row:
!      det ( A(N) ) = 2 * det ( A(N-1) ) - (-1) * (-1) * det ( A(N-2) )
!                = 2 * N - (N-1)
!                = N + 1
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Gregory, David Karney,
!    A Collection of Matrices for Testing Computational Algorithms,
!    Wiley, 1969,
!    ISBN: 0882756494,
!    LC: QA263.68
!
!    Morris Newman, John Todd,
!    Example A8,
!    The evaluation of matrix inversion programs,
!    Journal of the Society for Industrial and Applied Mathematics,
!    Volume 6, Number 4, pages 466-476, 1958.
!
!    John Todd,
!    Basic Numerical Mathematics,
!    Volume 2: Numerical Algebra,
!    Birkhauser, 1980,
!    ISBN: 0817608117,
!    LC: QA297.T58.
!
!    Joan Westlake,
!    A Handbook of Numerical Matrix Inversion and Solution of 
!    Linear Equations,
!    John Wiley, 1968,
!    ISBN13: 978-0471936756,
!    LC: QA263.W47.
!
!  Parameters:
!
!    Input, integer M, N, the order of the matrix.
!
!    Output, real ( kind = rk ) A(3), the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(3)
  integer m
  integer n

  a(1) = -1.0D+00
  a(2) =  2.0D+00
  a(3) = -1.0D+00

  return
end
subroutine r83s_gs_sl ( n, a, b, x, tol, it_max, it, diff )

!*****************************************************************************80
!
!! R83S_GS_SL solves an R83S system using Gauss-Seidel iteration.
!
!  Discussion:
!
!    The R83S storage format is used for a tridiagonal scalar matrix.
!    The vector A(3) contains the subdiagonal, diagonal, and superdiagonal
!    values that occur on every row.
!    RGE A(I,J) = R83S A(I-J+2).
!
!    This routine simply applies a given number of steps of the
!    iteration to an input approximate solution.  On first call, you can
!    simply pass in the zero vector as an approximate solution.  If
!    the returned value is not acceptable, you may call again, using
!    it as the starting point for additional iterations.
!
!  Example:
!
!    Here is how an R83S matrix of order 5, stored as (A1,A2,A3), would
!    be interpreted:
!
!      A2  A1   0   0   0
!      A3  A2  A1   0   0
!       0  A3  A2  A1   0 
!       0   0  A3  A2  A1
!       0   0   0  A3  A2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real ( kind = rk ) A(3), the R83S matrix.
!
!    Input, real ( kind = rk ) B(N), the right hand side of the linear system.
!
!    Input/output, real ( kind = rk ) X(N), an approximate solution to 
!    the system.
!
!    Input, real ( kind = rk ) TOL, a tolerance.  If the maximum change in
!    the solution is less than TOL, the iteration is terminated early.
!
!    Input, integer IT_MAX, the maximum number of iterations.
!
!    Output, integer IT, the number of iterations taken.
!
!    Output, real ( kind = rk ) DIFF, the maximum change in the solution
!    on the last iteration.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(3)
  real ( kind = rk ) b(n)
  real ( kind = rk ) diff
  integer i
  integer it
  integer it_max
  integer it_num
  real ( kind = rk ) tol
  real ( kind = rk ) x(n)
  real ( kind = rk ) x_norm
  real ( kind = rk ) x_old(n)
!
!  No diagonal matrix entry can be zero.
!
  if ( a(2) == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R83S_GS_SL - Fatal error!'
    write ( *, '(a,i8)' ) '  Zero diagonal entry.'
    stop 1
  end if

  do it_num = 1, it_max

    it = it_num

    x_old(1:n) = x(1:n)

    x(1) =   ( b(1)                 - a(3) * x(2)   ) / a(2)
    do i = 2, n - 1
      x(i) = ( b(i) - a(1) * x(i-1) - a(3) * x(i+1) ) / a(2)
    end do
    x(n) =   ( b(n) - a(1) * x(n-1)                 ) / a(2)

    x_norm = maxval ( abs ( x(1:n) ) )
    diff = maxval ( abs ( x(1:n) - x_old(1:n) ) )

    if ( diff <= tol * ( x_norm + 1.0D+00 ) ) then
      exit
    end if

  end do

  return
end
subroutine r83s_indicator ( m, n, a )

!*****************************************************************************80
!
!! R83S_INDICATOR sets an R83S indicator matrix.
!
!  Discussion:
!
!    The R83S storage format is used for a tridiagonal scalar matrix.
!    The vector A(3) contains the subdiagonal, diagonal, and superdiagonal
!    values that occur on every row.
!    RGE A(I,J) = R83S A(I-J+2).
!
!  Example:
!
!    Here is how an R83S matrix of order 5, stored as (A1,A2,A3), would
!    be interpreted:
!
!      A2  A1   0   0   0
!      A3  A2  A1   0   0
!       0  A3  A2  A1   0 
!       0   0  A3  A2  A1
!       0   0   0  A3  A2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 September 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the order of the matrix.
!
!    Output, real ( kind = rk ) A(3), the R83S matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(3)
  integer m
  integer n

  a(1) = 3.0D+00
  a(2) = 2.0D+00
  a(3) = 1.0D+00

  return
end
subroutine r83s_jac_sl ( n, a, b, x, tol, it_max, it, diff )

!*****************************************************************************80
!
!! R83S_JAC_SL solves an R83S system using Jacobi iteration.
!
!  Discussion:
!
!    The R83S storage format is used for a tridiagonal scalar matrix.
!    The vector A(3) contains the subdiagonal, diagonal, and superdiagonal
!    values that occur on every row.
!    RGE A(I,J) = R83S A(I-J+2).
!
!    This routine simply applies a given number of steps of the
!    iteration to an input approximate solution.  On first call, you can
!    simply pass in the zero vector as an approximate solution.  If
!    the returned value is not acceptable, you may call again, using
!    it as the starting point for additional iterations.
!
!  Example:
!
!    Here is how an R83S matrix of order 5, stored as (A1,A2,A3), would
!    be interpreted:
!
!      A2  A1   0   0   0
!      A3  A2  A1   0   0
!       0  A3  A2  A1   0 
!       0   0  A3  A2  A1
!       0   0   0  A3  A2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real ( kind = rk ) A(3), the R83S matrix.
!
!    Input, real ( kind = rk ) B(N), the right hand side of the linear system.
!
!    Input/output, real ( kind = rk ) X(N), an approximate solution 
!    to the system.
!
!    Input, real ( kind = rk ) TOL, a tolerance.  If the maximum change in
!    the solution is less than TOL, the iteration is terminated early.
!
!    Input, integer IT_MAX, the maximum number of iterations.
!
!    Output, integer IT, the number of iterations taken.
!
!    Output, real ( kind = rk ) DIFF, the maximum change in the solution
!    on the last iteration.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(3)
  real ( kind = rk ) b(n)
  real ( kind = rk ) diff
  integer i
  integer it
  integer it_max
  integer it_num
  real ( kind = rk ) tol
  real ( kind = rk ) x(n)
  real ( kind = rk ) x_new(n)
  real ( kind = rk ) x_norm
!
!  No diagonal matrix entry can be zero.
!
  if ( a(2) == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R83S_JAC_SL - Fatal error!'
    write ( *, '(a,i8)' ) '  Zero diagonal entry.'
    stop 1
  end if

  do it_num = 1, it_max

    it = it_num

    x_new(1) =   b(1)                 - a(3) * x(2)
    do i = 2, n - 1
      x_new(i) = b(i) - a(1) * x(i-1) - a(3) * x(i+1)
    end do
    x_new(n) =   b(n) - a(1) * x(n-1)
!
!  Divide by diagonal terms.
!
    x_new(1:n) = x_new(1:n) / a(2)
!
!  Measure norms of solution and change:
!
    x_norm = maxval ( abs ( x(1:n) ) )
    diff = maxval ( abs ( x_new(1:n) - x(1:n) ) )
!
!  Update.
!
    x(1:n) = x_new(1:n)
!
!  Test for early termination.
!
    if ( diff <= tol * ( x_norm + 1.0D+00 ) ) then
      exit
    end if

  end do

  return
end
subroutine r83s_mtv ( m, n, a, x, b )

!*****************************************************************************80
!
!! R83S_MTV computes b=A'*x, where A is an R83S matrix.
!
!  Discussion:
!
!    The R83S storage format is used for a tridiagonal scalar matrix.
!    The vector A(3) contains the subdiagonal, diagonal, and superdiagonal
!    values that occur on every row.
!    RGE A(I,J) = R83S A(I-J+2).
!
!  Example:
!
!    Here is how an R83S matrix of order 5, stored as (A1,A2,A3), would
!    be interpreted:
!
!      A2  A1   0   0   0
!      A3  A2  A1   0   0
!       0  A3  A2  A1   0 
!       0   0  A3  A2  A1
!       0   0   0  A3  A2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the order of the matrix.
!
!    Input, real ( kind = rk ) A(3), the R83S matrix.
!
!    Input, real ( kind = rk ) X(M), the vector to be multiplied by A'.
!
!    Output, real ( kind = rk ) B(N), the product A' * x.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(3)
  real ( kind = rk ) b(n)
  integer i
  integer j
  real ( kind = rk ) x(m)

  b(1:n) = 0.0D+00

  do j = 1, n
    do i = max ( 1, j - 1 ), min ( m, j + 1 )
      b(j) = b(j) + a(i-j+2) * x(i)
    end do
  end do

  return
end
subroutine r83s_mv ( m, n, a, x, b )

!*****************************************************************************80
!
!! R83S_MV multiplies an R83S matrix times an R8VEC.
!
!  Discussion:
!
!    The R83S storage format is used for a tridiagonal scalar matrix.
!    The vector A(3) contains the subdiagonal, diagonal, and superdiagonal
!    values that occur on every row.
!    RGE A(I,J) = R83S A(I-J+2).
!
!  Example:
!
!    Here is how an R83S matrix of order 5, stored as (A1,A2,A3), would
!    be interpreted:
!
!      A2  A1   0   0   0
!      A3  A2  A1   0   0
!       0  A3  A2  A1   0 
!       0   0  A3  A2  A1
!       0   0   0  A3  A2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the order of the matrix.
!
!    Input, real ( kind = rk ) A(3), the R83S matrix.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = rk ) B(M), the product A * x.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(3)
  real ( kind = rk ) b(m)
  integer i
  integer j
  real ( kind = rk ) x(n)

  b(1:m) = 0.0D+00

  do j = 1, n
    do i = max ( 1, j - 1 ), min ( m, j + 1 )
      b(i) = b(i) + a(i-j+2) * x(j)
    end do
  end do

  return
end
subroutine r83s_print ( m, n, a, title )

!*****************************************************************************80
!
!! R83S_PRINT prints an R83S matrix.
!
!  Discussion:
!
!    The R83S storage format is used for a tridiagonal scalar matrix.
!    The vector A(3) contains the subdiagonal, diagonal, and superdiagonal
!    values that occur on every row.
!    RGE A(I,J) = R83S A(I-J+2).
!
!  Example:
!
!    Here is how an R83S matrix of order 5, stored as (A1,A2,A3), would
!    be interpreted:
!
!      A2  A1   0   0   0
!      A3  A2  A1   0   0
!       0  A3  A2  A1   0 
!       0   0  A3  A2  A1
!       0   0   0  A3  A2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 September 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the order of the matrix.
!
!    Input, real ( kind = rk ) A(3), the R83S matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(3)
  integer m
  integer n
  character ( len = * ) title

  call r83s_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r83s_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R83S_PRINT_SOME prints some of an R83S matrix.
!
!  Discussion:
!
!    The R83S storage format is used for a tridiagonal scalar matrix.
!    The vector A(3) contains the subdiagonal, diagonal, and superdiagonal
!    values that occur on every row.
!    RGE A(I,J) = R83S A(I-J+2).
!
!  Example:
!
!    Here is how an R83S matrix of order 5, stored as (A1,A2,A3), would
!    be interpreted:
!
!      A2  A1   0   0   0
!      A3  A2  A1   0   0
!       0  A3  A2  A1   0 
!       0   0  A3  A2  A1
!       0   0   0  A3  A2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 September 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the order of the matrix.
!
!    Input, real ( kind = rk ) A(3), the R83S matrix.
!
!    Input, integer ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column, to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: incx = 5

  real ( kind = rk ) a(3)
  character ( len = 14 ) ctemp(incx)
  integer i
  integer i2hi
  integer i2lo
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  integer m
  integer n
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(a,5a14)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2lo = max ( i2lo, j2lo - 1 )
    i2hi = min ( ihi, m )
    i2hi = min ( i2hi, j2hi + 1 )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( i - j + 2 < 1 .or. 3 < i - j + 2 ) then
          ctemp(j2) = '              '
        else
          write ( ctemp(j2), '(g14.6)' ) a(i-j+2)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r83s_random ( m, n, seed, a )

!*****************************************************************************80
!
!! R83S_RANDOM randomizes an R83S matrix.
!
!  Discussion:
!
!    The R83S storage format is used for a tridiagonal scalar matrix.
!    The vector A(3) contains the subdiagonal, diagonal, and superdiagonal
!    values that occur on every row.
!    RGE A(I,J) = R83S A(I-J+2).
!
!  Example:
!
!    Here is how an R83S matrix of order 5, stored as (A1,A2,A3), would
!    be interpreted:
!
!      A2  A1   0   0   0
!      A3  A2  A1   0   0
!       0  A3  A2  A1   0 
!       0   0  A3  A2  A1
!       0   0   0  A3  A2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 September 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the order of the linear system.
!
!    Input/output, integer SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = rk ) A(3), the R83S matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(3)
  integer m
  integer n
  integer seed

  call r8vec_uniform_01 ( 3, seed, a )

  return
end
subroutine r83s_res ( m, n, a, x, b, r )

!*****************************************************************************80
!
!! R83S_RES computes the residual R = B-A*X for R83S matrices.
!
!  Discussion:
!
!    The R83S storage format is used for a tridiagonal scalar matrix.
!    The vector A(3) contains the subdiagonal, diagonal, and superdiagonal
!    values that occur on every row.
!    RGE A(I,J) = R83S A(I-J+2).
!
!  Example:
!
!    Here is how an R83S matrix of order 5, stored as (A1,A2,A3), would
!    be interpreted:
!
!      A2  A1   0   0   0
!      A3  A2  A1   0   0
!       0  A3  A2  A1   0 
!       0   0  A3  A2  A1
!       0   0   0  A3  A2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the order of the matrix.
!
!    Input, real ( kind = rk ) A(3), the matrix.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied by A.
!
!    Input, real ( kind = rk ) B(M), the desired result A * x.
!
!    Output, real ( kind = rk ) R(M), the residual R = B - A * X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(3)
  real ( kind = rk ) b(m)
  real ( kind = rk ) r(m)
  real ( kind = rk ) x(n)

  call r83s_mv ( m, n, a, x, r )

  r(1:m) = b(1:m) - r(1:m)

  return
end
subroutine r83s_to_r8ge ( m, n, a_83s, a_ge )

!*****************************************************************************80
!
!! R83S_TO_R8GE copies an R83S matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R83S storage format is used for a tridiagonal scalar matrix.
!    The vector A(3) contains the subdiagonal, diagonal, and superdiagonal
!    values that occur on every row.
!    RGE A(I,J) = R83S A(I-J+2).
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!  Example:
!
!    Here is how an R83S matrix of order 5, stored as (A1,A2,A3), would
!    be interpreted:
!
!      A2  A1   0   0   0
!      A3  A2  A1   0   0
!       0  A3  A2  A1   0 
!       0   0  A3  A2  A1
!       0   0   0  A3  A2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 September 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the order of the matrix.
!
!    Input, real ( kind = rk ) A_83S(3), the R83S matrix.
!
!    Output, real ( kind = rk ) A_GE(M,N), the R8GE matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a_83s(3)
  real ( kind = rk ) a_ge(m,n)
  integer i
  integer j

  a_ge(1:m,1:n) = 0.0D+00

  do j = 1, n
    do i = max ( 1, j - 1 ), min ( m, j + 1 )
      a_ge(i,j) = a_83s(i-j+2)
    end do
  end do

  return
end
subroutine r83s_zeros ( m, n, a )

!*****************************************************************************80
!
!! R83S_ZEROS zeros an R83S matrix.
!
!  Discussion:
!
!    The R83S storage format is used for a tridiagonal scalar matrix.
!    The vector A(3) contains the subdiagonal, diagonal, and superdiagonal
!    values that occur on every row.
!    RGE A(I,J) = R83S A(I-J+2).
!
!  Example:
!
!    Here is how an R83S matrix of order 5, stored as (A1,A2,A3), would
!    be interpreted:
!
!      A2  A1   0   0   0
!      A3  A2  A1   0   0
!       0  A3  A2  A1   0 
!       0   0  A3  A2  A1
!       0   0   0  A3  A2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 September 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the order of the matrix.
!
!    Output, real ( kind = rk ) A(3), the R83S matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(3)
  integer m
  integer n

  a(1:3) = 0.0D+00

  return
end
subroutine r8ge_mtv ( m, n, a, x, b )

!*****************************************************************************80
!
!! R8GE_MTV multiplies an R8VEC by an R8GE matrix.
!
!  Discussion:
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!    R8GE storage is used by LINPACK and LAPACK.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, real ( kind = rk ) A(M,N), the R8GE matrix.
!
!    Input, real ( kind = rk ) X(M), the vector to be multiplied by A.
!
!    Output, real ( kind = rk ) B(N), the product A' * x.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) b(n)
  real ( kind = rk ) x(m)

  b(1:n) = matmul ( transpose ( a(1:m,1:n) ), x(1:m) )

  return
end
subroutine r8ge_mv ( m, n, a, x, b )

!*****************************************************************************80
!
!! R8GE_MV multiplies an R8GE matrix by an R8VEC.
!
!  Discussion:
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!    R8GE storage is used by LINPACK and LAPACK.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, real ( kind = rk ) A(M,N), the R8GE matrix.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = rk ) B(M), the product A * x.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) b(m)
  real ( kind = rk ) x(n)

  b(1:m) = matmul ( a(1:m,1:n), x(1:n) )

  return
end
subroutine r8ge_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8GE_PRINT prints an R8GE matrix.
!
!  Discussion:
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!    R8GE storage is used by LINPACK and LAPACK.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, real ( kind = rk ) A(M,N), the R8GE matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  character ( len = * ) title

  call r8ge_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8ge_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8GE_PRINT_SOME prints some of an R8GE matrix.
!
!  Discussion:
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!    R8GE storage is used by LINPACK and LAPACK.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, real ( kind = rk ) A(M,N), the R8GE matrix.
!
!    Input, integer ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: incx = 5
  integer m
  integer n

  real ( kind = rk ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer i
  integer i2hi
  integer i2lo
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(''  Col:  '',5a14)' ) ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        write ( ctemp(j2), '(g14.6)' ) a(i,j)

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8vec_indicator1 ( n, a )

!*****************************************************************************80
!
!! R8VEC_INDICATOR1 sets an R8VEC to the indicator1 vector.
!
!  Discussion:
!
!    A(1:N) = (/ 1 : N /)
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements of A.
!
!    Output, real ( kind = rk ) A(N), the array to be initialized.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  integer i

  do i = 1, n
    a(i) = real ( i, kind = rk )
  end do

  return
end
function r8vec_norm ( n, a )

!*****************************************************************************80
!
!! R8VEC_NORM returns the L2 norm of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The vector L2 norm is defined as:
!
!      R8VEC_NORM = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in A.
!
!    Input, real ( kind = rk ) A(N), the vector whose L2 norm is desired.
!
!    Output, real ( kind = rk ) R8VEC_NORM, the L2 norm of A.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) r8vec_norm

  r8vec_norm = sqrt ( sum ( a(1:n)**2 ) )

  return
end
function r8vec_norm_affine ( n, v0, v1 )

!*****************************************************************************80
!
!! R8VEC_NORM_AFFINE returns the affine norm of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The affine vector L2 norm is defined as:
!
!      R8VEC_NORM_AFFINE(V0,V1)
!        = sqrt ( sum ( 1 <= I <= N ) ( V1(I) - V0(I) )^2
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the vectors.
!
!    Input, real ( kind = rk ) V0(N), the base vector.
!
!    Input, real ( kind = rk ) V1(N), the vector whose affine norm is desired.
!
!    Output, real ( kind = rk ) R8VEC_NORM_AFFINE, the L2 norm of V1-V0.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) r8vec_norm_affine
  real ( kind = rk ) v0(n)
  real ( kind = rk ) v1(n)

  r8vec_norm_affine = sqrt ( sum ( ( v0(1:n) - v1(1:n) )**2 ) )

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, real ( kind = rk ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  integer i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i8,g14.6)' ) i, a(i)
  end do

  return
end
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of real ( kind = rk ) values.
!
!    For now, the input quantity SEED is an integer variable.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer N, the number of entries 
!    in the vector.
!
!    Input/output, integer SEED, the "seed" value, 
!    which should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = rk ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer k
  integer seed
  real ( kind = rk ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + 2147483647
    end if

    r(i) = real ( seed, kind = rk ) * 4.656612875D-10

  end do

  return
end
subroutine r8vec2_print ( n, a1, a2, title )

!*****************************************************************************80
!
!! R8VEC2_PRINT prints an R8VEC2.
!
!  Discussion:
!
!    An R8VEC2 is a dataset consisting of N pairs of R8's, stored
!    as two separate vectors A1 and A2.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, real ( kind = rk ) A1(N), A2(N), the vectors to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a1(n)
  real ( kind = rk ) a2(n)
  integer i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i4,2x,g14.6,2x,g14.6)' ) i, a1(i), a2(i)
  end do

  return
end

