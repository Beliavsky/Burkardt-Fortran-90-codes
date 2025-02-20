function i4_log_10 ( i )

!*****************************************************************************80
!
!! i4_log_10() returns the integer part of the logarithm base 10 of an I4.
!
!  Discussion:
!
!    I4_LOG_10 ( I ) + 1 is the number of decimal digits in I.
!
!    An I4 is an integer value.
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
!    Input/output, integer SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = rk ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: i4_huge = 2147483647
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
    seed = seed + i4_huge
  end if

  r8_uniform_01 = real ( seed, kind = rk ) * 4.656612875D-10

  return
end
subroutine r83t_cg ( n, a, b, x )

!*****************************************************************************80
!
!! R83T_CG uses the conjugate gradient method on an R83T system.
!
!  Discussion:
!
!    The R83T storage format is used for an MxN tridiagonal matrix.
!    The superdiagonal is stored in entries (1:M-1,3), the diagonal in
!    entries (1:M,2), and the subdiagonal in (2:M,1).  Thus, the
!    the rows of the original matrix slide horizontally to form an
!    Mx3 stack of data.
!
!    An R83T matrix of order 3x5 would be stored:
!
!       *  A11 A12
!      A21 A22 A23
!      A32 A33 A34
!
!    An R83T matrix of order 5x5 would be stored:
!
!       *  A11 A12
!      A21 A22 A23
!      A32 A33 A34
!      A43 A44 A45
!      A54 A55  *
!
!    An R83T matrix of order 5x3 would be stored:
!
!       *  A11 A12
!      A21 A22 A23
!      A32 A33  *
!      A43  *   *
!       *   *   *
!
!    The matrix A must be a positive definite symmetric band matrix.
!
!    The method is designed to reach the solution after N computational
!    steps.  However, roundoff may introduce unacceptably large errors for
!    some problems.  In such a case, calling the routine again, using
!    the computed solution as the new starting estimate, should improve
!    the results.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 June 2014
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
!    N must be positive.
!
!    Input, real ( kind = rk ) A(N,3), the matrix.
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

  real ( kind = rk ) a(n,3)
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
  call r83t_mv ( n, n, a, x, ap )

  r(1:n) = b(1:n) - ap(1:n)
  p(1:n) = b(1:n) - ap(1:n)
!
!  Do the N steps of the conjugate gradient method.
!
  do it = 1, n
!
!  Compute the matrix*vector product AP=A*P.
!
    call r83t_mv ( n, n, a, p, ap )
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
subroutine r83t_dif2 ( m, n, a )

!*****************************************************************************80
!
!! R83T_DIF2 returns the DIF2 matrix in R83T format.
!
!  Discussion:
!
!    The R83T storage format is used for an MxN tridiagonal matrix.
!    The superdiagonal is stored in entries (1:M-1,3), the diagonal in
!    entries (1:M,2), and the subdiagonal in (2:M,1).  Thus, the
!    the rows of the original matrix slide horizontally to form an
!    Mx3 stack of data.
!
!    An R83T matrix of order 3x5 would be stored:
!
!       *  A11 A12
!      A21 A22 A23
!      A32 A33 A34
!
!    An R83T matrix of order 5x5 would be stored:
!
!       *  A11 A12
!      A21 A22 A23
!      A32 A33 A34
!      A43 A44 A45
!      A54 A55  *
!
!    An R83T matrix of order 5x3 would be stored:
!
!       *  A11 A12
!      A21 A22 A23
!      A32 A33  *
!      A43  *   *
!       *   *   *
!
!  Example:
!
!    N = 5
!
!    2 -1  .  .  .
!   -1  2 -1  .  .
!    . -1  2 -1  .
!    .  . -1  2 -1
!    .  .  . -1  2
!
!  Properties:
!
!    A is banded, with bandwidth 3.
!
!    A is tridiagonal.
!
!    Because A is tridiagonal, it has property A (bipartite).
!
!    A is a special case of the TRIS or tridiagonal scalar matrix.
!
!    A is integral, therefore det ( A ) is integral, and 
!    det ( A ) * inverse ( A ) is integral.
!
!    A is Toeplitz: constant along diagonals.
!
!    A is symmetric: A' = A.
!
!    Because A is symmetric, it is normal.
!
!    Because A is normal, it is diagonalizable.
!
!    A is persymmetric: A(I,J) = A(N+1-J,N+1-I).
!
!    A is positive definite.
!
!    A is an M matrix.
!
!    A is weakly diagonally dominant, but not strictly diagonally dominant.
!
!    A has an LU factorization A = L * U, without pivoting.
!
!      The matrix L is lower bidiagonal with subdiagonal elements:
!
!        L(I+1,I) = -I/(I+1)
!
!      The matrix U is upper bidiagonal, with diagonal elements
!
!        U(I,I) = (I+1)/I
!
!      and superdiagonal elements which are all -1.
!
!    A has a Cholesky factorization A = L * L', with L lower bidiagonal.
!
!      L(I,I) =    sqrt ( (I+1) / I )
!      L(I,I-1) = -sqrt ( (I-1) / I )
!
!    The eigenvalues are
!
!      LAMBDA(I) = 2 + 2 * COS(I*PI/(N+1))
!                = 4 SIN^2(I*PI/(2*N+2))
!
!    The corresponding eigenvector X(I) has entries
!
!       X(I)(J) = sqrt(2/(N+1)) * sin ( I*J*PI/(N+1) ).
!
!    Simple linear systems:
!
!      x = (1,1,1,...,1,1),   A*x=(1,0,0,...,0,1)
!
!      x = (1,2,3,...,n-1,n), A*x=(0,0,0,...,0,n+1)
!
!    det ( A ) = N + 1.
!
!    The value of the determinant can be seen by induction,
!    and expanding the determinant across the first row:
!
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
!    18 June 2014
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
!    Output, real ( kind = rk ) A(M,3), the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,3)
  integer mn

  a(1:m,1:3) = 0.0D+00

  mn = min ( m, n )

  a(2:mn,1)   = -1.0D+00
  a(1:mn,2)   =  2.0D+00
  a(1:mn-1,3) = -1.0D+00

  if ( m < n ) then
    a(mn,3) = -1.0D+00
  else if ( n < m ) then
    a(mn+1,1) = -1.0D+00
  end if
  
  return
end
subroutine r83t_gs_sl ( n, a, b, x, it_max )

!*****************************************************************************80
!
!! R83T_GS_SL solves an R83T system using Gauss-Seidel iteration.
!
!  Discussion:
!
!    The R83T storage format is used for an MxN tridiagonal matrix.
!    The superdiagonal is stored in entries (1:M-1,3), the diagonal in
!    entries (1:M,2), and the subdiagonal in (2:M,1).  Thus, the
!    the rows of the original matrix slide horizontally to form an
!    Mx3 stack of data.
!
!    An R83T matrix of order 3x5 would be stored:
!
!       *  A11 A12
!      A21 A22 A23
!      A32 A33 A34
!
!    An R83T matrix of order 5x5 would be stored:
!
!       *  A11 A12
!      A21 A22 A23
!      A32 A33 A34
!      A43 A44 A45
!      A54 A55  *
!
!    An R83T matrix of order 5x3 would be stored:
!
!       *  A11 A12
!      A21 A22 A23
!      A32 A33  *
!      A43  *   *
!       *   *   *
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 May 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 2.
!
!    Input, real ( kind = rk ) A(N,3), the R83T matrix.
!
!    Input, real ( kind = rk ) B(N), the right hand side of the linear system.
!
!    Input/output, real ( kind = rk ) X(N), an approximate solution to 
!    the system.
!
!    Input, integer IT_MAX, the maximum number of iterations.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,3)
  real ( kind = rk ) b(n)
  integer i
  integer it_max
  integer it_num
  real ( kind = rk ) x(n)
!
!  No diagonal matrix entry can be zero.
!
  do i = 1, n
    if ( a(i,2) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R83_GS_SL - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero diagonal entry, index = ', i
      stop 1
    end if
  end do

  do it_num = 1, it_max

    x(1) = ( b(1) - a(1,3) * x(2) ) / a(1,2)

    do i = 2, n - 1
      x(i) = ( b(i) - a(i,1) * x(i-1) - a(i,3) * x(i+1) ) / a(i,2)
    end do
    x(n) = ( b(n) - a(n,1) * x(n-1) ) / a(n,2)

  end do

  return
end
subroutine r83t_indicator ( m, n, a )

!*****************************************************************************80
!
!! R83T_INDICATOR sets the indicator matrix in R83T format.
!
!  Discussion:
!
!    The R83T storage format is used for an MxN tridiagonal matrix.
!    The superdiagonal is stored in entries (1:M-1,3), the diagonal in
!    entries (1:M,2), and the subdiagonal in (2:M,1).  Thus, the
!    the rows of the original matrix slide horizontally to form an
!    Mx3 stack of data.
!
!    An R83T matrix of order 3x5 would be stored:
!
!       *  A11 A12
!      A21 A22 A23
!      A32 A33 A34
!
!    An R83T matrix of order 5x5 would be stored:
!
!       *  A11 A12
!      A21 A22 A23
!      A32 A33 A34
!      A43 A44 A45
!      A54 A55  *
!
!    An R83T matrix of order 5x3 would be stored:
!
!       *  A11 A12
!      A21 A22 A23
!      A32 A33  *
!      A43  *   *
!       *   *   *
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 May 2016
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
!    Output, real ( kind = rk ) A(M,3), the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,3)
  real ( kind = rk ) fac
  integer i
  integer i4_log_10
  integer j
  integer k

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  do i = 1, m
    do k = 1, 3
      j = i - 2 + k
      if ( 1 <= j .and. j <= n ) then
        a(i,k) = real ( fac * i + j, kind = rk )
      else
        a(i,k) = 0.0D+00
      end if
    end do
  end do

  return
end
subroutine r83t_jac_sl ( n, a, b, x, it_max )

!*****************************************************************************80
!
!! R83T_JAC_SL solves an R83T system using Jacobi iteration.
!
!  Discussion:
!
!    The R83T storage format is used for an MxN tridiagonal matrix.
!    The superdiagonal is stored in entries (1:M-1,3), the diagonal in
!    entries (1:M,2), and the subdiagonal in (2:M,1).  Thus, the
!    the rows of the original matrix slide horizontally to form an
!    Mx3 stack of data.
!
!    An R83T matrix of order 3x5 would be stored:
!
!       *  A11 A12
!      A21 A22 A23
!      A32 A33 A34
!
!    An R83T matrix of order 5x5 would be stored:
!
!       *  A11 A12
!      A21 A22 A23
!      A32 A33 A34
!      A43 A44 A45
!      A54 A55  *
!
!    An R83T matrix of order 5x3 would be stored:
!
!       *  A11 A12
!      A21 A22 A23
!      A32 A33  *
!      A43  *   *
!       *   *   *
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 May 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 2.
!
!    Input, real ( kind = rk ) A(N,3), the R83T matrix.
!
!    Input, real ( kind = rk ) B(N), the right hand side of the linear system.
!
!    Input/output, real ( kind = rk ) X(N), an approximate solution 
!    to the system.
!
!    Input, integer IT_MAX, the maximum number of iterations.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,3)
  real ( kind = rk ) b(n)
  integer i
  integer it_max
  integer it_num
  real ( kind = rk ) x(n)
  real ( kind = rk ) x_new(n)
!
!  No diagonal matrix entry can be zero.
!
  do i = 1, n
    if ( a(i,2) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R83T_JAC_SL - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero diagonal entry, index = ', i
      stop 1
    end if
  end do

  do it_num = 1, it_max

    x_new(1) = b(1) - a(1,3) * x(2)
    do i = 2, n - 1
      x_new(i) = b(i) - a(i,1) * x(i-1) - a(i,3) * x(i+1)
    end do
    x_new(n) = b(n) - a(n,1) * x(n-1)
!
!  Divide by diagonal terms.
!
    x_new(1:n) = x_new(1:n) / a(1:n,2)
!
!  Update.
!
    x(1:n) = x_new(1:n)

  end do

  return
end
subroutine r83t_mtv ( m, n, a, x, b )

!*****************************************************************************80
!
!! R83T_MTV multiplies an R83T matrix transposed times an R8VEC.
!
!  Discussion:
!
!    The R83T storage format is used for an MxN tridiagonal matrix.
!    The superdiagonal is stored in entries (1:M-1,3), the diagonal in
!    entries (1:M,2), and the subdiagonal in (2:M,1).  Thus, the
!    the rows of the original matrix slide horizontally to form an
!    Mx3 stack of data.
!
!    An R83T matrix of order 3x5 would be stored:
!
!       *  A11 A12
!      A21 A22 A23
!      A32 A33 A34
!
!    An R83T matrix of order 5x5 would be stored:
!
!       *  A11 A12
!      A21 A22 A23
!      A32 A33 A34
!      A43 A44 A45
!      A54 A55  *
!
!    An R83T matrix of order 5x3 would be stored:
!
!       *  A11 A12
!      A21 A22 A23
!      A32 A33  *
!      A43  *   *
!       *   *   *
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 May 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = rk ) A(M,3), the matrix.
!
!    Input, real ( kind = rk ) X(M), the vector to be multiplied by A.
!
!    Output, real ( kind = rk ) B(N), the product A' * x.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,3)
  real ( kind = rk ) b(n)
  integer i
  integer j
  integer k
  real ( kind = rk ) x(m)

  b(1:n) = 0.0D+00

  do i = 1, m
    do k = 1, 3
      j = i - 2 + k
      if ( 1 <= j .and. j <= n ) then
        b(j) = b(j) + x(i) * a(i,k)
      end if
    end do
  end do

  return
end
subroutine r83t_mv ( m, n, a, x, b )

!*****************************************************************************80
!
!! R83T_MV multiplies an R83T matrix times an R8VEC.
!
!  Discussion:
!
!    The R83T storage format is used for an MxN tridiagonal matrix.
!    The superdiagonal is stored in entries (1:M-1,3), the diagonal in
!    entries (1:M,2), and the subdiagonal in (2:M,1).  Thus, the
!    the rows of the original matrix slide horizontally to form an
!    Mx3 stack of data.
!
!    An R83T matrix of order 3x5 would be stored:
!
!       *  A11 A12
!      A21 A22 A23
!      A32 A33 A34
!
!    An R83T matrix of order 5x5 would be stored:
!
!       *  A11 A12
!      A21 A22 A23
!      A32 A33 A34
!      A43 A44 A45
!      A54 A55  *
!
!    An R83T matrix of order 5x3 would be stored:
!
!       *  A11 A12
!      A21 A22 A23
!      A32 A33  *
!      A43  *   *
!       *   *   *
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 June 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = rk ) A(M,3), the matrix.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = rk ) B(M), the product A * x.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,3)
  real ( kind = rk ) b(m)
  integer mn
  real ( kind = rk ) x(n)

  b(1:m) = 0.0D+00

  mn = min ( m, n )

  if ( n == 1 ) then
    b(1) = a(1,2) * x(1)
    if ( 1 < m ) then
      b(2) = a(2,1) * x(1)
    end if
    return
  end if

  b(1)      = a(1,2)      * x(1) &
            + a(1,3)      * x(2)

  b(2:mn-1) = a(2:mn-1,1) * x(1:mn-2) &
            + a(2:mn-1,2) * x(2:mn-1) &
            + a(2:mn-1,3) * x(3:mn)

  b(mn)     = a(mn,1)     * x(mn-1) &
            + a(mn,2)     * x(mn)

  if ( n < m ) then
    b(mn+1) = b(mn+1) + a(mn+1,1) * x(mn)
  else if ( m < n ) then
    b(mn) = b(mn) + a(mn,3) * x(mn+1)
  end if

  return
end
subroutine r83t_print ( m, n, a, title )

!*****************************************************************************80
!
!! R83T_PRINT prints an R83T matrix.
!
!  Discussion:
!
!    The R83T storage format is used for an MxN tridiagonal matrix.
!    The superdiagonal is stored in entries (1:M-1,3), the diagonal in
!    entries (1:M,2), and the subdiagonal in (2:M,1).  Thus, the
!    the rows of the original matrix slide horizontally to form an
!    Mx3 stack of data.
!
!    An R83T matrix of order 3x5 would be stored:
!
!       *  A11 A12
!      A21 A22 A23
!      A32 A33 A34
!
!    An R83T matrix of order 5x5 would be stored:
!
!       *  A11 A12
!      A21 A22 A23
!      A32 A33 A34
!      A43 A44 A45
!      A54 A55  *
!
!    An R83T matrix of order 5x3 would be stored:
!
!       *  A11 A12
!      A21 A22 A23
!      A32 A33  *
!      A43  *   *
!       *   *   *
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 May 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the order of the matrix.
!    N must be positive.
!
!    Input, real ( kind = rk ) A(M,3), the R83T matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m

  real ( kind = rk ) a(m,3)
  integer n
  character ( len = * ) title

  call r83t_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r83t_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R83T_PRINT_SOME prints some of an R83T matrix.
!
!  Discussion:
!
!  Discussion:
!
!    The R83T storage format is used for an MxN tridiagonal matrix.
!    The superdiagonal is stored in entries (1:M-1,3), the diagonal in
!    entries (1:M,2), and the subdiagonal in (2:M,1).  Thus, the
!    the rows of the original matrix slide horizontally to form an
!    Mx3 stack of data.
!
!    An R83T matrix of order 3x5 would be stored:
!
!       *  A11 A12
!      A21 A22 A23
!      A32 A33 A34
!
!    An R83T matrix of order 5x5 would be stored:
!
!       *  A11 A12
!      A21 A22 A23
!      A32 A33 A34
!      A43 A44 A45
!      A54 A55  *
!
!    An R83T matrix of order 5x3 would be stored:
!
!       *  A11 A12
!      A21 A22 A23
!      A32 A33  *
!      A43  *   *
!       *   *   *
!
!  Example:
!
!    An R83 matrix of order 3x5 would be stored:
!
!       *  A12 A23 A34  *
!      A11 A22 A33  *   *
!      A21 A32  *   *   *
!
!    An R83 matrix of order 5x5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!    An R83 matrix of order 5x3 would be stored:
!
!       *  A12 A23
!      A11 A22 A33
!      A21 A32 A43
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 May 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the order of the matrix.
!    N must be positive.
!
!    Input, real ( kind = rk ) A(M,3), the R83T matrix.
!
!    Input, integer ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column, to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: incx = 5
  integer m

  real ( kind = rk ) a(m,3)
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
  integer n
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if
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
          write ( ctemp(j2), '(g14.6)' ) a(i,j-i+2)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r83t_random ( m, n, seed, a )

!*****************************************************************************80
!
!! R83T_RANDOM returns a random R83T matrix.
!
!  Discussion:
!
!    The R83T storage format is used for an MxN tridiagonal matrix.
!    The superdiagonal is stored in entries (1:M-1,3), the diagonal in
!    entries (1:M,2), and the subdiagonal in (2:M,1).  Thus, the
!    the rows of the original matrix slide horizontally to form an
!    Mx3 stack of data.
!
!    An R83T matrix of order 3x5 would be stored:
!
!       *  A11 A12
!      A21 A22 A23
!      A32 A33 A34
!
!    An R83T matrix of order 5x5 would be stored:
!
!       *  A11 A12
!      A21 A22 A23
!      A32 A33 A34
!      A43 A44 A45
!      A54 A55  *
!
!    An R83T matrix of order 5x3 would be stored:
!
!       *  A11 A12
!      A21 A22 A23
!      A32 A33  *
!      A43  *   *
!       *   *   *
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 August 2015
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
!    Input/output, integer SEED, a seed for the random number
!    generator.
!
!    Output, real ( kind = rk ) A(M,3), the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,3)
  integer i
  integer j
  integer k
  real ( kind = rk ) r8_uniform_01
  integer seed

  do i = 1, m
    do k = 1, 3
      j = i - 2 + k
      if ( 1 <= j .and. j <= n ) then
        a(i,k) = r8_uniform_01 ( seed )
      else
        a(i,k) = 0.0D+00
      end if
    end do
  end do

  return
end
subroutine r83t_res ( m, n, a, x, b, r )

!*****************************************************************************80
!
!! R83T_RES computes the residual R = B-A*X for R83T matrices.
!
!  Discussion:
!
!    The R83T storage format is used for an MxN tridiagonal matrix.
!    The superdiagonal is stored in entries (1:M-1,3), the diagonal in
!    entries (1:M,2), and the subdiagonal in (2:M,1).  Thus, the
!    the rows of the original matrix slide horizontally to form an
!    Mx3 stack of data.
!
!    An R83T matrix of order 3x5 would be stored:
!
!       *  A11 A12
!      A21 A22 A23
!      A32 A33 A34
!
!    An R83T matrix of order 5x5 would be stored:
!
!       *  A11 A12
!      A21 A22 A23
!      A32 A33 A34
!      A43 A44 A45
!      A54 A55  *
!
!    An R83T matrix of order 5x3 would be stored:
!
!       *  A11 A12
!      A21 A22 A23
!      A32 A33  *
!      A43  *   *
!       *   *   *
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 June 2014
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
!    Input, real ( kind = rk ) A(M,3), the matrix.
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

  real ( kind = rk ) a(m,3)
  real ( kind = rk ) b(m)
  real ( kind = rk ) r(m)
  real ( kind = rk ) x(n)

  call r83t_mv ( m, n, a, x, r )

  r(1:m) = b(1:m) - r(1:m)

  return
end
subroutine r83t_to_r8ge ( m, n, a_r83t, a_r8ge )

!*****************************************************************************80
!
!! R83T_TO_R8GE copies an R83T matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R83T storage format is used for an MxN tridiagonal matrix.
!    The superdiagonal is stored in entries (1:M-1,3), the diagonal in
!    entries (1:M,2), and the subdiagonal in (2:M,1).  Thus, the
!    the rows of the original matrix slide horizontally to form an
!    Mx3 stack of data.
!
!    An R83T matrix of order 3x5 would be stored:
!
!       *  A11 A12
!      A21 A22 A23
!      A32 A33 A34
!
!    An R83T matrix of order 5x5 would be stored:
!
!       *  A11 A12
!      A21 A22 A23
!      A32 A33 A34
!      A43 A44 A45
!      A54 A55  *
!
!    An R83T matrix of order 5x3 would be stored:
!
!       *  A11 A12
!      A21 A22 A23
!      A32 A33  *
!      A43  *   *
!       *   *   *
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 May 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the order of the matrix.
!    N must be at least 2.
!
!    Input, real ( kind = rk ) A_R83T(M,3), the R83T matrix.
!
!    Output, real ( kind = rk ) A_R8GE(M,N), the R8GE matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a_r83t(m,3)
  real ( kind = rk ) a_r8ge(m,n)
  integer i
  integer j
  integer k

  a_r8ge(1:m,1:n) = 0.0D+00

  do i = 1, m
    do k = 1, 3
      j = i - 2 + k
      if ( 1 <= j .and. j <= n ) then
        a_r8ge(i,j) = a_r83t(i,k)
      end if
    end do
  end do

  return
end
subroutine r83t_zeros ( m, n, a )

!*****************************************************************************80
!
!! R83T_ZEROS zeros an R83T matrix.
!
!  Discussion:
!
!    The R83T storage format is used for an MxN tridiagonal matrix.
!    The superdiagonal is stored in entries (1:M-1,3), the diagonal in
!    entries (1:M,2), and the subdiagonal in (2:M,1).  Thus, the
!    the rows of the original matrix slide horizontally to form an
!    Mx3 stack of data.
!
!    An R83T matrix of order 3x5 would be stored:
!
!       *  A11 A12
!      A21 A22 A23
!      A32 A33 A34
!
!    An R83T matrix of order 5x5 would be stored:
!
!       *  A11 A12
!      A21 A22 A23
!      A32 A33 A34
!      A43 A44 A45
!      A54 A55  *
!
!    An R83T matrix of order 5x3 would be stored:
!
!       *  A11 A12
!      A21 A22 A23
!      A32 A33  *
!      A43  *   *
!       *   *   *
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 August 2015
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
!    Output, real ( kind = rk ) A(M,3), the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,3)

  a(1:m,1:3) = 0.0D+00

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
!! R8VEC_INDICATOR1 sets an R8VEC to the indicator vector (1,2,3,...).
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 September 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements of A.
!
!    Output, real ( kind = rk ) A(N), the array.
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
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 August 2000
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
    write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
  end do

  return
end

