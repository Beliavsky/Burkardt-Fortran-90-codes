subroutine i4_fake_use ( n )

!*****************************************************************************80
!
!! i4_fake_use() pretends to use a variable.
!
!  Discussion:
!
!    Some compilers will issue a warning if a variable is unused.
!    Sometimes there's a good reason to include a variable in a program,
!    but not to use it.  Calling this function with that variable as
!    the argument will shut the compiler up.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 April 2020
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the variable to be "used".
!
  implicit none

  integer n

  if ( n /= n ) then
    write ( *, '(a)' ) '  i4_fake_use: variable is NAN.'
  end if

  return
end
subroutine orth_random ( n, a )

!*****************************************************************************80
!
!! ORTH_RANDOM returns the ORTH_RANDOM matrix.
!
!  Discussion:
!
!    The matrix is a random orthogonal matrix.
!
!  Properties:
!
!    The inverse of A is equal to A'.
!    A is orthogonal: A * A' = A' * A = I.
!    Because A is orthogonal, it is normal: A' * A = A * A'.
!    Columns and rows of A have unit Euclidean norm.
!    Distinct pairs of columns of A are orthogonal.
!    Distinct pairs of rows of A are orthogonal.
!    The L2 vector norm of A*x = the L2 vector norm of x for any vector x.
!    The L2 matrix norm of A*B = the L2 matrix norm of B for any matrix B.
!    det ( A ) = +1 or -1.
!    A is unimodular.
!    All the eigenvalues of A have modulus 1.
!    All singular values of A are 1.
!    All entries of A are between -1 and 1.
!
!  Discussion:
!
!    Thanks to Eugene Petrov, B I Stepanov Institute of Physics,
!    National Academy of Sciences of Belarus, for convincingly
!    pointing out the severe deficiencies of an earlier version of
!    this routine.
!
!    Essentially, the computation involves saving the Q factor of the
!    QR factorization of a matrix whose entries are normally distributed.
!    However, it is only necessary to generate this matrix a column at
!    a time, since it can be shown that when it comes time to annihilate
!    the subdiagonal elements of column K, these (transformed) elements of
!    column K are still normally distributed random values.  Hence, there
!    is no need to generate them at the beginning of the process and
!    transform them K-1 times.
!
!    For computational efficiency, the individual Householder transformations
!    could be saved, as recommended in the reference, instead of being
!    accumulated into an explicit matrix format.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Pete Stewart,
!    Efficient Generation of Random Orthogonal Matrices With an Application
!    to Condition Estimators,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 3, June 1980, pages 403-409.
!
!  Input:
!
!    integer N, the order of the matrix.
!
!  Output:
!
!    real ( kind = rk ) A(N,N), the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,n)
  integer i
  integer j
  real ( kind = rk ) r8_normal_01
  real ( kind = rk ) v(n)
  real ( kind = rk ) x(n)
!
!  Start with A = the identity matrix.
!
  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        a(i,j) = 1.0D+00
      else
        a(i,j) = 0.0D+00
      end if
    end do
  end do
!
!  Now behave as though we were computing the QR factorization of
!  some other random matrix.  Generate the N elements of the first column,
!  compute the Householder matrix H1 that annihilates the subdiagonal elements,
!  and set A := A * H1' = A * H.
!
!  On the second step, generate the lower N-1 elements of the second column,
!  compute the Householder matrix H2 that annihilates them,
!  and set A := A * H2' = A * H2 = H1 * H2.
!
!  On the N-1 step, generate the lower 2 elements of column N-1,
!  compute the Householder matrix HN-1 that annihilates them, and
!  and set A := A * H(N-1)' = A * H(N-1) = H1 * H2 * ... * H(N-1).
!  This is our random orthogonal matrix.
!
  do j = 1, n - 1
!
!  Set the vector that represents the J-th column to be annihilated.
!
    x(1:j-1) = 0.0D+00

    do i = j, n
      x(i) = r8_normal_01 ( )
    end do
!
!  Compute the vector V that defines a Householder transformation matrix
!  H(V) that annihilates the subdiagonal elements of X.
!
    call r8vec_house_column ( n, x, j, v )
!
!  Postmultiply the matrix A by H'(V) = H(V).
!
    call r8mat_house_axh ( n, a, v, a )

  end do

  return
end
function r8_normal_01 ( )

!*****************************************************************************80
!
!! r8_normal_01() returns a unit pseudonormal R8.
!
!  Discussion:
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!    Because this routine uses the Box Muller method, it requires pairs
!    of uniform random values to generate a pair of normal random values.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Output:
!
!    real ( kind = rk ) R8_NORMAL_01, a sample of the standard normal PDF.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r1
  real ( kind = rk ) r2
  real ( kind = rk ) r8_normal_01
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  integer, save :: used = 0
  real ( kind = rk ) x
  real ( kind = rk ), save :: y = 0.0D+00
!
!  On odd numbered calls, generate two uniforms, create two normals,
!  return the first normal.
!
  if ( mod ( used, 2 ) == 0 ) then

    call random_number ( harvest = r1 )
    call random_number ( harvest = r2 )

    x = sqrt ( -2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * r8_pi * r2 )
    y = sqrt ( -2.0D+00 * log ( r1 ) ) * sin ( 2.0D+00 * r8_pi * r2 )
!
!  On odd calls, return the second normal.
!
  else

    x = y

  end if

  used = used + 1

  r8_normal_01 = x

  return
end
subroutine r83_cg ( n, a, b, x )

!*****************************************************************************80
!
!! r83_cg() uses the conjugate gradient method on an R83 system.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
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
!    Here is how an R83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 June 2014
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
!    Input, real ( kind = rk ) A(3,N), the matrix.
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

  real ( kind = rk ) a(3,n)
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
  call r83_mv ( n, n, a, x, ap )

  r(1:n) = b(1:n) - ap(1:n)
  p(1:n) = b(1:n) - ap(1:n)
!
!  Do the N steps of the conjugate gradient method.
!
  do it = 1, n
!
!  Compute the matrix*vector product AP=A*P.
!
    call r83_mv ( n, n, a, p, ap )
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
subroutine r83_dif2 ( m, n, a )

!*****************************************************************************80
!
!! R83_DIF2 returns the DIF2 matrix in R83 format.
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
!    01 July 2000
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
!    Output, real ( kind = rk ) A(3,N), the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(3,n)
  integer mn

  a(1:3,1:n) = 0.0D+00

  mn = min ( m, n )

  a(2,1:mn)   = +2.0D+00
  a(1,2:mn)   = -1.0D+00

  if ( m <= n ) then
    a(3,1:mn-1) = -1.0D+00
  else if ( n < m ) then
    a(3,1:mn) = -1.0D+00
  end if
  
  return
end
subroutine r83_mv ( m, n, a, x, b )

!*****************************************************************************80
!
!! R83_MV multiplies an R83 matrix times an R8VEC.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!  Example:
!
!    Here is how an R83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 June 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = rk ) A(3,N), the R83 matrix.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = rk ) B(M), the product A * x.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(3,n)
  real ( kind = rk ) b(m)
  integer mn
  real ( kind = rk ) x(n)

  b(1:m) = 0.0D+00

  mn = min ( m, n )

  if ( n == 1 ) then
    b(1) = a(2,1) * x(1)
    if ( 1 < m ) then
      b(2) = a(3,1) * x(1)
    end if
    return
  end if

  b(1)      = a(2,1)      * x(1) &
            + a(1,2)      * x(2)

  b(2:mn-1) = a(3,1:mn-2) * x(1:mn-2) &
            + a(2,2:mn-1) * x(2:mn-1) &
            + a(1,3:mn)   * x(3:mn)

  b(mn)     = a(3,mn-1)   * x(mn-1) &
            + a(2,mn)     * x(mn)

  if ( n < m ) then
    b(n+1) = b(n+1) + a(3,n) * x(n)
  else if ( m < n ) then
    b(m) = b(m) + a(1,m+1) * x(m+1)
  end if

  return
end
subroutine r83_res ( m, n, a, x, b, r )

!*****************************************************************************80
!
!! R83_RES computes the residual R = B-A*X for R83 matrices.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 June 2014
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
!    Input, real ( kind = rk ) A(3,N), the matrix.
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

  real ( kind = rk ) a(3,n)
  real ( kind = rk ) b(m)
  real ( kind = rk ) r(m)
  real ( kind = rk ) x(n)

  call r83_mv ( m, n, a, x, r )

  r(1:m) = b(1:m) - r(1:m)

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
!      A2  A3   0   0   0
!      A1  A2  A3   0   0
!       0  A1  A2  A3   0 
!       0   0  A1  A2  A3
!       0   0   0  A1  A2
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
!    N must be positive.
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

  integer m
  integer n

  real ( kind = rk ) a(3)

  call i4_fake_use ( m )
  call i4_fake_use ( n )

  a(1) = -1.0D+00
  a(2) =  2.0D+00
  a(3) = -1.0D+00

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
!
!  Example:
!
!    Here is how an R83S matrix of order 5, stored as (A1,A2,A3), would
!    be interpreted:
!
!      A2  A3   0   0   0
!      A1  A2  A3   0   0
!       0  A1  A2  A3   0 
!       0   0  A1  A2  A3
!       0   0   0  A1  A2
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
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = rk ) A(3), the matrix.
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
  integer ihi
  integer ilo
  real ( kind = rk ) x(n)

  b(1:m) = 0.0D+00

  ilo = 2
  ihi = min ( m, n + 1 )
  b(ilo:ihi) = b(ilo:ihi) + a(1) * x(ilo-1:ihi-1)

  ilo = 1
  ihi = min ( m, n )
  b(ilo:ihi) = b(ilo:ihi) + a(2) * x(ilo:ihi)

  ilo = 1
  ihi = min ( m, n - 1 )
  b(ilo:ihi) = b(ilo:ihi) + a(3) * x(ilo+1:ihi+1)

  return
end
subroutine r83s_res ( m, n, a, x, b, r )

!*****************************************************************************80
!
!! R83S_RES computes the residual R = B-A*X for R83S matrices.
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
!  Parameters:
!
!    Input, integer M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
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
subroutine r83t_cg ( n, a, b, x )

!*****************************************************************************80
!
!! R83T_CG uses the conjugate gradient method on an R83T system.
!
!  Discussion:
!
!    The R83T storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1:N-1,3), the diagonal in
!    entries (1:N,2), and the subdiagonal in (2:N,1).  Thus, the
!    original matrix is "collapsed" horizontally into the array.
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
subroutine r83t_mv ( m, n, a, x, b )

!*****************************************************************************80
!
!! R83T_MV multiplies an R83T matrix times an R8VEC.
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
subroutine r83t_res ( m, n, a, x, b, r )

!*****************************************************************************80
!
!! R83T_RES computes the residual R = B-A*X for R83T matrices.
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
subroutine r8ge_cg ( n, a, b, x )

!*****************************************************************************80
!
!! R8GE_CG uses the conjugate gradient method on an R8GE system.
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
!    02 June 2014
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
!    Input, real ( kind = rk ) A(N,N), the matrix.
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

  real ( kind = rk ) a(n,n)
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
  call r8ge_mv ( n, n, a, x, ap )

  r(1:n) = b(1:n) - ap(1:n)
  p(1:n) = b(1:n) - ap(1:n)
!
!  Do the N steps of the conjugate gradient method.
!
  do it = 1, n
!
!  Compute the matrix*vector product AP=A*P.
!
    call r8ge_mv ( n, n, a, p, ap )
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
subroutine r8ge_dif2 ( m, n, a )

!*****************************************************************************80
!
!! R8GE_DIF2 returns the DIF2 matrix in R8GE format.
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
!    01 July 2000
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
!    Output, real ( kind = rk ) A(M,N), the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  integer i
  integer j

  do j = 1, n
    do i = 1, m

      if ( j == i - 1 ) then
        a(i,j) = -1.0D+00
      else if ( j == i ) then
        a(i,j) = 2.0D+00
      else if ( j == i + 1 ) then
        a(i,j) = -1.0D+00
      else
        a(i,j) = 0.0D+00
      end if

    end do
  end do

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
subroutine r8ge_res ( m, n, a, x, b, r )

!*****************************************************************************80
!
!! R8GE_RES computes the residual R = B-A*X for R8GE matrices.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 June 2014
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
!    Input, real ( kind = rk ) A(M,N), the matrix.
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

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) b(m)
  real ( kind = rk ) r(m)
  real ( kind = rk ) x(n)

  r(1:m) = b(1:m) - matmul ( a(1:m,1:n), x(1:n) )

  return
end
subroutine r8mat_house_axh ( n, a, v, ah )

!*****************************************************************************80
!
!! R8MAT_HOUSE_AXH computes A*H where H is a compact Householder matrix.
!
!  Discussion:
!
!    An R8MAT is a matrix of real ( kind = rk ) values.
!
!    The Householder matrix H(V) is defined by
!
!      H(V) = I - 2 * v * v' / ( v' * v )
!
!    This routine is not particularly efficient.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real ( kind = rk ) A(N,N), the matrix.
!
!    Input, real ( kind = rk ) V(N), a vector defining a Householder matrix.
!
!    Output, real ( kind = rk ) AH(N,N), the product A*H.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) ah(n,n)
  real ( kind = rk ) ah_temp(n,n)
  integer i
  integer j
  integer k
  real ( kind = rk ) v(n)
  real ( kind = rk ) v_normsq

  v_normsq = sum ( v(1:n) ** 2 )
!
!  Compute A*H' = A*H
!
  do i = 1, n
    do j = 1, n
      ah_temp(i,j) = a(i,j)
      do k = 1, n
        ah_temp(i,j) = ah_temp(i,j) - 2.0D+00 * a(i,k) * v(k) * v(j) / v_normsq
      end do
    end do
  end do
!
!  Copy the temporary result into AH.
!  Doing it this way means the user can identify the input arguments A and AH.
!
  ah(1:n,1:n) = ah_temp(1:n,1:n)

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows in A.
!
!    Input, integer N, the number of columns in A.
!
!    Input, real ( kind = rk ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  character ( len = * ) title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_PRINT_SOME prints some of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = rk ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ILO, JLO, the first row and column to print.
!
!    Input, integer IHI, JHI, the last row and column to print.
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

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)'
    return
  end if

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == real ( int ( a(i,j) ), kind = rk ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine r8pbu_cg ( n, mu, a, b, x )

!*****************************************************************************80
!
!! R8PBU_CG uses the conjugate gradient method on an R8PBU system.
!
!  Discussion:
!
!    The R8PBU storage format is for a symmetric positive definite band matrix.
!
!    To save storage, only the diagonal and upper triangle of A is stored,
!    in a compact diagonal format that preserves columns.
!
!    The diagonal is stored in row MU+1 of the array.
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
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
!    15 October 1998
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
!    Input, integer MU, the number of superdiagonals.
!    MU must be at least 0, and no more than N-1.
!
!    Input, real ( kind = rk ) A(MU+1,N), the R8PBU matrix.
!
!    Input, real ( kind = rk ) B(N), the right hand side vector.
!
!    Input/output, real ( kind = rk ) X(N).
!    On input, an estimate for the solution, which may be 0.
!    On output, the approximate solution vector.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer mu
  integer n

  real ( kind = rk ) a(mu+1,n)
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
  call r8pbu_mv ( n, n, mu, a, x, ap )

  r(1:n) = b(1:n) - ap(1:n)
  p(1:n) = b(1:n) - ap(1:n)
!
!  Do the N steps of the conjugate gradient method.
!
  do it = 1, n
!
!  Compute the matrix*vector product AP=A*P.
!
    call r8pbu_mv ( n, n, mu, a, p, ap )
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
subroutine r8pbu_dif2 ( m, n, mu, a )

!*****************************************************************************80
!
!! R8PBU_DIF2 returns the DIF2 matrix in R8PBU format.
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
!    01 July 2000
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
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer MU, the number of superdiagonals.
!    MU must be at least 0, and no more than N-1.
!
!    Output, real ( kind = rk ) A(MU+1,N), the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer mu
  integer n

  real ( kind = rk ) a(mu+1,n)
  integer m

  call i4_fake_use ( m )

  a(1:mu+1,1:n) = 0.0D+00

  a(mu,  2:n) = -1.0D+00
  a(mu+1,1:n) =  +2.0D+00
 
  return
end
subroutine r8pbu_mv ( m, n, mu, a, x, b )

!*****************************************************************************80
!
!! R8PBU_MV multiplies an R8PBU matrix by an R8VEC.
!
!  Discussion:
!
!    The R8PBU storage format is for a symmetric positive definite band matrix.
!
!    To save storage, only the diagonal and upper triangle of A is stored,
!    in a compact diagonal format that preserves columns.
!
!    The diagonal is stored in row MU+1 of the array.
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer MU, the number of superdiagonals in the matrix.
!    MU must be at least 0 and no more than N-1.
!
!    Input, real ( kind = rk ) A(MU+1,N), the R8PBU matrix.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = rk ) B(N), the result vector A * x.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer mu
  integer n

  real ( kind = rk ) a(mu+1,n)
  real ( kind = rk ) b(n)
  integer i
  integer ieqn
  integer j
  integer m
  real ( kind = rk ) x(n)

  call i4_fake_use ( m )
!
!  Multiply X by the diagonal of the matrix.
!
  b(1:n) = a(mu+1,1:n) * x(1:n)
!
!  Multiply X by the superdiagonals of the matrix.
!
  do i = mu, 1, -1
    do j = mu + 2 - i, n
      ieqn = i + j - mu - 1
      b(ieqn) = b(ieqn) + a(i,j) * x(j)
      b(j) = b(j) + a(i,j) * x(ieqn)
    end do
  end do

  return
end
subroutine r8pbu_res ( m, n, mu, a, x, b, r )

!*****************************************************************************80
!
!! R8PBU_RES computes the residual R = B-A*X for R8PBU matrices.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 June 2014
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
!    Input, integer MU, the number of superdiagonals in the matrix.
!    MU must be at least 0 and no more than N-1.
!
!    Input, real ( kind = rk ) A(MU+1,N), the matrix.
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
  integer mu
  integer n

  real ( kind = rk ) a(mu+1,n)
  real ( kind = rk ) b(m)
  real ( kind = rk ) r(m)
  real ( kind = rk ) x(n)

  call r8pbu_mv ( m, n, mu, a, x, r )

  r(1:m) = b(1:m) - r(1:m)

  return
end
subroutine r8sd_cg ( n, ndiag, offset, a, b, x )

!*****************************************************************************80
!
!! R8SD_CG uses the conjugate gradient method on an R8SD linear system.
!
!  Discussion:
!
!    The R8SD storage format is for symmetric matrices whose only nonzero
!    entries occur along a few diagonals, but for which these diagonals are 
!    not all close enough to the main diagonal for band storage to be efficient.
!
!    In that case, we assign the main diagonal the offset value 0, and 
!    each successive superdiagonal gets an offset value 1 higher, until
!    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
!
!    Assuming there are NDIAG nonzero diagonals (ignoring subdiagonals!),
!    we then create an array B that has N rows and NDIAG columns, and simply
!    "collapse" the matrix A to the left:
!
!    For the conjugate gradient method to be applicable, the matrix A must 
!    be a positive definite symmetric matrix.
!
!    The method is designed to reach the solution to the linear system
!      A * x = b
!    after N computational steps.  However, roundoff may introduce
!    unacceptably large errors for some problems.  In such a case,
!    calling the routine a second time, using the current solution estimate
!    as the new starting guess, should result in improved results.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 October 1998
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
!    Input, integer NDIAG, the number of diagonals that are stored.
!    NDIAG must be at least 1 and no more than N.
!
!    Input, integer OFFSET(NDIAG), the offsets for the diagonal
!    storage.
!
!    Input, real ( kind = rk ) A(N,NDIAG), the R8SD matrix.
!
!    Input, real ( kind = rk ) B(N), the right hand side vector.
!
!    Input/output, real ( kind = rk ) X(N).
!    On input, an estimate for the solution, which may be 0.
!    On output, the approximate solution vector.  Note that repeated
!    calls to this routine, using the value of X output on the previous
!    call, MAY improve the solution.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer ndiag

  real ( kind = rk ) a(n,ndiag)
  real ( kind = rk ) alpha
  real ( kind = rk ) ap(n)
  real ( kind = rk ) b(n)
  real ( kind = rk ) beta
  integer it
  integer offset(ndiag)
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
  call r8sd_mv ( n, n, ndiag, offset, a, x, ap )

  r(1:n) = b(1:n) - ap(1:n)
  p(1:n) = b(1:n) - ap(1:n)
!
!  Do the N steps of the conjugate gradient method.
!
  do it = 1, n
!
!  Compute the matrix*vector product AP = A*P.
!
    call r8sd_mv ( n, n, ndiag, offset, a, p, ap )
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
subroutine r8sd_dif2 ( m, n, ndiag, offset, a )

!*****************************************************************************80
!
!! R8SD_DIF2 returns the DIF2 matrix in R8SD format.
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
!    05 July 2000
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
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer NDIAG, the number of diagonals that are stored.
!    NDIAG must be at least 2.
!
!    Input, integer OFFSET(NDIAG), the offsets for the diagonal
!    storage.  It is simply presumed that OFFSET(1) = 0 and OFFSET(2) = 1.
!
!    Output, real ( kind = rk ) A(N,NDIAG), the R8SD matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer ndiag

  real ( kind = rk ) a(n,ndiag)
  integer m
  integer offset(ndiag)

  call i4_fake_use ( m )
  call i4_fake_use ( offset(1) )

  a(1:n,1:ndiag) = 0.0D+00

  a(1:n,  1) =  2.0D+00
  a(1:n-1,2) = -1.0D+00
 
  return
end
subroutine r8sd_mv ( m, n, ndiag, offset, a, x, b )

!*****************************************************************************80
!
!! R8SD_MV multiplies an R8SD matrix by an R8VEC.
!
!  Discussion:
!
!    The R8SD storage format is for symmetric matrices whose only nonzero 
!    entries occur along a few diagonals, but for which these diagonals are not 
!    all close enough to the main diagonal for band storage to be efficient.
!
!    In that case, we assign the main diagonal the offset value 0, and 
!    each successive superdiagonal gets an offset value 1 higher, until
!    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
!
!    Assuming there are NDIAG nonzero diagonals (ignoring subdiagonals!),
!    we then create an array B that has N rows and NDIAG columns, and simply
!    "collapse" the matrix A to the left:
!
!  Example:
!
!    The "offset" value is printed above each column.
!
!    Original matrix               New Matrix
!
!       0   1   2   3   4   5       0   1   3   5
!
!      11  12   0  14   0  16      11  12  14  16
!      21  22  23   0  25   0      22  23  25  --
!       0  32  33  34   0  36      33  34  36  --
!      41   0  43  44  45   0      44  45  --  --
!       0  52   0  54  55  56      55  56  --  --
!      61   0  63   0  65  66      66  --  --  --
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer NDIAG, the number of diagonals that are stored.
!    NDIAG must be at least 1 and no more than N.
!
!    Input, integer OFFSET(NDIAG), the offsets for the diagonal
!    storage.
!
!    Input, real ( kind = rk ) A(N,NDIAG), the R8SD matrix.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = rk ) B(N), the product A * x.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n
  integer ndiag

  real ( kind = rk ) a(n,ndiag)
  real ( kind = rk ) b(n)
  integer i
  integer j
  integer jdiag
  integer offset(ndiag)
  real ( kind = rk ) x(n)

  call i4_fake_use ( m )

  b(1:n) = 0.0D+00

  do i = 1, n
    do jdiag = 1, ndiag
      if ( 0 <= offset(jdiag) ) then
        j = i + offset(jdiag)
        if ( 1 <= j .and. j <= n ) then
          b(i) = b(i) + a(i,jdiag) * x(j)
          if ( offset(jdiag) /= 0 ) then
            b(j) = b(j) + a(i,jdiag) * x(i)
          end if
        end if
      end if
    end do
  end do

  return
end
subroutine r8sd_res ( m, n, ndiag, offset, a, x, b, r )

!*****************************************************************************80
!
!! R8SD_RES computes the residual R = B-A*X for R8SD matrices.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 June 2014
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
!    Input, integer NDIAG, the number of diagonals that are stored.
!    NDIAG must be at least 1 and no more than N.
!
!    Input, integer OFFSET(NDIAG), the offsets for the diagonal
!    storage.
!
!    Input, real ( kind = rk ) A(N,NDIAG), the R8SD matrix.
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
  integer ndiag

  real ( kind = rk ) a(n,ndiag)
  real ( kind = rk ) b(m)
  integer offset(ndiag)
  real ( kind = rk ) r(m)
  real ( kind = rk ) x(n)

  call r8sd_mv ( m, n, ndiag, offset, a, x, r )

  r(1:m) = b(1:m) - r(1:m)

  return
end
subroutine r8sp_cg ( n, nz_num, row, col, a, b, x )

!*****************************************************************************80
!
!! R8SP_CG uses the conjugate gradient method on an R8SP system.
!
!  Discussion:
!
!    The R8SP storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.
!
!    It is possible that a pair of indices (I,J) may occur more than
!    once.  Presumably, in this case, the intent is that the actual value
!    of A(I,J) is the sum of all such entries.  This is not a good thing
!    to do, but I seem to have come across this in MATLAB.
!
!    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP 
!    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
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
!    02 June 2014
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
!    Input, integer NZ_NUM, the number of nonzero elements in
!    the matrix.
!
!    Input, integer ROW(NZ_NUM), COL(NZ_NUM), the row and 
!    column indices of the nonzero elements.
!
!    Input, real ( kind = rk ) A(NZ_NUM), the nonzero elements of the matrix.
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
  integer nz_num

  real ( kind = rk ) a(nz_num)
  real ( kind = rk ) alpha
  real ( kind = rk ) ap(n)
  real ( kind = rk ) b(n)
  real ( kind = rk ) beta
  integer col(nz_num)
  integer it
  real ( kind = rk ) p(n)
  real ( kind = rk ) pap
  real ( kind = rk ) pr
  real ( kind = rk ) r(n)
  integer row(nz_num)
  real ( kind = rk ) rap
  real ( kind = rk ) x(n)
!
!  Initialize
!    AP = A * x,
!    R  = b - A * x,
!    P  = b - A * x.
!
  call r8sp_mv ( n, n, nz_num, row, col, a, x, ap )

  r(1:n) = b(1:n) - ap(1:n)
  p(1:n) = b(1:n) - ap(1:n)
!
!  Do the N steps of the conjugate gradient method.
!
  do it = 1, n
!
!  Compute the matrix*vector product AP=A*P.
!
    call r8sp_mv ( n, n, nz_num, row, col, a, p, ap )
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
subroutine r8sp_dif2 ( m, n, nz_num, row, col, a )

!*****************************************************************************80
!
!! R8SP_DIF2 returns the DIF2 matrix in R8SP format.
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
!    01 July 2000
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
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer NZ_NUM, the number of nonzero elements in
!    the matrix.
!
!    Output, integer ROW(NZ_NUM), COL(NZ_NUM), the row and 
!    column indices of the nonzero elements.
!
!    Output, real ( kind = rk ) A(NZ_NUM), the nonzero elements of the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nz_num

  real ( kind = rk ) a(nz_num)
  integer col(nz_num)
  integer i
  integer k
  integer m
  integer mn
  integer row(nz_num)

  mn = min ( m, n );

  k = 0
  do i = 1, mn

    if ( 0 < i - 1 ) then
      k = k + 1
      row(k) = i
      col(k) = i - 1
      a(k) = -1.0D+00
    end if

    k = k + 1
    row(k) = i
    col(k) = i
    a(k) = 2.0D+00

    if ( i < n ) then
      k = k + 1
      row(k) = i
      col(k) = i + 1
      a(k) = -1.0D+00
    end if

  end do

  return
end
subroutine r8sp_mv ( m, n, nz_num, row, col, a, x, b )

!*****************************************************************************80
!
!! R8SP_MV multiplies an R8SP matrix by an R8VEC.
!
!  Discussion:
!
!    The R8SP storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.
!
!    It is possible that a pair of indices (I,J) may occur more than
!    once.  Presumably, in this case, the intent is that the actual value
!    of A(I,J) is the sum of all such entries.  This is not a good thing
!    to do, but I seem to have come across this in MATLAB.
!
!    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP 
!    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns of 
!    the matrix.
!
!    Input, integer NZ_NUM, the number of nonzero elements in
!    the matrix.
!
!    Input, integer ROW(NZ_NUM), COL(NZ_NUM), the row and 
!    column indices of the nonzero elements.
!
!    Input, real ( kind = rk ) A(NZ_NUM), the nonzero elements of the matrix.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = rk ) B(M), the product vector A*X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n
  integer nz_num

  real ( kind = rk ) a(nz_num)
  real ( kind = rk ) b(m)
  integer col(nz_num)
  integer i
  integer j
  integer k
  integer row(nz_num)
  real ( kind = rk ) x(n)

  b(1:m) = 0.0D+00

  do k = 1, nz_num

    i = row(k)
    j = col(k)
    b(i) = b(i) + a(k) * x(j)

  end do

  return
end
subroutine r8sp_res ( m, n, nz_num, row, col, a, x, b, r )

!*****************************************************************************80
!
!! R8SP_RES computes the residual R = B-A*X for R8SP matrices.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 June 2014
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
!    Input, integer NZ_NUM, the number of nonzero elements in
!    the matrix.
!
!    Input, integer ROW(NZ_NUM), COL(NZ_NUM), the row and 
!    column indices of the nonzero elements.
!
!    Input, real ( kind = rk ) A(NZ_NUM), the nonzero elements of the matrix.
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
  integer nz_num

  real ( kind = rk ) a(nz_num)
  real ( kind = rk ) b(m)
  integer col(nz_num)
  real ( kind = rk ) r(m)
  integer row(nz_num)
  real ( kind = rk ) x(n)

  call r8sp_mv ( m, n, nz_num, row, col, a, x, r )

  r(1:m) = b(1:m) - r(1:m)

  return
end
subroutine r8vec_house_column ( n, a, k, v )

!*****************************************************************************80
!
!! R8VEC_HOUSE_COLUMN defines a Householder premultiplier that "packs" a column.
!
!  Discussion:
!
!    An R8VEC is a vector of real ( kind = rk ) values.
!
!    The routine returns a vector V that defines a Householder
!    premultiplier matrix H(V) that zeros out the subdiagonal entries of
!    column K of the matrix A.
!
!       H(V) = I - 2 * v * v'
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix A.
!
!    Input, real ( kind = rk ) A(N), column K of the matrix A.
!
!    Input, integer K, the column of the matrix to be modified.
!
!    Output, real ( kind = rk ) V(N), a vector of unit L2 norm which defines an
!    orthogonal Householder premultiplier matrix H with the property
!    that the K-th column of H*A is zero below the diagonal.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  integer k
  real ( kind = rk ) s
  real ( kind = rk ) v(n)

  v(1:n) = 0.0D+00

  if ( k < 1 .or. n <= k ) then
    return
  end if

  s = sqrt ( dot_product ( a(k:n), a(k:n) ) )

  if ( s == 0.0D+00 ) then
    return
  end if

  v(k) = a(k) + sign ( s, a(k) )
  v(k+1:n) = a(k+1:n)

  v(k:n) = v(k:n) / sqrt ( dot_product ( v(k:n), v(k:n) ) )

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
subroutine spd_random ( n, a )

!*****************************************************************************80
!
!! spd_random() returns a random symmetric positive definite matrix.
!
!  Discussion:
!
!    The matrix is a "random" symmetric positive definite matrix.
!
!    The matrix returned will have eigenvalues in the range [0,1].
!
!  Properties:
!
!    A is symmetric: A' = A.
!
!    A is positive definite: 0 < x'*A*x for nonzero x.
!
!    The eigenvalues of A will be real.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 May 2002
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the order of the matrix.
!
!  Output:
!
!    real ( kind = rk ) A(N,N), the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,n)
  integer i
  integer j
  integer k
  real ( kind = rk ) lambda(n)
  real ( kind = rk ) q(n,n)
!
!  Get a random set of eigenvalues.
!
  call random_number ( harvest = lambda(1:n) )
!
!  Get a random orthogonal matrix Q.
!
  call orth_random ( n, q )
!
!  Set A = Q * Lambda * Q'.
!
  do i = 1, n
    do j = 1, n
      a(i,j) = 0.0D+00
      do k = 1, n
        a(i,j) = a(i,j) + q(i,k) * lambda(k) * q(j,k)
      end do
    end do
  end do

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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
 
