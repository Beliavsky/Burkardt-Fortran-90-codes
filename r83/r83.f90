subroutine i4_fake_use ( n )

!*****************************************************************************80
!
!! i4_fake_use pretends to use a variable.
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

  integer, parameter :: rk = kind ( 1.0D+00 )

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
subroutine r83_cg ( n, a, b, x )

!*****************************************************************************80
!
!! R83_CG uses the conjugate gradient method on an R83 system.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:min(M+1,N)).
!    The diagonal in entries (2,1:min(M,N)).
!    The subdiagonal in (3,min(M-1,N)).
!    R8GE A(I,J) = R83 A(I-J+2,J).
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
subroutine r83_cr_fa ( n, a, a_cr )

!*****************************************************************************80
!
!! R83_CR_FA decomposes an R83 matrix using cyclic reduction.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:min(M+1,N)).
!    The diagonal in entries (2,1:min(M,N)).
!    The subdiagonal in (3,min(M-1,N)). 
!    R8GE A(I,J) = R83 A(I-J+2,J).
!
!    Once R83_CR_FA has decomposed a matrix A, then R83_CR_SL may be used 
!    to solve linear systems A * x = b.
!
!    R83_CR_FA does not employ pivoting.  Hence, the results can be more
!    sensitive to ill-conditioning than standard Gauss elimination.  In
!    particular, R83_CR_FA will fail if any diagonal element of the matrix
!    is zero.  Other matrices may also cause R83_CR_FA to fail.
!
!    R83_CR_FA can be guaranteed to work properly if the matrix is strictly
!    diagonally dominant, that is, if the absolute value of the diagonal
!    element is strictly greater than the sum of the absolute values of
!    the offdiagonal elements, for each equation.
!
!    The algorithm may be illustrated by the following figures:
!
!    The initial matrix is given by:
!
!          D1 U1
!          L1 D2 U2
!             L2 D3 U3
!                L3 D4 U4
!                   L4 D5 U5
!                      L5 D6
!
!    Rows and columns are permuted in an odd/even way to yield:
!
!          D1       U1
!             D3    L2 U3
!                D5    L4 U5
!          L1 U2    D2
!             L3 U4    D4
!                L5       D6
!
!    A block LU decomposition is performed to yield:
!
!          D1      |U1
!             D3   |L2 U3
!                D5|   L4 U5
!          --------+--------
!                  |D2'F3
!                  |F1 D4'F4
!                  |   F2 D6'
!
!    For large systems, this reduction is repeated on the lower right hand
!    tridiagonal subsystem until a completely upper triangular system
!    is obtained.  The system has now been factored into the product of a
!    lower triangular system and an upper triangular one, and the information
!    defining this factorization may be used by R83_CR_SL to solve linear
!    systems.
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
!    23 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Roger Hockney,
!    A fast direct solution of Poisson's equation using Fourier Analysis,
!    Journal of the ACM,
!    Volume 12, Number 1, pages 95-113, January 1965.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real ( kind = rk ) A(3,N), the R83 matrix.
!
!    Output, real ( kind = rk ) A_CR(3,0:2*N), factorization information 
!    needed by R83_CR_SL.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(3,n)
  real ( kind = rk ) a_cr(3,0:2*n)
  integer iful
  integer ifulp
  integer ihaf
  integer il
  integer ilp
  integer inc
  integer incr
  integer ipnt
  integer ipntp

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R83_CR_FA - Fatal error!'
    write ( *, '(a,i8)' ) '  Nonpositive N = ', n
    stop 1
  end if

  a_cr(1:3,0:2*n) = 0.0D+00
!
!  Copy the R83 matrix.
!
  if ( n == 1 ) then
    a_cr(2,1) = 1.0D+00 / a(2,1)
    return
  end if

  a_cr(1,1:n-1) = a(1,2:n)
  a_cr(2,1:n) = a(2,1:n)
  a_cr(3,1:n-1) = a(3,1:n-1)

  il = n
  ipntp = 0

  do while ( 1 < il )

    ipnt = ipntp
    ipntp = ipntp + il
    if ( mod ( il, 2 ) == 1 ) then
      inc = il + 1
    else
      inc = il
    end if

    incr = inc / 2
    il = il / 2
    ihaf = ipntp + incr + 1
    ifulp = ipnt + inc + 2

    do ilp = incr, 1, -1
      ifulp = ifulp - 2
      iful = ifulp - 1
      ihaf = ihaf - 1
      a_cr(2,iful) = 1.0D+00 / a_cr(2,iful)
      a_cr(3,iful)  = a_cr(3,iful)  * a_cr(2,iful)
      a_cr(1,ifulp) = a_cr(1,ifulp) * a_cr(2,ifulp+1)
      a_cr(2,ihaf)  = a_cr(2,ifulp) - a_cr(1,iful)  * a_cr(3,iful) &
                                  - a_cr(1,ifulp) * a_cr(3,ifulp)
      a_cr(3,ihaf) = - a_cr(3,ifulp) * a_cr(3,ifulp+1)
      a_cr(1,ihaf) = - a_cr(1,ifulp) * a_cr(1,ifulp+1)
    end do

  end do

  a_cr(2,ipntp+1) = 1.0D+00 / a_cr(2,ipntp+1)

  return
end
subroutine r83_cr_sl ( n, a_cr, b, x )

!*****************************************************************************80
!
!! R83_CR_SL solves a linear systems factored by R83_CR_FA.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:min(M+1,N)).
!    The diagonal in entries (2,1:min(M,N)).
!    The subdiagonal in (3,min(M-1,N)).
!    R8GE A(I,J) = R83 A(I-J+2,J).
!
!    The matrix A must be tridiagonal.  R83_CR_FA is called to compute the
!    LU factors of A.  It does so using a form of cyclic reduction.  If
!    the factors computed by R83_CR_FA are passed to R83_CR_SL, then 
!    a linear system involving the matrix A may be solved.
!
!    Note that R83_CR_FA does not perform pivoting, and so the solution 
!    produced by R83_CR_SL may be less accurate than a solution produced 
!    by a standard Gauss algorithm.  However, such problems can be 
!    guaranteed not to occur if the matrix A is strictly diagonally 
!    dominant, that is, if the absolute value of the diagonal coefficient 
!    is greater than the sum of the absolute values of the two off diagonal 
!    coefficients, for each row of the matrix.
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
!    06 May 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Roger Hockney,
!    A fast direct solution of Poisson's equation using Fourier Analysis,
!    Journal of the ACM,
!    Volume 12, Number 1, pages 95-113, January 1965.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real ( kind = rk ) A_CR(3,0:2*N), factorization information 
!    computed by R83_CR_FA.
!
!    Input, real ( kind = rk ) B(N), the right hand side.
!
!    Output, real ( kind = rk ) X(N), the solution of the linear systems.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a_cr(3,0:2*n)
  real ( kind = rk ) b(n)
  integer iful
  integer ifulm
  integer ihaf
  integer il
  integer ipnt
  integer ipntp
  integer ndiv
  real ( kind = rk ) rhs(0:2*n)
  real ( kind = rk ) x(n)

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R83_CR_SL - Fatal error!'
    write ( *, '(a,i8)' ) '  Nonpositive N = ', n
    stop 1
  end if

  if ( n == 1 ) then
    x(1) = a_cr(2,1) * b(1)
    return
  end if
!
!  Set up RHS.
!
  rhs(0) = 0.0D+00
  rhs(1:n) = b(1:n)
  rhs(n+1:2*n) = 0.0D+00

  il = n
  ndiv = 1
  ipntp = 0

  do while ( 1 < il )

    ipnt = ipntp
    ipntp = ipntp + il
    il = il / 2
    ndiv = ndiv * 2
    ihaf = ipntp

    do iful = ipnt + 2, ipntp, 2
      ihaf = ihaf + 1
      rhs(ihaf) = rhs(iful) &
        - a_cr(3,iful-1) * rhs(iful-1) &
        - a_cr(1,iful)   * rhs(iful+1)
    end do

  end do

  rhs(ihaf) = a_cr(2,ihaf) * rhs(ihaf)

  ipnt = ipntp

  do while ( 0 < ipnt )

    ipntp = ipnt
    ndiv = ndiv / 2
    il = n / ndiv
    ipnt = ipnt - il
    ihaf = ipntp

    do ifulm = ipnt + 1, ipntp, 2
      iful = ifulm + 1
      ihaf = ihaf + 1
      rhs(iful) = rhs(ihaf)
      rhs(ifulm) = a_cr(2,ifulm) &
        * (                     rhs(ifulm) &
            - a_cr(3,ifulm-1) * rhs(ifulm-1) &
            - a_cr(1,ifulm)   * rhs(iful) )
    end do

  end do

  x(1:n) = rhs(1:n)

  return
end
subroutine r83_cr_sls ( n, a_cr, nb, b, x )

!*****************************************************************************80
!
!! R83_CR_SLS solves several linear systems factored by R83_CR_FA.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:min(M+1,N)).
!    The diagonal in entries (2,1:min(M,N)).
!    The subdiagonal in (3,min(M-1,N)).
!    R8GE A(I,J) = R83 A(I-J+2,J).
!
!    The matrix A must be tridiagonal.  R83_CR_FA is called to compute the
!    LU factors of A.  It does so using a form of cyclic reduction.  If
!    the factors computed by R83_CR_FA are passed to R83_CR_SLS, then one or 
!    many linear systems involving the matrix A may be solved.
!
!    Note that R83_CR_FA does not perform pivoting, and so the solution 
!    produced by R83_CR_SLS may be less accurate than a solution produced 
!    by a standard Gauss algorithm.  However, such problems can be 
!    guaranteed not to occur if the matrix A is strictly diagonally 
!    dominant, that is, if the absolute value of the diagonal coefficient 
!    is greater than the sum of the absolute values of the two off diagonal 
!    coefficients, for each row of the matrix.
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
!    29 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Roger Hockney,
!    A fast direct solution of Poisson's equation using Fourier Analysis,
!    Journal of the ACM,
!    Volume 12, Number 1, pages 95-113, January 1965.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real ( kind = rk ) A_CR(3,0:2*N), factorization information 
!    computed by R83_CR_FA.
!
!    Input, integer NB, the number of right hand sides.
!
!    Input, real ( kind = rk ) B(N,NB), the right hand sides.
!
!    Output, real ( kind = rk ) X(N,NB), the solutions of the linear systems.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nb

  real ( kind = rk ) a_cr(3,0:2*n)
  real ( kind = rk ) b(n,nb)
  integer iful
  integer ifulm
  integer ihaf
  integer il
  integer ipnt
  integer ipntp
  integer ndiv
  real ( kind = rk ) rhs(0:2*n,nb)
  real ( kind = rk ) x(n,nb)

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R83_CR_SLS - Fatal error!'
    write ( *, '(a,i8)' ) '  Nonpositive N = ', n
    stop 1
  end if

  if ( n == 1 ) then
    x(1,1:nb) = a_cr(2,1) * b(1,1:nb)
    return
  end if
!
!  Set up RHS.
!
  rhs(0,1:nb) = 0.0D+00
  rhs(1:n,1:nb) = b(1:n,1:nb)
  rhs(n+1:2*n,1:nb) = 0.0D+00

  il = n
  ndiv = 1
  ipntp = 0

  do while ( 1 < il )

    ipnt = ipntp
    ipntp = ipntp + il
    il = il / 2
    ndiv = ndiv * 2
    ihaf = ipntp

    do iful = ipnt + 2, ipntp, 2
      ihaf = ihaf + 1
      rhs(ihaf,1:nb) = rhs(iful,1:nb) &
        - a_cr(3,iful-1) * rhs(iful-1,1:nb) &
        - a_cr(1,iful)   * rhs(iful+1,1:nb)
    end do

  end do

  rhs(ihaf,1:nb) = rhs(ihaf,1:nb) * a_cr(2,ihaf)
  ipnt = ipntp

  do while ( 0 < ipnt )

    ipntp = ipnt
    ndiv = ndiv / 2
    il = n / ndiv
    ipnt = ipnt - il
    ihaf = ipntp

    do ifulm = ipnt + 1, ipntp, 2
      iful = ifulm + 1
      ihaf = ihaf + 1
      rhs(iful,1:nb) = rhs(ihaf,1:nb)
      rhs(ifulm,1:nb) = a_cr(2,ifulm) &
        * (                     rhs(ifulm,1:nb) &
            - a_cr(3,ifulm-1) * rhs(ifulm-1,1:nb) &
            - a_cr(1,ifulm)   * rhs(iful,1:nb) )
    end do

  end do

  x(1:n,1:nb) = rhs(1:n,1:nb)

  return
end
subroutine r83_dif2 ( m, n, a )

!*****************************************************************************80
!
!! R83_DIF2 returns the DIF2 matrix in R83 format.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:min(M+1,N)).
!    The diagonal in entries (2,1:min(M,N)).
!    The subdiagonal in (3,min(M-1,N)).
!    R8GE A(I,J) = R83 A(I-J+2,J).
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
!    05 September 2015
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
  integer i
  integer j

  a(1:3,1:n) = 0.0D+00

  do j = 1, n
    do i = max ( 1, j - 1 ), min ( m, j + 1 )
      if ( i - j + 2 == 1 ) then
        a(i-j+2,j) = -1.0D+00
      else if ( i - j + 2 == 2 ) then
        a(i-j+2,j) = +2.0D+00
      else if ( i - j + 2 == 3 ) then
        a(i-j+2,j) = -1.0D+00
      end if
    end do
  end do
  
  return
end
subroutine r83_gs_sl ( n, a, b, x, it_max )

!*****************************************************************************80
!
!! R83_GS_SL solves an R83 system using Gauss-Seidel iteration.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:min(M+1,N)).
!    The diagonal in entries (2,1:min(M,N)).
!    The subdiagonal in (3,min(M-1,N)).
!    R8GE A(I,J) = R83 A(I-J+2,J).
!
!    This routine simply applies a given number of steps of the
!    iteration to an input approximate solution.  On first call, you can
!    simply pass in the zero vector as an approximate solution.  If
!    the returned value is not acceptable, you may call again, using
!    it as the starting point for additional iterations.
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
!    10 September 2015
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
!    Input, real ( kind = rk ) A(3,N), the R83 matrix.
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

  real ( kind = rk ) a(3,n)
  real ( kind = rk ) b(n)
  integer i
  integer it_max
  integer it_num
  real ( kind = rk ) x(n)
  real ( kind = rk ) x_old(n)
!
!  No diagonal matrix entry can be zero.
!
  do i = 1, n
    if ( a(2,i) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R83_GS_SL - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero diagonal entry, index = ', i
      stop 1
    end if
  end do

  do it_num = 1, it_max

    x_old(1:n) = x(1:n)

    x(1) = ( b(1) - a(1,2) * x(2) ) / a(2,1)
    do i = 2, n - 1
      x(i) = ( b(i) - a(3,i-1) * x(i-1) - a(1,i+1) * x(i+1) ) / a(2,i)
    end do
    x(n) = ( b(n) - a(3,n-1) * x(n-1) ) / a(2,n)

  end do

  return
end
subroutine r83_indicator ( m, n, a )

!*****************************************************************************80
!
!! R83_INDICATOR sets up an R83 indicator matrix.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:min(M+1,N)).
!    The diagonal in entries (2,1:min(M,N)).
!    The subdiagonal in (3,min(M-1,N)).
!    R8GE A(I,J) = R83 A(I-J+2,J).
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
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
!    Output, real ( kind = rk ) A(3,N), the R83 indicator matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(3,n)
  integer fac
  integer i
  integer i4_log_10
  integer j
  integer m

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  a(1:3,1:n) = 0.0D+00

  do j = 1, n
    do i = max ( 1, j - 1 ), min ( m, j + 1 )
      a(i-j+2,j) = real ( fac * i + j, kind = rk )
    end do
  end do

  return
end
subroutine r83_jac_sl ( n, a, b, x, it_max )

!*****************************************************************************80
!
!! R83_JAC_SL solves an R83 system using Jacobi iteration.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:min(M+1,N)).
!    The diagonal in entries (2,1:min(M,N)).
!    The subdiagonal in (3,min(M-1,N)).
!    R8GE A(I,J) = R83 A(I-J+2,J).
!
!    This routine simply applies a given number of steps of the
!    iteration to an input approximate solution.  On first call, you can
!    simply pass in the zero vector as an approximate solution.  If
!    the returned value is not acceptable, you may call again, using
!    it as the starting point for additional iterations.
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
!    05 September 2015
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
!    Input, real ( kind = rk ) A(3,N), the R83 matrix.
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

  real ( kind = rk ) a(3,n)
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
    if ( a(2,i) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R83_JAC_SL - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero diagonal entry, index = ', i
      stop 1
    end if
  end do

  do it_num = 1, it_max

    x_new(1) = b(1) - a(1,2) * x(2)
    do i = 2, n - 1
      x_new(i) = b(i) - a(3,i-1) * x(i-1) - a(1,i+1) * x(i+1)
    end do
    x_new(n) = b(n) - a(3,n-1) * x(n-1)
!
!  Divide by diagonal terms.
!
    x_new(1:n) = x_new(1:n) / a(2,1:n)
!
!  Update.
!
    x(1:n) = x_new(1:n)

  end do

  return
end
subroutine r83_mtv ( m, n, a, x, b )

!*****************************************************************************80
!
!! R83_MTV computes A'*x=b, where A is an R83 matrix.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:min(M+1,N)).
!    The diagonal in entries (2,1:min(M,N)).
!    The subdiagonal in (3,min(M-1,N)).
!    R8GE A(I,J) = R83 A(I-J+2,J).
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
!    05 September 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the order of the linear system.
!
!    Input, real ( kind = rk ) A(3,N), the MxN R83 matrix.
!
!    Input, real ( kind = rk ) X(M), the vector to be multiplied by A'.
!
!    Output, real ( kind = rk ) B(N), the product A' * x.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(3,n)
  real ( kind = rk ) b(n)
  integer i
  integer j
  real ( kind = rk ) x(m)

  b(1:n) = 0.0D+00
!
!  Find each nonzero A(I,J), multiply by X(I), add to B(J).
!
  do j = 1, n
    do i = max ( 1, j - 1 ), min ( m, j + 1 )
      b(j) = b(j) + x(i) * a(i-j+2,j)
    end do
  end do

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
!    The superdiagonal is stored in entries (1,2:min(M+1,N)).
!    The diagonal in entries (2,1:min(M,N)).
!    The subdiagonal in (3,min(M-1,N)).
!    R8GE A(I,J) = R83 A(I-J+2,J).
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
!    05 September 2015
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
  integer i
  integer j
  real ( kind = rk ) x(n)

  b(1:m) = 0.0D+00

  do j = 1, n
    do i = max ( 1, j - 1 ), min ( m, j + 1 )
      b(i) = b(i) + a(i-j+2,j) * x(j)
    end do
  end do

  return
end
subroutine r83_print ( m, n, a, title )

!*****************************************************************************80
!
!! R83_PRINT prints an R83 matrix.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:min(M+1,N)).
!    The diagonal in entries (2,1:min(M,N)).
!    The subdiagonal in (3,min(M-1,N)).
!    R8GE A(I,J) = R83 A(I-J+2,J).
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
!    31 August 2015
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
!    Input, real ( kind = rk ) A(3,N), the R83 matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(3,n)
  integer m
  character ( len = * ) title

  call r83_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r83_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R83_PRINT_SOME prints some of an R83 matrix.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:min(M+1,N)).
!    The diagonal in entries (2,1:min(M,N)).
!    The subdiagonal in (3,min(M-1,N)).
!    R8GE A(I,J) = R83 A(I-J+2,J).
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
!    05 September 2015
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
!    Input, real ( kind = rk ) A(3,N), the R83 matrix.
!
!    Input, integer ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column, to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: incx = 5
  integer n

  real ( kind = rk ) a(3,n)
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
          write ( ctemp(j2), '(g14.6)' ) a(i-j+2,j)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r83_random ( m, n, a )

!*****************************************************************************80
!
!! R83_RANDOM randomizes an R83 matrix.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:min(M+1,N)).
!    The diagonal in entries (2,1:min(M,N)).
!    The subdiagonal in (3,min(M-1,N)).
!    R8GE A(I,J) = R83 A(I-J+2,J).
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
!    05 September 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the linear system.
!
!    Output, real ( kind = rk ) A(3,N), the R83 matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(3,n)
  integer i
  integer j
  integer m

  a(1:3,1:n) = 0.0D+00

  do j = 1, n
    do i = max ( 1, j - 1 ), min ( m, j + 1 )
      call random_number ( harvest = a(i-j+2,j) )
    end do
  end do

  return
end
subroutine r83_res ( m, n, a, x, b, r )

!*****************************************************************************80
!
!! R83_RES computes the residual R = B-A*X for R83 matrices.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:min(M+1,N)).
!    The diagonal in entries (2,1:min(M,N)).
!    The subdiagonal in (3,min(M-1,N)).
!    R8GE A(I,J) = R83 A(I-J+2,J).
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
subroutine r83_to_r8ge ( m, n, a, b )

!*****************************************************************************80
!
!! R83_TO_R8GE copies an R83 matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:min(M+1,N)).
!    The diagonal in entries (2,1:min(M,N)).
!    The subdiagonal in (3,min(M-1,N)).
!    R8GE A(I,J) = R83 A(I-J+2,J).
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
!    05 September 2015
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
!    Input, real ( kind = rk ) A(3,N), the R83 matrix.
!
!    Output, real ( kind = rk ) B(M,N), the R8GE matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(3,n)
  real ( kind = rk ) b(m,n)
  integer i
  integer j

  b(1:m,1:n) = 0.0D+00

  do j = 1, n
    do i = max ( 1, j - 1 ), min ( m, j + 1 )
      b(i,j) = a(i-j+2,j)
    end do
  end do

  return
end
subroutine r83_zeros ( m, n, a )

!*****************************************************************************80
!
!! R83_ZEROS zeroes an R83 matrix.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:min(M+1,N)).
!    The diagonal in entries (2,1:min(M,N)).
!    The subdiagonal in (3,min(M-1,N)).
!    R8GE A(I,J) = R83 A(I-J+2,J).
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
!    26 January 2013
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, N, the order of the linear system.
!
!  Output:
!
!    real ( kind = rk ) A(3,N), the R83 matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(3,n)
  integer m

  call i4_fake_use ( m )

  a(1:3,1:n) = 0.0D+00

  return
end
subroutine r83np_fs ( n, a, b, x )

!*****************************************************************************80
!
!! R83NP_FS factors and solves an R83NP system.
!
!  Discussion:
!
!    The R83NP storage format is used for a tridiagonal matrix.
!    The subdiagonal   is in entries (1,2:N), 
!    the diagonal      is in entries (2,1:N), 
!    the superdiagonal is in entries (3,1:N-1). 
!
!    This algorithm requires that each diagonal entry be nonzero.
!    It does not use pivoting, and so can fail on systems that
!    are actually nonsingular.
!
!    The "R83NP" format used for this routine is different from the R83 format.
!    Here, we insist that the nonzero entries
!    for a given row now appear in the corresponding column of the
!    packed array.
!
!  Example:
!
!    Here is how an R83NP matrix of order 5 would be stored:
!
!       *  A21 A32 A43 A54
!      A11 A22 A33 A44 A55
!      A12 A23 A34 A45  *
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the linear system.
!
!    Input/output, real ( kind = rk ) A(3,N).
!    On input, the tridiagonal matrix.
!    On output, the data in these vectors has been overwritten
!    by factorization information.
!
!    Input, real ( kind = rk ) B(N), the right hand side of the linear system.
!
!    Output, real ( kind = rk ) X(N), the solution of the linear system.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(3,n)
  real ( kind = rk ) b(n)
  integer i
  real ( kind = rk ) x(n)
!
!  The diagonal entries can't be zero.
!
  do i = 1, n
    if ( a(2,i) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R83NP_FS - Fatal error!'
      write ( *, '(a,i8,a)' ) '  A(2,', i, ') = 0.'
      stop 1
    end if
  end do

  x(1:n) = b(1:n)

  do i = 2, n
    a(2,i) = a(2,i) - a(3,i-1) * a(1,i) / a(2,i-1)
    x(i)   = x(i)   - x(i-1)   * a(1,i) / a(2,i-1)
  end do

  x(n) = x(n) / a(2,n)
  do i = n - 1, 1, -1
    x(i) = ( x(i) - a(3,i) * x(i+1) ) / a(2,i)
  end do

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
!  Parameters:
!
!    None
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

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

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
