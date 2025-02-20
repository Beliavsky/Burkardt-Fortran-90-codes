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
subroutine r8ge_det ( n, a_lu, pivot, det )

!*****************************************************************************80
!
!! R8GE_DET: determinant of a matrix factored by R8GE_FA or R8GE_TRF.
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
!    29 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real ( kind = rk ) A_LU(N,N), the LU factors from R8GE_FA 
!    or R8GE_TRF.
!
!    Input, integer PIVOT(N), as computed by R8GE_FA or R8GE_TRF.
!
!    Output, real ( kind = rk ) DET, the determinant of the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a_lu(n,n)
  real ( kind = rk ) det
  integer i
  integer pivot(n)

  det = 1.0D+00

  do i = 1, n
    det = det * a_lu(i,i)
    if ( pivot(i) /= i ) then
      det = - det
    end if
  end do

  return
end
subroutine r8ge_fa ( n, a, pivot, info )

!*****************************************************************************80
!
!! R8GE_FA performs a LINPACK style PLU factorization of an R8GE matrix.
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
!    R8GE_FA is a simplified version of the LINPACK routine SGEFA.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input/output, real ( kind = rk ) A(N,N), the matrix to be factored.
!    On output, A contains an upper triangular matrix and the multipliers
!    which were used to obtain it.  The factorization can be written
!    A = L * U, where L is a product of permutation and unit lower
!    triangular matrices and U is upper triangular.
!
!    Output, integer PIVOT(N), a vector of pivot indices.
!
!    Output, integer INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,n)
  integer i
  integer info
  integer pivot(n)
  integer j
  integer k
  integer l
  real ( kind = rk ) t

  info = 0

  do k = 1, n - 1
!
!  Find L, the index of the pivot row.
!
    l = k
    do i = k + 1, n
      if ( abs ( a(l,k) ) < abs ( a(i,k) ) ) then
        l = i
      end if
    end do

    pivot(k) = l
!
!  If the pivot index is zero, the algorithm has failed.
!
    if ( a(l,k) == 0.0D+00 ) then
      info = k
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8GE_FA - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      stop 1
    end if
!
!  Interchange rows L and K if necessary.
!
    if ( l /= k ) then
      t      = a(l,k)
      a(l,k) = a(k,k)
      a(k,k) = t
    end if
!
!  Normalize the values that lie below the pivot entry A(K,K).
!
    a(k+1:n,k) = -a(k+1:n,k) / a(k,k)
!
!  Row elimination with column indexing.
!
    do j = k + 1, n

      if ( l /= k ) then
        t      = a(l,j)
        a(l,j) = a(k,j)
        a(k,j) = t
      end if

      a(k+1:n,j) = a(k+1:n,j) + a(k+1:n,k) * a(k,j)

    end do

  end do

  pivot(n) = n

  if ( a(n,n) == 0.0D+00 ) then
    info = n
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_FA - Fatal error!'
    write ( *, '(a,i8)' ) '  Zero pivot on step ', info
    stop 1
  end if

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
subroutine r8pbu_det ( n, mu, a_lu, det )

!*****************************************************************************80
!
!! R8PBU_DET computes the determinant of a matrix factored by R8PBU_FA.
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
!    29 October 1998
!
!  Author:
!
!    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer MU, the number of superdiagonals of the matrix.
!    MU must be at least 0 and no more than N-1.
!
!    Input, real ( kind = rk ) A_LU(MU+1,N), the LU factors from R8PBU_FA.
!
!    Output, real ( kind = rk ) DET, the determinant of the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer mu
  integer n

  real ( kind = rk ) a_lu(mu+1,n)
  real ( kind = rk ) det

  det = product ( a_lu(mu+1,1:n)**2 )

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

  a(1:mu+1,1:n) = 0.0D+00

  a(mu,  2:n) = -1.0D+00
  a(mu+1,1:n) =  +2.0D+00
 
  return
end
subroutine r8pbu_fa ( n, mu, a, info )

!*****************************************************************************80
!
!! R8PBU_FA factors an R8PBU matrix.
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
!    Once factored, linear systems A*x=b involving the matrix can be solved
!    by calling R8PBU_SL.  No pivoting is performed.  Pivoting is not necessary
!    for positive definite symmetric matrices.  If the matrix is not positive
!    definite, the algorithm may behave correctly, but it is also possible
!    that an illegal divide by zero will occur.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    31 October 1998
!
!  Author:
!
!    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer MU, the number of superdiagonals of the matrix.
!    MU must be at least 0, and no more than N-1.
!
!    Input/output, real ( kind = rk ) A(MU+1,N), the N by N matrix, stored 
!    in LINPACK positive definite symmetric band matrix storage.
!    On output, A contains information describing a factored form
!    of the matrix, that can be used to solve linear systems
!    A*x=b, using R8PBU_SL.
!
!    Output, integer INFO, singularity flag.
!    0, the matrix is nonsingular.
!    nonzero, the matrix is singular.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer mu
  integer n

  real ( kind = rk ) a(mu+1,n)
  integer ik
  integer info
  integer j
  integer jk
  integer k
  integer mm
  real ( kind = rk ) s

  info = 0

  do j = 1, n

    ik = mu + 1
    jk = max ( j - mu, 1 )
    mm = max ( mu + 2 - j, 1 )

    s = 0.0D+00

    do k = mm, mu

      a(k,j) = ( a(k,j) - sum ( a(ik:ik+k-mm-1,jk) * a(mm:k-1,j) ) ) &
        / a(mu+1,jk)

      s = s + a(k,j) * a(k,j)

      ik = ik - 1
      jk = jk + 1

    end do

    s = a(mu+1,j) - s

    if ( s <= 0.0D+00 ) then
      info = j
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8PBU_FA - Fatal error!'
      write ( *, '(a,i8)' ) '  Nonpositive pivot on step ', info
      stop 1
    end if

    a(mu+1,j) = sqrt ( s )

  end do

  return
end
subroutine r8pbu_indicator ( n, mu, a )

!*****************************************************************************80
!
!! R8PBU_INDICATOR sets up an R8PBU indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
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
!    14 January 2004
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
!    Output, real ( kind = rk ) A(MU+1,N), the R8PBU matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer mu
  integer n

  real ( kind = rk ) a(mu+1,n)
  integer fac
  integer i
  integer i4_log_10
  integer j

  fac = 10 ** ( i4_log_10 ( n ) + 1 )
!
!  Zero out the "junk" entries.
!
  do j = 1, mu
    do i = 1, mu + 1 - j
      a(i,j) = 0.0D+00
    end do
  end do
!
!  Set the meaningful values.
!
  do i = 1, n
    do j = i, min ( i + mu, n )
      a(mu+1+i-j,j) = real ( fac * i + j, kind = rk )
    end do
  end do

  return
end
subroutine r8pbu_ml ( n, mu, a_lu, x, b )

!*****************************************************************************80
!
!! R8PBU_ML multiplies an R8VEC times a matrix that was factored by R8PBU_FA.
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
!    31 October 1998
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
!    Input, integer MU, the number of superdiagonals of the matrix.
!    MU must be at least 0 and no more than N-1.
!
!    Input, real ( kind = rk ) A_LU(MU+1,N), the LU factors from R8PBU_FA.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = rk ) B(N), the product A * x.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer mu
  integer n

  real ( kind = rk ) a_lu(mu+1,n)
  real ( kind = rk ) b(n)
  integer i
  integer ilo
  integer j
  integer jhi
  integer k
  real ( kind = rk ) x(n)

  b(1:n) = x(1:n)
!
!  Multiply U * X = Y.
!
  do k = 1, n

    ilo = max ( 1, k - mu )
    do i = ilo, k - 1
      b(i) = b(i) + a_lu(mu+1+i-k,k) * b(k)
    end do

    b(k) = a_lu(mu+1,k) * b(k)

  end do
!
!  Multiply L * Y = B.
!
  do k = n, 1, -1

    jhi = min ( k + mu, n )
    do j = k + 1, jhi
      b(j) = b(j) + a_lu(mu+1+k-j,j) * b(k)
    end do

    b(k) = a_lu(mu+1,k) * b(k)

  end do

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
subroutine r8pbu_print ( n, mu, a, title )

!*****************************************************************************80
!
!! R8PBU_PRINT prints an R8PBU matrix.
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
!    07 May 2000
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
!    Input, integer MU, the upper (and lower) bandwidth.
!    MU must be nonnegative, and no greater than N-1.
!
!    Input, real ( kind = rk ) A(MU+1,N), the R8PBU matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer mu
  integer n

  real ( kind = rk ) a(mu+1,n)
  character ( len = * ) title

  call r8pbu_print_some ( n, mu, a, 1, 1, n, n, title )

  return
end
subroutine r8pbu_print_some ( n, mu, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8PBU_PRINT_SOME prints some of an R8PBU matrix.
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
!    21 March 2001
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
!    Input, integer MU, the upper (and lower) bandwidth.
!    MU must be nonnegative, and no greater than N-1.
!
!    Input, real ( kind = rk ) A(MU+1,N), the R8PBU matrix.
!
!    Input, integer ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: incx = 5
  integer mu
  integer n

  real ( kind = rk ) a(mu+1,n)
  real ( kind = rk ) aij
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

    write ( *, '(a,5a14)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2lo = max ( i2lo, j2lo - mu )
    i2hi = min ( ihi, n )
    i2hi = min ( i2hi, j2hi + mu )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( i <= j .and. j <= i + mu ) then
          aij = a(mu+1+i-j,j)
        else if ( i - mu <= j .and. j <= i ) then
          aij = a(mu+1+j-i,i)
        else
          aij = 0.0D+00
        end if

        if ( mu < i-j .or. mu < j-i ) then
          ctemp(j2) = '              '
        else
          write ( ctemp(j2), '(g14.6)' ) aij
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8pbu_random ( n, mu, seed, a )

!*****************************************************************************80
!
!! R8PBU_RANDOM randomizes an R8PBU matrix.
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
!    The matrix returned will be positive definite, but of limited
!    randomness.  The off diagonal elements are random values between
!    0 and 1, and the diagonal element of each row is selected to
!    ensure strict diagonal dominance.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 May 2003
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
!    Input/output, integer SEED, a seed for the random number 
!    generator.
!
!    Output, real ( kind = rk ) A(MU+1,N), the R8PBU matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer mu
  integer n

  real ( kind = rk ) a(mu+1,n)
  real ( kind = rk ) r8_uniform_01
  integer i
  integer j
  integer jhi
  integer jlo
  real ( kind = rk ) r
  integer seed
  real ( kind = rk ) sum2
!
!  Zero out the "junk" entries.
!
  do j = 1, mu
    a(1:mu+1-j,j) = 0.0D+00
  end do
!
!  Set the off diagonal values.
!
  do i = 1, n
    do j = i + 1, min ( i + mu, n )
      a(mu+1+i-j,j) = r8_uniform_01 ( seed )
    end do
  end do
!
!  Set the diagonal values.
!
  do i = 1, n

    sum2 = 0.0D+00

    jlo = max ( 1, i - mu )
    do j = jlo, i - 1
      sum2 = sum2 + abs ( a(mu+1+j-i,i) )
    end do

    jhi = min ( i + mu, n )
    do j = i + 1, jhi
      sum2 = sum2 + abs ( a(mu+1+i-j,j) )
    end do

    r = r8_uniform_01 ( seed )

    a(mu+1,i) = ( 1.0D+00 + r ) * ( sum2 + 0.01D+00 )

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
subroutine r8pbu_sl ( n, mu, a_lu, b )

!*****************************************************************************80
!
!! R8PBU_SL solves an R8PBU system factored by R8PBU_FA.
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
!    31 October 1998
!
!  Author:
!
!    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer MU, the number of superdiagonals of the matrix.
!    MU must be at least 0 and no more than N-1.
!
!    Input, real ( kind = rk ) A_LU(MU+1,N), the LU factors from R8PBU_FA.
!
!    Input/output, real ( kind = rk ) B(N).
!    On input, B contains the right hand side of the linear system
!    to be solved.
!    On output, B contains X, the solution vector.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer mu
  integer n

  real ( kind = rk ) a_lu(mu+1,n)
  real ( kind = rk ) b(n)
  integer i
  integer ilo
  integer k
!
!  Solve L * Y = B.
!
  do k = 1, n
    ilo = max ( 1, k - mu )
    b(k) = ( b(k) - sum ( b(ilo:k-1) * a_lu(mu+1+ilo-k:mu,k) ) ) &
      / a_lu(mu+1,k)
  end do
!
!  Solve U * X = Y.
!
  do k = n, 1, -1

    b(k) = b(k) / a_lu(mu+1,k)

    ilo = max ( 1, k - mu )
    do i = ilo, k - 1
      b(i) = b(i) - b(k) * a_lu(mu+1+i-k,k)
    end do

  end do

  return
end
subroutine r8pbu_sor ( n, mu, a, b, eps, itchk, itknt, itmax, omega, x )

!*****************************************************************************80
!
!! R8PBU_SOR uses SOR iteration to solve an R8PBU linear system.
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
!    A relaxation factor OMEGA may be used.
!
!    The iteration will proceed until a convergence test is met,
!    or the iteration limit is reached.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 October 1998
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
!    Input, integer MU, the number of superdiagonals in the 
!    matrix.  MU must be at least 0, and no more than N-1.
!
!    Input, real ( kind = rk ) A(MU+1,N), the R8PBU matrix.
!
!    Input, real ( kind = rk ) B(N), the right hand side of the system.
!
!    Input, real ( kind = rk ) EPS, convergence tolerance for the system. 
!    The vector b - A * x is computed every ITCHK iterations, and if the 
!    maximum entry of this vector is of norm less than EPS, the program
!    will return.
!
!    Input, integer ITCHK, the interval between convergence checks.
!    ITCHK steps will be taken before any check is made on whether the iteration
!    has converged.  ITCHK should be at least 1 and no greater
!    than ITMAX.
!
!    Output, integer ITKNT, the number of iterations taken.
!
!    Input, integer ITMAX, the maximum number of iterations 
!    allowed.  The program will return to the user if this many iterations 
!    are taken without convergence.
!
!    Input, real ( kind = rk ) OMEGA, the relaxation factor.  OMEGA must be
!    strictly between 0 and 2.  Use OMEGA = 1 for no relaxation, classical
!    Jacobi iteration.
!
!    Input/output, real ( kind = rk ) X(N).
!    On input, a starting vector for the iteration.
!    On output, the current approximation to the solution.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer mu
  integer n

  real ( kind = rk ) a(mu+1,n)
  real ( kind = rk ) b(n)
  real ( kind = rk ) eps
  real ( kind = rk ) err
  integer it
  integer itchk
  integer itknt
  integer itmax
  real ( kind = rk ) omega
  real ( kind = rk ) x(n)
  real ( kind = rk ) xtemp(n)

  if ( itchk <= 0 .or. itmax < itchk ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8PBU_SOR - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal ITCHK= ', itchk
    stop 1
  end if

  if ( itmax <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8PBU_SOR - Fatal error!'
    write ( *, '(a,i8)' ) '  Nonpositive ITMAX =', itmax
    stop 1
  end if

  if ( omega <= 0.0D+00 .or. 2.0D+00 <= omega ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8PBU_SOR - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Illegal value of OMEGA = ', omega
    stop 1
  end if 

  itknt = 0
!
!  Take ITCHK steps of the iteration before doing a convergence check.
!
  do while ( itknt <= itmax )

    do it = 1, itchk
!
!  Compute XTEMP(I) = B(I) + A(I,I) * X(I) - SUM ( J=1 to N ) A(I,J) * X(J).
!
      call r8pbu_mv ( n, n, mu, a, x, xtemp )

      xtemp(1:n) = x(1:n) + ( b(1:n) - xtemp(1:n) ) / a(mu+1,1:n)
!
!  Compute the next iterate as a weighted combination of the
!  old iterate and the just computed standard Jacobi iterate.
!
      if ( omega /= 1.0D+00 ) then
        xtemp(1:n) = ( 1.0D+00 - omega ) * x(1:n) + omega * xtemp(1:n)
      end if

      itknt = itknt + 1
!
!  Copy the new result into the old result vector.
!
      x(1:n) = xtemp(1:n)

    end do
!
!  Compute the maximum residual, the greatest entry in the vector
!  RESID(I) = B(I) - A(I,J) * X(J).
!
    call r8pbu_mv ( n, n, mu, a, x, xtemp )

    err = maxval ( abs ( b(1:n) - xtemp(1:n) ) )
!
!  Test to see if we can quit because of convergence,
!
    if ( err <= eps ) then
      return
    end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8PBU_SOR - Warning!'
  write ( *, '(a)' ) '  The iteration did not converge.'

  return
end
subroutine r8pbu_to_r8ge ( n, mu, a, b )

!*****************************************************************************80
!
!! R8PBU_TO_R8GE copies an R8PBU matrix to an R8GE matrix.
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
!    14 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrices.
!    N must be positive.
!
!    Input, integer MU, the upper bandwidth of A1.
!    MU must be nonnegative, and no greater than N-1.
!
!    Input, real ( kind = rk ) A(MU+1,N), the R8PBU matrix.
!
!    Output, real ( kind = rk ) B(N,N), the R8GE matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer mu
  integer n

  real ( kind = rk ) a(mu+1,n)
  real ( kind = rk ) b(n,n)
  integer i
  integer j

  do i = 1, n
    do j = 1, n
      if ( i <= j .and. j <= i + mu ) then
        b(i,j) = a(mu+1+i-j,j)
      else if ( i - mu <= j .and. j < i ) then
        b(i,j) = a(mu+1+j-i,i)
      else
        b(i,j) = 0.0D+00
      end if
    end do
  end do

  return
end
subroutine r8pbu_zeros ( n, mu, a )

!*****************************************************************************80
!
!! R8PBU_ZEROS zeroes an R8PBU matrix.
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
!    The matrix returned will be positive definite, but of limited
!    randomness.  The off diagonal elements are random values between
!    0 and 1, and the diagonal element of each row is selected to
!    ensure strict diagonal dominance.
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
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer MU, the number of superdiagonals in the matrix.
!    MU must be at least 0 and no more than N-1.
!
!    Output, real ( kind = rk ) A(MU+1,N), the R8PBU matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer mu
  integer n

  real ( kind = rk ) a(mu+1,n)

  a(1:mu+1,1:n) = 0.0D+00
 
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
subroutine r8vec_print_some ( n, a, max_print, title )

!*****************************************************************************80
!
!! R8VEC_PRINT_SOME prints "some" of an R8VEC.
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vector, is no more than MAX_PRINT, then
!    the entire vector is printed, one entry per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries of the vector.
!
!    Input, real ( kind = rk ) A(N), the vector to be printed.
!
!    Input, integer MAX_PRINT, the maximum number of lines 
!    to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  integer i
  integer max_print
  character ( len = * ) title

  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  if ( n <= max_print ) then

    if ( all ( a(1:n) == aint ( a(1:n) ) ) ) then
      do i = 1, n
        write ( *, '(i8,2x,i8)' ) i, int ( a(i) )
      end do
    else if ( all ( abs ( a(1:n) ) < 1000000.0D+00 ) ) then
      do i = 1, n
        write ( *, '(i8,2x,f14.6)' ) i, a(i)
      end do
    else
      do i = 1, n
        write ( *, '(i8,2x,g14.6)' ) i, a(i)
      end do
    end if

  else if ( 3 <= max_print ) then

    if ( all ( a(1:max_print-2) == aint ( a(1:max_print-2) ) ) ) then
      do i = 1, max_print - 2
        write ( *, '(i8,2x,i8)' ) i, int ( a(i) )
      end do
    else if ( all ( abs ( a(1:max_print-2) ) < 1000000.0D+00 ) ) then
      do i = 1, max_print - 2
        write ( *, '(i8,2x,f14.6)' ) i, a(i)
      end do
    else
      do i = 1, max_print - 2
        write ( *, '(i8,2x,g14.6)' ) i, a(i)
      end do
    end if

    write ( *, '(a)' ) '......  ..............'
    i = n

    if ( a(i) == real ( int ( a(i) ), kind = rk ) ) then
      write ( *, '(i8,2x,i8)' ) i, int ( a(i) )
    else if (  abs ( a(i) ) < 1000000.0D+00 ) then
      write ( *, '(i8,2x,f14.6)' ) i, a(i)
    else
      write ( *, '(i8,2x,g14.6)' ) i, a(i)
    end if

  else

    if ( all ( a(1:max_print-1) == aint ( a(1:max_print-1) ) ) ) then
      do i = 1, max_print - 1
        write ( *, '(i8,2x,i8)' ) i, int ( a(i) )
      end do
    else if ( all ( abs ( a(1:max_print-1) ) < 1000000.0D+00 ) ) then
      do i = 1, max_print - 1
        write ( *, '(i8,2x,f14.6)' ) i, a(i)
      end do
    else
      do i = 1, max_print - 1
        write ( *, '(i8,2x,g14.6)' ) i, a(i)
      end do
    end if

    i = max_print

    if ( a(i) == aint ( a(i) ) ) then
      write ( *, '(i8,2x,i8,a)' ) i, int ( a(i) ), '...more entries...'
    else if (  abs ( a(i) ) < 1000000.0D+00 ) then
      write ( *, '(i8,2x,f14.6,a)' ) i, a(i), '...more entries...'
    else
      write ( *, '(i8,2x,g14.6,a)' ) i, a(i), '...more entries...'
    end if

  end if

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
subroutine r8vec2_print_some ( n, x1, x2, max_print, title )

!*****************************************************************************80
!
!! R8VEC2_PRINT_SOME prints "some" of a pair of R8VEC's.
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vectors, is no more than MAX_PRINT, then
!    the entire vectors are printed, one entry of each per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries of the vectors.
!
!    Input, real ( kind = rk ) X1(N), X2(N), the vector to be printed.
!
!    Input, integer MAX_PRINT, the maximum number of lines 
!    to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer max_print
  character ( len = * ) title
  real ( kind = rk ) x1(n)
  real ( kind = rk ) x2(n)

  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  if ( n <= max_print ) then

    do i = 1, n
      write ( *, '(i8,2x,g14.6,2x,g14.6)' ) i, x1(i), x2(i)
    end do

  else if ( 3 <= max_print ) then

    do i = 1, max_print - 2
      write ( *, '(i8,2x,g14.6,2x,g14.6)' ) i, x1(i), x2(i)
    end do
    write ( *, '(a)' ) '......  ..............  ..............'
    i = n
    write ( *, '(i8,2x,g14.6,2x,g14.6)' ) i, x1(i), x2(i)

  else

    do i = 1, max_print - 1
      write ( *, '(i8,2x,g14.6,2x,g14.6)' ) i, x1(i), x2(i)
    end do
    i = max_print
    write ( *, '(i8,2x,g14.6,2x,g14.6,2x,a)' ) i, x1(i), x2(i), &
      '...more entries...'

  end if

  return
end

