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

  a(1:n,1:ndiag) = 0.0D+00

  a(1:n,  1) =  2.0D+00
  a(1:n-1,2) = -1.0D+00
 
  return
end
subroutine r8sd_indicator ( n, ndiag, offset, a )

!*****************************************************************************80
!
!! R8SD_INDICATOR sets up an R8SD indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
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
!    14 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer NDIAG, the number of diagonals that are stored.
!    NDIAG must be at least 1 and no more than N.
!
!    Input, integer OFFSET(NDIAG), the offsets for the diagonal
!    storage.
!
!    Output, real ( kind = rk ) A(N,NDIAG), the R8SD matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer ndiag

  real ( kind = rk ) a(n,ndiag)
  integer diag
  integer fac
  integer i
  integer i4_log_10
  integer j
  integer offset(ndiag)

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  do i = 1, n
    do diag = 1, ndiag
      j = i + offset(diag)
      if ( 1 <= j .and. j <= n ) then
        a(i,diag) = real ( fac * i + j, kind = rk )
      else
        a(i,diag) = 0.0D+00
      end if
    end do
  end do

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
subroutine r8sd_print ( n, ndiag, offset, a, title )

!*****************************************************************************80
!
!! R8SD_PRINT prints an R8SD matrix.
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
!    07 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of columns of the matrix.
!
!    Input, integer NDIAG, the number of diagonals of the matrix
!    that are stored in the array.
!    NDIAG must be at least 1, and no more than N.
!
!    Input, integer OFFSET(NDIAG), the offsets for the 
!    diagonal storage.
!
!    Input, real ( kind = rk ) A(N,NDIAG), the R8SD matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer ndiag

  real ( kind = rk ) a(n,ndiag)
  integer offset(ndiag)
  character ( len = * ) title

  call r8sd_print_some ( n, ndiag, offset, a, 1, 1, n, n, title )

  return
end
subroutine r8sd_print_some ( n, ndiag, offset, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8SD_PRINT_SOME prints some of an R8SD matrix.
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
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of columns of the matrix.
!
!    Input, integer NDIAG, the number of diagonals of the matrix
!    that are stored in the array.
!    NDIAG must be at least 1, and no more than N.
!
!    Input, integer OFFSET(NDIAG), the offsets for the diagonal
!    storage.
!
!    Input, real ( kind = rk ) A(N,NDIAG), the R8SD matrix.
!
!    Input, integer ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: incx = 5
  integer n
  integer ndiag

  real ( kind = rk ) a(n,ndiag)
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
  integer jdiag
  integer jhi
  integer jlo
  integer off
  integer offset(ndiag)
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
    i2hi = min ( ihi, n )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        aij = 0.0D+00
        off = j - i
        do jdiag = 1, ndiag
          if ( off == offset(jdiag) ) then
            aij = a(i,jdiag)
          else if ( off == - offset(jdiag) ) then
            aij = a(j,jdiag)
          end if
        end do

        write ( ctemp(j2), '(g14.6)' ) aij

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8sd_random ( n, ndiag, offset, seed, a )

!*****************************************************************************80
!
!! R8SD_RANDOM randomizes an R8SD matrix.
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
!    17 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer NDIAG, the number of diagonals that are stored.
!    NDIAG must be at least 1 and no more than N.
!
!    Input, integer OFFSET(NDIAG), the offsets for the diagonal
!    storage.
!
!    Input/output, integer SEED, a seed for the random number
!    generator.
!
!    Output, real ( kind = rk ) A(N,NDIAG), the R8SD matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer ndiag

  real ( kind = rk ) a(n,ndiag)
  real ( kind = rk ) r8_uniform_01
  integer i
  integer j
  integer jj
  integer offset(ndiag)
  integer seed

  do i = 1, n
    do j = 1, ndiag
      jj = i + offset(j)
      if ( 1 <= jj .and. jj <= n ) then
        a(i,j) = r8_uniform_01 ( seed )
      else
        a(i,j) = 0.0D+00
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
!
!    Input, integer N, the number of columns of the matrix.
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
subroutine r8sd_to_r8ge ( n, ndiag, offset, a, b )

!*****************************************************************************80
!
!! R8SD_TO_R8GE copies an R8SD matrix to an R8GE matrix.
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
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
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
!    30 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer NDIAG, the number of diagonals that are stored.
!    NDIAG must be at least 1 and no more than N.
!
!    Input, integer OFFSET(NDIAG), the offsets for the diagonal
!    storage.
!
!    Input, real ( kind = rk ) A(N,NDIAG), the R8SD matrix.
!
!    Output, real ( kind = rk ) B(N,N), the R8GE matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer ndiag

  real ( kind = rk ) a(n,ndiag)
  real ( kind = rk ) b(n,n)
  integer i
  integer j
  integer jj
  integer offset(ndiag)

  b(1:n,1:n) = 0.0D+00

  do i = 1, n
    do j = 1, ndiag
      jj = i + offset(j)
      if ( 1 <= jj .and. jj <= n ) then
        b(i,jj) = a(i,j)
        if ( i /= jj ) then
          b(jj,i) = a(i,j)
        end if
      end if
    end do
  end do

  return
end
subroutine r8sd_zeros ( n, ndiag, offset, a )

!*****************************************************************************80
!
!! R8SD_ZEROS zeroes an R8SD matrix.
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
!    26 January 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer NDIAG, the number of diagonals that are stored.
!    NDIAG must be at least 1 and no more than N.
!
!    Input, integer OFFSET(NDIAG), the offsets for the diagonal
!    storage.
!
!    Output, real ( kind = rk ) A(N,NDIAG), the R8SD matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer ndiag

  real ( kind = rk ) a(n,ndiag)
  integer offset(ndiag)

  a(1:n,1:ndiag) = 0.0D+00

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
!    19 December 2001
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

    do i = 1, n
      write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)
    end do

  else if ( 3 <= max_print ) then

    do i = 1, max_print - 2
      write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)
    end do
    write ( *, '(a)' ) '  ........  ..............'
    i = n
    write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)

  else

    do i = 1, max_print - 1
      write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)
    end do
    i = max_print
    write ( *, '(2x,i8,a,1x,g14.6,2x,a)' ) i, ':', a(i), '...more entries...'

  end if

  return
end

