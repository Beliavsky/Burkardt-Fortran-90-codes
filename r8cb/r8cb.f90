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
subroutine r8cb_det ( n, ml, mu, a_lu, det )

!*****************************************************************************80
!
!! R8CB_DET computes the determinant of an R8CB matrix factored by R8CB_NP_FA.
!
!  Discussion:
!
!    The R8CB storage format is used for a compact banded matrix.
!    It is assumed that the matrix has lower and upper bandwidths ML and MU,
!    respectively.  The matrix is stored in a way similar to that used
!    by LINPACK and LAPACK for a general banded matrix, except that in
!    this mode, no extra rows are set aside for possible fillin during pivoting.
!    Thus, this storage mode is suitable if you do not intend to factor
!    the matrix, or if you can guarantee that the matrix can be factored
!    without pivoting.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    25 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input, real ( kind = rk ) A_LU(ML+MU+1,N), the LU factors from R8CB_NP_FA.
!
!    Output, real ( kind = rk ) DET, the determinant of the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer n

  real ( kind = rk ) a_lu(ml+mu+1,n)
  real ( kind = rk ) det

  det = product ( a_lu(mu+1,1:n) )

  return
end
subroutine r8cb_dif2 ( m, n, ml, mu, a )

!*****************************************************************************80
!
!! R8CB_DIF2 sets up an R8CB second difference matrix.
!
!  Discussion:
!
!    The R8CB storage format is used for a compact banded matrix.
!    It is assumed that the matrix has lower and upper bandwidths ML and MU,
!    respectively.  The matrix is stored in a way similar to that used
!    by LINPACK and LAPACK for a general banded matrix, except that in
!    this mode, no extra rows are set aside for possible fillin during pivoting.
!    Thus, this storage mode is suitable if you do not intend to factor
!    the matrix, or if you can guarantee that the matrix can be factored
!    without pivoting.
!
!    The original M by N matrix is "collapsed" downward, so that diagonals
!    become rows of the storage array, while columns are preserved.  The
!    collapsed array is logically ML+MU+1 by N.  
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 July 2016
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
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than min(M,N)-1.
!
!    Output, real ( kind = rk ) A(ML+MU+1,N), the R8CB matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer m
  integer n

  real ( kind = rk ) a(ml+mu+1,n)
  integer diag
  integer i
  integer j

  a(1:ml+mu+1,1:n) = 0.0D+00

  do j = 1, n

    do diag = 1, ml + mu + 1

      i = diag + j - mu - 1

      if ( i == j ) then
        a(diag,j) = 2.0D+00
      else if ( i == j + 1 .or. i == j - 1 ) then
        a(diag,j) = -1.0D+00
      end if

    end do
  end do

  return
end
subroutine r8cb_indicator ( m, n, ml, mu, a )

!*****************************************************************************80
!
!! R8CB_INDICATOR sets up an R8CB indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R8CB storage format is used for a compact banded matrix.
!    It is assumed that the matrix has lower and upper bandwidths ML and MU,
!    respectively.  The matrix is stored in a way similar to that used
!    by LINPACK and LAPACK for a general banded matrix, except that in
!    this mode, no extra rows are set aside for possible fillin during pivoting.
!    Thus, this storage mode is suitable if you do not intend to factor
!    the matrix, or if you can guarantee that the matrix can be factored
!    without pivoting.
!
!    The original M by N matrix is "collapsed" downward, so that diagonals
!    become rows of the storage array, while columns are preserved.  The
!    collapsed array is logically ML+MU+1 by N.  
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 January 2004
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
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than min(M,N)-1.
!
!    Output, real ( kind = rk ) A(ML+MU+1,N), the R8CB matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer m
  integer n

  real ( kind = rk ) a(ml+mu+1,n)
  integer diag
  integer fac
  integer i
  integer i4_log_10
  integer j
  integer k
  real ( kind = rk ) value

  a(1:ml+mu+1,1:n) = 0.0D+00

  fac = 10 ** ( i4_log_10 ( n ) + 1 )
  k = 0

  do j = 1, n
    do diag = 1, ml + mu + 1

      i = diag + j - mu - 1

      if ( 1 <= i .and. i <= m .and. i - ml <= j .and. j <= i + mu ) then
        value = real ( fac * i + j, kind = rk )
      else
        k = k + 1
        value = - real ( k, kind = rk )
      end if

      a(diag,j) = value

    end do
  end do


  return
end
subroutine r8cb_ml ( n, ml, mu, a_lu, x, b, job )

!*****************************************************************************80
!
!! R8CB_ML computes A * x or A' * X, using R8CB_NP_FA factors.
!
!  Discussion:
!
!    The R8CB storage format is used for a compact banded matrix.
!    It is assumed that the matrix has lower and upper bandwidths ML and MU,
!    respectively.  The matrix is stored in a way similar to that used
!    by LINPACK and LAPACK for a general banded matrix, except that in
!    this mode, no extra rows are set aside for possible fillin during pivoting.
!    Thus, this storage mode is suitable if you do not intend to factor
!    the matrix, or if you can guarantee that the matrix can be factored
!    without pivoting.
!
!    It is assumed that R8CB_NP_FA has overwritten the original matrix
!    information by LU factors.  R8CB_ML is able to reconstruct the
!    original matrix from the LU factor data.
!
!    R8CB_ML allows the user to check that the solution of a linear
!    system is correct, without having to save an unfactored copy
!    of the matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input, real ( kind = rk ) A_LU(ML+MU+1,N), the LU factors from R8CB_NP_FA.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied.
!
!    Output, real ( kind = rk ) B(N), the result of the multiplication.
!
!    Input, integer JOB, specifies the operation to be done:
!    JOB = 0, compute A * x.
!    JOB nonzero, compute A' * x.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer n

  real ( kind = rk ) a_lu(ml+mu+1,n)
  real ( kind = rk ) b(n)
  integer i
  integer ihi
  integer ilo
  integer j
  integer jhi
  integer job
  real ( kind = rk ) x(n)

  b(1:n) = x(1:n)

  if ( job == 0 ) then
!
!  Y = U * X.
!
    do j = 1, n
      ilo = max ( 1, j - mu )
      do i = ilo, j - 1
        b(i) = b(i) + a_lu(i-j+mu+1,j) * b(j)
      end do
      b(j) = a_lu(j-j+mu+1,j) * b(j)
    end do
!
!  B = PL * Y = PL * U * X = A * x.
!
    do j = n - 1, 1, -1
      ihi = min ( n, j + ml )
      b(j+1:ihi) = b(j+1:ihi) - a_lu(mu+2:ihi-j+mu+1,j) * b(j)
    end do

  else
!
!  Y = ( PL )' * X.
!
    do j = 1, n - 1

      ihi = min ( n, j + ml )
      do i = j + 1, ihi
        b(j) = b(j) - b(i) * a_lu(i-j+mu+1,j)
      end do

    end do
!
!  B = U' * Y = ( PL * U )' * X = A' * X.
!
    do i = n, 1, -1
      jhi = min ( n, i + mu )
      do j = i + 1, jhi
        b(j) = b(j) + b(i) * a_lu(i-j+mu+1,j)
      end do
      b(i) = b(i) * a_lu(i-i+mu+1,i)
    end do

  end if

  return
end
subroutine r8cb_mtv ( m, n, ml, mu, a, x, b )

!*****************************************************************************80
!
!! R8CB_MTV computes b=A'*x, where A is an R8CB matrix.
!
!  Discussion:
!
!    The R8CB storage format is used for a compact banded matrix.
!    It is assumed that the matrix has lower and upper bandwidths ML and MU,
!    respectively.  The matrix is stored in a way similar to that used
!    by LINPACK and LAPACK for a general banded matrix, except that in
!    this mode, no extra rows are set aside for possible fillin during pivoting.
!    Thus, this storage mode is suitable if you do not intend to factor
!    the matrix, or if you can guarantee that the matrix can be factored
!    without pivoting.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the order of the matrix.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input, real ( kind = rk ) A(ML+MU+1,N), the R8CB matrix.
!
!    Input, real ( kind = rk ) X(M), the vector to be multiplied by A.
!
!    Output, real ( kind = rk ) B(N), the product A'*x.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer ml
  integer mu
  integer n

  real ( kind = rk ) a(ml+mu+1,n)
  real ( kind = rk ) b(n)
  integer i
  integer j
  integer jhi
  integer jlo
  real ( kind = rk ) x(m)

  b(1:n) = 0.0D+00

  do i = 1, m
    jlo = max ( 1, i - ml )
    jhi = min ( n, i + mu )
    do j = jlo, jhi
      b(j) = b(j) + x(i) * a(i-j+mu+1,j)
    end do
  end do

  return
end
subroutine r8cb_mv ( m, n, ml, mu, a, x, b )

!*****************************************************************************80
!
!! R8CB_MV computes b=A*x, where A is an R8CB matrix.
!
!  Discussion:
!
!    The R8CB storage format is used for a compact banded matrix.
!    It is assumed that the matrix has lower and upper bandwidths ML and MU,
!    respectively.  The matrix is stored in a way similar to that used
!    by LINPACK and LAPACK for a general banded matrix, except that in
!    this mode, no extra rows are set aside for possible fillin during pivoting.
!    Thus, this storage mode is suitable if you do not intend to factor
!    the matrix, or if you can guarantee that the matrix can be factored
!    without pivoting.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 July 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the order of the matrix.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input, real ( kind = rk ) A(ML+MU+1,N), the R8CB matrix.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = rk ) B(M), the product A * x.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer ml
  integer mu
  integer n

  real ( kind = rk ) a(ml+mu+1,n)
  real ( kind = rk ) b(m)
  integer i
  integer j
  integer jhi
  integer jlo
  real ( kind = rk ) x(n)

  b(1:m) = 0.0D+00

  do i = 1, m
    jlo = max ( 1, i - ml )
    jhi = min ( n, i + mu )
    do j = jlo, jhi
      b(i) = b(i) + a(i-j+mu+1,j) * x(j)
    end do
  end do

  return
end
subroutine r8cb_np_fa ( n, ml, mu, a, info )

!*****************************************************************************80
!
!! R8CB_NP_FA factors an R8CB matrix by Gaussian elimination.
!
!  Discussion:
!
!    The R8CB storage format is appropriate for a compact banded matrix.
!    It is assumed that the matrix has lower and upper bandwidths ML and MU,
!    respectively.  The matrix is stored in a way similar to that used
!    by LINPACK and LAPACK for a general banded matrix, except that in
!    this mode, no extra rows are set aside for possible fillin during pivoting.
!    Thus, this storage mode is suitable if you do not intend to factor
!    the matrix, or if you can guarantee that the matrix can be factored
!    without pivoting.
!
!    R8CB_NP_FA is a version of the LINPACK routine R8GBFA, modifed to use
!    no pivoting, and to be applied to the R8CB compressed band matrix storage
!    format.  It will fail if the matrix is singular, or if any zero
!    pivot is encountered.
!
!    If R8CB_NP_FA successfully factors the matrix, R8CB_NP_SL may be called
!    to solve linear systems involving the matrix.
!
!    The matrix is stored in a compact version of LINPACK general
!    band storage, which does not include the fill-in entires.
!    The following program segment will store the entries of a banded
!    matrix in the compact format used by this routine:
!
!      m = mu+1
!      do j = 1, n
!        i1 = max ( 1, j - mu )
!        i2 = min ( n, j + ml )
!        do i = i1, i2
!          k = i-j+m
!          a(k,j) = afull(i,j)
!        end do
!      end do
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input/output, real ( kind = rk ) A(ML+MU+1,N), the compact band matrix.
!    On input, the coefficient matrix of the linear system.
!    On output, the LU factors of the matrix.
!
!    Output, integer INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer n

  real ( kind = rk ) a(ml+mu+1,n)
  integer info
  integer j
  integer ju
  integer k
  integer lm
  integer m
  integer mm
!
!  The value of M is MU + 1 rather than ML + MU + 1.
!
  m = mu + 1
  info = 0
  ju = 0

  do k = 1, n - 1
!
!  If our pivot entry A(MU+1,K) is zero, then we must give up.
!
    if ( a(m,k) == 0.0D+00 ) then
      info = k
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8CB_NP_FA - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      stop 1
    end if
!
!  LM counts the number of nonzero elements that lie below the current
!  diagonal entry, A(K,K).
!
!  Multiply the LM entries below the diagonal by -1/A(K,K), turning
!  them into the appropriate "multiplier" terms in the L matrix.
!
    lm = min ( ml, n - k )
    a(m+1:m+lm,k) = - a(m+1:m+lm,k) / a(m,k)
!
!  MM points to the row in which the next entry of the K-th row is, A(K,J).
!  We then add L(I,K)*A(K,J) to A(I,J) for rows I = K+1 to K+LM.
!
    ju = max ( ju, mu + k )
    ju = min ( ju, n )
    mm = m

    do j = k + 1, ju
      mm = mm - 1
      a(mm+1:mm+lm,j) = a(mm+1:mm+lm,j) + a(mm,j) * a(m+1:m+lm,k)
    end do

  end do

  if ( a(m,n) == 0.0D+00 ) then
    info = n
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8CB_NP_FA - Fatal error!'
    write ( *, '(a,i8)' ) '  Zero pivot on step ', info
    stop 1
  end if

  return
end
subroutine r8cb_np_sl ( n, ml, mu, a_lu, b, job )

!*****************************************************************************80
!
!! R8CB_NP_SL solves an R8CB system factored by R8CB_NP_FA.
!
!  Discussion:
!
!    The R8CB storage format is used for a compact banded matrix.
!    It is assumed that the matrix has lower and upper bandwidths ML and MU,
!    respectively.  The matrix is stored in a way similar to that used
!    by LINPACK and LAPACK for a general banded matrix, except that in
!    this mode, no extra rows are set aside for possible fillin during pivoting.
!    Thus, this storage mode is suitable if you do not intend to factor
!    the matrix, or if you can guarantee that the matrix can be factored
!    without pivoting.
!
!    R8CB_NP_SL can also solve the related system A' * x = b.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input, real ( kind = rk ) A_LU(ML+MU+1,N), the LU factors from R8CB_NP_FA.
!
!    Input/output, real ( kind = rk ) B(N).
!    On input, B contains the right hand side of the linear system, B.
!    On output, B contains the solution of the linear system, X.
!
!    Input, integer JOB.
!    If JOB is zero, the routine will solve A * x = b.
!    If JOB is nonzero, the routine will solve A' * x = b.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer n

  real ( kind = rk ) a_lu(ml+mu+1,n)
  real ( kind = rk ) b(n)
  integer job
  integer k
  integer la
  integer lb
  integer lm
  integer m
!
!  The value of M is ML + 1, rather than MU + ML + 1.
!
  m = mu + 1
!
!  Solve A * x = b.
!
  if ( job == 0 ) then
!
!  Solve PL * Y = B.
!
    if ( 0 < ml ) then
      do k = 1, n - 1
        lm = min ( ml, n - k )
        b(k+1:k+lm) = b(k+1:k+lm) + b(k) * a_lu(m+1:m+lm,k)
      end do
    end if
!
!  Solve U * X = Y.
!
    do k = n, 1, -1

      b(k) = b(k) / a_lu(m,k)
      lm = min ( k, m ) - 1
      la = m - lm
      lb = k - lm

      b(lb:lb+lm-1) = b(lb:lb+lm-1) - b(k) * a_lu(la:la+lm-1,k)

    end do
!
!  Solve A' * X = B.
!
  else
!
!  Solve U' * Y = B.
!
    do k = 1, n
      lm = min ( k, m ) - 1
      la = m - lm
      lb = k - lm

      b(k) = ( b(k) - sum ( a_lu(la:la+lm-1,k) * b(lb:lb+lm-1) ) ) &
        / a_lu(m,k)

    end do
!
!  Solve ( PL )' * X = Y.
!
    if ( 0 < ml ) then

      do k = n - 1, 1, -1
        lm = min ( ml, n - k )
        b(k) = b(k) + sum ( a_lu(m+1:m+lm,k) * b(k+1:k+lm) )
      end do

    end if

  end if

  return
end
subroutine r8cb_print ( m, n, ml, mu, a, title )

!*****************************************************************************80
!
!! R8CB_PRINT prints an R8CB matrix.
!
!  Discussion:
!
!    The R8CB storage format is used for a compact banded matrix.
!    It is assumed that the matrix has lower and upper bandwidths ML and MU,
!    respectively.  The matrix is stored in a way similar to that used
!    by LINPACK and LAPACK for a general banded matrix, except that in
!    this mode, no extra rows are set aside for possible fillin during pivoting.
!    Thus, this storage mode is suitable if you do not intend to factor
!    the matrix, or if you can guarantee that the matrix can be factored
!    without pivoting.
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
!    Input, integer M, N, the number of rows and columns of 
!    the matrix.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than min(M,N)-1.
!
!    Input, real ( kind = rk ) A(ML+MU+1,N), the R8CB matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer n

  real ( kind = rk ) a(ml+mu+1,n)
  integer m
  character ( len = * ) title

  call r8cb_print_some ( m, n, ml, mu, a, 1, 1, m, n, title )

  return
end
subroutine r8cb_print_some ( m, n, ml, mu, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8CB_PRINT_SOME prints some of an R8CB matrix.
!
!  Discussion:
!
!    The R8CB storage format is used for a compact banded matrix.
!    It is assumed that the matrix has lower and upper bandwidths ML and MU,
!    respectively.  The matrix is stored in a way similar to that used
!    by LINPACK and LAPACK for a general banded matrix, except that in
!    this mode, no extra rows are set aside for possible fillin during pivoting.
!    Thus, this storage mode is suitable if you do not intend to factor
!    the matrix, or if you can guarantee that the matrix can be factored
!    without pivoting.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    13 January 2004
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
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than min(M,N)-1.
!
!    Input, real ( kind = rk ) A(ML+MU+1,N), the R8CB matrix.
!
!    Input, integer ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: incx = 5
  integer ml
  integer mu
  integer n

  real ( kind = rk ) a(ml+mu+1,n)
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
    i2hi = min ( ihi, m )
    i2hi = min ( i2hi, j2hi + ml )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( ml < i - j .or. mu < j - i ) then
          ctemp(j2) = '              '
        else
          write ( ctemp(j2), '(g14.6)' ) a(i-j+mu+1,j)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8cb_random ( m, n, ml, mu, seed, a )

!*****************************************************************************80
!
!! R8CB_RANDOM randomizes an R8CB matrix.
!
!  Discussion:
!
!    The R8CB storage format is used for a compact banded matrix.
!    It is assumed that the matrix has lower and upper bandwidths ML and MU,
!    respectively.  The matrix is stored in a way similar to that used
!    by LINPACK and LAPACK for a general banded matrix, except that in
!    this mode, no extra rows are set aside for possible fillin during pivoting.
!    Thus, this storage mode is suitable if you do not intend to factor
!    the matrix, or if you can guarantee that the matrix can be factored
!    without pivoting.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    24 July 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the order of the matrix.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input/output, integer SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = rk ) A(ML+MU+1,N), the R8CB matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer ml
  integer mu
  integer n

  real ( kind = rk ) a(ml+mu+1,n)
  real ( kind = rk ) r8_uniform_01
  integer i
  integer ihi
  integer ilo
  integer j
  integer seed

  a(1:ml+mu+1,1:n) = 0.0D+00
!
!  Set the entries that correspond to matrix elements.
!
  do j = 1, n

    ilo = max ( 1, j - mu )
    ihi = min ( m, j + ml )

    do i = ilo, ihi
      a(i-j+mu+1,j) = r8_uniform_01 ( seed )
    end do

  end do

  return
end
subroutine r8cb_to_r8vec ( m, n, ml, mu, a, x )

!*****************************************************************************80
!
!! R8CB_TO_R8VEC copies an R8CB matrix to an R8VEC.
!
!  Discussion:
!
!    In C++ and FORTRAN, this routine is not really needed.  In MATLAB,
!    a data item carries its dimensionality implicitly, and so cannot be
!    regarded sometimes as a vector and sometimes as an array.
!
!    The R8CB storage format is used for a compact banded matrix.
!    It is assumed that the matrix has lower and upper bandwidths ML and MU,
!    respectively.  The matrix is stored in a way similar to that used
!    by LINPACK and LAPACK for a general banded matrix, except that in
!    this mode, no extra rows are set aside for possible fillin during pivoting.
!    Thus, this storage mode is suitable if you do not intend to factor
!    the matrix, or if you can guarantee that the matrix can be factored
!    without pivoting.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 July 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns in 
!    the array.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!
!    Input, real ( kind = rk ) A(ML+MU+1,N), the array to be copied.
!
!    Output, real ( kind = rk ) X((ML+MU+1)*N), the vector.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer ml
  integer mu
  integer n

  real ( kind = rk ) a(ml+mu+1,n)
  integer i
  integer ihi
  integer ilo
  integer j
  real ( kind = rk ) x((ml+mu+1)*n)

  x(1:(ml+mu+1)*n) = 0.0D+00

  do j = 1, n

    ilo = max ( mu + 2 - j, 1 )
    ihi = mu + 1 + min ( ml, m - j )
    do i = ilo, ihi
      x(i+(j-1)*(ml+mu+1)) = a(i,j)
    end do

  end do

  return
end
subroutine r8cb_to_r8ge ( m, n, ml, mu, a, b )

!*****************************************************************************80
!
!! R8CB_TO_R8GE copies an R8CB matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R8CB storage format is used for a compact banded matrix.
!    It is assumed that the matrix has lower and upper bandwidths ML and MU,
!    respectively.  The matrix is stored in a way similar to that used
!    by LINPACK and LAPACK for a general banded matrix, except that in
!    this mode, no extra rows are set aside for possible fillin during pivoting.
!    Thus, this storage mode is suitable if you do not intend to factor
!    the matrix, or if you can guarantee that the matrix can be factored
!    without pivoting.
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
!    16 July 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the order of the matrices.
!
!    Input, integer ML, MU, the lower and upper bandwidths of A.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input, real ( kind = rk ) A(ML+MU+1,N), the R8CB matrix.
!
!    Output, real ( kind = rk ) B(M,N), the R8GE matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer ml
  integer mu
  integer n

  real ( kind = rk ) a(ml+mu+1,n)
  real ( kind = rk ) b(m,n)
  integer i
  integer j

  b(1:m,1:n) = 0.0D+00

  do i = 1, m
    do j = 1, n
      if ( j - mu <= i .and. i <= j + ml ) then
        b(i,j) = a(mu+1+i-j,j)
      end if
    end do
  end do

  return
end
subroutine r8cb_zeros ( n, ml, mu, a )

!*****************************************************************************80
!
!! R8CB_ZEROS zeroes an R8CB matrix.
!
!  Discussion:
!
!    The R8CB storage format is used for a compact banded matrix.
!    It is assumed that the matrix has lower and upper bandwidths ML and MU,
!    respectively.  The matrix is stored in a way similar to that used
!    by LINPACK and LAPACK for a general banded matrix, except that in
!    this mode, no extra rows are set aside for possible fillin during pivoting.
!    Thus, this storage mode is suitable if you do not intend to factor
!    the matrix, or if you can guarantee that the matrix can be factored
!    without pivoting.
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
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Output, real ( kind = rk ) A(ML+MU+1,N), the R8CB matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer n

  real ( kind = rk ) a(ml+mu+1,n)

  a(1:ml+mu+1,1:n) = 0.0D+00

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
subroutine r8ge_np_fa ( n, a, info )

!*****************************************************************************80
!
!! R8GE_NP_FA factors an R8GE matrix by nonpivoting Gaussian elimination.
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
!    R8GE_NP_FA is a version of the LINPACK routine SGEFA, but uses no
!    pivoting.  It will fail if the matrix is singular, or if any zero
!    pivot is encountered.
!
!    If R8GE_NP_FA successfully factors the matrix, R8GE_NP_SL may be called
!    to solve linear systems involving the matrix.
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
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real ( kind = rk ) A(N,N).
!    On input, A contains the matrix to be factored.
!    On output, A contains information about the factorization,
!    which must be passed unchanged to R8GE_NP_SL for solutions.
!
!    Output, integer INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,n)
  integer info
  integer j
  integer k

  info = 0

  do k = 1, n - 1

    if ( a(k,k) == 0.0D+00 ) then
      info = k
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8GE_NP_FA - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      stop 1
    end if

    a(k+1:n,k) = - a(k+1:n,k) / a(k,k)
    do j = k + 1, n
      a(k+1:n,j) = a(k+1:n,j) + a(k+1:n,k) * a(k,j)
    end do

  end do

  if ( a(n,n) == 0.0D+00 ) then
    info = n
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_NP_FA - Fatal error!'
    write ( *, '(a,i8)' ) '  Zero pivot on step ', info
    stop 1
  end if

  return
end
subroutine r8ge_np_sl ( n, a_lu, b, job )

!*****************************************************************************80
!
!! R8GE_NP_SL solves a system factored by R8GE_NP_FA.
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
!    16 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real ( kind = rk ) A_LU(N,N), the LU factors from R8GE_NP_FA.
!
!    Input/output, real ( kind = rk ) B(N).
!
!    On input, B contains the right hand side vector B.
!    On output, B contains the solution X.
!
!    Input, integer JOB.
!    If JOB is zero, the routine will solve A * x = b.
!    If JOB is nonzero, the routine will solve A' * x = b.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a_lu(n,n)
  real ( kind = rk ) b(n)
  integer job
  integer k
!
!  Solve A * x = b.
!
  if ( job == 0 ) then

    do k = 1, n - 1
      b(k+1:n) = b(k+1:n) + a_lu(k+1:n,k) * b(k)
    end do

    do k = n, 1, -1
      b(k) = b(k) / a_lu(k,k)
      b(1:k-1) = b(1:k-1) - a_lu(1:k-1,k) * b(k)
    end do
!
!  Solve A' * X = B.
!
  else

    do k = 1, n
      b(k) = ( b(k) - sum ( b(1:k-1) * a_lu(1:k-1,k) ) ) / a_lu(k,k)
    end do

    do k = n - 1, 1, -1
      b(k) = b(k) + sum ( b(k+1:n) * a_lu(k+1:n,k) )
    end do

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
subroutine r8vec_to_r8cb ( m, n, ml, mu, x, a )

!*****************************************************************************80
!
!! R8VEC_TO_R8CB copies an R8VEC into an R8CB matrix.
!
!  Discussion:
!
!    In C++ and FORTRAN, this routine is not really needed.  In MATLAB,
!    a data item carries its dimensionality implicitly, and so cannot be
!    regarded sometimes as a vector and sometimes as an array.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns 
!    in the array.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!
!    Input, real ( kind = rk ) X((ML+MU+1)*N), the vector to be copied
!    into the array.
!
!    Output, real ( kind = rk ) A(ML+MU+1,N), the array.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer ml
  integer mu
  integer n

  real ( kind = rk ) a(ml+mu+1,n)
  integer i
  integer j
  real ( kind = rk ) x((ml+mu+1)*n)

  do j = 1, n
    do i = 1, ml + mu + 1

      if ( 1 <= i + j - mu - 1 .and. i + j - mu - 1 <= m ) then
        a(i,j) = x(i+(ml+mu+1)*(j-1))
      else
        a(i,j) = 0.0D+00
      end if

    end do
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

