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
function i4_uniform_ab ( a, b, seed )

!*****************************************************************************80
!
!! I4_UNIFORM_AB returns a scaled pseudorandom I4 between A and B.
!
!  Discussion:
!
!    An I4 is an integer value.
!
!    The pseudorandom number will be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer A, B, the limits of the interval.
!
!    Input/output, integer SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer I4_UNIFORM_AB, a number between A and B.
!
  implicit none

  integer a
  integer b
  integer, parameter :: i4_huge = 2147483647
  integer i4_uniform_ab
  integer k
  real r
  integer seed
  integer value

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_UNIFORM_AB - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r = real ( seed ) * 4.656612875E-10
!
!  Scale R to lie between A-0.5 and B+0.5.
!
  r = ( 1.0E+00 - r ) * ( real ( min ( a, b ) ) - 0.5E+00 ) & 
    +             r   * ( real ( max ( a, b ) ) + 0.5E+00 )
!
!  Use rounding to convert R to an integer between A and B.
!
  value = nint ( r )

  value = max ( value, min ( a, b ) )
  value = min ( value, max ( a, b ) )

  i4_uniform_ab = value

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
subroutine r8cbb_add ( n1, n2, ml, mu, a, i, j, value )

!*****************************************************************************80
!
!! R8CBB_ADD adds a value to an entry of an R8CBB matrix.
!
!  Discussion:
!
!    The R8CBB storage format is for a compressed border banded matrix.  
!    Such a matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.  
!
!    The R8CBB format is the same as the R8BB format, except that the banded
!    matrix A1 is stored in compressed band form rather than standard
!    banded form.  In other words, we do not include the extra room
!    set aside for fill in during pivoting.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (ML+MU+1)*N1 entries of A, using the obvious variant
!    of the LINPACK general band format.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+MU+1)+(J-1)*(ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, the order of the banded and dense 
!    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input/output, real ( kind = rk ) A((ML+MU+1)*N1 + 2*N1*N2 + N2*N2),
!    the R8CBB matrix.
!
!    Input, integer I, J, the indices of the entry to be 
!    incremented.
!
!    Input, real ( kind = rk ) VALUE, the value to be added to the (I,J) entry.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer n1
  integer n2

  real ( kind = rk ) a((ml+mu+1)*n1+2*n1*n2+n2*n2)
  integer i
  integer ij
  integer j
  real ( kind = rk ) value

  if ( value == 0.0D+00 ) then
    return
  end if
!
!  Check for I or J out of bounds.
!
  if ( i <= 0 .or. n1 + n2 < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8CBB_ADD - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal input value of row index I = ', i
    stop 1
  end if

  if ( j <= 0 .or. n1 + n2 < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8CBB_ADD - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal input value of column index J = ',j
    stop 1
  end if
!
!  The A1 block of the matrix.
!
!  Check for out of band problems.
!
  if ( i <= n1 .and. j <= n1 ) then
    if ( mu < ( j - i ) .or. ml < ( i - j ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8CBB_ADD - Warning!'
      write ( *, '(a,i8,a,i8,a)' ) '  Unable to add to entry (', i, ',', j, ').'
      return
    else
      ij = ( i - j + mu + 1 ) + ( j - 1 ) * ( ml + mu + 1 )
    end if
!
!  The A2 block of the matrix:
!
  else if ( i <= n1 .and. n1 < j ) then
    ij = ( ml + mu + 1 ) * n1 + ( j - n1 - 1 ) * n1 + i
!
!  The A3 and A4 blocks of the matrix.
!
  else if ( n1 < i ) then
    ij = ( ml + mu + 1 ) * n1 + n2 * n1 + ( j - 1 ) * n2 + ( i - n1 )
  end if

  a(ij) = a(ij) + value

  return
end
subroutine r8cbb_dif2 ( n1, n2, ml, mu, a )

!*****************************************************************************80
!
!! R8CBB_DIF2 sets up an R8CBB second difference matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 July 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, the order of the banded and dense 
!    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    1 <= ML, 1 <= MU.
!
!    Output, real ( kind = rk ) A((ML+MU+1)*N1+2*N1*N2+N2*N2), the R8CBB matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer n1
  integer n2

  real ( kind = rk ) a((ml+mu+1)*n1+2*n1*n2+n2*n2)
  integer i
  integer j
  real ( kind = rk ) value

  a(1:(ml+mu+1)*n1+2*n1*n2+n2*n2) = 0.0D+00

  if ( ml < 1 .or. mu < 1 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'R8CBB_DIF2 - Fatal error!'
    write ( *, '(a)' ) '  1 <= ML and 1 <= MU required.'
    stop 1
  end if

  do i = 2, n1 + n2
    j = i - 1
    value = - 1.0D+00
    call r8cbb_set ( n1, n2, ml, mu, a, i, j, value )
  end do

  do i = 1, n1 + n2
    j = i
    value = 2.0D+00
    call r8cbb_set ( n1, n2, ml, mu, a, i, j, value )
  end do

  do i = 1, n1 + n2 - 1
    j = i + 1
    value = - 1.0D+00
    call r8cbb_set ( n1, n2, ml, mu, a, i, j, value )
  end do

  return
end
subroutine r8cbb_fa ( n1, n2, ml, mu, a, info )

!*****************************************************************************80
!
!! R8CBB_FA factors an R8CBB matrix.
!
!  Discussion:
!
!    The R8CBB storage format is for a compressed border banded matrix.  
!    Such a matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.  
!
!    The R8CBB format is the same as the R8BB format, except that the banded
!    matrix A1 is stored in compressed band form rather than standard
!    banded form.  In other words, we do not include the extra room
!    set aside for fill in during pivoting.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (ML+MU+1)*N1 entries of A, using the obvious variant
!    of the LINPACK general band format.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+MU+1)+(J-1)*(ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
!
!    Once the matrix has been factored by R8CBB_FA, R8CBB_SL may be called
!    to solve linear systems involving the matrix.
!
!    R8CBB_FA uses special non-pivoting versions of LINPACK routines to
!    carry out the factorization.  The special version of the banded
!    LINPACK solver also results in a space saving, since no entries
!    need be set aside for fill in due to pivoting.
!
!    The linear system must be border banded, of the form:
!
!      ( A1 A2 ) (X1) = (B1)
!      ( A3 A4 ) (X2)   (B2)
!
!    where A1 is a (usually big) banded square matrix, A2 and A3 are
!    column and row strips which may be nonzero, and A4 is a dense
!    square matrix.
!
!    The algorithm rewrites the system as:
!
!         X1 + inverse(A1) A2 X2 = inverse(A1) B1
!
!      A3 X1 +             A4 X2 = B2
!
!    and then rewrites the second equation as
!
!      ( A4 - A3 inverse(A1) A2 ) X2 = B2 - A3 inverse(A1) B1
!
!    The algorithm will certainly fail if the matrix A1 is singular,
!    or requires pivoting.  The algorithm will also fail if the A4 matrix,
!    as modified during the process, is singular, or requires pivoting.
!    All these possibilities are in addition to the failure that will
!    if the total matrix A is singular.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, the order of the banded and dense 
!    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input/output, real ( kind = rk ) A( (ML+MU+1)*N1 + 2*N1*N2 + N2*N2).
!    On input, A contains the compact border-banded coefficient matrix.
!    On output, A contains information describing a partial factorization
!    of the original coefficient matrix.  
!
!    Output, integer INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer n1
  integer n2

  real ( kind = rk ) a((ml+mu+1)*n1+2*n1*n2+n2*n2)
  integer i
  integer ij
  integer ik
  integer info
  integer j
  integer jk
  integer job
  integer k
  integer nband

  nband = ( ml + mu + 1 ) * n1
!
!  Factor the A1 band matrix, overwriting A1 by its factors.
!
  if ( 0 < n1 ) then

    call r8cb_np_fa ( n1, ml, mu, a, info )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8CBB_FA - Fatal error!'
      write ( *, '(a,i8)' ) '  R8CB_NP_FA returned INFO = ', info
      write ( *, '(a)' ) '  Factoring failed for column INFO.'
      write ( *, '(a)' ) '  The band matrix A1 is singular.'
      write ( *, '(a)' ) '  This algorithm cannot continue!'
      stop 1
    end if

  end if

  if ( 0 < n1 .and. 0 < n2 ) then
!
!  Set A2 := -inverse(A1) * A2.
!
    a(nband+1:nband+n1*n2) = - a(nband+1:nband+n1*n2)
    job = 0

    do j = 1, n2
      call r8cb_np_sl ( n1, ml, mu, a, a(nband+(j-1)*n1+1), job )
    end do
!
!  Set A4 := A4 + A3*A2
!
    do i = 1, n2
      do j = 1, n1
        ij = nband + n1 * n2 + ( j - 1 ) * n2 + i
        do k = 1, n2
          ik = nband + 2 * n1 * n2 + ( k - 1 ) * n2 + i
          jk = nband + ( k - 1 ) * n1 + j
          a(ik) = a(ik) + a(ij) * a(jk)
        end do
      end do
    end do

  end if
!
!  Factor A4.
!
  if ( 0 < n2 ) then

    call r8ge_np_fa ( n2, a(nband+2*n1*n2+1), info )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8CBB_FA - Fatal error!'
      write ( *, '(a,i8)' ) '  R8GE_NP_FA returned INFO = ',info
      write ( *, '(a)' ) '  This indicates singularity in column INFO'
      info = n1 + info
      write ( *, '(a,i8)' ) '  of the A4 submatrix, which is column ',info
      write ( *, '(a)' ) '  of the full matrix.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  It is possible that the full matrix is '
      write ( *, '(a)' ) '  nonsingular, but the algorithm R8CBB_FA may'
      write ( *, '(a)' ) '  not be used for this matrix.'
      stop 1
    end if

  end if

  return
end
subroutine r8cbb_get ( n1, n2, ml, mu, a, i, j, value )

!*****************************************************************************80
!
!! R8CBB_GET returns the value of an entry of an R8CBB matrix.
!
!  Discussion:
!
!    The R8CBB storage format is for a compressed border banded matrix.  
!    Such a matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.  
!
!    The R8CBB format is the same as the R8BB format, except that the banded
!    matrix A1 is stored in compressed band form rather than standard
!    banded form.  In other words, we do not include the extra room
!    set aside for fill in during pivoting.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (ML+MU+1)*N1 entries of A, using the obvious variant
!    of the LINPACK general band format.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+MU+1)+(J-1)*(ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
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
!    Input, integer N1, N2, the order of the banded and dense 
!    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input, real ( kind = rk ) A((ML+MU+1)*N1 + 2*N1*N2 + N2*N2),
!    the R8CBB matrix.
!
!    Input, integer I, J, the row and column of the entry to 
!    retrieve.
!
!    Output, real ( kind = rk ) VALUE, the value of the (I,J) entry.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer n1
  integer n2

  real ( kind = rk ) a((ml+mu+1)*n1+2*n1*n2+n2*n2)
  integer i
  integer ij
  integer j
  real ( kind = rk ) value
!
!  Check for I or J out of bounds.
!
  if ( i <= 0 .or. n1 + n2 < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8CBB_GET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal input value of row index I = ', i
    stop 1
  end if

  if ( j <= 0 .or. n1 + n2 < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8CBB_GET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal input value of column index J = ', j
    stop 1
  end if
!
!  The A1 block of the matrix.
!
!  Check for out of band problems.
!
  if ( i <= n1 .and. j <= n1 ) then
    if ( mu < ( j - i ) .or. ml < ( i - j ) ) then
      value = 0.0D+00
      return
    else
      ij = ( i - j + mu + 1 ) + ( j - 1 ) * ( ml + mu + 1 )
    end if
!
!  The A2 block of the matrix:
!
  else if ( i <= n1 .and. n1 < j ) then
    ij = ( ml + mu + 1 ) * n1 + ( j - n1 - 1 ) * n1 + i
!
!  The A3 and A4 blocks of the matrix.
!
  else if ( n1 < i ) then
    ij = ( ml + mu + 1 ) * n1 + n2 * n1 + ( j - 1 ) * n2 + ( i - n1 )
  end if

  value = a(ij)

  return
end
subroutine r8cbb_indicator ( n1, n2, ml, mu, a )

!*****************************************************************************80
!
!! R8CBB_INDICATOR sets up an R8CBB indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R8CBB storage format is for a compressed border banded matrix.  
!    Such a matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.  
!
!    The R8CBB format is the same as the R8BB format, except that the banded
!    matrix A1 is stored in compressed band form rather than standard
!    banded form.  In other words, we do not include the extra room
!    set aside for fill in during pivoting.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (ML+MU+1)*N1 entries of A, using the obvious variant
!    of the LINPACK general band format.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+MU+1)+(J-1)*(ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, the order of the banded and dense 
!    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative and no greater than N1-1.
!
!    Output, real ( kind = rk ) A((ML+MU+1)*N1+2*N1*N2+N2*N2), the R8CBB matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer n1
  integer n2

  real ( kind = rk ) a((ml+mu+1)*n1+2*n1*n2+n2*n2)
  integer base
  integer fac
  integer i
  integer i4_log_10
  integer j
  integer row

  a(1:(ml+mu+1)*n1+2*n1*n2+n2*n2) = 0.0D+00

  fac = 10 ** ( i4_log_10 ( n1 + n2 ) + 1 )
!
!  Set the banded matrix A1.
!
  do j = 1, n1
    do row = 1, ml + mu + 1
      i = row + j - mu - 1
      if ( 1 <= i .and. i <= n1 ) then
        a(row+(j-1)*(ml+mu+1)) = real ( fac * i + j, kind = rk )
      end if
    end do
  end do
!
!  Set the N1 by N2 rectangular strip A2.
!
  base = ( ml + mu + 1 ) * n1

  do i = 1, n1
    do j = n1 + 1, n1 + n2
      a(base + i + (j-n1-1)*n1 ) = real ( fac * i + j, kind = rk )
    end do
  end do
!
!  Set the N2 by N1 rectangular strip A3.
!
  base = ( ml + mu + 1 ) * n1 + n1 * n2

  do i = n1 + 1, n1 + n2
    do j = 1, n1    
      a(base + i-n1 + (j-1)*n2 ) = real ( fac * i + j, kind = rk )
    end do
  end do
!
!  Set the N2 by N2 square A4.
!
  base = ( ml + mu + 1 ) * n1 + n1 * n2 + n2 * n1

  do i = n1 + 1, n1 + n2
    do j = n1 + 1, n1 + n2
      a(base + i-n1 + (j-n1-1)*n2 ) = real ( fac * i + j, kind = rk )
    end do
  end do

  return
end
subroutine r8cbb_mtv ( n1, n2, ml, mu, a, x, b )

!*****************************************************************************80
!
!! R8CBB_MTV multiplies an R8VEC by an R8CBB matrix.
!
!  Discussion:
!
!    The R8CBB storage format is for a compressed border banded matrix.  
!    Such a matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.  
!
!    The R8CBB format is the same as the R8BB format, except that the banded
!    matrix A1 is stored in compressed band form rather than standard
!    banded form.  In other words, we do not include the extra room
!    set aside for fill in during pivoting.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (ML+MU+1)*N1 entries of A, using the obvious variant
!    of the LINPACK general band format.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+MU+1)+(J-1)*(ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
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
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input, integer N1, N2, the order of the banded and dense 
!    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, real ( kind = rk ) A((ML+MU+1)*N1 + 2*N1*N2 + N2*N2),
!    the R8CBB matrix.
!
!    Input, real ( kind = rk ) X(N1+N2), the vector to multiply the matrix.
!
!    Output, real ( kind = rk ) B(N1+N2), the product X * A.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer n1
  integer n2

  real ( kind = rk ) a((ml+mu+1)*n1+2*n1*n2+n2*n2)
  real ( kind = rk ) b(n1+n2)
  integer i
  integer ihi
  integer ij
  integer ilo
  integer j
  real ( kind = rk ) x(n1+n2)
!
!  Set B to zero.
!
  b(1:n1+n2) = 0.0D+00
!
!  Multiply by A1.
!
  do j = 1, n1
    ilo = max ( 1, j - mu )
    ihi = min ( n1, j + ml )
    ij = ( j - 1 ) * ( ml + mu + 1 ) - j + mu + 1
    b(j) = b(j) + sum ( x(ilo:ihi) * a(ij+ilo:ij+ihi) )
  end do
!
!  Multiply by A2.
!
  do j = n1 + 1, n1 + n2
    ij = ( ml + mu + 1 ) * n1 + ( j - n1 - 1 ) * n1
    b(j) = b(j) + sum ( x(1:n1) * a(ij+1:ij+n1) )
  end do
!
!  Multiply by A3 and A4.
!
  do j = 1, n1 + n2
    ij = ( ml + mu + 1 ) * n1 + n1 * n2 + ( j - 1 ) * n2 - n1
    do i = n1 + 1, n1 + n2
      b(j) = b(j) + x(i) * a(ij+i)
    end do
  end do

  return
end
subroutine r8cbb_mv ( n1, n2, ml, mu, a, x, b )

!*****************************************************************************80
!
!! R8CBB_MV multiplies an R8CBB matrix by an R8VEC.
!
!  Discussion:
!
!    The R8CBB storage format is for a compressed border banded matrix.  
!    Such a matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.  
!
!    The R8CBB format is the same as the R8BB format, except that the banded
!    matrix A1 is stored in compressed band form rather than standard
!    banded form.  In other words, we do not include the extra room
!    set aside for fill in during pivoting.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (ML+MU+1)*N1 entries of A, using the obvious variant
!    of the LINPACK general band format.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+MU+1)+(J-1)*(ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
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
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input, integer N1, N2, the order of the banded and dense 
!    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, real ( kind = rk ) A((ML+MU+1)*N1 + 2*N1*N2 + N2*N2),
!    the R8CBB matrix.
!
!    Input, real ( kind = rk ) X(N1+N2), the vector to be multiplied by A.
!
!    Output, real ( kind = rk ) B(N1+N2), the result of multiplying A by X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer n1
  integer n2

  real ( kind = rk ) a((ml+mu+1)*n1+2*n1*n2+n2*n2)
  real ( kind = rk ) b(n1+n2)
  integer ihi
  integer ij
  integer ilo
  integer j
  real ( kind = rk ) x(n1+n2)
!
!  Set B to zero.
!
  b(1:n1+n2) = 0.0D+00
!
!  Multiply by A1.
!
  do j = 1, n1
    ilo = max ( 1, j - mu )
    ihi = min ( n1, j + ml )
    ij = ( j - 1 ) * ( ml + mu + 1 ) - j + mu + 1
    b(ilo:ihi) = b(ilo:ihi) + a(ij+ilo:ij+ihi) * x(j)
  end do
!
!  Multiply by A2.
!
  do j = n1 + 1, n1 + n2
    ij = ( ml + mu + 1 ) * n1 + ( j - n1 - 1 ) * n1
    b(1:n1) = b(1:n1) + a(ij+1:ij+n1) * x(j)
  end do
!
!  Multiply by A3 and A4.
!
  do j = 1, n1 + n2
    ij = ( ml + mu + 1 ) * n1 + n1 * n2 + ( j - 1 ) * n2 - n1
    b(n1+1:n1+n2) = b(n1+1:n1+n2) + a(ij+n1+1:ij+n1+n2) * x(j)
  end do

  return
end
subroutine r8cbb_print ( n1, n2, ml, mu, a, title )

!*****************************************************************************80
!
!! R8CBB_PRINT prints an R8CBB matrix.
!
!  Discussion:
!
!    The R8CBB storage format is for a compressed border banded matrix.  
!    Such a matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.  
!
!    The R8CBB format is the same as the R8BB format, except that the banded
!    matrix A1 is stored in compressed band form rather than standard
!    banded form.  In other words, we do not include the extra room
!    set aside for fill in during pivoting.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (ML+MU+1)*N1 entries of A, using the obvious variant
!    of the LINPACK general band format.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+MU+1)+(J-1)*(ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
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
!    Input, integer N1, N2, the order of the banded and dense 
!    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input, real ( kind = rk ) A((ML+MU+1)*N1+2*N1*N2+N2*N2), the R8CBB matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer n1
  integer n2

  real ( kind = rk ) a((ml+mu+1)*n1+2*n1*n2+n2*n2)
  character ( len = * ) title

  call r8cbb_print_some ( n1, n2, ml, mu, a, 1, 1, n1 + n2, n1 + n2, title )

  return
end
subroutine r8cbb_print_some ( n1, n2, ml, mu, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8CBB_PRINT_SOME prints some of an R8CBB matrix.
!
!  Discussion:
!
!    The R8CBB storage format is for a compressed border banded matrix.  
!    Such a matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.  
!
!    The R8CBB format is the same as the R8BB format, except that the banded
!    matrix A1 is stored in compressed band form rather than standard
!    banded form.  In other words, we do not include the extra room
!    set aside for fill in during pivoting.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (ML+MU+1)*N1 entries of A, using the obvious variant
!    of the LINPACK general band format.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+MU+1)+(J-1)*(ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, the order of the banded and dense 
!    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input, real ( kind = rk ) A((ML+MU+1)*N1+2*N1*N2+N2*N2), the R8CBB matrix.
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
  integer n1
  integer n2

  real ( kind = rk ) a((ml+mu+1)*n1+2*n1*n2+n2*n2)
  real ( kind = rk ) aij
  character ( len = 14 ) ctemp(incx)
  integer i
  integer i2hi
  integer i2lo
  integer ihi
  integer ij
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
    j2hi = min ( j2hi, n1 + n2 )
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
    i2hi = min ( ihi, n1 + n2 )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        aij = 0.0D+00

        if ( i <= n1 .and. j <= n1 ) then
          if ( j - i <= mu .and. i - j <= ml ) then
            ij = ( i - j + mu + 1 ) + ( j - 1 ) * ( ml + mu + 1 )
            aij = a(ij)
          end if
        else if ( i <= n1 .and. n1 < j ) then
          ij = ( ml + mu + 1 ) * n1 + ( j - n1 - 1 ) * n1 + i
          aij = a(ij)
        else if ( n1 < i ) then
          ij = ( ml + mu + 1 ) * n1 + n2 * n1 + ( j - 1 ) * n2 + ( i - n1 )
          aij = a(ij)
        end if

        write ( ctemp(j2), '(g14.6)' ) aij

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8cbb_random ( n1, n2, ml, mu, seed, a )

!*****************************************************************************80
!
!! R8CBB_RANDOM randomizes an R8CBB matrix.
!
!  Discussion:
!
!    The R8CBB storage format is for a compressed border banded matrix.  
!    Such a matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.  
!
!    The R8CBB format is the same as the R8BB format, except that the banded
!    matrix A1 is stored in compressed band form rather than standard
!    banded form.  In other words, we do not include the extra room
!    set aside for fill in during pivoting.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (ML+MU+1)*N1 entries of A, using the obvious variant
!    of the LINPACK general band format.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+MU+1)+(J-1)*(ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, the order of the banded and dense 
!    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative and no greater than N1-1.
!
!    Input/output, integer SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = rk ) A((ML+MU+1)*N1 + 2*N1*N2 + N2*N2), 
!    the R8CBB matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer n1
  integer n2

  real ( kind = rk ) a((ml+mu+1)*n1+2*n1*n2+n2*n2)
  integer i
  integer j
  integer k
  integer khi
  integer klo
  real ( kind = rk ) r8_uniform_01
  integer row
  integer seed

  a(1:(ml+mu+1)*n1+2*n1*n2+n2*n2) = 0.0D+00
!
!  Randomize the banded matrix A1.
!
  do j = 1, n1
    do row = 1, ml + mu + 1
      i = row + j - mu - 1
      if ( 1 <= i .and. i <= n1 ) then
        a(row+(j-1)*(ml+mu+1)) = r8_uniform_01 ( seed )
      end if
    end do
  end do
!
!  Randomize the rectangular strips A2+A3+A4.
!
  klo = ( ml + mu + 1 ) * n1 + 1
  khi = ( ml + mu + 1 ) * n1 + n1 * n2 + n2 * n1 + n2 * n2

  do k = klo, khi
    a(k) = r8_uniform_01 ( seed )
  end do

  return
end
subroutine r8cbb_set ( n1, n2, ml, mu, a, i, j, value )

!*****************************************************************************80
!
!! R8CBB_SET sets the value of an entry in an R8CBB matrix.
!
!  Discussion:
!
!    The R8CBB storage format is for a compressed border banded matrix.  
!    Such a matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.  
!
!    The R8CBB format is the same as the R8BB format, except that the banded
!    matrix A1 is stored in compressed band form rather than standard
!    banded form.  In other words, we do not include the extra room
!    set aside for fill in during pivoting.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (ML+MU+1)*N1 entries of A, using the obvious variant
!    of the LINPACK general band format.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+MU+1)+(J-1)*(ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
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
!    Input, integer N1, N2, the order of the banded and dense 
!    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input/output, real ( kind = rk ) A((ML+MU+1)*N1 + 2*N1*N2 + N2*N2),
!    the R8CBB matrix.
!
!    Input, integer I, J, the row and column of the entry to set.
!
!    Input, real ( kind = rk ) VALUE, the value to be assigned to the
!    (I,J) entry.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer n1
  integer n2

  real ( kind = rk ) a((ml+mu+1)*n1+2*n1*n2+n2*n2)
  integer i
  integer ij
  integer j
  real ( kind = rk ) value
!
!  Check for I or J out of bounds.
!
  if ( i <= 0 .or. n1 + n2 < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8CBB_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal input value of row index I = ', i
    stop 1
  end if

  if ( j <= 0 .or. n1 + n2 < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8CBB_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal input value of column index J = ', j
    stop 1
  end if
!
!  The A1 block of the matrix.
!
!  Check for out of band problems.
!
  if ( i <= n1 .and. j <= n1 ) then
    if ( mu < ( j - i ) .or. ml < ( i - j ) ) then
      if ( value /= 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8CBB_SET - Warning!'
        write ( *, '(a,i8,a,i8,a)' ) '  Unable to set entry (', i, ',', j, ').'
      end if
      return
    else
      ij = ( i - j + mu + 1 ) + ( j - 1 ) * ( ml + mu + 1 )
    end if
!
!  The A2 block of the matrix:
!
  else if ( i <= n1 .and. n1 < j ) then
    ij = ( ml + mu + 1 ) * n1 + ( j - n1 - 1 ) * n1 + i
!
!  The A3 and A4 blocks of the matrix.
!
  else if ( n1 < i ) then
    ij = ( ml + mu + 1 ) * n1 + n2 * n1 + ( j - 1 ) * n2 + ( i - n1 )
  end if

  a(ij) = value

  return
end
subroutine r8cbb_sl ( n1, n2, ml, mu, a_lu, b )

!*****************************************************************************80
!
!! R8CBB_SL solves an R8CBB system factored by R8CBB_FA.
!
!  Discussion:
!
!    The R8CBB storage format is for a compressed border banded matrix.  
!    Such a matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.  
!
!    The R8CBB format is the same as the R8BB format, except that the banded
!    matrix A1 is stored in compressed band form rather than standard
!    banded form.  In other words, we do not include the extra room
!    set aside for fill in during pivoting.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (ML+MU+1)*N1 entries of A, using the obvious variant
!    of the LINPACK general band format.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+MU+1)+(J-1)*(ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
!
!
!    The linear system A * x = b is decomposable into the block system:
!
!      ( A1 A2 ) * (X1) = (B1)
!      ( A3 A4 )   (X2)   (B2)
!
!    where A1 is a (usually big) banded square matrix, A2 and A3 are
!    column and row strips which may be nonzero, and A4 is a dense
!    square matrix.
!
!    All the arguments except B are input quantities only, which are
!    not changed by the routine.  They should have exactly the same values
!    they had on exit from R8CBB_FA.
!
!    If more than one right hand side is to be solved, with the same
!    matrix, R8CBB_SL should be called repeatedly.  However, R8CBB_FA only
!    needs to be called once to create the factorization.
!
!    See the documentation of R8CBB_FA for details on the matrix storage.
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
!    Input, integer N1, N2, the order of the banded and dense 
!    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input, real ( kind = rk ) A_LU( (ML+MU+1)*N1 + 2*N1*N2 + N2*N2).
!    the LU factors from R8CBB_FA.
!
!    Input/output, real ( kind = rk ) B(N1+N2).
!    On input, B contains the right hand side of the linear system.
!    On output, B contains the solution.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer n1
  integer n2

  real ( kind = rk ) a_lu((ml+mu+1)*n1+2*n1*n2+n2*n2)
  real ( kind = rk ) b(n1+n2)
  integer i
  integer ij
  integer j
  integer job
  integer nband

  nband = ( ml + mu + 1 ) * n1
!
!  Set B1 := inverse(A1) * B1.
!
  if ( 0 < n1 ) then
    job = 0
    call r8cb_np_sl ( n1, ml, mu, a_lu, b, job )
  end if
!
!  Modify the right hand side of the second linear subsystem.
!  Replace B2 by B2-A3*B1.
!
  do j = 1, n1
    do i = 1, n2
      ij = nband + n1 * n2 + ( j - 1 ) * n2 + i
      b(n1+i) = b(n1+i) - a_lu(ij) * b(j)
    end do
  end do
!
!  Solve A4*B2 = B2.
!
  if ( 0 < n2 ) then
    job = 0
    call r8ge_np_sl ( n2, a_lu(nband+2*n1*n2+1), b(n1+1), job )
  end if
!
!  Modify the first subsolution.
!  Set B1 = B1+A2*B2.
!
  do i = 1, n1
    do j = 1, n2
      ij = nband + ( j - 1 ) * n1 + i
      b(i) = b(i) + a_lu(ij) * b(n1+j)
    end do
  end do

  return
end
subroutine r8cbb_to_r8ge ( n1, n2, ml, mu, a, b )

!*****************************************************************************80
!
!! R8CBB_TO_R8GE copies an R8CBB matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R8CBB storage format is for a compressed border banded matrix.  
!    Such a matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.  
!
!    The R8CBB format is the same as the R8BB format, except that the banded
!    matrix A1 is stored in compressed band form rather than standard
!    banded form.  In other words, we do not include the extra room
!    set aside for fill in during pivoting.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (ML+MU+1)*N1 entries of A, using the obvious variant
!    of the LINPACK general band format.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+MU+1)+(J-1)*(ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
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
!    28 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, the order of the banded and dense 
!    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input, real ( kind = rk ) A((ML+MU+1)*N1+2*N1*N2+N2*N2), the R8CBB matrix.
!
!    Output, real ( kind = rk ) B(N1+N2,N1+N2), the R8GE matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer n1
  integer n2

  real ( kind = rk ) a((ml+mu+1)*n1+2*n1*n2+n2*n2)
  real ( kind = rk ) b(n1+n2,n1+n2)
  integer i
  integer ij
  integer j

  b(1:n1+n2,1:n1+n2) = 0.0D+00

  do i = 1, n1
    do j = 1, n1

      if ( i - ml <= j .and. j <= i + mu ) then
        ij = ( i - j + mu + 1 ) + ( j - 1 ) * ( ml + mu + 1 )
        b(i,j) = a(ij)
      end if

    end do
  end do

  do i = 1, n1
    do j = n1 + 1, n1 + n2
      ij = ( ml + mu + 1 ) * n1 + ( j - n1 - 1 ) * n1 + i
      b(i,j) = a(ij)
    end do
  end do

  do i = n1 + 1, n1 + n2
    do j = 1, n1 + n2
      ij = ( ml + mu + 1 ) * n1 + n2 * n1 + ( j - 1 ) * n2 + ( i - n1 )
      b(i,j) = a(ij)
    end do
  end do

  return
end
subroutine r8cbb_zeros ( n1, n2, ml, mu, a )

!*****************************************************************************80
!
!! R8CBB_ZEROS zeroes an R8CBB matrix.
!
!  Discussion:
!
!    The R8CBB storage format is for a compressed border banded matrix.  
!    Such a matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.  
!
!    The R8CBB format is the same as the R8BB format, except that the banded
!    matrix A1 is stored in compressed band form rather than standard
!    banded form.  In other words, we do not include the extra room
!    set aside for fill in during pivoting.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (ML+MU+1)*N1 entries of A, using the obvious variant
!    of the LINPACK general band format.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+MU+1)+(J-1)*(ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
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
!    Input, integer N1, N2, the order of the banded and dense 
!    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative and no greater than N1-1.
!
!    Output, real ( kind = rk ) A((ML+MU+1)*N1 + 2*N1*N2 + N2*N2), 
!    the R8CBB matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer n1
  integer n2

  real ( kind = rk ) a((ml+mu+1)*n1+2*n1*n2+n2*n2)

  a(1:(ml+mu+1)*n1+2*n1*n2+n2*n2) = 0.0D+00

  return
end
subroutine r8gb_fa ( n, ml, mu, a, pivot, info )

!*****************************************************************************80
!
!! R8GB_FA performs a LINPACK-style PLU factorization of an R8GB matrix.
!
!  Discussion:
!
!    The R8GB storage format is for an M by N banded matrix, with lower 
!    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
!    extra superdiagonals, which may be required to store nonzero entries 
!    generated during Gaussian elimination.
!
!    The original M by N matrix is "collapsed" downward, so that diagonals
!    become rows of the storage array, while columns are preserved.  The
!    collapsed array is logically 2*ML+MU+1 by N.  
!
!    This routine is based on the LINPACK routine SGBFA.
!
!    R8GB storage is used by LINPACK and LAPACK.
!
!    The following program segment will set up the input.
!
!      m = ml + mu + 1
!      do j = 1, n
!        i1 = max ( 1, j - mu )
!        i2 = min ( n, j + ml )
!        do i = i1, i2
!          k = i - j + m
!          a(k,j) = afull(i,j)
!        end do
!      end do
!
!    This uses rows ML+1 through 2*ML+MU+1 of the array A.
!    In addition, the first ML rows in the array are used for
!    elements generated during the triangularization.
!
!    The ML+MU by ML+MU upper left triangle and the
!    ML by ML lower right triangle are not referenced.
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
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input/output, real ( kind = rk ) A(2*ML+MU+1,N), on input, 
!    the matrix in band storage, on output, information about 
!    the LU factorization.
!
!    Output, integer PIVOT(N), the pivot vector.
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

  real ( kind = rk ) a(2*ml+mu+1,n)
  integer i0
  integer info
  integer pivot(n)
  integer j
  integer j0
  integer j1
  integer ju
  integer jz
  integer k
  integer l
  integer lm
  integer m
  integer mm
  real ( kind = rk ) t

  m = ml + mu + 1
  info = 0
!
!  Zero out the initial fill-in columns.
!
  j0 = mu + 2
  j1 = min ( n, m ) - 1

  do jz = j0, j1
    i0 = m + 1 - jz
    a(i0:ml,jz) = 0.0D+00
  end do

  jz = j1
  ju = 0

  do k = 1, n - 1
!
!  Zero out the next fill-in column.
!
    jz = jz + 1
    if ( jz <= n ) then
      a(1:ml,jz) = 0.0D+00
    end if
!
!  Find L = pivot index.
!
    lm = min ( ml, n - k )

    l = m
    do j = m + 1, m + lm
      if ( abs ( a(l,k) ) < abs ( a(j,k) ) ) then
        l = j
      end if
    end do

    pivot(k) = l + k - m
!
!  Zero pivot implies this column already triangularized.
!
    if ( a(l,k) == 0.0D+00 ) then
      info = k
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8GB_FA - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      stop 1
    end if
!
!  Interchange if necessary.
!
    t      = a(l,k)
    a(l,k) = a(m,k)
    a(m,k) = t
!
!  Compute multipliers.
!
    a(m+1:m+lm,k) = - a(m+1:m+lm,k) / a(m,k)
!
!  Row elimination with column indexing.
!
    ju = max ( ju, mu + pivot(k) )
    ju = min ( ju, n )
    mm = m

    do j = k + 1, ju

      l = l - 1
      mm = mm - 1

      if ( l /= mm ) then
        t       = a(l,j)
        a(l,j)  = a(mm,j)
        a(mm,j) = t
      end if

      a(mm+1:mm+lm,j) = a(mm+1:mm+lm,j) + a(mm,j) * a(m+1:m+lm,k)

    end do

  end do

  pivot(n) = n

  if ( a(m,n) == 0.0D+00 ) then
    info = n
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GB_FA - Fatal error!'
    write ( *, '(a,i8)' ) '  Zero pivot on step ', info
    stop 1
  end if

  return
end
subroutine r8gb_sl ( n, ml, mu, a_lu, pivot, b, job )

!*****************************************************************************80
!
!! R8GB_SL solves a system factored by R8GB_FA.
!
!  Discussion:
!
!    The R8GB storage format is for an M by N banded matrix, with lower 
!    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
!    extra superdiagonals, which may be required to store nonzero entries 
!    generated during Gaussian elimination.
!
!    The original M by N matrix is "collapsed" downward, so that diagonals
!    become rows of the storage array, while columns are preserved.  The
!    collapsed array is logically 2*ML+MU+1 by N.  
!
!    R8GB storage is used by LINPACK and LAPACK.
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
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input, real ( kind = rk ) A_LU(2*ML+MU+1,N), the LU factors from R8GB_FA.
!
!    Input, integer PIVOT(N), the pivot vector from R8GB_FA.
!
!    Input/output, real ( kind = rk ) B(N).
!    On input, the right hand side vector.
!    On output, the solution.
!
!    Input, integer JOB.
!    0, solve A * x = b.
!    nonzero, solve A' * x = b.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer n

  real ( kind = rk ) a_lu(2*ml+mu+1,n)
  real ( kind = rk ) b(n)
  integer pivot(n)
  integer job
  integer k
  integer l
  integer la
  integer lb
  integer lm
  integer m
  real ( kind = rk ) t

  m = mu + ml + 1
!
!  Solve A * x = b.
!
  if ( job == 0 ) then
!
!  Solve L * Y = B.
!
    if ( 1 <= ml ) then

      do k = 1, n - 1

        lm = min ( ml, n - k )
        l = pivot(k)

        if ( l /= k ) then
          t    = b(l)
          b(l) = b(k)
          b(k) = t
        end if

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
!  Solve L' * X = Y.
!
    if ( 1 <= ml ) then

      do k = n - 1, 1, -1

        lm = min ( ml, n - k )
        b(k) = b(k) + sum ( a_lu(m+1:m+lm,k) * b(k+1:k+lm) )
        l = pivot(k)

        if ( l /= k ) then
          t    = b(l)
          b(l) = b(k)
          b(k) = t
        end if

      end do

    end if

  end if

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

