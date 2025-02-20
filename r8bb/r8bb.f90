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
subroutine r8bb_add ( n1, n2, ml, mu, a, i, j, value )

!*****************************************************************************80
!
!! R8BB_ADD adds a value to an entry in an R8BB matrix.
!
!  Discussion:
!
!    The R8BB storage format is for a border banded matrix.  Such a
!    matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
!    general band format.  The reason for the factor of 2 in front of
!    ML is to allocate space that may be required if pivoting occurs.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 September 1999
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
!    Input/output, real ( kind = rk ) A((2*ML+MU+1)*N1+2*N1*N2+N2*N2),
!    the R8BB matrix.
!
!    Input, integer I, J, the row and column of the entry to be
!    incremented.  Some combinations of I and J are illegal.
!
!    Input, real ( kind = rk ) VALUE, the value to be added to the 
!    (I,J)-th entry.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer n1
  integer n2

  real ( kind = rk ) a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
  integer i
  integer ij
  integer j
  real ( kind = rk ) value

  if ( value == 0.0D+00 ) then
    return
  end if

  if ( i <= 0 .or. n1 + n2 < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8BB_ADD - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal input value of row index I = ', i
    stop 1
  end if

  if ( j <= 0 .or. n1 + n2 < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8BB_ADD - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal input value of column index J = ', j
    stop 1
  end if
!
!  The A1 block of the matrix.
!
!  Check for out of band problems.
!
!  Normally, we would check the condition MU < (J-I), but the storage
!  format requires extra entries be set aside in case of pivoting, which
!  means that the condition becomes MU+ML < (J-I).
!
  if ( i <= n1 .and. j <= n1 ) then
    if ( ( mu + ml ) < ( j - i ) .or. ml < ( i - j ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8BB_ADD - Warning!'
      write ( *, '(a,i8,a,i8,a)' ) '  Unable to add to entry (', i, ',', j, ').'
    else
      ij = ( i - j + ml + mu + 1 ) + ( j - 1 ) * ( 2 * ml + mu + 1 )
    end if
!
!  The A2 block of the matrix.
!
  else if ( i <= n1 .and. n1 < j ) then
    ij = ( 2 * ml + mu + 1 ) * n1 + ( j - n1 - 1 ) * n1 + i
!
!  The A3 and A4 blocks of the matrix.
!
  else if ( n1 < i ) then
    ij = ( 2 * ml + mu + 1 ) * n1 + n2 * n1 + ( j - 1 ) * n2 + ( i - n1 )
  end if

  a(ij) = a(ij) + value

  return
end
subroutine r8bb_dif2 ( n1, n2, ml, mu, a )

!*****************************************************************************80
!
!! R8BB_DIF2 sets up an R8BB second difference matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 July 2016
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
!    Output, real ( kind = rk ) A((2*ML+MU+1)*N1+2*N1*N2+N2*N2), the R8BB matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer n1
  integer n2

  real ( kind = rk ) a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
  integer i
  integer j
  real ( kind = rk ) value

  a(1:(2*ml+mu+1)*n1+2*n1*n2+n2*n2) = 0.0D+00

  if ( ml < 1 .or. mu < 1 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'R8BB_DIF2 - Fatal error!'
    write ( *, '(a)' ) '  1 <= ML and 1 <= MU required.'
    stop 1
  end if

  do i = 2, n1 + n2
    j = i - 1
    value = - 1.0D+00
    call r8bb_set ( n1, n2, ml, mu, a, i, j, value )
  end do

  do i = 1, n1 + n2
    j = i
    value = 2.0D+00
    call r8bb_set ( n1, n2, ml, mu, a, i, j, value )
  end do

  do i = 1, n1 + n2 - 1
    j = i + 1
    value = - 1.0D+00
    call r8bb_set ( n1, n2, ml, mu, a, i, j, value )
  end do

  return
end
subroutine r8bb_fa ( n1, n2, ml, mu, a, pivot, info )

!*****************************************************************************80
!
!! R8BB_FA factors an R8BB matrix.
!
!  Discussion:
!
!    The R8BB storage format is for a border banded matrix.  Such a
!    matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
!    general band format.  The reason for the factor of 2 in front of
!    ML is to allocate space that may be required if pivoting occurs.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
!
!    Once the matrix has been factored by R8BB_FA, R8BB_SL may be called
!    to solve linear systems involving the matrix.
!
!    R8BB_FA uses LINPACK routines to carry out the factorization.
!
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
!         X1 + inv(A1) A2 X2 = inv(A1) B1
!
!      A3 X1 +         A4 X2 = B2
!
!    and then rewrites the second equation as
!
!      ( A4 - A3 inv(A1) A2 ) X2 = B2 - A3 inv(A1) B1
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 September 1999
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
!    Input/output, real ( kind = rk ) A( (2*ML+MU+1)*N1 + 2*N1*N2 + N2*N2 ).
!    On input, the border-banded matrix to be factored.
!    On output, information describing a partial factorization
!    of the original coefficient matrix.  This information is required
!    by R8BB_SL in order to solve linear systems associated with that
!    matrix.
!
!    Output, integer PIVOT(N1+N2), contains pivoting information.
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

  real ( kind = rk ) a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
  integer i
  integer ij
  integer ik
  integer info
  integer pivot(n1+n2)
  integer j
  integer jk
  integer job
  integer k
  integer nband

  nband = ( 2 * ml + mu + 1 ) * n1
!
!  Factor the A1 band matrix, overwriting A1 by its factors.
!
  if ( 0 < n1 ) then

    call r8gb_fa ( n1, ml, mu, a, pivot, info )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8BB_FA - Fatal error!'
      write ( *, '(a,i8)' ) '  R8GB_FA returned INFO = ', info
      write ( *, '(a)' ) '  Factoring failed for column INFO.'
      write ( *, '(a)' ) '  The band matrix A1 is singular.'
      write ( *, '(a)' ) '  This algorithm cannot continue!'
      stop 1
    end if

  end if

  if ( 0 < n1 .and. 0 < n2 ) then
!
!  Solve A1 * x = -A2 for x, and overwrite A2 by the results.
!
    do i = nband + 1, nband + n1 * n2
      a(i) = - a(i)
    end do

    job = 0
    do j = 1, n2
      call r8gb_sl ( n1, ml, mu, a, pivot, a(nband+(j-1)*n1+1), job )
    end do
!
!  A4 := A4 + A3 * A2.
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

    call r8ge_fa ( n2, a(nband+2*n1*n2+1), pivot(n1+1), info )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8BB_FA - Fatal error!'
      write ( *, '(a,i8)' ) '  R8GE_FA returned INFO = ',info
      write ( *, '(a)' ) '  This indicates singularity in column INFO.'
      write ( *, '(a,i8)' ) '  of the A4 submatrix, which is column ', n1+info
      write ( *, '(a)' ) '  of the full matrix.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  It is possible that the full matrix is '
      write ( *, '(a)' ) '  nonsingular, but the algorithm R8BB_FA may'
      write ( *, '(a)' ) '  not be used for this matrix.'
      stop 1
    end if
  end if

  return
end
subroutine r8bb_get ( n1, n2, ml, mu, a, i, j, value )

!*****************************************************************************80
!
!! R8BB_GET returns an entry of an R8BB matrix.
!
!  Discussion:
!
!    The R8BB storage format is for a border banded matrix.  Such a
!    matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
!    general band format.  The reason for the factor of 2 in front of
!    ML is to allocate space that may be required if pivoting occurs.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
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
!    blocks. N1 and N2 must be nonnegative, and at least one must be positive.
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input, real ( kind = rk ) A((2*ML+MU+1)*N1+2*N1*N2+N2*N2), the R8BB matrix.
!
!    Input, integer I, J, the row and column of the entry to 
!    be retrieved.
!
!    Output, real ( kind = rk ) VALUE, the value of the (I,J) entry.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer n1
  integer n2

  real ( kind = rk ) a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
  integer i
  integer ij
  integer j
  real ( kind = rk ) value

  if ( i <= 0 .or. n1 + n2 < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8BB_GET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal input value of row index I = ', i
    stop 1
  end if

  if ( j <= 0 .or. n1 + n2 < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8BB_GET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal input value of column index J = ', j
    stop 1
  end if
!
!  The A1 block of the matrix.
!
!  Check for out of band problems.
!
!  Normally, we would check the condition MU < (J-I), but the storage
!  format requires extra entries be set aside in case of pivoting, which
!  means that the condition becomes MU+ML < (J-I).
!
  if ( i <= n1 .and. j <= n1 ) then
    if ( mu + ml < ( j - i ) .or. ml < ( i - j ) ) then
      value = 0.0D+00
      return
    else
      ij = ( i - j + ml + mu + 1 ) + ( j - 1 ) * ( 2 * ml + mu + 1 )
    end if
!
!  The A2 block of the matrix.
!
  else if ( i <= n1 .and. n1 < j ) then
    ij = ( 2 * ml + mu + 1 ) * n1 + ( j - n1 - 1 ) * n1 + i
!
!  The A3 and A4 blocks of the matrix.
!
  else if ( n1 < i ) then
    ij = ( 2 * ml + mu + 1 ) * n1 + n2 * n1 + ( j - 1 ) * n2 + ( i - n1 )
  end if

  value = a(ij)

  return
end
subroutine r8bb_indicator ( n1, n2, ml, mu, a )

!*****************************************************************************80
!
!! R8BB_INDICATOR sets up an R8BB indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R8BB storage format is for a border banded matrix.  Such a
!    matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.
!
!  Example:
!
!    With N1 = 4, N2 = 1, ML = 1, MU = 2, the matrix entries would be:
!
!       00
!       00  00
!       00  00  00 --- ---
!      A11 A12 A13  00 ---  A16 A17
!      A21 A22 A23 A24  00  A26 A27
!      --- A32 A33 A34 A35  A36 A37
!      --- --- A43 A44 A45  A46 A47
!      --- --- --- A54 A55  A56 A57
!                       00
!
!      A61 A62 A63 A64 A65  A66 A67
!      A71 A72 A73 A74 A75  A76 A77
!
!    The matrix is actually stored as a vector, and we will simply suggest
!    the structure and values of the indicator matrix as:
!
!      00 00 00 00 00
!      00 00 13 24 35     16 17     61 62 63 64 65     66 67
!      00 12 23 34 45  +  26 27  +  71 72 73 74 75  +  76 77
!      11 22 33 44 55     36 37     
!      21 32 43 54 00     46 47     
!                         56 57     
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    08 January 2004
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
!    Output, real ( kind = rk ) A((2*ML+MU+1)*N1+2*N1*N2+N2*N2), the R8BB matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer n1
  integer n2

  real ( kind = rk ) a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
  integer base
  integer fac
  integer i
  integer i4_log_10
  integer j
  integer row

  a(1:(2*ml+mu+1)*n1+2*n1*n2+n2*n2) = 0.0D+00

  fac = 10 ** ( i4_log_10 ( n1 + n2 ) + 1 )
!
!  Set the banded matrix A1.
!
  do j = 1, n1
    do row = 1, 2 * ml + mu + 1
      i = row + j - ml - mu - 1
      if ( ml < row .and. 1 <= i .and. i <= n1 ) then
        a(row+(j-1)*(2*ml+mu+1)) = real ( fac * i + j, kind = rk )
      end if
    end do
  end do
!
!  Set the N1 by N2 rectangular strip A2.
!
  base = ( 2 * ml + mu + 1 ) * n1

  do i = 1, n1
    do j = n1 + 1, n1 + n2
      a(base + i + (j-n1-1)*n1 ) = real ( fac * i + j, kind = rk )
    end do
  end do
!
!  Set the N2 by N1 rectangular strip A3.
!
  base = ( 2 * ml + mu + 1 ) * n1 + n1 * n2

  do i = n1 + 1, n1 + n2
    do j = 1, n1    
      a(base + i-n1 + (j-1)*n2 ) = real ( fac * i + j, kind = rk )
    end do
  end do
!
!  Set the N2 by N2 square A4.
!
  base = ( 2 * ml + mu + 1 ) * n1 + n1 * n2 + n2 * n1

  do i = n1 + 1, n1 + n2
    do j = n1 + 1, n1 + n2
      a(base + i-n1 + (j-n1-1)*n2 ) = real ( fac * i + j, kind = rk )
    end do
  end do

  return
end
subroutine r8bb_mtv ( n1, n2, ml, mu, a, x, b )

!*****************************************************************************80
!
!! R8BB_MTV multiplies an R8VEC by an R8BB matrix.
!
!  Discussion:
!
!    The R8BB storage format is for a border banded matrix.  Such a
!    matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
!    general band format.  The reason for the factor of 2 in front of
!    ML is to allocate space that may be required if pivoting occurs.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
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
!    Input, integer N1, N2, the order of the banded and dense 
!    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative and no greater than N1-1.
!
!    Input, real ( kind = rk ) A((2*ML+MU+1)*N1 + 2*N1*N2 + N2*N2),
!    the R8BB matrix.
!
!    Input, real ( kind = rk ) X(N1+N2), the vector to multiply A.
!
!    Output, real ( kind = rk ) B(N1+N2), the product X times A.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer n1
  integer n2

  real ( kind = rk ) a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
  real ( kind = rk ) b(n1+n2)
  integer i
  integer ihi
  integer ij
  integer ilo
  integer j
  real ( kind = rk ) x(n1+n2)
!
!  Initialize B.
!
  b(1:n1+n2) = 0.0D+00
!
!  Multiply by A1.
!
  do j = 1, n1
    ilo = max ( 1, j - mu - ml )
    ihi = min ( n1, j + ml )
    ij = ( j - 1 ) * ( 2 * ml + mu + 1 ) - j + ml + mu + 1
    do i = ilo, ihi
      b(j) = b(j) + x(i) * a(ij+i)
    end do
  end do
!
!  Multiply by A2.
!
  do j = n1 + 1, n1 + n2
    ij = ( 2 * ml + mu + 1 ) * n1 + ( j - n1 - 1 ) * n1
    do i = 1, n1
      b(j) = b(j) + x(i) * a(ij+i)
    end do
  end do
!
!  Multiply by A3 and A4.
!
  do j = 1, n1 + n2
    ij = ( 2 * ml + mu + 1 ) * n1 + n1 * n2 + ( j - 1 ) * n2 - n1
    do i = n1 + 1, n1 + n2
      b(j) = b(j) + x(i) * a(ij+i)
    end do
  end do

  return
end
subroutine r8bb_mv ( n1, n2, ml, mu, a, x, b )

!*****************************************************************************80
!
!! R8BB_MV multiplies an R8BB matrix by an R8VEC.
!
!  Discussion:
!
!    The R8BB storage format is for a border banded matrix.  Such a
!    matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
!    general band format.  The reason for the factor of 2 in front of
!    ML is to allocate space that may be required if pivoting occurs.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
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
!    blocks  N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative and no greater than N1-1.
!
!    Input, real ( kind = rk ) A((2*ML+MU+1)*N1+2*N1*N2+N2*N2), the R8BB matrix.
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

  real ( kind = rk ) a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
  real ( kind = rk ) b(n1+n2)
  integer ihi
  integer ij
  integer ilo
  integer j
  real ( kind = rk ) x(n1+n2)
!
!  Initialize B.
!
  b(1:n1+n2) = 0.0D+00
!
!  Multiply by A1.
!
  do j = 1, n1

    ilo = max ( 1, j - mu - ml )
    ihi = min ( n1, j + ml )
    ij = ( j - 1 ) * ( 2 * ml + mu + 1 ) - j + ml + mu + 1

    b(ilo:ihi) = b(ilo:ihi) + a(ij+ilo:ij+ihi) * x(j)

  end do
!
!  Multiply by A2.
!
  do j = n1 + 1, n1 + n2
    ij = ( 2 * ml + mu + 1 ) * n1 + ( j - n1 - 1 ) * n1

    b(1:n1) = b(1:n1) + a(ij+1:ij+n1) * x(j)

  end do
!
!  Multiply by A3 and A4.
!
  do j = 1, n1 + n2
    ij = ( 2 * ml + mu + 1 ) * n1 + n1 * n2 + ( j - 1 ) * n2 - n1

    b(n1+1:n1+n2) = b(n1+1:n1+n2) + a(ij+n1+1:ij+n1+n2) * x(j)

  end do

  return
end
subroutine r8bb_print ( n1, n2, ml, mu, a, title )

!*****************************************************************************80
!
!! R8BB_PRINT prints an R8BB matrix.
!
!  Discussion:
!
!    The R8BB storage format is for a border banded matrix.  Such a
!    matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
!    general band format.  The reason for the factor of 2 in front of
!    ML is to allocate space that may be required if pivoting occurs.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
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
!    Input, real ( kind = rk ) A((2*ML+MU+1)*N1+2*N1*N2+N2*N2), the R8BB matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer n1
  integer n2

  real ( kind = rk ) a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
  character ( len = * ) title

  call r8bb_print_some ( n1, n2, ml, mu, a, 1, 1, n1 + n2, n1 + n2, title )

  return
end
subroutine r8bb_print_some ( n1, n2, ml, mu, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8BB_PRINT_SOME prints some of an R8BB matrix.
!
!  Discussion:
!
!    The R8BB storage format is for a border banded matrix.  Such a
!    matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
!    general band format.  The reason for the factor of 2 in front of
!    ML is to allocate space that may be required if pivoting occurs.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
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
!    Input, integer N1, N2, the order of the banded and dense 
!    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N1-1.
!
!    Input, real ( kind = rk ) A((2*ML+MU+1)*N1+2*N1*N2+N2*N2), the R8BB matrix.
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

  real ( kind = rk ) a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
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
          if ( ( j - i ) <= mu + ml .and. ( i - j ) <= ml ) then
            ij = ( i - j + ml + mu + 1 ) + ( j - 1 ) * ( 2 * ml + mu + 1 )
            aij = a(ij)
          end if
        else if ( i <= n1 .and. n1 < j ) then
          ij = ( 2 * ml + mu + 1 ) * n1 + ( j - n1 - 1 ) * n1 + i
          aij = a(ij)
        else if ( n1 < i ) then
          ij = ( 2 * ml + mu + 1 ) * n1 + n2 * n1 + ( j - 1 ) * n2 + ( i - n1 )
          aij = a(ij)
        end if

        write ( ctemp(j2), '(g14.6)' ) aij

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8bb_random ( n1, n2, ml, mu, seed, a )

!*****************************************************************************80
!
!! R8BB_RANDOM randomizes an R8BB matrix.
!
!  Discussion:
!
!    The R8BB storage format is for a border banded matrix.  Such a
!    matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
!    general band format.  The reason for the factor of 2 in front of
!    ML is to allocate space that may be required if pivoting occurs.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
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
!    Output, real ( kind = rk ) A((2*ML+MU+1)*N1+2*N1*N2+N2*N2), the R8BB matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer n1
  integer n2

  real ( kind = rk ) a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
  real ( kind = rk ) r8_uniform_01
  integer i
  integer ilo
  integer j
  integer row
  integer seed

  a(1:(2*ml+mu+1)*n1+2*n1*n2+n2*n2) = 0.0D+00
!
!  Randomize the banded matrix A1.
!
  do j = 1, n1
    do row = 1, 2 * ml + mu + 1
      i = row + j - ml - mu - 1
      if ( ml < row .and. 1 <= i .and. i <= n1 ) then
        a(row+(j-1)*(2*ml+mu+1)) = r8_uniform_01 ( seed )
      end if
    end do
  end do
!
!  Randomize the rectangular strips A2+A3+A4.
!
  ilo = ( 2 * ml + mu + 1 ) * n1 + 1

  call r8vec_uniform_01 ( n1 * n2 + n2 * n1 + n2 * n2, seed, a(ilo:) )

  return
end
subroutine r8bb_set ( n1, n2, ml, mu, a, i, j, value )

!*****************************************************************************80
!
!! R8BB_SET sets an entry of an R8BB matrix.
!
!  Discussion:
!
!    The R8BB storage format is for a border banded matrix.  Such a
!    matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
!    general band format.  The reason for the factor of 2 in front of
!    ML is to allocate space that may be required if pivoting occurs.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
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
!    Input/output, real ( kind = rk ) A((2*ML+MU+1)*N1+2*N1*N2+N2*N2),
!    the R8BB matrix.
!
!    Input, integer I, J, the row and column of the entry to 
!    be set.
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

  real ( kind = rk ) a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
  integer i
  integer ij
  integer j
  real ( kind = rk ) value

  if ( i <= 0 .or. n1 + n2 < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8BB_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal input value of row index I = ', i
    stop 1
  end if

  if ( j <= 0 .or. n1 + n2 < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8BB_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal input value of column index J = ', j
    stop 1
  end if
!
!  The A1 block of the matrix.
!
!  Check for out of band problems.
!
!  Normally, we would check the condition MU < (J-I), but the storage
!  format requires extra entries be set aside in case of pivoting, which
!  means that the condition becomes MU+ML < (J-I).
!
  if ( i <= n1 .and. j <= n1 ) then
    if ( mu + ml < ( j - i ) .or. ml < ( i - j ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8BB_SET - Warning!'
      write ( *, '(a,i8,a,i8,a)' ) '  Unable to set entry (', i, ',', j, ').'
      return
    else
      ij = ( i - j + ml + mu + 1 ) + ( j - 1 ) * ( 2 * ml + mu + 1 )
    end if
!
!  The A2 block of the matrix.
!
  else if ( i <= n1 .and. n1 < j ) then
    ij = ( 2 * ml + mu + 1 ) * n1 + ( j - n1 - 1 ) * n1 + i
!
!  The A3 and A4 blocks of the matrix.
!
  else if ( n1 < i ) then
    ij = ( 2 * ml + mu + 1 ) * n1 + n2 * n1 + ( j - 1 ) * n2 + ( i - n1 )
  end if

  a(ij) = value

  return
end
subroutine r8bb_sl ( n1, n2, ml, mu, a_lu, pivot, b )

!*****************************************************************************80
!
!! R8BB_SL solves an R8BB system factored by R8BB_FA.
!
!  Discussion:
!
!    The R8BB storage format is for a border banded matrix.  Such a
!    matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
!    general band format.  The reason for the factor of 2 in front of
!    ML is to allocate space that may be required if pivoting occurs.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!      (same formula used for A3).
!
!    The linear system A * x = b is decomposable into the block system:
!
!      ( A1 A2 ) * (X1) = (B1)
!      ( A3 A4 )   (X2)   (B2)
!
!    All the arguments except B are input quantities only, which are
!    not changed by the routine.  They should have exactly the same values
!    they had on exit from R8BB_FA.
!
!    If more than one right hand side is to be solved, with the same matrix,
!    R8BB_SL should be called repeatedly.  However, R8BB_FA only needs to be
!    called once to create the factorization.
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
!    ML and MU must be nonnegative and no greater than N1-1.
!
!    Input, real ( kind = rk ) A_LU( (2*ML+MU+1)*N1 + 2*N1*N2 + N2*N2), 
!    the LU factors from R8BB_FA.
!
!    Input, integer PIVOT(N1+N2), the pivoting information 
!    from R8BB_FA.
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

  real ( kind = rk ) a_lu((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
  real ( kind = rk ) b(n1+n2)
  integer i
  integer ij
  integer pivot(n1+n2)
  integer j
  integer job
  integer nband

  nband = ( 2 * ml + mu + 1 ) * n1
!
!  Set B1 := inverse(A1) * B1.
!
  if ( 0 < n1 ) then
    job = 0
    call r8gb_sl ( n1, ml, mu, a_lu, pivot, b, job )
  end if
!
!  Modify the right hand side of the second linear subsystem.
!  Set B2 := B2 - A3*B1.
!
  do i = 1, n2
    do j = 1, n1
      ij = nband + n1 * n2 + ( j - 1 ) * n2 + i
      b(n1+i) = b(n1+i) - a_lu(ij) * b(j)
    end do
  end do
!
!  Set B2 := inverse(A4) * B2.
!
  if ( 0 < n2 ) then
    job = 0
    call r8ge_sl ( n2, a_lu(nband+2*n1*n2+1), pivot(n1+1), b(n1+1), job )
  end if
!
!  Modify the first subsolution.
!  Set B1 := B1 + A2*B2.
!
  do i = 1, n1
    do j = 1, n2
      ij = nband + ( j - 1 ) * n1 + i
      b(i) = b(i) + a_lu(ij) * b(n1+j)
    end do
  end do

  return
end
subroutine r8bb_to_r8ge ( n1, n2, ml, mu, a, b )

!*****************************************************************************80
!
!! R8BB_TO_R8GE copies an R8BB matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R8BB storage format is for a border banded matrix.  Such a
!    matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
!    general band format.  The reason for the factor of 2 in front of
!    ML is to allocate space that may be required if pivoting occurs.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
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
!    10 July 2016
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
!    Input, real ( kind = rk ) A((2*ML+MU+1)*N1+2*N1*N2+N2*N2), the R8BB matrix.
!
!    Output, real ( kind = rk ) B(N1+N2,N1+N2), the R8GE matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer n1
  integer n2

  real ( kind = rk ) a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
  real ( kind = rk ) b(n1+n2,n1+n2)
  integer i
  integer ij
  integer j

  b(1:n1+n2,1:n1+n2) = 0.0D+00

  do i = 1, n1
    do j = 1, n1

      if ( mu + ml < ( j - i ) .or. ml < ( i - j ) ) then
        b(i,j) = 0.0D+00
      else
        ij = ( i - j + ml + mu + 1 ) + ( j - 1 ) * ( 2 * ml + mu + 1 )
        b(i,j) = a(ij)
      end if

    end do
  end do

  do i = 1, n1
    do j = n1 + 1, n1 + n2
      ij = ( 2 * ml + mu + 1 ) * n1 + ( j - n1 - 1 ) * n1 + i
      b(i,j) = a(ij)
    end do
  end do

  do i = n1 + 1, n1 + n2
    do j = 1, n1 + n2
      ij = ( 2 * ml + mu + 1 ) * n1 + n2 * n1 + ( j - 1 ) * n2 + ( i - n1 )
      b(i,j) = a(ij)
    end do
  end do

  return
end
subroutine r8bb_zeros ( n1, n2, ml, mu, a )

!*****************************************************************************80
!
!! R8BB_ZEROS zeroes an R8BB matrix.
!
!  Discussion:
!
!    The R8BB storage format is for a border banded matrix.  Such a
!    matrix has the logical form:
!
!      A1 | A2
!      ---+---
!      A3 | A4
!
!    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
!    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
!    respectively.
!
!    A should be defined as a vector.  The user must then store
!    the entries of the four blocks of the matrix into the vector A.
!    Each block is stored by columns.
!
!    A1, the banded portion of the matrix, is stored in
!    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
!    general band format.  The reason for the factor of 2 in front of
!    ML is to allocate space that may be required if pivoting occurs.
!
!    The following formulas should be used to determine how to store
!    the entry corresponding to row I and column J in the original matrix:
!
!    Entries of A1:
!
!      1 <= I <= N1, 1 <= J <= N1, (J-I) <= MU and (I-J) <= ML.
!
!      Store the I, J entry into location
!      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
!
!    Entries of A2:
!
!      1 <= I <= N1, N1+1 <= J <= N1+N2.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
!
!    Entries of A3:
!
!      N1+1 <= I <= N1+N2, 1 <= J <= N1.
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
!
!    Entries of A4:
!
!      N1+1 <= I <= N1+N2, N1+1 <= J <= N1+N2
!
!      Store the I, J entry into location
!      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
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
!    Output, real ( kind = rk ) A((2*ML+MU+1)*N1+2*N1*N2+N2*N2), the R8BB matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer n1
  integer n2

  real ( kind = rk ) a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)

  a(1:(2*ml+mu+1)*n1+2*n1*n2+n2*n2) = 0.0D+00

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
!        i1 = max ( 1, j-mu )
!        i2 = min ( n, j+ml )
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
!    N must be positive.
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
    lm = min ( ml, n-k )

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
!    N must be positive.
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

        lm = min ( ml, n-k )
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
subroutine r8ge_sl ( n, a_lu, pivot, b, job )

!*****************************************************************************80
!
!! R8GE_SL solves a system factored by R8GE_FA.
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
!    R8GE_SL is a simplified version of the LINPACK routine SGESL.
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
!    N must be positive.
!
!    Input, real ( kind = rk ) A_LU(N,N), the LU factors from R8GE_FA.
!
!    Input, integer PIVOT(N), the pivot vector from R8GE_FA.
!
!    Input/output, real ( kind = rk ) B(N).
!    On input, the right hand side vector.
!    On output, the solution vector.
!
!    Input, integer JOB, specifies the operation.
!    0, solve A * x = b.
!    nonzero, solve A' * x = b.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a_lu(n,n)
  real ( kind = rk ) b(n)
  integer pivot(n)
  integer job
  integer k
  integer l
  real ( kind = rk ) t
!
!  Solve A * x = b.
!
  if ( job == 0 ) then
!
!  Solve PL * Y = B.
!
    do k = 1, n - 1

      l = pivot(k)

      if ( l /= k ) then
        t    = b(l)
        b(l) = b(k)
        b(k) = t
      end if

      b(k+1:n) = b(k+1:n) + a_lu(k+1:n,k) * b(k)

    end do
!
!  Solve U * X = Y.
!
    do k = n, 1, -1
      b(k) = b(k) / a_lu(k,k)
      b(1:k-1) = b(1:k-1) - a_lu(1:k-1,k) * b(k)
    end do
!
!  Solve A' * X = B.
!
  else
!
!  Solve U' * Y = B.
!
    do k = 1, n
      b(k) = ( b(k) - sum ( b(1:k-1) * a_lu(1:k-1,k) ) ) / a_lu(k,k)
    end do
!
!  Solve ( PL )' * X = Y.
!
    do k = n - 1, 1, -1

      b(k) = b(k) + sum ( b(k+1:n) * a_lu(k+1:n,k) )

      l = pivot(k)

      if ( l /= k ) then
        t    = b(l)
        b(l) = b(k)
        b(k) = t
      end if

    end do

  end if

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
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
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
!    13 August 2014
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
!    Input, integer N, the number of entries in the vector.
!
!    Input/output, integer SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = rk ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer, parameter :: i4_huge = 2147483647
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
      seed = seed + i4_huge
    end if

    r(i) = real ( seed, kind = rk ) * 4.656612875D-10

  end do

  return
end

