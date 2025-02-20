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
subroutine r8gb_det ( n, ml, mu, a_lu, pivot, det )

!*****************************************************************************80
!
!! R8GB_DET: determinant of a matrix factored by R8GB_FA or R8GB_TRF.
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
!    28 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
!    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
!    Sven Hammarling, Alan McKenney, Danny Sorensen,
!    LAPACK User's Guide,
!    Second Edition,
!    SIAM, 1995.
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
!    Input, real ( kind = rk ) A_LU(2*ML+MU+1,N), the LU factors from 
!    R8GB_FA or R8GB_TRF.
!
!    Input, integer PIVOT(N), the pivot vector, as computed 
!    by R8GB_FA or R8GB_TRF.
!
!    Output, real ( kind = rk ) DET, the determinant of the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer n

  real ( kind = rk ) a_lu(2*ml+mu+1,n)
  real ( kind = rk ) det
  integer i
  integer pivot(n)

  det = product ( a_lu(ml+mu+1,1:n) )

  do i = 1, n
    if ( pivot(i) /= i ) then
      det = -det
    end if
  end do

  return
end
subroutine r8gb_dif2 ( m, n, ml, mu, a )

!*****************************************************************************80
!
!! R8GB_DIF2 sets up an R8GB second difference matrix.
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
!    18 June 2016
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
!    Output, real ( kind = rk ) A(2*ML+MU+1,N), the R8GB matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer m
  integer n

  real ( kind = rk ) a(2*ml+mu+1,n)
  integer diag
  integer i
  integer j

  a(1:2*ml+mu+1,1:n) = 0.0D+00

  do j = 1, n

    do diag = 1, 2 * ml + mu + 1

      i = diag + j - ml - mu - 1

      if ( i == j ) then
        a(diag,j) = 2.0D+00
      else if ( i == j + 1 .or. i == j - 1 ) then
        a(diag,j) = -1.0D+00
      end if

    end do
  end do

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
subroutine r8gb_indicator ( m, n, ml, mu, a )

!*****************************************************************************80
!
!! R8GB_INDICATOR sets up an R8GB indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    Note that the R8GB storage format includes extra room for
!    fillin entries that occur during Gauss elimination.  These entries
!    are not normally seen or used by the user.  This routine will
!    set those values to zero.
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
!    16 March 2005
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
!    Output, real ( kind = rk ) A(2*ML+MU+1,N), the R8GB matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer m
  integer n

  real ( kind = rk ) a(2*ml+mu+1,n)
  integer diag
  integer fac
  integer i
  integer i4_log_10
  integer j
  integer k
  real ( kind = rk ) value

  fac = 10 ** ( i4_log_10 ( n ) + 1 )
  k = 0

  do j = 1, n
    do diag = 1, 2 * ml + mu + 1

      i = diag + j - ml - mu - 1

      if ( 1 <= i .and. i <= m .and. i - ml <= j .and. j <= i + mu ) then
        value = real ( fac * i + j, kind = rk )
      else if ( 1 <= i .and. i <= m .and. &
        i - ml <= j .and. j <= i + mu + ml ) then
        value = 0.0D+00
      else
        k = k + 1
        value = - real ( k, kind = rk )
      end if

      a(diag,j) = value

    end do
  end do

  return
end
subroutine r8gb_ml ( n, ml, mu, a_lu, pivot, x, b, job )

!*****************************************************************************80
!
!! R8GB_ML computes A * x or A' * X, using R8GB_FA factors.
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
!    It is assumed that R8GB_FA has overwritten the original matrix
!    information by LU factors.  R8GB_ML is able to reconstruct the
!    original matrix from the LU factor data.
!
!    R8GB_ML allows the user to check that the solution of a linear
!    system is correct, without having to save an unfactored copy
!    of the matrix.
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
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input, real ( kind = rk ) A_LU(2*ML+MU+1,N), the LU factors from R8GB_FA.
!
!    Input, integer PIVOT(N), the pivot vector computed by R8GB_FA.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied.
!
!    Output, real ( kind = rk ) B(N), the result of the multiplication.
!
!    Input, integer JOB, specifies the operation to be done:
!    JOB = 0, compute A * x.
!    JOB nonzero, compute A' * X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer n

  real ( kind = rk ) a_lu(2*ml+mu+1,n)
  real ( kind = rk ) b(n)
  integer i
  integer ihi
  integer ilo
  integer pivot(n)
  integer j
  integer jhi
  integer job
  integer k
  real ( kind = rk ) t
  real ( kind = rk ) x(n)

  b(1:n) = x(1:n)

  if ( job == 0 ) then
!
!  Y = U * X.
!
    do j = 1, n
      ilo = max ( 1, j - ml - mu )
      do i = ilo, j - 1
        b(i) = b(i) + a_lu(i-j+ml+mu+1,j) * b(j)
      end do
      b(j) = a_lu(j-j+ml+mu+1,j) * b(j)
    end do
!
!  B = PL * Y = PL * U * X = A * x.
!
    do j = n - 1, 1, -1

      ihi = min ( n, j + ml )
      do i = j + 1, ihi
        b(i) = b(i) - a_lu(i-j+ml+mu+1,j) * b(j)
      end do

      k = pivot(j)

      if ( k /= j ) then
        t    = b(k)
        b(k) = b(j)
        b(j) = t
      end if

    end do

  else
!
!  Y = ( PL )' * X.
!
    do j = 1, n - 1

      k = pivot(j)

      if ( k /= j ) then
        t    = b(k)
        b(k) = b(j)
        b(j) = t
      end if

      jhi = min ( n, j + ml )
      do i = j + 1, jhi
        b(j) = b(j) - b(i) * a_lu(i-j+ml+mu+1,j)
      end do

    end do
!
!  B = U' * Y = ( PL * U )' * X = A' * X.
!
    do i = n, 1, -1

      jhi = min ( n, i + ml + mu )
      do j = i + 1, jhi
        b(j) = b(j) + b(i) * a_lu(i-j+ml+mu+1,j)
      end do
      b(i) = b(i) * a_lu(i-i+ml+mu+1,i)
    end do

  end if

  return
end
subroutine r8gb_mtv ( m, n, ml, mu, a, x, b )

!*****************************************************************************80
!
!! R8GB_MTV multiplies an R8VEC by an R8GB matrix.
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
!    LINPACK and LAPACK storage of general band matrices requires
!    an extra ML upper diagonals for possible fill in entries during
!    Gauss elimination.  This routine does not access any entries
!    in the fill in diagonals, because it assumes that the matrix
!    has NOT had Gauss elimination applied to it.  If the matrix
!    has been Gauss eliminated, then the routine R8GB_MU must be
!    used instead.
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
!    Input, integer M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than min(M,N)-1.
!
!    Input, real ( kind = rk ) A(2*ML+MU+1,N), the R8GB matrix.
!
!    Input, real ( kind = rk ) X(M), the vector to be multiplied by A.
!
!    Output, real ( kind = rk ) B(N), the product X*A.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer m
  integer n

  real ( kind = rk ) a(2*ml+mu+1,n)
  real ( kind = rk ) b(n)
  integer i
  integer ihi
  integer ilo
  integer j
  real ( kind = rk ) x(m)

  do j = 1, n
    b(j) = 0.0D+00
    ilo = max ( 1, j - mu )
    ihi = min ( m, j + ml )
    do i = ilo, ihi
      b(j) = b(j) + x(i) * a(i-j+ml+mu+1,j)
    end do
  end do

  return
end
subroutine r8gb_mu ( n, ml, mu, a_lu, pivot, x, b, job )

!*****************************************************************************80
!
!! R8GB_MU computes A * x or A' * X, using R8GB_TRF factors.
!
!  Warning:
!
!    This routine needs to be updated to allow for rectangular matrices.
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
!    It is assumed that R8GB_TRF has overwritten the original matrix
!    information by LU factors.  R8GB_MU is able to reconstruct the
!    original matrix from the LU factor data.
!
!    R8GB_MU allows the user to check that the solution of a linear
!    system is correct, without having to save an unfactored copy
!    of the matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
!    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
!    Sven Hammarling, Alan McKenney, Danny Sorensen,
!    LAPACK User's Guide,
!    Second Edition,
!    SIAM, 1995.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input, real ( kind = rk ) A_LU(2*ML+MU+1,N), the LU factors from R8GB_TRF.
!
!    Input, integer PIVOT(N), the pivot vector computed 
!    by R8GB_TRF.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied.
!
!    Output, real ( kind = rk ) B(N), the result of the multiplication.
!
!    Input, integer JOB, specifies the operation to be done:
!    JOB = 0, compute A * x.
!    JOB nonzero, compute A' * X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer n

  real ( kind = rk ) a_lu(2*ml+mu+1,n)
  real ( kind = rk ) b(n)
  integer i
  integer ihi
  integer ilo
  integer pivot(n)
  integer j
  integer jhi
  integer job
  integer k
  real ( kind = rk ) t
  real ( kind = rk ) x(n)

  b(1:n) = x(1:n)

  if ( job == 0 ) then
!
!  Y = U * X.
!
    do j = 1, n
      ilo = max ( 1, j - ml - mu )
      do i = ilo, j - 1
        b(i) = b(i) + a_lu(i-j+ml+mu+1,j) * b(j)
      end do
      b(j) = a_lu(j-j+ml+mu+1,j) * b(j)
    end do
!
!  B = PL * Y = PL * U * X = A * x.
!
    do j = n - 1, 1, -1

      ihi = min ( n, j + ml )
      do i = j + 1, ihi
        b(i) = b(i) + a_lu(i-j+ml+mu+1,j) * b(j)
      end do

      k = pivot(j)

      if ( k /= j ) then
        t    = b(k)
        b(k) = b(j)
        b(j) = t
      end if

    end do

  else
!
!  Y = ( PL )' * X.
!
    do j = 1, n - 1

      k = pivot(j)

      if ( k /= j ) then
        t    = b(k)
        b(k) = b(j)
        b(j) = t
      end if

      jhi = min ( n, j + ml )
      do i = j + 1, jhi
        b(j) = b(j) + b(i) * a_lu(i-j+ml+mu+1,j)
      end do

    end do
!
!  B = U' * Y = ( PL * U )' * X = A' * X.
!
    do i = n, 1, -1

      jhi = min ( n, i + ml + mu )
      do j = i + 1, jhi
        b(j) = b(j) + b(i) * a_lu(i-j+ml+mu+1,j)
      end do
      b(i) = b(i) * a_lu(i-i+ml+mu+1,i)
    end do

  end if

  return
end
subroutine r8gb_mv ( m, n, ml, mu, a, x, b )

!*****************************************************************************80
!
!! R8GB_MV multiplies an R8GB matrix by an R8VEC.
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
!    LINPACK and LAPACK storage of general band matrices requires
!    an extra ML upper diagonals for possible fill in entries during
!    Gauss elimination.  This routine does not access any entries
!    in the fill in diagonals, because it assumes that the matrix
!    has NOT had Gauss elimination applied to it.  If the matrix
!    has been Gauss eliminated, then the routine R8GB_MU must be
!    used instead.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 January 1999
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
!    Input, integer M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than min(M,N)-1.
!
!    Input, real ( kind = rk ) A(2*ML+MU+1,N), the R8GB matrix.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = rk ) B(M), the product A * x.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer m
  integer n

  real ( kind = rk ) a(2*ml+mu+1,n)
  real ( kind = rk ) b(m)
  integer i
  integer j
  integer jhi
  integer jlo
  real ( kind = rk ) x(n)

  do i = 1, m
    b(i) = 0.0D+00
    jlo = max ( 1, i - ml )
    jhi = min ( n, i + mu )
    do j = jlo, jhi
      b(i) = b(i) + a(i-j+ml+mu+1,j) * x(j)
    end do
  end do

  return
end
subroutine r8gb_nz_num ( m, n, ml, mu, a, nz_num )

!*****************************************************************************80
!
!! R8GB_NZ_NUM counts the nonzeroes in an R8GB matrix.
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
!    LINPACK and LAPACK band storage requires that an extra ML
!    superdiagonals be supplied to allow for fillin during Gauss
!    elimination.  Even though a band matrix is described as
!    having an upper bandwidth of MU, it effectively has an
!    upper bandwidth of MU+ML.  This routine will examine
!    values it finds in these extra bands, so that both unfactored
!    and factored matrices can be handled.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 October 2003
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
!    Input, real ( kind = rk ) A(2*ML+MU+1,N), the R8GB matrix.
!
!    Output, integer NZ_NUM, the number of nonzero entries in A.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer m
  integer n

  real ( kind = rk ) a(2*ml+mu+1,n)
  integer i
  integer j
  integer jhi
  integer jlo
  integer nz_num

  nz_num = 0

  do i = 1, m
    jlo = max ( 1, i - ml )
    jhi = min ( n, i + mu + ml )
    do j = jlo, jhi
      if ( a(i-j+ml+mu+1,j) /= 0.0D+00 ) then
        nz_num = nz_num + 1
      end if
    end do
  end do

  return
end
subroutine r8gb_print ( m, n, ml, mu, a, title )

!*****************************************************************************80
!
!! R8GB_PRINT prints an R8GB matrix.
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
!    04 April 2006
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
!    ML and MU must be nonnegative, and no greater than min(M,N)-1..
!
!    Input, real ( kind = rk ) A(2*ML+MU+1,N), the R8GB matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer n

  real ( kind = rk ) a(2*ml+mu+1,n)
  integer m
  character ( len = * ) title

  call r8gb_print_some ( m, n, ml, mu, a, 1, 1, m, n, title )

  return
end
subroutine r8gb_print_some ( m, n, ml, mu, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8GB_PRINT_SOME prints some of an R8GB matrix.
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
!    21 June 2016
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
!    ML and MU must be nonnegative, and no greater than min(M,N)-1..
!
!    Input, real ( kind = rk ) A(2*ML+MU+1,N), the R8GB matrix.
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

  real ( kind = rk ) a(2*ml+mu+1,n)
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

        if ( mu < j - i .or. ml < i - j ) then
          ctemp(j2) = '              '
        else
          write ( ctemp(j2), '(g14.6)' ) a(i-j+ml+mu+1,j)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8gb_random ( m, n, ml, mu, seed, a )

!*****************************************************************************80
!
!! R8GB_RANDOM randomizes an R8GB matrix.
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
!    LINPACK and LAPACK band storage requires that an extra ML
!    superdiagonals be supplied to allow for fillin during Gauss
!    elimination.  Even though a band matrix is described as
!    having an upper bandwidth of MU, it effectively has an
!    upper bandwidth of MU+ML.  This routine assumes it is setting
!    up an unfactored matrix, so it only uses the first MU upper bands,
!    and does not place nonzero values in the fillin bands.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 October 2004
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
!    Input/output, integer SEED, a seed for the random number 
!    generator.
!
!    Output, real ( kind = rk ) A(2*ML+MU+1,N), the R8GB matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer n

  real ( kind = rk ) a(2*ml+mu+1,n)
  real ( kind = rk ) r8_uniform_01
  integer i
  integer j
  integer m
  integer row
  integer seed

  do j = 1, n
    do row = 1, 2 * ml + mu + 1
      i = row + j - ml - mu - 1
      if ( ml < row .and. 1 <= i .and. i <= m ) then
        a(row,j) = r8_uniform_01 ( seed )
      else
        a(row,j) = 0.0D+00
      end if
    end do
  end do

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
subroutine r8gb_to_r8ge ( m, n, ml, mu, a, b )

!*****************************************************************************80
!
!! R8GB_TO_R8GE copies an R8GB matrix to an R8GE matrix.
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
!    LINPACK and LAPACK band storage requires that an extra ML
!    superdiagonals be supplied to allow for fillin during Gauss
!    elimination.  Even though a band matrix is described as
!    having an upper bandwidth of MU, it effectively has an
!    upper bandwidth of MU+ML.  This routine will copy nonzero
!    values it finds in these extra bands, so that both unfactored
!    and factored matrices can be handled.
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
!    17 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows of the matrices.
!    M must be positive.
!
!    Input, integer N, the number of columns of the matrices.
!    N must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths of A1.
!    ML and MU must be nonnegative, and no greater than min(M,N)-1.
!
!    Input, real ( kind = rk ) A(2*ML+MU+1,N), the R8GB matrix.
!
!    Output, real ( kind = rk ) B(M,N), the R8GE matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer ml
  integer mu
  integer n

  real ( kind = rk ) a(2*ml+mu+1,n)
  real ( kind = rk ) b(m,n)
  integer i
  integer j

  do i = 1, m
    do j = 1, n
      if ( i - ml <= j .and. j <= i + mu + ml ) then
        b(i,j) = a(ml+mu+1+i-j,j)
      else
        b(i,j) = 0.0D+00
      end if
    end do
  end do

  return
end
subroutine r8gb_to_r8s3 ( m, n, ml, mu, a, nz_num, sym, row, col, b )

!*****************************************************************************80
!
!! R8GB_TO_R8S3 copies an R8GB matrix to an R8S3 matrix.
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
!    LINPACK and LAPACK band storage requires that an extra ML
!    superdiagonals be supplied to allow for fillin during Gauss
!    elimination.  Even though a band matrix is described as
!    having an upper bandwidth of MU, it effectively has an
!    upper bandwidth of MU+ML.  This routine will copy nonzero
!    values it finds in these extra bands, so that both unfactored
!    and factored matrices can be handled.
!
!    The R8S3 storage format corresponds to the SLAP Triad format.
!
!    The R8S3 storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.  The entries may be given in any order.  No
!    check is made for the erroneous case in which a given matrix entry is
!    specified more than once.
!
!    There is a symmetry option for square matrices.  If the symmetric storage
!    option is used, the format specifies that only nonzeroes on the diagonal
!    and lower triangle are stored.  However, this routine makes no attempt
!    to enforce this.  The only thing it does is to "reflect" any nonzero
!    offdiagonal value.  Moreover, no check is made for the erroneous case
!    in which both A(I,J) and A(J,I) are specified, but with different values.
!
!    This routine reorders the entries of A so that the first N entries
!    are exactly the diagonal entries of the matrix, in order.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows of the matrices.
!    M must be positive.
!
!    Input, integer N, the number of columns of the matrices.
!    N must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths of A1.
!    ML and MU must be nonnegative, and no greater than min(M,N)-1.
!
!    Input, real ( kind = rk ) A(2*ML+MU+1,N), the R8GB matrix.
!
!    Input, integer NZ_NUM, the number of nonzero entries in A.
!    This number can be obtained by calling R8GB_NZ_NUM.
!
!    Output, integer SYM, is 0 if the matrix is not symmetric,
!    and 1 if the matrix is symmetric.  If the matrix is symmetric, then
!    only the nonzeroes on the diagonal and in the lower triangle are stored.
!    For this routine, SYM is always output 0.
!
!    Output, integer ROW(NZ_NUM), the row indices.
!
!    Output, integer COL(NZ_NUM), the column indices.
!
!    Output, real ( kind = rk ) B(NZ_NUM), the R8S3 matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer ml
  integer mu
  integer n
  integer nz_num

  real ( kind = rk ) a(2*ml+mu+1,n)
  real ( kind = rk ) b(nz_num)
  integer col(nz_num)
  integer i
  integer sym
  integer j
  integer nz
  integer row(nz_num)

  sym = 0
  nz = 0

  do i = 1, m
    do j = 1, n
      if ( i - ml <= j .and. j <= i + mu + ml ) then
        if ( a(ml+mu+1+i-j,j) /= 0.0D+00 ) then

          if ( nz_num <= nz ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'R8GB_TO_R8S3 - Fatal error!'
            write ( *, '(a,i8)' ) '  NZ_NUM = ', nz_num
            write ( *, '(a)' ) '  But the matrix has more nonzeros than that!'
            stop 1
          end if

          nz = nz + 1
          row(nz) = i
          col(nz) = j
          b(nz) = a(ml+mu+1+i-j,j)

        end if
      end if
    end do
  end do

  if ( nz < nz_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GB_TO_R8S3 - Warning!'
    write ( *, '(a,i8)' ) '  NZ_NUM = ', nz_num
    write ( *, '(a,i8)' ) '  But the number of nonzeros is ', nz
  end if

  return
end
subroutine r8gb_to_r8sp ( m, n, ml, mu, a, nz_num, row, col, b )

!*****************************************************************************80
!
!! R8GB_TO_R8SP copies an R8GB matrix to an R8SP matrix.
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
!    LINPACK and LAPACK band storage requires that an extra ML
!    superdiagonals be supplied to allow for fillin during Gauss
!    elimination.  Even though a band matrix is described as
!    having an upper bandwidth of MU, it effectively has an
!    upper bandwidth of MU+ML.  This routine will copy nonzero
!    values it finds in these extra bands, so that both unfactored
!    and factored matrices can be handled.
!
!    The R8SP storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.
!
!    It is possible that a pair of indices (I,J) may occur more than
!    once.  Presumably, in this case, the intent is that the actual value
!    of A(I,J) is the sum of all such entries.  This is not a good thing
!    to do, but I seem to have come across this in MATLAB.
!
!    The R8SP format is used by CSPARSE ("sparse triplet"), SLAP 
!    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows of the matrices.
!    M must be positive.
!
!    Input, integer N, the number of columns of the matrices.
!    N must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths of A1.
!    ML and MU must be nonnegative, and no greater than min(M,N)-1.
!
!    Input, real ( kind = rk ) A(2*ML+MU+1,N), the R8GB matrix.
!
!    Input, integer NZ_NUM, the number of nonzero entries in A.
!    This number can be obtained by calling R8GB_NZ_NUM.
!
!    Output, integer ROW(NZ_NUM), the row indices.
!
!    Output, integer COL(NZ_NUM), the column indices.
!
!    Output, real ( kind = rk ) B(NZ_NUM), the R8SP matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer ml
  integer mu
  integer n
  integer nz_num

  real ( kind = rk ) a(2*ml+mu+1,n)
  real ( kind = rk ) b(nz_num)
  integer col(nz_num)
  integer i
  integer j
  integer jhi
  integer jlo
  integer nz
  integer row(nz_num)

  nz = 0

  do i = 1, m

    jlo = max ( 1, i - ml )
    jhi = min ( n, i + mu + ml )

    do j = jlo, jhi

      if ( a(ml+mu+1+i-j,j) == 0.0D+00 ) then
        cycle
      end if

      if ( nz_num <= nz ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8GB_TO_R8SP - Fatal error!'
        write ( *, '(a,i8)' ) '  NZ_NUM = ', nz_num
        write ( *, '(a)' ) '  But the matrix has more nonzeros than that!'
        stop 1
      end if

      nz = nz + 1
      row(nz) = i
      col(nz) = j
      b(nz) = a(ml+mu+1+i-j,j)

    end do
  end do

  if ( nz < nz_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GB_TO_R8SP - Warning!'
    write ( *, '(a,i8)' ) '  NZ_NUM = ', nz_num
    write ( *, '(a,i8)' ) '  But the number of nonzeros is ', nz
  end if

  return
end
subroutine r8gb_to_r8vec ( m, n, ml, mu, a, x )

!*****************************************************************************80
!
!! R8GB_TO_R8VEC copies an R8GB matrix to an R8VEC.
!
!  Discussion:
!
!    In C++ and FORTRAN, this routine is not really needed.  In MATLAB,
!    a data item carries its dimensionality implicitly, and so cannot be
!    regarded sometimes as a vector and sometimes as an array.
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
!    17 March 2004
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
!    Input, real ( kind = rk ) A(2*ML+MU+1,N), the array to be copied.
!
!    Output, real ( kind = rk ) X((2*ML+MU+1)*N), the vector.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer ml
  integer mu
  integer n

  real ( kind = rk ) a(2*ml+mu+1,n)
  integer i
  integer ihi
  integer ilo
  integer j
  real ( kind = rk ) x((2*ml+mu+1)*n)

  do j = 1, n

    ihi = min ( ml + mu, ml + mu + 1 - j )
    do i = 1, ihi
      x(i+(j-1)*(2*ml+mu+1)) = 0.0D+00
    end do

    ilo = max ( ihi + 1, 1 )
    ihi = min ( 2 * ml + mu + 1, ml + mu + m + 1 - j )
    do i = ilo, ihi
      x(i+(j-1)*(2*ml+mu+1)) = a(i,j)
    end do

    ilo = ihi + 1
    ihi = 2 * ml + mu + 1
    do i = ilo, ihi
      x(i+(j-1)*(2*ml+mu+1)) = 0.0D+00
    end do

  end do

  return
end
subroutine r8gb_trf ( m, n, ml, mu, a, pivot, info )

!*****************************************************************************80
!
!! R8GB_TRF performs a LAPACK-style PLU factorization of an R8GB matrix.
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
!    This is a simplified, standalone version of the LAPACK
!    routine R8GBTRF.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 January 1999
!
!  Author:
!
!    Original FORTRAN77 version by the LAPACK group.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
!    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
!    Sven Hammarling, Alan McKenney, Danny Sorensen,
!    LAPACK User's Guide,
!    Second Edition,
!    SIAM, 1995.
!
!  Parameters:
!
!    Input, integer M, the number of rows of the matrix A.  
!    0 <= M.
!
!    Input, integer N, the number of columns of the matrix A.  
!    0 <= N.
!
!    Input, integer ML, the number of subdiagonals within the
!    band of A.  0 <= ML.
!
!    Input, integer MU, the number of superdiagonals within 
!    the band of A.  0 <= MU.
!
!    Input/output, real ( kind = rk ) A(2*ML+MU+1,N).  On input, the matrix A
!    in band storage, and on output, information about the PLU factorization.
!
!    Output, integer PIVOT(min(M,N)), the pivot indices;
!    for 1 <= i <= min(M,N), row i of the matrix was interchanged with
!    row IPIV(i).
!
!    Output, integer INFO, error flag.
!    = 0: successful exit;
!    < 0: an input argument was illegal;
!    > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
!         has been completed, but the factor U is exactly
!         singular, and division by zero will occur if it is used
!         to solve a system of equations.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer m
  integer n

  real ( kind = rk ) a(2*ml+mu+1,n)
  integer i
  integer info
  integer pivot(*)
  integer j
  integer jp
  integer ju
  integer k
  integer km
  integer kv
  real ( kind = rk ) piv
  real ( kind = rk ) t
!
  info = 0
!
!  KV is the number of superdiagonals in the factor U, allowing for fill-in.
!
  kv = mu + ml
!
!  Set fill-in elements in columns MU+2 to KV to zero.
!
  do j = mu + 2, min ( kv, n )
    do i = kv - j + 2, ml
      a(i,j) = 0.0D+00
    end do
  end do
!
!  JU is the index of the last column affected by the current stage
!  of the factorization.
!
  ju = 1

  do j = 1, min ( m, n )
!
!  Set the fill-in elements in column J+KV to zero.
!
    if ( j + kv <= n ) then
      a(1:ml,j+kv) = 0.0D+00
    end if
!
!  Find the pivot and test for singularity.
!  KM is the number of subdiagonal elements in the current column.
!
    km = min ( ml, m-j )

    piv = abs ( a(kv+1,j) )
    jp = kv + 1

    do i = kv + 2, kv + km + 1
      if ( piv < abs ( a(i,j) ) ) then
        piv = abs ( a(i,j) )
        jp = i
      end if
    end do

    jp = jp - kv

    pivot(j) = jp + j - 1

    if ( a(kv+jp,j) /= 0.0D+00 ) then

      ju = max ( ju, min ( j + mu + jp - 1, n ) )
!
!  Apply interchange to columns J to JU.
!
      if ( jp /= 1 ) then

        do i = 0, ju - j
          t = a(kv+jp-i,j+i)
          a(kv+jp-i,j+i) = a(kv+1-i,j+i)
          a(kv+1-i,j+i) = t
        end do

      end if
!
!  Compute the multipliers.
!
      if ( 0 < km ) then

        a(kv+2:kv+km+1,j) = a(kv+2:kv+km+1,j) / a(kv+1,j)
!
!  Update the trailing submatrix within the band.
!
        if ( j < ju ) then

          do k = 1, ju - j

            if ( a(kv+1-k,j+k) /= 0.0D+00 ) then

              do i = 1, km
                a(kv+i+1-k,j+k) = a(kv+i+1-k,j+k) - a(kv+i+1,j) * a(kv+1-k,j+k)
              end do

            end if

          end do

        end if

      end if

    else
!
!  If pivot is zero, set INFO to the index of the pivot
!  unless a zero pivot has already been found.
!
      if ( info == 0 ) then
        info = j
      end if

    end if

  end do

  return
end
subroutine r8gb_trs ( n, ml, mu, nrhs, trans, a, pivot, b, info )

!*****************************************************************************80
!
!! R8GB_TRS solves an R8GB linear system factored by R8GB_TRF.
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
!    19 January 1999
!
!  Author:
!
!    Original FORTRAN77 version by the LAPACK group.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
!    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
!    Sven Hammarling, Alan McKenney, Danny Sorensen,
!    LAPACK User's Guide,
!    Second Edition,
!    SIAM, 1995.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix A.
!    N must be positive.
!
!    Input, integer ML, the number of subdiagonals within the 
!    band of A.  ML must be at least 0, and no greater than N - 1.
!
!    Input, integer MU, the number of superdiagonals within the 
!    band of A.  MU must be at least 0, and no greater than N - 1.
!
!    Input, integer NRHS, the number of right hand sides and the 
!    number of columns of the matrix B.  NRHS must be positive.
!
!    Input, character TRANS, specifies the form of the system.
!    'N':  A * x = b  (No transpose)
!    'T':  A'* X = B  (Transpose)
!    'C':  A'* X = B  (Conjugate transpose = Transpose)
!
!    Input, real ( kind = rk ) A(2*ML+MU+1,N), the LU factorization of the 
!    band matrix A, computed by R8GB_TRF.  
!
!    Input, integer PIVOT(N), the pivot indices; for 1 <= I <= N, 
!    row I of the matrix was interchanged with row PIVOT(I).
!
!    Input/output, real ( kind = rk ) B(N,NRHS),
!    On entry, the right hand side vectors B.
!    On exit, the solution vectors, X.
!
!    Output, integer INFO, error flag.
!    = 0:  successful exit
!    < 0: if INFO = -K, the K-th argument had an illegal value
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer n
  integer nrhs

  real ( kind = rk ) a(2*ml+mu+1,n)
  real ( kind = rk ) b(n,nrhs)
  integer i
  integer info
  integer pivot(*)
  integer j
  integer k
  integer kd
  integer l
  integer lm
  real ( kind = rk ) t(nrhs)
  real ( kind = rk ) temp
  character trans
!
!  Test the input parameters.
!
  info = 0

  if ( trans /= 'N' .and. trans /= 'n' .and. &
       trans /= 'T' .and. trans /= 't' .and. &
       trans /= 'C' .and. trans /= 'c' ) then
    info = -1
  else if ( n <= 0 ) then
    info = -2
  else if ( ml < 0 ) then
    info = -3
  else if ( mu < 0 ) then
    info = -4
  else if ( nrhs <= 0 ) then
    info = -5
  end if

  if ( info /= 0 ) then
    return
  end if

  kd = mu + ml + 1
!
!  Solve A * x = b.
!
!  Solve L * x = b, overwriting b with x.
!
!  L is represented as a product of permutations and unit lower
!  triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
!  where each transformation L(i) is a rank-one modification of
!  the identity matrix.
!
  if ( trans == 'N' .or. trans == 'n' ) then

    if ( 0 < ml ) then

      do j = 1, n - 1

        lm = min ( ml, n-j )
        l = pivot(j)

        t(1:nrhs)   = b(l,1:nrhs)
        b(l,1:nrhs) = b(j,1:nrhs)
        b(j,1:nrhs) = t(1:nrhs)

        do k = 1, nrhs
          if ( b(j,k) /= 0.0D+00 ) then
            b(j+1:j+lm,k) = b(j+1:j+lm,k) - a(kd+1:kd+lm,j) * b(j,k)
          end if
        end do

      end do

    end if
!
!  Solve U * x = b, overwriting b with x.
!
    do i = 1, nrhs

      do j = n, 1, -1
        if ( b(j,i) /= 0.0D+00 ) then
          l = ml + mu + 1 - j
          b(j,i) = b(j,i) / a(ml+mu+1,j)
          do k = j - 1, max ( 1, j - ml - mu ), -1
            b(k,i) = b(k,i) - a(l+k,j) * b(j,i)
          end do
        end if
      end do

    end do

  else
!
!  Solve A' * x = b.
!
!  Solve U' * x = b, overwriting b with x.
!
    do i = 1, nrhs

      do j = 1, n
        temp = b(j,i)
        l = ml + mu + 1 - j
        do k = max ( 1, j - ml - mu ), j - 1
          temp = temp - a(l+k,j) * b(k,i)
        end do
        b(j,i) = temp / a(ml+mu+1,j)
      end do

    end do
!
!  Solve L' * x = b, overwriting b with x.
!
    if ( 0 < ml ) then

      do j = n - 1, 1, -1

        lm = min ( ml, n - j )

        do k = 1, nrhs
          b(j,k) = b(j,k) - sum ( b(j+1:j+lm,k) * a(kd+1:kd+lm,j) )
        end do

        l = pivot(j)

        t(1:nrhs)   = b(l,1:nrhs)
        b(l,1:nrhs) = b(j,1:nrhs)
        b(j,1:nrhs) = t(1:nrhs)

      end do

    end if

  end if

  return
end
subroutine r8gb_zeros ( m, n, ml, mu, a )

!*****************************************************************************80
!
!! R8GB_ZEROS zeroes an R8GB matrix.
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
!    LINPACK and LAPACK band storage requires that an extra ML
!    superdiagonals be supplied to allow for fillin during Gauss
!    elimination.  Even though a band matrix is described as
!    having an upper bandwidth of MU, it effectively has an
!    upper bandwidth of MU+ML.  This routine assumes it is setting
!    up an unfactored matrix, so it only uses the first MU upper bands,
!    and does not place nonzero values in the fillin bands.
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
!    Input, integer M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than min(M,N)-1.
!
!    Output, real ( kind = rk ) A(2*ML+MU+1,N), the R8GB matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer n

  real ( kind = rk ) a(2*ml+mu+1,n)
  integer m

  a(1:2*ml+mu+1,1:n) = 0.0D+00

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
subroutine r8ge_random ( m, n, seed, a )

!*****************************************************************************80
!
!! R8GE_RANDOM randomizes an R8GE matrix.
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
!    11 August 2004
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
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns in
!    the array.
!
!    Input/output, integer SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = rk ) A(M,N), the array.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  integer i
  integer, parameter :: i4_huge = 2147483647
  integer j
  integer k
  integer seed

  do j = 1, n

    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge
      end if

      a(i,j) = real ( seed, kind = rk ) * 4.656612875D-10

    end do
  end do

  return
end
subroutine r8ge_to_r8gb ( m, n, ml, mu, a, b )

!*****************************************************************************80
!
!! R8GE_TO_R8GB copies an R8GE matrix to an R8GB matrix.
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
!    It usually doesn't make sense to try to store a general matrix
!    in a band matrix format.  You can always do it, but it will take
!    more space, unless the general matrix is actually banded.
!
!    The purpose of this routine is to allow a user to set up a
!    banded matrix in the easy-to-use general format, and have this
!    routine take care of the compression of the data into general
!    format.  All the user has to do is specify the bandwidths.
!
!    Note that this routine "believes" what the user says about the
!    bandwidth.  It will assume that all entries in the general matrix
!    outside of the bandwidth are zero.
!
!    The original M by N matrix is "collapsed" downward, so that diagonals
!    become rows of the storage array, while columns are preserved.  The
!    collapsed array is logically 2*ML+MU+1 by N.
!
!    The R8GB storage format is for an M by N banded matrix, with lower
!    bandwidth ML and upper bandwidth MU.  Storage includes room for ML
!    extra superdiagonals, which may be required to store nonzero entries
!    generated during Gaussian elimination.
!
!    LINPACK and LAPACK band storage requires that an extra ML
!    superdiagonals be supplied to allow for fillin during Gauss
!    elimination.  Even though a band matrix is described as
!    having an upper bandwidth of MU, it effectively has an
!    upper bandwidth of MU+ML.  This routine will copy nonzero
!    values it finds in these extra bands, so that both unfactored
!    and factored matrices can be handled.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
!    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
!    Sven Hammarling, Alan McKenney, Danny Sorensen,
!    LAPACK User's Guide,
!    Second Edition,
!    SIAM, 1995.
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!  Parameters:
!
!    Input, integer M, the number of rows of the matrices.
!    M must be positive.
!
!    Input, integer N, the number of columns of the matrices.
!    N must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths of 
!    the matrix.  ML and MU must be nonnegative, and no greater than min(M,N)-1.
!
!    Input, real ( kind = rk ) A(M,N), the R8GE matrix.
!
!    Output, real ( kind = rk ) B(2*ML+MU+1,N), the R8GB matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer ml
  integer mu
  integer n

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) b(2*ml+mu+1,n)
  integer i
  integer j
  integer jhi
  integer jlo

  b(1:2*ml+mu+1,1:n) = 0.0D+00

  do i = 1, m
    jlo = max ( i - ml, 1 )
    jhi = min ( i + mu, n )
    do j = jlo, jhi
      b(ml+mu+1+i-j,j) = a(i,j)
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
subroutine r8vec_to_r8gb ( m, n, ml, mu, x, a )

!*****************************************************************************80
!
!! R8VEC_TO_R8GB copies an R8VEC into an R8GB matrix.
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
!    Input, real ( kind = rk ) X((2*ML+MU+1)*N), the vector to be copied
!    into the array.
!
!    Output, real ( kind = rk ) A(2*ML+MU+1,N), the array.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer ml
  integer mu
  integer n

  real ( kind = rk ) a(2*ml+mu+1,n)
  integer i
  integer j
  real ( kind = rk ) x((2*ml+mu+1)*n)

  do j = 1, n
    do i = 1, 2 * ml + mu + 1

      if ( 1 <= i + j - ml - mu - 1 .and. i + j - ml - mu - 1 <= m ) then
        a(i,j) = x(i+(2*ml+mu+1)*(j-1))
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

