program main

!*****************************************************************************80
!
!! linpack_bench() runs the real double precision linpack() benchmark.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 February 2025
!
!  Author:
!
!    Original Fortran77 version by Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart.
!    This version by John Burkardt.
!
!  Reference:
!
!    Jack Dongarra,
!    Performance of Various Computers Using Standard Linear Equations Software,
!    Technical Report CS-89-85,
!    Electrical Engineering and Computer Science Department,
!    University of Tennessee, 2008.
!
!  Local:
!
!    integer n: the problem size.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )
  
  integer, parameter :: n = 1000
  integer, parameter :: lda = n + 1

  real ( kind = rk8 ), dimension ( lda, n ) :: A
  real ( kind = rk8 ) A_norm
  real ( kind = rk8 ), dimension ( lda, n ) :: ALU
  real ( kind = rk8 ) b(n)
  real ( kind = rk8 ) cray
  real ( kind = rk8 ) cray_ratio
  real ( kind = rk8 ) eps
  integer info
  integer ipvt(n)
  integer job
  real ( kind = rk8 ) mflops
  real ( kind = rk8 ) ops
  real ( kind = rk8 ) r(n)
  real ( kind = rk8 ) r_norm
  real ( kind = rk8 ) r8mat_norm_li 
  real ( kind = rk8 ) ratio
  real ( kind = rk8 ) t
  real ( kind = rk8 ) t1
  real ( kind = rk8 ) t2
  real ( kind = rk8 ) unit
  real ( kind = rk8 ) x(n)
  real ( kind = rk8 ) x_exact(n)
  real ( kind = rk8 ) x_norm

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'linpack_bench():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Datatype: Real Double Precision'
  write ( *, '(a,i8)' ) '  Matrix order N =               ', n
  write ( *, '(a,i8)' ) '  Leading matrix dimension LDA = ', lda

  call d_matgen ( A, lda, n )
  ALU(1:lda,1:n) = A(1:lda,1:n)
  x_exact(1:n) = 1.0D+00
  b(1:n) = matmul ( A(1:n,1:n), x_exact(1:n) )

  call cpu_time ( t1 )
  call dgefa ( ALU, lda, n, ipvt, info )
  call cpu_time ( t2 )
  t = t2 - t1

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'linpack_bench(): Fatal error!'
    write ( *, '(a)' ) '  The matrix A is apparently singular.'
    write ( *, '(a)' ) '  Abnormal end of execution.'
    stop 1
  end if

  x(1:n) = b(1:n)
  job = 0
  call cpu_time ( t1 )
  call dgesl ( ALU, lda, n, ipvt, x, job )
  call cpu_time ( t2 )

  t = t + t2 - t1
!
!  Compute a residual to verify results.
!
  r = b - matmul ( A(1:n,1:n), x(1:n) )
!
!  Compute infinity norms.
!
  A_norm = r8mat_norm_li ( lda, n, n, A )
  r_norm = maxval ( abs ( r(1:n) ) )
  x_norm = maxval ( abs ( x(1:n) ) )
!
!  Report.
!
  eps = epsilon ( eps )
  ratio = r_norm / A_norm / x_norm / eps

  ops = real ( 2 * n * n * n, kind = rk8 ) / 3.0D+00 &
    + 2.0D+00 * real ( n * n, kind = rk8 )
  mflops = ops / ( 1000000 * t )
  unit = 2.0 / mflops
  cray = 0.056;
  cray_ratio = t / cray;

  write ( *, '(a)' ) ''
  write ( *, '(a,f14.6)' ) '  Normalized residual = ', ratio
  write ( *, '(a,f14.6)' ) '  Residual norm       = ', r_norm
  write ( *, '(a,f14.6)' ) '  Machine epsilon     = ', eps
  write ( *, '(a,f14.6)' ) '  First X[]           = ', x(1)
  write ( *, '(a,f14.6)' ) '  Last X[]            = ', x(n)
  write ( *, '(a,f14.6)' ) '  Time in seconds     = ', t
  write ( *, '(a,f14.6)' ) '  MegaFLOPS           = ', mflops
  write ( *, '(a,f14.6)' ) '  Unit                = ', unit
  write ( *, '(a,f14.6)' ) '  Cray ratio          = ', cray_ratio
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'linpack_bench():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine d_matgen ( a, lda, n )

!*****************************************************************************80
!
!! d_matgen() generates a random matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 April 2003
!
!  Author:
!
!    Original Fortran77 version by Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart.
!    This version by John Burkardt.
!
!  Parameters:
!
!    Output, real ( kind = rk8 ) A(LDA,N), the N by N matrix.
!
!    Input, integer LDA, the leading dimension of the matrix.
!
!    Input, integer N, the order of the matrix.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer lda
  integer n

  real ( kind = rk8 ) a(lda,n)
  real ( kind = rk8 ) d_random
  integer i
  integer init(4)
  integer j

  init(1:4) = (/ 1, 2, 3, 1325 /)

  do j = 1, n
     do i = 1, n
       a(i,j) = d_random ( init ) - 0.5D+00
    end do
  end do

  return
end
function d_random ( iseed )

!*****************************************************************************80
!
!! d_random() returns a uniformly distributed random number between 0 and 1.
!
!  Discussion:
!
!    This routine uses a multiplicative congruential method with modulus
!    2**48 and multiplier 33952834046453.
!
!    48-bit integers are stored in 4 integer array elements with 12 bits
!    per element.  Hence the routine is portable across machines with
!    integers of 32 bits or more.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 February 2025
!
!  Author:
!
!    Original Fortran77 version by Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart.
!    This version by John Burkardt.
!
!  Reference:
!
!    GS Fishman,
!    Multiplicative congruential random number generators with modulus
!    2**b: an exhaustive analysis for b = 32 and a partial analysis for b = 48,
!    Mathematics of Computation,
!    Volume 189, 1990, pages 331-344.
!
!  Parameters:
!
!    Input/output, integer ISEED(4).
!    On entry, the seed of the random number generator; the array
!    elements must be between 0 and 4095, and ISEED(4) must be odd.
!    On exit, the seed is updated.
!
!    Output, real ( kind = rk8 ) D_RANDOM, the next pseudorandom number.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) d_random
  integer, parameter :: ipw2 = 4096
  integer iseed(4)
  integer it1
  integer it2
  integer it3
  integer it4
  integer, parameter :: m1 = 494
  integer, parameter :: m2 = 322
  integer, parameter :: m3 = 2508
  integer, parameter :: m4 = 2549
  real ( kind = rk8 ) , parameter :: one = 1.0D+00
  real ( kind = rk8 ) , parameter :: r = 1.0D+00 / 4096.0D+00
!
!  Multiply the seed by the multiplier modulo 2**48.
!
  it4 = iseed(4) * m4
  it3 = it4 / ipw2
  it4 = it4 - ipw2 * it3
  it3 = it3 + iseed(3) * m4 + iseed(4) * m3
  it2 = it3 / ipw2
  it3 = it3 - ipw2 * it2
  it2 = it2 + iseed(2) * m4 + iseed(3) * m3 + iseed(4) * m2
  it1 = it2 / ipw2
  it2 = it2 - ipw2 * it1
  it1 = it1 + iseed(1) * m4 + iseed(2) * m3 + iseed(3) * m2 + iseed(4) * m1
  it1 = mod ( it1, ipw2 )
!
!  Return updated seed
!
  iseed(1) = it1
  iseed(2) = it2
  iseed(3) = it3
  iseed(4) = it4
!
!  Convert 48-bit integer to a real number in the interval (0,1)
!
  d_random = &
      r * ( real ( it1, kind = rk8 ) &
    + r * ( real ( it2, kind = rk8 ) &
    + r * ( real ( it3, kind = rk8 ) &
    + r * ( real ( it4, kind = rk8 ) ) ) ) )

  return
end
subroutine daxpy ( n, sa, x, incx, y, incy )

!*****************************************************************************80
!
!! daxpy() adds a constant times one vector to another.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 April 1999
!
!  Author:
!
!    Original Fortran77 version by Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart.
!    This version by John Burkardt.
!
!  Reference:
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Algorithm 539: 
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real ( kind = rk8 ) SA, the multiplier.
!
!    Input, real ( kind = rk8 ) X(*), the vector to be scaled and added to Y.
!
!    Input, integer INCX, the increment between successive entries of X.
!
!    Input/output, real ( kind = rk8 ) Y(*), the vector to which a 
!    multiple of X is to be added.
!
!    Input, integer INCY, the increment between successive entries of Y.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer i
  integer incx
  integer incy
  integer ix
  integer iy
  integer n
  real ( kind = rk8 ) sa
  real ( kind = rk8 ) x(*)
  real ( kind = rk8 ) y(*)

  if ( n <= 0 ) then

  else if ( sa == 0.0D+00 ) then

  else if ( incx == 1 .and. incy == 1 ) then

    y(1:n) = y(1:n) + sa * x(1:n)

  else

    if ( 0 <= incx ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    if ( 0 <= incy ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if

    do i = 1, n
      y(iy) = y(iy) + sa * x(ix)
      ix = ix + incx
      iy = iy + incy
    end do

  end if

  return
end
subroutine dgefa ( a, lda, n, ipvt, info )

!*****************************************************************************80
!
!! dgefa() factors a real matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 March 2001
!
!  Author:
!
!    Original Fortran77 version by Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart.
!    This version by John Burkardt.
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
!    Input/output, real ( kind = rk8 ) A(LDA,N).
!    On intput, the matrix to be factored.
!    On output, an upper triangular matrix and the multipliers used to obtain
!    it.  The factorization can be written A=L*U, where L is a product of
!    permutation and unit lower triangular matrices, and U is upper triangular.
!
!    Input, integer LDA, the leading dimension of A.
!
!    Input, integer N, the order of the matrix A.
!
!    Output, integer IPVT(N), the pivot indices.
!
!    Output, integer INFO, singularity indicator.
!    0, normal value.
!    K, if U(K,K) == 0.  This is not an error condition for this subroutine,
!    but it does indicate that DGESL or DGEDI will divide by zero if called.
!    Use RCOND in DGECO for a reliable indication of singularity.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer lda
  integer n

  real ( kind = rk8 ) a(lda,n)
  integer info
  integer ipvt(n)
  integer idamax
  integer j
  integer k
  integer l
  real ( kind = rk8 ) t
!
!  Gaussian elimination with partial pivoting.
!
  info = 0

  do k = 1, n - 1
!
!  Find L = pivot index.
!
    l = idamax ( n-k+1, a(k,k), 1 ) + k - 1
    ipvt(k) = l
!
!  Zero pivot implies this column already triangularized.
!
    if ( a(l,k) == 0.0D+00 ) then
      info = k
      cycle
    end if
!
!  Interchange if necessary.
!
    if ( l /= k ) then
      t      = a(l,k)
      a(l,k) = a(k,k)
      a(k,k) = t
    end if
!
!  Compute multipliers.
!
    a(k+1:n,k) = - a(k+1:n,k) / a(k,k)
!
!  Row elimination with column indexing.
!
    do j = k+1, n
      t = a(l,j)
      if ( l /= k ) then
        a(l,j) = a(k,j)
        a(k,j) = t
      end if
      call daxpy ( n-k, t, a(k+1,k), 1, a(k+1,j), 1 )
    end do

  end do

  ipvt(n) = n

  if ( a(n,n) == 0.0D+00 ) then
    info = n
  end if

  return
end
subroutine dgesl ( a, lda, n, ipvt, b, job )

!*****************************************************************************80
!
!! dgesl() solves a real general linear system A * X = B.
!
!  Discussion:
!
!    DGESL can solve either of the systems A * X = B or A' * X = B.
!
!    The system matrix must have been factored by DGECO or DGEFA.
!
!    A division by zero will occur if the input factor contains a
!    zero on the diagonal.  Technically this indicates singularity
!    but it is often caused by improper arguments or improper
!    setting of LDA.  It will not occur if the subroutines are
!    called correctly and if DGECO has set 0.0 < RCOND 
!    or DGEFA has set INFO == 0.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 March 2001
!
!  Author:
!
!    Original Fortran77 version by Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart.
!    This version by John Burkardt.
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
!    Input, real ( kind = rk8 ) A(LDA,N), the output from DGECO or DGEFA.
!
!    Input, integer LDA, the leading dimension of A.
!
!    Input, integer N, the order of the matrix A.
!
!    Input, integer IPVT(N), the pivot vector from DGECO or DGEFA.
!
!    Input/output, real ( kind = rk8 ) B(N).
!    On input, the right hand side vector.
!    On output, the solution vector.
!
!    Input, integer JOB.
!    0, solve A * X = B;
!    nonzero, solve A' * X = B.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer lda
  integer n

  real ( kind = rk8 ) a(lda,n)
  real ( kind = rk8 ) b(n)
  integer ipvt(n)
  integer job
  integer k
  integer l
  real ( kind = rk8 ) t
!
!  Solve A * X = B.
!
  if ( job == 0 ) then

    do k = 1, n - 1

      l = ipvt(k)
      t = b(l)

      if ( l /= k ) then
        b(l) = b(k)
        b(k) = t
      end if

      call daxpy ( n-k, t, a(k+1,k), 1, b(k+1), 1 )

    end do

    do k = n, 1, -1
      b(k) = b(k) / a(k,k)
      t = -b(k)
      call daxpy ( k-1, t, a(1,k), 1, b(1), 1 )
    end do

  else
!
!  Solve A' * X = B.
!
    do k = 1, n
      t = dot_product ( a(1:k-1,k), b(1:k-1) )
      b(k) = ( b(k) - t ) / a(k,k)
    end do

    do k = n - 1, 1, -1

      b(k) = b(k) + dot_product ( a(k+1:n,k), b(k+1:n) )
      l = ipvt(k)

      if ( l /= k ) then
        t    = b(l)
        b(l) = b(k)
        b(k) = t
      end if

    end do

  end if

  return
end
function idamax ( n, x, incx )

!*****************************************************************************80
!
!! idamax() finds the index of the vector element of maximum absolute value.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 April 1999
!
!  Author:
!
!    Original Fortran77 version by Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart.
!    This version by John Burkardt.
!
!  Reference:
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Algorithm 539: 
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!    Dongarra, Moler, Bunch and Stewart,
!    LINPACK User's Guide,
!    SIAM, (Society for Industrial and Applied Mathematics),
!    3600 University City Science Center,
!    Philadelphia, PA, 19104-2688.
!    ISBN 0-89871-172-X
!
!  Input:
!
!    integer N, the number of entries in the vector.
!
!    real ( kind = rk8 ) X(*), the vector to be examined.
!
!    integer INCX, the increment between successive entries of SX.
!
!  Output:
!
!    integer IDAMAX, the index of the element of SX of maximum absolute value.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer i
  integer incx
  integer idamax
  integer ix
  integer n
  real ( kind = rk8 ) samax
  real ( kind = rk8 ) x(*)

  if ( n <= 0 ) then

    idamax = 0

  else if ( n == 1 ) then

    idamax = 1

  else if ( incx == 1 ) then

    idamax = 1
    samax = abs ( x(1) )

    do i = 2, n

      if ( samax < abs ( x(i) ) ) then
        idamax = i
        samax = abs ( x(i) )
      end if

    end do

  else

    if ( 0 <= incx ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    idamax = 1
    samax = abs ( x(ix) )

    ix = ix + incx

    do i = 2, n
      if ( samax < abs ( x(ix) ) ) then
        idamax = i
        samax = abs ( x(ix) )
      end if
      ix = ix + incx
    end do

  end if

  return
end
function r8mat_norm_li ( lda, m, n, a )

!*****************************************************************************80
!
!! r8mat_norm_li() returns the matrix Loo norm of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    The matrix Loo norm is defined as:
!
!      R8MAT_NORM_LI =  max ( 1 <= I <= M ) sum ( 1 <= J <= N ) abs ( A(I,J) ).
!
!    The matrix Loo norm is derived from the vector Loo norm,
!    and satisifies:
!
!      r8vec_norm_li ( A * x ) <= r8mat_norm_li ( A ) * r8vec_norm_li ( x ).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, the number of rows in A.
!
!    integer N, the number of columns in A.
!
!    real ( kind = rk8 ) A(LDA,N), the matrix whose Loo
!    norm is desired.
!
!  Output:
!
!    real ( kind = rk8 ) R8MAT_NORM_LI, the Loo norm of A.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer lda
  integer m
  integer n

  real ( kind = rk8 ) a(lda,n)
  integer i
  real ( kind = rk8 ) r8mat_norm_li
  real ( kind = rk8 ) row_sum

  r8mat_norm_li = 0.0D+00

  do i = 1, m
    row_sum = sum ( abs ( a(i,1:n) ) )
    r8mat_norm_li = max ( r8mat_norm_li, row_sum )
  end do

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! timestamp() prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 8 ) ampm
  integer d
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = 10 )  time
  integer values(8)
  integer y
  character ( len = 5 ) zone

  call date_and_time ( date, time, zone, values )

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

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
