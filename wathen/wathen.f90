subroutine bandwidth ( m, n, a, b, l, d, u )

!*****************************************************************************80
!
!! bandwidth() returns the bandwidth of a matrix.
!
!  Discussion:
!
!    If the nonzeros of a matrix only occur in entries that are "close"
!    to the main diagonal, we say the matrix is banded.
!
!    Roughly speaking, the bandwidth B of a matrix is the number of 
!    diagonals containing nonzeros.  More precisely, it is the minimum number
!    of contiguous diagonals that contain all the nonzeros.  It is presumed
!    that the main diagonal is nonzero.
!
!    We can also measure U and L, the upper and lower "half-bandwidths" which
!    count the number of contiguous diagonals above or below the main
!    diagonal.
!
!    We may write
!      B = L + D + U
!    where D is presumably 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 June 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = rk ) A(M,N), the matrix.
!
!    Output, integer B, the total bandwidth.
!
!    Output, integer L, D, U, the lower, diagonal, and upper 
!    bandwidths.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  integer b
  integer d
  integer i
  integer j
  integer l
  integer u

  l = 0
  d = 0
  u = 0

  do i = 1, n

    j = 1
    do while ( l < i - j )
      if ( a(i,j) /= 0.0D+00 ) then
        l = i - j
        exit
      end if
      j = j + 1
    end do

    if ( a(i,i) /= 0.0D+00 ) then
      d = 1
    end if

    j = n
    do while ( u < j - i )
      if ( a(i,j) /= 0.0D+00 ) then
        u = j - i
        exit
      end if
      j = j - 1
    end do

  end do

  b = l + d + u

  return
end
subroutine cg_gb ( n, ml, mu, a, b, x )

!*****************************************************************************80
!
!! CG_GB uses the conjugate gradient method for a general banded (GB) matrix.
!
!  Discussion:
!
!    The linear system has the form A*x=b, where A is a positive-definite
!    symmetric matrix.
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
!    04 June 2014
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
!    Input, integer ML, MU, the lower and upper bandwidths.
!
!    Input, real ( kind = rk ) A(2*ML+MU+1,N), the band matrix.
!
!    Input, real ( kind = rk ) B(N), the right hand side vector.
!
!    Input/output, real ( kind = rk ) X(N).
!    On input, an estimate for the solution, which may be 0.
!    On output, the approximate solution vector.  
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer mu
  integer n

  real ( kind = rk ) a(2*ml+mu+1,n)
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
  call mv_gb ( n, n, ml, mu, a, x, ap )

  r(1:n) = b(1:n) - ap(1:n)
  p(1:n) = b(1:n) - ap(1:n)
!
!  Do the N steps of the conjugate gradient method.
!
  do it = 1, n
!
!  Compute the matrix*vector product AP = A*P.
!
    call mv_gb ( n, n, ml, mu, a, p, ap )
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
      exit
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
subroutine cg_ge ( n, a, b, x )

!*****************************************************************************80
!
!! CG_GE uses the conjugate gradient method for a general storage (GE) matrix.
!
!  Discussion:
!
!    The linear system has the form A*x=b, where A is a positive-definite
!    symmetric matrix, stored as a full storage matrix.
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
!    01 June 2014
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
!    Input, real ( kind = rk ) A(N,N), the matrix.
!
!    Input, real ( kind = rk ) B(N), the right hand side vector.
!
!    Input/output, real ( kind = rk ) X(N).
!    On input, an estimate for the solution, which may be 0.
!    On output,  the approximate solution vector.  
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
  ap = matmul ( a, x )
 
  r(1:n) = b(1:n) - ap(1:n)
  p(1:n) = b(1:n) - ap(1:n)
!
!  Do the N steps of the conjugate gradient method.
!
  do it = 1, n
!
!  Compute the matrix*vector product AP = A*P.
!
    ap = matmul ( a, p )
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
    rap = dot_product ( r(1:n), ap(1:n) )

    beta = - rap / pap
!
!  Update the perturbation vector
!    P = R + BETA * P.
!
    p(1:n) = r(1:n) + beta * p(1:n)

  end do

  return
end
subroutine cg_st ( n, nz_num, row, col, a, b, x )

!*****************************************************************************80
!
!! CG_ST uses the conjugate gradient method for a sparse triplet (ST) matrix.
!
!  Discussion:
!
!    The linear system has the form A*x=b, where A is a positive-definite
!    symmetric matrix, stored as a full storage matrix.
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
!    04 June 2014
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
!    Input, integer NZ_NUM, the number of nonzeros.
!
!    Input, integer ROW(NZ_NUM), COL(NZ_NUM), the row and column 
!    indices of the nonzero entries.
!
!    Input, real ( kind = rk ) A(NZ_NUM), the nonzero entries.
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
  real ( kind = rk ) rap
  integer row(nz_num)
  real ( kind = rk ) x(n)
!
!  Initialize
!    AP = A * x,
!    R  = b - A * x,
!    P  = b - A * x.
!
  call mv_st ( n, n, nz_num, row, col, a, x, ap )

  r(1:n) = b(1:n) - ap(1:n)
  p(1:n) = b(1:n) - ap(1:n)
!
!  Do the N steps of the conjugate gradient method.
!
  do it = 1, n
!
!  Compute the matrix*vector product AP = A*P.
!
    call mv_st ( n, n, nz_num, row, col, a, p, ap )
!
!  Compute the dot products
!    PAP = P*AP,
!    PR  = P*R
!  Set
!    ALPHA = PR / PAP.
!
    pap = dot_product ( p, ap )
    pr =  dot_product ( p, r )

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
subroutine daxpy ( n, da, dx, incx, dy, incy )

!*****************************************************************************80
!
!! DAXPY computes constant times a vector plus a vector.
!
!  Discussion:
!
!    This routine uses double precision real arithmetic.
!
!    This routine uses unrolled loops for increments equal to one.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 May 2005
!
!  Author:
!
!    Original FORTRAN77 version by Charles Lawson, Richard Hanson, 
!    David Kincaid, Fred Krogh.
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
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Algorithm 539, 
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    ACM Transactions on Mathematical Software, 
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer N, the number of elements in DX and DY.
!
!    Input, real ( kind = rk ) DA, the multiplier of DX.
!
!    Input, real ( kind = rk ) DX(*), the first vector.
!
!    Input, integer INCX, the increment between successive 
!    entries of DX.
!
!    Input/output, real ( kind = rk ) DY(*), the second vector.
!    On output, DY(*) has been replaced by DY(*) + DA * DX(*).
!
!    Input, integer INCY, the increment between successive 
!    entries of DY.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) da
  real ( kind = rk ) dx(*)
  real ( kind = rk ) dy(*)
  integer i
  integer incx
  integer incy
  integer ix
  integer iy
  integer m
  integer n

  if ( n <= 0 ) then
    return
  end if

  if ( da == 0.0D+00 ) then
    return
  end if
!
!  Code for unequal increments or equal increments
!  not equal to 1.
!
  if ( incx /= 1 .or. incy /= 1 ) then

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
      dy(iy) = dy(iy) + da * dx(ix)
      ix = ix + incx
      iy = iy + incy
    end do
!
!  Code for both increments equal to 1.
!
  else

    m = mod ( n, 4 )

    dy(1:m) = dy(1:m) + da * dx(1:m)

    do i = m + 1, n, 4
      dy(i  ) = dy(i  ) + da * dx(i  )
      dy(i+1) = dy(i+1) + da * dx(i+1)
      dy(i+2) = dy(i+2) + da * dx(i+2)
      dy(i+3) = dy(i+3) + da * dx(i+3)
    end do

  end if

  return
end
function ddot ( n, dx, incx, dy, incy )

!*****************************************************************************80
!
!! DDOT forms the dot product of two vectors.
!
!  Discussion:
!
!    This routine uses double precision real arithmetic.
!
!    This routine uses unrolled loops for increments equal to one.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 May 2005
!
!  Author:
!
!    Original FORTRAN77 version by Charles Lawson, Richard Hanson, 
!    David Kincaid, Fred Krogh.
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
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Algorithm 539, 
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    ACM Transactions on Mathematical Software, 
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vectors.
!
!    Input, real ( kind = rk ) DX(*), the first vector.
!
!    Input, integer INCX, the increment between successive 
!    entries in DX.
!
!    Input, real ( kind = rk ) DY(*), the second vector.
!
!    Input, integer INCY, the increment between successive 
!    entries in DY.
!
!    Output, real ( kind = rk ) DDOT, the sum of the product of the 
!    corresponding entries of DX and DY.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) ddot
  real ( kind = rk ) dtemp
  real ( kind = rk ) dx(*)
  real ( kind = rk ) dy(*)
  integer i
  integer incx
  integer incy
  integer ix
  integer iy
  integer m
  integer n

  ddot = 0.0D+00
  dtemp = 0.0D+00

  if ( n <= 0 ) then
    return
  end if
!
!  Code for unequal increments or equal increments
!  not equal to 1.
!
  if ( incx /= 1 .or. incy /= 1 ) then

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
      dtemp = dtemp + dx(ix) * dy(iy)
      ix = ix + incx
      iy = iy + incy
    end do
!
!  Code for both increments equal to 1.
!
  else

    m = mod ( n, 5 )

    do i = 1, m
      dtemp = dtemp + dx(i) * dy(i)
    end do

    do i = m + 1, n, 5

      dtemp = dtemp + dx(i  ) * dy(i  ) &
                    + dx(i+1) * dy(i+1) &
                    + dx(i+2) * dy(i+2) &
                    + dx(i+3) * dy(i+3) &
                    + dx(i+4) * dy(i+4)
    end do

  end if

  ddot = dtemp

  return
end
subroutine dgbfa ( abd, lda, n, ml, mu, ipvt, info )

!*****************************************************************************80
!
!! DGBFA factors a real band matrix by elimination.
!
!  Discussion:
!
!    DGBFA is usually called by DGBCO, but it can be called
!    directly with a saving in time if RCOND is not needed.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 May 2005
!
!  Author:
!
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
!    Input/output, real ( kind = rk ) ABD(LDA,N).  On input, the matrix in band
!    storage.  The columns of the matrix are stored in the columns of ABD
!    and the diagonals of the matrix are stored in rows ML+1 through
!    2*ML+MU+1 of ABD.  On output, an upper triangular matrix in band storage
!    and the multipliers which were used to obtain it.  The factorization
!    can be written A = L*U where L is a product of permutation and unit lower
!    triangular matrices and U is upper triangular.
!
!    Input, integer LDA, the leading dimension of the array ABD.
!    2*ML + MU + 1 <= LDA is required.
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer ML, MU, the number of diagonals below and above
!    the main diagonal.  0 <= ML < N, 0 <= MU < N.
!
!    Output, integer IPVT(N), the pivot indices.
!
!    Output, integer INFO, error flag.
!    0, normal value.
!    K, if U(K,K) == 0.0D+00.  This is not an error condition for this
!      subroutine, but it does indicate that DGBSL will divide by zero if
!      called.  Use RCOND in DGBCO for a reliable indication of singularity.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer lda
  integer n

  real ( kind = rk ) abd(lda,n)
  integer i
  integer i0
  integer info
  integer ipvt(n)
  integer idamax
  integer j
  integer j0
  integer j1
  integer ju
  integer jz
  integer k
  integer l
  integer lm
  integer m
  integer ml
  integer mm
  integer mu
  real ( kind = rk ) t

  m = ml + mu + 1
  info = 0
!
!  Zero initial fill-in columns.
!
  j0 = mu + 2
  j1 = min ( n, m ) - 1

  do jz = j0, j1
    i0 = m + 1 - jz
    do i = i0, ml
      abd(i,jz) = 0.0D+00
    end do
  end do

  jz = j1
  ju = 0
!
!  Gaussian elimination with partial pivoting.
!
  do k = 1, n-1
!
!  Zero out the next fill-in column.
!
    jz = jz + 1
    if ( jz <= n ) then
      abd(1:ml,jz) = 0.0D+00
    end if
!
!  Find L = pivot index.
!
    lm = min ( ml, n-k )
    l = idamax ( lm+1, abd(m,k), 1 ) + m - 1
    ipvt(k) = l + k - m
!
!  Zero pivot implies this column already triangularized.
!
    if ( abd(l,k) == 0.0D+00 ) then

      info = k
!
!  Interchange if necessary.
!
    else

      if ( l /= m ) then
        t = abd(l,k)
        abd(l,k) = abd(m,k)
        abd(m,k) = t
      end if
!
!  Compute multipliers.
!
      t = -1.0D+00 / abd(m,k)
      call dscal ( lm, t, abd(m+1,k), 1 )
!
!  Row elimination with column indexing.
!
      ju = min ( max ( ju, mu + ipvt(k) ), n )
      mm = m

      do j = k+1, ju
        l = l - 1
        mm = mm - 1
        t = abd(l,j)
        if ( l /= mm ) then
          abd(l,j) = abd(mm,j)
          abd(mm,j) = t
        end if
        call daxpy ( lm, t, abd(m+1,k), 1, abd(mm+1,j), 1 )
      end do

    end if

  end do

  ipvt(n) = n

  if ( abd(m,n) == 0.0D+00 ) then
    info = n
  end if

  return
end
subroutine dgbsl ( abd, lda, n, ml, mu, ipvt, b, job )

!*****************************************************************************80
!
!! DGBSL solves a real banded system factored by DGBCO or DGBFA.
!
!  Discussion:
!
!    DGBSL can solve either A * X = B  or  A' * X = B.
!
!    A division by zero will occur if the input factor contains a
!    zero on the diagonal.  Technically this indicates singularity
!    but it is often caused by improper arguments or improper
!    setting of LDA.  It will not occur if the subroutines are
!    called correctly and if DGBCO has set 0.0 < RCOND
!    or DGBFA has set INFO == 0.
!
!    To compute inverse(A) * C  where C is a matrix with P columns:
!
!      call dgbco ( abd, lda, n, ml, mu, ipvt, rcond, z )
!
!      if ( rcond is too small ) then
!        exit
!      end if
!
!      do j = 1, p
!        call dgbsl ( abd, lda, n, ml, mu, ipvt, c(1,j), 0 )
!      end do
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 May 2005
!
!  Author:
!
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
!    Input, real ( kind = rk ) ABD(LDA,N), the output from DGBCO or DGBFA.
!
!    Input, integer LDA, the leading dimension of the array ABD.
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer ML, MU, the number of diagonals below and above
!    the main diagonal.  0 <= ML < N, 0 <= MU < N.
!
!    Input, integer IPVT(N), the pivot vector from DGBCO or DGBFA.
!
!    Input/output, real ( kind = rk ) B(N).  On input, the right hand side.
!    On output, the solution.
!
!    Input, integer JOB, job choice.
!    0, solve A*X=B.
!    nonzero, solve A'*X=B.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer lda
  integer n

  real ( kind = rk ) abd(lda,n)
  real ( kind = rk ) b(n)
  integer ipvt(n)
  integer job
  integer k
  integer l
  integer la
  integer lb
  integer lm
  integer m
  integer ml
  integer mu
  real ( kind = rk ) ddot
  real ( kind = rk ) t

  m = mu + ml + 1
!
!  JOB = 0, Solve A * x = b.
!
!  First solve L * y = b.
!
  if ( job == 0 ) then

    if ( 0 < ml ) then

      do k = 1, n-1
        lm = min ( ml, n-k )
        l = ipvt(k)
        t = b(l)
        if ( l /= k ) then
          b(l) = b(k)
          b(k) = t
        end if
        call daxpy ( lm, t, abd(m+1,k), 1, b(k+1), 1 )
      end do

    end if
!
!  Now solve U * x = y.
!
    do k = n, 1, -1
      b(k) = b(k) / abd(m,k)
      lm = min ( k, m ) - 1
      la = m - lm
      lb = k - lm
      t = -b(k)
      call daxpy ( lm, t, abd(la,k), 1, b(lb), 1 )
    end do
!
!  JOB nonzero, solve A' * x = b.
!
!  First solve U' * y = b.
!
  else

    do k = 1, n
      lm = min ( k, m ) - 1
      la = m - lm
      lb = k - lm
      t = ddot ( lm, abd(la,k), 1, b(lb), 1 )
      b(k) = ( b(k) - t ) / abd(m,k)
    end do
!
!  Now solve L' * x = y.
!
    if ( 0 < ml ) then

      do k = n - 1, 1, -1
        lm = min ( ml, n - k )
        b(k) = b(k) + ddot ( lm, abd(m+1,k), 1, b(k+1), 1 )
        l = ipvt(k)
        if ( l /= k ) then
          t = b(l)
          b(l) = b(k)
          b(k) = t
        end if
      end do

    end if

  end if

  return
end
subroutine dgefa ( a, lda, n, ipvt, info )

!*****************************************************************************80
!
!! DGEFA factors a real general matrix.
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
!    Input/output, real ( kind = rk ) A(LDA,N).
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

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer lda
  integer n

  real ( kind = rk ) a(lda,n)
  integer info
  integer ipvt(n)
  integer idamax
  integer j
  integer k
  integer l
  real ( kind = rk ) t
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
      t = a(l,k)
      a(l,k) = a(k,k)
      a(k,k) = t
    end if
!
!  Compute multipliers.
!
    t = -1.0D+00 / a(k,k)
    a(k+1:n,k) = a(k+1:n,k) * t
!
!  Row elimination with column indexing.
!
    do j = k+1, n
      t = a(l,j)
      if ( l /= k ) then
        a(l,j) = a(k,j)
        a(k,j) = t
      end if
      a(k+1:n,j) = a(k+1:n,j) + t * a(k+1:n,k)
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
!! DGESL solves a real general linear system A * X = B.
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
!    Input, real ( kind = rk ) A(LDA,N), the output from DGECO or DGEFA.
!
!    Input, integer LDA, the leading dimension of A.
!
!    Input, integer N, the order of the matrix A.
!
!    Input, integer IPVT(N), the pivot vector from DGECO or DGEFA.
!
!    Input/output, real ( kind = rk ) B(N).
!    On input, the right hand side vector.
!    On output, the solution vector.
!
!    Input, integer JOB.
!    0, solve A * X = B;
!    nonzero, solve A' * X = B.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer lda
  integer n

  real ( kind = rk ) a(lda,n)
  real ( kind = rk ) b(n)
  integer ipvt(n)
  integer job
  integer k
  integer l
  real ( kind = rk ) t
!
!  Solve A * X = B.
!
  if ( job == 0 ) then

    do k = 1, n-1

      l = ipvt(k)
      t = b(l)

      if ( l /= k ) then
        b(l) = b(k)
        b(k) = t
      end if

      b(k+1:n) = b(k+1:n) + t * a(k+1:n,k)

    end do

    do k = n, 1, -1
      b(k) = b(k) / a(k,k)
      t = -b(k)
      b(1:k-1) = b(1:k-1) + t * a(1:k-1,k)
    end do

  else
!
!  Solve A' * X = B.
!
    do k = 1, n
      t = dot_product ( a(1:k-1,k), b(1:k-1) )
      b(k) = ( b(k) - t ) / a(k,k)
    end do

    do k = n-1, 1, -1

      b(k) = b(k) + dot_product ( a(k+1:n,k), b(k+1:n) )
      l = ipvt(k)

      if ( l /= k ) then
        t = b(l)
        b(l) = b(k)
        b(k) = t
      end if

    end do

  end if

  return
end
subroutine dscal ( n, sa, x, incx )

!*****************************************************************************80
!
!! DSCAL scales a vector by a constant.
!
!  Discussion:
!
!    This routine uses double precision real arithmetic.
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
!    Original FORTRAN77 version by Charles Lawson, Richard Hanson, 
!    David Kincaid, Fred Krogh.
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
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Algorithm 539, 
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real ( kind = rk ) SA, the multiplier.
!
!    Input/output, real ( kind = rk ) X(*), the vector to be scaled.
!
!    Input, integer INCX, the increment between successive 
!    entries of X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer incx
  integer ix
  integer m
  integer n
  real ( kind = rk ) sa
  real ( kind = rk ) x(*)

  if ( n <= 0 ) then

  else if ( incx == 1 ) then

    m = mod ( n, 5 )

    x(1:m) = sa * x(1:m)

    do i = m+1, n, 5
      x(i)   = sa * x(i)
      x(i+1) = sa * x(i+1)
      x(i+2) = sa * x(i+2)
      x(i+3) = sa * x(i+3)
      x(i+4) = sa * x(i+4)
    end do

  else

    if ( 0 <= incx ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    do i = 1, n
      x(ix) = sa * x(ix)
      ix = ix + incx
    end do

  end if

  return
end
function idamax ( n, dx, incx )

!*****************************************************************************80
!
!! IDAMAX indexes the array element of maximum absolute value.
!
!  Discussion:
!
!    This routine uses double precision real arithmetic.
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
!    Original FORTRAN77 version by Charles Lawson, Richard Hanson, 
!    David Kincaid, Fred Krogh.
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
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Algorithm 539, 
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real ( kind = rk ) X(*), the vector to be examined.
!
!    Input, integer INCX, the increment between successive 
!    entries of SX.
!
!    Output, integer IDAMAX, the index of the element of SX of 
!    maximum absolute value.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) dmax
  real ( kind = rk ) dx(*)
  integer i
  integer idamax
  integer incx
  integer ix
  integer n

  idamax = 0

  if ( n < 1 .or. incx <= 0 ) then
    return
  end if

  idamax = 1

  if ( n == 1 ) then
    return
  end if

  if ( incx == 1 ) then

    dmax = abs ( dx(1) )

    do i = 2, n
      if ( dmax < abs ( dx(i) ) ) then
        idamax = i
        dmax = abs ( dx(i) )
      end if
    end do

  else

    ix = 1
    dmax = abs ( dx(1) )
    ix = ix + incx

    do i = 2, n
      if ( dmax < abs ( dx(ix) ) ) then
        idamax = i
        dmax = abs ( dx(ix) )
      end if
      ix = ix + incx
    end do

  end if

  return
end
subroutine mv_gb ( m, n, ml, mu, a, x, b )

!*****************************************************************************80
!
!! MV_GB multiplies a banded matrix by an R8VEC.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 June 2014
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
!
!    Input, real ( kind = rk ) A(2*ML+MU+1,N), the matrix.
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

  real ( kind = rk ) a(2*ml+mu+1,n)
  real ( kind = rk ) b(m)
  integer i
  integer j
  integer jhi
  integer jlo
  real ( kind = rk ) x(n)

  b(1:m) = 0.0D+00

  do i = 1, n
    jlo = max ( 1, i - ml )
    jhi = min ( n, i + mu )
    do j = jlo, jhi
      b(i) = b(i) + a(i-j+ml+mu+1,j) * x(j)
    end do
  end do

  return
end
subroutine mv_ge ( m, n, a, x, b )

!*****************************************************************************80
!
!! MV_GE multiplies an R8GE matrix by an R8VEC.
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
subroutine mv_st ( m, n, nz_num, row, col, a, x, b )

!*****************************************************************************80
!
!! MV_ST multiplies a sparse triple matrix times a vector.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 June 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer NZ_NUM, the number of nonzero values.
!
!    Input, integer ROW(NZ_NUM), COL(NZ_NUM), the row and 
!    column indices.
!
!    Input, real ( kind = rk ) A(NZ_NUM), the nonzero values in the matrix.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied.
!
!    Output, real ( kind = rk ) B(M), the product A*X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n
  integer nz_num

  real ( kind = rk ) a(nz_num)
  real ( kind = rk ) b(m)
  integer col(nz_num)
  integer k
  integer row(nz_num)
  real ( kind = rk ) x(n)

  b(1:m) = 0.0D+00
  do k = 1, nz_num
    b(row(k)) = b(row(k)) + a(k) * x(col(k))
  end do

  return
end
subroutine nonzeros ( m, n, a, nnz )

!*****************************************************************************80
!
!! NONZEROS counts the nonzeros in a matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 June 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = rk ) A(M,N), the matrix.
!
!    Output, integer NNZ, the number of nonzero entries.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  integer i
  integer j
  integer nnz

  nnz = 0
  do j = 1, n
    do i = 1, m
      if ( a(i,j) /= 0.0D+00 ) then
        nnz = nnz + 1
      end if
    end do
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
subroutine wathen_bandwidth ( nx, ny, l, d, u )

!*****************************************************************************80
!
!! WATHEN_BANDWIDTH returns the bandwidth of the WATHEN matrix.
!
!  Discussion:
!
!    The bandwidth measures the minimal number of contiguous diagonals,
!    including the central diagonal, which contain all the nonzero elements
!    of a matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 June 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Nicholas Higham,
!    Algorithm 694: A Collection of Test Matrices in MATLAB,
!    ACM Transactions on Mathematical Software,
!    Volume 17, Number 3, September 1991, pages 289-305.
!
!    Andrew Wathen,
!    Realistic eigenvalue bounds for the Galerkin mass matrix,
!    IMA Journal of Numerical Analysis,
!    Volume 7, 1987, pages 449-457.
!
!  Parameters:
!
!    Input, integer NX, NY, values which determine the size of A.
!
!    Output, integer L, D, U, the lower, diagonal, and upper 
!    bandwidths of the matrix,
!
  implicit none

  integer d
  integer l
  integer nx
  integer ny
  integer u

  l = 3 * nx + 4
  d = 1
  u = 3 * nx + 4

  return
end
subroutine wathen_gb ( nx, ny, n, a )

!*****************************************************************************80
!
!! WATHEN_GB returns the Wathen matrix, using general banded (GB) storage.
!
!  Discussion:
!
!    The Wathen matrix is a finite element matrix which is sparse.
!
!    The entries of the matrix depend in part on a physical quantity
!    related to density.  That density is here assigned random values between
!    0 and 100.
!
!    The matrix order N is determined by the input quantities NX and NY,
!    which would usually be the number of elements in the X and Y directions.
!    The value of N is
!
!      N = 3*NX*NY + 2*NX + 2*NY + 1,
!
!    The matrix is the consistent mass matrix for a regular NX by NY grid
!    of 8 node serendipity elements.
!
!    The local element numbering is
!
!      3--2--1
!      |     |
!      4     8
!      |     |
!      5--6--7
!
!    Here is an illustration for NX = 3, NY = 2:
!
!     23-24-25-26-27-28-29
!      |     |     |     |
!     19    20    21    22
!      |     |     |     |
!     12-13-14-15-16-17-18
!      |     |     |     |
!      8     9    10    11
!      |     |     |     |
!      1--2--3--4--5--6--7
!
!    For this example, the total number of nodes is, as expected,
!
!      N = 3 * 3 * 2 + 2 * 2 + 2 * 3 + 1 = 29
!
!    The matrix is symmetric positive definite for any positive values of the
!    density RHO(X,Y).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Nicholas Higham,
!    Algorithm 694: A Collection of Test Matrices in MATLAB,
!    ACM Transactions on Mathematical Software,
!    Volume 17, Number 3, September 1991, pages 289-305.
!
!    Andrew Wathen,
!    Realistic eigenvalue bounds for the Galerkin mass matrix,
!    IMA Journal of Numerical Analysis,
!    Volume 7, Number 4, October 1987, pages 449-457.
!
!  Parameters:
!
!    Input, integer NX, NY, values which determine the size
!    of the matrix.
!
!    Input, integer N, the number of rows and columns.
!
!    Output, real ( kind = rk ) A(9*NX+13,N), the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nx
  integer ny

  real ( kind = rk ) a(9*nx+13,n)
  real ( kind = rk ), dimension ( 8, 8 ), save :: em =  reshape ( (/ &
     6.0, -6.0,  2.0, -8.0,  3.0, -8.0,  2.0, -6.0, &
    -6.0, 32.0, -6.0, 20.0, -8.0, 16.0, -8.0, 20.0, &
     2.0, -6.0,  6.0, -6.0,  2.0, -8.0,  3.0, -8.0, &
    -8.0, 20.0, -6.0, 32.0, -6.0, 20.0, -8.0, 16.0, &
     3.0, -8.0,  2.0, -6.0,  6.0, -6.0,  2.0, -8.0, &
    -8.0, 16.0, -8.0, 20.0, -6.0, 32.0, -6.0, 20.0, &
     2.0, -8.0,  3.0, -8.0,  2.0, -6.0,  6.0, -6.0, &
    -6.0, 20.0, -8.0, 16.0, -8.0, 20.0, -6.0, 32.0 /), &
    (/ 8, 8 /) )
  integer i
  integer ii
  integer j
  integer jj
  integer kcol
  integer krow
  integer ml
  integer mu
  integer node(8)
  real ( kind = rk ) rho

  ml = 3 * nx + 4
  mu = 3 * nx + 4

  a(1:9*nx+13,1:n) = 0.0D+00

  do j = 1, ny
    do i = 1, nx

      node(1) = 3 * j * nx + 2 * j + 2 * i + 1
      node(2) = node(1) - 1
      node(3) = node(1) - 2
      node(4) = ( 3 * j - 1 ) * nx + 2 * j + i - 1
      node(5) = ( 3 * j - 3 ) * nx + 2 * j + 2 * i - 3
      node(6) = node(5) + 1
      node(7) = node(5) + 2
      node(8) = node(4) + 1

      call random_number ( harvest = rho )
      rho = 100.0D+00 * rho

      do krow = 1, 8
        do kcol = 1, 8
          ii = node(krow);
          jj = node(kcol);
          a(ii-jj+ml+mu+1,jj) = a(ii-jj+ml+mu+1,jj) &
            + rho * em(krow,kcol)
        end do
      end do

    end do
  end do

  return
end
subroutine wathen_ge ( nx, ny, n, a )

!*****************************************************************************80
!
!! WATHEN_GE returns the Wathen matrix as a general storage (GE) matrix.
!
!  Discussion:
!
!    The Wathen matrix is a finite element matrix which is sparse.
!
!    The entries of the matrix depend in part on a physical quantity
!    related to density.  That density is here assigned random values between
!    0 and 100.
!
!    The matrix order N is determined by the input quantities NX and NY,
!    which would usually be the number of elements in the X and Y directions.
!    The value of N is
!
!      N = 3*NX*NY + 2*NX + 2*NY + 1,
!
!    The matrix is the consistent mass matrix for a regular NX by NY grid
!    of 8 node serendipity elements.
!
!    The local element numbering is
!
!      3--2--1
!      |     |
!      4     8
!      |     |
!      5--6--7
!
!    Here is an illustration for NX = 3, NY = 2:
!
!     23-24-25-26-27-28-29
!      |     |     |     |
!     19    20    21    22
!      |     |     |     |
!     12-13-14-15-16-17-18
!      |     |     |     |
!      8     9    10    11
!      |     |     |     |
!      1--2--3--4--5--6--7
!
!    For this example, the total number of nodes is, as expected,
!
!      N = 3 * 3 * 2 + 2 * 2 + 2 * 3 + 1 = 29
!
!    The matrix is symmetric positive definite for any positive values of the
!    density RHO(X,Y).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Nicholas Higham,
!    Algorithm 694: A Collection of Test Matrices in MATLAB,
!    ACM Transactions on Mathematical Software,
!    Volume 17, Number 3, September 1991, pages 289-305.
!
!    Andrew Wathen,
!    Realistic eigenvalue bounds for the Galerkin mass matrix,
!    IMA Journal of Numerical Analysis,
!    Volume 7, Number 4, October 1987, pages 449-457.
!
!  Parameters:
!
!    Input, integer NX, NY, values which determine the size 
!    of the matrix.
!
!    Input, integer N, the number of rows and columns.
!
!    Output, real ( kind = rk ) A(N,N), the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nx
  integer ny

  real ( kind = rk ) a(n,n)
  real ( kind = rk ), dimension ( 8, 8 ), save :: em =  reshape ( (/ &
     6.0, -6.0,  2.0, -8.0,  3.0, -8.0,  2.0, -6.0, &
    -6.0, 32.0, -6.0, 20.0, -8.0, 16.0, -8.0, 20.0, &
     2.0, -6.0,  6.0, -6.0,  2.0, -8.0,  3.0, -8.0, &
    -8.0, 20.0, -6.0, 32.0, -6.0, 20.0, -8.0, 16.0, &
     3.0, -8.0,  2.0, -6.0,  6.0, -6.0,  2.0, -8.0, &
    -8.0, 16.0, -8.0, 20.0, -6.0, 32.0, -6.0, 20.0, &
     2.0, -8.0,  3.0, -8.0,  2.0, -6.0,  6.0, -6.0, &
    -6.0, 20.0, -8.0, 16.0, -8.0, 20.0, -6.0, 32.0 /), &
    (/ 8, 8 /) )
  integer i
  integer j
  integer kcol
  integer krow
  integer node(8)
  real ( kind = rk ) rho

  a(1:n,1:n) = 0.0D+00

  do j = 1, ny
    do i = 1, nx

      node(1) = 3 * j * nx + 2 * j + 2 * i + 1
      node(2) = node(1) - 1
      node(3) = node(1) - 2
      node(4) = ( 3 * j - 1 ) * nx + 2 * j + i - 1
      node(5) = ( 3 * j - 3 ) * nx + 2 * j + 2 * i - 3
      node(6) = node(5) + 1
      node(7) = node(5) + 2
      node(8) = node(4) + 1

      call random_number ( harvest = rho )
      rho = 100.0D+00 * rho

      do krow = 1, 8
        do kcol = 1, 8
          a(node(krow),node(kcol)) = a(node(krow),node(kcol)) &
            + rho * em(krow,kcol)
        end do
      end do

    end do
  end do

  return
end
subroutine wathen_order ( nx, ny, n )

!*****************************************************************************80
!
!! WATHEN_ORDER returns the order of the WATHEN matrix.
!
!  Discussion:
!
!    N = 3 * 3 * 2 + 2 * 2 + 2 * 3 + 1 = 29
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 January 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Nicholas Higham,
!    Algorithm 694: A Collection of Test Matrices in MATLAB,
!    ACM Transactions on Mathematical Software,
!    Volume 17, Number 3, September 1991, pages 289-305.
!
!    Andrew Wathen,
!    Realistic eigenvalue bounds for the Galerkin mass matrix,
!    IMA Journal of Numerical Analysis,
!    Volume 7, 1987, pages 449-457.
!
!  Parameters:
!
!    Input, integer NX, NY, values which determine the size of A.
!
!    Output, integer N, the order of the matrix,
!    as determined by NX and NY.
!
  implicit none

  integer n
  integer nx
  integer ny

  n = 3 * nx * ny + 2 * nx + 2 * ny + 1

  return
end
subroutine wathen_st ( nx, ny, nz_num, row, col, a )

!*****************************************************************************80
!
!! WATHEN_ST: Wathen matrix stored in sparse triplet (ST) format.
!
!  Discussion:
!
!    When dealing with sparse matrices in MATLAB, it can be much more efficient
!    to work first with a triple of I, J, and X vectors, and only once
!    they are complete, convert to MATLAB's sparse format.
!
!    The Wathen matrix is a finite element matrix which is sparse.
!
!    The entries of the matrix depend in part on a physical quantity
!    related to density.  That density is here assigned random values between
!    0 and 100.
!
!    The matrix order N is determined by the input quantities NX and NY,
!    which would usually be the number of elements in the X and Y directions.
!
!    The value of N is
!
!      N = 3*NX*NY + 2*NX + 2*NY + 1,
!
!    The matrix is the consistent mass matrix for a regular NX by NY grid
!    of 8 node serendipity elements.
!
!    The local element numbering is
!
!      3--2--1
!      |     |
!      4     8
!      |     |
!      5--6--7
!
!    Here is an illustration for NX = 3, NY = 2:
!
!     23-24-25-26-27-28-29
!      |     |     |     |
!     19    20    21    22
!      |     |     |     |
!     12-13-14-15-16-17-18
!      |     |     |     |
!      8     9    10    11
!      |     |     |     |
!      1--2--3--4--5--6--7
!
!    For this example, the total number of nodes is, as expected,
!
!      N = 3 * 3 * 2 + 2 * 2 + 2 * 3 + 1 = 29
!
!    The matrix is symmetric positive definite for any positive values of the
!    density RHO(X,Y).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 July 2014
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Nicholas Higham,
!    Algorithm 694: A Collection of Test Matrices in MATLAB,
!    ACM Transactions on Mathematical Software,
!    Volume 17, Number 3, September 1991, pages 289-305.
!
!    Andrew Wathen,
!    Realistic eigenvalue bounds for the Galerkin mass matrix,
!    IMA Journal of Numerical Analysis,
!    Volume 7, Number 4, October 1987, pages 449-457.
!
!  Parameters:
!
!    Input, integer NX, NY, values which determine the size of 
!    the matrix.
!
!    Input, integer NZ_NUM, the number of values used to 
!    describe the matrix.
!
!    Output, integer ROW(NZ_NUM), COL(NZ_NUM), the row and 
!    column indices of the nonzero entries.
!
!    Output, real ( kind = rk ) A(NZ_NUM), the nonzero entries of the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nx
  integer ny
  integer nz_num

  real ( kind = rk ) a(nz_num)
  integer col(nz_num)
  real ( kind = rk ), dimension ( 8, 8 ), save :: em =  reshape ( (/ &
     6.0, -6.0,  2.0, -8.0,  3.0, -8.0,  2.0, -6.0, &
    -6.0, 32.0, -6.0, 20.0, -8.0, 16.0, -8.0, 20.0, &
     2.0, -6.0,  6.0, -6.0,  2.0, -8.0,  3.0, -8.0, &
    -8.0, 20.0, -6.0, 32.0, -6.0, 20.0, -8.0, 16.0, &
     3.0, -8.0,  2.0, -6.0,  6.0, -6.0,  2.0, -8.0, &
    -8.0, 16.0, -8.0, 20.0, -6.0, 32.0, -6.0, 20.0, &
     2.0, -8.0,  3.0, -8.0,  2.0, -6.0,  6.0, -6.0, &
    -6.0, 20.0, -8.0, 16.0, -8.0, 20.0, -6.0, 32.0 /), &
    (/ 8, 8 /) )
  integer i
  integer j
  integer k
  integer kcol
  integer krow
  integer node(8)
  real ( kind = rk ) rho
  integer row(nz_num)

  row(1:nz_num) = 0
  col(1:nz_num) = 0
  a(1:nz_num) = 0.0D+00
  
  k = 0

  do j = 1, ny
    do i = 1, nx

      node(1) = 3 * j * nx + 2 * j + 2 * i + 1
      node(2) = node(1) - 1
      node(3) = node(1) - 2
      node(4) = ( 3 * j - 1 ) * nx + 2 * j + i - 1
      node(5) = ( 3 * j - 3 ) * nx + 2 * j + 2 * i - 3
      node(6) = node(5) + 1
      node(7) = node(5) + 2
      node(8) = node(4) + 1

      call random_number ( harvest = rho )
      rho = 100.0D+00 * rho

      do krow = 1, 8
        do kcol = 1, 8
          k = k + 1
          row(k) = node(krow)
          col(k) = node(kcol)
          a(k) = rho * em(krow,kcol)
        end do
      end do

    end do
  end do

  return
end
subroutine wathen_st_size ( nx, ny, nz_num )

!*****************************************************************************80
!
!! WATHEN_ST_SIZE: Size of Wathen matrix stored in sparse triplet format.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 June 2014
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Nicholas Higham,
!    Algorithm 694: A Collection of Test Matrices in MATLAB,
!    ACM Transactions on Mathematical Software,
!    Volume 17, Number 3, September 1991, pages 289-305.
!
!    Andrew Wathen,
!    Realistic eigenvalue bounds for the Galerkin mass matrix,
!    IMA Journal of Numerical Analysis,
!    Volume 7, Number 4, October 1987, pages 449-457.
!
!  Parameters:
!
!    Input, integer NX, NY, values which determine the size of the matrix.
!
!    Output, integer NZ_NUM, the number of items of data used to describe
!    the matrix.
!
  implicit none

  integer nx
  integer ny
  integer nz_num

  nz_num = nx * ny * 64

  return
end
subroutine wathen_xy ( nx, ny, n, x, y )

!*****************************************************************************80
!
!! wathen_xy: coordinates of Wathen nodes.
!
!  Discussion:
!
!    We will take the region to be the unit square.
!
!    The grid uses quadratic serendipity elements.
!
!    Here is an illustration of the node numbering for NX = 3, NY = 2:
!
!     23-24-25-26-27-28-29
!      |     |     |     |
!     19    20    21    22
!      |     |     |     |
!     12-13-14-15-16-17-18
!      |     |     |     |
!      8     9    10    11
!      |     |     |     |
!      1--2--3--4--5--6--7
!
!    For this example, the total number of nodes is, as expected,
!
!      N = 3 * 3 * 2 + 2 * 2 + 2 * 3 + 1 = 29
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 February 2020
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Nicholas Higham,
!    Algorithm 694: A Collection of Test Matrices in MATLAB,
!    ACM Transactions on Mathematical Software,
!    Volume 17, Number 3, September 1991, pages 289-305.
!
!    Andrew Wathen,
!    Realistic eigenvalue bounds for the Galerkin mass matrix,
!    IMA Journal of Numerical Analysis,
!    Volume 7, Number 4, October 1987, pages 449-457.
!
!  Input:
!
!    integer NX, NY: values which determine the size of the matrix.
!
!    integer N: the number of variables.
!
!  Output:
!
!    real ( kind = rk ) X(N), Y(N): the node coordinates.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer j
  integer k
  integer nx
  integer ny
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)

  k = 0
  do j = 1, 2 * ny + 1

    if ( mod ( j, 2 ) == 1 ) then
      do i = 1, 2 * nx + 1
        x(k+i) = real ( i - 1, kind = rk ) / real ( 2 * nx, kind = rk ) 
      end do
      y(k+1:k+2*nx+1) = real ( j - 1, kind = rk ) / real ( 2 * ny, kind = rk ) 
      k = k + 2 * nx + 1
    else
      do i = 1, nx + 1
        x(k+i) = real ( i - 1, kind = rk ) / real ( nx, kind = rk )
      end do
      y(k+1:k+nx+1) = real ( j - 1, kind = rk ) / real ( 2 * ny, kind = rk ) 
      k = k + nx + 1
    end if

  end do

  return
end

