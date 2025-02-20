subroutine daxpy ( n, da, dx, incx, dy, incy )

!*****************************************************************************80
!
!! daxpy() computes constant times a vector plus a vector.
!
!  Discussion:
!
!    Uses unrolled loops for increments equal to one.
!
!  Author:
!
!    Jack Dongarra
!
!  Reference:
!
!    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
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

  if ( da  == 0.0D+00 ) then
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

    do i = 1, m
      dy(i) = dy(i) + da * dx(i)
    end do

    do i = m+1, n, 4
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
!    This routine uses unrolled loops for increments equal to one.
!
!  Author:
!
!    Jack Dongarra
!
!  Reference:
!
!    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vectors.
!
!    Input, real ( kind = rk ) DX(*), the first vector.
!
!    Input, integer INCX, the increment between successive entries 
!    in X.
!
!    Input, real ( kind = rk ) DY(*), the second vector.
!
!    Input, integer INCY, the increment between successive entries 
!    in Y.
!
!    Output, real DDOT, the sum of the product of the corresponding
!    entries of X and Y.
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

    do i = m+1, n, 5

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
function dnrm2 ( n, x, incx )

!*****************************************************************************80
!
!! DNRM2 returns the euclidean norm of a vector.
!
!  Discussion:
!
!     DNRM2 ( X ) = sqrt ( X' * X )
!
!  Author:
!
!    Sven Hammarling
!
!  Reference:
!
!    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real ( kind = rk ) X(*), the vector whose norm is to be computed.
!
!    Input, integer INCX, the increment between successive entries 
!    of X.
!
!    Output, real ( kind = rk ) DNRM2, the Euclidean norm of X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) absxi
  real ( kind = rk ) dnrm2
  integer incx
  integer ix
  integer n
  real ( kind = rk ) norm
  real ( kind = rk ) scale
  real ( kind = rk ) ssq
  real ( kind = rk ) x(*)

  if ( n < 1 .or. incx < 1 ) then

    norm  = 0.0D+00

  else if ( n == 1 ) then

    norm  = abs ( x(1) )

  else

    scale = 0.0D+00
    ssq = 1.0D+00

    do ix = 1, 1 + ( n - 1 )*incx, incx
      if ( x(ix) /= 0.0D+00 ) then
        absxi = abs ( x(ix) )
        if ( scale < absxi ) then
          ssq = 1.0D+00 + ssq * ( scale / absxi )**2
          scale = absxi
        else
          ssq = ssq + ( absxi / scale )**2
        end if
      end if
    end do
    norm  = scale * sqrt( ssq )
  end if

  dnrm2 = norm

  return
end
subroutine drot ( n, x, incx, y, incy, c, s )

!*****************************************************************************80
!
!! DROT applies a plane rotation.
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
!    Input/output, real ( kind = rk ) X(*), one of the vectors to be rotated.
!
!    Input, integer INCX, the increment between successive 
!    entries of X.
!
!    Input/output, real ( kind = rk ) Y(*), one of the vectors to be rotated.
!
!    Input, integer INCY, the increment between successive
!    elements of Y.
!
!    Input, real ( kind = rk ) C, S, parameters (presumably the cosine and
!    sine of some angle) that define a plane rotation.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) c
  integer i
  integer incx
  integer incy
  integer ix
  integer iy
  integer n
  real ( kind = rk ) s
  real ( kind = rk ) stemp
  real ( kind = rk ) x(*)
  real ( kind = rk ) y(*)

  if ( n <= 0 ) then

  else if ( incx == 1 .and. incy == 1 ) then

    do i = 1, n
      stemp = c * x(i) + s * y(i)
      y(i) = c * y(i) - s * x(i)
      x(i) = stemp
    end do

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
      stemp = c * x(ix) + s * y(iy)
      y(iy) = c * y(iy) - s * x(ix)
      x(ix) = stemp
      ix = ix + incx
      iy = iy + incy
    end do

  end if

  return
end
subroutine drotg ( sa, sb, c, s )

!*****************************************************************************80
!
!! DROTG constructs a Givens plane rotation.
!
!  Discussion:
!
!    This routine uses double precision real arithmetic.
!
!    Given values A and B, this routine computes
!
!    SIGMA = sign ( A ) if abs ( A ) >  abs ( B )
!          = sign ( B ) if abs ( A ) <= abs ( B );
!
!    R     = SIGMA * ( A * A + B * B );
!
!    C = A / R if R is not 0
!      = 1     if R is 0;
!
!    S = B / R if R is not 0,
!        0     if R is 0.
!
!    The computed numbers then satisfy the equation
!
!    (  C  S ) ( A ) = ( R )
!    ( -S  C ) ( B ) = ( 0 )
!
!    The routine also computes
!
!    Z = S     if abs ( A ) > abs ( B ),
!      = 1 / C if abs ( A ) <= abs ( B ) and C is not 0,
!      = 1     if C is 0.
!
!    The single value Z encodes C and S, and hence the rotation:
!
!    If Z = 1, set C = 0 and S = 1;
!    If abs ( Z ) < 1, set C = sqrt ( 1 - Z * Z ) and S = Z;
!    if abs ( Z ) > 1, set C = 1/ Z and S = sqrt ( 1 - C * C );
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 May 2006
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
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Algorithm 539, 
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input/output, real ( kind = rk ) SA, SB.  On input, SA and SB are the values
!    A and B.  On output, SA is overwritten with R, and SB is
!    overwritten with Z.
!
!    Output, real ( kind = rk ) C, S, the cosine and sine of the
!    Givens rotation.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) c
  real ( kind = rk ) r
  real ( kind = rk ) roe
  real ( kind = rk ) s
  real ( kind = rk ) sa
  real ( kind = rk ) sb
  real ( kind = rk ) scale
  real ( kind = rk ) z

  if ( abs ( sb ) < abs ( sa ) ) then
    roe = sa
  else
    roe = sb
  end if

  scale = abs ( sa ) + abs ( sb )

  if ( scale == 0.0D+00 ) then
    c = 1.0D+00
    s = 0.0D+00
    r = 0.0D+00
  else
    r = scale * sqrt ( ( sa / scale )**2 + ( sb / scale )**2 )
    r = sign ( 1.0D+00, roe ) * r
    c = sa / r
    s = sb / r
  end if

  if ( 0.0D+00 < abs ( c ) .and. abs ( c ) <= s ) then
    z = 1.0D+00 / c
  else
    z = s
  end if

  sa = r
  sb = z

  return
end
subroutine dscal ( n, sa, x, incx )

!*****************************************************************************80
!
!! DSCAL scales a vector by a constant.
!
!  Modified:
!
!    08 April 1999
!
!  Author:
!
!    Jack Dongarra
!
!  Reference:
!
!    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
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
subroutine dsvdc ( a, lda, m, n, s, e, u, ldu, v, ldv, work, job, info )

!*****************************************************************************80
!
!! DSVDC computes the singular value decomposition of a real rectangular matrix.
!
!  Discussion:
!
!    This routine reduces an M by N matrix A to diagonal form by orthogonal
!    transformations U and V.  The diagonal elements S(I) are the singular
!    values of A.  The columns of U are the corresponding left singular
!    vectors, and the columns of V the right singular vectors.
!
!    The form of the singular value decomposition is then
!
!      A(MxN) = U(MxM) * S(MxN) * V(NxN)'
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 September 2006
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
!    Input/output, real ( kind = rk ) A(LDA,N).  On input, the M by N
!    matrix whose singular value decomposition is to be computed.
!    On output, the matrix has been destroyed.  Depending on the user's
!    requests, the matrix may contain other useful information.
!
!    Input, integer LDA, the leading dimension of the array A.
!    LDA must be at least N.
!
!    Input, integer M, the number of rows of the matrix.
!
!    Input, integer N, the number of columns of the matrix A.
!
!    Output, real ( kind = rk ) S(MM), where MM = max(M+1,N).  The first
!    min(M,N) entries of S contain the singular values of A arranged in
!    descending order of magnitude.
!
!    Output, real ( kind = rk ) E(MM), where MM = max(M+1,N).  Ordinarily
!    contains zeros.  However see the discussion of INFO for exceptions.
!
!    Output, real ( kind = rk ) U(LDU,K).  If JOBA = 1 then K = M;
!    if 2 <= JOBA, then K = min(M,N).  U contains the M by M matrix of
!    left singular vectors.  U is not referenced if JOBA = 0.  If M <= N
!    or if JOBA = 2, then U may be identified with A in the subroutine call.
!
!    Input, integer LDU, the leading dimension of the array U.
!    LDU must be at least M.
!
!    Output, real ( kind = rk ) V(LDV,N), the N by N matrix of right singular
!    vectors.  V is not referenced if JOB is 0.  If N <= M, then V may be
!    identified with A in the subroutine call.
!
!    Input, integer LDV, the leading dimension of the array V.
!    LDV must be at least N.
!
!    Workspace, real ( kind = rk ) WORK(M).
!
!    Input, integer JOB, controls the computation of the singular
!    vectors.  It has the decimal expansion AB with the following meaning:
!      A =  0, do not compute the left singular vectors.
!      A =  1, return the M left singular vectors in U.
!      A >= 2, return the first min(M,N) singular vectors in U.
!      B =  0, do not compute the right singular vectors.
!      B =  1, return the right singular vectors in V.
!
!    Output, integer INFO, status indicator.
!    The singular values (and their corresponding singular vectors)
!    S(INFO+1), S(INFO+2),...,S(MN) are correct.  Here MN = min ( M, N ).
!    Thus if INFO is 0, all the singular values and their vectors are
!    correct.  In any event, the matrix B = U' * A * V is the bidiagonal
!    matrix with the elements of S on its diagonal and the elements of E on
!    its superdiagonal.  Thus the singular values of A and B are the same.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer lda
  integer ldu
  integer ldv
  integer m
  integer n

  real ( kind = rk ) a(lda,n)
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) cs
  real ( kind = rk ) e(*)
  real ( kind = rk ) el
  real ( kind = rk ) emm1
  real ( kind = rk ) f
  real ( kind = rk ) g
  integer info
  integer iter
  integer j
  integer job
  integer jobu
  integer k
  integer kase
  integer kk
  integer l
  integer ll
  integer lls
  integer ls
  integer lu
  integer, parameter :: maxit = 30
  integer mm
  integer mm1
  integer mn
  integer nct
  integer nctp1
  integer ncu
  integer nrt
  integer nrtp1
  real ( kind = rk ) s(*)
  real ( kind = rk ) scale
  real ( kind = rk ) ddot
  real ( kind = rk ) shift
  real ( kind = rk ) sl
  real ( kind = rk ) sm
  real ( kind = rk ) smm1
  real ( kind = rk ) sn
  real ( kind = rk ) dnrm2
  real ( kind = rk ) t
  real ( kind = rk ) t1
  real ( kind = rk ) test
  real ( kind = rk ) u(ldu,m)
  real ( kind = rk ) v(ldv,n)
  logical wantu
  logical wantv
  real ( kind = rk ) work(m)
  real ( kind = rk ) ztest
!
!  Determine what is to be computed.
!
  wantu = .false.
  wantv = .false.
  jobu = mod ( job, 100 ) / 10

  if ( 1 < jobu ) then
    ncu = min ( m, n )
  else
    ncu = m
  end if

  if ( jobu /= 0 ) then
    wantu = .true.
  end if

  if ( mod ( job, 10 ) /= 0 ) then
    wantv = .true.
  end if
!
!  Reduce A to bidiagonal form, storing the diagonal elements
!  in S and the super-diagonal elements in E.
!
  info = 0
  nct = min ( m-1, n )
  nrt = max ( 0, min ( m, n-2 ) )
  lu = max ( nct, nrt )

  do l = 1, lu
!
!  Compute the transformation for the L-th column and
!  place the L-th diagonal in S(L).
!
    if ( l <= nct ) then

      s(l) = dnrm2 ( m-l+1, a(l,l), 1 )

      if ( s(l) /= 0.0D+00 ) then
        if ( a(l,l) /= 0.0D+00 ) then
          s(l) = sign ( s(l), a(l,l) )
        end if
        call dscal ( m-l+1, 1.0D+00 / s(l), a(l,l), 1 )
        a(l,l) = 1.0D+00 + a(l,l)
      end if

      s(l) = -s(l)

    end if

    do j = l+1, n
!
!  Apply the transformation.
!
      if ( l <= nct .and. s(l) /= 0.0D+00 ) then
        t = -ddot ( m-l+1, a(l,l), 1, a(l,j), 1 ) / a(l,l)
        call daxpy ( m-l+1, t, a(l,l), 1, a(l,j), 1 )
      end if
!
!  Place the L-th row of A into E for the
!  subsequent calculation of the row transformation.
!
      e(j) = a(l,j)

    end do
!
!  Place the transformation in U for subsequent back multiplication.
!
    if ( wantu .and. l <= nct ) then
      u(l:m,l) = a(l:m,l)
    end if
!
!  Compute the L-th row transformation and place the
!  L-th superdiagonal in E(L).
!
    if ( l <= nrt ) then

      e(l) = dnrm2 ( n-l, e(l+1), 1 )

      if ( e(l) /= 0.0D+00 ) then
        if ( e(l+1) /= 0.0D+00 ) then
          e(l) = sign ( e(l), e(l+1) )
        end if
        call dscal ( n-l, 1.0D+00 / e(l), e(l+1), 1 )
        e(l+1) = 1.0D+00 + e(l+1)
      end if

      e(l) = -e(l)
!
!  Apply the transformation.
!
      if ( l + 1 <= m .and. e(l) /= 0.0D+00 ) then

        work(l+1:m) = 0.0D+00

        do j = l+1, n
          call daxpy ( m-l, e(j), a(l+1,j), 1, work(l+1), 1 )
        end do

        do j = l+1, n
          call daxpy ( m-l, -e(j)/e(l+1), work(l+1), 1, a(l+1,j), 1 )
        end do

      end if
!
!  Place the transformation in V for subsequent back multiplication.
!
      if ( wantv ) then
        v(l+1:n,l) = e(l+1:n)
      end if

    end if

  end do
!
!  Set up the final bidiagonal matrix of order MN.
!
  mn = min ( m + 1, n )
  nctp1 = nct + 1
  nrtp1 = nrt + 1

  if ( nct < n ) then
    s(nctp1) = a(nctp1,nctp1)
  end if

  if ( m < mn ) then
    s(mn) = 0.0D+00
  end if

  if ( nrtp1 < mn ) then
    e(nrtp1) = a(nrtp1,mn)
  end if

  e(mn) = 0.0D+00
!
!  If required, generate U.
!
  if ( wantu ) then

    u(1:m,nctp1:ncu) = 0.0D+00

    do j = nctp1, ncu
      u(j,j) = 1.0D+00
    end do

    do ll = 1, nct

      l = nct - ll + 1

      if ( s(l) /= 0.0D+00 ) then

        do j = l+1, ncu
          t = -ddot ( m-l+1, u(l,l), 1, u(l,j), 1 ) / u(l,l)
          call daxpy ( m-l+1, t, u(l,l), 1, u(l,j), 1 )
        end do

        u(l:m,l) = -u(l:m,l)
        u(l,l) = 1.0D+00 + u(l,l)
        u(1:l-1,l) = 0.0D+00

      else

        u(1:m,l) = 0.0D+00
        u(l,l) = 1.0D+00

      end if

    end do

  end if
!
!  If it is required, generate V.
!
  if ( wantv ) then

    do ll = 1, n

      l = n - ll + 1

      if ( l <= nrt .and. e(l) /= 0.0D+00 ) then

        do j = l + 1, n
          t = -ddot ( n-l, v(l+1,l), 1, v(l+1,j), 1 ) / v(l+1,l)
          call daxpy ( n-l, t, v(l+1,l), 1, v(l+1,j), 1 )
        end do

      end if

      v(1:n,l) = 0.0D+00
      v(l,l) = 1.0D+00

    end do

  end if
!
!  Main iteration loop for the singular values.
!
  mm = mn
  iter = 0

  do while ( 0 < mn )
!
!  If too many iterations have been performed, set flag and return.
!
    if ( maxit <= iter ) then
      info = mn
      return
    end if
!
!  This section of the program inspects for
!  negligible elements in the S and E arrays.
!
!  On completion the variables KASE and L are set as follows:
!
!  KASE = 1     if S(MN) and E(L-1) are negligible and L < MN
!  KASE = 2     if S(L) is negligible and L < MN
!  KASE = 3     if E(L-1) is negligible, L < MN, and
!               S(L), ..., S(MN) are not negligible (QR step).
!  KASE = 4     if E(MN-1) is negligible (convergence).
!
    do ll = 1, mn

      l = mn - ll

      if ( l == 0 ) then
        exit
      end if

      test = abs ( s(l) ) + abs ( s(l+1) )
      ztest = test + abs ( e(l) )

      if ( ztest == test ) then
        e(l) = 0.0D+00
        exit
      end if

    end do

    if ( l == mn - 1 ) then

      kase = 4

    else

      do lls = l + 1, mn + 1

        ls = mn - lls + l + 1

        if ( ls == l ) then
          exit
        end if

        test = 0.0D+00
        if ( ls /= mn ) then
          test = test + abs ( e(ls) )
        end if

        if ( ls /= l + 1 ) then
          test = test + abs ( e(ls-1) )
        end if

        ztest = test + abs ( s(ls) )

        if ( ztest == test ) then
          s(ls) = 0.0D+00
          exit
        end if

      end do

      if ( ls == l ) then
        kase = 3
      else if ( ls == mn ) then
        kase = 1
      else
        kase = 2
        l = ls
      end if

    end if

    l = l + 1
!
!  Deflate negligible S(MN).
!
    if ( kase == 1 ) then

      mm1 = mn - 1
      f = e(mn-1)
      e(mn-1) = 0.0D+00

      do kk = l, mm1

        k = mm1 - kk + l
        t1 = s(k)
        call drotg ( t1, f, cs, sn )
        s(k) = t1

        if ( k /= l ) then
          f = -sn * e(k-1)
          e(k-1) = cs * e(k-1)
        end if

        if ( wantv ) then
          call drot ( n, v(1,k), 1, v(1,mn), 1, cs, sn )
        end if

      end do
!
!  Split at negligible S(L).
!
    else if ( kase == 2 ) then

      f = e(l-1)
      e(l-1) = 0.0D+00

      do k = l, mn

        t1 = s(k)
        call drotg ( t1, f, cs, sn )
        s(k) = t1
        f = -sn * e(k)
        e(k) = cs * e(k)
        if ( wantu ) then
          call drot ( m, u(1,k), 1, u(1,l-1), 1, cs, sn )
        end if

      end do
!
!  Perform one QR step.
!
    else if ( kase == 3 ) then
!
!  Calculate the shift.
!
      scale = max ( abs ( s(mn) ), abs ( s(mn-1) ), abs ( e(mn-1) ), &
                    abs ( s(l) ), abs ( e(l) ) )

      sm = s(mn) / scale
      smm1 = s(mn-1) / scale
      emm1 = e(mn-1) / scale
      sl = s(l) / scale
      el = e(l) / scale
      b = ( ( smm1 + sm ) * ( smm1 - sm ) + emm1 * emm1 ) / 2.0D+00
      c = sm  * sm * emm1 * emm1
      shift = 0.0D+00

      if ( b /= 0.0D+00 .or. c /= 0.0D+00 ) then
        shift = sqrt ( b * b + c )
        if ( b < 0.0D+00 ) then
          shift = -shift
        end if
        shift = c / ( b + shift )
      end if

      f = ( sl + sm ) * ( sl - sm ) + shift
      g = sl * el
!
!  Chase zeros.
!
      mm1 = mn - 1

      do k = l, mm1

        call drotg ( f, g, cs, sn )

        if ( k /= l ) then
          e(k-1) = f
        end if

        f = cs * s(k) + sn * e(k)
        e(k) = cs * e(k) - sn * s(k)
        g = sn * s(k+1)
        s(k+1) = cs * s(k+1)

        if ( wantv ) then
          call drot ( n, v(1,k), 1, v(1,k+1), 1, cs, sn )
        end if

        call drotg ( f, g, cs, sn )
        s(k) = f
        f = cs * e(k) + sn * s(k+1)
        s(k+1) = -sn * e(k) + cs * s(k+1)
        g = sn * e(k+1)
        e(k+1) = cs * e(k+1)

        if ( wantu .and. k < m ) then
          call drot ( m, u(1,k), 1, u(1,k+1), 1, cs, sn )
        end if

      end do

      e(mn-1) = f
      iter = iter + 1
!
!  Convergence.
!
    else if ( kase == 4 ) then
!
!  Make the singular value nonnegative.
!
      if ( s(l) < 0.0D+00 ) then
        s(l) = -s(l)
        if ( wantv ) then
          v(1:n,l) = -v(1:n,l)
        end if
      end if
!
!  Order the singular value.
!
      do

        if ( l == mm ) then
          exit
        end if

        if ( s(l+1) <= s(l) ) then
          exit
        end if

        t = s(l)
        s(l) = s(l+1)
        s(l+1) = t

        if ( wantv .and. l < n ) then
          call dswap ( n, v(1,l), 1, v(1,l+1), 1 )
        end if

        if ( wantu .and. l < m ) then
          call dswap ( m, u(1,l), 1, u(1,l+1), 1 )
        end if

        l = l + 1

      end do

      iter = 0
      mn = mn - 1

    end if

  end do

  return
end
subroutine dswap ( n, x, incx, y, incy )

!*****************************************************************************80
!
!! DSWAP interchanges two vectors.
!
!  Modified:
!
!    08 April 1999
!
!  Reference:
!
!    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vectors.
!
!    Input/output, real ( kind = rk ) X(*), one of the vectors to swap.
!
!    Input, integer INCX, the increment between successive 
!    entries of X.
!
!    Input/output, real ( kind = rk ) Y(*), one of the vectors to swap.
!
!    Input, integer INCY, the increment between successive 
!    elements of Y.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer incx
  integer incy
  integer ix
  integer iy
  integer m
  integer n
  real ( kind = rk ) temp
  real ( kind = rk ) x(*)
  real ( kind = rk ) y(*)

  if ( n <= 0 ) then

  else if ( incx == 1 .and. incy == 1 ) then

    m = mod ( n, 3 )

    do i = 1, m
      temp = x(i)
      x(i) = y(i)
      y(i) = temp
    end do

    do i = m + 1, n, 3

      temp = x(i)
      x(i) = y(i)
      y(i) = temp

      temp = x(i+1)
      x(i+1) = y(i+1)
      y(i+1) = temp

      temp = x(i+2)
      x(i+2) = y(i+2)
      y(i+2) = temp

    end do

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
      temp = x(ix)
      x(ix) = y(iy)
      y(iy) = temp
      ix = ix + incx
      iy = iy + incy
    end do

  end if

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IUNIT, the free unit number.
!
  implicit none

  integer i
  integer ios
  integer iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
function i4_uniform_ab ( a, b )

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
!    Output, integer I4_UNIFORM_AB, a number between A and B.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer a
  integer b
  integer i4_uniform_ab
  real r
  integer value

  call random_number ( harvest = r )
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
function i4vec_max ( n, a )

!*****************************************************************************80
!
!! I4VEC_MAX computes the maximum element of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input, integer A(N), the array.
!
!    Output, integer I4VEC_MAX, the value of the largest entry.
!
  implicit none

  integer n

  integer a(n)
  integer i4vec_max

  i4vec_max = maxval ( a(1:n) )

  return
end
subroutine imtqlx ( n, d, e, z )

!*****************************************************************************80
!
!! IMTQLX diagonalizes a symmetric tridiagonal matrix.
!
!  Discussion:
!
!    This routine is a slightly modified version of the EISPACK routine to 
!    perform the implicit QL algorithm on a symmetric tridiagonal matrix. 
!
!    The authors thank the authors of EISPACK for permission to use this
!    routine. 
!
!    It has been modified to produce the product Q' * Z, where Z is an input 
!    vector and Q is the orthogonal matrix diagonalizing the input matrix.  
!    The changes consist (essentially) of applying the orthogonal 
!    transformations directly to Z as they are generated.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 December 2009
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!    Roger Martin, James Wilkinson,
!    The Implicit QL Algorithm,
!    Numerische Mathematik,
!    Volume 12, Number 5, December 1968, pages 377-383.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real ( kind = rk ) D(N), the diagonal entries of the matrix.
!    On output, the information in D has been overwritten.
!
!    Input/output, real ( kind = rk ) E(N), the subdiagonal entries of the 
!    matrix, in entries E(1) through E(N-1).  On output, the information in
!    E has been overwritten.
!
!    Input/output, real ( kind = rk ) Z(N).  On input, a vector.  On output,
!    the value of Q' * Z, where Q is the matrix that diagonalizes the
!    input symmetric tridiagonal matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) d(n)
  real ( kind = rk ) e(n)
  real ( kind = rk ) f
  real ( kind = rk ) g
  integer i
  integer ii
  integer, parameter :: itn = 30
  integer j
  integer k
  integer l
  integer m
  integer mml
  real ( kind = rk ) p
  real ( kind = rk ) prec
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) z(n)

  prec = epsilon ( prec )

  if ( n == 1 ) then
    return
  end if

  e(n) = 0.0D+00

  do l = 1, n

    j = 0

    do

      do m = l, n

        if ( m == n ) then
          exit
        end if

        if ( abs ( e(m) ) <= prec * ( abs ( d(m) ) + abs ( d(m+1) ) ) ) then
          exit
        end if

      end do

      p = d(l)

      if ( m == l ) then
        exit
      end if

      if ( itn <= j ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IMTQLX - Fatal error!'
        write ( *, '(a)' ) '  Iteration limit exceeded.'
        write ( *, '(a,i8)' ) '  J = ', j
        write ( *, '(a,i8)' ) '  L = ', l
        write ( *, '(a,i8)' ) '  M = ', m
        write ( *, '(a,i8)' ) '  N = ', n
        stop 1
      end if

      j = j + 1
      g = ( d(l+1) - p ) / ( 2.0D+00 * e(l) )
      r =  sqrt ( g * g + 1.0D+00 )
      g = d(m) - p + e(l) / ( g + sign ( r, g ) )
      s = 1.0D+00
      c = 1.0D+00
      p = 0.0D+00
      mml = m - l

      do ii = 1, mml

        i = m - ii
        f = s * e(i)
        b = c * e(i)

        if ( abs ( g ) <= abs ( f ) ) then
          c = g / f
          r =  sqrt ( c * c + 1.0D+00 )
          e(i+1) = f * r
          s = 1.0D+00 / r
          c = c * s
        else
          s = f / g
          r =  sqrt ( s * s + 1.0D+00 )
          e(i+1) = g * r
          c = 1.0D+00 / r
          s = s * c
        end if

        g = d(i+1) - p
        r = ( d(i) - g ) * s + 2.0D+00 * c * b
        p = s * r
        d(i+1) = g + p
        g = c * r - b
        f = z(i+1)
        z(i+1) = s * z(i) + c * f
        z(i) = c * z(i) - s * f

      end do

      d(l) = d(l) - p
      e(l) = g
      e(m) = 0.0D+00

    end do

  end do
!
!  Sorting.
!
  do ii = 2, n

    i = ii - 1
    k = i
    p = d(i)

    do j = ii, n
      if ( d(j) < p ) then
        k = j
        p = d(j)
      end if
    end do

    if ( k /= i ) then
      d(k) = d(i)
      d(i) = p
      p = z(i)
      z(i) = z(k)
      z(k) = p
    end if

  end do

  return
end
function r8_choose ( n, k )

!*****************************************************************************80
!
!! R8_CHOOSE computes the binomial coefficient C(N,K) as an R8.
!
!  Discussion:
!
!    The value is calculated in such a way as to avoid overflow and
!    roundoff.  The calculation is done in R8 arithmetic.
!
!    The formula used is:
!
!      C(N,K) = N! / ( K! * (N-K)! )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    ML Wolfson, HV Wright,
!    Algorithm 160:
!    Combinatorial of M Things Taken N at a Time,
!    Communications of the ACM,
!    Volume 6, Number 4, April 1963, page 161.
!
!  Parameters:
!
!    Input, integer N, K, are the values of N and K.
!
!    Output, real ( kind = rk ) R8_CHOOSE, the number of combinations of N
!    things taken K at a time.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer k
  integer mn
  integer mx
  integer n
  real ( kind = rk ) r8_choose
  real ( kind = rk ) value

  mn = min ( k, n - k )

  if ( mn < 0 ) then

    value = 0.0D+00

  else if ( mn == 0 ) then

    value = 1.0D+00

  else

    mx = max ( k, n - k )
    value = real ( mx + 1, kind = rk )

    do i = 2, mn
      value = ( value * real ( mx + i, kind = rk ) ) / real ( i, kind = rk )
    end do

  end if

  r8_choose = value

  return
end
function r8_factorial ( n )

!*****************************************************************************80
!
!! R8_FACTORIAL computes the factorial of N.
!
!  Discussion:
!
!    factorial ( N ) = product ( 1 <= I <= N ) I
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the argument of the factorial function.
!    If N is less than 1, the function value is returned as 1.
!
!    Output, real ( kind = rk ) R8_FACTORIAL, the factorial of N.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r8_factorial
  integer i
  integer n
  real ( kind = rk ) value

  value = 1.0D+00

  do i = 1, n
    value = value * real ( i, kind = rk )
  end do

  r8_factorial = value

  return
end
function r8_hyper_2f1 ( a_input, b_input, c_input, x_input )

!*****************************************************************************80
!
!! R8_HYPER_2F1 evaluates the hypergeometric function F(A,B,C,X).
!
!  Discussion:
!
!    A minor bug was corrected.  The HW variable, used in several places as
!    the "old" value of a quantity being iteratively improved, was not
!    being initialized.  JVB, 11 February 2008.
!
!    The original version of this program allowed the input arguments to
!    be modified, although they were restored to their input values before exit.
!    This is unacceptable if the input arguments are allowed to be constants.
!    The code has been modified so that the input arguments are never modified.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 October 2008
!
!  Author:
!
!    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
!    FORTRAN90 version by John Burkardt.
!
!    The original FORTRAN77 version of this routine is copyrighted by
!    Shanjie Zhang and Jianming Jin.  However, they give permission to
!    incorporate this routine into a user program provided that the copyright
!    is acknowledged.
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45
!
!  Parameters:
!
!    Input, real ( kind = rk ) A_INPUT, B_INPUT, C_INPUT, X_INPUT, 
!    the arguments of the function.  The user is allowed to pass these
!    values as constants or variables.
!    C_INPUT must not be equal to a nonpositive integer.
!    X_INPUT < 1.
!
!    Output, real ( kind = rk ) R8_HYPER_2F1, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) a_input
  real ( kind = rk ) a0
  real ( kind = rk ) aa
  real ( kind = rk ) b
  real ( kind = rk ) b_input
  real ( kind = rk ) bb
  real ( kind = rk ) c
  real ( kind = rk ) c_input
  real ( kind = rk ) c0
  real ( kind = rk ) c1
  real ( kind = rk ), parameter :: el = 0.5772156649015329D+00
  real ( kind = rk ) eps
  real ( kind = rk ) f0
  real ( kind = rk ) f1
  real ( kind = rk ) g0
  real ( kind = rk ) g1
  real ( kind = rk ) g2
  real ( kind = rk ) g3
  real ( kind = rk ) ga
  real ( kind = rk ) gabc
  real ( kind = rk ) gam
  real ( kind = rk ) gb
  real ( kind = rk ) gbm
  real ( kind = rk ) gc
  real ( kind = rk ) gca
  real ( kind = rk ) gcab
  real ( kind = rk ) gcb
  real ( kind = rk ) gm
  real ( kind = rk ) hw
  integer j
  integer k
  logical l0
  logical l1
  logical l2
  logical l3
  logical l4
  logical l5
  integer m
  integer nm
  real ( kind = rk ) pa
  real ( kind = rk ) pb
  real ( kind = rk ) r
  real ( kind = rk ) r0
  real ( kind = rk ) r1
  real ( kind = rk ) r8_hyper_2f1
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) r8_psi
  real ( kind = rk ) rm
  real ( kind = rk ) rp
  real ( kind = rk ) sm
  real ( kind = rk ) sp
  real ( kind = rk ) sp0
  real ( kind = rk ) value
  real ( kind = rk ) x
  real ( kind = rk ) x_input
  real ( kind = rk ) x1
!
!  Immediately copy the input arguments!
!
  a = a_input
  b = b_input
  c = c_input
  x = x_input

  l0 = ( c == aint ( c ) ) .and. ( c < 0.0D+00 )
  l1 = ( 1.0D+00 - x < 1.0D-15 ) .and. ( c - a - b <= 0.0D+00 )
  l2 = ( a == aint ( a ) ) .and. ( a < 0.0D+00 )
  l3 = ( b == aint ( b ) ) .and. ( b < 0.0D+00 )
  l4 = ( c - a == aint ( c - a ) ) .and. ( c - a <= 0.0D+00 )
  l5 = ( c - b == aint ( c - b ) ) .and. ( c - b <= 0.0D+00 )

  if ( l0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_HYPER_2F1 - Fatal error!'
    write ( *, '(a)' ) '  Integral C < 0.'
    write ( *, '(a)' ) '  The hypergeometric series is divergent.'
    stop 1
  end if

  if ( l1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_HYPER_2F1 - Fatal error!'
    write ( *, '(a)' ) '  The hypergeometric series is divergent.'
    write ( *, '(a)' ) '  1 - X < 0, C - A - B < 0.'
    stop 1
  end if

  if ( 0.95D+00 < x ) then
    eps = 1.0D-08
  else
    eps = 1.0D-15
  end if

  if ( x == 0.0D+00 .or. a == 0.0D+00 .or. b == 0.0D+00 ) then

    value = 1.0D+00
    r8_hyper_2f1 = value
    return

  else if ( 1.0D+00 - x == eps .and. 0.0D+00 < c - a - b ) then

    gc = gamma ( c )
    gcab = gamma ( c - a - b )
    gca = gamma ( c - a )
    gcb = gamma ( c - b )
    value = gc * gcab / ( gca * gcb )
    r8_hyper_2f1 = value
    return

  else if ( 1.0D+00 + x <= eps .and. abs ( c - a + b - 1.0D+00 ) <= eps ) then

    g0 = sqrt ( r8_pi ) * 2.0D+00**( - a )
    g1 = gamma ( c )
    g2 = gamma ( 1.0D+00 + a / 2.0D+00 - b )
    g3 = gamma ( 0.5D+00 + 0.5D+00 * a )
    value = g0 * g1 / ( g2 * g3 )
    r8_hyper_2f1 = value
    return

  else if ( l2 .or. l3 ) then

    if ( l2 ) then
      nm = int ( abs ( a ) )
    end if

    if ( l3 ) then
      nm = int ( abs ( b ) )
    end if

    value = 1.0D+00
    r = 1.0D+00

    do k = 1, nm
      r = r * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
        / ( k * ( c + k - 1.0D+00 ) ) * x
      value = value + r
    end do

    r8_hyper_2f1 = value
    return

  else if ( l4 .or. l5 ) then

    if ( l4 ) then
      nm = int ( abs ( c - a ) )
    end if

    if ( l5 ) then
      nm = int ( abs ( c - b ) )
    end if

    value = 1.0D+00
    r  = 1.0D+00
    do k = 1, nm
      r = r * ( c - a + k - 1.0D+00 ) * ( c - b + k - 1.0D+00 ) &
        / ( k * ( c + k - 1.0D+00 ) ) * x
      value = value + r
    end do
    value = ( 1.0D+00 - x )**( c - a - b ) * value
    r8_hyper_2f1 = value
    return

  end if

  aa = a
  bb = b
  x1 = x

  if ( x < 0.0D+00 ) then
    x = x / ( x - 1.0D+00 )
    if ( a < c .and. b < a .and. 0.0D+00 < b ) then
      a = bb
      b = aa
    end if
    b = c - b
  end if

  if ( 0.75D+00 <= x ) then

    gm = 0.0D+00

    if ( abs ( c - a - b - aint ( c - a - b ) ) < 1.0D-15 ) then

      m = int ( c - a - b )
      ga = gamma ( a )
      gb = gamma ( b )
      gc = gamma ( c )
      gam = gamma ( a + m )
      gbm = gamma ( b + m )

      pa = r8_psi ( a )
      pb = r8_psi ( b )

      if ( m /= 0 ) then
        gm = 1.0D+00
      end if

      do j = 1, abs ( m ) - 1
        gm = gm * j
      end do

      rm = 1.0D+00
      do j = 1, abs ( m )
        rm = rm * j
      end do

      f0 = 1.0D+00
      r0 = 1.0D+00
      r1 = 1.0D+00
      sp0 = 0.0D+00
      sp = 0.0D+00

      if ( 0 <= m ) then

        c0 = gm * gc / ( gam * gbm )
        c1 = - gc * ( x - 1.0D+00 ) ** m / ( ga * gb * rm )

        do k = 1, m - 1
          r0 = r0 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
            / ( k * ( k - m ) ) * ( 1.0D+00 - x )
          f0 = f0 + r0
        end do

        do k = 1, m
          sp0 = sp0 + 1.0D+00 / ( a + k - 1.0D+00 ) &
            + 1.0D+00 / ( b + k - 1.0D+00 ) - 1.0D+00 / real ( k, kind = rk )
        end do

        f1 = pa + pb + sp0 + 2.0D+00 * el + log ( 1.0D+00 - x )
        hw = f1

        do k = 1, 250

          sp = sp + ( 1.0D+00 - a ) / ( k * ( a + k - 1.0D+00 ) ) &
            + ( 1.0D+00 - b ) / ( k * ( b + k - 1.0D+00 ) )

          sm = 0.0D+00
          do j = 1, m
            sm = sm + ( 1.0D+00 - a ) &
              / ( ( j + k ) * ( a + j + k - 1.0D+00 ) ) &
              + 1.0D+00 / ( b + j + k - 1.0D+00 )
          end do

          rp = pa + pb + 2.0D+00 * el + sp + sm + log ( 1.0D+00 - x )

          r1 = r1 * ( a + m + k - 1.0D+00 ) * ( b + m + k - 1.0D+00 ) &
            / ( k * ( m + k ) ) * ( 1.0D+00 - x )

          f1 = f1 + r1 * rp

          if ( abs ( f1 - hw ) < abs ( f1 ) * eps ) then
            exit
          end if

          hw = f1

        end do

        value = f0 * c0 + f1 * c1

      else if ( m < 0 ) then

        m = - m
        c0 = gm * gc / ( ga * gb * ( 1.0D+00 - x )**m )
        c1 = - ( - 1 )**m * gc / ( gam * gbm * rm )

        do k = 1, m - 1
          r0 = r0 * ( a - m + k - 1.0D+00 ) * ( b - m + k - 1.0D+00 ) &
            / ( k * ( k - m ) ) * ( 1.0D+00 - x )
          f0 = f0 + r0
        end do

        do k = 1, m
          sp0 = sp0 + 1.0D+00 / real ( k, kind = rk )
        end do

        f1 = pa + pb - sp0 + 2.0D+00 * el + log ( 1.0D+00 - x )
        hw = f1

        do k = 1, 250

          sp = sp + ( 1.0D+00 - a ) &
            / ( k * ( a + k - 1.0D+00 ) ) &
            + ( 1.0D+00 - b ) / ( k * ( b + k - 1.0D+00 ) )

          sm = 0.0D+00
          do j = 1, m
            sm = sm + 1.0D+00 / real ( j + k, kind = rk )
          end do

          rp = pa + pb + 2.0D+00 * el + sp - sm + log ( 1.0D+00 - x )

          r1 = r1 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
            / ( k * ( m + k ) ) * ( 1.0D+00 - x )

          f1 = f1 + r1 * rp

          if ( abs ( f1 - hw ) < abs ( f1 ) * eps ) then
            exit
          end if

          hw = f1

        end do

        value = f0 * c0 + f1 * c1

      end if

    else

      ga = gamma ( a )
      gb = gamma ( b )
      gc = gamma ( c )
      gca = gamma ( c - a )
      gcb = gamma ( c - b )
      gcab = gamma ( c - a - b )
      gabc = gamma ( a + b - c )
      c0 = gc * gcab / ( gca * gcb )
      c1 = gc * gabc / ( ga * gb ) * ( 1.0D+00 - x )**( c - a - b )
      value = 0.0D+00
      hw = value
      r0 = c0
      r1 = c1

      do k = 1, 250

        r0 = r0 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
          / ( k * ( a + b - c + k ) ) * ( 1.0D+00 - x )

        r1 = r1 * ( c - a + k - 1.0D+00 ) * ( c - b + k - 1.0D+00 ) &
          / ( k * ( c - a - b + k ) ) * ( 1.0D+00 - x )

        value = value + r0 + r1

        if ( abs ( value - hw ) < abs ( value ) * eps ) then
          exit
        end if

        hw = value

      end do

      value = value + c0 + c1

    end if

  else

    a0 = 1.0D+00

    if ( a < c .and. c < 2.0D+00 * a .and. b < c .and. c < 2.0D+00 * b ) then

      a0 = ( 1.0D+00 - x )**( c - a - b )
      a = c - a
      b = c - b

    end if

    value = 1.0D+00
    hw = value
    r = 1.0D+00

    do k = 1, 250

      r = r * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
        / ( k * ( c + k - 1.0D+00 ) ) * x

      value = value + r

      if ( abs ( value - hw ) <= abs ( value ) * eps ) then
        exit
      end if

      hw = value

    end do

    value = a0 * value

  end if

  if ( x1 < 0.0D+00 ) then
    x = x1
    c0 = 1.0D+00 / ( 1.0D+00 - x ) ** aa
    value = c0 * value
  end if

  a = aa
  b = bb

  if ( 120 < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_HYPER_2F1 - Warning!'
    write ( *, '(a)' ) '  A large number of iterations were needed.'
    write ( *, '(a)' ) '  The accuracy of the results should be checked.'
  end if

  r8_hyper_2f1 = value

  return
end
function r8_mop ( i )

!*****************************************************************************80
!
!! R8_MOP returns the I-th power of -1 as an R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = rk ) value.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, the power of -1.
!
!    Output, real ( kind = rk ) R8_MOP, the I-th power of -1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  real ( kind = rk ) r8_mop
  real ( kind = rk ) value

  if ( mod ( i, 2 ) == 0 ) then
    value = + 1.0D+00
  else
    value = - 1.0D+00
  end if

  r8_mop = value

  return
end
function r8_psi ( xx )

!*****************************************************************************80
!
!! R8_PSI evaluates the function Psi(X).
!
!  Discussion:
!
!    This routine evaluates the logarithmic derivative of the
!    Gamma function,
!
!      PSI(X) = d/dX ( GAMMA(X) ) / GAMMA(X)
!             = d/dX LN ( GAMMA(X) )
!
!    for real X, where either
!
!      - XMAX1 < X < - XMIN, and X is not a negative integer,
!
!    or
!
!      XMIN < X.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody, Anthony Strecok, Henry Thacher,
!    Chebyshev Approximations for the Psi Function,
!    Mathematics of Computation,
!    Volume 27, Number 121, January 1973, pages 123-127.
!
!  Parameters:
!
!    Input, real ( kind = rk ) XX, the argument of the function.
!
!    Output, real ( kind = rk ) R8_PSI, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) aug
  real ( kind = rk ) den
  integer i
  integer n
  integer nq
  real ( kind = rk ), dimension ( 9 ) :: p1 = (/ &
   4.5104681245762934160D-03, &
   5.4932855833000385356D+00, &
   3.7646693175929276856D+02, &
   7.9525490849151998065D+03, &
   7.1451595818951933210D+04, &
   3.0655976301987365674D+05, &
   6.3606997788964458797D+05, &
   5.8041312783537569993D+05, &
   1.6585695029761022321D+05 /)
  real ( kind = rk ), dimension ( 7 ) :: p2 = (/ &
  -2.7103228277757834192D+00, &
  -1.5166271776896121383D+01, &
  -1.9784554148719218667D+01, &
  -8.8100958828312219821D+00, &
  -1.4479614616899842986D+00, &
  -7.3689600332394549911D-02, &
  -6.5135387732718171306D-21 /)
  real ( kind = rk ), parameter :: piov4 = 0.78539816339744830962D+00
  real ( kind = rk ), dimension ( 8 ) :: q1 = (/ &
   9.6141654774222358525D+01, &
   2.6287715790581193330D+03, &
   2.9862497022250277920D+04, &
   1.6206566091533671639D+05, &
   4.3487880712768329037D+05, &
   5.4256384537269993733D+05, &
   2.4242185002017985252D+05, &
   6.4155223783576225996D-08 /)
  real ( kind = rk ), dimension ( 6 ) :: q2 = (/ &
   4.4992760373789365846D+01, &
   2.0240955312679931159D+02, &
   2.4736979003315290057D+02, &
   1.0742543875702278326D+02, &
   1.7463965060678569906D+01, &
   8.8427520398873480342D-01 /)
  real ( kind = rk ) r8_psi
  real ( kind = rk ) sgn
  real ( kind = rk ) upper
  real ( kind = rk ) w
  real ( kind = rk ) x
  real ( kind = rk ), parameter :: x01 = 187.0D+00
  real ( kind = rk ), parameter :: x01d = 128.0D+00
  real ( kind = rk ), parameter :: x02 = 6.9464496836234126266D-04
  real ( kind = rk ), parameter :: xinf = 1.70D+38
  real ( kind = rk ), parameter :: xlarge = 2.04D+15
  real ( kind = rk ), parameter :: xmax1 = 3.60D+16
  real ( kind = rk ), parameter :: xmin1 = 5.89D-39
  real ( kind = rk ), parameter :: xsmall = 2.05D-09
  real ( kind = rk ) xx
  real ( kind = rk ) z

  x = xx
  w = abs ( x )
  aug = 0.0D+00
!
!  Check for valid arguments, then branch to appropriate algorithm.
!
  if ( xmax1 <= - x .or. w < xmin1 ) then

    if ( 0.0D+00 < x ) then
      r8_psi = - xinf
    else
      r8_psi = xinf
    end if

    return
  end if

  if ( x < 0.5D+00 ) then
!
!  X < 0.5, use reflection formula: psi(1-x) = psi(x) + pi * cot(pi*x)
!  Use 1/X for PI*COTAN(PI*X)  when  XMIN1 < |X| <= XSMALL.
!
    if ( w <= xsmall ) then

      aug = - 1.0D+00 / x
!
!  Argument reduction for cotangent.
!
    else

      if ( x < 0.0D+00 ) then
        sgn = piov4
      else
        sgn = - piov4
      end if

      w = w - real ( int ( w ), kind = rk )
      nq = int ( w * 4.0D+00 )
      w = 4.0D+00 * ( w - real ( nq, kind = rk ) * 0.25D+00 )
!
!  W is now related to the fractional part of 4.0 * X.
!  Adjust argument to correspond to values in the first
!  quadrant and determine the sign.
!
      n = nq / 2

      if ( n + n /= nq ) then
        w = 1.0D+00 - w
      end if

      z = piov4 * w

      if ( mod ( n, 2 ) /= 0 ) then
        sgn = - sgn
      end if
!
!  Determine the final value for  -pi * cotan(pi*x).
!
      n = ( nq + 1 ) / 2
      if ( mod ( n, 2 ) == 0 ) then
!
!  Check for singularity.
!
        if ( z == 0.0D+00 ) then

          if ( 0.0D+00 < x ) then
            r8_psi = - xinf
          else
            r8_psi = xinf
          end if

          return
        end if

        aug = sgn * ( 4.0D+00 / tan ( z ) )

      else

        aug = sgn * ( 4.0D+00 * tan ( z ) )

      end if

    end if

    x = 1.0D+00 - x

  end if
!
!  0.5 <= X <= 3.0.
!
  if ( x <= 3.0D+00 ) then

    den = x
    upper = p1(1) * x
    do i = 1, 7
      den = ( den + q1(i) ) * x
      upper = ( upper + p1(i+1) ) * x
    end do
    den = ( upper + p1(9) ) / ( den + q1(8) )
    x = ( x - x01 / x01d ) - x02
    r8_psi = den * x + aug
    return

  end if
!
!  3.0 < X.
!
  if ( x < xlarge ) then
    w = 1.0D+00 / ( x * x )
    den = w
    upper = p2(1) * w
    do i = 1, 5
      den = ( den + q2(i) ) * w
      upper = ( upper + p2(i+1) ) * w
    end do
    aug = ( upper + p2(7) ) / ( den + q2(6) ) - 0.5D+00 / x + aug
  end if

  r8_psi = aug + log ( x )

  return
end
function r8_sign ( x )

!*****************************************************************************80
!
!! R8_SIGN returns the sign of an R8.
!
!  Discussion:
!
!    value = -1 if X < 0;
!    value = +1 if X => 0.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, the number whose sign is desired.
!
!    Output, real ( kind = rk ) R8_SIGN, the sign of X:
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r8_sign
  real ( kind = rk ) x

  if ( x < 0.0D+00 ) then
    r8_sign = -1.0D+00
  else
    r8_sign = +1.0D+00
  end if

  return
end
function r8_uniform_ab ( a, b )

!*****************************************************************************80
!
!! R8_UNIFORM_AB returns a scaled pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = rk ) value.
!
!    The pseudorandom number should be uniformly distributed
!    between A and B.
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
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, the limits of the interval.
!
!    Output, real ( kind = rk ) R8_UNIFORM_AB, a number strictly between A and B.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) r
  real ( kind = rk ) r8_uniform_ab

  call random_number ( harvest = r )
  r8_uniform_ab = a + ( b - a ) * r

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
subroutine r8poly_print ( n, a, title )

!*****************************************************************************80
!
!! R8POLY_PRINT prints out a polynomial.
!
!  Discussion:
!
!    The power sum form is:
!
!      p(x) = a(0) + a(1) * x + ... + a(n-1) * x^(n-1) + a(n) * x^(n)
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of A.
!
!    Input, real ( kind = rk ) A(0:N), the polynomial coefficients.
!    A(0) is the constant term and
!    A(N) is the coefficient of X^N.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(0:n)
  integer i
  real ( kind = rk ) mag
  character plus_minus
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  if ( n < 0 ) then
    write ( *, '( ''  p(x) = 0'' )' )
    return
  end if

  if ( a(n) < 0.0D+00 ) then
    plus_minus = '-'
  else
    plus_minus = ' '
  end if

  mag = abs ( a(n) )

  if ( 2 <= n ) then
    write ( *, '( ''  p(x) = '', a1, g14.6, '' * x ^ '', i3 )' ) &
      plus_minus, mag, n
  else if ( n == 1 ) then
    write ( *, '( ''  p(x) = '', a1, g14.6, '' * x'' )' ) &
      plus_minus, mag
  else if ( n == 0 ) then
    write ( *, '( ''  p(x) = '', a1, g14.6 )' ) plus_minus, mag
  end if

  do i = n - 1, 0, -1

    if ( a(i) < 0.0D+00 ) then
      plus_minus = '-'
    else
      plus_minus = '+'
    end if

    mag = abs ( a(i) )

    if ( mag /= 0.0D+00 ) then

      if ( 2 <= i ) then
        write ( *, ' ( ''         '', a1, g14.6, '' * x ^ '', i3 )' ) &
          plus_minus, mag, i
      else if ( i == 1 ) then
        write ( *, ' ( ''         '', a1, g14.6, '' * x'' )' ) plus_minus, mag
      else if ( i == 0 ) then
        write ( *, ' ( ''         '', a1, g14.6 )' ) plus_minus, mag
      end if
    end if

  end do

  return
end
function r8vec_in_ab ( n, x, a, b )

!*****************************************************************************80
!
!! R8VEC_IN_AB is TRUE if the entries of an R8VEC are in the range [A,B].
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
!    15 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in X.
!
!    Input, real ( kind = rk ) X(N), the vector.
!
!    Input, real A, B, the limits of the range.
!
!    Output, logical R8VEC_IN_AB, is TRUE if every entry of A is
!    between A and B.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a
  real ( kind = rk ) b
  logical r8vec_in_ab
  real ( kind = rk ) x(n)

  if ( any ( x(1:n) < a .or. b < x(1:n) ) ) then
    r8vec_in_ab = .false.
  else
    r8vec_in_ab = .true.
  end if

  return
end
subroutine r8vec_linspace ( n, a_first, a_last, a )

!*****************************************************************************80
!
!! R8VEC_LINSPACE creates a vector of linearly spaced values.
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
!    14 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real ( kind = rk ) A_FIRST, A_LAST, the first and last entries.
!
!    Output, real ( kind = rk ) A(N), a vector of linearly spaced data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) a_first
  real ( kind = rk ) a_last
  integer i

  if ( n == 1 ) then

    a(1) = ( a_first + a_last ) / 2.0D+00

  else

    do i = 1, n
      a(i) = ( real ( n - i,     kind = rk ) * a_first &
             + real (     i - 1, kind = rk ) * a_last ) &
             / real ( n     - 1, kind = rk )
    end do

  end if

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
subroutine r8vec_uniform_ab ( n, a, b, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_AB returns a scaled pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    Each dimension ranges from A to B.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    31 May 2007
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
!    Input, integer N, the number of entries in the vector.
!
!    Input, real ( kind = rk ) A, B, the lower and upper limits.
!
!    Output, real ( kind = rk ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) r(n)

  call random_number ( harvest = r(1:n) )
  r(1:n) = a + ( b - a ) * r(1:n)

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
subroutine svd_solve ( m, n, a, b, x )

!*****************************************************************************80
!
!! SVD_SOLVE solves a linear system in the least squares sense.
!
!  Discussion:
!
!    The vector X returned by this routine should always minimize the 
!    Euclidean norm of the residual ||A*x-b||.
!
!    If the matrix A does not have full column rank, then there are multiple
!    vectors that attain the minimum residual.  In that case, the vector
!    X returned by this routine is the unique such minimizer that has the 
!    the minimum possible Euclidean norm, that is, ||A*x-b|| and ||x||
!    are both minimized.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 April 2012
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer M, the number of rows of A.
!
!    Input, integer N, the number of columns of A.
!
!    Input, real ( kind = rk ) A(M,N), the matrix.
!
!    Input, real ( kind = rk ) B(M), the right hand side.
!
!    Output, real ( kind = rk ) X(N), the least squares solution.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) a_copy(m,n)
  real ( kind = rk ) b(m)
  real ( kind = rk ) e(max(m+1,n))
  integer i
  integer info
  integer lda
  integer ldu
  integer ldv
  integer job
  real ( kind = rk ) sdiag(max(m+1,n))
  real ( kind = rk ) smax
  real ( kind = rk ) stol
  real ( kind = rk ) sub(n)
  real ( kind = rk ) u(m,m)
  real ( kind = rk ) ub(m)
  real ( kind = rk ) v(n,n)
  real ( kind = rk ), allocatable, dimension ( : ) :: work
  real ( kind = rk ) x(n)
!
!  Get the SVD.
!
  a_copy(1:m,1:n) = a(1:m,1:n)
  lda = m
  ldu = m
  ldv = n
  allocate ( work(1:m) )
  job = 11

  call dsvdc ( a_copy, lda, m, n, sdiag, e, u, ldu, v, ldv, work, job, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SVD_SOLVE - Failure!'
    write ( *, '(a)' ) '  The SVD could not be calculated.'
    write ( *, '(a)' ) '  LINPACK routine DSVDC returned a nonzero'
    write ( *, '(a,i8)' ) '  value of the error flag, INFO = ', info
    stop 1
  end if

  ub(1:m) = matmul ( transpose ( u(1:m,1:m) ), b(1:m) )

  sub(1:n) = 0.0D+00
!
!  For singular problems, there may be tiny but nonzero singular values
!  that should be ignored.  This is a reasonable attempt to avoid such 
!  problems, although in general, the user might wish to control the tolerance.
!
  smax = maxval ( sdiag(1:n) )
  if ( smax <= epsilon ( smax ) ) then
    smax = 1.0D+00
  end if

  stol = epsilon ( smax ) * smax

  do i = 1, n
    if ( i <= m ) then
      if ( stol <= sdiag(i) ) then
        sub(i) = ub(i) / sdiag(i)
      end if
    end if
  end do

  x(1:n) = matmul ( v(1:n,1:n), sub(1:n) )

  return
end
subroutine t_mass_matrix ( n, a )

!*****************************************************************************80
!
!! T_MASS_MATRIX computes the mass matrix for the Chebyshev T polynomial.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(0:n,0:n)
  integer i
  real ( kind = rk ), allocatable :: phi(:,:)
  real ( kind = rk ), allocatable :: phiw(:,:)
  real ( kind = rk ), allocatable :: w(:)
  real ( kind = rk ), allocatable :: x(:)

  allocate ( x(0:n) )
  allocate ( w(0:n) )

  call t_quadrature_rule ( n + 1, x, w )

  allocate ( phi(0:n,0:n) )

  call t_polynomial ( n + 1, n, x, phi )

  allocate ( phiw(0:n,0:n) )

  do i = 0, n
    phiw(0:n,i) = w(i) * phi(i,0:n)
  end do

  a(0:n,0:n) = matmul ( phiw(0:n,0:n), phi(0:n,0:n) )

  deallocate ( phi )
  deallocate ( phiw )
  deallocate ( w )
  deallocate ( x )

  return
end
function t_moment ( e )

!*****************************************************************************80
!
!! T_MOMENT: integral ( -1 <= x <= +1 ) x^e / sqrt ( 1 - x^2 ) dx.
!
!  Discussion:
!
!    Set 
!      x = cos ( theta ), 
!      dx = - sin ( theta ) d theta = - sqrt ( 1 - x^2 ) d theta
!    to transform the integral to
!      integral ( 0 <= theta <= pi ) - ( cos ( theta ) )^e d theta
!    which becomes
!      0 if E is odd,
!      (1/2^e) * choose ( e, e/2 ) * pi if E is even.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer E, the exponent of X.
!    0 <= E.
!
!    Output, real ( kind = rk ) T_MOMENT, the value of the integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer e
  real ( kind = rk ) r8_choose
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) t_moment
  real ( kind = rk ) value

  if ( mod ( e, 2 ) == 1 ) then

    value = 0.0D+00

  else

    value = r8_choose ( e, e / 2 ) * r8_pi / 2.0D+00 ** e

  end if

  t_moment = value

  return
end
subroutine t_polynomial ( m, n, x, v )

!*****************************************************************************80
!
!! T_POLYNOMIAL evaluates Chebyshev polynomials T(n,x).
!
!  Discussion:
!
!    Chebyshev polynomials are useful as a basis for representing the
!    approximation of functions since they are well conditioned, in the sense
!    that in the interval [-1,1] they each have maximum absolute value 1.
!    Hence an error in the value of a coefficient of the approximation, of
!    size epsilon, is exactly reflected in an error of size epsilon between
!    the computed approximation and the theoretical approximation.
!
!    Typical usage is as follows, where we assume for the moment
!    that the interval of approximation is [-1,1].  The value
!    of N is chosen, the highest polynomial to be used in the
!    approximation.  Then the function to be approximated is
!    evaluated at the N+1 points XJ which are the zeroes of the N+1-th
!    Chebyshev polynomial.  Let these values be denoted by F(XJ).
!
!    The coefficients of the approximation are now defined by
!
!      C(I) = 2/(N+1) * sum ( 1 <= J <= N+1 ) F(XJ) T(I,XJ)
!
!    except that C(0) is given a value which is half that assigned
!    to it by the above formula,
!
!    and the representation is
!
!    F(X) approximated by sum ( 0 <= J <= N ) C(J) T(J,X)
!
!    Now note that, again because of the fact that the Chebyshev polynomials
!    have maximum absolute value 1, if the higher order terms of the
!    coefficients C are small, then we have the option of truncating
!    the approximation by dropping these terms, and we will have an
!    exact value for maximum perturbation to the approximation that
!    this will cause.
!
!    It should be noted that typically the error in approximation
!    is dominated by the first neglected basis function (some multiple of
!    T(N+1,X) in the example above).  If this term were the exact error,
!    then we would have found the minimax polynomial, the approximating
!    polynomial of smallest maximum deviation from the original function.
!    The minimax polynomial is hard to compute, and another important
!    feature of the Chebyshev approximation is that it tends to behave
!    like the minimax polynomial while being easy to compute.
!
!    To evaluate a sum like 
!
!      sum ( 0 <= J <= N ) C(J) T(J,X), 
!
!    Clenshaw's recurrence formula is recommended instead of computing the
!    polynomial values, forming the products and summing.
!
!    Assuming that the coefficients C(J) have been computed
!    for J = 0 to N, then the coefficients of the representation of the
!    indefinite integral of the function may be computed by
!
!      B(I) = ( C(I-1) - C(I+1))/2*(I-1) for I=1 to N+1, 
!
!    with
! 
!      C(N+1)=0
!      B(0) arbitrary.  
!
!    Also, the coefficients of the representation of the derivative of the 
!    function may be computed by:
!
!      D(I) = D(I+2)+2*I*C(I) for I=N-1, N-2, ..., 0, 
!
!    with
!
!      D(N+1) = D(N)=0.
!
!    Some of the above may have to adjusted because of the irregularity of C(0).
!
!    The formula is:
!
!      T(N,X) = COS(N*ARCCOS(X))
!
!  Differential equation:
!
!    (1-X*X) Y'' - X Y' + N N Y = 0
!
!  First terms:
!
!    T(0,X) =  1
!    T(1,X) =  1 X
!    T(2,X) =  2 X^2 -   1
!    T(3,X) =  4 X^3 -   3 X
!    T(4,X) =  8 X^4 -   8 X^2 +  1
!    T(5,X) = 16 X^5 -  20 X^3 +  5 X
!    T(6,X) = 32 X^6 -  48 X^4 + 18 X^2 - 1
!    T(7,X) = 64 X^7 - 112 X^5 + 56 X^3 - 7 X
!
!  Inequality:
!
!    abs ( T(N,X) ) <= 1 for -1 <= X <= 1
!
!  Orthogonality:
!
!    For integration over [-1,1] with weight
!
!      W(X) = 1 / sqrt(1-X*X), 
!
!    if we write the inner product of T(I,X) and T(J,X) as
!
!      < T(I,X), T(J,X) > = integral ( -1 <= X <= 1 ) W(X) T(I,X) T(J,X) dX
!
!    then the result is:
!
!      < T(I,X), T(J,X) > = 0    if I /= J
!      < T(I,X), T(J,X) > = PI/2 if I == J /= 0
!      < T(I,X), T(J,X) > = PI   if I == J == 0
!
!    A discrete orthogonality relation is also satisfied at each of
!    the N zeroes of T(N,X):  sum ( 1 <= K <= N ) T(I,X) * T(J,X)
!                              = 0 if I /= J
!                              = N/2 if I == J /= 0
!                              = N if I == J == 0
!
!  Recursion:
!
!    T(0,X) = 1,
!    T(1,X) = X,
!    T(N,X) = 2 * X * T(N-1,X) - T(N-2,X)
!
!    T'(N,X) = N * ( -X * T(N,X) + T(N-1,X) ) / ( 1 - X^2 )
!
!  Special values:
!
!    T(N,1) = 1
!    T(N,-1) = (-1)^N
!    T(2N,0) = (-1)^N
!    T(2N+1,0) = 0
!    T(N,X) = (-1)^N * T(N,-X)
!
!  Zeroes:
!
!    M-th zero of T(N,X) is X = cos((2*M-1)*PI/(2*N)), M = 1 to N.
!
!  Extrema:
!
!    M-th extremum of T(N,X) is X = cos(PI*M/N), M = 0 to N.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer M, the number of evaluation points.
!
!    Input, integer N, the highest polynomial to compute.
!
!    Input, real ( kind = rk ) X(1:M), the evaluation points.
!
!    Output, real ( kind = rk ) V(1:M,0:N), the values of the polynomials.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer j
  real ( kind = rk ) v(1:m,0:n)
  real ( kind = rk ) x(1:m)

  if ( n < 0 ) then
    return
  end if

  v(1:m,0) = 1.0D+00

  if ( n < 1 ) then
    return
  end if

  v(1:m,1) = x(1:m)
 
  do j = 2, n
    v(1:m,j) = 2.0D+00 * x(1:m) * v(1:m,j-1) - v(1:m,j-2)
  end do
 
  return
end
subroutine t_polynomial_01_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! T_POLYNOMIAL_01_VALUES: values of shifted Chebyshev polynomials T(n,x).
!
!  Discussion:
!
!    The shifted Chebyshev polynomial T01(n,x) = T(n,2*x-1).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer N, the order of the function.
!
!    Output, real ( kind = rk ) X, the point where the function is evaluated.
!
!    Output, real ( kind = rk ) FX, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 25

  real ( kind = rk ) fx
  real ( kind = rk ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.0000000000000000D+00, &
     1.0000000000000000D+00, &
     0.7000000000000000D+00, &
    -0.0200000000000000D+00, &
    -0.7280000000000000D+00, &
    -0.9992000000000000D+00, &
    -0.6708800000000000D+00, &
     0.0599680000000000D+00, &
     0.7548352000000000D+00, &
     0.9968012800000000D+00, &
     0.6406865920000000D+00, &
    -0.0998400512000000D+00, &
    -0.7804626636800000D+00, &
    -0.9928076779520000D+00, &
    -1.0000000000000000D+00, &
     0.2063872000000000D+00, &
    -0.9784704000000000D+00, &
     0.2580224000000000D+00, &
     0.9870208000000000D+00, &
     0.0000000000000000D+00, &
    -0.9870208000000000D+00, &
    -0.2580224000000000D+00, &
     0.9784704000000000D+00, &
    -0.2063872000000000D+00, &
     1.0000000000000000D+00 /)
  integer n
  integer n_data
  integer, save, dimension ( n_max ) :: n_vec = (/ &
    -1, &
     0,  1,  2, &
     3,  4,  5, &
     6,  7,  8, &
     9, 10, 11, &
    12,  7,  7, &
     7,  7,  7, &
     7,  7,  7, &
     7,  7,  7 /)
  real ( kind = rk ) x
  real ( kind = rk ), save, dimension ( n_max ) :: x_vec = (/ &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.00D+00, &
    0.10D+00, &
    0.20D+00, &
    0.30D+00, &
    0.40D+00, &
    0.50D+00, &
    0.60D+00, &
    0.70D+00, &
    0.80D+00, &
    0.90D+00, &
    1.00D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine t_polynomial_ab ( a, b, m, n, xab, v )

!*****************************************************************************80
!
!! T_POLYNOMIAL_AB: evaluates Chebyshev polynomials TAB(n,x) in [A,B].
!
!  Discussion:
!
!    TAB(n,x) = T(n,(2*x-a-b)/(b-a))
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, the domain of definition.
!
!    Input, integer M, the number of evaluation points.
!
!    Input, integer N, the highest polynomial to compute.
!
!    Input, real ( kind = rk ) XAB(M), the evaluation points.
!    It must be the case that A <= XAB(*) <= B.
!
!    Output, real ( kind = rk ) V(M,N+1), the values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) v(1:m,0:n)
  real ( kind = rk ) x(1:m)
  real ( kind = rk ) xab(1:m)

  x(1:m) = ( 2.0D+00 * xab(1:m) - a - b ) / ( b - a )

  call t_polynomial ( m, n, x, v )
 
  return
end
function t_polynomial_ab_value ( a, b, n, xab )

!*****************************************************************************80
!
!! T_POLYNOMIAL_AB_VALUE: evaluates Chebyshev polynomials TAB(n,x) in [A,B].
!
!  Discussion:
!
!    TAB(n,x) = T(n,(2*x-a-b)/(b-a))
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, the domain of definition.
!
!    Input, integer N, the highest polynomial to compute.
!
!    Input, real ( kind = rk ) XAB, the evaluation point.
!    It must be the case that A <= XAB <= B.
!
!    Output, real ( kind = rk ) T_POLYNOMIAL_AB_VALUE, the value.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) t_polynomial_ab_value
  real ( kind = rk ) t_polynomial_value
  real ( kind = rk ) value
  real ( kind = rk ) x
  real ( kind = rk ) xab

  x = ( 2.0D+00 * xab - a - b ) / ( b - a )

  value = t_polynomial_value ( n, x )
 
  t_polynomial_ab_value = value

  return
end
subroutine t_polynomial_coefficients ( n, c )

!*****************************************************************************80
!
!! T_POLYNOMIAL_COEFFICIENTS: coefficients of the Chebyshev polynomial T(n,x).
!
!  First terms:
!
!    N/K     0     1      2      3       4     5      6    7      8    9   10
!
!     0      1
!     1      0     1
!     2     -1     0      2
!     3      0    -3      0      4
!     4      1     0     -8      0       8
!     5      0     5      0    -20       0    16
!     6     -1     0     18      0     -48     0     32
!     7      0    -7      0     56       0  -112      0    64
!
!  Recursion:
!
!    T(0,X) = 1,
!    T(1,X) = X,
!    T(N,X) = 2 * X * T(N-1,X) - T(N-2,X)
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Output, real ( kind = rk ) C(0:N,0:N), the coefficients of the Chebyshev T
!    polynomials.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) c(0:n,0:n)
  integer i

  if ( n < 0 ) then
    return
  end if

  c(0:n,0:n) = 0.0D+00

  c(0,0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  c(1,1) = 1.0D+00
 
  do i = 2, n
    c(i,0)     =                        - c(i-2,0)
    c(i,1:i-2) = 2.0D+00 * c(i-1,0:i-3) - c(i-2,1:i-2)
    c(i,  i-1) = 2.0D+00 * c(i-1,  i-2)
    c(i,  i  ) = 2.0D+00 * c(i-1,  i-1)
  end do
 
  return
end
subroutine t_polynomial_plot ( n_num, n_val, output_filename )

!*****************************************************************************80
!
!! T_POLYNOMIAL_PLOT plots Chebyshev polynomials T(n,x).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N_NUM, the number of polynomials to be plotted.
!
!    Input, integer N_VAL(N_NUM), the degrees of 1 or more 
!    Chebyshev polynomials to be plotted together.
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the name into which the 
!    graphics information is to be stored.  Note that the PNG format will 
!    be used.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 501
  integer n_num

  real ( kind = rk ) a
  real ( kind = rk ) b
  integer column
  character ( len = 255 ) command_filename
  integer command_unit
  character ( len = 255 ) data_filename
  integer data_unit
  integer i
  integer i4vec_max
  integer j
  integer n
  integer n_max
  integer n_val(n_num)
  character ( len = * ) output_filename
  real ( kind = rk ), allocatable :: v(:,:)
  real ( kind = rk ) x(m)

  a = -1.0D+00
  b = +1.0D+00

  call r8vec_linspace ( m, a, b, x )
!
!  Compute all the data.
!
  n_max = i4vec_max ( n_num, n_val )
  allocate ( v(m,0:n_max) )
  call t_polynomial ( m, n_max, x, v )
!
!  Create the data file.
!
  data_filename = 't_polynomial_data.txt'
  call get_unit ( data_unit )
  open ( unit = data_unit, file = data_filename, status = 'replace' )
  do i = 1, m
    write ( data_unit, '(2x,g14.6)', advance = 'no' ) x(i)
    do j = 1, n_num
      n = n_val(j)
      write ( data_unit, '(2x,g14.6)', advance = 'no' ) v(i,n)
    end do
    write ( data_unit, '(a)' ) ''
  end do
  close ( unit = data_unit )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) &
    '  Created graphics data file "' // trim ( data_filename ) // '".'
!
!  Plot the selected data.
!
  command_filename = 't_polynomial_commands.txt'
  call get_unit ( command_unit )
  open ( unit = command_unit, file = command_filename, status = 'replace' )

  write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) 'set nokey'
  write ( command_unit, '(a)' ) &
    'set output "' // trim ( output_filename ) // '"'
  write ( command_unit, '(a)' ) 'set xlabel "<---X--->"'
  write ( command_unit, '(a)' ) 'set ylabel "<---T(n,x)--->"'
  write ( command_unit, '(a)' ) &
    'set title "Chebyshev Polynomials T(n,x)"'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set style data lines'
  do j = 1, n_num
    column = n_val(j) + 1
    if ( j == 1 ) then
      write ( command_unit, '(a)', advance = 'no' ) 'plot '
    else
      write ( command_unit, '(a)', advance = 'no' ) '     '
    end if
    write ( command_unit, '(a,i2,a)', advance = 'no' ) &
      '"' // trim ( data_filename ) // &
      '" using 1:', column, ' lw 3 linecolor rgb "red"'
    if ( j < n_num ) then
      write ( command_unit, '(a)' ) ', \'
    else
      write ( command_unit, '(a)' ) ''
    end if
  end do

  close ( unit = command_unit )
  write ( *, '(a)' ) &
    '  Created graphics command file "' // trim ( command_filename ) // '".'

  deallocate ( v )

  return
end
function t_polynomial_value ( n, x )

!*****************************************************************************80
!
!! T_POLYNOMIAL_VALUE: returns the single value T(n,x).
!
!  Discussion:
!
!    In cases where calling T_POLYNOMIAL is inconvenient, because it returns
!    a vector of values for multiple arguments X, this simpler interface
!    may be appropriate.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the polynomial.
!
!    Input, real ( kind = rk ) X, the argument of the polynomial.
!
!    Output, real ( kind = rk ) T_POLYNOMIAL_VALUE, the value of T(n,x).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n
  real ( kind = rk ) t_polynomial_value
  real ( kind = rk ) value
  real ( kind = rk ), allocatable :: vec(:)
  real ( kind = rk ) x
  real ( kind = rk ) x_vec(1)

  if ( n < 0 ) then

    value = 0.0D+00

  else

    m = 1
    allocate ( vec(0:n) )

    x_vec(1) = x
    call t_polynomial ( m, n, x_vec, vec )

    value = vec(n)
    deallocate ( vec )

  end if

  t_polynomial_value = value

  return
end
subroutine t_polynomial_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! T_POLYNOMIAL_VALUES returns values of Chebyshev polynomials T(n,x).
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      ChebyshevT[n,x]
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer N, the order of the function.
!
!    Output, real ( kind = rk ) X, the point where the function is evaluated.
!
!    Output, real ( kind = rk ) FX, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 14

  real ( kind = rk ) fx
  real ( kind = rk ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.0000000000000000D+00, &
     0.1000000000000000D+01, &
     0.8000000000000000D+00, &
     0.2800000000000000D+00, &
    -0.3520000000000000D+00, &
    -0.8432000000000000D+00, &
    -0.9971200000000000D+00, &
    -0.7521920000000000D+00, &
    -0.2063872000000000D+00, &
     0.4219724800000000D+00, &
     0.8815431680000000D+00, &
     0.9884965888000000D+00, &
     0.7000513740800000D+00, &
     0.1315856097280000D+00 /)
  integer n
  integer n_data
  integer, save, dimension ( n_max ) :: n_vec = (/ &
    -1, &
     0,  1,  2, &
     3,  4,  5, &
     6,  7,  8, &
     9, 10, 11, &
    12 /)
  real ( kind = rk ) x
  real ( kind = rk ), save, dimension ( n_max ) :: x_vec = (/ &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine t_polynomial_zeros ( n, z )

!*****************************************************************************80
!
!! T_POLYNOMIAL_ZEROS returns zeroes of the Chebyshev polynomial T(n,x).
!
!  Discussion:
!
!    The I-th zero of T(N,X) is cos((2*I-1)*PI/(2*N)), I = 1 to N
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 March 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the polynomial.
!
!    Output, real ( kind = rk ) Z(N), the zeroes of T(N,X).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) angle
  integer i
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) z(n)

  do i = 1, n
    angle = real ( 2 * i - 1, kind = rk ) * r8_pi / real ( 2 * n, kind = rk );
    z(i) = cos ( angle );
  end do

  return
end
subroutine t_project_coefficients ( n, f, c )

!*****************************************************************************80
!
!! T_PROJECT_COEFFICIENTS: function projected onto Chebyshev polynomials T(n,x).
!
!  Discussion:
!
!    It is assumed that the interval of definition is -1 <= x <= +1.
!
!    Over this interval, f(x) will be well approximated by
!
!      f(x) approx sum ( 0 <= i <= n ) c(i) * T(i,x)
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the highest order polynomial to compute.
!
!    Input, external real ( kind = rk ) function F ( X ), evaluates the function.
!
!    Output, real ( kind = rk ) C(0:N), the projection coefficients of f(x) onto
!    T(0,x) through T(n,x).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) c(0:n)
  real ( kind = rk ) d(0:n)
  real ( kind = rk ), external :: f
  real ( kind = rk ) fac
  integer j
  integer k
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) total
  real ( kind = rk ) y

  do k = 0, n
    y = cos ( r8_pi * ( real ( k, kind = rk ) + 0.5D+00 ) &
      / real ( n + 1, kind = rk ) )
    d(k) = f ( y )
  end do

  fac = 2.0D+00 / real ( n + 1, kind = rk )

  do j = 0, n
    total = 0.0D+00
    do k = 0, n
      total = total + d(k) * cos ( ( r8_pi * real ( j, kind = rk ) ) &
        * ( ( real ( k, kind = rk ) + 0.5D+00 ) / real ( n + 1, kind = rk ) ) )
    end do
    c(j) = fac * total
  end do

  c(0) = c(0) / 2.0D+00

  return
end
subroutine t_project_coefficients_ab ( n, f, a, b, c )

!*****************************************************************************80
!
!! T_PROJECT_COEFFICIENTS_AB: function projected onto TAB(n,x) over [a,b].
!
!  Discussion:
!
!    TAB(n,x) = T(n,(2*x-a-b)/(b-a))
!
!    It is assumed that the interval of definition is a <= x <= b.
!
!    Over this interval, f(x) will be well approximated by
!
!      f(x) approx sum ( 0 <= i <= n ) c(i) * T(i,(2x-a-b)/(b-a))
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the highest order polynomial to compute.
!
!    Input, external real ( kind = rk ) function F ( X ), evaluates the function.
!
!    Input, real ( kind = rk ) A, B, the interval of definition.
!
!    Output, real ( kind = rk ) C(0:N), the projection coefficients of f(x) onto
!    T(0,x) through T(n,x).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c(0:n)
  real ( kind = rk ) d(0:n)
  real ( kind = rk ), external :: f
  real ( kind = rk ) fac
  integer j
  integer k
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) t
  real ( kind = rk ) total
  real ( kind = rk ) y

  do k = 0, n

    t = cos ( r8_pi * ( real ( k, kind = rk ) - 0.5D+00 ) &
      / real ( n + 1, kind = rk ) )

    y = ( ( 1.0D+00 + t ) * b   &
        + ( 1.0D+00 - t ) * a ) &
        /   2.0D+00

    d(k) = f ( y )

  end do

  fac = 2.0D+00 / real ( n + 1, kind = rk )

  do j = 0, n
    total = 0.0D+00
    do k = 0, n
      total = total + d(k) * cos ( ( r8_pi * real ( j, kind = rk ) ) &
        * ( ( real ( k, kind = rk ) + 0.5D+00 ) / real ( n + 1, kind = rk ) ) )
    end do
    c(j) = fac * total
  end do

  c(0) = c(0) / 2.0D+00

  return
end
subroutine t_project_coefficients_data ( a, b, m, n, x, d, c )

!*****************************************************************************80
!
!! T_PROJECT_COEFFICIENTS_DATA: project data onto Chebyshev polynomials T(n,x).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, the domain of definition.
!
!    Input, integer M, the number of data values.
!
!    Input, integer N, the desired order of the Chebyshev 
!    expansion.
!
!    Input, real ( kind = rk ) X(M), the data abscissas.  These need not 
!    be sorted.  It must be the case that A <= X() <= B.
!
!    Input, real ( kind = rk ) D(M), the data values.
!
!    Output, real ( kind = rk ) C(0:N), the approximate Chebshev coefficients.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c(0:n)
  real ( kind = rk ) d(m)
  logical r8vec_in_ab
  real ( kind = rk ) v(m,0:n)
  real ( kind = rk ) x(m)

  if ( .not. r8vec_in_ab ( m, x, a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) ' T_PROJECT_COEFFICIENTS_DATA- Fatal error!'
    write ( *, '(a)' ) '  Some X not in [A,B].'
    stop 1
  end if
!
!  Compute the M by N+1 Chebyshev Vandermonde matrix V.
!
  call t_polynomial_ab ( a, b, m, n, x, v )
!
!  Compute the least-squares solution C.
!
  call svd_solve ( m, n + 1, v, d, c )

  return
end
subroutine t_project_value ( m, n, x, c, v )

!*****************************************************************************80
!
!! T_PROJECT_VALUE evaluates an expansion in Chebyshev polynomials T(n,x).
!
!  Discussion:
!
!    The projection is assumed to be based on the interval [-1,+1].
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of evaluation points.
!
!    Input, integer N, the highest order polynomial to compute.
!
!    Input, real ( kind = rk ) X(M), the evaluation points.
!
!    Input, real ( kind = rk ) C(0:N), the expansion coefficients.
!
!    Output, real ( kind = rk ) V(M), the value of the Chebyshev function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n 

  real ( kind = rk ) b0(m)
  real ( kind = rk ) b1(m)
  real ( kind = rk ) b2(m)
  real ( kind = rk ) c(0:n)
  integer j
  real ( kind = rk ) v(m)
  real ( kind = rk ) x(m)

  b1(1:m) = 0.0D+00
  b0(1:m) = 0.0D+00

  do j = n, 0, -1
    b2(1:m) = b1(1:m)
    b1(1:m) = b0(1:m)
    b0(1:m) = c(j) + 2.0D+00 * x(1:m) * b1(1:m) - b2(1:m)
  end do

  v(1:m) = 0.5D+00 * ( c(0) + b0(1:m) - b2(1:m) )

  return
end
subroutine t_project_value_ab ( m, n, x, c, a, b, v )

!*****************************************************************************80
!
!! T_PROJECT_VALUE_AB evaluates an expansion in Chebyshev polynomials TAB(n,x).
!
!  Discussion:
!
!    TAB(n,x) = T(n,(2*x-a-b)/(b-a))
!
!    The projection is assumed to be based on the interval [A,B].
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of evaluation points.
!
!    Input, integer N, the highest order polynomial to compute.
!
!    Input, real ( kind = rk ) X(M), the evaluation points.
!
!    Input, real ( kind = rk ) C(0:N), the expansion coefficients.
!
!    Input, real ( kind = rk ) A, B, the interval of definition.
!
!    Output, real ( kind = rk ) V(M), the value of the Chebyshev function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n 

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) b0(m)
  real ( kind = rk ) b1(m)
  real ( kind = rk ) b2(m)
  real ( kind = rk ) c(0:n)
  integer j
  real ( kind = rk ) v(m)
  real ( kind = rk ) x(m)

  b1(1:m) = 0.0D+00
  b0(1:m) = 0.0D+00

  do j = n, 0, -1
    b2(1:m) = b1(1:m)
    b1(1:m) = b0(1:m)
    b0(1:m) = c(j) + 2.0D+00 / ( b - a ) * ( 2.0D+00 * x(1:m) - a - b ) &
      * b1(1:m) - b2(1:m)
  end do

  v(1:m) = 0.5D+00 * ( c(0) + b0(1:m) - b2(1:m) )

  return
end
subroutine t_quadrature_rule ( n, t, w )

!*****************************************************************************80
!
!! T_QUADRATURE_RULE: quadrature rule for T(n,x).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the rule.
!
!    Output, real ( kind = rk ) T(N), W(N), the points and weights of the rule.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) bj(n)
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) t(n)
  real ( kind = rk ) w(n)

  t(1:n) = 0.0D+00

  bj(1) = sqrt ( 0.5D+00 )
  bj(2:n) = 0.5D+00

  w(1) = sqrt ( r8_pi )
  w(2:n) = 0.0D+00

  call imtqlx ( n, t, bj, w )

  w(1:n) = w(1:n) ** 2

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
function tt_product ( i, j, x )

!*****************************************************************************80
!
!! TT_PRODUCT: evaluate T(i,x)*T(j,x)
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, J, the indices.
!
!    Input, real ( kind = rk ) X, the argument.
!
!    Output, real ( kind = rk ) TT_PRODUCT, the value.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer imj
  integer ipj
  integer j
  real ( kind = rk ) t_polynomial_value
  real ( kind = rk ) timj
  real ( kind = rk ) tipj
  real ( kind = rk ) tt_product
  real ( kind = rk ) value
  real ( kind = rk ) x

  if ( i < 0 .or. j < 0 ) then
    value = 0.0D+00
  else
    ipj = i + j
    tipj = t_polynomial_value ( ipj, x )
    imj = abs ( i - j )
    timj = t_polynomial_value ( imj, x )
    value = 0.5D+00 * ( tipj + timj )
  end if

  tt_product = value

  return
end
function tt_product_integral ( i, j )

!*****************************************************************************80
!
!! TT_PRODUCT_INTEGRAL: integral (-1<=x<=1) T(i,x)*T(j,x)/sqrt(1-x^2) dx
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, J, the polynomial indices.
!    0 <= I, J.
!
!    Output, real ( kind = rk ) TT_PRODUCT_INTEGRAL, the integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer j
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) tt_product_integral
  real ( kind = rk ) value

  if ( i < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TT_PRODUCT_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  0 <= I is required.'
    stop 1
  end if

  if ( j < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TT_PRODUCT_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  0 <= J is required.'
    stop 1
  end if

  if ( i /= j ) then
    value = 0.0D+00
  elseif ( i == 0 ) then
    value = r8_pi
  elseif ( 0 < i ) then
    value = r8_pi / 2.0D+00
  end if

  tt_product_integral = value

  return
end
function ttt_product_integral ( i, j, k )

!*****************************************************************************80
!
!! TTT_PRODUCT_INTEGRAL: int (-1<=x<=1) T(i,x)*T(j,x)*T(k,x)/sqrt(1-x^2) dx
!  
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    John Mason, David Handscomb,
!    Chebyshev Polynomials,
!    CRC Press, 2002,
!    ISBN: 0-8493-035509,
!    LC: QA404.5.M37.
!
!  Parameters:
!
!    Input, integer I, J, K, the polynomial indices.
!    0 <= I, J, K.
!
!    Output, real ( kind = rk ) TTT_PRODUCT_INTEGRAL, the integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer j
  integer k
  real ( kind = rk ) tt_product_integral
  real ( kind = rk ) ttt_product_integral
  real ( kind = rk ) value

  if ( i < 0 ) then
    value = 0.0D+00
  else if ( j < 0 ) then
    value = 0.0D+00
  else if ( k < 0 ) then
    value = 0.0D+00
  else
    value = 0.5D+00 * ( &
        tt_product_integral (       i + j,   k ) &
      + tt_product_integral ( abs ( i - j ), k ) )
  end if

  ttt_product_integral = value

  return
end
function tu_product ( i, j, x )

!*****************************************************************************80
!
!! TU_PRODUCT: evaluate T(i,x)*U(j,x)
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, J, the indices.
!
!    Input, real ( kind = rk ) X, the argument.
!
!    Output, real ( kind = rk ) TU_PRODUCT, the value.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer j
  real ( kind = rk ) tu_product
  real ( kind = rk ) u_polynomial_value
  real ( kind = rk ) uu_product
  real ( kind = rk ) value
  real ( kind = rk ) x

  if ( i < 0 ) then
    value = 0.0D+00
  else if ( j < 0 ) then
    value = 0.0D+00
  else if ( i == 0 ) then
    value = u_polynomial_value ( j, x )
  else
    value = 0.5D+00 * ( uu_product ( i, j, x ) - uu_product ( i - 2, j, x ) )
  end if

  tu_product = value

  return
end
subroutine u_mass_matrix ( n, a )

!*****************************************************************************80
!
!! U_MASS_MATRIX computes the mass matrix for the Chebyshev U polynomial.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(0:n,0:n)
  integer i
  real ( kind = rk ), allocatable :: phi(:,:)
  real ( kind = rk ), allocatable :: phiw(:,:)
  real ( kind = rk ), allocatable :: w(:)
  real ( kind = rk ), allocatable :: x(:)

  allocate ( x(0:n) )
  allocate ( w(0:n) )

  call u_quadrature_rule ( n + 1, x, w )

  allocate ( phi(0:n,0:n) )

  call u_polynomial ( n + 1, n, x, phi )

  allocate ( phiw(0:n,0:n) )

  do i = 0, n
    phiw(0:n,i) = w(i) * phi(i,0:n)
  end do

  a(0:n,0:n) = matmul ( phiw(0:n,0:n), phi(0:n,0:n) )

  deallocate ( phi )
  deallocate ( phiw )
  deallocate ( w )
  deallocate ( x )

  return
end
function u_moment ( e )

!*****************************************************************************80
!
!! U_MOMENT: integral ( -1 <= x <= +1 ) x^e sqrt ( 1 - x^2 ) dx.
!
!  Discussion:
!
!     E    U_INTEGRAL
!    --    --------------
!     0         pi /    2   
!     2         pi /    8
!     4         pi /   16
!     6     5 * pi /  128
!     8     7 * pi /  256
!    10    21 * pi / 1024
!    12    33 * pi / 2048
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer E, the exponent of X.
!    0 <= E.
!
!    Output, real ( kind = rk ) U_MOMENT, the value of the integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) arg1
  real ( kind = rk ) arg2
  integer e
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) u_moment
  real ( kind = rk ) value

  if ( mod ( e, 2 ) == 1 ) then

    value = 0.0D+00

  else

    arg1 = 0.5D+00 * real ( 1 + e, kind = rk )
    arg2 = 2.0D+00 + 0.5D+00 * real ( e, kind = rk )
    value = 0.5D+00 * sqrt ( r8_pi ) * gamma ( arg1 ) / gamma ( arg2 )

  end if

  u_moment = value

  return
end
subroutine u_polynomial ( m, n, x, v )

!*****************************************************************************80
!
!! U_POLYNOMIAL evaluates Chebyshev polynomials U(n,x).
!
!  Discussion:
!
!    The formula is:
!
!      If |X| <= 1, then
!
!        U(N,X) = sin ( (N+1) * arccos(X) ) / sqrt ( 1 - X^2 )
!               = sin ( (N+1) * arccos(X) ) / sin ( arccos(X) )
!
!      else
!
!        U(N,X) = sinh ( (N+1) * arccosh(X) ) / sinh ( arccosh(X) )
!
!  Differential equation:
!
!    (1-X*X) Y'' - 3 X Y' + N (N+2) Y = 0
!
!  First terms:
!
!    U(0,X) =   1
!    U(1,X) =   2 X
!    U(2,X) =   4 X^2 -   1
!    U(3,X) =   8 X^3 -   4 X
!    U(4,X) =  16 X^4 -  12 X^2 +  1
!    U(5,X) =  32 X^5 -  32 X^3 +  6 X
!    U(6,X) =  64 X^6 -  80 X^4 + 24 X^2 - 1
!    U(7,X) = 128 X^7 - 192 X^5 + 80 X^3 - 8X
!
!  Orthogonality:
!
!    For integration over [-1,1] with weight
!
!      W(X) = sqrt(1-X*X), 
!
!    we have
!
!      < U(I,X), U(J,X) > = integral ( -1 <= X <= 1 ) W(X) U(I,X) U(J,X) dX 
!
!    then the result is:
!
!      < U(I,X), U(J,X) >  =  0    if I /= J
!      < U(I,X), U(J,X) >  =  PI/2 if I == J
!
!  Recursion:
!
!    U(0,X) = 1,
!    U(1,X) = 2 * X,
!    U(N,X) = 2 * X * U(N-1,X) - U(N-2,X)
!
!  Special values:
!
!    U(N,1) = N + 1
!    U(2N,0) = (-1)^N
!    U(2N+1,0) = 0
!    U(N,X) = (-1)^N * U(N,-X)
!
!  Zeroes:
!
!    M-th zero of U(N,X) is X = cos( M*PI/(N+1)), M = 1 to N
!
!  Extrema:
!
!    M-th extremum of U(N,X) is X = cos( M*PI/N), M = 0 to N
!
!  Norm:
!
!    Integral ( -1 <= X <= 1 ) ( 1 - X^2 ) * U(N,X)^2 dX = PI/2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer M, the number of evaluation points.
!
!    Input, integer N, the highest polynomial to compute.
!
!    Input, real ( kind = rk ) X(M), the evaluation points.
!
!    Output, real ( kind = rk ) V(M,0:N), the values of the N+1 Chebyshev
!    polynomials.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer i
  real ( kind = rk ) v(m,0:n)
  real ( kind = rk ) x(m)

  if ( n < 0 ) then
    return
  end if

  v(1:m,0) = 1.0D+00

  if ( n < 1 ) then
    return
  end if

  v(1:m,1) = 2.0D+00 * x(1:m)

  do i = 2, n
    v(1:m,i) = 2.0D+00 * x(1:m) * v(1:m,i-1) - v(1:m,i-2)
  end do
 
  return
end
subroutine u_polynomial_01_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! U_POLYNOMIAL_01_VALUES: values of shifted Chebyshev polynomials U01(n,x).
!
!  Discussion:
!
!    The shifted Chebyshev polynomial U01(n,x) = U(n,2*x-1).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer N, the order of the function.
!
!    Output, real ( kind = rk ) X, the point where the function is evaluated.
!
!    Output, real ( kind = rk ) FX, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 25

  real ( kind = rk ) fx
  real ( kind = rk ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.000000000000000D+00, &
     1.000000000000000D+00, &
     1.400000000000000D+00, &
     0.9600000000000000D+00, &
    -0.05600000000000000D+00, &
    -1.038400000000000D+00, &
    -1.397760000000000D+00, &
    -0.9184640000000000D+00, &
     0.1119104000000000D+00, &
     1.075138560000000D+00, &
     1.393283584000000D+00, &
     0.8754584576000000D+00, &
    -0.1676417433600000D+00, &
    -1.110156898304000D+00, &
    -8.000000000000000D+00, &
     1.511014400000000D+00, &
    -1.133260800000000D+00, &
    -0.1636352000000000D+00, &
     1.019801600000000D+00, &
     0.000000000000000D+00, &
    -1.019801600000000D+00, &
     0.1636352000000000D+00, &
     1.133260800000000D+00, &
    -1.511014400000000D+00, &
     8.000000000000000D+00 /)
  integer n
  integer n_data
  integer, save, dimension ( n_max ) :: n_vec = (/ &
    -1, &
     0,  1,  2, &
     3,  4,  5, &
     6,  7,  8, &
     9, 10, 11, &
    12,  7,  7, &
     7,  7,  7, &
     7,  7,  7, &
     7,  7,  7 /)
  real ( kind = rk ) x
  real ( kind = rk ), save, dimension ( n_max ) :: x_vec = (/ &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.00D+00, &
    0.10D+00, &
    0.20D+00, &
    0.30D+00, &
    0.40D+00, &
    0.50D+00, &
    0.60D+00, &
    0.70D+00, &
    0.80D+00, &
    0.90D+00, &
    1.00D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine u_polynomial_ab ( a, b, m, n, xab, v )

!*****************************************************************************80
!
!! U_POLYNOMIAL_AB: evaluates Chebyshev polynomials UAB(n,x) in [A,B].
!
!  Discussion:
!
!    UAB(n,x) = U(n,(2*x-a-b)/(b-a))
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, the domain of definition.
!
!    Input, integer M, the number of evaluation points.
!
!    Input, integer N, the highest polynomial to compute.
!
!    Input, real ( kind = rk ) XAB(M), the evaluation points.
!    It must be the case that A <= XAB(*) <= B.
!
!    Output, real ( kind = rk ) V(M,N+1), the values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) v(1:m,0:n)
  real ( kind = rk ) x(1:m)
  real ( kind = rk ) xab(1:m)

  x(1:m) = ( 2.0D+00 * xab(1:m) - a - b ) / ( b - a )

  call u_polynomial ( m, n, x, v )
 
  return
end
function u_polynomial_ab_value ( a, b, n, xab )

!*****************************************************************************80
!
!! U_POLYNOMIAL_AB_VALUE: evaluates Chebyshev polynomials UAB(n,x) in [A,B].
!
!  Discussion:
!
!    UAB(n,x) = U(n,(2*x-a-b)/(b-a))
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, the domain of definition.
!
!    Input, integer N, the highest polynomial to compute.
!
!    Input, real ( kind = rk ) XAB, the evaluation point.
!    It must be the case that A <= XAB <= B.
!
!    Output, real ( kind = rk ) U_POLYNOMIAL_AB_VALUE, the value.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) u_polynomial_ab_value
  real ( kind = rk ) u_polynomial_value
  real ( kind = rk ) value
  real ( kind = rk ) x
  real ( kind = rk ) xab

  x = ( 2.0D+00 * xab - a - b ) / ( b - a )

  value = u_polynomial_value ( n, x )
 
  u_polynomial_ab_value = value

  return
end
subroutine u_polynomial_coefficients ( n, c )

!*****************************************************************************80
!
!! U_POLYNOMIAL_COEFFICIENTS: coefficients of Chebyshev polynomials U(n,x).
!
!  First terms:
!
!    N/K     0     1      2      3       4     5      6    7      8    9   10
!
!     0      1
!     1      0     2
!     2     -1     0      4
!     3      0    -4      0      8
!     4      1     0    -12      0      16
!     5      0     6      0    -32       0    32
!     6     -1     0     24      0     -80     0     64
!     7      0    -8      0     80       0  -192      0   128
!
!  Recursion:
!
!    U(0,X) = 1,
!    U(1,X) = 2*X,
!    U(N,X) = 2 * X * U(N-1,X) - U(N-2,X)
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Output, real ( kind = rk ) C(0:N,0:N), the coefficients.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) c(0:n,0:n)
  integer i

  if ( n < 0 ) then
    return
  end if

  c(0:n,0:n) = 0.0D+00

  c(0,0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  c(1,1) = 2.0D+00
 
  do i = 2, n
    c(i,0)     =                        - c(i-2,0)
    c(i,1:i-2) = 2.0D+00 * c(i-1,0:i-3) - c(i-2,1:i-2)
    c(i,  i-1) = 2.0D+00 * c(i-1,  i-2)
    c(i,  i  ) = 2.0D+00 * c(i-1,  i-1)
  end do
 
  return
end
subroutine u_polynomial_plot ( n_num, n_val, output_filename )

!*****************************************************************************80
!
!! U_POLYNOMIAL_PLOT plots Chebyshev polynomials U(n,x).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N_NUM, the number of polynomials to be plotted.
!
!    Input, integer N_VAL(N_NUM), the degrees of 1 or more 
!    Chebyshev polynomials to be plotted together.
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the name into which the 
!    graphics information is to be stored.  Note that the PNG format will 
!    be used.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 501
  integer n_num

  real ( kind = rk ) a
  real ( kind = rk ) b
  integer column
  character ( len = 255 ) command_filename
  integer command_unit
  character ( len = 255 ) data_filename
  integer data_unit
  integer i
  integer i4vec_max
  integer j
  integer n
  integer n_max
  integer n_val(n_num)
  character ( len = * ) output_filename
  real ( kind = rk ), allocatable :: v(:,:)
  real ( kind = rk ) x(m)

  a = -1.0D+00
  b = +1.0D+00

  call r8vec_linspace ( m, a, b, x )
!
!  Compute all the data.
!
  n_max = i4vec_max ( n_num, n_val )
  allocate ( v(m,0:n_max) )
  call u_polynomial ( m, n_max, x, v )
!
!  Create the data file.
!
  data_filename = 'u_polynomial_data.txt'
  call get_unit ( data_unit )
  open ( unit = data_unit, file = data_filename, status = 'replace' )
  do i = 1, m
    write ( data_unit, '(2x,g14.6)', advance = 'no' ) x(i)
    do j = 1, n_num
      n = n_val(j)
      write ( data_unit, '(2x,g14.6)', advance = 'no' ) v(i,n)
    end do
    write ( data_unit, '(a)' ) ''
  end do
  close ( unit = data_unit )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Created graphics data file "' // trim ( data_filename ) // '".'
!
!  Plot the selected data.
!
  command_filename = 'u_polynomial_commands.txt'
  call get_unit ( command_unit )
  open ( unit = command_unit, file = command_filename, status = 'replace' )

  write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) 'set nokey'
  write ( command_unit, '(a)' ) &
    'set output "' // trim ( output_filename ) // '"'
  write ( command_unit, '(a)' ) 'set xlabel "<---X--->"'
  write ( command_unit, '(a)' ) 'set ylabel "<---U(n,x)--->"'
  write ( command_unit, '(a)' ) &
    'set title "Chebyshev Polynomials U(n,x)"'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set style data lines'
  do j = 1, n_num
    column = n_val(j) + 1
    if ( j == 1 ) then
      write ( command_unit, '(a)', advance = 'no' ) 'plot '
    else
      write ( command_unit, '(a)', advance = 'no' ) '     '
    end if
    write ( command_unit, '(a,i2,a)', advance = 'no' ) &
      '"' // trim ( data_filename ) // &
      '" using 1:', column, ' lw 3 linecolor rgb "red"'
    if ( j < n_num ) then
      write ( command_unit, '(a)' ) ', \'
    else
      write ( command_unit, '(a)' ) ''
    end if
  end do

  close ( unit = command_unit )
  write ( *, '(a)' ) &
    '  Created graphics command file "' // trim ( command_filename ) // '".'

  deallocate ( v )

  return
end
function u_polynomial_value ( n, x )

!*****************************************************************************80
!
!! U_POLYNOMIAL_VALUE: returns the single value U(n,x).
!
!  Discussion:
!
!    In cases where calling U_POLYNOMIAL is inconvenient, because it returns
!    a vector of values for multiple arguments X, this simpler interface
!    may be appropriate.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the polynomial.
!
!    Input, real ( kind = rk ) X, the argument of the polynomial.
!
!    Output, real ( kind = rk ) U_POLYNOMIAL_VALUE, the value of U(n,x).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n
  real ( kind = rk ) u_polynomial_value
  real ( kind = rk ) value
  real ( kind = rk ), allocatable :: vec(:)
  real ( kind = rk ) x
  real ( kind = rk ) x_vec(1)

  if ( n < 0 ) then

    value = 0.0D+00

  else

    m = 1
    allocate ( vec(0:n) )

    x_vec(1) = x
    call u_polynomial ( m, n, x_vec, vec )

    value = vec(n)
    deallocate ( vec )

  end if

  u_polynomial_value = value

  return
end
subroutine u_polynomial_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! U_POLYNOMIAL_VALUES returns values of Chebyshev polynomials U(n,x).
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      ChebyshevU[n,x]
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer N, the order of the function.
!
!    Output, real ( kind = rk ) X, the point where the function is evaluated.
!
!    Output, real ( kind = rk ) FX, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 14

  real ( kind = rk ) fx
  real ( kind = rk ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.0000000000000000D+00, &
     0.1000000000000000D+01, &
     0.1600000000000000D+01, &
     0.1560000000000000D+01, &
     0.8960000000000000D+00, &
    -0.1264000000000000D+00, &
    -0.1098240000000000D+01, &
    -0.1630784000000000D+01, &
    -0.1511014400000000D+01, &
    -0.7868390400000000D+00, &
     0.2520719360000000D+00, &
     0.1190154137600000D+01, &
     0.1652174684160000D+01, &
     0.1453325357056000D+01 /)
  integer n
  integer n_data
  integer, save, dimension ( n_max ) :: n_vec = (/ &
    -1, &
     0,  1,  2, &
     3,  4,  5, &
     6,  7,  8, &
     9, 10, 11, &
    12 /)
  real ( kind = rk ) x
  real ( kind = rk ), save, dimension ( n_max ) :: x_vec = (/ &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine u_polynomial_zeros ( n, z )

!*****************************************************************************80
!
!! U_POLYNOMIAL_ZEROS returns zeroes of Chebyshev polynomials U(n,x).
!
!  Discussion:
!
!    The I-th zero of U(N,X) is cos((I-1)*PI/(N-1)), I = 1 to N
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the polynomial.
!
!    Output, real ( kind = rk ) Z(N), the zeroes of U(N,X).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) angle
  integer i
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) z(n)

  do i = 1, n
    angle = real ( i, kind = rk ) * r8_pi / real ( n + 1, kind = rk )
    z(i) = cos ( angle )
  end do

  return
end
subroutine u_quadrature_rule ( n, t, w )

!*****************************************************************************80
!
!! U_QUADRATURE_RULE: quadrature rule for U(n,x).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the rule.
!
!    Output, real ( kind = rk ) T(N), W(N), the points and weights of the rule.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) bj(n)
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) t(n)
  real ( kind = rk ) w(n)

  t(1:n) = 0.0D+00

  bj(1:n) = 0.5D+00

  w(1) = sqrt ( r8_pi / 2.0D+00 )
  w(2:n) = 0.0D+00

  call imtqlx ( n, t, bj, w )

  w(1:n) = w(1:n) ** 2

  return
end
function uu_product ( i, j, x )

!*****************************************************************************80
!
!! UU_PRODUCT: evaluate U(i,x)*U(j,x)
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, J, the indices.
!
!    Input, real ( kind = rk ) X, the argument.
!
!    Output, real ( kind = rk ) UU_PRODUCT, the value.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer j
  integer k
  real ( kind = rk ) u_polynomial_value
  real ( kind = rk ) uu_product
  real ( kind = rk ) value
  real ( kind = rk ) x

  value = 0.0D+00
  do k = abs ( i - j ), i + j, 2
    value = value + u_polynomial_value ( k, x )
  end do

  uu_product = value

  return
end
function uu_product_integral ( i, j )

!*****************************************************************************80
!
!! UU_PRODUCT_INTEGRAL: integral (-1<=x<=1) U(i,x)*U(j,x)*sqrt(1-x^2) dx
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, J, the polynomial indices.
!    0 <= I, J.
!
!    Output, real ( kind = rk ) UU_PRODUCT_INTEGRAL, the integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer j
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) uu_product_integral
  real ( kind = rk ) value

  if ( i < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'UU_PRODUCT_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  0 <= I is required.'
    stop 1
  end if

  if ( j < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'UU_PRODUCT_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  0 <= J is required.'
    stop 1
  end if

  if ( i /= j ) then
    value = 0.0D+00
  else
    value = r8_pi / 2.0D+00
  end if

  uu_product_integral = value

  return
end
subroutine v_mass_matrix ( n, a )

!*****************************************************************************80
!
!! V_MASS_MATRIX computes the mass matrix for the Chebyshev V polynomial.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 July 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(0:n,0:n)
  integer i
  real ( kind = rk ), allocatable :: phi(:,:)
  real ( kind = rk ), allocatable :: phiw(:,:)
  real ( kind = rk ), allocatable :: w(:)
  real ( kind = rk ), allocatable :: x(:)

  allocate ( x(0:n) )
  allocate ( w(0:n) )

  call v_quadrature_rule ( n + 1, x, w )

  allocate ( phi(0:n,0:n) )

  call v_polynomial ( n + 1, n, x, phi )

  allocate ( phiw(0:n,0:n) )

  do i = 0, n
    phiw(0:n,i) = w(i) * phi(i,0:n)
  end do

  a(0:n,0:n) = matmul ( phiw(0:n,0:n), phi(0:n,0:n) )

  deallocate ( phi )
  deallocate ( phiw )
  deallocate ( w )
  deallocate ( x )

  return
end
function v_moment ( e )

!*****************************************************************************80
!
!! V_MOMENT: integral ( -1 <= x <= +1 ) x^e sqrt(1+x) / sqrt(1-x) dx.
!
!  Discussion:
!
!     E    V_MOMENT
!    --    --------------
!     0      pi
!     1      pi / 2
!     2      pi / 2
!     3    3 pi / 8
!     4    3 pi / 8
!     5    5 pi / 16
!     6    5 pi / 16
!     7   35 pi / 128
!     8   35 pi / 128
!     9   63 pi / 256
!    10   63 pi / 256
!    11  231 pi / 1024
!    12  231 pi / 1024
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer E, the exponent of X.
!    0 <= E.
!
!    Output, real ( kind = rk ) V_MOMENT, the value of the integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) f1
  real ( kind = rk ) f2
  real ( kind = rk ) f3
  real ( kind = rk ) f4
  real ( kind = rk ) f5
  real ( kind = rk ) f6
  real ( kind = rk ) f7
  real ( kind = rk ) f8
  integer e
  real ( kind = rk ) r8_e
  real ( kind = rk ) r8_factorial
  real ( kind = rk ) r8_hyper_2f1
  real ( kind = rk ) r8_mop
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) v_moment
  real ( kind = rk ) value

  r8_e = real ( e, kind = rk )

  f1 = 1.0D+00 / gamma ( 1.5D+00 + r8_e )
  f2 = r8_mop ( e )
  f3 = r8_pi * gamma ( 1.5D+00 + r8_e )
  f4 = 2.0D+00 * r8_hyper_2f1 ( 0.5D+00, -r8_e, 1.0D+00, 2.0D+00 )
  f5 = ( -1.0D+00 + r8_mop ( e ) ) &
    * r8_hyper_2f1 ( 0.5D+00, -r8_e, 2.0D+00, 2.0D+00 )
  f6 = sqrt ( r8_pi ) * r8_factorial ( e )
  f7 = ( -1.0D+00 + r8_mop ( e ) ) &
    * r8_hyper_2f1 ( -0.5D+00, 1.0D+00 + r8_e, 1.5D+00 + r8_e, - 1.0D+00 )
  f8 = 2.0D+00 &
    * r8_hyper_2f1 ( 0.5D+00, 1.0D+00 + r8_e, 1.5D+00 + r8_e, -1.0D+00 )

  value = f1 * f2 * ( f3 * ( f4 + f5 ) - f6 * ( f7 + f8 ) )

  v_moment = value

  return
end
subroutine v_polynomial ( m, n, x, v )

!*****************************************************************************80
!
!! V_POLYNOMIAL evaluates Chebyshev polynomials V(n,x).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer M, the number of evaluation points.
!
!    Input, integer N, the highest polynomial to compute.
!
!    Input, real ( kind = rk ) X(M), the evaluation points.
!
!    Output, real ( kind = rk ) V(M,0:N), the values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer i
  real ( kind = rk ) v(m,0:n)
  real ( kind = rk ) x(m)

  if ( n < 0 ) then
    return
  end if

  v(1:m,0) = 1.0D+00

  if ( n < 1 ) then
    return
  end if

  v(1:m,1) = 2.0D+00 * x(1:m) - 1.0D+00

  do i = 2, n
    v(1:m,i) = 2.0D+00 * x(1:m) * v(1:m,i-1) - v(1:m,i-2)
  end do
 
  return
end
subroutine v_polynomial_01_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! V_POLYNOMIAL_01_VALUES: values of shifted Chebyshev polynomials V01(n,x).
!
!  Discussion:
!
!    The shifted Chebyshev polynomial V01(n,x) = V(n,2*x-1).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer N, the order of the function.
!
!    Output, real ( kind = rk ) X, the point where the function is evaluated.
!
!    Output, real ( kind = rk ) FX, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 25

  real ( kind = rk ) fx
  real ( kind = rk ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.0000000000000000D+00, &
     1.0000000000000000D+00, &
     0.4000000000000000D+00, &
    -0.4400000000000000D+00, &
    -1.0160000000000000D+00, &
    -0.9824000000000000D+00, &
    -0.3593600000000000D+00, &
     0.4792960000000000D+00, &
     1.0303744000000000D+00, &
     0.9632281600000000D+00, &
     0.3181450240000000D+00, &
    -0.5178251264000000D+00, &
    -1.0431002009600000D+00, &
    -0.9425151549440000D+00, &
    -15.000000000000000D+00, &
     3.1417984000000000D+00, &
    -1.3912448000000000D+00, &
    -1.2177792000000000D+00, &
     1.1837056000000000D+00, &
     1.0000000000000000D+00, &
    -0.8558976000000000D+00, &
    -0.8905088000000000D+00, &
     0.8752768000000000D+00, &
     0.1197696000000000D+00, &
     1.0000000000000000D+00 /)
  integer n
  integer n_data
  integer, save, dimension ( n_max ) :: n_vec = (/ &
    -1, &
     0,  1,  2, &
     3,  4,  5, &
     6,  7,  8, &
     9, 10, 11, &
    12,  7,  7, &
     7,  7,  7, &
     7,  7,  7, &
     7,  7,  7 /)
  real ( kind = rk ) x
  real ( kind = rk ), save, dimension ( n_max ) :: x_vec = (/ &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.00D+00, &
    0.10D+00, &
    0.20D+00, &
    0.30D+00, &
    0.40D+00, &
    0.50D+00, &
    0.60D+00, &
    0.70D+00, &
    0.80D+00, &
    0.90D+00, &
    1.00D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine v_polynomial_ab ( a, b, m, n, xab, v )

!*****************************************************************************80
!
!! V_POLYNOMIAL_AB: evaluates Chebyshev polynomials VAB(n,x) in [A,B].
!
!  Discussion:
!
!    VAB(n,x) = V(n,(2*x-a-b)/(b-a))
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, the domain of definition.
!
!    Input, integer M, the number of evaluation points.
!
!    Input, integer N, the highest polynomial to compute.
!
!    Input, real ( kind = rk ) XAB(M), the evaluation points.
!    It must be the case that A <= XAB(*) <= B.
!
!    Output, real ( kind = rk ) V(M,N+1), the values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) v(1:m,0:n)
  real ( kind = rk ) x(1:m)
  real ( kind = rk ) xab(1:m)

  x(1:m) = ( 2.0D+00 * xab(1:m) - a - b ) / ( b - a )

  call v_polynomial ( m, n, x, v )
 
  return
end
function v_polynomial_ab_value ( a, b, n, xab )

!*****************************************************************************80
!
!! V_POLYNOMIAL_AB_VALUE: evaluates Chebyshev polynomials VAB(n,x) in [A,B].
!
!  Discussion:
!
!    VAB(n,x) = V(n,(2*x-a-b)/(b-a))
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, the domain of definition.
!
!    Input, integer N, the highest polynomial to compute.
!
!    Input, real ( kind = rk ) XAB, the evaluation point.
!    It must be the case that A <= XAB <= B.
!
!    Output, real ( kind = rk ) V_POLYNOMIAL_AB_VALUE, the value.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) v_polynomial_ab_value
  real ( kind = rk ) v_polynomial_value
  real ( kind = rk ) value
  real ( kind = rk ) x
  real ( kind = rk ) xab

  x = ( 2.0D+00 * xab - a - b ) / ( b - a )

  value = v_polynomial_value ( n, x )
 
  v_polynomial_ab_value = value

  return
end
subroutine v_polynomial_coefficients ( n, c )

!*****************************************************************************80
!
!! V_POLYNOMIAL_COEFFICIENTS: coefficients of Chebyshev polynomials V(n,x).
!
!  First terms:
!
!    N/K     0     1      2      3       4     5      6    7      8    9   10
!
!     0      1
!     1     -1     2
!     2     -1    -2      4
!     3      1    -4     -4      8
!     4      1    +4    -12     -8      16
!     5     -1     6    +12    -32     -16    32
!     6     -1    -6     24    +32     -80   -32     64
!     7     +1    -8    -24     80     +80  -192    -64   128
!
!  Recursion:
!
!    V(0,X) = 1,
!    V(1,X) = 2 * X - 1,
!    V(N,X) = 2 * X * V(N-1,X) - V(N-2,X)
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Output, real ( kind = rk ) C(0:N,0:N), the coefficients.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) c(0:n,0:n)
  integer i

  if ( n < 0 ) then
    return
  end if

  c(0:n,0:n) = 0.0D+00

  c(0,0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  c(1,0) = -1.0D+00
  c(1,1) =  2.0D+00
 
  do i = 2, n
    c(i,0)     =                        - c(i-2,0)
    c(i,1:i-2) = 2.0D+00 * c(i-1,0:i-3) - c(i-2,1:i-2)
    c(i,  i-1) = 2.0D+00 * c(i-1,  i-2)
    c(i,  i  ) = 2.0D+00 * c(i-1,  i-1)
  end do
 
  return
end
subroutine v_polynomial_plot ( n_num, n_val, output_filename )

!*****************************************************************************80
!
!! V_POLYNOMIAL_PLOT plots Chebyshev polynomials V(n,x).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N_NUM, the number of polynomials to be plotted.
!
!    Input, integer N_VAL(N_NUM), the degrees of 1 or more 
!    Chebyshev polynomials to be plotted together.
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the name into which the 
!    graphics information is to be stored.  Note that the PNG format will 
!    be used.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 501
  integer n_num

  real ( kind = rk ) a
  real ( kind = rk ) b
  integer column
  character ( len = 255 ) command_filename
  integer command_unit
  character ( len = 255 ) data_filename
  integer data_unit
  integer i
  integer i4vec_max
  integer j
  integer n
  integer n_max
  integer n_val(n_num)
  character ( len = * ) output_filename
  real ( kind = rk ), allocatable :: v(:,:)
  real ( kind = rk ) x(m)

  a = -1.0D+00
  b = +1.0D+00

  call r8vec_linspace ( m, a, b, x )
!
!  Compute all the data.
!
  n_max = i4vec_max ( n_num, n_val )
  allocate ( v(m,0:n_max) )
  call v_polynomial ( m, n_max, x, v )
!
!  Create the data file.
!
  data_filename = 'v_polynomial_data.txt'
  call get_unit ( data_unit )
  open ( unit = data_unit, file = data_filename, status = 'replace' )
  do i = 1, m
    write ( data_unit, '(2x,g14.6)', advance = 'no' ) x(i)
    do j = 1, n_num
      n = n_val(j)
      write ( data_unit, '(2x,g14.6)', advance = 'no' ) v(i,n)
    end do
    write ( data_unit, '(a)' ) ''
  end do
  close ( unit = data_unit )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Created graphics data file "' // trim ( data_filename ) // '".'
!
!  Plot the selected data.
!
  command_filename = 'v_polynomial_commands.txt'
  call get_unit ( command_unit )
  open ( unit = command_unit, file = command_filename, status = 'replace' )

  write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) 'set nokey'
  write ( command_unit, '(a)' ) &
    'set output "' // trim ( output_filename ) // '"'
  write ( command_unit, '(a)' ) 'set xlabel "<---X--->"'
  write ( command_unit, '(a)' ) 'set ylabel "<---T(n,x)--->"'
  write ( command_unit, '(a)' ) &
    'set title "Chebyshev Polynomials V(n,x)"'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set style data lines'
  do j = 1, n_num
    column = n_val(j) + 1
    if ( j == 1 ) then
      write ( command_unit, '(a)', advance = 'no' ) 'plot '
    else
      write ( command_unit, '(a)', advance = 'no' ) '     '
    end if
    write ( command_unit, '(a,i2,a)', advance = 'no' ) &
      '"' // trim ( data_filename ) // &
      '" using 1:', column, ' lw 3 linecolor rgb "red"'
    if ( j < n_num ) then
      write ( command_unit, '(a)' ) ', \'
    else
      write ( command_unit, '(a)' ) ''
    end if
  end do

  close ( unit = command_unit )
  write ( *, '(a)' ) &
    '  Created graphics command file "' // trim ( command_filename ) // '".'

  deallocate ( v )

  return
end
function v_polynomial_value ( n, x )

!*****************************************************************************80
!
!! V_POLYNOMIAL_VALUE: returns the single value V(n,x).
!
!  Discussion:
!
!    In cases where calling V_POLYNOMIAL is inconvenient, because it returns
!    a vector of values for multiple arguments X, this simpler interface
!    may be appropriate.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the polynomial.
!
!    Input, real ( kind = rk ) X, the argument of the polynomial.
!
!    Output, real ( kind = rk ) V_POLYNOMIAL_VALUE, the value of V(n,x).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n
  real ( kind = rk ) v_polynomial_value
  real ( kind = rk ) value
  real ( kind = rk ), allocatable :: vec(:)
  real ( kind = rk ) x
  real ( kind = rk ) x_vec(1)

  if ( n < 0 ) then

    value = 0.0D+00

  else

    m = 1
    allocate ( vec(0:n) )

    x_vec(1) = x
    call v_polynomial ( m, n, x_vec, vec )

    value = vec(n)
    deallocate ( vec )

  end if

  v_polynomial_value = value

  return
end
subroutine v_polynomial_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! V_POLYNOMIAL_VALUES returns values of Chebyshev polynomials V(n,x).
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      u = Sqrt[(x+1)/2],
!      ChebyshevT[2*n+1,u] / u
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer N, the order of the function.
!
!    Output, real ( kind = rk ) X, the point where the function is evaluated.
!
!    Output, real ( kind = rk ) FX, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 14

  real ( kind = rk ) fx
  real ( kind = rk ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.0000000000000000D+00, &
     1.0000000000000000D+00, &
     0.6000000000000000D+00, &
    -0.0400000000000000D+00, &
    -0.6640000000000000D+00, &
    -1.0224000000000000D+00, &
    -0.9718400000000000D+00, &
    -0.5325440000000000D+00, &
     0.1197696000000000D+00, &
     0.7241753600000000D+00, &
     1.0389109760000000D+00, &
     0.9380822016000000D+00, &
     0.4620205465600000D+00, &
    -0.1988493271040000D+00 /)
  integer n
  integer n_data
  integer, save, dimension ( n_max ) :: n_vec = (/ &
    -1, &
     0,  1,  2, &
     3,  4,  5, &
     6,  7,  8, &
     9, 10, 11, &
    12 /)
  real ( kind = rk ) x
  real ( kind = rk ), save, dimension ( n_max ) :: x_vec = (/ &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine v_polynomial_zeros ( n, z )

!*****************************************************************************80
!
!! V_POLYNOMIAL_ZEROS returns zeroes of Chebyshev polynomials V(n,x).
!
!  Discussion:
!
!    The I-th zero of U(N,X) is cos((I-1/2)*PI/(N+1/2)), I = 1 to N
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    25 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the polynomial.
!
!    Output, real ( kind = rk ) Z(N), the zeroes.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) angle
  integer i
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) z(n)

  do i = 1, n
    angle = real ( 2 * n - 2 * i + 1, kind = rk ) * r8_pi &
      / real ( 2 * n + 1, kind = rk )
    z(i) = cos ( angle )
  end do

  return
end
subroutine v_quadrature_rule ( n, t, w )

!*****************************************************************************80
!
!! V_QUADRATURE_RULE: quadrature rule for V(n,x).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the rule.
!
!    Output, real ( kind = rk ) T(N), W(N), the points and weights of the rule.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) bj(n)
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) t(n)
  real ( kind = rk ) w(n)

  t(1:n) = 0.0D+00
  t(1) = + 0.5D+00

  bj(1:n) = 0.5D+00

  w(1) = sqrt ( r8_pi )
  w(2:n) = 0.0D+00

  call imtqlx ( n, t, bj, w )

  w(1:n) = w(1:n) ** 2

  return
end
function vv_product_integral ( i, j )

!*****************************************************************************80
!
!! VV_PRODUCT_INTEGRAL: int (-1<x<1) V(i,x)*V(j,x)*sqrt(1+x)/sqrt(1-x) dx
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, J, the polynomial indices.
!    0 <= I, J.
!
!    Output, real ( kind = rk ) VV_PRODUCT_INTEGRAL, the integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer j
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) vv_product_integral
  real ( kind = rk ) value

  if ( i < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VV_PRODUCT_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  0 <= I is required.'
    stop 1
  end if

  if ( j < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VV_PRODUCT_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  0 <= J is required.'
    stop 1
  end if

  if ( i /= j ) then
    value = 0.0D+00
  else
    value = r8_pi
  end if

  vv_product_integral = value

  return
end
subroutine w_mass_matrix ( n, a )

!*****************************************************************************80
!
!! W_MASS_MATRIX computes the mass matrix for the Chebyshev W polynomial.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 July 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(0:n,0:n)
  integer i
  real ( kind = rk ), allocatable :: phi(:,:)
  real ( kind = rk ), allocatable :: phiw(:,:)
  real ( kind = rk ), allocatable :: w(:)
  real ( kind = rk ), allocatable :: x(:)

  allocate ( x(0:n) )
  allocate ( w(0:n) )

  call w_quadrature_rule ( n + 1, x, w )

  allocate ( phi(0:n,0:n) )

  call w_polynomial ( n + 1, n, x, phi )

  allocate ( phiw(0:n,0:n) )

  do i = 0, n
    phiw(0:n,i) = w(i) * phi(i,0:n)
  end do

  a(0:n,0:n) = matmul ( phiw(0:n,0:n), phi(0:n,0:n) )

  deallocate ( phi )
  deallocate ( phiw )
  deallocate ( w )
  deallocate ( x )

  return
end
function w_moment ( e )

!*****************************************************************************80
!
!! W_MOMENT: integral ( -1 <= x <= +1 ) x^e sqrt(1-x) / sqrt(1+x) dx.
!
!  Discussion:
!
!     E    W_MOMENT
!    --    --------------
!     0        pi
!     1  -     pi / 2
!     2        pi / 2
!     3  -   3 pi / 8
!     4      3 pi / 8
!     5  -   5 pi / 16
!     6      5 pi / 16
!     7  -  35 pi / 128
!     8     35 pi / 128
!     9  -  63 pi / 256
!    10     63 pi / 256
!    11  - 231 pi / 1024
!    12    231 pi / 1024
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer E, the exponent of X.
!    0 <= E.
!
!    Output, real ( kind = rk ) W_MOMENT, the value of the integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) f1
  real ( kind = rk ) f2
  real ( kind = rk ) f3
  real ( kind = rk ) f4
  real ( kind = rk ) f5
  real ( kind = rk ) f6
  real ( kind = rk ) f7
  real ( kind = rk ) f8
  integer e
  real ( kind = rk ) r8_e
  real ( kind = rk ) r8_factorial
  real ( kind = rk ) r8_hyper_2f1
  real ( kind = rk ) r8_mop
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) w_moment
  real ( kind = rk ) value

  r8_e = real ( e, kind = rk )

  f1 = 1.0D+00 / gamma ( 1.5D+00 + r8_e )
  f2 = r8_mop ( e )
  f3 = r8_pi * gamma ( 1.5D+00 + r8_e )
  f4 = 2.0D+00 * r8_mop ( e ) &
    * r8_hyper_2f1 ( 0.5D+00, -r8_e, 1.0D+00, 2.0D+00 )
  f5 = ( -1.0D+00 + r8_mop ( e ) ) &
    * r8_hyper_2f1 ( 0.5D+00, -r8_e, 2.0D+00, 2.0D+00 )
  f6 = sqrt ( r8_pi ) * r8_factorial ( e )
  f7 = ( -1.0D+00 + r8_mop ( e ) ) &
    * r8_hyper_2f1 ( -0.5D+00, 1.0D+00 + r8_e, 1.5D+00 + r8_e, - 1.0D+00 )
  f8 = 2.0D+00 * r8_mop ( e ) &
    * r8_hyper_2f1 ( 0.5D+00, 1.0D+00 + r8_e, 1.5D+00 + r8_e, -1.0D+00 )

  value = f1 * f2 * ( f3 * ( f4 - f5 ) + f6 * ( f7 - f8 ) )

  w_moment = value

  return
end
subroutine w_polynomial ( m, n, x, v )

!*****************************************************************************80
!
!! W_POLYNOMIAL evaluates Chebyshev polynomials W(n,x).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer M, the number of evaluation points.
!
!    Input, integer N, the highest polynomial to compute.
!
!    Input, real ( kind = rk ) X(M), the evaluation points.
!
!    Output, real ( kind = rk ) V(M,0:N), the values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer i
  real ( kind = rk ) v(m,0:n)
  real ( kind = rk ) x(m)

  if ( n < 0 ) then
    return
  end if

  v(1:m,0) = 1.0D+00

  if ( n < 1 ) then
    return
  end if

  v(1:m,1) = 2.0D+00 * x(1:m) + 1.0D+00

  do i = 2, n
    v(1:m,i) = 2.0D+00 * x(1:m) * v(1:m,i-1) - v(1:m,i-2)
  end do
 
  return
end
subroutine w_polynomial_01_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! W_POLYNOMIAL_01_VALUES: values of shifted Chebyshev polynomials W01(n,x).
!
!  Discussion:
!
!    The shifted Chebyshev polynomial W01(n,x) = W(n,2*x-1).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer N, the order of the function.
!
!    Output, real ( kind = rk ) X, the point where the function is evaluated.
!
!    Output, real ( kind = rk ) FX, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 25

  real ( kind = rk ) fx
  real ( kind = rk ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.000000000000000D+00, &
     1.000000000000000D+00, &
     2.400000000000000D+00, &
     2.360000000000000D+00, &
     0.904000000000000D+00, &
    -1.094400000000000D+00, &
    -2.436160000000000D+00, &
    -2.316224000000000D+00, &
    -0.806553600000000D+00, &
     1.187048960000000D+00, &
     2.468422144000000D+00, &
     2.268742041600000D+00, &
     0.707816714240000D+00, &
    -1.277798641664000D+00, &
    -1.000000000000000D+00, &
    -0.119769600000000D+00, &
    -0.875276800000000D+00, &
     0.890508800000000D+00, &
     0.855897600000000D+00, &
    -1.000000000000000D+00, &
    -1.183705600000000D+00, &
     1.217779200000000D+00, &
     1.391244800000000D+00, &
    -3.141798400000000D+00, &
     15.00000000000000D+00 /)
  integer n
  integer n_data
  integer, save, dimension ( n_max ) :: n_vec = (/ &
    -1, &
     0,  1,  2, &
     3,  4,  5, &
     6,  7,  8, &
     9, 10, 11, &
    12,  7,  7, &
     7,  7,  7, &
     7,  7,  7, &
     7,  7,  7 /)
  real ( kind = rk ) x
  real ( kind = rk ), save, dimension ( n_max ) :: x_vec = (/ &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.00D+00, &
    0.10D+00, &
    0.20D+00, &
    0.30D+00, &
    0.40D+00, &
    0.50D+00, &
    0.60D+00, &
    0.70D+00, &
    0.80D+00, &
    0.90D+00, &
    1.00D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine w_polynomial_ab ( a, b, m, n, xab, v )

!*****************************************************************************80
!
!! W_POLYNOMIAL_AB: evaluates Chebyshev polynomials WAB(n,x) in [A,B].
!
!  Discussion:
!
!    WAB(n,x) = W(n,(2*x-a-b)/(b-a))
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, the domain of definition.
!
!    Input, integer M, the number of evaluation points.
!
!    Input, integer N, the highest polynomial to compute.
!
!    Input, real ( kind = rk ) XAB(M), the evaluation points.
!    It must be the case that A <= XAB(*) <= B.
!
!    Output, real ( kind = rk ) V(M,N+1), the values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) v(1:m,0:n)
  real ( kind = rk ) x(1:m)
  real ( kind = rk ) xab(1:m)

  x(1:m) = ( 2.0D+00 * xab(1:m) - a - b ) / ( b - a )

  call w_polynomial ( m, n, x, v )
 
  return
end
function w_polynomial_ab_value ( a, b, n, xab )

!*****************************************************************************80
!
!! W_POLYNOMIAL_AB_VALUE: evaluates Chebyshev polynomials WAB(n,x) in [A,B].
!
!  Discussion:
!
!    WAB(n,x) = W(n,(2*x-a-b)/(b-a))
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, the domain of definition.
!
!    Input, integer N, the highest polynomial to compute.
!
!    Input, real ( kind = rk ) XAB, the evaluation point.
!    It must be the case that A <= XAB <= B.
!
!    Output, real ( kind = rk ) W_POLYNOMIAL_AB_VALUE, the value.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) value
  real ( kind = rk ) w_polynomial_ab_value
  real ( kind = rk ) w_polynomial_value
  real ( kind = rk ) x
  real ( kind = rk ) xab

  x = ( 2.0D+00 * xab - a - b ) / ( b - a )

  value = w_polynomial_value ( n, x )
 
  w_polynomial_ab_value = value

  return
end
subroutine w_polynomial_coefficients ( n, c )

!*****************************************************************************80
!
!! W_POLYNOMIAL_COEFFICIENTS: coefficients of Chebyshev polynomials W(n,x).
!
!  First terms:
!
!    N/K     0     1      2      3       4     5      6    7      8    9   10
!
!     0      1
!     1      1     2
!     2     -1     2      4
!     3     -1    -4      4      8
!     4      1    -4    -12      8      16
!     5      1     6    -12    -32     +16    32
!     6     -1     6     24    -32     -80    32     64
!     7     -1    -8    +24    +80     -80  -192     64   128
!
!  Recursion:
!
!    W(0,X) = 1,
!    W(1,X) = 2 * X + 1,
!    W(N,X) = 2 * X * W(N-1,X) - W(N-2,X)
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Output, real ( kind = rk ) C(0:N,0:N), the coefficients.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) c(0:n,0:n)
  integer i

  if ( n < 0 ) then
    return
  end if

  c(0:n,0:n) = 0.0D+00

  c(0,0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  c(1,0) = + 1.0D+00
  c(1,1) = + 2.0D+00
 
  do i = 2, n
    c(i,0)     =                        - c(i-2,0)
    c(i,1:i-2) = 2.0D+00 * c(i-1,0:i-3) - c(i-2,1:i-2)
    c(i,  i-1) = 2.0D+00 * c(i-1,  i-2)
    c(i,  i  ) = 2.0D+00 * c(i-1,  i-1)
  end do
 
  return
end
subroutine w_polynomial_plot ( n_num, n_val, output_filename )

!*****************************************************************************80
!
!! W_POLYNOMIAL_PLOT plots Chebyshev polynomials W(n,x).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N_NUM, the number of polynomials to be plotted.
!
!    Input, integer N_VAL(N_NUM), the degrees of 1 or more 
!    Chebyshev polynomials to be plotted together.
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the name into which the 
!    graphics information is to be stored.  Note that the PNG format will 
!    be used.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 501
  integer n_num

  real ( kind = rk ) a
  real ( kind = rk ) b
  integer column
  character ( len = 255 ) command_filename
  integer command_unit
  character ( len = 255 ) data_filename
  integer data_unit
  integer i
  integer i4vec_max
  integer j
  integer n
  integer n_max
  integer n_val(n_num)
  character ( len = * ) output_filename
  real ( kind = rk ), allocatable :: v(:,:)
  real ( kind = rk ) x(m)

  a = -1.0D+00
  b = +1.0D+00

  call r8vec_linspace ( m, a, b, x )
!
!  Compute all the data.
!
  n_max = i4vec_max ( n_num, n_val )
  allocate ( v(m,0:n_max) )
  call w_polynomial ( m, n_max, x, v )
!
!  Create the data file.
!
  data_filename = 'w_polynomial_data.txt'
  call get_unit ( data_unit )
  open ( unit = data_unit, file = data_filename, status = 'replace' )
  do i = 1, m
    write ( data_unit, '(2x,g14.6)', advance = 'no' ) x(i)
    do j = 1, n_num
      n = n_val(j)
      write ( data_unit, '(2x,g14.6)', advance = 'no' ) v(i,n)
    end do
    write ( data_unit, '(a)' ) ''
  end do
  close ( unit = data_unit )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Created graphics data file "' // trim ( data_filename ) // '".'
!
!  Plot the selected data.
!
  command_filename = 'w_polynomial_commands.txt'
  call get_unit ( command_unit )
  open ( unit = command_unit, file = command_filename, status = 'replace' )

  write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) 'set nokey'
  write ( command_unit, '(a)' ) &
    'set output "' // trim ( output_filename ) // '"'
  write ( command_unit, '(a)' ) 'set xlabel "<---X--->"'
  write ( command_unit, '(a)' ) 'set ylabel "<---T(n,x)--->"'
  write ( command_unit, '(a)' ) &
    'set title "Chebyshev Polynomials W(n,x)"'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set style data lines'
  do j = 1, n_num
    column = n_val(j) + 1
    if ( j == 1 ) then
      write ( command_unit, '(a)', advance = 'no' ) 'plot '
    else
      write ( command_unit, '(a)', advance = 'no' ) '     '
    end if
    write ( command_unit, '(a,i2,a)', advance = 'no' ) &
      '"' // trim ( data_filename ) // &
      '" using 1:', column, ' lw 3 linecolor rgb "red"'
    if ( j < n_num ) then
      write ( command_unit, '(a)' ) ', \'
    else
      write ( command_unit, '(a)' ) ''
    end if
  end do

  close ( unit = command_unit )
  write ( *, '(a)' ) &
    '  Created graphics command file "' // trim ( command_filename ) // '".'

  deallocate ( v )

  return
end
function w_polynomial_value ( n, x )

!*****************************************************************************80
!
!! W_POLYNOMIAL_VALUE: returns the single value W(n,x).
!
!  Discussion:
!
!    In cases where calling W_POLYNOMIAL is inconvenient, because it returns
!    a vector of values for multiple arguments X, this simpler interface
!    may be appropriate.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the polynomial.
!
!    Input, real ( kind = rk ) X, the argument of the polynomial.
!
!    Output, real ( kind = rk ) W_POLYNOMIAL_VALUE, the value of T(n,x).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n
  real ( kind = rk ) value
  real ( kind = rk ), allocatable :: vec(:)
  real ( kind = rk ) w_polynomial_value
  real ( kind = rk ) x
  real ( kind = rk ) x_vec(1)

  if ( n < 0 ) then

    value = 0.0D+00

  else

    m = 1
    allocate ( vec(0:n) )

    x_vec(1) = x
    call w_polynomial ( m, n, x_vec, vec )

    value = vec(n)
    deallocate ( vec )

  end if

  w_polynomial_value = value

  return
end
subroutine w_polynomial_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! W_POLYNOMIAL_VALUES returns values of Chebyshev polynomials W(n,x).
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      u = Sqrt[(x+1)/2],
!      ChebyshevU[2*n,u]
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer N, the order of the function.
!
!    Output, real ( kind = rk ) X, the point where the function is evaluated.
!
!    Output, real ( kind = rk ) FX, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 14

  real ( kind = rk ) fx
  real ( kind = rk ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.000000000000000D+00, &
     1.000000000000000D+00, &
     2.600000000000000D+00, &
     3.160000000000000D+00, &
     2.456000000000000D+00, &
     0.769600000000000D+00, &
    -1.224640000000000D+00, &
    -2.729024000000000D+00, &
    -3.141798400000000D+00, &
    -2.297853440000000D+00, &
    -0.534767104000000D+00, &
     1.442226073600000D+00, &
     2.842328821760000D+00, &
     3.105500041216000D+00 /)
  integer n
  integer n_data
  integer, save, dimension ( n_max ) :: n_vec = (/ &
    -1, &
     0,  1,  2, &
     3,  4,  5, &
     6,  7,  8, &
     9, 10, 11, &
    12 /)
  real ( kind = rk ) x
  real ( kind = rk ), save, dimension ( n_max ) :: x_vec = (/ &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine w_polynomial_zeros ( n, z )

!*****************************************************************************80
!
!! W_POLYNOMIAL_ZEROS returns zeroes of Chebyshev polynomials W(n,x).
!
!  Discussion:
!
!    The I-th zero of U(N,X) is cos(I*PI/(N+1/2)), I = 1 to N
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    25 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the polynomial.
!
!    Output, real ( kind = rk ) Z(N), the zeroes.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) angle
  integer i
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) z(n)

  do i = 1, n
    angle = real ( 2 * ( n - i + 1 ), kind = rk ) * r8_pi &
      / real ( 2 * n + 1, kind = rk )
    z(i) = cos ( angle )
  end do

  return
end
subroutine w_quadrature_rule ( n, t, w )

!*****************************************************************************80
!
!! W_QUADRATURE_RULE: quadrature rule for W(n,x).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the rule.
!
!    Output, real ( kind = rk ) T(N), W(N), the points and weights of the rule.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) bj(n)
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) t(n)
  real ( kind = rk ) w(n)

  t(1:n) = 0.0D+00
  t(1) = - 0.5D+00

  bj(1:n) = 0.5D+00

  w(1) = sqrt ( r8_pi )
  w(2:n) = 0.0D+00

  call imtqlx ( n, t, bj, w )

  w(1:n) = w(1:n) ** 2

  return
end
function ww_product_integral ( i, j )

!*****************************************************************************80
!
!! WW_PRODUCT_INTEGRAL: int (-1<x<1) W(i,x)*W(j,x)*sqrt(1-x)/sqrt(1+x) dx
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, J, the polynomial indices.
!    0 <= I, J.
!
!    Output, real ( kind = rk ) WW_PRODUCT_INTEGRAL, the integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer j
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) value
  real ( kind = rk ) ww_product_integral

  if ( i < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WW_PRODUCT_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  0 <= I is required.'
    stop 1
  end if

  if ( j < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WW_PRODUCT_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  0 <= J is required.'
    stop 1
  end if

  if ( i /= j ) then
    value = 0.0D+00
  else
    value = r8_pi
  end if

  ww_product_integral = value

  return
end
