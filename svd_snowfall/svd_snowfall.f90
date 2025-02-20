subroutine ch_cap ( c )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, the character to capitalize.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character c
  integer itemp

  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end
function ch_eqi ( c1, c2 )

!*****************************************************************************80
!
!! CH_EQI is a case insensitive comparison of two characters for equality.
!
!  Example:
!
!    CH_EQI ( 'A', 'a' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C1, C2, the characters to compare.
!
!    Output, logical CH_EQI, the result of the comparison.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  logical ch_eqi
  character c1
  character c1_cap
  character c2
  character c2_cap

  c1_cap = c1
  c2_cap = c2

  call ch_cap ( c1_cap )
  call ch_cap ( c2_cap )

  if ( c1_cap == c2_cap ) then
    ch_eqi = .true.
  else
    ch_eqi = .false.
  end if

  return
end
subroutine ch_to_digit ( c, digit )

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the value of a base 10 digit.
!
!  Example:
!
!     C   DIGIT
!    ---  -----
!    '0'    0
!    '1'    1
!    ...  ...
!    '9'    9
!    ' '    0
!    'X'   -1
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the decimal digit, '0' through '9' or blank
!    are legal.
!
!    Output, integer DIGIT, the corresponding value.
!    If C was 'illegal', then DIGIT is -1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character c
  integer digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

    digit = ichar ( c ) - 48

  else if ( c == ' ' ) then

    digit = 0

  else

    digit = -1

  end if

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
!    This routine uses double precision real arithmetic.
!
!     DNRM2 ( X ) = sqrt ( X' * X )
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
!    Input, integer N, the number of entries in the vector.
!
!    Input, real ( kind = rk ) X(*), the vector whose norm is to be computed.
!
!    Input, integer INCX, the increment between successive 
!    entries of X.
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
    norm  = scale * sqrt ( ssq )
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

        do j = l+1, n
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

      do lls = l+1, mn+1

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
subroutine file_column_count ( input_filename, column_num )

!*****************************************************************************80
!
!! FILE_COLUMN_COUNT counts the number of columns in the first line of a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.
!
!    Most lines of the file is presumed to consist of COLUMN_NUM words,
!    separated by spaces.  There may also be some blank lines, and some
!    comment lines,
!    which have a "#" in column 1.
!
!    The routine tries to find the first non-comment non-blank line and
!    counts the number of words in that line.
!
!    If all lines are blanks or comments, it goes back and tries to analyze
!    a comment line.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the file.
!
!    Output, integer COLUMN_NUM, the number of columns in the file.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer column_num
  logical got_one
  character ( len = * ) input_filename
  integer input_status
  integer input_unit
  character ( len = 255 ) line
!
!  Open the file.
!
  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    form = 'formatted', access = 'sequential', iostat = input_status )

  if ( input_status /= 0 ) then
    column_num = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' &
      // trim ( input_filename ) // '" on unit ', input_unit
    return
  end if
!
!  Read one line, but skip blank lines and comment lines.
!
  got_one = .false.

  do

    read ( input_unit, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
      exit
    end if

    if ( len_trim ( line ) == 0 ) then
      cycle
    end if

    if ( line(1:1) == '#' ) then
      cycle
    end if

    got_one = .true.
    exit

  end do

  if ( .not. got_one ) then

    rewind ( input_unit )

    do

      read ( input_unit, '(a)', iostat = input_status ) line

      if ( input_status /= 0 ) then
        exit
      end if

      if ( len_trim ( line ) == 0 ) then
        cycle
      end if

      got_one = .true.
      exit

    end do

  end if

  close ( unit = input_unit )

  if ( .not. got_one ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Warning!'
    write ( *, '(a)' ) '  The file does not seem to contain any data.'
    column_num = -1
    return
  end if

  call s_word_count ( line, column_num )

  return
end
subroutine file_row_count ( input_filename, row_num )

!*****************************************************************************80
!
!! FILE_ROW_COUNT counts the number of row records in a file.
!
!  Discussion:
!
!    It does not count lines that are blank, or that begin with a
!    comment symbol '#'.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Output, integer ROW_NUM, the number of rows found.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer bad_num
  integer comment_num
  integer ierror
  character ( len = * ) input_filename
  integer input_status
  integer input_unit
  character ( len = 255 ) line
  integer record_num
  integer row_num

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = input_status )

  if ( input_status /= 0 ) then
    row_num = -1;
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_ROW_COUNT - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' // &
      trim ( input_filename ) // '" on unit ', input_unit
    stop
  end if

  comment_num = 0
  row_num = 0
  record_num = 0
  bad_num = 0

  do

    read ( input_unit, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
      ierror = record_num
      exit
    end if

    record_num = record_num + 1

    if ( line(1:1) == '#' ) then
      comment_num = comment_num + 1
      cycle
    end if

    if ( len_trim ( line ) == 0 ) then
      comment_num = comment_num + 1
      cycle
    end if

    row_num = row_num + 1

  end do

  close ( unit = input_unit )

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

  integer, parameter :: rk = kind ( 1.0D+00 )

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
subroutine r8col_normalize_li ( m, n, a )

!*****************************************************************************80
!
!! R8COL_NORMALIZE_LI normalizes an R8COL with the column infinity norm.
!
!  Discussion:
!
!    Each column is scaled so that the entry of maximum norm has the value 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input/output, real ( kind = rk ) A(M,N), the array to be normalized.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) c
  integer i
  integer j

  do j = 1, n

    c = a(1,j)

    do i = 2, m
      if ( abs ( c ) < abs ( a(i,j) ) ) then
        c = a(i,j)
      end if
    end do

    if ( c /= 0.0D+00 ) then
      a(1:m,j) = a(1:m,j) / c
    end if

  end do

  return
end
subroutine r8col_reverse ( m, n, a )

!*****************************************************************************80
!
!! R8COL_REVERSE reverses the order of columns in an R8COL.
!
!  Discussion:
!
!    To reverse the columns is to start with something like
!
!      11 12 13 14 15
!      21 22 23 24 25
!      31 32 33 34 35
!      41 42 43 44 45
!      51 52 53 54 55
!
!    and return
!
!      15 14 13 12 11
!      25 24 23 22 21
!      35 34 33 32 31
!      45 44 43 42 41
!      55 54 53 52 51
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input/output, real ( kind = rk ) A(M,N), the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  integer j
  integer jhi
  real ( kind = rk ) t(m)

  jhi = n / 2

  do j = 1, jhi
    t(1:m)       = a(1:m,j)
    a(1:m,j)     = a(1:m,n+1-j)
    a(1:m,n+1-j) = t(1:m)
  end do

  return
end
subroutine r8mat_data_read ( input_filename, m, n, table )

!*****************************************************************************80
!
!! R8MAT_DATA_READ reads data from an R8MAT file.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Discussion:
!
!    The file may contain more than N points, but this routine will
!    return after reading N of them.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of points.
!
!    Output, real ( kind = rk ) TABLE(M,N), the data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer ierror
  character ( len = * ) input_filename
  integer input_status
  integer input_unit
  integer j
  character ( len = 255 ) line
  real ( kind = rk ) table(m,n)
  real ( kind = rk ) x(m)

  ierror = 0

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = input_status )

  if ( input_status /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_DATA_READ - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' // &
      trim ( input_filename ) // '" on unit ', input_unit
    stop
  end if

  j = 0

  do while ( j < n )

    read ( input_unit, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8MAT_DATA_READ - Fatal error!'
      write ( *, '(a)' ) '  Error while reading lines of data.'
      write ( *, '(a,i8)' ) '  Number of values expected per line M = ', m
      write ( *, '(a,i8)' ) '  Number of data lines read, J =         ', j
      write ( *, '(a,i8)' ) '  Number of data lines needed, N =       ', n
      stop
    end if

    if ( line(1:1) == '#' .or. len_trim ( line ) == 0 ) then
      cycle
    end if

    call s_to_r8vec ( line, m, x, ierror )

    if ( ierror /= 0 ) then
      cycle
    end if

    j = j + 1

    table(1:m,j) = x(1:m)

  end do

  close ( unit = input_unit )

  return
end
subroutine r8mat_header_read ( input_filename, m, n )

!*****************************************************************************80
!
!! R8MAT_HEADER_READ reads the header from an R8MAT file.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Output, integer M, spatial dimension.
!
!    Output, integer N, the number of points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = * ) input_filename
  integer m
  integer n

  call file_column_count ( input_filename, m )

  if ( m <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data columns in'
    write ( *, '(a)' ) '  the file "' // trim ( input_filename ) // '".'
    stop
  end if

  call file_row_count ( input_filename, n )

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data rows in'
    write ( *, '(a)' ) '  the file "' // trim ( input_filename ) // '".'
    stop
  end if

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
subroutine r8mat_svd_linpack ( m, n, a, u, s, v )

!*****************************************************************************80
!
!! R8MAT_SVD_LINPACK gets the SVD of a matrix using a call to LINPACK.
!
!  Discussion:
!
!    The singular value decomposition of a real MxN matrix A has the form:
!
!      A = U * S * V'
!
!    where
!
!      U is MxM orthogonal,
!      S is MxN, and entirely zero except for the diagonal;
!      V is NxN orthogonal.
!
!    Moreover, the nonzero entries of S are positive, and appear
!    in order, from largest magnitude to smallest.
!
!    This routine calls the LINPACK routine DSVDC to compute the
!    factorization.
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
!    Input, integer M, N, the number of rows and columns
!    in the matrix A.
!
!    Input, real ( kind = rk ) A(M,N), the matrix whose singular value
!    decomposition we are investigating.
!
!    Output, real ( kind = rk ) U(M,M), S(M,N), V(N,N), the factors
!    that form the singular value decomposition of A.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) a_copy(m,n)
  real ( kind = rk ) e(max(m+1,n))
  integer i
  integer info
  integer lda
  integer ldu
  integer ldv
  integer job
  real ( kind = rk ) s(m,n)
  real ( kind = rk ) sdiag(max(m+1,n))
  real ( kind = rk ) u(m,m)
  real ( kind = rk ) v(n,n)
  real ( kind = rk ), allocatable, dimension ( : ) :: work

  allocate ( work(1:m) )
!
!  Compute the eigenvalues and eigenvectors.
!
  job = 11
  lda = m
  ldu = m
  ldv = n
!
!  The input matrix is destroyed by the routine.  Since we need to keep
!  it around, we only pass a copy to the routine.
!
  a_copy(1:m,1:n) = a(1:m,1:n)

  call dsvdc ( a_copy, lda, m, n, sdiag, e, u, ldu, v, ldv, work, job, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_SVD_LINPACK - Failure!'
    write ( *, '(a)' ) '  The SVD could not be calculated.'
    write ( *, '(a)' ) '  LINPACK routine DSVDC returned a nonzero'
    write ( *, '(a,i8)' ) '  value of the error flag, INFO = ', info
    return
  end if
!
!  Make the MxN matrix S from the diagonal values in SDIAG.
!
  s(1:m,1:n) = 0.0D+00
  do i = 1, min ( m, n )
    s(i,i) = sdiag(i)
  end do
!
!  Note that we do NOT need to transpose the V that comes out of LINPACK!
!
  deallocate ( work )

  return
end
subroutine r8mat_svd_low_rank ( m, n, r, u, s, v, a )

!*****************************************************************************80
!
!! R8MAT_SVD_LOW_RANK forms a rank R approximation using the SVD.
!
!  Discussion:
!
!    The singular value decomposition of a real MxN matrix A has the form:
!
!      A = U * S * V'
!
!    where 
!
!      U is MxM orthogonal,
!      S is MxN, and entirely zero except for the diagonal;
!      V is NxN orthogonal.
!
!    For R <= min ( M, N ), we can form a rank R approximation to A:
!
!      Ar = Ur * Sr * Vr'
!
!    where 
!
!      Ur is MxR orthogonal,
!      Sr is RxR, and entirely zero except for the diagonal;
!      Vr is NxR orthogonal.
!
!    Strictly speaking, we use RM = min ( R, M ), RN = min ( R, N ):
!
!      Ur is M  x RM,
!      Sr is RM x RN,
!      Vr is N  x RN
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns 
!    in the matrix A.
!
!    Input, integer R, the desired rank of the approximation.
!
!    Input, real ( kind = rk ) U(M,M), S(M,N), V(N,N), the factors
!    that form the singular value decomposition of A.
!
!    Output, real ( kind = rk ) A(M,N), the rank R SVD approximation to A.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  integer r
  integer rm
  integer rn
  real ( kind = rk ) s(m,n)
  real ( kind = rk ) u(m,m)
  real ( kind = rk ) v(n,n)

  rm = min ( r, m )
  rn = min ( r, n )

  a = matmul ( u(1:m,1:rm), matmul ( s(1:rm,1:rn), transpose ( v(1:n,1:rn) ) ) )

  return
end
subroutine r8row_reverse ( m, n, a )

!****************************************************************************80
!
!! R8ROW_REVERSE reverses the order of the rows of an R8ROW.
!
!  Discussion:
!
!    To reverse the rows is to start with something like
!
!      11 12 13 14 15
!      21 22 23 24 25
!      31 32 33 34 35
!      41 42 43 44 45
!      51 52 53 54 55
!
!    and return
!
!      51 52 53 54 55
!      41 42 43 44 45
!      31 32 33 34 35
!      21 22 23 24 25
!      11 12 13 14 15
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input/output, real ( kind = rk ) A(M,N), the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  integer i
  integer ihi
  real ( kind = rk ) t(n)

  ihi = m / 2

  do i = 1, ihi
    t(1:n) = a(i,1:n)
    a(i,1:n) = a(m+1-i,1:n)
    a(m+1-i,1:n) = t(1:n)
  end do

  return
end
subroutine r8vec_cum0 ( n, a, a_cum )

!*****************************************************************************80
!
!! R8VEC_CUM0 computes the cumulutive sums of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    Input:
!
!      A = (/ 1.0, 2.0, 3.0, 4.0 /)
!
!    Output:
!
!      A_CUM = (/ 0.0, 1.0, 3.0, 6.0, 10.0 /)
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real ( kind = rk ) A(N), the vector to be summed.
!
!    Output, real ( kind = rk ) A_CUM(0:N), the cumulative sums.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) a_cum(0:n)
  integer i

  a_cum(0) = 0.0D+00

  do i = 1, n
    a_cum(i) = a_cum(i-1) + a(i)
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
subroutine s_to_r8 ( s, dval, ierror, length )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 from a string.
!
!  Discussion:
!
!    The routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 blanks
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semicolon,
!
!    with most quantities optional.
!
!  Example:
!
!    S                 DVAL
!
!    '1'               1.0
!    '     1   '       1.0
!    '1A'              1.0
!    '12,34,56'        12.0
!    '  34 7'          34.0
!    '-1E2ABCD'        -100.0
!    '-1X2ABCD'        -1.0
!    ' 2E-1'           0.2
!    '23.45'           23.45
!    '-4.2E+2'         -420.0
!    '17d2'            1700.0
!    '-14e-2'         -0.14
!    'e2'              100.0
!    '-12.73e-9.23'   -12.73 * 10.0^(-9.23)
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string containing the
!    data to be read.  Reading will begin at position 1 and
!    terminate at the end of the string, or when no more
!    characters can be read to form a legal real.  Blanks,
!    commas, or other nonnumeric data will, in particular,
!    cause the conversion to halt.
!
!    Output, real ( kind = rk ) DVAL, the value read from the string.
!
!    Output, integer IERROR, error flag.
!    0, no errors occurred.
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer LENGTH, the number of characters read
!    to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character c
  logical ch_eqi
  real ( kind = rk ) dval
  integer ierror
  integer ihave
  integer isgn
  integer iterm
  integer jbot
  integer jsgn
  integer jtop
  integer length
  integer nchar
  integer ndig
  real ( kind = rk ) rbot
  real ( kind = rk ) rexp
  real ( kind = rk ) rtop
  character ( len = * ) s

  nchar = len_trim ( s )

  ierror = 0
  dval = 0.0D+00
  length = -1
  isgn = 1
  rtop = 0
  rbot = 1
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

  do

    length = length + 1

    if ( nchar < length+1 ) then
      exit
    end if

    c = s(length+1:length+1)
!
!  Blank character.
!
    if ( c == ' ' ) then

      if ( ihave == 2 ) then

      else if ( ihave == 6 .or. ihave == 7 ) then
        iterm = 1
      else if ( 1 < ihave ) then
        ihave = 11
      end if
!
!  Comma.
!
    else if ( c == ',' .or. c == ';' ) then

      if ( ihave /= 1 ) then
        iterm = 1
        ihave = 12
        length = length + 1
      end if
!
!  Minus sign.
!
    else if ( c == '-' ) then

      if ( ihave == 1 ) then
        ihave = 2
        isgn = -1
      else if ( ihave == 6 ) then
        ihave = 7
        jsgn = -1
      else
        iterm = 1
      end if
!
!  Plus sign.
!
    else if ( c == '+' ) then

      if ( ihave == 1 ) then
        ihave = 2
      else if ( ihave == 6 ) then
        ihave = 7
      else
        iterm = 1
      end if
!
!  Decimal point.
!
    else if ( c == '.' ) then

      if ( ihave < 4 ) then
        ihave = 4
      else if ( 6 <= ihave .and. ihave <= 8 ) then
        ihave = 9
      else
        iterm = 1
      end if
!
!  Scientific notation exponent marker.
!
    else if ( ch_eqi ( c, 'E' ) .or. ch_eqi ( c, 'D' ) ) then

      if ( ihave < 6 ) then
        ihave = 6
      else
        iterm = 1
      end if
!
!  Digit.
!
    else if (  ihave < 11 .and. lle ( '0', c ) .and. lle ( c, '9' ) ) then

      if ( ihave <= 2 ) then
        ihave = 3
      else if ( ihave == 4 ) then
        ihave = 5
      else if ( ihave == 6 .or. ihave == 7 ) then
        ihave = 8
      else if ( ihave == 9 ) then
        ihave = 10
      end if

      call ch_to_digit ( c, ndig )

      if ( ihave == 3 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = rk )
      else if ( ihave == 5 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = rk )
        rbot = 10.0D+00 * rbot
      else if ( ihave == 8 ) then
        jtop = 10 * jtop + ndig
      else if ( ihave == 10 ) then
        jtop = 10 * jtop + ndig
        jbot = 10 * jbot
      end if
!
!  Anything else is regarded as a terminator.
!
    else
      iterm = 1
    end if
!
!  If we haven't seen a terminator, and we haven't examined the
!  entire string, go get the next character.
!
    if ( iterm == 1 ) then
      exit
    end if

  end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LENGTH is equal to NCHAR.
!
  if ( iterm /= 1 .and. length + 1 == nchar ) then
    length = nchar
  end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
  if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then
    ierror = ihave
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'S_TO_R8 - Serious error!'
    write ( *, '(a)' ) '  Illegal or nonnumeric input:'
    write ( *, '(a)' ) '    ' // trim ( s )
    return
  end if
!
!  Number seems OK.  Form it.
!
  if ( jtop == 0 ) then
    rexp = 1.0D+00
  else
    if ( jbot == 1 ) then
      rexp = 10.0D+00 ** ( jsgn * jtop )
    else
      rexp = 10.0D+00 ** ( real ( jsgn * jtop, kind = rk ) &
        / real ( jbot, kind = rk ) )
    end if
  end if

  dval = real ( isgn, kind = rk ) * rexp * rtop / rbot

  return
end
subroutine s_to_r8vec ( s, n, rvec, ierror )

!*****************************************************************************80
!
!! S_TO_R8VEC reads an R8VEC from a string.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be read.
!
!    Input, integer N, the number of values expected.
!
!    Output, real ( kind = rk ) RVEC(N), the values read from the string.
!
!    Output, integer IERROR, error flag.
!    0, no errors occurred.
!    -K, could not read data for entries -K through N.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer ierror
  integer ilo
  integer lchar
  real ( kind = rk ) rvec(n)
  character ( len = * ) s

  i = 0
  ierror = 0
  ilo = 1

  do while ( i < n )

    i = i + 1

    call s_to_r8 ( s(ilo:), rvec(i), ierror, lchar )

    if ( ierror /= 0 ) then
      ierror = -i
      exit
    end if

    ilo = ilo + lchar

  end do

  return
end
subroutine s_word_count ( s, nword )

!*****************************************************************************80
!
!! S_WORD_COUNT counts the number of "words" in a string.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be examined.
!
!    Output, integer NWORD, the number of "words" in the string.
!    Words are presumed to be separated by one or more blanks.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  logical blank
  integer i
  integer lens
  integer nword
  character ( len = * ) s

  nword = 0
  lens = len ( s )

  if ( lens <= 0 ) then
    return
  end if

  blank = .true.

  do i = 1, lens

    if ( s(i:i) == ' ' ) then
      blank = .true.
    else if ( blank ) then
      nword = nword + 1
      blank = .false.
    end if

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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
