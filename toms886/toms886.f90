subroutine cheb ( deg, pt, tcheb )

!*****************************************************************************80
!
!! cheb() computes normalized Chebyshev polynomials.
!
!  Discussion:
!
!    This subroutine computes the array TCHEB of normalized Chebyshev 
!    polynomials from degree 0 to DEG:
!      T_0(x)=1, 
!      T_j(x) = sqrt(2) * cos ( j * acos(x) ) 
!    at the point x = PT.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!  
!  Modified:
!
!    12 February 2014
!
!  Author:
!
!    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi, 
!    Marco Vianello.
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Marco Caliari, Stefano de Marchi, Marco Vianello,
!    Algorithm 886:
!    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
!    ACM Transactions on Mathematical Software,
!    Volume 35, Number 3, October 2008, Article 21, 11 pages.
!
!  Parameters:
!
!    Input, integer DEG, the degree.
!    0 <= DEG.
!
!    Input, real ( kind = rk ) PT, the evaluation point.
!
!    Output, real ( kind = rk ) TCHEB(0:DEG), the value of the normalized
!    Chebyshev polynomials of degrees 0 through DEG at the point PT.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer deg

  integer j
  real ( kind = rk ) pt
  real ( kind = rk ), parameter :: sqrt2 = 1.4142135623730951D+00
  real ( kind = rk ) tcheb(0:deg)

  if ( deg < 0 ) then
    return
  end if

  tcheb(0) = 1.0D+00

  if ( deg < 1 ) then
    return
  end if

  tcheb(1) = sqrt2 * pt
 
  if ( deg < 2 ) then
    return
  end if

  tcheb(2) = 2.0D+00 * pt * tcheb(1) - sqrt2 * tcheb(0)
!
!  Chebyshev recurrence.
!
  do j = 3, deg
    tcheb(j) = 2.0D+00 * pt * tcheb(j-1) - tcheb(j-2)
  end do

  return
end
subroutine dgemm ( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, &
  ldc )

!*****************************************************************************80
!
!! DGEMM computes C = alpha * A * B and related operations.
!
!  Discussion:
!
!    DGEMM performs one of the matrix-matrix operations
!
!     C := alpha * op ( A ) * op ( B ) + beta * C,
!
!    where op ( X ) is one of
!
!      op ( X ) = X   or   op ( X ) = X',
!
!    ALPHA and BETA are scalars, and A, B and C are matrices, with op ( A )
!    an M by K matrix, op ( B ) a K by N matrix and C an N by N matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!    
!  Modified:
!
!    12 February 2014
!
!  Author:
!
!    Original FORTRAN77 version by Jack Dongarra.
!    This FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, character * 1 TRANSA, specifies the form of op( A ) to be used in
!    the matrix multiplication as follows:
!    'N' or 'n', op ( A ) = A.
!    'T' or 't', op ( A ) = A'.
!    'C' or 'c', op ( A ) = A'.
!
!    Input, character * 1 TRANSB, specifies the form of op ( B ) to be used in
!    the matrix multiplication as follows:
!    'N' or 'n', op ( B ) = B.
!    'T' or 't', op ( B ) = B'.
!    'C' or 'c', op ( B ) = B'.
!
!    Input, integer M, the number of rows of the  matrix op ( A ) 
!    and of the matrix C.  0 <= M.
!
!    Input, integer N, the number  of columns of the matrix 
!    op ( B ) and the number of columns of the matrix C.  0 <= N.
!
!    Input, integer K, the number of columns of the matrix op ( A )
!    and the number of rows of the matrix op ( B ).  0 <= K.
!
!    Input, real ( kind = rk ) ALPHA, the scalar multiplier 
!    for op ( A ) * op ( B ).
!
!    Input, real ( kind = rk ) A(LDA,KA), where:
!    if TRANSA is 'N' or 'n', KA is equal to K, and the leading M by K
!    part of the array contains A;
!    if TRANSA is not 'N' or 'n', then KA is equal to M, and the leading
!    K by M part of the array must contain the matrix A.
!
!    Input, integer LDA, the first dimension of A as declared in 
!    the calling routine.  When TRANSA = 'N' or 'n' then LDA must be at least 
!    max ( 1, M ), otherwise LDA must be at least max ( 1, K ).
!
!    Input, real ( kind = rk ) B(LDB,KB), where:
!    if TRANSB is 'N' or 'n', kB is N, and the leading K by N 
!    part of the array contains B;
!    if TRANSB is not 'N' or 'n', then KB is equal to K, and the leading
!    n by k  part of the array must contain the matrix B.
!
!    Input, integer LDB, the first dimension of B as declared in 
!    the calling routine.  When TRANSB = 'N' or 'n' then LDB must be at least 
!    max ( 1, K ), otherwise LDB must be at least max ( 1, N ).
!
!    Input, real ( kind = rk ) BETA, the scalar multiplier for C.
!
!    Input, real ( kind = rk ) C(LDC,N).
!    Before entry, the leading M by N part of this array must contain the 
!    matrix C, except when BETA is zero, in which case C need not be set 
!    on entry.
!    On exit, the array C is overwritten by the M by N matrix
!      alpha * op ( A ) * op ( B ) + beta * C.
!
!    Input, integer LDC, the first dimension of C as declared in 
!    the calling routine.  max ( 1, M ) <= LDC.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer lda
  integer ldb
  integer ldc

  real ( kind = rk ) a(lda,*)
  real ( kind = rk ) alpha
  real ( kind = rk ) b(ldb,*)
  real ( kind = rk ) beta
  real ( kind = rk ) c(ldc,*)
  integer i
  integer info
  integer j
  integer k
  integer l
  integer m
  integer n
  integer ncola
  integer nrowa
  integer nrowb
  logical nota
  logical notb
  real ( kind = rk ) temp
  character * 1 transa
  character * 1 transb
!
!  Set NOTA and NOTB as true if A and B respectively are not
!  transposed and set NROWA, NCOLA and NROWB as the number of rows
!  and columns of A and the number of rows of B respectively.
!
  nota = ( transa == 'N' .or. transa == 'n' )

  if ( nota ) then
    nrowa = m
    ncola = k
  else
    nrowa = k
    ncola = m
  end if

  notb = ( transb == 'N' .or. transb == 'n' )

  if ( notb ) then
    nrowb = k
  else
    nrowb = n
  end if
!
!  Test the input parameters.
!
  info = 0

  if ( transa /= 'N' .and. transa /= 'n' .and. &
       transa /= 'C' .and. transa /= 'c' .and. &
       transa /= 'T' .and. transa /= 't' ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'DGEMM - Fatal error!'
    write ( *, '(a)' ) '  Input TRANSA had illegal value.'
    stop 1
  end if

  if ( transb /= 'N' .and. transb /= 'n' .and. &
       transb /= 'C' .and. transb /= 'c' .and. &
       transb /= 'T' .and. transb /= 't' ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'DGEMM - Fatal error!'
    write ( *, '(a)' ) '  Input TRANSB had illegal value.'
    stop 1
  end if

  if ( m < 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'DGEMM - Fatal error!'
    write ( *, '(a)' ) '  Input M had illegal value.'
    stop 1
  end if

  if ( n < 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'DGEMM - Fatal error!'
    write ( *, '(a)' ) '  Input N had illegal value.'
    stop 1
  end if

  if ( k < 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'DGEMM - Fatal error!'
    write ( *, '(a)' ) '  Input K had illegal value.'
    stop 1
  end if

  if ( lda < max ( 1, nrowa ) ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'DGEMM - Fatal error!'
    write ( *, '(a)' ) '  Input LDA had illegal value.'
    stop 1
  end if

  if ( ldb < max ( 1, nrowb ) ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'DGEMM - Fatal error!'
    write ( *, '(a)' ) '  Input LDB had illegal value.'
    stop 1
  end if

  if ( ldc < max ( 1, m ) ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'DGEMM - Fatal error!'
    write ( *, '(a)' ) '  Input LDC had illegal value.'
    stop 1
  end if
!
!  Quick return if possible.
!
  if ( m == 0 ) then
    return
  end if

  if ( n == 0 ) then
    return
  end if

  if ( ( alpha == 0.0D+00 .or. k == 0 ) .and. beta == 1.0D+00 ) then
    return
  end if
!
!  And if alpha is zero.
!
  if ( alpha == 0.0D+00 ) then
    if ( beta == 0.0D+00 ) then
      c(1:m,1:n) = 0.0D+00
    else
      c(1:m,1:n) = beta * c(1:m,1:n)
    end if
    return
  end if
!
!  Start the operations.
!
  if ( notb ) then
!
!  Form  C := alpha*A*B + beta*C.
!
    if ( nota ) then

      do j = 1, n

        if ( beta == 0.0D+00 ) then
          c(1:m,j) = 0.0D+00
        else if ( beta /= 1.0D+00 ) then
          c(1:m,j) = beta * c(1:m,j)
        end if

        do l = 1, k
          if ( b(l,j) /= 0.0D+00 ) then
            c(1:m,j) = c(1:m,j) + alpha * b(l,j) * a(1:m,l)
          end if
        end do

      end do
!
!  Form  C := alpha*A'*B + beta*C
!
    else

      do j = 1, n
        do i = 1, m

          temp = dot_product ( a(1:k,i), b(1:k,j) )

          if ( beta == 0.0D+00 ) then
            c(i,j) = alpha * temp
          else
            c(i,j) = alpha * temp + beta * c(i,j)
          end if

        end do
      end do

    end if
!
!  Form  C := alpha*A*B' + beta*C
!
  else

    if ( nota ) then

      do j = 1, n

        if ( beta == 0.0D+00 ) then
          c(1:m,j) = 0.0D+00
        else if ( beta /= 1.0D+00 ) then
          c(1:m,j) = beta * c(1:m,j)
        end if

        do l = 1, k
          if ( b(j,l) /= 0.0D+00 ) then
            c(1:m,j) = c(1:m,j) + alpha * b(j,l) * a(1:m,l)
          end if
        end do

      end do
!
!  Form  C := alpha*A'*B' + beta*C
!
    else

      do j = 1, n
        do i = 1, m

          temp = dot_product ( a(1:k,i), b(j,1:k) )
          if ( beta == 0.0D+00 ) then
            c(i,j) = alpha * temp
          else
            c(i,j) = alpha * temp + beta * c(i,j)
          end if
        end do
      end do

    end if

  end if

  return
end
function franke ( x, y )

!*****************************************************************************80
!
!! FRANKE returns the value of the Franke function #1.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Richard Franke,
!    Scattered Data Interpolation: Tests of Some Methods,
!    Mathematics of Computation,
!    Volume 38, Number 157, January 1982, pages 181-200.
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, Y, the evalution points.
!
!    Output, real ( kind = rk ) FRANKE, the function values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) franke
  real ( kind = rk ) x
  real ( kind = rk ) y

  franke = &
      0.75D+00 * exp ( &
        - ( ( 9.0D+00 * x - 2.0D+00 )**2 &
          + ( 9.0D+00 * y - 2.0D+00 )**2 ) / 4.0D+00 ) &
    + 0.75D+00 * exp ( &
        - ( ( 9.0D+00 * x + 1.0D+00 )**2 ) / 49.0D+00 &
          - ( 9.0D+00 * y + 1.0D+00 )      / 10.0D+00 ) &      
    + 0.5D+00  * exp ( &
        - ( ( 9.0D+00 * x - 7.0D+00 )**2 &
          + ( 9.0D+00 * y - 3.0D+00 )**2 ) / 4.0D+00 ) &
    - 0.2D+00  * exp ( &
          - ( 9.0D+00 * x - 4.0D+00 )**2 &
          - ( 9.0D+00 * y - 7.0D+00 )**2 )

  return
end
subroutine padua2 ( deg, degmax, npd, wpd, fpd, raux1, raux2, c0, esterr )

!*****************************************************************************80
!
!! PADUA2 computes the Padua interpolation coefficient matrix.
!
!  Discussion:
!
!    This subroutine computes the coefficient matrix C0, in the 
!    orthonormal Chebyshev basis T_j(x)T_{k-j}(y), 0 <= j <= k <= DEG, 
!    T_0(x)=1, T_j(x) = sqrt(2) * cos(j * acos(x)), of the 
!    interpolation polynomial of degree DEG of the function values FPD 
!    at the set of NPD Padua points (PD1,PD2) in the square [-1,1]^2. 
!
!    The interpolant may be evaluated at an arbitrary point by the 
!    function PD2VAL. PD1, PD2 and WPD are the Padua points and weights 
!    computed by PDPTS.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!  
!  Modified:
!
!    12 February 2014
!
!  Author:
!
!    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi, 
!    Marco Vianello.
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Marco Caliari, Stefano de Marchi, Marco Vianello,
!    Algorithm 886:
!    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
!    ACM Transactions on Mathematical Software,
!    Volume 35, Number 3, October 2008, Article 21, 11 pages.
!
!  Parameters:
!
!    Input, integer DEG, the degree of approximation.
!
!    Input, integer DEGMAX, the maximum degree allowed.
!
!    Input, integer NPD, the number of Padua points.
!
!    Input, real ( kind = rk ) WPD(NPD), the weights.
!
!    Input, real ( kind = rk ) FPD(NPD), the value at the Padua points
!    of the function to be interpolated.
!
!    Workspace, real ( kind = rk ) RAUX1(0:DEGMAX,DEG+2).
!
!    Workspace, real ( kind = rk ) RAUX2(0:DEGMAX,DEG+2).
!
!    Output, real ( kind = rk ) C0(0:DEG,0:DEG), the coefficient matrix.
!
!    Output, real ( kind = rk ) ESTERR, the estimated error.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer deg
  integer degmax
  integer npd

  real ( kind = rk ) angle
  real ( kind = rk ) c0(0:degmax+1,0:deg+1)
  real ( kind = rk ) esterr
  real ( kind = rk ) fpd(npd)
  integer i
  integer j
  integer k
  real ( kind = rk ), parameter :: pi = 3.1415926535897931D+00
  real ( kind = rk ) pt
  real ( kind = rk ) raux1(0:degmax,deg+2)
  real ( kind = rk ) raux2(0:degmax,deg+2)
  real ( kind = rk ) wpd(npd)
!
!  Build the matrix P_2 and store it in RAUX2.
!
  do i = 0, deg + 1
    angle = real ( i, kind = rk ) * pi / real ( deg + 1, kind = rk )
    pt = - cos ( angle )
    call cheb ( deg, pt, raux2(0:deg,i+1) )
  end do
!
!  Build the matrix G(f) and store it in C0.
!
  c0(1:degmax+1,0:deg+1) = 0.0D+00

  k = 0
  do j = 0, deg + 1
    do i = 0, deg
      if ( mod ( i + j, 2 ) == 0 ) then
        k = k + 1
        c0(i,j) = fpd(k) * wpd(k)
      else
        c0(i,j) = 0.0D+00
      end if
    end do
  end do
!
!  Compute the matrix-matrix product G(f)*P_2' and store it in RAUX1.
!
  call dgemm ( 'n', 't', deg + 1, deg + 1, deg + 2, 1.0D+00, &
    c0, degmax + 2, raux2, degmax + 1, 0.0D+00, raux1, degmax + 1 )
!
!  Build the matrix P_1 and store it in RAUX2.
!
  do i = 0, deg
    angle = real ( i, kind = rk ) * pi / real ( deg, kind = rk )
    pt = - cos ( angle )
    call cheb ( deg, pt, raux2(0:deg,i+1) )
  end do
!
!  Compute the matrix-matrix product C(f) = P_1 * ( G(f) * P_2' ) 
!  and store it in C0.
!
  call dgemm ( 'n', 'n', deg + 1, deg + 1, deg + 1, 1.0D+00, &
    raux2, degmax + 1, raux1, degmax + 1, 0.0D+00, c0, degmax + 2 )

  c0(deg,0) = c0(deg,0) / 2.0D+00
!
!  Estimate the error.
!
  esterr = 0.0D+00
  do j = 0, 2
    do i = 0, deg - j
      esterr = esterr + abs ( c0(i,deg-i-j) )
    end do
  end do
  esterr = 2.0D+00 * esterr

  return
end
function pd2val ( deg, degmax, c0, tg1, tg2 )

!*****************************************************************************80
!
!! PD2VAL evaluates the Padua2 interpolant.
!
!  Discussion:
!
!    This function returns the value of the interpolant at (TG1,TG2).
!    C0 is the matrix of the coefficients computed by PADUA2.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!  
!  Modified:
!
!    13 February 2014
!
!  Author:
!
!    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi, 
!    Marco Vianello.
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Marco Caliari, Stefano de Marchi, Marco Vianello,
!    Algorithm 886:
!    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
!    ACM Transactions on Mathematical Software,
!    Volume 35, Number 3, October 2008, Article 21, 11 pages.
!
!  Parameters:
!
!    Input, integer DEG, the degree of approximation.
!
!    Input, integer DEGMAX, the maximum degree allowed.         
!
!    Input, real ( kind = rk ) C0(0:DEGMAX+1,0:DEG), the coefficient matrix.
!
!    Input, real ( kind = rk ) TG1, TG2, the first and second coordinates of
!    the target point.
!
!    Output, real ( kind = rk ) PD2VAL, the value of the interpolant at
!    the target point.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer deg
  integer degmax

  real ( kind = rk ) c0(0:degmax+1,0:deg)
  integer i
  real ( kind = rk ) pd2val
  real ( kind = rk ) tg1
  real ( kind = rk ) tg2
  real ( kind = rk ) ttg1(0:deg)
  real ( kind = rk ) ttg2(0:deg)
!
!  Compute the normalized Chebyshev polynomials at the target point.
!     
  call cheb ( deg, tg1, ttg1 )
  call cheb ( deg, tg2, ttg2 )
! 
!  Evaluate the interpolant
!
  pd2val = 0.0D+00
  do i = deg, 0, -1
    pd2val = pd2val + ttg2(deg-i) * dot_product ( ttg1(0:i), c0(0:i,deg-i) )
  end do

  return
end
subroutine pdpts ( deg, pd1, pd2, wpd, npd )

!*****************************************************************************80
!
!! PDPTS returns the points and weights for Padua interpolation.
!
!  Discussion:
!
!    This subroutine computes the first family of Padua points and 
!    weights corresponding to degree DEG.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!  
!  Modified:
!
!    12 February 2014
!
!  Author:
!
!    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi, 
!    Marco Vianello.
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Marco Caliari, Stefano de Marchi, Marco Vianello,
!    Algorithm 886:
!    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
!    ACM Transactions on Mathematical Software,
!    Volume 35, Number 3, October 2008, Article 21, 11 pages.
!
!  Parameters:
!
!    Input, integer DEG, the degree of approximation.
!
!    Output, real ( kind = rk ) PD1(NPD), PD2(NPD), the first and second
!    coordinates of the Padua points
!
!    Output, real ( kind = rk ) WPD(NPD), the weights.
!
!    Output, integer NPD, the number of Padua points.
!    NPD = ( DEG + 1 ) * ( DEG + 2 ) / 2.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer deg
  integer itemp0
  integer j
  integer k
  integer npd
  real ( kind = rk ) pd1(*)
  real ( kind = rk ) pd2(*)
  real ( kind = rk ), parameter :: pi = 3.1415926535897931D+00
  real ( kind = rk ) rtemp0
  real ( kind = rk ) wpd(*)
!
!  Compute the Padua points of the first family at degree DEG.
!
  if ( deg == 0 ) then    
    pd1(1) = -1.0D+00
    pd2(1) = -1.0D+00
    wpd(1) = 2.0D+00
    npd = 1
    return
  end if
   
  npd = 0
  itemp0 = deg * ( deg + 1 )
  rtemp0 = pi / itemp0

  do j = 0, deg + 1
    do k = mod ( j, 2 ), deg, 2

      npd = npd + 1
      pd1(npd) = - cos ( ( deg + 1 ) * k * rtemp0 )
      pd2(npd) = - cos ( deg * j * rtemp0 )
      wpd(npd) = 2.0D+00 / itemp0

      if ( k == 0 .or. k == deg ) then
        wpd(npd) = wpd(npd) / 2.0D+00
      end if

      if ( j == 0 .or. j == deg + 1 ) then
        wpd(npd) = wpd(npd) / 2.0D+00
      end if

    end do
  end do

  return
end
function r8_huge ( )

!*****************************************************************************80
!
!! R8_HUGE returns a very large R8.
!
!  Discussion:
!
!    The value returned by this function is NOT required to be the
!    maximum representable R8.  This value varies from machine to machine,
!    from compiler to compiler, and may cause problems when being printed.
!    We simply want a "very large" but non-infinite number.
!
!    FORTRAN90 provides a built-in routine HUGE ( X ) that
!    can return the maximum representable number of the same datatype
!    as X, if that is what is really desired.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 October 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) R8_HUGE, a "huge" value.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r8_huge

  r8_huge = 1.0D+30

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
