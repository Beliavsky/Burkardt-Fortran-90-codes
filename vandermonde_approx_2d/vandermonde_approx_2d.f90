function triangle_num ( n )

!*****************************************************************************80
!
!! TRIANGLE_NUM returns the N-th triangular number.
!
!  Definition:
!
!    The N-th triangular number T(N) is formed by the sum of the first
!    N integers:
!
!      T(N) = sum ( 1 <= I <= N ) I
!
!    By convention, T(0) = 0.
!
!  Formula:
!
!    T(N) = ( N * ( N + 1 ) ) / 2
!
!  First Values:
!
!     0
!     1
!     3
!     6
!    10
!    15
!    21
!    28
!    36
!    45
!    55
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the index of the desired number, 
!    which must be at least 0.
!
!    Output, integer TRIANGLE_NUM, the N-th triangular number.
!
  implicit none

  integer n
  integer triangle_num

  triangle_num = ( n * ( n + 1 ) ) / 2

  return
end
subroutine vandermonde_approx_2d_coef ( n, m, x, y, z, c )

!*****************************************************************************80
!
!! VANDERMONDE_APPROX_2D_COEF computes a 2D polynomial approximant.
!
!  Discussion:
!
!    We assume the approximating function has the form of a polynomial
!    in X and Y of total degree M.
!
!      p(x,y) = c00 
!             + c10 * x                + c01 *  y
!             + c20 * x^2   + c11 * xy + c02 * y^2
!             + ...
!             + cm0 * x^(m) + ...      + c0m * y^m.
!
!    If we let T(K) = the K-th triangular number 
!            = sum ( 1 <= I <= K ) I
!    then the number of coefficients in the above polynomial is T(M+1).
!
!    We have n data locations (x(i),y(i)) and values z(i) to approximate:
!
!      p(x(i),y(i)) = z(i)
!
!    This can be cast as an NxT(M+1) linear system for the polynomial
!    coefficients:
!
!      [ 1 x1 y1  x1^2 ... y1^m ] [ c00 ] = [  z1 ]
!      [ 1 x2 y2  x2^2 ... y2^m ] [ c10 ] = [  z2 ]
!      [ 1 x3 y3  x3^2 ... y3^m ] [ c01 ] = [  z3 ]
!      [ ...................... ] [ ... ] = [ ... ]
!      [ 1 xn yn  xn^2 ... yn^m ] [ c0m ] = [  zn ]
!
!    In the typical case, N is greater than T(M+1) (we have more data and 
!    equations than degrees of freedom) and so a least squares solution is 
!    appropriate, in which case the computed polynomial will be a least squares
!    approximant to the data.
!
!    The polynomial defined by the T(M+1) coefficients C could be evaluated 
!    at the Nx2-vector x by the command
!
!      pval = r8poly_value_2d ( m, c, n, x )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data points.
!
!    Input, integer M, the maximum degree of the polynomial.
!
!    Input, real ( kind = rk ) X(N), Y(N) the data locations.
!
!    Input, real ( kind = rk ) Z(N), the data values.
!
!    Output, real ( kind = rk ) C(T(M+1)), the coefficients of the approximating
!    polynomial.  C(1) is the constant term, and C(T(M+1)) multiplies Y^M.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ), allocatable :: a(:,:)
  real ( kind = rk ) c(*)
  integer m
  integer tm
  integer triangle_num
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)
  real ( kind = rk ) z(n)

  tm = triangle_num ( m + 1 )

  allocate ( a(1:n,1:tm) )
  call vandermonde_approx_2d_matrix ( n, m, tm, x, y, a )

  call qr_solve ( n, tm, a, z, c )

  deallocate ( a )

  return
end
subroutine vandermonde_approx_2d_matrix ( n, m, tm, x, y, a )

!*****************************************************************************80
!
!! VANDERMONDE_APPROX_2D_MATRIX computes a Vandermonde 2D approximation matrix.
!
!  Discussion:
!
!    We assume the approximating function has the form of a polynomial
!    in X and Y of total degree M.
!
!      p(x,y) = c00 
!             + c10 * x                + c01 * y
!             + c20 * x^2   + c11 * xy + c02 * y^2
!             + ...
!             + cm0 * x^(m) + ...      + c0m * y^m.
!
!    If we let T(K) = the K-th triangular number 
!            = sum ( 1 <= I <= K ) I
!    then the number of coefficients in the above polynomial is T(M+1).
!
!    We have n data locations (x(i),y(i)) and values z(i) to approximate:
!
!      p(x(i),y(i)) = z(i)
!
!    This can be cast as an NxT(M+1) linear system for the polynomial
!    coefficients:
!
!      [ 1 x1 y1  x1^2 ... y1^m ] [ c00 ] = [  z1 ]
!      [ 1 x2 y2  x2^2 ... y2^m ] [ c10 ] = [  z2 ]
!      [ 1 x3 y3  x3^2 ... y3^m ] [ c01 ] = [  z3 ]
!      [ ...................... ] [ ... ] = [ ... ]
!      [ 1 xn yn  xn^2 ... yn^m ] [ c0m ] = [  zn ]
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data points.
!
!    Input, integer M, the degree of the polynomial.
!
!    Input, integer TM, the M+1st triangular number.
!
!    Input, real ( kind = rk ) X(N), Y(N), the data locations.
!
!    Output, real ( kind = rk ) A(N,TM), the Vandermonde matrix for X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer tm

  real ( kind = rk ) a(n,tm)
  integer ex
  integer ey
  integer j
  integer m
  integer s
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)

  j = 0

  do s = 0, m
    do ex = s, 0, -1
      ey = s - ex
      j = j + 1
      a(1:n,j) = x(1:n) ** ex * y(1:n) ** ey 
    end do
  end do
 
  return
end
