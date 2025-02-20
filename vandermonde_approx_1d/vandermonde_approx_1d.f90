subroutine vandermonde_approx_1d_coef ( n, m, x, y, c )

!*****************************************************************************80
!
!! VANDERMONDE_APPROX_1D_COEF computes a 1D polynomial approximant.
!
!  Discussion:
!
!    We assume the approximating function has the form
!
!      p(x) = c0 + c1 * x + c2 * x^2 + ... + cm * x^m.
!
!    We have n data values (x(i),y(i)) which must be approximated:
!
!      p(x(i)) = c0 + c1 * x(i) + c2 * x(i)^2 + ... + cm * x(i)^m = y(i)
!
!    This can be cast as an Nx(M+1) linear system for the polynomial
!    coefficients:
!
!      [ 1 x1 x1^2 ... x1^m ] [  c0 ] = [  y1 ]
!      [ 1 x2 x2^2 ... x2^m ] [  c1 ] = [  y2 ]
!      [ .................. ] [ ... ] = [ ... ]
!      [ 1 xn xn^2 ... xn^m ] [  cm ] = [  yn ]
!
!    In the typical case, N is greater than M+1 (we have more data and equations
!    than degrees of freedom) and so a least squares solution is appropriate,
!    in which case the computed polynomial will be a least squares approximant
!    to the data.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 September 2012
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
!    Input, real ( kind = rk ) X(N), Y(N), the data values.
!
!    Output, real ( kind = rk ) C(0:M), the coefficients of the approximating
!    polynomial.  C(0) is the constant term, and C(M) multiplies X^M.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(1:n,0:m)
  real ( kind = rk ) c(0:m)
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)

  call vandermonde_approx_1d_matrix ( n, m, x, a )

  call qr_solve ( n, m + 1, a, y, c )

  return
end
subroutine vandermonde_approx_1d_matrix ( n, m, x, a )

!*****************************************************************************80
!
!! VANDERMONDE_APPROX_1D_MATRIX computes a Vandermonde 1D approximation matrix.
!
!  Discussion:
!
!    We assume the approximant has the form
!
!      p(x) = c0 + c1 * x + c2 * x^2 + ... + cm * x^m.
!
!    We have n data values (x(i),y(i)) which must be approximated:
!
!      p(x(i)) = c0 + c1 * x(i) + c2 * x(i)^2 + ... + cm * x(i)^m = y(i)
!
!    This can be cast as an Nx(M+1) linear system for the polynomial
!    coefficients:
!
!      [ 1 x1 x1^2 ... x1^m ] [  c0 ] = [  y1 ]
!      [ 1 x2 x2^2 ... x2^m ] [  c1 ] = [  y2 ]
!      [ .................. ] [ ... ] = [ ... ]
!      [ 1 xn xn^2 ... xn^m ] [  cm ] = [  yn ]
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 September 2012
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
!    Input, real ( kind = rk ) X(N), the data values.
!
!    Output, real ( kind = rk ) A(N,0:M), the Vandermonde matrix for X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(n,0:m)
  integer j
  real ( kind = rk ) x(n)

  a(1:n,0) = 1.0D+00
  do j = 1, m
    a(1:n,j) = a(1:n,j-1) * x(1:n)
  end do

  return
end
