subroutine imtqlx ( n, d, e, z )

!*****************************************************************************80
!
!! imtqlx() diagonalizes a symmetric tridiagonal matrix.
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
!    This version by John Burkardt.
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
        stop
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
function j_double_product_integral ( i, j, a, b )

!*****************************************************************************80
!
!! j_double_product_integral(): integral of J(i,x)*J(j,x)*(1-x)^a*(1+x)^b.
!
!  Discussion:
!
!    VALUE = integral ( -1 <= x <= +1 ) J(i,x)*J(j,x)*(1-x)^a*(1+x)^b dx
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, J, the polynomial indices.
!
!    Input, real ( kind = rk ) A, B, the parameters.
!    -1 < A, B.
!
!    Output, real ( kind = rk ) VALUE, the value of the integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  integer i
  real ( kind = rk ) i_r8
  integer j
  real ( kind = rk ) j_double_product_integral
  real ( kind = rk ) r8_factorial
  real ( kind = rk ) value

  if ( i /= j ) then

    value = 0.0D+00

  else

    i_r8 = real ( i, kind = rk )

    value = 2.0D+00 ** ( a + b + 1.0D+00 ) &
      / ( 2.0D+00 * i_r8 + a + b + 1.0D+00 ) &
      * gamma ( i_r8 + a + 1.0D+00 ) &
      * gamma ( i_r8 + b + 1.0D+00 ) &
      / r8_factorial ( i ) &
      / gamma ( i_r8 + a + b + 1.0D+00 )

  end if

  j_double_product_integral = value

  return
end
function j_integral ( n )

!*****************************************************************************80
!
!! j_integral() evaluates a monomial integral associated with J(n,a,b,x).
!
!  Discussion:
!
!    The integral:
!
!      integral ( -1 <= x < +1 ) x^n dx
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the exponent.
!    0 <= N.
!
!    Output, real ( kind = rk ) J_INTEGRAL, the value of the integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) j_integral
  integer n
  real ( kind = rk ) value

  if ( mod ( n, 2 ) == 1 ) then
    value = 0.0D+00
  else
    value = 2.0D+00 / real ( n + 1, kind = rk )
  end if

  j_integral = value

  return
end
subroutine j_polynomial ( m, n, alpha, beta, x, cx )

!*****************************************************************************80
!
!! j_polynomial() evaluates the Jacobi polynomials J(n,a,b,x).
!
!  Differential equation:
!
!    (1-X*X) Y'' + (BETA-ALPHA-(ALPHA+BETA+2) X) Y' + N (N+ALPHA+BETA+1) Y = 0
!
!  Recursion:
!
!    P(0,ALPHA,BETA,X) = 1,
!
!    P(1,ALPHA,BETA,X) = ( (2+ALPHA+BETA)*X + (ALPHA-BETA) ) / 2
!
!    P(N,ALPHA,BETA,X)  = 
!      ( 
!        (2*N+ALPHA+BETA-1) 
!        * ((ALPHA^2-BETA^2)+(2*N+ALPHA+BETA)*(2*N+ALPHA+BETA-2)*X) 
!        * P(N-1,ALPHA,BETA,X)
!        -2*(N-1+ALPHA)*(N-1+BETA)*(2*N+ALPHA+BETA) * P(N-2,ALPHA,BETA,X)
!      ) / 2*N*(N+ALPHA+BETA)*(2*N-2+ALPHA+BETA)
!
!  Restrictions:
!
!    -1 < ALPHA
!    -1 < BETA
!
!  Norm:
!
!    Integral ( -1 <= X <= 1 ) ( 1 - X )^ALPHA * ( 1 + X )^BETA 
!      * P(N,ALPHA,BETA,X)^2 dX 
!    = 2^(ALPHA+BETA+1) * Gamma ( N + ALPHA + 1 ) * Gamma ( N + BETA + 1 ) /
!      ( 2 * N + ALPHA + BETA ) * N! * Gamma ( N + ALPHA + BETA + 1 )
!
!  Special values:
!
!    P(N,ALPHA,BETA,1) = (N+ALPHA)!/(N!*ALPHA!) for integer ALPHA.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 April 2012
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
!    Input, integer N, the highest order polynomial to compute.  
!    Note that polynomials 0 through N will be computed.
!
!    Input, real ( kind = rk ) ALPHA, one of the parameters defining the Jacobi
!    polynomials, ALPHA must be greater than -1.
!
!    Input, real ( kind = rk ) BETA, the second parameter defining the Jacobi
!    polynomials, BETA must be greater than -1.
!
!    Input, real ( kind = rk ) X(M), the point at which the polynomials are 
!    to be evaluated.
!
!    Output, real ( kind = rk ) CX(M,0:N), the values of the first N+1 Jacobi
!    polynomials at the point X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) alpha
  real ( kind = rk ) beta
  real ( kind = rk ) cx(1:m,0:n)
  real ( kind = rk ) c1
  real ( kind = rk ) c2
  real ( kind = rk ) c3
  real ( kind = rk ) c4
  integer i
  real ( kind = rk ) r_i
  real ( kind = rk ) x(m)

  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'J_POLYNOMIAL - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Illegal input value of ALPHA = ', alpha
    write ( *, '(a)' ) '  But ALPHA must be greater than -1.'
    stop
  end if
 
  if ( beta <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'J_POLYNOMIAL - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Illegal input value of BETA = ', beta
    write ( *, '(a)' ) '  But BETA must be greater than -1.'
    stop
  end if
  
  if ( n < 0 ) then
    return
  end if

  cx(1:m,0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  cx(1:m,1) = ( 1.0D+00 + 0.5D+00 * ( alpha + beta ) ) * x(1:m) &
    + 0.5D+00 * ( alpha - beta )
 
  do i = 2, n

    r_i = real ( i, kind = rk ) 

    c1 = 2.0D+00 * r_i * ( r_i + alpha + beta ) &
       * ( 2.0D+00 * r_i - 2.0D+00 + alpha + beta )

    c2 = ( 2.0D+00 * r_i - 1.0D+00 + alpha + beta ) &
       * ( 2.0D+00 * r_i  + alpha + beta ) &
       * ( 2.0D+00 * r_i - 2.0D+00 + alpha + beta )

    c3 = ( 2.0D+00 * r_i - 1.0D+00 + alpha + beta ) &
      * ( alpha + beta ) * ( alpha - beta )

    c4 = - 2.0D+00 * ( r_i - 1.0D+00 + alpha ) &
      * ( r_i - 1.0D+00 + beta )  * ( 2.0D+00 * r_i + alpha + beta )

    cx(1:m,i) = ( ( c3 + c2 * x(1:m) ) * cx(1:m,i-1) + c4 * cx(1:m,i-2) ) / c1

  end do

  return
end
subroutine j_polynomial_values ( n_data, n, a, b, x, fx )

!*****************************************************************************80
!
!! j_polynomial_values() returns some values of the Jacobi polynomial.
!
!  Discussion:
!
!    In Mathematica, the function
!
!      JacobiP[ n, a, b, x ]
!
!    returns the value of the Jacobi polynomial.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 July 2004
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
!    Output, integer N, the degree of the polynomial.
!
!    Output, real ( kind = rk ) A, B, parameters of the function.
!
!    Output, real ( kind = rk ) X, the argument of the function.
!
!    Output, real ( kind = rk ) FX, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 26

  real ( kind = rk ) a
  real ( kind = rk ), save, dimension ( n_max ) :: a_vec = (/ &
     0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, &
     0.0D+00, 0.0D+00, 1.0D+00, 2.0D+00, &
     3.0D+00, 4.0D+00, 5.0D+00, 0.0D+00, &
     0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, &
     0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, &
     0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, &
     0.0D+00, 0.0D+00 /)
  real ( kind = rk ) b
  real ( kind = rk ), save, dimension ( n_max ) :: b_vec = (/ &
    1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00, 2.0D+00, &
    3.0D+00, 4.0D+00, 5.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00 /)
  real ( kind = rk ) fx
  real ( kind = rk ), save, dimension ( n_max ) :: fx_vec = (/ &
      1.000000000000000D+00, &
      0.250000000000000D+00, &
     -0.375000000000000D+00, &
     -0.484375000000000D+00, &
     -0.132812500000000D+00, &
      0.275390625000000D+00, &
     -0.164062500000000D+00, &
     -1.174804687500000D+00, &
     -2.361328125000000D+00, &
     -2.616210937500000D+00, &
      0.117187500000000D+00, &
      0.421875000000000D+00, &
      0.504882812500000D+00, &
      0.509765625000000D+00, &
      0.430664062500000D+00, &
     -6.000000000000000D+00, &
      0.038620000000000D+00, &
      0.811840000000000D+00, &
      0.036660000000000D+00, &
     -0.485120000000000D+00, &
     -0.312500000000000D+00, &
      0.189120000000000D+00, &
      0.402340000000000D+00, &
      0.012160000000000D+00, &
     -0.439620000000000D+00, &
      1.000000000000000D+00 /)
  integer n
  integer n_data
  integer, save, dimension ( n_max ) :: n_vec = (/ &
     0, 1, 2, 3, &
     4, 5, 5, 5, &
     5, 5, 5, 5, &
     5, 5, 5, 5, &
     5, 5, 5, 5, &
     5, 5, 5, 5, &
     5, 5 /)
  real ( kind = rk ) x
  real ( kind = rk ), save, dimension ( n_max ) :: x_vec = (/ &
     0.5D+00,  0.5D+00,  0.5D+00,  0.5D+00, &
     0.5D+00,  0.5D+00,  0.5D+00,  0.5D+00, &
     0.5D+00,  0.5D+00,  0.5D+00,  0.5D+00, &
     0.5D+00,  0.5D+00,  0.5D+00, -1.0D+00, &
    -0.8D+00, -0.6D+00, -0.4D+00, -0.2D+00, &
     0.0D+00,  0.2D+00,  0.4D+00,  0.6D+00, &
     0.8D+00,  1.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    a = 0.0D+00
    b = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    a = a_vec(n_data)
    b = b_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine j_polynomial_zeros ( n, alpha, beta, x )

!*****************************************************************************80
!
!! j_polynomial_zeros(): zeros of Jacobi polynomial J(n,a,b,x).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 April 2012
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer, N, the order.
!
!    Input, real ( kind = rk ), ALPHA, BETA, the parameters.
!    -1 < ALPHA, BETA.
!
!    Output, real ( kind = rk ) X(N), the zeros.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a2b2
  real ( kind = rk ) ab
  real ( kind = rk ) abi
  real ( kind = rk ) alpha
  real ( kind = rk ) beta
  real ( kind = rk ) bj(n)
  integer i
  real ( kind = rk ) i_r8
  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n)
  real ( kind = rk ) zemu

  ab = alpha + beta
  abi = 2.0D+00 + ab
!
!  Define the zero-th moment.
!
  zemu = 2.0D+00 ** ( ab + 1.0D+00 ) * gamma ( alpha + 1.0D+00 ) &
    * gamma ( beta + 1.0D+00 ) / gamma ( abi )
!
!  Define the Jacobi matrix.
!
  x(1) = ( beta - alpha ) / abi
  x(2:n) = 0.0D+00

  bj(1) = 4.0D+00 * ( 1.0D+00 + alpha ) * ( 1.0D+00 + beta ) &
    / ( ( abi + 1.0D+00 ) * abi * abi )
  bj(2:n) = 0.0D+00

  a2b2 = beta * beta - alpha * alpha

  do i = 2, n
    i_r8 = real ( i, kind = rk )
    abi = 2.0D+00 * i_r8 + ab
    x(i) = a2b2 / ( ( abi - 2.0D+00 ) * abi )
    abi = abi ** 2
    bj(i) = 4.0D+00 * i_r8 * ( i_r8 + alpha ) * ( i_r8 + beta ) &
      * ( i_r8 + ab ) / ( ( abi - 1.0D+00 ) * abi )
  end do
  bj(1:n) =  sqrt ( bj(1:n) )

  w(1) = sqrt ( zemu )
  w(2:n) = 0.0D+00
!
!  Diagonalize the Jacobi matrix.
!
  call imtqlx ( n, x, bj, w )

  return
end
subroutine j_quadrature_rule ( n, alpha, beta, x, w )

!*****************************************************************************80
!
!! j_quadrature_rule(): Gauss-Jacobi quadrature based on J(n,a,b,x).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 April 2012
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer, N, the order.
!
!    Input, real ( kind = rk ), ALPHA, BETA, the parameters.
!    -1 < ALPHA, BETA.
!
!    Output, real ( kind = rk ) X(N), the abscissas.
!
!    Output, real ( kind = rk ) W(N), the weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a2b2
  real ( kind = rk ) ab
  real ( kind = rk ) abi
  real ( kind = rk ) alpha
  real ( kind = rk ) beta
  real ( kind = rk ) bj(n)
  integer i
  real ( kind = rk ) i_r8
  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n)
  real ( kind = rk ) zemu

  ab = alpha + beta
  abi = 2.0D+00 + ab
!
!  Define the zero-th moment.
!
  zemu = 2.0D+00 ** ( ab + 1.0D+00 ) * gamma ( alpha + 1.0D+00 ) &
    * gamma ( beta + 1.0D+00 ) / gamma ( abi )
!
!  Define the Jacobi matrix.
!
  x(1) = ( beta - alpha ) / abi
  x(2:n) = 0.0D+00

  bj(1) = 4.0D+00 * ( 1.0D+00 + alpha ) * ( 1.0D+00 + beta ) &
    / ( ( abi + 1.0D+00 ) * abi * abi )
  bj(2:n) = 0.0D+00

  a2b2 = beta * beta - alpha * alpha

  do i = 2, n
    i_r8 = real ( i, kind = rk )
    abi = 2.0D+00 * i_r8 + ab
    x(i) = a2b2 / ( ( abi - 2.0D+00 ) * abi )
    abi = abi ** 2
    bj(i) = 4.0D+00 * i_r8 * ( i_r8 + alpha ) * ( i_r8 + beta ) &
      * ( i_r8 + ab ) / ( ( abi - 1.0D+00 ) * abi )
  end do
  bj(1:n) =  sqrt ( bj(1:n) )

  w(1) = sqrt ( zemu )
  w(2:n) = 0.0D+00
!
!  Diagonalize the Jacobi matrix.
!
  call imtqlx ( n, x, bj, w )

  w(1:n) = w(1:n) ** 2

  return
end
function r8_factorial ( n )

!*****************************************************************************80
!
!! r8_factorial() computes the factorial of N.
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

  r8_factorial = 1.0D+00

  do i = 1, n
    r8_factorial = r8_factorial * real ( i, kind = rk )
  end do

  return
end
function r8_sign ( x )

!*****************************************************************************80
!
!! r8_sign() returns the sign of an R8.
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
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! r8mat_print() prints an R8MAT.
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
!! r8mat_print_some() prints some of an R8MAT.
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
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! r8vec_print() prints an R8VEC.
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
subroutine r8vec2_print ( n, a1, a2, title )

!*****************************************************************************80
!
!! r8vec2_print() prints an R8VEC2.
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

