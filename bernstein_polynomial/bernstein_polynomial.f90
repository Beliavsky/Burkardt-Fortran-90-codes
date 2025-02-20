subroutine bernstein_matrix ( n, a )

!*****************************************************************************80
!
!! bernstein_matrix() returns the Bernstein matrix.
!
!  Discussion:
!
!    The Bernstein matrix of order N is an NxN matrix A which can be used to
!    transform a vector of power basis coefficients C representing a polynomial 
!    P(X) to a corresponding Bernstein basis coefficient vector B:
!
!      B = A * C
!
!    The N power basis vectors are ordered as (1,X,X^2,...X^(N-1)) and the N 
!    Bernstein basis vectors as ((1-X)^(N-1), X*(1-X)^(N-2),...,X^(N-1)).
!
!    For N = 5, the matrix has the form:
!
!      1 -4   6  -4  1
!      0  4 -12  12 -4
!      0  0   6 -12  6
!      0  0   0   4 -4
!      0  0   0   0  1
!
!    and the numbers in each column represent the coefficients in the power
!    series expansion of a Bernstein polynomial, so that 
!
!      B(5,4) = - 4 x^4 + 12 x^3 - 12 x^2 + 4 x
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the order of the matrix.
!
!  Output:
!
!    real ( kind = rk ) A(N,N), the Bernstein matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,n)
  integer i0
  integer j0
  integer n0
  real ( kind = rk ) r8_choose
  real ( kind = rk ) r8_mop

  a(1:n,1:n) = 0.0D+00

  n0 = n - 1

  do j0 = 0, n0
    do i0 = 0, j0
      a(i0+1,j0+1) = r8_mop ( j0 - i0 ) * r8_choose ( n0 - i0, j0 - i0 ) &
        * r8_choose ( n0, i0 )
    end do
  end do

  return
end
subroutine bernstein_matrix_determinant ( n, value )

!*****************************************************************************80
!
!! bernstein_matrix_determinant() returns the determinant of the BERNSTEIN matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the order of the matrix.
!
!  Output:
!
!    real ( kind = rk ) VALUE, the determinant.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer n
  real ( kind = rk ) r8_choose
  real ( kind = rk ) value

  value = 1.0D+00
  do i = 0, n - 1
    value = value * r8_choose ( n - 1, i )
  end do

  return
end
subroutine bernstein_matrix_inverse ( n, a )

!*****************************************************************************80
!
!! bernstein_matrix_inverse() returns the inverse Bernstein matrix.
!
!  Discussion:
!
!    The inverse Bernstein matrix of order N is an NxN matrix A which can 
!    be used to transform a vector of Bernstein basis coefficients B
!    representing a polynomial P(X) to a corresponding power basis 
!    coefficient vector C:
!
!      C = A * B
!
!    The N power basis vectors are ordered as (1,X,X^2,...X^(N-1)) and the N 
!    Bernstein basis vectors as ((1-X)^(N-1), X*(1-X)^(N-2),...,X^(N-1)).
!
!    For N = 5, the matrix has the form:
!
!      1   1    1    1   1
!      0  1/4  1/2  3/4  1
!      0   0   1/6  1/2  1
!      0   0    0   1/4  1
!      0   0    0    0   1
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the order of the matrix.
!
!  Output:
!
!    real ( kind = rk ) A(N,N), the inverse Bernstein matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,n)
  integer i0
  integer j0
  integer n0
  real ( kind = rk ) r8_choose

  a(1:n,1:n) = 0.0D+00

  n0 = n - 1

  do j0 = 0, n0
    do i0 = 0, j0
      a(i0+1,j0+1) = r8_choose ( j0, i0 ) / r8_choose ( n0, i0 )
    end do
  end do

  return
end
subroutine bernstein_poly_01 ( n, x, bern )

!*****************************************************************************80
!
!! bernstein_poly_01() evaluates the Bernstein polynomials based in [0,1].
!
!  Discussion:
!
!    The Bernstein polynomials are assumed to be based on [0,1].
!
!    The formula is:
!
!      B(N,I)(X) = [N!/(I!*(N-I)!)] * (1-X)^(N-I) * X^I
!
!  First values:
!
!    B(0,0)(X) = 1
!
!    B(1,0)(X) =      1-X
!    B(1,1)(X) =                X
!
!    B(2,0)(X) =     (1-X)^2
!    B(2,1)(X) = 2 * (1-X)    * X
!    B(2,2)(X) =                X^2
!
!    B(3,0)(X) =     (1-X)^3
!    B(3,1)(X) = 3 * (1-X)^2 * X
!    B(3,2)(X) = 3 * (1-X)   * X^2
!    B(3,3)(X) =               X^3
!
!    B(4,0)(X) =     (1-X)^4
!    B(4,1)(X) = 4 * (1-X)^3 * X
!    B(4,2)(X) = 6 * (1-X)^2 * X^2
!    B(4,3)(X) = 4 * (1-X)   * X^3
!    B(4,4)(X) =               X^4
!
!  Special values:
!
!    B(N,I)(X) has a unique maximum value at X = I/N.
!
!    B(N,I)(X) has an I-fold zero at 0 and and N-I fold zero at 1.
!
!    B(N,I)(1/2) = C(N,K) / 2^N
!
!    For a fixed X and N, the polynomials add up to 1:
!
!      Sum ( 0 <= I <= N ) B(N,I)(X) = 1
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the degree of the Bernstein polynomials 
!    to be used.  For any N, there is a set of N+1 Bernstein polynomials,
!    each of degree N, which form a basis for polynomials on [0,1].
!
!    real ( kind = rk ) X, the evaluation point.
!
!  Output:
!
!    real ( kind = rk ) BERN(0:N), the values of the N+1 
!    Bernstein polynomials at X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) bern(0:n)
  integer i
  integer j
  real ( kind = rk ) x

  if ( n == 0 ) then
 
    bern(0) = 1.0D+00
 
  else if ( 0 < n ) then
 
    bern(0) = 1.0D+00 - x
    bern(1) = x
 
    do i = 2, n
      bern(i) = x * bern(i-1)
      do j = i - 1, 1, -1
        bern(j) =             x   * bern(j-1) &
                + ( 1.0D+00 - x ) * bern(j)
      end do
      bern(0) = ( 1.0D+00 - x ) * bern(0)
    end do
 
  end if
 
  return
end
subroutine bernstein_poly_01_matrix ( m, n, x, b )

!*****************************************************************************80
!
!! bernstein_poly_01_matrix() evaluates Bernstein polynomials on [0,1].
!
!  Discussion:
!
!    The Bernstein polynomials are assumed to be based on [0,1].
!
!  Formula:
!
!    B(N,I)(X) = [N!/(I!*(N-I)!)] * (1-X)^(N-I) * X^I
!
!  First values:
!
!    B(0,0)(X) = 1
!
!    B(1,0)(X) =      1-X
!    B(1,1)(X) =                X
!
!    B(2,0)(X) =     (1-X)^2
!    B(2,1)(X) = 2 * (1-X)   * X
!    B(2,2)(X) =                X^2
!
!    B(3,0)(X) =     (1-X)^3
!    B(3,1)(X) = 3 * (1-X)^2 * X
!    B(3,2)(X) = 3 * (1-X)   * X^2
!    B(3,3)(X) =               X^3
!
!    B(4,0)(X) =     (1-X)^4
!    B(4,1)(X) = 4 * (1-X)^3 * X
!    B(4,2)(X) = 6 * (1-X)^2 * X^2
!    B(4,3)(X) = 4 * (1-X)   * X^3
!    B(4,4)(X) =               X^4
!
!  Special values:
!
!    B(N,I)(X) has a unique maximum value at X = I/N.
!
!    B(N,I)(X) has an I-fold zero at 0 and and N-I fold zero at 1.
!
!    B(N,I)(1/2) = C(N,K) / 2^N
!
!    For a fixed X and N, the polynomials add up to 1:
!
!      Sum ( 0 <= I <= N ) B(N,I)(X) = 1
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, the number of evaluation points.
!
!    integer N, the degree of the Bernstein polynomials to 
!    be used.  For any N, there is a set of N+1 Bernstein polynomials,
!    each of degree N, which form a basis for polynomials on [0,1].
!
!    real ( kind = rk ) X(M), the evaluation points.
!
!  Output:
!
!    real ( kind = rk ) B(M,1:N+1), the values of the N+1 Bernstein 
!    polynomials at X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) b(m,n+1)
  integer i
  integer j
  real ( kind = rk ) x(m)

  if ( n == 0 ) then
 
    b(1:m,1) = 1.0D+00
 
  else if ( 0 < n ) then
 
    b(1:m,1) = 1.0D+00 - x(1:m)
    b(1:m,2) = x(1:m)
 
    do i = 2, n
      b(1:m,i+1) = x(1:m) * b(1:m,i)
      do j = i - 1, 1, -1
        b(1:m,j+1) = x(1:m) * b(1:m,j) + ( 1.0D+00 - x(1:m) ) * b(1:m,j+1)
      end do
      b(1:m,1) = ( 1.0D+00 - x(1:m) ) * b(1:m,1)
    end do
 
  end if
 
  return
end
subroutine bernstein_poly_01_values ( n_data, n, k, x, b )

!*****************************************************************************80
!
!! bernstein_poly_01_values() returns some values of the Bernstein polynomials.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N_DATA: the user sets N_DATA to 0 on the first call.
!
!  Output:
!
!    integer N_DATA: N_DATA is incremented and that test data is returned. 
!    When there is no more test data, N_DATA is set to 0.
!
!    integer N, the degree of the polynomial.
!
!    integer K, the index of the polynomial.
! 
!    real ( kind = rk ) X, the argument of the polynomial.
!
!    real ( kind = rk ) B, the value of the polynomial B(N,K)(X).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nmax = 15

  real ( kind = rk ) b
  real ( kind = rk ), save, dimension ( nmax ) :: b_vec = (/ &
    1.0D+00, &
    0.75D+00,       0.25D+00, &
    0.5625D+00,     0.3750D+00,   0.0625D+00, &
    0.421875D+00,   0.421875D+00, 0.140625D+00,  0.015625D+00, &
    0.31640625D+00, 0.421875D+00, 0.2109375D+00, 0.046875D+00, 0.00390625D+00 /)
  integer k
  integer, save, dimension ( nmax ) :: k_vec = (/ &
    0, &
    0, 1, &
    0, 1, 2, &
    0, 1, 2, 3, &
    0, 1, 2, 3, 4 /)
  integer n
  integer n_data
  integer, save, dimension ( nmax ) :: n_vec = (/ &
    0, &
    1, 1, &
    2, 2, 2, &
    3, 3, 3, 3, &
    4, 4, 4, 4, 4 /)
  real ( kind = rk ) x
  real ( kind = rk ), save, dimension ( nmax ) :: x_vec = (/ &
    0.25D+00, &
    0.25D+00, 0.25D+00, &
    0.25D+00, 0.25D+00, 0.25D+00, &
    0.25D+00, 0.25D+00, 0.25D+00, 0.25D+00, &
    0.25D+00, 0.25D+00, 0.25D+00, 0.25D+00, 0.25D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( nmax < n_data ) then
    n_data = 0
    n = 0
    k = 0
    x = 0.0D+00
    b = 0.0D+00
  else
    n = n_vec(n_data)
    k = k_vec(n_data)
    x = x_vec(n_data)
    b = b_vec(n_data)
  end if

  return
end
subroutine bernstein_poly_ab ( n, a, b, x, bern )

!*****************************************************************************80
!
!! bernstein_poly_ab() evaluates the Bernstein polynomials based in [A,B].
!
!  Discussion:
!
!    The formula is:
!
!      BERN(N,I)(X) = [N!/(I!*(N-I)!)] * (B-X)^(N-I) * (X-A)^I / (B-A)^N
!
!  First values:
!
!    B(0,0)(X) =   1
!
!    B(1,0)(X) = (      B-X                ) / (B-A)
!    B(1,1)(X) = (                 X-A     ) / (B-A)
!
!    B(2,0)(X) = (     (B-X)^2             ) / (B-A)^2
!    B(2,1)(X) = ( 2 * (B-X)    * (X-A)    ) / (B-A)^2
!    B(2,2)(X) = (                (X-A)^2  ) / (B-A)^2
!
!    B(3,0)(X) = (     (B-X)^3             ) / (B-A)^3
!    B(3,1)(X) = ( 3 * (B-X)^2  * (X-A)    ) / (B-A)^3
!    B(3,2)(X) = ( 3 * (B-X)    * (X-A)^2  ) / (B-A)^3
!    B(3,3)(X) = (                (X-A)^3  ) / (B-A)^3
!
!    B(4,0)(X) = (     (B-X)^4             ) / (B-A)^4
!    B(4,1)(X) = ( 4 * (B-X)^3  * (X-A)    ) / (B-A)^4
!    B(4,2)(X) = ( 6 * (B-X)^2  * (X-A)^2  ) / (B-A)^4
!    B(4,3)(X) = ( 4 * (B-X)    * (X-A)^3  ) / (B-A)^4
!    B(4,4)(X) = (                (X-A)^4  ) / (B-A)^4
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the degree of the Bernstein polynomials 
!    to be used.  For any N, there is a set of N+1 Bernstein polynomials, 
!    each of degree N, which form a basis for polynomials on [A,B].
!
!    real ( kind = rk ) A, B, the endpoints of the interval on which the
!    polynomials are to be based.  A and B should not be equal.
!
!    real ( kind = rk ) X, the point at which the polynomials 
!    are to be evaluated.
!
!  Output:
!
!    real ( kind = rk ) BERN(0:N), the values of the N+1
!    Bernstein polynomials at X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) bern(0:n)
  integer i
  integer j
  real ( kind = rk ) x

  if ( b == a ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'bernstein_poly_ab - Fatal error!'
    write ( *, '(a,g14.6)' ) '  A = B = ', a
    stop 1
  end if

  if ( n == 0 ) then
 
    bern(0) = 1.0D+00
 
  else if ( 0 < n ) then
 
    bern(0) = ( b - x ) / ( b - a )
    bern(1) = ( x - a ) / ( b - a )
 
    do i = 2, n
      bern(i) = ( x - a ) * bern(i-1) / ( b - a )
      do j = i - 1, 1, -1
        bern(j) = ( ( b - x     ) * bern(j)     &
                  + (     x - a ) * bern(j-1) ) &
                  / ( b     - a )
      end do
      bern(0) = ( b - x ) * bern(0) / ( b - a )
    end do
 
  end if
 
  return
end
subroutine bernstein_poly_ab_approx ( n, a, b, ydata, nval, xval, yval )

!*****************************************************************************80
!
!! bernstein_poly_ab_approx() evaluates Bernstein approximant to F(X) on [A,B].
!
!  Formula:
!
!    BPAB(F)(X) = sum ( 0 <= I <= N ) F(X(I)) * B_BASE(I,X)
!
!    where
!
!      X(I) = ( ( N - I ) * A + I * B ) / N
!      B_BASE(I,X) is the value of the I-th Bernstein basis polynomial at X.
!
!  Discussion:
!
!    The Bernstein polynomial BPAB(F) for F(X) over [A,B] is an approximant, 
!    not an interpolant; in other words, its value is not guaranteed to equal
!    that of F at any particular point.  However, for a fixed interval
!    [A,B], if we let N increase, the Bernstein polynomial converges
!    uniformly to F everywhere in [A,B], provided only that F is continuous.
!    Even if F is not continuous, but is bounded, the polynomial converges
!    pointwise to F(X) at all points of continuity.  On the other hand,
!    the convergence is quite slow compared to other interpolation
!    and approximation schemes.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Input:
!
!    integer N, the degree of the Bernstein polynomial
!    to be used.  N must be at least 0.
!
!    real ( kind = rk ) A, B, the endpoints of the interval on which the
!    approximant is based.  A and B should not be equal.
!
!    real ( kind = rk ) YDATA(0:N), the data values at N+1 equally
!    spaced points in [A,B].  If N = 0, then the evaluation point should
!    be 0.5 * ( A + B).  Otherwise, evaluation point I should be
!    ( (N-I)*A + I*B ) / N ).
!
!    integer NVAL, the number of points at which the
!    approximant is to be evaluated.
!
!    real ( kind = rk ) XVAL(NVAL), the point at which the Bernstein 
!    polynomial approximant is to be evaluated.  The entries of XVAL do not 
!    have to lie in the interval [A,B].
!
!  Output:
!
!    real ( kind = rk ) YVAL(NVAL), the values of the Bernstein 
!    polynomial approximant for F, based in [A,B], evaluated at XVAL.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nval

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) bvec(0:n)
  integer i
  real ( kind = rk ) xval(nval)
  real ( kind = rk ) ydata(0:n)
  real ( kind = rk ) yval(nval)

  do i = 1, nval
!
!  Evaluate the Bernstein basis polynomials at XVAL.
!
    call bernstein_poly_ab ( n, a, b, xval(i), bvec )
!
!  Now compute the sum of YDATA(I) * BVEC(I).
!
    yval(i) = dot_product ( ydata(0:n), bvec(0:n) )

  end do

  return
end
subroutine bernstein_to_legendre ( n, a )

!*****************************************************************************80
!
!! bernstein_to_legendre() returns the Bernstein-to-Legendre matrix.
!
!  Discussion:
!
!    The Legendre polynomials are often defined on [-1,+1], while the
!    Bernstein polynomials are defined on [0,1].  For this function,
!    the Legendre polynomials have been shifted to share the [0,1]
!    interval of definition.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the maximum degree of the polynomials.
!
!
!    real ( kind = rk ) A(0:N,0:N), the Bernstein-to-Legendre matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(0:n,0:n)
  integer i
  integer j
  integer k
  real ( kind = rk ) r8_choose
  real ( kind = rk ) r8_mop

  a(0:n,0:n) = 0.0D+00

  do i = 0, n
    do j = 0, n
      do k = 0, i
        a(i,j) = a(i,j) &
          + r8_mop ( i + k ) * r8_choose ( i, k ) ** 2 &
          / r8_choose ( n + i, j + k )
      end do
      a(i,j) = a(i,j) * r8_choose ( n, j ) &
        * real ( 2 * i + 1, kind = rk ) / real ( n + i + 1, kind = rk )
    end do
  end do

  return
end
subroutine bernstein_to_power ( n, a )

!*****************************************************************************80
!
!! bernstein_to_power() returns the Bernstein-to-Power matrix.
!
!  Discussion:
!
!    The Bernstein-to-Power matrix of degree N is an N+1xN+1 matrix A which can 
!    be used to transform the N+1 coefficients of a polynomial of degree N
!    from a vector B of Bernstein basis polynomial coefficients ((1-x)^n,...,x^n).
!    to a vector P of coefficients of the power basis (1,x,x^2,...,x^n).
!
!    If we are using N=4-th degree polynomials, the matrix has the form:
!
!      1   0   0   0  0
!     -4   4   0   0  0
!      6 -12   6   0  0
!     -4  12 -12   4  0
!      1  -4   6  -4  1
!
!   and a polynomial with the Bernstein basis representation
!     p(x) = 3/4 * b(4,1) + 1/2 b(4,2)
!   whose Bernstein coefficient vector is
!     B = ( 0, 3/4, 1/2, 0, 0 )
!   will have the Bernstein basis coefficients 
!     P = A * B = ( 0, 3, -6, 3, 0 ).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the degree of the polynomials.
!
!  Output:
!
!    real ( kind = rk ) A(0:N,0:N), the Bernstein-to-Power matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(0:n,0:n)
  integer i
  integer j
  real ( kind = rk ) r8_choose
  real ( kind = rk ) r8_mop

  a(0:n,0:n) = 0.0D+00

  do j = 0, n
    do i = 0, j
      a(n-i,n-j) = r8_mop ( j - i ) * r8_choose ( n - i, j - i ) &
        * r8_choose ( n, i )
    end do
  end do

  return
end
subroutine bernstein_vandermonde ( n, v )

!*****************************************************************************80
!
!! bernstein_vandermonde() returns the Bernstein Vandermonde matrix.
!
!  Discussion:
!
!    The Bernstein Vandermonde matrix of order N is constructed by
!    evaluating the N Bernstein polynomials of degree N-1 at N equally
!    spaced points between 0 and 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the order of the matrix.
!
!  Output:
!
!    real ( kind = rk ) A(N,N), the Bernstein Vandermonde matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) b(n)
  integer i
  real ( kind = rk ) v(n,n)
  real ( kind = rk ) x

  if ( n == 1 ) then
    v(1,1) = 1.0D+00
    return
  end if

  do i = 1, n
    x = real ( i - 1, kind = rk ) / real ( n - 1, kind = rk )
    call bernstein_poly_01 ( n - 1, x, b );
    v(i,1:n) = b(1:n)
  end do

  return
end
subroutine legendre_to_bernstein ( n, a )

!*****************************************************************************80
!
!! legendre_to_bernstein() returns the Legendre-to-Bernstein matrix.
!
!  Discussion:
!
!    The Legendre polynomials are often defined on [-1,+1], while the
!    Bernstein polynomials are defined on [0,1].  For this function,
!    the Legendre polynomials have been shifted to share the [0,1]
!    interval of definition.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the maximum degree of the polynomials.
!
!  Output:
!
!    real A(0:N,0:N), the Legendre-to-Bernstein matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(0:n,0:n)
  integer i
  integer j
  integer k
  real ( kind = rk ) r8_choose
  real ( kind = rk ) r8_mop

  a(0:n,0:n) = 0.0D+00

  do i = 0, n
    do j = 0, n
      do k = max ( 0, i + j - n ), min ( i, j )
        a(i,j) = a(i,j) &
          + r8_mop ( j + k ) * r8_choose ( j, k ) ** 2 &
          * r8_choose ( n - j, i - k )
      end do
      a(i,j) = a(i,j) / r8_choose ( n, i )
    end do
  end do

  return
end
subroutine power_to_bernstein ( n, a )

!*****************************************************************************80
!
!! power_to_bernstein() returns the Power-to-Bernstein matrix.
!
!  Discussion:
!
!    The Power-to-Bernstein matrix of degree N is an N+1xN+1 matrix A which can 
!    be used to transform the N+1 coefficients of a polynomial of degree N
!    from a vector P of coefficients of the power basis (1,x,x^2,...,x^n)
!    to a vector B of Bernstein basis polynomial coefficients ((1-x)^n,...,x^n).
!
!    If we are using N=4-th degree polynomials, the matrix has the form:
!
!          1   0    0    0   0
!          1  1/4   0    0   0
!      A = 1  1/2  1/6   0   0
!          1  3/4  1/2  1/4  1
!          1   1    1    1   1
!
!   and a polynomial 
!     p(x) = 3x - 6x^2 + 3x^3
!   whose power coefficient vector is
!     P = ( 0, 3, -6, 3, 0 )
!   will have the Bernstein basis coefficients 
!     B = A * P = ( 0, 3/4, 1/2, 0, 0 ).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the degree of the polynomials.
!
!  Output:
!
!    real ( kind = rk ) A(0:N,0:N), the Power-to-Bernstein matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(0:n,0:n)
  integer i
  integer j
  real ( kind = rk ) r8_choose

  a(0:n,0:n) = 0.0D+00

  do j = 0, n
    do i = 0, j
      a(n-i,n-j) = r8_choose ( j, i ) / r8_choose ( n, i )
    end do
  end do

  return
end
function r8_choose ( n, k )

!*****************************************************************************80
!
!! r8_choose() computes the binomial coefficient C(N,K) as an R8.
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
!    02 September 2021
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
!  Input:
!
!    integer N, K, are the values of N and K.
!
!  Output:
!
!    real ( kind = rk ) r8_choose, the number of combinations of N
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
function r8_mop ( i )

!*****************************************************************************80
!
!! r8_mop() returns the I-th power of -1 as an R8.
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
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer I, the power of -1.
!
!  Output:
!
!    real ( kind = rk ) r8_mop, the I-th power of -1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  real ( kind = rk ) r8_mop

  if ( mod ( i, 2 ) == 0 ) then
    r8_mop = + 1.0D+00
  else
    r8_mop = - 1.0D+00
  end if

  return
end
subroutine r8mat_is_identity ( n, a, error_frobenius )

!*****************************************************************************80
!
!! r8mat_is_identity() determines if an R8MAT is the identity.
!
!  Discussion:
!
!    An R8MAT is a matrix of real ( kind = rk ) values.
!
!    The routine returns the Frobenius norm of A - I.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the order of the matrix.
!
!    real ( kind = rk ) A(N,N), the matrix.
!
!  Output:
!
!    real ( kind = rk ) ERROR_FROBENIUS, the Frobenius norm
!    of the difference matrix A - I, which would be exactly zero
!    if A were the identity matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) error_frobenius
  integer i
  integer j

  error_frobenius = 0.0D+00

  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        error_frobenius = error_frobenius + ( a(i,j) - 1.0D+00 )**2
      else
        error_frobenius = error_frobenius + a(i,j)**2
      end if
    end do 
  end do

  error_frobenius = sqrt ( error_frobenius )

  return
end
function r8mat_norm_fro ( m, n, a )

!*****************************************************************************80
!
!! r8mat_norm_fro() returns the Frobenius norm of an M by N R8MAT.
!
!  Discussion:
!
!    An R8MAT is a matrix of real ( kind = rk ) values.
!
!    The Frobenius norm is defined as
!
!      r8mat_norm_fro = sqrt (
!        sum ( 1 <= I <= M ) Sum ( 1 <= J <= N ) A(I,J)^2 )
!
!    The matrix Frobenius-norm is not derived from a vector norm, but
!    is compatible with the vector L2 norm, so that:
!
!      r8vec_norm_l2 ( A*x ) <= r8mat_norm_fro ( A ) * r8vec_norm_l2 ( x ).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, N, the order of the matrix.
!
!    real ( kind = rk ) A(M,N), the matrix.
!
!  Output:
!
!    real ( kind = rk ) r8mat_norm_fro, the Frobenius norm of A.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) r8mat_norm_fro

  r8mat_norm_fro = sqrt ( sum ( a(1:m,1:n)**2 ) )

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! r8mat_print() prints an R8MAT.
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
!    02 September 2021
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
!    real ( kind = rk ) A(M,N), the matrix.
!
!    character ( len = * ) TITLE, a title.
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
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, N, the number of rows and columns.
!
!    real ( kind = rk ) A(M,N), an M by N matrix to be printed.
!
!    integer ILO, JLO, the first row and column to print.
!
!    integer IHI, JHI, the last row and column to print.
!
!    character ( len = * ) TITLE, a title.
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
subroutine r8vec_linspace ( n, a_first, a_last, a )

!*****************************************************************************80
!
!! r8vec_linspace() creates a vector of linearly spaced values.
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
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of entries in the vector.
!
!    real ( kind = rk ) A_FIRST, A_LAST, the first and last entries.
!
!  Output:
!
!    real ( kind = rk ) A(N), a vector of linearly spaced data.
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
subroutine timestamp ( )

!*****************************************************************************80
!
!! timestamp() prints the current YMDHMS date as a time stamp.
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
!    02 September 2021
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
 
