subroutine bernstein_to_legendre01 ( n, bcoef, lcoef )

!*****************************************************************************80
!
!! bernstein_to_legendre01() converts from Bernstein to Legendre01 form.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 February 2024
!
!  Author:
!
!    John Burkardt.
!
!  Input:
!
!    integer n: the order of the polynomial.
!
!    real bcoef(0:n): the Bernstein coefficients of the polynomial.
! 
!  Output:
!
!    real lcoef(0:n): the Legendre01 coefficients of the polynomial.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ), allocatable :: A(:,:)
  real ( kind = rk ) bcoef(0:n)
  real ( kind = rk ) lcoef(0:n)

  if ( n < 0 ) then
    return
  end if
!
!  Get the matrix.
!
  allocate ( A(0:n,0:n) )

  call bernstein_to_legendre01_matrix ( n, A )
!
!  Apply the linear transformation.
!
  lcoef = matmul ( A, bcoef )

  deallocate ( A )

  return
end
subroutine bernstein_to_legendre01_matrix ( n, A )

!*****************************************************************************80
!
!! bernstein_to_legendre01_matrix() returns the Bernstein-to-Legendre matrix.
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
!    25 February 2024
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
!    real ( kind = rk ) A(0:N,0:N), the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) A(0:n,0:n)
  integer i
  integer i4_choose
  integer j
  integer k

  A(0:n,0:n) = 0.0D+00

  do i = 0, n
    do j = 0, n
      do k = 0, i
        A(i,j) = A(i,j) &
          + ( -1.0D+00 ) ** ( i + k ) * i4_choose ( i, k ) ** 2 &
          / i4_choose ( n + i, j + k )
      end do
      A(i,j) = A(i,j) * i4_choose ( n, j ) &
        * real ( 2 * i + 1, kind = rk ) / real ( n + i + 1, kind = rk )
    end do
  end do

  return
end
subroutine bernstein_to_monomial ( n, bcoef, mcoef )

!*****************************************************************************80
!
!! bernstein_to_monomial() converts from Bernstein to monomial form.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 February 2024
!
!  Author:
!
!    John Burkardt.
!
!  Input:
!
!    integer n: the order of the polynomial.
!
!    real bcoef(0:n): the Bernstein coefficients of the polynomial.
! 
!  Output:
!
!    real mcoef(0:n): the monomial coefficients of the polynomial.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ), allocatable :: A(:,:)
  real ( kind = rk ) bcoef(0:n)
  real ( kind = rk ) mcoef(0:n)

  if ( n < 0 ) then
    return
  end if
!
!  Get the matrix.
!
  allocate ( A(0:n,0:n) )

  call bernstein_to_monomial_matrix ( n + 1, A )
!
!  Apply the linear transformation.
!
  mcoef = matmul ( A, bcoef )

  deallocate ( A )

  return
end
subroutine bernstein_to_monomial_matrix ( n, A )

!*****************************************************************************80
!
!! bernstein_to_monomial_matrix() converts Bernstein to monomial form.
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
!    Bernstein basis vectors as ((1-X)^(N-1), X*(1_X)^(N-2),...,X^(N-1)).
!
!  Example:
!
!    N = 5
!
!    1    -4     6    -4     1
!    0     4   -12    12    -4
!    0     0     6   -12     6
!    0     0     0     4    -4
!    0     0     0     0     1
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 July 2011
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
!    real ( kind = rk ) A(N,N), the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) A(n,n)
  integer i0
  integer i4_choose
  integer j0
  integer n0

  A(1:n,1:n) = 0.0D+00

  n0 = n - 1

  do j0 = 0, n0
    do i0 = 0, j0
      A(i0+1,j0+1) = ( -1 ) ** ( j0 - i0 ) * i4_choose ( n0 - i0, j0 - i0 ) &
        * i4_choose ( n0, i0 )
    end do
  end do

  return
end
subroutine chebyshev_to_monomial ( n, ccoef, mcoef )

!*****************************************************************************80
!
!! chebyshev_to_monomial() converts from Chebyshev to monomial form.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 February 2024
!
!  Author:
!
!    Original Fortran77 version by Fred Krogh.
!    This version by John Burkardt.
!
!  Input:
!
!    integer n: the order of the polynomial.
!
!    real ( kind = rk ) ccoef(0:n): the Chebyshev coefficients of the polynomial.
! 
!  Output:
!
!    real ( kind = rk ) mcoef(0:n): the monomial coefficients of the polynomial.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) ccoef(0:n)
  integer i
  integer j
  real ( kind = rk ) mcoef(0:n)
  real ( kind = rk ) tp

  if ( n < 0 ) then
    return
  end if

  tp = 1.0D+00

  mcoef(0:n) = ccoef(0:n)

  do j = 0, n - 2
    do i = n - 2, j, -1
      mcoef(i) = mcoef(i) - mcoef(i+2)
    end do
    mcoef(j+1) = 0.5D+00 * mcoef(j+1)
    mcoef(j) = tp * mcoef(j)
    tp = 2.0D+00 * tp
  end do

  mcoef(n) = tp * mcoef(n)
  if ( 0 < n ) then
    mcoef(n-1) = tp * mcoef(n-1)
  end if

  return
end
subroutine chebyshev_to_monomial_matrix ( n, A )

!*****************************************************************************80
!
!! chebyshev_to_monomial_matrix() converts a Chebyshev T polynomial to monomial.
!
!  Discussion
!
!    This is the Chebyshev T matrix, associated with the Chebyshev
!    "T" polynomials, or Chebyshev polynomials of the first kind.
!
!  Example:
!
!    N = 11
!
!    1  .  -1    .    1    .   -1    .     1    .    -1
!    .  1   .   -3    .    5    .   -7     .    9     .
!    .  .   2    .   -8    .   18    .   -32    .    50
!    .  .   .    4    .  -20    .   56     . -120     .
!    .  .   .    .    8    .  -48    .   160    .  -400
!    .  .   .    .    .   16    . -112     .  432     .
!    .  .   .    .    .    .   32    .  -256    .  1120
!    .  .   .    .    .    .    .   64     . -576     .
!    .  .   .    .    .    .    .    .   128    . -1280
!    .  .   .    .    .    .    .    .     .  256     .
!    .  .   .    .    .    .    .    .     .    .   512
!
!  Properties:
!
!    A is generally not symmetric: A' /= A.
!
!    A is integral, therefore det ( A ) is integral, and 
!    det ( A ) * inverse ( A ) is integral.
!
!    A is reducible.
!
!    A is lower triangular.
!
!    Each row of A sums to 1.
!
!    det ( A ) = 2^( (N-1) * (N-2) / 2 )
!
!    A is not normal: A' * A /= A * A'.
!
!    For I = 1:
!
!      LAMBDA(1) = 1
!
!    For 1 < I
!
!      LAMBDA(I) = 2^(I-2)
!
!    The family of matrices is nested as a function of N.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 April 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N: the order of the matrix.
!
!  Output:
!
!    real ( kind = rk ) A(N,N): the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) A(n,n)
  integer i
  integer j

  if ( n <= 0 ) then
    return
  end if

  A(1:n,1:n) = 0.0D+00

  A(1,1) = 1.0D+00

  if ( n == 1 ) then
    return
  end if

  A(2,2) = 1.0D+00

  do j = 3, n
    do i = 1, n
      if ( i == 1 ) then
        A(i,j) = - A(i,j-2)
      else
        A(i,j) = 2.0D+00 * A(i-1,j-1) - A(i,j-2)
      end if
    end do
  end do

  return
end
subroutine gegenbauer_to_monomial ( n, alpha, gcoef, mcoef )

!*****************************************************************************80
!
!! gegenbauer_to_monomial() converts from Gegenbauer to monomial form.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 February 2024
!
!  Author:
!
!    John Burkardt.
!
!  Input:
!
!    integer n: the order of the polynomial.
!
!    real alpha: the Gegenbauer parameter.
!
!    real gcoef(0:n): the Gegenbauer coefficients of the polynomial.
! 
!  Output:
!
!    real mcoef(0:n): the monomial coefficients of the polynomial.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ), allocatable :: A(:,:)
  real ( kind = rk ) alpha
  real ( kind = rk ) gcoef(0:n)
  real ( kind = rk ) mcoef(0:n)

  if ( n < 0 ) then
    return
  end if
!
!  Get the matrix.
!
  allocate ( A(0:n,0:n) )

  call gegenbauer_to_monomial_matrix ( n + 1, alpha, A )
!
!  Apply the linear transformation.
!
  mcoef = matmul ( A, gcoef )

  deallocate ( A )

  return
end
subroutine gegenbauer_to_monomial_matrix ( n, alpha, A )

!*****************************************************************************80
!
!! gegenbauer_to_monomial_matrix(): Gegenbauer to monomial conversion matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 April 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N: the order of A.
!
!    real alpha: the parameter.
!
!  Output:
!
!    real ( kind = rk ) A(N,N): the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) A(n,n)
  real ( kind = rk ) alpha
  real ( kind = rk ) c1
  real ( kind = rk ) c2
  integer i
  integer j
  integer nn

  if ( n <= 0 ) then
    return
  end if

  do j = 1, n
    do i = 1, n
      A(i,j) = 0.0D+00
    end do
  end do

  A(1,1) = 1.0D+00

  if ( n == 1 ) then
    return
  end if

  A(2,2) = 2.0D+00 * alpha
!
!  Perform convex sum.
!  Translating "(n+1) C(n+1) = 2 (n+alpha) x C(n) - ( n + 2 alpha - 1 ) C(n-1)"
!  drove me nuts, between indexing at 1 rather than 0, and dealing with
!  the interpretation of "n+1", because I now face the rare "off by 2" error!
!
  do j = 3, n
    nn = j - 2
    c1 = ( 2.0 * nn + 2.0 * alpha       ) / real ( nn + 1, kind = rk )
    c2 = (     - nn - 2.0 * alpha + 1.0 ) / real ( nn + 1, kind = rk )
    A(2:j,j)   =              c1 * A(1:j-1,j-1)
    A(1:j-2,j) = A(1:j-2,j) + c2 * A(1:j-2,j-2)
  end do

  return
end
subroutine hermite_to_monomial ( n, hcoef, mcoef )

!*****************************************************************************80
!
!! hermite_to_monomial() converts from Hermite to monomial form.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 February 2024
!
!  Author:
!
!    John Burkardt.
!
!  Input:
!
!    integer n: the order of the polynomial.
!
!    real hcoef(0:n): the Hermite coefficients of the polynomial.
! 
!  Output:
!
!    real mcoef(0:n): the monomial coefficients of the polynomial.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ), allocatable :: A(:,:)
  real ( kind = rk ) hcoef(0:n)
  real ( kind = rk ) mcoef(0:n)

  if ( n < 0 ) then
    return
  end if
!
!  Get the matrix.
!
  allocate ( A(0:n,0:n) )

  call hermite_to_monomial_matrix ( n + 1, A )
!
!  Apply the linear transformation.
!
  mcoef = matmul ( A, hcoef )

  deallocate ( A )

  return
end
subroutine hermite_to_monomial_matrix ( n, A )

!*****************************************************************************80
!
!! hermite_to_monomial_matrix() converts from Hermite to monomial form.
!
!  Example:
!
!    N = 11
!
!      1   .   -2      .      12     .   -120     .   1680      .  -30240
!      .   2    .    -12       .   120      . -1680      .  30240       .
!      .   .    4      .     -48     .    720     . -13440      .  302400
!      .   .    .      8       .  -160      .  3360      . -80640       .
!      .   .    .      .      16     .   -480     .  13440      . -403200
!      .   .    .      .       .    32      . -1344      .  48384       .
!      .   .    .      .       .     .     64     .  -3584      .  161280
!      .   .    .      .       .     .      .   128      .  -9216       .
!      .   .    .      .       .     .      .     .    256      .  -23040
!      .   .    .      .       .     .      .     .      .    512       .
!      .   .    .      .       .     .      .     .      .      .    1024
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 July 2000
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
!    real ( kind = rk ) A(N,N), the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) A(n,n)
  integer i
  integer j

  if ( n <= 0 ) then
    return
  end if

  A(1:n,1:n) = 0.0D+00

  A(1,1) = 1.0D+00

  if ( n == 1 ) then
    return
  end if

  A(2,2) = 2.0D+00

  if ( n == 2 ) then
    return
  end if

  do j = 3, n
    do i = 1, n
      if ( i == 1 ) then
        A(i,j) =                      - 2.0D+00 &
          * real ( j - 2, kind = rk ) * A(i,j-2)
      else
        A(i,j) = 2.0D+00 * a(i-1,j-1) - 2.0D+00 &
          * real ( j - 2, kind = rk ) * A(i,j-2)
      end if
    end do
  end do

  return
end
function i4_choose ( n, k )

!*****************************************************************************80
!
!! i4_choose() computes the binomial coefficient C(N,K) as an I4.
!
!  Discussion:
!
!    The value is calculated in such a way as to avoid overflow and
!    roundoff.  The calculation is done in integer arithmetic.
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
!    30 October 2014
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
!    integer I4_CHOOSE, the number of combinations of N
!    things taken K at a time.
!
  implicit none

  integer i
  integer i4_choose
  integer k
  integer mn
  integer mx
  integer n
  integer value

  mn = min ( k, n - k )
  mx = max ( k, n - k )

  if ( mn < 0 ) then

    value = 0

  else if ( mn == 0 ) then

    value = 1

  else

    value = mx + 1

    do i = 2, mn
      value = ( value * ( mx + i ) ) / i
    end do

  end if

  i4_choose = value

  return
end
subroutine laguerre_to_monomial ( n, lcoef, mcoef )

!*****************************************************************************80
!
!! laguerre_to_monomial() converts from Laguerre to monomial form.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 February 2024
!
!  Author:
!
!    John Burkardt.
!
!  Input:
!
!    integer n: the order of the polynomial.
!
!    real lcoef(0:n): the Laguerre coefficients of the polynomial.
! 
!  Output:
!
!    real mcoef(0:n): the monomial coefficients of the polynomial.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ), allocatable :: A(:,:)
  real ( kind = rk ) lcoef(0:n)
  real ( kind = rk ) mcoef(0:n)

  if ( n < 0 ) then
    return
  end if
!
!  Get the matrix.
!
  allocate ( A(0:n,0:n) )

  call laguerre_to_monomial_matrix ( n + 1, A )
!
!  Apply the linear transformation.
!
  mcoef = matmul ( A, lcoef )

  deallocate ( A )

  return
end
subroutine laguerre_to_monomial_matrix ( n, A )

!*****************************************************************************80
!
!! laguerre_to_monomial_matrix() converts from laguerre to monomial form.
!
!  Example:
!
!    N = 8 (each column must be divided by the factor below it.)
!
!      1      1      2      6     24    120    720   5040
!      .     -1     -4    -18    -96   -600  -4320 -35280
!      .      .      1      9     72    600   5400  52920
!      .      .      .      1    -16   -200  -2400 -29400
!      .      .      .      .      1     25    450   7350
!      .      .      .      .      .     -1    -36   -882
!      .      .      .      .      .      .      1     49
!      .      .      .      .      .      .      .     -1
!
!     /1     /1     /2     /6    /24   /120   /720  /5040
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 February 2024
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    US Department of Commerce, 1964.
!
!  Input:
!
!    integer N, the order of the matrix.
!
!  Output:
!
!    real ( kind = rk ) A(N,N), the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) A(n,n)
  integer i
  integer j

  if ( n <= 0 ) then
    return
  end if

  A(1:n,1:n) = 0.0D+00

  A(1,1) = 1.0D+00

  if ( n == 1 ) then
    return
  end if

  A(1,2) = 1.0D+00
  A(2,2) = -1.0D+00

  if ( n == 2 ) then
    return
  end if

  do j = 3, n
    do i = 1, n
      if ( i == 1 ) then
        A(i,j) = ( real ( 2 * j - 3, kind = rk ) * A(i,j-1) &
                 + real (   - j + 2, kind = rk ) * A(i,j-2) ) &
                 / real (     j - 1, kind = rk )
      else
        A(i,j) = ( real ( 2 * j - 3, kind = rk ) * A(i,j-1) &
                 - real (         1, kind = rk ) * A(i-1,j-1) &
                 + real (   - j + 2, kind = rk ) * A(i,j-2) ) &
                 / real (     j - 1, kind = rk )
      end if
    end do
  end do

  return
end
subroutine legendre_to_monomial ( n, lcoef, mcoef )

!*****************************************************************************80
!
!! legendre_to_monomial() converts from Legendre to monomial form.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 February 2024
!
!  Author:
!
!    John Burkardt.
!
!  Input:
!
!    integer n: the order of the polynomial.
!
!    real lcoef(0:n): the Legendre coefficients of the polynomial.
! 
!  Output:
!
!    real mcoef(0:n): the monomial coefficients of the polynomial.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ), allocatable :: A(:,:)
  real ( kind = rk ) lcoef(0:n)
  real ( kind = rk ) mcoef(0:n)

  if ( n < 0 ) then
    return
  end if
!
!  Get the matrix.
!
  allocate ( A(0:n,0:n) )

  call legendre_to_monomial_matrix ( n + 1, A )
!
!  Apply the linear transformation.
!
  mcoef = matmul ( A, lcoef )

  deallocate ( A )

  return
end
subroutine legendre_to_monomial_matrix ( n, A )

!*****************************************************************************80
!
!! legendre_to_monomial_matrix() returns the LEGENDRE matrix.
!
!  Example:
!
!    N = 11 (each column must be divided by factor at bottom)
!
!     1    .    -1     .      3     .     -5      .      35     .   -63
!     .    1     .    -3      .    15      .    -25       .   315     .
!     .    .     3     .    -30     .    105      .   -1260     .  3465
!     .    .     .     5      .   -70      .    315       . -4620     .
!     .    .     .     .     35     .   -315      .    6930     .-30030
!     .    .     .     .      .    63      .   -693       . 18018     .
!     .    .     .     .      .     .    231      .  -12012     . 90090
!     .    .     .     .      .     .      .    429       .-25740     .
!     .    .     .     .      .     .      .      .    6435     -109395
!     .    .     .     .      .     .      .      .       . 12155     .
!     .    .     .     .      .     .      .      .       .     . 46189
!
!    /1   /1    /2    /2     /8    /8    /16    /16    /128  /128  /256
!
!  Properties:
!
!    A is generally not symmetric: A' /= A.
!
!    A is lower triangular.
!
!    The elements of each row sum to 1.
!
!    Because it has a constant row sum of 1,
!    A has an eigenvalue of 1, and
!    a (right) eigenvector of ( 1, 1, 1, ..., 1 ).
!
!    A is reducible.
!
!    The diagonals form a pattern of zero, positive, zero, negative.
!
!    The family of matrices is nested as a function of N.
!
!    A is not diagonally dominant.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 February 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the order of A.
!
!  Output:
!
!    real A(N,N), the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) A(n,n)
  integer i
  integer j

  A(1:n,1:n) = 0.0D+00

  if ( n <= 0 ) then
    return
  end if

  A(1,1) = 1.0D+00

  if ( n == 1 ) then
    return
  end if

  A(2,2) = 1.0D+00

  if ( n == 2 ) then
    return
  end if

  do j = 3, n
    do i = 1, n
      if ( i == 1 ) then
        A(i,j) = - ( j - 2 ) * A(i,j-2) &
                 / ( j - 1 )
      else
        A(i,j) = ( ( 2 * j - 3 ) * A(i-1,j-1) &
                 + (   - j + 2 ) * A(i,j-2) ) &
                 / (     j - 1 )
      end if
    end do
  end do

  return
end
subroutine legendre01_to_bernstein ( n, lcoef, bcoef )

!*****************************************************************************80
!
!! legendre01_to_bernstein() converts from Legendre01 to Bernstein form.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 February 2024
!
!  Author:
!
!    John Burkardt.
!
!  Input:
!
!    integer n: the order of the polynomial.
!
!    real lcoef(1:n+1): the Legendre01 coefficients of the polynomial.
! 
!  Output:
!
!    real bcoef(1:n+1): the Bernstein coefficients of the polynomial.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ), allocatable :: A(:,:)
  real ( kind = rk ) bcoef(0:n)
  real ( kind = rk ) lcoef(0:n)

  if ( n < 0 ) then
    return
  end if
!
!  Get the matrix.
!
  allocate ( A(0:n,0:n) )
  call legendre01_to_bernstein_matrix ( n, A )
!
!  Apply the linear transformation.
!
  bcoef = matmul ( A, lcoef )

  deallocate ( A )

  return
end
subroutine legendre01_to_bernstein_matrix ( n, a )

!*****************************************************************************80
!
!! legendre01_to_bernstein_matrix() returns the Legendre01-to-Bernstein matrix.
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
!    25 February 2024
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
!    real A(0:N,0:N), the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) A(0:n,0:n)
  integer i
  integer i4_choose
  integer j
  integer k

  A(0:n,0:n) = 0.0D+00

  do i = 0, n
    do j = 0, n
      do k = max ( 0, i + j - n ), min ( i, j )
        A(i,j) = A(i,j) &
          + ( -1.0D+00 ) ** ( j + k ) * i4_choose ( j, k ) ** 2 &
          * i4_choose ( n - j, i - k )
      end do
      A(i,j) = A(i,j) / i4_choose ( n, i )
    end do
  end do

  return
end
subroutine monomial_to_bernstein ( n, mcoef, bcoef )

!*****************************************************************************80
!
!! monomial_to_bernstein() converts from monomial to Bernstein form.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 February 2024
!
!  Author:
!
!    John Burkardt.
!
!  Input:
!
!    integer n: the order of the polynomial.
!
!    real mcoef(1:n+1): the monomial coefficients of the polynomial.
! 
!  Output:
!
!    real bcoef(1:n+1): the Bernstein coefficients of the polynomial.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ), allocatable :: A(:,:)
  real ( kind = rk ) bcoef(0:n)
  real ( kind = rk ) mcoef(0:n)

  if ( n < 0 ) then
    return
  end if
!
!  Get the matrix.
!
  allocate ( A(0:n,0:n) )
  call monomial_to_bernstein_matrix ( n + 1, A )
!
!  Apply the linear transformation.
!
  bcoef = matmul ( A, mcoef )

  deallocate ( A )

  return
end
subroutine monomial_to_bernstein_matrix ( n, A )

!*****************************************************************************80
!
!! monomial_to_bernstein_matrix() converts monomial to Bernsteim form.
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
!  Example:
!
!    N = 5
!
!   1.0000    1.0000    1.0000    1.0000    1.0000
!        0    0.2500    0.5000    0.7500    1.0000
!        0         0    0.1667    0.5000    1.0000
!        0         0         0    0.2500    1.0000
!        0         0         0         0    1.0000
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 July 2011
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
!    real ( kind = rk ) A(N,N), the inverse matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) A(n,n)
  integer i0
  integer i4_choose
  integer j0
  integer n0

  A(1:n,1:n) = 0.0D+00

  n0 = n - 1

  do j0 = 0, n0
    do i0 = 0, j0
      A(i0+1,j0+1) = real ( i4_choose ( j0, i0 ), kind = rk ) & 
                   / real ( i4_choose ( n0, i0 ), kind = rk )
    end do
  end do

  return
end
subroutine monomial_to_chebyshev ( n, mcoef, ccoef )

!*****************************************************************************80
!
!! monomial_to_chebyshev() converts from monomial to Chebyshev form.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 February 2024
!
!  Author:
!
!    Original Fortran77 version by Fred Krogh.
!    This version by John Burkardt.
!
!  Input:
!
!    integer n: the order of the polynomial.
!
!    real ( kind = rk ) mcoef(0:n): the monomial coefficients of the polynomial.
! 
!  Output:
!
!    real ( kind = rk ) ccoef(0:n): the Chebyshev coefficients of the polynomial.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) ccoef(0:n)
  integer i
  integer j
  real ( kind = rk ) mcoef(0:n)
  real ( kind = rk ) tp

  if ( n < 0 ) then
    return
  end if

  ccoef(0:n) = mcoef(0:n)

  tp = 0.5D+00 ** ( n - 1 )
  ccoef(n) = tp * ccoef(n)

  if ( n == 0 ) then
    return
  end if

  ccoef(n-1) = tp * ccoef(n-1)
  do j = n - 2, 0, -1
    tp = 2.0D+00 * tp
    ccoef(j) = tp * ccoef(j)
    ccoef(j+1) = 2.0D+00 * ccoef(j+1)
    do i = j, n - 2
      ccoef(i) = ccoef(i) + ccoef(i+2)
    end do
  end do

  return
end
subroutine monomial_to_chebyshev_matrix ( n, A )

!*****************************************************************************80
!
!! monomial_to_chebyshev_matrix() converts from monomial to Chebyshev T polynomial.
!
!  Example:
!
!    N = 11
!
!    Each column must be divided by the divisor below it.
!
!      1   .   1   .   3   .  10   .   35    .  126
!      .   1   .   3   .  10   .  35    .  126    .
!      .   .   1   .   4   .  15   .   56    .  210
!      .   .   .   1   .   5   .  21    .   84    .
!      .   .   .   .   1   .   6   .   28    .  120
!      .   .   .   .   .   1   .   7    .   36    .
!      .   .   .   .   .   .   1   .    8    .   45
!      .   .   .   .   .   .   .   1    .    9    .
!      .   .   .   .   .   .   .   .    1    .   10
!      .   .   .   .   .   .   .   .    .    1    .
!      .   .   .   .   .   .   .   .    .    .    1  
!     /1  /1  /2  /4  /8 /16 /32 /64 /128 /256 /512 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 April 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N: the order of the matrix.
!
!  Output:
!
!    real ( kind = rk ) A(N,N): the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) A(n,n)
  integer i
  integer j

  if ( n <= 0 ) then
    return
  end if

  A(1:n,1:n) = 0.0D+00

  A(1,1) = 1.0D+00

  if ( n == 1 ) then
    return
  end if

  A(2,2) = 1.0D+00

  do j = 3, n
    do i = 1, n
      if ( i == 1 ) then
        A(i,j) =                          A(i+1,j-1)   / 2.0D+00
      else if ( i == 2 ) then
        A(i,j) = ( 2.0D+00 * A(i-1,j-1) + A(i+1,j-1) ) / 2.0D+00
      else if ( i < n ) then
        A(i,j) = (           A(i-1,j-1) + A(i+1,j-1) ) / 2.0D+00
      else
        A(i,j) =             A(i-1,j-1)                / 2.0D+00
      end if
    end do
  end do

  return
end
subroutine monomial_to_gegenbauer ( n, alpha, mcoef, gcoef )

!*****************************************************************************80
!
!! monomial_to_gegenbauer() converts from monomial to Gegenbauer form.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 February 2024
!
!  Author:
!
!    John Burkardt.
!
!  Input:
!
!    integer n: the order of the polynomial.
!
!    real alpha: the Gegenbauer parameter.
!
!    real mcoef(1:n+1): the monomial coefficients of the polynomial.
! 
!  Output:
!
!    real gcoef(1:n+1): the Gegenbauer coefficients of the polynomial.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ), allocatable :: A(:,:)
  real ( kind = rk ) alpha
  real ( kind = rk ) gcoef(0:n)
  real ( kind = rk ) mcoef(0:n)

  if ( n < 0 ) then
    return
  end if
!
!  Get the matrix.
!
  allocate ( A(0:n,0:n) )
  call monomial_to_gegenbauer_matrix ( n + 1, alpha, A )
!
!  Apply the linear transformation.
!
  gcoef = matmul ( A, mcoef )

  deallocate ( A )

  return
end
subroutine monomial_to_gegenbauer_matrix ( n, alpha, A )

!*****************************************************************************80
!
!! monomial_to_gegenbauer_matrix(): monomial to Gegenbauer conversion matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 April 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer n: the order of A.
!
!    real ( kind = rk ) alpha: the parameter.
!
!  Output:
!
!    real ( kind = rk ) A[N,N]: the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) A(0:n-1,0:n-1)
  real ( kind = rk ) alpha
  real ( kind = rk ) bot
  integer i
  integer ilo
  integer j
  real ( kind = rk ) r8_factorial
  real ( kind = rk ) top

  A(0:n-1,0:n-1) = 0.0D+00

  if ( n <= 0 ) then
    return
  end if

  do j = 0, n - 1

    ilo = mod ( j, 2 )

    do i = ilo, j, 2

      top = ( i + alpha ) * r8_factorial ( j ) * gamma ( alpha )

      bot = 2.0D+00 ** j * r8_factorial ( ( j - i ) / 2 ) &
        * gamma ( ( j + i ) / 2 + alpha + 1.0D+00 )

      A(i,j) = top / bot

    end do
  end do

  return
end
subroutine monomial_to_hermite ( n, mcoef, hcoef )

!*****************************************************************************80
!
!! monomial_to_hermite() converts from monomial to Hermite form.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 February 2024
!
!  Author:
!
!    John Burkardt.
!
!  Input:
!
!    integer n: the order of the polynomial.
!
!    real mcoef(1:n+1): the monomial coefficients of the polynomial.
! 
!  Output:
!
!    real hcoef(1:n+1): the Hermite coefficients of the polynomial.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ), allocatable :: A(:,:)
  real ( kind = rk ) hcoef(0:n)
  real ( kind = rk ) mcoef(0:n)

  if ( n < 0 ) then
    return
  end if
!
!  Get the matrix.
!
  allocate ( A(0:n,0:n) )
  call monomial_to_hermite_matrix ( n + 1, A )
!
!  Apply the linear transformation.
!
  hcoef = matmul ( A, mcoef )

  deallocate ( A )

  return
end
subroutine monomial_to_hermite_matrix ( n, A )

!*****************************************************************************80
!
!! monomial_to_hermite_matrix() converts monomial to Hermite form.
!
!  Example:
!
!    N = 11 (Each column must be divided by the factor below it).
!
!    1     .     2     .    12     .   120     .  1680     . 30240
!    .     1     .     6     .    60     .   840     . 15120     .
!    .     .     1     .    12     .   180     .  3360     . 75600
!    .     .     .     1     .    20     .   420     . 10080     .
!    .     .     .     .     1     .    30     .   840     . 25200
!    .     .     .     .     .     1     .    42     .  1512     .
!    .     .     .     .     .     .     1     .    56     .  2520
!    .     .     .     .     .     .     .     1     .    72     .
!    .     .     .     .     .     .     .     .     1     .    90
!    .     .     .     .     .     .     .     .     .     1     .
!    .     .     .     .     .     .     .     .     .     .     1
!
!   /1    /2    /4    /8   /16   /32   /64  /128  /256  /512 /1024
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 February 2024
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
!    real ( kind = rk ) A(N,N), the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) A(n,n)
  integer i
  integer j

  if ( n <= 0 ) then
    return
  end if

  A(1:n,1:n) = 0.0D+00

  A(1,1) = 1.0D+00

  if ( n == 1 ) then
    return
  end if

  A(2,2) = 0.5D+00

  if ( n == 2 ) then
    return
  end if

  do j = 3, n
    do i = 1, n
      if ( i == 1 ) then
        A(i,j) = ( real ( j - 2, kind = rk ) * A(i,j-2)              ) / 2.0D+00
      else
        A(i,j) = ( real ( j - 2, kind = rk ) * A(i,j-2) + A(i-1,j-1) ) / 2.0D+00
      end if
    end do
  end do

  return
end
subroutine monomial_to_laguerre ( n, mcoef, lcoef )

!*****************************************************************************80
!
!! monomial_to_laguerre() converts from monomial to Laguerre form.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 February 2024
!
!  Author:
!
!    John Burkardt.
!
!  Input:
!
!    integer n: the order of the polynomial.
!
!    real mcoef(1:n+1): the monomial coefficients of the polynomial.
! 
!  Output:
!
!    real lcoef(1:n+1): the Laguerre coefficients of the polynomial.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ), allocatable :: A(:,:)
  real ( kind = rk ) lcoef(0:n)
  real ( kind = rk ) mcoef(0:n)

  if ( n < 0 ) then
    return
  end if
!
!  Get the matrix.
!
  allocate ( A(0:n,0:n) )
  call monomial_to_laguerre_matrix ( n + 1, A )
!
!  Apply the linear transformation.
!
  lcoef = matmul ( A, mcoef )

  deallocate ( A )

  return
end
subroutine monomial_to_laguerre_matrix ( n, A )

!*****************************************************************************80
!
!! monomial_to_laguerre_matrix() converts from monomial to Laguerre form.
!
!  Example:
!
!    N = 9
!
!        1       1       2       6      24     120     720    5040    40320
!        .      -1      -4     -18     -96    -600   -4320  -35280  -322560
!        .       .       2      18     144    1200   10800  105840  1128960
!        .       .       .      -6     -96   -1200  -14400 -176400 -2257920
!        .       .       .       .      24     600   10800  176400  2822400
!        .       .       .       .       .    -120   -4320 -105840 -2257920
!        .       .       .       .       .       .     720   35280  1128960
!        .       .       .       .       .       .       .   -5040  -322560
!        .       .       .       .       .       .       .       .    40320
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 February 2024
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
!    real ( kind = rk ) A(N,N), the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) A(n,n)
  integer i
  integer j

  if ( n <= 0 ) then
    return
  end if

  A(1:n,1:n) = 0.0D+00

  A(1,1) = 1.0D+00

  if ( n == 1 ) then
    return
  end if

  A(1,2) = 1.0D+00
  A(2,2) = -1.0D+00

  if ( n == 2 ) then
    return
  end if

  do j = 3, n
    do i = 1, n
      if ( i == 1 ) then
        A(i,j) = real ( j - 1, kind = rk ) * ( A(i,j-1)              )
      else
        A(i,j) = real ( j - 1, kind = rk ) * ( A(i,j-1) - A(i-1,j-1) )
      end if
    end do
  end do

  return
end
subroutine monomial_to_legendre ( n, mcoef, lcoef )

!*****************************************************************************80
!
!! monomial_to_legendre() converts from monomial to Legendre form.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 February 2024
!
!  Author:
!
!    John Burkardt.
!
!  Input:
!
!    integer n: the order of the polynomial.
!
!    real mcoef(1:n+1): the monomial coefficients of the polynomial.
! 
!  Output:
!
!    real lcoef(1:n+1): the Legendre coefficients of the polynomial.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ), allocatable :: A(:,:)
  real ( kind = rk ) lcoef(0:n)
  real ( kind = rk ) mcoef(0:n)

  if ( n < 0 ) then
    return
  end if
!
!  Get the matrix.
!
  allocate ( A(0:n,0:n) )
  call monomial_to_legendre_matrix ( n + 1, A )
!
!  Apply the linear transformation.
!
  lcoef = matmul ( A, mcoef )

  deallocate ( A )

  return
end
subroutine monomial_to_legendre_matrix ( n, A )

!*****************************************************************************80
!
!! monomial_to_legendre_matrix(): convert monomial to Legendre basis.
!
!  Discussion:
!
!    If PM(x) is a linear combination of monomials
!    with coefficients CM, then PL(x) is a linear combination of
!    Legendre polynomials with coefficients CL = A * CM.
!    
!    Note that we assume the coefficients are ordered such that
!    the constant term is first.
!
!  Example:
!
!    N = 11 (each column must be divided by the underlying factor).
!
!       1     .      1     .      7     .    33    .   715    . 4199
!       .     1      .     3      .    27     .  143     . 3315    .
!       .     .      2     .     20     .   110    .  2600    .16150
!       .     .      .     2      .    28     .  182     . 4760    .
!       .     .      .     .      8     .    72    .  2160    .15504
!       .     .      .     .      .     8     .   88     . 2992    .
!       .     .      .     .      .     .    16    .   832    . 7904
!       .     .      .     .      .     .     .   16     .  960    .
!       .     .      .     .      .     .     .    .   128    . 2176
!       .     .      .     .      .     .     .    .     .  128    .
!       .     .      .     .      .     .     .    .     .    .  256
!
!      /1    /1     /3    /5    /35   /63  /231 /429 /6435/12155/46189  
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 February 2024
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
!    real A(N,N), the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) A(n,n)
  integer i
  integer j

  A(1:n,1:n) = 0.0D+00

  if ( n <= 0 ) then
    return
  end if

  A(1,1) = 1.0D+00

  if ( n == 1 ) then
    return
  end if

  A(2,2) = 1.0D+00

  if ( n == 2 ) then
    return
  end if

  do j = 3, n
    do i = 1, n
      if ( i == 1 ) then

        A(i,j) = (     i     ) * A(i+1,j-1) / ( 2 * i + 1 )

      else if ( i < n ) then

        A(i,j) = (     i - 1 ) * A(i-1,j-1) / ( 2 * i - 3 ) &
               + (     i     ) * A(i+1,j-1) / ( 2 * i + 1 )

      else

        A(i,j) = (     i - 1 ) * A(i-1,j-1) / ( 2 * i - 3 )

      end if
    end do
  end do

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
!  Input:
!
!    integer N: the argument of the factorial function.
!    If N is less than 1, the function value is returned as 1.
!
!  Output:
!
!    real ( kind = rk ) R8_FACTORIAL: the factorial of N.
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

