function nchoosek ( n, k )

!*****************************************************************************80
!
!! nchoosek() computes the binomial coefficient C(N,K) as an I4.
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
!    integer nchoosek, the number of combinations of N
!    things taken K at a time.
!
  implicit none

  integer i

  integer k
  integer mn
  integer mx
  integer n
  integer nchoosek
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

  nchoosek = value

  return
end
subroutine sigmoid_derivative_coef ( n, coef )

!*****************************************************************************80
!
!! sigmoid_derivative_coef(): sigmoid derivative expansion coefficients.
!
!  Discussion:
!
!    s(x) = 1 / ( 1 + exp ( - x ) )
!
!    s(x)^(n) = sum ( 1 <= j <= n + 1 ) coef(j) * s(x)^j 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 May 2024
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joe McKenna,
!    Derivatives of the sigmoid function,
!    https://joepatmckenna.github.io/calculus/derivative/sigmoid!20function/linear!20albegra/2018/01/20/sigmoid-derivs/
!
!  Input:
!
!    integer n: the order of the derivative.
!    0 <= n.
!
!  Output:
!
!    real ( kind = rk8 ) coef(n+2): the expansion coefficients.
!    coef(1) is always 0, and coef(2) is always 1.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  real ( kind = rk8 ) cnk
  real ( kind = rk8 ) coef(n+2)
  integer j
  integer k
  real ( kind = rk8 ) mop
  integer nchoosek
 
  coef(1:n+2) = 0.0D+00

  do k = 0, n
    cnk = 0.0D+00
    mop = -1
    do j = 0, k
      mop = - mop
      cnk = cnk + mop * ( j + 1 )**n * nchoosek ( k, j )
    end do
    coef(k+2) = cnk
  end do

  return
end
function sigmoid_derivative_value ( n, x )

!*****************************************************************************80
!
!! sigmoid_derivative_value(): evaluate sigmoid derivative at x.
!
!  Discussion:
!
!    s(x) = 1 / ( 1 + exp ( - x ) )
!
!    s(x)^(n) = sum ( 1 <= j <= n + 1 ) coef(j) * s(x)^j 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 May 2024
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joe McKenna,
!    Derivatives of the sigmoid function,
!    https://joepatmckenna.github.io/calculus/derivative/sigmoid!20function/linear!20albegra/2018/01/20/sigmoid-derivs/
!
!  Input:
!
!    integer n: the order of the derivative.
!
!    real ( kind = rk8 ) x: the evaluation point.
!
!  Output:
!
!    real ( kind = rk8 ) sigmoid_derivative_value: the value of s^(n)(x).
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) cnk
  integer j
  integer k
  integer mop
  integer n
  integer nchoosek
  real ( kind = rk8 ) sig
  real ( kind = rk8 ) sigk
  real ( kind = rk8 ) sigmoid_derivative_value
  real ( kind = rk8 ) value
  real ( kind = rk8 ) x

  sig = 1.0D+00 / ( 1.0D+00 + exp ( - x ) )
  value = 0.0D+00
  sigk = 1.0D+00

  do k = 0, n
    sigk = sigk * sig
    cnk = 0.0D+00
    mop = -1
    do j = 0, k
      mop = - mop
      cnk = cnk + mop * ( j + 1 )**n * nchoosek ( k, j )
    end do
    value = value + cnk * sigk
  end do

  sigmoid_derivative_value = value

  return 
end

