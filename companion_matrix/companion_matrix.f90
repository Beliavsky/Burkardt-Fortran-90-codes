subroutine companion_chebyshev ( n, p, C )

!*****************************************************************************80
!
!! companion_chebyshev() returns the Chebyshev basis companion matrix for a polynomial.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 April 2024
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    John Boyd,
!    Solving Transcendental Equations,
!    The Chebyshev Polynomial Proxy and other Numerical Rootfinders,
!    Perturbation Series, and Oracles,
!    SIAM, 2014,
!    ISBN: 978-1-611973-51-8,
!    LC: QA:353.T7B69
!
!  Input:
!
!    real p(n+1): the polynomial coefficients, in order of increasing degree.
!
!  Output:
!
!    real C(n,n): the companion matrix.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  real ( kind = rk8 ) C(n,n)
  integer i
  integer j
  real ( kind = rk8 ) p(n+1)

  C(1:n,1:n) = 0.0

  C(1,2) = 1.0

  do i = 2, n - 1
    C(i,i-1) = 0.5
    C(i,i+1) = 0.5
  end do

  do j = 1, n
    C(n,j) = - 0.5 * p(j) / p(n+1)
  end do

  C(n,n-1) = C(n,n-1) + 0.5

  return
end
subroutine companion_gegenbauer ( n, p, C, alpha )

!*****************************************************************************80
!
!! companion_gegenbauer() returns the Gegenbauer basis companion matrix for a polynomial.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 April 2024
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    John Boyd,
!    Solving Transcendental Equations,
!    The gegenbauer Polynomial Proxy and other Numerical Rootfinders,
!    Perturbation Series, and Oracles,
!    SIAM, 2014,
!    ISBN: 978-1-611973-51-8,
!    LC: QA:353.T7B69
!
!  Input:
!
!    real p(n+1): the polynomial coefficients, in order of increasing degree.
!
!    real alpha: the Gegenbauer parameter.
!
!  Output:
!
!    real C(n,n): the companion matrix.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  real ( kind = rk8 ) C(n,n)
  real ( kind = rk8 ) alpha
  integer i
  integer j
  real ( kind = rk8 ) p(n+1)

  C(1:n,1:n) = 0.0

  C(1,2) = 0.5 / alpha

  do i = 2, n - 1
    C(i,i-1) = 0.5 * ( i - 1 + 2 * alpha - 1 ) / ( i - 1 + alpha )
    C(i,i+1) = 0.5 * ( i                     ) / ( i - 1 + alpha )
  end do

  do j = 1, n
    C(n,j) = - 0.5 * p(j) * n / ( n - 1 + alpha ) / p(n+1)
  end do

  C(n,n-1) = C(n,n-1) + 0.5 * ( n - 2 + 2 * alpha ) / ( n - 1 + alpha )

  return
end
subroutine companion_hermite ( n, p, C )

!*****************************************************************************80
!
!! companion_hermite() returns the Hermite basis companion matrix for a polynomial.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 April 2024
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    John Boyd,
!    Solving Transcendental Equations,
!    The legendre Polynomial Proxy and other Numerical Rootfinders,
!    Perturbation Series, and Oracles,
!    SIAM, 2014,
!    ISBN: 978-1-611973-51-8,
!    LC: QA:353.T7B69
!
!  Input:
!
!    real p(n+1): the polynomial coefficients, in order of increasing degree.
!
!  Output:
!
!    real C(n,n): the companion matrix.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  real ( kind = rk8 ) C(n,n)
  integer i
  integer j
  real ( kind = rk8 ) p(n+1)

  C(1:n,1:n) = 0.0

  do i = 1, n - 1
    C(i,i+1) = 0.5
  end do

  do i = 2, n
    C(i,i-1) = i - 1
  end do

  do j = 1, n
    C(n,j) = C(n,j) - p(j) / 2.0 * p(n+1)
  end do

  return
end
subroutine companion_laguerre ( n, p, C )

!*****************************************************************************80
!
!! companion_laguerre() returns the Laguerre basis companion matrix for a polynomial.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 April 2024
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    John Boyd,
!    Solving Transcendental Equations,
!    The gegenbauer Polynomial Proxy and other Numerical Rootfinders,
!    Perturbation Series, and Oracles,
!    SIAM, 2014,
!    ISBN: 978-1-611973-51-8,
!    LC: QA:353.T7B69
!
!  Input:
!
!    real p(n+1): the polynomial coefficients, in order of increasing degree.
!
!  Output:
!
!    real C(n,n): the companion matrix.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  real ( kind = rk8 ) C(n,n)
  integer i
  integer j
  real ( kind = rk8 ) p(n+1)

  C(1:n,1:n) = 0.0

  do i = 1, n
    C(i,i) = 2 * i - 1
  end do

  do i = 2, n
    C(i,i-1) = - i + 1
  end do

  do i = 1, n - 1
    C(i,i+1) = - i
  end do

  do j = 1, n
    C(n,j) = C(n,j) + n * p(j) / p(n+1)
  end do

  return
end
subroutine companion_legendre ( n, p, C )

!*****************************************************************************80
!
!! companion_legendre() returns the Legendre basis companion matrix for a polynomial.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 April 2024
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    John Boyd,
!    Solving Transcendental Equations,
!    The legendre Polynomial Proxy and other Numerical Rootfinders,
!    Perturbation Series, and Oracles,
!    SIAM, 2014,
!    ISBN: 978-1-611973-51-8,
!    LC: QA:353.T7B69
!
!  Input:
!
!    real p(n+1): the polynomial coefficients, in order of increasing degree.
!
!  Output:
!
!    real C(n,n): the companion matrix.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  real ( kind = rk8 ) C(n,n)
  integer i
  integer j
  real ( kind = rk8 ) p(n+1)

  C(1:n,1:n) = 0.0

  C(1,2) = 1.0

  do i = 2, n - 1
    C(i,i-1) = real ( i - 1, kind = rk8 ) / real ( 2 * i - 1, kind = rk8 )
    C(i,i+1) = real ( i,     kind = rk8 ) / real ( 2 * i - 1, kind = rk8 )
  end do

  do j = 1, n
    C(n,j) = - p(j) / p(n+1) * real ( n, kind = rk8 ) &
      / real ( 2 * n - 1, kind = rk8 )
  end do

  C(n,n-1) = C(n,n-1) + real ( n - 1, kind = rk8 ) &
    / real ( 2 * n - 1, kind = rk8 )

  return
end
subroutine companion_monomial ( n, p, C )

!*****************************************************************************80
!
!! companion_monomial() returns the monomial basis companion matrix for a polynomial.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 April 2024
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    John Boyd,
!    Solving Transcendental Equations,
!    The Chebyshev Polynomial Proxy and other Numerical Rootfinders,
!    Perturbation Series, and Oracles,
!    SIAM, 2014,
!    ISBN: 978-1-611973-51-8,
!    LC: QA:353.T7B69
!
!  Input:
!
!    real p(n+1): the polynomial coefficients, in order of increasing degree.
!
!  Output:
!
!    real C(n,n): the companion matrix.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  real ( kind = rk8 ) C(n,n)
  integer i
  integer j
  real ( kind = rk8 ) p(n+1)

  C(1:n,1:n) = 0.0

  do i = 1, n - 1
    C(i,i+1) = 1.0
  end do

  do j = 1, n
    C(n,j) = - p(j) / p(n+1)
  end do

  return
end

