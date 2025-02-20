program main

!*****************************************************************************80
!
!! companion_matrix_test() tests companion_matrix().
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
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'companion_matrix_test():'
  write ( *, '(a)' ) '  Fortran90 version '
  write ( *, '(a)' ) '  companion_matrix() computes the companion matrix'
  write ( *, '(a)' ) '  of a polynomial, in various bases.'

  call companion_chebyshev_test ( )
  call companion_gegenbauer_test ( )
  call companion_hermite_test ( )
  call companion_laguerre_test ( )
  call companion_legendre_test ( )
  call companion_monomial_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'companion_matrix_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop
end
subroutine companion_chebyshev_test ( )

!*****************************************************************************80
!
!! companion_chebyshev_test() tests companion_chebyshev().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 June 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ck8 = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk8 ) A(n+1,n+1)
  real ( kind = rk8 ) C(n,n)
  real ( kind = rk8 ) p(0:n)
  real ( kind = rk8 ) q(0:n)
  complex ( kind = ck8 ) r(n)
  complex ( kind = ck8 ) s(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'companion_chebyshev_test():'
  write ( *, '(a)' ) '  companion_chebyshev() computes the companion matrix'
  write ( *, '(a)' ) '  of a polynomial p(x) in the Chebyshev basis.'

  p = (/ 6.0, 5.0, 4.0, 3.0, 2.0, 1.0 /)

  call polynomial_chebyshev_print ( n, p, '  p(x)' )
!
!  Compute monomial form so we can get roots directly.
!
  call chebyshev_to_monomial_matrix ( n + 1, A )
  q = matmul ( A, p )
  call polynomial_monomial_print ( n, q, '  Monomial q(x)' )

  call r8poly_roots ( n, q, r )
  call c8vec_print ( n, r, '  Roots of q(x):' )

  call companion_chebyshev ( n, p, C )
  call r8mat_print ( n, n, C, '  Chebyshev companion matrix C(p):' )

  call eigs ( n, C, s )
  call c8vec_print ( n, s, '  Eigenvalues of C(p):' )

  return
end
subroutine companion_gegenbauer_test ( )

!*****************************************************************************80
!
!! companion_gegenbauer_test() tests companion_gegenbauer().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 June 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ck8 = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk8 ) A(n+1,n+1)
  real ( kind = rk8 ) alpha
  real ( kind = rk8 ) C(n,n)
  integer i
  real ( kind = rk8 ) p(0:n)
  real ( kind = rk8 ) q(0:n)
  complex ( kind = ck8 ) r(n)
  complex ( kind = ck8 ) s(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'companion_gegenbauer_test():'
  write ( *, '(a)' ) '  companion_gegenbauer() computes the companion matrix'
  write ( *, '(a)' ) '  of a polynomial p(x) in the Gegenbauer basis.'

  p = (/ 6.0, 5.0, 4.0, 3.0, 2.0, 1.0 /)

  call polynomial_gegenbauer_print ( n, p, '  p(x)' )
!
!  Try several parameter values.
!
  do i = 1, 2

    alpha = i / 2.0D+00

    write ( *, '(a)' ) ''
    write ( *, '(a,g14.6)' ) '  Gegenbauer parameter alpha will be ', alpha
!
!  Compute monomial form so we can get roots directly.
!
    call gegenbauer_to_monomial_matrix ( n + 1, A, alpha )
    q = matmul ( A, p )
    call polynomial_monomial_print ( n, q, '  Monomial q(x)' )

    call r8poly_roots ( n, q, r )
    call c8vec_print ( n, r, '  Roots of q(x):' )

    call companion_gegenbauer ( n, p, C, alpha )
    call r8mat_print ( n, n, C, '  Gegenbauer companion matrix C(p)' )

    call eigs ( n, C, s )
    call c8vec_print ( n, s, '  Eigenvalues of C(p):' )

  end do

  return
end
subroutine companion_hermite_test ( )

!*****************************************************************************80
!
!! companion_hermite_test() tests companion_hermite().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 June 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ck8 = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk8 ) A(n+1,n+1)
  real ( kind = rk8 ) C(n,n)
  real ( kind = rk8 ) p(0:n)
  real ( kind = rk8 ) q(0:n)
  complex ( kind = ck8 ) r(n)
  complex ( kind = ck8 ) s(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'companion_hermite_test():'
  write ( *, '(a)' ) '  companion_hermite() computes the companion matrix'
  write ( *, '(a)' ) '  of a polynomial p(x) in the Hermite basis.'

  p = (/ 6.0, 5.0, 4.0, 3.0, 2.0, 1.0 /)

  call polynomial_hermite_print ( n, p, '  Hermite p(x)' )
!
!  Compute monomial form so we can get roots directly.
!
  call h_to_monomial_matrix ( n + 1, A )
  q = matmul ( A, p )
  call polynomial_monomial_print ( n, q, '  Monomial q(x)' )

  call r8poly_roots ( n, q, r )
  call c8vec_print ( n, r, '  Roots of q(x):' )
 
  call companion_hermite ( n, p, C )
  call r8mat_print ( n, n, C, '  Hermite companion matrix C(p):' )

  call eigs ( n, C, s )
  call c8vec_print ( n, s, '  Eigenvalues of C(p):' )

  return
end
subroutine companion_laguerre_test ( )

!*****************************************************************************80
!
!! companion_laguerre_test() tests companion_laguerre().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 June 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ck8 = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk8 ) A(n+1,n+1)
  real ( kind = rk8 ) C(n,n)
  real ( kind = rk8 ) p(0:n)
  real ( kind = rk8 ) q(0:n)
  complex ( kind = ck8 ) r(n)
  complex ( kind = ck8 ) s(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'companion_laguerre_test():'
  write ( *, '(a)' ) '  companion_laguerre() computes the companion matrix'
  write ( *, '(a)' ) '  of a polynomial p(x) in the Laguerre basis.'

  p = (/ 6.0, 5.0, 4.0, 3.0, 2.0, 1.0 /)

  call polynomial_laguerre_print ( n, p, '  Laguerre p(x)' )
!
!  Compute monomial form so we can get roots directly.
!
  call laguerre_to_monomial_matrix ( n + 1, A )
  q = matmul ( A, p )
  call polynomial_monomial_print ( n, q, '  Monomial q(x)' )

  call r8poly_roots ( n, q, r )
  call c8vec_print ( n, r, '  Roots of q(x):' )
 
  call companion_laguerre ( n, p, C )
  call r8mat_print ( n, n, C, '  Laguerre companion matrix C(p):' )

  call eigs ( n, C, s )
  call c8vec_print ( n, s, '  Eigenvalues of C(p):' )

  return
end
subroutine companion_legendre_test ( )

!*****************************************************************************80
!
!! companion_legendre_test() tests companion_legendre().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 June 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ck8 = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk8 ) A(n+1,n+1)
  real ( kind = rk8 ) C(n,n)
  real ( kind = rk8 ) p(0:n)
  real ( kind = rk8 ) q(0:n)
  complex ( kind = ck8 ) r(n)
  complex ( kind = ck8 ) s(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'companion_legendre_test():'
  write ( *, '(a)' ) '  companion_legendre() computes the companion matrix'
  write ( *, '(a)' ) '  of a polynomial p(x) in the Legendre basis.'

  p = (/ 6.0, 5.0, 4.0, 3.0, 2.0, 1.0 /)

  call polynomial_legendre_print ( n, p, '  Legendre p(x)' )
!
!  Compute monomial form so we can get roots directly.
!
  call legendre_to_monomial_matrix ( n + 1, A )
  q = matmul ( A, p )
  call polynomial_monomial_print ( n, q, '  Monomial q(x)' )

  call r8poly_roots ( n, q, r )
  call c8vec_print ( n, r, '  Roots of q(x):' )
 
  call companion_legendre ( n, p, C )
  call r8mat_print ( n, n, C, 'Legendre companion matrix C(p):' )

  call eigs ( n, C, s )
  call c8vec_print ( n, s, '  Eigenvalues of C(p):' )

  return
end
subroutine companion_monomial_test ( )

!*****************************************************************************80
!
!! companion_monomial_test() tests companion_monomial().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 June 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ck8 = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk8 ) C(n,n)
  real ( kind = rk8 ) p(0:n)
  complex ( kind = ck8 ) r(n)
  complex ( kind = ck8 ) s(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'companion_monomial_test():'
  write ( *, '(a)' ) '  companion_monomial() computes the companion matrix'
  write ( *, '(a)' ) '  of a polynomial p(x) in the monomial basis.'

  p = (/ 6.0, 5.0, 4.0, 3.0, 2.0, 1.0 /)

  call polynomial_monomial_print ( n, p, '  p(x)' )
!
!  Get roots directly.
!
  call r8poly_roots ( n, p, r )
  call c8vec_print ( n, r, '  Roots of q(x):' )

  call companion_monomial ( n, p, C )
  call r8mat_print ( n, n, C, '  Monomial companion matrix C(p):' )

  call eigs ( n, C, s )
  call c8vec_print ( n, s, '  Eigenvalues of C(p):' )

  return
end
subroutine chebyshev_to_monomial_matrix ( n, A )

!*****************************************************************************80
!
!! chebyshev_to_monomial_matrix(): convert Chebyshev polynomial to monomial form.
!
!  Discussion:
!
!    1     0     -1      0       1     0     -1    0
!    0     1      0     -3       0     5      0   -7
!    0     0      2      0      -8     0     18    0
!    0     0      0      4       0   -20      0   56
!    0     0      0      0       8     0    -48    0
!    0     0      0      0       0    16      0 -112
!    0     0      0      0       0     0     32    0
!    0     0      0      0       0     0      0   64
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
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical subroutines,
!    US Department of Commerce, 1964.
!
!  Input:
!
!    integer N, the order of the matrix.
!
!  Output:
!
!    real A(n,n), the matrix.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  real ( kind = rk8 ) A(n,n)
  integer i

  A(1:n,1:n) = 0.0

  A(1,1) = 1.0

  if ( n == 1 ) then
    return
  end if

  A(2,2) = 1.0
 
  do i = 2, n - 1
    A(1,i+1)     =                  - A(1,i-1)
    A(2:i-1,i+1) = 2.0 * A(1:i-2,i) - A(2:i-1,i-1)
    A(i,i+1) =     2.0 * A(i-1,i)
    A(i+1,i+1) =   2.0 * A(i,i)
  end do
 
  return
end
subroutine gegenbauer_to_monomial_matrix ( n, A, alpha )

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
!    04 April 2024
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
!    real A(N,N): the matrix.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  real ( kind = rk8 ) A(n,n)
  real ( kind = rk8 ) alpha
  real ( kind = rk8 ) c1
  real ( kind = rk8 ) c2
  integer j
  integer nn

  A(1:n,1:n) = 0.0

  A(1,1) = 1.0

  if ( n == 1 ) then
    return
  end if

  A(2,2) = 2.0 * alpha
!
!  Perform convex sum.
!  Translating "(n+1) C(n+1) = 2 (n+alpha) x C(n) - ( n + 2 alpha - 1 ) C(n-1)"
!  drove me nuts, between indexing at 1 rather than 0, and dealing with
!  the interpretation of "n+1", because I now face the rare "off by 2" error!
!
  do j = 3, n
    nn = j - 2 
    c1 = ( 2 * nn + 2 * alpha     ) / ( nn + 1 )
    c2 = (   - nn - 2 * alpha + 1 ) / ( nn + 1 )
    A(2:j,j)   =              c1 * A(1:j-1,j-1)
    A(1:j-2,j) = A(1:j-2,j) + c2 * A(1:j-2,j-2)
  end do

  return
end
subroutine h_to_monomial_matrix ( n, A )

!*****************************************************************************80
!
!! h_to_monomial_matrix(): physicist's Hermite to monomial conversion matrix.
!
!  Example:
!
!    N = 11
!
!      1  .  -2   .    12     .    -120     .   1680     .   -30240
!      .  2   .  12     .   120      .  -1680      . 30240        .
!      .  .   4   .   -48     .     720     . -13440     .   302400
!      .  .   .   8     .  -160       .  3360   .   -80640        .
!      .  .   .   .    16     .    -480     .  13440     .  -403200
!      .  .   .   .     .    32       . -1344      . 48384        .
!      .  .   .   .     .     .      64     .  -3584     .   161280
!      .  .   .   .     .     .       .   128      . -9216        .
!      .  .   .   .     .     .       .     .    256     .   -23040
!      .  .   .   .     .     .       .     .      .   512        .
!      .  .   .   .     .     .       .     .      .     .     1024
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
!  Input:
!
!    integer N, the order of A.
!
!  Output:
!
!    real A(N,N), the matrix.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  real ( kind = rk8 ) A(n,n)
  integer i
  integer j

  A(1:n,1:n) = 0.0

  A(1,1) = 1.0

  if ( n == 1 ) then
    return
  end if

  A(2,2) = 2.0

  if ( n == 2 ) then
    return
  end if

  do j = 3, n
    do i = 1, n
      if ( i == 1 ) then
        A(i,j) =                  - 2.0 * ( j - 2 ) * A(i,j-2)
      else
        A(i,j) = 2.0 * A(i-1,j-1) - 2.0 * ( j - 2 ) * A(i,j-2)
      end if
    end do
  end do

  return
end
subroutine laguerre_to_monomial_matrix ( n, A )

!*****************************************************************************80
!
!! laguerre_to_monomial_matrix() converts from Laguerre to monomial form.
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
!    01 April 2024
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical subroutines,
!    US Department of Commerce, 1964.
!
!  Input:
!
!    integer N: the order of the matrix.
!
!  Output:
!
!    real A(N,N): the matrix.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  real ( kind = rk8 ) A(n,n)
  integer i
  integer j

  A(1:n,1:n) = 0.0

  A(1,1) = 1.0

  if ( n == 1 ) then
    return
  end if

  A(1,2) = 1.0
  A(2,2) = -1.0

  do j = 3, n
    do i = 1, n
      if ( i == 1 ) then
        A(i,j) = (  ( 2 * j - 3 ) * A(i,j-1) &
                 +  (   - j + 2 ) * A(i,j-2) ) &
                 /  (     j - 1 )
      else
        A(i,j) = (  ( 2 * j - 3 ) * A(i,j-1) &
                 -  (         1 ) * A(i-1,j-1) &
                 +  (   - j + 2 ) * A(i,j-2) ) &
                 /  (     j - 1 )
      end if
    end do
  end do

  return
end
subroutine legendre_to_monomial_matrix ( n, A )

!*****************************************************************************80
!
!! legendre_to_monomial_matrix(): Legendre coefficient conversion matrix.
!
!  Discussion:
!
!    If PL(x) is a linear combination of Legendre polynomials
!    with coefficients CL, then PM(x) is a linear combination of
!    monomials with coefficients CM = A * CL.
!    
!    The coefficients are ordered so the constant term is first.
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
!  Input:
!
!    integer N, the order of A.
!
!  Output:
!
!    real A(N,N), the matrix.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  real ( kind = rk8 ) A(n,n)
  integer i
  integer j

  A(1:n,1:n) = 0.0

  A(1,1) = 1.0

  if ( n == 1 ) then
    return
  end if

  A(2,2) = 1.0

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
subroutine polynomial_chebyshev_print ( d, c, label )

!*****************************************************************************80
!
!! polynomial_chebyshev_print() prints a polynomial in the Chebyshev basis.
!
!  Discussion:
!
!    The form of a Chebyshev polynomial is:
!
!      p(x) = c(0)*T0(x) + c(1)*T(1)(x) + ... + c(n)*T(n)(x)
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
!  Input:
!
!    integer D: the polynomial degree.
!
!    real C(0:d): the polynomial.
!
!    character LABEL: an optional title.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer d

  real ( kind = rk8 ) c(0:d)
  integer i
  character ( len = * ) label

  write ( *, '(a)' ) ''
  if ( 0 < len_trim ( label ) ) then
    write ( *, '(a)' ) trim ( label ), ' = '
  end if

  if ( d < 0 ) then
    write ( *, '(a)' ) '  Zero polynomial'
    return
  end if

  if ( all ( c == 0.0 ) ) then
    write ( *, '(a)' ) '  0'
    return
  end if

  do i = d, 0, -1

    if ( c(i) /= 0.0 ) then
      write ( *, '(a,g14.6,a,i2,a)' ) '  +', c(i), ' * T', i, '(x)'
    end if

  end do

  return
end
subroutine polynomial_gegenbauer_print ( d, c, label )

!*****************************************************************************80
!
!! polynomial_gegenbauer_print() prints a polynomial in the Gegenbauer basis.
!
!  Discussion:
!
!    The form of a Gegenbauer polynomial is:
!
!      p(x) = c(0)*C0(x) + c(1)*C1(x) + ... + c(n)*C(n)(x)
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
!  Input:
!
!    integer D: the polynomial degree.
!
!    real C(1:D+1): the polynomial.
!
!    character LABEL: an optional title.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer d

  real ( kind = rk8 ) c(0:d)
  integer i
  character ( len = * ) label

  write ( *, '(a)' ) ''
  if ( 0 < len_trim ( label ) ) then
    write ( *, '(a)' ) trim ( label ), ' = '
  end if

  if ( d < 0 ) then
    write ( *, '(a)' ) '  Zero polynomial'
    return
  end if

  if ( all ( c == 0.0 ) ) then
    write ( *, '(a)' ) '  0'
    return
  end if

  do i = d, 0, -1

    if ( c(i) /= 0.0 ) then
      write ( *, '(a,g14.6,a,i2,a)' ) '  +', c(i), ' * C', i, '(x)'
    end if

  end do

  return
end
subroutine polynomial_hermite_print ( d, c, label )

!*****************************************************************************80
!
!! polynomial_hermite_print() prints a polynomial in the Hermite basis.
!
!  Discussion:
!
!    The form of a Hermite polynomial is:
!
!      p(x) = c(0)*H0(x) + c(1)*H(1)(x) + ... + c(n)*H(n)(x)
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
!  Input:
!
!    integer D: the polynomial degree.
!
!    real C(1:D+1): the polynomial.
!
!    character LABEL: an optional title.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer d

  real ( kind = rk8 ) c(0:d)
  integer i
  character ( len = * ) label

  write ( *, '(a)' ) ''
  if ( 0 < len_trim ( label ) ) then
    write ( *, '(a)' ) trim ( label ), ' = '
  end if

  if ( d < 0 ) then
    write ( *, '(a)' ) '  Zero polynomial'
    return
  end if

  if ( all ( c == 0.0 ) ) then
    write ( *, '(a)' ) '  0'
    return
  end if

  do i = d, 0, -1

    if ( c(i) /= 0.0 ) then
      write ( *, '(a,g14.6,a,i2,a)' ) '  +', c(i), ' * H', i, '(x)'
    end if

  end do

  return
end
subroutine polynomial_laguerre_print ( d, c, label )

!*****************************************************************************80
!
!! polynomial_laguerre_print() prints a polynomial in the Laguerre basis.
!
!  Discussion:
!
!    The form of a Laguerre polynomial is:
!
!      p(x) = c(0)*L0(x) + c(1)*L(1)(x) + ... + c(n)*L(n)(x)
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
!  Input:
!
!    integer D: the polynomial degree.
!
!    real C(1:D+1): the polynomial coefficients.
!
!    character LABEL: an optional title.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer d

  real ( kind = rk8 ) c(0:d)
  integer i
  character ( len = * ) label

  write ( *, '(a)' ) ''
  if ( 0 < len_trim ( label ) ) then
    write ( *, '(a)' ) trim ( label ), ' = '
  end if

  if ( d < 0 ) then
    write ( *, '(a)' ) '  Zero polynomial'
    return
  end if

  if ( all ( c == 0.0 ) ) then
    write ( *, '(a)' ) '  0'
    return
  end if

  do i = d, 0, -1

    if ( c(i) /= 0.0 ) then
      write ( *, '(a,g14.6,a,i2,a)' ) '  +', c(i), ' * L', i, '(x)'
    end if

  end do

  return
end
subroutine polynomial_legendre_print ( d, c, label )

!*****************************************************************************80
!
!! polynomial_legendre_print() prints a polynomial in the Legendre basis.
!
!  Discussion:
!
!    The form of a Legendre polynomial is:
!
!      p(x) = c(0)*P0(x) + c(1)*P(1)(x) + ... + c(n)*P(n)(x)
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
!  Input:
!
!    integer D: the polynomial degree.
!
!    real C(1:D+1): the polynomial coefficients.
!
!    character LABEL: an optional title.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer d

  real ( kind = rk8 ) c(0:d)
  integer i
  character ( len = * ) label

  write ( *, '(a)' ) ''
  if ( 0 < len_trim ( label ) ) then
    write ( *, '(a)' ) trim ( label ), ' = '
  end if

  if ( d < 0 ) then
    write ( *, '(a)' ) '  Zero polynomial'
    return
  end if

  if ( all ( c == 0.0 ) ) then
    write ( *, '(a)' ) '  0'
    return
  end if

  do i = d, 0, -1

    if ( c(i) /= 0.0 ) then
      write ( *, '(a,g14.6,a,i2,a)' ) '  +', c(i), ' * P', i, '(x)'
    end if

  end do

  return
end
subroutine polynomial_monomial_print ( d, p, label )

!*****************************************************************************80
!
!! polynomial_monomial_print() prints a polynomial.
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
!  Input:
!
!    integer d: the degree of the polynomial.
!
!    real p(0:d), the polynomial coefficients.
!    p(0) is the constant term and
!    p(N) is the coefficient of X^N.
!
!    character ( len = * ) label: a title.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer d

  integer i
  real ( kind = rk8 ) mag
  real ( kind = rk8 ) p(0:d)
  character plus_minus
  character ( len = * ) label

  if ( 0 < len_trim ( label ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( label )
  end if

  write ( *, '(a)' ) ' '

  if ( d < 0 ) then
    write ( *, '( ''  p(x) = 0'' )' )
    return
  end if

  if ( p(d) < 0.0D+00 ) then
    plus_minus = '-'
  else
    plus_minus = ' '
  end if

  mag = abs ( p(d) )

  if ( 2 <= d ) then
    write ( *, '( ''  p(x) = '', a1, g14.6, '' * x ^ '', i3 )' ) &
      plus_minus, mag, d
  else if ( d == 1 ) then
    write ( *, '( ''  p(x) = '', a1, g14.6, '' * x'' )' ) &
      plus_minus, mag
  else if ( d == 0 ) then
    write ( *, '( ''  p(x) = '', a1, g14.6 )' ) plus_minus, mag
  end if

  do i = d - 1, 0, -1

    if ( p(i) < 0.0D+00 ) then
      plus_minus = '-'
    else
      plus_minus = '+'
    end if

    mag = abs ( p(i) )

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
subroutine c8vec_print ( n, a, title )

!*****************************************************************************80
!
!! c8vec_print() prints a C8VEC, with an optional title.
!
!  Discussion:
!
!    A C8VEC is a vector of C8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of components of the vector.
!
!    complex ( kind = ck8 ) A(N), the vector to be printed.
!
!    character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: ck8 = kind ( ( 1.0D+00, 1.0D+00 ) )

  integer n

  complex ( kind = ck8 ) a(n)
  integer i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,2x,2g14.6)' ) i, a(i)
  end do

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! r8mat_print() prints a real matrix.
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
!  Input:
!
!    integer M, the number of rows in A.
!
!    integer N, the number of columns in A.
!
!    real ( kind = rk8 ) A(M,N), the matrix.
!
!    character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk8 ) a(m,n)
  character ( len = * ) title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! r8mat_print_some prints some of a real matrix.
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
!  Input:
!
!    integer M, N, the number of rows and columns.
!
!    real ( kind = rk8 ) A(M,N), an M by N matrix to be printed.
!
!    integer ILO, JLO, the first row and column to print.
!
!    integer IHI, JHI, the last row and column to print.
!
!    character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: incx = 5
  integer m
  integer n

  real ( kind = rk8 ) a(m,n)
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

        write ( ctemp(j2), '(g14.6)' ) a(i,j)

      end do

      write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine r8poly_roots ( n, p, r )

!*****************************************************************************80
!
!! r8poly_roots() returns the roots of a polynomial.
!
!  Discussion:
!
!    The monomial form of a polynomial is:
!
!      p(x) = c(0)*x^0 + c(1)*x + ... + c(n)*x^(n)
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    08 June 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer n: the polynomial degree.
!
!    real p(n+1): the polynomial coefficients, constant first.
!
!  Output:
!
!    complex r(n): the polynomial roots.
!
  implicit none

  integer, parameter :: ck8 = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk8 = kind ( 1.0D+00 )

  logical fail
  integer i
  integer n
  real ( kind = rk8 ) op(n+1)
  real ( kind = rk8 ) p(n+1)
  complex ( kind = ck8 ) r(n)
  real ( kind = rk8 ) zeroi(n)
  real ( kind = rk8 ) zeror(n)
      
  op(1:n+1) = p(n+1:1:-1)

  call rpoly ( op, n, zeror, zeroi, fail )

  if ( fail ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'r8poly_roots(): Fatal error!'
    write ( *, '(a)' ) '  rpoly() from toms493() failed.'
    stop 1
  end if

  do i = 1, n
    r(i) = dcmplx ( zeror(i), zeroi(i) )
  end do
      
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
!    15 August 2021
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

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
