program polynomial_conversion_test

!*****************************************************************************80
!
!! polynomial_conversion_test() tests polynomial_conversion().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 April 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'polynomial_conversion_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test polynomial_conversion().'

  call bernstein_to_legendre01_test ( )
  call bernstein_to_legendre01_matrix_test ( )
  call legendre01_to_bernstein_test ( )
  call legendre01_to_bernstein_matrix_test ( )
  call bernstein_legendre01_bernstein_test ( )

  call bernstein_to_monomial_test ( )
  call bernstein_to_monomial_matrix_test ( )
  call monomial_to_bernstein_test ( )
  call monomial_to_bernstein_matrix_test ( )
  call bernstein_monomial_bernstein_test ( )

  call chebyshev_to_monomial_test ( )
  call monomial_to_chebyshev_test ( )
  call chebyshev_monomial_chebyshev_test ( )

  call gegenbauer_to_monomial_test ( )
  call gegenbauer_to_monomial_matrix_test ( )
  call monomial_to_gegenbauer_test ( )
  call monomial_to_gegenbauer_matrix_test ( )
  call gegenbauer_monomial_gegenbauer_test ( )

  call hermite_to_monomial_test ( )
  call hermite_to_monomial_matrix_test ( )
  call monomial_to_hermite_test ( )
  call monomial_to_hermite_matrix_test ( )
  call hermite_monomial_hermite_test ( )

  call laguerre_to_monomial_test ( )
  call laguerre_to_monomial_matrix_test ( )
  call monomial_to_laguerre_test ( )
  call monomial_to_laguerre_matrix_test ( )
  call laguerre_monomial_laguerre_test ( )

  call legendre_to_monomial_test ( )
  call legendre_to_monomial_matrix_test ( )
  call monomial_to_legendre_test ( )
  call monomial_to_legendre_matrix_test ( )
  call legendre_monomial_legendre_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'polynomial_conversion_test():'
  write ( *, '(a)' ) '  Normal end of execution.'

  call timestamp ( )

  stop
end
subroutine bernstein_legendre01_bernstein_test ( )

!*****************************************************************************80
!
!! bernstein_legendre01_bernstein_test() tests accuracy.
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 10

  real ( kind = rk ) bcoef(0:n)
  real ( kind = rk ) bcoef2(0:n)
  real ( kind = rk ) e
  real ( kind = rk ) lcoef(0:n)
  real ( kind = rk ) r8vec_diff_norm_l2

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'bernstein_legendre01_bernstein_test ( ):'
  write ( *, '(a)' ) '  Convert a polynomial from Bernstein form'
  write ( *, '(a)' ) '  to legendre01 form and back.'

  call random_number ( bcoef(0:n) )
  call bernstein_to_legendre01 ( n, bcoef, lcoef )
  call legendre01_to_bernstein ( n, lcoef, bcoef2 )

  e = r8vec_diff_norm_l2 ( n, bcoef, bcoef2 )
  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  L2 difference = ', e

  return
end
subroutine bernstein_monomial_bernstein_test ( )

!*****************************************************************************80
!
!! bernstein_monomial_bernstein_test() tests accuracy.
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 10

  real ( kind = rk ) bcoef(0:n)
  real ( kind = rk ) bcoef2(0:n)
  real ( kind = rk ) e
  real ( kind = rk ) mcoef(0:n)
  real ( kind = rk ) r8vec_diff_norm_l2

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'bernstein_monomial_bernstein_test ( ):'
  write ( *, '(a)' ) '  Convert a polynomial from Bernstein form'
  write ( *, '(a)' ) '  to monomial form and back.'

  call random_number ( bcoef(0:n) )
  call bernstein_to_monomial ( n, bcoef, mcoef )
  call monomial_to_bernstein ( n, mcoef, bcoef2 )

  e = r8vec_diff_norm_l2 ( n, bcoef, bcoef2 )
  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  L2 difference = ', e

  return
end
subroutine bernstein_to_legendre01_test ( )

!*****************************************************************************80
!
!! bernstein_to_legendre01_test() tests bernstein_to_legendre01().
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nmax = 6

  real ( kind = rk ) bcoef(0:nmax)
  integer k
  real ( kind = rk ) lcoef(0:nmax)
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'bernstein_to_legendre01_test ( ):'
  write ( *, '(a)' ) '  bernstein_to_legendre01() converts a'
  write ( *, '(a)' ) '  polynomial from Bernstein form'
  write ( *, '(a)' ) '  to Legendre01 form.'

  write ( *, '(a)' ) ''
  write ( *, '(10x,7(''P01'',i1,''    ''))' ) ( k, k = 0, nmax )
  write ( *, '(a)' ) ''
  do n = 0, nmax
    do k = 0, n - 1
      bcoef(k) = 0.0D+00
    end do
    bcoef(n) = 1.0D+00
    call bernstein_to_legendre01 ( n, bcoef, lcoef )
    write ( *, '(a,i1,a,f7.3,8f8.3)' ) 'B', n, '(x) = ', ( lcoef(k), k = 0, n )
  end do

  return
end
subroutine bernstein_to_legendre01_matrix_test ( )

!*****************************************************************************80
!
!! bernstein_to_legendre01_matrix_test() tests bernstein_to_legendre01_matrix().
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: A(:,:)
  integer i
  integer j
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'bernstein_to_legendre01_matrix_test ( ):'
  write ( *, '(a)' ) '  bernstein_to_legendre01_matrix() returns the matrix'
  write ( *, '(a)' ) '  which converts a polynomial from Bernstein form'
  write ( *, '(a)' ) '  to Legendre01 form.'

  n = 5

  allocate ( A(0:n,0:n) )
  call bernstein_to_legendre01_matrix ( n, A )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  A:'
  write ( *, '(a)' ) ''
  do i = 0, n
    do j = 0, n
      write ( *, '(f9.4)', advance = 'no' ) A(i,j)
    end do
    write ( *, '(a)' ) ''
  end do

  deallocate ( A )

  return
end
subroutine bernstein_to_monomial_test ( )

!*****************************************************************************80
!
!! bernstein_to_monomial_test() tests bernstein_to_monomial().
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nmax = 6

  real ( kind = rk ) bcoef(0:nmax)
  integer k
  real ( kind = rk ) mcoef(0:nmax)
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'bernstein_to_monomial_test ( ):'
  write ( *, '(a)' ) '  bernstein_to_monomial() converts a'
  write ( *, '(a)' ) '  polynomial from Bernstein form'
  write ( *, '(a)' ) '  to monomial form.'

  write ( *, '(a)' ) ''
  write ( *, '(10x,7(''X**'',i1,''    ''))' ) ( k, k = 0, nmax )
  write ( *, '(a)' ) ''
  do n = 0, nmax
    do k = 0, n - 1
      bcoef(k) = 0.0D+00
    end do
    bcoef(n) = 1.0D+00
    call bernstein_to_monomial ( n, bcoef, mcoef )
    write ( *, '(a,i1,a,f7.3,8f8.3)' ) 'B', n, '(x) = ', ( mcoef(k), k = 0, n )
  end do

  return
end
subroutine bernstein_to_monomial_matrix_test ( )

!*****************************************************************************80
!
!! bernstein_to_monomial_matrix_test() tests bernstein_to_monomial_matrix().
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: A(:,:)
  integer i
  integer j
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'bernstein_to_monomial_matrix_test ( ):'
  write ( *, '(a)' ) '  bernstein_to_monomial_matrix() returns the matrix'
  write ( *, '(a)' ) '  which converts a polynomial from Bernstein form'
  write ( *, '(a)' ) '  to monomial form.'

  n = 4

  allocate ( A(0:n,0:n) )
  call bernstein_to_monomial_matrix ( n + 1, A )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  A:'
  write ( *, '(a)' ) ''
  do i = 0, n
    do j = 0, n
      write ( *, '(f9.4)', advance = 'no' ) A(i,j)
    end do
    write ( *, '(a)' ) ''
  end do

  deallocate ( A )

  return
end
subroutine chebyshev_monomial_chebyshev_test ( )

!*****************************************************************************80
!
!! chebyshev_monomial_chebyshev_test() tests accuracy.
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
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 10

  real ( kind = rk ) ccoef(0:n)
  real ( kind = rk ) ccoef2(0:n)
  real ( kind = rk ) e
  real ( kind = rk ) mcoef(0:n)
  real ( kind = rk ) r8vec_diff_norm_l2

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'chebyshev_monomial_chebyshev_test ( ):'
  write ( *, '(a)' ) '  Convert a polynomial from Chebyshev form'
  write ( *, '(a)' ) '  to monomial form and back.'

  call random_number ( ccoef(0:n) )
  call chebyshev_to_monomial ( n, ccoef, mcoef )
  call monomial_to_chebyshev ( n, mcoef, ccoef2 )

  e = r8vec_diff_norm_l2 ( n, ccoef, ccoef2 )
  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  L2 difference = ', e

  return
end
subroutine chebyshev_to_monomial_test ( )

!*****************************************************************************80
!
!! chebyshev_to_monomial_test() tests chebyshev_to_monomial().
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
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nmax = 6

  real ( kind = rk ) ccoef(0:nmax)
  integer k
  real ( kind = rk ) mcoef(0:nmax)
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'chebyshev_to_monomial_test ( ):'
  write ( *, '(a)' ) '  chebyshev_to_monomial() converts a'
  write ( *, '(a)' ) '  polynomial from Chebyshev form'
  write ( *, '(a)' ) '  to monomial form.'

  write ( *, '(a)' ) ''
  write ( *, '(10x,7(''X**'',i1,''    ''))' ) ( k, k = 0, nmax )
  write ( *, '(a)' ) ''
  do n = 0, nmax
    do k = 0, n - 1
      ccoef(k) = 0.0D+00
    end do
    ccoef(n) = 1.0D+00
    call chebyshev_to_monomial ( n, ccoef, mcoef )
    write ( *, '(a,i1,a,f7.3,8f8.3)' ) 'T', n, '(x) = ', ( mcoef(k), k = 0, n )
  end do

  return
end
subroutine gegenbauer_monomial_gegenbauer_test ( )

!*****************************************************************************80
!
!! gegenbauer_monomial_gegenbauer_test() tests accuracy.
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
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 10

  real ( kind = rk ) alpha
  real ( kind = rk ) e
  real ( kind = rk ) gcoef(0:n)
  real ( kind = rk ) gcoef2(0:n)
  real ( kind = rk ) mcoef(0:n)
  real ( kind = rk ) r8vec_diff_norm_l2

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'gegenbauer_monomial_gegenbauer_test ( ):'
  write ( *, '(a)' ) '  Convert a polynomial from Gegenbauer form'
  write ( *, '(a)' ) '  to monomial form and back.'

  alpha = 0.5D+00
  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Gegenbauer parameter = ', alpha

  call random_number ( gcoef(0:n) )
  call gegenbauer_to_monomial ( n, alpha, gcoef, mcoef )
  call monomial_to_gegenbauer ( n, alpha, mcoef, gcoef2 )

  e = r8vec_diff_norm_l2 ( n, gcoef, gcoef2 )
  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  L2 difference = ', e

  return
end
subroutine gegenbauer_to_monomial_test ( )

!*****************************************************************************80
!
!! gegenbauer_to_monomial_test() tests gegenbauer_to_monomial().
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
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nmax = 6

  real ( kind = rk ) alpha
  real ( kind = rk ) gcoef(0:nmax)
  integer k
  real ( kind = rk ) mcoef(0:nmax)
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'gegenbauer_to_monomial_test ( ):'
  write ( *, '(a)' ) '  gegenbauer_to_monomial() converts a'
  write ( *, '(a)' ) '  polynomial from Gegenbauer form'
  write ( *, '(a)' ) '  to monomial form.'

  alpha = 0.5D+00
  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Gegenbauer parameter = ', alpha
  write ( *, '(a)' ) ''
  write ( *, '(10x,7(''X**'',i1,''    ''))' ) ( k, k = 0, nmax )
  write ( *, '(a)' ) ''
  do n = 0, nmax
    do k = 0, n - 1
      gcoef(k) = 0.0D+00
    end do
    gcoef(n) = 1.0D+00
    call gegenbauer_to_monomial ( n, alpha, gcoef, mcoef )
    write ( *, '(a,i1,a,8f8.3)' ) 'C', n, '(x) = ', ( mcoef(k), k = 0, n )
  end do

  return
end
subroutine gegenbauer_to_monomial_matrix_test ( )

!*****************************************************************************80
!
!! gegenbauer_to_monomial_matrix_test() tests gegenbauer_to_monomial_matrix().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 April 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: A(:,:)
  real ( kind = rk ) alpha
  integer i
  integer j
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'gegenbauer_to_monomial_matrix_test ( ):'
  write ( *, '(a)' ) '  gegenbauer_to_monomial_matrix() returns the matrix'
  write ( *, '(a)' ) '  which converts a polynomial from Gegenbauer form'
  write ( *, '(a)' ) '  to monomial form.'

  n = 4
  alpha = 0.5

  allocate ( A(0:n,0:n) )
  call gegenbauer_to_monomial_matrix ( n + 1, alpha, A )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  A:'
  write ( *, '(a)' ) ''
  do i = 0, n
    do j = 0, n
      write ( *, '(f9.4)', advance = 'no' ) A(i,j)
    end do
    write ( *, '(a)' ) ''
  end do

  deallocate ( A )

  return
end
subroutine hermite_monomial_hermite_test ( )

!*****************************************************************************80
!
!! hermite_monomial_hermite_test() tests accuracy.
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 10

  real ( kind = rk ) e
  real ( kind = rk ) hcoef(0:n)
  real ( kind = rk ) hcoef2(0:n)
  real ( kind = rk ) mcoef(0:n)
  real ( kind = rk ) r8vec_diff_norm_l2

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'hermite_monomial_hermite_test ( ):'
  write ( *, '(a)' ) '  Convert a polynomial from Hermite form'
  write ( *, '(a)' ) '  to monomial form and back.'

  call random_number ( hcoef(0:n) )
  call hermite_to_monomial ( n, hcoef, mcoef )
  call monomial_to_hermite ( n, mcoef, hcoef2 )

  e = r8vec_diff_norm_l2 ( n, hcoef, hcoef2 )
  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  L2 difference = ', e

  return
end
subroutine hermite_to_monomial_test ( )

!*****************************************************************************80
!
!! hermite_to_monomial_test() tests hermite_to_monomial().
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nmax = 6

  real ( kind = rk ) hcoef(0:nmax)
  integer k
  real ( kind = rk ) mcoef(0:nmax)
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'hermite_to_monomial_test ( ):'
  write ( *, '(a)' ) '  hermite_to_monomial() converts a'
  write ( *, '(a)' ) '  polynomial from Hermite form'
  write ( *, '(a)' ) '  to monomial form.'

  write ( *, '(a)' ) ''
  write ( *, '(10x,7(''X**'',i1,''    ''))' ) ( k, k = 0, nmax )
  write ( *, '(a)' ) ''
  do n = 0, nmax
    do k = 0, n - 1
      hcoef(k) = 0.0D+00
    end do
    hcoef(n) = 1.0D+00
    call hermite_to_monomial ( n, hcoef, mcoef )
    write ( *, '(a,i1,a,8f8.3)' ) 'H', n, '(x) = ', ( mcoef(k), k = 0, n )
  end do

  return
end
subroutine hermite_to_monomial_matrix_test ( )

!*****************************************************************************80
!
!! hermite_to_monomial_matrix_test() tests hermite_to_monomial_matrix().
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: A(:,:)
  integer i
  integer j
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'hermite_to_monomial_matrix_test ( ):'
  write ( *, '(a)' ) '  hermite_to_monomial_matrix() returns the matrix'
  write ( *, '(a)' ) '  which converts a polynomial from Hermite form'
  write ( *, '(a)' ) '  to monomial form.'

  n = 4

  allocate ( A(0:n,0:n) )
  call hermite_to_monomial_matrix ( n + 1, A )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  A:'
  write ( *, '(a)' ) ''
  do i = 0, n
    do j = 0, n
      write ( *, '(f9.4)', advance = 'no' ) A(i,j)
    end do
    write ( *, '(a)' ) ''
  end do

  deallocate ( A )

  return
end
subroutine laguerre_monomial_laguerre_test ( )

!*****************************************************************************80
!
!! laguerre_monomial_laguerre_test() tests accuracy.
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 10

  real ( kind = rk ) e
  real ( kind = rk ) lcoef(0:n)
  real ( kind = rk ) lcoef2(0:n)
  real ( kind = rk ) mcoef(0:n)
  real ( kind = rk ) r8vec_diff_norm_l2

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'laguerre_monomial_laguerre_test ( ):'
  write ( *, '(a)' ) '  Convert a polynomial from Laguerre form'
  write ( *, '(a)' ) '  to monomial form and back.'

  call random_number ( lcoef(0:n) )
  call laguerre_to_monomial ( n, lcoef, mcoef )
  call monomial_to_laguerre ( n, mcoef, lcoef2 )

  e = r8vec_diff_norm_l2 ( n, lcoef, lcoef2 )
  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  L2 difference = ', e

  return
end
subroutine laguerre_to_monomial_test ( )

!*****************************************************************************80
!
!! laguerre_to_monomial_test() tests laguerre_to_monomial().
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nmax = 6

  integer k
  real ( kind = rk ) lcoef(0:nmax)
  real ( kind = rk ) mcoef(0:nmax)
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'laguerre_to_monomial_test ( ):'
  write ( *, '(a)' ) '  laguerre_to_monomial() converts a'
  write ( *, '(a)' ) '  polynomial from Laguerre form'
  write ( *, '(a)' ) '  to monomial form.'

  write ( *, '(a)' ) ''
  write ( *, '(10x,7(''X**'',i1,''    ''))' ) ( k, k = 0, nmax )
  write ( *, '(a)' ) ''
  do n = 0, nmax
    do k = 0, n - 1
      lcoef(k) = 0.0D+00
    end do
    lcoef(n) = 1.0D+00
    call laguerre_to_monomial ( n, lcoef, mcoef )
    write ( *, '(a,i1,a,8f8.3)' ) 'L', n, '(x) = ', ( mcoef(k), k = 0, n )
  end do

  return
end
subroutine laguerre_to_monomial_matrix_test ( )

!*****************************************************************************80
!
!! laguerre_to_monomial_matrix_test() tests laguerre_to_monomial_matrix().
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: A(:,:)
  integer i
  integer j
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'laguerre_to_monomial_matrix_test ( ):'
  write ( *, '(a)' ) '  laguerre_to_monomial_matrix() returns the matrix'
  write ( *, '(a)' ) '  which converts a polynomial from Laguerre form'
  write ( *, '(a)' ) '  to monomial form.'

  n = 4

  allocate ( A(0:n,0:n) )
  call laguerre_to_monomial_matrix ( n + 1, A )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  A:'
  write ( *, '(a)' ) ''
  do i = 0, n
    do j = 0, n
      write ( *, '(f9.4)', advance = 'no' ) A(i,j)
    end do
    write ( *, '(a)' ) ''
  end do

  deallocate ( A )

  return
end
subroutine legendre_monomial_legendre_test ( )

!*****************************************************************************80
!
!! legendre_monomial_legendre_test() tests accuracy.
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 10

  real ( kind = rk ) e
  real ( kind = rk ) lcoef(0:n)
  real ( kind = rk ) lcoef2(0:n)
  real ( kind = rk ) mcoef(0:n)
  real ( kind = rk ) r8vec_diff_norm_l2

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'legendre_monomial_legendre_test ( ):'
  write ( *, '(a)' ) '  Convert a polynomial from Legendre form'
  write ( *, '(a)' ) '  to monomial form and back.'

  call random_number ( lcoef(0:n) )
  call legendre_to_monomial ( n, lcoef, mcoef )
  call monomial_to_legendre ( n, mcoef, lcoef2 )

  e = r8vec_diff_norm_l2 ( n, lcoef, lcoef2 )
  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  L2 difference = ', e

  return
end
subroutine legendre_to_monomial_test ( )

!*****************************************************************************80
!
!! legendre_to_monomial_test() tests legendre_to_monomial().
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nmax = 6

  integer k
  real ( kind = rk ) lcoef(0:nmax)
  real ( kind = rk ) mcoef(0:nmax)
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'legendre_to_monomial_test ( ):'
  write ( *, '(a)' ) '  legendre_to_monomial() converts a'
  write ( *, '(a)' ) '  polynomial from Legendre form'
  write ( *, '(a)' ) '  to monomial form.'

  write ( *, '(a)' ) ''
  write ( *, '(10x,7(''X**'',i1,''    ''))' ) ( k, k = 0, nmax )
  write ( *, '(a)' ) ''
  do n = 0, nmax
    do k = 0, n - 1
      lcoef(k) = 0.0D+00
    end do
    lcoef(n) = 1.0D+00
    call legendre_to_monomial ( n, lcoef, mcoef )
    write ( *, '(a,i1,a,f7.3,8f8.3)' ) 'P', n, '(x) = ', ( mcoef(k), k = 0, n )
  end do

  return
end
subroutine legendre_to_monomial_matrix_test ( )

!*****************************************************************************80
!
!! legendre_to_monomial_matrix_test() tests legendre_to_monomial_matrix().
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: A(:,:)
  integer i
  integer j
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'legendre_to_monomial_matrix_test ( ):'
  write ( *, '(a)' ) '  legendre_to_monomial_matrix() returns the matrix'
  write ( *, '(a)' ) '  which converts a polynomial from Legendre form'
  write ( *, '(a)' ) '  to monomial form.'

  n = 4

  allocate ( A(0:n,0:n) )
  call legendre_to_monomial_matrix ( n + 1, A )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  A:'
  write ( *, '(a)' ) ''
  do i = 0, n
    do j = 0, n
      write ( *, '(f12.4)', advance = 'no' ) A(i,j)
    end do
    write ( *, '(a)' ) ''
  end do

  deallocate ( A )

  return
end
subroutine legendre01_to_bernstein_test ( )

!*****************************************************************************80
!
!! legendre01_to_bernstein_test() tests legendre01_to_bernstein().
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nmax = 6

  real ( kind = rk ) bcoef(0:nmax)
  integer k
  real ( kind = rk ) lcoef(0:nmax)
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'legendre01_to_bernstein_test ( ):'
  write ( *, '(a)' ) '  legendre01_to_bernstein() converts a'
  write ( *, '(a)' ) '  polynomial from Legendre01 form to '
  write ( *, '(a)' ) '  Bernstein form.'

  write ( *, '(a)' ) ''
  write ( *, '(8x,7(''B'',i1,''(x)   ''))' ) ( k, k = 0, nmax )
  write ( *, '(a)' ) ''
  do n = 0 , nmax
    do k = 0, n - 1
      lcoef(k) = 0.0D+00
    end do
    lcoef(n) = 1.0D+00
    call legendre01_to_bernstein ( n, lcoef, bcoef )
    write ( *, '(a,i1,a,8f8.2)' ) 'P01', n, ' = ', ( bcoef(k), k = 0, n )
  end do

  return
end
subroutine legendre01_to_bernstein_matrix_test ( )

!*****************************************************************************80
!
!! legendre01_to_bernstein_matrix_test() tests legendre01_to_bernstein_matrix().
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: A(:,:)
  integer i
  integer j
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'legendre01_to_bernstein_matrix_test ( ):'
  write ( *, '(a)' ) '  legendre01_to_bernstein_matrix() returns the matrix'
  write ( *, '(a)' ) '  which converts a polynomial from Legendre01 form'
  write ( *, '(a)' ) '  to Bernstein form.'

  n = 5

  allocate ( A(0:n,0:n) )
  call legendre01_to_bernstein_matrix ( n, A )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  A:'
  write ( *, '(a)' ) ''
  do i = 0, n
    do j = 0, n
      write ( *, '(f9.4)', advance = 'no' ) A(i,j)
    end do
    write ( *, '(a)' ) ''
  end do

  deallocate ( A )

  return
end
subroutine monomial_to_bernstein_test ( )

!*****************************************************************************80
!
!! monomial_to_bernstein_test() tests monomial_to_bernstein().
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nmax = 6

  real ( kind = rk ) bcoef(0:nmax)
  integer k
  real ( kind = rk ) mcoef(0:nmax)
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'monomial_to_bernstein_test ( ):'
  write ( *, '(a)' ) '  monomial_to_bernstein() converts a'
  write ( *, '(a)' ) '  polynomial from monomial form to '
  write ( *, '(a)' ) '  Bernstein form.'

  write ( *, '(a)' ) ''
  write ( *, '(8x,7(''B'',i1,''(x)   ''))' ) ( k, k = 0, nmax )
  write ( *, '(a)' ) ''
  do n = 0 , nmax
    do k = 0, n - 1
      mcoef(k) = 0.0D+00
    end do
    mcoef(n) = 1.0D+00
    call monomial_to_bernstein ( n, mcoef, bcoef )
    write ( *, '(a,i1,a,8f8.5)' ) 'X**', n, ' = ', ( bcoef(k), k = 0, n )
  end do

  return
end
subroutine monomial_to_bernstein_matrix_test ( )

!*****************************************************************************80
!
!! monomial_to_bernstein_matrix_test() tests monomial_to_bernstein_matrix().
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: A(:,:)
  integer i
  integer j
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'monomial_to_bernstein_matrix_test ( ):'
  write ( *, '(a)' ) '  monomial_to_bernstein_matrix() returns the matrix'
  write ( *, '(a)' ) '  which converts a polynomial from monomial form'
  write ( *, '(a)' ) '  to Bernstein form.'

  n = 4

  allocate ( A(0:n,0:n) )
  call monomial_to_bernstein_matrix ( n + 1, A )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  A:'
  write ( *, '(a)' ) ''
  do i = 0, n
    do j = 0, n
      write ( *, '(f9.4)', advance = 'no' ) A(i,j)
    end do
    write ( *, '(a)' ) ''
  end do

  deallocate ( A )

  return
end
subroutine monomial_to_chebyshev_test ( )

!*****************************************************************************80
!
!! monomial_to_chebyshev_test() tests monomial_to_chebyshev().
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
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nmax = 6

  real ( kind = rk ) ccoef(0:nmax)
  integer k
  real ( kind = rk ) mcoef(0:nmax)
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'monomial_to_chebyshev_test ( ):'
  write ( *, '(a)' ) '  monomial_to_chebyshev() converts a'
  write ( *, '(a)' ) '  polynomial from monomial form to '
  write ( *, '(a)' ) '  Chebyshev form.'

  write ( *, '(a)' ) ''
  write ( *, '(8x,7(''T'',i1,''(x)   ''))' ) ( k, k = 0, nmax )
  write ( *, '(a)' ) ''
  do n = 0 , nmax
    do k = 0, n - 1
      mcoef(k) = 0.0D+00
    end do
    mcoef(n) = 1.0D+00
    call monomial_to_chebyshev ( n, mcoef, ccoef )
    write ( *, '(a,i1,a,8f8.5)' ) 'X**', n, ' = ', ( ccoef(k), k = 0, n )
  end do

  return
end
subroutine monomial_to_gegenbauer_test ( )

!*****************************************************************************80
!
!! monomial_to_gegenbauer_test() tests monomial_to_gegenbauer().
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
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nmax = 6

  real ( kind = rk ) alpha
  real ( kind = rk ) gcoef(0:nmax)
  integer k
  real ( kind = rk ) mcoef(0:nmax)
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'monomial_to_gegenbauer_test ( ):'
  write ( *, '(a)' ) '  monomial_to_gegenbauer() converts a'
  write ( *, '(a)' ) '  polynomial from monomial form to '
  write ( *, '(a)' ) '  Gegenbauer form.'

  alpha = 0.5D+00
  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Gegenbauer parameter = ', alpha

  write ( *, '(a)' ) ''
  write ( *, '(8x,7(''C'',i1,''(x)   ''))' ) ( k, k = 0, nmax )
  write ( *, '(a)' ) ''
  do n = 0 , nmax
    do k = 0, n - 1
      mcoef(k) = 0.0D+00
    end do
    mcoef(n) = 1.0D+00
    call monomial_to_gegenbauer ( n, alpha, mcoef, gcoef )
    write ( *, '(a,i1,a,8f8.5)' ) 'X**', n, ' = ', ( gcoef(k), k = 0, n )
  end do

  return
end
subroutine monomial_to_gegenbauer_matrix_test ( )

!*****************************************************************************80
!
!! monomial_to_gegenbauer_matrix_test() tests monomial_to_gegenbauer_matrix().
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
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: A(:,:)
  real ( kind = rk ) alpha
  integer i
  integer j
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'monomial_to_gegenbauer_matrix_test ( ):'
  write ( *, '(a)' ) '  monomial_to_gegenbauer_matrix() returns the matrix'
  write ( *, '(a)' ) '  which converts a polynomial from monomial form'
  write ( *, '(a)' ) '  to Gegenbauer form.'

  n = 4
  alpha = 0.5

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Gegenbauer parameter = ', alpha

  allocate ( A(0:n,0:n) )
  call monomial_to_gegenbauer_matrix ( n + 1, alpha, A )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  A:'
  write ( *, '(a)' ) ''
  do i = 0, n
    do j = 0, n
      write ( *, '(f9.4)', advance = 'no' ) A(i,j)
    end do
    write ( *, '(a)' ) ''
  end do

  deallocate ( A )

  return
end
subroutine monomial_to_hermite_test ( )

!*****************************************************************************80
!
!! monomial_to_hermite_test() tests monomial_to_hermite().
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nmax = 6

  real ( kind = rk ) hcoef(0:nmax)
  integer k
  real ( kind = rk ) mcoef(0:nmax)
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'monomial_to_hermite_test ( ):'
  write ( *, '(a)' ) '  monomial_to_hermite() converts a'
  write ( *, '(a)' ) '  polynomial from monomial form to '
  write ( *, '(a)' ) '  Hermite form.'

  write ( *, '(a)' ) ''
  write ( *, '(8x,7(''H'',i1,''(x)   ''))' ) ( k, k = 0, nmax )
  write ( *, '(a)' ) ''
  do n = 0 , nmax
    do k = 0, n - 1
      mcoef(k) = 0.0D+00
    end do
    mcoef(n) = 1.0D+00
    call monomial_to_hermite ( n, mcoef, hcoef )
    write ( *, '(a,i1,a,8f8.5)' ) 'X**', n, ' = ', ( hcoef(k), k = 0, n )
  end do

  return
end
subroutine monomial_to_hermite_matrix_test ( )

!*****************************************************************************80
!
!! monomial_to_hermite_matrix_test() tests monomial_to_hermite_matrix().
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: A(:,:)
  integer i
  integer j
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'monomial_to_hermite_matrix_test ( ):'
  write ( *, '(a)' ) '  monomial_to_hermite_matrix() returns the matrix'
  write ( *, '(a)' ) '  which converts a polynomial from monomial form'
  write ( *, '(a)' ) '  to Hermite form.'

  n = 4

  allocate ( A(0:n,0:n) )
  call monomial_to_hermite_matrix ( n + 1, A )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  A:'
  write ( *, '(a)' ) ''
  do i = 0, n
    do j = 0, n
      write ( *, '(f9.4)', advance = 'no' ) A(i,j)
    end do
    write ( *, '(a)' ) ''
  end do

  deallocate ( A )

  return
end
subroutine monomial_to_laguerre_test ( )

!*****************************************************************************80
!
!! monomial_to_laguerre_test() tests monomial_to_laguerre().
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nmax = 6

  integer k
  real ( kind = rk ) lcoef(0:nmax)
  real ( kind = rk ) mcoef(0:nmax)
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'monomial_to_laguerre_test ( ):'
  write ( *, '(a)' ) '  monomial_to_laguerre() converts a'
  write ( *, '(a)' ) '  polynomial from monomial form to '
  write ( *, '(a)' ) '  Laguerre form.'

  write ( *, '(a)' ) ''
  write ( *, '(8x,7(''L'',i1,''(x)   ''))' ) ( k, k = 0, nmax )
  write ( *, '(a)' ) ''
  do n = 0 , nmax
    do k = 0, n - 1
      mcoef(k) = 0.0D+00
    end do
    mcoef(n) = 1.0D+00
    call monomial_to_laguerre ( n, mcoef, lcoef )
    write ( *, '(a,i1,a,8g12.4)' ) 'X**', n, ' = ', ( lcoef(k), k = 0, n )
  end do

  return
end
subroutine monomial_to_laguerre_matrix_test ( )

!*****************************************************************************80
!
!! monomial_to_laguerre_matrix_test() tests monomial_to_laguerre_matrix().
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: A(:,:)
  integer i
  integer j
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'monomial_to_laguerre_matrix_test ( ):'
  write ( *, '(a)' ) '  monomial_to_laguerre_matrix() returns the matrix'
  write ( *, '(a)' ) '  which converts a polynomial from monomial form'
  write ( *, '(a)' ) '  to Laguerre form.'

  n = 4

  allocate ( A(0:n,0:n) )
  call monomial_to_laguerre_matrix ( n + 1, A )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  A:'
  write ( *, '(a)' ) ''
  do i = 0, n
    do j = 0, n
      write ( *, '(f9.4)', advance = 'no' ) A(i,j)
    end do
    write ( *, '(a)' ) ''
  end do

  deallocate ( A )

  return
end
subroutine monomial_to_legendre_test ( )

!*****************************************************************************80
!
!! monomial_to_legendre_test() tests monomial_to_legendre().
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nmax = 6

  integer k
  real ( kind = rk ) lcoef(0:nmax)
  real ( kind = rk ) mcoef(0:nmax)
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'monomial_to_legendre_test ( ):'
  write ( *, '(a)' ) '  monomial_to_legendre() converts a'
  write ( *, '(a)' ) '  polynomial from monomial form to '
  write ( *, '(a)' ) '  Legendre form.'

  write ( *, '(a)' ) ''
  write ( *, '(8x,7(''P'',i1,''(x)   ''))' ) ( k, k = 0, nmax )
  write ( *, '(a)' ) ''
  do n = 0 , nmax
    do k = 0, n - 1
      mcoef(k) = 0.0D+00
    end do
    mcoef(n) = 1.0D+00
    call monomial_to_legendre ( n, mcoef, lcoef )
    write ( *, '(a,i1,a,8f8.5)' ) 'X**', n, ' = ', ( lcoef(k), k = 0, n )
  end do

  return
end
subroutine monomial_to_legendre_matrix_test ( )

!*****************************************************************************80
!
!! monomial_to_legendre_matrix_test() tests monomial_to_legendre_matrix().
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: A(:,:)
  integer i
  integer j
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'monomial_to_legendre_matrix_test ( ):'
  write ( *, '(a)' ) '  monomial_to_legendre_matrix() returns the matrix'
  write ( *, '(a)' ) '  which converts a polynomial from monomial form'
  write ( *, '(a)' ) '  to Legendre_matrix form.'

  n = 4

  allocate ( A(0:n,0:n) )
  call monomial_to_legendre_matrix ( n + 1, A )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  A:'
  write ( *, '(a)' ) ''
  do i = 0, n
    do j = 0, n
      write ( *, '(f9.4)', advance = 'no' ) A(i,j)
    end do
    write ( *, '(a)' ) ''
  end do

  deallocate ( A )

  return
end
function r8vec_diff_norm_l2 ( n, a, b )

!*****************************************************************************80
!
!! r8vec_diff_norm_l2() returns the L2 norm of the difference of R8VEC's.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The vector L2 norm is defined as:
!
!      R8VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 April 2010
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of entries in A.
!
!    real ( kind = rk ) A(N), B(N), the vectors
!
!  Output:
!
!    real ( kind = rk ) R8VEC_DIFF_NORM_L2, the L2 norm of A - B.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) b(n)
  real ( kind = rk ) r8vec_diff_norm_l2

  r8vec_diff_norm_l2 = sqrt ( sum ( ( a(1:n) - b(1:n) ) ** 2 ) )

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
