program main

!*****************************************************************************80
!
!! toms726_test() tests toms726().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 April 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer n

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'toms726_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test toms726().'

  n = 1
  call test01 ( n )
  n = 4
  call test01 ( n )
  n = 16
  call test01 ( n )
  n = 64
  call test01 ( n )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'toms726_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( n )

!*****************************************************************************80
!
!! TEST01 calls QLAG_R8 to compute a Gauss-Laguerre quadrature rule.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 April 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the rule to be generated.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) exact
  integer i
  integer ierr
  integer j
  real ( kind = rk ) quad
  real ( kind = rk ) r8_factorial
  real ( kind = rk ), allocatable :: w(:)
  real ( kind = rk ), allocatable :: x(:)
  real ( kind = rk ), allocatable :: xj(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  QLAG_R8 computes points and weights for'
  write ( *, '(a)' ) '  a Gauss-Laguerre quadrature rule.'
!
!  Compute the rule.
!
  allocate ( x(1:n) )
  allocate ( xj(1:n) )
  allocate ( w(1:n) )

  call qlag_r8 ( n, x, w, ierr )
!
!  Print the rule.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       I      X(i)                      W(i)'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i6,2x,g25.18,2x,g25.18)' ) i, x(i), w(i)
  end do
!
!  Test the rule on a few monomials.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       J      Integral(X^J)             Quad(X^J)'
  write ( *, '(a)' ) ' '
  do j = 0, 5
    exact = r8_factorial ( j )
    xj(1:n) = x(1:n) ** j
    quad = dot_product ( w, xj )
    write ( *, '(2x,i6,2x,g25.18,2x,g25.18)' ) j, exact, quad
  end do

  deallocate ( w )
  deallocate ( x )
  deallocate ( xj )

  return
end
function r8_factorial ( n )

!*****************************************************************************80
!
!! R8_FACTORIAL computes the factorial of N.
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
subroutine qcheb_r8 ( n, x, w, i, ierr )

!*****************************************************************************80
!
!! QCHEB_R8 returns a Gauss-Chebyshev rule.
!
!  Modified:
!
!    29 April 2013
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer ierr
  integer k
  real ( kind = rk ) om2
  real ( kind = rk ) pi
  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n)

  common / common08 / om2

  pi = 4.0D+00 * atan ( 1.0D+00 )

  do k = 1, n
    x(k) = cos ( real ( 2 * k - 1, kind = rk ) * pi / real ( 2 * n, kind = rk ) )
    w(k) = pi / ( real ( n, kind = rk ) * sqrt ( 1.0D+00 - om2 * x(k)**2 ) )
  end do

  ierr = 0

  return
end
subroutine qjac_r8 ( n, x, w, i, ierr )

!*****************************************************************************80
!
!! QJAC_R8 returns a Gauss-Jacobi rule.
!
!  Modified:
!
!    29 April 2013
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(41)
  real ( kind = rk ) b(41)
  real ( kind = rk ) e(41)
  real ( kind = rk ) epsma
  integer i
  integer ierr
  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n)

  common / common05 / a, b, e, epsma

  call gauss_r8 ( n, a, b, epsma, x, w, ierr, e )

  w(1:n) = w(1:n) / b(1)

  return
end
subroutine qlag_r8 ( n, x, w, ierr )

!*****************************************************************************80
!
!! QLAG_R8 returns a Gauss-Laguerre rule.
!
!  Modified:
!
!    24 April 2013
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) alpha
  real ( kind = rk ) b(n)
  real ( kind = rk ) beta
  real ( kind = rk ) e(n)
  real ( kind = rk ) epsma
  integer ierr
  integer ipoly
  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n)
!
!  Get the recursion coefficients.
!
  ipoly = 7
  alpha = 0.0D+00
  beta = 0.0D+00

  call recur_r8 ( n, ipoly, alpha, beta, a, b, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'QLAGR8 - Fatal error!'
    write ( *, '(a,i6)' ) '  RECUR_R8 returned IERR = ', ierr
    stop
  end if
!
!  Generate the Gaussian quadrature formula.
!
  epsma = epsilon ( epsma )
  call gauss_r8 ( n, a, b, epsma, x, w, ierr, e )

  if ( ierr == 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'QLAGR8 - Warning'
    write ( *, '(a,i6)' ) '  The accuracy request was not achieved.'
  else if ( ierr /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'QLAGR8 - Fatal error!'
    write ( *, '(a,i6)' ) '  GAUSS_R8 returned IERR = ', ierr
    stop
  end if

  return
end
subroutine qchle_r8 ( n, x, w, i, ierr )

!*****************************************************************************80
!
!! QCHLE_R8 returns Gauss-Chebyshev or Gauss-Legendre rules.
!
!  Modified:
!
!    29 April 2013
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(81)
  real ( kind = rk ) b(81)
  real ( kind = rk ) c
  real ( kind = rk ) e(81)
  real ( kind = rk ) epsma
  integer i
  integer ierr
  integer k
  real ( kind = rk ) pi
  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n)

  pi = 4.0D+00 * atan ( 1.0D+00 )
!
!  Gauss-Chebyshev rule.
!
  if ( i == 1 ) then

    do k = 1, n
      x(k) = cos ( real ( 2 * k - 1, kind = rk ) * pi / real ( 2 * n, kind = rk ) )
      w(k) = pi / real ( k, kind = rk )
    end do
!
!  Gauss-Legendre rule.
!
  else if ( i == 2 ) then

    call recur_r8 ( n, 1, 0.0D+00, 0.0D+00, a, b, ierr )

    call gauss_r8 ( n, a, b, epsma, x, w, ierr, e )
!
!  Missing value for C here.
!
!   w(1:n) = c * w(1:n)

  end if

  return
end
function wf_r8 ( dx, i )

!*****************************************************************************80
!
!! wf_r8() ???
!
  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) dx
  integer i
  real ( kind = rk ) wf_r8
  wf_r8 = 1.0D+00

  return
end
