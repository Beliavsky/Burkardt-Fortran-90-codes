program main

!*****************************************************************************80
!
!! FEM1D_LAGRANGE_TEST tests FEM1D_LAGRANGE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 November 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer q_num
  integer x_num

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'FEM1D_LAGRANGE_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) '  Test the FEM1D_LAGRANGE library.'

  call legendre_set_test ( )
  call lagrange_value_test ( )
  call lagrange_derivative_test ( )

  x_num = 11
  q_num = 5
  call fem1d_lagrange_stiffness_test ( x_num, q_num )

  x_num = 11
  q_num = 10
  call fem1d_lagrange_stiffness_test ( x_num, q_num )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'FEM1D_LAGRANGE_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine legendre_set_test ( )

!*****************************************************************************80
!
!! LEGENDRE_SET_TEST tests LEGENDRE_SET.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 November 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) e1
  real ( kind = rk ) e2
  real ( kind = rk ) e3
  integer n
  real ( kind = rk ), allocatable :: w(:)
  real ( kind = rk ), allocatable :: x(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'LEGENDRE_SET_TEST'
  write ( *, '(a)' ) '  LEGENDRE_SET returns points and weights of'
  write ( *, '(a)' ) '  Gauss-Legendre quadrature rules.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '   N               1             X^4           Runge'
  write ( *, '(a)' ) ''

  do n = 1, 10

    allocate ( x(1:n) )
    allocate ( w(1:n) )

    call legendre_set ( n, x, w )
    e1 = sum ( w(1:n) )
    e2 = dot_product ( w(1:n), x(1:n) ** 4 )
    e3 = dot_product ( w(1:n), 1.0D+00 / ( 1.0D+00 + 25.0D+00 * x(1:n) ** 2 ) )
    write ( *, '(2x,i2,2x,g14.6,2x,g14.6,2x,g14.6)' ) n, e1, e2, e3

    deallocate ( w )
    deallocate ( x )

  end do

  return
end
subroutine lagrange_value_test ( )

!*****************************************************************************80
!
!! LAGRANGE_VALUE_TEST tests LAGRANGE_VALUE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 November 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer j
  real ( kind = rk ), allocatable :: li(:,:)
  integer nd
  integer ni
  real ( kind = rk ), allocatable :: xd(:)
  real ( kind = rk ) xhi
  real ( kind = rk ), allocatable :: xi(:)
  real ( kind = rk ) xlo

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'LAGRANGE_VALUE_TEST'
  write ( *, '(a)' ) '  LAGRANGE_VALUE evaluates the Lagrange basis polynomials.'

  nd = 5
  xlo = 0.0D+00
  xhi = real ( nd - 1, kind = rk ) 
  allocate ( xd(1:nd) )
  call r8vec_linspace ( nd, xlo, xhi, xd )

  call r8vec_print ( nd, xd, '  Lagrange basis points:' )
!
!  Evaluate the polynomials.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) &
    '   I      X          L1(X)       L2(X)       L3(X)       L4(X)       L5(X)'
  write ( *, '(a)' ) ''
 
  ni = 2 * nd - 1
  allocate ( xi(1:ni) )
  call r8vec_linspace ( ni, xlo, xhi, xi )

  allocate ( li(1:ni,1:nd) )
  call lagrange_value ( nd, xd, ni, xi, li )

  do i = 1, ni
 
    write ( *, '(2x,i2,f10.4)', advance = 'no' ) i, xi(i)
    do j = 1, nd
      write ( *, '(2x,f10.4)', advance = 'no' ) li(i,j)
    end do
    write ( *, '(a)' ) ''
 
  end do

  deallocate ( li )
  deallocate ( xd )
  deallocate ( xi )

  return
end
subroutine lagrange_derivative_test ( )

!*****************************************************************************80
!
!! LAGRANGE_DERIVATIVE_TEST tests LAGRANGE_DERIVATIVE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 November 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer j
  real ( kind = rk ), allocatable :: lpi(:,:)
  integer nd
  integer ni
  real ( kind = rk ), allocatable :: xd(:)
  real ( kind = rk ) xhi
  real ( kind = rk ), allocatable :: xi(:)
  real ( kind = rk ) xlo

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'LAGRANGE_DERIVATIVE_TEST'
  write ( *, '(a)' ) '  LAGRANGE_DERIVATIVE evaluates the Lagrange basis derivative.'

  nd = 5
  xlo = 0.0D+00
  xhi = real ( nd - 1, kind = rk ) 
  allocate ( xd(1:nd) )
  call r8vec_linspace ( nd, xlo, xhi, xd )

  call r8vec_print ( nd, xd, '  Lagrange basis points:' )
!
!  Evaluate the polynomials.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) &
    '   I      X         L1''(X)      L2''(X)      L3''(X)      L4''(X)      L5''(X)'
  write ( *, '(a)' ) ''
 
  ni = 2 * nd - 1
  allocate ( xi(1:ni) )
  call r8vec_linspace ( ni, xlo, xhi, xi )

  allocate ( lpi(1:ni,1:nd) )
  call lagrange_derivative ( nd, xd, ni, xi, lpi )

  do i = 1, ni
 
    write ( *, '(2x,i2,f10.4)', advance = 'no' ) i, xi(i)
    do j = 1, nd
      write ( *, '(2x,f10.4)', advance = 'no' ) lpi(i,j)
    end do
    write ( *, '(a)' ) ''
 
  end do

  deallocate ( lpi )
  deallocate ( xd )
  deallocate ( xi )

  return
end
subroutine fem1d_lagrange_stiffness_test ( x_num, q_num )

!*****************************************************************************80
!
!! FEM1D_LAGRANGE_STIFFNESS_TEST tests FEM1D_LAGRANGE_STIFFNESS.
!
!  Discussion:
!
!    The results are very sensitive to the quadrature rule.
!
!    In particular, if X_NUM points are used, the mass matrix will
!    involve integrals of polynomials of degree 2*(X_NUM-1), so the
!    quadrature rule should use at least Q_NUM = X_NUM - 1 points.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 November 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer X_NUM, the number of nodes.
!
!    Input, integer Q_NUM, the number of quadrature points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: a(:,:)
  real ( kind = rk ), allocatable :: b(:)
  real ( kind = rk ), external :: f
  integer i
  integer info
  real ( kind = rk ), allocatable :: k(:,:)
  real ( kind = rk ), allocatable :: m(:,:)
  integer q_num
  real ( kind = rk ), allocatable :: u(:)
  real ( kind = rk ), allocatable :: u_e(:)
  real ( kind = rk ), allocatable :: x(:)
  real ( kind = rk ) x_hi
  real ( kind = rk ) x_lo
  integer x_num

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'FEM1D_LAGRANGE_STIFFNESS_TEST'
  write ( *, '(a)' ) '  FEM1D_LAGRANGE_STIFFNESS computes the stiffness matrix,'
  write ( *, '(a)' ) '  the mass matrix, and right hand side vector for a'
  write ( *, '(a)' ) '  finite element problem using Lagrange interpolation'
  write ( *, '(a)' ) '  basis polynomials.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Solving:'
  write ( *, '(a)' ) '    -u"+u=x on 0 < x < 1'
  write ( *, '(a)' ) '    u(0) = u(1) = 0'
  write ( *, '(a)' ) '  Exact solution:'
  write ( *, '(a)' ) '    u(x) = x - sinh(x)/sinh(1)'
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Number of mesh points = ', x_num
  write ( *, '(a,i4)' ) '  Number of quadrature points = ', q_num

  x_lo = 0.0D+00
  x_hi = 1.0D+00
  allocate ( x(1:x_num) )
  call r8vec_linspace ( x_num, x_lo, x_hi, x )

  allocate ( a(1:x_num,1:x_num) )
  allocate ( m(1:x_num,1:x_num) )
  allocate ( b(1:x_num) )

  call fem1d_lagrange_stiffness ( x_num, x, q_num, f, a, m, b )

  allocate ( k(1:x_num,1:x_num) )

  k = a + m

  k(1,1:x_num) = 0.0D+00
  k(1,1) = 1.0D+00
  b(1) = 0.0D+00

  k(x_num,1:x_num) = 0.0D+00
  k(x_num,x_num) = 1.0D+00
  b(x_num) = 0.0D+00

  allocate ( u(1:x_num) )

  u(1:x_num) = b(1:x_num)
  call r8mat_fs ( x_num, k, u, info )

  allocate ( u_e(1:x_num) )
  call exact ( x_num, x, u_e )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '   I      X             U              U(exact)         Error'
  write ( *, '(a)' ) ''

  do i = 1, x_num
    write ( *, '(2x,i2,2x,f8.4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      i, x(i), u(i), u_e(i), abs ( u(i) - u_e(i) )
  end do

  deallocate ( a )
  deallocate ( b )
  deallocate ( k )
  deallocate ( m )
  deallocate ( u )
  deallocate ( u_e )
  deallocate ( x )

  return
end
function f ( x )

!*****************************************************************************80
!
!! F evaluates the right hand side function.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 November 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, the evaluation point.
!
!    Output, real ( kind = rk ) F, the value of the right hand side at X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) f
  real ( kind = rk ) x

  f = x

  return
end
subroutine exact ( x_num, x, ue )

!*****************************************************************************80
!
!! EXACT returns the exact solution.
!
!  Discussion:
!
!    The results are very sensitive to the quadrature rule.
!
!    In particular, if X_NUM points are used, the mass matrix will
!    involve integrals of polynomials of degree 2*(X_NUM-1), so the
!    quadrature rule should use at least Q_NUM = X_NUM - 1 points.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 November 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer X_NUM, the number of nodes.
!
!    Input, real ( kind = rk ) X(X_NUM), the nodes.
!
!    Output, real ( kind = rk ) UE(X_NUM), the exact solution at the nodes.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer x_num

  real ( kind = rk ) ue(x_num)
  real ( kind = rk ) x(x_num)

  ue(1:x_num) = x(1:x_num) - sinh ( x(1:x_num) ) / sinh ( 1.0D+00 )

  return
end

