program main

!*****************************************************************************80
!
!! quadrature_weights_vandermonde_test() tests quadrature_weights_vandermonde().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 February 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'quadrature_weights_vandermonde_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the quadrature_weights_vandermonde library.'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'quadrature_weights_vandermonde_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests quadrature_weights_vandermonde for a Newton-Cotes rule.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 April 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: x(:)

  a =  0.0D+00
  b = +1.0D+00
  n = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  Use the Vandermonde procedure to compute the'
  write ( *, '(a)' ) '  quadrature weights for a Newton-Cotes rule.'
  write ( *, '(a,i4)' ) '  Order N = ', n
  write ( *, '(a,f10.4,a,f10.4,a)' ) '  Interval = [', a, ',', b, ']'
!
!  Set the points.
!
  allocate ( x(1:n) )
  call r8vec_even ( n, a, b, x )
  call r8vec_print ( n, x, '  Abscissas:' )
!
!  Compute the weights.
!
  allocate ( w(1:n) )

  call quadrature_weights_vandermonde ( n, a, b, x, w )

  call r8vec_print ( n, w, '  Weights:' )
!
!  Free memory.
!
  deallocate ( w )
  deallocate ( x )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests quadrature_weights_vandermonde for a Clenshaw-Curtis rule.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 April 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) theta
  real ( kind = 8 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: x(:)

  a =  -1.0D+00
  b = +1.0D+00
  n = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Use the Vandermonde procedure to compute the'
  write ( *, '(a)' ) '  quadrature weights for a Clenshaw-Curtis rule.'
  write ( *, '(a,i4)' ) '  Order N = ', n
  write ( *, '(a,f10.4,a,f10.4,a)' ) '  Interval is [', a, ',', b, ']'
!
!  Set the points.
!
  allocate ( x(1:n) )

  do i = 1, n

    theta = real ( n - i, kind = 8 ) * pi &
          / real ( n - 1, kind = 8 )

    x(i) = ( ( 1 - cos ( theta ) ) * a   &
           + ( 1 + cos ( theta ) ) * b ) &
           /   2.0D+00

  end do

  call r8vec_print ( n, x, '  Abscissas:' )
!
!  Determine the corresponding weights.
!
  allocate ( w(1:n) )

  call quadrature_weights_vandermonde ( n, a, b, x, w )

  call r8vec_print ( n, w, '  Weights:' )
!
!  Free memory.
!
  deallocate ( w )
  deallocate ( x )

  return
end
