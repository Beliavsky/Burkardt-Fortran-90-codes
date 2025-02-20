program main

!*****************************************************************************80
!
!! FEM2D_BVP_QUADRATIC_TEST tests the FEM2D_BVP_QUADRATIC library.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 June 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM2D_BVP_QUADRATIC_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the FEM2D_BVP_QUADRATIC library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM2D_BVP_QUADRATIC_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 carries out test case #1.
!
!  Discussion:
!
!    Use A1, C1, F1, EXACT1, EXACT_UX1, EXACT_UY1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 June 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nx = 3
  integer, parameter :: ny = 3

  real ( kind = rk ), external :: a1
  real ( kind = rk ), external :: c1
  real ( kind = rk ), external :: exact1
  real ( kind = rk ), external :: exact_ux1
  real ( kind = rk ), external :: exact_uy1
  real ( kind = rk ), external :: f1
  real ( kind = rk ) e1
  real ( kind = rk ) e2
  real ( kind = rk ) h1s
  integer i
  integer j
  real ( kind = rk ) u(nx,ny)
  real ( kind = rk ) uexact
  real ( kind = rk ) x(nx)
  real ( kind = rk ) x_first
  real ( kind = rk ) x_last
  real ( kind = rk ) y(nx)
  real ( kind = rk ) y_first
  real ( kind = rk ) y_last

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Solve - del ( A del U ) + C U = F '
  write ( *, '(a)' ) '  on the unit square with zero boundary conditions.'
  write ( *, '(a)' ) '  A1(X,Y) = 1.0'
  write ( *, '(a)' ) '  C1(X,Y) = 0.0'
  write ( *, '(a)' ) '  F1(X,Y) = 2*X*(1-X)+2*Y*(1-Y)'
  write ( *, '(a)' ) '  U1(X,Y) = X * ( 1 - X ) * Y * ( 1 - Y )'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of X grid values NX = ', nx
  write ( *, '(a,i8)' ) '  Number of Y grid values NY = ', ny
!
!  Geometry definitions.
!
  x_first = 0.0D+00
  x_last = 1.0D+00
  call r8vec_even ( nx, x_first, x_last, x )

  y_first = 0.0D+00
  y_last = 1.0D+00
  call r8vec_even ( ny, y_first, y_last, y )

  call fem2d_bvp_quadratic ( nx, ny, a1, c1, f1, x, y, u )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I     J    X         Y         U         Uexact    Error'
  write ( *, '(a)' ) ' '

  do j = 1, ny
    do i = 1, nx
      uexact = exact1 ( x(i), y(j) )
      write ( *, '(2x,i4,2x,i4,2x,f8.2,2x,f8.2,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
        i, j, x(i), y(j), u(i,j), uexact, abs ( u(i,j) - uexact )
    end do
  end do

  call fem2d_l1_error ( nx, ny, x, y, u, exact1, e1 )
  call fem2d_l2_error_quadratic ( nx, ny, x, y, u, exact1, e2 )
  call fem2d_h1s_error_quadratic ( nx, ny, x, y, u, exact_ux1, exact_uy1, h1s )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  l1 norm of error  = ', e1
  write ( *, '(a,g14.6)' ) '  L2 norm of error  = ', e2
  write ( *, '(a,g14.6)' ) '  Seminorm of error = ', h1s

  return
end
function a1 ( x, y )

!*****************************************************************************80
!
!! A1 evaluates A function #1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 June 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, Y, the evaluation point.
!
!    Output, real ( kind = rk ) A1, the value of A(X,Y).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a1
  real ( kind = rk ) x
  real ( kind = rk ) y

  call r8_fake_use ( x )
  call r8_fake_use ( y )

  a1 = 1.0D+00

  return
end
function c1 ( x, y )

!*****************************************************************************80
!
!! C1 evaluates C function #1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 June 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, Y, the evaluation point.
!
!    Output, real ( kind = rk ) C1, the value of C(X,Y).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) c1
  real ( kind = rk ) x
  real ( kind = rk ) y

  call r8_fake_use ( x )
  call r8_fake_use ( y )

  c1 = 0.0D+00

  return
end
function exact1 ( x, y )

!*****************************************************************************80
!
!! EXACT1 evaluates exact solution #1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 June 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, Y, the evaluation point.
!
!    Output, real ( kind = rk ) EXACT1, the value of U(X,Y).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) exact1
  real ( kind = rk ) x
  real ( kind = rk ) y

  exact1 = x * ( 1.0D+00 - x ) * y * ( 1.0D+00 - y )

  return
end
function exact_ux1 ( x, y )

!*****************************************************************************80
!
!! EXACT_UX1 evaluates the X derivative of the exact solution #1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 June 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, Y, the evaluation point.
!
!    Output, real ( kind = rk ) EXACT_UX1, the value of dUdX(X,Y).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) exact_ux1
  real ( kind = rk ) x
  real ( kind = rk ) y

  exact_ux1 = ( 1.0D+00 - 2.0D+00 * x ) * ( y - y * y )

  return
end
function exact_uy1 ( x, y )

!*****************************************************************************80
!
!! EXACT_UY1 evaluates the Y derivative of the exact solution #1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 June 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, Y, the evaluation point.
!
!    Output, real ( kind = rk ) EXACT_UY1, the value of dUdY(X,Y).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) exact_uy1
  real ( kind = rk ) x
  real ( kind = rk ) y

  exact_uy1 = ( x - x * x ) * ( 1.0D+00 - 2.0D+00 * y )

  return
end
function f1 ( x, y )

!*****************************************************************************80
!
!! F1 evaluates right hand side function #1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 June 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, Y, the evaluation point.
!
!    Output, real ( kind = rk ) F1, the value of F(X,Y).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) f1
  real ( kind = rk ) x
  real ( kind = rk ) y

  f1 = 2.0D+00 * x * ( 1.0D+00 - x ) &
     + 2.0D+00 * y * ( 1.0D+00 - y )

  return
end

