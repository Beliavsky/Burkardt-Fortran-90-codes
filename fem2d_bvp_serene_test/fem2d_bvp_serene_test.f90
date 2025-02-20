program main

!*****************************************************************************80
!
!! fem2d_bvp_serene_test() tests fem2d_bvp_serene().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 July 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'FEM2D_BVP_SERENE_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the FEM2D_BVP_SERENE library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'FEM2D_BVP_SERENE_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
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
!    01 July 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nx = 5
  integer, parameter :: ny = 5

  real ( kind = rk ), external :: a1
  real ( kind = rk ), external :: c1
  real ( kind = rk ) e1
  real ( kind = rk ) e2
  real ( kind = rk ), external :: exact1
  real ( kind = rk ), external :: exact_ux1
  real ( kind = rk ), external :: exact_uy1
  real ( kind = rk ), external :: f1
  integer fem2d_bvp_serene_node_num
  real ( kind = rk ) h1s
  integer i
  integer inc
  integer j
  integer k
  integer node_num
  real ( kind = rk ), parameter :: one = 1.0D+00
  logical show11
  real ( kind = rk ), allocatable :: u(:)
  real ( kind = rk ) uexact
  real ( kind = rk ) x(nx)
  real ( kind = rk ) y(ny)
  real ( kind = rk ), parameter :: zero = 0.0D+00

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) ' TEST01'
  write ( *, '(a)' ) '  Solve - del ( A del U ) + C U = F '
  write ( *, '(a)' ) '  on the unit square with zero boundary conditions.'
  write ( *, '(a)' ) '  A1(X,Y) = 1.0'
  write ( *, '(a)' ) '  C1(X,Y) = 0.0'
  write ( *, '(a)' ) '  F1(X,Y) = 2*X*(1-X)+2*Y*(1-Y).'
  write ( *, '(a)' ) '  U1(X,Y) = X * ( 1 - X ) * Y * ( 1 - Y )'
  write ( *, '(a)' ) ''
  write ( *, '(a,i4,a,i4,a)' ) '  The grid uses ', nx, ' by ', ny, ' nodes.'

  node_num = fem2d_bvp_serene_node_num ( nx, ny )
  allocate ( u(1:node_num) )

  write ( *, '(a,i4)' ) '  The number of nodes is ', node_num
!
!  Geometry definitions.
!
  call r8vec_linspace ( nx, zero, one, x )
  call r8vec_linspace ( ny, zero, one, y )

  show11 = .false.
  call fem2d_bvp_serene ( nx, ny, a1, c1, f1, x, y, show11, u )

  if ( nx * ny <= 25 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) '     I     J      X         Y         U         Uexact    Error'
    write ( *, '(a)' ) ''

    k = 0

    do j = 1, ny

      if ( mod ( j, 2 ) == 1 ) then
        inc = 1
      else
        inc = 2
      end if

      do i = 1, nx, inc
        k = k + 1
        uexact = exact1 ( x(i), y(j) )
        write ( *, '(2x,i4,2x,i4,2x,f8.4,2x,f8.4,2x,f8.4,2x,f8.4,2x,e10.2)' ) &
          i, j, x(i), y(j), u(k), uexact, abs ( u(k) - uexact )
      end do
    end do

  end if

  call fem2d_l1_error_serene ( nx, ny, x, y, u, exact1, e1 )
  call fem2d_l2_error_serene ( nx, ny, x, y, u, exact1, e2 )
  call fem2d_h1s_error_serene ( nx, ny, x, y, u, exact_ux1, exact_uy1, h1s )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  l1 error   = ', e1
  write ( *, '(a,g14.6)' ) '  L2 error   = ', e2
  write ( *, '(a,g14.6)' ) '  H1S error  = ', h1s

  deallocate ( u )

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
!    01 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, Y, the evaluation point.
!
!    Output, real A1, the value of A(X).
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
!    01 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, Y, the evaluation point.
!
!    Output, real C1, the value of C(X).
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
!    01 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, Y, the evaluation point.
!
!    Output, real EXACT1, the value of the solution.
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
!! EXACT_UX1 evaluates the derivative dUdX of exact solution #1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, Y, the evaluation point.
!
!    Output, real EXACT_UX1, the value of dUdX.
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
!! EXACT_UY1 evaluates the derivative dUdY of exact solution #1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, Y, the evaluation point.
!
!    Output, real EXACT_UY1, the value of dUdX.
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
!    01 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, Y, the evaluation point.
!
!    Output, real F1, the value of the right hand side.
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
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 checks the basis functions.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 July 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer j
  real ( kind = rk ) v(8,8)
  real ( kind = rk ) vx(8)
  real ( kind = rk ) vy(8)
  real ( kind = rk ) xe
  real ( kind = rk ) xq
  real ( kind = rk ) xw
  real ( kind = rk ) xx(8)
  real ( kind = rk ) yn
  real ( kind = rk ) yq
  real ( kind = rk ) ys
  real ( kind = rk ) yy(8)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Basis function checks.'
!
!  Check that V is identity matrix at nodes.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  The matrix Aij = V(i)(X(j),Y(j)) should be the identity.'
  write ( *, '(a)' ) ''

  xw = 0.0D+00
  ys = 3.0D+00
  xe = 2.0D+00
  yn = 5.0D+00
  xx(1:8) = (/ 2.0, 1.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0 /)
  yy(1:8) = (/ 5.0, 5.0, 5.0, 4.0, 3.0, 3.0, 3.0, 4.0 /)

  do j = 1, 8
    xq = xx(j)
    yq = yy(j)
    call basis_serene ( xq, yq, xw, ys, xe, yn, xx, yy, v(1:8,j) )
  end do

  call r8mat_print ( 8, 8, v, '  V(i)(X(j),Y(j))' )
!
!  Check that VX and VY sum to zero anywhere.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  The vectors dVdX(1:8)(X,Y) and dVdY(1:8)(X,Y)'
  write ( *, '(a)' ) '  should both sum to zero for any (X,Y).'

  call random_number ( harvest = xq )
  xq = 2.0 * xq
  call random_number ( harvest = yq )
  yq = 3.0 + 2.0 * yq
  xw = 0.0
  ys = 3.0
  xe = 2.0
  yn = 5.0
  xx(1:8) = (/ 2.0, 1.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0 /)
  yy(1:8) = (/ 5.0, 5.0, 5.0, 4.0, 3.0, 3.0, 3.0, 4.0 /)
  call basisd_serene ( xq, yq, xw, ys, xe, yn, xx, yy, vx, vy )
  
  write ( *, '(a)' ) ''
  write ( *, '(a,g10.4,a,g10.4,a)' ) &
    '  Random evaluation point is (', xq, ',', yq, ')'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '              dVdX        dVdY'
  write ( *, '(a)' ) ''
  do i = 1, 8
    write ( *, '(2x,i4,2x,g10.4,2x,g10.4)' ) i, vx(i), vy(i)
  end do
  write ( *, '(a)' ) ''
  write ( *, '(a,g10.4,2x,g10.4)' ) '  Sum:  ', sum ( vx(1:8) ), sum ( vy(1:8) )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 carries out test case #3.
!
!  Discussion:
!
!    Use A3, C3, F3, EXACT3, EXACT_UX3, EXACT_UY3.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 July 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nx = 5
  integer, parameter :: ny = 5

  real ( kind = rk ), external :: a3
  real ( kind = rk ), allocatable :: amat(:,:)
  real ( kind = rk ), external :: c3
  real ( kind = rk ) e1
  real ( kind = rk ) e2
  real ( kind = rk ), external :: exact3
  real ( kind = rk ), external :: exact_ux3
  real ( kind = rk ), external :: exact_uy3
  real ( kind = rk ), external :: f3
  integer fem2d_bvp_serene_node_num
  real ( kind = rk ) h1s
  integer i
  integer inc
  integer j
  integer k
  integer node_num
  real ( kind = rk ), parameter :: one = 1.0D+00
  real ( kind = rk ) scale
  logical show11
  real ( kind = rk ), allocatable :: u(:)
  real ( kind = rk ) uexact
  real ( kind = rk ) x(nx)
  real ( kind = rk ) y(ny)
  real ( kind = rk ), parameter :: zero = 0.0D+00

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Solve - del ( A del U ) + C U = F '
  write ( *, '(a)' ) '  on the unit square with zero boundary conditions.'
  write ( *, '(a)' ) '  A1(X,Y) = 0.0'
  write ( *, '(a)' ) '  C1(X,Y) = 1.0'
  write ( *, '(a)' ) '  F1(X,Y) = X * ( 1 - X ) * Y * ( 1 - Y ).'
  write ( *, '(a)' ) '  U1(X,Y) = X * ( 1 - X ) * Y * ( 1 - Y )'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  This example is contrived so that the system matrix'
  write ( *, '(a)' ) '  is the WATHEN matrix.'
  write ( *, '(a)' ) ''
  write ( *, '(a,i4,a,i4,a)' ) '  The grid uses ', nx, ' by ', ny, ' nodes.'

  node_num = fem2d_bvp_serene_node_num ( nx, ny )
  allocate ( u(1:node_num) )

  write ( *, '(a,i4)' ) '  The number of nodes is ', node_num
!
!  Geometry definitions.
!
  call r8vec_linspace ( nx, zero, one, x )
  call r8vec_linspace ( ny, zero, one, y )

  show11 = .true.

  call fem2d_bvp_serene ( nx, ny, a3, c3, f3, x, y, show11, u )

  if ( nx * ny <= 25 ) then

    write ( *, '(a)' ) ''
    write ( *, '(a)' ) '     I     J      X         Y         U         Uexact    Error'
    write ( *, '(a)' ) ''

    k = 0

    do j = 1, ny

      if ( mod ( j, 2 ) == 1 ) then
        inc = 1
      else
        inc = 2
      end if

      do i = 1, nx, inc
        k = k + 1
        uexact = exact3 ( x(i), y(j) )
        write ( *, '(2x,i4,2x,i4,2x,f8.4,2x,f8.4,2x,f8.4,2x,f8.4,2x,e10.2)' ) &
          i, j, x(i), y(j), u(k), uexact, abs ( u(k) - uexact )
      end do
    end do

  end if

  call fem2d_l1_error_serene ( nx, ny, x, y, u, exact3, e1 )
  call fem2d_l2_error_serene ( nx, ny, x, y, u, exact3, e2 )
  call fem2d_h1s_error_serene ( nx, ny, x, y, u, exact_ux3, exact_uy3, h1s )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  l1 error   = ', e1
  write ( *, '(a,g14.6)' ) '  L2 error   = ', e2
  write ( *, '(a,g14.6)' ) '  H1S error  = ', h1s
!
!  Pull out the Wathen matrix from MATLAB.
!  It will have been multiplied by a random scale factor.
!  While my numbering scheme is
!    3  2  1
!    4     8
!    5  6  7
!  the numbering scheme used here is 
!    1  2  3
!    4     5
!    6  7  8
!    
  allocate ( amat(1:8,1:8) )

  call wathen ( 1, 1, 8, amat )
 
  scale = 0.5 * amat(1,3)
  amat(1:8,1:8) = amat(1:8,1:8) / scale

  call r8mat_print ( 8, 8, amat, '  WATHEN matrix (permuted)' )

  deallocate ( amat )
  deallocate ( u )

  return
end
function a3 ( x, y )

!*****************************************************************************80
!
!! A3 evaluates A function #3.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, Y, the evaluation point.
!
!    Output, real A3, the value of A(X).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a3
  real ( kind = rk ) x
  real ( kind = rk ) y

  call r8_fake_use ( x )
  call r8_fake_use ( y )

  a3 = 0.0D+00

  return
end
function c3 ( x, y )

!*****************************************************************************80
!
!! C3 evaluates C function #3.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, Y, the evaluation point.
!
!    Output, real C3, the value of C(X).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) c3
  real ( kind = rk ) x
  real ( kind = rk ) y

  call r8_fake_use ( x )
  call r8_fake_use ( y )

  c3 = 1.0D+00

  return
end
function exact3 ( x, y )

!*****************************************************************************80
!
!! EXACT3 evaluates exact solution #3.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, Y, the evaluation point.
!
!    Output, real EXACT3, the value of the solution.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) exact3
  real ( kind = rk ) x
  real ( kind = rk ) y

  exact3 = x * ( 1.0D+00 - x ) * y * ( 1.0D+00 - y )

  return
end
function exact_ux3 ( x, y )

!*****************************************************************************80
!
!! EXACT_UX3 evaluates the derivative dUdX of exact solution #3.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, Y, the evaluation point.
!
!    Output, real EXACT_UX3, the value of dUdX.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) exact_ux3
  real ( kind = rk ) x
  real ( kind = rk ) y

  exact_ux3 = ( 1.0D+00 - 2.0D+00 * x ) * ( y - y * y )

  return
end
function exact_uy3 ( x, y )

!*****************************************************************************80
!
!! EXACT_UY3 evaluates the derivative dUdY of exact solution #3.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, Y, the evaluation point.
!
!    Output, real EXACT_UY3, the value of dUdX.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) exact_uy3
  real ( kind = rk ) x
  real ( kind = rk ) y

  exact_uy3 = ( x - x * x ) * ( 1.0D+00 - 2.0D+00 * y )

  return
end
function f3 ( x, y )

!*****************************************************************************80
!
!! F3 evaluates right hand side function #3.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, Y, the evaluation point.
!
!    Output, real F3, the value of the right hand side.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) f3
  real ( kind = rk ) x
  real ( kind = rk ) y

  f3 = x * ( 1.0D+00 - x ) * y * ( 1.0D+00 - y )

  return
end

