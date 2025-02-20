subroutine fem2d_bvp_quadratic ( nx, ny, a, c, f, x, y, u )

!*****************************************************************************80
!
!! fem2d_bvp_quadratic() solves a boundary value problem on a rectangle.
!
!  Discussion:
!
!    The procedure uses the finite element method, with piecewise quadratic
!    basis functions to solve a 2D boundary value problem over a rectangle.
!
!    The following differential equation is imposed inside the region:
!
!      - d/dx a(x,y) du/dx - d/dy a(x,y) du/dy + c(x,y) * u(x,y) = f(x,y)
!
!    where a(x,y), c(x,y), and f(x,y) are given functions.
!
!    On the boundary, the solution is constrained to have the value 0,
!    known as "zero Dirichlet" boundary conditions.
!
!    The finite element method will use a regular grid of NX nodes in X, and 
!    NY nodes in Y.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 June 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NX, NY, the number of X and Y grid values.
!    NX and NY must be odd and at least 3.
!
!    Input, function A ( X ), evaluates a(x);
!
!    Input, function C ( X ), evaluates c(x);
!
!    Input, function F ( X ), evaluates f(x);
!
!    Input, real ( kind = rk ) X(NX), Y(NY), the mesh points.
!
!  Output:
!
!    Output, real ( kind = rk ) U(NX,NY), the finite element coefficients, which 
!    are also the value of the computed solution at the mesh points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nx
  integer ny
  integer, parameter :: quad_num = 3

  real ( kind = rk ), external :: a
  real ( kind = rk ), dimension ( quad_num ) :: abscissa = (/ &
    -0.774596669241483377035853079956D+00, &
     0.000000000000000000000000000000D+00, &
     0.774596669241483377035853079956D+00 /)
  real ( kind = rk ) amat(nx*ny,nx*ny)
  real ( kind = rk ) aq
  real ( kind = rk ) b(nx*ny)
  real ( kind = rk ), external :: c
  integer cc
  real ( kind = rk ) cq
  integer e
  integer ex
  integer ex_num
  integer ey
  integer ey_num
  real ( kind = rk ), external :: f
  real ( kind = rk ) fq
  integer i
  integer ierror
  integer ii
  integer il
  integer il2
  integer il3
  integer j
  integer jj
  integer jl
  integer jl2
  integer jl3
  integer k
  integer mm
  integer mn
  integer n
  integer node(9)
  integer qx
  integer qy
  integer s
  real ( kind = rk ) t
  real ( kind = rk ) u(nx*ny)
  real ( kind = rk ) v(9)
  real ( kind = rk ) vx(9)
  real ( kind = rk ) vy(9)
  integer w
  real ( kind = rk ), dimension ( quad_num ) :: weight = (/ &
    0.555555555555555555555555555556D+00, &
    0.888888888888888888888888888889D+00, &
    0.555555555555555555555555555556D+00 /)
  real ( kind = rk ) wq
  real ( kind = rk ) x(nx)
  real ( kind = rk ) xq
  real ( kind = rk ) xx(3)
  real ( kind = rk ) y(ny)
  real ( kind = rk ) yq
  real ( kind = rk ) yy(3)
!
!  Zero out the matrix and right hand side.
!
  mn = nx * ny
  amat(1:mn,1:mn) = 0.0D+00
  b(1:mn) = 0.0D+00

  ex_num = ( nx - 1 ) / 2
  ey_num = ( ny - 1 ) / 2

  do ex = 1, ex_num

    w = 2 * ex - 1
    cc = 2 * ex
    e = 2 * ex + 1

    xx(1) = x(w)
    xx(2) = x(cc)
    xx(3) = x(e)

    do ey = 1, ey_num

      s = 2 * ey - 1
      mm = 2 * ey
      n = 2 * ey + 1

      yy(1) = y(s)
      yy(2) = y(mm)
      yy(3) = y(n)
!
!  Node indices
!
!  7  8  9   wn cn en
!  4  5  6   wm cm em
!  1  2  3   ws cs es
!
      node(1) = ( 2 * ey - 2 ) * nx + ( ex - 1 ) * 2 + 1
      node(2) = ( 2 * ey - 2 ) * nx + ( ex - 1 ) * 2 + 2
      node(3) = ( 2 * ey - 2 ) * nx + ( ex - 1 ) * 2 + 3
      node(4) = ( 2 * ey - 1 ) * nx + ( ex - 1 ) * 2 + 1
      node(5) = ( 2 * ey - 1 ) * nx + ( ex - 1 ) * 2 + 2
      node(6) = ( 2 * ey - 1 ) * nx + ( ex - 1 ) * 2 + 3
      node(7) = ( 2 * ey     ) * nx + ( ex - 1 ) * 2 + 1
      node(8) = ( 2 * ey     ) * nx + ( ex - 1 ) * 2 + 2
      node(9) = ( 2 * ey     ) * nx + ( ex - 1 ) * 2 + 3

      do qx = 1, quad_num

        xq = ( ( 1.0D+00 - abscissa(qx) ) * xx(1)  &
             + ( 1.0D+00 + abscissa(qx) ) * xx(3) ) &
               / 2.0D+00

        do qy = 1, quad_num

          yq = ( ( 1.0D+00 - abscissa(qy) ) * yy(1)   &
               + ( 1.0D+00 + abscissa(qy) ) * yy(3) ) &
                 / 2.0D+00

          wq = weight(qx) * ( xx(3) - xx(1) ) / 2.0D+00 &
             * weight(qy) * ( yy(3) - yy(1) ) / 2.0D+00


          v(1:9) = 1.0D+00
          vx(1:9) = 0.0D+00
          vy(1:9) = 0.0D+00

          k = 0

          do jl = 1, 3
            do il = 1, 3

              k = k + 1

              do il2 = 1, 3
                if ( il2 /= il ) then
                  v(k) = v(k) * ( xq - xx(il2) ) / ( xx(il) - xx(il2) )
                  t = 1.0D+00 / ( xx(il) - xx(il2 ) )
                  do il3 = 1, 3
                    if ( il3 /= il .and. il3 /= il2 ) then
                      t = t * ( xq - xx(il3) ) / ( xx(il) - xx(il3) )
                    end if
                  end do
                  do jl2 = 1, 3
                    if ( jl2 /= jl ) then
                      t = t * ( yq - yy(jl2) ) / ( yy(jl) - yy(jl2) )
                    end if
                  end do
                  vx(k) = vx(k) + t
                end if
              end do

              do jl2 = 1, 3
                if ( jl2 /= jl ) then
                  v(k) = v(k) * ( yq - yy(jl2) ) / ( yy(jl) - yy(jl2) )
                  t = 1.0D+00 / ( yy(jl) - yy(jl2 ) )
                  do il2 = 1, 3
                    if ( il2 /= il ) then
                      t = t * ( xq - xx(il2) ) / ( xx(il) - xx(il2) )
                    end if
                  end do
                  do jl3 = 1, 3
                    if ( jl3 /= jl .and. jl3 /= jl2 ) then
                      t = t * ( yq - yy(jl3) ) / ( yy(jl) - yy(jl3) )
                    end if
                  end do
                  vy(k) = vy(k) + t
                end if
              end do

            end do
          end do

          aq = a ( xq, yq )
          cq = c ( xq, yq )
          fq = f ( xq, yq )

          do i = 1, 9
            ii = node(i)
            do j = 1, 9
              jj = node(j)
              amat(ii,jj) = amat(ii,jj) + wq * ( vx(i) * aq * vx(j) &
                                               + vy(i) * aq * vy(j) &
                                               + v(i)  * cq * v(j) )
            end do
            b(ii) = b(ii) + wq * ( v(i) * fq )
          end do
 
        end do
      end do
    end do
  end do
!
!  Where a node is on the boundary, 
!  replace the finite element equation by a boundary condition.
!
  k = 0
  do j = 1, ny
    do i = 1, nx
      k = k + 1
      if ( i == 1 .or. i == nx .or. j == 1 .or. j == ny ) then
        amat(k,1:mn) = 0.0D+00
        amat(1:mn,k) = 0.0D+00
        amat(k,k) = 1.0D+00
        b(k) = 0.0D+00
      end if
    end do
  end do
!
!  Solve the linear system.
!
  call r8mat_solve2 ( mn, amat, b, u, ierror )

  return
end
subroutine fem2d_h1s_error_quadratic ( nx, ny, x, y, u, exact_ux, &
  exact_uy, h1s )

!*****************************************************************************80
!
!! FEM2D_H1S_ERROR_QUADRATIC: seminorm error of a finite element solution.
!
!  Discussion:
!
!    The finite element method has been used, over a rectangle,
!    involving a grid of NX*NY nodes, with piecewise quadratic elements used 
!    for the basis.
!
!    The finite element solution U(x,y) has been computed, and formulas for the
!    exact derivatives Vx(x,y) and Vy(x,y) are known.
!
!    This function estimates the H1 seminorm of the error:
!
!      H1S = sqrt ( integral ( x, y )   ( Ux(x,y) - Vx(x,y) )^2 
!                                     + ( Uy(x,y) - Vy(x,y) )^2 dx dy )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 June 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NX, NY, the number of nodes in the X 
!    and Y directions.
!
!    Input, real ( kind = rk ) X(NX), Y(NY), the grid coordinates.
!
!    Input, real ( kind = rk ) U(NX,NY), the finite element coefficients.
!
!    Input, function EXACT_UX(X,Y), EXACT_UY(X,Y), return the exact
!    derivatives with respect to X and Y, at the point (X,Y).
!
!    Output, real ( kind = rk ) H1S, the estimated seminorm of 
!    the error.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nx
  integer ny
  integer, parameter :: quad_num = 3

  real ( kind = rk ), dimension ( quad_num ) :: abscissa = (/ &
    -0.774596669241483377035853079956D+00, &
     0.000000000000000000000000000000D+00, &
     0.774596669241483377035853079956D+00 /)
  integer cc
  integer e
  integer ex
  integer ex_num
  integer ey
  integer ey_num
  real ( kind = rk ), external :: exact_ux
  real ( kind = rk ), external :: exact_uy
  real ( kind = rk ) exq
  real ( kind = rk ) eyq
  real ( kind = rk ) h1s
  integer il
  integer il2
  integer il3
  integer jl
  integer jl2
  integer jl3
  integer k
  integer mm
  integer n
  integer node(9)
  integer qx
  integer qy
  integer s
  real ( kind = rk ) t
  real ( kind = rk ) u(nx*ny)
  real ( kind = rk ) uxq
  real ( kind = rk ) uyq
  real ( kind = rk ) vx(9)
  real ( kind = rk ) vy(9)
  integer w
  real ( kind = rk ), dimension ( quad_num ) :: weight = (/ &
    0.555555555555555555555555555556D+00, &
    0.888888888888888888888888888889D+00, &
    0.555555555555555555555555555556D+00 /)
  real ( kind = rk ) wq
  real ( kind = rk ) x(nx)
  real ( kind = rk ) xq
  real ( kind = rk ) xx(3)
  real ( kind = rk ) y(ny)
  real ( kind = rk ) yq
  real ( kind = rk ) yy(3)

  h1s = 0.0D+00

  ex_num = ( nx - 1 ) / 2
  ey_num = ( ny - 1 ) / 2

  do ex = 1, ex_num

    w = 2 * ex - 1
    cc = 2 * ex
    e = 2 * ex + 1

    xx(1) = x(w)
    xx(2) = x(cc)
    xx(3) = x(e)

    do ey = 1, ey_num

      s = 2 * ey - 1
      mm = 2 * ey
      n = 2 * ey + 1

      yy(1) = y(s)
      yy(2) = y(mm)
      yy(3) = y(n)
!
!  Node indices
!
!  7  8  9   wn cn en
!  4  5  6   wm cm em
!  1  2  3   ws cs es
!
      node(1) = ( 2 * ey - 2 ) * nx + ( ex - 1 ) * 2 + 1
      node(2) = ( 2 * ey - 2 ) * nx + ( ex - 1 ) * 2 + 2
      node(3) = ( 2 * ey - 2 ) * nx + ( ex - 1 ) * 2 + 3
      node(4) = ( 2 * ey - 1 ) * nx + ( ex - 1 ) * 2 + 1
      node(5) = ( 2 * ey - 1 ) * nx + ( ex - 1 ) * 2 + 2
      node(6) = ( 2 * ey - 1 ) * nx + ( ex - 1 ) * 2 + 3
      node(7) = ( 2 * ey     ) * nx + ( ex - 1 ) * 2 + 1
      node(8) = ( 2 * ey     ) * nx + ( ex - 1 ) * 2 + 2
      node(9) = ( 2 * ey     ) * nx + ( ex - 1 ) * 2 + 3

      do qx = 1, quad_num

        xq = ( ( 1.0D+00 - abscissa(qx) ) * xx(1)  &
             + ( 1.0D+00 + abscissa(qx) ) * xx(3) ) &
               / 2.0D+00

        do qy = 1, quad_num

          yq = ( ( 1.0D+00 - abscissa(qy) ) * yy(1)   &
               + ( 1.0D+00 + abscissa(qy) ) * yy(3) ) &
                 / 2.0D+00

          wq = weight(qx) * ( xx(3) - xx(1) ) / 2.0D+00 &
             * weight(qy) * ( yy(3) - yy(1) ) / 2.0D+00

          vx(1:9) = 0.0D+00
          vy(1:9) = 0.0D+00

          uxq = 0.0D+00
          uyq = 0.0D+00

          k = 0

          do jl = 1, 3
            do il = 1, 3

              k = k + 1

              do il2 = 1, 3
                if ( il2 /= il ) then
                  t = 1.0D+00 / ( xx(il) - xx(il2 ) )
                  do il3 = 1, 3
                    if ( il3 /= il .and. il3 /= il2 ) then
                      t = t * ( xq - xx(il3) ) / ( xx(il) - xx(il3) )
                    end if
                  end do
                  do jl2 = 1, 3
                    if ( jl2 /= jl ) then
                      t = t * ( yq - yy(jl2) ) / ( yy(jl) - yy(jl2) )
                    end if
                  end do
                  vx(k) = vx(k) + t
                end if
              end do

              do jl2 = 1, 3
                if ( jl2 /= jl ) then

                  t = 1.0D+00 / ( yy(jl) - yy(jl2 ) )
                  do il2 = 1, 3
                    if ( il2 /= il ) then
                      t = t * ( xq - xx(il2) ) / ( xx(il) - xx(il2) )
                    end if
                  end do
                  do jl3 = 1, 3
                    if ( jl3 /= jl .and. jl3 /= jl2 ) then
                      t = t * ( yq - yy(jl3) ) / ( yy(jl) - yy(jl3) )
                    end if
                  end do
                  vy(k) = vy(k) + t
                end if
              end do

              uxq = uxq + u(node(k)) * vx(k)
              uyq = uyq + u(node(k)) * vy(k)

            end do
          end do

          exq = exact_ux ( xq, yq )
          eyq = exact_uy ( xq, yq )

          h1s = h1s + wq * ( ( uxq - exq ) ** 2 + ( uyq - eyq ) ** 2 )

        end do
      end do
    end do
  end do

  h1s = sqrt ( h1s )

  return
end
subroutine fem2d_l1_error ( nx, ny, x, y, u, exact, e1 )

!*****************************************************************************80
!
!! FEM2D_L1_ERROR estimates the l1 error norm of a finite element solution.
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
!    Input, integer NX, NY, the number of nodes in the X 
!    and Y directions.
!
!    Input, real ( kind = rk ) X(NX), Y(NY), the grid coordinates.
!
!    Input, real ( kind = rk ) U(NX,NY), the finite element coefficients.
!
!    Input, function EQ = EXACT ( X, Y ), returns the value of the exact
!    solution at the point (X,Y).
!
!    Output, real ( kind = rk ) E1, the estimated L2 norm of the error.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nx
  integer ny

  real ( kind = rk ) e1
  real ( kind = rk ), external :: exact
  integer i
  integer j
  real ( kind = rk ) u(nx,ny)
  real ( kind = rk ) x(nx)
  real ( kind = rk ) y(ny)

  e1 = 0.0D+0
  do j = 1, ny
    do i = 1, nx
      e1 = e1 + abs ( u(i,j) - exact ( x(i), y(j) ) )
    end do
  end do

  e1 = e1 / real ( nx, kind = rk ) / real ( ny, kind = rk )

  return
end
subroutine fem2d_l2_error_quadratic ( nx, ny, x, y, u, exact, e2 )

!*****************************************************************************80
!
!! FEM2D_L2_ERROR_LINEAR: L2 error norm of a finite element solution.
!
!  Discussion:
!
!    The finite element method has been used, over a rectangle,
!    involving a grid of NX*NY nodes, with piecewise quadratic elements used 
!    for the basis.
!
!    The finite element coefficients have been computed, and a formula for the
!    exact solution is known.
!
!    This function estimates E2, the L2 norm of the error:
!
!      E2 = Integral ( X, Y ) ( U(X,Y) - EXACT(X,Y) )^2 dX dY
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 June 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NX, NY, the number of nodes in the X 
!    and Y directions.
!
!    Input, real ( kind = rk ) X(NX), Y(NY), the grid coordinates.
!
!    Input, real ( kind = rk ) U(NX,NY), the finite element coefficients.
!
!    Input, real ( kind = rk ), external EXACT(X,Y), returns the value of 
!    the exact solution at the point (X,Y).
!
!    Output, real ( kind = rk ) E2, the estimated L2 norm of the error.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nx
  integer ny
  integer, parameter :: quad_num = 3

  real ( kind = rk ), dimension ( quad_num ) :: abscissa = (/ &
    -0.774596669241483377035853079956D+00, &
     0.000000000000000000000000000000D+00, &
     0.774596669241483377035853079956D+00 /)
  integer cc
  integer e
  real ( kind = rk ) e2
  real ( kind = rk ) eq
  integer ex
  integer ex_num
  real ( kind = rk ), external :: exact
  integer ey
  integer ey_num
  integer il
  integer il2
  integer jl
  integer jl2
  integer k
  integer mm
  integer n
  integer node(9)
  integer qx
  integer qy
  integer s
  real ( kind = rk ) u(nx*ny)
  real ( kind = rk ) uq
  real ( kind = rk ) v(9)
  integer w
  real ( kind = rk ), dimension ( quad_num ) :: weight = (/ &
    0.555555555555555555555555555556D+00, &
    0.888888888888888888888888888889D+00, &
    0.555555555555555555555555555556D+00 /)
  real ( kind = rk ) wq
  real ( kind = rk ) x(nx)
  real ( kind = rk ) xq
  real ( kind = rk ) xx(3)
  real ( kind = rk ) y(ny)
  real ( kind = rk ) yq
  real ( kind = rk ) yy(3)

  e2 = 0.0D+00

  ex_num = ( nx - 1 ) / 2
  ey_num = ( ny - 1 ) / 2

  do ex = 1, ex_num

    w = 2 * ex - 1
    cc = 2 * ex
    e = 2 * ex + 1

    xx(1) = x(w)
    xx(2) = x(cc)
    xx(3) = x(e)

    do ey = 1, ey_num

      s = 2 * ey - 1
      mm = 2 * ey
      n = 2 * ey + 1

      yy(1) = y(s)
      yy(2) = y(mm)
      yy(3) = y(n)
!
!  Node indices
!
!  7  8  9   wn cn en
!  4  5  6   wm cm em
!  1  2  3   ws cs es
!
      node(1) = ( 2 * ey - 2 ) * nx + ( ex - 1 ) * 2 + 1
      node(2) = ( 2 * ey - 2 ) * nx + ( ex - 1 ) * 2 + 2
      node(3) = ( 2 * ey - 2 ) * nx + ( ex - 1 ) * 2 + 3
      node(4) = ( 2 * ey - 1 ) * nx + ( ex - 1 ) * 2 + 1
      node(5) = ( 2 * ey - 1 ) * nx + ( ex - 1 ) * 2 + 2
      node(6) = ( 2 * ey - 1 ) * nx + ( ex - 1 ) * 2 + 3
      node(7) = ( 2 * ey     ) * nx + ( ex - 1 ) * 2 + 1
      node(8) = ( 2 * ey     ) * nx + ( ex - 1 ) * 2 + 2
      node(9) = ( 2 * ey     ) * nx + ( ex - 1 ) * 2 + 3

      do qx = 1, quad_num

        xq = ( ( 1.0D+00 - abscissa(qx) ) * xx(1)  &
             + ( 1.0D+00 + abscissa(qx) ) * xx(3) ) &
               / 2.0D+00

        do qy = 1, quad_num

          yq = ( ( 1.0D+00 - abscissa(qy) ) * yy(1)   &
               + ( 1.0D+00 + abscissa(qy) ) * yy(3) ) &
                 / 2.0D+00

          wq = weight(qx) * ( xx(3) - xx(1) ) / 2.0D+00 &
             * weight(qy) * ( yy(3) - yy(1) ) / 2.0D+00

          uq = 0.0D+00

          v(1:9) = 1.0D+00

          k = 0

          do jl = 1, 3
            do il = 1, 3

              k = k + 1

              do il2 = 1, 3
                if ( il2 /= il ) then
                  v(k) = v(k) * ( xq - xx(il2) ) / ( xx(il) - xx(il2) )
                end if
              end do

              do jl2 = 1, 3
                if ( jl2 /= jl ) then
                  v(k) = v(k) * ( yq - yy(jl2) ) / ( yy(jl) - yy(jl2) )
                end if
              end do

              uq = uq + u(node(k)) * v(k)

            end do
          end do

          eq = exact ( xq, yq )

          e2 = e2 + wq * ( uq - eq ) ** 2

        end do
      end do
    end do
  end do

  e2 = sqrt ( e2 )

  return
end
subroutine r8_fake_use ( x )

!*****************************************************************************80
!
!! r8_fake_use() pretends to use an R8 variable.
!
!  Discussion:
!
!    Some compilers will issue a warning if a variable is unused.
!    Sometimes there's a good reason to include a variable in a program,
!    but not to use it.  Calling this function with that variable as
!    the argument will shut the compiler up.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 April 2020
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) X, the variable to be "used".
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) x

  if ( x /= x ) then
    write ( *, '(a)' ) '  r8_fake_use: variable is NAN.'
  end if

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
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
!  Parameters:
!
!    Input, integer M, the number of rows in A.
!
!    Input, integer N, the number of columns in A.
!
!    Input, real ( kind = rk ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  character ( len = * ) title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_PRINT_SOME prints some of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
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
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = rk ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ILO, JLO, the first row and column to print.
!
!    Input, integer IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: incx = 5
  integer m
  integer n

  real ( kind = rk ) a(m,n)
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

        if ( a(i,j) == real ( int ( a(i,j) ), kind = rk ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine r8mat_solve2 ( n, a, b, x, ierror )

!*****************************************************************************80
!
!! R8MAT_SOLVE2 computes the solution of an N by N linear system.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!    The linear system may be represented as
!
!      A*X = B
!
!    If the linear system is singular, but consistent, then the routine will
!    still produce a solution.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    29 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of equations.
!
!    Input/output, real ( kind = rk ) amat(N,N).
!    On input, A is the coefficient matrix to be inverted.
!    On output, A has been overwritten.
!
!    Input/output, real ( kind = rk ) B(N).
!    On input, B is the right hand side of the system.
!    On output, B has been overwritten.
!
!    Output, real ( kind = rk ) X(N), the solution of the linear system.
!
!    Output, integer IERROR.
!    0, no error detected.
!    1, consistent singularity.
!    2, inconsistent singularity.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) amax
  real ( kind = rk ) b(n)
  integer i
  integer ierror
  integer imax
  integer ipiv(n)
  integer j
  integer k
  real ( kind = rk ) x(n)

  ierror = 0

  ipiv(1:n) = 0
  x(1:n) = 0.0D+00
!
!  Process the matrix.
!
  do k = 1, n
!
!  In column K:
!    Seek the row IMAX with the properties that:
!      IMAX has not already been used as a pivot;
!      amat(IMAX,K) is larger in magnitude than any other candidate.
!
    amax = 0.0D+00
    imax = 0
    do i = 1, n
      if ( ipiv(i) == 0 ) then
        if ( amax < abs ( a(i,k) ) ) then
          imax = i
          amax = abs ( a(i,k) )
        end if
      end if
    end do
!
!  If you found a pivot row IMAX, then,
!    eliminate the K-th entry in all rows that have not been used for pivoting.
!
    if ( imax /= 0 ) then

      ipiv(imax) = k
      a(imax,k+1:n) = a(imax,k+1:n) / a(imax,k)
      b(imax) = b(imax) / a(imax,k)
      a(imax,k) = 1.0D+00

      do i = 1, n

        if ( ipiv(i) == 0 ) then
          a(i,k+1:n) = a(i,k+1:n) - a(i,k) * a(imax,k+1:n)
          b(i) = b(i) - a(i,k) * b(imax)
          a(i,k) = 0.0D+00
        end if

      end do

    end if

  end do
!
!  Now, every row with nonzero IPIV begins with a 1, and
!  all other rows are all zero.  Begin solution.
!
  do j = n, 1, -1

    imax = 0
    do k = 1, n
      if ( ipiv(k) == j ) then
        imax = k
      end if
    end do

    if ( imax == 0 ) then

      x(j) = 0.0D+00

      if ( b(j) == 0.0D+00 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8MAT_SOLVE2 - Warning:'
        write ( *, '(a,i8)' ) '  Consistent singularity, equation = ', j
      else
        ierror = 2
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8MAT_SOLVE2 - Error:'
        write ( *, '(a,i8)' ) '  Inconsistent singularity, equation = ', j
      end if

    else

      x(j) = b(imax)

      do i = 1, n
        if ( i /= imax ) then
          b(i) = b(i) - a(i,j) * x(j)
        end if
      end do

    end if

  end do

  return
end
subroutine r8vec_even ( n, alo, ahi, a )

!*****************************************************************************80
!
!! R8VEC_EVEN returns an R8VEC of evenly spaced values.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    If N is 1, then the midpoint is returned.
!
!    Otherwise, the two endpoints are returned, and N-2 evenly
!    spaced points between them.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of values.
!
!    Input, real ( kind = rk ) ALO, AHI, the low and high values.
!
!    Output, real ( kind = rk ) amat(N), N evenly spaced values.
!    Normally, amat(1) = ALO and amat(N) = AHI.
!    However, if N = 1, then amat(1) = 0.5*(ALO+AHI).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) ahi
  real ( kind = rk ) alo
  integer i

  if ( n == 1 ) then

    a(1) = 0.5D+00 * ( alo + ahi )

  else

    do i = 1, n
      a(i) = ( real ( n - i,     kind = rk ) * alo   &
             + real (     i - 1, kind = rk ) * ahi ) &
             / real ( n     - 1, kind = rk )
    end do

  end if

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
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
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
