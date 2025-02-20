subroutine fem1d_bvp_quadratic ( n, a, c, f, x, u )

!*****************************************************************************80
!
!! FEM1D_BVP_QUADRATIC solves a two point boundary value problem.
!
!  Discussion:
!
!    The finite element method is used, with a mesh of piecewise quadratic
!    elements.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 June 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of nodes.
!
!    Input, function A ( X ), evaluates a(x);
!
!    Input, function C ( X ), evaluates c(x);
!
!    Input, function F ( X ), evaluates f(x);
!
!    Input, real ( kind = rk ) X(N), the mesh points.
!
!    Output, real ( kind = rk ) U(N), the finite element coefficients, which 
!    are also the value of the computed solution at the mesh points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: quad_num = 3

  real ( kind = rk ), external :: a
  real ( kind = rk ) abscissa(quad_num)
  real ( kind = rk ) amat(n,n)
  real ( kind = rk ) axq
  real ( kind = rk ) b(n)
  real ( kind = rk ), external :: c
  real ( kind = rk ) cxq
  integer e
  integer e_num
  real ( kind = rk ), external :: f
  real ( kind = rk ) fxq
  integer ierror
  integer l
  integer m
  integer q
  integer r
  real ( kind = rk ) weight(quad_num)
  real ( kind = rk ) wq
  real ( kind = rk ) u(n)
  real ( kind = rk ) vl
  real ( kind = rk ) vlp
  real ( kind = rk ) vm
  real ( kind = rk ) vmp
  real ( kind = rk ) vr
  real ( kind = rk ) vrp
  real ( kind = rk ) x(n)
  real ( kind = rk ) xl
  real ( kind = rk ) xm
  real ( kind = rk ) xq
  real ( kind = rk ) xr
!
!  Quadrature definitions.
!
  abscissa(1) = -0.774596669241483377035853079956D+00
  abscissa(2) = 0.000000000000000000000000000000D+00
  abscissa(3) = 0.774596669241483377035853079956D+00
  weight(1) = 0.555555555555555555555555555556D+00
  weight(2) = 0.888888888888888888888888888889D+00
  weight(3) = 0.555555555555555555555555555556D+00
!
!  Zero out the matrix and right hand side.
!
  amat(1:n,1:n) = 0.0D+00
  b(1:n) = 0.0D+00
!
!  Integrate over element E.
!
  e_num = ( n - 1 ) / 2

  do e = 1, e_num
!
!  Element E uses nodes
!    L = 2 * E - 1
!    M = 2 * E
!    R = 2 * E + 1
!
    l = 2 * e - 1
    m = 2 * e
    r = 2 * e + 1

    xl = x(l)
    xm = x(m)
    xr = x(r)

    do q = 1, quad_num

      xq = ( ( 1.0 - abscissa(q) ) * xl   &
           + ( 1.0 + abscissa(q) ) * xr ) &
           / 2.0;

      wq = weight(q) * ( xr - xl ) / 2.0

      axq = a ( xq );
      cxq = c ( xq );
      fxq = f ( xq );

      vl = ( ( xq - xm ) / ( xl - xm ) ) &
         * ( ( xq - xr ) / ( xl - xr ) )

      vm = ( ( xq - xl ) / ( xm - xl ) ) &
         * ( ( xq - xr ) / ( xm - xr ) )

      vr = ( ( xq - xl ) / ( xr - xl ) ) &
         * ( ( xq - xm ) / ( xr - xm ) )

      vlp = (     1.0D+00 / ( xl - xm ) ) &
          * ( ( xq - xr ) / ( xl - xr ) ) &
          + ( ( xq - xm ) / ( xl - xm ) ) &
          * (     1.0D+00 / ( xl - xr ) )

      vmp = (     1.0D+00 / ( xm - xl ) ) &
          * ( ( xq - xr ) / ( xm - xr ) ) &
          + ( ( xq - xl ) / ( xm - xl ) ) &
          * (     1.0D+00 / ( xm - xr ) )

      vrp = (     1.0D+00 / ( xr - xl ) ) &
          * ( ( xq - xm ) / ( xr - xm ) ) &
          + ( ( xq - xl ) / ( xr - xl ) ) &
          * (     1.0D+00 / ( xr - xm ) )

      amat(l,l) = amat(l,l) + wq * ( vlp * axq * vlp + vl * cxq * vl )
      amat(l,m) = amat(l,m) + wq * ( vlp * axq * vmp + vl * cxq * vm )
      amat(l,r) = amat(l,r) + wq * ( vlp * axq * vrp + vl * cxq * vr )
      b(l)   = b(l)   + wq * ( vl * fxq )

      amat(m,l) = amat(m,l) + wq * ( vmp * axq * vlp + vm * cxq * vl )
      amat(m,m) = amat(m,m) + wq * ( vmp * axq * vmp + vm * cxq * vm )
      amat(m,r) = amat(m,r) + wq * ( vmp * axq * vrp + vm * cxq * vr )
      b(m) =   b(m)   + wq * ( vm * fxq )

      amat(r,l) = amat(r,l) + wq * ( vrp * axq * vlp + vr * cxq * vl )
      amat(r,m) = amat(r,m) + wq * ( vrp * axq * vmp + vr * cxq * vm )
      amat(r,r) = amat(r,r) + wq * ( vrp * axq * vrp + vr * cxq * vr )
      b(r) =   b(r)   + wq * ( vr * fxq )

    end do

  end do
!
!  Equation 1 is the left boundary condition, U(0.0) = 0.0;
!
  amat(1,1:n) = 0.0D+00
  amat(1,1) = 1.0D+00
  b(1) = 0.0D+00
!
!  Equation N is the right boundary condition, U(1.0) = 0.0;
!
  amat(n,1:n) = 0.0D+00
  amat(n,n) = 1.0D+00
  b(n) = 0.0D+00
!
!  Solve the linear system.
!
  call r8mat_solve2 ( n, amat, b, u, ierror )

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IUNIT, the free unit number.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer ios
  integer iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
subroutine h1s_error_quadratic ( n, x, u, exact_ux, h1s )

!*****************************************************************************80
!
!! H1S_ERROR_QUADRATIC: seminorm error of a finite element solution.
!
!  Discussion:
!
!    We assume the finite element method has been used, over an interval [A,B]
!    involving N nodes, with piecewise quadratic elements used for the basis.
!    The coefficients U(1:N) have been computed, and a formula for the
!    exact derivative is known.
!
!    This function estimates the seminorm of the error:
!
!      SEMINORM = Integral ( A <= X <= B ) ( dU(X)/dx - EXACT_UX(X) )^2 dX
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 June 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of nodes.
!
!    Input, real ( kind = rk ) X(N), the mesh points.
!
!    Input, real ( kind = rk ) U(N), the finite element coefficients.
!
!    Input, function EQ = EXACT_UX ( X ), returns the value of the exact
!    derivative at the point X.
!
!    Output, real ( kind = rk ) H1S the estimated seminorm of 
!    the error.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: quad_num = 3

  real ( kind = rk ) abscissa(quad_num)
  integer e
  integer e_num
  real ( kind = rk ), external :: exact_ux
  real ( kind = rk ) exq
  real ( kind = rk ) h1s
  integer l
  integer m
  integer q
  integer r
  real ( kind = rk ) u(n)
  real ( kind = rk ) uxq
  real ( kind = rk ) vlp
  real ( kind = rk ) vmp
  real ( kind = rk ) vrp
  real ( kind = rk ) weight(quad_num)
  real ( kind = rk ) wq
  real ( kind = rk ) x(n)
  real ( kind = rk ) xl
  real ( kind = rk ) xm
  real ( kind = rk ) xq
  real ( kind = rk ) xr

  h1s = 0.0D+00
!
!  Quadrature definitions.
!
  abscissa(1) = -0.774596669241483377035853079956D+00
  abscissa(2) = 0.000000000000000000000000000000D+00
  abscissa(3) = 0.774596669241483377035853079956D+00
  weight(1) = 0.555555555555555555555555555556D+00
  weight(2) = 0.888888888888888888888888888889D+00
  weight(3) = 0.555555555555555555555555555556D+00
!
!  Integrate over element E.
!
  e_num = ( n - 1 ) / 2

  do e = 1, e_num
!
!  Element E uses nodes
!    L = 2 * E - 1
!    M = 2 * E
!    R = 2 * E + 1
!
    l = 2 * e - 1
    m = 2 * e
    r = 2 * e + 1

    xl = x(l)
    xm = x(m)
    xr = x(r)

    do q = 1, quad_num

      xq = ( ( 1.0D+00 - abscissa(q) ) * xl   &
           + ( 1.0D+00 + abscissa(q) ) * xr ) &
           /   2.0D+00

      wq = weight(q) * ( xr - xl ) / 2.0

      vlp = (     1.0D+00 / ( xl - xm ) ) &
          * ( ( xq - xr ) / ( xl - xr ) ) &
          + ( ( xq - xm ) / ( xl - xm ) ) &
          * (     1.0D+00 / ( xl - xr ) )

      vmp = (     1.0D+00 / ( xm - xl ) ) &
          * ( ( xq - xr ) / ( xm - xr ) ) &
          + ( ( xq - xl ) / ( xm - xl ) ) &
          * (     1.0D+00 / ( xm - xr ) )

      vrp = (     1.0D+00 / ( xr - xl ) ) &
          * ( ( xq - xm ) / ( xr - xm ) ) &
          + ( ( xq - xl ) / ( xr - xl ) ) &
          * (     1.0D+00 / ( xr - xm ) )

      uxq = u(l) * vlp + u(m) * vmp + u(r) * vrp

      exq = exact_ux ( xq )

      h1s = h1s + wq * ( uxq - exq ) ** 2

    end do

  end do

  h1s = sqrt ( h1s )

  return
end
subroutine l1_error ( n, x, u, exact, e1 )

!*****************************************************************************80
!
!! L1_ERROR estimates the l1 error norm of a finite element solution.
!
!  Discussion:
!
!    We assume the finite element method has been used, over an interval [A,B]
!    involving N nodes.
!
!    The coefficients U(1:N) have been computed, and a formula for the
!    exact solution is known.
!
!    This function estimates the little l1 norm of the error:
!      L1_NORM = sum ( 1 <= I <= N ) abs ( U(i) - EXACT(X(i)) )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 June 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of nodes.
!
!    Input, real ( kind = rk ) X(N), the mesh points.
!
!    Input, real ( kind = rk ) U(N), the finite element coefficients.
!
!    Input, function EQ = EXACT ( X ), returns the value of the exact
!    solution at the point X.
!
!    Output, real ( kind = rk ) E1, the estimated L2 norm of the error.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) e1
  real ( kind = rk ), external :: exact
  integer i
  real ( kind = rk ) u(n)
  real ( kind = rk ) x(n)

  e1 = 0.0D+0
  do i = 1, n
    e1 = e1 + abs ( u(i) - exact ( x(i) ) )
  end do

  e1 = e1 / real ( n, kind = rk )

  return
end
subroutine l2_error_quadratic ( n, x, u, exact, e2 )

!*****************************************************************************80
!
!! L2_ERROR_QUADRATIC estimates the L2 error norm of a finite element solution.
!
!  Discussion:
!
!    We assume the finite element method has been used, over an interval [A,B]
!    involving N nodes, with piecewise quadratic elements used for the basis.
!
!    The coefficients U(1:N) have been computed, and a formula for the
!    exact solution is known.
!
!    This function estimates the L2 norm of the error:
!
!      L2_NORM = Integral ( A <= X <= B ) ( U(X) - EXACT(X) )^2 dX
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 June 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of nodes.
!
!    Input, real ( kind = rk ) X(N), the mesh points.
!
!    Input, real ( kind = rk ) U(N), the finite element coefficients.
!
!    Input, function EQ = EXACT ( X ), returns the value of the exact
!    solution at the point X.
!
!    Output, real ( kind = rk ) E2, the estimated L2 norm of the error.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: quad_num = 3

  real ( kind = rk ) abscissa(quad_num)
  integer e
  integer e_num
  real ( kind = rk ) e2
  real ( kind = rk ) eq
  real ( kind = rk ), external :: exact
  integer l
  integer m
  integer q
  integer r
  real ( kind = rk ) u(n)
  real ( kind = rk ) vl
  real ( kind = rk ) vm
  real ( kind = rk ) vr
  real ( kind = rk ) uq
  real ( kind = rk ) weight(quad_num)
  real ( kind = rk ) wq
  real ( kind = rk ) x(n)
  real ( kind = rk ) xl
  real ( kind = rk ) xm
  real ( kind = rk ) xq
  real ( kind = rk ) xr

  e2 = 0.0D+00
!
!  Quadrature definitions.
!
  abscissa(1) = -0.774596669241483377035853079956D+00
  abscissa(2) = 0.000000000000000000000000000000D+00
  abscissa(3) = 0.774596669241483377035853079956D+00
  weight(1) = 0.555555555555555555555555555556D+00
  weight(2) = 0.888888888888888888888888888889D+00
  weight(3) = 0.555555555555555555555555555556D+00
!
!  Integrate over element E.
!
  e_num = ( n - 1 ) / 2

  do e = 1, e_num
!
!  Element E uses nodes
!    L = 2 * E - 1
!    M = 2 * E
!    R = 2 * E + 1
!
    l = 2 * e - 1
    m = 2 * e
    r = 2 * e + 1

    xl = x(l)
    xm = x(m)
    xr = x(r)

    do q = 1, quad_num

      xq = ( ( 1.0D+00 - abscissa(q) ) * xl   &
           + ( 1.0D+00 + abscissa(q) ) * xr ) &
           /   2.0D+00

      wq = weight(q) * ( xr - xl ) / 2.0D+00

      vl = ( ( xq - xm ) / ( xl - xm ) ) &
         * ( ( xq - xr ) / ( xl - xr ) )

      vm = ( ( xq - xl ) / ( xm - xl ) ) &
         * ( ( xq - xr ) / ( xm - xr ) )

      vr = ( ( xq - xl ) / ( xr - xl ) ) &
         * ( ( xq - xm ) / ( xr - xm ) )

      uq = u(l) * vl + u(m) * vm + u(r) * vr
      eq = exact ( xq )

      e2 = e2 + wq * ( uq - eq ) ** 2

    end do

  end do

  e2 = sqrt ( e2 )

  return
end
subroutine max_error_quadratic ( n, x, u, exact, value )

!*****************************************************************************80
!
!! MAX_ERROR_QUADRATIC estimates max error norm of a finite element solution.
!
!  Discussion:
!
!    We assume the finite element method has been used, over an interval [A,B]
!    involving N nodes, with piecewise quadratic elements used for the basis.
!    The coefficients U(1:N) have been computed, and a formula for the
!    exact solution is known.
!
!    This function estimates the max norm of the error:
!
!      MAX_NORM = Integral ( A <= X <= B ) max ( abs ( U(X) - EXACT(X) ) ) dX
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of nodes.
!
!    Input, real ( kind = rk ) X(N), the mesh points.
!
!    Input, real ( kind = rk ) U(N), the finite element coefficients.
!
!    Input, function EQ = EXACT ( X ), returns the value of the exact
!    solution at the point X.
!
!    Output, real ( kind = rk ) VALUE, the estimated max norm of the error.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer e
  integer e_num
  real ( kind = rk ) eq
  real ( kind = rk ), external :: exact
  integer l
  integer m
  integer q
  integer quad_num
  integer r
  real ( kind = rk ) u(n)
  real ( kind = rk ) ul
  real ( kind = rk ) um
  real ( kind = rk ) ur
  real ( kind = rk ) uq
  real ( kind = rk ) value
  real ( kind = rk ) vl
  real ( kind = rk ) vm
  real ( kind = rk ) vr
  real ( kind = rk ) x(n)
  real ( kind = rk ) xl
  real ( kind = rk ) xm
  real ( kind = rk ) xq
  real ( kind = rk ) xr

  quad_num = 8
  value = 0.0D+00

  e_num = ( n - 1 ) / 2

  do e = 1, e_num

    l = 2 * e - 1
    xl = x(l)
    ul = u(l)

    m = 2 * e
    xm = x(m)
    um = u(m)

    r = 2 * e + 1
    xr = x(r)
    ur = u(r)

    do q = 0, quad_num - 1

      xq = ( real ( quad_num - q, kind = rk ) * xl   &
           + real (            q, kind = rk ) * xr ) &
         /   real ( quad_num, kind = rk )

      vl = ( ( xq - xm ) / ( xl - xm ) ) &
         * ( ( xq - xr ) / ( xl - xr ) )

      vm = ( ( xq - xl ) / ( xm - xl ) ) &
         * ( ( xq - xr ) / ( xm - xr ) )

      vr = ( ( xq - xl ) / ( xr - xl ) ) &
         * ( ( xq - xm ) / ( xr - xm ) )

      uq = u(l) * vl + u(m) * vm + u(r) * vr

      eq = exact ( xq )

      value = max ( value, abs ( uq - eq ) )

    end do

  end do
!
!  For completeness, check last node.
!
  xq = x(n)
  uq = u(n)
  eq = exact ( xq )

  value = max ( value, abs ( uq - eq ) )
!
!  Integral approximation requires multiplication by interval length.
!
  value = value * ( x(n) - x(1) )

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
!    Input/output, real ( kind = rk ) A(N,N).
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
!      A(IMAX,K) is larger in magnitude than any other candidate.
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
subroutine r8vec_linspace ( n, a, b, x )

!*****************************************************************************80
!
!! R8VEC_LINSPACE creates a vector of linearly spaced values.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    4 points evenly spaced between 0 and 12 will yield 0, 4, 8, 12.
!
!    In other words, the interval is divided into N-1 even subintervals,
!    and the endpoints of intervals are used as the points.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real ( kind = rk ) A, B, the first and last entries.
!
!    Output, real ( kind = rk ) X(N), a vector of linearly spaced data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a
  real ( kind = rk ) b
  integer i
  real ( kind = rk ) x(n)

  if ( n == 1 ) then

    x(1) = ( a + b ) / 2.0D+00

  else

    do i = 1, n
      x(i) = ( real ( n - i,     kind = rk ) * a   &
             + real (     i - 1, kind = rk ) * b ) &
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
