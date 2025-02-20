subroutine diffusivity_2d_bnt ( dc0, omega, n, x, y, dc )

!*****************************************************************************80
!
!! diffusivity_2d_bnt() evaluates a 2D stochastic diffusivity function.
!
!  Discussion:
!
!    The 2D diffusion equation has the form
!
!      - Del ( DC(X,Y) Del U(X,Y) ) = F(X,Y)
!
!    where DC(X,Y) is a function called the diffusivity.
!
!    In the stochastic version of the problem, the diffusivity function
!    includes the influence of stochastic parameters:
!
!      - Del ( DC(X,Y;OMEGA) Del U(X,Y;OMEGA) ) = F(X,Y).
!
!    In this function, the domain is the rectangle [-1.5,0]x[-0.4,0.8].
!
!    The four stochastic parameters OMEGA(1:4) are assumed to be independent
!    identically distributed random variables with mean value zero and 
!    variance 1.  The distribution is typically taken to be Gaussian or
!    uniform.
!
!    A collocation approach to this problem would then use the roots of
!    Hermite or Legendre polynomials.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 July 2013
!
!  Author:
!
!    John Burkardt

!
!  Reference:
!
!    Ivo Babuska, Fabio Nobile, Raul Tempone,
!    A stochastic collocation method for elliptic partial differential equations
!    with random input data,
!    SIAM Journal on Numerical Analysis,
!    Volume 45, Number 3, 2007, pages 1005-1034.
!
!  Parameters:
!
!    Input, real ( kind = rk8 ) DC0, the constant term in the expansion of the 
!    diffusion coefficient.  Take DC0 = 10.
!
!    Input, real ( kind = rk8 ) OMEGA(4), the stochastic parameters.
!
!    Input, integer N, the number of evaluation points.
!
!    Input, real ( kind = rk8 ) X(N), Y(N), the points where the diffusion 
!    coefficient is to be evaluated.
!
!    Output, real ( kind = rk8 ) DC(N), the value of the diffusion coefficient 
!    at (X,Y).
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  real ( kind = rk8 ) arg(n)
  real ( kind = rk8 ) dc(n)
  real ( kind = rk8 ) dc0
  real ( kind = rk8 ) omega(4)
  real ( kind = rk8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk8 ) x(n)
  real ( kind = rk8 ) y(n)

  arg(1:n) = omega(1) * cos ( pi * x(1:n) ) &
           + omega(2) * sin ( pi * x(1:n) ) &
           + omega(3) * cos ( pi * y(1:n) ) &
           + omega(4) * sin ( pi * y(1:n) )

  arg(1:n)= exp ( - 0.125D+00 ) * arg(1:n)

  dc(1:n) = dc0 + exp ( arg(1:n) )

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! get_unit() returns a free FORTRAN unit number.
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
!    This code is distributed under the GNU LGPL license.
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
subroutine interior ( omega, nx, ny, x, y, f, n, a, rhs )

!*****************************************************************************80
!
!! interior() sets up the matrix and right hand side at interior nodes.
!
!  Discussion:
!
!    Nodes are assigned a single index K, which increases as:
!
!    (NY-1)*NX+1  (NY-1)*NX+2  ...  NY * NX
!           ....         ....  ...    .....
!           NX+1         NX+2  ...   2 * NX
!              1            2  ...       NX
!
!    Therefore, the neighbors of an interior node numbered C are
!
!             C+NY
!              |
!      C-1 --- C --- C+1
!              |
!             C-NY
!
!    If we number rows from bottom I = 1 to top I = NY
!    and columns from left J = 1 to right J = NX, then the relationship
!    between the single index K and the row and column indices I and J is:
!      K = ( I - 1 ) * NX + J
!    and
!      J = 1 + mod ( K - 1, NX )
!      I = 1 + ( K - J ) / NX
!      
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk8 ) OMEGA(4), the stochastic coefficients.
!
!    Input, integer NX, NY, the number of grid points in X and Y.
!
!    Input, real ( kind = rk8 ) X(NX), Y(NY), the coordinates of grid lines.
!
!    Input, real ( kind = rk8 ) function F ( X, Y ), evaluates the heat 
!    source term.
!
!    Input, integer N, the number of nodes.
!
!    Output, real ( kind = rk8 ) A(N,N), the system matrix, with the entries for 
!    the interior nodes filled in.
!
!    Output, real ( kind = rk8 ) RHS(N), the system right hand side, with the 
!    entries for the interior nodes filled in.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n
  integer nx
  integer ny

  real ( kind = rk8 ) a(n,n)
  real ( kind = rk8 ) dc0
  real ( kind = rk8 ) dce(1)
  real ( kind = rk8 ) dcn(1)
  real ( kind = rk8 ) dcs(1)
  real ( kind = rk8 ) dcw(1)
  real ( kind = rk8 ) dx
  real ( kind = rk8 ) dy
  real ( kind = rk8 ), external :: f
  integer ic
  integer in
  integer is
  integer jc
  integer je
  integer jw
  integer kc
  integer ke
  integer kn
  integer ks
  integer kw
  real ( kind = rk8 ) omega(4)
  real ( kind = rk8 ) rhs(n)
  real ( kind = rk8 ) x(nx)
  real ( kind = rk8 ) xce(1)
  real ( kind = rk8 ) xcw(1)
  real ( kind = rk8 ) y(ny)
  real ( kind = rk8 ) ycn(1)
  real ( kind = rk8 ) ycs(1)

  dc0 = 1.0D+00
!
!  For now, assume X and Y are equally spaced.
!
  dx = x(2) - x(1)
  dy = y(2) - y(1)

  do ic = 2, ny - 1
    do jc = 2, nx - 1

      in = ic + 1
      is = ic - 1
      je = jc + 1
      jw = jc - 1

      kc = ( ic - 1 ) * nx + jc
      ke = kc + 1
      kw = kc - 1
      kn = kc + nx
      ks = kc - nx

      xce(1) = 0.5D+00 * ( x(jc) + x(je) )
      call diffusivity_2d_bnt ( dc0, omega, 1, xce(1), y(ic),  dce(1) )
      xcw(1) = 0.5D+00 * ( x(jc) + x(jw) )
      call diffusivity_2d_bnt ( dc0, omega, 1, xcw(1), y(ic),  dcw(1) )
      ycn(1) = 0.5D+00 * ( y(ic) + y(in) )
      call diffusivity_2d_bnt ( dc0, omega, 1, x(jc),  ycn(1), dcn(1) )
      ycs(1) = 0.5D+00 * ( y(ic) + y(is) )
      call diffusivity_2d_bnt ( dc0, omega, 1, x(jc),  ycs(1), dcs(1) )

      a(kc,kc) = ( dce(1) + dcw(1) ) / dx / dx + ( dcn(1) + dcs(1) ) / dy / dy
      a(kc,ke) = - dce(1)            / dx / dx
      a(kc,kw) =          - dcw(1)   / dx / dx
      a(kc,kn) =                                 - dcn(1)            / dy / dy
      a(kc,ks) =                                          - dcs(1)   / dy / dy

      rhs(kc) = f ( x(jc), y(ic) )

    end do
  end do

  return
end
subroutine r8mat_fs ( n, a, b, info )

!*****************************************************************************80
!
!! r8mat_fs() factors and solves a system with one right hand side.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    This routine differs from R8MAT_FSS in two ways:
!    * only one right hand side is allowed;
!    * the input matrix A is not modified.
!
!    This routine uses partial pivoting, but no pivot vector is required.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 January 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input/output, real ( kind = rk8 ) A(N,N).
!    On input, A is the coefficient matrix of the linear system.
!    On output, A is in unit upper triangular form, and
!    represents the U factor of an LU factorization of the
!    original coefficient matrix.
!
!    Input/output, real ( kind = rk8 ) B(N).
!    On input, the right hand side of the linear system.
!    On output, the solution of the linear systems.
!
!    Output, integer INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  real ( kind = rk8 ) a(n,n)
  real ( kind = rk8 ) a2(n,n)
  real ( kind = rk8 ) b(n)
  integer i
  integer info
  integer ipiv
  integer jcol
  real ( kind = rk8 ) piv
  real ( kind = rk8 ) row(n)
  real ( kind = rk8 ) t
  real ( kind = rk8 ) temp

  a2(1:n,1:n) = a(1:n,1:n)

  info = 0

  do jcol = 1, n
!
!  Find the maximum element in column I.
!
    piv = abs ( a2(jcol,jcol) )
    ipiv = jcol
    do i = jcol + 1, n
      if ( piv < abs ( a2(i,jcol) ) ) then
        piv = abs ( a2(i,jcol) )
        ipiv = i
      end if
    end do

    if ( piv == 0.0D+00 ) then
      info = jcol
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8MAT_FS - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      stop
    end if
!
!  Switch rows JCOL and IPIV, and B.
!
    if ( jcol /= ipiv ) then

      row(1:n) = a2(jcol,1:n)
      a2(jcol,1:n) = a2(ipiv,1:n)
      a2(ipiv,1:n) = row(1:n)

      t       = b(jcol)
      b(jcol) = b(ipiv)
      b(ipiv) = t

    end if
!
!  Scale the pivot row.
!
    a2(jcol,jcol+1:n) = a2(jcol,jcol+1:n) / a2(jcol,jcol)
    b(jcol) = b(jcol) / a2(jcol,jcol)
    a2(jcol,jcol) = 1.0D+00
!
!  Use the pivot row to eliminate lower entries in that column.
!
    do i = jcol + 1, n
      if ( a2(i,jcol) /= 0.0D+00 ) then
        temp = - a2(i,jcol)
        a2(i,jcol) = 0.0D+00
        a2(i,jcol+1:n) = a2(i,jcol+1:n) + temp * a2(jcol,jcol+1:n)
        b(i) = b(i) + temp * b(jcol)
      end if
    end do

  end do
!
!  Back solve.
!
  do jcol = n, 2, -1
    b(1:jcol-1) = b(1:jcol-1) - a2(1:jcol-1,jcol) * b(jcol)
  end do

  return
end
subroutine r8vec_linspace ( n, a, b, x )

!*****************************************************************************80
!
!! r8vec_linspace() creates a vector of linearly spaced values.
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
!    This code is distributed under the GNU LGPL license.
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
!    Input, real ( kind = rk8 ) A_FIRST, A_LAST, the first and last entries.
!
!    Output, real ( kind = rk8 ) X(N), a vector of linearly spaced data.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  integer i
  real ( kind = rk8 ) x(n)

  if ( n == 1 ) then

    x(1) = ( a + b ) / 2.0D+00

  else

    do i = 1, n
      x(i) = ( real ( n - i,     kind = rk8 ) * a   &
             + real (     i - 1, kind = rk8 ) * b ) &
             / real ( n     - 1, kind = rk8 )
    end do

  end if

  return
end
subroutine r8vec_mesh_2d ( nx, ny, xvec, yvec, xmat, ymat )

!*****************************************************************************80
!
!! r8vec_mesh_2d() creates a 2D mesh from X and Y vectors.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    NX = 2
!    XVEC = ( 1, 2, 3 )
!    NY = 3
!    YVEC = ( 4, 5 )
!
!    XMAT = (
!      1, 2, 3
!      1, 2, 3 )
!
!    YMAT = (
!      4, 4, 4
!      5, 5, 5 ) 
!    
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 July 2013
!
!  Parameters:
!
!    Input, integer NX, NY, the number of X and Y values.
!
!    Input, real ( kind = rk8 ) XVEC(NX), YVEC(NY), the X and Y coordinate
!    values.
!
!    Output, real ( kind = rk8 ) XMAT(NX,NY), YMAT(NX,NY), the coordinate
!    values of points on an NX by NY mesh.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer nx
  integer ny

  integer j
  real ( kind = rk8 ) xmat(nx,ny)
  real ( kind = rk8 ) xvec(nx)
  real ( kind = rk8 ) ymat(nx,ny)
  real ( kind = rk8 ) yvec(ny)

  do j = 1, ny
    xmat(1:nx,j) = xvec(1:nx)
  end do

  do j = 1, ny
    ymat(1:nx,j) = yvec(j)
  end do

 return
end
subroutine r8vec_normal_01 ( n, x )

!*****************************************************************************80
!
!! r8vec_normal_01() returns a unit pseudonormal R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!    This routine can generate a vector of values on one call.  It
!    has the feature that it should provide the same results
!    in the same order no matter how we break up the task.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of values desired.  
!
!  Output:
!
!    real ( kind = rk8 ) X(N), a sample of the standard normal PDF.
!
!  Local:
!
!    real ( kind = rk8 ) R(N+1), is used to store some uniform
!    random values.  Its dimension is N+1, but really it is only needed
!    to be the smallest even number greater than or equal to N.
!
!    integer X_LO_INDEX, X_HI_INDEX, records the range of entries of
!    X that we need to compute.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  integer m
  real ( kind = rk8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk8 ) r(n+1)
  real ( kind = rk8 ) x(n)
  integer x_hi_index
  integer x_lo_index
!
!  Record the range of X we need to fill in.
!
  x_lo_index = 1
  x_hi_index = n
!
!  Maybe we don't need any more values.
!
  if ( x_hi_index - x_lo_index + 1 == 1 ) then

    call random_number ( harvest = r(1:2) )

    x(x_hi_index) = &
             sqrt ( -2.0D+00 * log ( r(1) ) ) * cos ( 2.0D+00 * pi * r(2) )
!
!  If we require an even number of values, that's easy.
!
  else if ( mod ( x_hi_index - x_lo_index + 1, 2 ) == 0 ) then

    m = ( x_hi_index - x_lo_index + 1 ) / 2

    call random_number ( harvest = r(1:2*m) )
 
    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * cos ( 2.0D+00 * pi * r(2:2*m:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * sin ( 2.0D+00 * pi * r(2:2*m:2) )
!
!  If we require an odd number of values, we generate an even number,
!  and handle the last pair specially, storing one in X(N), and
!  saving the other for later.
!
  else

    x_hi_index = x_hi_index - 1

    m = ( x_hi_index - x_lo_index + 1 ) / 2 + 1

    call random_number ( harvest = r(1:2*m) )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * cos ( 2.0D+00 * pi * r(2:2*m-2:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * sin ( 2.0D+00 * pi * r(2:2*m-2:2) )

    x(n) = sqrt ( -2.0D+00 * log ( r(2*m-1) ) ) &
      * cos ( 2.0D+00 * pi * r(2*m) )

  end if

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! r8vec_print() prints an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, real ( kind = rk8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  real ( kind = rk8 ) a(n)
  integer i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
  end do

  return
end
subroutine stochastic_heat2d ( omega, nx, ny, x, y, f, u )

!*****************************************************************************80
!
!! stochastic_heat2d() solves the steady 2D heat equation.
!
!  Discussion:
!
!    Nodes are assigned a singled index K, which increases as:
!
!    (NY-1)*NX+1  (NY-1)*NX+2  ...  NY * NX
!           ....         ....  ...    .....
!           NX+1         NX+2  ...   2 * NX
!              1            2  ...       NX
!
!    Therefore, the neighbors of an interior node numbered C are
!
!             C+NY
!              |
!      C-1 --- C --- C+1
!              |
!             C-NY
!
!    Nodes on the lower boundary satisfy:
!      1 <= K <= NX
!    Nodes on the upper boundary satisfy:
!      (NY-1)*NX+1 <= K <= NY * NX
!    Nodes on the left boundary satisfy:
!      mod ( K, NX ) = 1
!    Nodes on the right boundary satisfy:
!      mod ( K, NX ) = 0
!
!    If we number rows from bottom I = 1 to top I = NY
!    and columns from left J = 1 to right J = NX, we have
!      K = ( I - 1 ) * NX + J
!    and
!      J = 1 + mod ( K - 1, NX )
!      I = 1 + ( K - J ) / NX
!      
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk8 ) OMEGA(4), the stochastic coefficients.
!
!    Input, integer NX, NY, the number of grid points in X and Y.
!
!    Input, real ( kind = rk8 ) X(NX), Y(NY), the coordinates of grid lines.
!
!    Input, real ( kind = rk8 ) function F ( X, Y ), evaluates the heat 
!    source term.
!
!    Output, real ( kind = rk8 ) U(NX,NY), the approximation to the solution at 
!    the grid points.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer nx
  integer ny

  real ( kind = rk8 ), allocatable :: a(:,:)
  real ( kind = rk8 ), external :: f
  integer info
  integer n
  real ( kind = rk8 ) omega(4)
  real ( kind = rk8 ), allocatable :: rhs(:)
  real ( kind = rk8 ) u(nx,ny)
  real ( kind = rk8 ) x(nx)
  real ( kind = rk8 ) y(ny)
!
!  Set the total number of unknowns.
!
  n = nx * ny
!
!  Set up the matrix and right hand side.
!
  allocate ( a(1:n,1:n) )
  allocate ( rhs(1:n) )
!
!  Define the matrix at interior points.
!
  call interior ( omega, nx, ny, x, y, f, n, a, rhs )
!
!  Handle boundary conditions.
!
  call boundary ( nx, ny, x, y, n, a, rhs )
!
!  Solve the linear system.
!
  call r8mat_fs ( n, a, rhs, info )

  u = reshape ( rhs, (/ nx, ny /) )
!
!  Free memory.
!
  deallocate ( a )
  deallocate ( rhs )

  return
end

