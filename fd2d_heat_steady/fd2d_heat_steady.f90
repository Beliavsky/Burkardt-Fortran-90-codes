subroutine fd2d_heat_steady ( nx, ny, x, y, d, f, u )

!*****************************************************************************80
!
!! FD2D_HEAT_STEADY solves the steady 2D heat equation.
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
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NX, NY, the number of grid points in X and Y.
!
!    Input, real ( kind = rk ) X(NX), Y(NY), the coordinates of grid lines.
!
!    Input, real ( kind = rk ) function D ( X, Y ), evaluates the thermal
!    conductivity.
!
!    Input, real ( kind = rk ) function F ( X, Y ), evaluates the heat 
!    source term.
!
!    Output, real ( kind = rk ) U(NX,NY), the approximation to the solution at 
!    the grid points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nx
  integer ny

  real ( kind = rk ), allocatable :: a(:,:)
  real ( kind = rk ), external :: d
  real ( kind = rk ), external :: f
  integer info
  integer n
  real ( kind = rk ), allocatable :: rhs(:)
  real ( kind = rk ) u(nx,ny)
  real ( kind = rk ) x(nx)
  real ( kind = rk ) y(ny)
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
  call interior ( nx, ny, x, y, d, f, n, a, rhs )
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
subroutine interior ( nx, ny, x, y, d, f, n, a, rhs )

!*****************************************************************************80
!
!! INTERIOR sets up the matrix and right hand side at interior nodes.
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
!    This code is distributed under the MIT license.
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
!    Input, integer NX, NY, the number of grid points in X and Y.
!
!    Input, real ( kind = rk ) X(NX), Y(NY), the coordinates of grid lines.
!
!    Input, real ( kind = rk ) function D ( X, Y ), evaluates the thermal
!    conductivity.
!
!    Input, real ( kind = rk ) function F ( X, Y ), evaluates the heat 
!    source term.
!
!    Input, integer N, the number of nodes.
!
!    Output, real ( kind = rk ) A(N,N), the system matrix, with the entries for 
!    the interior nodes filled in.
!
!    Output, real ( kind = rk ) RHS(N), the system right hand side, with the 
!    entries for the interior nodes filled in.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nx
  integer ny

  real ( kind = rk ) a(n,n)
  real ( kind = rk ), external :: d
  real ( kind = rk ) dc0
  real ( kind = rk ) dce
  real ( kind = rk ) dcn
  real ( kind = rk ) dcs
  real ( kind = rk ) dcw
  real ( kind = rk ) dx
  real ( kind = rk ) dy
  real ( kind = rk ), external :: f
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
  real ( kind = rk ) rhs(n)
  real ( kind = rk ) x(nx)
  real ( kind = rk ) y(ny)

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

      dce = d ( 0.5D+00 * ( x(jc) + x(je) ),             y(ic) )
      dcw = d ( 0.5D+00 * ( x(jc) + x(jw) ),             y(ic) )
      dcn = d (             x(jc),           0.5D+00 * ( y(ic) + y(in) ) )
      dcs = d (             x(jc),           0.5D+00 * ( y(ic) + y(is) ) )

      a(kc,kc) = ( dce + dcw ) / dx / dx + ( dcn + dcs ) / dy / dy
      a(kc,ke) = - dce         / dx / dx
      a(kc,kw) =       - dcw   / dx / dx
      a(kc,kn) =                           - dcn         / dy / dy
      a(kc,ks) =                                 - dcs   / dy / dy

      rhs(kc) = f ( x(jc), y(ic) )

    end do
  end do

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
subroutine r8mat_fs ( n, a, b, info )

!*****************************************************************************80
!
!! R8MAT_FS factors and solves a system with one right hand side.
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
!    This code is distributed under the MIT license. 
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
!    Input/output, real ( kind = rk ) A(N,N).
!    On input, A is the coefficient matrix of the linear system.
!    On output, A is in unit upper triangular form, and
!    represents the U factor of an LU factorization of the
!    original coefficient matrix.
!
!    Input/output, real ( kind = rk ) B(N).
!    On input, the right hand side of the linear system.
!    On output, the solution of the linear systems.
!
!    Output, integer INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) a2(n,n)
  real ( kind = rk ) b(n)
  integer i
  integer info
  integer ipiv
  integer jcol
  real ( kind = rk ) piv
  real ( kind = rk ) row(n)
  real ( kind = rk ) t
  real ( kind = rk ) temp

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
!    Input, real ( kind = rk ) A_FIRST, A_LAST, the first and last entries.
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
subroutine r8vec_mesh_2d ( nx, ny, xvec, yvec, xmat, ymat )

!*****************************************************************************80
!
!! R8VEC_MESH_2D creates a 2D mesh from X and Y vectors.
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
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 July 2013
!
!  Parameters:
!
!    Input, integer NX, NY, the number of X and Y values.
!
!    Input, real ( kind = rk ) XVEC(NX), YVEC(NY), the X and Y coordinate
!    values.
!
!    Output, real ( kind = rk ) XMAT(NX,NY), YMAT(NX,NY), the coordinate
!    values of points on an NX by NY mesh.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nx
  integer ny

  integer j
  real ( kind = rk ) xmat(nx,ny)
  real ( kind = rk ) xvec(nx)
  real ( kind = rk ) ymat(nx,ny)
  real ( kind = rk ) yvec(ny)

  do j = 1, ny
    xmat(1:nx,j) = xvec(1:nx)
  end do

  do j = 1, ny
    ymat(1:nx,j) = yvec(j)
  end do

 return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
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
!    Input, real ( kind = rk ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
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
