subroutine poisson_2d ( nx, ny, rhs, u_exact )

!*****************************************************************************80
!
!! poisson_2d() solves the Poisson equation in a 2D region (a square).
!
!  Discussion:
!
!    The Poisson equation
!
!      - DEL^2 U(x,y) = F(x,y)
!
!    is solved on the unit square [0,1] x [0,1] using a grid of NX by
!    NY evenly spaced points.  The first and last points in each direction
!    are boundary points.
!
!    The boundary conditions and F are set so that the exact solution is
!
!      U(x,y) = sin ( pi * x * y )
!
!    so that
!
!      - DEL^2 U(x,y) = pi^2 * ( x^2 + y^2 ) * sin ( pi * x * y )
!
!    The Jacobi iteration is repeatedly applied until convergence is detected.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 October 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer nx
  integer ny

  logical converged
  real ( kind = rk8 ) diff
  real ( kind = rk8 ) dx
  real ( kind = rk8 ) dy
  real ( kind = rk8 ) error
  real ( kind = rk8 ) f(nx,ny)
  integer i
  integer itnew
  integer itold
  integer j
  integer, parameter :: jt_max = 20
  real ( kind = rk8 ) rms_norm
  real ( kind = rk8 ), parameter :: tolerance = 0.000001D+00
  real ( kind = rk8 ) u(nx,ny)
  real ( kind = rk8 ), external :: u_exact
  real ( kind = rk8 ) u_norm
  real ( kind = rk8 ) udiff(nx,ny)
  real ( kind = rk8 ) uexact(nx,ny)
  real ( kind = rk8 ) unew(nx,ny)
  real ( kind = rk8 ) unew_norm
  real ( kind = rk8 ) x
  real ( kind = rk8 ) y

  dx = 1.0D+00 / real ( nx - 1, kind = rk8 )
  dy = 1.0D+00 / real ( ny - 1, kind = rk8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'poisson_openmp():'
  write ( *, '(a)' ) '  A program for solving the Poisson equation.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  -DEL^2 U = F(X,Y)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  on the rectangle 0 <= X <= 1, 0 <= Y <= 1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  F(X,Y) = pi^2 * ( x^2 + y^2 ) * sin ( pi * x * y )'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of interior X grid points is ', nx
  write ( *, '(a,i8)' ) '  The number of interior Y grid points is ', ny
  write ( *, '(a,f10.4)' ) '  The X grid spacing is ', dx
  write ( *, '(a,f10.4)' ) '  The Y grid spacing is ', dy
!
!  Set the right hand side array F.
!
  call rhs ( nx, ny, f )
!
!  Set the initial solution estimate UNEW.
!  We are "allowed" to pick up the boundary conditions exactly.
!
  unew(1:nx,1:ny) = 0.0D+00
  unew(1,   1:ny) = f(1,   1:ny)
  unew(  nx,1:ny) = f(  nx,1:ny)
  unew(1:nx,1)    = f(1:nx,1)
  unew(1:nx,  ny) = f(1:nx,  ny)

  unew_norm = rms_norm ( nx, ny, unew )
!
!  Set up the exact solution UEXACT.
!
  do j = 1, ny 
    y = real ( j - 1, kind = rk8 ) / real ( ny - 1, kind = rk8 )
    do i = 1, nx
      x = real ( i - 1, kind = rk8 ) / real ( nx - 1, kind = rk8 )
      uexact(i,j) = u_exact ( x, y )
    end do
  end do
  u_norm = rms_norm ( nx, ny, uexact )
  write ( *, '(a,g14.6)' ) '  RMS of exact solution = ', u_norm
!
!  Do the iteration.
!
  converged = .false.

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Step    ||Unew||     ||Unew-U||     ||Unew-Exact||'
  write ( *, '(a)' ) ' '

  udiff = unew(1:nx,1:ny) - uexact(1:nx,1:ny)
  error = rms_norm ( nx, ny, udiff )
  write ( *, '(2x,i6,2x,g14.6,2x,14x,2x,g14.6)' ) 0, unew_norm, error

  itnew = 0

  do

    itold = itnew
    itnew = itold + 500
!
!  Carry out 500 Jacobi steps in parallel before we come
!  back to check for convergence.
!
    call jacobi ( nx, ny, dx, dy, f, itold, itnew, u, unew )
!
!  We declare convergence if successive iterates change very little.
!
    u_norm = unew_norm
    unew_norm = rms_norm ( nx, ny, unew )
    udiff = unew(1:nx,1:ny) - u(1:nx,1:ny)
    diff = rms_norm ( nx, ny, udiff )
    udiff = unew(1:nx,1:ny) - uexact(1:nx,1:ny)
    error = rms_norm ( nx, ny, udiff )

    write ( *, '(2x,i6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      itnew, unew_norm, diff, error

    if ( diff <= tolerance ) then
      converged = .true.
      exit
    end if

  end do

  if ( converged ) then
    write ( *, '(a)' ) '  The iteration has converged.'
  else
    write ( *, '(a)' ) '  The iteration has NOT converged.'
  end if

  return
end
function rms_norm ( m, n, a )

!*****************************************************************************80
!
!! rms_norm() returns the root mean square of data stored as an R8MAT.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, N, the number of rows and columns in A.
!
!    real ( kind = rk8 ) A(M,N), the data whose RMS is desired.
!
!  Output:
!
!    real ( kind = rk8 ) rms_norm, the root mean square of A.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk8 ) a(m,n)
  real ( kind = rk8 ) rms_norm

  rms_norm = sqrt ( sum ( a(1:m,1:n)**2 ) / real ( m * n, kind = rk8 ) )

  return
end
subroutine jacobi ( nx, ny, dx, dy, f, itold, itnew, u, unew )

!*****************************************************************************80
!
!! jacobi() carries out several steps of the Jacobi iteration.
!
!  Discussion:
!
!    Assuming DX = DY, we can approximate
!
!      - ( d/dx d/dx + d/dy d/dy ) U(X,Y) 
!
!    by
!
!      ( U(i-1,j) + U(i+1,j) + U(i,j-1) + U(i,j+1) - 4*U(i,j) ) / dx / dy
!
!    The discretization employed below will not be correct in the general
!    case where DX and DY are not equal.  It's only a little more complicated
!    to allow DX and DY to be different, but we're not going to worry about 
!    that right now.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer NX, NY, the X and Y grid dimensions.
!
!    real ( kind = rk8 ) DX, DY, the spacing between grid points.
!
!    real ( kind = rk8 ) F(NX,NY), the right hand side data.
!
!    integer ITOLD, the iteration index on input.
!
!    integer ITNEW, the desired iteration index on output.
!
!    real ( kind = rk8 ) UNEW(NX,NY); the solution estimate on iteration ITOLD.
!
!  Output:
!
!    real ( kind = rk8 ) U(NX,NY), the solution estimate on iteration ITNEW-1.
!
!    real ( kind = rk8 ) UNEW(NX,NY); the solution estimate on iteration ITNEW.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer nx
  integer ny

  real ( kind = rk8 ) dx
  real ( kind = rk8 ) dy
  real ( kind = rk8 ) f(nx,ny)
  integer i
  integer it
  integer itnew
  integer itold
  integer j
  real ( kind = rk8 ) u(nx,ny)
  real ( kind = rk8 ) unew(nx,ny)
!
!  Carry out iterations ITOLD+1 through ITNEW.
!
  do it = itold + 1, itnew

!  Save the current estimate.
!
    do j = 1, ny
      do i = 1, nx
        u(i,j) = unew(i,j)
      end do
    end do
!
!  Compute a new estimate.
!
    do j = 1, ny
      do i = 1, nx

        if ( i == 1 .or. i == nx .or. j == 1 .or. j == ny ) then
          unew(i,j) = f(i,j)
        else
          unew(i,j) = 0.25D+00 * ( &
            u(i-1,j) + u(i,j+1) + u(i,j-1) + u(i+1,j) + f(i,j) * dx * dy )
        end if

      end do
    end do

  end do

  return
end

