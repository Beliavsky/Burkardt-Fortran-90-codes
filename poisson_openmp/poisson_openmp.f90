program main

!*****************************************************************************80
!
!! poisson_openmp() solves the Poisson equation, with OpenMP for parallel execution.
!
!  Discussion:
!
!    This program runs in parallel, using OpenMP in the SWEEP function.
!
!    The Poisson equation
!
!      - DEL^2 U(x,y) = F(x,y)
!
!    is solved on the unit square [0,1] x [0,1] using a grid of NX by
!    NX evenly spaced points.  The first and last points in each direction
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
!    For convenience in writing the discretized equations, we assume NX = NY.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 December 2011
!
!  Author:
!
!    John Burkardt
!
  use omp_lib

  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nx = 161
  integer, parameter :: ny = 161

  logical converged
  real ( kind = rk ) diff
  real ( kind = rk ) dx
  real ( kind = rk ) dy
  real ( kind = rk ) error
  real ( kind = rk ) f(nx,ny)
  integer i
  integer id
  integer itnew
  integer itold
  integer j
  integer, parameter :: jt_max = 20
  real ( kind = rk ) r8mat_rms
  real ( kind = rk ), parameter :: tolerance = 0.000001D+00
  real ( kind = rk ) u(nx,ny)
  real ( kind = rk ) u_exact
  real ( kind = rk ) u_norm
  real ( kind = rk ) udiff(nx,ny)
  real ( kind = rk ) uexact(nx,ny)
  real ( kind = rk ) unew(nx,ny)
  real ( kind = rk ) unew_norm
  real ( kind = rk ) wtime
  real ( kind = rk ) x
  real ( kind = rk ) y

  dx = 1.0D+00 / real ( nx - 1, kind = rk )
  dy = dx

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POISSON_OPENMP:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  A program for solving the Poisson equation.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Use OpenMP for parallel execution.'
  write ( *, '(a,i4)' ) &
    '  The number of processors available is ', omp_get_num_procs ( )
!$omp parallel
  id = omp_get_thread_num ( )
  if ( id == 0 ) then
    write ( *, '(a,i4)' ) &
      '  The number of threads available is ', omp_get_num_threads ( )
  end if
!$omp end parallel
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

  unew_norm = r8mat_rms ( nx, ny, unew )
!
!  Set up the exact solution UEXACT.
!
  do j = 1, ny 
    y = real ( j - 1, kind = rk ) / real ( ny - 1, kind = rk )
    do i = 1, nx
      x = real ( i - 1, kind = rk ) / real ( nx - 1, kind = rk )
      uexact(i,j) = u_exact ( x, y )
    end do
  end do
  u_norm = r8mat_rms ( nx, ny, uexact )
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
  error = r8mat_rms ( nx, ny, udiff )
  write ( *, '(2x,i6,2x,g14.6,2x,14x,2x,g14.6)' ) 0, unew_norm, error

  wtime = omp_get_wtime ( )

  itnew = 0

  do

    itold = itnew
    itnew = itold + 500
!
!  SWEEP carries out 500 Jacobi steps in parallel before we come
!  back to check for convergence.
!
    call sweep ( nx, ny, dx, dy, f, itold, itnew, u, unew )
!
!  We declare convergence if successive iterates change very little.
!
    u_norm = unew_norm
    unew_norm = r8mat_rms ( nx, ny, unew )
    udiff = unew(1:nx,1:ny) - u(1:nx,1:ny)
    diff = r8mat_rms ( nx, ny, udiff )
    udiff = unew(1:nx,1:ny) - uexact(1:nx,1:ny)
    error = r8mat_rms ( nx, ny, udiff )

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

  wtime = omp_get_wtime ( ) - wtime
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Elapsed seconds = ', wtime
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POISSON_OPENMP:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
function r8mat_rms ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_RMS returns the root mean square of data stored as an R8MAT.
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
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns in A.
!
!    Input, real ( kind = rk ) A(M,N), the data whose RMS is desired.
!
!    Output, real ( kind = rk ) R8MAT_RMS, the root mean square of A.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) r8mat_rms

  r8mat_rms = sqrt ( sum ( a(1:m,1:n)**2 ) / real ( m * n, kind = rk ) )

  return
end
subroutine rhs ( nx, ny, f )

!*****************************************************************************80
!
!! RHS initializes the right hand side "vector".
!
!  Discussion:
!
!    It is convenient for us to set up RHS as a 2D array.  However, each
!    entry of RHS is really the right hand side of a linear system of the
!    form
!
!      A * U = F
!
!    In cases where U(I,J) is a boundary value, then the equation is simply
!
!      U(I,J) = F(i,j)
!
!    and F(I,J) holds the boundary data.
!
!    Otherwise, the equation has the form
!
!      (1/DX^2) * ( U(I+1,J)+U(I-1,J)+U(I,J-1)+U(I,J+1)-4*U(I,J) ) = F(I,J)
!
!    where DX is the spacing and F(I,J) is the value at X(I), Y(J) of
!
!      pi^2 * ( x^2 + y^2 ) * sin ( pi * x * y )
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
!  Parameters:
!
!    Input, integer NX, NY, the X and Y grid dimensions.
!
!    Output, real ( kind = rk ) F(NX,NY), the right hand side data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nx
  integer ny

  real ( kind = rk ) f(nx,ny)
  real ( kind = rk ) fnorm
  integer i
  integer j
  real ( kind = rk ) r8mat_rms
  real ( kind = rk ) u_exact
  real ( kind = rk ) uxxyy_exact
  real ( kind = rk ) x
  real ( kind = rk ) y
!
!  The "boundary" entries of F will store the boundary values of the solution.
!
  x = 0.0D+00
  do j = 1, ny
    y = real ( j - 1, kind = rk ) / real ( ny - 1, kind = rk )
    f(1,j) = u_exact ( x, y )
  end do

  x = 1.0D+00
  do j = 1, ny
    y = real ( j - 1, kind = rk ) / real ( ny - 1, kind = rk )
    f(nx,j) = u_exact ( x, y )
  end do

  y = 0.0D+00
  do i = 1, nx
    x = real ( i - 1, kind = rk ) / real ( nx - 1, kind = rk )
    f(i,1) = u_exact ( x, y )
  end do

  y = 1.0D+00
  do i = 1, nx
    x = real ( i - 1, kind = rk ) / real ( nx - 1, kind = rk )
    f(i,ny) = u_exact ( x, y )
  end do
!
!  The "interior" entries of F store the right hand sides 
!  of the Poisson equation.
!
  do j = 2, ny - 1
    y = real ( j - 1, kind = rk ) / real ( ny - 1, kind = rk )
    do i = 2, nx - 1
      x = real ( i - 1, kind = rk ) / real ( nx - 1, kind = rk )
      f(i,j) = - uxxyy_exact ( x, y )
    end do
  end do

  fnorm = r8mat_rms ( nx, ny, f ) 

  write ( *, '(a,g14.6)' ) '  RMS of F = ', fnorm

  return
end
subroutine sweep ( nx, ny, dx, dy, f, itold, itnew, u, unew )

!*****************************************************************************80
!
!! SWEEP carries out several steps of the Jacobi iteration.
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
!  Parameters:
!
!    Input, integer NX, NY, the X and Y grid dimensions.
!
!    Input, real ( kind = rk ) DX, DY, the spacing between grid points.
!
!    Input, real ( kind = rk ) F(NX,NY), the right hand side data.
!
!    Input, integer ITOLD, the iteration index on input.
!
!    Input, integer ITNEW, the desired iteration index
!    on output.
!
!    Output, real ( kind = rk ) U(NX,NY), the solution estimate on 
!    iteration ITNEW-1.
!
!    Input/output, real ( kind = rk ) UNEW(NX,NY); on input, the solution 
!    estimate on iteration ITOLD.  On output, the solution estimate on 
!    iteration ITNEW.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nx
  integer ny

  real ( kind = rk ) dx
  real ( kind = rk ) dy
  real ( kind = rk ) f(nx,ny)
  integer i
  integer it
  integer itnew
  integer itold
  integer j
  real ( kind = rk ) u(nx,ny)
  real ( kind = rk ) unew(nx,ny)
!
!  Carry out iterations ITOLD+1 through ITNEW.
!
!$omp parallel &
!$omp   shared ( dx, dy, f, itnew, itold, nx, ny, u, unew ) &
!$omp   private ( i, it, j )
!
  do it = itold + 1, itnew

!  Save the current estimate.
!
!$omp do
    do j = 1, ny
      do i = 1, nx
        u(i,j) = unew(i,j)
      end do
    end do
!$omp end do
!
!  Compute a new estimate.
!
!$omp do
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
!$omp end do

  end do
!$omp end parallel

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
  implicit none

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
function u_exact ( x, y )

!*****************************************************************************80
!
!! U_EXACT evaluates the exact solution.
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
!  Parameters:
!
!    Input, real ( kind = rk ) X, Y, the coordinates of a point.
!
!    Output, real ( kind = rk ) U_EXACT, the value of the exact solution 
!    at (X,Y).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) u_exact
  real ( kind = rk ) x
  real ( kind = rk ) y

  u_exact = sin ( pi * x * y )

  return
end
function uxxyy_exact ( x, y )

!*****************************************************************************80
!
!! UXXYY_EXACT evaluates ( d/dx d/dx + d/dy d/dy ) of the exact solution.
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
!  Parameters:
!
!    Input, real ( kind = rk ) X, Y, the coordinates of a point.
!
!    Output, real ( kind = rk ) UXXYY_EXACT, the value of 
!    ( d/dx d/dx + d/dy d/dy ) of the exact solution at (X,Y).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) uxxyy_exact
  real ( kind = rk ) x
  real ( kind = rk ) y

  uxxyy_exact = - pi * pi * ( x * x + y * y ) * sin ( pi * x * y )

  return
end
