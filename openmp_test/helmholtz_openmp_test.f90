program main

!*****************************************************************************80
!
!! helmholtz_openmp_test() solves the Helmholtz equation with OpenMP.
!
!  Discussion:
!
!    The two dimensional region given is:
!
!      -1 <= X <= +1
!      -1 <= Y <= +1
!
!    The region is discretized by a set of M by N nodes:
!
!      P(I,J) = ( X(I), Y(J) )
!
!    where, for 1 <= I <= M, 1 <= J <= N, (FORTRAN convention)
!
!      X(I) = ( 2 * I - M - 1 ) / ( M - 1 )
!      Y(J) = ( 2 * J - N - 1 ) / ( N - 1 )
!
!    The Helmholtz equation for the scalar function U(X,Y) is
!
!      - Uxx(X,Y) -Uyy(X,Y) + ALPHA * U(X,Y) = F(X,Y)
!
!    where ALPHA is a positive constant.  We suppose that Dirichlet
!    boundary conditions are specified, that is, that the value of
!    U(X,Y) is given for all points along the boundary.
!
!    We suppose that the right hand side function F(X,Y) is specified in 
!    such a way that the exact solution is
!
!      U(X,Y) = ( 1 - X**2 ) * ( 1 - Y**2 )
!
!    Using standard finite difference techniques, the second derivatives
!    of U can be approximated by linear combinations of the values
!    of U at neighboring points.  Using this fact, the discretized
!    differential equation becomes a set of linear equations of the form:
!
!      A * U = F
!
!    These linear equations are then solved using a form of the Jacobi 
!    iterative method with a relaxation factor.
!
!    Directives are used in this code to achieve parallelism.
!    All do loops are parallized with default 'static' scheduling.
!
!    Note that the use of the data types "INTEGER ( KIND = 4 )" and
!    "REAL ( kind = rk )" is somewhat nonstandard.  If these declarations
!    are rejected by your compiler, or cause problems in computation,
!    they may be replaced by the simpler "INTEGER" and "DOUBLE PRECISION".
!
!  Modified:
!
!    19 April 2009
!
!  Author:
!
!    Joseph Robicheaux, Sanjiv Shah.
!
  use omp_lib

  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), parameter :: alpha = 0.25D+00
  integer, parameter :: it_max = 100
  integer, parameter :: m = 500
  integer, parameter :: n = 500
  real ( kind = rk ), parameter :: omega = 1.1D+00
  real ( kind = rk ), parameter :: tol = 1.0D-08
  real ( kind = rk ) wtime

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HELMHOLTZ():'
  write ( *, '(a)' ) '  FORTRAN90/OpenMP version'
  write ( *, '(a)' ) '  A program which solves the 2D Helmholtz equation.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This program is being run in parallel.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of processors available = ', omp_get_num_procs ( )
  write ( *, '(a,i8)' ) '  The number of threads available    = ', omp_get_max_threads ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The region is [-1,1] x [-1,1].'
  write ( *, '(a,i8)' ) '  The number of nodes in the X direction is M = ', m
  write ( *, '(a,i8)' ) '  The number of nodes in the Y direction is N = ', n
  write ( *, '(a,i8)' ) '  Number of variables in linear system M * N  = ', m * n
  write ( *, '(a,g14.6)' ) &
    '  The scalar coefficient in the Helmholtz equation is ALPHA = ', alpha
  write ( *, '(a,g14.6)' ) &
    '  The relaxation value is OMEGA = ', omega
  write ( *, '(a,g14.6)' ) &
    '  The error tolerance is TOL = ', tol
  write ( *, '(a,i8)' ) &
    '  The maximum number of Jacobi iterations is IT_MAX = ', it_max
!
!  Call the driver routine.
!
  wtime = omp_get_wtime ( )

  call driver ( m, n, it_max, alpha, omega, tol )

  wtime = omp_get_wtime ( ) - wtime

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Elapsed wall clock time = ', wtime
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HELMHOLTZ():'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
subroutine driver ( m, n, it_max, alpha, omega, tol )

!*****************************************************************************80
!
!! DRIVER allocates arrays and solves the problem.
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    Joseph Robicheaux, Sanjiv Shah.
!
!  Parameters:
!
!    Input, integer M, N, the number of grid points in the 
!    X and Y directions.
!
!    Input, integer IT_MAX, the maximum number of Jacobi 
!    iterations allowed.
!
!    Input, real ( kind = rk ) ALPHA, the scalar coefficient in the
!    Helmholtz equation.
!
!    Input, real ( kind = rk ) OMEGA, the relaxation parameter, which
!    should be strictly between 0 and 2.  For a pure Jacobi method,
!    use OMEGA = 1.
!
!    Input, real ( kind = rk ) TOL, an error tolerance for the linear
!    equation solver.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) alpha
  real ( kind = rk ), allocatable, dimension ( :, : ) :: f
  integer it_max
  real ( kind = rk ) omega
  real ( kind = rk ) tol
  real ( kind = rk ), allocatable, dimension ( :, : ) :: u
!
!  Initialize the data.
!
  allocate ( f(1:m,1:n) )

  call rhs_set ( m, n, alpha, f )
!
!  Solve the Helmholtz equation.
!
  allocate ( u(1:m,1:n) )

!$omp parallel
!$omp workshare

  u(1:m,1:n) = 0.0D+00

!$omp end workshare
!$omp end parallel

  call jacobi ( m, n, alpha, omega, u, f, tol, it_max )
!
!  Determine the error.
!
  call error_check ( m, n, alpha, u, f )

  deallocate ( f )
  deallocate ( u )

  return
end
subroutine error_check ( m, n, alpha, u, f )

!*****************************************************************************80
!
!! ERROR_CHECK determines the error in the numerical solution.
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    Joseph Robicheaux, Sanjiv Shah.
!
!  Parameters:
!
!    Input, integer M, N, the number of grid points in the 
!    X and Y directions.
!
!    Input, real ( kind = rk ) ALPHA, the scalar coefficient in the
!    Helmholtz equation.  ALPHA should be positive.
!
!    Input, real ( kind = rk ) U(M,N), the solution of the Helmholtz equation 
!    at the grid points.
!
!    Input, real ( kind = rk ) F(M,N), values of the right hand side function 
!    for the Helmholtz equation at the grid points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) alpha
  real ( kind = rk ) error_norm
  real ( kind = rk ) f(m,n)
  integer i
  integer j
  real ( kind = rk ) u(m,n)
  real ( kind = rk ) u_exact
  real ( kind = rk ) u_norm
  real ( kind = rk ) u_true
  real ( kind = rk ) u_true_norm
  real ( kind = rk ) x
  real ( kind = rk ) y

!$omp parallel
!$omp workshare

  u_norm = sqrt ( sum ( u(1:m,1:n)**2 ) )

!$omp end workshare
!$omp end parallel

  u_true_norm = 0.0D+00
  error_norm = 0.0D+00

!$omp parallel &
!$omp   shared ( m, n, u ) &
!$omp   private ( i, j, x, y, u_true )

!$omp do reduction ( + : error_norm, u_true_norm )
  do j = 1, n
    do i = 1, m
      x = real ( 2 * i - m - 1, kind = rk ) / real ( m - 1, kind = rk )
      y = real ( 2 * j - n - 1, kind = rk ) / real ( n - 1, kind = rk )
      u_true = u_exact ( x, y )
      error_norm = error_norm + ( u(i,j) - u_true )**2
      u_true_norm = u_true_norm + u_true**2
    end do
  end do

!$omp end do
!$omp end parallel

  error_norm = sqrt ( error_norm )
  u_true_norm = sqrt ( u_true_norm )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Computed U l2 norm :       ', u_norm
  write ( *, '(a,g14.6)' ) '  Computed U_EXACT l2 norm : ', u_true_norm
  write ( *, '(a,g14.6)' ) '  Error l2 norm:             ', error_norm

  return
end
subroutine jacobi ( m, n, alpha, omega, u, f, tol, it_max )

!*****************************************************************************80
!
!! JACOBI applies the Jacobi iterative method to solve the linear system.
!
!  Modified:
!
!    17 November 2007
!
!  Author:
!
!    Joseph Robicheaux, Sanjiv Shah.
!
!  Parameters:
!
!    Input, integer M, N, the number of grid points in the 
!    X and Y directions.
!
!    Input, real ( kind = rk ) ALPHA, the scalar coefficient in the
!    Helmholtz equation.  ALPHA should be positive.
!
!    Input, real ( kind = rk ) OMEGA, the relaxation parameter, which
!    should be strictly between 0 and 2.  For a pure Jacobi method,
!    use OMEGA = 1.
!
!    Input/output, real ( kind = rk ) U(M,N), the solution of the Helmholtz
!    equation at the grid points.
!
!    Input, real ( kind = rk ) F(M,N), values of the right hand side function 
!    for the Helmholtz equation at the grid points.
!
!    Input, real ( kind = rk ) TOL, an error tolerance for the linear
!    equation solver.
!
!    Input, integer IT_MAX, the maximum number of Jacobi 
!    iterations allowed.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) alpha
  real ( kind = rk ) ax
  real ( kind = rk ) ay
  real ( kind = rk ) b
  real ( kind = rk ) dx
  real ( kind = rk ) dy
  real ( kind = rk ) error
  real ( kind = rk ) error_norm
  real ( kind = rk ) f(m,n)
  integer i
  integer it
  integer it_max
  integer j
  real ( kind = rk ) omega
  real ( kind = rk ) tol
  real ( kind = rk ) u(m,n)
  real ( kind = rk ), allocatable, dimension ( :, : ) :: u_old

  allocate ( u_old(1:m,1:n) )
!
!  Initialize the coefficients.
!
  dx = 2.0D+00 / real ( m - 1, kind = rk )
  dy = 2.0D+00 / real ( n - 1, kind = rk )

  ax = -1.0D+00 / dx / dx
  ay = -1.0D+00 / dy / dy
  b  = +2.0D+00 / dx / dx + 2.0D+00 / dy / dy + alpha

  do it = 1, it_max
!
!  Copy new solution into old.
!
!$omp parallel
!$omp workshare

    u_old(1:m,1:n) = u(1:m,1:n)

!$omp end workshare
!$omp end parallel
!
!  Compute stencil, residual, and update.
!
    error_norm = 0.0D+00

!$omp parallel &
!$omp   shared ( ax, ay, b, f, m, n, omega, u, u_old ) &
!$omp   private ( error, i, j )

!$omp do reduction ( + : error_norm )
    do j = 1, n
      do i = 1, m
!
!  Evaluate the residual.
!
        if ( i == 1 .or. i == m .or. j == 1 .or. j == n ) then
          error = u_old(i,j) - f(i,j)
        else
          error = ( ax * ( u_old(i-1,j) + u_old(i+1,j) ) &
                  + ay * ( u_old(i,j-1) + u_old(i,j+1) ) &
            + b * u_old(i,j) - f(i,j) ) / b
        end if
!
!  Update the solution.
!
        u(i,j) = u_old(i,j) - omega * error
!
!  Accumulate the residual error.
!
        error_norm = error_norm + error * error

      end do
    end do

!$omp end do
!$omp end parallel
!
!  Error check.
!
    error_norm = sqrt ( error_norm ) / real ( m * n, kind = rk )

    write ( *, '(2x,i4,a,g14.6)' ) it, '  Residual RMS ', error_norm

    if ( error_norm <= tol ) then
      exit
    end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) 'Total number of iterations ', it

  deallocate ( u_old )

  return
end
subroutine rhs_set ( m, n, alpha, f )

!*****************************************************************************80
!
!! RHS_SET initializes the right hand side.
!
!  Discussion:
!
!    The routine assumes that the exact solution and its second
!    derivatives are given by the routine EXACT.
!
!    The appropriate Dirichlet boundary conditions are determined
!    by getting the value of U returned by EXACT.
!
!    The appropriate right hand side function is determined by
!    having EXACT return the values of U, UXX and UYY, and setting
!
!      F = -UXX - UYY + ALPHA * U
!
!  Modified:
!
!    20 March 2002
!
!  Author:
!
!    Joseph Robicheaux, Sanjiv Shah.
!
!  Parameters:
!
!    Input, integer M, N, the number of grid points in the 
!    X and Y directions.
!
!    Input, real ( kind = rk ) ALPHA, the scalar coefficient in the
!    Helmholtz equation.  ALPHA should be positive.
!
!    Output, real ( kind = rk ) F(M,N), values of the right hand side function 
!    for the Helmholtz equation at the grid points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) alpha
  real ( kind = rk ) f(m,n)
  real ( kind = rk ) f_norm
  integer i
  integer j
  real ( kind = rk ) u_exact
  real ( kind = rk ) u_xx_exact
  real ( kind = rk ) u_yy_exact
  real ( kind = rk ) x
  real ( kind = rk ) y

!$omp parallel
!$omp workshare
  f(1:m,1:n) = 0.0D+00
!$omp end workshare
!$omp end parallel
!
!  Set the boundary conditions.
!

!$omp parallel &
!$omp shared ( alpha, f, m, n ) &
!$omp private ( i, j, x, y )

!$omp do
  do i = 1, m
    j = 1
    y = real ( 2 * j - n - 1, kind = rk ) / real ( n - 1, kind = rk )
    x = real ( 2 * i - m - 1, kind = rk ) / real ( m - 1, kind = rk )
    f(i,j) = u_exact ( x, y )
  end do
!$omp end do

!$omp do
  do i = 1, m
    j = n
    y = real ( 2 * j - n - 1, kind = rk ) / real ( n - 1, kind = rk )
    x = real ( 2 * i - m - 1, kind = rk ) / real ( m - 1, kind = rk )
    f(i,j) = u_exact ( x, y )
  end do
!$omp end do

!$omp do
  do j = 1, n
    i = 1
    x = real ( 2 * i - m - 1, kind = rk ) / real ( m - 1, kind = rk )
    y = real ( 2 * j - n - 1, kind = rk ) / real ( n - 1, kind = rk )
    f(i,j) = u_exact ( x, y )
  end do
!$omp end do

!$omp do
  do j = 1, n
    i = m
    x = real ( 2 * i - m - 1, kind = rk ) / real ( m - 1, kind = rk )
    y = real ( 2 * j - n - 1, kind = rk ) / real ( n - 1, kind = rk )
    f(i,j) = u_exact ( x, y )
  end do
!$omp end do

!$omp do
  do j = 2, n - 1
    do i = 2, m - 1
      x = real ( 2 * i - m - 1, kind = rk ) / real ( m - 1, kind = rk )
      y = real ( 2 * j - n - 1, kind = rk ) / real ( n - 1, kind = rk )
      f(i,j) = - u_xx_exact ( x, y ) - u_yy_exact ( x, y ) + alpha * u_exact ( x, y )
    end do
  end do
!$omp end do

!$omp end parallel

!$omp parallel
!$omp workshare
  f_norm = sqrt ( sum ( f(1:m,1:n)**2 ) )
!$omp end workshare
!$omp end parallel

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Right hand side l2 norm = ', f_norm

  return
end
function u_exact ( x, y )

!*****************************************************************************80
!
!! U_EXACT returns the exact value of the solution.
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, Y, the point at which the values are needed.
!
!    Output, real ( kind = rk ) U_EXACT, the exact value of the solution.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) u_exact
  real ( kind = rk ) x
  real ( kind = rk ) y

  u_exact = ( 1.0D+00 - x**2 ) * ( 1.0D+00 - y**2 )

  return
end
function u_xx_exact ( x, y )

!*****************************************************************************80
!
!! U_XX_EXACT returns the exact value of the second X derivative of the solution.
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, Y, the point at which the values are needed.
!
!    Output, real ( kind = rk ) U_XX_EXACT, the second X derivative of the solution.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) u_xx_exact
  real ( kind = rk ) x
  real ( kind = rk ) y

  u_xx_exact = -2.0D+00 * ( 1.0D+00 + y ) * ( 1.0D+00 - y )

  return
end
function u_yy_exact ( x, y )

!*****************************************************************************80
!
!! U_YY_EXACT returns the exact value of the second Y derivative of the solution.
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, Y, the point at which the values are needed.
!
!    Output, real ( kind = rk ) U_YY_EXACT, the second Y derivative of the solution.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) u_yy_exact
  real ( kind = rk ) x
  real ( kind = rk ) y

  u_yy_exact = -2.0D+00 * ( 1.0D+00 + x ) * ( 1.0D+00 - x )

  return
end
