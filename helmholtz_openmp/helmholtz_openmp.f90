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
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 November 2024
!
!  Author:
!
!    Original Fortran77 version by Joseph Robicheaux, Sanjiv Shah.
!    This version by John Burkardt.
!
  use omp_lib

  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ), parameter :: alpha = 0.5D+00
  integer, parameter :: it_max = 100
  integer, parameter :: m = 64
  integer, parameter :: n = 64
  real ( kind = rk8 ), parameter :: omega = 0.9D+00
  real ( kind = rk8 ), parameter :: tol = 1.0D-08
  real ( kind = rk8 ) wtime

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'helmholtz_openmp():'
  write ( *, '(a)' ) '  Fortran90/OpenMP version'
  write ( *, '(a)' ) '  A program which solves the 2D Helmholtz equation.'
  write ( *, '(a)' ) '  This program is being run in parallel under OpenMP.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of processors = ', omp_get_num_procs ( )
  write ( *, '(a,i8)' ) '  Number of threads    = ', omp_get_max_threads ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The region is [-1,1] x [-1,1].'
  write ( *, '(a,i8)' ) '  The number of nodes in the X direction is M = ', m
  write ( *, '(a,i8)' ) '  The number of nodes in the Y direction is N = ', n
  write ( *, '(a,i8)' ) '  Number of variables in linear system M * N  = ', m * n
  write ( *, '(a,g14.6)' ) &
    '  The scalar coefficient in the Helmholtz equation is ALPHA = ', alpha
  write ( *, '(a,g14.6)' ) &
    '  The relaxation value is OMEGA = ', omega
  write ( *, '(a,g14.6)' ) '  The error tolerance is TOL = ', tol
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
  write ( *, '(a)' ) 'helmholtz_openmp():'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
subroutine driver ( m, n, it_max, alpha, omega, tol )

!*****************************************************************************80
!
!! driver() allocates arrays and solves the problem.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 November 2024
!
!  Author:
!
!    Original Fortran77 version by Joseph Robicheaux, Sanjiv Shah.
!    This version by John Burkardt.
!
!  Input:
!
!    integer M, N, the number of X and Y grid points.
!
!    integer IT_MAX, the maximum number of Jacobi iterations.
!
!    real ( kind = rk8 ) ALPHA, the scalar coefficient in the
!    Helmholtz equation.
!
!    real ( kind = rk8 ) OMEGA, the relaxation parameter, which
!    should be strictly between 0 and 2.  For a pure Jacobi method,
!    use OMEGA = 1.
!
!    real ( kind = rk8 ) TOL, an error tolerance for the linear
!    equation solver.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk8 ) alpha
  real ( kind = rk8 ), allocatable, dimension ( :, : ) :: f
  integer it_max
  real ( kind = rk8 ) omega
  real ( kind = rk8 ) tol
  real ( kind = rk8 ), allocatable, dimension ( :, : ) :: u
!
!  Initialize the right hand side F.
!
  allocate ( f(1:m,1:n) )

  call rhs_set ( m, n, alpha, f )
!
!  Initialize the solution estimate U.
!
  allocate ( u(1:m,1:n) )
!$omp parallel
!$omp workshare
  u(1:m,1:n) = 0.0D+00
!$omp end workshare
!$omp end parallel
!
!  Solve the Helmholtz equation.
!
  call jacobi ( m, n, alpha, omega, u, f, tol, it_max )
!
!  Determine the error ||U-UEXACT||.
!
  call error_check ( m, n, u )

  deallocate ( f )
  deallocate ( u )

  return
end
subroutine error_check ( m, n, u )

!*****************************************************************************80
!
!! error_check() determines the error in the numerical solution.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 November 2024
!
!  Author:
!
!    Original Fortran77 version by Joseph Robicheaux, Sanjiv Shah.
!    This version by John Burkardt.
!
!  Input:
!
!    integer M, N, the number of X and Y grid points.
!
!    real ( kind = rk8 ) U(M,N), the estimated solution.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk8 ) error_norm
  integer i
  integer j
  real ( kind = rk8 ) u(m,n)
  real ( kind = rk8 ) u_exact
  real ( kind = rk8 ) u_norm
  real ( kind = rk8 ) u_true
  real ( kind = rk8 ) u_true_norm
  real ( kind = rk8 ) x
  real ( kind = rk8 ) y
!
!  1: Compute U norm
!
!$omp parallel
!$omp workshare
  u_norm = sqrt ( sum ( u(1:m,1:n)**2 ) / m / n )
!$omp end workshare
!$omp end parallel
!
!  2: Compute U_True norm.
!
  u_true_norm = 0.0D+00

!$omp parallel &
!$omp   shared ( m, n, u ) &
!$omp   private ( i, j, x, y, u_true )

!$omp do reduction ( + : u_true_norm )
  do j = 1, n
    do i = 1, m
      x = real ( 2 * i - m - 1, kind = rk8 ) / real ( m - 1, kind = rk8 )
      y = real ( 2 * j - n - 1, kind = rk8 ) / real ( n - 1, kind = rk8 )
      u_true = u_exact ( x, y )
      u_true_norm = u_true_norm + u_true**2
    end do
  end do

!$omp end do
!$omp end parallel

  u_true_norm = sqrt ( u_true_norm / m / n )
!
!  3: Rescale U
!
  u = u * u_true_norm / u_norm
  u_norm = u_true_norm
!
!  4: Compute U-Utrue norm
!
  error_norm = 0.0D+00

!$omp parallel &
!$omp   shared ( m, n, u ) &
!$omp   private ( i, j, x, y, u_true )

!$omp do reduction ( + : error_norm )
  do j = 1, n
    do i = 1, m
      x = real ( 2 * i - m - 1, kind = rk8 ) / real ( m - 1, kind = rk8 )
      y = real ( 2 * j - n - 1, kind = rk8 ) / real ( n - 1, kind = rk8 )
      u_true = u_exact ( x, y )
      error_norm = error_norm + ( u(i,j) - u_true )**2
    end do
  end do

!$omp end do
!$omp end parallel

  error_norm = sqrt ( error_norm / m / n )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  RMS U :       ', u_norm
  write ( *, '(a,g14.6)' ) '  RMS UEXACT :  ', u_true_norm
  write ( *, '(a,g14.6)' ) '  RMS U-UEXACT: ', error_norm

  return
end
subroutine jacobi ( m, n, alpha, omega, u, f, tol, it_max )

!*****************************************************************************80
!
!! jacobi() applies the Jacobi iterative method to solve the linear system.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 November 2024
!
!  Author:
!
!    Original Fortran77 version by Joseph Robicheaux, Sanjiv Shah.
!    This version by John Burkardt.
!
!  Input:
!
!    integer M, N, the number of X and Y grid points.
!
!    real ( kind = rk8 ) ALPHA, the scalar coefficient in the
!    Helmholtz equation.  ALPHA should be positive.
!
!    real ( kind = rk8 ) OMEGA, the relaxation parameter, which
!    should be strictly between 0 and 2.  For a pure Jacobi method,
!    use OMEGA = 1.
!
!    real ( kind = rk8 ) U(M,N), the solution estimate.
!
!    real ( kind = rk8 ) F(M,N), values of the right hand side function 
!    for the Helmholtz equation at the grid points.
!
!    real ( kind = rk8 ) TOL, an error tolerance for the linear
!    equation solver.
!
!    integer IT_MAX, the maximum number of Jacobi iterations allowed.
!
!  Output:
!
!    real ( kind = rk8 ) U(M,N), the updated solution estimate.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk8 ) alpha
  real ( kind = rk8 ) ax
  real ( kind = rk8 ) ay
  real ( kind = rk8 ) b
  real ( kind = rk8 ) dx
  real ( kind = rk8 ) dy
  real ( kind = rk8 ) error
  real ( kind = rk8 ) error_norm
  real ( kind = rk8 ) f(m,n)
  integer i
  integer it
  integer it_max
  integer j
  real ( kind = rk8 ) omega
  real ( kind = rk8 ) tol
  real ( kind = rk8 ) u(m,n)
  real ( kind = rk8 ), allocatable, dimension ( :, : ) :: u_old

  allocate ( u_old(1:m,1:n) )
!
!  Initialize the coefficients.
!
  dx = 2.0D+00 / real ( m - 1, kind = rk8 )
  dy = 2.0D+00 / real ( n - 1, kind = rk8 )

  ax = 1.0D+00 / dx / dx
  ay = 1.0D+00 / dy / dy
  b  = -2.0D+00 / dx / dx - 2.0D+00 / dy / dy - alpha

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
    error_norm = sqrt ( error_norm / m / n )

    write ( *, '(2x,i4,a,g14.6)' ) it, '  RMS Residual:', error_norm

    if ( error_norm <= tol ) then
      exit
    end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) 'Number of iterations ', it

  deallocate ( u_old )

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
!    real ( kind = rk8 ) X, the variable to be "used".
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) x

  if ( x /= x ) then
    write ( *, '(a)' ) '  r8_fake_use(): variable is NAN.'
  end if

  return
end
subroutine rhs_set ( m, n, alpha, f )

!*****************************************************************************80
!
!! rhs_set() initializes the right hand side.
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
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 November 2024
!
!  Author:
!
!    Original Fortran77 version by Joseph Robicheaux, Sanjiv Shah.
!    This version by John Burkardt.
!
!  Input:
!
!    integer M, N, the number of X and Y grid points.
!
!    real ( kind = rk8 ) ALPHA, the scalar coefficient in the
!    Helmholtz equation.  ALPHA should be positive.
!
!  Output:
!
!    real ( kind = rk8 ) F(M,N), values of the right hand side function.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk8 ) alpha
  real ( kind = rk8 ) f(m,n)
  real ( kind = rk8 ) f_norm
  integer i
  integer j
  real ( kind = rk8 ) u_exact
  real ( kind = rk8 ) u_xx_exact
  real ( kind = rk8 ) u_yy_exact
  real ( kind = rk8 ) x
  real ( kind = rk8 ) y

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
    y = real ( 2 * j - n - 1, kind = rk8 ) / real ( n - 1, kind = rk8 )
    x = real ( 2 * i - m - 1, kind = rk8 ) / real ( m - 1, kind = rk8 )
    f(i,j) = u_exact ( x, y )
  end do
!$omp end do

!$omp do
  do i = 1, m
    j = n
    y = real ( 2 * j - n - 1, kind = rk8 ) / real ( n - 1, kind = rk8 )
    x = real ( 2 * i - m - 1, kind = rk8 ) / real ( m - 1, kind = rk8 )
    f(i,j) = u_exact ( x, y )
  end do
!$omp end do

!$omp do
  do j = 1, n
    i = 1
    x = real ( 2 * i - m - 1, kind = rk8 ) / real ( m - 1, kind = rk8 )
    y = real ( 2 * j - n - 1, kind = rk8 ) / real ( n - 1, kind = rk8 )
    f(i,j) = u_exact ( x, y )
  end do
!$omp end do

!$omp do
  do j = 1, n
    i = m
    x = real ( 2 * i - m - 1, kind = rk8 ) / real ( m - 1, kind = rk8 )
    y = real ( 2 * j - n - 1, kind = rk8 ) / real ( n - 1, kind = rk8 )
    f(i,j) = u_exact ( x, y )
  end do
!$omp end do

!$omp do
  do j = 2, n - 1
    do i = 2, m - 1
      x = real ( 2 * i - m - 1, kind = rk8 ) / real ( m - 1, kind = rk8 )
      y = real ( 2 * j - n - 1, kind = rk8 ) / real ( n - 1, kind = rk8 )
      f(i,j) = + u_xx_exact ( x, y ) + u_yy_exact ( x, y ) - alpha * u_exact ( x, y )
    end do
  end do
!$omp end do

!$omp end parallel

!$omp parallel
!$omp workshare
  f_norm = sqrt ( sum ( f(1:m,1:n)**2 ) / m / n )
!$omp end workshare
!$omp end parallel

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  RMS right hand side = ', f_norm

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! timestamp() prints the current YMDHMS date as a time stamp.
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
!    15 August 2021
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

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
function u_exact ( x, y )

!*****************************************************************************80
!
!! u_exact() returns the exact value of the solution.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 November 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk8 ) X, Y, the point at which the values are needed.
!
!  Output:
!
!    real ( kind = rk8 ) U_EXACT, the exact value of the solution.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) u_exact
  real ( kind = rk8 ) x
  real ( kind = rk8 ) y

  u_exact = ( 1.0D+00 - x**2 ) * ( 1.0D+00 - y**2 )

  return
end
function u_xx_exact ( x, y )

!*****************************************************************************80
!
!! u_xx_exact() returns the exact value of the second X derivative of the solution.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 November 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk8 ) X, Y, the point at which the values are needed.
!
!  Output:
!
!    real ( kind = rk8 ) U_XX_EXACT, the second X derivative of the solution.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) u_xx_exact
  real ( kind = rk8 ) x
  real ( kind = rk8 ) y

  call r8_fake_use ( x )

  u_xx_exact = -2.0D+00 * ( 1.0D+00 + y ) * ( 1.0D+00 - y )

  return
end
function u_yy_exact ( x, y )

!*****************************************************************************80
!
!! u_yy_exact() returns the exact value of the second Y derivative of the solution.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 November 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk8 ) X, Y, the point at which the values are needed.
!
!  Output:
!
!    real ( kind = rk8 ) U_YY_EXACT, the second Y derivative of the solution.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) u_yy_exact
  real ( kind = rk8 ) x
  real ( kind = rk8 ) y

  call r8_fake_use ( y )

  u_yy_exact = -2.0D+00 * ( 1.0D+00 + x ) * ( 1.0D+00 - x )

  return
end
