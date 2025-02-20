subroutine hermite_compute ( order, xtab, weight )

!*****************************************************************************80
!
!! hermite_compute() computes a Gauss-Hermite quadrature rule.
!
!  Discussion:
!
!    The abscissas are the zeros of the N-th order Hermite polynomial.
!
!    The integration interval is ( -oo, +oo ).
!
!    The weight function is w(x) = exp ( - x * x ).
!
!    The integral to approximate:
!
!      integral ( -oo < X < +oo ) exp ( - X * X ) * F(X) dX
!
!    The quadrature rule:
!
!      sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input, integer ORDER, the order of the formula.
!
!    Output, real ( kind = rk ) XTAB(ORDER), the abscissas.
!
!    Output, real ( kind = rk ) WEIGHT(ORDER), the weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer order

  real ( kind = rk ) cc
  real ( kind = rk ) dp2
  integer i
  real ( kind = rk ) p1
  real ( kind = rk ) s
  real ( kind = rk ) temp
  real ( kind = rk ) weight(order)
  real ( kind = rk ) x
  real ( kind = rk ) xtab(order)

  cc = 1.7724538509D+00 * gamma ( real ( order, kind = rk ) ) &
    / ( 2.0D+00 ** ( order - 1 ) )

  s = ( 2.0D+00 * real ( order, kind = rk ) + 1.0D+00 ) ** ( 1.0D+00 / 6.0D+00 )

  do i = 1, ( order + 1 ) / 2

    if ( i == 1 ) then

      x = s ** 3 - 1.85575D+00 / s

    else if ( i == 2 ) then

      x = x - 1.14D+00 * ( ( real ( order, kind = rk ) ) ** 0.426D+00 ) / x

    else if ( i == 3 ) then

      x = 1.86D+00 * x - 0.86D+00 * xtab(1)

    else if ( i == 4 ) then

      x = 1.91D+00 * x - 0.91D+00 * xtab(2)

    else

      x = 2.0D+00 * x - xtab(i-2)

    end if

    call hermite_root ( x, order, dp2, p1 )

    xtab(i) = x
    weight(i) = ( cc / dp2 ) / p1

    xtab(order-i+1) = - x
    weight(order-i+1) = weight(i)

  end do
!
!  Reverse the order.
!
  do i = 1, order / 2
    temp            = xtab(i)
    xtab(i)         = xtab(order+1-i)
    xtab(order+1-i) = temp
  end do

  return
end
subroutine hermite_integral ( n, value )

!*****************************************************************************80
!
!! HERMITE_INTEGRAL returns the value of a Hermite polynomial integral.
!
!  Discussion:
!
!    H(n) = Integral ( -oo < x < +oo ) x^n exp(-x^2) dx
!
!    H(n) is 0 for n odd.
!
!    H(n) = (n-1)!! * sqrt(pi) / 2^(n/2) for n even.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the integral.  
!    0 <= N.
!
!    Output, real ( kind = rk ) VALUE, the value of the integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i4_factorial2
  integer n
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) value

  if ( n < 0 ) then

    value = - huge ( value )

  else if ( mod ( n, 2 ) == 1 ) then

    value = 0.0D+00

  else

    value = real ( i4_factorial2 ( n - 1 ), kind = rk ) * sqrt ( r8_pi ) &
      / 2.0D+00 ** ( n / 2 )

  end if

  return
end
subroutine hermite_recur ( p2, dp2, p1, x, order )

!*****************************************************************************80
!
!! HERMITE_RECUR finds the value and derivative of a Hermite polynomial.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 September 1998
!
!  Author:
!
!    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Output, real ( kind = rk ) P2, the value of H(ORDER)(X).
!
!    Output, real ( kind = rk ) DP2, the value of H'(ORDER)(X).
!
!    Output, real ( kind = rk ) P1, the value of H(ORDER-1)(X).
!
!    Input, real ( kind = rk ) X, the point at which polynomials are evaluated.
!
!    Input, integer ORDER, the order of the polynomial 
!    to be computed.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  real ( kind = rk ) dp0
  real ( kind = rk ) dp1
  real ( kind = rk ) dp2
  integer order
  real ( kind = rk ) p0
  real ( kind = rk ) p1
  real ( kind = rk ) p2
  real ( kind = rk ) x

  p1 = 1.0D+00
  dp1 = 0.0D+00

  p2 = x
  dp2 = 1.0D+00

  do i = 2, order

    p0 = p1
    dp0 = dp1

    p1 = p2
    dp1 = dp2

    p2  = x * p1 - 0.5D+00 * ( real ( i, kind = rk ) - 1.0D+00 ) * p0
    dp2 = x * dp1 + p1 - 0.5D+00 * ( real ( i, kind = rk ) - 1.0D+00 ) * dp0

  end do

  return
end
subroutine hermite_root ( x, order, dp2, p1 )

!*****************************************************************************80
!
!! HERMITE_ROOT improves an approximate root of a Hermite polynomial.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 September 1998
!
!  Author:
!
!    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input/output, real ( kind = rk ) X, the approximate root, which
!    should be improved on output.
!
!    Input, integer ORDER, the order of the Hermite polynomial.
!
!    Output, real ( kind = rk ) DP2, the value of H'(ORDER)(X).
!
!    Output, real ( kind = rk ) P1, the value of H(ORDER-1)(X).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) d
  real ( kind = rk ) dp2
  real ( kind = rk ), parameter :: eps = 1.0D-12
  integer order
  real ( kind = rk ) p1
  real ( kind = rk ) p2
  integer step
  integer, parameter :: step_max = 10
  real ( kind = rk ) x

  do step = 1, step_max

    call hermite_recur ( p2, dp2, p1, x, order )

    d = p2 / dp2
    x = x - d

    if ( abs ( d ) <= eps * ( abs ( x ) + 1.0D+00 ) ) then
      return
    end if

  end do

  return
end
function i4_factorial2 ( n )

!*****************************************************************************80
!
!! I4_FACTORIAL2 computes the double factorial function.
!
!  Discussion:
!
!    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
!                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
!
!  Example:
!
!     N    N!!
!
!     0     1
!     1     1
!     2     2
!     3     3
!     4     8
!     5    15
!     6    48
!     7   105
!     8   384
!     9   945
!    10  3840
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the argument of the double factorial 
!    function.  If N is less than 1, I4_FACTORIAL2 is returned as 1.
!
!    Output, integer I4_FACTORIAL2, the value of the double
!    factorial.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i4_factorial2
  integer n
  integer n_copy

  if ( n < 1 ) then
    i4_factorial2 = 1
    return
  end if

  n_copy = n
  i4_factorial2 = 1

  do while ( 1 < n_copy )
    i4_factorial2 = i4_factorial2 * n_copy
    n_copy = n_copy - 2
  end do

  return
end
subroutine p00_exact ( problem, exact )

!*****************************************************************************80
!
!! P00_EXACT returns the exact integral for any problem.
!
!  Discussion:
!
!    This routine provides a "generic" interface to the exact integral
!    routines for the various problems, and allows a problem to be called
!    by index (PROBLEM) rather than by name.
!
!    In most cases, the "exact" value of the integral is not given;
!    instead a "respectable" approximation is available.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer PROBLEM, the index of the problem.
!
!    Output, real ( kind = rk ) EXACT, the exact value of the integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) exact
  integer problem

  if ( problem == 1 ) then
    call p01_exact ( exact )
  else if ( problem == 2 ) then
    call p02_exact ( exact )
  else if ( problem == 3 ) then
    call p03_exact ( exact )
  else if ( problem == 4 ) then
    call p04_exact ( exact )
  else if ( problem == 5 ) then
    call p05_exact ( exact )
  else if ( problem == 6 ) then
    call p06_exact ( exact )
  else if ( problem == 7 ) then
    call p07_exact ( exact )
  else if ( problem == 8 ) then
    call p08_exact ( exact )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_EXACT - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal problem number = ', problem
    stop
  end if

  return
end
subroutine p00_fun ( problem, option, n, x, f )

!*****************************************************************************80
!
!! P00_FUN evaluates the integrand for any problem.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer PROBLEM, the index of the problem.
!
!    Input, integer OPTION:
!    0, integrand is f(x).
!    1, integrand is exp(-x*x) * f(x);
!    2, integrand is exp(-x*x/2) * f(x);
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) X(N), the evaluation points.
!
!    Output, real ( kind = rk ) F(N), the function values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) f(n)
  integer option
  integer problem
  real ( kind = rk ) x(n)

  if ( problem == 1 ) then
    call p01_fun ( option, n, x, f )
  else if ( problem == 2 ) then
    call p02_fun ( option, n, x, f )
  else if ( problem == 3 ) then
    call p03_fun ( option, n, x, f )
  else if ( problem == 4 ) then
    call p04_fun ( option, n, x, f )
  else if ( problem == 5 ) then
    call p05_fun ( option, n, x, f )
  else if ( problem == 6 ) then
    call p06_fun ( option, n, x, f )
  else if ( problem == 7 ) then
    call p07_fun ( option, n, x, f )
  else if ( problem == 8 ) then
    call p08_fun ( option, n, x, f )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_FUN - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal problem number = ', problem
    stop
  end if

  return
end
subroutine p00_gauss_hermite ( problem, order, result )

!*****************************************************************************80
!
!! P00_GAUSS_HERMITE applies a Gauss-Hermite quadrature rule.
!
!  Discussion:
!
!    The Gauss-Hermite rule that is applied has been defined for
!      I(f) = integral ( -oo < x < +oo ) f(x) exp(-x^2) dx
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer PROBLEM, the index of the problem.
!
!    Input, integer ORDER, the order of the Gauss-Laguerre rule
!    to apply.
!
!    Output, real ( kind = rk ) RESULT, the approximate integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable, dimension ( : ) :: f_vec
  integer option
  integer order
  integer problem
  real ( kind = rk ) result
  real ( kind = rk ), allocatable, dimension ( : ) :: weight
  real ( kind = rk ), allocatable, dimension ( : ) :: xtab

  allocate ( f_vec(1:order) )
  allocate ( weight(1:order) )
  allocate ( xtab(1:order) )

  call hermite_compute ( order, xtab, weight )

  option = 1
  call p00_fun ( problem, option, order, xtab, f_vec )

  result = dot_product ( weight(1:order), f_vec(1:order) )

  deallocate ( f_vec )
  deallocate ( weight )
  deallocate ( xtab )

  return
end
subroutine p00_monte_carlo ( problem, order, result )

!*****************************************************************************80
!
!! P00_MONTE_CARLO applies a Monte Carlo procedure to Hermite integrals.
!
!  Discussion:
!
!    We wish to estimate the integral:
!
!      I(f) = integral ( -oo < x < +oo ) f(x) exp ( - x * x ) dx
!
!    We do this by a Monte Carlo sampling procedure, in which 
!    we select N points X(1:N) from a standard normal distribution,
!    and estimate
!
!      Q(f) = sum ( 1 <= I <= N ) f(x(i)) / sqrt ( pi )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 May 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer PROBLEM, the index of the problem.
!
!    Input, integer ORDER, the order of the Gauss-Laguerre rule
!    to apply.
!
!    Output, real ( kind = rk ) RESULT, the approximate integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable, dimension ( : ) :: f_vec
  integer option
  integer order
  integer problem
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) result
  real ( kind = rk ) weight
  real ( kind = rk ), allocatable, dimension ( : ) :: x_vec

  allocate ( f_vec(1:order) )
  allocate ( x_vec(1:order) )

  call r8vec_normal_01 ( order, x_vec )

  option = 2
  call p00_fun ( problem, option, order, x_vec, f_vec )

  weight = real ( order, kind = rk ) / sqrt ( r8_pi ) / sqrt ( 2.0D+00 )

  result = sum ( f_vec(1:order) ) / weight

  deallocate ( f_vec )
  deallocate ( x_vec )

  return
end
subroutine p00_problem_num ( problem_num )

!*****************************************************************************80
!
!! P00_PROBLEM_NUM returns the number of test integration problems.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer PROBLEM_NUM, the number of test problems.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer problem_num

  problem_num = 8

  return
end
subroutine p00_title ( problem, title )

!*****************************************************************************80
!
!! P00_TITLE returns the title for any problem.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer PROBLEM, the index of the problem.
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer problem
  character ( len = * ) title

  if ( problem == 1 ) then
    call p01_title ( title )
  else if ( problem == 2 ) then
    call p02_title ( title )
  else if ( problem == 3 ) then
    call p03_title ( title )
  else if ( problem == 4 ) then
    call p04_title ( title )
  else if ( problem == 5 ) then
    call p05_title ( title )
  else if ( problem == 6 ) then
    call p06_title ( title )
  else if ( problem == 7 ) then
    call p07_title ( title )
  else if ( problem == 8 ) then
    call p08_title ( title )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_TITLE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal problem number = ', problem
    stop
  end if

  return
end
subroutine p00_turing ( problem, h, tol, n, result )

!*****************************************************************************80
!
!! P00_TURING applies the Turing quadrature rule.
!
!  Discussion:
!
!    We consider the approximation:
!
!      Integral ( -oo < x < +oo ) f(x) dx
!
!      = h * Sum ( -oo < i < +oo ) f(nh) + error term
!
!    Given H and a tolerance TOL, we start summing at I = 0, and
!    adding one more term in the positive and negative I directions,
!    until the absolute value of the next two terms being added 
!    is less than TOL.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan Turing,
!    A Method for the Calculation of the Zeta Function,
!    Proceedings of the London Mathematical Society,
!    Volume 48, 1943, pages 180-197.
!
!  Parameters:
!
!    Input, integer PROBLEM, the index of the problem.
!
!    Input, real ( kind = rk ) H, the spacing to use.
!
!    Input, real ( kind = rk ) TOL, the tolerance.  
!
!    Output, integer N, the number of pairs of steps taken.
!    The actual number of function evaluations is 2*N+1.
!
!    Output, real ( kind = rk ) RESULT, the approximate integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) f_vec(2)
  real ( kind = rk ) h
  integer n
  integer, parameter :: n_too_many = 100000
  integer option
  integer order
  integer problem
  real ( kind = rk ) result
  real ( kind = rk ) tol
  real ( kind = rk ) xtab(2)

  option = 0
  n = 0

  result = 0.0D+00
  order = 1
  xtab(1) = 0.0D+00
  call p00_fun ( problem, option, order, xtab, f_vec )
  result = result + h * f_vec(1)

  do

    n = n + 1

    xtab(1) =   real ( n, kind = rk ) * h
    xtab(2) = - real ( n, kind = rk ) * h

    order = 2
    call p00_fun ( problem, option, order, xtab, f_vec )

    result = result + h * ( f_vec(1) + f_vec(2) )
!
!  Just do a simple-minded absolute error tolerance check to start with.
!
    if ( abs ( f_vec(1) ) + abs ( f_vec(2) ) <= tol ) then
      exit
    end if
!
!  Just in case things go crazy.
!
    if ( n_too_many <= n ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P00_TURING - Warning!'
      write ( *, '(a,i8)' ) &
        '  Number of steps exceeded N_TOO_MANY = ', n_too_many
      exit
    end if

  end do

  return
end
subroutine p01_exact ( exact )

!*****************************************************************************80
!
!! P01_EXACT returns the exact integral for problem 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) EXACT, the value of the integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) exact
  real ( kind = rk ), parameter :: omega = 1.0D+00
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00

  exact = sqrt ( r8_pi ) * exp ( - omega * omega )

  return
end
subroutine p01_fun ( option, n, x, f )

!*****************************************************************************80
!
!! P01_FUN evaluates the integrand for problem 1.
!
!  Discussion:
!
!    Squire gives exact value as sqrt(pi) * exp(-w*w).
!
!    Integral ( -oo < x < +oo ) exp(-x*x) cos(2*w*x) dx
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Squire,
!    Comparison of Gauss-Hermite and Midpoint Quadrature with Application
!    to the Voigt Function,
!    in Numerical Integration: 
!    Recent Developments, Software and Applications,
!    edited by Patrick Keast, Graeme Fairweather,
!    Reidel, 1987, pages 337-340,
!    ISBN: 9027725144,
!    LC: QA299.3.N38.
!
!  Parameters:
!
!    Input, integer OPTION:
!    0, integrand is f(x).
!    1, integrand is exp(-x*x) * f(x);
!    2, integrand is exp(-x*x/2) * f(x);
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) X(N), the evaluation points.
!
!    Output, real ( kind = rk ) F(N), the function values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) f(n)
  real ( kind = rk ), parameter :: omega = 1.0D+00
  integer option
  real ( kind = rk ) x(n)

  if ( option == 0 ) then
    f(1:n) = cos ( 2.0D+00 * omega * x(1:n) ) &
      * exp ( - x(1:n) * x(1:n) )
  else if ( option == 1 ) then
    f(1:n) = cos ( 2.0D+00 * omega * x(1:n) )
  else if ( option == 2 ) then
    f(1:n) = cos ( 2.0D+00 * omega * x(1:n) ) &
      * exp ( - 0.5D+00 * x(1:n) * x(1:n) )
  end if

  return
end
subroutine p01_title ( title )

!*****************************************************************************80
!
!! P01_TITLE returns the title for problem 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = * ) title

  title = 'cos(2*omega*x)'

  return
end
subroutine p02_exact ( exact )

!*****************************************************************************80
!
!! P02_EXACT returns the exact integral for problem 2.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) EXACT, the value of the integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) exact
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00

  exact = sqrt ( r8_pi )

  return
end
subroutine p02_fun ( option, n, x, f )

!*****************************************************************************80
!
!! P02_FUN evaluates the integrand for problem 2.
!
!  Discussion:
!
!    The exact value is sqrt(pi).
!
!    Integral ( -oo < x < +oo ) exp(-x*x) dx
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer OPTION:
!    0, integrand is f(x).
!    1, integrand is exp(-x*x) * f(x);
!    2, integrand is exp(-x*x/2) * f(x);
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) X(N), the evaluation points.
!
!    Output, real ( kind = rk ) F(N), the function values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) f(n)
  integer option
  real ( kind = rk ) x(n)

  if ( option == 0 ) then
    f(1:n) = 1.0D+00 &
      * exp ( - x(1:n) * x(1:n) ) 
  else if ( option == 1 ) then
    f(1:n) = 1.0D+00
  else if ( option == 2 ) then
    f(1:n) = 1.0D+00 &
      * exp ( - 0.5D+00 * x(1:n) * x(1:n) )
  end if

  return
end
subroutine p02_title ( title )

!*****************************************************************************80
!
!! P02_TITLE returns the title for problem 2.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = * ) title

  title = 'exp(-x*x)'

  return
end
subroutine p03_exact ( exact )

!*****************************************************************************80
!
!! P03_EXACT returns the exact integral for problem 3.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) EXACT, the value of the integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) exact
  real ( kind = rk ), parameter :: p = 1.0D+00
  real ( kind = rk ), parameter :: q = 3.0D+00
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00

  exact = r8_pi / ( q * sin ( r8_pi * p / q ) )

  return
end
subroutine p03_fun ( option, n, x, f )

!*****************************************************************************80
!
!! P03_FUN evaluates the integrand for problem 3.
!
!  Discussion:
!
!    The exact value is pi / (q*sin(pi*p/q) ), assuming 0 < p < q.
!
!    Integral ( -oo < x < +oo ) exp(-px) / ( 1 + exp ( -qx) ) dx
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer OPTION:
!    0, integrand is f(x).
!    1, integrand is exp(-x*x) * f(x);
!    2, integrand is exp(-x*x/2) * f(x);
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) X(N), the evaluation points.
!
!    Output, real ( kind = rk ) F(N), the function values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) f(n)
  integer option
  real ( kind = rk ), parameter :: p = 1.0D+00
  real ( kind = rk ), parameter :: q = 3.0D+00
  real ( kind = rk ) x(n)

  if ( option == 0 ) then
    f(1:n) = exp ( - p * x(1:n) ) / ( 1.0D+00 + exp ( -q * x(1:n) ) )
  else if ( option == 1 ) then
    f(1:n) = exp ( - p * x(1:n) ) / ( 1.0D+00 + exp ( -q * x(1:n) ) ) &
      * exp ( + x(1:n) * x(1:n) )
  else if ( option == 2 ) then
    f(1:n) = exp ( - p * x(1:n) ) / ( 1.0D+00 + exp ( -q * x(1:n) ) ) &
      * exp ( + 0.5D+00 * x(1:n) * x(1:n) )
  end if

  return
end
subroutine p03_title ( title )

!*****************************************************************************80
!
!! P03_TITLE returns the title for problem 3.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = * ) title

  title = 'exp(-px) / ( 1 + exp(-qx) )'

  return
end
subroutine p04_exact ( exact )

!*****************************************************************************80
!
!! P04_EXACT returns the exact integral for problem 4.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) EXACT, the value of the integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) exact
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00

  exact = sqrt ( r8_pi / 2.0D+00 )

  return
end
subroutine p04_fun ( option, n, x, f )

!*****************************************************************************80
!
!! P04_FUN evaluates the integrand for problem 4.
!
!  Discussion:
!
!    The exact value is sqrt ( pi / 2 )
!
!    Integral ( -oo < x < +oo ) sin ( x*x ) dx
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer OPTION:
!    0, integrand is f(x).
!    1, integrand is exp(-x*x) * f(x);
!    2, integrand is exp(-x*x/2) * f(x);
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) X(N), the evaluation points.
!
!    Output, real ( kind = rk ) F(N), the function values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) f(n)
  integer option
  real ( kind = rk ) x(n)

  if ( option == 0 ) then
    f(1:n) = sin ( x(1:n)**2 )
  else if ( option == 1 ) then
    f(1:n) = sin ( x(1:n)**2 ) & 
      * exp ( + x(1:n) * x(1:n) )
  else if ( option == 2 ) then
    f(1:n) = sin ( x(1:n)**2 ) &
      * exp ( + 0.5D+00 * x(1:n) * x(1:n) )
  end if

  return
end
subroutine p04_title ( title )

!*****************************************************************************80
!
!! P04_TITLE returns the title for problem 4.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = * ) title

  title = 'sin(x^2)'

  return
end
subroutine p05_exact ( exact )

!*****************************************************************************80
!
!! P05_EXACT returns the exact integral for problem 5.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) EXACT, the value of the integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) exact
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00

  exact = r8_pi / 3.0D+00

  return
end
subroutine p05_fun ( option, n, x, f )

!*****************************************************************************80
!
!! P05_FUN evaluates the integrand for problem 5.
!
!  Discussion:
!
!    The exact value is pi / 3.
!
!    Integral ( -oo < x < +oo ) dx / ( (1+x^2) sqrt(4+3x^2) )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer OPTION:
!    0, integrand is f(x).
!    1, integrand is exp(-x*x) * f(x);
!    2, integrand is exp(-x*x/2) * f(x);
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) X(N), the evaluation points.
!
!    Output, real ( kind = rk ) F(N), the function values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) f(n)
  integer option
  real ( kind = rk ) x(n)

  if ( option == 0 ) then
    f(1:n) = 1.0D+00 / ( ( 1.0D+00 + x(1:n)**2 ) &
      * sqrt ( 4.0D+00 + 3.0D+00 * x(1:n)**2 ) )
  else if ( option == 1 ) then
    f(1:n) = 1.0D+00 / ( ( 1.0D+00 + x(1:n)**2 ) &
      * sqrt ( 4.0D+00 + 3.0D+00 * x(1:n)**2 ) ) &
      * exp ( x(1:n)**2 )
  else if ( option == 2 ) then
    f(1:n) = 1.0D+00 / ( ( 1.0D+00 + x(1:n)**2 ) &
      * sqrt ( 4.0D+00 + 3.0D+00 * x(1:n)**2 ) ) &
      * exp ( 0.5D+00 * x(1:n)**2 )
  end if

  return
end
subroutine p05_title ( title )

!*****************************************************************************80
!
!! P05_TITLE returns the title for problem 5.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = * ) title

  title = '1/( (1+x^2) sqrt(4+3x^2) )'

  return
end
subroutine p06_exact ( exact )

!*****************************************************************************80
!
!! P06_EXACT returns the exact integral for problem 6.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) EXACT, the value of the integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) exact
  integer i4_factorial2
  integer m
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00

  call p06_param ( 'G', 'M', m )

  if ( m <= -1 ) then

    exact = - huge ( exact )

  else if ( mod ( m, 2 ) == 1 ) then

    exact = 0.0D+00

  else

    exact = real ( i4_factorial2 ( m - 1 ), kind = rk ) * sqrt ( r8_pi ) &
      / 2.0D+00**( m / 2 )

  end if

  return
end
subroutine p06_fun ( option, n, x, f )

!*****************************************************************************80
!
!! P06_FUN evaluates the integrand for problem 6.
!
!  Discussion:
!
!    The exact value is (m-1)!! * sqrt ( pi ) / sqrt ( 2^m ).
!
!    Integral ( -oo < x < +oo ) x^m exp (-x*x) dx
!
!    The parameter M is set by calling P06_PARAM.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer OPTION:
!    0, integrand is f(x).
!    1, integrand is exp(-x*x) * f(x);
!    2, integrand is exp(-x*x/2) * f(x);
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) X(N), the evaluation points.
!
!    Output, real ( kind = rk ) F(N), the function values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) f(n)
  integer m
  integer option
  real ( kind = rk ) x(n)

  call p06_param ( 'G', 'M', m )

  if ( option == 0 ) then
    f(1:n) = x(1:n) ** m &
      * exp ( - x(1:n)**2 )
  else if ( option == 1 ) then
    f(1:n) = x(1:n) ** m
  else if ( option == 2 ) then
    f(1:n) = x(1:n) ** m &
      * exp ( - 0.5D+00 * x(1:n)**2 )
  end if

  return
end
subroutine p06_param ( action, name, value )

!*****************************************************************************80
!
!! P06_PARAM gets or sets parameters for problem 6.
!
!  Discussion:
!
!    The parameter is named "M", and it represents the value of the exponent
!    in the integrand function:
!
!    Integral ( -oo < x < +oo ) x^m exp (-x*x) dx
!
!    M must be greater than -1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ACTION, the action.
!    'S' to set the value,
!    'G' to get the value.
!
!    Input, character NAME, the parameter name.
!    'M', the exponent.
!
!    Input/output, integer VALUE, the parameter value.
!    If ACTION = 'S', then VALUE is an input quantity, and M is set to VALUE.
!    If ACTION = 'G', then VALUE is an output quantity, and VALUE is set to M.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character action
  character name
  integer, save :: m = 0
  integer value

  if ( action == 'S' .or. action == 's' ) then

    if ( value <= -1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P06_PARAM - Fatal error!'
      write ( *, '(a)' ) '  Parameter M must be greater than -1.'
      stop
    end if

    m = value

  else if ( action == 'G' .or. action == 'g' ) then

    value = m

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P06_PARAM - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized value of ACTION = "' // action // '".'
    stop

  end if

  return
end
subroutine p06_title ( title )

!*****************************************************************************80
!
!! P06_TITLE returns the title for problem 6.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = * ) title

  title = 'x^m exp(-x*x)'

  return
end
subroutine p07_exact ( exact )

!*****************************************************************************80
!
!! P07_EXACT returns the exact integral for problem 7.
!
!  Discussion:
!
!    The 20 digit value of e^(1/4) was computed by Mathematica.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) EXACT, the value of the integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), parameter :: e_sqrt_sqrt = 1.2840254166877414841D+00
  real ( kind = rk ) exact
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00

  exact = 0.25D+00 * sqrt ( r8_pi ) / e_sqrt_sqrt

  return
end
subroutine p07_fun ( option, n, x, f )

!*****************************************************************************80
!
!! P07_FUN evaluates the integrand for problem 7.
!
!  Discussion:
!
!    The exact value is (1/4) sqrt(pi) / sqrt(sqrt(e)).
!
!    Integral ( -oo < x < +oo ) x^2 cos(x) e^(-x^2) dx
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Prem Kythe, Michael Schaeferkotter,
!    Handbook of Computational Methods for Integration,
!    Chapman and Hall, 2004,
!    ISBN: 1-58488-428-2,
!    LC: QA299.3.K98.
!
!  Parameters:
!
!    Input, integer OPTION:
!    0, integrand is f(x).
!    1, integrand is exp(-x*x) * f(x);
!    2, integrand is exp(-x*x/2) * f(x);
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) X(N), the evaluation points.
!
!    Output, real ( kind = rk ) F(N), the function values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) f(n)
  integer option
  real ( kind = rk ) x(n)

  if ( option == 0 ) then
    f(1:n) = x(1:n)**2 * cos ( x(1:n) ) * exp ( - x(1:n)**2 )
  else if ( option == 1 ) then
    f(1:n) = x(1:n)**2 * cos ( x(1:n) )
  else if ( option == 2 ) then
    f(1:n) = x(1:n)**2 * cos ( x(1:n) ) * exp ( - x(1:n)**2 / 2.0D+00 )
  end if

  return
end
subroutine p07_title ( title )

!*****************************************************************************80
!
!! P07_TITLE returns the title for problem 7.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = * ) title

  title = 'x^2 cos ( x ) exp(-x*x)'

  return
end
subroutine p08_exact ( exact )

!*****************************************************************************80
!
!! P08_EXACT returns the exact integral for problem 8.
!
!  Discussion:
!
!    The 20 digit value of the answer was computed by Mathematica.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 July 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) EXACT, the value of the integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) exact

  exact = 3.0088235661136433510D+00

  return
end
subroutine p08_fun ( option, n, x, f )

!*****************************************************************************80
!
!! P08_FUN evaluates the integrand for problem 8.
!
!  Discussion:
!
!    The exact value is sqrt ( 2 pi ) * HypergeometricU ( -1/2, 0, 1 ).
!
!    Integral ( -oo < x < +oo ) sqrt(1+x*x/2) * exp(-x*x/2) dx
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 July 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Prem Kythe, Michael Schaeferkotter,
!    Handbook of Computational Methods for Integration,
!    Chapman and Hall, 2004,
!    ISBN: 1-58488-428-2,
!    LC: QA299.3.K98.
!
!  Parameters:
!
!    Input, integer OPTION:
!    0, integrand is f(x).
!    1, integrand is exp(-x*x) * f(x);
!    2, integrand is exp(-x*x/2) * f(x);
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) X(N), the evaluation points.
!
!    Output, real ( kind = rk ) F(N), the function values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) f(n)
  integer option
  real ( kind = rk ) x(n)

  if ( option == 0 ) then
    f(1:n) = sqrt ( 1.0D+00 + 0.5D+00 * x(1:n)**2 ) &
      * exp ( - 0.5D+00 * x(1:n)**2 )
  else if ( option == 1 ) then
    f(1:n) = sqrt ( 1.0D+00 + 0.5D+00 * x(1:n)**2 ) &
      * exp ( + 0.5D+00 * x(1:n)**2 )
  else if ( option == 2 ) then
    f(1:n) = sqrt ( 1.0D+00 + 0.5D+00 * x(1:n)**2 )
  end if

  return
end
subroutine p08_title ( title )

!*****************************************************************************80
!
!! P08_TITLE returns the title for problem 8.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 July 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = * ) title

  title = 'sqrt(1+x*x/2) * exp(-x*x/2)'

  return
end
subroutine r8vec_normal_01 ( n, x )

!*****************************************************************************80
!
!! R8VEC_NORMAL_01 returns a unit pseudonormal R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 May 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of values desired.
!
!    Output, real ( kind = rk ) X(N), a sample of the standard normal PDF.
!
!  Local:
!
!    Local, real ( kind = rk ) R(N+1), is used to store some uniform
!    random values.  Its dimension is N+1, but really it is only needed
!    to be the smallest even number greater than or equal to N.
!
!    Local, integer X_LO_INDEX, X_HI_INDEX, records the range
!    of entries of X that we need to compute
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer m
  real ( kind = rk ) r(n+1)
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) x(n)
  integer x_hi_index
  integer x_lo_index
!
!  Record the range of X we need to fill in.
!
  x_lo_index = 1
  x_hi_index = n
!
!  If we need just one new value, do that here to avoid null arrays.
!
  if ( x_hi_index - x_lo_index + 1 == 1 ) then

    call random_number ( harvest = r(1:2) )

    x(x_hi_index) = &
      sqrt ( - 2.0D+00 * log ( r(1) ) ) * cos ( 2.0D+00 * r8_pi * r(2) )
!
!  If we require an even number of values, that's easy.
!
  else if ( mod ( x_hi_index - x_lo_index, 2 ) == 1 ) then

    m = ( x_hi_index - x_lo_index + 1 ) / 2

    call random_number ( harvest = r(1:2*m) )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( - 2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * cos ( 2.0D+00 * r8_pi * r(2:2*m:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( - 2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * sin ( 2.0D+00 * r8_pi * r(2:2*m:2) )
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
      sqrt ( - 2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * cos ( 2.0D+00 * r8_pi * r(2:2*m-2:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( - 2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * sin ( 2.0D+00 * r8_pi * r(2:2*m-2:2) )

    x(n) = sqrt ( - 2.0D+00 * log ( r(2*m-1) ) ) &
      * cos ( 2.0D+00 * r8_pi * r(2*m) )

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

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
