function polynomial_root_bound ( n, c )

!*****************************************************************************80
!
!! polynomial_root_bound() bounds the roots of a polynomial.
!
!  Discussion:
!
!    The method is due to Cauchy, and the bound is called the Cauchy bound.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 December 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    John D Cook,
!    Bounding complex roots by a positive root,
!    https://www.johndcook.com/blog/2023/12/02/cauchy-radius/
!    02 December 2023.
!
!  Input:
!
!    integer n: the degree of the polynomial.
!
!    complex c(0:n): the coefficients of the polynomial.
!    c(n) multiplies x^n and c(0) is the constant term (multiplying x^0).
!
!  Output:
!
!    real b: a bound on the norm of all roots of the polynomial.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  complex ( kind = ck ) c(0:n)
  real ( kind = rk ) fx
  integer it
  integer job
  real ( kind = rk ) polynomial_root_bound
  real ( kind = rk ) polyval
  real ( kind = rk ), allocatable :: q(:)
  real ( kind = rk ) qneg
  real ( kind = rk ) qpos
  real ( kind = rk ) x
  real ( kind = rk ) xneg
  real ( kind = rk ) xpos
!
!  If the constant term is zero, the problem reduces to finding roots of p(z)/z.
!  Shift coefficients, repeatedly if necessary.
!
  do while ( .true. )
!
!  Polynomial could be identically zero.
!
    if ( n + 1 <= 0 ) then
      polynomial_root_bound = 0.0
      write ( *, '(a)' ) 'polynomial_root_bound(): p(z) = 0.'
      return
    end if

    if ( c(0) /= 0.0 ) then
      exit
    end if
!
!  Drop the first coefficient, which is zero and reduce the order.
!
    n = n - 1
    c(0:n) = c(1:n+1)

  end do
!
!  Define associated real polynomial q(x):
!    q(x) = |cn| x^(n) - |cn-1| x^(n-1) - |cn-2| x^(n-2) - ... - |c0|
!
  allocate ( q(0:n) )

  q(0:n) = 0.0

  q(n) = abs ( c(n) )
  q(0:n-1) = - abs ( c(0:n-1) )
!
!  Determine change of sign interval for Q(x).
!
  xneg = 0.0
  qneg = q(0)

  xpos = 1.0
  qpos = polyval ( n, q, xpos )
  do while ( qpos <= 0.0 )
    xpos = xpos * 2.0
    qpos = polyval ( n, q, xpos )
  end do
!
!  Use bisection to find X such that Q(X)=0.
!
  job = 0
  it = 0
  do while ( .true. )
    call bisection_rc ( xneg, xpos, x, fx, job )
    fx = polyval ( n, q, x )
    it = it + 1
    if ( 100 < it ) then
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'polynomial_root_bound(): Warning!'
      write ( *, '(a)' ) '  Too many bisection steps!'
      exit
    end if
    if ( abs ( fx ) < 1.0E-10 ) then
      write ( *, '(a,i3,a)' ) 'Convergence after ', it, ' iterations'
      exit
    end if
  end do

  polynomial_root_bound = 0.5 * ( xneg + xpos )

  deallocate ( q )

  return
end
subroutine bisection_rc ( a, b, x, fx, job )

!*****************************************************************************80
!
!! bisection_rc() seeks a zero of f(x) in a change of sign interval.
!
!  Discussion:
!
!    The bisection method is used.
!
!    This routine uses reverse communication, so that the function is always
!    evaluated in the calling program.
!
!    On the first call, the user sets JOB = 0, and the values of A and B.
!    Thereafter, the user checks the returned value of JOB and follows 
!    directions.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) A, B: the change of sign interval.
!
!    real ( kind = rk ) FX:
!    on first call, with JOB = 0, a value of FX is not needed.
!    Thereafter, FX should be set to f(x), where x is output of previous call.
!
!    integer JOB: set to 0 on first call, to request initialization.
!
!  Output:
!
!    real ( kind = rk ) A, B: the updated change of sign interval.
!
!    real ( kind = rk ) X: the next point at which a function value is needed.
!
!    integer JOB: reset to 1, indicating that initialization was completed.
!  
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ), save :: fa
  real ( kind = rk ), save :: fb
  real ( kind = rk ) fx
  integer job
  integer, save :: state
  real ( kind = rk ) x
!
!  Initialize.  
!    Accept value of a.
!    Request value of f(a).
!
  if ( job == 0 ) then

    fa = 0.0D+00
    fb = 0.0D+00
    state = 1
    x = a
    job = 1
!
!  Accept value of f(a).
!  Request value of f(b).
!
  else if ( state == 1 ) then

    fa = fx
    x = b
    state = 2
!
!  Accept value of f(b).
!  Request value of f at midpoint.
!
  else if ( state == 2 ) then

    if ( fx == 0.0 ) then
      a = b
      state = 3
      return
    end if

    fb = fx

    if ( 0.0 < fa * fb ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'bisection_rc(): Fatal error!'
      write ( *, '(a)' ) '  f(a) and f(b) have the same sign.'
      stop ( 1 )
    end if

    x = ( a + b ) / 2.0D+00
    state = 3
!
!  Use sign of function value to bisect current interval.
!  Determine new midpoint and request function value.
!
  else

    if ( fx == 0.0 ) then
      a = x
      b = x
      return
    end if

    if ( 0.0 < fx * fa ) then
      a = x
      fa = fx
    else
      b = x
      fb = fx
    end if
    x = ( a + b ) / 2.0D+00
    state = 3

  end if

  return
end
function polyval ( d, c, x )

!*****************************************************************************80
!
!! polyval() evaluates a polynomial using a naive method.
!
!  Discussion:
!
!    The polynomial
!
!      p(x) = c0 + c1 * x + c2 * x^2 + ... + cd * x^d
!
!    is to be evaluated at X.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 August 2017
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer d: the degree.
!
!    real ( kind = rk ) c(0:d): the polynomial coefficients.
!    C(I) is the coefficient of X^I.
!
!    real ( kind = rk ) X: the evaluation point.
!
!  Output:
!
!    real ( kind = rk ) polyval: the polynomial value.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer d

  real ( kind = rk ) c(0:d)
  integer i
  real ( kind = rk ) polyval
  real ( kind = rk ) value
  real ( kind = rk ) x
  real ( kind = rk ) xi

  value = c(0)
  xi = 1.0D+00
  do i = 1, d
    xi = xi * x
    value = value + c(i) * xi
  end do

  polyval = value

  return
end
