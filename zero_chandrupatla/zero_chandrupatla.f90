subroutine zero_chandrupatla ( f, x1, x2, xm, fm, calls )

!*****************************************************************************80
!
!! zero_chandrupatla() seeks a zero of a function using Chandrupatla's algorithm.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 March 2024
!
!  Author:
!
!    Original QBASIC version by Tirupathi Chandrupatla.
!    This version by John Burkardt.
!
!  Reference:
!
!    Tirupathi Chandrupatla,
!    A new hybrid quadratic/bisection algorithm for finding the zero of a
!    nonlinear function without using derivatives,
!    Advances in Engineering Software,
!    Volume 28, Number 3, pages 145-149, 1997.
!
!  Input:
!
!    function f ( x ): the name of the user-supplied function.
!
!    real ( kind = rk8 ) a, b: the endpoints of the change of sign interval.
!
!  Output:
!
!    real ( kind = rk8 ) z, fz: the estimated root and its function value.
!
!    integer calls: the number of function calls.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) al
  real ( kind = rk8 ) b
  real ( kind = rk8 ) c
  integer calls
  real ( kind = rk8 ) d
  real ( kind = rk8 ) delta
  real ( kind = rk8 ) eps
  real ( kind = rk8 ), external :: f 
  real ( kind = rk8 ) f0
  real ( kind = rk8 ) f1
  real ( kind = rk8 ) f2
  real ( kind = rk8 ) f3
  real ( kind = rk8 ) fh
  real ( kind = rk8 ) fl
  real ( kind = rk8 ) fm
  real ( kind = rk8 ) ph
  real ( kind = rk8 ) t
  real ( kind = rk8 ) tl
  real ( kind = rk8 ) tol
  real ( kind = rk8 ) x0
  real ( kind = rk8 ) x1
  real ( kind = rk8 ) x2
  real ( kind = rk8 ) x3
  real ( kind = rk8 ) xi
  real ( kind = rk8 ) xm

  eps = 1.0D-10
  delta = 0.00001

  f1 = f ( x1 )
  f2 = f ( x2 )
  calls = 2

  t = 0.5

  do while ( .true. )

    x0 = x1 + t * ( x2 - x1 )
    f0 = f ( x0 )
    calls = calls + 1
!
!  Arrange 2-1-3: 2-1 Interval, 1 Middle, 3 Discarded point.
!
    if ( ( 0.0 < f0 ) .eqv. ( 0.0 < f1 ) ) then
      x3 = x1
      f3 = f1
      x1 = x0
      f1 = f0
    else
      x3 = x2
      f3 = f2
      x2 = x1
      f2 = f1
      x1 = x0
      f1 = f0
    end if
!
!  Identify the one that approximates zero.
!
    if ( abs ( f2 ) < abs ( f1 ) ) then
      xm = x2
      fm = f2
    else
      xm = x1
      fm = f1
    end if

    tol = 2.0 * eps * abs ( xm ) + 0.5 * delta
    tl = tol / abs ( x2 - x1 )

    if ( 0.5 < tl .or. fm == 0.0 ) then
      exit
    end if
!
!  If inverse quadratic interpolation holds, use it.
!
    xi = ( x1 - x2 ) / ( x3 - x2 )
    ph = ( f1 - f2 ) / ( f3 - f2 )
    fl = 1.0 - sqrt ( 1.0 - xi )
    fh = sqrt ( xi )

    if ( fl < ph .and. ph < fh ) then
      al = ( x3 - x1 ) / ( x2 - x1 )
      a = f1 / ( f2 - f1 )
      b = f3 / ( f2 - f3 )
      c = f1 / ( f3 - f1 )
      d = f2 / ( f3 - f2 )
      t = a * b + c * d * al
    else
      t = 0.5
    end if
!
!  Adjust T away from the interval boundary.
!
    t = max ( t, tl )
    t = min ( t, 1.0 - tl )

  end do

  return
end

