subroutine laguerre ( x0, degree, abserr, kmax, f, x, ierror, k )

!*****************************************************************************80
!
!! laguerre() implements the Laguerre rootfinding method for polynomials.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 March 2024
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Eldon Hansen, Merrell Patrick,
!    A Family of Root Finding Methods,
!    Numerische Mathematik,
!    Volume 27, 1977, pages 257 - 269.
!
!  Input:
!
!    real ( kind = rk ) X0.
!    An initial estimate for the root of the equation.
!
!    integer degree: the polynomial degree of the function, at least 2.
!
!    real ( kind = rk ) ABSERR: an error tolerance.
!
!    integer KMAX: the maximum number of iterations allowed.
!
!    real ( kind = rk ) external F: the name of the routine that
!    evaluates the function or its derivatives, of the form
!      function f ( x, ider )
!
!  Output:
!
!    real ( kind = rk ) X: the estimated solution, if IERROR=0.
!
!    integer IERROR: error indicator.
!    0, no error occurred.
!    nonzero, an error occurred, and the iteration was halted.
!
!    integer K: the number of steps taken.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) abserr
  real ( kind = rk ) beta
  real ( kind = rk ) bot
  real ( kind = rk ) d2fx
  integer degree
  real ( kind = rk ) dfx
  real ( kind = rk ) dx
  real ( kind = rk ), external :: f
  real ( kind = rk ) fx
  integer ierror
  integer k
  integer kmax
  real ( kind = rk ) x
  real ( kind = rk ) x0
  real ( kind = rk ) z
!
!  Initialization.
!
  ierror = 0
  x = x0

  beta = 1.0D+00 / real ( degree - 1, kind = rk )

  k = 0
  fx = f ( x, 0 )
  dfx = f ( x, 1 )
  d2fx = f ( x, 2 )
!
!  Iteration loop:
!
  do
!
!  If the error tolerance is satisfied, then exit.
!
    if ( abs ( fx ) <= abserr ) then
      exit
    end if

    k = k + 1

    if ( kmax < k ) then
      ierror = 2
      return
    end if

    z = dfx**2 - ( beta + 1.0D+00 ) * fx * d2fx
    z = max ( z, 0.0D+00 )

    bot = beta * dfx + sqrt ( z )

    if ( bot == 0.0D+00 ) then
      ierror = 3
      return
    end if
!
!  Set the increment.
!
    dx = - ( beta + 1.0D+00 ) * fx / bot
!
!  Update the iterate and function values.
!
    x = x + dx
    fx = f ( x, 0 )
    dfx = f ( x, 1 )
    d2fx = f ( x, 2 )

  end do

  return
end
