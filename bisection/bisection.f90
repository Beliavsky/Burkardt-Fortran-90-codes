subroutine bisection ( a, b, tol, f, it )

!*****************************************************************************80
!
!! bisection() carries out the bisection method to seek a root of F(X) = 0.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 December 2023
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) A, B: the endpoints of a change of sign interval.
!
!    real ( kind = rk ) TOL: the interval size tolerance.  Once |B-A| < TOL,
!    the iteration will stop.
!
!    function f(x): the name of the routine that evaluates the function.
!
!  Output:
!
!    real ( kind = rk ) A, B: the new endpoints that constitute an
!    change of sign interval no larger than TOL.
!
!    integer IT: the number of bisections.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ), external :: f
  real ( kind = rk ) fa
  real ( kind = rk ) fb
  real ( kind = rk ) fc
  integer it
  real ( kind = rk ) tol

  it = 0

  fa = f ( a )
  if ( fa == 0.0 ) then
    b = a
    return
  end if

  fb = f ( b )
  if ( fb == 0.0 ) then
    a = b
    return
  end if

  if ( 0.0 < fa * fb ) then
    write ( *, '(a)' ) 'bisection(): [A,B] not a change of sign interval!'
    stop ( 1 )
  end if

  do while ( tol < abs ( b - a ) )

    c = ( a + b ) / 2.0
    fc = f(c)
    it = it + 1

    if ( 100 < it ) then
      write ( *, '(a)' ) 'bisection(): Too many iterations!'
      stop ( 1 )
    end if

    if ( fc == 0.0 ) then
      a = c
      b = c
      return
    else if ( 0.0 < fc * fa ) then
      a = c
      fa = fc
    else
      b = c
      fb = fc
    end if

  end do

  return
end
