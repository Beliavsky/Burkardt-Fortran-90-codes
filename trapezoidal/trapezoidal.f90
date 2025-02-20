subroutine trapezoidal ( dydt, tspan, y0, n, m, t, y )

!*****************************************************************************80
!
!! trapezoidal() uses the trapezoidal method + fsolve_tr() to solve an ODE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 November 2023
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    external dydt: a subroutine that evaluates the right
!    hand side of the ODE, of the form
!    subroutine dydt ( t, y, dy )
!
!    real ( kind = rk ) tspan(2): the initial and final times.
!
!    real ( kind = rk ) y0(m): the initial condition.
!
!    integer n: the number of steps to take.
!
!    integer m: the number of variables.
!
!  Output:
!
!    real ( kind = rk ) t(n+1), y(n+1,m): the times and solution values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) dt
  external dydt
  real ( kind = rk ) fn(m)
  integer i
  integer info
  real ( kind = rk ) t(n+1)
  real ( kind = rk ) tn
  real ( kind = rk ) to
  real ( kind = rk ) tol
  real ( kind = rk ) tspan(2)
  real ( kind = rk ) y(n+1,m)
  real ( kind = rk ) y0(m)
  real ( kind = rk ) yn(m)
  real ( kind = rk ) yo(m)

  dt = ( tspan(2) - tspan(1) ) / n

  tol = 1.0D-05

  do i = 0, n

    if ( i == 0 ) then

      t(i+1) = tspan(1)
      y(i+1,1:m) = y0(1:m)

    else

      to = t(i)
      yo = y(i,1:m)

      tn = t(i) + dt 
      yn(1:m) = y(i,1:m)
!
!  Call fsolve_tr() to compute yn.
!
      call fsolve_tr ( dydt, m, to, yo, tn, yn, fn, tol, info )

      if ( info /= 1 ) then
        write ( *, '(a)' ) ''
        write ( *, '(a)' ) 'trapezoidal(): Fatal error!'
        write ( *, '(a,i10)' ) '  info = ', info
        stop 1
      end if

      t(i+1) = tn
      y(i+1,1:m) = yn(1:m)

    end if

  end do

  return
end

