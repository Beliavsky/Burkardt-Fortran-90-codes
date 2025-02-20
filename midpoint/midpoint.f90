subroutine midpoint ( dydt, tspan, y0, n, m, t, y )

!*****************************************************************************80
!
!! midpoint() uses the midpoint method + fsolve_be() to solve an ODE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 November 2023
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
!    real ( kind = rk ) tspan(2): contains the initial and final times.
!
!    real ( kind = rk ) y0(m): a column vector containing the initial condition.
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
  real ( kind = rk ) fm(m)
  integer i
  integer info
  real ( kind = rk ) t(n+1)
  real ( kind = rk ) theta
  real ( kind = rk ) tm
  real ( kind = rk ) to
  real ( kind = rk ) tol
  real ( kind = rk ) tspan(2)
  real ( kind = rk ) y(n+1,m)
  real ( kind = rk ) y0(m)
  real ( kind = rk ) ym(m)
  real ( kind = rk ) yo(m)

  dt = ( tspan(2) - tspan(1) ) / n

  theta = 0.5D+00
  tol = 1.0D-05

  do i = 0, n

    if ( i == 0 ) then

      t(i+1) = tspan(1)
      y(i+1,1:m) = y0(1:m)

    else

      to = t(i)
      yo = y(i,1:m)

      tm = t(i) + theta * dt 
      ym(1:m) = y(i,1:m)
!
!  Call fsolve_be() to compute ym.
!
      call fsolve_be ( dydt, m, to, yo, tm, ym, fm, tol, info )

      if ( info /= 1 ) then
        write ( *, '(a)' ) ''
        write ( *, '(a)' ) 'midpoint(): Fatal error!'
        write ( *, '(a,i10)' ) '  info = ', info
        stop 1
      end if

      t(i+1) = t(i) + dt
      y(i+1,1:m) = (           1.0D+00 / theta ) * ym(1:m) &
                 + ( 1.0D+00 - 1.0D+00 / theta ) * y(i,1:m)

    end if

  end do

  return
end

