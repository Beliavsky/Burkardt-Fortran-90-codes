subroutine bdf2 ( dydt, tspan, y0, n, m, t, y )

!*****************************************************************************80
!
!! bdf2() uses the BDF2 method + fsolve_bdf2() to solve an ODE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 November 2023
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

  real ( kind = rk ) dely(m)
  real ( kind = rk ) dt
  external dydt
  real ( kind = rk ) fn(m)
  integer i
  integer info
  real ( kind = rk ) t(n+1)
  real ( kind = rk ) t1
  real ( kind = rk ) t2
  real ( kind = rk ) t3
  real ( kind = rk ) tol
  real ( kind = rk ) tspan(2)
  real ( kind = rk ) y(n+1,m)
  real ( kind = rk ) y0(m)
  real ( kind = rk ) y1(m)
  real ( kind = rk ) y2(m)
  real ( kind = rk ) y3(m)

  dt = ( tspan(2) - tspan(1) ) / n

  tol = 1.0D-05

  do i = 0, n

    if ( i == 0 ) then

      t(i+1) = tspan(1)
      y(i+1,1:m) = y0(1:m)

    else if ( i == 1 ) then

      t1 = t(i)
      y1 = y(i,1:m)

      t2 = t1 + dt
      call dydt ( t1, y1, dely )
      y2(1:m) = y1(1:m) + dt * dely(1:m)
      call fsolve_be ( dydt, m, t1, y1, t2, y2, fn, tol, info )

      t(i+1) = t2
      y(i+1,1:m) = y2(1:m)

    else

      t1 = t(i-1)
      y1(1:m) = y(i-1,1:m)
      t2 = t1 + dt
      y2(1:m) = y(i,1:m);
      t3 = t2 + dt 
      call dydt ( t2, y2, dely )
      y3(1:m) = y(i,1:m) + dt * dely(1:m)

      call fsolve_bdf2 ( dydt, m, t1, y1, t2, y2, t3, y3, fn, tol, info )

      if ( info /= 1 ) then
        write ( *, '(a)' ) ''
        write ( *, '(a)' ) 'bdf2(): Fatal error!'
        write ( *, '(a,i10)' ) '  info = ', info
        stop 1
      end if

      t(i+1) = t3
      y(i+1,1:m) = y3(1:m)

    end if

  end do

  return
end

