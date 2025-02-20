program main

!*****************************************************************************80
!
!! fsolve_test() tests fsolve().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 November 2023
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'fsolve_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test fsolve(), which solves systems of nonlinear equations.'

  call fsolve_test1 ( )
  call fsolve_test2 ( )
  call fsolve_test3 ( )
  call fsolve_test4 ( )

  call predator_prey_be_test ( )
  call predator_prey_tr_test ( )

  call stiff_bdf2_test ( )
  call stiff_be_test ( )
  call stiff_tr_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'fsolve_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine fsolve_test1 ( )

!*****************************************************************************80
!
!! fsolve_test1() tests fsolve() on a system of 1 equation.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 April 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 1

  external f1
  real ( kind = rk ) fx(n)
  integer info
  real ( kind = rk ) tol
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'fsolve_test1():'
  write ( *, '(a)' ) '  fsolve() solves a nonlinear system of 1 equation.'

  x(1:1) = (/ 0.0D+00 /)

  call f1 ( n, x, fx )      
  call r8vec2_print ( n, x, fx, '  Initial X and F(X)' )

  tol = 0.00001D+00

  call fsolve ( f1, n, x, fx, tol, info )

  write ( *, '(a)' ) ' '
  if ( info == 1 ) then
    write ( *, '(a)' ) '  Satisfactory computation.'
  else
    write ( *, '(a,i6)' ) '  Returned value of INFO = ', info
  end if

  call r8vec2_print ( n, x, fx, '  Final X and F(X)' )

  return
end
subroutine f1 ( n, x, fx )

!*****************************************************************************80
!
!! f1() evaluates a nonlinear system of 1 equation.
!
!  Discussion:
!
!    This is Kepler's equation.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 April 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of variables.
!
!    real ( kind = rk ) X(N), the variable values.
!
!  Output:
!
!    real ( kind = rk ) FX(N), the function values at X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) e
  real ( kind = rk ) fx(n)
  real ( kind = rk ) m
  real ( kind = rk ) x(n)

  e = 0.8D+00
  m = 5.0D+00
  fx(1) = x(1) - m - e * sin ( x(1) )

  return
end
subroutine fsolve_test2 ( )

!*****************************************************************************80
!
!! fsolve_test2() tests fsolve() on a system of 2 equations.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 April 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 2

  external f2
  real ( kind = rk ) fx(n)
  integer info
  real ( kind = rk ) tol
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'fsolve_test2():'
  write ( *, '(a)' ) '  fsolve() solves a nonlinear system of 2 equations.'

  x(1:2) = (/ 3.0D+00, 0.0D+00 /)

  call f2 ( n, x, fx )      
  call r8vec2_print ( n, x, fx, '  Initial X and F(X)' )

  tol = 0.00001D+00

  call fsolve ( f2, n, x, fx, tol, info )

  write ( *, '(a)' ) ' '
  if ( info == 1 ) then
    write ( *, '(a)' ) '  Satisfactory computation.'
  else
    write ( *, '(a,i6)' ) '  Returned value of INFO = ', info
  end if

  call r8vec2_print ( n, x, fx, '  Final X and F(X)' )

  return
end
subroutine f2 ( n, x, fx )

!*****************************************************************************80
!
!! f2() evaluates a nonlinear system of 2 equations.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 August 2016
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of variables.
!
!    real ( kind = rk ) X(N), the variable values.
!
!  Output:
!
!    real ( kind = rk ) FX(N), the function values at X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) fx(n)
  real ( kind = rk ) x(n)

  fx(1) = x(1) * x(1) - 10.0D+00 * x(1) + x(2) * x(2) + 8.0D+00
  fx(2) = x(1) * x(2) * x(2) + x(1) - 10.0D+00 * x(2) + 8.0D+00

  return
end
subroutine fsolve_test3 ( )

!*****************************************************************************80
!
!! fsolve_test3() tests fsolve() on a system of 4 equations.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 April 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 4

  external f3
  real ( kind = rk ) fx(n)
  integer info
  real ( kind = rk ) tol
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'fsolve_test3():'
  write ( *, '(a)' ) '  fsolve() solves a nonlinear system of 4 equations.'

  x(1:n) = 0.0D+00

  call f3 ( n, x, fx )      
  call r8vec2_print ( n, x, fx, '  Initial X and F(X)' )

  tol = 0.00001D+00

  call fsolve ( f3, n, x, fx, tol, info )

  write ( *, '(a)' ) ' '
  if ( info == 1 ) then
    write ( *, '(a)' ) '  Satisfactory computation.'
  else
    write ( *, '(a,i6)' ) '  Returned value of INFO = ', info
  end if

  call r8vec2_print ( n, x, fx, '  Final X and F(X)' )

  return
end
subroutine f3 ( n, x, fx )

!*****************************************************************************80
!
!! f3() evaluates a nonlinear system of 4 equations.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 April 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of variables.
!
!    real ( kind = rk ) X(N), the variable values.
!
!  Output:
!
!    real ( kind = rk ) FX(N), the function values at X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) fx(n)
  integer i
  real ( kind = rk ) x(n)

  do i = 1, n
    fx(i) = ( x(i) - i )**2
  end do

  return
end
subroutine fsolve_test4 ( )

!*****************************************************************************80
!
!! fsolve_test4() tests fsolve() on a system of 8 equations.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 April 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 8

  external f4
  real ( kind = rk ) fx(n)
  integer info
  real ( kind = rk ) tol
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'fsolve_test4():'
  write ( *, '(a)' ) '  fsolve() solves a nonlinear system of 8 equations.'

  x(1:n) = 0.0D+00

  call f4 ( n, x, fx )      
  call r8vec2_print ( n, x, fx, '  Initial X and F(X)' )

  tol = 0.00001D+00

  call fsolve ( f4, n, x, fx, tol, info )

  write ( *, '(a)' ) ' '
  if ( info == 1 ) then
    write ( *, '(a)' ) '  Satisfactory computation.'
  else
    write ( *, '(a,i6)' ) '  Returned value of INFO = ', info
  end if

  call r8vec2_print ( n, x, fx, '  Final X and F(X)' )

  return
end
subroutine f4 ( n, x, fx )

!*****************************************************************************80
!
!! f4() evaluates a nonlinear system of 8 equations.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 April 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of variables.
!
!    real ( kind = rk ) X(N), the variable values.
!
!  Output:
!
!    real ( kind = rk ) FX(N), the function values at X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) fx(n)
  integer i
  real ( kind = rk ) x(n)

  do i = 1, n

    fx(i) = ( 3.0D+00 - 2.0D+00 * x(i) ) * x(i) + 1.0D+00

    if ( 1 < i ) then
      fx(i) = fx(i) - x(i-1)
    end if

    if ( i < n ) then
      fx(i) = fx(i) - 2.0D+00 * x(i+1)
    end if

  end do

  return
end
subroutine predator_prey_be_test ( )

!*****************************************************************************80
!
!! predator_prey_be_test() tests fsolve_be() on the predator prey ODE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 November 2023
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 2

  external predator_prey_dydt
  real ( kind = rk ) fm(m)
  integer info
  real ( kind = rk ) tm
  real ( kind = rk ) to
  real ( kind = rk ) tol
  real ( kind = rk ) ym(m)
  real ( kind = rk ) yo(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'predator_prey_be_test():'

  to = 0.0D+00
  yo = (/ 5000.0D+00, 100.0D+00 /)
  tm = 1.0D+00
  ym = yo
  tol = 1.0D-05

  call backward_euler_residual ( predator_prey_dydt, m, to, yo, tm, ym, fm )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Initial ||ym-yo-(tm-to)*dydt(tm,ym)||:'
  write ( *, '(g14.6)' ) sqrt ( sum ( fm**2 ) )

  call fsolve_be ( predator_prey_dydt, m, to, yo, tm, ym, fm, tol, info )
  if ( info /= 1 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a,i6)' ) '  fsolve_be() returned error flag info = ', info
  end if

  call backward_euler_residual ( predator_prey_dydt, m, to, yo, tm, ym, fm )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Final ||ym-yo-(tm-to)*dydt(tm,ym)||:'
  write ( *, '(g14.6)' ) sqrt ( sum ( fm**2 ) )

  return
end
subroutine predator_prey_tr_test ( )

!*****************************************************************************80
!
!! predator_prey_tr_test() tests fsolve_tr() on the predator prey ODE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 November 2023
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 2

  external predator_prey_dydt
  real ( kind = rk ) ft(m)
  integer info
  real ( kind = rk ) tn
  real ( kind = rk ) to
  real ( kind = rk ) tol
  real ( kind = rk ) yn(m)
  real ( kind = rk ) yo(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'predator_prey_tr_test():'

  to = 0.0D+00
  yo = (/ 5000.0D+00, 100.0D+00 /)
  tn = 1.0D+00
  yn = yo
  tol = 1.0D-05

  call trapezoidal_residual ( predator_prey_dydt, m, to, yo, tn, yn, ft )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Initial residual norm:'
  write ( *, '(g14.6)' ) sqrt ( sum ( ft**2 ) )

  call fsolve_tr ( predator_prey_dydt, m, to, yo, tn, yn, ft, tol, info )
  if ( info /= 1 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a,i6)' ) '  fsolve_tr() returned error flag info = ', info
  end if

  call trapezoidal_residual ( predator_prey_dydt, m, to, yo, tn, yn, ft )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Final residual norm:'
  write ( *, '(g14.6)' ) sqrt ( sum ( ft**2 ) )

  return
end
subroutine predator_prey_dydt ( t, y, dydt )

!*****************************************************************************80
!
!! predator_prey_dydt() evaluates the right hand side of the system.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 February 2020
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    George Lindfield, John Penny,
!    Numerical Methods Using MATLAB,
!    Second Edition,
!    Prentice Hall, 1999,
!    ISBN: 0-13-012641-1,
!    LC: QA297.P45.
!
!  Input:
!
!    real ( kind = rk ) t: the current time.
!
!    real ( kind = rk ) y(2): the current solution variables, rabbits and foxes.
!
!  Output:
!
!    real ( kind = rk ) dydt(2): the right hand side of the 2 ODE's.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) dydt(2)
  real ( kind = rk ) t
  real ( kind = rk ) y(2)

  call r8_fake_use ( t )

  dydt(1) =    2.0 * y(1) - 0.001 * y(1) * y(2)
  dydt(2) = - 10.0 * y(2) + 0.002 * y(1) * y(2)

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
!    real ( kind = rk ) X, the variable to be "used".
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) x

  if ( x /= x ) then
    write ( *, '(a)' ) '  r8_fake_use(): variable is NAN.'
  end if

  return
end
subroutine r8vec2_print ( n, a1, a2, title )

!*****************************************************************************80
!
!! r8vec2_print prints an R8VEC2.
!
!  Discussion:
!
!    An R8VEC2 is a dataset consisting of N pairs of R8's, stored
!    as two separate vectors A1 and A2.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of components of the vector.
!
!    real ( kind = rk ) A1(N), A2(N), the vectors to be printed.
!
!    character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a1(n)
  real ( kind = rk ) a2(n)
  integer i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i4,2x,g14.6,2x,g14.6)' ) i, a1(i), a2(i)
  end do

  return
end
subroutine stiff_be_test ( )

!*****************************************************************************80
!
!! stiff_be_test() tests fsolve_be() on the stiff ODE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 November 2023
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 1

  external stiff_dydt
  real ( kind = rk ) fm(m)
  integer info
  real ( kind = rk ) tm
  real ( kind = rk ) to
  real ( kind = rk ) tol
  real ( kind = rk ) ym(m)
  real ( kind = rk ) yo(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'stiff_be_test():'

  to = 0.0D+00
  yo = (/ 0.0D+00 /)
  tm = 0.25
  ym = yo
  tol = 1.0D-05

  call backward_euler_residual ( stiff_dydt, m, to, yo, tm, ym, fm )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Initial ||ym-yo-(tm-to)*dydt(tm,ym)||:'
  write ( *, '(g14.6)' ) sqrt ( sum ( fm**2 ) )

  call fsolve_be ( stiff_dydt, m, to, yo, tm, ym, fm, tol, info )

  if ( info /= 1 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a,i6)' ) '  fsolve_be() returned error flag info = ', info
  end if

  call backward_euler_residual ( stiff_dydt, m, to, yo, tm, ym, fm )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Final ||ym-yo-(tm-to)*dydt(tm,ym)||:'
  write ( *, '(g14.6)' ) sqrt ( sum ( fm**2 ) )

  return
end
subroutine stiff_tr_test ( )

!*****************************************************************************80
!
!! stiff_tr_test() tests fsolve_tr() on the stiff ODE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 November 2023
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 1

  external stiff_dydt
  real ( kind = rk ) ft(m)
  integer info
  real ( kind = rk ) tn
  real ( kind = rk ) to
  real ( kind = rk ) tol
  real ( kind = rk ) yn(m)
  real ( kind = rk ) yo(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'stiff_tr_test():'

  to = 0.0D+00
  yo = (/ 0.0D+00 /)
  tn = 0.25
  yn = yo
  tol = 1.0D-05

  call trapezoidal_residual ( stiff_dydt, m, to, yo, tn, yn, ft )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Initial residual norm:'
  write ( *, '(g14.6)' ) sqrt ( sum ( ft**2 ) )

  call fsolve_tr ( stiff_dydt, m, to, yo, tn, yn, ft, tol, info )

  if ( info /= 1 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a,i6)' ) '  fsolve_tr() returned error flag info = ', info
  end if

  call trapezoidal_residual ( stiff_dydt, m, to, yo, tn, yn, ft )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Final residual norm:'
  write ( *, '(g14.6)' ) sqrt ( sum ( ft**2 ) )

  return
end
subroutine stiff_dydt ( t, y, dydt )

!*****************************************************************************80
!
!! stiff_dydt() evaluates the right hand side of the stiff ODE.
!
!  Discussion:
!
!    y' = 50 * ( cos(t) - y )
!    y(0) = 0
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 November 2023
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) T, Y(1): the time and solution value.
!
!  Output:
!
!    real ( kind = rk ) DYDT(1): the derivative value.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) dydt(1)
  real ( kind = rk ) t
  real ( kind = rk ) y(1)

  dydt(1) = 50.0D+00 * ( cos ( t ) - y(1) )

  return
end
subroutine stiff_bdf2_test ( )

!*****************************************************************************80
!
!! stiff_bdf2_test() tests fsolve_bdf2() on the stiff ODE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 November 2023
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 1

  external stiff_dydt
  real ( kind = rk ) fm(m)
  integer info
  real ( kind = rk ) t1
  real ( kind = rk ) t2
  real ( kind = rk ) t3
  real ( kind = rk ) tol
  real ( kind = rk ) y1(m)
  real ( kind = rk ) y2(m)
  real ( kind = rk ) y3(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'stiff_bdf2_test():'

  t1 = 0.0D+00
  y1 = (/ 0.0D+00 /)
  t2 = 0.25
  y2(1) = 50.0D+00 * ( sin ( t2 ) + 50.0D+00 * cos ( t2 ) &
        - 50.0D+00 * exp ( - 50.0D+00 * t2 ) ) / 2501.0D+00
  t3 = 0.5
  y3 = y2
  tol = 1.0D-05

  call bdf2_residual ( stiff_dydt, m, t1, y1, t2, y2, t3, y3, fm )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Initial residual:'
  write ( *, '(g14.6)' ) sqrt ( sum ( fm**2 ) )

  call fsolve_bdf2 ( stiff_dydt, m, t1, y1, t2, y2, t3, y3, fm, tol, info )

  if ( info /= 1 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a,i6)' ) '  fsolve_bdf2() returned error flag info = ', info
  end if

  call bdf2_residual ( stiff_dydt, m, t1, y1, t2, y2, t3, y3, fm )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Final residual:'
  write ( *, '(g14.6)' ) sqrt ( sum ( fm**2 ) )

  return
end
subroutine stiff_exact ( n, t, y )

!*****************************************************************************80
!
!! stiff_exact() evaluates the exact solution of the stiff ODE.
!
!  Discussion:
!
!    y' = 50 * ( cos(t) - y )
!    y(0) = 0
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 April 2020
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N: the number of values.
!
!    real ( kind = rk ) T(N): the evaluation times.
!
!  Output:
!
!    real ( kind = rk ) Y(N): the exact solution values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) t(n)
  real ( kind = rk ) y(n)

  y(1:n) = 50.0D+00 * ( sin ( t ) + 50.0D+00 * cos(t) &
    - 50.0D+00 * exp ( - 50.0D+00 * t ) ) / 2501.0D+00

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

