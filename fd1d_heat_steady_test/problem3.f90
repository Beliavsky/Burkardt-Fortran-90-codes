program main

!*****************************************************************************80
!
!! MAIN is the main program for problem 3.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 May 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 21

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ), external :: f3
  real ( kind = rk ), external :: k3
  real ( kind = rk ) u(n)
  character ( len = 80 ) u_file
  real ( kind = rk ) ua
  real ( kind = rk ) ub
  real ( kind = rk ) x(n)
  character ( len = 80 ) x_file

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PROBLEM3:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  A test problem for FD1D_HEAT_STEADY.'
  write ( *, '(a)' ) '  Interior source term.'

  a = 0.0D+00
  b = 1.0D+00

  ua = 0.0D+00
  ub = 100.0D+00

  call fd1d_heat_steady ( n, a, b, ua, ub, k3, f3, x, u )

  x_file = 'problem3_nodes.txt'
  call r8mat_write ( x_file, 1, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X data written to "' // trim ( x_file ) // '".'

  u_file = 'problem3_values.txt'
  call r8mat_write ( u_file, 1, n, u )

  write ( *, '(a)' ) '  U data written to "' // trim ( u_file ) // '".'
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PROBLEM3:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
function k3 ( x )

!*****************************************************************************80
!
!! K3 evaluates the heat transfer coefficient K(X).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, the position.
!
!    Output, real ( kind = rk ) K3, the value of K(X).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) k3
  real ( kind = rk ) x

  call r8_fake_use ( x )

  k3 = 1.0D+00

  return
end
function f3 ( x )

!*****************************************************************************80
!
!! F3 evaluates the right hand side of the steady state heat equation.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, the position.
!
!    Output, real ( kind = rk ) F3, the value of F(X).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) f3
  real ( kind = rk ) x

  if ( x < 0.15D+00 ) then
    f3 = 0.0D+00
  else if ( x < 0.45D+00 ) then
    f3 = 200.0D+00
  else
    f3 = 0.0D+00
  end if

  return
end


