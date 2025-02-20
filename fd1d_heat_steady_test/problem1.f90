program main

!*****************************************************************************80
!
!! MAIN is the main program for problem 1.
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

  integer, parameter :: n = 11

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ), external :: f1
  real ( kind = rk ), external :: k1
  real ( kind = rk ) u(n)
  character ( len = 80 ) u_file
  real ( kind = rk ) ua
  real ( kind = rk ) ub
  real ( kind = rk ) x(n)
  character ( len = 80 ) x_file

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PROBLEM1:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  A test problem for FD1D_HEAT_STEADY.'

  a = 0.0D+00
  b = 1.0D+00

  ua = 0.0D+00
  ub = 1.0D+00

  call fd1d_heat_steady ( n, a, b, ua, ub, k1, f1, x, u )

  x_file = 'problem1_nodes.txt'
  call r8mat_write ( x_file, 1, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X data written to "' // trim ( x_file ) // '".'

  u_file = 'problem1_values.txt'
  call r8mat_write ( u_file, 1, n, u )

  write ( *, '(a)' ) '  U data written to "' // trim ( u_file ) // '".'
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PROBLEM1:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
function k1 ( x )

!*****************************************************************************80
!
!! K1 evaluates the heat transfer coefficient K(X).
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
!    Output, real ( kind = rk ) K1, the value of K(X).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) k1
  real ( kind = rk ) x

  call r8_fake_use ( x )

  k1 = 1.0D+00

  return
end
function f1 ( x )

!*****************************************************************************80
!
!! F1 evaluates the right hand side of the steady state heat equation.
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
!    Output, real ( kind = rk ) F1, the value of F(X).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) f1
  real ( kind = rk ) x

  call r8_fake_use ( x )

  f1 = 0.0D+00

  return
end


