program main

!*****************************************************************************80
!
!! MAIN is the main program for problem 2.
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
  real ( kind = rk ), external :: f2
  real ( kind = rk ), external :: k2
  real ( kind = rk ) u(n)
  character ( len = 80 ) u_file
  real ( kind = rk ) ua
  real ( kind = rk ) ub
  real ( kind = rk ) x(n)
  character ( len = 80 ) x_file

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PROBLEM2:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  A test problem for FD1D_HEAT_STEADY.'
  write ( *, '(a)' ) '  Low K, then high K, then moderate K.'

  a = 0.0D+00
  b = 1.0D+00

  ua = 0.0D+00
  ub = 1.0D+00

  call fd1d_heat_steady ( n, a, b, ua, ub, k2, f2, x, u )

  x_file = 'problem2_nodes.txt'
  call r8mat_write ( x_file, 1, n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X data written to "' // trim ( x_file ) // '".'

  u_file = 'problem2_values.txt'
  call r8mat_write ( u_file, 1, n, u )

  write ( *, '(a)' ) '  U data written to "' // trim ( u_file ) // '".'
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PROBLEM2:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
function k2 ( x )

!*****************************************************************************80
!
!! K2 evaluates the heat transfer coefficient K(X).
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
!    Output, real ( kind = rk ) K2, the value of K(X).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) k2
  real ( kind = rk ) x

  if ( x < 0.5D+00 ) then
    k2 = 0.25D+00
  else if ( x < 0.75D+00 ) then
    k2 = 4.0D+00
  else
    k2 = 1.0D+00
  end if

  return
end
function f2 ( x )

!*****************************************************************************80
!
!! F2 evaluates the right hand side of the steady state heat equation.
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
!    Output, real ( kind = rk ) F2, the value of F(X).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) f2
  real ( kind = rk ) x

  call r8_fake_use ( x )
  f2 = 0.0D+00

  return
end


