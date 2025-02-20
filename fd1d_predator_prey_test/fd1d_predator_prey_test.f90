subroutine u_init ( n, x, u0 )

!*****************************************************************************80
!
!! U_INIT supplies the initial value of U at each node.
!
!  Modified:
!
!    02 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of nodes.
!
!    Input, real ( kind = rk ) X(N), the X coordinates of the nodes.
!
!    Output, real ( kind = rk ) U0(N), the initial value of U at the nodes.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) u0(n)
  real ( kind = rk ) x(n)

  u0(1:n) = exp ( - ( x(1:n) - real ( 100.0, kind = rk ) )**2 ) &
    / real ( 5.0, kind = rk )

  return
end
subroutine v_init ( n, x, v0 )

!*****************************************************************************80
!
!! V_INIT supplies the initial value of V at each node.
!
!  Modified:
!
!    02 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of nodes.
!
!    Input, real ( kind = rk ) X(N), the X coordinates of the nodes.
!
!    Output, real ( kind = rk ) U0(N), the initial value of U at the nodes.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) v0(n)
  real ( kind = rk ) x(n)

  call r8_fake_use ( x(1) )

  v0(1:n) = real ( 2.0, kind = rk ) / real ( 5.0, kind = rk )

  return
end
subroutine r8_fake_use ( x )

!*****************************************************************************80
!
!! r8_fake_use pretends to use a variable.
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
    write ( *, '(a)' ) '  r8_fake_use: variable is NAN.'
  end if

  return
end

