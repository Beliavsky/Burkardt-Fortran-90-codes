program main

!*****************************************************************************80
!
!! toms178_test() tests toms178().
!
!  Modified:
!
!    12 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'toms178_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test TOMS178().'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS178_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests HOOKE with the Rosenbrock function.
!
!  Modified:
!
!    12 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nvars = 2

  real ( kind = rk ) endpt(nvars)
  real ( kind = rk ) eps
  integer hooke
  integer i
  integer it
  integer itermax
  real ( kind = rk ) rho
  real ( kind = rk ), external :: rosenbrock
  real ( kind = rk ) startpt(nvars)
  real ( kind = rk ) value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'toms178_test01'
  write ( *, '(a)' ) '  hooke seeks a minimizer of F(X).'
  write ( *, '(a)' ) '  Here we use the Rosenbrock function.'
!
!  Starting guess for Rosenbrock.
!
  startpt(1) = -1.2D+00
  startpt(2) = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Initial estimate X = '
  write ( *, '(a)' ) ' '
  do i = 1, nvars
    write ( *, '(2x,i8,2x,g14.6)' ) i, startpt(i)
  end do

  value = rosenbrock ( startpt, nvars )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X) = ', value
!
!  Call HOOKE.
!
  itermax = 5000
  rho = 0.5D+00
  eps = 1.0D-06

  it = hooke ( nvars, startpt, endpt, rho, eps, itermax, rosenbrock )
!
!  Results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of iterations taken = ', it
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X* = '
  write ( *, '(a)' ) ' '
  do i = 1, nvars
    write ( *, '(2x,i8,2x,g14.6)' ) i, endpt(i)
  end do

  value = rosenbrock ( endpt, nvars )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X*) = ', value

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests HOOKE with the WOODS function.
!
!  Discussion:
!
!    The Hooke and Jeeves algorithm works well when RHO = 0.5, but
!    does poorly when RHO = 0.6, and better when RHO = 0.8
!
!  Modified:
!
!    12 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nvars = 4

  real ( kind = rk ) endpt(nvars)
  real ( kind = rk ) eps
  integer hooke
  integer i
  integer it
  integer itermax
  real ( kind = rk ) rho
  real ( kind = rk ) startpt(nvars)
  real ( kind = rk ) value
  real ( kind = rk ), external :: woods

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  HOOKE seeks a minimizer of F(X).'
  write ( *, '(a)' ) '  Here we use the Woods function.'
!
!  Starting guess.
!
  startpt(1) = -3.0D+00
  startpt(2) = -1.0D+00
  startpt(3) = -3.0D+00
  startpt(4) = -1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Initial estimate X = '
  write ( *, '(a)' ) ' '
  do i = 1, nvars
    write ( *, '(2x,i8,2x,g14.6)' ) i, startpt(i)
  end do

  value = woods ( startpt, nvars )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X) = ', value
!
!  Call HOOKE.
!
  itermax = 5000
  rho = 0.5D+00
  eps = 1.0D-06

  it = hooke ( nvars, startpt, endpt, rho, eps, itermax, woods )
!
!  Results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of iterations taken = ', it
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X* = '
  write ( *, '(a)' ) ' '
  do i = 1, nvars
    write ( *, '(2x,i8,2x,g14.6)' ) i, endpt(i)
  end do

  value = woods ( endpt, nvars )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X*) = ', value

  return
end
function rosenbrock ( x, n )

!*****************************************************************************80
!
!! ROSENBROCK evaluates the Rosenbrock function.
!
!  Discussion:
!
!    The Hooke and Jeeves algorithm works reasonably well on
!    Rosenbrock's test function, depending on the value of RHO chosen.
!
!  Modified:
!
!    12 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) X(N), the argument of the function.
!
!    Input, integer N, the spatial dimension.
!
!    Output, real ( kind = rk ) ROSENBROCK, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) rosenbrock
  real ( kind = rk ) x(n)

  rosenbrock = 100.0 * ( x(2) - x(1) * x(1) )**2 &
             +         ( 1.0D+00 - x(1) )**2

  return
end
function woods ( x, n )

!*****************************************************************************80
!
!! WOODS evaluates the Woods function.
!
!  Modified:
!
!    12 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) X(N), the argument of the function.
!
!    Input, integer N, the spatial dimension.
!
!    Output, real ( kind = rk ) WOODS, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) s1
  real ( kind = rk ) s2
  real ( kind = rk ) s3
  real ( kind = rk ) t1
  real ( kind = rk ) t2
  real ( kind = rk ) t3
  real ( kind = rk ) t4
  real ( kind = rk ) t5
  real ( kind = rk ) woods
  real ( kind = rk ) x(n)

  s1 = x(2) - x(1) * x(1)
  s2 = 1.0D+00 - x(1)
  s3 = x(2) - 1.0D+00
  t1 = x(4) - x(3) * x(3)
  t2 = 1.0D+00 - x(3)
  t3 = x(4) - 1.0D+00
  t4 = s3 + t3
  t5 = s3 - t3

  woods = 100.0D+00 * s1**2 &
        +             s2**2 &
        +  90.0D+00 * t1**2 &
        +             t2**2 &
        +  10.0D+00 * t4**2 &
        +   0.1D+00 * t5**2

  return
end

