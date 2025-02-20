program main

!*****************************************************************************80
!
!! bisection_test() tests bisection().
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ), external :: fcos
  real ( kind = rk ), external :: fpoly
  real ( kind = rk ), external :: kepler

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'bisection_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test bisection()'

  a = 0.0
  b = 8.0
  call bisection_example ( a, b, fpoly, 'x^2 - 2*x - 15' )

  a = 0.0
  b = 1.0
  call bisection_example ( a, b, fcos, 'cos(x) - x' )

  a = 0.0
  b = 10.0
  call bisection_example ( a, b, kepler, 'Kepler function'  )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'bisection_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  call timestamp ( )

  return
end
subroutine bisection_example ( a, b, f, f_string )

!*****************************************************************************80
!
!! bisection_example() applies bisection() to a particular example.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 December 2023
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ), external :: f
  character * ( * ) f_string
  real ( kind = rk ) fa
  real ( kind = rk ) fb
  real ( kind = rk ) fx
  integer it
  real ( kind = rk ) tol
  real ( kind = rk ) x

  tol = 10.0 * epsilon ( tol ) * ( b - a )
  call bisection ( a, b, tol, f, it )
  x = ( a + b ) / 2.0
  fa = f(a)
  fb = f(b)
  fx = f(x)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Function = "' // f_string // '"'
  write ( *, '(a,g24.16,a,g24.16)' ) '  a = ', a, ', f(a) = ', fa
  write ( *, '(a,g24.16,a,g24.16)' ) '  b = ', b, ', f(b) = ', fb
  write ( *, '(a,g10.4)' ) '  Interval tolerance = ', tol
  write ( *, '(a,i6)' ) '  Number of bisections = ', it
  write ( *, '(a,g24.16,a,g24.16)' ) '  x = ', x, ', f(x) = ', fx

  return
end
function fcos ( x )

!*****************************************************************************80
!
!! fcos() evaluates the function cos(x)-x.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 December 2023
!
!  Input:
!
!    real ( kind = rk ) x, the argument.
!
!  Output:
!
!    real ( kind = rk ) fcos: the function value.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) fcos
  real ( kind = rk ) x

  fcos = cos ( x ) - x

  return
end
function fpoly ( x )

!*****************************************************************************80
!
!! fpoly() evaluates a polynomial function.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 December 2023
!
!  Input:
!
!    real ( kind = rk ) x, the argument.
!
!  Output:
!
!    real ( kind = rk ) fpoly: the function value.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) fpoly
  real ( kind = rk ) x

  fpoly = x * x - 2.0 * x - 15.0

  return
end
function kepler ( x )

!*****************************************************************************80
!
!! kepler() evaluates a version of Kepler's equation.
!
!  Discussion:
!
!    Kepler's equation relates the mean anomaly M, the eccentric anomaly E,
!    andthe eccentricity e of a planetary orbit.
!
!    Typically, e is a fixed feature of the orbit, the value of M is determined
!    by observation, and the value of E is desired.
!
!    Kepler's equation states that:
!      M = E - e sin(E)
!
!    Suppose we have an orbit with e = 2, and we have observed M = 5.  What is
!    the value of E?  The equation becomes:
!      5 = E - 2 sin ( E ).
!
!    To solve for E, we need to rewrite this as a function:
!      F(E) = 5 - E + 2 sin ( E )
!    and then use a nonlinear equation solver to solve for the value of E
!    such that F(E)=0.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 December 2023
!
!  Input:
!
!    real ( kind = rk ) x, the current estimate for the value of E.
!
!  Output:
!
!    real ( kind = rk ) kepler: the Kepler equation residual F(E).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) kepler
  real ( kind = rk ) x

  kepler = 5.0 - x + 2.0 * sin ( x )

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
!    15 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

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

