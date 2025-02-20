program main

!*****************************************************************************80
!
!! asa047_test() tests asa047().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'asa047_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test asa047().'

  call rosenbrock_test ( )
  call powell_test ( )
  call helical_test ( )
  call quartic_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'asa047_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine rosenbrock_test ( )

!*****************************************************************************80
!
!! rosenbrock_test() tests nelmin() on rosenbrock().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 2

  integer i
  integer icount
  integer ifault
  integer kcount
  integer konvge
  integer numres
  real ( kind = rk ) reqmin
  real ( kind = rk ), external :: rosenbrock
  real ( kind = rk ) start(n)
  real ( kind = rk ) step(n)
  real ( kind = rk ) xmin(n)
  real ( kind = rk ) ynewlo

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'rosenbrock_test():'
  write ( *, '(a)' ) '  Apply nelmin() to the ROSENBROCK function.'

  start(1:n) = (/ -1.2D+00, 1.0D+00 /)

  reqmin = 1.0D-08

  step(1:n) = (/ 1.0D+00, 1.0D+00 /)

  konvge = 10
  kcount = 500

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Starting point X:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,g14.6)' ) start(i)
  end do

  ynewlo = rosenbrock ( start )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X) = ', ynewlo

  call nelmin ( rosenbrock, n, start, xmin, ynewlo, reqmin, step, &
    konvge, kcount, icount, numres, ifault )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return code IFAULT = ', ifault
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Estimate of minimizing value X*:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,g14.6)' ) xmin(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X*) = ', ynewlo

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of iterations = ', icount
  write ( *, '(a,i8)' ) '  Number of restarts =   ', numres

  return
end
function rosenbrock ( x )

!*****************************************************************************80
!
!! rosenbrock() evaluates the Rosenbrock parabolic value function.
!
!  Discussion:
!
!    Thanks to Vivek Rao for pointing out a discrepancy in the local
!    dimensioning of X.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 August 2021
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    R ONeill,
!    Algorithm AS 47:
!    Function Minimization Using a Simplex Procedure,
!    Applied Statistics,
!    Volume 20, Number 3, 1971, pages 338-345.
!
!  Input:
!
!    real ( kind = rk ) X(2), the argument.
!
!  Output:
!
!    real ( kind = rk ) ROSENBROCK, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) fx
  real ( kind = rk ) fx1
  real ( kind = rk ) fx2
  real ( kind = rk ) rosenbrock
  real ( kind = rk ) x(2)

  fx1 = x(2) - x(1) * x(1)
  fx2 = 1.0D+00 - x(1)

  fx = 100.0D+00 * fx1 * fx1 &
     +             fx2 * fx2

  rosenbrock = fx

  return
end
subroutine powell_test ( )

!*****************************************************************************80
!
!! powell_test() tests nelmin() on powell().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 4

  integer i
  integer icount
  integer ifault
  integer kcount
  integer konvge
  integer numres
  real ( kind = rk ), external :: powell
  real ( kind = rk ) reqmin
  real ( kind = rk ) start(n)
  real ( kind = rk ) step(n)
  real ( kind = rk ) xmin(n)
  real ( kind = rk ) ynewlo

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'powell_test():'
  write ( *, '(a)' ) '  Apply nelmin() to POWELL quartic function.'

  start(1:n) = (/  3.0D+00, - 1.0D+00,   0.0D+00,   1.0D+00 /)

  reqmin = 1.0D-08

  step(1:n) = (/ 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00 /)

  konvge = 10
  kcount = 500

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Starting point X:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,g14.6)' ) start(i)
  end do

  ynewlo = powell ( start )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X) = ', ynewlo

  call nelmin ( powell, n, start, xmin, ynewlo, reqmin, step, &
    konvge, kcount, icount, numres, ifault )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return code IFAULT = ', ifault
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Estimate of minimizing value X*:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,g14.6)' ) xmin(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X*) = ', ynewlo

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of iterations = ', icount
  write ( *, '(a,i8)' ) '  Number of restarts =   ', numres

  return
end
function powell ( x )

!*****************************************************************************80
!
!! powell() evaluates the Powell quartic function.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 August 2021
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    R ONeill,
!    Algorithm AS 47:
!    Function Minimization Using a Simplex Procedure,
!    Applied Statistics,
!    Volume 20, Number 3, 1971, pages 338-345.
!
!  Input:
!
!    real ( kind = rk ) X(4), the argument.
!
!  Output:
!
!    real ( kind = rk ) POWELL, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) fx
  real ( kind = rk ) fx1
  real ( kind = rk ) fx2
  real ( kind = rk ) fx3
  real ( kind = rk ) fx4
  real ( kind = rk ) powell
  real ( kind = rk ) x(4)

  fx1 = x(1) + 10.0D+00 * x(2)
  fx2 = x(3) - x(4)
  fx3 = x(2) - 2.0D+00 * x(3)
  fx4 = x(1) - x(4)

  fx =            fx1 * fx1 &
     +  5.0D+00 * fx2 * fx2 &
     +            fx3 * fx3 * fx3 * fx3 &
     + 10.0D+00 * fx4 * fx4 * fx4 * fx4

  powell = fx

  return
end
subroutine helical_test ( )

!*****************************************************************************80
!
!! helical_test() tests nelmin() on helical().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 3

  real ( kind = rk ), external :: helical
  integer i
  integer icount
  integer ifault
  integer kcount
  integer konvge
  integer numres
  real ( kind = rk ) reqmin
  real ( kind = rk ) start(n)
  real ( kind = rk ) step(n)
  real ( kind = rk ) xmin(n)
  real ( kind = rk ) ynewlo

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'helical_test():'
  write ( *, '(a)' ) '  Test nelmin() on the HELICAL function.'

  start(1:n) = (/ - 1.0D+00,   0.0D+00,   0.0D+00 /)

  reqmin = 1.0D-08

  step(1:n) = (/ 1.0D+00, 1.0D+00, 1.0D+00 /)

  konvge = 10
  kcount = 500

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Starting point X:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,g14.6)' ) start(i)
  end do

  ynewlo = helical ( start )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X) = ', ynewlo

  call nelmin ( helical, n, start, xmin, ynewlo, reqmin, step, &
    konvge, kcount, icount, numres, ifault )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return code IFAULT = ', ifault
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Estimate of minimizing value X*:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,g14.6)' ) xmin(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X*) = ', ynewlo

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of iterations = ', icount
  write ( *, '(a,i8)' ) '  Number of restarts =   ', numres

  return
end
function helical ( x )

!*****************************************************************************80
!
!! helical() evaluates the Fletcher-Powell helical valley function.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 August 2021
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    R ONeill,
!    Algorithm AS 47:
!    Function Minimization Using a Simplex Procedure,
!    Applied Statistics,
!    Volume 20, Number 3, 1971, pages 338-345.
!
!  Input:
!
!    real ( kind = rk ) X(3), the argument.
!
!  Output:
!
!    real ( kind = rk ) HELICAL, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) fx
  real ( kind = rk ) fx1
  real ( kind = rk ) fx2
  real ( kind = rk ) fx3
  real ( kind = rk ) helical
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) theta
  real ( kind = rk ) x(3)

  if ( 0.0D+00 < x(1) ) then
    theta = atan2 ( x(2), x(1) ) / 2.0D+00 / r8_pi
  else if ( x(1) < 0.0D+00 ) then
    theta = 0.5D+00 + atan2 ( x(2), x(1) ) / 2.0D+00 / r8_pi
  else if ( x(1) == 0.0D+00 ) then
    theta = 0.25D+00
  end if

  fx1 = x(3) - 10.0D+00 * theta
  fx2 = sqrt ( x(1) * x(1) + x(2) * x(2) )
  fx3 = x(3)

  fx = 100.0D+00 * fx1 * fx1 &
     +             fx2 * fx2 &
     +             fx3 * fx3

  helical = fx

  return
end
subroutine quartic_test ( )

!*****************************************************************************80
!
!! quartic_test() tests nelmin() on quartic().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 10

  integer i
  integer icount
  integer ifault
  integer kcount
  integer konvge
  integer numres
  real ( kind = rk ), external :: quartic
  real ( kind = rk ) reqmin
  real ( kind = rk ) start(n)
  real ( kind = rk ) step(n)
  real ( kind = rk ) xmin(n)
  real ( kind = rk ) ynewlo

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'quartic_test():'
  write ( *, '(a)' ) '  Apply nelmin() to the QUARTIC function.'

  start(1:n) = 1.0D+00

  reqmin = 1.0D-08

  step(1:n) = 1.0D+00

  konvge = 10
  kcount = 2000

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Starting point X:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,g14.6)' ) start(i)
  end do

  ynewlo = quartic ( start )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X) = ', ynewlo

  call nelmin ( quartic, n, start, xmin, ynewlo, reqmin, step, &
    konvge, kcount, icount, numres, ifault )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return code IFAULT = ', ifault
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Estimate of minimizing value X*:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,g14.6)' ) xmin(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  F(X*) = ', ynewlo

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of iterations = ', icount
  write ( *, '(a,i8)' ) '  Number of restarts =   ', numres

  return
end
function quartic ( x )

!*****************************************************************************80
!
!! quartic() evaluates a function defined by a sum of fourth powers.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 August 2021
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    R ONeill,
!    Algorithm AS 47:
!    Function Minimization Using a Simplex Procedure,
!    Applied Statistics,
!    Volume 20, Number 3, 1971, pages 338-345.
!
!  Input:
!
!    real ( kind = rk ) X(10), the argument.
!
!  Output:
!
!    real ( kind = rk ) QUARTIC, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) quartic
  real ( kind = rk ) x(10)

  quartic = sum ( x(1:10) ** 4 )

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
!    27 August 2021
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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
 
