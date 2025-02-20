program main

!*****************************************************************************80
!
!! poisson_1d_test() tests poisson_1d().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 October 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'poisson_1d_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test poisson_1d().'

  call test01 ( ) 
  call test02 ( ) 
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'poisson_1d_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( ) 

!*****************************************************************************80
!
!! test01() tests poisson_1d() on test case #1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  real ( kind = rk8 ) difmax
  real ( kind = rk8 ), external :: exact1
  real ( kind = rk8 ), external :: force1
  integer i
  integer it_num
  integer k
  integer n
  real ( kind = rk8 ), allocatable :: u(:)
  real ( kind = rk8 ) ua
  real ( kind = rk8 ) ub
  real ( kind = rk8 ), allocatable :: x(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test01():'
  write ( *, '(a)' ) '  poisson_1d_monogrid() solves a 1D Poisson BVP'
  write ( *, '(a)' ) '  using the Gauss-Seidel method.'

  a = 0.0D+00
  b = 1.0D+00
  ua = 0.0D+00
  ub = 0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  -u"(x) = 1, for 0 < x < 1'
  write ( *, '(a)' ) '  u(0) = u(1) = 0.'
  write ( *, '(a)' ) '  Solution is u(x) = ( -x^2 + x ) / 2'

  do k = 5, 5

    n = 2**k
    allocate ( u(1:n+1) )
    allocate ( x(1:n+1) )
    do i = 1, n + 1
      x(i) = ( ( n + 1 - i ) * a + ( i - 1 ) * b ) / n
    end do

    write ( *, '(a)' ) ' '
    write ( *, '(a,i4)' ) '  Mesh index K = ', k
    write ( *, '(a,i6)' ) '  Number of intervals N=2^K = ', n
    write ( *, '(a,i6)' ) '  Number of nodes = 2^K+1 =   ', n + 1

    call poisson_1d ( n, a, b, ua, ub, force1, it_num, u )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '     I        X(I)      U(I)         U Exact(X(I))'
    write ( *, '(a)' ) ' '
    do i = 1, n + 1
      write ( *, '(2x,i4,2x,f10.4,2x,g14.6,2x,g14.6)' ) &
        i, x(i), u(i), exact1 ( x(i) )
    end do

    write ( *, '(a)' ) ' '

    difmax = 0.0D+00
    do i = 1, n + 1
      difmax = max ( difmax, abs ( u(i) - exact1 ( x(i) ) ) )
    end do 
    write ( *, '(a,g14.6)' ) '  Maximum error = ', difmax
    write ( *, '(a,i6)' ) '  Number of Gauss-Seidel iterations = ', it_num

    deallocate ( u )
    deallocate ( x )

  end do

  return
end
function exact1 ( x )

!*****************************************************************************80
!                                                    
!! exact1() evaluates the exact solution for test case #1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Hager,
!    Applied Numerical Linear Algebra,
!    Prentice-Hall, 1988,
!    ISBN13: 978-0130412942,
!    LC: QA184.H33.
!
!  Input:
!
!    real ( kind = rk8 ) X, the evaluation point.
!
!  Parameters:
!
!    real ( kind = rk8 ) EXACT1, the value of the exact solution at X.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) exact1
  real ( kind = rk8 ) x

  exact1 = 0.5D+00 * ( - x * x + x )

  return
end
function force1 ( x )

!*****************************************************************************80
!                                                    
!! force1() evaluates the forcing function for test case #1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Hager,
!    Applied Numerical Linear Algebra,
!    Prentice-Hall, 1988,
!    ISBN13: 978-0130412942,
!    LC: QA184.H33.
!
!  Input:
!
!    real ( kind = rk8 ) X, the evaluation point.
!
!  Output:
!
!    real ( kind = rk8 ) FORCE1, the value of the forcing function at X.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) force1
  real ( kind = rk8 ) x

  call r8_fake_use ( x )

  force1 = 1.0D+00

  return
end
subroutine test02 ( ) 

!*****************************************************************************80
!
!! test02() tests poisson_1d() on test case #2.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  real ( kind = rk8 ) difmax
  real ( kind = rk8 ), external :: exact2
  real ( kind = rk8 ), external :: force2
  integer i
  integer it_num
  integer k
  integer n
  real ( kind = rk8 ), allocatable :: u(:)
  real ( kind = rk8 ) ua
  real ( kind = rk8 ) ub
  real ( kind = rk8 ), allocatable :: x(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test02():'
  write ( *, '(a)' ) '  poisson_1d() solves a 1D Poisson BVP'
  write ( *, '(a)' ) '  using the Gauss-Seidel method.'

  a = 0.0D+00
  b = 1.0D+00
  ua = 0.0D+00
  ub = 0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  -u"(x) = - x * (x+3) * exp(x), for 0 < x < 1'
  write ( *, '(a)' ) '  u(0) = u(1) = 0.'
  write ( *, '(a)' ) '  Solution is u(x) = x * (x-1) * exp(x)'

  do k = 5, 5

    n = 2**k
    allocate ( u(1:n+1) )
    allocate ( x(1:n+1) )
    do i = 1, n + 1
      x(i) = ( ( n + 1 - i ) * a + ( i - 1 ) * b ) / n
    end do

    write ( *, '(a)' ) ' '
    write ( *, '(a,i4)' ) '  Mesh index K = ', k
    write ( *, '(a,i6)' ) '  Number of intervals N=2^K = ', n
    write ( *, '(a,i6)' ) '  Number of nodes = 2^K+1 =   ', n + 1

    call poisson_1d ( n, a, b, ua, ub, force2, it_num, u )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '     I        X(I)      U(I)         U Exact(X(I))'
    write ( *, '(a)' ) ' '
    do i = 1, n + 1
      write ( *, '(2x,i4,2x,f10.4,2x,g14.6,2x,g14.6)' ) &
        i, x(i), u(i), exact2 ( x(i) )
    end do

    write ( *, '(a)' ) ' '

    difmax = 0.0D+00
    do i = 1, n + 1
      difmax = max ( difmax, abs ( u(i) - exact2 ( x(i) ) ) )
    end do 
    write ( *, '(a,g14.6)' ) '  Maximum error = ', difmax
    write ( *, '(a,i6)' ) '  Number of Gauss-Seidel iterations = ', it_num

    deallocate ( u )
    deallocate ( x )

  end do

  return
end
function exact2 ( x )

!*****************************************************************************80
!                                                    
!! exact2() evaluates the exact solution for test case #2.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Hager,
!    Applied Numerical Linear Algebra,
!    Prentice-Hall, 1988,
!    ISBN13: 978-0130412942,
!    LC: QA184.H33.
!
!  Input:
!
!    real ( kind = rk8 ) X, the evaluation point.
!
!  Output:
!
!    real ( kind = rk8 ) EXACT2, the value of the exact solution at X.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) exact2
  real ( kind = rk8 ) x

  exact2 = x * ( x - 1.0D+00 ) * exp ( x )

  return
end
function force2 ( x )

!*****************************************************************************80
!                                                    
!! force2() evaluates the forcing function for test case #2.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Hager,
!    Applied Numerical Linear Algebra,
!    Prentice-Hall, 1988,
!    ISBN13: 978-0130412942,
!    LC: QA184.H33.
!
!  Input:
!
!    real ( kind = rk8 ) X, the evaluation point.
!
!  Output:
!
!    real ( kind = rk8 ) FORCE2, the value of the forcing function at X.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) force2
  real ( kind = rk8 ) x

  force2 =  - x * ( x + 3.0D+00 ) * exp ( x )

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
!    real ( kind = rk8 ) X, the variable to be "used".
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) x

  if ( x /= x ) then
    write ( *, '(a)' ) '  r8_fake_use(): variable is NAN.'
  end if

  return
end
