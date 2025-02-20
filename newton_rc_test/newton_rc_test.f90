program main

!*****************************************************************************80
!
!! newton_rc_test tests newton_rc().
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

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'newton_rc_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  newton_rc() solves nonlinear equations'
  write ( *, '(a)' ) '  using reverse communication.'

  call newton_rc_test01 ( )
  call newton_rc_test02 ( )
  call newton_rc_test03 ( )
  call newton_rc_test04 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'newton_rc_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine newton_rc_test01 ( )

!*****************************************************************************80
!
!! newton_rc_test01 calls newton_rc() for 1 equation.
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

  real ( kind = rk ) fx(n)
  integer ido
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'newton_rc_test01'
  write ( *, '(a)' ) '  Use newton_rc() to solve a system of 1 nonlinear equation.'
!
!  Initialization.
!
  ido = 0
  x(1:n) = 0.0D+00
  call f1 ( n, x, fx )

  call r8vec2_print ( n, x, fx, '  Initial X and F(X)' )

  do

    call newton_rc ( ido, n, x, fx )

    if ( ido == 0 ) then
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) '  Convergence:'
      exit
    else if ( ido == 1 .or. ido == 2 ) then
      call f1 ( n, x, fx )
    else
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) '  Convergence failure:'
      exit
    end if

  end do

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
subroutine newton_rc_test02 ( )

!*****************************************************************************80
!
!! newton_rc_test02 calls newton_rc() for 2 equations.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 April 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 2

  real ( kind = rk ) fx(n)
  integer ido
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'newton_rc_test02'
  write ( *, '(a)' ) '  Use newton_rc() to solve a system of 2 nonlinear equations.'
!
!  Initialization.
!
  ido = 0
  x(1) = 3.0D+00
  x(2) = 0.0D+00
  call f2 ( n, x, fx )

  call r8vec2_print ( n, x, fx, '  Initial X and F(X)' )

  do

    call newton_rc ( ido, n, x, fx )

    if ( ido == 0 ) then
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) '  Convergence:'
      exit
    else if ( ido == 1 .or. ido == 2 ) then
      call f2 ( n, x, fx )
    else
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) '  Convergence failure:'
      exit
    end if

  end do

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
subroutine newton_rc_test03 ( )

!*****************************************************************************80
!
!! newton_rc_test03 calls newton_rc() for 4 equations.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 April 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 4

  real ( kind = rk ) fx(n)
  integer ido
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'newton_rc_test03'
  write ( *, '(a)' ) '  Use newton_rc() to solve a system of 4 nonlinear equations.'
!
!  Initialization.
!
  ido = 0
  x(1:n) = 0.0D+00
  call f3 ( n, x, fx )

  call r8vec2_print ( n, x, fx, '  Initial X and F(X)' )

  do

    call newton_rc ( ido, n, x, fx )

    if ( ido == 0 ) then
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) '  Convergence:'
      exit
    else if ( ido == 1 .or. ido == 2 ) then
      call f3 ( n, x, fx )
    else
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) '  Convergence failure:'
      exit
    end if

  end do

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
subroutine newton_rc_test04 ( )

!*****************************************************************************80
!
!! newton_rc_test04 calls newton_rc() for 8 equations.
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

  real ( kind = rk ) fx(n)
  integer ido
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'newton_rc_test04'
  write ( *, '(a)' ) '  Use newton_rc() to solve a system of 8 nonlinear equations.'
!
!  Initialization.
!
  ido = 0
  x(1:n) = 0.0D+00
  call f4 ( n, x, fx )

  call r8vec2_print ( n, x, fx, '  Initial X and F(X)' )

  do

    call newton_rc ( ido, n, x, fx )

    if ( ido == 0 ) then
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) '  Convergence:'
      exit
    else if ( ido == 1 .or. ido == 2 ) then
      call f4 ( n, x, fx )
    else
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) '  Convergence failure:'
      exit
    end if

  end do

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

