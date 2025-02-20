program main

!*****************************************************************************80
!
!! zero_muller_test() tests zero_muller().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 March 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'zero_muller_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test zero_muller(), which uses Muller''s method,'
  write ( *, '(a)' ) '  with complex arithmeic, to solve a nonlinear equation.'

  call test01 ( )
  call test02 ( )
  call test03 ()
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'zero_muller_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! test01() tests zero_muller() on F(X) = X*X+9.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )
  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  real ( kind = rk ) fatol
  external func01
  complex ( kind = ck ) fxnew
  integer itmax
  complex ( kind = ck ) x1
  complex ( kind = ck ) x2
  complex ( kind = ck ) x3
  complex ( kind = ck ) xnew
  real ( kind = rk ) xatol
  real ( kind = rk ) xrtol

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test01():'
  write ( *, '(a)' ) '  Demonstrate zero_muller() on F(X) = X*X+9.'

  fatol = 1.0D-05
  itmax = 10
  x1 = cmplx ( 1.0D+00, 0.0D+00, kind = ck )
  x2 = cmplx ( 0.0D+00, 1.0D+00, kind = ck )
  x3 = cmplx ( 0.5D+00, 0.5D+00, kind = ck )
  xatol = 1.0D-05
  xrtol = 1.0D-05

  call zero_muller ( func01, fatol, itmax, x1, x2, x3, xatol, xrtol, &
    xnew, fxnew )

  write ( *, '(a)' ) ' '
  write ( *, '(a,f20.10,f20.10)' ) '     X   = ', xnew
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  with function value F(X):'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f20.10,f20.10)' ) '    FX   = ', fxnew
  write ( *, '(a,f20.10)' )        '  ||FX|| = ', abs ( fxnew )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! test02() tests zero_muller() on F(X) = (X*X+4) * (X-10) * (X+20)
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )
  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  real ( kind = rk ) fatol
  external func02
  complex ( kind = ck ) fxnew
  integer itmax
  complex ( kind = ck ) x1
  complex ( kind = ck ) x2
  complex ( kind = ck ) x3
  complex ( kind = ck ) xnew
  real ( kind = rk ) xatol
  real ( kind = rk ) xrtol

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test02():'
  write ( *, '(a)' ) '  Demonstrate zero_muller() on F(X) = (X*X+4)*(X-10)*(X+20).'

  fatol = 1.0D-05
  itmax = 10
  x1 = cmplx ( 1.0D+00, 0.0D+00, kind = ck )
  x2 = cmplx ( 0.0D+00, 1.0D+00, kind = ck )
  x3 = cmplx ( 0.5D+00, 0.5D+00, kind = ck )
  xatol = 1.0D-05
  xrtol = 1.0D-05

  call zero_muller ( func02, fatol, itmax, x1, x2, x3, xatol, xrtol, &
    xnew, fxnew )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a,f20.10,f20.10)' ) '     X   = ', xnew
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  with function value F(X):'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f20.10,f20.10)' ) '    FX   = ', fxnew
  write ( *, '(a,f20.10)' )        '  ||FX|| = ', abs ( fxnew )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! test03() tests zero_muller() on Zhelyazkov's function
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )
  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  real ( kind = rk ) fatol
  external func03
  complex ( kind = ck ) fxnew
  integer itmax
  integer test
  complex ( kind = ck ) x1
  complex ( kind = ck ) x2
  complex ( kind = ck ) x3
  complex ( kind = ck ) xnew
  real ( kind = rk ) xatol
  real ( kind = rk ) xrtol

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test03():'
  write ( *, '(a)' ) '  Demonstrate zero_muller() on Zhelyazkov''s function.'

  fatol = 1.0D-07
  itmax = 10

  do test = 1, 2
!
!  First set of starting points.
!  Result is X = ( 1.5705798926, 0.0 )
!
    if ( test == 1 ) then
      x1 = cmplx ( 1.0D+00, 0.0D+00, kind = ck )
      x2 = cmplx ( 0.0D+00, 1.0D+00, kind = ck )
      x3 = cmplx ( 0.5D+00, 0.5D+00, kind = ck )
!
!  Second set of starting points.
!  Result is X = ( -0.5802520567, 0.0 ).
!
    else if ( test == 2 ) then
      x1 = cmplx (  0.0D+00, 1.0D+00, kind = ck )
      x2 = cmplx (  1.0D+00, 2.0D+00, kind = ck )
      x3 = cmplx ( -1.0D+00, 2.0D+00, kind = ck )
    end if

    xatol = 1.0D-07
    xrtol = 1.0D-07

    call zero_muller ( func03, fatol, itmax, x1, x2, x3, xatol, xrtol, &
      xnew, fxnew )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) ' '
    write ( *, '(a,f20.10,f20.10)' ) '     X   = ', xnew
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  with function value F(X):'
    write ( *, '(a)' ) ' '
    write ( *, '(a,f20.10,f20.10)' ) '    FX   = ', fxnew
    write ( *, '(a,f20.10)' )        '  ||FX|| = ', abs ( fxnew )

  end do

  return
end
subroutine func01 ( x, fx )

!*****************************************************************************80
!
!! func01() evaluates F(X) = X*X+9.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) X, the point at which the function is to
!    be evaluated.
!
!  Output:
!
!    complex ( kind = ck ) FX, the function value at X.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  complex ( kind = ck ) fx
  complex ( kind = ck ) x

  fx = x * x + 9.0D+00
 
  return
end
subroutine func02 ( x, fx )

!*****************************************************************************80
!
!! func02() evaluates F(X) = (X*X+4)*(X-1)*(X+2).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) X, the point at which the function is to
!    be evaluated.
!
!  Output:
!
!    complex ( kind = ck ) FX, the function value at X.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  complex ( kind = ck ) fx
  complex ( kind = ck ) x

  fx = ( x * x + 4.0D+00 ) * ( x - 10.0D+00 ) * ( x + 20.0D+00 )
 
  return
end
subroutine func03 ( z, fz )

!*****************************************************************************80
!
!! func03() evaluates Zhelyazkov's function.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) Z, the point at which the function is to
!    be evaluated.
!
!  Output:
!
!    complex ( kind = ck ) FZ, the function value at Z.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )
  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  complex ( kind = ck ) a
  complex ( kind = ck ) b
  complex ( kind = ck ), parameter :: eps = ( 0.4D+00, 0.0D+00 )
  complex ( kind = ck ), parameter :: eta = ( 0.64D+00, 0.0D+00 )
  complex ( kind = ck ) fz
  real ( kind = rk ), parameter :: me = 0.384D+00
  real ( kind = rk ), parameter :: mo = 0.5D+00
  complex ( kind = ck ) of
  complex ( kind = ck ) ok
  complex ( kind = ck ), parameter :: one = ( 1.0D+00, 0.0D+00 )
  real ( kind = rk ), parameter :: x = 0.5D+00
  complex ( kind = ck ) z

  ok = z - me / sqrt ( eta )
  of = z - mo

  a = of * of + ( ok * ok ) * eta * dtanh ( x )

  b = ( of - ok * eta ) / ( of - ok * eta * eta )

  fz = of * of - one + ( eta * ok * ok - one ) * &
    dtanh ( x ) - x * x * eps * eps * a * b

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

