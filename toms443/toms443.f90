subroutine lambert_w_values ( n_data, x, fx )

!*****************************************************************************80
!
!! lambert_w_values() returns some values of the Lambert W function.
!
!  Discussion:
!
!    The function W(X) is defined implicitly by:
!
!      W(X) * e^W(X) = X
!
!    The function is also known as the "Omega" function.
!
!    In Mathematica, the function can be evaluated by:
!
!      W = ProductLog [ X ]
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Corless, Gaston Gonnet, David Hare, David Jeffrey,
!    Donald Knuth,
!    On the Lambert W Function,
!    Advances in Computational Mathematics,
!    Volume 5, 1996, pages 329-359.
!
!    Brian Hayes,
!    "Why W?",
!    The American Scientist,
!    Volume 93, March-April 2005, pages 104-108.
!
!    Eric Weisstein,
!    CRC Concise Encyclopedia of Mathematics,
!    CRC Press, 2002,
!    Second edition,
!    ISBN: 1584883472,
!    LC: QA5.W45.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = rk ) X, the argument of the function.
!
!    Output, real ( kind = rk ) FX, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 22

  real ( kind = rk ) fx
  real ( kind = rk ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.0000000000000000D+00, &
    0.3517337112491958D+00, &
    0.5671432904097839D+00, &
    0.7258613577662263D+00, &
    0.8526055020137255D+00, &
    0.9585863567287029D+00, &
    0.1000000000000000D+01, &
    0.1049908894964040D+01, &
    0.1130289326974136D+01, &
    0.1202167873197043D+01, &
    0.1267237814307435D+01, &
    0.1326724665242200D+01, &
    0.1381545379445041D+01, &
    0.1432404775898300D+01, &
    0.1479856830173851D+01, &
    0.1524345204984144D+01, &
    0.1566230953782388D+01, &
    0.1605811996320178D+01, &
    0.1745528002740699D+01, &
    0.3385630140290050D+01, &
    0.5249602852401596D+01, &
    0.1138335808614005D+02 /)
  integer n_data
  real ( kind = rk ) x
  real ( kind = rk ), save, dimension ( n_max ) :: x_vec = (/ &
    0.0000000000000000D+00, &
    0.5000000000000000D+00, &
    0.1000000000000000D+01, &
    0.1500000000000000D+01, &
    0.2000000000000000D+01, &
    0.2500000000000000D+01, &
    0.2718281828459045D+01, &
    0.3000000000000000D+01, &
    0.3500000000000000D+01, &
    0.4000000000000000D+01, &
    0.4500000000000000D+01, &
    0.5000000000000000D+01, &
    0.5500000000000000D+01, &
    0.6000000000000000D+01, &
    0.6500000000000000D+01, &
    0.7000000000000000D+01, &
    0.7500000000000000D+01, &
    0.8000000000000000D+01, &
    0.1000000000000000D+02, &
    0.1000000000000000D+03, &
    0.1000000000000000D+04, &
    0.1000000000000000D+07 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
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
function wew_a ( x, en )

!*****************************************************************************80
!
!! WEW_A estimates Lambert's W function.
!
!  Discussion:
!
!    For a given X, this routine estimates the solution W of Lambert's 
!    equation:
!
!      X = W * EXP ( W )
!
!    This routine has higher accuracy than WEW_B.
!
!  Modified:
!
!    08 June 2014
!
!  Reference:
!
!    Fred Fritsch, R Shafer, W Crowley,
!    Algorithm 443: Solution of the transcendental equation w e^w = x,
!    Communications of the ACM,
!    October 1973, Volume 16, Number 2, pages 123-124.
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, the argument of W(X)
!
!    Output, real ( kind = rk ) WEW_A, the estimated value of W(X).
!
!    Output, real ( kind = rk ) EN, the last relative correction to W(X).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), parameter :: c1 = 4.0D+00 / 3.0D+00
  real ( kind = rk ), parameter :: c2 = 7.0D+00 / 3.0D+00
  real ( kind = rk ), parameter :: c3 = 5.0D+00 / 6.0D+00
  real ( kind = rk ), parameter :: c4 = 2.0D+00 / 3.0D+00
  real ( kind = rk ) en
  real ( kind = rk ) f
  real ( kind = rk ) temp
  real ( kind = rk ) temp2
  real ( kind = rk ) wew_a
  real ( kind = rk ) wn
  real ( kind = rk ) x
  real ( kind = rk ) y
  real ( kind = rk ) zn
!
!  Initial guess.
!
  f = log ( x )

  if ( x <= 6.46D+00 ) then

    wn = x * ( 1.0D+00 + c1 * x ) / ( 1.0D+00 + x * ( c2 + c3 * x ) )
    zn = f - wn - log ( wn )

  else

    wn = f
    zn = - log ( wn )

  end if
!
!  Iteration 1.
!
  temp = 1.0D+00 + wn
  y = 2.0D+00 * temp * ( temp + c4 * zn ) - zn
  wn = wn * ( 1.0D+00 + zn * y / ( temp * ( y - zn ) ) )
!
!  Iteration 2.
!
  zn = f - wn - log ( wn )
  temp = 1.0D+00 + wn
  temp2 = temp + c4 * zn
  en = zn * temp2 / ( temp * temp2 - 0.5D+00 * zn )
  wn = wn * ( 1.0D+00 + en )

  wew_a = wn

  return
end
function wew_b ( x, en )

!*****************************************************************************80
!
!! WEW_B estimates Lambert's W function.
!
!  Discussion:
!
!    For a given X, this routine estimates the solution W of Lambert's 
!    equation:
!
!      X = W * EXP ( W )
!
!    This routine has lower accuracy than WEW_A.
!
!  Modified:
!
!    08 June 2014
!
!  Reference:
!
!    Fred Fritsch, R Shafer, W Crowley,
!    Algorithm 443: Solution of the transcendental equation w e^w = x,
!    Communications of the ACM,
!    October 1973, Volume 16, Number 2, pages 123-124.
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, the argument of W(X)
!
!    Output, real ( kind = rk ) WEW_B, the estimated value of W(X).
!
!    Output, real ( kind = rk ) EN, the last relative correction to W(X).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), parameter :: c1 = 4.0D+00 / 3.0D+00
  real ( kind = rk ), parameter :: c2 = 7.0D+00 / 3.0D+00
  real ( kind = rk ), parameter :: c3 = 5.0D+00 / 6.0D+00
  real ( kind = rk ), parameter :: c4 = 2.0D+00 / 3.0D+00
  real ( kind = rk ) en
  real ( kind = rk ) f
  real ( kind = rk ) temp
  real ( kind = rk ) wew_b
  real ( kind = rk ) wn
  real ( kind = rk ) x
  real ( kind = rk ) y
  real ( kind = rk ) zn
!
!  Initial guess.
!
  f = log ( x )

  if ( x <= 0.7385D+00 ) then
    wn = x * ( 1.0D+00 + c1 * x ) / ( 1.0D+00 + x * ( c2 + c3 * x ) )
  else
    wn = f - 24.0D+00 * ( ( f + 2.0D+00 ) * f - 3.0D+00 ) &
      / ( ( 0.7D+00 * f + 58.0D+00 ) * f + 127.0D+00 )
  end if
!
!  Iteration 1.
!
  zn = f - wn - log ( wn )
  temp = 1.0D+00 + wn
  y = 2.0D+00 * temp * ( temp + c4 * zn ) - zn
  en = zn * y / ( temp * ( y - zn ) )
  wn = wn * ( 1.0D+00 + en )

  wew_b = wn

  return
end
