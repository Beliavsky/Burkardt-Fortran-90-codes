program main

!*****************************************************************************80
!
!! lambert_w_test() tests lambert_w().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 June 2023
!
!  Author:
!
!    John Burkardt.
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'lambert_w_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test lambert_w().'

  call lambert_w_test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'lambert_w_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine lambert_w_test01 ( )

!*****************************************************************************80
!
!! lambert_w_test01() compares lambert_w() to stored values, and a build in function.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 June 2023
!
!  Author:
!
!    John Burkardt.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer branch
  real ( kind = rk ) fx1
  real ( kind = rk ) fx2
  real ( kind = rk ) lambert_w
  integer n_data
  real ( kind = rk ) x

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'lambert_w_test01():'
  write ( *, '(a)' ) '  Compare stored, computed, and system values for LambertW(x).'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '      X           Stored value       lambert_w(x)'
  write ( *, '(a)' ) ''

  n_data = 0

  do

    call lambert_w_values ( n_data, x, fx1, branch )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = lambert_w ( x, branch, 0 )

    write ( *, '(2x,g12.6,2x,f25.16,2x,f25.16)' ) x, fx1, fx2

  end do

  return
end
subroutine lambert_w_values ( n_data, x, fx, b )

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
!    The function is defined for -1/e <= x.
!
!    There are two branches, joining at -1/e = x.
!    The lower branch extends from -1/e <= x < 0
!    The upper branch extends from -1/e <= x
!
!    The function is also known as the "Omega" function.
!
!    In Mathematica, the function can be evaluated by:
!      W = ProductLog [ X ]
!
!    In MATLAB,
!      W = lambertw ( b, x )
!    where b = -1 for lower branch, 0 for upper branch.
!
!    In Python,
!      W = scipy.special.lambertw ( x, b )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 June 2023
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
!  Input:
!
!    integer n_data: The user sets N_DATA to 0 before the first call.
!
!  Output:
!
!    integer n_data: On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    real ( kind = rk ) x: the argument of the function.
!
!    real ( kind = rk ) fx: the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 41

  integer b
  integer, save, dimension ( n_max ) :: branch_vec = (/ &
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
     0 /)
  real ( kind = rk ) fx
  real ( kind = rk ), save, dimension ( n_max ) :: fx_vec = (/ &
   -4.889720169867429D+00, &
   -3.994308347002122D+00, &
   -3.439216483280204D+00, &
   -3.022313245324657D+00, &
   -2.678346990016661D+00, &
   -2.376421342062887D+00, &
   -2.097349210703492D+00, &
   -1.824388309032984D+00, &
   -1.531811608389612D+00, &
   -1.000000000000000D+00, &
   -0.608341284733432D+00, &
   -0.471671909743522D+00, &
   -0.374493134019498D+00, &
   -0.297083462446424D+00, &
   -0.231960952986534D+00, &
   -0.175356500529299D+00, &
   -0.125066982982524D+00, &
   -0.079678160511477D+00, &
   -0.038221241746799D+00, &
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
   -0.036787944117144D+00, &
   -0.073575888234288D+00, &
   -0.110363832351433D+00, &
   -0.147151776468577D+00, &
   -0.183939720585721D+00, &
   -0.220727664702865D+00, &
   -0.257515608820010D+00, &
   -0.294303552937154D+00, &
   -0.331091497054298D+00, &
   -0.367879441171442D+00, &
   -0.331091497054298D+00, &
   -0.294303552937154D+00, &
   -0.257515608820010D+00, &
   -0.220727664702865D+00, &
   -0.183939720585721D+00, &
   -0.147151776468577D+00, &
   -0.110363832351433D+00, &
   -0.073575888234288D+00, &
   -0.036787944117144D+00, &
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
    b = 0
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
    b = branch_vec(n_data)
  end if

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

