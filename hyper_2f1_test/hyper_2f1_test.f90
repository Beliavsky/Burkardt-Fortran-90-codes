program main

!*****************************************************************************80
!
!! hyper_2f1_test() tests hyper_2f1().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 December 2023
!
!  Author:
!
!    John Burkardt
!
  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'hyper_2f1_test():'
  write ( *, '(a)' ) '  Fortran90 version.'
  write ( *, '(a)' ) '  Test hyper_2f1().'

  call hyper_2f1_real_test ( )
  call hyper_2f1_complex_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'hyper_2f1_test()'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  return
end
subroutine hyper_2f1_real_test ( )

!*****************************************************************************80
!
!! hyper_2f1_real_test() tests hyper_2f1_real().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 December 2023
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  complex ( kind = ck ) ac
  real ( kind = rk ) b
  complex ( kind = ck ) bc
  real ( kind = rk ) c
  complex ( kind = ck ) cc
  real ( kind = rk ) fx1
  real ( kind = rk ) fx2
  real ( kind = rk ) hyper_2f1
  integer n_data
  real ( kind = rk ) x
  complex ( kind = ck ) z
  
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'hyper_2f1_real_test():'
  write ( *, '(a)' ) '  Test hyper_2f1().'

  n_data = 0

  do while ( .true. )

    call hyper_2f1_real_values ( n_data, a, b, c, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    ac = a
    bc = b
    cc = c
    z = x
    fx2 = hyper_2f1 ( ac, bc, cc, z )

    write ( *, '(a)' ) ''
    write ( *, '(a,f8.4,2x,f8.4,2x,f8.4,2x,f8.4)' ) &
      '  A,B,C,X:      ', a, b, c, x
    write ( *, '(a,g14.6)' ) '  FX exact:     ', fx1
    write ( *, '(a,g14.6)' ) '  FX computed:  ', fx2
    write ( *, '(a,g14.6)' ) '  Error:        ', abs ( fx1 - fx2 )

  end do


  return
end
subroutine hyper_2f1_real_values ( n_data, a, b, c, x, fx )

!*****************************************************************************80
!
!! hyper_2f1_real_values() returns some values of the hypergeometric function 2F1.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      fx = Hypergeometric2F1 [ a, b, c, x ]
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 September 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996,
!    ISBN: 0-8493-2479-3,
!    LC: QA47.M315.
!
!  Input:
!
!    integer N_DATA: The user sets N_DATA to 0 before the first call.
!
!  Output:
!
!    integer N_DATA: On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    real ( kind = rk ) A, B, C: the parameters.
!
!    real ( kind = rk ) X: the argument.
!
!    real ( kind = rk ) FX, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 24

  real ( kind = rk ) a
  real ( kind = rk ), save, dimension ( n_max ) :: a_vec = (/ &
   -2.5D+00, &
   -0.5D+00, &
    0.5D+00, &
    2.5D+00, &
   -2.5D+00, &
   -0.5D+00, &
    0.5D+00, &
    2.5D+00, &
   -2.5D+00, &
   -0.5D+00, &
    0.5D+00, &
    2.5D+00, &
    3.3D+00, &
    1.1D+00, &
    1.1D+00, &
    3.3D+00, &
    3.3D+00, &
    1.1D+00, &
    1.1D+00, &
    3.3D+00, &
    3.3D+00, &
    1.1D+00, &
    1.1D+00, &
    3.3D+00 /)
  real ( kind = rk ) b
  real ( kind = rk ), save, dimension ( n_max ) :: b_vec = (/ &
    3.3D+00, &
    1.1D+00, &
    1.1D+00, &
    3.3D+00, &
    3.3D+00, &
    1.1D+00, &
    1.1D+00, &
    3.3D+00, &
    3.3D+00, &
    1.1D+00, &
    1.1D+00, &
    3.3D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00 /)
  real ( kind = rk ) c
  real ( kind = rk ), save, dimension ( n_max ) :: c_vec = (/ &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
    6.7D+00, &
   -5.5D+00, &
   -0.5D+00, &
    0.5D+00, &
    4.5D+00, &
   -5.5D+00, &
   -0.5D+00, &
    0.5D+00, &
    4.5D+00, &
   -5.5D+00, &
   -0.5D+00, &
    0.5D+00, &
    4.5D+00 /)
  real ( kind = rk ) fx
  real ( kind = rk ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.72356129348997784913D+00, &
    0.97911109345277961340D+00, &
    1.0216578140088564160D+00, &
    1.4051563200112126405D+00, &
    0.46961431639821611095D+00, &
    0.95296194977446325454D+00, &
    1.0512814213947987916D+00, &
    2.3999062904777858999D+00, &
    0.29106095928414718320D+00, &
    0.92536967910373175753D+00, &
    1.0865504094806997287D+00, &
    5.7381565526189046578D+00, &
    15090.669748704606754D+00, &
   -104.31170067364349677D+00, &
    21.175050707768812938D+00, &
    4.1946915819031922850D+00, &
    1.0170777974048815592D+10, &
   -24708.635322489155868D+00, &
    1372.2304548384989560D+00, &
    58.092728706394652211D+00, &
    5.8682087615124176162D+18, &
   -4.4635010147295996680D+08, &
    5.3835057561295731310D+06, &
    20396.913776019659426D+00 /)
  integer n_data
  real ( kind = rk ) x
  real ( kind = rk ), save, dimension ( n_max ) :: x_vec = (/ &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.55D+00, &
    0.55D+00, &
    0.55D+00, &
    0.55D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.55D+00, &
    0.55D+00, &
    0.55D+00, &
    0.55D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00, &
    0.85D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0D+00
    b = 0.0D+00
    c = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    b = b_vec(n_data)
    c = c_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine hyper_2f1_complex_test ( )

!*****************************************************************************80
!
!! hyper_2f1_complex_test() tests hyper_2f1_complex().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 December 2023
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  complex ( kind = ck ) ac
  real ( kind = rk ) b
  complex ( kind = ck ) bc
  real ( kind = rk ) c
  complex ( kind = ck ) cc
  complex ( kind = ck ) fz1
  complex ( kind = ck ) fz2
  complex ( kind = ck ) hyper_2f1
  integer n_data
  complex ( kind = ck ) z
  
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'hyper_2f1_complex_test():'
  write ( *, '(a)' ) '  Test hyper_2f1().'

  n_data = 0

  do while ( .true. )

    call hyper_2f1_complex_values ( n_data, a, b, c, z, fz1 )

    if ( n_data == 0 ) then
      exit
    end if

    ac = a
    bc = b
    cc = c
    fz2 = hyper_2f1 ( ac, bc, cc, z )

    write ( *, '(a)' ) ''
    write ( *, '(a,f8.4,2x,f8.4,2x,f8.4,2x,f8.4,f8.4)' ) &
      '  A,B,C,Z:      ', a, b, c, z
    write ( *, '(a,g24.16,g24.16)' ) '  FX exact:     ', fz1
    write ( *, '(a,g24.16,g24.16)' ) '  FX computed:  ', fz2
    write ( *, '(a,g14.6)' ) '  Error:        ', abs ( fz1 - fz2 )

  end do

  return
end
subroutine hyper_2f1_complex_values ( n_data, a, b, c, z, fz )

!*****************************************************************************80
!
!! hyper_2f1_complex_values() returns some values of the hypergeometric function 2F1.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      fz = Hypergeometric2F1 [ a, b, c, z ]
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 December 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996,
!    ISBN: 0-8493-2479-3,
!    LC: QA47.M315.
!
!  Input:
!
!    integer N_DATA: The user sets N_DATA to 0 before the first call.
!
!  Output:
!
!    integer N_DATA: On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    real ( kind = rk ) A, B, C: the parameters.
!
!    complex ( kind = ck ) Z: the argument.
!
!    complex ( kind = ck ) FZ: the value of the function.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 15

  real ( kind = rk ) a
  real ( kind = rk ), save, dimension ( n_max ) :: a_vec = (/ &
    3.2, &
    3.2, &
   -5.0, &
    3.3, &
   -7.0, &
    4.3, &
    3.3, &
    3.5, &
    3.3, &
    7.0, &
    5.0, &
    3.5, &
    2.1, &
    8.7, &
    8.7 /)
  real ( kind = rk ) b
  real ( kind = rk ), save, dimension ( n_max ) :: b_vec = (/ &
   1.8, &
  -1.8, &
   3.3, &
  -6.0, &
   3.3, &
  -8.0, &
   5.8, &
  -2.4, &
   4.3, &
   5.0, &
   7.0, &
   1.2, &
   5.4, &
   3.2, &
   2.7 /)
  real ( kind = rk ) c
  real ( kind = rk ), save, dimension ( n_max ) :: c_vec = (/ &
   6.7, &
   6.7, &
   6.7, &
   3.7, &
  -3.7, &
  -3.7, &
   6.7, &
   6.7, &
   6.7, &
   4.1, &
   4.1, &
   9.7, &
   9.7, &
   6.7, &
   6.7 /)
  complex ( kind = ck ) fz
  complex ( kind = ck ), save, dimension ( n_max ) :: fz_vec = (/ &
   (  5.468999154361234D+00,    +0.00000000D+00 ), &
   (  0.3375063477462785D+00,   +0.00000000D+00 ), &
   (  116.8274991533609D+00,    +603.8909562709345D+00 ), &
   (  17620.41819334182D+00,    +38293.80901310932D+00 ), &
   ( -11772775115.27448D+00,    -14382285977.20268D+00 ), &
   (  1316118577866.058D+00,    -101298889382.4362D+00 ), &
   (  1.733055678355656D+00,    +0.6340102904953357D+00 ), &
   (  0.6476224071999852D+00,   -0.5211050690999773D+00 ), &
   ( -1.483008322270093D+00,    +8.374426179451589D+00 ), &
   ( -0.004037609523971226D+00, -0.002956632645480181D+00 ), &
   ( -0.004037609523971226D+00, -0.002956632645480181D+00 ), &
   (  1.034313610729953D+00,    +0.5447389238499308D+00 ), &
   (  0.6885043978280027D+00,   +1.227418679098749D+00 ), &
   ( -0.9004649679297319D+00,   -1.11988994714304D+00 ), &
   ( -0.4608388640599718D+00,   -0.5457569650549665D+00 ) /)
  integer n_data
  complex ( kind = ck ) z
  complex ( kind = ck ), save, dimension ( n_max ) :: z_vec = (/ &
   ( 1.0, +0.0 ), &
   ( 1.0, +0.0 ), &
   ( 5.2, +4.8 ), &
   ( 5.2, -4.8 ), &
   ( 5.2, -4.8 ), &
   ( 5.2, +4.8 ), &
   ( 0.2, +0.1 ), &
   ( 0.2, +0.5 ), &
   ( 0.8, +0.3 ), &
   ( 3.0, -1.0 ), &
   ( 3.0, -1.0 ), &
   ( 0.6, +0.9 ), &
   ( 0.5, +0.7 ), &
   ( 0.5, +0.7 ), &
   ( 0.6, +0.9 ) /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0D+00
    b = 0.0D+00
    c = 0.0D+00
    z = 0.0D+00
    fz = ( 0.0D+00, 0.0D+00 )
  else
    a = a_vec(n_data)
    b = b_vec(n_data)
    c = c_vec(n_data)
    z = z_vec(n_data)
    fz = fz_vec(n_data)
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
