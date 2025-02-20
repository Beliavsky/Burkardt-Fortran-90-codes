function rc ( x, y, errtol, ierr )

!*****************************************************************************80
!
!! rc() computes the elementary integral RC(X,Y).
!
!  Discussion:
!
!    This function computes the elementary integral
!
!      RC(X,Y) = Integral ( 0 <= T < oo )
!
!                  -1/2     -1
!        (1/2)(T+X)    (T+Y)  DT,
!
!    where X is nonnegative and Y is positive.  The duplication
!    theorem is iterated until the variables are nearly equal,
!    and the function is then expanded in Taylor series to fifth
!    order.
!
!    Logarithmic, inverse circular, and inverse hyperbolic
!    functions can be expressed in terms of RC.
!
!    Check by addition theorem:
!
!      RC(X,X+Z) + RC(Y,Y+Z) = RC(0,Z),
!      where X, Y, and Z are positive and X * Y = Z * Z.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 May 2018
!
!  Author:
!
!    Original FORTRAN77 version by Bille Carlson, Elaine Notis.
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Bille Carlson,
!    Computing Elliptic Integrals by Duplication,
!    Numerische Mathematik,
!    Volume 33, 1979, pages 1-16.
!
!    Bille Carlson, Elaine Notis,
!    Algorithm 577, Algorithms for Incomplete Elliptic Integrals,
!    ACM Transactions on Mathematical Software,
!    Volume 7, Number 3, pages 398-403, September 1981.
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, Y, the arguments in the integral.
!
!    Input, real ( kind = rk ) ERRTOL, the error tolerance.
!    Relative error due to truncation is less than
!      16 * ERRTOL ^ 6 / (1 - 2 * ERRTOL).
!    Sample choices:
!      ERRTOL   Relative truncation error less than
!      1.D-3    2.D-17
!      3.D-3    2.D-14
!      1.D-2    2.D-11
!      3.D-2    2.D-8
!      1.D-1    2.D-5
!
!    Output, integer IERR, the error flag.
!    0, no error occurred.
!    1, abnormal termination.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) c1
  real ( kind = rk ) c2
  real ( kind = rk ) errtol
  integer ierr
  real ( kind = rk ) lamda
  real ( kind = rk ) lolim
  real ( kind = rk ) mu
  real ( kind = rk ) rc
  real ( kind = rk ) s
  real ( kind = rk ) sn
  real ( kind = rk ) uplim
  real ( kind = rk ) x
  real ( kind = rk ) xn
  real ( kind = rk ) y
  real ( kind = rk ) yn
!
!  LOLIM AND UPLIM DETERMINE THE RANGE OF VALID ARGUMENTS.
!  LOLIM IS NOT LESS THAN THE MACHINE MINIMUM MULTIPLIED BY 5.
!  UPLIM IS NOT GREATER THAN THE MACHINE MAXIMUM DIVIDED BY 5.
!
  save lolim
  save uplim

  data lolim /3.D-78/
  data uplim /1.D+75/

  if ( &
    x < 0.0d0 .or. &
    y <= 0.0d0 .or. &
    ( x + y ) < lolim .or. &
    uplim < x .or. &
    uplim < y ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'RC - Error!'
    write ( *, '(a)' ) '  Invalid input arguments.'
    write ( *, '(a,d23.16)' ) '  X = ', x
    write ( *, '(a,d23.16)' ) '  Y = ', y
    write ( *, '(a)' ) ''
    ierr = 1
    rc = 0.0D+00
    return
  end if

  ierr = 0
  xn = x
  yn = y

  do

    mu = ( xn + yn + yn ) / 3.0d0
    sn = ( yn + mu ) / mu - 2.0d0

    if ( abs ( sn ) < errtol ) then
      c1 = 1.0d0 / 7.0d0
      c2 = 9.0d0 / 22.0d0
      s = sn * sn * ( 0.3d0 &
    + sn * ( c1 + sn * ( 0.375d0 + sn * c2 ) ) )
      rc = ( 1.0d0 + s ) / sqrt ( mu )
      return
    end if

    lamda = 2.0d0 * sqrt ( xn ) * sqrt ( yn ) + yn
    xn = ( xn + lamda ) * 0.25d0
    yn = ( yn + lamda ) * 0.25d0

  end do

end
function rd ( x, y, z, errtol, ierr )

!*****************************************************************************80
!
!! RD computes an incomplete elliptic integral of the second kind, RD(X,Y,Z).
!
!  Discussion:
!
!    This function computes an incomplete elliptic integral of the second kind.
!
!    RD(X,Y,Z) = Integral ( 0 <= T < oo )
!
!                    -1/2     -1/2     -3/2
!          (3/2)(T+X)    (T+Y)    (T+Z)    DT,
!
!    where X and Y are nonnegative, X + Y is positive, and Z is positive.
!
!    If X or Y is zero, the integral is complete.
!
!    The duplication theorem is iterated until the variables are
!    nearly equal, and the function is then expanded in Taylor
!    series to fifth order.
!
!    Check:
!
!      RD(X,Y,Z) + RD(Y,Z,X) + RD(Z,X,Y) = 3 / sqrt ( X * Y * Z ),
!      where X, Y, and Z are positive.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 May 2018
!
!  Author:
!
!    Original FORTRAN77 version by Bille Carlson, Elaine Notis.
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Bille Carlson,
!    Computing Elliptic Integrals by Duplication,
!    Numerische Mathematik,
!    Volume 33, 1979, pages 1-16.
!
!    Bille Carlson, Elaine Notis,
!    Algorithm 577, Algorithms for Incomplete Elliptic Integrals,
!    ACM Transactions on Mathematical Software,
!    Volume 7, Number 3, pages 398-403, September 1981.
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, Y, Z, the arguments in the integral.
!
!    Input, real ( kind = rk ) ERRTOL, the error tolerance.
!    The relative error due to truncation is less than
!      3 * ERRTOL ^ 6 / (1-ERRTOL) ^ 3/2.
!    Sample choices:
!      ERRTOL   Relative truncation error less than
!      1.D-3    4.D-18
!      3.D-3    3.D-15
!      1.D-2    4.D-12
!      3.D-2    3.D-9
!      1.D-1    4.D-6
!
!    Output, integer IERR, the error flag.
!    0, no error occurred.
!    1, abnormal termination.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) c1
  real ( kind = rk ) c2
  real ( kind = rk ) c3
  real ( kind = rk ) c4
  real ( kind = rk ) ea
  real ( kind = rk ) eb
  real ( kind = rk ) ec
  real ( kind = rk ) ed
  real ( kind = rk ) ef
  real ( kind = rk ) epslon
  real ( kind = rk ) errtol
  integer ierr
  real ( kind = rk ) lamda
  real ( kind = rk ) lolim
  real ( kind = rk ) mu
  real ( kind = rk ) power4
  real ( kind = rk ) rd
  real ( kind = rk ) sigma
  real ( kind = rk ) s1
  real ( kind = rk ) s2
  real ( kind = rk ) uplim
  real ( kind = rk ) x
  real ( kind = rk ) xn
  real ( kind = rk ) xndev
  real ( kind = rk ) xnroot
  real ( kind = rk ) y
  real ( kind = rk ) yn
  real ( kind = rk ) yndev
  real ( kind = rk ) ynroot
  real ( kind = rk ) z
  real ( kind = rk ) zn
  real ( kind = rk ) zndev
  real ( kind = rk ) znroot
!
!  LOLIM AND UPLIM DETERMINE THE RANGE OF VALID ARGUMENTS.
!  LOLIM IS NOT LESS THAN 2 / (MACHINE MAXIMUM) ^ (2/3).
!  UPLIM IS NOT GREATER THAN (0.1 * ERRTOL / MACHINE
!  MINIMUM) ^ (2/3), WHERE ERRTOL IS DESCRIBED BELOW.
!  IN THE FOLLOWING TABLE IT IS ASSUMED THAT ERRTOL WILL
!  NEVER BE CHOSEN SMALLER THAN 1.D-5.
!
  save lolim
  save uplim

  data lolim /6.D-51/
  data uplim /1.D+48/

  if ( &
    x < 0.0D+00 .or. &
    y < 0.0D+00 .or. &
    x + y < lolim .or. &
    z < lolim .or. &
    uplim < x .or. &
    uplim < y .or. &
    uplim < z ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'RD - Error!'
    write ( *, '(a)' ) '  Invalid input arguments.'
    write ( *, '(a,d23.16)' ) '  X = ', x
    write ( *, '(a,d23.16)' ) '  Y = ', y
    write ( *, '(a,d23.16)' ) '  Z = ', z
    write ( *, '(a)' ) ''
    ierr = 1
    rd = 0.0D+00
    return
  end if

  ierr = 0
  xn = x
  yn = y
  zn = z
  sigma = 0.0d0
  power4 = 1.0d0

  do

    mu = ( xn + yn + 3.0d0 * zn ) * 0.2d0
    xndev = ( mu - xn ) / mu
    yndev = ( mu - yn ) / mu
    zndev = ( mu - zn ) / mu
    epslon = max ( abs ( xndev ), abs ( yndev ), abs ( zndev ) )

    if ( epslon < errtol ) then
      c1 = 3.0d0 / 14.0d0
      c2 = 1.0d0 / 6.0d0
      c3 = 9.0d0 / 22.0d0
      c4 = 3.0d0 / 26.0d0
      ea = xndev * yndev
      eb = zndev * zndev
      ec = ea - eb
      ed = ea - 6.0d0 * eb
      ef = ed + ec + ec
      s1 = ed * ( - c1 + 0.25d0 * c3 * ed - 1.5d0 * c4 * zndev * ef )
      s2 = zndev  * ( c2 * ef + zndev * ( - c3 * ec + zndev * c4 * ea ) )
      rd = 3.0d0 * sigma  + power4 * ( 1.0d0 + s1 + s2 ) / ( mu * sqrt ( mu ) )

      return
    end if

    xnroot = sqrt ( xn )
    ynroot = sqrt ( yn )
    znroot = sqrt ( zn )
    lamda = xnroot * ( ynroot + znroot ) + ynroot * znroot
    sigma = sigma + power4 / ( znroot * ( zn + lamda ) )
    power4 = power4 * 0.25d0
    xn = ( xn + lamda ) * 0.25d0
    yn = ( yn + lamda ) * 0.25d0
    zn = ( zn + lamda ) * 0.25d0

  end do

end
function rf ( x, y, z, errtol, ierr )

!*****************************************************************************80
!
!! RF computes an incomplete elliptic integral of the first kind, RF(X,Y,Z).
!
!  Discussion:
!
!    This function computes the incomplete elliptic integral of the first kind.
!
!    RF(X,Y,Z) = Integral ( 0 <= T < oo )
!
!                    -1/2     -1/2     -1/2
!          (1/2)(T+X)    (T+Y)    (T+Z)    DT,
!
!    where X, Y, and Z are nonnegative and at most one of them is zero.
!
!    If X or Y or Z is zero, the integral is complete.
!
!    The duplication theorem is iterated until the variables are
!    nearly equal, and the function is then expanded in Taylor
!    series to fifth order.
!
!    Check by addition theorem:
!
!      RF(X,X+Z,X+W) + RF(Y,Y+Z,Y+W) = RF(0,Z,W),
!      where X, Y, Z, W are positive and X * Y = Z * W.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 May 2018
!
!  Author:
!
!    Original FORTRAN77 version by Bille Carlson, Elaine Notis.
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Bille Carlson,
!    Computing Elliptic Integrals by Duplication,
!    Numerische Mathematik,
!    Volume 33, 1979, pages 1-16.
!
!    Bille Carlson, Elaine Notis,
!    Algorithm 577, Algorithms for Incomplete Elliptic Integrals,
!    ACM Transactions on Mathematical Software,
!    Volume 7, Number 3, pages 398-403, September 1981.
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, Y, Z, the arguments in the integral.
!
!    Input, real ( kind = rk ) ERRTOL, the error tolerance.
!    Relative error due to truncation is less than
!      ERRTOL ^ 6 / (4 * (1 - ERRTOL)).
!    Sample choices:
!      ERRTOL   Relative truncation error less than
!      1.D-3    3.D-19
!      3.D-3    2.D-16
!      1.D-2    3.D-13
!      3.D-2    2.D-10
!      1.D-1    3.D-7
!
!    Output, integer IERR, the error flag.
!    0, no error occurred.
!    1, abnormal termination.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) c1
  real ( kind = rk ) c2
  real ( kind = rk ) c3
  real ( kind = rk ) e2
  real ( kind = rk ) e3
  real ( kind = rk ) epslon
  real ( kind = rk ) errtol
  integer ierr
  real ( kind = rk ) lamda
  real ( kind = rk ) lolim
  real ( kind = rk ) mu
  real ( kind = rk ) rf
  real ( kind = rk ) s
  real ( kind = rk ) uplim
  real ( kind = rk ) x
  real ( kind = rk ) xn
  real ( kind = rk ) xndev
  real ( kind = rk ) xnroot
  real ( kind = rk ) y
  real ( kind = rk ) yn
  real ( kind = rk ) yndev
  real ( kind = rk ) ynroot
  real ( kind = rk ) z
  real ( kind = rk ) zn
  real ( kind = rk ) zndev
  real ( kind = rk ) znroot
!
!  LOLIM AND UPLIM DETERMINE THE RANGE OF VALID ARGUMENTS.
!  LOLIM IS NOT LESS THAN THE MACHINE MINIMUM MULTIPLIED BY 5.
!  UPLIM IS NOT GREATER THAN THE MACHINE MAXIMUM DIVIDED BY 5.
!
  save lolim
  save uplim

  data lolim /3.D-78/
  data uplim /1.D+75/

  if ( &
    x < 0.0D+00 .or. &
    y < 0.0D+00 .or. &
    z < 0.0D+00 .or. &
    x + y < lolim .or. &
    x + z < lolim .or. &
    y + z < lolim .or. &
    uplim <= x .or. &
    uplim <= y .or. &
    uplim <= z ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'RF - Error!'
    write ( *, '(a)' ) '  Invalid input arguments.'
    write ( *, '(a,d23.16)' ) '  X = ', x
    write ( *, '(a,d23.16)' ) '  Y = ', y
    write ( *, '(a,d23.16)' ) '  Z = ', z
    write ( *, '(a)' ) ''
    ierr = 1
    rf = 0.0D+00
    return
  end if

  ierr = 0
  xn = x
  yn = y
  zn = z

  do

    mu = ( xn + yn + zn ) / 3.0d0
    xndev = 2.0d0 - ( mu + xn ) / mu
    yndev = 2.0d0 - ( mu + yn ) / mu
    zndev = 2.0d0 - ( mu + zn ) / mu
    epslon = max ( abs ( xndev ), abs ( yndev ), abs ( zndev ) )

    if ( epslon < errtol ) then
      c1 = 1.0d0 / 24.0d0
      c2 = 3.0d0 / 44.0d0
      c3 = 1.0d0 / 14.0d0
      e2 = xndev * yndev - zndev * zndev
      e3 = xndev * yndev * zndev
      s = 1.0d0 + ( c1 * e2 - 0.1d0 - c2 * e3 ) * e2 + c3 * e3
      rf = s / sqrt ( mu )
      return
    end if

    xnroot = sqrt ( xn )
    ynroot = sqrt ( yn )
    znroot = sqrt ( zn )
    lamda = xnroot * ( ynroot + znroot ) + ynroot * znroot
    xn = ( xn + lamda ) * 0.25d0
    yn = ( yn + lamda ) * 0.25d0
    zn = ( zn + lamda ) * 0.25d0

  end do

end
function rj ( x, y, z, p, errtol, ierr )

!*****************************************************************************80
!
!! RJ computes an incomplete elliptic integral of the third kind, RJ(X,Y,Z,P).
!
!  Discussion:
!
!    This function computes an incomplete elliptic integral of the third kind.
!
!    RJ(X,Y,Z,P) = Integral ( 0 <= T < oo )
!
!                  -1/2     -1/2     -1/2     -1
!        (3/2)(T+X)    (T+Y)    (T+Z)    (T+P)  DT,
!
!    where X, Y, and Z are nonnegative, at most one of them is
!    zero, and P is positive.
!
!    If X or Y or Z is zero, then the integral is complete.
!
!    The duplication theorem is iterated until the variables are nearly equal,
!    and the function is then expanded in Taylor series to fifth order.
!
!    Check by addition theorem:
!
!      RJ(X,X+Z,X+W,X+P)
!      + RJ(Y,Y+Z,Y+W,Y+P) + (A-B) * RJ(A,B,B,A) + 3 / sqrt ( A)
!      = RJ(0,Z,W,P), where X,Y,Z,W,P are positive and X * Y
!      = Z * W,  A = P * P * (X+Y+Z+W),  B = P * (P+X) * (P+Y),
!      and B - A = P * (P-Z) * (P-W).
!
!    The sum of the third and fourth terms on the left side is 3 * RC(A,B).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 May 2018
!
!  Author:
!
!    Original FORTRAN77 version by Bille Carlson, Elaine Notis.
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Bille Carlson,
!    Computing Elliptic Integrals by Duplication,
!    Numerische Mathematik,
!    Volume 33, 1979, pages 1-16.
!
!    Bille Carlson, Elaine Notis,
!    Algorithm 577, Algorithms for Incomplete Elliptic Integrals,
!    ACM Transactions on Mathematical Software,
!    Volume 7, Number 3, pages 398-403, September 1981.
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, Y, Z, P, the arguments in the integral.
!
!    Input, real ( kind = rk ) ERRTOL, the error tolerance.
!    Relative error due to truncation of the series for rj
!    is less than 3 * ERRTOL ^ 6 / (1 - ERRTOL) ^ 3/2.
!    An error tolerance (ETOLRC) will be passed to the subroutine
!    for RC to make the truncation error for RC less than for RJ.
!    Sample choices:
!      ERRTOL   Relative truncation error less than
!      1.D-3    4.D-18
!      3.D-3    3.D-15
!      1.D-2    4.D-12
!      3.D-2    3.D-9
!      1.D-1    4.D-6
!
!    Output, integer IERR, the error flag.
!    0, no error occurred.
!    1, abnormal termination.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) alfa
  real ( kind = rk ) beta
  real ( kind = rk ) c1
  real ( kind = rk ) c2
  real ( kind = rk ) c3
  real ( kind = rk ) c4
  real ( kind = rk ) ea
  real ( kind = rk ) eb
  real ( kind = rk ) ec
  real ( kind = rk ) e2
  real ( kind = rk ) e3
  real ( kind = rk ) epslon
  real ( kind = rk ) errtol
  real ( kind = rk ) etolrc
  integer ierr
  real ( kind = rk ) lamda
  real ( kind = rk ) lolim
  real ( kind = rk ) mu
  real ( kind = rk ) p
  real ( kind = rk ) pn
  real ( kind = rk ) pndev
  real ( kind = rk ) power4
  real ( kind = rk ) rc
  real ( kind = rk ) rj
  real ( kind = rk ) sigma
  real ( kind = rk ) s1
  real ( kind = rk ) s2
  real ( kind = rk ) s3
  real ( kind = rk ) uplim
  real ( kind = rk ) x
  real ( kind = rk ) xn
  real ( kind = rk ) xndev
  real ( kind = rk ) xnroot
  real ( kind = rk ) y
  real ( kind = rk ) yn
  real ( kind = rk ) yndev
  real ( kind = rk ) ynroot
  real ( kind = rk ) z
  real ( kind = rk ) zn
  real ( kind = rk ) zndev
  real ( kind = rk ) znroot
!
!  LOLIM AND UPLIM DETERMINE THE RANGE OF VALID ARGUMENTS.
!  LOLIM IS NOT LESS THAN THE CUBE ROOT OF THE VALUE
!  OF LOLIM USED IN THE SUBROUTINE FOR RC.
!  UPLIM IS NOT GREATER THAN 0.3 TIMES THE CUBE ROOT OF
!  THE VALUE OF UPLIM USED IN THE SUBROUTINE FOR RC.
!
  save lolim
  save uplim

  data lolim /2.D-26/
  data uplim /3.D+24/

  if ( &
    x < 0.0D+00 .or. &
    y < 0.0D+00 .or. &
    z < 0.0D+00 .or. &
    x + y < lolim .or. &
    x + z < lolim .or. &
    y + z < lolim .or. &
    p < lolim .or. &
    uplim < x .or. &
    uplim < y .or. &
    uplim < z .or. &
    uplim < p ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'RJ - Error!'
    write ( *, '(a)' ) '  Invalid input arguments.'
    write ( *, '(a,d23.16)' ) '  X = ', x
    write ( *, '(a,d23.16)' ) '  Y = ', y
    write ( *, '(a,d23.16)' ) '  Z = ', z
    write ( *, '(a,d23.16)' ) '  P = ', p
    write ( *, '(a)' ) ''
    ierr = 1
    rj = 0.0D+00
    return
  end if

  ierr = 0
  xn = x
  yn = y
  zn = z
  pn = p
  sigma = 0.0d0
  power4 = 1.0d0
  etolrc = 0.5d0 * errtol

  do

    mu = ( xn + yn + zn + pn + pn ) * 0.2d0
    xndev = ( mu - xn ) / mu
    yndev = ( mu - yn ) / mu
    zndev = ( mu - zn ) / mu
    pndev = ( mu - pn ) / mu
    epslon = max ( abs ( xndev ), abs ( yndev ), abs ( zndev ), abs ( pndev ) )

    if ( epslon < errtol ) then
      c1 = 3.0d0 / 14.0d0
      c2 = 1.0d0 / 3.0d0
      c3 = 3.0d0 / 22.0d0
      c4 = 3.0d0 / 26.0d0
      ea = xndev * ( yndev + zndev ) + yndev * zndev
      eb = xndev * yndev * zndev
      ec = pndev * pndev
      e2 = ea - 3.0d0 * ec
      e3 = eb + 2.0d0 * pndev * ( ea - ec )
      s1 = 1.0d0 + e2 * ( - c1 + 0.75d0 * c3 * e2 - 1.5d0 * c4 * e3 )
      s2 = eb * ( 0.5d0 * c2 + pndev * ( - c3 - c3 + pndev * c4 ) )
      s3 = pndev * ea * ( c2 - pndev * c3 ) - c2 * pndev * ec
      rj = 3.0d0 * sigma + power4 * ( s1 + s2 + s3 ) / ( mu * sqrt ( mu ) )
      return
    end if

    xnroot = sqrt ( xn )
    ynroot = sqrt ( yn )
    znroot = sqrt ( zn )
    lamda = xnroot * ( ynroot + znroot ) + ynroot * znroot
    alfa = pn * ( xnroot + ynroot + znroot ) &
      + xnroot * ynroot * znroot
    alfa = alfa * alfa
    beta = pn * ( pn + lamda ) * ( pn + lamda )
    sigma = sigma + power4 * rc ( alfa, beta, etolrc, ierr )

    if ( ierr /= 0 ) then
      rj = 0.0D+00
      return
    end if

    power4 = power4 * 0.25d0
    xn = ( xn + lamda ) * 0.25d0
    yn = ( yn + lamda ) * 0.25d0
    zn = ( zn + lamda ) * 0.25d0
    pn = ( pn + lamda ) * 0.25d0

  end do

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

