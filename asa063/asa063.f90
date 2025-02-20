function alogam ( x, ifault )

!*****************************************************************************80
!
!! alogam() computes the logarithm of the Gamma function.
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
!    Original FORTRAN77 version by Malcolm Pike, David Hill.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Malcolm Pike, David Hill,
!    Algorithm 291:
!    Logarithm of Gamma Function,
!    Communications of the ACM,
!    Volume 9, Number 9, September 1966, page 684.
!
!  Input:
!
!    real ( kind = rk ) X, the argument of the Gamma function.
!    X should be greater than 0.
!
!  Output:
!
!    integer IFAULT, error flag.
!    0, no error.
!    1, X <= 0.
!
!    real ( kind = rk ) ALOGAM, the logarithm of the Gamma
!    function of X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) alogam
  real ( kind = rk ) f
  integer ifault
  real ( kind = rk ) x
  real ( kind = rk ) y
  real ( kind = rk ) z

  if ( x <= 0.0D+00 ) then
    ifault = 1
    alogam = 0.0D+00
    return
  end if

  ifault = 0
  y = x

  if ( x < 7.0D+00 ) then

    f = 1.0D+00
    z = y

    do while ( z < 7.0D+00 )
      f = f * z
      z = z + 1.0D+00
    end do

    y = z
    f = - log ( f )

  else

    f = 0.0D+00

  end if

  z = 1.0D+00 / y / y

  alogam = f + ( y - 0.5D+00 ) * log ( y ) - y &
    + 0.918938533204673D+00 + &
    ((( &
    - 0.000595238095238D+00   * z &
    + 0.000793650793651D+00 ) * z &
    - 0.002777777777778D+00 ) * z &
    + 0.083333333333333D+00 ) / y

  return
end
subroutine beta_inc_values ( n_data, a, b, x, fx )

!*****************************************************************************80
!
!! beta_inc_values() returns some values of the incomplete Beta function.
!
!  Discussion:
!
!    The incomplete Beta function may be written
!
!      BETA_INC(A,B,X) = integral (0 <= t <= X) T^(A-1) * (1-T)^(B-1) dT
!                      / integral (0 <= t <= 1) T^(A-1) * (1-T)^(B-1) dT
!
!    Thus,
!
!      BETA_INC(A,B,0.0) = 0.0
!      BETA_INC(A,B,1.0) = 1.0
!
!    The incomplete Beta function is also sometimes called the
!    "modified" Beta function, or the "normalized" Beta function
!    or the Beta CDF (cumulative density function).
!
!    In Mathematica, the function can be evaluated by:
!
!      BETA[X,A,B] / BETA[A,B]
!
!    The function can also be evaluated by using the Statistics package:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = BetaDistribution [ a, b ]
!      CDF [ dist, x ]
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
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Karl Pearson,
!    Tables of the Incomplete Beta Function,
!    Cambridge University Press, 1968,
!    LC: QA351.P38.
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
!    integer N_DATA.  The user sets N_DATA to 0 before the first call.  
!
!  Output:
!
!    integer N_DATA.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    real ( kind = rk ) A, B, the parameters of the function.
!
!    real ( kind = rk ) X, the argument of the function.
!
!    real ( kind = rk ) FX, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 45

  real ( kind = rk ) a
  real ( kind = rk ), save, dimension ( n_max ) :: a_vec = (/ &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     1.0D+00, &
     1.0D+00, &
     1.0D+00, &
     1.0D+00, &
     1.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     5.5D+00, &
    10.0D+00, &
    10.0D+00, &
    10.0D+00, &
    10.0D+00, &
    20.0D+00, &
    20.0D+00, &
    20.0D+00, &
    20.0D+00, &
    20.0D+00, &
    30.0D+00, &
    30.0D+00, &
    40.0D+00, &
     1.0D+00, &
     1.0D+00, &
     1.0D+00, &
     1.0D+00, &
     1.0D+00, &
     1.0D+00, &
     1.0D+00, &
     1.0D+00, &
     2.0D+00, &
     3.0D+00, &
     4.0D+00, &
     5.0D+00, &
     1.30625D+00, &
     1.30625D+00, &
     1.30625D+00 /)
  real ( kind = rk ) b
  real ( kind = rk ), save, dimension ( n_max ) :: b_vec = (/ &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     1.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     5.0D+00, &
     0.5D+00, &
     5.0D+00, &
     5.0D+00, &
    10.0D+00, &
     5.0D+00, &
    10.0D+00, &
    10.0D+00, &
    20.0D+00, &
    20.0D+00, &
    10.0D+00, &
    10.0D+00, &
    20.0D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     0.5D+00, &
     2.0D+00, &
     3.0D+00, &
     4.0D+00, &
     5.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
     2.0D+00, &
    11.7562D+00, &
    11.7562D+00, & 
    11.7562D+00 /)
  real ( kind = rk ) fx
  real ( kind = rk ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.6376856085851985D-01, &
    0.2048327646991335D+00, &
    0.1000000000000000D+01, &
    0.0000000000000000D+00, &
    0.5012562893380045D-02, &
    0.5131670194948620D-01, &
    0.2928932188134525D+00, &
    0.5000000000000000D+00, &
    0.2800000000000000D-01, &
    0.1040000000000000D+00, &
    0.2160000000000000D+00, &
    0.3520000000000000D+00, &
    0.5000000000000000D+00, &
    0.6480000000000000D+00, &
    0.7840000000000000D+00, &
    0.8960000000000000D+00, &
    0.9720000000000000D+00, &
    0.4361908850559777D+00, &
    0.1516409096347099D+00, &
    0.8978271484375000D-01, &
    0.1000000000000000D+01, &
    0.5000000000000000D+00, &
    0.4598773297575791D+00, &
    0.2146816102371739D+00, &
    0.9507364826957875D+00, &
    0.5000000000000000D+00, &
    0.8979413687105918D+00, &
    0.2241297491808366D+00, &
    0.7586405487192086D+00, &
    0.7001783247477069D+00, &
    0.5131670194948620D-01, &
    0.1055728090000841D+00, &
    0.1633399734659245D+00, &
    0.2254033307585166D+00, &
    0.3600000000000000D+00, &
    0.4880000000000000D+00, &
    0.5904000000000000D+00, &
    0.6723200000000000D+00, &
    0.2160000000000000D+00, &
    0.8370000000000000D-01, &
    0.3078000000000000D-01, &
    0.1093500000000000D-01, &
    0.918884684620518D+00, &
    0.21052977489419D+00, &
    0.1824130512500673D+00 /)
  integer n_data
  real ( kind = rk ) x
  real ( kind = rk ), save, dimension ( n_max ) :: x_vec = (/ &
    0.01D+00, &
    0.10D+00, &
    1.00D+00, &
    0.00D+00, &
    0.01D+00, &
    0.10D+00, &
    0.50D+00, &
    0.50D+00, &
    0.10D+00, &
    0.20D+00, &
    0.30D+00, &
    0.40D+00, &
    0.50D+00, &
    0.60D+00, &
    0.70D+00, &
    0.80D+00, &
    0.90D+00, &
    0.50D+00, &
    0.90D+00, &
    0.50D+00, &
    1.00D+00, &
    0.50D+00, &
    0.80D+00, &
    0.60D+00, &
    0.80D+00, &
    0.50D+00, &
    0.60D+00, &
    0.70D+00, &
    0.80D+00, &
    0.70D+00, &
    0.10D+00, &
    0.20D+00, &
    0.30D+00, &
    0.40D+00, &
    0.20D+00, &
    0.20D+00, &
    0.20D+00, &
    0.20D+00, &
    0.30D+00, &
    0.30D+00, &
    0.30D+00, &
    0.30D+00, &
    0.225609D+00, &
    0.0335568D+00, &
    0.0295222D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0D+00
    b = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    b = b_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function betain ( x, p, q, beta, ifault )

!*****************************************************************************80
!
!! betain() computes the incomplete Beta function ratio.
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
!    Original FORTRAN77 version by KL Majumder, GP Bhattacharjee.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    KL Majumder, GP Bhattacharjee,
!    Algorithm AS 63:
!    The incomplete Beta Integral,
!    Applied Statistics,
!    Volume 22, Number 3, 1973, pages 409-411.
!
!  Input:
!
!    real ( kind = rk ) X, the argument, between 0 and 1.
!
!    real ( kind = rk ) P, Q, the parameters, which
!    must be positive.
!
!    real ( kind = rk ) BETA, the logarithm of the complete
!    beta function.
!
!  Output:
!
!    integer IFAULT, error flag.
!    0, no error.
!    nonzero, an error occurred.
!
!    real ( kind = rk ) BETAIN, the value of the incomplete
!    Beta function ratio.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), parameter :: acu = 0.1D-14
  real ( kind = rk ) ai
  real ( kind = rk ) beta
  real ( kind = rk ) betain
  real ( kind = rk ) cx
  integer ifault
  logical indx
  integer ns
  real ( kind = rk ) p
  real ( kind = rk ) pp
  real ( kind = rk ) psq
  real ( kind = rk ) q
  real ( kind = rk ) qq
  real ( kind = rk ) rx
  real ( kind = rk ) temp
  real ( kind = rk ) term
  real ( kind = rk ) x
  real ( kind = rk ) xx

  betain = x
  ifault = 0
!
!  Check the input arguments.
!
  if ( p <= 0.0D+00 .or. q <= 0.0D+00 ) then
    ifault = 1
    return
  end if

  if ( x < 0.0D+00 .or. 1.0D+00 < x ) then
    ifault = 2
    return
  end if
!
!  Special cases.
!
  if ( x == 0.0D+00 .or. x == 1.0D+00 ) then
    return
  end if
!
!  Change tail if necessary and determine S.
!
  psq = p + q
  cx = 1.0D+00 - x

  if ( p < psq * x ) then
    xx = cx
    cx = x
    pp = q
    qq = p
    indx = .true.
  else
    xx = x
    pp = p
    qq = q
    indx = .false.
  end if

  term = 1.0D+00
  ai = 1.0D+00
  betain = 1.0D+00
  ns = int ( qq + cx * psq )
!
!  Use Soper's reduction formula.
!
  rx = xx / cx
  temp = qq - ai
  if ( ns == 0 ) then
    rx = xx
  end if

  do

    term = term * temp * rx / ( pp + ai )
    betain = betain + term
    temp = abs ( term )

    if ( temp <= acu .and. temp <= acu * betain ) then

      betain = betain * exp ( pp * log ( xx ) &
        + ( qq - 1.0D+00 ) * log ( cx ) - beta ) / pp

      if ( indx ) then
        betain = 1.0D+00 - betain
      end if

      exit

    end if

    ai = ai + 1.0D+00
    ns = ns - 1

    if ( 0 <= ns ) then
      temp = qq - ai
      if ( ns == 0 ) then
        rx = xx
      end if
    else
      temp = psq
      psq = psq + 1.0D+00
    end if

  end do

  return
end

