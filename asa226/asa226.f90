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
!    28 March 1999
!
!  Author:
!
!    Original FORTRAN77 version by Malcolm Pike, David Hill.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Malcolm Pike, David Hill,
!    Algorithm 291:
!    Logarithm of Gamma Function,
!    Communications of the ACM,
!    Volume 9, Number 9, September 1966, page 684.
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, the argument of the Gamma function.
!    X should be greater than 0.
!
!    Output, integer IFAULT, error flag.
!    0, no error.
!    1, X <= 0.
!
!    Output, real ( kind = rk ) ALOGAM, the logarithm of the Gamma
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
subroutine beta_noncentral_cdf_values ( n_data, a, b, lambda, x, fx )

!*****************************************************************************80
!
!! beta_noncentral_cdf_values() returns some values of the noncentral Beta CDF.
!
!  Discussion:
!
!    The values presented here are taken from the reference, where they
!    were given to a limited number of decimal places.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    R Chattamvelli, R Shanmugam,
!    Algorithm AS 310:
!    Computing the Non-central Beta Distribution Function,
!    Applied Statistics,
!    Volume 46, Number 1, 1997, pages 146-156.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = rk ) A, B, the shape parameters.
!
!    Output, real ( kind = rk ) LAMBDA, the noncentrality parameter.
!
!    Output, real ( kind = rk ) X, the argument of the function.
!
!    Output, real ( kind = rk ) FX, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 25

  real ( kind = rk ) a
  real ( kind = rk ), save, dimension ( n_max ) :: a_vec = (/ &
        5.0D+00, &
        5.0D+00, &
        5.0D+00, &
       10.0D+00, &
       10.0D+00, &
       10.0D+00, &
       20.0D+00, &
       20.0D+00, &
       20.0D+00, &
       10.0D+00, &
       10.0D+00, &
       15.0D+00, &
       20.0D+00, &
       20.0D+00, &
       20.0D+00, &
       30.0D+00, &
       30.0D+00, &
       10.0D+00, &
       10.0D+00, &
       10.0D+00, &
       15.0D+00, &
       10.0D+00, &
       12.0D+00, &
       30.0D+00, &
       35.0D+00 /)
  real ( kind = rk ) b
  real ( kind = rk ), save, dimension ( n_max ) :: b_vec = (/ &
        5.0D+00, &
        5.0D+00, &
        5.0D+00, &
       10.0D+00, &
       10.0D+00, &
       10.0D+00, &
       20.0D+00, &
       20.0D+00, &
       20.0D+00, &
       20.0D+00, &
       10.0D+00, &
        5.0D+00, &
       10.0D+00, &
       30.0D+00, &
       50.0D+00, &
       20.0D+00, &
       40.0D+00, &
        5.0D+00, &
       10.0D+00, &
       30.0D+00, &
       20.0D+00, &
        5.0D+00, &
       17.0D+00, &
       30.0D+00, &
       30.0D+00 /)
  real ( kind = rk ) fx
  real ( kind = rk ), save, dimension ( n_max ) :: fx_vec = (/ &
       0.4563021D+00, &
       0.1041337D+00, &
       0.6022353D+00, &
       0.9187770D+00, &
       0.6008106D+00, &
       0.0902850D+00, &
       0.9998655D+00, &
       0.9925997D+00, &
       0.9641112D+00, &
       0.9376626573D+00, &
       0.7306817858D+00, &
       0.1604256918D+00, &
       0.1867485313D+00, &
       0.6559386874D+00, &
       0.9796881486D+00, &
       0.1162386423D+00, &
       0.9930430054D+00, &
       0.0506899273D+00, &
       0.1030959706D+00, &
       0.9978417832D+00, &
       0.2555552369D+00, &
       0.0668307064D+00, &
       0.0113601067D+00, &
       0.7813366615D+00, &
       0.8867126477D+00 /)
  real ( kind = rk ) lambda
  real ( kind = rk ), save, dimension ( n_max ) :: lambda_vec = (/ &
        54.0D+00, &
       140.0D+00, &
       170.0D+00, &
        54.0D+00, &
       140.0D+00, &
       250.0D+00, &
        54.0D+00, &
       140.0D+00, &
       250.0D+00, &
       150.0D+00, &
       120.0D+00, &
        80.0D+00, &
       110.0D+00, &
        65.0D+00, &
       130.0D+00, &
        80.0D+00, &
       130.0D+00, &
        20.0D+00, &
        54.0D+00, &
        80.0D+00, &
       120.0D+00, &
        55.0D+00, &
        64.0D+00, &
       140.0D+00, &
        20.0D+00 /)
  integer n_data
  real ( kind = rk ) x
  real ( kind = rk ), save, dimension ( n_max ) :: x_vec = (/ &
       0.8640D+00, &
       0.9000D+00, &
       0.9560D+00, &
       0.8686D+00, &
       0.9000D+00, &
       0.9000D+00, &
       0.8787D+00, &
       0.9000D+00, &
       0.9220D+00, &
       0.868D+00, &
       0.900D+00, &
       0.880D+00, &
       0.850D+00, &
       0.660D+00, &
       0.720D+00, &
       0.720D+00, &
       0.800D+00, &
       0.644D+00, &
       0.700D+00, &
       0.780D+00, &
       0.760D+00, &
       0.795D+00, &
       0.560D+00, &
       0.800D+00, &
       0.670D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0D+00
    b = 0.0D+00
    lambda = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    b = b_vec(n_data)
    lambda = lambda_vec(n_data)
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
!    12 January 2008
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
!  Parameters:
!
!    Input, real ( kind = rk ) X, the argument, between 0 and 1.
!
!    Input, real ( kind = rk ) P, Q, the parameters, which
!    must be positive.
!
!    Input, real ( kind = rk ) BETA, the logarithm of the complete
!    beta function.
!
!    Output, integer IFAULT, error flag.
!    0, no error.
!    nonzero, an error occurred.
!
!    Output, real ( kind = rk ) BETAIN, the value of the incomplete
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
function betanc ( x, a, b, lambda, ifault )

!*****************************************************************************80
!
!! betanc() computes the tail of the noncentral Beta distribution.
!
!  Discussion:
!
!    This routine returns the cumulative probability of X for the non-central
!    Beta distribution with parameters A, B and non-centrality LAMBDA.
!
!    Note that if LAMBDA = 0, the standard Beta distribution is defined.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by Russell Lenth.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Russell Lenth,
!    Algorithm AS 226:
!    Computing Noncentral Beta Probabilities,
!    Applied Statistics,
!    Volume 36, Number 2, 1987, pages 241-244.
!
!    H Frick,
!    Algorithm AS R84:
!    A Remark on Algorithm AS 226:
!    Computing Noncentral Beta Probabilities,
!    Applied Statistics,
!    Volume 39, Number 2, 1990, pages 311-312.
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, the value defining the cumulative
!    probability lower tail.  Normally, 0 <= X <= 1, but any value
!    is allowed.
!
!    Input, real ( kind = rk ) A, B, the parameters of the distribution.
!    0 < A, 0 < B.
!
!    Input, real ( kind = rk ) LAMBDA, the noncentrality parameter
!    of the distribution.  0 <= LAMBDA.  The program can produce reasonably
!    accurate results for values of LAMBDA up to about 100.
!
!    Output, integer IFAULT, error flag.
!    0, no error occurred.
!    nonzero, an error occurred.
!
!    Output, real ( kind = rk ) BETANC, the cumulative probability
!    of X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) alogam
  real ( kind = rk ) ax
  real ( kind = rk ) b
  real ( kind = rk ) beta
  real ( kind = rk ) betain
  real ( kind = rk ) betanc
  real ( kind = rk ) c
  real ( kind = rk ) errbd
  real ( kind = rk ), parameter :: errmax = 1.0D-07
  real ( kind = rk ) gx
  integer ifault
  integer, parameter :: itrmax = 150
  real ( kind = rk ) lambda
  real ( kind = rk ) q
  real ( kind = rk ) sumq
  real ( kind = rk ) temp
  real ( kind = rk ), parameter :: ualpha = 5.0D+00
  real ( kind = rk ) x
  real ( kind = rk ) xj

  ifault = 0

  if ( lambda < 0.0D+00 .or. &
       a <= 0.0D+00 .or. &
       b <= 0.0D+00 ) then
    ifault = 2
    betanc = -1.0D+00
    return
  end if

  if ( x <= 0.0D+00 ) then
    betanc = 0.0D+00
    return
  end if

  if ( 1.0D+00 <= x ) then
    betanc = 1.0D+00
    return
  end if

  c = 0.5D+00 * lambda
!
!  Initialize the series.
!
  beta = alogam ( a, ifault ) &
       + alogam ( b, ifault ) &
       - alogam ( a + b, ifault )

  temp = betain ( x, a, b, beta, ifault )

  gx = exp ( a * log ( x ) + b * log ( 1.0D+00 - x ) &
    - beta - log ( a ) )

  q = exp ( - c )

  xj = 0.0D+00
  ax = q * temp
  sumq = 1.0D+00 - q
  betanc = ax
!
!  Recur over subsequent terms until convergence is achieved.
!
  ifault = 1

  do

    xj = xj + 1.0D+00
    temp = temp - gx
    gx = x * ( a + b + xj - 1.0D+00 ) * gx / ( a + xj )
    q = q * c / xj
    sumq = sumq - q
    ax = temp * q
    betanc = betanc + ax
!
!  Check for convergence and act accordingly.
!
    errbd = abs ( ( temp - gx ) * sumq )

    if ( errbd <= errmax ) then
      ifault = 0
      exit
    end if

    if (  itrmax < int ( xj ) ) then
      exit
    end if

  end do

  return
end

