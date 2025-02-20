function alngam ( xvalue, ifault )

!*****************************************************************************80
!
!! alngam() computes the logarithm of the gamma function.
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
!    Original FORTRAN77 version by Allan Macleod.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Allan Macleod,
!    Algorithm AS 245,
!    A Robust and Reliable Algorithm for the Logarithm of the Gamma Function,
!    Applied Statistics,
!    Volume 38, Number 2, 1989, pages 397-402.
!
!  Parameters:
!
!    Input, real ( kind = rk ) XVALUE, the argument of the Gamma function.
!
!    Output, integer IFAULT, error flag.
!    0, no error occurred.
!    1, XVALUE is less than or equal to 0.
!    2, XVALUE is too big.
!
!    Output, real ( kind = rk ) ALNGAM, the logarithm of the gamma function of X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) alngam
  real ( kind = rk ), parameter :: alr2pi = 0.918938533204673D+00
  integer ifault
  real ( kind = rk ), dimension ( 9 ) :: r1 = (/ &
    -2.66685511495D+00, &
    -24.4387534237D+00, &
    -21.9698958928D+00, &
     11.1667541262D+00, &
     3.13060547623D+00, &
     0.607771387771D+00, &
     11.9400905721D+00, &
     31.4690115749D+00, &
     15.2346874070D+00 /)
  real ( kind = rk ), dimension ( 9 ) :: r2 = (/ &
    -78.3359299449D+00, &
    -142.046296688D+00, &
     137.519416416D+00, &
     78.6994924154D+00, &
     4.16438922228D+00, &
     47.0668766060D+00, &
     313.399215894D+00, &
     263.505074721D+00, &
     43.3400022514D+00 /)
  real ( kind = rk ), dimension ( 9 ) :: r3 = (/ &
    -2.12159572323D+05, &
     2.30661510616D+05, &
     2.74647644705D+04, &
    -4.02621119975D+04, &
    -2.29660729780D+03, &
    -1.16328495004D+05, &
    -1.46025937511D+05, &
    -2.42357409629D+04, &
    -5.70691009324D+02 /)
  real ( kind = rk ), dimension ( 5 ) :: r4 = (/ &
     0.279195317918525D+00, &
     0.4917317610505968D+00, &
     0.0692910599291889D+00, &
     3.350343815022304D+00, &
     6.012459259764103D+00 /)
  real ( kind = rk ) x
  real ( kind = rk ) x1
  real ( kind = rk ) x2
  real ( kind = rk ), parameter :: xlge = 5.10D+05
  real ( kind = rk ), parameter :: xlgst = 1.0D+30
  real ( kind = rk ) xvalue
  real ( kind = rk ) y

  x = xvalue
  alngam = 0.0D+00
!
!  Check the input.
!
  if ( xlgst <= x ) then
    ifault = 2
    return
  end if

  if ( x <= 0.0D+00 ) then
    ifault = 1
    return
  end if

  ifault = 0
!
!  Calculation for 0 < X < 0.5 and 0.5 <= X < 1.5 combined.
!
  if ( x < 1.5D+00 ) then

    if ( x < 0.5D+00 ) then

      alngam = - log ( x )
      y = x + 1.0D+00
!
!  Test whether X < machine epsilon.
!
      if ( y == 1.0D+00 ) then
        return
      end if

    else

      alngam = 0.0D+00
      y = x
      x = ( x - 0.5D+00 ) - 0.5D+00

    end if

    alngam = alngam + x * (((( &
        r1(5)   * y &
      + r1(4) ) * y &
      + r1(3) ) * y &
      + r1(2) ) * y &
      + r1(1) ) / (((( &
                  y &
      + r1(9) ) * y &
      + r1(8) ) * y &
      + r1(7) ) * y &
      + r1(6) )

    return

  end if
!
!  Calculation for 1.5 <= X < 4.0.
!
  if ( x < 4.0D+00 ) then

    y = ( x - 1.0D+00 ) - 1.0D+00

    alngam = y * (((( &
        r2(5)   * x &
      + r2(4) ) * x &
      + r2(3) ) * x &
      + r2(2) ) * x &
      + r2(1) ) / (((( &
                  x &
      + r2(9) ) * x &
      + r2(8) ) * x &
      + r2(7) ) * x &
      + r2(6) )
!
!  Calculation for 4.0 <= X < 12.0.
!
  else if ( x < 12.0D+00 ) then

    alngam = (((( &
        r3(5)   * x &
      + r3(4) ) * x &
      + r3(3) ) * x &
      + r3(2) ) * x &
      + r3(1) ) / (((( &
                  x &
      + r3(9) ) * x &
      + r3(8) ) * x &
      + r3(7) ) * x &
      + r3(6) )
!
!  Calculation for 12.0 <= X.
!
  else

    y = log ( x )
    alngam = x * ( y - 1.0D+00 ) - 0.5D+00 * y + alr2pi

    if ( x <= xlge ) then

      x1 = 1.0D+00 / x
      x2 = x1 * x1

      alngam = alngam + x1 * ( ( &
             r4(3)   * &
        x2 + r4(2) ) * &
        x2 + r4(1) ) / ( ( &
        x2 + r4(5) ) * &
        x2 + r4(4) )

    end if

  end if

  return
end
function alnorm ( x, upper )

!*****************************************************************************80
!
!! alnorm() computes the cumulative density of the standard normal distribution.
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
!    Original FORTRAN77 version by David Hill.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Hill,
!    Algorithm AS 66:
!    The Normal Integral,
!    Applied Statistics,
!    Volume 22, Number 3, 1973, pages 424-427.
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, is one endpoint of the semi-infinite interval
!    over which the integration takes place.
!
!    Input, logical UPPER, determines whether the upper or lower
!    interval is to be integrated:
!    .TRUE.  => integrate from X to + Infinity;
!    .FALSE. => integrate from - Infinity to X.
!
!    Output, real ( kind = rk ) ALNORM, the integral of the standard normal
!    distribution over the desired interval.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), parameter :: a1 = 5.75885480458D+00
  real ( kind = rk ), parameter :: a2 = 2.62433121679D+00
  real ( kind = rk ), parameter :: a3 = 5.92885724438D+00
  real ( kind = rk ) alnorm
  real ( kind = rk ), parameter :: b1 = -29.8213557807D+00
  real ( kind = rk ), parameter :: b2 = 48.6959930692D+00
  real ( kind = rk ), parameter :: c1 = -0.000000038052D+00
  real ( kind = rk ), parameter :: c2 = 0.000398064794D+00
  real ( kind = rk ), parameter :: c3 = -0.151679116635D+00
  real ( kind = rk ), parameter :: c4 = 4.8385912808D+00
  real ( kind = rk ), parameter :: c5 = 0.742380924027D+00
  real ( kind = rk ), parameter :: c6 = 3.99019417011D+00
  real ( kind = rk ), parameter :: con = 1.28D+00
  real ( kind = rk ), parameter :: d1 = 1.00000615302D+00
  real ( kind = rk ), parameter :: d2 = 1.98615381364D+00
  real ( kind = rk ), parameter :: d3 = 5.29330324926D+00
  real ( kind = rk ), parameter :: d4 = -15.1508972451D+00
  real ( kind = rk ), parameter :: d5 = 30.789933034D+00
  real ( kind = rk ), parameter :: ltone = 7.0D+00
  real ( kind = rk ), parameter :: p = 0.398942280444D+00
  real ( kind = rk ), parameter :: q = 0.39990348504D+00
  real ( kind = rk ), parameter :: r = 0.398942280385D+00
  logical up
  logical upper
  real ( kind = rk ), parameter :: utzero = 18.66D+00
  real ( kind = rk ) x
  real ( kind = rk ) y
  real ( kind = rk ) z

  up = upper
  z = x

  if ( z < 0.0D+00 ) then
    up = .not. up
    z = - z
  end if

  if ( ltone < z .and. ( ( .not. up ) .or. utzero < z ) ) then

    if ( up ) then
      alnorm = 0.0D+00
    else
      alnorm = 1.0D+00
    end if

    return

  end if

  y = 0.5D+00 * z * z

  if ( z <= con ) then

    alnorm = 0.5D+00 - z * ( p - q * y &
      / ( y + a1 + b1 &
      / ( y + a2 + b2 &
      / ( y + a3 ))))

  else

    alnorm = r * exp ( - y ) &
      / ( z + c1 + d1 &
      / ( z + c2 + d2 &
      / ( z + c3 + d3 &
      / ( z + c4 + d4 &
      / ( z + c5 + d5 &
      / ( z + c6 ))))))

  end if

  if ( .not. up ) then
    alnorm = 1.0D+00 - alnorm
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
subroutine student_noncentral_cdf_values ( n_data, df, lambda, x, fx )

!*****************************************************************************80
!
!! student_noncentral_cdf_values() returns values of the noncentral Student CDF.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = NoncentralStudentTDistribution [ df, lambda ]
!      CDF [ dist, x ]
!
!    Mathematica seems to have some difficulty computing this function
!    to the desired number of digits.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 September 2004
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
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer DF, real ( kind = rk ) LAMBDA, the parameters
!    of the function.
!
!    Output, real ( kind = rk ) X, the argument of the function.
!
!    Output, real ( kind = rk ) FX, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 30

  integer df
  integer, save, dimension ( n_max ) :: df_vec = (/ &
     1,  2,  3, &
     1,  2,  3, &
     1,  2,  3, &
     1,  2,  3, &
     1,  2,  3, &
    15, 20, 25, &
     1,  2,  3, &
    10, 10, 10, &
    10, 10, 10, &
    10, 10, 10 /)
  real ( kind = rk ) fx
  real ( kind = rk ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.8975836176504333D+00, &
    0.9522670169D+00, &
    0.9711655571887813D+00, &
    0.8231218864D+00, &
    0.9049021510D+00, &
    0.9363471834D+00, &
    0.7301025986D+00, &
    0.8335594263D+00, &
    0.8774010255D+00, &
    0.5248571617D+00, &
    0.6293856597D+00, &
    0.6800271741D+00, &
    0.20590131975D+00, &
    0.2112148916D+00, &
    0.2074730718D+00, &
    0.9981130072D+00, &
    0.9994873850D+00, &
    0.9998391562D+00, &
    0.168610566972D+00, &
    0.16967950985D+00, &
    0.1701041003D+00, &
    0.9247683363D+00, &
    0.7483139269D+00, &
    0.4659802096D+00, &
    0.9761872541D+00, &
    0.8979689357D+00, &
    0.7181904627D+00, &
    0.9923658945D+00, &
    0.9610341649D+00, &
    0.8688007350D+00 /)
  real ( kind = rk ) lambda
  real ( kind = rk ), save, dimension ( n_max ) :: lambda_vec = (/ &
    0.0D+00, &
    0.0D+00, &
    0.0D+00, &
    0.5D+00, &
    0.5D+00, &
    0.5D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    2.0D+00, &
    2.0D+00, &
    2.0D+00, &
    4.0D+00, &
    4.0D+00, &
    4.0D+00, &
    7.0D+00, &
    7.0D+00, &
    7.0D+00, &
    1.0D+00, &
    1.0D+00, &
    1.0D+00, &
    2.0D+00, &
    3.0D+00, &
    4.0D+00, &
    2.0D+00, &
    3.0D+00, &
    4.0D+00, &
    2.0D+00, &
    3.0D+00, &
    4.0D+00 /)
  integer n_data
  real ( kind = rk ) x
  real ( kind = rk ), save, dimension ( n_max ) :: x_vec = (/ &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
     3.00D+00, &
    15.00D+00, &
    15.00D+00, &
    15.00D+00, &
     0.05D+00, &
     0.05D+00, &
     0.05D+00, &
     4.00D+00, &
     4.00D+00, &
     4.00D+00, &
     5.00D+00, &
     5.00D+00, &
     5.00D+00, &
     6.00D+00, &
     6.00D+00, &
     6.00D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    df = 0
    lambda = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    df = df_vec(n_data)
    lambda = lambda_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function tnc ( t, df, delta, ifault )

!*****************************************************************************80
!
!! tnc() computes the tail of the noncentral T distribution.
!
!  Discussion:
!
!    This routine computes the cumulative probability at T of the
!    non-central T-distribution with DF degrees of freedom (which may
!    be fractional) and non-centrality parameter DELTA.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by Russell Lenth.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Russell Lenth,
!    Algorithm AS 243:
!    Cumulative Distribution Function of the Non-Central T Distribution,
!    Applied Statistics,
!    Volume 38, Number 1, 1989, pages 185-189.
!
!    William Guenther,
!    Evaluation of probabilities for the noncentral distributions and
!    difference of two T-variables with a desk calculator,
!    Journal of Statistical Computation and Simulation,
!    Volume 6, Number 3-4, 1978, pages 199-206.
!
!  Parameters:
!
!    Input, real ( kind = rk ) T, the point whose cumulative probability
!    is desired.
!
!    Input, real ( kind = rk ) DF, the number of degrees of freedom.
!
!    Input, real ( kind = rk ) DELTA, the noncentrality parameter.
!
!    Output, integer IFAULT, error flag.
!    0, no error.
!    nonzero, an error occcurred.
!
!    Output, real ( kind = rk ) TNC, the tail of the noncentral
!    T distribution.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) albeta
  real ( kind = rk ) alngam
  real ( kind = rk ) alnorm
  real ( kind = rk ), parameter :: alnrpi = 0.57236494292470008707D+00
  real ( kind = rk ) b
  real ( kind = rk ) betain
  real ( kind = rk ) del
  real ( kind = rk ) delta
  real ( kind = rk ) df
  real ( kind = rk ) en
  real ( kind = rk ) errbd
  real ( kind = rk ), parameter :: errmax = 1.0D-10
  real ( kind = rk ) geven
  real ( kind = rk ) godd
  integer ifault
  integer, parameter :: itrmax = 100
  real ( kind = rk ) lambda
  logical negdel
  real ( kind = rk ) p
  real ( kind = rk ) q
  real ( kind = rk ), parameter :: r2pi = 0.79788456080286535588D+00
  real ( kind = rk ) rxb
  real ( kind = rk ) s
  real ( kind = rk ) t
  real ( kind = rk ) tnc
  real ( kind = rk ) tt
  real ( kind = rk ) x
  real ( kind = rk ) xeven
  real ( kind = rk ) xodd

  tnc = 0.0D+00

  if ( df <= 0.0D+00 ) then
    ifault = 2
    return
  end if

  ifault = 0

  tt = t
  del = delta
  negdel = .false.

  if ( t < 0.0D+00 ) then
    negdel = .true.
    tt = - tt
    del = - del
  end if
!
!  Initialize twin series.
!
  en = 1.0D+00
  x = t * t / ( t * t + df )

  if ( x <= 0.0D+00 ) then

    ifault = 0
    tnc = tnc + alnorm ( del, .true. )

    if ( negdel ) then
      tnc = 1.0D+00 - tnc
    end if

    return

  end if

  lambda = del * del
  p = 0.5D+00 * exp ( - 0.5D+00 * lambda )
  q = r2pi * p * del
  s = 0.5D+00 - p
  a = 0.5D+00
  b = 0.5D+00 * df
  rxb = ( 1.0D+00 - x )**b
  albeta = alnrpi + alngam ( b, ifault ) - alngam ( a + b, ifault )
  xodd = betain ( x, a, b, albeta, ifault )
  godd = 2.0D+00 * rxb * exp ( a * log ( x ) - albeta )
  xeven = 1.0D+00 - rxb
  geven = b * x * rxb
  tnc = p * xodd + q * xeven
!
!  Repeat until convergence.
!
  do

    a = a + 1.0D+00
    xodd = xodd - godd
    xeven = xeven - geven
    godd = godd * x * ( a + b - 1.0D+00 ) / a
    geven = geven * x * ( a + b - 0.5D+00 ) / ( a + 0.5D+00 )
    p = p * lambda / ( 2.0D+00 * en )
    q = q * lambda / ( 2.0D+00 * en + 1.0D+00 )
    s = s - p
    en = en + 1.0D+00
    tnc = tnc + p * xodd + q * xeven
    errbd = 2.0D+00 * s * ( xodd - godd )

    if ( errbd <= errmax ) then
      ifault = 0
      exit
    end if

    if ( itrmax < en ) then
      ifault = 1
      exit
    end if

  end do

  tnc = tnc + alnorm ( del, .true. )

  if ( negdel ) then
    tnc = 1.0D+00 - tnc
  end if

  return
end
