subroutine normal_01_cdf_values ( n_data, x, fx )

!*****************************************************************************80
!
!! normal_01_cdf_values() returns some values of the Normal 01 CDF.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Needs["Statistics`ContinuousDistributions`"]
!      dist = NormalDistribution [ 0, 1 ]
!      CDF [ dist, x ]
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz and Irene Stegun,
!    Handbook of Mathematical Functions,
!    US Department of Commerce, 1964.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
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

  integer, parameter :: n_max = 17

  real ( kind = rk ) fx
  real ( kind = rk ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.5000000000000000D+00, &
    0.5398278372770290D+00, &
    0.5792597094391030D+00, &
    0.6179114221889526D+00, &
    0.6554217416103242D+00, &
    0.6914624612740131D+00, &
    0.7257468822499270D+00, &
    0.7580363477769270D+00, &
    0.7881446014166033D+00, &
    0.8159398746532405D+00, &
    0.8413447460685429D+00, &
    0.9331927987311419D+00, &
    0.9772498680518208D+00, &
    0.9937903346742239D+00, &
    0.9986501019683699D+00, &
    0.9997673709209645D+00, &
    0.9999683287581669D+00 /)
  integer n_data
  real ( kind = rk ) x
  real ( kind = rk ), save, dimension ( n_max ) :: x_vec = (/ &
    0.0000000000000000D+00, &
    0.1000000000000000D+00, &
    0.2000000000000000D+00, &
    0.3000000000000000D+00, &
    0.4000000000000000D+00, &
    0.5000000000000000D+00, &
    0.6000000000000000D+00, &
    0.7000000000000000D+00, &
    0.8000000000000000D+00, &
    0.9000000000000000D+00, &
    0.1000000000000000D+01, &
    0.1500000000000000D+01, &
    0.2000000000000000D+01, &
    0.2500000000000000D+01, &
    0.3000000000000000D+01, &
    0.3500000000000000D+01, &
    0.4000000000000000D+01 /)

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
function r8_normal_01_cdf_inverse ( p )

!*****************************************************************************80
!
!! r8_normal_01_cdf_inverse() inverts the standard normal CDF.
!
!  Discussion:
!
!    The result is accurate to about 1 part in 10^16.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 December 2004
!
!  Author:
!
!    Original FORTRAN77 version by Michael Wichura.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Michael Wichura,
!    The Percentage Points of the Normal Distribution,
!    Algorithm AS 241,
!    Applied Statistics,
!    Volume 37, Number 3, pages 477-484, 1988.
!
!  Parameters:
!
!    Input, real ( kind = rk ) P, the value of the cumulative probability 
!    densitity function.  0 < P < 1.  If P is outside this range,
!    an "infinite" value will be returned.
!
!    Output, real ( kind = rk ) D_NORMAL_01_CDF_INVERSE, the normal deviate 
!    value with the property that the probability of a standard normal 
!    deviate being less than or equal to the value is P.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), parameter, dimension ( 8 ) :: a = (/ &
    3.3871328727963666080D+00, &
    1.3314166789178437745D+02, &
    1.9715909503065514427D+03, &
    1.3731693765509461125D+04, &
    4.5921953931549871457D+04, &
    6.7265770927008700853D+04, &
    3.3430575583588128105D+04, &
    2.5090809287301226727D+03 /)
  real ( kind = rk ), parameter, dimension ( 8 ) :: b = (/ &
    1.0D+00, &
    4.2313330701600911252D+01, &
    6.8718700749205790830D+02, &
    5.3941960214247511077D+03, &
    2.1213794301586595867D+04, &
    3.9307895800092710610D+04, &
    2.8729085735721942674D+04, &
    5.2264952788528545610D+03 /)
  real ( kind = rk ), parameter, dimension ( 8 ) :: c = (/ &
    1.42343711074968357734D+00, &
    4.63033784615654529590D+00, &
    5.76949722146069140550D+00, &
    3.64784832476320460504D+00, &
    1.27045825245236838258D+00, &
    2.41780725177450611770D-01, &
    2.27238449892691845833D-02, &
    7.74545014278341407640D-04 /)
  real ( kind = rk ), parameter :: const1 = 0.180625D+00
  real ( kind = rk ), parameter :: const2 = 1.6D+00
  real ( kind = rk ), parameter, dimension ( 8 ) :: d = (/ &
    1.0D+00, &
    2.05319162663775882187D+00, &
    1.67638483018380384940D+00, &
    6.89767334985100004550D-01, &
    1.48103976427480074590D-01, &
    1.51986665636164571966D-02, &
    5.47593808499534494600D-04, &
    1.05075007164441684324D-09 /)
  real ( kind = rk ), parameter, dimension ( 8 ) :: e = (/ &
    6.65790464350110377720D+00, &
    5.46378491116411436990D+00, &
    1.78482653991729133580D+00, &
    2.96560571828504891230D-01, &
    2.65321895265761230930D-02, &
    1.24266094738807843860D-03, &
    2.71155556874348757815D-05, &
    2.01033439929228813265D-07 /)
  real ( kind = rk ), parameter, dimension ( 8 ) :: f = (/ &
    1.0D+00, &
    5.99832206555887937690D-01, &
    1.36929880922735805310D-01, &
    1.48753612908506148525D-02, &
    7.86869131145613259100D-04, &
    1.84631831751005468180D-05, &
    1.42151175831644588870D-07, &
    2.04426310338993978564D-15 /)
  real ( kind = rk ) p
  real ( kind = rk ) q
  real ( kind = rk ) r
  real ( kind = rk ) r8_normal_01_cdf_inverse 
  real ( kind = rk ) r8poly_value
  real ( kind = rk ), parameter :: split1 = 0.425D+00
  real ( kind = rk ), parameter :: split2 = 5.0D+00

  if ( p <= 0.0D+00 ) then
    r8_normal_01_cdf_inverse = - huge ( p )
    return
  end if

  if ( 1.0D+00 <= p ) then
    r8_normal_01_cdf_inverse = huge ( p )
    return
  end if

  q = p - 0.5D+00

  if ( abs ( q ) <= split1 ) then

    r = const1 - q * q
    r8_normal_01_cdf_inverse = q * r8poly_value ( 8, a, r ) &
                                 / r8poly_value ( 8, b, r )

  else

    if ( q < 0.0D+00 ) then
      r = p
    else
      r = 1.0D+00 - p
    end if

    if ( r <= 0.0D+00 ) then
      r8_normal_01_cdf_inverse = - 1.0D+00
      stop
    end if

    r = sqrt ( -log ( r ) )

    if ( r <= split2 ) then

      r = r - const2
      r8_normal_01_cdf_inverse = r8poly_value ( 8, c, r ) &
                               / r8poly_value ( 8, d, r )

    else

      r = r - split2
      r8_normal_01_cdf_inverse = r8poly_value ( 8, e, r ) &
                               / r8poly_value ( 8, f, r )
   
    end if

    if ( q < 0.0D+00 ) then
      r8_normal_01_cdf_inverse = - r8_normal_01_cdf_inverse
    end if

  end if

  return
end
function r8poly_value ( n, a, x )

!*****************************************************************************80
!
!! r8poly_value() evaluates an R8POLY
!
!  Discussion:
!
!    For sanity's sake, the value of N indicates the NUMBER of 
!    coefficients, or more precisely, the ORDER of the polynomial,
!    rather than the DEGREE of the polynomial.  The two quantities
!    differ by 1, but cause a great deal of confusion.
!
!    Given N and A, the form of the polynomial is:
!
!      p(x) = a(1) + a(2) * x + ... + a(n-1) * x^(n-2) + a(n) * x^(n-1)
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the polynomial.
!
!    Input, real ( kind = rk ) A(N), the coefficients of the polynomial.
!    A(1) is the constant term.
!
!    Input, real ( kind = rk ) X, the point at which the polynomial is 
!    to be evaluated.
!
!    Output, real ( kind = rk ) R8POLY_VALUE, the value of the polynomial at X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  integer i
  real ( kind = rk ) r8poly_value
  real ( kind = rk ) x

  r8poly_value = 0.0D+00
  do i = n, 1, -1
    r8poly_value = r8poly_value * x + a(i)
  end do

  return
end

