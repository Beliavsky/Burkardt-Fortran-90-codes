subroutine bivariate_normal_cdf_values ( n_data, x, y, r, fxy )

!*****************************************************************************80
!
!! bivariate_normal_cdf_values() returns some values of the bivariate normal CDF.
!
!  Discussion:
!
!    FXY is the probability that two variables A and B, which are
!    related by a bivariate normal distribution with correlation R,
!    respectively satisfy A <= X and B <= Y.
!
!    Mathematica can evaluate the bivariate normal CDF via the commands:
!
!      <<MultivariateStatistics`
!      cdf = CDF[MultinormalDistribution[{0,0}{{1,r},{r,1}}],{x,y}]
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 November 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    National Bureau of Standards,
!    Tables of the Bivariate Normal Distribution and Related Functions,
!    Applied Mathematics Series, Number 50, 1959.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = rk ) X, Y, the parameters of the function.
!
!    Output, real ( kind = rk ) R, the correlation value.
!
!    Output, real ( kind = rk ) FXY, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 41

  real ( kind = rk ) fxy
  real ( kind = rk ), save, dimension ( n_max ) :: fxy_vec = (/ &
  0.02260327218569867D+00, &
  0.1548729518584100D+00, &
  0.4687428083352184D+00, &
  0.7452035868929476D+00, &
  0.8318608306874188D+00, &
  0.8410314261134202D+00, &
  0.1377019384919464D+00, &
  0.1621749501739030D+00, &
  0.1827411243233119D+00, &
  0.2010067421506235D+00, &
  0.2177751155265290D+00, &
  0.2335088436446962D+00, &
  0.2485057781834286D+00, &
  0.2629747825154868D+00, &
  0.2770729823404738D+00, &
  0.2909261168683812D+00, &
  0.3046406378726738D+00, &
  0.3183113449213638D+00, &
  0.3320262544108028D+00, &
  0.3458686754647614D+00, &
  0.3599150462310668D+00, &
  0.3742210899871168D+00, &
  0.3887706405282320D+00, &
  0.4032765198361344D+00, &
  0.4162100291953678D+00, &
  0.6508271498838664D+00, &
  0.8318608306874188D+00, &
  0.0000000000000000D+00, &
  0.1666666666539970D+00, &
  0.2500000000000000D+00, &
  0.3333333333328906D+00, &
  0.5000000000000000D+00, &
  0.7452035868929476D+00, &
  0.1548729518584100D+00, &
  0.1548729518584100D+00, &
  0.06251409470431653D+00, &
  0.7452035868929476D+00, &
  0.1548729518584100D+00, &
  0.1548729518584100D+00, &
  0.06251409470431653D+00, &
  0.6337020457912916D+00 /)
  integer n_data
  real ( kind = rk ) r
  real ( kind = rk ), save, dimension ( n_max ) :: r_vec = (/ &
     0.500D+00,  0.500D+00,  0.500D+00,  0.500D+00,  0.500D+00, &
     0.500D+00, -0.900D+00, -0.800D+00, -0.700D+00, -0.600D+00, &
    -0.500D+00, -0.400D+00, -0.300D+00, -0.200D+00, -0.100D+00, &
     0.000D+00,  0.100D+00,  0.200D+00,  0.300D+00,  0.400D+00, &
     0.500D+00,  0.600D+00,  0.700D+00,  0.800D+00,  0.900D+00, &
     0.673D+00,  0.500D+00, -1.000D+00, -0.500D+00,  0.000D+00, &
     0.500D+00,  1.000D+00,  0.500D+00,  0.500D+00,  0.500D+00, &
     0.500D+00,  0.500D+00,  0.500D+00,  0.500D+00,  0.500D+00, &
     0.500D+00 /)
  real ( kind = rk ) x
  real ( kind = rk ), save, dimension ( n_max ) :: x_vec = (/ &
    -2.0D+00, -1.0D+00,  0.0D+00,  1.0D+00,  2.0D+00, &
     3.0D+00, -0.2D+00, -0.2D+00, -0.2D+00, -0.2D+00, &
    -0.2D+00, -0.2D+00, -0.2D+00, -0.2D+00, -0.2D+00, &
    -0.2D+00, -0.2D+00, -0.2D+00, -0.2D+00, -0.2D+00, &
    -0.2D+00, -0.2D+00, -0.2D+00, -0.2D+00, -0.2D+00, &
     1.0D+00,  2.0D+00,  0.0D+00,  0.0D+00,  0.0D+00, &
     0.0D+00,  0.0D+00,  1.0D+00,  1.0D+00, -1.0D+00, &
    -1.0D+00,  1.0D+00,  1.0D+00, -1.0D+00, -1.0D+00, &
     0.7071067811865475D+00 /)
  real ( kind = rk ) y
  real ( kind = rk ), save, dimension ( n_max ) :: y_vec = (/ &
     1.0D+00,  1.0D+00,  1.0D+00,  1.0D+00,  1.0D+00, &
     1.0D+00,  0.5D+00,  0.5D+00,  0.5D+00,  0.5D+00, &
     0.5D+00,  0.5D+00,  0.5D+00,  0.5D+00,  0.5D+00, &
     0.5D+00,  0.5D+00,  0.5D+00,  0.5D+00,  0.5D+00, &
     0.5D+00,  0.5D+00,  0.5D+00,  0.5D+00,  0.5D+00, &
     0.5D+00,  1.0D+00,  0.0D+00,  0.0D+00,  0.0D+00, &
     0.0D+00,  0.0D+00,  1.0D+00, -1.0D+00,  1.0D+00, &
    -1.0D+00,  1.0D+00, -1.0D+00,  1.0D+00, -1.0D+00, &
     0.7071067811865475D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    r = 0.0D+00
    x = 0.0D+00
    y = 0.0D+00
    fxy = 0.0D+00
  else
    r = r_vec(n_data)
    x = x_vec(n_data)
    y = y_vec(n_data)
    fxy = fxy_vec(n_data)
  end if

  return
end
function bivnor ( ah, ak, r )

!*****************************************************************************80
!
!! BIVNOR computes the bivariate normal CDF.
!
!  Discussion:
!
!    BIVNOR computes the probability for two normal variates X and Y
!    whose correlation is R, that AH <= X and AK <= Y.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 April 2012
!
!  Author:
!
!    Original FORTRAN77 version by Thomas Donnelly.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Thomas Donnelly,
!    Algorithm 462: Bivariate Normal Distribution,
!    Communications of the ACM,
!    October 1973, Volume 16, Number 10, page 638.
!
!  Parameters:
!
!    Input, real ( kind = rk ) AH, AK, the lower limits of integration.
!
!    Input, real ( kind = rk ) R, the correlation between X and Y.
!
!    Output, real ( kind = rk ) BIVNOR, the bivariate normal CDF.
!
!  Local Parameters:
!
!    Local, integer IDIG, the number of significant digits
!    to the right of the decimal point desired in the answer.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a2
  real ( kind = rk ) ah
  real ( kind = rk ) ak
  real ( kind = rk ) ap
  real ( kind = rk ) b
  real ( kind = rk ) bivnor
  real ( kind = rk ) cn
  real ( kind = rk ) con
  real ( kind = rk ) conex
  real ( kind = rk ) ex
  real ( kind = rk ) g2
  real ( kind = rk ) gauss
  real ( kind = rk ) gh
  real ( kind = rk ) gk
  real ( kind = rk ) gw
  real ( kind = rk ) h2
  real ( kind = rk ) h4
  integer i
  integer, parameter :: idig = 15
  integer is
  real ( kind = rk ) r
  real ( kind = rk ) rr
  real ( kind = rk ) s1
  real ( kind = rk ) s2
  real ( kind = rk ) sgn
  real ( kind = rk ) sn
  real ( kind = rk ) sp
  real ( kind = rk ) sqr
  real ( kind = rk ) t
  real ( kind = rk ), parameter :: twopi = 6.283185307179587D+00
  real ( kind = rk ) w2
  real ( kind = rk ) wh
  real ( kind = rk ) wk

  gauss ( t ) = ( 1.0D+00 + erf ( t / sqrt ( 2.0D+00 ) ) ) / 2.0D+00

!  GAUSS is a univariate lower normal tail area calculated here from the
!  central error function ERF.
!
  b = 0.0D+00

  gh = gauss ( - ah ) / 2.0D+00
  gk = gauss ( - ak ) / 2.0D+00

  if ( r == 0.0D+00 ) then
    b = 4.0D+000 * gh * gk
    b = max ( b, 0.0D+00 )
    b = min ( b, 1.0D+00 )
    bivnor = b
    return
  end if

  rr = ( 1.0D+00 + r ) * ( 1.0D+00 - r )

  if ( rr < 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BIVNOR - Fatal error!'
    write ( *, '(a)' ) '  1 < |R|.'
    stop
  end if

  if ( rr == 0.0D+00 ) then

    if ( r < 0.0D+00 ) then

      if ( ah + ak < 0.0D+00 ) then
        b = 2.0D+00 * ( gh + gk ) - 1.0D+00
      end if

    else

      if ( ah - ak < 0.0D+00 ) then
        b = 2.0D+00 * gk
      else
        b = 2.0D+00 * gh
      end if

    end if

    b = max ( b, 0.0D+00 )
    b = min ( b, 1.0D+00 )
    bivnor = b
    return

  end if

  sqr = sqrt ( rr )

  if ( idig == 15 ) then
    con = twopi * 1.0D-15 / 2.0D+00
  else
    con = twopi / 2.0D+00
    do i = 1, idig
      con = con / 10.0D+00
    end do
  end if
!
!  (0,0)
!
  if ( ah == 0.0D+00 .and. ak == 0.0D+00 ) then
    b = 0.25D+00 + asin ( r ) / twopi
    b = max ( b, 0.0D+00 )
    b = min ( b, 1.0D+00 )
    bivnor = b
    return
  end if
!
!  (0,nonzero)
!
  if ( ah == 0.0D+00 .and. ak /= 0.0D+00 ) then

    b = gk
    wh = -ak
    wk = ( ah / ak - r ) / sqr
    gw = 2.0D+00 * gk
    is = 1
!
!  (nonzero,0)
!
  else if ( ah /= 0.0D+00 .and. ak == 0.0D+00 ) then

    b = gh
    wh = -ah
    wk = ( ak / ah - r ) / sqr
    gw = 2.0D+00 * gh
    is = -1
!
!  (nonzero,nonzero)
!
  else if ( ah /= 0.0 .and. ak /= 0.0 ) then

    b = gh + gk
    if ( ah * ak < 0.0D+00 ) then
      b = b - 0.5D+00
    end if
    wh = - ah
    wk = ( ak / ah - r ) / sqr
    gw = 2.0D+00 * gh
    is = -1

  end if

  do

    sgn = -1.0D+00
    t = 0.0D+00

    if ( wk /= 0.0D+00 ) then

      if ( abs ( wk ) == 1.0D+00 ) then

        t = wk * gw * ( 1.0D+00 - gw ) / 2.0D+00
        b = b + sgn * t

      else

        if ( 1.0D+00 < abs ( wk ) ) then

          sgn = -sgn
          wh = wh * wk
          g2 = gauss ( wh )
          wk = 1.0D+00 / wk

          if ( wk < 0.0D+00 ) then
            b = b + 0.5D+00
          end if

          b = b - ( gw + g2 ) / 2.0D+00 + gw * g2

        end if

        h2 = wh * wh
        a2 = wk * wk
        h4 = h2 / 2.0D+00
        ex = exp ( - h4 )
        w2 = h4 * ex
        ap = 1.0D+00
        s2 = ap - ex
        sp = ap
        s1 = 0.0D+00
        sn = s1
        conex = abs ( con / wk )

        do

          cn = ap * s2 / ( sn + sp )
          s1 = s1 + cn

          if ( abs ( cn ) <= conex ) then
            exit
          end if

          sn = sp
          sp = sp + 1.0D+00
          s2 = s2 - w2
          w2 = w2 * h4 / sp
          ap = - ap * a2

        end do

        t = ( atan ( wk ) - wk * s1 ) / twopi
        b = b + sgn * t

      end if

    end if

    if ( 0 <= is ) then
      exit
    end if

    if ( ak == 0.0D+00 ) then
      exit
    end if

    wh = -ak
    wk = ( ah / ak - r ) / sqr
    gw = 2.0D+00 * gk
    is = 1

  end do

  b = max ( b, 0.0D+00 )
  b = min ( b, 1.0D+00 )
  bivnor = b

  return
end
subroutine bivnor_test ( )

!*****************************************************************************80
!
!! BIVNOR_TEST demonstrates the use of BIVNOR.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 April 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) bivnor
  real ( kind = rk ) fxy1
  real ( kind = rk ) fxy2
  integer n_data
  real ( kind = rk ) r
  real ( kind = rk ) x
  real ( kind = rk ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BIVNOR_TEST:'
  write ( *, '(a)' ) '  BIVNOR computes the bivariate normal probability.'
  write ( *, '(a)' ) '  Compare with tabulated values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The value R is the correlation between X and Y.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
  '          X               Y               ', &
  'R           P                         P                       DIFF'
  write ( *, '(a,a)' ) &
  '                                          ', &
  '           (Tabulated)               (BIVNOR)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bivariate_normal_cdf_values ( n_data, x, y, r, fxy1 )

    if ( n_data == 0 ) then
      exit
    end if
!
!  BIVNOR computes the "tail" of the probability, and we want the
!  initial part!  To get that value, negate X and Y.
!
    fxy2 = bivnor (  - x, - y, r )

    write ( *, '(2x,f14.6,2x,f14.6,2x,f14.6,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    x, y, r, fxy1, fxy2, abs ( fxy1 - fxy2 )

  end do

  return
end
function bivprb ( h, k, r )

!*****************************************************************************80
!
!! BIVPRB computes a bivariate normal CDF for correlated X and Y.
!
!  Discussion:
!
!    This routine computes P( H < X, K < Y ) for X and Y standard normal
!    variables with correlation R.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by Mike Patefield, David Tandy.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Mike Patefield, David Tandy,
!    Fast and Accurate Calculation of Owen's T Function,
!    Journal of Statistical Software,
!    Volume 5, Number 5, 2000, pages 1-25.
!
!  Parameters:
!
!    Input, real ( kind = rk ) H, K, the lower limits for X and Y.
!    0 <= H, 0 <= K.
!
!    Input, real ( kind = rk ) R, the correlation between X and Y.
!
!    Output, real ( kind = rk ) BIVPRB, the requested probability.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) bivprb
  real ( kind = rk ) h
  real ( kind = rk ) k
  real ( kind = rk ) q
  real ( kind = rk ) r
  real ( kind = rk ) ri
  real ( kind = rk ) rr
  real ( kind = rk ), parameter :: rroot2 = 0.70710678118654752440D+00
  real ( kind = rk ), parameter :: rtwopi = 0.15915494309189533577D+00
  real ( kind = rk ) znorm2

  if ( h < 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BIVPRB - Fatal error!'
    write ( *, '(a,g14.6)' ) '  H < 0, H = ', h
    stop
  end if

  if ( k < 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BIVPRB - Fatal error!'
    write ( *, '(a,g14.6)' ) '  K < 0, K = ', k
    stop
  end if

  if ( r == 0.0D+00 ) then

    bivprb = znorm2 ( h ) * znorm2 ( k )

  else

    rr = 1.0D+00 - r * r

    if ( 0.0D+00 < rr ) then

      ri = 1.0D+00 / sqrt ( rr )

      if ( 0.0D+00 < h .and. 0.0D+00 < k ) then

        bivprb = q ( h, ( k / h - r ) * ri ) &
               + q ( k, ( h / k - r ) * ri )
      else if ( 0.0D+00 < h ) then
        bivprb = q ( h, - r * ri )
      else if ( 0.0D+00 < k ) then
        bivprb = q ( k, - r * ri )
      else if ( h == 0.0D+00 .and. k == 0.0D+00 ) then
        bivprb = 0.25D+00 + rtwopi * asin ( r )
      end if

    else if ( 1.0D+00 <= r ) then

      if ( k <= h ) then
        bivprb = znorm2 ( h )
      else
        bivprb = znorm2 ( k )
      end if

    else

      bivprb = 0.0D+00

    end if

  end if

  return
end
subroutine normal_01_cdf_values ( n_data, x, fx )

!*****************************************************************************80
!
!! NORMAL_01_CDF_VALUES returns some values of the Normal 01 CDF.
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
subroutine owen_values ( n_data, h, a, t )

!*****************************************************************************80
!
!! OWEN_VALUES returns some values of Owen's T function.
!
!  Discussion:
!
!    Owen's T function is useful for computation of the bivariate normal
!    distribution and the distribution of a skewed normal distribution.
!
!    Although it was originally formulated in terms of the bivariate
!    normal function, the function can be defined more directly as
!
!      T(H,A) = 1 / ( 2 * pi ) *
!        Integral ( 0 <= X <= A ) e^(-H^2*(1+X^2)/2) / (1+X^2) dX
!
!    In Mathematica, the function can be evaluated by:
!
!      fx = 1/(2*Pi) * Integrate [ E^(-h^2*(1+x^2)/2)/(1+x^2), {x,0,a} ]
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Mike Patefield, David Tandy,
!    Fast and Accurate Calculation of Owen's T Function,
!    Journal of Statistical Software,
!    Volume 5, Number 5, 2000, pages 1-25.
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
!    Output, real ( kind = rk ) H, a parameter.
!
!    Output, real ( kind = rk ) A, the upper limit of the integral.
!
!    Output, real ( kind = rk ) T, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 28

  real ( kind = rk ) a
  real ( kind = rk ), save, dimension ( n_max ) :: a_vec = (/ &
    0.2500000000000000D+00, &
    0.4375000000000000D+00, &
    0.9687500000000000D+00, &
    0.0625000000000000D+00, &
    0.5000000000000000D+00, &
    0.9999975000000000D+00, &
    0.5000000000000000D+00, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.5000000000000000D+00, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.5000000000000000D+00, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.5000000000000000D+00, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.5000000000000000D+00, &
    0.1000000000000000D+01, &
    0.2000000000000000D+01, &
    0.3000000000000000D+01, &
    0.1000000000000000D+02, &
    0.1000000000000000D+03 /)
  real ( kind = rk ) h
  real ( kind = rk ), save, dimension ( n_max ) :: h_vec = (/ &
    0.0625000000000000D+00, &
    6.5000000000000000D+00, &
    7.0000000000000000D+00, &
    4.7812500000000000D+00, &
    2.0000000000000000D+00, &
    1.0000000000000000D+00, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.5000000000000000D+00, &
    0.2500000000000000D+00, &
    0.2500000000000000D+00, &
    0.2500000000000000D+00, &
    0.2500000000000000D+00, &
    0.1250000000000000D+00, &
    0.1250000000000000D+00, &
    0.1250000000000000D+00, &
    0.1250000000000000D+00, &
    0.7812500000000000D-02, &
    0.7812500000000000D-02, &
    0.7812500000000000D-02, &
    0.7812500000000000D-02, &
    0.7812500000000000D-02, &
    0.7812500000000000D-02 /)
  integer n_data
  real ( kind = rk ) t
  real ( kind = rk ), save, dimension ( n_max ) :: t_vec = (/ &
    3.8911930234701366D-02, &
    2.0005773048508315D-11, &
    6.3990627193898685D-13, &
    1.0632974804687463D-07, &
    8.6250779855215071D-03, &
    6.6741808978228592D-02, &
    0.4306469112078537D-01, &
    0.6674188216570097D-01, &
    0.7846818699308410D-01, &
    0.7929950474887259D-01, &
    0.6448860284750376D-01, &
    0.1066710629614485D+00, &
    0.1415806036539784D+00, &
    0.1510840430760184D+00, &
    0.7134663382271778D-01, &
    0.1201285306350883D+00, &
    0.1666128410939293D+00, &
    0.1847501847929859D+00, &
    0.7317273327500385D-01, &
    0.1237630544953746D+00, &
    0.1737438887583106D+00, &
    0.1951190307092811D+00, &
    0.7378938035365546D-01, &
    0.1249951430754052D+00, &
    0.1761984774738108D+00, &
    0.1987772386442824D+00, &
    0.2340886964802671D+00, &
    0.2479460829231492D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    h = 0.0D+00
    a = 0.0D+00
    t = 0.0D+00
  else
    h = h_vec(n_data)
    a = a_vec(n_data)
    t = t_vec(n_data)
  end if

  return
end
function q ( h, ah )

!*****************************************************************************80
!
!! Q computes (1/2) * p(H<Z) - T(H,AH).
!
!  Discussion:
!
!    The routine computes Q = (1/2) * P( H < Z ) - T ( H, AH ).
!
!    The result for Q is non-negative.
!
!    Warning : Q is computed as the difference between two terms;
!    When the two terms are of similar value this may produce
!    error in Q.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by Mike Patefield, David Tandy.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Mike Patefield, David Tandy,
!    Fast and Accurate Calculation of Owen's T Function,
!    Journal of Statistical Software,
!    Volume 5, Number 5, 2000, pages 1-25.
!
!  Parameters:
!
!    Input, real ( kind = rk ) H, the lower limit for Z.
!    0 < H.
!
!    Input, real ( kind = rk ) AH, one of the arguments for the T function.
!
!    Output, real ( kind = rk ) Q, the desired quantity.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) ah
  real ( kind = rk ) ahh
  real ( kind = rk ) h
  real ( kind = rk ) q
  real ( kind = rk ), parameter :: rroot2 = 0.70710678118654752440D+00
  real ( kind = rk ) t
  real ( kind = rk ) tfun
  real ( kind = rk ) znorm1
  real ( kind = rk ) znorm2

  if ( h <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Q - Fatal error!'
    write ( *, '(a,g14.6)' ) '  H <= 0.0, H = ', h
    stop
  end if

  if ( 1.0D+00 < ah ) then
    ahh = ah * h
    q = tfun ( ahh, 1.0D+00 / ah, h ) - znorm2 ( ahh ) * znorm1 ( h )
  else
    q = 0.5D+00 * znorm2 ( h ) - t ( h, ah )
  end if

  return
end
function t ( h, a )

!*****************************************************************************80
!
!! T computes Owen's T function for arbitrary H and A.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 April 2012
!
!  Author:
!
!    Original FORTRAN77 version by Mike Patefield, David Tandy.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Mike Patefield, David Tandy,
!    Fast and Accurate Calculation of Owen's T Function,
!    Journal of Statistical Software,
!    Volume 5, Number 5, 2000, pages 1-25.
!
!  Parameters:
!
!    Input, real ( kind = rk ) H, A, the arguments.
!
!    Output, real ( kind = rk ) T, the value of Owen's T function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) absa
  real ( kind = rk ) absh
  real ( kind = rk ) ah
  real ( kind = rk ), parameter :: cut = 0.67D+00
  real ( kind = rk ) h
  real ( kind = rk ) normah
  real ( kind = rk ) normh
  real ( kind = rk ), parameter :: rroot2 = 0.70710678118654752440D+00
  real ( kind = rk ) t
  real ( kind = rk ) tfun
  real ( kind = rk ) znorm1
  real ( kind = rk ) znorm2

  absh = abs ( h )
  absa = abs ( a )
  ah = absa * absh

  if ( absa <= 1.0D+00 ) then

    t = tfun ( absh, absa, ah )

  else if ( absh <= cut ) then

    t = 0.25D+00 - znorm1 ( absh ) * znorm1 ( ah ) &
      - tfun ( ah, 1.0D+00 / absa, absh )

  else

    normh = znorm2 ( absh )
    normah = znorm2 ( ah )
    t = 0.5D+00 * ( normh + normah ) - normh * normah &
    - tfun ( ah, 1.0D+00 / absa, absh )

  end if

  if ( a < 0.0D+00 ) then
    t = - t
  end if

  return
end
subroutine t_test ( )

!*****************************************************************************80
!
!! T_TEST demonstrates the use of T.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) h
  integer n_data
  real ( kind = rk ) t
  real ( kind = rk ) t1
  real ( kind = rk ) t2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'T_TEST:'
  write ( *, '(a)' ) '  T computes the Owen T function.'
  write ( *, '(a)' ) '  Compare with tabulated values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) '          H               A           ', &
  'T                         T                       DIFF'
  write ( *, '(a,a)' ) '                                     ', &
  '(Tabulated)               (TFN)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call owen_values ( n_data, h, a, t1 )

    if ( n_data == 0 ) then
      exit
    end if

    t2 = t ( h, a )

    write ( *, '(2x,f14.6,2x,f14.6,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    h, a, t1, t2, abs ( t1 - t2 )

  end do

  return
end
function tfun ( h, a, ah )

!*****************************************************************************80
!
!! TFUN computes Owen's T function for a restricted range of parameters.
!
!  Discussion:
!
!    This routine computes Owen's T-function of H and A.
!
!    This routine was originally named "TF", but was renamed "TFUN" to
!    avoid a conflict with a built in MATLAB function.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 December 2011
!
!  Author:
!
!    Original FORTRAN77 version by Mike Patefield, David Tandy.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Mike Patefield, David Tandy,
!    Fast and Accurate Calculation of Owen's T Function,
!    Journal of Statistical Software,
!    Volume 5, Number 5, 2000, pages 1-25.
!
!  Parameters:
!
!    Input, real ( kind = rk ) H, the H argument of the function.
!    0 <= H.
!
!    Input, real ( kind = rk ) A, the A argument of the function.
!    0 <= A <= 1.
!
!    Input, real ( kind = rk ) AH, the value of A*H.
!
!    Output, real ( kind = rk ) TF, the value of Owen's T function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) ah
  real ( kind = rk ) ai
  real ( kind = rk ) aj
  real ( kind = rk ), dimension ( 7 ) :: arange = (/ &
    0.025D+00, 0.09D+00, 0.15D+00, 0.36D+00, 0.5D+00, &
    0.9D+00, 0.99999D+00 /)
  real ( kind = rk ) as
  real ( kind = rk ), dimension ( 21 ) :: c2 = (/ &
                                   0.99999999999999987510D+00, &
     -0.99999999999988796462D+00,  0.99999999998290743652D+00, &
     -0.99999999896282500134D+00,  0.99999996660459362918D+00, &
     -0.99999933986272476760D+00,  0.99999125611136965852D+00, &
     -0.99991777624463387686D+00,  0.99942835555870132569D+00, &
     -0.99697311720723000295D+00,  0.98751448037275303682D+00, &
     -0.95915857980572882813D+00,  0.89246305511006708555D+00, &
     -0.76893425990463999675D+00,  0.58893528468484693250D+00, &
     -0.38380345160440256652D+00,  0.20317601701045299653D+00, &
     -0.82813631607004984866D-01,  0.24167984735759576523D-01, &
     -0.44676566663971825242D-02,  0.39141169402373836468D-03 /)
  real ( kind = rk ) dhs
  real ( kind = rk ) dj
  real ( kind = rk ) gj
  real ( kind = rk ) h
  real ( kind = rk ), dimension ( 14 ) :: hrange = (/ &
    0.02D+00, 0.06D+00, 0.09D+00, 0.125D+00, 0.26D+00, &
    0.4D+00,  0.6D+00,  1.6D+00,  1.7D+00,   2.33D+00, &
    2.4D+00,  3.36D+00, 3.4D+00,  4.8D+00 /)
  real ( kind = rk ) hs
  integer i
  integer iaint
  integer icode
  integer ihint
  integer ii
  integer j
  integer jj
  integer m
  integer maxii
  integer, dimension ( 18 ) :: meth = (/ &
    1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 4, 4, 4, 4, 5, 6 /)
  real ( kind = rk ) normh
  integer, dimension ( 18 ) :: ord = (/ &
    2, 3, 4, 5, 7,10,12,18,10,20,30,20, 4, 7, 8,20,13, 0 /)
  real ( kind = rk ), dimension ( 13 ) :: pts = (/ &
                                   0.35082039676451715489D-02, &
      0.31279042338030753740D-01,  0.85266826283219451090D-01, &
      0.16245071730812277011D+00,  0.25851196049125434828D+00, &
      0.36807553840697533536D+00,  0.48501092905604697475D+00, &
      0.60277514152618576821D+00,  0.71477884217753226516D+00, &
      0.81475510988760098605D+00,  0.89711029755948965867D+00, &
      0.95723808085944261843D+00,  0.99178832974629703586D+00 /)
  real ( kind = rk ) r
  real ( kind = rk ), parameter :: rroot2 = 0.70710678118654752440D+00
  real ( kind = rk ), parameter :: rrtpi = 0.39894228040143267794D+00
  real ( kind = rk ), parameter :: rtwopi = 0.15915494309189533577D+00
  integer, dimension ( 15, 8 ) :: select = reshape ( (/&
    1, 1, 2,13,13,13,13,13,13,13,13,16,16,16, 9, &
    1, 2, 2, 3, 3, 5, 5,14,14,15,15,16,16,16, 9, &
    2, 2, 3, 3, 3, 5, 5,15,15,15,15,16,16,16,10, &
    2, 2, 3, 5, 5, 5, 5, 7, 7,16,16,16,16,16,10, &
    2, 3, 3, 5, 5, 6, 6, 8, 8,17,17,17,12,12,11, &
    2, 3, 5, 5, 5, 6, 6, 8, 8,17,17,17,12,12,12, &
    2, 3, 4, 4, 6, 6, 8, 8,17,17,17,17,17,12,12, &
    2, 3, 4, 4, 6, 6,18,18,18,18,17,17,17,12,12 /), (/ 15, 8 /) )
  real ( kind = rk ) tfun
  real ( kind = rk ) vi
  real ( kind = rk ), dimension ( 13 ) :: wts = (/ &
                                   0.18831438115323502887D-01, &
      0.18567086243977649478D-01,  0.18042093461223385584D-01, &
      0.17263829606398753364D-01,  0.16243219975989856730D-01, &
      0.14994592034116704829D-01,  0.13535474469662088392D-01, &
      0.11886351605820165233D-01,  0.10070377242777431897D-01, &
      0.81130545742299586629D-02,  0.60419009528470238773D-02, &
      0.38862217010742057883D-02,  0.16793031084546090448D-02 /)
  real ( kind = rk ) y
  real ( kind = rk ) yi
  real ( kind = rk ) z
  real ( kind = rk ) zi
  real ( kind = rk ) znorm1
  real ( kind = rk ) znorm2
!
!  Determine appropriate method from t1...t6
!
  ihint = 15

  do i = 1, 14
    if ( h <= hrange(i) ) then
      ihint = i
      exit
    end if
  end do

  iaint = 8

  do i = 1, 7
    if ( a <= arange(i) ) then
      iaint = i
      exit
    end if
  end do

  icode = select(ihint,iaint)
  m = ord(icode)
!
!  t1(h, a, m) ; m = 2, 3, 4, 5, 7, 10, 12 or 18
!  jj = 2j - 1 ; gj = exp(-h*h/2) * (-h*h/2)**j / j!
!  aj = a**(2j-1) / (2*pi)
!
  if ( meth(icode) == 1 ) then

    hs = - 0.5D+00 * h * h
    dhs = exp ( hs )
    as = a * a
    j = 1
    jj = 1
    aj = rtwopi * a
    tfun = rtwopi * atan ( a )
    dj = dhs - 1.0D+00
    gj = hs * dhs

    do

      tfun = tfun + dj * aj / real ( jj, kind = rk )

      if ( m <= j ) then
        return
      end if

      j = j + 1
      jj = jj + 2
      aj = aj * as
      dj = gj - dj
      gj = gj * hs / real ( j, kind = rk )

    end  do
!
!  t2(h, a, m) ; m = 10, 20 or 30
!  z = (-1)**(i-1) * zi ; ii = 2i - 1
!  vi = (-1)**(i-1) * a**(2i-1) * exp[-(a*h)**2/2] / sqrt(2*pi)
!
  else if ( meth(icode) == 2 ) then

    maxii = m + m + 1
    ii = 1
    tfun = 0.0D+00
    hs = h * h
    as = - a * a
    vi = rrtpi * a * exp ( - 0.5D+00 * ah * ah )
    z = znorm1 ( ah ) / h
    y = 1.0D+00 / hs

    do

      tfun = tfun + z

      if ( maxii <= ii ) then
        tfun = tfun * rrtpi * exp ( - 0.5D+00 * hs )
        return
      end if

      z = y * ( vi - real ( ii, kind = rk ) * z )
      vi = as * vi
      ii = ii + 2

    end do
!
!  t3(h, a, m) ; m = 20
!  ii = 2i - 1
!  vi = a**(2i-1) * exp[-(a*h)**2/2] / sqrt(2*pi)
!
  else if ( meth(icode) == 3 ) then

    i = 1
    ii = 1
    tfun = 0.0D+00
    hs = h * h
    as = a * a
    vi = rrtpi * a * exp ( - 0.5D+00 * ah * ah )
    zi = znorm1 ( ah ) / h
    y = 1.0D+00 / hs

    do

      tfun = tfun + zi * c2(i)

      if ( m < i ) then
        tfun = tfun * rrtpi * exp ( - 0.5D+00 * hs )
        return
      end if

      zi = y  * ( real ( ii, kind = rk ) * zi - vi )
      vi = as * vi
      i = i + 1
      ii = ii + 2

    end do
!
!  t4(h, a, m) ; m = 4, 7, 8 or 20;  ii = 2i + 1
!  ai = a * exp[-h*h*(1+a*a)/2] * (-a*a)**i / (2*pi)
!
  else if ( meth(icode) == 4 ) then

    maxii = m + m + 1
    ii = 1
    hs = h * h
    as = - a * a
    tfun = 0.0D+00
    ai = rtwopi * a * exp ( - 0.5D+00 * hs * ( 1.0D+00 - as ) )
    yi = 1.0D+00

    do

      tfun = tfun + ai * yi

      if ( maxii <= ii ) then
        return
      end if

      ii = ii + 2
      yi = ( 1.0D+00 - hs * yi ) / real ( ii, kind = rk )
      ai = ai * as

    end do
!
!  t5(h, a, m) ; m = 13
!  2m - point gaussian quadrature
!
  else if ( meth(icode) == 5 ) then

    tfun = 0.0D+00
    as = a * a
    hs = - 0.5D+00 * h * h
    do i = 1, m
      r = 1.0D+00 + as * pts(i)
      tfun = tfun + wts(i) * exp ( hs * r ) / r
    end do
    tfun = a * tfun
!
!  t6(h, a);  approximation for a near 1, (a<=1)
!
  else if ( meth(icode) == 6 ) then

    normh = znorm2 ( h )
    tfun = 0.5D+00 * normh * ( 1.0D+00 - normh )
    y = 1.0D+00 - a
    r = atan ( y / ( 1.0D+00 + a ) )

    if ( r /= 0.0D+00 ) then
      tfun = tfun - rtwopi * r * exp ( - 0.5D+00 * y * h * h / r )
    end if

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
function znorm1 ( z )

!*****************************************************************************80
!
!! ZNORM1 evaluates the normal CDF from 0 to Z.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) Z, the upper limit.
!
!    Output, real ( kind = rk ) ZNORM1, the probability that a standard
!    normal variable will lie between 0 and Z.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) z
  real ( kind = rk ) znorm1

  znorm1 = 0.5D+00 * erf ( z / sqrt ( 2.0D+00 ) )

  return
end
subroutine znorm1_test ( )

!*****************************************************************************80
!
!! ZNORM1_TEST tests ZNORM1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 May 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) fx1
  real ( kind = rk ) fx2
  integer n_data
  real ( kind = rk ) x
  real ( kind = rk ) znorm1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ZNORM1_TEST:'
  write ( *, '(a)' ) '  ZNORM1 computes the normal CDF starting at 0.'
  write ( *, '(a)' ) '  Compare with tabulated values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
  '          X           P                         P                       DIFF'
  write ( *, '(a,a)' ) &
  '                     (Tabulated)               (ZNORM1)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call normal_01_cdf_values ( n_data, x, fx1 )

    fx1 = fx1 - 0.5D+00

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = znorm1 ( x )

    write ( *, '(2x,f14.6,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    x, fx1, fx2, abs ( fx1 - fx2 )

  end do

  return
end
function znorm2 ( z )

!*****************************************************************************80
!
!! ZNORM2 evaluates the normal CDF from Z to +oo.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) Z, the lower limit.
!
!    Output, real ( kind = rk ) ZNORM2, the probability that a standard
!    normal variable will lie between Z and +oo.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) z
  real ( kind = rk ) znorm2

  znorm2 = 0.5D+00 * erfc ( z / sqrt ( 2.0D+00 ) )

  return
end
subroutine znorm2_test ( )

!*****************************************************************************80
!
!! ZNORM2_TEST tests ZNORM2.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 November 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) fx1
  real ( kind = rk ) fx2
  integer n_data
  real ( kind = rk ) x
  real ( kind = rk ) znorm2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ZNORM2_TEST:'
  write ( *, '(a)' ) '  ZNORM2 computes the complementary normal CDF.'
  write ( *, '(a)' ) '  Compare with tabulated values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
  '          X           P                         P                       DIFF'
  write ( *, '(a,a)' ) &
  '                     (Tabulated)               (ZNORM2)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call normal_01_cdf_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx1 = 1.0D+00 - fx1
    fx2 = znorm2 ( x )

    write ( *, '(2x,f14.6,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    x, fx1, fx2, abs ( fx1 - fx2 )

  end do

  return
end

