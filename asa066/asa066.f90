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
!    27 August 2021
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
!  Input:
!
!    real ( kind = rk ) X, is one endpoint of the semi-infinite interval
!    over which the integration takes place.
!
!    logical UPPER, determines whether the upper or lower
!    interval is to be integrated:
!    .TRUE.  => integrate from X to + Infinity;
!    .FALSE. => integrate from - Infinity to X.
!
!  Output:
!
!    real ( kind = rk ) ALNORM, the integral of the standard normal
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
!    real ( kind = rk ) X, the argument of the function.
!
!    real ( kind = rk ) FX, the value of the function.
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
subroutine normp ( z, p, q, pdf )

!*****************************************************************************80
!
!! normp() computes the cumulative density of the standard normal distribution.
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
!    Original FORTRAN77 version by Alan Miller.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thacher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Input:
!
!    real ( kind = rk ) Z, divides the real line into two
!    semi-infinite intervals, over each of which the standard normal
!    distribution is to be integrated.
!
!  Output:
!
!    real ( kind = rk ) P, Q, the integrals of the standard normal
!    distribution over the intervals ( - Infinity, Z] and
!    [Z, + Infinity ), respectively.
!
!    real ( kind = rk ) PDF, the value of the standard normal
!    distribution at Z.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) :: cutoff = 7.071D+00
  real ( kind = rk ) expntl
  real ( kind = rk ) p
  real ( kind = rk ) :: p0 = 220.2068679123761D+00
  real ( kind = rk ) :: p1 = 221.2135961699311D+00
  real ( kind = rk ) :: p2 = 112.0792914978709D+00
  real ( kind = rk ) :: p3 = 33.91286607838300D+00
  real ( kind = rk ) :: p4 = 6.373962203531650D+00
  real ( kind = rk ) :: p5 = 0.7003830644436881D+00
  real ( kind = rk ) :: p6 = 0.03526249659989109D+00
  real ( kind = rk ) pdf
  real ( kind = rk ) q
  real ( kind = rk ) :: q0 = 440.4137358247522D+00
  real ( kind = rk ) :: q1 = 793.8265125199484D+00
  real ( kind = rk ) :: q2 = 637.3336333788311D+00
  real ( kind = rk ) :: q3 = 296.5642487796737D+00
  real ( kind = rk ) :: q4 = 86.78073220294608D+00
  real ( kind = rk ) :: q5 = 16.06417757920695D+00
  real ( kind = rk ) :: q6 = 1.755667163182642D+00
  real ( kind = rk ) :: q7 = 0.08838834764831844D+00
  real ( kind = rk ) :: root2pi = 2.506628274631001D+00
  real ( kind = rk ) z
  real ( kind = rk ) zabs

  zabs = abs ( z )
!
!  37 < |Z|.
!
  if ( 37.0D+00 < zabs ) then

    pdf = 0.0D+00
    p = 0.0D+00
!
!  |Z| <= 37.
!
  else

    expntl = exp ( - 0.5D+00 * zabs * zabs )
    pdf = expntl / root2pi
!
!  |Z| < CUTOFF = 10 / sqrt(2).
!
    if ( zabs < cutoff ) then

      p = expntl * (((((( &
          p6   * zabs &
        + p5 ) * zabs &
        + p4 ) * zabs &
        + p3 ) * zabs &
        + p2 ) * zabs &
        + p1 ) * zabs &
        + p0 ) / ((((((( &
          q7   * zabs &
        + q6 ) * zabs &
        + q5 ) * zabs &
        + q4 ) * zabs &
        + q3 ) * zabs &
        + q2 ) * zabs &
        + q1 ) * zabs &
      + q0 )
!
!  CUTOFF <= |Z|.
!
    else

      p = pdf / ( &
        zabs + 1.0D+00 / ( &
        zabs + 2.0D+00 / ( &
        zabs + 3.0D+00 / ( &
        zabs + 4.0D+00 / ( &
        zabs + 0.65D+00 )))))

    end if

  end if

  if ( z < 0.0D+00 ) then
    q = 1.0D+00 - p
  else
    q = p
    p = 1.0D+00 - q
  end if

  return
end
subroutine nprob ( z, p, q, pdf )

!*****************************************************************************80
!
!! nprob() computes the cumulative density of the standard normal distribution.
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
!    Original FORTRAN77 version by AG Adams.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    AG Adams,
!    Algorithm 39:
!    Areas Under the Normal Curve,
!    Computer Journal,
!    Volume 12, Number 2, May 1969, pages 197-198.
!
!  Input:
!
!    real ( kind = rk ) Z, divides the real line into
!    two semi-infinite intervals, over each of which the standard normal
!    distribution is to be integrated.
!
!  Output:
!
!    real ( kind = rk ) P, Q, the integrals of the standard normal
!    distribution over the intervals ( - Infinity, Z] and
!    [Z, + Infinity ), respectively.
!
!    real ( kind = rk ) PDF, the value of the standard normal
!    distribution at Z.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), parameter :: a0 = 0.5D+00
  real ( kind = rk ), parameter :: a1 = 0.398942280444D+00
  real ( kind = rk ), parameter :: a2 = 0.399903438504D+00
  real ( kind = rk ), parameter :: a3 = 5.75885480458D+00
  real ( kind = rk ), parameter :: a4 = 29.8213557808D+00
  real ( kind = rk ), parameter :: a5 = 2.62433121679D+00
  real ( kind = rk ), parameter :: a6 = 48.6959930692D+00
  real ( kind = rk ), parameter :: a7 = 5.92885724438D+00
  real ( kind = rk ), parameter :: b0 = 0.398942280385D+00
  real ( kind = rk ), parameter :: b1 = 0.000000038052D+00
  real ( kind = rk ), parameter :: b2 = 1.00000615302D+00
  real ( kind = rk ), parameter :: b3 = 0.000398064794D+00
  real ( kind = rk ), parameter :: b4 = 1.98615381364D+00
  real ( kind = rk ), parameter :: b5 = 0.151679116635D+00
  real ( kind = rk ), parameter :: b6 = 5.29330324926D+00
  real ( kind = rk ), parameter :: b7 = 4.8385912808D+00
  real ( kind = rk ), parameter :: b8 = 15.1508972451D+00
  real ( kind = rk ), parameter :: b9 = 0.742380924027D+00
  real ( kind = rk ), parameter :: b10 = 30.789933034D+00
  real ( kind = rk ), parameter :: b11 = 3.99019417011D+00
  real ( kind = rk ) p
  real ( kind = rk ) pdf
  real ( kind = rk ) q
  real ( kind = rk ) y
  real ( kind = rk ) z
  real ( kind = rk ) zabs

  zabs = abs ( z )
!
!  |Z| between 0 and 1.28
!
  if ( abs ( z ) <= 1.28D+00 ) then

    y = a0 * z * z
    pdf = exp ( - y ) * b0

    q = a0 - zabs * ( a1 - a2 * y &
      / ( y + a3 - a4 &
      / ( y + a5 + a6 &
      / ( y + a7 ))))
!
!  |Z| between 1.28 and 12.7
!
  else if ( abs ( z ) <= 12.7D+00 ) then

    y = a0 * z * z
    pdf = exp ( - y ) * b0

    q = pdf &
      / ( zabs - b1 + b2 &
      / ( zabs + b3 + b4 &
      / ( zabs - b5 + b6 &
      / ( zabs + b7 - b8 &
      / ( zabs + b9 + b10 &
      / ( zabs + b11 ))))))
!
!  Z far out in tail.
!
  else

    q = 0.0D+00
    pdf = 0.0D+00

  end if

  if ( z < 0.0D+00 ) then
    p = q
    q = 1.0D+00 - p
  else
    p = 1.0D+00 - q
  end if

  return
end

