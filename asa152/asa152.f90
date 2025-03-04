function alnfac ( n )

!*****************************************************************************80
!
!! alnfac() computes the logarithm of the factorial of N.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 August 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the argument of the factorial.
!
!  Output:
!
!    real ( kind = rk ) ALNFAC, the logarithm of the factorial of N.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) alnfac
  real ( kind = rk ) alngam
  integer ier
  integer n

  alnfac = alngam ( real ( n + 1, kind = rk ), ier )

  return
end
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
!    28 August 2021
!
!  Author:
!
!    Original FORTRAN77 version by Allan Macleod
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Allan Macleod,
!    Algorithm AS 245,
!    A Robust and Reliable Algorithm for the Logarithm of the Gamma Function,
!    Applied Statistics,
!    Volume 38, Number 2, 1989, pages 397-402.
!
!  Input:
!
!    real ( kind = rk ) XVALUE, the argument of the Gamma function.
!
!  Output:
!
!    integer IFAULT, error flag.
!    0, no error occurred.
!    1, XVALUE is less than or equal to 0.
!    2, XVALUE is too big.
!
!    real ( kind = rk ) ALNGAM, the logarithm of the gamma function of X.
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
!    28 August 2021
!
!  Author:
!
!    Original FORTRAN77 version by David Hill
!    FORTRAN90 version by John Burkardt
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
function chyper ( point, kk, ll, mm, nn, ifault )

!*****************************************************************************80
!
!! chyper() computes point or cumulative hypergeometric probabilities.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 August 2021
!
!  Author:
!
!    Original FORTRAN77 version by Richard Lund
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    PR Freeman,
!    Algorithm AS 59:
!    Hypergeometric Probabilities,
!    Applied Statistics,
!    Volume 22, Number 1, 1973, pages 130-133.
!
!    Richard Lund,
!    Algorithm AS 152:
!    Cumulative hypergeometric probabilities,
!    Applied Statistics,
!    Volume 29, Number 2, 1980, pages 221-223.
!
!    BL Shea,
!    Remark AS R77:
!    A Remark on Algorithm AS 152: Cumulative hypergeometric probabilities,
!    Applied Statistics,
!    Volume 38, Number 1, 1989, pages 199-204.
!
!  Input:
!
!    logical POINT, is TRUE if the point probability is desired,
!    and FALSE if the cumulative probability is desired.
!
!    integer KK, the sample size.
!    0 <= KK <= MM.
!
!    integer LL, the number of successes in the sample.
!    0 <= LL <= KK.
!
!    integer MM, the population size that was sampled.
!    0 <= MM.
!
!    integer NN, the number of "successes" in 
!    the population.  0 <= NN <= MM.
!
!  Output:
!
!    integer IFAULT, error flag.
!    0, no error occurred.
!    nonzero, an error occurred.
!
!    real ( kind = rk ) CHYPER, the PDF (point probability) of
!    exactly LL successes out of KK samples, or the CDF (cumulative
!    probability) of up to LL successes out of KK samples.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) alnfac
  real ( kind = rk ) alnorm
  real ( kind = rk ) arg
  real ( kind = rk ) chyper
  logical dir
  real ( kind = rk ), parameter :: elimit = - 88.0D+00
  integer i
  integer ifault
  integer j
  integer k
  integer kk
  integer kl
  integer l
  integer ll
  integer m
  integer, parameter :: mbig = 600
  real ( kind = rk ) mean
  integer mm
  integer mnkl
  integer, parameter :: mvbig = 1000
  integer n
  integer nl
  integer nn
  real ( kind = rk ) p
  logical point
  real ( kind = rk ) pt
  real ( kind = rk ), parameter :: rootpi = 2.506628274631001D+00
  real ( kind = rk ), parameter :: scale = 1.0D+35
  real ( kind = rk ) sig

  ifault = 0

  k = kk + 1
  l = ll + 1
  m = mm + 1
  n = nn + 1

  dir = .true.
!
!  Check arguments are within permitted limits.
!
  chyper = 0.0D+00

  if ( n < 1 .or. m < n .or. k < 1 .or. m < k ) then
    ifault = 1
    return
  end if

  if ( l < 1 .or. m - n < k - l ) then
    ifault = 2
    return
  end if

  if ( .not. point ) then
    chyper = 1.0D+00
  end if

  if ( n < l .or. k < l ) then
    ifault = 2
    return
  end if

  ifault = 0
  chyper = 1.0D+00

  if ( k == 1 .or. k == m .or. n == 1 .or. n == m ) then
    return
  end if

  if ( .not. point .and. ll == min ( kk, nn ) ) then
    return
  end if

  p = real ( nn, kind = rk ) / real ( mm - nn, kind = rk )

  if ( 16.0D+00 * max ( p, 1.0D+00 / p ) &
    < real ( min ( kk, mm - kk ), kind = rk ) .and. &
    mvbig < mm .and. - 100.0D+00 < elimit ) then
!
!  Use a normal approximation.
!
    mean = real ( kk * nn, kind = rk ) / real ( mm, kind = rk )

    sig = sqrt ( mean * ( real ( mm - nn, kind = rk ) / real ( mm, kind = rk ) ) &
    * ( real ( mm - kk, kind = rk ) / ( real ( mm - 1, kind = rk ) ) ) )

    if ( point ) then

      arg = - 0.5D+00 * ((( real ( ll, kind = rk ) - mean ) / sig )**2 )
      if ( elimit <= arg ) then
        chyper = exp ( arg ) / ( sig * rootpi )
      else
        chyper = 0.0D+00
      end if

    else

      chyper = alnorm ( ( real ( ll, kind = rk ) + 0.5D+00 - mean ) / sig, &
        .false. )

    end if

  else
!
!  Calculate exact hypergeometric probabilities.
!  Interchange K and N if this saves calculations.
!
    if ( min ( n - 1, m - n ) < min ( k - 1, m - k ) ) then
      i = k
      k = n
      n = i
    end if

    if ( m - k < k - 1 ) then
      dir = .not. dir
      l = n - l + 1
      k = m - k + 1
    end if

    if ( mbig < mm ) then
!
!  Take logarithms of factorials.
!
      p = alnfac ( nn ) &
        - alnfac ( mm ) &
        + alnfac ( mm - kk ) &
        + alnfac ( kk ) &
        + alnfac ( mm - nn ) &
        - alnfac ( ll ) &
        - alnfac ( nn - ll ) &
        - alnfac ( kk - ll ) &
        - alnfac ( mm - nn - kk + ll )

      if ( elimit <= p ) then
        chyper = exp ( p )
      else
        chyper = 0.0D+00
      end if

    else
!
!  Use Freeman/Lund algorithm.
!
      do i = 1, l-1
        chyper = chyper * real ( ( k - i ) * ( n - i ), kind = rk ) &
        / real ( ( l - i ) * ( m - i ), kind = rk )
      end do

      if ( l /= k ) then
        j = m - n + l
        do i = l, k-1
          chyper = chyper * real ( j - i, kind = rk ) / real ( m - i, kind = rk )
        end do

      end if

    end if

    if ( point ) then
      return
    end if

    if ( chyper == 0.0D+00 ) then
!
!  We must recompute the point probability since it has underflowed.
!
      if ( mm <= mbig ) then
        p = alnfac ( nn ) &
          - alnfac ( mm ) &
          + alnfac ( kk ) &
          + alnfac ( mm - nn ) &
          - alnfac ( ll ) &
          - alnfac ( nn - ll ) &
          - alnfac ( kk - ll ) &
          - alnfac ( mm - nn - kk + ll ) &
          + alnfac ( mm - kk )
      end if

      p = p + log ( scale )

      if ( p < elimit ) then
        ifault = 3
        if ( real ( nn * kk + nn + kk + 1, kind = rk ) &
          / real ( mm + 2, kind = rk ) < real ( ll, kind = rk ) ) then
          chyper = 1.0D+00
        end if
        return
      else
        p = exp ( p )
      end if

    else
!
!  Scale up at this point.
!
      p = chyper * scale

    end if

    pt = 0.0D+00
    nl = n - l
    kl = k - l
    mnkl = m - n - kl + 1

    if ( l <= kl ) then

      do i = 1, l - 1
        p = p * real ( ( l - i ) * ( mnkl - i ), kind = rk ) / &
        real ( ( nl + i ) * ( kl + i ), kind = rk )
        pt = pt + p
      end do

    else

      dir = .not. dir
      do j = 0, kl - 1
        p = p * real ( ( nl - j ) * ( kl - j ), kind = rk ) &
        / real ( ( l + j ) * ( mnkl + j ), kind = rk )
        pt = pt + p
      end do

    end if

    if ( p == 0.0D+00 ) then
      ifault = 3
    end if

    if ( dir ) then
      chyper = chyper + ( pt / scale )
    else
      chyper = 1.0D+00 - ( pt / scale )
    end if

  end if

  return
end
subroutine hypergeometric_cdf_values ( n_data, sam, suc, pop, n, fx )

!*****************************************************************************80
!
!! hypergeometric_cdf_values() returns some values of the hypergeometric CDF.
!
!  Discussion:
!
!    CDF(X)(A,B) is the probability of at most X successes in A trials,
!    given that the probability of success on a single trial is B.
!
!    In Mathematica, the function can be evaluated by:
!
!      dist = HypergeometricDistribution [ sam, suc, pop ]
!      CDF [ dist, n ]
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 August 2021
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
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996,
!    ISBN: 0-8493-2479-3,
!    LC: QA47.M315.
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
!    integer SAM, SUC, POP, the sample size,
!    success size, and population parameters of the function.
!
!    integer N, the argument of the function.
!
!    real ( kind = rk ) FX, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 16

  real ( kind = rk ) fx
  real ( kind = rk ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.6001858177500578D-01, &
    0.2615284665839845D+00, &
    0.6695237889132748D+00, &
    0.1000000000000000D+01, &
    0.1000000000000000D+01, &
    0.5332595856827856D+00, &
    0.1819495964117640D+00, &
    0.4448047017527730D-01, &
    0.9999991751316731D+00, &
    0.9926860896560750D+00, &
    0.8410799901444538D+00, &
    0.3459800113391901D+00, &
    0.0000000000000000D+00, &
    0.2088888139634505D-02, &
    0.3876752992448843D+00, &
    0.9135215248834896D+00 /)
  integer n
  integer n_data
  integer, save, dimension ( n_max ) :: n_vec = (/ &
     7,  8,  9, 10, &
     6,  6,  6,  6, &
     6,  6,  6,  6, &
     0,  0,  0,  0 /)
  integer pop
  integer, save, dimension ( n_max ) :: pop_vec = (/ &
    100, 100, 100, 100, &
    100, 100, 100, 100, &
    100, 100, 100, 100, &
    90,  200, 1000, 10000 /)
  integer sam
  integer, save, dimension ( n_max ) :: sam_vec = (/ &
    10, 10, 10, 10, &
     6,  7,  8,  9, &
    10, 10, 10, 10, &
    10, 10, 10, 10 /)
  integer suc
  integer, save, dimension ( n_max ) :: suc_vec = (/ &
    90, 90, 90, 90, &
    90, 90, 90, 90, &
    10, 30, 50, 70, &
    90, 90, 90, 90 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    sam = 0
    suc = 0
    pop = 0
    n = 0
    fx = 0.0D+00
  else
    sam = sam_vec(n_data)
    suc = suc_vec(n_data)
    pop = pop_vec(n_data)
    n = n_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine hypergeometric_pdf_values ( n_data, sam, suc, pop, n, fx )

!*****************************************************************************80
!
!! hypergeometric_pdf_values() returns some values of the hypergeometric PDF.
!
!  Discussion:
!
!    PDF(X)(A,B) is the probability of X successes in A trials,
!    given that the probability of success on a single trial is B.
!
!    In Mathematica, the function can be evaluated by:
!
!      dist = HypergeometricDistribution [ sam, suc, pop ]
!      PDF [ dist, n ]
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 August 2021
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
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996,
!    ISBN: 0-8493-2479-3,
!    LC: QA47.M315.
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
!    integer SAM, SUC, POP, the sample size,
!    success size, and population parameters of the function.
!
!    integer N, the argument of the function.
!
!    real ( kind = rk ) FX, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 16

  real ( kind = rk ) fx
  real ( kind = rk ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.05179370533242827D+00, &
    0.2015098848089788D+00, &
    0.4079953223292903D+00, &
    0.3304762110867252D+00, &
    0.5223047493549780D+00, &
    0.3889503452643453D+00, &
    0.1505614239732950D+00, &
    0.03927689321042477D+00, &
    0.00003099828465518108D+00, &
    0.03145116093938197D+00, &
    0.2114132170316862D+00, &
    0.2075776621999210D+00, &
    0.0000000000000000D+00, &
    0.002088888139634505D+00, &
    0.3876752992448843D+00, &
    0.9135215248834896D+00 /)
  integer n
  integer n_data
  integer, save, dimension ( n_max ) :: n_vec = (/ &
     7,  8,  9, 10, &
     6,  6,  6,  6, &
     6,  6,  6,  6, &
     0,  0,  0,  0 /)
  integer pop
  integer, save, dimension ( n_max ) :: pop_vec = (/ &
    100, 100, 100, 100, &
    100, 100, 100, 100, &
    100, 100, 100, 100, &
    90,  200, 1000, 10000 /)
  integer sam
  integer, save, dimension ( n_max ) :: sam_vec = (/ &
    10, 10, 10, 10, &
     6,  7,  8,  9, &
    10, 10, 10, 10, &
    10, 10, 10, 10 /)
  integer suc
  integer, save, dimension ( n_max ) :: suc_vec = (/ &
    90, 90, 90, 90, &
    90, 90, 90, 90, &
    10, 30, 50, 70, &
    90, 90, 90, 90 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    sam = 0
    suc = 0
    pop = 0
    n = 0
    fx = 0.0D+00
  else
    sam = sam_vec(n_data)
    suc = suc_vec(n_data)
    pop = pop_vec(n_data)
    n = n_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end

