subroutine log_normal_cdf ( x, mu, sigma, cdf )

!*****************************************************************************80
!
!! LOG_NORMAL_CDF evaluates the Log Normal CDF.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, the argument of the PDF.
!    0.0 < X.
!
!    Input, real ( kind = rk ) MU, SIGMA, the parameters of the PDF.
!    0.0 < SIGMA.
!
!    Output, real ( kind = rk ) CDF, the value of the CDF.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) cdf
  real ( kind = rk ) logx
  real ( kind = rk ) mu
  real ( kind = rk ) sigma
  real ( kind = rk ) x

  if ( x <= 0.0D+00 ) then

    cdf = 0.0D+00

  else

    logx = log ( x )

    call normal_cdf ( logx, mu, sigma, cdf )

  end if

  return
end
subroutine log_normal_cdf_inv ( cdf, mu, sigma, x )

!*****************************************************************************80
!
!! LOG_NORMAL_CDF_INV inverts the Log Normal CDF.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!    Input, real ( kind = rk ) MU, SIGMA, the parameters of the PDF.
!    0.0 < SIGMA.
!
!    Input, real ( kind = rk ) X, the corresponding argument.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) cdf
  real ( kind = rk ) logx
  real ( kind = rk ) mu
  real ( kind = rk ) sigma
  real ( kind = rk ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LOG_NORMAL_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop 1
  end if

  call normal_cdf_inv ( cdf, mu, sigma, logx )

  x = exp ( logx )

  return
end
subroutine log_normal_pdf ( x, mu, sigma, pdf )

!*****************************************************************************80
!
!! LOG_NORMAL_PDF evaluates the Log Normal PDF.
!
!  Discussion:
!
!    PDF(A,BX)
!      = exp ( - 0.5 * ( ( log ( X ) - MU ) / SIGMA )^2 )
!        / ( SIGMA * X * sqrt ( 2 * PI ) )
!
!    The Log Normal PDF is also known as the Cobb-Douglas PDF,
!    and as the Antilog_normal PDF.
!
!    The Log Normal PDF describes a variable X whose logarithm
!    is normally distributed.
!
!    The special case MU = 0, SIGMA = 1 is known as Gilbrat's PDF.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, the argument of the PDF.
!    0.0 < X
!
!    Input, real ( kind = rk ) MU, SIGMA, the parameters of the PDF.
!    0.0 < SIGMA.
!
!    Output, real ( kind = rk ) PDF, the value of the PDF.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) mu
  real ( kind = rk ) pdf
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) sigma
  real ( kind = rk ) x

  if ( x <= 0.0D+00 ) then
    pdf = 0.0D+00
  else
    pdf = exp ( - 0.5D+00 * ( ( log ( x ) - mu ) / sigma ) ** 2 ) &
      / ( sigma * x * sqrt ( 2.0D+00 * r8_pi ) )
  end if

  return
end
subroutine log_normal_truncated_ab_cdf ( x, mu, sigma, a, b, cdf )

!*****************************************************************************80
!
!! LOG_NORMAL_TRUNCATED_AB_CDF evaluates the truncated AB Lognormal CDF.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 March 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, the argument of the PDF.
!    0.0 < X.
!
!    Input, real ( kind = rk ) MU, SIGMA, the parameters of the log normal PDF.
!    0.0 < SIGMA.
!
!    Input, real ( kind = rk ) A, B, the lower and upper truncation limits.
!    A < B.
!
!    Output, real ( kind = rk ) CDF, the value of the CDF.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  logical ( kind = 4 ) check
  real ( kind = rk ) cdf
  real ( kind = rk ) lncdf_a
  real ( kind = rk ) lncdf_b
  real ( kind = rk ) lncdf_x
  real ( kind = rk ) mu
  real ( kind = rk ) sigma
  real ( kind = rk ) x

  call log_normal_truncated_ab_check ( mu, sigma, a, b, check )

  if ( .not. check ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'LOG_NORMAL_TRUNCATED_AB_CDF - Fatal error!'
    write ( *, '(a)' ) '  Parameters are not legal.'
  end if

  if ( x <= a ) then

    cdf = 0.0D+00

  else if ( b <= x ) then

    cdf = 1.0D+00

  else

    call log_normal_cdf ( a, mu, sigma, lncdf_a )
    call log_normal_cdf ( b, mu, sigma, lncdf_b )
    call log_normal_cdf ( x, mu, sigma, lncdf_x )

    cdf = ( lncdf_x - lncdf_a ) / ( lncdf_b - lncdf_a )

  end if

  return
end
subroutine log_normal_truncated_ab_cdf_inv ( cdf, mu, sigma, a, b, x )

!*****************************************************************************80
!
!! LOG_NORMAL_TRUNCATED_AB_CDF_INV inverts the truncated AB Lognormal CDF.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 March 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) CDF, the value of the CDF.
!    0 < CDF < 1.
!
!    Input, real ( kind = rk ) MU, SIGMA, the parameters of the log normal PDF.
!    0.0 < SIGMA.
!
!    Input, real ( kind = rk ) A, B, the lower and upper truncation limits.
!    A < B.
!
!    Output, real ( kind = rk ) X, the argument of the PDF with the given CDF.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) cdf
  logical ( kind = 4 ) check
  real ( kind = rk ) lncdf_a
  real ( kind = rk ) lncdf_b
  real ( kind = rk ) lncdf_x
  real ( kind = rk ) mu
  real ( kind = rk ) sigma
  real ( kind = rk ) x

  call log_normal_truncated_ab_check ( mu, sigma, a, b, check )

  if ( .not. check ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'LOG_NORMAL_TRUNCATED_AB_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  Parameters are not legal.'
  end if

  if ( cdf <= 0.0D+00 ) then

    x = a

  else if ( 1.0D+00 <= cdf ) then

    x = b

  else

    call log_normal_cdf ( a, mu, sigma, lncdf_a )
    call log_normal_cdf ( b, mu, sigma, lncdf_b )

    lncdf_x = lncdf_a + cdf * ( lncdf_b - lncdf_a )
    call log_normal_cdf_inv ( lncdf_x, mu, sigma, x )

  end if

  return
end
subroutine log_normal_truncated_ab_check ( mu, sigma, a, b, check )

!*****************************************************************************80
!
!! LOG_NORMAL_TRUNCATED_AB_CHECK: check truncated normal AB Lognormal PDF.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 March 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) MU, SIGMA, the parameters of the log normal PDF.
!    0.0 < SIGMA.
!
!    Input, real ( kind = rk ) A, B, the lower and upper truncation limits.
!    A < B.
!
!    Output, logical CHECK, is true if the parameters are legal.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  logical ( kind = 4 ) check
  real ( kind = rk ) mu
  real ( kind = rk ) sigma

  check = .true.

  call r8_fake_use ( mu )

  if ( sigma <= 0.0D+00 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'LOG_NORMAL_TRUNCATED_AB_CHECK - Fatal error!'
    write ( *, '(a)' ) '  SIGMA <= 0.'
    check = .false.
  end if

  if ( b <= a ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'LOG_NORMAL_TRUNCATED_AB_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B <= A.'
    check = .false.
  end if

  return
end
subroutine log_normal_truncated_ab_mean ( mu, sigma, a, b, mean )

!*****************************************************************************80
!
!! LOG_NORMAL_TRUNCATED_AB_MEAN: mean of the Log Normal Truncated AB PDF.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 March 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) MU, SIGMA, the parameters of the log normal PDF.
!    0.0 < SIGMA.
!
!    Input, real ( kind = rk ) A, B, the lower and upper truncation limits.
!    A < B.
!
!    Output, real ( kind = rk ) MEAN, the mean of the PDF.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) a0
  real ( kind = rk ) b
  real ( kind = rk ) b0
  real ( kind = rk ) c1
  real ( kind = rk ) c2
  real ( kind = rk ) c3
  real ( kind = rk ) c4
  logical ( kind = 4 ) check
  real ( kind = rk ) ln_mean
  real ( kind = rk ) mean
  real ( kind = rk ) mu
  real ( kind = rk ) sigma

  call log_normal_truncated_ab_check ( mu, sigma, a, b, check )

  if ( .not. check ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'LOG_NORMAL_TRUNCATED_AB_MEAN - Fatal error!'
    write ( *, '(a)' ) '  Parameters are not legal.'
  end if

  a0 = ( log ( a ) - mu ) / sigma
  b0 = ( log ( b ) - mu ) / sigma

  call normal_01_cdf ( sigma - a0, c1 )
  call normal_01_cdf ( sigma - b0, c2 )
  call normal_01_cdf ( + a0, c3 )
  call normal_01_cdf ( + b0, c4 )

  ln_mean = exp ( mu + 0.5D+00 * sigma ** 2 )

  mean = ln_mean * ( c1 - c2 ) / ( c4 - c3 )

  return
end
subroutine log_normal_truncated_ab_pdf ( x, mu, sigma, a, b, pdf )

!*****************************************************************************80
!
!! LOG_NORMAL_TRUNCATED_AB_PDF evaluates the truncated AB Lognormal PDF.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 March 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, the argument of the PDF.
!    0.0 < X.
!
!    Input, real ( kind = rk ) MU, SIGMA, the parameters of the log normal PDF.
!    0.0 < SIGMA.
!
!    Input, real ( kind = rk ) A, B, the lower and upper truncation limits.
!    A < B.
!
!    Output, real ( kind = rk ) PDF, the value of the PDF.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  logical ( kind = 4 ) check
  real ( kind = rk ) lncdf_a
  real ( kind = rk ) lncdf_b
  real ( kind = rk ) lnpdf_x
  real ( kind = rk ) mu
  real ( kind = rk ) pdf
  real ( kind = rk ) sigma
  real ( kind = rk ) x

  call log_normal_truncated_ab_check ( mu, sigma, a, b, check )

  if ( .not. check ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'LOG_NORMAL_TRUNCATED_AB_PDF - Fatal error!'
    write ( *, '(a)' ) '  Parameters are not legal.'
  end if

  if ( x <= a ) then

    pdf = 0.0D+00

  else if ( b <= x ) then

    pdf = 0.0D+00

  else

    call log_normal_cdf ( a, mu, sigma, lncdf_a )
    call log_normal_cdf ( b, mu, sigma, lncdf_b )
    call log_normal_pdf ( x, mu, sigma, lnpdf_x )

    pdf = lnpdf_x / ( lncdf_b - lncdf_a )

  end if

  return
end
subroutine log_normal_truncated_ab_sample ( mu, sigma, a, b, seed, x )

!*****************************************************************************80
!
!! LOG_NORMAL_TRUNCATED_AB_SAMPLE samples the truncated AB Lognormal CDF.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 March 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) MU, SIGMA, the parameters of the log normal PDF.
!    0.0 < SIGMA.
!
!    Input, real ( kind = rk ) A, B, the lower and upper truncation limits.
!    A < B.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, real ( kind = rk ) X, a random sample
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) cdf
  logical ( kind = 4 ) check
  real ( kind = rk ) lncdf_a
  real ( kind = rk ) lncdf_b
  real ( kind = rk ) mu
  real ( kind = rk ) r8_uniform_ab
  integer ( kind = 4 ) seed
  real ( kind = rk ) sigma
  real ( kind = rk ) x

  call log_normal_truncated_ab_check ( mu, sigma, a, b, check )

  if ( .not. check ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'LOG_NORMAL_TRUNCATED_AB_SAMPLE - Fatal error!'
    write ( *, '(a)' ) '  Parameters are not legal.'
  end if

  call log_normal_cdf ( a, mu, sigma, lncdf_a )
  call log_normal_cdf ( b, mu, sigma, lncdf_b )

  cdf = r8_uniform_ab ( lncdf_a, lncdf_b, seed )

  call log_normal_cdf_inv ( cdf, mu, sigma, x )

  return
end
subroutine log_normal_truncated_ab_variance ( mu, sigma, a, b, variance )

!*****************************************************************************80
!
!! LOG_NORMAL_TRUNCATED_AB_VARIANCE: variance of Log Normal Truncated AB PDF.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 March 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) MU, SIGMA, the parameters of the log normal PDF.
!    0.0 < SIGMA.
!
!    Input, real ( kind = rk ) A, B, the lower and upper truncation limits.
!    A < B.
!
!    Output, real ( kind = rk ) VARIANCE, the mean of the PDF.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) a0
  real ( kind = rk ) b
  real ( kind = rk ) b0
  real ( kind = rk ) c1
  real ( kind = rk ) c2
  real ( kind = rk ) c3
  real ( kind = rk ) c4
  logical ( kind = 4 ) check
  real ( kind = rk ) ln_xsquared
  real ( kind = rk ) lntab_xsquared
  real ( kind = rk ) mean
  real ( kind = rk ) mu
  real ( kind = rk ) sigma
  real ( kind = rk ) variance

  call log_normal_truncated_ab_check ( mu, sigma, a, b, check )

  if ( .not. check ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'LOG_NORMAL_TRUNCATED_AB_VARIANCE - Fatal error!'
    write ( *, '(a)' ) '  Parameters are not legal.'
  end if

  call log_normal_truncated_ab_mean ( mu, sigma, a, b, mean )

  a0 = ( log ( a ) - mu ) / sigma
  b0 = ( log ( b ) - mu ) / sigma

  call normal_01_cdf ( 2.0 * sigma - a0, c1 )
  call normal_01_cdf ( 2.0 * sigma - b0, c2 )
  call normal_01_cdf ( + a0, c3 )
  call normal_01_cdf ( + b0, c4 )

  ln_xsquared = exp ( 2.0D+00 * mu + 2.0D+00 * sigma ** 2 )

  lntab_xsquared = ln_xsquared * ( c1 - c2 ) / ( c4 - c3 )

  variance = lntab_xsquared - mean ** 2

  return
end
subroutine normal_01_cdf ( x, cdf )

!*****************************************************************************80
!
!! NORMAL_01_CDF evaluates the Normal 01 CDF.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    AG Adams,
!    Algorithm 39,
!    Areas Under the Normal Curve,
!    Computer Journal,
!    Volume 12, pages 197-198, 1969.
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, the argument of the CDF.
!
!    Output, real ( kind = rk ) CDF, the value of the CDF.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), parameter :: a1 = 0.398942280444D+00
  real ( kind = rk ), parameter :: a2 = 0.399903438504D+00
  real ( kind = rk ), parameter :: a3 = 5.75885480458D+00
  real ( kind = rk ), parameter :: a4 = 29.8213557808D+00
  real ( kind = rk ), parameter :: a5 = 2.62433121679D+00
  real ( kind = rk ), parameter :: a6 = 48.6959930692D+00
  real ( kind = rk ), parameter :: a7 = 5.92885724438D+00
  real ( kind = rk ), parameter :: b0 = 0.398942280385D+00
  real ( kind = rk ), parameter :: b1 = 3.8052D-08
  real ( kind = rk ), parameter :: b2 = 1.00000615302D+00
  real ( kind = rk ), parameter :: b3 = 3.98064794D-04
  real ( kind = rk ), parameter :: b4 = 1.98615381364D+00
  real ( kind = rk ), parameter :: b5 = 0.151679116635D+00
  real ( kind = rk ), parameter :: b6 = 5.29330324926D+00
  real ( kind = rk ), parameter :: b7 = 4.8385912808D+00
  real ( kind = rk ), parameter :: b8 = 15.1508972451D+00
  real ( kind = rk ), parameter :: b9 = 0.742380924027D+00
  real ( kind = rk ), parameter :: b10 = 30.789933034D+00
  real ( kind = rk ), parameter :: b11 = 3.99019417011D+00
  real ( kind = rk ) cdf
  real ( kind = rk ) q
  real ( kind = rk ) x
  real ( kind = rk ) y
!
!  |X| <= 1.28.
!
  if ( abs ( x ) <= 1.28D+00 ) then

    y = 0.5D+00 * x * x

    q = 0.5D+00 - abs ( x ) * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5 &
      + a6 / ( y + a7 ) ) ) )
!
!  1.28 < |X| <= 12.7
!
  else if ( abs ( x ) <= 12.7D+00 ) then

    y = 0.5D+00 * x * x

    q = exp ( - y ) * b0 / ( abs ( x ) - b1 &
      + b2 / ( abs ( x ) + b3 &
      + b4 / ( abs ( x ) - b5 &
      + b6 / ( abs ( x ) + b7 &
      - b8 / ( abs ( x ) + b9 &
      + b10 / ( abs ( x ) + b11 ) ) ) ) ) )
!
!  12.7 < |X|
!
  else

    q = 0.0D+00

  end if
!
!  Take account of negative X.
!
  if ( x < 0.0D+00 ) then
    cdf = q
  else
    cdf = 1.0D+00 - q
  end if

  return
end
subroutine normal_01_cdf_inv ( p, x )

!*****************************************************************************80
!
!! NORMAL_01_CDF_INV inverts the standard normal CDF.
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
!    05 June 2007
!
!  Author:
!
!    Original FORTRAN77 version by Michael Wichura.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Michael Wichura,
!    Algorithm AS241:
!    The Percentage Points of the Normal Distribution,
!    Applied Statistics,
!    Volume 37, Number 3, pages 477-484, 1988.
!
!  Parameters:
!
!    Input, real ( kind = rk ) P, the value of the cumulative probability
!    densitity function.  0 < P < 1.  If P is outside this range, an
!    "infinite" value will be returned.
!
!    Output, real ( kind = rk ) X, the normal deviate value
!    with the property that the probability of a standard normal deviate being
!    less than or equal to the value is P.
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
  real ( kind = rk ) r8poly_value_horner
  real ( kind = rk ), parameter :: split1 = 0.425D+00
  real ( kind = rk ), parameter :: split2 = 5.0D+00
  real ( kind = rk ) x

  if ( p <= 0.0D+00 ) then
    x = - huge ( x )
    return
  end if

  if ( 1.0D+00 <= p ) then
    x = huge ( x )
    return
  end if

  q = p - 0.5D+00

  if ( abs ( q ) <= split1 ) then

    r = const1 - q * q
    x = q * r8poly_value_horner ( 7, a, r ) / r8poly_value_horner ( 7, b, r )

  else

    if ( q < 0.0D+00 ) then
      r = p
    else
      r = 1.0D+00 - p
    end if

    if ( r <= 0.0D+00 ) then

      x = huge ( x )

    else

      r = sqrt ( - log ( r ) )

      if ( r <= split2 ) then

        r = r - const2
        x = r8poly_value_horner ( 7, c, r ) / r8poly_value_horner ( 7, d, r )

      else

        r = r - split2
        x = r8poly_value_horner ( 7, e, r ) / r8poly_value_horner ( 7, f, r )

      end if

    end if

    if ( q < 0.0D+00 ) then
      x = -x
    end if

  end if

  return
end
subroutine normal_cdf ( x, mu, sigma, cdf )

!*****************************************************************************80
!
!! NORMAL_CDF evaluates the Normal CDF.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, the argument of the CDF.
!
!    Input, real ( kind = rk ) MU, SIGMA, the parameters of the PDF.
!    0.0 < SIGMA.
!
!    Output, real ( kind = rk ) CDF, the value of the CDF.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) cdf
  real ( kind = rk ) mu
  real ( kind = rk ) sigma
  real ( kind = rk ) x
  real ( kind = rk ) y

  y = ( x - mu ) / sigma

  call normal_01_cdf ( y, cdf )

  return
end
subroutine normal_cdf_inv ( cdf, mu, sigma, x )

!*****************************************************************************80
!
!! NORMAL_CDF_INV inverts the Normal CDF.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!    Input, real ( kind = rk ) MU, SIGMA, the parameters of the PDF.
!    0.0 < SIGMA.
!
!    Output, real ( kind = rk ) X, the corresponding argument.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) cdf
  real ( kind = rk ) mu
  real ( kind = rk ) sigma
  real ( kind = rk ) x
  real ( kind = rk ) x2

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NORMAL_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop 1
  end if

  call normal_01_cdf_inv ( cdf, x2 )

  x = mu + sigma * x2

  return
end
subroutine r8_fake_use ( x )

!*****************************************************************************80
!
!! r8_fake_use() pretends to use an R8 variable.
!
!  Discussion:
!
!    Some compilers will issue a warning if a variable is unused.
!    Sometimes there's a good reason to include a variable in a program,
!    but not to use it.  Calling this function with that variable as
!    the argument will shut the compiler up.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 April 2020
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) X, the variable to be "used".
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) x

  if ( x /= x ) then
    write ( *, '(a)' ) '  r8_fake_use: variable is NAN.'
  end if

  return
end
function r8_uniform_ab ( a, b, seed )

!*****************************************************************************80
!
!! R8_UNIFORM_AB returns a scaled pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = rk ) value.
!
!    For now, the input quantity SEED is an integer variable.
!
!    The pseudorandom number should be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = rk ) R8_UNIFORM_AB, a number strictly between A and B.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  real ( kind = rk ) r8_uniform_ab
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM_AB - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r8_uniform_ab = a + ( b - a ) * real ( seed, kind = rk ) * 4.656612875D-10

  return
end
function r8poly_value_horner ( m, c, x )

!*****************************************************************************80
!
!! R8POLY_VALUE_HORNER evaluates a polynomial using Horner's method.
!
!  Discussion:
!
!    The polynomial 
!
!      p(x) = c0 + c1 * x + c2 * x^2 + ... + cm * x^m
!
!    is to be evaluated at the value X.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 January 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the degree.
!
!    Input, real ( kind = rk ) C(0:M), the polynomial coefficients.  
!    C(I) is the coefficient of X^I.
!
!    Input, real ( kind = rk ) X, the evaluation point.
!
!    Output, real ( kind = rk ) R8POLY_VALUE_HORNER, the polynomial value.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ( kind = 4 ) m

  real ( kind = rk ) c(0:m)
  integer ( kind = 4 ) i
  real ( kind = rk ) r8poly_value_horner
  real ( kind = rk ) value
  real ( kind = rk ) x

  value = c(m)
  do i = m - 1, 0, -1
    value = value * x + c(i)
  end do

  r8poly_value_horner = value

  return
end
subroutine r8vec_max ( n, a, amax )

!*****************************************************************************80
!
!! R8VEC_MAX returns the maximum value in an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = rk ) A(N), the array.
!
!    Output, real ( kind = rk ) AMAX, the value of the largest entry.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ( kind = 4 ) n

  real ( kind = rk ) a(n)
  real ( kind = rk ) amax

  amax = maxval ( a(1:n) )

  return
end
subroutine r8vec_mean ( n, x, mean )

!*****************************************************************************80
!
!! R8VEC_MEAN returns the mean of an R8VEC.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = rk ) X(N), the vector whose mean is desired.
!
!    Output, real ( kind = rk ) MEAN, the mean, or average,
!    of the vector entries.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ( kind = 4 ) n

  real ( kind = rk ) mean
  real ( kind = rk ) x(n)

  mean = sum ( x(1:n) ) / real ( n, kind = rk )

  return
end
subroutine r8vec_min ( n, a, amin )

!*****************************************************************************80
!
!! R8VEC_MIN returns the minimum value of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = rk ) A(N), the array.
!
!    Output, real ( kind = rk ) AMIN, the value of the smallest entry.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ( kind = 4 ) n

  real ( kind = rk ) a(n)
  real ( kind = rk ) amin

  amin = minval ( a(1:n) )

  return
end
subroutine r8vec_variance ( n, x, variance )

!*****************************************************************************80
!
!! R8VEC_VARIANCE returns the variance of an R8VEC.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = rk ) X(N), the vector whose variance is desired.
!
!    Output, real ( kind = rk ) VARIANCE, the variance of the vector entries.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ( kind = 4 ) n

  real ( kind = rk ) mean
  real ( kind = rk ) variance
  real ( kind = rk ) x(n)

  call r8vec_mean ( n, x, mean )

  variance = sum ( ( x(1:n) - mean ) ** 2 )

  if ( 1 < n ) then
    variance = variance / real ( n - 1, kind = rk )
  else
    variance = 0.0D+00
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
!  Parameters:
!
!    None
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
