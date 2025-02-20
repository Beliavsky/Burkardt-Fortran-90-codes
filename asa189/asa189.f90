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
!    This version by John Burkardt.
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
subroutine bbl ( mu, theta, rl, mrl, lm, rnl )

!*****************************************************************************80
!
!! BBL calculates the beta binomial log likelihood.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 March 1999
!
!  Author:
!
!    Original FORTRAN77 version by D Smith.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    D Smith,
!    Algorithm AS 189,
!    Maximum Likelihood Estimation of the Parameters of the Beta
!    Binomial Distribution,
!    Applied Statistics,
!    Volume 32, Number 2, 1983, pages 196-204.
!
!  Parameters:
!
!    Input, real ( kind = rk ) MU, the estimated value of MU.
!
!    Input, real ( kind = rk ) THETA, the estimated value of THETA.
!
!    Input, integer RL(MRL,3), array of coefficients of
!    (MU + R * THETA), ( 1 - MU + R * THETA) and ( 1 + R * THETA ) terms.
!
!    Input, integer MRL, the first dimension of the RL array,
!    which must be at least the maximum of the values in IN(*).
!
!    Input, integer LM(3), contain the values Max ( IX(J) - 1 ),
!    Max ( IN(J) - IX(J) - 1 ), and Max ( IN(J) - 1 ).
!
!    Output, real ( kind = rk ) RNL, the log likelihood.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer mrl

  real ( kind = rk ) a
  integer i
  integer lm(3)
  real ( kind = rk ) mu
  integer mlm
  integer rl(mrl,3)
  real ( kind = rk ) rnl
  real ( kind = rk ) theta

  rnl = 0.0D+00
  mlm = lm(3)

  do i = 1, mlm

    a = real ( i - 1, kind = rk ) * theta

    if ( i <= lm(1) ) then
      rnl = rnl + real ( rl(i,1), kind = rk ) * log ( mu + a )
    end if

    if ( i <= lm(2) ) then
      rnl = rnl + real ( rl(i,2), kind = rk ) * log ( 1.0D+00 - mu + a )
    end if

    rnl = rnl - real ( rl(i,3), kind = rk ) * log ( 1.0D+00 + a )

  end do

  return
end
subroutine bbme ( n, ix, in, inf, mu, theta )

!*****************************************************************************80
!
!! BBME estimates MU and THETA of the beta binomial distribution.
!
!  Discussion:
!
!    The method of moments is used.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 February 2003
!
!  Author:
!
!    Original FORTRAN77 version by D Smith.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    D Smith,
!    Algorithm AS 189,
!    Maximum Likelihood Estimation of the Parameters of the Beta
!    Binomial Distribution,
!    Applied Statistics,
!    Volume 32, Number 2, 1983, pages 196-204.
!
!  Parameters:
!
!    Input, integer N, the number of observations or trials.
!
!    Input, integer IX(N), contains the number of successes for
!    each trial.
!
!    Input, integer IN(N), the number tested on each trial.
!
!    Input, real ( kind = rk ) INF, the maximum acceptable value for THETA.
!
!    Output, real ( kind = rk ) MU, the estimate for MU.
!
!    Output, real ( kind = rk ) THETA, the estimate for THETA.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) d1
  real ( kind = rk ) d2
  integer i
  integer in(n)
  real ( kind = rk ) inf
  integer ix(n)
  logical j
  real ( kind = rk ) mu
  real ( kind = rk ) p(n)
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) theta
  real ( kind = rk ) tp
  real ( kind = rk ) w(n)
  real ( kind = rk ) w_sum

  j = .false.
  w(1:n) = real ( in(1:n), kind = rk )
  p(1:n) = real ( ix(1:n), kind = rk ) / real ( in(1:n), kind = rk )

  do

    w_sum = sum ( w(1:n) )
    tp = dot_product ( w(1:n), p(1:n) )

    tp = tp / w_sum

    s = 0.0D+00
    d1 = 0.0D+00
    d2 = 0.0D+00
    do i = 1, n
      r = p(i) - tp
      s = s + w(i) * r * r
      r = w(i) * ( 1.0D+00 - w(i) / w_sum )
      d1 = d1 + r / real ( in(i), kind = rk )
      d2 = d2 + r
    end do

    s = real ( n - 1, kind = rk ) * s / real ( n, kind = rk )
    r = tp * ( 1.0D+00 - tp )

    if ( r == 0.0D+00 ) then
      exit
    end if

    r = ( s - r * d1 ) / ( r * ( d2 - d1 ) )
    r = max ( r, 0.0D+00 )

    if ( j ) then
      exit
    end if

    w(1:n) = w(1:n) / ( 1.0D+00 + r * ( w(1:n) - 1.0D+00 ) )

    j = .true.

  end do
!
!  Set the estimates.
!
  mu = tp

  if ( r < 1.0D+00 ) then

    theta = r / ( 1.0D+00 - r )

    theta = min ( theta, inf )

  else

    theta = inf

  end if

  return
end
subroutine bbml ( n, ix, in, rl, mrl, iter, ccrit, mu, theta, mu_se, theta_se, &
  rnl, ifault )

!*****************************************************************************80
!
!! BBML estimates the parameters of a beta binomial distribution.
!
!  Discussion:
!
!    The beta binomial probability density function for X successes
!    out of N trials is
!
!      PDF(X) ( N, MU, THETA ) =
!        C(N,X) * Product ( 0 <= R <= X - 1 ) ( MU + R * THETA )
!               * Product ( 0 <= R <= N - X - 1 ) ( 1 - MU + R * THETA )
!               / Product ( 0 <= R <= N - 1 )  ( 1 + R * THETA )
!
!    where
!
!      C(N,X) is the combinatorial coefficient;
!      MU is the expectation of the underlying Beta distribution;
!      THETA is a shape parameter.
!
!    A THETA value of 0 results in a PDF equivalent to the binomial
!    distribution:
!
!      PDF(X) ( N, MU, 0 ) = C(N,X) * MU^X * ( 1 - MU )^(N-X)
!
!    This PDF can be reformulated as:
!
!      PDF2(X)(A,B,C) = Beta(A+X,B+C-X)
!        / ( (C+1) * Beta(X+1,C-X+1) * Beta(A,B) )  for 0 <= X <= C.
!
!    Given A, B, C for PDF2, the equivalent PDF has:
!
!      N     = C
!      MU    = A / ( A + B )
!      THETA = 1 / ( A + B )
!
!    Given N, MU, THETA for PDF, the equivalent PDF2 has:
!
!      A = MU / THETA
!      B = ( 1 - MU ) / THETA
!      C = N
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 August 2012
!
!  Author:
!
!    Original FORTRAN77 version by D Smith.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    D Smith,
!    Algorithm AS 189,
!    Maximum Likelihood Estimation of the Parameters of the Beta
!    Binomial Distribution,
!    Applied Statistics,
!    Volume 32, Number 2, 1983, pages 196-204.
!
!  Parameters:
!
!    Input, integer N, the number of observations or trials.
!
!    Input, integer IX(N), contains the number of successes for
!    each trial.
!
!    Input, integer IN(N), the number tested on each trial.
!
!    Workspace, integer RL(MRL,3), array of coefficients of
!    (MU + R * THETA), ( 1 - MU + R * THETA) and ( 1 + R * THETA ) terms.
!
!    Input, integer MRL, the first dimension of the RL array,
!    which must be at least the maximum of the values in IN(*).
!
!    Input/output, integer ITER;
!    On input, the maximum number of iterations allowed.
!    On output, the number of iterations taken.
!
!    Input, real ( kind = rk ) CCRIT, the convergence criterion.  The iteration 
!    is judged to have converged when abs ( delta MU ) and
!    abs ( delta THETA) are less than or equal to CCRIT.
!
!    Output, real ( kind = rk ) MU, the maximum likelihood estimate of MU, the 
!    mean of the beta binomial distribution.
!
!    Output, real ( kind = rk ) THETA, the maximum likelihood estimate of THETA, 
!    the shape parameter of the beta binomial distribution.
!
!    Output, real ( kind = rk ) MU_SE, the standard error of the estimate of MU;
!    returned as -1.0 if it cannot be calculated.
!
!    Output, real ( kind = rk ) THETA_SE, the standard error of the estimate of
!    THETA; returned as -1.0 if it cannot be calculated.
!
!    Output, real ( kind = rk ) RNL, the log likelihood for the maximum 
!    likelihood estimates.
!
!    Output, integer IFAULT, error flag.
!    0, no error.
!    1, N <= 1;
!    2, IX(I) = 0 for all I;
!    3, IX(I) = IN(I) for all I;
!    4, max ( IN(I) ) > MRL;
!    5, either IX(I) < 0 or IN(I) < IX(I) for some I;
!    6, MU went outside the range of [0,1], or THETA went outside the
!       range [0,INF], where INF represents Infinity;
!    7, if the maximum number of iterations was exceeded;
!    8, if the damped Newton-Raphson iteration failed.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) ccrit
  real ( kind = rk ) d
  real ( kind = rk ) del
  real ( kind = rk ) dum
  real ( kind = rk ) e
  real ( kind = rk ) eps
  real ( kind = rk ) f
  real ( kind = rk ) fd(2)
  integer ifault
  integer in(n)
  real ( kind = rk ), parameter :: inf = 1.0D+06
  integer iter
  integer iter_max
  integer ix(n)
  integer lm(3)
  logical mc
  real ( kind = rk ) mu
  real ( kind = rk ) mu_se
  integer mrl
  integer nnd
  integer rd1(2,2)
  integer rd2(2,3)
  integer rd3(2,4)
  integer rl(mrl,3)
  real ( kind = rk ) rnl
  real ( kind = rk ) sd(3)
  real ( kind = rk ) td(4)
  real ( kind = rk ) theta
  real ( kind = rk ) theta_se
  real ( kind = rk ) ub(2)

  rd1(1,1) =   1
  rd1(2,1) = - 1
  rd1(1,2) =   1
  rd1(2,2) =   1

  rd2(1,1) = - 1
  rd2(2,1) = - 1
  rd2(1,2) = - 1
  rd2(2,2) =   1
  rd2(1,3) = - 1
  rd2(2,3) = - 1

  rd3(1,1) =   2
  rd3(2,1) = - 2
  rd3(1,2) =   2
  rd3(2,2) =   2
  rd3(1,3) =   2
  rd3(2,3) = - 2
  rd3(1,4) =   2
  rd3(2,4) =   2

  iter_max = iter
  iter = 0
  mc = .true.
  ub(1) = 0.01D+00
  ub(2) = 0.01D+00
!
!  Set the arrays RL and LM.
!
  call set ( n, ix, in, rl, mrl, lm, ifault )

  if ( ifault /= 0 ) then
    return
  end if

  mu_se = - 1.0D+00
  theta_se = - 1.0D+00
  nnd = 0
!
!  Calculate initial estimates by the method of moments.
!
  call bbme ( n, ix, in, inf, mu, theta )

  if ( theta == inf ) then
    ifault = 6
    call bbl ( mu, theta, rl, mrl, lm, rnl )
    return
  end if
!
!  Newton-Raphson iteration on first derivatives.
!
  do

    do

      if ( iter_max < iter ) then
        ifault = 7
        call bbl ( mu, theta, rl, mrl, lm, rnl )
        return
      end if
!
!  Calculate the first derivatives of the log likelihood.
!
      call gder ( mu, theta, rl, mrl, lm, 2, rd1, fd )
!
!  Calculate the second derivatives of the log likelihood.
!
      call gder ( mu, theta, rl, mrl, lm, 3, rd2, sd )
!
!  Calculate the third derivatives of the log likelihood.
!
      call gder ( mu, theta, rl, mrl, lm, 4, rd3, td )
!
!  Calculate the increments.
!
      dum = sd(1) * sd(3) - sd(2) * sd(2)

      if ( sd(1) < 0.0D+00 .and. 0.0D+00 < dum ) then
        exit
      end if
!
!  Non negative definite matrix.
!
      nnd = nnd + 1
!
!  SD(1) is always negative so a gradient step is made on MU.
!
      a = mu - fd(1) / sd(1)
      b = theta

      if ( fd(2) /= 0.0D+00 ) then
        b = b + sign ( ub(2), fd(2) )
      end if

      if ( a <= 0.0D+00 ) then
        a = 0.0001D+00
      else if ( 1.0D+00 <= a ) then
        a = 0.9999D+00
      end if

      b = max ( b, 0.0D+00 )
      b = min ( b, inf )

      call bbl ( mu, theta, rl, mrl, lm, c )
      call bbl ( a, b, rl, mrl, lm, d )

      if ( 10 < nnd .or. d <= c ) then
        ifault = 8
        call bbl ( mu, theta, rl, mrl, lm, rnl )
        return
      end if

      iter = iter + 1
      mu = a
      theta = b

    end do

    del = ( fd(2) * sd(2) - fd(1) * sd(3) ) / dum
    eps = ( fd(1) * sd(2) - fd(2) * sd(1) ) / dum
!
!  Check to see if the Lipschitz condition is satisfied.
!
    a = sd(2) * td(2) - td(1) * sd(3)
    b = sd(2) * td(3) - td(2) * sd(3)
    c = td(1) * sd(2) - td(2) * sd(1)
    d = sd(2) * td(2) - sd(1) * td(3)
    e = sd(2) * td(4) - td(3) * sd(3)
    f = td(3) * sd(2) - td(4) * sd(1)

    a = del * a + eps * b
    c = del * c + eps * d
    e = del * b + eps * e
    f = del * d + eps * f

    dum = ( a * a + c * c + e * e + f * f ) / dum / dum
!
!  Failure of the Lipschitz condition.
!  A step in the direction of the gradient is made.
!
    if ( 1.0D+00 <= dum ) then

      a = fd(1)**2
      b = fd(2)**2
      c = a * sd(1) + 2.0D+00 * sd(2) * fd(1) * fd(2) + b * sd(3)

      if ( c /= 0.0D+00 ) then

        c = - ( a + b ) / c
        del = c * fd(1)
        eps = c * fd(2)

        if ( ub(1) < abs ( del ) ) then
          del = sign ( ub(1), del )
        end if

        ub(1) = 2.0D+00 * abs ( del )
 
        if ( ub(2) < abs ( eps ) ) then
          eps = sign ( ub(2), eps )
        end if

        ub(2) = 2.0D+00 * abs ( eps )

      else

        if ( fd(1) /= 0.0D+00 ) then
          del = sign ( ub(1), fd(1) )
        else
          del = 0.0D+00
        end if

        if ( fd(2) /= 0.0D+00 ) then
          eps = sign ( ub(2), fd(2) )
        else
          eps = 0.0D+00
        end if

      end if

      call bbl ( mu, theta, rl, mrl, lm, c )
!
!  Begin loop.
!
      do

        a = mu + del
        b = theta + eps
  
        if ( a <= 0.0D+00 ) then
          a = 0.0001D+00
        end if

        if ( 1.0D+00 <= a ) then
          a = 0.9999D+00
        end if

        del = a - mu

        b = max ( b, 0.0D+00 )
        b = min ( b, inf )

        eps = b - theta
        call bbl ( a, b, rl, mrl, lm, d )
!
!  Check to see if gradient step has increased log likelihood.
!
        if ( c < d ) then
          exit
        end if

        del = del / 2.0D+00
        eps = eps / 2.0D+00

        if ( abs ( del ) <= ccrit .and. abs ( eps ) <= ccrit ) then
          ifault = 8
          call bbl ( mu, theta, rl, mrl, lm, rnl )
          return
        end if

      end do

    else if ( abs ( del ) <= ccrit .and. abs ( eps ) <= ccrit ) then

      mc = .false.

    end if

    iter = iter + 1
    a = mu + del
    b = theta + eps

    if ( a <= 0.0D+00 .or. &
         1.0D+00 <= a .or. &
         b < 0.0D+00 .or. &
         inf < b ) then

      if ( a <= 0.0D+00 ) then
        mu = 0.0D+00
      end if

      if ( 1.0D+00 <= a ) then
        mu = 1.0D+00
      end if

      if ( b < 0.0D+00 ) then
        theta = 0.0D+00
      end if

      if ( inf < b ) then
        theta = inf
      end if

      ifault = 6
      call bbl ( mu, theta, rl, mrl, lm, rnl )
      return

    end if

    mu = a
    theta = b

    if ( .not. mc ) then
      exit
    end if

  end do
!
!  Calculate the standard errors.
!
  if ( sd(1) < 0.0D+00 ) then
    mu_se = sqrt ( - 1.0D+00 / sd(1) )
  end if

  if ( sd(3) < 0.0D+00 ) then
    theta_se = sqrt ( - 1.0D+00 / sd(3) )
  end if
!
!  Calculate the log likelihood.
!
  call bbl ( mu, theta, rl, mrl, lm, rnl )

  return
end
function beta ( a, b )

!*****************************************************************************80
!
!! BETA returns the value of the Beta function.
!
!  Discussion:
!
!    BETA(A,B) = ( GAMMA ( A ) * GAMMA ( B ) ) / GAMMA ( A + B )
!              = Integral ( 0 <= T <= 1 ) T^(A-1) (1-T)^(B-1) dT.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, the parameters of the function.
!    0.0 < A,
!    0.0 < B.
!
!    Output, real ( kind = rk ) BETA, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) alngam
  real ( kind = rk ) b
  real ( kind = rk ) beta
  integer ifault
!
!  Check.
!
  if ( a <= 0.0D+00 .or. b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BETA - Fatal error!'
    write ( *, '(a)' ) '  Both A and B must be greater than 0.'
    stop
  end if

  beta = exp ( alngam ( a, ifault ) &
             + alngam ( b, ifault ) &
             - alngam ( a + b, ifault ) )

  return
end
subroutine beta_binomial_cdf_inv ( cdf, a, b, c, x )

!*****************************************************************************80
!
!! BETA_BINOMIAL_CDF_INV inverts the Beta Binomial CDF.
!
!  Discussion:
!
!    A simple-minded discrete approach is used.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) CDF, the value of the CDF.
!
!    Input, real ( kind = rk ) A, B, parameters of the PDF.
!    0.0 < A,
!    0.0 < B.
!
!    Input, integer C, a parameter of the PDF.
!    0 <= C.
!
!    Output, integer X, the smallest X whose cumulative density
!    function is greater than or equal to CDF.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) beta
  integer c
  real ( kind = rk ) cdf
  real ( kind = rk ) cum
  real ( kind = rk ) pdf
  integer x
  integer y

  if ( cdf <= 0.0D+00 ) then

    x = 0

  else

    cum = 0.0D+00

    do y = 0, c

      pdf = beta ( a + real ( y, kind = rk ), &
        b + real ( c - y, kind = rk ) ) / ( real ( c + 1, kind = rk ) &
        * beta ( real ( y + 1, kind = rk ), &
        real ( c - y + 1, kind = rk ) ) * beta ( a, b ) )

      cum = cum + pdf

      if ( cdf <= cum ) then
        x = y
        return
      end if

    end do

    x = c

  end if

  return
end
subroutine beta_binomial_check ( a, b, c )

!*****************************************************************************80
!
!! BETA_BINOMIAL_CHECK() checks the parameters of the Beta Binomial PDF.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, parameters of the PDF.
!    0.0 < A,
!    0.0 < B.
!
!    Input, integer C, a parameter of the PDF.
!    0 <= C.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  integer c

  if ( a <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BETA_BINOMIAL_CHECK - Fatal error!'
    write ( *, '(a)' ) '  A <= 0.'
    stop
  end if

  if ( b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BETA_BINOMIAL_CHECK - Fatal error!'
    write ( *, '(a)' ) '  B <= 0.'
    stop
  end if

  if ( c < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BETA_BINOMIAL_CHECK - Fatal error!'
    write ( *, '(a)' ) '  C < 0.'
    stop
  end if

  return
end
subroutine beta_binomial_sample ( a, b, c, x )

!*****************************************************************************80
!
!! beta_binomial_sample() samples the Beta Binomial CDF.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) A, B, parameters of the PDF.
!    0.0 < A,
!    0.0 < B.
!
!    integer C, a parameter of the PDF.
!    0 <= C.
!
!  Output:
!
!    integer( kind = 4 )  X, a sample of the PDF.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  integer c
  real ( kind = rk ) cdf
  integer x

  call random_number ( harvest = cdf )

  call beta_binomial_cdf_inv ( cdf, a, b, c, x )

  return
end
subroutine gder ( mu, theta, rl, mrl, lm, ider, rd, pd )

!*****************************************************************************80
!
!! GDER computes derivatives of the log likelihood.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 March 1999
!
!  Author:
!
!    Original FORTRAN77 version by D Smith.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    D Smith,
!    Algorithm AS 189,
!    Maximum Likelihood Estimation of the Parameters of the Beta
!    Binomial Distribution,
!    Applied Statistics,
!    Volume 32, Number 2, 1983, pages 196-204.
!
!  Parameters:
!
!    Input, real ( kind = rk ) MU, the estimated value of MU.
!
!    Input, real ( kind = rk ) THETA, the estimated value of THETA.
!
!    Input, integer RL(MRL,3), array of coefficients of
!    (MU + R * THETA), ( 1 - MU + R * THETA) and ( 1 + R * THETA ) terms.
!
!    Input, integer MRL, the first dimension of the RL array,
!    which must be at least the maximum of the values in IN(*).
!
!    Input, integer LM(3), contain the values Max ( IX(J) - 1 ),
!    Max ( IN(J) - IX(J) - 1 ), and Max ( IN(J) - 1 ).
!
!    Input, integer IDER, 1 plus the order of the derivative
!    desired.  IDER can be 2, 3 or 4.
!
!    Input, integer RD(2,IDER), an array of coefficients.
!
!    Output, real ( kind = rk ) PD(IDER), the derivatives of the log likelihood.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ider
  integer mrl

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) d
  integer i
  integer j
  integer k
  integer kk
  integer lm(3)
  real ( kind = rk ) mu
  integer mlm
  real ( kind = rk ) pd(ider)
  integer rd(2,ider)
  integer rl(mrl,3)
  real ( kind = rk ) theta

  mlm = lm(3)
  kk = ider - 1

  pd(1:ider) = 0.0D+00

  do i = 1, mlm

    c = real ( i - 1, kind = rk )
    a = c * theta

    do j = 1, 3

      if ( i <= lm(j) ) then

        if ( j == 1 ) then
          d = mu + a
        else if ( j == 2 ) then
          d = 1.0D+00 - mu + a
        else if ( j == 3 ) then
          d = 1.0D+00 + a
        end if

        b = real ( rl(i,j) ) / d**kk

        if ( j /= 3 ) then

          do k = 1, ider
            pd(k) = pd(k) + real ( rd(j,k), kind = rk ) * b
            b = b * c
          end do

        else

          d = - real ( rd(1,1),  kind = rk ) * b * c**kk
          pd(ider) = pd(ider) + d

        end if

      end if

    end do

  end do

  return
end
subroutine set ( n, ix, in, rl, mrl, lm, ifault )

!*****************************************************************************80
!
!! set() sets up the arrays RL and LM.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 March 1999
!
!  Author:
!
!    Original FORTRAN77 version by D Smith.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    D Smith,
!    Algorithm AS 189,
!    Maximum Likelihood Estimation of the Parameters of the Beta
!    Binomial Distribution,
!    Applied Statistics,
!    Volume 32, Number 2, 1983, pages 196-204.
!
!  Parameters:
!
!    Input, integer N, the number of observations or trials.
!
!    Input, integer IX(N), contains the number of successes for
!    each trial.
!
!    Input, integer IN(N), the number tested on each trial.
!
!    Output, integer RL(MRL,3), array of coefficients of
!    (MU + R * THETA), ( 1 - MU + R * THETA) and ( 1 + R * THETA ) terms.
!
!    Input, integer MRL, the first dimension of the RL array,
!    which must be at least the maximum of the values in IN(*).
!
!    Output, integer LM(3), contain the values Max ( IX(J) - 1 ),
!    Max ( IN(J) - IX(J) - 1 ), and Max ( IN(J) - 1 ).
!
!    Output, integer IFAULT, error flag.
!    0, no error.
!    1, N <= 1;
!    2, IX(I) = 0 for all I;
!    3, IX(I) = IN(I) for all I;
!    4, max ( IN(I) ) > MRL;
!    5, either IX(I) < 0 or IN(I) < IX(I) for some I.
!
  implicit none

  integer mrl
  integer n

  integer i
  integer ifault
  integer in(n)
  integer ix(n)
  integer j
  integer jj
  integer k
  integer lm(3)
  integer mar
  integer rl(mrl,3)
!
!  Check the input data.
!
  if ( n <= 1 ) then
    ifault = 1
    return
  end if

  ifault = 2
  do i = 1, n
    if ( 0 < ix(i) ) then
      ifault = 0
    end if
  end do

  if ( ifault == 2 ) then
    return
  end if

  ifault = 3
  do i = 1, n
    if ( ix(i) < in(i) ) then
      ifault = 0
    end if
  end do

  if ( ifault == 3 ) then
    return
  end if
!
!  Form the matrix of counts.
!
  ifault = 4

  lm(1:3) = 0
  rl(1:mrl,1:3) = 0

  do i = 1, n

    jj = ix(i)
    mar = 1

    do

      if ( jj < 0 ) then
        ifault = 5
        return
      end if

      if ( 0 < jj ) then

        if ( mrl < jj ) then
          return
        end if

        if ( lm(mar) < jj ) then
          lm(mar) = jj
        end if

        rl(jj,mar) = rl(jj,mar) + 1

      end if

      if ( mar == 1 ) then
        jj = in(i) - ix(i)
        mar = 2
      else if ( mar == 2 ) then
        jj = in(i)
        mar = 3
      else
        exit
      end if

    end do

  end do

  ifault = 0
!
!  Evaluate number of calls to different terms of likelihood function.
!
  do i = 1, 3

    jj = lm(i) - 1
    k = jj
    do j = 1, jj
      rl(k,i) = rl(k,i) + rl(k+1,i)
      k = k - 1
    end do

  end do

  return
end

