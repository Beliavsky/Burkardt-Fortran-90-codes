program main

!*****************************************************************************80
!
!! MAIN is the main program for QUADPACK_DOUBLE_TEST.
!
!  Discussion:
!
!    QUADPACK_DOUBLE_TEST tests the QUADPACK_DOUBLE library.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QUADPACK_DOUBLE_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the QUADPACK_DOUBLE library.'

  call dqag_test ( )
  call dqagi_test ( )
  call dqagp_test ( )
  call dqags_test ( )
  call dqawc_test ( )
  call dqawf_test ( )
  call dqawo_test ( )
  call dqaws_test ( )
  call dqk15_test ( )
  call dqk21_test ( )
  call dqk31_test ( )
  call dqk41_test ( )
  call dqk51_test ( )
  call dqk61_test ( )
  call dqng_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QUADPACK_DOUBLE_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )
 
  stop 0
end
subroutine dqag_test ( )

!*****************************************************************************80
!
!! DQAG_TEST tests DQAG.
!
!  Discussion:
!
!    DQAG is an adaptive automatic integrator using a Gauss-Kronrod rule.
!
!    Integrate cos(100*sin(x)) from 0 to pi.
!
!    The exact answer is pi * j0(100), or roughly 0.06278740.
!
!    KEY chooses the order of the integration rule, from 1 to 6.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: limit = 500

  integer ( kind = 4 ), parameter :: lenw = 4 * limit

  real ( kind = 8 ), parameter :: a = 0.0D+00
  real ( kind = 8 ) abserr
  real ( kind = 8 ) b
  real ( kind = 8 ), parameter :: epsabs = 0.0D+00
  real ( kind = 8 ), parameter :: epsrel = 0.001D+00
  real ( kind = 8 ), external :: f02
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) iwork(limit)
  integer ( kind = 4 ), parameter :: key = 6
  integer ( kind = 4 ) last
  integer ( kind = 4 ) neval
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = 8 ) result
  real ( kind = 8 ), parameter :: true = 0.06278740D+00
  real ( kind = 8 ) work(lenw)

  b = r8_pi

  call dqag ( f02, a, b, epsabs, epsrel, key, result, abserr, neval, ier, &
    limit, lenw, last, iwork, work )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DQAG_TEST'
  write ( *, '(a)' ) '  Test DQAG'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is COS(100*SIN(X))'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,i8)' ) '  Number of function evaluations, NEVAL = ', neval
  write ( *, '(a,i8)' ) '  Error return code IER = ', ier

  return
end
subroutine dqagi_test ( )

!*****************************************************************************80
!
!! DQAGI_TEST tests DQAGI.
!
!  Discussion:
!
!    DQAGI is an adaptive quadrature routine for infinite intervals.
!
!    Integrate log(x)/(1+100*x*x) from 0 to infinity.
!
!    The exact answer is -pi*log(10)/20 = -0.3616892
!
!    give the type of infinity
!
!    inf=1 means a to infinity
!       -1      -infinity to a
!        2      -infinity to infinity
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: limit = 500

  integer ( kind = 4 ), parameter :: lenw = limit * 4

  real ( kind = 8 ), parameter :: a = 0.0D+00
  real ( kind = 8 ) abserr
  real ( kind = 8 ), parameter :: epsabs = 0.0D+00
  real ( kind = 8 ), parameter :: epsrel = 0.001D+00
  real ( kind = 8 ), external :: f05
  integer ( kind = 4 ) ier
  integer ( kind = 4 ), parameter :: inf = 1
  integer ( kind = 4 ) iwork(limit)
  integer ( kind = 4 ) last
  integer ( kind = 4 ) neval
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = 8 ) result
  real ( kind = 8 ) true
  real ( kind = 8 ) work(lenw)

  call dqagi ( f05, a, inf, epsabs, epsrel, result, abserr, neval, &
    ier, limit, lenw, last, iwork, work )

  true = - r8_pi * log ( 10.0D+00 ) / 20.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DQAGI_TEST'
  write ( *, '(a)' ) '  Test DQAGI'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is log(x)/(1+100*x*x)'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =    Infinity'
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,i8)' ) '  Number of function evaluations, NEVAL = ', neval
  write ( *, '(a,i8)' ) '  Error return code IER = ', ier

  return
end
subroutine dqagp_test ( )

!*****************************************************************************80
!
!! DQAGP_TEST tests DQAGP.
!
!  Discussion:
!
!    DQAGP is an adaptive integrator that can handle singularities
!    of the integrand at user specified points,
!
!    Integrate 
!
!      x**3 * log(abs( (x*x-1)*(x*x-2) )) 
!
!    from 0 to 3.
!
!    The exact answer is 61*log(2)+77*log(7)/4 - 27.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: limit = 500
  integer ( kind = 4 ), parameter :: npts = 2
  integer ( kind = 4 ), parameter :: npts2 = 2 * npts

  integer ( kind = 4 ), parameter :: leniw = 2 * limit + npts2
  integer ( kind = 4 ), parameter :: lenw = leniw * 2 - npts2

  real ( kind = 8 ), parameter :: a = 0.0D+00
  real ( kind = 8 ) abserr
  real ( kind = 8 ), parameter :: b = 3.0D+00
  real ( kind = 8 ), parameter :: epsabs = 0.0D+00
  real ( kind = 8 ), parameter :: epsrel = 0.001D+00
  real ( kind = 8 ), external :: f04
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) iwork(leniw)
  integer ( kind = 4 ) last
  integer ( kind = 4 ) neval
  real ( kind = 8 ) points(npts2)
  real ( kind = 8 ) result
  real ( kind = 8 ) true
  real ( kind = 8 ) work(lenw)
!
!  Singularity points:
!
  points(1) = 1.0D+00
  points(2) = sqrt ( 2.0D+00 )

  call dqagp ( f04, a, b, npts2, points, epsabs, epsrel, result, abserr, &
    neval, ier, leniw, lenw, last, iwork, work )

  true = 61.0D+00 * log ( 2.0D+00 ) &
    + 77.0D+00 * log ( 7.0D+00 ) / 4.0D+00 - 27.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QAGP_TEST'
  write ( *, '(a)' ) '  Test QAGP'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is x**3 * log(abs((x*x-1)*(x*x-2)))'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,i8)' ) '  Number of function evaluations, NEVAL = ', neval
  write ( *, '(a,i8)' ) '  Error return code IER = ', ier

  return
end
subroutine dqags_test ( )

!*****************************************************************************80
!
!! DQAGS_TEST tests DQAGS.
!
!  Discussion:
!
!    DQAGS is an adaptive integrator for endpoint singularities.
!
!    integrate log(x)/sqrt(x) from 0 to 1.
!
!    The exact answer is -4.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: limit = 500

  integer ( kind = 4 ), parameter :: lenw = limit * 4

  real ( kind = 8 ), parameter :: a = 0.0D+00
  real ( kind = 8 ) abserr
  real ( kind = 8 ), parameter :: b = 1.0D+00
  real ( kind = 8 ), parameter :: epsabs = 0.0D+00
  real ( kind = 8 ), parameter :: epsrel = 0.001D+00
  real ( kind = 8 ), external :: f03
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) iwork(limit)
  integer ( kind = 4 ) last
  integer ( kind = 4 ) neval
  real ( kind = 8 ) result
  real ( kind = 8 ), parameter :: true = -4.0D+00
  real ( kind = 8 ) work(lenw)

  call dqags ( f03, a, b, epsabs, epsrel, result, abserr, neval, ier, &
    limit, lenw, last, iwork, work )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DQAGS_TEST'
  write ( *, '(a)' ) '  Test DQAGS'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is LOG(X)/SQRT(X)'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,i8)' ) '  Number of function evaluations, NEVAL = ', neval
  write ( *, '(a,i8)' ) '  Error return code IER = ', ier

  return
end
subroutine dqawc_test ( )

!*****************************************************************************80
!
!! DQAWC_TEST tests DQAWC.
!
!  Discussion:
!
!    DQAWC is an adaptive integrator for finding the Cauchy
!    principal value of the integral of f(x)*w(x) over (a,b)
!    where w(x)=1/(x-c), c between a and b.
!
!    Integrate 1/(x*(5*x*x*x+6)) from -1 to 5
!
!    The exact answer is log(125/631) / 18 = -0.08994401
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: limit = 500

  integer ( kind = 4 ), parameter :: lenw = limit * 4

  real ( kind = 8 ), parameter :: a = -1.0D+00
  real ( kind = 8 ) abserr
  real ( kind = 8 ), parameter :: b = 5.0D+00
  real ( kind = 8 ), parameter :: c = 0.0D+00
  real ( kind = 8 ), parameter :: epsabs = 0.0D+00
  real ( kind = 8 ), parameter :: epsrel = 0.001D+00
  real ( kind = 8 ), external :: f09
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) iwork(limit)
  integer ( kind = 4 ) last
  integer ( kind = 4 ) neval
  real ( kind = 8 ) result
  real ( kind = 8 ) true
  real ( kind = 8 ) work(lenw)

  call dqawc ( f09, a, b, c, epsabs, epsrel, result, abserr, neval, ier, &
    limit, lenw, last, iwork, work )

  true = log ( 125.0D+00 / 631.0D+00 ) / 18.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DQAWC_TEST'
  write ( *, '(a)' ) '  Test DQAWC'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is 1/(x*(5*x**3+6)'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Point of singularity c =      ', c
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,i8)' ) '  Number of function evaluations, NEVAL = ', neval
  write ( *, '(a,i8)' ) '  Error return code IER = ', ier

  return
end
subroutine dqawf_test ( )

!*****************************************************************************80
!
!! DQAWF_TEST tests DQAWF.
!
!  Discussion:
!
!    DQAWF handles fourier integration of f(x)*w(x) from
!    a to infinity, with w(x)=cos(omega*x) or sine(omega*x)
!
!    integrate cos(pi*x/2) /sqrt(x) from 0 to infinity.
!
!    The exact answer is 1.0
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: limit = 500
  integer ( kind = 4 ), parameter :: limlst = 50

  integer ( kind = 4 ), parameter :: leniw = 2 * limit + limlst
  integer ( kind = 4 ), parameter :: maxp1 = 21
  integer ( kind = 4 ), parameter :: lenw = leniw * 2 + maxp1 * 25

  real ( kind = 8 ), parameter :: a = 0.0D+00
  real ( kind = 8 ) abserr
  real ( kind = 8 ), parameter :: epsabs = 0.001D+00
  real ( kind = 8 ), external :: f07
  integer ( kind = 4 ) ier
  integer ( kind = 4 ), parameter :: integr = 1
  integer ( kind = 4 ) iwork(leniw)
  integer ( kind = 4 ) lst
  integer ( kind = 4 ) neval
  real ( kind = 8 ) omega
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = 8 ) result
  real ( kind = 8 ), parameter :: true = 1.0D+00
  real ( kind = 8 ) work(lenw)
!
!  set argument of sine or cosine
!  set integr=1 for cosine, 2 for sine
!
  omega = 0.5D+00 * r8_pi

  call dqawf ( f07, a, omega, integr, epsabs, result, abserr, neval, ier, &
    limlst, lst, leniw, maxp1, lenw, iwork, work )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DQAWF_TEST'
  write ( *, '(a)' ) '  Test DQAWF'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is cos(pi*x/2)/sqrt(x)'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,i8)' ) '  Number of function evaluations, NEVAL = ', neval
  write ( *, '(a,i8)' ) '  Error return code IER = ', ier

  return
end
subroutine dqawo_test ( )

!*****************************************************************************80
!
!! DQAWO_TEST tests DQAWO.
!
!  Discussion:
!
!    DQAWO integrates functions of the form 
!      f(x) * sin(omega*x)
!    or 
!      f(x) * cos(omega*x)
!
!    Here, we estimate 
!
!      Integral ( 0 <= x <= 1 ) log(x) sin(10*pi*x) dx
!
!    The exact answer is
!
!      exact = - ( gamma + log(10*pi) - ci(10*pi) ) / (10*pi)
!            = - 0.1281316...
!
!    Here Gamma is Euler's constant.
!
!    ci is the cosine integral:
!
!      ci(x) = integral ( x <= v < +oo ) - cos ( v ) / v dv.
!
!    We specify 
!      * INTEGR=1 for integrands with a cosine factor;
!      * INTEGR=2 for integrands with a sine factor.
!
!    Thanks to William Gandler for pointing out errors in the documentation
!    and text of this example, 29 October 2010.
!
!  Modified:
!
!    12 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: leniw = 500
  integer ( kind = 4 ), parameter :: maxp1 = 20
  integer ( kind = 4 ), parameter :: lenw = leniw * 2 + maxp1 * 25

  real ( kind = 8 ), parameter :: a = 0.0D+00
  real ( kind = 8 ) abserr
  real ( kind = 8 ), parameter :: b = 1.0D+00
  real ( kind = 8 ), parameter :: ci = -0.001007D+00
  real ( kind = 8 ), parameter :: epsabs = 0.0D+00
  real ( kind = 8 ), parameter :: epsrel = 0.001D+00
  real ( kind = 8 ), external :: f06
  real ( kind = 8 ), parameter :: gamma = 0.5772156649D+00
  integer ( kind = 4 ) ier
  integer ( kind = 4 ), parameter :: integr = 2
  integer ( kind = 4 ) iwork(leniw)
  integer ( kind = 4 ) last
  integer ( kind = 4 ) neval
  real ( kind = 8 ) omega
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = 8 ) result
  real ( kind = 8 ) true
  real ( kind = 8 ) work(lenw)

  omega = 10.0D+00 * r8_pi

  call dqawo ( f06, a, b, omega, integr, epsabs, epsrel, result, abserr, &
    neval, ier, leniw, maxp1, lenw, last, iwork, work )

  true = - ( gamma + log ( 10.0D+00 * r8_pi ) - ci ) / ( 10.0D+00 * r8_pi )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DQAWO_TEST'
  write ( *, '(a)' ) '  Test DQAWO'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is log(x)*sin(10*pi*x)'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,i8)' ) '  Number of function evaluations, NEVAL = ', neval
  write ( *, '(a,i8)' ) '  Error return code IER = ', ier

  return
end
subroutine dqaws_test ( )

!*****************************************************************************80
!
!! DQAWS_TEST tests DQAWS.
!
!  Discussion:
!
!    DQAWS is an adaptive integrator for integrands with
!    algebraic or logarithmic singularities at the endpoints.
!
!    Here we estimate:
!
!      Integral ( 0 <= x <= 1 ) log(x) / ( 1 + log(x)^2 )^2 dx
!
!    The exact answer is 
!
!      exact = 0.5 * ( ci(1.0) * sin(1.0) - si(G&R)(1.0) * cos(1.0) - 1.0 )
!
!    Numerically:
!
!      exact = -0.18927518788209332118...
!
!    ci is the cosine integral:
!
!      ci(x) = integral ( x <= v < oo ) - cos(v) / v dv
!
!    si (according to Mathematica, for instance) is the sine integral:
!
!      si(x) = integral ( 0 <= v <= x )    sin(v) / v dv
!
!    Note that when Gradshteyn and Rizhik refer to si(x), they
!    mean, instead the complementary form:
!
!      si(G&R)(x) = integral ( x <= v < oo ) sin(v) / v dv
!                 = ( pi / 2 ) - si(x)
!
!    ci(1.0)      = 0.33740392290096813466
!    si(1.0)      = 0.94608307036718301494
!    si(G&R)(1.0) = 0.62471325642771360429
!
!    Thanks to William Gandler for questioning a previous version of
!    the documentation and text of this example, which led to the 
!    clarification of the difference between the Gradshteyn & Rizhik
!    convention versus the Mathematica convention for the sine integral
!    function, 02 November 2010.
!
!    Note that the original QUADPACK documentation lists the answer as
!      (Ci(1)sin(1)+(pi/2-Si(1))*cos(1))/pi
!    which has an incorrect final divisor of pi (it should be 2) and
!    which uses the Mathematica convention for Si.
!
!  Modified:
!
!    12 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: limit = 500

  integer ( kind = 4 ), parameter :: lenw = limit * 4

  real ( kind = 8 ), parameter :: a = 0.0D+00
  real ( kind = 8 ) abserr
  real ( kind = 8 ), parameter :: alfa = 0.0D+00
  real ( kind = 8 ), parameter :: b = 1.0D+00
  real ( kind = 8 ), parameter :: beta = 0.0D+00
  real ( kind = 8 ), parameter :: ci = 0.33740392290096813466D+00
  real ( kind = 8 ), parameter :: epsabs = 0.0D+00
  real ( kind = 8 ), parameter :: epsrel = 0.001D+00
  real ( kind = 8 ), external :: f08
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) integr
  integer ( kind = 4 ) iwork(limit)
  integer ( kind = 4 ) last
  integer ( kind = 4 ) neval
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = 8 ) result
  real ( kind = 8 ), parameter :: si = 0.94608307036718301494D+00
  real ( kind = 8 ) true
  real ( kind = 8 ) work(lenw)
!
!  INTEGR = 2 means the weight function is:
!
!    (x-a)**alfa * (b-x)**beta * log ( x - a )
!
  integr = 2

  call dqaws ( f08, a, b, alfa, beta, integr, epsabs, epsrel, result, &
    abserr, neval, ier, limit, lenw, last, iwork, work )

  true = 0.5D+00 * ( ci * sin ( 1.0D+00 ) &
    + ( r8_pi/ 2.0D+00 - si ) * cos ( 1.0D+00 ) - 1.0D+00 ) 

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DQAWS_TEST'
  write ( *, '(a)' ) '  Test DQAWS'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is log(x)/(1+(log(x))**2)**2'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,i8)' ) '  Number of function evaluations, NEVAL = ', neval
  write ( *, '(a,i8)' ) '  Error return code IER = ', ier

  return
end
subroutine dqk15_test ( )

!*****************************************************************************80
!
!! DQK15_TEST tests DQK15.
!
!  Discussion:
!
!    DQK15 is a Gauss-Kronrod quadrature rule.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), parameter :: a = 0.0D+00
  real ( kind = 8 ) abserr
  real ( kind = 8 ), parameter :: b = 1.0D+00
  real ( kind = 8 ), external :: f01
  real ( kind = 8 ) resabs
  real ( kind = 8 ) resasc
  real ( kind = 8 ) result
  real ( kind = 8 ), parameter :: true = - 4.0D+00 / 9.0D+00

  call dqk15 ( f01, a, b, result, abserr, resabs, resasc )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DQK15_TEST'
  write ( *, '(a)' ) '  Test DQK15'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is SQRT(X)*LOG(X)'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,g14.6)' ) '  RESABS =                      ', resabs
  write ( *, '(a,g14.6)' ) '  RESASC =                      ', resasc

  return
end
subroutine dqk21_test ( )

!*****************************************************************************80
!
!! DQK21_TEST tests DQK21.
!
!  Discussion:
!
!    DQK21 is a Gauss-Kronrod quadrature rule.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), parameter :: a = 0.0D+00
  real ( kind = 8 ) abserr
  real ( kind = 8 ), parameter :: b = 1.0D+00
  real ( kind = 8 ), external :: f01
  real ( kind = 8 ) resabs
  real ( kind = 8 ) resasc
  real ( kind = 8 ) result
  real ( kind = 8 ), parameter :: true = - 4.0D+00 / 9.0D+00

  call dqk21 ( f01, a, b, result, abserr, resabs, resasc )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DQK21_TEST'
  write ( *, '(a)' ) '  Test DQK21'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is SQRT(X)*LOG(X)'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,g14.6)' ) '  RESABS =                      ', resabs
  write ( *, '(a,g14.6)' ) '  RESASC =                      ', resasc

  return
end
subroutine dqk31_test ( )

!*****************************************************************************80
!
!! DQK31_TEST tests DQK31.
!
!  Discussion:
!
!    QK31 is a Gauss-Kronrod quadrature rule.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), parameter :: a = 0.0D+00
  real ( kind = 8 ) abserr
  real ( kind = 8 ), parameter :: b = 1.0D+00
  real ( kind = 8 ), external :: f01
  real ( kind = 8 ) resabs
  real ( kind = 8 ) resasc
  real ( kind = 8 ) result
  real ( kind = 8 ), parameter :: true = - 4.0D+00 / 9.0D+00

  call dqk31 ( f01, a, b, result, abserr, resabs, resasc )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DQK31_TEST'
  write ( *, '(a)' ) '  Test DQK31'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is SQRT(X)*LOG(X)'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,g14.6)' ) '  RESABS =                      ', resabs
  write ( *, '(a,g14.6)' ) '  RESASC =                      ', resasc

  return
end
subroutine dqk41_test ( )

!*****************************************************************************80
!
!! DQK41_TEST tests DQK41.
!
!  Discussion:
!
!    DQK41 is a Gauss-Kronrod quadrature rule.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), parameter :: a = 0.0D+00
  real ( kind = 8 ) abserr
  real ( kind = 8 ), parameter :: b = 1.0D+00
  real ( kind = 8 ), external :: f01
  real ( kind = 8 ) resabs
  real ( kind = 8 ) resasc
  real ( kind = 8 ) result
  real ( kind = 8 ), parameter :: true = - 4.0D+00 / 9.0D+00

  call dqk41 ( f01, a, b, result, abserr, resabs, resasc )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DQK41_TEST'
  write ( *, '(a)' ) '  Test DQK41'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is SQRT(X)*LOG(X)'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,g14.6)' ) '  RESABS =                      ', resabs
  write ( *, '(a,g14.6)' ) '  RESASC =                      ', resasc

  return
end
subroutine dqk51_test ( )

!*****************************************************************************80
!
!! DQK51_TEST tests DQK51.
!
!  Discussion:
!
!    DQK51 is a Gauss-Kronrod quadrature rule.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), parameter :: a = 0.0D+00
  real ( kind = 8 ) abserr
  real ( kind = 8 ), parameter :: b = 1.0D+00
  real ( kind = 8 ), external :: f01
  real ( kind = 8 ) resabs
  real ( kind = 8 ) resasc
  real ( kind = 8 ) result
  real ( kind = 8 ), parameter :: true = - 4.0D+00 / 9.0D+00

  call dqk51 ( f01, a, b, result, abserr, resabs, resasc )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DQK51_TEST'
  write ( *, '(a)' ) '  Test DQK51'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is SQRT(X)*LOG(X)'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,g14.6)' ) '  RESABS =                      ', resabs
  write ( *, '(a,g14.6)' ) '  RESASC =                      ', resasc

  return
end
subroutine dqk61_test ( )

!*****************************************************************************80
!
!! DQK61_TEST tests DQK61.
!
!  Discussion:
!
!    QK61 is a Gauss-Kronrod quadrature rule.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), parameter :: a = 0.0D+00
  real ( kind = 8 ) abserr
  real ( kind = 8 ), parameter :: b = 1.0D+00
  real ( kind = 8 ), external :: f01
  real ( kind = 8 ) resabs
  real ( kind = 8 ) resasc
  real ( kind = 8 ) result
  real ( kind = 8 ), parameter :: true = - 4.0D+00 / 9.0D+00

  call dqk61 ( f01, a, b, result, abserr, resabs, resasc )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DQK61_TEST'
  write ( *, '(a)' ) '  Test DQK61'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is SQRT(X)*LOG(X)'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,g14.6)' ) '  RESABS =                      ', resabs
  write ( *, '(a,g14.6)' ) '  RESASC =                      ', resasc

  return
end
subroutine dqng_test ( )

!*****************************************************************************80
!
!! DQNG_TEST tests DQNG.
!
!  Discussion:
!
!    DQNG is a nonadaptive automatic integrator using a Gauss-Kronrod or 
!    Patterson rule.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), parameter :: a = 0.0D+00
  real ( kind = 8 ) abserr
  real ( kind = 8 ), parameter :: b = 1.0D+00
  real ( kind = 8 ), parameter :: epsabs = 0.0D+00
  real ( kind = 8 ), parameter :: epsrel = 0.001D+00
  real ( kind = 8 ), external :: f01
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) neval
  real ( kind = 8 ) result
  real ( kind = 8 ), parameter :: true = - 4.0D+00 / 9.0D+00

  call dqng ( f01, a, b, epsabs, epsrel, result, abserr, neval, ier )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DQNG_TEST'
  write ( *, '(a)' ) '  Test DQNG'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is SQRT(X)*LOG(X)'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,i8)' ) '  Number of function evaluations, NEVAL = ', neval
  write ( *, '(a,i8)' ) '  Error return code IER = ', ier

  return
end
function f01 ( x )

!*****************************************************************************80
!
!! F01 is the integrand function SQRT(X) * LOG(X).
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) f01
  real ( kind = 8 ) x

  if ( x <= 0.0D+00 )then
    f01 = 0.0D+00
  else
    f01 = sqrt ( x ) * log ( x )
  end if

  return
end
function f02 ( x )

!*****************************************************************************80
!
!! F02 is the integrand function COS(100*SIN(X)).
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) f02
  real ( kind = 8 ) x

  f02 = cos ( 100.0D+00 * sin ( x ) )

  return
end
function f03 ( x )

!*****************************************************************************80
!
!! F03 is the integrand function LOG(X)/SQRT(X).
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) f03
  real ( kind = 8 ) x

  if ( x <= 0.0D+00 ) then
    f03 = 0.0D+00
  else
    f03 = log ( x ) / sqrt ( x )
  end if

  return
end
function f04 ( x )

!*****************************************************************************80
!
!! F04 is the integrand function X^3 LOG((X^2-1)*(X^2-2))
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) f04
  real ( kind = 8 ) x

  f04 = x**3 * log ( abs ( ( x - 1.0D+00 ) * ( x + 1.0D+00 ) * ( x**2 - 2.0D+00 ) ) )

  return
end
function f05 ( x )

!*****************************************************************************80
!
!! F05 is the integrand function LOG(X)/(1+100X^2).
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) f05
  real ( kind = 8 ) x

  f05 = log ( x ) / ( 1.0D+00 + 100.0D+00 * x**2 )

  return
end
function f06 ( x )

!*****************************************************************************80
!
!! F06 is the integrand function LOG(X).
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) f06
  real ( kind = 8 ) x

  if ( x <= 0.0D+00 ) then
    f06 = 0.0D+00
  else
    f06 = log ( x )
  end if

  return
end
function f07 ( x )

!*****************************************************************************80
!
!! F07 is the integrand function 1/SQRT(X).
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) f07
  real ( kind = 8 ) x

  if ( x <= 0.0D+00 ) then
    f07 = 0.0D+00
  else
    f07 = 1.0D+00 / sqrt ( x )
  end if

  return
end
function f08 ( x )

!*****************************************************************************80
!
!! F08 is the integrand function 1/(1+LOG(X)**2)**2
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) f08
  real ( kind = 8 ) x

  if ( 0.0D+00 < x ) then
    f08 = 1.0D+00 / ( 1.0D+00 + log ( x )**2 )**2
  else
    f08 = 0.0D+00
  end if

  return
end
function f09 ( x )

!*****************************************************************************80
!
!! F09 is the integrand function 1 / ( 5 X^3 + 6 ).
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) f09
  real ( kind = 8 ) x

  f09 = 1.0D+00 / ( 5.0D+00 * x**3 + 6.0D+00 )

  return
end
