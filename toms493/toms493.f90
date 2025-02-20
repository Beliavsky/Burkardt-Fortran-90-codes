subroutine rpoly ( op, degree, zeror, zeroi, fail )

!*****************************************************************************80
!
!! rpoly() finds the zeros of a real polynomial.
!
!  Discussion:
!
!    In the original code, the input quantity DEGREE could be altered
!    before output.  This is an unrighteous idea and has been quashed.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 June 2024
!
!  Author:
!
!    Original Fortran77 version by Michael Jenkins.
!    This version by John Burkardt.
!
!  Reference:
!
!    Michael Jenkins,
!    Algorithm 493: Zeros of a Real Polynomial,
!    ACM Transactions on Mathematical Software,
!    Volume 1, Number 2, June 1975, pages 178-189.
!
!  Input:
!
!    real op(degree+1): coefficients in order of decreasing powers.
!
!    integer degree: the degree of the polynomial.
!
!  Output:
!
!    real zeror(degree), zeroi(degree): real and imaginary parts 
!    of the zeros.
!
!    logical fail: True only if leading coefficient is zero or if RPOLY
!    has found fewer than DEGREE zeros.  
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer degree

  real ( kind = rk8 ) a
  real ( kind = rk8 ) a1
  real ( kind = rk8 ) a3
  real ( kind = rk8 ) a7
  real ( kind = rk8 ) aa
  real ( kind = rk8 ) are
  real ( kind = rk8 ) b
  real ( kind = rk8 ) base
  real ( kind = rk8 ) bb
  real ( kind = rk8 ) bnd
  real ( kind = rk8 ) c
  real ( kind = rk8 ) cc
  integer cnt
  real ( kind = rk8 ) cosr
  real ( kind = rk8 ) d
  real ( kind = rk8 ) df
  real ( kind = rk8 ) dx
  real ( kind = rk8 ) e
  real ( kind = rk8 ) eta
  real ( kind = rk8 ) f
  real ( kind = rk8 ) factor
  logical fail
  real ( kind = rk8 ) ff
  real ( kind = rk8 ) g
  real ( kind = rk8 ) h
  integer i
  real ( kind = rk8 ) infin
  integer j
  integer jj
  real ( kind = rk8 ) k(degree+1)
  integer l
  integer l2
  real ( kind = rk8 ) lo
  real ( kind = rk8 ) lzi
  real ( kind = rk8 ) lzr
  real ( kind = rk8 ) mre
  integer n
  integer nn
  integer nz
  real ( kind = rk8 ) op(degree+1)
  real ( kind = rk8 ) p(degree+1)
  real ( kind = rk8 ) pt(degree+1)
  real ( kind = rk8 ) qk(degree+1)
  real ( kind = rk8 ) qp(degree+1)
  real ( kind = rk8 ) smalno
  real ( kind = rk8 ) svk(degree+1)
  real ( kind = rk8 ) t
  real ( kind = rk8 ) temp(degree+1)
  real ( kind = rk8 ) sc
  real ( kind = rk8 ) si
  real ( kind = rk8 ) sinr
  real ( kind = rk8 ) sr
  real ( kind = rk8 ) szi
  real ( kind = rk8 ) szr
  real ( kind = rk8 ) u
  real ( kind = rk8 ) v
  real ( kind = rk8 ) x
  real ( kind = rk8 ) xm
  real ( kind = rk8 ) xmax
  real ( kind = rk8 ) xmin
  real ( kind = rk8 ) xx
  real ( kind = rk8 ) xxx
  real ( kind = rk8 ) yy
  real ( kind = rk8 ) zeroi(degree)
  logical zerok
  real ( kind = rk8 ) zeror(degree)
!
! THE FOLLOWING STATEMENTS SET MACHINE CONSTANTS USED
! IN VARIOUS PARTS OF THE PROGRAM. THE MEANING OF THE
! FOUR CONSTANTS ARE...
! ETA     THE MAXIMUM RELATIVE REPRESENTATION ERROR
!         WHICH CAN BE DESCRIBED AS THE SMALLEST
!         POSITIVE FLOATING POINT NUMBER SUCH THAT
!         1.0+ETA IS GREATER THAN 1.
! INFINY  THE LARGEST FLOATING-POINT NUMBER.
! SMALNO  THE SMALLEST POSITIVE FLOATING-POINT NUMBER
!         IF THE EXPONENT RANGE DIFFERS IN SINGLE AND
!         real ( kind = rk8 ) THEN SMALNO AND INFIN
!         SHOULD INDICATE THE SMALLER RANGE.
!
! BASE    THE BASE OF THE FLOATING-POINT NUMBER SYSTEM USED.
!
  base = 2.0
  eta = epsilon ( eta )
  infin = huge ( infin )
  smalno = tiny ( smalno )
!
!  ARE and MRE are the unit error in addition and multiplication.
!
  are = eta
  mre = eta
  lo = smalno / eta
!
!  Initialization of constants for shift rotation.
!
  xx = 0.70710678
  yy = - xx
  cosr = -0.069756474
  sinr = 0.99756405
  fail = .false.
  n = degree
  nn = n + 1
!
!  Algorithm fails if the leading coefficient is zero.
!
  if ( op(1) == 0.0 ) then
    fail = .true.
    return
  end if
!
!  Remove the zeros at the origin if any.
!
  do while ( op(nn) == 0.0 )
    j = degree - n + 1
    zeror(j) = 0.0
    zeroi(j) = 0.0
    nn = nn - 1
    n = n - 1
  end do
!
!  Copy the coefficients.
!
  do i = 1, nn
    p(i) = op(i)
  end do
!
!  Start the algorithm for one zero.
!
   40 continue

  if ( n < 1 ) then
    return
  else if ( n == 1 ) then
    zeror(degree) = - p(2) / p(1)
    zeroi(degree) = 0.0
    return
  else if ( n == 2 ) then
    call quad ( p(1), p(2), p(3), zeror(degree-1), &
      zeroi(degree-1), zeror(degree), zeroi(degree) )
    return
  end if
!
!  Find largest and smallest moduli of coefficients.
!
  xmax = 0.0
  xmin = infin
  do i = 1, nn
    x = abs ( p(i) )
    xmax = max ( xmax, x )
    if ( x /= 0.0 ) then
      xmin = min ( xmin, x )
    end if
  end do
!
!  Scale if there are large or very small coefficients.
!
  sc = lo / xmin

  if ( 1.0 < sc ) then
    if ( xmax <= infin ) then
      l = int ( log ( sc ) / log ( base ) + 0.5 )
      factor = ( base * 1.0 ) ** l

      if ( factor /= 1.0 ) then
        do i = 1, nn
          p(i) = factor * p(i)
        end do
      end if
    end if
    go to 110
  end if

  if ( xmax < 10.0 ) then
    go to 110
  end if

  if ( sc == 0.0 ) then
    sc = smalno
  end if

  l = int ( log ( sc ) / log ( base ) + 0.5 )
  factor = ( base * 1.0 ) ** l

  if ( factor /= 1.0 ) then
    do i = 1, nn
      p(i) = factor * p(i)
    end do
  end if
!
!  Compute lower bound on moduli of zeros.
!
110   continue

  do i = 1, nn
    pt(i) = abs ( p(i) )
  end do
  pt(nn) = - pt(nn)
!
!  Compute upper estimate of bound.
!
  x = exp (( log ( -pt(nn) )- log ( pt(1) ) ) / dble ( n ) )
!
!  If Newton step at the origin is better, use it.
!
  if ( pt(n) /= 0.0 ) then
    xm = - pt(nn) / pt(n)
    x = min ( x, xm )
  end if
!
!  Chop the interval (0,X) until FF <= 0
!
  do

    xm = x * 0.1D+00
    ff = pt(1)
    do i = 2, nn
      ff = ff * xm + pt(i)
    end do

    if ( ff <= 0.0 ) then
      exit
    end if
    x = xm

  end do

  dx = x
!
!  Do Newton iteration until X converges to two decimal places.
!
  do while ( 0.005 < abs ( dx / x ) )
    ff = pt(1)
    df = ff
    do i = 2, n
      ff = ff * x + pt(i)
      df = df * x + ff
    end do
    ff = ff * x + pt(nn)
    dx = ff / df
    x = x - dx
  end do

  bnd = x
!
!  Compute the derivative as the initial K polynomial.
!
  do i = 2, n
    k(i) = dble ( nn - i ) * p(i) / dble ( n )
  end do

  k(1) = p(1)
  aa = p(nn)
  bb = p(n)
  zerok = k(n) == 0.0
!
!  Do 5 steps with no shift.
!
  do jj = 1, 5

    cc = k(n)
!
!  Use scaled form of recurrence if value of k at 0 is nonzero.
!
    if ( .not. zerok ) then

      t = - aa / cc
      do i = 1, n - 1
        j = nn - i
        k(j) = t * k(j-1) + p(j)
      end do

      k(1) = p(1)
      zerok = abs ( k(n) ) <= abs ( bb ) * eta * 10.0
!
!  Use unscaled form of recurrence.
!
    else

      do i = 1, n - 1
        j = nn - i
        k(j) = k(j-1)
      end do

      k(1) = 0.0
      zerok = k(n) == 0.0

    end if

  end do
!
!  Save K for restarts with new shifts.
!
  do i = 1, n
    temp(i) = k(i)
  end do
!
!  Loop to select the quadratic corresponding to each new shift.
!
  do cnt = 1, 20
!
!  Quadratic corresponds to a double shift to a
!  non-real point and its complex conjugate. 
!  The point has modulus BND and amplitude rotated by 94 degrees
!  from the previous shift.
!
    xxx = cosr * xx - sinr * yy
    yy =  sinr * xx + cosr * yy
    xx = xxx
    sr = bnd * xx
    si = bnd * yy
    u = -2.0 * sr
    v = bnd
!
!  Second stage calculation, fixed quadratic.
!
    l2 = 20 * cnt
    call fxshfr ( a, a1, a3, a7, are, b, c, d, e, eta, &
      f, g, h, k, l2, lzr, lzi, mre, n, nn, p, qk, qp, sr, svk, &
      szr, szi, u, v, nz )
!
!  The second stage jumps directly to one of the third stage iterations 
!  and returns here if successful.
!
!  Deflate the polynomial, store the zero or zeros and
!  return to the main algorithm.
!
    if ( 0 < nz ) then

      j = degree - n + 1
      zeror(j) = szr
      zeroi(j) = szi
      nn = nn - nz
      n = nn - 1
      do i = 1, nn
        p(i) = qp(i)
      end do

      if ( 1 < nz ) then
        zeror(j+1) = lzr
        zeroi(j+1) = lzi
      end if

      go to 40

    end if
!
!  If the iteration is unsuccessful another quadratic
!  is chosen after restoring K.
!
    do i = 1, n
      k(i) = temp(i)
    end do

  end do
!
!  Return with failure if no convergence with 20 shifts.
!
  fail = .true.

  return
end
subroutine fxshfr ( a, a1, a3, a7, are, b, c, d, e, eta, &
  f, g, h, k, l2, lzr, lzi, mre, n, nn, p, qk, qp, sr, svk, &
  szr, szi, u, v, nz )

!*****************************************************************************80
!
!! fxshfr() computes up to L2 fixed shift K-polynomials.
!
!  Discussion:
!
!    In the linear or quadratic cases, a convergence test is made.
!
!    The code initiates one of the variable shift iterations and 
!    returns the number of zeros found.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 June 2024
!
!  Author:
!
!    Original Fortran77 version by Michael Jenkins.
!    This version by John Burkardt.
!
!  Reference:
!
!    Michael Jenkins,
!    Algorithm 493: Zeros of a Real Polynomial,
!    ACM Transactions on Mathematical Software,
!    Volume 1, Number 2, June 1975, pages 178-189.
!
!  Input:
!
!    integer L2: limit of fixed shift steps
!
!  Output:
!
!    integer NZ: the number of zeros found.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer nn

  real ( kind = rk8 ) a
  real ( kind = rk8 ) a1
  real ( kind = rk8 ) a3
  real ( kind = rk8 ) a7
  real ( kind = rk8 ) are
  real ( kind = rk8 ) b
  real ( kind = rk8 ) betas
  real ( kind = rk8 ) betav
  real ( kind = rk8 ) c
  real ( kind = rk8 ) d
  logical doit
  real ( kind = rk8 ) e
  real ( kind = rk8 ) eta
  real ( kind = rk8 ) f
  real ( kind = rk8 ) g
  real ( kind = rk8 ) h
  integer i
  integer iflag
  integer j
  real ( kind = rk8 ) k(nn)
  integer l2
  real ( kind = rk8 ) lzi
  real ( kind = rk8 ) lzr
  real ( kind = rk8 ) mre
  integer n
  integer nz
  real ( kind = rk8 ) oss
  real ( kind = rk8 ) ots
  real ( kind = rk8 ) otv
  real ( kind = rk8 ) ovv
  real ( kind = rk8 ) p(nn)
  real ( kind = rk8 ) qk(nn)
  real ( kind = rk8 ) qp(nn)
  real ( kind = rk8 ) s
  logical spass
  real ( kind = rk8 ) sr
  real ( kind = rk8 ) ss
  logical stry
  real ( kind = rk8 ) svk(nn)
  real ( kind = rk8 ) svu
  real ( kind = rk8 ) svv
  real ( kind = rk8 ) szi
  real ( kind = rk8 ) szr
  real ( kind = rk8 ) ts
  real ( kind = rk8 ) tss
  real ( kind = rk8 ) tv
  real ( kind = rk8 ) tvv
  integer type
  real ( kind = rk8 ) u
  real ( kind = rk8 ) ui
  real ( kind = rk8 ) v
  real ( kind = rk8 ) vi
  logical vpass
  logical vtry
  real ( kind = rk8 ) vv

  nz = 0
  betav = 0.25
  betas = 0.25
  oss = sr
  ovv = v
!
!  Evaluate polynomial by synthetic division.
!
  call quadsd ( nn, u, v, p, qp, a, b )

  call calcsc ( nn, n, a, b, eta, k, type, u, v, a1, a3, a7, &
    c, d, e, f, g, h, qk )

  do j = 1, l2
!
!  Calculate next K polynomial and estimate V.
!
    call nextk ( nn, type, a, a1, a3, a7, b, eta, qk, qp, k )

    call calcsc ( nn, n, a, b, eta, k, type, u, v, a1, a3, a7, &
      c, d, e, f, g, h, qk )

    call newest ( nn, n, a, a1, a3, a7, b, c, d, f, g, h, &
      k, p, type, u, v, ui, vi )
    vv = vi
!
!  Estimate S.
!
    if ( k(n) == 0.0 ) then
      ss = 0.0
    else
      ss = - p(nn) / k(n)
    end if

    tv = 1.0
    ts = 1.0

    if ( j == 1 .or. type == 3 ) then
      ovv = vv
      oss = ss
      otv = tv
      ots = ts
      cycle
    end if
!
!  Compute relative measures of convergence of S and V sequences.
!
    if ( vv /= 0.0 ) then
      tv = abs ( ( vv - ovv ) / vv )
    end if

    if ( ss /= 0.0 ) then
      ts = abs ( ( ss - oss ) / ss )
    end if
!
!  If decreasing, multiply two most recent convergence measures.
!
    if ( tv < otv ) then
      tvv = tv * otv
    else
      tvv = 1.0
    end if

    if ( ts < ots ) then
      tss = ts * ots
    else
      tss = 1.0
    end if
!
!  Compare with convergence criteria.
!
    vpass = tvv < betav
    spass = tss < betas

    if ( .not. ( spass .or. vpass ) ) then
      ovv = vv
      oss = ss
      otv = tv
      ots = ts
      cycle
    end if
!
!  At least one sequence has passed the convergence test. 
!  Store variables before iterating
!
    svu = u
    svv = v
    do i = 1, n
      svk(i) = k(i)
    end do
    s = ss
!
!  Choose iteration according to the fastest converging sequence.
!
    vtry = .false.
    stry = .false.

    if ( spass .and. ( ( .not. vpass ) .or. tss < tvv ) ) then
      doit = .false.
    else
      doit = .true.
    end if

   do while ( .true. )

      if ( .not. doit ) then

        doit = .true.

      else

        call quadit ( a, a1, a3, a7, are, b, c, d, e, eta, f, g, &
          h, k, lzr, lzi, mre, n, nn, p, qk, qp, szr, szi, &
          u, v, ui, vi, nz )

        if ( 0 < nz ) then
          return
        end if
!
!  Quadratic iteration has failed. Flag that it has
!  been tried and decrease the convergence criterion.
!
        vtry = .true.
        betav = betav * 0.25
!
!  Try linear iteration if it has not been tried and
!  the S sequence is converging.
!
        if ( stry .or. ( .not. spass ) ) then

          u = svu
          v = svv
          do i = 1, n
            k(i) = svk(i)
          end do
!
!  Try quadratic iteration if it has not been tried
!  and the v sequence is converging.
!
          if ( .not. vpass .or. vtry ) then
            exit
          end if

          cycle

        end if

        do i = 1, n
          k(i) = svk(i)
        end do

      end if

      call realit ( are, eta, k, mre, n, nn, p, qk, qp, szr, szi, &
        s, nz, iflag )

      if ( 0 < nz ) then
        return
      end if
!
!  Linear iteration has failed. Flag that it has been
!  tried and decrease the convergence criterion.
!
      stry = .true.
      betas = betas * 0.25
!
!  If linear iteration signals an almost double real
!  zero, attempt quadratic interation.
!
      if ( iflag /= 0 ) then
        ui = - ( s + s )
        vi = s * s
        cycle
      end if
!
!  Restore variables.
!
      u = svu
      v = svv
      do i = 1, n
        k(i) = svk(i)
      end do
!
!  Try quadratic iteration if it has not been tried
!  and the v sequence is converging.
!
      if ( .not. vpass .or. vtry ) then
        exit
      end if

    end do
!
!  Recompute QP and scalar values to continue the second stage.
!
    call quadsd ( nn, u, v, p, qp, a, b )

    call calcsc ( nn, n, a, b, eta, k, type, u, v, a1, a3, a7, &
      c, d, e, f, g, h, qk )

    ovv = vv
    oss = ss
    otv = tv
    ots = ts

  end do

  return
end
subroutine quadit ( a, a1, a3, a7, are, b, c, d, e, eta, f, g, &
  h, k, lzr, lzi, mre, n, nn, p, qk, qp, szr, szi, u, v, uu, &
  vv, nz )

!*****************************************************************************80
!
!! quadit() applies an iteration for a quadratic factor.
!
!  Discussion:
!
!    Variable-shift K-polynomial iteration for a quadratic factor.
!
!    Converges only if the zeros are equimodular or nearly so.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 June 2024
!
!  Author:
!
!    Original Fortran77 version by Michael Jenkins.
!    This version by John Burkardt.
!
!  Reference:
!
!    Michael Jenkins,
!    Algorithm 493: Zeros of a Real Polynomial,
!    ACM Transactions on Mathematical Software,
!    Volume 1, Number 2, June 1975, pages 178-189.
!
!  Input:
!
!    real UU, VV: coefficients of starting quadratic.
!
!  Output:
!
!    integer nz: the number of zeros found.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer nn

  real ( kind = rk8 ) a
  real ( kind = rk8 ) a1
  real ( kind = rk8 ) a3
  real ( kind = rk8 ) a7
  real ( kind = rk8 ) are
  real ( kind = rk8 ) b
  real ( kind = rk8 ) c
  real ( kind = rk8 ) d
  real ( kind = rk8 ) e
  real ( kind = rk8 ) ee
  real ( kind = rk8 ) eta
  real ( kind = rk8 ) f
  real ( kind = rk8 ) g
  real ( kind = rk8 ) h
  integer i
  integer j
  real ( kind = rk8 ) k(nn)
  real ( kind = rk8 ) lzi
  real ( kind = rk8 ) lzr
  real ( kind = rk8 ) mp
  real ( kind = rk8 ) mre
  integer n
  integer nz
  real ( kind = rk8 ) omp
  real ( kind = rk8 ) p(nn)
  real ( kind = rk8 ) qk(nn)
  real ( kind = rk8 ) qp(nn)
  real ( kind = rk8 ) relstp
  real ( kind = rk8 ) szi
  real ( kind = rk8 ) szr
  real ( kind = rk8 ) t
  logical tried
  integer type
  real ( kind = rk8 ) u
  real ( kind = rk8 ) ui
  real ( kind = rk8 ) uu
  real ( kind = rk8 ) v
  real ( kind = rk8 ) vi
  real ( kind = rk8 ) vv
  real ( kind = rk8 ) zm

  nz = 0
  tried = .false.
  u = uu
  v = vv
  j = 0

  do

    a = 1.0
    call quad ( a, u, v, szr, szi, lzr, lzi )
!
!  Return if roots of the quadratic are real and not
!  close to multiple or nearly equal and of opposite sign.
!
    if ( 0.01 * abs ( lzr ) &
      < abs ( abs ( szr ) - abs ( lzr ) ) ) then
      exit
    end if
!
!  Evaluate polynomial by quadratic synthetic division.
!
    call quadsd ( nn, u, v, p, qp, a, b )
    mp = abs ( a - szr * b ) + abs ( szi * b )
!
!  Compute a rigorous bound on the rounding error in evaluating P.
!
    zm = sqrt ( abs ( v ) )
    ee = 2.0 * abs ( qp(1) )
    t = - szr * b
    do i = 2, n
      ee = ee * zm + abs ( qp(i) )
    end do
    ee = ee * zm + abs ( a + t )
    ee = ( 5.0 * mre + 4.0 * are ) * ee &
      - ( 5.0 * mre + 2.0 * are ) * &
      ( abs ( a + t ) + abs ( b ) * zm ) + 2.0 * are * abs ( t )
!
!  Iteration has converged sufficiently if the
!  polynomial value is less than 20 times this bound.
!
    if ( mp <= 20.0 * ee ) then
      nz = 2
      exit
    end if

    j = j + 1
!
!  Stop iteration after 20 steps
!
    if ( 20 < j ) then
      exit
    end if
!
!  A cluster appears to be stalling the convergence.
!  Five fixed shift steps are taken with a U, V close to the cluster.
!
    if ( 2 <= j .and. &
      relstp <= 0.01 .and. &
      omp <= mp .and. &
      .not. tried ) then

      relstp = max ( relstp, eta )
      relstp = sqrt ( relstp )
      u = u - u * relstp
      v = v + v * relstp

      call quadsd ( nn, u, v, p, qp, a, b )

      do i = 1, 5

        call calcsc ( nn, n, a, b, eta, k, type, u, v, a1, a3, a7, &
          c, d, e, f, g, h, qk )

        call nextk ( nn, type, a, a1, a3, a7, b, eta, qk, qp, k )

      end do

      tried = .true.
      j = 0

    end if

    omp = mp
!
!  Calculate next K polynomial and new U and V.
!
    call calcsc ( nn, n, a, b, eta, k, type, u, v, a1, a3, a7, &
      c, d, e, f, g, h, qk )

    call nextk ( nn, type, a, a1, a3, a7, b, eta, qk, qp, k )

    call calcsc ( nn, n, a, b, eta, k, type, u, v, a1, a3, a7, &
      c, d, e, f, g, h, qk )

    call newest ( nn, n, a, a1, a3, a7, b, c, d, f, g, h, &
      k, p, type, u, v, ui, vi )
!
!  If VI is zero, the iteration is not converging.
!
    if ( vi == 0.0 ) then
      exit
    end if

    relstp = abs ( ( vi - v ) / vi )
    u = ui
    v = vi
  
  end do

  return
end
subroutine realit ( are, eta, k, mre, n, nn, p, qk, qp, szr, szi, &
  sss, nz, iflag )

!*****************************************************************************80
!
!! realit() is a variable-shift H polynomial iteration for a real zero.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 June 2024
!
!  Author:
!
!    Original Fortran77 version by Michael Jenkins.
!    This version by John Burkardt.
!
!  Reference:
!
!    Michael Jenkins,
!    Algorithm 493: Zeros of a Real Polynomial,
!    ACM Transactions on Mathematical Software,
!    Volume 1, Number 2, June 1975, pages 178-189.
!
!  Input:
!
!    real sss: the starting iterate.
!
!  Output:
!
!    real sss: the final iterate.
!
!    integer nz: the number of zeros found.
!
!    integer iflag: indicates a pair of zeros near real axis.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n
  integer nn

  real ( kind = rk8 ) are
  real ( kind = rk8 ) ee
  real ( kind = rk8 ) eta
  integer i
  integer iflag
  integer j
  real ( kind = rk8 ) k(nn)
  real ( kind = rk8 ) kv
  real ( kind = rk8 ) mp
  real ( kind = rk8 ) mre
  real ( kind = rk8 ) ms
  integer nz
  real ( kind = rk8 ) omp
  real ( kind = rk8 ) p(nn)
  real ( kind = rk8 ) pv
  real ( kind = rk8 ) qk(nn)
  real ( kind = rk8 ) qp(nn)
  real ( kind = rk8 ) s
  real ( kind = rk8 ) sss
  real ( kind = rk8 ) szi
  real ( kind = rk8 ) szr
  real ( kind = rk8 ) t

  nz = 0
  s = sss
  iflag = 0
  j = 0

  do

    pv = p(1)
!
!  Evaluate P(S).
!
    qp(1) = pv
    do i = 2, nn
      pv = pv * s + p(i)
      qp(i) = pv
    end do
    mp = abs ( pv )
!
!  Compute a rigorous bound on the error in evaluating P.
!
    ms = abs ( s )
    ee = ( mre / ( are + mre ) ) * abs ( qp(1) )
    do i = 2, nn
      ee = ee * ms + abs ( qp(i) )
    end do
!
!  Iteration has converged sufficiently if the
!  polynomial value is less than 20 times this bound.
!
    if ( mp <= 20.0 * ( ( are + mre ) * ee - mre * mp ) ) then
      nz = 1
      szr = s
      szi = 0.0
      exit
    end if

    j = j + 1
!
!  Stop iteration after 10 steps.
!
    if ( 10 < j ) then
      exit
    end if
!
!  A cluster of zeros near the real axis has been encountered.
!  Return with IFLAG set to initiate a quadratic iteration.
!
    if ( 2 <= j .and. &
      abs ( t ) < 0.001 * abs ( s - t ) .and. &
      omp < mp ) then
      iflag = 1
      sss = s
      exit
    end if
!
!  Return if the polynomial value has increased significantly.
!
    omp = mp
!
!  Compute T, the next polynomial, and the new iterate.
!
    kv = k(1)
    qk(1) = kv
    do i = 2, n
      kv = kv * s + k(i)
      qk(i) = kv
    end do
!
!  Use the scaled form of the recurrence if the value
!  of K at S is nonzero.
!
    if ( abs ( k(n) ) * 10.0 * eta < abs ( kv ) ) then

      t = - pv / kv
      k(1) = qp(1)
      do i = 2, n
        k(i) = t * qk(i-1) + qp(i)
      end do

    else

      k(1) = 0.0
      do i = 2, n
        k(i) = qk(i-1)
      end do

    end if

    kv = k(1)
    do i = 2, n
      kv = kv * s + k(i)
    end do

    if ( abs ( kv ) <= abs ( k(n) ) * 10.0 * eta ) then
      t = 0.0
    else
      t = - pv / kv
    end if

    s = s + t

  end do

  return
end
subroutine calcsc ( nn, n, a, b, eta, k, type, u, v, a1, a3, a7, &
  c, d, e, f, g, h, qk )

!*****************************************************************************80
!
!! calcsc() calculates scalar quantities for the computation.
!
!  Discussion:
!
!    The quantities are used to compute the next K polynomial and 
!    new estimates of the quadratic coefficients.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 June 2024
!
!  Author:
!
!    Original Fortran77 version by Michael Jenkins.
!    This version by John Burkardt.
!
!  Reference:
!
!    Michael Jenkins,
!    Algorithm 493: Zeros of a Real Polynomial,
!    ACM Transactions on Mathematical Software,
!    Volume 1, Number 2, June 1975, pages 178-189.
!
!  Input:
!
!    integer TYPE: indicates how calculations are normalized 
!    to avoid overflow.
!
!    real k(nn):
!
!  Output:
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer nn

  real ( kind = rk8 ) a
  real ( kind = rk8 ) a1
  real ( kind = rk8 ) a3
  real ( kind = rk8 ) a7
  real ( kind = rk8 ) b
  real ( kind = rk8 ) c
  real ( kind = rk8 ) d
  real ( kind = rk8 ) e
  real ( kind = rk8 ) eta
  real ( kind = rk8 ) f
  real ( kind = rk8 ) g
  real ( kind = rk8 ) h
  real ( kind = rk8 ) k(nn)
  integer n
  real ( kind = rk8 ) qk(nn)
  integer type
  real ( kind = rk8 ) u
  real ( kind = rk8 ) v
!
!  Synthetic division of K by the quadratic 1,U,V.
!
  call quadsd ( n, u, v, k, qk, c, d )

  if ( abs ( c ) <= abs ( k(n) ) * 100.0 * eta .and. &
    abs ( d ) <= abs ( k(n-1) ) * 100.0 * eta ) then
    type = 3
!
!  TYPE=3 indicates the quadratic is almost a factor of K.
!
  else if ( abs ( c ) <= abs ( d ) ) then
    type = 2
!
!  TYPE=2 indicates that all formulas are divided by D.
!
    e = a / d
    f = c / d
    g = u * b
    h = v * b
    a3 = ( a + g ) * e + h * ( b / d )
    a1 = b * f - a
    a7 = ( f + u ) * a + h
!
!  TYPE=1 indicates that all formulas are divided by C
!
  else

    type = 1
    e = a / c
    f = d / c
    g = u * e
    h = v * b
    a3 = a * e + ( h / c + g ) * b
    a1 = b - a * ( d / c )
    a7 = a + g * d + h * f

  end if

  return
end
subroutine nextk ( nn, type, a, a1, a3, a7, b, eta, qk, qp, k )

!*****************************************************************************80
!
!! nextk() computes the next K polynomials using scalars from calcsc().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 June 2024
!
!  Author:
!
!    Original Fortran77 version by Michael Jenkins.
!    This version by John Burkardt.
!
!  Reference:
!
!    Michael Jenkins,
!    Algorithm 493: Zeros of a Real Polynomial,
!    ACM Transactions on Mathematical Software,
!    Volume 1, Number 2, June 1975, pages 178-189.
!
!  Input:
!
!    integer nn:
!
!    integer type:
!
!    real a:
!
!    real a1, a3, a7:
!
!    real b:
!
!    real eta:
!
!    real qk(nn):
!
!    real qp(nn):
!
!  Output:
!
!    real k(nn):
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer nn

  real ( kind = rk8 ) a
  real ( kind = rk8 ) a1
  real ( kind = rk8 ) a3
  real ( kind = rk8 ) a7
  real ( kind = rk8 ) b
  real ( kind = rk8 ) eta
  integer i
  real ( kind = rk8 ) k(nn)
  real ( kind = rk8 ) qk(nn)
  real ( kind = rk8 ) qp(nn)
  real ( kind = rk8 ) temp
  integer type

  if ( type == 1 ) then
    temp = b
  else if ( type == 2 ) then
    temp = a
  end if

  if ( type == 3 ) then

    k(1) = 0.0
    k(2) = 0.0
    do i = 3, nn
      k(i) = qk(i-2)
    end do

  else if ( abs ( temp ) * eta * 10.0 <= abs ( a1 ) ) then

    a7 = a7 / a1
    a3 = a3 / a1
    k(1) = qp(1)
    k(2) = qp(2) - a7 * qp(1)
    do i = 3, nn - 1
      k(i) = a3 * qk(i-2) - a7 * qp(i-1) + qp(i)
    end do

  else

    k(1) = 0.0
    k(2) = - a7 * qp(1)
    do i = 3, nn - 1
      k(i) = a3 * qk(i-2) - a7 * qp(i-1)
    end do

  end if

  return
end
subroutine newest ( nn, n, a, a1, a3, a7, b, c, d, f, g, h, &
  k, p, type, u, v, ui, vi )

!*****************************************************************************80
!
!! newest() computes new estimates of the quadratic coefficients.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 June 2024
!
!  Author:
!
!    Original Fortran77 version by Michael Jenkins.
!    This version by John Burkardt.
!
!  Reference:
!
!    Michael Jenkins,
!    Algorithm 493: Zeros of a Real Polynomial,
!    ACM Transactions on Mathematical Software,
!    Volume 1, Number 2, June 1975, pages 178-189.
!
!  Input:
!
!    integer type:
!
!  Output:
!
!    real ui, vi:
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) a1
  real ( kind = rk8 ) a3
  real ( kind = rk8 ) a4
  real ( kind = rk8 ) a5
  real ( kind = rk8 ) a7
  real ( kind = rk8 ) b
  real ( kind = rk8 ) b1
  real ( kind = rk8 ) b2
  real ( kind = rk8 ) c
  real ( kind = rk8 ) c1
  real ( kind = rk8 ) c2
  real ( kind = rk8 ) c3
  real ( kind = rk8 ) c4
  real ( kind = rk8 ) d
  real ( kind = rk8 ) f
  real ( kind = rk8 ) g
  real ( kind = rk8 ) h
  real ( kind = rk8 ) k(nn)
  integer n
  integer nn
  real ( kind = rk8 ) p(nn)
  real ( kind = rk8 ) temp
  integer type
  real ( kind = rk8 ) u
  real ( kind = rk8 ) ui
  real ( kind = rk8 ) v
  real ( kind = rk8 ) vi
!
!  Use formulas appropriate to setting of TYPE.
!
  if ( type == 3 ) then
    ui = 0.d0
    vi = 0.d0
    return
  end if

  if ( type == 1 ) then
    a4 = a + u * b + h * f
    a5 = c + ( u + v * f ) * d
  else
    a4 = ( a + g ) * f + h
    a5 = ( f + u ) * c + v * d
  end if
!
!  Evaluate new quadratic coefficients.
!
  b1 = - k(n) / p(nn)
  b2 = - ( k(n-1) + b1 * p(n) ) / p(nn)
  c1 = v * b2 * a1
  c2 = b1 * a7
  c3 = b1 * b1 * a3
  c4 = c1 - c2 - c3
  temp = a5 + b1 * a4 - c4
  if ( temp == 0.0 ) then
    ui = 0.0
    vi = 0.0
  else
    ui = u - ( u * ( c3 + c2 ) &
      + v * ( b1 * a1 + b2 * a7 ) ) / temp
    vi = v * ( 1.0 + c4 / temp )
  end if

  return
end
subroutine quadsd ( nn, u, v, p, q, a, b )

!*****************************************************************************80
!
!! quadsd() divides p by the quadratic 1,u,v using synthetic division.
!
!  Discussion:
!
!    The quotient is in Q and the remainder in A, B.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 June 2024
!
!  Author:
!
!    Original Fortran77 version by Michael Jenkins.
!    This version by John Burkardt.
!
!  Reference:
!
!    Michael Jenkins,
!    Algorithm 493: Zeros of a Real Polynomial,
!    ACM Transactions on Mathematical Software,
!    Volume 1, Number 2, June 1975, pages 178-189.
!
!  Input:
!
!    integer nn:
!
!    real u, v:
!
!    real p(nn): the polynomial coefficients.
!
!  Output:
!
!    real q(nn): the quotient.
!
!    real a, b: the remainder.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer nn

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  real ( kind = rk8 ) c
  integer i
  real ( kind = rk8 ) p(nn)
  real ( kind = rk8 ) q(nn)
  real ( kind = rk8 ) u
  real ( kind = rk8 ) v

  b = p(1)
  q(1) = b
  a = p(2) - u * b
  q(2) = a
  do i = 3, nn
    c = p(i) - u * a - v * b
    q(i) = c
    b = a
    a = c
  end do

  return
end
subroutine quad ( a, b1, c, sr, si, lr, li )

!*****************************************************************************80
!
!! quad() calculates the zeros of A*Z^2+B1*Z+C.
!
!  Discussion:
!
!    The quadratic formula, adjusted to avoid
!    overflow, is used to find the larger zero if the
!    zeros are real and both zeros are complex.
!
!    The smaller real zero is found directly from the
!    product of the zeros c/a.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 June 2024
!
!  Author:
!
!    Original Fortran77 version by Michael Jenkins.
!    This version by John Burkardt.
!
!  Reference:
!
!    Michael Jenkins,
!    Algorithm 493: Zeros of a Real Polynomial,
!    ACM Transactions on Mathematical Software,
!    Volume 1, Number 2, June 1975, pages 178-189.
!
!  Input:
!
!    real a, b1, c: the coefficients of a quadratic polynomial.
!
!  Output:
!
!    real sr, si: the real and imaginary parts of the smaller zero.
!
!    real lr, li: the real and imaginary parts of the larger zero.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  real ( kind = rk8 ) b1
  real ( kind = rk8 ) c
  real ( kind = rk8 ) d
  real ( kind = rk8 ) e
  real ( kind = rk8 ) li
  real ( kind = rk8 ) lr
  real ( kind = rk8 ) si
  real ( kind = rk8 ) sr

  if ( a == 0.0 ) then
    if ( b1 == 0.0 ) then
      sr = 0.0
    else
      sr = - c / b1
    end if
    lr = 0.0
    si = 0.0
    li = 0.0
    return
  end if

  if ( c == 0.0 ) then
    sr = 0.0
    lr = - b1 / a
    si = 0.0
    li = 0.0
    return
  end if
!
!  Compute discriminant avoiding overflow.
!
  b = b1 / 2.0

  if ( abs ( c ) <= abs ( b ) ) then
    e = 1.0 - ( a / b ) * ( c / b )
    d = sqrt ( abs ( e ) ) * abs ( b )
  else
    if ( c < 0.0 ) then
      e = b * ( b / abs ( c ) ) + a
    else
      e = b * ( b / abs ( c ) ) - a
    end if
    d = sqrt ( abs ( e ) ) * sqrt ( abs ( c ) )
  end if
!
!  Real zeros.
!
  if ( 0.0 <= e ) then

    if ( 0.0 <= b ) then
      d = - d
    end if
    lr = ( - b + d ) / a
    sr = 0.0
    if ( lr /= 0.0 ) then
      sr = ( c / lr ) / a
    end if
    si = 0.0
    li = 0.0
!
!  Complex conjugate zeros.
!
  else

    sr = - b / a
    lr = sr
    si = abs ( d / a )
    li = - si

  end if

  return
  end

