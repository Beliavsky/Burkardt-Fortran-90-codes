subroutine cpoly ( opr, opi, degree, zeror, zeroi, fail )

!*****************************************************************************80
!
!! cpoly() finds the zeros of a complex polynomial.
!
!  Discussion:
!
!    The program has been written to reduce the chance of overflow
!    occurring.  If it does occur, there is still a possibility that
!    the zerofinder will work provided the overflowed quantity is
!    replaced by a large number.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 March 2024
!
!  Author:
!
!    Original Fortran77 version by Michael Jenkins, Joseph Traub.
!    This version by John Burkardt.
!
!  Reference:
!
!    Michael Jenkins, Joseph Traub,
!    Algorithm 419: Zeros of a Complex Polynomial,
!    Communications of the Association for Computing Machinery,
!    Volume 15, Number 2, February 1972, pages 97-99.
!
!  Input:
!
!    real ( kind = rk8 ) opr(degree+1), opi(degree+1): the real and imaginary
!    parts of the coefficients in order of decreasing powers.
!
!    integer degree: the degree of the polynomial.
!
!  Output:
!
!    real ( kind = rk8 ) zeror(degree), zeroi(degree): the real and imaginary parts 
!    of the zeros.
!
!    logical fail:  true if the leading coefficient is zero or if cpoly
!    has found fewer than degree zeros.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer degree

  real ( kind = rk8 ) are
  real ( kind = rk8 ) base
  real ( kind = rk8 ) bnd
  real ( kind = rk8 ) cmod
  real ( kind = rk8 ) cosr
  real ( kind = rk8 ) eta
  logical fail
  real ( kind = rk8 ) hi(degree+1)
  real ( kind = rk8 ) hr(degree+1)
  integer i
  integer idnn2
  real ( kind = rk8 ) infin
  real ( kind = rk8 ) mre
  integer n
  integer nn
  real ( kind = rk8 ) opi(degree+1)
  real ( kind = rk8 ) opr(degree+1)
  real ( kind = rk8 ) pi(degree+1)
  real ( kind = rk8 ) pr(degree+1)
  real ( kind = rk8 ) pvi
  real ( kind = rk8 ) pvr
  real ( kind = rk8 ) qhi(degree+1)
  real ( kind = rk8 ) qhr(degree+1)
  real ( kind = rk8 ) qpi(degree+1)
  real ( kind = rk8 ) qpr(degree+1)
  real ( kind = rk8 ) rescale
  real ( kind = rk8 ) shr(degree+1)
  real ( kind = rk8 ) si
  real ( kind = rk8 ) sinr
  real ( kind = rk8 ) smalno
  real ( kind = rk8 ) sr
  real ( kind = rk8 ) ti
  real ( kind = rk8 ) tr
  real ( kind = rk8 ) xx
  real ( kind = rk8 ) yy
  real ( kind = rk8 ) zeroi(degree)
  real ( kind = rk8 ) zeror(degree)
  real ( kind = rk8 ) zi
  real ( kind = rk8 ) zr
!
!  Set constants.
!
  call mcon ( eta, infin, smalno, base )
  are = eta
  mre = 2.0d0 * sqrt ( 2.0d0 ) * eta
  xx = 0.70710678
  yy = -xx
  cosr = -0.060756474
  sinr = 0.99756405
  fail = .false.
  nn = degree + 1
  n = degree
!
!  Algorithm fails if the leading coefficient is zero.
!
  if ( opr(1) == 0.0d0 .and. opi(1) == 0.0d0 ) then
    fail = .true.
    return
  end if
!
!  Remove the zeros at the origin if any.
!
  do while ( opr(nn) == 0.0d0 .and. opi(nn) == 0.0d0 )
    idnn2 = degree - nn + 2
    zeror(idnn2) = 0.0d0
    zeroi(idnn2) = 0.0d0
    nn = nn - 1
  end do
!
!  Copy the coefficients.
!
  do i = 1, nn
    pr(i) = opr(i)
    pi(i) = opi(i)
  end do
!
!  Rescale the polynomial.
!
  do i = 1, nn
    shr(i) = cmod ( pr(i), pi(i) )
  end do

  bnd = rescale ( nn, shr )

  do i = 1, nn
    pr(i) = bnd * pr(i)
    pi(i) = bnd * pi(i)
  end do
!
!  Start the algorithm for one zero.
!
  fail = .false.

  do while ( .not. fail )

    if ( nn <= 2 ) then
      call cdivid ( -pr(2), -pi(2), pr(1), pi(1), zeror(degree), zeroi(degree) )
      return
    end if

    call loopy ( degree, nn, hr, hi, qhr, qhi, pr, pi, qpr, qpi, shr, pvr, &
      pvi, sr, si, tr, ti, zr, zi, cosr, sinr, xx, yy, are, mre, zeror, &
      zeroi, fail )

  end do

  return
end
subroutine loopy ( degree, nn, hr, hi, qhr, qhi, pr, pi, qpr, qpi, shr, pvr, &
  pvi, sr, si, tr, ti, zr, zi, cosr, sinr, xx, yy, are, mre, zeror, &
  zeroi, fail )

!*****************************************************************************80
!
!! loopy() applies several shifts to the polynomial to search for zeros.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 March 2024
!
!  Author:
!
!    Original Fortran77 version by Michael Jenkins, Joseph Traub.
!    This version by John Burkardt.
!
!  Reference:
!
!    Michael Jenkins, Joseph Traub,
!    Algorithm 419: Zeros of a Complex Polynomial,
!    Communications of the Association for Computing Machinery,
!    Volume 15, Number 2, February 1972, pages 97-99.
!
!  Input:
!
!  Output:
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer nn

  real ( kind = rk8 ) are
  real ( kind = rk8 ) bnd
  real ( kind = rk8 ) cauchy
  real ( kind = rk8 ) cmod
  integer cnt1
  integer cnt2
  logical conv
  real ( kind = rk8 ) cosr
  integer degree
  logical fail
  real ( kind = rk8 ) hi(nn)
  real ( kind = rk8 ) hr(nn)
  integer i
  integer idnn2
  real ( kind = rk8 ) mre
  real ( kind = rk8 ) pi(nn)
  real ( kind = rk8 ) pr(nn)
  real ( kind = rk8 ) pvi
  real ( kind = rk8 ) pvr
  real ( kind = rk8 ) qhi(nn)
  real ( kind = rk8 ) qhr(nn)
  real ( kind = rk8 ) qpi(nn)
  real ( kind = rk8 ) qpr(nn)
  real ( kind = rk8 ) shr(nn)
  real ( kind = rk8 ) si
  real ( kind = rk8 ) sinr
  real ( kind = rk8 ) sr
  real ( kind = rk8 ) ti
  real ( kind = rk8 ) tr
  real ( kind = rk8 ) xx
  real ( kind = rk8 ) xxx
  real ( kind = rk8 ) yy
  real ( kind = rk8 ) zeroi(nn)
  real ( kind = rk8 ) zeror(nn)
  real ( kind = rk8 ) zi
  real ( kind = rk8 ) zr

  fail = .false.
!
!  Calculate a lower bound on the modulus of the zeros.
!
  do i = 1, nn
    shr(i) = cmod ( pr(i), pi(i) )
  end do

  bnd = cauchy ( nn, shr )
!
!  Outer loop to control 2 major passes with different sequences of shifts.
!
  do cnt1 = 1, 2
!
!  First stage calculation, no shift.
!
    call noshft ( nn, 5, pr, pi, tr, ti, hr, hi )
!
!  Inner loop to select a shift.
!
    do cnt2 = 1, 9
!
!  Shift is chosen with modulus bnd and amplitude rotated by
!  94 degrees from the previous shift.
!
      xxx = cosr * xx - sinr * yy
      yy  = sinr * xx + cosr * yy
      xx = xxx
      sr = bnd * xx
      si = bnd * yy
!
!  Second stage calculation, fixed shift.
!
      call fxshft ( nn, 10 * cnt2, are, mre, pr, pi, hr, hi, qpr, qpi, &
        qhr, qhi, sr, si, tr, ti, pvr, pvi, zr, zi, conv )
!
!  The second stage jumps directly to the third stage iteration.
!  If successful the zero is stored and the polynomial deflated.
!
      if ( conv ) then
        idnn2 = degree - nn + 2
        zeror(idnn2) = zr
        zeroi(idnn2) = zi
        nn = nn - 1
        do i = 1, nn
          pr(i) = qpr(i)
          pi(i) = qpi(i)
        end do
        return
      end if
!
!  If the iteration is unsuccessful another shift is chosen.
!
    end do
!
!  If 9 shifts fail, the outer loop is repeated with another
!  sequence of shifts.
!
  end do
!
!  The zerofinder has failed on two major passes.
!
  fail = .true.

  return
end
subroutine noshft ( nn, l1, pr, pi, tr, ti, hr, hi )

!*****************************************************************************80
!
!! noshft() computes the derivative polynomial as the initial h polynomial.
!
!  Discussion:
!
!    It also computes l1 no-shift h polynomials.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 March 2024
!
!  Author:
!
!    Original Fortran77 version by Michael Jenkins, Joseph Traub.
!    This version by John Burkardt.
!
!  Reference:
!
!    Michael Jenkins, Joseph Traub,
!    Algorithm 419: Zeros of a Complex Polynomial,
!    Communications of the Association for Computing Machinery,
!    Volume 15, Number 2, February 1972, pages 97-99.
!
!  Input:
!
!    integer l1: the number of no-shift h polynomials to compute.
!
!  Output:
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer nn

  real ( kind = rk8 ) base
  real ( kind = rk8 ) cmod
  real ( kind = rk8 ) eta
  real ( kind = rk8 ) hi(nn)
  real ( kind = rk8 ) hr(nn)
  integer i
  real ( kind = rk8 ) infin
  integer j
  integer jj
  integer l1
  integer n
  integer nm1
  real ( kind = rk8 ) pi(nn)
  real ( kind = rk8 ) pr(nn)
  real ( kind = rk8 ) smalno
  real ( kind = rk8 ) t1
  real ( kind = rk8 ) t2
  real ( kind = rk8 ) ti
  real ( kind = rk8 ) tr
  real ( kind = rk8 ) xni

  call mcon ( eta, infin, smalno, base )

  n = nn - 1
  nm1 = n - 1
  do i = 1, n
    xni = nn - i
    hr(i) = xni * pr(i) / float ( n )
    hi(i) = xni * pi(i) / float ( n )
  end do

  do jj = 1, l1

    if ( eta * 10.0d0 * cmod ( pr(n), pi(n) ) < cmod ( hr(n), hi(n) ) ) then
 
      call cdivid ( -pr(nn), -pi(nn), hr(n), hi(n), tr, ti )

      do i = 1, nm1
        j = nn - i
        t1 = hr(j-1)
        t2 = hi(j-1)
        hr(j) = tr * t1 - ti * t2 + pr(j)
        hi(j) = tr * t2 + ti * t1 + pi(j)
      end do
      hr(1) = pr(1)
      hi(1) = pi(1)
!
!  If the constant term is essentially zero, shift h coefficients.
!
    else

      do i = 1, nm1
    j = nn - i
    hr(j) = hr(j-1)
    hi(j) = hi(j-1)
      end do
      hr(1) = 0.0d0
      hi(1) = 0.0d0

    end if

  end do

  return
end
subroutine fxshft ( nn, l2, are, mre, pr, pi, hr, hi, qpr, qpi, qhr, qhi, &
  sr, si, tr, ti, pvr, pvi, zr, zi, conv )

!*****************************************************************************80
!
!! fxshft() computes l2 fixed-shift h polynomials and tests for convergence.
!
!  Discussion:
!
!    It initiates a variable-shift iteration and returns with the
!    approximate zero if successful.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 March 2024
!
!  Author:
!
!    Original Fortran77 version by Michael Jenkins, Joseph Traub.
!    This version by John Burkardt.
!
!  Reference:
!
!    Michael Jenkins, Joseph Traub,
!    Algorithm 419: Zeros of a Complex Polynomial,
!    Communications of the Association for Computing Machinery,
!    Volume 15, Number 2, February 1972, pages 97-99.
!
!  Input:
!
!    integer nn: ?
!
!    integer l2: the maximum number of fixed shift steps to take.
!
!    real ( kind = rk8 ) are: error bound on complex addition.
!
!    real ( kind = rk8 ) mre: error bound on complex multiplication.
!
!    real ( kind = rk8 ) pr(nn), pi(nn): ?
!
!    real ( kind = rk8 ) hr(nn), hi(nn): ?
!
!    real ( kind = rk8 ) qpr(nn), qpi(nn): ?
!
!    real ( kind = rk8 ) qhr(nn), qhi(nn): ?
!
!    real ( kind = rk8 ) sr, si: ?
!
!    real ( kind = rk8 ) tr, ti: ?
!
!    real ( kind = rk8 ) pvr, pvi: ?
!
!  Output:
!
!    real ( kind = rk8 ) zr, zi: approximate zero if conv is true.
!
!    logical conv: true if the stage 3 iteration converged.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n
  integer nn

  real ( kind = rk8 ) are
  logical bool
  real ( kind = rk8 ) cmod
  logical conv
  real ( kind = rk8 ) hi(nn)
  real ( kind = rk8 ) hr(nn)
  integer i
  integer j
  integer l2
  real ( kind = rk8 ) mre
  real ( kind = rk8 ) oti
  real ( kind = rk8 ) otr
  logical pasd
  real ( kind = rk8 ) pi(nn)
  real ( kind = rk8 ) pr(nn)
  real ( kind = rk8 ) pvi
  real ( kind = rk8 ) pvr
  real ( kind = rk8 ) qhi(nn)
  real ( kind = rk8 ) qhr(nn)
  real ( kind = rk8 ) qpi(nn)
  real ( kind = rk8 ) qpr(nn)
  real ( kind = rk8 ) shi(nn)
  real ( kind = rk8 ) shr(nn)
  real ( kind = rk8 ) si
  real ( kind = rk8 ) sr
  real ( kind = rk8 ) svsi
  real ( kind = rk8 ) svsr
  logical test
  real ( kind = rk8 ) ti
  real ( kind = rk8 ) tr
  real ( kind = rk8 ) zi
  real ( kind = rk8 ) zr

  n = nn - 1
!
!  Evaluate polynomial at s.
!
  call polyev ( nn, sr, si, pr, pi, qpr, qpi, pvr, pvi )
  test = .true.
  pasd = .false.
!
!  Calculate first t = -p(s)/h(s).
!
  call calct ( nn, are, bool, sr, si, hr, hi, qhr, qhi, pvr, pvi, tr, ti )
!
!  Main loop for one second stage step.
!
  do j = 1, l2

    otr = tr
    oti = ti
!
!  Compute next h polynomial and new t.
!
    call nexth ( bool, nn, qhr, qhi, qpr, qpi, tr, ti, hr, hi )
    call calct ( nn, are, bool, sr, si, hr, hi, qhr, qhi, pvr, pvi, tr, ti )
    zr = sr + tr
    zi = si + ti
!
!  Test for convergence unless stage 3 has failed once or this
!  is the last h polynomial .
!
    if ( bool .or. .not. test .or. j == l2 ) then
      cycle
    end if

    if ( 0.5d0 * cmod ( zr, zi ) <= cmod ( tr - otr, ti - oti ) ) then
      pasd = .false.
      cycle
    end if

    if ( .not. pasd ) then
      pasd = .true.
      cycle
    end if
!
!  The weak convergence test has been passed twice.
!  Save the current h polynomial and shift.
!
    do i = 1, n
      shr(i) = hr(i)
      shi(i) = hi(i)
    end do
    svsr = sr
    svsi = si
!
!  Start the third stage iteration. 
!
    call vrshft ( nn, 10, are, mre, sr, si, pr, pi, qpr, qpi, pvr, pvi, &
      hr, hi, qhr, qhi, tr, ti, zr, zi, conv )

    if ( conv ) then
      return
    end if
!
!  The iteration failed to converge. 
!  Turn off testing and restore h, s, pv and t.
!
    test = .false.

    do i = 1, n
      hr(i) = shr(i)
      hi(i) = shi(i)
    end do
    sr = svsr
    si = svsi
    call polyev ( nn, sr, si, pr, pi, qpr, qpi, pvr, pvi )
    call calct ( nn, are, bool, sr, si, hr, hi, qhr, qhi, pvr, pvi, tr, ti )

  end do
!
!  Attempt an iteration with final h polynomial from second stage.
!
  call vrshft ( nn, 10, are, mre, sr, si, pr, pi, qpr, qpi, pvr, pvi, &
    hr, hi, qhr, qhi, tr, ti, zr, zi, conv )

  return
end
subroutine vrshft ( nn, l3, are, mre, sr, si, pr, pi, qpr, qpi, pvr, pvi, &
  hr, hi, qhr, qhi, tr, ti, zr, zi, conv )

!*****************************************************************************80
!
!! vrshft() carries out the third stage iteration.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 March 2024
!
!  Author:
!
!    Original Fortran77 version by Michael Jenkins, Joseph Traub.
!    This version by John Burkardt.
!
!  Reference:
!
!    Michael Jenkins, Joseph Traub,
!    Algorithm 419: Zeros of a Complex Polynomial,
!    Communications of the Association for Computing Machinery,
!    Volume 15, Number 2, February 1972, pages 97-99.
!
!  Input:
!
!    integer l3: the maximum number of steps to take in stage 3.
!
!    real ( kind = rk8 ) are: error bound on complex addition.
!
!    real ( kind = rk8 ) mre: error bound on complex multiplication.
!
!    real ( kind = rk8 ) zr, zi: the real and imaginary parts of the 
!    initial iterate.
!
!  Output:
!
!    real ( kind = rk8 ) zr, zi: the real and imaginary parts of the
!    final iterate.
!
!    logical conv: true if the iteration converged.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer nn

  real ( kind = rk8 ) are
  logical b
  real ( kind = rk8 ) base
  logical bool
  real ( kind = rk8 ) cmod
  logical conv
  real ( kind = rk8 ) errev
  real ( kind = rk8 ) eta
  real ( kind = rk8 ) hi(nn)
  real ( kind = rk8 ) hr(nn)
  integer i
  real ( kind = rk8 ) infin
  integer j
  integer l3
  real ( kind = rk8 ) mp
  real ( kind = rk8 ) mre
  real ( kind = rk8 ) ms
  real ( kind = rk8 ) omp
  real ( kind = rk8 ) pi(nn)
  real ( kind = rk8 ) pr(nn)
  real ( kind = rk8 ) pvi
  real ( kind = rk8 ) pvr
  real ( kind = rk8 ) qhi(nn)
  real ( kind = rk8 ) qhr(nn)
  real ( kind = rk8 ) qpi(nn)
  real ( kind = rk8 ) qpr(nn)
  real ( kind = rk8 ) r1
  real ( kind = rk8 ) r2
  real ( kind = rk8 ) relstp
  real ( kind = rk8 ) si
  real ( kind = rk8 ) smalno
  real ( kind = rk8 ) sr
  real ( kind = rk8 ) tp
  real ( kind = rk8 ) ti
  real ( kind = rk8 ) tr
  real ( kind = rk8 ) zi
  real ( kind = rk8 ) zr

  call mcon ( eta, infin, smalno, base )

  conv = .false.
  b = .false.
  sr = zr
  si = zi
!
!  Main loop for stage three
! 
  do i = 1, l3
!
!  Evaluate polynomial at s and test for convergence.
!
    call polyev ( nn, sr, si, pr, pi, qpr, qpi, pvr, pvi )
    mp = cmod ( pvr, pvi )
    ms = cmod ( sr, si )
!
!  Polynomial value is smaller in value than a bound on the error
!  in evaluating p, terminate the iteration.
!
    if ( mp < 20.0d0 * errev ( nn, qpr, qpi, ms, mp, are, mre ) ) then
      conv = .true.
      zr = sr
      zi = si
      return
    end if

    if ( i == 1 ) then
      omp = mp
      cycle
    end if

    if ( b .or. mp < omp .or. 0.05d0 <= relstp ) then

      if ( omp < mp * 0.1d0 ) then
    return
      end if
      omp = mp
!
!  Iteration has stalled. probably a cluster of zeros. do 5 fixed
!  shift steps into the cluster to force one zero to dominate.
!
    else

      tp = relstp
      b = .true.
      if ( relstp < eta ) then
    tp = eta
      end if
      r1 = sqrt ( tp )
      r2 = sr * ( 1.0d0 + r1 ) - si * r1
      si = sr * r1 + si * ( 1.0d0 + r1 )
      sr = r2
      call polyev ( nn, sr, si, pr, pi, qpr, qpi, pvr, pvi )

      do j = 1, 5
        call calct ( nn, are, bool, sr, si, hr, hi, qhr, qhi, pvr, pvi, tr, ti )
        call nexth ( bool, nn, qhr, qhi, qpr, qpi, tr, ti, hr, hi )
      end do

      omp = infin

    end if

    call calct ( nn, are, bool, sr, si, hr, hi, qhr, qhi, pvr, pvi, tr, ti )
    call nexth ( bool, nn, qhr, qhi, qpr, qpi, tr, ti, hr, hi )
    call calct ( nn, are, bool, sr, si, hr, hi, qhr, qhi, pvr, pvi, tr, ti )

    if ( .not. bool ) then
      relstp = cmod ( tr, ti ) / cmod ( sr, si )
      sr = sr + tr
      si = si + ti
    end if

  end do

  return
end
subroutine calct ( nn, are, bool, sr, si, hr, hi, qhr, qhi, pvr, pvi, tr, ti )

!*****************************************************************************80
!
!! calct() computes t = -p(s)/h(s).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 March 2024
!
!  Author:
!
!    Original Fortran77 version by Michael Jenkins, Joseph Traub.
!    This version by John Burkardt.
!
!  Reference:
!
!    Michael Jenkins, Joseph Traub,
!    Algorithm 419: Zeros of a Complex Polynomial,
!    Communications of the Association for Computing Machinery,
!    Volume 15, Number 2, February 1972, pages 97-99.
!
!  Input:
!
!    real ( kind = rk8 ) are: error bound on complex addition.
!
!  Output:
!
!    logical bool: set true if h(s) is essentially zero.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer nn

  real ( kind = rk8 ) are
  logical bool
  real ( kind = rk8 ) cmod
  real ( kind = rk8 ) hi(nn)
  real ( kind = rk8 ) hr(nn)
  real ( kind = rk8 ) hvi
  real ( kind = rk8 ) hvr
  integer n
  real ( kind = rk8 ) pvi
  real ( kind = rk8 ) pvr
  real ( kind = rk8 ) qhi(nn)
  real ( kind = rk8 ) qhr(nn)
  real ( kind = rk8 ) si
  real ( kind = rk8 ) sr
  real ( kind = rk8 ) ti
  real ( kind = rk8 ) tr

  n = nn - 1
!
!  Evaluate h(s).
!
  call polyev ( n, sr, si, hr, hi, qhr, qhi, hvr, hvi )

  bool = cmod ( hvr, hvi ) <= are * 10.0d0 * cmod ( hr(n), hi(n) )

  if ( .not. bool ) then
    call cdivid ( -pvr, -pvi, hvr, hvi, tr, ti )
  else
    tr = 0.0d0
    ti = 0.0d0
  end if

  return
end
subroutine nexth ( bool, nn, qhr, qhi, qpr, qpi, tr, ti, hr, hi )

!*****************************************************************************80
!
!! nexth() calculates the next shifted h polynomial.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 March 2024
!
!  Author:
!
!    Original Fortran77 version by Michael Jenkins, Joseph Traub.
!    This version by John Burkardt.
!
!  Reference:
!
!    Michael Jenkins, Joseph Traub,
!    Algorithm 419: Zeros of a Complex Polynomial,
!    Communications of the Association for Computing Machinery,
!    Volume 15, Number 2, February 1972, pages 97-99.
!
!  Input:
!
!    logical bool: if true, h(s) is essentially zero.
!
!  Output:
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer nn

  logical bool
  real ( kind = rk8 ) hr(nn)
  real ( kind = rk8 ) hi(nn)
  integer j
  integer n
  integer nm1
  real ( kind = rk8 ) qhr(nn)
  real ( kind = rk8 ) qhi(nn)
  real ( kind = rk8 ) qpr(nn)
  real ( kind = rk8 ) qpi(nn)
  real ( kind = rk8 ) t1
  real ( kind = rk8 ) t2
  real ( kind = rk8 ) ti
  real ( kind = rk8 ) tr

  n = nn - 1
  nm1 = n - 1

  if ( .not. bool ) then

    do j = 2, n
      t1 = qhr(j-1)
      t2 = qhi(j-1)
      hr(j) = tr * t1 - ti * t2 + qpr(j)
      hi(j) = tr * t2 + ti * t1 + qpi(j)
    end do
    hr(1) = qpr(1)
    hi(1) = qpi(1)
!
!  If h(s) is zero, replace h with qh.
!
  else

    do j = 2, n
      hr(j) = qhr(j-1)
      hi(j) = qhi(j-1)
    end do
    hr(1) = 0.0d0
    hi(1) = 0.0d0

  end if

  return
end
subroutine polyev ( nn, sr, si, pr, pi, qr, qi, pvr, pvi )

!*****************************************************************************80
!
!! polyev() evaluates a polynomial by the Horner recurrence.
!
!  Discussion:
!
!    It places the partial sums in q and the computed value in pv.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 March 2024
!
!  Author:
!
!    Original Fortran77 version by Michael Jenkins, Joseph Traub.
!    This version by John Burkardt.
!
!  Reference:
!
!    Michael Jenkins, Joseph Traub,
!    Algorithm 419: Zeros of a Complex Polynomial,
!    Communications of the Association for Computing Machinery,
!    Volume 15, Number 2, February 1972, pages 97-99.
!
!  Input:
!
!    integer nn: the polynomial order.
!
!    double precison sr, si: the real and imaginary parts of the
!    polynomial argument.
!
!    double precison pr(nn), pi(nn): the real and imaginary parts of the
!    polynomial coefficients.
!
!  Output:
!
!    real ( kind = rk8 ) qr(nn), qi(nn): the partial sums.
!
!    real ( kind = rk8 ) pvr, pvi: the real and imaginary parts of the
!    polynomial value.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer nn

  integer i
  real ( kind = rk8 ) pi(nn)
  real ( kind = rk8 ) pr(nn)
  real ( kind = rk8 ) pvi
  real ( kind = rk8 ) pvr
  real ( kind = rk8 ) qi(nn)
  real ( kind = rk8 ) qr(nn)
  real ( kind = rk8 ) si
  real ( kind = rk8 ) sr
  real ( kind = rk8 ) t

  qr(1) = pr(1)
  qi(1) = pi(1)
  pvr = qr(1)
  pvi = qi(1)
  do i = 2, nn
    t = pvr * sr - pvi * si + pr(i)
    pvi = pvr * si + pvi * sr + pi(i)
    pvr = t
    qr(i) = pvr
    qi(i) = pvi
  end do

  return
end
function errev ( nn, qr, qi, ms, mp, are, mre )

!*****************************************************************************80
!
!! errev() bounds the error evaluating the polynomial by the Horner recurrence.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 March 2024
!
!  Author:
!
!    Original Fortran77 version by Michael Jenkins, Joseph Traub.
!    This version by John Burkardt.
!
!  Reference:
!
!    Michael Jenkins, Joseph Traub,
!    Algorithm 419: Zeros of a Complex Polynomial,
!    Communications of the Association for Computing Machinery,
!    Volume 15, Number 2, February 1972, pages 97-99.
!
!  Input:
!
!    integer nn: the order of the polynomial.
!
!    real ( kind = rk8 ) qr(nn), qi(nn): the partial sums.
!
!    real ( kind = rk8 ) ms: the modulus of the points.
!
!    real ( kind = rk8 ) mp: the modulus of the polynomial value.
!
!    real ( kind = rk8 ) are, mre: error bounds on complex addition 
!    and multiplication.
!
!  Output:
!
!    real ( kind = rk8 ) errev: the polynomial evaluation error bound.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer nn

  real ( kind = rk8 ) are
  real ( kind = rk8 ) cmod
  real ( kind = rk8 ) e
  real ( kind = rk8 ) errev
  integer i
  real ( kind = rk8 ) mp
  real ( kind = rk8 ) mre
  real ( kind = rk8 ) ms
  real ( kind = rk8 ) qi(nn)
  real ( kind = rk8 ) qr(nn)

  e = cmod ( qr(1), qi(1) ) * mre / ( are + mre )
  do i = 1, nn
    e = e * ms + cmod ( qr(i), qi(i) )
  end do

  errev = e * ( are + mre ) - mp * mre

  return
end
function cauchy ( nn, pt )

!*****************************************************************************80
!
!! cauchy() computes a lower bound on the moduli of the zeros of a polynomial.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 March 2024
!
!  Author:
!
!    Original Fortran77 version by Michael Jenkins, Joseph Traub.
!    This version by John Burkardt.
!
!  Reference:
!
!    Michael Jenkins, Joseph Traub,
!    Algorithm 419: Zeros of a Complex Polynomial,
!    Communications of the Association for Computing Machinery,
!    Volume 15, Number 2, February 1972, pages 97-99.
!
!  Input:
!
!    integer nn: the order of the polynomial.
!
!    real ( kind = rk8 ) pt(nn): the modulus of the coefficients.
!
!  Output:
!
!    real ( kind = rk8 ) cauchy: a lower bound on the moduli of the zeros.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer nn

  real ( kind = rk8 ) cauchy
  real ( kind = rk8 ) df
  real ( kind = rk8 ) dx
  real ( kind = rk8 ) f
  integer i
  integer n
  real ( kind = rk8 ) pt(nn)
  real ( kind = rk8 ) q(nn)
  real ( kind = rk8 ) x
  real ( kind = rk8 ) xm

  pt(nn) = - pt(nn)
!
!  Compute upper estimate of bound.
!
  n = nn - 1
  x = exp ( ( log ( - pt(nn) ) - log ( pt(1) ) ) / float ( n ) )
!
!  If Newton step at the origin is better, use it.
!
  if ( pt(n) /= 0.0d0 ) then
    xm = - pt(nn) / pt(n)
    x = min ( x, xm )
  end if
!
!  Chop the interval (0,x) unitl f <= 0.
!
  do

    xm = x * 0.1d0

    f = pt(1)
    do i = 2, nn
      f = f * xm + pt(i)
    end do

    if ( f <= 0.0d0 ) then
      exit
    end if

    x = xm
 
  end do

  dx = x
!
!  Do Newton iteration until x converges to two decimal places.
!
  do while ( 0.005D0 < abs ( dx / x ) )

    q(1) = pt(1)
    do i = 2, nn
      q(i) = q(i-1) * x + pt(i)
    end do
    f = q(nn)

    df = q(1)
    do i = 2, n
      df = df * x + q(i)
    end do

    dx = f / df
    x = x - dx

  end do

  cauchy = x

  return
end
function rescale ( nn, pt )

!*****************************************************************************80
!
!! rescale() returns a scale factor for the coefficients of the polynomial. 
!
!  Discussion:
!
!    The scaling is done to avoid overflow and to avoid undetected underflow 
!    interfering with the convergence criterion.  
!
!    The factor is a power of the base.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 March 2024
!
!  Author:
!
!    Original Fortran77 version by Michael Jenkins, Joseph Traub.
!    This version by John Burkardt.
!
!  Reference:
!
!    Michael Jenkins, Joseph Traub,
!    Algorithm 419: Zeros of a Complex Polynomial,
!    Communications of the Association for Computing Machinery,
!    Volume 15, Number 2, February 1972, pages 97-99.
!
!  Input:
!
!    integer nn: the order of the polynomial.
!
!    real ( kind = rk8 ) pt(nn): modulus of coefficients of the polynomial.
!
!  Output:
!
!    real ( kind = rk8 ) rescale: the scale factor.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer nn

  real ( kind = rk8 ) base
  real ( kind = rk8 ) eta
  real ( kind = rk8 ) hi
  integer i
  real ( kind = rk8 ) infin
  integer l
  real ( kind = rk8 ) lo
  real ( kind = rk8 ) pt(nn)
  real ( kind = rk8 ) rescale
  real ( kind = rk8 ) sc
  real ( kind = rk8 ) smalno
  real ( kind = rk8 ) x
  real ( kind = rk8 ) xmax
  real ( kind = rk8 ) xmin
!
!  Get machine constants.
!
  call mcon ( eta, infin, smalno, base )
!
!  Find largest and smallest moduli of coefficients.
!
  hi = sqrt ( infin )
  lo = smalno / eta
  xmax = 0.0d0
  xmin = infin

  do i = 1, nn
    x = pt(i)
    xmax = max ( xmax, x )
    if ( x /= 0.0d0 ) then
      xmin = min ( xmin, x )
    end if
  end do
!
!  Scale only if there are very large or very small components.
!
  rescale = 1.0d0
  if ( lo <= xmin .and. xmax <= hi ) then
    return
  end if

  x = lo / xmin
  if ( x <= 1.0d0 ) then
    sc = 1.0d0 / ( sqrt ( xmax ) * sqrt ( xmin ) )
  else
    sc = x
    if ( xmax < infin / sc ) then
      sc = 1.0d0
    end if
  end if

  l = int ( log ( sc ) / log ( base ) + 0.5d0 )
  rescale = base ** l

  return
end
subroutine cdivid ( ar, ai, br, bi, cr, ci )

!*****************************************************************************80
!
!! cdivid() performs complex division c = a/b, avoiding overflow.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 March 2024
!
!  Author:
!
!    Original Fortran77 version by Michael Jenkins, Joseph Traub.
!    This version by John Burkardt.
!
!  Reference:
!
!    Michael Jenkins, Joseph Traub,
!    Algorithm 419: Zeros of a Complex Polynomial,
!    Communications of the Association for Computing Machinery,
!    Volume 15, Number 2, February 1972, pages 97-99.
!
!  Input:
!
!    real ( kind = rk8 ) ar, ai: the real and imaginary part of the numerator.
!
!    real ( kind = rk8 ) br, bi: the real and imaginary part of the denominator.
!
!  Output:
!
!    real ( kind = rk8 ) cr, ci: the real and imaginary part of the result.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) ar
  real ( kind = rk8 ) ai
  real ( kind = rk8 ) base
  real ( kind = rk8 ) br
  real ( kind = rk8 ) bi
  real ( kind = rk8 ) cr
  real ( kind = rk8 ) ci
  real ( kind = rk8 ) d
  real ( kind = rk8 ) eta
  real ( kind = rk8 ) infin
  real ( kind = rk8 ) r
  real ( kind = rk8 ) smalno
!
!  Division by zero, c = infinity.
!
  if ( br == 0.0d0 .and. bi == 0.0d0 ) then

    call mcon ( eta, infin, smalno, base )

    cr = infin
    ci = infin

  else if ( abs ( br ) < abs ( bi ) ) then

    r = br / bi
    d = bi + r * br
    cr = ( ar * r + ai ) / d
    ci = ( ai * r - ar ) / d

  else

    r = bi / br
    d = br + r * bi
    cr = ( ar + ai * r ) / d
    ci = ( ai - ar * r ) / d

  end if

  return
end
function cmod ( r, i )

!*****************************************************************************80
!
!! cmod() computes the modulus of a complex number, avoiding overflow.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 March 2024
!
!  Author:
!
!    Original Fortran77 version by Michael Jenkins, Joseph Traub.
!    This version by John Burkardt.
!
!  Reference:
!
!    Michael Jenkins, Joseph Traub,
!    Algorithm 419: Zeros of a Complex Polynomial,
!    Communications of the Association for Computing Machinery,
!    Volume 15, Number 2, February 1972, pages 97-99.
!
!  Input:
!
!    real ( kind = rk8 ) r, i: the real and imaginary parts of a complex number.
!
!  Output:
!
!    real ( kind = rk8 ) cmod: the modulus of the complex number.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) ai
  real ( kind = rk8 ) ar
  real ( kind = rk8 ) cmod
  real ( kind = rk8 ) i
  real ( kind = rk8 ) r

  ar = abs ( r )
  ai = abs ( i )

  if ( ar < ai ) then
    cmod = ai * sqrt ( 1.0d0 + ( ar / ai )**2 )
  else if ( ai < ar ) then
    cmod = ar * sqrt ( 1.0d0 + ( ai / ar )**2 )
  else
    cmod = ar * sqrt ( 2.0d0 )
  end if

  return
end
subroutine mcon ( eta, infin, smalno, base )

!*****************************************************************************80
!
!! mcon() provides machine constants.
!
!  Discussion:
!
!    The user may either set them directly or use the
!    statements below to compute them. 
!
!    Let t be the number of base-digits in each floating-point
!    number(real ( kind = rk8 )). Then eta is either .5*b**(1-t)
!    or b**(1-t) depending on whether rounding or truncation
!    is used.
!
!    Let m be the largest exponent and n the smallest exponent
!    in the number system. Then infiny is (1-base**(-t))*base**m
!    and smalno is base**n.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 March 2024
!
!  Author:
!
!    Original Fortran77 version by Michael Jenkins, Joseph Traub.
!    This version by John Burkardt.
!
!  Reference:
!
!    Michael Jenkins, Joseph Traub,
!    Algorithm 419: Zeros of a Complex Polynomial,
!    Communications of the Association for Computing Machinery,
!    Volume 15, Number 2, February 1972, pages 97-99.
!
!  Output:
!
!    real ( kind = rk8 ) eta: the smallest positive floating-point number such 
!    that 1.0d0 + eta is greater than 1.0d0.
!
!    real ( kind = rk8 ) infin: the largest floating-point number.
!
!    real ( kind = rk8 ) smalno: the smallest positive floating-point number.
!
!    real ( kind = rk8 ) base: the base of the floating-point number system used.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) base
  real ( kind = rk8 ) eta
  real ( kind = rk8 ) infin
  integer m
  integer n
  real ( kind = rk8 ) smalno
  integer t

  base = 16.0d0
  t = 14
  m = 63
  n = -65
  eta = base ** ( 1 - t )
  infin = base * ( 1.0d0 - base ** ( - t ) ) * base ** ( m - 1 )
  smalno = ( base ** ( n + 3 ) ) / base ** 3

  return
  end
