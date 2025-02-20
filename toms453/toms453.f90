subroutine bromin ( n, s, tol, xr, xi, wr, wi, eps, ier )

!*****************************************************************************80
!
!! bromin() calculates a Bromwich quadrature rule.
!
!  this subroutine calculates abscissas and weights of the
!  gaussian quadrature formula of order n for the bromwich
!  integral.  only the abscissas of the first quadrant of
!  the complex plane, the real abscissa (if n is odd) and
!  the corresponding weights are calculated.  the other
!  abscissas and weights are complex conjugates.
!
!  input parameters
!
!    n, order of the quadrature formula.
!      n must be greater than 2.
!    tol, requested relative accuracy of the abscissas.
!    s, parameters of the weight function.
!
!  output parameters
!
!    xr and xi contain the real and imaginary parts of
!      the abscissas.  if n is odd, the real abscissa
!      is xr(1).
!    wr and wi contain the real and imaginary parts of
!      the corresponding weights.
!    eps is a crude estimation of the obtained relative
!      accuracy of the abscissas.
!    ier is an error code.
!      if ier = 0, the computation is carried out to
!        the requested accuracy.
!      if ier.gt.0, the ier-th abscissa is not found.
!      if ier = -1, the computations are carried out,
!        but the requested accuracy is not
!        achieved.
!      if ier = -2, n is less than 3.
!
!  function programs required
!    function gamma(x), which evaluates the gamma
!      function for positive x.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) ak
  real ( kind = rk ) an
  real ( kind = rk ) arg
  real ( kind = rk ) ci
  real ( kind = rk ) cr
  real ( kind = rk ) d
  real ( kind = rk ) d1
  real ( kind = rk ) d2
  real ( kind = rk ) e
  real ( kind = rk ) eps
  real ( kind = rk ) fac
  real ( kind = rk ) facti
  real ( kind = rk ) factr
  integer ier
  integer ignal
  integer j
  integer k
  integer l
  integer n1
  integer num
  integer nup
  real ( kind = rk ) pi
  real ( kind = rk ) pr
  real ( kind = rk ) qi
  real ( kind = rk ) qr
  real ( kind = rk ) ri
  real ( kind = rk ) rr
  real ( kind = rk ) s
  real ( kind = rk ) t1
  real ( kind = rk ) t2
  real ( kind = rk ) tol
  real ( kind = rk ) u
  real ( kind = rk ) v
  real ( kind = rk ) wi(n)
  real ( kind = rk ) wr(n)
  real ( kind = rk ) xi(n)
  real ( kind = rk ) xr(n)
  real ( kind = rk ) yi
  real ( kind = rk ) yr
  real ( kind = rk ) z

  if ( n < 3 ) then
    ier = - 2
    return
  end if

  n1 = ( n + 1 ) / 2
  l = n - 1
  an = n
  ier = 0
  eps = tol
  arg = 0.034D+00 * ( 30.0D+00 + an + an ) / ( an - 1.0D+00 )
  factr = cos ( arg )
  facti = sin ( arg )
  fac = 1.0D+00
  ak = 0.0D+00
  do k = 1, l
    ak = ak + 1.0D+00
    fac = - fac * ak
  end do
  fac = fac * ( an + an + s - 2.0D+00 )**2 / &
    ( an * gamma ( an + s - 1.0D+00 ) )
!
!  Approximate the first abscissa.
!
  yr = 1.333D+00 * an + s - 1.5D+00

  if ( mod ( n, 2 ) .eq. 0 ) then
    yi = 1.6D+00 + 0.07D+00 * s
  else
    yi = 0.0D+00
  end if

  do k = 1, n1

    e = tol
    ignal = 0
    num = 0
    nup = 0
!
!  Newton-Raphson method.
!
    d = yr * yr + yi * yi
    yr = yr / d
    yi = - yi / d

    do

      qr = s * yr - 1.0D+00
      qi = s * yi
      pr = ( s + 1.0D+00 ) * ( ( s + 2.0D+00 ) * ( yr * yr - yi * yi ) &
        - 2.0D+00 * yr ) + 1.0D+00
      pi = 2.0D+00 * ( s + 1.0D+00 ) * yi * ( ( s + 2.0D+00 ) * yr - 1.0D+00 )
      z = 2.0D+00

      do j = 3, n
        rr = qr
        ri = qi
        qr = pr
        qi = pi
        z = z + 1.0D+00
        u = z + s - 2.0D+00
        v = u + z
        d = ( v * yr + ( 2.0D+00 - s ) / ( v - 2.0D+00 ) ) / u
        d1 = ( z - 1.0D+00 ) * v / ( u * ( v - 2.0D+00 ) )
        d2 = v * yi / u
        pr = ( v - 1.0D+00 ) * ( qr * d - qi * d2 ) + d1 * rr
        pi = ( v - 1.0D+00 ) * ( qi * d + qr * d2 ) + d1 * ri
      end do

      if ( ignal == 1 ) then
        exit
      end if

      d = ( yr * yr + yi * yi ) * v
      d1 = ( ( pr + qr ) * yr + ( pi + qi ) * yi ) / d + pr
      d2 = ( ( pi + qi ) * yr - ( pr + qr ) * yi ) / d + pi
      d = ( d1 * d1 + d2 * d2 ) * an
      t1 = pr * yr - pi * yi
      t2 = pi * yr + pr * yi
      cr = ( t1 * d1 + t2 * d2 ) / d
      ci = ( t2 * d1 - t1 * d2 ) / d
      yr = yr - cr
      yi = yi - ci
      num = num + 1
!
!  Test of convergence of iteration process.
!
      if ( cr * cr + ci * ci <= e * e * ( yr * yr + yi * yi ) ) then

        ignal = 1
!
!  Test of number of iteration steps.
!
      else if ( 10 < num ) then

        e = e * 10.0D+00
        ier = - 1
        nup = nup + 1

        if ( 5 < nup ) then
          ier = k
          return
        end if

      end if

    end do
!
!  Calculation of weights.
!
    if ( eps < e ) then
      eps = e
    end if

    d = qr * qr + qi * qi
    d = d * d
    d1 = yr * qr + yi * qi
    d2 = yi * qr - yr * qi
    wr(k) = fac * ( d1 * d1 - d2 * d2 ) / d
    wi(k) = 2.0D+00 * fac * d2 * d1 / d
    d = yr * yr + yi * yi
    xr(k) =   yr / d
    xi(k) = - yi / d

    if ( n1 < k + 1 ) then
      exit
    end if

    if ( n1 == k + 1 ) then
      factr = cos ( 1.5D+00 * arg )
      facti = sin ( 1.5D+00 * arg )
    end if
!
!  Approximate the (K+1)-th abscissa.
!
    yr = ( xr(k) + 0.67D+00 * an ) * factr - xi(k) * facti - 0.67D+00 * an
    yi = ( xr(k) + 0.67D+00 * an ) * facti + xi(k) * factr

  end do

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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
