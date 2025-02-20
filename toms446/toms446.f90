subroutine binom ( x, xx, npl, m, nt, xa )

!*****************************************************************************80
!
!! binom(): binomial expansion series for the (-1/M) power of a Chebyshev series.
!
!  Discussion:
!
!    This routine uses a certain number of terms of the binomial expansion 
!    series to estimate the (-1/M) power of a given Chebyshev series. 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Roger Broucke.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Roger Broucke,
!    Algorithm 446:
!    Ten Subroutines for the Manipulation of Chebyshev Series,
!    Communications of the ACM,
!    October 1973, Volume 16, Number 4, pages 254-256.
!
!  Parameters:
!
!    Input, real ( kind = rk ) X(NPL), the given Chebyshev series.
!
!    Input, real ( kind = rk ) XX(NPL), an initial estimate for
!    the Chebyshev series for the input function raised to the (-1/M) power.
!
!    Input, integer NPL, the number of terms in the 
!    Chebyshev series.
!
!    Input, integer M, defines the exponent, (-1/M).
!    0 < M.
!
!    Input, integer NT, the number of terms of the binomial
!    series to be used.
!
!    Output, real ( kind = rk ) XA(NPL), the estimated Chebyshev series
!    for the input function raised to the (-1/M) power.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer npl

  real ( kind = rk ) alfa
  real ( kind = rk ) coef
  real ( kind = rk ) dkm2
  real ( kind = rk ) dkmm
  real ( kind = rk ) dm
  integer k
  integer m
  integer nt
  real ( kind = rk ) w2(npl)
  real ( kind = rk ) w3(npl)
  real ( kind = rk ) ww(npl)
  real ( kind = rk ) x(npl)
  real ( kind = rk ) xa(npl)
  real ( kind = rk ) xx(npl)

  dm = real ( m, kind = rk )
  alfa = -1.0D+00 / dm

  ww(1:npl) = x(1:npl)

  do k = 1, m

    call mltply ( ww, xx, npl, w2 )

    ww(1:npl) = w2(1:npl)

  end do

  ww(1) = ww(1) - 2.0D+00

  xa(1) = 2.0D+00
  xa(2:npl) = 0.0D+00

  w3(1) = 2.0D+00
  w3(2:npl) = 0.0D+00

  do k = 2, nt

    dkmm = real ( k - 1, kind = rk )
    dkm2 = real ( k - 2, kind = rk )
    coef = ( alfa - dkm2 ) / dkmm

    call mltply ( w3, ww, npl, w2 )

    w3(1:npl) = w2(1:npl) * coef

    xa(1:npl) = xa(1:npl) + w3(1:npl)

  end do

  call mltply ( xa, xx, npl, w2 )

  xa(1:npl) = w2(1:npl)

  return
end
subroutine cheby ( nf, npl, functn, x )

!*****************************************************************************80
!
!! CHEBY carries out the Chebyshev analysis of one or more functions.
!
!  Discussion:
!
!    This routine carries out the simultaneous Chebyshev analysis of 
!    NF functions.
!
!    The output is a matrix containing one Chebyshev series per column.
!
!    An example of a routine to compute the function values is:
!
!      subroutine functn ( a, val )
!      real ( kind = rk ) a, val(2)
!      val(1) = sin(a)
!      val(2) = cos(a)
!      return
!      end
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Roger Broucke.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Roger Broucke,
!    Algorithm 446:
!    Ten Subroutines for the Manipulation of Chebyshev Series,
!    Communications of the ACM,
!    October 1973, Volume 16, Number 4, pages 254-256.
!
!  Parameters:
!
!    Input, integer NF, the number of functions to be analyzed.
!
!    Input, integer NPL, the number of terms in the 
!    Chebyshev series.
!
!    Input, integer NPLMAX, the leading dimension of X.
!
!    Input, external FUNCTN, the name of a routine which computes
!    the function values at any given point.
!
!    Output, real ( kind = rk ) X(NPL,NF), the Chebyshev series.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nf
  integer npl

  real ( kind = rk ) enn
  real ( kind = rk ) fk
  real ( kind = rk ) fxj(nf)
  real ( kind = rk ) gc(2*(npl-1))
  integer j
  integer k
  integer l
  integer lm
  integer n
  real ( kind = rk ) pen
  real ( kind = rk ) x(npl,nf)
  real ( kind = rk ) xj

  x(1:npl,1:nf) = 0.0D+00

  n = npl - 1
  enn = real ( n, kind = rk )
  pen = 3.1415926535897932D+00 / enn

  do k = 1, 2 * n
    fk = real ( k - 1, kind = rk )
    gc(k) = cos ( fk * pen )
  end do

  do j = 1, npl

    xj = gc(j)
    call functn ( xj, fxj )

    if ( j == 1 .or. j == npl ) then
      fxj(1:nf) = 0.5D+00 * fxj(1:nf)
    end if

    do l = 1, npl
      lm = mod ( ( l - 1 ) * ( j - 1 ), 2 * n ) + 1
      x(l,1:nf) = x(l,1:nf) + fxj(1:nf) * gc(lm)
    end do

  end do

  x(1:npl,1:nf) = 2.0D+00 * x(1:npl,1:nf) / enn

  return
end
subroutine dfrnt ( xx, npl, x2 )

!*****************************************************************************80
!
!! DFRNT determines the derivative of a Chebyshev series.
!
!  Discussion:
!
!    This routine computes the Chebyshev series of the derivative of a 
!    function whose Chebyshev series is given.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Roger Broucke.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Roger Broucke,
!    Algorithm 446:
!    Ten Subroutines for the Manipulation of Chebyshev Series,
!    Communications of the ACM,
!    October 1973, Volume 16, Number 4, pages 254-256.
!
!  Parameters:
!
!    Input, real ( kind = rk ) XX(NPL), the given Chebyshev series.
!
!    Input, integer NPL, the number of terms in the 
!    Chebyshev series.
!
!    Output, real ( kind = rk ) X2(NPL), the Chebyshev series for the
!    derivative.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer npl

  real ( kind = rk ) dl
  real ( kind = rk ) dn
  integer k
  integer l
  real ( kind = rk ) x2(npl)
  real ( kind = rk ) xx(npl)
  real ( kind = rk ) xxl
  real ( kind = rk ) xxn

  dn = real ( npl - 1, kind = rk )
  xxn = xx(npl-1)
  x2(npl-1) = 2.0D+00 * xx(npl) * dn
  x2(npl) = 0.0D+00

  do k = 3, npl
    l = npl - k + 1
    dl = real ( l, kind = rk )
    xxl = xx(l)
    x2(l) = x2(l+2) + 2.0D+00 * xxn * dl
    xxn = xxl
  end do

  return
end
subroutine echeb ( x, coef, npl, fx )

!*****************************************************************************80
!
!! ECHEB evaluates a Chebyshev series at a point.
!
!  Discussion:
!
!  Ê This routine evaluates a Chebyshev series at a point in [-1,+1].
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Roger Broucke.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Roger Broucke,
!    Algorithm 446:
!    Ten Subroutines for the Manipulation of Chebyshev Series,
!    Communications of the ACM,
!    October 1973, Volume 16, Number 4, pages 254-256.
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, the evaluation point.
!    -1 <= X <= +1.
!
!    Input, real ( kind = rk ) COEF(NPL), the Chebyshev series.
!
!    Input, integer NPL, the number of terms in the 
!    Chebyshev series.
!
!    Output, real ( kind = rk ) FX, the value of the Chebyshev series at X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer npl

  real ( kind = rk ) br
  real ( kind = rk ) brp2
  real ( kind = rk ) brpp
  real ( kind = rk ) coef(npl)
  real ( kind = rk ) fx
  integer j
  integer k
  real ( kind = rk ) x

  br = 0.0D+00
  brpp = 0.0D+00

  do k = 1, npl
    j = npl - k + 1
    brp2 = brpp
    brpp = br
    br = 2.0D+00 * x * brpp - brp2 + coef(j)
  end do

  fx = 0.5D+00 * ( br - brp2 )

  return
end
subroutine edcheb ( x, coef, npl, fx )

!*****************************************************************************80
!
!! EDCHEB evaluates the derivative of a Chebyshev series at a point.
!
!  Discussion:
!
!    This routine evaluates the derivative of a Chebyshev series 
!    at a point in [-1,+1].
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Roger Broucke.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Roger Broucke,
!    Algorithm 446:
!    Ten Subroutines for the Manipulation of Chebyshev Series,
!    Communications of the ACM,
!    October 1973, Volume 16, Number 4, pages 254-256.
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, the evaluation point.
!    -1 <= X <= +1.
!
!    Input, real ( kind = rk ) COEF(NPL), the Chebyshev series.
!
!    Input, integer NPL, the number of terms in the 
!    Chebyshev series.
!
!    Output, real ( kind = rk ) FX, the value of the derivative of the
!    Chebyshev series at X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer npl

  real ( kind = rk ) bf
  real ( kind = rk ) bj
  real ( kind = rk ) bjp2
  real ( kind = rk ) bjpl
  real ( kind = rk ) coef(npl)
  real ( kind = rk ) dj
  real ( kind = rk ) fx
  integer j
  integer k
  integer n
  real ( kind = rk ) x
  real ( kind = rk ) xj
  real ( kind = rk ) xjp2
  real ( kind = rk ) xjpl

  xjp2 = 0.0D+00
  xjpl = 0.0D+00
  bjp2 = 0.0D+00
  bjpl = 0.0D+00
  n = npl - 1

  do k = 1, n
    j = npl - k
    dj = real ( j, kind = rk )
    xj = 2.0D+00 * coef(j+1) * dj + xjp2
    bj = 2.0D+00 * x * bjpl - bjp2 + xj
    bf = bjp2
    bjp2 = bjpl
    bjpl = bj
    xjp2 = xjpl
    xjpl = xj
  end do

  fx = 0.5D+00 * ( bj - bf )

  return
end
subroutine invert ( x, xx, npl, net, xnvse )

!*****************************************************************************80
!
!! INVERT computes the inverse Chebyshev series.
!
!  Discussion:
!
!    This routine computes the inverse of a Chebyshev series, starting with
!    an initial approximation XX. 
!
!    The routine uses the Euler method and computes all powers EPS^K 
!    up to K=2^(NET+1), where EPS = 1 - X * ( XX inverse ).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Roger Broucke.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Roger Broucke,
!    Algorithm 446:
!    Ten Subroutines for the Manipulation of Chebyshev Series,
!    Communications of the ACM,
!    October 1973, Volume 16, Number 4, pages 254-256.
!
!  Parameters:
!
!    Input, real ( kind = rk ) X(NPL), the Chebyshev series.
!
!    Input, real ( kind = rk ) XX(NPL), an initial approximation for the
!    inverse Chebyshev series.
!
!    Input, integer NPL, the number of terms in the 
!    Chebyshev series.
!
!    Input, integer NET, the number of iterations to take.
!
!    Output, real ( kind = rk ) XNVSE(NPL), the estimated Chebyshev
!    series of the inverse function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer npl

  integer k
  integer net
  real ( kind = rk ) w2(npl)
  real ( kind = rk ) ww(npl)
  real ( kind = rk ) x(npl)
  real ( kind = rk ) xnvse(npl)
  real ( kind = rk ) xx(npl)

  call mltply ( x, xx, npl, ww )

  ww(1) = 2.0D+00 - ww(1)

  ww(2:npl) = - ww(2:npl)

  call mltply ( ww, ww, npl, w2 )

  ww(1) = 2.0D+00 * ww(1)

  do k = 1, net

    call mltply ( ww, w2, npl, xnvse )

    ww(1:npl) = ww(1:npl) + xnvse(1:npl)

    call mltply ( w2, w2, npl, xnvse )

    w2(1:npl) = xnvse(1:npl)

  end do

  call mltply ( ww, xx, npl, xnvse )

  return
end
subroutine mltply ( xx, x2, npl, x3 )

!*****************************************************************************80
!
!! MLTPLY multiplies two Chebyshev series.
!
!  Discussion:
!
!    This routine multiplies two given Chebyshev series, XX and X2,
!    to produce an output Chebyshev series, X3.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Roger Broucke.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Roger Broucke,
!    Algorithm 446:
!    Ten Subroutines for the Manipulation of Chebyshev Series,
!    Communications of the ACM,
!    October 1973, Volume 16, Number 4, pages 254-256.
!
!  Parameters:
!
!    Input, real ( kind = rk ) XX(NPL), the first Chebyshev series.
!
!    Input, real ( kind = rk ) X2(NPL), the second Chebyshev series.
!
!    Input, integer NPL, the number of terms in the 
!    Chebyshev series.
!
!    Output, real ( kind = rk ) X3(NPL), the Chebyshev series of the
!    product.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer npl

  real ( kind = rk ) ex
  integer k
  integer l
  integer m
  integer mm
  real ( kind = rk ) x2(npl)
  real ( kind = rk ) x3(npl)
  real ( kind = rk ) xx(npl)

  x3(1:npl) = 0.0D+00

  do k = 1, npl
    ex = 0.0D+00
    mm = npl - k + 1
    do m = 1, mm
      l = m + k - 1
      ex = ex + xx(m) * x2(l) + xx(l) * x2(m)
    end do
    x3(k) = 0.5D+00 * ex
  end do

  x3(1) = x3(1) - 0.5D+00 * xx(1) * x2(1)

  do k = 3, npl
    ex = 0.0D+00
    mm = k - 1
    do m = 2, mm
      l = k - m + 1
      ex = ex + xx(m) * x2(l)
    end do
    x3(k) = 0.5D+00 * ex + x3(k)
  end do

  return
end
subroutine ntgrt ( xx, npl, x2 )

!*****************************************************************************80
!
!! NTGRT determines the integral of a Chebyshev series.
!
!  Discussion:
!
!    This routine computes the Chebyshev series for the integral of a 
!    function whose Chebyshev series is given.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Roger Broucke.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Roger Broucke,
!    Algorithm 446:
!    Ten Subroutines for the Manipulation of Chebyshev Series,
!    Communications of the ACM,
!    October 1973, Volume 16, Number 4, pages 254-256.
!
!  Parameters:
!
!    Input, real ( kind = rk ) XX(NPL), the Chebyshev series.
!
!    Input, integer NPL, the number of terms in the 
!    Chebyshev series.
!
!    Output, real ( kind = rk ) X2(NPL), the Chebyshev series for the
!    integral of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer npl

  real ( kind = rk ) dk
  integer k
  integer n
  real ( kind = rk ) term
  real ( kind = rk ) x2(npl)
  real ( kind = rk ) xpr
  real ( kind = rk ) xx(npl)

  xpr = xx(1)
  x2(1) = 0.0D+00
  n = npl - 1

  do k = 2, n
    dk = real ( k - 1, kind = rk )
    term = ( xpr - xx(k+1) ) / ( 2.0D+00 * dk )
    xpr = xx(k)
    x2(k) = term
  end do

  dk = real ( n, kind = rk )
  x2(npl) = xpr / ( 2.0D+00 * dk )

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
subroutine xalfa2 ( x, xx, npl, m, maxet, epsln, net )

!*****************************************************************************80
!
!! XALFA2 computes a Chebyshev series raised to the (-1/M) power.
!
!  Discussion:
!
!    This routine estimates the Chebyshev series for a function raised
!    to the (-1/M) power, given the Chebyshev series for the function,
!    and a starting estimate for the desired series.
!
!    The convergence is quadratic.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Roger Broucke.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Roger Broucke,
!    Algorithm 446:
!    Ten Subroutines for the Manipulation of Chebyshev Series,
!    Communications of the ACM,
!    October 1973, Volume 16, Number 4, pages 254-256.
!
!  Parameters:
!
!    Input, real ( kind = rk ) X(NPL), the Chebyshev series of the function.
!
!    Input/output, real ( kind = rk ) XX(NPL).  On input, an initial
!    approximation to the Chebyshev series of the function raised to the
!    (-1/M) power.  On output, an improved approximation.
!
!    Input, integer NPL, the number of terms in the 
!    Chebyshev series.
!
!    Input, integer M, determines the exponent (-1/M).
!
!    Input, integer MAXET, the maximum number of iterations.
!
!    Input, real ( kind = rk ) EPSLN, the required precision.
!
!    Output, integer NET, the actual number of iterations.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer npl

  real ( kind = rk ) dalfa
  real ( kind = rk ) dm
  real ( kind = rk ) epsln
  integer jx
  integer k
  integer m
  integer maxet
  integer net
  real ( kind = rk ) s
  real ( kind = rk ) tdmm
  real ( kind = rk ) w2(npl)
  real ( kind = rk ) ww(npl)
  real ( kind = rk ) x(npl)
  real ( kind = rk ) xx(npl)

  dm = real ( m, kind = rk )
  dalfa = 1.0D+00 / dm
  tdmm = 2.0D+00 * ( dm + 1.0D+00 )

  do jx = 1, maxet

    ww(1:npl) = x(1:npl)
 
    do  k = 1, m

      call mltply ( ww, xx, npl, w2 )

      ww(1:npl) = w2(1:npl)

    end do

    s = - 2.0D+00 + sum ( abs ( ww(1:npl) ) )

    ww(1:npl) = - ww(1:npl)

    ww(1) = ww(1) + tdmm
    call mltply ( ww, xx, npl, w2 )

    xx(1:npl) = w2(1:npl) * dalfa

    net = jx

    if ( abs ( s ) < epsln ) then
      return
    end if

  end do

  return
end
subroutine xalfa3 ( x, xx, npl, m, maxet, epsln, net )

!*****************************************************************************80
!
!! XALFA3 computes a Chebyshev series raised to the (-1/M) power.
!
!  Discussion:
!
!    This routine estimates the Chebyshev series for a function raised
!    to the (-1/M) power, given the Chebyshev series for the function,
!    and a starting estimate for the desired series.
!
!    The convergence is of order three.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 September 2011
!
!  Author:
!
!    Original FORTRAN77 version by Roger Broucke.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Roger Broucke,
!    Algorithm 446:
!    Ten Subroutines for the Manipulation of Chebyshev Series,
!    Communications of the ACM,
!    October 1973, Volume 16, Number 4, pages 254-256.
!
!  Parameters:
!
!    Input, real ( kind = rk ) X(NPL), the Chebyshev series of the function.
!
!    Input/output, real ( kind = rk ) XX(NPL).  On input, an initial
!    approximation to the Chebyshev series of the function raised to the
!    (-1/M) power.  On output, an improved approximation.
!
!    Input, integer NPL, the number of terms in the 
!    Chebyshev series.
!
!    Input, integer M, determines the exponent (-1/M).
!
!    Input, integer MAXET, the maximum number of iterations.
!
!    Input, real ( kind = rk ) EPSLN, the required precision.
!
!    Output, integer NET, the actual number of iterations.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer npl

  real ( kind = rk ) dalfa
  real ( kind = rk ) dm
  real ( kind = rk ) epsln
  integer jx
  integer k
  integer m
  integer maxet
  integer net
  real ( kind = rk ) p5dml
  real ( kind = rk ) s
  real ( kind = rk ) tdmm
  real ( kind = rk ) w2(npl)
  real ( kind = rk ) ww(npl)
  real ( kind = rk ) x(npl)
  real ( kind = rk ) xx(npl)

  dm = real ( m, kind = rk )
  dalfa = 1.0D+00 / dm
  tdmm = 2.0D+00 * ( dm + 1.0D+00 )
  p5dml = 0.5D+00 * ( dm + 1.0D+00 )

  do jx = 1, maxet

    ww(1:npl) = x(1:npl)

    do k = 1, m
      call mltply ( ww, xx, npl, w2 )
      ww(1:npl) = w2(1:npl)
    end do

    s = - 2.0D+00 + sum ( abs ( ww(1:npl) ) )

    ww(1) = ww(1) - 2.0D+00
    ww(1:npl) = ww(1:npl) * dalfa

    call mltply ( ww, ww, npl, w2 )

    ww(1:npl) = - ww(1:npl)

    w2(1:npl) = w2(1:npl) * p5dml

    ww(1) = ww(1) + 2.0D+00

    w2(1:npl) = w2(1:npl) + ww(1:npl)

    call mltply ( w2, xx, npl, ww )

    xx(1:npl) = ww(1:npl)

    net = jx

    if ( abs ( s ) < epsln ) then
      return
    end if

  end do

  return
end
