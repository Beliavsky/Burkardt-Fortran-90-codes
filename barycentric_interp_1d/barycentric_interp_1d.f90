subroutine lagcheby1_interp_1d ( nd, xd, yd, ni, xi, yi )

!*****************************************************************************80
!
!! lagcheby1_interp_1d() evaluates the Lagrange Chebyshev 1 interpolant.
!
!  Discussion:
!
!    The weight vector WD computed below is only valid if the data points
!    XD are, as expected, the Chebyshev Type 1 points for [-1,+1], or a linearly
!    mapped version for [A,B].  The XD values may be computed by:
!
!      xd = r8vec_cheby1space ( nd, a, b )
!
!    for instance.
!
!    Thanks to John Ferrier for point out that DENOM and NUMER were not
!    being initialized, 16 September 2013.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jean-Paul Berrut, Lloyd Trefethen,
!    Barycentric Lagrange Interpolation,
!    SIAM Review,
!    Volume 46, Number 3, September 2004, pages 501-517.
!
!  Input:
!
!    integer ND, the number of data points.
!    ND must be at least 1.
!
!    real ( kind = rk ) XD(ND), the data points.
!
!    real ( kind = rk ) YD(ND), the data values.
!
!    integer NI, the number of interpolation points.
!
!    real ( kind = rk ) XI(NI), the interpolation points.
!
!  Output:
!
!    real ( kind = rk ) YI(NI), the interpolated values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nd
  integer ni

  real ( kind = rk ) denom(ni)
  integer exact(ni)
  integer i
  integer j
  real ( kind = rk ) numer(ni)
  real ( kind = rk ) r8_mop
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) t
  real ( kind = rk ) theta
  real ( kind = rk ) wd
  real ( kind = rk ) xd(nd)
  real ( kind = rk ) xi(ni)
  real ( kind = rk ) yd(nd)
  real ( kind = rk ) yi(ni)

  exact(1:ni) = 0
  numer(1:ni) = 0.0D+00
  denom(1:ni) = 0.0D+00

  do j = 1, nd

    theta = real ( 2 * j - 1, kind = rk ) * r8_pi / real ( 2 * nd, kind = rk )
    wd = r8_mop ( j + 1 ) * sin ( theta )

    do i = 1, ni

      if ( xi(i) == xd(j) ) then
        exact(i) = j
        numer(i) = yd(j)
        denom(i) = 1.0D+00
      end if

      if ( exact(i) == 0 ) then
        t = wd / ( xi(i) - xd(j) )
        numer(i) = numer(i) + t * yd(j)
        denom(i) = denom(i) + t
      end if

    end do
  end do

  yi(1:ni) = numer(1:ni) / denom(1:ni)

  return
end
subroutine lagcheby2_interp_1d ( nd, xd, yd, ni, xi, yi )

!*****************************************************************************80
!
!! lagcheby2_interp_1d() evaluates the Lagrange Chebyshev 2 interpolant.
!
!  Discussion:
!
!    The weight vector WD computed below is only valid if the data points
!    XD are, as expected, the Chebyshev Type 2 points for [-1,+1], or a
!    linearly mapped version for [A,B].  The XD values may be computed by:
!
!      xd = r8vec_cheby2space ( nd, a, b )
!
!    for instance.
!
!    Thanks to John Ferrier for point out that DENOM and NUMER were not
!    being initialized, 16 September 2013.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jean-Paul Berrut, Lloyd Trefethen,
!    Barycentric Lagrange Interpolation,
!    SIAM Review,
!    Volume 46, Number 3, September 2004, pages 501-517.
!
!  Input:
!
!    integer ND, the number of data points.
!    ND must be at least 1.
!
!    real ( kind = rk ) XD(ND), the data points.
!
!    real ( kind = rk ) YD(ND), the data values.
!
!    integer NI, the number of interpolation points.
!
!    real ( kind = rk ) XI(NI), the interpolation points.
!
!  Output:
!
!    real ( kind = rk ) YI(NI), the interpolated values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nd
  integer ni

  real ( kind = rk ) denom(ni)
  integer exact(ni)
  integer i
  integer j
  real ( kind = rk ) numer(ni)
  real ( kind = rk ) r8_mop
  real ( kind = rk ) t
  real ( kind = rk ) wd
  real ( kind = rk ) xd(nd)
  real ( kind = rk ) xi(ni)
  real ( kind = rk ) yd(nd)
  real ( kind = rk ) yi(ni)

  exact(1:ni) = 0
  numer(1:ni) = 0.0D+00
  denom(1:ni) = 0.0D+00

  do j = 1, nd

    wd = r8_mop ( j + 1 )
    if ( j == 1 .or. j == nd ) then
      wd = 0.5D+00 * wd
    end if

    do i = 1, ni

      if ( xi(i) == xd(j) ) then
        exact(i) = j
        numer(i) = yd(j)
        denom(i) = 1.0D+00
      end if

      if ( exact(i) == 0 ) then
        t = wd / ( xi(i) - xd(j) )
        numer(i) = numer(i) + t * yd(j)
        denom(i) = denom(i) + t
      end if

    end do
  end do

  yi(1:ni) = numer(1:ni) / denom(1:ni)

  return
end
subroutine lageven_interp_1d ( nd, xd, yd, ni, xi, yi )

!*****************************************************************************80
!
!! lageven_value_1d() evaluates the Lagrange evenly-spaced interpolant.
!
!  Discussion:
!
!    The weight vector WD computed below is only valid if the data points
!    XD are, as expected, evenly spaced in an interval [A,B] with
!    spacing (B-A)/N.  The XD values might be computed by:
!
!      xd(i) = ( ( 2 * nd - 2 * i + 1 ) * a 
!              + (          2 * i - 1 ) * b ) 
!              / ( 2 * nd             )
!
!    for instance.
!
!    Thanks to John Ferrier for point out that DENOM and NUMER were not
!    being initialized, 16 September 2013.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jean-Paul Berrut, Lloyd Trefethen,
!    Barycentric Lagrange Interpolation,
!    SIAM Review,
!    Volume 46, Number 3, September 2004, pages 501-517.
!
!  Input:
!
!    integer ND, the number of data points.
!    ND must be at least 1.
!
!    real ( kind = rk ) XD(ND), the data points.
!
!    real ( kind = rk ) YD(ND), the data values.
!
!    integer NI, the number of interpolation points.
!
!    real ( kind = rk ) XI(NI), the interpolation points.
!
!  Output:
!
!    real ( kind = rk ) YI(NI), the interpolated values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nd
  integer ni

  real ( kind = rk ) denom(ni)
  integer exact(ni)
  integer i
  integer j
  real ( kind = rk ) numer(ni)
  real ( kind = rk ) r8_choose
  real ( kind = rk ) r8_mop
  real ( kind = rk ) t
  real ( kind = rk ) wd
  real ( kind = rk ) xd(nd)
  real ( kind = rk ) xi(ni)
  real ( kind = rk ) yd(nd)
  real ( kind = rk ) yi(ni)

  exact(1:ni) = 0
  numer(1:ni) = 0.0D+00
  denom(1:ni) = 0.0D+00

  do j = 1, nd

    wd = r8_mop ( j ) * r8_choose ( nd, j )

    do i = 1, ni

      if ( xi(i) == xd(j) ) then
        exact(i) = j
        numer(i) = yd(j)
        denom(i) = 1.0D+00
      end if

      if ( exact(i) == 0 ) then
        t = wd / ( xi(i) - xd(j) )
        numer(i) = numer(i) + t * yd(j)
        denom(i) = denom(i) + t
      end if

    end do
  end do

  yi(1:ni) = numer(1:ni) / denom(1:ni)

  return
end
 
