subroutine nearest_interp_1d ( nd, xd, yd, ni, xi, yi )

!*****************************************************************************80
!
!! NEAREST_INTERP_1D evaluates the nearest neighbor interpolant.
!
!  Discussion:
!
!    The nearest neighbor interpolant L(ND,XD,YD)(X) is the piecewise
!    constant function which interpolates the data (XD(I),YD(I)) for I = 1
!    to ND.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ND, the number of data points.
!    ND must be at least 1.
!
!    Input, real ( kind = rk ) XD(ND), the data points.
!
!    Input, real ( kind = rk ) YD(ND), the data values.
!
!    Input, integer NI, the number of interpolation points.
!
!    Input, real ( kind = rk ) XI(NI), the interpolation points.
!
!    Output, real ( kind = rk ) YI(NI), the interpolated values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nd
  integer ni

  real ( kind = rk ) d
  real ( kind = rk ) d2
  integer i
  integer j
  integer k
  real ( kind = rk ) xd(nd)
  real ( kind = rk ) xi(ni)
  real ( kind = rk ) yd(nd)
  real ( kind = rk ) yi(ni)

  do i = 1, ni

    k = 1
    d = abs ( xi(i) - xd(k) )

    do j = 2, nd

      d2 = abs ( xi(i) - xd(j) )

      if ( d2 < d ) then
        k = j
        d = d2
      end if

    end do

    yi(i) = yd(k)

  end do

  return
end
