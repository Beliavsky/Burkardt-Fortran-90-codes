subroutine chebyshev_coef_1d ( nd, xd, yd, c, xmin, xmax )

!*****************************************************************************80
!
!! CHEBYSHEV_COEF_1D determines the Chebyshev interpolant coefficients.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 September 2012
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
!    Input, real ( kind = rk ) XD(ND), the data locations.
!
!    Input, real ( kind = rk ) YD(ND), the data values.
!
!    Output, real ( kind = rk ) C(ND), the Chebyshev coefficients.
!
!    Output, real ( kind = rk ) XMIN, XMAX, the interpolation interval.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nd

  real ( kind = rk ) a(nd,nd)
  real ( kind = rk ) c(nd)
  integer i
  integer j
  real ( kind = rk ) x(nd)
  real ( kind = rk ) xd(nd)
  real ( kind = rk ) xmax
  real ( kind = rk ) xmin
  real ( kind = rk ) yd(nd)

  if ( nd == 1 ) then
    xmin = xd(1)
    xmax = xd(1)
    c(1) = 1.0D+00
    return
  end if

  xmin = minval ( xd(1:nd) )
  xmax = maxval ( xd(1:nd) )
!
!  Map XD to [-1,+1].
!
  x(1:nd) = ( 2.0D+00 * xd(1:nd) - xmin - xmax ) / ( xmax - xmin )
!
!  Form the Chebyshev Vandermonde matrix.
!
  do j = 1, nd
    do i = 1, nd
      a(i,j) = cos ( acos ( x(i) ) * real ( j - 1, kind = rk ) )
    end do
  end do 
!
!  Solve for the expansion coefficients.
!
  call qr_solve ( nd, nd, a, yd, c )

  return
end
subroutine chebyshev_interp_1d ( nd, xd, yd, ni, xi, yi )

!*****************************************************************************80
!
!! CHEBYSHEV_INTERP_1D determines and evaluates the Chebyshev interpolant.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 September 2012
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
!    Input, real ( kind = rk ) XD(ND), the data locations.
!
!    Input, real ( kind = rk ) YD(ND), the data values.
!
!    Input, integer NI, the number of interpolation points.
!
!    Input, real ( kind = rk ) XI(NI), the interpolation points, which
!    must be each be in the interval [ min(XD), max(XD)].
!
!    Output, real ( kind = rk ) YI(NI), the interpolated values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nd
  integer ni

  real ( kind = rk ) c(nd)
  real ( kind = rk ) xd(nd)
  real ( kind = rk ) xi(ni)
  real ( kind = rk ) xmax
  real ( kind = rk ) xmin
  real ( kind = rk ) yd(nd)
  real ( kind = rk ) yi(ni)

  call chebyshev_coef_1d ( nd, xd, yd, c, xmin, xmax )

  call chebyshev_value_1d ( nd, c, xmin, xmax, ni, xi, yi )

  return
end
subroutine chebyshev_value_1d ( nd, c, xmin, xmax, ni, xi, yi )

!*****************************************************************************80
!
!! CHEBYSHEV_VALUE_1D evaluates a Chebyshev interpolant, given its coefficients.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 September 2012
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
!    Input, real ( kind = rk ) C(ND), the Chebyshev coefficients.
!
!    Input, real ( kind = rk ) XMIN, XMAX, the interpolation interval.
!
!    Input, integer NI, the number of interpolation points.
!
!    Input, real ( kind = rk ) XI(NI), the interpolation points, which
!    must be each be in the interval [XMIN,XMAX].
!
!    Output, real ( kind = rk ) YI(NI), the interpolated values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nd
  integer ni

  real ( kind = rk ) a(ni,nd)
  real ( kind = rk ) c(nd)
  integer i
  integer j
  real ( kind = rk ) x(ni)
  real ( kind = rk ) xi(ni)
  real ( kind = rk ) xmax
  real ( kind = rk ) xmin
  real ( kind = rk ) yi(ni)

  if ( nd == 1 ) then
    yi(1) = c(1)
    return
  end if
!
!  Map XI to [-1,+1].
!
  x(1:ni) = ( 2.0D+00 * xi(1:ni) - xmin - xmax ) / ( xmax - xmin )

  do j = 1, nd
    do i = 1, ni
      a(i,j) = cos ( acos ( x(i) ) * real ( j - 1, kind = rk ) )
    end do
  end do 

  yi(1:ni) = matmul ( a(1:ni,1:nd), c(1:nd) )

  return
end
