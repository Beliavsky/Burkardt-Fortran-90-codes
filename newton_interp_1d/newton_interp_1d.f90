subroutine newton_coef_1d ( nd, xd, yd, cd )

!*****************************************************************************80
!
!! newton_coef_1d() computes coefficients of a Newton 1D interpolant.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Input:
!
!    integer ND, the number of data points.
!
!    real ( kind = rk ) XD(ND), the X values at which data was taken.
!
!    real ( kind = rk ) YD(ND), the corresponding Y values.
!
!  Output:
!
!    real ( kind = rk ) CD(ND), the divided difference coefficients.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nd

  real ( kind = rk ) cd(nd)
  integer i
  integer j
  real ( kind = rk ) xd(nd)
  real ( kind = rk ) yd(nd)
!
!  Copy the data values.
!
  cd(1:nd) = yd(1:nd)
!
!  Compute the divided differences.
!
  do i = 2, nd
    do j = nd, i, -1

      cd(j) = ( cd(j) - cd(j-1) ) / ( xd(j) - xd(j+1-i) )

    end do
  end do

  return
end
subroutine newton_value_1d ( nd, xd, cd, ni, xi, yi )

!*****************************************************************************80
!
!! newton_value_1d() evaluates a Newton 1D interpolant.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Input:
!
!    integer ND, the order of the difference table.
!
!    real ( kind = rk ) XD(ND), the X values of the difference table.
!
!    real ( kind = rk ) CD(ND), the divided differences.
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

  real ( kind = rk ) cd(nd)
  integer i
  real ( kind = rk ) xd(nd)
  real ( kind = rk ) xi(ni)
  real ( kind = rk ) yi(ni)

  yi(1:ni) = cd(nd)
  do i = nd - 1, 1, -1
    yi(1:ni) = cd(i) + ( xi(1:ni) - xd(i) ) * yi(1:ni)
  end do

  return
end

