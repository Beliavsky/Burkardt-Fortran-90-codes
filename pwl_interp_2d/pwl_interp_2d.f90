subroutine pwl_interp_2d ( nxd, nyd, xd, yd, zd, ni, xi, yi, zi )

!*****************************************************************************80
!
!! PWL_INTERP_2D: piecewise linear interpolant to data defined on a 2D grid.
!
!  Discussion:
!
!    Thanks to Adam Hirst for pointing out an error in the formula that
!    chooses the interpolation triangle, 04 February 2018.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 February 2018
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NXD, NYD, the number of X and Y data values.
!
!    Input, real ( kind = rk ) XD(NXD), YD(NYD), the sorted X and Y data.
!
!    Input, real ( kind = rk ) ZD(NXD,NYD), the Z data.
!
!    Input, integer NI, the number of interpolation points.
!
!    Input, real ( kind = rk ) XI(NI), YI(NI), the coordinates of the
!    interpolation points.
!
!    Output, real ( kind = rk ) ZI(NI), the value of the interpolant.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ni
  integer nxd
  integer nyd

  real ( kind = rk ) alpha
  real ( kind = rk ) beta
  real ( kind = rk ) det
  real ( kind = rk ) dxa
  real ( kind = rk ) dxb
  real ( kind = rk ) dxi
  real ( kind = rk ) dya
  real ( kind = rk ) dyb
  real ( kind = rk ) dyi
  real ( kind = rk ) gamma
  integer i
  integer j
  integer k
  real ( kind = rk ) r8_huge
  integer r8vec_bracket5
  real ( kind = rk ) xd(nxd)
  real ( kind = rk ) xi(ni)
  real ( kind = rk ) yd(nyd)
  real ( kind = rk ) yi(ni)
  real ( kind = rk ) zd(nxd,nyd)
  real ( kind = rk ) zi(ni)

  do k = 1, ni
!
!  For interpolation point (xi(k),yi(k)), find data intervals I and J so that:
!
!    xd(i) <= xi(k) <= xd(i+1),
!    yd(j) <= yi(k) <= yd(j+1).
!
!  But if the interpolation point is not within a data interval, 
!  assign the dummy interpolant value zi(k) = infinity.
!
    i = r8vec_bracket5 ( nxd, xd, xi(k) )
    if ( i == -1 ) then
      zi(k) = r8_huge ( )
      cycle
    end if

    j = r8vec_bracket5 ( nyd, yd, yi(k) )
    if ( j == -1 ) then
      zi(k) = r8_huge ( )
      cycle
    end if
!
!  The rectangular cell is arbitrarily split into two triangles.
!  The linear interpolation formula depends on which triangle 
!  contains the data point.
!
!    (I,J+1)--(I+1,J+1)
!      |\       |
!      | \      |
!      |  \     |
!      |   \    |
!      |    \   |
!      |     \  |
!    (I,J)---(I+1,J)
!
    if ( yi(k) < yd(j+1) &
      + ( yd(j) - yd(j+1) ) * ( xi(k) - xd(i) ) / ( xd(i+1) - xd(i) ) ) then

      dxa = xd(i+1) - xd(i)
      dya = yd(j)   - yd(j)

      dxb = xd(i)   - xd(i)
      dyb = yd(j+1) - yd(j)

      dxi = xi(k)   - xd(i)
      dyi = yi(k)   - yd(j)

      det = dxa * dyb - dya * dxb

      alpha = ( dxi * dyb - dyi * dxb ) / det
      beta =  ( dxa * dyi - dya * dxi ) / det
      gamma = 1.0D+00 - alpha - beta

      zi(k) = alpha * zd(i+1,j) + beta * zd(i,j+1) + gamma * zd(i,j)

    else

      dxa = xd(i)   - xd(i+1)
      dya = yd(j+1) - yd(j+1)

      dxb = xd(i+1) - xd(i+1)
      dyb = yd(j)   - yd(j+1)

      dxi = xi(k)   - xd(i+1)
      dyi = yi(k)   - yd(j+1)

      det = dxa * dyb - dya * dxb

      alpha = ( dxi * dyb - dyi * dxb ) / det
      beta =  ( dxa * dyi - dya * dxi ) / det
      gamma = 1.0D+00 - alpha - beta

      zi(k) = alpha * zd(i,j+1) + beta * zd(i+1,j) + gamma * zd(i+1,j+1)

    end if

  end do

  return
end
