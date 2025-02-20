subroutine pwl_basis_1d ( nd, xd, ni, xi, bk )

!*****************************************************************************80
!
!! PWL_BASIS_1D evaluates a 1D piecewise linear basis function.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ND, the number of data points.
!
!    Input, real ( kind = rk ) XD(ND), the data points.
!
!    Input, integer NI, the number of interpolation points.
!
!    Input, real ( kind = rk ) XI(NI), the interpolation points.
!
!    Output, real ( kind = rk ) BK(NI,ND), the basis functions at the 
!    interpolation points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nd
  integer ni

  real ( kind = rk ) bk(ni,nd)
  integer i
  integer j
  real ( kind = rk ) t
  real ( kind = rk ) xd(nd)
  real ( kind = rk ) xi(ni)

  bk(1:ni,1:nd) = 0.0D+00

  if ( nd == 1 ) then
    bk(1:ni,1:nd) = 1.0D+00
    return
  end if

  do i = 1, ni

    do j = 1, nd

      if ( j == 1 .and. xi(i) <= xd(j) ) then

        t = ( xi(i) - xd(j) ) / ( xd(j+1) - xd(j) )
        bk(i,j) = 1.0D+00 - t

      else if ( j == nd .and. xd(j) <= xi(i) ) then

        t = ( xi(i) - xd(j-1) ) / ( xd(j) - xd(j-1) )
        bk(i,j) = t

      else if ( xd(j-1) < xi(i) .and. xi(i) <= xd(j) ) then

        t = ( xi(i) - xd(j-1) ) / ( xd(j) - xd(j-1) )
        bk(i,j) = t

      else if ( xd(j) <= xi(i) .and. xi(i) < xd(j+1) ) then

        t = ( xi(i) - xd(j) ) / ( xd(j+1) - xd(j) )
        bk(i,j) = 1.0D+00 - t

      end if

    end do

  end do

  return
end
subroutine pwl_value_1d ( nd, xd, yd, ni, xi, yi )

!*****************************************************************************80
!
!! PWL_VALUE_1D evaluates the piecewise linear interpolant.
!
!  Discussion:
!
!    The piecewise linear interpolant L(ND,XD,YD)(X) is the piecewise
!    linear function which interpolates the data (XD(I),YD(I)) for I = 1
!    to ND.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 September 2012
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

  integer i
  integer k
  real ( kind = rk ) t
  real ( kind = rk ) xd(nd)
  real ( kind = rk ) yd(nd)
  real ( kind = rk ) xi(ni)
  real ( kind = rk ) yi(ni)

  yi(1:ni) = 0.0D+00

  if ( nd == 1 ) then
    yi(1:ni) = yd(1)
    return
  end if

  do i = 1, ni

    if ( xi(i) <= xd(1) ) then

      t = ( xi(i) - xd(1) ) / ( xd(2) - xd(1) )
      yi(i) = ( 1.0D+00 - t ) * yd(1) + t * yd(2)

    else if ( xd(nd) <= xi(i) ) then

      t = ( xi(i) - xd(nd-1) ) / ( xd(nd) - xd(nd-1) )
      yi(i) = ( 1.0D+00 - t ) * yd(nd-1) + t * yd(nd)

    else

      do k = 2, nd

        if ( xd(k-1) <= xi(i) .and. xi(i) <= xd(k) ) then

          t = ( xi(i) - xd(k-1) ) / ( xd(k) - xd(k-1) )
          yi(i) = ( 1.0D+00 - t ) * yd(k-1) + t * yd(k)
          exit

        end if

      end do

    end if

  end do
  
  return
end
