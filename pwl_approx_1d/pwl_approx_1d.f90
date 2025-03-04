subroutine pwl_approx_1d ( nd, xd, yd, nc, xc, yc )

!*****************************************************************************80
!
!! PWL_APPROX_1D determines the control values for a PWL approximant.
!
!  Discussion:
!
!    The piecewise linear approximant is defined by NC control pairs 
!    (XC(I),YC(I)) and approximates ND data pairs (XD(I),YD(I)).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 October 2012
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
!    Input, integer NC, the number of control points.
!    NC must be at least 1.
!
!    Input, real ( kind = rk ) XC(NC), the control points.
!
!    Output, real ( kind = rk ) YC(NC), the control values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nc
  integer nd

  real ( kind = rk ) a(nd,nc)
  real ( kind = rk ) xc(nc)
  real ( kind = rk ) xd(nd)
  real ( kind = rk ) yc(nc)
  real ( kind = rk ) yd(nd)
!
!  Define the NDxNC linear system that determines the control values.
!
  call pwl_approx_1d_matrix ( nd, xd, yd, nc, xc, a )
!
!  Solve the system.
!
  call qr_solve ( nd, nc, a, yd, yc )

  return
end
subroutine pwl_approx_1d_matrix ( nd, xd, yd, nc, xc, a )

!*****************************************************************************80
!
!! PWL_APPROX_1D_MATRIX returns the matrix for the PWL approximant controls.
!
!  Discussion:
!
!    The value of the piecewise linear approximant, using control points XC
!    and control values YC, evaluated at the point XD, can be represented by
!
!      YD = A * YC
!
!    where A is a matrix whose values depend on XC and XD.
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
!    Input, integer NC, the number of control points.
!    NC must be at least 1.
!
!    Input, real ( kind = rk ) XC(NC), the control points.
!
!    Output, real ( kind = rk ) A(ND,NC), the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nc
  integer nd

  real ( kind = rk ) a(nd,nc)
  integer i
  integer j
  integer k
  real ( kind = rk ) t
  real ( kind = rk ) xc(nc)
  real ( kind = rk ) xd(nd)
  real ( kind = rk ) yd(nd)

  call r8_fake_use ( yd(1) )

  a(1:nd,1:nc) = 0.0D+00

  do i = 1, nd
    k = nc - 1
    do j = 2, nc - 1
      if ( xd(i) < xc(j) ) then
        k = j - 1
        exit
      end if
    end do
    t = ( xd(i) - xc(k) ) / ( xc(k+1) - xc(k) )
    a(i,k)   = 1.0D+00 - t
    a(i,k+1) =           t
  end do

  return
end
subroutine pwl_interp_1d ( nd, xd, yd, ni, xi, yi )

!*****************************************************************************80
!
!! PWL_INTERP_1D evaluates the piecewise linear interpolant.
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
!    06 September 2012
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

