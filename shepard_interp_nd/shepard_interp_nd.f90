subroutine shepard_interp_nd ( m, nd, xd, zd, p, ni, xi, zi )

!*****************************************************************************80
!
!! SHEPARD_INTERP_ND evaluates a Shepard interpolant in M dimensions.
!
!  Discussion:
!
!    This code should be vectorized.
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
!  Reference:
!
!    Donald Shepard,
!    A two-dimensional interpolation function for irregularly spaced data,
!    ACM '68: Proceedings of the 1968 23rd ACM National Conference,
!    ACM, pages 517-524, 1969.
!
!  Parameters:
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer ND, the number of data points.
!
!    Input, real ( kind = rk ) XD(M,ND), the data points.
!
!    Input, real ( kind = rk ) ZD(ND), the data values.
!
!    Input, real ( kind = rk ) P, the power.
!
!    Input, integer NI, the number of interpolation points.
!
!    Input, real ( kind = rk ) XI(M,NI), the interpolation points.
!
!    Output, real ( kind = rk ) ZI(NI), the interpolated values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer nd
  integer ni

  integer i
  integer j
  real ( kind = rk ) p
  real ( kind = rk ) r8vec_norm_affine
  real ( kind = rk ) s
  real ( kind = rk ) w(nd)
  real ( kind = rk ) xd(m,nd)
  real ( kind = rk ) xi(m,nd)
  integer z
  real ( kind = rk ) zd(nd)
  real ( kind = rk ) zi(ni)

  do i = 1, ni

    if ( p == 0.0D+00 ) then

      w(1:nd) = 1.0D+00 / real ( nd, kind = rk )

    else

      z = -1
      do j = 1, nd
        w(j) = r8vec_norm_affine ( m, xi(1:m,i), xd(1:m,j) )
        if ( w(j) == 0.0D+00 ) then
          z = j
          exit
        end if
      end do

      if ( z /= -1 ) then
        w(1:nd) = 0.0D+00
        w(z) = 1.0D+00
      else
        w(1:nd) = 1.0D+00 / w(1:nd) ** p
        s = sum ( w )
        w(1:nd) = w(1:nd) / s
      end if

    end if

    zi(i) = dot_product ( w, zd )

  end do

  return
end
