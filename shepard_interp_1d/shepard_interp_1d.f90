subroutine shepard_basis_1d ( nd, xd, p, ni, xi, bk )

!*****************************************************************************80
!
!! SHEPARD_BASIS_1D evaluates a 1D Shepard basis function.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 July 2015
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
!    Input, integer ND, the number of data points.
!
!    Input, real ( kind = rk ) XD(ND), the data points.
!
!    Input, real ( kind = rk ) P, the power.
!
!    Input, integer NI, the number of interpolation points.
!
!    Input, real ( kind = rk ) XI(NI), the interpolation points.
!
!    Output, real ( kind = rk ) BK(NI,ND), the basis function at the 
!    interpolation points.
! 
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nd
  integer ni

  real ( kind = rk ) bk(ni,nd)
  integer i
  integer j
  real ( kind = rk ) p
  real ( kind = rk ) s
  real ( kind = rk ) w(nd)
  real ( kind = rk ) xd(nd)
  real ( kind = rk ) xi(ni)
  integer z

  do i = 1, ni

    if ( p == 0.0D+00 ) then

      w(1:nd) = 1.0D+00 / real ( nd, kind = rk )

    else

      z = -1
      do j = 1, nd
        w(j) = abs ( xi(i) - xd(j) )
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
        s = sum ( w(1:nd) )
        w(1:nd) = w(1:nd) / s
      end if

    end if

    bk(i,1:nd) = w(1:nd)

  end do


  return
end
subroutine shepard_value_1d ( nd, xd, yd, p, ni, xi, yi )

!*****************************************************************************80
!
!! SHEPARD_VALUE_1D evaluates a 1D Shepard interpolant.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 September 2012
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
!    Input, integer ND, the number of data points.
!
!    Input, real ( kind = rk ) XD(ND), the data points.
!
!    Input, real ( kind = rk ) YD(ND), the data values.
!
!    Input, real ( kind = rk ) P, the power.
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
  integer j
  real ( kind = rk ) p
  real ( kind = rk ) s
  real ( kind = rk ) w(nd)
  real ( kind = rk ) xd(nd)
  real ( kind = rk ) xi(ni)
  real ( kind = rk ) yd(nd)
  real ( kind = rk ) yi(ni)
  integer z

  do i = 1, ni

    if ( p == 0.0D+00 ) then

      w(1:nd) = 1.0D+00 / real ( nd, kind = rk )

    else

      z = 0
      do j = 1, nd
        w(j) = abs ( xi(i) - xd(j) )
        if ( w(j) == 0.0D+00 ) then
          z = j
          exit
        end if
      end do

      if ( z /= 0 ) then
        w(1:nd) = 0.0D+00
        w(z) = 1.0D+00
      else
        w(1:nd) = 1.0D+00 / w(1:nd) ** p
        s = sum ( w )
        w(1:nd) = w(1:nd) / s
      end if

    end if

    yi(i) = dot_product ( w(1:nd), yd(1:nd) )

  end do

  return
end
