subroutine lagrange_basis_function_1d ( mx, xd, i, xi, yi ) 

!*****************************************************************************80
!
!! LAGRANGE_BASIS_FUNCTION_1D evaluates one 1D Lagrange basis function.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MX, the degree of the basis function.
!
!    Input, real ( kind = rk ) XD(MX+1), the interpolation nodes.
!
!    Input, integer ( kind = 4 ) I, the index of the basis function.
!    1 <= I <= MX+1.
!
!    Input, real ( kind = rk ) XI, the evaluation point.
!
!    Output, real ( kind = rk ) YI, the value of the I-th Lagrange 1D basis 
!    function for the nodes XD, evaluated at XI.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ( kind = 4 ) mx

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = rk ) xd(mx+1)
  real ( kind = rk ) xi
  real ( kind = rk ) yi

  yi = 1.0D+00

  if ( xi /= xd(i) ) then
    do j = 1, mx + 1
      if ( j /= i ) then
        yi = yi * ( xi - xd(j) ) / ( xd(i) - xd(j) )
      end if
    end do
  end if

  return
end
subroutine lagrange_interp_2d ( mx, my, xd_1d, yd_1d, zd, ni, xi, yi, zi )

!*****************************************************************************80
!
!! LAGRANGE_INTERP_2D evaluates the Lagrange interpolant for a product grid.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MX, MY, the polynomial degree in X and Y.
!
!    Input, real ( kind = rk ) XD_1D(MX+1), YD_1D(MY+1), the 1D data locations.
!
!    Input, real ( kind = rk ) ZD((MX+1),(MY+1)), the 2D array of data values.
!
!    Input, integer ( kind = 4 ) NI, the number of 2D interpolation points.
!
!    Input, real ( kind = rk ) XI(NI), YI(NI), the 2D interpolation points.
!
!    Output, real ( kind = rk ) ZI(NI), the interpolated values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ( kind = 4 ) mx
  integer ( kind = 4 ) my
  integer ( kind = 4 ) ni

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = rk ) lx
  real ( kind = rk ) ly
  real ( kind = rk ) xd_1d(mx+1)
  real ( kind = rk ) xi(ni)
  real ( kind = rk ) yd_1d(my+1)
  real ( kind = rk ) yi(ni)
  real ( kind = rk ) zd(mx+1,my+1)
  real ( kind = rk ) zi(ni)

  do k = 1, ni
    l = 0
    zi(k) = 0.0D+00
    do i = 1, mx + 1
      do j = 1, my + 1
        l = l + 1
        call lagrange_basis_function_1d ( mx, xd_1d, i, xi(k), lx )
        call lagrange_basis_function_1d ( my, yd_1d, j, yi(k), ly )
        zi(k) = zi(k) + zd(i,j) * lx * ly
      end do
    end do
  end do

  return
end
