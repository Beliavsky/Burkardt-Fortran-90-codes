
subroutine closest_point_brute ( m, nr, r, s, near_index, near_dist )

!*****************************************************************************80
!
!! closest_point_brute() finds the nearest data point to a given point.
!
!  Discussion:
!
!    We are given R, a set of NR points in M dimensions.
!
!    We are given S, a point in M dimensions.
!
!    We seek the index J of the point R(J)
!    which is nearest to S over all points in R.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 August 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, the spatial dimension.
!
!    integer NR, the number of data points.
!
!    real ( kind = rk8 ) R(M,NR), the data points.
!
!    real ( kind = rk8 ) S(M), the sample points.
!
!  Output:
!
!    integer near_index, the index of the nearest data point.
!
!    real ( kind = rk8 ) near_dist, the distance of the nearest data point.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer m
  integer nr

  real ( kind = rk8 ) dist_sq
  integer jr
  real ( kind = rk8 ) near_dist
  integer near_index
  real ( kind = rk8 ) r(m,nr)
  real ( kind = rk8 ) s(m)

  near_dist = huge ( near_dist )
  near_index = -1

  do jr = 1, nr

    dist_sq = sum ( ( r(1:m,jr) - s(1:m) ) **2 )

    if ( dist_sq < near_dist ) then
      near_dist = dist_sq
      near_index = jr
    end if

  end do

  near_dist = sqrt ( near_dist )

  return
end

