subroutine closest_pair_brute ( n, xy, d, i, j )

!*****************************************************************************80
!
!! closest_pair_brute() finds closest pair of points in 2d using brute force.
!
!  License:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 April 2024
!
!  Author:
!
!    Original MATLAB code by Cleve Moler.
!    This version by John Burkardt.
!
!  Reference:
!
!    Cleve Moler,
!    https://blogs.mathworks.com/cleve/2024/03/28/closest-pair-of-points-problem/
!    Closest pair of points problem,
!    28 March 2024.
!
!  Input:
!
!    integer n: the number of points.
!
!    real ( kind = rk8 ) xy(n,2): the X and Y coordinates of the points.
!
!  Output:
!
!    real ( kind = rk8 ) d: the distance of the closest pair.
!
!    integer i, j: the indices of the closest pair.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  real ( kind = rk8 ) d
  real ( kind = rk8 ) d2
  integer i
  integer i2 
  integer j
  integer j2
  real ( kind = rk8 ) xy(n,2)

  d = huge ( d )
  i = 0
  j = 0

  do i2 = 1, n
    do j2 = i2 + 1, n
      d2 = sqrt ( sum ( ( xy(i2,:) - xy(j2,:) )**2 ) )
      if ( d2 < d ) then
        d = d2
        i = i2
        j = j2
      end if
    end do
  end do

  return
end

