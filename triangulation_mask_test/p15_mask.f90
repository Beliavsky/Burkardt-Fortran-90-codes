subroutine triangle_mask ( dim_num, triangle_order, nodes, coord, mask )

!*****************************************************************************80
!
!! triangle_mask() is a user routine which masks triangles.
!
!  Discussion:
!
!    The region to be considered is the union of two rectangles.
!    The first is  -8 <= X <= 2, -1 <= Y <= 0,
!    the second is -2 <= X <= 8,  0 <= Y <= 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer TRIANGLE_ORDER, the number of nodes in the triangle.
!
!    Input, integer NODES(TRIANGLE_ORDER), the indices of the nodes.
!
!    Input, real ( kind = rk ) COORD(DIM_NUM,TRIANGLE_ORDER), the coordinates
!    of the nodes.
!
!    Output, logical MASK, is TRUE if the triangle should be discarded,
!    and FALSE if the triangle should be retained.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num
  integer triangle_order

  real ( kind = rk ) centroid(dim_num)
  real ( kind = rk ) coord(dim_num,triangle_order)
  integer dim
  logical mask
  integer nodes(triangle_order)

  call i4_fake_use ( nodes(1) )
!
!  Compute the centroid.
!
  do dim = 1, dim_num
    centroid(dim) = sum ( coord(dim,1:triangle_order) ) &
                  / real ( triangle_order, kind = rk )
  end do
!
!  MASK = The centroid is outside the region.
!
  if (      -8.0D+00 <= centroid(1)            .and. &
                        centroid(1) <= 2.0D+00 .and. &
            -1.0D+00 <= centroid(2)            .and. &
                        centroid(2) <= 0.0D+00 ) then

    mask = .false.

  else if ( -2.0D+00 <= centroid(1)            .and. &
                        centroid(1) <= 8.0D+00 .and. &
             0.0D+00 <= centroid(2)            .and. &
                        centroid(2) <= 1.0D+00 ) then

    mask = .false.

  else

    mask = .true.

  end if

  return
end
subroutine i4_fake_use ( n )

!*****************************************************************************80
!
!! i4_fake_use() pretends to use a variable.
!
!  Discussion:
!
!    Some compilers will issue a warning if a variable is unused.
!    Sometimes there's a good reason to include a variable in a program,
!    but not to use it.  Calling this function with that variable as
!    the argument will shut the compiler up.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 April 2020
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the variable to be "used".
!
  implicit none

  integer n

  if ( n /= n ) then
    write ( *, '(a)' ) '  i4_fake_use(): variable is NAN.'
  end if

  return
end
