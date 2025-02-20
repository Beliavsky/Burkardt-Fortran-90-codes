subroutine triangle_mask ( dim_num, triangle_order, nodes, coord, mask )

!*****************************************************************************80
!
!! triangle_mask() is a user routine which masks triangles.
!
!  Discussion:
!
!    The region to be considered is the [0,4]x[0,4] square.
!
!    We want to remove the lower left triangular corner,
!    and part of the upper right triangular corner.
!
!    The following diagram of the 25 nodes indicates by "O" the
!    nodes that should end up being deleted, although the deletion
!    is actually done by triangles.
!
!    Before masking:
!
!      X - X - X - X - X
!      | \ | \ | \ | \ |
!      X - X - X - X - X
!      | \ | \ | \ | \ |
!      X - X - X - X - X
!      | \ | \ | \ | \ |
!      X - X - X - X - X
!      | \ | \ | \ | \ |
!      X - X - X - X - X
!
!    After masking:
!
!      X - X   O   O   O
!      | \ | \          
!      X - X - X   O   O
!      | \ | \ | \      
!      X - X - X - X - X
!        \ | \ | \ | \ |
!      O   X - X - X - X
!            \ | \ | \ |
!      O   O   X - X - X
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 January 2007
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
!  Remove the lower left corner
!
  if ( centroid(1) + centroid(2) < 2.0D+00 ) then

    mask = .true.
!
!  Remove the upper right section.
!
  else if ( 5.0D+00 < centroid(1) + centroid(2) .and. &
            2.0D+00 < centroid(2) ) then

    mask = .true.
!
!  Keep everything else.
!
  else

    mask = .false.

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

