program main

!*****************************************************************************80
!
!! bezier_surface_test() tests bezier_surface().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'bezier_surface_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test bezier_surface().'

  call test01 ( )
  call test02 ( )
  call test03 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'bezier_surface_test()'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! test01() tests routines to read a Bezier surface definition.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = 100 ) :: node_file_name = 'teapot_nodes.txt'
  integer node_num
  real ( kind = rk ), allocatable, dimension ( :, : ) :: node_xyz
  character ( len = 100 ) :: rectangle_file_name = 'teapot_rectangles.txt'
  integer rectangle_num
  integer, allocatable, dimension ( :, : ) :: rectangle_node

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test01()'
  write ( *, '(a)' ) '  bezier_surface_node_size() determines the number of'
  write ( *, '(a)' ) '  nodes in a Bezier surface node file.'
  write ( *, '(a)' ) '  bezier_surface_node_read() reads the'
  write ( *, '(a)' ) '  nodes in a Bezier surface node file.'
  write ( *, '(a)' ) '  bezier_surface_rectangles_size() determines the number of'
  write ( *, '(a)' ) '  rectangles in a Bezier surface rectangle file.'
  write ( *, '(a)' ) '  bezier_surface_rectangles_read() reads the'
  write ( *, '(a)' ) '  rectangles in a Bezier surface rectangle file.'
!
!  Get the number of nodes, allocate space for them, and read them in.
!
  call bezier_surface_node_size ( node_file_name, node_num )
  allocate ( node_xyz(1:3,1:node_num) )
  call bezier_surface_node_read ( node_file_name, node_num, node_xyz )
  call bezier_surface_node_print ( node_num, node_xyz )
!
!  Get the number of rectangles, allocate space for them, and read them in.
!
  call bezier_surface_rectangle_size ( rectangle_file_name, rectangle_num )
  allocate ( rectangle_node(1:16,1:rectangle_num) )
  call bezier_surface_rectangle_read ( rectangle_file_name, rectangle_num, &
    rectangle_node )
  call bezier_surface_rectangle_print ( rectangle_num, rectangle_node )
!
!  Free the memory.
!
  deallocate ( node_xyz )
  deallocate ( rectangle_node )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! test02() tests bezier_surface_neighbors().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = 100 ) :: rectangle_file_name = 'teapot_rectangles.txt'
  integer rectangle_num
  integer, allocatable, dimension ( :, : ) :: rectangle_neighbor
  integer, allocatable, dimension ( :, : ) :: rectangle_node

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test02():'
  write ( *, '(a)' ) '  bezier_surface_neighbors() determines patch neighbors.'
  write ( *, '(a)' ) '  Note that, for this example, the teapot, there are'
  write ( *, '(a)' ) '  cases where more than two patches meet at a'
  write ( *, '(a)' ) '  (degenerate) side.  This routine will not handle'
  write ( *, '(a)' ) '  such cases completely.'
!
!  Get the number of rectangles, allocate space for them, and read them in.
!
  call bezier_surface_rectangle_size ( rectangle_file_name, rectangle_num )
  allocate ( rectangle_node(1:16,1:rectangle_num) )
  call bezier_surface_rectangle_read ( rectangle_file_name, rectangle_num, &
    rectangle_node )
! call bezier_surface_rectangle_print ( rectangle_num, rectangle_node )
!
!  Compute and print the neighbor array.
!
  allocate ( rectangle_neighbor(1:4,1:rectangle_num) )
  call bezier_surface_neighbors ( rectangle_num, rectangle_node, &
    rectangle_neighbor )
  call i4mat_transpose_print ( 4, rectangle_num, rectangle_neighbor, &
    '  Bezier patch neighbors:' )
!
!  Free the memory.
!
  deallocate ( rectangle_neighbor )
  deallocate ( rectangle_node )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! test03() tests bezier_patch_evaluate().
!
!  Discussion:
!
!    For simplicity, we set up a Bezier surface of a single patch.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: node_num = 16
  integer, parameter :: point_num = 16
  integer, parameter :: rectangle_num = 1

  integer i
  integer j
  integer node
  real ( kind = rk ) node_xyz(3,node_num)
  integer patch
  integer point
  real ( kind = rk ) point_uv(2,point_num)
  real ( kind = rk ) point_xyz(3,point_num)
  integer, dimension(16,rectangle_num) :: rectangle_node = reshape ( (/ &
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16 /), (/ 16, 1 /) )
  real ( kind = rk ) x
  real ( kind = rk ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test03():'
  write ( *, '(a)' ) '  bezier_patch_evaluate() evaluates points in one'
  write ( *, '(a)' ) '  patch of a Bezier surface.'

  node = 0
  do j = 1, 4
    y = real ( j - 1, kind = rk ) / 3.0D+00
    do i = 1, 4
      x = real ( i - 1, kind = rk ) / 3.0D+00
      node = node + 1
      node_xyz(1,node) = x
      node_xyz(2,node) = y
      node_xyz(3,node) = x * ( 1.0D+00 - x ) * y * y
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Nodal coordinates:'
  write ( *, '(a)' ) ' '

  do node = 1, node_num
    write ( *, '(2x,i2,2x,5g12.4)' ) node, node_xyz(1:3,node)
  end do

  patch = 1

  point = 0
  do j = 1, 4
    do i = 1, 4
      point = point + 1
      point_uv(1,point) = real ( i - 1, kind = rk ) / 3.0D+00
      point_uv(2,point) = real ( j - 1, kind = rk ) / 3.0D+00
    end do
  end do

  call bezier_patch_evaluate ( node_num, node_xyz, rectangle_num, &
    rectangle_node, patch, point_num, point_uv, point_xyz )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  (U,V) --> (X,Y,Z) coordinates:'
  write ( *, '(a)' ) ' '

  do point = 1, point_num
    write ( *, '(2x,i2,2x,5g12.4)' ) point, point_uv(1:2,point), &
      point_xyz(1:3,point)
  end do

  return
end
