program main

!*****************************************************************************80
!
!! TRI_SURFACE_IO_TEST() tests TRI_SURFACE_IO().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 February 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRI_SURFACE_IO_TEST():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test TRI_SURFACE_IO().'

  call test01 ( 'sphere_nodes.txt', 'sphere_elements.txt' )
  call test02 ( 'sphere_nodes.txt', 'sphere_elements.txt' )
  call test03 ( 'cube_nodes.txt', 'cube_elements.txt' )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRI_SURFACE_IO_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( node_file_name, triangle_file_name )

!*****************************************************************************80
!
!! TEST01 tests TRI_SURFACE_SIZE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 September 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer dim_num
  character ( len = *  ) node_file_name
  integer node_num
  integer order_num
  integer triangle_num
  character ( len = *  ) triangle_file_name

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  TRI_SURFACE_SIZE determines the size of various objects'
  write ( *, '(a)' ) '  in a TRI_SURFACE file.'

  call tri_surface_size ( node_file_name, triangle_file_name, dim_num, node_num, &
    order_num, triangle_num )

  call tri_surface_size_print ( node_file_name, triangle_file_name, dim_num, &
    node_num, order_num, triangle_num )

  return
end
subroutine test02 ( node_file_name, triangle_file_name )

!*****************************************************************************80
!
!! TEST02 tests TRI_SURFACE_READ.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 September 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num
  character ( len = *  ) node_file_name
  integer node_num
  real ( kind = rk ), allocatable, dimension(:,:) :: node_xyz
  integer order_num
  integer triangle_num
  character ( len = *  ) triangle_file_name
  integer, allocatable, dimension(:,:) :: triangle_node

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  TRI_SURFACE_READ reads data from a TRI_SURFACE file.'

  call tri_surface_size ( node_file_name, triangle_file_name, dim_num, node_num, &
    order_num, triangle_num )

  call tri_surface_size_print ( node_file_name, triangle_file_name, dim_num, &
    node_num, order_num, triangle_num )

  allocate ( node_xyz(dim_num,node_num) )
  allocate ( triangle_node(order_num,triangle_num) )

  call tri_surface_read ( node_file_name, triangle_file_name, dim_num, node_num, &
    order_num, triangle_num, node_xyz, triangle_node )

  call tri_surface_print ( node_file_name, triangle_file_name, dim_num, node_num, &
    order_num, triangle_num, node_xyz, triangle_node )

  deallocate ( node_xyz )
  deallocate ( triangle_node )

  return
end
subroutine test03 ( node_file_name, triangle_file_name )

!*****************************************************************************80
!
!! TEST03 tests TRI_SURFACE_WRITE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 September 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer, parameter :: node_num = 8
  integer, parameter :: order_num = 3
  integer, parameter :: triangle_num = 12

  character ( len = *  ) node_file_name
  real ( kind = rk ), dimension ( dim_num, node_num ) :: node_xyz = reshape ( (/ &
    0.0D+00, 0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00, 0.0D+00, &
    1.0D+00, 1.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, 1.0D+00, &
    1.0D+00, 0.0D+00, 1.0D+00, &
    0.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00 /), (/ dim_num, node_num /) )
  character ( len = *  ) triangle_file_name
  integer, dimension ( order_num, triangle_num ) :: triangle_node = reshape ( (/ &
    1, 3, 2, &
    2, 3, 4, &
    1, 6, 5, &
    1, 2, 6, &
    3, 7, 4, &
    4, 7, 8, &
    5, 6, 8, &
    5, 8, 7, &
    1, 5, 7, &
    1, 7, 3, &
    2, 4, 6, &
    6, 4, 8 /), (/ order_num, triangle_num /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  TRI_SURFACE_WRITE writes TRI_SURFACE data to two files.'

  call tri_surface_write ( node_file_name, triangle_file_name, dim_num, &
    node_num, order_num, triangle_num, node_xyz, triangle_node )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Graphics data was written to:'
  write ( *, '(a)' ) '    Node file:     "' // trim ( node_file_name ) // '".'
  write ( *, '(a)' ) '    Triangle file: "' // trim ( triangle_file_name ) // '".'

  return
end
