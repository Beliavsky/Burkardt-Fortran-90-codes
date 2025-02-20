program main

!*****************************************************************************80
!
!! gmsh_io_test() tests gmsh_io().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 October 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'gmsh_io_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test gmsh_io().'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'GMSH_IO_TEST():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 gets the example 2D data and writes it to a file.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 October 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, allocatable :: element_node(:,:)
  integer element_num
  integer element_order
  character ( len = 255 ) gmsh_filename
  integer m
  integer node_num
  real ( kind = rk ), allocatable :: node_x(:,:)

  gmsh_filename = 'example_2d.msh'

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  Get example 2D data, write to a file.'
!
!  Get sizes.
!
  call gmsh_mesh2d_node_size_example ( node_num, m )

  call gmsh_mesh2d_element_size_example ( element_num, element_order  )
!
!  Print the sizes.
!
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Number of nodes = ', node_num
  write ( *, '(a,i4)' ) '  Spatial dimension = ', m
  write ( *, '(a,i4)' ) '  Number of elements = ', element_num
  write ( *, '(a,i4)' ) '  Order of elements = ', element_order
!
!  Allocate memory.
!
  allocate ( node_x(1:m,1:node_num) )
  allocate ( element_node(1:element_order,1:element_num) )
!
!  Get the data.
!
  call gmsh_mesh2d_node_data_example ( node_num, m, node_x )

  call gmsh_mesh2d_element_data_example ( element_num, element_order, &
    element_node )
!
!  Print some of the data.
!
  call r8mat_transpose_print_some ( m, node_num, node_x, &
    1, 1, m, 10, '  Coordinates for first 10 nodes:' )

  call i4mat_transpose_print_some ( element_order, element_num, element_node, &
    1, 1, element_order, 10, '  Node connectivity of first 10 elements:' )
!
!  Write the GMSH file.
!
  call gmsh_mesh2d_write ( gmsh_filename, m, node_num, node_x, &
    element_order, element_num, element_node )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Wrote example data to file "' // &
    trim ( gmsh_filename ) // '"'
!
!  Clean up.
!
  deallocate ( element_node )
  deallocate ( node_x )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 reads the example data from a file.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 October 2014
!
!  Author:
!
!   John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, allocatable :: element_node(:,:)
  integer element_num
  integer element_order
  character ( len = 255 ) gmsh_filename
  integer m
  integer node_num
  real ( kind = rk ), allocatable :: node_x(:,:)

  gmsh_filename = 'example_2d.msh'

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST02:'
  write ( *, '(a)' ) '  Read data from a file.'
!
!  Get the data size.
!
  call gmsh_size_read ( gmsh_filename, node_num, m, element_num, &
    element_order )
!
!  Print the sizes.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Node data read from file "' // trim ( gmsh_filename ) // '"'
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Number of nodes = ', node_num
  write ( *, '(a,i4)' ) '  Spatial dimension = ', m
  write ( *, '(a,i4)' ) '  Number of elements = ', element_num
  write ( *, '(a,i4)' ) '  Element order = ', element_order
!
!  Allocate memory.
!
  allocate ( node_x(1:m,1:node_num) )
  allocate ( element_node(1:element_order,1:element_num) )
!
!  Get the data.
!
  call gmsh_data_read ( gmsh_filename, m, node_num, node_x, &
    element_order, element_num, element_node )
!
!  Print some of the data.
!
  call r8mat_transpose_print_some ( m, node_num, node_x, &
    1, 1, m, 10, '  Coordinates for first 10 nodes:' )
!
  call i4mat_transpose_print_some ( element_order, element_num, element_node, &
    1, 1, element_order, 10, '  Connectivity for first 10 elements:' )
!
!  Clean up.
!
  deallocate ( element_node )
  deallocate ( node_x )

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 8 ) ampm
  integer d
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  integer values(8)
  integer y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
