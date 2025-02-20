program main

!*****************************************************************************80
!
!! GRF_IO_TEST tests GRF_IO.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GRF_IO_TEST():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test GRF_IO().'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GRF_IO_TEST():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests GRF_HEADER_WRITE and GRF_DATA_WRITE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer edge_num
  integer, allocatable :: edge_data(:)
  integer, allocatable :: edge_pointer(:)
  integer node_num
  character ( len = 80 ) :: output_filename = 'coxeter.grf'
  real ( kind = rk ), allocatable :: xy(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01():'
  write ( *, '(a)' ) '  GRF_HEADER_WRITE writes the header of a GRF file.'
  write ( *, '(a)' ) '  GRF_DATAWRITE writes the data of a GRF file.'

  call grf_example_size ( node_num, edge_num )

  call grf_header_print ( node_num, edge_num )

  allocate ( edge_pointer(1:node_num+1) )
  allocate ( edge_data(1:edge_num) )
  allocate ( xy(1:2,1:node_num) )

  call grf_example ( node_num, edge_num, edge_pointer, edge_data, xy )

  call grf_write ( output_filename, node_num, edge_num, edge_pointer, &
    edge_data, xy )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Wrote the GRF file "' // trim ( output_filename ) // '".'

  deallocate ( edge_data )
  deallocate ( edge_pointer )
  deallocate ( xy )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests GRF_HEADER_READ and GRF_DATA_READ.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer edge_num
  character ( len = 80 ) :: input_filename = 'coxeter.grf'
  integer, allocatable :: edge_data(:)
  integer, allocatable :: edge_pointer(:)
  integer node_num
  real ( kind = rk ), allocatable :: xy(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  GRF_HEADER_READ reads the header of a GRF file.'
  write ( *, '(a)' ) '  GRF_DATA_READ reads the data of a GRF file.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Reading the GRF file "' // trim ( input_filename ) // '".'

  call grf_header_read ( input_filename, node_num, edge_num )

  call grf_header_print ( node_num, edge_num )

  allocate ( edge_pointer(1:node_num+1) )
  allocate ( edge_data(1:edge_num) )
  allocate ( xy(1:2,1:node_num) )

  call grf_data_read ( input_filename, node_num, edge_num, edge_pointer, &
    edge_data, xy )

  call grf_data_print ( node_num, edge_num, edge_pointer, edge_data, xy )

  deallocate ( edge_data )
  deallocate ( edge_pointer )
  deallocate ( xy )

  return
end
