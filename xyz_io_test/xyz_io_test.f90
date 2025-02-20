program main

!*****************************************************************************80
!
!! xyz_io_test() tests xyz_io().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'xyz_io_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test xyz_io().'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'xyz_io_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests XYZ_EXAMPLE, XYZ_WRITE.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  character ( len = 80 ) :: file_name = 'helix.xyz'
  integer point_num
  real ( kind = rk8 ), allocatable :: xyz(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  XYZ_EXAMPLE sets up sample XYZ data.'
  write ( *, '(a)' ) '  XYZ_WRITE writes an XYZ file.'

  call xyz_example_size ( point_num )

  write ( *, '(a,i8)' ) '  Example dataset size is ', point_num

  allocate ( xyz(3,1:point_num) )

  call xyz_example ( point_num, xyz )

  call xyz_write ( file_name, point_num, xyz )

  deallocate ( xyz )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  XYZ_WRITE wrote the header and data for "' &
    // trim ( file_name ) //'".'
  write ( *, '(a,i8)' ) '  Number of points = ', point_num

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests XYZ_HEADER_READ, XYZ_DATA_READ.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  character ( len = 80 ) :: file_name = 'xyz_io_test02.xyz'
  integer i
  integer k
  integer point_num
  real ( kind = rk8 ), allocatable :: xyz(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  XYZ_HEADER_READ reads the header of an XYZ file.'
  write ( *, '(a)' ) '  XYZ_DATA_READ reads the data of an XYZ file.'

  call xyz_write_test ( file_name )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  XYZ_WRITE_TEST created some data.'

  call xyz_header_read ( file_name, point_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  XYZ_HEADER_READ has read the header.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of points = ', point_num

  allocate ( xyz(3,1:point_num) )

  call xyz_data_read ( file_name, point_num, xyz )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  XYZ_DATA_READ has read the data.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sample data:'
  write ( *, '(a)' ) ' '

  do k = 1, 11
    i = ( ( 11 - k     ) * 1 &
        + (      k - 1 ) * point_num ) &
        / ( 11     - 1 )
    write ( *, '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4)' ) i, xyz(1:3,i)
  end do

  deallocate ( xyz )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests XYZL_EXAMPLE.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, allocatable, dimension ( : ) :: line_data
  integer, allocatable, dimension ( : ) :: line_pointer
  integer line_data_num
  integer line_num
  integer point_num
  real ( kind = rk8 ), allocatable :: xyz(:,:)
  character ( len = 80 ) :: xyz_filename = 'cube.xyz'
  character ( len = 80 ) :: xyzl_filename = 'cube.xyzl'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  XYZL_EXAMPLE sets up XYZ and XYZL data.'

  call xyzl_example_size ( point_num, line_num, line_data_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Example has:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of points     = ', point_num
  write ( *, '(a,i8)' ) '  Number of lines      = ', line_num
  write ( *, '(a,i8)' ) '  Number of line items = ', line_data_num

  allocate ( line_data(line_data_num) )
  allocate ( line_pointer(line_num+1) )
  allocate ( xyz(3,1:point_num) )

  call xyzl_example ( point_num, line_num, line_data_num, xyz, line_pointer, &
    line_data )

  call xyz_write ( xyz_filename, point_num, xyz )

  call xyzl_write ( xyzl_filename, point_num, line_num, line_data_num, &
    line_pointer, line_data )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Wrote the XYZ file "' // trim ( xyz_filename ) // '",'
  write ( *, '(a)' ) '  and the XYZL file "' // trim ( xyzl_filename ) // '".'

  deallocate ( line_data )
  deallocate ( line_pointer )
  deallocate ( xyz )

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests XYZL_READ.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer i
  integer j
  integer line
  integer, allocatable :: line_data(:)
  integer, allocatable :: line_pointer(:)
  integer line_data_num
  integer line_num
  integer point_num
  real ( kind = rk8 ), allocatable :: xyz(:,:)
  character ( len = 80 ) :: xyz_filename = 'cube.xyz'
  character ( len = 80 ) :: xyzl_filename = 'cube.xyzl'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  XYZ_HEADER_READ  reads the header of an XYZ  file.'
  write ( *, '(a)' ) '  XYZ_DATA_READ    reads the data   of an XYZ  file.'
  write ( *, '(a)' ) '  XYZL_HEADER_READ reads the header of an XYZL file.'
  write ( *, '(a)' ) '  XYZL_DATA_READ   reads the data   of an XYZL file.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Examine XYZ file "' // trim ( xyz_filename ) // '".'

  call xyz_header_read ( xyz_filename, point_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of points = ', point_num

  allocate ( xyz(3,1:point_num) )

  call xyz_data_read ( xyz_filename, point_num, xyz )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Point data:'
  write ( *, '(a)' ) ' '

  do i = 1, point_num
    write ( *, '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4)' ) i, xyz(1:3,i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Examine XYZL file "' // trim ( xyzl_filename ) // '".'

  call xyzl_header_read ( xyzl_filename, line_num, line_data_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of lines      = ', line_num
  write ( *, '(a,i8)' ) '  Number of line items = ', line_data_num

  allocate ( line_data(line_data_num) )
  allocate ( line_pointer(line_num+1) )

  call xyzl_data_read ( xyzl_filename, line_num, line_data_num, line_pointer, line_data )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Line pointers:'
  write ( *, '(a)' ) ' '

  do line = 1, line_num
    write ( *, '(2x,i4,2x,i8,2x,i8)' ) line, line_pointer(line), line_pointer(line+1) - 1
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Line data:'
  write ( *, '(a)' ) ' '

  do line = 1, line_num
    write ( *, '(2x,i4,4x)', advance = 'no' ) line
    do j = line_pointer(line), line_pointer(line+1) - 1
      write ( *, '(2x,i8)', advance = 'no' ) line_data(j)
    end do
    write ( *, '(a)', advance = 'yes' ) ' '
  end do

  deallocate ( line_data )
  deallocate ( line_pointer )
  deallocate ( xyz )

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests XYZF_EXAMPLE.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, allocatable, dimension ( : ) :: face_data
  integer, allocatable, dimension ( : ) :: face_pointer
  integer face_data_num
  integer face_num
  integer point_num
  real ( kind = rk8 ), allocatable :: xyz(:,:)
  character ( len = 80 ) :: xyz_filename = 'cube.xyz'
  character ( len = 80 ) :: xyzf_filename = 'cube.xyzf'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  XYZF_EXAMPLE sets up XYZ and XYZF data.'

  call xyzf_example_size ( point_num, face_num, face_data_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Example has:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of points     = ', point_num
  write ( *, '(a,i8)' ) '  Number of faces      = ', face_num
  write ( *, '(a,i8)' ) '  Number of face items = ', face_data_num

  allocate ( face_data(face_data_num) )
  allocate ( face_pointer(face_num+1) )
  allocate ( xyz(3,1:point_num) )

  call xyzf_example ( point_num, face_num, face_data_num, xyz, face_pointer, &
    face_data )

  call xyz_write ( xyz_filename, point_num, xyz )

  call xyzf_write ( xyzf_filename, point_num, face_num, face_data_num, &
    face_pointer, face_data )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Wrote the XYZ file "' // trim ( xyz_filename ) // '",'
  write ( *, '(a)' ) '  and the XYZF file "' // trim ( xyzf_filename ) // '".'

  deallocate ( face_data )
  deallocate ( face_pointer )
  deallocate ( xyz )

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests XYZF_READ.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer i
  integer j
  integer face
  integer, allocatable :: face_data(:)
  integer, allocatable :: face_pointer(:)
  integer face_data_num
  integer face_num
  integer point_num
  real ( kind = rk8 ), allocatable :: xyz(:,:)
  character ( len = 80 ) :: xyz_filename = 'cube.xyz'
  character ( len = 80 ) :: xyzf_filename = 'cube.xyzf'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  XYZ_HEADER_READ  reads the header of an XYZ  file.'
  write ( *, '(a)' ) '  XYZ_DATA_READ    reads the data   of an XYZ  file.'
  write ( *, '(a)' ) '  XYZF_HEADER_READ reads the header of an XYZF file.'
  write ( *, '(a)' ) '  XYZF_DATA_READ   reads the data   of an XYZF file.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Examine XYZ file "' // trim ( xyz_filename ) // '".'

  call xyz_header_read ( xyz_filename, point_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of points = ', point_num

  allocate ( xyz(3,1:point_num) )

  call xyz_data_read ( xyz_filename, point_num, xyz )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Point data:'
  write ( *, '(a)' ) ' '

  do i = 1, point_num
    write ( *, '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4)' ) i, xyz(1:3,i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Examine XYZF file "' // trim ( xyzf_filename ) // '".'

  call xyzf_header_read ( xyzf_filename, face_num, face_data_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of faces      = ', face_num
  write ( *, '(a,i8)' ) '  Number of face items = ', face_data_num

  allocate ( face_data(face_data_num) )
  allocate ( face_pointer(face_num+1) )

  call xyzf_data_read ( xyzf_filename, face_num, face_data_num, face_pointer, face_data )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Face pointers:'
  write ( *, '(a)' ) ' '

  do face = 1, face_num
    write ( *, '(2x,i4,2x,i8,2x,i8)' ) face, face_pointer(face), face_pointer(face+1) - 1
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Face data:'
  write ( *, '(a)' ) ' '

  do face = 1, face_num
    write ( *, '(2x,i4,4x)', advance = 'no' ) face
    do j = face_pointer(face), face_pointer(face+1) - 1
      write ( *, '(2x,i8)', advance = 'no' ) face_data(j)
    end do
    write ( *, '(a)', advance = 'yes' ) ' '
  end do

  deallocate ( face_data )
  deallocate ( face_pointer )
  deallocate ( xyz )

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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end

