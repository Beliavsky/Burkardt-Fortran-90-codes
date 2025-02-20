program main

!*****************************************************************************80
!
!! xy_io_test() tests xy_io().
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
  write ( *, '(a)' ) 'xy_io_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test xy_io().'

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
  write ( *, '(a)' ) 'xy_io_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests XY_EXAMPLE, XY_WRITE.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 June 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: point_num = 300

  character ( len = 80 ) :: file_name = 'xy_io_test01.xy'
  real ( kind = rk8 ) xy(2,point_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  XY_EXAMPLE sets up sample XY data.'
  write ( *, '(a)' ) '  XY_WRITE writes an XY file.'

  call xy_example ( point_num, xy )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  XY_EXAMPLE has created the data.'

  call xy_write ( file_name, point_num, xy )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  XY_WRITE wrote the header and data for "' &
    // trim ( file_name ) //'".'
  write ( *, '(a,i8)' ) '  Number of points = ', point_num

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests XY_READ.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 June 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  character ( len = 80 ) :: file_name = 'xy_io_test02.xy'
  integer i
  integer k
  integer point_num
  real ( kind = rk8 ), allocatable :: xy(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  XY_READ reads an XY file.'

  call xy_write_test ( file_name )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  XY_WRITE_TEST created some data.'

  call xy_header_read ( file_name, point_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  XY_HEADER_READ has read the header.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of points = ', point_num

  allocate ( xy(2,point_num) )

  call xy_data_read ( file_name, point_num, xy )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  XY_DATA_READ has read the data.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sample data:'
  write ( *, '(a)' ) ' '

  do k = 1, 10
    i = ( ( 10 - k ) * 1 + ( k - 1 ) * point_num ) / ( 10 - 1 )
    write ( *, '(i4,2x,f10.4,2x,f10.4)' ) i, xy(1:2,i)
  end do

  deallocate ( xy )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests XYL_EXAMPLE, XYL_WRITE.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    31 December 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, allocatable :: line_data(:)
  integer, allocatable :: line_pointer(:)
  integer line_data_num
  integer line_num
  integer point_num
  real ( kind = rk8 ), allocatable :: xy(:,:)
  character ( len = 80 ) :: xy_filename = 'house.xy'
  character ( len = 80 ) :: xyl_filename = 'house.xyl'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  XYL_EXAMPLE sets up XY and XYL data.'

  call xyl_example_size ( point_num, line_num, line_data_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Example has:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of points     = ', point_num
  write ( *, '(a,i8)' ) '  Number of lines      = ', line_num
  write ( *, '(a,i8)' ) '  Number of line items = ', line_data_num

  allocate ( line_data(line_data_num) )
  allocate ( line_pointer(line_num+1) )
  allocate ( xy(2,point_num) )

  call xyl_example ( point_num, line_num, line_data_num, xy, &
    line_pointer, line_data )

  call xy_write ( xy_filename, point_num, xy )

  call xyl_write ( xyl_filename, point_num, line_num, line_data_num, &
    line_pointer, line_data )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Wrote the XY file "' // trim ( xy_filename ) // '",'
  write ( *, '(a)' ) '  and the XYL file "' // trim ( xyl_filename ) // '".'

  deallocate ( line_data )
  deallocate ( line_pointer )
  deallocate ( xy )

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! test04() tests xyl_header_read(), xyl_data_read().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    31 December 2008
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
  real ( kind = rk8 ), allocatable :: xy(:,:)
  character ( len = 80 ) :: xy_filename = 'house.xy'
  character ( len = 80 ) :: xyl_filename = 'house.xyl'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  XY_HEADER_READ  reads the header of an XY  file.'
  write ( *, '(a)' ) '  XY_DATA_READ    reads the data   of an XY  file.'
  write ( *, '(a)' ) '  XYL_HEADER_READ reads the header of an XYL file.'
  write ( *, '(a)' ) '  XYL_DATA_READ   reads the data   of an XYL file.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Examine XY file "' // trim ( xy_filename ) // '".'

  call xy_header_read ( xy_filename, point_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of points     = ', point_num

  allocate ( xy(2,point_num) )

  call xy_data_read ( xy_filename, point_num, xy )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Point data:'
  write ( *, '(a)' ) ' '

  do i = 1, point_num
    write ( *, '(2x,i4,2x,f10.4,2x,f10.4)' ) i, xy(1:2,i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Examine XYL file "' // trim ( xyl_filename ) // '".'

  call xyl_header_read ( xyl_filename, line_num, line_data_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of lines      = ', line_num
  write ( *, '(a,i8)' ) '  Number of line items = ', line_data_num

  allocate ( line_data(line_data_num) )
  allocate ( line_pointer(line_num+1) )

  call xyl_data_read ( xyl_filename, line_num, line_data_num, line_pointer, line_data )

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
  deallocate ( xy )

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests XYF_EXAMPLE, XYF_WRITE.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, allocatable :: face_data(:)
  integer, allocatable :: face_pointer(:)
  integer face_data_num
  integer face_num
  integer point_num
  real ( kind = rk8 ), allocatable :: xy(:,:)
  character ( len = 80 ) :: xy_filename = 'annulus.xy'
  character ( len = 80 ) :: xyf_filename = 'annulus.xyf'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  XYF_EXAMPLE sets up XY and XYF data.'

  call xyf_example_size ( point_num, face_num, face_data_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Example has:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of points     = ', point_num
  write ( *, '(a,i8)' ) '  Number of faces      = ', face_num
  write ( *, '(a,i8)' ) '  Number of face items = ', face_data_num

  allocate ( face_data(face_data_num) )
  allocate ( face_pointer(face_num+1) )
  allocate ( xy(2,point_num) )

  call xyf_example ( point_num, face_num, face_data_num, xy, &
    face_pointer, face_data )

  call xy_write ( xy_filename, point_num, xy )

  call xyf_write ( xyf_filename, point_num, face_num, face_data_num, &
    face_pointer, face_data )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Wrote the XY file "' // trim ( xy_filename ) // '",'
  write ( *, '(a)' ) '  and the XYF file "' // trim ( xyf_filename ) // '".'

  deallocate ( face_data )
  deallocate ( face_pointer )
  deallocate ( xy )

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests XYF_HEADER_READ, XYF_DATA_READ.
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

  integer face
  integer, allocatable :: face_data(:)
  integer, allocatable :: face_pointer(:)
  integer face_data_num
  integer face_num
  integer i
  integer j
  integer point_num
  real ( kind = rk8 ), allocatable :: xy(:,:)
  character ( len = 80 ) :: xy_filename = 'annulus.xy'
  character ( len = 80 ) :: xyf_filename = 'annulus.xyf'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  XY_HEADER_READ  reads the header of an XY  file.'
  write ( *, '(a)' ) '  XY_DATA_READ    reads the data   of an XY  file.'
  write ( *, '(a)' ) '  XYF_HEADER_READ reads the header of an XYF file.'
  write ( *, '(a)' ) '  XYF_DATA_READ   reads the data   of an XYF file.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Examine XY file "' // trim ( xy_filename ) // '".'

  call xy_header_read ( xy_filename, point_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of points     = ', point_num

  allocate ( xy(2,point_num) )

  call xy_data_read ( xy_filename, point_num, xy )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Point data:'
  write ( *, '(a)' ) ' '

  do i = 1, point_num
    write ( *, '(2x,i4,2x,f10.4,2x,f10.4)' ) i, xy(1:2,i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Examine XYF file "' // trim ( xyf_filename ) // '".'

  call xyf_header_read ( xyf_filename, face_num, face_data_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of faces      = ', face_num
  write ( *, '(a,i8)' ) '  Number of face items = ', face_data_num

  allocate ( face_data(face_data_num) )
  allocate ( face_pointer(face_num+1) )

  call xyf_data_read ( xyf_filename, face_num, face_data_num, face_pointer, face_data )

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
  deallocate ( xy )

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

