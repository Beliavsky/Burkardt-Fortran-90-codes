program main

!*****************************************************************************80
!
!! vtk_io_test() tests vtk_io().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'vtk_io_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test vtk_io().'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'vtk_io_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests VTK_PUVW_WRITE.
!
!  Discussion:
!
!    VTK_PUVW_WRITE writes pressure (P) and velocity (UVW) for a 3D
!    fluid flow calculation.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: element_num = 6
  integer, parameter :: element_order = 8
  integer, parameter :: node_num = 24

  integer :: element_node(element_order,element_num) = reshape ( &
    (/ &
      1,  2,  5,  6, 13, 14, 17, 18, &
      2,  3,  6,  7, 14, 15, 18, 19, &
      3,  4,  7,  8, 15, 16, 19, 20, &
      5,  6,  9, 10, 17, 18, 21, 22, &
      6,  7, 10, 11, 18, 19, 22, 23, &
      7,  8, 11, 12, 19, 20, 23, 24  &
    /), (/ element_order, element_num /) )
  integer i
  integer j
  integer k
  integer node
  character ( len = 80 ) output_filename
  integer output_unit
  real ( kind = rk ) p(node_num)
  character ( len = 80 ) title
  real ( kind = rk ) uvw(3,node_num)
  real ( kind = rk ) x
  real ( kind = rk ) xyz(3,node_num)
  real ( kind = rk ) y
  real ( kind = rk ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  VTK_PUVW_WRITE writes 3d fluid data, pressure and '
  write ( *, '(a)' ) '  velocity, to a VTK file.'

  output_filename = 'puvw_data.txt'
  title = 'Sample data for VTK_PUVW_WRITE.'

  node = 0
  do k = 1, 2
    z = real ( k - 1, kind = rk )
    do j = 1, 3
      y = real ( j - 1, kind = rk )
      do i = 1, 4
        x = real ( i - 1, kind = rk )
        node = node + 1
        xyz(1:3,node) = (/ x, y, z /)
        p(node_num) = 10.0D+00 * x
        uvw(1,node_num) = 2.0D+00 * x
        uvw(2,node_num) = 3.0D+00 * y
        uvw(3,node_num) = 4.0D+00 * z
      end do
    end do
  end do

  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, status = 'replace' )

  call vtk_puvw_write ( output_unit, title, node_num, element_num, &
    element_order, xyz, element_node, p, uvw )

  close (  unit = output_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  VTK_PUVW_WRITE created the file "' &
    // trim ( output_filename ) // '"'

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! timestamp() prints the current YMDHMS date as a time stamp.
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
!  Parameters:
!
!    None
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

