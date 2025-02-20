program main

!*****************************************************************************80
!
!! pgma_io_test() tests pgma_io().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'pgma_io_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test PGMA_IO().'

  call test01 ( )
  call test02 ( )
  call test03 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PGMA_IO_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests PGMA_EXAMPLE, PGMA_WRITE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ncol = 300
  integer, parameter :: nrow = 300

  character ( len = 80 ) :: file_name = 'pgma_io_test_01.ascii.pgm'
  integer g(nrow,ncol)
  integer ierror

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  PGMA_EXAMPLE sets up ASCII PGM data.'
  write ( *, '(a)' ) '  PGMA_WRITE writes an ASCII PGM file.'

  call pgma_example ( nrow, ncol, g )

  call pgma_write ( file_name, nrow, ncol, g, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) 'PGMA_WRITE returns IERROR = ', ierror
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Wrote the header and data for "' &
    // trim ( file_name ) //'".'
  write ( *, '(a,i8)' ) '  Number of rows of data =    ', nrow
  write ( *, '(a,i8)' ) '  Number of columns of data = ', ncol

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests PGMA_READ_DATA, PGMA_READ_HEADER.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 80 ) :: file_name = 'pgma_io_test_02.ascii.pgm'
  integer file_unit
  integer, allocatable, dimension ( :, : ) :: g
  integer i
  integer ierror
  integer ios
  integer j
  integer k
  integer maxg
  integer ncol
  integer nrow

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  PGMA_READ reads an ASCII PGM file.'

  call pgma_write_test ( file_name )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  PGMA_WRITE_TEST created some data.'

  call get_unit ( file_unit )

  open ( unit = file_unit, file = file_name, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST02 - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file.'
    return
  end if

  call pgma_read_header ( file_unit, nrow, ncol, maxg )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  PGMA_READ_HEADER read the header.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of rows of data =    ', nrow
  write ( *, '(a,i8)' ) '  Number of columns of data = ', ncol
  write ( *, '(a,i8)' ) '  Maximum G value =           ', maxg

  allocate ( g(nrow,ncol) )

  call pgma_read_data ( file_unit, nrow, ncol, g )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  PGMA_READ_DATA read the data.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sample data:'
  write ( *, '(a)' ) ' '

  do k = 1, 10
    i = ( ( 10 - k ) * 1 + ( k - 1 ) * nrow ) / ( 10 - 1 )
    j = ( ( 10 - k ) * 1 + ( k - 1 ) * ncol ) / ( 10 - 1 )
    write ( *, '(i4,2x,i4,2x,i6)' ) i, j, g(i,j)
  end do

  call pgma_check_data ( nrow, ncol, maxg, g, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST02'
    write ( *, '(a,i8)' ) '  The data was not accepted by PGMA_CHECK_DATA.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The data was accepted by PGMA_CHECK_DATA.'

  deallocate ( g )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests PGMA_WRITE.
!
!  Discussion:
!
!    This example makes a sort of grayscale checkerboard.
!
!    The gray scale values were computed by the routine
!    GRAYSCALE_RGB in the COLORS library.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 June 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: ncol = 300
  integer, parameter :: ngray = 11
  integer, parameter :: nrow = 300

  character ( len = 80 ) :: file_name = 'pgma_io_test_03.ascii.pgm'
  integer g(nrow,ncol)
  real ( kind = rk ), dimension ( ngray ) :: gray = (/ &
    0.000D+00, 0.291D+00, 0.434D+00, 0.540D+00, 0.629D+00, &
    0.706D+00, 0.774D+00, 0.837D+00, 0.895D+00, 0.949D+00, &
    1.000D+00 /)
  integer i
  integer ierror
  integer j
  integer k

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  PGMA_WRITE writes an ASCII PGM file.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this example, we make a sort of grayscale'
  write ( *, '(a)' ) '  checkerboard.'

  do i = 1, nrow
    do j = 1, ncol
      k = ( i - 1 + j - 1 ) * ngray / min ( nrow, ncol )
      k = 1 + mod ( k, ngray )
      g(i,j) = int ( 255.0D+00 * gray(k) )
    end do
  end do

  call pgma_write ( file_name, nrow, ncol, g, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) 'PGMA_WRITE returns IERROR = ', ierror
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Wrote the header and data for "' &
      // trim ( file_name ) //'".'
    write ( *, '(a,i8)' ) '  Number of rows of data =    ', nrow
    write ( *, '(a,i8)' ) '  Number of columns of data = ', ncol
  end if

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

