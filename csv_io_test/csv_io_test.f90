program main

!*****************************************************************************80
!
!! csv_io_test() tests csv_io().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 November 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CSV_IO_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test CSV_IO.'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CSV_IO_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 writes a variety of data items to a CSV file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 November 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 10

  integer, dimension ( n ) :: column2 = (/ &
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10 /)
  character ( len = 10 ), dimension ( n ) :: column3 = (/ &
    'one       ', 'two       ', 'three     ', 'four      ', 'five      ', &
    'six       ', 'seven     ', 'eight     ', 'nine      ', 'ten       ' /)
  real ( kind = rk ), dimension ( n ) :: column4 = (/ &
    1.1D+00, 2.2D+00, 3.3D+00, 4.4D+00,  5.5D+00, &
    6.6D+00, 7.7D+00, 8.8D+00, 9.9D+00, 10.10D+00 /)
  character ( len = 80 ) :: csv_file_name = 'csv_4col_5row.csv'
  integer csv_file_unit
  integer i
  character ( len = 80 ) record

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Test routines to create a CSV file.'

  call csv_file_open_write ( csv_file_name, csv_file_unit )
!
!  Write the header line.
!
  record = ' '
  call csv_record_append_s ( 'I4', record )
  call csv_record_append_s ( 'S', record )
  call csv_record_append_s ( 'R8', record )

  call csv_file_header_write ( csv_file_unit, record )

  do i = 1, n

    record = ' '

    call csv_record_append_i4 ( column2(i), record )
    call csv_record_append_s  ( column3(i), record )
    call csv_record_append_r8 ( column4(i), record )

    call csv_file_record_write ( csv_file_unit, record )

  end do

  call csv_file_close_write ( csv_file_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to "' // trim ( csv_file_name ) // '".'

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 reads a CSV file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 November 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 255 ) :: csv_file_name = 'csv_4col_5row.csv'
  integer csv_file_status
  integer csv_file_unit
  integer csv_record_status
  integer i
  integer line_num
  character ( len = 255 ) record
  integer value_count

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Read data from CSV file created by TEST01.'

  call csv_file_line_count ( csv_file_name, line_num )

  write ( *, '(a,i8,a)' ) '  File contains ', line_num, ' lines.'

  call csv_file_open_read ( csv_file_name, csv_file_unit )

  do i = 1, line_num
    read ( csv_file_unit, '(a)', iostat = csv_file_status ) record
    write ( *, '(i4, a)' ) i, trim ( record )
    call csv_value_count ( record, csv_record_status, value_count )
    write ( *, '(i4,2x,i4,a)' ) i, value_count, ' values read.'
  end do

  call csv_file_close_read ( csv_file_unit )

  return
end

