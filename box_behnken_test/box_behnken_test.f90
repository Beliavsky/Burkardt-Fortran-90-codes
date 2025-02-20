program main

!*****************************************************************************80
!
!! box_behnken_test() tests box_behnken().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'box_behnken_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test box_behnken().'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'box_behnken_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! test01() tests box_behnken().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ), dimension (dim_num,2) :: range = reshape ( (/ &
    0.0D+00, 10.0D+00,  5.0D+00, &
    1.0D+00, 11.0D+00, 15.0D+00 /), (/ dim_num, 2 /) )
  integer x_num
  real ( kind = rk ), allocatable, dimension(:,:) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test01():'
  write ( *, '(a)' ) '  box_behnken() computes a Box-Behnken dataset.'

  call r8mat_transpose_print ( dim_num, 2, range, &
    '  The ranges:' )

  call box_behnken_size ( dim_num, x_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a,i6)' ) '  For dimension DIM_NUM = ', dim_num, &
    ' the Box-Behnken design is of size ', x_num

  allocate ( x(1:dim_num,1:x_num) )

  call box_behnken ( dim_num, x_num, range, x )

  call r8mat_transpose_print ( dim_num, x_num, x, &
    '  The Box-Behnken design:' )

  deallocate ( x )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! test02() writes a box_behnken() dataset to a file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 4

  character ( len = 80 ) file_out_name
  real ( kind = rk ), dimension (dim_num,2) :: range = reshape ( (/ &
    0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00 /), (/ dim_num, 2 /) )
  integer x_num
  real ( kind = rk ), allocatable, dimension(:,:) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test02()'
  write ( *, '(a)' ) '  r8mat_write() writes a Box-Behnken dataset'
  write ( *, '(a)' ) '  to a file.'

  call r8mat_transpose_print ( dim_num, 2, range, &
    '  The ranges:' )

  call box_behnken_size ( dim_num, x_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a,i6)' ) '  For dimension DIM_NUM = ', dim_num, &
    ' the Box-Behnken design is of size ', x_num

  allocate ( x(1:dim_num,1:x_num) )

  call box_behnken ( dim_num, x_num, range, x )

  file_out_name = 'box_behnken_04_33.txt'

  call r8mat_write ( file_out_name, dim_num, x_num, x )

  deallocate ( x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The data was written to the file "' // &
    trim ( file_out_name ) // '".'

  return
end

