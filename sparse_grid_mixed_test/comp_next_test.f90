program main

!****************************************************************************80
!
!! comp_next_test() tests comp_next().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 December 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Fabio Nobile, Raul Tempone, Clayton Webster,
!    A Sparse Grid Stochastic Collocation Method for Partial Differential
!    Equations with Random Input Data,
!    SIAM Journal on Numerical Analysis,
!    Volume 46, Number 5, 2008, pages 2309-2345.
!
  implicit none

  integer, parameter :: test_num = 12

  integer dim_num
  integer :: dim_num_array(test_num) = (/ &
    2, 2, 2, 2, 2, &
    3, 3, 3, 3, 3, &
    4, 4 /)
  integer level_max
  integer :: level_max_array(test_num) = (/ &
    0, 1, 2, 3, 4, &
    0, 1, 2, 3, 4, &
    2, 3 /)
  integer test

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COMP_NEXT_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Call COMP_NEXT_TEST with various arguments..'

  do test = 1, test_num

    dim_num = dim_num_array(test)
    level_max = level_max_array(test)
    call comp_next_test ( dim_num, level_max )

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COMP_NEXT_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine comp_next_test ( dim_num, level_max )

!****************************************************************************80
!
!! COMP_NEXT_TEST tests COMP_NEXT, which computes 1D level vectors.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer LEVEL_MAX, the maximum level.
!
  implicit none

  integer dim_num

  integer h
  integer i
  integer level
  integer level_1d(dim_num)
  integer level_max
  integer level_min
  logical more_grids
  integer t

  level_min = max ( 0, level_max + 1 - dim_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COMP_NEXT_TEST'
  write ( *, '(a)' ) '  COMP_NEXT generates, one at a time, vectors'
  write ( *, '(a)' ) '  LEVEL_1D(1:DIM_NUM) whose components add up to LEVEL.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a,i8)' ) '  We call with:'
  write ( *, '(a,i8)' ) '  DIM_NUM = ', dim_num
  write ( * , '(2x,i8,a,i8)' ) level_min, &
    ' = LEVEL_MIN <= LEVEL <= LEVEL_MAX = ', level_max
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     LEVEL     INDEX  LEVEL_1D Vector'
!
!  The outer loop generates values of LEVEL from LEVEL_MIN to LEVEL_MAX.
!
  do level = level_min, level_max

    write ( *, '(a)' ) ' '
!
!  The inner loop generates vectors LEVEL_1D(1:DIM_NUM) whose components
!  add up to LEVEL.
!
    more_grids = .false.
    h = 0
    t = 0
    i = 0

    do

      call comp_next ( level, dim_num, level_1d, more_grids, h, t )

      i = i + 1
      write ( *, '(2x,i8,2x,i8,6(2x,i8))'  ) level, i, level_1d(1:dim_num)

      if ( .not. more_grids ) then
        exit
      end if

    end do

  end do

  return
end
