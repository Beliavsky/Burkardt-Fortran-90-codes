program main

!*****************************************************************************80
!
!! niederreiter2_test() tests niederreiter2().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'niederreiter2_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test NIEDERREITER2().'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NIEDERREITER2_TEST():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests NIEDERREITER2.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 August 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_max = 4

  integer dim_num
  integer i
  real ( kind = rk ) r(dim_max)
  integer seed
  integer seed_in
  integer seed_out

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  NIEDERREITER2 computes the next element of '
  write ( *, '(a)' ) '  a Niederreiter quasirandom sequence using base 2.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we call NIEDERREITER2 repeatedly.'

  do dim_num = 2, dim_max

    seed = 0

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Using dimension DIM_NUM =   ', dim_num
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '    Seed    Seed     Niederreiter2'
    write ( *, '(a)' ) '      In     Out'
    write ( *, '(a)' ) ' '
    do i = 0, 110
      seed_in = seed
      call niederreiter2 ( dim_num, seed, r )
      seed_out = seed
      if ( i <= 11 .or. 95 <= i ) then
        write ( *, '(i8,i8,4f10.4)' ) seed_in, seed_out, r(1:dim_num)
      else if ( i == 12 ) then
        write ( *, '(a)' ) '......................'
      end if
    end do

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests NIEDERREITER2.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 August 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  integer i
  real ( kind = rk ) r(dim_num)
  integer seed
  integer seed_in
  integer seed_out

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  NIEDERREITER2 computes the next element of '
  write ( *, '(a)' ) '  a Niederreiter quasirandom sequence using base 2.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we demonstrate how the SEED can be'
  write ( *, '(a)' ) '  manipulated to skip ahead in the sequence, or'
  write ( *, '(a)' ) '  to come back to any part of the sequence.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Using dimension DIM_NUM =   ', dim_num

  seed = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Seed    Seed     Niederreiter2'
  write ( *, '(a)' ) '      In     Out'
  write ( *, '(a)' ) ' '
  do i = 0, 10
    seed_in = seed
    call niederreiter2 ( dim_num, seed, r )
    seed_out = seed
    write ( *, '(2i8,4f10.4)' ) seed_in, seed_out, r(1:dim_num)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Jump ahead by increasing SEED:'
  write ( *, '(a)' ) ' '

  seed = 100

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Seed    Seed     Niederreiter2'
  write ( *, '(a)' ) '      In     Out'
  write ( *, '(a)' ) ' '
  do i = 1, 5
    seed_in = seed
    call niederreiter2 ( dim_num, seed, r )
    seed_out = seed
    write ( *, '(2i8,4f10.4)' ) seed_in, seed_out, r(1:dim_num)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Jump back by decreasing SEED:'
  write ( *, '(a)' ) ' '

  seed = 3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Seed    Seed     Niederreiter2'
  write ( *, '(a)' ) '      In     Out'
  write ( *, '(a)' ) ' '
  do i = 0, 10
    seed_in = seed
    call niederreiter2 ( dim_num, seed, r )
    seed_out = seed
    write ( *, '(2i8,4f10.4)' ) seed_in, seed_out, r(1:dim_num)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Jump ahead by increasing SEED:'
  write ( *, '(a)' ) ' '

  seed = 98

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Seed    Seed     Niederreiter2'
  write ( *, '(a)' ) '      In     Out'
  write ( *, '(a)' ) ' '
  do i = 1, 5
    seed_in = seed
    call niederreiter2 ( dim_num, seed, r )
    seed_out = seed
    write ( *, '(2i8,4f10.4)' ) seed_in, seed_out, r(1:dim_num)
  end do

  return
end
