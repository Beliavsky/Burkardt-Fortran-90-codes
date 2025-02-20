program main

!*****************************************************************************80
!
!! latin_center_test() tests latin_center().
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

  integer seed
  integer seed_save

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'latin_center_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test latin_center().'

  call test00 ( seed )

  seed_save = seed
  call test01 ( seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Repeat TEST01 with a different seed from the first run.'

  call test01 ( seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Repeat TEST01 with the same seed as the first run.'

  seed = seed_save
  call test01 ( seed )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'latin_center_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test00 ( seed )

!*****************************************************************************80
!
!! TEST00 tests GET_SEED and RANDOM_INITIALIZE.
!
!  Discussion:
!
!    If the LATIN routines are set to use the system random number
!    generator, then setting SEED here is enough.
!
!    If they are set to use the portable UNIFORM routine, then
!    SEED must be passed to those routines.
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

  integer seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST00'
  write ( *, '(a)' ) '  GET_SEED returns a seed for the random number'
  write ( *, '(a)' ) '  generator, based on the current time.'
  write ( *, '(a)' ) '  RANDOM_INITIALIZE uses that seed to initialize'
  write ( *, '(a)' ) '  the FORTRAN90 random number generator.'

  call get_seed ( seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  GET_SEED returns SEED = ', seed

  call random_initialize ( seed )

  return
end
subroutine test01 ( seed )

!*****************************************************************************80
!
!! TEST01 tests LATIN_CENTER.
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

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2
  integer, parameter :: point_num = 10

  integer j
  integer seed
  real ( kind = rk ) x(dim_num,point_num)

  write ( *, '(a)'     ) ' '
  write ( *, '(a)'     ) 'TEST01'
  write ( *, '(a)'     ) '  LATIN_CENTER chooses a Latin cell arrangement,'
  write ( *, '(a)'     ) '  and returns the centers of those cells.'
  write ( *, '(a)'     ) ' '
  write ( *, '(a,i8)'  ) '  Spatial dimension =  ', dim_num
  write ( *, '(a,i8)'  ) '  Number of points =   ', point_num
  write ( *, '(a,i12)' ) '  Random number SEED = ', seed

  call latin_center ( dim_num, point_num, seed, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The Latin Center Square points:'
  write ( *, '(a)' ) ' '

  do j = 1, point_num
    write ( *, '(2x,f10.4,2x,f10.4)' ) x(1:dim_num,j)
  end do

  return
end
