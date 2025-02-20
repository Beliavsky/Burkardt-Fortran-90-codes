program main

!*****************************************************************************80
!
!! weekday_test() tests weekday().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 July 2017
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'weekday_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test weekday().'

  call jed_to_weekday_test ( )
  call ymd_to_weekday_common_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'weekday_test():'
  write ( *, '(a)' ) '  Noraml end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
