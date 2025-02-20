program main

!*****************************************************************************80
!
!! CUBE_FELIPPA_RULE_TEST() tests CUBE_FELIPPA_RULE().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 July 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer degree_max

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CUBE_FELIPPA_RULE_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test CUBE_FELIPPA_RULE().'

  degree_max = 4
  call cube_monomial_test ( degree_max )

  degree_max = 6
  call cube_quad_test ( degree_max )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CUBE_FELIPPA_RULE_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end

