program main

!*****************************************************************************80
!
!! square_felippa_rule_test() tests square_felippa_rule().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 September 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer degree_max

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'square_felippa_rule_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test square_felippa_rule().'

  degree_max = 4
  call square_monomial_test ( degree_max )

  degree_max = 5
  call square_quad_test ( degree_max )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'square_felippa_rule_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end

