program main

!*****************************************************************************80
!
!! line_felippa_rule_test() tests line_felippa_rule().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 September 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer degree_max

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'line_felippa_rule_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the LINE_FELIPPA_RULE library.'

  degree_max = 4
  call line_monomial_test ( degree_max )

  degree_max = 11
  call line_quad_test ( degree_max )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'line_felippa_rule_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end

