program main

!*****************************************************************************80
!
!! CONTINUED_FRACTION_TEST tests CONTINUED_FRACTION.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 August 2017
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'CONTINUED_FRACTION_TEST:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  CONTINUED_FRACTION is a library for dealing with'
  write ( *, '(a)' ) '  expresssions representing a continued fraction.'

  call i4cf_evaluate_test ( )
  call i4scf_evaluate_test ( )
  call r8_to_i4scf_test ( )
  call r8cf_evaluate_test ( )
  call r8scf_evaluate_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'CONTINUED_FRACTION_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  call timestamp ( )

  stop
end

