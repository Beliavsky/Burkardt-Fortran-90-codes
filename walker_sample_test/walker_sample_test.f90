program main

!*****************************************************************************80
!
!! walker_sample_test() tests walker_sample().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'walker_sample_test():'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) '  Test walker_sample().'
  
  call i4_choose_test ( )
  call i4_uniform_ab_test ( )
  call normalize_test ( )
  call r8vec_print_test ( )
  call random_permutation_test ( )
  call walker_build_test ( )
  call walker_sampler_test ( )
  call walker_verify_test ( )
  call zipf_probability_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'walker_sample_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
