program main

!*****************************************************************************80
!
!! test_partial_digest_test() tests TEST_PARTIAL_DIGEST().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 January 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST_PARTIAL_DIGEST_TEST:'
  write ( *, '(a)' ) '  FORTRAN90 version.'

  call i4_uniform_ab_test ( )
  call i4vec_distances_test ( )
  call i4vec_heap_d_test ( )
  call i4vec_print_test ( )
  call i4vec_sort_heap_a_test ( )
  call ksub_random_test ( )
  call example_partial_digest_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST_PARTIAL_DIGEST_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
