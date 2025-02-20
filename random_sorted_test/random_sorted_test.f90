program main

!*****************************************************************************80
!
!! random_sorted_test() tests random_sorted().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 March 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'random_sorted_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test random_sorted().'

  call r8vec_normal_01_sorted_test ( )
  call r8vec_uniform_01_sorted1_test ( )
  call r8vec_uniform_01_sorted2_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'random_sorted_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine r8vec_normal_01_sorted_test ( )

!*****************************************************************************80
!
!! R8VEC_NORMAL_01_SORTED_TEST tests R8VEC_NORMAL_01_SORTED_TEST,
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 March 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  integer i
  real ( kind = rk ) r8vec(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8VEC_NORMAL_01_SORTED_TEST:'
  write ( *, '(a)' ) '  R8VEC_NORMAL_01_SORTED generates a vector of N normal 01'
  write ( *, '(a)' ) '  random values in ascending sorted order.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Generate several examples:'
  write ( *, '(a)' ) ''

  do i = 1, 10
    call r8vec_normal_01_sorted ( n, r8vec )
    call r8vec_transpose_print ( n, r8vec, '  R8VEC:' )
  end do

  return
end
subroutine r8vec_uniform_01_sorted1_test ( )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01_SORTED1_TEST tests R8VEC_UNIFORM_01_SORTED1_TEST,
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 March 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  integer i
  real ( kind = rk ) r8vec(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8VEC_UNIFORM_01_SORTED1_TEST:'
  write ( *, '(a)' ) '  R8VEC_UNIFORM_01_SORTED1 generates a vector of N random'
  write ( *, '(a)' ) '  values in ascending sorted order.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Generate several examples:'
  write ( *, '(a)' ) ''

  do i = 1, 10
    call r8vec_uniform_01_sorted1 ( n, r8vec )
    call r8vec_transpose_print ( n, r8vec, '  R8VEC:' )
  end do

  return
end
subroutine r8vec_uniform_01_sorted2_test ( )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01_SORTED2_TEST tests R8VEC_UNIFORM_01_SORTED2_TEST.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 March 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  integer i
  real ( kind = rk ) r8vec(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8VEC_UNIFORM_01_SORTED2_TEST:'
  write ( *, '(a)' ) '  R8VEC_UNIFORM_01_SORTED2 generates a vector of N random'
  write ( *, '(a)' ) '  values in ascending sorted order.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Generate several examples:'
  write ( *, '(a)' ) ''

  do i = 1, 10
    call r8vec_uniform_01_sorted2 ( n, r8vec )
    call r8vec_transpose_print ( n, r8vec, '  R8VEC:' )
  end do

  return
end
