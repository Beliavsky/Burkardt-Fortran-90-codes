program main

!*****************************************************************************80
!
!! colored_noise_test() tests colored_noise().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 June 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) alpha
  integer i
  integer n
  real ( kind = rk ) q_d

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COLORED_NOISE_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test COLORED_NOISE().'

  call r8vec_sftf_test ( )

  n = 128
  q_d = 1.0D+00
  alpha = 0.00D+00

  do i = 0, 8
    alpha = 0.25D+00 * real ( i, kind = rk )
    call colored_noise_test01 ( n, q_d, alpha )
  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COLORED_NOISE_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, * ) ' '
  call timestamp ( )

  stop 0
end
