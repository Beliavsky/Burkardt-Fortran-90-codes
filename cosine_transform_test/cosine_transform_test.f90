program main

!*****************************************************************************80
!
!! cosine_transform_test() tests cosine_transform().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 August 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'COSINE_TRANSFORM_TEST():'
  write ( *, '(a)' ) '  Fortran90 version.'
  write ( *, '(a)' ) '  Test COSINE_TRANSFORM().'

  call cosine_transform_test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'COSINE_TRANSFORM_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine cosine_transform_test01 ( )

!*****************************************************************************80
!
!! COSINE_TRANSFORM_TEST01 does a DCT and its inverse.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 August 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 10

  integer i
  real ( kind = rk ) r(n)
  real ( kind = rk ) s(n)
  real ( kind = rk ) t(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'COSINE_TRANSFORM_TEST01:'
  write ( *, '(a)' ) '  COSINE_TRANSFORM_DATA does a cosine transform of data'
  write ( *, '(a)' ) '  defined by a vector.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Apply the transform, then its inverse.'
  write ( *, '(a)' ) '  Let R be a random N vector.'
  write ( *, '(a)' ) '  Let S be the transform of D.'
  write ( *, '(a)' ) '  Let T be the transform of E.'
  write ( *, '(a)' ) '  Then R and T will be equal.'

  call random_number ( harvest = r(1:n) )
  call cosine_transform_data ( n, r, s )
  call cosine_transform_inverse ( n, s, t )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '     I      R(I)        S(I)        T(I)'
  write ( *, '(a)' ) ''

  do i = 1, n
    write ( *, '(2x,i4,2x,f10.6,2x,f10.6,2x,f10.6)' ) i, r(i), s(i), t(i)
  end do

  return
end
