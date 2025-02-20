program main

!*****************************************************************************80
!
!! hankel_cholesky_test() tests hankel_cholesky().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 January 2017
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'hankel_cholesky_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test HANKEL_CHOLESKY().'

  call hankel_cholesky_upper_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'HANKEL_CHOLESKY_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine hankel_cholesky_upper_test ( )

!*****************************************************************************80
!
!! HANKEL_CHOLESKY_UPPER_TEST tests HANKEL_CHOLESKY_UPPER.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 January 2017
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  integer flag
  real ( kind = rk ) h(n,n)
  real ( kind = rk ) hanti(2*n-1)
  real ( kind = rk ) l(n,n)
  real ( kind = rk ) lii(n)
  real ( kind = rk ) liim1(n-1)
  real ( kind = rk ) r1(n,n)
  real ( kind = rk ) r2(n,n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'hankel_cholesky_upper_test():'
  write ( *, '(a)' ) '  HANKEL_CHOLESKY_UPPER() is given a Hankel matrix H and'
  write ( *, '(a)' ) '  computes an upper triangular matrix R such that'
  write ( *, '(a)' ) '  H = R'' * R'
!
!  Get a Hankel matrix that is positive definite.
!
  call random_number ( harvest = lii(1:n) )
  call random_number ( harvest = liim1(1:n-1) )
  call hankel_pds_cholesky_lower ( n, lii, liim1, l )
  h = matmul ( l, transpose ( l ) )
  call r8mat_print ( n, n, h, '  The Hankel matrix H:' )
!
!  Compute R using R8MAT_CHOLESKY_FACTOR_UPPER.
!
  call r8mat_cholesky_factor_upper ( n, h, r1, flag )
  if ( flag /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) ' R8MAT_CHOLESKY_FACTOR_UPPER says H is not positive definite.'
  else
    call r8mat_print ( n, n, r1, '  R computed by R8MAT_CHOLESKY_FACTOR_UPPER:' )
  end if
!
!  Compute R using HANKEL_CHOLESKY.
!
  hanti(1:n) = h(1:n,1)
  hanti(n+1:2*n-1) = h(n,2:n)

  call hankel_cholesky_upper ( n, hanti, r2 )
  call r8mat_print ( n, n, r2, '  R computed by HANKEL_CHOLESKY:' )

  return
end
