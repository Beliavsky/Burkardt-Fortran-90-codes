program main

!*****************************************************************************80
!
!! MAIN is the main program for TOEPLITZ_CHOLESKY_TEST.
!
!  Discussion:
!
!    TOEPLITZ_CHOLESKY_TEST tests the TOEPLITZ_CHOLESKY library.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 January 2017
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TOEPLITZ_CHOLESKY_TEST:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TOEPLITZ_CHOLESKY library.'

  call t_cholesky_lower_test ( )
  call toep_cholesky_lower_test ( )
  call toeplitz_cholesky_lower_test ( )

  call t_cholesky_upper_test ( )
  call toep_cholesky_upper_test ( )
  call toeplitz_cholesky_upper_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TOEPLITZ_CHOLESKY_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine t_cholesky_lower_test ( )

!*****************************************************************************80
!
!! T_CHOLESKY_LOWER_TEST tests T_CHOLESKY_LOWER.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 January 2017
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 3

  real ( kind = rk ) b(n,n)
  real ( kind = rk ) l(n,n)
  real ( kind = rk ) t(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'T_CHOLESKY_LOWER_TEST'
  write ( *, '(a)' ) '  T_CHOLESKY_LOWER produces the lower Cholesky'
  write ( *, '(a)' ) '  factor L of a positive definite symmetric'
  write ( *, '(a)' ) '  Toeplitz matrix T, so that T=L*L''.'
  write ( *, '(a)' ) '  The first row of T is input..'

  t = (/ 1.0D+00, 0.5D+00, -0.375D+00 /)

  call r8vec_print ( n, t, '  First row of Toeplitz matrix T:' )

  call t_cholesky_lower ( n, t, l )
  call r8mat_print ( n, n, l, '  Computed lower Cholesky factor L:' )

  b = matmul ( l, transpose ( l ) )
  call r8mat_print ( n, n, b, '  Product LL'':' )

  return
end
subroutine toep_cholesky_lower_test ( )

!*****************************************************************************80
!
!! TOEP_CHOLESKY_LOWER_TEST tests TOEP_CHOLESKY_LOWER.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 January 2017
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 3

  real ( kind = rk ) b(n,n)
  real ( kind = rk ) g(2,n)
  real ( kind = rk ) l(n,n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TOEP_CHOLESKY_LOWER_TEST'
  write ( *, '(a)' ) '  TOEP_CHOLESKY_LOWER produces the lower Cholesky'
  write ( *, '(a)' ) '  factor L of a positive definite symmetric'
  write ( *, '(a)' ) '  Toeplitz matrix T, so that T=L*L''.'
  write ( *, '(a)' ) '  T is input in a compressed (2,N) array.'

  g = reshape ( (/ &
    1.0D+00, 0.0D+00,  &
    0.5D+00, 0.5D+00,  &
   -0.375D+00, -0.375D+00 /), (/ 2, 3 /) )

  call r8mat_print ( 2, n, g, '  Compressed Toeplitz matrix G:' )

  call toep_cholesky_lower ( n, g, l )
  call r8mat_print ( n, n, l, '  Computed lower Cholesky factor L:' )

  b = matmul ( l, transpose ( l ) )
  call r8mat_print ( n, n, b, '  Product LL'':' )

  return
end
subroutine toeplitz_cholesky_lower_test ( )

!*****************************************************************************80
!
!! TOEPLITZ_CHOLESKY_LOWER_TEST tests TOEPLITZ_CHOLESKY_LOWER.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 January 2017
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 3

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) b(n,n)
  real ( kind = rk ) l(n,n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TOEPLITZ_CHOLESKY_LOWER_TEST'
  write ( *, '(a)' ) '  TOEPLITZ_CHOLESKY_LOWER produces the lower Cholesky'
  write ( *, '(a)' ) '  factor L of a positive definite symmetric'
  write ( *, '(a)' ) '  Toeplitz matrix T, so that T=L*L''.'
  write ( *, '(a)' ) '  T is input as an NxN array.'

  a = reshape ( (/ &
        1.0D+00,   0.5D+00, -0.375D+00, &
        0.5D+00,   1.0D+00,  0.5D+00, &
       -0.375D+00, 0.5D+00,  1.0D+00 /), (/ 3, 3 /) )

  call r8mat_print ( n, n, a, '  Toeplitz matrix A:' )

  call toeplitz_cholesky_lower ( n, a, l )
  call r8mat_print ( n, n, l, '  Computed lower Cholesky factor L:' )

  b = matmul ( l, transpose ( l ) )
  call r8mat_print ( n, n, b, '  Product LL'':' )

  return
end
subroutine t_cholesky_upper_test ( )

!*****************************************************************************80
!
!! T_CHOLESKY_UPPER_TEST tests T_CHOLESKY_UPPER.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 January 2017
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 3

  real ( kind = rk ) b(n,n)
  real ( kind = rk ) r(n,n)
  real ( kind = rk ) t(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'T_CHOLESKY_UPPER_TEST'
  write ( *, '(a)' ) '  T_CHOLESKY_UPPER produces the upper Cholesky'
  write ( *, '(a)' ) '  factor R of a positive definite symmetric'
  write ( *, '(a)' ) '  Toeplitz matrix T, so that T=R''*R.'
  write ( *, '(a)' ) '  The first row of T is input.'
  
  t = (/  1.0D+00, 0.5D+00, -0.375D+00 /)

  call r8vec_print ( n, t, '  First row of Toeplitz matrix T:' )

  call t_cholesky_upper ( n, t, r )
  call r8mat_print ( n, n, r, '  Computed upper Cholesky factor R:' )

  b = matmul ( transpose ( r ), r )
  call r8mat_print ( n, n, b, '  Product R''R:' )

  return
end
subroutine toep_cholesky_upper_test ( )

!*****************************************************************************80
!
!! TOEP_CHOLESKY_UPPER_TEST tests TOEP_CHOLESKY_UPPER.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 January 2017
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 3

  real ( kind = rk ) b(n,n)
  real ( kind = rk ) g(2,n)
  real ( kind = rk ) r(n,n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TOEP_CHOLESKY_UPPER_TEST'
  write ( *, '(a)' ) '  TOEP_CHOLESKY_UPPER produces the upper Cholesky'
  write ( *, '(a)' ) '  factor R of a positive definite symmetric'
  write ( *, '(a)' ) '  Toeplitz matrix T, so that T=R''*R.'
  write ( *, '(a)' ) '  T is input in a compressed (2,N) array.'
  
  g = reshape ( (/ &
    1.0D+00,    0.0D+00,  &
    0.5D+00,    0.5D+00,  &
   -0.375D+00, -0.375D+00 /), (/ 2, 3 /) )

  call r8mat_print ( 2, n, g, '  Compressed Toeplitz matrix G:' )

  call toep_cholesky_upper ( n, g, r )
  call r8mat_print ( n, n, r, '  Computed upper Cholesky factor R:' )

  b = matmul ( transpose ( r ), r )
  call r8mat_print ( n, n, b, '  Product R''R:' )

  return
end
subroutine toeplitz_cholesky_upper_test ( )

!*****************************************************************************80
!
!! TOEPLITZ_CHOLESKY_UPPER_TEST tests TOEPLITZ_CHOLESKY_UPPER.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 January 2017
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 3

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) b(n,n)
  real ( kind = rk ) r(n,n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TOEPLITZ_CHOLESKY_UPPER_TEST'
  write ( *, '(a)' ) '  TOEPLITZ_CHOLESKY_UPPER produces the upper Cholesky'
  write ( *, '(a)' ) '  factor R of a positive definite symmetric'
  write ( *, '(a)' ) '  Toeplitz matrix T, so that T=R''*R.'
  write ( *, '(a)' ) '  T is input as an NxN array.'

  a = reshape ( (/ &
        1.0D+00,   0.5D+00, -0.375D+00, &
        0.5D+00,   1.0D+00,  0.5D+00, &
       -0.375D+00, 0.5D+00,  1.0D+00 /), (/ 3, 3 /) )

  call r8mat_print ( n, n, a, '  Toeplitz matrix A:' )

  call toeplitz_cholesky_upper ( n, a, r )
  call r8mat_print ( n, n, r, '  Computed upper Cholesky factor R:' )

  b = matmul ( transpose ( r ), r )
  call r8mat_print ( n, n, b, '  Product R''R:' )

  return
end




