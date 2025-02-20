program main

!*****************************************************************************80
!
!! SIMPLEX_COORDINATES_TEST tests SIMPLEX_COORDINATES.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 September 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SIMPLEX_COORDINATES_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the SIMPLEX_COORDINATES library.'

  n = 3
  call simplex_coordinates1_test ( n )
  call simplex_coordinates2_test ( n )

  n = 4
  call simplex_coordinates1_test ( n )
  call simplex_coordinates2_test ( n )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SIMPLEX_COORDINATES_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine simplex_coordinates1_test ( n )

!*****************************************************************************80
!
!! SIMPLEX_COORDINATES1_TEST calls SIMPLEX_COORDINATES1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the spatial dimension.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) r8_factorial
  real ( kind = rk ) side
  real ( kind = rk ) volume
  real ( kind = rk ) volume2
  real ( kind = rk ) x(n,n+1)
  real ( kind = rk ) xtx(n+1,n+1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SIMPLEX_COORDINATES1_TEST'
  write ( *, '(a)' ) '  Call SIMPLEX_COORDINATES1'

  call simplex_coordinates1 ( n, x )

  call r8mat_transpose_print ( n, n + 1, x, '  Simplex vertex coordinates:' )

  side = sqrt ( sum ( ( x(1:n,1) - x(1:n,2) )**2 ) )

  call simplex_volume ( n, x, volume )

  volume2 = sqrt ( real ( n + 1, kind = rk ) ) / r8_factorial ( n ) &
    / sqrt ( 2.0D+00**n ) * side**n

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Side length =     ', side
  write ( *, '(a,g14.6)' ) '  Volume =          ', volume
  write ( *, '(a,g14.6)' ) '  Expected volume = ', volume2

  xtx = matmul ( transpose ( x ), x )

  call r8mat_transpose_print ( n + 1, n + 1, xtx, '  Dot product matrix:' )

  return
end
subroutine simplex_coordinates2_test ( n )

!*****************************************************************************80
!
!! SIMPLEX_COORDINATES2_TEST calls SIMPLEX_COORDINATES2.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the spatial dimension.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) r8_factorial
  real ( kind = rk ) side
  real ( kind = rk ) volume
  real ( kind = rk ) volume2
  real ( kind = rk ) x(n,n+1)
  real ( kind = rk ) xtx(n+1,n+1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SIMPLEX_COORDINATES2_TEST'
  write ( *, '(a)' ) '  Call SIMPLEX_COORDINATES2'

  call simplex_coordinates2 ( n, x )

  call r8mat_transpose_print ( n, n + 1, x, '  Simplex vertex coordinates:' )

  side = sqrt ( sum ( ( x(1:n,1) - x(1:n,2) )**2 ) )

  call simplex_volume ( n, x, volume )

  volume2 = sqrt ( real ( n + 1, kind = rk ) ) / r8_factorial ( n ) &
    / sqrt ( 2.0D+00**n ) * side**n

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Side length =     ', side
  write ( *, '(a,g14.6)' ) '  Volume =          ', volume
  write ( *, '(a,g14.6)' ) '  Expected volume = ', volume2

  xtx = matmul ( transpose ( x ), x )

  call r8mat_transpose_print ( n + 1, n + 1, xtx, '  Dot product matrix:' )

  return
end
