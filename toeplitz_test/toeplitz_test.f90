program main

!*****************************************************************************80
!
!! toeplitz_test() tests toeplitz().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 October 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOEPLITZ_TEST():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test TOEPLITZ().'

  call test03 ( )
  call test04 ( )
  call test08 ( )
  call test09 ( )
  call test10 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOEPLITZ_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests C8CI_SL.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 October 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ck8 = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: n = 4

  complex ( kind = ck8 ) a(n)
  integer seed
  complex ( kind = ck8 ) x(n)
  complex ( kind = ck8 ) x2(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  C8CI_SL solves a complex circulant system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  seed = 123456789
  call c8ci_random ( n, seed, a )

  call c8ci_print ( n, a, '  The circulant matrix:' )
!
!  Set the desired solution.
!
  call c8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
  call c8ci_mxv ( n, a, x, x2 )
!
!  Solve the linear system.
!
  call c8ci_sl ( n, a, x2 )

  call c8vec_print ( n, x2, '  Solution:' )

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests C8TO_SL.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 October 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ck8 = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: n = 4

  complex ( kind = ck8 ) a(2*n-1)
  complex ( kind = ck8 ) b(n)
  integer job
  integer seed
  complex ( kind = ck8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  C8TO_SL solves a complex Toeplitz system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  seed = 123456789
  call c8to_random ( n, seed, a )

  call c8to_print ( n, a, '  The Toeplitz matrix:' )

  do job = 0, 1
!
!  Set the desired solution.
!
    call c8vec_indicator ( n, x )

    if ( job == 0 ) then
      call c8vec_print ( n, x, '  Desired solution:' )
    else
      call c8vec_print ( n, x, '  Desired solution to transposed system:' )
    end if
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call c8to_mxv ( n, a, x, b )
    else
      call c8to_vxm ( n, a, x, b )
    end if

    if ( job == 0 ) then
      call c8vec_print ( n, b, '  Right Hand Side:' )
    else
      call c8vec_print ( n, b, '  Right Hand Side of transposed system:' )
    end if
!
!  Solve the linear system.
!
    call c8to_sl ( n, a, b, x, job )

    if ( job == 0 ) then
      call c8vec_print ( n, x, '  Solution:' )
    else
      call c8vec_print ( n, x, '  Solution to transposed system:' )
    end if

  end do

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests R8TO_SL.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 October 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ck8 = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: n = 4

  real ( kind = rk8 ) a(2*n-1)
  real ( kind = rk8 ) b(n)
  integer job
  integer seed
  real ( kind = rk8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  R8TO_SL solves a real Toeplitz system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  seed = 123456789
  call r8to_random ( n, seed, a )

  call r8to_print ( n, a, '  The Toeplitz matrix:' )

  do job = 0, 1
!
!  Set the desired solution.
!
    call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call r8to_mxv ( n, a, x, b )
    else
      call r8to_vxm ( n, a, x, b )
    end if
!
!  Solve the linear system.
!
    call r8to_sl ( n, a, b, x, job )

    if ( job == 0 ) then
      call r8vec_print ( n, x, '  Solution:' )
    else
      call r8vec_print ( n, x, '  Solution to transposed system:' )
    end if

  end do

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests R8BTO_MXV, R8BTO_VXM.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 October 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ck8 = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: l = 3
  integer, parameter :: m = 2

  real ( kind = rk8 ), dimension ( m, m, l ) ::  a1 = reshape ( (/ &
    1.0D+00, 5.0D+00, 2.0D+00, 5.0D+00, &
    3.0D+00, 6.0D+00, 4.0D+00, 6.0D+00, &
    5.0D+00, 7.0D+00, 6.0D+00, 7.0D+00 /), (/ m, m, l /) )

  real ( kind = rk8 ), dimension ( m, m, l-1 ) :: a2 = reshape ( (/ &
    7.0D+00, 8.0D+00, 8.0D+00, 8.0D+00, &
    9.0D+00, 9.0D+00, 0.0D+00, 9.0D+00 /), (/ m, m, l-1 /) )
  real ( kind = rk8 ) b(m,l)
  real ( kind = rk8 ) x(m,l)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  For a real block Toeplitz matrix,'
  write ( *, '(a)' ) '  R8BTO_MXV computes A * x.'
  write ( *, '(a)' ) '  R8BTO_VXM computes x * A.'

  call r8bto_print ( m, l, a1, a2, '  The block Toeplitz matrix:' )

  call r8vec_indicator ( m*l, x )

  call r8vec_print ( m*l, x, '  The vector x:' )

  call r8bto_mxv ( m, l, a1, a2, x, b )

  call r8vec_print ( m*l, b, '  The product A*x:' )

  call r8bto_vxm ( m, l, a1, a2, x, b )

  call r8vec_print ( m*l, b, '  The product x*A:' )

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests R8BTO_SL.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 October 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ck8 = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: m = 2
  integer, parameter :: l = 3
  integer, parameter :: n = m * l

  real ( kind = rk8 ), dimension ( m, m, l ) ::  a1 = reshape ( (/ &
    9.0D+00, 2.0D+00, 1.0D+00, 8.0D+00, &
    3.0D+00, 6.0D+00, 4.0D+00, 6.0D+00, &
    5.0D+00, 7.0D+00, 6.0D+00, 7.0D+00 /), (/ m, m, l /) )

  real ( kind = rk8 ), dimension ( m, m, l-1 ) :: a2 = reshape ( (/ &
    7.0D+00, 8.0D+00, 8.0D+00, 8.0D+00, &
    9.0D+00, 9.0D+00, 0.0D+00, 9.0D+00 /), (/ m, m, l-1 /) )
  real ( kind = rk8 ) b(n)
  real ( kind = rk8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  R8BTO_SL solves a block Toeplitz system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8bto_print ( m, l, a1, a2, '  The block Toeplitz matrix:' )
!
!  Set the desired solution.
!
  call r8vec_indicator ( n, x )
!
!  Compute the right hand side.
!
  call r8bto_mxv ( m, l, a1, a2, x, b )

  call r8vec_print ( n, b, '  Right hand side:' )

  call r8bto_sl ( m, l, a1, a2, b, x )

  call r8vec_print ( n, x, '  Computed solution:' )

  return
end
