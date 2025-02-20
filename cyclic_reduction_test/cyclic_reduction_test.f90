program main

!*****************************************************************************80
!
!! cyclic_reduction_test() tests cyclic_reduction().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    29 May 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CYCLIC_REDUCTION_TEST():'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test CYCLIC_REDUCTION().'

  call test01 ( )
  call test02 ( )
  call test03 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CYCLIC_REDUCTION_TEST():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end 
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests C83_CR_FA, C83_CR_SL.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 May 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 10

  complex ( kind = ck ) a(3,n)
  complex ( kind = ck ) a_cr(3,0:2*n)
  complex ( kind = ck ) b(n)
  integer i
  integer j
  real ( kind = rk ) t
  complex ( kind = ck ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  C83_CR_FA factors a complex tridiagonal matrix;'
  write ( *, '(a)' ) '  C83_CR_SL solves a factored system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix values.
!
  a(1,1) = 0.0D+00
  do j = 2, n
    t = real ( - j + 1, kind = rk )
    a(1,j) = cmplx ( -1.0D+00, t, kind = ck )
  end do

  do j = 1, n
    t = 2.0D+00 * real ( j, kind = rk )
    a(2,j) = cmplx ( 2.0D+00, t, kind = ck )
  end do

  do j = 1, n - 1
    t = real ( - j - 1, kind = rk )
    a(3,j) = cmplx ( -1.0D+00, t, kind = ck )
  end do
  a(3,n) = 0.0D+00
!
!  Set the desired solution.
!
  do i = 1, n
    t = real ( i, kind = rk )
    x(i) = cmplx ( t, t, kind = ck )
  end do
!
!  Compute the corresponding right hand side.
!
  call c83_mxv ( n, a, x, b )
!
!  Factor the matrix.
!
  call c83_cr_fa ( n, a, a_cr )
!
!  Solve the linear system.
!
  call c83_cr_sl ( n, a_cr, b, x )

  call c8vec_print ( n, x, '  Solution:' )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests R83_CR_FA, R83_CR_SLS.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    29 May 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5
  integer, parameter :: nb = 2

  real ( kind = rk ) a(3,n)
  real ( kind = rk ) a_cr(3,0:2*n)
  real ( kind = rk ) b(n,nb)
  logical, parameter :: debug = .true.
  real ( kind = rk ) x(n,nb)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  R83_CR_FA factors a real tridiagonal matrix;'
  write ( *, '(a)' ) '  R83_CR_SLS solves 1 or more systems.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
  write ( *, '(a)' ) '  Demonstrate multiple system solution method.'
!
!  Set the matrix.
!
  a(1,1) = 0.0D+00
  a(1,2:n) = -1.0D+00

  a(2,1:n) = 2.0D+00

  a(3,1:n-1) = -1.0D+00
  a(3,n) = 0.0D+00

  if ( debug ) then
    call r83_print ( n, a, '  Input matrix:' )
  end if
!
!  Factor the matrix once.
!
  call r83_cr_fa ( n, a, a_cr )

  if ( debug ) then
    call r83_print ( 2 * n + 1, a_cr, '  Cyclic reduction factor information:' )
  end if
!
!  Solve 2 systems simultaneously.
!
  b(1:n-1,1) = 0.0D+00
  b(n,1) = real ( n + 1, kind = rk )

  b(1,2) = 1.0D+00
  b(2:n-1,2) = 0.0D+00
  b(n,2) = 1.0D+00
!
!  Solve the linear systems.
!
  call r83_cr_sls ( n, a_cr, nb, b, x )

  call r8mat_print_some ( n, nb, x, 1, 1, 10, nb, '  Solutions:' )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests R83_CR_FA, R83_CR_SL.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    29 May 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 10

  real ( kind = rk ) a(3,n)
  real ( kind = rk ) a_cr(3,0:2*n)
  real ( kind = rk ) b(n)
  logical, parameter :: debug = .false.
  integer j
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  For a real tridiagonal matrix,'
  write ( *, '(a)' ) '  using CYCLIC REDUCTION,'
  write ( *, '(a)' ) '  R83_CR_FA factors;'
  write ( *, '(a)' ) '  R83_CR_SL solves 1 system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
  write ( *, '(a)' ) '  The matrix is NOT symmetric.'
!
!  Set the matrix values.
!
  a(1,1) = 0.0D+00
  do j = 2, n
    a(1,j) = real ( j, kind = rk )
  end do

  do j = 1, n
    a(2,j) = 4.0D+00 * real ( j, kind = rk )
  end do

  do j = 1, n - 1
    a(3,j) = real ( j, kind = rk )
  end do
  a(3,n) = 0.0D+00

  if ( debug ) then
    call r83_print ( n, a, '  The matrix:' )
  end if
!
!  Set the desired solution.
!
  call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r83_mxv ( n, a, x, b )

  if ( debug ) then
    call r8vec_print ( n, b, '  The right hand side:' )
  end if

  x(1:n) = 0.0D+00
!
!  Factor the matrix.
!
  call r83_cr_fa ( n, a, a_cr )

  if ( debug ) then
    call r83_print ( 2 * n + 1, a_cr, '  The factor information:' )
  end if
!
!  Solve the linear system.
!
  call r83_cr_sl ( n, a_cr, b, x )

  call r8vec_print ( n, x, '  The solution:' )

  return
end
