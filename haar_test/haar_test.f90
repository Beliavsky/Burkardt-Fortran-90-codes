program main

!*****************************************************************************80
!
!! haar_test() tests haar().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 March 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'haar_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test haar().'

  call test01 ( )
  call test02 ( )
  call test03 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'haar_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests HAAR_1D.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 March 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a_first
  real ( kind = rk ) a_last
  real ( kind = rk ) err
  integer i
  integer n
  real ( kind = rk ) r8vec_diff_norm
  real ( kind = rk ), allocatable :: u(:)
  real ( kind = rk ), allocatable :: v(:)
  real ( kind = rk ), allocatable :: w(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  HAAR_1D computes the Haar transform of a vector.'
!
!  Random data.
!
  n = 16
  allocate ( u(1:n) )
  allocate ( v(1:n) )
  allocate ( w(1:n) )
  call random_number ( harvest = u(1:n) )
  v(1:n) = u(1:n)

  call haar_1d ( n, v )

  w(1:n) = v(1:n)
  call haar_1d_inverse ( n, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U(i)        H(U)(i)  Hinv(H(U))(i)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i2,2x,f10.4,2x,f10.4,2x,f10.4)' ) &
      i, u(i), v(i), w(i)
  end do

  deallocate ( u )
  deallocate ( v )
  deallocate ( w )
!
!  Constant signal.
!
  n = 8
  allocate ( u(1:n) )
  allocate ( v(1:n) )
  allocate ( w(1:n) )
  u(1:n) = 1.0D+00
  v(1:n) = u(1:n)

  call haar_1d ( n, v )

  w(1:n) = v(1:n)
  call haar_1d_inverse ( n, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U(i)        H(U)(i)  Hinv(H(U))(i)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i2,2x,f10.4,2x,f10.4,2x,f10.4)' ) &
      i, u(i), v(i), w(i)
  end do

  deallocate ( u )
  deallocate ( v )
  deallocate ( w )
!
!  Linear signal.
!
  n = 16
  allocate ( u(1:n) )
  allocate ( v(1:n) )
  allocate ( w(1:n) )
  a_first = 1.0D+00
  a_last = dble ( n )
  call r8vec_linspace ( n, a_first, a_last, u )
  v(1:n) = u(1:n)

  call haar_1d ( n, v )

  w(1:n) = v(1:n)
  call haar_1d_inverse ( n, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U(i)        H(U)(i)  Hinv(H(U))(i)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i2,2x,f10.4,2x,f10.4,2x,f10.4)' ) &
      i, u(i), v(i), w(i)
  end do

  deallocate ( u )
  deallocate ( v )
  deallocate ( w )
!
!  Quadratic data.
!
  n = 8
  allocate ( u(1:n) )
  allocate ( v(1:n) )
  allocate ( w(1:n) )
  u(1) = 25.0D+00
  u(2) = 16.0D+00
  u(3) = 9.0D+00
  u(4) = 4.0D+00
  u(5) = 1.0D+00
  u(6) = 0.0D+00
  u(7) = 1.0D+00
  u(8) = 4.0D+00
  v(1:n) = u(1:n)

  call haar_1d ( n, v )

  w(1:n) = v(1:n)
  call haar_1d_inverse ( n, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   i      U(i)        H(U)(i)  Hinv(H(U))(i)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i2,2x,f10.4,2x,f10.4,2x,f10.4)' ) &
      i, u(i), v(i), w(i)
  end do

  deallocate ( u )
  deallocate ( v )
  deallocate ( w )
!
!  N not a power of 2.
!
  n = 99
  allocate ( u(1:n) )
  allocate ( v(1:n) )
  allocate ( w(1:n) )

  call random_number ( harvest = u(1:n) )

  v(1:n) = u(1:n)
  call haar_1d ( n, v )

  w(1:n) = v(1:n)
  call haar_1d_inverse ( n, w )

  err = r8vec_diff_norm ( n, u, w )

  write ( *, '(a)' ) ''
  write ( *, '(a,i4,a,g14.6)' ) &
    '  For N = ', n, &
    ', ||u-haar_1d_inverse(haar_1d(u))|| = ', err

  deallocate ( u )
  deallocate ( v )
  deallocate ( w )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests HAAR_2D and HAAR_2D_INVERSE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 March 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) err
  integer m
  integer n
  real ( kind = rk ) r8mat_diff_frobenius
  real ( kind = rk ), allocatable :: u(:,:)
  real ( kind = rk ), allocatable :: v(:,:)
  real ( kind = rk ), allocatable :: w(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  HAAR_2D computes the Haar transform of an array.'
  write ( *, '(a)' ) '  HAAR_2D_INVERSE inverts the transform.'
!
!  Demonstrate successful inversion.
!
  m = 16
  n = 4

  allocate ( u(1:m,1:n) )
  allocate ( v(1:m,1:n) )
  allocate ( w(1:m,1:n) )

  call random_number ( harvest = u(1:m,1:n) )

  call r8mat_print ( m, n, u, '  Input array U:' )

  v(1:m,1:n) = u(1:m,1:n)
  call haar_2d ( m, n, v )

  call r8mat_print ( m, n, v, '  Transformed array V:' )

  w(1:m,1:n) = v(1:m,1:n)
  call haar_2d_inverse ( m, n, w )

  call r8mat_print ( m, n, w, '  Recovered array W:' )

  deallocate ( u )
  deallocate ( v )
  deallocate ( w )
!
!  M, N not powers of 2.
!
  m = 37
  n = 53
  allocate ( u(1:m,1:n) )
  allocate ( v(1:m,1:n) )
  allocate ( w(1:m,1:n) )

  call random_number ( harvest = u(1:m,1:n) )

  v(1:m,1:n) = u(1:m,1:n)
  call haar_2d ( m, n, v )

  w(1:m,1:n) = v(1:m,1:n)
  call haar_2d_inverse ( m, n, w )

  err = r8mat_diff_frobenius ( m, n, u, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4,a,i4,a,g14.6)' ) &
    '  M = ', m, &
    ', N = ', n, &
    ', ||haar_2d_inverse(haar_2d(u))-u|| = ', err

  deallocate ( u )
  deallocate ( v )
  deallocate ( w )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests HAAR_2D and HAAR_2D_INVERSE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 May 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 128
  integer, parameter :: n = 128

  integer i
  real ( kind = rk ) u(m,n)
  real ( kind = rk ) v(m,n)
  real ( kind = rk ) w(m,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  HAAR_2D computes the Haar transform of an array.'
  write ( *, '(a)' ) '  HAAR_2D_INVERSE inverts the transform.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Apply this to a 128x128 matrix of 0 and 1 values'
  write ( *, '(a)' ) '  which is actually a bit map of the Sierpinski triangle.'
!
!  Demonstrate successful inversion.
!
  open ( unit = 1, file = 'sierpinski.txt', status = 'old' )
  do i = 1, n
    read ( 1, * ) u(i,1:n)
  end do
  close ( unit = 1 )

  call r8mat_print_some ( m, n, u, 1,1,10,10,'  Input array U:' )

  v(1:m,1:n) = u(1:m,1:n)
  call haar_2d ( m, n, v )

  call r8mat_print_some ( m, n, v, 1,1,10,10,'  Transformed array V:' )

  w(1:m,1:n) = v(1:m,1:n)
  call haar_2d_inverse ( m, n, w )
!
!  Merely for neatness, round the data (this gets rid of tiny nonzeros).
!
  w(1:m,1:n) = floor ( w(1:m,1:n) + 0.5 )

  call r8mat_print_some ( m, n, w, 1,1,10,10,'  Recovered array W:' )

  return
end
