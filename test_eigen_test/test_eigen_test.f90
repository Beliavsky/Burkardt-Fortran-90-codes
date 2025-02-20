program main

!*****************************************************************************80
!
!! test_eigen_test() tests test_eigen().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 June 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test_eigen_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test test_eigen().'

  call r8symm_gen_test ( )
  call r8nsymm_gen_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test_eigen_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine r8symm_gen_test ( )

!*****************************************************************************80
!
!! r8symm_gen_test() tests r8symm_gen().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 June 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk8 ) a(n,n)
  real ( kind = rk8 ) aq(n,n)
  integer j
  real ( kind = rk8 ) lambda(n)
  real ( kind = rk8 ) lambda2(n)
  real ( kind = rk8 ), parameter :: lambda_dev = 5.0D+00
  real ( kind = rk8 ), parameter :: lambda_mean = 10.0D+00
  real ( kind = rk8 ) q(n,n)
  integer :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8symm_gen_test():'
  write ( *, '(a)' ) '  r8symm_gen() generates an arbitrary size symmetric matrix'
  write ( *, '(a)' ) '  with known eigenvalues and eigenvectors.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Real data is declared as "REAL ( kind = rk8 )".'

  call r8symm_gen ( n, lambda_mean, lambda_dev, seed, a, q, lambda )

  call r8mat_print ( n, n, A, '  The symmetric matrix A:' )
  call r8mat_print ( n, n, Q, '  The eigenvector matrix Q:' )

  aq(1:n,1:n) = matmul ( a(1:n,1:n), q(1:n,1:n) )

  do j = 1, n
    lambda2(j) = sqrt ( sum ( aq(1:n,j)**2 ) )
  end do

  call r8vec2_print ( n, lambda, lambda2, &
    '  LAMBDA versus the column norms of A*Q:' )

  return
end
subroutine r8nsymm_gen_test ( )

!*****************************************************************************80
!
!! r8nsymm_gen_test() tests r8nsymm_gen().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 June 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk8 ) A(n,n)
  integer i
  real ( kind = rk8 ) lambda(n)
  real ( kind = rk8 ), parameter :: lambda_dev = 5.0D+00
  real ( kind = rk8 ), parameter :: lambda_mean = 10.0D+00
  real ( kind = rk8 ) Q(n,n)
  integer :: seed = 123456789
  real ( kind = rk8 ) T(n,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8nsymm_gen_test():'
  write ( *, '(a)' ) '  r8nsymm_gen() generates an arbitrary size nonsymmetric'
  write ( *, '(a)' ) '  matrix with known eigenvalues and eigenvectors.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Real data is declared as "REAL ( kind = rk8 )".'

  call r8nsymm_gen ( n, lambda_mean, lambda_dev, seed, A, Q, T )

  do i = 1, n
    lambda(i) = T(i,i)
  end do

  call r8vec_sort_bubble_a ( n, lambda )

  call r8mat_print ( n, n, A, '  The nonsymmetric matrix A:' )
  call r8mat_print ( n, n, Q, '  The eigenvector matrix Q:' )
  call r8mat_print ( n, n, T, '  The upper triangular matrix T:' )
  call r8vec_print ( n, lambda, '  The sorted eigenvalues LAMBDA:' )

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! timestamp() prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = 8 ) ampm
  integer d
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  integer values(8)
  integer y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
