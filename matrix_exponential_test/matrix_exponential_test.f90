program main

!*****************************************************************************80
!
!! MATRIX_EXPONENTIAL_TEST() tests MATRIX_EXPONENTIAL().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MATRIX_EXPONENTIAL_TEST:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test MATRIX_EXPONENTIAL.'
  write ( *, '(a)' ) '  The C8LIB and R8LIB libraries are needed.'
  write ( *, '(a)' ) '  This test needs the TEST_MATRIX_EXPONENTIAL library.'

  call matrix_exponential_test01 ( )
  call matrix_exponential_test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MATRIX_EXPONENTIAL_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine matrix_exponential_test01 ( )

!*****************************************************************************80
!
!! MATRIX_EXPONENTIAL_TEST01 compares matrix exponential algorithms.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: a(:,:)
  real ( kind = rk ), allocatable :: a_exp(:,:)
  integer n
  integer test
  integer test_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MATRIX_EXPONENTIAL_TEST01:'
  write ( *, '(a)' ) '  R8MAT_EXPM1 is based on MATLAB''s matrix exponential function;'
  write ( *, '(a)' ) '  R8MAT_EXPM2 uses a Taylor series approach;'

  call r8mat_exp_test_num ( test_num )

  do test = 1, test_num

    write ( *, '(a)' ) ' '
    write ( *, '(a,i4)' ) '  Test #', test

    call r8mat_exp_story ( test )

    call r8mat_exp_n ( test, n )

    write ( *, '(a,i4)' ) '  Matrix order N = ', n

    allocate ( a(1:n,1:n) )

    call r8mat_exp_a ( test, n, a )

    call r8mat_print ( n, n, a, '  Matrix:' )

    allocate ( a_exp(1:n,1:n) )

    call r8mat_expm1 ( n, a, a_exp )
    call r8mat_print ( n, n, a_exp, '  EXPM1(A):' )

    call r8mat_expm2 ( n, a, a_exp )
    call r8mat_print ( n, n, a_exp, '  EXPM2(A):' )

    call r8mat_exp_expa ( test, n, a_exp )
    call r8mat_print ( n, n, a_exp, '  Exact Exponential:' )

    deallocate ( a )
    deallocate ( a_exp )

  end do

  return
end
subroutine matrix_exponential_test02 ( )

!*****************************************************************************80
!
!! MATRIX_EXPONENTIAL_TEST02 compares matrix exponential algorithms.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 March 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  complex ( kind = ck ), allocatable :: a(:,:)
  complex ( kind = ck ), allocatable :: a_exp(:,:)
  integer n
  integer test
  integer test_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MATRIX_EXPONENTIAL_TEST02:'
  write ( *, '(a)' ) '  C8MAT_EXPM1 is based on MATLAB''s matrix exponential function;'
  write ( *, '(a)' ) '  C8MAT_EXPM2 uses a Taylor series approach;'

  call c8mat_exp_test_num ( test_num )

  do test = 1, test_num

    write ( *, '(a)' ) ' '
    write ( *, '(a,i4)' ) '  Test #', test

    call c8mat_exp_story ( test )

    call c8mat_exp_n ( test, n )

    write ( *, '(a,i4)' ) '  Matrix order N = ', n

    allocate ( a(1:n,1:n) )

    call c8mat_exp_a ( test, n, a )

    call c8mat_print ( n, n, a, '  Matrix:' )

    allocate ( a_exp(1:n,1:n) )

    call c8mat_expm1 ( n, a, a_exp )
    call c8mat_print ( n, n, a_exp, '  EXPM1(A):' )

    call c8mat_exp_expa ( test, n, a_exp )
    call c8mat_print ( n, n, a_exp, '  Exact Exponential:' )

    deallocate ( a )
    deallocate ( a_exp )

  end do

  return
end
