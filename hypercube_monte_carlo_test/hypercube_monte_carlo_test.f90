program main

!*****************************************************************************80
!
!! HYPERCUBE_MONTE_CARLO_TEST() tests HYPERCUBE_MONTE_CARLO().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 January 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HYPERCUBE_MONTE_CARLO_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test HYPERCUBE_MONTE_CARLO().'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HYPERCUBE_MONTE_CARLO_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 estimates integrals over the unit hypercube in 3D.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 January 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 3

  integer e(m)
  integer :: e_test(m,10) = reshape ( (/ &
    0, 0, 0, &
    1, 0, 0, &
    0, 1, 0, &
    0, 0, 1, &
    2, 0, 0, &
    1, 1, 0, &
    1, 0, 1, &
    0, 2, 0, &
    0, 1, 1, &
    0, 0, 2 /), (/ m, 10 /) )
  real ( kind = rk ) hypercube01_volume
  integer j
  integer n
  real ( kind = rk ) result(10)
  real ( kind = rk ), allocatable :: value(:)
  real ( kind = rk ), allocatable :: x(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Use HYPERCUBE01_SAMPLE to estimate integrals '
  write ( *, '(a)' ) '  over the interior of the unit hypercube in 3D.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '         N' // &
    '        1' // &
    '               X' // &
    '               Y ' // &
    '              Z' // &
    '               X^2' // &
    '              XY' // &
    '             XZ' // &
    '              Y^2' // &
    '             YZ' // &
    '               Z^2'
  write ( *, '(a)' ) ' '

  n = 1

  do while ( n <= 65536 )

    allocate ( value(1:n) )
    allocate ( x(1:m,1:n) )

    call hypercube01_sample ( m, n, x )

    do j = 1, 10

      e(1:m) = e_test(1:m,j)

      call monomial_value ( m, n, e, x, value )

      result(j) = hypercube01_volume ( m ) * sum ( value(1:n) ) &
        / real ( n, kind = rk )

    end do

    write ( *, '(2x,i8,10(2x,g14.6))' ) n, result(1:10)

    deallocate ( value )
    deallocate ( x )

    n = 2 * n

  end do

  write ( *, '(a)' ) ' '

  do j = 1, 10

    e(1:m) = e_test(1:m,j)

    call hypercube01_monomial_integral ( m, e, result(j) )

  end do

  write ( *, '(2x,a8,10(2x,g14.6))' ) '   Exact', result(1:10)

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 estimates integrals over the unit hypercube in 6D.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 January 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 6


  integer e(m)
  integer :: e_test(m,7) = reshape ( (/ &
    0, 0, 0, 0, 0, 0, &
    1, 0, 0, 0, 0, 0, &
    0, 2, 0, 0, 0, 0, &
    0, 2, 2, 0, 0, 0, &
    0, 0, 0, 4, 0, 0, &
    2, 0, 0, 0, 2, 2, &
    0, 0, 0, 0, 0, 6 /), (/ m, 7 /) )
  real ( kind = rk ) hypercube01_volume
  integer j
  integer n
  real ( kind = rk ) result(7)
  real ( kind = rk ), allocatable :: value(:)
  real ( kind = rk ), allocatable :: x(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Use HYPERCUBE01_SAMPLE to estimate integrals '
  write ( *, '(a)' ) '  over the interior of the unit hypercube in 6D.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '         N' // &
    '        1      ' // &
    '        U      ' // &
    '         V^2   ' // &
    '         V^2W^2' // &
    '         X^4   ' // &
    '         Y^2Z^2' // &
    '         Z^6'
  write ( *, '(a)' ) ' '

  n = 1

  do while ( n <= 65536 )

    allocate ( value(1:n) )
    allocate ( x(1:m,1:n) )

    call hypercube01_sample ( m, n, x )

    do j = 1, 7

      e(1:m) = e_test(1:m,j)

      call monomial_value ( m, n, e, x, value )

      result(j) = hypercube01_volume ( m ) * sum ( value(1:n) ) &
        / real ( n, kind = rk )

    end do

    write ( *, '(2x,i8,7(2x,g14.6))' ) n, result(1:7)

    deallocate ( value )
    deallocate ( x )

    n = 2 * n

  end do

  write ( *, '(a)' ) ' '

  do j = 1, 7

    e(1:m) = e_test(1:m,j)

    call hypercube01_monomial_integral ( m, e, result(j) )

  end do

  write ( *, '(2x,a8,7(2x,g14.6))' ) '   Exact', result(1:7)

  return
end
