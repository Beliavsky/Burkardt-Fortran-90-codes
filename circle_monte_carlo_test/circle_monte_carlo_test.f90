program main

!*****************************************************************************80
!
!! circle_monte_carlo_test() tests circle_monte_carlo().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'circle_monte_carlo_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test circle_monte_carlo().'

  call circle01_sample_random_test ( )
  call circle01_sample_ergodic_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'circle_monte_carlo_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine circle01_sample_random_test ( )

!*****************************************************************************80
!
!! circle01_sample_random_test() tests circle01_sample_random().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) circle01_length
  integer e(2)
  integer :: e_test(2,7) = reshape ( (/ &
    0, 0, &
    2, 0, &
    0, 2, &
    4, 0, &
    2, 2, &
    0, 4, &
    6, 0 /), (/ 2, 7 /) )
  integer j
  integer n
  real ( kind = rk ) result(7)
  real ( kind = rk ), allocatable :: value(:)
  real ( kind = rk ), allocatable :: x(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'circle01_sample_random_test():'
  write ( *, '(a)' ) '  circle01_sample_random() randomly samples the unit circle.'
  write ( *, '(a)' ) '  Use it to estimate integrals.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N        1              X^2             Y^2' // &
    '             X^4           X^2Y^2             Y^4           X^6'
  write ( *, '(a)' ) ' '

  n = 1

  do while ( n <= 65536 )

    allocate ( value(1:n) )
    allocate ( x(1:2,1:n) )

    call circle01_sample_random ( n, x )

    do j = 1, 7

      e(1:2) = e_test(1:2,j)

      call monomial_value ( 2, n, e, x, value )

      result(j) = circle01_length ( ) * sum ( value(1:n) ) &
        / real ( n, kind = rk )

    end do

    write ( *, '(2x,i8,7(2x,g14.6))' ) n, result(1:7)

    deallocate ( value )
    deallocate ( x )

    n = 2 * n

  end do

  write ( *, '(a)' ) ' '

  do j = 1, 7

    e(1:2) = e_test(1:2,j)

    call circle01_monomial_integral ( e, result(j) )

  end do

  write ( *, '(2x,a8,7(2x,g14.6))' ) '   Exact', result(1:7)

  return
end
subroutine circle01_sample_ergodic_test ( )

!*****************************************************************************80
!
!! circle01_sample_ergodic_test() tests circle01_sample_ergodic().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) angle
  real ( kind = rk ) circle01_length
  integer e(2)
  integer :: e_test(2,7) = reshape ( (/ &
    0, 0, &
    2, 0, &
    0, 2, &
    4, 0, &
    2, 2, &
    0, 4, &
    6, 0 /), (/ 2, 7 /) )
  integer j
  integer n
  real ( kind = rk ) result(7)
  real ( kind = rk ), allocatable :: value(:)
  real ( kind = rk ), allocatable :: x(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'circle01_sample_ergodic_test():'
  write ( *, '(a)' ) '  circle01_sample_ergodic() ergodically samples the unit circle'
  write ( *, '(a)' ) '  Use it to estimate integrals.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N        1              X^2             Y^2' // &
    '             X^4           X^2Y^2             Y^4           X^6'
  write ( *, '(a)' ) ' '

  n = 1

  do while ( n <= 65536 )

    allocate ( value(1:n) )
    allocate ( x(1:2,1:n) )

    angle = 0.0D+00
    call circle01_sample_ergodic ( n, angle, x )

    do j = 1, 7

      e(1:2) = e_test(1:2,j)

      call monomial_value ( 2, n, e, x, value )

      result(j) = circle01_length ( ) * sum ( value(1:n) ) &
        / real ( n, kind = rk )

    end do

    write ( *, '(2x,i8,7(2x,g14.6))' ) n, result(1:7)

    deallocate ( value )
    deallocate ( x )

    n = 2 * n

  end do

  write ( *, '(a)' ) ' '

  do j = 1, 7

    e(1:2) = e_test(1:2,j)

    call circle01_monomial_integral ( e, result(j) )

  end do

  write ( *, '(2x,a8,7(2x,g14.6))' ) '   Exact', result(1:7)

  return
end
 
