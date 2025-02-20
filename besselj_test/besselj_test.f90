program main

!*****************************************************************************80
!
!! besselj_test() tests besselj().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'besselj_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test besselj().'

  call dbesj_test ( )
  call rjbesl_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'besselj_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine dbesj_test ( )

!*****************************************************************************80
!
!! dbesj_test() tests dbesj().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) alpha
  real ( kind = rk ) fx
  integer n
  integer n_data
  integer nz
  real ( kind = rk ) x
  real ( kind = rk ) y(1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'dbesj_test():'
  write ( *, '(a)' ) '  dbesj() evaluates the Bessel J function with NONINTEGER'
  write ( *, '(a)' ) '  order.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        ALPHA           X                     FX                        FX'
  write ( *, '(a)' ) '                                              exact                     computed'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_jx_values ( n_data, alpha, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    n = 1
    call dbesj ( x, alpha, n, y, nz )

    write ( *, '(2x,f12.6,2x,f24.16,2x,g24.16,2x,g24.16)' ) &
      alpha, x, fx, y(1)

  end do

  return
end
subroutine rjbesl_test ( )

!*****************************************************************************80
!
!! rjbesl_test() tests rjbesl().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) alpha
  real ( kind = rk ), allocatable :: b(:)
  real ( kind = rk ) fx
  integer n
  integer n_data
  integer nb
  integer ncalc
  real ( kind = rk ) order
  real ( kind = rk ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'rjbesl_test():'
  write ( *, '(a)' ) '  rjbesl() evaluates the Bessel J function with NONINTEGER'
  write ( *, '(a)' ) '  order.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        ORDER           X                     FX                        FX'
  write ( *, '(a)' ) '                                              exact                     computed'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_jx_values ( n_data, order, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    n = int ( order )
    alpha = order - real ( n, kind = rk )
    allocate ( b(0:n) )
    nb = n + 1

    call rjbesl ( x, alpha, nb, b, ncalc )

    write ( *, '(2x,f12.6,2x,f24.16,2x,g24.16,2x,g24.16)' ) &
      order, x, fx, b(n)

    deallocate ( b )

  end do

  return
end
