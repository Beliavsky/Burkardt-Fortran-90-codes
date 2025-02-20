program main

!*****************************************************************************80
!
!! toms179_test() tests toms179().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 January 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS179_TEST():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test TOMS179().'

  call alogam_test ( )
  call mdbeta_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS179_TEST():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine alogam_test ( )

!*****************************************************************************80
!
!! ALOGAM_TEST tests ALOGAM.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 January 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) alogam
  real ( kind = rk ) fx
  real ( kind = rk ) fx2
  integer ifault
  integer n_data
  real ( kind = rk ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ALOGAM_TEST'
  write ( *, '(a)' ) '  ALOGAM estimates the logarithm of the Gamma function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) '      X         Exact Value               ', &
    'Computed                Diff'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gamma_log_values ( n_data, x, fx )

    if ( n_data <= 0 ) then
      exit
    end if

    fx2 = alogam ( x, ifault )

    write ( *, '(2x,f8.4,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
      x, fx, fx2, abs ( fx - fx2 )

  end do

  return
end
subroutine mdbeta_test ( )

!*****************************************************************************80
!
!! MDBETA_TEST tests MDBETA.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 January 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) fx
  real ( kind = rk ) fx2
  integer ier
  integer n_data
  real ( kind = rk ) p
  real ( kind = rk ) q
  real ( kind = rk ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MDBETA_TEST'
  write ( *, '(a)' ) '  MDBETA_TEST estimates the modified Beta function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) '      X         P         Q         ', &
    'Exact Value               Computed                Diff'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call beta_cdf_values ( n_data, p, q, x, fx )

    if ( n_data <= 0 ) then
      exit
    end if

    call mdbeta ( x, p, q, fx2, ier )

    write ( *, &
      '(2x,f8.4,2x,f8.4,2x,f8.4,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
      x, p, q, fx, fx2, abs ( fx - fx2 )

  end do

  return
end

