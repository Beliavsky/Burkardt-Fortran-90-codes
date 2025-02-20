program main

!*****************************************************************************80
!
!! toms291_test() tests toms291().
!
!  Modified:
!
!    28 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS291_TEST():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test TOMS291().'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS291_TEST():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 compares results from ALOGAM with tabulated data.
!
!  Modified:
!
!    07 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) alogam
  real ( kind = rk ) fx1
  real ( kind = rk ) fx2
  integer ifault
  integer n_data
  real ( kind = rk ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) &
  '  ALOGAM computes the logarithm of the Gamma function.'
  write ( *, '(a)' ) '  Compare against tabulated data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) '          X        ', &
    '   FX                        FX'
  write ( *, '(a,a)' ) '                   ', &
    'Tabulated                  ALOGAM'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gamma_log_values ( n_data, x, fx1 )

    if ( n_data .eq. 0 ) then
       exit
    end if

    fx2 = alogam ( x, ifault )

    write ( *, '(2x,f14.6,2x,g24.16,2x,g24.16)' ) x, fx1, fx2

  end do

  return
end
