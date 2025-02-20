program main

!*****************************************************************************80
!
!! clausen_test() tests clausen().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'clausen_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test clausen().'

  call clausen_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'clausen_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine clausen_test ( )

!*****************************************************************************80
!
!! clausen_test() tests clausen().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) clausen
  real ( kind = rk ) diff
  real ( kind = rk ) fx1
  real ( kind = rk ) fx2
  integer n_data
  real ( kind = rk ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'clausen_test():'
  write ( *, '(a)' ) '  clausen() evaluates the Clausen function.'
  write ( *, '(a)' ) '  Compare its results to tabulated data.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '                               Tabulated               Computed'
  write ( *, '(a)' ) '           X                          FX                     FX        Diff'
  write ( *, '(a)' ) ''

  n_data = 0

  do

    call clausen_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = clausen ( x )

    diff = abs ( fx1 - fx2 )

    write ( *, '(2x,g14.6,2x,g24.16,2x,g24.16,2x,e7.1)' ) x, fx1, fx2, diff

  end do

  return
end
