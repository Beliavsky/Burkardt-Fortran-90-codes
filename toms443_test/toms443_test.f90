program main

!*****************************************************************************80
!
!! toms443_test() tests toms443().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 June 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'toms443_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test TOMS443().'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS433_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests WEW_A
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 June 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) en
  integer n_data
  real ( kind = rk ) w1
  real ( kind = rk ) w2
  real ( kind = rk ) wew_a
  real ( kind = rk ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Test WEW_A to evaluate'
  write ( *, '(a)' ) '  Lambert''s W function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X      Exact    Computed      Error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call lambert_w_values ( n_data, x, w1 )

    if ( n_data <= 0 ) then
      exit
    end if

    if ( x == 0.0D+00 ) then
      w2 = 0.0D+00
    else
      w2 = wew_a ( x, en )
    end if

    write ( *, '(2x,f12.4,2x,g16.8,2x,g16.8,2x,e10.2)' ) &
      x, w1, w2, abs ( w1 - w2 )

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests WEW_B
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 June 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) en
  integer n_data
  real ( kind = rk ) w1
  real ( kind = rk ) w2
  real ( kind = rk ) wew_b
  real ( kind = rk ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Test WEW_B to evaluate'
  write ( *, '(a)' ) '  Lambert''s W function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X      Exact    Computed      Error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call lambert_w_values ( n_data, x, w1 )

    if ( n_data <= 0 ) then
      exit
    end if

    if ( x == 0.0D+00 ) then
      w2 = 0.0D+00
    else
      w2 = wew_b ( x, en )
    end if

    write ( *, '(2x,f12.4,2x,g16.8,2x,g16.8,2x,e10.2)' ) &
      x, w1, w2, abs ( w1 - w2 )

  end do

  return
end

