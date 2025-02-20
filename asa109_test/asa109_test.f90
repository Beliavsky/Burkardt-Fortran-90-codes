program main

!*****************************************************************************80
!
!! asa109_test() tests asa109().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'asa109_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test asa109().'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'asa109_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! test01() tests xinbta().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) beta_log
  real ( kind = rk ) fx
  integer ifault
  integer n_data
  real ( kind = rk ) x
  real ( kind = rk ) x2
  real ( kind = rk ) xinbta

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test01():'
  write ( *, '(a)' ) '  xinbta() inverts the incomplete Beta function.'
  write ( *, '(a)' ) '  Compare with tabulated values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) '      A       B           FX      ', &
    '    X                         X                       DIFF'
  write ( *, '(a,a)' ) '                                  ', &
    '   (tabulated)               (XINBTA)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call beta_inc_values ( n_data, a, b, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    beta_log = lgamma ( a ) &
             + lgamma ( b ) &
             - lgamma ( a + b )

    x2 = xinbta ( a, b, beta_log, fx, ifault )

    write ( *, &
    '(2x,f6.2,2x,f6.2,2x,f14.6,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    a, b, fx, x, x2, abs ( x - x2 )

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! test02() tests beta_inc_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) beta_log
  real ( kind = rk ) betain
  real ( kind = rk ) fx
  real ( kind = rk ) fx2
  integer ifault
  integer n_data
  real ( kind = rk ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test02():'
  write ( *, '(a)' ) '  beta_inc_values() returns values of '
  write ( *, '(a)' ) '  the incomplete Beta function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '          A               B               X           BETA_INC(A,B)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call beta_inc_values ( n_data, a, b, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    beta_log = lgamma ( a ) &
             + lgamma ( b ) &
             - lgamma ( a + b )

    fx2 = betain ( x, a, b, beta_log, ifault )

    write ( *, &
      '(2x,f12.8,2x,f12.8,2x,f14.6,2x,g24.16,2x,g24.16,2x,e10.4)' ) &
      a, b, x, fx, fx2, abs ( fx - fx2 )

  end do

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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end

