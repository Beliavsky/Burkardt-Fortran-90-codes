program main

!*****************************************************************************80
!
!! asa241_test() tests asa241().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'asa241_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test asa241().'

  call test02 ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'asa241_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test02 ( )

!*****************************************************************************80
!
!! test02() tests r8_normal_01_cdf_inverse().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) fx
  integer n_data
  real ( kind = rk ) r8_normal_01_cdf_inverse
  real ( kind = rk ) x
  real ( kind = rk ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test02():'
  write ( *, '(a)' ) '  Let FX = NormalCDF ( X ).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  normal_01_cdf_values() returns some values of ( X, FX ).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  r8_normal_01_cdf_inverse() takes the value of FX, and'
  write ( *, '(a)' ) '  computes an estimate X2, of the corresponding input, '
  write ( *, '(a)' ) '  argument, accurate to about 16 decimal places.'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) &
    '      FX                         X                         X2' &
    // '                     DIFF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call normal_01_cdf_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    x2 = r8_normal_01_cdf_inverse ( fx )

    write ( *, '(2x,g24.16,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    fx, x, x2, abs ( x - x2 )

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

