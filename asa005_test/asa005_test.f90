program main

!*****************************************************************************80
!
!! asa005_test() tests asa005().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'asa005_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test asa005().'

  call prncst_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'asa005_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine prncst_test ( )

!*****************************************************************************80
!
!! prncst_test() tests prncst().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer df
  real ( kind = rk ) fx
  real ( kind = rk ) fx2
  integer ifault
  real ( kind = rk ) lambda
  integer n_data
  real ( kind = rk ) prncst
  real ( kind = rk ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'prncst_test():'
  write ( *, '(a)' ) '  prncst() computes the noncentral Student T '
  write ( *, '(a)' ) '  Cumulative Density Function.'
  write ( *, '(a)' ) '  Compare with tabulated values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) '      X   LAMBDA  DF     ', &
    ' CDF                       CDF                     DIFF'
  write ( *, '(a,a)' ) '                         ', &
    ' Tabulated                 PRNCST'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call student_noncentral_cdf_values ( n_data, df, lambda, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = prncst ( x, df, lambda, ifault )

    write ( *, &
    '(2x,f6.2,2x,f6.2,2x,i2,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    x, lambda, df, fx, fx2, abs ( fx - fx2 )

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
!    26 August 2021
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

