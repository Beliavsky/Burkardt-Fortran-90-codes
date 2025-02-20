program main

!*****************************************************************************80
!
!! asa310_test() tests asa310().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA310_TEST():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test ASA310().'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA310_TEST():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests NCBETA against tabulated values.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 February
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) errmax
  real ( kind = rk ) fx
  real ( kind = rk ) fx2
  integer ifault
  real ( kind = rk ) lambda
  integer n_data
  real ( kind = rk ) ncbeta
  real ( kind = rk ) x

  errmax = 1.0D-10

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  Compare tabulated values of the noncentral'
  write ( *, '(a)' ) '  incomplete Beta Function against values'
  write ( *, '(a)' ) '  computed by NCBETA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      A        B     LAMBDA        X       ', &
    ' CDF               CDF          DIFF'
  write ( *, '(a,a)' ) &
    '                                           ', &
    '(tabulated)       (NCBETA)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call beta_noncentral_cdf_values ( n_data, a, b, lambda, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = ncbeta ( a, b, lambda, x, errmax, ifault )

    write ( *, '(2x,f7.1,2x,f7.1,2x,f7.1,2x,f10.4, 2x,g14.6,2x,g14.6,2x,g10.4)' ) &
    a, b, lambda, x, fx, fx2, abs ( fx - fx2 )

  end do

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
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

