program main

!*****************************************************************************80
!
!! asa076_test() tests asa076().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'asa076_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test asa076().'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'asa076_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! test01() tests tfn().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) h
  integer n_data
  real ( kind = rk ) t1
  real ( kind = rk ) t2
  real ( kind = rk ) tfn

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test01():'
  write ( *, '(a)' ) '  tfn() computes the Owen T function.'
  write ( *, '(a)' ) '  Compare with tabulated values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) '          H               A           ', &
  'T                         T                       DIFF'
  write ( *, '(a,a)' ) '                                     ', &
  '(Tabulated)               (TFN)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call owen_values ( n_data, h, a, t1 )

    if ( n_data == 0 ) then
      exit
    end if

    t2 = tfn ( h, a )

    write ( *, '(2x,f14.6,2x,f14.6,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    h, a, t1, t2, abs ( t1 - t2 )

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! test02() tests tha().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) h
  integer n_data
  real ( kind = rk ) t1
  real ( kind = rk ) t2
  real ( kind = rk ) tha

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test02():'
  write ( *, '(a)' ) '  tha() computes the Owen T function.'
  write ( *, '(a)' ) '  Compare with tabulated values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) '          H               A           ', &
 'T                         T                       DIFF'
  write ( *, '(a,a)' ) '                                     ', &
 '(Tabulated)               (THA)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call owen_values ( n_data, h, a, t1 )

    if ( n_data == 0 ) then
      exit
    end if

    t2 = tha ( h, 1.0D+00, a, 1.0D+00 )

    write ( *, '(2x,f14.6,2x,f14.6,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    h, a, t1, t2, abs ( t1 - t2 )

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
!    27 August 2021
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

