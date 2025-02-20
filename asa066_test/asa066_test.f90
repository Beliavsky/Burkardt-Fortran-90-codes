program main

!*****************************************************************************80
!
!! asa066_test() tests asa066().
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
  write ( *, '(a)' ) 'asa066_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test asa066().'

  call test01 ( )
  call test02 ( )
  call test03 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'asa066_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! test01() compares alnorm() against tabulated values.
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

  real ( kind = rk ) alnorm
  real ( kind = rk ) fx
  real ( kind = rk ) fx2
  integer n_data
  logical upper
  real ( kind = rk ) x

  upper = .false.

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test01():'
  write ( *, '(a)' ) '  Compare tabulated values of the normal '
  write ( *, '(a)' ) '  Cumulative Density Function against values'
  write ( *, '(a)' ) '  computed by alnorm().'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '         X        CDF                       CDF', &
    '                    DIFF'
  write ( *, '(a)' ) &
    '               (tabulated)                 (ALNORM)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do
    call normal_01_cdf_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = alnorm ( x, upper )

    write ( *, '(2x,f10.4,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    x, fx, fx2, abs ( fx - fx2 )

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! test02() compares normp() against tabulated values.
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

  real ( kind = rk ) fx
  real ( kind = rk ) fx2
  real ( kind = rk ) fx3
  integer n_data
  real ( kind = rk ) pdf
  logical upper
  real ( kind = rk ) x

  upper = .false.

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test02():'
  write ( *, '(a)' ) '  Compare tabulated values of the normal '
  write ( *, '(a)' ) '  Cumulative Density Function against values'
  write ( *, '(a)' ) '  computed by normp().'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '         X        CDF                       CDF', &
    '                    DIFF'
  write ( *, '(a)' ) &
    '               (tabulated)                 (NORMP)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call normal_01_cdf_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    call normp ( x, fx2, fx3, pdf )

    write ( *, '(2x,f10.4,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    x, fx, fx2, abs ( fx - fx2 )

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! test03() compares nprob() against tabulated values.
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

  real ( kind = rk ) fx
  real ( kind = rk ) fx2
  real ( kind = rk ) fx3
  integer n_data
  real ( kind = rk ) pdf
  logical upper
  real ( kind = rk ) x

  upper = .false.

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test03():'
  write ( *, '(a)' ) '  Compare tabulated values of the normal '
  write ( *, '(a)' ) '  Cumulative Density Function against values'
  write ( *, '(a)' ) '  computed by nprob().'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '         X        CDF                       CDF', &
    '                    DIFF'
  write ( *, '(a)' ) &
    '               (tabulated)                 (NPROB)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call normal_01_cdf_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    call nprob ( x, fx2, fx3, pdf )

    write ( *, '(2x,f10.4,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    x, fx, fx2, abs ( fx - fx2 )

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
 
