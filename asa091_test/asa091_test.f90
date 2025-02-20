program main

!*****************************************************************************80
!
!! asa091_test() tests asa091().
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
  write ( *, '(a)' ) 'asa091_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test asa091().'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'asa091_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! test01() makes a single simple calculation with ppchi2().
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

  real ( kind = rk ) g
  integer ifault
  real ( kind = rk ) p
  real ( kind = rk ) ppchi2
  real ( kind = rk ) v
  real ( kind = rk ) value
  real ( kind = rk ), parameter :: value_correct = 0.4D+00

  p = 0.017523D+00
  v = 4.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test01():'
  write ( *, '(a)' ) '  ppchi2() inverts the Chi-Squared CDF.'

  g = lgamma ( v / 2.0D+00 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g24.16)' ) '  P =                  ', p
  write ( *, '(a,g24.16)' ) '  V =                  ', v
  write ( *, '(a,g24.16)' ) '  G Log(Gamma(V/2)) =  ', g

  value = ppchi2 ( p, v, g, ifault )

  write ( *, '(a,g24.16)' ) '  VALUE =              ', value
  write ( *, '(a,g24.16)' ) '  VALUE (correct) =    ', value_correct

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'  ) '  Error flag IFAULT = ', ifault

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! test02() compares ppchi2() against tabulated values.
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

  integer a
  real ( kind = rk ) fx
  integer n_data
  real ( kind = rk ) g
  real ( kind = rk ) ppchi2
  integer ifault
  real ( kind = rk ) v
  real ( kind = rk ) x
  real ( kind = rk ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test02():'
  write ( *, '(a)' ) '  Compare tabulated values of the Chi-Squared '
  write ( *, '(a)' ) '  Cumulative Density Function against values'
  write ( *, '(a)' ) '  computed by ppchi2().'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
  '         N        CDF       X                        ', &
  ' X2                    DIFF'
  write ( *, '(a,a)' ) &
  '                           (tabulated)               ', &
  '(PPCHI2)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call chi_square_cdf_values ( n_data, a, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    v = dble ( a )

    g = lgamma ( v / 2.0D+00 )

    x2 = ppchi2 ( fx, v, g, ifault )

    write ( *, '(2x,i8,2x,f10.4,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    a, fx, x, x2, abs ( x - x2 )

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
 
