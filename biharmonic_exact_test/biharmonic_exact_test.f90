program main

!*****************************************************************************80
!
!! biharmonic_exact_test() tests biharmonic_exact().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 August 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  real ( kind = rk8 ) c
  real ( kind = rk8 ) d
  real ( kind = rk8 ) e
  real ( kind = rk8 ) f
  real ( kind = rk8 ) g
  integer i
  real ( kind = rk8 ) R
  real ( kind = rk8 ) X(n)
  real ( kind = rk8 ) Y(n)

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'biharmonic_exact_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test biharmonic_exact().'

  call random_number ( harvest = X )
  call random_number ( harvest = Y )
  a = 1.0
  b = 2.0
  c = 3.0
  d = 4.0
  e = 5.0
  f = 6.0
  g = 7.0
!
!  Exact function 1.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Exact solution 1 residual at random points:'

  do i = 1, n
    call biharmonic_exact_r1 ( X(i), Y(i), a, b, c, d, e, f, g, R )
    write ( *, '(2x,i2,2x,f10.4,2x,f10.4,2x,g14.6)' ) i, X(i), Y(i), R
  end do
!
!  Exact function 2.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Exact solution 2 residual at random points:'

  do i = 1, n
    call biharmonic_exact_r2 ( X(i), Y(i), a, b, c, d, e, f, g, R )
    write ( *, '(2x,i2,2x,f10.4,2x,f10.4,2x,g14.6)' ) i, X(i), Y(i), R
  end do
!
!  Exact function 3.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Exact solution 3 residual at random points:'

  do i = 1, n
    call biharmonic_exact_r3 ( X(i), Y(i), a, b, c, d, e, f, R )
    write ( *, '(2x,i2,2x,f10.4,2x,f10.4,2x,g14.6)' ) i, X(i), Y(i), R
  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'biharmonic_exact_test()'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop
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
!    15 August 2021
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

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end

