program main

!*****************************************************************************80
!
!! toms659_test() tests toms659().
!
!  Modified:
!
!    14 March 2021
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: maxdim = 1111

  integer atmost
  integer dimen
  real ( kind = rk ) f
  integer i
  integer j
  real ( kind = rk ) quasi(maxdim)
  real ( kind = rk ) sum
  real ( kind = rk ) t1
  real ( kind = rk ) t2

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'toms659_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test TOMS659().'

  do

    read (*,*) dimen, atmost

    if ( dimen == 0 ) then
      exit
    end if

    write ( *, '(a)' ) ''
    write ( *, '(a,i8)' ) '  DIMENSION = ', DIMEN
    write ( *, '(a,i8)' ) '  ATMOST = ', ATMOST

    call cpu_time ( t1 )

    call insobl ( dimen, atmost )

    write ( *, * ) '  Start time = ',T1
    write ( *, '(a)' ) '  I = Iteration number'
    write ( *, '(a)' ) '  EI = Estimate of integral'
    write ( *, '(a)' ) ''

    sum = 0.0D+00
    do i = 1, atmost

      call i4_sobol ( dimen, quasi )

      f = 1.0
      do j = 1, dimen
        f = f * abs ( 4.0D+00 * quasi(j) - 2.0D+00 )
      end do
      sum = sum + f

      if ( mod ( i, 5000 ) == 0 ) then
        write ( *, '(a,i6)' ) '  i = ', i
        write ( *,*) '  ei = ', sum / i
        call cpu_time ( t2 )
        write ( *,*) '  TIMEI = ', t2 - t1
        write ( *,'(1h )')
      end if

    end do

    write ( *, * ) ' EI = ', sum / atmost

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'toms659_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  call timestamp ( )

  stop 0
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

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end

