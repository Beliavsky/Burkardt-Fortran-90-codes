program main

!*****************************************************************************80
!
!! asa136_test() tests asa136().
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
  write ( *, '(a)' ) 'asa136_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test asa136().'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'asa136_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ()

!*****************************************************************************80
!
!! test01() tests asa136().
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

  integer, parameter :: k = 5
  integer, parameter :: m = 100
  integer, parameter :: n = 2

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) an1(k)
  real ( kind = rk ) an2(k)
  real ( kind = rk ) c(k,n)
  real ( kind = rk ) d(m)
  integer i
  integer ic1(m)
  integer ic2(m)
  integer ifault
  integer iter
  integer itran(k)
  integer j
  integer live(k)
  integer nc(k)
  integer nc_sum
  integer ncp(k)
  real ( kind = rk ) wss(k)
  real ( kind = rk ) wss_sum

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test01():'
  write ( *, '(a)' ) '  asa136() implements the K-means algorithm.'
!
!  Read the data.
!
  open ( unit = 1, file = 'points_100.txt', status = 'old' )

  do i = 1, m
    read ( 1, * ) ( a(i,j), j = 1, n )
  end do

  close ( unit = 1 )
!
!  Print a few data values.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  First 5 data values:'
  write ( *, '(a)' ) ' '

  do i = 1, 5
    write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) i, ( a(i,j), j = 1, n )
  end do
!
!  Initialize the cluster centers.
!  Here, we arbitrarily make the first K data points cluster centers.
!
  do i = 1, k
    do j = 1, n
      c(i,j) = a(i,j)
    end do
  end do

  iter = 50
!
!  Compute the clusters.
!
  call kmns ( a, m, n, c, k, ic1, ic2, nc, an1, an2, ncp, d, &
    itran, live, iter, wss, ifault )

  if ( ifault /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST01 - Fatal error!'
    write ( *, '(a,i8)' ) '  KMNS returned IFAULT = ', ifault
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Cluster  Population  Energy'
  write ( *, '(a)' ) ' '

  nc_sum = 0
  wss_sum = 0.0D+00

  do i = 1, k
    write ( *, '(2x,i8,2x,i8,2x,g14.6)' ) i, nc(i), wss(i)
    nc_sum = nc_sum + nc(i)
    wss_sum = wss_sum + wss(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(2x,a8,2x,i8,2x,g14.6)' ) '   Total', nc_sum, wss_sum

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

