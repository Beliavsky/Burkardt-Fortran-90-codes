program main

!*****************************************************************************80
!
!! asa058_test() tests asa058().
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
  write ( *, '(a)' ) 'asa058_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test asa058().'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'asa058_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! test01() tests clustr().
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

  integer, parameter :: k = 5
  integer, parameter :: m = 2
  integer, parameter :: n = 100

  integer b(n)
  real ( kind = rk ) d(k,m)
  real ( kind = rk ) dev(k)
  real ( kind = rk ) dev_sum
  integer e(k)
  integer e_sum
  real ( kind = rk ) f(n)
  integer i
  integer j
  integer k2
  integer nz
  real ( kind = rk ) x(n,m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test01():'
  write ( *, '(a)' ) '  clustr() uses K-means for clustering data.'
  write ( *, '(a)' ) '  Applied Statistics Algorithm 58'
!
!  Read the data.
!
  write (  *, '(a)' ) ' '
  write (  *, '(a)' ) '  Reading the data.'

  open ( unit = 1, file = 'points_100.txt', status = 'old' )

  do i = 1, n
    read ( 1, * ) ( x(i,j), j = 1, m )
  end do

  close ( unit = 1 )
!
!  Print a few data values.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  First 5 data values:'
  write ( *, '(a)' ) ' '

  do i = 1, 5
    write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) i, ( x(i,j), j = 1, m )
  end do
!
!  Initialize the cluster centers arbitrarily.
!
  d(1:k,1:m) = x(1:k,1:m)
!
!  Compute the clusters.
!
  nz = 1
  k2 = k

  call clustr ( x, d, dev, b, f, e, n, m, k, nz, k2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Cluster  Population  Energy'
  write ( *, '(a)' ) ' '

  do i = 1, k
    write ( *, '(2x,i8,2x,i8,2x,g14.6)' ) i, e(i), dev(i)
  end do

  e_sum = sum ( e(1:k) )
  dev_sum = sum ( dev(1:k) )

  write ( *, '(a)' ) ' '
  write ( *, '(2x,a8,2x,i8,2x,g14.6)' ) '   Total', e_sum, dev_sum

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
 
