program main

!*****************************************************************************80
!
!! closest_point_brute_test() tests closest_point_brute().
!
!  Discussion:
!
!    We are given R, a set of NR points in M dimensions.
!
!    We are given S, a points in M dimensions.
!
!    We seek the index J of the point R(J)
!    which is nearest to S over all points in R.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 August 2024
!
!  Author:
!
!    John Burkardt
!
!  Local:
!
!    integer M, the spatial dimension.
!
!    integer NR, the number of data points.
!
!    real ( kind = rk8 ) R(M,NR), the data points.
!
!    real ( kind = rk8 ) RT(NR,M), the transposed data points.
!
!    real ( kind = rk8 ) S(M), the sample point. 
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer m
  real ( kind = rk8 ) near_dist
  integer near_index
  integer nr
  real ( kind = rk8 ), allocatable :: r(:,:)
  real ( kind = rk8 ), allocatable :: s(:)
  real ( kind = rk8 ) t1
  real ( kind = rk8 ) t2

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'closest_point_brute_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test closest_point_brute().'

  m = 2
  write ( *, '(a,i2)' ) '  Spatial dimension m = ', m
  write ( *, '(a)' ) ''

  nr = 1

  do while ( nr <= 65536 )

    allocate ( r(1:nr,1:m) )
    call random_number ( harvest = r )
 
    allocate ( s(1:m) )
    call random_number ( harvest = s )

    call cpu_time ( t1 )
    call closest_point_brute ( m, nr, r, s, near_index, near_dist )
    call cpu_time ( t2 )

    write ( *, '(i5,2x,i5,2x,g8.2,2x,g8.2)' ) &
      nr, near_index, near_dist, ( t2 - t1 )

    deallocate ( r )
    deallocate ( s )

    nr = nr * 2

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'closest_point_brute_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

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
