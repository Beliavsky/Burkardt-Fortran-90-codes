subroutine i4mat_shortest_path ( n, m )

!*****************************************************************************80
!
!! i4mat_shortest_path() computes the shortest distance between pairs of points.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 March 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Floyd,
!    Algorithm 97, Shortest Path,
!    Communications of the ACM,
!    Volume 5, Number 6, June 1962, page 345.
!
!  Parameters:
!
!    Input, integer N, the number of points.
!
!    Input/output, integer M(N,N).
!    On input, M(I,J) contains the length of the direct link between 
!    nodes I and J, or HUGE if there is no direct link.
!    On output, M(I,J) contains the distance between nodes I and J,
!    that is, the length of the shortest path between them.  If there
!    is no such path, then M(I,J) will remain HUGE.
!
  implicit none

  integer n

  integer i
  integer i4_inf
  integer j
  integer k
  integer m(n,n)
  integer s

  i4_inf = huge ( i4_inf )

  do i = 1, n
    do j = 1, n
      if ( m(j,i) < i4_inf ) then
        do k = 1, n
          if ( m(i,k) < i4_inf ) then
            s = m(j,i) + m(i,k)
            if ( s < m(j,k) ) then
              m(j,k) = s
            end if
          end if
        end do
      end if
    end do
  end do

  return
end
subroutine r8mat_shortest_path ( n, m )

!*****************************************************************************80
!
!! R8MAT_SHORTEST_PATH computes the shortest distance between pairs of points.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 March 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Floyd,
!    Algorithm 97, Shortest Path,
!    Communications of the ACM,
!    Volume 5, Number 6, June 1962, page 345.
!
!  Parameters:
!
!    Input, integer N, the number of points.
!
!    Input/output, real ( kind = rk ) M(N,N).
!    On input, M(I,J) contains the length of the direct link between 
!    nodes I and J, or HUGE if there is no direct link.
!    On output, M(I,J) contains the distance between nodes I and J,
!    that is, the length of the shortest path between them.  If there
!    is no such path, then M(I,J) will remain HUGE.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer j
  integer k
  real ( kind = rk ) m(n,n)
  real ( kind = rk ) r8_inf
  real ( kind = rk ) s

  r8_inf = huge ( r8_inf )

  do i = 1, n
    do j = 1, n
      if ( m(j,i) < r8_inf ) then
        do k = 1, n
          if ( m(i,k) < r8_inf ) then
            s = m(j,i) + m(i,k)
            if ( s < m(j,k) ) then
              m(j,k) = s
            end if
          end if
        end do
      end if
    end do
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
!  Parameters:
!
!    None
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
