subroutine cosine_transform_data ( n, d, c )

!*****************************************************************************80
!
!! cosine_transform_data() does a cosine transform on a vector of data.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 August 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data points.
!
!    Input, real ( kind = rk ) D(N), the vector of data.
!
!    Output, real ( kind = rk ) C(N), the cosine transform coefficients.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) c(n)
  real ( kind = rk ) d(n)
  integer i
  integer j
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00

  c(1:n) = 0.0D+00

  do i = 1, n
    do j = 1, n
      c(i) = c(i) + cos ( r8_pi * real ( 2 * j - 1, kind = rk ) &
        * real ( i - 1, kind = rk ) / 2.0D+00 / real ( n, kind = rk ) ) * d(j);
    end do
  end do

  c(1:n) = c(1:n) * sqrt ( 2.0D+00 / real ( n, kind = rk ) )

  return
end
subroutine cosine_transform_inverse ( n, c, d )

!*****************************************************************************80
!
!! COSINE_TRANSFORM_INVERSE does an inverse cosine transform.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 August 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data points.
!
!    Input, real ( kind = rk ) C(N), the cosine transform coefficients.
!
!    Output, real ( kind = rk ) D(N), the vector of data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) c(n)
  real ( kind = rk ) d(n)
  integer i
  integer j
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00

  do i = 1, n
    d(i) = c(1) / 2.0D+00
    do j = 2, n
      d(i) = d(i) + cos ( r8_pi * real ( 2 * i - 1, kind = rk ) &
        * real ( j - 1, kind = rk ) / 2.0D+00 / real ( n, kind = rk ) ) * c(j)
    end do
  end do

  d(1:n) = d(1:n) * sqrt ( 2.0D+00 / real ( n, kind = rk ) )

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

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
