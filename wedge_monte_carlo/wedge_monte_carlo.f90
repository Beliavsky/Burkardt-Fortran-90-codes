subroutine monomial_value ( m, n, e, x, v )

!*****************************************************************************80
!
!! monomial_value() evaluates a monomial.
!
!  Discussion:
!
!    This routine evaluates a monomial of the form
!
!      product ( 1 <= i <= m ) x(i)^e(i)
!
!    where the exponents are nonnegative integers.  Note that
!    if the combination 0^0 is encountered, it should be treated
!    as 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 April 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of points at which the
!    monomial is to be evaluated.
!
!    Input, integer E(M), the exponents.
!
!    Input, real ( kind = rk ) X(M,N), the point coordinates.
!
!    Output, real ( kind = rk ) V(N), the value of the monomial.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer e(m)
  integer i
  real ( kind = rk ) v(n)
  real ( kind = rk ) x(m,n)

  v(1:n) = 1.0D+00

  do i = 1, m
    if ( 0 /= e(i) ) then
      v(1:n) = v(1:n) * x(i,1:n) ** e(i)
    end if
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
function wedge01_volume ( )

!*****************************************************************************80
!
!! WEDGE01_VOLUME returns the volume of the unit wedge in 3D.
!
!  Discussion:
!
!    The unit wedge is:
!
!      0 <= X
!      0 <= Y
!      X + Y <= 1
!      -1 <= Z <= 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 August 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) WEDGE01_VOLUME, the volume of the unit wedge.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) wedge01_volume

  wedge01_volume = 1.0D+00

  return
end
subroutine wedge01_integral ( e, value )

!*****************************************************************************80
!
!! WEDGE01_INTEGRAL returns the integral of a monomial in the unit wedge in 3D.
!
!  Discussion:
!
!    This routine returns the integral of
!
!      product ( 1 <= I <= 3 ) X(I)^E(I)
!
!    over the unit wedge.
!
!    The integration region is:
!
!      0 <= X
!      0 <= Y
!      X + Y <= 1
!      -1 <= Z <= 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 August 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer E(3), the exponents.
!
!    Output, real ( kind = rk ) VALUE, the integral of the monomial.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer e(3)
  integer i
  integer k
  real ( kind = rk ) value

  value = 1.0D+00

  k = e(1)

  do i = 1, e(2)
    k = k + 1
    value = value * real ( i, kind = rk ) / real ( k, kind = rk )
  end do

  k = k + 1
  value = value / real ( k, kind = rk )

  k = k + 1
  value = value / real ( k, kind = rk )
!
!  Now account for integration in Z.
!
  if ( e(3) == - 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WEDGE01_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  E(3) = -1 is not a legal input.'
    stop 1
  else if ( mod ( e(3), 2 ) == 1 ) then
    value = 0.0D+00
  else
    value = value * 2.0D+00 / real ( e(3) + 1, kind = rk )
  end if

  return
end
subroutine wedge01_sample ( n, x )

!*****************************************************************************80
!
!! WEDGE01_SAMPLE samples points uniformly from the unit wedge in 3D.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 August 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Reuven Rubinstein,
!    Monte Carlo Optimization, Simulation, and Sensitivity 
!    of Queueing Networks,
!    Krieger, 1992,
!    ISBN: 0894647644,
!    LC: QA298.R79.
!
!  Parameters:
!
!    Input, integer N, the number of points.
!
!    Output, real ( kind = rk ) X(3,N), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 3
  integer n

  real ( kind = rk ) e(4)
  real ( kind = rk ) e_sum
  integer j
  real ( kind = rk ) x(m,n)

  do j = 1, n

    call random_number ( harvest = e(1:4) )

    e(1:3) = - log ( e(1:3) )

    e_sum = sum ( e(1:3) )

    x(1:2,j) = e(1:2) / e_sum
    x(3,j) = 2.0D+00 * e(4) - 1.0D+00

  end do

  return
end
