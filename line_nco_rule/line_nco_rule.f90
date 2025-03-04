subroutine line_nco_rule ( n, a, b, x, w )

!*****************************************************************************80
!
!! line_nco_rule() computes a Newton-Cotes Open (NCO) quadrature rule.
!
!  Discussion:
!
!    The integral:
!
!      Integral ( A <= X <= B ) F(X) dx
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= N ) W(I) * F ( X(I) ).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order.
!
!    Input, real ( kind = rk ) A, B, the endpoints of the interval.
!
!    Input, real ( kind = rk ) X(N), the abscissas.
!
!    Output, real ( kind = rk ) W(N), the weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) d(n)
  integer i
  integer j
  integer k
  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n)
  real ( kind = rk ) y_a
  real ( kind = rk ) y_b
!
!  Define the points X.
!
  call r8vec_linspace2 ( n, a, b, x )
!
!  Compute the Lagrange basis polynomial which is 1 at X(I),
!  and zero at the other nodes.
!
  do i = 1, n

    d(1:n) = 0.0D+00
    d(i) = 1.0D+00

    do j = 2, n
      do k = j, n
        d(n+j-k) = ( d(n+j-k-1) - d(n+j-k) ) / ( x(n+1-k) - x(n+j-k) )
      end do
    end do

    do j = 1, n - 1
      do k = 1, n - j
        d(n-k) = d(n-k) - x(n-k-j+1) * d(n-k+1)
      end do
    end do
!
!  Evaluate the antiderivative of the polynomial at the endpoints.
!
    y_a = d(n) / real ( n, kind = rk )
    do j = n - 1, 1, -1
      y_a = y_a * a + d(j) / real ( j, kind = rk )
    end do
    y_a = y_a * a

    y_b = d(n) / real ( n, kind = rk )
    do j = n - 1, 1, -1
      y_b = y_b * b + d(j) / real ( j, kind = rk )
    end do
    y_b = y_b * b

    w(i) = y_b - y_a

  end do

  return
end
subroutine r8vec_linspace2 ( n, a, b, x )

!*****************************************************************************80
!
!! R8VEC_LINSPACE2 creates a vector of linearly spaced values.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    4 points evenly spaced between 0 and 12 will yield 2, 4, 6, 8, 10.
!
!    In other words, the interval is divided into N+1 even subintervals,
!    and the endpoints of internal intervals are used as the points.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real ( kind = rk ) A, B, the first and last entries.
!
!    Output, real ( kind = rk ) X(N), a vector of linearly spaced data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a
  real ( kind = rk ) b
  integer i
  real ( kind = rk ) x(n)

  do i = 1, n
    x(i) = ( real ( n  - i + 1, kind = rk ) * a &
           + real (      i,     kind = rk ) * b ) &
           / real ( n      + 1, kind = rk )
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
