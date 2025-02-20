subroutine moment ( n, x, y, p, q, nu_pq )

!*****************************************************************************80
!
!! MOMENT computes an unnormalized moment of a polygon.
!
!  Discussion:
!
!    Nu(P,Q) = Integral ( x, y in polygon ) x^p y^q dx dy
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carsten Steger,
!    On the calculation of arbitrary moments of polygons,
!    Technical Report FGBV-96-05,
!    Forschungsgruppe Bildverstehen, Informatik IX,
!    Technische Universitaet Muenchen, October 1996.
!
!  Parameters:
!
!    Input, integer N, the number of vertices of the polygon.
!
!    Input, real ( kind = rk ) X(N), Y(N), the vertex coordinates.
!
!    Input, integer P, Q, the indices of the moment.
!
!    Output, real ( kind = rk ) NU_PQ, the unnormalized moment Nu(P,Q).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer k
  integer l
  real ( kind = rk ) nu_pq
  integer p
  integer q
  real ( kind = rk ) r8_choose
  real ( kind = rk ) s_pq
  real ( kind = rk ) x(n)
  real ( kind = rk ) xi
  real ( kind = rk ) xj
  real ( kind = rk ) y(n)
  real ( kind = rk ) yi
  real ( kind = rk ) yj

  nu_pq = 0.0D+00

  xj = x(n)
  yj = y(n)

  do i = 1, n

    xi = x(i)
    yi = y(i)

    s_pq = 0.0D+00
    do k = 0, p
      do l = 0, q
        s_pq = s_pq &
          + r8_choose ( k + l, l ) * r8_choose ( p + q - k - l, q - l ) &
          * xi ** k * xj ** ( p - k ) &
          * yi ** l * yj ** ( q - l )
      end do
    end do

    nu_pq = nu_pq + ( xj * yi - xi * yj ) * s_pq

    xj = xi
    yj = yi

  end do

  nu_pq = nu_pq / real ( p + q + 2, kind = rk ) &
    / real ( p + q + 1, kind = rk ) &
    / r8_choose ( p + q, p )

  return
end
subroutine moment_central ( n, x, y, p, q, mu_pq )

!*****************************************************************************80
!
!! MOMENT_CENTRAL computes central moments of a polygon.
!
!  Discussion:
!
!    The central moment Mu(P,Q) is defined by
!
!      Mu(P,Q) = Integral ( polygon ) (x-Alpha(1,0))^p (y-Alpha(0,1))^q dx dy
!              / Area ( polygon )
!
!    where 
!
!      Alpha(1,0) = Integral ( polygon ) x dx dy / Area ( polygon )
!      Alpha(0,1) = Integral ( polygon ) y dx dy / Area ( polygon )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carsten Steger,
!    On the calculation of arbitrary moments of polygons,
!    Technical Report FGBV-96-05,
!    Forschungsgruppe Bildverstehen, Informatik IX,
!    Technische Universitaet Muenchen, October 1996.
!
!  Parameters:
!
!    Input, integer N, the number of vertices of the polygon.
!
!    Input, real ( kind = rk ) X(N), Y(N), the vertex coordinates.
!
!    Input, integer P, Q, the indices of the moment.
!
!    Output, real ( kind = rk ) MU_PQ, the unnormalized moment Mu(P,Q).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) alpha_01
  real ( kind = rk ) alpha_10
  real ( kind = rk ) alpha_ij
  integer i
  integer j
  real ( kind = rk ) mu_pq
  integer p
  integer q
  real ( kind = rk ) r8_choose
  real ( kind = rk ) r8_mop
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)

  call moment_normalized ( n, x, y, 1, 0, alpha_10 )
  call moment_normalized ( n, x, y, 0, 1, alpha_01 )

  mu_pq = 0.0D+00

  do i = 0, p
    do j = 0, q

      call moment_normalized ( n, x, y, i, j, alpha_ij )

      mu_pq = mu_pq + r8_mop ( p + q - i - j ) &
        * r8_choose ( p, i ) * r8_choose ( q, j ) &
        * alpha_10 ** ( p - i ) * alpha_01 ** ( q - j ) * alpha_ij

    end do
  end do

  return
end
subroutine moment_normalized ( n, x, y, p, q, alpha_pq )

!*****************************************************************************80
!
!! MOMENT_NORMALIZED computes a normalized moment of a polygon.
!
!  Discussion:
!
!    Alpha(P,Q) = Integral ( x, y in polygon ) x^p y^q dx dy / Area ( polygon )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carsten Steger,
!    On the calculation of arbitrary moments of polygons,
!    Technical Report FGBV-96-05,
!    Forschungsgruppe Bildverstehen, Informatik IX,
!    Technische Universitaet Muenchen, October 1996.
!
!  Parameters:
!
!    Input, integer N, the number of vertices of the polygon.
!
!    Input, real ( kind = rk ) X(N), Y(N), the vertex coordinates.
!
!    Input, integer P, Q, the indices of the moment.
!
!    Output, real ( kind = rk ) ALPHA_PQ, the normalized moment Alpha(P,Q).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) alpha_pq
  real ( kind = rk ) nu_00
  real ( kind = rk ) nu_pq
  integer p
  integer q
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)

  call moment ( n, x, y, p, q, nu_pq )
  call moment ( n, x, y, 0, 0, nu_00 )

  alpha_pq = nu_pq / nu_00

  return
end
function r8_choose ( n, k )

!*****************************************************************************80
!
!! R8_CHOOSE computes the binomial coefficient C(N,K) as an R8.
!
!  Discussion:
!
!    The value is calculated in such a way as to avoid overflow and
!    roundoff.  The calculation is done in R8 arithmetic.
!
!    The formula used is:
!
!      C(N,K) = N! / ( K! * (N-K)! )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    ML Wolfson, HV Wright,
!    Algorithm 160:
!    Combinatorial of M Things Taken N at a Time,
!    Communications of the ACM,
!    Volume 6, Number 4, April 1963, page 161.
!
!  Parameters:
!
!    Input, integer N, K, are the values of N and K.
!
!    Output, real ( kind = rk ) R8_CHOOSE, the number of combinations of N
!    things taken K at a time.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer k
  integer mn
  integer mx
  integer n
  real ( kind = rk ) r8_choose
  real ( kind = rk ) value

  mn = min ( k, n - k )

  if ( mn < 0 ) then

    value = 0.0D+00

  else if ( mn == 0 ) then

    value = 1.0D+00

  else

    mx = max ( k, n - k )
    value = real ( mx + 1, kind = rk )

    do i = 2, mn
      value = ( value * real ( mx + i, kind = rk ) ) / real ( i, kind = rk )
    end do

  end if

  r8_choose = value

  return
end
function r8_mop ( i )

!*****************************************************************************80
!
!! R8_MOP returns the I-th power of -1 as an R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = rk ) value.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, the power of -1.
!
!    Output, real ( kind = rk ) R8_MOP, the I-th power of -1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  real ( kind = rk ) r8_mop

  if ( mod ( i, 2 ) == 0 ) then
    r8_mop = + 1.0D+00
  else
    r8_mop = - 1.0D+00
  end if

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

  integer, parameter :: rk = kind ( 1.0D+00 )

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
