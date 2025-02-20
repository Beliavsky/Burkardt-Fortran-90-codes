subroutine pwl_product_integral ( a, b, f_num, f_x, f_v, g_num, &
  g_x, g_v, integral )

!*****************************************************************************80
!
!! pwl_PRODUCT_INTEGRAL: piecewise linear product integral.
!
!  Discussion:
!
!    We are given two piecewise linear functions F(X) and G(X) and we wish
!    to compute the exact value of the integral
!
!      INTEGRAL = Integral ( A <= X <= B ) F(X) * G(X) dx
!
!    The functions F(X) and G(X) are defined as tables of coordinates X and
!    values V.  A piecewise linear function is evaluated at a point X by 
!    evaluating the interpolant to the data at the endpoints of the interval 
!    containing X.  
!
!    It must be the case that A <= B.
!
!    It must be the case that the node coordinates F_X(*) and G_X(*) are
!    given in ascending order.
!
!    It must be the case that:
!
!      F_X(1) <= A and B <= F_X(F_NUM)
!      G_X(1) <= A and B <= G_X(G_NUM)
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, the limits of integration.
!
!    Input, integer F_NUM, the number of nodes for F.
!
!    Input, real ( kind = rk ) F_X(F_NUM), the node coordinates for F.
!
!    Input, real ( kind = rk ) F_V(F_NUM), the nodal values for F.
!
!    Input, integer G_NUM, the number of nodes for G.
!
!    Input, real ( kind = rk ) G_X(G_NUM), the node coordinates for G.
!
!    Input, real ( kind = rk ) G_V(G_NUM), the nodal values for G.
!
!    Output, real ( kind = rk ) INTEGRAL, the integral of F(X) * G(X)
!    from A to B.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer f_num
  integer g_num

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) bit
  integer f_left
  real ( kind = rk ) f_v(f_num)
  real ( kind = rk ) f_x(f_num)
  real ( kind = rk ) f0
  real ( kind = rk ) f1
  real ( kind = rk ) fl
  real ( kind = rk ) fr
  integer g_left
  real ( kind = rk ) g_v(g_num)
  real ( kind = rk ) g_x(g_num)
  real ( kind = rk ) g0
  real ( kind = rk ) g1
  real ( kind = rk ) gl
  real ( kind = rk ) gr
  real ( kind = rk ) h0
  real ( kind = rk ) h1
  real ( kind = rk ) h2
  integer i
  real ( kind = rk ) integral
  real ( kind = rk ) xl
  real ( kind = rk ) xr
  real ( kind = rk ) xr_max

  integral = 0.0D+00

  if ( f_x(f_num) <= a .or. g_x(g_num) <= a ) then
    return
  end if

  if ( f_num < 2 .or. g_num < 2 ) then
    return
  end if

  xr = a

  f_left = 1
  call r8vec_bracket3 ( f_num, f_x, xr, f_left )
  fr = f_v(f_left) + ( xr - f_x(f_left) ) * ( f_v(f_left+1) - f_v(f_left) ) &
    / ( f_x(f_left+1) - f_x(f_left) )

  g_left = 1
  call r8vec_bracket3 ( g_num, g_x, xr, g_left )
  gr = g_v(g_left) + ( xr - g_x(g_left) ) * ( g_v(g_left+1) - g_v(g_left) ) &
    / ( g_x(g_left+1) - g_x(g_left) )

  xr_max = b
  xr_max = min ( xr_max, f_x(f_num) )
  xr_max = min ( xr_max, g_x(g_num) )

  do while ( xr < xr_max )
!
!  Shift right values to left.
!
    xl = xr
    fl = fr
    gl = gr
!
!  Determine the new right values.
!  The hard part is figuring out how to advance XR some, but not too much.
!
    xr = xr_max

    do i = 1, 2
      if ( f_left + i <= f_num ) then
        if ( xl < f_x(f_left+i) .and. f_x(f_left+i) < xr ) then
          xr = f_x(f_left+i)
          exit
        end if
      end if
    end do

    do i = 1, 2
      if ( g_left + i <= g_num ) then
        if ( xl < g_x(g_left+i) .and. g_x(g_left+i) < xr ) then
          xr = g_x(g_left+i)
          exit
        end if
      end if
    end do

    call r8vec_bracket3 ( f_num, f_x, xr, f_left )
    fr = f_v(f_left) + ( xr - f_x(f_left) ) * ( f_v(f_left+1) - f_v(f_left) ) &
      / ( f_x(f_left+1) - f_x(f_left) )

    call r8vec_bracket3 ( g_num, g_x, xr, g_left )
    gr = g_v(g_left) + ( xr - g_x(g_left) ) * ( g_v(g_left+1) - g_v(g_left) ) &
      / ( g_x(g_left+1) - g_x(g_left) )
!
!  Form the linear polynomials for F(X) and G(X) over [XL,XR],
!  then the product H(X), integrate H(X) and add to the running total.
!
    if ( epsilon ( xl - xr ) <= abs ( xr - xl ) ) then

      f1 = fl - fr
      f0 = fr * xl - fl * xr

      g1 = gl - gr
      g0 = gr * xl - gl * xr

      h2 = f1 * g1
      h1 = f1 * g0 + f0 * g1
      h0 = f0 * g0

      h2 = h2 / 3.0D+00
      h1 = h1 / 2.0D+00

      bit = ( ( h2 * xr + h1 ) * xr + h0 ) * xr &
          - ( ( h2 * xl + h1 ) * xl + h0 ) * xl

      integral = integral + bit / ( xr - xl ) / ( xr - xl )

    end if

  end do

  return
end
subroutine pwl_product_quad ( a, b, f_num, f_x, f_v, g_num, &
  g_x, g_v, quad_num, quad )

!*****************************************************************************80
!
!! pwl_PRODUCT_QUAD: estimate piecewise linear product integral.
!
!  Discussion:
!
!    We are given two piecewise linear functions F(X) and G(X) and we wish
!    to estimate the value of the integral
!
!      INTEGRAL = Integral ( A <= X <= B ) F(X) * G(X) dx
!
!    The functions F(X) and G(X) are defined as tables of coordinates X and
!    values V.  A piecewise linear function is evaluated at a point X by 
!    evaluating the interpolant to the data at the endpoints of the interval 
!    containing X.  
!
!    It must be the case that A <= B.
!
!    It must be the case that the node coordinates F_X(*) and G_X(*) are
!    given in ascending order.
!
!    It must be the case that:
!
!      F_X(1) <= A and B <= F_X(F_NUM)
!      G_X(1) <= A and B <= G_X(G_NUM)
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, the limits of integration.
!
!    Input, integer F_NUM, the number of nodes for F.
!
!    Input, real ( kind = rk ) F_X(F_NUM), the node coordinates for F.
!
!    Input, real ( kind = rk ) F_V(F_NUM), the nodal values for F.
!
!    Input, integer G_NUM, the number of nodes for G.
!
!    Input, real ( kind = rk ) G_X(G_NUM), the node coordinates for G.
!
!    Input, real ( kind = rk ) G_V(G_NUM), the nodal values for G.
!
!    Input, integer QUAD_NUM, the number of quadrature points.
!
!    Output, real ( kind = rk ) QUAD, an estimate for the integral of F(X) * G(X)
!    from A to B.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer f_num
  integer g_num

  real ( kind = rk ) a
  real ( kind = rk ) a2
  real ( kind = rk ) b
  real ( kind = rk ) b2
  integer f_left
  real ( kind = rk ) f_v(f_num)
  real ( kind = rk ) f_x(f_num)
  real ( kind = rk ) fq
  integer g_left
  real ( kind = rk ) g_v(g_num)
  real ( kind = rk ) g_x(g_num)
  real ( kind = rk ) gq
  integer i
  real ( kind = rk ) quad
  integer quad_num
  real ( kind = rk ) xq

  quad = 0.0D+00

  f_left = 1
  g_left = 1

  a2 = a
  a2 = max ( a2, f_x(1) )
  a2 = max ( a2, g_x(1) )

  b2 = b
  b2 = min ( b2, f_x(f_num) )
  b2 = min ( b2, g_x(g_num) )

  do i = 1, quad_num

    xq =  ( real (                2 * i - 1, kind = rk ) * b2 & 
          + real ( 2 * quad_num - 2 * i + 1, kind = rk ) * a2 )  & 
          / real ( 2 * quad_num,             kind = rk )

    call r8vec_bracket3 ( f_num, f_x, xq, f_left )

    fq = f_v(f_left) + ( xq - f_x(f_left) ) * ( f_v(f_left+1) - f_v(f_left) ) &
      / ( f_x(f_left+1) - f_x(f_left) )

    call r8vec_bracket3 ( g_num, g_x, xq, g_left )

    gq = g_v(g_left) + ( xq - g_x(g_left) ) * ( g_v(g_left+1) - g_v(g_left) ) &
      / ( g_x(g_left+1) - g_x(g_left) )

    quad = quad + fq * gq

  end do

  quad = quad * ( b - a ) / real ( quad_num, kind = rk )

  return
end
subroutine r8vec_bracket3 ( n, t, tval, left )

!*****************************************************************************80
!
!! R8VEC_BRACKET3 finds the interval containing or nearest a given value.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The routine always returns the index LEFT of the sorted array
!    T with the property that either
!    *  T is contained in the interval [ T(LEFT), T(LEFT+1) ], or
!    *  T < T(LEFT) = T(1), or
!    *  T > T(LEFT+1) = T(N).
!
!    The routine is useful for interpolation problems, where
!    the abscissa must be located within an interval of data
!    abscissas for interpolation, or the "nearest" interval
!    to the (extreme) abscissa must be found so that extrapolation
!    can be carried out.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, length of the input array.
!
!    Input, real ( kind = rk ) T(N), an array that has been sorted
!    into ascending order.
!
!    Input, real ( kind = rk ) TVAL, a value to be bracketed by entries of T.
!
!    Input/output, integer LEFT.
!    On input, if 1 <= LEFT <= N-1, LEFT is taken as a suggestion for the
!    interval [ T(LEFT), T(LEFT+1) ] in which TVAL lies.  This interval
!    is searched first, followed by the appropriate interval to the left
!    or right.  After that, a binary search is used.
!    On output, LEFT is set so that the interval [ T(LEFT), T(LEFT+1) ]
!    is the closest to TVAL; it either contains TVAL, or else TVAL
!    lies outside the interval [ T(1), T(N) ].
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer high
  integer left
  integer low
  integer mid
  real ( kind = rk ) t(n)
  real ( kind = rk ) tval
!
!  Check the input data.
!
  if ( n < 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_BRACKET3 - Fatal error!'
    write ( *, '(a)' ) '  N must be at least 2.'
    stop
  end if
!
!  If LEFT is not between 1 and N-1, set it to the middle value.
!
  if ( left < 1 .or. n - 1 < left ) then
    left = ( n + 1 ) / 2
  end if
!
!  CASE 1: TVAL < T(LEFT):
!  Search for TVAL in [T(I), T(I+1)] for intervals I = 1 to LEFT-1.
!
  if ( tval < t(left) ) then

    if ( left == 1 ) then
      return
    else if ( left == 2 ) then
      left = 1
      return
    else if ( t(left-1) <= tval ) then
      left = left - 1
      return
    else if ( tval <= t(2) ) then
      left = 1
      return
    end if
!
!  ...Binary search for TVAL in [T(I), T(I+1)] for intervals I = 2 to LEFT-2.
!
    low = 2
    high = left - 2

    do

      if ( low == high ) then
        left = low
        return
      end if

      mid = ( low + high + 1 ) / 2

      if ( t(mid) <= tval ) then
        low = mid
      else
        high = mid - 1
      end if

    end do
!
!  CASE2: T(LEFT+1) < TVAL:
!  Search for TVAL in [T(I),T(I+1)] for intervals I = LEFT+1 to N-1.
!
  else if ( t(left+1) < tval ) then

    if ( left == n - 1 ) then
      return
    else if ( left == n - 2 ) then
      left = left + 1
      return
    else if ( tval <= t(left+2) ) then
      left = left + 1
      return
    else if ( t(n-1) <= tval ) then
      left = n - 1
      return
    end if
!
!  ...Binary search for TVAL in [T(I), T(I+1)] for intervals I = LEFT+2 to N-2.
!
    low = left + 2
    high = n - 2

    do

      if ( low == high ) then
        left = low
        return
      end if

      mid = ( low + high + 1 ) / 2

      if ( t(mid) <= tval ) then
        low = mid
      else
        high = mid - 1
      end if

    end do
!
!  CASE3: T(LEFT) <= TVAL <= T(LEFT+1):
!  T is in [T(LEFT), T(LEFT+1)], as the user said it might be.
!
  else

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
