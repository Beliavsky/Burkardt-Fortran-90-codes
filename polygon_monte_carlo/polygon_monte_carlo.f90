function angle_degree ( x1, y1, x2, y2, x3, y3 )

!*****************************************************************************80
!
!! angle_degree() returns the degree angle defined by three points.
!
!  Discussion:
!
!        P1
!        /
!       /
!      /
!     /
!    P2--------->P3
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 August 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, X3, Y3, the coordinates of the points
!    P1, P2, P3.
!
!    Output, real VALUE, the angle swept out by the rays, measured
!    in degrees.  0 <= VALUE < 360.  If either ray has zero length,
!    then VALUE is set to 0.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) angle_degree
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) value
  real ( kind = rk ) x
  real ( kind = rk ) x1
  real ( kind = rk ) x2
  real ( kind = rk ) x3
  real ( kind = rk ) y
  real ( kind = rk ) y1
  real ( kind = rk ) y2
  real ( kind = rk ) y3

  x = ( x3 - x2 ) * ( x1 - x2 ) + ( y3 - y2 ) * ( y1 - y2 )

  y = ( x3 - x2 ) * ( y1 - y2 ) - ( y3 - y2 ) * ( x1 - x2 )

  if ( x == 0.0D+00 .and. y == 0.0D+00 ) then

    value = 0.0D+00

  else

    value = atan2 ( y, x )

    if ( value < 0.0D+00 ) then
      value = value + 2.0D+00 * r8_pi
    end if

    value = 180.0D+00 * value / r8_pi

  end if

  angle_degree = value

  return
end
function between ( xa, ya, xb, yb, xc, yc )

!*****************************************************************************80
!
!! between() is TRUE if vertex C is between vertices A and B.
!
!  Discussion:
!
!    The points must be (numerically) collinear.
!
!    Given that condition, we take the greater of XA - XB and YA - YB
!    as a "scale" and check where C's value lies.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 May 2014
!
!  Author:
!
!    Original C version by Joseph ORourke.
!    This version by John Burkardt.
!
!  Reference:
!
!    Joseph ORourke,
!    Computational Geometry in C,
!    Cambridge, 1998,
!    ISBN: 0521649765,
!    LC: QA448.D38.
!
!  Parameters:
!
!    Input, real ( kind = rk ) XA, YA, XB, YB, XC, YC, the coordinates of 
!    the vertices.
!
!    Output, logical BETWEEN, is TRUE if C is between A and B.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  logical between
  logical collinear
  logical value
  real ( kind = rk ) xa
  real ( kind = rk ) xb
  real ( kind = rk ) xc
  real ( kind = rk ) xmax
  real ( kind = rk ) xmin
  real ( kind = rk ) ya
  real ( kind = rk ) yb
  real ( kind = rk ) yc
  real ( kind = rk ) ymax
  real ( kind = rk ) ymin

  if ( .not. collinear ( xa, ya, xb, yb, xc, yc ) ) then
    value = .false.
  else if ( abs ( ya - yb ) < abs ( xa - xb ) ) then
    xmax = max ( xa, xb )
    xmin = min ( xa, xb )
    value = ( xmin <= xc .and. xc <= xmax )
  else
    ymax = max ( ya, yb )
    ymin = min ( ya, yb )
    value = ( ymin <= yc .and. yc <= ymax )
  end if

  between = value

  return
end
function collinear ( xa, ya, xb, yb, xc, yc )

!*****************************************************************************80
!
!! collinear() returns a measure of collinearity for three points.
!
!  Discussion:
!
!    In order to deal with collinear points whose coordinates are not
!    numerically exact, we compare the area of the largest square
!    that can be created by the line segment between two of the points
!    to (twice) the area of the triangle formed by the points.
!
!    If the points are collinear, their triangle has zero area.
!    If the points are close to collinear, then the area of this triangle
!    will be small relative to the square of the longest segment.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 September 2016
!
!  Author:
!
!    Original C version by Joseph ORourke.
!    This version by John Burkardt.
!
!  Reference:
!
!    Joseph ORourke,
!    Computational Geometry in C,
!    Cambridge, 1998,
!    ISBN: 0521649765,
!    LC: QA448.D38.
!
!  Parameters:
!
!    Input, real ( kind = rk ) XA, YA, XB, YB, XC, YC, the coordinates of 
!    the vertices.
!
!    Output, logical COLLINEAR, is TRUE if the points are judged 
!    to be collinear.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) area
  logical collinear
  real ( kind = rk ), parameter :: r8_eps = 2.220446049250313D-016
  real ( kind = rk ) side_ab_sq
  real ( kind = rk ) side_bc_sq
  real ( kind = rk ) side_ca_sq
  real ( kind = rk ) side_max_sq
  real ( kind = rk ) triangle_area
  logical value
  real ( kind = rk ) xa
  real ( kind = rk ) xb
  real ( kind = rk ) xc
  real ( kind = rk ) ya
  real ( kind = rk ) yb
  real ( kind = rk ) yc

  area = triangle_area ( xa, ya, xb, yb, xc, yc )

  side_ab_sq = ( xa - xb ) ** 2 + ( ya - yb ) ** 2
  side_bc_sq = ( xb - xc ) ** 2 + ( yb - yc ) ** 2
  side_ca_sq = ( xc - xa ) ** 2 + ( yc - ya ) ** 2

  side_max_sq = max ( side_ab_sq, max ( side_bc_sq, side_ca_sq ) )

  if ( side_max_sq <= r8_eps ) then
    value = .true.
  else if ( 2.0D+00 * abs ( area ) <= r8_eps * side_max_sq ) then
    value = .true.
  else
    value = .false.
  end if

  collinear = value

  return
end
function diagonal ( im1, ip1, n, prev_node, next_node, x, y )

!*****************************************************************************80
!
!! diagonal(): VERTEX(IM1) to VERTEX(IP1) is a proper internal diagonal.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 May 2014
!
!  Author:
!
!    Original C version by Joseph ORourke.
!    This version by John Burkardt.
!
!  Reference:
!
!    Joseph ORourke,
!    Computational Geometry in C,
!    Cambridge, 1998,
!    ISBN: 0521649765,
!    LC: QA448.D38.
!
!  Parameters:
!
!    Input, integer IM1, IP1, the indices of two vertices.
!
!    Input, integer N, the number of vertices.
!
!    Input, integer PREV_NODE(N), the previous neighbor of 
!    each vertex.
!
!    Input, integer NEXT_NODE(N), the next neighbor of each vertex.
!
!    Input, real ( kind = rk ) X(N), Y(N), the coordinates of each vertex.
!
!    Output, logical DIAGONAL, the value of the test.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  logical diagonal
  logical diagonalie
  integer im1
  logical in_cone
  integer ip1
  integer next_node(n)
  integer prev_node(n)
  logical value1
  logical value2
  logical value3
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)

  value1 = in_cone ( im1, ip1, n, prev_node, next_node, x, y )
  value2 = in_cone ( ip1, im1, n, prev_node, next_node, x, y )
  value3 = diagonalie ( im1, ip1, n, next_node, x, y )

  diagonal = ( value1 .and. value2 .and. value3 )

  return
end
function diagonalie ( im1, ip1, n, next_node, x, y )

!*****************************************************************************80
!
!! diagonalie() is true if VERTEX(IM1):VERTEX(IP1) is a proper diagonal.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 May 2014
!
!  Author:
!
!    Original C version by Joseph ORourke.
!    This version by John Burkardt.
!
!  Reference:
!
!    Joseph ORourke,
!    Computational Geometry in C,
!    Cambridge, 1998,
!    ISBN: 0521649765,
!    LC: QA448.D38.
!
!  Parameters:
!
!    Input, integer IM1, IP1, the indices of two vertices.
!
!    Input, integer N, the number of vertices.
!
!    Input, integer NEXT_NODE(N), the next neighbor of each vertex.
!
!    Input, real ( kind = rk ) X(N), Y(N), the coordinates of each vertex.
!
!    Output, logical DIAGONALIE, the value of the test.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  logical diagonalie
  integer first
  integer im1
  logical intersect
  integer ip1
  integer j
  integer jp1
  integer next_node(n)
  logical value
  logical value2
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)

  first = im1
  j = first
  jp1 = next_node(first)

  value = .true.
!
!  For each edge VERTEX(J):VERTEX(JP1) of the polygon:
!
  do
!
!  Skip any edge that includes vertex IM1 or IP1.
!
    if ( j == im1 .or. j == ip1 .or. jp1 == im1 .or. jp1 == ip1 ) then

    else

      value2 = intersect ( x(im1), y(im1), x(ip1), y(ip1), x(j), y(j), &
        x(jp1), y(jp1) )

      if ( value2 ) then
        value = .false.
        exit
      end if

    end if

    j = jp1
    jp1 = next_node(j)

    if ( j == first ) then
      exit
    end if

  end do

  diagonalie = value

  return
end
function in_cone ( im1, ip1, n, prev_node, next_node, x, y )

!*****************************************************************************80
!
!! in_cone() is TRUE if the diagonal VERTEX(IM1):VERTEX(IP1) is strictly internal.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 May 2014
!
!  Author:
!
!    Original C version by Joseph ORourke.
!    This version by John Burkardt.
!
!  Reference:
!
!    Joseph ORourke,
!    Computational Geometry in C,
!    Cambridge, 1998,
!    ISBN: 0521649765,
!    LC: QA448.D38.
!
!  Parameters:
!
!    Input, integer IM1, IP1, the indices of two vertices.
!
!    Input, integer N, the number of vertices.
!
!    Input, integer PREV_NODE(N), the previous neighbor of 
!    each vertex.
!
!    Input, integer NEXT_NODE(N), the next neighbor of each vertex.
!
!    Input, real ( kind = rk ) X(N), Y(N), the coordinates of each vertex.
!
!    Output, logical IN_CONE, the value of the test.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer im1
  integer im2
  logical in_cone
  integer ip1
  integer next_node(n)
  integer prev_node(n)
  real ( kind = rk ) t1
  real ( kind = rk ) t2
  real ( kind = rk ) t3
  real ( kind = rk ) t4
  real ( kind = rk ) t5
  real ( kind = rk ) triangle_area
  logical value
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)

  im2 = prev_node(im1)
  i = next_node(im1)

  t1 = triangle_area ( x(im1), y(im1), x(i), y(i), x(im2), y(im2) )

  if ( 0.0D+00 <= t1 ) then

    t2 = triangle_area ( x(im1), y(im1), x(ip1), y(ip1), x(im2), y(im2) )
    t3 = triangle_area ( x(ip1), y(ip1), x(im1), y(im1), x(i), y(i) )
    value = ( ( 0.0D+00 < t2 ) .and. ( 0.0D+00 < t3 ) )

  else

    t4 = triangle_area ( x(im1), y(im1), x(ip1), y(ip1), x(i), y(i) )
    t5 = triangle_area ( x(ip1), y(ip1), x(im1), y(im1), x(im2), y(im2) )
    value = .not. ( ( 0.0D+00 <= t4 ) .and. ( 0.0D+00 <= t5 ) )

  end if

  in_cone = value

  return
end
function intersect ( xa, ya, xb, yb, xc, yc, xd, yd )

!*****************************************************************************80
!
!! intersect() is true if lines VA:VB and VC:VD intersect.
!
!  Discussion:
!
!    Thanks to Gene Dial for correcting the call to intersect_prop(),
!    08 September 2016.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 September 2016
!
!  Author:
!
!    Original C version by Joseph ORourke.
!    This version by John Burkardt.
!
!  Reference:
!
!    Joseph ORourke,
!    Computational Geometry in C,
!    Cambridge, 1998,
!    ISBN: 0521649765,
!    LC: QA448.D38.
!
!  Parameters:
!
!    Input, real ( kind = rk ) XA, YA, XB, YB, XC, YC, XD, YD, the X and Y 
!    coordinates of the four vertices.
!
!    Output, logical VALUE, the value of the test.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  logical between
  logical intersect
  logical intersect_prop
  logical value
  real ( kind = rk ) xa
  real ( kind = rk ) xb
  real ( kind = rk ) xc
  real ( kind = rk ) xd
  real ( kind = rk ) ya
  real ( kind = rk ) yb
  real ( kind = rk ) yc
  real ( kind = rk ) yd

  if ( intersect_prop ( xa, ya, xb, yb, xc, yc, xd, yd ) ) then
    value = .true.
  else if ( between ( xa, ya, xb, yb, xc, yc ) ) then
    value = .true.
  else if ( between ( xa, ya, xb, yb, xd, yd ) ) then
    value = .true.
  else if ( between ( xc, yc, xd, yd, xa, ya ) ) then
    value = .true.
  else if ( between ( xc, yc, xd, yd, xb, yb ) ) then
    value = .true.
  else
    value = .false.
  end if

  intersect = value

  return
end
function intersect_prop ( xa, ya, xb, yb, xc, yc, xd, yd )

!*****************************************************************************80
!
!! intersect_prop() is TRUE if lines VA:VB and VC:VD have a proper intersection.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 May 2014
!
!  Author:
!
!    Original C version by Joseph ORourke.
!    This version by John Burkardt.
!
!  Reference:
!
!    Joseph ORourke,
!    Computational Geometry in C,
!    Cambridge, 1998,
!    ISBN: 0521649765,
!    LC: QA448.D38.
!
!  Parameters:
!
!    Input, real ( kind = rk ) XA, YA, XB, YB, XC, YC, XD, YD, the X and Y 
!    coordinates of the four vertices.
!
!    Output, logical INTERSECT_PROP, the result of the test.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  logical collinear
  logical intersect_prop
  logical l4_xor
  real ( kind = rk ) t1
  real ( kind = rk ) t2
  real ( kind = rk ) t3
  real ( kind = rk ) t4
  real ( kind = rk ) triangle_area
  logical value
  logical value1
  logical value2
  logical value3
  logical value4
  real ( kind = rk ) xa
  real ( kind = rk ) xb
  real ( kind = rk ) xc
  real ( kind = rk ) xd
  real ( kind = rk ) ya
  real ( kind = rk ) yb
  real ( kind = rk ) yc
  real ( kind = rk ) yd

  if ( collinear ( xa, ya, xb, yb, xc, yc ) ) then
    value = .false.
  else if ( collinear ( xa, ya, xb, yb, xd, yd ) ) then
    value = .false.
  else if ( collinear ( xc, yc, xd, yd, xa, ya ) ) then
    value = .false.
  else if ( collinear ( xc, yc, xd, yd, xb, yb ) ) then
    value = .false.
  else

    t1 = triangle_area ( xa, ya, xb, yb, xc, yc )
    t2 = triangle_area ( xa, ya, xb, yb, xd, yd )
    t3 = triangle_area ( xc, yc, xd, yd, xa, ya )
    t4 = triangle_area ( xc, yc, xd, yd, xb, yb )

    value1 = ( 0.0D+00 < t1 )
    value2 = ( 0.0D+00 < t2 )
    value3 = ( 0.0D+00 < t3 )
    value4 = ( 0.0D+00 < t4 )

    value = ( l4_xor ( value1, value2 ) ) .and. ( l4_xor ( value3, value4 ) )
 
  end if

  intersect_prop = value

  return
end
function l4_xor ( l1, l2 )

!*****************************************************************************80
!
!! l4_xor() returns the exclusive OR of two L4's.
!
!  Discussion:
!
!    An L4 is a logical value.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 May 2014
!
!  Author:
!
!   John Burkardt
!
!  Parameters:
!
!    Input, logical L1, L2, two values whose exclusive OR 
!    is needed.
!
!    Output, logical L4_XOR, the exclusive OR of L1 and L2.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  logical l1
  logical l2
  logical l4_xor
  logical value1
  logical value2

  value1 = (         l1   .and. ( .not. l2 ) )
  value2 = ( ( .not. l1 ) .and.         l2   )

  l4_xor = ( value1 .or. value2 )

  return
end
subroutine monomial_value ( m, n, e, x, value )

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
!    Output, real ( kind = rk ) VALUE(N), the value of the monomial.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer e(m)
  integer i
  real ( kind = rk ) value(n)
  real ( kind = rk ) x(m,n)

  value(1:n) = 1.0D+00

  do i = 1, m
    if ( 0 /= e(i) ) then
      value(1:n) = value(1:n) * x(i,1:n) ** e(i)
    end if
  end do

  return
end
function polygon_area ( nv, v )

!*****************************************************************************80
!
!! polygon_area() determines the area of a polygon.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 May 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NV, the number of vertices of the polygon.
!
!    Input, real ( kind = rk ) V(2,NV), the vertex coordinates.
!
!    Output, real ( kind = rk ) POLYGON_AREA, the area of the polygon.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nv

  real ( kind = rk ) area
  integer e(2)
  real ( kind = rk ) polygon_area
  real ( kind = rk ) v(2,nv)

  e(1) = 0
  e(2) = 0

  call polygon_monomial_integral ( nv, v, e, area )

  polygon_area = area

  return
end
function polygon_area2 ( n, x, y )

!*****************************************************************************80
!
!! polygon_area2() returns the area of a polygon.
!
!  Discussion:
!
!    The vertices should be listed in counter-clockwise order so that
!    the area will be positive.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 September 2016
!
!  Author:
!
!    John Burkardt.
!
!  Parameters:
!
!    Input, integer N, the number of vertices.
!
!    Input, real ( kind = rk ) X(N), Y(N), the vertex coordinates.
!
!    Output, real ( kind = rk ) POLYGON_AREA2, the area of the polygon.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) area
  integer i
  integer im1
  real ( kind = rk ) polygon_area2
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)

  area = 0.0D+00
  im1 = n

  do i = 1, n
    area = area + x(im1) * y(i) - x(i) * y(im1)
    im1 = i
  end do

  area = 0.5D+00 * area

  polygon_area2 = area

  return
end
subroutine polygon_monomial_integral ( nv, v, e, nu_pq )

!*****************************************************************************80
!
!! polygon_monomial_integral() integrates a monomial over a polygon.
!
!  Discussion:
!
!    Nu(P,Q) = Integral ( x, y in polygon ) x^e(1) y^e(2) dx dy
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
!    Input, integer NV, the number of vertices of the polygon.
!
!    Input, real ( kind = rk ) V(2,NV), the vertex coordinates.
!
!    Input, integer E(2), the exponents of the monomial.
!
!    Output, real ( kind = rk ) NU_PQ, the unnormalized moment Nu(P,Q).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nv

  integer e(2)
  integer i
  integer k
  integer l
  real ( kind = rk ) nu_pq
  integer p
  integer q
  real ( kind = rk ) r8_choose
  real ( kind = rk ) s_pq
  real ( kind = rk ) v(2,nv)
  real ( kind = rk ) xi
  real ( kind = rk ) xj
  real ( kind = rk ) yi
  real ( kind = rk ) yj

  p = e(1)
  q = e(2)

  nu_pq = 0.0D+00

  xj = v(1,nv)
  yj = v(2,nv)

  do i = 1, nv

    xi = v(1,i)
    yi = v(2,i)

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
subroutine polygon_sample ( nv, v, n, s )

!*****************************************************************************80
!
!! polygon_sample() uniformly samples a polygon.
!
!  Discussion:
!
!    The indexing of triangles() was changed to (nv-2,3), 15 January 2025.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 January 2025
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer NV, the number of vertices.
!
!    real ( kind = rk ) V(2,NV), the vertices of the polygon, listed in
!    counterclockwise order.
!
!    integer N, the number of points to create.
!
!  Output:
!
!    real ( kind = rk ) S(2,N), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nv

  real ( kind = rk ) area_cumulative(nv-2)
  real ( kind = rk ) area_polygon
  real ( kind = rk ) area_relative(nv-2)
  real ( kind = rk ) area_triangle(nv-2)
  real ( kind = rk ) area_percent
  integer i
  integer j
  integer k
  real ( kind = rk ) r(2)
  real ( kind = rk ) triangle_area
  integer triangles(nv-2,3)
  real ( kind = rk ) s(2,n)
  real ( kind = rk ) v(2,nv)
!
!  Triangulate the polygon.
!
  call polygon_triangulate ( nv, v(1,1:nv), v(2,1:nv), triangles )
!
!  Determine the areas of each triangle.
!
  do i = 1, nv - 2
    area_triangle(i) = triangle_area ( &
      v(1,triangles(i,1)), v(2,triangles(i,1)), &
      v(1,triangles(i,2)), v(2,triangles(i,2)), &
      v(1,triangles(i,3)), v(2,triangles(i,3)) )
  end do
!
!  Normalize the areas.
!
  area_polygon = sum ( area_triangle(1:nv-2) )
  area_relative(1:nv-2) = area_triangle(1:nv-2) / area_polygon
!
!  Replace each area by the sum of itself and all previous ones.
!
  area_cumulative(1) = area_relative(1)
  do i = 2, nv - 2
    area_cumulative(i) = area_relative(i) + area_cumulative(i-1)
  end do

  do j = 1, n
!
!  Choose triangle I at random, based on areas.
!
    call random_number ( harvest = area_percent )

    do k = 1, nv - 2

      i = k

      if ( area_percent <= area_cumulative(k) ) then
        exit
      end if

    end do
!
!  Now choose a point at random in triangle I.
!
    call random_number ( harvest = r(1:2) )

    if ( 1.0D+00 < sum ( r(1:2) ) ) then
      r(1:2) = 1.0D+00 - r(1:2)
    end if

    s(1:2,j) = ( 1.0D+00 - r(1) - r(2) ) * v(1:2,triangles(i,1)) &
                         + r(1)          * v(1:2,triangles(i,2)) &
                                + r(2)   * v(1:2,triangles(i,3))
  end do

  return
end
subroutine polygon_triangulate ( n, x, y, triangles )

!*****************************************************************************80
!
!! polygon_triangulate() determines a triangulation of a polygon.
!
!  Discussion:
!
!    There are N-3 triangles in the triangulation.
!
!    For the first N-2 triangles, the first edge listed is always an
!    internal diagonal.
!
!    Thanks to Gene Dial for pointing out a mistake in the area calculation,
!    10 September 2016.
!
!    Gene Dial requested an angle tolerance of about 1 millionth radian or 
!    5.7E-05 degrees, 26 June 2018.
!
!    The indexing of triangles() was changed to (nv-2,3), 15 January 2025.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 January 2025
!
!  Author:
!
!    Original C version by Joseph ORourke.
!    This version by John Burkardt.
!
!  Reference:
!
!    Joseph ORourke,
!    Computational Geometry in C,
!    Cambridge, 1998,
!    ISBN: 0521649765,
!    LC: QA448.D38.
!
!  Input:
!
!    integer N, the number of vertices.
!
!    real ( kind = rk ) X(N), Y(N), the coordinates of each vertex.
!
!  Output:
!
!    Ointeger TRIANGLES(N-2,3), the triangles of the triangulation.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) angle
  real ( kind = rk ) angle_degree
  real ( kind = rk ), parameter :: angle_tol = 5.7D-05
  real ( kind = rk ) area
  logical diagonal
  logical ear(n)
  integer i
  integer i0
  integer i1
  integer i2
  integer i3
  integer i4
  integer next_node(n)
  integer node
  integer node_m1
  integer node1
  integer node2
  integer node3
  real ( kind = rk ) polygon_area2
  integer prev_node(n)
  integer triangle_num
  integer triangles(n-2,3)
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)
!
!  We must have at least 3 vertices.
!
  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POLYGON_TRIANGULATE - Fatal error!'
    write ( *, '(a)' ) '  N < 3.'
    stop 1
  end if
!
!  Consecutive vertices cannot be equal.
!
  node_m1 = n
  do node = 1, n
    if ( x(node_m1) == x(node) .and. y(node_m1) == y(node) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POLYGON_TRIANGULATE - Fatal error!'
      write ( *, '(a)' ) '  Two consecutive nodes are identical.'
      stop 1
    end if
    node_m1 = node
  end do
!
!  No node can be the vertex of an angle less than 1 degree 
!  in absolute value.
!
  node1 = n

  do node2 = 1, n

    node3 = mod ( node2, n ) + 1

    angle = angle_degree ( &
      x(node1), y(node1), &
      x(node2), y(node2), &
      x(node3), y(node3) )

    if ( abs ( angle ) <= angle_tol ) then
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'POLYGON_TRIANGULATE - Fatal error!'
      write ( *, '(a,g14.6)' ) '  Polygon has an angle smaller than ', angle_tol
      write ( *, '(a,i4)' ) '  occurring at node ', node2
      stop 1
    end if

    node1 = node2

  end do
!
!  Area must be positive.
!
  area = polygon_area2 ( n, x, y )

  if ( area <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POLYGON_TRIANGULATE - Fatal error!'
    write ( *, '(a)' ) '  Polygon has zero or negative area.'
    stop 1
  end if
!
!  PREV and NEXT point to the previous and next nodes.
!
  i = 1
  prev_node(i) = n
  next_node(i) = i + 1

  do i = 2, n - 1
    prev_node(i) = i - 1
    next_node(i) = i + 1
  end do

  i = n
  prev_node(i) = i - 1
  next_node(i) = 1
!
!  EAR indicates whether the node and its immediate neighbors form an ear
!  that can be sliced off immediately.
!
  do i = 1, n
    ear(i) = diagonal ( prev_node(i), next_node(i), n, prev_node, next_node, &
      x, y )
  end do

  triangle_num = 0

  i2 = 1

  do while ( triangle_num < n - 3 )
!
!  If I2 is an ear, gather information necessary to carry out
!  the slicing operation and subsequent "healing".
!
    if ( ear(i2) ) then

      i3 = next_node(i2)
      i4 = next_node(i3)
      i1 = prev_node(i2)
      i0 = prev_node(i1)
!
!  Make vertex I2 disappear.
!
      next_node(i1) = i3
      prev_node(i3) = i1
!
!  Update the earity of I1 and I3, because I2 disappeared.
!
      ear(i1) = diagonal ( i0, i3, n, prev_node, next_node, x, y )
      ear(i3) = diagonal ( i1, i4, n, prev_node, next_node, x, y )
!
!  Add the diagonal [I3, I1, I2] to the list.
!
      triangle_num = triangle_num + 1
      triangles(triangle_num,1) = i3
      triangles(triangle_num,2) = i1
      triangles(triangle_num,3) = i2

    end if
!
!  Try the next vertex.
!
    i2 = next_node(i2)

  end do
!
!  The last triangle is formed from the three remaining vertices.
!
  i3 = next_node(i2)
  i1 = prev_node(i2)

  triangle_num = triangle_num + 1
  triangles(triangle_num,1) = i3
  triangles(triangle_num,2) = i1
  triangles(triangle_num,3) = i2

  return
end
function r8_choose ( n, k )

!*****************************************************************************80
!
!! r8_choose() computes the binomial coefficient C(N,K) as an R8.
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
function triangle_area ( xa, ya, xb, yb, xc, yc )

!*****************************************************************************80
!
!! triangle_area() computes the signed area of a triangle.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 May 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) XA, YA, XB, YB, XC, YC, the coordinates of
!    the vertices of the triangle, given in counterclockwise order.
!
!    Output, real ( kind = rk ) TRIANGLE_AREA, the signed area of the triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) triangle_area
  real ( kind = rk ) value
  real ( kind = rk ) xa
  real ( kind = rk ) xb
  real ( kind = rk ) xc
  real ( kind = rk ) ya
  real ( kind = rk ) yb
  real ( kind = rk ) yc

  value = 0.5D+00 * ( &
      ( xb - xa ) * ( yc - ya ) &
    - ( xc - xa ) * ( yb - ya ) )

  triangle_area = value

  return
end

