subroutine angle_half ( p1, p2, p3, p4 )

!*****************************************************************************80
!
!! angle_half() finds half an angle.
!
!  Discussion:
!
!    The original angle is defined by the sequence of points P1, P2 and P3.
!
!    The point P4 is calculated so that:
!
!      (P1,P2,P4) = (P1,P2,P3) / 2
!
!        P1
!        /
!       /   P4
!      /  .  
!     / .
!    P2--------->P3
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) P1(2), P2(2), P3(2), points defining the angle. 
!
!    real ( kind = rk ) P4(2), a point defining the half angle.
!    The vector P4 - P2 will have unit norm.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) p1(2)
  real ( kind = rk ) p2(2)
  real ( kind = rk ) p3(2)
  real ( kind = rk ) p4(2)

  p4(1:2) = 0.5D+00 * ( &
      ( p1(1:2) - p2(1:2) ) / sqrt ( sum ( ( p1(1:2) - p2(1:2) )**2 ) ) &
    + ( p3(1:2) - p2(1:2) ) / sqrt ( sum ( ( p3(1:2) - p2(1:2) )**2 ) ) )

   p4(1:2) = p2(1:2) + p4(1:2) / sqrt ( sum ( p4(1:2)**2 ) )

  return
end
function angle_rad ( p1, p2, p3 )

!*****************************************************************************80
!
!! angle_rad() returns the angle in radians swept out between two rays.
!
!  Discussion:
!
!    Except for the zero angle case, it should be true that
!
!      ANGLE_RAD ( P1, P2, P3 ) + ANGLE_RAD ( P3, P2, P1 ) = 2 * PI
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
!    15 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) P1(2), P2(2), P3(2), define the rays
!    P1 - P2 and P3 - P2 which define the angle.
!
!  Output:
!
!    real ( kind = rk ) ANGLE_RAD, the angle swept out by the rays,
!    in radians.  0 <= ANGLE_RAD < 2 * PI.  If either ray has zero
!    length, then ANGLE_RAD is set to 0.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) angle_rad
  real ( kind = rk ) p(2)
  real ( kind = rk ) p1(2)
  real ( kind = rk ) p2(2)
  real ( kind = rk ) p3(2)
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00

  p(1) = ( p3(1) - p2(1) ) * ( p1(1) - p2(1) ) &
       + ( p3(2) - p2(2) ) * ( p1(2) - p2(2) )

  p(2) = ( p3(1) - p2(1) ) * ( p1(2) - p2(2) ) &
       - ( p3(2) - p2(2) ) * ( p1(1) - p2(1) )

  if ( all ( p(1:2) == 0.0D+00)  ) then
    angle_rad = 0.0D+00
    return
  end if

  angle_rad = atan2 ( p(2), p(1) )

  if ( angle_rad < 0.0D+00 ) then
    angle_rad = angle_rad + 2.0D+00 * r8_pi
  end if

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
!  Input:
!
!    real ( kind = rk ) XA, YA, XB, YB, XC, YC, the coordinates of 
!    the vertices.
!
!  Output:
!
!    logical BETWEEN, is TRUE if C is between A and B.
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
!  Input:
!
!    real ( kind = rk ) XA, YA, XB, YB, XC, YC, the coordinates of 
!    the vertices.
!
!  Output:
!
!    logical COLLINEAR, is TRUE if the points are judged 
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
  logical value
  real ( kind = rk ) xa
  real ( kind = rk ) xb
  real ( kind = rk ) xc
  real ( kind = rk ) ya
  real ( kind = rk ) yb
  real ( kind = rk ) yc

  area = 0.5D+00 * ( &
      ( xb - xa ) * ( yc - ya ) &
    - ( xc - xa ) * ( yb - ya ) )

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
function diagonal ( im1, ip1, n, prev, next, x, y )

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
!  Input:
!
!    integer IM1, IP1, the indices of two vertices.
!
!    integer N, the number of vertices.
!
!    integer PREV(N), the previous neighbor of each vertex.
!
!    integer NEXT(N), the next neighbor of each vertex.
!
!    real ( kind = rk ) X(N), Y(N), the coordinates of each vertex.
!
!  Output:
!
!    logical DIAGONAL, the value of the test.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  logical diagonal
  logical diagonalie
  integer im1
  logical in_cone
  integer ip1
  integer next(n)
  integer prev(n)
  logical value1
  logical value2
  logical value3
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)

  value1 = in_cone ( im1, ip1, n, prev, next, x, y )
  value2 = in_cone ( ip1, im1, n, prev, next, x, y )
  value3 = diagonalie ( im1, ip1, n, next, x, y )

  diagonal = ( value1 .and. value2 .and. value3 )

  return
end
function diagonalie ( im1, ip1, n, next, x, y )

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
!  Input:
!
!    integer IM1, IP1, the indices of two vertices.
!
!    integer N, the number of vertices.
!
!    integer NEXT(N), the next neighbor of each vertex.
!
!    real ( kind = rk ) X(N), Y(N), the coordinates of each vertex.
!
!  Output:
!
!    logical DIAGONALIE, the value of the test.
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
  integer next(n)
  logical value
  logical value2
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)

  first = im1
  j = first
  jp1 = next(first)

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
    jp1 = next(j)

    if ( j == first ) then
      exit
    end if

  end do

  diagonalie = value

  return
end
function i4_modp ( i, j )

!*****************************************************************************80
!
!! i4_modp() returns the nonnegative remainder of I4 division.
!
!  Discussion:
!
!    If
!      NREM = I4_MODP ( I, J )
!      NMULT = ( I - NREM ) / J
!    then
!      I = J * NMULT + NREM
!    where NREM is always nonnegative.
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
!
!    An I4 is an integer value.
!
!  Example:
!
!        I     J     MOD I4_MODP    Factorization
!
!      107    50       7       7    107 =  2 *  50 + 7
!      107   -50       7       7    107 = -2 * -50 + 7
!     -107    50      -7      43   -107 = -3 *  50 + 43
!     -107   -50      -7      43   -107 =  3 * -50 + 43
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer I, the number to be divided.
!
!    integer J, the number that divides I.
!
!  Output:
!
!    integer I4_MODP, the nonnegative remainder when I is divided by J.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer i4_modp
  integer j
  integer value

  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MODP(): Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal divisor J = ', j
    stop 1
  end if

  value = mod ( i, j )

  if ( value < 0 ) then
    value = value + abs ( j )
  end if

  i4_modp = value

  return
end
function i4_wrap ( ival, ilo, ihi )

!*****************************************************************************80
!
!! i4_wrap() forces an I4 to lie between given limits by wrapping.
!
!  Discussion:
!
!    An I4 is an integer value.
!
!    There appears to be a bug in the GFORTRAN compiler which can lead to
!    erroneous results when the first argument of I4_WRAP is an expression.
!    In particular:
!
!    do i = 1, 3
!      if ( test ) then
!        i4 = i4_wrap ( i + 1, 1, 3 )
!      end if
!    end do
!
!    was, when I = 3, returning I4 = 3.  So I had to replace this with
!
!    do i = 1, 3
!      if ( test ) then
!        i4 = i + 1
!        i4 = i4_wrap ( i4, 1, 3 )
!      end if
!    end do
!
!  Example:
!
!    ILO = 4, IHI = 8
!
!    I  Value
!
!    -2     8
!    -1     4
!     0     5
!     1     6
!     2     7
!     3     8
!     4     4
!     5     5
!     6     6
!     7     7
!     8     8
!     9     4
!    10     5
!    11     6
!    12     7
!    13     8
!    14     4
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer IVAL, a value.
!
!    integer ILO, IHI, the desired bounds.
!
!  Output:
!
!    integer I4_WRAP, a "wrapped" version of the value.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i4_modp
  integer i4_wrap
  integer ihi
  integer ilo
  integer ival
  integer jhi
  integer jlo
  integer value
  integer wide

  jlo = min ( ilo, ihi )
  jhi = max ( ilo, ihi )

  wide = jhi - jlo + 1

  if ( wide == 1 ) then
    value = jlo
  else
    value = jlo + i4_modp ( ival - jlo, wide )
  end if

  i4_wrap = value

  return
end
function in_cone ( im1, ip1, n, prev, next, x, y )

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
!  Input:
!
!    integer IM1, IP1, the indices of two vertices.
!
!    integer N, the number of vertices.
!
!    integer PREV(N), the previous neighbor of each vertex.
!
!    integer NEXT(N), the next neighbor of each vertex.
!
!    real ( kind = rk ) X(N), Y(N), the coordinates of each vertex.
!
!  Output:
!
!    logical IN_CONE, the value of the test.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer im1
  integer im2
  logical in_cone
  integer ip1
  integer next(n)
  integer prev(n)
  real ( kind = rk ) t1
  real ( kind = rk ) t2
  real ( kind = rk ) t3
  real ( kind = rk ) t4
  real ( kind = rk ) t5
  real ( kind = rk ) triangle_area
  logical value
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)

  im2 = prev(im1)
  i = next(im1)

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
!  Input:
!
!    real ( kind = rk ) XA, YA, XB, YB, XC, YC, XD, YD, the X and Y 
!    coordinates of the four vertices.
!
!  Output:
!
!    logical VALUE, the value of the test.
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

  if ( intersect_prop ( xa, ya, xb, yb, xc, yc, xc, yd ) ) then
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
!  Input:
!
!    real ( kind = rk ) XA, YA, XB, YB, XC, YC, XD, YD, the X and Y 
!    coordinates of the four vertices.
!
!  Output:
!
!    logical INTERSECT_PROP, the result of the test.
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
!  Input:
!
!    logical L1, L2, two values whose exclusive OR 
!    is needed.
!
!  Output:
!
!    logical L4_XOR, the exclusive OR of L1 and L2.
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
subroutine polygon_angles ( n, v, angle )

!*****************************************************************************80
!
!! polygon_angles() computes the interior angles of a polygon.
!
!  Discussion:
!
!    The vertices should be listed in counter clockwise order.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of vertices of the polygon.
!
!    real ( kind = rk ) V(2,N), the vertices.
!
!  Output:
!
!    real ( kind = rk ) ANGLE(N), the angles of the polygon,
!    in radians.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) angle(n)
  real ( kind = rk ) angle_rad
  integer i
  integer i4_wrap
  integer im1
  integer ip1
  real ( kind = rk ) v(2,n)

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'polygon_angles(): Fatal error!'
    write ( *, '(a)' ) '  The number of vertices must be at least 3.'
    write ( *, '(a,i8)' ) '  The input value of N = ', n
    stop 1
  end if
 
  do i = 1, n

    im1 = i4_wrap ( i - 1, 1, n )
    ip1 = i4_wrap ( i + 1, 1, n )

    angle(i) = angle_rad ( v(1:2,im1), v(1:2,i), v(1:2,ip1) )

  end do

  return
end
subroutine polygon_area ( n, v, area )

!*****************************************************************************80
!
!! polygon_area() computes the area of a polygon.
!
!  Discussion:
!
!    AREA = 1/2 * abs ( sum ( 1 <= I <= N ) X(I) * ( Y(I+1) - Y(I-1) ) )
!    where Y(0) should be replaced by Y(N), and Y(N+1) by Y(1).
!
!    If the vertices are given in counter clockwise order, the area
!    will be positive.  If the vertices are given in clockwise order,
!    the area will be negative.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of vertices of the polygon.
!
!    real ( kind = rk ) V(2,N), the vertices.
!
!  Output:
!
!    real ( kind = rk ) AREA, the absolute area of the polygon.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) area
  integer i
  integer i4_wrap
  integer im1
  integer ip1
  real ( kind = rk ) v(2,n)

  area = 0.0D+00

  do i = 1, n

    im1 = i4_wrap ( i-1, 1, n )
    ip1 = i4_wrap ( i+1, 1, n )

    area = area + v(1,i) * ( v(2,ip1) - v(2,im1) )

  end do

  area = 0.5D+00 * area

  return
end
subroutine polygon_area_2 ( n, v, area )

!*****************************************************************************80
!
!! polygon_area_2() computes the area of a polygon.
!
!  Discussion:
!
!    The area is the sum of the areas of the triangles formed by
!    node N with consecutive pairs of nodes.
!
!    If the vertices are given in counter clockwise order, the area
!    will be positive.  If the vertices are given in clockwise order,
!    the area will be negative.
!
!    Thanks to Martin Pineault for noticing that an earlier version
!    of this routine would not correctly compute the area of a nonconvex
!    polygon.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Adrian Bowyer, John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983,
!    ISBN: 0408012420.
!
!  Input:
!
!    integer N, the number of vertices of the polygon.
!
!    real ( kind = rk ) V(2,N), the vertices.
!
!  Output:
!
!    real ( kind = rk ) AREA, the absolute area of the polygon.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) area
  real ( kind = rk ) area_triangle
  integer i
  real ( kind = rk ) triangle_area
  real ( kind = rk ) v(2,n)

  area = 0.0D+00

  do i = 1, n - 2

    area_triangle = triangle_area ( &
      v(1,i),   v(2,i), &
      v(1,i+1), v(2,i+1), &
      v(1,n),   v(2,n) )

    area = area + area_triangle

  end do

  return
end
subroutine polygon_area_lattice ( i, b, area )

!*****************************************************************************80
!
!! polygon_area_lattice() computes the area of a lattice polygon.
!
!  Discussion:
!
!    We define a lattice to be the 2D plane, in which the points
!    whose (X,Y) coordinates are both integers are given a special
!    status as "lattice points".
!
!    A lattice polygon is a polygon whose vertices are lattice points.
!
!    The area of a lattice polygon can be computed by Pick's Theorem:
!
!      Area = I + B / 2 - 1
!
!    where
!
!      I = the number of lattice points contained strictly inside the polygon;
!
!      B = the number of lattice points that lie exactly on the boundary.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Branko Gruenbaum, Geoffrey Shephard,
!    Pick's Theorem,
!    The American Mathematical Monthly,
!    Volume 100, Number 2, February 1993, pages 150-161.
!
!  Input:
!
!    integer I, the number of interior lattice points.
!
!    integer B, the number of boundary lattice points.
!
!  Output:
!
!    real ( kind = rk ) AREA, the area of the lattice polygon.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) area
  integer b
  integer i

  area = real ( i, kind = rk ) + real ( b, kind = rk ) / 2.0D+00 - 1.0D+00

  return
end
subroutine polygon_centroid ( n, v, centroid )

!*****************************************************************************80
!
!! polygon_centroid() computes the centroid of a polygon.
!
!  Discussion:
!
!    Denoting the centroid coordinates by CENTROID, then
!
!      CENTROID(1) = Integral ( Polygon interior ) x dx dy / Area ( Polygon )
!      CENTROID(2) = Integral ( Polygon interior ) y dx dy / Area ( Polygon ).
!
!    Green's theorem states that for continuously differentiable functions
!    M(x,y) and N(x,y),
!
!      Integral ( Polygon boundary ) ( M dx + N dy ) =
!      Integral ( Polygon interior ) ( dN/dx - dM/dy ) dx dy.
!
!    Using M(x,y) = 0 and N(x,y) = x*x/2, we get:
!
!      CENTROID(1) = 0.5 * Integral ( Polygon boundary ) x*x dy 
!                  / Area ( Polygon ),
!
!    which becomes
!
!      CENTROID(1) = 1/6 sum ( 1 <= I <= N )
!        ( X(I+1) + X(I) ) * ( X(I) * Y(I+1) - X(I+1) * Y(I))
!        / Area ( Polygon )
!
!    where, when I = N, the index "I+1" is replaced by 1.
!
!    A similar calculation gives us a formula for CENTROID(2).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 July 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Gerard Bashein, Paul Detmer,
!    Centroid of a Polygon,
!    in Graphics Gems IV, 
!    edited by Paul Heckbert,
!    AP Professional, 1994,
!    T385.G6974.
!
!  Input:
!
!    integer N, the number of sides of the polygon.
!
!    real ( kind = rk ) V(2,N), the coordinates of the vertices.
!
!  Output:
!
!    real ( kind = rk ) CENTROID(2), the coordinates of the centroid.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) area
  real ( kind = rk ) centroid(2)
  integer i
  integer ip1
  real ( kind = rk ) temp
  real ( kind = rk ) v(2,n)

  area = 0.0D+00
  centroid(1:2) = 0.0D+00

  do i = 1, n

    if ( i < n ) then
      ip1 = i + 1
    else
      ip1 = 1
    end if

    temp = ( v(1,i) * v(2,ip1) - v(1,ip1) * v(2,i) )

    area = area + temp

    centroid(1:2) = centroid(1:2) + ( v(1:2,ip1) + v(1:2,i) ) * temp

  end do

  area = area / 2.0D+00

  if ( area == 0.0D+00 ) then
    centroid(1:2) = v(1:2,1)
  else
    centroid(1:2) = centroid(1:2) / ( 6.0D+00 * area )
  end if

  return
end
subroutine polygon_centroid_2 ( n, v, centroid )

!*****************************************************************************80
!
!! polygon_centroid_2() computes the centroid of a polygon.
!
!  Discussion:
!
!    The centroid is the area-weighted sum of the centroids of
!    disjoint triangles that make up the polygon.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 July 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Adrian Bowyer, John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983,
!    ISBN: 0408012420.
!
!  Input:
!
!    integer N, the number of vertices of the polygon.
!
!    real ( kind = rk ) V(2,N), the coordinates of the vertices.
!
!  Output:
!
!    real ( kind = rk ) CENTROID(2), the coordinates of the centroid.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) area_polygon
  real ( kind = rk ) area_triangle
  real ( kind = rk ) centroid(2)
  integer i
  real ( kind = rk ) triangle_area
  real ( kind = rk ) v(2,n)

  area_polygon = 0.0D+00
  centroid(1:2) = 0.0D+00

  do i = 1, n - 2

    area_triangle = triangle_area ( &
      v(1,i),   v(2,i), &
      v(1,i+1), v(2,i+1), &
      v(1,n),   v(2,n) )

    area_polygon = area_polygon + area_triangle

    centroid(1:2) = centroid(1:2) + area_triangle &
      * ( v(1:2,i) + v(1:2,i+1) + v(1:2,n) ) / 3.0D+00

  end do

  if ( area_polygon == 0.0D+00 ) then
    centroid(1:2) = v(1:2,1)
  else
    centroid(1:2) = centroid(1:2) / area_polygon
  end if

  return
end
subroutine polygon_contains_point ( n, v, p, inside )

!*****************************************************************************80
!
!! polygon_contains_point() finds if a point is inside a polygon.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 November 2016
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of nodes or vertices in 
!    the polygon.  N must be at least 3.
!
!    real ( kind = rk ) V(2,N), the vertices of the polygon.
!
!    real ( kind = rk ) P(2), the coordinates of the point to be tested.
!
!  Output:
!
!    logical INSIDE, is TRUE if the point is inside 
!    the polygon.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  logical inside
  real ( kind = rk ) p(2)
  real ( kind = rk ) px1
  real ( kind = rk ) px2
  real ( kind = rk ) py1
  real ( kind = rk ) py2
  real ( kind = rk ) v(2,n)
  real ( kind = rk ) xints

  inside = .false.

  px1 = v(1,1)
  py1 = v(2,1)
  xints = p(1) - 1.0D+00

  do i = 1, n

    px2 = v(1,mod(i,n)+1)
    py2 = v(2,mod(i,n)+1)

    if ( min ( py1, py2 ) < p(2) ) then
      if ( p(2) <= max ( py1, py2 ) ) then
        if ( p(1) <= max ( px1, px2 ) ) then
          if ( py1 /= py2 ) then
            xints = ( p(2) - py1 ) * ( px2 - px1 ) / ( py2 - py1 ) + px1
          end if
          if ( px1 == px2 .or. p(1) <= xints ) then
            inside = .not. inside
          end if
        end if
      end if
    end if

    px1 = px2
    py1 = py2

  end do

  return
end
subroutine polygon_convex_contains_point ( n, v, p, inside )

!*****************************************************************************80
!
!! polygon_convex_contains_point(): is a point inside a convex polygon.
!
!  Discussion:
!
!    The result is not reliable if the polygon is not convex.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of nodes or vertices in the 
!    polygon.  N must be at least 3.
!
!    real ( kind = rk ) V(2,N), the vertices of the polygon.
!
!    real ( kind = rk ) P(2), the coordinates of the point to be tested.
!
!  Output:
!
!    logical INSIDE, is TRUE if the point is inside
!    the polygon or on its boundary.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  logical inside
  real ( kind = rk ) p(2)
  real ( kind = rk ) t(2,3)
  real ( kind = rk ) v(2,n)

  inside = .false.
!
!  A point is inside a convex polygon if and only if it is inside
!  one of the triangles formed by X(1),Y(1) and any two consecutive
!  points on the polygon's circumference.
!
  t(1:2,1) = v(1:2,1)

  do i = 2, n - 1

    t(1:2,2) = v(1:2,i)
    t(1:2,3) = v(1:2,i+1)

    call triangle_contains_point_1 ( t, p, inside )

    if ( inside ) then
      return
    end if

  end do

  return
end
subroutine polygon_contains_point_3 ( n, v, p, inside )

!*****************************************************************************80
!
!! polygon_contains_point_3() finds if a point is inside a simple polygon.
!
!  Discussion:
!
!    A simple polygon is one whose boundary never crosses itself.
!    The polygon does not need to be convex.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Moshe Shimrat,
!    ACM Algorithm 112,
!    Position of Point Relative to Polygon,
!    Communications of the ACM,
!    Volume 5, Number 8, page 434, August 1962.
!
!  Input:
!
!    integer N, the number of nodes or vertices in 
!    the polygon.  N must be at least 3.
!
!    real ( kind = rk ) V(2,N), the vertices of the polygon.
!
!    real ( kind = rk ) P(2), the coordinates of the point to be tested.
!
!  Output:
!
!    logical INSIDE, is TRUE if the point is inside the polygon.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  logical inside
  integer ip1
  real ( kind = rk ) p(2)
  real ( kind = rk ) v(2,n)

  inside = .false.

  do i = 1, n

    if ( i < n ) then
      ip1 = i + 1
    else
      ip1 = 1
    end if

    if ( ( v(2,i)   <  p(2) .and. p(2) <= v(2,ip1)   ) .or. &
         ( p(2) <= v(2,i)   .and. v(2,ip1)   < p(2) ) ) then
      if ( ( p(1) - v(1,i) ) - ( p(2) - v(2,i) ) &
         * ( v(1,ip1) - v(1,i) ) / ( v(2,ip1) - v(2,i) ) < 0.0D+00 ) then
        inside = .not. inside
      end if
    end if

  end do

  return
end
subroutine polygon_data_inrad ( n, radin, area, radout, side )

!*****************************************************************************80
!
!! polygon_data_inrad() determines polygonal data from its inner radius.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of sides of the polygon.
!    N must be at least 3.
!
!    real ( kind = rk ) RADIN, the inner radius of the polygon, that is,
!    the radius of the largest circle that can be inscribed within
!    the polygon.
!
!  Output:
!
!    real ( kind = rk ) AREA, the area of the regular polygon.
!
!    real ( kind = rk ) RADOUT, the outer radius of the polygon, that is,
!    the radius of the smallest circle that can be described about
!    the polygon.
!
!    real ( kind = rk ) SIDE, the length of one side of the polygon.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) angle
  real ( kind = rk ) area
  integer n
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) radin
  real ( kind = rk ) radout
  real ( kind = rk ) side

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'polygon_data_inrad(): Fatal error!'
    write ( *, '(a)' ) '  Input value of N must be at least 3'
    write ( *, '(a,i8)' ) '  but your input value was N = ', n
    stop 1
  end if

  angle = r8_pi / real ( n, kind = rk )
  area = real ( n, kind = rk ) * radin * radin * tan ( angle )
  side = 2.0D+00 * radin * tan ( angle )
  radout = 0.5D+00 * side / sin ( angle )

  return
end
subroutine polygon_data_outrad ( n, radout, area, radin, side )

!*****************************************************************************80
!
!! polygon_data_outrad() determines polygonal data from its outer radius.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of sides of the polygon.
!    N must be at least 3.
!
!    real ( kind = rk ) RADOUT, the outer radius of the polygon, that is,
!    the radius of the smallest circle that can be described
!    around the polygon.
!
!  Output:
!
!    real ( kind = rk ) AREA, the area of the regular polygon.
!
!    real ( kind = rk ) RADIN, the inner radius of the polygon, that is,
!    the radius of the largest circle that can be inscribed
!    within the polygon.
!
!    real ( kind = rk ) SIDE, the length of one side of the polygon.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) angle
  real ( kind = rk ) area
  integer n
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) radin
  real ( kind = rk ) radout
  real ( kind = rk ) side

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'polygon_data_outrad(): Fatal error!'
    write ( *, '(a)' ) '  Input value of N must be at least 3'
    write ( *, '(a,i8)' ) '  but your input value was N = ', n
    stop 1
  end if

  angle = r8_pi / real ( n, kind = rk )
  area = 0.5D+00 * real ( n, kind = rk ) * radout * radout &
    * sin ( 2.0D+00 * angle )
  side = 2.0D+00 * radout * sin ( angle )
  radin = 0.5D+00 * side / tan ( angle )

  return
end
subroutine polygon_data_side ( n, side, area, radin, radout )

!*****************************************************************************80
!
!! polygon_data_side() determines polygonal data from its side length.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    29 June 2005
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of sides of the polygon.
!    N must be at least 3.
!
!    real ( kind = rk ) SIDE, the length of one side of the polygon.
!
!  Output:
!
!    real ( kind = rk ) AREA, the area of the regular polygon.
!
!    real ( kind = rk ) RADIN, the inner radius of the polygon, that is,
!    the radius of the largest circle that can be inscribed within
!    the polygon.
!
!    real ( kind = rk ) RADOUT, the outer radius of the polygon, that is,
!    the radius of the smallest circle that can be described about
!    the polygon.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) angle
  real ( kind = rk ) area
  integer n
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) radin
  real ( kind = rk ) radout
  real ( kind = rk ) side

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'polygon_data_side(): Fatal error!'
    write ( *, '(a)' ) '  Input value of N must be at least 3'
    write ( *, '(a,i8)' ) '  but your input value was N = ', n
    stop 1
  end if

  angle = r8_pi / real ( n, kind = rk )
  area = 0.25D+00 * real ( n, kind = rk ) * side * side / tan ( angle )
  radin = 0.5D+00 * side / tan ( angle )
  radout = 0.5D+00 * side / sin ( angle )

  return
end
subroutine polygon_diameter ( n, v, diameter )

!*****************************************************************************80
!
!! polygon_diameter() computes the diameter of a polygon.
!
!  Discussion:
!
!    The diameter of a polygon is the maximum distance between any
!    two points on the polygon.  It is guaranteed that this maximum
!    distance occurs between two vertices of the polygon.  It is
!    sufficient to check the distance between all pairs of vertices.
!    This is an N^2 algorithm.  There is an algorithm by Shamos which
!    can compute this quantity in order N time instead.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of vertices of the polygon.
!
!    real ( kind = rk ) V(2,N), the vertices.
!
!  Output:
!
!    real ( kind = rk ) DIAMETER, the diameter of the polygon.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) diameter
  integer i
  integer j
  real ( kind = rk ) v(2,n)

  diameter = 0.0D+00

  do i = 1, n

    do j = i + 1, n
      diameter = max ( diameter, &
        sqrt ( ( v(1,i) - v(1,j) ) ** 2 + ( v(2,i) - v(2,j) ) ** 2 ) )
    end do

  end do

  return
end
subroutine polygon_expand ( n, v, h, w )

!*****************************************************************************80
!
!! polygon_expand() expands a polygon.
!
!  Discussion:
!
!    This routine simple moves each vertex of the polygon outwards
!    in such a way that the sides of the polygon advance by H.  
!
!    This approach should always work if the polygon is convex, or 
!    star-shaped.  But for general polygons, it is possible
!    that this procedure, for large enough H, will create a polygon
!    whose sides intersect.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of sides of the polygon.
!
!    real ( kind = rk ) V(2,N), the coordinates of the vertices.
!
!    real ( kind = rk ) H, the expansion amount.
!
!  Output:
!
!    real ( kind = rk ) W(2,N), the "expanded" coordinates.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) angle
  real ( kind = rk ) angle_rad
  real ( kind = rk ) h
  real ( kind = rk ) h2
  integer i
  integer i4_wrap
  integer im1
  integer ip1
  real ( kind = rk ) p4(2)
  real ( kind = rk ) v(2,n)
  real ( kind = rk ) w(2,n)
!
!  Consider each angle, formed by the nodes P(I-1), P(I), P(I+1).
!
  do i = 1, n

    im1 = i4_wrap ( i-1, 1, n )
    ip1 = i4_wrap ( i+1, 1, n )
!
!        P1
!        /
!       /   P4
!      /  .  
!     / .
!    P2--------->P3
!
    call angle_half ( v(1:2,im1), v(1:2,i), v(1:2,ip1), p4 )
!
!  Compute the value of the half angle.
!
    angle = angle_rad ( v(1:2,im1), v(1:2,i), p4(1:2) )
!
!  The stepsize along the ray must be adjusted so that the sides
!  move out by H.
!
    h2 = h / sin ( angle )

    w(1:2,i) = v(1:2,i) - h2 * ( p4(1:2) - v(1:2,i) )

  end do

  return
end
subroutine polygon_integral_1 ( n, v, result )

!*****************************************************************************80
!
!! polygon_integral_1() integrates the function 1 over a polygon.
!
!  Discussion:
!
!    The polygon is bounded by the points (X(1:N), Y(1:N)).
!
!    INTEGRAL = 0.5 * sum ( 1 <= I <= N )
!      ( X(I) + X(I-1) ) * ( Y(I) - Y(I-1) )
!
!    where X(0) and Y(0) should be replaced by X(N) and Y(N).
!
!    The integral of 1 over a polygon is the area of the polygon.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    SF Bockman,
!    Generalizing the Formula for Areas of Polygons to Moments,
!    American Mathematical Society Monthly,
!    1989, pages 131-132.
!
!  Input:
!
!    integer N, the number of vertices of the polygon.
!    N should be at least 3 for a nonzero result.
!
!    real ( kind = rk ) V(2,N), the coordinates of the vertices
!    of the polygon.  These vertices should be given in counter clockwise order.
!
!  Output:
!
!    real ( kind = rk ) RESULT, the value of the integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer im1
  real ( kind = rk ) result
  real ( kind = rk ) v(2,n)

  result = 0.0D+00

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'polygon_integral_1(): Fatal error!'
    write ( *, '(a)' ) '  The number of vertices must be at least 3.'
    write ( *, '(a,i8)' ) '  The input value of N = ', n
    stop 1
  end if

  do i = 1, n

    if ( i == 1 ) then
      im1 = n
    else
      im1 = i - 1
    end if

    result = result + 0.5D+00 * ( v(1,i) + v(1,im1) ) * ( v(2,i) - v(2,im1) )

  end do

  return
end
subroutine polygon_integral_x ( n, v, result )

!*****************************************************************************80
!
!! polygon_integral_x() integrates the function X over a polygon.
!
!  Discussion:
!
!    The polygon is bounded by the points (X(1:N), Y(1:N)).
!
!    INTEGRAL = (1/6) * sum ( 1 <= I <= N )
!      ( X(I)*X(I) + X(I) * X(I-1) + X(I-1)*X(I-1) ) * ( Y(I) - Y(I-1) )
!
!    where X(0) and Y(0) should be replaced by X(N) and Y(N).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    SF Bockman,
!    Generalizing the Formula for Areas of Polygons to Moments,
!    American Mathematical Society Monthly,
!    1989, pages 131-132.
!
!  Input:
!
!    integer N, the number of vertices of the polygon.
!    N should be at least 3 for a nonzero result.
!
!    real ( kind = rk ) V(2,N), the coordinates of the vertices
!    of the polygon.  These vertices should be given in counter clockwise order.
!
!  Output:
!
!    real ( kind = rk ) RESULT, the value of the integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer im1
  real ( kind = rk ) result
  real ( kind = rk ) v(2,n)

  result = 0.0D+00

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'polygon_integral_x(): Fatal error!'
    write ( *, '(a)' ) '  The number of vertices must be at least 3.'
    write ( *, '(a,i8)' ) '  The input value of N = ', n
    stop 1
  end if

  do i = 1, n

    if ( i == 1 ) then
      im1 = n
    else
      im1 = i - 1
    end if

    result = result + ( v(1,i)**2 + v(1,i) * v(1,im1) + v(1,im1)**2 ) &
      * ( v(2,i) - v(2,im1) )

  end do

  result = result / 6.0D+00

  return
end
subroutine polygon_integral_xx ( n, v, result )

!*****************************************************************************80
!
!! polygon_integral_xx() integrates the function X*X over a polygon.
!
!  Discussion:
!
!    The polygon is bounded by the points (X(1:N), Y(1:N)).
!
!    INTEGRAL = (1/12) * sum ( 1 <= I <= N )
!      ( X(I)^3 + X(I)^2 * X(I-1) + X(I) * X(I-1)^2 + X(I-1)^3 )
!      * ( Y(I) - Y(I-1) )
!
!    where X(0) and Y(0) should be replaced by X(N) and Y(N).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    SF Bockman,
!    Generalizing the Formula for Areas of Polygons to Moments,
!    American Mathematical Society Monthly,
!    1989, pages 131-132.
!
!  Input:
!
!    integer N, the number of vertices of the polygon.
!    N should be at least 3 for a nonzero result.
!
!    real ( kind = rk ) V(2,N), the coordinates of the vertices
!    of the polygon.  These vertices should be given in
!    counter clockwise order.
!
!  Output:
!
!    real ( kind = rk ) RESULT, the value of the integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer im1
  real ( kind = rk ) result
  real ( kind = rk ) v(2,n)

  result = 0.0D+00

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'polygon_integral_xx(): Fatal error!'
    write ( *, '(a)' ) '  The number of vertices must be at least 3.'
    write ( *, '(a,i8)' ) '  The input value of N = ', n
    stop 1
  end if

  do i = 1, n

    if ( i == 1 ) then
      im1 = n
    else
      im1 = i - 1
    end if

    result = result + ( v(1,i) ** 3 + v(1,i) ** 2 * v(1,im1) &
      + v(1,i) * v(1,im1) ** 2 + v(1,im1) ** 3 ) * ( v(2,i) - v(2,im1) )

  end do

  result = result / 12.0D+00

  return
end
subroutine polygon_integral_xy ( n, v, result )

!*****************************************************************************80
!
!! polygon_integral_xy() integrates the function X*Y over a polygon.
!
!  Discussion:
!
!    The polygon is bounded by the points (X(1:N), Y(1:N)).
!
!    INTEGRAL = (1/24) * sum ( 1 <= I <= N )
!      ( Y(I)   * ( 3 * X(I)^2 + 2 * X(I) * X(I-1) +     X(I-1)^2 )
!      + Y(I-1) * (     X(I)^2 + 2 * X(I) * X(I-1) + 3 * X(I-1)^2 ) )
!      * ( Y(I) - Y(I-1) )
!
!    where X(0) and Y(0) should be replaced by X(N) and Y(N).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    SF Bockman,
!    Generalizing the Formula for Areas of Polygons to Moments,
!    American Mathematical Society Monthly,
!    1989, pages 131-132.
!
!  Input:
!
!    integer N, the number of vertices of the polygon.
!    N should be at least 3 for a nonzero result.
!
!    real ( kind = rk ) V(2,N), the coordinates of the vertices
!    of the polygon.  These vertices should be given in
!    counter clockwise order.
!
!  Output:
!
!    real ( kind = rk ) RESULT, the value of the integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer im1
  real ( kind = rk ) result
  real ( kind = rk ) v(2,n)

  result = 0.0D+00

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'polygon_integral_xy(): Fatal error!'
    write ( *, '(a)' ) '  The number of vertices must be at least 3.'
    write ( *, '(a,i8)' ) '  The input value of N = ', n
    stop 1
  end if

  do i = 1, n

    if ( i == 1 ) then
      im1 = n
    else
      im1 = i - 1
    end if

    result = result + ( &
      v(2,i) * ( 3.0D+00 * v(1,i)**2 + 2.0D+00 * v(1,i) * v(1,im1) &
      + v(1,im1)**2 ) + v(2,im1) * ( v(1,i)**2 + 2.0D+00 * v(1,i) * v(1,im1) &
      + 3.0D+00 * v(1,im1)**2 ) ) * ( v(2,i) - v(2,im1) )

  end do

  result = result / 24.0D+00

  return
end
subroutine polygon_integral_y ( n, v, result )

!*****************************************************************************80
!
!! polygon_integral_y() integrates the function Y over a polygon.
!
!  Discussion:
!
!    The polygon is bounded by the points (X(1:N), Y(1:N)).
!
!    INTEGRAL = (1/6) * sum ( 1 <= I <= N )
!      - ( Y(I)^2 + Y(I) * Y(I-1) + Y(I-1)^2 ) * ( X(I) - X(I-1) )
!
!    where X(0) and Y(0) should be replaced by X(N) and Y(N).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    SF Bockman,
!    Generalizing the Formula for Areas of Polygons to Moments,
!    American Mathematical Society Monthly,
!    1989, pages 131-132.
!
!  Input:
!
!    integer N, the number of vertices of the polygon.
!    N should be at least 3 for a nonzero result.
!
!    real ( kind = rk ) V(2,N), the coordinates of the vertices
!    of the polygon.  These vertices should be given in
!    counter clockwise order.
!
!  Output:
!
!    real ( kind = rk ) RESULT, the value of the integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer im1
  real ( kind = rk ) result
  real ( kind = rk ) v(2,n)

  result = 0.0D+00

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'polygon_integral_y(): Fatal error!'
    write ( *, '(a)' ) '  The number of vertices must be at least 3.'
    write ( *, '(a,i8)' ) '  The input value of N = ', n
    stop 1
  end if

  do i = 1, n

    if ( i == 1 ) then
      im1 = n
    else
      im1 = i - 1
    end if

    result = result - ( v(2,i)**2 + v(2,i) * v(2,im1) + v(2,im1)**2 ) &
      * ( v(1,i) - v(1,im1) )

  end do

  result = result / 6.0D+00

  return
end
subroutine polygon_integral_yy ( n, v, result )

!*****************************************************************************80
!
!! polygon_integral_yy() integrates the function Y*Y over a polygon.
!
!  Discussion:
!
!    The polygon is bounded by the points (X(1:N), Y(1:N)).
!
!    INTEGRAL = (1/12) * sum ( 1 <= I <= N )
!      - ( Y(I)^3 + Y(I)^2 * Y(I-1) + Y(I) * Y(I-1)^2 + Y(I-1)^3 )
!      * ( X(I) - X(I-1) )
!
!    where X(0) and Y(0) should be replaced by X(N) and Y(N).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    SF Bockman,
!    Generalizing the Formula for Areas of Polygons to Moments,
!    American Mathematical Society Monthly,
!    1989, pages 131-132.
!
!  Input:
!
!    integer N, the number of vertices of the polygon.
!    N should be at least 3 for a nonzero result.
!
!    real ( kind = rk ) V(2,N), the coordinates of the vertices
!    of the polygon.  These vertices should be given in
!    counter clockwise order.
!
!  Output:
!
!    real ( kind = rk ) RESULT, the value of the integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer im1
  real ( kind = rk ) result
  real ( kind = rk ) v(2,n)

  result = 0.0D+00

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'polygon_integral_yy(): Fatal error!'
    write ( *, '(a)' ) '  The number of polygonal vertices must be '
    write ( *, '(a,i8)' ) '  at least 3, but the input polygon has N = ', n
    stop 1
  end if

  do i = 1, n

    if ( i == 1 ) then
      im1 = n
    else
      im1 = i - 1
    end if

    result = result - ( v(2,i)**3 + v(2,i)**2 * v(2,im1) &
      + v(2,i) * v(2,im1)**2 + v(2,im1)**3 ) * ( v(1,i) - v(1,im1) )

  end do

  result = result / 12.0D+00

  return
end
subroutine polygon_is_ccw ( n, v, value )

!*****************************************************************************80
!
!! polygon_is_ccw() determines whether a polygon is counter clockwise.
!
!  Discussion:
!
!    Notice that the polygon is dimensioned V(N,2), rather than V(2,N).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 December 2022
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer n: the number of vertices in the polygon.
!
!    real V(N,2), the coordinates of the vertices of the polygon.  
!
!  Output:
!
!    logical VALUE: true if the polygon vertices are listed in CCW order.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) angle
  real ( kind = rk ) cross
  real ( kind = rk ) dot
  integer i
  integer j
  integer k
  real ( kind = rk ) v(n,2)
  logical value

  value = .true.

  do i = 1, n

    j = mod ( i, n ) + 1
    k = mod ( i + 1, n ) + 1

    dot =   ( v(k,1) - v(j,1) ) * ( v(i,1) - v(j,1) ) &
          + ( v(k,2) - v(j,2) ) * ( v(i,2) - v(j,2) )

    cross =   ( v(k,1) - v(j,1) ) * ( v(i,2) - v(j,2) ) &
            - ( v(i,1) - v(j,1) ) * ( v(k,2) - v(j,2) )

    angle = atan2 ( cross, dot )

    if ( angle < 0.0D+00 ) then
      value = .false.
      return
    end if

  end do

  return
end
function polygon_is_convex ( n, v )

!*****************************************************************************80
!
!! polygon_is_convex() determines whether a polygon is convex.
!
!  Discussion:
!
!    If the polygon has less than 3 distinct vertices, it is
!    classified as convex degenerate.
!
!    If the polygon "goes around" more than once, it is classified
!    as NOT convex.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Peter Schorn, Frederick Fisher,
!    Testing the Convexity of a Polygon,
!    in Graphics Gems IV, 
!    edited by Paul Heckbert,
!    AP Professional, 1994,
!    T385.G6974.
!
!  Parameters
!
!    integer N, the number of vertices.
!
!    real ( kind = rk ) V(2,N), the coordinates of the vertices 
!    of the polygon.
!
!  Output:
!
!    integer POLYGON_IS_CONVEX:
!    -1, the polygon is not convex;
!     0, the polygon has less than 3 vertices; it is "degenerately" convex;
!     1, the polygon is convex and counter clockwise;
!     2, the polygon is convex and clockwise.
!
!    real ( kind = rk ) V(2,N), duplicate consecutive points have been 
!    deleted, and the vertices have been reordered so that the 
!    lexicographically least point comes first.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ), parameter :: RAD_TO_DEG = 180.0D+00 / r8_pi

  integer n

  real ( kind = rk ) angle
  integer, parameter :: CONVEX_CCW = 1
  integer, parameter :: CONVEX_CW = 2
  real ( kind = rk ) cross
  integer, parameter :: DEGENERATE_CONVEX = 0
  real ( kind = rk ) dot
  real ( kind = rk ) exterior_total
  integer i
  integer ip1
  integer ip2
  integer, parameter :: NOT_CONVEX = -1
  integer polygon_is_convex
  real ( kind = rk ) sense
  real ( kind = rk ), parameter :: tol = 1.0D+00
  real ( kind = rk ) v(2,n)

  exterior_total = 0.0D+00
!
!  If there are not at least 3 distinct vertices, we are done.
!
  if ( n < 3 ) then
    polygon_is_convex = DEGENERATE_CONVEX
    return
  end if

  sense = 0.0D+00
!
!  Consider each polygonal vertex I.
!
  do i = 1, n

    ip1 = i + 1
    if ( n < ip1 ) then
      ip1 = ip1 - n
    end if

    ip2 = i + 2
    if ( n < ip2 ) then
      ip2 = ip2 - n
    end if

    dot =   ( v(1,ip2) - v(1,ip1) ) * ( v(1,i) - v(1,ip1) ) &
          + ( v(2,ip2) - v(2,ip1) ) * ( v(2,i) - v(2,ip1) )

    cross =   ( v(1,ip2) - v(1,ip1) ) * ( v(2,i) - v(2,ip1) ) &
            - ( v(1,i)   - v(1,ip1) ) * ( v(2,ip2) - v(2,ip1) )

    angle = atan2 ( cross, dot )
!
!  See if the turn defined by this vertex is our first indication of
!  the "sense" of the polygon, or if it disagrees with the previously
!  defined sense.
!
    if ( sense == 0.0D+00 ) then

      if ( angle < 0.0D+00 ) then
        sense = -1.0D+00
      else if ( 0.0D+00 < angle ) then
        sense = +1.0D+00
      end if

    else if ( sense == 1.0D+00 ) then

      if ( angle < 0.0D+00 ) then
        polygon_is_convex = NOT_CONVEX
        return
      end if

    else if ( sense == -1.0D+00 ) then

      if ( 0.0D+00 < angle ) then
        polygon_is_convex = NOT_CONVEX
        return
      end if

    end if
!
!  If the exterior total is greater than 360, then the polygon is
!  going around again.
!
    angle = atan2 ( -cross, -dot )

    exterior_total = exterior_total + angle

    if ( 360.0D+00 + tol < abs ( exterior_total ) * RAD_TO_DEG ) then
      polygon_is_convex = NOT_CONVEX
      return
    end if

  end do

  if ( sense == +1.0D+00 ) then
    polygon_is_convex = CONVEX_CCW
  else if ( sense == -1.0D+00 ) then
    polygon_is_convex = CONVEX_CW
  end if

  return
end
subroutine polygon_perimeter ( n, v, perimeter )

!*****************************************************************************80
!
!! polygon_perimeter() computes the perimeter of a polygon.
!
!  Discussion:
!
!    The perimeter is simply the sum of the side lengths.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 October 2015
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of vertices of the polygon.
!
!    real ( kind = rk ) V(2,N), the vertices.
!
!  Output:
!
!    real ( kind = rk ) PERIMETER, the perimeter.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer im1
  real ( kind = rk ) l
  real ( kind = rk ) perimeter
  real ( kind = rk ) v(2,n)

  perimeter = 0.0D+00

  im1 = n

  do i = 1, n
    l = sqrt ( ( v(1,im1) - v(1,i) ) ** 2 + ( v(2,im1) - v(2,i) ) ** 2 )
    perimeter = perimeter + l
    im1 = i
  end do

  return
end
subroutine polygon_perimeter_quad ( n, v, hmax, f, value )

!*****************************************************************************80
!
!! polygon_perimeter_quad() estimates an integral over the perimeter of a polygon.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 October 2015
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of vertices of the polygon.
!
!    real ( kind = rk ) V(2,N), the vertices.
!
!    real ( kind = rk ) HMAX, the maximum length of a quadrature interval.
!
!    real ( kind = rk ), external F ( X, Y ), a function whose integral 
!    over the perimeter is desired.
!
!  Output:
!
!    real ( kind = rk ) VALUE, the estimated integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) dxy
  real ( kind = rk ), external :: f
  real ( kind = rk ) fxy
  real ( kind = rk ) hmax
  integer i
  integer i4_wrap
  integer ip1
  integer j
  real ( kind = rk ) l
  integer m
  real ( kind = rk ) v(2,n)
  real ( kind = rk ) value
  real ( kind = rk ) x
  real ( kind = rk ) y

  value = 0.0D+00

  do i = 1, n

    ip1 = i4_wrap ( i + 1, 1, n )
    l = sqrt ( ( v(1,ip1) - v(1,i) ) ** 2 + ( v(2,ip1) - v(2,i) ) ** 2 )
    m = ceiling ( l / hmax )
    dxy = l / real ( m, kind = rk )

    do j = 1, 2 * m - 1, 2
      x = ( real ( 2 * m - j, kind = rk ) * v(1,i) &
          + real (         j, kind = rk ) * v(1,ip1) ) &
          / real ( 2 * m,     kind = rk )
      y = ( real ( 2 * m - j, kind = rk ) * v(2,i) &
          + real (         j, kind = rk ) * v(2,ip1) ) &
          / real ( 2 * m,     kind = rk )
      fxy = f ( x, y )
      value = value + fxy * dxy
    end do

  end do

  return
end
subroutine polygon_point_dist ( n, v, p, dist )

!*****************************************************************************80
!
!! polygon_point_dist(): distance ( polygon, point ).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    28 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of vertices.
!
!    real ( kind = rk ) V(2,N), the triangle vertices.
!
!    real ( kind = rk ) P(2), the point to be checked.
!
!  Output:
!
!    real ( kind = rk ) DIST, the distance from the point to the
!    polygon.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) dist
  real ( kind = rk ) dist2
  integer i4_wrap
  integer j
  integer jp1
  real ( kind = rk ) p(2)
  real ( kind = rk ) v(2,n)
!
!  Find the distance to each of the line segments.
!
  dist = huge ( dist )

  do j = 1, n

    jp1 = i4_wrap ( j+1, 1, n )

    call segment_point_dist ( v(1:2,j), v(1:2,jp1), p, dist2 )

    if ( dist2 < dist ) then
      dist = dist2
    end if

  end do

  return
end
subroutine polygon_point_near ( n, v, p, pn, dist )

!*****************************************************************************80
!
!! polygon_point_near() computes the nearest point on a polygon.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    28 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) V(2,N), the polygon vertices.
!
!    real ( kind = rk ) P(2), the point whose nearest polygon point
!    is to be determined.
!
!  Output:
!
!    real ( kind = rk ) PN(2), the nearest point to P.
!
!    real ( kind = rk ) DIST, the distance from the point to the
!    polygon.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) dist
  real ( kind = rk ) dist2
  integer i4_wrap
  integer j
  integer jp1
  real ( kind = rk ) p(2)
  real ( kind = rk ) pn(2)
  real ( kind = rk ) pn2(2)
  real ( kind = rk ) tval
  real ( kind = rk ) v(2,n)
!
!  Find the distance to each of the line segments that make up the edges
!  of the polygon.
!
  dist = huge ( dist )
  pn(1:2) = 0.0D+00

  do j = 1, n

    jp1 = i4_wrap ( j+1, 1, n )

    call segment_point_near ( v(1:2,j), v(1:2,jp1), p, &
      pn2, dist2, tval )

    if ( dist2 < dist ) then
      dist = dist2
      pn(1:2) = pn2(1:2)
    end if

  end do

  return
end
subroutine polygon_sample ( nv, v, n,  s )

!*****************************************************************************80
!
!! polygon_sample() uniformly samples a polygon.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 August 2005
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
  integer triangles(3,nv-2)
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
      v(1,triangles(1,i)), v(2,triangles(1,i)), &
      v(1,triangles(2,i)), v(2,triangles(2,i)), &
      v(1,triangles(3,i)), v(2,triangles(3,i)) )
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

    s(1:2,j) = ( 1.0D+00 - r(1) - r(2) ) * v(1:2,triangles(1,i)) &
                         + r(1)          * v(1:2,triangles(2,i)) &
                                + r(2)   * v(1:2,triangles(3,i))
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
!  Input:
!
!    integer N, the number of vertices.
!
!    real ( kind = rk ) X(N), Y(N), the coordinates of each vertex.
!
!  Output:
!
!    integer TRIANGLES(3,N-2), the triangles of the 
!    triangulation.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) area
  logical diagonal
  logical ear(n)
  integer i
  integer i0
  integer i1
  integer i2
  integer i3
  integer i4
  integer next(n)
  integer node
  integer node_m1
  integer prev(n)
  integer triangle_num
  integer triangles(3,n-2)
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)
!
!  We must have at least 3 vertices.
!
  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'polygon_triangulate(): Fatal error!'
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
      write ( *, '(a)' ) 'polygon_triangulate(): Fatal error!'
      write ( *, '(a)' ) '  Two consecutive nodes are identical.'
      stop 1
    end if
    node_m1 = node
  end do
!
!  Area must be positive.
!
  area = 0.0D+00
  do node = 1, n - 2
    area = area + 0.5D+00 * &
    ( &
        ( x(node+1) - x(node) ) * ( y(node+2) - y(node) ) &
      - ( x(node+2) - x(node) ) * ( y(node+1) - y(node) ) &
    )
  end do

  if ( area <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'polygon_triangulate(): Fatal error!'
    write ( *, '(a)' ) '  Polygon has zero or negative area.'
    stop 1
  end if
!
!  PREV and NEXT point to the previous and next nodes.
!
  i = 1
  prev(i) = n
  next(i) = i + 1

  do i = 2, n - 1
    prev(i) = i - 1
    next(i) = i + 1
  end do

  i = n
  prev(i) = i - 1
  next(i) = 1
!
!  EAR indicates whether the node and its immediate neighbors form an ear
!  that can be sliced off immediately.
!
  do i = 1, n
    ear(i) = diagonal ( prev(i), next(i), n, prev, next, x, y )
  end do

  triangle_num = 0

  i2 = 1

  do while ( triangle_num < n - 3 )
!
!  If I2 is an ear, gather information necessary to carry out
!  the slicing operation and subsequent "healing".
!
    if ( ear(i2) ) then

      i3 = next(i2)
      i4 = next(i3)
      i1 = prev(i2)
      i0 = prev(i1)
!
!  Make vertex I2 disappear.
!
      next(i1) = i3
      prev(i3) = i1
!
!  Update the earity of I1 and I3, because I2 disappeared.
!
      ear(i1) = diagonal ( i0, i3, n, prev, next, x, y )
      ear(i3) = diagonal ( i1, i4, n, prev, next, x, y )
!
!  Add the diagonal [I3, I1, I2] to the list.
!
      triangle_num = triangle_num + 1
      triangles(1,triangle_num) = i3
      triangles(2,triangle_num) = i1
      triangles(3,triangle_num) = i2

    end if
!
!  Try the next vertex.
!
    i2 = next(i2)

  end do
!
!  The last triangle is formed from the three remaining vertices.
!
  i3 = next(i2)
  i1 = prev(i2)

  triangle_num = triangle_num + 1
  triangles(1,triangle_num) = i3
  triangles(2,triangle_num) = i1
  triangles(3,triangle_num) = i2

  return
end
function r8_degrees ( radians )

!*****************************************************************************80
!
!! r8_degrees() converts an angle from radian to degree measure.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) RADIANS, the angle measurement in radians.
!
!  Output:
!
!    real ( kind = rk ) R8_DEGREES, the angle measurement in degrees.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r8_degrees
  real ( kind = rk ), parameter :: r8_pi = 3.1415926535897932384626434D+00
  real ( kind = rk ) radians

  r8_degrees = radians * 180.0D+00 / r8_pi

  return
end
subroutine r8_fake_use ( x )

!*****************************************************************************80
!
!! r8_fake_use() pretends to use an R8 variable.
!
!  Discussion:
!
!    Some compilers will issue a warning if a variable is unused.
!    Sometimes there's a good reason to include a variable in a program,
!    but not to use it.  Calling this function with that variable as
!    the argument will shut the compiler up.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 April 2020
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) X, the variable to be "used".
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) x

  if ( x /= x ) then
    write ( *, '(a)' ) '  r8_fake_use: variable is NAN.'
  end if

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! r8mat_print() prints a real matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, the number of rows in A.
!
!    integer N, the number of columns in A.
!
!    real ( kind = rk ) A(M,N), the matrix.
!
!    character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  character ( len = * ) title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! r8mat_print_some prints some of a real matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, N, the number of rows and columns.
!
!    real ( kind = rk ) A(M,N), an M by N matrix to be printed.
!
!    integer ILO, JLO, the first row and column to print.
!
!    integer IHI, JHI, the last row and column to print.
!
!    character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: incx = 5
  integer m
  integer n

  real ( kind = rk ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer i
  integer i2hi
  integer i2lo
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)'
    return
  end if

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        write ( ctemp(j2), '(g14.6)' ) a(i,j)

      end do

      write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine r8mat_solve ( n, rhs_num, a, info )

!*****************************************************************************80
!
!! r8mat_solve() uses Gauss-Jordan elimination to solve an N by N linear system.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the order of the matrix.
!
!    integer RHS_NUM, the number of right hand sides.
!    RHS_NUM must be at least 0.
!
!    real ( kind = rk ) A(N,N+RHS_NUM), contains in rows and
!    columns 1 to N the coefficient matrix, and in columns N+1 through
!    N+RHS_NUM, the right hand sides.
!
!  Output:
!
!    integer INFO, singularity flag.
!    0, the matrix was not singular, the solutions were computed;
!    J, factorization failed on step J, and the solutions could not
!    be computed.
!
!    real ( kind = rk ) A(N,N+RHS_NUM), the coefficient matrix
!    area has been destroyed, while the right hand sides have
!    been overwritten with the corresponding solutions.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer rhs_num

  real ( kind = rk ) a(n,n+rhs_num)
  real ( kind = rk ) apivot
  real ( kind = rk ) factor
  integer i
  integer info
  integer ipivot
  integer j
  real ( kind = rk ) t(n+rhs_num)

  info = 0

  do j = 1, n
!
!  Choose a pivot row.
!
    ipivot = j
    apivot = a(j,j)

    do i = j + 1, n
      if ( abs ( apivot ) < abs ( a(i,j) ) ) then
        apivot = a(i,j)
        ipivot = i
      end if
    end do

    if ( apivot == 0.0D+00 ) then
      info = j
      return
    end if
!
!  The pivot row moves into the J-th row.
!
    if ( ipivot /= j ) then
      t(       1:n+rhs_num) = a(ipivot,1:n+rhs_num)
      a(ipivot,1:n+rhs_num) = a(j,     1:n+rhs_num)
      a(j,     1:n+rhs_num) = t(       1:n+rhs_num)
    end if
!
!  A(J,J) becomes 1.
!
    a(j,j) = 1.0D+00
    a(j,j+1:n+rhs_num) = a(j,j+1:n+rhs_num) / apivot
!
!  A(I,J) becomes 0.
!
    do i = 1, n

      if ( i /= j ) then
        factor = a(i,j)
        a(i,j) = 0.0D+00
        a(i,j+1:n+rhs_num) = a(i,j+1:n+rhs_num) - factor * a(j,j+1:n+rhs_num)
      end if

    end do

  end do

  return
end
subroutine r8mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! r8mat_transpose_print() prints an R8MAT, transposed.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, N, the number of rows and columns.
!
!    real ( kind = rk ) A(M,N), an M by N matrix to be printed.
!
!    character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  character ( len = * ) title

  call r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! r8mat_transpose_print_some() prints some of an R8MAT, transposed.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, N, the number of rows and columns.
!
!    real ( kind = rk ) A(M,N), an M by N matrix to be printed.
!
!    integer ILO, JLO, the first row and column to print.
!
!    integer IHI, JHI, the last row and column to print.
!
!    character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: incx = 5
  integer m
  integer n

  real ( kind = rk ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer i
  integer i2
  integer i2hi
  integer i2lo
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)'
    return
  end if

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(i8,6x)' ) i
    end do

    write ( *, '(''  Row   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc
        i = i2lo - 1 + i2
        write ( ctemp(i2), '(g14.6)' ) a(i,j)
      end do

      write ( *, '(i5,a,5a14)' ) j, ':', ( ctemp(i), i = 1, inc )

    end do

  end do

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! r8vec_print() prints an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of components of the vector.
!
!    real ( kind = rk ) A(N), the vector to be printed.
!
!    character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  integer i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
  end do

  return
end
subroutine segment_point_dist ( p1, p2, p, dist )

!*****************************************************************************80
!
!! segment_point_dist(): distance ( line segment, point ).
!
!  Discussion:
!
!    A line segment is the finite portion of a line that lies between
!    two points P1 and P2.
!
!    The nearest point will satisfy the condition
!
!      PN = (1-T) * P1 + T * P2.
!
!    T will always be between 0 and 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 May 2006
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) P1(2), P2(2), the endpoints of the line segment.
!
!    real ( kind = rk ) P(2), the point whose nearest neighbor on the line
!    segment is to be determined.
!
!  Output:
!
!    real ( kind = rk ) DIST, the distance from the point to the
!    line segment.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) bot
  real ( kind = rk ) dist
  real ( kind = rk ) p(2)
  real ( kind = rk ) p1(2)
  real ( kind = rk ) p2(2)
  real ( kind = rk ) pn(2)
  real ( kind = rk ) t
!
!  If the line segment is actually a point, then the answer is easy.
!
  if ( all ( p1(1:2) == p2(1:2) ) ) then

    t = 0.0D+00

  else

    bot = sum ( ( p2(1:2) - p1(1:2) )**2 )

    t = sum ( ( p(1:2)  - p1(1:2) ) &
            * ( p2(1:2) - p1(1:2) ) ) / bot

    t = max ( t, 0.0D+00 )
    t = min ( t, 1.0D+00 )

  end if

  pn(1:2) = p1(1:2) + t * ( p2(1:2) - p1(1:2) )

  dist = sqrt ( sum ( ( p(1:2) - pn(1:2) )**2 ) )

  return
end
subroutine segment_point_near ( p1, p2, p, pn, dist, t )

!*****************************************************************************80
!
!! segment_point_near(): nearest point on line segment to point.
!
!  Discussion:
!
!    A line segment is the finite portion of a line that lies between
!    two points P1 and P2.
!
!    The nearest point will satisfy the condition
!
!      PN = (1-T) * P1 + T * P2.
!
!    T will always be between 0 and 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 May 2006
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) P1(2), P2(2), the endpoints of the line segment.
!
!    real ( kind = rk ) P(2), the point whose nearest neighbor
!    on the line segment is to be determined.
!
!  Output:
!
!    real ( kind = rk ) PN(2), the point on the line segment which is
!    nearest the point P.
!
!    real ( kind = rk ) DIST, the distance from the point to the 
!    nearest point on the line segment.
!
!    real ( kind = rk ) T, the relative position of the point PN
!    to the points P1 and P2.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) bot
  real ( kind = rk ) dist
  real ( kind = rk ) p(2)
  real ( kind = rk ) p1(2)
  real ( kind = rk ) p2(2)
  real ( kind = rk ) pn(2)
  real ( kind = rk ) t
!
!  If the line segment is actually a point, then the answer is easy.
!
  if ( all ( p1(1:2) == p2(1:2) ) ) then

    t = 0.0D+00

  else

    bot = sum ( ( p2(1:2) - p1(1:2) )**2 )

    t = sum ( ( p(1:2)  - p1(1:2) ) &
            * ( p2(1:2) - p1(1:2) ) ) / bot

    t = max ( t, 0.0D+00 )
    t = min ( t, 1.0D+00 )

  end if

  pn(1:2) = p1(1:2) + t * ( p2(1:2) - p1(1:2) )

  dist = sqrt ( sum ( ( p(1:2) - pn(1:2) )**2 ) )

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
!  Input:
!
!    real ( kind = rk ) XA, YA, XB, YB, XC, YC, the coordinates of
!    the vertices of the triangle, given in counterclockwise order.
!
!  Output:
!
!    real ( kind = rk ) TRIANGLE_AREA, the signed area of the triangle.
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
subroutine triangle_barycentric ( t, p, xsi )

!*****************************************************************************80
!
!! triangle_barycentric() finds the barycentric coordinates of a point.
!
!  Discussion:
!
!    The barycentric coordinate of point P related to vertex A can be
!    interpreted as the ratio of the area of the triangle with 
!    vertex A replaced by vertex P to the area of the original 
!    triangle.
!
!    This routine assumes that the triangle vertices are given in
!    counter clockwise order.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    28 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) T(2,3), the triangle vertices.
!    The vertices should be given in counter clockwise order.
!
!    real ( kind = rk ) P(2), the point to be checked.
!
!  Output:
!
!    real ( kind = rk ) XSI(3), the barycentric coordinates of P
!    with respect to the triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(2,3)
  integer info
  real ( kind = rk ) p(2)
  real ( kind = rk ) t(2,3)
  real ( kind = rk ) xsi(3)
!
!  Set up the linear system
!
!    ( X2-X1  X3-X1 ) XSI(1)  = X-X1
!    ( Y2-Y1  Y3-Y1 ) XSI(2)    Y-Y1
!
!  which is satisfied by the barycentric coordinates of P.
!
  a(1,1) = t(1,2) - t(1,1)
  a(1,2) = t(1,3) - t(1,1)
  a(1,3) = p(1)   - t(1,1)

  a(2,1) = t(2,2) - t(2,1)
  a(2,2) = t(2,3) - t(2,1)
  a(2,3) = p(2)   - t(2,1)
!
!  Solve the linear system.
!
  call r8mat_solve ( 2, 1, a, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'triangle_barycentric(): Fatal error!'
    write ( *, '(a)' ) '  The linear system is singular.'
    write ( *, '(a)' ) '  The input data does not form a proper triangle.'
    stop 1
  end if

  xsi(1) = a(1,3)
  xsi(2) = a(2,3)
  xsi(3) = 1.0D+00 - xsi(1) - xsi(2)

  return
end
subroutine triangle_contains_point_1 ( t, p, inside )

!*****************************************************************************80
!
!! triangle_contains_point_1() finds if a point is inside a triangle.
!
!  Discussion:
!
!    It is conventional to list the triangle vertices in counter clockwise
!    order.  However, this routine does not require a particular order
!    for the vertices.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) T(2,3), the triangle vertices.
!
!    real ( kind = rk ) P(2), the point to be checked.
!
!  Output:
!
!    logical INSIDE, is TRUE if the point is inside the triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  logical inside
  real ( kind = rk ) p(2)
  real ( kind = rk ) t(2,3)
  real ( kind = rk ) xsi(3)

  call triangle_barycentric ( t, p, xsi )

  if ( any ( xsi(1:3) < 0.0D+00 ) ) then
    inside = .false.
  else
    inside = .true.
  end if

  return
end
 
