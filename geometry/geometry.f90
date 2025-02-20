subroutine angle_box_2d ( dist, p1, p2, p3, p4, p5 )

!*****************************************************************************80
!
!! angle_box_2d() "boxes" an angle defined by three points in 2D.
!
!  Discussion:
!
!    The routine is given points P1, P2 and P3, determining the two lines:
!      P1 to P2
!    and
!      P2 to P3
!    and a nonnegative distance
!      DIST.
!
!    The routine returns a pair of "corner" points
!      P4 and P5
!    both of which are a distance DIST from both lines, and in fact,
!    both of which are a distance DIST from P2.
!
!                         /  P3
!                        /   /   /
!     - - - - - - - - -P4 - / -P6 - - -
!                      /   /   /
!    P1---------------/--P2-----------------
!                    /   /   /
!     - - - - - - -P7 - / -P5 - - - - -
!                  /   /   /
!
!    In the illustration, P1, P2 and P3 are the points defining the lines.
!
!    P4 and P5 represent the desired "corner points", which
!    are on the positive or negative sides of both lines.
!
!    P6 and P7 represent the undesired points, which
!    are on the positive side of one line and the negative of the other.
!
!    Special cases:
!
!    if P1 = P2, this is the same as extending the line from
!    P3 through P2 without a bend.
!
!    if P3 = P2, this is the same as extending the line from
!    P1 through P2 without a bend.
!
!    if P1 = P2 = P3 this is an error.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) DIST, the nonnegative distance from P1
!    to the computed points P4 and P5.
!
!    real ( kind = rk ) P1(2), P2(2), P3(2).
!    P1 and P2 are distinct points that define a line.
!    P2 and P3 are distinct points that define a line.
!
!  Output:
!
!    real ( kind = rk ) P4(2), P5(2), points which lie DIST units from
!    the line between P1 and P2, and from the line between P2 and P3.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) dist
  real ( kind = rk ) p1(2)
  real ( kind = rk ) p2(2)
  real ( kind = rk ) p3(2)
  real ( kind = rk ) p4(2)
  real ( kind = rk ) p5(2)
  real ( kind = rk ) stheta
  real ( kind = rk ) temp1
  real ( kind = rk ) temp2
  real ( kind = rk ) u(2)
  real ( kind = rk ) u1(2)
  real ( kind = rk ) u2(2)
!
!  If DIST = 0, assume the user knows best.
!
  if ( dist == 0.0D+00 ) then
    p4(1:2) = p2(1:2)
    p5(1:2) = p2(1:2)
    return
  end if
!
!  Fail if all three points are equal.
!
  if ( all ( p1(1:2) == p2(1:2) ) .and. &
       all ( p2(1:2) == p3(1:2) ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ANGLE_BOX_2D - Fatal error!'
    write ( *, '(a)' ) '  Input points P1 = P2 = P3.'
    write ( *, '(a,2g14.6)' ) '  P1 = ', p1(1:2)
    stop 1
  end if
!
!  If P1 = P2, extend the line through the doubled point.
!
  if ( all ( p1(1:2) == p2(1:2) ) ) then
    u2(1) = p3(2) - p2(2)
    u2(2) = p2(1) - p3(1)
    temp1 = sqrt ( sum ( u2(1:2)**2 ) )
    u2(1:2) = u2(1:2) / temp1
    p4(1:2) = p2(1:2) + dist * u2(1:2)
    p5(1:2) = p2(1:2) - dist * u2(1:2)
    return
  end if
!
!  If P2 = P3, extend the line through the doubled point.
!
  if ( all ( p2(1:2) == p3(1:2) ) ) then
    u1(1) = p1(2) - p2(2)
    u1(2) = p2(1) - p1(1)
    temp1 = sqrt ( sum ( u1(1:2)**2 ) )
    u1(1:2) = u1(1:2) / temp1
    p4(1:2) = p2(1:2) + dist * u1(1:2)
    p5(1:2) = p2(1:2) - dist * u1(1:2)
    return
  end if
!
!  Compute the unit normal vectors to each line.
!  We choose the sign so that the unit normal to line 1 has
!  a positive dot product with line 2.
!
  u1(1) = p1(2) - p2(2)
  u1(2) = p2(1) - p1(1)
  temp1 = sqrt ( sum ( u1(1:2)**2 ) )
  u1(1:2) = u1(1:2) / temp1

  temp1 = dot_product ( u1(1:2), p3(1:2) - p2(1:2) )

  if ( temp1 < 0.0D+00 ) then
    u1(1:2) = -u1(1:2)
  end if

  u2(1) = p3(2) - p2(2)
  u2(2) = p2(1) - p3(1)
  temp1 = sqrt ( sum ( u2(1:2)**2 ) )
  u2(1:2) = u2(1:2) / temp1

  temp1 = dot_product ( u2(1:2), p1(1:2) - p2(1:2) )
  if ( temp1 < 0.0D+00 ) then
    u2(1:2) = -u2(1:2)
  end if
!
!  Try to catch the case where we can't determine the
!  sign of U1, because both U1 and -U1 are perpendicular
!  to (P3-P2)...and similarly for U2 and (P1-P2).
!
  temp1 = dot_product ( u1(1:2), p3(1:2) - p2(1:2) )
  temp2 = dot_product ( u2(1:2), p1(1:2) - p2(1:2) )

  if ( temp1 == 0.0D+00 .or. temp2 == 0.0D+00 ) then

    if ( dot_product ( u1(1:2), u2(1:2) ) < 0.0D+00 ) then
      u1(1:2) = -u1(1:2)
    end if

  end if
!
!  Try to catch a line turning back on itself, evidenced by
!    Cos(theta) = (P3-P2) dot (P2-P1) / ( norm(P3-P2) * norm(P2-P1) )
!  being -1, or very close to -1.
!
  temp1 = dot_product ( p3(1:2) - p2(1:2), p2(1:2) - p1(1:2) ) 

  temp1 = temp1 / &
        ( sqrt ( sum ( ( p3(1:2) - p2(1:2) )**2 ) ) &
        * sqrt ( sum ( ( p2(1:2) - p1(1:2) )**2 ) ) )

  if ( temp1 < -0.99D+00 ) then
    temp1 = sqrt ( sum ( ( p2(1:2) - p1(1:2) )**2 ) )
    p4(1:2) = p2(1:2) + dist * ( p2(1:2) - p1(1:2) ) &
      / temp1 + dist * u1(1:2)
    p5(1:2) = p2(1:2) + dist * ( p2(1:2) - p1(1:2) ) &
      / temp1 - dist * u1(1:2)
    return
  end if
!
!  Compute the "average" unit normal vector.
!
!  The average of the unit normals could be zero, but only when
!  the second line has the same direction and opposite sense
!  of the first, and we've already checked for that case.
!
!  Well, check again!  This problem "bit" me in the case where
!  P1 = P2, which I now treat specially just to guarantee I
!  avoid this problem!
!
  if ( dot_product ( u1(1:2), u2(1:2) ) < 0.0D+00 ) then
    u2(1:2) = -u2(1:2)
  end if

  u(1:2) = 0.5D+00 * ( u1(1:2) + u2(1:2) )
  temp1 = sqrt ( sum ( u(1:2)**2 ) )
  u(1:2) = u(1:2) / temp1
!
!  You must go DIST/STHETA units along this unit normal to
!  result in a distance DIST from line1 (and line2).
!
  stheta = dot_product ( u(1:2), u1(1:2) )

  p4(1:2) = p2(1:2) + dist * u(1:2) / stheta
  p5(1:2) = p2(1:2) - dist * u(1:2) / stheta

  return
end
subroutine angle_contains_point_2d ( p1, p2, p3, p, inside )

!*****************************************************************************80
!
!! angle_contains_point_2d determines if an angle contains a point, in 2D.
!
!  Discussion:
!
!    The angle is defined by the sequence of points P1, P2 and P3.
!
!    The point is "contained" by the angle if the ray P - P2
!    is between (in a counter clockwise sense) the rays P1 - P2
!    and P3 - P2.
!
!        P1
!        /
!       /   P
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
!    15 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) P1(2), P2(2), P3(2), the coordinates of
!    three points that define the angle.  The order of these points matters!
!
!    real ( kind = rk ) P(2), the point to be checked.
!
!  Output:
!
!    logical INSIDE, is TRUE if the point is inside the angle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) angle_rad_2d
  logical inside
  real ( kind = rk ) p(2)
  real ( kind = rk ) p1(2)
  real ( kind = rk ) p2(2)
  real ( kind = rk ) p3(2)

  if ( angle_rad_2d ( p1, p2, p ) <= angle_rad_2d ( p1, p2, p3 ) ) then
    inside = .true.
  else
    inside = .false.
  end if

  return
end
function angle_degrees_2d ( p1, p2, p3 )

!*****************************************************************************80
!
!! angle_degrees_2d returns the angle swept out between two rays in 2D.
!
!  Discussion:
!
!    Except for the zero angle case, it should be true that
!
!      angle_degrees_2d ( P1, P2, P3 ) + angle_degrees_2d ( P3, P2, P1 ) = 360.0
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
!    14 January 2005
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
!    real ( kind = rk ) angle_degrees_2d, the angle swept out by the 
!    rays, measured in degrees.  0 <= angle_degrees_2d < 360.  If either ray 
!    has zero length, then angle_degrees_2d is set to 0.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) angle_degrees_2d
  real ( kind = rk ) angle_rad_2d
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) radians_to_degrees
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) p3(dim_num)

  p(1) = ( p3(1) - p2(1) ) * ( p1(1) - p2(1) ) &
       + ( p3(2) - p2(2) ) * ( p1(2) - p2(2) )

  p(2) = ( p3(1) - p2(1) ) * ( p1(2) - p2(2) ) &
       - ( p3(2) - p2(2) ) * ( p1(1) - p2(1) )

  if ( p(1) == 0.0D+00 .and. p(2) == 0.0D+00 ) then
    angle_degrees_2d = 0.0D+00
    return
  end if

  angle_rad_2d = atan2 ( p(2), p(1) )

  if ( angle_rad_2d < 0.0D+00 ) then
    angle_rad_2d = angle_rad_2d + 2.0D+00 * r8_pi
  end if

  angle_degrees_2d = radians_to_degrees ( angle_rad_2d )

  return
end
subroutine angle_half_2d ( p1, p2, p3, p4 )

!*****************************************************************************80
!
!! ANGLE_HALF_2D finds half an angle in 2D.
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
!  Output:
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
function angle_rad_2d ( p1, p2, p3 )

!*****************************************************************************80
!
!! ANGLE_RAD_2D returns the angle in radians swept out between two rays in 2D.
!
!  Discussion:
!
!    Except for the zero angle case, it should be true that
!
!      ANGLE_RAD_2D ( P1, P2, P3 ) + ANGLE_RAD_2D ( P3, P2, P1 ) = 2 * PI
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
!    real ( kind = rk ) ANGLE_RAD_2D, the angle swept out by the rays,
!    in radians.  0 <= ANGLE_RAD_2D < 2 * PI.  If either ray has zero
!    length, then ANGLE_RAD_2D is set to 0.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) angle_rad_2d
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
    angle_rad_2d = 0.0D+00
    return
  end if

  angle_rad_2d = atan2 ( p(2), p(1) )

  if ( angle_rad_2d < 0.0D+00 ) then
    angle_rad_2d = angle_rad_2d + 2.0D+00 * r8_pi
  end if

  return
end
function angle_rad_3d ( p1, p2, p3 )

!*****************************************************************************80
!
!! ANGLE_RAD_3D returns the angle in radians between two rays in 3D.
!
!  Discussion:
!
!    The routine always computes the SMALLER of the two angles between
!    two rays.  Thus, if the rays make an (exterior) angle of
!    1.5 pi radians, the (interior) angle of 0.5 pi radians will be reported.
!
!    X dot Y = Norm(X) * Norm(Y) * Cos ( Angle(X,Y) )
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) P1(3), P2(3), P3(3), points defining an angle.
!    The rays are P1 - P2 and P3 - P2.
!
!  Output:
!
!    real ( kind = rk ) ANGLE_RAD_3D, the angle between the two rays,
!    in radians.  This value will always be between 0 and PI.  If either ray has
!    zero length, then the angle is returned as zero.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) angle_rad_3d
  real ( kind = rk ) dot
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) p3(dim_num)
  real ( kind = rk ) r8_acos
  real ( kind = rk ) v1norm
  real ( kind = rk ) v2norm

  v1norm = sqrt ( sum ( ( p1(1:dim_num) - p2(1:dim_num) )**2 ) )

  if ( v1norm == 0.0D+00 ) then
    angle_rad_3d = 0.0D+00
    return
  end if

  v2norm = sqrt ( sum ( ( p3(1:dim_num) - p2(1:dim_num) )**2 ) )

  if ( v2norm == 0.0D+00 ) then
    angle_rad_3d = 0.0D+00
    return
  end if

  dot = sum ( ( p1(1:dim_num) - p2(1:dim_num) ) &
            * ( p3(1:dim_num) - p2(1:dim_num) ) )

  angle_rad_3d = r8_acos ( dot / ( v1norm * v2norm ) )

  return
end
function angle_rad_nd ( dim_num, v1, v2 )

!*****************************************************************************80
!
!! ANGLE_RAD_ND returns the angle in radians between two rays in ND.
!
!  Discussion:
!
!    This routine always computes the SMALLER of the two angles between
!    two rays.  Thus, if the rays make an (exterior) angle of 1.5 PI,
!    then the (interior) angle of 0.5 PI is reported.
!
!    X dot Y = Norm(X) * Norm(Y) * Cos( Angle(X,Y) )
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer DIM_NUM, the spatial dimension.
!
!    real ( kind = rk ) V1(DIM_NUM), V2(DIM_NUM), the two rays.
!
!  Output:
!
!    real ( kind = rk ) ANGLE_RAD_ND, the angle between the rays,
!    in radians.  This value will always be between 0 and PI.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num

  real ( kind = rk ) angle_rad_nd
  real ( kind = rk ) dot
  real ( kind = rk ) r8_acos
  real ( kind = rk ) v1(dim_num)
  real ( kind = rk ) v1norm
  real ( kind = rk ) v2(dim_num)
  real ( kind = rk ) v2norm

  dot = dot_product ( v1(1:dim_num), v2(1:dim_num) )

  v1norm = sqrt ( sum ( v1(1:dim_num)**2 ) )

  if ( v1norm == 0.0D+00 ) then
    angle_rad_nd = 0.0D+00
    return
  end if

  v2norm = sqrt ( sum ( v2(1:dim_num)**2 ) )

  if ( v2norm == 0.0D+00 ) then
    angle_rad_nd = 0.0D+00
    return
  end if

  angle_rad_nd = r8_acos ( dot / ( v1norm * v2norm ) )

  return
end
subroutine angle_turn_2d ( p1, p2, p3, turn )

!*****************************************************************************80
!
!! ANGLE_TURN_2D computes a turning angle in 2D.
!
!  Discussion:
!
!    This routine is most useful when considering the vertices of a 
!    polygonal shape.  We wish to distinguish between angles that "turn
!    in" to the shape, (between 0 and 180 degrees) and angles that 
!    "turn out" (between 180 and 360 degrees), as we traverse the boundary.
!
!    If we compute the interior angle and subtract 180 degrees, we get the
!    supplementary angle, which has the nice property that it is
!    negative for "in" angles and positive for "out" angles, and is zero if
!    the three points actually lie along a line.
!
!    Assuming P1, P2 and P3 define an angle, the TURN can be
!    defined to be either:
!
!    * the supplementary angle to the angle formed by P1=P2=P3, or
!
!    * the angle between the vector ( P3-P2) and the vector -(P1-P2), 
!      where -(P1-P2) can be understood as the vector that continues
!      through P2 from the direction P1.
!
!    The turning will be zero if P1, P2 and P3 lie along a straight line.
!
!    It will be a positive angle if the turn from the previous direction
!    is counter clockwise, and negative if it is clockwise.
!
!    The turn is given in radians, and will lie between -PI and PI.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    13 August 2014
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) P1(2), P2(2), P3(2), the points that form
!    the angle.
!
!  Output:
!
!    real ( kind = rk ) TURN, the turn angle, between -PI and PI.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) p3(dim_num)
  real ( kind = rk ) r8_atan
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) turn

  p(1) = ( p3(1) - p2(1) ) * ( p1(1) - p2(1) ) &
       + ( p3(2) - p2(2) ) * ( p1(2) - p2(2) )

  p(2) = ( p3(1) - p2(1) ) * ( p1(2) - p2(2) ) &
       - ( p3(2) - p2(2) ) * ( p1(1) - p2(1) )

  if ( p(1) == 0.0D+00 .and. p(2) == 0.0D+00 ) then
    turn = 0.0D+00
  else
    turn = r8_pi - r8_atan ( p(2), p(1) )
  end if

  return
end
subroutine annulus_area_2d ( r1, r2, area )

!*****************************************************************************80
!
!! ANNULUS_AREA_2D computes the area of a circular annulus in 2D.
!
!  Discussion:
!
!    A circular annulus with center (XC,YC), inner radius R1 and
!    outer radius R2, is the set of points (X,Y) so that
!
!      R1^2 <= (X-XC)^2 + (Y-YC)^2 <= R2^2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) R1, R2, the inner and outer radii.
!
!  Output:
!
!    real ( kind = rk ) AREA, the area.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) area
  real ( kind = rk ) r1
  real ( kind = rk ) r2
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00

  area = r8_pi * ( r2 + r1 ) * ( r2 - r1 )

  return
end
subroutine annulus_sector_area_2d ( r1, r2, theta1, theta2, area )

!*****************************************************************************80
!
!! ANNULUS_SECTOR_AREA_2D computes the area of an annular sector in 2D.
!
!  Discussion:
!
!    An annular sector with center PC, inner radius R1 and
!    outer radius R2, and angles THETA1, THETA2, is the set of points 
!    P so that
!
!      R1^2 <= (P(1)-PC(1))^2 + (P(2)-PC(2))^2 <= R2^2
!
!    and
!
!      THETA1 <= THETA ( P - PC ) <= THETA2
!
!  Modified:
!
!    02 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) R1, R2, the inner and outer radii.
!
!    real ( kind = rk ) THETA1, THETA2, the angles.
!
!  Output:
!
!    real ( kind = rk ) AREA, the area.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) area
  real ( kind = rk ) r1
  real ( kind = rk ) r2
  real ( kind = rk ) theta1
  real ( kind = rk ) theta2

  area = 0.5D+00 * ( theta2 - theta1 ) * ( r2 + r1 ) * ( r2 - r1 )

  return
end
subroutine annulus_sector_centroid_2d ( pc, r1, r2, theta1, theta2, centroid )

!*****************************************************************************80
!
!! ANNULUS_SECTOR_CENTROID_2D computes the centroid of an annular sector in 2D.
!
!  Discussion:
!
!    An annular sector with center PC, inner radius R1 and
!    outer radius R2, and angles THETA1, THETA2, is the set of points 
!    P so that
!
!      R1^2 <= (P(1)-PC(1))^2 + (P(2)-PC(2))^2 <= R2^2
!
!    and
!
!      THETA1 <= THETA ( P - PC ) <= THETA2
!
!    Thanks to Ed Segall for pointing out a mistake in the computation
!    of the angle THETA associated with the centroid.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 December 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    John Harris, Horst Stocker,
!    Handbook of Mathematics and Computational Science,
!    Springer, 1998, QA40.S76
!
!  Input:
!
!    real ( kind = rk ) PC(2), the center.
!
!    real ( kind = rk ) R1, R2, the inner and outer radii.
!
!    real ( kind = rk ) THETA1, THETA2, the angles.
!
!  Output:
!
!    real ( kind = rk ) CENTROID(2), the centroid.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) centroid(2)
  real ( kind = rk ) pc(2)
  real ( kind = rk ) r
  real ( kind = rk ) r1
  real ( kind = rk ) r2
  real ( kind = rk ) theta
  real ( kind = rk ) theta1
  real ( kind = rk ) theta2

  theta = theta2 - theta1

  r = 4.0D+00 * sin ( theta / 2.0D+00 ) / ( 3.0D+00 * theta ) &
    * ( r1 * r1 + r1 * r2 + r2 * r2 ) / ( r1 + r2 )

  centroid(1) = pc(1) + r * cos ( theta1 + theta / 2.0D+00 )
  centroid(2) = pc(2) + r * sin ( theta1 + theta / 2.0D+00 )

  return
end
subroutine ball01_sample_2d ( p )

!*****************************************************************************80
!
!! BALL01_SAMPLE_2D picks a random point in the unit ball in 2D.
!
!  Discussion:
!
!    The unit ball is the set of points P such that
!
!      P(1) * P(1) + P(2) * P(2) <= 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Output:
!
!    real ( kind = rk ) P(2), a random point in the unit ball.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) r
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) theta
  real ( kind = rk ) u(dim_num)

  call random_number ( harvest = u(1:dim_num) )

  r = sqrt ( u(1) )
  theta = 2.0D+00 * r8_pi * u(2)

  p(1) = r * cos ( theta )
  p(2) = r * sin ( theta )

  return
end
subroutine ball01_sample_3d ( p )

!*****************************************************************************80
!
!! BALL01_SAMPLE_3D picks a random point in the unit ball in 3D.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Output:
!
!    real ( kind = rk ) P(3), the sample point.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) phi
  real ( kind = rk ) r
  real ( kind = rk ) r8_acos
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) theta
  real ( kind = rk ) u(dim_num)
  real ( kind = rk ) vdot

  call random_number ( harvest = u(1:dim_num) )
!
!  Pick a uniformly random VDOT, which must be between -1 and 1.
!  This represents the dot product of the random vector with the Z unit vector.
!
!  Note: this works because the surface area of the sphere between
!  Z and Z + dZ is independent of Z.  So choosing Z uniformly chooses
!  a patch of area uniformly.
!
  vdot = 2.0D+00 * u(1) - 1.0D+00

  phi = r8_acos ( vdot )
!
!  Pick a uniformly random rotation between 0 and 2 Pi around the
!  axis of the Z vector.
!
  theta = 2.0D+00 * r8_pi * u(2)
!
!  Pick a random radius R.
!
  r = u(3) ** ( 1.0D+00 / 3.0D+00 )

  p(1) = r * cos ( theta ) * sin ( phi )
  p(2) = r * sin ( theta ) * sin ( phi )
  p(3) = r * cos ( phi )

  return
end
subroutine ball01_sample_nd ( dim_num, p )

!*****************************************************************************80
!
!! BALL01_SAMPLE_ND picks a random point in the unit ball in ND.
!
!  Discussion:
!
!    N-1 random Givens rotations are applied to the point ( 1, 0, 0, ..., 0 ).
!
!    The I-th Givens rotation is in the plane of coordinate axes I and I+1,
!    and has the form:
!
!     [ cos ( theta )  - sin ( theta ) ] * x(i)      = x'(i)
!     [ sin ( theta )    cos ( theta ) ]   x(i+1)      x'(i+1)
!
!    Finally, a scaling is applied to set the point at a distance R
!    from the origin, in a way that results in a uniform distribution.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer DIM_NUM, the spatial dimension.
!
!  Output:
!
!    real ( kind = rk ) P(N), the random point.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num

  integer i
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) pi
  real ( kind = rk ) r
  real ( kind = rk ) random_cosine
  real ( kind = rk ) random_sign
  real ( kind = rk ) random_sine
 
  p(1) = 1.0D+00
  p(2:dim_num) = 0.0D+00

  do i = 1, dim_num - 1

    call random_number ( harvest = r )
    random_cosine = 2.0D+00 * r - 1.0D+00
    call random_number ( harvest = r )
    random_sign = real ( 2 * int ( 2.0D+00 * r ) - 1, kind = rk )
    call random_number ( harvest = r )
    random_sine = random_sign * sqrt ( 1.0D+00 - random_cosine * random_cosine )

    pi = p(i)
    p(i  ) = random_cosine * pi
    p(i+1) = random_sine   * pi

  end do

  call random_number ( harvest = r )

  r = r ** ( 1.0D+00 / real ( dim_num, kind = rk ) )

  p(1:dim_num) = r * p(1:dim_num)

  return
end
function ball01_volume ( )

!*****************************************************************************80
!
!! BALL01_VOLUME returns the volume of the unit ball in 3D.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 January 2018
!
!  Author:
!
!    John Burkardt
!
!  Output:
!
!    real ( kind = rk ) BALL01_VOLUME_3D, the volume.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) ball01_volume
  real ( kind = rk ) r
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  
  r = 1.0D+00
  ball01_volume = 4.0D+00 * r8_pi * r ** 3 / 3.0D+00

  return
end
subroutine basis_map_3d ( u, v, a, ierror )

!*****************************************************************************80
!
!! BASIS_MAP_3D computes the matrix which maps one basis to another in 3D.
!
!  Discussion:
!
!    As long as the column vectors U1, U2 and U3 are linearly independent,
!    a matrix A will be computed that maps U1 to V1, U2 to V2, and
!    U3 to V3, where V1, V2 and V3 are the columns of V.
!
!    Depending on the values of the vectors, A may represent a
!    rotation, reflection, dilation, projection, or a combination of these
!    basic linear transformations.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) U(3,3), the columns of U are the three 
!    "domain" or "preimage" vectors, which should be linearly independent.
!
!    real ( kind = rk ) V(3,3), the columns of V are the three 
!    "range" or "image" vectors.
!
!  Output:
!
!    real ( kind = rk ) A(3,3), a matrix with the property that 
!    A * U1 = V1, A * U2 = V2 and A * U3 = V3.
!
!    integer IERROR, error flag.
!    0, no error occurred.
!    nonzero, the matrix [ U1 | U2 | U3 ] is exactly singular.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(3,3)
  real ( kind = rk ) b(3,3)
  real ( kind = rk ) c(3,3)
  real ( kind = rk ) det
  integer ierror
  real ( kind = rk ) u(3,3)
  real ( kind = rk ) v(3,3)

  ierror = 0
!
!  Compute C = the inverse of [ U1 | U2 | U3 ].
!
  b(1:3,1:3) = u(1:3,1:3)

  call r8mat_inverse_3d ( b, c, det )

  if ( det == 0.0D+00 ) then
    ierror = 1
    return
  end if
!
!  A = [ V1 | V2 | V3 ] * inverse [ U1 | U2 | U3 ].
!
  a(1:3,1:3) = matmul ( v(1:3,1:3), c(1:3,1:3) )

  return
end
function box_contains_point_2d ( p1, p2, p )

!*****************************************************************************80
!
!! BOX_CONTAINS_POINT_2D determines if a point is inside a box in 2D.
!
!  Discussion:
!
!    A box in 2D is a rectangle with sides aligned on coordinate
!    axes.  It can be described by its low and high corners, P1 and P2
!    as the set of points P satisfying:
!
!      P1(1:2) <= P(1:2) <= P2(1:2).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 June 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(2), P2(2), the low and high 
!    corners of the box.
!
!    Input, real ( kind = rk ) P(2), the point to be checked.
!
!    Output, logical BOX_CONTAINS_POINT_2D, is TRUE if the point 
!    is inside the box.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  logical box_contains_point_2d
  real ( kind = rk ) p(2)
  real ( kind = rk ) p1(2)
  real ( kind = rk ) p2(2)

  if ( p(1)  < p1(1) .or. &
       p2(1) < p(1)  .or. &
       p(2)  < p1(2) .or. &
       p2(2) < p(2) ) then
    box_contains_point_2d = .false.
  else
    box_contains_point_2d = .true.
  end if

  return
end
function box_contains_point_nd ( dim_num, p1, p2, p )

!*****************************************************************************80
!
!! BOX_CONTAINS_POINT_ND determines if a point is inside a box in ND.
!
!  Discussion:
!
!    A box is a rectangle with sides aligned on coordinate
!    axes.  It can be described by its low and high corners, P1 and P2
!    as the set of points P satisfying:
!
!      P1(1:DIM_NUM) <= P(1:DIM_NUM) <= P2(1:DIM_NUM).
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
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = rk ) P1(DIM_NUM), P2(DIM_NUM), the low and high 
!    corners of the box.
!
!    Input, real ( kind = rk ) P(DIM_NUM), the point to be checked.
!
!    Output, logical BOX_CONTAINS_POINT_ND, is TRUE if the point 
!    is inside the box.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num

  logical box_contains_point_nd
  integer i
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)

  box_contains_point_nd = .false.

  do i = 1, dim_num
    if ( p(i) < p1(i) .or. p2(i) < p(i) ) then
      return
    end if
  end do

  box_contains_point_nd = .true.

  return
end
function box_contains_segment_nd ( dim_num, p1, p2, pa, pb  )

!*****************************************************************************80
!
!! BOX_CONTAINS_SEGMENT_ND reports if a box contains a line segment in ND.
!
!  Discussion:
!
!    A box is assumed to be a rectangle with sides aligned on coordinate
!    axes.  It can be described by its low and high corners, P1 and P2
!    as the set of points P satisfying:
!
!      P1(1:DIM_NUM) <= P(1:DIM_NUM) <= P2(1:DIM_NUM).
!
!    A line segment is the finite portion of a line that lies between
!    two points.
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
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = rk ) P1(DIM_NUM), P2(DIM_NUM), the low and high corners
!    of the box.
!
!    Input, real ( kind = rk ) PA(DIM_NUM), PB(DIM_NUM), the endpoints of the 
!    line segment.
!
!    Output, logical BOX_CONTAINS_SEGMENT_ND, is TRUE if the box 
!    contains the line segment.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num

  logical box_contains_segment_nd
  logical box_contains_point_nd
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) pa(dim_num)
  real ( kind = rk ) pb(dim_num)

  box_contains_segment_nd = .false.

  if ( .not. box_contains_point_nd ( dim_num, p1, p2, pa ) ) then
    return
  end if

  if ( .not. box_contains_point_nd ( dim_num, p1, p2, pb ) ) then
    return
  end if

  box_contains_segment_nd = .true.

  return
end
subroutine box_ray_int_2d ( p1, p2, pa, pb, pint )

!*****************************************************************************80
!
!! BOX_RAY_INT_2D: intersection ( box, ray ) in 2D.
!
!  Discussion:
!
!    A box in 2D is a rectangle with sides aligned on coordinate
!    axes.  It can be described by its low and high corners, P1 and P2
!    as the set of points P satisfying:
!
!      P1(1:2) <= P(1:2) <= P2(1:2).
!
!    The origin of the ray is assumed to be inside the box.  This
!    guarantees that the ray will intersect the box in exactly one point.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(2), P2(2), the low and high corners of the box.
!
!    Input, real ( kind = rk ) PA(2), the origin of the ray, which should be
!    inside the box.
!
!    Input, real ( kind = rk ) PB(2), a second point on the ray.
!
!    Output, real ( kind = rk ) PINT(2), the point on the box intersected 
!    by the ray.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  logical inside
  integer ival
  real ( kind = rk ) p1(2)
  real ( kind = rk ) p2(2)
  real ( kind = rk ) pa(2)
  real ( kind = rk ) pb(2)
  real ( kind = rk ) pc(2)
  real ( kind = rk ) pd(2)
  real ( kind = rk ) pint(2)
  integer side

  do side = 1, 4

    if ( side == 1 ) then
      pd(1:2) = (/ p1(1), p1(2) /)
      pc(1:2) = (/ p2(1), p1(2) /)
    else if ( side == 2 ) then
      pd(1:2) = (/ p2(1), p1(2) /)
      pc(1:2) = (/ p2(1), p2(2) /)
    else if ( side == 3 ) then
      pd(1:2) = (/ p2(1), p2(2) /)
      pc(1:2) = (/ p1(1), p2(2) /)
    else if ( side == 4 ) then
      pd(1:2) = (/ p1(1), p2(2) /)
      pc(1:2) = (/ p1(1), p1(2) /)
    end if

    call angle_contains_point_2d ( pc, pa, pd, pb, inside )

    if ( inside ) then
      exit
    end if

    if ( side == 4 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BOX_RAY_INT_2D - Fatal error!'
      write ( *, '(a)' ) '  No intersection could be found.'
      stop 1
    end if

  end do

  call lines_exp_int_2d ( pa, pb, pc, pd, ival, pint )

  return
end
subroutine box_segment_clip_2d ( p1, p2, pa, pb, ival )

!*****************************************************************************80
!
!! BOX_SEGMENT_CLIP_2D uses a box to clip a line segment in 2D.
!
!  Discussion:
!
!    A box in 2D is a rectangle with sides aligned on coordinate
!    axes.  It can be described by its low and high corners, P1 and P2
!    as the set of points P satisfying:
!
!      P1(1:2) <= P(1:2) <= P2(1:2).
!
!    A line segment is the finite portion of a line that lies between
!    two points.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(2), P2(2), the low and high corners of the box.
!
!    Input/output, real ( kind = rk ) PA(2), PB(2); on input, the endpoints 
!    of a line segment.  On output, the endpoints of the portion of the
!    line segment that lies inside the box.  However, if no part of the
!    initial line segment lies inside the box, the output value is the
!    same as the input value.
!
!    Output, integer IVAL:
!    -1, no part of the line segment is within the box.
!     0, no clipping was necessary.
!     1, PA was clipped.
!     2, PB was clipped.
!     3, PA and PB were clipped.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  logical clip_a
  logical clip_b
  integer ival
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) pa(dim_num)
  real ( kind = rk ) pb(dim_num)
  real ( kind = rk ) q(dim_num)

  clip_a = .false.
  clip_b = .false.
!
!  Require that XMIN <= X.
!
  if ( pa(1) < p1(1) .and. pb(1) < p1(1) ) then
    ival = -1
    return
  end if

  if ( pa(1) < p1(1) .and. p1(1) <= pb(1) ) then
    q(1) = p1(1)
    q(2) = pa(2) + ( pb(2) - pa(2) ) * ( q(1) - pa(1) ) / ( pb(1) - pa(1) )
    pa(1:2) = q(1:2)
    clip_a = .true.
  else if ( p1(1) <= pa(1) .and. pb(1) < p1(1) ) then
    q(1) = p1(1)
    q(2) = pa(2) + ( pb(2) - pa(2) ) * ( q(1) - pa(1) ) / ( pb(1) - pa(1) )
    pb(1:2) = q(1:2)
    clip_b = .true.
  end if
!
!  Require that X <= XMAX.
!
  if ( p2(1) < pa(1) .and. p2(1) < pb(1) ) then
    ival = -1
    return
  end if

  if ( p2(1) < pa(1) .and. pb(1) <= p2(1) ) then
    q(1) = p2(1)
    q(2) = pa(2) + ( pb(2) - pa(2) ) * ( q(1) - pa(1) ) / ( pb(1) - pa(1) )
    pa(1:2) = q(1:2)
    clip_a = .true.
  else if ( pa(1) <= p2(1) .and. p2(1) < pb(1) ) then
    q(1) = p2(1)
    q(2) = pa(2) + ( pb(2) - pa(2) ) * ( q(1) - pa(1) ) / ( pb(1) - pa(1) )
    pb(1:2) = q(1:2)
    clip_b = .true.
  end if
!
!  Require that YMIN <= Y.
!
  if ( pa(2) < p1(2) .and. pb(2) < p1(2) ) then
    ival = -1
    return
  end if

  if ( pa(2) < p1(2) .and. p1(2) <= pb(2) ) then
    q(2) = p1(2)
    q(1) = pa(1) + ( pb(1) - pa(1) ) * ( q(2) - pa(2) ) / ( pb(2) - pa(2) )
    pa(1:2) = q(1:2)
    clip_a = .true.
  else if ( p1(2) <= pa(2) .and. pb(2) < p1(2) ) then
    q(2) = p1(2)
    q(1) = pa(1) + ( pb(1) - pa(1) ) * ( q(2) - pa(2) ) / ( pb(2) - pa(2) )
    pb(1:2) = q(1:2)
    clip_b = .true.
  end if
!
!  Require that Y <= YMAX.
!
  if ( p2(2) < pa(2) .and. p2(2) < pb(2) ) then
    ival = -1
    return
  end if

  if ( p2(2) < pa(2) .and. pb(2) <= p2(2) ) then
    q(2) = p2(2)
    q(1) = pa(1) + ( pb(1) - pa(1) ) * ( q(2) - pa(2) ) / ( pb(2) - pa(2) )
    pa(1:2) = q(1:2)
    clip_a = .true.
  else if ( pa(2) <= p2(2) .and. p2(2) < pb(2) ) then
    q(2) = p2(2)
    q(1) = pa(1) + ( pb(1) - pa(1) ) * ( p2(2) - pa(2) ) / ( pb(2) - pa(2) )
    pb(1:2) = q(1:2)
    clip_b = .true.
  end if

  ival = 0

  if ( clip_a ) then
    ival = ival + 1
  end if

  if ( clip_b ) then
    ival = ival + 2
  end if

  return
end
function box01_contains_point_2d ( p )

!*****************************************************************************80
!
!! BOX01_CONTAINS_POINT_2D determines if a point is inside the unit box in 2D.
!
!  Discussion:
!
!    A unit box is assumed to be a rectangle with sides aligned on coordinate
!    axes.  It can be described as the set of points P satisfying:
!
!      0.0 <= P(1:DIM_NUM) <= 1.0
!
!      0.0 <= P(1:2) <= 1.0 
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P(2), the point to be checked.
!
!    Output, logical BOX01_CONTAINS_POINT_2D, is TRUE if the
!    point is  inside the box.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  logical box01_contains_point_2d
  real ( kind = rk ) p(2)

  box01_contains_point_2d = &
    all ( 0.0D+00 <= p(1:2) ) .and. all ( p(1:2) <= 1.0D+00 )

  return
end
function box01_contains_point_nd ( dim_num, p )

!*****************************************************************************80
!
!! BOX01_CONTAINS_POINT_ND determines if a point is inside the unit box in ND.
!
!  Discussion:
!
!    A unit box is assumed to be a rectangle with sides aligned on coordinate
!    axes.  It can be described as the set of points P satisfying:
!
!      0.0 <= P(1:DIM_NUM) <= 1.0
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 June 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = rk ) P(DIM_NUM), the point to be checked.
!
!    Output, logical BOX_01_CONTAINS_POINT_ND, is TRUE 
!    if the point is inside the box.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num

  logical box01_contains_point_nd
  real ( kind = rk ) p(dim_num)

  box01_contains_point_nd = &
    all ( 0.0D+00 <= p(1:dim_num) ) .and. all ( p(1:dim_num) <= 1.0D+00 )

  return
end
subroutine circle_arc_point_near_2d ( r, pc, theta1, theta2, p, pn, &
  dist )

!*****************************************************************************80
!
!! CIRCLE_ARC_POINT_NEAR_2D : nearest point on a circular arc.
!
!  Discussion:
!
!    A circular arc is defined by the portion of a circle (R,C)
!    between two angles (THETA1,THETA2).
!
!    Thus, a point P on a circular arc satisfies
!
!      ( P(1) - PC(1) ) * ( P(1) - PC(1) ) 
!    + ( P(2) - PC(2) ) * ( P(2) - PC(2) ) = R * R
!
!    and
!
!      Theta1 <= Theta <= Theta2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) R, the radius of the circle.
!
!    real ( kind = rk ) PC(2), the center of the circle.
!
!    real ( kind = rk ) THETA1, THETA2, the angles defining the arc,
!    in radians.  Normally, THETA1 < THETA2.
!
!    real ( kind = rk ) P(2), the point to be checked.
!
!  Output:
!
!    real ( kind = rk ) PN(2), a point on the circular arc which is
!    nearest to the point.
!
!    real ( kind = rk ) DIST, the distance to the nearest point.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) dist
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) pn(dim_num)
  real ( kind = rk ) r
  real ( kind = rk ) r2
  real ( kind = rk ) r8_atan
  real ( kind = rk ) r8_modp
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) theta
  real ( kind = rk ) theta1
  real ( kind = rk ) theta2
!
!  Special case, the zero circle.
!
  if ( r == 0.0D+00 ) then
    pn(1:dim_num) = pc(1:dim_num)
    dist = sqrt ( sum ( ( p(1:dim_num) - pn(1:dim_num) )**2 ) )
    return
  end if
!
!  Determine the angle made by the point.
!
  theta = r8_atan ( p(2) - pc(2), p(1) - pc(1) )
!
!  If the angle is between THETA1 and THETA2, then you can
!  simply project the point onto the arc.
!
  if ( r8_modp ( theta  - theta1,  2.0D+00 * r8_pi ) <= &
       r8_modp ( theta2 - theta1,  2.0D+00 * r8_pi ) ) then

    r2 = sqrt ( sum ( ( p(1:dim_num) - pc(1:dim_num) )**2 ) )

    pn(1:dim_num) = pc(1:dim_num) + ( p(1:dim_num) - pc(1:dim_num) ) * r / r2
!
!  Otherwise, if the angle is less than the negative of the
!  average of THETA1 and THETA2, it's on the side of the arc
!  where the endpoint associated with THETA2 is closest.
!
  else if ( r8_modp ( theta - 0.5D+00 * ( theta1 + theta2 ), 2.0D+00 * r8_pi ) &
    <= r8_pi ) then

    pn(1:dim_num) = pc(1:dim_num) + r * (/ cos ( theta2 ), sin ( theta2 ) /)
!
!  Otherwise, the endpoint associated with THETA1 is closest.
!
  else

    pn(1:dim_num) = pc(1:dim_num) + r * (/ cos ( theta1 ), sin ( theta1 ) /)

  end if

  dist = sqrt ( sum ( ( p(1:dim_num) - pn(1:dim_num) )**2 ) )

  return
end
subroutine circle_area_2d ( r, area )

!*****************************************************************************80
!
!! CIRCLE_AREA_2D computes the area of a circle in 2D.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) R, the radius of the circle.
!
!  Output:
!
!    real ( kind = rk ) AREA, the area of the circle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) area
  real ( kind = rk ) r
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00

  area = r8_pi * r * r

  return
end
subroutine circle_dia2imp_2d ( p1, p2, r, pc )

!*****************************************************************************80
!
!! CIRCLE_DIA2IMP_2D converts a diameter to an implicit circle in 2D.
!
!  Discussion:
!
!    The diameter form of a circle is:
!
!      P1 and P2 are the endpoints of a diameter.
!
!    Points P on an implicit circle in 2D satisfy the equation:
!
!      ( P(1) - PC(1) )^2 + ( P(2) - PC(2) )^2 = R^2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) P1(2), P2(2), two points that are the 
!    endpoints of a diameter of the circle.
!
!  Output:
!
!    real ( kind = rk ) R, the radius of the circle.
!
!    real ( kind = rk ) PC(2), the center of the circle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) r

  r = 0.5D+00 * sqrt ( sum ( ( p2(1:dim_num) - p1(1:dim_num) )**2 ) )

  pc(1:dim_num) = 0.5D+00 * ( p1(1:dim_num) + p2(1:dim_num) )

  return
end
subroutine circle_exp_contains_point_2d ( p1, p2, p3, p, inside )

!*****************************************************************************80
!
!! CIRCLE_EXP_CONTAINS_POINT_2D: explicit circle contains a point in 2D.
!
!  Discussion:
!
!    The explicit form of a circle in 2D is:
!
!      The circle passing through points P1, P2 and P3.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) P1(2), P2(2), P3(2), three points on a circle.
!
!    real ( kind = rk ) P(2), the point to test.
!
!  Output:
!
!    integer INSIDE, reports the result:
!   -1, the three points are distinct and noncolinear,
!    and P lies inside the circle.
!    0, the three points are distinct and noncolinear,
!    and P lies on the circle.
!    1, the three points are distinct and noncolinear,
!    and P lies outside the circle.
!    2, the three points are distinct and colinear,
!    and P lies on the line.
!    3, the three points are distinct and colinear,
!    and P does not lie on the line.
!    4, two points are distinct, and P lies on the line.
!    5, two points are distinct, and P does not lie on the line.
!    6, all three points are equal, and P is equal to them,
!    7, all three points are equal, and P is not equal to them.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) a(4,4)
  real ( kind = rk ) det
  real ( kind = rk ) r8mat_det_4d
  integer inside
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) p3(dim_num)
!
!  P1 = P2?
!
  if ( all ( p1(1:dim_num) == p2(1:dim_num) ) ) then

    if ( all ( p1(1:dim_num) == p3(1:dim_num) ) ) then

      if ( all ( p1(1:dim_num) == p(1:dim_num) ) ) then
        inside = 6
      else
        inside = 7
      end if

    else

      det = ( p1(1) - p3(1) ) * ( p(2)  - p3(2) ) &
          - ( p(1)  - p3(1) ) * ( p1(2) - p3(2) )

      if ( det == 0.0D+00 ) then
        inside = 4
      else
        inside = 5
      end if
    end if

    return

  end if
!
!  P1 does not equal P2.  Does P1 = P3?
!
  if ( all ( p1(1:dim_num) == p3(1:dim_num) ) ) then

    det = ( p1(1) - p2(1) ) * ( p(2)  - p2(2) ) &
        - ( p(1)  - p2(1) ) * ( p1(2) - p2(2) )

    if ( det == 0.0D+00 ) then
      inside = 4
    else
      inside = 5
    end if

    return

  end if
!
!  The points are distinct.  Are they colinear?
!
  det = ( p1(1) - p2(1) ) * ( p3(2) - p2(2) ) &
      - ( p3(1) - p2(1) ) * ( p1(2) - p2(2) )

  if ( det == 0.0D+00 ) then

    det = ( p1(1) - p2(1) ) * ( p(2)  - p2(2) ) &
        - ( p(1)  - p2(1) ) * ( p1(2) - p2(2) )

    if ( det == 0.0D+00 ) then
      inside = 2
    else
      inside = 3
    end if

    return

  end if
!
!  The points are distinct and non-colinear.
!
!  Compute the determinant
!
  a(1,1) = p1(1)
  a(1,2) = p1(2)
  a(1,3) = p1(1) * p1(1) + p1(2) * p1(2)
  a(1,4) = 1.0D+00

  a(2,1) = p2(1)
  a(2,2) = p2(2)
  a(2,3) = p2(1) * p2(1) + p2(2) * p2(2)
  a(2,4) = 1.0D+00

  a(3,1) = p3(1)
  a(3,2) = p3(2)
  a(3,3) = p3(1) * p3(1) + p3(2) * p3(2)
  a(3,4) = 1.0D+00

  a(4,1) = p(1)
  a(4,2) = p(2)
  a(4,3) = p(1) * p(1) + p(2) * p(2)
  a(4,4) = 1.0D+00

  det = r8mat_det_4d ( a )

  if ( det < 0.0D+00 ) then
    inside = 1
  else if ( det == 0.0D+00 ) then
    inside = 0
  else
    inside = -1
  end if

  return
end
subroutine circle_exp2imp_2d ( p1, p2, p3, r, pc )

!*****************************************************************************80
!
!! circle_exp2imp_2d converts a circle from explicit to implicit form in 2D.
!
!  Discussion:
!
!    The explicit form of a circle in 2D is:
!
!      The circle passing through points P1, P2 and P3.
!
!    Points P on an implicit circle in 2D satisfy the equation:
!
!      ( P(1) - PC(1) )^2 + ( P(2) - PC(2) )^2 = R^2
!
!    Any three distinct points define a circle, as long as they don't lie 
!    on a straight line.  (If the points do lie on a straight line, we 
!    could stretch the definition of a circle to allow an infinite radius 
!    and a center at some infinite point.)
!
!    The diameter of the circle can be found by solving a 2 by 2 linear system.
!    This is because the vectors P2 - P1 and P3 - P1 are secants of the circle,
!    and each forms a right triangle with the diameter.  Hence, the dot product
!    of P2 - P1 with the diameter is equal to the square of the length
!    of P2 - P1, and similarly for P3 - P1.  These two equations determine the
!    diameter vector originating at P1.
!
!    If all three points are equal, return a circle of radius 0 and 
!    the obvious center.
!
!    If two points are equal, return a circle of radius half the distance
!    between the two distinct points, and center their average.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    08 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph ORourke,
!    Computational Geometry,
!    Second Edition,
!    Cambridge, 1998,
!    ISBN: 0521649765,
!    LC: QA448.D38.
!
!  Input:
!
!    real ( kind = rk ) P1(2), P2(2), P3(2), three points on the circle.
!
!  Output:
!
!    real ( kind = rk ) R, the radius of the circle.  Normally, R will
!    be positive.  R will be (meaningfully) zero if all three points are 
!    equal.  If two points are equal, R is returned as the distance between
!    two nonequal points.  R is returned as -1 in the unlikely event that 
!    the points are numerically collinear; philosophically speaking, R 
!    should actually be "infinity" in this case.
!
!    real ( kind = rk ) PC(2), the center of the circle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) e
  real ( kind = rk ) f
  real ( kind = rk ) g
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) p3(dim_num)
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) r
!
!  If all three points are equal, then the
!  circle of radius 0 and center P1 passes through the points.
!
  if ( all ( p1(1:dim_num) == p2(1:dim_num) ) .and. &
       all ( p1(1:dim_num) == p3(1:dim_num) ) ) then
    r = 0.0D+00
    pc(1:dim_num) = p1(1:dim_num)
    return
  end if
!
!  If exactly two points are equal, then the circle is defined as
!  having the obvious radius and center.
!
       if ( all ( p1(1:dim_num) == p2(1:dim_num) ) ) then

    r = 0.5D+00 * sqrt ( sum ( ( p1(1:dim_num) - p3(1:dim_num) )**2 ) )
    pc(1:dim_num) = 0.5D+00 * ( p1(1:dim_num) + p3(1:dim_num)  )
    return

  else if ( all ( p1(1:dim_num) == p3(1:dim_num) ) ) then

    r = 0.5D+00 * sqrt ( sum ( ( p1(1:dim_num) - p2(1:dim_num) )**2 ) )
    pc(1:dim_num) = 0.5D+00 * ( p1(1:dim_num) + p2(1:dim_num)  )
    return

  else if ( all ( p2(1:dim_num) == p3(1:dim_num) ) ) then

    r = 0.5D+00 * sqrt ( sum ( ( p1(1:dim_num) - p2(1:dim_num) )**2 ) )
    pc(1:dim_num) = 0.5D+00 * ( p1(1:dim_num) + p2(1:dim_num)  )
    return

  end if
!
!  We check for collinearity.  A more useful check would compare the
!  absolute value of G to a small quantity.
!
  e = ( p2(1) - p1(1) ) * ( p1(1) + p2(1) ) &
    + ( p2(2) - p1(2) ) * ( p1(2) + p2(2) )

  f = ( p3(1) - p1(1) ) * ( p1(1) + p3(1) ) &
    + ( p3(2) - p1(2) ) * ( p1(2) + p3(2) )

  g = ( p2(1) - p1(1) ) * ( p3(2) - p2(2) ) &
    - ( p2(2) - p1(2) ) * ( p3(1) - p2(1) )

  if ( g == 0.0D+00 ) then
    pc(1:2) = (/ 0.0D+00, 0.0D+00 /)
    r = -1.0D+00
    return
  end if
!
!  The center is halfway along the diameter vector from P1.
!
  pc(1) = 0.5D+00 * ( ( p3(2) - p1(2) ) * e - ( p2(2) - p1(2) ) * f ) / g
  pc(2) = 0.5D+00 * ( ( p2(1) - p1(1) ) * f - ( p3(1) - p1(1) ) * e ) / g
!
!  Knowing the center, the radius is now easy to compute.
!
  r = sqrt ( sum ( ( p1(1:dim_num) - pc(1:dim_num) )**2 ) )

  return
end
subroutine circle_imp_contains_point_2d ( r, pc, p, inside )

!*****************************************************************************80
!
!! circle_imp_contains_point_2d: implicit circle contains a point in 2D?
!
!  Discussion:
!
!    Points P on an implicit circle in 2D satisfy the equation:
!
!      ( P(1) - PC(1) )^2 + ( P(2) - PC(2) )^2 = R^2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) R, the radius of the circle.
!
!    real ( kind = rk ) PC(2), the center of the circle.
!
!    real ( kind = rk ) P(2), the point to be checked.
!
!  Output:
!
!    logical INSIDE, is TRUE if the point is inside or
!    on the circle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  logical inside
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) r

  if ( ( p(1) - pc(1) ) * ( p(1) - pc(1) ) &
     + ( p(2) - pc(2) ) * ( p(2) - pc(2) ) <= r * r ) then
    inside = .true.
  else
    inside = .false.
  end if

  return
end
subroutine circle_imp_line_exp_dist_2d ( r, pc, p1, p2, dist )

!*****************************************************************************80
!
!! circle_imp_line_exp_dist_2d: distance ( implicit circle, explicit line ) in 2D.
!
!  Discussion:
!
!    The distance is zero if the line intersects the circle.
!
!    Points P on an implicit circle in 2D satisfy the equation:
!
!      ( P(1) - PC(1) )^2 + ( P(2) - PC(2) )^2 = R^2
!
!    The explicit form of a line in 2D is:
!
!      the line through the points P1 and P2.
!
!    The distance between the circle and the line is zero if
!    and only if they intersect.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) R, the radius of the circle.
!
!    real ( kind = rk ) PC(2), the center of the circle.
!
!    real ( kind = rk ) P1(2), P2(2), two points on the line.
!
!  Output:
!
!    real ( kind = rk ) DIST, the distance of the line to the circle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) dist
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) r

  call line_exp_point_dist_2d ( p1, p2, pc, dist )

  dist = dist - r

  if ( dist < 0.0D+00 ) then
    dist = 0.0D+00
  end if

  return
end
subroutine circle_imp_line_exp_intersect ( r, pc, p0, p1, int_num, p )

!*****************************************************************************80
!
!! circle_imp_line_exp_intersect: ( imp circle, exp line ) intersection.
!
!  Discussion:
!
!    Points P on an implicit circle in 2D satisfy the equation:
!
!      ( P(1) - PC(1) )^2 + ( P(2) - PC(2) )^2 = R^2
!
!    An explicit line in 2D is:
!
!      P(t) = (1-t) * P0 + t * P1
!
!    where P0 and P1 are two distinct points on the line.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 September 2020
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Eric W Weisstein, 
!    "Circle-Line Intersection." 
!    From MathWorld--A Wolfram Web Resource. 
!    https://mathworld.wolfram.com/Circle-LineIntersection.html
!
!  Input:
!
!    real ( kind = rk ) R, the radius of the circle.
!
!    real ( kind = rk ) PC(2), the center of the circle.
!
!    real ( kind = rk ) P0(2), P1(2), two points on the explicit line.
!
!  Output:
!
!    integer INT_NUM, the number of intersections.
!    INT_NUM will be 0, 1 or 2.
!
!    real ( kind = rk ) P(2,INT_NUM), the intersecting points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) d
  real ( kind = rk ) disc
  real ( kind = rk ) dr
  real ( kind = rk ) dx
  real ( kind = rk ) dy
  integer int_num
  real ( kind = rk ) p(2,2)
  real ( kind = rk ) p0(2)
  real ( kind = rk ) p1(2)
  real ( kind = rk ) pc(2)
  real ( kind = rk ) r
  real ( kind = rk ) x1
  real ( kind = rk ) x2
  real ( kind = rk ) xm
  real ( kind = rk ) xp
  real ( kind = rk ) y1
  real ( kind = rk ) y2
  real ( kind = rk ) ym
  real ( kind = rk ) yp
!
!  Solve the problem in the coordinate system where the circle 
!  has center (0,0).
!
  x1 = p0(1) - pc(1)
  y1 = p0(2) - pc(2)

  x2 = p1(1) - pc(1)
  y2 = p1(2) - pc(2)

  dx = x2 - x1
  dy = y2 - y1

  dr = sqrt ( dx * dx + dy * dy )
  d = x1 * y2 - x2 * y1

  disc = r**2 * dr**2 - d**2

  if ( disc < 0.0 ) then
    int_num = 0
    p(1,1) = pc(1)
    p(2,1) = pc(2)
    p(1,2) = pc(1)
    p(2,2) = pc(2)
  else if ( disc == 0.0 ) then
    int_num = 1
    xp = (   d * dy ) / dr**2
    yp = ( - d * dx ) / dr**2
    p(1,1) = xp + pc(1)
    p(2,1) = yp + pc(2)
    p(1,2) = pc(1)
    p(2,2) = pc(2)
  else
    int_num = 2
    xp = (   d * dy + sign ( 1.0D+00, dy ) * dx * sqrt ( disc ) ) / dr**2
    xm = (   d * dy - sign ( 1.0D+00, dy ) * dx * sqrt ( disc ) ) / dr**2
    yp = ( - d * dx + abs ( dy )                * sqrt ( disc ) ) / dr**2
    ym = ( - d * dx - abs ( dy )                * sqrt ( disc ) ) / dr**2
    p(1,1) = xp + pc(1)
    p(2,1) = yp + pc(2)
    p(1,2) = xm + pc(1)
    p(2,2) = ym + pc(2)
  end if

  return
end
subroutine circle_imp_line_par_int_2d ( r, pc, x0, y0, f, g, int_num, p )

!*****************************************************************************80
!
!! circle_imp_line_par_int_2d: ( imp circle, param line ) intersection in 2D.
!
!  Discussion:
!
!    Points P on an implicit circle in 2D satisfy the equation:
!
!      ( P(1) - PC(1) )^2 + ( P(2) - PC(2) )^2 = R^2
!
!    The parametric form of a line in 2D is:
!
!      X = X0 + F * T
!      Y = Y0 + G * T
!
!    We normalize by always choosing F^2 + G^2 = 1, and F nonnegative.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) R, the radius of the circle.
!
!    real ( kind = rk ) PC(2), the center of the circle.
!
!    real ( kind = rk ) F, G, X0, Y0, the parametric parameters of the line.
!
!  Output:
!
!    integer INT_NUM, the number of intersecting 
!    points found.  INT_NUM will be 0, 1 or 2.
!
!    real ( kind = rk ) P(2,INT_NUM), the intersecting points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) f
  real ( kind = rk ) g
  integer int_num
  real ( kind = rk ) p(dim_num,2)
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) r
  real ( kind = rk ) root
  real ( kind = rk ) t
  real ( kind = rk ) x0
  real ( kind = rk ) y0

  root = r * r * ( f * f + g * g ) - ( f * ( pc(2) - y0 ) &
    - g * ( pc(1) - x0 ) )**2

  if ( root < 0.0D+00 ) then

    int_num = 0

  else if ( root == 0.0D+00 ) then

    int_num = 1

    t = ( f * ( pc(1) - x0 ) + g * ( pc(2) - y0 ) ) / ( f * f + g * g )
    p(1,1) = x0 + f * t
    p(2,1) = y0 + g * t

  else if ( 0.0D+00 < root ) then

    int_num = 2

    t = ( ( f * ( pc(1) - x0 ) + g * ( pc(2) - y0 ) ) &
      - sqrt ( root ) ) / ( f * f + g * g )

    p(1,1) = x0 + f * t
    p(2,1) = y0 + g * t

    t = ( ( f * ( pc(1) - x0 ) + g * ( pc(2) - y0 ) ) &
      + sqrt ( root ) ) / ( f * f + g * g )

    p(1,2) = x0 + f * t
    p(2,2) = y0 + g * t

  end if

  return
end
subroutine circle_imp_point_dist_2d ( r, pc, p, dist )

!*****************************************************************************80
!
!! CIRCLE_IMP_POINT_DIST_2D: distance ( implicit circle, point ) in 2D.
!
!  Discussion:
!
!    The distance is zero if the point is on the circle.
!
!    Points P on an implicit circle in 2D satisfy the equation:
!
!      ( P(1) - PC(1) )^2 + ( P(2) - PC(2) )^2 = R^2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) R, the radius of the circle.
!
!    real ( kind = rk ) PC(2), the center of the circle.
!
!    real ( kind = rk ) P(2), the point to be checked.
!
!  Output:
!
!    real ( kind = rk ) DIST, the distance of the point to the circle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) dist
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) r
  real ( kind = rk ) r2

  r2 = sqrt ( sum ( ( p(1:2) - pc(1:2) )**2 ) )

  dist = abs ( r2 - r )

  return
end
subroutine circle_imp_point_dist_signed_2d ( r, pc, p, dist )

!*****************************************************************************80
!
!! circle_imp_point_dist_signed_2d: signed distance ( imp circle, point ) in 2D.
!
!  Discussion:
!
!    The signed distance is zero if the point is on the circle.
!    The signed distance is negative if the point is inside the circle.
!
!    Points P on an implicit circle in 2D satisfy the equation:
!
!      ( P(1) - PC(1) )^2 + ( P(2) - PC(2) )^2 = R^2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) R, the radius of the circle.
!
!    real ( kind = rk ) PC(2), the center of the circle.
!
!    real ( kind = rk ) P(2), the point to be checked.
!
!  Output:
!
!    real ( kind = rk ) DIST, the signed distance of the point
!    to the circle.  If the point is inside the circle, the signed distance
!    is negative.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) dist
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) r
  real ( kind = rk ) r2

  r2 = sqrt ( sum ( ( p(1:dim_num) - pc(1:dim_num) )**2 ) )

  dist = r2 - r

  return
end
subroutine circle_imp_point_near_2d ( r, pc, p, pn, dist )

!*****************************************************************************80
!
!! circle_imp_point_near_2d: nearest ( implicit circle, point ) in 2D.
!
!  Discussion:
!
!    This routine finds the distance from a point to an implicitly
!    defined circle, and returns the point on the circle that is
!    nearest to the given point.
!
!    If the given point is the center of the circle, than any point
!    on the circle is "the" nearest.
!
!    Points P on an implicit circle in 2D satisfy the equation:
!
!      ( P(1) - PC(1) )^2 + ( P(2) - PC(2) )^2 = R^2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) R, the radius of the circle.
!
!    real ( kind = rk ) PC(2), the center of the circle.
!
!    real ( kind = rk ) P(2), the point to be checked.
!
!  Output:
!
!    real ( kind = rk ) PN(2), the nearest point on the circle.
!
!    real ( kind = rk ) DIST, the distance of the point to the circle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) dist
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) pn(dim_num)
  real ( kind = rk ) r
  real ( kind = rk ) r2

  if ( all ( p(1:dim_num) == pc(1:dim_num) ) ) then
    dist = r
    pn(1:dim_num) = pc(1:dim_num) + r / sqrt ( real ( dim_num, kind = rk ) )
    return
  end if

  r2 = sqrt ( sum ( ( p(1:dim_num) - pc(1:dim_num) )**2 ) )

  dist = abs (  r2 - r )

  pn(1:dim_num) = pc(1:dim_num) + r * ( p(1:dim_num) - pc(1:dim_num) ) / r2

  return
end
subroutine circle_imp_points_2d ( r, pc, n, p )

!*****************************************************************************80
!
!! circle_imp_points_2d returns points on an implicit circle in 2D.
!
!  Discussion:
!
!    The first point is always ( PC(1) + R, PC(2) ), and subsequent
!    points proceed counter clockwise around the circle.
!
!    Points P on an implicit circle in 2D satisfy the equation:
!
!      ( P(1) - PC(1) )^2 + ( P(2) - PC(2) )^2 = R^2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) R, the radius of the circle.
!
!    real ( kind = rk ) PC(2), the center of the circle.
!
!    integer N, the number of points desired.  
!    N must be at least 1.
!
!  Output:
!
!    real ( kind = rk ) P(2,N), the coordinates of points on the circle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: dim_num = 2

  integer j
  real ( kind = rk ) p(2,n)
  real ( kind = rk ) pc(2)
  real ( kind = rk ) r
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) theta

  do j = 1, n
    theta = 2.0D+00 * r8_pi * real ( j - 1, kind = rk ) / real ( n, kind = rk )
    p(1:dim_num,j) = pc(1:dim_num) + r * (/ cos ( theta ), sin ( theta ) /)
  end do

  return
end
subroutine circle_imp_points_3d ( r, pc, nc, n, p )

!*****************************************************************************80
!
!! circle_imp_points_3d returns points on an implicit circle in 3D.
!
!  Discussion:
!
!    Points P on an implicit circle in 3D satisfy the equations:
!
!      ( P(1) - PC(1) )^2 
!    + ( P(2) - PC(2) )^2 
!    + ( P(3) - PC(3) )^2 = R^2
!
!    and
!
!      ( P(1) - PC(1) ) * NC(1) 
!    + ( P(2) - PC(2) ) * NC(2) 
!    + ( P(3) - PC(3) ) * NC(3) = 0
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) R, the radius of the circle.
!
!    Input, real ( kind = rk ) PC(3), the center of the circle.
!
!    Input, real ( kind = rk ) NC(3), a nonzero vector that is normal to
!    the plane of the circle.  It is customary, but not necessary,
!    that this vector have unit norm.
!
!    Input, integer N, the number of points desired.  
!    N must be at least 1.
!
!    Output, real ( kind = rk ) P(3,N), the coordinates of points 
!    on the circle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: dim_num = 3

  integer j
  real ( kind = rk ) n1(dim_num)
  real ( kind = rk ) n2(dim_num)
  real ( kind = rk ) nc(dim_num)
  real ( kind = rk ) p(dim_num,n)
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) r
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) theta
!
!  Get two unit vectors N1 and N2 which are orthogonal to each other,
!  and to NC.
!
  call plane_normal_basis_3d ( pc, nc, n1, n2 )
!
!  Rotate R units away from PC in the plane of N1 and N2.
!
  do j = 1, n

    theta = ( 2.0D+00 * r8_pi * real ( j - 1, kind = rk ) ) &
      / real ( n, kind = rk )

    p(1:dim_num,j) = pc(1:dim_num) &
      + r * ( cos ( theta ) * n1(1:dim_num) &
            + sin ( theta ) * n2(1:dim_num) )

  end do

  return
end
subroutine circle_imp_points_arc_2d ( r, pc, theta1, theta2, n, p )

!*****************************************************************************80
!
!! circle_imp_points_arc_2d: N points on an arc of an implicit circle in 2D.
!
!  Discussion:
!
!    The first point is 
!      ( PC(1) + R * COS ( THETA1 ), PC(2) + R * SIN ( THETA1 ) );
!    The last point is
!      ( PC(1) + R * COS ( THETA2 ), PC(2) + R * SIN ( THETA2 ) );
!    and the intermediate points are evenly spaced in angle between these,
!    and in counter clockwise order.
!
!    Points P on an implicit circle in 2D satisfy the equation:
!
!      ( P(1) - PC(1) )^2 + ( P(2) - PC(2) )^2 = R^2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) R, the radius of the circle.
!
!    real ( kind = rk ) PC(2), the center of the circle.
!
!    real ( kind = rk ) THETA1, THETA2, the angular coordinates of 
!    the first and last points to be drawn, in radians.
!
!    integer N, the number of points desired.  
!    N must be at least 1.
!
!  Output:
!
!    real ( kind = rk ) P(2,N), the points on the circle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: dim_num = 2

  integer i
  real ( kind = rk ) p(dim_num,n)
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) r
  real ( kind = rk ) r8_modp
  real ( kind = rk ) theta
  real ( kind = rk ) theta1
  real ( kind = rk ) theta2
  real ( kind = rk ) theta3
!
!  THETA3 is the smallest angle, no less than THETA1, which
!  coincides with THETA2.
!
  theta3 = theta1 + r8_modp ( theta2 - theta1, 2.0D+00 * r8_pi )

  do i = 1, n

    if ( 1 < n ) then
      theta = ( real ( n - i,     kind = rk ) * theta1   &
              + real (     i - 1, kind = rk ) * theta3 ) &
              / real ( n     - 1, kind = rk )
    else
      theta = 0.5D+00 * ( theta1 + theta3 )
    end if

    p(1:dim_num,i) = pc(1:dim_num) + r * (/ cos ( theta ), sin ( theta ) /)

  end do

  return
end
subroutine circle_imp_print_2d ( r, pc, title )

!*****************************************************************************80
!
!! circle_imp_print_2d prints an implicit circle in 2D.
!
!  Discussion:
!
!    Points P on an implicit circle in 2D satisfy the equation:
!
!      ( P(1) - PC(1) )^2 + ( P(2) - PC(2) )^2 = R^2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) R, the radius of the circle.
!
!    real ( kind = rk ) PC(2), the center of the circle.
!
!    character ( length = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) r
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)'        ) ' '
  write ( *, '(a,g14.6)'  ) '  Radius = ', r
  write ( *, '(a,2g14.6)' ) '  Center = ', pc(1:dim_num)

  return
end
subroutine circle_imp_print_3d ( r, pc, nc, title )

!*****************************************************************************80
!
!! circle_imp_print_3d prints an implicit circle in 3D.
!
!  Discussion:
!
!    Points P on an implicit circle in 3D satisfy the equations:
!
!      ( P(1) - PC(1) )^2 
!    + ( P(2) - PC(2) )^2 
!    + ( P(3) - PC(3) )^2 = R^2
!
!    and
!
!      ( P(1) - PC(1) ) * NC(1) 
!    + ( P(2) - PC(2) ) * NC(2) 
!    + ( P(3) - PC(3) ) * NC(3) = 0
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) R, the radius of the circle.
!
!    Input, real ( kind = rk ) PC(3), the center of the circle.
!
!    Input, real ( kind = rk ) NC(3), the normal vector to the circle.
!
!    Input, character ( length = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) nc(dim_num)
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) r
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)'        ) ' '
  write ( *, '(a,g14.6)'  ) '  Radius = ', r
  write ( *, '(a,3g14.6)' ) '  Center = ', pc(1:dim_num)
  write ( *, '(a,3g14.6)' ) '  Normal = ', nc(1:dim_num)

  return
end
subroutine circle_imp_segment_intersect ( r, pc, p0, p1, int_num, p )

!*****************************************************************************80
!
!! circle_imp_segment_intersect: ( imp circle, line segment ) intersection.
!
!  Discussion:
!
!    Points P on an implicit circle in 2D satisfy the equation:
!
!      ( P(1) - PC(1) )^2 + ( P(2) - PC(2) )^2 = R^2
!
!    A line segment in 2D is:
!
!      P(t) = (1-t) * P0 + t * P1
!
!    where P0 and P1 are distinct points and 0 <= t <= 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 July 2020
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) R, the radius of the circle.
!
!    real ( kind = rk ) PC(2), the center of the circle.
!
!    real ( kind = rk ) P0(2), P1(2), the endpoints of the line segment.
!
!  Output:
!
!    integer INT_NUM, the number of intersections.
!    INT_NUM will be 0, 1 or 2.
!
!    real ( kind = rk ) P(2,INT_NUM), the intersecting points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer int_num
  integer int_num_line
  integer j
  real ( kind = rk ) p(2,2)
  real ( kind = rk ) p0(2)
  real ( kind = rk ) p1(2)
  real ( kind = rk ) p2(2)
  real ( kind = rk ) pc(2)
  real ( kind = rk ) p_line(2,2)
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) t
!
!  Find intersections of the full line with the circle.
!
  call circle_imp_line_exp_intersect ( r, pc, p0, p1, int_num_line, p_line )
!
!  Only accept intersection points that lie between P0 and P1.
!
  int_num = 0
  p(1:2,1:2) = 0.0D+00

  do j = 1, int_num_line

    p2(1:2) = p_line(1:2,j)

    call segment_point_coords_2d ( p0, p1, p2, s, t )

    if ( 0.0D+00 <= t .and. t <= 1.0D+00 ) then
      int_num = int_num + 1
      p(1:2,int_num) = p2(1:2)
    end if

  end do
  
  return
end
subroutine circle_imp2exp_2d ( r, pc, p1, p2, p3 )

!*****************************************************************************80
!
!! circle_imp2exp_2d converts a circle from implicit to explicit form in 2D.
!
!  Discussion:
!
!    Points P on an implicit circle in 2D satisfy the equation:
!
!      ( P(1) - PC(1) )^2 + ( P(2) - PC(2) )^2 = R^2
!
!    The explicit form of a circle in 2D is:
!
!      The circle passing through points P1, P2 and P3.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph ORourke,
!    Computational Geometry,
!    Second Edition,
!    Cambridge, 1998,
!    ISBN: 0521649765,
!    LC: QA448.D38.
!
!  Input:
!
!    real ( kind = rk ) R, PC(2), the radius and center of the circle.
!
!  Output:
!
!    real ( kind = rk ) P1(2), P2(2), P3(2), three points on the circle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) p3(dim_num)
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) r
  real ( kind = rk ) theta

  theta = 0.0D+00
  p1(1) = pc(1) + r * cos ( theta )
  p1(2) = pc(2) + r * sin ( theta )

  theta = 2.0D+00 * r8_pi / 3.0D+00
  p2(1) = pc(1) + r * cos ( theta )
  p2(2) = pc(2) + r * sin ( theta )

  theta = 4.0D+00 * r8_pi / 3.0D+00
  p3(1) = pc(1) + r * cos ( theta )
  p3(2) = pc(2) + r * sin ( theta )

  return
end
subroutine circle_llr2imp_2d ( p1, p2, q1, q2, r, pc )

!*****************************************************************************80
!
!! circle_llr2imp_2d converts a circle from LLR to implicit form in 2D.
!
!  Discussion:
!
!    The LLR form of a circle in 2D is:
!
!      The circle of radius R tangent to the lines L1 and L2.
!
!    Points P on an implicit circle in 2D satisfy the equation:
!
!      ( P(1) - PC(1) )^2 + ( P(2) - PC(2) )^2 = R^2
!
!    Let S be the scaled distance of a point on L1 from P1 to P2, 
!    and let N1 be a unit normal vector to L1.  Then a point P that is
!    R units from L1 satisfies:
!
!      P = P1 + s * ( P2 - P1 ) + R * N1.
!
!    Let t be the scaled distance of a point on L2 from Q1 to Q2,
!    and let N2 be a unit normal vector to L2.  Then a point Q that is
!    R units from L2 satisfies:
!
!      Q = Q1 + t * ( Q2 - Q1 ) + R * N2.
!
!    For the center of the circle, then, we have P = Q, that is
!
!      ( P2 - P1 ) * s - ( Q2 - Q1 ) * t = - P1 + Q1 - R * N1 + R * N2 )
!
!    This is a linear system for ( s and t ) from which we can compute
!    the points of tangency, and the center.
!
!    Note that we have four choices for the circle based on the use
!    of plus or minus N1 and plus or minus N2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(2), P2(2), two points on line 1.
!
!    Input, real ( kind = rk ) Q1(2), Q2(2), two points on line 2.
!
!    Input, real ( kind = rk ) R, the radius of the circle.  
!
!    Output, real ( kind = rk ) PC(2,4), the centers of the circles.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) a(2,2)
  real ( kind = rk ) b(2)
  real ( kind = rk ) det
  real ( kind = rk ) n1(dim_num)
  real ( kind = rk ) n2(dim_num)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) pc(dim_num,4)
  real ( kind = rk ) q1(dim_num)
  real ( kind = rk ) q2(dim_num)
  real ( kind = rk ) r
  real ( kind = rk ) x(dim_num)
!
!  Compute the normals N1 and N2.
!
  call line_exp_normal_2d ( p1, p2, n1 )

  call line_exp_normal_2d ( q1, q2, n2 )
!
!  Set the linear system.
!
  a(1:2,1) =   p2(1:2) - p1(1:2)
  a(1:2,2) = - q2(1:2) + q1(1:2)
!
!  Solve the 4 linear systems, using every combination of 
!  signs on the normal vectors.
!
  b(1:2) = - p1(1:2) + q1(1:2) + r * n1(1:2) + r * n2(1:2)

  call r8mat_solve_2d ( a, b, det, x )

  pc(1:2,1) = p1(1:2) + ( p2(1:2) - p1(1:2) ) * x(1) - r * n1(1:2) 

  b(1:2) = - p1(1:2) + q1(1:2) + r * n1(1:2) - r * n2(1:2)

  call r8mat_solve_2d ( a, b, det, x )

  pc(1:2,2) = p1(1:2) + ( p2(1:2) - p1(1:2) ) * x(1) - r * n1(1:2) 

  b(1:2) = - p1(1:2) + q1(1:2) - r * n1(1:2) + r * n2(1:2)

  call r8mat_solve_2d ( a, b, det, x )

  pc(1:2,3) = p1(1:2) + ( p2(1:2) - p1(1:2) ) * x(1) + r * n1(1:2) 

  b(1:2) = - p1(1:2) + q1(1:2) - r * n1(1:2) - r * n2(1:2)

  call r8mat_solve_2d ( a, b, det, x )

  pc(1:2,4) = p1(1:2) + ( p2(1:2) - p1(1:2) ) * x(1) + r * n1(1:2) 

  return
end
subroutine circle_lune_angle_by_height_2d ( r, h, angle )

!*****************************************************************************80
!
!! CIRCLE_LUNE_ANGLE_BY_HEIGHT_2D computes the angle of a circular lune.
!
!  Discussion:
!
!    Draw the chord connecting two points on the circumference of a circle.
!    The region between the chord and the circumference is a "lune".
!    We wish to know the angle subtended by the lune.
!
!    The distance from the center of the circle to the midpoint of the chord
!    is the "height" H of the lune.  It is natural to expect 0 <= H <= R.
!    However, if we allow -R <= H < 0 as well, this allows us to include
!    lunes which involve more than half the circle's area.
!
!    If H < -R or R < H, then no lune is formed, and we return a zero angle.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 January 2018
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) R, the radius of the circle.
!
!    Input, real ( kind = rk ) H, the height of the lune.
!
!    Output, real ( kind = rk ) ANGLE, the angle of the lune.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) angle
  real ( kind = rk ) h
  real ( kind = rk ) r

  if ( -r <= h .and. h <= r ) then
    angle = 2.0D+00 * acos ( h / r );
  else
    angle = 0.0D+00
  end if

  return
end
subroutine circle_lune_area_by_angle_2d ( r, pc, theta1, theta2, area )

!*****************************************************************************80
!
!! CIRCLE_LUNE_AREA_BY_ANGLE_2D returns the area of a circular lune in 2D.
!
!  Discussion:
!
!    A lune is formed by drawing a circular arc, and joining its endpoints.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) R, the radius of the circle.
!
!    Input, real ( kind = rk ) PC(2), the center of the circle.
!
!    Input, real ( kind = rk ) THETA1, THETA2, the angles defining the arc,
!    in radians.  Normally, THETA1 < THETA2.
!
!    Output, real ( kind = rk ) AREA, the area of the lune.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) area
  real ( kind = rk ) area_sector
  real ( kind = rk ) area_triangle
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) r
  real ( kind = rk ) theta1
  real ( kind = rk ) theta2

  call circle_sector_area_2d ( r, pc, theta1, theta2, area_sector )
  call circle_triangle_area_2d ( r, pc, theta1, theta2, area_triangle )

  area = area_sector - area_triangle

  return
end
subroutine circle_lune_area_by_height_2d ( r, h )

!*****************************************************************************80
!
!! CIRCLE_LUNE_AREA_BY_ANGLE_2D returns the area of a circular lune in 2D.
!
!  Discussion:
!
!    Draw the chord connecting two points on the circumference of a circle.
!    The region between the chord and the circumference is a "lune".
!    We wish to know the area of this region.
!
!    The distance from the center of the circle to the midpoint of the chord
!    is the "height" H of the lune.  It is natural to expect 0 <= H <= R.
!    However, if we allow -R <= H < 0 as well, this allows us to include
!    lunes which involve more than half the circle's area.
!
!    If H < -R or R < H, then no lune is formed and we have zero area.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 January 2018
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) R, the radius of the circle.
!
!    Input, real ( kind = rk ) H, the height of the lune.
!
!    Output, real ( kind = rk ) AREA, the area of the lune.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) area
  real ( kind = rk ) h
  real ( kind = rk ) r

  if ( -r <= h .and. h <= r ) then
    area = r ** 2 * acos ( h / r ) - h * sqrt ( r ** 2 - h ** 2 )
  else
    area = 0.0D+00
  end if

  return
end
subroutine circle_lune_centroid_2d ( r, pc, theta1, theta2, centroid )

!*****************************************************************************80
!
!! CIRCLE_LUNE_CENTROID_2D returns the centroid of a circular lune in 2D.
!
!  Discussion:
!
!    A lune is formed by drawing a circular arc, and joining its endpoints.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 October 2004
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
!  Parameters:
!
!    Input, real ( kind = rk ) R, the radius of the circle.
!
!    Input, real ( kind = rk ) PC(2), the center of the circle.
!
!    Input, real ( kind = rk ) THETA1, THETA2, the angles defining the arc,
!    in radians.  Normally, THETA1 < THETA2.
!
!    Output, real ( kind = rk ) CENTROID(2), the coordinates of the centroid
!    of the lune.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) centroid(dim_num)
  real ( kind = rk ) d
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) r
  real ( kind = rk ) theta
  real ( kind = rk ) theta1
  real ( kind = rk ) theta2

  theta = theta2 - theta1

  if ( theta == 0.0D+00 ) then
    d = r
  else
    d = 4.0D+00 * r * ( sin ( 0.5D+00 * theta ) )**3 / &
      ( 3.0D+00 * ( theta - sin ( theta ) ) )
  end if

  centroid(1:2) = (/ pc(1) + d * cos ( theta ), &
                     pc(2) + d * sin ( theta ) /)

  return
end
subroutine circle_lune_height_by_angle_2d ( r, angle, height )

!*****************************************************************************80
!
!! CIRCLE_LUNE_HEIGHT_BY_ANGLE_2D computes the height of a circular lune.
!
!  Discussion:
!
!    Draw the chord connecting two points on the circumference of a circle.
!    The region between the chord and the circumference is a "lune".
!    The lune subtends a given angle between 0 and 2 pi.
!
!    The distance from the center of the circle to the midpoint of the chord
!    is the "height" H of the lune and we wish to determine this value.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 January 2018
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) R, the radius of the circle.
!
!    Input, real ( kind = rk ) ANGLE, the angle subtended by the lune.
!
!    Output, real ( kind = rk ) HEIGHT, the height of the lune
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) angle
  real  ( kind = rk ) height
  real  ( kind = rk ) r

  height = r * cos ( angle / 2.0D+00 )

  return
end
subroutine circle_pppr2imp_3d ( p1, p2, p3, r, pc, normal )

!*****************************************************************************80
!
!! CIRCLE_PPPR2IMP_3D converts a circle from PPPR to implicit form in 3D.
!
!  Discussion:
!
!    The PPPR form of a circle in 3D is:
!
!      The circle of radius R passing through points P1 and P2,
!      and lying in the plane of P1, P2 and P3.
!
!    Points P on an implicit circle in 2D satisfy the equations:
!
!        ( P(1) - PC(1) )^2 + ( P(2) - PC(2) )^2 + ( P(3) - PC(3) )^2 = R^2
!      and 
!        ( P - PC ) dot NORMAL = 0.
!
!    There may be zero, one, or two circles that satisfy the
!    requirements of the PPPR form.
!
!    If there is no such circle, then PC(1:2,1) and PC(1:2,2)
!    are set to the midpoint of (P1,P2).
!
!    If there is one circle, PC(1:2,1) and PC(1:2,2) will be equal.
!
!    If there are two circles, then PC(1:2,1) is the first center,
!    and PC(1:2,2) is the second.
!
!    This calculation is equivalent to finding the intersections of
!    spheres of radius R at points P1 and P2, which lie in the plane
!    defined by P1, P2 and P3.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) P1(3), P2(3), two points on the circle.
!
!    real ( kind = rk ) P3(3), a third point.
!
!    real ( kind = rk ) R, the radius of the circle.
!
!  Output:
!
!    real ( kind = rk ) PC(3,2), the centers of the two circles.
!
!    real ( kind = rk ) NORMAL(3), the normal to the circles.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) dist
  real ( kind = rk ) dot
  real ( kind = rk ) h
  integer j
  real ( kind = rk ) length
  real ( kind = rk ) normal(dim_num)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) p3(dim_num)
  real ( kind = rk ) pc(dim_num,2)
  real ( kind = rk ) r
  real ( kind = rk ) v(dim_num)
!
!  Compute the distance from P1 to P2.
!
  dist = sqrt ( sum ( ( p2(1:dim_num) - p1(1:dim_num) )**2 ) )
!
!  If R is smaller than DIST, we don't have a circle.
!
  if ( 2.0D+00 * r < dist ) then
    do j = 1, 2
      pc(1:dim_num,j) = 0.5D+00 * ( p1(1:dim_num) + p2(1:dim_num) )
    end do
    return
  end if
!
!  H is the distance from the midpoint of (P1,P2) to the center.
!
  h = sqrt ( ( r + 0.5D+00 * dist ) * ( r - 0.5D+00 * dist ) )
!
!  Define a unit direction V that is normal to P2-P1, and lying
!  in the plane (P1,P2,P3).
!
!  To do this, subtract from P3-P1 the component in the direction P2-P1.
!
  v(1:dim_num) = p3(1:dim_num) - p1(1:dim_num)
  dot = dot_product ( v(1:dim_num), p2(1:dim_num) - p1(1:dim_num) )
  dot = dot / dist

  v(1:dim_num) = v(1:dim_num) - dot * ( p2(1:dim_num) - p1(1:dim_num) ) / dist

  length = sqrt ( sum ( v(1:dim_num)**2 ) )

  v(1:dim_num) = v(1:dim_num) / length
!
!  We can go with or against the given normal direction.
!
  pc(1:dim_num,1) = 0.5D+00 * ( p2(1:dim_num) + p1(1:dim_num) ) &
    + h * v(1:dim_num)

  pc(1:dim_num,2) = 0.5D+00 * ( p2(1:dim_num) + p1(1:dim_num) ) &
    - h * v(1:dim_num)

  call plane_exp_normal_3d ( p1, p2, p3, normal )

  return
end
subroutine circle_ppr2imp_2d ( p1, p2, r, pc )

!*****************************************************************************80
!
!! CIRCLE_PPR2IMP_2D converts a circle from PPR to implicit form in 2D.
!
!  Discussion:
!
!    The PPR form of a circle in 2D is:
!
!      The circle of radius R passing through points P1 and P2.
!
!    Points P on an implicit circle in 2D satisfy the equation:
!
!      ( P(1) - PC(1) )^2 + ( P(2) - PC(2) )^2 = R^2
!
!    There may be zero, one, or two circles that satisfy the 
!    requirements of the PPR form.
!
!    If there is no such circle, then PC(1:2,1) and PC(1:2,2)
!    are set to the midpoint of (P1,P2).
!
!    If there is one circle, PC(1:2,1) and PC(1:2,2) will be equal.
!
!    If there are two circles, then PC(1:2,1) is the first center,
!    and PC(1:2,2) is the second.
!
!    This calculation is equivalent to finding the intersections of
!    circles of radius R at points P1 and P2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) P1(2), P2(2), two points on the circle.
!
!    real ( kind = rk ) R, the radius of the circle.  
!
!  Output:
!
!    real ( kind = rk ) PC(2,2), the centers of the two circles.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) dist
  real ( kind = rk ) h
  integer j
  real ( kind = rk ) normal(dim_num)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) pc(dim_num,2)
  real ( kind = rk ) r
!
!  Compute the distance from P1 to P2.
!
  dist = sqrt ( sum ( ( p2(1:dim_num) - p1(1:dim_num) )**2 ) )
!
!  If R is smaller than DIST, we don't have a circle.
!
  if ( 2.0D+00 * r < dist ) then
    do j = 1, 2
      pc(1:dim_num,j) = 0.5D+00 * ( p1(1:dim_num) + p2(1:dim_num) )
    end do
    return
  end if
!
!  H is the distance from the midpoint of (P1,P2) to the center.
!
  h = sqrt ( ( r + 0.5D+00 * dist ) * ( r - 0.5D+00 * dist ) )
!
!  Determine the unit normal direction.
!
  normal(1) =   ( p2(2) - p1(2) ) / dist
  normal(2) = - ( p2(1) - p1(1) ) / dist
!
!  We can go with or against the given normal direction.
!
  pc(1:dim_num,1) = 0.5D+00 * ( p2(1:dim_num) + p1(1:dim_num) ) &
    + h * normal(1:dim_num)

  pc(1:dim_num,2) = 0.5D+00 * ( p2(1:dim_num) + p1(1:dim_num) ) &
    - h * normal(1:dim_num)

  return
end
subroutine circle_sector_area_2d ( r, pc, theta1, theta2, area )

!*****************************************************************************80
!
!! CIRCLE_SECTOR_AREA_2D computes the area of a circular sector in 2D.
!
!  Discussion:
!
!    A circular sector is formed by a circular arc, and the two straight line 
!    segments that join its ends to the center of the circle.
!
!    A circular sector is defined by the two conditions
!
!      ( P(1) - PC(1) )^2 + ( P(2) - PC(2) )^2 = R^2
!
!    and
!
!      Theta1 <= Theta <= Theta2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) R, the radius of the circle.
!
!    real ( kind = rk ) PC(2), the center of the circle.
!
!    real ( kind = rk ) THETA1, THETA2, the two angles defining the
!    sector, in radians.  Normally, THETA1 < THETA2.
!
!  Output:
!
!    real ( kind = rk ) AREA, the area of the circle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) area
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) r
  real ( kind = rk ) theta1
  real ( kind = rk ) theta2

  call r8_fake_use ( pc(1) )

  area = 0.5D+00 * r * r * ( theta2 - theta1 )

  return
end
subroutine circle_sector_centroid_2d ( r, pc, theta1, theta2, centroid )

!*****************************************************************************80
!
!! CIRCLE_SECTOR_CENTROID_2D returns the centroid of a circular sector in 2D.
!
!  Discussion:
!
!    A circular sector is formed by a circular arc, and the two straight line 
!    segments that join its ends to the center of the circle.
!
!    A circular sector is defined by
!
!      ( P(1) - PC(1) )^2 + ( P(2) - PC(2) )^2 = R^2
!
!    and
!
!      Theta1 <= Theta <= Theta2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 October 2004
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
!    real ( kind = rk ) R, the radius of the circle.
!
!    real ( kind = rk ) PC(2), the center of the circle.
!
!    real ( kind = rk ) THETA1, THETA2, the angles defining the arc,
!    in radians.  Normally, THETA1 < THETA2.
!
!  Output:
!
!    real ( kind = rk ) CENTROID(2), the coordinates of the centroid
!    of the sector.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) centroid(dim_num)
  real ( kind = rk ) d
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) r
  real ( kind = rk ) theta
  real ( kind = rk ) theta1
  real ( kind = rk ) theta2

  theta = theta2 - theta1

  if ( theta == 0.0D+00 ) then
    d = 2.0D+00 * r / 3.0D+00
  else
    d = 4.0D+00 * r * sin ( 0.5D+00 * theta ) / &
      ( 3.0D+00 * theta )
  end if

  centroid(1:2) = (/ pc(1) + d * cos ( theta ), &
                     pc(2) + d * sin ( theta ) /)

  return
end
subroutine circle_sector_contains_point_2d ( r, pc, theta1, theta2, &
  p, inside )

!*****************************************************************************80
!
!! CIRCLE_SECTOR_CONTAINS_POINT_2D : is a point inside a circular sector?
!
!  Discussion:
!
!    A circular sector is formed by a circular arc, and the two straight line 
!    segments that join its ends to the center of the circle.
!
!    A circular sector is defined by
!
!      ( P(1) - PC(1) )^2 + ( P(2) - PC(2) )^2 = R^2
!
!    and
!
!      Theta1 <= Theta <= Theta2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) R, the radius of the circle.
!
!    real ( kind = rk ) PC(2), the center of the circle.
!
!    real ( kind = rk ) THETA1, THETA2, the angles defining the arc,
!    in radians.  Normally, THETA1 < THETA2.
!
!    real ( kind = rk ) P(2), the point to be checked.
!
!  Output:
!
!    logical INSIDE, is TRUE if the point is inside or 
!    on the circular sector.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  logical inside
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) r
  real ( kind = rk ) r8_atan
  real ( kind = rk ) r8_modp
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) theta
  real ( kind = rk ) theta1
  real ( kind = rk ) theta2

  inside = .false.
!
!  Is the point inside the (full) circle?
!
  if ( ( p(1) - pc(1) ) * ( p(1) - pc(1) ) &
     + ( p(2) - pc(2) ) * ( p(2) - pc(2) ) <= r * r ) then
!
!  Is the point's angle within the arc's range?
!  Try to force the angles to lie between 0 and 2 * PI.
!
    theta = r8_atan ( p(2) - pc(2), p(1) - pc(1) )

    if ( r8_modp ( theta  - theta1,  2.0D+00 * r8_pi ) <= &
         r8_modp ( theta2 - theta1,  2.0D+00 * r8_pi ) ) then

      inside = .true.

    end if

  end if

  return
end
subroutine circle_sector_print_2d ( r, pc, theta1, theta2 )

!*****************************************************************************80
!
!! CIRCLE_SECTOR_PRINT_2D prints a circular sector in 2D.
!
!  Discussion:
!
!    A circular sector is formed by a circular arc, and the two straight line 
!    segments that join its ends to the center of the circle.
!
!    A circular sector is defined by
!
!      ( P(1) - PC(1) )^2 + ( P(2) - PC(2) )^2 = R^2
!
!    and
!
!      Theta1 <= Theta <= Theta2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) R, the radius of the circle.
!
!    real ( kind = rk ) PC(2), the center of the circle.
!
!    real ( kind = rk ) THETA1, THETA2, the angles defining the arc,
!    in radians.  Normally, THETA1 < THETA2.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) r
  real ( kind = rk ) theta1
  real ( kind = rk ) theta2

  write ( *, '(a)'        ) ' '
  write ( *, '(a)'        ) '  Circular sector definition:'
  write ( *, '(a)'        ) ' '
  write ( *, '(a,g14.6)'  ) '    Radius = ', r
  write ( *, '(a,2g14.6)' ) '    Center = ', pc(1:2)
  write ( *, '(a,2g14.6)' ) '    Theta  = ', theta1, theta2

  return
end
subroutine circle_triangle_area_2d ( r, pc, theta1, theta2, area )

!*****************************************************************************80
!
!! circle_triangle_area_2d returns the area of a circle triangle in 2D.
!
!  Discussion:
!
!    A circle triangle is formed by drawing a circular arc, and considering
!    the triangle formed by the endpoints of the arc plus the center of
!    the circle.
!
!    Note that for angles greater than PI, the triangle will actually
!    have NEGATIVE area.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) R, the radius of the circle.
!
!    real ( kind = rk ) PC(2), the center of the circle.
!
!    real ( kind = rk ) THETA1, THETA2, the angles defining the arc,
!    in radians.  Normally, THETA1 < THETA2.
!
!  Output:
!
!    real ( kind = rk ) AREA, the (signed) area of the triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) area
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) r
  real ( kind = rk ) theta1
  real ( kind = rk ) theta2

  call r8_fake_use ( pc(1) )

  area = 0.5D+00 * r * r * sin ( theta2 - theta1 )

  return
end
subroutine circle_triple_angles_2d ( r1, r2, r3, angle1, angle2, angle3 )

!*****************************************************************************80
!
!! circle_triple_angles_2d returns the angles formed by three cotangent circles in 2D.
!
!  Discussion:
!
!    A circle triple is a set of three tangent circles.  We assume
!    that no circle is contained in another.
!
!    We consider the triangle formed by joining the centers of the circles.
!
!    We report the angles in that triangle.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 September 2020
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Kenneth Stephenson,
!    Circle Packing, The Theory of Discrete Analytic Functions,
!    Cambridge, 2005.
!
!  Input:
!
!    real ( kind = rk ) R1, R2, R3, the radii of the circles.
!
!  Output:
!
!    real ( kind = rk ) ANGLE1, ANGLE2, ANGLE3, the angles in the triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) angle1
  real ( kind = rk ) angle2
  real ( kind = rk ) angle3
  real ( kind = rk ) r1
  real ( kind = rk ) r2
  real ( kind = rk ) r3

  angle1 = acos &
  ( &
    ( &
      ( r1 + r2 )**2 + ( r1 + r3 )**2 - ( r2 + r3 )**2 &
    ) / &
    ( &
      2.0D+00 * ( r1 + r2 ) * ( r1 + r3 ) &
    ) &
  ) 

  angle2 = acos &
  ( &
    ( &
      ( r2 + r3 )**2 + ( r2 + r1 )**2 - ( r3 + r1 )**2 &
    ) / &
    ( &
      2.0D+00 * ( r2 + r3 ) * ( r2 + r1 ) &
    ) &
  )

  angle3 = acos &
  ( &
    ( &
      ( r3 + r1 )**2 + ( r3 + r2 )**2 - ( r1 + r2 )**2 &
    ) / &
    ( &
      2.0D+00 * ( r3 + r1 ) * ( r3 + r2 ) & 
    ) &
  ) 

  return
end
subroutine circles_intersect_points_2d ( r1, pc1, r2, pc2, int_num, p )

!*****************************************************************************80
!
!! CIRCLES_INTERSECT_POINTS_2D: intersection points of two circles in 2D.
!
!  Discussion:
!
!    Two circles can intersect in 0, 1, 2 or infinitely many points.
!
!    The 0 and 2 intersection cases are numerically robust; the 1 and
!    infinite intersection cases are numerically fragile.  The routine
!    uses a tolerance to try to detect the 1 and infinite cases.
!
!    Points P on an implicit circle in 2D satisfy the equation:
!
!      ( P(1) - PC(1) )^2 + ( P(2) - PC(2) )^2 = R^2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) R1, the radius of the first circle.
!
!    real ( kind = rk ) PC1(2), the center of the first circle.
!
!    real ( kind = rk ) R2, the radius of the second circle.
!
!    real ( kind = rk ) PC2(2), the center of the second circle.
!
!  Output:
!
!    integer INT_NUM, the number of intersecting points 
!    found.  INT_NUM will be 0, 1, 2 or 3.  3 indicates that there are an 
!    infinite number of intersection points.
!
!    real ( kind = rk ) P(2,2), if INT_NUM is 1 or 2,
!    the coordinates of the intersecting points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) distsq
  integer int_num
  real ( kind = rk ) p(dim_num,2)
  real ( kind = rk ) pc1(dim_num)
  real ( kind = rk ) pc2(dim_num)
  real ( kind = rk ) r1
  real ( kind = rk ) r2
  real ( kind = rk ) root
  real ( kind = rk ) sc1
  real ( kind = rk ) sc2
  real ( kind = rk ) t1
  real ( kind = rk ) t2
  real ( kind = rk ) tol

  tol = epsilon ( tol )

  p(1:dim_num,1:2) = 0.0D+00
!
!  Take care of the case in which the circles have the same center.
!
  t1 = ( abs ( pc1(1) - pc2(1) ) &
       + abs ( pc1(2) - pc2(2) ) ) / 2.0D+00

  t2 = ( abs ( pc1(1) ) + abs ( pc2(1) ) &
       + abs ( pc1(2) ) + abs ( pc2(2) ) + 1.0D+00 ) / 5.0D+00

  if ( t1 <= tol * t2 ) then

    t1 = abs ( r1 - r2 )
    t2 = ( abs ( r1 ) + abs ( r2 ) + 1.0D+00 ) / 3.0D+00

    if ( t1 <= tol * t2 ) then
      int_num = 3
    else
      int_num = 0
    end if

    return

  end if

  distsq = ( pc1(1) - pc2(1) )**2 + ( pc1(2) - pc2(2) )**2

  root = 2.0D+00 * ( r1**2 + r2**2 ) * distsq - distsq**2 &
    - ( r1 - r2 )**2 * ( r1 + r2 )**2

  if ( root < -tol ) then
    int_num = 0
    return
  end if

  sc1 = ( distsq - ( r2**2 - r1**2 ) ) / distsq

  if ( root < tol ) then
    int_num = 1
    p(1:dim_num,1) = pc1(1:dim_num) &
      + 0.5D+00 * sc1 * ( pc2(1:dim_num) - pc1(1:dim_num) )
    return
  end if

  sc2 = sqrt ( root ) / distsq

  int_num = 2

  p(1,1) = pc1(1) + 0.5D+00 * sc1 * ( pc2(1) - pc1(1) ) &
                  - 0.5D+00 * sc2 * ( pc2(2) - pc1(2) )
  p(2,1) = pc1(2) + 0.5D+00 * sc1 * ( pc2(2) - pc1(2) ) &
                  + 0.5D+00 * sc2 * ( pc2(1) - pc1(1) )

  p(1,2) = pc1(1) + 0.5D+00 * sc1 * ( pc2(1) - pc1(1) ) &
                  + 0.5D+00 * sc2 * ( pc2(2) - pc1(2) )
  p(2,2) = pc1(2) + 0.5D+00 * sc1 * ( pc2(2) - pc1(2) ) &
                  - 0.5D+00 * sc2 * ( pc2(1) - pc1(1) )

  return
end
subroutine combin2 ( n, k, icnk )

!*****************************************************************************80
!
!! COMBIN2 computes the binomial coefficient C(N,K).
!
!  Discussion:
!
!    The value is calculated in such a way as to avoid overflow and
!    roundoff.  The calculation is done in integer arithmetic.
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
!    17 January 1999
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
!    April, 1963.
!
!  Input:
!
!    integer N, K, are the values of N and K.
!
!  Output:
!
!    integer ICNK, the number of combinations of N
!    things taken K at a time.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer icnk
  integer k
  integer mn
  integer mx
  integer n

  mn = min ( k, n - k )

  if ( mn < 0 ) then

    icnk = 0

  else if ( mn == 0 ) then

    icnk = 1

  else

    mx = max ( k, n - k )
    icnk = mx + 1

    do i = 2, mn
      icnk = ( icnk * ( mx + i ) ) / i
    end do

  end if

  return
end
subroutine cone_area_3d ( h, r, area )

!*****************************************************************************80
!
!! CONE_AREA_3D computes the surface area of a right circular cone in 3D.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) H, R, the height of the cone, and the radius
!    of the circle that forms the base of the cone.
!
!    Output, real ( kind = rk ) AREA, the surface area of the cone.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) area
  real ( kind = rk ) h
  real ( kind = rk ) r
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00

  area = r8_pi * r * sqrt ( h * h + r * r )

  return
end
subroutine cone_centroid_3d ( r, pc, pt, centroid )

!*****************************************************************************80
!
!! CONE_CENTROID_3D returns the centroid of a cone in 3D.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 January 2005
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
!  Parameters:
!
!    Input, real ( kind = rk ) R, the radius of the circle at the base of
!    the cone.
!
!    Input, real ( kind = rk ) PC(3), the center of the circle.
!
!    Input, real ( kind = rk ) PT(3), the coordinates of the tip of the cone.
!
!    Output, real ( kind = rk ) CENTROID(3), the coordinates of the centroid
!    of the cone.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) centroid(dim_num)
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) pt(dim_num)
  real ( kind = rk ) r

  call r8_fake_use ( r )

  centroid(1:dim_num) = 0.75D+00 * pc(1:dim_num) + 0.25D+00 * pt(1:dim_num)

  return
end
subroutine cone_volume_3d ( h, r, volume )

!*****************************************************************************80
!
!! CONE_VOLUME_3D computes the volume of a right circular cone in 3D.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) H, R, the height of the cone, and the radius
!    of the circle that forms the base of the cone.
!
!    Output, real ( kind = rk ) VOLUME, the volume of the cone.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) h
  real ( kind = rk ) r
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) volume

  volume = r8_pi * r * r * h / 3.0D+00

  return
end
subroutine conv3d ( axis, theta, n, cor3, cor2 )

!*****************************************************************************80
!
!! CONV3D converts 3D data to a 2D projection.
!
!  Discussion:
!
!    A "presentation angle" THETA is used to project the 3D point
!    (X3D, Y3D, Z3D) to the 2D projection (XVAL,YVAL).
!
!    If AXIS = 'X':
!
!      X2D = Y3D - sin ( THETA ) * X3D
!      Y2D = Z3D - sin ( THETA ) * X3D
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    24 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character AXIS, the coordinate axis to be projected.
!    AXIS should be 'X', 'Y', or 'Z'.
!
!    Input, real ( kind = rk ) THETA, the presentation angle in degrees.
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) COR3(3,N), the 3D points.
!
!    Output, real ( kind = rk ) COR2(2,N), the 2D projections.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  character axis
  real ( kind = rk ) cor2(2,n)
  real ( kind = rk ) cor3(3,n)
  real ( kind = rk ) degrees_to_radians
  real ( kind = rk ) stheta
  real ( kind = rk ) theta

  stheta = sin ( degrees_to_radians ( theta ) )

  if ( axis == 'X' .or. axis == 'x' ) then

    cor2(1,1:n) = cor3(2,1:n) - stheta * cor3(1,1:n)
    cor2(2,1:n) = cor3(3,1:n) - stheta * cor3(1,1:n)

  else if ( axis == 'Y' .or. axis == 'y' ) then

    cor2(1,1:n) = cor3(1,1:n) - stheta * cor3(2,1:n)
    cor2(2,1:n) = cor3(3,1:n) - stheta * cor3(2,1:n)

  else if ( axis == 'Z' .or. axis == 'z' ) then

    cor2(1,1:n) = cor3(1,1:n) - stheta * cor3(3,1:n)
    cor2(2,1:n) = cor3(2,1:n) - stheta * cor3(3,1:n)

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CONV3D - Fatal error!'
    write ( *, '(a)' ) '  Illegal coordinate index = "' // axis // '".'
    stop 1

  end if

  return
end
function cot_rad ( angle_rad )

!*****************************************************************************80
!
!! COT_RAD returns the cotangent of an angle.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) ANGLE_RAD, the angle, in radians.
!
!    Output, real ( kind = rk ) COT_RAD, the cotangent of the angle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) angle_rad
  real ( kind = rk ) cot_rad

  cot_rad  = cos ( angle_rad ) / sin ( angle_rad )

  return
end
subroutine cube_shape_3d ( point_num, face_num, face_order_max, &
  point_coord, face_order, face_point )

!*****************************************************************************80
!
!! CUBE_SHAPE_3D describes a cube in 3D.
!
!  Discussion:
!
!    The vertices lie on the unit sphere.
!
!    The dual of the cube is the octahedron.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer POINT_NUM, the number of points.
!
!    Input, integer FACE_NUM, the number of faces.
!
!    Input, integer FACE_ORDER_MAX, the maximum number of vertices
!    in a face.
!
!    Output, real ( kind = rk ) POINT_COORD(3,POINT_NUM),
!    the vertices.
!
!    Output, integer FACE_ORDER(FACE_NUM), the number of vertices
!    per face.
!
!    Output, integer FACE_POINT(FACE_ORDER_MAX,FACE_NUM); 
!    FACE_POINT(I,J) contains the index of the I-th point in the J-th face.  The
!    points are listed in the counter clockwise direction defined
!    by the outward normal at the face.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer face_num
  integer face_order_max
  integer, parameter :: dim_num = 3
  integer point_num

  real ( kind = rk ) a
  integer face_order(face_num)
  integer face_point(face_order_max,face_num)
  real ( kind = rk ) point_coord(dim_num,point_num)
!
!  Set point coordinates.
!
  a = sqrt ( 1.0D+00 / 3.0D+00 )

  point_coord(1:dim_num,1:point_num) = reshape ( (/ &
     -a, -a, -a, &
      a, -a, -a, &
      a,  a, -a, &
     -a,  a, -a, &
     -a, -a,  a, &
      a, -a,  a, &
      a,  a,  a, &
     -a,  a,  a /), (/ dim_num, point_num /) )
!
!  Set the face orders.
!
  face_order(1:face_num) = (/ &
    4, 4, 4, 4, 4, 4 /)
!
!  Set the faces.
!
  face_point(1:face_order_max,1:face_num) = reshape ( (/ &
     1, 4, 3, 2, &
     1, 2, 6, 5, &
     2, 3, 7, 6, &
     3, 4, 8, 7, &
     1, 5, 8, 4, &
     5, 6, 7, 8 /), (/ face_order_max, face_num /) )

  return
end
subroutine cube_size_3d ( point_num, edge_num, face_num, face_order_max )

!*****************************************************************************80
!
!! CUBE_SIZE_3D gives "sizes" for a cube in 3D.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer POINT_NUM, the number of points.
!
!    Output, integer EDGE_NUM, the number of edges.
!
!    Output, integer FACE_NUM, the number of faces.
!
!    Output, integer FACE_ORDER_MAX, the maximum order of any face.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer edge_num
  integer face_num
  integer face_order_max
  integer point_num

  point_num = 8
  edge_num = 12
  face_num = 6
  face_order_max = 4

  return
end
function cube01_volume ( )

!*****************************************************************************80
!
!! CUBE01_VOLUME returns the volume of the unit cube in 3D.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 January 2018
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) CUBE01_VOLUME, the volume.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) cube01_volume

  cube01_volume = 1.0D+00

  return
end
subroutine cylinder_point_dist_3d ( p1, p2, r, p, distance )

!*****************************************************************************80
!
!! CYLINDER_POINT_DIST_3D: distance from a cylinder to a point in 3D.
!
!  Discussion:
!
!    We are computing the distance to the SURFACE of the cylinder.
!
!    The surface of a (right) (finite) cylinder in 3D is defined by an axis,
!    which is the line segment from point P1 to P2, and a radius R.  The points 
!    on the surface of the cylinder are:
!    * points at a distance R from the line through P1 and P2, and whose nearest
!      point on the line through P1 and P2 is strictly between P1 and P2, 
!    PLUS
!    * points at a distance less than or equal to R from the line through P1
!      and P2, whose nearest point on the line through P1 and P2 is either 
!      P1 or P2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(3), P2(3), the first and last points
!    on the axis line of the cylinder.
!
!    Input, real ( kind = rk ) R, the radius of the cylinder.
!
!    Input, real ( kind = rk ) P(3), the point.
!
!    Output, real ( kind = rk ) DISTANCE, the distance from the point 
!    to the cylinder.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) axis(dim_num)
  real ( kind = rk ) axis_length
  real ( kind = rk ) distance
  real ( kind = rk ) r8vec_norm
  real ( kind = rk ) off_axis_component
  real ( kind = rk ) p(3)
  real ( kind = rk ) p_dot_axis
  real ( kind = rk ) p_length
  real ( kind = rk ) p1(3)
  real ( kind = rk ) p2(3)
  real ( kind = rk ) r

  axis(1:dim_num) = p2(1:dim_num) - p1(1:dim_num)
  axis_length = r8vec_norm ( dim_num, axis )

  if ( axis_length == 0.0D+00 ) then
    distance = -huge ( distance )
    return
  end if

  axis(1:dim_num) = axis(1:dim_num) / axis_length

  p_dot_axis = dot_product ( p(1:dim_num) - p1(1:dim_num), axis )
!
!  Case 1: Below bottom cap.
!
  if ( p_dot_axis <= 0.0D+00 ) then

    call disk_point_dist_3d ( p1, r, axis, p, distance )
!
!  Case 2: between cylinder planes.
!
  else if ( p_dot_axis <= axis_length ) then

    p_length = r8vec_norm ( dim_num, p(1:dim_num) - p1(1:dim_num) )
    off_axis_component = sqrt ( p_length**2 - p_dot_axis**2 )

    distance = abs ( off_axis_component - r )

    if ( off_axis_component < r ) then
      distance = min ( distance, axis_length - p_dot_axis )
      distance = min ( distance, p_dot_axis )
    end if
!
!  Case 3: Above the top cap.
!  
  else if ( axis_length < p_dot_axis ) then

    call disk_point_dist_3d ( p2, r, axis, p, distance )

  end if

  return
end
subroutine cylinder_point_dist_signed_3d ( p1, p2, r, p, distance )

!*****************************************************************************80
!
!! CYLINDER_POINT_DIST_SIGNED_3D: signed distance from cylinder to point in 3D.
!
!  Discussion:
!
!    We are computing the signed distance to the SURFACE of the cylinder.
!
!    The surface of a (right) (finite) cylinder in 3D is defined by an axis,
!    which is the line segment from point P1 to P2, and a radius R.  The points 
!    on the surface of the cylinder are:
!    * points at a distance R from the line through P1 and P2, and whose nearest
!      point on the line through P1 and P2 is strictly between P1 and P2, 
!    PLUS
!    * points at a distance less than or equal to R from the line through P1
!      and P2, whose nearest point on the line through P1 and P2 is either 
!      P1 or P2.
!
!    Points inside the surface have a negative distance.
!    Points on the surface have a zero distance.
!    Points outside the surface have a positive distance.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(3), P2(3), the first and last points
!    on the axis line of the cylinder.
!
!    Input, real ( kind = rk ) R, the radius of the cylinder.
!
!    Input, real ( kind = rk ) P(3), the point.
!
!    Output, real ( kind = rk ) DISTANCE, the signed distance from the point 
!    to the cylinder.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) axis(dim_num)
  real ( kind = rk ) axis_length
  real ( kind = rk ) distance
  real ( kind = rk ) r8vec_norm
  real ( kind = rk ) off_axis_component
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) p_dot_axis
  real ( kind = rk ) p_length
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) r

  axis(1:dim_num) = p2(1:dim_num) - p1(1:dim_num)
  axis_length = r8vec_norm ( dim_num, axis )

  if ( axis_length == 0.0D+00 ) then
    distance = -huge ( distance )
    return
  end if

  axis(1:dim_num) = axis(1:dim_num) / axis_length

  p_dot_axis = dot_product ( p(1:dim_num) - p1(1:dim_num), axis )
!
!  Case 1: Below bottom cap.
!
  if ( p_dot_axis <= 0.0D+00 ) then

    call disk_point_dist_3d ( p1, r, axis, p, distance )
!
!  Case 2: between cylinder planes.
!
  else if ( p_dot_axis <= axis_length ) then

    p_length = r8vec_norm ( dim_num, p(1:dim_num) - p1(1:dim_num) )
    off_axis_component = sqrt ( p_length**2 - p_dot_axis**2 )

    distance = off_axis_component - r 

    if ( distance < 0.0D+00 ) then
      distance = max ( distance, p_dot_axis - axis_length )
      distance = max ( distance, -p_dot_axis )
    end if
!
!  Case 3: Above the top cap.
!  
  else if ( axis_length < p_dot_axis ) then

    call disk_point_dist_3d ( p2, r, axis, p, distance )

  end if

  return
end
subroutine cylinder_point_inside_3d ( p1, p2, r, p, inside )

!*****************************************************************************80
!
!! CYLINDER_POINT_INSIDE_3D determines if a cylinder contains a point in 3D.
!
!  Discussion:
!
!    The surface and interior of a (right) (finite) cylinder in 3D is defined 
!    by an axis, which is the line segment from point P1 to P2, and a 
!    radius R.  The points contained in the volume include:
!    * points at a distance less than or equal to R from the line through P1
!      and P2, whose nearest point on the line through P1 and P2 is, in fact,
!      P1, P2, or any point between them.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(3), P2(3), the first and last points
!    on the axis line of the cylinder.
!
!    Input, real ( kind = rk ) R, the radius of the cylinder.
!
!    Input, real ( kind = rk ) P(3), the point.
!
!    Output, logical INSIDE, is TRUE if the point is 
!    inside the cylinder.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) axis(dim_num)
  real ( kind = rk ) axis_length
  logical inside
  real ( kind = rk ) off_axis_component
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) p_dot_axis
  real ( kind = rk ) p_length
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) r
  real ( kind = rk ) r8vec_norm

  axis(1:dim_num) = p2(1:dim_num) - p1(1:dim_num)
  axis_length = r8vec_norm ( dim_num, axis )

  if ( axis_length == 0.0D+00 ) then
    inside = .false.
    return
  end if

  axis(1:dim_num) = axis(1:dim_num) / axis_length

  p_dot_axis = dot_product ( p(1:dim_num) - p1(1:dim_num), axis )
!
!  If the point lies below or above the "caps" of the cylinder, we're done.
!
  if ( p_dot_axis < 0.0D+00 .or. axis_length < p_dot_axis ) then

    inside = .false.
!
!  Otherwise, determine the distance from P to the axis.
!
  else

    p_length = r8vec_norm ( dim_num, p(1:dim_num) - p1(1:dim_num) )

    off_axis_component = sqrt ( p_length**2 - p_dot_axis**2 )

    if ( off_axis_component <= r ) then
      inside = .true.
    else
      inside = .false.
    end if

  end if

  return
end
subroutine cylinder_point_near_3d ( p1, p2, r, p, pn )

!*****************************************************************************80
!
!! CYLINDER_POINT_NEAR_3D: nearest point on a cylinder to a point in 3D.
!
!  Discussion:
!
!    We are computing the nearest point on the SURFACE of the cylinder.
!
!    The surface of a (right) (finite) cylinder in 3D is defined by an axis,
!    which is the line segment from point P1 to P2, and a radius R.  The points 
!    on the surface of the cylinder are:
!    * points at a distance R from the line through P1 and P2, and whose nearest
!      point on the line through P1 and P2 is strictly between P1 and P2, 
!    PLUS
!    * points at a distance less than or equal to R from the line through P1
!      and P2, whose nearest point on the line through P1 and P2 is either 
!      P1 or P2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(3), P2(3), the first and last points
!    on the axis line of the cylinder.
!
!    Input, real ( kind = rk ) R, the radius of the cylinder.
!
!    Input, real ( kind = rk ) P(3), the point.
!
!    Output, real ( kind = rk ) PN(3), the nearest point on the cylinder.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) axial_component
  real ( kind = rk ) axis(dim_num)
  real ( kind = rk ) axis_length
  real ( kind = rk ) distance
  real ( kind = rk ) r8vec_norm
  real ( kind = rk ) off_axis(dim_num)
  real ( kind = rk ) off_axis_component
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) pn(dim_num)
  real ( kind = rk ) r

  axis(1:dim_num) = p2(1:dim_num) - p1(1:dim_num)
  axis_length = r8vec_norm ( dim_num, axis )
  axis(1:dim_num) = axis(1:dim_num) / axis_length

  axial_component = dot_product ( p(1:dim_num) - p1(1:dim_num), axis )

  off_axis(1:dim_num) = p(1:dim_num) - p1(1:dim_num) &
    - axial_component * axis(1:dim_num)

  off_axis_component = r8vec_norm ( dim_num, off_axis )
!
!  Case 1: Below bottom cap.
!
  if ( axial_component <= 0.0D+00 ) then

    if ( off_axis_component <= r ) then
      pn(1:dim_num) = p1(1:dim_num) + off_axis(1:dim_num)
    else
      pn(1:dim_num) = p1(1:dim_num) &
        + ( r / off_axis_component ) * off_axis(1:dim_num)
    end if
!
!  Case 2: between cylinder planes.
!
  else if ( axial_component <= axis_length ) then

    if ( off_axis_component == 0.0D+00 ) then

      call r8vec_any_normal ( dim_num, axis, off_axis )
      
      pn(1:dim_num) = p(1:dim_num) + r * off_axis(1:dim_num)

    else

      distance = abs ( off_axis_component - r )

      pn(1:dim_num) = p1(1:dim_num) + axial_component * axis(1:dim_num) &
        + ( r / off_axis_component ) * off_axis(1:dim_num)

      if ( off_axis_component < r ) then

        if ( axis_length - axial_component < distance ) then
          distance = axis_length - axial_component
          pn(1:dim_num) = p2(1:dim_num) + off_axis(1:dim_num)
        end if

        if ( axial_component < distance ) then
          distance = axial_component
          pn(1:dim_num) = p1(1:dim_num) + off_axis(1:dim_num)
        end if

      end if

    end if
!
!  Case 3: Above the top cap.
!  
  else if ( axis_length < axial_component ) then

    if ( off_axis_component <= r ) then
      pn(1:dim_num) = p2(1:dim_num) + off_axis(1:dim_num)
    else
      pn(1:dim_num) = p2(1:dim_num) &
        + ( r / off_axis_component ) * off_axis(1:dim_num)
    end if

  end if

  return
end
subroutine cylinder_sample_3d ( p1, p2, r, n, p )

!*****************************************************************************80
!
!! CYLINDER_SAMPLE_3D samples a cylinder in 3D.
!
!  Discussion:
!
!    We are sampling the interior of a right finite cylinder in 3D.
!
!    The interior of a (right) (finite) cylinder in 3D is defined by an axis,
!    which is the line segment from point P1 to P2, and a radius R.  The points 
!    on or inside the cylinder are:
!    * points whose distance from the line through P1 and P2 is less than
!      or equal to R, and whose nearest point on the line through P1 and P2
!      lies (nonstrictly) between P1 and P2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(3), P2(3), the first and last points
!    on the axis line of the cylinder.
!
!    Input, real ( kind = rk ) R, the radius of the cylinder.
!
!    Input, integer N, the number of sample points to compute.
!
!    Input, real ( kind = rk ) P(3,N), the sample points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer n

  real ( kind = rk ) axis(dim_num)
  real ( kind = rk ) axis_length
  real ( kind = rk ) r8vec_norm
  integer i
  real ( kind = rk ) p(dim_num,n)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) r
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) radius(n)
  real ( kind = rk ) theta(n)
  real ( kind = rk ) v2(dim_num)
  real ( kind = rk ) v3(dim_num)
  real ( kind = rk ) z(n)
!
!  Compute the axis vector.
!
  axis(1:dim_num) = p2(1:dim_num) - p1(1:dim_num)
  axis_length = r8vec_norm ( dim_num, axis )
  axis(1:dim_num) = axis(1:dim_num) / axis_length
!
!  Compute vectors V2 and V3 that form an orthogonal triple with AXIS.
!
  call plane_normal_basis_3d ( p1, axis, v2, v3 )
!
!  Assemble the randomized information.
!
  call random_number ( harvest = radius(1:n) )
  radius(1:n) = r * sqrt ( radius(1:n) )

  call random_number ( harvest = theta(1:n) )
  theta(1:n) = 2.0D+00 * r8_pi * theta(1:n)

  call random_number ( harvest = z(1:n) )
  z(1:n) = axis_length * z(1:n)

  do i = 1, dim_num

    p(i,1:n) =                                       p1(i)   &
              + z(1:n)                             * axis(i) &
              + radius(1:n) * cos ( theta(1:n) )   * v2(i)   &
              + radius(1:n) * sin ( theta(1:n) )   * v3(i)

  end do


  return
end
subroutine cylinder_volume_3d ( p1, p2, r, volume )

!*****************************************************************************80
!
!! CYLINDER_VOLUME_3D determines the volume of a cylinder in 3D.
!
!  Discussion:
!
!    The surface and interior of a (right) (finite) cylinder in 3D is defined 
!    by an axis, which is the line segment from point P1 to P2, and a radius R.
!    The points contained in the volume include:
!    * points at a distance less than or equal to R from the line through P1 
!      and P2, whose nearest point on the line through P1 and P2 is, in fact,
!      P1, P2, or any point between them.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(3), P2(3), the first and last points
!    on the axis line of the cylinder.
!
!    Input, real ( kind = rk ) R, the radius of the cylinder.
!
!    Output, real ( kind = rk ) VOLUME, the volume of the cylinder.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) h
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) r
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) volume

  h = sqrt ( sum ( ( p1(1:dim_num) - p2(1:dim_num) )**2 ) )

  volume = r8_pi * r * r * h

  return
end
function degrees_to_radians ( angle_degrees )

!*****************************************************************************80
!
!! DEGREES_TO_RADIANS converts an angle from degrees to radians.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) angle_degrees, an angle in degrees.
!
!  Output:
!
!    real ( kind = rk ) DEGREES_TO_RADIANS, the equivalent angle
!    in radians.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) angle_degrees
  real ( kind = rk ) degrees_to_radians
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00

  degrees_to_radians = ( angle_degrees / 180.0D+00 ) * r8_pi

  return
end
subroutine direction_pert_3d ( sigma, vbase, vran )

!*****************************************************************************80
!
!! DIRECTION_PERT_3D randomly perturbs a direction vector in 3D.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) SIGMA, determines the strength of the
!    perturbation.
!    SIGMA <= 0 results in a completely random direction.
!    1 <= SIGMA results in VBASE.
!    0 < SIGMA < 1 results in a perturbation from VBASE, which is
!    large when SIGMA is near 0, and small when SIGMA is near 1.
!
!    Input, real ( kind = rk ) VBASE(3), the base direction vector, which
!    should have unit norm.
!
!    Output, real ( kind = rk ) VRAN(3), the perturbed vector, which will
!    have unit norm.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) dphi
  real ( kind = rk ) phi
  real ( kind = rk ) psi
  real ( kind = rk ) r
  real ( kind = rk ) r8_acos
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) sigma
  real ( kind = rk ) theta
  real ( kind = rk ) v(dim_num)
  real ( kind = rk ) vbase(dim_num)
  real ( kind = rk ) vdot
  real ( kind = rk ) vran(dim_num)
  real ( kind = rk ) x
!
!  1 <= SIGMA, just use the base vector.
!
  if ( 1.0D+00 <= sigma ) then

    vran(1:dim_num) = vbase(1:dim_num)

  else if ( sigma <= 0.0D+00 ) then

    call random_number ( harvest = vdot )
    vdot = 2.0D+00 * vdot - 1.0D+00

    phi = r8_acos ( vdot )

    call random_number ( harvest = theta )
    theta = 2.0D+00 * r8_pi * theta

    vran(1) = cos ( theta ) * sin ( phi )
    vran(2) = sin ( theta ) * sin ( phi )
    vran(3) = cos ( phi )

  else

    phi = r8_acos ( vbase(3) )
    theta = atan2 ( vbase(2), vbase(1) )
!
!  Pick VDOT, which must be between -1 and 1.  This represents
!  the dot product of the perturbed vector with the base vector.
!
!  Return a uniformly random value between 0 and 1.
!  The operations we perform on this quantity tend to bias it
!  out towards 1, as SIGMA grows from 0 to 1.
!
!  VDOT, in turn, is a value between -1 and 1, which, for large
!  SIGMA, we want biased towards 1.
!
    call random_number ( harvest = r )
    x = exp ( ( 1.0D+00 - sigma ) * log ( r ) )
    dphi = r8_acos ( 2.0D+00 * x - 1.0D+00 )
!
!  Now we know enough to write down a vector that is rotated DPHI
!  from the base vector.
!
    v(1) = cos ( theta ) * sin ( phi + dphi )
    v(2) = sin ( theta ) * sin ( phi + dphi )
    v(3) = cos ( phi + dphi )
!
!  Pick a uniformly random rotation between 0 and 2 Pi around the
!  axis of the base vector.
!
    call random_number ( harvest = psi )
    psi = 2.0D+00 * r8_pi * psi
!
!  Carry out the rotation.
!
    call rotation_axis_vector_3d ( vbase, psi, v, vran )

  end if

  return
end
subroutine direction_uniform_2d ( vran )

!*****************************************************************************80
!
!! DIRECTION_UNIFORM_2D picks a random direction vector in 2D.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) VRAN(2), the random direction vector, with
!    unit norm.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) theta
  real ( kind = rk ) vran(dim_num)

  call random_number ( harvest = theta )
  theta = 2.0D+00 * r8_pi * theta

  vran(1) = cos ( theta )
  vran(2) = sin ( theta )

  return
end
subroutine direction_uniform_3d ( vran )

!*****************************************************************************80
!
!! DIRECTION_UNIFORM_3D picks a random direction vector in 3D.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) VRAN(3), the random direction vector,
!    with unit norm.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) phi
  real ( kind = rk ) r8_acos
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) theta
  real ( kind = rk ) vdot
  real ( kind = rk ) vran(dim_num)
!
!  Pick a uniformly random VDOT, which must be between -1 and 1.
!  This represents the dot product of the random vector with the Z unit vector.
!
!  Note: this works because the surface area of the sphere between
!  Z and Z + dZ is independent of Z.  So choosing Z uniformly chooses
!  a patch of area uniformly.
!
  call random_number ( harvest = vdot )
  vdot = 2.0D+00 * vdot - 1.0D+00

  phi = r8_acos ( vdot )
!
!  Pick a uniformly random rotation between 0 and 2 Pi around the
!  axis of the Z vector.
!
  call random_number ( harvest = theta )
  theta = 2.0D+00 * r8_pi * theta

  vran(1) = cos ( theta ) * sin ( phi )
  vran(2) = sin ( theta ) * sin ( phi )
  vran(3) = cos ( phi )

  return
end
subroutine direction_uniform_nd ( dim_num, w )

!*****************************************************************************80
!
!! DIRECTION_UNIFORM_ND generates a random direction vector in ND.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    13 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = rk ) W(DIM_NUM), a random direction vector,
!    with unit norm.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num

  real ( kind = rk ) norm
  real ( kind = rk ) w(dim_num)
!
!  Get N values from a standard normal distribution.
!
  call r8vec_normal_01 ( dim_num, w )
!
!  Compute the length of the vector.
!
  norm = sqrt ( sum ( w(1:dim_num)**2 ) )
!
!  Normalize the vector.
!
  w(1:dim_num) = w(1:dim_num) / norm

  return
end
subroutine disk_point_dist_3d ( pc, r, axis, p, dist )

!*****************************************************************************80
!
!! DISK_POINT_DIST_3D determines the distance from a disk to a point in 3D.
!
!  Discussion:
!
!    A disk in 3D satisfies the equations:
!
!      ( P(1) - PC(1) )^2 + ( P(2) - PC(2) )^2 + ( P(3) - PC(3) <= R^2
!
!    and
!
!      P(1) * AXIS(1) + P(2) * AXIS(2) + P(3) * AXIS(3) = 0
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) PC(3), the center of the disk.
!
!    Input, real ( kind = rk ) R, the radius of the disk.
!
!    Input, real ( kind = rk ) AXIS(3), the axis vector.
!
!    Input, real ( kind = rk ) P(3), the point to be checked.
!
!    Output, real ( kind = rk ) DIST, the distance of the point to the disk.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) axial_component
  real ( kind = rk ) axis(dim_num)
  real ( kind = rk ) axis_length
  real ( kind = rk ) dist
  real ( kind = rk ) r8vec_norm
  real ( kind = rk ) off_axis_component
  real ( kind = rk ) off_axis(dim_num)
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) r
!
!  Special case: the point is the center.
!
  if ( all ( p(1:dim_num) == pc(1:dim_num) ) ) then
    dist = 0.0D+00
    return
  end if

  axis_length = r8vec_norm ( dim_num, axis(1:dim_num) )

  if ( axis_length == 0.0D+00 ) then
    dist = -huge ( dist )
    return
  end if

  axial_component = dot_product ( p(1:dim_num) - pc(1:dim_num), &
    axis(1:dim_num) ) / axis_length
!
!  Special case: the point satisfies the disk equation exactly.
!
  if ( sum ( p(1:dim_num) - pc(1:dim_num) )**2 <= r * r .and. &
        axial_component == 0.0D+00 ) then
    dist = 0.0D+00
    return
  end if
!
!  Decompose P-PC into axis component and off-axis component.
!
  off_axis(1:dim_num) = p(1:dim_num) - pc(1:dim_num) &
    - axial_component * axis(1:dim_num) / axis_length

  off_axis_component = r8vec_norm ( dim_num, off_axis )
!
!  If the off-axis component has norm less than R, the nearest point is
!  the projection to the disk along the axial direction, and the distance 
!  is just the dot product of P-PC with unit AXIS.
!
  if ( off_axis_component <= r ) then
    dist = abs ( axial_component )
    return
  end if
!
!  Otherwise, the nearest point is along the perimeter of the disk.
!
  dist = sqrt ( axial_component**2 + ( off_axis_component - r )**2 )

  return
end
subroutine dms_to_radians ( degrees, minutes, seconds, radians )

!*****************************************************************************80
!
!! DMS_TO_RADIANS converts an angle from degrees/minutes/seconds to radians.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DEGREES, MINUTES, SECONDS, an angle in 
!    degrees, minutes, and seconds.
!
!    Output, real ( kind = rk ) RADIANS, the equivalent angle in radians.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) angle
  integer degrees
  integer minutes
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) radians
  integer seconds

  angle =   real ( degrees, kind = rk ) &
        + ( real ( minutes, kind = rk ) &
        + ( real ( seconds, kind = rk ) / 60.0D+00 ) ) / 60.0D+00

  radians = ( angle / 180.0D+00 ) * r8_pi

  return
end
subroutine dodec_shape_3d ( point_num, face_num, face_order_max, &
  point_coord, face_order, face_point )

!*****************************************************************************80
!
!! DODEC_SHAPE_3D describes a dodecahedron in 3D.
!
!  Discussion:
!
!    The vertices lie on the unit sphere.
!
!    The dual of a dodecahedron is an icosahedron.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer POINT_NUM, the number of points.
!
!    Input, integer FACE_NUM, the number of faces.
!
!    Input, integer FACE_ORDER_MAX, the maximum number of vertices
!    per face.
!
!    Output, real ( kind = rk ) POINT_COORD(3,POINT_NUM), the vertices.
!
!    Output, integer FACE_ORDER[FACE_NUM], the number of vertices
!    per face.
!
!    Output, integer FACE_POINT(FACE_ORDER_MAX,POINT_NUM); 
!    FACE_POINT(I,J) contains the index of the I-th point in the J-th face.  
!    The points are listed in the counter clockwise direction defined
!    by the outward normal at the face.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer face_num
  integer face_order_max
  integer point_num

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  integer face_order(face_num)
  integer face_point(face_order_max,face_num)
  real ( kind = rk ) phi
  real ( kind = rk ) point_coord(dim_num,point_num)
  real ( kind = rk ) z
!
!  Set point coordinates.
!
  phi = 0.5D+00 * ( sqrt ( 5.0D+00 ) + 1.0D+00 )

  a = 1.0D+00 / sqrt ( 3.0D+00 )
  b = phi / sqrt ( 3.0D+00 )
  c = ( phi - 1.0D+00 ) / sqrt ( 3.0D+00 )
  z = 0.0D+00

  point_coord(1:dim_num,1:point_num) = reshape ( (/ &
      a,  a,  a, &
      a,  a, -a, &
      a, -a,  a, &
      a, -a, -a, &
     -a,  a,  a, &
     -a,  a, -a, &
     -a, -a,  a, &
     -a, -a, -a, &
      c,  b,  z, &
     -c,  b,  z, &
      c, -b,  z, &
     -c, -b,  z, &
      b,  z,  c, &
      b,  z, -c, &
     -b,  z,  c, &
     -b,  z, -c, &
      z,  c,  b, &
      z, -c,  b, &
      z,  c, -b, &
      z, -c, -b /), (/ dim_num, point_num /) )
!
!  Set the face orders.
!
  face_order(1:face_num) = (/ &
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5 /)
!
!  Set the faces.
!
  face_point(1:face_order_max,1:face_num) = reshape ( (/ &
      2,  9,  1, 13, 14, &
      5, 10,  6, 16, 15, &
      3, 11,  4, 14, 13, &
      8, 12,  7, 15, 16, &
      3, 13,  1, 17, 18, &
      2, 14,  4, 20, 19, &
      5, 15,  7, 18, 17, &
      8, 16,  6, 19, 20, &
      5, 17,  1,  9, 10, &
      3, 18,  7, 12, 11, &
      2, 19,  6, 10,  9, &
      8, 20,  4, 11, 12 /), (/ face_order_max, face_num /) )

  return
end
subroutine dodec_size_3d ( point_num, edge_num, face_num, face_order_max )

!*****************************************************************************80
!
!! DODEC_SIZE_3D gives "sizes" for a dodecahedron in 3D.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer POINT_NUM, the number of points.
!
!    Output, integer EDGE_NUM, the number of edges.
!
!    Output, integer FACE_NUM, the number of faces.
!
!    Output, integer FACE_ORDER_MAX, the maximum order of any face.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer edge_num
  integer face_num
  integer face_order_max
  integer point_num

  point_num = 20
  edge_num = 30
  face_num = 12
  face_order_max = 5

  return
end
subroutine dual_shape_3d ( point_num, face_num, face_order_max, &
  point_coord, face_order, face_point, point_num2, face_num2, &
  face_order_max2, point_coord2, face_order2, face_point2 )

!*****************************************************************************80
!
!! DUAL_SHAPE_3D constructs the dual of a shape in 3D.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer POINT_NUM, the number of points.
!
!    Input, integer FACE_NUM, the number of faces.
!
!    Input, integer FACE_ORDER_MAX, the maximum number of vertices
!    per face.
!
!    Input, real ( kind = rk ) POINT_COORD(3,POINT_NUM), the points.
!
!    Input, integer FACE_ORDER(FACE_NUM), the number of vertices
!    per face.
!
!    Input, integer FACE_POINT(FACE_ORDER_MAX,FACE_NUM); 
!    FACE_POINT(I,J) is the index of the I-th point in the J-th face.  The
!    points are listed in the counter clockwise direction defined
!    by the outward normal at the face.
!
!    Input, integer POINT_NUM2, the number of points in the dual.
!
!    Input, integer FACE_NUM2, the number of faces in the dual.
!
!    Input, integer FACE_ORDER_MAX2, the maximum number of 
!    vertices per face in the dual.
!
!    Output, real ( kind = rk ) POINT_COORD2(3,POINT_NUM2), the point 
!    coordinates of the dual.
!
!    Output, integer FACE_ORDER2(FACE_NUM2), the number of 
!    vertices per face.
!
!    Output, integer FACE_POINT2(FACE_ORDER_MAX2,FACE_NUM2), 
!    the vertices of each face in the dual.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer face_num
  integer face_num2
  integer face_order_max
  integer face_order_max2
  integer, parameter :: dim_num = 3
  integer point_num
  integer point_num2

  integer col
  integer face
  integer face_order(face_num)
  integer face_order2(face_num2)
  integer face_point(face_order_max,face_num)
  integer face_point2(face_order_max2,face_num2)
  integer i
  integer inext
  integer iprev
  integer istop
  integer j
  integer k
  real ( kind = rk ) norm
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) point_coord(dim_num,point_num)
  real ( kind = rk ) point_coord2(dim_num,point_num2)
  integer row
!
!  This computation should really compute the center of gravity
!  of the face, in the general case.
!
!  We'll also assume the vertices of the original and the dual
!  are to lie on the unit sphere, so we can normalize the
!  position vector of the vertex.
!
  do face = 1, face_num

    p(1:dim_num) = 0.0D+00

    do j = 1, face_order(face)
      k = face_point(j,face)
      p(1:dim_num) = p(1:dim_num) + point_coord(1:dim_num,k)
    end do

    norm = sqrt ( sum ( p(1:dim_num)**2 ) )

    point_coord2(1:dim_num,face) = p(1:dim_num) / norm

  end do
!
!  Now build the face in the dual associated with each node FACE.
!
  do face = 1, face_num2
!
!  Initialize the order.
!
    face_order2(face) = 0
!
!  Find the first occurrence of FACE in an edge of polyhedron.
!
    call i4col_find_item ( face_order_max, face_num, face_point, &
      face, row, col )

    if ( row <= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DUAL_SHAPE_3D - Fatal error!'
      write ( *, '(a,i8)' ) '  Could not find an edge using node ', face
      stop 1
    end if
!
!  Save the following node as ISTOP.
!  When we encounter ISTOP again, this will mark the end of our search.
!
    i = row + 1
    if ( face_order(col) < i ) then
      i = 1
    end if

    istop = face_point(i,col)
!
!  Save the previous node as INEXT.
!
    do

      i = row - 1
      if ( i < 1 ) then
        i = i + face_order(col)
      end if

      inext = face_point(i,col)

      face_order2(face) = face_order2(face) + 1

      face_point2(face_order2(face),face) = col
!
!  If INEXT =/= ISTOP, continue.
!
      if ( inext == istop ) then
        exit
      end if
!
!  Set IPREV:= INEXT.
!
      iprev = inext
!
!  Search for the occurrence of the edge FACE-IPREV.
!
      call i4col_find_pair_wrap ( face_order_max, face_num, face_point, &
        face, iprev, row, col )

      if ( row <= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DUAL_SHAPE_3D - Fatal error!'
        write ( *, '(a,i8)' ) '  No edge from node ', iprev
        write ( *, '(a,i8)' ) '  to node ', face
        stop 1
      end if

    end do

  end do

  return
end
subroutine dual_size_3d ( point_num, edge_num, face_num, face_order_max, &
  point_coord, face_order, face_point, point_num2, edge_num2, face_num2, &
  face_order_max2 )

!*****************************************************************************80
!
!! DUAL_SIZE_3D determines sizes for a dual of a shape in 3D.
!
!  Discussion:
!
!    We don't actually need FACE_POINT as input here.  But since the
!    three arrays occur together everywhere else, it seems unnecessarily
!    user-confusing to vary the usage here!
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer POINT_NUM, the number of points.
!
!    Input, integer EDGE_NUM, the number of edges.
!
!    Input, integer FACE_NUM, the number of faces.
!
!    Input, integer FACE_ORDER_MAX, the maximum number of vertices
!    per face.
!
!    Input, real ( kind = rk ) POINT_COORD(3,POINT_NUM), the points.
!
!    Input, integer FACE_ORDER(FACE_NUM), the number of vertices
!    per face.
!
!    Input, integer FACE_POINT(FACE_ORDER_MAX,FACE_NUM); 
!    FACE_POINT(I,J) is the index of the I-th point in the J-th face.  The
!    points are listed in the counter clockwise direction defined
!    by the outward normal at the face.
!
!    Output, integer POINT_NUM2, the number of points in the dual.
!
!    Output, integer EDGE_NUM2, the number of edges in the dual.
!
!    Output, integer FACE_NUM2, the number of faces in the dual.
!
!    Output, integer FACE_ORDER_MAX2, the maximum number of 
!    vertices per face in the dual.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer face_num
  integer face_order_max
  integer point_num

  integer edge_num
  integer edge_num2
  integer face
  integer face_num2
  integer face_order(face_num)
  integer face_order2(point_num)
  integer face_order_max2
  integer face_point(face_order_max,face_num)
  integer face2
  integer i
  integer point_num2
  real ( kind = rk ) point_coord(dim_num,point_num)

  call r8_fake_use ( point_coord(1,1) )
!
!  These values are easy to compute:
!
  point_num2 = face_num
  edge_num2 = edge_num
  face_num2 = point_num
!
!  To determine FACE_ORDER_MAX2 is not so easy.
!  You have to construct the FACE_ORDER array for the dual shape.
!  The order of a dual face is the number of edges that the vertex occurs in.
!  But then all we have to do is count how many times each item shows up
!  in the FACE_POINT array.
!
  face_order_max2 = 0
  face_order2(1:face_num2) = 0

  do face = 1, face_num
    do i = 1, face_order(face)
      face2 = face_point(i,face)
      face_order2(face2) = face_order2(face2) + 1
    end do
  end do

  face_order_max2 = maxval ( face_order2(1:face_num2) )

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! get_unit() returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
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
!  Parameters:
!
!    Output, integer IUNIT.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5 and 6).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer ios
  integer iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
subroutine glob2loc_3d ( cospitch, cosroll, cosyaw, sinpitch, sinroll, sinyaw, &
  globas, glopts, locpts )

!*****************************************************************************80
!
!! GLOB2LOC_3D converts from a global to a local coordinate system in 3D.
!
!  Discussion:
!
!    A global coordinate system is given.
!
!    A local coordinate system has been translated to the point with
!    global coordinates GLOBAS, and rotated through a yaw, a pitch, and
!    a roll.
!
!    A point has global coordinates GLOPTS, and it is desired to know
!    the point's local coordinates LOCPTS.
!
!    The transformation may be written as
!
!      LOC = M_ROLL * M_PITCH * M_YAW * ( GLOB - GLOBAS )
!
!    where
!
!               (       1            0            0      )
!    M_ROLL =   (       0        cos(Roll)    sin(Roll)  )
!               (       0      - sin(Roll)    cos(Roll)  )
!
!               (   cos(Pitch)       0      - sin(Pitch) )
!    M_PITCH =  (       0            1            0      )
!               (   sin(Pitch)       0        cos(Pitch) )
!
!               (   cos(Yaw)     sin(Yaw)         0      )
!    M_YAW    = ( - sin(Yaw)     cos(Yaw)         0      )
!               (       0            0            1      )
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) COSPITCH, COSROLL, COSYAW, the cosines of
!    the pitch, roll and yaw angles.
!
!    Input, real ( kind = rk ) SINPITCH, SINROLL, SINYAW, the sines of the pitch,
!    roll and yaw angles.
!
!    Input, real ( kind = rk ) GLOBAS(3), the global base vector.
!
!    Input, real ( kind = rk ) GLOPTS(3), the global coordinates
!    of the point whose coordinates are to be transformed.
!
!    Output, real ( kind = rk ) LOCPTS(3), the local coordinates of the point
!    whose global coordinates were given in GLOPTS.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) cospitch
  real ( kind = rk ) cosroll
  real ( kind = rk ) cosyaw
  real ( kind = rk ) globas(dim_num)
  real ( kind = rk ) glopts(dim_num)
  real ( kind = rk ) locpts(dim_num)
  real ( kind = rk ) sinpitch
  real ( kind = rk ) sinroll
  real ( kind = rk ) sinyaw

  locpts(1) = ( cosyaw * cospitch ) * ( glopts(1) - globas(1) ) &
            + ( sinyaw * cospitch ) * ( glopts(2) - globas(2) ) &
            -   sinpitch * ( glopts(3) - globas(3) )

  locpts(2) = ( cosyaw * sinpitch * sinroll - sinyaw * cosroll ) &
    * ( glopts(1) - globas(1) ) &
    + ( sinyaw * sinpitch * sinroll + cosyaw * cosroll ) &
    * ( glopts(2) - globas(2) ) &
    +   cospitch * sinroll * ( glopts(3) - globas(3) )

  locpts(3) = ( cosyaw * sinpitch * cosroll + sinyaw * sinroll ) &
    * ( glopts(1) - globas(1) ) &
    + ( sinyaw * sinpitch * cosroll - cosyaw * sinroll  ) &
    * ( glopts(2) - globas(2) ) &
    + ( cospitch * cosroll ) * ( glopts(3) - globas(3) )

  return
end
function halfplane_contains_point_2d ( p1, p2, p )

!*****************************************************************************80
!
!! HALFPLANE_CONTAINS_POINT_2D reports if a half-plane contains a point in 2d.
!
!  Discussion:
!
!    The halfplane is assumed to be all the points "to the left" of the
!    line that passes from P1 through P2.  Thus, one way to
!    understand where the point P is, is to compute the signed
!    area of the triangle ( P1, P2, P ).
!
!    If this area is
!      positive, the point is strictly inside the halfplane,
!      zero, the point is on the boundary of the halfplane,
!      negative, the point is strictly outside the halfplane.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(2), P2(2), two distinct points
!    on the line defining the half plane.
!
!    Input, real ( kind = rk ) P(2), the point to be checked.
!
!    Output, logical HALFPLANE_CONTAINS_POINT_2D, is TRUE if 
!    the halfplane contains the point.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) area_signed
  logical halfplane_contains_point_2d
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)

  area_signed = 0.5D+00 *       &
    ( p1(1) * ( p2(2) - p(2)  ) &
    + p2(1) * ( p(2)  - p1(2) ) &
    + p(1)  * ( p1(2) - p2(2) ) )

  halfplane_contains_point_2d = ( 0.0D+00 <= area_signed )

  return
end
subroutine halfspace_imp_triangle_int_3d ( a, b, c, d, t, int_num, pint )

!*****************************************************************************80
!
!! HALFSPACE_IMP_TRIANGLE_INT_3D: intersection ( imp halfspace, triangle ).
!
!  Discussion:
!
!    The implicit form of a half-space in 3D may be described as the set
!    of points P on or "above" an implicit plane:
!
!      0 <= A * P(1) + B * P(2) + C * P(3) + D
!
!    The triangle is specified by listing its three vertices.
!
!    The intersection may be described by the number of vertices of the
!    triangle that are included in the halfspace, and by the location of
!    points between vertices that separate a side of the triangle into
!    an included part and an unincluded part.
!
!    0 vertices, 0 separators    (no intersection)
!    1 vertex, 0 separators      (point intersection)
!    2 vertices, 0 separators    (line intersection)
!    3 vertices, 0 separators    (triangle intersection)
!
!    1 vertex, 2 separators,     (intersection is a triangle)
!    2 vertices, 2 separators,   (intersection is a quadrilateral).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, C, D, the parameters that define the
!    implicit plane, which in turn define the implicit halfspace.
!
!    Input, real ( kind = rk ) T(3,3), the vertices of the triangle.
!
!    Output, integer INT_NUM, the number of intersection points 
!    returned, which will always be between 0 and 4.
!
!    Output, real ( kind = rk ) PINT(3,4), the coordinates of the INT_NUM
!    intersection points.  The points will lie in sequence on the triangle.
!    Some points will be vertices, and some may be separators.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) d
  real ( kind = rk ) dist1
  real ( kind = rk ) dist2
  real ( kind = rk ) dist3
  integer int_num
  real ( kind = rk ) pint(dim_num,4)
  real ( kind = rk ) t(dim_num,3)
!
!  Compute the signed distances between the vertices and the plane.
!
  dist1 = a * t(1,1) + b * t(2,1) + c * t(3,1) + d
  dist2 = a * t(1,2) + b * t(2,2) + c * t(3,2) + d
  dist3 = a * t(1,3) + b * t(2,2) + c * t(3,3) + d
!
!  Now we can find the intersections.
!
  call halfspace_triangle_int_3d ( dist1, dist2, dist3, t, int_num, pint )

  return
end
subroutine halfspace_normal_triangle_int_3d ( pp, normal, t, int_num, pint )

!*****************************************************************************80
!
!! HALFSPACE_NORMAL_TRIANGLE_INT_3D: intersection ( norm halfspace, triangle ).
!
!  Discussion:
!
!    The normal form of a halfspace in 3D may be described as the set
!    of points P on or "above" a plane described in normal form:
!
!      PP is a point on the plane,
!      NORMAL is the unit normal vector, pointing "out" of the
!      halfspace.
!
!    The triangle is specified by listing its three vertices.
!
!    The intersection may be described by the number of vertices of the
!    triangle that are included in the halfspace, and by the location of
!    points between vertices that separate a side of the triangle into
!    an included part and an unincluded part.
!
!    0 vertices, 0 separators    (no intersection)
!    1 vertex, 0 separators      (point intersection)
!    2 vertices, 0 separators    (line intersection)
!    3 vertices, 0 separators    (triangle intersection)
!
!    1 vertex, 2 separators,     (intersection is a triangle)
!    2 vertices, 2 separators,   (intersection is a quadrilateral).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) PP(3), a point on the bounding plane
!    that defines the halfspace.
!
!    Input, real ( kind = rk ) NORMAL(3), the components of the normal vector
!    to the bounding plane that defines the halfspace.  By convention, the
!    normal vector points "outwards" from the halfspace.
!
!    Input, real ( kind = rk ) T(3,3), the vertices of the triangle.
!
!    Output, integer INT_NUM, the number of intersection points 
!    returned, which will always be between 0 and 4.
!
!    Output, real ( kind = rk ) PINT(3,4), the coordinates of the INT_NUM
!    intersection points.  The points will lie in sequence on the triangle.
!    Some points will be vertices, and some may be separators.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) d
  real ( kind = rk ) dist1
  real ( kind = rk ) dist2
  real ( kind = rk ) dist3
  real ( kind = rk ) normal(dim_num)
  integer int_num
  real ( kind = rk ) pp(dim_num)
  real ( kind = rk ) pint(dim_num,4)
  real ( kind = rk ) t(dim_num,3)
!
!  Compute the signed distances between the vertices and the plane.
!
  d = - dot_product ( normal(1:dim_num), pp(1:dim_num) )
!
!  Compute the signed distances between the vertices and the plane.
!
  dist1 = d + dot_product ( normal(1:dim_num), t(1:dim_num,1) )
  dist2 = d + dot_product ( normal(1:dim_num), t(1:dim_num,2) )
  dist3 = d + dot_product ( normal(1:dim_num), t(1:dim_num,3) )
!
!  Now we can find the intersections.
!
  call halfspace_triangle_int_3d ( dist1, dist2, dist3, t, int_num, pint )

  return
end
subroutine halfspace_triangle_int_3d ( dist1, dist2, dist3, t, int_num, pint )

!*****************************************************************************80
!
!! HALFSPACE_TRIANGLE_INT_3D: intersection ( halfspace, triangle ) in 3D.
!
!  Discussion:
!
!    The triangle is specified by listing its three vertices.
!
!    The halfspace is not described in the input data.  Rather, the
!    distances from the triangle vertices to the halfspace are given.
!
!    The intersection may be described by the number of vertices of the
!    triangle that are included in the halfspace, and by the location of
!    points between vertices that separate a side of the triangle into
!    an included part and an unincluded part.
!
!    0 vertices, 0 separators    (no intersection)
!    1 vertex, 0 separators      (point intersection)
!    2 vertices, 0 separators    (line intersection)
!    3 vertices, 0 separators    (triangle intersection)
!
!    1 vertex, 2 separators,     (intersection is a triangle)
!    2 vertices, 2 separators,   (intersection is a quadrilateral).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) DIST1, DIST2, DIST3, the distances from each of
!    the three vertices of the triangle to the halfspace.  The distance is
!    zero if a vertex lies within the halfspace, or on the plane that
!    defines the boundary of the halfspace.  Otherwise, it is the
!    distance from that vertex to the bounding plane.
!
!    Input, real ( kind = rk ) T(3,3), the vertices of the triangle.
!
!    Output, integer INT_NUM, the number of intersection points
!    returned, which will always be between 0 and 4.
!
!    Output, real ( kind = rk ) PINT(3,4), the coordinates of the INT_NUM
!    intersection points.  The points will lie in sequence on the triangle.
!    Some points will be vertices, and some may be separators.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) dist1
  real ( kind = rk ) dist2
  real ( kind = rk ) dist3
  integer int_num
  real ( kind = rk ) pint(dim_num,4)
  real ( kind = rk ) t(dim_num,3)
!
!  Walk around the triangle, looking for vertices that are included,
!  and points of separation.
!
  int_num = 0

  if ( dist1 <= 0.0D+00 ) then

    int_num = int_num + 1
    pint(1:dim_num,int_num) = t(1:dim_num,1)

  end if

  if ( dist1 * dist2 < 0.0D+00 ) then

    int_num = int_num + 1
    pint(1:dim_num,int_num) = &
      ( dist1 * t(1:dim_num,2) - dist2 * t(1:dim_num,1) ) &
      / ( dist1 - dist2 )

  end if

  if ( dist2 <= 0.0D+00 ) then

    int_num = int_num + 1
    pint(1:dim_num,int_num) = t(1:dim_num,2)

  end if

  if ( dist2 * dist3 < 0.0D+00 ) then

    int_num = int_num + 1

    pint(1:dim_num,int_num) = &
      ( dist2 * t(1:dim_num,3) - dist3 * t(1:dim_num,2) ) &
      / ( dist2 - dist3 )

  end if

  if ( dist3 <= 0.0D+00 ) then

    int_num = int_num + 1
    pint(1:dim_num,int_num) = t(1:dim_num,3)

  end if

  if ( dist3 * dist1 < 0.0D+00 ) then

    int_num = int_num + 1
    pint(1:dim_num,int_num) = &
      ( dist3 * t(1:dim_num,1) - dist1 * t(1:dim_num,3) ) &
      / ( dist3 - dist1 )

  end if

  return
end
function haversine ( a )

!*****************************************************************************80
!
!! HAVERSINE computes the haversine of an angle.
!
!  Discussion:
!
!    haversine(A) = ( 1 - cos ( A ) ) / 2
!
!    The haversine is useful in spherical trigonometry.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, the angle.
!
!    Output, real ( kind = rk ) HAVERSINE, the haversine of the angle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) haversine

  haversine = ( 1.0D+00 - cos ( a ) ) / 2.0D+00

  return
end
subroutine helix_shape_3d ( a, n, r, theta1, theta2, p )

!*****************************************************************************80
!
!! HELIX_SHAPE_3D computes points on a helix in 3D.
!
!  Discussion:
!
!    The user specifies the parameters A and R, the first and last
!    THETA values, and the number of equally spaced THETA values
!    at which point values are to be computed.
!
!    X = R * COS ( THETA )
!    Y = R * SIN ( THETA )
!    Z = A * THETA
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, the rate at which Z advances with THETA.
!
!    Input, integer N, the number of points to compute on 
!    the helix.
!
!    Input, real ( kind = rk ) R, the radius of the helix.
!
!    Input, real ( kind = rk ) THETA1, THETA2, the first and last THETA values at
!    which to compute points on the helix.  THETA is measured in
!    radians.
!
!    Output, real ( kind = rk ) P(3,N), the coordinates of points on the helix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: dim_num = 3

  real ( kind = rk ) a
  integer i
  real ( kind = rk ) p(dim_num,n)
  real ( kind = rk ) r
  real ( kind = rk ) theta
  real ( kind = rk ) theta1
  real ( kind = rk ) theta2

  do i = 1, n

    if ( n == 1 ) then
      theta = 0.5D+00 * ( theta1 + theta2 )
    else
      theta = ( real ( n - i,     kind = rk ) * theta1 &
              + real (     i - 1, kind = rk ) * theta2 ) &
              / real ( n     - 1, kind = rk )
    end if

    p(1,i) = r * cos ( theta )
    p(2,i) = r * sin ( theta )
    p(3,i) = a * theta

  end do

  return
end
function hexagon_area_2d ( r )

!*****************************************************************************80
!
!! HEXAGON_AREA_2D returns the area of a regular hexagon in 2D.
!
!  Discussion:
!
!    The radius of a regular hexagon is the distance from the center
!    of the hexagon to any vertex.  This happens also to equal the 
!    length of any side.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) R, the radius of the hexagon.
!
!    Output, real ( kind = rk ) HEXAGON_AREA_2D, the area of the hexagon.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) hexagon_area_2d
  real ( kind = rk ) hexagon01_area_2d
  real ( kind = rk ) r

  hexagon_area_2d = r * r * hexagon01_area_2d ( )

  return
end
function hexagon_contains_point_2d ( v, p )

!*****************************************************************************80
!
!! HEXAGON_CONTAINS_POINT_2D finds if a point is inside a hexagon in 2D.
!
!  Discussion:
!
!    This test is only valid if the hexagon is convex.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) V(2,6), the vertices, in counter clockwise order.
!
!    real ( kind = rk ) P(2), the point to be tested.
!
!  Output:
!
!    logical HEXAGON_CONTAINS_POINT_2D, is TRUE 
!    if X is in the hexagon.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 6
  integer, parameter :: dim_num = 2

  logical hexagon_contains_point_2d
  integer i
  integer j
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) v(dim_num,n)
!
!  A point is inside a convex hexagon if and only if it is "inside"
!  each of the 6 halfplanes defined by lines through consecutive
!  vertices.
!
  do i = 1, n

    j = mod ( i, n ) + 1

    if (  v(1,i) * ( v(2,j) - p(2  ) ) &
        + v(1,j) * ( p(2  ) - v(2,i) ) &
        + p(1  ) * ( v(2,i) - v(2,j) ) < 0.0D+00 ) then

      hexagon_contains_point_2d = .false.
      return

    end if

  end do

  hexagon_contains_point_2d = .true.

  return
end
subroutine hexagon_shape_2d ( angle_degrees, p )

!*****************************************************************************80
!
!! HEXAGON_SHAPE_2D returns points on the unit regular hexagon in 2D.
!
!  Diagram:
!
!      120_____60
!        /     \
!    180/       \0
!       \       /
!        \_____/
!      240     300
!
!  Discussion:
!
!    The unit regular hexagon has radius 1.  The radius is the distance from
!    the center to any vertex, and it is also the length of any side.  
!    An example of a unit hexagon is the convex hull of the points:
!
!      (   1,              0 ),
!      (   0.5,   sqrt (3)/2 ),
!      ( - 0.5,   sqrt (3)/2 ),
!      ( - 1,              0 ),
!      ( - 0.5, - sqrt (3)/2 ),
!      (   0.5, - sqrt (3)/2 ).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) angle_degrees, the angle, in degrees, of the point.
!
!  Output:
!
!    real ( kind = rk ) P(2), the coordinates of the point.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) angle_degrees
  real ( kind = rk ) angle2
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) r8_cotd
  real ( kind = rk ) r8_modp
  real ( kind = rk ) r8_tand
!
!  Ensure that 0 <= ANGLE < 360.
!
  angle2 = r8_modp ( angle_degrees, 360.0D+00 )
!
!  y = - sqrt(3) * x + sqrt(3)
!
  if ( 0.0D+00 <= angle2 .and. angle2 <= 60.0D+00 ) then

    p(1) = sqrt ( 3.0D+00 ) / ( r8_tand ( angle2 ) + sqrt ( 3.0D+00 ) )
    p(2) = r8_tand ( angle2 ) * p(1)
!
!  y = sqrt(3) / 2
!
  else if ( angle2 <= 120.0D+00 ) then

    p(2) = sqrt ( 3.0D+00 ) / 2.0D+00
    p(1) = r8_cotd ( angle2 ) * p(2)
!
!  y = sqrt(3) * x + sqrt(3)
!
  else if ( angle2 <= 180.0D+00 ) then

    p(1) = sqrt ( 3.0D+00 ) / ( r8_tand ( angle2 ) - sqrt ( 3.0D+00 ) )
    p(2) = r8_tand ( angle2 ) * p(1)
!
!  y = - sqrt(3) * x - sqrt(3)
!
  else if ( angle2 <= 240.0D+00 ) then

    p(1) = - sqrt ( 3.0D+00 ) / ( r8_tand ( angle2 ) + sqrt ( 3.0D+00 ) )
    p(2) = r8_tand ( angle2 ) * p(1)
!
!  y = - sqrt(3) / 2
!
  else if ( angle2 <= 300.0D+00 ) then

    p(2) = - sqrt ( 3.0D+00 ) / 2.0D+00
    p(1) = r8_cotd ( angle2 ) * p(2)
!
!  y = sqrt(3) * x - sqrt(3)
!
  else if ( angle2 <= 360.0D+00 ) then

    p(1) = - sqrt ( 3.0D+00 ) / ( r8_tand ( angle2 ) - sqrt ( 3.0D+00 ) )
    p(2) = r8_tand ( angle2 ) * p(1)

  end if

  return
end
subroutine hexagon_vertices_2d ( p )

!*****************************************************************************80
!
!! HEXAGON_VERTICES_2D returns the vertices of the unit hexagon in 2D.
!
!  Discussion:
!
!    The unit hexagon has maximum radius 1, and is the hull of the points
!
!      (   1,              0 ),
!      (   0.5,   sqrt (3)/2 ),
!      ( - 0.5,   sqrt (3)/2 ),
!      ( - 1,              0 ),
!      ( - 0.5, - sqrt (3)/2 ),
!      (   0.5, - sqrt (3)/2 ).
!
!  Diagram:
!
!      120_____60
!        /     \
!    180/       \0
!       \       /
!        \_____/
!      240     300
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) P(2,6), the coordinates of the vertices.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ), parameter :: a = 0.8660254037844386D+00
  real ( kind = rk ) p(dim_num,6)

  p(1:2,1:6) = reshape ( (/ &
     1.0D+00,  0.0D+00, &
     0.5D+00,  a, &
    -0.5D+00,  a, &
    -1.0D+00,  0.0D+00, &
    -0.5D+00, -a, &
     0.5D+00, -a /), (/ dim_num, 6 /) )

  return
end
function hexagon01_area_2d ( )

!*****************************************************************************80
!
!! HEXAGON01_AREA_2D returns the area of a unit regular hexagon in 2D.
!
!  Discussion:
!
!    A "unit" regular hexagon has both a "radius" of 1 (distance
!    from the center to any vertex), and a side length of 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) HEXAGON01_AREA_2D, the area of the hexagon.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) hexagon01_area_2d

  hexagon01_area_2d = 3.0D+00 * sqrt ( 3.0D+00 ) / 2.0D+00

  return
end
function hyperball01_volume ( m )

!*****************************************************************************80
!
!! HYPERBALL01_VOLUME returns the volume of the unit hyperball in M dimensions.
!
!  Discussion:
!
!     M  Volume
!
!     1    2
!     2    1        * PI
!     3  ( 4 /   3) * PI
!     4  ( 1 /   2) * PI^2
!     5  ( 8 /  15) * PI^2
!     6  ( 1 /   6) * PI^3
!     7  (16 / 105) * PI^3
!     8  ( 1 /  24) * PI^4
!     9  (32 / 945) * PI^4
!    10  ( 1 / 120) * PI^5
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 January 2018
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the spatial dimension.
!
!    Output, real ( kind = rk ) HYPERBALL01_VOLUME, the volume of the unit ball.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) hyperball01_volume
  integer i
  integer m
  integer m_half
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) volume

  if ( mod ( m, 2 ) == 0 ) then
    m_half = ( m / 2 )
    volume = r8_pi ** m_half
    do i = 1, m_half
      volume = volume / real ( i, kind = rk )
    end do
  else
    m_half = ( ( m - 1 ) / 2 )
    volume = r8_pi ** m_half * 2.0D+00 ** m
    do i = m_half + 1, 2 * m_half + 1
      volume = volume / real ( i, kind = rk )
    end do
  end if

  hyperball01_volume = volume

  return
end
function i4_dedekind_factor ( p, q )

!*****************************************************************************80
!
!! I4_DEDEKIND_FACTOR computes a function needed for a Dedekind sum. 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 July 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Hans Rademacher, Emil Grosswald,
!    Dedekind Sums,
!    Mathematics Association of America, 1972,
!    LC: QA241.R2.
!
!  Parameters:
!
!    Input, integer P, Q, two positive integers.
!
!    Input, real ( kind = rk ) I4_DEDEKIND_FACTOR, the Dedekind factor of P / Q.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) i4_dedekind_factor
  integer p
  integer q

  if ( mod ( p, q ) == 0 ) then
    i4_dedekind_factor = 0.0D+00
  else
    i4_dedekind_factor = real ( p, kind = rk ) / real ( q, kind = rk ) &
      - real ( ( p / q ), kind = rk ) - 0.5D+00
  end if

  return
end
subroutine i4_dedekind_sum ( p, q, s )

!*****************************************************************************80
!
!! I4_DEDEKIND_SUM computes the Dedekind sum of two I4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 July 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Hans Rademacher, Emil Grosswald,
!    Dedekind Sums,
!    Mathematics Association of America, 1972,
!    LC: QA241.R2.
!
!  Parameters:
!
!    Input, integer P, Q, two positive integers.
!
!    Output, real ( kind = rk ) S, the Dedekind sum of P and Q.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  real ( kind = rk ) i4_dedekind_factor
  integer p
  integer q
  real ( kind = rk ) s

  s = 0.0D+00

  do i = 1, q
    s = s + i4_dedekind_factor ( i, q ) * i4_dedekind_factor ( p * i, q )
  end do

  return
end
function i4_factorial2 ( n )

!*****************************************************************************80
!
!! I4_FACTORIAL2 computes the double factorial function.
!
!  Discussion:
!
!    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
!                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the argument of the double factorial 
!    function.  If N is less than 1, I4_FACTORIAL2 is returned as 1.
!
!    Output, integer I4_FACTORIAL2, the value of N!!.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i4_factorial2
  integer n
  integer n_copy

  if ( n < 1 ) then
    i4_factorial2 = 1
    return
  end if

  n_copy = n
  i4_factorial2 = 1

  do while ( 1 < n_copy )
    i4_factorial2 = i4_factorial2 * n_copy
    n_copy = n_copy - 2
  end do

  return
end
subroutine i4_fake_use ( n )

!*****************************************************************************80
!
!! i4_fake_use pretends to use a variable.
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
!    integer N, the variable to be "used".
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  if ( n /= n ) then
    write ( *, '(a)' ) '  i4_fake_use: variable is NAN.'
  end if

  return
end
function i4_gcd ( i, j )

!*****************************************************************************80
!
!! I4_GCD finds the greatest common divisor of two I4's.
!
!  Discussion:
!
!    Note that only the absolute values of I and J are
!    considered, so that the result is always nonnegative.
!
!    If I or J is 0, I4_GCD is returned as max ( 1, abs ( I ), abs ( J ) ).
!
!    If I and J have no common factor, I4_GCD is returned as 1.
!
!    Otherwise, using the Euclidean algorithm, I4_GCD is the
!    greatest common divisor of I and J.
!
!    An I4 is an integer value.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    25 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, J, two numbers whose GCD is desired.
!
!    Output, integer I4_GCD, the greatest common divisor
!    of I and J.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer i4_gcd
  integer j
  integer p
  integer q
  integer r

  i4_gcd = 1
!
!  Return immediately if either I or J is zero.
!
  if ( i == 0 ) then
    i4_gcd = max ( 1, abs ( j ) )
    return
  else if ( j == 0 ) then
    i4_gcd = max ( 1, abs ( i ) )
    return
  end if
!
!  Set P to the larger of I and J, Q to the smaller.
!  This way, we can alter P and Q as we go.
!
  p = max ( abs ( i ), abs ( j ) )
  q = min ( abs ( i ), abs ( j ) )
!
!  Carry out the Euclidean algorithm.
!
  do

    r = mod ( p, q )

    if ( r == 0 ) then
      exit
    end if

    p = q
    q = r

  end do

  i4_gcd = q

  return
end
function i4_huge ( )

!*****************************************************************************80
!
!! I4_HUGE returns a "huge" I4.
!
!  Discussion:
!
!    On an IEEE 32 bit machine, I4_HUGE should be 2**31 - 1, and its
!    bit pattern should be
!
!     01111111111111111111111111111111
!
!    In this case, its numerical value is 2147483647.
!
!    Using the Dec/Compaq/HP Alpha FORTRAN compiler FORT, I could
!    use I4_HUGE() and HUGE interchangeably.
!
!    An I4 is an integer value.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer I4_HUGE, a "huge" I4.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i4_huge

  i4_huge = 2147483647

  return
end
function i4_lcm ( i, j )

!*****************************************************************************80
!
!! I4_LCM computes the least common multiple of two I4's.
!
!  Discussion:
!
!    The least common multiple may be defined as
!
!      LCM(I,J) = ABS( I * J ) / GCD(I,J)
!
!    where GCD(I,J) is the greatest common divisor of I and J.
!
!    An I4 is an integer value.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, J, the integers whose I4_LCM is desired.
!
!    Output, integer I4_LCM, the least common multiple of I and J.
!    I4_LCM is never negative.  I4_LCM is 0 if either I or J is zero.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer i4_gcd
  integer j
  integer i4_lcm

  i4_lcm = abs ( i * ( j / i4_gcd ( i, j ) ) )

  return
end
function i4_modp ( i, j )

!*****************************************************************************80
!
!! I4_MODP returns the nonnegative remainder of integer division.
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
!  Example:
!
!        I     J     MOD  I4_MODP    Factorization
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
!  Parameters:
!
!    Input, integer I, the number to be divided.
!
!    Input, integer J, the number that divides I.
!
!    Output, integer I4_MODP, the nonnegative remainder when I is
!    divided by J.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer i4_modp
  integer j

  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MODP - Fatal error!'
    write ( *, '(a,i8)' ) '  I4_MODP ( I, J ) called with J = ', j
    stop 1
  end if

  i4_modp = mod ( i, j )

  if ( i4_modp < 0 ) then
    i4_modp = i4_modp + abs ( j )
  end if

  return
end
subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! I4_SWAP switches two I4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer I, J.  On output, the values of I and
!    J have been interchanged.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer j
  integer k

  k = i
  i = j
  j = k

  return
end
function i4_uniform_ab ( a, b )

!*****************************************************************************80
!
!! I4_UNIFORM_AB returns a scaled pseudorandom I4.
!
!  Discussion:
!
!    An I4 is an integer value.
!
!    The pseudorandom number will be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer A, B, the limits of the interval.
!
!    Output, integer I4_UNIFORM_AB, a number between A and B.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer a
  integer b
  integer i4_uniform_ab
  real r
  integer value

  call random_number ( harvest = r )
!
!  Scale R to lie between A-0.5 and B+0.5.
!
  r = ( 1.0E+00 - r ) * ( real ( min ( a, b ) ) - 0.5E+00 ) & 
    +             r   * ( real ( max ( a, b ) ) + 0.5E+00 )
!
!  Use rounding to convert R to an integer between A and B.
!
  value = nint ( r )

  value = max ( value, min ( a, b ) )
  value = min ( value, max ( a, b ) )

  i4_uniform_ab = value

  return
end
function i4_wrap ( ival, ilo, ihi )

!*****************************************************************************80
!
!! I4_WRAP forces an I4 to lie between given limits by wrapping.
!
!  Example:
!
!    ILO = 4, IHI = 8
!
!    I  I4_WRAP
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
!    19 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer IVAL, an integer value.
!
!    Input, integer ILO, IHI, the desired bounds for the integer
!    value.
!
!    Output, integer I4_WRAP, a "wrapped" version of IVAL.
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
  integer wide

  jlo = min ( ilo, ihi )
  jhi = max ( ilo, ihi )

  wide = jhi - jlo + 1

  if ( wide == 1 ) then
    i4_wrap = jlo
  else
    i4_wrap = jlo + i4_modp ( ival - jlo, wide )
  end if

  return
end
subroutine i4col_compare ( m, n, a, i, j, isgn )

!*****************************************************************************80
!
!! I4COL_COMPARE compares columns I and J of an I4COL.
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, I = 2, J = 4
!
!      A = (
!        1  2  3  4
!        5  6  7  8
!        9 10 11 12 )
!
!    Output:
!
!      ISGN = -1
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer A(M,N), an array of N columns of vectors 
!    of length M.
!
!    Input, integer I, J, the columns to be compared.
!    I and J must be between 1 and N.
!
!    Output, integer ISGN, the results of the comparison:
!    -1, column I < column J,
!     0, column I = column J,
!    +1, column J < column I.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer a(m,n)
  integer i
  integer isgn
  integer j
  integer k
!
!  Check.
!
  if ( i < 1 .or. n < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Column index I is out of bounds.'
    stop 1
  end if

  if ( j < 1 .or. n < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Column index J is out of bounds.'
    stop 1
  end if

  isgn = 0

  if ( i == j ) then
    return
  end if

  k = 1

  do while ( k <= m )

    if ( a(k,i) < a(k,j) ) then
      isgn = -1
      return
    else if ( a(k,j) < a(k,i) ) then
      isgn = +1
      return
    end if

    k = k + 1

  end do

  return
end
subroutine i4col_find_item ( m, n, a, item, row, col )

!*****************************************************************************80
!
!! I4COL_FIND_ITEM searches a table by columns for a given scalar value.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns in
!    the table.
!
!    Input, integer A(M,N), an array of N columns of vectors
!    of length M.
!
!    Input, integer ITEM, the value to search for.
!
!    Output, integer ROW, COL, the row and column indices
!    of the first occurrence of the value ITEM.  The search
!    is conducted by columns.  If the item is not found, then
!    ROW = COL = -1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer a(m,n)
  integer col
  integer i
  integer item
  integer j
  integer row

  do j = 1, n
    do i = 1, m
      if ( a(i,j) == item ) then
        row = i
        col = j
        return
      end if
    end do
  end do

  row = -1
  col = -1

  return
end
subroutine i4col_find_pair_wrap ( m, n, a, item1, item2, row, col )

!*****************************************************************************80
!
!! I4COL_FIND_PAIR_WRAP searches a table by columns for a pair of items.
!
!  Discussion:
!
!    The items (ITEM1, ITEM2) must occur consecutively.
!    However, wrapping is allowed, that is, if ITEM1 occurs
!    in the last row, and ITEM2 "follows" it in the first row
!    of the same column, a match is declared. 
!
!    If the pair of items is not found, then ROW = COL = -1.  
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns in 
!    the array.
!
!    Input, integer A(M,N), the array to search.
!
!    Input, integer ITEM1, ITEM2, the values to search for.
!
!    Output, integer ROW, COL, the row and column indices
!    of the first occurrence of the value ITEM1 followed immediately
!    by ITEM2.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer a(m,n)
  integer col
  integer i
  integer i2
  integer item1
  integer item2
  integer j
  integer row

  do j = 1, n
    do i = 1, m

      if ( a(i,j) == item1 ) then

        i2 = i + 1

        if ( m < i2 ) then
          i2 = 1
        end if

        if ( a(i2,j) == item2 ) then
          row = i
          col = j
          return
        end if

      end if

    end do
  end do

  row = -1
  col = -1

  return
end
subroutine i4col_sort_a ( m, n, a )

!*****************************************************************************80
!
!! I4COL_SORT_A ascending sorts an integer array of columns.
!
!  Discussion:
!
!    In lexicographic order, the statement "X < Y", applied to two real
!    vectors X and Y of length M, means that there is some index I, with
!    1 <= I <= M, with the property that
!
!      X(J) = Y(J) for J < I,
!    and
!      X(I) < Y(I).
!
!    In other words, the first time they differ, X is smaller.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    25 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows of A, and the length of
!    a vector of data.
!
!    Input, integer N, the number of columns of A.
!
!    Input/output, integer A(M,N).
!    On input, the array of N columns of M-vectors.
!    On output, the columns of A have been sorted in ascending
!    lexicographic order.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer a(m,n)
  integer i
  integer indx
  integer isgn
  integer j

  if ( m <= 0 ) then
    return
  end if

  if ( n <= 1 ) then
    return
  end if
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      call i4col_swap ( m, n, a, i, j )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4col_compare ( m, n, a, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine i4col_sorted_unique_count ( m, n, a, unique_num )

!*****************************************************************************80
!
!! I4COL_SORTED_UNIQUE_COUNT counts unique elements in an I4COL.
!
!  Discussion:
!
!    The columns of the array may be ascending or descending sorted.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer A(M,N), a sorted array, containing
!    N columns of data.
!
!    Output, integer UNIQUE_NUM, the number of unique columns.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer a(m,n)
  integer j1
  integer j2
  integer unique_num

  if ( n <= 0 ) then
    unique_num = 0
    return
  end if

  unique_num = 1
  j1 = 1

  do j2 = 2, n

    if ( any ( a(1:m,j1) /= a(1:m,j2) ) ) then
      unique_num = unique_num + 1
      j1 = j2
    end if

  end do

  return
end
subroutine i4col_swap ( m, n, a, i, j )

!*****************************************************************************80
!
!! I4COL_SWAP swaps columns I and J of a integer array of column data.
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, I = 2, J = 4
!
!      A = (
!        1  2  3  4
!        5  6  7  8
!        9 10 11 12 )
!
!    Output:
!
!      A = (
!        1  4  3  2
!        5  8  7  6
!        9 12 11 10 )
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns in 
!    the array.
!
!    Input/output, integer A(M,N), an array of N columns of 
!    length M.
!
!    Input, integer I, J, the columns to be swapped.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer a(m,n)
  integer col(m)
  integer i
  integer j

  if ( i < 1 .or. n < i .or. j < 1 .or. n < j ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_SWAP - Fatal error!'
    write ( *, '(a)' ) '  I or J is out of bounds.'
    write ( *, '(a,i8)' ) '  I =    ', i
    write ( *, '(a,i8)' ) '  J =    ', j
    write ( *, '(a,i8)' ) '  N =    ', n
    stop 1

  end if

  if ( i == j ) then
    return
  end if

  col(1:m) = a(1:m,i)
  a(1:m,i) = a(1:m,j)
  a(1:m,j) = col(1:m)

  return
end
subroutine i4mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! I4MAT_PRINT prints an integer matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows in A.
!
!    Input, integer N, the number of columns in A.
!
!    Input, integer A(M,N), the matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer a(m,n)
  integer ihi
  integer ilo
  integer jhi
  integer jlo
  character ( len = * )  title

  ilo = 1
  ihi = m
  jlo = 1
  jhi = n

  call i4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

  return
end
subroutine i4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! I4MAT_PRINT_SOME prints some of an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4 values.
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
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer A(M,N), an M by N matrix to be printed.
!
!    Input, integer ILO, JLO, the first row and column to print.
!
!    Input, integer IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: incx = 10
  integer m
  integer n

  integer a(m,n)
  character ( len = 8 )  ctemp(incx)
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
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8)' ) j
    end do

    write ( *, '(''  Col '',10a8)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        write ( ctemp(j2), '(i8)' ) a(i,j)

      end do

      write ( *, '(i5,a,10a8)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine i4mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! I4MAT_TRANSPOSE_PRINT prints an I4MAT, transposed.
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
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer A(M,N), an M by N matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer a(m,n)
  character ( len = * ) title

  call i4mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine i4mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! I4MAT_TRANSPOSE_PRINT_SOME prints some of the transpose of an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4 values.
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
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer A(M,N), an M by N matrix to be printed.
!
!    Input, integer ILO, JLO, the first row and column to print.
!
!    Input, integer IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: incx = 10
  integer m
  integer n

  integer a(m,n)
  character ( len = 8 )  ctemp(incx)
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
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(i8)' ) i
    end do

    write ( *, '(''  Row '',10a8)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc

        i = i2lo - 1 + i2

        write ( ctemp(i2), '(i8)' ) a(i,j)

      end do

      write ( *, '(i5,a,10a8)' ) j, ':', ( ctemp(i), i = 1, inc )

    end do

  end do

  return
end
subroutine i4row_compare ( m, n, a, i, j, isgn )

!*****************************************************************************80
!
!! I4ROW_COMPARE compares two rows of a integer array.
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, I = 2, J = 3
!
!      A = (
!        1  2  3  4
!        5  6  7  8
!        9 10 11 12 )
!
!    Output:
!
!      ISGN = -1
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer A(M,N), an array of M rows of vectors of 
!    length N.
!
!    Input, integer I, J, the rows to be compared.
!    I and J must be between 1 and M.
!
!    Output, integer ISGN, the results of the comparison:
!    -1, row I < row J,
!     0, row I = row J,
!    +1, row J < row I.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer a(m,n)
  integer i
  integer isgn
  integer j
  integer k
!
!  Check that I and J are legal.
!
  if ( i < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4ROW_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Row index I is less than 1.'
    write ( *, '(a,i8)' ) '  I = ', i
    stop 1
  else if ( m < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4ROW_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Row index I is out of bounds.'
    write ( *, '(a,i8)' ) '  I = ', i
    write ( *, '(a,i8)' ) '  Maximum legal value is M = ', m
    stop 1
  end if

  if ( j < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4ROW_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Row index J is less than 1.'
    write ( *, '(a,i8)' ) '  J = ', j
    stop 1
  else if ( m < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4ROW_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Row index J is out of bounds.'
    write ( *, '(a,i8)' ) '  J = ', j 
    write ( *, '(a,i8)' ) '  Maximum legal value is M = ', m
    stop 1
  end if

  isgn = 0

  if ( i == j ) then
    return
  end if

  k = 1

  do while ( k <= n )

    if ( a(i,k) < a(j,k) ) then
      isgn = -1
      return
    else if ( a(j,k) < a(i,k) ) then
      isgn = +1
      return
    end if

    k = k + 1

  end do

  return
end
subroutine i4row_sort_a ( m, n, a )

!*****************************************************************************80
!
!! I4ROW_SORT_A ascending sorts the rows of an integer array.
!
!  Discussion:
!
!    In lexicographic order, the statement "X < Y", applied to two
!    vectors X and Y of length M, means that there is some index I, with
!    1 <= I <= M, with the property that
!
!      X(J) = Y(J) for J < I,
!    and
!      X(I) < Y(I).
!
!    In other words, X is less than Y if, at the first index where they
!    differ, the X value is less than the Y value.
!
!  Example:
!
!    Input:
!
!      M = 5, N = 3
!
!      A =
!        3  2  1
!        2  4  3
!        3  1  8
!        2  4  2
!        1  9  9
!
!    Output:
!
!      A =
!        1  9  9
!        2  4  2
!        2  4  3
!        3  1  8
!        3  2  1
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows of A.
!
!    Input, integer N, the number of columns of A.
!
!    Input/output, integer A(M,N).
!    On input, the array of M rows of N-vectors.
!    On output, the rows of A have been sorted in ascending
!    lexicographic order.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer a(m,n)
  integer i
  integer indx
  integer isgn
  integer j
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( m, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      call i4row_swap ( m, n, a, i, j )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4row_compare ( m, n, a, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine i4row_sorted_unique_count ( m, n, a, unique_num )

!*****************************************************************************80
!
!! I4ROW_SORTED_UNIQUE_COUNT counts unique elements in an IROW array.
!
!  Discussion:
!
!    The rows of the array may be ascending or descending sorted.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer A(M,N), a sorted array, containing
!    M rows of data.
!
!    Output, integer UNIQUE_NUM, the number of unique rows.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer a(m,n)
  integer i1
  integer i2
  integer unique_num

  if ( n <= 0 ) then
    unique_num = 0
    return
  end if

  unique_num = 1
  i1 = 1

  do i2 = 2, m

    if ( any ( a(i1,1:n) /= a(i2,1:n) ) ) then
      unique_num = unique_num + 1
      i1 = i2
    end if

  end do

  return
end
subroutine i4row_swap ( m, n, a, irow1, irow2 )

!*****************************************************************************80
!
!! I4ROW_SWAP swaps two rows of an integer array.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input/output, integer A(M,N), an array of data.
!
!    Input, integer IROW1, IROW2, the two rows to swap.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer a(m,n)
  integer irow1
  integer irow2
  integer row(n)
!
!  Check.
!
  if ( irow1 < 1 .or. m < irow1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4ROW_SWAP - Fatal error!'
    write ( *, '(a)' ) '  IROW1 is out of range.'
    stop 1
  end if

  if ( irow2 < 1 .or. m < irow2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4ROW_SWAP - Fatal error!'
    write ( *, '(a)' ) '  IROW2 is out of range.'
    stop 1
  end if

  if ( irow1 == irow2 ) then
    return
  end if

  row(1:n) = a(irow1,1:n)
  a(irow1,1:n) = a(irow2,1:n)
  a(irow2,1:n) = row(1:n)

  return
end
subroutine i4vec_heap_d ( n, a )

!*****************************************************************************80
!
!! I4VEC_HEAP_D reorders an array of integers into a descending heap.
!
!  Discussion:
!
!    A descending heap is an array A with the property that, for every index J,
!    A(2*J) <= A(J) and A(2*J+1) <= A(J), (as long as the indices
!    2*J and 2*J+1 are legal).
!
!  Diagram:
!
!                  A(1)
!                /      \
!            A(2)         A(3)
!          /     \        /  \
!      A(4)       A(5)  A(6) A(7)
!      /  \       /   \
!    A(8) A(9) A(10) A(11)
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer N, the size of the input array.
!
!    Input/output, integer A(N).
!    On input, an unsorted array.
!    On output, the array has been reordered into a heap.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer a(n)
  integer i
  integer ifree
  integer key
  integer m
!
!  Only nodes N/2 down to 1 can be "parent" nodes.
!
  do i = n/2, 1, -1
!
!  Copy the value out of the parent node.
!  Position IFREE is now "open".
!
    key = a(i)
    ifree = i

    do
!
!  Positions 2*IFREE and 2*IFREE + 1 are the descendants of position
!  IFREE.  (One or both may not exist because they exceed N.)
!
      m = 2 * ifree
!
!  Does the first position exist?
!
      if ( n < m ) then
        exit
      end if
!
!  Does the second position exist?
!
      if ( m + 1 <= n ) then
!
!  If both positions exist, take the larger of the two values,
!  and update M if necessary.
!
        if ( a(m) < a(m+1) ) then
          m = m + 1
        end if

      end if
!
!  If the large descendant is larger than KEY, move it up,
!  and update IFREE, the location of the free position, and
!  consider the descendants of THIS position.
!
      if ( a(m) <= key ) then
        exit
      end if

      a(ifree) = a(m)
      ifree = m

    end do
!
!  Once there is no more shifting to do, KEY moves into the free spot IFREE.
!
    a(ifree) = key

  end do

  return
end
subroutine i4vec_indicator ( n, a )

!*****************************************************************************80
!
!! I4VEC_INDICATOR sets an integer vector to the indicator vector.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements of A.
!
!    Output, integer A(N), the array to be initialized.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer a(n)
  integer i

  do i = 1, n
    a(i) = i
  end do

  return
end
function i4vec_lcm ( n, v )

!*****************************************************************************80
!
!! I4VEC_LCM returns the least common multiple of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    The value LCM returned has the property that it is the smallest integer
!    which is evenly divisible by every element of V.
!
!    The entries in V may be negative.
!
!    If any entry of V is 0, then LCM is 0.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 July 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of V.
!
!    Input, integer V(N), the vector.
!
!    Output, integer I4VEC_LCM, the least common multiple of V.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer i4_lcm
  integer i4vec_lcm
  integer lcm
  integer v(n)

  lcm = 1

  do i = 1, n

    if ( v(i) == 0 ) then
      lcm = 0
      i4vec_lcm = lcm
      return
    end if

    lcm = i4_lcm ( lcm, v(i) )

  end do

  i4vec_lcm = lcm

  return
end
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 May 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, integer A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer a(n)
  integer i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,a,2x,i12)' ) i, ':', a(i)
  end do

  return
end
subroutine i4vec_sort_heap_a ( n, a )

!*****************************************************************************80
!
!! I4VEC_SORT_HEAP_A ascending sorts an integer array using heap sort.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input/output, integer A(N).
!    On input, the array to be sorted;
!    On output, the array has been sorted.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer a(n)
  integer n1

  if ( n <= 1 ) then
    return
  end if
!
!  1: Put A into descending heap form.
!
  call i4vec_heap_d ( n, a )
!
!  2: Sort A.
!
!  The largest object in the heap is in A(1).
!  Move it to position A(N).
!
  call i4_swap ( a(1), a(n) )
!
!  Consider the diminished heap of size N1.
!
  do n1 = n - 1, 2, -1
!
!  Restore the heap structure of A(1) through A(N1).
!
    call i4vec_heap_d ( n1, a )
!
!  Take the largest object from A(1) and move it to A(N1).
!
    call i4_swap ( a(1), a(n1) )

  end do

  return
end
subroutine i4vec_sorted_unique ( n, a, unique_num )

!*****************************************************************************80
!
!! I4VEC_SORTED_UNIQUE gets the unique elements in a sorted integer array.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements in A.
!
!    Input/output, integer A(N).  On input, the sorted
!    integer array.  On output, the unique elements in A.
!
!    Output, integer UNIQUE_NUM, the number of unique elements
!    in A.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer a(n)
  integer itest
  integer unique_num

  unique_num = 0

  if ( n <= 0 ) then
    return
  end if

  unique_num = 1

  do itest = 2, n

    if ( a(itest) /= a(unique_num) ) then
      unique_num = unique_num + 1
      a(unique_num) = a(itest)
    end if

  end do

  return
end
subroutine i4vec_uniform_ab ( n, a, b, x )

!*****************************************************************************80
!
!! I4VEC_UNIFORM_AB returns a scaled pseudorandom I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    The pseudorandom numbers should be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of the vector.
!
!    Input, integer A, B, the limits of the interval.
!
!    Output, integer X(N), a vector of numbers between A and B.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer a
  integer b
  integer i
  real r
  integer value
  integer x(n)

  do i = 1, n

    call random_number ( harvest = r )
!
!  Scale R to lie between A-0.5 and B+0.5.
!
    r = ( 1.0E+00 - r ) * ( real ( min ( a, b ) ) - 0.5E+00 ) & 
      +             r   * ( real ( max ( a, b ) ) + 0.5E+00 )
!
!  Use rounding to convert R to an integer between A and B.
!
    value = nint ( r )

    value = max ( value, min ( a, b ) )
    value = min ( value, max ( a, b ) )

    x(i) = value

  end do

  return
end
subroutine i4vec2_compare ( n, a1, a2, i, j, isgn )

!*****************************************************************************80
!
!! I4VEC2_COMPARE compares pairs of integers stored in two vectors.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data items.
!
!    Input, integer A1(N), A2(N), contain the two components 
!    of each item.
!
!    Input, integer I, J, the items to be compared.
!
!    Output, integer ISGN, the results of the comparison:
!    -1, item I < item J,
!     0, item I = item J,
!    +1, item J < item I.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer a1(n)
  integer a2(n)
  integer i
  integer isgn
  integer j

  isgn = 0

       if ( a1(i) < a1(j) ) then

    isgn = -1

  else if ( a1(i) == a1(j) ) then

         if ( a2(i) < a2(j) ) then
      isgn = -1
    else if ( a2(i) < a2(j) ) then
      isgn = 0
    else if ( a2(j) < a2(i) ) then
      isgn = +1
    end if

  else if ( a1(j) < a1(i) ) then

    isgn = +1

  end if

  return
end
subroutine i4vec2_sort_a ( n, a1, a2 )

!*****************************************************************************80
!
!! I4VEC2_SORT_A ascending sorts a vector of pairs of integers.
!
!  Discussion:
!
!    Each item to be sorted is a pair of integers (I,J), with the I
!    and J values stored in separate vectors A1 and A2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of items of data.
!
!    Input/output, integer A1(N), A2(N), the data to be sorted..
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer a1(n)
  integer a2(n)
  integer i
  integer indx
  integer isgn
  integer j
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      call i4_swap ( a1(i), a1(j) )
      call i4_swap ( a2(i), a2(j) )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4vec2_compare ( n, a1, a2, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine i4vec2_sorted_unique ( n, a1, a2, unique_num )

!*****************************************************************************80
!
!! I4VEC2_SORTED_UNIQUE gets the unique elements in a sorted I4VEC2.
!
!  Discussion:
!
!    Item I is stored as the pair A1(I), A2(I).
!
!    The items must have been sorted, or at least it must be the
!    case that equal items are stored in adjacent vector locations.
!
!    If the items were not sorted, then this routine will only
!    replace a string of equal values by a single representative.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of items.
!
!    Input/output, integer A1(N), A2(N).
!    On input, the array of N items.
!    On output, an array of unique items.
!
!    Output, integer UNIQUE_NUM, the number of unique items.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer a1(n)
  integer a2(n)
  integer itest
  integer unique_num

  unique_num = 0

  if ( n <= 0 ) then
    return
  end if

  unique_num = 1

  do itest = 2, n

    if ( a1(itest) /= a1(unique_num) .or. a2(itest) /= a2(unique_num) ) then

      unique_num = unique_num + 1

      a1(unique_num) = a1(itest)
      a2(unique_num) = a2(itest)

    end if

  end do

  return
end
subroutine icos_shape ( point_num, edge_num, face_num, face_order_max, &
  point_coord, edge_point, face_order, face_point )

!*****************************************************************************80
!
!! ICOS_SHAPE describes an icosahedron.
!
!  Discussion:
!
!    The input data required for this routine can be retrieved from
!    ICOS_SIZE.
!
!    The vertices lie on the unit sphere.
!
!    The dual of an icosahedron is a dodecahedron.
!
!    The data has been rearranged from a previous assignment.  
!    The STRIPACK program refuses to triangulate data if the first
!    three nodes are "collinear" on the sphere.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer POINT_NUM, the number of points (12).
!
!    Input, integer EDGE_NUM, the number of edges (30).
!
!    Input, integer FACE_NUM, the number of faces (20).
!
!    Input, integer FACE_ORDER_MAX, the maximum number of 
!    vertices per face (3).
!
!    Output, real ( kind = rk ) POINT_COORD(3,POINT_NUM), the points.
!
!    Output, integer EDGE_POINT(2,EDGE_NUM), the points that 
!    make up each edge, listed in ascending order of their indexes.
!
!    Output, integer FACE_ORDER(FACE_NUM), the number of vertices
!    per face.
!
!    Output, integer FACE_POINT(FACE_ORDER_MAX,FACE_NUM); 
!    FACE_POINT(I,J) is the index of the I-th point in the J-th face.  The
!    points are listed in the counter clockwise direction defined
!    by the outward normal at the face.  The nodes of each face are ordered 
!    so that the lowest index occurs first.  The faces are then sorted by
!    nodes.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer edge_num
  integer, parameter :: edge_order = 2
  integer face_num
  integer face_order_max
  integer, parameter :: dim_num = 3
  integer point_num

  real ( kind = rk ) a
  real ( kind = rk ) b
  integer edge_point(edge_order,edge_num)
  integer face_order(face_num)
  integer face_point(face_order_max,face_num)
  real ( kind = rk ) phi
  real ( kind = rk ) point_coord(dim_num,point_num)
  real ( kind = rk ) z
!
!  Set the point coordinates.
!
  phi = 0.5D+00 * ( sqrt ( 5.0D+00 ) + 1.0D+00 )

  a = phi / sqrt ( 1.0D+00 + phi * phi )
  b = 1.0D+00 / sqrt ( 1.0D+00 + phi * phi )
  z = 0.0D+00
!
!  A*A + B*B + Z*Z = 1.
!
  point_coord(1:dim_num,1:point_num) = reshape ( (/ &
      a,  b,  z, &
      a, -b,  z, &
      b,  z,  a, &
      b,  z, -a, &
      z,  a,  b, &
      z,  a, -b, &
      z, -a,  b, &
      z, -a, -b, &
     -b,  z,  a, &
     -b,  z, -a, &
     -a,  b,  z, &
     -a, -b,  z /), (/ dim_num, point_num /) )
!
!  Set the edges.
!
  edge_point(1:edge_order,1:edge_num) = reshape ( (/ &
     1,  2, &
     1,  3, &
     1,  4, &
     1,  5, &
     1,  6, &
     2,  3, &
     2,  4, &
     2,  7, &
     2,  8, &
     3,  5, &
     3,  7, &
     3,  9, &
     4,  6, &
     4,  8, &
     4, 10, &
     5,  6, &
     5,  9, &
     5, 11, &
     6, 10, &
     6, 11, &
     7,  8, &
     7,  9, &
     7, 12, &
     8, 10, &
     8, 12, &
     9, 11, &
     9, 12, &
    10, 11, &
    10, 12, &
    11, 12 /), (/ edge_order, edge_num /) )
!
!  Set the face orders.
!
  face_order(1:face_num) = (/ &
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3 /)
!
!  Set the faces.
!
  face_point(1:face_order_max,1:face_num) = reshape ( (/ &
     1,  2,  4, &
     1,  3,  2, &
     1,  4,  6, &
     1,  5,  3, &
     1,  6,  5, &
     2,  3,  7, &
     2,  7,  8, &
     2,  8,  4, &
     3,  5,  9, &
     3,  9,  7, &
     4,  8, 10, &
     4, 10,  6, &
     5,  6, 11, &
     5, 11,  9, &
     6, 10, 11, &
     7,  9, 12, &
     7, 12,  8, &
     8, 12, 10, &
     9, 11, 12, &
    10, 12, 11 /), (/ face_order_max, face_num /) )

  return
end
subroutine icos_size ( point_num, edge_num, face_num, face_order_max )

!*****************************************************************************80
!
!! ICOS_SIZE gives "sizes" for an icosahedron.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer POINT_NUM, the number of points.
!
!    Output, integer EDGE_NUM, the number of edges.
!
!    Output, integer FACE_NUM, the number of faces.
!
!    Output, integer FACE_ORDER_MAX, the maximum order of any face.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer edge_num
  integer face_num
  integer face_order_max
  integer point_num

  point_num = 12
  edge_num = 30
  face_num = 20
  face_order_max = 3

  return
end
function line_exp_is_degenerate_nd ( dim_num, p1, p2 )

!*****************************************************************************80
!
!! line_exp_is_degenerate_nd finds if an explicit line is degenerate in ND.
!
!  Discussion:
!
!    The explicit form of a line in ND is:
!
!      the line through the points P1 and P2.
!
!    An explicit line is degenerate if the two defining points are equal.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = rk ) P1(DIM_NUM), P2(DIM_NUM), two points on the line.
!
!    Output, logical LINE_EXP_IS_DEGENERATE_ND, is TRUE if the line
!    is degenerate.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num

  logical line_exp_is_degenerate_nd
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)

  line_exp_is_degenerate_nd = ( all ( p1(1:dim_num) == p2(1:dim_num) ) )

  return
end
subroutine line_exp_normal_2d ( p1, p2, normal )

!*****************************************************************************80
!
!! LINE_EXP_NORMAL_2D computes a unit normal vector to a line in 2D.
!
!  Discussion:
!
!    The explicit form of a line in 2D is:
!
!      the line through the points P1 and P2.
!
!    The sign of the normal vector N is chosen so that the normal vector
!    points "to the left" of the direction of the line.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 October 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(2), P2(2), two distinct points on the line.
!
!    Output, real ( kind = rk ) NORMAL(2), a unit normal vector to the line.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  logical line_exp_is_degenerate_nd
  real ( kind = rk ) norm
  real ( kind = rk ) normal(dim_num)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)

  if ( line_exp_is_degenerate_nd ( dim_num, p1, p2 ) ) then
    normal(1:dim_num) = sqrt ( 2.0D+00 )
    return
  end if

  norm = sqrt ( ( p2(1) - p1(1) )**2 + ( p2(2) - p1(2) )**2 )

  normal(1) = - ( p2(2) - p1(2) ) / norm
  normal(2) =   ( p2(1) - p1(1) ) / norm

  return
end
subroutine line_exp_perp_2d ( p1, p2, p3, p4, flag )

!*****************************************************************************80
!
!! LINE_EXP_PERP_2D computes a line perpendicular to a line and through a point.
!
!  Discussion:
!
!    The explicit form of a line in 2D is:
!
!      the line through the points P1 and P2.
!
!    The input point P3 should NOT lie on the line (P1,P2).  If it
!    does, then the output value P4 will equal P3.
!
!    P1-----P4-----------P2
!            |
!            |
!           P3
!
!    P4 is also the nearest point on the line (P1,P2) to the point P3.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    31 July 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(2), P2(2), two points on the line.
!
!    Input, real ( kind = rk ) P3(2), a point (presumably not on the 
!    line (P1,P2)), through which the perpendicular must pass.
!
!    Output, real ( kind = rk ) P4(2), a point on the line (P1,P2),
!    such that the line (P3,P4) is perpendicular to the line (P1,P2).
!
!    Output, logical FLAG, is TRUE if the value could 
!    not be computed.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) bot
  logical flag
  logical line_exp_is_degenerate_nd
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) p3(dim_num)
  real ( kind = rk ) p4(dim_num)
  real ( kind = rk ) r8_huge
  real ( kind = rk ) t

  flag = .false.

  if ( line_exp_is_degenerate_nd ( dim_num, p1, p2 ) ) then
    flag = .true.
    p4(1:2) = r8_huge ( )
    return
  end if

  bot = sum ( ( p2(1:dim_num) - p1(1:dim_num) )**2 )
!
!  (P3-P1) dot (P2-P1) = Norm(P3-P1) * Norm(P2-P1) * Cos(Theta).
!
!  (P3-P1) dot (P2-P1) / Norm(P3-P1)^2 = normalized coordinate T
!  of the projection of (P3-P1) onto (P2-P1).
!
  t = sum ( ( p1(1:dim_num) - p3(1:dim_num) ) &
          * ( p1(1:dim_num) - p2(1:dim_num) ) ) / bot

  p4(1:dim_num) = p1(1:dim_num) + t * ( p2(1:dim_num) - p1(1:dim_num) )

  return
end
subroutine line_exp_point_dist_2d ( p1, p2, p, dist )

!*****************************************************************************80
!
!! LINE_EXP_POINT_DIST_2D: distance ( explicit line, point ) in 2D.
!
!  Discussion:
!
!    The explicit form of a line in 2D is:
!
!      the line through the points P1 and P2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(2), P2(2), two points on the line.
!
!    Input, real ( kind = rk ) P(2), the point whose distance from the line is
!    to be measured.
!
!    Output, real ( kind = rk ) DIST, the distance from the point to the line.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) bot
  real ( kind = rk ) dist
  real ( kind = rk ) dot
  logical line_exp_is_degenerate_nd
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) pn(dim_num)
  real ( kind = rk ) t

  if ( line_exp_is_degenerate_nd ( dim_num, p1, p2 ) ) then

    pn(1:dim_num) = p1(1:dim_num)
!
!  (P-P1) dot (P2-P1) = Norm(P-P1) * Norm(P2-P1) * Cos(Theta).
!
!  (P-P1) dot (P2-P1) / Norm(P-P1)^2 = normalized coordinate T
!  of the projection of (P-P1) onto (P2-P1).
!
  else

    dot = sum ( ( p(1:dim_num) - p1(1:dim_num) ) &
              * ( p2(1:dim_num) - p1(1:dim_num) ) )

    bot = sum ( ( p2(1:dim_num) - p1(1:dim_num) )**2 )

    t = dot / bot

    pn(1:dim_num) = p1(1:dim_num) + t * ( p2(1:dim_num) - p1(1:dim_num) )

  end if

  dist = sqrt ( sum ( ( p(1:dim_num) - pn(1:dim_num) )**2 ) )

  return
end
subroutine line_exp_point_dist_3d ( p1, p2, p, dist )

!*****************************************************************************80
!
!! LINE_EXP_POINT_DIST_3D: distance ( explicit line, point ) in 3D.
!
!  Discussion:
!
!    The explicit form of a line in 3D is:
!
!      the line through the points P1 and P2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(3), P2(3), two points on the line.
!
!    Input, real ( kind = rk ) P(3), the point whose distance from the line is
!    to be measured.
!
!    Output, real ( kind = rk ) DIST, the distance from the point to the line.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) bot
  real ( kind = rk ) dist
  logical line_exp_is_degenerate_nd
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) pn(dim_num)
  real ( kind = rk ) t

  if ( line_exp_is_degenerate_nd ( dim_num, p1, p2 ) ) then

    pn(1:dim_num) = p1(1:dim_num)
!
!  (P-P1) dot (P2-P1) = Norm(P-P1) * Norm(P2-P1) * Cos(Theta).
!
!  (P-P1) dot (P2-P1) / Norm(P-P1)^2 = normalized coordinate T
!  of the projection of (P-P1) onto (P2-P1).
!
  else

    bot = sum ( ( p2(1:dim_num) - p1(1:dim_num) )**2 )

    t = sum ( ( p(1:dim_num) - p1(1:dim_num) ) &
            * ( p2(1:dim_num) - p1(1:dim_num) ) ) / bot

    pn(1:dim_num) = p1(1:dim_num) + t * ( p2(1:dim_num) - p1(1:dim_num) )

  end if
!
!  Now compute the distance between the projection point and P.
!
  dist = sqrt ( sum ( ( p(1:dim_num) - pn(1:dim_num) )**2 ) )

  return
end
subroutine line_exp_point_dist_signed_2d ( p1, p2, p, dist_signed )

!*****************************************************************************80
!
!! LINE_EXP_POINT_DIST_SIGNED_2D: signed distance ( exp line, point ) in 2D.
!
!  Discussion:
!
!    The explicit form of a line in 2D is:
!
!      the line through the points P1 and P2.
!
!    The signed distance has two interesting properties:
!
!    *  The absolute value of the signed distance is the
!        usual (Euclidean) distance.
!
!    *  Points with signed distance 0 lie on the line,
!       points with a negative signed distance lie on one side
!         of the line,
!       points with a positive signed distance lie on the
!         other side of the line.
!
!    Assuming that C is nonnegative, then if a point is a positive
!    distance away from the line, it is on the same side of the
!    line as the point (0,0), and if it is a negative distance
!    from the line, it is on the opposite side from (0,0).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(2), P2(2), two points on the line.
!
!    Input, real ( kind = rk ) P(2), the point whose signed distance is desired.
!
!    Output, real ( kind = rk ) DIST_SIGNED, the signed distance from the
!    point to the line.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) dist_signed
  logical line_exp_is_degenerate_nd
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
!
!  If the explicit line degenerates to a point, the computation is easy.
!
  if ( line_exp_is_degenerate_nd ( dim_num, p1, p2 ) ) then

    dist_signed = sqrt ( sum ( ( p1(1:dim_num) - p(1:dim_num) )**2 ) )
!
!  Convert the explicit line to the implicit form A * P(1) + B * P(2) + C = 0.
!  This makes the computation of the signed distance to (X,Y) easy.
!
  else

    a = p2(2) - p1(2)
    b = p1(1) - p2(1)
    c = p2(1) * p1(2) - p1(1) * p2(2)

    dist_signed = ( a * p(1) + b * p(2) + c ) / sqrt ( a * a + b * b )

  end if

  return
end
subroutine line_exp_point_near_2d ( p1, p2, p, pn, dist, t )

!*****************************************************************************80
!
!! LINE_EXP_POINT_NEAR_2D: point on an explicit line nearest a point in 2D.
!
!  Discussion:
!
!    The explicit form of a line in 2D is:
!
!      the line through the points P1 and P2.
!
!    The nearest point PN = (XN,YN) has the form:
!
!      PN = (1-T) * P1 + T * P2.
!
!    If T is less than 0, PN is furthest from P2.
!    If T is between 0 and 1, PN is between P1 and P2.
!    If T is greater than 1, PN is furthest from P1.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(2), P2(2), two points on the line.
!
!    Input, real ( kind = rk ) P(2), the point whose nearest neighbor on the
!    line is to be determined.
!
!    Output, real ( kind = rk ) PN(2), the nearest point on the line to P.
!
!    Output, real ( kind = rk ) DIST, the distance from the point to the line.
!
!    Output, real ( kind = rk ) T, the relative position of the point
!    PN to the points P1 and P2.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) bot
  real ( kind = rk ) dist
  logical line_exp_is_degenerate_nd
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) pn(dim_num)
  real ( kind = rk ) t

  if ( line_exp_is_degenerate_nd ( dim_num, p1, p2 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LINE_POINT_NEAR_2D - Fatal error!'
    write ( *, '(a)' ) '  The line is degenerate.'
    stop 1
  end if
!
!  (P-P1) dot (P2-P1) = Norm(P-P1) * Norm(P2-P1) * Cos(Theta).
!
!  (P-P1) dot (P2-P1) / Norm(P-P1)^2 = normalized coordinate T
!  of the projection of (P-P1) onto (P2-P1).
!
  bot = sum ( ( p2(1:dim_num) - p1(1:dim_num) )**2 )

  t = sum ( ( p1(1:dim_num) - p(1:dim_num) ) &
          * ( p1(1:dim_num) - p2(1:dim_num) ) ) / bot

  pn(1:dim_num) = p1(1:dim_num) + t * ( p2(1:dim_num) - p1(1:dim_num) )

  dist = sqrt ( sum ( ( pn(1:dim_num) - p(1:dim_num) )**2 ) )

  return
end
subroutine line_exp_point_near_3d ( p1, p2, p, pn, dist, t )

!*****************************************************************************80
!
!! LINE_EXP_POINT_NEAR_3D: nearest point on explicit line to point in 3D.
!
!  Discussion:
!
!    The explicit form of a line in 3D is:
!
!      the line through the points P1 and P2.
!
!    The nearest point PN has the form:
!
!      PN = ( 1 - T ) * P1 + T * P2.
!
!    If T is less than 0, PN is furthest away from P2.
!    If T is between 0 and 1, PN is between P1 and P2.
!    If T is greater than 1, PN is furthest away from P1.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(3), P2(3), two points on the line.
!
!    Input, real ( kind = rk ) P(3), the point whose nearest neighbor on
!    the line is to be determined.
!
!    Output, real ( kind = rk ) PN(3), the point which is the nearest
!    point on the line to P.
!
!    Output, real ( kind = rk ) DIST, the distance from the point to the 
!    nearest point on the line.
!
!    Output, real ( kind = rk ) T, the relative position of the point
!    PN to P1 and P2.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) bot
  real ( kind = rk ) dist
  logical line_exp_is_degenerate_nd
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) pn(dim_num)
  real ( kind = rk ) t

  if ( line_exp_is_degenerate_nd ( dim_num, p1, p2 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LINE_EXP_POINT_NEAR_3D - Fatal error!'
    write ( *, '(a)' ) '  The line is degenerate.'
    stop 1
  end if
!
!  (P-P1) dot (P2-P1) = Norm(P-P1) * Norm(P2-P1) * Cos(Theta).
!
!  (P-P1) dot (P2-P1) / Norm(P-P1)^2 = normalized coordinate T
!  of the projection of (P-P1) onto (P2-P1).
!
  bot = sum ( ( p2(1:dim_num) - p1(1:dim_num) )**2 )

  t = sum ( ( p(1:dim_num) - p1(1:dim_num) ) &
          * ( p2(1:dim_num) - p1(1:dim_num) ) ) / bot
!
!  Now compute the location of the projection point.
!
  pn(1:dim_num) = p1(1:dim_num) + t * ( p2(1:dim_num) - p1(1:dim_num) )
!
!  Now compute the distance between the projection point and P.
!
  dist = sqrt ( sum ( ( pn(1:dim_num) - p(1:dim_num) )**2 ) )

  return
end
subroutine line_exp_point_projection_2d ( p1, p2, p3, t )

!*****************************************************************************80
!
!! line_exp_point_projection_2d projects a point onto an explicit line.
!
!  Discussion:
!
!    The explicit form of a line in 2D is:
!
!      P(t) = (1-t)*P1 + t * P2
!
!    The projection of a point P3 is the nearest point on that line.
!    We can characterize that point by T, its "coordinate".
!
!    If P3 is actually on the line, then we have the following useful facts:
!    t < 0: P3 is "to the left of" P1.
!    t = 0: P3 = P1
!    0 < t < 1: P3 is between P1 and P2.
!    t = 1: P3 = P2.
!    1 < t: P3 is "to the right of" P2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 July 2020
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) P1(2), P2(2), two points on the line.
!
!    real ( kind = rk ) P3(2), the point whose projection is needed.
!
!  Output:
!
!    real ( kind = rk ) T, the projection coordinate of P3.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) d(2)
  real ( kind = rk ) dn
  real ( kind = rk ) e(2)
  real ( kind = rk ) p1(2)
  real ( kind = rk ) p2(2)
  real ( kind = rk ) p3(2)
  real ( kind = rk ) t

  d(1) = p2(1) - p1(1)
  d(2) = p2(2) - p1(2)
  dn = sqrt ( d(1)**2 + d(2)**2 )
  d(1) = d(1) / dn
  d(2) = d(2) / dn

  e(1) = p3(1) - p1(1)
  e(2) = p3(2) - p1(2)

  t = d(1) * e(1) + e(2) * d(2)

  return
end
subroutine line_exp_print_2d ( p1, p2, label )

!*****************************************************************************80
!
!! line_exp_print_2d: prints an explicit line in 2D.
!
!  Discussion:
!
!    The explicit form of a line in 3D is:
!
!      the line through the distinct points P1 and P2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 July 2020
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) P1(2), P2(2), two points on the line.
!
!    character ( len = * ) label: a label for the line.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = * ) label
  real ( kind = rk ) p1(2)
  real ( kind = rk ) p2(2)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) label
  write ( *, '(a,2f8.4)' ) '    P1 = ', p1(1:2)
  write ( *, '(a,2f8.4)' ) '    P2 = ', p2(1:2)

  return
end
subroutine line_exp2imp_2d ( p1, p2, a, b, c )

!*****************************************************************************80
!
!! LINE_EXP2IMP_2D converts an explicit line to implicit form in 2D.
!
!  Discussion:
!
!    The explicit form of a line in 2D is:
!
!      the line through the points P1 and P2.
!
!    The implicit form of a line in 2D is:
!
!      A * X + B * Y + C = 0
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) P1(2), P2(2), two points on the line.
!
!  Output:
!
!    real ( kind = rk ) A, B, C, the implicit form of the line.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  logical line_exp_is_degenerate_nd
  real ( kind = rk ) norm
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
!
!  Take care of degenerate cases.
!
  if ( line_exp_is_degenerate_nd ( dim_num, p1, p2 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LINE_EXP2IMP_2D - Warning!'
    write ( *, '(a)' ) '  The line is degenerate.'
  end if

  a = p2(2) - p1(2)
  b = p1(1) - p2(1)
  c = p2(1) * p1(2) - p1(1) * p2(2)

  norm = a * a + b * b + c * c

  if ( 0.0D+00 < norm ) then
    a = a / norm
    b = b / norm
    c = c / norm
  end if

  if ( a < 0.0D+00 ) then
    a = -a
    b = -b
    c = -c
  end if

  return
end
subroutine line_exp2par_2d ( p1, p2, f, g, x0, y0 )

!*****************************************************************************80
!
!! LINE_EXP2PAR_2D converts a line from explicit to parametric form in 2D.
!
!  Discussion:
!
!    The explicit form of a line in 2D is:
!
!      the line through the points P1 and P2.
!
!    The parametric form of a line in 2D is:
!
!      X = X0 + F * T
!      Y = Y0 + G * T
!
!    We normalize by always choosing F^2 + G^2 = 1, and F nonnegative.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) P1(2), P2(2), two points on the line.
!
!  Output:
!
!    real ( kind = rk ) F, G, X0, Y0, the parametric parameters
!    of the line.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) f
  real ( kind = rk ) g
  logical line_exp_is_degenerate_nd
  real ( kind = rk ) norm
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) x0
  real ( kind = rk ) y0

  if ( line_exp_is_degenerate_nd ( dim_num, p1, p2 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LINE_EXP2PAR_2D - Warning!'
    write ( *, '(a)' ) '  The line is degenerate.'
  end if

  x0 = p1(1)
  y0 = p1(2)

  f = p2(1) - p1(1)
  g = p2(2) - p1(2)

  norm = sqrt ( f * f + g * g )

  if ( norm /= 0.0D+00 ) then
    f = f / norm
    g = g / norm
  end if

  if ( f < 0.0D+00 ) then
    f = -f
    g = -g
  end if

  return
end
subroutine line_exp2par_3d ( p1, p2, f, g, h, x0, y0, z0 )

!*****************************************************************************80
!
!! LINE_EXP2PAR_3D converts a line from explicit to parametric form in 3D.
!
!  Discussion:
!
!    The explicit form of a line in 3D is:
!
!      the line through the points P1 and P2.
!
!    The parametric form of a line in 3D is:
!
!      X = X0 + F * T
!      Y = Y0 + G * T
!      Z = Z0 + H * T
!
!    We normalize by always choosing F^2 + G^2 + H^2 = 1, and F nonnegative.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(3), P2(3), two points on the line.
!
!    Output, real ( kind = rk ) F, G, H, X0, Y0, Z0, the parametric parameters
!    of the line.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) f
  real ( kind = rk ) g
  real ( kind = rk ) h
  logical line_exp_is_degenerate_nd
  real ( kind = rk ) norm
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) x0
  real ( kind = rk ) y0
  real ( kind = rk ) z0

  if ( line_exp_is_degenerate_nd ( dim_num, p1, p2 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LINE_EXP2PAR_3D - Warning!'
    write ( *, '(a)' ) '  The line is degenerate.'
  end if

  x0 = p1(1)
  y0 = p1(2)
  z0 = p1(3)

  f = p2(1) - p1(1)
  g = p2(2) - p1(2)
  h = p2(3) - p1(3)

  norm = sqrt ( f * f + g * g + h * h )

  if ( norm /= 0.0D+00 ) then
    f = f / norm
    g = g / norm
    h = h / norm
  end if

  if ( f < 0.0D+00 ) then
    f = -f
    g = -g
    h = -h
  end if

  return
end
function line_imp_is_degenerate_2d ( a, b, c )

!*****************************************************************************80
!
!! LINE_IMP_IS_DEGENERATE_2D finds if an implicit point is degenerate in 2D.
!
!  Discussion:
!
!    The implicit form of a line in 2D is:
!
!      A * X + B * Y + C = 0
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, C, the implicit line parameters.
!
!    Output, logical LINE_IMP_IS_DEGENERATE_2D, is true if the
!    line is degenerate.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  logical line_imp_is_degenerate_2d

  call r8_fake_use ( c )

  line_imp_is_degenerate_2d = ( a * a + b * b == 0.0D+00 )

  return
end
subroutine line_imp_point_dist_2d ( a, b, c, p, dist )

!*****************************************************************************80
!
!! LINE_IMP_POINT_DIST_2D: distance ( implicit line, point ) in 2D.
!
!  Discussion:
!
!    The implicit form of a line in 2D is:
!
!      A * X + B * Y + C = 0
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, C, the implicit line parameters.
!
!    Input, real ( kind = rk ) P(2), the point whose distance from the line is
!    to be measured.
!
!    Output, real ( kind = rk ) DIST, the distance from the point to the line.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) dist
  logical line_imp_is_degenerate_2d
  real ( kind = rk ) p(dim_num)

  if ( line_imp_is_degenerate_2d ( a, b, c ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LINE_IMP_POINT_DIST_2D - Fatal error!'
    write ( *, '(a)' ) '  The line is degenerate.'
    stop 1
  end if

  dist = abs ( a * p(1) + b * p(2) + c ) / sqrt ( a * a + b * b )

  return
end
subroutine line_imp_point_dist_signed_2d ( a, b, c, p, dist_signed )

!*****************************************************************************80
!
!! LINE_IMP_POINT_DIST_SIGNED_2D: signed distance ( imp line, point ) in 2D.
!
!  Discussion:
!
!    The implicit form of a line in 2D is:
!
!      A * X + B * Y + C * Z + D = 0
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, C, the implicit line parameters.
!
!    Input, real ( kind = rk ) P(2), the coordinates of the point.
!
!    Output, real ( kind = rk ) DIST_SIGNED, the signed distance from the
!    point to the line.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) dist_signed
  logical line_imp_is_degenerate_2d
  real ( kind = rk ) p(dim_num)

  if ( line_imp_is_degenerate_2d ( a, b, c ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LINE_IMP_POINT_DIST_SIGNED_2D - Fatal error!'
    write ( *, '(a)' ) '  The line is degenerate.'
    stop 1
  end if

  dist_signed = - sign ( 1.0D+00, c ) * ( a * p(1) + b * p(2) + c ) / &
    sqrt ( a * a + b * b )

  return
end
subroutine line_imp2exp_2d ( a, b, c, p1, p2 )

!*****************************************************************************80
!
!! LINE_IMP2EXP_2D converts an implicit line to explicit form in 2D.
!
!  Discussion:
!
!    The implicit form of line in 2D is:
!
!      A * X + B * Y + C = 0
!
!    The explicit form of a line in 2D is:
!
!      the line through the points P1 and P2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 January 2005
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
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, C, the implicit line parameters.
!
!    Output, real ( kind = rk ) P1(2), P2(2), two points on the line.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  logical line_imp_is_degenerate_2d
  real ( kind = rk ) normsq
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)

  if ( line_imp_is_degenerate_2d ( a, b, c ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LINE_IMP2EXP_2D - Fatal error!'
    write ( *, '(a)' ) '  The line is degenerate.'
    stop 1
  end if

  normsq = a * a + b * b

  p1(1) = - a * c / normsq
  p1(2) = - b * c / normsq

  if ( abs ( b ) < abs ( a ) ) then
    p2(1) = - ( a - b / a ) * c / normsq
    p2(2) = - ( b + 1.0D+00 ) * c / normsq
  else
    p2(1) = - ( a + 1.0D+00 ) * c / normsq
    p2(2) = - ( b - a / b ) * c / normsq
  end if

  return
end
subroutine line_imp2par_2d ( a, b, c, f, g, x0, y0 )

!*****************************************************************************80
!
!! LINE_IMP2PAR_2D converts an implicit line to parametric form in 2D.
!
!  Discussion:
!
!    The implicit form of line in 2D is:
!
!      A * X + B * Y + C = 0
!
!    The parametric form of a line in 2D is:
!
!      X = X0 + F * T
!      Y = Y0 + G * T
!
!    We may normalize by choosing F*F + G*G = 1, and F nonnegative.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    24 January 1999
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
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, C, the implicit line parameters.
!
!    Output, real ( kind = rk ) F, G, X0, Y0, the parametric parameters of
!    the line.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) f
  logical line_imp_is_degenerate_2d
  real ( kind = rk ) g
  real ( kind = rk ) x0
  real ( kind = rk ) y0

  if ( line_imp_is_degenerate_2d ( a, b, c ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LINE_IMP2PAR_2D - Fatal error!'
    write ( *, '(a)' ) '  The line is degenerate.'
    stop 1
  end if

  x0 = - a * c / ( a * a + b * b )
  y0 = - b * c / ( a * a + b * b )

  f =   b / sqrt ( a * a + b * b )
  g = - a / sqrt ( a * a + b * b )

  if ( f < 0.0D+00 ) then
    f = -f
    g = -g
  end if

  return
end
subroutine line_par_point_dist_2d ( f, g, x0, y0, p, dist )

!*****************************************************************************80
!
!! LINE_PAR_POINT_DIST_2D: distance ( parametric line, point ) in 2D.
!
!  Discussion:
!
!    The parametric form of a line in 2D is:
!
!      X = X0 + F * T
!      Y = Y0 + G * T
!
!    We normalize by always choosing F*F + G*G = 1, and F nonnegative.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 January 2005
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
!  Parameters:
!
!    Input, real ( kind = rk ) F, G, X0, Y0, the parametric line parameters.
!
!    Input, real ( kind = rk ) P(2), the point whose distance from the line is
!    to be measured.
!
!    Output, real ( kind = rk ) DIST, the distance from the point to the line.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) dist
  real ( kind = rk ) dx
  real ( kind = rk ) dy
  real ( kind = rk ) f
  real ( kind = rk ) g
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) x0
  real ( kind = rk ) y0

  dx =   g * g * ( p(1) - x0 ) - f * g * ( p(2) - y0 )
  dy = - f * g * ( p(1) - x0 ) + f * f * ( p(2) - y0 )

  dist = sqrt ( dx * dx + dy * dy ) / ( f * f + g * g )

  return
end
subroutine line_par_point_dist_3d ( f, g, h, x0, y0, z0, p, dist )

!*****************************************************************************80
!
!! LINE_PAR_POINT_DIST_3D: distance ( parametric line, point ) in 3D.
!
!  Discussion:
!
!    The parametric form of a line in 3D is:
!
!      X = X0 + F * T
!      Y = Y0 + G * T
!      Z = Z0 + H * T
!
!    We normalize by always choosing F*F + G*G + H*H = 1, and F nonnegative.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 January 2005
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
!  Parameters:
!
!    Input, real ( kind = rk ) F, G, H, X0, Y0, Z0, the parametric line
!    parameters.
!
!    Input, real ( kind = rk ) P(3), the point whose distance from the line is
!    to be measured.
!
!    Output, real ( kind = rk ) DIST, the distance from the point to the line.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) dist
  real ( kind = rk ) dx
  real ( kind = rk ) dy
  real ( kind = rk ) dz
  real ( kind = rk ) f
  real ( kind = rk ) g
  real ( kind = rk ) h
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) x0
  real ( kind = rk ) y0
  real ( kind = rk ) z0

  dx =   g * ( f * ( p(2) - y0 ) - g * ( p(1) - x0 ) ) &
       + h * ( f * ( p(3) - z0 ) - h * ( p(1) - x0 ) )

  dy =   h * ( g * ( p(3) - z0 ) - h * ( p(2) - y0 ) )  &
       - f * ( f * ( p(2) - y0 ) - g * ( p(1) - x0 ) )

  dz = - f * ( f * ( p(3) - z0 ) - h * ( p(1) - x0 ) ) &
       - g * ( g * ( p(3) - z0 ) - h * ( p(2) - y0 ) )

  dist = sqrt ( dx * dx + dy * dy + dz * dz ) &
    / ( f * f + g * g + h * h )

  return
end
subroutine line_par_point_near_2d ( f, g, x0, y0, p, pn )

!*****************************************************************************80
!
!! LINE_PAR_POINT_NEAR_2D: nearest point on parametric line to given point, 2D.
!
!  Discussion:
!
!    The parametric form of a line in 2D is:
!
!      X = X0 + F * T
!      Y = Y0 + G * T
!
!    We may normalize by choosing F*F + G*G = 1, and F nonnegative.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 April 2013
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
!  Parameters:
!
!    Input, real ( kind = rk ) F, G, X0, Y0, the parametric line parameters.
!
!    Input, real ( kind = rk ) P(2), the point whose distance from the line is
!    to be measured.
!
!    Output, real ( kind = rk ) PN(2), the point on the parametric line which
!    is nearest to P.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) f
  real ( kind = rk ) g
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) pn(dim_num)
  real ( kind = rk ) t
  real ( kind = rk ) x0
  real ( kind = rk ) y0

  t = ( f * ( p(1) - x0 ) + g * ( p(2) - y0 ) ) / ( f * f + g * g )

  pn(1) = x0 + t * f
  pn(2) = y0 + t * g

  return
end
subroutine line_par_point_near_3d ( f, g, h, x0, y0, z0, p, pn )

!*****************************************************************************80
!
!! LINE_PAR_POINT_NEAR_3D: nearest point on parametric line to given point, 3D.
!
!  Discussion:
!
!    The parametric form of a line in 3D is:
!
!      X = X0 + F * T
!      Y = Y0 + G * T
!      Z = Z0 + H * T
!
!    We may normalize by choosing F*F + G*G + H*H = 1, and F nonnegative.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 April 2013
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
!  Parameters:
!
!    Input, real ( kind = rk ) F, G, H, X0, Y0, Z0, the parametric 
!    line parameters.
!
!    Input, real ( kind = rk ) P(3), the point whose distance from the line is
!    to be measured.
!
!    Output, real ( kind = rk ) PN(3), the point on the parametric line which
!    is nearest to P.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) f
  real ( kind = rk ) g
  real ( kind = rk ) h
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) pn(dim_num)
  real ( kind = rk ) t
  real ( kind = rk ) x0
  real ( kind = rk ) y0
  real ( kind = rk ) z0

  t = ( f * ( p(1) - x0 ) + g * ( p(2) - y0 ) + h * ( p(3) - z0 ) ) &
    / ( f * f + g * g + h * h )

  pn(1) = x0 + t * f
  pn(2) = y0 + t * g
  pn(3) = z0 + t * h

  return
end
subroutine line_par_print_2d ( f, g, x0, y0, label )

!*****************************************************************************80
!
!! line_par_print_2d: prints a parametric line in 2D.
!
!  Discussion:
!
!    The parametric form of a line in 2D is:
!
!      X = X0 + F * T
!      Y = Y0 + G * T
!
!    We may normalize by choosing F*F + G*G = 1, and F nonnegative.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 July 2020
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
!    real ( kind = rk ) F, G, X0, Y0, the parametric line parameters.
!
!    character ( len = * ) label: a label for the line.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) f
  real ( kind = rk ) g
  character ( len = * ) label
  real ( kind = rk ) x0
  real ( kind = rk ) y0

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) label
  write ( *, '(a,2f8.4)' ) '    F,  G =  ', f, g
  write ( *, '(a,2f8.4)' ) '    X0, Y0 = ', x0, y0

  return
end
subroutine line_par2exp_2d ( f, g, x0, y0, p1, p2 )

!*****************************************************************************80
!
!! LINE_PAR2EXP_2D converts a parametric line to explicit form in 2D.
!
!  Discussion:
!
!    The parametric form of a line in 2D is:
!
!      X = X0 + F * T
!      Y = Y0 + G * T
!
!    We normalize by always choosing F*F + G*G = 1, and F nonnegative.
!
!    The explicit form of a line in 2D is:
!
!      the line through the points P1 and P2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 January 2005
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
!  Parameters:
!
!    Input, real ( kind = rk ) F, G, X0, Y0, the parametric line parameters.
!
!    Output, real ( kind = rk ) P1(2), P2(2), two points on the line.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) f
  real ( kind = rk ) g
  real ( kind = rk ) x0
  real ( kind = rk ) y0
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)

  p1(1) = x0
  p1(2) = y0

  p2(1) = p1(1) + f
  p2(2) = p1(2) + g

  return
end
subroutine line_par2exp_3d ( f, g, h, x0, y0, z0, p1, p2 )

!*****************************************************************************80
!
!! LINE_PAR2EXP_3D converts a parametric line to explicit form in 3D.
!
!  Discussion:
!
!    The parametric form of a line in 3D is:
!
!      X = X0 + F * T
!      Y = Y0 + G * T
!      Z = Z0 + H * T
!
!    We may normalize by choosing F*F + G*G + H*H = 1, and F nonnegative.
!
!    The explicit form of a line in 3D is:
!
!      the line through the points P1 and P2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    13 April 2013
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
!  Parameters:
!
!    Input, real ( kind = rk ) F, G, H, X0, Y0, Z0, the parametric 
!    line parameters.
!
!    Output, real ( kind = rk ) P1(3), P2(3), two points on the line.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) f
  real ( kind = rk ) g
  real ( kind = rk ) h
  real ( kind = rk ) x0
  real ( kind = rk ) y0
  real ( kind = rk ) z0
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)

  p1(1) = x0
  p1(2) = y0
  p1(3) = z0

  p2(1) = p1(1) + f
  p2(2) = p1(2) + g
  p2(3) = p1(3) + h

  return
end
subroutine line_par2imp_2d ( f, g, x0, y0, a, b, c )

!*****************************************************************************80
!
!! LINE_PAR2IMP_2D converts a parametric line to implicit form in 2D.
!
!  Discussion:
!
!    The parametric form of a line in 2D is:
!
!      X = X0 + F * T
!      Y = Y0 + G * T
!
!    We normalize by always choosing F*F + G*G = 1, and F nonnegative.
!
!    The implicit form of a line in 2D is:
!
!      A * X + B * Y + C = 0
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    24 January 1999
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
!  Parameters:
!
!    Input, real ( kind = rk ) F, G, X0, Y0, the parametric line parameters.
!
!    Output, real ( kind = rk ) A, B, C, the implicit line parameters.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) f
  real ( kind = rk ) g
  real ( kind = rk ) x0
  real ( kind = rk ) y0

  a = -g
  b = f
  c = g * x0 - f * y0

  return
end
subroutine lines_exp_angle_3d ( p1, p2, q1, q2, angle )

!*****************************************************************************80
!
!! LINES_EXP_ANGLE_3D finds the angle between two explicit lines in 3D.
!
!  Discussion:
!
!    The explicit form of a line in 3D is:
!
!      the line through the points P1 and P2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(3), P2(3), two points on the first line.
!
!    Input, real ( kind = rk ) Q1(3), Q2(3), two points on the second line.
!
!    Output, real ( kind = rk ) ANGLE, the angle in radians between the two
!    lines.  The angle is computed using the ACOS function, and so lies between
!    0 and PI.  But if one of the lines is degenerate, the angle is 
!    returned as -1.0.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) angle
  real ( kind = rk ) ctheta
  logical line_exp_is_degenerate_nd
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) pdotq
  real ( kind = rk ) pnorm
  real ( kind = rk ) q1(dim_num)
  real ( kind = rk ) q2(dim_num)
  real ( kind = rk ) qnorm
  real ( kind = rk ) r8_acos

  if ( line_exp_is_degenerate_nd ( dim_num, p1, p2 ) ) then
!   write ( *, '(a)' ) ' '
!   write ( *, '(a)' ) 'LINES_EXP_ANGLE_3D - Fatal error!'
!   write ( *, '(a)' ) '  The line (P1,P2) is degenerate!'
    angle = -1.0D+00
    return
  end if

  if ( line_exp_is_degenerate_nd ( dim_num, q1, q2 ) ) then
!   write ( *, '(a)' ) ' '
!   write ( *, '(a)' ) 'LINES_EXP_ANGLE_3D - Warning!'
!   write ( *, '(a)' ) '  The line (Q1,Q2) is degenerate!'
    angle = -1.0D+00
    return
  end if

  pnorm = sqrt ( sum ( ( p2(1:dim_num) - p1(1:dim_num) )**2 ) )

  qnorm = sqrt ( sum ( ( q2(1:dim_num) - q1(1:dim_num) )**2 ) )

  pdotq = sum ( ( p2(1:dim_num) - p1(1:dim_num) ) &
              * ( q2(1:dim_num) - q1(1:dim_num) ) )

  ctheta = pdotq / ( pnorm * qnorm )

  angle = r8_acos ( ctheta )

  return
end
subroutine lines_exp_angle_nd ( dim_num, p1, p2, q1, q2, angle )

!*****************************************************************************80
!
!! LINES_EXP_ANGLE_ND returns the angle between two explicit lines in ND.
!
!  Discussion:
!
!    The explicit form of a line in ND is:
!
!      the line through the points P1 and P2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = rk ) P1(DIM_NUM), P2(DIM_NUM), two points 
!    on the first line.
!
!    Input, real ( kind = rk ) Q1(DIM_NUM), Q2(DIM_NUM), two points 
!    on the second line.
!
!    Output, real ( kind = rk ) ANGLE, the angle in radians between the two
!    lines.  The angle is computed using the ACOS function, and so lies
!    between 0 and PI.  But if one of the lines is degenerate, the angle
!    is returned as -1.0.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num

  real ( kind = rk ) angle
  real ( kind = rk ) ctheta
  logical line_exp_is_degenerate_nd
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) pdotq
  real ( kind = rk ) pnorm
  real ( kind = rk ) q1(dim_num)
  real ( kind = rk ) q2(dim_num)
  real ( kind = rk ) qnorm
  real ( kind = rk ) r8_acos

  if ( line_exp_is_degenerate_nd ( dim_num, p1, p2 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LINES_EXP_ANGLE_3D - Fatal error!'
    write ( *, '(a)' ) '  The line (P1,P2) is degenerate!'
    angle = -1.0D+00
    stop 1
  end if

  if ( line_exp_is_degenerate_nd ( dim_num, q1, q2 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LINES_EXP_ANGLE_3D - Fatal error!'
    write ( *, '(a)' ) '  The line (Q1,Q2) is degenerate!'
    angle = -1.0D+00
    stop 1
  end if

  pnorm = sqrt ( sum ( ( p2(1:dim_num) - p1(1:dim_num) )**2 ) )
  qnorm = sqrt ( sum ( ( q2(1:dim_num) - q1(1:dim_num) )**2 ) )

  pdotq = sum ( ( p2(1:dim_num) - p1(1:dim_num) ) &
              * ( q2(1:dim_num) - q1(1:dim_num) ) )

  ctheta = pdotq / ( pnorm * qnorm )
  angle = r8_acos ( ctheta )

  return
end
subroutine lines_exp_dist_3d ( p1, p2, q1, q2, dist )

!*****************************************************************************80
!
!! LINES_EXP_DIST_3D computes the distance between two explicit lines in 3D.
!
!  Discussion:
!
!    The explicit form of a line in 3D is:
!
!      the line through the points P1 and P2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(3), P2(3), two points on the first line.
!
!    Input, real ( kind = rk ) Q1(3), Q2(3), two points on the second line.  
!
!    Output, real ( kind = rk ) DIST, the distance between the lines.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) a11
  real ( kind = rk ) a12
  real ( kind = rk ) a13
  real ( kind = rk ) a21
  real ( kind = rk ) a22
  real ( kind = rk ) a23
  real ( kind = rk ) a31
  real ( kind = rk ) a32
  real ( kind = rk ) a33
  real ( kind = rk ) bot
  real ( kind = rk ) cr1
  real ( kind = rk ) cr2
  real ( kind = rk ) cr3
  real ( kind = rk ) dist
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) q1(dim_num)
  real ( kind = rk ) q2(dim_num)
  real ( kind = rk ) top
!
!  The distance is found by computing the volume of a parallelipiped,
!  and dividing by the area of its base.
!
!  But if the lines are parallel, we compute the distance by
!  finding the distance between the first line and any point
!  on the second line.
!
  a11 = q1(1) - p1(1)
  a12 = q1(2) - p1(2)
  a13 = q1(3) - p1(3)

  a21 = p2(1) - p1(1)
  a22 = p2(2) - p1(2)
  a23 = p2(3) - p1(3)

  a31 = q2(1) - q1(1)
  a32 = q2(2) - q1(2)
  a33 = q2(3) - q1(3)
!
!  Compute the cross product.
!
  cr1 = a22 * a33 - a23 * a32
  cr2 = a23 * a31 - a21 * a33
  cr3 = a21 * a32 - a22 * a31

  bot = sqrt ( cr1 * cr1 + cr2 * cr2 + cr3 * cr3 )

  if ( bot == 0.0D+00 ) then

    call line_exp_point_dist_3d ( p1, p2, q1, dist )

  else

    top = abs (   a11 * ( a22 * a33 - a23 * a32 ) &
                - a12 * ( a21 * a33 - a23 * a31 ) &
                + a13 * ( a21 * a32 - a22 * a31 ) )

    dist = top / bot

  end if

  return
end
subroutine lines_exp_dist_3d_2 ( p1, p2, q1, q2, dist )

!*****************************************************************************80
!
!! LINES_EXP_DIST_3D_2 computes the distance between two explicit lines in 3D.
!
!  Discussion:
!
!    The explicit form of a line in 3D is:
!
!      the line through the points P1 and P2.
!
!    This routine uses a method that is essentially independent of dimension.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(3), P2(3), two points on the first line.
!
!    Input, real ( kind = rk ) Q1(3), Q2(3), two points on the second line.  
!
!    Output, real ( kind = rk ) DIST, the distance between the lines.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) d
  real ( kind = rk ) det
  real ( kind = rk ) dist
  real ( kind = rk ) e
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) pn(dim_num)
  real ( kind = rk ) q1(dim_num)
  real ( kind = rk ) q2(dim_num)
  real ( kind = rk ) qn(dim_num)
  real ( kind = rk ) sn
  real ( kind = rk ) tn
  real ( kind = rk ) u(dim_num)
  real ( kind = rk ) v(dim_num)
  real ( kind = rk ) w0(dim_num)
!
!  Let U = (P2-P1) and V = (Q2-Q1) be the direction vectors on
!  the two lines.
!
  u(1:dim_num) = p2(1:dim_num) - p1(1:dim_num)
  v(1:dim_num) = q2(1:dim_num) - q1(1:dim_num)
!
!  Let SN be the unknown coordinate of the nearest point PN on line 1,
!  so that PN = P(SN) = P1 + SN * (P2-P1).
!
!  Let TN be the unknown coordinate of the nearest point QN on line 2,
!  so that QN = Q(TN) = Q1 + TN * (Q2-Q1).
!
!  Let W0 = (P1-Q1).
!
  w0(1:dim_num) = p1(1:dim_num) - q1(1:dim_num)
!
!  The vector direction WC = P(SN) - Q(TC) is unique (among directions)
!  perpendicular to both U and V, so
!
!    U dot WC = 0
!    V dot WC = 0
!
!  or, equivalently:
!
!    U dot ( P1 + SN * (P2 - P1) - Q1 - TN * (Q2 - Q1) ) = 0
!    V dot ( P1 + SN * (P2 - P1) - Q1 - TN * (Q2 - Q1) ) = 0
!
!  or, equivalently:
!
!    (u dot u ) * sn - (u dot v ) tc = -u * w0
!    (v dot u ) * sn - (v dot v ) tc = -v * w0
!
!  or, equivalently:
!
!   ( a  -b ) * ( sn ) = ( -d )
!   ( b  -c )   ( tc )   ( -e )
!
  a = dot_product ( u, u )
  b = dot_product ( u, v )
  c = dot_product ( v, v )
  d = dot_product ( u, w0 )
  e = dot_product ( v, w0 )
!
!  Check the determinant.
!
  det = - a * c + b * b

  if ( det == 0.0D+00 ) then
    sn = 0.0D+00
    if ( abs ( b ) < abs ( c ) ) then
      tn = e / c
    else
      tn = d / b
    end if
  else
    sn = ( c * d - b * e ) / det
    tn = ( b * d - a * e ) / det
  end if

  pn(1:dim_num) = p1(1:dim_num) + sn * ( p2(1:dim_num) - p1(1:dim_num) )
  qn(1:dim_num) = q1(1:dim_num) + tn * ( q2(1:dim_num) - q1(1:dim_num) )

  dist = sqrt ( sum ( ( pn(1:dim_num) - qn(1:dim_num) )**2 ) )

  return
end
function lines_exp_equal_2d ( p1, p2, q1, q2 )

!*****************************************************************************80
!
!! LINES_EXP_EQUAL_2D determines if two explicit lines are equal in 2D.
!
!  Discussion:
!
!    The explicit form of a line in 2D is:
!
!      the line through the points P1 and P2.
!
!    It is essentially impossible to accurately determine whether two
!    explicit lines are equal in 2D.  However, for form's sake, and
!    because occasionally the correct result can be determined, we
!    provide this routine.  Since divisions are avoided, if the
!    input data is exactly representable, the result should be
!    correct.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(2), P2(2), two points on the first line.
!
!    Input, real ( kind = rk ) Q1(2), Q2(2), two points on the second line.
!
!    Output, logical LINES_EXP_EQUAL_2D, is TRUE if the two lines 
!    are determined to be identical.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  logical lines_exp_equal_2d
  real ( kind = rk ) p1(2)
  real ( kind = rk ) p2(2)
  real ( kind = rk ) q1(2)
  real ( kind = rk ) q2(2)
  real ( kind = rk ) test1
  real ( kind = rk ) test2
  real ( kind = rk ) test3
  real ( kind = rk ) test4
!
!  Slope (P1,P2) = Slope (P2,Q1).
!
  test1 = ( p2(2) - p1(2) ) * ( q1(1) - p2(1) ) &
        - ( p2(1) - p1(1) ) * ( q1(2) - p2(2) )

  if ( test1 /= 0.0D+00 ) then
    lines_exp_equal_2d = .false.
    return
  end if
!
!  Slope (Q1,Q2) = Slope (P2,Q1).
!
  test2 = ( q2(2) - q1(2) ) * ( q1(1) - p2(1) ) &
        - ( q2(1) - q1(1) ) * ( q1(2) - p2(2) ) 

  if ( test2 /= 0.0D+00 ) then
    lines_exp_equal_2d = .false.
    return
  end if
!
!  Slope (P1,P2) = Slope (P1,Q2).
!
  test3 = ( p2(2) - p1(2) ) * ( q2(1) - p1(1) ) &
        - ( p2(1) - p1(1) ) * ( q2(2) - p1(2) )

  if ( test3 /= 0.0D+00 ) then
    lines_exp_equal_2d = .false.
    return
  end if
!
!  Slope (Q1,Q2) = Slope (P1,Q2).
!
  test4 = ( q2(2) - q1(2) ) * ( q2(1) - p1(1) ) &
        - ( q2(1) - q1(1) ) * ( q2(2) - p1(2) ) 

  if ( test4 /= 0.0D+00 ) then
    lines_exp_equal_2d = .false.
    return
  end if

  lines_exp_equal_2d = .true.

  return
end
subroutine lines_exp_int_2d ( p1, p2, q1, q2, ival, p )

!*****************************************************************************80
!
!! LINES_EXP_INT_2D determines where two explicit lines intersect in 2D.
!
!  Discussion:
!
!    The explicit form of a line in 2D is:
!
!      the line through the points P1 and P2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(2), P2(2), two points on the first line.
!
!    Input, real ( kind = rk ) Q1(2), Q2(2), two points on the second line.
!
!    Output, integer IVAL, reports on the intersection:
!    0, no intersection, the lines may be parallel or degenerate.
!    1, one intersection point, returned in P.
!    2, infinitely many intersections, the lines are identical.
!
!    Output, real ( kind = rk ) P(2), if IVAl = 1, P is
!    the intersection point.  Otherwise, P = 0.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) a1
  real ( kind = rk ) a2
  real ( kind = rk ) b1
  real ( kind = rk ) b2
  real ( kind = rk ) c1
  real ( kind = rk ) c2
  integer ival
  logical point_1
  logical point_2
  real ( kind = rk ) p(2)
  real ( kind = rk ) p1(2)
  real ( kind = rk ) p2(2)
  real ( kind = rk ) q1(2)
  real ( kind = rk ) q2(2)

  ival = 0
  p(1:dim_num) = 0.0D+00
!
!  Check whether either line is a point.
!
  if ( all ( p1(1:dim_num) == p2(1:dim_num) ) ) then
    point_1 = .true.
  else
    point_1 = .false.
  end if

  if ( all ( q1(1:dim_num) == q2(1:dim_num) ) ) then
    point_2 = .true.
  else
    point_2 = .false.
  end if
!
!  Convert the lines to ABC format.
!
  if ( .not. point_1 ) then
    call line_exp2imp_2d ( p1, p2, a1, b1, c1 )
  end if

  if ( .not. point_2 ) then
    call line_exp2imp_2d ( q1, q2, a2, b2, c2 )
  end if
!
!  Search for intersection of the lines.
!
  if ( point_1 .and. point_2 ) then
    if ( all ( p1(1:dim_num) == q1(1:dim_num) ) ) then
      ival = 1
      p(1:dim_num) = p1(1:dim_num)
    end if
  else if ( point_1 ) then
    if ( a2 * p1(1) + b2 * p1(2) == c2 ) then
      ival = 1
      p(1:dim_num) = p1(1:dim_num)
    end if
  else if ( point_2 ) then
    if ( a1 * q1(1) + b1 * q1(2) == c1 ) then
      ival = 1
      p(1:dim_num) = q1(1:dim_num)
    end if
  else
    call lines_imp_int_2d ( a1, b1, c1, a2, b2, c2, ival, p )
  end if

  return
end
subroutine lines_exp_near_3d ( p1, p2, q1, q2, pn, qn )

!*****************************************************************************80
!
!! LINES_EXP_NEAR_3D computes the nearest points on two explicit lines in 3D.
!
!  Discussion:
!
!    The explicit form of a line in 3D is:
!
!      the line through the points P1 and P2.
!
!    This routine uses a method that is essentially independent of dimension.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(3), P2(3), two points on the first line.
!
!    Input, real ( kind = rk ) Q1(3), Q2(3), two points on the second line.  
!
!    Output, real ( kind = rk ) PN(3), QN(3), the points on the first and
!    second lines that are nearest.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) d
  real ( kind = rk ) det
  real ( kind = rk ) e
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) pn(dim_num)
  real ( kind = rk ) q1(dim_num)
  real ( kind = rk ) q2(dim_num)
  real ( kind = rk ) qn(dim_num)
  real ( kind = rk ) sn
  real ( kind = rk ) tn
  real ( kind = rk ) u(dim_num)
  real ( kind = rk ) v(dim_num)
  real ( kind = rk ) w0(dim_num)
!
!  Let U = (P2-P1) and V = (Q2-Q1) be the direction vectors on
!  the two lines.
!
  u(1:dim_num) = p2(1:dim_num) - p1(1:dim_num)
  v(1:dim_num) = q2(1:dim_num) - q1(1:dim_num)
!
!  Let SN be the unknown coordinate of the nearest point PN on line 1,
!  so that PN = P(SN) = P1 + SN * (P2-P1).
!
!  Let TN be the unknown coordinate of the nearest point QN on line 2,
!  so that QN = Q(TN) = Q1 + TN * (Q2-Q1).
!
!  Let W0 = (P1-Q1).
!
  w0(1:dim_num) = p1(1:dim_num) - q1(1:dim_num)
!
!  The vector direction WC = P(SN) - Q(TC) is unique (among directions)
!  perpendicular to both U and V, so
!
!    U dot WC = 0
!    V dot WC = 0
!
!  or, equivalently:
!
!    U dot ( P1 + SN * (P2 - P1) - Q1 - TN * (Q2 - Q1) ) = 0
!    V dot ( P1 + SN * (P2 - P1) - Q1 - TN * (Q2 - Q1) ) = 0
!
!  or, equivalently:
!
!    (u dot u ) * sn - (u dot v ) tc = -u * w0
!    (v dot u ) * sn - (v dot v ) tc = -v * w0
!
!  or, equivalently:
!
!   ( a  -b ) * ( sn ) = ( -d )
!   ( b  -c )   ( tc )   ( -e )
!
  a = dot_product ( u, u )
  b = dot_product ( u, v )
  c = dot_product ( v, v )
  d = dot_product ( u, w0 )
  e = dot_product ( v, w0 )
!
!  Check the determinant.
!
  det = - a * c + b * b

  if ( det == 0.0D+00 ) then
    sn = 0.0D+00
    if ( abs ( b ) < abs ( c ) ) then
      tn = e / c
    else
      tn = d / b
    end if
  else
    sn = ( c * d - b * e ) / det
    tn = ( b * d - a * e ) / det
  end if

  pn(1:dim_num) = p1(1:dim_num) + sn * ( p2(1:dim_num) - p1(1:dim_num) )
  qn(1:dim_num) = q1(1:dim_num) + tn * ( q2(1:dim_num) - q1(1:dim_num) )

  return
end
function lines_exp_parallel_2d ( p1, p2, q1, q2 )

!*****************************************************************************80
!
!! LINES_EXP_PARALLEL_2D determines if two lines are parallel in 2D.
!
!  Discussion:
!
!    The explicit form of a line in 2D is:
!
!      the line through the points P1 and P2.
!
!    The test is essentially a comparison of slopes, but should be
!    more accurate than an explicit slope comparison, and unfazed
!    by degenerate cases.
!
!    On the other hand, there is NO tolerance for error.  If the
!    slopes differ by a single digit in the last place, then the
!    lines are judged to be nonparallel.  A more robust test would
!    be to compute the angle between the lines, because then it makes
!    sense to say the lines are "almost" parallel: the angle is small.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(2), P2(2), two points on the first line.
!
!    Input, real ( kind = rk ) Q1(2), Q2(2), two points on the second line.
!
!    Output, logical LINES_EXP_PARALLEL_2D is TRUE if the 
!    lines are parallel.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  logical lines_exp_parallel_2d
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) q1(dim_num)
  real ( kind = rk ) q2(dim_num)

  lines_exp_parallel_2d = ( p2(1) - p1(1) ) * ( q2(2) - q1(2) ) == &
                          ( q2(1) - q1(1) ) * ( p2(2) - p1(2) )

  return
end
function lines_exp_parallel_3d ( p1, p2, q1, q2 )

!*****************************************************************************80
!
!! LINES_EXP_PARALLEL_3D determines if two lines are parallel in 3D.
!
!  Discussion:
!
!    The explicit form of a line in 3D is:
!
!      the line through the points P1 and P2.
!
!    The points P1, P2 define a direction (P2-P1).  Similarly, the
!    points (Q1,Q2) define a direction (Q2-Q1).  The quantity
!
!      (P2-P1) dot (Q2-Q1) = norm(P2-P1) * norm(Q2-Q1) * cos ( angle )
!
!    Therefore, the following value is between 0 and 1;
!
!      abs ( (P2-P1) dot (Q2-Q1) / ( norm(P2-P1) * norm(Q2-Q1) ) )
!
!    and the lines are parallel if
!
!      abs ( (P2-P1) dot (Q2-Q1) / ( norm(P2-P1) * norm(Q2-Q1) ) ) = 1
!
!    We can rephrase this as requiring:
!
!      ( (P2-P1)dot(Q2-Q1) )^2 = (P2-P1)dot(P2-P1) * (Q2-Q1)dot(Q2-Q1)
!
!    which avoids division and square roots.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(3), P2(3), two points on the first line.
!
!    Input, real ( kind = rk ) Q1(3), Q2(3), two points on the second line.
!
!    Output, logical LINES_EXP_PARALLEL_3D is TRUE if the lines 
!    are parallel.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  logical lines_exp_parallel_3d
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) pdotp
  real ( kind = rk ) pdotq
  real ( kind = rk ) q1(dim_num)
  real ( kind = rk ) q2(dim_num)
  real ( kind = rk ) qdotq

  pdotq = dot_product ( p2(1:dim_num) - p1(1:dim_num), &
                        q2(1:dim_num) - q1(1:dim_num) )

  pdotp = dot_product ( p2(1:dim_num) - p1(1:dim_num), &
                        p2(1:dim_num) - p1(1:dim_num) )

  qdotq = dot_product ( q2(1:dim_num) - q1(1:dim_num), &
                        q2(1:dim_num) - q1(1:dim_num) )

  lines_exp_parallel_3d = ( pdotq * pdotq == pdotp * qdotq )

  return
end
subroutine lines_imp_angle_2d ( a1, b1, c1, a2, b2, c2, theta )

!*****************************************************************************80
!
!! LINES_IMP_ANGLE_2D finds the angle between two implicit lines in 2D.
!
!  Discussion:
!
!    The implicit form of a line in 2D is:
!
!      A * X + B * Y + C = 0
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    24 January 1999
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
!  Parameters:
!
!    Input, real ( kind = rk ) A1, B1, C1, the implicit parameters of the 
!    first line.
!
!    Input, real ( kind = rk ) A2, B2, C2, the implicit parameters of the
!    second line.
!
!    Output, real ( kind = rk ) THETA, the angle between the two lines.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) a1
  real ( kind = rk ) a2
  real ( kind = rk ) b1
  real ( kind = rk ) b2
  real ( kind = rk ) c1
  real ( kind = rk ) c2
  real ( kind = rk ) pdotq
  real ( kind = rk ) pnorm
  real ( kind = rk ) qnorm
  real ( kind = rk ) r8_acos
  real ( kind = rk ) theta

  call r8_fake_use ( c1 )
  call r8_fake_use ( c2 )

  pdotq = a1 * a2 + b1 * b2
  pnorm = sqrt ( a1 * a1 + b1 * b1 )
  qnorm = sqrt ( a2 * a2 + b2 * b2 )

  theta = r8_acos ( pdotq / ( pnorm * qnorm ) )

  return
end
subroutine lines_imp_dist_2d ( a1, b1, c1, a2, b2, c2, dist )

!*****************************************************************************80
!
!! LINES_IMP_DIST_2D determines the distance between two implicit lines in 2D.
!
!  Discussion:
!
!    If the lines intersect, then their distance is zero.
!    If the two lines are parallel, then they have a nonzero distance.
!
!    The implicit form of a line in 2D is:
!
!      A * X + B * Y + C = 0
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A1, B1, C1, define the first line.
!    At least one of A1 and B1 must be nonzero.
!
!    Input, real ( kind = rk ) A2, B2, C2, define the second line.
!    At least one of A2 and B2 must be nonzero.
!
!    Output, real ( kind = rk ) DIST, the distance between the two lines.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) a1
  real ( kind = rk ) a2
  real ( kind = rk ) b1
  real ( kind = rk ) b2
  real ( kind = rk ) c1
  real ( kind = rk ) c2
  real ( kind = rk ) dist
  logical line_imp_is_degenerate_2d
!
!  Refuse to handle degenerate lines.
!
  if ( line_imp_is_degenerate_2d ( a1, b1, c1 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LINES_IMP_DIST_2D - Fatal error!'
    write ( *, '(a)' ) '  Line 1 is degenerate.'
    stop 1
  end if

  if ( line_imp_is_degenerate_2d ( a2, b2, c2 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LINES_IMP_DIST_2D - Fatal error!'
    write ( *, '(a)' ) '  Line 2 is degenerate.'
    stop 1
  end if
!
!  Determine if the lines intersect.
!
  if ( a1 * b2 /= a2 * b1 ) then
    dist = 0.0D+00
    return
  end if
!
!  Determine the distance between the parallel lines.
!
  dist = abs ( c2 / sqrt ( a2 * a2 + b2 * b2 ) &
             - c1 / sqrt ( a1 * a1 + b1 * b1 ) )

  return
end
subroutine lines_imp_int_2d ( a1, b1, c1, a2, b2, c2, ival, p )

!*****************************************************************************80
!
!! LINES_IMP_INT_2D determines where two implicit lines intersect in 2D.
!
!  Discussion:
!
!    The implicit form of a line in 2D is:
!
!      A * X + B * Y + C = 0
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    25 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A1, B1, C1, define the first line.
!    At least one of A1 and B1 must be nonzero.
!
!    Input, real ( kind = rk ) A2, B2, C2, define the second line.
!    At least one of A2 and B2 must be nonzero.
!
!    Output, integer IVAL, reports on the intersection.
!
!    -1, both A1 and B1 were zero.
!    -2, both A2 and B2 were zero.
!     0, no intersection, the lines are parallel.
!     1, one intersection point, returned in P.
!     2, infinitely many intersections, the lines are identical.
!
!    Output, real ( kind = rk ) P(2), if IVAL = 1, then P is
!    the intersection point.  Otherwise, P = 0.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) a(dim_num,dim_num+1)
  real ( kind = rk ) a1
  real ( kind = rk ) a2
  real ( kind = rk ) b1
  real ( kind = rk ) b2
  real ( kind = rk ) c1
  real ( kind = rk ) c2
  integer info
  integer ival
  logical line_imp_is_degenerate_2d
  real ( kind = rk ) p(dim_num)

  p(1:dim_num) = 0.0D+00
!
!  Refuse to handle degenerate lines.
!
  if ( line_imp_is_degenerate_2d ( a1, b1, c1 ) ) then
    ival = -1
    return
  end if

  if ( line_imp_is_degenerate_2d ( a2, b2, c2 ) ) then
    ival = -2
    return
  end if
!
!  Set up and solve a linear system.
!
  a(1,1) = a1
  a(1,2) = b1
  a(1,3) = -c1

  a(2,1) = a2
  a(2,2) = b2
  a(2,3) = -c2

  call r8mat_solve ( 2, 1, a, info )
!
!  If the inverse exists, then the lines intersect at the solution point.
!
  if ( info == 0 ) then

    ival = 1
    p(1:dim_num) = a(1:dim_num,3)
!
!  If the inverse does not exist, then the lines are parallel
!  or coincident.  Check for parallelism by seeing if the
!  C entries are in the same ratio as the A or B entries.
!
  else

    ival = 0

    if ( a1 == 0.0D+00 ) then
      if ( b2 * c1 == c2 * b1 ) then
        ival = 2
      end if
    else
      if ( a2 * c1 == c2 * a1 ) then
        ival = 2
      end if
    end if

  end if

  return
end
subroutine lines_par_angle_2d ( f1, g1, x01, y01, f2, g2, x02, y02, theta )

!*****************************************************************************80
!
!! LINES_PAR_ANGLE_2D finds the angle between two parametric lines in 2D.
!
!  Discussion:
!
!    The parametric form of a line in 2D is:
!
!      X = X0 + F * T
!      Y = Y0 + G * T
!
!    We normalize by always choosing F*F + G*G = 1, and F nonnegative.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    24 January 1999
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
!  Parameters:
!
!    Input, real ( kind = rk ) F1, G1, X01, Y01, the parametric parameters of the
!    first line.
!
!    Input, real ( kind = rk ) F2, G2, X02, Y02, the parametric parameters of the
!    second line.
!
!    Output, real ( kind = rk ) THETA, the angle between the two lines.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) f1
  real ( kind = rk ) f2
  real ( kind = rk ) g1
  real ( kind = rk ) g2
  real ( kind = rk ) pdotq
  real ( kind = rk ) pnorm
  real ( kind = rk ) qnorm
  real ( kind = rk ) r8_acos
  real ( kind = rk ) theta
  real ( kind = rk ) x01
  real ( kind = rk ) x02
  real ( kind = rk ) y01
  real ( kind = rk ) y02

  call r8_fake_use ( x01 )
  call r8_fake_use ( x02 )
  call r8_fake_use ( y01 )
  call r8_fake_use ( y02 )

  pdotq = f1 * f2 + g1 * g2
  pnorm = sqrt ( f1 * f1 + g1 * g1 )
  qnorm = sqrt ( f2 * f2 + g2 * g2 )

  theta = r8_acos ( pdotq / ( pnorm * qnorm ) )

  return
end
subroutine lines_par_angle_3d ( f1, g1, h1, x01, y01, z01, f2, g2, h2, &
  x02, y02, z02, theta )

!*****************************************************************************80
!
!! LINES_PAR_ANGLE_3D finds the angle between two parametric lines in 3D.
!
!  Discussion:
!
!    The parametric form of a line in 3D is:
!
!      X = X0 + F * T
!      Y = Y0 + G * T
!      Z = Z0 + H * T
!
!    We normalize by always choosing F*F + G*G + H*H = 1, and F nonnegative.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 February 1999
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
!  Parameters:
!
!    Input, real ( kind = rk ) F1, G1, H1, X01, Y01, Z01, the parametric
!    parameters of the first line.
!
!    Input, real ( kind = rk ) F2, G2, H2, X02, Y02, Z02, the parametric
!    parameters of the second line.
!
!    Output, real ( kind = rk ) THETA, the angle between the two lines.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) f1
  real ( kind = rk ) f2
  real ( kind = rk ) g1
  real ( kind = rk ) g2
  real ( kind = rk ) h1
  real ( kind = rk ) h2
  real ( kind = rk ) pdotq
  real ( kind = rk ) pnorm
  real ( kind = rk ) qnorm
  real ( kind = rk ) r8_acos
  real ( kind = rk ) theta
  real ( kind = rk ) x01
  real ( kind = rk ) x02
  real ( kind = rk ) y01
  real ( kind = rk ) y02
  real ( kind = rk ) z01
  real ( kind = rk ) z02

  call r8_fake_use ( x01 )
  call r8_fake_use ( x02 )
  call r8_fake_use ( y01 )
  call r8_fake_use ( y02 )
  call r8_fake_use ( z01 )
  call r8_fake_use ( z02 )

  pdotq = f1 * f2 + g1 * g2 + h1 * h2
  pnorm = sqrt ( f1 * f1 + g1 * g1 + h1 * h1 )
  qnorm = sqrt ( f2 * f2 + g2 * g2 + h2 * h2 )

  theta = r8_acos ( pdotq / ( pnorm * qnorm ) )

  return
end
subroutine lines_par_dist_3d ( f1, g1, h1, x01, y01, z01, f2, g2, h2, &
  x02, y02, z02, dist )

!*****************************************************************************80
!
!! LINES_PAR_DIST_3D finds the distance between two parametric lines in 3D.
!
!  Discussion:
!
!    The parametric form of a line in 3D is:
!
!      X = X0 + F * T
!      Y = Y0 + G * T
!      Z = Z0 + H * T
!
!    We normalize by always choosing F*F + G*G + H*H = 1, and F nonnegative.
!
!    This code does not work for parallel or near parallel lines.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 February 1999
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
!  Parameters:
!
!    Input, real ( kind = rk ) F1, G1, H1, X01, Y01, Z01, the parametric
!    parameters of the first line.
!
!    Input, real ( kind = rk ) F2, G2, H2, X02, Y02, Z02, the parametric
!    parameters of the second line.
!
!    Output, real ( kind = rk ) DIST, the distance between the two lines.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) dist
  real ( kind = rk ) f1
  real ( kind = rk ) f2
  real ( kind = rk ) g1
  real ( kind = rk ) g2
  real ( kind = rk ) h1
  real ( kind = rk ) h2
  real ( kind = rk ) x01
  real ( kind = rk ) x02
  real ( kind = rk ) y01
  real ( kind = rk ) y02
  real ( kind = rk ) z01
  real ( kind = rk ) z02

  dist = abs ( ( x02 - x01 ) * ( g1 * h2 - g2 * h1 ) &
             + ( y02 - y01 ) * ( h1 * f2 - h2 * f1 ) &
             + ( z02 - z01 ) * ( f1 * g2 - f2 * g1 ) )  / &
             ( ( f1 * g2 - f2 * g1 ) ** 2 &
             + ( g1 * h2 - g2 * h1 ) ** 2 &
             + ( h1 * f2 - h2 * f1 ) ** 2 )

  return
end
subroutine lines_par_int_2d ( f1, g1, x1, y1, f2, g2, x2, y2, t1, t2, pint )

!*****************************************************************************80
!
!! LINES_PAR_INT_2D determines where two parametric lines intersect in 2D.
!
!  Discussion:
!
!    The parametric form of a line in 2D is:
!
!      X = X0 + F * T
!      Y = Y0 + G * T
!
!    We normalize by always choosing F*F + G*G = 1, and F nonnegative.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 January 1999
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
!  Parameters:
!
!    Input, real ( kind = rk ) F1, G1, X1, Y1, define the first parametric line.
!
!    Input, real ( kind = rk ) F2, G2, X2, Y2, define the second parametric line.
!
!    Output, real ( kind = rk ) T1, T2, the T parameters on the first and second
!    lines of the intersection point.
!
!    Output, real ( kind = rk ) PINT(2), the intersection point.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) det
  real ( kind = rk ) f1
  real ( kind = rk ) f2
  real ( kind = rk ) g1
  real ( kind = rk ) g2
  real ( kind = rk ) pint(dim_num)
  real ( kind = rk ) t1
  real ( kind = rk ) t2
  real ( kind = rk ) x1
  real ( kind = rk ) x2
  real ( kind = rk ) y1
  real ( kind = rk ) y2

  det = f2 * g1 - f1 * g2

  if ( det == 0.0D+00 ) then
    t1 = 0.0D+00
    t2 = 0.0D+00
    pint(1:dim_num) = 0.0D+00
  else
    t1 = ( f2 * ( y2 - y1 ) - g2 * ( x2 - x1 ) ) / det
    t2 = ( f1 * ( y2 - y1 ) - g1 * ( x2 - x1 ) ) / det
    pint(1) = x1 + f1 * t1
    pint(2) = y1 + g1 * t1
  end if

  return
end
subroutine loc2glob_3d ( cospitch, cosroll, cosyaw, sinpitch, sinroll, sinyaw, &
  globas, locpts, glopts )

!*****************************************************************************80
!
!! LOC2GLOB_3D converts from a local to global coordinate system in 3D.
!
!  Discussion:
!
!    A global coordinate system is given.
!
!    A local coordinate system has been translated to the point with
!    global coordinates GLOBAS, and rotated through a yaw, a pitch, and
!    a roll.
!
!    A point has local coordinates LOCPTS, and it is desired to know
!    the point's global coordinates GLOPTS.
!
!    The transformation may be written as
!
!      GLOB = GLOBAS + N_YAW * N_PITCH * N_ROLL * LOC
!
!    where
!
!               (  cos(Yaw)   -sin(Yaw)        0      )
!    N_YAW    = (  sin(Yaw)    cos(Yaw)        0      )
!               (      0           0           1      )
!
!               (  cos(Pitch)      0       sin(Pitch) )
!    N_PITCH =  (      0           1           0      )
!               ( -sin(Pitch)      0       cos(Pitch) )
!
!               (      1           0           0      )
!    N_ROLL =   (      0       cos(Roll)  -sin(Roll)  )
!               (      0       sin(Roll)   cos(Roll)  )
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) COSPITCH, COSROLL, COSYAW, the cosines of the
!    pitch, roll and yaw angles.
!
!    Input, real ( kind = rk ) SINPITCH, SINROLL, SINYAW, the sines of the pitch,
!    roll and yaw angles.
!
!    Input, real ( kind = rk ) GLOBAS(3), the global coordinates of the base
!    vector.
!
!    Input, real ( kind = rk ) LOCPTS(3), the local coordinates of the point.
!
!    Output, real ( kind = rk ) GLOPTS(3), the global coordinates of the point.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) cospitch
  real ( kind = rk ) cosroll
  real ( kind = rk ) cosyaw
  real ( kind = rk ) globas(dim_num)
  real ( kind = rk ) glopts(dim_num)
  real ( kind = rk ) locpts(dim_num)
  real ( kind = rk ) sinpitch
  real ( kind = rk ) sinroll
  real ( kind = rk ) sinyaw

  glopts(1) = globas(1) + (  cosyaw * cospitch ) * locpts(1) &
    + (  cosyaw * sinpitch * sinroll - sinyaw * cosroll ) * locpts(2) &
    + (  cosyaw * sinpitch * cosroll + sinyaw * sinroll ) * locpts(3)

  glopts(2) = globas(2) + (  sinyaw * cospitch ) * locpts(1) &
    + (  sinyaw * sinpitch * sinroll + cosyaw * cosroll ) * locpts(2) &
    + (  sinyaw * sinpitch * cosroll - cosyaw * sinroll ) * locpts(3)

  glopts(3) = globas(3) + ( -sinpitch ) * locpts(1) &
    + (  cospitch * sinroll ) * locpts(2) &
    + (  cospitch * cosroll ) * locpts(3)

  return
end
subroutine l4vec_print ( n, a, title )

!*****************************************************************************80
!
!! L4VEC_PRINT prints an L4VEC.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, logical A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  logical a(n)
  integer i
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,a,1x,l1)' ) i, ':', a(i)
  end do

  return
end
subroutine minabs ( x1, y1, x2, y2, x3, y3, xmin, ymin )

!*****************************************************************************80
!
!! MINABS finds a local minimum of F(X) = A * abs ( X ) + B.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    28 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = rk ) X1, Y1, X2, Y2, X3, Y3, are three sets of
!    data of the form ( X, F(X) ).  The three X values must be distinct.
!    On output, the data has been sorted so that X1 < X2 < X3,
!    and the Y values have been rearranged accordingly.
!
!    Output, real ( kind = rk ) XMIN, YMIN.  XMIN is a point within the interval
!    spanned by X1, X2 and X3, at which F takes its local minimum
!    value YMIN.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) slope
  real ( kind = rk ) slope12
  real ( kind = rk ) slope13
  real ( kind = rk ) slope23
  real ( kind = rk ) x1
  real ( kind = rk ) x2
  real ( kind = rk ) x3
  real ( kind = rk ) xmin
  real ( kind = rk ) y1
  real ( kind = rk ) y2
  real ( kind = rk ) y3
  real ( kind = rk ) ymin
!
!  Refuse to deal with coincident data.
!
  if ( x1 == x2 .or. x2 == x3 .or. x3 == x1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MINABS - Fatal error!'
    write ( *, '(a)' ) '  X values are equal.'
    stop 1
  end if
!
!  Sort the data.
!
  if ( x2 < x1 ) then
    call r8_swap ( x1, x2 )
    call r8_swap ( y1, y2 )
  end if

  if ( x3 < x1 ) then
    call r8_swap ( x1, x3 )
    call r8_swap ( y1, y3 )
  end if

  if ( x3 < x2 ) then
    call r8_swap ( x2, x3 )
    call r8_swap ( y2, y3 )
  end if
!
!  Now determine the slopes.
!
  slope12 = ( y2 - y1 ) / ( x2 - x1 )
  slope23 = ( y3 - y2 ) / ( x3 - x2 )
  slope13 = ( y3 - y1 ) / ( x3 - x1 )
!
!  Case 1: Minimum must be at an endpoint.
!
  if ( slope13 <= slope12 .or. 0.0D+00 <= slope12 ) then

    if ( y1 < y3 ) then
      xmin = x1
      ymin = y1
    else
      xmin = x3
      ymin = y3
    end if
!
!  Case 2: The curve decreases, and decreases faster than the line
!  joining the endpoints.
!
!  Whichever of SLOPE12 and SLOPE23 is the greater in magnitude
!  represents the actual slope of the underlying function.
!  Find where two lines of that slope, passing through the
!  endpoint data, intersect.
!
  else

    slope = max ( abs ( slope12 ), slope23 )

    xmin = 0.5D+00 * ( x1 + x3 + ( y1 - y3 ) / slope )
    ymin = y1 - slope * ( xmin - x1 )

  end if

  return
end
subroutine minquad ( x1, y1, x2, y2, x3, y3, xmin, ymin )

!*****************************************************************************80
!
!! MINQUAD finds a local minimum of F(X) = A * X * X + B * X + C.
!
!  Discussion:
!
!    MINQUAD is primarily intended as a utility routine.
!    The square of the distance function between a point
!    and a line segment has the form of F(X).  Hence, we can seek
!    the line on the second segment which minimizes the square of
!    the distance to the other line segment.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = rk ) X1, Y1, X2, Y2, X3, Y3, three sets of data
!    of the form ( X, F(X) ).  The three X values must be distinct.
!    On output, the data has been sorted so that X1 < X2 < X3,
!    and the Y values have been rearranged accordingly.
!
!    Output, real ( kind = rk ) XMIN, YMIN.  XMIN is a point within the interval
!    spanned by X1, X2 and X3, at which F takes its local minimum value YMIN.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ierror
  real ( kind = rk ) x
  real ( kind = rk ) x1
  real ( kind = rk ) x2
  real ( kind = rk ) x3
  real ( kind = rk ) xleft
  real ( kind = rk ) xmin
  real ( kind = rk ) xrite
  real ( kind = rk ) y
  real ( kind = rk ) y1
  real ( kind = rk ) y2
  real ( kind = rk ) y3
  real ( kind = rk ) ymin
!
!  Refuse to deal with coincident data.
!
  if ( x1 == x2 .or. x2 == x3 .or. x3 == x1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MINQUAD - Fatal error!'
    write ( *, '(a)' ) '  X values are equal.'
    stop 1
  end if
!
!  Find the interval endpoints.
!
  xleft = min ( x1, x2, x3 )
  xrite = max ( x1, x2, x3 )
!
!  Find the minimizer and its function value, over the three input points.
!
  if ( y1 <= y2 .and. y1 <= y3 ) then
    xmin = x1
    ymin = y1
  else if ( y2 <= y1 .and. y2 <= y3 ) then
    xmin = x2
    ymin = y2
  else
    xmin = x3
    ymin = y3
  end if
!
!  Find the minimizer and its function value over the real line.
!
  call parabola_ex ( x1, y1, x2, y2, x3, y3, x, y, ierror )
!
!  If F is linear, then take the already computed min.
!
  if ( ierror == 2 ) then
!
!  If F has a maximum, then take the already computed min.
!
  else if ( ymin < y ) then
!
!  If the minimizer is to the left, take the already computed min.
!
  else if ( x < xleft ) then
!
!  If the minimizer is to the right, take the already computed min.
!
  else if ( xrite < x ) then

  else

    xmin = x
    ymin = y

  end if

  return
end
subroutine octahedron_shape_3d ( point_num, face_num, face_order_max, &
  point_coord, face_order, face_point )

!*****************************************************************************80
!
!! OCTAHEDRON_SHAPE_3D describes an octahedron in 3D.
!
!  Discussion:
!
!    The vertices lie on the unit sphere.
!
!    The dual of the octahedron is the cube.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    08 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer POINT_NUM, the number of points.
!
!    Input, integer FACE_NUM, the number of faces.
!
!    Input, integer FACE_ORDER_MAX, the maximum number of vertices 
!    per face.
!
!    Output, real ( kind = rk ) POINT_COORD(3,POINT_NUM), the points.
!
!    Output, integer FACE_ORDER(FACE_NUM), the number of 
!    vertices per face.
!
!    Output, integer FACE_POINT(FACE_ORDER_MAX,FACE_NUM); 
!    FACE_POINT(I,J) is the index of the I-th point in the J-th face.  The
!    points are listed in the counter clockwise direction defined
!    by the outward normal at the face.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer face_num
  integer face_order_max
  integer, parameter :: dim_num = 3
  integer point_num

  integer face_order(face_num)
  integer face_point(face_order_max,face_num)
  real ( kind = rk ) point_coord(dim_num,point_num)
!
!  Set point coordinates.
!
  point_coord(1:dim_num,1:point_num) = reshape ( (/ &
     0.0D+00,  0.0D+00, -1.0D+00, &
     0.0D+00, -1.0D+00,  0.0D+00, &
     1.0D+00,  0.0D+00,  0.0D+00, &
     0.0D+00,  1.0D+00,  0.0D+00, &
    -1.0D+00,  0.0D+00,  0.0D+00, &
     0.0D+00,  0.0D+00,  1.0D+00 /), (/ dim_num, point_num /) )
!
!  Set the face orders.
!
  face_order(1:face_num) = (/ &
    3, 3, 3, 3, 3, 3, 3, 3 /)
!
!  Set the faces.
!
  face_point(1:face_order_max,1:face_num) = reshape ( (/ &
     1, 3, 2, &
     1, 4, 3, &
     1, 5, 4, &
     1, 2, 5, &
     2, 3, 6, &
     3, 4, 6, &
     4, 5, 6, &
     5, 2, 6 /), (/ face_order_max, face_num /) )

  return
end
subroutine octahedron_size_3d ( point_num, edge_num, face_num, face_order_max )

!*****************************************************************************80
!
!! OCTAHEDRON_SIZE_3D returns size information for an octahedron in 3D.
!
!  Discussion:
!
!    This routine can be called before calling OCTAHEDRON_SHAPE_3D,
!    so that space can be allocated for the arrays.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer POINT_NUM, the number of points.
!
!    Output, integer EDGE_NUM, the number of edges.
!
!    Output, integer FACE_NUM, the number of faces.
!
!    Output, integer FACE_ORDER_MAX, the maximum number of 
!    vertices per face.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer edge_num
  integer face_num
  integer face_order_max
  integer point_num

  point_num = 6
  edge_num = 12
  face_num = 8
  face_order_max = 3

  return
end
subroutine parallelogram_area_2d ( p, area )

!*****************************************************************************80
!
!! PARALLELOGRAM_AREA_2D computes the area of a parallelogram in 2D.
!
!  Discussion:
!
!    A parallelogram is a polygon having four sides, with the property
!    that each pair of opposite sides is paralell.
!
!    Given the first three vertices of the parallelogram, 
!    P1, P2, and P3, the fourth vertex must satisfy
!
!      P4 = P1 + ( P3 - P2 )
!
!    This routine uses the fact that the norm of the cross product 
!    of two vectors is the area of the parallelogram they form:
!
!      Area = ( P3 - P2 ) x ( P1 - P2 ).
!
!        P4<-----P3
!        /       /
!       /       /   
!      P1----->P2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 May 2010
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) P(2,4), the parallelogram vertices,
!    given in counterclockwise order.  The fourth vertex is ignored.
!
!  Output:
!
!    real ( kind = rk ) AREA, the (signed) area.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) area
  real ( kind = rk ) p(2,4)
!
!  Compute the cross product vector, which only has a single
!  nonzero component.  
!
  area = ( p(1,2) - p(1,1) ) * ( p(2,3) - p(2,1) ) &
       - ( p(2,2) - p(2,1) ) * ( p(1,3) - p(1,1) )

  return
end
subroutine parallelogram_area_3d ( p, area )

!*****************************************************************************80
!
!! PARALLELOGRAM_AREA_3D computes the area of a parallelogram in 3D.
!
!  Discussion:
!
!    A parallelogram is a polygon having four sides, with the property
!    that each pair of opposite sides is paralell.
!
!    A parallelogram in 3D must have the property that it is "really"
!    a 2D object, that is, that the four vertices that define it lie
!    in some plane.  
!
!    Given the first three vertices of the parallelogram (in 2D or 3D), 
!    P1, P2, and P3, the fourth vertex must satisfy
!
!      P4 = P1 + ( P3 - P2 )
!
!    This routine uses the fact that the norm of the cross product 
!    of two vectors is the area of the parallelogram they form:
!
!      Area = ( P3 - P2 ) x ( P1 - P2 ).
!
!        P4<-----P3
!        /       /
!       /       /   
!      P1----->P2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 May 2010
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) P(3,4), the parallelogram vertices,
!    given in counterclockwise order.  The fourth vertex is ignored.
!
!  Output:
!
!    real ( kind = rk ) AREA, the area
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) area
  real ( kind = rk ) cross(3)
  real ( kind = rk ) p(3,4)
!
!  Compute the cross product vector.
!
  cross(1) = ( p(2,2) - p(2,1) ) * ( p(3,3) - p(3,1) ) &
           - ( p(3,2) - p(3,1) ) * ( p(2,3) - p(2,1) )

  cross(2) = ( p(3,2) - p(3,1) ) * ( p(1,3) - p(1,1) ) &
           - ( p(1,2) - p(1,1) ) * ( p(3,3) - p(3,1) )

  cross(3) = ( p(1,2) - p(1,1) ) * ( p(2,3) - p(2,1) ) &
           - ( p(2,2) - p(2,1) ) * ( p(1,3) - p(1,1) )

  area = sqrt ( sum ( cross(1:3) ** 2 ) )

  return
end
function parallelogram_contains_point_2d ( p1, p2, p3, p )

!*****************************************************************************80
!
!! PARALLELOGRAM_CONTAINS_POINT_2D: is point inside a parallelogram in 2D.
!
!  Discussion:
!
!       P2..............
!       /              .
!      /              .
!     /              .
!    P1----------->P3
!
!    The algorithm used here essentially computes the barycentric
!    coordinates of the point P, and accepts it if both coordinates
!    are between 0 and 1.  ( For a triangle, they must be positive,
!    and sum to no more than 1.)  The same trick works for a parallelepiped.
!
!    05 August 2005: Thanks to Gernot Grabmair for pointing out that a previous
!    version of this routine was incorrect.
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
!  Parameters:
!
!    Input, real ( kind = rk ) P1(2), P2(2), P3(2), three corners of the 
!    parallelogram, with P1 between P2 and P3.
!
!    Input, real ( kind = rk ) P(2), the point to be checked.
!
!    Output, logical PARALLELOGRAM_CONTAINS_POINT_2D, 
!    is TRUE if P is inside the parallelogram.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) a(dim_num,dim_num+1)
  integer info
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) p3(dim_num)
  logical parallelogram_contains_point_2d
!
!  Set up the linear system
!
!    ( X2-X1  X3-X1 ) XSI(1)  = X-X1
!    ( Y2-Y1  Y3-Y1 ) XSI(2)    Y-Y1
!
!  which is satisfied by the barycentric coordinates of P.
!
  a(1,1) = p2(1) - p1(1)
  a(1,2) = p3(1) - p1(1)
  a(1,3) = p(1)  - p1(1)

  a(2,1) = p2(2) - p1(2)
  a(2,2) = p3(2) - p1(2)
  a(2,3) = p(2)  - p1(2)
!
!  Solve the linear system.
!
  call r8mat_solve ( dim_num, 1, a, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PARALLELOGRAM_CONTAINS_CONTAIN_2D - Fatal error!'
    write ( *, '(a)' ) '  The linear system is singular.'
    write ( *, '(a)' ) '  The input data does not form a proper triangle.'
    stop 1
  end if

  if ( a(1,3) < 0.0D+00 .or. 1.0D+00 < a(1,3) ) then
    parallelogram_contains_point_2d = .false.
  else if ( a(2,3) < 0.0D+00 .or. 1.0D+00 < a(2,3) ) then
    parallelogram_contains_point_2d = .false.
  else
    parallelogram_contains_point_2d = .true.
  end if

  return
end
function parallelogram_contains_point_3d ( p1, p2, p3, p )

!*****************************************************************************80
!
!! PARALLELOGRAM_CONTAINS_POINT_3D: point "inside" parallelogram in 3D.
!
!  Discussion:
!
!    The parallelogram is a 2-dimensional object in a 3D space.
!    For a point to be "inside" the parallelogram, it should
!    lie in the plane defined by the sides of the parallelogram,
!    and, within that plane, lie inside the parallelogram.
!
!    The algorithm constructs an auxilliary point P4, such that
!    P4-P1 is normal to P2-P1 and P3-P1.  The barycentric coordinates
!    of the point P can be used to determine if the point lies in
!    the plane, and within the parallelogram.
!
!       P2..............
!       /              .
!      /              .
!     /              .
!    P1----------->P3
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(3), P2(3), P3(3), three corners of the 
!    parallelogram, with P1 between P2 and P3.
!
!    Input, real ( kind = rk ) P(3), the point to be checked.
!
!    Output, logical PARALLELOGRAM_CONTAINS_POINT_3D, 
!    is TRUE if P is inside the parallelogram, or on its boundary.
!    A slight amount of leeway is allowed for error, since a three
!    dimensional point may lie exactly in the plane of the parallelogram,
!    and yet be computationally slightly outside it.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) a(dim_num,dim_num+1)
  real ( kind = rk ) r8vec_norm
  integer info
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) p3(dim_num)
  real ( kind = rk ) p4(dim_num)
  logical parallelogram_contains_point_3d
  real ( kind = rk ), parameter :: tol = 0.0001D+00
!
!  Turn the triangle into a tetrahedron by computing the normal to
!  P2-P1 and P3-P1.
!
  call r8vec_cross_product_3d ( p2(1:dim_num) - p1(1:dim_num), &
                  p3(1:dim_num) - p1(1:dim_num), p4 )

  p4(1:dim_num) = p4(1:dim_num) / r8vec_norm ( dim_num, p4 )
!
!  Set up the linear system
!
!    ( X2-X1  X3-X1 X4-X1 ) XSI(1)  = X-X1
!    ( Y2-Y1  Y3-Y1 Y4-Y1 ) XSI(2)    Y-Y1
!    ( Z2-Z1  Z3-Z1 Z4-Z1 ) XSI(3)    Z-Z1
!
!  which is satisfied by the barycentric coordinates of P.
!
  a(1,1) = p2(1) - p1(1)
  a(1,2) = p3(1) - p1(1)
  a(1,3) = p4(1) - p1(1)
  a(1,4) = p(1)  - p1(1)

  a(2,1) = p2(2) - p1(2)
  a(2,2) = p3(2) - p1(2)
  a(2,3) = p4(2) - p1(2)
  a(2,4) = p(2)  - p1(2)

  a(3,1) = p2(3) - p1(3)
  a(3,2) = p3(3) - p1(3)
  a(3,3) = p4(3) - p1(3)
  a(3,4) = p(3)  - p1(3)
!
!  Solve the linear system.
!
  call r8mat_solve ( dim_num, 1, a, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PARALLELOGRAM_CONTAINS_CONTAIN_3D - Fatal error!'
    write ( *, '(a)' ) '  The linear system is singular.'
    write ( *, '(a)' ) '  The input data does not form a proper triangle.'
    stop 1
  end if

  if ( a(1,4) < 0.0D+00 .or. 1.0D+00 < a(1,4) ) then
    parallelogram_contains_point_3d = .false.
  else if ( a(2,4) < 0.0D+00 .or. 1.0D+00 < a(2,4) ) then
    parallelogram_contains_point_3d = .false.
  else if ( tol < abs ( a(3,4) ) ) then
    parallelogram_contains_point_3d = .false.
  else
    parallelogram_contains_point_3d = .true.
  end if

  return
end
subroutine parallelogram_point_dist_3d ( p1, p2, p3, p, dist )

!*****************************************************************************80
!
!! PARALLELOGRAM_POINT_DIST_3D: distance ( parallelogram, point ) in 3D.
!
!  Discussion:
!
!       P2..............
!       /              .
!      /              .
!     /              .
!    P1----------->P3
!
!    Note that we are asking for the distance, in 3D, to a parallelogram,
!    which is a 2D object.
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
!  Parameters:
!
!    Input, real ( kind = rk ) P1(3), P2(3), P3(3), three corners of the 
!    parallelogram, with P1 between P2 and P3.
!
!    Input, real ( kind = rk ) P(3), the point which is to be checked.
!
!    Output, real ( kind = rk ) DIST, the distance from the point to the
!    parallelogram.  DIST is zero if the point lies exactly on the
!    parallelogram.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) dis13
  real ( kind = rk ) dis21
  real ( kind = rk ) dis34
  real ( kind = rk ) dis42
  real ( kind = rk ) dist
  logical inside
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) p3(dim_num)
  real ( kind = rk ) p4(dim_num)
  logical parallelogram_contains_point_3d
  real ( kind = rk ) pn(dim_num)
  real ( kind = rk ) pp(dim_num)
  real ( kind = rk ) t
  real ( kind = rk ) temp
!
!  Compute PP, the unit normal to X2-X1 and X3-X1:
!
  pp(1) = ( p2(2) - p1(2) ) * ( p3(3) - p1(3) ) &
        - ( p2(3) - p1(3) ) * ( p3(2) - p1(2) )
  pp(2) = ( p2(3) - p1(3) ) * ( p3(1) - p1(1) ) &
        - ( p2(1) - p1(1) ) * ( p3(3) - p1(3) )
  pp(3) = ( p2(1) - p1(1) ) * ( p3(2) - p1(2) ) &
        - ( p2(2) - p1(2) ) * ( p3(1) - p1(1) )

  temp = sqrt ( sum ( pp(1:dim_num) ** 2 ) )

  if ( temp == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PARALLELOGRAM_POINT_DIST_3D - Fatal error!'
    write ( *, '(a)' ) '  The normal vector is zero.'
    stop 1
  end if

  pp(1:dim_num) = pp(1:dim_num) / temp
!
!  Find PN, the nearest point to P in the plane.
!
  t = dot_product ( pp(1:dim_num), p(1:dim_num) - p1(1:dim_num) ) 

  pn(1:dim_num) = p(1:dim_num) - pp(1:dim_num) * t
!
!  If P lies WITHIN the parallelogram, we're done.
!
  inside = parallelogram_contains_point_3d ( p1, p2, p3, p )

  if ( inside ) then
    dist = sqrt ( sum ( ( pn(1:dim_num) - p(1:dim_num) ) ** 2 ) )
    return
  end if
!
!  Otherwise, find the distance between P and each of the
!  four line segments that make up the boundary of the parallelogram.
!
  p4(1:dim_num) = p2(1:dim_num) + p3(1:dim_num) - p1(1:dim_num)

  call segment_point_dist_3d ( p1, p3, p, dis13 )
  call segment_point_dist_3d ( p3, p4, p, dis34 )
  call segment_point_dist_3d ( p4, p2, p, dis42 )
  call segment_point_dist_3d ( p2, p1, p, dis21 )

  dist = min ( dis13, dis34, dis42, dis21 )

  return
end
subroutine parabola_ex ( x1, y1, x2, y2, x3, y3, x, y, ierror )

!*****************************************************************************80
!
!! PARABOLA_EX: extremal point of a parabola determined by three points.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) X1, Y1, X2, Y2, X3, Y3, the coordinates of 
!    three points on the parabola.  X1, X2 and X3 must be distinct.
!
!    Output, real ( kind = rk ) X, Y, the X coordinate of the extremal point
!    of the parabola, and the value of the parabola at that point.
!
!    Output, integer IERROR, error flag.
!    0, no error.
!    1, two of the X values are equal.
!    2, the data lies on a straight line; there is no finite extremal point.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) bot
  integer ierror
  real ( kind = rk ) x
  real ( kind = rk ) x1
  real ( kind = rk ) x2
  real ( kind = rk ) x3
  real ( kind = rk ) y
  real ( kind = rk ) y1
  real ( kind = rk ) y2
  real ( kind = rk ) y3

  ierror = 0

  if ( x1 == x2 .or. x2 == x3 .or. x3 == x1 ) then
    ierror = 1
    return
  end if

  if ( y1 == y2 .and. y2 == y3 .and. y3 == y1 ) then
    x = x1
    y = y1
    return
  end if

  bot = ( x2 - x3 ) * y1 - ( x1 - x3 ) * y2 + ( x1 - x2 ) * y3

  if ( bot == 0.0D+00 ) then
    ierror = 2
    return
  end if

  x = 0.5D+00 * ( x1 * x1 * ( y3 - y2 ) &
                + x2 * x2 * ( y1 - y3 ) &
                + x3 * x3 * ( y2 - y1 ) ) / bot

  y = (    ( x  - x2 ) * ( x  - x3 ) * ( x2 - x3 ) * y1 &
         - ( x  - x1 ) * ( x  - x3 ) * ( x1 - x3 ) * y2 &
         + ( x  - x1 ) * ( x  - x2 ) * ( x1 - x2 ) * y3 ) / &
         ( ( x1 - x2 ) * ( x2 - x3 ) * ( x1 - x3 ) )

  return
end
subroutine parabola_ex2 ( x1, y1, x2, y2, x3, y3, x, y, a, b, c, ierror )

!*****************************************************************************80
!
!! PARABOLA_EX2: extremal point of a parabola determined by three points.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    29 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) X1, Y1, X2, Y2, X3, Y3, the coordinates of 
!    three points on the parabola.  X1, X2 and X3 must be distinct.
!
!    Output, real ( kind = rk ) X, Y, the X coordinate of the extremal point
!    of the parabola, and the value of the parabola at that point.
!
!    Output, real ( kind = rk ) A, B, C, the coefficients that define the
!    parabola: P(X) = A * X * X + B * X + C.
!
!    Output, integer IERROR, error flag.
!    0, no error.
!    1, two of the X values are equal.
!    2, the data lies on a straight line; there is no finite extremal
!    point.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) det
  integer ierror
  real ( kind = rk ) v(3,3)
  real ( kind = rk ) w(3,3)
  real ( kind = rk ) x
  real ( kind = rk ) x1
  real ( kind = rk ) x2
  real ( kind = rk ) x3
  real ( kind = rk ) y
  real ( kind = rk ) y1
  real ( kind = rk ) y2
  real ( kind = rk ) y3

  ierror = 0

  if ( x1 == x2 .or. x2 == x3 .or. x3 == x1 ) then
    ierror = 1
    return
  end if

  if ( y1 == y2 .and. y2 == y3 .and. y3 == y1 ) then
    x = x1
    y = y1
    return
  end if
!
!  Set up the Vandermonde matrix.
!
  v(1,1) = 1.0D+00
  v(1,2) = x1
  v(1,3) = x1 * x1

  v(2,1) = 1.0D+00
  v(2,2) = x2
  v(2,3) = x2 * x2

  v(3,1) = 1.0D+00
  v(3,2) = x3
  v(3,3) = x3 * x3
!
!  Get the inverse.
!
  call r8mat_inverse_3d ( v, w, det )
!
!  Compute the parabolic coefficients.
!
  c = w(1,1) * y1 + w(1,2) * y2 + w(1,3) * y3
  b = w(2,1) * y1 + w(2,2) * y2 + w(2,3) * y3
  a = w(3,1) * y1 + w(3,2) * y2 + w(3,3) * y3
!
!  Determine the extremal point.
!
  if ( a == 0.0D+00 ) then
    ierror = 2
    return
  end if

  x = - b / ( 2.0D+00 * a )
  y = a * x * x + b * x + c

  return
end
function parallelepiped_contains_point_3d ( p1, p2, p3, p4, p )

!*****************************************************************************80
!
!! PARALLELEPIPED_CONTAINS_POINT_3D: point inside parallelepiped in 3D.
!
!  Discussion:
!
!    A parallelepiped is a "slanted box", that is, opposite
!    sides are parallel planes.
!
!         *------------------*
!        / .                / \
!       /   .              /   \
!      /     .            /     \
!    P4------------------*       \
!      \        .         \       \
!       \        .         \       \
!        \        .         \       \
!         \       P2.........\.......\
!          \     .            \     /
!           \   .              \   /
!            \ .                \ /
!            P1-----------------P3
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(3), P2(3), P3(3), P4(3), four corners 
!    of the parallelepiped.  It is assumed that P2, P3 and P4 are
!    immediate neighbors of P1.
!
!    Input, real ( kind = rk ) P(3), the point to be checked.
!
!    Output, logical PARALLELEPIPED_CONTAINS_POINT_3D, 
!    is true if P is inside the parallelepiped, or on its boundary.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) dot
  logical parallelepiped_contains_point_3d
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) p3(dim_num)
  real ( kind = rk ) p4(dim_num)
 
  parallelepiped_contains_point_3d = .false.

  dot = dot_product ( p(1:dim_num)  - p1(1:dim_num), &
                      p2(1:dim_num) - p1(1:dim_num) )

  if ( dot < 0.0D+00 ) then
    return
  end if

  if ( sum ( ( p2(1:dim_num) - p1(1:dim_num) ) ** 2 ) < dot ) then
    return
  end if

  dot = dot_product ( p(1:dim_num)  - p1(1:dim_num), &
                      p3(1:dim_num) - p1(1:dim_num) )

  if ( dot < 0.0D+00 ) then
    return
  end if

  if ( sum ( ( p3(1:dim_num) - p1(1:dim_num) ) ** 2 ) < dot ) then
    return
  end if

  dot = dot_product ( p(1:dim_num)  - p1(1:dim_num), &
                      p4(1:dim_num) - p1(1:dim_num) )

  if ( dot < 0.0D+00 ) then
    return
  end if

  if ( sum ( ( p4(1:dim_num) - p1(1:dim_num) ) ** 2 ) < dot ) then
    return
  end if

  parallelepiped_contains_point_3d = .true.

  return
end
subroutine parallelepiped_point_dist_3d ( p1, p2, p3, p4, p, dist )

!*****************************************************************************80
!
!! PARALLELEPIPED_POINT_DIST_3D: distance ( parallelepiped, point ) in 3D.
!
!  Discussion:
!
!    A parallelepiped is a "slanted box", that is, opposite
!    sides are parallel planes.
!
!         *------------------*
!        / .                / \
!       /   .              /   \
!      /     .            /     \
!    P4------------------*       \
!      \        .         \       \
!       \        .         \       \
!        \        .         \       \
!         \       P2.........\.......\
!          \     .            \     /
!           \   .              \   /
!            \ .                \ /
!            P1-----------------P3
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(3), P2(3), P3(3), P4(3), 
!    half of the corners of the box, from which the other corners can be
!    deduced.  The corners should be chosen so that the first corner
!    is directly connected to the other three.  The locations of
!    corners 5, 6, 7 and 8 will be computed by the parallelogram
!    relation.
!
!    Input, real ( kind = rk ) P(3), the point which is to be checked.
!
!    Output, real ( kind = rk ) DIST, the distance from the point to the box. 
!    DIST is zero if the point lies exactly on the box.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) dis
  real ( kind = rk ) dist
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) p3(dim_num)
  real ( kind = rk ) p4(dim_num)
  real ( kind = rk ) p5(dim_num)
  real ( kind = rk ) p6(dim_num)
  real ( kind = rk ) p7(dim_num)
  real ( kind = rk ) p8(dim_num)
!
!  Fill in the other corners
!
  p5(1:dim_num) = p2(1:dim_num) + p3(1:dim_num) - p1(1:dim_num)
  p6(1:dim_num) = p2(1:dim_num) + p4(1:dim_num) - p1(1:dim_num)
  p7(1:dim_num) = p3(1:dim_num) + p4(1:dim_num) - p1(1:dim_num)
  p8(1:dim_num) = p2(1:dim_num) + p3(1:dim_num) + p4(1:dim_num) &
    - 2.0D+00 * p1(1:dim_num)
!
!  Compute the distance from the point P to each of the six
!  parallelogram faces.
!
  call parallelogram_point_dist_3d ( p1, p2, p3, p, dis )

  dist = dis

  call parallelogram_point_dist_3d ( p1, p2, p4, p, dis )

  dist = min ( dist, dis )

  call parallelogram_point_dist_3d ( p1, p3, p4, p, dis )

  dist = min ( dist, dis )

  call parallelogram_point_dist_3d ( p8, p5, p6, p, dis )

  dist = min ( dist, dis )

  call parallelogram_point_dist_3d ( p8, p5, p7, p, dis )

  dist = min ( dist, dis )

  call parallelogram_point_dist_3d ( p8, p6, p7, p, dis )

  dist = min ( dist, dis )

  return
end
subroutine perm_inverse ( n, p )

!*****************************************************************************80
!
!! PERM_INVERSE inverts a permutation "in place".
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    25 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of objects being permuted.
!
!    Input/output, integer P(N), the permutation, in standard 
!    index form.  On output, P describes the inverse permutation
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer i0
  integer i1
  integer i2
  integer is
  integer p(n)

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_INVERSE - Fatal error!'
    write ( *, '(a,i8)' ) '  Input value of N = ', n
    stop 1
  end if

  is = 1

  do i = 1, n

    i1 = p(i)

    do while ( i < i1 )
      i2 = p(i1)
      p(i1) = -i2
      i1 = i2
    end do

    is = -sign ( 1, p(i) )
    p(i) = sign ( p(i), is )

  end do

  do i = 1, n

    i1 = -p(i)

    if ( 0 <= i1 ) then

      i0 = i

      do

        i2 = p(i1)
        p(i1) = i0

        if ( i2 < 0 ) then
          exit
        end if

        i0 = i1
        i1 = i2

      end do

    end if

  end do

  return
end
subroutine plane_exp_grid_3d ( p1, p2, p3, ncor3, line_num, cor3, lines, &
  maxcor3, line_max, ierror )

!*****************************************************************************80
!
!! PLANE_EXP_GRID_3D computes points and lines making up a planar grid in 3D.
!
!  Discussion:
!
!    The explicit form of a plane in 3D is:
!
!      the plane through P1, P2 and P3.
!
!    The data format used is that of SGI Inventor.
!
!    On input, if NCOR3 is zero (or negative), then the data computed by 
!    this routine will be stored normally in COR3.  But if NCOR3 is 
!    positive, it is assumed that COR3 already contains NCOR3 items
!    of useful data.  The new data is appended to COR3.  On output, NCOR3 
!    is increased by the number of points computed by this routine.
!
!    On input, if LINE_NUM is zero (or negative), then the data computed by 
!    this routine will be stored normally in LINES.  But if LINE_NUM is
!    positive, it is assumed that LINES already contains some useful data.  The
!    new data is appended to LINES.  On output, LINE_NUM is increased by the 
!    number of lines computed by this routine.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(3), P2(3), P3(3), three points on the plane.
!
!    Input/output, integer NCOR3, the number of points stored 
!    in COR3.
!
!    Input/output, integer LINE_NUM, the number of line data items.
!
!    Input/output, real ( kind = rk ) COR3(3,MAXCOR3), the grid points.
!
!    Input/output, integer LINES(LINE_MAX), the indices of 
!    points used in the lines of the grid.  Successive entries of LINES are 
!    joined by a line, unless an entry equals -1.  Note that indices begin
!    with 0.
!
!    Input, integer MAXCOR3, the maximum number of points.
!
!    Input, integer LINE_MAX, the maximum number of lines.
!
!    Output, integer IERROR, error indicator.
!    0, no error.
!    1, more space for point coordinates is needed.
!    2, more space for line data is needed.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer maxcor3
  integer line_max

  real ( kind = rk ) a
  real ( kind = rk ) amax
  real ( kind = rk ) amin
  real ( kind = rk ) b
  real ( kind = rk ) bmax
  real ( kind = rk ) bmin
  real ( kind = rk ) cor3(dim_num,maxcor3)
  real ( kind = rk ) dot
  integer i
  integer ierror
  integer j
  integer line_num
  integer lines(line_max)
  integer nbase
  integer ncor3
  integer, parameter :: nx = 5
  integer, parameter :: ny = 5
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) p3(dim_num)
  real ( kind = rk ) v1(dim_num)
  real ( kind = rk ) v2(dim_num)

  ierror = 0

  if ( ncor3 <= 0 ) then
    ncor3 = 0
  end if

  if ( line_num <= 0 ) then
    line_num = 0
  end if

  nbase = ncor3
!
!  Compute the two basis vectors for the affine plane.
!
  v1(1:dim_num) = p2(1:dim_num) - p1(1:dim_num)

  call vector_unit_nd ( dim_num, v1 )

  v2(1:dim_num) = p3(1:dim_num) - p1(1:dim_num)

  dot = dot_product ( v1(1:dim_num), v2(1:dim_num) )
!
!  Remove the component of V1 from V2, and give the
!  resulting vector unit norm.  V1 and V2 are now orthogonal
!  and of unit length, and represent the two direction vectors
!  of our plane.
!
  v2(1:dim_num) = v2(1:dim_num) - dot * v1(1:dim_num)

  call vector_unit_nd ( dim_num, v2 )
!
!  Compute the (V1,V2) coordinate range of the input data, if any.
!
  if ( ncor3 == 0 ) then

    amin = 0.0D+00
    amax = 1.0D+00
    bmin = 0.0D+00
    bmax = 1.0D+00

  else

    do i = 1, ncor3

      a = dot_product ( v1(1:dim_num), cor3(1:dim_num,i) )
      b = dot_product ( v2(1:dim_num), cor3(1:dim_num,i) )

      if ( i == 1 ) then
        amin = a
        amax = a
        bmin = b
        bmax = b
      else
        amin = min ( amin, a )
        amax = max ( amax, a )
        bmin = min ( bmin, b )
        bmax = max ( bmax, b )
      end if

    end do

  end if
!
!  Generate the points we will use.
!
  if ( maxcor3 < ncor3 + nx * ny ) then
    ierror = 1
    return
  end if

  do j = 1, ny

    b = ( real ( ny - j,     kind = rk ) * bmin &
        + real (      j - 1, kind = rk ) * bmax ) &
        / real ( ny     - 1, kind = rk )

    do i = 1, nx

      a = ( real ( nx - i,     kind = rk ) * amin &
          + real (      i - 1, kind = rk ) * amax ) &
          / real ( nx     - 1, kind = rk )

      ncor3 = ncor3 + 1
      cor3(1:dim_num,ncor3) = a * v1(1:dim_num) + b * v2(1:dim_num)

    end do

  end do
!
!  Do the "horizontals".
!
  do i = 1, nx

    do j = 1, ny

      if ( line_max <= line_num ) then
        ierror = 2
        return
      end if

      line_num = line_num + 1
      lines(line_num) = nbase + ( j - 1 ) * nx + i

    end do

    if ( line_max <= line_num ) then
      ierror = 2
      return
    end if

    line_num = line_num + 1
    lines(line_num) = 0

  end do
!
!  Do the "verticals".
!
  do j = 1, ny

    do i = 1, nx

      if ( line_max <= line_num ) then
        ierror = 2
        return
      end if

      line_num = line_num + 1
      lines(line_num) = nbase + ( j - 1 ) * nx + i

    end do

    if ( line_max <= line_num ) then
      ierror = 2
      return
    end if

    line_num = line_num + 1
    lines(line_num) = 0

  end do

  return
end
subroutine plane_exp_normal_3d ( p1, p2, p3, normal )

!*****************************************************************************80
!
!! PLANE_EXP_NORMAL_3D finds the normal to an explicit plane in 3D.
!
!  Discussion:
!
!    The explicit form of a plane in 3D is
!
!      the plane through P1, P2 and P3.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(3), P2(3), P3(3), three points on the plane.
!
!    Output, real ( kind = rk ) NORMAL(3), the coordinates of the unit normal
!    vector to the plane containing the three points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) normal(dim_num)
  real ( kind = rk ) normal_norm
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) p3(dim_num)
!
!  The cross product (P2-P1) x (P3-P1) is normal to (P2-P1) and (P3-P1).
!
  normal(1) = ( p2(2) - p1(2) ) * ( p3(3) - p1(3) ) &
            - ( p2(3) - p1(3) ) * ( p3(2) - p1(2) )

  normal(2) = ( p2(3) - p1(3) ) * ( p3(1) - p1(1) ) &
            - ( p2(1) - p1(1) ) * ( p3(3) - p1(3) )

  normal(3) = ( p2(1) - p1(1) ) * ( p3(2) - p1(2) ) &
            - ( p2(2) - p1(2) ) * ( p3(1) - p1(1) )

  normal_norm = sqrt ( sum ( normal(1:dim_num) ** 2 ) )

  if ( normal_norm == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLANE_EXP_NORMAL_3D - Fatal error!'
    write ( *, '(a)' ) '  The plane is poorly defined.'
    stop 1
  end if

  normal(1:dim_num) = normal(1:dim_num) / normal_norm

  return
end
subroutine plane_exp_point_dist_3d ( p1, p2, p3, p, dist )

!*****************************************************************************80
!
!! PLANE_EXP_POINT_DIST_3D: distance ( explicit plane, point ) in 3D.
!
!  Discussion:
!
!    The explicit form of a plane in 3D is
!
!      the plane through P1, P2 and P3.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(3), P2(3), P3(3), three points on the plane.
!
!    Input, real ( kind = rk ) P(3), the coordinates of the point.
!
!    Output, real ( kind = rk ) DIST, the distance from the point to the plane.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) d
  real ( kind = rk ) dist
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) p3(dim_num)

  call plane_exp2imp_3d ( p1, p2, p3, a, b, c, d )

  call plane_imp_point_dist_3d ( a, b, c, d, p, dist )

  return
end
subroutine plane_exp_pro2 ( p1, p2, p3, n, p, pp )

!*****************************************************************************80
!
!! PLANE_EXP_PRO2 produces 2D coordinates of points that lie in a plane, in 3D.
!
!  Discussion:
!
!    The explicit form of a plane in 3D is:
!
!      the plane through P1, P2 and P3.
!
!    The first thing to do is to compute two orthonormal vectors V1 and
!    V2, so that any point P that lies in the plane may be written as
!
!      P = P1 + alpha * V1 + beta * V2
!
!    The vector V1 lies in the direction P2-P1, and V2 lies in
!    the plane, is orthonormal to V1, and has a positive component
!    in the direction of P3-P1.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(3), P2(3), P3(3), three points on the plane.
!
!    Input, integer N, the number of points to project.
!
!    Input, real ( kind = rk ) P(3,N), are the Cartesian
!    coordinates of points which lie on the plane spanned by the
!    three points.  These points are not checked to ensure that
!    they lie on the plane.
!
!    Output, real ( kind = rk ) PP(2,N), the "in-plane"
!    coordinates of the points.  
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: dim_num = 3

  real ( kind = rk ) dot
  integer i
  real ( kind = rk ) p(dim_num,n)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) p3(dim_num)
  real ( kind = rk ) pp(2,dim_num)
  real ( kind = rk ) v1(dim_num)
  real ( kind = rk ) v2(dim_num)
!
!  Compute the two basis vectors for the affine plane.
!
  v1(1:dim_num) = p2(1:dim_num) - p1(1:dim_num)

  call vector_unit_nd ( dim_num, v1 )

  v2(1:dim_num) = p3(1:dim_num) - p1(1:dim_num)

  dot = dot_product ( v1(1:dim_num), v2(1:dim_num) )

  v2(1:dim_num) = v2(1:dim_num) - dot * v1(1:dim_num)

  call vector_unit_nd ( dim_num, v2 )
!
!  Now decompose each point.
!
  do i = 1, n
    pp(1,i) = dot_product ( p(1:dim_num,i) - p1(1:dim_num), v1(1:dim_num) )
    pp(2,i) = dot_product ( p(1:dim_num,i) - p2(1:dim_num), v2(1:dim_num) )
  end do

  return
end
subroutine plane_exp_pro3 ( p1, p2, p3, n, p, pp )

!*****************************************************************************80
!
!! PLANE_EXP_PRO3 projects points orthographically onto a plane, in 3D.
!
!  Discussion:
!
!    The explicit form of a plane in 3D is:
!
!      the plane through P1, P2 and P3.
!
!    PP may share the same memory as PO, in
!    which case the projections will overwrite the original data.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(3), P2(3), P3(3), three points on the plane.
!
!    Input, integer N, the number of points to project.
!
!    Input, real ( kind = rk ) P(3,N), the points.
!
!    Output, real ( kind = rk ) PP(3,N), the projections of the points through 
!    the focus point onto the plane.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: dim_num = 3

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) d
  integer i
  real ( kind = rk ) p(dim_num,n)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) p3(dim_num)
  real ( kind = rk ) pp(dim_num,n)
!
!  Put the plane into ABCD form.
!
  call plane_exp2imp_3d ( p1, p2, p3, a, b, c, d )
!
!  For each point, its image in the plane is the nearest point
!  in the plane.
!
  do i = 1, n

    call plane_imp_point_near_3d ( a, b, c, d, p(1:dim_num,i), pp(1:dim_num,i) )

  end do

  return
end
subroutine plane_exp_project_3d ( p1, p2, p3, pf, n, po, pp, ivis )

!*****************************************************************************80
!
!! PLANE_EXP_PROJECT_3D projects points through a point onto a plane in 3D.
!
!  Discussion:
!
!    The explicit form of a plane in 3D is
!
!      the plane through P1, P2 and P3.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(3), P2(3), P3(3), three points on the plane.
!
!    Input, real ( kind = rk ) PF(3), the focus point.
!
!    Input, integer N, the number of points to project.
!
!    Input, real ( kind = rk ) PO(3,N), the object points.
!
!    Output, real ( kind = rk ) PP(3,N), are the 
!    coordinates of the projections of the object points through the focus
!    point onto the plane.  PP may share the same memory as PO,
!    in which case the projections will overwrite the original data.
!
!    Output, integer IVIS(N), visibility indicator:
!    3, the object was behind the plane;
!    2, the object was already on the plane;
!    1, the object was between the focus and the plane;
!    0, the line from the object to the focus is parallel to the plane,
!    so the object is "invisible".
!    -1, the focus is between the object and the plane.  The object
!    might be considered invisible.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: dim_num = 3

  real ( kind = rk ) a
  real ( kind = rk ) alpha
  real ( kind = rk ) angle_rad_3d
  real ( kind = rk ) b
  real ( kind = rk ) beta
  real ( kind = rk ) c
  real ( kind = rk ) d
  real ( kind = rk ) disfo
  real ( kind = rk ) disfn
  integer i
  integer ivis(n)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) p3(dim_num)
  real ( kind = rk ) pf(dim_num)
  real ( kind = rk ) pn(dim_num)
  real ( kind = rk ) po(dim_num,n)
  real ( kind = rk ) pp(dim_num,n)
!
!  Put the plane into ABCD form.
!
  call plane_exp2imp_3d ( p1, p2, p3, a, b, c, d )
!
!  Get the nearest point on the plane to the focus.
!
  call plane_imp_point_near_3d ( a, b, c, d, pf, pn )
!
!  Get the distance from the focus to the plane.
!
  disfn = sqrt ( sum ( ( pf(1:dim_num) - pn(1:dim_num) ) ** 2 ) )
!
!  If the focus lies in the plane, this is bad.  We could still
!  project points that actually lie in the plane, but we'll
!  just bail out.
!
  if ( disfn == 0.0D+00 ) then
    ivis(1:n) = 0
    do i = 1, dim_num
      pp(i,1:n) = pf(i)
    end do
    return
  end if
!
!  Process the points.
!
  do i = 1, n
!
!  Get the distance from the focus to the object.
!
    disfo = sqrt ( sum ( ( po(1:dim_num,i) - pf(1:dim_num) ) ** 2 ) )

    if ( disfo == 0.0D+00 ) then

      ivis(i) = 0
      pp(1:dim_num,i) = pn(1:dim_num)

    else
!
!  Compute ALPHA, the angle between (OBJECT-FOCUS) and (NEAREST-FOCUS).
!
      alpha = angle_rad_3d ( po(1:3,i), pf(1:3), pn(1:3) )

      if ( cos ( alpha ) == 0.0D+00 ) then

        ivis(i) = 0
        pp(1:dim_num,i) = pn(1:dim_num)

      else
!
!  BETA is Dist(NEAREST-FOCUS) / ( Cos(ALPHA)*Dist(OBJECT-FOCUS) )
!
        beta = disfn / ( cos ( alpha ) * disfo )

        if ( 1.0D+00 < beta ) then
          ivis(i) = 1
        else if ( beta == 1.0D+00 ) then
          ivis(i) = 2
        else if ( 0.0D+00 < beta ) then
          ivis(i) = 3
        else
          ivis(i) = -1
        end if
!
!  Set the projected point.
!
        pp(1:dim_num,i) = pf(1:dim_num) &
          + beta * ( po(1:dim_num,i) - pf(1:dim_num) )

      end if

    end if

  end do

  return
end
subroutine plane_exp2imp_3d ( p1, p2, p3, a, b, c, d )

!*****************************************************************************80
!
!! PLANE_EXP2IMP_3D converts an explicit plane to implicit form in 3D.
!
!  Discussion:
!
!    The explicit form of a plane in 3D is
!
!      the plane through P1, P2 and P3.
!
!    The implicit form of a plane in 3D is
!
!      A * X + B * Y + C * Z + D = 0
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 February 2005
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
!  Parameters:
!
!    Input, real ( kind = rk ) P1(3), P2(3), P3(3), three points on the plane.
!
!    Output, real ( kind = rk ) A, B, C, D, coefficients which describe 
!    the plane.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) d
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) p3(dim_num)

  a = ( p2(2) - p1(2) ) * ( p3(3) - p1(3) ) &
    - ( p2(3) - p1(3) ) * ( p3(2) - p1(2) )

  b = ( p2(3) - p1(3) ) * ( p3(1) - p1(1) ) &
    - ( p2(1) - p1(1) ) * ( p3(3) - p1(3) )

  c = ( p2(1) - p1(1) ) * ( p3(2) - p1(2) ) &
    - ( p2(2) - p1(2) ) * ( p3(1) - p1(1) )

  d = - p2(1) * a - p2(2) * b - p2(3) * c

  return
end
subroutine plane_exp2normal_3d ( p1, p2, p3, pp, normal )

!*****************************************************************************80
!
!! PLANE_EXP2NORMAL_3D converts an explicit plane to normal form in 3D.
!
!  Discussion:
!
!    The explicit form of a plane in 3D is
!
!      the plane through P1, P2 and P3.
!
!    The normal form of a plane in 3D is
!
!      PP, a point on the plane, and
!      N, the unit normal to the plane.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(3), P2(3), P3(3), three points on the plane.
!
!    Output, real ( kind = rk ) PP(3), a point on the plane.
!
!    Output, real ( kind = rk ) NORMAL(3), a unit normal vector to the plane.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) norm
  real ( kind = rk ) normal(dim_num)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) p3(dim_num)
  real ( kind = rk ) pp(dim_num)

  pp(1:dim_num) = p1(1:dim_num)

  normal(1) = ( p2(2) - p1(2) ) * ( p3(3) - p1(3) ) &
            - ( p2(3) - p1(3) ) * ( p3(2) - p1(2) )

  normal(2) = ( p2(3) - p1(3) ) * ( p3(1) - p1(1) ) &
            - ( p2(1) - p1(1) ) * ( p3(3) - p1(3) )

  normal(3) = ( p2(1) - p1(1) ) * ( p3(2) - p1(2) ) &
            - ( p2(2) - p1(2) ) * ( p3(1) - p1(1) )

  norm = sqrt ( sum ( normal(1:dim_num) ** 2 ) )

  if ( norm == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLANE_EXP2NORMAL_3D - Fatal error!'
    write ( *, '(a)' ) '  The normal vector is null.'
    write ( *, '(a)' ) '  Two points coincide, or nearly so.'
    stop 1
  end if

  normal(1:dim_num) = normal(1:dim_num) / norm

  return
end
function plane_imp_is_degenerate_3d ( a, b, c )

!*****************************************************************************80
!
!! PLANE_IMP_IS_DEGENERATE_3D is TRUE if an implicit plane is degenerate.
!
!  Discussion:
!
!    The implicit form of a plane in 3D is:
!
!      A * X + B * Y + C * Z + D = 0
!
!    The implicit plane is degenerate if A = B = C = 0.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, C, the implicit plane parameters.
!
!    Output, logical PLANE_IMP_IS_DEGENERATE_3D, 
!    is TRUE if the plane is degenerate.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  logical plane_imp_is_degenerate_3d

  if ( a == 0.0D+00 .and. b == 0.0D+00 .and. c == 0.0D+00 ) then
    plane_imp_is_degenerate_3d = .true.
  else
    plane_imp_is_degenerate_3d = .false.
  end if

  return
end
subroutine plane_imp_line_par_int_3d ( a, b, c, d, x0, y0, z0, f, g, h, &
  intersect, p )

!*****************************************************************************80
!
!! PLANE_IMP_LINE_PAR_INT_3D: intersection ( impl plane, param line ) in 3D.
!
!  Discussion:
!
!    The implicit form of a plane in 3D is:
!
!      A * X + B * Y + C * Z + D = 0
!
!    The parametric form of a line in 3D is:
!
!      X = X0 + F * T
!      Y = Y0 + G * T
!      Z = Z0 + H * T
!
!    We normalize by always choosing F*F + G*G + H*H = 1, 
!    and F nonnegative.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 January 2005
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
!    ISBN: 0408012420,
!    page 111.
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, C, D, the implicit plane parameters.
!
!    Input, real ( kind = rk ) X0, Y0, Z0, F, G, H, parameters that define the
!    parametric line.
!
!    Output, logical INTERSECT, is TRUE if the line and the plane
!    intersect.
!
!    Output, real ( kind = rk ) P(3), is a point of intersection of the line
!    and the plane, if INTERSECT is TRUE.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) d
  real ( kind = rk ) denom
  real ( kind = rk ) f
  real ( kind = rk ) g
  real ( kind = rk ) h
  logical intersect
  real ( kind = rk ) norm1
  real ( kind = rk ) norm2
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) t
  real ( kind = rk ), parameter :: tol = 0.00001D+00
  real ( kind = rk ) x0
  real ( kind = rk ) y0
  real ( kind = rk ) z0
!
!  Check.
!
  norm1 = sqrt ( a * a + b * b + c * c )

  if ( norm1 == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLANE_IMP_LINE_PAR_INT_3D - Fatal error!'
    write ( *, '(a)' ) '  The plane normal vector is null.'
    stop 1
  end if

  norm2 = sqrt ( f * f + g * g + h * h )

  if ( norm2 == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLANE_IMP_LINE_PAR_INT_3D - Fatal error!'
    write ( *, '(a)' ) '  The line direction vector is null.'
    stop 1
  end if

  denom = a * f + b * g + c * h
!
!  The line and the plane may be parallel.
!
  if ( abs ( denom ) < tol * norm1 * norm2 ) then

    if ( a * x0 + b * y0 + c * z0 + d == 0.0D+00 ) then
      intersect = .true.
      p(1) = x0
      p(2) = y0
      p(3) = z0
    else
      intersect = .false.
      p(1:dim_num) = 0.0D+00
    end if
!
!  If they are not parallel, they must intersect.
!
  else

    intersect = .true.
    t = - ( a * x0 + b * y0 + c * z0 + d ) / denom
    p(1) = x0 + t * f
    p(2) = y0 + t * g
    p(3) = z0 + t * h

  end if

  return
end
subroutine plane_imp_point_dist_3d ( a, b, c, d, p, dist )

!*****************************************************************************80
!
!! PLANE_IMP_POINT_DIST_3D: distance ( implicit plane, point ) in 3D.
!
!  Discussion:
!
!    The implicit form of a plane in 3D is:
!
!      A * X + B * Y + C * Z + D = 0
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 January 2005
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
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, C, D, the implicit plane parameters.
!
!    Input, real ( kind = rk ) P(3), the coordinates of the point.
!
!    Output, real ( kind = rk ) DIST, the distance from the point to the plane.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) d
  real ( kind = rk ) dist
  real ( kind = rk ) norm
  real ( kind = rk ) p(dim_num)

  norm = sqrt ( a * a + b * b + c * c )

  if ( norm == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLANE_IMP_POINT_DIST_3D - Fatal error!'
    write ( *, '(a)' ) '  The plane normal vector is null.'
    stop 1
  end if

  dist = abs ( a * p(1) + b * p(2) + c * p(3) + d ) / norm

  return
end
subroutine plane_imp_point_dist_signed_3d ( a, b, c, d, p, dist_signed )

!*****************************************************************************80
!
!! PLANE_IMP_POINT_DIST_SIGNED_3D: signed distance ( imp plane, point) in 3D.
!
!  Discussion:
!
!    The implicit form of a plane in 3D is:
!
!      A * X + B * Y + C * Z + D = 0
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Priamos Georgiades,
!    Signed Distance From Point To Plane,
!    in Graphics Gems III,
!    edited by David Kirk,
!    Academic Press, 1992, pages 233-235, T385.G6973.
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, C, D, the implicit plane parameters.
!
!    Input, real ( kind = rk ) P(3), the coordinates of the point.
!
!    Output, real ( kind = rk ) DIST_SIGNED, the signed distance from 
!    the point to the plane.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) d
  real ( kind = rk ) dist_signed
  real ( kind = rk ) norm
  real ( kind = rk ) p(dim_num)

  norm = sqrt ( a * a + b * b + c * c )

  if ( norm == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLANE_IMP_POINT_DIST_SIGNED_3D - Fatal error!'
    write ( *, '(a)' ) '  The plane normal vector is null.'
    stop 1
  end if

  dist_signed = - sign ( 1.0D+00, d ) &
    * ( a * p(1) + b * p(2) + c * p(3) + d ) / norm

  return
end
subroutine plane_imp_point_near_3d ( a, b, c, d, p, pn )

!*****************************************************************************80
!
!! PLANE_IMP_POINT_NEAR_3D: nearest point on a implicit plane to a point in 3D.
!
!  Discussion:
!
!    The implicit form of a plane in 3D is:
!
!      A * X + B * Y + C * Z + D = 0
!
!    A normal vector to the plane is (A,B,C).
!
!    The line defined by (XN-P(1))/A = (YN-P(2))/B = (ZN-P(3))/C = T
!    goes through P and is parallel to N.
!
!    Solving for the point (XN,YN,ZN) we get
!
!      XN = A*T+P(1)
!      YN = B*T+P(2)
!      ZN = C*T+P(3)
!
!    Now place these values in the equation for the plane:
!
!      A*(A*T+P(1)) + B*(B*T+P(2)) + C*(C*T+P(3)) + D = 0
!
!    and solve for T:
!
!      T = (-A*P(1)-B*P(2)-C*P(3)-D) / (A * A + B * B + C * C )
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, C, D, the implicit plane parameters.
!
!    Input, real ( kind = rk ) P(3), the coordinates of the point.
!
!    Output, real ( kind = rk ) PN(3), the nearest point on the plane.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) d
  real ( kind = rk ) p(dim_num)
  logical plane_imp_is_degenerate_3d
  real ( kind = rk ) pn(dim_num)
  real ( kind = rk ) t

  if ( plane_imp_is_degenerate_3d ( a, b, c ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLANE_IMP_POINT_NEAR_3D - Fatal error!'
    write ( *, '(a)' ) '  A = B = C = 0.'
    stop 1
  end if

  t = - ( a * p(1) + b * p(2) + c * p(3) + d ) / ( a * a + b * b + c * c )

  pn(1) = p(1) + a * t
  pn(2) = p(2) + b * t
  pn(3) = p(3) + c * t

  return
end
subroutine plane_imp_segment_near_3d ( p1, p2, a, b, c, d, dist, p, pn )

!*****************************************************************************80
!
!! PLANE_IMP_SEGMENT_NEAR_3D: nearest ( implicit plane, line segment ) in 3D.
!
!  Discussion:
!
!    The implicit form of a plane in 3D is:
!
!      A * X + B * Y + C * Z + D = 0
!
!    A line segment is the finite portion of a line that lies between
!    two points.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(3), P2(3), the endpoints of the line
!    segment.
!
!    Input, real ( kind = rk ) A, B, C, D, the implicit plane parameters.
!
!    Output, real ( kind = rk ) DIST, the distance between the line segment and
!    the plane.
!
!    Output, real ( kind = rk ) P(3), the nearest point on the plane.
!
!    Output, real ( kind = rk ) PN(3), the nearest point on the line
!    segment to the plane.  If DIST is zero, the PN is a point of
!    intersection of the plane and the line segment.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) a
  real ( kind = rk ) alpha
  real ( kind = rk ) an
  real ( kind = rk ) b
  real ( kind = rk ) bn
  real ( kind = rk ) c
  real ( kind = rk ) cn
  real ( kind = rk ) d
  real ( kind = rk ) dist
  real ( kind = rk ) dn
  real ( kind = rk ) dot1
  real ( kind = rk ) dot2
  real ( kind = rk ) norm
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) pn(dim_num)

  pn(1:dim_num) = 0.0D+00
  p(1:dim_num) = 0.0D+00

  norm = sqrt ( a * a + b * b + c * c )

  if ( norm == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLANE_IMP_SEGMENT_NEAR_3D - Fatal error!'
    write ( *, '(a)' ) '  Plane normal vector is null.'
    stop 1
  end if
!
!  The normalized coefficients allow us to compute the (signed) distance.
!
  an = a / norm
  bn = b / norm
  cn = c / norm
  dn = d / norm
!
!  If the line segment is actually a point, then the answer is easy.
!
  if ( all ( p1(1:dim_num) == p2(1:dim_num) ) ) then

    dot1 = an * p1(1) + bn * p1(2) + cn * p1(3) + dn
    dist = abs ( dot1 )
    pn(1:dim_num) = p1(1:dim_num)
    p(1) = pn(1) - an * dot1
    p(2) = pn(2) - bn * dot1
    p(3) = pn(3) - cn * dot1
    return

  end if
!
!  Compute the projections of the two points onto the normal vector.
!
  dot1 = an * p1(1) + bn * p1(2) + cn * p1(3) + dn
  dot2 = an * p2(1) + bn * p2(2) + cn * p2(3) + dn
!
!  If these have the same sign, then the line segment does not
!  cross the plane, and one endpoint is the nearest point.
!
  if ( ( 0.0D+00 < dot1 .and. 0.0D+00 < dot2 ) .or. &
       ( dot1 < 0.0D+00 .and. dot2 < 0.0D+00 ) ) then

    dot1 = abs ( dot1 )
    dot2 = abs ( dot2 )

    if ( dot1 < dot2 ) then
      pn(1:dim_num) = p1(1:dim_num)
      p(1) = pn(1) - an * dot1
      p(2) = pn(2) - bn * dot1
      p(3) = pn(3) - cn * dot1
      dist = dot1
    else
      pn(1:dim_num) = p2(1:dim_num)
      dist = dot2
      p(1) = pn(1) - an * dot2
      p(2) = pn(2) - bn * dot2
      p(3) = pn(3) - cn * dot2
    end if
!
!  If the projections differ in sign, the line segment crosses the plane.
!
  else

    if ( dot1 == 0.0D+00 ) then
      alpha = 0.0D+00
    else if ( dot2 == 0.0D+00 ) then
      alpha = 1.0D+00
    else
      alpha = dot2 / ( dot2 - dot1 )
    end if

    pn(1:dim_num) =             alpha   * p1(1:dim_num) &
                  + ( 1.0D+00 - alpha ) * p2(1:dim_num)

    p(1:dim_num) = pn(1:dim_num)

    dist = 0.0D+00

  end if
 
  return
end
subroutine plane_imp_triangle_int_3d ( a, b, c, d, t, int_num, pint )

!*****************************************************************************80
!
!! PLANE_IMP_TRIANGLE_INT_3D: intersection ( implicit plane, triangle ) in 3D.
!
!  Discussion:
!
!    The implicit form of a plane in 3D is:
!
!      A * X + B * Y + C * Z + D = 0
!
!    There may be 0, 1, 2 or 3 points of intersection returned.
!
!    If two intersection points are returned, then the entire line
!    between them comprises points of intersection.
!
!    If three intersection points are returned, then all points of
!    the triangle intersect the plane.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, C, D, the implicit plane parameters.
!
!    Input, real ( kind = rk ) T(3,3), the vertices of the triangle.
!
!    Output, integer INT_NUM, the number of intersection points 
!    returned.
!
!    Output, real ( kind = rk ) PINT(3,3), the intersection points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) d
  real ( kind = rk ) dist1
  real ( kind = rk ) dist2
  real ( kind = rk ) dist3
  integer int_num
  real ( kind = rk ) pint(dim_num,3)
  real ( kind = rk ) t(dim_num,3)

  int_num = 0
!
!  Compute the signed distances between the vertices and the plane.
!
  dist1 = a * t(1,1) + b * t(2,1) + c * t(3,1) + d
  dist2 = a * t(1,2) + b * t(2,2) + c * t(3,2) + d
  dist3 = a * t(1,3) + b * t(2,3) + c * t(3,3) + d
!
!  Consider any zero distances.
!
  if ( dist1 == 0.0D+00 ) then
    int_num = int_num + 1
    pint(1:dim_num,int_num) = t(1:dim_num,1)
  end if

  if ( dist2 == 0.0D+00 ) then
    int_num = int_num + 1
    pint(1:dim_num,int_num) = t(1:dim_num,2)
  end if

  if ( dist3 == 0.0D+00 ) then
    int_num = int_num + 1
    pint(1:dim_num,int_num) = t(1:dim_num,3)
  end if
!
!  If 2 or 3 of the nodes intersect, we're already done.
!
  if ( 2 <= int_num ) then
    return
  end if
!
!  If one node intersects, then we're done unless the other two
!  are of opposite signs.
!
  if ( int_num == 1 ) then

    if ( dist1 == 0.0D+00 ) then

      call plane_imp_triangle_int_add_3d ( t(1:dim_num,2), t(1:dim_num,3), &
        dist2, dist3, int_num, pint )

    else if ( dist2 == 0.0D+00 ) then

      call plane_imp_triangle_int_add_3d ( t(1:dim_num,1), t(1:dim_num,3), &
        dist1, dist3, int_num, pint )

    else if ( dist3 == 0.0D+00 ) then

      call plane_imp_triangle_int_add_3d ( t(1:dim_num,1), t(1:dim_num,2), &
        dist1, dist2, int_num, pint )

    end if

    return

  end if
!
!  All nodal distances are nonzero, and there is at least one
!  positive and one negative.
!
  if ( dist1 * dist2 < 0.0D+00 .and. dist1 * dist3 < 0.0D+00 ) then

    call plane_imp_triangle_int_add_3d ( t(1:dim_num,1), t(1:dim_num,2), &
      dist1, dist2, int_num, pint )

    call plane_imp_triangle_int_add_3d ( t(1:dim_num,1), t(1:dim_num,3), &
      dist1, dist3, int_num, pint )

  else if ( dist2 * dist1 < 0.0D+00 .and. dist2 * dist3 < 0.0D+00 ) then

    call plane_imp_triangle_int_add_3d ( t(1:dim_num,2), t(1:dim_num,1), &
      dist2, dist1, int_num, pint )

    call plane_imp_triangle_int_add_3d ( t(1:dim_num,2), t(1:dim_num,3), &
      dist2, dist3, int_num, pint )

  else if ( dist3 * dist1 < 0.0D+00 .and. dist3 * dist2 < 0.0D+00 ) then

    call plane_imp_triangle_int_add_3d ( t(1:dim_num,3), t(1:dim_num,1), &
      dist3, dist1, int_num, pint )

    call plane_imp_triangle_int_add_3d ( t(1:dim_num,3), t(1:dim_num,2), &
      dist3, dist2, int_num, pint )

  end if

  return
end
subroutine plane_imp_triangle_int_add_3d ( p1, p2, dist1, dist2, int_num, pint )

!*****************************************************************************80
!
!! PLANE_IMP_TRIANGLE_INT_ADD_3D is a utility for plane/triangle intersections.
!
!  Discussion:
!
!    This routine is called to consider the value of the signed distance
!    from a plane of two nodes of a triangle.  If the two values
!    have opposite signs, then there is a point of intersection between
!    them.  The routine computes this point and adds it to the list.
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
!  Parameters:
!
!    Input, real ( kind = rk ) P1(3), P2(3), the coordinates of two vertices 
!    of a triangle.
!
!    Input, real ( kind = rk ) DIST1, DIST2, the signed distances of the 
!    two vertices from a plane.
!
!    Input/output, integer INT_NUM, the number of intersection 
!    points.
!
!    Input/output, real ( kind = rk ) PINT(3,INT_NUM), the intersection points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) alpha
  real ( kind = rk ) dist1
  real ( kind = rk ) dist2
  integer int_num
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) pint(dim_num,3)

  if ( dist1 == 0.0D+00 ) then
    int_num = int_num + 1
    pint(1:dim_num,int_num) = p1(1:dim_num)
  else if ( dist2 == 0.0D+00 ) then
    int_num = int_num + 1
    pint(1:dim_num,int_num) = p2(1:dim_num)
  else if ( dist1 * dist2 < 0.0D+00 ) then
    alpha = dist2 / ( dist2 - dist1 )
    int_num = int_num + 1
    pint(1:dim_num,int_num) =             alpha   * p1(1:dim_num) &
                            + ( 1.0D+00 - alpha ) * p2(1:dim_num)
  end if

  return
end
subroutine plane_imp_triangle_near_3d ( t, a, b, c, d, dist, near_num, pn )

!*****************************************************************************80
!
!! PLANE_IMP_TRIANGLE_NEAR_3D: nearest ( implicit plane, triangle ) in 3D.
!
!  Discussion:
!
!    The implicit form of a plane in 3D is:
!
!      A * X + B * Y + C * Z + D = 0
!
!    If DIST = 0, then each point is a point of intersection, and there
!    will be at most 3 such points returned.
!
!    If 0 < DIST, then the points are listed in pairs, with the first
!    being on the triangle, and the second on the plane.  Two points will
!    be listed in the most common case, but possibly 4 or 6.
!
!    I should see to it that the underlying distance routine always returns
!    one of the endpoints if the entire line segment is at zero distance.
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
!  Parameters:
!
!    Input, real ( kind = rk ) T(3,3), the vertices of the triangle.
!
!    Input, real ( kind = rk ) A, B, C, D, the implicit plane parameters.
!
!    Output, real ( kind = rk ) DIST, the distance between the triangle
!    and the plane.
!
!    Output, integer NEAR_NUM, the number of nearest points 
!    returned.
!
!    Output, real ( kind = rk ) PN(3,6), a collection of nearest points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) d
  real ( kind = rk ) dist
  real ( kind = rk ) dist12
  real ( kind = rk ) dist23
  real ( kind = rk ) dist31
  integer near_num
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) pn(dim_num,6)
  real ( kind = rk ) pt(dim_num)
  real ( kind = rk ) t(dim_num,3)

  near_num = 0
!
!  Consider the line segment P1 - P2.
!
  call plane_imp_segment_near_3d ( t(1:dim_num,1), t(1:dim_num,2), &
    a, b, c, d, dist12, p, pt )

  dist = dist12

  near_num = near_num + 1
  pn(1:dim_num,near_num) = pt(1:dim_num)

  if ( 0.0D+00 < dist12 ) then
    near_num = near_num + 1
    pn(1:dim_num,near_num) = p(1:dim_num)
  end if
!
!  Consider the line segment P2 - P3.
!
  call plane_imp_segment_near_3d ( t(1:dim_num,2), t(1:dim_num,3), &
    a, b, c, d, dist23, p, pt )

  if ( dist23 < dist ) then

    near_num = 0
    dist = dist23

    near_num = near_num + 1
    pn(1:dim_num,near_num) = pt(1:dim_num)

    if ( 0.0D+00 < dist23 ) then
      near_num = near_num + 1
      pn(1:dim_num,near_num) = p(1:dim_num)
    end if

  else if ( dist23 == dist ) then

    near_num = near_num + 1
    pn(1:dim_num,near_num) = pt(1:dim_num)

    if ( 0.0D+00 < dist23 ) then
      near_num = near_num + 1
      pn(1:dim_num,near_num) = p(1:dim_num)
    end if

  end if
!
!  Consider the line segment P3 - P1.
!
  call plane_imp_segment_near_3d ( t(1:dim_num,3), t(1:dim_num,1), &
    a, b, c, d, dist31, p, pt )

  if ( dist31 < dist ) then

    near_num = 0
    dist = dist31

    near_num = near_num + 1
    pn(1:dim_num,near_num) = pt(1:dim_num)

    if ( 0.0D+00 < dist31 ) then
      near_num = near_num + 1
      pn(1:dim_num,near_num) = p(1:dim_num)
    end if

  else if ( dist31 == dist ) then

    near_num = near_num + 1
    pn(1:dim_num,near_num) = pt(1:dim_num)

    if ( 0.0D+00 < dist31 ) then
      near_num = near_num + 1
      pn(1:dim_num,near_num) = p(1:dim_num)
    end if

  end if

  return
end
subroutine plane_imp2exp_3d ( a, b, c, d, p1, p2, p3 )

!*****************************************************************************80
!
!! PLANE_IMP2EXP_3D converts an implicit plane to explicit form in 3D.
!
!  Discussion:
!
!    The implicit form of a plane in 3D is
!
!      A * X + B * Y + C * Z + D = 0.
!
!    The explicit form of a plane in 3D is
!
!      the plane through P1, P2 and P3.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, C, D, the implicit plane parameters.
!
!    Output, real ( kind = rk ) P1(3), P2(3), P3(3), three points on the plane.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) d
  real ( kind = rk ) normal(dim_num)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) p3(dim_num)
  real ( kind = rk ) pp(dim_num)

  call plane_imp2normal_3d ( a, b, c, d, pp, normal )

  call plane_normal2exp_3d ( pp, normal, p1, p2, p3 )

  return
end
subroutine plane_imp2normal_3d ( a, b, c, d, pp, normal )

!*****************************************************************************80
!
!! PLANE_IMP2NORMAL_3D converts an implicit plane to normal form in 3D.
!
!  Discussion:
!
!    The implicit form of a plane in 3D is
!
!      A * X + B * Y + C * Z + D = 0.
!
!    The normal form of a plane in 3D is
!
!      PP, a point on the plane, and
!      N, the unit normal to the plane.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, C, D, the implicit plane parameters.
!
!    Output, real ( kind = rk ) PP(3), a point on the plane.
!
!    Output, real ( kind = rk ) NORMAL(3), the unit normal vector to the plane.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) d
  real ( kind = rk ) norm
  real ( kind = rk ) normal(dim_num)
  real ( kind = rk ) pp(dim_num)

  norm = sqrt ( a * a + b * b + c * c )

  if ( norm == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLANE_IMP2NORMAL_3D - Fatal error!'
    write ( *, '(a)' ) '  The plane (A,B,C) has zero norm.'
    stop 1
  end if

  normal(1) = a / norm
  normal(2) = b / norm
  normal(3) = c / norm

  if ( a /= 0.0D+00 ) then
    pp(1) = - d / a
    pp(2) = 0.0D+00
    pp(3) = 0.0D+00
  else if ( b /= 0.0D+00 ) then
    pp(1) = 0.0D+00
    pp(2) = - d / b
    pp(3) = 0.0D+00
  else if ( c /= 0.0D+00 ) then
    pp(1) = 0.0D+00
    pp(2) = 0.0D+00
    pp(3) = - d / c
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLANE_IMP2NORMAL_3D - Fatal error!'
    write ( *, '(a)' ) '  The (A,B,C) vector is null.'
    stop 1
  end if

  return
end
subroutine plane_normal_basis_3d ( pp, normal, pq, pr )

!*****************************************************************************80
!
!! PLANE_NORMAL_BASIS_3D finds two perpendicular vectors in a plane in 3D.
!
!  Discussion:
!
!    The normal form of a plane in 3D is:
!
!      PP is a point on the plane,
!      N is a normal vector to the plane.
!
!    The two vectors to be computed, PQ and PR, can be regarded as
!    the basis of a Cartesian coordinate system for points in the plane.
!    Any point in the plane can be described in terms of the "origin" 
!    point PP plus a weighted sum of the two vectors PQ and PR:
!
!      P = PP + a * PQ + b * PR.
!
!    The vectors PQ and PR have unit length, and are perpendicular to N
!    and to each other.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) PP(3), a point on the plane.  (Actually,
!    we never need to know these values to do the calculation!)
!
!    Input, real ( kind = rk ) NORMAL(3), a normal vector N to the plane.  The
!    vector must not have zero length, but it is not necessary for N
!    to have unit length.
!
!    Output, real ( kind = rk ) PQ(3), a vector of unit length,
!    perpendicular to the vector N and the vector PR.
!
!    Output, real ( kind = rk ) PR(3), a vector of unit length,
!    perpendicular to the vector N and the vector PQ.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) r8vec_norm
  real ( kind = rk ) normal(dim_num)
  real ( kind = rk ) normal_norm
  real ( kind = rk ) pp(dim_num)
  real ( kind = rk ) pq(dim_num)
  real ( kind = rk ) pr(dim_num)
  real ( kind = rk ) pr_norm

  call r8_fake_use ( pp(1) )
!
!  Compute the length of NORMAL.
!
  normal_norm = r8vec_norm ( dim_num, normal )

  if ( normal_norm == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLANE_NORMAL_BASIS_3D - Fatal error!'
    write ( *, '(a)' ) '  The normal vector is 0.'
    stop 1
  end if
!
!  Find a vector PQ that is normal to NORMAL and has unit length.
!
  call r8vec_any_normal ( dim_num, normal, pq )
!
!  Now just take the cross product NORMAL x PQ to get the PR vector.
!
  call r8vec_cross_product_3d ( normal, pq, pr )

  pr_norm = r8vec_norm ( dim_num, pr )

  pr(1:dim_num) = pr(1:dim_num) / pr_norm

  return
end
subroutine plane_normal_line_exp_int_3d ( pp, normal, p1, p2, ival, pint )

!*****************************************************************************80
!
!! PLANE_NORMAL_LINE_EXP_INT_3D: intersection of plane and line in 3D.
!
!  Discussion:
!
!    The normal form of a plane in 3D is:
!
!      PP is a point on the plane,
!      N is a normal vector to the plane.
!
!    The explicit form of a line in 3D is:
!
!      P1, P2 are two points on the line.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) PP(3), a point on the plane.
!
!    Input, real ( kind = rk ) NORMAL(3), a normal vector to the plane.
!
!    Input, real ( kind = rk ) P1(3), P2(3), two distinct points on the line.
!
!    Output, integer IVAL, the kind of intersection;
!    0, the line and plane seem to be parallel and separate;
!    1, the line and plane intersect at a single point;
!    2, the line and plane seem to be parallel and joined.
!
!    Output, real ( kind = rk ) PINT(3), the coordinates of a
!    common point of the plane and line, when IVAL is 1 or 2.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) direction(dim_num)
  integer ival
  logical line_exp_is_degenerate_nd
  real ( kind = rk ) normal(dim_num)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) pint(dim_num)
  real ( kind = rk ) pp(dim_num)
  real ( kind = rk ) temp
!
!  Make sure the line is not degenerate.
!
  if ( line_exp_is_degenerate_nd ( dim_num, p1, p2 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLANE_NORMAL_LINE_EXP_INT_3D - Fatal error!'
    write ( *, '(a)' ) '  The line is degenerate.'
    stop 1
  end if
!
!  Make sure the plane normal vector is a unit vector.
!
  temp = sqrt ( sum ( normal(1:dim_num) ** 2 ) )

  if ( temp == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLANE_NORMAL_LINE_EXP_INT_3D - Fatal error!'
    write ( *, '(a)' ) '  The normal vector of the plane is degenerate.'
    stop 1
  end if

  normal(1:dim_num) = normal(1:dim_num) / temp
!
!  Determine the unit direction vector of the line.
!
  direction(1:dim_num) = p2(1:dim_num) - p1(1:dim_num)
  temp = sqrt ( sum ( direction(1:dim_num) ** 2 ) )
  direction(1:dim_num) = direction(1:dim_num) / temp
!
!  If the normal and direction vectors are orthogonal, then
!  we have a special case to deal with.
!
  if ( dot_product ( normal(1:dim_num), direction(1:dim_num) ) == 0.0D+00 ) then

    temp = dot_product ( normal(1:dim_num), p1(1:dim_num) - pp(1:dim_num) )

    if ( temp == 0.0D+00 ) then
      ival = 2
      pint(1:dim_num) = p1(1:dim_num)
    else
      ival = 0
      pint(1:dim_num) = huge ( temp )
    end if

    return
  end if
!
!  Determine the distance along the direction vector to the intersection point.
!
  temp = dot_product ( pp(1:dim_num) - p1(1:dim_num), normal(1:dim_num) ) &
    / dot_product ( direction(1:dim_num), normal(1:dim_num) )

  ival = 1
  pint(1:dim_num) = p1(1:dim_num) + temp * direction(1:dim_num)

  return
end
subroutine plane_normal_qr_to_xyz ( pp, normal, pq, pr, n, qr, xyz )

!*****************************************************************************80
!
!! PLANE_NORMAL_QR_TO_XYZ: QR_TO_XYZ coordinates for a normal form plane.
!
!  Discussion:
!
!    The normal form of a plane in 3D is:
!
!      PP is a point on the plane,
!      NORMAL is a normal vector to the plane.
!
!    Two vectors PQ and PR can be computed with the properties that
!    * NORMAL, PQ and PR are pairwise orthogonal;
!    * PQ and PR have unit length;
!    * every point P in the plane has a "QR" representation
!      as P = PP + q * PQ + r * PR.
!
!    This function is given the QR coordinates of a set of points on the
!    plane, and returns the XYZ coordinates.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 November 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) PP(3), a point on the plane.
!
!    Input, real ( kind = rk ) NORMAL(3), a normal vector N to the plane.  The
!    vector must not have zero length, but it is not necessary for N
!    to have unit length.
!
!    Input, real ( kind = rk ) PQ(3), a vector of unit length,
!    perpendicular to the vector N and the vector PR.
!
!    Input, real ( kind = rk ) PR(3), a vector of unit length,
!    perpendicular to the vector N and the vector PQ.
!
!    Input, integer N, the number of points on the plane.
!
!    Input, real ( kind = rk ) QR(2,N), the QR coordinates of the points.
!
!    Output, real ( kind = rk ) XYZ(3,N), the XYZ coordinates of the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )
  
  integer n

  real ( kind = rk ) normal(3)
  real ( kind = rk ) pp(3)
  real ( kind = rk ) pq(3)
  real ( kind = rk ) pqpr(3,2)
  real ( kind = rk ) pr(3)
  real ( kind = rk ) qr(2,n)
  real ( kind = rk ) xyz(3,n)

  call r8_fake_use ( normal(1) )

  xyz(1,1:n) = pp(1)
  xyz(2,1:n) = pp(2)
  xyz(3,1:n) = pp(3)

  pqpr(1:3,1) = pq(1:3)
  pqpr(1:3,2) = pr(1:3)

  xyz(1:3,1:n) = xyz(1:3,1:n) + matmul ( pqpr(1:3,1:2), qr(1:2,1:n) )

  return
end
subroutine plane_normal_tetrahedron_intersect ( pp, normal, t, int_num, pint )

!*****************************************************************************80
!
!! PLANE_NORMAL_TETRAHEDRON_INTERSECT intersects a plane and a tetrahedron.
!
!  Discussion:
!
!    The intersection of a plane and a tetrahedron is one of:
!    0) empty
!    1) a single point
!    2) a single line segment
!    3) a triangle
!    4) a quadrilateral.
!
!    In each case, the region of intersection can be described by the
!    corresponding number of points.  In particular, cases 2, 3 and 4
!    are described by the vertices that bound the line segment, triangle,
!    or quadrilateral.
!
!    The normal form of a plane is:
!
!      PP is a point on the plane,
!      N is a normal vector to the plane.
!
!    The form of a tetrahedron is
!
!      T(1:3,1:4) contains the coordinates of the vertices.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) PP(3), a point on the plane.
!
!    Input, real ( kind = rk ) NORMAL(3), a normal vector to the plane.
!
!    Input, real ( kind = rk ) T(3,4), the tetrahedron vertices.
!
!    Output, integer INT_NUM, the number of intersection 
!    points returned.  This will be 0, 1, 2, 3 or 4.
!
!    Output, real ( kind = rk ) PINT(3,4), the coordinates of the
!    intersection points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) area1
  real ( kind = rk ) area2
  real ( kind = rk ) d(4)
  real ( kind = rk ) dn
  real ( kind = rk ) dpp
  integer int_num
  integer j1
  integer j2
  real ( kind = rk ) normal(3)
  real ( kind = rk ) pint(3,4)
  real ( kind = rk ) pp(3)
  logical r8_sign_opposite_strict
  real ( kind = rk ) t(3,4)
  real ( kind = rk ) temp(3)

  int_num = 0
  pint(1:3,1:4) = 0.0D+00
!
!  DN is the length of the normal vector.
!
  dn = dot_product ( normal(1:3), normal(1:3) )
!
!  DPP is the distance between the origin and the projection of the
!  point PP onto the normal vector.
!
  dpp = dn - dot_product ( normal(1:3), pp(1:3) )
!
!  D(I) is positive, zero, or negative if vertex I is above,
!  on, or below the plane.
!
  d(1:4) = dn - matmul ( normal(1:3), t(1:3,1:4) ) - dpp
!
!  If all D are positive or negative, no intersection.
!
  if ( all ( d(1:4) < 0.0D+00 ) .or. all ( 0.0D+00 < d(1:4) ) ) then
    int_num = 0
    return
  end if
!
!  Points with zero distance are automatically added to the list.
!
!  For each point with nonzero distance, seek another point
!  with opposite sign and higher index, and compute the intersection
!  of the line between those points and the plane.
!
  do j1 = 1, 4

    if ( d(j1) == 0.0D+00 ) then
      int_num = int_num + 1
      pint(1:3,int_num) = t(1:3,j1)
    else
      do j2 = j1 + 1, 4
        if ( r8_sign_opposite_strict ( d(j1), d(j2) ) ) then
          int_num = int_num + 1
          pint(1:3,int_num) = ( d(j1)         * t(1:3,j2)   &
                                      - d(j2) * t(1:3,j1) ) &
                            / ( d(j1) - d(j2) )
        end if
      end do
    end if
  end do
!
!  If four points were found, try to order them properly.
!
  if ( int_num == 4 ) then
    call quad_area_3d ( pint, area1 )
    temp(1:3) = pint(1:3,3)
    pint(1:3,3) = pint(1:3,4)
    pint(1:3,4) = temp(1:3)
    call quad_area_3d ( pint, area2 )
    if ( area2 < area1 ) then
      temp(1:3) = pint(1:3,3)
      pint(1:3,3) = pint(1:3,4)
      pint(1:3,4) = temp(1:3)
    end if
  end if

  return
end
subroutine plane_normal_triangle_int_3d ( pp, normal, t, int_num, pint )

!*****************************************************************************80
!
!! PLANE_NORMAL_TRIANGLE_INT_3D: intersection ( normal plane, triangle ) in 3D.
!
!  Discussion:
!
!    The normal form of a plane in 3D is:
!
!      PP is a point on the plane,
!      N is a normal vector to the plane.
!
!    There may be 0, 1, 2 or 3 points of intersection returned.
!
!    If two intersection points are returned, then the entire line
!    between them comprises points of intersection.
!
!    If three intersection points are returned, then all points of
!    the triangle intersect the plane.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) PP(3), a point on the plane.
!
!    Input, real ( kind = rk ) NORMAL(3), a normal vector to the plane.
!
!    Input, real ( kind = rk ) T(3,3), the vertices of the triangle.
!
!    Output, integer INT_NUM, the number of intersection 
!    points returned.
!
!    Output, real ( kind = rk ) PINT(3,3), the coordinates of the
!    intersection points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) d
  real ( kind = rk ) dist1
  real ( kind = rk ) dist2
  real ( kind = rk ) dist3
  real ( kind = rk ) normal(dim_num)
  integer int_num
  real ( kind = rk ) pint(dim_num,3)
  real ( kind = rk ) pp(dim_num)
  real ( kind = rk ) t(dim_num,3)

  int_num = 0
!
!  Compute the signed distances between the vertices and the plane.
!
  d = - dot_product ( normal(1:dim_num), pp(1:dim_num) )

  dist1 = dot_product ( normal(1:dim_num), t(1:dim_num,1) ) + d
  dist2 = dot_product ( normal(1:dim_num), t(1:dim_num,2) ) + d
  dist3 = dot_product ( normal(1:dim_num), t(1:dim_num,3) ) + d
!
!  Consider any zero distances.
!
  if ( dist1 == 0.0D+00 ) then

    int_num = int_num + 1
    pint(1:dim_num,int_num) = t(1:dim_num,1)

  end if

  if ( dist2 == 0.0D+00 ) then

    int_num = int_num + 1
    pint(1:dim_num,int_num) = t(1:dim_num,2)

  end if

  if ( dist3 == 0.0D+00 ) then

    int_num = int_num + 1
    pint(1:dim_num,int_num) = t(1:dim_num,3)

  end if
!
!  If 2 or 3 of the nodes intersect, we're already done.
!
  if ( 2 <= int_num ) then
    return
  end if
!
!  If one node intersects, then we're done unless the other two
!  are of opposite signs.
!
  if ( int_num == 1 ) then

    if ( dist1 == 0.0D+00 ) then

      call plane_imp_triangle_int_add_3d ( t(1:dim_num,2), t(1:dim_num,3), &
        dist2, dist3, int_num, pint )

    else if ( dist2 == 0.0D+00 ) then

      call plane_imp_triangle_int_add_3d ( t(1:dim_num,1), t(1:dim_num,3), &
        dist1, dist3, int_num, pint )

    else if ( dist3 == 0.0D+00 ) then

      call plane_imp_triangle_int_add_3d ( t(1:dim_num,1), t(1:dim_num,2), &
        dist1, dist2, int_num, pint )

    end if

    return

  end if
!
!  All nodal distances are nonzero, and there is at least one
!  positive and one negative.
!
  if ( dist1 * dist2 < 0.0D+00 .and. dist1 * dist3 < 0.0D+00 ) then

    call plane_imp_triangle_int_add_3d ( t(1:dim_num,1), t(1:dim_num,2), &
      dist1, dist2, int_num, pint )

    call plane_imp_triangle_int_add_3d ( t(1:dim_num,1), t(1:dim_num,3), &
      dist1, dist3, int_num, pint )

  else if ( dist2 * dist1 < 0.0D+00 .and. dist2 * dist3 < 0.0D+00 ) then

    call plane_imp_triangle_int_add_3d ( t(1:dim_num,2), t(1:dim_num,1), &
      dist2, dist1, int_num, pint )

    call plane_imp_triangle_int_add_3d ( t(1:dim_num,2), t(1:dim_num,3), &
      dist2, dist3, int_num, pint )

  else if ( dist3 * dist1 < 0.0D+00 .and. dist3 * dist2 < 0.0D+00 ) then

    call plane_imp_triangle_int_add_3d ( t(1:dim_num,3), t(1:dim_num,1), &
      dist3, dist1, int_num, pint )

    call plane_imp_triangle_int_add_3d ( t(1:dim_num,3), t(1:dim_num,2), &
      dist3, dist2, int_num, pint )

  end if

  return
end
subroutine plane_normal_uniform_3d ( pp, normal )

!*****************************************************************************80
!
!! PLANE_NORMAL_UNIFORM_3D generates a random normal plane in 3D.
!
!  Discussion:
!
!    The normal form of a plane is:
!
!      PP is a point on the plane,
!      N is a normal vector to the plane.
!
!    The point PP will be chosen at random inside the unit sphere.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) PP(3), a point on the plane.
!
!    Output, real ( kind = rk ) NORMAL(3), the unit normal vector.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) norm
  real ( kind = rk ) normal(dim_num)
  real ( kind = rk ) pp(dim_num)
!
!  Pick PP as a random point inside the unit sphere in ND.
!
  call ball01_sample_3d ( pp )
!
!  Get values from a standard normal distribution.
!
  call r8vec_normal_01 ( dim_num, normal )
!
!  Compute the length of the vector.
!
  norm = sqrt ( sum ( normal(1:dim_num) ** 2 ) )
!
!  Normalize the vector.
!
  normal(1:dim_num) = normal(1:dim_num) / norm

  return
end
subroutine plane_normal_uniform_nd ( dim_num, pp, normal )

!*****************************************************************************80
!
!! PLANE_NORMAL_UNIFORM_ND generates a random normal plane in ND.
!
!  Discussion:
!
!    The normal form of a plane is:
!
!      PP is a point on the plane,
!      N is a normal vector to the plane.
!
!    The point PP will be chosen at random inside the unit sphere.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = rk ) PP(DIM_NUM), a point on the plane.
!
!    Output, real ( kind = rk ) NORMAL(DIM_NUM), the unit normal vector.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num

  real ( kind = rk ) norm
  real ( kind = rk ) normal(dim_num)
  real ( kind = rk ) pp(dim_num)
!
!  Pick PP as a random point inside the unit sphere in ND.
!
  call ball01_sample_nd ( dim_num, pp )
!
!  Get values from a standard normal distribution.
!
  call r8vec_normal_01 ( dim_num, normal )
!
!  Compute the length of the vector.
!
  norm = sqrt ( sum ( normal(1:dim_num) ** 2 ) )
!
!  Normalize the vector.
!
  normal(1:dim_num) = normal(1:dim_num) / norm

  return
end
subroutine plane_normal_xyz_to_qr ( pp, normal, pq, pr, n, xyz, qr )

!*****************************************************************************80
!
!! PLANE_NORMAL_XYZ_TO_QR: XYZ to QR coordinates for a normal form plane.
!
!  Discussion:
!
!    The normal form of a plane in 3D is:
!
!      PP is a point on the plane,
!      NORMAL is a normal vector to the plane.
!
!    Two vectors PQ and PR can be computed with the properties that
!    * NORMAL, PQ and PR are pairwise orthogonal;
!    * PQ and PR have unit length;
!    * every point P in the plane has a "QR" representation
!      as P = PP + q * PQ + r * PR.
!
!    This function is given the XYZ coordinates of a set of points on the
!    plane, and returns the QR coordinates.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 May 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) PP(3), a point on the plane.
!
!    Input, real ( kind = rk ) NORMAL(3), a normal vector N to the plane.  The
!    vector must not have zero length, but it is not necessary for N
!    to have unit length.
!
!    Input, real ( kind = rk ) PQ(3), a vector of unit length,
!    perpendicular to the vector N and the vector PR.
!
!    Input, real ( kind = rk ) PR(3), a vector of unit length,
!    perpendicular to the vector N and the vector PQ.
!
!    Input, integer N, the number of points on the plane.
!
!    Input, real ( kind = rk ) XYZ(3,N), the XYZ coordinates of the points.
!
!    Output, real ( kind = rk ) QR(2,N), the QR coordinates of the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )
  
  integer n

  real ( kind = rk ) normal(3)
  real ( kind = rk ) pp(3)
  real ( kind = rk ) pq(3)
  real ( kind = rk ) pr(3)
  real ( kind = rk ) qr(2,n)
  real ( kind = rk ) rpqp(2,3)
  real ( kind = rk ) xyz(3,n)

  call r8_fake_use ( normal(1) )

  rpqp(1,1:3) = pq(1:3)
  rpqp(2,1:3) = pr(1:3)

  qr(1:2,1:n) = matmul ( rpqp(1:2,1:3), xyz(1:3,1:n) )

  qr(1,1:n) = qr(1,1:n) - dot_product ( pq(1:3), pp(1:3) )
  qr(2,1:n) = qr(2,1:n) - dot_product ( pr(1:3), pp(1:3) )

  return
end
subroutine plane_normal2exp_3d ( pp, normal, p1, p2, p3 )

!*****************************************************************************80
!
!! PLANE_NORMAL2EXP_3D converts a normal plane to explicit form in 3D.
!
!  Discussion:
!
!    The normal form of a plane in 3D is
!
!      PP, a point on the plane, and
!      N, the unit normal to the plane.
!
!    The explicit form of a plane in 3D is
!
!      the plane through P1, P2 and P3.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) PP(3), a point on the plane.
!
!    Input, real ( kind = rk ) NORMAL(3), a normal vector N to the plane.  The
!    vector must not have zero length, but it is not necessary for N
!    to have unit length.
!
!    Output, real ( kind = rk ) P1(3), P2(3), P3(3), three points on the plane.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) normal(dim_num)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) p3(dim_num)
  real ( kind = rk ) pp(dim_num)
  real ( kind = rk ) pq(dim_num)
  real ( kind = rk ) pr(dim_num)

  call plane_normal_basis_3d ( pp, normal, pq, pr )

  p1(1:dim_num) = pp(1:dim_num)
  p2(1:dim_num) = pp(1:dim_num) + pq(1:dim_num)
  p3(1:dim_num) = pp(1:dim_num) + pr(1:dim_num)

  return
end
subroutine plane_normal2imp_3d ( pp, normal, a, b, c, d )

!*****************************************************************************80
!
!! PLANE_NORMAL2IMP_3D converts a normal form plane to implicit form in 3D.
!
!  Discussion:
!
!    The normal form of a plane in 3D is
!
!      PP, a point on the plane, and
!      N, the unit normal to the plane.
!
!    The implicit form of a plane in 3D is
!
!      A * X + B * Y + C * Z + D = 0.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) PP(3), a point on the plane.
!
!    Input, real ( kind = rk ) NORMAL(3), the unit normal vector to the plane.
!
!    Output, real ( kind = rk ) A, B, C, D, the implicit plane parameters.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) d
  real ( kind = rk ) normal(dim_num)
  real ( kind = rk ) pp(dim_num)

  a = normal(1)
  b = normal(2)
  c = normal(3)
  d = - a * pp(1) - b * pp(2) - c * pp(3)

  return
end
subroutine planes_imp_angle_3d ( a1, b1, c1, d1, a2, b2, c2, d2, angle )

!*****************************************************************************80
!
!! PLANES_IMP_ANGLE_3D: dihedral angle between implicit planes in 3D.
!
!  Discussion:
!
!    The implicit form of a plane in 3D is:
!
!      A * X + B * Y + C * Z + D = 0
!
!    If two planes P1 and P2 intersect in a nondegenerate way, then there is a
!    line of intersection L0.  Consider any plane perpendicular to L0.  The
!    dihedral angle of P1 and P2 is the angle between the lines L1 and L2, where
!    L1 is the intersection of P1 and P0, and L2 is the intersection of P2 
!    and P0.
!
!    The dihedral angle may also be calculated as the angle between the normal
!    vectors of the two planes.  Note that if the planes are parallel or
!    coincide, the normal vectors are identical, and the dihedral angle is 0.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    08 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Daniel Zwillinger, editor,
!    CRC Standard Math Tables and Formulae, 30th edition,
!    Section 4.13, "Planes",
!    CRC Press, 1996, pages 305-306.
!
!  Parameters:
!
!    Input, real ( kind = rk ) A1, B1, C1, D1, coefficients that define the
!    first plane.
!
!    Input, real ( kind = rk ) A2, B2, C2, D2, coefficients that define
!    the second plane.
!
!    Output, real ( kind = rk ) ANGLE, the dihedral angle, in radians,
!    defined by the two planes.  If either plane is degenerate, or they
!    do not intersect, or they coincide, then the angle is set to HUGE(1.0).
!    Otherwise, the angle is between 0 and PI.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) a1
  real ( kind = rk ) a2
  real ( kind = rk ) angle
  real ( kind = rk ) b1
  real ( kind = rk ) b2
  real ( kind = rk ) c1
  real ( kind = rk ) c2
  real ( kind = rk ) cosine
  real ( kind = rk ) d1
  real ( kind = rk ) d2
  real ( kind = rk ) norm1
  real ( kind = rk ) norm2
  real ( kind = rk ) r8_acos

  call r8_fake_use ( d1 )
  call r8_fake_use ( d2 )

  norm1 = sqrt ( a1 * a1 + b1 * b1 + c1 * c1 )

  if ( norm1 == 0.0D+00 ) then
    angle = huge ( angle )
    return
  end if

  norm2 = sqrt ( a2 * a2 + b2 * b2 + c2 * c2 )

  if ( norm2 == 0.0D+00 ) then
    angle = huge ( angle )
    return
  end if

  cosine = ( a1 * a2 + b1 * b2 + c1 * c2 ) / ( norm1 * norm2 )

  angle = r8_acos ( cosine )

  return
end
function points_avoid_point_naive_2d ( n, p_set, p )

!*****************************************************************************80
!
!! POINTS_AVOID_POINT_NAIVE_2D: is a point "far" from a set of points in 2D?
!
!  Discussion:
!
!    The routine discards points that are too close to other points.
!    The method used to check this is quadratic in the number of points,
!    and may take an inordinate amount of time if there are a large
!    number of points.  But in that case, what do you want?  If you want
!    lots of points, you don't want to delete any because it won't matter.
!
!    The test point is "far enough" from an accepted point if
!    the Euclidean distance is at least 100 times EPSILON.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    24 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of accepted points.
!
!    Input, real ( kind = rk ) P_SET(2,N), the accepted points.
!
!    Input, real ( kind = rk ) P(2), a point to be tested.
!
!    Output, logical POINTS_AVOID_POINT_NAIVE_2D, is TRUE if 
!    XY_TEST is "far enough" from all the accepted points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: dim_num = 2

  integer j
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) p_set(dim_num,n)
  logical points_avoid_point_naive_2d
  real ( kind = rk ) tol

  tol = 100.0D+00 * epsilon ( tol )

  points_avoid_point_naive_2d = .true.

  do j = 1, n

    if ( sqrt ( sum ( ( p_set(1:dim_num,j) - p(1:dim_num) )**2 ) ) < tol ) then
      points_avoid_point_naive_2d = .false.
      return
    end if

  end do

  return
end
subroutine points_bisect_line_imp_2d ( p1, p2, a, b, c )

!*****************************************************************************80
!
!! POINTS_BISECT_LINE_IMP_2D: implicit bisector line between two points in 2D.
!
!  Discussion:
!
!    This routine finds, in implicit form, the equation of the line
!    that is equidistant from two points.
!
!    The implicit form of a line in 2D is:
!
!      A * X + B * Y + C = 0
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 January 1999
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
!  Parameters:
!
!    Input, real ( kind = rk ) P1(2), P2(2), the coordinates of two points.
!
!    Output, real ( kind = rk ) A, B, C, the parameters of the implicit line
!    equidistant from both points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)

  a = p1(1) - p2(1)
  b = p1(2) - p2(2)
  c = - 0.5D+00 * ( ( p1(1) * p1(1) + p1(2) * p1(2) ) &
                  - ( p2(1) * p2(1) + p2(2) * p2(2) ) )

  return
end
subroutine points_bisect_line_par_2d ( p1, p2, f, g, x, y )

!*****************************************************************************80
!
!! POINTS_BISECT_LINE_PAR_2D: parametric bisector line between points in 2D.
!
!  Discussion:
!
!    This routine finds, in parametric form, the equation of the line
!    that is equidistant from two points.
!
!    The parametric form of a line in 2D is:
!
!      X = X0 + F * T
!      Y = Y0 + G * T
!
!    We normalize by always choosing F*F + G*G = 1, and F nonnegative.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 January 2005
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
!  Parameters:
!
!    Input, real ( kind = rk ) P1(2), P2(2), two points.
!
!    Output, real ( kind = rk ) F, G, X, Y, the parameters of the parametric line
!    equidistant from both points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) f
  real ( kind = rk ) g
  real ( kind = rk ) norm
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) x
  real ( kind = rk ) y

  f = 0.5D+00 * ( p1(1) + p2(1) )
  g = 0.5D+00 * ( p1(2) + p2(2) )

  norm = f * f + g * g

  if ( norm /= 0.0D+00 ) then
    f = f / norm
    g = g / norm
  end if

  if ( f < 0.0D+00 ) then
    f = -f
    g = -g
  end if

  x = - ( p2(2) - p1(2) )
  y = + ( p2(1) - p1(1) )

  return
end
subroutine points_centroid_2d ( n, p, centroid_index )

!*****************************************************************************80
!
!! POINTS_CENTROID_2D computes the discrete centroid of a point set in 2D.
!
!  Discussion:
!
!    Given a discrete set of points S, the discrete centroid z is defined by
!
!                           sum ( x in S ) ( x - z )^2
!        = min ( y in S ) { sum ( x in S ) ( x - y )^2
!
!    In other words, the discrete centroid is a point in the set whose distance
!    to the other points is minimized.  The discrete centroid of a point set
!    need not be unique.  Consider a point set that comprises the
!    vertices of an equilateral triangle.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) P(2,N), the points.
!
!    Output, integer CENTROID_INDEX, the index of a discrete
!    centroid of the set, between 1 and N.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: dim_num = 2

  integer centroid_index
  real ( kind = rk ) dist
  real ( kind = rk ) dist_min
  integer i
  integer j
  real ( kind = rk ) p(dim_num,n)

  dist_min = 0.0D+00
  centroid_index = -1

  do i = 1, n

    dist = 0.0D+00
    do j = 1, n
      dist = dist + sum ( ( p(1:dim_num,i) - p(1:dim_num,j) )**2 )
    end do

    if ( i == 1 ) then
      dist_min = dist
      centroid_index = i
    else if ( dist < dist_min ) then
      dist_min = dist
      centroid_index = i
    end if

  end do

  return
end
subroutine points_colin_2d ( p1, p2, p3, colin )

!*****************************************************************************80
!
!! POINTS_COLIN_2D estimates the colinearity of 3 points in 2D.
!
!  Discussion:
!
!    The estimate of colinearity that is returned is the ratio
!    of the area of the triangle spanned by the points to the area
!    of the equilateral triangle with the same perimeter.
!
!    This estimate is 1 if the points are maximally noncolinear, 0 if the
!    points are exactly colinear, and otherwise is closer to 1 or 0 depending
!    on whether the points are far or close to colinearity.
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
!  Parameters:
!
!    Input, real ( kind = rk ) P1(2), P2(2), P3(2), the points.
!
!    Output, real ( kind = rk ) COLIN, the colinearity estimate.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) area_triangle
  real ( kind = rk ) area2
  real ( kind = rk ) colin
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) p3(dim_num)
  real ( kind = rk ) perim
  real ( kind = rk ) side
  real ( kind = rk ) t(dim_num,3)

  t(1:dim_num,1:3) = reshape ( (/ &
    p1(1:dim_num), p2(1:dim_num), p3(1:dim_num) /), (/ dim_num, 3 /) )

  call triangle_area_2d ( t, area_triangle )

  if ( area_triangle == 0.0D+00 ) then

    colin = 0.0D+00

  else

    perim = sqrt ( sum ( ( p2(1:dim_num) - p1(1:dim_num) )**2 ) ) &
          + sqrt ( sum ( ( p3(1:dim_num) - p2(1:dim_num) )**2 ) ) &
          + sqrt ( sum ( ( p1(1:dim_num) - p3(1:dim_num) )**2 ) )

    side = perim / 3.0D+00

    area2 = 0.25D+00 * sqrt ( 3.0D+00 ) * side * side

    colin = abs ( area_triangle ) / area2

  end if

  return
end
subroutine points_colin_3d ( p1, p2, p3, colin )

!*****************************************************************************80
!
!! POINTS_COLIN_3D estimates the colinearity of 3 points in 3D.
!
!  Discussion:
!
!    The estimate of colinearity that is returned is the ratio
!    of the area of the triangle spanned by the points to the area
!    of the equilateral triangle with the same perimeter.
!
!    This estimate is 1 if the points are maximally noncolinear, 0 if the
!    points are exactly colinear, and otherwise is closer to 1 or 0 depending
!    on whether the points are far or close to colinearity.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(3), P2(3), P3(3), the points.
!
!    Output, real ( kind = rk ) COLIN, the colinearity estimate. 
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) area_triangle
  real ( kind = rk ) area2
  real ( kind = rk ) colin
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) p3(dim_num)
  real ( kind = rk ) perim
  real ( kind = rk ) side
  real ( kind = rk ) t(dim_num,3)

  t(1:dim_num,1:3) = reshape ( (/ &
    p1(1:dim_num), p2(1:dim_num), p3(1:dim_num) /), (/ dim_num, 3 /) )

  call triangle_area_3d ( t, area_triangle )

  if ( area_triangle == 0.0D+00 ) then

    colin = 0.0D+00

  else

    perim = sqrt ( sum ( ( p2(1:dim_num) - p1(1:dim_num) )**2 ) ) &
          + sqrt ( sum ( ( p3(1:dim_num) - p2(1:dim_num) )**2 ) ) &
          + sqrt ( sum ( ( p1(1:dim_num) - p3(1:dim_num) )**2 ) )

    side = perim / 3.0D+00

    area2 = 0.25D+00 * sqrt ( 3.0D+00 ) * side * side

    colin = abs ( area_triangle ) / area2

  end if

  return
end
subroutine points_dist_nd ( dim_num, p1, p2, dist )

!*****************************************************************************80
!
!! POINTS_DIST_ND finds the distance between two points in ND.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = rk ) P1(DIM_NUM), P2(DIM_NUM), the coordinates 
!    of two points.
!
!    Output, real ( kind = rk ) DIST, the distance between the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num

  real ( kind = rk ) dist
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)

  dist = sqrt ( sum ( ( p1(1:dim_num) - p2(1:dim_num) )**2 ) )

  return
end
subroutine points_hull_2d ( node_num, node_xy, hull_num, hull )

!*****************************************************************************80
!
!! POINTS_HULL_2D computes the convex hull of 2D points.
!
!  Discussion:
!
!    The work involved is N*log(H), where N is the number of points, and H is
!    the number of points that are on the hull.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer NODE_NUM, the number of nodes.
!
!    real ( kind = rk ) NODE_XY(2,NODE_NUM), the coordinates of the nodes.
!
!  Output:
!
!    integer HULL_NUM, the number of nodes that lie on 
!    the convex hull.
!
!    integer HULL(NODE_NUM).  Entries 1 through HULL_NUM 
!    contain the indices of the nodes that form the convex hull, in order.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer node_num

  real ( kind = rk ) angle
  real ( kind = rk ) angle_max
  real ( kind = rk ) angle_rad_2d
  real ( kind = rk ) di
  real ( kind = rk ) dr
  integer first
  integer hull(node_num)
  integer hull_num
  integer i
  real ( kind = rk ) node_xy(2,node_num)
  real ( kind = rk ) p_xy(2)
  integer q
  real ( kind = rk ) q_xy(2)
  integer r
  real ( kind = rk ) r_xy(2)

  if ( node_num < 1 ) then
    hull_num = 0
    return
  end if
!
!  If NODE_NUM = 1, the hull is the point.
!
  if ( node_num == 1 ) then
    hull_num = 1
    hull(1) = 1
    return
  end if
!
!  If NODE_NUM = 2, then the convex hull is either the two distinct points,
!  or possibly a single (repeated) point.
!
  if ( node_num == 2 ) then

    if ( node_xy(1,1) /= node_xy(1,2) .or. node_xy(2,1) /= node_xy(2,2) ) then
      hull_num = 2
      hull(1) = 1
      hull(2) = 2
    else
      hull_num = 1
      hull(1) = 1
    end if

    return

  end if
!
!  Find the leftmost point and call it "Q".
!  In case of ties, take the bottom-most.
!
  q = 1
  do i = 2, node_num
    if ( node_xy(1,i) < node_xy(1,q) .or. &
       ( node_xy(1,i) == node_xy(1,q) .and. node_xy(2,i) < node_xy(2,q) ) ) then
      q = i
    end if
  end do

  q_xy(1:2) = node_xy(1:2,q)
!
!  Remember the starting point, so we know when to stop!
!
  first = q
  hull_num = 1
  hull(1) = q
!
!  For the first point, make a dummy previous point, 1 unit south,
!  and call it "P".
!
  p_xy(1) = q_xy(1)
  p_xy(2) = q_xy(2) - 1.0D+00
!
!  Now, having old point P, and current point Q, find the new point R
!  so the angle PQR is maximal.
!
!  Watch out for the possibility that the two nodes are identical.
!
  do

    r = 0
    angle_max = 0.0D+00

    do i = 1, node_num

      if ( i /= q .and. &
           ( node_xy(1,i) /= q_xy(1) .or. node_xy(2,i) /= q_xy(2) ) ) then

        angle = angle_rad_2d ( p_xy, q_xy, node_xy(1:2,i) )

        if ( r == 0 .or. angle_max < angle ) then

          r = i
          r_xy(1:2) = node_xy(1:2,r)
          angle_max = angle
!
!  In case of ties, choose the nearer point.
!
        else if ( r /= 0 .and. angle == angle_max ) then

          di = ( node_xy(1,i) - q_xy(1) )**2 + ( node_xy(2,i) - q_xy(2) )**2
          dr = ( r_xy(1)      - q_xy(1) )**2 + ( r_xy(2)      - q_xy(2) )**2

          if ( di < dr ) then
            r = i
            r_xy(1:2) = node_xy(1:2,r)
            angle_max = angle
          end if

        end if

      end if

    end do
!
!  We are done when we have returned to the first point on the convex hull.
!
    if ( r == first ) then
      exit
    end if

    hull_num = hull_num + 1

    if ( node_num < hull_num ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POINTS_HULL_2D - Fatal error!'
      write ( *, '(a)' ) '  The algorithm has failed.'
      stop 1
    end if
!
!  Add point R to convex hull.
!
    hull(hull_num) = r
!
!  Set P := Q, Q := R, and prepare to search for next point R.
!
    q = r

    p_xy(1:2) = q_xy(1:2)
    q_xy(1:2) = r_xy(1:2)

  end do

  return
end
subroutine points_plot ( file_name, node_num, node_xy, node_label )

!*****************************************************************************80
!
!! POINTS_PLOT plots a pointset.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the output file.
!
!    Input, integer NODE_NUM, the number of points.
!
!    Input, real ( kind = rk ) NODE_XY(2,NODE_NUM), the nodes.
!
!    Input, logical NODE_LABEL, is TRUE if the nodes should 
!    be labeled.
!
!  Local parameters:
!
!    Local, integer CIRCLE_SIZE, controls the size of the circles depicting
!    the nodes, measured in PostScript points (1/72 of an inch).  
!    Currently set to 5.  3 is pretty small, and 1 is barely visible.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )
  
  integer node_num

  integer, parameter :: circle_size = 5
  integer delta
  character ( len = * ) file_name
  integer file_unit
  integer ios
  integer node
  logical node_label
  real ( kind = rk ) node_xy(2,node_num)
  character ( len = 40 ) string
  real ( kind = rk ) x_max
  real ( kind = rk ) x_min
  integer x_ps
  integer :: x_ps_max = 576
  integer :: x_ps_max_clip = 594
  integer :: x_ps_min = 36
  integer :: x_ps_min_clip = 18
  real ( kind = rk ) x_scale
  real ( kind = rk ) y_max
  real ( kind = rk ) y_min
  integer y_ps
  integer :: y_ps_max = 666
  integer :: y_ps_max_clip = 684
  integer :: y_ps_min = 126
  integer :: y_ps_min_clip = 108
  real ( kind = rk ) y_scale
!
!  We need to do some figuring here, so that we can determine
!  the range of the data, and hence the height and width
!  of the piece of paper.
!
  x_max = maxval ( node_xy(1,1:node_num) )
  x_min = minval ( node_xy(1,1:node_num) )
  x_scale = x_max - x_min

  x_max = x_max + 0.05D+00 * x_scale
  x_min = x_min - 0.05D+00 * x_scale
  x_scale = x_max - x_min

  y_max = maxval ( node_xy(2,1:node_num) )
  y_min = minval ( node_xy(2,1:node_num) )
  y_scale = y_max - y_min

  y_max = y_max + 0.05D+00 * y_scale
  y_min = y_min - 0.05D+00 * y_scale
  y_scale = y_max - y_min

  if ( x_scale < y_scale ) then

    delta = nint ( real ( x_ps_max - x_ps_min, kind = rk ) &
      * ( y_scale - x_scale ) / ( 2.0D+00 * y_scale ) )

    x_ps_max = x_ps_max - delta
    x_ps_min = x_ps_min + delta

    x_ps_max_clip = x_ps_max_clip - delta
    x_ps_min_clip = x_ps_min_clip + delta

    x_scale = y_scale

  else if ( y_scale < x_scale ) then

    delta = nint ( real ( y_ps_max - y_ps_min, kind = rk ) &
      * ( x_scale - y_scale ) / ( 2.0D+00 * x_scale ) )

    y_ps_max = y_ps_max - delta
    y_ps_min = y_ps_min + delta

    y_ps_max_clip = y_ps_max_clip - delta
    y_ps_min_clip = y_ps_min_clip + delta

    y_scale = x_scale

  end if

  call get_unit ( file_unit )

  open ( unit = file_unit, file = file_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POINTS_PLOT - Fatal error!'
    write ( *, '(a)' ) '  Can not open output file.'
    stop 1
  end if

  write ( file_unit, '(a)' ) '%!PS-Adobe-3.0 EPSF-3.0'
  write ( file_unit, '(a)' ) '%%Creator: points_plot.f90'
  write ( file_unit, '(a)' ) '%%Title: ' // trim ( file_name )
  write ( file_unit, '(a)' ) '%%Pages: 1'
  write ( file_unit, '(a,i3,2x,i3,2x,i3,2x,i3)' ) '%%BoundingBox: ', &
    x_ps_min, y_ps_min, x_ps_max, y_ps_max
  write ( file_unit, '(a)' ) '%%Document-Fonts: Times-Roman'
  write ( file_unit, '(a)' ) '%%LanguageLevel: 1'
  write ( file_unit, '(a)' ) '%%EndComments'
  write ( file_unit, '(a)' ) '%%BeginProlog'
  write ( file_unit, '(a)' ) '/inch {72 mul} def'
  write ( file_unit, '(a)' ) '%%EndProlog'
  write ( file_unit, '(a)' ) '%%Page: 1 1'
  write ( file_unit, '(a)' ) 'save'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Set the RGB line color to very light gray.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '0.900  0.900  0.900 setrgbcolor'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Draw a gray border around the page.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) 'newpath'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_min, y_ps_min, ' moveto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_max, y_ps_min, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_max, y_ps_max, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_min, y_ps_max, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_min, y_ps_min, ' lineto'
  write ( file_unit, '(a)' ) 'stroke'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Set the RGB line color to black.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '0.000  0.000  0.000 setrgbcolor'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Set the font and its size.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '/Times-Roman findfont'
  write ( file_unit, '(a)' ) '0.50 inch scalefont'
  write ( file_unit, '(a)' ) 'setfont'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Print a title.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  210  702  moveto'
  write ( file_unit, '(a)' ) '%  (Pointset)  show'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Define a clipping polygon.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) 'newpath'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
    x_ps_min_clip, y_ps_min_clip, ' moveto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
    x_ps_max_clip, y_ps_min_clip, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
    x_ps_max_clip, y_ps_max_clip, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
    x_ps_min_clip, y_ps_max_clip, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
    x_ps_min_clip, y_ps_min_clip, ' lineto'
  write ( file_unit, '(a)' ) 'clip newpath'
!
!  Draw the nodes.
!
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Draw filled dots at each node.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Set the RGB color to blue.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '0.000  0.150  0.750 setrgbcolor'
  write ( file_unit, '(a)' ) '%'

  do node = 1, node_num

    x_ps = int ( &
      ( ( x_max - node_xy(1,node)         ) * real ( x_ps_min, kind = rk )   &
      + (         node_xy(1,node) - x_min ) * real ( x_ps_max, kind = rk ) ) &
      / ( x_max                   - x_min ) )

    y_ps = int ( &
      ( ( y_max - node_xy(2,node)         ) * real ( y_ps_min, kind = rk )   &
      + (         node_xy(2,node) - y_min ) * real ( y_ps_max, kind = rk ) ) &
      / ( y_max                   - y_min ) )

    write ( file_unit, '(a,i4,2x,i4,2x,i4,2x,a)' ) 'newpath ', x_ps, y_ps, &
      circle_size, '0 360 arc closepath fill'

  end do
!
!  Label the nodes.
!
  if ( node_label ) then

    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Label the nodes:'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Set the RGB color to darker blue.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '0.000  0.250  0.850 setrgbcolor'
    write ( file_unit, '(a)' ) '/Times-Roman findfont'
    write ( file_unit, '(a)' ) '0.20 inch scalefont'
    write ( file_unit, '(a)' ) 'setfont'

    do node = 1, node_num

      x_ps = int ( &
        ( ( x_max - node_xy(1,node)         ) * real ( x_ps_min, kind = rk )   &
        + (       + node_xy(1,node) - x_min ) * real ( x_ps_max, kind = rk ) ) &
        / ( x_max                   - x_min ) )

      y_ps = int ( &
        ( ( y_max - node_xy(2,node)         ) * real ( y_ps_min, kind = rk )   &
        + (         node_xy(2,node) - y_min ) * real ( y_ps_max, kind = rk ) ) &
        / ( y_max                   - y_min ) )

      write ( string, '(i4)' ) node
      string = adjustl ( string )

      write ( file_unit, '(i4,2x,i4,a)' ) x_ps, y_ps+5, &
        ' moveto (' // trim ( string ) // ') show'

    end do

  end if

  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) 'restore  showpage'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  End of page.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%%Trailer'
  write ( file_unit, '(a)' ) '%%EOF'
  close ( unit = file_unit )

  return
end
subroutine points_point_near_naive_nd ( dim_num, set_num, pset, p, i_min, &
  dist_min )

!*****************************************************************************80
!
!! POINTS_POINT_NEAR_NAIVE_ND finds the nearest point to a given point in ND.
!
!  Discussion:
!
!    A naive algorithm is used.  The distance to every point is calculated,
!    in order to determine the smallest.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer SET_NUM, the number of points in the set.
!
!    Input, real ( kind = rk ) PSET(DIM_NUM,SET_NUM), the points in the set.
!
!    Input, real ( kind = rk ) P(DIM_NUM), the point whose nearest neighbor
!    is sought.
!
!    Output, integer I_MIN, the index of the nearest point in 
!    PSET to P.
!
!    Output, real ( kind = rk ) DIST_MIN, the distance between P(*) 
!    and PSET(*,I_MIN).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num
  integer set_num

  real ( kind = rk ) d
  real ( kind = rk ) dist_min
  integer i
  integer i_min
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) pset(dim_num,set_num)

  dist_min = huge ( dist_min )
  i_min = -1

  do i = 1, set_num
    d = sum ( ( p(1:dim_num) - pset(1:dim_num,i) )**2 )
    if ( d < dist_min ) then
      dist_min = d
      i_min = i
    end if
  end do

  dist_min = sqrt ( dist_min )

  return
end
subroutine polar_to_xy ( r, t, xy )

!*****************************************************************************80
!
!! POLAR_TO_XY converts polar coordinates to XY coordinates.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) R, T, the radius and angle (in radians).
!
!    Output, real ( kind = rk ) XY(2), the Cartesian coordinates.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r
  real ( kind = rk ) t
  real ( kind = rk ) xy(2)

  xy(1) = r * cos ( t )
  xy(2) = r * sin ( t )

  return
end
subroutine polygon_area_2d ( n, v, area )

!*****************************************************************************80
!
!! polygon_area_2d() computes the area of a polygon in 2D.
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
  integer, parameter :: dim_num = 2

  real ( kind = rk ) area
  integer i
  integer i4_wrap
  integer im1
  integer ip1
  real ( kind = rk ) v(dim_num,n)

  area = 0.0D+00

  do i = 1, n

    im1 = i4_wrap ( i-1, 1, n )
    ip1 = i4_wrap ( i+1, 1, n )

    area = area + v(1,i) * ( v(2,ip1) - v(2,im1) )

  end do

  area = 0.5D+00 * area

  return
end
subroutine polygon_area_3d ( n, v, area, normal )

!*****************************************************************************80
!
!! polygon_area_3d() computes the area of a polygon in 3D.
!
!  Discussion:
!
!    The computation is not valid unless the vertices of the polygon
!    lie in a plane, so that the polygon that is defined is "flat".
!
!    The polygon does not have to be "regular", that is, neither its
!    sides nor its angles need to be equal.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Allen Van Gelder,
!    Efficient Computation of Polygon Area and Polyhedron Volume,
!    Graphics Gems V, 
!    edited by Alan Paeth,
!    AP Professional, 1995, T385.G6975.
!
!  Parameters:
!
!    Input, integer N, the number of vertices.
!
!    Input, real ( kind = rk ) V(3,N), the coordinates of the vertices.
!    The vertices should be listed in neighboring order.
!
!    Output, real ( kind = rk ) AREA, the area of the polygon.
!
!    Output, real ( kind = rk ) NORMAL(3), the unit normal vector to the polygon.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer, parameter :: dim_num = 3

  real ( kind = rk ) area
  real ( kind = rk ) cross(dim_num)
  integer i
  integer ip1
  real ( kind = rk ) normal(dim_num)
  real ( kind = rk ) v(dim_num,n)

  normal(1:dim_num) = 0.0D+00

  do i = 1, n

    if ( i < n ) then
      ip1 = i + 1
    else
      ip1 = 1
    end if
!
!  Compute the cross product vector.
!
    cross(1) = v(2,i) * v(3,ip1) - v(3,i) * v(2,ip1)
    cross(2) = v(3,i) * v(1,ip1) - v(1,i) * v(3,ip1)
    cross(3) = v(1,i) * v(2,ip1) - v(2,i) * v(1,ip1)

    normal(1:dim_num) = normal(1:dim_num) + cross(1:dim_num)

  end do

  area = sqrt ( sum ( normal(1:dim_num)**2 ) )

  if ( area /= 0.0D+00 ) then
    normal(1:dim_num) = normal(1:dim_num) / area
  else
    normal(1:dim_num) = 1.0D+00 / sqrt ( real ( dim_num, kind = rk ) )
  end if

  area = 0.5D+00 * area

  return
end
subroutine polygon_centroid_3d ( n, v, centroid )

!*****************************************************************************80
!
!! polygon_centroid_3d() computes the centroid of a polygon in 3D.
!
!  Discussion:
!
!    The polygon is described by its vertices.  In many applications,
!    these vertices will lie in a common plane, and the polygon will
!    be "flat".  However, that is not required for this formula.
!
!    This formula triangulates the polygon, computes the area of
!    each triangle and its centroid, and then computes the centroid
!    of the polygon as the weight-averaged sum of the triangle centroids.
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
!    Adrian Bowyer, John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983,
!    ISBN: 0408012420.
!
!  Parameters:
!
!    Input, integer N, the number of vertices of the polygon.
!
!    Input, real ( kind = rk ) V(3,N), the coordinates of the vertices.
!
!    Output, real ( kind = rk ) CENTROID(3), the coordinates of the centroid.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: dim_num = 3

  real ( kind = rk ) area_polygon
  real ( kind = rk ) area_triangle
  real ( kind = rk ) centroid(dim_num)
  integer i
  real ( kind = rk ) t(dim_num,3)
  real ( kind = rk ) v(dim_num,n)

  area_polygon = 0.0D+00
  centroid(1:dim_num) = 0.0D+00

  do i = 1, n - 2

    t(1:dim_num,1:3) = reshape ( (/ &
      v(1:dim_num,i), v(1:dim_num,i+1), v(1:dim_num,n) /), &
      (/ dim_num, 3 /) )

    call triangle_area_3d ( t, area_triangle )

    area_polygon = area_polygon + area_triangle

    centroid(1:dim_num) = centroid(1:dim_num) + area_triangle &
      * ( v(1:dim_num,i) + v(1:dim_num,i+1) + v(1:dim_num,n) ) / 3.0D+00

  end do

  if ( area_polygon == 0.0D+00 ) then
    centroid(1:dim_num) = v(1:dim_num,1)
  else
    centroid(1:dim_num) = centroid(1:dim_num) / area_polygon
  end if

  return
end
subroutine polygon_normal_3d ( n, v, normal ) 

!*****************************************************************************80
!
!! polygon_normal_3d() computes the normal vector to a polygon in 3D.
!
!  Discussion:
!
!    If the polygon is planar, then this calculation is correct.
!
!    Otherwise, the normal vector calculated is the simple average
!    of the normals defined by the planes of successive triples
!    of vertices.
!
!    If the polygon is "almost" planar, this is still acceptable.
!    But as the polygon is less and less planar, so this averaged normal
!    vector becomes more and more meaningless.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paulo Cezar Pinto Carvalho, Paulo Roma Cavalcanti,
!    Point in Polyhedron Testing Using Spherical Polygons,
!    in Graphics Gems V,
!    edited by Alan Paeth,
!    Academic Press, 1995,
!    ISBN: 0125434553,
!    LC: T385.G6975.
!
!  Parameters:
!
!    Input, integer N, the number of vertices.
!
!    Input, real ( kind = rk ) V(3,N), the coordinates of the vertices.
!
!    Output, real ( kind = rk ) NORMAL(3), the averaged normal vector
!    to the polygon. 
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer n

  real ( kind = rk ) r8vec_norm
  integer j
  real ( kind = rk ) normal(dim_num)
  real ( kind = rk ) normal_norm
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) v(dim_num,n)
  real ( kind = rk ) v1(dim_num)
  real ( kind = rk ) v2(dim_num)

  normal(1:dim_num) = 0.0D+00

  v1(1:dim_num) = v(1:dim_num,2) - v(1:dim_num,1)

  do j = 3, n

    v2(1:dim_num) = v(1:dim_num,j) - v(1:dim_num,1)

    call r8vec_cross_product_3d ( v1, v2, p )

    normal(1:dim_num) = normal(1:dim_num) + p(1:dim_num)

    v1(1:dim_num) = v2(1:dim_num)

  end do
!
!  Normalize.
!
  normal_norm = r8vec_norm ( dim_num, normal )

  if ( normal_norm == 0.0D+00 ) then
    return
  end if

  normal(1:dim_num) = normal(1:dim_num) / normal_norm

  return
end
subroutine polygon_solid_angle_3d ( n, v, p, solid_angle )

!*****************************************************************************80
!
!! POLYGON_SOLID_ANGLE_3D: projected solid angle of a 3D plane polygon.
!
!  Discussion:
!
!    A point P is at the center of a unit sphere.  A planar polygon
!    is to be projected onto the surface of this sphere, by drawing
!    the ray from P to each polygonal vertex, and noting where this ray
!    intersects the sphere.
!
!    We compute the area on the sphere of the projected polygon.
!
!    Since we are projecting the polygon onto a unit sphere, the area 
!    of the projected polygon is equal to the solid angle subtended by 
!    the polygon.
!
!    The value returned by this routine will include a sign.  The
!    angle subtended will be NEGATIVE if the normal vector defined by
!    the polygon points AWAY from the viewing point, and will be
!    POSITIVE if the normal vector points towards the viewing point.
!
!    If the orientation of the polygon is of no interest to you,
!    then you can probably simply take the absolute value of the
!    solid angle as the information you want.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    29 October 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paulo Cezar Pinto Carvalho, Paulo Roma Cavalcanti,
!    Point in Polyhedron Testing Using Spherical Polygons,
!    in Graphics Gems V,
!    edited by Alan Paeth,
!    Academic Press, 1995,
!    ISBN: 0125434553,
!    LC: T385.G6975.
!
!  Parameters:
!
!    Input, integer N, the number of vertices.
!
!    Input, real ( kind = rk ) V(3,N), the coordinates of the vertices.
!
!    Input, real ( kind = rk ) P(3), the point at the center of the unit sphere.
!
!    Output, double SOLID_ANGLE, the solid angle subtended
!    by the polygon, as projected onto the unit sphere around the point P.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer n

  real ( kind = rk ) a(dim_num)
  real ( kind = rk ) angle
  real ( kind = rk ) area
  real ( kind = rk ) b(dim_num)
  real ( kind = rk ) r8vec_norm
  real ( kind = rk ) r8vec_scalar_triple_product
  integer i4_wrap
  integer j
  integer jp1
  real ( kind = rk ) normal1(dim_num)
  real ( kind = rk ) normal1_norm
  real ( kind = rk ) normal2(dim_num)
  real ( kind = rk ) normal2_norm
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) plane(dim_num)
  real ( kind = rk ) r1(dim_num)
  real ( kind = rk ) r8_acos
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) s
  real ( kind = rk ) solid_angle
  real ( kind = rk ) v(dim_num,n)

  if ( n < 3 ) then
    solid_angle = 0.0D+00
    return
  end if

  call polygon_normal_3d ( n, v, plane )
 
  a(1:dim_num) = v(1:dim_num,n) - v(1:dim_num,1)

  area = 0.0D+00

  do j = 1, n

    r1(1:dim_num) = v(1:dim_num,j) - p(1:dim_num)

    jp1 = i4_wrap ( j + 1, 1, n )

    b(1:dim_num) = v(1:dim_num,jp1) - v(1:dim_num,j)

    call r8vec_cross_product_3d ( a, r1, normal1 )

    normal1_norm = r8vec_norm ( dim_num, normal1 )

    call r8vec_cross_product_3d ( r1, b, normal2 )

    normal2_norm = r8vec_norm ( dim_num, normal2 )
    
    s = dot_product ( normal1(1:dim_num), normal2(1:dim_num) ) &
      / ( normal1_norm * normal2_norm )

    angle = r8_acos ( s )

    s = r8vec_scalar_triple_product ( b, a, plane )

    if ( 0.0D+00 < s ) then
      area = area + r8_pi - angle
    else
      area = area + r8_pi + angle
    end if

    a(1:dim_num) = -b(1:dim_num)

  end do

  area = area - r8_pi * real ( n - 2, kind = rk )

  if ( 0.0D+00 < dot_product ( plane(1:dim_num), r1(1:dim_num) ) ) then
    solid_angle = -area
  else
    solid_angle = area
  end if

  return
end
subroutine polyhedron_area_3d ( coord, order_max, face_num, node, &
  node_num, order, area )

!*****************************************************************************80
!
!! POLYHEDRON_AREA_3D computes the surface area of a polyhedron in 3D.
!
!  Discussion:
!
!    The computation is not valid unless the faces of the polyhedron
!    are planar polygons.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Allen Van Gelder,
!    Efficient Computation of Polygon Area and Polyhedron Volume,
!    in Graphics Gems V, 
!    edited by Alan Paeth,
!    AP Professional, 1995, T385.G6975
!
!  Parameters:
!
!    Input, real ( kind = rk ) COORD(3,NODE_NUM), the coordinates of the
!    vertices.  The vertices may be listed in any order.
!
!    Input, integer ORDER_MAX, the maximum number of vertices 
!    that make up a face of the polyhedron.
!
!    Input, integer FACE_NUM, the number of faces of the 
!    polyhedron.
!
!    Input, integer NODE(FACE_NUM,ORDER_MAX).  Face I is defined 
!    by the vertices NODE(I,1) through NODE(I,ORDER(I)).  These vertices
!    are listed in neighboring order.
!
!    Input, integer NODE_NUM, the number of points stored in COORD.
!
!    Input, integer ORDER(FACE_NUM), the number of vertices 
!    making up each face.
!
!    Output, real ( kind = rk ) AREA, the total surface area of the polyhedron.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer face_num
  integer order_max
  integer, parameter :: dim_num = 3
  integer node_num

  real ( kind = rk ) ainc
  real ( kind = rk ) area
  real ( kind = rk ) coord(dim_num,node_num)
  integer face
  integer j
  integer k1
  integer k2
  integer node(face_num,order_max)
  integer order(face_num)
  real ( kind = rk ) v(dim_num)

  area = 0.0D+00
!
!  For each face
!
  do face = 1, face_num

    v(1:dim_num) = 0.0D+00
!
!  For each triangle in the face, compute the normal vector.
!
    do j = 1, order(face)

      k1 = node(face,j)

      if ( j < order(face) ) then
        k2 = node(face,j+1)
      else
        k2 = node(face,1)
      end if
!
!  Compute the cross product.
!
      v(1) = v(1) + coord(2,k1) * coord(3,k2) - coord(3,k1) * coord(2,k2)
      v(2) = v(2) + coord(3,k1) * coord(1,k2) - coord(1,k1) * coord(3,k2)
      v(3) = v(3) + coord(1,k1) * coord(2,k2) - coord(2,k1) * coord(1,k2)

    end do
!
!  Add the magnitude of the normal vector to the sum.
!
    ainc = sqrt ( sum ( v(1:dim_num)**2 ) )
    area = area + ainc

  end do

  area = 0.5D+00 * area

  return
end
subroutine polyhedron_centroid_3d ( coord, order_max, face_num, node, &
  node_num, order, centroid )

!*****************************************************************************80
!
!! POLYHEDRON_CENTROID_3D computes the centroid of a polyhedron in 3D.
!
!  Discussion:
!
!    The centroid can be computed as the volume-weighted average of
!    the centroids of the tetrahedra defined by choosing a point in
!    the interior of the polyhedron, and using as a base every triangle
!    created by triangulating the faces of the polyhedron.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) COORD(3,NODE_NUM), the vertices.
!    The vertices may be listed in any order.
!
!    Input, integer ORDER_MAX, the maximum number of vertices 
!    that make up a face of the polyhedron.
!
!    Input, integer FACE_NUM, the number of faces of the 
!    polyhedron.
!
!    Input, integer NODE(FACE_NUM,ORDER_MAX).  Face I is defined 
!    by the vertices NODE(I,1) through NODE(I,ORDER(I)).  These vertices
!    are listed in neighboring order.
!
!    Input, integer NODE_NUM, the number of points stored in COORD.
!
!    Input, integer ORDER(FACE_NUM), the number of vertices making 
!    up each face.
!
!    Output, real ( kind = rk ) CENTROID(3), the centroid of the polyhedron.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer face_num
  integer order_max
  integer, parameter :: dim_num = 3
  integer node_num

  real ( kind = rk ) area
  real ( kind = rk ) centroid(dim_num)
  real ( kind = rk ) coord(dim_num,node_num)
  integer face
  integer n1
  integer n2
  integer n3
  integer node(face_num,order_max)
  real ( kind = rk ) normal(dim_num)
  integer order(face_num)
  real ( kind = rk ) point(dim_num)
  real ( kind = rk ) polygon_area
  real ( kind = rk ) polygon_centroid(dim_num)
  real ( kind = rk ) tetra(dim_num,4)
  real ( kind = rk ) tetra_centroid(dim_num)
  real ( kind = rk ) tetra_volume
  integer vert
  integer vert_num
  real ( kind = rk ) volume
  real ( kind = rk ) v(dim_num,order_max)
!
!  Compute a point in the interior.
!  We take the area-weighted centroid of each face.
!
  point(1:dim_num) = 0.0D+00
  area = 0.0D+00

  do face = 1, face_num

    vert_num = order(face)

    v(1:dim_num,1:vert_num) = coord(1:dim_num,node(face,1:vert_num))

    call polygon_area_3d ( vert_num, v, polygon_area, normal )

    call polygon_centroid_3d ( vert_num, v, polygon_centroid )

    point(1:dim_num) = point(1:dim_num) &
      + polygon_area * polygon_centroid(1:dim_num)

    area = area + polygon_area

  end do

  point(1:dim_num) = point(1:dim_num) / area
!
!  Now triangulate each face.
!  For each triangle, consider the tetrahedron created by including POINT.
!
  centroid(1:dim_num) = 0.0D+00
  volume = 0.0D+00

  do face = 1, face_num

    n3 = node(face,order(face))

    do vert = 1, order(face) - 2

      n1 = node(face,vert)
      n2 = node(face,vert+1)

      tetra(1:dim_num,1:4) = reshape ( (/ &
        coord(1:dim_num,n1), coord(1:dim_num,n2), coord(1:dim_num,n3), &
        point(1:dim_num) /), (/ dim_num, 4 /) )

      call tetrahedron_volume_3d ( tetra, tetra_volume )

      call tetrahedron_centroid_3d ( tetra, tetra_centroid )

      centroid(1:dim_num) = centroid(1:dim_num) &
        + tetra_volume * tetra_centroid(1:dim_num)

      volume = volume + tetra_volume

    end do
  end do

  centroid(1:dim_num) = centroid(1:dim_num) / volume

  return
end
subroutine polyhedron_contains_point_3d ( node_num, face_num, &
  face_order_max, v, face_order, face_point, p, inside )

!*****************************************************************************80
!
!! POLYHEDRON_CONTAINS_POINT_3D determines if a point is inside a polyhedron.
!
!  Discussion:
!
!    The reference states that the polyhedron should be simple (that
!    is, the faces should form a single connected surface), and that 
!    the individual faces should be consistently oriented.
!
!    However, the polyhedron does not, apparently, need to be convex.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paulo Cezar Pinto Carvalho, Paulo Roma Cavalcanti,
!    Point in Polyhedron Testing Using Spherical Polygons,
!    in Graphics Gems V,
!    edited by Alan Paeth,
!    Academic Press, 1995,
!    ISBN: 0125434553,
!    LC: T385.G6975.
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of vertices.
!
!    Input, integer FACE_NUM, the number of faces.
!
!    Input, integer FACE_ORDER_MAX, the maximum order of any face.
!
!    Input, real ( kind = rk ) V(3,NODE_NUM), the coordinates of the vertices.
!
!    Input, integer FACE_ORDER(FACE_NUM), the order of each face.
!
!    Input, integer FACE_POINT(FACE_ORDER_MAX,FACE_NUM), the 
!    indices of the nodes that make up each face.
!
!    Input, real ( kind = rk ) P(3), the point to be tested.
!
!    Output, logical INSIDE, is true if the point 
!    is inside the polyhedron.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer face_num
  integer face_order_max
  integer node_num

  real ( kind = rk ) area
  integer face
  integer face_order(face_num)
  integer face_point(face_order_max,face_num)
  logical inside
  integer k
  integer node
  integer node_num_face
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) solid_angle
  real ( kind = rk ) v(dim_num,node_num)
  real ( kind = rk ) v_face(dim_num,face_order_max)

  area = 0.0D+00

  do face = 1, face_num

    node_num_face = face_order(face)

    do k = 1, node_num_face

      node = face_point(k,face)

      v_face(1:dim_num,k) = v(1:dim_num,node)

    end do

    call polygon_solid_angle_3d ( node_num_face, v_face, p, solid_angle )

    area = area + solid_angle

  end do
!
!  AREA should be -4*PI, 0, or 4*PI.
!  So this test should be quite safe!
!
  if ( area < -2.0D+00 * r8_pi .or. 2.0D+00 * r8_pi < area ) then
    inside = .true.
  else
    inside = .false.
  end if

  return
end
subroutine polyhedron_volume_3d ( coord, order_max, face_num, node, &
  node_num, order, volume )

!*****************************************************************************80
!
!! POLYHEDRON_VOLUME_3D computes the volume of a polyhedron in 3D.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) COORD(3,NODE_NUM), the coordinates of 
!    the vertices.  The vertices may be listed in any order.
!
!    Input, integer ORDER_MAX, the maximum number of vertices 
!    that make up a face of the polyhedron.
!
!    Input, integer FACE_NUM, the number of faces of the 
!    polyhedron.
!
!    Input, integer NODE(FACE_NUM,ORDER_MAX).  Face I is defined by
!    the vertices NODE(I,1) through NODE(I,ORDER(I)).  These vertices
!    are listed in neighboring order.
!
!    Input, integer NODE_NUM, the number of points stored in COORD.
!
!    Input, integer ORDER(FACE_NUM), the number of vertices making 
!    up each face.
!
!    Output, real ( kind = rk ) VOLUME, the volume of the polyhedron.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer face_num
  integer order_max
  integer, parameter :: dim_num = 3
  integer node_num

  real ( kind = rk ) coord(dim_num,node_num)
  integer face
  integer n1
  integer n2
  integer n3
  integer node(face_num,order_max)
  integer order(face_num)
  integer v
  real ( kind = rk ) volume

  volume = 0.0D+00
!
!  Triangulate each face.
!
  do face = 1, face_num

    n3 = node(face,order(face))

    do v = 1, order(face) - 2

      n1 = node(face,v)
      n2 = node(face,v+1)

      volume = volume &
        + coord(1,n1) &
        * ( coord(2,n2) * coord(3,n3) - coord(2,n3) * coord(3,n2) ) &
        + coord(1,n2) &
        * ( coord(2,n3) * coord(3,n1) - coord(2,n1) * coord(3,n3) ) &
        + coord(1,n3) &
        * ( coord(2,n1) * coord(3,n2) - coord(2,n2) * coord(3,n1) )

    end do

  end do

  volume = volume / 6.0D+00

  return
end
subroutine polyhedron_volume_3d_2 ( coord, order_max, face_num, node, &
  node_num, order, volume )

!*****************************************************************************80
!
!! POLYHEDRON_VOLUME_3D_2 computes the volume of a polyhedron in 3D.
!
!  Discussion:
!
!    The computation is not valid unless the faces of the polyhedron
!    are planar polygons.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Allen Van Gelder,
!    Efficient Computation of Polygon Area and Polyhedron Volume,
!    in Graphics Gems V,
!    edited by Alan Paeth,
!    AP Professional, 1995, T385.G6975.
!
!  Parameters:
!
!    Input, real ( kind = rk ) COORD(3,NODE_NUM), the vertices.
!    The vertices may be listed in any order.
!
!    Input, integer ORDER_MAX, the maximum number of vertices 
!    that make up a face of the polyhedron.
!
!    Input, integer FACE_NUM, the number of faces of the 
!    polyhedron.
!
!    Input, integer NODE(FACE_NUM,ORDER_MAX).  Face I is defined 
!    by the vertices NODE(I,1) through NODE(I,ORDER(I)).  These vertices
!    are listed in neighboring order.
!
!    Input, integer NODE_NUM, the number of points stored in COORD.
!
!    Input, integer ORDER(FACE_NUM), the number of vertices making 
!    up each face.
!
!    Output, real ( kind = rk ) VOLUME, the volume of the polyhedron.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer face_num
  integer order_max
  integer, parameter :: dim_num = 3
  integer node_num

  real ( kind = rk ) coord(dim_num,node_num)
  integer face
  integer j
  integer k
  integer k1
  integer k2
  integer node(face_num,order_max)
  real ( kind = rk ) normal(dim_num)
  integer order(face_num)
  real ( kind = rk ) v(dim_num)
  real ( kind = rk ) volume

  volume = 0.0D+00

  do face = 1, face_num

    v(1:dim_num) = 0.0D+00
!
!  Compute the area vector for this face.
!
    do j = 1, order(face)

      k1 = node(face,j)

      if ( j < order(face) ) then
        k2 = node(face,j+1)
      else
        k2 = node(face,1)
      end if
!
!  Compute the cross product.
!
      normal(1) = coord(2,k1) * coord(3,k2) - coord(3,k1) * coord(2,k2)
      normal(2) = coord(3,k1) * coord(1,k2) - coord(1,k1) * coord(3,k2)
      normal(3) = coord(1,k1) * coord(2,k2) - coord(2,k1) * coord(1,k2)

      v(1:dim_num) = v(1:dim_num) + normal(1:dim_num)

    end do
!
!  Area vector dot any vertex.
!
    k = node(face,1)
    volume = volume + dot_product ( v(1:dim_num), coord(1:dim_num,k) )

  end do

  volume = volume / 6.0D+00

  return
end
subroutine polyline_arclength_nd ( dim_num, n, p, s )

!*****************************************************************************80
!
!! POLYLINE_ARCLENGTH_ND computes the arclength of points on a polyline in ND.
!
!  Discussion:
!
!    A polyline of order N is the geometric structure consisting of
!    the N-1 line segments that lie between successive elements of a list
!    of N points.
!
!    An ordinary line segment is a polyline of order 2.
!    The letter "V" is a polyline of order 3.
!    The letter "N" is a polyline of order 4, and so on.
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
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer N, the number of points defining the polyline.
!
!    Input, real ( kind = rk ) P(DIM_NUM,N), the points defining the polyline.
!
!    Output, real ( kind = rk ) S(N), the arclength coordinates
!    of each point.  The first point has S(1) = 0 and the 
!    last point has S(N) = arclength of the entire polyline.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num
  integer n

  integer i
  real ( kind = rk ) p(dim_num,n)
  real ( kind = rk ) s(n)

  s(1) = 0.0D+00

  do i = 2, n

    s(i) = s(i-1) + sqrt ( sum ( ( p(1:dim_num,i) - p(1:dim_num,i-1) )**2 ) )

  end do

  return
end
subroutine polyline_index_point_nd ( dim_num, n, p, t, pt )

!*****************************************************************************80
!
!! POLYLINE_INDEX_POINT_ND evaluates a polyline at a given arclength in ND.
!
!  Discussion:
!
!    The polyline is defined as the set of N-1 line segments lying
!    between a sequence of N points.  The arclength of a point lying
!    on the polyline is simply the length of the broken line from the
!    initial point.  Any point on the polyline can be found by
!    specifying its arclength.
!
!    If the given arclength coordinate is less than 0, or greater
!    than the arclength coordinate of the last given point, then
!    extrapolation is used, that is, the first and last line segments
!    are extended as necessary.
!
!    The arclength coordinate system measures the distance between
!    any two points on the polyline as the length of the segment of the
!    line that joins them.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer N, the number of points defining the polyline.
!
!    Input, real ( kind = rk ) P(DIM_NUM,N), the points defining the polyline.
!
!    Input, real ( kind = rk ) T, the desired arclength coordinate.
!
!    Output, real ( kind = rk ) PT(DIM_NUM), the point corresponding to the
!    arclength.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer dim_num

  integer i
  real ( kind = rk ) p(dim_num,n)
  real ( kind = rk ) pt(dim_num)
  real ( kind = rk ) t
  real ( kind = rk ) t1
  real ( kind = rk ) t2

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POLYLINE_INDEX_POINT_ND - Fatal error!'
    write ( *, '(a)' ) '  The input quantity N is nonpositive.'
    write ( *, '(a,i8)' ) '  N = ', n
    stop 1
  end if

  if ( n == 1 ) then

    pt(1:dim_num) = p(1:dim_num,1)

  else

    t2 = 0.0D+00

    do i = 1, n - 1
!
!  Find the distance between points I and I+1.
!
      t1 = t2
      t2 = t1 + sqrt ( sum ( ( p(1:dim_num,i+1) - p(1:dim_num,i) )**2 ) )
!
!  Interpolate or extrapolate in an interval.
!
      if ( t <= t2 .or. i == n - 1 ) then

        pt(1:dim_num) = ( ( t2 - t      ) * p(1:dim_num,i)     &
                        + (      t - t1 ) * p(1:dim_num,i+1) ) &
                        / ( t2     - t1 )

        return
      end if
    end do
  end if

  return
end
subroutine polyline_length_nd ( dim_num, nk, pk, length )

!*****************************************************************************80
!
!! POLYLINE_LENGTH_ND computes the length of a polyline in ND.
!
!  Discussion:
!
!    A polyline of order NK is the geometric structure consisting of
!    the NK-1 line segments that lie between successive elements of a list
!    of NK points.
!
!    An ordinary line segment is a polyline of order 2.
!    The letter "V" is a polyline of order 3.
!    The letter "N" is a polyline of order 4, and so on.
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
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer NK, the number of points defining the polyline.
!
!    Input, real ( kind = rk ) PK(DIM_NUM,NK), the points defining the polyline.
!
!    Output, real ( kind = rk ) LENGTH, the length of the polyline.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num
  integer nk

  integer i
  real ( kind = rk ) length
  real ( kind = rk ) pk(dim_num,nk)

  length = 0.0D+00

  do i = 2, nk

    length = length &
      + sqrt ( sum ( ( pk(1:dim_num,i) - pk(1:dim_num,i-1) )**2 ) )

  end do

  return
end
subroutine polyline_points_nd ( dim_num, n, p, nt, pt )

!*****************************************************************************80
!
!! POLYLINE_POINTS_ND computes equally spaced points on a polyline in ND.
!
!  Discussion:
!
!    A polyline of order N is the geometric structure consisting of
!    the N-1 line segments that lie between successive elements of a list
!    of N points.
!
!    An ordinary line segment is a polyline of order 2.
!    The letter "V" is a polyline of order 3.
!    The letter "N" is a polyline of order 4, and so on.
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
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer N, the number of points defining the polyline.
!
!    Input, real ( kind = rk ) P(DIM_NUM,N), the points defining the polyline.
!
!    Input, integer NT, the number of points to be sampled.
!
!    Output, real ( kind = rk ) PT(DIM_NUM,NT), equally spaced points
!    on the polyline.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num
  integer n
  integer nt

  integer it
  integer j
  real ( kind = rk ) p(dim_num,n)
  real ( kind = rk ) pt(dim_num,nt)
  real ( kind = rk ) s(n)
  real ( kind = rk ) st

  call polyline_arclength_nd ( dim_num, n, p, s )

  j = 1

  do it = 1,  nt

    st = ( real ( nt - it,     kind = rk ) * 0.0D+00 + &
           real (      it - 1, kind = rk ) * s(n) ) &
         / real ( nt      - 1, kind = rk )

    do

      if ( s(j) <= st .and. st <= s(j+1) ) then
        exit
      end if

      if ( n - 1 <= j ) then
        exit
      end if

      j = j + 1

    end do

    pt(1:dim_num,it) = ( ( s(j+1) - st        ) * p(1:dim_num,j) &
                       + (          st - s(j) ) * p(1:dim_num,j+1) ) &
                       / ( s(j+1)      - s(j) )

  end do

  return
end
subroutine polyloop_arclength_nd ( dim_num, nk, pk, sk )

!*****************************************************************************80
!
!! POLYLOOP_ARCLENGTH_ND computes the arclength of points on a polyloop in ND.
!
!  Discussion:
!
!    A polyloop of order NK is the geometric structure consisting of
!    the NK line segments that lie between successive elements of a list
!    of NK points, with the last point joined to the first.
!
!    Warning: I just made up the word "polyloop", so don't go repeating it!
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
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer NK, the number of points defining the polyloop.
!
!    Input, real ( kind = rk ) PK(DIM_NUM,NK), the points defining the polyloop.
!
!    Output, real ( kind = rk ) SK(NK+1), the arclength coordinates
!    of each point.  The first point has two arc length values,
!    namely SK(1) = 0 and SK(NK+1) = LENGTH.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num
  integer nk

  integer i
  integer j
  real ( kind = rk ) pk(dim_num,nk)
  real ( kind = rk ) sk(nk+1)

  sk(1) = 0.0D+00

  do i = 2, nk + 1

    if ( i <= nk ) then
      j = i
    else
      j = 1
    end if

    sk(i) = sk(i-1) &
      + sqrt ( sum ( ( pk(1:dim_num,j) - pk(1:dim_num,i-1) )**2 ) )

  end do

  return
end
subroutine polyloop_length_nd ( dim_num, nk, pk, length )

!*****************************************************************************80
!
!! POLYLOOP_LENGTH_ND computes the length of a polyloop in ND.
!
!  Discussion:
!
!    A polyloop of order NK is the geometric structure consisting of
!    the NK line segments that lie between successive elements of a list
!    of NK points, with the last point joined to the first.
!
!    Warning: I just made up the word "polyloop", so don't go repeating it!
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
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer NK, the number of points defining the polyloop.
!
!    Input, real ( kind = rk ) PK(DIM_NUM,NK), the points defining the polyloop.
!
!    Output, real ( kind = rk ) LENGTH, the length of the polyloop.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num
  integer nk

  integer i
  integer j
  real ( kind = rk ) length
  real ( kind = rk ) pk(dim_num,nk)

  length = 0.0D+00

  do i = 2, nk + 1

    if ( i <= nk ) then
      j = i
    else
      j = 1
    end if

    length = length &
      + sqrt ( sum ( ( pk(1:dim_num,j) - pk(1:dim_num,i-1) )**2 ) )

  end do

  return
end
subroutine polyloop_points_nd ( dim_num, nk, pk, nt, pt )

!*****************************************************************************80
!
!! POLYLOOP_POINTS_ND computes equally spaced points on a polyloop in ND.
!
!  Discussion:
!
!    A polyloop of order NK is the geometric structure consisting of
!    the NK line segments that lie between successive elements of a list
!    of NK points, including a segment from the last point to the first.
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
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer NK, the number of points defining the polyloop.
!
!    Input, real ( kind = rk ) PK(DIM_NUM,NK), the points defining the polyloop.
!
!    Input, integer NT, the number of points to be sampled.
!
!    Input, real ( kind = rk ) PT(DIM_NUM,NT), equally spaced points
!    on the polyloop.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num
  integer nk
  integer nt

  integer it
  integer i4_wrap
  integer j
  integer jp1
  real ( kind = rk ) pk(dim_num,nk)
  real ( kind = rk ) pt(dim_num,nt)
  real ( kind = rk ) sk(nk+1)
  real ( kind = rk ) st

  call polyloop_arclength_nd ( dim_num, nk, pk, sk )

  j = 1

  do it = 1,  nt

    st = ( real ( nt - it,     kind = rk ) * 0.0D+00 + &
           real (      it - 1, kind = rk ) * sk(nk+1) ) &
         / real ( nt      - 1, kind = rk )

    do

      if ( sk(j) <= st .and. st <= sk(j+1) ) then
        exit
      end if

      if ( nk <= j ) then
        exit
      end if

      j = j + 1

    end do

    jp1 = i4_wrap ( j + 1, 1, nk )

    pt(1:dim_num,it) = ( ( sk(j+1) - st         ) * pk(1:dim_num,j) &
                       + (           st - sk(j) ) * pk(1:dim_num,jp1) ) &
                       / ( sk(j+1)      - sk(j) )

  end do

  return
end
subroutine provec ( m, n, base, vecm, vecn, vecnm )

!*****************************************************************************80
!
!! PROVEC projects a vector from M space into N space.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the dimension of the higher order space.
!
!    Input, integer N, the dimension of the lower order space.
!
!    Input, real ( kind = rk ) BASE(M,N).  The columns of BASE contain
!    N vectors, each of length M, which form the basis for
!    a space of dimension N.
!
!    Input, real ( kind = rk ) VECM(M), is an M dimensional vector.
!
!    Output, real ( kind = rk ) VECN(N), the projection of VECM into the
!    lower dimensional space.  These values represent
!    coordinates in the lower order space.
!
!    Output, real ( kind = rk ) VECNM(M), the projection of VECM into the
!    lower dimensional space, but using coordinates in
!    the higher dimensional space.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) base(m,n)
  integer i
  integer j
  real ( kind = rk ) temp
  real ( kind = rk ) vecm(m)
  real ( kind = rk ) vecn(n)
  real ( kind = rk ) vecnm(m)
!
!  For each vector, remove all projections onto previous vectors,
!  and then normalize.  This should result in a matrix BASE
!  whose columns are orthonormal.
!
  do j = 1, n

    do i = 1, j-1

      temp = dot_product ( base(1:m,i), base(1:m,j) )

      base(1:m,j) = base(1:m,j) - temp * base(1:m,i)

    end do

    temp = sqrt ( sum ( base(1:m,j)**2 ) )

    if ( 0.0D+00 < temp ) then
      base(1:m,j) = base(1:m,j) / temp
    end if

  end do
!
!  Compute the coordinates of the projection of the vector
!  simply by taking dot products.
!
  do j = 1, n
    vecn(j) = dot_product ( vecm(1:m), base(1:m,j) )
  end do
!
!  Compute the coordinates of the projection in terms of
!  the original space.
!
  do i = 1, m
    vecnm(i) = dot_product ( base(i,1:n), vecn(1:n) )
  end do

  return
end
subroutine pyramid_volume_3d ( h, s, volume )

!*****************************************************************************80
!
!! PYRAMID_VOLUME_3D computes the volume of a pyramid with square base in 3D.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) H, S, the height of the pyramid, and the 
!    length of one side of the square base.
!
!    Output, real ( kind = rk ) VOLUME, the volume of the pyramid.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) h
  real ( kind = rk ) s
  real ( kind = rk ) volume

  volume = s * s * h / 3.0D+00

  return
end
function pyramid01_volume ( )

!*****************************************************************************80
!
!! PYRAMID01_VOLUME returns the volume of a unit pyramid.
!
!  Discussion:
!
!    A pyramid with square base can be regarded as the upper half of a
!    3D octahedron.
!
!    The integration region:
!
!      - ( 1 - Z ) <= X <= 1 - Z
!      - ( 1 - Z ) <= Y <= 1 - Z
!                0 <= Z <= 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 January 2018
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) VALUE, the volume of the pyramid.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) pyramid01_volume
  real ( kind = rk ) volume

  volume = 4.0D+00 / 3.0D+00

  pyramid01_volume = volume

  return
end
subroutine quad_area_2d ( q, area )

!*****************************************************************************80
!
!! QUAD_AREA_2D computes the area of a quadrilateral in 2D.
!
!  Discussion:
!
!    A quadrilateral is a polygon defined by 4 vertices.
!
!    This algorithm should be able to handle nonconvex quadrilaterals.
!
!    The vertices of the quadrilateral should be listed in counter clockwise
!    order, so that the area is positive.
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
!  Parameters:
!
!    Input, real ( kind = rk ) Q(2,4), the vertices, specified in
!    counter clockwise order.
!
!    Output, real ( kind = rk ) AREA, the area of the quadrilateral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) area
  real ( kind = rk ) area_triangle
  real ( kind = rk ) q(dim_num,4)
  real ( kind = rk ) t(dim_num,3)

  area = 0.0D+00

  t(1:dim_num,1:3) = reshape ( (/ &
    q(1:2,1), q(1:2,2), q(1:2,3) /), (/ dim_num, 3 /) )

  call triangle_area_2d ( t, area_triangle )

  area = area + area_triangle

  t(1:dim_num,1:3) = reshape ( (/ &
    q(1:2,3), q(1:2,4), q(1:2,1) /), (/ dim_num, 3 /) )

  call triangle_area_2d ( t, area_triangle )

  area = area + area_triangle

  return
end
subroutine quad_area2_2d ( q, area )

!*****************************************************************************80
!
!! QUAD_AREA2_2D computes the area of a quadrilateral in 2D.
!
!  Discussion:
!
!    A quadrilateral is a polygon defined by 4 vertices.
!
!    This algorithm computes the area of the related 
!    Varignon parallelogram first.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 May 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) Q(2,4), the vertices, specified in
!    counter clockwise order.
!
!    Output, real ( kind = rk ) AREA, the area of the quadrilateral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) area
  real ( kind = rk ) p(2,4)
  real ( kind = rk ) q(2,4)
!
!  Define a parallelogram by averaging consecutive vertices.
!
  p(1:2,1:3) = ( q(1:2,1:3) + q(1:2,2:4) ) / 2.0D+00
  p(1:2,  4) = ( q(1:2,  4) + q(1:2,1  ) ) / 2.0D+00
!
!  Compute the area.
!
  call parallelogram_area_2d ( p, area )
!
!  The quadrilateral's area is twice that of the parallelogram.
!
  area = 2.0D+00 * area

  return
end
subroutine quad_area_3d ( q, area )

!*****************************************************************************80
!
!! QUAD_AREA_3D computes the area of a quadrilateral in 3D.
!
!  Discussion:
!
!    A quadrilateral is a polygon defined by 4 vertices.
!
!    It is assumed that the four vertices of the quadrilateral
!    are coplanar.
!
!    This algorithm computes the area of the related 
!    Varignon parallelogram first.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 May 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) Q(3,4), the vertices, specified in
!    counter clockwise order.
!
!    Output, real ( kind = rk ) AREA, the area of the quadrilateral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) area
  real ( kind = rk ) p(3,4)
  real ( kind = rk ) q(3,4)
!
!  Define a parallelogram by averaging consecutive vertices.
!
  p(1:3,1:3) = ( q(1:3,1:3) + q(1:3,2:4) ) / 2.0D+00
  p(1:3,  4) = ( q(1:3,  4) + q(1:3,1  ) ) / 2.0D+00
!
!  Compute the area.
!
  call parallelogram_area_3d ( p, area )
!
!  The quadrilateral's area is twice that of the parallelogram.
!
  area = 2.0D+00 * area

  return
end
subroutine quad_contains_point_2d ( q, p, inside )

!*****************************************************************************80
!
!! QUAD_CONTAINS_POINT_2D: is point inside a convex quadrilateral in 2D.
!
!  Discussion:
!
!    A quadrilateral is a polygon defined by 4 vertices.
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
!  Parameters:
!
!    Input, real ( kind = rk ) Q(2,4), the vertices of the quadrilateral.
!
!    Input, real ( kind = rk ) P(2), the point to be checked.
!
!    Output, logical INSIDE, is TRUE if the point is in the
!    quadrilateral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) angle_1
  real ( kind = rk ) angle_2
  real ( kind = rk ) angle_rad_2d
  logical inside
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) q(dim_num,4)
!
!  This will only handle convex quadrilaterals.
!
  inside = .false.

  angle_1 = angle_rad_2d ( q(1:2,1), q(1:2,2), q(1:2,3) )
  angle_2 = angle_rad_2d ( q(1:2,1), q(1:2,2), p(1:2) )

  if ( angle_1 < angle_2 ) then
    return
  end if

  angle_1 = angle_rad_2d ( q(1:2,2), q(1:2,3), q(1:2,4) )
  angle_2 = angle_rad_2d ( q(1:2,2), q(1:2,3), p(1:2) )

  if ( angle_1 < angle_2 ) then
    return
  end if

  angle_1 = angle_rad_2d ( q(1:2,3), q(1:2,4), q(1:2,1) )
  angle_2 = angle_rad_2d ( q(1:2,3), q(1:2,4), p(1:2) )

  if ( angle_1 < angle_2 ) then
    return
  end if

  angle_1 = angle_rad_2d ( q(1:2,4), q(1:2,1), q(1:2,2) )
  angle_2 = angle_rad_2d ( q(1:2,4), q(1:2,1), p(1:2) )

  if ( angle_1 < angle_2 ) then
    return
  end if

  inside = .true.

  return
end
subroutine quad_convex_random ( xy )

!*****************************************************************************80
!
!! QUAD_CONVEX_RANDOM returns a random convex quadrilateral.
!
!  Description:
!
!    A quadrilateral is a polygon defined by 4 vertices.
!
!    The quadrilateral is constrained in that the vertices must all lie
!    with the unit square.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 June 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) XY(2,NODE_NUM), the coordinates of the 
!    nodes of the quadrilateral, given in counterclockwise order.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: node_num = 4

  integer hull(node_num)
  integer hull_num
  integer j
  real ( kind = rk ) xy(2,node_num)
  real ( kind = rk ) xy_random(2,node_num)

  do
!
!  Generate 4 random points.
!
    call random_number ( harvest = xy_random(1:2,1:node_num) )
!
!  Determine the convex hull.
!
    call points_hull_2d ( node_num, xy_random, hull_num, hull )
!
!  If HULL_NUM < NODE_NUM, then our convex hull is a triangle.
!  Try again.
!
    if ( hull_num == node_num ) then
      exit
    end if

  end do
!
!  Make an ordered copy of the random points.
!
  do j = 1, node_num
    xy(1:2,j) = xy_random(1:2,hull(j))
  end do

  return
end
subroutine quad_point_dist_2d ( q, p, dist )

!*****************************************************************************80
!
!! QUAD_POINT_DIST_2D: distance ( quadrilateral, point ) in 2D.
!
!  Discussion:
!
!    A quadrilateral is a polygon defined by 4 vertices.
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
!  Parameters:
!
!    Input, real ( kind = rk ) Q(2,4), the quadrilateral vertices.
!
!    Input, real ( kind = rk ) P(2), the point to be checked.
!
!    Output, real ( kind = rk ) DIST, the distance from the point to the
!    quadrilateral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2
  integer, parameter :: side_num = 4

  real ( kind = rk ) dist
  real ( kind = rk ) dist2
  integer i4_wrap
  integer j
  integer jp1
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) q(dim_num,side_num)
!
!  Find the distance to each of the line segments.
!
  dist = huge ( dist )

  do j = 1, side_num

    jp1 = i4_wrap ( j+1, 1, side_num )

    call segment_point_dist_2d ( q(1:dim_num,j), q(1:dim_num,jp1), p, dist2 )

    if ( dist2 < dist ) then
      dist = dist2
    end if

  end do

  return
end
subroutine quad_point_dist_signed_2d ( q, p, dist_signed )

!*****************************************************************************80
!
!! QUAD_POINT_DIST_SIGNED_2D: signed distance ( quadrilateral, point ) in 2D.
!
!  Discussion:
!
!    A quadrilateral is a polygon defined by 4 vertices.
!
!    The quadrilateral must be convex.  DIST_SIGNED is actually the maximum 
!    of the signed distances from the point to each of the four lines that 
!    make up the quadrilateral.
!
!    Essentially, if the point is outside the convex quadrilateral,
!    only one of the signed distances can be positive, or two can
!    be positive and equal.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) Q(2,4), the vertices of the quadrilateral.
!
!    Input, real ( kind = rk ) P(2), the point which is to be checked.
!
!    Output, real ( kind = rk ) DIST_SIGNED, the signed distance from the 
!    point to the convex quadrilateral.  If DIST_SIGNED is
!    0.0, the point is on the boundary;
!    negative, the point is in the interior;
!    positive, the point is in the exterior.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) dis
  real ( kind = rk ) dis12
  real ( kind = rk ) dis23
  real ( kind = rk ) dis34
  real ( kind = rk ) dis41
  real ( kind = rk ) dist_signed
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) pm(dim_num)
  real ( kind = rk ) q(dim_num,4)
!
!  Compare the signed distance from each line segment to the point,
!  with the signed distance to the midpoint of the opposite line.
!
!  The signed distances should all be negative if the point is inside.
!
!  Side 12
!
  call line_exp_point_dist_signed_2d ( q(1:2,1), q(1:2,2), p, dis12 )

  pm(1:dim_num) = 0.5D+00 * ( q(1:dim_num,3) + q(1:dim_num,4) )

  call line_exp_point_dist_signed_2d ( q(1:2,1), q(1:2,2), pm, dis )

  if ( 0.0D+00 < dis ) then
    dis = -dis
    dis12 = -dis12
  end if
!
!  Side 23
!
  call line_exp_point_dist_signed_2d ( q(1:2,2), q(1:2,3), p, dis23 )

  pm(1:dim_num) = 0.5D+00 * ( q(1:dim_num,4) + q(1:dim_num,1) )

  call line_exp_point_dist_signed_2d ( q(1:2,2), q(1:2,3), pm, dis )

  if ( 0.0D+00 < dis ) then
    dis = -dis
    dis23 = -dis23
  end if
!
!  Side 34
!
  call line_exp_point_dist_signed_2d ( q(1:2,3), q(1:2,4), p, dis34 )

  pm(1:dim_num) = 0.5D+00 * ( q(1:dim_num,1) + q(1:dim_num,2) )

  call line_exp_point_dist_signed_2d ( q(1:2,3), q(1:2,4), pm, dis )

  if ( 0.0D+00 < dis ) then
    dis = -dis
    dis34 = -dis34
  end if
!
!  Side 41
!
  call line_exp_point_dist_signed_2d ( q(1:2,4), q(1:2,1), p, dis41 )

  pm(1:dim_num) = 0.5D+00 * ( q(1:dim_num,2) + q(1:dim_num,3) )

  call line_exp_point_dist_signed_2d ( q(1:2,4), q(1:2,1), pm, dis )

  if ( 0.0D+00 < dis ) then
    dis = -dis
    dis41 = -dis41
  end if

  dist_signed = max ( dis12, dis23, dis34, dis41 )

  return
end
subroutine quad_point_near_2d ( q, p, pn, dist )

!*****************************************************************************80
!
!! QUAD_POINT_NEAR_2D computes the nearest point on a quadrilateral in 2D.
!
!  Discussion:
!
!    A quadrilateral is a polygon defined by 4 vertices.
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
!  Parameters:
!
!    Input, real ( kind = rk ) Q(2,4), the quadrilateral vertices.
!
!    Input, real ( kind = rk ) P(2), the point whose nearest quadrilateral point
!    is to be determined.
!
!    Output, real ( kind = rk ) PN(2), the nearest point to P.
!
!    Output, real ( kind = rk ) DIST, the distance from the point to the
!    quadrilateral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2
  integer, parameter :: side_num = 4

  real ( kind = rk ) dist
  real ( kind = rk ) dist2
  integer i4_wrap
  integer j
  integer jp1
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) pn(dim_num)
  real ( kind = rk ) pn2(dim_num)
  real ( kind = rk ) q(dim_num,side_num)
  real ( kind = rk ) tval
!
!  Find the distance to each of the line segments that make up the edges
!  of the quadrilateral.
!
  dist = huge ( dist )
  pn(1:dim_num) = 0.0D+00

  do j = 1, side_num

    jp1 = i4_wrap ( j+1, 1, side_num )

    call segment_point_near_2d ( q(1:dim_num,j), q(1:dim_num,jp1), p, &
      pn2, dist2, tval )

    if ( dist2 < dist ) then
      dist = dist2
      pn(1:dim_num) = pn2(1:dim_num)
    end if

  end do

  return
end
function r8_acos ( c )

!*****************************************************************************80
!
!! R8_ACOS computes the arc cosine function, with argument truncation.
!
!  Discussion:
!
!    If you call your system ACOS routine with an input argument that is
!    even slightly outside the range [-1.0, 1.0 ], you may get an unpleasant
!    surprise (I did).
!
!    This routine simply truncates arguments outside the range.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) C, the argument.
!
!    Output, real ( kind = rk ) R8_ACOS, an angle whose cosine is C.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) c
  real ( kind = rk ) c2
  real ( kind = rk ) r8_acos

  c2 = c
  c2 = max ( c2, -1.0D+00 )
  c2 = min ( c2, +1.0D+00 )

  r8_acos = acos ( c2 )

  return
end
function r8_asin ( s )

!*****************************************************************************80
!
!! R8_ASIN computes the arc sine function, with argument truncation.
!
!  Discussion:
!
!    If you call your system ASIN routine with an input argument that is
!    even slightly outside the range [-1.0, 1.0 ], you may get an unpleasant 
!    surprise (I did).
!
!    This routine simply truncates arguments outside the range.
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
!  Parameters:
!
!    Input, real ( kind = rk ) S, the argument.
!
!    Output, real ( kind = rk ) R8_ASIN, an angle whose sine is S.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r8_asin
  real ( kind = rk ) s
  real ( kind = rk ) s2

  s2 = s
  s2 = max ( s2, -1.0D+00 )
  s2 = min ( s2, +1.0D+00 )

  r8_asin = asin ( s2 )

  return
end
function r8_atan ( y, x )

!*****************************************************************************80
!
!! R8_ATAN computes the inverse tangent of the ratio Y / X.
!
!  Discussion:
!
!    R8_ATAN returns an angle whose tangent is ( Y / X ), a job which
!    the built in functions ATAN and ATAN2 already do.
!
!    However:
!
!    * R8_ATAN always returns a positive angle, between 0 and 2 PI,
!      while ATAN and ATAN2 return angles in the interval [-PI/2,+PI/2]
!      and [-PI,+PI] respectively;
!
!    * R8_ATAN accounts for the signs of X and Y, (as does ATAN2).  The ATAN
!     function by contrast always returns an angle in the first or fourth
!     quadrants.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) Y, X, two quantities which represent the
!    tangent of an angle.  If Y is not zero, then the tangent is (Y/X).
!
!    Output, real ( kind = rk ) R8_ATAN, an angle between 0 and 2 * PI, whose
!    tangent is (Y/X), and which lies in the appropriate quadrant so that
!    the signs of its cosine and sine match those of X and Y.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) abs_x
  real ( kind = rk ) abs_y
  real ( kind = rk ) r8_atan
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) theta
  real ( kind = rk ) theta_0
  real ( kind = rk ) x
  real ( kind = rk ) y
!
!  Special cases:
!
  if ( x == 0.0D+00 ) then

    if ( 0.0D+00 < y ) then
      theta = r8_pi / 2.0D+00
    else if ( y < 0.0D+00 ) then
      theta = 3.0D+00 * r8_pi / 2.0D+00
    else if ( y == 0.0D+00 ) then
      theta = 0.0D+00
    end if

  else if ( y == 0.0D+00 ) then

    if ( 0.0D+00 < x ) then
      theta = 0.0D+00
    else if ( x < 0.0D+00 ) then
      theta = r8_pi
    end if
!
!  We assume that ATAN2 is correct when both arguments are positive.
!
  else

    abs_y = abs ( y )
    abs_x = abs ( x )

    theta_0 = atan2 ( abs_y, abs_x )

    if ( 0.0D+00 < x .and. 0.0D+00 < y ) then
      theta = theta_0
    else if ( x < 0.0D+00 .and. 0.0D+00 < y ) then
      theta = r8_pi - theta_0
    else if ( x < 0.0D+00 .and. y < 0.0D+00 ) then
      theta = r8_pi + theta_0
    else if ( 0.0D+00 < x .and. y < 0.0D+00 ) then
      theta = 2.0D+00 * r8_pi - theta_0
    end if

  end if

  r8_atan = theta

  return
end
function r8_cosd ( degrees )

!*****************************************************************************80
!
!! R8_COSD returns the cosine of an angle given in degrees.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) DEGREES, the angle in degrees.
!
!    Output, real ( kind = rk ) R8_COSD, the cosine of the angle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) degrees
  real ( kind = rk ) r8_cosd
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) radians

  radians = r8_pi * ( degrees / 180.0D+00 )
  r8_cosd = cos ( radians )

  return
end
function r8_cotd ( degrees )

!*****************************************************************************80
!
!! R8_COTD returns the cotangent of an angle given in degrees.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) DEGREES, the angle in degrees.
!
!    Output, real ( kind = rk ) R8_COTD, the cotangent of the angle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) degrees
  real ( kind = rk ) r8_cotd
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) radians

  radians = r8_pi * ( degrees / 180.0D+00 )
  r8_cotd = cos ( radians ) / sin ( radians )

  return
end
function r8_cscd ( degrees )

!*****************************************************************************80
!
!! R8_CSCD returns the cosecant of an angle given in degrees.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) DEGREES, the angle in degrees.
!
!    Output, real ( kind = rk ) R8_CSCD, the cosecant of the angle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) degrees
  real ( kind = rk ) r8_cscd
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) radians

  radians = r8_pi * ( degrees / 180.0D+00 )
  r8_cscd  = 1.0D+00 / sin ( radians )

  return
end
subroutine r8_fake_use ( x )

!*****************************************************************************80
!
!! r8_fake_use pretends to use a variable.
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
function r8_huge ( )

!*****************************************************************************80
!
!! R8_HUGE returns a very large R8.
!
!  Discussion:
!
!    The value returned by this function is NOT required to be the
!    maximum representable R8.  This value varies from machine to machine,
!    from compiler to compiler, and may cause problems when being printed.
!    We simply want a "very large" but non-infinite number.
!
!    FORTRAN90 provides a built-in routine HUGE ( X ) that
!    can return the maximum representable number of the same datatype
!    as X, if that is what is really desired.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 October 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) R8_HUGE, a "huge" value.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r8_huge

  r8_huge = 1.0D+30

  return
end
function r8_is_int ( r )

!*****************************************************************************80
!
!! R8_IS_INT determines if a real number represents an integer value.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) R, the number to be checked.
!
!    Output, logical R8_IS_INT, is TRUE if R is an integer value.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  real ( kind = rk ) r
  logical r8_is_int

  if ( real ( huge ( i ), kind = rk ) < r ) then
    r8_is_int = .false.
  else if ( r < - real ( huge ( i ), kind = rk ) ) then
    r8_is_int = .false.
  else if ( r == real ( int ( r ), kind = rk ) ) then
    r8_is_int = .true.
  else
    r8_is_int = .false.
  end if

  return
end
function r8_modp ( x, y )

!*****************************************************************************80
!
!! r8_modp() returns the nonnegative remainder of real division.
!
!  Discussion:
!
!    If
!      REM = R8_MODP ( X, Y )
!      RMULT = ( X - REM ) / Y
!    then
!      X = Y * RMULT + REM
!    where REM is always nonnegative.
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360.0) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, R8_MODP(A,360.0) is between 0 and 360, always.
!
!  Example:
!
!        I         J     MOD  R8_MODP  R8_MODP Factorization
!
!      107        50       7       7    107 =  2 *  50 + 7
!      107       -50       7       7    107 = -2 * -50 + 7
!     -107        50      -7      43   -107 = -3 *  50 + 43
!     -107       -50      -7      43   -107 =  3 * -50 + 43
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    24 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) X, the number to be divided.
!
!    real ( kind = rk ) Y, the number that divides X.
!
!  Output:
!
!    real ( kind = rk ) R8_MODP, the nonnegative remainder 
!    when X is divided by Y.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r8_modp
  real ( kind = rk ) x
  real ( kind = rk ) y

  if ( y == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_MODP - Fatal error!'
    write ( *, '(a,g14.6)' ) '  R8_MODP ( X, Y ) called with Y = ', y
    stop 1
  end if

  r8_modp = mod ( x, y )

  if ( r8_modp < 0.0D+00 ) then
    r8_modp = r8_modp + abs ( y )
  end if

  return
end
function r8_normal_01 ( )

!*****************************************************************************80
!
!! R8_NORMAL_01 returns a unit pseudonormal R8.
!
!  Discussion:
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!    Because this routine uses the Box Muller method, it requires pairs
!    of uniform random values to generate a pair of normal random values.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) R8_NORMAL_01, a sample of the standard 
!    normal PDF.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r1
  real ( kind = rk ) r2
  real ( kind = rk ) r8_normal_01
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  integer, save :: used = 0
  real ( kind = rk ) x
  real ( kind = rk ), save :: y = 0.0D+00
!
!  On odd numbered calls, generate two uniforms, create two normals,
!  return the first normal.
!
  if ( mod ( used, 2 ) == 0 ) then

    call random_number ( harvest = r1 )
    call random_number ( harvest = r2 )

    x = sqrt ( -2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * r8_pi * r2 )
    y = sqrt ( -2.0D+00 * log ( r1 ) ) * sin ( 2.0D+00 * r8_pi * r2 )
!
!  On odd calls, return the second normal.
!
  else

    x = y

  end if

  used = used + 1

  r8_normal_01 = x

  return
end
function r8_pi ( )

!*****************************************************************************80
!
!! R8_PI returns the value of pi.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) R8_PI, the value of pi.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r8_pi

  r8_pi = 3.141592653589793D+00

  return
end
function r8_sign_opposite_strict ( r1, r2 )

!*****************************************************************************80
!
!! R8_SIGN_OPPOSITE_STRICT is TRUE if two R8's are strictly of opposite sign.
!
!  Discussion:
!
!    This test could be coded numerically as
!
!      if ( r1 * r2 < 0.0 ) then ...
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) R1, R2, the values to check.
!
!    Output, logical R8_SIGN_OPPOSITE_STRICT, is TRUE 
!    if ( R1 < 0 and 0 < R2 ) or ( R2 < 0 and 0 < R1 ).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r1
  real ( kind = rk ) r2
  logical r8_sign_opposite_strict

  r8_sign_opposite_strict = ( r1 < 0.0D+00 .and. 0.0D+00 < r2 ) .or. &
                            ( r2 < 0.0D+00 .and. 0.0D+00 < r1 )

  return
end
function r8_sind ( degrees )

!*****************************************************************************80
!
!! R8_SIND returns the sine of an angle given in degrees.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) DEGREES, the angle in degrees.
!
!    Output, real ( kind = rk ) R8_SIND, the sine of the angle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) degrees
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) r8_sind
  real ( kind = rk ) radians

  radians = r8_pi * ( degrees / 180.0D+00 )
  r8_sind  = sin ( radians )

  return
end
function r8_secd ( degrees )

!*****************************************************************************80
!
!! R8_SECD returns the secant of an angle given in degrees.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) DEGREES, the angle in degrees.
!
!    Output, real ( kind = rk ) R8_SECD, the secant of the angle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) degrees
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) r8_secd
  real ( kind = rk ) radians

  radians = r8_pi * ( degrees / 180.0D+00 )
  r8_secd = 1.0D+00 / cos ( radians )

  return
end
subroutine r8_swap ( x, y )

!*****************************************************************************80
!
!! R8_SWAP switches two R8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = rk ) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) x
  real ( kind = rk ) y
  real ( kind = rk ) z

  z = x
  x = y
  y = z

  return
end
function r8_tand ( degrees )

!*****************************************************************************80
!
!! R8_TAND returns the tangent of an angle given in degrees.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) DEGREES, the angle in degrees.
!
!    Output, real ( kind = rk ) R8_TAND, the tangent of the angle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) degrees
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) r8_tand
  real ( kind = rk ) radians

  radians = r8_pi * ( degrees / 180.0D+00 )
  r8_tand  = tan ( radians )

  return
end
function r8_uniform_ab ( a, b )

!*****************************************************************************80
!
!! R8_UNIFORM_AB returns a scaled pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = rk ) value.
!
!    The pseudorandom number should be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, the limits of the interval.
!
!    Output, real ( kind = rk ) R8_UNIFORM_AB, a number strictly between A and B.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) r
  real ( kind = rk ) r8_uniform_ab

  call random_number ( harvest = r )
  r8_uniform_ab = a + ( b - a ) * r

  return
end
subroutine r82vec_permute ( n, p, a )

!*****************************************************************************80
!
!! R82VEC_PERMUTE permutes a R82 vector in place.
!
!  Discussion:
!
!    This routine permutes an array of real "objects", but the same
!    logic can be used to permute an array of objects of any arithmetic
!    type, or an array of objects of any complexity.  The only temporary
!    storage required is enough to store a single object.  The number
!    of data movements made is N + the number of cycles of order 2 or more,
!    which is never more than N + N/2.
!
!  Example:
!
!    Input:
!
!      N = 5
!      P = (   2,    4,    5,    1,    3 )
!      A = ( 1.0,  2.0,  3.0,  4.0,  5.0 )
!          (11.0, 22.0, 33.0, 44.0, 55.0 )
!
!    Output:
!
!      A    = (  2.0,  4.0,  5.0,  1.0,  3.0 )
!             ( 22.0, 44.0, 55.0, 11.0, 33.0 ).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    28 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of objects.
!
!    Input, integer P(N), the permutation.  P(I) = J means
!    that the I-th element of the output array should be the J-th
!    element of the input array.  P must be a legal permutation
!    of the integers from 1 to N, otherwise the algorithm will
!    fail catastrophically.
!
!    Input/output, real ( kind = rk ) A(2,N), the array to be permuted.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(2,n)
  real ( kind = rk ) a_temp(2)
  integer iget
  integer iput
  integer istart
  integer p(n)
!
!  Search for the next element of the permutation that has not been used.
!
  do istart = 1, n

    if ( p(istart) < 0 ) then

      cycle

    else if ( p(istart) == istart ) then

      p(istart) = - p(istart)
      cycle

    else

      a_temp(1:2) = a(1:2,istart)
      iget = istart
!
!  Copy the new value into the vacated entry.
!
      do

        iput = iget
        iget = p(iget)

        p(iput) = - p(iput)

        if ( iget < 1 .or. n < iget ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R82VEC_PERMUTE - Fatal error!'
          stop 1
        end if

        if ( iget == istart ) then
          a(1:2,iput) = a_temp(1:2)
          exit
        end if

        a(1:2,iput) = a(1:2,iget)

      end do

    end if

  end do
!
!  Restore the signs of the entries.
!
  p(1:n) = -p(1:n)

  return
end
subroutine r82vec_sort_heap_index_a ( n, a, indx )

!*****************************************************************************80
!
!! R82VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R82VEC.
!
!  Discussion:
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      A(1:2,INDX(I)), I = 1 to N is sorted,
!
!    or explicitly, by the call
!
!      call R82VEC_PERMUTE ( N, A, INDX )
!
!    after which A(1:2,I), I = 1 to N is sorted.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input, real ( kind = rk ) A(2,N), an array to be index-sorted.
!
!    Output, integer INDX(N), the sort index.  The
!    I-th element of the sorted array is A(1:2,INDX(I)).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(2,n)
  real ( kind = rk ) aval(2)
  integer i
  integer indx(n)
  integer indxt
  integer ir
  integer j
  integer l

  if ( n < 1 ) then
    return
  end if

  if ( n == 1 ) then
    indx(1) = 1
    return
  end if

  call i4vec_indicator ( n, indx )

  l = n / 2 + 1
  ir = n

  do

    if ( 1 < l ) then

      l = l - 1
      indxt = indx(l)
      aval(1:2) = a(1:2,indxt)

    else

      indxt = indx(ir)
      aval(1:2) = a(1:2,indxt)
      indx(ir) = indx(1)
      ir = ir - 1

      if ( ir == 1 ) then
        indx(1) = indxt
        exit
      end if

    end if

    i = l
    j = l + l

    do while ( j <= ir )

      if ( j < ir ) then
        if (   a(1,indx(j)) <  a(1,indx(j+1)) .or. &
             ( a(1,indx(j)) == a(1,indx(j+1)) .and. &
               a(2,indx(j)) <  a(2,indx(j+1)) ) ) then
          j = j + 1
        end if
      end if

      if (   aval(1) <  a(1,indx(j)) .or. &
           ( aval(1) == a(1,indx(j)) .and. &
             aval(2) <  a(2,indx(j)) ) ) then
        indx(i) = indx(j)
        i = j
        j = j + j
      else
        j = ir + 1
      end if

    end do

    indx(i) = indxt

  end do

  return
end
subroutine r8ge_det ( n, a, pivot, det )

!*****************************************************************************80
!
!! R8GE_DET computes the determinant of a matrix factored by R8GE_FA.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    29 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jack Dongarra, James Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real ( kind = rk ) A(N,N), the LU factors computed by R8GE_FA.
!
!    Input, integer PIVOT(N), as computed by R8GE_FA.
!
!    Output, real ( kind = rk ) DET, the determinant of the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) det
  integer i
  integer pivot(n)

  det = 1.0D+00

  do i = 1, n
    det = det * a(i,i)
    if ( pivot(i) /= i ) then
      det = - det
    end if
  end do

  return
end
subroutine r8ge_fa ( n, a, pivot, info )

!*****************************************************************************80
!
!! R8GE_FA factors a general matrix.
!
!  Discussion:
!
!    R8GE_FA is a simplified version of the LINPACK routine DGEFA.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jack Dongarra, James Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input/output, real ( kind = rk ) A(N,N), the matrix to be factored.
!    On output, A contains an upper triangular matrix and the multipliers
!    which were used to obtain it.  The factorization can be written
!    A = L * U, where L is a product of permutation and unit lower
!    triangular matrices and U is upper triangular.
!
!    Output, integer PIVOT(N), a vector of pivot indices.
!
!    Output, integer INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,n)
  integer i
  integer info
  integer pivot(n)
  integer j
  integer k
  integer l
  real ( kind = rk ) t

  info = 0

  do k = 1, n - 1
!
!  Find L, the index of the pivot row.
!
    l = k
    do i = k+1, n
      if ( abs ( a(l,k) ) < abs ( a(i,k) ) ) then
        l = i
      end if
    end do

    pivot(k) = l
!
!  If the pivot index is zero, the algorithm has failed.
!
    if ( a(l,k) == 0.0D+00 ) then
      info = k
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8GE_FA - Warning!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      return
    end if
!
!  Interchange rows L and K if necessary.
!
    if ( l /= k ) then
      t      = a(l,k)
      a(l,k) = a(k,k)
      a(k,k)= t
    end if
!
!  Normalize the values that lie below the pivot entry A(K,K).
!
    a(k+1:n,k) = - a(k+1:n,k) / a(k,k)
!
!  Row elimination with column indexing.
!
    do j = k+1, n

      if ( l /= k ) then
        t      = a(l,j)
        a(l,j) = a(k,j)
        a(k,j) = t
      end if

      a(k+1:n,j) = a(k+1:n,j) + a(k+1:n,k) * a(k,j)

    end do

  end do

  pivot(n) = n

  if ( a(n,n) == 0.0D+00 ) then
    info = n
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_FA - Warning!'
    write ( *, '(a,i8)' ) '  Zero pivot on step ', info
  end if

  return
end
subroutine r8ge_sl ( n, a, pivot, b, job )

!*****************************************************************************80
!
!! R8GE_SL solves a system factored by R8GE_FA.
!
!  Discussion:
!
!    R8GE_SL is a simplified version of the LINPACK routine DGESL.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real ( kind = rk ) A(N,N), the LU factors from R8GE_FA.
!
!    Input, integer PIVOT(N), the pivot vector from R8GE_FA.
!
!    Input/output, real ( kind = rk ) B(N).
!    On input, the right hand side vector.
!    On output, the solution vector.
!
!    Input, integer JOB, specifies the operation.
!    0, solve A * x = b.
!    nonzero, solve A' * x = b.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) b(n)
  integer pivot(n)
  integer job
  integer k
  integer l
!
!  Solve A * x = b.
!
  if ( job == 0 ) then
!
!  Solve PL * Y = B.
!
    do k = 1, n - 1

      l = pivot(k)

      if ( l /= k ) then
        call r8_swap ( b(l), b(k) )
      end if

      b(k+1:n) = b(k+1:n) + a(k+1:n,k) * b(k)

    end do
!
!  Solve U * X = Y.
!
    do k = n, 1, -1
      b(k) = b(k) / a(k,k)
      b(1:k-1) = b(1:k-1) - a(1:k-1,k) * b(k)
    end do
!
!  Solve A' * X = B.
!
  else
!
!  Solve U' * Y = B.
!
    do k = 1, n
      b(k) = ( b(k) - sum ( b(1:k-1) * a(1:k-1,k) ) ) / a(k,k)
    end do
!
!  Solve ( PL )' * X = Y.
!
    do k = n - 1, 1, -1

      b(k) = b(k) + sum ( b(k+1:n) * a(k+1:n,k) )

      l = pivot(k)

      if ( l /= k ) then
        call r8_swap ( b(l), b(k) )
      end if

    end do

  end if

  return
end
function r8mat_det_2d ( a )

!*****************************************************************************80
!
!! R8MAT_DET_2D computes the determinant of a 2 by 2 matrix.
!
!  Discussion:
!
!    The determinant is the area spanned by the vectors making up the rows
!    or columns of the matrix.
!
!    R8MAT_DET_2D = A(1,1) * A(2,2) - A(1,2) * A(2,1).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A(2,2), the matrix whose determinant is desired.
!
!    Output, real ( kind = rk ) R8MAT_DET_2D, the determinant of the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(2,2)
  real ( kind = rk ) r8mat_det_2d

  r8mat_det_2d = a(1,1) * a(2,2) - a(1,2) * a(2,1)

  return
end
function r8mat_det_3d ( a )

!*****************************************************************************80
!
!! R8MAT_DET_3D computes the determinant of a 3 by 3 matrix.
!
!  Discussion:
!
!    The determinant is the volume of the shape spanned by the vectors
!    making up the rows or columns of the matrix.
!
!    det = a11 * a22 * a33 - a11 * a23 * a32
!        + a12 * a23 * a31 - a12 * a21 * a33
!        + a13 * a21 * a32 - a13 * a22 * a31
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A(3,3), the matrix whose determinant is desired.
!
!    Output, real ( kind = rk ) R8MAT_DET_3D, the determinant of the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(3,3)
  real ( kind = rk ) r8mat_det_3d

  r8mat_det_3d =   a(1,1) * ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) &
              + a(1,2) * ( a(2,3) * a(3,1) - a(2,1) * a(3,3) ) &
              + a(1,3) * ( a(2,1) * a(3,2) - a(2,2) * a(3,1) )

  return
end
function r8mat_det_4d ( a )

!*****************************************************************************80
!
!! R8MAT_DET_4D computes the determinant of a 4 by 4 matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A(4,4), the matrix whose determinant is desired.
!
!    Output, real ( kind = rk ) R8MAT_DET_4D, the determinant of the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(4,4)
  real ( kind = rk ) r8mat_det_4d

  r8mat_det_4d = &
      a(1,1) * ( &
        a(2,2) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
      - a(2,3) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
      + a(2,4) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) ) &
    - a(1,2) * ( &
        a(2,1) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
      - a(2,3) * ( a(3,1) * a(4,4) - a(3,4) * a(4,1) ) &
      + a(2,4) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) ) &
    + a(1,3) * ( &
        a(2,1) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
      - a(2,2) * ( a(3,1) * a(4,4) - a(3,4) * a(4,1) ) &
      + a(2,4) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) ) &
    - a(1,4) * ( &
        a(2,1) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) &
      - a(2,2) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) &
      + a(2,3) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) )

  return
end
function r8mat_det_5d ( a )

!*****************************************************************************80
!
!! R8MAT_DET_5D computes the determinant of a 5 by 5 matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A(5,5), the matrix whose determinant is desired.
!
!    Output, real ( kind = rk ) R8MAT_DET_5D, the determinant of the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(5,5)
  real ( kind = rk ) b(4,4)
  real ( kind = rk ) r8mat_det_4d
  real ( kind = rk ) r8mat_det_5d
  integer i
  integer inc
  integer j
  integer k
!
!  Expand the determinant into the sum of the determinants of the
!  five 4 by 4 matrices created by dropping row 1, and column k.
!
  r8mat_det_5d = 0.0D+00

  do k = 1, 5

    do i = 1, 4
      do j = 1, 4

        if ( j < k ) then
          inc = 0
        else
          inc = 1
        end if

        b(i,j) = a(i+1,j+inc)

      end do
    end do

    r8mat_det_5d = r8mat_det_5d + (-1)**( k + 1 ) * a(1,k) * r8mat_det_4d ( b )

  end do

  return
end
subroutine r8mat_inverse_2d ( a, b, det )

!*****************************************************************************80
!
!! R8MAT_INVERSE_2D inverts a 2 by 2 real matrix using Cramer's rule.
!
!  Discussion:
!
!    If DET is zero, then A is singular, and does not have an
!    inverse.  In that case, B is simply set to zero, and a
!    message is printed.
!
!    If DET is nonzero, then its value is roughly an estimate
!    of how nonsingular the matrix A is.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A(2,2), the matrix to be inverted.
!
!    Output, real ( kind = rk ) B(2,2), the inverse of the matrix A.
!
!    Output, real ( kind = rk ) DET, the determinant of the matrix A.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(2,2)
  real ( kind = rk ) b(2,2)
  real ( kind = rk ) det
!
!  Compute the determinant.
!
  det = a(1,1) * a(2,2) - a(1,2) * a(2,1)
!
!  If the determinant is zero, bail out.
!
  if ( det == 0.0D+00 ) then

    b(1:2,1:2) = 0.0D+00

    return
  end if
!
!  Compute the entries of the inverse matrix using an explicit formula.
!
  b(1,1) = + a(2,2) / det
  b(1,2) = - a(1,2) / det
  b(2,1) = - a(2,1) / det
  b(2,2) = + a(1,1) / det

  return
end
subroutine r8mat_inverse_3d ( a, b, det )

!*****************************************************************************80
!
!! R8MAT_INVERSE_3D inverts a 3 by 3 real matrix using Cramer's rule.
!
!  Discussion:
!
!    If DET is zero, then A is singular, and does not have an
!    inverse.  In that case, B is simply set to zero, and a
!    message is printed.
!
!    If DET is nonzero, then its value is roughly an estimate
!    of how nonsingular the matrix A is.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A(3,3), the matrix to be inverted.
!
!    Output, real ( kind = rk ) B(3,3), the inverse of the matrix A.
!
!    Output, real ( kind = rk ) DET, the determinant of the matrix A.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(3,3)
  real ( kind = rk ) b(3,3)
  real ( kind = rk ) det
!
!  Compute the determinant of A
!
  det =   a(1,1) * ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) &
        + a(1,2) * ( a(2,3) * a(3,1) - a(2,1) * a(3,3) ) &
        + a(1,3) * ( a(2,1) * a(3,2) - a(2,2) * a(3,1) )
!
!  If the determinant is zero, bail out.
!
  if ( det == 0.0D+00 ) then

    b(1:3,1:3) = 0.0D+00

    return
  end if
!
!  Compute the entries of the inverse matrix using an explicit
!  formula.
!
  b(1,1) = + ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) / det
  b(1,2) = - ( a(1,2) * a(3,3) - a(1,3) * a(3,2) ) / det
  b(1,3) = + ( a(1,2) * a(2,3) - a(1,3) * a(2,2) ) / det

  b(2,1) = - ( a(2,1) * a(3,3) - a(2,3) * a(3,1) ) / det
  b(2,2) = + ( a(1,1) * a(3,3) - a(1,3) * a(3,1) ) / det
  b(2,3) = - ( a(1,1) * a(2,3) - a(1,3) * a(2,1) ) / det

  b(3,1) = + ( a(2,1) * a(3,2) - a(2,2) * a(3,1) ) / det
  b(3,2) = - ( a(1,1) * a(3,2) - a(1,2) * a(3,1) ) / det
  b(3,3) = + ( a(1,1) * a(2,2) - a(1,2) * a(2,1) ) / det

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints a real matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows in A.
!
!    Input, integer N, the number of columns in A.
!
!    Input, real ( kind = rk ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
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
!! R8MAT_PRINT_SOME prints some of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
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
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = rk ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ILO, JLO, the first row and column to print.
!
!    Input, integer IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
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
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

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

        if ( a(i,j) == real ( int ( a(i,j) ), kind = rk ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine r8mat_solve ( n, rhs_num, a, info )

!*****************************************************************************80
!
!! R8MAT_SOLVE uses Gauss-Jordan elimination to solve an N by N linear system.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    29 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer RHS_NUM, the number of right hand sides.  
!    RHS_NUM must be at least 0.
!
!    Input/output, real ( kind = rk ) A(N,N+rhs_num), contains in rows and
!    columns 1 to N the coefficient matrix, and in columns N+1 through
!    N+rhs_num, the right hand sides.  On output, the coefficient matrix
!    area has been destroyed, while the right hand sides have
!    been overwritten with the corresponding solutions.
!
!    Output, integer INFO, singularity flag.
!    0, the matrix was not singular, the solutions were computed;
!    J, factorization failed on step J, and the solutions could not
!    be computed.
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

  info = 0

  do j = 1, n
!
!  Choose a pivot row.
!
    ipivot = j
    apivot = a(j,j)

    do i = j+1, n
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
!  Interchange.
!
    do i = 1, n + rhs_num
      call r8_swap ( a(ipivot,i), a(j,i) )
    end do
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
subroutine r8mat_solve_2d ( a, b, det, x )

!*****************************************************************************80
!
!! R8MAT_SOLVE_2D solves a 2 by 2 linear system using Cramer's rule.
!
!  Discussion:
!
!    If DET is zero, then A is singular, and does not have an
!    inverse.  In that case, X is simply set to zero, and a
!    message is printed.
!
!    If DET is nonzero, then its value is roughly an estimate
!    of how nonsingular the matrix A is.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A(2,2), the matrix.
!
!    Input, real ( kind = rk ) B(2), the right hand side.
!
!    Output, real ( kind = rk ) DET, the determinant of the matrix A.
!
!    Output, real ( kind = rk ) X(2), the solution of the system, 
!    if DET is nonzero.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(2,2)
  real ( kind = rk ) b(2)
  real ( kind = rk ) det
  real ( kind = rk ) x(2)
!
!  Compute the determinant.
!
  det = a(1,1) * a(2,2) - a(1,2) * a(2,1)
!
!  If the determinant is zero, bail out.
!
  if ( det == 0.0D+00 ) then
    x(1:2) = 0.0D+00
    return
  end if
!
!  Compute the solution.
!
  x(1) = (  a(2,2) * b(1) - a(1,2) * b(2) ) / det
  x(2) = ( -a(2,1) * b(1) + a(1,1) * b(2) ) / det

  return
end
subroutine r8mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT prints a R8MAT, transposed.
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
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = rk ) A(M,N), an M by N matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
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
!! R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
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
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = rk ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ILO, JLO, the first row and column to print.
!
!    Input, integer IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
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
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

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
subroutine r8mat_uniform_ab ( m, n, a, b, r )

!*****************************************************************************80
!
!! R8MAT_UNIFORM_AB fills scaled pseudorandom R8MAT.
!
!  Discussion:
!
!    An R8MAT is an array of real ( kind = rk ) values.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns
!    in the array.
!
!    Input, real ( kind = rk ) A, B, the lower and upper limits.
!
!    Output, real ( kind = rk ) R(M,N), the array of pseudorandom values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) r(m,n)

  call random_number ( harvest = r(1:m,1:n) )

  r(1:m,1:n) = a + ( b - a ) * r(1:m,1:n)

  return
end
subroutine r8vec_angle_3d ( u, v, angle )

!*****************************************************************************80
!
!! R8VEC_ANGLE_3D computes the angle between two vectors in 3D.
!
!  Modified:
!
!    07 July 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) U(3), V(3), the vectors.
!
!    Output, real ( kind = rk ) ANGLE, the angle between the two vectors.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) angle
  real ( kind = rk ) angle_cos
  real ( kind = rk ) r8_acos
  real ( kind = rk ) u(3)
  real ( kind = rk ) u_norm
  real ( kind = rk ) uv_dot
  real ( kind = rk ) v(3)
  real ( kind = rk ) v_norm

  uv_dot = dot_product ( u(1:3), v(1:3) )

  u_norm = sqrt ( dot_product ( u(1:3), u(1:3) ) )

  v_norm = sqrt ( dot_product ( v(1:3), v(1:3) ) )

  angle_cos = uv_dot / u_norm / v_norm

  angle = r8_acos ( angle_cos )

  return
end
subroutine r8vec_any_normal ( dim_num, v1, v2 )

!*****************************************************************************80
!
!! R8VEC_ANY_NORMAL returns some normal vector to V1.
!
!  Discussion:
!
!    If DIM_NUM < 2, then no normal vector can be returned.
!
!    If V1 is the zero vector, then any unit vector will do.
!
!    No doubt, there are better, more robust algorithms.  But I will take
!    just about ANY reasonable unit vector that is normal to V1.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = rk ) V1(DIM_NUM), the vector.
!
!    Output, real ( kind = rk ) V2(DIM_NUM), a vector that is
!    normal to V2, and has unit Euclidean length.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num

  real ( kind = rk ) r8vec_norm
  integer i
  integer j
  integer k
  real ( kind = rk ) v1(dim_num)
  real ( kind = rk ) v2(dim_num)
  real ( kind = rk ) vj
  real ( kind = rk ) vk

  if ( dim_num < 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_ANY_NORMAL - Fatal error!'
    write ( *, '(a)' ) '  Called with DIM_NUM < 2.'
    stop 1
  end if

  if ( r8vec_norm ( dim_num, v1 ) == 0.0D+00 ) then
    v2(1) = 1.0D+00
    v2(2:dim_num) = 0.0D+00
    return
  end if
!
!  Seek the largest entry in V1, VJ = V1(J), and the
!  second largest, VK = V1(K).
!
!  Since V1 does not have zero norm, we are guaranteed that
!  VJ, at least, is not zero.
!
  j = -1
  vj = 0.0D+00

  k = -1
  vk = 0.0D+00

  do i = 1, dim_num

    if ( abs ( vk ) < abs ( v1(i) ) .or. k < 1 ) then

      if ( abs ( vj ) < abs ( v1(i) ) .or. j < 1 ) then
        k = j
        vk = vj
        j = i
        vj = v1(i)
      else
        k = i
        vk = v1(i)
      end if

    end if

  end do
!
!  Setting V2 to zero, except that V2(J) = -VK, and V2(K) = VJ,
!  will just about do the trick.
!
  v2(1:dim_num) = 0.0D+00

  v2(j) = -vk / sqrt ( vk * vk + vj * vj )
  v2(k) =  vj / sqrt ( vk * vk + vj * vj )

  return
end
subroutine r8vec_bracket ( n, x, xval, left, right )

!*****************************************************************************80
!
!! R8VEC_BRACKET searches a sorted array for successive brackets of a value.
!
!  Discussion:
!
!    If the values in the vector are thought of as defining intervals
!    on the real line, then this routine searches for the interval
!    nearest to or containing the given value.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, length of input array.
!
!    Input, real ( kind = rk ) X(N), an array that has been sorted into
!    ascending order.
!
!    Input, real ( kind = rk ) XVAL, a value to be bracketed.
!
!    Output, integer LEFT, RIGHT, the results of the search.
!    Either:
!      XVAL < X(1), when LEFT = 1, RIGHT = 2;
!      X(N) < XVAL, when LEFT = N-1, RIGHT = N;
!    or
!      X(LEFT) <= XVAL <= X(RIGHT).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer left
  integer right
  real ( kind = rk ) x(n)
  real ( kind = rk ) xval

  do i = 2, n - 1

    if ( xval < x(i) ) then
      left = i - 1
      right = i
      return
    end if

   end do

  left = n - 1
  right = n

  return
end
function r8vec_cross_product_2d ( v1, v2 )

!*****************************************************************************80
!
!! R8VEC_CROSS_PRODUCT_2D finds the cross product of a pair of vectors in 2D.
!
!  Discussion:
!
!    Strictly speaking, the vectors V1 and V2 should be considered
!    to lie in a 3D space, both having Z coordinate zero.  The cross 
!    product value V3 then represents the standard cross product vector 
!    (0,0,V3).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) V1(2), V2(2), the vectors.
!
!    Output, real ( kind = rk ) R8VEC_CROSS_PRODUCT_2D, the cross product.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r8vec_cross_product_2d
  real ( kind = rk ) v1(2)
  real ( kind = rk ) v2(2)

  r8vec_cross_product_2d = v1(1) * v2(2) - v1(2) * v2(1)

  return
end
function r8vec_cross_product_affine_2d ( v0, v1, v2 )

!*****************************************************************************80
!
!! R8VEC_CROSS_PRODUCT_AFFINE_2D finds the affine cross product in 2D.
!
!  Discussion:
!
!    Strictly speaking, the vectors V1 and V2 should be considered
!    to lie in a 3D space, both having Z coordinate zero.  The cross 
!    product value V3 then represents the standard cross product vector 
!    (0,0,V3).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) V0(2), the base vector.
!
!    Input, real ( kind = rk ) V1(2), V2(2), the vectors.
!
!    Output, real ( kind = rk ) R8VEC_CROSS_PRODUCT_AFFINE_2D, 
!    the cross product (V1-V0) x (V2-V0).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r8vec_cross_product_affine_2d
  real ( kind = rk ) v0(2)
  real ( kind = rk ) v1(2)
  real ( kind = rk ) v2(2)

  r8vec_cross_product_affine_2d = &
      ( v1(1) - v0(1) ) * ( v2(2) - v0(2) ) &
    - ( v2(1) - v0(1) ) * ( v1(2) - v0(2) )

  return
end
subroutine r8vec_cross_product_3d ( v1, v2, v3 )

!*****************************************************************************80
!
!! R8VEC_CROSS_PRODUCT_3D computes the cross product of two vectors in 3D.
!
!  Discussion:
!
!    The cross product in 3D can be regarded as the determinant of the
!    symbolic matrix:
!
!          |  i  j  k |
!      det | x1 y1 z1 |
!          | x2 y2 z2 |
!
!      = ( y1 * z2 - z1 * y2 ) * i
!      + ( z1 * x2 - x1 * z2 ) * j
!      + ( x1 * y2 - y1 * x2 ) * k
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) V1(3), V2(3), the two vectors.
!
!    Output, real ( kind = rk ) V3(3), the cross product vector.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) v1(3)
  real ( kind = rk ) v2(3)
  real ( kind = rk ) v3(3)

  v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
  v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
  v3(3) = v1(1) * v2(2) - v1(2) * v2(1)

  return
end
subroutine r8vec_cross_product_affine_3d ( v0, v1, v2, v3 )

!*****************************************************************************80
!
!! R8VEC_CROSS_PRODUCT_AFFINE_3D computes the affine cross product in 3D.
!
!  Discussion:
!
!    The cross product in 3D can be regarded as the determinant of the
!    symbolic matrix:
!
!          |  i  j  k |
!      det | x1 y1 z1 |
!          | x2 y2 z2 |
!
!      = ( y1 * z2 - z1 * y2 ) * i
!      + ( z1 * x2 - x1 * z2 ) * j
!      + ( x1 * y2 - y1 * x2 ) * k
!
!    Here, we use V0 as the base of an affine system so we compute
!    the cross product of (V1-V0) and (V2-V0).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) V0(3), the base vector.
!
!    Input, real ( kind = rk ) V1(3), V2(3), the two vectors.
!
!    Output, real ( kind = rk ) V3(3), the cross product vector
!    ( V1-V0) x (V2-V0).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) v0(3)
  real ( kind = rk ) v1(3)
  real ( kind = rk ) v2(3)
  real ( kind = rk ) v3(3)

  v3(1) = ( v1(2) - v0(2) ) * ( v2(3) - v0(3) ) &
        - ( v2(2) - v0(2) ) * ( v1(3) - v0(3) )

  v3(2) = ( v1(3) - v0(3) ) * ( v2(1) - v0(1) ) &
        - ( v2(3) - v0(3) ) * ( v1(1) - v0(1) )

  v3(3) = ( v1(1) - v0(1) ) * ( v2(2) - v0(2) ) &
        - ( v2(1) - v0(1) ) * ( v1(2) - v0(2) )

  return
end
function r8vec_distance ( dim_num, v1, v2 )

!*****************************************************************************80
!
!! R8VEC_DISTANCE returns the Euclidean distance between two vectors.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = rk ) V1(DIM_NUM), V2(DIM_NUM), the vectors.
!
!    Output, real ( kind = rk ) R8VEC_DISTANCE, the Euclidean distance 
!    between the vectors.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num

  real ( kind = rk ) r8vec_distance
  real ( kind = rk ) v1(dim_num)
  real ( kind = rk ) v2(dim_num)

  r8vec_distance = sqrt ( sum ( ( v1(1:dim_num) - v2(1:dim_num) )**2 ) )

  return
end
function r8vec_dot_product ( dim_num, v1, v2 )

!*****************************************************************************80
!
!! R8VEC_DOT_PRODUCT finds the dot product of a pair of vectors in ND.
!
!  Discussion:
!
!    In FORTRAN, the system routine DOT_PRODUCT should be called
!    directly.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = rk ) V1(DIM_NUM), V2(DIM_NUM), the vectors.
!
!    Output, real ( kind = rk ) R8VEC_DOT_PRODUCT, the dot product.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num

  real ( kind = rk ) r8vec_dot_product
  real ( kind = rk ) v1(dim_num)
  real ( kind = rk ) v2(dim_num)

  r8vec_dot_product = dot_product ( v1(1:dim_num), v2(1:dim_num) )

  return
end
function r8vec_dot_product_affine ( n, v0, v1, v2 )

!*****************************************************************************80
!
!! R8VEC_DOT_PRODUCT_AFFINE computes the affine dot product V1-V0 * V2-V0.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the spatial dimension.
!
!    Input, real ( kind = rk ) V0(N), the base vector.
!
!    Input, real ( kind = rk ) V1(N), V2(N), the vectors.
!
!    Output, real ( kind = rk ) R8VEC_DOT_PRODUCT_AFFINE, the dot product.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) r8vec_dot_product_affine
  real ( kind = rk ) v0(n)
  real ( kind = rk ) v1(n)
  real ( kind = rk ) v2(n)

  r8vec_dot_product_affine = dot_product ( &
    v1(1:n) - v0(1:n),  &
    v2(1:n) - v0(1:n) )

  return
end
function r8vec_eq ( n, a1, a2 )

!*****************************************************************************80
!
!! R8VEC_EQ is true if every pair of entries in two vectors is equal.
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
!  Parameters:
!
!    Input, integer N, the number of entries in the vectors.
!
!    Input, real ( kind = rk ) A1(N), A2(N), two vectors to compare.
!
!    Output, logical R8VEC_EQ.
!    R8VEC_EQ is TRUE if every pair of elements A1(I) and A2(I) are equal.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a1(n)
  real ( kind = rk ) a2(n)
  logical r8vec_eq

  r8vec_eq = ( all ( a1(1:n) == a2(1:n) ) )

  return
end
function r8vec_gt ( n, a1, a2 )

!*****************************************************************************80
!
!! R8VEC_GT == ( A1 > A2 ) for real vectors.
!
!  Discussion:
!
!    The comparison is lexicographic.
!
!    A1 > A2  <=>                              A1(1) > A2(1) or
!                 ( A1(1)     == A2(1)     and A1(2) > A2(2) ) or
!                 ...
!                 ( A1(1:N-1) == A2(1:N-1) and A1(N) > A2(N)
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of the vectors.
!
!    Input, real ( kind = rk ) A1(N), A2(N), the vectors to be compared.
!
!    Output, logical R8VEC_GT, is TRUE if and only if A1 > A2.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a1(n)
  real ( kind = rk ) a2(n)
  integer i
  logical r8vec_gt

  r8vec_gt = .false.

  do i = 1, n

    if ( a2(i) < a1(i) ) then
      r8vec_gt = .true.
      exit
    else if ( a1(i) < a2(i) ) then
      r8vec_gt = .false.
      exit
    end if

  end do

  return
end
function r8vec_lt ( n, a1, a2 )

!*****************************************************************************80
!
!! R8VEC_LT == ( A1 < A2 ) for real vectors.
!
!  Discussion:
!
!    The comparison is lexicographic.
!
!    A1 < A2  <=>                              A1(1) < A2(1) or
!                 ( A1(1)     == A2(1)     and A1(2) < A2(2) ) or
!                 ...
!                 ( A1(1:N-1) == A2(1:N-1) and A1(N) < A2(N)
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of the vectors.
!
!    Input, real ( kind = rk ) A1(N), A2(N), the vectors to be compared.
!
!    Output, logical R8VEC_LT, is TRUE if and only if A1 < A2.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a1(n)
  real ( kind = rk ) a2(n)
  integer i
  logical r8vec_lt

  r8vec_lt = .false.

  do i = 1, n

    if ( a1(i) < a2(i) ) then
      r8vec_lt = .true.
      exit
    else if ( a2(i) < a1(i) ) then
      r8vec_lt = .false.
      exit
    end if

  end do

  return
end
function r8vec_norm ( n, a )

!*****************************************************************************80
!
!! R8VEC_NORM returns the L2 norm of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The vector L2 norm is defined as:
!
!      R8VEC_NORM = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in A.
!
!    Input, real ( kind = rk ) A(N), the vector whose L2 norm is desired.
!
!    Output, real ( kind = rk ) R8VEC_NORM, the L2 norm of A.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) r8vec_norm

  r8vec_norm = sqrt ( sum ( a(1:n)**2 ) )

  return
end
function r8vec_norm_affine ( n, v0, v1 )

!*****************************************************************************80
!
!! R8VEC_NORM_AFFINE returns the affine norm of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The affine vector L2 norm is defined as:
!
!      R8VEC_NORM_AFFINE(V0,V1) 
!        = sqrt ( sum ( 1 <= I <= N ) ( V1(I) - V0(I) )^2 )
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the vectors.
!
!    Input, real ( kind = rk ) V0(N), the base vector.
!
!    Input, real ( kind = rk ) V1(N), the vector whose affine norm is desired.
!
!    Output, real ( kind = rk ) R8VEC_NORM_AFFINE, the L2 norm of V1-V0.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) r8vec_norm_affine
  real ( kind = rk ) v0(n)
  real ( kind = rk ) v1(n)

  r8vec_norm_affine = sqrt ( sum ( ( v0(1:n) - v1(1:n) )**2 ) )

  return
end
subroutine r8vec_normal_01 ( n, x )

!*****************************************************************************80
!
!! R8VEC_NORMAL_01 samples the unit normal probability distribution.
!
!  Discussion:
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!    This routine can generate a vector of values on one call.  It
!    has the feature that it should provide the same results
!    in the same order no matter how we break up the task.
!
!    The Box-Muller method is used, which is efficient, but
!    generates an even number of values each time.  On any call
!    to this routine, an even number of new values are generated.
!    Depending on the situation, one value may be left over.
!    In that case, it is saved for the next call.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 January 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of values desired.  If N is
!    negative, then the code will flush its internal memory; in particular,
!    if there is a saved value to be used on the next call, it is
!    instead discarded. 
!
!    Output, real ( kind = rk ) X(N), a sample of the standard normal PDF.
!
!  Local:
!
!    Local, integer MADE, records the number of values that have
!    been computed.  On input with negative N, this value overwrites
!    the return value of N, so the user can get an accounting of
!    how much work has been done.
!
!    Local, real ( kind = rk ) R(N+1), is used to store some uniform
!    random values.  Its dimension is N+1, but really it is only needed
!    to be the smallest even number greater than or equal to N.
!
!    Local, integer SAVED, is 0 or 1 depending on whether there is a
!    single saved value left over from the previous call.
!
!    Local, integer X_LO_INDEX, X_HI_INDEX, records the range of entries of
!    X that we need to compute.  This starts off as 1:N, but is adjusted
!    if we have a saved value that can be immediately stored in X(1),
!    and so on.
!
!    Local, real ( kind = rk ) Y, the value saved from the previous call, if
!    SAVED is 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer m
  integer, save :: made = 0
  real ( kind = rk ) r(n+1)
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  integer, save :: saved = 0
  real ( kind = rk ) x(n)
  integer x_hi_index
  integer x_lo_index
  real ( kind = rk ), save :: y = 0.0D+00
!
!  I'd like to allow the user to reset the internal data.
!  But this won't work properly if we have a saved value Y.
!  I'm making a crock option that allows the user to signal
!  explicitly that any internal memory should be flushed,
!  by passing in a negative value for N.
!
  if ( n < 0 ) then
    n = made
    made = 0
    saved = 0
    y = 0.0D+00
    return
  else if ( n == 0 ) then
    return
  end if
!
!  Record the range of X we need to fill in.
!
  x_lo_index = 1
  x_hi_index = n
!
!  Use up the old value, if we have it.
!
  if ( saved == 1 ) then
    x(1) = y
    saved = 0
    x_lo_index = 2
  end if
!
!  Maybe we don't need any more values.
!
  if ( x_hi_index - x_lo_index + 1 == 0 ) then
!
!  If we need just one new value, do that here to avoid null arrays.
!
  else if ( x_hi_index - x_lo_index + 1 == 1 ) then

    call random_number ( harvest = r(1:2) )

    x(x_hi_index) = &
             sqrt ( -2.0D+00 * log ( r(1) ) ) * cos ( 2.0D+00 * r8_pi * r(2) )
    y =      sqrt ( -2.0D+00 * log ( r(1) ) ) * sin ( 2.0D+00 * r8_pi * r(2) )

    saved = 1

    made = made + 2
!
!  If we require an even number of values, that's easy.
!
  else if ( mod ( x_hi_index - x_lo_index + 1, 2 ) == 0 ) then

    m = ( x_hi_index - x_lo_index + 1 ) / 2

    call random_number ( harvest = r(1:2*m) )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * cos ( 2.0D+00 * r8_pi * r(2:2*m:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * sin ( 2.0D+00 * r8_pi * r(2:2*m:2) )

    made = made + x_hi_index - x_lo_index + 1
!
!  If we require an odd number of values, we generate an even number,
!  and handle the last pair specially, storing one in X(N), and
!  saving the other for later.
!
  else

    x_hi_index = x_hi_index - 1

    m = ( x_hi_index - x_lo_index + 1 ) / 2 + 1

    call random_number ( harvest = r(1:2*m) )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * cos ( 2.0D+00 * r8_pi * r(2:2*m-2:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * sin ( 2.0D+00 * r8_pi * r(2:2*m-2:2) )

    x(n) = sqrt ( -2.0E+00 * log ( r(2*m-1) ) ) &
      * cos ( 2.0D+00 * r8_pi * r(2*m) )

    y = sqrt ( -2.0D+00 * log ( r(2*m-1) ) ) &
      * sin ( 2.0D+00 * r8_pi * r(2*m) )

    saved = 1

    made = made + x_hi_index - x_lo_index + 2

  end if

  return
end
function r8vec_normsq ( n, v )

!*****************************************************************************80
!
!! R8VEC_NORMSQ returns the square of the L2 norm of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The square of the vector L2 norm is defined as:
!
!      R8VEC_NORMSQ = sum ( 1 <= I <= N ) V(I)^2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    28 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the vector dimension.
!
!    Input, real ( kind = rk ) V(N), the vector.
!
!    Output, real ( kind = rk ) R8VEC_NORMSQ, the squared L2 norm.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) r8vec_normsq
  real ( kind = rk ) v(n)

  r8vec_normsq = sum ( v(1:n)**2 )

  return
end
function r8vec_normsq_affine ( n, v0, v1 )

!*****************************************************************************80
!
!! R8VEC_NORMSQ_AFFINE returns the affine squared norm of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The affine squared vector L2 norm is defined as:
!
!      R8VEC_NORMSQ_AFFINE(V0,V1) 
!        = sum ( 1 <= I <= N ) ( V1(I) - V0(I) )^2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    28 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the vector dimension.
!
!    Input, real ( kind = rk ) V0(N), the base vector.
!
!    Input, real ( kind = rk ) V1(N), the vector.
!
!    Output, real ( kind = rk ) R8VEC_NORMSQ_AFFINE, the affine squared L2 norm.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) r8vec_normsq_affine
  real ( kind = rk ) v0(n)
  real ( kind = rk ) v1(n)

  r8vec_normsq_affine = sum ( ( v0(1:n) - v1(1:n) )**2 )

  return
end
subroutine r8vec_polarize ( n, a, p, a_normal, a_parallel )

!*****************************************************************************80
!
!! R8VEC_POLARIZE decomposes an R8VEC into normal and parallel components.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The (nonzero) vector P defines a direction.
!
!    The vector A can be written as the sum
!
!      A = A_normal + A_parallel
!
!    where A_parallel is a linear multiple of P, and A_normal
!    is perpendicular to P.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    08 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input, real ( kind = rk ) A(N), the vector to be polarized.
!
!    Input, real ( kind = rk ) P(N), the polarizing direction.
!
!    Output, real ( kind = rk ) A_NORMAL(N), A_PARALLEL(N), the normal
!    and parallel components of A.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) a_dot_p
  real ( kind = rk ) a_normal(n)
  real ( kind = rk ) a_parallel(n)
  real ( kind = rk ) p(n)
  real ( kind = rk ) p_norm

  p_norm = sqrt ( sum ( p(1:n)**2 ) )

  if ( p_norm == 0.0D+00 ) then
    a_normal(1:n) = a(1:n)
    a_parallel(1:n) = 0.0D+00
    return
  end if

  a_dot_p = dot_product ( a(1:n), p(1:n) ) / p_norm

  a_parallel(1:n) = a_dot_p * p(1:n) / p_norm

  a_normal(1:n) = a(1:n) - a_parallel(1:n)

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
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
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, real ( kind = rk ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  integer i
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
  end do

  return
end
function r8vec_scalar_triple_product ( v1, v2, v3 )

!*****************************************************************************80
!
!! R8VEC_SCALAR_TRIPLE_PRODUCT finds the scalar triple product in 3D.
!
!  Discussion:
!
!    [A,B,C] = A dot ( B cross C ) 
!            = B dot ( C cross A )
!            = C dot ( A cross B )
!
!    The volume of a parallelepiped, whose sides are given by
!    vectors A, B, and C, is abs ( A dot ( B cross C ) ).
!
!    Three vectors are coplanar if and only if their scalr triple 
!    product vanishes.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Eric Weisstein,
!    "Scalar Triple Product",
!    CRC Concise Encyclopedia of Mathematics, 
!    CRC, 1999
!
!  Parameters:
!
!    Input, real ( kind = rk ) V1(3), V2(3), V3(3), the vectors.
!
!    Output, real ( kind = rk ) R8VEC_SCALAR_TRIPLE_PRODUCT, the scalar
!    triple product.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r8vec_scalar_triple_product
  real ( kind = rk ) v1(3)
  real ( kind = rk ) v2(3)
  real ( kind = rk ) v3(3)
  real ( kind = rk ) v4(3)

  call r8vec_cross_product_3d ( v2, v3, v4 )

  r8vec_scalar_triple_product = dot_product ( v1(1:3), v4(1:3) )

  return
end
subroutine r8vec_swap ( n, a1, a2 )

!*****************************************************************************80
!
!! R8VEC_SWAP swaps the entries of two real vectors.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the arrays.
!
!    Input/output, real ( kind = rk ) A1(N), A2(N), the vectors to swap.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a1(n)
  real ( kind = rk ) a2(n)
  real ( kind = rk ) a3(n)

  a3(1:n) = a1(1:n)
  a1(1:n) = a2(1:n)
  a2(1:n) = a3(1:n)

  return
end
subroutine r8vec_uniform_ab ( n, a, b, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_AB returns a scaled pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of real ( kind = rk ) values.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer M, the number of entries in the vector.
!
!    Input, real ( kind = rk ) A, B, the lower and upper limits.
!
!    Output, real ( kind = rk ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) r(n)

  call random_number ( harvest = r(1:n) )
  r(1:n) = a + ( b - a ) * r(1:n)

  return
end
subroutine r8vec_uniform_unit ( m, w )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_UNIT generates a uniformly random unit vector.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the spatial dimension.
!
!    Output, real ( kind = rk ) W(M), a random direction vector,
!    with unit norm.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m

  real ( kind = rk ) norm
  real ( kind = rk ) w(m)
!
!  Get M values from a standard normal distribution.
!
  call r8vec_normal_01 ( m, w )
!
!  Compute the length of the vector.
!
  norm = sqrt ( sum ( w(1:m)**2 ) )
!
!  Normalize the vector.
!
  w(1:m) = w(1:m) / norm

  return
end
subroutine radec_distance_3d ( ra1, dec1, ra2, dec2, theta )

!*****************************************************************************80
!
!! RADEC_DISTANCE_3D - angular distance, astronomical units, sphere in 3D.
!
!  Discussion:
!
!    Right ascension is measured in hours, between 0 and 24, and
!    essentially measures longitude.
!
!    Declination measures the angle from the equator towards the north pole,
!    and ranges from -90 (South Pole) to 90 (North Pole).
!
!    On the unit sphere, the angular separation between two points is 
!    equal to their geodesic or great circle distance.  On any other
!    sphere, multiply the angular separation by the radius of the
!    sphere to get the geodesic or great circle distance.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) RA1, DEC1, RA2, DEC2, the right ascension and
!    declination of the two points.
!
!    Output, real ( kind = rk ) THETA, the angular separation between the points,
!    in radians.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) cos_theta
  real ( kind = rk ) dec1
  real ( kind = rk ) dec2
  real ( kind = rk ) degrees_to_radians
  real ( kind = rk ) norm_v1
  real ( kind = rk ) norm_v2
  real ( kind = rk ) phi1
  real ( kind = rk ) phi2
  real ( kind = rk ) r8_acos
  real ( kind = rk ) ra1
  real ( kind = rk ) ra2
  real ( kind = rk ) theta
  real ( kind = rk ) theta1
  real ( kind = rk ) theta2
  real ( kind = rk ) v1(dim_num)
  real ( kind = rk ) v2(dim_num)

  theta1 = degrees_to_radians ( 15.0D+00 * ra1 )
  phi1 = degrees_to_radians ( dec1 )

  v1(1:dim_num) = (/ cos ( theta1 ) * cos ( phi1 ), &
                     sin ( theta1 ) * cos ( phi1 ), &
                                      sin ( phi1 ) /)

  norm_v1 = sqrt ( sum ( v1(1:dim_num)**2 ) )

  theta2 = degrees_to_radians ( 15.0D+00 * ra2 )
  phi2 = degrees_to_radians ( dec2 )

  v2(1:dim_num) = (/ cos ( theta2 ) * cos ( phi2 ), &
                     sin ( theta2 ) * cos ( phi2 ), &
                                      sin ( phi2 ) /)

  norm_v2 = sqrt ( sum ( v2(1:dim_num)**2 ) )

  cos_theta = dot_product ( v1(1:dim_num), v2(1:dim_num) ) &
    / ( norm_v1 * norm_v2 )

  theta = r8_acos ( cos_theta )

  return
end
subroutine radec_to_xyz ( ra, dec, p )

!*****************************************************************************80
!
!! RADEC_TO_XYZ converts right ascension/declination to (X,Y,Z) coordinates.
!
!  Discussion:
!
!    Right ascension is measured in hours, between 0 and 24, and
!    essentially measures longitude.
!
!    Declination measures the angle from the equator towards the north pole,
!    and ranges from -90 (South Pole) to 90 (North Pole).
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
!  Parameters:
!
!    Input, real ( kind = rk ) RA, DEC, the right ascension and declination
!    of a point.
!
!    Output, real ( kind = rk ) P(3), the corresponding coordinates of
!    a point with radius 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) dec
  real ( kind = rk ) degrees_to_radians
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) phi
  real ( kind = rk ) ra
  real ( kind = rk ) theta

  theta = degrees_to_radians ( 15.0D+00 * ra )
  phi = degrees_to_radians ( dec )

  p(1) = cos ( theta ) * cos ( phi )
  p(2) = sin ( theta ) * cos ( phi )
  p(3) = sin ( phi )

  return
end
function radians_to_degrees ( angle_rad )

!*****************************************************************************80
!
!! RADIANS_TO_DEGREES converts an angle from radians to degrees.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) ANGLE_RAD, an angle in radians.
!
!    Output, real ( kind = rk ) RADIANS_TO_DEGREES, the equivalent angle
!    in degrees.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) angle_rad
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) radians_to_degrees

  radians_to_degrees = ( angle_rad / r8_pi ) * 180.0D+00

  return
end
subroutine radians_to_dms ( angle_rad, degrees, minutes, seconds )

!*****************************************************************************80
!
!! RADIANS_TO_DMS converts an angle from radians to degrees/minutes/seconds.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) ANGLE_RAD, the angle in radians.
!
!    Output, integer DEGREES, MINUTES, SECONDS, the equivalent 
!    angle in degrees, minutes, and seconds.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) angle_degrees
  real ( kind = rk ) angle_rad
  integer degrees
  integer minutes
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  integer seconds

  angle_degrees = 180.0D+00 * abs ( angle_rad ) / r8_pi

  degrees = int ( angle_degrees )
  angle_degrees = ( angle_degrees - real ( degrees, kind = rk ) ) * 60.0D+00
  minutes = int ( angle_degrees )
  angle_degrees = ( angle_degrees - real ( minutes, kind = rk ) ) * 60.0D+00
  seconds = nint ( angle_degrees )

  if ( angle_rad < 0.0D+00 ) then
    degrees = - degrees
    minutes = - minutes
    seconds = - seconds
  end if

  return
end
subroutine rotation_axis_vector_3d ( axis, angle, v, w )

!*****************************************************************************80
!
!! ROTATION_AXIS_VECTOR_3D rotates a vector around an axis vector in 3D.
!
!  Discussion:
!
!    Thanks to Cody Farnell for correcting some mistakes in an earlier
!    version of this routine.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) AXIS(3), the axis vector for the rotation.
!
!    Input, real ( kind = rk ) ANGLE, the angle, in radians, of the rotation.
!
!    Input, real ( kind = rk ) V(3), the vector to be rotated.
!
!    Output, real ( kind = rk ) W(3), the rotated vector.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) angle
  real ( kind = rk ) axis(dim_num)
  real ( kind = rk ) axis_norm
  real ( kind = rk ) dot
  real ( kind = rk ) norm
  real ( kind = rk ) normal(dim_num)
  real ( kind = rk ) normal_component
  real ( kind = rk ) normal2(dim_num)
  real ( kind = rk ) parallel(dim_num)
  real ( kind = rk ) rot(dim_num)
  real ( kind = rk ) u(dim_num)
  real ( kind = rk ) v(dim_num)
  real ( kind = rk ) w(dim_num)
!
!  Compute the length of the rotation axis.
!
  u(1:dim_num) = axis(1:dim_num)

  axis_norm = sqrt ( sum ( u(1:dim_num)**2 ) )

  if ( axis_norm == 0.0D+00 ) then
    w(1:dim_num) = 0.0D+00
    return
  end if

  u(1:dim_num) = u(1:dim_num) / axis_norm
!
!  Compute the dot product of the vector and the unit rotation axis.
!
  dot = dot_product ( u(1:dim_num), v(1:dim_num) )
!
!  Compute the parallel component of the vector.
!
  parallel(1:dim_num) = dot * u(1:dim_num)
!
!  Compute the normal component of the vector.
!
  normal(1:dim_num) = v(1:dim_num) - parallel(1:dim_num)

  normal_component = sqrt ( sum ( normal(1:dim_num)**2 ) )

  if ( normal_component == 0.0D+00 ) then
    w(1:dim_num) = parallel(1:dim_num)
    return
  end if

  normal(1:dim_num) = normal(1:dim_num) / normal_component
!
!  Compute a second vector, lying in the plane, perpendicular
!  to V, and forming a right-handed system, as the cross product
!  of the first two vectors.
!
  normal2(1) = u(2) * normal(3) - u(3) * normal(2)
  normal2(2) = u(3) * normal(1) - u(1) * normal(3)
  normal2(3) = u(1) * normal(2) - u(2) * normal(1)

  norm = sqrt ( sum ( normal2(1:dim_num)**2 ) )

  normal2(1:dim_num) = normal2(1:dim_num) / norm
!
!  Rotate the normal component by the angle.
!
  rot(1:dim_num) = normal_component * ( &
      cos ( angle ) * normal(1:dim_num) &
    + sin ( angle ) * normal2(1:dim_num) )
!
!  The rotated vector is the parallel component plus the rotated component.
!
  w(1:dim_num) = parallel(1:dim_num) + rot(1:dim_num)

  return
end
subroutine rtp_to_xyz ( r, theta, phi, xyz )

!*****************************************************************************80
!
!! RTP_TO_XYZ converts (R,Theta,Phi) to (X,Y,Z) coordinates.
!
!  Discussion:
!
!    R measures the distance of the point to the origin.
!
!    Theta measures the "longitude" of the point, between 0 and 2 PI.
!
!    PHI measures the angle from the "north pole", between 0 and PI.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) R, THETA, PHI, the radius, longitude, and
!    declination of a point.
!
!    Output, real ( kind = rk ) XYZ(3), the corresponding Cartesian coordinates. 
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) phi
  real ( kind = rk ) r
  real ( kind = rk ) theta
  real ( kind = rk ) xyz(3)

  xyz(1) = r * cos ( theta ) * sin ( phi )
  xyz(2) = r * sin ( theta ) * sin ( phi )
  xyz(3) = r *                 cos ( phi )

  return
end
subroutine segment_contains_point_1d ( p1, p2, p, t )

!*****************************************************************************80
!
!! SEGMENT_CONTAINS_POINT_1D reports if a line segment contains a point in 1D.
!
!  Discussion:
!
!    A line segment is the finite portion of a line that lies between
!    two points P1 and P2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    28 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1, P2, two points defining a line segment.
!    The line segment has T = 0 at P1, and T = 1 at P2.
!
!    Input, real ( kind = rk ) P, a point to be tested.
!
!    Output, real ( kind = rk ) T, the coordinate of P3 in units of (P2-P1).
!    The point P3 is contained in the line segment if 0 <= T <= 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) p
  real ( kind = rk ) p1
  real ( kind = rk ) p2
  real ( kind = rk ) t
  real ( kind = rk ) unit

  unit = p2 - p1

  if ( unit == 0.0D+00 ) then

    if ( p == p1 ) then
      t = 0.5D+00
    else if ( p < p1 ) then
      t = - huge ( t )
    else if ( p1 < p ) then
      t = huge ( t )
    end if

  else

    t = ( p - p1 ) / unit

  end if

  return
end
subroutine segment_contains_point_2d ( p1, p2, p, u )

!*****************************************************************************80
!
!! SEGMENT_CONTAINS_POINT_2D reports if a line segment contains a point in 2D.
!
!  Discussion:
!
!    A line segment is the finite portion of a line that lies between
!    two points P1 and P2.
!
!    In exact arithmetic, point P is on the line segment between
!    P1 and P2 if and only if 0 <= U <= 1 and V = 0.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(2), P2(2), the endpoints of the line segment.
!
!    Input, real ( kind = rk ) P(2), a point to be tested.
!
!    Output, real ( kind = rk ) U(2), the components of P, with the first
!    component measured along the axis with origin at P1 and unit at P2, 
!    and second component the magnitude of the off-axis portion of the
!    vector P-P1, measured in units of (P2-P1).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) normsq
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) u(dim_num)

  normsq = sum ( ( p2(1:dim_num) - p1(1:dim_num) )**2 )

  if ( normsq == 0.0D+00 ) then

    if ( all ( p(1:dim_num) == p1(1:dim_num) ) ) then
      u(1) = 0.5D+00
      u(2) = 0.0D+00
    else
      u(1) = 0.5D+00
      u(2) = huge ( u(2) )
    end if

  else

    u(1) = sum ( ( p(1:dim_num)  - p1(1:dim_num) ) &
               * ( p2(1:dim_num) - p1(1:dim_num) ) ) / normsq

    u(2) = sqrt ( ( ( u(1) - 1.0D+00 ) * p1(1) - u(1) * p2(1) + p(1) )**2 &
                + ( ( u(1) - 1.0D+00 ) * p1(2) - u(1) * p2(2) + p(2) )**2 ) &
                / sqrt ( normsq )

  end if

  return
end
subroutine segment_point_coords_2d ( p1, p2, p, s, t )

!*****************************************************************************80
!
!! SEGMENT_POINT_COORDS_2D: coordinates of a point on a line segment in 2D.
!
!  Discussion:
!
!    A line segment is the finite portion of a line that lies between
!    two points P1 and P2.
!
!    By the coordinates of a point P with respect to a line segment [P1,P2]
!    we mean numbers S and T such that S gives us the distance from the
!    point P to the nearest point PN on the line (not the line segment!), 
!    and T gives us the position of PN relative to P1 and P2.
!
!    If S is zero, then P lies on the line.
!
!    If 0 <= T <= 1, then PN lies on the line segment.
!
!    If both conditions hold, then P lies on the line segment.
!
!    If E is the length of the line segment, then the distance of the 
!    point to the line segment is:
!
!      sqrt ( S^2 +  T^2    * E^2 )     if      T <= 0;
!             S                         if 0 <= T <= 1
!      sqrt ( S^2 + (T-1)^2 * E^2 )     if 1 <= T
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    28 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(2), P2(2), the endpoints of the line segment.
!
!    Input, real ( kind = rk ) P(2), the point to be considered.
!
!    Output, real ( kind = rk ) S, the distance of P to the nearest point PN
!    on the line through P1 and P2.  (S will always be nonnegative.)
!
!    Output, real ( kind = rk ) T, the relative position of the point PN
!    to the points P1 and P2.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) bot
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) pn(dim_num)
  real ( kind = rk ) s
  real ( kind = rk ) t
!
!  If the line segment is actually a point, then the answer is easy.
!
  if ( all ( p1(1:dim_num) == p2(1:dim_num) ) ) then

    t = 0.0D+00

  else

    bot = sum ( ( p2(1:dim_num) - p1(1:dim_num) )**2 )

    t = sum ( ( p(1:dim_num)  - p1(1:dim_num) ) &
            * ( p2(1:dim_num) - p1(1:dim_num) ) ) / bot

  end if

  pn(1:dim_num) = p1(1:dim_num) + t * ( p2(1:dim_num) - p1(1:dim_num) )

  s = sqrt ( sum ( ( p(1:dim_num) - pn(1:dim_num) )**2 ) )

  return
end
subroutine segment_point_coords_3d ( p1, p2, p, s, t )

!*****************************************************************************80
!
!! SEGMENT_POINT_COORDS_3D: coordinates of a point on a line segment in 3D.
!
!  Discussion:
!
!    A line segment is the finite portion of a line that lies between
!    two points P1 and P2.
!
!    By the coordinates of a point P with respect to a line segment [P1,P2]
!    we mean numbers S and T such that S gives us the distance from the
!    point P to the nearest point PN on the line (not the line segment!), 
!    and T gives us the position of PN relative to P1 and P2.
!
!    If S is zero, then P lies on the line.
!
!    If 0 <= T <= 1, then PN lies on the line segment.
!
!    If both conditions hold, then P lies on the line segment.
!
!    If E is the length of the line segment, then the distance of the 
!    point to the line segment is:
!
!      sqrt ( S^2 +  T^2    * E^2 )     if      T <= 0;
!             S                         if 0 <= T <= 1
!      sqrt ( S^2 + (T-1)^2 * E^2 )     if 1 <= T
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(3), P2(3), the endpoints of the line segment.
!
!    Input, real ( kind = rk ) P(3), the point to be considered.
!
!    Output, real ( kind = rk ) S, the distance of P to the nearest point PN
!    on the line through P1 and P2.  (S will always be nonnegative.)
!
!    Output, real ( kind = rk ) T, the relative position of the point PN
!    to the points P1 and P2.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) bot
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) pn(dim_num)
  real ( kind = rk ) s
  real ( kind = rk ) t
!
!  If the line segment is actually a point, then the answer is easy.
!
  if ( all ( p1(1:dim_num) == p2(1:dim_num) ) ) then

    t = 0.0D+00

  else

    bot = sum ( ( p2(1:dim_num) - p1(1:dim_num) )**2 )

    t = sum ( ( p(1:dim_num)  - p1(1:dim_num) ) &
            * ( p2(1:dim_num) - p1(1:dim_num) ) ) / bot

  end if

  pn(1:dim_num) = p1(1:dim_num) + t * ( p2(1:dim_num) - p1(1:dim_num) )

  s = sqrt ( sum ( ( p(1:dim_num) - pn(1:dim_num) )**2 ) )

  return
end
subroutine segment_point_dist_2d ( p1, p2, p, dist )

!*****************************************************************************80
!
!! SEGMENT_POINT_DIST_2D: distance ( line segment, point ) in 2D.
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
!  Parameters:
!
!    Input, real ( kind = rk ) P1(2), P2(2), the endpoints of the line segment.
!
!    Input, real ( kind = rk ) P(2), the point whose nearest neighbor on the line
!    segment is to be determined.
!
!    Output, real ( kind = rk ) DIST, the distance from the point to the
!    line segment.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) bot
  real ( kind = rk ) dist
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) pn(dim_num)
  real ( kind = rk ) t
!
!  If the line segment is actually a point, then the answer is easy.
!
  if ( all ( p1(1:dim_num) == p2(1:dim_num) ) ) then

    t = 0.0D+00

  else

    bot = sum ( ( p2(1:dim_num) - p1(1:dim_num) )**2 )

    t = sum ( ( p(1:dim_num)  - p1(1:dim_num) ) &
            * ( p2(1:dim_num) - p1(1:dim_num) ) ) / bot

    t = max ( t, 0.0D+00 )
    t = min ( t, 1.0D+00 )

  end if

  pn(1:dim_num) = p1(1:dim_num) + t * ( p2(1:dim_num) - p1(1:dim_num) )

  dist = sqrt ( sum ( ( p(1:dim_num) - pn(1:dim_num) )**2 ) )

  return
end
subroutine segment_point_dist_3d ( p1, p2, p, dist )

!*****************************************************************************80
!
!! SEGMENT_POINT_DIST_3D: distance ( line segment, point ) in 3D.
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
!  Parameters:
!
!    Input, real ( kind = rk ) P1(3), P2(3), the endpoints of the segment.
!
!    Input, real ( kind = rk ) P(3), the point whose nearest neighbor on
!    the line segment is to be determined.
!
!    Output, real ( kind = rk ) DIST, the distance from the point to the
!    line segment.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) bot
  real ( kind = rk ) dist
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) pn(dim_num)
  real ( kind = rk ) t
!
!  If the line segment is actually a point, then the answer is easy.
!
  if ( all ( p1(1:dim_num) == p2(1:dim_num) ) ) then

    t = 0.0D+00

  else

    bot = sum ( ( p2(1:dim_num) - p1(1:dim_num) )**2 )

    t = sum ( ( p(1:dim_num)  - p1(1:dim_num) ) &
            * ( p2(1:dim_num) - p1(1:dim_num) ) ) / bot

    t = max ( t, 0.0D+00 )
    t = min ( t, 1.0D+00 )

  end if

  pn(1:dim_num) = p1(1:dim_num) + t * ( p2(1:dim_num) - p1(1:dim_num) )

  dist = sqrt ( sum ( ( p(1:dim_num) - pn(1:dim_num) )**2 ) )

  return
end
subroutine segment_point_near_2d ( p1, p2, p, pn, dist, t )

!*****************************************************************************80
!
!! SEGMENT_POINT_NEAR_2D: nearest point on line segment to point in 2D.
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
!  Parameters:
!
!    Input, real ( kind = rk ) P1(2), P2(2), the endpoints of the line segment.
!
!    Input, real ( kind = rk ) P(2), the point whose nearest neighbor
!    on the line segment is to be determined.
!
!    Output, real ( kind = rk ) PN(2), the point on the line segment which is
!    nearest the point P.
!
!    Output, real ( kind = rk ) DIST, the distance from the point to the 
!    nearest point on the line segment.
!
!    Output, real ( kind = rk ) T, the relative position of the point PN
!    to the points P1 and P2.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) bot
  real ( kind = rk ) dist
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) pn(dim_num)
  real ( kind = rk ) t
!
!  If the line segment is actually a point, then the answer is easy.
!
  if ( all ( p1(1:dim_num) == p2(1:dim_num) ) ) then

    t = 0.0D+00

  else

    bot = sum ( ( p2(1:dim_num) - p1(1:dim_num) )**2 )

    t = sum ( ( p(1:dim_num)  - p1(1:dim_num) ) &
            * ( p2(1:dim_num) - p1(1:dim_num) ) ) / bot

    t = max ( t, 0.0D+00 )
    t = min ( t, 1.0D+00 )

  end if

  pn(1:dim_num) = p1(1:dim_num) + t * ( p2(1:dim_num) - p1(1:dim_num) )

  dist = sqrt ( sum ( ( p(1:dim_num) - pn(1:dim_num) )**2 ) )

  return
end
subroutine segment_point_near_3d ( p1, p2, p, pn, dist, t )

!*****************************************************************************80
!
!! SEGMENT_POINT_NEAR_3D: nearest point on line segment to point in 3D.
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
!  Parameters:
!
!    Input, real ( kind = rk ) P1(3), P2(3), the endpoints of the segment.
!
!    Input, real ( kind = rk ) P(3), the point whose nearest neighbor
!    on the line segment is to be determined.
!
!    Output, real ( kind = rk ) PN(3), the point on the line segment
!    nearest to P.
!
!    Output, real ( kind = rk ) DIST, the distance from the point to the
!    nearest point on the line segment.
!
!    Output, real ( kind = rk ) T, the relative position of the nearest point
!    P to P1 and P2, that is PN = (1-T)*P1 + T*P2.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) bot
  real ( kind = rk ) dist
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) pn(dim_num)
  real ( kind = rk ) t
!
!  If the line segment is actually a point, then the answer is easy.
!
  if ( all ( p1(1:dim_num) == p2(1:dim_num) ) ) then

    t = 0.0D+00

  else

    bot = sum ( ( p2(1:dim_num) - p1(1:dim_num) )**2 )

    t = sum ( ( p(1:dim_num)  - p1(1:dim_num) ) &
            * ( p2(1:dim_num) - p1(1:dim_num) ) ) / bot

    t = max ( t, 0.0D+00 )
    t = min ( t, 1.0D+00 )

  end if

  pn(1:dim_num) = p1(1:dim_num) + t * ( p2(1:dim_num) - p1(1:dim_num) )

  dist = sqrt ( sum ( ( p(1:dim_num) - pn(1:dim_num) )**2 ) )

  return
end
subroutine segments_curvature_2d ( p1, p2, p3, curvature )

!*****************************************************************************80
!
!! SEGMENTS_CURVATURE_2D computes the curvature of two line segments in 2D.
!
!  Discussion:
!
!    A line segment is the finite portion of a line that lies between
!    two points P1 and P2.
!
!    We assume that the segments are [P1,P2] and [P2,P3].
!
!    We compute the circle that passes through P1, P2 and P3.
!
!    The inverse of the radius of this circle is the local "curvature"
!    associated with the three points.
!
!    If curvature is 0, the two line segments have the same slope,
!    and the three points are collinear.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(2), P2(2), P3(2), the points.
!
!    Output, real ( kind = rk ) CURVATURE, the local curvature.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) curvature
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) p3(dim_num)
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) r

  call circle_exp2imp_2d ( p1, p2, p3, r, pc )

  if ( 0.0D+00 < r ) then
    curvature = 1.0D+00 / r
  else
    curvature = 0.0D+00
  end if

  return
end
subroutine segments_dist_2d ( p1, p2, q1, q2, dist )

!*****************************************************************************80
!
!! SEGMENTS_DIST_2D computes the distance between two line segments in 2D.
!
!  Discussion:
!
!    A line segment is the finite portion of a line that lies between
!    two points P1 and P2.
!
!    If the lines through [P1,P2] and [Q1,Q2] intersect, and both
!    line segments include the point of intersection, then the distance
!    is zero and we are done.
!
!    Therefore, we compute the intersection of the two lines, and
!    find the coordinates of that intersection point on each line.
!    This will tell us if the zero distance case has occurred.
!
!    Otherwise, let PN and QN be points in [P1,P2] and [Q1,Q2] for which 
!    the distance is minimal.  If the lines do not intersect, then it 
!    cannot be the case that both PN and QN are strictly interior to their 
!    line segments, aside from the exceptional singular case when
!    the line segments overlap or are parallel.  Even then, one of PN
!    and QN may be taken to be a segment endpoint.
!
!    Therefore, our second computation finds the minimum of:
!
!      Distance ( P1, [Q1,Q2] );
!      Distance ( P2, [Q1,Q2] );
!      Distance ( Q1, [P1,P2] );
!      Distance ( Q2, [P1,P2] );
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(2), P2(2), the endpoints of the first
!    segment.
!
!    Input, real ( kind = rk ) Q1(2), Q2(2), the endpoints of the second
!    segment.
!
!    Output, real ( kind = rk ) DIST, the distance between the line segments.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) dist
  real ( kind = rk ) dist2
  integer ival
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) q1(dim_num)
  real ( kind = rk ) q2(dim_num)
  real ( kind = rk ) r(dim_num)
  real ( kind = rk ) rps
  real ( kind = rk ) rpt
  real ( kind = rk ) rqs
  real ( kind = rk ) rqt
!
!  Determine whether and where the underlying lines intersect.
!
  call lines_exp_int_2d ( p1, p2, q1, q2, ival, r )
!
!  If there is exactly one intersection point part of both lines, 
!  check that it is part of both line segments.
!
  if ( ival == 1 ) then

    call segment_point_coords_2d ( p1, p2, r, rps, rpt )
    call segment_point_coords_2d ( q1, q2, r, rqs, rqt )

    if ( 0.0D+00 <= rpt .and. rpt <= 1.0D+00 .and. &
         0.0D+00 <= rqt .and. rqt <= 1.0D+00 ) then
      dist = 0.0D+00
      return
    end if

  end if
!
!  If there is no intersection, or the intersection point is
!  not part of both line segments, then an endpoint of one
!  line segment achieves the minimum distance.
!
  call segment_point_dist_2d ( q1, q2, p1, dist2 )
  dist = dist2
  call segment_point_dist_2d ( q1, q2, p2, dist2 )
  dist = min ( dist, dist2 )
  call segment_point_dist_2d ( p1, p2, q1, dist2 )
  dist = min ( dist, dist2 )
  call segment_point_dist_2d ( p1, p2, q2, dist2 )
  dist = min ( dist, dist2 )

  return
end
subroutine segments_dist_3d ( p1, p2, q1, q2, dist )

!*****************************************************************************80
!
!! SEGMENTS_DIST_3D computes the distance between two line segments in 3D.
!
!  Discussion:
!
!
!    NOTE: The special cases for identical and parallel lines have not been
!    worked out yet; those cases are exceptional, and so this code
!    is made available in a slightly unfinished form!
!
!
!    A line segment is the finite portion of a line that lies between
!    two points P1 and P2.
!
!    Given two line segments, consider the underlying lines on which
!    they lie.
!   
!    A) If the lines are identical, then the distance between the line segments
!    is 0, if the segments overlap, or otherwise is attained by the
!    minimum of the distances between each endpoint and the opposing
!    line segment.
!
!    B) If the lines are parallel, then the distance is either the distance
!    between the lines, if the projection of one line segment onto
!    the other overlaps, or otherwise is attained by the
!    minimum of the distances between each endpoint and the opposing
!    line segment.
!
!    C) If the lines are not identical, and not parallel, then there are 
!    unique points PN and QN which are the closest pair of points on the lines.
!    If PN is interior to [P1,P2] and QN is interior to [Q1,Q2],
!    then the distance between the two line segments is the distance
!    between PN and QN.  Otherwise, the nearest distance can be computed
!    by taking the minimum of the distance from each endpoing to the
!    opposing line segment.
!
!    Therefore, our computation first checks whether the lines are
!    identical, parallel, or other, and checks for the special case
!    where the minimum occurs in the interior.
!
!    If that case is ruled out, it computes and returns the minimum of:
!
!      Distance ( P1, [Q1,Q2] );
!      Distance ( P2, [Q1,Q2] );
!      Distance ( Q1, [P1,P2] );
!      Distance ( Q2, [P1,P2] );
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(3), P2(3), the endpoints of the first
!    segment.
!
!    Input, real ( kind = rk ) Q1(3), Q2(3), the endpoints of the second
!    segment.
!
!    Output, real ( kind = rk ) DIST, the distance between the line segments.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) d
  real ( kind = rk ) det
  real ( kind = rk ) dist
  real ( kind = rk ) dist2
  real ( kind = rk ) e
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) pn(dim_num)
  real ( kind = rk ) q1(dim_num)
  real ( kind = rk ) q2(dim_num)
  real ( kind = rk ) qn(dim_num)
  real ( kind = rk ) sn
  real ( kind = rk ) tn
  real ( kind = rk ) u(dim_num)
  real ( kind = rk ) v(dim_num)
  real ( kind = rk ) w0(dim_num)
!
!  The lines are identical.
!  THIS CASE NOT SET UP YET
!
! if ( lines_exp_equal_3d ( p1, p2, q1, q2 ) ) then
! end if
!
!  The lines are not identical, but parallel
!  THIS CASE NOT SET UP YET.
!
! if ( lines_exp_parallel_3d ( p1, p2, q1, q2 ) ) then
! end if
!
!  C: The lines are not identical, not parallel.
!

!
!  Let U = (P2-P1) and V = (Q2-Q1) be the direction vectors on
!  the two lines.
!
  u(1:dim_num) = p2(1:dim_num) - p1(1:dim_num)
  v(1:dim_num) = q2(1:dim_num) - q1(1:dim_num)
!
!  Let SN be the unknown coordinate of the nearest point PN on line 1,
!  so that PN = P(SN) = P1 + SN * (P2-P1).
!
!  Let TN be the unknown coordinate of the nearest point QN on line 2,
!  so that QN = Q(TN) = Q1 + TN * (Q2-Q1).
!
!  Let W0 = (P1-Q1).
!
  w0(1:dim_num) = p1(1:dim_num) - q1(1:dim_num)
!
!  The vector direction WC = P(SN) - Q(TC) is unique (among directions)
!  perpendicular to both U and V, so
!
!    U dot WC = 0
!    V dot WC = 0
!
!  or, equivalently:
!
!    U dot ( P1 + SN * (P2 - P1) - Q1 - TN * (Q2 - Q1) ) = 0
!    V dot ( P1 + SN * (P2 - P1) - Q1 - TN * (Q2 - Q1) ) = 0
!
!  or, equivalently:
!
!    (u dot u ) * sn - (u dot v ) tc = -u * w0
!    (v dot u ) * sn - (v dot v ) tc = -v * w0
!
!  or, equivalently:
!
!   ( a  -b ) * ( sn ) = ( -d )
!   ( b  -c )   ( tc )   ( -e )
!
  a = dot_product ( u, u )
  b = dot_product ( u, v )
  c = dot_product ( v, v )
  d = dot_product ( u, w0 )
  e = dot_product ( v, w0 )
!
!  Check the determinant.
!
  det = - a * c + b * b

  if ( det == 0.0D+00 ) then
    sn = 0.0D+00
    if ( abs ( b ) < abs ( c ) ) then
      tn = e / c
    else
      tn = d / b
    end if
  else
    sn = ( c * d - b * e ) / det
    tn = ( b * d - a * e ) / det
  end if
!
!  Now if both nearest points on the lines
!  also happen to lie inside their line segments,
!  then we have found the nearest points on the line segments.
!
  if ( 0.0D+00 <= sn .and. sn <= 1.0D+00 .and. &
       0.0D+00 <= tn .and. tn <= 1.0D+00 ) then
    pn(1:dim_num) = p1(1:dim_num) + sn * ( p2(1:dim_num) - p1(1:dim_num) )
    qn(1:dim_num) = q1(1:dim_num) + tn * ( q2(1:dim_num) - q1(1:dim_num) )
    dist = sqrt ( sum ( ( pn(1:dim_num) - qn(1:dim_num) )**2 ) )
    return
  end if
!
!  The nearest point did not occur in the interior.
!  Therefore it must be achieved at an endpoint.
!
  call segment_point_dist_3d ( q1, q2, p1, dist2 )
  dist = dist2
  call segment_point_dist_3d ( q1, q2, p2, dist2 )
  dist = min ( dist, dist2 )
  call segment_point_dist_3d ( p1, p2, q1, dist2 )
  dist = min ( dist, dist2 )
  call segment_point_dist_3d ( p1, p2, q2, dist2 )
  dist = min ( dist, dist2 )

  return
end
subroutine segments_dist_3d_old ( p1, p2, q1, q2, dist )

!*****************************************************************************80
!
!! SEGMENTS_DIST_3D_OLD computes the distance between two line segments in 3D.
!
!  Discussion:
!
!    A line segment is the finite portion of a line that lies between
!    two points P1 and P2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(3), P2(3), the endpoints of the
!    first segment.
!
!    Input, real ( kind = rk ) Q1(3), Q2(3), the endpoints of the
!    second segment.
!
!    Output, real ( kind = rk ) DIST, the distance between the line segments.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) d1
  real ( kind = rk ) d2
  real ( kind = rk ) dist
  real ( kind = rk ) dl
  real ( kind = rk ) dm
  real ( kind = rk ) dr
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) pm(dim_num)
  real ( kind = rk ) pn1(dim_num)
  real ( kind = rk ) pn2(dim_num)
  real ( kind = rk ) q1(dim_num)
  real ( kind = rk ) q2(dim_num)
  real ( kind = rk ) t1
  real ( kind = rk ) t2
  real ( kind = rk ) tl
  real ( kind = rk ) tm
  real ( kind = rk ) tmin
  real ( kind = rk ) tr
!
!  Find the nearest points on line 2 to the endpoints of line 1.
!
  call segment_point_near_3d ( q1, q2, p1, pn1, d1, t1 )
  call segment_point_near_3d ( q1, q2, p2, pn2, d2, t2 )

  if ( t1 == t2 ) then
    call segment_point_dist_3d ( p1, p2, pn1, dist )
    return
  end if

  pm(1:dim_num) = 0.5D+00 * ( pn1(1:dim_num) + pn2(1:dim_num) )
!
!  On line 2, over the interval between the points nearest to line 1,
!  the square of the distance of any point to line 1 is a quadratic function.
!  Evaluate it at three points, and seek its local minimum.
!
  call segment_point_dist_3d ( p1, p2, pn1, dl )
  call segment_point_dist_3d ( p1, p2, pm, dm )
  call segment_point_dist_3d ( p1, p2, pn2, dr )

  tl = 0.0D+00
  tm = 0.5D+00
  tr = 1.0D+00

  dl = dl * dl
  dm = dm * dm
  dr = dr * dr

  call minquad ( tl, dl, tm, dm, tr, dr, tmin, dist )

  dist = sqrt ( dist )

  return
end
subroutine segments_int_1d ( p1, p2, q1, q2, dist, r1, r2 )

!*****************************************************************************80
!
!! SEGMENTS_INT_1D computes the intersection of two line segments in 1D.
!
!  Discussion:
!
!    A line segment is the finite portion of a line that lies between
!    two points P1 and P2.
!
!    In 1D, two line segments "intersect" if they overlap.
!
!    Using a real number DIST to report overlap is preferable to 
!    returning a TRUE/FALSE flag, since DIST is better able to 
!    handle cases where the segments "almost" interlap.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1, P2, the endpoints of the first segment.
!
!    Input, real ( kind = rk ) Q1, Q2, the endpoints of the second segment.
!
!    Output, real ( kind = rk ) DIST, the "distance" between the segments.
!    < 0, the segments overlap, and the overlap is DIST units long;
!    = 0, the segments overlap at a single point;
!    > 0, the segments do not overlap.  The distance between the nearest
!    points is DIST units.
!
!    Output, real ( kind = rk ) R1, R2, the endpoints of the intersection
!    segment.  
!    If DIST < 0, then the interval [R1,R2] is the common intersection
!    of the two segments.
!    If DIST = 0, then R1 = R2 is the single common point of the two segments.
!    If DIST > 0, then (R1,R2) is an open interval separating the two
!    segments, which do not overlap at all.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) dist
  real ( kind = rk ) p1
  real ( kind = rk ) p2
  real ( kind = rk ) q1
  real ( kind = rk ) q2
  real ( kind = rk ) r1
  real ( kind = rk ) r2

  r1 = max ( min ( p1, p2 ), &
             min ( q1, q2 ) )

  r2 = min ( max ( p1, p2 ), &
             max ( q1, q2 ) )

  dist = r1 - r2

  return
end
subroutine segments_int_2d ( p1, p2, q1, q2, flag, r )

!*****************************************************************************80
!
!! SEGMENTS_INT_2D computes the intersection of two line segments in 2D.
!
!  Discussion:
!
!    A line segment is the finite portion of a line that lies between
!    two points P1 and P2.
!
!    In 2D, two line segments might not intersect, even though the
!    lines, of which they are portions, intersect.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(2), P2(2), the endpoints of the first
!    segment.
!
!    Input, real ( kind = rk ) Q1(2), Q2(2), the endpoints of the second
!    segment.
!
!    Output, integer FLAG, records the results.
!    0, the line segments do not intersect.
!    1, the line segments intersect.
!
!    Output, real ( kind = rk ) R(2), an intersection point, if there is one.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  integer flag
  integer ival
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) q1(dim_num)
  real ( kind = rk ) q2(dim_num)
  real ( kind = rk ) r(dim_num)
  real ( kind = rk ), parameter :: tol = 0.001D+00
  real ( kind = rk ) u(dim_num)
!
!  Find the intersection of the two lines.
!
  r(1:dim_num) = (/ 0.0D+00, 0.0D+00 /)

  call lines_exp_int_2d ( p1, p2, q1, q2, ival, r )

  if ( ival == 0 ) then
    flag = 0
    return
  end if
!
!  Is the intersection point part of the first line segment?
!
  call segment_contains_point_2d ( p1, p2, r, u )

  if ( u(1) < 0.0D+00 .or. 1.0D+00 < u(1) .or. tol < u(2) ) then
    flag = 0
    return
  end if
!
!  Is the intersection point part of the second line segment?
!
  call segment_contains_point_2d ( q1, q2, r, u )

  if ( u(1) < 0.0D+00 .or. 1.0D+00 < u(1) .or. tol < u(2) ) then
    flag = 0
    return
  end if

  flag = 1

  return
end
subroutine shape_point_dist_2d ( pc, p1, side_num, p, dist )

!*****************************************************************************80
!
!! SHAPE_POINT_DIST_2D: distance ( regular shape, point ) in 2D.
!
!  Discussion:
!
!    The "regular shape" is assumed to be an equilateral and equiangular
!    polygon, such as the standard square, pentagon, hexagon, and so on.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) PC(2), the center of the shape.
!
!    Input, real ( kind = rk ) P1(2), the first vertex of the shape.
!
!    Input, integer SIDE_NUM, the number of sides.
!
!    Input, real ( kind = rk ) P(2), the point to be checked.
!
!    Output, real ( kind = rk ) DIST, the distance from the point to the shape.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) angle
  real ( kind = rk ) angle_degrees_2d
  real ( kind = rk ) angle2
  real ( kind = rk ) degrees_to_radians
  real ( kind = rk ) dist
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) pa(dim_num)
  real ( kind = rk ) pb(dim_num)
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) radius
  real ( kind = rk ) sector_angle
  integer sector_index
  integer side_num
!
!  Determine the angle subtended by a single side.
!
  sector_angle = 360.0D+00 / real ( side_num, kind = rk )
!
!  How long is the half-diagonal?
!
  radius = sqrt ( sum ( ( p1(1:dim_num) - pc(1:dim_num) )**2 ) )
!
!  If the radius is zero, then the shape is a point and the computation is easy.
!
  if ( radius == 0.0D+00 ) then
    dist = sqrt ( sum ( ( p(1:dim_num) - pc(1:dim_num) )**2 ) )
    return
  end if
!
!  If the test point is at the pc, then the computation is easy.
!  The angle subtended by any side is ( 2 * PI / SIDE_NUM ) and the
!  nearest distance is the midpoint of any such side.
!
  if ( all ( p(1:dim_num) == pc(1:dim_num) ) ) then
    dist = radius * cos ( r8_pi / real ( side_num, kind = rk ) )
    return
  end if
!
!  Determine the angle between the ray to the first corner,
!  and the ray to the test point.
!
  angle = angle_degrees_2d ( p1(1:2), pc(1:2), p(1:2) )
!
!  Determine the sector of the point.
!
  sector_index = int ( angle / sector_angle ) + 1
!
!  Generate the two corner points that terminate the SECTOR-th side.
!
  angle2 = real ( sector_index - 1, kind = rk ) * sector_angle
  angle2 = degrees_to_radians ( angle2 )

  call vector_rotate_base_2d ( p1, pc, angle2, pa )

  angle2 = real ( sector_index, kind = rk ) * sector_angle
  angle2 = degrees_to_radians ( angle2 )

  call vector_rotate_base_2d ( p1, pc, angle2, pb )
!
!  Determine the distance from the test point to the line segment that
!  is the SECTOR-th side.
!
  call segment_point_dist_2d ( pa, pb, p, dist )

  return
end
subroutine shape_point_near_2d ( pc, p1, side_num, p, pn, dist )

!*****************************************************************************80
!
!! SHAPE_POINT_NEAR_2D: nearest point ( regular shape, point ) in 2D.
!
!  Discussion:
!
!    The "regular shape" is assumed to be an equilateral and equiangular
!    polygon, such as the standard square, pentagon, hexagon, and so on.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) PC(2), the center of the shape.
!
!    Input, real ( kind = rk ) P1(2), the first vertex of the shape.
!
!    Input, integer SIDE_NUM, the number of sides.
!
!    Input, real ( kind = rk ) P(2), the point to be checked.
!
!    Output, real ( kind = rk ) PN(2), the point on the shape that is nearest
!    to the given point.
!
!    Output, real ( kind = rk ) DIST, the distance between the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) angle
  real ( kind = rk ) angle_degrees_2d
  real ( kind = rk ) angle2
  real ( kind = rk ) degrees_to_radians
  real ( kind = rk ) dist
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) pa(dim_num)
  real ( kind = rk ) pb(dim_num)
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) pd(dim_num)
  real ( kind = rk ) pn(dim_num)
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) radius
  real ( kind = rk ) sector_angle
  integer sector_index
  integer side_num
  real ( kind = rk ) t
!
!  Determine the angle subtended by a single side.
!
  sector_angle = 360.0D+00 / real ( side_num, kind = rk )
!
!  How long is the half-diagonal?
!
  radius = sqrt ( sum ( ( p1(1:dim_num) - pc(1:dim_num) )**2 ) )
!
!  If the radius is zero, then the shape is a point and the computation is easy.
!
  if ( radius == 0.0D+00 ) then
    pn(1:dim_num) = pc(1:dim_num)
    dist = sqrt ( sum ( ( p(1:dim_num) - pn(1:dim_num) )**2 ) )
    return
  end if
!
!  If the test point is at the pc, then the computation is easy.
!  The angle subtended by any side is ( 2 * PI / SIDE_NUM ) and the
!  nearest distance is the midpoint of any such side.
!
  if ( all ( p(1:dim_num) == pc(1:dim_num) ) ) then
    angle = r8_pi / real ( side_num, kind = rk )
    pd(1) =   ( p(1) - pc(1) ) * cos ( angle ) &
            + ( p(2) - pc(2) ) * sin ( angle )
    pd(2) = - ( p(1) - pc(1) ) * sin ( angle ) &
            + ( p(2) - pc(2) ) * cos ( angle )
    pn(1) = pc(1) + pd(1) * cos ( angle )
    pn(2) = pc(2) + pd(2) * sin ( angle )
    dist = radius * cos ( angle )
    return
  end if
!
!  Determine the angle between the ray to the first corner,
!  and the ray to the test point.
!
  angle = angle_degrees_2d ( p1(1:2), pc(1:2), p(1:2) )
!
!  Determine the sector of the point.
!
  sector_index = int ( angle / sector_angle ) + 1
!
!  Generate the two corner points that terminate the SECTOR-th side.
!
  angle2 = real ( sector_index - 1, kind = rk ) * sector_angle
  angle2 = degrees_to_radians ( angle2 )

  call vector_rotate_base_2d ( p1, pc, angle2, pa )

  angle2 = real ( sector_index, kind = rk ) * sector_angle
  angle2 = degrees_to_radians ( angle2 )

  call vector_rotate_base_2d ( p1, pc, angle2, pb )
!
!  Determine the point on the SECTOR-th side of the shape which is
!  nearest.
! 
  call segment_point_near_2d ( pa, pb, p, pn, dist, t )

  return
end
subroutine shape_print_3d ( point_num, face_num, face_order_max, &
  point_coord, face_order, face_point )

!*****************************************************************************80
!
!! SHAPE_PRINT_3D prints information about a polyhedron in 3D.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer POINT_NUM, the number of points.
!
!    Input, integer FACE_NUM, the number of faces.
!
!    Input, integer FACE_ORDER_MAX, the number of vertices 
!    per face.
!
!    Input, real ( kind = rk ) POINT_COORD(3,POINT_NUM), the vertices.
!
!    Input, integer FACE_ORDER(FACE_NUM), the number of vertices
!    per face.
!
!    Input, integer FACE_POINT(FACE_ORDER_MAX,FACE_NUM); 
!    FACE_POINT(I,J) contains the index of the I-th point in the J-th face.  
!    The points are listed in the counter clockwise direction defined
!    by the outward normal at the face.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer face_num
  integer face_order_max
  integer, parameter :: dim_num = 3
  integer point_num

  integer face_order(face_num)
  integer face_point(face_order_max,face_num)
  integer i
  real ( kind = rk ) point_coord(dim_num,point_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SHAPE_PRINT_3D'
  write ( *, '(a)' ) '  Information about a polytope.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of vertices is ', point_num
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Vertices:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Index          X               Y               Z'
  write ( *, '(a)' ) ' '

  do i = 1, point_num
    write ( *, '(2x,i8,2x,3f16.8)' ) i, point_coord(1:dim_num,i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of faces is ', face_num
  write ( *, '(a,i8)' ) '  The maximum order of any face is ', face_order_max
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Index     Order         Indices of Nodes in Face'
  write ( *, '(22x,10i8)' ) ( i, i = 1, face_order_max )
  write ( *, '(a)' ) ' '
  
  do i = 1, face_num
    write ( *, '(2x,i8,2x,i8,2x,10i8)' ) i, face_order(i), &
     face_point(1:face_order(i),i)
  end do

  return
end
subroutine shape_ray_int_2d ( pc, p1, side_num, pa, pb, pint )

!*****************************************************************************80
!
!! SHAPE_RAY_INT_2D: intersection ( regular shape, ray ) in 2D.
!
!  Discussion:
!
!    The "regular shape" is assumed to be an equilateral and equiangular
!    polygon, such as the standard square, pentagon, hexagon, and so on.
!
!    The origin of the ray is assumed to be inside the shape.  This
!    guarantees that the ray will intersect the shape in exactly one point.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    29 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) PC(2), the center of the shape.
!
!    Input, real ( kind = rk ) P1(2), the first vertex of the shape.
!
!    Input, integer SIDE_NUM, the number of sides.
!
!    Input, real ( kind = rk ) PA(2), the origin of the ray.
!
!    Input, real ( kind = rk ) PB(2), a second point on the ray.
!
!    Output, real ( kind = rk ) PINT(2), the point on the shape intersected
!    by the ray.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) angle2
  real ( kind = rk ) degrees_to_radians
  logical inside
  integer ival
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) pa(dim_num)
  real ( kind = rk ) pb(dim_num)
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) pint(dim_num)
  real ( kind = rk ) radius
  real ( kind = rk ) sector_angle
  integer sector_index
  integer side_num
  real ( kind = rk ) v1(dim_num)
  real ( kind = rk ) v2(dim_num)
!
!  Warning!
!  No check is made to ensure that the ray origin is inside the shape.
!  These calculations are not valid if that is not true!
!
!  Determine the angle subtended by a single side.
!
  sector_angle = 360.0D+00 / real ( side_num, kind = rk )
!
!  How long is the half-diagonal?
!
  radius = sqrt ( sum ( ( p1(1:dim_num) - pc(1:dim_num) )**2 ) )
!
!  If the radius is zero, refuse to continue.
!
  if ( radius == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SHAPE_RAY_INT_2D - Fatal error!'
    write ( *, '(a)' ) '  The shape has radius zero.'
    stop 1
  end if
!
!  Determine which sector side intersects the ray.
!
  v2(1:dim_num) = (/ 0.0D+00, 0.0D+00 /)

  do sector_index = 1, side_num
!
!  Determine the two vertices that define this sector.
!
    if ( sector_index == 1 ) then

      angle2 = real ( sector_index - 1, kind = rk ) * sector_angle
      angle2 = degrees_to_radians ( angle2 )

      call vector_rotate_base_2d ( p1, pc, angle2, v1 )

    else

      v1(1:dim_num) = v2(1:dim_num)

    end if

    angle2 = real ( sector_index, kind = rk ) * sector_angle
    angle2 = degrees_to_radians ( angle2 )

    call vector_rotate_base_2d ( p1, pc, angle2, v2 )
!
!  Draw the angle from one vertex to the ray origin to the next vertex,
!  and see if that angle contains the ray.  If so, then the ray
!  must intersect the shape side of that sector.
!   
    call angle_contains_point_2d ( v1, pa, v2, pb, inside )
!
!  Determine the intersection of the lines defined by the ray and the
!  sector side.  (We're already convinced that the ray and sector line
!  segment intersect, so we can use the simpler code that treats them
!  as full lines).
!
    if ( inside ) then

      call lines_exp_int_2d ( pa, pb, v1, v2, ival, pint )

      return

    end if

  end do
!
!  If the calculation fell through the loop, then something's wrong.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SHAPE_RAY_INT_2D - Fatal error!'
  write ( *, '(a)' ) '  Cannot find intersection of ray and shape.'
  stop 1
end
subroutine simplex_lattice_layer_point_next ( n, c, v, more )

!*****************************************************************************80
!
!! SIMPLEX_LATTICE_LAYER_POINT_NEXT: next simplex lattice layer point.
!
!  Discussion:
!
!    The simplex lattice layer L is bounded by the lines
!
!      0 <= X(1:N),
!      L - 1 < sum X(1:N) / C(1:N)  <= L.
!
!    In particular, layer L = 0 always contains just the origin.
!
!    This function returns, one at a time, the points that lie within 
!    a given simplex lattice layer.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    08 July 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the spatial dimension.
!
!    Input, integer C(N+1), coefficients defining the 
!    lattice layer in entries 1 to N, and the laver index in C(N+1).  
!    The coefficients should be positive, and C(N+1) must be nonnegative.
!
!    Input/output, integer V(N).  On first call for a given layer,
!    the input value of V is not important.  On a repeated call for the same
!    layer, the input value of V should be the output value from the previous 
!    call.  On output, V contains the next lattice layer point.
!
!    Input/output, logical MORE.  On input, set MORE to FALSE 
!    to indicate that this is the first call for a given layer.  Thereafter,
!    the input value should be the output value from the previous call.  On 
!    output, MORE is TRUE if the returned value V is a new point.
!    If the output value is FALSE, then no more points were found,
!    and V was reset to 0, and the lattice layer has been exhausted.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer c(n+1)
  integer c1n
  integer i
  integer i4vec_lcm
  integer j
  integer lhs
  logical more
  integer rhs1
  integer rhs2
  integer v(n)
!
!  Treat layer C(N+1) = 0 specially.
!
  if ( c(n+1) == 0 ) then
    if ( .not. more ) then
      v(1:n) = 0
      more = .true.
    else
      more = .false.
    end if
    return
  end if
!
!  Compute the first point.
!
  if ( .not. more ) then

    v(1) = ( c(n+1) - 1 ) * c(1) + 1
    v(2:n) = 0   
    more = .true.

  else

    c1n = i4vec_lcm ( n, c )

    rhs1 = c1n * ( c(n+1) - 1 )
    rhs2 = c1n *   c(n+1)
!
!  Try to increment component I.
!
    do i = 1, n

      v(i) = v(i) + 1

      v(1:i-1) = 0

      if ( 1 < i ) then
        v(1) = rhs1
        do j = 2, n
          v(1) = v(1) - ( c1n / c(j) ) * v(j)
        end do
        v(1) = ( c(1) * v(1) ) / c1n
        v(1) = max ( v(1), 0 )
      end if

      lhs = 0
      do j = 1, n
        lhs = lhs + ( c1n / c(j) ) * v(j)
      end do

      if ( lhs <= rhs1 ) then
        v(1) = v(1) + 1
        lhs = lhs + c1n / c(1)
      end if

      if ( lhs <= rhs2 ) then
        return
      end if

    end do

    v(1:n) = 0
    more = .false.

  end if

  return
end
subroutine simplex_lattice_point_next ( n, c, v, more )

!*****************************************************************************80
!
!! SIMPLEX_LATTICE_POINT_NEXT returns the next simplex lattice point.
!
!  Discussion:
!
!    The lattice simplex is defined by the vertices:
!
!      (0,0,...,0), (C(N+1)/C(1),0,...,0), (0,C(N+1)/C(2),...,0) ... 
!      (0,0,...C(N+1)/C(N))
!
!    The lattice simplex is bounded by the lines
!
!      0 <= V(1:N),
!      V(1) / C(1) + V(2) / C(2) + ... + V(N) / C(N) <= C(N+1)
!
!    Lattice points are listed one at a time, starting at the origin,
!    with V(1) increasing first.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    08 July 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the spatial dimension.
!
!    Input, integer C(N+1), coefficients defining the 
!    lattice simplex.  These should be positive.
!
!    Input/output, integer V(N).  On first call, the input
!    value is not important.  On a repeated call, the input value should
!    be the output value from the previous call.  On output, V contains
!    the next lattice point.
!
!    Input/output, logical MORE.  On input, set MORE to FALSE to 
!    indicate that this is the first call for a given simplex.  Thereafter, 
!    the input value should be the output value from the previous call.  On 
!    output, MORE is TRUE if not only is the returned value V a lattice point,
!    but the routine can be called again for another lattice point.
!    If the output value is FALSE, then no more lattice points were found,
!    and V was reset to 0, and the routine should not be called further
!    for this simplex.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer c(n+1)
  integer c1n
  integer i
  integer i4vec_lcm
  integer j
  integer lhs
  logical more
  integer rhs
  integer term
  integer v(n)

  if ( .not. more ) then

    v(1:n) = 0
    more = .true.

  else

    c1n = i4vec_lcm ( n, c )
    rhs = c1n * c(n+1)

    lhs = 0
    do i = 1, n
      term = 1
      do j = 1, n
        if ( i == j ) then
          term = term * v(j)
        else
          term = term * c(j)
        end if
      end do
      lhs = lhs + term
    end do

    do i = 1, n
      if ( lhs + c1n / c(i) <= rhs ) then
        v(i) = v(i) + 1
        more = .true.
        return
      end if
      lhs = lhs - c1n * v(i) / c(i)
      v(i) = 0
    end do

    more = .false.

  end if

  return
end
subroutine simplex01_lattice_point_num_nd ( d, s, n )

!*****************************************************************************80
!
!! SIMPLEX01_LATTICE_POINT_NUM_ND: count lattice points.
!
!  Discussion:
!
!    The simplex is assumed to be the unit D-dimensional simplex:
!
!    ( (0,0,...,0), (1,0,...,0), (0,1,...,0), ... (0,,0,...,1) )
!
!    or a copy of this simplex scaled by an integer S:
!
!    ( (0,0,...,0), (S,0,...,0), (0,S,...,0), ... (0,,0,...,S) )
!
!    The routine returns the number of integer lattice points that appear
!    inside the simplex or on its boundary.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 July 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Matthias Beck, Sinai Robins,
!    Computing the Continuous Discretely,
!    Springer, 2006,
!    ISBN13: 978-0387291390,
!    LC: QA640.7.B43.
!
!  Parameters:
!
!    Input, integer D, the spatial dimension.
!
!    Input, integer S, the scale factor.
!
!    Output, integer N, the number of lattice points.
!
  implicit  none

  integer d
  integer i
  integer n
  integer s

  n = 1
  do i = 1, d
    n = ( n * ( s + i ) ) / i
  end do

  return
end
subroutine simplex01_volume_nd ( dim_num, volume )

!*****************************************************************************80
!
!! SIMPLEX01_VOLUME_ND computes the volume of the unit simplex in ND.
!
!  Discussion:
!
!    The formula is simple: volume = 1/N!.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    29 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = rk ) VOLUME, the volume.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer dim_num
  real ( kind = rk ) volume

  volume = 1.0D+00
  do i = 1, dim_num
    volume = volume / real ( i, kind = rk )
  end do

  return
end
subroutine simplex_volume_nd ( dim_num, a, volume )

!*****************************************************************************80
!
!! SIMPLEX_VOLUME_ND computes the volume of a simplex in ND.
!
!  Discussion:
!
!    The formula is: 
!
!      volume = 1/N! * det ( A )
!
!    where A is the N by N matrix obtained by subtracting one
!    vector from all the others.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    29 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the dimension of the space.
!
!    Input, real ( kind = rk ) A(DIM_NUM,DIM_NUM+1), the vertices.
!
!    Output, real ( kind = rk ) VOLUME, the volume of the simplex.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num

  real ( kind = rk ) a(dim_num,dim_num+1)
  real ( kind = rk ) b(dim_num,dim_num)
  real ( kind = rk ) det
  integer i
  integer info
  integer j
  integer pivot(dim_num)
  real ( kind = rk ) volume

  b(1:dim_num,1:dim_num) = a(1:dim_num,1:dim_num)
  do j = 1, dim_num
    b(1:dim_num,j) = b(1:dim_num,j) - a(1:dim_num,dim_num+1)
  end do

  call r8ge_fa ( dim_num, b, pivot, info )

  if ( info /= 0 ) then

    volume = -1.0D+00

  else

    call r8ge_det ( dim_num, b, pivot, det )

    volume = abs ( det )
    do i = 1, dim_num
      volume = volume / real ( i, kind = rk )
    end do

  end if

  return
end
function sin_power_int ( a, b, n )

!*****************************************************************************80
!
!! SIN_POWER_INT evaluates the sine power integral.
!
!  Discussion:
!
!    The function is defined by
!
!      SIN_POWER_INT(A,B,N) = Integral ( A <= T <= B ) ( sin ( t ))^n dt
!
!    The algorithm uses the following fact:
!
!      Integral sin^n ( t ) = (1/n) * (
!        sin^(n-1)(t) * cos(t) + ( n-1 ) * Integral sin^(n-2) ( t ) dt )
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, real ( kind = rk ) A, B, the limits of integration.
!
!    Input, integer N, the power of the sine function.
!
!    Output, real ( kind = rk ) SIN_POWER_INT, the value of the integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) ca
  real ( kind = rk ) cb
  integer m
  integer mlo
  integer n
  real ( kind = rk ) sa
  real ( kind = rk ) sb
  real ( kind = rk ) sin_power_int
  real ( kind = rk ) value

  if ( n < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SIN_POWER_INT - Fatal error!'
    write ( *, '(a)' ) '  Power N < 0.'
    value = 0.0D+00
    stop 1
  end if

  sa = sin ( a )
  sb = sin ( b )
  ca = cos ( a )
  cb = cos ( b )

  if ( mod ( n, 2 ) == 0 ) then

    value = b - a
    mlo = 2
  else
    value = ca - cb
    mlo = 3
  end if

  do m = mlo, n, 2
    value = ( real ( m - 1, kind = rk ) * value &
              + sa**( m - 1 ) * ca - sb**( m - 1 ) * cb ) &
      / real ( m, kind = rk )
  end do

  sin_power_int = value

  return
end
subroutine soccer_shape_3d ( point_num, face_num, face_order_max, point_coord, &
  face_order, face_point )

!*****************************************************************************80
!
!! SOCCER_SHAPE_3D describes a truncated icosahedron in 3D.
!
!  Discussion:
!
!    The shape is a truncated icosahedron, which is the design used
!    on a soccer ball.  There are 12 pentagons and 20 hexagons.
!
!    Call SOCCER_SIZE_3D to get the values of POINT_NUM, FACE_NUM, and 
!    FACE_ORDER_MAX, so you can allocate space for the arrays.
!
!    For each face, the face list must be of length FACE_ORDER_MAX.
!    In cases where a face is of lower than maximum order (the
!    12 pentagons, in this case), the extra entries are listed as
!    "-1".
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    08 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    http://mathworld.wolfram.com/TruncatedIcosahedron.html
!
!  Parameters:
!
!    Input, integer POINT_NUM, the number of points (60).
!
!    Input, integer FACE_NUM, the number of faces (32).
!
!    Input, integer FACE_ORDER_MAX, the maximum order of any 
!    face (6).
!
!    Output, real ( kind = rk ) POINT_COORD(3,POINT_NUM), the vertices.
!
!    Output, integer FACE_ORDER(FACE_NUM), the number of
!    vertices per face.
!
!    Output, integer FACE_POINT(FACE_ORDER_MAX,FACE_NUM);
!    FACE_POINT(I,J) contains the index of the I-th point in the J-th face.
!    The points are listed in the counter clockwise direction defined
!    by the outward normal at the face.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer face_num
  integer face_order_max
  integer, parameter :: dim_num = 3
  integer point_num

  integer face_order(face_num)
  integer face_point(face_order_max,face_num)
  real ( kind = rk ) point_coord(dim_num,point_num)
!
!  Set the point coordinates.
!
  point_coord(1:dim_num,1:point_num) = reshape ( (/ &
       -0.100714D+01,   0.153552D+00,   0.067258D+00, &
       -0.960284D+00,   0.0848813D+00, -0.336290D+00, &
       -0.951720D+00,  -0.153552D+00,   0.336290D+00, &
       -0.860021D+00,   0.529326D+00,   0.150394D+00, &
       -0.858000D+00,  -0.290893D+00,  -0.470806D+00, &
       -0.849436D+00,  -0.529326D+00,   0.201774D+00, &
       -0.802576D+00,  -0.597996D+00,  -0.201774D+00, &
       -0.784200D+00,   0.418215D+00,  -0.502561D+00, &
       -0.749174D+00,  -0.0848813D+00,  0.688458D+00, &
       -0.722234D+00,   0.692896D+00,  -0.201774D+00, &
       -0.657475D+00,   0.597996D+00,   0.502561D+00, &
       -0.602051D+00,   0.290893D+00,   0.771593D+00, &
       -0.583675D+00,  -0.692896D+00,   0.470806D+00, &
       -0.579632D+00,  -0.333333D+00,  -0.771593D+00, &
       -0.521710D+00,  -0.418215D+00,   0.771593D+00, &
       -0.505832D+00,   0.375774D+00,  -0.803348D+00, &
       -0.489955D+00,  -0.830237D+00,  -0.336290D+00, &
       -0.403548D+00,   0.000000D+00,  -0.937864D+00, &
       -0.381901D+00,   0.925138D+00,  -0.201774D+00, &
       -0.352168D+00,  -0.666667D+00,  -0.688458D+00, &
       -0.317142D+00,   0.830237D+00,   0.502561D+00, &
       -0.271054D+00,  -0.925138D+00,   0.336290D+00, &
       -0.227464D+00,   0.333333D+00,   0.937864D+00, &
       -0.224193D+00,  -0.993808D+00,  -0.067258D+00, &
       -0.179355D+00,   0.993808D+00,   0.150394D+00, &
       -0.165499D+00,   0.608015D+00,  -0.803348D+00, &
       -0.147123D+00,  -0.375774D+00,   0.937864D+00, &
       -0.103533D+00,   0.882697D+00,  -0.502561D+00, &
       -0.513806D-01,   0.666667D+00,   0.771593D+00, &
        0.000000D+00,   0.000000D+00,   1.021000D+00, &
        0.000000D+00,   0.000000D+00,  -1.021000D+00, &
        0.513806D-01,  -0.666667D+00,  -0.771593D+00, &
        0.103533D+00,  -0.882697D+00,   0.502561D+00, &
        0.147123D+00,   0.375774D+00,  -0.937864D+00, &
        0.165499D+00,  -0.608015D+00,   0.803348D+00, &
        0.179355D+00,  -0.993808D+00,  -0.150394D+00, &
        0.224193D+00,   0.993808D+00,   0.067258D+00, &
        0.227464D+00,  -0.333333D+00,  -0.937864D+00, &
        0.271054D+00,   0.925138D+00,  -0.336290D+00, &
        0.317142D+00,  -0.830237D+00,  -0.502561D+00, &
        0.352168D+00,   0.666667D+00,   0.688458D+00, &
        0.381901D+00,  -0.925138D+00,   0.201774D+00, &
        0.403548D+00,   0.000000D+00,   0.937864D+00, &
        0.489955D+00,   0.830237D+00,   0.336290D+00, &
        0.505832D+00,  -0.375774D+00,   0.803348D+00, &
        0.521710D+00,   0.418215D+00,  -0.771593D+00, &
        0.579632D+00,   0.333333D+00,   0.771593D+00, &
        0.583675D+00,   0.692896D+00,  -0.470806D+00, &
        0.602051D+00,  -0.290893D+00,  -0.771593D+00, &
        0.657475D+00,  -0.597996D+00,  -0.502561D+00, &
        0.722234D+00,  -0.692896D+00,   0.201774D+00, &
        0.749174D+00,   0.0848813D+00, -0.688458D+00, &
        0.784200D+00,  -0.418215D+00,   0.502561D+00, &
        0.802576D+00,   0.597996D+00,   0.201774D+00, &
        0.849436D+00,   0.529326D+00,  -0.201774D+00, &
        0.858000D+00,   0.290893D+00,   0.470806D+00, &
        0.860021D+00,  -0.529326D+00,  -0.150394D+00, &
        0.951720D+00,   0.153552D+00,  -0.336290D+00, &
        0.960284D+00,  -0.0848813D+00,  0.336290D+00, &
        1.007140D+00,  -0.153552D+00,  -0.067258D+00 /), &
    (/ dim_num, point_num /) )
!
!  Set the face orders.
!
  face_order(1:face_num) = (/ &
    6, 6, 5, 6, 5, 6, 5, 6, 6, 6, &
    5, 6, 5, 6, 5, 6, 6, 6, 5, 6, &
    5, 5, 6, 6, 6, 5, 6, 5, 6, 6, &
    5, 6 /)
!
!  Set faces.
!
  face_point(1:face_order_max,1:face_num) = reshape ( (/ &
       30, 43, 47, 41, 29, 23, &
       30, 23, 12,  9, 15, 27, &
       30, 27, 35, 45, 43, -1, &
       43, 45, 53, 59, 56, 47, &
       23, 29, 21, 11, 12, -1, &
       27, 15, 13, 22, 33, 35, &
       47, 56, 54, 44, 41, -1, &
       45, 35, 33, 42, 51, 53, &
       12, 11,  4,  1,  3,  9, &
       29, 41, 44, 37, 25, 21, &
       15,  9,  3,  6, 13, -1, &
       56, 59, 60, 58, 55, 54, &
       53, 51, 57, 60, 59, -1, &
       11, 21, 25, 19, 10,  4, &
       33, 22, 24, 36, 42, -1, &
       13,  6,  7, 17, 24, 22, &
       54, 55, 48, 39, 37, 44, &
       51, 42, 36, 40, 50, 57, &
        4, 10,  8,  2,  1, -1, &
        3,  1,  2,  5,  7,  6, &
       25, 37, 39, 28, 19, -1, &
       55, 58, 52, 46, 48, -1, &
       60, 57, 50, 49, 52, 58, &
       10, 19, 28, 26, 16,  8, &
       36, 24, 17, 20, 32, 40, &
        7,  5, 14, 20, 17, -1, &
       48, 46, 34, 26, 28, 39, &
       50, 40, 32, 38, 49, -1, &
        8, 16, 18, 14,  5,  2, &
       46, 52, 49, 38, 31, 34, &
       16, 26, 34, 31, 18, -1, &
       32, 20, 14, 18, 31, 38 /), (/ face_order_max, face_num /) )

  return
end
subroutine soccer_size_3d ( point_num, edge_num, face_num, face_order_max )

!*****************************************************************************80
!
!! SOCCER_SIZE_3D gives "sizes" for a truncated icosahedron in 3D.
!
!  Discussion:
!
!    The shape is a truncated icosahedron, which is the design used
!    on a soccer ball.  There are 12 pentagons and 20 hexagons.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    http://mathworld.wolfram.com/TruncatedIcosahedron.html
!
!  Parameters:
!
!    Output, integer POINT_NUM, the number of points.
!
!    Output, integer EDGE_NUM, the number of edges.
!
!    Output, integer FACE_NUM, the number of faces.
!
!    Output, integer FACE_ORDER_MAX, the maximum order of any face.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer edge_num
  integer face_num
  integer face_order_max
  integer point_num

  point_num = 60
  edge_num = 90
  face_num = 32
  face_order_max = 6

  return
end
subroutine sort_heap_external ( n, indx, i, j, isgn )

!*****************************************************************************80
!
!! SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
!
!  Discussion:
!
!    The actual list of data is not passed to the routine.  Hence this
!    routine may be used to sort integers, real ( kind = rk )s, numbers, names,
!    dates, shoe sizes, and so on.  After each call, the routine asks
!    the user to compare or interchange two items, until a special
!    return value signals that the sorting is completed.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 February 2004
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis and Herbert Wilf
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer N, the number of items to be sorted.
!
!    Input/output, integer INDX, the main communication signal.
!
!    The user must set INDX to 0 before the first call.
!    Thereafter, the user should not change the value of INDX until
!    the sorting is done.
!
!    On return, if INDX is
!
!      greater than 0,
!      * interchange items I and J;
!      * call again.
!
!      less than 0,
!      * compare items I and J;
!      * set ISGN = -1 if I < J, ISGN = +1 if J < I;
!      * call again.
!
!      equal to 0, the sorting is done.
!
!    Output, integer I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ISGN, results of comparison of elements I 
!    and J.  (Used only when the previous call returned INDX less than 0).
!    ISGN <= 0 means I is less than or equal to J;
!    0 <= ISGN means I is greater than or equal to J.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer, save :: i_save = 0
  integer indx
  integer isgn
  integer j
  integer, save :: j_save = 0
  integer, save :: k = 0
  integer, save :: k1 = 0
  integer n
  integer, save :: n1 = 0
!
!  INDX = 0: This is the first call.
!
  if ( indx == 0 ) then

    i_save = 0
    j_save = 0
    k = n / 2
    k1 = k
    n1 = n
!
!  INDX < 0: The user is returning the results of a comparison.
!
  else if ( indx < 0 ) then

    if ( indx == -2 ) then

      if ( isgn < 0 ) then
        i_save = i_save + 1
      end if

      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return

    end if

    if ( 0 < isgn ) then
      indx = 2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then

      if ( n1 == 1 ) then
        i_save = 0
        j_save = 0
        indx = 0
      else
        i_save = n1
        n1 = n1 - 1
        j_save = 1
        indx = 1
      end if

      i = i_save
      j = j_save
      return

    end if

    k = k - 1
    k1 = k
!
!  0 < INDX, the user was asked to make an interchange.
!
  else if ( indx == 1 ) then

    k1 = k

  end if

  do

    i_save = 2 * k1

    if ( i_save == n1 ) then
      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return
    else if ( i_save <= n1 ) then
      j_save = i_save + 1
      indx = -2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then
      exit
    end if

    k = k - 1
    k1 = k

  end do

  if ( n1 == 1 ) then
    i_save = 0
    j_save = 0
    indx = 0
    i = i_save
    j = j_save
  else
    i_save = n1
    n1 = n1 - 1
    j_save = 1
    indx = 1
    i = i_save
    j = j_save
  end if

  return
end
subroutine sphere_cap_area_2d ( r, h, area )

!*****************************************************************************80
!
!! SPHERE_CAP_AREA_2D computes the surface area of a spherical cap in 2D.
!
!  Discussion:
!
!    Draw any radius of the sphere and note the point P where the radius
!    intersects the sphere.  Consider the point on the radius line which is
!    H units from P.  Draw the circle that lies in the plane perpendicular to
!    the radius, and which intersects the sphere.  The circle divides the sphere
!    into two pieces, and the corresponding disk divides the solid sphere into
!    two pieces.  The spherical cap is the part of the solid sphere that
!    includes the point P.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) R, the radius of the sphere.
!
!    Input, real ( kind = rk ) H, the "height" of the spherical cap. 
!    H must be between 0 and 2 * R.
!
!    Output, real ( kind = rk ) AREA, the area of the spherical cap.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) area
  real ( kind = rk ) h
  real ( kind = rk ) r
  real ( kind = rk ) r8_asin
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) theta

  if ( h <= 0.0D+00 ) then
    area = 0.0D+00
  else if ( 2.0D+00 * r <= h ) then
    area = 2.0D+00 * r8_pi * r
  else

    theta = 2.0D+00 * r8_asin ( sqrt ( r * r - ( r - h )**2 ) / r )
    area = r * theta

    if ( r <= h ) then
      area = 2.0D+00 * r8_pi * r - area
    end if

  end if

  return
end
subroutine sphere_cap_area_3d ( r, h, area )

!*****************************************************************************80
!
!! SPHERE_CAP_AREA_3D computes the surface area of a spherical cap in 3D.
!
!  Discussion:
!
!    Draw any radius of the sphere and note the point P where the radius
!    intersects the sphere.  Consider the point on the radius line which is
!    H units from P.  Draw the circle that lies in the plane perpendicular to
!    the radius, and which intersects the sphere.  The circle divides the sphere
!    into two pieces, and the corresponding disk divides the solid sphere into
!    two pieces.  The spherical cap is the part of the solid sphere that
!    includes the point P.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) R, the radius of the sphere.
!
!    Input, real ( kind = rk ) H, the "height" of the spherical cap. 
!    H must be between 0 and 2 * R.
!
!    Output, real ( kind = rk ) AREA, the area of the spherical cap.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) area
  real ( kind = rk ) h
  real ( kind = rk ) r
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00

  if ( h <= 0.0D+00 ) then
    area = 0.0D+00
  else if ( 2.0D+00 * r <= h ) then
    area = 4.0D+00 * r8_pi * r * r
  else
    area = 2.0D+00 * r8_pi * r * h
  end if

  return
end
subroutine sphere_cap_area_nd ( dim_num, r, h, area )

!*****************************************************************************80
!
!! SPHERE_CAP_AREA_ND computes the area of a spherical cap in ND.
!
!  Discussion:
!
!    The spherical cap is a portion of the surface of the sphere:
!
!      sum ( X(1:N)^2 ) = R^2
!
!    which is no more than H units from the uppermost point on the sphere.
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
!  Reference:
!
!    Thomas Ericson, Victor Zinoviev,
!    Codes on Euclidean Spheres,
!    Elsevier, 2001, pages 439-441.
!    QA166.7 E75
!
!  Parameters:
!
!    Input, integer DIM_NUM, the dimension of the space.
!
!    Input, real ( kind = rk ) R, the radius of the sphere.
!
!    Input, real ( kind = rk ) H, the "thickness" of the spherical cap,
!    which is normally between 0 and 2 * R.
!
!    Output, real ( kind = rk ) AREA, the area of the spherical cap.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) area
  real ( kind = rk ) area2
  real ( kind = rk ) h
  real ( kind = rk ) haver_sine
  integer i
  integer dim_num
  real ( kind = rk ) r
  real ( kind = rk ) r8_asin
  real ( kind = rk ) sphere_k
  real ( kind = rk ) theta
  real ( kind = rk ) ti
  real ( kind = rk ) tj
  real ( kind = rk ) tk

  if ( h <= 0.0D+00 ) then
    area = 0.0D+00
    return
  end if

  if ( 2.0D+00 * r <= h ) then
    call sphere_imp_area_nd ( dim_num, r, area )
    return
  end if
!
!  For cases where R < H < 2 * R, work with the complementary region.
!
  haver_sine = sqrt ( ( 2.0D+00 * r - h ) * h )

  theta = r8_asin ( haver_sine / r )

  if ( dim_num < 1 ) then

    area = -1.0D+00
    return

  else if ( dim_num == 1 ) then

    area = 0.0D+00

  else if ( dim_num == 2 ) then

    area = 2.0D+00 * theta * r

  else

    ti = theta

    tj = ti
    ti = 1.0D+00 - cos ( theta )

    do i = 2, dim_num-2
      tk = tj
      tj = ti
      ti = ( real ( i - 1, kind = rk ) * tk &
        - cos ( theta ) * sin ( theta ) ** ( i - 1 ) ) &
        / real ( i, kind = rk )
    end do

    area = sphere_k ( dim_num-1 ) * ti * r ** ( dim_num - 1 )

  end if
!
!  Adjust for cases where R < H < 2R.
!
  if ( r < h ) then
    call sphere_imp_area_nd ( dim_num, r, area2 )
    area = area2 - area
  end if

  return
end
subroutine sphere_cap_volume_2d ( r, h, volume )

!*****************************************************************************80
!
!! SPHERE_CAP_VOLUME_2D computes the volume of a spherical cap in 2D.
!
!  Discussion:
!
!    Draw any radius R of the circle and denote as P the point where the
!    radius intersects the circle.  Now consider the point Q which lies
!    on the radius and which is H units from P.  The line which is
!    perpendicular to the radius R and passes through Q divides the
!    circle into two pieces.  The piece including the point P is the
!    spherical (circular) cap of height (or thickness) H.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) R, the radius of the sphere.
!
!    Input, real ( kind = rk ) H, the "height" of the spherical cap.  H must
!    be between 0 and 2 * R.
!
!    Output, real ( kind = rk ) VOLUME, the volume (area) of the spherical cap.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) h
  real ( kind = rk ) r
  real ( kind = rk ) r8_asin
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) theta
  real ( kind = rk ) volume

  if ( h <= 0.0D+00 ) then

    volume = 0.0D+00

  else if ( 2.0D+00 * r <= h ) then

    volume = r8_pi * r * r

  else

    theta = 2.0D+00 * r8_asin ( sqrt ( r * r - ( r - h )**2 ) / r )
    volume = r * r * ( theta - sin ( theta ) ) / 2.0D+00

    if ( r < h ) then
      volume = r8_pi * r * r - volume
    end if

  end if

  return
end
subroutine sphere_cap_volume_3d ( r, h, volume )

!*****************************************************************************80
!
!! SPHERE_CAP_VOLUME_3D computes the volume of a spherical cap in 3D.
!
!  Discussion:
!
!    Draw any radius of the sphere and note the point P where the radius
!    intersects the sphere.  Consider the point on the radius line which is
!    H units from P.  Draw the circle that lies in the plane perpendicular to
!    the radius, and which intersects the sphere.  The circle divides the sphere
!    into two pieces, and the corresponding disk divides the solid sphere into
!    two pieces.  The part of the solid sphere that includes the point P
!    is the spherical cap of height (or thickness) H.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) R, the radius of the sphere.
!
!    Input, real ( kind = rk ) H, the "height" of the spherical cap.  H must
!    be between 0 and 2 * R.
!
!    Output, real ( kind = rk ) VOLUME, the volume of the spherical cap.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) h
  real ( kind = rk ) r
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) volume

  if ( h <= 0.0D+00 ) then
    volume = 0.0D+00
  else if ( 2.0D+00 * r <= h ) then
    volume = ( 4.0D+00 / 3.0D+00 ) * r8_pi * r * r * r
  else
    volume = ( 1.0D+00 / 3.0D+00 ) * r8_pi * h * h * ( 3.0D+00 * r - h )
  end if

  return
end
subroutine sphere_cap_volume_nd ( dim_num, r, h, volume )

!*****************************************************************************80
!
!! SPHERE_CAP_VOLUME_ND computes the volume of a spherical cap in ND.
!
!  Discussion:
!
!    The spherical cap is a portion of the surface and interior of the sphere:
!
!      sum ( X(1:N)^2 ) <= R^2
!
!    which is no more than H units from some point P on the sphere.
!
!
!    The algorithm proceeds from the observation that the N-dimensional
!    sphere can be parameterized by a quantity RC that runs along the
!    radius from the center to the point P.  The value of RC at the
!    base of the spherical cap is (R-H) and at P it is R.  We intend to
!    use RC as our integration parameeter.
!
!    The volume of the spherical cap is then the integral, as RC goes
!    from (R-H) to R, of the N-1 dimensional volume of the sphere
!    of radius RS, where RC^2 + RS^2 = R^2.
!
!    The volume of the N-1 dimensional sphere of radius RS is simply 
!    some constants times RS^(N-1).
! 
!    After factoring out the constant terms, and writing RC = R * cos ( T ),
!    and RS = R * sin ( T ), and letting 
!      T_MAX = arc_sine ( sqrt ( ( 2.0D+00 * r - h ) * h / r ) ),
!    the "interesting part" of our integral becomes
!
!      constants * R^N * Integral ( T = 0 to T_MAX ) sin^N ( T ) dT
!
!    The integral of sin^N ( T ) dT can be handled by recursion.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the dimension of the space.
!
!    Input, real ( kind = rk ) R, the radius of the sphere.
!
!    Input, real ( kind = rk ) H, the "thickness" of the spherical cap,
!    which is normally between 0 and 2 * R.
!
!    Output, real ( kind = rk ) VOLUME, the volume of the spherical cap.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) angle
  real ( kind = rk ) factor1
  real ( kind = rk ) factor2
  real ( kind = rk ) h
  integer dim_num
  real ( kind = rk ) r
  real ( kind = rk ) r8_asin
  real ( kind = rk ) sin_power_int
  real ( kind = rk ) sphere01_volume_nd
  real ( kind = rk ) volume
  real ( kind = rk ) volume2

  if ( h <= 0.0D+00 ) then
    volume = 0.0D+00
    return
  end if

  if ( 2.0D+00 * r <= h ) then
    call sphere_imp_volume_nd ( dim_num, r, volume )
    return
  end if

  if ( dim_num < 1 ) then

    volume = -1.0D+00

  else if ( dim_num == 1 ) then

    volume = h

  else

    factor1 = sphere01_volume_nd ( dim_num - 1 )

    angle = r8_asin ( sqrt ( ( 2.0D+00 * r - h ) * h / r ) )

    factor2 = sin_power_int ( 0.0D+00, angle, dim_num )

    volume = factor1 * factor2 * r**dim_num

    if ( r < h ) then
      call sphere_imp_volume_nd ( dim_num, r, volume2 )
      volume = volume2 - volume
    end if

  end if

  return
end
subroutine sphere_dia2imp_3d ( p1, p2, r, pc )

!*****************************************************************************80
!
!! SPHERE_DIA2IMP_3D converts a diameter to an implicit sphere in 3D.
!
!  Discussion:
!
!    An implicit sphere in 3D satisfies the equation:
!
!      ( P(1) - PC(1) )^2 + ( P(2) - PC(2) )^2 + ( P(3) - PC(3) )^2 = R^2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(3), P2(3), are two points which form a 
!    diameter of the sphere.
!
!    Output, real ( kind = rk ) R, the computed radius of the sphere.
!
!    Output, real ( kind = rk ) PC(3), the computed center of the sphere.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) p1(3)
  real ( kind = rk ) p2(3)
  real ( kind = rk ) pc(3)
  real ( kind = rk ) r
  real ( kind = rk ) r8vec_norm_affine

  r = 0.5D+00 * r8vec_norm_affine ( 3, p1, p2 )

  pc(1:3) = 0.5D+00 * ( p1(1:3) + p2(1:3) )

  return
end
subroutine sphere_distance_xyz ( xyz1, xyz2, dist )

!*****************************************************************************80
!
!! SPHERE_DISTANCE_XYZ computes great circle distances on a sphere.
!
!  Discussion:
!
!    XYZ coordinates are used.
!
!    We assume the points XYZ1 and XYZ2 lie on the same sphere.
!
!    This computation is a special form of the Vincenty formula.
!    It should be less sensitive to errors associated with very small 
!    or very large angular separations.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    "Great-circle distance",
!    Wikipedia.
!
!  Parameters:
!
!    Input, real ( kind = rk ) XYZ1(3), the coordinates of the first point.
!
!    Input, real ( kind = rk ) XYZ2(3), the coordinates of the second point.
!
!    Output, real ( kind = rk ) DIST, the great circle distance between
!    the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) bot
  real ( kind = rk ) dist
  real ( kind = rk ) lat1
  real ( kind = rk ) lat2
  real ( kind = rk ) lon1
  real ( kind = rk ) lon2
  real ( kind = rk ) r
  real ( kind = rk ) r8_asin
  real ( kind = rk ) r8_atan
  real ( kind = rk ) r8vec_norm
  real ( kind = rk ) top
  real ( kind = rk ) xyz1(3)
  real ( kind = rk ) xyz2(3)

  r = r8vec_norm ( 3, xyz1 )

  lat1 = r8_asin ( xyz1(3) )
  lon1 = r8_atan ( xyz1(2), xyz1(1) )

  lat2 = r8_asin ( xyz2(3) )
  lon2 = r8_atan ( xyz2(2), xyz2(1) )

  top = ( cos ( lat2 ) * sin ( lon1 - lon2 ) )**2 &
      + ( cos ( lat1 ) * sin ( lat2 ) &
      -   sin ( lat1 ) * cos ( lat2 ) * cos ( lon1 - lon2 ) )**2

  top = sqrt ( top )

  bot = sin ( lat1 ) * sin ( lat2 ) &
      + cos ( lat1 ) * cos ( lat2 ) * cos ( lon1 - lon2 )

  dist = r * atan2 ( top, bot )

  return
end
subroutine sphere_distance1 ( lat1, lon1, lat2, lon2, r, dist )

!*****************************************************************************80
!
!! SPHERE_DISTANCE1 computes great circle distances on a sphere.
!
!  Discussion:
!
!    This computation is based on the law of cosines for spheres.
!    This formula can suffer from rounding errors when the angular
!    distances are small.
!
!    Here we assume that latitude is measured in radians, and goes
!    from -PI at the south pole to +PI at the north pole.  
!
!    Similarly, longitude is measured in radians, and goes from 0 to 2 PI.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 February 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    "Great-circle distance",
!    Wikipedia.
!
!  Parameters:
!
!    Input, real ( kind = rk ) LAT1, LON1, the latitude and longitude of 
!    the first point.
!
!    Input, real ( kind = rk ) LAT2, LON2, the latitude and longitude of 
!    the second point.
!
!    Input, real ( kind = rk ) R, the radius of the sphere.
!
!    Output, real ( kind = rk ) DIST, the great circle distance between
!    the points, measured in the same units as R.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) c
  real ( kind = rk ) dist
  real ( kind = rk ) lat1
  real ( kind = rk ) lat2
  real ( kind = rk ) lon1
  real ( kind = rk ) lon2
  real ( kind = rk ) r

  c = cos ( lat1 ) * cos ( lat2 ) * cos ( lon1 - lon2 ) &
    + sin ( lat1 ) * sin ( lat2 )

  dist = r * acos ( c )

  return
end
subroutine sphere_distance2 ( lat1, lon1, lat2, lon2, r, dist )

!*****************************************************************************80
!
!! SPHERE_DISTANCE2 computes great circle distances on a sphere.
!
!  Discussion:
!
!    This computation is written in terms of haversines, and can be more
!    accurate when measuring small angular distances.  It can be somewhat
!    inaccurate when the two points are antipodal.
!
!    Here we assume that latitude is measured in radians, and goes
!    from -PI at the south pole to +PI at the north pole.  
!
!    Similarly, longitude is measured in radians, and goes from 0 to 2 PI.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 February 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    "Great-circle distance",
!    Wikipedia.
!
!  Parameters:
!
!    Input, real ( kind = rk ) LAT1, LON1, the latitude and longitude of 
!    the first point.
!
!    Input, real ( kind = rk ) LAT2, LON2, the latitude and longitude of 
!    the second point.
!
!    Input, real ( kind = rk ) R, the radius of the sphere.
!
!    Output, real ( kind = rk ) DIST, the great circle distance between
!    the points, measured in the same units as R.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) dist
  real ( kind = rk ) lat1
  real ( kind = rk ) lat2
  real ( kind = rk ) lon1
  real ( kind = rk ) lon2
  real ( kind = rk ) r
  real ( kind = rk ) s

  s = ( sin ( ( lat1 - lat2 ) / 2.0D+00 ) )**2 &
    + cos ( lat1 ) * cos ( lat2 ) * ( sin ( ( lon1 - lon2 ) / 2.0D+00 ) )**2
  s = sqrt ( s )

  dist = 2.0D+00 * r * asin ( s )

  return
end
subroutine sphere_distance3 ( lat1, lon1, lat2, lon2, r, dist )

!*****************************************************************************80
!
!! SPHERE_DISTANCE3 computes great circle distances on a sphere.
!
!  Discussion:
!
!    This computation is a special form of the Vincenty formula.
!    It should be less sensitive to errors associated with very small 
!    or very large angular separations.
!
!    Here we assume that latitude is measured in radians, and goes
!    from -PI at the south pole to +PI at the north pole.  
!
!    Similarly, longitude is measured in radians, and goes from 0 to 2 PI.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 February 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    "Great-circle distance",
!    Wikipedia.
!
!  Parameters:
!
!    Input, real ( kind = rk ) LAT1, LON1, the latitude and longitude of 
!    the first point.
!
!    Input, real ( kind = rk ) LAT2, LON2, the latitude and longitude of 
!    the second point.
!
!    Input, real ( kind = rk ) R, the radius of the sphere.
!
!    Output, real ( kind = rk ) DIST, the great circle distance between
!    the points, measured in the same units as R.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) bot
  real ( kind = rk ) dist
  real ( kind = rk ) lat1
  real ( kind = rk ) lat2
  real ( kind = rk ) lon1
  real ( kind = rk ) lon2
  real ( kind = rk ) r
  real ( kind = rk ) top

  top = ( cos ( lat2 ) * sin ( lon1 - lon2 ) )**2 &
      + ( cos ( lat1 ) * sin ( lat2 ) &
      -   sin ( lat1 ) * cos ( lat2 ) * cos ( lon1 - lon2 ) )**2

  top = sqrt ( top )

  bot = sin ( lat1 ) * sin ( lat2 ) &
      + cos ( lat1 ) * cos ( lat2 ) * cos ( lon1 - lon2 )

  dist = r * atan2 ( top, bot )

  return
end
subroutine sphere_exp_contains_point_3d ( p1, p2, p3, p4, p, inside )

!*****************************************************************************80
!
!! SPHERE_EXP_CONTAINS_POINT_3D: does an explicit sphere contain a point in 3D.
!
!  Discussion:
!
!    An explicit sphere in 3D is determined by four points,
!    which should be distinct, and not coplanar.
!
!    The computation checks the determinant of the 5 by 5 matrix:
!
!      x1  y1  z1  x1^2+y1^2+z1^2  1
!      x2  y2  z2  x2^2+y2^2+z2^2  1
!      x3  y3  z3  x3^2+y3^2+z3^2  1
!      x4  y4  z4  x4^2+y4^2+z4^2  1
!      x   y   z   x^2 +y^2 +z^2   1
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(3), P2(3), P3(3), P4(3),
!    four distinct noncoplanar points on the sphere.
!
!    Input, real ( kind = rk ) P(3), the coordinates of a point, whose
!    position relative to the sphere is desired.
!
!    Output, logical INSIDE, is TRUE if the point is in the sphere.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) a(5,5)
  real ( kind = rk ) det
  logical inside
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) p3(dim_num)
  real ( kind = rk ) p4(dim_num)
  real ( kind = rk ) r8mat_det_5d
!
!  Compute the determinant.
!
  a(1,1:dim_num) = p1(1:dim_num)
  a(1,4) = sum ( p1(1:dim_num)**2 )
  a(1,5) = 1.0D+00

  a(2,1:dim_num) = p2(1:dim_num)
  a(2,4) = sum ( p2(1:dim_num)**2 )
  a(2,5) = 1.0D+00

  a(3,1:dim_num) = p3(1:dim_num)
  a(3,4) = sum ( p3(1:dim_num)**2 )
  a(3,5) = 1.0D+00

  a(4,1:dim_num) = p4(1:dim_num)
  a(4,4) = sum ( p4(1:dim_num)**2 )
  a(4,5) = 1.0D+00

  a(5,1:dim_num) = p(1:dim_num)
  a(5,4) = sum ( p(1:dim_num)**2 )
  a(5,5) = 1.0D+00

  det = r8mat_det_5d ( a )

  if ( det < 0.0D+00 ) then
    inside = .false.
  else if ( 0.0D+00 <= det ) then
    inside = .true.
  end if

  return
end
subroutine sphere_exp_point_near_3d ( p1, p2, p3, p4, p, pn )

!*****************************************************************************80
!
!! SPHERE_EXP_POINT_NEAR_3D: nearest point on explicit sphere to a point in 3D.
!
!  Discussion:
!
!    An explicit sphere in 3D is determined by four points,
!    which should be distinct, and not coplanar.
!
!    If the center of the sphere is PC, and the point is P, then
!    the desired point lies at a positive distance R along the vector 
!    P-PC unless P = PC in which case any point on the sphere is "nearest".
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(3), P2(3), P3(3), P4(3),
!    four distinct noncoplanar points on the sphere.
!
!    Input, real ( kind = rk ) P(3), a point whose nearest point on the 
!    sphere is desired.
!
!    Output, real ( kind = rk ) PN(3), the nearest point on the sphere.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) norm
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) p3(dim_num)
  real ( kind = rk ) p4(dim_num)
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) pn(dim_num)
  real ( kind = rk ) r
!
!  Find the center.
!
  call sphere_exp2imp_3d ( p1, p2, p3, p4, r, pc )
!
!  If P = PC, bail out now.
!
  norm = sqrt ( sum ( ( p(1:dim_num) - pc(1:dim_num) )**2 ) )

  if ( norm == 0.0D+00 ) then
    pn(1:dim_num) = pc(1:dim_num) + r / sqrt ( real ( dim_num, kind = rk ) )
    return
  end if
!
!  Compute the nearest point.
!
  pn(1:dim_num) = pc(1:dim_num) + r * ( p(1:dim_num) - pc(1:dim_num) ) / norm

  return
end
subroutine sphere_exp2imp_3d ( p1, p2, p3, p4, r, pc )

!*****************************************************************************80
!
!! SPHERE_EXP2IMP_3D converts a sphere from explicit to implicit form in 3D.
!
!  Discussion:
!
!    An explicit sphere in 3D is determined by four points,
!    which should be distinct, and not coplanar.
!
!    An implicit sphere in 3D satisfies the equation:
!
!      sum ( ( P(1:DIM_NUM) - PC(1:DIM_NUM) )^2 ) = R^2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 February 2005
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
!  Parameters:
!
!    Input, real ( kind = rk ) P1(3), P2(3), P3(3), P4(3),
!    four distinct noncoplanar points on the sphere.
!
!    Output, real ( kind = rk ) R, PC(3), the radius and the center
!    of the sphere.  If the linear system is
!    singular, then R = -1, PC(1:3) = 0.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) p3(dim_num)
  real ( kind = rk ) p4(dim_num)
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) r
  real ( kind = rk ) tetra(dim_num,4)

  tetra(1:dim_num,1:4) = reshape ( (/ &
    p1(1:dim_num), p2(1:dim_num), p3(1:dim_num), p4(1:dim_num) /), &
    (/ dim_num, 4 /) )

  call tetrahedron_circumsphere_3d ( tetra, r, pc )

  return
end
subroutine sphere_exp2imp_nd ( n, p, r, pc )

!*****************************************************************************80
!
!! SPHERE_EXP2IMP_ND finds an N-dimensional sphere through N+1 points.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    25 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the spatial dimension.
!
!    Input, real ( kind = rk ) P(N,N+1), the points.
!
!    Output, real ( kind = rk ) R, PC(N), the radius and center of the
!    sphere.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,n+1)
  integer i
  integer info
  integer j
  real ( kind = rk ) pc(n)
  real ( kind = rk ) r
  real ( kind = rk ) p(n,n+1)
!
!  Set up the linear system.
!
  a(1:n,1:n) = transpose ( p(1:n,2:n+1) )

  do j = 1, n
    a(1:n,j) = a(1:n,j) - p(j,1)
  end do

  do i = 1, n
    a(i,n+1) = sum ( a(i,1:n)**2 )
  end do
!
!  Solve the linear system.
!
  call r8mat_solve ( n, 1, a, info )
!
!  If the system was singular, return a consolation prize.
!
  if ( info /= 0 ) then
    r = -1.0D+00
    pc(1:n) = 0.0D+00
    return
  end if
!
!  Compute the radius and center.
!
  r = 0.5D+00 * sqrt ( sum ( a(1:n,n+1)**2 ) )

  pc(1:n) = p(1:n,1) + 0.5D+00 * a(1:n,n+1)

  return
end
subroutine sphere_imp_area_3d ( r, area )

!*****************************************************************************80
!
!! SPHERE_IMP_AREA_3D computes the surface area of an implicit sphere in 3D.
!
!  Discussion:
!
!    An implicit sphere in 3D satisfies the equation:
!
!      sum ( ( P(1:DIM_NUM) - PC(1:DIM_NUM) )^2 ) = R^2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    13 August 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) R, the radius of the sphere.
!
!    Output, real ( kind = rk ) AREA, the area of the sphere.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) area
  real ( kind = rk ) r
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00

  area = 4.0D+00 * r8_pi * r * r

  return
end
subroutine sphere_imp_area_nd ( dim_num, r, area )

!*****************************************************************************80
!
!! SPHERE_IMP_AREA_ND computes the surface area of an implicit sphere in ND.
!
!  Discussion:
!
!    An implicit sphere in ND satisfies the equation:
!
!      sum ( ( P(1:DIM_NUM) - pc(1:DIM_NUM) )^2 ) = R^2
!
!    DIM_NUM   Area
!
!    2      2       * PI   * R
!    3      4       * PI   * R^2
!    4      2       * PI^2 * R^3
!    5      (8/3)   * PI^2 * R^4
!    6                PI^3 * R^5
!    7      (16/15) * PI^3 * R^6
!
!    Sphere_Area ( DIM_NUM, R ) = 
!      2 * PI^(DIM_NUM/2) * R^(DIM_NUM-1) / Gamma ( DIM_NUM / 2 )
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the dimension of the space.
!
!    Input, real ( kind = rk ) R, the radius of the sphere.
!
!    Output, real ( kind = rk ) AREA, the area of the sphere.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) area
  integer dim_num
  real ( kind = rk ) r
  real ( kind = rk ) sphere01_area_nd

  area = r**( dim_num - 1  ) * sphere01_area_nd ( dim_num )

  return
end
subroutine sphere_imp_contains_point_3d ( r, pc, p, inside )

!*****************************************************************************80
!
!! SPHERE_IMP_CONTAINS_POINT_3D: point in implicit sphere in 3D?
!
!  Discussion:
!
!    An implicit sphere in 3D satisfies the equation:
!
!      sum ( ( P(1:DIM_NUM) - PC(1:DIM_NUM) )^2 ) = R^2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) R, the radius of the sphere.
!
!    Input, real ( kind = rk ) PC(3), the center of the sphere.
!
!    Input, real ( kind = rk ) P(3), the point to be checked.
!
!    Output, logical INSIDE, is TRUE if the point is
!    inside the sphere.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  logical inside
  real ( kind = rk ) p(3)
  real ( kind = rk ) pc(3)
  real ( kind = rk ) r

  if ( sum ( ( p(1:3) - pc(1:3) ) ** 2 ) <= r * r ) then
    inside = .true.
  else
    inside = .false.
  end if

  return
end
subroutine sphere_imp_line_project_3d ( r, pc, n, p, maxpnt2, n2, pp, &
  theta_min, theta_max )

!*****************************************************************************80
!
!! SPHERE_IMP_LINE_PROJECT_3D projects a line onto an implicit sphere in 3D.
!
!  Discussion:
!
!    An implicit sphere in 3D satisfies the equation:
!
!      sum ( ( P(1:DIM_NUM) - PC(1:DIM_NUM) )^2 ) = R^2
!
!    The line to be projected is specified as a sequence of points.
!    If two successive points subtend a small angle, then the second
!    point is essentially dropped.  If two successive points subtend
!    a large angle, then intermediate points are inserted, so that
!    the projected line stays closer to the sphere.
!
!    Note that if any P coincides with the center of the sphere, then
!    its projection is mathematically undefined.  PP will
!    be returned as the center.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) R, the radius of the sphere.  If R is
!    zero, PP will be returned as the pc, and if R is
!    negative, points will end up diametrically opposite from where
!    you would expect them for a positive R.
!
!    Input, real ( kind = rk ) PC(3), the center of the sphere.
!
!    Input, integer N, the number of points on the line that is
!    to be projected.
!
!    Input, real ( kind = rk ) P(3,N), the coordinates of
!    the points on the line that is to be projected.
!
!    Input, integer MAXPNT2, the maximum number of points on the
!    projected line.  Even if the routine thinks that more points are needed,
!    no more than MAXPNT2 will be generated.
!
!    Output, integer N2, the number of points on the projected
!    line.  N2 can be zero, if the line has an angular projection of less
!    than THETA_MIN radians.
!
!    Output, real ( kind = rk ) PP(3,N2), the coordinates
!    of the points representing the projected line.  These points lie on the
!    sphere.  Successive points are separated by at least THETA_MIN
!    radians, and by no more than THETA_MAX radians.
!
!    Input, real ( kind = rk ) THETA_MIN, THETA_MAX, the minimum and maximum
!    angular projections allowed between successive projected points.
!    If two successive points on the original line have projections
!    separated by more than THETA_MAX radians, then intermediate points
!    will be inserted, in an attempt to keep the line closer to the
!    sphere.  If two successive points are separated by less than
!    THETA_MIN radians, then the second point is dropped, and the
!    line from the first point to the next point is considered.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxpnt2
  integer, parameter :: dim_num = 3
  integer n

  real ( kind = rk ) alpha
  real ( kind = rk ) ang3d
  real ( kind = rk ) dot
  integer i
  integer j
  integer nfill
  integer n2
  real ( kind = rk ) p(dim_num,n)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) pd(dim_num)
  real ( kind = rk ) pp(dim_num,maxpnt2)
  real ( kind = rk ) r
  real ( kind = rk ) r8_acos
  real ( kind = rk ) theta_max
  real ( kind = rk ) theta_min
  real ( kind = rk ) tnorm
!
!  Check the input.
!
  if ( r == 0.0D+00 ) then
    n2 = 0
    return
  end if

  p1(1:dim_num) = pc(1:dim_num)
  p2(1:dim_num) = pc(1:dim_num)

  n2 = 0

  do i = 1, n

    if ( all ( p(1:dim_num,i) == pc(1:dim_num) ) ) then

    else

      p1(1:dim_num) = p2(1:dim_num)

      alpha = sqrt ( sum ( ( p(1:dim_num,i) - pc(1:dim_num) )**2 ) )

      p2(1:dim_num) = pc(1:dim_num) &
        + r * ( p(1:dim_num,i) - pc(1:dim_num) ) / alpha
!
!  If we haven't gotten any points yet, take this point as our start.
!
      if ( n2 == 0 ) then

        n2 = n2 + 1
        pp(1:dim_num,n2) = p2(1:dim_num)
!
!  Compute the angular projection of P1 to P2.
!
      else if ( 1 <= n2 ) then

        dot = sum ( ( p1(1:dim_num) - pc(1:dim_num) ) &
                  * ( p2(1:dim_num) - pc(1:dim_num) ) )

        ang3d = r8_acos (  dot / ( r * r ) )
!
!  If the angle is at least THETA_MIN, (or it's the last point),
!  then we will draw a line segment.
!
        if ( theta_min < abs ( ang3d ) .or. i == n ) then
!
!  Now we check to see if the line segment is too long.
!
          if ( theta_max < abs ( ang3d ) ) then

            nfill = int ( abs ( ang3d ) / theta_max )

            do j = 1, nfill-1

              pd(1:dim_num) = &
                ( real ( nfill - j, kind = rk ) &
                * ( p1(1:dim_num) - pc(1:dim_num) ) &
                + real (         j, kind = rk ) &
                * ( p2(1:dim_num) - pc(1:dim_num) ) )

              tnorm = sqrt ( sum ( pd(1:dim_num)**2 ) )

              if ( tnorm /= 0.0D+00 ) then
                pd(1:dim_num) = pc(1:dim_num) + r * pd(1:dim_num) / tnorm
                n2 = n2 + 1
                pp(1:dim_num,n2) = pd(1:dim_num)
              end if

            end do

          end if
!
!  Now tack on the projection of point 2.
!
          n2 = n2 + 1
          pp(1:dim_num,n2) = p2(1:dim_num)

        end if

      end if

    end if

  end do

  return
end
subroutine sphere_imp_local2xyz_3d ( r, pc, theta, phi, p )

!*****************************************************************************80
!
!! SPHERE_IMP_LOCAL2XYZ_3D: local to XYZ coordinates on implicit sphere in 3D.
!
!  Discussion:
!
!    An implicit sphere in 3D satisfies the equation:
!
!      sum ( ( P(1:DIM_NUM) - PC(1:DIM_NUM) )^2 ) = R^2
!
!    The "local" spherical coordinates of a point are two angles, THETA and PHI.
!    PHI measures the angle that the vector from the origin to the point
!    makes with the positive Z axis.  THETA measures the angle that the
!    projection of the vector onto the XY plane makes with the positive X axis.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) R, the radius of the sphere.
!
!    Input, real ( kind = rk ) PC(3), the center of the sphere.
!
!    Input, real ( kind = rk ) THETA, PHI, the local (THETA,PHI) spherical
!    coordinates of a point on the sphere.  THETA and PHI are angles,
!    measured in radians.  Usually, 0 <= THETA < 2 * PI, and 0 <= PHI <= PI.
!
!    Output, real ( kind = rk ) P(3), the XYZ coordinates of the point.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) phi
  real ( kind = rk ) r
  real ( kind = rk ) theta

  p(1) = pc(1) + r * sin ( phi ) * cos ( theta )
  p(2) = pc(2) + r * sin ( phi ) * sin ( theta )
  p(3) = pc(3) + r * cos ( phi )

  return
end
subroutine sphere_imp_point_near_3d ( r, pc, p, pn )

!*****************************************************************************80
!
!! SPHERE_IMP_POINT_NEAR_3D: nearest point on implicit sphere to a point in 3D.
!
!  Discussion:
!
!    An implicit sphere in 3D satisfies the equation:
!
!      sum ( ( P(1:DIM_NUM) - PC(1:DIM_NUM) )^2 ) = R^2
!
!    If the center of the sphere is PC, and the point is P, then
!    the desired point lies at a positive distance R along the vector 
!    P-PC unless P = PC, in which case any point on the sphere is "nearest".
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) R, the radius of the sphere.
!
!    Input, real ( kind = rk ) PC(3), the center of the sphere.
!
!    Input, real ( kind = rk ) P(3), a point whose
!    nearest point on the sphere is desired.
!
!    Output, real ( kind = rk ) PN(3), the nearest point on the sphere.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) norm
  real ( kind = rk ) p(3)
  real ( kind = rk ) pc(3)
  real ( kind = rk ) pn(3)
  real ( kind = rk ) r
!
!  If P = PC, bail out now.
!
  norm = sqrt ( sum ( ( p(1:dim_num) - pc(1:dim_num) )**2 ) )

  if ( norm == 0.0D+00 ) then
    pn(1:dim_num) = pc(1:dim_num) + r / sqrt ( real ( dim_num, kind = rk ) )
    return
  end if
!
!  Compute the nearest point.
!
  pn(1:dim_num) = pc(1:dim_num) + r * ( p(1:dim_num) - pc(1:dim_num) ) / norm

  return
end
subroutine sphere_imp_point_project_3d ( r, pc, p, pp )

!*****************************************************************************80
!
!! SPHERE_IMP_POINT_PROJECT_3D projects a point onto an implicit sphere in 3D.
!
!  Discussion:
!
!    An implicit sphere in 3D satisfies the equation:
!
!      sum ( ( P(1:DIM_NUM) - PC(1:DIM_NUM) )^2 ) = R^2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) R, the radius of the sphere.
!
!    Input, real ( kind = rk ) PC(3), the center of the sphere.
!
!    Input, real ( kind = rk ) P(3), a point.
!
!    Output, real ( kind = rk ) PP(3), the projected point.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) norm
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) pp(dim_num)
  real ( kind = rk ) r

  if ( r == 0.0D+00 ) then

    pp(1:dim_num) = pc(1:dim_num)

  else if ( all ( p(1:dim_num) == pc(1:dim_num) ) ) then

    pp(1:dim_num) = pc(1:dim_num) + r / sqrt ( real ( dim_num, kind = rk ) )

  else

    norm = sqrt ( sum ( ( p(1:dim_num) - pc(1:dim_num) )**2 ) )
 
    pp(1:dim_num) = pc(1:dim_num) + r * ( p(1:dim_num) - pc(1:dim_num) ) / norm

  end if

  return
end
subroutine sphere_imp_volume_3d ( r, volume )

!*****************************************************************************80
!
!! SPHERE_IMP_VOLUME_3D computes the volume of an implicit sphere in 3D.
!
!  Discussion:
!
!    An implicit sphere in 3D satisfies the equation:
!
!      sum ( ( P(1:DIM_NUM) - pc(1:DIM_NUM) )^2 ) = R^2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) R, the radius of the sphere.
!
!    Output, real ( kind = rk ) VOLUME, the volume of the sphere.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) volume

  volume = ( 4.0D+00 / 3.0D+00 ) * r8_pi * r * r * r

  return
end
subroutine sphere_imp_volume_nd ( dim_num, r, volume )

!*****************************************************************************80
!
!! SPHERE_IMP_VOLUME_ND computes the volume of an implicit sphere in ND.
!
!  Discussion:
!
!    An implicit sphere in ND satisfies the equation:
!
!      sum ( ( X(1:N) - PC(1:N) )^2 ) = R^2
!
!    where R is the radius and PC is the center.
!
!    Results for the first few values of N are:
!
!    DIM_NUM  Volume
!    -     -----------------------
!    2                PI   * R^2
!    3     (4/3)    * PI   * R^3
!    4     (1/2)    * PI^2 * R^4
!    5     (8/15)   * PI^2 * R^5
!    6     (1/6)    * PI^3 * R^6
!    7     (16/105) * PI^3 * R^7
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the dimension of the space.
!
!    Input, real ( kind = rk ) R, the radius of the sphere.
!
!    Output, real ( kind = rk ) VOLUME, the volume of the sphere.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num
  real ( kind = rk ) r
  real ( kind = rk ) sphere01_volume_nd
  real ( kind = rk ) volume

  volume = r**dim_num * sphere01_volume_nd ( dim_num )

  return
end
subroutine sphere_imp_zone_area_3d ( r, h1, h2, area  )

!*****************************************************************************80
!
!! SPHERE_IMP_ZONE_AREA_3D computes the surface area of a spherical zone in 3D.
!
!  Discussion:
!
!    An implicit sphere in 3D satisfies the equation:
!
!      sum ( ( P(1:DIM_NUM) - PC(1:DIM_NUM) )^2 ) = R^2
!
!    Draw any radius of the sphere and note the point P where the radius
!    intersects the sphere.  Now choose two points on the radius line, a
!    distance H1 and H2 from the point P.  Consider all the points on or within
!    the sphere whose projection onto the radius lies between these two points.
!    These points constitute the spherical zone, which can also be considered
!    the difference of two spherical caps.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) R, the radius of the sphere.
!
!    Input, real ( kind = rk ) H1, H2, the distances that define the 
!    thickness of the zone.  H1 and H2 must be between 0 and 2 * R.
!
!    Output, real ( kind = rk ) AREA, the area of the spherical zone.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) area
  real ( kind = rk ) h
  real ( kind = rk ) h1
  real ( kind = rk ) h2
  real ( kind = rk ) r
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00

  h = abs ( h1 - h2 )

  if ( h <= 0.0D+00 ) then
    area = 0.0D+00
  else if ( 2.0D+00 * r <= h ) then
    area = 4.0D+00 * r8_pi * r * r
  else
    area = 2.0D+00 * r8_pi * r * h
  end if

  return
end
subroutine sphere_imp_zone_volume_3d ( r, h1, h2, volume )

!*****************************************************************************80
!
!! SPHERE_IMP_ZONE_VOLUME_3D computes the volume of a spherical zone in 3D.
!
!  Discussion:
!
!    An implicit sphere in 3D satisfies the equation:
!
!      sum ( ( P(1:DIM_NUM) - pc(1:DIM_NUM) )^2 ) = R^2
!
!    Draw any radius of the sphere and note the point P where the radius
!    intersects the sphere.  Now choose two points on the radius line, a
!    distance H1 and H2 from the point P.  Consider all the points on or within
!    the sphere whose projection onto the radius lies between these two points.
!    These points constitute the spherical zone, which can also be considered
!    the difference of two spherical caps.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) R, the radius of the sphere.
!
!    Input, real ( kind = rk ) H1, H2, the distances that define the 
!    thickness of the zone.  H1 and H2 must be between 0 and 2 * R.
!
!    Output, real ( kind = rk ) VOLUME, the volume of the spherical zone
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) h1
  real ( kind = rk ) h11
  real ( kind = rk ) h2
  real ( kind = rk ) h22
  real ( kind = rk ) r
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) volume

  h11 = min ( h1, h2 )
  h11 = max ( h11, 0.0D+00 )

  if ( 2.0D+00 * r <= h11 ) then
    volume = 0.0D+00
    return
  end if

  h22 = max ( h1, h2 )
  h22 = min ( h22, 2.0D+00 * r )

  if ( h22 <= 0.0D+00 ) then
    volume = 0.0D+00
    return
  end if

  volume = ( 1.0D+00 / 3.0D+00 ) * r8_pi * ( &
      h22 * h22 * ( 3.0D+00 * r - h22 ) &
    - h11 * h11 * ( 3.0D+00 * r - h11 ) )

  return
end
subroutine sphere_imp2exp_3d ( r, pc, p1, p2, p3, p4 )

!*****************************************************************************80
!
!! SPHERE_IMP2EXP_3D converts a sphere from implicit to explicit form in 3D.
!
!  Discussion:
!
!    An implicit sphere in 3D satisfies the equation:
!
!      sum ( ( P(1:DIM_NUM) - PC(1:DIM_NUM) )^2 ) = R^2
!
!    An explicit sphere in 3D is determined by four points,
!    which should be distinct, and not coplanar.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 February 2005
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
!  Parameters:
!
!    Input, real ( kind = rk ) R, PC(3), the radius and center of the sphere.
!
!    Output, real ( kind = rk ) P1(3), P2(3), P3(3), P4(3),
!    four distinct noncoplanar points on the sphere.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) p3(dim_num)
  real ( kind = rk ) p4(dim_num)
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) phi
  real ( kind = rk ) r
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) theta

  theta = 0.0D+00
  phi = 0.0D+00

  p1(1) = pc(1) + r * cos ( theta ) * sin ( phi )
  p1(2) = pc(2) + r * sin ( theta ) * sin ( phi )
  p1(3) = pc(3) + r                 * cos ( phi )

  theta = 0.0D+00
  phi = 2.0D+00 * r8_pi / 3.0D+00

  p2(1) = pc(1) + r * cos ( theta ) * sin ( phi )
  p2(2) = pc(2) + r * sin ( theta ) * sin ( phi )
  p2(3) = pc(3) + r                 * cos ( phi )

  theta = 2.0D+00 * r8_pi / 3.0D+00
  phi = 2.0D+00 * r8_pi / 3.0D+00

  p3(1) = pc(1) + r * cos ( theta ) * sin ( phi )
  p3(2) = pc(2) + r * sin ( theta ) * sin ( phi )
  p3(3) = pc(3) + r                 * cos ( phi )

  theta = 4.0D+00 * r8_pi / 3.0D+00
  phi = 2.0D+00 * r8_pi / 3.0D+00

  p4(1) = pc(1) + r * cos ( theta ) * sin ( phi )
  p4(2) = pc(2) + r * sin ( theta ) * sin ( phi )
  p4(3) = pc(3) + r                 * cos ( phi )

  return
end
function sphere_k ( dim_num )

!*****************************************************************************80
!
!! SPHERE_K computes a factor useful for spherical computations.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Thomas Ericson, Victor Zinoviev,
!    Codes on Euclidean Spheres,
!    Elsevier, 2001, pages 439-441.
!    QA166.7 E75
!
!  Parameters:
!
!    Input, integer DIM_NUM, the dimension of the space.
!
!    Output, real ( kind = rk ) SPHERE_K, the factor.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i4_factorial2
  integer dim_num
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) sphere_k

  if ( mod ( dim_num, 2 ) == 0 ) then
    sphere_k = ( 2.0D+00 * r8_pi ) ** ( dim_num / 2 )
  else
    sphere_k = 2.0D+00 * ( 2.0D+00 * r8_pi ) ** ( ( dim_num - 1 ) / 2 )
  end if

  sphere_k = sphere_k / real ( i4_factorial2 ( dim_num - 2 ), kind = rk )

  return
end
subroutine sphere_triangle_angles_to_area ( r, a, b, c, area )

!*****************************************************************************80
!
!! SPHERE_TRIANGLE_ANGLES_TO_AREA computes the area of a spherical triangle.
!
!  Discussion:
!
!    A sphere centered at 0 in 3D satisfies the equation:
!
!      X*X + Y*Y + Z*Z = R*R
!
!    A spherical triangle is specified by three points on the surface
!    of the sphere.
!
!    The area formula is known as Girard's formula.
!
!    The area of a spherical triangle is:
!
!      AREA = ( A + B + C - PI ) * R*R
!
!    where A, B and C are the (surface) angles of the triangle.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) R, the radius of the sphere.
!
!    Input, real ( kind = rk ) A, B, C, the angles of the triangle.
!
!    Output, real ( kind = rk ) AREA, the area of the spherical triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) area
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) r
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
!
!  Apply Girard's formula.
!
  area = r * r * ( a + b + c - r8_pi )

  return
end
subroutine sphere_triangle_sides_to_angles ( r, as, bs, cs, a, b, c )

!*****************************************************************************80
!
!! SPHERE_TRIANGLE_SIDES_TO_ANGLES computes spherical triangle angles.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) R, the radius of the sphere.
!
!    Input, real ( kind = rk ) AS, BS, CS, the (geodesic) length of the 
!    sides of the triangle.
!
!    Output, real ( kind = rk ) A, B, C, the spherical angles of the triangle.
!    Angle A is opposite the side of length AS, and so on.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) as
  real ( kind = rk ) asu
  real ( kind = rk ) b
  real ( kind = rk ) bs
  real ( kind = rk ) bsu
  real ( kind = rk ) c
  real ( kind = rk ) cs
  real ( kind = rk ) csu
  real ( kind = rk ) r
  real ( kind = rk ) ssu
  real ( kind = rk ) tan_a2
  real ( kind = rk ) tan_b2
  real ( kind = rk ) tan_c2

  asu = as / r
  bsu = bs / r
  csu = cs / r
  ssu = ( asu + bsu + csu ) / 2.0D+00

  tan_a2 = sqrt ( ( sin ( ssu - bsu ) * sin ( ssu - csu ) ) / & 
                  ( sin ( ssu ) * sin ( ssu - asu )     ) )

  a = 2.0D+00 * atan ( tan_a2 )

  tan_b2 = sqrt ( ( sin ( ssu - asu ) * sin ( ssu - csu ) ) / & 
                  ( sin ( ssu ) * sin ( ssu - bsu )     ) )

  b = 2.0D+00 * atan ( tan_b2 )

  tan_c2 = sqrt ( ( sin ( ssu - asu ) * sin ( ssu - bsu ) ) / & 
                  ( sin ( ssu ) * sin ( ssu - csu )     ) )

  c = 2.0D+00 * atan ( tan_c2 )

  return
end
subroutine sphere_triangle_vertices_to_angles ( r, v1, v2, v3, a, b, c )

!*****************************************************************************80
!
!! SPHERE_TRIANGLE_VERTICES_TO_ANGLES: spherical triangle angles from vertices.
!
!  Discussion:
!
!    A sphere centered at 0 in 3D satisfies the equation:
!
!      X * X + Y * Y + Z * Z = R * R
!
!    A spherical triangle is specified by three points on the surface
!    of the sphere.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    24 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) R, the radius of the sphere.
!
!    Input, real ( kind = rk ) V1(3), V2(3), V3(3), the vertices of the triangle.
!
!    Output, real ( kind = rk ) A, B, C, the angles of the spherical triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) as
  real ( kind = rk ) b
  real ( kind = rk ) bs
  real ( kind = rk ) c
  real ( kind = rk ) cs
  real ( kind = rk ) r
  real ( kind = rk ) v1(3)
  real ( kind = rk ) v2(3)
  real ( kind = rk ) v3(3)
!
!  Compute the lengths of the sides of the spherical triangle.
!
  call sphere_triangle_vertices_to_sides ( r, v1, v2, v3, as, bs, cs )
!
!  Get the spherical angles.
!
  call sphere_triangle_sides_to_angles ( r, as, bs, cs, a, b, c )

  return
end
subroutine sphere_triangle_vertices_to_area ( r, v1, v2, v3, area )

!*****************************************************************************80
!
!! SPHERE_TRIANGLE_VERTICES_TO_AREA computes the area of a spherical triangle.
!
!  Discussion:
!
!    A sphere centered at 0 in 3D satisfies the equation:
!
!      X * X + Y * Y + Z * Z = R * R
!
!    A spherical triangle is specified by three points on the surface
!    of the sphere.
!
!    The area formula is known as Girard's formula.
!
!    The area of a spherical triangle is:
!
!      AREA = ( A + B + C - PI ) * R*R
!
!    where A, B and C are the (surface) angles of the triangle.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) R, the radius of the sphere.
!
!    Input, real ( kind = rk ) V1(3), V2(3), V3(3), the vertices of the triangle.
!
!    Output, real ( kind = rk ) AREA, the area of the spherical triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) area
  real ( kind = rk ) as
  real ( kind = rk ) b
  real ( kind = rk ) bs
  real ( kind = rk ) c
  real ( kind = rk ) cs
  real ( kind = rk ) r
  real ( kind = rk ) v1(3)
  real ( kind = rk ) v2(3)
  real ( kind = rk ) v3(3)
!
!  Compute the lengths of the sides of the spherical triangle.
!
  call sphere_triangle_vertices_to_sides ( r, v1, v2, v3, as, bs, cs )
!
!  Get the spherical angles.
!
  call sphere_triangle_sides_to_angles ( r, as, bs, cs, a, b, c )
!
!  Get the area.
!
  call sphere_triangle_angles_to_area ( r, a, b, c, area )

  return
end
subroutine sphere_triangle_vertices_to_centroid ( r, v1, v2, v3, vs )

!*****************************************************************************80
!
!! SPHERE_TRIANGLE_VERTICES_TO_CENTROID gets a spherical triangle centroid.
!
!  Discussion:
!
!    A sphere centered at 0 in 3D satisfies the equation:
!
!      X*X + Y*Y + Z*Z = R*R
!
!    A spherical triangle is specified by three points on the sphere.
!
!    The (true) centroid of a spherical triangle is the point
!
!      VT = (XT,YT,ZT) = Integral ( X, Y, Z ) dArea / Integral 1 dArea
!
!    Note that the true centroid does NOT, in general, lie on the sphere.  
!
!    The "flat" centroid VF is the centroid of the planar triangle defined by
!    the vertices of the spherical triangle.
!
!    The "spherical" centroid VS of a spherical triangle is computed by
!    the intersection of the geodesic bisectors of the triangle angles.
!    The spherical centroid lies on the sphere.
!
!    VF, VT and VS lie on a line through the center of the sphere.  We can
!    easily calculate VF by averaging the vertices, and from this determine
!    VS by normalizing.
!
!    Of course, we still will not have actually computed VT, which lies
!    somewhere between VF and VS!
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) R, the radius of the sphere.
!
!    Input, real ( kind = rk ) V1(3), V2(3), V3(3), the vertices of the triangle.
!
!    Output, real ( kind = rk ) VS(3), the coordinates of the "spherical
!    centroid" of the spherical triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) norm
  real ( kind = rk ) r
  real ( kind = rk ) r8vec_norm
  real ( kind = rk ) v1(3)
  real ( kind = rk ) v2(3)
  real ( kind = rk ) v3(3)
  real ( kind = rk ) vs(3)

  vs(1:3) = ( v1(1:3) + v2(1:3) + v3(1:3) ) / 3.0D+00

  norm = r8vec_norm ( 3, vs )

  vs(1:3) = r * vs(1:3) / norm

  return
end
subroutine sphere_triangle_vertices_to_orientation ( a, b, c, o )

!*****************************************************************************80
!
!! SPHERE_TRIANGLE_VERTICES_TO_ORIENTATION: orientation of a spherical triangle.
!
!  Discussion:
!
!    Three points on a sphere actually compute two triangles; typically
!    we are interested in the smaller of the two.
!
!    As long as our triangle is "small", we can define an orientation
!    by comparing the direction of the centroid against the normal
!    vector (C-B) x (A-B).  If the dot product of these vectors
!    is positive, we say the triangle has positive orientation.
!
!    By using information from the triangle orientation, we can correctly
!    determine the area of a Voronoi polygon by summing up the pieces
!    of Delaunay triangles, even in the case when the Voronoi vertex
!    lies outside the Delaunay triangle.  In that case, the areas of
!    some of the Delaunay triangle pieces must be formally negative.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 May 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A(3), B(3), C(3), three points on a sphere.
!
!    Output, integer O, is +1 if the spherical triangle is 
!    judged to have positive orientation, and -1 otherwise.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(3)
  real ( kind = rk ) b(3)
  real ( kind = rk ) c(3)
  real ( kind = rk ) cd(3)
  real ( kind = rk ) cp(3)
  integer o
  real ( kind = rk ) v1(3)
  real ( kind = rk ) v2(3)
!
!  Centroid.
!
  cd(1:3) = ( a(1:3) + b(1:3) + c(1:3) ) / 3.0D+00
!
!  Cross product ( C - B ) x ( A - B );
!
  v1(1:3) = c(1:3) - b(1:3)
  v2(1:3) = a(1:3) - b(1:3)

  cp(1) = v1(2) * v2(3) - v1(3) * v2(2)
  cp(2) = v1(3) * v2(1) - v1(1) * v2(3)
  cp(3) = v1(1) * v2(2) - v1(2) * v2(1)
!
!  Compare the directions.
!
  if ( dot_product ( cp, cd ) < 0.0D+00 ) then
    o = - 1
  else
    o = + 1
  end if

  return
end
subroutine sphere_triangle_vertices_to_sides ( r, v1, v2, v3, as, bs, cs )

!*****************************************************************************80
!
!! SPHERE_TRIANGLE_VERTICES_TO_SIDES computes spherical triangle sides.
!
!  Discussion:
!
!    We can use the ACOS system call here, but the ARC_COSINE routine
!    will automatically take care of cases where the input argument is
!    (usually slightly) out of bounds.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) R, the radius of the sphere.
!
!    Input, real ( kind = rk ) V1(3), V2(3), V3(3), the vertices of the spherical
!    triangle.
!
!    Output, real ( kind = rk ) AS, BS, CS, the (geodesic) length of the sides
!    of the triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) as
  real ( kind = rk ) bs
  real ( kind = rk ) cs
  real ( kind = rk ) r
  real ( kind = rk ) r8_acos
  real ( kind = rk ) v1(3)
  real ( kind = rk ) v2(3)
  real ( kind = rk ) v3(3)

  as = r * r8_acos ( dot_product ( v2(1:3), v3(1:3) ) / r**2 )
  bs = r * r8_acos ( dot_product ( v3(1:3), v1(1:3) ) / r**2 )
  cs = r * r8_acos ( dot_product ( v1(1:3), v2(1:3) ) / r**2 )

  return
end
function sphere01_area_nd ( dim_num )

!*****************************************************************************80
!
!! SPHERE01_AREA_ND computes the surface area of a unit sphere in ND.
!
!  Discussion:
!
!    The unit sphere in ND satisfies:
!
!      sum ( 1 <= I <= DIM_NUM ) X(I) * X(I) = 1
!
!    Results for the first few values of N are:
!
!    DIM_NUM   Area
!
!     2    2        * PI
!     3    4        * PI
!     4  ( 2 /   1) * PI^2
!     5  ( 8 /   3) * PI^2
!     6  ( 1 /   1) * PI^3
!     7  (16 /  15) * PI^3
!     8  ( 1 /   3) * PI^4
!     9  (32 / 105) * PI^4
!    10  ( 1 /  12) * PI^5
!
!    For the unit sphere, Area(DIM_NUM) = DIM_NUM * Volume(DIM_NUM)
!
!    Sphere_Unit_Area ( DIM_NUM ) = 2 * PI^(DIM_NUM/2) / Gamma ( DIM_NUM / 2 )
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the dimension of the space.
!
!    Output, real ( kind = rk ) SPHERE01_AREA_ND, the area of the sphere.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) area
  integer dim_num
  integer i
  integer m
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) sphere01_area_nd

  if ( mod ( dim_num, 2 ) == 0 ) then
    m = dim_num / 2
    area = 2.0D+00 * ( r8_pi ) ** m
    do i = 1, m - 1
      area = area / real ( i, kind = rk )
    end do
  else
    m = ( dim_num - 1 ) / 2
    area = ( r8_pi ) ** m * 2.0D+00 ** dim_num
    do i = m + 1, 2 * m
      area = area / real ( i,  kind = rk )
    end do
  end if

  sphere01_area_nd = area

  return
end
subroutine sphere01_area_values ( n_data, n, area )

!*****************************************************************************80
!
!! SPHERE01_AREA_VALUES returns some areas of the unit sphere in ND.
!
!  Discussion:
!
!    The formula for the surface area of the unit sphere in N dimensions is:
!
!      Sphere_Unit_Area ( N ) = 2 * pi^(N/2) / Gamma ( N / 2 )
!
!    Some values of the function include:
!
!       N   Area
!
!       2    2        * PI
!       3  ( 4 /    ) * PI
!       4  ( 2 /   1) * PI^2
!       5  ( 8 /   3) * PI^2
!       6  ( 1 /   1) * PI^3
!       7  (16 /  15) * PI^3
!       8  ( 1 /   3) * PI^4
!       9  (32 / 105) * PI^4
!      10  ( 1 /  12) * PI^5
!
!    For the unit sphere, Area(N) = N * Volume(N)
!
!    In Mathematica, the function can be evaluated by:
!
!      2 * Pi^(n/2) / Gamma[n/2]
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer N, the spatial dimension.
!
!    Output, real ( kind = rk ) AREA, the area of the unit sphere
!    in that dimension.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 20

  real ( kind = rk ) area
  real ( kind = rk ), save, dimension ( n_max ) :: area_vec = (/ &
    0.2000000000000000D+01, &
    0.6283185307179586D+01, &
    0.1256637061435917D+02, &
    0.1973920880217872D+02, &
    0.2631894506957162D+02, &
    0.3100627668029982D+02, &
    0.3307336179231981D+02, &
    0.3246969701133415D+02, &
    0.2968658012464836D+02, &
    0.2550164039877345D+02, &
    0.2072514267328890D+02, &
    0.1602315322625507D+02, &
    0.1183817381218268D+02, &
    0.8389703410491089D+01, &
    0.5721649212349567D+01, &
    0.3765290085742291D+01, &
    0.2396678817591364D+01, &
    0.1478625959000308D+01, &
    0.8858104195716824D+00, &
    0.5161378278002812D+00 /)
  integer n_data
  integer n
  integer, save, dimension ( n_max ) :: n_vec = (/ &
     1, &
     2, &
     3, &
     4, &
     5, &
     6, &
     7, &
     8, &
     9, &
    10, &
    11, &
    12, &
    13, &
    14, &
    15, &
    16, &
    17, &
    18, &
    19, &
    20 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    area = 0.0D+00
  else
    n = n_vec(n_data)
    area = area_vec(n_data)
  end if

  return
end
subroutine sphere01_sample_2d ( x )

!*****************************************************************************80
!
!! SPHERE01_SAMPLE_2D picks a random point on the unit sphere (circle) in 2D.
!
!  Discussion:
!
!    The unit sphere in 2D satisfies:
!
!      X * X + Y * Y = 1
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) X(2), a random point on the unit circle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) u
  real ( kind = rk ) x(dim_num)

  call random_number ( harvest = u )

  x(1) = cos ( 2.0D+00 * r8_pi * u )
  x(2) = sin ( 2.0D+00 * r8_pi * u )

  return
end
subroutine sphere01_sample_3d ( x )

!*****************************************************************************80
!
!! SPHERE01_SAMPLE_3D picks a random point on the unit sphere in 3D.
!
!  Discussion:
!
!    The unit sphere in 3D satisfies:
!
!      X * X + Y * Y + Z * Z = 1
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) X(3), the sample point.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) r8_acos
  real ( kind = rk ) phi
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) theta
  real ( kind = rk ) vdot
  real ( kind = rk ) x(dim_num)
!
!  Pick a uniformly random VDOT, which must be between -1 and 1.
!  This represents the dot product of the random vector with the Z unit vector.
!
!  This works because the surface area of the sphere between
!  Z and Z + dZ is independent of Z.  So choosing Z uniformly chooses
!  a patch of area uniformly.
!
  call random_number ( harvest = vdot )
  vdot = 2.0D+00 * vdot - 1.0D+00

  phi = r8_acos ( vdot )
!
!  Pick a uniformly random rotation between 0 and 2 Pi around the
!  axis of the Z vector.
!
  call random_number ( harvest = theta )
  theta = 2.0D+00 * r8_pi * theta

  x(1) = cos ( theta ) * sin ( phi )
  x(2) = sin ( theta ) * sin ( phi )
  x(3) = cos ( phi )

  return
end
subroutine sphere01_sample_3d_2 ( x )

!*****************************************************************************80
!
!! SPHERE01_SAMPLE_3D_2 is a BAD method for sampling the unit sphere in 3D.
!
!  Discussion:
!
!    The unit sphere in 3D satisfies:
!
!      X * X + Y * Y + Z * Z = 1
!
!    Points on the unit sphere have coordinates ( PHI, THETA ) where
!    PHI varies from 0 to PI, and THETA from 0 to 2 PI, so that:
!
!    x = cos ( theta ) * sin ( phi )
!    y = sin ( theta ) * sin ( phi )
!    z =                 cos ( phi )
!
!    This routine implements a sampling of the sphere that simply
!    picks PHI and THETA uniformly at random from their ranges.
!    This is a uniform sampling on the cylinder, but it is NOT
!    a uniform sampling on the sphere.  I implement it here just
!    so I can run some tests against the code in SPHERE_UNIT_SAMPLE_3D.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) X(3), the sample point.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) phi
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) theta
  real ( kind = rk ) x(dim_num)

  call random_number ( harvest = phi )
  phi = r8_pi * phi

  call random_number ( harvest = theta )
  theta = 2.0D+00 * r8_pi * theta

  x(1) = cos ( theta ) * sin ( phi )
  x(2) = sin ( theta ) * sin ( phi )
  x(3) = cos ( phi )

  return
end
subroutine sphere01_sample_nd ( dim_num, x )

!*****************************************************************************80
!
!! SPHERE01_SAMPLE_ND picks a random point on the unit sphere in ND.
!
!  Discussion:
!
!    The unit sphere in ND satisfies:
!
!      sum ( 1 <= I <= DIM_NUM ) X(I) * X(I) = 1
!
!    DIM_NUM-1 random Givens rotations are applied to the point 
!    ( 1, 0, 0, ..., 0 ).
!
!    The I-th Givens rotation is in the plane of coordinate axes I and I+1,
!    and has the form:
!
!     [ cos ( theta )  - sin ( theta ) ] * x(i)      = x'(i)
!     [ sin ( theta )    cos ( theta ) ]   x(i+1)      x'(i+1)
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = rk ) X(DIM_NUM), the random point.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num

  integer i
  real ( kind = rk ) random_cosine
  real ( kind = rk ) random_sign
  real ( kind = rk ) random_sine
  real ( kind = rk ) x(dim_num)
  real ( kind = rk ) xi

  x(1) = 1.0D+00
  x(2:dim_num) = 0.0D+00

  do i = 1, dim_num-1
    call random_number ( harvest = random_cosine )
    random_cosine = 2.0D+00 * random_cosine - 1.0D+00
    call random_number ( harvest = random_sign )
    random_sign = real ( 2 * int ( 2.0D+00 * random_sign ) - 1,  kind = rk )
    random_sine = random_sign * sqrt ( 1.0D+00 - random_cosine**2 )
    xi = x(i)
    x(i  ) = random_cosine * xi
    x(i+1) = random_sine   * xi
  end do

  return
end
subroutine sphere01_sample_nd_2 ( dim_num, x )

!*****************************************************************************80
!
!! SPHERE01_SAMPLE_ND_2 picks a random point on the unit sphere in ND.
!
!  Discussion:
!
!    The unit sphere in ND satisfies:
!
!      sum ( 1 <= I <= DIM_NUM ) X(I) * X(I) = 1
!
!    DIM_NUM independent normally distributed random numbers are generated,
!    and then scaled to have unit norm.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = rk ) X(DIM_NUM), the random point.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num

  real ( kind = rk ) norm
  real ( kind = rk ) x(dim_num)

  call r8vec_normal_01 ( dim_num, x )

  norm = sqrt ( sum ( x(1:dim_num)**2 ) )

  x(1:dim_num) = x(1:dim_num) / norm

  return
end
subroutine sphere01_sample_nd_3 ( dim_num, x )

!*****************************************************************************80
!
!! SPHERE01_SAMPLE_ND_3 picks a random point on the unit sphere in ND.
!
!  Discussion:
!
!    The unit sphere in ND satisfies:
!
!      sum ( 1 <= I <= DIM_NUM ) X(I) * X(I) = 1
!
!    Points in the [-1,1] cube are generated.  Points lying outside
!    the sphere are rejected.  Points inside the unit sphere are normalized
!    to lie on the sphere.
!
!    Because the volume of the unit sphere
!    relative to the unit cube decreases drastically in higher dimensions,
!    this routine becomes increasingly inefficient at higher DIM_NUM.  
!    Above DIM_NUM = 5, this problem will become significant.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = rk ) X(DIM_NUM), the random point.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num

  real ( kind = rk ) norm
  real ( kind = rk ) x(dim_num)

  do

    call random_number ( harvest = x(1:dim_num) )

    x(1:dim_num) = 2.0D+00 * x(1:dim_num) - 1.0D+00

    norm = sqrt ( sum ( x(1:dim_num)**2 ) )

    if ( norm <= 1.0D+00 ) then
      x(1:dim_num) = x(1:dim_num) / norm
      exit
    end if

  end do

  return
end
function sphere01_volume_nd ( dim_num )

!*****************************************************************************80
!
!! SPHERE01_VOLUME_ND computes the volume of a unit sphere in ND.
!
!  Discussion:
!
!    The unit sphere in ND satisfies:
!
!      sum ( 1 <= I <= DIM_NUM ) X(I) * X(I) = 1
!
!    Results for the first few values of DIM_NUM are:
!
!     DIM_NUM  Volume
!
!     1    2
!     2    1        * PI
!     3  ( 4 /   3) * PI
!     4  ( 1 /   2) * PI^2
!     5  ( 8 /  15) * PI^2
!     6  ( 1 /   6) * PI^3
!     7  (16 / 105) * PI^3
!     8  ( 1 /  24) * PI^4
!     9  (32 / 945) * PI^4
!    10  ( 1 / 120) * PI^5
!
!    For the unit sphere, Volume(DIM_NUM) = 2 * PI * Volume(DIM_NUM-2)/ DIM_NUM
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Output, real ( kind = rk ) SPHERE_UNIT_VOLUME_ND, the volume of the sphere.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num
  integer i
  integer m
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) sphere01_volume_nd
  real ( kind = rk ) volume

  if ( mod ( dim_num, 2 ) == 0 ) then
    m = dim_num / 2
    volume = r8_pi ** m
    do i = 1, m
      volume = volume / real ( i, kind = rk )
    end do
  else
    m = ( dim_num - 1 ) / 2
    volume = r8_pi ** m * 2.0D+00**dim_num
    do i = m + 1, 2 * m + 1
      volume = volume / real ( i, kind = rk )
    end do
  end if

  sphere01_volume_nd = volume

  return
end
subroutine sphere01_volume_values ( n_data, n, volume )

!*****************************************************************************80
!
!! SPHERE01_VOLUME_VALUES returns some volumes of the unit sphere in ND.
!
!  Discussion:
!
!    The formula for the volume of the unit sphere in N dimensions is
!
!      Volume(N) = 2 * pi^(N/2) / ( N * Gamma ( N / 2 ) )
!
!    This function satisfies the relationships:
!
!      Volume(N) = 2 * pi * Volume(N-2) / N
!      Volume(N) = Area(N) / N
!
!    Some values of the function include:
!
!       N  Volume
!
!       1    1
!       2    1        * PI
!       3  ( 4 /   3) * PI
!       4  ( 1 /   2) * PI^2
!       5  ( 8 /  15) * PI^2
!       6  ( 1 /   6) * PI^3
!       7  (16 / 105) * PI^3
!       8  ( 1 /  24) * PI^4
!       9  (32 / 945) * PI^4
!      10  ( 1 / 120) * PI^5
!
!    In Mathematica, the function can be evaluated by:
!
!      2 * Pi^(n/2) / ( n * Gamma[n/2] )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer N, the spatial dimension.
!
!    Output, real ( kind = rk ) VOLUME, the volume of the unit
!    sphere in that dimension.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 20

  integer n_data
  integer n
  integer, save, dimension ( n_max ) :: n_vec = (/ &
     1,  2, &
     3,  4, &
     5,  6, &
     7,  8, &
     9, 10, &
    11, 12, &
    13, 14, &
    15, 16, &
    17, 18, &
    19, 20 /)
  real ( kind = rk ) volume
  real ( kind = rk ), save, dimension ( n_max ) :: volume_vec = (/ &
    0.2000000000000000D+01, &
    0.3141592653589793D+01, &
    0.4188790204786391D+01, &
    0.4934802200544679D+01, &
    0.5263789013914325D+01, &
    0.5167712780049970D+01, &
    0.4724765970331401D+01, &
    0.4058712126416768D+01, &
    0.3298508902738707D+01, &
    0.2550164039877345D+01, &
    0.1884103879389900D+01, &
    0.1335262768854589D+01, &
    0.9106287547832831D+00, &
    0.5992645293207921D+00, &
    0.3814432808233045D+00, &
    0.2353306303588932D+00, &
    0.1409811069171390D+00, &
    0.8214588661112823D-01, &
    0.4662160103008855D-01, &
    0.2580689139001406D-01 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    volume = 0.0D+00
  else
    n = n_vec(n_data)
    volume = volume_vec(n_data)
  end if

  return
end
subroutine sphere01_distance_xyz ( xyz1, xyz2, dist )

!*****************************************************************************80
!
!! SPHERE01_DISTANCE_XYZ computes great circle distances on a unit sphere.
!
!  Discussion:
!
!    XYZ coordinates are used.
!
!    We assume the points XYZ1 and XYZ2 lie on the unit sphere.
!
!    This computation is a special form of the Vincenty formula.
!    It should be less sensitive to errors associated with very small 
!    or very large angular separations.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    "Great-circle distance",
!    Wikipedia.
!
!  Parameters:
!
!    Input, real ( kind = rk ) XYZ1(3), the coordinates of the first point.
!
!    Input, real ( kind = rk ) XYZ2(3), the coordinates of the second point.
!
!    Output, real ( kind = rk ) DIST, the great circle distance between
!    the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) bot
  real ( kind = rk ) dist
  real ( kind = rk ) lat1
  real ( kind = rk ) lat2
  real ( kind = rk ) lon1
  real ( kind = rk ) lon2
  real ( kind = rk ) r8_asin
  real ( kind = rk ) r8_atan
  real ( kind = rk ) top
  real ( kind = rk ) xyz1(3)
  real ( kind = rk ) xyz2(3)

  lat1 = r8_asin ( xyz1(3) )
  lon1 = r8_atan ( xyz1(2), xyz1(1) )

  lat2 = r8_asin ( xyz2(3) )
  lon2 = r8_atan ( xyz2(2), xyz2(1) )

  top = ( cos ( lat2 ) * sin ( lon1 - lon2 ) )**2 &
      + ( cos ( lat1 ) * sin ( lat2 ) &
      -   sin ( lat1 ) * cos ( lat2 ) * cos ( lon1 - lon2 ) )**2

  top = sqrt ( top )

  bot = sin ( lat1 ) * sin ( lat2 ) &
      + cos ( lat1 ) * cos ( lat2 ) * cos ( lon1 - lon2 )

  dist = atan2 ( top, bot )

  return
end
function sphere01_polygon_area ( n, lat, lon )

!*****************************************************************************80
!
!! sphere01_polygon_area() returns the area of a spherical polygon.
!
!  Discussion:
!
!    On a unit sphere, the area of a spherical polygon with N sides
!    is equal to the spherical excess:
!
!      E = sum ( interior angles ) - ( N - 2 ) * pi.
!
!    On a sphere with radius R, the area is the spherical excess multiplied
!    by R * R.
!
!    The code was revised in accordance with suggestions in Carvalho and
!    Cavalcanti.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 August 2005
!
!  Author:
!
!    Original C version by Robert Miller.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paulo Cezar Pinto Carvalho, Paulo Roma Cavalcanti,
!    Point in Polyhedron Testing Using Spherical Polygons,
!    in Graphics Gems V,
!    edited by Alan Paeth,
!    Academic Press, 1995,
!    ISBN: 0125434553,
!    LC: T385.G6975.
!
!    Robert Miller,
!    Computing the Area of a Spherical Polygon,
!    Graphics Gems, Volume IV, pages 132-138,
!    Edited by Paul Heckbert,
!    Academic Press, 1994, T385.G6974.
!
!    Eric Weisstein,
!    "Spherical Polygon",
!    CRC Concise Encyclopedia of Mathematics,
!    CRC Press, 1999.
!
!  Parameters:
!
!    Input, integer N, the number of vertices.
!
!    Input, real ( kind = rk ) LAT[N], LON[N], the latitudes and longitudes 
!    of the vertices of the spherical polygon.
!
!    Output, real ( kind = rk ) SPHERE01_POLYGON_AREA, the area of the 
!    spherical polygon, measured in spherical radians.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a
  real ( kind = rk ) area
  real ( kind = rk ) b
  real ( kind = rk ) beta1
  real ( kind = rk ) beta2
  real ( kind = rk ) c
  real ( kind = rk ) cos_b1
  real ( kind = rk ) cos_b2
  real ( kind = rk ) excess
  real ( kind = rk ) hav_a
  real ( kind = rk ) haversine
  integer j
  integer k
  real ( kind = rk ) lam
  real ( kind = rk ) lam1
  real ( kind = rk ) lam2
  real ( kind = rk ) lat(n)
  real ( kind = rk ) lon(n)
  real ( kind = rk ), parameter :: r8_pi_half = 1.5707963267948966192313D+00
  real ( kind = rk ) s
  real ( kind = rk ) sphere01_polygon_area
  real ( kind = rk ) t

  area = 0.0D+00

  do j = 1, n + 1

    if ( j == 1 ) then
      lam1 = lon(j)
      beta1 = lat(j)
      lam2 = lon(j+1)
      beta2 = lat(j+1)
      cos_b1 = cos ( beta1 )
      cos_b2 = cos ( beta2 )
    else
      k = mod ( j + 1, n + 1 )
      lam1 = lam2
      beta1 = beta2
      lam2 = lon(k)
      beta2 = lat(k)
      cos_b1 = cos_b2
      cos_b2 = cos ( beta2 )
    end if

    if ( lam1 /= lam2 ) then

      hav_a = haversine ( beta2 - beta1 ) &
        + cos_b1 * cos_b2 * haversine ( lam2 - lam1 )
      a = 2.0D+00 * asin ( sqrt ( hav_a ) )

      b = r8_pi_half - beta2
      c = r8_pi_half - beta1
      s = 0.5D+00 * ( a + b + c )
!
!  Given the three sides of a spherical triangle, we can use a formula
!  to find the spherical excess.
!
      t = tan ( s / 2.0D+00 ) * tan ( ( s - a ) / 2.0D+00 ) &
        * tan ( ( s - b ) / 2.0D+00 ) * tan ( ( s - c ) / 2.0D+00 )

      excess = abs ( 4.0D+00 * atan ( sqrt ( abs ( t ) ) ) )

      if ( lam1 < lam2 ) then
        lam = lam2 - lam1
      else
        lam = lam2 - lam1 + 4.0D+00 * r8_pi_half
      end if

      if ( 2.0D+00 * r8_pi_half < lam ) then
        excess = -excess 
      end if

      area = area + excess

    end if

  end do

  sphere01_polygon_area = abs ( area )

  return
end
subroutine sphere01_triangle_angles_to_area ( a, b, c, area )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_ANGLES_TO_AREA computes the area of a spherical triangle.
!
!  Discussion:
!
!    A unit sphere in 3D satisfies the equation:
!
!      X^2 + Y^2 + Z^2 = 1
!
!    A spherical triangle is specified by three points on the surface
!    of the sphere.
!
!    The area formula is known as Girard's formula.
!
!    The area of a spherical triangle on a unit sphere is:
!
!      AREA = ( A + B + C - PI )
!
!    where A, B and C are the (surface) angles of the triangle.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, C, the angles of the triangle.
!
!    Output, real ( kind = rk ) AREA, the area of the sphere.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) area
  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
!
!  Apply Girard's formula.
!
  area = a + b + c - r8_pi

  return
end
subroutine sphere01_triangle_sides_to_angles ( as, bs, cs, a, b, c )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_SIDES_TO_ANGLES computes spherical triangle angles.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) AS, BS, CS, the (geodesic) length of the 
!    sides of the triangle.
!
!    Output, real ( kind = rk ) A, B, C, the spherical angles of the triangle.
!    Angle A is opposite the side of length AS, and so on.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) as
  real ( kind = rk ) asu
  real ( kind = rk ) b
  real ( kind = rk ) bs
  real ( kind = rk ) bsu
  real ( kind = rk ) c
  real ( kind = rk ) cs
  real ( kind = rk ) csu
  real ( kind = rk ) ssu
  real ( kind = rk ) tan_a2
  real ( kind = rk ) tan_b2
  real ( kind = rk ) tan_c2

  asu = as
  bsu = bs
  csu = cs
  ssu = ( asu + bsu + csu ) / 2.0D+00

  tan_a2 = sqrt ( ( sin ( ssu - bsu ) * sin ( ssu - csu ) ) / &
                  ( sin ( ssu ) * sin ( ssu - asu )     ) )

  a = 2.0D+00 * atan ( tan_a2 )

  tan_b2 = sqrt ( ( sin ( ssu - asu ) * sin ( ssu - csu ) ) / &
                  ( sin ( ssu ) * sin ( ssu - bsu )     ) )

  b = 2.0D+00 * atan ( tan_b2 )

  tan_c2 = sqrt ( ( sin ( ssu - asu ) * sin ( ssu - bsu ) ) / &
                  ( sin ( ssu ) * sin ( ssu - csu )     ) )

  c = 2.0D+00 * atan ( tan_c2 )

  return
end
subroutine sphere01_triangle_vertices_to_angles ( v1, v2, v3, a, b, c )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_VERTICES_TO_ANGLES: spherical triangle angles by vertices.
!
!  Discussion:
!
!    A unit sphere centered at 0 in 3D satisfies the equation:
!
!      X * X + Y * Y + Z * Z = 1
!
!    A spherical triangle is specified by three points on the surface
!    of the sphere.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    25 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) V1(3), V2(3), V3(3), the vertices of the triangle.
!
!    Output, real ( kind = rk ) A, B, C, the angles of the spherical triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) as
  real ( kind = rk ) b
  real ( kind = rk ) bs
  real ( kind = rk ) c
  real ( kind = rk ) cs
  real ( kind = rk ) v1(3)
  real ( kind = rk ) v2(3)
  real ( kind = rk ) v3(3)
!
!  Compute the lengths of the sides of the spherical triangle.
!
  call sphere01_triangle_vertices_to_sides ( v1, v2, v3, as, bs, cs )
!
!  Get the spherical angles.
!
  call sphere01_triangle_sides_to_angles ( as, bs, cs, a, b, c )

  return
end
subroutine sphere01_triangle_vertices_to_area ( v1, v2, v3, area )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_VERTICES_TO_AREA computes the area of a spherical triangle.
!
!  Discussion:
!
!    A unit sphere in 3D satisfies the equation:
!
!      X^2 + Y^2 + Z^2 = 1
!
!    A spherical triangle is specified by three points on the surface
!    of the sphere.
!
!    The area formula is known as Girard's formula.
!
!    The area of a spherical triangle on a unit sphere is:
!
!      AREA = ( A + B + C - PI )
!
!    where A, B and C are the (surface) angles of the triangle.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) V1(3), V2(3), V3(3), the vertices of the triangle.
!
!    Output, real ( kind = rk ) AREA, the area of the sphere.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) area
  real ( kind = rk ) a
  real ( kind = rk ) as
  real ( kind = rk ) b
  real ( kind = rk ) bs
  real ( kind = rk ) c
  real ( kind = rk ) cs
  real ( kind = rk ) v1(3)
  real ( kind = rk ) v2(3)
  real ( kind = rk ) v3(3)
!
!  Compute the lengths of the sides of the spherical triangle.
!
  call sphere01_triangle_vertices_to_sides ( v1, v2, v3, as, bs, cs )
!
!  Get the spherical angles.
!
  call sphere01_triangle_sides_to_angles ( as, bs, cs, a, b, c )
!
!  Get the area.
!
  call sphere01_triangle_angles_to_area ( a, b, c, area )

  return
end
subroutine sphere01_triangle_vertices_to_centroid ( v1, v2, v3, vs )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_VERTICES_TO_CENTROID gets a spherical triangle "centroid".
!
!  Discussion:
!
!    A unit sphere in 3D satisfies the equation:
!
!      X^2 + Y^2 + Z^2 = 1
!
!    A spherical triangle is specified by three points on the sphere.
!
!    The (true) centroid of a spherical triangle is the point
!
!      VT = (XT,YT,ZT) = Integral ( X, Y, Z ) dArea / Integral 1 dArea
!
!    Note that the true centroid does NOT, in general, lie on the sphere.
!
!    The "flat" centroid VF is the centroid of the planar triangle defined by
!    the vertices of the spherical triangle.
!
!    The "spherical" centroid VS of a spherical triangle is computed by
!    the intersection of the geodesic bisectors of the triangle angles.
!    The spherical centroid lies on the sphere.
!
!    VF, VT and VS lie on a line through the center of the sphere.  We can
!    easily calculate VF by averaging the vertices, and from this determine
!    VS by normalizing.
!
!    (Of course, we still will not have actually computed VT, which lies
!    somewhere between VF and VS!)
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) V1(3), V2(3), V3(3), the vertices of the triangle.
!
!    Output, real ( kind = rk ) VS(3), the coordinates of the "spherical
!    centroid" of the spherical triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) norm
  real ( kind = rk ) v1(3)
  real ( kind = rk ) v2(3)
  real ( kind = rk ) v3(3)
  real ( kind = rk ) vs(3)

  vs(1:3) = ( v1(1:3) + v2(1:3) + v3(1:3) ) / 3.0D+00

  norm = sqrt ( sum ( vs(1:3)**2 ) )

  vs(1:3) = vs(1:3) / norm

  return
end
subroutine sphere01_triangle_vertices_to_midpoints ( v1, v2, v3, m1, m2, m3 )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_VERTICES_TO_MIDPOINTS: midsides of a spherical triangle.
!
!  Discussion:
!
!    The points are assumed to lie on the unit sphere.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) V1(3), V2(3), V3(3), the vertices of the triangle.
!
!    Output, real ( kind = rk ) M1(3), M2(3), M3(3), the coordinates of 
!    the midpoints of the sides of the spherical triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) m1(3)
  real ( kind = rk ) m2(3)
  real ( kind = rk ) m3(3)
  real ( kind = rk ) norm
  real ( kind = rk ) v1(3)
  real ( kind = rk ) v2(3)
  real ( kind = rk ) v3(3)

  m1(1:3) = ( v1(1:3) + v2(1:3) ) / 2.0D+00
  norm = sqrt ( sum ( m1(1:3)**2 ) )
  m1(1:3) = m1(1:3) / norm

  m2(1:3) = ( v2(1:3) + v3(1:3) ) / 2.0D+00
  norm = sqrt ( sum ( m2(1:3)**2 ) )
  m2(1:3) = m2(1:3) / norm

  m3(1:3) = ( v3(1:3) + v1(1:3) ) / 2.0D+00
  norm = sqrt ( sum ( m3(1:3)**2 ) )
  m3(1:3) = m3(1:3) / norm

  return
end
subroutine sphere01_triangle_vertices_to_sides ( v1, v2, v3, as, bs, cs )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_VERTICES_TO_SIDES computes spherical triangle sides.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) V1(3), V2(3), V3(3), the vertices of the spherical
!    triangle.
!
!    Output, real ( kind = rk ) AS, BS, CS, the (geodesic) length of the 
!    sides of the triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) as
  real ( kind = rk ) bs
  real ( kind = rk ) cs
  real ( kind = rk ) v1(3)
  real ( kind = rk ) v2(3)
  real ( kind = rk ) v3(3)

  as = acos ( dot_product ( v2(1:3), v3(1:3) ) )
  bs = acos ( dot_product ( v3(1:3), v1(1:3) ) )
  cs = acos ( dot_product ( v1(1:3), v2(1:3) ) )

  return
end
subroutine string_2d ( nvec, p1, p2, string_num, order, string )

!*****************************************************************************80
!
!! STRING_2D groups line segments into connected lines in 2D.
!
!  Discussion:
!
!    The routine receives an unordered set of line segments, described by
!    pairs of coordinates P1 and P2, and tries to group them
!    into ordered lists that constitute connected jagged lines.
!
!    This routine will not match two endpoints unless they are exactly equal.
!
!    On input, line segment I has endpoints P1(I), P2(I).
!
!    On output, the order of the components may have been switched.
!    That is, for some I, P1(I) and P2(I) may have been swapped.
!
!    More importantly, both points P1(I) and P2(I) may have been swapped
!    with another pair P1(J), P2(J).
!
!    The resulting coordinates will have been sorted in order
!    of the string to which they belong, and then by the order
!    of their traversal within that string.
!
!    The array STRING(I) identifies the string to which segment I belongs.
!
!    If two segments I and J have the same value of STRING, then
!    ORDER(I) and ORDER(J) give the relative order of the two segments 
!    in the string.  Thus if ORDER(I) = -3 and ORDER(J) = 2, then when 
!    the string is traversed, segment I is traversed first, then four other
!    segments are traversed, and then segment J is traversed.
!
!    For each string, the segment with ORDER(I) = 0 is the initial segment 
!    from which the entire string was "grown" (with growth possible to both the
!    left and the right).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NVEC, the number of line segments to be 
!    analyzed.
!
!    Input/output, real ( kind = rk ) P1(2,NVEC), P2VEC(2,NVEC), the 
!    line segments.
!
!    Output, integer ORDER(NVEC), the order vector.
!
!    Output, integer STRING(NVEC), the string to which each 
!    segment belongs.
!
!    Output, integer STRING_NUM, the number of strings created.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2
  integer nvec

  integer i
  integer indx
  integer isgn
  integer j
  integer jval
  integer kval
  integer match
  integer order(nvec)
  real ( kind = rk ) p1(dim_num,nvec)
  real ( kind = rk ) p2(dim_num,nvec)
  integer seed
  integer string(nvec)
  integer string_num
  real ( kind = rk ) x1val
  real ( kind = rk ) x2val
  real ( kind = rk ) y1val
  real ( kind = rk ) y2val
!
!  Mark STRING so that each segment is alone.
!
  order(1:nvec) = 0
  string(1:nvec) = nvec + i
!
!  Starting with the lowest numbered group of line segments,
!  see if any higher numbered groups belong.
!
  seed = 1
  string_num = 1
  string(seed) = string_num

  do

    x1val = p1(1,seed)
    y1val = p1(2,seed)

    x2val = p2(1,seed)
    y2val = p2(2,seed)

    jval = order(seed)
    kval = order(seed)

    do

      match = 0

      do j = 1, nvec

        if ( string_num < string(j) ) then

          if ( x1val == p1(1,j) .and. y1val == p1(2,j) ) then

            jval = jval - 1
            order(j) = jval
            string(j) = string_num
            x1val = p2(1,j)
            y1val = p2(2,j)
            match = match + 1

            call r8_swap ( p1(1,j), p2(1,j) )
            call r8_swap ( p1(2,j), p2(2,j) )

          else if ( x1val == p2(1,j) .and. y1val == p2(2,j) ) then

            jval = jval - 1
            order(j) = jval
            string(j) = string_num
            x1val = p1(1,j)
            y1val = p1(2,j)
            match = match + 1

          else if ( x2val == p1(1,j) .and. y2val == p1(2,j) ) then

            kval = kval + 1
            order(j) = kval
            string(j) = string_num
            x2val = p2(1,j)
            y2val = p2(2,j)
            match = match + 1

          else if ( x2val == p2(1,j) .and. y2val == p2(2,j) ) then

            kval = kval + 1
            order(j) = kval
            string(j) = string_num
            x2val = p1(1,j)
            y2val = p1(2,j)
            match = match + 1

            call r8_swap ( p1(1,j), p2(1,j) )
            call r8_swap ( p1(2,j), p2(2,j) )

          end if

        end if

      end do
!
!  If the string has closed on itself, then we don't want to
!  look for any more matches for this string.
!
      if ( x1val == x2val .and. y1val == y2val ) then
        exit
      end if
!
!  If we made no matches this pass, we're done.
!
      if ( match <= 0 ) then
        exit
      end if

    end do
!
!  This string is "exhausted".  Are there any line segments we
!  haven't looked at yet?
!
    seed = 0

    do i = 1, nvec
      if ( string_num < string(i) ) then
        seed = i
        string_num = string_num + 1
        string(i) = string_num
        exit
      end if
    end do

    if ( seed == 0 ) then
      exit
    end if

  end do
!
!  There are no more line segments to look at.  Renumber the
!  isolated segments.
!
!  Question: Can this ever happen?
!
  do i = 1, nvec
    if ( nvec < string(i) ) then
      string_num = string_num + 1
      string(i) = string_num
    end if
  end do
!
!  Now sort the line segments by string and by order of traversal.
!
  i = 0
  isgn = 0
  j = 0

  indx = 0

  do

    call sort_heap_external ( nvec, indx, i, j, isgn )

    if ( 0 < indx ) then

      call i4_swap ( order(i), order(j) )
      call i4_swap ( string(i), string(j) )
      call r8_swap ( p1(1,i), p1(1,j) )
      call r8_swap ( p1(2,i), p1(2,j) )
      call r8_swap ( p2(1,i), p2(1,j) )
      call r8_swap ( p2(2,i), p2(2,j) )

    else if ( indx < 0 ) then

      if ( ( string(i) < string(j) ) .or. &
           ( string(i) == string(j) .and. order(i) < order(j) ) ) then

        isgn = -1

      else

        isgn = + 1

      end if

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine super_ellipse_points_2d ( pc, r1, r2, expo, psi, n, p )

!*****************************************************************************80
!
!! SUPER_ELLIPSE_POINTS_2D returns N points on a tilted superellipse in 2D.
!
!  Discussion:
!
!    The points are "equally spaced" in the angular sense.  They are
!    not equally spaced along the perimeter.
!
!    The parametric formula of the (untilted) superellipse is:
!
!      X = R1 * cos^EXPO ( THETA )
!      Y = R2 * sin^EXPO ( THETA )
!
!    An implicit form of the (untilted) superellipse is:
!
!      (X/R1)^(2/EXPO) + (Y/R2)^(2/EXPO) = 1
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Martin Gardner,
!    The Mathematical Carnival,
!    Knopf, 1975, pages 240-254.
! 
!  Parameters:
!
!    Input, real ( kind = rk ) PC(2), the center of the superellipse.
!
!    Input, real ( kind = rk ) R1, R2, the "radius" of the superellipse
!    in the major and minor axis directions.  A circle has these values equal.
!
!    Input, real ( kind = rk ) EXPO, the exponent of the superellipse. 
!    0 = a rectangle;
!    between 0 and 1, a "rounded" rectangle;
!    1.0 = an ellipse;
!    2.0 = a diamond;
!    > 2.0 a pinched shape.
!
!    Input, real ( kind = rk ) PSI, the angle that the major axis of the
!    superellipse makes with the X axis.  A value of 0.0 means that the
!    major and minor axes of the superellipse will be the X and Y 
!    coordinate axes.
!
!    Input, integer N, the number of points desired.  N must
!    be at least 1.
!
!    Output, real ( kind = rk ) P(2,N), the coordinates of points 
!    on the superellipse.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: dim_num = 2

  real ( kind = rk ) act
  real ( kind = rk ) ast
  integer i
  real ( kind = rk ) expo
  real ( kind = rk ) p(dim_num,n)
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) psi
  real ( kind = rk ) r1
  real ( kind = rk ) r2
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) sct
  real ( kind = rk ) sst
  real ( kind = rk ) theta

  do i = 1, n

    theta = ( 2.0D+00 * r8_pi * real ( i - 1, kind = rk ) ) &
      / real ( n, kind = rk )

    act = abs ( cos ( theta ) )
    sct = sign ( 1.0D+00, cos ( theta ) )
    ast = abs ( sin ( theta ) )
    sst = sign ( 1.0D+00, sin ( theta ) )

    p(1,i) = pc(1) + r1 * cos ( psi ) * sct * ( act ) ** expo &
                   - r2 * sin ( psi ) * sst * ( ast ) ** expo
 
    p(2,i) = pc(2) + r1 * sin ( psi ) * sct * ( act ) ** expo &
                   + r2 * cos ( psi ) * sst * ( ast ) ** expo

  end do

  return
end
subroutine tetrahedron_barycentric_3d ( tetra, p, c )

!*****************************************************************************80
!
!! tetrahedron_barycentric_3d(): barycentric coordinates of a point in 3D.
!
!  Discussion:
!
!    The barycentric coordinates of a point P with respect to
!    a tetrahedron are a set of four values C(1:4), each associated
!    with a vertex of the tetrahedron.  The values must sum to 1.
!    If all the values are between 0 and 1, the point is contained
!    within the tetrahedron.
!
!    The barycentric coordinate of point P related to vertex A can be
!    interpreted as the ratio of the volume of the tetrahedron with 
!    vertex A replaced by vertex P to the volume of the original 
!    tetrahedron.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) TETRA(3,4) the tetrahedron vertices.
!
!    Input, real ( kind = rk ) P(3), the point to be checked.
!
!    Output, real ( kind = rk ) C(4), the barycentric coordinates of P with
!    respect to the tetrahedron.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer, parameter :: rhs_num = 1

  real ( kind = rk ) a(dim_num,dim_num+rhs_num)
  real ( kind = rk ) c(dim_num+1)
  integer i
  integer info
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) tetra(dim_num,4)
!
!  Set up the linear system
!
!    ( X2-X1  X3-X1  X4-X1 ) C2    X - X1
!    ( Y2-Y1  Y3-Y1  Y4-Y1 ) C3  = Y - Y1
!    ( Z2-Z1  Z3-Z1  Z4-Z1 ) C4    Z - Z1
!
!  which is satisfied by the barycentric coordinates of P.
!
  a(1:dim_num,1:3) = tetra(1:dim_num,2:4)
  a(1:dim_num,4) = p(1:dim_num)

  do i = 1, dim_num
    a(i,1:4) = a(i,1:4) - tetra(i,1)
  end do
!
!  Solve the linear system.
!
  call r8mat_solve ( dim_num, rhs_num, a, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TETRAHEDRON_BARYCENTRIC_3D - Fatal error!'
    write ( *, '(a)' ) '  The linear system is singular.'
    write ( *, '(a)' ) '  The input data does not form a proper tetrahedron.'
    stop 1
  end if

  c(2:4) = a(1:dim_num,4)

  c(1) = 1.0D+00 - sum ( c(2:4) )

  return
end
subroutine tetrahedron_centroid_3d ( tetra, centroid )

!*****************************************************************************80
!
!! TETRAHEDRON_CENTROID_3D computes the centroid of a tetrahedron in 3D.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) TETRA(3,4) the tetrahedron vertices.
!
!    Output, real ( kind = rk ) CENTROID(3), the coordinates of the centroid.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) centroid(dim_num)
  integer i
  real ( kind = rk ) tetra(dim_num,4)

  do i = 1, dim_num
    centroid(i) = sum ( tetra(i,1:4) ) / 4.0D+00
  end do

  return
end
subroutine tetrahedron_circumsphere_3d ( tetra, r, pc )

!*****************************************************************************80
!
!! tetrahedron_circumsphere_3d() computes the circumsphere of a tetrahedron in 3D.
!
!  Discussion:
!
!    The circumsphere, or circumscribed sphere, of a tetrahedron is the 
!    sphere that passes through the four vertices.  The circumsphere is
!    not necessarily the smallest sphere that contains the tetrahedron.
!
!    Surprisingly, the diameter of the sphere can be found by solving
!    a 3 by 3 linear system.  This is because the vectors P2 - P1,
!    P3 - P1 and P4 - P1 are secants of the sphere, and each forms a
!    right triangle with the diameter through P1.  Hence, the dot product of
!    P2 - P1 with that diameter is equal to the square of the length
!    of P2 - P1, and similarly for P3 - P1 and P4 - P1.  This determines
!    the diameter vector originating at P1, and hence the radius and
!    center.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 August 2005
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
!  Parameters:
!
!    Input, real ( kind = rk ) TETRA(3,4) the tetrahedron vertices.
!
!    Output, real ( kind = rk ) R, PC(3), the center of the
!    circumscribed sphere, and its radius.  If the linear system is
!    singular, then R = -1, PC(1:3) = 0.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer, parameter :: rhs_num = 1

  real ( kind = rk ) a(dim_num,dim_num+rhs_num)
  integer i
  integer info
  integer j
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) r
  real ( kind = rk ) tetra(dim_num,4)
!
!  Set up the linear system.
!
  a(1:dim_num,1:3) = transpose ( tetra(1:dim_num,2:4) )

  do j = 1, dim_num
    a(1:dim_num,j) = a(1:dim_num,j) - tetra(j,1)
  end do

  do i = 1, 3
    a(i,4) = sum ( a(i,1:3)**2 )
  end do
!
!  Solve the linear system.
!
  call r8mat_solve ( dim_num, rhs_num, a, info )
!
!  If the system was singular, return a consolation prize.
!
  if ( info /= 0 ) then
    r = -1.0D+00
    pc(1:dim_num) = 0.0D+00
    return
  end if
!
!  Compute the radius and center.
!
  r = 0.5D+00 * sqrt ( sum ( a(1:dim_num,4)**2 ) )

  pc(1:dim_num) = tetra(1:dim_num,1) + 0.5D+00 * a(1:dim_num,4)

  return
end
subroutine tetrahedron_volume_3d ( tetra, volume )

!*****************************************************************************80
!
!! TETRAHEDRON_VOLUME_3D computes the volume of a tetrahedron in 3D.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) TETRA(3,4), the vertices of the tetrahedron.
!
!    Output, real ( kind = rk ) VOLUME, the volume of the tetrahedron.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) a(4,4)
  real ( kind = rk ) r8mat_det_4d
  real ( kind = rk ) tetra(dim_num,4)
  real ( kind = rk ) volume

  a(1:dim_num,1:4) = tetra(1:dim_num,1:4)
  a(4,1:4) = 1.0D+00

  volume = abs ( r8mat_det_4d ( a ) ) / 6.0D+00

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
!    06 August 2005
!
!  Author:
!
!    John Burkardt
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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine tmat_init ( a )

!*****************************************************************************80
!
!! TMAT_INIT initializes the geometric transformation matrix.
!
!  Discussion:
!
!    The geometric transformation matrix can be thought of as a 4 by 4
!    matrix "A" having components:
!
!      r11 r12 r13 t1
!      r21 r22 r23 t2
!      r31 r32 r33 t3
!        0   0   0  1
!
!    This matrix encodes the rotations, scalings and translations that
!    are applied to graphical objects.
!
!    A point P = (x,y,z) is rewritten in "homogeneous coordinates" as
!    PH = (x,y,z,1).  Then to apply the transformations encoded in A to
!    the point P, we simply compute A * PH.
!
!    Individual transformations, such as a scaling, can be represented
!    by simple versions of the transformation matrix.  If the matrix
!    A represents the current set of transformations, and we wish to
!    apply a new transformation B, then the original points are
!    transformed twice:  B * ( A * PH ).  The new transformation B can
!    be combined with the original one A, to give a single matrix C that
!    encodes both transformations: C = B * A.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Foley, Andries van Dam, Steven Feiner, John Hughes,
!    Computer Graphics, Principles and Practice,
!    Second Edition,
!    Addison Wesley, 1990.
!
!  Parameters:
!
!    Input, real ( kind = rk ) A(4,4), the geometric transformation matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(4,4)
  integer i
  integer j

  do i = 1, 4
    do j = 1, 4
      if ( i == j ) then
        a(i,j) = 1.0D+00
      else
        a(i,j) = 0.0D+00
      end if
    end do
  end do

  return
end
subroutine tmat_mxm ( a, b, c )

!*****************************************************************************80
!
!! TMAT_MXM multiplies two geometric transformation matrices.
!
!  Discussion:
!
!    The product is accumulated in a temporary array, and then assigned
!    to the result.  Therefore, it is legal for any two, or all three,
!    of the arguments to share memory.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Foley, Andries van Dam, Steven Feiner, John Hughes,
!    Computer Graphics, Principles and Practice,
!    Second Edition,
!    Addison Wesley, 1990.
!
!  Parameters:
!
!    Input, real ( kind = rk ) A(4,4), the first geometric transformation matrix.
!
!    Input, real ( kind = rk ) B(4,4), the second geometric transformation
!    matrix.
!
!    Output, real ( kind = rk ) C(4,4), the product A * B.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(4,4)
  real ( kind = rk ) b(4,4)
  real ( kind = rk ) c(4,4)

  c(1:4,1:4) = matmul ( a(1:4,1:4), b(1:4,1:4) )

  return
end
subroutine tmat_mxp ( a, x, y )

!*****************************************************************************80
!
!! TMAT_MXP multiplies a geometric transformation matrix times a point.
!
!  Discussion:
!
!    The matrix will normally have the form
!
!      xx xy xz tx
!      yx yy yz ty
!      zx zy zz tz
!       0  0  0  1
!
!    where the 3x3 initial block controls rotations and scalings,
!    and the values [ tx, ty, tz ] implement a translation.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Foley, Andries van Dam, Steven Feiner, John Hughes,
!    Computer Graphics, Principles and Practice,
!    Second Edition,
!    Addison Wesley, 1990.
!
!  Parameters:
!
!    Input, real ( kind = rk ) A(4,4), the geometric transformation matrix.
!
!    Input, real ( kind = rk ) X(3), the point to be multiplied.  The fourth
!    component of X is implicitly assigned the value of 1.
!
!    Output, real ( kind = rk ) Y(3), the result of A*X.  The product is
!    accumulated in a temporary vector, and then assigned to the result.
!    Therefore, it is legal for X and Y to share memory.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(4,4)
  real ( kind = rk ) x(3)
  real ( kind = rk ) y(3)

  y(1:3) = a(1:3,4) + matmul ( a(1:3,1:3), x(1:3) )

  return
end
subroutine tmat_mxp2 ( a, n, x, y )

!*****************************************************************************80
!
!! TMAT_MXP2 multiplies a geometric transformation matrix times N points.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Foley, Andries van Dam, Steven Feiner, John Hughes,
!    Computer Graphics, Principles and Practice,
!    Second Edition,
!    Addison Wesley, 1990.
!
!  Parameters:
!
!    Input, real ( kind = rk ) A(4,4), the geometric transformation matrix.
!
!    Input, integer N, the number of points to be multiplied.
!
!    Input, real ( kind = rk ) X(3,N), the points to be multiplied.
!
!    Output, real ( kind = rk ) Y(3,N), the transformed points.  Each product is
!    accumulated in a temporary vector, and then assigned to the
!    result.  Therefore, it is legal for X and Y to share memory.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(4,4)
  integer i
  real ( kind = rk ) x(3,n)
  real ( kind = rk ) y(3,n)

  do i = 1, 3
    y(i,1:n) = a(i,4)
  end do

  y(1:3,1:n) = y(1:3,1:n) + matmul ( a(1:3,1:3), x(1:3,1:n) )

  return
end
subroutine tmat_mxv ( a, x, y )

!*****************************************************************************80
!
!! TMAT_MXV multiplies a geometric transformation matrix times a vector.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Foley, Andries van Dam, Steven Feiner, John Hughes,
!    Computer Graphics, Principles and Practice,
!    Second Edition,
!    Addison Wesley, 1990.
!
!  Parameters:
!
!    Input, real ( kind = rk ) A(4,4), the geometric transformation matrix.
!
!    Input, real ( kind = rk ) X(3), the vector to be multiplied.  The fourth
!    component of X is implicitly assigned the value of 1.
!
!    Output, real ( kind = rk ) Y(3), the result of A*X.  The product is
!    accumulated in a temporary vector, and then assigned to the result. 
!    Therefore, it is legal for X and Y to share memory.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(4,4)
  real ( kind = rk ) x(3)
  real ( kind = rk ) y(3)

  y(1:3) = a(1:3,4) + matmul ( a(1:3,1:3), x(1:3) )

  return
end
subroutine tmat_rot_axis ( a, angle, axis, b )

!*****************************************************************************80
!
!! TMAT_ROT_AXIS: coordinate axis rotation to geometric transformation matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Foley, Andries van Dam, Steven Feiner, John Hughes,
!    Computer Graphics, Principles and Practice,
!    Second Edition,
!    Addison Wesley, 1990.
!
!  Parameters:
!
!    Input, real ( kind = rk ) A(4,4), the current geometric transformation
!    matrix.
!
!    Input, real ( kind = rk ) ANGLE, the angle, in degrees, of the rotation.
!
!    Input, character AXIS, is 'X', 'Y' or 'Z', specifying the coordinate
!    axis about which the rotation occurs.
!
!    Output, real ( kind = rk ) B(4,4), the modified geometric 
!    transformation matrix.
!    A and B may share the same memory.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real      ( kind = rk ) a(4,4)
  real      ( kind = rk ) angle
  real      ( kind = rk ) angle_rad
  character              axis
  real      ( kind = rk ) b(4,4)
  real      ( kind = rk ) c(4,4)
  real      ( kind = rk ) degrees_to_radians

  angle_rad = degrees_to_radians ( angle )

  call tmat_init ( c )

  if ( axis == 'X' .or. axis == 'x' ) then
    c(2,2) =   cos ( angle_rad )
    c(2,3) = - sin ( angle_rad )
    c(3,2) =   sin ( angle_rad )
    c(3,3) =   cos ( angle_rad )
  else if ( axis == 'Y' .or. axis == 'y' ) then
    c(1,1) =   cos ( angle_rad )
    c(1,3) =   sin ( angle_rad )
    c(3,1) = - sin ( angle_rad )
    c(3,3) =   cos ( angle_rad )
  else if ( axis == 'Z' .or. axis == 'z' ) then
    c(1,1) =   cos ( angle_rad )
    c(1,2) = - sin ( angle_rad )
    c(2,1) =   sin ( angle_rad )
    c(2,2) =   cos ( angle_rad )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TMAT_ROT_AXIS - Fatal error!'
    write ( *, '(a)' ) '  Illegal rotation axis: ' // axis
    write ( *, '(a)' ) '  Legal choices are ''X'', ''Y'', or ''Z''.'
    stop 1
  end if

  b(1:4,1:4) = matmul ( c(1:4,1:4), a(1:4,1:4) )

  return
end
subroutine tmat_rot_vector ( a, angle, axis, b )

!*****************************************************************************80
!
!! TMAT_ROT_VECTOR: arbitrary axis rotation to geometric transformation matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Foley, Andries van Dam, Steven Feiner, John Hughes,
!    Computer Graphics, Principles and Practice,
!    Second Edition,
!    Addison Wesley, 1990.
!
!  Parameters:
!
!    Input, real ( kind = rk ) A(4,4), the current geometric transformation
!    matrix.
!
!    Input, real ( kind = rk ) ANGLE, the angle, in degrees, of the rotation.
!
!    Input, real ( kind = rk ) AXIS(3), the axis vector about which 
!    rotation occurs.  AXIS may not be the zero vector.
!
!    Output, real ( kind = rk ) B(4,4), the modified geometric 
!    transformation matrix.
!    A and B may share the same memory.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(4,4)
  real ( kind = rk ) angle
  real ( kind = rk ) angle_rad
  real ( kind = rk ) axis(3)
  real ( kind = rk ) b(4,4)
  real ( kind = rk ) c(4,4)
  real ( kind = rk ) ca
  real ( kind = rk ) degrees_to_radians
  real ( kind = rk ) norm
  real ( kind = rk ) sa
  real ( kind = rk ) v1
  real ( kind = rk ) v2
  real ( kind = rk ) v3

  v1 = axis(1)
  v2 = axis(2)
  v3 = axis(3)

  norm = sqrt ( v1 * v1 + v2 * v2 + v3 * v3 )

  if ( norm == 0.0D+00 ) then
    return
  end if

  v1 = v1 / norm
  v2 = v2 / norm
  v3 = v3 / norm

  angle_rad = degrees_to_radians ( angle )
  ca = cos ( angle_rad )
  sa = sin ( angle_rad )

  call tmat_init ( c )

  c(1,1) =                    v1 * v1 + ca * ( 1.0D+00 - v1 * v1 )
  c(1,2) = ( 1.0D+00 - ca ) * v1 * v2 - sa * v3
  c(1,3) = ( 1.0D+00 - ca ) * v1 * v3 + sa * v2

  c(2,1) = ( 1.0D+00 - ca ) * v2 * v1 + sa * v3
  c(2,2) =                    v2 * v2 + ca * ( 1.0D+00 - v2 * v2 )
  c(2,3) = ( 1.0D+00 - ca ) * v2 * v3 - sa * v1

  c(3,1) = ( 1.0D+00 - ca ) * v3 * v1 - sa * v2
  c(3,2) = ( 1.0D+00 - ca ) * v3 * v2 + sa * v1
  c(3,3) =                    v3 * v3 + ca * ( 1.0D+00 - v3 * v3 )

  b(1:4,1:4) = matmul ( c(1:4,1:4), a(1:4,1:4) )

  return
end
subroutine tmat_scale ( a, s, b )

!*****************************************************************************80
!
!! TMAT_SCALE applies a scaling to the geometric transformation matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Foley, Andries van Dam, Steven Feiner, John Hughes,
!    Computer Graphics, Principles and Practice,
!    Second Edition,
!    Addison Wesley, 1990.
!
!  Parameters:
!
!    Input, real ( kind = rk ) A(4,4), the current geometric transformation
!    matrix.
!
!    Input, real ( kind = rk ) S(3), the scalings to be applied to the 
!    X, Y and Z coordinates.
!
!    Output, real ( kind = rk ) B(4,4), the modified geometric transformation
!    matrix.  A and B may share the same memory.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(4,4)
  real ( kind = rk ) b(4,4)
  real ( kind = rk ) c(4,4)
  real ( kind = rk ) s(3)

  call tmat_init ( c )

  c(1,1) = s(1)
  c(2,2) = s(2)
  c(3,3) = s(3)

  b(1:4,1:4) = matmul ( c(1:4,1:4), a(1:4,1:4) )

  return
end
subroutine tmat_shear ( a, axis, s, b )

!*****************************************************************************80
!
!! TMAT_SHEAR applies a shear to the geometric transformation matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Foley, Andries van Dam, Steven Feiner, John Hughes,
!    Computer Graphics, Principles and Practice,
!    Second Edition,
!    Addison Wesley, 1990.
!
!  Parameters:
!
!    Input, real ( kind = rk ) A(4,4), the current geometric transformation
!    matrix.
!
!    Input, character ( len = 2 ) AXIS, is 'XY', 'XZ', 'YX', 'YZ', 'ZX' or 'ZY',
!    specifying the shear equation:
!
!      XY:  x' = x + s * y;
!      XZ:  x' = x + s * z;
!      YX:  y' = y + s * x;
!      YZ:  y' = y + s * z;
!      ZX:  z' = z + s * x;
!      ZY:  z' = z + s * y.
!
!    Input, real ( kind = rk ) S, the shear coefficient.
!
!    Output, real ( kind = rk ) B(4,4), the modified geometric transformation
!    matrix.  A and B may share the same memory.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(4,4)
  character ( len = 2 ) axis
  real ( kind = rk ) b(4,4)
  real ( kind = rk ) c(4,4)
  real ( kind = rk ) s

  call tmat_init ( c )

  if ( axis == 'XY' .or. axis == 'xy' ) then
    c(1,2) = s
  else if ( axis == 'XZ' .or. axis == 'xz' ) then
    c(1,3) = s
  else if ( axis == 'YX' .or. axis == 'yx' ) then
    c(2,1) = s
  else if ( axis == 'YZ' .or. axis == 'yz' ) then
    c(2,3) = s
  else if ( axis == 'ZX' .or. axis == 'zx' ) then
    c(3,1) = s
  else if ( axis == 'ZY' .or. axis == 'zy' ) then
    c(3,2) = s
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TMAT_SHEAR - Fatal error!'
    write ( *, '(a)' ) '  Illegal shear axis: "' // axis // '".'
    write ( *, '(a)' ) '  Legal choices are XY, XZ, YX, YZ, ZX, or ZY.'
    stop 1
  end if

  b(1:4,1:4) = matmul ( c(1:4,1:4), a(1:4,1:4) )

  return
end
subroutine tmat_trans ( a, t, b )

!*****************************************************************************80
!
!! TMAT_TRANS applies a translation to the geometric transformation matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Foley, Andries van Dam, Steven Feiner, John Hughes,
!    Computer Graphics, Principles and Practice,
!    Second Edition,
!    Addison Wesley, 1990.
!
!  Parameters:
!
!    Input, real ( kind = rk ) A(4,4), the current geometric transformation
!    matrix.
!
!    Input, real ( kind = rk ) T(3), the translation.  This may be thought
!    of as the point that the origin moves to under the translation.
!
!    Output, real ( kind = rk ) B(4,4), the modified transformation matrix.
!    A and B may share the same memory.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(4,4)
  real ( kind = rk ) b(4,4)
  real ( kind = rk ) t(3)

  b(1:4,1:4) = a(1:4,1:4)

  b(1:3,4) = b(1:3,4) + t(1:3)

  return
end
function torus_area_3d ( r1, r2 )

!*****************************************************************************80
!
!! TORUS_AREA_3D returns the area of a torus in 3D.
!
!  Discussion:
!
!    A torus with radii R1 and R2 is the set of points P satisfying:
!
!    ( sqrt ( P(1)^2 + P(2)^2 ) - R1 )^2 + P(3)^2 <= R2^2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) R1, R2, the two radii that define the torus.
!
!    Output, real ( kind = rk ) TORUS_AREA_3D, the area of the torus.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r1
  real ( kind = rk ) r2
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) torus_area_3d

  torus_area_3d = 4.0D+00 * r8_pi * r8_pi * r1 * r2

  return
end
subroutine torus_volume_3d ( r1, r2, volume )

!*****************************************************************************80
!
!! TORUS_VOLUME_3D computes the volume of a torus in 3D.
!
!  Discussion:
!
!    A torus with radii R1 and R2 is the set of points P satisfying:
!
!    ( sqrt ( P(1)*^ + P(2)^2 ) - R1 )^2 + P(3)^2 <= R2^2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) R1, R2, the "inner" and "outer" radii of the
!    torus.
!
!    Output, real ( kind = rk ) VOLUME, the volume of the torus.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r1
  real ( kind = rk ) r2
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) volume

  volume = 2.0D+00 * r8_pi * r8_pi * r1 * r2 * r2

  return
end
subroutine tp_to_xyz ( theta, phi, v )

!*****************************************************************************80
!
!! TP_TO_XYZ converts unit spherical TP coordinates to XYZ coordinates.
!
!  Discussion:
!
!    The point is assume to lie on the unit sphere centered at the origin.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) THETA, PHI, the angular coordinates of a point
!    on the unit sphere.
!
!    Output, real ( kind = rk ) V(3), the XYZ coordinates.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) phi
  real ( kind = rk ) theta
  real ( kind = rk ) v(3)

  v(1) = cos ( theta ) * sin ( phi )
  v(2) = sin ( theta ) * sin ( phi )
  v(3) =                 cos ( phi )

  return
end
subroutine triangle_angles_2d ( t, angle )

!*****************************************************************************80
!
!! TRIANGLE_ANGLES_2D computes the angles of a triangle in 2D.
!
!  Discussion:
!
!    The law of cosines is used:
!
!      C^2 = A^2 + B^2 - 2 * A * B * COS ( GAMMA )
!
!    where GAMMA is the angle opposite side C.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) T(2,3), the triangle vertices.
!
!    Output, real ( kind = rk ) ANGLE(3), the angles opposite
!    sides P1-P2, P2-P3 and P3-P1, in radians.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) a
  real ( kind = rk ) angle(3)
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) r8_acos
  real ( kind = rk ) t(dim_num,3)
!
!  Compute the length of each side.
!
  a = sqrt ( sum ( ( t(1:dim_num,1) - t(1:dim_num,2) )**2 ) )
  b = sqrt ( sum ( ( t(1:dim_num,2) - t(1:dim_num,3) )**2 ) )
  c = sqrt ( sum ( ( t(1:dim_num,3) - t(1:dim_num,1) )**2 ) )
!
!  Take care of ridiculous special cases.
!
  if ( a == 0.0D+00 .and. b == 0.0D+00 .and. c == 0.0D+00 ) then
    angle(1:3) = 2.0D+00 * r8_pi / 3.0D+00
    return
  end if

  if ( c == 0.0D+00 .or. a == 0.0D+00 ) then
    angle(1) = r8_pi
  else
    angle(1) = r8_acos ( ( c * c + a * a - b * b ) / ( 2.0D+00 * c * a ) )
  end if

  if ( a == 0.0D+00 .or. b == 0.0D+00 ) then
    angle(2) = r8_pi
  else
    angle(2) = r8_acos ( ( a * a + b * b - c * c ) / ( 2.0D+00 * a * b ) )
  end if

  if ( b == 0.0D+00 .or. c == 0.0D+00 ) then
    angle(3) = r8_pi
  else
    angle(3) = r8_acos ( ( b * b + c * c - a * a ) / ( 2.0D+00 * b * c ) )
  end if

  return
end
subroutine triangle_angles_3d ( t, angle )

!*****************************************************************************80
!
!! TRIANGLE_ANGLES_3D computes the angles of a triangle in 3D.
!
!  Discussion:
!
!    The law of cosines is used:
!
!      C * C = A * A + B * B - 2 * A * B * COS ( GAMMA )
!
!    where GAMMA is the angle opposite side C.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) T(3,3), the triangle vertices.
!
!    Output, real ( kind = rk ) ANGLE(3), the angles opposite
!    sides P1-P2, P2-P3 and P3-P1, in radians.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) a
  real ( kind = rk ) angle(3)
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) r8_acos
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) t(dim_num,3)
!
!  Compute the length of each side.
!
  a = sqrt ( sum ( ( t(1:dim_num,1) - t(1:dim_num,2) )**2 ) )
  b = sqrt ( sum ( ( t(1:dim_num,2) - t(1:dim_num,3) )**2 ) )
  c = sqrt ( sum ( ( t(1:dim_num,3) - t(1:dim_num,1) )**2 ) )
!
!  Take care of a ridiculous special case.
!
  if ( a == 0.0D+00 .and. b == 0.0D+00 .and. c == 0.0D+00 ) then
    angle(1:3) = 2.0D+00 * r8_pi / 3.0D+00
    return
  end if

  if ( c == 0.0D+00 .or. a == 0.0D+00 ) then
    angle(1) = r8_pi
  else
    angle(1) = r8_acos ( ( c * c + a * a - b * b ) / ( 2.0D+00 * c * a ) )
  end if

  if ( a == 0.0D+00 .or. b == 0.0D+00 ) then
    angle(2) = r8_pi
  else
    angle(2) = r8_acos ( ( a * a + b * b - c * c ) / ( 2.0D+00 * a * b ) )
  end if

  if ( b == 0.0D+00 .or. c == 0.0D+00 ) then
    angle(3) = r8_pi
  else
    angle(3) = r8_acos ( ( b * b + c * c - a * a ) / ( 2.0D+00 * b * c ) )
  end if

  return
end
subroutine triangle_area_2d ( t, area )

!*****************************************************************************80
!
!! TRIANGLE_AREA_2D computes the area of a triangle in 2D.
!
!  Discussion:
!
!    If the triangle's vertices are given in counter clockwise order,
!    the area will be positive.  If the triangle's vertices are given
!    in clockwise order, the area will be negative!
!
!    An earlier version of this routine always returned the absolute
!    value of the computed area.  I am convinced now that that is
!    a less useful result!  For instance, by returning the signed 
!    area of a triangle, it is possible to easily compute the area 
!    of a nonconvex polygon as the sum of the (possibly negative) 
!    areas of triangles formed by node 1 and successive pairs of vertices.
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
!  Parameters:
!
!    Input, real ( kind = rk ) T(2,3), the triangle vertices.
!
!    Output, real ( kind = rk ) AREA, the area of the triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) area
  real ( kind = rk ) t(dim_num,3)

  area = 0.5D+00 * ( &
      t(1,1) * ( t(2,2) - t(2,3) ) &
    + t(1,2) * ( t(2,3) - t(2,1) ) &
    + t(1,3) * ( t(2,1) - t(2,2) ) )

  return
end
subroutine triangle_area_3d ( t, area )

!*****************************************************************************80
!
!! TRIANGLE_AREA_3D computes the area of a triangle in 3D.
!
!  Discussion:
!
!    This routine uses the fact that the norm of the cross product 
!    of two vectors is the area of the parallelogram they form.  
!
!    Therefore, the area of the triangle is half of that value.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 December 2004
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
!  Parameters:
!
!    Input, real ( kind = rk ) T(3,3), the triangle vertices.
!
!    Output, real ( kind = rk ) AREA, the area of the triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) area
  real ( kind = rk ) cross(dim_num)
  real ( kind = rk ) t(dim_num,3)
!
!  Compute the cross product vector.
!
  cross(1) = ( t(2,2) - t(2,1) ) * ( t(3,3) - t(3,1) ) &
           - ( t(3,2) - t(3,1) ) * ( t(2,3) - t(2,1) )

  cross(2) = ( t(3,2) - t(3,1) ) * ( t(1,3) - t(1,1) ) &
           - ( t(1,2) - t(1,1) ) * ( t(3,3) - t(3,1) )

  cross(3) = ( t(1,2) - t(1,1) ) * ( t(2,3) - t(2,1) ) &
           - ( t(2,2) - t(2,1) ) * ( t(1,3) - t(1,1) )

  area = 0.5D+00 * sqrt ( sum ( cross(1:3)**2 ) )

  return
end
subroutine triangle_area_3d_2 ( t, area )

!*****************************************************************************80
!
!! TRIANGLE_AREA_3D_2 computes the area of a triangle in 3D.
!
!  Discussion:
!
!    This routine computes the area "the hard way".
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) T(3,3), the triangle vertices.
!
!    Output, real ( kind = rk ) AREA, the area of the triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) alpha
  real ( kind = rk ) area
  real ( kind = rk ) base
  real ( kind = rk ) dot
  real ( kind = rk ) height
  real ( kind = rk ) t(dim_num,3)
!
!  Find the projection of (P3-P1) onto (P2-P1).
!
  dot = ( t(1,2) - t(1,1) ) * ( t(1,3) - t(1,1) ) &
      + ( t(2,2) - t(2,1) ) * ( t(2,3) - t(2,1) ) &
      + ( t(3,2) - t(3,1) ) * ( t(3,3) - t(3,1) )
!
!  Find the length of (P2-P1).
!
  base = sqrt ( ( t(1,2) - t(1,1) )**2 &
              + ( t(2,2) - t(2,1) )**2 &
              + ( t(3,2) - t(3,1) )**2 )
!
!  The height of the triangle is the length of (P3-P1) after its
!  projection onto (P2-P1) has been subtracted.
!
  if ( base == 0.0D+00 ) then

    height = 0.0D+00

  else

    alpha = dot / ( base * base )

    height = sqrt ( &
        ( t(1,1) + alpha * ( t(1,2) - t(1,1) ) - t(1,3) )**2 &
      + ( t(2,1) + alpha * ( t(2,2) - t(2,1) ) - t(2,3) )**2 &
      + ( t(3,1) + alpha * ( t(3,2) - t(3,1) ) - t(3,3) )**2 )

  end if

  area = 0.5D+00 * base * height

  return
end
subroutine triangle_area_3d_3 ( t, area )

!*****************************************************************************80
!
!! TRIANGLE_AREA_3D_3 computes the area of a triangle in 3D.
!
!  Discussion:
!
!    This routine uses Heron's formula
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 December 2004
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
!  Parameters:
!
!    Input, real ( kind = rk ) T(3,3), the triangle vertices.
!
!    Output, real ( kind = rk ) AREA, the area of the triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) area
  integer i
  integer j
  integer jp1
  real ( kind = rk ) s(3)
  real ( kind = rk ) t(dim_num,3)

  do j = 1, 3
    jp1 = mod ( j, 3 ) + 1
    s(j) = 0.0D+00
    do i = 1, dim_num
      s(j) = s(j) + ( t(i,j) - t(i,jp1) )**2
    end do
    s(j) = sqrt ( s(j) )
  end do

  area = (   s(1) + s(2) + s(3) ) &
       * ( - s(1) + s(2) + s(3) ) &
       * (   s(1) - s(2) + s(3) ) &
       * (   s(1) + s(2) - s(3) )

  if ( area < 0.0D+00 ) then
    area = -1.0D+00
    return
  end if

  area = 0.25D+00 * sqrt ( area )

  return
end
subroutine triangle_area_heron ( s, area )

!*****************************************************************************80
!
!! TRIANGLE_AREA_HERON computes the area of a triangle using Heron's formula.
!
!  Discussion:
!
!    The formula is valid for any spatial dimension, depending only
!    on the lengths of the sides, and not the coordinates of the vertices.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) S(3), the lengths of the three sides.
!
!    Output, real ( kind = rk ) AREA, the area of the triangle, or -1.0 if the
!    sides cannot constitute a triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) area
  real ( kind = rk ) s(3)

  area = (   s(1) + s(2) + s(3) ) &
       * ( - s(1) + s(2) + s(3) ) &
       * (   s(1) - s(2) + s(3) ) &
       * (   s(1) + s(2) - s(3) )

  if ( area < 0.0D+00 ) then
    area = -1.0D+00
    return
  end if

  area = 0.25D+00 * sqrt ( area )

  return
end
subroutine triangle_area_vector_3d ( t, area_vector )

!*****************************************************************************80
!
!! TRIANGLE_AREA_VECTOR_3D computes the area vector of a triangle in 3D.
!
!  Discussion:
!
!    The "area vector" of a triangle is simply a cross product of,
!    for instance, the vectors (V2-V1) and (V3-V1), where V1, V2
!    and V3 are the vertices of the triangle.
!    
!    The norm of the cross product vector of two vectors is the area 
!    of the parallelogram they form.  
!
!    Therefore, the area of the triangle is half of the norm of the
!    area vector:
!
!      area = 0.5 * sqrt ( sum ( area_vector(1:3)^2 ) )
!
!    The reason for looking at the area vector rather than the area
!    is that this makes it possible to compute the area of a flat
!    polygon in 3D by summing the areas of the triangles that form
!    a decomposition of the polygon, while allowing for both positive
!    and negative areas.  (Sum the vectors, THEN take the norm and
!    multiply by 1/2).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 October 2005
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
!  Parameters:
!
!    Input, real ( kind = rk ) T(3,3), the triangle vertices.
!
!    Output, real ( kind = rk ) AREA_VECTOR(3), the area vector of the triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) area_vector(dim_num)
  real ( kind = rk ) t(dim_num,3)

  area_vector(1) = ( t(2,2) - t(2,1) ) * ( t(3,3) - t(3,1) ) &
                 - ( t(3,2) - t(3,1) ) * ( t(2,3) - t(2,1) )

  area_vector(2) = ( t(3,2) - t(3,1) ) * ( t(1,3) - t(1,1) ) &
                 - ( t(1,2) - t(1,1) ) * ( t(3,3) - t(3,1) )

  area_vector(3) = ( t(1,2) - t(1,1) ) * ( t(2,3) - t(2,1) ) &
                 - ( t(2,2) - t(2,1) ) * ( t(1,3) - t(1,1) )

  return
end
subroutine triangle_barycentric_2d ( t, p, xsi )

!*****************************************************************************80
!
!! TRIANGLE_BARYCENTRIC_2D finds the barycentric coordinates of a point in 2D.
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
!  Parameters:
!
!    Input, real ( kind = rk ) T(2,3), the triangle vertices.
!    The vertices should be given in counter clockwise order.
!
!    Input, real ( kind = rk ) P(2), the point to be checked.
!
!    Output, real ( kind = rk ) XSI(3), the barycentric coordinates of P
!    with respect to the triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2
  integer, parameter :: rhs_num = 1

  real ( kind = rk ) a(dim_num,dim_num+rhs_num)
  integer info
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) t(dim_num,3)
  real ( kind = rk ) xsi(dim_num+1)
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
  call r8mat_solve ( dim_num, rhs_num, a, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGLE_BARYCENTRIC_2D - Fatal error!'
    write ( *, '(a)' ) '  The linear system is singular.'
    write ( *, '(a)' ) '  The input data does not form a proper triangle.'
    stop 1
  end if

  xsi(1) = a(1,3)
  xsi(2) = a(2,3)
  xsi(3) = 1.0D+00 - xsi(1) - xsi(2)

  return
end
subroutine triangle_centroid_2d ( t, centroid )

!*****************************************************************************80
!
!! TRIANGLE_CENTROID_2D computes the centroid of a triangle in 2D.
!
!  Discussion:
!
!    The centroid of a triangle can also be considered the 
!    center of gravity, or center of mass, assuming that the triangle 
!    is made of a thin uniform sheet of massy material.
!
!    The centroid of a triangle is the intersection of the medians.
!
!    A median of a triangle is a line connecting a vertex to the
!    midpoint of the opposite side.
!
!    In barycentric coordinates, in which the vertices of the triangle
!    have the coordinates (1,0,0), (0,1,0) and (0,0,1), the centroid
!    has coordinates (1/3,1/3,1/3).
!
!    In geometry, the centroid of a triangle is often symbolized by "G".
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
!  Reference:
!
!    Adrian Bowyer, John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983,
!    ISBN: 0408012420.
!
!  Parameters:
!
!    Input, real ( kind = rk ) T(2,3), the triangle vertices.
!
!    Output, real ( kind = rk ) CENTROID(2), the coordinates of the centroid.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) centroid(dim_num)
  integer i
  real ( kind = rk ) t(dim_num,3)

  do i = 1, dim_num
    centroid(i) = sum ( t(i,1:3) ) / 3.0D+00
  end do

  return
end
subroutine triangle_centroid_3d ( t, centroid )

!*****************************************************************************80
!
!! TRIANGLE_CENTROID_3D computes the centroid of a triangle in 3D.
!
!  Discussion:
!
!    The centroid of a triangle can also be considered the 
!    center of gravity or center of mass, assuming that the triangle 
!    is made of a thin uniform sheet of massy material.
!
!    The centroid of a triangle is the intersection of the medians.
!    A median of a triangle is a line connecting any vertex to the
!    midpoint of the opposite side.
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
!  Parameters:
!
!    Input, real ( kind = rk ) T(3,3), the triangle vertices.
!
!    Output, real ( kind = rk ) CENTROID(3), the coordinates of the centroid.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) centroid(dim_num)
  integer i
  real ( kind = rk ) t(dim_num,3)

  do i = 1, dim_num
    centroid(i) = sum ( t(i,1:3) ) / 3.0D+00
  end do

  return
end
subroutine triangle_circumcenter_2d ( t, pc )

!*****************************************************************************80
!
!! TRIANGLE_CIRCUMCENTER_2D computes the circumcenter of a triangle in 2D.
!
!  Discussion:
!
!    The circumcenter of a triangle is the center of the circumcircle, the
!    circle that passes through the three vertices of the triangle.
!
!    The circumcircle contains the triangle, but it is not necessarily the
!    smallest triangle to do so.
!
!    If all angles of the triangle are no greater than 90 degrees, then
!    the center of the circumscribed circle will lie inside the triangle.
!    Otherwise, the center will lie outside the triangle.
!
!    The circumcenter is the intersection of the perpendicular bisectors
!    of the sides of the triangle.
!
!    In geometry, the circumcenter of a triangle is often symbolized by "O".
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) T(2,3), the triangle vertices.
!
!    Output, real ( kind = rk ) PC(2), the circumcenter of the triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) det
  real ( kind = rk ) f(2)
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) t(dim_num,3)
  real ( kind = rk ) top(dim_num)

  f(1) = ( t(1,2) - t(1,1) )**2 + ( t(2,2) - t(2,1) )**2
  f(2) = ( t(1,3) - t(1,1) )**2 + ( t(2,3) - t(2,1) )**2
  
  top(1) =    ( t(2,3) - t(2,1) ) * f(1) - ( t(2,2) - t(2,1) ) * f(2)
  top(2) =  - ( t(1,3) - t(1,1) ) * f(1) + ( t(1,2) - t(1,1) ) * f(2)

  det  =    ( t(2,3) - t(2,1) ) * ( t(1,2) - t(1,1) ) &
          - ( t(2,2) - t(2,1) ) * ( t(1,3) - t(1,1) ) 

  pc(1:2) = t(1:2,1) + 0.5D+00 * top(1:2) / det

  return
end
subroutine triangle_circumcenter_2d_2 ( t, pc )

!*****************************************************************************80
!
!! TRIANGLE_CIRCUMCENTER_2D_2 computes the circumcenter of a triangle in 2D.
!
!  Discussion:
!
!    The circumcenter of a triangle is the center of the circumcircle, the
!    circle that passes through the three vertices of the triangle.
!
!    The circumcircle contains the triangle, but it is not necessarily the
!    smallest triangle to do so.
!
!    If all angles of the triangle are no greater than 90 degrees, then
!    the center of the circumscribed circle will lie inside the triangle.
!    Otherwise, the center will lie outside the triangle.
!
!    The circumcenter is the intersection of the perpendicular bisectors
!    of the sides of the triangle.
!
!    Surprisingly, the diameter of the circle can be found by solving
!    a 2 by 2 linear system.  If we label the vertices of the triangle
!    P1, P2 and P3, then the vectors P2 - P1 and P3 - P1 are secants of
!    the circle, and each forms a right triangle with the diameter 
!    vector through P1. 
! 
!    Hence, the dot product of P2 - P1 with the diameter vector is equal 
!    to the square of the length of P2 - P1, and similarly for P3 - P1.
!    This determines the diameter vector originating at P1.
!
!    In geometry, the circumcenter of a triangle is often symbolized by "O".
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
!  Reference:
!
!    Adrian Bowyer, John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983,
!    ISBN: 0408012420.
!
!  Parameters:
!
!    Input, real ( kind = rk ) T(2,3), the triangle vertices.
!
!    Output, real ( kind = rk ) PC(2), the circumcenter of the triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2
  integer, parameter :: rhs_num = 1

  real ( kind = rk ) a(dim_num,dim_num+rhs_num)
  integer info
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) t(dim_num,3)
!
!  Set up the linear system.
!
  a(1,1) = t(1,2) - t(1,1)
  a(1,2) = t(2,2) - t(2,1)
  a(1,3) = ( t(1,2) - t(1,1) )**2 + ( t(2,2) - t(2,1) )**2

  a(2,1) = t(1,3) - t(1,1)
  a(2,2) = t(2,3) - t(2,1)
  a(2,3) = ( t(1,3) - t(1,1) )**2 + ( t(2,3) - t(2,1) )**2
!
!  Solve the linear system.
!
  call r8mat_solve ( dim_num, rhs_num, a, info )
!
!  Compute the center
!
  if ( info /= 0 ) then
    pc(1:dim_num) = 0.0D+00
  else
    pc(1:dim_num) = t(1:dim_num,1) + 0.5D+00 * a(1:dim_num,dim_num+1)
  end if

  return
end
subroutine triangle_circumcenter ( n, t, p )

!*****************************************************************************80
!
!! TRIANGLE_CIRCUMCENTER computes the circumcenter of a triangle in ND.
!
!  Discussion:
!
!    Three ND points A, B and C lie on a circle.  
!
!    The circumcenter P has the formula
!
!      P = ( Area ( PBC ) * A + Area ( APC) * B + Area ( ABP ) * C )
!        / ( Area ( PBC )     + Area ( APC )    + Area ( ABP ) )
!
!    The details of the formula rely on information supplied
!    by Oscar Lanzi III.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    28 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the spatial dimension.
!
!    Input, real ( kind = rk ) T(N,3), the triangle vertices.
!
!    Output, real ( kind = rk ) P(N), the circumcenter of the triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a
  real ( kind = rk ) abp
  real ( kind = rk ) apc
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) p(n)
  real ( kind = rk ) pbc
  real ( kind = rk ) r8vec_normsq_affine
  real ( kind = rk ) t(n,3)

  a = r8vec_normsq_affine ( n, t(1:n,2), t(1:n,3) )
  b = r8vec_normsq_affine ( n, t(1:n,3), t(1:n,1) )
  c = r8vec_normsq_affine ( n, t(1:n,1), t(1:n,2) )

  pbc = a * ( - a + b + c )
  apc = b * (   a - b + c ) 
  abp = c * (   a + b - c )

  p(1:n) = ( pbc * t(1:n,1) + apc * t(1:n,2) + abp * t(1:n,3) ) &
         / ( pbc            + apc            + abp )

  return
end
subroutine triangle_circumcircle_2d ( t, r, pc )

!*****************************************************************************80
!
!! TRIANGLE_CIRCUMCIRCLE_2D computes the circumcircle of a triangle in 2D.
!
!  Discussion:
!
!    The circumcenter of a triangle is the center of the circumcircle, the
!    circle that passes through the three vertices of the triangle.
!
!    The circumcircle contains the triangle, but it is not necessarily the
!    smallest triangle to do so.
!
!    If all angles of the triangle are no greater than 90 degrees, then
!    the center of the circumscribed circle will lie inside the triangle.
!    Otherwise, the center will lie outside the triangle.
!
!    The circumcenter is the intersection of the perpendicular bisectors
!    of the sides of the triangle.
!
!    In geometry, the circumcenter of a triangle is often symbolized by "O".
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) T(2,3), the triangle vertices.
!
!    Output, real ( kind = rk ) R, PC(2), the circumradius and circumcenter
!    of the triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) bot
  real ( kind = rk ) c
  real ( kind = rk ) det
  real ( kind = rk ) f(2)
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) r
  real ( kind = rk ) top(dim_num)
  real ( kind = rk ) t(dim_num,3)
!
!  Circumradius.
!
  a = sqrt ( ( t(1,2) - t(1,1) )**2 + ( t(2,2) - t(2,1) )**2 )
  b = sqrt ( ( t(1,3) - t(1,2) )**2 + ( t(2,3) - t(2,2) )**2 )
  c = sqrt ( ( t(1,1) - t(1,3) )**2 + ( t(2,1) - t(2,3) )**2 )

  bot = ( a + b + c ) * ( - a + b + c ) * (   a - b + c ) * (   a + b - c )

  if ( bot <= 0.0D+00 ) then
    r = -1.0D+00
    pc(1:2) = 0.0D+00
    return
  end if

  r = a * b * c / sqrt ( bot )
!
!  Circumcenter.
! 
  f(1) = ( t(1,2) - t(1,1) )**2 + ( t(2,2) - t(2,1) )**2
  f(2) = ( t(1,3) - t(1,1) )**2 + ( t(2,3) - t(2,1) )**2
  
  top(1) =    ( t(2,3) - t(2,1) ) * f(1) - ( t(2,2) - t(2,1) ) * f(2)
  top(2) =  - ( t(1,3) - t(1,1) ) * f(1) + ( t(1,2) - t(1,1) ) * f(2)

  det  =    ( t(2,3) - t(2,1) ) * ( t(1,2) - t(1,1) ) &
          - ( t(2,2) - t(2,1) ) * ( t(1,3) - t(1,1) ) 

  pc(1:2) = t(1:2,1) + 0.5D+00 * top(1:2) / det

  return
end
subroutine triangle_circumcircle_2d_2 ( t, r, pc )

!*****************************************************************************80
!
!! TRIANGLE_CIRCUMCIRCLE_2D_2 computes the circumcircle of a triangle in 2D.
!
!  Discussion:
!
!    The circumscribed circle of a triangle is the circle that passes through
!    the three vertices of the triangle.  The circumscribed circle contains
!    the triangle, but it is not necessarily the smallest triangle to do so.
!
!    Surprisingly, the diameter of the circle can be found by solving
!    a 2 by 2 linear system.  This is because the vectors P2 - P1
!    and P3 - P1 are secants of the circle, and each forms a right
!    triangle with the diameter.  Hence, the dot product of
!    P2 - P1 with the diameter is equal to the square of the length
!    of P2 - P1, and similarly for P3 - P1.  This determines the
!    diameter vector originating at P1.
!
!    If all angles of the triangle are no greater than 90 degrees, then
!    the center of the circumscribed circle will lie inside the triangle.
!    Otherwise, the center will lie outside the triangle.
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
!  Reference:
!
!    Adrian Bowyer, John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983,
!    ISBN: 0408012420.
!
!  Parameters:
!
!    Input, real ( kind = rk ) T(2,3), the triangle vertices.
!
!    Output, real ( kind = rk ) R, PC(2), the circumradius and circumcenter.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2
  integer, parameter :: rhs_num = 1

  real ( kind = rk ) a(dim_num,dim_num+rhs_num)
  integer info
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) r
  real ( kind = rk ) t(dim_num,3)
!
!  Set up the linear system.
!
  a(1,1) = t(1,2) - t(1,1)
  a(1,2) = t(2,2) - t(2,1)
  a(1,3) = ( t(1,2) - t(1,1) )**2 + ( t(2,2) - t(2,1) )**2

  a(2,1) = t(1,3) - t(1,1)
  a(2,2) = t(2,3) - t(2,1)
  a(2,3) = ( t(1,3) - t(1,1) )**2 + ( t(2,3) - t(2,1) )**2
!
!  Solve the linear system.
!
  call r8mat_solve ( dim_num, rhs_num, a, info )

  if ( info /= 0 ) then
    r = -1.0D+00
    pc(1:dim_num) = 0.0D+00
  end if

  r = 0.5D+00 * sqrt ( a(1,dim_num+1) * a(1,dim_num+1) &
                     + a(2,dim_num+1) * a(2,dim_num+1) )
  pc(1:dim_num) = t(1:dim_num,1) + 0.5D+00 * a(1:dim_num,dim_num+1)

  return
end
subroutine triangle_circumradius_2d ( t, r )

!*****************************************************************************80
!
!! TRIANGLE_CIRCUMRADIUS_2D computes the circumradius of a triangle in 2D.
!
!  Discussion:
!
!    The circumscribed circle of a triangle is the circle that passes through
!    the three vertices of the triangle.  The circumscribed circle contains
!    the triangle, but it is not necessarily the smallest triangle to do so.
!
!    The circumradius of a triangle is the radius of the circumscribed
!    circle.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) T(2,3), the triangle vertices.
!
!    Output, real ( kind = rk ) R, the circumradius of the circumscribed circle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) bot
  real ( kind = rk ) c
  real ( kind = rk ) r
  real ( kind = rk ) t(dim_num,3)
!
!  Compute the length of each side.
!
  a = sqrt ( sum ( ( t(1:dim_num,1) - t(1:dim_num,2) )**2 ) )
  b = sqrt ( sum ( ( t(1:dim_num,2) - t(1:dim_num,3) )**2 ) )
  c = sqrt ( sum ( ( t(1:dim_num,3) - t(1:dim_num,1) )**2 ) )

  bot = ( a + b + c ) * ( - a + b + c ) * (   a - b + c ) * (   a + b - c )

  if ( bot <= 0.0D+00 ) then
    r = -1.0D+00
    return
  end if

  r = a * b * c / sqrt ( bot )

  return
end
subroutine triangle_contains_line_exp_3d ( t, p1, p2, inside, pint )

!*****************************************************************************80
!
!! TRIANGLE_CONTAINS_LINE_EXP_3D finds if a line is inside a triangle in 3D.
!
!  Discussion:
!
!    A line will "intersect" the plane of a triangle in 3D if 
!    * the line does not lie in the plane of the triangle
!      (there would be infinitely many intersections), AND 
!    * the line does not lie parallel to the plane of the triangle
!      (there are no intersections at all).
!
!    Therefore, if a line intersects the plane of a triangle, it does so
!    at a single point.  We say the line is "inside" the triangle if, 
!    regarded as 2D objects, the intersection point of the line and the plane
!    is inside the triangle.
!
!    A triangle in 3D is determined by three points:
!
!      T(1:3,1), T(1:3,2) and T(1:3,3).
!
!    The explicit form of a line in 3D is:
!
!      the line through the points P1(1:3), P2(1:3).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Steve Marschner, Cornell University,
!    CS465 Notes: Simple Ray-Triangle Intersection.
!
!  Parameters:
!
!    Input, real ( kind = rk ) T(3,3), the triangle vertices.
!
!    Input, real ( kind = rk ) P1(3), P2(3), two points on the line.
!
!    Output, logical INSIDE, is TRUE if (the intersection point of)
!    the line is inside the triangle.
!
!    Output, real ( kind = rk ) PINT(3), the point where the line
!    intersects the plane of the triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  logical inside
  integer ival
  logical line_exp_is_degenerate_nd
  real ( kind = rk ) normal(dim_num)
  real ( kind = rk ) normal2(dim_num)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) pint(dim_num)
  real ( kind = rk ) t(dim_num,3)
  real ( kind = rk ) temp
  logical triangle_is_degenerate_nd
  real ( kind = rk ) v1(dim_num)
  real ( kind = rk ) v2(dim_num)
!
!  Make sure the line is not degenerate.
!
  if ( line_exp_is_degenerate_nd ( dim_num, p1, p2 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGLE_CONTAINS_LINE_EXP_3D - Fatal error!'
    write ( *, '(a)' ) '  The explicit line is degenerate.'
    stop 1
  end if
!
!  Make sure the triangle is not degenerate.
!
  if ( triangle_is_degenerate_nd ( dim_num, t ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGLE_CONTAINS_LINE_EXP_3D - Fatal error!'
    write ( *, '(a)' ) '  The triangle is degenerate.'
    stop 1
  end if
!
!  Determine a unit normal vector associated with the plane of
!  the triangle.
!
  v1(1:dim_num) = t(1:dim_num,2) - t(1:dim_num,1)
  v2(1:dim_num) = t(1:dim_num,3) - t(1:dim_num,1)

  normal(1) = v1(2) * v2(3) - v1(3) * v2(2)
  normal(2) = v1(3) * v2(1) - v1(1) * v2(3)
  normal(3) = v1(1) * v2(2) - v1(2) * v2(1)

  temp = sqrt ( sum ( normal(1:dim_num)**2 ) )
  normal(1:dim_num) = normal(1:dim_num) / temp
!
!  Find the intersection of the plane and the line.
!
  call plane_normal_line_exp_int_3d ( t(1:dim_num,1), normal, p1, p2, &
    ival, pint )

  if ( ival == 0 ) then
    inside = .false.
    pint(1:dim_num) = huge ( temp )
    return
  else if ( ival == 2 ) then
    inside = .false.
    pint(1:dim_num) = p1(1:dim_num)
    return
  end if
!
!  Now, check that all three triangles made by two vertices and
!  the intersection point have the same "clock sense" as the
!  triangle's normal vector.
!
  v1(1:dim_num) = t(1:dim_num,2)  - t(1:dim_num,1)
  v2(1:dim_num) = pint(1:dim_num) - t(1:dim_num,1)

  normal2(1) = v1(2) * v2(3) - v1(3) * v2(2)
  normal2(2) = v1(3) * v2(1) - v1(1) * v2(3)
  normal2(3) = v1(1) * v2(2) - v1(2) * v2(1)

  if ( dot_product ( normal(1:dim_num), normal2(1:dim_num) ) < 0.0D+00 ) then
    inside = .false.
    return
  end if

  v1(1:dim_num) = t(1:dim_num,3)  - t(1:dim_num,2)
  v2(1:dim_num) = pint(1:dim_num) - t(1:dim_num,2)

  normal2(1) = v1(2) * v2(3) - v1(3) * v2(2)
  normal2(2) = v1(3) * v2(1) - v1(1) * v2(3)
  normal2(3) = v1(1) * v2(2) - v1(2) * v2(1)

  if ( dot_product ( normal(1:dim_num), normal2(1:dim_num) ) < 0.0D+00 ) then
    inside = .false.
    return
  end if

  v1(1:dim_num) = t(1:dim_num,1)  - t(1:dim_num,3)
  v2(1:dim_num) = pint(1:dim_num) - t(1:dim_num,3)

  normal2(1) = v1(2) * v2(3) - v1(3) * v2(2)
  normal2(2) = v1(3) * v2(1) - v1(1) * v2(3)
  normal2(3) = v1(1) * v2(2) - v1(2) * v2(1)

  if ( dot_product ( normal(1:dim_num), normal2(1:dim_num) ) < 0.0D+00 ) then
    inside = .false.
    return
  end if

  inside = .true.

  return
end
subroutine triangle_contains_line_par_3d ( t, p0, pd, inside, p )

!*****************************************************************************80
!
!! TRIANGLE_CONTAINS_LINE_PAR_3D: finds if a line is inside a triangle in 3D.
!
!  Discussion:
!
!    A line will "intersect" the plane of a triangle in 3D if 
!    * the line does not lie in the plane of the triangle
!      (there would be infinitely many intersections), AND 
!    * the line does not lie parallel to the plane of the triangle
!      (there are no intersections at all).
!
!    Therefore, if a line intersects the plane of a triangle, it does so
!    at a single point.  We say the line is "inside" the triangle if, 
!    regarded as 2D objects, the intersection point of the line and the plane
!    is inside the triangle.
!
!    A triangle in 3D is determined by three points:
!
!      T(1:3,1), T(1:3,2) and T(1:3,3).
!
!    The parametric form of a line in 3D is:
!
!      P(1:3) = P0(1:3) + PD(1:3) * T
!
!    We can normalize by requiring PD to have euclidean norm 1, 
!    and the first nonzero entry positive.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 February 2007
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
!    ISBN: 0408012420,
!    page 111.
!
!  Parameters:
!
!    Input, real ( kind = rk ) T(3,3), the three points that define
!    the triangle.
!
!    Input, real ( kind = rk ) P0(3), PD(3), parameters that define the
!    parametric line.
!
!    Output, logical INSIDE, is TRUE if (the intersection point of)
!    the line is inside the triangle.
!
!    Output, real ( kind = rk ) P(3), is the point of intersection of the line
!    and the plane of the triangle, unless they are parallel.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) a
  real ( kind = rk ) angle_sum
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) d
  real ( kind = rk ) denom
  logical inside
  logical intersect
  real ( kind = rk ) norm
  real ( kind = rk ) norm1
  real ( kind = rk ) norm2
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) p0(dim_num)
  real ( kind = rk ) pd(dim_num)
  real ( kind = rk ) r8_acos
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) t(dim_num,3)
  real ( kind = rk ) t_int
  real ( kind = rk ), parameter :: tol = 0.00001D+00
  real ( kind = rk ) v1(dim_num)
  real ( kind = rk ) v2(dim_num)
  real ( kind = rk ) v3(dim_num)
!
!  Determine the implicit form (A,B,C,D) of the plane containing the
!  triangle.
!
  a = ( t(2,2) - t(2,1) ) * ( t(3,3) - t(3,1) ) &
    - ( t(3,2) - t(3,1) ) * ( t(2,3) - t(2,1) )

  b = ( t(3,2) - t(3,1) ) * ( t(1,3) - t(1,1) ) &
    - ( t(1,2) - t(1,1) ) * ( t(3,3) - t(3,1) )

  c = ( t(1,2) - t(1,1) ) * ( t(2,3) - t(2,1) ) &
    - ( t(2,2) - t(2,1) ) * ( t(1,3) - t(1,1) )

  d = - t(1,2) * a - t(2,2) * b - t(3,2) * c
!
!  Make sure the plane is well-defined.
!
  norm1 = sqrt ( a * a + b * b + c * c )

  if ( norm1 == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGLE_LINE_PAR_INT_3D - Fatal error!'
    write ( *, '(a)' ) '  The plane normal vector is null.'
    inside = .false.
    p(1:dim_num) = 0.0D+00
    stop 1
  end if
!
!  Make sure the implicit line is well defined.
!
  norm2 = sqrt ( sum ( pd(1:dim_num)**2 ) )

  if ( norm2 == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGLE_LINE_PAR_INT_3D - Fatal error!'
    write ( *, '(a)' ) '  The line direction vector is null.'
    inside = .false.
    p(1:dim_num) = 0.0D+00
    stop 1
  end if
!
!  Determine the denominator of the parameter in the
!  implicit line definition that determines the intersection
!  point.
!
  denom = a * pd(1) + b * pd(2) + c * pd(3)
!
!  If DENOM is zero, or very small, the line and the plane may be 
!  parallel or almost so.
!
  if ( abs ( denom ) < tol * norm1 * norm2 ) then
!
!  The line may actually lie in the plane.  We're not going
!  to try to address this possibility.
!
    if ( a * p0(1) + b * p0(2) + c * p0(3) + d == 0.0D+00 ) then

      intersect = .true.
      inside = .false.
      p(1:dim_num) = p0(1:dim_num)
!
!  The line and plane are parallel and disjoint.
!
    else

      intersect = .false.
      inside = .false.
      p(1:dim_num) = 0.0D+00

    end if
!
!  The line and plane intersect at a single point P.
!
  else

    intersect = .true.
    t_int = - ( a * p0(1) + b * p0(2) + c * p0(3) + d ) / denom
    p(1:dim_num) = p0(1:dim_num) + t_int * pd(1:dim_num)
!
!  To see if P is included in the triangle, sum the angles
!  formed by P and pairs of the vertices.  If the point is in the
!  triangle, we get a total 360 degree view.  Otherwise, we
!  get less than 180 degrees.
!
    v1(1:dim_num) = t(1:dim_num,1) - p(1:dim_num)
    v2(1:dim_num) = t(1:dim_num,2) - p(1:dim_num)
    v3(1:dim_num) = t(1:dim_num,3) - p(1:dim_num)

    norm = sqrt ( sum ( v1(1:dim_num)**2 ) )

    if ( norm == 0.0D+00 ) then
      inside = .true.
      return
    end if

    v1(1:dim_num) = v1(1:dim_num) / norm

    norm = sqrt ( sum ( v2(1:dim_num)**2 ) )

    if ( norm == 0.0D+00 ) then
      inside = .true.
      return
    end if

    v2(1:dim_num) = v2(1:dim_num) / norm

    norm = sqrt ( sum ( v3(1:dim_num)**2 ) )

    if ( norm == 0.0D+00 ) then
      inside = .true.
      return
    end if

    v3(1:dim_num) = v3(1:dim_num) / norm

    angle_sum = r8_acos ( dot_product ( v1(1:3), v2(1:3) ) ) &
              + r8_acos ( dot_product ( v2(1:3), v3(1:3) ) ) &
              + r8_acos ( dot_product ( v3(1:3), v1(1:3) ) )

    if ( nint ( angle_sum / r8_pi ) == 2 ) then
      inside = .true.
    else
      inside = .false.
    end if

  end if

  return
end
subroutine triangle_contains_point_2d_1 ( t, p, inside )

!*****************************************************************************80
!
!! TRIANGLE_CONTAINS_POINT_2D_1 finds if a point is inside a triangle in 2D.
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
!  Parameters:
!
!    Input, real ( kind = rk ) T(2,3), the triangle vertices.
!
!    Input, real ( kind = rk ) P(2), the point to be checked.
!
!    Output, logical INSIDE, is TRUE if the point is inside 
!    the triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  logical inside
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) t(dim_num,3)
  real ( kind = rk ) xsi(dim_num+1)

  call triangle_barycentric_2d ( t, p, xsi )

  if ( any ( xsi(1:3) < 0.0D+00 ) ) then
    inside = .false.
  else
    inside = .true.
  end if

  return
end
subroutine triangle_contains_point_2d_2 ( t, p, inside )

!*****************************************************************************80
!
!! TRIANGLE_CONTAINS_POINT_2D_2 finds if a point is inside a triangle in 2D.
!
!  Discussion:
!
!    The routine assumes that the vertices are given in counter clockwise
!    order.  If the triangle vertices are actually given in clockwise 
!    order, this routine will behave as though the triangle contains
!    no points whatsoever!
!
!    The routine determines if a point P is "to the right of" each of the lines
!    that bound the triangle.  It does this by computing the cross product
!    of vectors from a vertex to its next vertex, and to P.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) T(2,3), the triangle vertices.
!    The vertices should be given in counter clockwise order.
!
!    Input, real ( kind = rk ) P(2), the point to be checked.
!
!    Output, logical INSIDE, is TRUE if the point is 
!    inside the triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  logical inside
  integer j
  integer k
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) t(dim_num,3)

  do j = 1, 3

    k = mod ( j, 3 ) + 1

    if ( 0.0D+00 < ( p(1) - t(1,j) ) * ( t(2,k) - t(2,j) ) &
                 - ( p(2) - t(2,j) ) * ( t(1,k) - t(1,j) ) ) then
      inside = .false.
      return
    end if

  end do

  inside = .true.

  return
end
subroutine triangle_contains_point_2d_3 ( t, p, inside )

!*****************************************************************************80
!
!! TRIANGLE_CONTAINS_POINT_2D_3 finds if a point is inside a triangle in 2D.
!
!  Discussion:
!
!    This routine is the same as TRIANGLE_CONTAINS_POINT_2D_2, except
!    that it does not assume an ordering of the points.  It should
!    work correctly whether the vertices of the triangle are listed
!    in clockwise or counter clockwise order.
!
!    The routine determines if a point P is "to the right of" each of the lines
!    that bound the triangle.  It does this by computing the cross product
!    of vectors from a vertex to its next vertex, and to P.
!
!    The point is inside the triangle if it is to the right of all
!    the lines, or to the left of all the lines.
!
!    This version was suggested by Paulo Ernesto of Maptek Brasil.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) T(2,3), the triangle vertices.
!
!    Input, real ( kind = rk ) P(2), the point to be checked.
!
!    Output, logical INSIDE, is TRUE if the point is 
!    inside the triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) dir_new
  real ( kind = rk ) dir_old
  logical inside
  integer j
  integer k
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) t(dim_num,3)

  dir_old = 0.0D+00

  do j = 1, 3

    k = mod ( j, 3 ) + 1

    dir_new = ( p(1) - t(1,j) ) * ( t(2,k) - t(2,j) ) &
            - ( p(2) - t(2,j) ) * ( t(1,k) - t(1,j) )
    
    if ( dir_new * dir_old < 0.0D+00 ) then
      inside = .false.
      return
    end if

    if ( dir_new /= 0.0D+00 ) then
      dir_old = dir_new
    end if

  end do

  inside = .true.

  return
end
subroutine triangle_diameter_2d ( t, diameter )

!*****************************************************************************80
!
!! TRIANGLE_DIAMETER_2D computes the diameter of a triangle in 2D.
!
!  Discussion:
!
!    The diameter of a triangle is the diameter of the smallest circle
!    that can be drawn around the triangle.  At least two of the vertices
!    of the triangle will intersect the circle, but not necessarily
!    all three!
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) T(2,3), the triangle vertices.
!
!    Output, real ( kind = rk ) DIAMETER, the diameter of the triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) a
  real ( kind = rk ) asq
  real ( kind = rk ) b
  real ( kind = rk ) bsq
  real ( kind = rk ) c
  real ( kind = rk ) csq
  real ( kind = rk ) diameter
  real ( kind = rk ) t(dim_num,3)
!
!  Compute the squared length of each side.
!
  asq = sum ( t(1:dim_num,1) - t(1:dim_num,2) )**2
  bsq = sum ( t(1:dim_num,2) - t(1:dim_num,3) )**2
  csq = sum ( t(1:dim_num,3) - t(1:dim_num,1) )**2
!
!  Take care of a zero side.
!
  if ( asq == 0.0D+00 ) then
    diameter = sqrt ( bsq )
    return
  else if ( bsq == 0.0D+00 ) then
    diameter = sqrt ( csq )
    return
  else if ( csq == 0.0D+00 ) then
    diameter = sqrt ( asq )
    return
  end if
!
!  Make ASQ the largest.
!
  if ( asq < bsq ) then
    call r8_swap ( asq, bsq )
  end if

  if ( asq < csq ) then
    call r8_swap ( asq, csq )
  end if
!
!  If ASQ is very large...
!
  if ( bsq + csq < asq ) then

    diameter = sqrt ( asq )

  else

    a = sqrt ( asq )
    b = sqrt ( bsq )
    c = sqrt ( csq )

    diameter = 2.0D+00 * a * b * c / sqrt ( ( a + b + c ) * ( - a + b + c ) &
      * ( a - b + c ) * ( a + b - c ) )

  end if

  return
end
subroutine triangle_edge_length_2d ( t, edge_length )

!*****************************************************************************80
!
!! TRIANGLE_EDGE_LENGTH_2D returns edge lengths of a triangle in 2D.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) T(2,3), the triangle vertices.
!
!    Output, real ( kind = rk ) EDGE_LENGTH(3), the length of the edges.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) edge_length(3)
  integer i4_wrap
  integer j1
  integer j2
  real ( kind = rk ) r8vec_norm
  real ( kind = rk ) t(dim_num,3)

  do j1 = 1, 3
    j2 = i4_wrap ( j1 + 1, 1, 3 )
    edge_length(j1) = &
      r8vec_norm ( dim_num, t(1:dim_num,j2) - t(1:dim_num,j1) )
  end do

  return
end
subroutine triangle_gridpoints_2d ( t, sub_num, grid_max, grid_num, g )

!*****************************************************************************80
!
!! TRIANGLE_GRIDPOINTS_2D computes gridpoints within a triangle in 2D.
!
!  Discussion:
!
!    The gridpoints are computed by repeated halving of the triangle.
!    The 0-th set of grid points is the vertices themselves.
!    The first set of grid points is the midpoints of the sides.
!    These points can be used to draw 4 triangles that make up the original
!    triangle.  The second set of grid points is the side midpoints and centers
!    of these four triangles.
!
!     SUB_NUM                     GRID_NUM
!    -----                        -----
!        0      1                  =  1  (centroid)
!        1      1 + 2              =  3  (vertices)
!        2      1 + 2 + 3          =  6
!        3      1 + 2 + 3 + 4      = 10
!        4      1 + 2 + 3 + 4 + 5  = 15
!
!    GRID_NUM is the sum of the integers from 1 to SUB_NUM+1 or
!
!      GRID_NUM = (SUB_NUM+1) * (SUB_NUM+2) / 2
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
!  Parameters:
!
!    Input, real ( kind = rk ) T(2,3), the triangle vertices.
!
!    Input, integer SUB_NUM, the number of subdivisions.
!
!    Input, integer GRID_MAX, the maximum number of grid points.
!
!    Output, integer GRID_NUM, the number of grid points returned.
!
!    Output, real ( kind = rk ) G(2,GRID_MAX), the grid points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer grid_max
  integer, parameter :: dim_num = 2

  real ( kind = rk ) g(dim_num,grid_max)
  integer i
  integer j
  integer grid_num
  integer sub_num
  real ( kind = rk ) t(dim_num,3)

  grid_num = 0
!
!  Special case, SUB_NUM = 0.
!
  if ( sub_num == 0 ) then
    if ( 1 <= grid_max ) then
      grid_num = 1
      g(1,1) = ( t(1,1) + t(1,2) + t(1,3) ) / 3.0D+00
      g(2,1) = ( t(2,1) + t(2,2) + t(2,3) ) / 3.0D+00
    end if
    return
  end if

  do i = 0, sub_num

    do j = 0, sub_num - i

      if ( grid_num < grid_max ) then

        grid_num = grid_num + 1

        g(1,grid_num) = ( real (           i,     kind = rk ) * t(1,1) &
                        + real (               j, kind = rk ) * t(1,2) &
                        + real ( sub_num - i - j, kind = rk ) * t(1,3) ) &
                        / real ( sub_num,         kind = rk )

        g(2,grid_num) = ( real (           i,     kind = rk ) * t(2,1) &
                        + real (               j, kind = rk ) * t(2,2) &
                        + real ( sub_num - i - j, kind = rk ) * t(2,3) ) &
                        / real ( sub_num,         kind = rk )
      end if

    end do
  end do

  return
end
subroutine triangle_incenter_2d ( t, pc )

!*****************************************************************************80
!
!! TRIANGLE_INCENTER_2D computes the incenter of a triangle in 2D.
!
!  Discussion:
!
!    The incenter of a triangle is the center of the inscribed circle.
!
!    The inscribed circle of a triangle is the largest circle that can
!    be drawn inside the triangle.
!
!    The inscribed circle is tangent to all three sides of the triangle.
!
!    The angle bisectors of the triangle intersect at the center of the
!    inscribed circle.
!
!    In geometry, the incenter is often represented by "I".
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 August 2003
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
!  Parameters:
!
!    Input, real ( kind = rk ) T(2,3), the triangle vertices.
!
!    Output, real ( kind = rk ) PC(2), the incenter.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) perimeter
  real ( kind = rk ) t(dim_num,3)
!
!  Compute the length of each side.
!
  a = sqrt ( sum ( ( t(1:dim_num,1) - t(1:dim_num,2) )**2 ) )
  b = sqrt ( sum ( ( t(1:dim_num,2) - t(1:dim_num,3) )**2 ) )
  c = sqrt ( sum ( ( t(1:dim_num,3) - t(1:dim_num,1) )**2 ) )

  perimeter = a + b + c

  if ( perimeter == 0.0D+00 ) then
    pc(1:dim_num) = t(1:dim_num,1)
  else
    pc(1:dim_num) = ( b * t(1:dim_num,1) &
                    + c * t(1:dim_num,2) &
                    + a * t(1:dim_num,3) ) / perimeter
  end if

  return
end
subroutine triangle_incircle_2d ( t, r, pc )

!*****************************************************************************80
!
!! TRIANGLE_INCIRCLE_2D computes the inscribed circle of a triangle in 2D.
!
!  Discussion:
!
!    The inscribed circle of a triangle is the largest circle that can
!    be drawn inside the triangle.  It is tangent to all three sides,
!    and the lines from its center to the vertices bisect the angles
!    made by each vertex.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 December 2004
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
!  Parameters:
!
!    Input, real ( kind = rk ) T(2,3), the triangle vertices.
!
!    Output, real ( kind = rk ) R, PC(2), the radius and center of the
!    inscribed circle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) perimeter
  real ( kind = rk ) r
  real ( kind = rk ) t(dim_num,3)
!
!  Compute the length of each side.
!
  a = sqrt ( sum ( ( t(1:dim_num,1) - t(1:dim_num,2) )**2 ) )
  b = sqrt ( sum ( ( t(1:dim_num,2) - t(1:dim_num,3) )**2 ) )
  c = sqrt ( sum ( ( t(1:dim_num,3) - t(1:dim_num,1) )**2 ) )

  perimeter = a + b + c

  if ( perimeter == 0.0D+00 ) then
    pc(1:dim_num) = t(1:dim_num,1)
    r = 0.0D+00
    return
  end if

  pc(1:dim_num) = (  &
      b * t(1:dim_num,1) &
    + c * t(1:dim_num,2) &
    + a * t(1:dim_num,3) ) / perimeter

  r = 0.5D+00 * sqrt ( &
      ( - a + b + c )  &
    * ( + a - b + c )  &
    * ( + a + b - c ) / perimeter )

  return
end
subroutine triangle_inradius_2d ( t, r )

!*****************************************************************************80
!
!! TRIANGLE_INRADIUS_2D: radius of the inscribed circle of a triangle in 2D.
!
!  Discussion:
!
!    The inscribed circle of a triangle is the largest circle that can
!    be drawn inside the triangle.  It is tangent to all three sides,
!    and the lines from its center to the vertices bisect the angles
!    made by each vertex.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    13 April 2005
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
!  Parameters:
!
!    Input, real ( kind = rk ) T(2,3), the triangle vertices.
!
!    Output, real ( kind = rk ) R, the radius of the inscribed circle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) perimeter
  real ( kind = rk ) r
  real ( kind = rk ) t(dim_num,3)
!
!  Compute the length of each side.
!
  a = sqrt ( sum ( ( t(1:dim_num,1) - t(1:dim_num,2) )**2 ) )
  b = sqrt ( sum ( ( t(1:dim_num,2) - t(1:dim_num,3) )**2 ) )
  c = sqrt ( sum ( ( t(1:dim_num,3) - t(1:dim_num,1) )**2 ) )

  perimeter = a + b + c

  if ( perimeter == 0.0D+00 ) then
    r = 0.0D+00
    return
  end if

  r = 0.5D+00 * sqrt ( &
      ( - a + b + c )  &
    * ( + a - b + c )  &
    * ( + a + b - c ) / perimeter )

  return
end
function triangle_is_degenerate_nd ( dim_num, t )

!*****************************************************************************80
!
!! TRIANGLE_IS_DEGENERATE_ND finds if a triangle is degenerate in ND.
!
!  Discussion:
!
!    A triangle in ND is described by the coordinates of its 3 vertices.
!
!    A triangle in ND is degenerate if any two vertices are equal.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = rk ) T(DIM_NUM,3), the triangle vertices.
!
!    Output, logical TRIANGLE_IS_DEGENERATE_ND, is TRUE if the
!    triangle is degenerate.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num

  real ( kind = rk ) t(dim_num,3)
  logical triangle_is_degenerate_nd

  triangle_is_degenerate_nd = &
     ( all ( t(1:dim_num,1) == t(1:dim_num,2) ) .or. &
       all ( t(1:dim_num,2) == t(1:dim_num,3) ) .or. &
       all ( t(1:dim_num,3) == t(1:dim_num,1) ) )

  return
end
subroutine triangle_lattice_layer_point_next ( c, v, more )

!*****************************************************************************80
!
!! TRIANGLE_LATTICE_LAYER_POINT_NEXT: next triangle lattice layer point.
!
!  Discussion:
!
!    The triangle lattice layer L is bounded by the lines
!
!      0 <= X,
!      0 <= Y,
!      L - 1 < X / C(1) + Y / C(2) <= L.
!
!    In particular, layer L = 0 always contains the single point (0,0).
!
!    This function returns, one at a time, the points that lie within 
!    a given triangle lattice layer.
!
!    Thus, if we set C(1) = 2, C(2) = 3, then we get the following layers:
!
!    L = 0: (0,0)
!    L = 1: (1,0), (2,0), (0,1), (1,1), (0,2), (0,3)
!    L = 2: (3,0), (4,0), (2,1), (3,1), (1,2), (2,2), (1,3), (2,3),
!           (0,4), (1,4), (0,5), (0,6).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    08 July 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer C(3), coefficients defining the 
!    lattice layer.  Entry C(3) contains the layer index.
!    C(1) and C(2) should be positive, and C(3) must be nonnegative.
!
!    Input/output, integer V(2).  On first call for a given layer,
!    the input value of V is not important.  On a repeated call for the same
!    layer, the input value of V should be the output value from the previous 
!    call.  On output, V contains the next lattice layer point.
!
!    Input/output, logical MORE.  On input, set MORE to FALSE to 
!    indicate that this is the first call for a given layer.  Thereafter, the 
!    input value should be the output value from the previous call.  On output,
!    MORE is TRUE if the returned value V is a new point.
!    If the output value is FALSE, then no more points were found,
!    and V was reset to 0, and the lattice layer has been exhausted.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer c(3)
  integer c1n
  integer i4vec_lcm
  logical more
  integer, parameter :: n = 2
  integer rhs1
  integer rhs2
  integer v(2)
!
!  Treat layer C(N+1) = 0 specially.
!
  if ( c(n+1) == 0 ) then
    if ( .not. more ) then
      v(1:n) = 0
      more = .true.
    else
      more = .false.
    end if
    return
  end if
!
!  Compute first point.
!
  if ( .not. more ) then

    v(1) = ( c(n+1) - 1 ) * c(1) + 1
    v(2) = 0      
    more = .true.

  else

    c1n = i4vec_lcm ( n, c )

    rhs1 = c1n * ( c(n+1) - 1 )
    rhs2 = c1n *   c(n+1)

    if ( c(2) * ( v(1) + 1 ) + c(1) * v(2) <= rhs2 ) then
      v(1) = v(1) + 1
    else
      v(1) = ( rhs1 - c(1) * ( v(2) + 1 ) ) / c(2)
      v(1) = max ( v(1), 0 )
      v(2) = v(2) + 1
      if ( c(2) * v(1) + c(1) * v(2) <= rhs1 ) then
        v(1) = v(1) + 1
      end if
      if ( c(2) * v(1) + c(1) * v(2) <= rhs2 ) then

      else
        v(1:n) = 0
        more = .false.
      end if
    end if    
  end if

  return
end
subroutine triangle_lattice_point_next ( c, v, more )

!*****************************************************************************80
!
!! TRIANGLE_LATTICE_POINT_NEXT returns the next triangle lattice point.
!
!  Discussion:
!
!    The lattice triangle is defined by the vertices:
!
!      (0,0), (C(3)/C(1), 0) and (0,C(3)/C(2))
!
!    The lattice triangle is bounded by the lines
!
!      0 <= X,
!      0 <= Y
!      X / C(1) + Y / C(2) <= C(3)
!
!    Lattice points are listed one at a time, starting at the origin,
!    with X increasing first.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    08 July 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer C(3), coefficients defining the 
!    lattice triangle.  These should be positive.
!
!    Input/output, integer V(2).  On first call, the input
!    value is not important.  On a repeated call, the input value should
!    be the output value from the previous call.  On output, V contains
!    the next lattice point.
!
!    Input/output, logical MORE.  On input, set MORE to FALSE to 
!    indicate that this is the first call for a given triangle.  Thereafter, 
!    the input value should be the output value from the previous call.  On 
!    output, MORE is TRUE if the returned value V is a new lattice point.
!    If the output value is FALSE, then no more lattice points were found,
!    and V was reset to 0, and the routine should not be called further
!    for this triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer c(3)
  integer c1n
  integer i4vec_lcm
  logical more
  integer, parameter :: n = 2
  integer rhs
  integer v(2)

  if ( .not. more ) then

    v(1:n) = 0
    more = .true.

  else

    c1n = i4vec_lcm ( n, c )

    rhs = c1n * c(n+1)

    if ( c(2) * ( v(1) + 1 ) + c(1) * v(2) <= rhs ) then
      v(1) = v(1) + 1
    else
      v(1) = 0
      if ( c(2) * v(1) + c(1) * ( v(2) + 1 ) <= rhs ) then
        v(2) = v(2) + 1
      else
        v(2) = 0
        more = .false.
      end if
    end if    
  end if

  return
end
subroutine triangle_line_imp_int_2d ( t, a, b, c, int_num, pint )

!*****************************************************************************80
!
!! TRIANGLE_LINE_IMP_INT_2D: implicit line intersects a triangle in 2D.
!
!  Discussion:
!
!    An implicit line is the set of points ( X, Y ) satisfying
!
!      A * X + B * Y + C = 0
!
!    where at least one of A and B is not zero.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) T(2,3), the triangle vertices.
!
!    Input, real ( kind = rk ) A, B, C, determine the equation of the line:
!    A*X + B*Y + C = 0.
!
!    Output, integer INT_NUM, the number of points of intersection
!    of the line with the triangle.  INT_NUM may be 0, 1, 2 or 3.
!
!    Output, real ( kind = rk ) PINT(2,3), contains the intersection points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) a
  real ( kind = rk ) a1
  real ( kind = rk ) b
  real ( kind = rk ) b1
  real ( kind = rk ) c
  real ( kind = rk ) c1
  integer i
  integer i4_wrap
  integer int_num
  integer ival
  integer j
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) pint(dim_num,3)
  real ( kind = rk ) t(dim_num,3)
  real ( kind = rk ) test1
  real ( kind = rk ) test2

  int_num = 0

  do i = 1, 3

    j = i4_wrap ( i+1, 1, 3 )
!
!  Get the implicit form of the line through vertices I and I+1.
!
    call line_exp2imp_2d ( t(1:2,i), t(1:2,j), a1, b1, c1 )
!
!  Seek an intersection with the original line.
!
    call lines_imp_int_2d ( a, b, c, a1, b1, c1, ival, p )
!
!  If there is an intersection, determine if it lies between the two vertices.
!
    if ( ival == 1 ) then

      test1 = sum ( ( p(1:dim_num)   - t(1:dim_num,i) ) &
                  * ( t(1:dim_num,j) - t(1:dim_num,i) ) )
      test2 = sum ( ( t(1:dim_num,j) - t(1:dim_num,i) ) &
                  * ( t(1:dim_num,j) - t(1:dim_num,i) ) )

      if ( 0 <= test1 .and. test1 <= test2 ) then
        int_num = int_num + 1
        pint(1:dim_num,int_num) = p(1:dim_num)
      end if

    end if

  end do

  return
end
function triangle_orientation_2d ( t )

!*****************************************************************************80
!
!! TRIANGLE_ORIENTATION_2D determines the orientation of a triangle in 2D.
!
!  Discussion:
!
!    Three distinct non-colinear points in the plane define a circle.
!    If the points are visited in the order P1, P2, and then
!    P3, this motion defines a clockwise or counter clockwise
!    rotation along the circle.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) T(2,3), the triangle vertices.
!
!    Output, integer TRIANGLE_ORIENTATION_2D, reports if the 
!    three points lie clockwise on the circle that passes through them.  
!    The possible return values are:
!    0, the points are distinct, noncolinear, and lie counter clockwise
!    on their circle.
!    1, the points are distinct, noncolinear, and lie clockwise
!    on their circle.
!    2, the points are distinct and colinear.
!    3, at least two of the points are identical.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) det
  integer triangle_orientation_2d
  real ( kind = rk ) t(dim_num,3)

  if ( all ( t(1:dim_num,1) == t(1:dim_num,2) ) .or. &
       all ( t(1:dim_num,2) == t(1:dim_num,3) ) .or. &
       all ( t(1:dim_num,3) == t(1:dim_num,1) ) ) then
    triangle_orientation_2d = 3
    return
  end if

  det = ( t(1,1) - t(1,3) ) * ( t(2,2) - t(2,3) ) &
      - ( t(1,2) - t(1,3) ) * ( t(2,1) - t(2,3) )

  if ( det == 0.0D+00 ) then
    triangle_orientation_2d = 2
  else if ( det < 0.0D+00 ) then
    triangle_orientation_2d = 1
  else if ( 0.0D+00 < det ) then
    triangle_orientation_2d = 0
  end if

  return
end
subroutine triangle_orthocenter_2d ( t, pc )

!*****************************************************************************80
!
!! TRIANGLE_ORTHOCENTER_2D computes the orthocenter of a triangle in 2D.
!
!  Discussion:
!
!    The orthocenter is defined as the intersection of the three altitudes
!    of a triangle.
!
!    An altitude of a triangle is the line through a vertex of the triangle
!    and perpendicular to the opposite side.
!
!    In geometry, the orthocenter of a triangle is often symbolized by "H".
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    31 July 2009
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
!  Parameters:
!
!    Input, real ( kind = rk ) T(2,3), the triangle vertices.
!
!    Output, real ( kind = rk ) PC(2), the orthocenter of the triangle.
!
!    Output, logical FLAG, is TRUE if the value could not 
!    be computed.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  logical flag
  integer ival
  real ( kind = rk ) p23(dim_num)
  real ( kind = rk ) p31(dim_num)
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) r8_huge
  real ( kind = rk ) t(dim_num,3)
!
!  Determine a point P23 common to the line (P2,P3) and
!  its perpendicular through P1.
!
  call line_exp_perp_2d ( t(1:2,2), t(1:2,3), t(1:2,1), p23, flag )

  if ( flag ) then
    pc(1:2) = r8_huge ( )
    return
  end if
!
!  Determine a point P31 common to the line (P3,P1) and
!  its perpendicular through P2.
!
  call line_exp_perp_2d ( t(1:2,3), t(1:2,1), t(1:2,2), p31, flag )

  if ( flag ) then
    pc(1:2) = r8_huge ( )
    return
  end if
!
!  Determine PC, the intersection of the lines (P1,P23) and (P2,P31).
!
  call lines_exp_int_2d ( t(1:2,1), p23(1:2), t(1:2,2), p31(1:2), ival, pc )

  if ( ival /= 1 ) then
    pc(1:2) = r8_huge ( )
    flag = .true.
    return
  end if

  return
end
subroutine triangle_point_dist_2d ( t, p, dist )

!*****************************************************************************80
!
!! TRIANGLE_POINT_DIST_2D: distance ( triangle, point ) in 2D.
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
!  Parameters:
!
!    Input, real ( kind = rk ) T(2,3), the triangle vertices.
!
!    Input, real ( kind = rk ) P(2), the point to be checked.
!
!    Output, real ( kind = rk ) DIST, the distance from the point to the
!    triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2
  integer, parameter :: side_num = 3

  real ( kind = rk ) dist
  real ( kind = rk ) dist2
  integer i4_wrap
  integer j
  integer jp1
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) t(dim_num,side_num)
!
!  Find the distance to each of the line segments.
!
  dist = huge ( dist )

  do j = 1, side_num

    jp1 = i4_wrap ( j+1, 1, side_num )

    call segment_point_dist_2d ( t(1:dim_num,j), t(1:dim_num,jp1), p, dist2 )

    if ( dist2 < dist ) then
      dist = dist2
    end if

  end do

  return
end
subroutine triangle_point_dist_3d ( t, p, dist )

!*****************************************************************************80
!
!! TRIANGLE_POINT_DIST_3D: distance ( triangle, point ) in 3D.
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
!  Parameters:
!
!    Input, real ( kind = rk ) T(3,3), the triangle vertices.
!
!    Input, real ( kind = rk ) P(3), the point which is to be checked.
!
!    Output, real ( kind = rk ) DIST, the distance from the point to the
!    triangle.  DIST is zero if the point lies exactly on the triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) dist
  real ( kind = rk ) dist2
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) t(dim_num,3)
!
!  Compute the distances from the point to each of the sides.
!
  call segment_point_dist_3d ( t(1:dim_num,1), t(1:dim_num,2), p, dist2 )

  dist = dist2

  call segment_point_dist_3d ( t(1:dim_num,2), t(1:dim_num,3), p, dist2 )

  dist = min ( dist, dist2 )

  call segment_point_dist_3d ( t(1:dim_num,3), t(1:dim_num,1), p, dist2 )

  dist = min ( dist, dist2 )

  return
end
subroutine triangle_point_dist_signed_2d ( t, p, dist_signed )

!*****************************************************************************80
!
!! TRIANGLE_POINT_DIST_SIGNED_2D: signed distance ( triangle, point ) in 2D.
!
!  Discussion:
!
!    If the signed distance is:
!    0, the point is on the boundary of the triangle;
!    negative, the point is in the triangle;
!    positive, the point is outside the triangle.
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
!  Parameters:
!
!    Input, real ( kind = rk ) T(2,3), the triangle vertices.
!    These should be given in counter clockwise order.
!
!    Input, real ( kind = rk ) P(2), the point which is to be checked.
!
!    Output, real ( kind = rk ) DIST_SIGNED, the signed distance from the
!    point to the triangle.  
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) dis12
  real ( kind = rk ) dis23
  real ( kind = rk ) dis31
  real ( kind = rk ) dist_signed
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) t(dim_num,3)
!
!  Compute the signed line distances to the point.
!
  call line_exp_point_dist_signed_2d ( t(1:2,1), t(1:2,2), p, dis12 )

  call line_exp_point_dist_signed_2d ( t(1:2,2), t(1:2,3), p, dis23 )

  call line_exp_point_dist_signed_2d ( t(1:2,3), t(1:2,1), p, dis31 )
!
!  If the point is inside the triangle, all the line distances are negative.
!  The largest (negative) line distance has the smallest magnitude,
!  and is the signed triangle distance.
!
  if ( dis12 <= 0.0D+00 .and. dis23 <= 0.0D+00 .and. dis31 <= 0.0D+00 ) then
    dist_signed = max ( dis12, dis23, dis31 )
!
!  If the point is outside the triangle, then we have to compute
!  the (positive) line segment distances and take the minimum.
!
  else

    call segment_point_dist_2d ( t(1:2,1), t(1:2,2), p, dis12 )
    call segment_point_dist_2d ( t(1:2,2), t(1:2,3), p, dis23 )
    call segment_point_dist_2d ( t(1:2,3), t(1:2,1), p, dis31 )

    dist_signed = min ( dis12, dis23, dis31 )

  end if

  return
end
subroutine triangle_point_near_2d ( t, p, pn, dist )

!*****************************************************************************80
!
!! TRIANGLE_POINT_NEAR_2D computes the nearest point on a triangle in 2D.
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
!  Parameters:
!
!    Input, real ( kind = rk ) T(2,3), the triangle vertices.
!
!    Input, real ( kind = rk ) P(2), the point whose nearest triangle point
!    is to be determined.
!
!    Output, real ( kind = rk ) PN(2), the nearest point to P.
!
!    Output, real ( kind = rk ) DIST, the distance from the point to the
!    triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2
  integer, parameter :: side_num = 3

  real ( kind = rk ) dist
  real ( kind = rk ) dist2
  integer i4_wrap
  integer j
  integer jp1
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) pn(dim_num)
  real ( kind = rk ) pn2(dim_num)
  real ( kind = rk ) t(dim_num,side_num)
  real ( kind = rk ) tval
!
!  Find the distance to each of the line segments that make up the edges
!  of the triangle.
!
  dist = huge ( dist )
  pn(1:dim_num) = 0.0D+00

  do j = 1, side_num

    jp1 = i4_wrap ( j+1, 1, side_num )

    call segment_point_near_2d ( t(1:dim_num,j), t(1:dim_num,jp1), p, &
      pn2, dist2, tval )

    if ( dist2 < dist ) then
      dist = dist2
      pn(1:dim_num) = pn2(1:dim_num)
    end if

  end do

  return
end
subroutine triangle_quality_2d ( t, quality )

!*****************************************************************************80
!
!! TRIANGLE_QUALITY_2D: "quality" of a triangle in 2D.
!
!  Discussion:
!
!    The quality of a triangle is 2.0 times the ratio of the radius of 
!    the inscribed circle divided by that of the circumscribed circle.  
!    An equilateral triangle achieves the maximum possible quality of 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    31 July 2009
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
!  Parameters:
!
!    Input, real ( kind = rk ) T(2,3), the triangle vertices.
!
!    Output, real ( kind = rk ) QUALITY, the quality of the triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) quality
  real ( kind = rk ) t(dim_num,3)
!
!  Compute the length of each side.
!
  a = sqrt ( sum ( ( t(1:dim_num,1) - t(1:dim_num,2) )**2 ) )
  b = sqrt ( sum ( ( t(1:dim_num,2) - t(1:dim_num,3) )**2 ) )
  c = sqrt ( sum ( ( t(1:dim_num,3) - t(1:dim_num,1) )**2 ) )

  if ( a * b * c == 0.0D+00 ) then
    quality = 0.0D+00
  else
    quality = ( - a + b + c ) * ( a - b + c ) * ( a + b - c ) &
      / ( a * b * c )
  end if

  return
end
subroutine triangle_right_lattice_point_num_2d ( a, b, n )

!*****************************************************************************80
!
!! TRIANGLE_RIGHT_LATTICE_POINT_NUM_2D: count lattice points.
!
!  Discussion:
!
!    The triangle is assumed to be a right triangle which, without loss
!    of generality, has the coordinates:
!
!    ( (0,0), (a,0), (0,b) )
!
!    The routine returns the number of integer lattice points that appear
!    inside the triangle or on its edges or vertices.
!
!    The formula for this function occurred to me (JVB) after some thought,
!    on 06 July 2009.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 July 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer A, B, define the vertices.
!
!    Output, integer N, the number of lattice points.
!
  implicit  none

  integer a
  integer b
  integer i4_gcd
  integer n

  n = ( ( a + 1 ) * ( b + 1 ) + i4_gcd ( a, b ) + 1 ) / 2

  return
end
subroutine triangle_sample ( t, n, p )

!*****************************************************************************80
!
!! TRIANGLE_SAMPLE returns random points in a triangle.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) T(2,3), the triangle vertices.
!
!    Input, integer N, the number of points to generate.
!
!    Output, real ( kind = rk ) P(2,N), random points in the triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2
  integer n

  real ( kind = rk ) alpha(n)
  integer dim
  real ( kind = rk ) p(dim_num,n)
  real ( kind = rk ) p12(dim_num,n)
  real ( kind = rk ) p13(dim_num,n)
  real ( kind = rk ) t(dim_num,3)

  call random_number ( harvest = alpha(1:n) )
!
!  Interpret R as a percentage of the triangle's area.
!
!  Imagine a line L, parallel to side 1, so that the area between
!  vertex 1 and line L is R percent of the full triangle's area.
!
!  The line L will intersect sides 2 and 3 at a fraction
!  ALPHA = SQRT ( R ) of the distance from vertex 1 to vertices 2 and 3.
!
  alpha(1:n) = sqrt ( alpha(1:n) )
!
!  Determine the coordinates of the points on sides 2 and 3 intersected
!  by line L.
!
  do dim = 1, dim_num

    p12(dim,1:n) = ( 1.0D+00 - alpha(1:n) ) * t(dim,1) &
                             + alpha(1:n)   * t(dim,2)

    p13(dim,1:n) = ( 1.0D+00 - alpha(1:n) ) * t(dim,1) &
                             + alpha(1:n)   * t(dim,3)

  end do
!
!  Now choose, uniformly at random, a point on the line L.
!
  call random_number ( harvest = alpha(1:n) )

  do dim = 1, dim_num

    p(dim,1:n) = ( 1.0D+00 - alpha(1:n) ) * p12(dim,1:n) &
                           + alpha(1:n)   * p13(dim,1:n)

  end do

  return
end
subroutine triangle01_lattice_point_num_2d ( s, n )

!*****************************************************************************80
!
!! TRIANGLE01_LATTICE_POINT_NUM_2D: count lattice points.
!
!  Discussion:
!
!    The triangle is assumed to be the unit triangle:
!
!    ( (0,0), (1,0), (0,1) )
!
!    or a copy of this triangle scaled by an integer S:
!
!    ( (0,0), (S,0), (0,S) ).
!
!    The routine returns the number of integer lattice points that appear
!    inside the triangle or on its edges or vertices.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 July 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Matthias Beck, Sinai Robins,
!    Computing the Continuous Discretely,
!    Springer, 2006,
!    ISBN13: 978-0387291390,
!    LC: QA640.7.B43.
!
!  Parameters:
!
!    Input, integer S, the scale factor.
!
!    Output, integer N, the number of lattice points.
!
  implicit  none

  integer n
  integer s

  n = ( ( s + 2 ) * ( s + 1 ) ) / 2

  return
end
subroutine triangle_xsi_to_xy_2d ( t, xsi, p )

!*****************************************************************************80
!
!! TRIANGLE_XSI_TO_XY_2D converts from barycentric to XY coordinates in 2D.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) T(2,3), the triangle vertices.
!
!    Input, real ( kind = rk ) XSI(3), the barycentric coordinates of a point.
!    XSI(1) + XSI(2) + XSI(3) should equal 1, but this is not checked.
!
!    Output, real ( kind = rk ) P(2), the XY coordinates of the point.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) t(dim_num,3)
  real ( kind = rk ) xsi(dim_num+1)

  p(1:dim_num) = matmul ( t(1:dim_num,1:3), xsi(1:dim_num+1) )

  return
end
subroutine triangle_xy_to_xsi_2d ( t, p, xsi )

!*****************************************************************************80
!
!! TRIANGLE_XY_TO_XSI_2D converts from XY to barycentric in 2D.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) T(2,3), the triangle vertices.
!
!    Input, real ( kind = rk ) P(2), the XY coordinates of a point.
!
!    Output, real ( kind = rk ) XSI(3), the barycentric coordinates of the point.
!    XSI1 + XSI2 + XSI3 should equal 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) det
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) t(dim_num,3)
  real ( kind = rk ) xsi(3)

  det = ( t(1,1) - t(1,3) ) * ( t(2,2) - t(2,3) ) &
      - ( t(1,2) - t(1,3) ) * ( t(2,1) - t(2,3) )

  xsi(1) = (   ( t(2,2) - t(2,3) ) * ( p(1) - t(1,3) ) &
             - ( t(1,2) - t(1,3) ) * ( p(2) - t(2,3) ) ) / det

  xsi(2) = ( - ( t(2,1) - t(2,3) ) * ( p(1) - t(1,3) ) &
             + ( t(1,1) - t(1,3) ) * ( p(2) - t(2,3) ) ) / det

  xsi(3) = 1.0D+00 - xsi(1) - xsi(2)

  return
end
subroutine truncated_octahedron_shape_3d ( point_num, face_num, &
  face_order_max, point_coord, face_order, face_point )

!*****************************************************************************80
!
!! TRUNCATED_OCTAHEDRON_SHAPE_3D describes a truncated octahedron in 3D.
!
!  Discussion:
!
!    The shape is a truncated octahedron.  There are 8 hexagons and 6
!    squares.
!
!    The truncated octahedron is an interesting shape because it
!    is "space filling".  In other words, all of 3D space can be
!    filled by a regular lattice of these shapes.
!
!    Call TRUNCATED_OCTAHEDRON_SIZE_3D to get the values of POINT_NUM, 
!    FACE_NUM, and FACE_ORDER_MAX, so you can allocate space for the arrays.
!
!    For each face, the face list must be of length FACE_ORDER_MAX.
!    In cases where a face is of lower than maximum order (the
!    squares, in this case), the extra entries are listed as "-1".
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer POINT_NUM, the number of points (24).
!
!    Input, integer FACE_NUM, the number of faces (14).
!
!    Input, integer FACE_ORDER_MAX, the maximum order of any 
!    face (6).
!
!    Output, real ( kind = rk ) POINT_COORD(3,POINT_NUM), the vertices.
!
!    Output, integer FACE_ORDER(FACE_NUM), the number of 
!    vertices per face.
!
!    Output, integer FACE_POINT(FACE_ORDER_MAX,FACE_NUM); 
!    FACE_POINT(I,J) contains the index of the I-th point in the J-th face.  
!    The points are listed in the counter clockwise direction defined
!    by the outward normal at the face.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer face_num
  integer face_order_max
  integer, parameter :: dim_num = 3
  integer point_num

  integer face_order(face_num)
  integer face_point(face_order_max,face_num)
  real ( kind = rk ) point_coord(dim_num,point_num)
!
!  Set the point coordinates.
!
  point_coord(1:dim_num,1:point_num) = reshape ( (/ &
    -1.5D+00, -0.5D+00,  0.0D+00,        &
    -1.5D+00,  0.5D+00,  0.0D+00,        &
    -1.0D+00, -1.0D+00, -0.70710677D+00, &
    -1.0D+00, -1.0D+00,  0.70710677D+00, &
    -1.0D+00,  1.0D+00, -0.70710677D+00, &
    -1.0D+00,  1.0D+00,  0.70710677D+00, &
    -0.5D+00, -1.5D+00,  0.0D+00,        &
    -0.5D+00, -0.5D+00, -1.4142135D+00,  &
    -0.5D+00, -0.5D+00,  1.4142135D+00,  &
    -0.5D+00,  0.5D+00, -1.4142135D+00,  &
    -0.5D+00,  0.5D+00,  1.4142135D+00,  &
    -0.5D+00,  1.5D+00,  0.0D+00,        &
     0.5D+00, -1.5D+00,  0.0D+00,        &
     0.5D+00, -0.5D+00, -1.4142135D+00,  &
     0.5D+00, -0.5D+00,  1.4142135D+00,  &
     0.5D+00,  0.5D+00, -1.4142135D+00,  &
     0.5D+00,  0.5D+00,  1.4142135D+00,  &
     0.5D+00,  1.5D+00,  0.0D+00,        &
     1.0D+00, -1.0D+00, -0.70710677D+00, &
     1.0D+00, -1.0D+00,  0.70710677D+00, &
     1.0D+00,  1.0D+00, -0.70710677D+00, &
     1.0D+00,  1.0D+00,  0.70710677D+00, &
     1.5D+00, -0.5D+00,  0.0D+00,        &
     1.5D+00,  0.5D+00,  0.0D+00 /), (/ dim_num, point_num /) )
!
!  Set the face orders.
!
  face_order(1:face_num) = (/ &
    4, 4, 4, 4, 4, 4, 6, 6, 6, 6, &
    6, 6, 6, 6 /)
!
!  Set faces.
!
  face_point(1:face_order_max,1:face_num) = reshape ( (/ &
    17, 11,  9, 15, -1, -1, &
    14,  8, 10, 16, -1, -1, &
    22, 24, 21, 18, -1, -1, &
    12,  5,  2,  6, -1, -1, &
    13, 19, 23, 20, -1, -1, &
     4,  1,  3,  7, -1, -1, &
    19, 13,  7,  3,  8, 14, &
    15,  9,  4,  7, 13, 20, &
    16, 10,  5, 12, 18, 21, &
    22, 18, 12,  6, 11, 17, &
    20, 23, 24, 22, 17, 15, &
    14, 16, 21, 24, 23, 19, &
     9, 11,  6,  2,  1,  4, &
     3,  1,  2,  5, 10,  8 /), (/ face_order_max, face_num /) )

  return
end
subroutine truncated_octahedron_size_3d ( point_num, edge_num, face_num, &
  face_order_max )

!*****************************************************************************80
!
!! TRUNCATED_OCTAHEDRON_SIZE_3D gives "sizes" for a truncated octahedron in 3D.
!
!  Discussion:
!
!    The truncated octahedron is "space-filling".
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer POINT_NUM, the number of points.
!
!    Output, integer EDGE_NUM, the number of edges.
!
!    Output, integer FACE_NUM, the number of faces.
!
!    Output, integer FACE_ORDER_MAX, the maximum order of any face.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer edge_num
  integer face_num
  integer face_order_max
  integer point_num

  point_num = 24
  edge_num = 36
  face_num = 14
  face_order_max = 6

  return
end
subroutine tube_2d ( dist, n, p, p1, p2 )

!*****************************************************************************80
!
!! TUBE_2D constructs a "tube" of given width around a path in 2D.
!
!  Discussion:
!
!    The routine is given a sequence of N points, and a distance DIST.
!
!    It returns the coordinates of the corners of the top and bottom
!    of a tube of width 2*DIST, which envelopes the line connecting
!    the points.
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
!  Parameters:
!
!    Input, real ( kind = rk ) DIST, the radius of the tube.
!
!    Input, integer N, the number of points defining the line.
!    N must be at least 2.
!
!    Input, real ( kind = rk ) P(2,N), the points which comprise the broken
!    line which is to be surrounded by the tube.  Points should
!    not be immediately repeated, that is, it should never be
!    the case that
!      P(1,I) = P(1,I+1) and P(2,I) = P(2,I+1).
!
!    Output, real ( kind = rk ) P1(2,N), P2(2,N), the points P1 form
!    one side of the tube, and P2 the other.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: dim_num = 2

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) dis1
  real ( kind = rk ) dis2
  real ( kind = rk ) dist
  integer i
  integer i4_wrap
  integer im1
  integer ip1
  real ( kind = rk ) p(dim_num,n)
  real ( kind = rk ) p1(dim_num,n)
  real ( kind = rk ) p2(dim_num,n)
  real ( kind = rk ) temp
!
!  Check that N is at least 3.
!
  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TUBE_2D - Fatal error!'
    write ( *, '(a)' ) '  N must be at least 3'
    write ( *, '(a,i8)' ) '  but your input value was N = ', n
    stop 1
  end if
!
!  Check that consecutive points are distinct.
!
  do i = 1, n - 1
    if ( all ( p(1:2,i) == p(1:2,i+1) ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TUBE_2D - Fatal error!'
      write ( *, '(a,i8)' ) '  P(1:2,I) = P(1:2,I+1) for I = ', i
      write ( *, '(a,2g14.6)' ) '  P(1:2,I) = ', p(1:2,i)
      stop 1
    end if
  end do

  do i = 1, n

    im1 = i4_wrap ( i-1, 1, n )
    ip1 = i4_wrap ( i+1, 1, n )

    call angle_box_2d ( dist, p(1:2,im1), p(1:2,i), &
      p(1:2,ip1), p1(1:2,i), p2(1:2,i) )
!
!  On the first and last steps, translate the corner points DIST units
!  along the line, to make an extra buffer.
!
    if ( i == 1 ) then

      temp = sqrt ( ( p(1,2) - p(1,1) )**2 + ( p(2,2) - p(2,1) )**2 )
      p1(1:2,1) = p1(1:2,1) - dist * ( p(1:2,2) - p(1:2,1) ) / temp
      p2(1:2,1) = p2(1:2,1) - dist * ( p(1:2,2) - p(1:2,1) ) / temp

    else if ( i == n ) then

      temp = sqrt ( ( p(1,n) - p(1,n-1) )**2 + ( p(2,n) - p(2,n-1) )**2 )
      p1(1:2,n) = p1(1:2,n) + dist * ( p(1:2,n) - p(1:2,n-1) ) / temp
      p2(1:2,n) = p2(1:2,n) + dist * ( p(1:2,n) - p(1:2,n-1) ) / temp

    end if
!
!  The new points P1 and P2 may need to be swapped.
!
!  Compute the signed distance from the points to the line.
!
    if ( 1 < i ) then

      a = p(2,i-1) - p(2,i)
      b = p(1,i) - p(1,i-1)
      c = p(1,i-1) * p(2,i) - p(1,i) * p(2,i-1)

      dis1 = ( a * p1(1,i-1) + b * p1(2,i-1) + c ) / sqrt ( a * a + b * b )

      dis2 = ( a * p1(1,i) + b * p1(2,i) + c ) / sqrt ( a * a + b * b )

      if ( sign ( 1.0D+00, dis1 ) /= sign ( 1.0D+00, dis2 ) ) then

        call r8_swap ( p1(1,i), p2(1,i) )
        call r8_swap ( p1(2,i), p2(2,i) )

      end if

    end if

  end do

  return
end
subroutine vector_directions_nd ( dim_num, v, angle )

!*****************************************************************************80
!
!! VECTOR_DIRECTIONS_ND returns the direction angles of a vector in ND.
!
!  Discussion:
!
!    Let V be the vector, and let E(I) be the I-th unit coordinate axis vector.
!    The I-th direction angle is the angle between V and E(I), which is
!    the angle whose cosine is equal to the direction cosine:
!
!      Direction_Cosine(I) = V dot E(I) / |V|.
!
!    If V is the null or zero vector, then the direction cosines and
!    direction angles are undefined, and this routine simply returns
!    zeroes.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = rk ) V(DIM_NUM), the vector.
!
!    Output, real ( kind = rk ) ANGLE(DIM_NUM), the direction angles, in radians,
!    that the vector V makes with the coordinate axes.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num

  real ( kind = rk ) angle(dim_num)
  real ( kind = rk ) v(dim_num)
  real ( kind = rk ) vnorm
!
!  Get the norm of the vector.
!
  vnorm = sqrt ( sum ( v(1:dim_num)**2 ) )

  if ( vnorm == 0.0D+00 ) then
    angle(1:dim_num) = 0.0D+00
    return
  end if

  angle(1:dim_num) = acos ( v(1:dim_num) / vnorm )

  return
end
subroutine vector_rotate_2d ( v, angle, w )

!*****************************************************************************80
!
!! VECTOR_ROTATE_2D rotates a vector around the origin in 2D.
!
!  Discussion:
!
!    To see why this formula is so, consider that the original point
!    has the form ( R cos Theta, R sin Theta ), and the rotated point
!    has the form ( R cos ( Theta + Angle ), R sin ( Theta + Angle ) ).
!    Now use the addition formulas for cosine and sine to relate
!    the new point to the old one:
!
!      ( W1 ) = ( cos Angle  - sin Angle ) * ( V1 )
!      ( W2 )   ( sin Angle    cos Angle )   ( V2 )
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) V(2), the components of the vector to be
!    rotated.
!
!    Input, real ( kind = rk ) ANGLE, the angle, in radians, of the rotation
!    to be carried out.  A positive angle rotates the vector in the
!    counter clockwise direction.
!
!    Output, real ( kind = rk ) W(2), the rotated vector.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) angle
  real ( kind = rk ) v(dim_num)
  real ( kind = rk ) w(dim_num)

  w(1) = cos ( angle ) * v(1) - sin ( angle ) * v(2)
  w(2) = sin ( angle ) * v(1) + cos ( angle ) * v(2)

  return
end
subroutine vector_rotate_3d ( v1, axis, angle, v2 )

!*****************************************************************************80
!
!! VECTOR_ROTATE_3D rotates a vector around an axis vector in 3D.
!
!  Discussion:
!
!    Thanks to Cody Farnell for correcting some errors in a previous
!    version of this routine!
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    28 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) V1(3), the vector to be rotated.
!
!    Input, real ( kind = rk ) AXIS(3), the vector about which the
!    rotation is to be carried out.
!
!    Input, real ( kind = rk ) ANGLE, the angle, in radians, of the rotation
!    to be carried out.
!
!    Output, real ( kind = rk ) V2(3), the rotated vector.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) angle
  real ( kind = rk ) axis(3)
  real ( kind = rk ) dot
  real ( kind = rk ) norm
  real ( kind = rk ) norm_vn
  real ( kind = rk ) normal2(3)
  real ( kind = rk ) r8vec_norm
  real ( kind = rk ) v1(3)
  real ( kind = rk ) v2(3)
  real ( kind = rk ) vn(3)
  real ( kind = rk ) vp(3)
  real ( kind = rk ) vr(3)
!
!  Compute the length of the rotation axis.
!
  norm = r8vec_norm ( 3, axis )

  if ( norm == 0.0D+00 ) then
    v2(1:3) = v1(1:3)
    return
  end if
!
!  Compute the dot product of the vector and the (unit) rotation axis.
!
  dot = dot_product ( v1(1:3), axis(1:3) ) / norm
!
!  Compute the parallel component of the vector.
!
  vp(1:3) = dot * axis(1:3) / norm
!
!  Compute the normal component of the vector.
!
  vn(1:3) = v1(1:3) - vp(1:3)

  norm_vn = r8vec_norm ( 3, vn )

  if ( norm_vn == 0.0D+00 ) then
    v2(1:3) = vp(1:3)
    return
  end if

  vn(1:3) = vn(1:3) / norm_vn
!
!  Compute a second vector, lying in the plane, perpendicular
!  to V1 and VN, and forming a right-handed system.
!
  normal2(1) = axis(2) * vn(3) - axis(3) * vn(2)
  normal2(2) = axis(3) * vn(1) - axis(1) * vn(3)
  normal2(3) = axis(1) * vn(2) - axis(2) * vn(1)

  norm = r8vec_norm ( 3, normal2 )
  if ( norm /= 0.0D+00 ) then
    normal2(1:3) = normal2(1:3) / norm
  end if
!
!  Rotate the normal component by the angle.
!
  vr(1:3) = norm_vn * ( cos ( angle ) * vn(1:3) + sin ( angle ) * normal2(1:3) )
!
!  The rotated vector is the parallel component plus the rotated component.
!
  v2(1:3) = vp(1:3) + vr(1:3)

  return
end
subroutine vector_rotate_base_2d ( p1, pb, angle, p2 )

!*****************************************************************************80
!
!! VECTOR_ROTATE_BASE_2D rotates a vector around a base point in 2D.
!
!  Discussion:
!
!    The original vector is assumed to be ( X1-XB, Y1-YB ), and the
!    rotated vector is ( X2-XB, Y2-YB ).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    29 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) P1(2), the endpoint of the original vector.
!
!    Input, real ( kind = rk ) PB(2), the location of the base point.
!
!    Input, real ( kind = rk ) ANGLE, the angle, in radians, of the rotation
!    to be carried out.  A positive angle rotates the vector in the
!    counter clockwise direction.
!
!    Output, real ( kind = rk ) P2(2), the endpoint of the rotated vector.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) angle
  real ( kind = rk ) p1(2)
  real ( kind = rk ) p2(2)
  real ( kind = rk ) pb(2)

  p2(1) = pb(1) + cos ( angle ) * ( p1(1) - pb(1) ) &
                - sin ( angle ) * ( p1(2) - pb(2) )

  p2(2) = pb(2) + sin ( angle ) * ( p1(1) - pb(1) ) &
                + cos ( angle ) * ( p1(2) - pb(2) )

  return
end
subroutine vector_separation_nd ( dim_num, v1, v2, theta )

!*****************************************************************************80
!
!! VECTOR_SEPARATION_ND finds the angular separation between vectors in ND.
!
!  Discussion:
!
!    Any two vectors lie in a plane, and are separated by a plane angle.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = rk ) V1(DIM_NUM), V2(DIM_NUM), the two vectors.
!
!    Output, real ( kind = rk ) THETA, the angle between the two vectors.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num

  real ( kind = rk ) cos_theta
  real ( kind = rk ) r8_acos
  real ( kind = rk ) theta
  real ( kind = rk ) v1(dim_num)
  real ( kind = rk ) v1_norm
  real ( kind = rk ) v2(dim_num)
  real ( kind = rk ) v2_norm

  v1_norm = sqrt ( sum ( v1(1:dim_num)**2 ) )

  v2_norm = sqrt ( sum ( v2(1:dim_num)**2 ) )

  cos_theta = dot_product ( v1(1:dim_num), v2(1:dim_num) ) &
    / ( v1_norm * v2_norm )

  theta = r8_acos ( cos_theta )

  return
end
subroutine vector_unit_nd ( dim_num, v )

!*****************************************************************************80
!
!! VECTOR_UNIT_ND normalizes a vector in ND.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input/output, real ( kind = rk ) V(DIM_NUM), the vector to be normalized.
!    On output, V should have unit Euclidean norm.  However, if the input vector
!    has zero Euclidean norm, it is not altered.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num

  real ( kind = rk ) norm
  real ( kind = rk ) v(dim_num)

  norm = sqrt ( sum ( v(1:dim_num)**2 ) )

  if ( norm /= 0.0D+00 ) then
    v(1:dim_num) = v(1:dim_num) / norm
  end if

  return
end
function voxels_dist_l1_nd ( dim_num, v1, v2 )

!*****************************************************************************80
!
!! VOXELS_DIST_L1_ND computes the L1 distance between voxels in ND.
!
!  Discussion:
!
!    A voxel is generally a point in 3D space with integer coordinates.
!    There's no reason to stick with 3D, so this routine will handle
!    any dimension.
!
!    We can imagine that, in traveling from V1 to V2, we are allowed to 
!    increment or decrement just one coordinate at a time.  The minimum number 
!    of such changes required is the L1 distance. 
!
!    More formally,
!
!      DIST_L1 ( V1, V2 ) = sum ( 1 <= I <= N ) | V1(I) - V2(I) |
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
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer V1(DIM_NUM), the voxel that begins the line.
!
!    Input, integer V2(DIM_NUM), the voxel that ends the line.
!
!    Output, integer VOXELS_DIST_L1_ND, the L1 distance 
!    between the voxels.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num

  integer v1(dim_num)
  integer v2(dim_num)
  integer voxels_dist_l1_nd

  voxels_dist_l1_nd = sum ( abs ( v1(1:dim_num) - v2(1:dim_num) ) )

  return
end
subroutine voxels_line_3d ( v1, v2, n, v )

!*****************************************************************************80
!
!! VOXELS_LINE_3D computes voxels along a line in 3D.
!
!  Discussion:
!
!    The line itself is defined by two voxels.  The line will begin
!    at the first voxel, and move towards the second.  If the value of
!    N is equal to the L1 distance between the two voxels, then the
!    line will "almost" reach the second voxel.  Depending on the
!    direction, 1, 2 or 3 more steps may be needed.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 March 2005
!
!  Author:
!
!    Daniel Cohen
!
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Daniel Cohen,
!    Voxel Traversal along a 3D Line,
!    in Graphics Gems IV,
!    edited by Paul Heckbert,
!    AP Professional, 1994,
!    T385.G6974.
!
!  Parameters:
!
!    Input, integer V1(3), the voxel that begins the line.
!
!    Input, integer V2(3), the voxel that ends the line.
!
!    Input, integer N, the number of voxels to compute.
!
!    Output, integer V(3,N), a sequence of voxels, whose
!    first value is V1 and which proceeds towards V2.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: dim_num = 3

  integer a(3)
  integer exy
  integer exz
  integer ezy
  integer i
  integer s(3)
  integer v(3,n)
  integer v1(3)
  integer v2(3)

  if ( n <= 0 ) then
    return
  end if
!
!  Determine the number of voxels on the line.
!
  s(1:dim_num) = sign ( 1, v2(1:dim_num) - v1(1:dim_num) )
  a(1:dim_num) = abs ( v2(1:dim_num) - v1(1:dim_num) )

  exy = a(2) - a(1)
  exz = a(3) - a(1)
  ezy = a(2) - a(3)
!
!  We start at the starting point.
!
  v(1:dim_num,1) = v1(1:dim_num)

  do i = 2, n

    v(1:dim_num,i) = v(1:dim_num,i-1)

    if ( exy < 0 ) then

      if ( exz < 0 ) then
        v(1,i) = v(1,i) + s(1)
        exy = exy + 2 * a(2)
        exz = exz + 2 * a(3)
      else
        v(3,i) = v(3,i) + s(3)
        exz = exz - 2 * a(1)
        ezy = ezy + 2 * a(2)
      end if

    else if ( ezy < 0 ) then

      v(3,i) = v(3,i) + s(3)
      exz = exz - 2 * a(1)
      ezy = ezy + 2 * a(2)

    else

      v(2,i) = v(2,i) + s(2)
      exy = exy - 2 * a(1)
      ezy = ezy - 2 * a(3)

    end if

  end do

  return
end
subroutine voxels_region_3d ( list_max, nx, ny, nz, ishow, list_num, list, &
  region_num )

!*****************************************************************************80
!
!! VOXELS_REGION_3D arranges contiguous voxels into regions in 3D.
!
!  Discussion:
!
!    On input, the ISHOW array contains zero and nonzero values.  The nonzero
!    values are taken to be active voxels.  On output, the zero voxels remain
!    zero, and all the active voxels have been assigned a value which now
!    indicates membership in a region, or group of contiguous voxels.
!
!    On output, the array LIST contains information about the regions.
!    The last used element of LIST is LIST_NUM.
!
!    The number of elements in region REGION_NUM is NELEM = LIST(LIST_NUM).  
!    The (I,J,K) indices of the last element in this region are in
!    LIST(LIST_NUM-3) through LIST(LIST_NUM-1), and the first element is
!    listed in LIST(LIST_NUM-3*NELEM), LIST(LIST_NUM-3*NELEM+1),
!    LIST(LIST_NUM-3*NELEM+2).
!
!    The number of elements in REGION_NUM-1 is listed in
!    LIST(LIST_NUM-3*NELEM-1), 
!    and the (I,J,K) indices of the these elements are listed there.
!
!    Thanks to Emre Evren for pointing out a hard-to-spot error involving
!    a DO loop that mistakenly read "DO 1 = 1, N".
!
!  Picture:
!
!    Input:
!
!      0  2  0  0 17  0  3
!      0  0  3  0  1  0  4
!      1  0  4  8  8  0  7
!      3  0  6 45  0  0  0
!      3 17  0  5  9  2  5
!
!    Output:
!
!      0  1  0  0  2  0  3
!      0  0  2  0  2  0  3
!      4  0  2  2  2  0  3
!      4  0  2  2  0  0  0
!      4  4  0  2  2  2  2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LIST_MAX, the maximum length of the array 
!    used to list the elements of the regions.
!
!    Input, integer NX, NY, NZ, the number of voxels in the X, Y 
!    and Z directions.
!
!    Input/output, integer ISHOW(NX,NY,NZ).  On input, the only 
!    significance to the entries is whether they are zero or nonzero.  On 
!    output, the nonzero entries have now been revalued so that contiguous 
!    entries have the same value, indicating a grouping into a region.
!
!    Output, integer LIST_NUM, the number of entries of LIST that 
!    were used.  However, if LIST_MAX < LIST_NUM, then there was not enough 
!    space in LIST to store the data properly, and LIST should not be used,
!    although the data in ISHOW should be correct.
!
!    Output, integer LIST(LIST_MAX), contains, in stack form, a 
!    list of the indices of the elements in each region.
!
!    Output, integer REGION_NUM, the number of regions discovered.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: maxstack = 100

  integer list_max
  integer nx
  integer ny
  integer nz

  integer i
  integer i2
  integer ibase
  integer ihi
  integer ilo
  integer ishow(nx,ny,nz)
  integer j
  integer j2
  integer jbase
  integer jhi
  integer jlo
  integer k
  integer k2
  integer kbase
  integer khi
  integer klo
  integer list(list_max)
  integer list_num
  integer nabes
  integer ncan
  integer nelements
  integer nstack
  integer region_num
  integer stack(maxstack)
!
!  Reset all nonzero entries of ISHOW to -1.
!
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx

        if ( ishow(i,j,k) /= 0 ) then
          ishow(i,j,k) = -1
        end if

      end do
    end do
  end do
!
!  Start the number of items in the region list at 0.
!
  list_num = 0
!
!  Start the number of regions at 0.
!
  region_num = 0
!
!  The stack begins empty.
!
  nstack = 0
!
!  Search for an unused "ON" voxel from which we can "grow" a new region.
!
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
!
!  We found a voxel that is "ON", and does not belong to any region.
!
        if ( ishow(i,j,k) == -1 ) then
!
!  Increase the number of regions.
!
          region_num = region_num + 1
!
!  Add this voxel to the region.
!
          ishow(i,j,k) = region_num
!
!  Add this voxel to the stack.
!
          if ( maxstack < nstack + 4 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'VOXELS_REGION - Fatal error!'
            write ( *, '(a)' ) '  The internal stack overflowed.'
            write ( *, '(a)' ) '  The algorithm has failed.'
            stop 1
          end if

          stack(nstack+1) = i
          stack(nstack+2) = j
          stack(nstack+3) = k
          stack(nstack+4) = 1

          nstack = nstack + 4
!
!  Add this voxel to the description of the region.
!
          nelements = 1

          if ( list_num + 3 <= list_max ) then
            list(list_num+1) = i
            list(list_num+2) = j
            list(list_num+3) = k
          end if

          list_num = list_num + 3

          do
!
!  Find all neighbors of BASE that are "ON" but unused.
!  Mark them as belonging to this region, and stack their indices.
!
            ibase = stack(nstack-3)
            jbase = stack(nstack-2)
            kbase = stack(nstack-1)

            ilo = max ( ibase-1, 1 )
            ihi = min ( ibase+1, nx )
            jlo = max ( jbase-1, 1 )
            jhi = min ( jbase+1, ny )
            klo = max ( kbase-1, 1 )
            khi = min ( kbase+1, nz )

            nabes = 0

            do i2 = ilo, ihi
              do j2 = jlo, jhi
                do k2 = klo, khi
!
!  We found a neighbor to our current search point, which is "ON" and unused.
!
                  if ( ishow(i2,j2,k2) == -1 ) then
!
!  Increase the number of neighbors.
!
                    nabes = nabes + 1
!
!  Mark the neighbor as belonging to the region.
!
                    ishow(i2,j2,k2) = region_num
!
!  Add the neighbor to the stack.
!
                    if ( maxstack < nstack + 3 ) then
                      write ( *, '(a)' ) ' '
                      write ( *, '(a)' ) 'VOXELS_REGION - Fatal error!'
                      write ( *, '(a)' ) '  The internal stack overflowed.'
                      write ( *, '(a)' ) '  The algorithm has failed.'
                      stop 1
                    end if

                    stack(nstack+1) = i2
                    stack(nstack+2) = j2
                    stack(nstack+3) = k2

                    nstack = nstack + 3
!
!  Add the neighbor to the description of the region.
!
                    nelements = nelements + 1

                    if ( list_num + 3 <= list_max ) then
                      list(list_num+1) = i2
                      list(list_num+2) = j2
                      list(list_num+3) = k2
                    end if

                    list_num = list_num + 3

                  end if

                end do
              end do
            end do
!
!  If any new neighbors were found, take the last one as the basis
!  for a deeper search.
!
            if ( 0 < nabes ) then

              if ( maxstack < nstack + 1 ) then
                write ( *, '(a)' ) ' '
                write ( *, '(a)' ) 'VOXELS_REGION - Fatal error!'
                write ( *, '(a)' ) '  The internal stack overflowed.'
                write ( *, '(a)' ) '  The algorithm has failed.'
                stop 1
              end if

              stack(nstack+1) = nabes
              nstack = nstack + 1
              cycle

            end if
!
!  If the current search point had no new neighbors, drop it from the stack.
!
            ncan = stack(nstack) - 1
            nstack = nstack - 3
            stack(nstack) = ncan
!
!  If there are still any unused candidates at this level, take the
!  last one as the basis for a deeper search.
!
            if ( 0 < stack(nstack) ) then
              cycle
            end if
!
!  If there are no more unused candidates at this level, then we need
!  to back up a level in the stack.  If there are any candidates at
!  that earlier level, then we can still do more searching.
!
            nstack = nstack - 1

            if ( nstack <= 0 ) then
              exit
            end if

          end do
!
!  If we have exhausted the stack, we have completed this region.
!  Tag the number of elements to the end of the region description list.
!
          list_num = list_num + 1
          if ( list_num <= list_max ) then
            list(list_num) = nelements
          end if

        end if

      end do
    end do
  end do
!
!  Print some warnings.
!
  if ( list_max < list_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VOXELS_REGION - Warning!'
    write ( *, '(a)' ) '  LIST_MAX was too small to list the regions.'
    write ( *, '(a)' ) '  Do not try to use the LIST array!'
    write ( *, '(a)' ) '  The ISHOW data is OK, however.'
  end if

  return
end
subroutine voxels_step_3d ( v1, v2, inc, jnc, knc, v3 )

!*****************************************************************************80
!
!! VOXELS_STEP_3D computes voxels along a line from a given point in 3D.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer V1(3), the coordinates of the base voxel from
!    which the line begins.
!
!    Input, integer V2(3), the coordinates of the current voxel
!    on the line.  For the first call, these might be equal to V1.
!
!    Input, integer INC, JNC, KNC, the increments to the voxels.
!    These values define the direction along which the line proceeds.
!    However, the voxels on the line will typically be incremented
!    by a fractional value of the vector (INC,JNC,KNC), and the
!    result is essentially rounded.
!
!    Output, integer V3(3), the coordinates of the next voxel along
!    the line.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) alpha
  real ( kind = rk ) alphai
  real ( kind = rk ) alphaj
  real ( kind = rk ) alphak
  integer inc
  integer jnc
  integer knc
  integer v1(3)
  integer v2(3)
  integer v3(3)

  v3(1:dim_num) = v2(1:dim_num)
!
!  Assuming for the moment that (I,J,K) can take on real values,
!  points on the line have the form:
!
!    I = V1(2) + alpha * inc
!    J = V1(2) + alpha * jnc
!    K = V1(3) + alpha * knc
!
  if ( inc == 0 .and. jnc == 0 .and. knc == 0 ) then
    return
  end if

  alpha = 0.0D+00
!
!  Compute the smallest ALPHA that will change one of V2(1:3) by +-0.5.
!
  if ( 0 < inc ) then
    alphai = ( real ( v2(1) - v1(1), kind = rk ) + 0.5D+00 ) &
             / real ( inc, kind = rk )
  else if ( inc < 0 ) then
    alphai = ( real ( v2(1) - v1(1), kind = rk ) - 0.5D+00 ) &
             / real ( inc, kind = rk )
  else
    alphai = huge ( alphai )
  end if

  if ( 0 < jnc ) then
    alphaj = ( real ( v2(2) - v1(2), kind = rk ) + 0.5D+00 ) &
             / real ( jnc, kind = rk )
  else if ( jnc < 0 ) then
    alphaj = ( real ( v2(2) - v1(2), kind = rk ) - 0.5D+00 ) &
             / real ( jnc, kind = rk )
  else
    alphaj = huge ( alphaj )
  end if

  if ( 0 < knc ) then
    alphak = ( real ( v2(3) - v1(3), kind = rk ) + 0.5D+00 ) &
             / real ( knc, kind = rk )
  else if ( knc < 0 ) then
    alphak = ( real ( v2(3) - v1(3), kind = rk ) - 0.5D+00 ) &
             / real ( knc, kind = rk )
  else
    alphaj = huge ( alphaj )
  end if
!
!  The ALPHA of smallest positive magnitude represents the closest next voxel.
!
  alpha = huge ( alpha )

  if ( 0.0D+00 < alphai ) then
    alpha = min ( alpha, alphai )
  end if

  if ( 0.0D+00 < alphaj ) then
    alpha = min ( alpha, alphaj )
  end if

  if ( 0.0D+00 < alphak ) then
    alpha = min ( alpha, alphak )
  end if
!
!  Move to the new voxel.  Whichever index just made the half
!  step must be forced to take a whole step.
!
  if ( alpha == alphai ) then
    v3(1) = v2(1) + sign ( 1, inc )
    v3(2) = v1(2) + nint ( alpha * real ( jnc, kind = rk ) )
    v3(3) = v1(3) + nint ( alpha * real ( knc, kind = rk ) )
  else if ( alpha == alphaj ) then
    v3(1) = v1(1) + nint ( alpha * real ( inc, kind = rk ) )
    v3(2) = v2(2) + sign ( 1, jnc )
    v3(3) = v1(3) + nint ( alpha * real ( knc, kind = rk ) )
  else if ( alpha == alphak ) then
    v3(1) = v1(1) + nint ( alpha * real ( inc, kind = rk ) )
    v3(2) = v1(2) + nint ( alpha * real ( jnc, kind = rk ) )
    v3(3) = v2(3) + sign ( 1, knc )
  end if

  return
end
function wedge01_volume ( )

!*****************************************************************************80
!
!! WEDGE01_VOLUME: volume of the unit wedge in 3D.
!
!  Discussion:
!
!    The integration region is defined as:
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
!    15 January 2018
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) WEDGE01_VOLUME, the volume.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) wedge01_volume

  wedge01_volume = 1.0D+00

  return
end
subroutine xy_to_polar ( xy, r, t )

!*****************************************************************************80
!
!! XY_TO_POLAR converts XY coordinates to polar coordinates.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) XY(2), the Cartesian coordinates.
!
!    Output, real ( kind = rk ) R, T, the radius and angle (in radians).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r
  real ( kind = rk ) r8_atan
  real ( kind = rk ) t
  real ( kind = rk ) xy(2)

  r = sqrt ( xy(1) * xy(1) + xy(2) * xy(2) )

  if ( r == 0.0D+00 ) then
    t = 0.0D+00
  else
    t = r8_atan ( xy(2), xy(1) )
  end if

  return
end
subroutine xyz_to_radec ( p, ra, dec )

!*****************************************************************************80
!
!! XYZ_TO_RADEC converts (X,Y,Z) to right ascension/declination coordinates.
!
!  Discussion:
!
!    Given an XYZ point, compute its distance R from the origin, and
!    regard it as lying on a sphere of radius R, whose axis is the Z
!    axis.
!
!    The right ascension of the point is the "longitude", measured in hours,
!    between 0 and 24, with the X axis having right ascension 0, and the
!    Y axis having right ascension 6.
!
!    Declination measures the angle from the equator towards the north pole,
!    and ranges from -90 (South Pole) to 90 (North Pole).
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
!  Parameters:
!
!    Input, real ( kind = rk ) P(3), the coordinates of a point in 3D.
!
!    Output, real ( kind = rk ) RA, DEC, the corresponding right ascension
!    and declination.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) dec
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) p_norm
  real ( kind = rk ) phi
  real ( kind = rk ) r8_asin
  real ( kind = rk ) r8_atan
  real ( kind = rk ) ra
  real ( kind = rk ) radians_to_degrees
  real ( kind = rk ) theta

  p_norm = sqrt ( sum ( p(1:dim_num)**2 )  )

  if ( p_norm == 0.0D+00 ) then
    dec = 0.0D+00
    ra = 0.0D+00
    return
  end if

  phi = r8_asin ( p(3) / p_norm )

  if ( cos ( phi ) == 0.0D+00 ) then
    theta = 0.0D+00
  else
    theta = r8_atan ( p(2), p(1) )
  end if

  dec = radians_to_degrees ( phi )
  ra = radians_to_degrees ( theta ) / 15.0D+00

  return
end
subroutine xyz_to_rtp ( xyz, r, theta, phi )

!*****************************************************************************80
!
!! XYZ_TO_RTP converts (X,Y,Z) to (R,Theta,Phi) coordinates.
!
!  Discussion:
!
!    Given an XYZ point, compute its distance R from the origin, and
!    regard it as lying on a sphere of radius R, whose axis is the Z
!    axis.
!
!    Theta measures the "longitude" of the point, between 0 and 2 PI.
!
!    PHI measures the angle from the "north pole", between 0 and PI.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) XYZ(3), the coordinates of a point in 3D.
!
!    Output, real ( kind = rk ) R, THETA, PHI, the radius, longitude and
!    declination of the point.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r
  real ( kind = rk ) r8_acos
  real ( kind = rk ) r8_atan
  real ( kind = rk ) phi
  real ( kind = rk ) theta
  real ( kind = rk ) xyz(3)

  r = sqrt ( sum ( xyz(1:3)**2 )  )

  if ( r == 0.0D+00 ) then
    theta = 0.0D+00
    phi = 0.0D+00
    return
  end if

  phi = r8_acos ( xyz(3) / r )

  theta = r8_atan ( xyz(2), xyz(1) )

  return
end
subroutine xyz_to_tp ( xyz, theta, phi )

!*****************************************************************************80
!
!! XYZ_TO_TP converts (X,Y,Z) to (Theta,Phi) coordinates.
!
!  Discussion:
!
!    Given an XYZ point, regard it as lying on a sphere of radius R, 
!    centered at the origin, whose axis is the Z axis.
!
!    We assume that the actual value of R is of no interest, and do
!    not report it.  This is especially appropriate if the point is
!    expected to lie on the unit sphere, for instance.
!
!    THETA measures the "longitude" of the point, between 0 and 2 PI.
!
!    PHI measures the angle from the "north pole", between 0 and PI.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) XYZ(3), the coordinates of a point in 3D.
!
!    Output, real ( kind = rk ) THETA, PHI, the longitude and declination 
!    of the point.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r
  real ( kind = rk ) r8_acos
  real ( kind = rk ) r8_atan
  real ( kind = rk ) phi
  real ( kind = rk ) theta
  real ( kind = rk ) xyz(3)

  r = sqrt ( sum ( xyz(1:3)**2 )  )

  if ( r == 0.0D+00 ) then
    theta = 0.0D+00
    phi = 0.0D+00
    return
  end if

  phi = r8_acos ( xyz(3) / r )

  theta = r8_atan ( xyz(2), xyz(1) )

  return
end
