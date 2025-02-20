function arc_cosine ( c )

!*****************************************************************************80
!
!! ARC_COSINE computes the arc cosine function, with argument truncation.
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
!    02 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) C, the argument.
!
!    Output, real ( kind = rk ) ARC_COSINE, an angle whose cosine is C.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) arc_cosine
  real ( kind = rk ) c
  real ( kind = rk ) c2

  c2 = c
  c2 = max ( c2, -1.0D+00 )
  c2 = min ( c2, +1.0D+00 )

  arc_cosine = acos ( c2 )

  return
end
subroutine ball_unit_sample_3d ( seed, p )

!*****************************************************************************80
!
!! BALL_UNIT_SAMPLE_3D picks a random point in the unit ball in 3D.
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
!    Input/output, integer SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = rk ) P(3), the sample point.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) arc_cosine
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) phi
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) r
  integer seed
  real ( kind = rk ) theta
  real ( kind = rk ) u(dim_num)
  real ( kind = rk ) vdot

  call r8vec_uniform_01 ( dim_num, seed, u )
!
!  Pick a uniformly random VDOT, which must be between -1 and 1.
!  This represents the dot product of the random vector with the Z unit vector.
!
!  Note: this works because the surface area of the sphere between
!  Z and Z + dZ is independent of Z.  So choosing Z uniformly chooses
!  a patch of area uniformly.
!
  vdot = 2.0D+00 * u(1) - 1.0D+00

  phi = arc_cosine ( vdot )
!
!  Pick a uniformly random rotation between 0 and 2 Pi around the
!  axis of the Z vector.
!
  theta = 2.0D+00 * pi * u(2)
!
!  Pick a random radius R.
!
  r = u(3)**( 1.0D+00 / 3.0D+00 )

  p(1) = r * cos ( theta ) * sin ( phi )
  p(2) = r * sin ( theta ) * sin ( phi )
  p(3) = r * cos ( phi )

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
!    17 May 2007
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
  real ( kind = rk ) off_axis_component
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) p_dot_axis
  real ( kind = rk ) p_length
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) r
  real ( kind = rk ) r8vec_length

  axis(1:dim_num) = p2(1:dim_num) - p1(1:dim_num)
  axis_length = r8vec_length ( dim_num, axis )

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

    p_length = r8vec_length ( dim_num, p(1:dim_num) - p1(1:dim_num) )
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
  real ( kind = rk ) off_axis_component
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) p_dot_axis
  real ( kind = rk ) p_length
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) r
  real ( kind = rk ) r8vec_length

  axis(1:dim_num) = p2(1:dim_num) - p1(1:dim_num)
  axis_length = r8vec_length ( dim_num, axis )

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

    p_length = r8vec_length ( dim_num, p(1:dim_num) - p1(1:dim_num) )
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
!    Output, logical INSIDE, is TRUE if the point is inside the cylinder.
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
  real ( kind = rk ) r8vec_length

  axis(1:dim_num) = p2(1:dim_num) - p1(1:dim_num)
  axis_length = r8vec_length ( dim_num, axis )

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

    p_length = r8vec_length ( dim_num, p(1:dim_num) - p1(1:dim_num) )

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
  real ( kind = rk ) axial_length
  real ( kind = rk ) axis(dim_num)
  real ( kind = rk ) axis_length
  real ( kind = rk ) distance
  real ( kind = rk ) off_axis(dim_num)
  real ( kind = rk ) off_axis_component
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) p1(dim_num)
  real ( kind = rk ) p2(dim_num)
  real ( kind = rk ) pn(dim_num)
  real ( kind = rk ) r
  real ( kind = rk ) r8vec_length

  axis(1:dim_num) = p2(1:dim_num) - p1(1:dim_num)
  axis_length = r8vec_length ( dim_num, axis )
  axis(1:dim_num) = axis(1:dim_num) / axis_length

  axial_component = dot_product ( p(1:dim_num) - p1(1:dim_num), axis )

  off_axis(1:dim_num) = p(1:dim_num) - p1(1:dim_num) &
    - axial_component * axis(1:dim_num)

  off_axis_component = r8vec_length ( dim_num, off_axis )
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
subroutine disk_point_dist_3d ( pc, r, axis, p, dist )

!*****************************************************************************80
!
!! DISK_POINT_DIST_3D determines the distance from a disk to a point in 3D.
!
!  Discussion:
!
!    A disk in 3D satisfies the equations:
!
!      ( P(1) - PC(1) )**2 + ( P(2) - PC(2) )**2 + ( P(3) - PC(3) <= R**2
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
  real ( kind = rk ) off_axis_component
  real ( kind = rk ) off_axis(dim_num)
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) r
  real ( kind = rk ) r8vec_length
!
!  Special case: the point is the center.
!
  if ( all ( p(1:dim_num) == pc(1:dim_num) ) ) then
    dist = 0.0D+00
    return
  end if

  axis_length = r8vec_length ( dim_num, axis(1:dim_num) )

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

  off_axis_component = r8vec_length ( dim_num, off_axis )
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
subroutine p00_boundary_edge_num ( test, boundary_edge_num )

!*****************************************************************************80
!
!! P00_BOUNDARY_EDGE_NUM counts the boundary edges in a problem.
!
!  Discussion:
!
!    Call this routine if you need to allocate memory for an
!    array in which to store the boundary edges.
!
!    We assume that the boundary is triangulated.  BOUNDARY_EDGE_NUM
!    is the number of edges associated with this triangulation.
!    Thus, for a cube, the value is NOT the expected 12, but rather
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer TEST, the index of the test problem
!
!    Output, integer BOUNDARY_EDGE_NUM, the number of
!    boundary edges.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer boundary_edge_num
  integer test

  if ( test == 1 ) then
    call p01_boundary_edge_num ( boundary_edge_num )
  else if ( test == 2 ) then
    call p02_boundary_edge_num ( boundary_edge_num )
  else if ( test == 3 ) then
    call p03_boundary_edge_num ( boundary_edge_num )
  else if ( test == 4 ) then
    call p04_boundary_edge_num ( boundary_edge_num )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_BOUNDARY_EDGE_NUM - Fatal error!'
    write ( *, '(a)' ) '  Input value of TEST is out of range.'
    stop
  end if

  return
end
subroutine p00_boundary_edges ( test, boundary_edge_num, boundary_edge )

!*****************************************************************************80
!
!! P00_BOUNDARY_EDGES returns the boundary edges for any problem.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer TEST, the index of the test problem
!
!    Input, integer BOUNDARY_EDGE_NUM, the number of 
!    boundary edges.
!
!    Output, integer BOUNDARY_EDGE(2,BOUNDARY_EDGE_NUM),
!    the boundary edges, described as pairs of node indices.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer boundary_edge_num

  integer boundary_edge(2,boundary_edge_num)
  integer test

  if ( test == 1 ) then
    call p01_boundary_edges ( boundary_edge_num, boundary_edge )
  else if ( test == 2 ) then
    call p02_boundary_edges ( boundary_edge_num, boundary_edge )
  else if ( test == 3 ) then
    call p03_boundary_edges ( boundary_edge_num, boundary_edge )
  else if ( test == 4 ) then
    call p04_boundary_edges ( boundary_edge_num, boundary_edge )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_BOUNDARY_EDGES - Fatal error!'
    write ( *, '(a)' ) '  Input value of TEST is out of range.'
    stop
  end if

  return
end
subroutine p00_boundary_face_num ( test, boundary_face_num )

!*****************************************************************************80
!
!! P00_BOUNDARY_FACE_NUM counts the boundary faces in a problem.
!
!  Discussion:
!
!    Call this routine if you need to allocate memory for an
!    array in which to store the boundary faces.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer TEST, the index of the test problem
!
!    Output, integer BOUNDARY_FACE_NUM, the number of 
!    boundary faces.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer boundary_face_num
  integer test

  if ( test == 1 ) then
    call p01_boundary_face_num ( boundary_face_num )
  else if ( test == 2 ) then
    call p02_boundary_face_num ( boundary_face_num )
  else if ( test == 3 ) then
    call p03_boundary_face_num ( boundary_face_num )
  else if ( test == 4 ) then
    call p04_boundary_face_num ( boundary_face_num )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_BOUNDARY_FACE_NUM - Fatal error!'
    write ( *, '(a)' ) '  Input value of TEST is out of range.'
    stop
  end if

  return
end
subroutine p00_boundary_faces ( test, boundary_face_num, boundary_face )

!*****************************************************************************80
!
!! P00_BOUNDARY_FACES returns the boundary faces in any problem.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer TEST, the index of the test problem
!
!    Input, integer BOUNDARY_FACE_NUM, the number of
!    boundary faces.
!
!    Output, integer BOUNDARY_FACE(3,BOUNDARY_FACE_NUM),
!    the boundary faces, described as triples of node indices.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer boundary_face_num

  integer boundary_face(3,boundary_face_num)
  integer test

  if ( test == 1 ) then
    call p01_boundary_faces ( boundary_face_num, boundary_face )
  else if ( test == 2 ) then
    call p02_boundary_faces ( boundary_face_num, boundary_face )
  else if ( test == 3 ) then
    call p03_boundary_faces ( boundary_face_num, boundary_face )
  else if ( test == 4 ) then
    call p04_boundary_faces ( boundary_face_num, boundary_face )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_BOUNDARY_FACES - Fatal error!'
    write ( *, '(a)' ) '  Input value of TEST is out of range.'
    stop
  end if

  return
end
subroutine p00_boundary_node_num ( test, boundary_node_num )

!*****************************************************************************80
!
!! P00_BOUNDARY_NODE_NUM counts the boundary nodes in a problem.
!
!  Discussion:
!
!    Call this routine if you need to allocate memory for an
!    array in which to store the boundary nodes.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer TEST, the index of the test problem
!
!    Output, integer BOUNDARY_NODE_NUM, the number of 
!    boundary nodes.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer boundary_node_num
  integer test

  if ( test == 1 ) then
    call p01_boundary_node_num ( boundary_node_num )
  else if ( test == 2 ) then
    call p02_boundary_node_num ( boundary_node_num )
  else if ( test == 3 ) then
    call p03_boundary_node_num ( boundary_node_num )
  else if ( test == 4 ) then
    call p04_boundary_node_num ( boundary_node_num )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_BOUNDARY_NODE_NUM - Fatal error!'
    write ( *, '(a)' ) '  Input value of TEST is out of range.'
    stop
  end if

  return
end
subroutine p00_boundary_nodes ( test, boundary_node_num, boundary_node )

!*****************************************************************************80
!
!! P00_BOUNDARY_NODES returns the boundary nodes in problem 01.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer TEST, the index of the test problem
!
!    Input, integer BOUNDARY_NODE_NUM, the number of 
!    boundary nodes.
!
!    Output, real ( kind = rk ) BOUNDARY_NODE(3,BOUNDARY_NODE_NUM),
!    the boundary nodes.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer boundary_node_num
  integer, parameter :: dim_num = 3

  real ( kind = rk ) boundary_node(dim_num,boundary_node_num)
  integer test

  if ( test == 1 ) then
    call p01_boundary_nodes ( boundary_node_num, boundary_node )
  else if ( test == 2 ) then
    call p02_boundary_nodes ( boundary_node_num, boundary_node )
  else if ( test == 3 ) then
    call p03_boundary_nodes ( boundary_node_num, boundary_node )
  else if ( test == 4 ) then
    call p04_boundary_nodes ( boundary_node_num, boundary_node )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_BOUNDARY_NODES - Fatal error!'
    write ( *, '(a)' ) '  Input value of TEST is out of range.'
    stop
  end if

  return
end
subroutine p00_boundary_project ( test, n, point )

!*****************************************************************************80
!
!! P00_BOUNDARY_PROJECT projects exterior points to the boundary.
!
!  Discussion:
!
!    Interior points are not changed.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer TEST, the index of the test problem
!
!    Input, integer N, the number of points.
!
!    Input/output, real ( kind = rk ) POINT(3,N), the coordinates
!    of the points.  Any input points that are exterior to the region
!    are replaced on output by the nearest boundary point.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer n

  real ( kind = rk ), dimension ( dim_num, n ) :: point
  integer test

  if ( test == 1 ) then
    call p01_boundary_project ( n, point )
  else if ( test == 2 ) then
    call p02_boundary_project ( n, point )
  else if ( test == 3 ) then
    call p03_boundary_project ( n, point )
  else if ( test == 4 ) then
    call p04_boundary_project ( n, point )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_BOUNDARY_PROJECT - Fatal error!'
    write ( *, '(a)' ) '  Input value of TEST is out of range.'
    stop
  end if

  return
end
subroutine p00_box ( test, lo, hi )

!*****************************************************************************80
!
!! P00_BOX returns a bounding box for a problem.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer TEST, the index of the test problem
!
!    Output, real ( kind = rk ) LO(DIMENSION), HI(DIMENSION), the lower and
!    upper corners of a bounding box.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) hi(dim_num)
  real ( kind = rk ) lo(dim_num)
  integer test

  if ( test == 1 ) then
    call p01_box ( lo, hi )
  else if ( test == 2 ) then
    call p02_box ( lo, hi )
  else if ( test == 3 ) then
    call p03_box ( lo, hi )
  else if ( test == 4 ) then
    call p04_box ( lo, hi )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_BOX - Fatal error!'
    write ( *, '(a)' ) '  Input value of TEST is out of range.'
    stop
  end if

  return
end
subroutine p00_fixed_num ( test, fixed_num )

!*****************************************************************************80
!
!! P00_FIXED_NUM returns the number of fixed points in a problem.
!
!  Discussion:
!
!    Call this routine if you need to allocate memory for an
!    array in which to store the fixed points.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer TEST, the index of the test problem
!
!    Output, integer FIXED_NUM, the number of fixed points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer fixed_num
  integer test

  if ( test == 1 ) then
    call p01_fixed_num ( fixed_num )
  else if ( test == 2 ) then
    call p02_fixed_num ( fixed_num )
  else if ( test == 3 ) then
    call p03_fixed_num ( fixed_num )
  else if ( test == 4 ) then
    call p04_fixed_num ( fixed_num )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_FIXED_NUM - Fatal error!'
    write ( *, '(a)' ) '  Input value of TEST is out of range.'
    stop
  end if

  return
end
subroutine p00_fixed_points ( test, fixed_num, fixed )

!*****************************************************************************80
!
!! P00_FIXED_POINTS returns the fixed points in a problem.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer TEST, the index of the test problem
!
!    Input, integer FIXED_NUM, the number of fixed points.
!
!    Output, real ( kind = rk ) FIXED(3,FIXED_NUM), the
!    coordinates of the fixed points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer fixed_num

  real ( kind = rk ) fixed(dim_num,fixed_num)
  integer test

  if ( test == 1 ) then
    call p01_fixed_points ( fixed_num, fixed )
  else if ( test == 2 ) then
    call p02_fixed_points ( fixed_num, fixed )
  else if ( test == 3 ) then
    call p03_fixed_points ( fixed_num, fixed )
  else if ( test == 4 ) then
    call p04_fixed_points ( fixed_num, fixed )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_FIXED_POINTS - Fatal error!'
    write ( *, '(a)' ) '  Input value of TEST is out of range.'
    stop
  end if

  return
end
subroutine p00_header ( test )

!*****************************************************************************80
!
!! P00_HEADER prints some information about a problem.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer TEST, the index of the test problem
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer test

  if ( test == 1 ) then
    call p01_header ( )
  else if ( test == 2 ) then
    call p02_header ( )
  else if ( test == 3 ) then
    call p03_header ( )
  else if ( test == 4 ) then
    call p04_header ( )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_HEADER - Fatal error!'
    write ( *, '(a)' ) '  Input value of TEST is out of range.'
    stop
  end if

  return
end
subroutine p00_inside ( test, n, point, inside )

!*****************************************************************************80
!
!! P00_INSIDE reports if a point is inside the region in a problem.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer TEST, the index of the test problem
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, logical INSIDE(N), is TRUE if the point is in the region.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer n

  logical inside(n)
  real ( kind = rk ) point(dim_num,n)
  integer test

  if ( test == 1 ) then
    call p01_inside ( n, point, inside )
  else if ( test == 2 ) then
    call p02_inside ( n, point, inside )
  else if ( test == 3 ) then
    call p03_inside ( n, point, inside )
  else if ( test == 4 ) then
    call p04_inside ( n, point, inside )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_INSIDE - Fatal error!'
    write ( *, '(a)' ) '  Input value of TEST is out of range.'
    stop
  end if

  return
end
subroutine p00_sample ( test, n, seed, point )

!*****************************************************************************80
!
!! P00_SAMPLE samples points from the region in a problem.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer TEST, the index of the test problem
!
!    Input, integer N, the number of points.
!
!    Input/output, integer SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = rk ) POINT(3,N), the coordinates
!    of the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer n

  real ( kind = rk ) point(dim_num,n)
  integer seed
  integer test

  if ( test == 1 ) then
    call p01_sample ( n, seed, point )
  else if ( test == 2 ) then
    call p02_sample ( n, seed, point )
  else if ( test == 3 ) then
    call p03_sample ( n, seed, point )
  else if ( test == 4 ) then
    call p04_sample ( n, seed, point )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_SAMPLE - Fatal error!'
    write ( *, '(a)' ) '  Input value of TEST is out of range.'
    stop
  end if

  return
end
subroutine p00_sample_h1 ( test, n, h, seed, point )

!*****************************************************************************80
!
!! P00_SAMPLE_H1 samples points from the enlarged region in a problem.
!
!  Discussion:
!
!    The region is enlarged by an amount H.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer TEST, the index of the test problem
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) H, the enlargement amount.
!
!    Input/output, integer SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = rk ) POINT(3,N), the coordinates
!    of the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer n

  real ( kind = rk ) h
  real ( kind = rk ) point(dim_num,n)
  integer seed
  integer test

  if ( test == 1 ) then
    call p01_sample_h1 ( n, h, seed, point )
  else if ( test == 2 ) then
    call p02_sample_h1 ( n, h, seed, point )
  else if ( test == 3 ) then
    call p03_sample_h1 ( n, h, seed, point )
  else if ( test == 4 ) then
    call p04_sample_h1 ( n, h, seed, point )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_SAMPLE_H1 - Fatal error!'
    write ( *, '(a)' ) '  Input value of TEST is out of range.'
    stop
  end if

  return
end
subroutine p00_sdist ( test, n, point, sdist )

!*****************************************************************************80
!
!! P00_SDIST returns the signed distance to the region in a problem.
!
!  Discussion:
!
!    A positive distance indicates the point is outside the region.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer TEST, the index of the test problem
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) POINT(M,N), the coordinates
!    of the points.
!
!    Output, real ( kind = rk ) SDIST(N), the signed distance
!    of each point to the region.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer n

  real ( kind = rk ) point(dim_num,n)
  real ( kind = rk ) sdist(n)
  integer test

  if ( test == 1 ) then
    call p01_sdist ( n, point, sdist )
  else if ( test == 2 ) then
    call p02_sdist ( n, point, sdist )
  else if ( test == 3 ) then
    call p03_sdist ( n, point, sdist )
  else if ( test == 4 ) then
    call p04_sdist ( n, point, sdist )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_SDIST - Fatal error!'
    write ( *, '(a)' ) '  Input value of TEST is out of range.'
    stop
  end if

  return
end
subroutine p00_test_num ( test_num )

!*****************************************************************************80
!
!! P00_TEST_NUM returns the number of available tests.
!
!  Discussion:
!
!    Call this routine if you need to cycle through all the tests.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer TEST_NUM, the number of tests.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer test_num

  test_num = 4

  return
end
subroutine p00_title ( test, title )

!*****************************************************************************80
!
!! P00_TITLE returns a title for a problem.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer TEST, the index of the test problem
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer test
  character ( len = * ) title

  if ( test == 1 ) then
    call p01_title ( title )
  else if ( test == 2 ) then
    call p02_title ( title )
  else if ( test == 3 ) then
    call p03_title ( title )
  else if ( test == 4 ) then
    call p04_title ( title )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_TITLE - Fatal error!'
    write ( *, '(a)' ) '  Input value of TEST is out of range.'
    stop
  end if

  return
end
subroutine p01_boundary_edge_num ( boundary_edge_num )

!*****************************************************************************80
!
!! P01_BOUNDARY_EDGE_NUM counts the boundary edges in problem 01.
!
!  Discussion:
!
!    Call this routine if you need to allocate memory for an
!    array in which to store the boundary edges.
!
!    We assume that the boundary is triangulated.  BOUNDARY_EDGE_NUM
!    is the number of edges associated with this triangulation.
!    Thus, for a cube, the value is NOT the expected 12, but rather
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
!    Output, integer BOUNDARY_EDGE_NUM, the number of 
!    boundary edges.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer boundary_edge_num

  boundary_edge_num = 18

  return
end
subroutine p01_boundary_edges ( boundary_edge_num, boundary_edge )

!*****************************************************************************80
!
!! P01_BOUNDARY_EDGES returns the boundary edges in problem 01.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer BOUNDARY_EDGE_NUM, the number of 
!    boundary edges.
!
!    Output, integer BOUNDARY_EDGE(2,BOUNDARY_EDGE_NUM),
!    the boundary edges, described as pairs of node indices.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer boundary_edge_num

  integer boundary_edge(2,boundary_edge_num)

  boundary_edge = reshape ( (/ &
    1, 2, &
    1, 3, &
    1, 5, &
    1, 6, &
    2, 3, &
    2, 4, &
    2, 5, &
    2, 7, &
    3, 4, &
    3, 6, &
    3, 8, &
    4, 7, &
    4, 8, &
    5, 6, &
    5, 7, &
    5, 8, &
    6, 8, &
    7, 8 /), (/ 2, boundary_edge_num /) )

  return
end
subroutine p01_boundary_face_num ( boundary_face_num )

!*****************************************************************************80
!
!! P01_BOUNDARY_FACE_NUM counts the boundary faces in problem 01.
!
!  Discussion:
!
!    Call this routine if you need to allocate memory for an
!    array in which to store the boundary faces.
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
!    Output, integer BOUNDARY_FACE_NUM, the number of 
!    boundary faces.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer boundary_face_num

  boundary_face_num = 12

  return
end
subroutine p01_boundary_faces ( boundary_face_num, boundary_face )

!*****************************************************************************80
!
!! P01_BOUNDARY_FACES returns the boundary faces in problem 01.
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
!    Input, integer BOUNDARY_FACE_NUM, the number of 
!    boundary faces.
!
!    Output, integer BOUNDARY_FACE(3,BOUNDARY_FACE_NUM),
!    the boundary faces, described as triples of node indices.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer boundary_face_num

  integer boundary_face(3,boundary_face_num)

  boundary_face = reshape ( (/ &
    5, 7, 2, &
    5, 2, 1, &
    8, 6, 3, &
    8, 4, 3, &
    8, 7, 5, &
    8, 5, 6, &
    6, 5, 1, &
    6, 1, 3, &
    4, 2, 7, &
    4, 7, 8, &
    3, 1, 2, &
    3, 2, 4 /), (/ 2, boundary_face_num /) )

  return
end
subroutine p01_boundary_node_num ( boundary_node_num )

!*****************************************************************************80
!
!! P01_BOUNDARY_NODE_NUM counts the boundary nodes in problem 01.
!
!  Discussion:
!
!    Call this routine if you need to allocate memory for an
!    array in which to store the boundary nodes.
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
!    Output, integer BOUNDARY_NODE_NUM, the number of 
!    boundary nodes.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer boundary_node_num

  boundary_node_num = 8

  return
end
subroutine p01_boundary_nodes ( boundary_node_num, boundary_node )

!*****************************************************************************80
!
!! P01_BOUNDARY_NODES returns the boundary nodes in problem 01.
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
!    Input, integer BOUNDARY_NODE_NUM, the number of 
!    boundary nodes.
!
!    Output, real ( kind = rk ) BOUNDARY_NODE(3,BOUNDARY_NODE_NUM),
!    the boundary nodes.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer boundary_node_num
  integer, parameter :: dim_num = 3

  real ( kind = rk ) boundary_node(dim_num,boundary_node_num)

  boundary_node = reshape ( (/ &
    0.0D+00, 0.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, 1.0D+00, &
    0.0D+00, 1.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00, 1.0D+00, &
    3.0D+00, 0.0D+00, 0.0D+00, &
    3.0D+00, 0.0D+00, 1.0D+00, &
    3.0D+00, 1.0D+00, 0.0D+00, &
    3.0D+00, 1.0D+00, 1.0D+00 /), (/ dim_num, boundary_node_num /) )

  return
end
subroutine p01_boundary_project ( n, point )

!*****************************************************************************80
!
!! P01_BOUNDARY_PROJECT projects exterior points to the boundary in problem 01.
!
!  Discussion:
!
!    Points in the interior are not changed.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!
!    Input/output, real ( kind = rk ) POINT(3,N), the coordinates
!    of the points.  On output, all exterior points have been
!    replaced by the nearest point on the boundary.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer n

  real ( kind = rk ), dimension ( dim_num, n ) :: point
  real ( kind = rk ), parameter :: xlo =  0.0D+00
  real ( kind = rk ), parameter :: xhi = +3.0D+00
  real ( kind = rk ), parameter :: ylo =  0.0D+00
  real ( kind = rk ), parameter :: yhi = +1.0D+00
  real ( kind = rk ), parameter :: zlo =  0.0D+00
  real ( kind = rk ), parameter :: zhi = +1.0D+00

  point(1,1:n) = max ( point(1,1:n), xlo )
  point(1,1:n) = min ( point(1,1:n), xhi )
  point(2,1:n) = max ( point(2,1:n), ylo )
  point(2,1:n) = min ( point(2,1:n), yhi )
  point(3,1:n) = max ( point(3,1:n), zlo )
  point(3,1:n) = min ( point(3,1:n), zhi )

  return
end
subroutine p01_box ( lo, hi )

!*****************************************************************************80
!
!! P01_BOX returns a bounding box for problem 01.
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
!  Parameters:
!
!    Output, real ( kind = rk ) LO(3), HI(3), coordinates of the
!    low and high corners of the box.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) hi(dim_num)
  real ( kind = rk ) lo(dim_num)

  lo(1:dim_num) = (/ 0.0D+00, 0.0D+00, 0.0D+00 /)
  hi(1:dim_num) = (/ 3.0D+00, 1.0D+00, 1.0D+00 /)

  return
end
subroutine p01_fixed_num ( fixed_num )

!*****************************************************************************80
!
!! P01_FIXED_NUM returns the number of fixed points in problem 01.
!
!  Discussion:
!
!    Call this routine if you need to allocate memory for an
!    array in which to store the fixed points.
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
!    Output, integer FIXED_NUM, the number of fixed points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer fixed_num

  fixed_num = 8

  return
end
subroutine p01_fixed_points ( fixed_num, fixed )

!*****************************************************************************80
!
!! P01_FIXED_POINTS returns the fixed points in problem 01.
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
!    Input, integer FIXED_NUM, the number of fixed points.
!
!    Output, real ( kind = rk ) FIXED(3,FIXED_NUM), the coordinates
!    of the fixed points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer fixed_num

  real ( kind = rk ) fixed(dim_num,fixed_num)

  fixed(1:dim_num,1) = (/ 0.0D+00, 0.0D+00, 0.0D+00 /)
  fixed(1:dim_num,2) = (/ 0.0D+00, 0.0D+00, 1.0D+00 /)
  fixed(1:dim_num,3) = (/ 0.0D+00, 1.0D+00, 0.0D+00 /)
  fixed(1:dim_num,4) = (/ 0.0D+00, 1.0D+00, 1.0D+00 /)
  fixed(1:dim_num,5) = (/ 3.0D+00, 0.0D+00, 0.0D+00 /)
  fixed(1:dim_num,6) = (/ 3.0D+00, 0.0D+00, 1.0D+00 /)
  fixed(1:dim_num,7) = (/ 3.0D+00, 1.0D+00, 0.0D+00 /)
  fixed(1:dim_num,8) = (/ 3.0D+00, 1.0D+00, 1.0D+00 /)

  return
end
subroutine p01_header ( )

!*****************************************************************************80
!
!! P01_HEADER prints some information about problem 01.
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
!  Parameters:
!
!    None.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'P01:'
  write ( *, '(a)' ) '  Sample problem #1'
  write ( *, '(a)' ) '  The channel.'
  write ( *, '(a)' ) '  3 x 1 x 1 box.'

  return
end
subroutine p01_inside ( n, point, inside )

!*****************************************************************************80
!
!! P01_INSIDE reports if a point is inside the region in problem 01.
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
!  Parameters:
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) POINT(3,N), the coordinates
!    of the points.
!
!    Output, logical INSIDE(N), is TRUE if the point is in the region.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer n

  logical inside(n)
  integer j
  real ( kind = rk ) point(dim_num,n)

  do j = 1, n
    inside(j) = 0.0D+00 <= point(1,j) .and. point(1,j) <= 3.0D+00 .and. &
                0.0D+00 <= point(2,j) .and. point(2,j) <= 1.0D+00 .and. &
                0.0D+00 <= point(3,j) .and. point(3,j) <= 1.0D+00

  end do

  return
end
subroutine p01_sample ( n, seed, point )

!*****************************************************************************80
!
!! P01_SAMPLE samples points from the region in problem 01.
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
!  Parameters:
!
!    Input, integer N, the number of points.
!
!    Input/output, integer SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = rk ) POINT(3,N), the coordinates
!    of the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer n

  real ( kind = rk ) point(dim_num,n)
  integer seed

  call random_number ( harvest = point(1:dim_num,1:n) )

  point(1,1:n) = 3.0D+00 * point(1,1:n)

  return
end
subroutine p01_sample_h1 ( n, h, seed, point )

!*****************************************************************************80
!
!! P01_SAMPLE_H1 samples points from the enlarged region in problem 01.
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
!  Parameters:
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) H, the enlargement amount.
!
!    Input/output, integer SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = rk ) POINT(3,N), the coordinates
!    of the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer n

  real ( kind = rk ) h
  real ( kind = rk ) point(dim_num,n)
  integer seed

  call random_number ( harvest = point(1:dim_num,1:n) )

  point(1,1:n) = ( 1.0D+00 - point(1,1:n) ) * (         - h ) &
               +             point(1,1:n)   * ( 3.0D+00 + h )

  point(2,1:n) = ( 1.0D+00 - point(2,1:n) ) * (         - h ) &
               +             point(2,1:n)   * ( 1.0D+00 + h )

  point(3,1:n) = ( 1.0D+00 - point(3,1:n) ) * (         - h ) &
               +             point(3,1:n)   * ( 1.0D+00 + h )

  return
end
subroutine p01_sdist ( n, point, sdist )

!*****************************************************************************80
!
!! P01_SDIST returns the signed distance to the region in problem 01.
!
!  Discussion:
!
!    A positive distance indicates the point is outside the region.
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
!  Parameters:
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) POINT(3,N), the coordinates
!    of the points.
!
!    Output, real ( kind = rk ) SDIST(N), the signed distance of
!    each point to the region.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer n

  integer j
  real ( kind = rk ) point(dim_num,n)
  real ( kind = rk ) sdist(n)

  do j = 1, n
    sdist(j) = min ( -point(1,j), 3.0D+00 - point(1,j), &
                     -point(2,j), 1.0D+00 - point(2,j), &
                     -point(3,j), 1.0D+00 - point(3,j) )

  end do

  return
end
subroutine p01_title ( title )

!*****************************************************************************80
!
!! P01_TITLE returns a title for problem 01.
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
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = * ) title

  title = '#1: The 3x1x1 box.'

  return
end
subroutine p02_boundary_edge_num ( boundary_edge_num )

!*****************************************************************************80
!
!! P02_BOUNDARY_EDGE_NUM counts the boundary edges in problem 02.
!
!  Discussion:
!
!    Call this routine if you need to allocate memory for an
!    array in which to store the boundary edges.
!
!    We assume that the boundary is triangulated.  BOUNDARY_EDGE_NUM
!    is the number of edges associated with this triangulation.
!    Thus, for a cube, the value is NOT the expected 12, but rather
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
!    Output, integer BOUNDARY_EDGE_NUM, the number of 
!    boundary edges.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer boundary_edge_num

  boundary_edge_num = 0

  return
end
subroutine p02_boundary_edges ( boundary_edge_num, boundary_edge )

!*****************************************************************************80
!
!! P02_BOUNDARY_EDGES returns the boundary edges in problem 02.
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
!    Input, integer BOUNDARY_EDGE_NUM, the number of 
!    boundary edges.
!
!    Output, integer BOUNDARY_EDGE(3,BOUNDARY_EDGE_NUM),
!    the boundary edges, described as pairs of node indices.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer boundary_edge_num

  integer boundary_edge(2,boundary_edge_num)

  return
end
subroutine p02_boundary_face_num ( boundary_face_num )

!*****************************************************************************80
!
!! P02_BOUNDARY_FACE_NUM counts the boundary faces in problem 02.
!
!  Discussion:
!
!    Call this routine if you need to allocate memory for an
!    array in which to store the boundary faces.
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
!    Output, integer BOUNDARY_FACE_NUM, the number of 
!    boundary faces.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer boundary_face_num

  boundary_face_num = 0

  return
end
subroutine p02_boundary_faces ( boundary_face_num, boundary_face )

!*****************************************************************************80
!
!! P02_BOUNDARY_FACES returns the boundary faces in problem 02.
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
!    Input, integer BOUNDARY_FACE_NUM, the number of 
!    boundary faces.
!
!    Output, integer BOUNDARY_FACE(3,BOUNDARY_FACE_NUM),
!    the boundary faces, described as triples of node indices.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer boundary_face_num

  integer boundary_face(3,boundary_face_num)

  return
end
subroutine p02_boundary_node_num ( boundary_node_num )

!*****************************************************************************80
!
!! P02_BOUNDARY_NODE_NUM counts the boundary nodes in problem 02.
!
!  Discussion:
!
!    Call this routine if you need to allocate memory for an
!    array in which to store the boundary nodes.
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
!    Output, integer BOUNDARY_NODE_NUM, the number of 
!    boundary nodes.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer boundary_node_num

  boundary_node_num = 0

  return
end
subroutine p02_boundary_nodes ( boundary_node_num, boundary_node )

!*****************************************************************************80
!
!! P02_BOUNDARY_NODES returns the boundary nodes in problem 02.
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
!    Input, integer BOUNDARY_NODE_NUM, the number of 
!    boundary nodes.
!
!    Output, real ( kind = rk ) BOUNDARY_NODE(3,BOUNDARY_NODE_NUM),
!    the boundary nodes.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer boundary_node_num
  integer, parameter :: dim_num = 3

  real ( kind = rk ) boundary_node(dim_num,boundary_node_num)

  return
end
subroutine p02_boundary_project ( n, point )

!*****************************************************************************80
!
!! P02_BOUNDARY_PROJECT projects exterior points to the boundary in problem 02.
!
!  Discussion:
!
!    Points in the interior are not changed.
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
!    Input, integer N, the number of points.
!
!    Input/output, real ( kind = rk ) POINT(3,N), the coordinates
!    of the points.  On output, all exterior points have been
!    replaced by the nearest point on the boundary.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer n
  real ( kind = rk ), parameter :: height = 4.0D+00

  integer j
  real ( kind = rk ), dimension(dim_num) :: p1 = (/ &
    0.0D+00, 0.0D+00, 0.0D+00 /)
  real ( kind = rk ), dimension(dim_num) :: p2 = (/ &
    0.0D+00, 0.0D+00, height /)
  real ( kind = rk ) pn(dim_num)
  real ( kind = rk ), dimension ( dim_num, n ) :: point
  real ( kind = rk ), parameter :: radius = 1.0D+00

  do j = 1, n

    call cylinder_point_near_3d ( p1, p2, radius, point(1:dim_num,j), &
      pn(1:dim_num) )

    point(1:dim_num,j) = pn(1:dim_num)

  end do

  return
end
subroutine p02_box ( lo, hi )

!*****************************************************************************80
!
!! P02_BOX returns a bounding box for problem 02.
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
!    Output, real ( kind = rk ) LO(3), HI(3), coordinates of the
!    low and high corners of the box.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ), parameter :: height = 4.0D+00
  real ( kind = rk ) hi(dim_num)
  real ( kind = rk ) lo(dim_num)
  real ( kind = rk ), parameter :: radius = 1.0D+00

  lo(1:dim_num) = (/ -radius, -radius, 0.0D+00 /)
  hi(1:dim_num) = (/  radius,  radius, height /)

  return
end
subroutine p02_fixed_num ( fixed_num )

!*****************************************************************************80
!
!! P02_FIXED_NUM returns the number of fixed points in problem 02.
!
!  Discussion:
!
!    Call this routine if you need to allocate memory for an
!    array in which to store the fixed points.
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
!    Output, integer FIXED_NUM, the number of fixed points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer fixed_num

  fixed_num = 0

  return
end
subroutine p02_fixed_points ( fixed_num, fixed )

!*****************************************************************************80
!
!! P02_FIXED_POINTS returns the fixed points in problem 02.
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
!    Input, integer FIXED_NUM, the number of fixed points.
!
!    Output, real ( kind = rk ) FIXED(3,FIXED_NUM), the coordinates
!    of the fixed points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer fixed_num

  real ( kind = rk ) fixed(dim_num,fixed_num)

  return
end
subroutine p02_header ( )

!*****************************************************************************80
!
!! P02_HEADER prints some information about problem 02.
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
!    None.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), parameter :: height = 4.0D+00
  real ( kind = rk ), parameter :: radius = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'P02:'
  write ( *, '(a)' ) '  Sample problem #2'
  write ( *, '(a)' ) '  The vertical cylinder.'
  write ( *, '(a,g14.6)' ) '  radius R = ', radius
  write ( *, '(a,g14.6)' ) '  height H = ', height

  return
end
subroutine p02_inside ( n, point, inside )

!*****************************************************************************80
!
!! P02_INSIDE reports if a point is inside the region in problem 02.
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
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) POINT(3,N), the coordinates
!    of the points.
!
!    Output, logical INSIDE(N), is TRUE if the point is in the region.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  real ( kind = rk ), parameter :: height = 4.0D+00
  integer n

  logical inside(n)
  integer j
  real ( kind = rk ), dimension(dim_num) :: p1 = (/ &
    0.0D+00, 0.0D+00, 0.0D+00 /)
  real ( kind = rk ), dimension(dim_num) :: p2 = (/ &
    0.0D+00, 0.0D+00, height /)
  real ( kind = rk ) point(dim_num,n)
  real ( kind = rk ), parameter :: radius = 1.0D+00

  do j = 1, n

    call cylinder_point_inside_3d ( p1, p2, radius, point(1:dim_num,j), &
      inside(j) )

  end do

  return
end
subroutine p02_sample ( n, seed, point )

!*****************************************************************************80
!
!! P02_SAMPLE samples points from the region in problem 02.
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
!    Input, integer N, the number of points.
!
!    Input/output, integer SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = rk ) POINT(3,N), the coordinates
!    of the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer n

  real ( kind = rk ), parameter :: height = 4.0D+00
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) point(dim_num,n)
  real ( kind = rk ) r(n)
  real ( kind = rk ), parameter :: radius = 1.0D+00
  integer seed
  real ( kind = rk ) theta(n)
  real ( kind = rk ) z(n)

  call random_number ( harvest = r(1:n) )
  r(1:n) = radius * sqrt ( r(1:n) )

  call random_number ( harvest = theta(1:n) )
  theta(1:n) = 2.0D+00 * pi * theta(1:n)

  call random_number ( harvest = z(1:n) )
  z(1:n) = height * z(1:n)

  point(1,1:n) = r(1:n) * cos ( theta(1:n) )
  point(2,1:n) = r(1:n) * sin ( theta(1:n) )
  point(3,1:n) = z(1:n)

  return
end
subroutine p02_sample_h1 ( n, h, seed, point )

!*****************************************************************************80
!
!! P02_SAMPLE_H1 samples points from the enlarged region in problem 02.
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
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) H, the enlargement amount.
!
!    Input/output, integer SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = rk ) POINT(3,N), the coordinates
!    of the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer n

  real ( kind = rk ) h
  real ( kind = rk ), parameter :: height = 4.0D+00
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) point(dim_num,n)
  real ( kind = rk ) r(n)
  real ( kind = rk ), parameter :: radius = 1.0D+00
  integer seed
  real ( kind = rk ) theta(n)
  real ( kind = rk ) z(n)

  call random_number ( harvest = r(1:n) )
  r(1:n) = radius * sqrt ( r(1:n) )

  call random_number ( harvest = theta(1:n) )
  theta(1:n) = 2.0D+00 * pi * theta(1:n)

  call random_number ( harvest = z(1:n) )
  z(1:n) = ( height + h ) *             z(1:n) &
                    - h   * ( 1.0D+00 - z(1:n) )

  point(1,1:n) = r(1:n) * cos ( theta(1:n) )
  point(2,1:n) = r(1:n) * sin ( theta(1:n) )
  point(3,1:n) = z(1:n)

  return
end
subroutine p02_sdist ( n, point, sdist )

!*****************************************************************************80
!
!! P02_SDIST returns the signed distance to the region in problem 02.
!
!  Discussion:
!
!    A positive distance indicates the point is outside the region.
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
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) POINT(3,N), the coordinates
!    of the points.
!
!    Output, real ( kind = rk ) SDIST(N), the signed distance of
!    each point to the region.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  real ( kind = rk ), parameter :: height = 4.0D+00
  integer n

  integer j
  real ( kind = rk ), dimension(dim_num) :: p1 = (/ &
    0.0D+00, 0.0D+00, 0.0D+00 /)
  real ( kind = rk ), dimension(dim_num) :: p2 = (/ &
    0.0D+00, 0.0D+00, height /)
  real ( kind = rk ) point(dim_num,n)
  real ( kind = rk ), parameter :: radius = 1.0D+00
  real ( kind = rk ) sdist(n)

  do j = 1, n

    call cylinder_point_dist_signed_3d ( p1, p2, radius, &
      point(1:dim_num,j), sdist(j) )

  end do

  return
end
subroutine p02_title ( title )

!*****************************************************************************80
!
!! P02_TITLE returns a title for problem 02.
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
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = * ) title

  title = '#2: The vertical cylinder.'

  return
end
subroutine p03_boundary_edge_num ( boundary_edge_num )

!*****************************************************************************80
!
!! P03_BOUNDARY_EDGE_NUM counts the boundary edges in problem 03.
!
!  Discussion:
!
!    Call this routine if you need to allocate memory for an
!    array in which to store the boundary edges.
!
!    We assume that the boundary is triangulated.  BOUNDARY_EDGE_NUM
!    is the number of edges associated with this triangulation.
!    Thus, for a cube, the value is NOT the expected 12, but rather
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer BOUNDARY_EDGE_NUM, the number of 
!    boundary edges.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer boundary_edge_num

  boundary_edge_num = 18

  return
end
subroutine p03_boundary_edges ( boundary_edge_num, boundary_edge )

!*****************************************************************************80
!
!! P03_BOUNDARY_EDGES returns the boundary edges in problem 03.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer BOUNDARY_EDGE_NUM, the number of 
!    boundary edges.
!
!    Output, integer BOUNDARY_EDGE(2,BOUNDARY_EDGE_NUM),
!    the boundary edges, described as pairs of node indices.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer boundary_edge_num

  integer boundary_edge(2,boundary_edge_num)

  boundary_edge = reshape ( (/ &
    1, 2, &
    1, 3, &
    1, 5, &
    1, 6, &
    2, 3, &
    2, 4, &
    2, 5, &
    2, 7, &
    3, 4, &
    3, 6, &
    3, 8, &
    4, 7, &
    4, 8, &
    5, 6, &
    5, 7, &
    5, 8, &
    6, 8, &
    7, 8 /), (/ 2, boundary_edge_num /) )

  return
end
subroutine p03_boundary_face_num ( boundary_face_num )

!*****************************************************************************80
!
!! P03_BOUNDARY_FACE_NUM counts the boundary faces in problem 03.
!
!  Discussion:
!
!    Call this routine if you need to allocate memory for an
!    array in which to store the boundary faces.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer BOUNDARY_FACE_NUM, the number of 
!    boundary faces.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer boundary_face_num

  boundary_face_num = 12

  return
end
subroutine p03_boundary_faces ( boundary_face_num, boundary_face )

!*****************************************************************************80
!
!! P03_BOUNDARY_FACES returns the boundary faces in problem 03.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer BOUNDARY_FACE_NUM, the number of 
!    boundary faces.
!
!    Output, integer BOUNDARY_FACE(3,BOUNDARY_FACE_NUM),
!    the boundary faces, described as triples of node indices.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer boundary_face_num

  integer boundary_face(3,boundary_face_num)

  boundary_face = reshape ( (/ &
    5, 7, 2, &
    5, 2, 1, &
    8, 6, 3, &
    8, 4, 3, &
    8, 7, 5, &
    8, 5, 6, &
    6, 5, 1, &
    6, 1, 3, &
    4, 2, 7, &
    4, 7, 8, &
    3, 1, 2, &
    3, 2, 4 /), (/ 2, boundary_face_num /) )

  return
end
subroutine p03_boundary_node_num ( boundary_node_num )

!*****************************************************************************80
!
!! P03_BOUNDARY_NODE_NUM counts the boundary nodes in problem 03.
!
!  Discussion:
!
!    Call this routine if you need to allocate memory for an
!    array in which to store the boundary nodes.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer BOUNDARY_NODE_NUM, the number of
!    boundary nodes.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer boundary_node_num

  boundary_node_num = 8

  return
end
subroutine p03_boundary_nodes ( boundary_node_num, boundary_node )

!*****************************************************************************80
!
!! P03_BOUNDARY_NODES returns the boundary nodes in problem 03.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer BOUNDARY_NODE_NUM, the number of 
!    boundary nodes.
!
!    Output, real ( kind = rk ) BOUNDARY_NODE(3,BOUNDARY_NODE_NUM),
!    the boundary nodes.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer boundary_node_num
  integer, parameter :: dim_num = 3

  real ( kind = rk ) boundary_node(dim_num,boundary_node_num)

  boundary_node = reshape ( (/ &
    0.0D+00, 0.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, 1.0D+00, &
    0.0D+00, 1.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 0.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00 /), (/ dim_num, boundary_node_num /) )

  return
end
subroutine p03_boundary_project ( n, point )

!*****************************************************************************80
!
!! P03_BOUNDARY_PROJECT projects exterior points to the boundary in problem 03.
!
!  Discussion:
!
!    Points in the interior are not changed.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!
!    Input/output, real ( kind = rk ) POINT(3,N), the coordinates
!    of the points.  On output, all exterior points have been
!    replaced by the nearest point on the boundary.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer n

  real ( kind = rk ), dimension ( dim_num, n ) :: point
  real ( kind = rk ), parameter :: xlo =  0.0D+00
  real ( kind = rk ), parameter :: xhi = +1.0D+00
  real ( kind = rk ), parameter :: ylo =  0.0D+00
  real ( kind = rk ), parameter :: yhi = +1.0D+00
  real ( kind = rk ), parameter :: zlo =  0.0D+00
  real ( kind = rk ), parameter :: zhi = +1.0D+00

  point(1,1:n) = max ( point(1,1:n), xlo )
  point(1,1:n) = min ( point(1,1:n), xhi )
  point(2,1:n) = max ( point(2,1:n), ylo )
  point(2,1:n) = min ( point(2,1:n), yhi )
  point(3,1:n) = max ( point(3,1:n), zlo )
  point(3,1:n) = min ( point(3,1:n), zhi )

  return
end
subroutine p03_box ( lo, hi )

!*****************************************************************************80
!
!! P03_BOX returns a bounding box for problem 03.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) LO(3), HI(3), coordinates of the
!    low and high corners of the box.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) hi(dim_num)
  real ( kind = rk ) lo(dim_num)

  lo(1:dim_num) = (/ 0.0D+00, 0.0D+00, 0.0D+00 /)
  hi(1:dim_num) = (/ 1.0D+00, 1.0D+00, 1.0D+00 /)

  return
end
subroutine p03_fixed_num ( fixed_num )

!*****************************************************************************80
!
!! P03_FIXED_NUM returns the number of fixed points in problem 03.
!
!  Discussion:
!
!    Call this routine if you need to allocate memory for an
!    array in which to store the fixed points.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer FIXED_NUM, the number of fixed points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer fixed_num

  fixed_num = 8

  return
end
subroutine p03_fixed_points ( fixed_num, fixed )

!*****************************************************************************80
!
!! P03_FIXED_POINTS returns the fixed points in problem 03.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer FIXED_NUM, the number of fixed points.
!
!    Output, real ( kind = rk ) FIXED(3,FIXED_NUM), the coordinates
!    of the fixed points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer fixed_num

  real ( kind = rk ) fixed(dim_num,fixed_num)

  fixed(1:dim_num,1) = (/ 0.0D+00, 0.0D+00, 0.0D+00 /)
  fixed(1:dim_num,2) = (/ 0.0D+00, 0.0D+00, 1.0D+00 /)
  fixed(1:dim_num,3) = (/ 0.0D+00, 1.0D+00, 0.0D+00 /)
  fixed(1:dim_num,4) = (/ 0.0D+00, 1.0D+00, 1.0D+00 /)
  fixed(1:dim_num,5) = (/ 1.0D+00, 0.0D+00, 0.0D+00 /)
  fixed(1:dim_num,6) = (/ 1.0D+00, 0.0D+00, 1.0D+00 /)
  fixed(1:dim_num,7) = (/ 1.0D+00, 1.0D+00, 0.0D+00 /)
  fixed(1:dim_num,8) = (/ 1.0D+00, 1.0D+00, 1.0D+00 /)

  return
end
subroutine p03_header ( )

!*****************************************************************************80
!
!! P03_HEADER prints some information about problem 03.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'P03:'
  write ( *, '(a)' ) '  Sample problem #3'
  write ( *, '(a)' ) '  The unit cube.'

  return
end
subroutine p03_inside ( n, point, inside )

!*****************************************************************************80
!
!! P03_INSIDE reports if a point is inside the region in problem 03.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) POINT(3,N), the coordinates
!    of the points.
!
!    Output, logical INSIDE(N), is TRUE if the point is in the region.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer n

  logical inside(n)
  integer j
  real ( kind = rk ) point(dim_num,n)

  do j = 1, n
    inside(j) = 0.0D+00 <= point(1,j) .and. point(1,j) <= 1.0D+00 .and. &
                0.0D+00 <= point(2,j) .and. point(2,j) <= 1.0D+00 .and. &
                0.0D+00 <= point(3,j) .and. point(3,j) <= 1.0D+00

  end do

  return
end
subroutine p03_sample ( n, seed, point )

!*****************************************************************************80
!
!! P03_SAMPLE samples points from the region in problem 03.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!
!    Input/output, integer SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = rk ) POINT(3,N), the coordinates
!    of the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer n

  real ( kind = rk ) point(dim_num,n)
  integer seed

  call random_number ( harvest = point(1:dim_num,1:n) )

  return
end
subroutine p03_sample_h1 ( n, h, seed, point )

!*****************************************************************************80
!
!! P03_SAMPLE_H1 samples points from the enlarged region in problem 03.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) H, the enlargement amount.
!
!    Input/output, integer SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = rk ) POINT(3,N), the coordinates
!    of the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer n

  real ( kind = rk ) h
  integer i
  real ( kind = rk ) point(dim_num,n)
  integer seed

  call random_number ( harvest = point(1:dim_num,1:n) )

  do i = 1, 3

    point(i,1:n) = ( 1.0D+00 - point(i,1:n) ) * (         - h ) &
                 +             point(i,1:n)   * ( 1.0D+00 + h )

  end do

  return
end
subroutine p03_sdist ( n, point, sdist )

!*****************************************************************************80
!
!! P03_SDIST returns the signed distance to the region in problem 03.
!
!  Discussion:
!
!    A positive distance indicates the point is outside the region.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) POINT(3,N), the coordinates
!    of the points.
!
!    Output, real ( kind = rk ) SDIST(N), the signed distance of
!    each point to the region.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer n

  integer j
  real ( kind = rk ) point(dim_num,n)
  real ( kind = rk ) sdist(n)

  do j = 1, n
    sdist(j) = min ( -point(1,j), 1.0D+00 - point(1,j), &
                     -point(2,j), 1.0D+00 - point(2,j), &
                     -point(3,j), 1.0D+00 - point(3,j) )

  end do

  return
end
subroutine p03_title ( title )

!*****************************************************************************80
!
!! P03_TITLE returns a title for problem 03.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = * ) title

  title = '#3: The unit cube.'

  return
end
subroutine p04_boundary_edge_num ( boundary_edge_num )

!*****************************************************************************80
!
!! P04_BOUNDARY_EDGE_NUM counts the boundary edges in problem 04.
!
!  Discussion:
!
!    Call this routine if you need to allocate memory for an
!    array in which to store the boundary edges.
!
!    We assume that the boundary is triangulated.  BOUNDARY_EDGE_NUM
!    is the number of edges associated with this triangulation.
!    Thus, for a cube, the value is NOT the expected 12, but rather
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer BOUNDARY_EDGE_NUM, the number of 
!    boundary edges.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer boundary_edge_num

  boundary_edge_num = 0

  return
end
subroutine p04_boundary_edges ( boundary_edge_num, boundary_edge )

!*****************************************************************************80
!
!! P04_BOUNDARY_EDGES returns the boundary edges in problem 04.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer BOUNDARY_EDGE_NUM, the number of 
!    boundary edges.
!
!    Output, integer BOUNDARY_EDGE(2,BOUNDARY_EDGE_NUM),
!    the boundary edges, described as pairs of node indices.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer boundary_edge_num

  integer boundary_edge(2,boundary_edge_num)

  return
end
subroutine p04_boundary_face_num ( boundary_face_num )

!*****************************************************************************80
!
!! P04_BOUNDARY_FACE_NUM counts the boundary faces in problem 04.
!
!  Discussion:
!
!    Call this routine if you need to allocate memory for an
!    array in which to store the boundary faces.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer BOUNDARY_FACE_NUM, the number of 
!    boundary faces.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer boundary_face_num

  boundary_face_num = 0

  return
end
subroutine p04_boundary_faces ( boundary_face_num, boundary_face )

!*****************************************************************************80
!
!! P04_BOUNDARY_FACES returns the boundary faces in problem 04.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer BOUNDARY_FACE_NUM, the number of 
!    boundary faces.
!
!    Output, integer BOUNDARY_FACE(3,BOUNDARY_FACE_NUM),
!    the boundary faces, described as triples of node indices.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer boundary_face_num

  integer boundary_face(3,boundary_face_num)

  return
end
subroutine p04_boundary_node_num ( boundary_node_num )

!*****************************************************************************80
!
!! P04_BOUNDARY_NODE_NUM counts the boundary nodes in problem 04.
!
!  Discussion:
!
!    Call this routine if you need to allocate memory for an
!    array in which to store the boundary nodes.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer BOUNDARY_NODE_NUM, the number of 
!    boundary nodes.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer boundary_node_num

  boundary_node_num = 0

  return
end
subroutine p04_boundary_nodes ( boundary_node_num, boundary_node )

!*****************************************************************************80
!
!! P04_BOUNDARY_NODES returns the boundary nodes in problem 04.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer BOUNDARY_NODE_NUM, the number of 
!    boundary nodes.
!
!    Output, real ( kind = rk ) BOUNDARY_NODE(3,BOUNDARY_NODE_NUM),
!    the boundary nodes.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer boundary_node_num
  integer, parameter :: dim_num = 3

  real ( kind = rk ) boundary_node(dim_num,boundary_node_num)

  return
end
subroutine p04_boundary_project ( n, point )

!*****************************************************************************80
!
!! P04_BOUNDARY_PROJECT projects exterior points to the boundary in problem 04.
!
!  Discussion:
!
!    Points in the interior are not changed.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!
!    Input/output, real ( kind = rk ) POINT(3,N), the coordinates
!    of the points.  On output, all exterior points have been
!    replaced by the nearest point on the boundary.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer n

  integer j
  real ( kind = rk ) norm
  real ( kind = rk ), dimension ( dim_num, n ) :: point

  do j = 1, n

    norm = sqrt ( sum ( point(1:3,j)**2 ) )

    if ( norm == 0.0D+00 ) then
      point(1:3,j) = 1.0D+00 / sqrt ( 3.0D+00 )
    else
      point(1:3,j) = point(1:3,j) / norm
    end if

  end do

  return
end
subroutine p04_box ( lo, hi )

!*****************************************************************************80
!
!! P04_BOX returns a bounding box for problem 04.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) LO(3), HI(3), coordinates of the
!    low and high corners of the box.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) hi(dim_num)
  real ( kind = rk ) lo(dim_num)

  lo(1:dim_num) = (/ -1.0D+00, -1.0D+00, -1.0D+00 /)
  hi(1:dim_num) = (/  1.0D+00,  1.0D+00,  1.0D+00 /)

  return
end
subroutine p04_fixed_num ( fixed_num )

!*****************************************************************************80
!
!! P04_FIXED_NUM returns the number of fixed points in problem 04.
!
!  Discussion:
!
!    Call this routine if you need to allocate memory for an
!    array in which to store the fixed points.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer FIXED_NUM, the number of fixed points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer fixed_num

  fixed_num = 0

  return
end
subroutine p04_fixed_points ( fixed_num, fixed )

!*****************************************************************************80
!
!! P04_FIXED_POINTS returns the fixed points in problem 04.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer FIXED_NUM, the number of fixed points.
!
!    Output, real ( kind = rk ) FIXED(3,FIXED_NUM), the coordinates
!    of the fixed points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer fixed_num

  real ( kind = rk ) fixed(dim_num,fixed_num)

  return
end
subroutine p04_header ( )

!*****************************************************************************80
!
!! P04_HEADER prints some information about problem 04.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'P04:'
  write ( *, '(a)' ) '  Sample problem #4'
  write ( *, '(a)' ) '  The unit sphere.'

  return
end
subroutine p04_inside ( n, point, inside )

!*****************************************************************************80
!
!! P04_INSIDE reports if a point is inside the region in problem 04.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) POINT(3,N), the coordinates
!    of the points.
!
!    Output, logical INSIDE(N), is TRUE if the point is in the region.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer n

  logical inside(n)
  integer j
  real ( kind = rk ) point(dim_num,n)

  do j = 1, n
    inside(j) = sum ( point(1:3,j)**2 ) < 1.0D+00
  end do

  return
end
subroutine p04_sample ( n, seed, point )

!*****************************************************************************80
!
!! P04_SAMPLE samples points from the region in problem 04.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!
!    Input/output, integer SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = rk ) POINT(3,N), the coordinates
!    of the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer n

  real ( kind = rk ) point(dim_num,n)
  integer seed

  call ball_unit_sample_3d ( seed, point )

  return
end
subroutine p04_sample_h1 ( n, h, seed, point )

!*****************************************************************************80
!
!! P04_SAMPLE_H1 samples points from the enlarged region in problem 04.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) H, the enlargement amount.
!
!    Input/output, integer SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = rk ) POINT(3,N), the coordinates
!    of the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer n

  real ( kind = rk ) h
  real ( kind = rk ) point(dim_num,n)
  integer seed

  call ball_unit_sample_3d ( seed, point )

  point(1:3,1:n) = ( 1.0D+00 + h ) * point(1:3,1:n)

  return
end
subroutine p04_sdist ( n, point, sdist )

!*****************************************************************************80
!
!! P04_SDIST returns the signed distance to the region in problem 04.
!
!  Discussion:
!
!    A positive distance indicates the point is outside the region.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) POINT(3,N), the coordinates
!    of the points.
!
!    Output, real ( kind = rk ) SDIST(N), the signed distance of
!    each point to the region.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer n

  integer j
  real ( kind = rk ) point(dim_num,n)
  real ( kind = rk ) sdist(n)

  do j = 1, n
    sdist(j) = sqrt ( sum ( point(1:3,j)**2 ) ) - 1.0D+00
  end do

  return
end
subroutine p04_title ( title )

!*****************************************************************************80
!
!! P04_TITLE returns a title for problem 04.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, a title for the problem.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = * ) title

  title = '#4: The unit sphere.'

  return
end
subroutine r8mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
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
!    Input, character ( len = * ) TITLE, an optional title.
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
!    Input, integer ILO, JLO, the first row and column to print.
!
!    Input, integer IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, an optional title.
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

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(i7,7x)') i
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

      write ( *, '(i5,1x,5a14)' ) j, ( ctemp(i), i = 1, inc )

    end do

  end do

  write ( *, '(a)' ) ' '

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

  integer i
  integer j
  integer k
  real ( kind = rk ) r8vec_length
  real ( kind = rk ) v1(dim_num)
  real ( kind = rk ) v2(dim_num)
  real ( kind = rk ) vj
  real ( kind = rk ) vk

  if ( dim_num < 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_ANY_NORMAL - Fatal error!'
    write ( *, '(a)' ) '  Called with DIM_NUM < 2.'
    stop
  end if

  if ( r8vec_length ( dim_num, v1 ) == 0.0D+00 ) then
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
function r8vec_length ( dim_num, x )

!*****************************************************************************80
!
!! R8VEC_LENGTH returns the Euclidean length of an R8VEC.
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
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = rk ) X(DIM_NUM), the vector.
!
!    Output, real ( kind = rk ) R8VEC_LENGTH, the Euclidean length of the vector.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num

  real ( kind = rk ) r8vec_length
  real ( kind = rk ) x(dim_num)

  r8vec_length = sqrt ( sum ( ( x(1:dim_num) )**2 ) )

  return
end
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, L E Schrage,
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
!    P A Lewis, A S Goodman, J M Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer M, the number of entries in the vector.
!
!    Input/output, integer SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = rk ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer k
  integer seed
  real ( kind = rk ) r(n)

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + 2147483647
    end if

    r(i) = real ( seed, kind = rk ) * 4.656612875D-10

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
