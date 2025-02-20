function arc_cosine ( c )

!*****************************************************************************80
!
!! arc_cosine() computes the arc cosine function, with argument truncation.
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
  c2 = max ( c2, - 1.0D+00 )
  c2 = min ( c2, + 1.0D+00 )

  arc_cosine = acos ( c2 )

  return
end
function arc_sine ( s )

!*****************************************************************************80
!
!! ARC_SINE computes the arc sine function, with argument truncation.
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
!    Output, real ( kind = rk ) ARC_SINE, an angle whose sine is S.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) arc_sine
  real ( kind = rk ) s
  real ( kind = rk ) s2

  s2 = s
  s2 = max ( s2, - 1.0D+00 )
  s2 = min ( s2, + 1.0D+00 )

  arc_sine = asin ( s2 )

  return
end
function atan4 ( y, x )

!*****************************************************************************80
!
!! ATAN4 computes the inverse tangent of the ratio Y / X.
!
!  Discussion:
!
!    ATAN4 returns an angle whose tangent is ( Y / X ), a job which
!    the built in functions ATAN and ATAN2 already do.
!
!    However:
!
!    * ATAN4 always returns a positive angle, between 0 and 2 PI,
!      while ATAN and ATAN2 return angles in the interval [-PI/2,+PI/2]
!      and [-PI,+PI] respectively;
!
!    * ATAN4 accounts for the signs of X and Y, (as does ATAN2).  The ATAN
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
!    Output, real ( kind = rk ) ATAN4, an angle between 0 and 2 * PI, 
!    whose tangent is (Y/X), and which lies in the appropriate quadrant so 
!    that the signs of its cosine and sine match those of X and Y.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) abs_x
  real ( kind = rk ) abs_y
  real ( kind = rk ) atan4
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) theta
  real ( kind = rk ) theta_0
  real ( kind = rk ) x
  real ( kind = rk ) y
!
!  Special cases:
!
  if ( x == 0.0D+00 ) then

    if ( 0.0D+00 < y ) then
      theta = pi / 2.0D+00
    else if ( y < 0.0D+00 ) then
      theta = 3.0D+00 * pi / 2.0D+00
    else if ( y == 0.0D+00 ) then
      theta = 0.0D+00
    end if

  else if ( y == 0.0D+00 ) then

    if ( 0.0D+00 < x ) then
      theta = 0.0D+00
    else if ( x < 0.0D+00 ) then
      theta = PI
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
      theta = pi - theta_0
    else if ( x < 0.0D+00 .and. y < 0.0D+00 ) then
      theta = pi + theta_0
    else if ( 0.0D+00 < x .and. y < 0.0D+00 ) then
      theta = 2.0D+00 * pi - theta_0
    end if

  end if

  atan4 = theta

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
!    The input data required for this routine can be retrieved from ICOS_SIZE.
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
  integer point_num

  real ( kind = rk ) a
  real ( kind = rk ) b
  integer edge_point(edge_order,edge_num)
  integer face_order(face_num)
  integer face_point(face_order_max,face_num)
  real ( kind = rk ) phi
  real ( kind = rk ) point_coord(3,point_num)
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
  point_coord(1:3,1:point_num) = reshape ( (/ &
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
     -a, -b,  z /), (/ 3, point_num /) )
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
!! ICOS_SIZE gives "sizes" for an icosahedron in 3D.
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
subroutine s_cat ( s1, s2, s3 )

!*****************************************************************************80
!
!! S_CAT concatenates two strings to make a third string.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S1, the "prefix" string.
!
!    Input, character ( len = * ) S2, the "postfix" string.
!
!    Output, character ( len = * ) S3, the string made by
!    concatenating S1 and S2, ignoring any trailing blanks.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = * ) s1
  character ( len = * ) s2
  character ( len = * ) s3

  if ( s1 == ' ' .and. s2 == ' ' ) then
    s3 = ' '
  else if ( s1 == ' ' ) then
    s3 = s2
  else if ( s2 == ' ' ) then
    s3 = s1
  else
    s3 = trim ( s1 ) // trim ( s2 )
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

  real ( kind = rk ) arc_sine
  real ( kind = rk ) atan4
  real ( kind = rk ) bot
  real ( kind = rk ) dist
  real ( kind = rk ) lat1
  real ( kind = rk ) lat2
  real ( kind = rk ) lon1
  real ( kind = rk ) lon2
  real ( kind = rk ) top
  real ( kind = rk ) xyz1(3)
  real ( kind = rk ) xyz2(3)

  lat1 = arc_sine ( xyz1(3) )
  lon1 = atan4 ( xyz1(2), xyz1(1) )

  lat2 = arc_sine ( xyz2(3) )
  lon2 = atan4 ( xyz2(2), xyz2(1) )

  top = ( cos ( lat2 ) * sin ( lon1 - lon2 ) )**2 &
      + ( cos ( lat1 ) * sin ( lat2 ) &
      -   sin ( lat1 ) * cos ( lat2 ) * cos ( lon1 - lon2 ) )**2

  top = sqrt ( top )

  bot = sin ( lat1 ) * sin ( lat2 ) &
      + cos ( lat1 ) * cos ( lat2 ) * cos ( lon1 - lon2 )

  dist = atan2 ( top, bot )

  return
end
subroutine sphere01_monomial_integral ( e, integral )

!*****************************************************************************80
!
!! SPHERE01_MONOMIAL_INTEGRAL returns monomial integrals on the unit sphere.
!
!  Discussion:
!
!    The integration region is 
!
!      X^2 + Y^2 + Z^2 = 1.
!
!    The monomial is F(X,Y,Z) = X^E(1) * Y^E(2) * Z^E(3).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    24 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Academic Press, 1984, page 263.
!
!  Parameters:
!
!    Input, integer E(3), the exponents of X, Y and Z in the 
!    monomial.  Each exponent must be nonnegative.
!
!    Output, real ( kind = rk ) INTEGRAL, the integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer e(3)
  integer i
  real ( kind = rk ) integral
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00

  if ( any ( e(1:3) < 0 ) ) then
    integral = - huge ( integral )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPHERE01_MONOMIAL_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  All exponents must be nonnegative.'
    write ( *, '(a,i8)' ) '  E(1) = ', e(1)
    write ( *, '(a,i8)' ) '  E(2) = ', e(2)
    write ( *, '(a,i8)' ) '  E(3) = ', e(3)
    stop
  end if

  if ( all ( e(1:3) == 0 ) ) then

    integral = 2.0D+00 * sqrt ( pi**3 ) / gamma ( 1.5D+00 )

  else if ( any ( mod ( e(1:3), 2 ) == 1 ) ) then

    integral = 0.0D+00

  else

    integral = 2.0D+00

    do i = 1, 3
      integral = integral * gamma ( 0.5D+00 * real ( e(i) + 1, kind = rk ) )
    end do

    integral = integral &
      / gamma ( 0.5D+00 * ( real ( sum ( e(1:3) + 1 ), kind = rk ) ) )

  end if

  return
end
subroutine sphere01_quad_icos1c ( factor, fun, node_num, result )

!*****************************************************************************80
!
!! SPHERE01_QUAD_ICOS1C: centroid rule, subdivide then project.
!
!  Discussion:
!
!    This function estimates an integral over the surface of the unit sphere.
!
!    This function sets up an icosahedral grid, and subdivides each
!    edge of the icosahedron into FACTOR subedges.  These edges define a grid
!    within each triangular icosahedral face.  The centroids of these
!    triangles can be determined.  All of these calculations are done,
!    essentially, on the FLAT faces of the icosahedron.  Only then are
!    the triangle vertices and centroids projected to the sphere.  
!
!    The resulting grid of spherical triangles and projected centroids
!    is used to apply a centroid quadrature rule over the surface of
!    the unit sphere.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer FACTOR, the subdivision factor, which must
!    be at least 1.
!
!    Input, external :: FUN, evaluates the integrand, of the form:
!      subroutine fun ( n, x, v )
!      integer n
!      real ( kind = rk ) v(n)
!      real ( kind = rk ) x(3,n)
!
!    Output, integer NODE_NUM, the number of evaluation points.
!
!    Output, real ( kind = rk ) RESULT, the estimated integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer node_num

  integer a
  real ( kind = rk ) a_xyz(3)
  real ( kind = rk ) a2_xyz(3)
  real ( kind = rk ) area
  real ( kind = rk ) area_total
  integer b
  real ( kind = rk ) b_xyz(3)
  real ( kind = rk ) b2_xyz(3)
  integer c
  real ( kind = rk ) c_xyz(3)
  real ( kind = rk ) c2_xyz(3)
  integer edge_num
  integer, allocatable, dimension ( :, : ) :: edge_point
  integer f1
  integer f2
  integer f3
  integer face
  integer face_num
  integer, allocatable, dimension ( : ) :: face_order
  integer, allocatable, dimension ( :, : ) :: face_point
  integer face_order_max
  integer factor
  external             fun
  real ( kind = rk ) node_xyz(3)
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ), allocatable, dimension ( :, : ) :: point_coord
  integer point_num
  real ( kind = rk ) result
  real ( kind = rk ) v
!
!  Size the icosahedron.
!
  call icos_size ( point_num, edge_num, face_num, face_order_max )
!
!  Set the icosahedron.
!
  allocate ( point_coord(1:3,1:point_num) )
  allocate ( edge_point(1:2,1:edge_num) )
  allocate ( face_order(1:face_num) )
  allocate ( face_point(1:face_order_max,1:face_num) )

  call icos_shape ( point_num, edge_num, face_num, face_order_max, &
    point_coord, edge_point, face_order, face_point )
!
!  Initialize the integral data.
!
  result = 0.0D+00
  area_total = 0.0D+00
  node_num = 0
!
!  Pick a face of the icosahedron, and identify its vertices as A, B, C.
!
  do face = 1, face_num

    a = face_point(1,face)
    b = face_point(2,face)
    c = face_point(3,face)

    a_xyz(1:3) = point_coord(1:3,a)
    b_xyz(1:3) = point_coord(1:3,b)
    c_xyz(1:3) = point_coord(1:3,c)
!
!  Some subtriangles will have the same direction as the face.
!  Generate each in turn, by determining the barycentric coordinates
!  of the centroid (F1,F2,F3), from which we can also work out the barycentric
!  coordinates of the vertices of the subtriangle.
!
    do f3 = 1, 3 * factor - 2, 3
      do f2 = 1, 3 * factor - f3 - 1, 3

        f1 = 3 * factor - f3 - f2

        call sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1, f2, f3, &
          node_xyz )

        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1 + 2, f2 - 1, f3 - 1, a2_xyz )
        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1 - 1, f2 + 2, f3 - 1, b2_xyz )
        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1 - 1, f2 - 1, f3 + 2, c2_xyz )

        call sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz, area )

        call fun ( 1, node_xyz, v )    

        node_num = node_num + 1
        result = result + area * v
        area_total = area_total + area

      end do
    end do
!
!  The other subtriangles have the opposite direction from the face.
!  Generate each in turn, by determining the barycentric coordinates
!  of the centroid (F1,F2,F3), from which we can also work out the barycentric
!  coordinates of the vertices of the subtriangle.
!
    do f3 = 2, 3 * factor - 4, 3
      do f2 = 2, 3 * factor - f3 - 2, 3

        f1 = 3 * factor - f3 - f2

        call sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1, f2, f3, &
          node_xyz )

        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1 - 2, f2 + 1, f3 + 1, a2_xyz )
        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1 + 1, f2 - 2, f3 + 1, b2_xyz )
        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1 + 1, f2 + 1, f3 - 2, c2_xyz )

        call sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz, area )

        call fun ( 1, node_xyz, v )  

        node_num = node_num + 1  
        result = result + area * v
        area_total = area_total + area

      end do
    end do

  end do
!
!  Discard allocated memory.
!
  deallocate ( edge_point )
  deallocate ( face_order )
  deallocate ( face_point )
  deallocate ( point_coord )

  return
end
subroutine sphere01_quad_icos1m ( factor, fun, node_num, result )

!*****************************************************************************80
!
!! SPHERE01_QUAD_ICOS1M: midside rule, subdivide then project.
!
!  Discussion:
!
!    This function estimates an integral over the surface of the unit sphere.
!
!    This function sets up an icosahedral grid, and subdivides each
!    edge of the icosahedron into FACTOR subedges.  These edges define a grid
!    within each triangular icosahedral face.  The midsides of these
!    triangles can be determined.  All of these calculations are done,
!    essentially, on the FLAT faces of the icosahedron.  Only then are
!    the triangle vertices and midsides projected to the sphere.  
!
!    The resulting grid of spherical triangles and projected midsides
!    is used to apply a midside quadrature rule over the surface of
!    the unit sphere.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer FACTOR, the subdivision factor, which must
!    be at least 1.
!
!    Input, external :: FUN, evaluates the integrand, of the form:
!      subroutine fun ( n, x, v )
!      integer n
!      real ( kind = rk ) v(n)
!      real ( kind = rk ) x(3,n)
!
!    Output, integer NODE_NUM, the number of evaluation points.
!
!    Output, real ( kind = rk ) RESULT, the estimated integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer node_num

  integer a
  real ( kind = rk ) a_xyz(3)
  real ( kind = rk ) a2_xyz(3)
  real ( kind = rk ) a3_xyz(3)
  real ( kind = rk ) area
  real ( kind = rk ) area_total
  integer b
  real ( kind = rk ) b_xyz(3)
  real ( kind = rk ) b2_xyz(3)
  real ( kind = rk ) b3_xyz(3)
  integer c
  real ( kind = rk ) c_xyz(3)
  real ( kind = rk ) c2_xyz(3)
  real ( kind = rk ) c3_xyz(3)
  integer edge_num
  integer, allocatable, dimension ( :, : ) :: edge_point
  integer f1
  integer f2
  integer f3
  integer face
  integer face_num
  integer, allocatable, dimension ( : ) :: face_order
  integer, allocatable, dimension ( :, : ) :: face_point
  integer face_order_max
  integer factor
  external fun
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ), allocatable, dimension ( :, : ) :: point_coord
  integer point_num
  real ( kind = rk ) result
  real ( kind = rk ) va
  real ( kind = rk ) vb
  real ( kind = rk ) vc
!
!  Size the icosahedron.
!
  call icos_size ( point_num, edge_num, face_num, face_order_max )
!
!  Set the icosahedron.
!
  allocate ( point_coord(1:3,1:point_num) )
  allocate ( edge_point(1:2,1:edge_num) )
  allocate ( face_order(1:face_num) )
  allocate ( face_point(1:face_order_max,1:face_num) )

  call icos_shape ( point_num, edge_num, face_num, face_order_max, &
    point_coord, edge_point, face_order, face_point )
!
!  Initialize the integral data.
!
  result = 0.0D+00
  node_num = 0
  area_total = 0.0D+00
!
!  Consider each face.
!
  do face = 1, face_num

    a = face_point(1,face)
    b = face_point(2,face)
    c = face_point(3,face)

    a_xyz(1:3) = point_coord(1:3,a)
    b_xyz(1:3) = point_coord(1:3,b)
    c_xyz(1:3) = point_coord(1:3,c)
!
!  Deal with subtriangles that have same orientation as face.
!
    do f1 = 0, factor - 1
      do f2 = 0, factor - f1 - 1
        f3 = factor - f1 - f2

        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1 + 1, f2,     f3 - 1, a2_xyz )
        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1,     f2 + 1, f3 - 1, b2_xyz )
        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1,     f2,     f3,     c2_xyz )

        call sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz, area )

        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, 2 * f1 + 1, 2 * f2 + 1, 2 * f3 - 2, a3_xyz )
        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, 2 * f1,     2 * f2 + 1, 2 * f3 - 1, b3_xyz )
        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, 2 * f1 + 1, 2 * f2,     2 * f3 - 1, c3_xyz )

        node_num = node_num + 3
        call fun ( 1, a3_xyz, va )
        call fun ( 1, b3_xyz, vb )   
        call fun ( 1, c3_xyz, vc )   
        result = result + area * ( va + vb + vc ) / 3.0D+00
        area_total = area_total + area

      end do
    end do
!
!  Deal with subtriangles that have opposite orientation as face.
!
    do f3 = 0, factor - 2
      do f2 = 1, factor - f3 - 1
        f1 = factor - f2 - f3

        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1 - 1, f2,     f3 + 1, a2_xyz )
        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1,     f2 - 1, f3 + 1, b2_xyz )
        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1,     f2,     f3,     c2_xyz )

        call sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz, area )

        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, 2 * f1 - 1, 2 * f2 - 1, 2 * f3 + 2, a3_xyz )
        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, 2 * f1,     2 * f2 - 1, 2 * f3 + 1, b3_xyz )
        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, 2 * f1 - 1, 2 * f2,     2 * f3 + 1, c3_xyz )

        node_num = node_num + 3
        call fun ( 1, a3_xyz, va )
        call fun ( 1, b3_xyz, vb )   
        call fun ( 1, c3_xyz, vc )   
        result = result + area * ( va + vb + vc ) / 3.0D+00
        area_total = area_total + area

      end do
    end do
  end do
!
!  Discard allocated memory.
!
  deallocate ( edge_point )
  deallocate ( face_order )
  deallocate ( face_point )
  deallocate ( point_coord )

  return
end
subroutine sphere01_quad_icos1v ( factor, fun, node_num, result )

!*****************************************************************************80
!
!! SPHERE01_QUAD_ICOS1V: vertex rule, subdivide then project.
!
!  Discussion:
!
!    This function estimates an integral over the surface of the unit sphere.
!
!    This function sets up an icosahedral grid, and subdivides each
!    edge of the icosahedron into FACTOR subedges.  These edges define a grid
!    within each triangular icosahedral face.  The vertices of these
!    triangles can be determined.  All of these calculations are done,
!    essentially, on the FLAT faces of the icosahedron.  Only then are
!    the triangle vertices projected to the sphere.  
!
!    The resulting grid of spherical triangles is used to apply a vertex
!    quadrature rule over the surface of the unit sphere.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer FACTOR, the subdivision factor, which must
!    be at least 1.
!
!    Input, external :: FUN, evaluates the integrand, of the form:
!      subroutine fun ( n, x, v )
!      integer n
!      real ( kind = rk ) v(n)
!      real ( kind = rk ) x(3,n)
!
!    Output, integer NODE_NUM, the number of evaluation points.
!
!    Output, real ( kind = rk ) RESULT, the estimated integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer node_num

  integer a
  real ( kind = rk ) a_xyz(3)
  real ( kind = rk ) a2_xyz(3)
  real ( kind = rk ) area
  real ( kind = rk ) area_total
  integer b
  real ( kind = rk ) b_xyz(3)
  real ( kind = rk ) b2_xyz(3)
  integer c
  real ( kind = rk ) c_xyz(3)
  real ( kind = rk ) c2_xyz(3)
  integer edge_num
  integer, allocatable, dimension ( :, : ) :: edge_point
  integer f1
  integer f2
  integer f3
  integer face
  integer face_num
  integer, allocatable, dimension ( : ) :: face_order
  integer, allocatable, dimension ( :, : ) :: face_point
  integer face_order_max
  integer factor
  external fun
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ), allocatable, dimension ( :, : ) :: point_coord
  integer point_num
  real ( kind = rk ) result
  real ( kind = rk ) va
  real ( kind = rk ) vb
  real ( kind = rk ) vc
!
!  Size the icosahedron.
!
  call icos_size ( point_num, edge_num, face_num, face_order_max )
!
!  Set the icosahedron.
!
  allocate ( point_coord(1:3,1:point_num) )
  allocate ( edge_point(1:2,1:edge_num) )
  allocate ( face_order(1:face_num) )
  allocate ( face_point(1:face_order_max,1:face_num) )

  call icos_shape ( point_num, edge_num, face_num, face_order_max, &
    point_coord, edge_point, face_order, face_point )
!
!  Initialize the integral data.
!
  result = 0.0D+00
  node_num = 0
  area_total = 0.0D+00
!
!  Consider each face.
!
  do face = 1, face_num

    a = face_point(1,face)
    b = face_point(2,face)
    c = face_point(3,face)

    a_xyz(1:3) = point_coord(1:3,a)
    b_xyz(1:3) = point_coord(1:3,b)
    c_xyz(1:3) = point_coord(1:3,c)
!
!  Deal with subtriangles that have same orientation as face.
!
    do f1 = 0, factor - 1
      do f2 = 0, factor - f1 - 1
        f3 = factor - f1 - f2

        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1 + 1, f2,     f3 - 1, a2_xyz )
        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1,     f2 + 1, f3 - 1, b2_xyz )
        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1,     f2,     f3,     c2_xyz )

        call sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz, area )

        node_num = node_num + 3
        call fun ( 1, a2_xyz, va )
        call fun ( 1, b2_xyz, vb )   
        call fun ( 1, c2_xyz, vc )   
        result = result + area * ( va + vb + vc ) / 3.0D+00
        area_total = area_total + area

      end do
    end do
!
!  Deal with subtriangles that have opposite orientation as face.
!
    do f3 = 0, factor - 2
      do f2 = 1, factor - f3 - 1
        f1 = factor - f2 - f3

        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1 - 1, f2,     f3 + 1, a2_xyz )
        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1,     f2 - 1, f3 + 1, b2_xyz )
        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1,     f2,     f3,     c2_xyz )

        call sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz, area )

        node_num = node_num + 3
        call fun ( 1, a2_xyz, va )
        call fun ( 1, b2_xyz, vb )   
        call fun ( 1, c2_xyz, vc )   
        result = result + area * ( va + vb + vc ) / 3.0D+00
        area_total = area_total + area

      end do
    end do
  end do
!
!  Discard allocated memory.
!
  deallocate ( edge_point )
  deallocate ( face_order )
  deallocate ( face_point )
  deallocate ( point_coord )

  return
end
subroutine sphere01_quad_icos2v ( factor, fun, node_num, result )

!*****************************************************************************80
!
!! SPHERE01_QUAD_ICOS2V: vertex rule, subdivide then project.
!
!  Discussion:
!
!    This function estimates an integral over the surface of the unit sphere.
!
!    This function sets up an icosahedral grid, and subdivides each
!    edge of the icosahedron into FACTOR subedges.  These edges define a grid
!    within each triangular icosahedral face.  The vertices of these
!    triangles can be determined.  All of these calculations are done,
!    essentially, on the FLAT faces of the icosahedron.  Only then are
!    the triangle vertices projected to the sphere.  
!
!    The resulting grid of spherical triangles is used to apply a vertex
!    quadrature rule over the surface of the unit sphere.
!
!    This is a revision of SPHERE01_QUAD_ICOS2V that attempted to use a more
!    sophisticated scheme to map points from the planar triangle to the surface
!    of the unit sphere.  Very little improvement to the estimated integral
!    was observed, so development of this scheme has been set aside for now.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer FACTOR, the subdivision factor, which must
!    be at least 1.
!
!    Input, external :: FUN, evaluates the integrand, of the form:
!      subroutine fun ( n, x, v )
!      integer n
!      real ( kind = rk ) v(n)
!      real ( kind = rk ) x(3,n)
!
!    Output, integer NODE_NUM, the number of evaluation points.
!
!    Output, real ( kind = rk ) RESULT, the estimated integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer node_num

  integer a
  real ( kind = rk ) a_xyz(3)
  real ( kind = rk ) a2_xyz(3)
  real ( kind = rk ) area
  real ( kind = rk ) area_total
  integer b
  real ( kind = rk ) b_xyz(3)
  real ( kind = rk ) b2_xyz(3)
  integer c
  real ( kind = rk ) c_xyz(3)
  real ( kind = rk ) c2_xyz(3)
  integer edge_num
  integer, allocatable, dimension ( :, : ) :: edge_point
  integer f1
  integer f2
  integer f3
  integer face
  integer face_num
  integer, allocatable, dimension ( : ) :: face_order
  integer, allocatable, dimension ( :, : ) :: face_point
  integer face_order_max
  integer factor
  external fun
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ), allocatable, dimension ( :, : ) :: point_coord
  integer point_num
  real ( kind = rk ) result
  real ( kind = rk ) va
  real ( kind = rk ) vb
  real ( kind = rk ) vc
!
!  Size the icosahedron.
!
  call icos_size ( point_num, edge_num, face_num, face_order_max )
!
!  Set the icosahedron.
!
  allocate ( point_coord(1:3,1:point_num) )
  allocate ( edge_point(1:2,1:edge_num) )
  allocate ( face_order(1:face_num) )
  allocate ( face_point(1:face_order_max,1:face_num) )

  call icos_shape ( point_num, edge_num, face_num, face_order_max, &
    point_coord, edge_point, face_order, face_point )
!
!  Initialize the integral data.
!
  result = 0.0D+00
  node_num = 0
  area_total = 0.0D+00
!
!  Consider each face.
!
  do face = 1, face_num

    a = face_point(1,face)
    b = face_point(2,face)
    c = face_point(3,face)

    a_xyz(1:3) = point_coord(1:3,a)
    b_xyz(1:3) = point_coord(1:3,b)
    c_xyz(1:3) = point_coord(1:3,c)
!
!  Deal with subtriangles that have same orientation as face.
!
    do f1 = 0, factor - 1
      do f2 = 0, factor - f1 - 1
        f3 = factor - f1 - f2

        call sphere01_triangle_project2 ( &
          a_xyz, b_xyz, c_xyz, f1 + 1, f2,     f3 - 1, a2_xyz )
        call sphere01_triangle_project2 ( &
          a_xyz, b_xyz, c_xyz, f1,     f2 + 1, f3 - 1, b2_xyz )
        call sphere01_triangle_project2 ( &
          a_xyz, b_xyz, c_xyz, f1,     f2,     f3,     c2_xyz )

        call sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz, area )

        node_num = node_num + 3
        call fun ( 1, a2_xyz, va )
        call fun ( 1, b2_xyz, vb )   
        call fun ( 1, c2_xyz, vc )   
        result = result + area * ( va + vb + vc ) / 3.0D+00
        area_total = area_total + area

      end do
    end do
!
!  Deal with subtriangles that have opposite orientation as face.
!
    do f3 = 0, factor - 2
      do f2 = 1, factor - f3 - 1
        f1 = factor - f2 - f3

        call sphere01_triangle_project2 ( &
          a_xyz, b_xyz, c_xyz, f1 - 1, f2,     f3 + 1, a2_xyz )
        call sphere01_triangle_project2 ( &
          a_xyz, b_xyz, c_xyz, f1,     f2 - 1, f3 + 1, b2_xyz )
        call sphere01_triangle_project2 ( &
          a_xyz, b_xyz, c_xyz, f1,     f2,     f3,     c2_xyz )

        call sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz, area )

        node_num = node_num + 3
        call fun ( 1, a2_xyz, va )
        call fun ( 1, b2_xyz, vb )   
        call fun ( 1, c2_xyz, vc )   
        result = result + area * ( va + vb + vc ) / 3.0D+00
        area_total = area_total + area

      end do
    end do
  end do
!
!  Discard allocated memory.
!
  deallocate ( edge_point )
  deallocate ( face_order )
  deallocate ( face_point )
  deallocate ( point_coord )

  return
end
subroutine sphere01_quad_llc ( f, h, n, result )

!*****************************************************************************80
!
!! SPHERE01_QUAD_LLC: Longitude/Latitude grid with centroid rule.
!
!  Discussion:
!
!    The sphere is broken up into spherical triangles, whose sides
!    do not exceed the length H.  Then a centroid rule is used on
!    each spherical triangle.
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
!    Input, external :: F, evaluates the integrand, of the form:
!      subroutine f ( n, x, v )
!      integer n
!      real ( kind = rk ) v(n)
!      real ( kind = rk ) x(3,n)
!
!    Input, real ( kind = rk ) H, the maximum length of a side of the spherical
!    quadrilaterals.
!
!    Output, integer N, the number of points used.
!
!    Output, real ( kind = rk ) RESULT, the approximate integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) area
  external f
  real ( kind = rk ) h
  integer i
  integer j
  integer n
  real ( kind = rk ) phi
  integer phi_num
  real ( kind = rk ) phi1
  real ( kind = rk ) phi2
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) result
  real ( kind = rk ) sector_area
  real ( kind = rk ) sphere_area
  real ( kind = rk ) theta
  integer theta_num
  real ( kind = rk ) theta1
  real ( kind = rk ) theta2
  real ( kind = rk ) v(1)
  real ( kind = rk ) x(3)
  real ( kind = rk ) x1(3)
  real ( kind = rk ) x11(3)
  real ( kind = rk ) x12(3)
  real ( kind = rk ) x2(3)
  real ( kind = rk ) x21(3)
  real ( kind = rk ) x22(3)
!
!  Choose PHI and THETA counts that make short sides.
!
  phi_num = int ( pi / h )

  if ( h * real ( phi_num, kind = rk ) < pi ) then
    phi_num = phi_num + 1
  end if

  theta_num = int ( 2.0D+00 * pi / h )

  if ( h * real ( theta_num, kind = rk ) < pi ) then
    theta_num = theta_num + 1
  end if

  n = 0
  result = 0.0D+00
!
!  Only one THETA (and hence, only one PHI.)
!
  if ( theta_num == 1 ) then

    sphere_area = 4.0D+00 * pi

    theta = 0.0D+00
    phi = pi / 2.0D+00
    call tp_to_xyz ( theta, phi, x )

    call f ( 1, x, v )
    n = n + 1
    result = sphere_area * v(1)
!
!  Several THETA's, one PHI.
!
  else if ( phi_num == 1 ) then

    sphere_area = 4.0D+00 * pi
    sector_area = sphere_area / real ( theta_num, kind = rk )

    result = 0.0D+00

    do j = 1, theta_num

      theta = real ( ( j - 1 ) * 2, kind = rk ) * pi &
            / real ( theta_num, kind = rk )
      phi = pi / 2.0D+00
      call tp_to_xyz ( theta, phi, x )
      call f ( 1, x, v )
      n = n + 1
      result = result + sector_area * v(1)

    end do
!
!  At least two PHI's.
!
  else

    result = 0.0D+00
!
!  Picture in top row, with V1 = north pole:
!
!        V1
!       /  \
!      /    \
!    V12----V22
!
    phi1 = 0.0D+00
    phi2 = pi / real ( phi_num, kind = rk )

    do j = 1, theta_num

      theta1 = real ( j - 1, kind = rk ) * 2.0D+00 * pi &
             / real ( theta_num, kind = rk )
      theta2 = real ( j    , kind = rk ) * 2.0D+00 * pi &
             / real ( theta_num, kind = rk )

      call tp_to_xyz ( theta1, phi1, x1 )
      call tp_to_xyz ( theta1, phi2, x12 )
      call tp_to_xyz ( theta2, phi2, x22 )

      call sphere01_triangle_vertices_to_area ( x1, x12, x22, area )
      call sphere01_triangle_vertices_to_centroid ( x1, x12, x22, x )
      call f ( 1, x, v )
      n = n + 1
      result = result + area * v(1)

    end do
!
!  Picture in all intermediate rows:
!
!    V11--V21
!     | \  |
!     |  \ |
!    V12--V22
!
    do i = 2, phi_num-1

      phi1 = real ( i - 1, kind = rk ) * pi / real ( phi_num, kind = rk )
      phi2 = real ( i,     kind = rk ) * pi / real ( phi_num, kind = rk )

      do j = 1, theta_num

        theta1 = real ( j - 1, kind = rk ) * 2.0D+00 * pi &
               / real ( theta_num, kind = rk )
        theta2 = real ( j,     kind = rk ) * 2.0D+00 * pi &
               / real ( theta_num, kind = rk )

        call tp_to_xyz ( theta1, phi1, x11 )
        call tp_to_xyz ( theta2, phi1, x21 )
        call tp_to_xyz ( theta1, phi2, x12 )
        call tp_to_xyz ( theta2, phi2, x22 )

        call sphere01_triangle_vertices_to_area ( x11, x12, x22, area )
        call sphere01_triangle_vertices_to_centroid ( x11, x12, x22, x )
        call f ( 1, x, v )
        n = n + 1
        result = result + area * v(1)

        call sphere01_triangle_vertices_to_area ( x22, x21, x11, area )
        call sphere01_triangle_vertices_to_centroid ( x22, x21, x11, x )
        call f ( 1, x, v )
        n = n + 1
        result = result + area * v(1)

      end do

    end do
!
!  Picture in last row, with V2 = south pole:
!
!    V11----V21
!      \    /
!       \  /
!        V2
!
    phi1 = real ( phi_num - 1, kind = rk ) * pi &
         / real ( phi_num, kind = rk )
    phi2 =                                  pi

    do j = 1, theta_num

      theta1 = real ( j - 1, kind = rk ) * 2.0D+00 * pi &
             / real ( theta_num, kind = rk )
      theta2 = real ( j,     kind = rk ) * 2.0D+00 * pi &
             / real ( theta_num, kind = rk )

      call tp_to_xyz ( theta1, phi1, x11 )
      call tp_to_xyz ( theta2, phi1, x21 )
      call tp_to_xyz ( theta2, phi2, x2 )

      call sphere01_triangle_vertices_to_area ( x11, x2, x21, area )
      call sphere01_triangle_vertices_to_centroid ( x11, x2, x21, x )
      call f ( 1, x, v )
      n = n + 1
      result = result + area * v(1)

    end do

  end if

  return
end
subroutine sphere01_quad_llm ( f, h, n, result )

!*****************************************************************************80
!
!! SPHERE01_QUAD_LLM: longitude/latitude grid plus midside rule.
!
!  Discussion:
!
!    The sphere is broken up into spherical triangles, whose sides
!    do not exceed the length H.  Then the function is evaluated
!    at the midsides, and the average is multiplied by the area.
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
!    Input, external :: F, evaluates the integrand, of the form:
!      subroutine f ( n, x, v )
!      integer n
!      real ( kind = rk ) v(n)
!      real ( kind = rk ) x(3,n)
!
!    Input, real ( kind = rk ) H, the maximum length of a side of the spherical
!    quadrilaterals.
!
!    Output, integer N, the number of points used.
!
!    Output, real ( kind = rk ) RESULT, the approximate integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) area
  external f
  real ( kind = rk ) h
  integer i
  integer j
  real ( kind = rk ) m1(3)
  real ( kind = rk ) m2(3)
  real ( kind = rk ) m3(3)
  integer n
  real ( kind = rk ) phi
  integer phi_num
  real ( kind = rk ) phi1
  real ( kind = rk ) phi2
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) result
  real ( kind = rk ) sector_area
  real ( kind = rk ) sphere_area
  real ( kind = rk ) theta
  integer theta_num
  real ( kind = rk ) theta1
  real ( kind = rk ) theta2
  real ( kind = rk ) v(1)
  real ( kind = rk ) x(3)
  real ( kind = rk ) x1(3)
  real ( kind = rk ) x11(3)
  real ( kind = rk ) x12(3)
  real ( kind = rk ) x2(3)
  real ( kind = rk ) x21(3)
  real ( kind = rk ) x22(3)
!
!  Choose PHI and THETA counts that make short sides.
!
  phi_num = int ( pi / h )

  if ( h * real ( phi_num, kind = rk ) < pi ) then
    phi_num = phi_num + 1
  end if

  theta_num = int ( 2.0D+00 * pi / h )

  if ( h * real ( theta_num, kind = rk ) < pi ) then
    theta_num = theta_num + 1
  end if

  n = 0
  result = 0.0D+00
!
!  Only one THETA (and hence, only one PHI.)
!
  if ( theta_num == 1 ) then

    sphere_area = 4.0D+00 * pi

    theta = 0.0D+00
    phi = pi / 2.0D+00
    call tp_to_xyz ( theta, phi, x )
    call f ( 1, x, v )
    n = n + 1
    result = sphere_area * v(1)
!
!  Several THETA's, one PHI.
!
  else if ( phi_num == 1 ) then

    sphere_area = 4.0D+00 * pi
    sector_area = sphere_area / real ( theta_num, kind = rk )

    result = 0.0D+00

    do j = 1, theta_num

      theta = real ( ( j - 1 ) * 2, kind = rk ) * pi &
            / real ( theta_num, kind = rk )
      phi = pi / 2.0D+00
      call tp_to_xyz ( theta, phi, x )
      call f ( 1, x, v )
      n = n + 1
      result = result + sector_area * v(1)

    end do
!
!  At least two PHI's.
!
  else

    result = 0.0D+00
!
!  Picture:
!
!        V1
!       /  \
!      /    \
!    V12----V22
!
    phi1 = 0.0D+00
    phi2 = pi / real ( phi_num, kind = rk )

    do j = 1, theta_num

      theta1 = real ( j - 1, kind = rk ) * 2.0D+00 * pi &
             / real ( theta_num, kind = rk )
      theta2 = real ( j,     kind = rk ) * 2.0D+00 * pi &
             / real ( theta_num, kind = rk )

      call tp_to_xyz ( theta1, phi1, x1 )
      call tp_to_xyz ( theta1, phi2, x12 )
      call tp_to_xyz ( theta2, phi2, x22 )

      call sphere01_triangle_vertices_to_area ( x1, x12, x22, area )

      call sphere01_triangle_vertices_to_midpoints ( x1, x12, x22, &
        m1, m2, m3 )

      call f ( 1, m1, v )
      n = n + 1
      result = result + area * v(1) / 3.0D+00
      call f ( 1, m2, v )
      n = n + 1
      result = result + area * v(1) / 3.0D+00
      call f ( 1, m3, v )
      n = n + 1
      result = result + area * v(1) / 3.0D+00

    end do
!
!  Picture:
!
!    V11--V21
!     | \  |
!     |  \ |
!    V12--V22
!
    do i = 2, phi_num-1

      phi1 = real ( i - 1, kind = rk ) * pi &
           / real ( phi_num, kind = rk )
      phi2 = real ( i,     kind = rk ) * pi &
           / real ( phi_num, kind = rk )

      do j = 1, theta_num

        theta1 = real ( j - 1, kind = rk ) * 2.0D+00 * pi &
               / real ( theta_num, kind = rk )
        theta2 = real ( j,     kind = rk ) * 2.0D+00 * pi &
               / real ( theta_num, kind = rk )

        call tp_to_xyz ( theta1, phi1, x11 )
        call tp_to_xyz ( theta2, phi1, x21 )
        call tp_to_xyz ( theta1, phi2, x12 )
        call tp_to_xyz ( theta2, phi2, x22 )

        call sphere01_triangle_vertices_to_area ( x11, x12, x22, area )

        call sphere01_triangle_vertices_to_midpoints ( x11, x12, x22, &
          m1, m2, m3 )

        call f ( 1, m1, v )
        n = n + 1
        result = result + area * v(1) / 3.0D+00
        call f ( 1, m2, v )
        n = n + 1
        result = result + area * v(1) / 3.0D+00
        call f ( 1, m3, v )
        n = n + 1
        result = result + area * v(1) / 3.0D+00

        call sphere01_triangle_vertices_to_area ( x22, x21, x11, area )

        call sphere01_triangle_vertices_to_midpoints ( x22, x21, x11, &
          m1, m2, m3 )

        call f ( 1, m1, v )
        n = n + 1
        result = result + area * v(1) / 3.0D+00
        call f ( 1, m2, v )
        n = n + 1
        result = result + area * v(1) / 3.0D+00
        call f ( 1, m3, v )
        n = n + 1
        result = result + area * v(1) / 3.0D+00

      end do

    end do
!
!  Picture:
!
!    V11----V21
!      \    /
!       \  /
!        V2
!
    phi1 = real ( phi_num - 1, kind = rk ) * pi &
         / real ( phi_num, kind = rk )
    phi2 =                                  pi

    do j = 1, theta_num

      theta1 = real ( j - 1, kind = rk ) * 2.0D+00 * pi &
             / real ( theta_num, kind = rk )
      theta2 = real ( j,     kind = rk ) * 2.0D+00 * pi &
             / real ( theta_num, kind = rk )

      call tp_to_xyz ( theta1, phi1, x11 )
      call tp_to_xyz ( theta2, phi1, x21 )
      call tp_to_xyz ( theta2, phi2, x2 )

      call sphere01_triangle_vertices_to_area ( x11, x2, x21, area )

      call sphere01_triangle_vertices_to_midpoints ( x11, x2, x21, &
        m1, m2, m3 )

      call f ( 1, m1, v )
      n = n + 1
      result = result + area * v(1) / 3.0D+00
      call f ( 1, m2, v )
      n = n + 1
      result = result + area * v(1) / 3.0D+00
      call f ( 1, m3, v )
      n = n + 1
      result = result + area * v(1) / 3.0D+00

    end do

  end if

  return
end
subroutine sphere01_quad_llv ( f, h, n, result )

!*****************************************************************************80
!
!! SPHERE01_QUAD_LLV: longitude/latitude grid with vertex rule.
!
!  Discussion:
!
!    The sphere is broken up into spherical triangles, whose sides
!    do not exceed the length H.  Then the function is evaluated
!    at the vertices, and the average is multiplied by the area.
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
!    Input, external :: F, evaluates the integrand, of the form:
!      subroutine f ( n, x, v )
!      integer n
!      real ( kind = rk ) v(n)
!      real ( kind = rk ) x(3,n)
!
!    Input, real ( kind = rk ) H, the maximum length of a side of the spherical
!    quadrilaterals.
!
!    Output, integer N, the number of points used.
!
!    Output, real ( kind = rk ) RESULT, the approximate integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) area
  external f
  real ( kind = rk ) h
  integer i
  integer j
  integer n
  real ( kind = rk ) phi
  integer phi_num
  real ( kind = rk ) phi1
  real ( kind = rk ) phi2
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) result
  real ( kind = rk ) sector_area
  real ( kind = rk ) sphere_area
  real ( kind = rk ) theta
  integer theta_num
  real ( kind = rk ) theta1
  real ( kind = rk ) theta2
  real ( kind = rk ) v(1)
  real ( kind = rk ) x(3)
  real ( kind = rk ) x1(3)
  real ( kind = rk ) x11(3)
  real ( kind = rk ) x12(3)
  real ( kind = rk ) x2(3)
  real ( kind = rk ) x21(3)
  real ( kind = rk ) x22(3)
!
!  Choose PHI and THETA counts that make short sides.
!
  phi_num = int ( pi / h )

  if ( h * real ( phi_num, kind = rk ) < pi ) then
    phi_num = phi_num + 1
  end if

  theta_num = int ( 2.0D+00 * pi / h )

  if ( h * real ( theta_num, kind = rk ) < pi ) then
    theta_num = theta_num + 1
  end if

  n = 0
  result = 0.0D+00
!
!  Only one THETA (and hence, only one PHI.)
!
  if ( theta_num == 1 ) then

    sphere_area = 4.0D+00 * pi

    theta = 0.0D+00
    phi = pi / 2.0D+00
    call tp_to_xyz ( theta, phi, x )
    call f ( 1, x, v )
    result = sphere_area * v(1)
!
!  Several THETA's, one PHI.
!
  else if ( phi_num == 1 ) then

    sphere_area = 4.0D+00 * pi
    sector_area = sphere_area / real ( theta_num, kind = rk )

    result = 0.0D+00

    do j = 1, theta_num

      theta = real ( ( j - 1 ) * 2, kind = rk ) * pi &
        / real ( theta_num, kind = rk )
      phi = pi / 2.0D+00
      call tp_to_xyz ( theta, phi, x )
      call f ( 1, x, v )
      n = n + 1
      result = result + sector_area * v(1)

    end do
!
!  At least two PHI's.
!
  else

    result = 0.0D+00
!
!  Picture:
!
!        V1
!       /  \
!      /    \
!    V12----V22
!
    phi1 = 0.0D+00
    phi2 = pi / real ( phi_num, kind = rk )

    do j = 1, theta_num

      theta1 = real ( j - 1, kind = rk ) * 2.0D+00 * pi &
             / real ( theta_num, kind = rk )
      theta2 = real ( j,     kind = rk ) * 2.0D+00 * pi &
             / real ( theta_num, kind = rk )

      call tp_to_xyz ( theta1, phi1, x1 )
      call tp_to_xyz ( theta1, phi2, x12 )
      call tp_to_xyz ( theta2, phi2, x22 )

      call sphere01_triangle_vertices_to_area ( x1, x12, x22, area )

      call f ( 1, x1, v )
      n = n + 1
      result = result + area * v(1) / 3.0D+00
      call f ( 1, x12, v )
      n = n + 1
      result = result + area * v(1) / 3.0D+00
      call f ( 1, x22, v )
      n = n + 1
      result = result + area * v(1) / 3.0D+00

    end do
!
!  Picture:
!
!    V11--V21
!     | \  |
!     |  \ |
!    V12--V22
!
    do i = 2, phi_num-1

      phi1 = real ( i - 1, kind = rk ) * pi &
           / real ( phi_num, kind = rk )
      phi2 = real ( i,     kind = rk ) * pi &
           / real ( phi_num, kind = rk )

      do j = 1, theta_num

        theta1 = real ( j - 1, kind = rk ) * 2.0D+00 * pi &
               / real ( theta_num, kind = rk )
        theta2 = real ( j,     kind = rk ) * 2.0D+00 * pi &
               / real ( theta_num, kind = rk )

        call tp_to_xyz ( theta1, phi1, x11 )
        call tp_to_xyz ( theta2, phi1, x21 )
        call tp_to_xyz ( theta1, phi2, x12 )
        call tp_to_xyz ( theta2, phi2, x22 )

        call sphere01_triangle_vertices_to_area ( x11, x12, x22, area )

        call f ( 1, x11, v )
        n = n + 1
        result = result + area * v(1) / 3.0D+00
        call f ( 1, x12, v )
        n = n + 1
        result = result + area * v(1) / 3.0D+00
        call f ( 1, x22, v )
        n = n + 1
        result = result + area * v(1) / 3.0D+00

        call sphere01_triangle_vertices_to_area ( x22, x21, x11, area )

        call f ( 1, x22, v )
        n = n + 1
        result = result + area * v(1) / 3.0D+00
        call f ( 1, x21, v )
        n = n + 1
        result = result + area * v(1) / 3.0D+00
        call f ( 1, x11, v )
        n = n + 1
        result = result + area * v(1) / 3.0D+00

      end do

    end do
!
!  Picture:
!
!    V11----V21
!      \    /
!       \  /
!        V2
!
    phi1 = real ( phi_num - 1, kind = rk ) * pi &
         / real ( phi_num, kind = rk )
    phi2 =                                  pi

    do j = 1, theta_num

      theta1 = real ( j - 1, kind = rk ) * 2.0D+00 * pi &
             / real ( theta_num, kind = rk )
      theta2 = real ( j,     kind = rk ) * 2.0D+00 * pi &
             / real ( theta_num, kind = rk )

      call tp_to_xyz ( theta1, phi1, x11 )
      call tp_to_xyz ( theta2, phi1, x21 )
      call tp_to_xyz ( theta2, phi2, x2 )

      call sphere01_triangle_vertices_to_area ( x11, x2, x21, area )

      call f ( 1, x11, v )
      n = n + 1
      result = result + area * v(1) / 3.0D+00
      call f ( 1, x2, v )
      n = n + 1
      result = result + area * v(1) / 3.0D+00
      call f ( 1, x21, v )
      n = n + 1
      result = result + area * v(1) / 3.0D+00

    end do

  end if

  return
end
subroutine sphere01_quad_mc ( f, h, n, result )

!*****************************************************************************80
!
!! SPHERE01_QUAD_MC uses the Monte Carlo rule for sphere quadrature.
!
!  Discussion:
!
!    A number of points N are chosen at random on the sphere, with N
!    being determined so that, if the points were laid out on a regular
!    grid, the average spacing would be no more than H.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    13 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, external :: F, evaluates the integrand, of the form:
!      subroutine f ( n, x, v )
!      integer n
!      real ( kind = rk ) v(n)
!      real ( kind = rk ) x(3,n)
!
!    Input, real ( kind = rk ) H, the maximum length of a side of the spherical
!    quadrilaterals.
!
!    Input, integer N, the number of points used.
!
!    Output, real ( kind = rk ) RESULT, the approximate integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  external f
  real ( kind = rk ) h
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) result
  real ( kind = rk ) sphere_area
  real ( kind = rk ) v(n)
  real ( kind = rk ) x(3,n)

  call r8_fake_use ( h )

  sphere_area = 4.0D+00 * pi

  call sphere01_sample_3d ( n, x )

  call f ( n, x, v )

  result = sphere_area * sum ( v(1:n) ) / real ( n, kind = rk )

  return
end
subroutine sphere01_quad_mc_size ( h, n )

!*****************************************************************************80
!
!! SPHERE01_QUAD_MC_SIZE sizes a Monte Carlo rule for sphere quadrature.
!
!  Discussion:
!
!    A number of points N are chosen at random on the sphere, with N
!    being determined so that, if the points were laid out on a regular
!    grid, the average spacing would be no more than H.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) H, the maximum length of a side of the spherical
!    quadrilaterals.
!
!    Output, integer N, the number of points to use.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) h
  integer n
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) sphere_area
!
!  The sphere's area is 4 * PI.
!  Choose N so that we divide this area into N subareas of PI * H * H.
!
  sphere_area = 4.0D+00 * pi

  n = int ( sphere_area / h**2 )
  n = max ( n, 1 )

  return
end
subroutine sphere01_sample_3d ( n, x )

!*****************************************************************************80
!
!! SPHERE01_SAMPLE_3D picks random points on a sphere in 3D.
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
!    Input, integer N, the number of samples.
!
!    Output, real ( kind = rk ) X(3,N), the sample points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer j
  real ( kind = rk ) phi
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) theta
  real ( kind = rk ) vdot
  real ( kind = rk ) x(3,n)

  do j = 1, n
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

    phi = acos ( vdot )
!
!  Pick a uniformly random rotation between 0 and 2 Pi around the
!  axis of the Z vector.
!
    call random_number ( harvest = theta )
    theta = 2.0D+00 * pi * theta

    x(1,j) = cos ( theta ) * sin ( phi )
    x(2,j) = sin ( theta ) * sin ( phi )
    x(3,j) =                 cos ( phi )

  end do

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
!    The area of a spherical triangle is:
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
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
!
!  Apply Girard's formula.
!
  area = a + b + c - pi

  return
end
subroutine sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1, f2, f3, &
  node_xyz )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_PROJECT projects from plane to spherical triangle.
!
!  Discussion:
!
!    We assume that points A, B and C lie on the unit sphere, and they
!    thus define a spherical triangle.
!
!    They also, of course, define a planar triangle.
!
!    Let (F1,F2,F3) be the barycentric coordinates of a point in this 
!    planar triangle.
!
!    This function determines the coordinates of the point in the planar
!    triangle identified by the barycentric coordinates, and returns the
!    coordinates of the projection of that point onto the unit sphere.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A_XYZ(3), B_XYZ(3), C_XYZ(3), the coordinates
!    of the points A, B, and C.
!
!    Input, integer F1, F2, F3, the barycentric coordinates
!    of a point in the triangle ABC.  Normally, these coordinates would
!    be real numbers, and would sum to 1.  For convenience, we allow these
!    to be integers which must be divided by F1+F2+F3.
!
!    Output, real ( kind = rk ) NODE_XYZ(3), the coordinates of the 
!    point on the unit sphere which is the projection of the point on the plane
!    whose barycentric coordinates with respect to A, B, and C is
!    (F1,F2,F3)/(F1+F2+F3).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a_xyz(3)
  real ( kind = rk ) b_xyz(3)
  real ( kind = rk ) c_xyz(3)
  integer f1
  integer f2
  integer f3
  real ( kind = rk ) node_norm
  real ( kind = rk ) node_xyz(3)
  real ( kind = rk ) r8vec_norm

  node_xyz(1:3) = &
    ( real ( f1,           kind = rk ) * a_xyz(1:3)   &
    + real (      f2,      kind = rk ) * b_xyz(1:3)   &
    + real (           f3, kind = rk ) * c_xyz(1:3) ) &
    / real ( f1 + f2 + f3, kind = rk )

  node_norm = r8vec_norm ( 3, node_xyz(1:3) )

  node_xyz(1:3) = node_xyz(1:3) / node_norm

  return
end
subroutine sphere01_triangle_project2 ( a_xyz, b_xyz, c_xyz, f1, f2, f3, &
  node_xyz )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_PROJECT2 projects from plane to spherical triangle.
!
!  Discussion:
!
!    We assume that points A, B and C lie on the unit sphere, and they
!    thus define a spherical triangle.
!
!    They also, of course, define a planar triangle.
!
!    Let (F1,F2,F3) be the barycentric coordinates of a point in this 
!    planar triangle.
!
!    This function determines the coordinates of the point in the planar
!    triangle identified by the barycentric coordinates, and returns the
!    coordinates of the projection of that point onto the unit sphere.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A_XYZ(3), B_XYZ(3), C_XYZ(3), the coordinates
!    of the points A, B, and C.
!
!    Input, integer F1, F2, F3, the barycentric coordinates
!    of a point in the triangle ABC.  Normally, these coordinates would
!    be real numbers, and would sum to 1.  For convenience, we allow these
!    to be integers which must be divided by F1+F2+F3.
!
!    Output, real ( kind = rk ) NODE_XYZ(3), the coordinates of the 
!    point on the unit sphere which is the projection of the point on the 
!    plane whose barycentric coordinates with respect to A, B, and C is
!    (F1,F2,F3)/(F1+F2+F3).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a_xyz(3)
  real ( kind = rk ) ab(3)
  real ( kind = rk ) ac(3)
  real ( kind = rk ) acn(3)
  real ( kind = rk ) acp(3)
  real ( kind = rk ) angle
  real ( kind = rk ) b_xyz(3)
  real ( kind = rk ) bn(3)
  real ( kind = rk ) bp(3)
  real ( kind = rk ) c_xyz(3)
  real ( kind = rk ) cn(3)
  real ( kind = rk ) cp(3)
  integer f1
  integer f2
  integer f3
  real ( kind = rk ) node_xyz(3)
  real ( kind = rk ) r8vec_norm
  real ( kind = rk ) theta_ab
  real ( kind = rk ) theta_ac
  real ( kind = rk ) theta_bc
!
!  This check avoids 0/0 calculations later.
!
  if ( f2 == 0 .and. f3 == 0 ) then
    node_xyz(1:3) = a_xyz(1:3)
    return
  else if ( f1 == 0 .and. f3 == 0 ) then
    node_xyz(1:3) = b_xyz(1:3)
    return
  else if ( f1 == 0 .and. f2 == 0 ) then
    node_xyz(1:3) = c_xyz(1:3)
    return
  end if
!
!  Determine the angular distances (A,B) and (A,C).
!
  call sphere01_distance_xyz ( a_xyz, b_xyz, theta_ab )

  call sphere01_distance_xyz ( a_xyz, c_xyz, theta_ac )
!
!  Polarize B = BP + BN
!  Normalize BN, 
!  Same for C.
!
  call r8vec_polarize ( 3, b_xyz, a_xyz, bn, bp )
  bn(1:3) = bn(1:3) / r8vec_norm ( 3, bn )

  call r8vec_polarize ( 3, c_xyz, a_xyz, cn, cp )
  cn(1:3) = cn(1:3) / r8vec_norm ( 3, cn )
!
!  Determine AB and AC that use cos ( ( F2 + F3 ) / ( F1 + F2 + F3 ) ) of A
!  and cos ( F1 / ( F1 + F2 + F3 ) ) of B or C.
!
  angle = ( real ( f2 + f3, kind = rk ) * theta_ab ) &
    / real ( f1 + f2 + f3, kind = rk )
  ab(1:3) = cos ( angle ) * a_xyz(1:3) + sin ( angle ) * bn(1:3)

  angle = ( real ( f2 + f3, kind = rk ) * theta_ac ) &
    / real ( f1 + f2 + f3, kind = rk )
  ac(1:3) = cos ( angle ) * a_xyz(1:3) + sin ( angle ) * cn(1:3)
!
!  Determine the angular distance between AB and AC.
!
  call sphere01_distance_xyz ( ab(1:3), ac(1:3), theta_bc )
!
!  Polarize AC = ACP + ACN, normalize ACN.
!
  call r8vec_polarize ( 3, ac, ab, acn, acp )
  acn(1:3) = acn(1:3) / r8vec_norm ( 3, acn )
!
!  The interval between AB and AC is marked by F2+F3+1 vertices 0 through F2+F3.
!
  angle = ( real ( f3, kind = rk ) * theta_bc ) / real ( f2 + f3, kind = rk )

  node_xyz(1:3) = cos ( angle ) * ab(1:3) + sin ( angle ) * acn(1:3)

  return
end
subroutine sphere01_triangle_sample ( n, v1, v2, v3, x )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_SAMPLE: sample points from triangle on unit sphere.
!
!  Discussion:
!
!    The sphere has center 0 and radius 1.
!
!    A spherical triangle on the surface of the unit sphere contains those 
!    points with radius R = 1, bounded by the vertices V1, V2, V3.
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
!  Reference:
!
!    James Arvo,
!    Stratified sampling of spherical triangles,
!    Computer Graphics Proceedings, Annual Conference Series, 
!    ACM SIGGRAPH '95, pages 437-438, 1995.
!
!  Parameters:
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) V1(3), V2(3), V3(3), the XYZ coordinates of
!    the vertices of the spherical triangle.
!
!    Output, real ( kind = rk ) X(3,N), the XYZ coordinates of the 
!    sample points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a
  real ( kind = rk ) alpha
  real ( kind = rk ) area
  real ( kind = rk ) area_hat
  real ( kind = rk ) b
  real ( kind = rk ) beta
  real ( kind = rk ) c
  real ( kind = rk ) gamma
  integer j
  real ( kind = rk ) q
  real ( kind = rk ) r8vec_norm
  real ( kind = rk ) s
  real ( kind = rk ) t
  real ( kind = rk ) u
  real ( kind = rk ) v
  real ( kind = rk ) v1(3)
  real ( kind = rk ) v2(3)
  real ( kind = rk ) v3(3)
  real ( kind = rk ) v31(3)
  real ( kind = rk ) v4(3)
  real ( kind = rk ) v42(3)
  real ( kind = rk ) w
  real ( kind = rk ) x(3,n)
  real ( kind = rk ) xsi1
  real ( kind = rk ) xsi2
  real ( kind = rk ) z

  call sphere01_triangle_vertices_to_sides ( v1, v2, v3, a, b, c )

  call sphere01_triangle_sides_to_angles ( a, b, c, alpha, beta, gamma )

  call sphere01_triangle_angles_to_area ( alpha, beta, gamma, area )

  do j = 1, n
!
!  Select the new area.
!
    call random_number ( harvest = xsi1 )
    area_hat = xsi1 * area
!
!  Compute the sine and cosine of the angle phi.
!
    s = sin ( area_hat - alpha )
    t = cos ( area_hat - alpha )
!
!  Compute the pair that determines beta_hat.
!
    u = t - cos ( alpha )
    v = s + sin ( alpha ) * cos ( c )
!
!  Q is the cosine of the new edge length b_hat.
!
    q = ( ( v * t - u * s ) * cos ( alpha ) - v ) &
      / ( ( v * s + u * t ) * sin ( alpha ) )
!
!  We very occasionally get a Q value out of bounds.
!
    q = max ( q, - 1.0D+00 )
    q = min ( q, + 1.0D+00 )
!
!  V31 = normalized ( V3 - ( V3 dot V1 ) * V1 )
!
    w = dot_product ( v3, v1 )
    v31(1:3) = v3(1:3) - w * v1(1:3)
    v31(1:3) = v31(1:3) / r8vec_norm ( 3, v31(1:3) )
!
!  V4 is the third vertex of the subtriangle V1, V2, V4.
!
    v4(1:3) = q * v1(1:3) + sqrt ( 1.0D+00 - q * q ) * v31(1:3)
!
!  Select cos theta, which will sample along the edge from V2 to V4.
!
    call random_number ( harvest = xsi2 )
    z = 1.0D+00 - xsi2 * ( 1.0D+00 - dot_product ( v4, v2 ) )
!
!  V42 = normalized ( V4 - ( V4 dot V2 ) * V2 )
!
    w = dot_product ( v4, v2 )
    v42(1:3) = v4(1:3) - w * v2(1:3)
    v42(1:3) = v42(1:3) / r8vec_norm ( 3, v42(1:3) )
!
!  Construct the point.
!
    x(1:3,j) = z * v2(1:3) + sqrt ( 1.0D+00 - z * z ) * v42(1:3)

  end do

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
subroutine sphere01_triangle_vertices_to_area ( v1, v2, v3, area )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_VERTICES_TO_AREA computes the area of a spherical triangle.
!
!  Discussion:
!
!    A sphere in 3D satisfies the equation:
!
!      X^2 + Y^2 + Z^2 = 1
!
!    A spherical triangle is specified by three points on the surface
!    of the sphere.
!
!    The area formula is known as Girard's formula.
!
!    The area of a spherical triangle is:
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
!    A sphere in 3D satisfies the equation:
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
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    31 May 2001
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
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = 10 ) time
  integer values(8)
  integer y
  character ( len = 5 ) zone

  call date_and_time ( date, time, zone, values )

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

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine tp_to_xyz ( theta, phi, v )

!*****************************************************************************80
!
!! TP_TO_XYZ converts spherical TP coordinates to XYZ coordinates.
!
!  Discussion:
!
!    The sphere is assumed to have radius 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) THETA, PHI, the spherical coordinates 
!    of a point.
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
