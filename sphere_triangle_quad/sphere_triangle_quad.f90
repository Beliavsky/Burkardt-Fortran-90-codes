function r8_acos ( c )

!*****************************************************************************80
!
!! r8_acos() computes the arc cosine function, with argument truncation.
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
function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = rk ) value.
!
!    For now, the input quantity SEED is an integer variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
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
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
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
!    Input/output, integer SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = rk ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: i4_huge = 2147483647
  integer k
  real ( kind = rk ) r8_uniform_01
  integer seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = rk ) * 4.656612875D-10

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
subroutine r8vec_transpose_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_TRANSPOSE_PRINT prints an R8VEC "transposed".
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Example:
!
!    A = (/ 1.0, 2.1, 3.2, 4.3, 5.4, 6.5, 7.6, 8.7, 9.8, 10.9, 11.0 /)
!    TITLE = 'My vector:  '
!
!    My vector:
!        1.0    2.1    3.2    4.3    5.4
!        6.5    7.6    8.7    9.8   10.9
!       11.0
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
  integer ihi
  integer ilo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  do ilo = 1, n, 5
    ihi = min ( ilo + 5 - 1, n )
    write ( *, '(5g14.6)' ) a(ilo:ihi)
  end do

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
subroutine sphere01_sample ( n, seed, x )

!*****************************************************************************80
!
!! SPHERE01_SAMPLE picks random points on the unit sphere in 3D.
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
!    Input/output, integer SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = rk ) X(3,N), the sample points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer j
  real ( kind = rk ) phi
  real ( kind = rk ) r8_acos
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) r8_uniform_01
  integer seed
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
    vdot = r8_uniform_01 ( seed )
    vdot = 2.0D+00 * vdot - 1.0D+00

    phi = r8_acos ( vdot )
!
!  Pick a uniformly random rotation between 0 and 2 Pi around the
!  axis of the Z vector.
!
    theta = r8_uniform_01 ( seed )
    theta = 2.0D+00 * r8_pi * theta

    x(1,j) = cos ( theta ) * sin ( phi )
    x(2,j) = sin ( theta ) * sin ( phi )
    x(3,j) = cos ( phi )

  end do

  return
end
subroutine sphere01_triangle_angles_to_area ( a, b, c, area )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_ANGLES_TO_AREA: area of a triangle on the unit sphere.
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
!    The area of a spherical triangle on the unit sphere is:
!
!      AREA = A + B + C - PI
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
subroutine sphere01_triangle_quad_00 ( n, v1, v2, v3, f, seed, result )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_QUAD_00: quadrature over a triangle on the unit sphere.
!
!  Discussion:
!
!    This is a Monte Carlo approach.
!
!    The integral is approximated by averaging the values at N random points,
!    multiplied by the area.
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
!    Input, integer N, the number of sample points.
!
!    Input, real ( kind = rk ) V1(3), V2(3), V3(3), the XYZ coordinates of
!    the vertices of the triangle.
!
!    Input, external, real ( kind = rk ) F, evaluates the integrand at X.
!
!    Input/output, integer SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = rk ) RESULT, the approximate integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) area
  real ( kind = rk ), external :: f
  integer j
  real ( kind = rk ) quad
  real ( kind = rk ) result
  integer seed
  real ( kind = rk ) v1(3)
  real ( kind = rk ) v2(3)
  real ( kind = rk ) v3(3)
  real ( kind = rk ) vc(3,n)

  call sphere01_triangle_vertices_to_area ( v1, v2, v3, area )

  call sphere01_triangle_sample ( n, v1, v2, v3, seed, vc )

  quad = 0.0D+00
  do j = 1, n
    quad = quad + f ( vc(1:3,j) )
  end do

  result = quad * area / real ( n, kind = rk )

  return
end
subroutine sphere01_triangle_quad_01 ( v1, v2, v3, f, result )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_QUAD_01: quadrature over a triangle on the unit sphere.
!
!  Discussion:
!
!    The integral is approximated by the value at the centroid,
!    multiplied by the area.
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
!    Input, real ( kind = rk ) V1(3), V2(3), V3(3), the XYZ coordinates of
!    the vertices of the triangle.
!
!    Input, external, real ( kind = rk ) F, evaluates the integrand at X.
!
!    Output, real ( kind = rk ) RESULT, the approximate integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) area
  real ( kind = rk ), external :: f
  real ( kind = rk ) quad
  real ( kind = rk ) result
  real ( kind = rk ) v1(3)
  real ( kind = rk ) v2(3)
  real ( kind = rk ) v3(3)
  real ( kind = rk ) vc(3)

  call sphere01_triangle_vertices_to_area ( v1, v2, v3, area )

  call sphere01_triangle_vertices_to_centroid ( v1, v2, v3, vc )

  quad = f ( vc )
  result = quad * area

  return
end
subroutine sphere01_triangle_quad_02 ( v1, v2, v3, f, result )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_QUAD_02: quadrature over a triangle on the unit sphere.
!
!  Discussion:
!
!    The integral is approximated by the average of the vertex values,
!    multiplied by the area.
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
!    Input, real ( kind = rk ) V1(3), V2(3), V3(3), the XYZ coordinates of
!    the vertices of the triangle.
!
!    Input, external, real ( kind = rk ) F, evaluates the integrand at X.
!
!    Output, real ( kind = rk ) RESULT, the approximate integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) area
  real ( kind = rk ), external :: f
  real ( kind = rk ) quad
  real ( kind = rk ) result
  real ( kind = rk ) v1(3)
  real ( kind = rk ) v2(3)
  real ( kind = rk ) v3(3)

  call sphere01_triangle_vertices_to_area ( v1, v2, v3, area )

  quad = ( f ( v1 ) + f ( v2 ) + f ( v3 ) ) / 3.0D+00

  result = quad * area

  return
end
subroutine sphere01_triangle_quad_03 ( v1, v2, v3, f, result )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_QUAD_03: quadrature over a triangle on the unit sphere.
!
!  Discussion:
!
!    The integral is approximated by the average of the midside values,
!    multiplied by the area.
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
!    Input, real ( kind = rk ) V1(3), V2(3), V3(3), the XYZ coordinates of
!    the vertices of the triangle.
!
!    Input, external, real ( kind = rk ) F, evaluates the integrand at X.
!
!    Output, real ( kind = rk ) RESULT, the approximate integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) area
  real ( kind = rk ), external :: f
  real ( kind = rk ) quad
  real ( kind = rk ) result
  real ( kind = rk ) v1(3)
  real ( kind = rk ) v2(3)
  real ( kind = rk ) v3(3)
  real ( kind = rk ) v4(3)
  real ( kind = rk ) v5(3)
  real ( kind = rk ) v6(3)

  call sphere01_triangle_vertices_to_area ( v1, v2, v3, area )

  call sphere01_triangle_vertices_to_midpoints ( v1, v2, v3, v4, v5, v6 )

  quad = ( f ( v4 ) + f ( v5 ) + f ( v6 ) ) / 3.0D+00

  result = quad * area

  return
end
subroutine sphere01_triangle_quad_icos1c ( a_xyz, b_xyz, c_xyz, factor, &
  fun, node_num, result )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_QUAD_ICOS1C: centroid rule, subdivide then project.
!
!  Discussion:
!
!    This function estimates an integral over a spherical triangle on the
!    unit sphere.
!
!    This function sets up an icosahedral grid, and subdivides each
!    edge of the icosahedron into FACTOR subedges.  These edges define a grid
!    within each triangular icosahedral face.  The centroids of these
!    triangles can be determined.  All of these calculations are done,
!    essentially, on the FLAT faces of the icosahedron.  Only then are
!    the triangle vertices and centroids projected to the sphere.  
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A_XYZ(3), B_XYZ(3), C_XYZ(3), the vertices
!    of the spherical triangle.
!
!    Input, integer FACTOR, the subdivision factor, which must
!    be at least 1.
!
!    Input, external, real ( kind = rk ) FUN, evaluates the integrand at X.
!
!    Output, integer NODE_NUM, the number of evaluation points.
!
!    Output, real ( kind = rk ) RESULT, the estimated integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a_xyz(3)
  real ( kind = rk ) a2_xyz(3)
  real ( kind = rk ) area
  real ( kind = rk ) area_total
  real ( kind = rk ) b_xyz(3)
  real ( kind = rk ) b2_xyz(3)
  real ( kind = rk ) c_xyz(3)
  real ( kind = rk ) c2_xyz(3)
  integer f1
  integer f2
  integer f3
  integer factor
  real ( kind = rk ), external :: fun
  integer node_num
  real ( kind = rk ) node_xyz(3)
  real ( kind = rk ) result
  real ( kind = rk ) v
!
!  Initialize the integral data.
!
  result = 0.0D+00
  area_total = 0.0D+00
  node_num = 0
!
!  Some subtriangles will have the same direction as the triangle.
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

      v = fun ( node_xyz )    

      node_num = node_num + 1
      result = result + area * v
      area_total = area_total + area

    end do
  end do
!
!  The other subtriangles have the opposite direction from the triangle.
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

      v = fun ( node_xyz )  

      node_num = node_num + 1  
      result = result + area * v
      area_total = area_total + area

    end do
  end do

  return
end
subroutine sphere01_triangle_quad_icos1m ( a_xyz, b_xyz, c_xyz, factor, &
  fun, node_num, result )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_QUAD_ICOS1M: midpoint rule, subdivide then project.
!
!  Discussion:
!
!    This function estimates an integral over a spherical triangle on the
!    unit sphere.
!
!    This function sets up an icosahedral grid, and subdivides each
!    edge of the icosahedron into FACTOR subedges.  These edges define a grid
!    within each triangular icosahedral face.  The midpoints of the edges
!    of these triangles can be determined.  All of these calculations are done,
!    essentially, on the FLAT faces of the icosahedron.  Only then are
!    the triangle vertices and midpoints projected to the sphere.  
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A_XYZ(3), B_XYZ(3), C_XYZ(3), the vertices
!    of the spherical triangle.
!
!    Input, integer FACTOR, the subdivision factor, which must
!    be at least 1.
!
!    Input, external, real ( kind = rk ) FUN, evaluates the integrand at X.
!
!    Output, integer NODE_NUM, the number of evaluation points.
!
!    Output, real ( kind = rk ) RESULT, the estimated integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a_xyz(3)
  real ( kind = rk ) a2_xyz(3)
  real ( kind = rk ) a3_xyz(3)
  real ( kind = rk ) area
  real ( kind = rk ) area_total
  real ( kind = rk ) b_xyz(3)
  real ( kind = rk ) b2_xyz(3)
  real ( kind = rk ) b3_xyz(3)
  real ( kind = rk ) c_xyz(3)
  real ( kind = rk ) c2_xyz(3)
  real ( kind = rk ) c3_xyz(3)
  integer f1
  integer f2
  integer f3
  integer factor
  real ( kind = rk ), external :: fun
  integer node_num
  real ( kind = rk ) result
  real ( kind = rk ) va
  real ( kind = rk ) vb
  real ( kind = rk ) vc
!
!  Initialize the integral data.
!
  result = 0.0D+00
  area_total = 0.0D+00
  node_num = 0
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
      va = fun ( a3_xyz )
      vb = fun ( b3_xyz )   
      vc = fun ( c3_xyz )   
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
      va = fun ( a3_xyz )
      vb = fun ( b3_xyz )   
      vc = fun ( c3_xyz )  
      result = result + area * ( va + vb + vc ) / 3.0D+00
      area_total = area_total + area

    end do
  end do

  return
end
subroutine sphere01_triangle_quad_icos1v ( a_xyz, b_xyz, c_xyz, factor, &
  fun, node_num, result )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_QUAD_ICOS1V: vertex rule, subdivide then project.
!
!  Discussion:
!
!    This function estimates an integral over a spherical triangle on the
!    unit sphere.
!
!    This function sets up an icosahedral grid, and subdivides each
!    edge of the icosahedron into FACTOR subedges.  These edges define a grid
!    within each triangular icosahedral face.   All of these calculations are 
!    done, essentially, on the FLAT faces of the icosahedron.  Only then are
!    the triangle vertices projected to the sphere.  
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A_XYZ(3), B_XYZ(3), C_XYZ(3), the vertices
!    of the spherical triangle.
!
!    Input, integer FACTOR, the subdivision factor, which must
!    be at least 1.
!
!    Input, external, real ( kind = rk ) FUN, evaluates the integrand at X.
!
!    Output, integer NODE_NUM, the number of evaluation points.
!
!    Output, real ( kind = rk ) RESULT, the estimated integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a_xyz(3)
  real ( kind = rk ) a2_xyz(3)
  real ( kind = rk ) area
  real ( kind = rk ) area_total
  real ( kind = rk ) b_xyz(3)
  real ( kind = rk ) b2_xyz(3)
  real ( kind = rk ) c_xyz(3)
  real ( kind = rk ) c2_xyz(3)
  integer f1
  integer f2
  integer f3
  integer factor
  real ( kind = rk ), external :: fun
  integer node_num
  real ( kind = rk ) result
  real ( kind = rk ) va
  real ( kind = rk ) vb
  real ( kind = rk ) vc
!
!  Initialize the integral data.
!
  result = 0.0D+00
  area_total = 0.0D+00
  node_num = 0
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
      va = fun ( a2_xyz )
      vb = fun ( b2_xyz )   
      vc = fun ( c2_xyz )   
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
      va = fun ( a2_xyz )
      vb = fun ( b2_xyz )   
      vc = fun ( c2_xyz )  
      result = result + area * ( va + vb + vc ) / 3.0D+00
      area_total = area_total + area

    end do
  end do

  return
end
subroutine sphere01_triangle_quad_icos2v ( a_xyz, b_xyz, c_xyz, factor, &
  fun, node_num, result )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_QUAD_ICOS2V: vertex rule, subdivide then project.
!
!  Discussion:
!
!    This function estimates an integral over a spherical triangle on the
!    unit sphere.
!
!    This function sets up an icosahedral grid, and subdivides each
!    edge of the icosahedron into FACTOR subedges.  These edges define a grid
!    within each triangular icosahedral face.   All of these calculations are 
!    done, essentially, on the FLAT faces of the icosahedron.  Only then are
!    the triangle vertices projected to the sphere.  
!
!    This function uses a more sophisticated projection scheme than
!    SPHERE01_TRIANGLE_QUAD_ICOS1V, but this does not seem to improve
!    the results significantly.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A_XYZ(3), B_XYZ(3), C_XYZ(3), the vertices
!    of the spherical triangle.
!
!    Input, integer FACTOR, the subdivision factor, which must
!    be at least 1.
!
!    Input, external, real ( kind = rk ) FUN, evaluates the integrand at X.
!
!    Output, integer NODE_NUM, the number of evaluation points.
!
!    Output, real ( kind = rk ) RESULT, the estimated integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a_xyz(3)
  real ( kind = rk ) a2_xyz(3)
  real ( kind = rk ) area
  real ( kind = rk ) area_total
  real ( kind = rk ) b_xyz(3)
  real ( kind = rk ) b2_xyz(3)
  real ( kind = rk ) c_xyz(3)
  real ( kind = rk ) c2_xyz(3)
  integer f1
  integer f2
  integer f3
  integer factor
  real ( kind = rk ), external :: fun
  integer node_num
  real ( kind = rk ) result
  real ( kind = rk ) va
  real ( kind = rk ) vb
  real ( kind = rk ) vc
!
!  Initialize the integral data.
!
  result = 0.0D+00
  area_total = 0.0D+00
  node_num = 0
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
      va = fun ( a2_xyz )
      vb = fun ( b2_xyz )   
      vc = fun ( c2_xyz )   
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
      va = fun ( a2_xyz )
      vb = fun ( b2_xyz )   
      vc = fun ( c2_xyz )  
      result = result + area * ( va + vb + vc ) / 3.0D+00
      area_total = area_total + area

    end do
  end do

  return
end
subroutine sphere01_triangle_sample ( n, v1, v2, v3, seed, x )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_SAMPLE: sample spherical triangle on unit sphere.
!
!  Discussion:
!
!    The unit sphere has center 0 and radius 1.
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
!    Input/output, integer SEED, a seed for the random 
!    number generator.
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
  real ( kind = rk ) r8_uniform_01
  real ( kind = rk ) r8vec_norm
  real ( kind = rk ) s
  integer seed
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
!
!  Compute the sides, angles, and area of the spherical triangle;
!
  call sphere01_triangle_vertices_to_sides ( v1, v2, v3, a, b, c )

  call sphere01_triangle_sides_to_angles ( a, b, c, alpha, beta, gamma )

  call sphere01_triangle_angles_to_area ( alpha, beta, gamma, area )

  do j = 1, n
!
!  Select the new area.
!
    xsi1 = r8_uniform_01 ( seed )

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
    t = r8vec_norm ( 3, v31(1:3) )
    if ( 0.0D+00 < t ) then
      v31(1:3) = v31(1:3) / t
    end if
!
!  V4 is the third vertex of the subtriangle V1, V2, V4.
!
    v4(1:3) = q * v1(1:3) + sqrt ( 1.0D+00 - q * q ) * v31(1:3)
!
!  Select cos theta, which will sample along the edge from V2 to V4.
!
    xsi2 = r8_uniform_01 ( seed )

    z = 1.0D+00 - xsi2 * ( 1.0D+00 - dot_product ( v4, v2 ) )
!
!  V42 = normalized ( V4 - ( V4 dot V2 ) * V2 )
!
    w = dot_product ( v4, v2 )

    v42(1:3) = v4(1:3) - w * v2(1:3)
    t = r8vec_norm ( 3, v42(1:3) )
    if ( 0.0D+00 < t ) then
      v42(1:3) = v42(1:3) / t
    end if
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
!! SPHERE01_TRIANGLE_SIDES_TO_ANGLES: angles of triangle on unit sphere.
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
!! SPHERE01_TRIANGLE_VERTICES_TO_AREA: area of triangle on unit sphere.
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
!    The area of a spherical triangle on the unit sphere is:
!
!      AREA = A + B + C - PI
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
!! SPHERE01_TRIANGLE_VERTICES_TO_CENTROID: centroid of triangle on unit sphere.
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
!! SPHERE01_TRIANGLE_VERTICES_TO_MIDPOINTS: midsides of triangle on unit sphere.
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
!! SPHERE01_TRIANGLE_VERTICES_TO_SIDES: sides of triangle on unit sphere.
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
  real ( kind = rk ) r8_acos
  real ( kind = rk ) v1(3)
  real ( kind = rk ) v2(3)
  real ( kind = rk ) v3(3)

  as = r8_acos ( dot_product ( v2(1:3), v3(1:3) ) )
  bs = r8_acos ( dot_product ( v3(1:3), v1(1:3) ) )
  cs = r8_acos ( dot_product ( v1(1:3), v2(1:3) ) )

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
  implicit none

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
