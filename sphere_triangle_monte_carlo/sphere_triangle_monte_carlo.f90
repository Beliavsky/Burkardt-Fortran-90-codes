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
subroutine r8vec_normalize ( n, a )

!*****************************************************************************80
!
!! R8VEC_NORMALIZE normalizes an R8VEC in the Euclidean norm.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The euclidean norm is also sometimes called the l2 or
!    least squares norm.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of the vector.
!
!    Input/output, real ( kind = rk ) A(N), the vector to be normalized.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) norm

  norm = sqrt ( sum ( a(1:n)**2 ) )

  if ( norm /= 0.0D+00 ) then
    a(1:n) = a(1:n) / norm
  end if

  return
end
subroutine sphere01_sample ( n, seed, x )

!*****************************************************************************80
!
!! SPHERE01_SAMPLE picks random points on the unit sphere.
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

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
