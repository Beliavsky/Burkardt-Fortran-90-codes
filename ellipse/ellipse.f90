function ellipse_area1 ( a, r )

!*****************************************************************************80
!
!! ellipse_area1() returns the area of an ellipse defined by a matrix.
!
!  Discussion:
!
!    The points X in the ellipse are described by a 2 by 2
!    positive definite symmetric matrix A, and a "radius" R, such that
!      X' * A * X <= R * R
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 April 2014
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) A(2,2), the matrix that describes
!    the ellipse.  A must be symmetric and positive definite.
!
!    real ( kind = rk ) R, the "radius" of the ellipse.
!
!  Output:
!
!    real ( kind = rk ) ELLIPSE_AREA1, the area of the ellipse.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(2,2)
  real ( kind = rk ) ellipse_area1
  real ( kind = rk ) r
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00

  ellipse_area1 = r ** 2 * r8_pi / sqrt ( a(1,1) * a(2,2) - a(2,1) * a(1,2) )

  return
end
function ellipse_area2 ( a, b, c, d )

!*****************************************************************************80
!
!! ellipse_area2() returns the area of an ellipse defined by an equation.
!
!  Discussion:
!
!    The ellipse is described by the formula
!      a x^2 + b xy + c y^2 = d
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 November 2016
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) A, B, C, coefficients on the left hand side.
!
!    real ( kind = rk ) D, the right hand side.
!
!  Output:
!
!    real ( kind = rk ) ELLIPSE_AREA2, the area of the ellipse.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) d
  real ( kind = rk ) ellipse_area2
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00

  ellipse_area2 = 2.0D+00 * d * d * r8_pi / sqrt ( 4.0D+00 * a * c - b * b )

  return
end
function ellipse_area3 ( r1, r2 )

!*****************************************************************************80
!
!! ellipse_area3() returns the area of an ellipse in 2D.
!
!  Discussion:
!
!    An ellipse in standard position has a center at the origin, and
!    axes aligned with the coordinate axes.  Any point P on the ellipse
!    satisfies
!
!      (  P(1) / R1 )^2 + ( P(2) / R2 )^2 == 1
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    28 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) R1, R2, the "radius" of the ellipse in the major
!    and minor axis directions.  A circle has these values equal.
!
!  Output:
!
!    real ( kind = rk ) ELLIPSE_AREA3, the area of the ellipse.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) ellipse_area3
  real ( kind = rk ) r1
  real ( kind = rk ) r2
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00

  ellipse_area3 = r8_pi * r1 * r2

  return
end
function ellipse_aspect_ratio ( a, b )

!*****************************************************************************80
!
!! ellipse_aspect_ratio() computes the aspect ratio of an ellipse.
!
!  Discussion:
!
!    The ellipse has major and minor semi-axes a and b.  In particular, it
!    could have the form:
!
!      (x/a)^2 + (y/b)^2 = 1
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 October 2022
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    John D Cook,
!    Eccentricity, Flattening, and Aspect Ratio,
!    https://www.johndcook.com/blog/2022/10/14/eccentricity-flatness-aspect/
!    Posted 14 October 2022.
!
!  Input:
!
!    real ( kind = rk ) A, B, the major and minor semi-axes.
!
!  Output:
!
!    real ( kind = rk ) ellipse_aspect_ratio: the aspect ratio.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) ellipse_aspect_ratio
  real ( kind = rk ) t

  a = abs ( a )
  b = abs ( b )

  if ( a < b ) then
    t = a
    a = b
    b = t
  end if

  ellipse_aspect_ratio = b / a

  return
end
function ellipse_eccentricity ( a, b )

!*****************************************************************************80
!
!! ellipse_eccentricity() computes the eccentricity of an ellipse.
!
!  Discussion:
!
!    The ellipse has major and minor semi-axes a and b.  In particular, it
!    could have the form:
!
!      (x/a)^2 + (y/b)^2 = 1
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 October 2022
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    John D Cook,
!    Eccentricity, Flattening, and Aspect Ratio,
!    https://www.johndcook.com/blog/2022/10/14/eccentricity-flatness-aspect/
!    Posted 14 October 2022.
!
!  Input:
!
!    real ( kind = rk ) A, B, the major and minor semi-axes.
!
!  Output:
!
!    real ( kind = rk ) ELLIPSE_ECCENTRICITY, the eccentricity.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) ellipse_eccentricity
  real ( kind = rk ) t

  a = abs ( a )
  b = abs ( b )

  if ( a < b ) then
    t = a
    a = b
    b = t
  end if

  ellipse_eccentricity = sqrt ( 1.0D+00 - ( b / a )**2 )

  return
end
function ellipse_flattening ( a, b )

!*****************************************************************************80
!
!! ellipse_flattening() computes the flattening of an ellipse.
!
!  Discussion:
!
!    The ellipse has major and minor semi-axes a and b.  In particular, it
!    could have the form:
!
!      (x/a)^2 + (y/b)^2 = 1
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 October 2022
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    John D Cook,
!    Eccentricity, Flattening, and Aspect Ratio,
!    https://www.johndcook.com/blog/2022/10/14/eccentricity-flatness-aspect/
!    Posted 14 October 2022.
!
!  Input:
!
!    real ( kind = rk ) A, B, the major and minor semi-axes.
!
!  Output:
!
!    real ( kind = rk ) ellipse_flattening: the flattening.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) ellipse_flattening
  real ( kind = rk ) t

  a = abs ( a )
  b = abs ( b )

  if ( a < b ) then
    t = a
    a = b
    b = t
  end if

  ellipse_flattening = ( a - b ) / a

  return
end
subroutine ellipse_point_dist_2d ( r1, r2, p, dist )

!*****************************************************************************80
!
!! ellipse_point_dist_2d() finds the distance from a point to an ellipse in 2D.
!
!  Discussion:
!
!    An ellipse in standard position has a center at the origin, and
!    axes aligned with the coordinate axes.  Any point P on the ellipse
!    satisfies
!
!      (  P(1) / R1 )^2 + ( P(2) / R2 )^2 == 1
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
!    Dianne O'Leary,
!    Elastoplastic Torsion: Twist and Stress,
!    Computing in Science and Engineering,
!    July/August 2004, pages 74-76.
!    September/October 2004, pages 63-65.
!
!  Input:
!
!    real ( kind = rk ) R1, R2, the ellipse parameters.  Normally,
!    these are both positive quantities.  Generally, they are also
!    distinct.
!
!    real ( kind = rk ) P(2), the point.
!
!  Output:
!
!    real ( kind = rk ) DIST, the distance to the ellipse.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) dist
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) pn(dim_num)
  real ( kind = rk ) r1
  real ( kind = rk ) r2

  call ellipse_point_near_2d ( r1, r2, p, pn )

  dist = sqrt ( sum ( ( p(1:dim_num) - pn(1:dim_num) )**2 ) )

  return
end
subroutine ellipse_point_near_2d ( r1, r2, p, pn )

!*****************************************************************************80
!
!! ellipse_point_near_2d() finds the nearest point on an ellipse in 2D.
!
!  Discussion:
!
!    An ellipse in standard position has a center at the origin, and
!    axes aligned with the coordinate axes.  Any point P on the ellipse
!    satisfies
!
!      (  P(1) / R1 )^2 + ( P(2) / R2 )^2 == 1
!
!    The nearest point PN on the ellipse has the property that the
!    line from PN to P is normal to the ellipse.  Points on the ellipse
!    can be parameterized by T, to have the form
!
!      ( R1 * cos ( T ), R2 * sin ( T ) ).
!
!    The tangent vector to the ellipse has the form
!
!      ( -R1 * sin ( T ), R2 * cos ( T ) ) 
!
!    At PN, the dot product of this vector with  ( P - PN ) must be
!    zero:
!
!      - R1 * sin ( T ) * ( X - R1 * cos ( T ) )
!      + R2 * cos ( T ) * ( Y - R2 * sin ( T ) ) = 0
!
!    This nonlinear equation for T can be solved by Newton's method.
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
!  Input:
!
!    real ( kind = rk ) R1, R2, the ellipse parameters.  Normally,
!    these are both positive quantities.  Generally, they are also
!    distinct.
!
!    real ( kind = rk ) P(2), the point.
!
!  Output:
!
!    real ( kind = rk ) PN(2), the point on the ellipse which
!    is closest to P.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) ct
  real ( kind = rk ) f
  real ( kind = rk ) fp
  integer iteration
  integer, parameter :: iteration_max = 100
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) pn(dim_num)
  real ( kind = rk ) r1
  real ( kind = rk ) r2
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) st
  real ( kind = rk ) t
  real ( kind = rk ) x
  real ( kind = rk ) y

  x = abs ( p(1) )
  y = abs ( p(2) )

  if ( y == 0.0D+00 .and. r1 * r1 - r2 * r2 <= r1 * x ) then

    t = 0.0D+00

  else if ( x == 0.0D+00 .and. r2 * r2 - r1 * r1 <= r2 * y ) then

    t = r8_pi / 2.0D+00

  else

    if ( y == 0.0D+00 ) then
      y = sqrt ( epsilon ( y ) ) * abs ( r2 )
    end if

    if ( x == 0.0D+00 ) then
      x = sqrt ( epsilon ( x ) ) * abs ( r1 )
    end if
!
!  Initial parameter T:
!
    t = atan2 ( y, x )

    iteration = 0

    do

      ct = cos ( t )
      st = sin ( t )

      f = ( x - abs ( r1 ) * ct ) * abs ( r1 ) * st &
        - ( y - abs ( r2 ) * st ) * abs ( r2 ) * ct

      if ( abs ( f ) <= 100.0D+00 * epsilon ( f ) ) then
        exit
      end if

      if ( iteration_max <= iteration ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'ELLIPSE_POINT_NEAR_2D(): Warning!'
        write ( *, '(a)' ) '  Reached iteration limit.'
        write ( *, '(a,f8.6)' ) '  T = ', t
        write ( *, '(a,g14.6)' ) '  F = ', f
        exit
      end if

      iteration = iteration + 1

      fp = r1 * r1 * st * st + r2 * r2 * ct * ct &
         + ( x - abs ( r1 ) * ct ) * abs ( r1 ) * ct &
         + ( y - abs ( r2 ) * st ) * abs ( r2 ) * st

      t = t - f / fp

    end do

  end if
!
!  From the T value, we get the nearest point.
!
  pn(1) = abs ( r1 ) * cos ( t )
  pn(2) = abs ( r2 ) * sin ( t )
!
!  Take care of case where the point was in another quadrant.
!
  pn(1) = sign ( 1.0D+00, p(1) ) * pn(1)
  pn(2) = sign ( 1.0D+00, p(2) ) * pn(2)

  return
end
subroutine ellipse_points_2d ( pc, r1, r2, psi, n, p )

!*****************************************************************************80
!
!! ellipse_points_2d() returns N points on an tilted ellipse in 2D.
!
!  Discussion:
!
!    An ellipse in standard position has a center at the origin, and
!    axes aligned with the coordinate axes.  Any point P on the ellipse
!    satisfies
!
!      (  P(1) / R1 )^2 + ( P(2) / R2 )^2 == 1
!
!    The points are "equally spaced" in the angular sense.  They are
!    not equally spaced along the perimeter of the ellipse.
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
!    real ( kind = rk ) PC(2), the center of the ellipse.
!
!    real ( kind = rk ) R1, R2, the "radius" of the ellipse in the major
!    and minor axis directions.  A circle has these values equal.
!
!    real ( kind = rk ) PSI, the angle that the major axis of the ellipse
!    makes with the X axis.  A value of 0.0 means that the major and
!    minor axes of the ellipse will be the X and Y coordinate axes.
!
!    integer N, the number of points desired.  N must 
!    be at least 1.
!
!  Output:
!
!    real ( kind = rk ) P(2,N), points on the ellipse.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer, parameter :: dim_num = 2

  integer i
  real ( kind = rk ) p(dim_num,n)
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) psi
  real ( kind = rk ) r1
  real ( kind = rk ) r2
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) theta

  do i = 1, n

    theta = ( 2.0D+00 * r8_pi * real ( i - 1, kind = rk ) ) &
      / real ( n, kind = rk )

    p(1,i) = pc(1) + r1 * cos ( psi ) * cos ( theta ) &
                   - r2 * sin ( psi ) * sin ( theta )

    p(2,i) = pc(2) + r1 * sin ( psi ) * cos ( theta ) &
                   + r2 * cos ( psi ) * sin ( theta )

  end do

  return
end
subroutine ellipse_points_arc_2d ( pc, r1, r2, psi, theta1, theta2, n, p )

!*****************************************************************************80
!
!! ellipse_points_arc_2d() returns N points on a tilted elliptical arc in 2D.
!
!  Discussion:
!
!    An ellipse in standard position has a center at the origin, and
!    axes aligned with the coordinate axes.  Any point P on the ellipse
!    satisfies
!
!      (  P(1) / R1 )^2 + ( P(2) / R2 )^2 == 1
!
!    The points are "equally spaced" in the angular sense.  They are
!    not equally spaced along the perimeter of the ellipse.
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
!    real ( kind = rk ) PC(2), the coordinates of the center of
!    the ellipse.
!
!    real ( kind = rk ) R1, R2, the "radius" of the ellipse in the major
!    and minor axis directions.  A circle has these values equal.
!
!    real ( kind = rk ) PSI, the angle that the major axis of the ellipse
!    makes with the X axis.  A value of 0.0 means that the major and
!    minor axes of the ellipse will be the X and Y coordinate axes.
!
!    real ( kind = rk ) THETA1, THETA2, the angular coordinates of
!    the first and last points to be drawn, in radians.  This angle is measured
!    with respect to the (possibly tilted) major axis.
!
!    integer N, the number of points desired.  N must 
!    be at least 1.
!
!  Output:
!
!    real ( kind = rk ) P(2,N), points on the ellipse.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: dim_num = 2

  real ( kind = rk ) r8_modp
  integer i
  real ( kind = rk ) p(dim_num,n)
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) psi
  real ( kind = rk ) r1
  real ( kind = rk ) r2
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
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
      theta = ( real ( n - i,     kind = rk ) * theta1 &
              + real (     i - 1, kind = rk ) * theta3 ) &
              / real ( n     - 1, kind = rk )
    else
      theta = 0.5D+00 * ( theta1 + theta3 )
    end if

    p(1,i) = pc(1) + r1 * cos ( psi ) * cos ( theta ) &
                   - r2 * sin ( psi ) * sin ( theta )

    p(2,i) = pc(2) + r1 * sin ( psi ) * cos ( theta ) &
                   + r2 * cos ( psi ) * sin ( theta )

  end do

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

