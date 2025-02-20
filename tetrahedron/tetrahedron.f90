function i4_gcd ( i, j )

!*****************************************************************************80
!
!! i4_gcd() finds the greatest common divisor of two I4's.
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
!  Input:
!
!    integer I, J, two numbers whose GCD is desired.
!
!  Output:
!
!    integer I4_GCD, the greatest common divisor
!    of I and J.
!
  implicit none

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
function i4_lcm ( i, j )

!*****************************************************************************80
!
!! i4_lcm() computes the least common multiple of two I4's.
!
!  Discussion:
!
!    The least common multiple may be defined as
!
!      LCM(I,J) = ABS ( I * J ) / GCD(I,J)
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
!  Input:
!
!    integer I, J, the integers whose I4_LCM is desired.
!
!  Output:
!
!    integer I4_LCM, the least common multiple of I and J.
!    I4_LCM is never negative.  I4_LCM is 0 if either I or J is zero.
!
  implicit none

  integer i
  integer i4_gcd
  integer j
  integer i4_lcm

  i4_lcm = abs ( i * ( j / i4_gcd ( i, j ) ) )

  return
end
function i4vec_lcm ( n, v )

!*****************************************************************************80
!
!! i4vec_lcm() returns the least common multiple of an I4VEC.
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
!    31 July 2009
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the order of V.
!
!    integer V(N), the vector.
!
!  Output:
!
!    integer I4VEC_LCM, the least common multiple of V.
!
  implicit none

  integer n

  integer i
  integer i4_lcm
  integer i4vec_lcm
  integer v(n)
  integer value

  value = 1

  do i = 1, n

    if ( v(i) == 0 ) then
      value = 0
      exit
    end if

    value = i4_lcm ( value, v(i) )

  end do

  i4vec_lcm = value

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
!  Input:
!
!    integer N, the number of vertices.
!
!    real ( kind = rk ) V(3,N), the coordinates of the vertices.
!    The vertices should be listed in neighboring order.
!
!  Output:
!
!    real ( kind = rk ) AREA, the area of the polygon.
!
!    real ( kind = rk ) NORMAL(3), the unit normal vector to the polygon.
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
!  Input:
!
!    real ( kind = rk ) C, the argument.
!
!  Output:
!
!    real ( kind = rk ) R8_ACOS, an angle whose cosine is C.
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
subroutine r8_swap ( x, y )

!*****************************************************************************80
!
!! r8_swap() swaps two R8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) X, Y, the values to be interchanged.
!
!  Output:
!
!    real ( kind = rk ) X, Y, the interchanged values.
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
function r8mat_det_4d ( a )

!*****************************************************************************80
!
!! r8mat_det_4d() computes the determinant of a 4 by 4 matrix.
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
!  Input:
!
!    real ( kind = rk ) A(4,4), the matrix whose determinant is desired.
!
!  Output:
!
!    real ( kind = rk ) R8MAT_DET_4D, the determinant of the matrix.
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
subroutine r8mat_solve ( n, rhs_num, a, info )

!*****************************************************************************80
!
!! r8mat_solve() uses Gauss-Jordan elimination to solve an N by N linear system.
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
!  Input:
!
!    integer N, the order of the matrix.
!
!    integer RHS_NUM, the number of right hand sides.
!    RHS_NUM must be at least 0.
!
!    real ( kind = rk ) A(N,N+rhs_num), contains in rows and
!    columns 1 to N the coefficient matrix, and in columns N+1 through
!    N+rhs_num, the right hand sides.
!
!  Output:
!
!    integer INFO, singularity flag.
!    0, the matrix was not singular, the solutions were computed;
!    J, factorization failed on step J, and the solutions could not
!    be computed.

!    real ( kind = rk ) A(N,N+rhs_num), the coefficient matrix
!    area has been destroyed, while the right hand sides have
!    been overwritten with the corresponding solutions.
!
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
subroutine r8mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! r8mat_transpose_print() prints an R8MAT, transposed.
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
!    character ( len = * ) TITLE, an optional title.
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
!    An R8MAT is an array of R8 values.
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
!    integer ILO, JLO, the first row and column to print.
!
!    integer IHI, JHI, the last row and column to print.
!
!    character ( len = * ) TITLE, an optional title.
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

      write ( *, '(i5,1x,5a14)' ) j, ( ctemp(i), i = 1, inc )

    end do

  end do

  return
end
subroutine r8vec_angle_3d ( u, v, angle )

!*****************************************************************************80
!
!! r8vec_angle_3d() computes the angle between two vectors in 3D.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 July 2009
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) U(3), V(3), the vectors.
!
!  Output:
!
!    real ( kind = rk ) ANGLE, the angle between the two vectors.
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
subroutine r8vec_cross_product_3d ( v1, v2, v3 )

!*****************************************************************************80
!
!! r8vec_cross_product_3d() computes the cross product of two vectors in 3D.
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
!  Input:
!
!    real ( kind = rk ) V1(3), V2(3), the two vectors.
!
!  Output:
!
!    real ( kind = rk ) V3(3), the cross product vector.
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
function r8vec_length ( dim_num, x )

!*****************************************************************************80
!
!! r8vec_length() returns the Euclidean length of a vector.
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
!  Input:
!
!    integer DIM_NUM, the spatial dimension.
!
!    real ( kind = rk ) X(DIM_NUM), the vector.
!
!  Output:
!
!    real ( kind = rk ) R8VEC_LENGTH, the Euclidean length of the vector.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num

  real ( kind = rk ) r8vec_length
  real ( kind = rk ) x(dim_num)

  r8vec_length = sqrt ( sum ( ( x(1:dim_num) )**2 ) )

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
!    character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  integer i
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,2x,g16.8)' ) i, a(i)
  end do

  return
end
subroutine r8vec_transpose_print ( n, a, title )

!*****************************************************************************80
!
!! r8vec_transpose_print() prints an R8VEC "transposed".
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
  integer ihi
  integer ilo
  character ( len = * ) title
  integer title_length

  title_length = len_trim ( title )

  do ilo = 1, n, 5
    if ( ilo == 1 ) then
      write ( *, '(a)', advance = 'NO' ) trim ( title )
    else
      write ( *, '(a)', advance = 'NO' ) ( ' ', i = 1, title_length )
    end if
    write ( *, '(2x)', advance = 'NO' )
    ihi = min ( ilo + 5 - 1, n )
    write ( *, '(5g14.6)' ) a(ilo:ihi)
  end do

  return
end
subroutine shape_print ( point_num, face_num, face_order_max, &
  point_coord, face_order, face_point )

!*****************************************************************************80
!
!! shape_print() prints information about a polyhedron.
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
!  Input:
!
!    integer POINT_NUM, the number of points.
!
!    integer FACE_NUM, the number of faces.
!
!    integer FACE_ORDER_MAX, the number of vertices 
!    per face.
!
!    real ( kind = rk ) POINT_COORD(3,POINT_NUM), the vertices.
!
!    integer FACE_ORDER(FACE_NUM), the number of vertices
!    per face.
!
!    integer FACE_POINT(FACE_ORDER_MAX,FACE_NUM); 
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
  write ( *, '(a)' ) 'shape_print():'
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
subroutine tetrahedron_barycentric ( tetra, p, c )

!*****************************************************************************80
!
!! tetrahedron_barycentric(): barycentric coordinates of a point in 3D.
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
!  Input:
!
!    real ( kind = rk ) TETRA(3,4) the tetrahedron vertices.
!
!    real ( kind = rk ) P(3), the point to be checked.
!
!  Output:
!
!    real ( kind = rk ) C(4), the barycentric coordinates of P with
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
    write ( *, '(a)' ) 'tetrahedron_barycentric(): Fatal error!'
    write ( *, '(a)' ) '  The linear system is singular.'
    write ( *, '(a)' ) '  The input data does not form a proper tetrahedron.'
    stop 1
  end if

  c(2:4) = a(1:dim_num,4)

  c(1) = 1.0D+00 - sum ( c(2:4) )

  return
end
subroutine tetrahedron_centroid ( tetra, centroid )

!*****************************************************************************80
!
!! tetrahedron_centroid() computes the centroid of a tetrahedron.
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
!  Input:
!
!    real ( kind = rk ) TETRA(3,4) the tetrahedron vertices.
!
!  Output:
!
!    real ( kind = rk ) CENTROID(3), the coordinates of the centroid.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) centroid(3)
  integer i
  real ( kind = rk ) tetra(3,4)

  do i = 1, 3
    centroid(i) = sum ( tetra(i,1:4) ) / 4.0D+00
  end do

  return
end
subroutine tetrahedron_circumsphere ( tetra, r, pc )

!*****************************************************************************80
!
!! tetrahedron_circumsphere() computes the circumsphere of a tetrahedron.
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
!    Butterworths, 1983.
!
!  Input:
!
!    real ( kind = rk ) TETRA(3,4) the tetrahedron vertices.
!
!  Output:
!
!    real ( kind = rk ) R, PC(3), the center of the
!    circumscribed sphere, and its radius.  If the linear system is
!    singular, then R = -1, PC(1:3) = 0.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(3,4)
  integer i
  integer info
  integer j
  real ( kind = rk ) pc(3)
  real ( kind = rk ) r
  real ( kind = rk ) tetra(3,4)
!
!  Set up the linear system.
!
  a(1:3,1:3) = transpose ( tetra(1:3,2:4) )

  do j = 1, 3
    a(1:3,j) = a(1:3,j) - tetra(j,1)
  end do

  do i = 1, 3
    a(i,4) = sum ( a(i,1:3)**2 )
  end do
!
!  Solve the linear system.
!
  call r8mat_solve ( 3, 1, a, info )
!
!  If the system was singular, return a consolation prize.
!
  if ( info /= 0 ) then
    r = -1.0D+00
    pc(1:3) = 0.0D+00
    return
  end if
!
!  Compute the radius and center.
!
  r = 0.5D+00 * sqrt ( sum ( a(1:3,4)**2 ) )

  pc(1:3) = tetra(1:3,1) + 0.5D+00 * a(1:3,4)

  return
end
subroutine tetrahedron_contains_point ( tetra, p, inside )

!*****************************************************************************80
!
!! tetrahedron_contains_point() finds if a point is inside a tetrahedron in 3D.
!
!  Discussion:
!
!    We compute the barycentric coordinates C(1:4) of the point, with respect
!    to the tetrahedron.  The point is inside the tetrahedron if and only
!    if each coordinate is nonnegative, and their sum is no greater than 1.
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
!  Input:
!
!    real ( kind = rk ) TETRA(3,4), the tetrahedron vertices.
!
!    real ( kind = rk ) P(3), the point to be checked.
!
!  Output:
!
!    logical INSIDE, is TRUE if P is inside the 
!    tetrahedron.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) c(dim_num+1)
  logical inside
  real ( kind = rk ) p(dim_num)
  real ( kind = rk ) tetra(dim_num,4)

  call tetrahedron_barycentric ( tetra, p, c )
!
!  If the point is in the tetrahedron, its barycentric coordinates
!  must be nonnegative.
!
  if ( any ( c(1:dim_num+1) < 0.0D+00 ) ) then
    inside = .false.
  else
    inside = .true.
  end if

  return
end
subroutine tetrahedron_dihedral_angles ( tetra, angle )

!*****************************************************************************80
!
!! tetrahedron_dihedral_angles() computes dihedral angles of a tetrahedron.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 July 2009
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) TETRA(3,4), the vertices of the tetrahedron.
!
!  Output:
!
!    real ( kind = rk ) ANGLE(6), the dihedral angles along the
!    axes AB, AC, AD, BC, BD and CD, respectively.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) ab(3)
  real ( kind = rk ) abc_normal(3)
  real ( kind = rk ) abd_normal(3)
  real ( kind = rk ) ac(3)
  real ( kind = rk ) acd_normal(3)
  real ( kind = rk ) ad(3)
  real ( kind = rk ) angle(6)
  real ( kind = rk ) bc(3)
  real ( kind = rk ) bcd_normal(3)
  real ( kind = rk ) bd(3)
  real ( kind = rk ) cd(3)
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) tetra(3,4)

  call tetrahedron_edges ( tetra, ab, ac, ad, bc, bd, cd )

  call r8vec_cross_product_3d ( ac, ab, abc_normal )
  call r8vec_cross_product_3d ( ab, ad, abd_normal )
  call r8vec_cross_product_3d ( ad, ac, acd_normal )
  call r8vec_cross_product_3d ( bc, bd, bcd_normal )

  call r8vec_angle_3d ( abc_normal, abd_normal, angle(1) )
  call r8vec_angle_3d ( abc_normal, acd_normal, angle(2) )
  call r8vec_angle_3d ( abd_normal, acd_normal, angle(3) )
  call r8vec_angle_3d ( abc_normal, bcd_normal, angle(4) )
  call r8vec_angle_3d ( abd_normal, bcd_normal, angle(5) )
  call r8vec_angle_3d ( acd_normal, bcd_normal, angle(6) )

  angle(1:6) = r8_pi - angle(1:6)

  return
end
subroutine tetrahedron_edges ( tetra, ab, ac, ad, bc, bd, cd )

!*****************************************************************************80
!
!! tetrahedron_edges() computes the edges of a tetrahedron.
!
!  Discussion:
!
!    The vertices are A, B, C, D.  The edge from A to B is denoted by AB.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 May 2014
!
!  Author:
!
!    John Burkardt.
!
!  Input:
!
!    real ( kind = rk ) TETRA(3,4), the vertices of the tetrahedron.
!
!  Output:
!
!    real ( kind = rk ) AB(3), AC(3), AD(3), BC(3), BD(3), CD(3), 
!    vectors that represent the edges of the tetrahedron.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) ab(3)
  real ( kind = rk ) ac(3)
  real ( kind = rk ) ad(3)
  real ( kind = rk ) bc(3)
  real ( kind = rk ) bd(3)
  real ( kind = rk ) cd(3)
  real ( kind = rk ) tetra(3,4)

  ab(1:3) = tetra(1:3,2) - tetra(1:3,1)
  ac(1:3) = tetra(1:3,3) - tetra(1:3,1)
  ad(1:3) = tetra(1:3,4) - tetra(1:3,1)
  bc(1:3) = tetra(1:3,3) - tetra(1:3,2)
  bd(1:3) = tetra(1:3,4) - tetra(1:3,2)
  cd(1:3) = tetra(1:3,4) - tetra(1:3,3)

  return
end
subroutine tetrahedron_edge_length ( tetra, edge_length )

!*****************************************************************************80
!
!! tetrahedron_edge_length() returns edge lengths of a tetrahedron.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) TETRA(3,4), the tetrahedron vertices.
!
!  Output:
!
!    real ( kind = rk ) EDGE_LENGTH(6), the length of the edges.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r8vec_length
  real ( kind = rk ) edge_length(6)
  integer j1
  integer j2
  integer k
  real ( kind = rk ) tetra(3,4)

  k = 0
  do j1 = 1, 3
    do j2 = j1 + 1, 4
      k = k + 1
      edge_length(k) = r8vec_length ( 3, tetra(1:3,j2) - tetra(1:3,j1) )
     end do
  end do

  return
end
subroutine tetrahedron_face_angles ( tetra, angles )

!*****************************************************************************80
!
!! tetrahedron_face_angles() returns the 12 face angles of a tetrahedron.
!
!  Discussion:
!
!    The tetrahedron has 4 triangular faces.  This routine computes the
!    3 planar angles associated with each face.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 July 2009
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) TETRA(3,4) the tetrahedron vertices.
!
!  Output:
!
!    real ( kind = rk ) ANGLES(3,4), the face angles.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) angles(3,4)
  real ( kind = rk ) tri(3,3)
  real ( kind = rk ) tetra(3,4)
!
!  Face 123
!
  tri(1:3,1:3) = tetra(1:3,1:3)
  call triangle_angles_3d ( tri, angles(1:3,1) )
!
!  Face 124
!
  tri(1:3,1:2) = tetra(1:3,1:2)
  tri(1:3,3) = tetra(1:3,4)
  call triangle_angles_3d ( tri, angles(1:3,2) )
!
!  Face 134
!
  tri(1:3,1) = tetra(1:3,1)
  tri(1:3,2:3) = tetra(1:3,3:4)
  call triangle_angles_3d ( tri, angles(1:3,3) )
!
!  Face 234
!
  tri(1:3,1:3) = tetra(1:3,2:4)
  call triangle_angles_3d ( tri, angles(1:3,4) )

  return
end
subroutine tetrahedron_face_areas ( tetra, areas )

!*****************************************************************************80
!
!! tetrahedron_face_areas() returns the 4 face areas of a tetrahedron.
!
!  Discussion:
!
!    The tetrahedron has 4 triangular faces.  This routine computes the
!    area of each face.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 July 2009
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) TETRA(3,4) the tetrahedron vertices.
!
!  Output:
!
!    real ( kind = rk ) AREAS(4), the face areas.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) areas(4)
  real ( kind = rk ) tri(3,3)
  real ( kind = rk ) tetra(3,4)
!
!  Face 123
!
  tri(1:3,1:3) = tetra(1:3,1:3)
  call triangle_area_3d ( tri, areas(1) )
!
!  Face 124
!
  tri(1:3,1:2) = tetra(1:3,1:2)
  tri(1:3,3) = tetra(1:3,4)
  call triangle_area_3d ( tri, areas(2) )
!
!  Face 134
!
  tri(1:3,1) = tetra(1:3,1)
  tri(1:3,2:3) = tetra(1:3,3:4)
  call triangle_area_3d ( tri, areas(3) )
!
!  Face 234
!
  tri(1:3,1:3) = tetra(1:3,2:4)
  call triangle_area_3d ( tri, areas(4) )

  return
end
subroutine tetrahedron_insphere ( tetra, r, pc )

!*****************************************************************************80
!
!! tetrahedron_insphere() finds the insphere of a tetrahedron.
!
!  Discussion:
!
!    The insphere of a tetrahedron is the inscribed sphere, which touches
!    each face of the tetrahedron at a single point.
!
!    The points of contact are the centroids of the triangular faces
!    of the tetrahedron.  Therefore, the point of contact for a face
!    can be computed as the average of the vertices of that face.
!
!    The sphere can then be determined as the unique sphere through
!    the four given centroids.
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
!    Philip Schneider, David Eberly,
!    Geometric Tools for Computer Graphics,
!    Elsevier, 2002,
!    ISBN: 1558605940,
!    LC: T385.G6974.
!
!  Input:
!
!    real ( kind = rk ) TETRA(3,4), the vertices of the tetrahedron.
!
!  Output:
!
!    real ( kind = rk ) R, PC(3), the radius and the center
!    of the sphere.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) b(4,4)
  real ( kind = rk ) r8mat_det_4d
  real ( kind = rk ) r8vec_length
  real ( kind = rk ) gamma
  real ( kind = rk ) l123
  real ( kind = rk ) l124
  real ( kind = rk ) l134
  real ( kind = rk ) l234
  real ( kind = rk ) n123(1:3)
  real ( kind = rk ) n124(1:3)
  real ( kind = rk ) n134(1:3)
  real ( kind = rk ) n234(1:3)
  real ( kind = rk ) pc(1:3)
  real ( kind = rk ) r
  real ( kind = rk ) tetra(1:3,4)
  real ( kind = rk ) v21(1:3)
  real ( kind = rk ) v31(1:3)
  real ( kind = rk ) v41(1:3)
  real ( kind = rk ) v32(1:3)
  real ( kind = rk ) v42(1:3)
  real ( kind = rk ) v43(1:3)

  call tetrahedron_edges ( tetra, v21, v31, v41, v32, v42, v43 )

  call r8vec_cross_product_3d ( v21, v31, n123 )
  call r8vec_cross_product_3d ( v41, v21, n124 )
  call r8vec_cross_product_3d ( v31, v41, n134 )
  call r8vec_cross_product_3d ( v42, v32, n234 )

  l123 = r8vec_length ( 3, n123 )
  l124 = r8vec_length ( 3, n124 )
  l134 = r8vec_length ( 3, n134 )
  l234 = r8vec_length ( 3, n234 )

  pc(1:3) = ( l234 * tetra(1:3,1)   &
            + l134 * tetra(1:3,2)   &
            + l124 * tetra(1:3,3)   &
            + l123 * tetra(1:3,4) ) &
            / ( l234 + l134 + l124 + l123 )

  b(1:3,1:4) = tetra(1:3,1:4)
  b(4,1:4) = 1.0D+00

  gamma = abs ( r8mat_det_4d ( b ) )

  r = gamma / ( l234 + l134 + l124 + l123 )

  return
end
subroutine tetrahedron_lattice_layer_point_next ( c, v, more )

!*****************************************************************************80
!
!! tetrahedron_lattice_layer_point_next(): next tetrahedron lattice layer point.
!
!  Discussion:
!
!    The tetrahedron lattice layer L is bounded by the lines
!
!      0 <= X,
!      0 <= Y,
!      0 <= Z,
!      L - 1 < X / C(1) + Y / C(2) + Z/C(3) <= L.
!
!    In particular, layer L = 0 always contains the single point (0,0).
!
!    This function returns, one at a time, the points that lie within 
!    a given tetrahedron lattice layer.
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
!  Input:
!
!    integer C(4), coefficients defining the 
!    lattice layer in entries 1 to 3, and the laver index in C(4).  
!    The coefficients should be positive, and C(4) must be nonnegative.
!
!    integer V(3).  On first call for a given layer,
!    the input value of V is not important.  On a repeated call for the same
!    layer, the input value of V should be the output value from the previous 
!    call.
!
!    logical MORE, set MORE to FALSE to 
!    indicate that this is the first call for a given layer.  Thereafter, the 
!    input value should be the output value from the previous call.
!
!  Output:
!
!    integer V(3), the next lattice layer point.
!
!    logical MORE, TRUE if the returned value V is a new point.
!    If the output value is FALSE, then no more points were found,
!    and V was reset to 0, and the lattice layer has been exhausted.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer c(4)
  integer c1n
  integer i4vec_lcm
  integer lhs
  logical more
  integer, parameter :: n = 3
  integer rhs1
  integer rhs2
  integer v(3)
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
!  Can we simply increase X?
!
    v(1) = v(1) + 1

    lhs = ( c1n / c(1) ) * v(1) &
        + ( c1n / c(2) ) * v(2) &
        + ( c1n / c(3) ) * v(3)

    if ( lhs <= rhs2 ) then
!
!  No.  Increase Y, and set X so we just exceed RHS1...if possible.
!
    else

      v(2) = v(2) + 1

      v(1) = ( c(1) * ( rhs1 - ( c1n / c(2) ) * v(2) &
                             - ( c1n / c(3) ) * v(3) ) ) / c1n
      v(1) = max ( v(1), 0 )

      lhs = ( c1n / c(1) ) * v(1) &
          + ( c1n / c(2) ) * v(2) &
          + ( c1n / c(3) ) * v(3)

      if ( lhs <= rhs1 ) then
        v(1) = v(1) + 1
        lhs = lhs + c1n / c(1)
      end if
!
!  We have increased Y by 1.  Have we stayed below the upper bound?
!
      if ( lhs <= rhs2 ) then

      else
!
!  No.  Increase Z, and set X so we just exceed RHS1...if possible.
!
        v(3) = v(3) + 1
        v(2) = 0
        v(1) = ( c(1) * ( rhs1 - ( c1n / c(2) ) * v(2) &
                               - ( c1n / c(3) ) * v(3) ) ) / c1n
        v(1) = max ( v(1), 0 )

        lhs = ( c1n / c(1) ) * v(1) &
            + ( c1n / c(2) ) * v(2) &
            + ( c1n / c(3) ) * v(3)

        if ( lhs <= rhs1 ) then
          v(1) = v(1) + 1
          lhs = lhs + c1n / c(1)
        end if

        if ( lhs <= rhs2 ) then

        else
          more = .false.
          v(1:n) = 0
        end if

      end if
    end if
  end if

  return
end
subroutine tetrahedron_lattice_point_next ( c, v, more )

!*****************************************************************************80
!
!! tetrahedron_lattice_point_next() returns the next tetrahedron lattice point.
!
!  Discussion:
!
!    The lattice tetrahedron is defined by the vertices:
!
!      (0,0,0), (C(4)/C(1),0,0), (0,C(4)/C(2),0) and (0,0,C(4)/C(3))
!
!    The lattice tetrahedron is bounded by the lines
!
!      0 <= X,
!      0 <= Y
!      0 <= Z,
!      X / C(1) + Y / C(2) + Z / C(3) <= C(4)
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
!  Input:
!
!    integer C(4), coefficients defining the 
!    lattice tetrahedron.  These should be positive.
!
!    integer V(3).  On first call, the input
!    value is not important.  On a repeated call, the input value should
!    be the output value from the previous call.
!
!    logical MORE, TRUE if not only is the returned value V a lattice point,
!    but the routine can be called again for another lattice point.
!    If the output value is FALSE, then no more lattice points were found,
!    and V was reset to 0, and the routine should not be called further
!    for this tetrahedron.
!
!  Output:
!
!    integer V(3), the next lattice point.
!
!    logical MORE, TRUE if not only is the returned value V a lattice point,
!    but the routine can be called again for another lattice point.
!    If the output value is FALSE, then no more lattice points were found,
!    and V was reset to 0, and the routine should not be called further
!    for this tetrahedron.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer c(4)
  integer c1n
  integer i4vec_lcm
  integer lhs
  logical more
  integer, parameter :: n = 3
  integer rhs
  integer v(3)

  if ( .not. more ) then

    v(1:n) = 0
    more = .true.

  else

    c1n = i4vec_lcm ( n, c )

    rhs = c1n * c(n+1)

    lhs =        c(2) * c(3) * v(1) &
        + c(1) *        c(3) * v(2) &
        + c(1) * c(2)        * v(3)

    if ( lhs +  c1n / c(1) <= rhs ) then

      v(1) = v(1) + 1

    else

      lhs = lhs - c1n * v(1) / c(1)
      v(1) = 0

      if ( lhs + c1n / c(2) <= rhs ) then

        v(2) = v(2) + 1

      else

        lhs = lhs - c1n * v(2) / c(2)
        v(2) = 0

        if ( lhs + c1n / c(3) <= rhs ) then

          v(3) = v(3) + 1

        else

          v(3) = 0
          more = .false.

        end if

      end if

    end if
      
  end if

  return
end
subroutine tetrahedron_quality1 ( tetra, quality )

!*****************************************************************************80
!
!! tetrahedron_quality1(): "quality" of a tetrahedron.
!
!  Discussion:
!
!    The quality of a tetrahedron is 3 times the ratio of the radius of
!    the inscribed sphere divided by that of the circumscribed sphere.
!
!    An equilateral tetrahredron achieves the maximum possible quality of 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) TETRA(3,4), the tetrahedron vertices.
!
!  Output:
!
!    real ( kind = rk ) QUALITY, the quality of the tetrahedron.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) pc(3)
  real ( kind = rk ) quality
  real ( kind = rk ) r_in
  real ( kind = rk ) r_out
  real ( kind = rk ) tetra(3,4)

  call tetrahedron_circumsphere ( tetra, r_out, pc )

  call tetrahedron_insphere ( tetra, r_in, pc )

  quality = 3.0D+00 * r_in / r_out

  return
end
subroutine tetrahedron_quality2 ( tetra, quality2 )

!*****************************************************************************80
!
!! tetrahedron_quality2(): "quality" of a tetrahedron.
!
!  Discussion:
!
!    The quality measure #2 of a tetrahedron is:
!
!      QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
!
!    where
!
!      RIN = radius of the inscribed sphere;
!      LMAX = length of longest side of the tetrahedron.
!
!    An equilateral tetrahredron achieves the maximum possible quality of 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Qiang Du, Desheng Wang,
!    The Optimal Centroidal Voronoi Tesselations and the Gersho's
!    Conjecture in the Three-Dimensional Space,
!    Computers and Mathematics with Applications,
!    Volume 49, 2005, pages 1355-1373.
!
!  Input:
!
!    real ( kind = rk ) TETRA(3,4), the tetrahedron vertices.
!
!  Output:
!
!    real ( kind = rk ) QUALITY2, the quality of the tetrahedron.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) edge_length(6)
  real ( kind = rk ) l_max
  real ( kind = rk ) pc(3)
  real ( kind = rk ) quality2
  real ( kind = rk ) r_in
  real ( kind = rk ) tetra(3,4)

  call tetrahedron_edge_length ( tetra, edge_length )

  l_max = maxval ( edge_length(1:6) )

  call tetrahedron_insphere ( tetra, r_in, pc )

  quality2 = 2.0D+00 * sqrt ( 6.0D+00 ) * r_in / l_max

  return
end
subroutine tetrahedron_quality3 ( tetra, quality3 )

!*****************************************************************************80
!
!! tetrahedron_quality3() computes the mean ratio of a tetrahedron.
!
!  Discussion:
!
!    This routine computes QUALITY3, the eigenvalue or mean ratio of
!    a tetrahedron.
!
!      QUALITY3 = 12 * ( 3 * volume )^(2/3) / (sum of squares of edge lengths).
!
!    This value may be used as a shape quality measure for the tetrahedron.
!
!    For an equilateral tetrahedron, the value of this quality measure
!    will be 1.  For any other tetrahedron, the value will be between
!    0 and 1.
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
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Input:
!
!    real ( kind = rk ) TETRA(3,4), the vertices of the tetrahedron.
!
!  Output:
!
!    real ( kind = rk ) QUALITY3, the mean ratio of the tetrahedron.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) ab(3)
  real ( kind = rk ) ac(3)
  real ( kind = rk ) ad(3)
  real ( kind = rk ) bc(3)
  real ( kind = rk ) bd(3)
  real ( kind = rk ) cd(3)
  real ( kind = rk ) denom
  real ( kind = rk ) lab
  real ( kind = rk ) lac
  real ( kind = rk ) lad
  real ( kind = rk ) lbc
  real ( kind = rk ) lbd
  real ( kind = rk ) lcd
  real ( kind = rk ) quality3
  real ( kind = rk ) tetra(3,4)
  real ( kind = rk ) volume
!
!  Compute the vectors representing the sides of the tetrahedron.
!
  call tetrahedron_edges ( tetra, ab, ac, ad, bc, bd, cd )
!
!  Compute the squares of the lengths of the sides.
!
  lab = sum ( ab(1:3)**2 )
  lac = sum ( ac(1:3)**2 )
  lad = sum ( ad(1:3)**2 )
  lbc = sum ( bc(1:3)**2 )
  lbd = sum ( bd(1:3)**2 )
  lcd = sum ( cd(1:3)**2 )
!
!  Compute the volume.
!
  volume = abs ( &
      ab(1) * ( ac(2) * ad(3) - ac(3) * ad(2) ) &
    + ab(2) * ( ac(3) * ad(1) - ac(1) * ad(3) ) &
    + ab(3) * ( ac(1) * ad(2) - ac(2) * ad(1) ) ) / 6.0D+00

  denom = lab + lac + lad + lbc + lbd + lcd

  if ( denom == 0.0D+00 ) then
    quality3 = 0.0D+00
  else
    quality3 = 12.0D+00 * ( 3.0D+00 * volume )**( 2.0D+00 / 3.0D+00 ) / denom
  end if

  return
end
subroutine tetrahedron_quality4 ( tetra, quality4 )

!*****************************************************************************80
!
!! tetrahedron_quality4() computes the minimum solid angle of a tetrahedron.
!
!  Discussion:
!
!    This routine computes a quality measure for a tetrahedron, based
!    on the sine of half the minimum of the four solid angles.
!
!    The quality measure for an equilateral tetrahedron should be 1,
!    since the solid angles of such a tetrahedron are each equal to pi.
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
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Input:
!
!    real ( kind = rk ) TETRA(3,4), the vertices of the tetrahedron.
!
!  Output:
!
!    real ( kind = rk ) QUALITY4, the value of the quality measure.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) ab(3)
  real ( kind = rk ) ac(3)
  real ( kind = rk ) ad(3)
  real ( kind = rk ) bc(3)
  real ( kind = rk ) bd(3)
  real ( kind = rk ) cd(3)
  real ( kind = rk ) denom
  real ( kind = rk ) l1
  real ( kind = rk ) l2
  real ( kind = rk ) l3
  real ( kind = rk ) lab
  real ( kind = rk ) lac
  real ( kind = rk ) lad
  real ( kind = rk ) lbc
  real ( kind = rk ) lbd
  real ( kind = rk ) lcd
  real ( kind = rk ) quality4
  real ( kind = rk ) tetra(3,4)
  real ( kind = rk ) volume
!
!  Compute the vectors that represent the sides.
!
  call tetrahedron_edges ( tetra, ab, ac, ad, bc, bd, cd )
!
!  Compute the lengths of the sides.
!
  lab = sqrt ( sum ( ab(1:3)**2 ) )
  lac = sqrt ( sum ( ac(1:3)**2 ) )
  lad = sqrt ( sum ( ad(1:3)**2 ) )
  lbc = sqrt ( sum ( bc(1:3)**2 ) )
  lbd = sqrt ( sum ( bd(1:3)**2 ) )
  lcd = sqrt ( sum ( cd(1:3)**2 ) )
!
!  Compute the volume
!
  volume = abs ( &
      ab(1) * ( ac(2) * ad(3) - ac(3) * ad(2) ) &
    + ab(2) * ( ac(3) * ad(1) - ac(1) * ad(3) ) &
    + ab(3) * ( ac(1) * ad(2) - ac(2) * ad(1) ) ) / 6.0D+00

  quality4 = 1.0D+00

  l1 = lab + lac
  l2 = lab + lad
  l3 = lac + lad

  denom = ( l1 + lbc ) * ( l1 - lbc ) &
        * ( l2 + lbd ) * ( l2 - lbd ) &
        * ( l3 + lcd ) * ( l3 - lcd )

  if ( denom <= 0.0D+00 ) then
    quality4 = 0.0D+00
  else
    quality4 = min ( quality4, 12.0D+00 * volume / sqrt ( denom ) )
  end if

  l1 = lab + lbc
  l2 = lab + lbd
  l3 = lbc + lbd

  denom = ( l1 + lac ) * ( l1 - lac ) &
        * ( l2 + lad ) * ( l2 - lad ) &
        * ( l3 + lcd ) * ( l3 - lcd )

  if ( denom <= 0.0D+00 ) then
    quality4 = 0.0D+00
  else
    quality4 = min ( quality4, 12.0D+00 * volume / sqrt ( denom ) )
  end if

  l1 = lac + lbc
  l2 = lac + lcd
  l3 = lbc + lcd

  denom = ( l1 + lab ) * ( l1 - lab ) &
        * ( l2 + lad ) * ( l2 - lad ) &
        * ( l3 + lbd ) * ( l3 - lbd )

  if ( denom <= 0.0D+00 ) then
    quality4 = 0.0D+00
  else
    quality4 = min ( quality4, 12.0D+00 * volume / sqrt ( denom ) )
  end if

  l1 = lad + lbd
  l2 = lad + lcd
  l3 = lbd + lcd

  denom = ( l1 + lab ) * ( l1 - lab ) &
        * ( l2 + lac ) * ( l2 - lac ) &
        * ( l3 + lbc ) * ( l3 - lbc )

  if ( denom <= 0.0D+00 ) then
    quality4 = 0.0D+00
  else
    quality4 = min ( quality4, 12.0D+00 * volume / sqrt ( denom ) )
  end if

  quality4 = quality4 * 1.5D+00 * sqrt ( 6.0D+00 )

  return
end
subroutine tetrahedron_rhombic_shape ( point_num, face_num, &
  face_order_max, point_coord, face_order, face_point )

!*****************************************************************************80
!
!! tetrahedron_rhombic_shape() describes a rhombic tetrahedron in 3D.
!
!  Discussion:
!
!    Call TETRAHEDRON_RHOMBIC_SIZE() first, to get dimension information.
!
!    The tetrahedron is described using 10 nodes.  If we label the vertices
!    P0, P1, P2 and P3, then the extra nodes lie halfway between vertices,
!    and have the labels P01, P02, P03, P12, P13 and P23.
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
!  Reference:
!
!    Anwei Liu, Barry Joe,
!    Quality Local Refinement of Tetrahedral Meshes Based
!    on 8-Subtetrahedron Subdivision,
!    Mathematics of Computation,
!    Volume 65, Number 215, July 1996, pages 1183-1200.
!
!  Input:
!
!    integer POINT_NUM, the number of points.
!
!    integer FACE_NUM, the number of faces.
!
!    integer FACE_ORDER_MAX, the maximum number of 
!    vertices per face.
!
!  Output:
!
!    real ( kind = rk ) POINT_COORD(3,POINT_NUM), the vertices.
!
!    integer FACE_ORDER(FACE_NUM), the number of vertices
!    for each face.
!
!    integer FACE_POINT(FACE_ORDER_MAX,FACE_NUM); 
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
  real ( kind = rk ) d
  integer face_order(face_num)
  integer face_point(face_order_max,face_num)
  real ( kind = rk ) point_coord(dim_num,point_num)
  real ( kind = rk ), parameter :: z = 0.0D+00

  a =        1.0D+00   / sqrt ( 3.0D+00 )
  b = sqrt ( 2.0D+00 ) / sqrt ( 3.0D+00 )
  c = sqrt ( 3.0D+00 ) /        6.0D+00
  d =        1.0D+00   / sqrt ( 6.0D+00 )
!
!  Set the point coordinates.
!
  point_coord(1:dim_num,1)  = (/ -b,  z,  z /)
  point_coord(1:dim_num,2)  = (/  z, -a,  z /)
  point_coord(1:dim_num,3)  = (/  z,  a,  z /)
  point_coord(1:dim_num,4)  = (/  z,  z,  b /)
  point_coord(1:dim_num,5)  = (/ -d, -c,  z /)
  point_coord(1:dim_num,6)  = (/ -d,  c,  z /)
  point_coord(1:dim_num,7)  = (/ -d,  z,  d /)
  point_coord(1:dim_num,8)  = (/  z,  z,  z /)
  point_coord(1:dim_num,9)  = (/  z, -c,  d /)
  point_coord(1:dim_num,10) = (/  z,  c,  d /)
!
!  Set the face orders.
!
  face_order(1:face_num) = (/ &
    6, 6, 6, 6 /)
!
!  Set faces.
!
  face_point(1:face_order_max,1:face_num) = reshape ( (/ &
     1,  5,  2,  9,  4,  7, &
     2,  8,  3, 10,  4,  9, &
     3,  6,  1,  7,  4, 10, &
     1,  6,  3,  8,  2,  5 /), (/ face_order_max, face_num /) )

  return
end
subroutine tetrahedron_rhombic_size ( point_num, edge_num, face_num, &
  face_order_max )

!*****************************************************************************80
!
!! tetrahedron_rhombic_size() gives "sizes" for a rhombic tetrahedron in 3D.
!
!  Discussion:
!
!    Call this routine first, in order to learn the required dimensions
!    of arrays to be set up by TETRAHEDRON_RHOMBIC_SHAPE().
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
!  Output:
!
!    integer POINT_NUM, the number of vertices.
!
!    integer EDGE_NUM, the number of edges.
!
!    integer FACE_NUM, the number of faces.
!
!    integer FACE_ORDER_MAX, the maximum order of any face.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer edge_num
  integer face_num
  integer face_order_max
  integer point_num

  point_num = 10
  edge_num = 6
  face_num = 4
  face_order_max = 6

  return
end
subroutine tetrahedron_sample ( t, n, p )

!*****************************************************************************80
!
!! tetrahedron_sample() returns random points in a tetrahedron.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) T(3,4), the tetrahedron vertices.
!
!    integer N, the number of points to sample.
!
!  Output:
!
!    real ( kind = rk ) P(3,N), random points in the tetrahedron.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer n

  real ( kind = rk ) alpha
  real ( kind = rk ) beta
  real ( kind = rk ) gamma
  integer j
  real ( kind = rk ) p(dim_num,n)
  real ( kind = rk ) p12(dim_num)
  real ( kind = rk ) p13(dim_num)
  real ( kind = rk ) r
  real ( kind = rk ) t(dim_num,dim_num+1)
  real ( kind = rk ) tr(dim_num,3)

  do j = 1, n

    call random_number ( harvest = r )
!
!  Interpret R as a percentage of the tetrahedron's volume.
!
!  Imagine a plane, parallel to face 1, so that the volume between
!  vertex 1 and the plane is R percent of the full tetrahedron volume.
!
!  The plane will intersect sides 12, 13, and 14 at a fraction
!  ALPHA = R^1/3 of the distance from vertex 1 to vertices 2, 3, and 4.
!  
    alpha = r**( 1.0D+00 / 3.0D+00 )
!
!  Determine the coordinates of the points on sides 12, 13 and 14 intersected
!  by the plane, which form a triangle TR.
!
    tr(1:dim_num,1) = ( 1.0D+00 - alpha ) * t(1:dim_num,1) &
                                + alpha   * t(1:dim_num,2)
    tr(1:dim_num,2) = ( 1.0D+00 - alpha ) * t(1:dim_num,1) &
                                + alpha   * t(1:dim_num,3)
    tr(1:dim_num,3) = ( 1.0D+00 - alpha ) * t(1:dim_num,1) &
                                + alpha   * t(1:dim_num,4)
!
!  Choose, uniformly at random, a point in this triangle.
!
    call random_number ( harvest = r )
!
!  Interpret R as a percentage of the triangle's area.
!
!  Imagine a line L, parallel to side 1, so that the area between
!  vertex 1 and line L is R percent of the full triangle's area.
!
!  The line L will intersect sides 2 and 3 at a fraction
!  ALPHA = SQRT ( R ) of the distance from vertex 1 to vertices 2 and 3.
!
    beta = sqrt ( r )
!
!  Determine the coordinates of the points on sides 2 and 3 intersected
!  by line L.
!
    p12(1:dim_num) = ( 1.0D+00 - beta ) * tr(1:dim_num,1) &
                               + beta   * tr(1:dim_num,2)
    p13(1:dim_num) = ( 1.0D+00 - beta ) * tr(1:dim_num,1) &
                               + beta   * tr(1:dim_num,3)
!
!  Now choose, uniformly at random, a point on the line L.
!
    call random_number ( harvest = gamma )

    p(1:dim_num,j) = ( 1.0D+00 - gamma ) * p12(1:dim_num) &
                   +             gamma   * p13(1:dim_num)

  end do

  return
end
subroutine tetrahedron_shape ( point_num, face_num, face_order_max, &
  point_coord, face_order, face_point )

!*****************************************************************************80
!
!! tetrahedron_shape() describes a tetrahedron in 3D.
!
!  Discussion:
!
!    Call tetrahedron_size() first, to get dimension information.
!
!    The vertices lie on the unit sphere.
!
!    The dual of the tetrahedron is the tetrahedron.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer POINT_NUM, the number of points.
!
!    integer FACE_NUM, the number of faces.
!
!    integer FACE_ORDER_MAX, the maximum number of 
!    vertices per face.
!
!  Output:
!
!    real ( kind = rk ) POINT_COORD(3,POINT_NUM), the vertices.
!
!    integer FACE_ORDER(FACE_NUM), the number of vertices
!    for each face.
!
!    integer FACE_POINT(FACE_ORDER_MAX,FACE_NUM); 
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

  integer face_order(face_num)
  integer face_point(face_order_max,face_num)
  real ( kind = rk ) point_coord(dim_num,point_num)
!
!  Set the point coordinates.
!
  point_coord(1:dim_num,1:point_num) = reshape ( (/ &
        0.942809D+00,    0.000000D+00,   -0.333333D+00, &
       -0.471405D+00,    0.816497D+00,   -0.333333D+00, &
       -0.471405D+00,   -0.816497D+00,   -0.333333D+00, &
        0.000000D+00,    0.000000D+00,    1.000000D+00 /), &
    (/ dim_num, point_num /) )
!
!  Set the face orders.
!
  face_order(1:face_num) = (/ &
    3, 3, 3, 3 /)
!
!  Set faces.
!
  face_point(1:face_order_max,1:face_num) = reshape ( (/ &
       1, 3, 2, &
       1, 2, 4, &
       1, 4, 3, &
       2, 3, 4 /), (/ face_order_max, face_num /) )

  return
end
subroutine tetrahedron_size ( point_num, edge_num, face_num, &
  face_order_max )

!*****************************************************************************80
!
!! tetrahedron_size() gives "sizes" for a tetrahedron in 3D.
!
!  Discussion:
!
!    Call this routine first, in order to learn the required dimensions
!    of arrays to be set up by tetrahedron_shape().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 July 2003
!
!  Author:
!
!    John Burkardt
!
!  Output:
!
!    integer POINT_NUM, the number of vertices.
!
!    integer EDGE_NUM, the number of edges.
!
!    integer FACE_NUM, the number of faces.
!
!    integer FACE_ORDER_MAX, the maximum order of any face.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer edge_num
  integer face_num
  integer face_order_max
  integer point_num

  point_num = 4
  edge_num = 6
  face_num = 4
  face_order_max = 3

  return
end
subroutine tetrahedron_solid_angles ( tetra, angle )

!*****************************************************************************80
!
!! tetrahedron_solid_angles() computes solid angles of a tetrahedron.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 July 2009
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) TETRA(3,4), the vertices of the tetrahedron.
!
!  Output:
!
!    real ( kind = rk ) ANGLE(4), the solid angles.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) angle(4)
  real ( kind = rk ) dihedral_angles(6)
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) tetra(3,4)

  call tetrahedron_dihedral_angles ( tetra, dihedral_angles )

  angle(1) = dihedral_angles(1) &
           + dihedral_angles(2) &
           + dihedral_angles(3) - r8_pi

  angle(2) = dihedral_angles(1) &
           + dihedral_angles(4) &
           + dihedral_angles(5) - r8_pi

  angle(3) = dihedral_angles(2) &
           + dihedral_angles(4) &
           + dihedral_angles(6) - r8_pi

  angle(4) = dihedral_angles(3) &
           + dihedral_angles(5) &
           + dihedral_angles(6) - r8_pi

  return
end
subroutine tetrahedron_volume ( tetra, volume )

!*****************************************************************************80
!
!! tetrahedron_volume() computes the volume of a tetrahedron.
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
!  Input:
!
!    real ( kind = rk ) TETRA(3,4), the vertices of the tetrahedron.
!
!  Output:
!
!    real ( kind = rk ) VOLUME, the volume of the tetrahedron.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(4,4)
  real ( kind = rk ) r8mat_det_4d
  real ( kind = rk ) tetra(3,4)
  real ( kind = rk ) volume

  a(1:3,1:4) = tetra(1:3,1:4)
  a(4,1:4) = 1.0D+00

  volume = abs ( r8mat_det_4d ( a ) ) / 6.0D+00

  return
end
subroutine triangle_angles_3d ( t, angle )

!*****************************************************************************80
!
!! triangle_angles_3d() computes the angles of a triangle in 3D.
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
!  Input:
!
!    real ( kind = rk ) T(3,3), the triangle vertices.
!
!  Output:
!
!    real ( kind = rk ) ANGLE(3), the angles opposite
!    sides P1-P2, P2-P3 and P3-P1, in radians.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) angle(3)
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) r8_acos
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) t(3,3)
!
!  Compute the length of each side.
!
  a = sqrt ( sum ( ( t(1:3,1) - t(1:3,2) )**2 ) )
  b = sqrt ( sum ( ( t(1:3,2) - t(1:3,3) )**2 ) )
  c = sqrt ( sum ( ( t(1:3,3) - t(1:3,1) )**2 ) )
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
subroutine triangle_area_3d ( t, area )

!*****************************************************************************80
!
!! triangle_area_3d() computes the area of a triangle in 3D.
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
!    Butterworths, 1983.
!
!  Input:
!
!    real ( kind = rk ) T(3,3), the triangle vertices.
!
!  Output:
!
!    real ( kind = rk ) AREA, the area of the triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) area
  real ( kind = rk ) cross(3)
  real ( kind = rk ) t(3,3)
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

