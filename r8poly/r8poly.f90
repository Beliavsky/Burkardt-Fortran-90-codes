function r8_sign ( x )

!*****************************************************************************80
!
!! r8_sign() returns the sign of an R8.
!
!  Discussion:
!
!    value = -1 if X < 0;
!    value = +1 if X => 0.
!
!    Note that the standard FORTRAN90 "sign" function is more complicated.
!    In particular,
!
!      Z = sign ( X, Y )
!
!    means that
!
!      Z =   |X| if 0 <= Y;
!          - |X| if Y < 0;
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, the number whose sign is desired.
!
!    Output, real ( kind = rk ) R8_SIGN, the sign of X:
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r8_sign
  real ( kind = rk ) value
  real ( kind = rk ) x

  if ( x < 0.0D+00 ) then
    value = -1.0D+00
  else
    value = +1.0D+00
  end if

  r8_sign = value

  return
end
subroutine r82poly2_print ( a, b, c, d, e, f )

!*****************************************************************************80
!
!! r82poly2_print() prints a second order polynomial in two variables.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, C, D, E, F, the coefficients.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) d
  real ( kind = rk ) e
  real ( kind = rk ) f

  write ( *, &
    '( 2x, f8.4, '' * x^2 + '', f8.4, '' * y^2 + '', f8.4, '' * xy  + '' )' ) &
    a, b, c

  write ( *, &
    '( 2x, f8.4, '' * x + '', f8.4, '' * y + '', f8.4, '' = 0 '' )' ) d, e, f

  return
end
subroutine r82poly2_type ( a, b, c, d, e, f, type )

!*****************************************************************************80
!
!! r82poly2_type() analyzes a second order polynomial in two variables.
!
!  Discussion:
!
!    The polynomial has the form
!
!      A x^2 + B y^2 + C xy + Dx + Ey + F = 0
!
!    The possible types of the solution set are:
!
!     1: a hyperbola;
!        9x^2 -  4y^2       -36x - 24y -  36 = 0
!     2: a parabola;
!        4x^2 +  1y^2 - 4xy + 3x -  4y +   1 = 0;
!     3: an ellipse;
!        9x^2 + 16y^2       +36x - 32y -  92 = 0;
!     4: an imaginary ellipse (no real solutions);
!         x^2 +   y^2       - 6x - 10y + 115 = 0;
!     5: a pair of intersecting lines;
!                        xy + 3x -   y -   3 = 0
!     6: one point;
!         x^2 +  2y^2       - 2x + 16y +  33 = 0;
!     7: a pair of distinct parallel lines;
!                 y^2            -  6y +   8 = 0
!     8: a pair of imaginary parallel lines (no real solutions);
!                 y^2            -  6y +  10 = 0
!     9: a pair of coincident lines.
!                 y^2            -  2y +   1 = 0
!    10: a single line;
!                             2x -   y +   1 = 0;
!    11; all space;
!                                          0 = 0;
!    12; no solutions;
!                                          1 = 0;
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    CRC Press, 30th Edition, 1996, pages 282-284.
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, C, D, E, F, the coefficients.
!
!    Output, integer TYPE, indicates the type of the solution set.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) d
  real ( kind = rk ) delta
  real ( kind = rk ) e
  real ( kind = rk ) f
  real ( kind = rk ) j
  real ( kind = rk ) k
  integer type
!
!  Handle the degenerate case.
!
  if ( a == 0.0D+00 .and. &
       b == 0.0D+00 .and. &
       c == 0.0D+00 ) then
    if ( d == 0.0D+00 .and. e == 0.0D+00 ) then
      if ( f == 0.0D+00 ) then
        type = 11
      else
        type = 12
      end if
    else
      type = 10
    end if
    return
  end if

  delta = &
      8.0D+00 * a * b * f &
    + 2.0D+00 * c * e * d &
    - 2.0D+00 * a * e * e &
    - 2.0D+00 * b * d * d &
    - 2.0D+00 * f * c * c

  j = 4.0D+00 * a * b - c * c

  if ( delta /= 0.0D+00 ) then
    if ( j < 0.0D+00 ) then
      type = 1
    else if ( j == 0.0D+00 ) then
      type = 2
    else if ( 0.0D+00 < j ) then
      if ( sign ( 1.0D+00, delta ) /= sign ( 1.0D+00, ( a + b ) ) ) then
        type = 3
      else if ( sign ( 1.0D+00, delta ) == sign ( 1.0D+00, ( a + b ) ) ) then
        type = 4
      end if
    end if
  else if ( delta == 0.0D+00 ) then
    if ( j < 0.0D+00 ) then
      type = 5
    else if ( 0.0D+00 < j ) then
      type = 6
    else if ( j == 0.0D+00 ) then

      k = 4.0D+00 * ( a + b ) * f - d * d - e * e

      if ( k < 0.0D+00 ) then
        type = 7
      else if ( 0.0D+00 < k ) then
        type = 8
      else if ( k == 0.0D+00 ) then
        type = 9
      end if

    end if
  end if

  return
end
subroutine r82poly2_type_print ( type )

!*****************************************************************************80
!
!! r82poly2_type_print() prints the meaning of the output from R82POLY2_TYPE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer TYPE, the type index returned by R82POLY2_TYPE.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer type

  if ( type == 1 ) then
    write ( *, '(a)' ) '  The set of solutions forms a hyperbola.'
  else if ( type == 2 ) then
    write ( *, '(a)' ) '  The set of solutions forms a parabola.'
  else if ( type == 3 ) then
    write ( *, '(a)' ) '  The set of solutions forms an ellipse.'
  else if ( type == 4 ) then
    write ( *, '(a)' ) '  The set of solutions forms an imaginary ellipse.'
    write ( *, '(a)' ) '  (There are no real solutions).'
  else if ( type == 5 ) then
    write ( *, '(a)' ) &
      '  The set of solutions forms a pair of intersecting lines.'
  else if ( type == 6 ) then
    write ( *, '(a)' ) '  The set of solutions is a single point.'
  else if ( type == 7 ) then
    write ( *, '(a)' ) &
      '  The set of solutions form a pair of distinct parallel lines.'
  else if ( type == 8 ) then
    write ( *, '(a)' ) &
      '  The set of solutions forms a pair of imaginary parallel lines.'
    write ( *, '(a)' ) '  (There are no real solutions).'
  else if ( type == 9 ) then
    write ( *, '(a)' ) &
      '  The set of solutions forms a pair of coincident lines.'
  else if ( type == 10 ) then
    write ( *, '(a)' ) '  The set of solutions forms a single line.'
  else if ( type == 11 ) then
    write ( *, '(a)' ) '  The set of solutions is all space.'
  else if ( type == 12 ) then
    write ( *, '(a)' ) '  The set of solutions is empty.'
  else
    write ( *, '(a)' ) '  This type index is unknown.'
  end if

  return
end
function r8mat_det_3d ( a )

!*****************************************************************************80
!
!! r8mat_det_3d() computes the determinant of a 3 by 3 R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    The formula for the determinant of a 3 by 3 matrix is
!
!        a11 * a22 * a33 - a11 * a23 * a32
!      + a12 * a23 * a31 - a12 * a21 * a33
!      + a13 * a21 * a32 - a13 * a22 * a31
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
!    Input, real ( kind = rk ) A(3,3), the matrix whose determinant is desired.
!
!    Output, real ( kind = rk ) R8MAT_DET_3D, the determinant of the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(3,3)
  real ( kind = rk ) r8mat_det_3d

  r8mat_det_3d = &
         a(1,1) * ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) &
       + a(1,2) * ( a(2,3) * a(3,1) - a(2,1) * a(3,3) ) &
       + a(1,3) * ( a(2,1) * a(3,2) - a(2,2) * a(3,1) )

  return
end
subroutine r8mat_inverse_3d ( a, b, det )

!*****************************************************************************80
!
!! r8mat_inverse_3d() inverts a 3 by 3 R8MAT using Cramer's rule.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    If the determinant is zero, then A is singular, and does not have an
!    inverse.  In that case, B is simply set to zero, and a
!    message is printed.
!
!    If the determinant is nonzero, then its value is roughly an estimate
!    of how nonsingular the matrix A is.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 November 1998
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
  real ( kind = rk ) r8mat_det_3d
!
!  Compute the determinant of A.
!
  det = r8mat_det_3d ( a )
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
  b(1,1) =  ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) / det
  b(1,2) = -( a(1,2) * a(3,3) - a(1,3) * a(3,2) ) / det
  b(1,3) =  ( a(1,2) * a(2,3) - a(1,3) * a(2,2) ) / det

  b(2,1) = -( a(2,1) * a(3,3) - a(2,3) * a(3,1) ) / det
  b(2,2) =  ( a(1,1) * a(3,3) - a(1,3) * a(3,1) ) / det
  b(2,3) = -( a(1,1) * a(2,3) - a(1,3) * a(2,1) ) / det

  b(3,1) =  ( a(2,1) * a(3,2) - a(2,2) * a(3,1) ) / det
  b(3,2) = -( a(1,1) * a(3,2) - a(1,2) * a(3,1) ) / det
  b(3,3) =  ( a(1,1) * a(2,2) - a(1,2) * a(2,1) ) / det

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! r8mat_print() prints an R8MAT.
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
!    12 September 2004
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
!! r8mat_print_some() prints some of an R8MAT.
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
subroutine r8poly2_ex ( x1, y1, x2, y2, x3, y3, x, y, ierror )

!*****************************************************************************80
!
!! r8poly2_ex() finds the extremal point of a parabola determined by three points.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 December 2004
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
!    2, the data lies on a straight line; there is no finite extremal
!    point.
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

  x = 0.5D+00 * ( &
          x1 ** 2 * ( y3 - y2 ) &
        + x2 ** 2 * ( y1 - y3 ) &
        + x3 ** 2 * ( y2 - y1 ) ) / bot

  y = ( &
         ( x - x2 ) * ( x - x3 ) * ( x2 - x3 ) * y1 &
       - ( x - x1 ) * ( x - x3 ) * ( x1 - x3 ) * y2 &
       + ( x - x1 ) * ( x - x2 ) * ( x1 - x2 ) * y3 ) / &
       ( ( x1 - x2 ) * ( x2 - x3 ) * ( x1 - x3 ) )

  return
end
subroutine r8poly2_ex2 ( x1, y1, x2, y2, x3, y3, x, y, a, b, c, ierror )

!*****************************************************************************80
!
!! r8poly2_ex2() finds extremal point of a parabola determined by three points.
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
!    Output, real ( kind = rk ) X, Y, the X coordinate of the extremal
!    point of the parabola, and the value of the parabola at that point.
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

  x = -b / ( 2.0D+00 * a )
  y = a * x * x + b * x + c

  return
end
subroutine r8poly2_root ( a, b, c, r1, r2 )

!*****************************************************************************80
!
!! r8poly2_root() returns the two roots of a quadratic polynomial.
!
!  Discussion:
!
!    The polynomial has the form:
!
!      A * X * X + B * X + C = 0
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, C, the coefficients of the polynomial.
!    A must not be zero.
!
!    Output, complex ( kind = ck8 ) R1, R2, the roots of the polynomial, which
!    might be real and distinct, real and equal, or complex conjugates.
!
  implicit none

  integer, parameter :: ck8 = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  complex ( kind = ck8 ) disc
  complex ( kind = ck8 ) q
  complex ( kind = ck8 ) r1
  complex ( kind = ck8 ) r2

  if ( a == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8POLY2_ROOT - Fatal error!'
    write ( *, '(a)' ) '  The coefficient A is zero.'
    stop 1
  end if

  disc = b * b - 4.0D+00 * a * c
  q = -0.5D+00 * ( b + sign ( 1.0D+00, b ) * sqrt ( disc ) )
  r1 = q / a
  r2 = c / q

  return
end
subroutine r8poly2_rroot ( a, b, c, r1, r2 )

!*****************************************************************************80
!
!! r8poly2_rroot() returns the real parts of the roots of a quadratic polynomial.
!
!  Example:
!
!     A    B    C       roots              R1   R2
!    --   --   --     ------------------   --   --
!     1   -4    3     1          3          1    3
!     1    0    4     2*i      - 2*i        0    0
!     1   -6   10     3 +   i    3 -   i    3    3
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 December 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, C, the coefficients of the quadratic
!    polynomial A * X * X + B * X + C = 0 whose roots are desired.
!    A must not be zero.
!
!    Output, real ( kind = rk ) R1, R2, the real parts of the roots
!    of the polynomial.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) disc
  real ( kind = rk ) q
  real ( kind = rk ) r1
  real ( kind = rk ) r2

  if ( a == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8POLY2_RROOT - Fatal error!'
    write ( *, '(a)' ) '  The coefficient A is zero.'
    stop 1
  end if

  disc = b * b - 4.0D+00 * a * c
  if ( 0.0D+00 <= disc ) then
    q = ( b + sign ( 1.0D+00, b ) * sqrt ( disc ) )
    r1 = -0.5D+00 * q / a
    r2 = -2.0D+00 * c / q
  else
    r1 = b / 2.0D+00 / a
    r2 = b / 2.0D+00 / a
  end if

  return
end
subroutine r8poly2_val ( x1, y1, x2, y2, x3, y3, x, y, yp, ypp )

!*****************************************************************************80
!
!! r8poly2_val() evaluates a parabola defined by three data values.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) X1, Y1, X2, Y2, X3, Y3, three pairs of data.
!    If the X values are distinct, then all the Y values represent
!    actual values of the parabola.
!
!    Three special cases are allowed:
!
!      X1 == X2 /= X3: Y2 is the derivative at X1;
!      X1 /= X2 == X3: Y3 is the derivative at X3;
!      X1 == X2 == X3: Y2 is the derivative at X1, and
!                      Y3 is the second derivative at X1.
!
!    Input, real ( kind = rk ) X, an abscissa at which the parabola is to be
!    evaluated.
!
!    Output, real ( kind = rk ) Y, YP, YPP, the values of the parabola and
!    its first and second derivatives at X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer distinct
  real ( kind = rk ) dif1
  real ( kind = rk ) dif2
  real ( kind = rk ) t
  real ( kind = rk ) x
  real ( kind = rk ) x1
  real ( kind = rk ) x2
  real ( kind = rk ) x3
  real ( kind = rk ) y
  real ( kind = rk ) y1
  real ( kind = rk ) y2
  real ( kind = rk ) y3
  real ( kind = rk ) yp
  real ( kind = rk ) ypp
!
!  If any X's are equal, put them and the Y data first.
!
  if ( x1 == x2 .and. x2 == x3 ) then
    distinct = 1
  else if ( x1 == x2 ) then
    distinct = 2
  else if ( x1 == x3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8POLY2_VAL - Fatal error!'
    write ( *, '(a)' ) '  X1 = X3 =/= X2.'
    write ( *, '(a,g14.6)' ) '  X1 = ', x1
    write ( *, '(a,g14.6)' ) '  X2 = ', x2
    write ( *, '(a,g14.6)' ) '  X3 = ', x3
    stop 1
  else if ( x2 == x3 ) then
    distinct = 2
    t  = x1
    x1 = x2
    x2 = t
    t  = x2
    x2 = x3
    x3 = t
    t  = y1
    y1 = y2
    y2 = t
    t  = y2
    y2 = y3
    y3 = t
  else
    distinct = 3
  end if
!
!  Set up the coefficients.
!
  if ( distinct == 1 ) then

    dif1 = y2
    dif2 = 0.5D+00 * y3

  else if ( distinct == 2 ) then

    dif1 = y2
    dif2 = ( ( y3 - y1 ) / ( x3 - x1 ) - y2 ) / ( x3 - x2 )

  else if ( distinct == 3 ) then

    dif1 = ( y2 - y1 ) / ( x2 - x1 )
    dif2 =  ( ( y3 - y1 ) / ( x3 - x1 ) &
            - ( y2 - y1 ) / ( x2 - x1 ) ) / ( x3 - x2 )

  end if
!
!  Evaluate.
!
  y = y1 + ( x - x1 ) * dif1 + ( x - x1 ) * ( x - x2 ) * dif2
  yp = dif1 + ( 2.0D+00 * x - x1 - x2 ) * dif2
  ypp = 2.0D+00 * dif2

  return
end
subroutine r8poly2_val2 ( dim_num, ndata, tdata, ydata, left, tval, yval )

!*****************************************************************************80
!
!! r8poly2_val2() evaluates a parabolic interpolant through tabular data.
!
!  Discussion:
!
!    This routine is a utility routine used by OVERHAUSER_SPLINE_VAL.
!    It constructs the parabolic interpolant through the data in
!    3 consecutive entries of a table and evaluates this interpolant
!    at a given abscissa value.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the dimension of a single data point.
!    DIM_NUM must be at least 1.
!
!    Input, integer NDATA, the number of data points.
!    NDATA must be at least 3.
!
!    Input, real ( kind = rk ) TDATA(NDATA), the abscissas of the data points.
!    The values in TDATA must be in strictly ascending order.
!
!    Input, real ( kind = rk ) YDATA(DIM_NUM,NDATA), the data points
!    corresponding to the abscissas.
!
!    Input, integer LEFT, the location of the first of the three
!    consecutive data points through which the parabolic interpolant
!    must pass.  1 <= LEFT <= NDATA - 2.
!
!    Input, real ( kind = rk ) TVAL, the value of T at which the parabolic
!    interpolant is to be evaluated.  Normally, TDATA(1) <= TVAL <= T(NDATA),
!    and the data will be interpolated.  For TVAL outside this range,
!    extrapolation will be used.
!
!    Output, real ( kind = rk ) YVAL(DIM_NUM), the value of the parabolic
!    interpolant at TVAL.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ndata
  integer dim_num

  real ( kind = rk ) dif1
  real ( kind = rk ) dif2
  integer i
  integer left
  real ( kind = rk ) t1
  real ( kind = rk ) t2
  real ( kind = rk ) t3
  real ( kind = rk ) tval
  real ( kind = rk ) tdata(ndata)
  real ( kind = rk ) ydata(dim_num,ndata)
  real ( kind = rk ) y1
  real ( kind = rk ) y2
  real ( kind = rk ) y3
  real ( kind = rk ) yval(dim_num)
!
!  Check.
!
  if ( left < 1 .or. ndata-2 < left ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8POLY2_VAL2 - Fatal error!'
    write ( *, '(a)' ) '  LEFT < 1 or NDATA-2 < LEFT.'
    write ( *, '(a,i8)' ) '  LEFT = ', left
    stop 1
  end if

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8POLY2_VAL2 - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1.'
    write ( *, '(a,i8)' ) '  DIM_NUM = ', dim_num
    stop 1
  end if
!
!  Copy out the three abscissas.
!
  t1 = tdata(left)
  t2 = tdata(left+1)
  t3 = tdata(left+2)

  if ( t2 <= t1 .or. t3 <= t2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8POLY2_VAL2 - Fatal error!'
    write ( *, '(a)' ) '  T2 <= T1 or T3 <= T2.'
    write ( *, '(a,g14.6)' ) '  T1 = ', t1
    write ( *, '(a,g14.6)' ) '  T2 = ', t2
    write ( *, '(a,g14.6)' ) '  T3 = ', t3
    stop 1
  end if
!
!  Construct and evaluate a parabolic interpolant for the data
!  in each dimension.
!
  do i = 1, dim_num

    y1 = ydata(i,left)
    y2 = ydata(i,left+1)
    y3 = ydata(i,left+2)

    dif1 = ( y2 - y1 ) / ( t2 - t1 )
    dif2 = ( ( y3 - y1 ) / ( t3 - t1 ) &
           - ( y2 - y1 ) / ( t2 - t1 ) ) / ( t3 - t2 )

    yval(i) = y1 + ( tval - t1 ) * ( dif1 + ( tval - t2 ) * dif2 )

  end do

  return
end
subroutine r8poly3_root ( a, b, c, d, r1, r2, r3 )

!*****************************************************************************80
!
!! r8poly3_root() returns the three roots of a cubic polynomial.
!
!  Discussion:
!
!    The polynomial has the form
!
!      A * X^3 + B * X^2 + C * X + D = 0
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 December 2004
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, C, D, the coefficients of the polynomial.
!    A must not be zero.
!
!    Output, complex ( kind = ck8 ) R1, R2, R3, the roots of the polynomial, which
!    will include at least one real root.
!
  implicit none

  integer, parameter :: ck8 = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) d
  complex ( kind = ck8 ) i
  complex ( kind = ck8 ) one
  real ( kind = rk ) q
  real ( kind = rk ) r
  complex ( kind = ck8 ) r1
  complex ( kind = ck8 ) r2
  complex ( kind = ck8 ) r3
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) s1
  real ( kind = rk ) s2
  real ( kind = rk ) temp
  real ( kind = rk ) theta

  if ( a == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8POLY3_ROOT - Fatal error!'
    write ( *, '(a)' ) '  A must not be zero!'
    stop 1
  end if

  one = cmplx ( 1.0D+00, 0.0D+00, kind = ck8 )
  i = sqrt ( - one )

  q = ( ( b / a ) ** 2 - 3.0D+00 * ( c / a ) ) / 9.0D+00

  r = ( 2.0D+00 * ( b / a ) ** 3 - 9.0D+00 * ( b / a ) * ( c / a ) &
      + 27.0D+00 * ( d / a ) ) / 54.0D+00

  if ( r * r < q * q * q ) then

    theta = acos ( r / sqrt ( q ** 3 ) )
    r1 = - 2.0D+00 * sqrt ( q ) * cos (   theta                     / 3.0D+00 )
    r2 = - 2.0D+00 * sqrt ( q ) * cos ( ( theta + 2.0D+00 * r8_pi ) / 3.0D+00 )
    r3 = - 2.0D+00 * sqrt ( q ) * cos ( ( theta + 4.0D+00 * r8_pi ) / 3.0D+00 )

  else if ( q * q * q <= r * r ) then

    temp = - r + sqrt ( r ** 2 - q ** 3 )
    s1 = sign ( 1.0D+00, temp ) * ( abs ( temp ) ) ** ( 1.0D+00 / 3.0D+00)

    temp = - r - sqrt ( r ** 2 - q ** 3 )
    s2 = sign ( 1.0D+00, temp ) * ( abs ( temp ) ) ** ( 1.0D+00 / 3.0D+00)

    r1 = s1 + s2
    r2 = - 0.5D+00 * ( s1 + s2 ) + i * 0.5D+00 * sqrt ( 3.0D+00 ) * ( s1 - s2 )
    r3 = - 0.5D+00 * ( s1 + s2 ) - i * 0.5D+00 * sqrt ( 3.0D+00 ) * ( s1 - s2 )

  end if

  r1 = r1 - b / ( 3.0D+00 * a )
  r2 = r2 - b / ( 3.0D+00 * a )
  r3 = r3 - b / ( 3.0D+00 * a )

  return
end
subroutine r8poly4_root ( a, b, c, d, e, r1, r2, r3, r4 )

!*****************************************************************************80
!
!! r8poly4_root() returns the four roots of a quartic polynomial.
!
!  Discussion:
!
!    The polynomial has the form:
!
!      A * X^4 + B * X^3 + C * X^2 + D * X + E = 0
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 December 2004
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, C, D, E, the coefficients of the polynomial.
!    A must not be zero.
!
!    Output, complex ( kind = ck8 ) R1, R2, R3, R4, the roots of the polynomial.
!
  implicit none

  integer, parameter :: ck8 = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) a3
  real ( kind = rk ) a4
  real ( kind = rk ) b
  real ( kind = rk ) b3
  real ( kind = rk ) b4
  real ( kind = rk ) c
  real ( kind = rk ) c3
  real ( kind = rk ) c4
  real ( kind = rk ) d
  real ( kind = rk ) d3
  real ( kind = rk ) d4
  real ( kind = rk ) e
  complex ( kind = ck8 ) p
  complex ( kind = ck8 ) q
  complex ( kind = ck8 ) r
  complex ( kind = ck8 ) r1
  complex ( kind = ck8 ) r2
  complex ( kind = ck8 ) r3
  complex ( kind = ck8 ) r4
  complex ( kind = ck8 ) zero

  zero = cmplx ( 0.0D+00, 0.0D+00, kind = ck8 )

  if ( a == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8POLY4_ROOT - Fatal error!'
    write ( *, '(a)') '  A must not be zero!'
    stop 1
  end if

  a4 = b / a
  b4 = c / a
  c4 = d / a
  d4 = e / a
!
!  Set the coefficients of the resolvent cubic equation.
!
  a3 = 1.0D+00
  b3 = -b4
  c3 = a4 * c4 - 4.0D+00 * d4
  d3 = -a4 * a4 * d4 + 4.0D+00 * b4 * d4 - c4 * c4
!
!  Find the roots of the resolvent cubic.
!
  call r8poly3_root ( a3, b3, c3, d3, r1, r2, r3 )
!
!  Choose one root of the cubic, here R1.
!
!  Set R = sqrt ( 0.25D+00 * A4 ** 2 - B4 + R1 )
!
  r = sqrt ( 0.25D+00 * a4 ** 2 - b4 + r1 )

  if ( r /= zero ) then

    p = sqrt ( 0.75D+00 * a4 ** 2 - r ** 2 - 2.0D+00 * b4 &
        + 0.25D+00 * ( 4.0D+00 * a4 * b4 - 8.0D+00 * c4 - a4 ** 3 ) / r )

    q = sqrt ( 0.75D+00 * a4 ** 2 - r ** 2 - 2.0D+00 * b4 &
        - 0.25D+00 * ( 4.0D+00 * a4 * b4 - 8.0D+00 * c4 - a4 ** 3 ) / r )

  else

    p = sqrt ( 0.75D+00 * a4 ** 2 - 2.0D+00 * b4 &
      + 2.0D+00 * sqrt ( r1 ** 2 - 4.0D+00 * d4 ) )

    q = sqrt ( 0.75D+00 * a4 ** 2 - 2.0D+00 * b4 &
      - 2.0D+00 * sqrt ( r1 ** 2 - 4.0D+00 * d4 ) )

  end if
!
!  Set the roots.
!
  r1 = -0.25D+00 * a4 + 0.5D+00 * r + 0.5D+00 * p
  r2 = -0.25D+00 * a4 + 0.5D+00 * r - 0.5D+00 * p
  r3 = -0.25D+00 * a4 - 0.5D+00 * r + 0.5D+00 * q
  r4 = -0.25D+00 * a4 - 0.5D+00 * r - 0.5D+00 * q

  return
end
subroutine r8poly_add ( na, a, nb, b, c )

!*****************************************************************************80
!
!! r8poly_add() adds two R8POLY's.
!
!  Discussion:
!
!    The polynomials are in power sum form.
!
!    The power sum form is:
!
!      p(x) = a(0) + a(1)*x + ... + a(n-1)*x^(n-1) + a(n)*x^(n)
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 November 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NA, the degree of polynomial A.
!
!    Input, real ( kind = rk ) A(0:NA), the coefficients of the first
!    polynomial factor.
!
!    Input, integer NB, the degree of polynomial B.
!
!    Input, real ( kind = rk ) B(0:NB), the coefficients of the
!    second polynomial factor.
!
!    Output, real ( kind = rk )) C(0:max(NA,NB)), the coefficients of A + B.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer na
  integer nb

  real ( kind = rk ) a(0:na)
  real ( kind = rk ) b(0:nb)
  real ( kind = rk ) c(0:max(na,nb))
  real ( kind = rk ) d(0:max(na,nb))

  if ( nb == na ) then
    d(0:na) = a(0:na) + b(0:na)
  else if ( nb < na ) then
    d(0:nb) = a(0:nb) + b(0:nb)
    d(nb+1:na) = a(nb+1:na)
  else if ( na < nb ) then
    d(0:na) = a(0:na) + b(0:na)
    d(na+1:nb) = b(na+1:nb)
  end if

  c(0:max(na,nb)) = d(0:max(na,nb))

  return
end
subroutine r8poly_ant_coef ( n, poly_cof, poly_cof2 )

!*****************************************************************************80
!
!! r8poly_ant_coef() integrates a polynomial in standard form.
!
!  Discussion:
!
!    The antiderivative of a polynomial P(X) is any polynomial Q(X)
!    with the property that d/dX Q(X) = P(X).
!
!    This routine chooses the antiderivative whose constant term is zero.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 October 2019
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the order of the polynomial.
!
!    real ( kind = rk ) POLY_COF(1:N+1), the polynomial coefficients.
!    POLY_COF(1) is the constant term, and POLY_COF(N+1) is the
!    coefficient of X^(N).
!
!  Output:
!
!    real ( kind = rk ) POLY_COF2(1:N+2), the coefficients of
!    the antiderivative polynomial, in standard form.  The constant
!    term is set to zero.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  real ( kind = rk ) poly_cof(n+1)
  real ( kind = rk ) poly_cof2(n+2)
!
!  Set the constant term.
!
  poly_cof2(1) = 0.0D+00
!
!  Integrate the polynomial.
!
  do i = 2, n + 2
    poly_cof2(i) = poly_cof(i-1) / real ( i - 1, kind = rk )
  end do

  return
end
function r8poly_ant_value ( n, poly_cof, xval )

!*****************************************************************************80
!
!! r8poly_ant_value() evaluates the antiderivative of a polynomial.
!
!  Discussion:
!
!    The constant term of the antiderivative is taken to be zero.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 October 2019
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the order of the polynomial.
!
!    real ( kind = rk ) POLY_COF(1:N+1), the polynomial coefficients.
!    POLY_COF(1) is the constant term, and POLY_COF(N+1) is the coefficient of
!    X^(N).
!
!    real ( kind = rk ) XVAL, the evaluation point.
!
!  Output:
!
!    real ( kind = rk ) r8poly_ant_value, the value of the antiderivative at XVAL.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  real ( kind = rk ) poly_cof(n+1)
  real ( kind = rk ) r8poly_ant_value
  real ( kind = rk ) xval
  real ( kind = rk ) yval

  yval = 0.0D+00
  do i = n + 1, 1, -1
    yval = ( yval + poly_cof(i) / real ( i, kind = rk ) ) * xval
  end do

  r8poly_ant_value = yval

  return
end
function r8poly_degree ( na, a )

!*****************************************************************************80
!
!! r8poly_degree() returns the degree of a polynomial.
!
!  Discussion:
!
!    The degree of a polynomial is the index of the highest power
!    of X with a nonzero coefficient.
!
!    The degree of a constant polynomial is 0.  The degree of the
!    zero polynomial is debatable, but this routine returns the
!    degree as 0.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 January 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NA, the dimension of A.
!
!    Input, real ( kind = rk ) A(0:NA), the coefficients of the polynomials.
!
!    Output, integer R8POLY_DEGREE, the degree of A.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer na

  real ( kind = rk ) a(0:na)
  integer r8poly_degree
  integer value

  value = na

  do while ( 0 < value )

    if ( a(value) /= 0.0D+00 ) then
      exit
    end if

    value = value - 1

  end do

  r8poly_degree = value

  return
end
subroutine r8poly_deriv ( n, c, p, cp )

!*****************************************************************************80
!
!! r8poly_deriv() returns the derivative of a polynomial.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the degree of the polynomial.
!
!    Input, real ( kind = rk ) C(0:N), the polynomial coefficients.
!    C(I) is the coefficient of X^I.
!
!    Input, integer P, the order of the derivative.
!    0 means no derivative is taken.
!    1 means first derivative,
!    2 means second derivative and so on.
!    Values of P less than 0 are meaningless.  Values of P greater
!    than N are meaningful, but the code will behave as though the
!    value of P was N+1.
!
!    Output, real ( kind = rk ) CP(0:N-P), the polynomial coefficients of
!    the derivative.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) c(0:n)
  real ( kind = rk ) cp(0:*)
  real ( kind = rk ) cp_temp(0:n)
  integer d
  integer i
  integer p

  if ( n < p ) then
    return
  end if

  cp_temp(0:n) = c(0:n)

  do d = 1, p
    do i = 0, n - d
      cp_temp(i) = real ( i + 1, kind = rk ) * cp_temp(i+1)
    end do
    cp_temp(n-d+1) = 0.0D+00
  end do

  cp(0:n-p) = cp_temp(0:n-p)

  return
end
subroutine r8poly_division ( na, a, nb, b, nq, q, nr, r )

!*****************************************************************************80
!
!! r8poly_division() computes the quotient and remainder of two polynomials.
!
!  Discussion:
!
!    The polynomials are assumed to be stored in power sum form.
!
!    The power sum form of a polynomial is:
!
!      p(x) = a(0) + a(1) * x + ... + a(n-1) * x^(n-1) + a(n) * x^(n)
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 October 2019
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer NA, the dimension of A.
!
!    real ( kind = rk ) A(0:NA), the coefficients of the polynomial
!    to be divided.
!
!    integer NB, the dimension of B.
!
!    real ( kind = rk ) B(0:NB), the coefficients of the divisor
!    polynomial.
!
!  Output:
!
!    integer NQ, the degree of Q.
!    If the divisor polynomial is zero, NQ is returned as -1.
!
!    real ( kind = rk ) Q(0:NA-NB), contains the quotient of A/B.
!    If A and B have full degree, Q should be dimensioned Q(0:NA-NB).
!    In any case, Q(0:NA) should be enough.
!
!    integer NR, the degree of R.
!    If the divisor polynomial is zero, NR is returned as -1.
!
!    real ( kind = rk ) R(0:NB-1), contains the remainder of A/B.
!    If B has full degree, R should be dimensioned R(0:NB-1).
!    Otherwise, R will actually require less space.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer na
  integer nb

  real ( kind = rk ) a(0:na)
  real ( kind = rk ) a2(0:na)
  real ( kind = rk ) b(0:nb)
  integer i
  integer na2
  integer nb2
  integer nq
  integer nr
  real ( kind = rk ) q(0:*)
  real ( kind = rk ) r(0:*)
  integer r8poly_degree

  na2 = r8poly_degree ( na, a )
  nb2 = r8poly_degree ( nb, b )

  if ( b(nb2) == 0.0D+00 ) then
    nq = -1
    nr = -1
    return
  end if

  a2(0:na) = a(0:na)

  nq = na2 - nb2
  nr = nb2 - 1

  do i = nq, 0, -1
    q(i) = a2(i+nb2) / b(nb2)
    a2(i+nb2) = 0.0D+00
    a2(i:i+nb2-1) = a2(i:i+nb2-1) - q(i) * b(0:nb2-1)
  end do

  r(0:nr) = a2(0:nr)

  return
end
subroutine r8poly_lagrange_0 ( npol, xpol, xval, wval )

!*****************************************************************************80
!
!! r8poly_lagrange_0() evaluates the Lagrange factor at a point.
!
!  Formula:
!
!    W(X) = Product ( 1 <= I <= NPOL ) ( X - XPOL(I) )
!
!  Discussion:
!
!    For a set of points XPOL(I), 1 <= I <= NPOL, the IPOL-th Lagrange basis
!    polynomial L(IPOL)(X), has the property:
!
!      L(IPOL)( XPOL(J) ) = delta ( IPOL, J )
!
!    and may be expressed as:
!
!      L(IPOL)(X) = W(X) / ( ( X - XPOL(IPOL) ) * W'(XPOL(IPOL)) )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NPOL, the number of abscissas.
!    NPOL must be at least 1.
!
!    Input, real ( kind = rk ) XPOL(NPOL), the abscissas, which
!    should be distinct.
!
!    Input, real ( kind = rk ) XVAL, the point at which the Lagrange
!    factor is to be evaluated.
!
!    Output, real ( kind = rk ) WVAL, the value of the Lagrange factor at XVAL.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer npol

  real ( kind = rk ) wval
  real ( kind = rk ) xpol(npol)
  real ( kind = rk ) xval

  wval = product ( xval - xpol(1:npol) )

  return
end
subroutine r8poly_lagrange_1 ( npol, xpol, xval, dwdx )

!*****************************************************************************80
!
!! r8poly_lagrange_1() evaluates the first derivative of the Lagrange factor.
!
!  Discussion:
!
!    W(XPOL(1:NPOL))(X) = Product ( 1 <= I <= NPOL ) ( X - XPOL(I) )
!
!    W'(XPOL(1:NPOL))(X)
!      = Sum ( 1 <= J <= NPOL ) Product ( I /= J ) ( X - XPOL(I) )
!
!    We also have the recursion:
!
!      W'(XPOL(1:NPOL))(X) = d/dX ( ( X - XPOL(NPOL) ) * W(XPOL(1:NPOL-1))(X) )
!                    = W(XPOL(1:NPOL-1))(X)
!                    + ( X - XPOL(NPOL) ) * W'(XPOL(1:NPOL-1))(X)
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NPOL, the number of abscissas.
!
!    Input, real ( kind = rk ) XPOL(NPOL), the abscissas, which should
!    be distinct.
!
!    Input, real ( kind = rk ) XVAL, the point at which the Lagrange
!    factor is to be evaluated.
!
!    Output, real ( kind = rk ) DWDX, the derivative of W with respect to X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer npol

  real ( kind = rk ) dwdx
  integer i
  real ( kind = rk ) w
  real ( kind = rk ) xpol(npol)
  real ( kind = rk ) xval

  dwdx = 0.0D+00
  w = 1.0D+00

  do i = 1, npol

    dwdx = w + ( xval - xpol(i) ) * dwdx
    w = w * ( xval - xpol(i) )

  end do

  return
end
subroutine r8poly_lagrange_2 ( npol, xpol, xval, dw2dx2 )

!*****************************************************************************80
!
!! r8poly_lagrange_2() evaluates the second derivative of the Lagrange factor.
!
!  Discussion:
!
!    W(X)  = Product ( 1 <= I <= NPOL ) ( X - XPOL(I) )
!
!    W'(X) = Sum ( 1 <= J <= NPOL )
!            Product ( I /= J ) ( X - XPOL(I) )
!
!    W"(X) = Sum ( 1 <= K <= NPOL )
!            Sum ( J =/ K )
!            Product ( I /= K, J ) ( X - XPOL(I) )
!
!    For a set of points XPOL(I), 1 <= I <= NPOL, the IPOL-th Lagrange basis
!    polynomial L(IPOL)(X), has the property:
!
!      L(IPOL)( XPOL(J) ) = delta ( IPOL, J )
!
!    and may be expressed as:
!
!      L(IPOL)(X) = W(X) / ( ( X - XPOL(IPOL) ) * W'(XPOL(IPOL)) )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NPOL, the number of abscissas.
!    NPOL must be at least 1.
!
!    Input, real ( kind = rk ) XPOL(NPOL), the abscissas, which should
!    be distinct.
!
!    Input, real ( kind = rk ) XVAL, the point at which the Lagrange
!    factor is to be evaluated.
!
!    Output, real ( kind = rk ) DW2DX2, the second derivative of W
!    with respect to XVAL.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer npol

  real ( kind = rk ) dw2dx2
  integer i
  integer j
  integer k
  real ( kind = rk ) term
  real ( kind = rk ) xpol(npol)
  real ( kind = rk ) xval

  dw2dx2 = 0.0D+00

  do k = 1, npol

    do j = 1, npol

      if ( j /= k ) then
        term = 1.0D+00

        do i = 1, npol
          if ( i /= j .and. i /= k ) then
            term = term * ( xval - xpol(i) )
          end if
        end do

        dw2dx2 = dw2dx2 + term

      end if

    end do

  end do

  return
end
subroutine r8poly_lagrange_coef ( npol, ipol, xpol, pcof )

!*****************************************************************************80
!
!! r8poly_lagrange_coef() returns the coefficients of a Lagrange polynomial.
!
!  Discussion:
!
!    Given distinct abscissas XPOL(1:NPOL), the IPOL-th Lagrange
!    polynomial L(IPOL)(X) is defined as the polynomial of degree
!    NPOL - 1 which is 1 at XPOL(IPOL) and 0 at the NPOL - 1 other
!    abscissas.
!
!    A formal representation is:
!
!      L(IPOL)(X) = Product ( 1 <= I <= NPOL, I /= IPOL )
!       ( X - X(I) ) / ( X(IPOL) - X(I) )
!
!    However sometimes it is desirable to be able to write down
!    the standard polynomial coefficients of L(IPOL)(X).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NPOL, the number of abscissas.
!    NPOL must be at least 1.
!
!    Input, integer IPOL, the index of the polynomial to evaluate.
!    IPOL must be between 1 and NPOL.
!
!    Input, real ( kind = rk ) XPOL(NPOL), the abscissas of the
!    Lagrange polynomials.  The entries in XPOL must be distinct.
!
!    Output, real ( kind = rk ) PCOF(0:NPOL-1), the standard polynomial
!    coefficients of the IPOL-th Lagrange polynomial:
!      L(IPOL)(X) = SUM ( 0 <= I <= NPOL-1 ) PCOF(I) * X^I
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer npol

  integer i
  integer indx
  integer ipol
  integer j
  real ( kind = rk ) pcof(0:npol-1)
  logical r8vec_is_distinct
  real ( kind = rk ) xpol(npol)
!
!  Make sure IPOL is legal.
!
  if ( ipol < 1 .or. npol < ipol ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8POLY_LAGRANGE_COEF - Fatal error!'
    write ( *, '(a)' ) '  1 <= IPOL <= NPOL is required.'
    write ( *, '(a,i8)' ) '  IPOL = ', ipol
    write ( *, '(a,i8)' ) '  NPOL = ', npol
    stop 1
  end if
!
!  Check that the abscissas are distinct.
!
  if ( .not. r8vec_is_distinct ( npol, xpol ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8POLY_LAGRANGE_COEF - Fatal error!'
    write ( *, '(a)' ) '  Two or more entries of XPOL are equal:'
    stop 1
  end if

  pcof(0) = 1.0D+00
  pcof(1:npol-1) = 0.0D+00

  indx = 0

  do i = 1, npol

    if ( i /= ipol ) then

      indx = indx + 1

      do j = indx, 0, -1

        pcof(j) = -xpol(i) * pcof(j) / ( xpol(ipol) - xpol(i) )

        if ( 0 < j ) then
          pcof(j) = pcof(j) + pcof(j-1) / ( xpol(ipol) - xpol(i) )
        end if

      end do

    end if

  end do

  return
end
subroutine r8poly_lagrange_factor ( npol, xpol, xval, wval, dwdx )

!*****************************************************************************80
!
!! r8poly_lagrange_factor() evaluates the polynomial Lagrange factor at a point.
!
!  Discussion:
!
!    W(X) = Product ( 1 <= I <= NPOL ) ( X - XPOL(I) )
!
!  Discussion:
!
!    Suppose F(X) is at least N times continuously differentiable in the
!    interval [A,B].  Pick NPOL distinct points XPOL(I) in [A,B] and compute
!    the interpolating polynomial P(X) of order NPOL ( and degree NPOL-1)
!    which passes through all the points ( XPOL(I), F(XPOL(I)) ).
!    Then in the interval [A,B], the maximum error
!
!      abs ( F(X) - P(X) )
!
!    is bounded by:
!
!      C * FNMAX * W(X)
!
!    where
!
!      C is a constant,
!      FNMAX is the maximum value of the NPOL-th derivative of F in [A,B],
!      W(X) is the Lagrange factor.
!
!    Thus, the value of W(X) is useful as part of an estimated bound
!    for the interpolation error.
!
!    Note that the Chebyshev abscissas have the property that they minimize
!    the value of W(X) over the interval [A,B].  Hence, if the abscissas may
!    be chosen arbitrarily, the Chebyshev abscissas have this advantage over
!    other choices.
!
!    For a set of points XPOL(I), 1 <= I <= NPOL, the IPOL-th Lagrange basis
!    polynomial L(IPOL)(X), has the property:
!
!      L(IPOL)( XPOL(J) ) = delta ( IPOL, J )
!
!    and may be expressed as:
!
!      L(IPOL)(X) = W(X) / ( ( X - XPOL(IPOL) ) * W'(XPOL(IPOL)) )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NPOL, the number of abscissas.
!    NPOL must be at least 1.
!
!    Input, real ( kind = rk ) XPOL(NPOL), the abscissas, which should
!    be distinct.
!
!    Input, real ( kind = rk ) XVAL, the point at which the Lagrange
!    factor is to be evaluated.
!
!    Output, real ( kind = rk ) WVAL, the value of the Lagrange factor at XVAL.
!
!    Output, real ( kind = rk ) DWDX, the derivative of W with respect to XVAL.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer npol

  real ( kind = rk ) dwdx
  integer i
  integer j
  real ( kind = rk ) term
  real ( kind = rk ) wval
  real ( kind = rk ) xpol(npol)
  real ( kind = rk ) xval

  wval = product ( xval - xpol(1:npol) )

  dwdx = 0.0D+00

  do i = 1, npol

    term = 1.0D+00

    do j = 1, npol
      if ( i /= j ) then
        term = term * ( xval - xpol(j) )
      end if
    end do

    dwdx = dwdx + term

  end do

  return
end
subroutine r8poly_lagrange_value ( npol, ipol, xpol, xval, pval, dpdx )

!*****************************************************************************80
!
!! r8poly_lagrange_value() evaluates the IPOL-th Lagrange polynomial.
!
!  Discussion:
!
!    Given NPOL distinct abscissas, XPOL(1:NPOL), the IPOL-th Lagrange
!    polynomial L(IPOL)(X) is defined as the polynomial of degree
!    NPOL - 1 which is 1 at XPOL(IPOL) and 0 at the NPOL - 1 other
!    abscissas.
!
!    A formal representation is:
!
!      L(IPOL)(X) = Product ( 1 <= I <= NPOL, I /= IPOL )
!       ( X - X(I) ) / ( X(IPOL) - X(I) )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NPOL, the number of abscissas.
!    NPOL must be at least 1.
!
!    Input, integer IPOL, the index of the polynomial to evaluate.
!    IPOL must be between 1 and NPOL.
!
!    Input, real ( kind = rk ) XPOL(NPOL), the abscissas of the Lagrange
!    polynomials.  The entries in XPOL must be distinct.
!
!    Input, real ( kind = rk ) XVAL, the point at which the IPOL-th
!    Lagrange polynomial is to be evaluated.
!
!    Output, real ( kind = rk ) PVAL, the value of the IPOL-th Lagrange
!    polynomial at XVAL.
!
!    Output, real ( kind = rk ) DPDX, the derivative of the IPOL-th
!    Lagrange polynomial at XVAL.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer npol

  real ( kind = rk ) dpdx
  integer i
  integer ipol
  integer j
  real ( kind = rk ) p2
  real ( kind = rk ) pval
  logical r8vec_is_distinct
  real ( kind = rk ) xpol(npol)
  real ( kind = rk ) xval
!
!  Make sure IPOL is legal.
!
  if ( ipol < 1 .or. npol < ipol ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'r8poly_lagrange_value - Fatal error!'
    write ( *, '(a)' ) '  1 <= IPOL <= NPOL is required.'
    write ( *, '(a,i8)' ) '  IPOL = ', ipol
    stop 1
  end if
!
!  Check that the abscissas are distinct.
!
  if ( .not. r8vec_is_distinct ( npol, xpol ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'r8poly_lagrange_value - Fatal error!'
    write ( *, '(a)' ) '  Two or more entries of XPOL are equal:'
    stop 1
  end if
!
!  Evaluate the polynomial.
!
  pval = 1.0D+00

  do i = 1, npol

    if ( i /= ipol ) then

      pval = pval * ( xval - xpol(i) ) / ( xpol(ipol) - xpol(i) )

    end if

  end do
!
!  Evaluate the derivative, which can be found by summing up the result
!  of differentiating one factor at a time, successively.
!
  dpdx = 0.0D+00

  do i = 1, npol

    if ( i /= ipol ) then

      p2 = 1.0D+00
      do j = 1, npol

        if ( j == i ) then
          p2 = p2                      / ( xpol(ipol) - xpol(j) )
        else if ( j /= ipol ) then
          p2 = p2 * ( xval - xpol(j) ) / ( xpol(ipol) - xpol(j) )
        end if

      end do

      dpdx = dpdx + p2

    end if

  end do

  return
end
subroutine r8poly_multiply ( na, a, nb, b, c )

!*****************************************************************************80
!
!! r8poly_multiply() computes the product of two real polynomials A and B.
!
!  Discussion:
!
!    The polynomials are in power sum form.
!
!    The power sum form is:
!
!      p(x) = a(0) + a(1) * x + ... + a(n-1) * x^(n-1) + a(n) * x^(n)
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NA, the dimension of A.
!
!    Input, real ( kind = rk ) A(0:NA), the coefficients of the first
!    polynomial factor.
!
!    Input, integer NB, the dimension of B.
!
!    Input, real ( kind = rk ) B(0:NB), the coefficients of the second
!    polynomial factor.
!
!    Output, real ( kind = rk ) C(0:NA+NB), the coefficients of A * B.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer na
  integer nb

  real ( kind = rk ) a(0:na)
  real ( kind = rk ) b(0:nb)
  real ( kind = rk ) c(0:na+nb)
  real ( kind = rk ) d(0:na+nb)
  integer i

  d(0:na+nb) = 0.0D+00

  do i = 0, na
    d(i:i+nb) = d(i:i+nb) + a(i) * b(0:nb)
  end do

  c(0:na+nb) = d(0:na+nb)

  return
end
subroutine r8poly_power ( na, a, p, b )

!*****************************************************************************80
!
!! r8poly_power() computes a positive integer power of a polynomial.
!
!  Discussion:
!
!    The power sum form is:
!
!      p(x) = a(0) + a(1)*x + ... + a(n-1)*x^(n-1) + a(n)*x^(n)
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 June 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NA, the dimension of A.
!
!    Input, real ( kind = rk ) A(0:NA), the polynomial to be raised to the power.
!
!    Input, integer P, the nonnegative power to which A is raised.
!
!    Output, real ( kind = rk ) B(0:P*NA), the power of the polynomial.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer na
  integer p

  real ( kind = rk ) a(0:na)
  real ( kind = rk ) b(0:p*na)
  integer i
  integer j
  integer nonzer
!
!  Zero out B.
!
  b(0:p*na) = 0.0D+00
!
!  Search for the first nonzero element in A.
!
  nonzer = -1

  do i = 0, na
    if ( a(i) /= 0.0D+00 ) then
      nonzer = i
      exit
    end if
  end do

  if ( nonzer == -1 ) then
    return
  end if

  b(0) = a(nonzer) ** p

  do i = 1, p * ( na - nonzer )

    if ( i + nonzer <= na ) then
      b(i) = real ( i * p, kind = rk ) * b(0) * a(i+nonzer)
    else
      b(i) = 0.0D+00
    end if

    do j = 1, i - 1

      if ( j + nonzer <= na ) then
        b(i) = b(i) - real ( i - j, kind = rk ) * a(j+nonzer) * b(i-j)
      end if

      if ( i - j + nonzer <= na ) then
        b(i) = b(i) + real ( i - j, kind = rk ) * real ( p, kind = rk ) &
          * b(j) * a(i-j+nonzer)
      end if

    end do

    b(i) = b(i) / ( real ( i, kind = rk ) * a(nonzer) )

  end do
!
!  Shift B up.
!
  do i = p * nonzer, p * na
    b(i) = b(i-p*nonzer)
  end do

  do i = 0, p * nonzer - 1
    b(i) = 0.0D+00
  end do

  return
end
subroutine r8poly_print ( n, a, title )

!*****************************************************************************80
!
!! r8poly_print() prints a polynomial.
!
!  Discussion:
!
!    The power sum form is:
!
!      p(x) = a(0) + a(1) * x + ... + a(n-1) * x^(n-1) + a(n) * x^(n)
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the dimension of A.
!
!    real ( kind = rk ) A(0:N), the polynomial coefficients.
!    A(0) is the constant term and
!    A(N) is the coefficient of X^N.
!
!    character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(0:n)
  integer i
  real ( kind = rk ) mag
  character plus_minus
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  if ( n < 0 ) then
    write ( *, '( ''  p(x) = 0'' )' )
    return
  end if

  if ( a(n) < 0.0D+00 ) then
    plus_minus = '-'
  else
    plus_minus = ' '
  end if

  mag = abs ( a(n) )

  if ( 2 <= n ) then
    write ( *, '( ''  p(x) = '', a1, g14.6, '' * x ^ '', i3 )' ) &
      plus_minus, mag, n
  else if ( n == 1 ) then
    write ( *, '( ''  p(x) = '', a1, g14.6, '' * x'' )' ) &
      plus_minus, mag
  else if ( n == 0 ) then
    write ( *, '( ''  p(x) = '', a1, g14.6 )' ) plus_minus, mag
  end if

  do i = n - 1, 0, -1

    if ( a(i) < 0.0D+00 ) then
      plus_minus = '-'
    else
      plus_minus = '+'
    end if

    mag = abs ( a(i) )

    if ( mag /= 0.0D+00 ) then

      if ( 2 <= i ) then
        write ( *, ' ( ''         '', a1, g14.6, '' * x ^ '', i3 )' ) &
          plus_minus, mag, i
      else if ( i == 1 ) then
        write ( *, ' ( ''         '', a1, g14.6, '' * x'' )' ) plus_minus, mag
      else if ( i == 0 ) then
        write ( *, ' ( ''         '', a1, g14.6 )' ) plus_minus, mag
      end if
    end if

  end do

  return
end
subroutine r8poly_shift ( scale, shift, n, poly_cof )

!*****************************************************************************80
!
!! r8poly_shift() adjusts the coefficients of a polynomial for a new argument.
!
!  Discussion:
!
!    Assuming P(X) is a polynomial in the argument X, of the form:
!
!      P(X) =
!          C(N) * X^N
!        + ...
!        + C(1) * X
!        + C(0),
!
!    and that Z is related to X by the formula:
!
!      Z = SCALE * X + SHIFT
!
!    then this routine computes coefficients C for the polynomial Q(Z):
!
!      Q(Z) =
!          C(N) * Z^N
!        + ...
!        + C(1) * Z
!        + C(0)
!
!    so that:
!
!      Q(Z(X)) = P(X)
!
!  Example:
!
!    P(X) = 2 * X^2 - X + 6
!
!    Z = 2.0 * X + 3.0
!
!    Q(Z) = 0.5 *         Z^2 -  3.5 * Z + 12
!
!    Q(Z(X)) = 0.5 * ( 4.0 * X^2 + 12.0 * X +  9 )
!            - 3.5 * (              2.0 * X +  3 )
!                                           + 12
!
!            = 2.0         * X^2 -  1.0 * X +  6
!
!            = P(X)
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 October 1999
!
!  Reference:
!
!    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
!    Numerical Recipes: The Art of Scientific Computing,
!    Cambridge University Press.
!
!  Parameters:
!
!    Input, real ( kind = rk ) SHIFT, SCALE, the shift and scale applied to X,
!    so that Z = SCALE * X + SHIFT.
!
!    Input, integer N, the number of coefficients.
!
!    Input/output, real ( kind = rk ) POLY_COF(0:N).
!    On input, the coefficient array in terms of the X variable.
!    On output, the coefficient array in terms of the Z variable.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer j
  real ( kind = rk ) poly_cof(0:n)
  real ( kind = rk ) scale
  real ( kind = rk ) shift

  do i = 1, n
    poly_cof(i:n) = poly_cof(i:n) / scale
  end do

  do i = 0, n - 1
    do j = n - 1, i, -1
      poly_cof(j) = poly_cof(j) - shift * poly_cof(j+1)
    end do
  end do

  return
end
function r8poly_value ( m, c, x )

!*****************************************************************************80
!
!! r8poly_value() evaluates a polynomial using a naive method.
!
!  Discussion:
!
!    The polynomial
!
!      p(x) = c0 + c1 * x + c2 * x^2 + ... + cm * x^m
!
!    is to be evaluated at the value X.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 August 2017
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, the degree.
!
!    real ( kind = rk ) C(0:M), the polynomial coefficients.
!    C(I) is the coefficient of X^I.
!
!    real ( kind = rk ) X, the evaluation point.
!
!  Output:
!
!    real ( kind = rk ) R8POLY_VALUE, the polynomial value.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m

  real ( kind = rk ) c(0:m)
  integer i
  real ( kind = rk ) r8poly_value
  real ( kind = rk ) value
  real ( kind = rk ) x
  real ( kind = rk ) xi

  value = c(0)
  xi = 1.0D+00
  do i = 1, m
    xi = xi * x
    value = value + c(i) * xi
  end do

  r8poly_value = value

  return
end
subroutine r8poly_value_2d ( m, c, n, x, y, p )

!*****************************************************************************80
!
!! r8poly_value_2d() evaluates a polynomial in 2 variables, X and Y.
!
!  Discussion:
!
!    We assume the polynomial is of total degree M, and has the form:
!
!      p(x,y) = c00
!             + c10 * x                + c01 * y
!             + c20 * x^2   + c11 * xy + c02 * y^2
!             + ...
!             + cm0 * x^(m) + ...      + c0m * y^m.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the degree of the polynomial.
!
!    Input, real ( kind = rk ) C(T(M+1)), the polynomial coefficients.
!    C(1) is the constant term.  T(M+1) is the M+1-th triangular number.
!    The coefficients are stored consistent with the following ordering
!    of monomials: 1, X, Y, X^2, XY, Y^2, X^3, X^2Y, XY^2, Y^3, X^4, ...
!
!    Input, integer N, the number of evaluation points.
!
!    Input, real ( kind = rk ) X(N), Y(N), the evaluation points.
!
!    Output, real ( kind = rk ) P(N), the value of the polynomial at the
!    evaluation points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) c(*)
  integer ex
  integer ey
  integer j
  integer m
  real ( kind = rk ) p(n)
  integer s
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)

  p(1:n) = 0.0D+00

  j = 0
  do s = 0, m
    do ex = s, 0, -1
      ey = s - ex
      j = j + 1
      p(1:n) = p(1:n) + c(j) * x(1:n) ** ex * y(1:n) ** ey
    end do
  end do

  return
end
function r8poly_value_horner ( m, c, x )

!*****************************************************************************80
!
!! r8poly_value_horner() evaluates a polynomial using Horner's method.
!
!  Discussion:
!
!    The polynomial
!
!      p(x) = c0 + c1 * x + c2 * x^2 + ... + cm * x^m
!
!    is to be evaluated at the value X.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 January 2014
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, the degree.
!
!    real ( kind = rk ) C(0:M), the polynomial coefficients.
!    C(I) is the coefficient of X^I.
!
!    real ( kind = rk ) X, the evaluation point.
!
!  Output:
!
!    real ( kind = rk ) R8POLY_VALUE_HORNER, the polynomial value.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m

  real ( kind = rk ) c(0:m)
  integer i
  real ( kind = rk ) r8poly_value_horner
  real ( kind = rk ) value
  real ( kind = rk ) x

  value = c(m)
  do i = m - 1, 0, -1
    value = value * x + c(i)
  end do

  r8poly_value_horner = value

  return
end
subroutine r8poly_values_horner ( m, c, n, x, p )

!*****************************************************************************80
!
!! r8poly_values_horner() evaluates a polynomial using Horner's method.
!
!  Discussion:
!
!    The polynomial
!
!      p(x) = c0 + c1 * x + c2 * x^2 + ... + cm * x^m
!
!    is to be evaluated at the vector of values X.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, the degree.
!
!    real ( kind = rk ) C(0:M), the polynomial coefficients.
!    C(I) is the coefficient of X^I.
!
!    integer N, the number of evaluation points.
!
!    real ( kind = rk ) X(N), the evaluation points.
!
!  Output:
!
!    real ( kind = rk ) P(N), the polynomial values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) c(0:m)
  integer i
  real ( kind = rk ) p(n)
  real ( kind = rk ) x(n)

  p(1:n) = c(m)
  do i = m - 1, 0, -1
    p(1:n) = p(1:n) * x(1:n) + c(i)
  end do

  return
end
subroutine r8r8_print ( a1, a2, title )

!*****************************************************************************80
!
!! r8r8_print() prints an R8R8.
!
!  Discussion:
!
!    An R8R8 is simply a pair of R8R8's, stored separately.
!
!    A format is used which suggests a coordinate pair:
!
!  Example:
!
!    Center : ( 1.23, 7.45 )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A1, A2, the coordinates of the vector.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a1
  real ( kind = rk ) a2
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '( 2x, a, a4, g14.6, a1, g14.6, a1 )' ) &
      trim ( title ), ' : (', a1, ',', a2, ')'
  else
    write ( *, '( 2x, a1, g14.6, a1, g14.6, a1 )' ) '(', a1, ',', a2, ')'
  end if

  return
end
subroutine r8vec_even ( n, alo, ahi, a )

!*****************************************************************************80
!
!! r8vec_even() returns an R8VEC of evenly spaced values.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    If N is 1, then the midpoint is returned.
!
!    Otherwise, the two endpoints are returned, and N-2 evenly
!    spaced points between them.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of values.
!
!    Input, real ( kind = rk ) ALO, AHI, the low and high values.
!
!    Output, real ( kind = rk ) A(N), N evenly spaced values.
!    Normally, A(1) = ALO and A(N) = AHI.
!    However, if N = 1, then A(1) = 0.5*(ALO+AHI).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) ahi
  real ( kind = rk ) alo
  integer i

  if ( n == 1 ) then

    a(1) = 0.5D+00 * ( alo + ahi )

  else

    do i = 1, n
      a(i) = ( real ( n - i,     kind = rk ) * alo   &
             + real (     i - 1, kind = rk ) * ahi ) &
             / real ( n     - 1, kind = rk )
    end do

  end if

  return
end
subroutine r8vec_even_select ( n, xlo, xhi, ival, xval )

!*****************************************************************************80
!
!! r8vec_even_select() returns the I-th of N evenly spaced values in [ XLO, XHI ].
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    XVAL = ( (N-IVAL) * XLO + (IVAL-1) * XHI ) / real ( N - 1 )
!
!    Unless N = 1, X(1) = XLO and X(N) = XHI.
!
!    If N = 1, then X(1) = 0.5*(XLO+XHI).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of values.
!
!    Input, real ( kind = rk ) XLO, XHI, the low and high values.
!
!    Input, integer IVAL, the index of the desired point.
!    IVAL is normally between 1 and N, but may be any integer value.
!
!    Output, real ( kind = rk ) XVAL, the IVAL-th of N evenly spaced values
!    between XLO and XHI.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer ival
  real ( kind = rk ) xhi
  real ( kind = rk ) xlo
  real ( kind = rk ) xval

  if ( n == 1 ) then

    xval = 0.5D+00 * ( xlo + xhi )

  else

    xval = ( real ( n - ival,     kind = rk ) * xlo   &
           + real (     ival - 1, kind = rk ) * xhi ) &
           / real ( n        - 1, kind = rk )

  end if

  return
end
subroutine r8vec_indicator1 ( n, a )

!*****************************************************************************80
!
!! r8vec_indicator1() sets an R8VEC to the indicator vector (1,2,3,...).
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
!    27 September 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements of A.
!
!    Output, real ( kind = rk ) A(N), the array.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  integer i

  do i = 1, n
    a(i) = real ( i, kind = rk )
  end do

  return
end
function r8vec_is_distinct ( n, a )

!*****************************************************************************80
!
!! r8vec_is_distinct() is true if the entries in an R8VEC are distinct.
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
!    31 March 2018
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real ( kind = rk ) A(N), the vector to be checked.
!
!    Output, logical R8VEC_IS_DISTINCT is TRUE if the entries
!    are distinct.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  integer i
  integer j
  logical r8vec_is_distinct
  logical value

  value = .true.

  do i = 2, n
    do j = 1, i - 1
      if ( a(i) == a(j) ) then
        value = .false.
        exit
      end if
    end do
  end do

  r8vec_is_distinct = value

  return
end
subroutine r8vec_linspace ( n, a, b, x )

!*****************************************************************************80
!
!! r8vec_linspace() creates a vector of linearly spaced values.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    4 points evenly spaced between 0 and 12 will yield 0, 4, 8, 12.
!
!    In other words, the interval is divided into N-1 even subintervals,
!    and the endpoints of intervals are used as the points.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real ( kind = rk ) A, B, the first and last entries.
!
!    Output, real ( kind = rk ) X(N), a vector of linearly spaced data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a
  real ( kind = rk ) b
  integer i
  real ( kind = rk ) x(n)

  if ( n == 1 ) then

    x(1) = ( a + b ) / 2.0D+00

  else

    do i = 1, n
      x(i) = ( real ( n - i,     kind = rk ) * a   &
             + real (     i - 1, kind = rk ) * b ) &
             / real ( n     - 1, kind = rk )
    end do

  end if

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
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
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
!    My vector:  1.0    2.1    3.2    4.3    5.4
!                6.5    7.6    8.7    9.8   10.9
!               11.0
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 August 2016
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
  integer ihi
  integer ilo
  integer j
  character ( len = * ) title
  integer title_length

  title_length = len_trim ( title )

  do ilo = 1, n, 5
    if ( ilo == 1 ) then
      write ( *, '(a)', advance = 'NO' ) trim ( title )
    else
      do i = 1, title_length
        write ( *, '(1x)', advance = 'NO' )
      end do
    end if
    write ( *, '(2x)', advance = 'NO' )
    ihi = min ( ilo + 5 - 1, n )
    do j = ilo, ihi
      write ( *, '(g14.6)', advance = 'NO' ) a(j)
    end do
    write ( *, '(a)' )

  end do

  return
end
subroutine r8vec2_print ( n, a1, a2, title )

!*****************************************************************************80
!
!! r8vec2_print() prints an R8VEC2.
!
!  Discussion:
!
!    An R8VEC2 is a dataset consisting of N pairs of R8's, stored
!    as two separate vectors A1 and A2.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, real ( kind = rk ) A1(N), A2(N), the vectors to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a1(n)
  real ( kind = rk ) a2(n)
  integer i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i4,2x,g14.6,2x,g14.6)' ) i, a1(i), a2(i)
  end do

  return
end
subroutine roots_to_r8poly ( n, x, c )

!*****************************************************************************80
!
!! roots_to_r8poly() converts polynomial roots to polynomial coefficients.
!
!  Discussion:
!
!    p(x) = c(0) x^n + x(1) x^(n-1) + ... + c(n)
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 December 2023
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of roots specified.
!
!    real ( kind = rk ) X(N), the roots.
!
!  Output:
!
!    real ( kind = rk ) C(0:N), the coefficients of the polynomial.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) c(0:n)
  integer i
  integer j
  real ( kind = rk ) x(n)
!
!  Initialize C to (0, 0, ..., 0, 1).
!  Essentially, we are setting up a divided difference table.
!
  c(0:n-1) = 0.0D+00
  c(n) = 1.0D+00
!
!  Convert to standard polynomial form by shifting the abscissas
!  of the divided difference table to 0.
!
  do j = 1, n
    do i = 1, n + 1 - j
      c(n-i) = c(n-i) - x(n+1-i-j+1) * c(n-i+1)
    end do
  end do

  return
end
subroutine r8mat_vand2 ( n, x, a )

!*****************************************************************************80
!
!! r8mat_vand2() returns the N by N row Vandermonde matrix A.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    The row Vandermonde matrix returned by this routine reads "across"
!    rather than down.  In particular, each row begins with a 1, followed by
!    some value X, followed by successive powers of X.
!
!  Formula:
!
!    A(I,J) = X(I)^(J-1)
!
!  Properties:
!
!    A is nonsingular if, and only if, the X values are distinct.
!
!    The determinant of A is
!
!      det(A) = product ( 2 <= I <= N ) (
!        product ( 1 <= J <= I-1 ) ( ( X(I) - X(J) ) ) ).
!
!    The matrix A is generally ill-conditioned.
!
!  Example:
!
!    N = 5, X = (2, 3, 4, 5, 6)
!
!    1 2  4   8   16
!    1 3  9  27   81
!    1 4 16  64  256
!    1 5 25 125  625
!    1 6 36 216 1296
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
!    Input, integer N, the order of the matrix desired.
!
!    Input, real ( kind = rk ) X(N), the values that define A.
!
!    Output, real ( kind = rk ) A(N,N), the N by N row Vandermonde matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,n)
  integer i
  integer j
  real ( kind = rk ) x(n)

  do i = 1, n
    do j = 1, n

      if ( j == 1 .and. x(i) == 0.0D+00 ) then
        a(i,j) = 1.0D+00
      else
        a(i,j) = x(i) ** ( j - 1 )
      end if

    end do
  end do

  return
end

