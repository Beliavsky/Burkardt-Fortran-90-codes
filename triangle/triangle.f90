function i4_modp ( i, j )

!*****************************************************************************80
!
!! i4_modp() returns the nonnegative remainder of integer division.
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

  integer, parameter :: rk8 = kind ( 1.0D+00 )

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
function i4_wrap ( ival, ilo, ihi )

!*****************************************************************************80
!
!! i4_wrap() forces an I4 to lie between given limits by wrapping.
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

  integer, parameter :: rk8 = kind ( 1.0D+00 )

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
function line_exp_is_degenerate_nd ( dim_num, p1, p2 )

!*****************************************************************************80
!
!! line_exp_is_degenerate_nd() finds if an explicit line is degenerate in ND.
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
!    Input, real ( kind = rk8 ) P1(DIM_NUM), P2(DIM_NUM), two points on the line.
!
!    Output, logical LINE_EXP_IS_DEGENERATE_ND, is TRUE if the line
!    is degenerate.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer dim_num

  logical line_exp_is_degenerate_nd
  real ( kind = rk8 ) p1(dim_num)
  real ( kind = rk8 ) p2(dim_num)

  line_exp_is_degenerate_nd = ( all ( p1(1:dim_num) == p2(1:dim_num) ) )

  return
end
subroutine line_exp_perp ( p1, p2, p3, p4, flag )

!*****************************************************************************80
!
!! line_exp_perp() computes a line perpendicular to a line and through a point.
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
!    01 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk8 ) P1(2), P2(2), two points on the line.
!
!    Input, real ( kind = rk8 ) P3(2), a point (presumably not on the
!    line (P1,P2)), through which the perpendicular must pass.
!
!    Output, real ( kind = rk8 ) P4(2), a point on the line (P1,P2),
!    such that the line (P3,P4) is perpendicular to the line (P1,P2).
!
!    Output, logical FLAG, is TRUE if the point could not be computed.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) bot
  logical flag
  logical line_exp_is_degenerate_nd
  real ( kind = rk8 ) p1(2)
  real ( kind = rk8 ) p2(2)
  real ( kind = rk8 ) p3(2)
  real ( kind = rk8 ) p4(2)
  real ( kind = rk8 ) r8_huge
  real ( kind = rk8 ) t

  flag = .false.

  if ( line_exp_is_degenerate_nd ( 2, p1, p2 ) ) then
    flag = .true.
    p4(1:2) = r8_huge ( )
    return
  end if

  bot = sum ( ( p2(1:2) - p1(1:2) ) ** 2 )
!
!  (P3-P1) dot (P2-P1) = Norm(P3-P1) * Norm(P2-P1) * Cos(Theta).
!
!  (P3-P1) dot (P2-P1) / Norm(P3-P1)**2 = normalized coordinate T
!  of the projection of (P3-P1) onto (P2-P1).
!
  t = sum ( ( p1(1:2) - p3(1:2) ) &
          * ( p1(1:2) - p2(1:2) ) ) / bot

  p4(1:2) = p1(1:2) + t * ( p2(1:2) - p1(1:2) )

  return
end
subroutine line_exp2imp ( p1, p2, a, b, c )

!*****************************************************************************80
!
!! line_exp2imp() converts an explicit line to implicit form in 2D.
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
!  Parameters:
!
!    Input, real ( kind = rk8 ) P1(2), P2(2), two points on the line.
!
!    Output, real ( kind = rk8 ) A, B, C, the implicit form of the line.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  real ( kind = rk8 ) c
  logical line_exp_is_degenerate_nd
  real ( kind = rk8 ) norm
  real ( kind = rk8 ) p1(2)
  real ( kind = rk8 ) p2(2)
!
!  Take care of degenerate cases.
!
  if ( line_exp_is_degenerate_nd ( 2, p1, p2 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LINE_EXP2IMP - Warning!'
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
function line_imp_is_degenerate ( a, b, c )

!*****************************************************************************80
!
!! line_imp_is_degenerate() finds if an implicit point is degenerate in 2D.
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
!    19 August 2023
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk8 ) A, B, C, the implicit line parameters.
!
!    Output, logical LINE_IMP_IS_DEGENERATE, is true if the
!    line is degenerate.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  real ( kind = rk8 ) c
  logical line_imp_is_degenerate

  call r8_fake_use ( c )

  line_imp_is_degenerate = ( a * a + b * b == 0.0D+00 )

  return
end
subroutine lines_exp_int ( p1, p2, q1, q2, ival, p )

!*****************************************************************************80
!
!! lines_exp_int() determines where two explicit lines intersect in 2D.
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
!    Input, real ( kind = rk8 ) P1(2), P2(2), two points on the first line.
!
!    Input, real ( kind = rk8 ) Q1(2), Q2(2), two points on the second line.
!
!    Output, integer IVAL, reports on the intersection:
!    0, no intersection, the lines may be parallel or degenerate.
!    1, one intersection point, returned in P.
!    2, infinitely many intersections, the lines are identical.
!
!    Output, real ( kind = rk8 ) P(2), if IVAl = 1, P is
!    the intersection point.  Otherwise, P = 0.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a1
  real ( kind = rk8 ) a2
  real ( kind = rk8 ) b1
  real ( kind = rk8 ) b2
  real ( kind = rk8 ) c1
  real ( kind = rk8 ) c2
  integer ival
  logical point_1
  logical point_2
  real ( kind = rk8 ) p(2)
  real ( kind = rk8 ) p1(2)
  real ( kind = rk8 ) p2(2)
  real ( kind = rk8 ) q1(2)
  real ( kind = rk8 ) q2(2)

  ival = 0
  p(1:2) = 0.0D+00
!
!  Check whether either line is a point.
!
  if ( all ( p1(1:2) == p2(1:2) ) ) then
    point_1 = .true.
  else
    point_1 = .false.
  end if

  if ( all ( q1(1:2) == q2(1:2) ) ) then
    point_2 = .true.
  else
    point_2 = .false.
  end if
!
!  Convert the lines to ABC format.
!
  if ( .not. point_1 ) then
    call line_exp2imp ( p1, p2, a1, b1, c1 )
  end if

  if ( .not. point_2 ) then
    call line_exp2imp ( q1, q2, a2, b2, c2 )
  end if
!
!  Search for intersection of the lines.
!
  if ( point_1 .and. point_2 ) then
    if ( all ( p1(1:2) == q1(1:2) ) ) then
      ival = 1
      p(1:2) = p1(1:2)
    end if
  else if ( point_1 ) then
    if ( a2 * p1(1) + b2 * p1(2) == c2 ) then
      ival = 1
      p(1:2) = p1(1:2)
    end if
  else if ( point_2 ) then
    if ( a1 * q1(1) + b1 * q1(2) == c1 ) then
      ival = 1
      p(1:2) = q1(1:2)
    end if
  else
    call lines_imp_int ( a1, b1, c1, a2, b2, c2, ival, p )
  end if

  return
end
subroutine lines_imp_int ( a1, b1, c1, a2, b2, c2, ival, p )

!*****************************************************************************80
!
!! lines_imp_int() determines where two implicit lines intersect in 2D.
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
!    Input, real ( kind = rk8 ) A1, B1, C1, define the first line.
!    At least one of A1 and B1 must be nonzero.
!
!    Input, real ( kind = rk8 ) A2, B2, C2, define the second line.
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
!    Output, real ( kind = rk8 ) P(2), if IVAL = 1, then P is
!    the intersection point.  Otherwise, P = 0.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a(2,3)
  real ( kind = rk8 ) a1
  real ( kind = rk8 ) a2
  real ( kind = rk8 ) b1
  real ( kind = rk8 ) b2
  real ( kind = rk8 ) c1
  real ( kind = rk8 ) c2
  integer info
  integer ival
  logical line_imp_is_degenerate
  real ( kind = rk8 ) p(2)

  p(1:2) = 0.0D+00
!
!  Refuse to handle degenerate lines.
!
  if ( line_imp_is_degenerate ( a1, b1, c1 ) ) then
    ival = -1
    return
  end if

  if ( line_imp_is_degenerate ( a2, b2, c2 ) ) then
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
    p(1:2) = a(1:2,3)
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
!    Input, real ( kind = rk8 ) C, the argument.
!
!    Output, real ( kind = rk8 ) R8_ACOS, an angle whose cosine is C.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) c
  real ( kind = rk8 ) c2
  real ( kind = rk8 ) r8_acos

  c2 = c
  c2 = max ( c2, -1.0D+00 )
  c2 = min ( c2, +1.0D+00 )

  r8_acos = acos ( c2 )

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
!    real ( kind = rk8 ) X, the variable to be "used".
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) x

  if ( x /= x ) then
    write ( *, '(a)' ) '  r8_fake_use(): variable is NAN.'
  end if

  return
end
function r8_huge ( )

!*****************************************************************************80
!
!! r8_huge() returns a very large R8.
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
!    Output, real ( kind = rk8 ) R8_HUGE, a "huge" value.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) r8_huge

  r8_huge = 1.0D+30

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
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer RHS_NUM, the number of right hand sides.
!    RHS_NUM must be at least 0.
!
!    Input/output, real ( kind = rk8 ) A(N,N+RHS_NUM), contains in rows and
!    columns 1 to N the coefficient matrix, and in columns N+1 through
!    N+RHS_NUM, the right hand sides.  On output, the coefficient matrix
!    area has been destroyed, while the right hand sides have
!    been overwritten with the corresponding solutions.
!
!    Output, integer INFO, singularity flag.
!    0, the matrix was not singular, the solutions were computed;
!    J, factorization failed on step J, and the solutions could not
!    be computed.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n
  integer rhs_num

  real ( kind = rk8 ) a(n,n+rhs_num)
  real ( kind = rk8 ) apivot
  real ( kind = rk8 ) factor
  integer i
  integer info
  integer ipivot
  integer j
  real ( kind = rk8 ) t(n+rhs_num)

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
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = rk8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk8 ) a(m,n)
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
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = rk8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ILO, JLO, the first row and column to print.
!
!    Input, integer IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: incx = 5
  integer m
  integer n

  real ( kind = rk8 ) a(m,n)
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
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = rk8 ) X(DIM_NUM), the vector.
!
!    Output, real ( kind = rk8 ) R8VEC_LENGTH, the Euclidean length of the vector.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer dim_num

  real ( kind = rk8 ) r8vec_length
  real ( kind = rk8 ) x(dim_num)

  r8vec_length = sqrt ( sum ( ( x(1:dim_num) ) ** 2 ) )

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
!    Input, real ( kind = rk8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  real ( kind = rk8 ) a(n)
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
subroutine segment_point_dist ( p1, p2, p, dist )

!*****************************************************************************80
!
!! segment_point_dist(): distance ( line segment, point ) in 2D.
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
!    Input, real ( kind = rk8 ) P1(2), P2(2), the endpoints of the line segment.
!
!    Input, real ( kind = rk8 ) P(2), the point whose nearest neighbor on the line
!    segment is to be determined.
!
!    Output, real ( kind = rk8 ) DIST, the distance from the point to the
!    line segment.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) bot
  real ( kind = rk8 ) dist
  real ( kind = rk8 ) p(2)
  real ( kind = rk8 ) p1(2)
  real ( kind = rk8 ) p2(2)
  real ( kind = rk8 ) pn(2)
  real ( kind = rk8 ) t
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
!! segment_point_near(): nearest point on line segment to point in 2D.
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
!    Input, real ( kind = rk8 ) P1(2), P2(2), the endpoints of the line segment.
!
!    Input, real ( kind = rk8 ) P(2), the point whose nearest neighbor
!    on the line segment is to be determined.
!
!    Output, real ( kind = rk8 ) PN(2), the point on the line segment which is
!    nearest the point P.
!
!    Output, real ( kind = rk8 ) DIST, the distance from the point to the 
!    nearest point on the line segment.
!
!    Output, real ( kind = rk8 ) T, the relative position of the point PN
!    to the points P1 and P2.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) bot
  real ( kind = rk8 ) dist
  real ( kind = rk8 ) p(2)
  real ( kind = rk8 ) p1(2)
  real ( kind = rk8 ) p2(2)
  real ( kind = rk8 ) pn(2)
  real ( kind = rk8 ) t
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
subroutine triangle_angles ( t, angle )

!*****************************************************************************80
!
!! triangle_angles() computes the angles of a triangle.
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
!    Input, real ( kind = rk8 ) T(2,3), the triangle vertices.
!
!    Output, real ( kind = rk8 ) ANGLE(3), the angles opposite
!    sides P1-P2, P2-P3 and P3-P1, in radians.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) angle(3)
  real ( kind = rk8 ) b
  real ( kind = rk8 ) c
  real ( kind = rk8 ) r8_acos
  real ( kind = rk8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk8 ) t(2,3)
!
!  Compute the length of each side.
!
  a = sqrt ( sum ( ( t(1:2,1) - t(1:2,2) ) ** 2 ) )
  b = sqrt ( sum ( ( t(1:2,2) - t(1:2,3) ) ** 2 ) )
  c = sqrt ( sum ( ( t(1:2,3) - t(1:2,1) ) ** 2 ) )
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
subroutine triangle_area ( t, area )

!*****************************************************************************80
!
!! triangle_area() computes the area of a triangle in 2D.
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
!    Input, real ( kind = rk8 ) T(2,3), the triangle vertices.
!
!    Output, real ( kind = rk8 ) AREA, the area of the triangle.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) area
  real ( kind = rk8 ) t(2,3)

  area = 0.5D+00 * ( &
      t(1,1) * ( t(2,2) - t(2,3) ) &
    + t(1,2) * ( t(2,3) - t(2,1) ) &
    + t(1,3) * ( t(2,1) - t(2,2) ) )

  return
end
subroutine triangle_centroid ( t, centroid )

!*****************************************************************************80
!
!! triangle_centroid() computes the centroid of a triangle in 2D.
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
!    Butterworths, 1983.
!
!  Parameters:
!
!    Input, real ( kind = rk8 ) T(2,3), the triangle vertices.
!
!    Output, real ( kind = rk8 ) CENTROID(2), the coordinates of the centroid.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) centroid(2)
  integer i
  real ( kind = rk8 ) t(2,3)

  do i = 1, 2
    centroid(i) = sum ( t(i,1:3) ) / 3.0D+00
  end do

  return
end
subroutine triangle_circumcircle ( t, r, pc )

!*****************************************************************************80
!
!! triangle_circumcircle() computes the circumcircle of a triangle in 2D.
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
!    Input, real ( kind = rk8 ) T(2,3), the triangle vertices.
!
!    Output, real ( kind = rk8 ) R, PC(2), the circumradius and circumcenter
!    of the triangle.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  real ( kind = rk8 ) bot
  real ( kind = rk8 ) c
  real ( kind = rk8 ) det
  real ( kind = rk8 ) f(2)
  real ( kind = rk8 ) pc(2)
  real ( kind = rk8 ) r
  real ( kind = rk8 ) top(2)
  real ( kind = rk8 ) t(2,3)
!
!  Circumradius.
!
  a = sqrt ( ( t(1,2) - t(1,1) ) ** 2 + ( t(2,2) - t(2,1) ) ** 2 )
  b = sqrt ( ( t(1,3) - t(1,2) ) ** 2 + ( t(2,3) - t(2,2) ) ** 2 )
  c = sqrt ( ( t(1,1) - t(1,3) ) ** 2 + ( t(2,1) - t(2,3) ) ** 2 )

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
  f(1) = ( t(1,2) - t(1,1) ) ** 2 + ( t(2,2) - t(2,1) ) ** 2
  f(2) = ( t(1,3) - t(1,1) ) ** 2 + ( t(2,3) - t(2,1) ) ** 2

  top(1) =    ( t(2,3) - t(2,1) ) * f(1) - ( t(2,2) - t(2,1) ) * f(2)
  top(2) =  - ( t(1,3) - t(1,1) ) * f(1) + ( t(1,2) - t(1,1) ) * f(2)

  det  =    ( t(2,3) - t(2,1) ) * ( t(1,2) - t(1,1) ) &
          - ( t(2,2) - t(2,1) ) * ( t(1,3) - t(1,1) )

  pc(1:2) = t(1:2,1) + 0.5D+00 * top(1:2) / det

  return
end
subroutine triangle_contains_point ( t, p, inside )

!*****************************************************************************80
!
!! triangle_contains_point() finds if a point is inside a triangle in 2D.
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
!    Input, real ( kind = rk8 ) T(2,3), the triangle vertices.
!    The vertices should be given in counter clockwise order.
!
!    Input, real ( kind = rk8 ) P(2), the point to be checked.
!
!    Output, logical INSIDE, is TRUE if the point is inside the triangle.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  logical inside
  integer j
  integer k
  real ( kind = rk8 ) p(2)
  real ( kind = rk8 ) t(2,3)

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
subroutine triangle_diameter ( t, diameter )

!*****************************************************************************80
!
!! triangle_diameter() computes the diameter of a triangle in 2D.
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
!    Input, real ( kind = rk8 ) T(2,3), the triangle vertices.
!
!    Output, real ( kind = rk8 ) DIAMETER, the diameter of the triangle.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) asq
  real ( kind = rk8 ) b
  real ( kind = rk8 ) bsq
  real ( kind = rk8 ) c
  real ( kind = rk8 ) csq
  real ( kind = rk8 ) diameter
  real ( kind = rk8 ) t(2,3)
  real ( kind = rk8 ) temp
!
!  Compute the squared length of each side.
!
  asq = sum ( t(1:2,1) - t(1:2,2) ) ** 2
  bsq = sum ( t(1:2,2) - t(1:2,3) ) ** 2
  csq = sum ( t(1:2,3) - t(1:2,1) ) ** 2
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
    temp = asq
    asq = bsq
    bsq = temp
  end if

  if ( asq < csq ) then
    temp = asq
    asq = csq
    csq = temp
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
subroutine triangle_edge_length ( t, edge_length )

!*****************************************************************************80
!
!! triangle_edge_length() returns edge lengths of a triangle in 2D.
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
!    Input, real ( kind = rk8 ) T(2,3), the triangle vertices.
!
!    Output, real ( kind = rk8 ) EDGE_LENGTH(3), the length of the edges.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) edge_length(3)
  integer i4_wrap
  integer j1
  integer j2
  real ( kind = rk8 ) r8vec_length
  real ( kind = rk8 ) t(2,3)

  do j1 = 1, 3
    j2 = i4_wrap ( j1 + 1, 1, 3 )
    edge_length(j1) = &
      r8vec_length ( 2, t(1:2,j2) - t(1:2,j1) )
  end do

  return
end
subroutine triangle_incircle ( t, r, pc )

!*****************************************************************************80
!
!! triangle_incircle() computes the inscribed circle of a triangle in 2D.
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
!    Butterworths, 1983.
!
!  Parameters:
!
!    Input, real ( kind = rk8 ) T(2,3), the triangle vertices.
!
!    Output, real ( kind = rk8 ) R, PC(2), the radius and center of the
!    inscribed circle.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  real ( kind = rk8 ) c
  real ( kind = rk8 ) pc(2)
  real ( kind = rk8 ) perimeter
  real ( kind = rk8 ) r
  real ( kind = rk8 ) t(2,3)
!
!  Compute the length of each side.
!
  a = sqrt ( sum ( ( t(1:2,1) - t(1:2,2) ) ** 2 ) )
  b = sqrt ( sum ( ( t(1:2,2) - t(1:2,3) ) ** 2 ) )
  c = sqrt ( sum ( ( t(1:2,3) - t(1:2,1) ) ** 2 ) )

  perimeter = a + b + c

  if ( perimeter == 0.0D+00 ) then
    pc(1:2) = t(1:2,1)
    r = 0.0D+00
    return
  end if

  pc(1:2) = (  &
      b * t(1:2,1) &
    + c * t(1:2,2) &
    + a * t(1:2,3) ) / perimeter

  r = 0.5D+00 * sqrt ( &
      ( - a + b + c )  &
    * ( + a - b + c )  &
    * ( + a + b - c ) / perimeter )

  return
end
function triangle_orientation ( t )

!*****************************************************************************80
!
!! triangle_orientation() determines the orientation of a triangle in 2D.
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
!    Input, real ( kind = rk8 ) T(2,3), the triangle vertices.
!
!    Output, integer TRIANGLE_ORIENTATION, reports if the
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

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) det
  integer triangle_orientation
  real ( kind = rk8 ) t(2,3)

  if ( all ( t(1:2,1) == t(1:2,2) ) .or. &
       all ( t(1:2,2) == t(1:2,3) ) .or. &
       all ( t(1:2,3) == t(1:2,1) ) ) then
    triangle_orientation = 3
    return
  end if

  det = ( t(1,1) - t(1,3) ) * ( t(2,2) - t(2,3) ) &
      - ( t(1,2) - t(1,3) ) * ( t(2,1) - t(2,3) )

  if ( det == 0.0D+00 ) then
    triangle_orientation = 2
  else if ( det < 0.0D+00 ) then
    triangle_orientation = 1
  else if ( 0.0D+00 < det ) then
    triangle_orientation = 0
  end if

  return
end
subroutine triangle_orthocenter ( t, pc, flag )

!*****************************************************************************80
!
!! triangle_orthocenter() computes the orthocenter of a triangle in 2D.
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
!    30 July 2009
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
!  Parameters:
!
!    Input, real ( kind = rk8 ) T(2,3), the triangle vertices.
!
!    Output, real ( kind = rk8 ) PC(2), the orthocenter of the triangle.
!
!    Output, logical FLAG, is TRUE if there was an error condition.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  logical flag
  integer ival
  real ( kind = rk8 ) p23(2)
  real ( kind = rk8 ) p31(2)
  real ( kind = rk8 ) pc(2)
  real ( kind = rk8 ) r8_huge
  real ( kind = rk8 ) t(2,3)
!
!  Determine a point P23 common to the line (P2,P3) and
!  its perpendicular through P1.
!
  call line_exp_perp ( t(1:2,2), t(1:2,3), t(1:2,1), p23, flag )

  if ( flag ) then
    pc(1:2) = r8_huge ( )
    return
  end if
!
!  Determine a point P31 common to the line (P3,P1) and
!  its perpendicular through P2.
!
  call line_exp_perp ( t(1:2,3), t(1:2,1), t(1:2,2), p31, flag )

  if ( flag ) then
    pc(1:2) = r8_huge ( )
    return
  end if
!
!  Determine PC, the intersection of the lines (P1,P23) and (P2,P31).
!
  call lines_exp_int ( t(1:2,1), p23(1:2), t(1:2,2), p31(1:2), ival, pc )

  if ( ival /= 1 ) then
    flag = .true.
    pc(1:2) = r8_huge ( )
    return
  end if

  return
end
subroutine triangle_point_dist ( t, p, dist )

!*****************************************************************************80
!
!! triangle_point_dist(): distance ( triangle, point ) in 2D.
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
!    Input, real ( kind = rk8 ) T(2,3), the triangle vertices.
!
!    Input, real ( kind = rk8 ) P(2), the point to be checked.
!
!    Output, real ( kind = rk8 ) DIST, the distance from the point to the
!    triangle.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) dist
  real ( kind = rk8 ) dist2
  integer i4_wrap
  integer j
  integer jp1
  real ( kind = rk8 ) p(2)
  real ( kind = rk8 ) t(2,3)
!
!  Find the distance to each of the line segments.
!
  dist = huge ( dist )

  do j = 1, 3

    jp1 = i4_wrap ( j + 1, 1, 3 )

    call segment_point_dist ( t(1:2,j), t(1:2,jp1), p, dist2 )

    if ( dist2 < dist ) then
      dist = dist2
    end if

  end do

  return
end
subroutine triangle_point_near ( t, p, pn, dist )

!*****************************************************************************80
!
!! triangle_point_near() computes the nearest point on a triangle in 2D.
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
!    Input, real ( kind = rk8 ) T(2,3), the triangle vertices.
!
!    Input, real ( kind = rk8 ) P(2), the point whose nearest triangle point
!    is to be determined.
!
!    Output, real ( kind = rk8 ) PN(2), the nearest point to P.
!
!    Output, real ( kind = rk8 ) DIST, the distance from the point to the
!    triangle.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) dist
  real ( kind = rk8 ) dist2
  integer i4_wrap
  integer j
  integer jp1
  real ( kind = rk8 ) p(2)
  real ( kind = rk8 ) pn(2)
  real ( kind = rk8 ) pn2(2)
  real ( kind = rk8 ) t(2,3)
  real ( kind = rk8 ) tval
!
!  Find the distance to each of the line segments that make up the edges
!  of the triangle.
!
  dist = huge ( dist )
  pn(1:2) = 0.0D+00

  do j = 1, 3

    jp1 = i4_wrap ( j+1, 1, 3 )

    call segment_point_near ( t(1:2,j), t(1:2,jp1), p, pn2, dist2, tval )

    if ( dist2 < dist ) then
      dist = dist2
      pn(1:2) = pn2(1:2)
    end if

  end do

  return
end
subroutine triangle_quality ( t, quality )

!*****************************************************************************80
!
!! triangle_quality(): "quality" of a triangle in 2D.
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
!    30 July 2009
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
!  Parameters:
!
!    Input, real ( kind = rk8 ) T(2,3), the triangle vertices.
!
!    Output, real ( kind = rk8 ) QUALITY, the quality of the triangle.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  real ( kind = rk8 ) c
  real ( kind = rk8 ) quality
  real ( kind = rk8 ) t(2,3)
!
!  Compute the length of each side.
!
  a = sqrt ( sum ( ( t(1:2,1) - t(1:2,2) ) ** 2 ) )
  b = sqrt ( sum ( ( t(1:2,2) - t(1:2,3) ) ** 2 ) )
  c = sqrt ( sum ( ( t(1:2,3) - t(1:2,1) ) ** 2 ) )

  if ( a * b * c == 0.0D+00 ) then
    quality = 0.0D+00
  else
    quality = ( - a + b + c ) * ( a - b + c ) * ( a + b - c ) &
      / ( a * b * c )
  end if

  return
end
subroutine triangle_reference_sample ( n, p )

!*****************************************************************************80
!
!! triangle_reference_sample() returns random points in the reference triangle.
!
!  Diagram:
!
!       3
!    s  |\
!    i  | \
!    d  |  \
!    e  |   \  side 2
!       |    \
!    3  |     \
!       |      \
!       1-------2
!
!         side 1
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points to generate.
!
!    Output, real ( kind = rk8 ) P(2,N), random points in the triangle.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  real ( kind = rk8 ) alpha
  real ( kind = rk8 ) beta
  integer j
  real ( kind = rk8 ) p(2,n)
  real ( kind = rk8 ) r

  do j = 1, n

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
    alpha = sqrt ( r )
!
!  Now choose, uniformly at random, a point on the line L.
!
    call random_number ( harvest = beta )

    p(1,j) = ( 1.0D+00 - beta ) * alpha
    p(2,j) =             beta   * alpha

  end do

  return
end
subroutine triangle_sample ( t, n, p )

!*****************************************************************************80
!
!! triangle_sample() returns random points in a triangle.
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
!    Input, real ( kind = rk8 ) T(2,3), the triangle vertices.
!
!    Input, integer N, the number of points to generate.
!
!    Output, real ( kind = rk8 ) P(2,N), random points in the triangle.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  real ( kind = rk8 ) alpha(n)
  integer dim
  real ( kind = rk8 ) p(2,n)
  real ( kind = rk8 ) p12(2,n)
  real ( kind = rk8 ) p13(2,n)
  real ( kind = rk8 ) t(2,3)
!
!  For comparison between F90, C++ and MATLAB codes, call R8VEC_UNIFORM_01.
!
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
  do dim = 1, 2

    p12(dim,1:n) = ( 1.0D+00 - alpha(1:n) ) * t(dim,1) &
                             + alpha(1:n)   * t(dim,2)

    p13(dim,1:n) = ( 1.0D+00 - alpha(1:n) ) * t(dim,1) &
                             + alpha(1:n)   * t(dim,3)

  end do
!
!  Now choose, uniformly at random, a point on the line L.
!
  call random_number ( harvest = alpha(1:n) )

  do dim = 1, 2

    p(dim,1:n) = ( 1.0D+00 - alpha(1:n) ) * p12(dim,1:n) &
                           + alpha(1:n)   * p13(dim,1:n)

  end do

  return
end
subroutine triangle_sample_reflect ( t, n, p )

!*****************************************************************************80
!
!! triangle_sample_reflect() returns random points in a triangle using reflection.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 January 2025
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk8 ) t(2,3): the triangle vertices.
!
!    integer n: the number of points to generate.
!
!  Output:
!
!    real ( kind = rk8 ) p(2,n): random points in the triangle.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  integer i
  real ( kind = rk8 ) p(2,n)
  real ( kind = rk8 ) t(2,3)
  real ( kind = rk8 ) u(n)
  real ( kind = rk8 ) v(n)

  call random_number ( harvest = u(1:n) )
  call random_number ( harvest = v(1:n) )

  do i = 1, n
    if ( 1.0D+00 < u(i) + v(i) ) then
      u(i) = 1.0D+00 - u(i)
      v(i) = 1.0D+00 - v(i)
    end if
  end do

  p(1,1:n) = ( 1.0D+00 - u(1:n) - v(1:n) ) * t(1,1) &
                       + u(1:n)            * t(1,2) &
                                + v(1:n)   * t(1,3)

  p(2,1:n) = ( 1.0D+00 - u(1:n) - v(1:n) ) * t(2,1) &
                       + u(1:n)            * t(2,2) &
                                + v(1:n)   * t(2,3)

  return
end
subroutine triangle_xsi_to_xy ( t, xsi, p )

!*****************************************************************************80
!
!! triangle_xsi_to_xy() converts from barycentric to XY coordinates in 2D.
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
!    Input, real ( kind = rk8 ) T(2,3), the triangle vertices.
!
!    Input, real ( kind = rk8 ) XSI(3), the barycentric coordinates of a point.
!    XSI(1) + XSI(2) + XSI(3) should equal 1, but this is not checked.
!
!    Output, real ( kind = rk8 ) P(2), the XY coordinates of the point.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk8 ) p(dim_num)
  real ( kind = rk8 ) t(dim_num,3)
  real ( kind = rk8 ) xsi(dim_num+1)

  p(1:dim_num) = matmul ( t(1:dim_num,1:3), xsi(1:dim_num+1) )

  return
end
subroutine triangle_xy_to_xsi ( t, p, xsi )

!*****************************************************************************80
!
!! triangle_xy_to_xsi() converts from XY to barycentric in 2D.
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
!    Input, real ( kind = rk8 ) T(2,3), the triangle vertices.
!
!    Input, real ( kind = rk8 ) P(2), the XY coordinates of a point.
!
!    Output, real ( kind = rk8 ) XSI(3), the barycentric coordinates of the point.
!    XSI1 + XSI2 + XSI3 should equal 1.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk8 ) det
  real ( kind = rk8 ) p(dim_num)
  real ( kind = rk8 ) t(dim_num,3)
  real ( kind = rk8 ) xsi(3)

  det = ( t(1,1) - t(1,3) ) * ( t(2,2) - t(2,3) ) &
      - ( t(1,2) - t(1,3) ) * ( t(2,1) - t(2,3) )

  xsi(1) = (   ( t(2,2) - t(2,3) ) * ( p(1) - t(1,3) ) &
             - ( t(1,2) - t(1,3) ) * ( p(2) - t(2,3) ) ) / det

  xsi(2) = ( - ( t(2,1) - t(2,3) ) * ( p(1) - t(1,3) ) &
             + ( t(1,1) - t(1,3) ) * ( p(2) - t(2,3) ) ) / det

  xsi(3) = 1.0D+00 - xsi(1) - xsi(2)

  return
end
