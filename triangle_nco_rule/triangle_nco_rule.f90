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

  integer i
  integer i4_modp
  integer j
  integer value

  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MODP - Fatal error!'
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
!! I4_WRAP forces an I4 to lie between given limits by wrapping.
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
!    Input, integer ILO, IHI, the desired bounds for the
!    integer value.
!
!    Output, integer I4_WRAP, a "wrapped" version of IVAL.
!
  implicit none

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
subroutine reference_to_physical_t3 ( node_xy, n, ref, phy )

!*****************************************************************************80
!
!! REFERENCE_TO_PHYSICAL_T3 maps T3 reference points to physical points.
!
!  Discussion:
!
!    Given the vertices of an order 3 physical triangle and a point
!    (XSI,ETA) in the reference triangle, the routine computes the value
!    of the corresponding image point (X,Y) in physical space.
!
!    This routine is also appropriate for an order 4 triangle,
!    as long as the fourth node is the centroid of the triangle.
!
!    This routine may also be appropriate for an order 6
!    triangle, if the mapping between reference and physical space
!    is linear.  This implies, in particular, that the sides of the
!    image triangle are straight and that the "midside" nodes in the
!    physical triangle are literally halfway along the sides of
!    the physical triangle.
!
!  Reference Element T3:
!
!    |
!    1  3
!    |  |\
!    |  | \
!    S  |  \
!    |  |   \
!    |  |    \
!    0  1-----2
!    |
!    +--0--R--1-->
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 June 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) NODE_XY(2,3), the coordinates of the vertices.
!    The vertices are assumed to be the images of (0,0), (1,0) and
!    (0,1) respectively.
!
!    Input, integer N, the number of objects to transform.
!
!    Input, real ( kind = rk ) REF(2,N), points in the reference triangle.
!
!    Output, real ( kind = rk ) PHY(2,N), corresponding points in the
!    physical triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  real ( kind = rk ) node_xy(2,3)
  real ( kind = rk ) phy(2,n)
  real ( kind = rk ) ref(2,n)

  do i = 1, 2
    phy(i,1:n) = node_xy(i,1) * ( 1.0D+00 - ref(1,1:n) - ref(2,1:n) ) &
               + node_xy(i,2) *             ref(1,1:n)                &
               + node_xy(i,3) *                          ref(2,1:n)
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
subroutine triangle_area ( node_xy, area )

!*****************************************************************************80
!
!! TRIANGLE_AREA computes the area of a triangle.
!
!  Discussion:
!
!    If the triangle's vertices are given in counterclockwise order,
!    the area will be positive.  If the triangle's vertices are given
!    in clockwise order, the area will be negative!
!
!    If you cannot guarantee counterclockwise order, and you need to
!    have the area positive, then you can simply take the absolute value
!    of the result of this routine.
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
!    Input, real ( kind = rk ) NODE_XY(2,3), the triangle vertices.
!
!    Output, real ( kind = rk ) AREA, the area of the triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) area
  real ( kind = rk ) node_xy(2,3)

  area = 0.5D+00 * ( &
      node_xy(1,1) * ( node_xy(2,2) - node_xy(2,3) ) &
    + node_xy(1,2) * ( node_xy(2,3) - node_xy(2,1) ) &
    + node_xy(1,3) * ( node_xy(2,1) - node_xy(2,2) ) )

  return
end
subroutine triangle_nco_degree ( rule, degree )

!*****************************************************************************80
!
!! TRIANGLE_NCO_DEGREE returns the degree of an NCO rule for the triangle.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Peter Silvester,
!    Symmetric Quadrature Formulae for Simplexes,
!    Mathematics of Computation,
!    Volume 24, Number 109, January 1970, pages 95-100.
!
!  Parameters:
!
!    Input, integer RULE, the index of the rule.
!
!    Output, integer DEGREE, the polynomial degree of exactness of
!    the rule.
!
  implicit none

  integer degree
  integer rule

  if ( 1 <= rule .and. rule <= 9 ) then

    degree = rule - 1

  else

    degree = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGLE_NCO_DEGREE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal RULE = ', rule
    stop 1

  end if

  return
end
subroutine triangle_nco_order_num ( rule, order_num )

!*****************************************************************************80
!
!! TRIANGLE_NCO_ORDER_NUM returns the order of an NCO rule for the triangle.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Peter Silvester,
!    Symmetric Quadrature Formulae for Simplexes,
!    Mathematics of Computation,
!    Volume 24, Number 109, January 1970, pages 95-100.
!
!  Parameters:
!
!    Input, integer RULE, the index of the rule.
!
!    Output, integer ORDER_NUM, the order (number of points)
!    of the rule.
!
  implicit none

  integer order_num
  integer rule
  integer, allocatable, dimension ( : ) :: suborder
  integer suborder_num

  call triangle_nco_suborder_num ( rule, suborder_num )

  allocate ( suborder(1:suborder_num) )

  call triangle_nco_suborder ( rule, suborder_num, suborder )

  order_num = sum ( suborder(1:suborder_num) )

  deallocate ( suborder )

  return
end
subroutine triangle_nco_rule ( rule, order_num, xy, w )

!*****************************************************************************80
!
!! TRIANGLE_NCO_RULE returns the points and weights of an NCO rule.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Peter Silvester,
!    Symmetric Quadrature Formulae for Simplexes,
!    Mathematics of Computation,
!    Volume 24, Number 109, January 1970, pages 95-100.
!
!  Parameters:
!
!    Input, integer RULE, the index of the rule.
!
!    Input, integer ORDER_NUM, the order (number of points)
!    of the rule.
!
!    Output, real ( kind = rk ) XY(2,ORDER_NUM), the points of the rule.
!
!    Output, real ( kind = rk ) W(ORDER_NUM), the weights of the rule.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer order_num

  integer i4_wrap
  integer k
  integer o
  integer rule
  integer s
  integer, allocatable, dimension ( : ) :: suborder
  integer suborder_num
  real ( kind = rk ), allocatable, dimension ( : ) :: suborder_w
  real ( kind = rk ), allocatable, dimension ( :, : ) :: suborder_xyz
  real ( kind = rk ) w(order_num)
  real ( kind = rk ) xy(2,order_num)
!
!  Get the suborder information.
!
  call triangle_nco_suborder_num ( rule, suborder_num )

  allocate ( suborder(suborder_num) )
  allocate ( suborder_xyz(3,suborder_num) )
  allocate ( suborder_w(suborder_num) )

  call triangle_nco_suborder ( rule, suborder_num, suborder )

  call triangle_nco_subrule ( rule, suborder_num, suborder_xyz, suborder_w )
!
!  Expand the suborder information to a full order rule.
!
  o = 0

  do s = 1, suborder_num

    if ( suborder(s) == 1 ) then

      o = o + 1
      xy(1:2,o) = suborder_xyz(1:2,s)
      w(o) = suborder_w(s)

    else if ( suborder(s) == 3 ) then

      do k = 1, 3
        o = o + 1
        xy(1,o) = suborder_xyz ( i4_wrap(k,   1, 3 ), s )
        xy(2,o) = suborder_xyz ( i4_wrap(k+1, 1, 3 ), s )
        w(o) = suborder_w(s)
      end do

    else if ( suborder(s) == 6 ) then

      do k = 1, 3
        o = o + 1
        xy(1,o) = suborder_xyz ( i4_wrap(k,  1,3), s )
        xy(2,o) = suborder_xyz ( i4_wrap(k+1,1,3), s )
        w(o) = suborder_w(s)
      end do

      do k = 1, 3
        o = o + 1
        xy(1,o) = suborder_xyz ( i4_wrap(k+1,1,3), s )
        xy(2,o) = suborder_xyz ( i4_wrap(k,  1,3), s )
        w(o) = suborder_w(s)
      end do

    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TRIANGLE_NCO_RULE - Fatal error!'
      write ( *, '(a,i8,a,i8)' ) '  Illegal SUBORDER(', s, ') = ', suborder(s)
      write ( *, '(a,i8)' ) '  RULE =    ', rule
      write ( *, '(a,i8)' ) '  ORDER_NUM = ', order_num
      stop 1

    end if

  end do

  deallocate ( suborder )
  deallocate ( suborder_xyz )
  deallocate ( suborder_w )

  return
end
subroutine triangle_nco_rule_num ( rule_num )

!*****************************************************************************80
!
!! TRIANGLE_NCO_RULE_NUM returns the number of NCO rules available.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Peter Silvester,
!    Symmetric Quadrature Formulae for Simplexes,
!    Mathematics of Computation,
!    Volume 24, Number 109, January 1970, pages 95-100.
!
!  Parameters:
!
!    Output, integer RULE_NUM, the number of rules available.
!
  implicit none

  integer rule_num

  rule_num = 9

  return
end
subroutine triangle_nco_suborder ( rule, suborder_num, suborder )

!*****************************************************************************80
!
!! TRIANGLE_NCO_SUBORDER returns the suborders for an NCO rule.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Peter Silvester,
!    Symmetric Quadrature Formulae for Simplexes,
!    Mathematics of Computation,
!    Volume 24, Number 109, January 1970, pages 95-100.
!
!  Parameters:
!
!    Input, integer RULE, the index of the rule.
!
!    Input, integer SUBORDER_NUM, the number of suborders of
!    the rule.
!
!    Output, integer SUBORDER(SUBORDER_NUM), the suborders of
!    the rule.
!
  implicit none

  integer suborder_num

  integer rule
  integer suborder(suborder_num)

  if ( rule == 1 ) then
    suborder(1:suborder_num) = (/ &
      1 /)
  else if ( rule == 2 ) then
    suborder(1:suborder_num) = (/ &
      3 /)
  else if ( rule == 3 ) then
    suborder(1:suborder_num) = (/ &
      3, 3 /)
  else if ( rule == 4 ) then
    suborder(1:suborder_num) = (/ &
      3, 6, 1 /)
  else if ( rule == 5 ) then
    suborder(1:suborder_num) = (/ &
      3, 6, 3, 3 /)
  else if ( rule == 6 ) then
    suborder(1:suborder_num) = (/ &
      3, 6, 6, 3, 3 /)
  else if ( rule == 7 ) then
    suborder(1:suborder_num) = (/ &
      3, 6, 6, 3, 3, 6, 1 /)
  else if ( rule == 8 ) then
    suborder(1:suborder_num) = (/ &
      3, 6, 6, 3, 6, 6, 3, 3 /)
  else if ( rule == 9 ) then
    suborder(1:suborder_num) = (/ &
      3, 6, 6, 3, 6, 6, 3, 6, 3, 3 /)
  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGLE_NCO_SUBORDER - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal RULE = ', rule
    stop 1

  end if

  return
end
subroutine triangle_nco_suborder_num ( rule, suborder_num )

!*****************************************************************************80
!
!! TRIANGLE_NCO_SUBORDER_NUM returns the number of suborders for an NCO rule.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Peter Silvester,
!    Symmetric Quadrature Formulae for Simplexes,
!    Mathematics of Computation,
!    Volume 24, Number 109, January 1970, pages 95-100.
!
!  Parameters:
!
!    Input, integer RULE, the index of the rule.
!
!    Output, integer SUBORDER_NUM, the number of suborders of
!    the rule.
!
  implicit none

  integer rule
  integer, dimension(1:9) :: suborder = (/ &
     1, 1, 2, 3, 4, 5, 7, 8, 10 /)

  integer suborder_num

  if ( 1 <= rule .and. rule <= 9 ) then
    suborder_num = suborder(rule)
  else
    suborder_num = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGLE_NCO_SUBORDER_NUM - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal RULE = ', rule
    stop 1
  end if

  return
end
subroutine triangle_nco_subrule ( rule, suborder_num, suborder_xyz, suborder_w )

!*****************************************************************************80
!
!! TRIANGLE_NCO_SUBRULE returns a compressed NCO rule.
!
!  Discussion:
!
!    The value "3757007" listed in the reference as the sixth numerator
!    for the open case with degree 8 (our rule 9) is incorrect.  It
!    should be "-3757007".
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Peter Silvester,
!    Symmetric Quadrature Formulae for Simplexes,
!    Mathematics of Computation,
!    Volume 24, Number 109, January 1970, pages 95-100.
!
!  Parameters:
!
!    Input, integer RULE, the index of the rule.
!
!    Input, integer SUBORDER_NUM, the number of suborders of
!    the rule.
!
!    Output, real ( kind = rk ) SUBORDER_XYZ(3,SUBORDER_NUM),
!    the barycentric coordinates of the abscissas.
!
!    Output, real ( kind = rk ) SUBORDER_W(SUBORDER_NUM), the
!    suborder weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer suborder_num

  integer rule
  real ( kind = rk ) suborder_w(suborder_num)
  integer suborder_w_n(suborder_num)
  integer suborder_w_d
  real ( kind = rk ) suborder_xyz(3,suborder_num)
  integer suborder_xyz_n(3,suborder_num)
  integer suborder_xyz_d

  if ( rule == 1 ) then

    suborder_xyz_n(1:3,1:suborder_num) = reshape ( (/ &
      0,  0, 0  &
      /), (/ 3, suborder_num /) )

    suborder_xyz_d = 0

    suborder_w_n(1:suborder_num) = (/ &
      1 /)

    suborder_w_d = 1

  else if ( rule == 2 ) then

    suborder_xyz_n(1:3,1:suborder_num) = reshape ( (/ &
      1, 0, 0  &
      /), (/ 3, suborder_num /) )

    suborder_xyz_d = 1

    suborder_w_n(1:suborder_num) = (/ &
      1 /)

    suborder_w_d = 3

  else if ( rule == 3 ) then

    suborder_xyz_n(1:3,1:suborder_num) = reshape ( (/ &
      2, 0, 0, &
      1, 1, 0  &
      /), (/ 3, suborder_num /) )

    suborder_xyz_d = 2

    suborder_w_n(1:suborder_num) = (/ &
      7, -3 /)

     suborder_w_d = 12

  else if ( rule == 4 ) then

    suborder_xyz_n(1:3,1:suborder_num) = reshape ( (/ &
      3, 0, 0,  &
      2, 1, 0,  &
      1, 1, 1   &
      /), (/ 3, suborder_num /) )

    suborder_xyz_d = 3

    suborder_w_n(1:suborder_num) = (/ &
      8, 3, -12 /)

     suborder_w_d = 30

  else if ( rule == 5 ) then

    suborder_xyz_n(1:3,1:suborder_num) = reshape ( (/ &
      4, 0, 0, &
      3, 1, 0,  &
      2, 2, 0,  &
      2, 1, 1   &
      /), (/ 3, suborder_num /) )

    suborder_xyz_d = 4

    suborder_w_n(1:suborder_num) = (/ &
      307, -316, 629, -64 /)

     suborder_w_d = 720

  else if ( rule == 6 ) then

    suborder_xyz_n(1:3,1:suborder_num) = reshape ( (/ &
      5, 0, 0,  &
      4, 1, 0,  &
      3, 2, 0,  &
      3, 1, 1,  &
      2, 2, 1  &
      /), (/ 3, suborder_num /) )

    suborder_xyz_d = 5

    suborder_w_n(1:suborder_num) = (/ &
      71, -13, 57, -167, 113 /)

     suborder_w_d = 315

  else if ( rule == 7 ) then

    suborder_xyz_n(1:3,1:suborder_num) = reshape ( (/ &
      6, 0, 0,  &
      5, 1, 0,  &
      4, 2, 0,  &
      4, 1, 1,  &
      3, 3, 0,  &
      3, 2, 1,  &
      2, 2, 2   &
      /), (/ 3, suborder_num /) )

    suborder_xyz_d = 6

    suborder_w_n(1:suborder_num) = (/ &
      767, -1257, 2901, 387, -3035, -915, 3509 /)

     suborder_w_d = 2240

  else if ( rule == 8 ) then

    suborder_xyz_n(1:3,1:suborder_num) = reshape ( (/ &
      7, 0, 0,  &
      6, 1, 0,  &
      5, 2, 0,  &
      5, 1, 1,  &
      4, 3, 0,  &
      4, 2, 1,  &
      3, 3, 1,  &
      3, 2, 2   &
      /), (/ 3, suborder_num /) )

    suborder_xyz_d = 7

    suborder_w_n(1:suborder_num) = (/ &
      898, -662, 1573, -2522, -191, 2989, -5726, 1444 /)

     suborder_w_d = 4536

  else if ( rule == 9 ) then

    suborder_xyz_n(1:3,1:suborder_num) = reshape ( (/ &
      8, 0, 0,  &
      7, 1, 0,  &
      6, 2, 0,  &
      6, 1, 1,  &
      5, 3, 0,  &
      5, 2, 1,  &
      4, 4, 0,  &
      4, 3, 1,  &
      4, 2, 2,  &
      3, 3, 2   &
      /), (/ 3, suborder_num /) )

    suborder_xyz_d = 8

    suborder_w_n(1:suborder_num) = (/ &
      1051445, -2366706, 6493915, 1818134, -9986439,-3757007, 12368047, &
       478257, 10685542, -6437608 /)

     suborder_w_d = 3628800

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGLE_NCO_SUBRULE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal RULE = ', rule
    stop 1

  end if

  suborder_xyz(1:3,1:suborder_num) = &
      real ( 1 + suborder_xyz_n(1:3,1:suborder_num), kind = rk ) &
    / real ( 3 + suborder_xyz_d,                     kind = rk )

  suborder_w(1:suborder_num) = &
      real ( suborder_w_n(1:suborder_num), kind = rk ) &
    / real ( suborder_w_d,                 kind = rk )

  return
end
