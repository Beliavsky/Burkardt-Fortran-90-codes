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
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IUNIT, the free unit number.
!
  implicit none

  integer i
  integer ios
  integer iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

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
subroutine reference_to_physical_t4 ( tetra, n, ref, phy )

!*****************************************************************************80
!
!! REFERENCE_TO_PHYSICAL_T4 maps tetrahedron reference to physical points.
!
!  Discussion:
!
!    Given the vertices of an order 4 physical tetrahedron and a point
!    (R,S,T) in the reference tetrahedron, the routine computes the value
!    of the corresponding image point (X,Y,Z) in physical space.
!
!    This routine will also be correct for an order 10 tetrahedron,
!    if the mapping between reference and physical space
!    is linear.  This implies, in particular, that the sides of the
!    image tetrahedron are straight, the faces are flat, and
!    the "midside" nodes in the physical tetrahedron are
!    halfway along the edges of the physical tetrahedron.
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
!  Parameters:
!
!    Input, real ( kind = rk ) TETRA(3,4), the coordinates of the vertices.
!    The vertices are assumed to be the images of (0,0,0), (1,0,0),
!    (0,1,0) and (0,0,1) respectively.
!
!    Input, integer N, the number of points to transform.
!
!    Input, real ( kind = rk ) REF(3,N), points in the reference tetrahedron.
!
!    Output, real ( kind = rk ) PHY(3,N), corresponding points in the
!    physical tetrahedron.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  real ( kind = rk ) phy(3,n)
  real ( kind = rk ) ref(3,n)
  real ( kind = rk ) tetra(3,4)

  do i = 1, 3
    phy(i,1:n) = &
      tetra(i,1) * ( 1.0D+00 - ref(1,1:n) - ref(2,1:n) - ref(3,1:n) ) &
    + tetra(i,2) *             ref(1,1:n)                             &
    + tetra(i,3) *                          ref(2,1:n)                &
    + tetra(i,4) *                                       ref(3,1:n)
  end do

  return
end
subroutine tetrahedron_ncc_degree ( rule, degree )

!*****************************************************************************80
!
!! TETRAHEDRON_NCC_DEGREE returns the degree of an NCC rule.
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

  if ( 1 <= rule .and. rule <= 7 ) then

    degree = rule - 1

  else

    degree = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TETRAHEDRON_NCC_DEGREE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal RULE = ', rule
    stop 1

  end if

  return
end
subroutine tetrahedron_ncc_order_num ( rule, order_num )

!*****************************************************************************80
!
!! TETRAHEDRON_NCC_ORDER_NUM returns the order of an NCC rule.
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

  call tetrahedron_ncc_suborder_num ( rule, suborder_num )

  allocate ( suborder(1:suborder_num) )

  call tetrahedron_ncc_suborder ( rule, suborder_num, suborder )

  order_num = sum ( suborder(1:suborder_num) )

  deallocate ( suborder )

  return
end
subroutine tetrahedron_ncc_rule ( rule, order_num, xyz, w )

!*****************************************************************************80
!
!! TETRAHEDRON_NCC_RULE returns the points and weights of an NCC rule.
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
!    Input, integer ORDER_NUM, the order (number of points)
!     of the rule.
!
!    Output, real ( kind = rk ) XYZ(3,ORDER_NUM), the points of the rule.
!
!    Output, real ( kind = rk ) W(ORDER_NUM), the weights of the rule.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer order_num

  integer o
  integer rule
  integer s
  integer, allocatable, dimension ( : ) :: suborder
  integer suborder_num
  real ( kind = rk ), allocatable, dimension ( : ) :: suborder_w
  real ( kind = rk ), allocatable, dimension ( :, : ) :: suborder_xyz
  real ( kind = rk ) w(order_num)
  real ( kind = rk ) xyz(3,order_num)
!
!  Get the suborder information.
!
  call tetrahedron_ncc_suborder_num ( rule, suborder_num )

  allocate ( suborder(suborder_num) )
  allocate ( suborder_xyz(4,suborder_num) )
  allocate ( suborder_w(suborder_num) )

  call tetrahedron_ncc_suborder ( rule, suborder_num, suborder )

  call tetrahedron_ncc_subrule ( rule, suborder_num, suborder_xyz, suborder_w )
!
!  Expand the suborder information to a full order rule.
!
  o = 0

  do s = 1, suborder_num

    if ( suborder(s) == 1 ) then

      xyz(1:3,o+1) = suborder_xyz(1:3,s)
      w(o+1) = suborder_w(s)

      o = o + 1
!
!  Fourfold symmetry on (A,A,A,B)
!
!    123 AAA
!    124 AAB
!    142 ABA
!    412 BAA
!
    else if ( suborder(s) == 4 ) then

      xyz(1:3,o+1) = suborder_xyz( (/ 1, 2, 3 /), s )
      xyz(1:3,o+2) = suborder_xyz( (/ 1, 2, 4 /), s )
      xyz(1:3,o+3) = suborder_xyz( (/ 1, 4, 2 /), s )
      xyz(1:3,o+4) = suborder_xyz( (/ 4, 1, 2 /), s )
      w(o+1:o+4) = suborder_w(s)

      o = o + 4
!
!  Sixfold symmetry on (A,A,B,B):
!
!    123 (A,A,B)
!    132 (A,B,A),
!    134 (A,B,B)
!    312 (B,A,A)
!    314 (B,A,B)
!    341 (B,B,A)
!
    else if ( suborder(s) == 6 ) then

      xyz(1:3,o+1) = suborder_xyz( (/ 1, 2, 3 /), s )
      xyz(1:3,o+2) = suborder_xyz( (/ 1, 3, 2 /), s )
      xyz(1:3,o+3) = suborder_xyz( (/ 1, 3, 4 /), s )
      xyz(1:3,o+4) = suborder_xyz( (/ 3, 1, 2 /), s )
      xyz(1:3,o+5) = suborder_xyz( (/ 3, 1, 4 /), s )
      xyz(1:3,o+6) = suborder_xyz( (/ 3, 4, 1 /), s )
      w(o+1:o+6) = suborder_w(s)

      o = o + 6
!
!  Twelvefold symmetry on (A,A,B,C):
!
!    123 (A,A,B)
!    124 (A,A,C)
!    132 (A,B,A)
!    134 (A,B,C)
!    142 (A,C,A)
!    143 (A,C,B)
!    312 (B,A,A)
!    314 (B,A,C)
!    341 (B,C,A)
!    412 (C,A,A)
!    413 (C,A,B)
!    431 (C,B,A)
!
    else if ( suborder(s) == 12 ) then

      xyz(1:3,o+1)  = suborder_xyz( (/ 1, 2, 3 /), s )
      xyz(1:3,o+2)  = suborder_xyz( (/ 1, 2, 4 /), s )
      xyz(1:3,o+3)  = suborder_xyz( (/ 1, 3, 2 /), s )
      xyz(1:3,o+4)  = suborder_xyz( (/ 1, 3, 4 /), s )
      xyz(1:3,o+5)  = suborder_xyz( (/ 1, 4, 2 /), s )
      xyz(1:3,o+6)  = suborder_xyz( (/ 1, 4, 3 /), s )
      xyz(1:3,o+7)  = suborder_xyz( (/ 3, 1, 2 /), s )
      xyz(1:3,o+8)  = suborder_xyz( (/ 3, 1, 4 /), s )
      xyz(1:3,o+9)  = suborder_xyz( (/ 3, 4, 1 /), s )
      xyz(1:3,o+10) = suborder_xyz( (/ 4, 1, 2 /), s )
      xyz(1:3,o+11) = suborder_xyz( (/ 4, 1, 3 /), s )
      xyz(1:3,o+12) = suborder_xyz( (/ 4, 3, 1 /), s )
      w(o+1:o+12) = suborder_w(s)

      o = o + 12
!
!  24 fold symmetry on (A,B,C,D):
!
!    123 (A,B,C)
!    124 (A,B,D)
!    132 (A,C,B)
!    134 (A,C,D)
!    142 (A,D,B)
!    143 (A,D,C)
!    213 (B,A,C)
!    214 (B,A,D)
!    231 (B,C,A)
!    234 (B,C,D)
!    241 (B,D,A)
!    243 (B,D,C)
!    312 (C,A,B)
!    314 (C,A,D)
!    321 (C,B,A)
!    324 (C,B,D)
!    341 (C,D,A)
!    342 (C,D,B)
!    412 (D,A,B)
!    413 (D,A,C)
!    421 (D,B,A)
!    423 (D,B,C)
!    431 (D,C,A)
!    432 (D,C,B)
!
    else if ( suborder(s) == 24 ) then

      xyz(1:3,o+1)  = suborder_xyz( (/ 1, 2, 3 /), s )
      xyz(1:3,o+2)  = suborder_xyz( (/ 1, 2, 4 /), s )
      xyz(1:3,o+3)  = suborder_xyz( (/ 1, 3, 2 /), s )
      xyz(1:3,o+4)  = suborder_xyz( (/ 1, 3, 4 /), s )
      xyz(1:3,o+5)  = suborder_xyz( (/ 1, 4, 2 /), s )
      xyz(1:3,o+6)  = suborder_xyz( (/ 1, 4, 3 /), s )

      xyz(1:3,o+7)  = suborder_xyz( (/ 2, 1, 3 /), s )
      xyz(1:3,o+8)  = suborder_xyz( (/ 2, 1, 4 /), s )
      xyz(1:3,o+9)  = suborder_xyz( (/ 2, 3, 1 /), s )
      xyz(1:3,o+10) = suborder_xyz( (/ 2, 3, 4 /), s )
      xyz(1:3,o+11) = suborder_xyz( (/ 2, 4, 1 /), s )
      xyz(1:3,o+12) = suborder_xyz( (/ 2, 4, 3 /), s )

      xyz(1:3,o+13) = suborder_xyz( (/ 3, 1, 2 /), s )
      xyz(1:3,o+14) = suborder_xyz( (/ 3, 1, 4 /), s )
      xyz(1:3,o+15) = suborder_xyz( (/ 3, 2, 1 /), s )
      xyz(1:3,o+16) = suborder_xyz( (/ 3, 2, 4 /), s )
      xyz(1:3,o+17) = suborder_xyz( (/ 3, 4, 1 /), s )
      xyz(1:3,o+18) = suborder_xyz( (/ 3, 4, 2 /), s )

      xyz(1:3,o+19) = suborder_xyz( (/ 4, 1, 2 /), s )
      xyz(1:3,o+20) = suborder_xyz( (/ 4, 1, 3 /), s )
      xyz(1:3,o+21) = suborder_xyz( (/ 4, 2, 1 /), s )
      xyz(1:3,o+22) = suborder_xyz( (/ 4, 2, 3 /), s )
      xyz(1:3,o+23) = suborder_xyz( (/ 4, 3, 1 /), s )
      xyz(1:3,o+24) = suborder_xyz( (/ 4, 3, 2 /), s )
      w(o+1:o+24) = suborder_w(s)

      o = o + 24

    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TETRAHEDRON_NCC_RULE - Fatal error!'
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
subroutine tetrahedron_ncc_rule_num ( rule_num )

!*****************************************************************************80
!
!! TETRAHEDRON_NCC_RULE_NUM returns the number of NCC rules available.
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

  rule_num = 7

  return
end
subroutine tetrahedron_ncc_suborder ( rule, suborder_num, suborder )

!*****************************************************************************80
!
!! TETRAHEDRON_NCC_SUBORDER returns the suborders for an NCC rule.
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
!    Input, integer SUBORDER_NUM, the number of suborders
!    of the rule.
!
!    Output, integer SUBORDER(SUBORDER_NUM), the suborders
!    of the rule.
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
      4 /)
  else if ( rule == 3 ) then
    suborder(1:suborder_num) = (/ &
      4, 6 /)
  else if ( rule == 4 ) then
    suborder(1:suborder_num) = (/ &
      4, 12, 4 /)
  else if ( rule == 5 ) then
    suborder(1:suborder_num) = (/ &
      4, 12, 6, 12, 1 /)
  else if ( rule == 6 ) then
    suborder(1:suborder_num) = (/ &
      4, 12, 12, 12, 12, 4 /)
  else if ( rule == 7 ) then
    suborder(1:suborder_num) = (/ &
      4, 12, 12, 12, 6, 24, 4, 4, 6 /)
  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TETRAHEDRON_NCC_SUBORDER - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal RULE = ', rule
    stop 1

  end if

  return
end
subroutine tetrahedron_ncc_suborder_num ( rule, suborder_num )

!*****************************************************************************80
!
!! TETRAHEDRON_NCC_SUBORDER_NUM returns the number of suborders for an NCC rule.
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
!    Output, integer SUBORDER_NUM, the number of suborders
!    of the rule.
!
  implicit none

  integer rule
  integer, dimension(1:7) :: suborder = (/ &
     1, 1, 2, 3, 5, 6, 9 /)

  integer suborder_num

  if ( 1 <= rule .and. rule <= 7 ) then
    suborder_num = suborder(rule)
  else
    suborder_num = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TETRAHEDRON_NCC_SUBORDER_NUM - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal RULE = ', rule
    stop 1

  end if

  return
end
subroutine tetrahedron_ncc_subrule ( rule, suborder_num, suborder_xyz, &
  suborder_w )

!*****************************************************************************80
!
!! TETRAHEDRON_NCC_SUBRULE returns a compressed NCC rule.
!
!  Discussion:
!
!    In order for these compressed rules to be "unwrapped" correctly,
!    it's necessary that the values in SUBORDER_XYZ_N be listed
!    in a particular order for each kind of symmetry.  Basically,
!    the repeated equal values must come first.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 January 2007
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
!    Input, integer SUBORDER_NUM, the number of suborders
!    of the rule.
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
  real ( kind = rk ) suborder_xyz(4,suborder_num)
  integer suborder_xyz_n(4,suborder_num)
  integer suborder_xyz_d

  if ( rule == 1 ) then

    suborder_xyz_n(1:4,1:suborder_num) = reshape ( (/ &
      1, 1, 1, 1  &
      /), (/ 4, suborder_num /) )

    suborder_xyz_d = 4

    suborder_w_n(1:suborder_num) = (/ &
      1 /)

    suborder_w_d = 1

  else if ( rule == 2 ) then

    suborder_xyz_n(1:4,1:suborder_num) = reshape ( (/ &
      0, 0, 0, 1  &
      /), (/ 4, suborder_num /) )

    suborder_xyz_d = 1

    suborder_w_n(1:suborder_num) = (/ &
      1 /)

    suborder_w_d = 4

  else if ( rule == 3 ) then

    suborder_xyz_n(1:4,1:suborder_num) = reshape ( (/ &
      0, 0, 0, 2, &
      1, 1, 0, 0  &
      /), (/ 4, suborder_num /) )

    suborder_xyz_d = 2

    suborder_w_n(1:suborder_num) = (/ &
      -1, 4 /)

     suborder_w_d = 20

  else if ( rule == 4 ) then

    suborder_xyz_n(1:4,1:suborder_num) = reshape ( (/ &
      0, 0, 0, 3,  &
      0, 0, 1, 2,  &
      1, 1, 1, 0   &
      /), (/ 4, suborder_num /) )

    suborder_xyz_d = 3

    suborder_w_n(1:suborder_num) = (/ &
      1, 0, 9 /)

     suborder_w_d = 40

  else if ( rule == 5 ) then

    suborder_xyz_n(1:4,1:suborder_num) = reshape ( (/ &
      0, 0, 0, 4,  &
      0, 0, 3, 1,  &
      2, 2, 0, 0,  &
      1, 1, 0, 2,  &
      1, 1, 1, 1   &
      /), (/ 4, suborder_num /) )

    suborder_xyz_d = 4

    suborder_w_n(1:suborder_num) = (/ &
      -5, 16, -12, 16, 128 /)

     suborder_w_d = 420

  else if ( rule == 6 ) then

    suborder_xyz_n(1:4,1:suborder_num) = reshape ( (/ &
      0, 0, 0, 5,  &
      0, 0, 4, 1,  &
      0, 0, 3, 2,  &
      1, 1, 0, 3,  &
      2, 2, 1, 0,  &
      1, 1, 1, 2   &
      /), (/ 4, suborder_num /) )

    suborder_xyz_d = 5

    suborder_w_n(1:suborder_num) = (/ &
      33, -35, 35, 275, -75, 375 /)

     suborder_w_d = 4032

  else if ( rule == 7 ) then

    suborder_xyz_n(1:4,1:suborder_num) = reshape ( (/ &
      0, 0, 0, 6,  &
      0, 0, 5, 1,  &
      0, 0, 4, 2,  &
      1, 1, 0, 4,  &
      3, 3, 0, 0,  &
      3, 2, 1, 0,  &
      1, 1, 1, 3,  &
      2, 2, 2, 0,  &
      2, 2, 1, 1   &
      /), (/ 4, suborder_num /) )

    suborder_xyz_d = 6

    suborder_w_n(1:suborder_num) = (/ &
      -7, 24, -30, 0, 40, 30, 180, -45, 0 /)

     suborder_w_d = 1400

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TETRAHEDRON_NCC_SUBRULE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal RULE = ', rule
    stop 1

  end if

  suborder_xyz(1:4,1:suborder_num) = &
      real ( suborder_xyz_n(1:4,1:suborder_num), kind = rk ) &
    / real ( suborder_xyz_d,                     kind = rk )

  suborder_w(1:suborder_num) = &
      real ( suborder_w_n(1:suborder_num), kind = rk ) &
    / real ( suborder_w_d,                 kind = rk )

  return
end
subroutine tetrahedron_volume ( tetra, volume )

!*****************************************************************************80
!
!! TETRAHEDRON_VOLUME computes the volume of a tetrahedron in 3D.
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
