program main

!*****************************************************************************80
!
!! triangle_to_medit() converts mesh data from TRIANGLE to MEDIT format.
!
!  Usage:
!
!    triangle_to_medit prefix
!
!    where 'prefix' is the common filename prefix:
!
!    * 'prefix'.node contains the triangle node coordinates,
!    * 'prefix'.ele contains the triangle element node connectivity.
!    * 'prefix'.mesh will be the MESH file created by the program.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 October 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer arg_num
  integer dim
  integer edge_label(0)
  integer edge_vertex(0,0)
  integer edges
  real ( kind = rk ), allocatable :: element_att(:,:)
  integer element_att_num
  integer, allocatable :: element_node(:,:)
  integer element_num
  integer element_order
  integer hexahedron_label(0)
  integer hexahedron_vertex(0,0)
  integer hexahedrons
  integer iarg
  integer iargc
  integer m
  character ( len = 255 ) medit_filename
  real ( kind = rk ), allocatable :: node_att(:,:)
  integer node_att_num
  integer, allocatable :: node_marker(:,:)
  integer node_marker_num
  integer node_num;
  real ( kind = rk ), allocatable :: node_x(:,:)
  character ( len = 255 ) prefix
  integer quadrilateral_label(0)
  integer quadrilateral_vertex(0,0)
  integer quadrilaterals
  integer tetrahedron_label(0)
  integer tetrahedron_vertex(0,0)
  integer tetrahedrons
  character ( len = 255 ) triangle_element_filename
  integer, allocatable :: triangle_label(:)
  character ( len = 255 ) triangle_node_filename
  integer triangles
  integer vertices

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TRIANGLE_TO_MEDIT'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Read a mesh description created by TRIANGLE:'
  write ( *, '(a)' ) '  * "prefix".node, node coordinates.'
  write ( *, '(a)' ) '  * "prefix".ele, element connectivity.'
  write ( *, '(a)' ) '  Write a corresponding MEDIT mesh file.'
  write ( *, '(a)' ) '  * "prefix".mesh'
!
!  Get the number of command line arguments.
!
  arg_num = iargc ( )
!
!  Get the filename prefix.
!
  if ( arg_num < 1 ) then

    write ( *, '(a)' ) ''
    write ( *, '(a)' ) '  Please enter the filename prefix.'

    read ( *, '(a)' ) prefix

  else 

    iarg = 1
    call getarg ( iarg, prefix )

  end if
!
!  Create the filenames.
!
  triangle_node_filename = trim ( prefix ) // '.node'
  triangle_element_filename = trim ( prefix ) // '.ele'
  medit_filename = trim ( prefix ) // '.mesh'
!
!  Read TRIANGLE sizes.
!
  call triangle_node_size_read ( triangle_node_filename, node_num, m, &
    node_att_num, node_marker_num )

  call triangle_element_size_read ( triangle_element_filename, element_num, &
    element_order, element_att_num )
!
!  Report sizes.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Size information from TRIANGLE files:'
  write ( *, '(a,i4)' ) '  Spatial dimension M = ', m
  write ( *, '(a,i4)' ) '  Number of nodes NODE_NUM = ', node_num
  write ( *, '(a,i4)' ) '  NODE_ATT_NUM = ', node_att_num
  write ( *, '(a,i4)' ) '  NODE_MARKER_NUM = ', node_marker_num
  write ( *, '(a,i4)' ) '  Number of elements ELEMENT_NUM = ', element_num
  write ( *, '(a,i4)' ) '  Element order ELEMENT_ORDER = ', element_order
  write ( *, '(a,i4)' ) '  ELEMENT_ATT_NUM = ', element_att_num
!
!  Allocate memory.
!
  allocate ( node_att(1:node_att_num,1:node_num) )
  allocate ( node_marker(1:node_marker_num,1:node_num) )
  allocate ( node_x(m,node_num) )
  allocate ( element_node(1:element_order,1:element_num) )
  allocate ( element_att(1:element_att_num,1:element_num) )
!
!  Read TRIANGLE data.
!
  call triangle_node_data_read ( triangle_node_filename, node_num, m, &
    node_att_num, node_marker_num, node_x, node_att, node_marker )

  call triangle_element_data_read ( triangle_element_filename, element_num, &
    element_order, element_att_num, element_node, element_att )
!
!  Write the MEDIT data.
!
  dim = 2
  vertices = node_num
  edges = 0
  triangles = element_num
  quadrilaterals = 0
  tetrahedrons = 0
  hexahedrons = 0
  allocate ( triangle_label(1:triangles) )
  triangle_label(1:triangles) = 0

  call mesh_write ( medit_filename, dim, vertices, edges, triangles, &
    quadrilaterals, tetrahedrons, hexahedrons, node_x, node_marker, &
    edge_vertex, edge_label, element_node, triangle_label, &
    quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex, &
    tetrahedron_label, hexahedron_vertex, hexahedron_label )
!
!  Free memory.
!
  deallocate ( element_att )
  deallocate ( element_node )
  deallocate ( node_att )
  deallocate ( node_marker )
  deallocate ( node_x )
  deallocate ( triangle_label )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TRIANGLE_TO_MEDIT:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop
end
subroutine ch_cap ( c )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, the character to capitalize.
!
  implicit none

  character c
  integer itemp

  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end
function ch_eqi ( c1, c2 )

!*****************************************************************************80
!
!! CH_EQI is a case insensitive comparison of two characters for equality.
!
!  Example:
!
!    CH_EQI ( 'A', 'a' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C1, C2, the characters to compare.
!
!    Output, logical CH_EQI, the result of the comparison.
!
  implicit none

  logical ch_eqi
  character c1
  character c1_cap
  character c2
  character c2_cap

  c1_cap = c1
  c2_cap = c2

  call ch_cap ( c1_cap )
  call ch_cap ( c2_cap )

  if ( c1_cap == c2_cap ) then
    ch_eqi = .true.
  else
    ch_eqi = .false.
  end if

  return
end
subroutine ch_to_digit ( c, digit )

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the integer value of a base 10 digit.
!
!  Example:
!
!     C   DIGIT
!    ---  -----
!    '0'    0
!    '1'    1
!    ...  ...
!    '9'    9
!    ' '    0
!    'X'   -1
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the decimal digit, '0' through '9' or blank
!    are legal.
!
!    Output, integer DIGIT, the corresponding integer value.
!    If C was 'illegal', then DIGIT is -1.
!
  implicit none

  character c
  integer digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

    digit = ichar ( c ) - 48

  else if ( c == ' ' ) then

    digit = 0

  else

    digit = -1

  end if

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 September 2005
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
subroutine mesh_write ( filename, dim, vertices, edges, triangles, &
  quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, &
  vertex_label, edge_vertex, edge_label, triangle_vertex, triangle_label, &
  quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex, &
  tetrahedron_label, hexahedron_vertex, hexahedron_label )

!*****************************************************************************80
!
!! MESH_WRITE writes sizes and data to a MESH file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 December 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Pascal Frey,
!    MEDIT: An interactive mesh visualization software,
!    Technical Report RT-0253,
!    Institut National de Recherche en Informatique et en Automatique,
!    03 December 2001.
!
!  Parameters:
!
!    Input, character ( len = * ) FILENAME, the name of the file to be created.
!    Ordinarily, the name should include the extension ".mesh".
!
!    Input, integer DIM, the spatial dimension, which should
!    be 2 or 3.
!
!    Input, integer VERTICES, the number of vertices.
!
!    Input, real ( kind = rk ) VERTEX_COORDINATE(DIM,VERTICES), the coordinates
!    of each vertex.
!
!    Input, integer VERTEX_LABEL(VERTICES), a label for
!    each vertex.
!
!    Input, integer EDGES, the number of edges (may be 0).
!
!    Input, integer EDGE_VERTEX(2,EDGES), the vertices that form
!    each edge.
!
!    Input, integer EDGE_LABEL(EDGES), a label for each edge.
!
!    Input, integer TRIANGLES, the number of triangles (may be 0).
!
!    Input, integer TRIANGLE_VERTEX(3,TRIANGLES), the vertices
!    that form each triangle.
!
!    Input, integer TRIANGLE_LABEL(TRIANGLES), a label for each
!    triangle.
!
!    Input, integer QUADRILATERALS, the number of quadrilaterals
!    (may be 0).
!
!    Input, integer QUADRILATERAL_VERTEX(4,QUADRILATERALS), the
!    vertices that form each quadrilateral.
!
!    Input, integer QUADRILATERAL_LABEL(QUADRILATERALS), a label
!    for each quadrilateral.
!
!    Input, integer TETRAHEDRONS, the number of tetrahedrons
!    (may be 0).
!
!    Input, integer TETRAHEDRON_VERTEX(4,TETRAHEDRONS), the
!    vertices that form each tetrahedron.
!
!    Input, integer TETRAHEDRON_LABEL(TETRAHEDRONS), a label for
!    each tetrahedron.
!
!    Input, integer HEXAHEDRONS, the number of hexahedrons
!    (may be 0).
!
!    Input, integer HEXAHEDRON_VERTEX(8,HEXAHEDRONS), the vertices
!    that form each hexahedron.
!
!    Input, integer HEXAHEDRON_LABEL(HEXAHEDRONS), a label for
!    each hexahedron.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim
  integer edges
  integer hexahedrons
  integer quadrilaterals
  integer tetrahedrons
  integer triangles
  integer vertices

  integer edge_label(edges)
  integer edge_vertex(2,edges)
  character ( len = * ) filename
  integer fileunit
  integer hexahedron_label(hexahedrons)
  integer hexahedron_vertex(8,hexahedrons)
  integer ios
  integer j
  integer quadrilateral_label(quadrilaterals)
  integer quadrilateral_vertex(4,quadrilaterals)
  integer tetrahedron_label(tetrahedrons)
  integer tetrahedron_vertex(4,tetrahedrons)
  integer triangle_label(triangles)
  integer triangle_vertex(3,triangles)
  real ( kind = rk ) vertex_coordinate(dim,vertices)
  integer vertex_label(vertices)
!
!  Open the file.
!
  call get_unit ( fileunit )

  open ( unit = fileunit, file = filename, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MESH_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open file.'
    stop
  end if

  write ( fileunit, '(a)' ) 'MeshVersionFormatted 1'
  write ( fileunit, '(a)' ) '#  Created by mesh_write.f90'
!
!  Dimension information.
!
  write ( fileunit, '(a)' ) ' '
  write ( fileunit, '(a)' ) 'Dimension'
  write ( fileunit, '(i8)' ) dim
!
!  Vertices.
!
  write ( fileunit, '(a)' ) ' '
  write ( fileunit, '(a)' ) 'Vertices'
  write ( fileunit, '(i8)' ) vertices
  if ( dim == 2 ) then
    do j = 1, vertices
      write ( fileunit, '(2(2x,f10.6),2x,i8)' ) &
        vertex_coordinate(1:dim,j), vertex_label(j)
    end do
  else if ( dim == 3 ) then
    do j = 1, vertices
      write ( fileunit, '(3(2x,f10.6),2x,i8)' ) &
        vertex_coordinate(1:dim,j), vertex_label(j)
    end do
  end if
!
!  Edges.
!
  if ( 0 < edges ) then
    write ( fileunit, '(a)' ) ' '
    write ( fileunit, '(a)' ) 'Edges'
    write ( fileunit, '(i8)' ) edges
    do j = 1, edges
      write ( fileunit, '(2(2x,i8),2x,i8)' ) &
        edge_vertex(1:2,j), edge_label(j)
    end do
  end if
!
!  Triangles.
!
  if ( 0 < triangles ) then
    write ( fileunit, '(a)' ) ' '
    write ( fileunit, '(a)' ) 'Triangles'
    write ( fileunit, '(i8)' ) triangles
    do j = 1, triangles
      write ( fileunit, '(3(2x,i8),2x,i8)' ) &
        triangle_vertex(1:3,j), triangle_label(j)
    end do
  end if
!
!  Quadrilaterals.
!
  if ( 0 < quadrilaterals ) then
    write ( fileunit, '(a)' ) ' '
    write ( fileunit, '(a)' ) 'Quadrilaterals'
    write ( fileunit, '(i8)' ) quadrilaterals
    do j = 1, quadrilaterals
      write ( fileunit, '(4(2x,i8),2x,i8)' ) &
        quadrilateral_vertex(1:4,j), quadrilateral_label(j)
    end do
  end if
!
!  Tetrahedron.
!
  if ( 0 < tetrahedrons ) then
    write ( fileunit, '(a)' ) ' '
    write ( fileunit, '(a)' ) 'Tetrahedra'
    write ( fileunit, '(i8)' ) tetrahedrons
    do j = 1, tetrahedrons
      write ( fileunit, '(4(2x,i8),2x,i8)' ) &
        tetrahedron_vertex(1:4,j), tetrahedron_label(j)
    end do
  end if
!
!  Hexahedron.
!
  if ( 0 < hexahedrons ) then
    write ( fileunit, '(a)' ) ' '
    write ( fileunit, '(a)' ) 'Hexahedra'
    write ( fileunit, '(i8)' ) hexahedrons
    do j = 1, hexahedrons
      write ( fileunit, '(8(2x,i8),2x,i8)' ) &
        hexahedron_vertex(1:8,j), hexahedron_label(j)
    end do
  end if
!
!  End.
!
  write ( fileunit, '(a)' ) ' '
  write ( fileunit, '(a)' ) 'End'

  close ( unit = fileunit )

  return
end
subroutine s_to_i4 ( s, ival, ierror, length )

!*****************************************************************************80
!
!! S_TO_I4 reads an I4 from a string.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string to be examined.
!
!    Output, integer IVAL, the integer value read from the string.
!    If the string is blank, then IVAL will be returned 0.
!
!    Output, integer IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer LENGTH, the number of characters of S
!    used to make IVAL.
!
  implicit none

  character c
  integer i
  integer ierror
  integer isgn
  integer istate
  integer ival
  integer length
  character ( len = * ) s

  ierror = 0
  istate = 0
  isgn = 1
  ival = 0

  do i = 1, len_trim ( s )

    c = s(i:i)
!
!  Haven't read anything.
!
    if ( istate == 0 ) then

      if ( c == ' ' ) then

      else if ( c == '-' ) then
        istate = 1
        isgn = -1
      else if ( c == '+' ) then
        istate = 1
        isgn = + 1
      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read the sign, expecting digits.
!
    else if ( istate == 1 ) then

      if ( c == ' ' ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read at least one digit, expecting more.
!
    else if ( istate == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        ival = 10 * ival + ichar ( c ) - ichar ( '0' )
      else
        ival = isgn * ival
        length = i - 1
        return
      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( istate == 2 ) then
    ival = isgn * ival
    length = len_trim ( s )
  else
    ierror = 1
    length = 0
  end if

  return
end
subroutine s_to_i4vec ( s, n, ivec, ierror )

!*****************************************************************************80
!
!! S_TO_I4VEC reads an I4VEC from a string.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be read.
!
!    Input, integer N, the number of values expected.
!
!    Output, integer IVEC(N), the values read from the string.
!
!    Output, integer IERROR, error flag.
!    0, no errors occurred.
!    -K, could not read data for entries -K through N.
!
  implicit none

  integer n

  integer i
  integer ierror
  integer ilo
  integer ivec(n)
  integer length
  character ( len = * ) s

  i = 0
  ierror = 0
  ilo = 1

  do while ( i < n )

    i = i + 1

    call s_to_i4 ( s(ilo:), ivec(i), ierror, length )

    if ( ierror /= 0 ) then
      ierror = -i
      exit
    end if

    ilo = ilo + length

  end do

  return
end
subroutine s_to_r8 ( s, dval, ierror, length )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 from a string.
!
!  Discussion:
!
!    The routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 blanks
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semicolon,
!
!    with most quantities optional.
!
!  Example:
!
!    S                 DVAL
!
!    '1'               1.0
!    '     1   '       1.0
!    '1A'              1.0
!    '12,34,56'        12.0
!    '  34 7'          34.0
!    '-1E2ABCD'        -100.0
!    '-1X2ABCD'        -1.0
!    ' 2E-1'           0.2
!    '23.45'           23.45
!    '-4.2E+2'         -420.0
!    '17d2'            1700.0
!    '-14e-2'         -0.14
!    'e2'              100.0
!    '-12.73e-9.23'   -12.73 * 10.0^(-9.23)
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string containing the
!    data to be read.  Reading will begin at position 1 and
!    terminate at the end of the string, or when no more
!    characters can be read to form a legal real.  Blanks,
!    commas, or other nonnumeric data will, in particular,
!    cause the conversion to halt.
!
!    Output, real ( kind = rk ) DVAL, the value read from the string.
!
!    Output, integer IERROR, error flag.
!    0, no errors occurred.
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer LENGTH, the number of characters read
!    to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character c
  logical ch_eqi
  real ( kind = rk ) dval
  integer ierror
  integer ihave
  integer isgn
  integer iterm
  integer jbot
  integer jsgn
  integer jtop
  integer length
  integer nchar
  integer ndig
  real ( kind = rk ) rbot
  real ( kind = rk ) rexp
  real ( kind = rk ) rtop
  character ( len = * ) s

  nchar = len_trim ( s )

  ierror = 0
  dval = 0.0D+00
  length = -1
  isgn = 1
  rtop = 0
  rbot = 1
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

  do

    length = length + 1

    if ( nchar < length+1 ) then
      exit
    end if

    c = s(length+1:length+1)
!
!  Blank character.
!
    if ( c == ' ' ) then

      if ( ihave == 2 ) then

      else if ( ihave == 6 .or. ihave == 7 ) then
        iterm = 1
      else if ( 1 < ihave ) then
        ihave = 11
      end if
!
!  Comma.
!
    else if ( c == ',' .or. c == ';' ) then

      if ( ihave /= 1 ) then
        iterm = 1
        ihave = 12
        length = length + 1
      end if
!
!  Minus sign.
!
    else if ( c == '-' ) then

      if ( ihave == 1 ) then
        ihave = 2
        isgn = -1
      else if ( ihave == 6 ) then
        ihave = 7
        jsgn = -1
      else
        iterm = 1
      end if
!
!  Plus sign.
!
    else if ( c == '+' ) then

      if ( ihave == 1 ) then
        ihave = 2
      else if ( ihave == 6 ) then
        ihave = 7
      else
        iterm = 1
      end if
!
!  Decimal point.
!
    else if ( c == '.' ) then

      if ( ihave < 4 ) then
        ihave = 4
      else if ( 6 <= ihave .and. ihave <= 8 ) then
        ihave = 9
      else
        iterm = 1
      end if
!
!  Scientific notation exponent marker.
!
    else if ( ch_eqi ( c, 'E' ) .or. ch_eqi ( c, 'D' ) ) then

      if ( ihave < 6 ) then
        ihave = 6
      else
        iterm = 1
      end if
!
!  Digit.
!
    else if (  ihave < 11 .and. lle ( '0', c ) .and. lle ( c, '9' ) ) then

      if ( ihave <= 2 ) then
        ihave = 3
      else if ( ihave == 4 ) then
        ihave = 5
      else if ( ihave == 6 .or. ihave == 7 ) then
        ihave = 8
      else if ( ihave == 9 ) then
        ihave = 10
      end if

      call ch_to_digit ( c, ndig )

      if ( ihave == 3 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = rk )
      else if ( ihave == 5 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = rk )
        rbot = 10.0D+00 * rbot
      else if ( ihave == 8 ) then
        jtop = 10 * jtop + ndig
      else if ( ihave == 10 ) then
        jtop = 10 * jtop + ndig
        jbot = 10 * jbot
      end if
!
!  Anything else is regarded as a terminator.
!
    else
      iterm = 1
    end if
!
!  If we haven't seen a terminator, and we haven't examined the
!  entire string, go get the next character.
!
    if ( iterm == 1 ) then
      exit
    end if

  end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LENGTH is equal to NCHAR.
!
  if ( iterm /= 1 .and. length + 1 == nchar ) then
    length = nchar
  end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
  if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then
    ierror = ihave
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'S_TO_R8 - Serious error!'
    write ( *, '(a)' ) '  Illegal or nonnumeric input:'
    write ( *, '(a)' ) '    ' // trim ( s )
    return
  end if
!
!  Number seems OK.  Form it.
!
  if ( jtop == 0 ) then
    rexp = 1.0D+00
  else
    if ( jbot == 1 ) then
      rexp = 10.0D+00 ** ( jsgn * jtop )
    else
      rexp = 10.0D+00 ** ( real ( jsgn * jtop, kind = rk ) &
        / real ( jbot, kind = rk ) )
    end if
  end if

  dval = real ( isgn, kind = rk ) * rexp * rtop / rbot

  return
end
subroutine s_to_r8vec ( s, n, rvec, ierror )

!*****************************************************************************80
!
!! S_TO_R8VEC reads an R8VEC from a string.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be read.
!
!    Input, integer N, the number of values expected.
!
!    Output, real ( kind = rk ) RVEC(N), the values read from the string.
!
!    Output, integer IERROR, error flag.
!    0, no errors occurred.
!    -K, could not read data for entries -K through N.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer ierror
  integer ilo
  integer lchar
  real ( kind = rk ) rvec(n)
  character ( len = * ) s

  i = 0
  ierror = 0
  ilo = 1

  do while ( i < n )

    i = i + 1

    call s_to_r8 ( s(ilo:), rvec(i), ierror, lchar )

    if ( ierror /= 0 ) then
      ierror = -i
      exit
    end if

    ilo = ilo + lchar

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
subroutine triangle_element_data_read ( element_filename, element_num, &
  element_order, element_att_num, element_node, element_att )

!*****************************************************************************80
!
!! TRIANGLE_ELEMENT_DATA_READ reads the data from an element file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 October 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ELEMENT_FILENAME, the name of the file.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer ELEMENT_ORDER, the order of the elements.
!
!    Input, integer ELEMENT_ATT_NUM, number of element attributes 
!    listed on each node record.
!
!    Output, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM), 
!    the indices of the nodes that make up each element.
!
!    Output, real ( kind = rk ) ELEMENT_ATT(ELEMENT_ATT_NUM,ELEMENT_NUM), the 
!    attributes of each element.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer element_att_num
  integer element_num
  integer element_order

  integer element
  real ( kind = rk ) element_att(element_att_num,element_num)
  character ( len = * ) element_filename
  integer element_node(element_order,element_num)
  integer i
  integer i1
  integer i2
  integer i3
  integer ierror
  integer input
  integer ios
  integer ival
  integer length
  character ( len = 255 ) text
  real ( kind = rk ) value

  element = 0

  call get_unit ( input )

  open ( unit = input, file = element_filename, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'TRIANGLE_ELEMENT_DATA_READ - Fatal error!'
    write ( *, '(a)' ) '  Unable to open file.'
    stop 1
  end if

  do

    read ( input, '(a)', iostat = ios ) text

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TRIANGLE_ELEMENT_DATA_READ - Fatal error!'
      write ( *, '(a)' ) '  Unexpected end of file while reading.'
      stop 1
    end if

    if ( len_trim ( text ) == 0 ) then
      cycle
    end if

    if ( text(1:1) == '#' ) then
      cycle
    end if

    if ( element == 0 ) then

      call s_to_i4 ( text, i1, ierror, length )
      text = text(length+1:)
      call s_to_i4 ( text, i2, ierror, length )
      text = text(length+1:)
      call s_to_i4 ( text, i3, ierror, length )
      text = text(length+1:)

    else

      call s_to_i4 ( text, ival, ierror, length )
      text = text(length+1:)

      do i = 1, element_order
        call s_to_i4 ( text, ival, ierror, length )
        text = text(length+1:)
        element_node(i,element) = ival
      end do

      do i = 1, element_att_num
        call s_to_r8 ( text, value, length, ierror )
        text = text(length+1:)
        element_att(i,element) = value
      end do

    end if

    element = element + 1

    if ( element_num < element ) then
      exit
    end if

  end do

  close ( unit = input )

  return
end
subroutine triangle_element_size_read ( element_filename, element_num, &
  element_order, element_att_num )

!*****************************************************************************80
!
!! TRIANGLE_ELEMENT_SIZE_READ reads the header information from an element file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 October 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ELEMENT_FILENAME, the name of the 
!    element file.
!
!    Output, integer ELEMENT_NUM, the number of elements.
!
!    Output, integer ELEMENT_ORDER, the order of the elements.
!
!    Output, integer ELEMENT_ATT_NUM, the number of 
!    element attributes.
!
  implicit none

  integer element_att_num
  character ( len = * ) element_filename
  integer element_num
  integer element_order
  integer ierror
  integer input
  integer ios
  integer length
  character ( len = 255 ) text

  call get_unit ( input )

  open ( unit = input, file = element_filename, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'TRIANGLE_ELEMENT_SIZE_READ - Fatal error!'
    write ( *, '(a)' ) '  Unable to open file.'
    stop 1
  end if

  do

    read ( input, '(a)', iostat = ios ) text

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TRIANGLE_ELEMENT_SIZE_READ - Fatal error!'
      write ( *, '(a)' ) '  Unexpected end of file while reading.'
      stop 1
    end if

    if ( len_trim ( text ) == 0 ) then
      cycle
    end if

    if ( text(1:1) == '#' ) then
      cycle
    end if

    call s_to_i4 ( text, element_num, ierror, length )
    text = text(length+1:)
    call s_to_i4 ( text, element_order, ierror, length )
    text = text(length+1:)
    call s_to_i4 ( text, element_att_num, ierror, length )
    text = text(length+1:)

    exit

  end do

  close ( unit = input )

  return
end
subroutine triangle_node_data_read ( node_filename, node_num, node_dim, &
  node_att_num, node_marker_num, node_coord, node_att, node_marker )

!*****************************************************************************80
!
!! TRIANGLE_NODE_DATA_READ reads the data from a node file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 October 2014
!
!  Author:
!
!    John Burkardt
!
! Parameters:
!
!   Input, character ( len = * ) NODE_FILENAME, the name of the node file.
!
!   Input, integer NODE_NUM, the number of nodes.
!
!   Input, integer NODE_DIM, the spatial dimension.
!
!   Input, integer NODE_ATT_NUM, number of node attributes 
!   listed on each node record.
!
!   Input, integer NODE_MARKER_NUM, 1 if every node record 
!   includes a final boundary marker value.
!
!   Output, real ( kind = rk ) NODE_COORD(NODE_DIM,NODE_NUM), the nodal 
!   coordinates.
!
!   Output, real ( kind = rk ) NODE_ATT(NODE_ATT_NUM,NODE_NUM), the nodal 
!   attributes.
!
!   Output, integer NODE_MARKER(NODE_MARKER_NUM,NODE_NUM), the 
!   node markers.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer node_att_num
  integer node_dim
  integer node_marker_num
  integer node_num

  integer i
  integer ierror
  integer input
  integer ios
  integer ival
  integer length
  integer node
  real ( kind = rk ) node_att(node_att_num,node_num)
  real ( kind = rk ) node_coord(node_dim,node_num)
  character ( len = * ) node_filename
  integer node_marker(node_marker_num,node_num)
  character ( len = 255 ) text
  real ( kind = rk ) value

  node = 0

  call get_unit ( input )

  open ( unit = input, file = node_filename, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'TRIANGLE_NODE_DATA_READ - Fatal error!'
    write ( *, '(a)' ) '  Unable to open file.'
    stop 1
  end if

  do

    read ( input, '(a)', iostat = ios ) text

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TRIANGLE_NODE_DATA_READ - Fatal error!'
      write ( *, '(a)' ) '  Unexpected end of file while reading.'
      stop 1
    end if

    if ( len_trim ( text ) == 0 ) then
      cycle
    end if

    if ( text(1:1) == '#' ) then
      cycle
    end if
!
!  Ignore the dimension line.
!
    if ( node == 0 ) then

    else

      call s_to_i4 ( text, ival, ierror, length )
      text = text(length+1:)

      do i = 1, node_dim
        call s_to_r8 ( text, value, ierror, length )
        text = text(length+1:)
        node_coord(i,node) = value
      end do

      do i = 1, node_att_num
        call s_to_r8 ( text, value, ierror, length )
        text = text(length+1:)
        node_att(i,node) = value;
      end do

      do i = 1, node_marker_num
        call s_to_i4 ( text, ival, ierror, length )
        text = text(length+1:)
        node_marker(i,node) = ival
      end do

    end if

    node = node + 1

    if ( node_num < node ) then
      exit
    end if

  end do

  close ( unit = input )

  return
end
subroutine triangle_node_size_read ( node_filename, node_num, node_dim, &
  node_att_num, node_marker_num )

!*****************************************************************************80
!
!! TRIANGLE_NODE_SIZE_READ reads the header information from a node file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 October 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) NODE_FILENAME, the name of the node file.
!
!    Output, integer NODE_NUM, the number of nodes.
!
!    Output, integer NODE_DIM, the spatial dimension.
!
!    Output, integer NODE_ATT_NUM, number of node attributes 
!    listed on each node record.
!
!    Output, integer NODE_MARKER_NUM, 1 if every node record 
!    includes a final boundary marker value.
!
  implicit none

  integer ierror
  integer input
  integer ios
  integer length
  integer node_att_num
  integer node_dim
  character ( len = * ) node_filename
  integer node_marker_num
  integer node_num
  character ( len = 255 ) text

  call get_unit ( input )

  open ( unit = input, file = node_filename, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'TRIANGLE_NODE_SIZE_READ - Fatal error!'
    write ( *, '(a)' ) '  Unable to open file.'
    stop 1
  end if

  do

    read ( input, '(a)', iostat = ios ) text

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TRIANGLE_NODE_SIZE_READ - Fatal error!'
      write ( *, '(a)' ) '  Unexpected end of file while reading.'
      stop 1
    end if

    if ( len_trim ( text ) == 0 ) then
      cycle
    end if

    if ( text(1:1) == '#' ) then
      cycle
    end if

    call s_to_i4 ( text, node_num, ierror, length )
    text = text(length+1:)
    call s_to_i4 ( text, node_dim, ierror, length )
    text = text(length+1:)
    call s_to_i4 ( text, node_att_num, ierror, length )
    text = text(length+1:)
    call s_to_i4 ( text, node_marker_num, ierror, length )
    text = text(length+1:)

    exit

  end do

  close ( unit = input )

  return
end
