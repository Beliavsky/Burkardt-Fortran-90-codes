program main

!*****************************************************************************80
!
!! medit_to_fem() converts mesh data from MEDIT to FEM format.
!
!  Usage:
!
!    medit_to_fem prefix
!
!    where 'prefix' is the common filename prefix:
!
!    * 'prefix'.mesh is the MEDIT filename.
!    * 'prefix'_nodes.txt will contain the FEM node coordinates.
!    * 'prefix'_elements.txt will contain the FEM element node connectivity.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 November 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer arg_num
  integer dim
  integer, allocatable :: edge_label(:)
  integer, allocatable :: edge_vertex(:,:)
  integer edges
  integer element_num
  integer element_order
  character ( len = 255 ) fem_element_filename
  character ( len = 255 ) fem_node_filename
  integer, allocatable :: hexahedron_label(:)
  integer, allocatable :: hexahedron_vertex(:,:)
  integer hexahedrons
  integer iarg
  integer iargc
  integer m
  character ( len = 255 ) medit_filename
  integer node_num;
  character ( len = 255 ) prefix
  integer, allocatable :: quadrilateral_label(:)
  integer, allocatable :: quadrilateral_vertex(:,:)
  integer quadrilaterals
  integer, allocatable :: tetrahedron_label(:)
  integer, allocatable :: tetrahedron_vertex(:,:)
  integer tetrahedrons
  integer, allocatable :: triangle_label(:)
  integer, allocatable :: triangle_vertex(:,:)
  integer triangles
  real ( kind = rk ), allocatable :: vertex_coordinate(:,:)
  integer, allocatable :: vertex_label(:)
  integer vertices

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'medit_to_fem():'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Read a MESH file created by MEDIT:'
  write ( *, '(a)' ) '  * "prefix".mesh.'
  write ( *, '(a)' ) '  Write two simple FEM files listing nodes and elements.'
  write ( *, '(a)' ) '  * "prefix"_nodes.txt, node coordinates.'
  write ( *, '(a)' ) '  * "prefix"_elements.txt, element connectivity.'
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
  medit_filename = trim ( prefix ) // '.mesh'
  fem_node_filename = trim ( prefix ) // '_nodes.txt'
  fem_element_filename = trim ( prefix ) // '_elements.txt'
!
!  Read MEDIT sizes.
!
  call mesh_size_read ( medit_filename, dim, vertices, edges, triangles, &
  quadrilaterals, tetrahedrons, hexahedrons )
!
!  Report sizes.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Size information from MEDIT file:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of dimensions = ', dim
  write ( *, '(a,i8)' ) '  Number of vertices = ', vertices
  write ( *, '(a,i8)' ) '  Number of edges = ', edges
  write ( *, '(a,i8)' ) '  Number of triangles = ', triangles
  write ( *, '(a,i8)' ) '  Number of quadrilaterals = ', quadrilaterals
  write ( *, '(a,i8)' ) '  Number of tetrahedrons = ', tetrahedrons
  write ( *, '(a,i8)' ) '  Number of hexahedrons = ', hexahedrons
!
!  Allocate memory.
!
  allocate ( edge_label(edges) )
  allocate ( edge_vertex(2,edges) )
  allocate ( hexahedron_label(hexahedrons) )
  allocate ( hexahedron_vertex(8,hexahedrons) )
  allocate ( quadrilateral_label(quadrilaterals) )
  allocate ( quadrilateral_vertex(4,quadrilaterals) )
  allocate ( tetrahedron_label(tetrahedrons) )
  allocate ( tetrahedron_vertex(4,tetrahedrons) )
  allocate ( triangle_label(triangles) )
  allocate ( triangle_vertex(3,triangles) )
  allocate ( vertex_coordinate(dim,vertices) )
  allocate ( vertex_label(vertices) )
!
!  Read MEDIT data.
!
  call mesh_data_read ( medit_filename, dim, vertices, edges, triangles, &
    quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, &
    vertex_label, edge_vertex, edge_label, triangle_vertex, triangle_label, &
    quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex, &
    tetrahedron_label, hexahedron_vertex, hexahedron_label )
!
!  Choose the FEM data.
!
!  We need to assume that there is only one element type.
!  If there are elements of multiple dimension, take the highest.
!
  m = dim
  node_num = vertices
  call r8mat_write ( fem_node_filename, dim, vertices, vertex_coordinate )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Created node coordinate file "' // &
    trim ( fem_node_filename ) // '"'

  if ( 0 < hexahedrons .and. dim == 3 ) then

    element_order = 8
    element_num = hexahedrons
    call i4mat_write ( fem_element_filename, element_order, element_num, &
      hexahedron_vertex )

  else if ( 0 < tetrahedrons .and. dim == 3 ) then

    element_order = 4
    element_num = tetrahedrons
    call i4mat_write ( fem_element_filename, element_order, element_num, &
      tetrahedron_vertex )

  else if ( 0 < quadrilaterals .and. dim == 2 ) then

    element_order = 4
    element_num = quadrilaterals
    call i4mat_write ( fem_element_filename, element_order, element_num, &
      quadrilateral_vertex )

  else if ( 0 < triangles .and. dim == 2 ) then

    element_order = 3
    element_num = triangles
    call i4mat_write ( fem_element_filename, element_order, element_num, &
      triangle_vertex )

  else

    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'medit_to_fem(): Fatal error!'
    write ( *, '(a)' ) '  Unexpected combination of spatial dimension'
    write ( *, '(a)' ) '  and number of nonzero objects.'
    stop 1

  end if

  write ( *, '(a)' ) '  Created element connectivity file "' // &
    trim ( fem_element_filename ) // '"'
!
!  Free memory.
!
  deallocate ( edge_label )
  deallocate ( edge_vertex )
  deallocate ( hexahedron_label )
  deallocate ( hexahedron_vertex )
  deallocate ( quadrilateral_label )
  deallocate ( quadrilateral_vertex )
  deallocate ( tetrahedron_label )
  deallocate ( tetrahedron_vertex )
  deallocate ( triangle_label )
  deallocate ( triangle_vertex )
  deallocate ( vertex_coordinate )
  deallocate ( vertex_label )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'medit_to_fem():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop
end
subroutine ch_cap ( c )

!*****************************************************************************80
!
!! CH_CAP() capitalizes a single character.
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
subroutine i4mat_write ( output_filename, m, n, table )

!*****************************************************************************80
!
!! I4MAT_WRITE writes an I4MAT file.
!
!  Discussion:
!
!    An I4MAT is an array of I4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the output file name.
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of points.
!
!    Input, integer TABLE(M,N), the data.
!
  implicit none

  integer m
  integer n

  integer j
  character ( len = * ) output_filename
  integer output_status
  integer output_unit
  character ( len = 30 ) string
  integer table(m,n)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4MAT_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_filename ) // '" on unit ', output_unit
    output_unit = -1
    stop 1
  end if
!
!  Create a format string.
!
  if ( 0 < m .and. 0 < n ) then

    write ( string, '(a1,i8,a4)' ) '(', m, 'i10)'
!
!  Write the data.
!
    do j = 1, n
      write ( output_unit, string ) table(1:m,j)
    end do

  end if
!
!  Close the file.
!
  close ( unit = output_unit )

  return
end
subroutine mesh_data_read ( filename, dim, vertices, edges, triangles, &
  quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, &
  vertex_label, edge_vertex, edge_label, triangle_vertex, triangle_label, &
  quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex, &
  tetrahedron_label, hexahedron_vertex, hexahedron_label )

!*****************************************************************************80
!
!! MESH_DATA_READ reads data from a MESH file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 October 2010
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
!    Input, character ( len = * ) FILENAME, the name of the MESH file.
!
!    Input, integer DIM, the spatial dimension, which should
!    be 2 or 3.
!
!    Input, integer VERTICES, the number of vertices.
!
!    Input, integer EDGES, the number of edges (may be 0).
!
!    Input, integer TRIANGLES, the number of triangles (may be 0).
!
!    Input, integer QUADRILATERALS, the number of quadrilaterals
!    (may be 0).
!
!    Input, integer TETRAHEDRONS, the number of tetrahedrons
!    (may be 0).
!
!    Input, integer HEXAHEDRONS, the number of hexahedrons
!    (may be 0).
!
!    Output, real VERTEX_COORDINATE(DIM,VERTICES), the coordinates
!    of each vertex.
!
!    Output, integer VERTEX_LABEL(VERTICES), a label for
!    each vertex.
!
!    Output, integer EDGE_VERTEX(2,EDGES), the vertices that form
!    each edge.
!
!    Output, integer EDGE_LABEL(EDGES), a label for each edge.
!
!    Output, integer TRIANGLE_VERTEX(3,TRIANGLES), the vertices
!    that form each triangle.
!
!    Output, integer TRIANGLE_LABEL(TRIANGLES), a label for each
!    triangle.
!
!    Output, integer QUADRILATERAL_VERTEX(4,QUADRILATERALS), the
!    vertices that form each quadrilateral.
!
!    Output, integer QUADRILATERAL_LABEL(QUADRILATERALS), a label
!    for each quadrilateral.
!
!    Output, integer TETRAHEDRON_VERTEX(4,TETRAHEDRONS), the
!    vertices that form each tetrahedron.
!
!    Output, integer TETRAHEDRON_LABEL(TETRAHEDRONS), a label for
!    each tetrahedron.
!
!    Output, integer HEXAHEDRON_VERTEX(8,HEXAHEDRONS), the vertices
!    that form each hexahedron.
!
!    Output, integer HEXAHEDRON_LABEL(HEXAHEDRONS), a label for
!    each hexahedron.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer edges
  integer hexahedrons
  integer quadrilaterals
  integer tetrahedrons
  integer triangles
  integer vertices

  integer dim
  integer edge
  integer edge_label(edges)
  integer edge_vertex(2,edges)
  character ( len = * ) filename
  integer fileunit
  integer hexahedron
  integer hexahedron_label(hexahedrons)
  integer hexahedron_vertex(8,hexahedrons)
  integer i4vec(9)
  integer ierror
  integer ios
  character ( len = 80 ) keyword
  integer length
  integer line_num
  integer quadrilateral
  integer quadrilateral_label(quadrilaterals)
  integer quadrilateral_vertex(4,quadrilaterals)
  real ( kind = rk ) r8vec(9)
  logical s_begin
  logical s_eqi
  integer tetrahedron
  integer tetrahedron_label(tetrahedrons)
  integer tetrahedron_vertex(4,tetrahedrons)
  character ( len = 255 ) text
  integer triangle
  integer triangle_label(triangles)
  integer triangle_vertex(3,triangles)
  integer vertex
  real ( kind = rk ) vertex_coordinate(dim,vertices)
  integer vertex_label(vertices)
!
!  Initialize everything to nothing.
!
  vertex_coordinate(1:dim,1:vertices) = 0.0D+00
  vertex_label(1:vertices) = 0
  edge_vertex(1:2,1:edges) = 0
  edge_label(1:edges) = 0
  triangle_vertex(1:3,1:triangles) = 0
  triangle_label(1:triangles) = 0
  quadrilateral_vertex(1:4,1:quadrilaterals) = 0
  quadrilateral_label(1:quadrilaterals) = 0
  tetrahedron_vertex(1:4,1:tetrahedrons) = 0
  tetrahedron_label(1:tetrahedrons) = 0
  hexahedron_vertex(1:8,1:hexahedrons) = 0
  hexahedron_label(1:hexahedrons) = 0
!
!  Open the file.
!
  call get_unit ( fileunit )

  open ( unit = fileunit, file = filename, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MESH_DATA_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open file.'
    stop
  end if
!
!  Read lines til you get alphanumerics and determine a "mode"
!
  line_num = 0
  keyword = 'NONE'

  do

    read ( fileunit, '(a)', iostat = ios ) text

    if ( ios /= 0 ) then
      exit
    end if

    line_num = line_num + 1

    if ( len_trim ( text ) == 0 ) then
      keyword = 'NONE'
      cycle
    end if

    if ( text(1:1) == '#' ) then
      cycle
    end if
!
!  Remove initial blanks.
!
    text = adjustl ( text )
!
!  Expecting a keyword.
!
        if ( s_eqi ( text, 'CORNERS' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'DIMENSION' ) ) then

      keyword = 'DIMENSION'

    else if ( s_eqi ( text, 'EDGES' ) ) then

      keyword = 'EDGES'

    else if ( s_eqi ( text, 'END' ) ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  END statement encountered.'
      exit

    else if ( s_eqi ( text, 'HEXAHEDRA' ) .or. &
              s_eqi ( text, 'HEXAHEDRONS' ) ) then

      keyword = 'HEXAHEDRONS'

    else if ( s_begin ( text, 'MESHVERSIONFORMATTED' ) ) then

    else if ( s_eqi ( text, 'NORMALATQUADRILATERALVERTICES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'NORMALATTRIANGLEVERTICES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'NORMALATVERTICES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'NORMALS' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'QUADRILATERALS' ) ) then

      keyword = 'QUADRILATERALS'

    else if ( s_eqi ( text, 'REQUIREDEDGES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'REQUIREDVERTICES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'RIDGES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'TANGENTATEDGES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'TANGENTS' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'TETRAHEDRA' ) .or. &
              s_eqi ( text, 'TETRAHEDRONS' ) ) then

      keyword = 'TETRAHEDRONS'

    else if ( s_eqi ( text, 'TRIANGLES' ) ) then

      keyword = 'TRIANGLES'

    else if ( s_eqi ( text, 'VERTICES' ) ) then

      keyword = 'VERTICES'
!
!  Presumably, numeric data to be processed by keyword.
!
    else if ( s_eqi ( keyword, 'DIMENSION' ) ) then

      call s_to_i4 ( text, dim, ierror, length )

      keyword = 'NONE'

    else if ( s_eqi ( keyword, 'EDGES' ) ) then

      call s_to_i4 ( text, edges, ierror, length )

      keyword = 'EDGE_VERTEX'
      edge = 0

    else if ( s_eqi ( keyword, 'EDGE_VERTEX' ) ) then

      call s_to_i4vec ( text, 3, i4vec, ierror )
      edge = edge + 1
      edge_vertex(1:2,edge) = i4vec(1:2)
      edge_label(edge) = i4vec(3)

    else if ( s_eqi ( keyword, 'HEXAHEDRONS' ) ) then

      call s_to_i4 ( text, hexahedrons, ierror, length )

      keyword = 'HEXAHEDRON_VERTEX'
      hexahedron = 0

    else if ( s_eqi ( keyword, 'HEXAHEDRON_VERTEX' ) ) then

      call s_to_i4vec ( text, 9, i4vec, ierror )
      hexahedron = hexahedron + 1
      hexahedron_vertex(1:8,hexahedron) = i4vec(1:8)
      hexahedron_label(hexahedron) = i4vec(9)

    else if ( s_eqi ( keyword, 'QUADRILATERALS' ) ) then

      call s_to_i4 ( text, quadrilaterals, ierror, length )

      keyword = 'QUADRILATERAL_VERTEX'
      quadrilateral = 0

    else if ( s_eqi ( keyword, 'QUADRILATERAL_VERTEX' ) ) then

      call s_to_i4vec ( text, 5, i4vec, ierror )
      quadrilateral = quadrilateral + 1
      quadrilateral_vertex(1:4,quadrilateral) = i4vec(1:4)
      quadrilateral_label(quadrilateral) = i4vec(5)

    else if ( s_eqi ( keyword, 'TETRAHEDRONS' ) ) then

      call s_to_i4 ( text, tetrahedrons, ierror, length )

      keyword = 'TETRAHEDRON_VERTEX'
      tetrahedron = 0

    else if ( s_eqi ( keyword, 'TETRAHEDRON_VERTEX' ) ) then

      call s_to_i4vec ( text, 5, i4vec, ierror )
      tetrahedron = tetrahedron + 1
      tetrahedron_vertex(1:4,tetrahedron) = i4vec(1:4)
      tetrahedron_label(tetrahedron) = i4vec(5)

    else if ( s_eqi ( keyword, 'TRIANGLES' ) ) then

      call s_to_i4 ( text, triangles, ierror, length )

      keyword = 'TRIANGLE_VERTEX'
      triangle = 0

    else if ( s_eqi ( keyword, 'TRIANGLE_VERTEX' ) ) then

      call s_to_i4vec ( text, 4, i4vec, ierror )
      triangle = triangle + 1
      triangle_vertex(1:3,triangle) = i4vec(1:3)
      triangle_label(triangle) = i4vec(4)

    else if ( s_eqi ( keyword, 'VERTICES' ) ) then

      call s_to_i4 ( text, vertices, ierror, length )

      keyword = 'VERTEX_COORDINATE'
      vertex = 0

    else if ( s_eqi ( keyword, 'VERTEX_COORDINATE' ) ) then

      call s_to_r8vec ( text, dim + 1, r8vec, ierror )
      vertex = vertex + 1
      vertex_coordinate(1:dim,vertex) = r8vec(1:dim)
      vertex_label(vertex) = int ( r8vec(dim+1) )

    else if ( s_eqi ( keyword, 'SKIP' ) ) then

    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MESH_DATA_READ - Fatal error!'
      write ( *, '(a,i8)' ) &
        '  Could not find keyword while reading line ', line_num
      write ( *, '(a)' ) '"' // trim ( text ) // '".'
      stop

    end if

  end do
!
!  Close the file.
!
  close ( unit = fileunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) '  Read ', line_num, &
    ' lines from "' // trim ( filename ) // '".'

  return
end
subroutine mesh_size_read ( filename, dim, vertices, edges, triangles, &
  quadrilaterals, tetrahedrons, hexahedrons )

!*****************************************************************************80
!
!! MESH_SIZE_READ reads sizes from a MESH file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 October 2010
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
!    Input, character ( len = * ) FILENAME, the name of the MESH file.
!
!    Output, integer DIM, the spatial dimension, 
!    which should be 2 or 3.
!
!    Output, integer VERTICES, the number of vertices.
!
!    Output, integer EDGES, the number of edges (may be 0).
!
!    Output, integer TRIANGLES, the number of triangles (may be 0).
!
!    Output, integer QUADRILATERALS, the number of quadrilaterals
!    (may be 0).
!
!    Output, integer TETRAHEDRONS, the number of tetrahedrons
!    (may be 0).
!
!    Output, integer HEXAHEDRONS, the number of hexahedrons
!    (may be 0).
!
  implicit none

  integer dim
  integer edges
  character ( len = * ) filename
  integer fileunit
  integer hexahedrons
  integer ierror
  integer ios
  character ( len = 80 ) keyword
  integer length
  integer line_num
  integer quadrilaterals
  logical s_begin
  logical s_eqi
  integer tetrahedrons
  character ( len = 255 ) text
  integer triangles
  integer vertices
!
!  Initialize everything to nothing.
!
  dim = 0
  vertices = 0
  edges = 0
  triangles = 0
  quadrilaterals = 0
  tetrahedrons = 0
  hexahedrons = 0
!
!  Open the file.
!
  call get_unit ( fileunit )

  open ( unit = fileunit, file = filename, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MESH_SIZE_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open file.'
    stop
  end if
!
!  Read lines til you get alphanumerics and determine a "mode"
!
  line_num = 0
  keyword = 'NONE'

  do

    read ( fileunit, '(a)', iostat = ios ) text

    if ( ios /= 0 ) then
      exit
    end if

    line_num = line_num + 1

    if ( len_trim ( text ) == 0 ) then
      keyword = 'NONE'
      cycle
    end if

    if ( text(1:1) == '#' ) then
      cycle
    end if
!
!  Remove initial blanks.
!
    text = adjustl ( text )
!
!  Expecting a keyword.
!
        if ( s_eqi ( text, 'CORNERS' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'DIMENSION' ) ) then

      keyword = 'DIMENSION'

    else if ( s_eqi ( text, 'EDGES' ) ) then

      keyword = 'EDGES'

    else if ( s_eqi ( text, 'END' ) ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  END statement encountered.'
      exit

    else if ( s_eqi ( text, 'HEXAHEDRA' ) .or. &
              s_eqi ( text, 'HEXAHEDRONS' ) ) then

      keyword = 'HEXAHEDRONS'

    else if ( s_begin ( text, 'MESHVERSIONFORMATTED' ) ) then

    else if ( s_eqi ( text, 'NORMALATQUADRILATERALVERTICES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'NORMALATTRIANGLEVERTICES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'NORMALATVERTICES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'NORMALS' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'QUADRILATERALS' ) ) then

      keyword = 'QUADRILATERALS'

    else if ( s_eqi ( text, 'REQUIREDEDGES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'REQUIREDVERTICES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'RIDGES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'TANGENTATEDGES' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'TANGENTS' ) ) then

      keyword = 'SKIP'

    else if ( s_eqi ( text, 'TETRAHEDRA' ) .or. &
              s_eqi ( text, 'TETRAHEDRONS' ) ) then

      keyword = 'TETRAHEDRONS'

    else if ( s_eqi ( text, 'TRIANGLES' ) ) then

      keyword = 'TRIANGLES'

    else if ( s_eqi ( text, 'VERTICES' ) ) then

      keyword = 'VERTICES'
!
!  Presumably, numeric data to be processed by keyword.
!
    else if ( s_eqi ( keyword, 'DIMENSION' ) ) then

      call s_to_i4 ( text, dim, ierror, length )

      keyword = 'NONE'

    else if ( s_eqi ( keyword, 'EDGES' ) ) then

      call s_to_i4 ( text, edges, ierror, length )

      keyword = 'EDGE_VERTEX'

    else if ( s_eqi ( keyword, 'EDGE_VERTEX' ) ) then

    else if ( s_eqi ( keyword, 'HEXAHEDRONS' ) ) then

      call s_to_i4 ( text, hexahedrons, ierror, length )

      keyword = 'HEXAHEDRON_VERTEX'

    else if ( s_eqi ( keyword, 'HEXAHEDRON_VERTEX' ) ) then

    else if ( s_eqi ( keyword, 'QUADRILATERALS' ) ) then

      call s_to_i4 ( text, quadrilaterals, ierror, length )

      keyword = 'QUADRILATERAL_VERTEX'

    else if ( s_eqi ( keyword, 'QUADRILATERAL_VERTEX' ) ) then

    else if ( s_eqi ( keyword, 'TETRAHEDRONS' ) ) then

      call s_to_i4 ( text, tetrahedrons, ierror, length )

      keyword = 'TETRAHEDRON_VERTEX'

    else if ( s_eqi ( keyword, 'TETRAHEDRON_VERTEX' ) ) then

    else if ( s_eqi ( keyword, 'TRIANGLES' ) ) then

      call s_to_i4 ( text, triangles, ierror, length )

      keyword = 'TRIANGLE_VERTEX'

    else if ( s_eqi ( keyword, 'TRIANGLE_VERTEX' ) ) then

    else if ( s_eqi ( keyword, 'VERTICES' ) ) then

      call s_to_i4 ( text, vertices, ierror, length )

      keyword = 'VERTEX_COORDINATE'

    else if ( s_eqi ( keyword, 'VERTEX_COORDINATE' ) ) then

    else if ( s_eqi ( keyword, 'SKIP' ) ) then

    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MESH_SIZE_READ - Fatal error!'
      write ( *, '(a,i8)' ) &
        '  Could not find keyword while reading line ', line_num
      write ( *, '(a)' ) '"' // trim ( text ) // '".'
      stop

    end if

  end do
!
!  Close the file.
!
  close ( unit = fileunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) '  Read ', line_num, &
    ' lines from "' // trim ( filename ) // '".'

  return
end
subroutine r8mat_write ( output_filename, m, n, table )

!*****************************************************************************80
!
!! R8MAT_WRITE writes an R8MAT file.
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
!    31 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the output file name.
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) TABLE(M,N), the data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer j
  character ( len = * ) output_filename
  integer output_status
  integer output_unit
  character ( len = 30 ) string
  real ( kind = rk ) table(m,n)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_filename ) // '" on unit ', output_unit
    output_unit = -1
    stop 1
  end if
!
!  Create a format string.
!
!  For less precision in the output file, try:
!
!                                            '(', m, 'g', 14, '.', 6, ')'
!
  if ( 0 < m .and. 0 < n ) then

    write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 24, '.', 16, ')'
!
!  Write the data.
!
    do j = 1, n
      write ( output_unit, string ) table(1:m,j)
    end do

  end if
!
!  Close the file.
!
  close ( unit = output_unit )

  return
end
function s_begin ( s1, s2 )

!*****************************************************************************80
!
!! S_BEGIN is TRUE if one string matches the beginning of the other.
!
!  Discussion:
!
!    The strings are compared, ignoring blanks, spaces and capitalization.
!
!  Example:
!
!     S1              S2      S_BEGIN
!
!    'Bob'          'BOB'     TRUE
!    '  B  o b '    ' bo b'   TRUE
!    'Bob'          'Bobby'   TRUE
!    'Bobo'         'Bobb'    FALSE
!    ' '            'Bob'     FALSE    (Do not allow a blank to match
!                                       anything but another blank string.)
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S1, S2, the strings to be compared.
!
!    Output, logical S_BEGIN, is TRUE if the strings match up to
!    the end of the shorter string, ignoring case.
!
  implicit none

  logical ch_eqi
  integer i1
  integer i2
  logical s_begin
  character ( len = * ) s1
  integer s1_length
  character ( len = * ) s2
  integer s2_length

  s1_length = len_trim ( s1 )
  s2_length = len_trim ( s2 )
!
!  If either string is blank, then both must be blank to match.
!  Otherwise, a blank string matches anything, which is not
!  what most people want.
!
  if ( s1_length == 0 .or. s2_length == 0 ) then

    if ( s1_length == 0 .and. s2_length == 0 ) then
      s_begin = .true.
    else
      s_begin = .false.
    end if

    return

  end if

  i1 = 0
  i2 = 0
!
!  Find the next nonblank in S1.
!
  do

    do

      i1 = i1 + 1

      if ( s1_length < i1 ) then
        s_begin = .true.
        return
      end if

      if ( s1(i1:i1) /= ' ' ) then
        exit
      end if

    end do
!
!  Find the next nonblank in S2.
!
    do

      i2 = i2 + 1

      if ( s2_length < i2 ) then
        s_begin = .true.
        return
      end if

      if ( s2(i2:i2) /= ' ' ) then
        exit
      end if

    end do
!
!  If the characters match, get the next pair.
!
    if ( .not. ch_eqi ( s1(i1:i1), s2(i2:i2) ) ) then
      exit
    end if

  end do

  s_begin = .false.

  return
end
function s_eqi ( s1, s2 )

!*****************************************************************************80
!
!! S_EQI is a case insensitive comparison of two strings for equality.
!
!  Discussion:
!
!    S_EQI ( 'Anjana', 'ANJANA' ) is TRUE.
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
!    Input, character ( len = * ) S1, S2, the strings to compare.
!
!    Output, logical S_EQI, the result of the comparison.
!
  implicit none

  character c1
  character c2
  integer i
  integer lenc
  logical s_eqi
  character ( len = * ) s1
  integer s1_length
  character ( len = * ) s2
  integer s2_length

  s1_length = len ( s1 )
  s2_length = len ( s2 )
  lenc = min ( s1_length, s2_length )

  s_eqi = .false.

  do i = 1, lenc

    c1 = s1(i:i)
    c2 = s2(i:i)
    call ch_cap ( c1 )
    call ch_cap ( c2 )

    if ( c1 /= c2 ) then
      return
    end if

  end do

  do i = lenc + 1, s1_length
    if ( s1(i:i) /= ' ' ) then
      return
    end if
  end do

  do i = lenc + 1, s2_length
    if ( s2(i:i) /= ' ' ) then
      return
    end if
  end do

  s_eqi = .true.

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

