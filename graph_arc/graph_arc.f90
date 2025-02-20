subroutine digraph_arc_edge_sort ( nedge, inode, jnode )

!*****************************************************************************80
!
!! digraph_arc_edge_sort() sorts the edge array of a graph.
!
!  Discussion:
!
!    The edges are sorted in dictionary order.
!
!  Example:
!
!    Input:
!
!      INODE  JNODE
!
!        3      2
!        2      4
!        4      3
!        2      1
!        1      4
!
!    Output:
!
!      INODE  JNODE
!
!        1      4
!        2      1
!        2      4
!        3      2
!        4      3
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    24 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NEDGE, the number of edges.
!
!    Input/output, integer INODE(NEDGE), JNODE(NEDGE), the edge
!    array.  The I-th edge goes from node INODE(I) to node JNODE(I).
!    On output, the INODE and JNODE arrays have been sorted as described.
!
  implicit none

  integer nedge

  integer iedge
  integer indx
  integer inode(nedge)
  integer isgn
  integer jedge
  integer jnode(nedge)

  if ( nedge <= 1 ) then
    return
  end if
!
!  Sort the edges using an external heap sort.
!
  iedge = 0
  jedge = 0
  indx = 0
  isgn = 0

  do

    call sort_heap_external ( nedge, indx, iedge, jedge, isgn )
!
!  Interchange edges IEDGE and JEDGE.
!
    if ( 0 < indx ) then

      call i4_swap ( inode(iedge), inode(jedge) )
      call i4_swap ( jnode(iedge), jnode(jedge) )
!
!  Compare edges IEDGE and JEDGE.
!
    else if ( indx < 0 ) then

      if ( ( inode(iedge) < inode(jedge) ) .or. &
        ( inode(iedge) == inode(jedge) .and. &
          jnode(iedge) < jnode(jedge) ) ) then
        isgn = -1
      else
        isgn = +1
      end if

    else if ( indx == 0 ) then

      exit

    end if

  end do
 
  return
end
subroutine digraph_arc_print ( nedge, inode, jnode, title )

!*****************************************************************************80
!
!! digraph_arc_print() prints a digraph from an edge list.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NEDGE, the number of edges.
!
!    Input, integer INODE(NEDGE), JNODE(NEDGE), the beginning and
!    end nodes of the edges.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer nedge

  integer i
  integer inode(nedge)
  integer jnode(nedge)
  character ( len = * ) title

  if ( len_trim ( title ) /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  do i = 1, nedge
    write ( *, '(i8,4x,2i8)' ) i, inode(i), jnode(i)
  end do

  return
end
subroutine edges_to_ps ( plotxmin2, plotymin2, alpha, iunit, inode, jnode, &
  nedge, nnode, x, y, xmin, ymin )

!*****************************************************************************80
!
!! edges_to_ps() writes subplot edges to a PostScript file.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer PLOTXMIN2, PLOTYMIN2, the Postscript origin.
!
!    Input, real ( kind = rk ) ALPHA, the physical-to-Postscript scale factor.
!
!    Input, integer IUNIT, the output FORTRAN unit.
!
!    Input, integer INODE(NEDGE), JNODE(NEDGE), the edge array.
!    The I-th edge connects nodes INODE(I) and JNODE(I).
!
!    Input, integer NEDGE, the number of edges.
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, real ( kind = rk ) X(NNODE), Y(NNODE), the X and Y components
!    of points.
!
!    Input, real ( kind = rk ) XMIN, YMIN, the physical origin.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nedge
  integer nnode

  real ( kind = rk ) alpha
  integer i
  integer inode(nedge)
  integer iunit
  integer jnode(nedge)
  integer node
  integer plotxmin2
  integer plotymin2
  integer px1
  integer px2
  integer py1
  integer py2
  real ( kind = rk ) x(nnode)
  real ( kind = rk ) xmin
  real ( kind = rk ) y(nnode)
  real ( kind = rk ) ymin
!
!  Draw lines.
!
  do i = 1, nedge

    node = inode(i)
    px1 = plotxmin2 + nint ( alpha * ( x(node) - xmin ) )
    py1 = plotymin2 + nint ( alpha * ( y(node) - ymin ) )

    node = jnode(i)
    px2 = plotxmin2 + nint ( alpha * ( x(node) - xmin ) )
    py2 = plotymin2 + nint ( alpha * ( y(node) - ymin ) )

    write ( iunit, '(2i4,a,2i4,a)' ) px1, py1, ' moveto ', px2, py2, &
      ' lineto stroke'

  end do

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! get_unit() returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
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
!    Output, integer IUNIT.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5 and 6).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
  implicit none

  integer i
  integer ios
  integer iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do
!
!  No free unit was found.
!
  iunit = 0

  return
end
subroutine graph_adj_edge_count ( adj, nnode, nedge )

!*****************************************************************************80
!
!! graph_adj_edge_count() counts the edges in a graph.
!
!  Discussion:
!
!    Self-edges are allowed.
!
!    The adjacency matrix is assumed to be symmetric, so all edges
!    except self-edges will appear twice.
!
!    The resulting sum must count the self-edges double, and then divide
!    the number of edges counted by 2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer adj(nnode,nnode), the adjacency information.
!    ADJ(I,J) is 1 if there is an edge from node I to node J.
!
!    integer NNODE, the number of nodes.
!
!  Output:
!
!    integer NEDGE, the number of edges in the graph.
!
  implicit none

  integer nnode

  integer adj(nnode,nnode)
  integer i
  integer j
  integer nedge

  nedge = 0

  do i = 1, nnode
    do j = 1, nnode

      if ( i == j ) then
        nedge = nedge + 2 * adj(i,j)
      else
        nedge = nedge + adj(i,j)
      end if

    end do
  end do

  nedge = nedge / 2

  return
end
subroutine graph_adj_edge_select ( adj, nnode, ni, nj )

!*****************************************************************************80
!
!! graph_adj_edge_select() returns one edge from a graph.
!
!  Discussion:
!
!    This function returns the first edge it encounters.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    25 February 2023
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer adj(nnode,nnode), the adjacency matrix for the graph.
!
!    integer NNODE, the number of nodes.
!
!  Output:
!
!    integer ni, nj: the endpoints of an edge of the graph.
!    If no edge was found, ni and nj are -1.
!
  implicit none

  integer nnode

  integer adj(nnode,nnode)
  integer i
  integer j
  integer ni
  integer nj

  ni = -1
  nj = -1

  do i = 1, nnode
    do j = 1, nnode
      if ( adj(i,j) /= 0 ) then
        ni = i
        nj = j
        return
      end if
    end do
  end do

  return
end
subroutine graph_adj_is_edge_connected ( adj, nnode, is_connected )

!*****************************************************************************80
!
!! graph_adj_is_edge_connected() determines if a graph is edgewise connected.
!
!  Definition:
!
!    A graph is edgewise connected if from any edge it is possible to reach
!    any other edge.  An edgewise connected graph may include isolated nodes.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    25 February 2023
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer adj(nnode,nnode), the adjacency matrix for the 
!    graph.  ADJ(I,J) is nonzero if there is an edge from node I to node J.
!
!    integer NNODE, the number of nodes.
!
!  Output:
!
!    logical is_connected: true if the graph is edgewise connected.
!
  implicit none

  integer nnode

  integer adj(nnode,nnode)
  integer found(nnode)
  integer i
  integer ihi
  integer ii
  integer ilo
  logical is_connected
  integer j
  integer jhi
  integer jlo
  integer list(nnode)
  integer ni
  integer nj
!
!  Select an arbitrary edge (ni,nj) from the graph.
!
  call graph_adj_edge_select ( adj, nnode, ni, nj )

  if ( ni == -1 ) then
    is_connected = .false.
    return
  end if
!
!  FOUND(I) is 1 if edge I has been reached.
!  LIST(I) contains a list of the nodes as they are reached.
!
  list(1:nnode) = 0
  found(1:nnode) = 0

  ilo = 1
  ihi = 1
  list(1) = ni
  found(ni) = 1
  if ( ni /= nj ) then
    ihi = 2
    list(2) = nj
    found(nj) = 1
  end if

  adj(ni,nj) = - adj(ni,nj)
  adj(nj,ni) = - adj(nj,ni)
!
!  From the batch of edge nodes found last time, LIST(ILO:IHI),
!  look for unfound neighbors, and store their indices in LIST(JLO:JHI).
!
  do while ( .true. )

    jlo = ihi + 1
    jhi = ihi

    do ii = ilo, ihi

      i = list(ii)

      do j = 1, nnode

        if ( 0 < adj(i,j) ) then

          adj(i,j) = - adj(i,j)
          if ( 0 < adj(j,i) ) then
            adj(j,i) = - adj(j,i)
          end if

          if ( found(j) == 0 ) then
            jhi = jhi + 1
            list(jhi) = j
            found(j) = 1
          end if

        end if

      end do

    end do

    if ( jhi < jlo ) then
      exit
    end if

    ilo = jlo
    ihi = jhi

  end do
!
!  If any edges were unvisited, then the graph is not edgewise connected.
!
  is_connected = .true.

  do i = 1, nnode
    do j = 1, nnode
      if ( 0 < adj(i,j) ) then
        is_connected = .false.
      end if
    end do
  end do
!
!  Restore the positive sign of ADJ.
!
  adj(1:nnode,1:nnode) = abs ( adj(1:nnode,1:nnode) )

  return
end
subroutine graph_adj_is_node_connected ( adj, nnode, result )

!*****************************************************************************80
!
!! graph_adj_is_node_connected() determines if a graph is nodewise connected.
!
!  Definition:
!
!    A graph is nodewise connected if, from every node, there is a path
!    to any other node.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 March 2022
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer adj(nnode,nnode), the adjacency matrix for the 
!    graph.  ADJ(I,J) is nonzero if there is an edge from node I to node J.
!
!    integer NNODE, the number of nodes.
!
!  Output:
!
!    integer RESULT.
!    0, the graph is not nodewise connected.
!    1, the graph is nodewise connected.
!
  implicit none

  integer nnode

  integer adj(nnode,nnode)
  integer found(nnode)
  integer i
  integer ihi
  integer ii
  integer ilo
  integer j
  integer jhi
  integer jlo
  integer list(nnode)
  integer result
!
!  FOUND(I) is 1 if node I has been reached.
!  LIST(I) contains a list of the nodes as they are reached.
!
  list(1:nnode) = 0
  found(1:nnode) = 0
!
!  Start at node 1.
!
  found(1) = 1
  list(1) = 1
  ilo = 1
  ihi = 1
!
!  From the batch of nodes found last time, LIST(ILO:IHI),
!  look for unfound neighbors, and store their indices in LIST(JLO:JHI).
!
  do

    jlo = ihi + 1
    jhi = ihi

    do ii = ilo, ihi

      i = list(ii)

      do j = 1, nnode

        if ( adj(i,j) /= 0 .or. adj(j,i) /= 0 ) then

          if ( found(j) == 0 ) then
          jhi = jhi + 1
            list(jhi) = j
            found(j) = 1
          end if

        end if

      end do

    end do
!
!  If any neighbors were found, go back and find THEIR neighbors.
!
    if ( jhi < jlo ) then
      exit
    end if

    ilo = jlo
    ihi = jhi

  end do
!
!  No more neighbors were found.  Have we reached all nodes?
!
  if ( ihi == nnode ) then
    result = 1
  else
    result = 0
  end if

  return
end
subroutine graph_adj_is_tree ( adj, nnode, result )

!*****************************************************************************80
!
!! graph_adj_is_tree() determines whether a graph is a tree.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer adj(nnode,nnode), the adjacency matrix for the 
!    graph.  ADJ(I,J) is nonzero if there is an edge from node I to node J.
!
!    Input, integer NNODE, the number of nodes.
!
!    Output, integer RESULT.
!    0, the graph is not a tree.
!    1, the graph is a tree.
!
  implicit none

  integer nnode

  integer adj(nnode,nnode)
  integer nedge
  integer result

  if ( nnode <= 1 ) then
    result = 1
    return
  end if
!
!  Every node must be connected to every other node.
!
  call graph_adj_is_node_connected ( adj, nnode, result )

  if ( result == 0 ) then
    return
  end if
!
!  There must be exactly NNODE-1 edges.
!
  call graph_adj_edge_count ( adj, nnode, nedge )

  if ( nedge == nnode - 1 ) then
    result = 1
  else
    result = 0
  end if

  return
end
subroutine graph_adj_print ( adj, nnode, title )

!*****************************************************************************80
!
!! graph_adj_print() prints out an adjacency matrix for a graph.
!
!  Discussion:
!
!    This routine actually allows the entries of ADJ to have ANY value.
!    Values between 0 and 9 will be printed as is.  Other values will
!    be printed as '*'.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ADJ(NNODEA,NNODE), the adjacency matrix of a 
!    graph.  ADJ(I,J) is 1 if there is a direct connection FROM node I TO 
!    node J, and is 0 otherwise.
!
!    Input, integer NNODE, the number of nodes.  
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer nnode

  integer adj(nnode,nnode)
  integer i
  integer j
  integer jhi
  character ( len = 80 ) string
  character ( len = * ) title

  if ( len_trim ( title ) /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  do i = 1, nnode

    jhi = min ( nnode, 74 )

    do j = 1, jhi

      if ( 0 <= adj(i,j) .and. adj(i,j) <= 9 ) then
        string(j:j) = char ( 48 + adj(i,j) )
      else
        string(j:j) = '*'
      end if

    end do

    write ( *, '(2x,i3,2x,a)' ) i, string(1:jhi)

  end do

  return
end
subroutine graph_arc_chromatic ( nnode, nedge, inode, jnode, iarray, jarray, &
  karray, stack, maxstack )

!*****************************************************************************80
!
!! graph_arc_chromatic() calculates the chromatic polynomial of a connected graph.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 March 2023
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, integer NEDGE, the number of edges.
!
!    Input, integer INODE(NEDGE), JNODE(NEDGE), the edge endpoints.
!
!    Output, integer IARRAY(NNODE).  Coefficients of the chromatic
!    polynomial in power form:
!
!      P(X) = 
!        IARRAY(N)   * X^NNODE
!      - IARRAY(N-1) * X^NNODE-1
!      + IARRAY(N-2) * X^NNODE-2
!      ...
!      +-IARRAY(1)   * X
!
!    Output, integer JARRAY(NNODE).  Coefficients of the chromatic
!    polynomial using the Tutte or tree form:
!
!      P(X) = SUM ( I = 1 TO NNODE ) 
!        (-1)^(NNODE-I) * IARRAY(I) * X * (X-1)^(I-1)
!
!    Output, integer KARRAY(NNODE).  The Stirling or factorial form
!    of the chromatic polynomial.
!
!      P(X) = SUM ( I = 1 TO NNODE ) KARRAY(I) * (X)(I)
!
!    Here (X)(I) is meant to represent X*(X-1)*(X-2)...*(X-I+1).
!
!    Workspace, integer STACK(2,MAXSTACK).
!
!    Input, integer MAXSTACK, dimension of working storage.  An 
!    upper estimate for the amount of storage required is
!    NNODE * ( IE - 0.5*(NNODE-1)).
!
  implicit none

  integer nedge
  integer nnode
  integer maxstack

  integer i
  integer iarray(nnode)
  integer nedge1
  integer ien(2)
  integer iendpt(2,nedge)
  integer inode(nedge)
  integer is
  integer iu
  integer iv
  integer j
  integer jarray(nnode)
  integer jhi
  integer jnode(nedge)
  integer k
  integer karray(nnode)
  integer l
  integer nnode1
  integer stack(2,maxstack)

  iendpt(1,1:nedge) = inode(1:nedge)
  iendpt(2,1:nedge) = jnode(1:nedge)

  is = 0
  jarray(1:nnode) = 0
  nedge1 = nedge
  nnode1 = nnode
 
10    continue
 
  if ( 0 < nnode1 .and. 0 < nedge1 ) then
    call span_forest ( nnode1, nedge1, iendpt, k, karray )
  else
    k = 0
  end if
 
  if ( k /= 1 ) then
    go to 50
  end if

  if ( nedge1 < nnode1 ) then
    go to 40
  end if
 
  if ( nedge1 == nnode1 ) then
 
    jarray(nnode1) = jarray(nnode1) + 1
 
  else
 
    do i = 1, nedge1
      is = is + 1
      stack(1,is) = iendpt(1,i)
      stack(2,is) = iendpt(2,i)
    end do
 
    stack(1,is) = nnode1
    stack(2,is) = nedge1 - 1
 
  end if
 
20    continue
 
  iarray(1:nnode) = 0
  iu = min ( iendpt(1,nedge1), iendpt(2,nedge1) )
  iv = iendpt(1,nedge1) + iendpt(2,nedge1) - iu
  jhi = nedge1 - 1
  nedge1 = 0
 
  do j = 1, jhi
 
    do l = 1, 2

      ien(l) = iendpt(l,j)

      if ( ien(l) == iv ) then
        ien(l) = iu
      end if

      if ( ien(l) == nnode1 ) then
        ien(l) = iv
      end if

    end do
 
    do l = 1, 2
 
      if ( ien(l) == iu ) then
        if ( iarray(ien(3-l)) /= 0 ) then
          go to 30
        end if
        iarray(ien(3-l)) = 1
      end if
 
    end do
 
    nedge1 = nedge1 + 1
 
    iendpt(1,nedge1) = ien(1)
    iendpt(2,nedge1) = ien(2)
 
30      continue
 
  end do
 
  nnode1 = nnode1 - 1
  go to 10
 
40    continue
 
  jarray(nnode1) = jarray(nnode1) + 1
 
  if ( is /= 0 ) then
 
    nnode1 = stack(1,is)
    nedge1 = stack(2,is)
    is = is - nedge1 - 1
 
    do i = 1, nedge1
      iendpt(1,i) = stack(1,is+i)
      iendpt(2,i) = stack(2,is+i)
    end do
 
    if ( nedge1 == nnode1 ) then
      jarray(nnode1) = jarray(nnode1) + 1
    else
      is = is + nedge1
      stack(1,is) = nnode1
      stack(2,is) = nedge1 - 1
    end if
 
    go to 20
 
  end if
 
50    continue
 
  do i = 1, nnode
    iarray(i) = jarray(i)
    karray(i) = ( 1 - 2 * mod ( nnode-i, 2 ) ) * jarray(i)
  end do
 
  call poly ( nnode, iarray, 1, nnode, iv )
  call poly ( nnode, karray, 0, -2, iv )
 
  return
end
subroutine graph_arc_complement ( inode, jnode, inode2, jnode2, maxedge, &
  nedge, nedge2, nnode )

!*****************************************************************************80
!
!! graph_arc_complement() returns the edge list of the complement of a graph.
!
!  Discussion:
!
!    This routine can also handle a directed graph.
!
!  Definition:
!
!    The complement of a graph G is a graph H with the property that
!    nodes U and V are connected in H if and only if they are not
!    connected in G.  However, edges from a node to itself are not
!    allowed.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    24 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer INODE(NEDGE), JNODE(NEDGE).  INODE(I) 
!    and JNODE(I) are the start and end nodes of the I-th edge of the graph G.
!    On output, the data in INODE and JNODE will have been sorted, but not
!    otherwise disrupted.
!
!    Output, integer INODE2(MAXEDGE), JNODE2(MAXEDGE).  INODE2(I) 
!    and JNODE2(I) are the start and end nodes of the I-th edge of the 
!    complement graph H.
!
!    Input, integer MAXEDGE, the amount of storage available in 
!    INODE2 and JNODE2.  MAXEDGE only needs to be as large as NEDGE2, and NEDGE2
!    can be precomputed, assuming that the input value of NEDGE does not
!    count any self edges (edges from a node to itself), and does not
!    count an edge twice (that is, counting the edge from I to J, and
!    the edge from J to I, as distinct).  If this is so, then you can
!    set MAXEDGE = NEDGE2 = 0.5 * ( NNODE * ( NNODE - 1 ) ) - NEDGE.
!
!    Input, integer NEDGE, the number of edges in the graph G.
!
!    Output, integer NEDGE2, the number of edges in the complement 
!    graph H.
!
!    Input, integer NNODE, the number of nodes.
!
  implicit none

  integer maxedge
  integer nedge

  integer i
  integer i1
  integer i2
  integer inedge
  integer inode(nedge)
  integer inode2(maxedge)
  integer j
  integer j1
  integer j2
  integer jnode(nedge)
  integer jnode2(maxedge)
  integer nedge2
  integer nnode
!
!  Sort the input edge array.
!
  call graph_arc_edge_sort ( nedge, inode, jnode )
!
!  Compute the complementary edges.
!
  nedge2 = 0
 
  inedge = 0
  i2 = 1
  j2 = 1

  do while ( inedge < nedge )
 
    inedge = inedge + 1
    i1 = i2
    j1 = j2

    if ( inedge <= nedge ) then
      i2 = inode(inedge)
      j2 = jnode(inedge)
    else
      i2 = nnode
      j2 = nnode
    end if
 
    if ( i1 == i2 ) then
 
      do j = j1+1, j2-1
        if ( i1 < j ) then
          nedge2 = nedge2 + 1
          inode2(nedge2) = i2
          jnode2(nedge2) = j
        end if
      end do
 
    else
 
      do j = j1+1, nnode
        if ( i1 < j ) then
          nedge2 = nedge2 + 1
          inode2(nedge2) = i1
          jnode2(nedge2) = j
        end if
      end do
 
      do i = i1+1, i2-1
        do j = 1, nnode
          if ( i < j ) then
            nedge2 = nedge2 + 1
            inode2(nedge2) = i
            jnode2(nedge2) = j
          end if
        end do
      end do

      do j = 1, j2-1
        if ( i2 < j ) then
          nedge2 = nedge2 + 1
          inode2(nedge2) = i2
          jnode2(nedge2) = j
        end if
      end do
 
    end if

  end do
 
  return
end
subroutine graph_arc_degree ( nnode, nedge, inode, jnode, degree )

!*****************************************************************************80
!
!! graph_arc_degree() determines the degree of the nodes of a graph.
!
!  Definition:
!
!    The degree of a node is the number of edges that include the node.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    24 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer NNODE, the number of nodes.
!
!    integer NEDGE, the number of edges.
!
!    integer INODE(NEDGE), JNODE(NEDGE), the pairs of nodes
!    that form the edges.
!
!  Output:
!
!    integer DEGREE(NNODE), the degree of each node, 
!    that is, the number of edges that include the node.
!
  implicit none

  integer nedge
  integer nnode

  integer degree(nnode)
  integer i
  integer inode(nedge)
  integer jnode(nedge)
  integer n

  degree(1:nnode) = 0

  do i = 1, nedge

    n = inode(i)
    if ( 1 <= n .and. n <= nnode ) then
      degree(n) = degree(n) + 1
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GRAPH_ARC_DEGREE(): Fatal error!'
      write ( *, '(a,i8)' ) '  Out-of-range node value = ', n
      stop 1
    end if

    n = jnode(i)
    if ( 1 <= n .and. n <= nnode ) then
      degree(n) = degree(n) + 1
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GRAPH_ARC_DEGREE(): Fatal error!'
      write ( *, '(a,i8)' ) '  Out-of-range node value = ', n
      stop 1
    end if

  end do

  return
end
subroutine graph_arc_edge_con2 ( nnode, nedge, inode, jnode, edge_con )

!*****************************************************************************80
!
!! graph_arc_edge_con2() finds the edge-connectivity of a connected graph.
!
!  Method:
!
!    A graph G has edge connectivity K if, given any pair of distinct nodes 
!    I and J, there are K paths from I to J, no two of which use a common edge.
!
!    Thus, in particular, if a graph G is Hamiltonian, it must have 
!    edge connectivity at least 2.  For we can simply take the Hamiltonian
!    circuit, and use the part from I to J as the first path, and the
!    part from J to I as the second, simply reversing the direction
!    of traversal.
!
!    To determine the edge connectivity, for each J from 2 to NNODE do 
!    the following:
!
!      Take node 1 as the source, node J as the sink in G, assign a unit 
!      capacity to all edges in both directions, and find the value of the
!      maximum flow G(J) in the resulting network.  
!
!    The edge-connectivity is then equal to the minimum of G(2:NNODE).
!
!    This routine finds the edge-connectivity of a given undirected graph with 
!    the help of a maximum flow algorithm.
!
!    The maximum network flow algorithm requires O(NNODE**3) operations.  The 
!    edge-connectivity of a graph will therefore be found in O(NNODE**4) 
!    operations.
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
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Input:
!
!    integer NNODE, the number of nodes.
!
!    integer NEDGE, the number of edges.
!
!    integer INODE(NEDGE), JNODE(NEDGE), the end nodes of 
!    the edges.
!
!  Output:
!
!    integer EDGE_CON, the edge-connectivity of the graph.
!
  implicit none

  integer nedge
  integer nnode

  integer capflo(2,2*nedge)
  integer edge_con
  integer i
  integer icut(nnode)
  integer iendpt(2,2*nedge)
  integer inode(nedge)
  integer isink
  integer isorce
  integer j
  integer jnode(nedge)
  integer node_flow(nnode)
!
!  Create the network from the graph.
!
  j = 0
 
  do i = 1, nedge

    j = j + 1
    iendpt(1,j) = inode(i)
    iendpt(2,j) = jnode(i)
    capflo(1,j) = 1
    capflo(2,j) = 0

    j = j + 1
    iendpt(1,j) = jnode(i)
    iendpt(2,j) = inode(i)
    capflo(1,j) = 1
    capflo(2,j) = 0

  end do
!
!  Call the network flow algorithm.
!
  edge_con = nnode
  isorce = 1

  do isink = 2, nnode

    call network_flow_max ( nnode, 2*nedge, iendpt, capflo, isorce, isink, &
      icut, node_flow )
 
    if ( node_flow(isorce) < edge_con ) then
      edge_con = node_flow(isorce)
    end if

  end do
 
  return
end
subroutine graph_arc_edge_sort ( nedge, inode, jnode )

!*****************************************************************************80
!
!! graph_arc_edge_sort() sorts the edge array of a graph.
!
!  Discussion:
!
!    The pair of nodes (I,J) representing an edge is reordered so
!    that the smaller node is listed first.
!
!    Then the edges are sorted in dictionary order.
!
!  Example:
!
!    Input:
!
!      INODE  JNODE
!
!        3      2
!        4      3
!        2      1
!        1      4
!
!    Output:
!
!      INODE  JNODE
!
!        1      2
!        1      4
!        2      3
!        3      4
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    24 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NEDGE, the number of edges.
!
!    Input/output, integer INODE(NEDGE), JNODE(NEDGE), the edge 
!    array of a graph.  The I-th edge of the graph connects nodes INODE(I) and
!    JNODE(I).
!
!    On output, the INODE and JNODE arrays have been sorted as described.
!
  implicit none

  integer nedge

  integer i
  integer iedge
  integer indx
  integer inode(nedge)
  integer isgn
  integer jedge
  integer jnode(nedge)

  if ( nedge <= 1 ) then
    return
  end if
!
!  Sort the node pairs.
!
  do i = 1, nedge
    if ( jnode(i) < inode(i) ) then
      call i4_swap ( inode(i), jnode(i) )
    end if
  end do
!
!  Sort the edges using an external heap sort.
!
  iedge = 0
  jedge = 0
  indx = 0
  isgn = 0

  do

    call sort_heap_external ( nedge, indx, iedge, jedge, isgn )
!
!  Interchange edges IEDGE and JEDGE.
!
    if ( 0 < indx ) then

      call i4_swap ( inode(iedge), inode(jedge) )
      call i4_swap ( jnode(iedge), jnode(jedge) )
!
!  Compare edges IEDGE and JEDGE.
!
    else if ( indx < 0 ) then

      if ( ( inode(iedge) < inode(jedge) ) .or. &
         ( inode(iedge) == inode(jedge) .and. &
           jnode(iedge) < jnode(jedge) ) ) then
        isgn = -1
      else
        isgn = +1
      end if

    else

      exit

    end if

  end do
 
  return
end
subroutine graph_arc_euler_circ ( nnode, nedge, inode, jnode, circuit )

!*****************************************************************************80
!
!! graph_arc_euler_circ() finds an Euler circuit in an Eulerian graph.
!
!  Discussion:
!
!    An Euler circuit of a graph is a path that uses each edge exactly once.
!
!    A graph is Eulerian if it has an Euler circuit.
!
!    An Eulerian graph may have many circuits; this routine only finds one.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 July 2000
!
!  Reference:
!
!    Hang Tong Lau,
!    Combinatorial Heuristic Algorithms in FORTRAN,
!    Springer Verlag, 1986.
!
!  Input:
!
!    integer NNODE, the number of nodes in the graph.
!
!    integer NEDGE, the number of edges in the graph.
!
!    integer INODE(NEDGE), JNODE(NEDGE), the two end nodes 
!    of each edge.
!
!  Output:
!
!    integer CIRCUIT(NEDGE), the Euler circuit, as a 
!    series of nodes.
!
  implicit none

  integer nedge
  integer nnode

  integer circuit(nedge)
  logical copyon
  logical found
  integer i
  integer ibase
  integer iforwd
  integer inode(nedge)
  integer insert
  integer ipivot
  integer iwork1(nedge)
  integer iwork2(nedge)
  integer iwork3(nnode)
  integer iwork4(nnode)
  integer iwork5(nnode)
  integer iwork6(nnode)
  integer j
  integer jnode(nedge)
  integer k
  integer l
  integer locbas
  integer nbreak
  integer ncopy
  integer numarc
  integer numnode
!
!  The number of times each node has been visited begins at 0.
!
  iwork3(1:nnode) = 0
  circuit(1:nedge) = 0
  iwork1(1:nedge) = 0
  iwork2(1:nedge) = 0
!
!  Begin the Euler circuit with the edge INODE(1), JNODE(1).
!
  numarc = 1
  iwork2(1) = 1

  numnode = 1
  i = inode(1)
  iwork1(numnode) = i
  iwork3(i) = 1

  numnode = numnode + 1
  j = jnode(1)
  iwork1(numnode) = j
  iwork3(j) = 1

  ibase = j
  nbreak = 0
!
!  Look for the next arc.
!
30    continue

  do i = 2, nedge

    if ( iwork2(i) == 0 ) then

      if ( inode(i) == ibase ) then

        found = .true.
        ibase = jnode(i)

      else if ( jnode(i) == ibase ) then

        found = .true.
        ibase = inode(i)

      else

        found = .false.

      end if

      if ( found ) then
        iwork2(i) = 1
        numarc = numarc + 1
        numnode = numnode + 1
        if ( numnode <= nedge ) then
          iwork1(numnode) = ibase
        end if
        iwork3(ibase) = 1
        go to 30
      end if

    end if

  end do
!
!  A cycle has been found.
!
  if ( 0 < nbreak ) then
    numnode = numnode - 1
    iwork5(nbreak) = numnode
  end if

  if ( numarc < nedge ) then

    iwork1(numnode) = ibase
!
!  Find a node in the current Euler circuit.
!
    do i = 2, nedge

      if ( iwork2(i) == 0 ) then

        if ( iwork3(inode(i)) /= 0 ) then

          found = .true.
          j = inode(i)
          k = jnode(i)

        else if ( iwork3(jnode(i)) /= 0 ) then

          found = .true.
          j = jnode(i)
          k = inode(i)

        else

          found = .false.

        end if
!
!  Identify a path which will be added to the circuit.
!
        if ( found ) then
          nbreak = nbreak + 1
          iwork6(nbreak) = j
          ibase = k
          iwork3(k) = 1
          numnode = numnode + 1
          iwork4(nbreak) = numnode
          iwork1(numnode) = ibase
          iwork2(i) = 1
          numarc = numarc + 1
          go to 30
        end if

      end if

    end do

  end if
!
!  Form the Euler circuit.
!
  if ( nbreak == 0 ) then
    numnode = numnode - 1
    circuit(1:numnode) = iwork1(1:numnode)
    return
  end if

  insert = 1
  ipivot = iwork6(insert)
  iforwd = 0

  do

    ncopy = 1
    ibase = iwork1(1)
    locbas = 1
    circuit(ncopy) = ibase
!
!  A path identified before is added to the circuit.
!
80  continue

    if ( ibase == ipivot ) then

      j = iwork4(insert) + iforwd
      k = iwork5(insert) + iforwd

      do l = j, k
        ncopy = ncopy + 1
        circuit(ncopy) = iwork1(l)
        iwork1(l) = 0
      end do

      ncopy = ncopy + 1
!
!  Add the intersecting node to the circuit.
!
      circuit(ncopy) = ibase
      iforwd = iforwd + 1

      if ( ncopy < numnode ) then

        do

          if ( nedge <= ncopy ) then
            exit
          end if

          locbas = locbas + 1

          if ( nedge <= locbas ) then
            exit
          end if

          ibase = iwork1(locbas)

          if ( ibase /= 0 ) then
            ncopy = ncopy + 1
            circuit(ncopy) = ibase
          end if

        end do

      end if

    else

      ncopy = ncopy + 1

       if ( ncopy <= numnode ) then
         locbas = locbas + 1
         ibase = iwork1(locbas)
         circuit(ncopy) = ibase
         go to 80
       end if

     end if
!
!  Check if more paths are to be added to the circuit.
!
    copyon = .false.
    insert = insert + 1

    if ( insert <= nbreak ) then
      copyon = .true.
      ipivot = iwork6(insert)
    end if

    if ( .not. copyon ) then
      exit
    end if

    iwork1(1:nedge) = circuit(1:nedge)

  end do

  return
end
subroutine graph_arc_euler_circ_cand ( nedge, inode, jnode, circuit, k, &
  nstack, stack, maxstack, ncan, iwork )

!*****************************************************************************80
!
!! graph_arc_euler_circ_cand(): candidates for the K-th edge of an Euler circuit.
!
!  Discussion:
!
!    This routine is used in conjunction with I4VEC_BACKTRACK, which directs 
!    the search for a complete Euler circuit.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 July 2000
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer NEDGE, the number of edges in the graph.
!
!    Input, integer INODE(NEDGE), JNODE(NEDGE), the edge array of 
!    the graph.  The I-th edge extends from node INODE(I) to JNODE(I).
!
!    Input, integer CIRCUIT(NEDGE), CIRCUIT(I) is the I-th edge 
!    in the circuit.  A full circuit will have NEDGE edges, but on input we 
!    only have K-1.
!
!    Input, integer K, the index of the next edge to be determined
!    in circuit.
!
!    Input/output, integer NSTACK, the current length of the stack.
!
!    Input, integer STACK(MAXSTACK).  As yet unused candidates for 
!    positions 1 to K-1.
!
!    Input, integer MAXSTACK, the dimension of STACK.
!
!    Input/output, integer NCAN(NEDGE), the number of candidates 
!    for each position.
!
!    Workspace, integer IWORK(NEDGE).
!
  implicit none

  integer nedge
  integer maxstack

  integer circuit(nedge)
  integer i
  integer inode(nedge)
  integer it
  integer iwork(nedge)
  integer jnode(nedge)
  integer k
  logical lwork(nedge)
  integer ncan(nedge)
  integer nstack
  integer stack(maxstack)

  ncan(k) = 0

  if ( k == 1 ) then
    iwork(1) = jnode(1)
    stack(1) = 1
    nstack = 1
    ncan(k) = 1
    return
  end if
 
  if ( 2 < k ) then
    iwork(k-1) = inode(circuit(k-1)) + jnode(circuit(k-1)) - iwork(k-2)
  end if
 
  it = iwork(k-1)
 
  do i = 1, nedge
    lwork(i) = it == inode(i) .or. it == jnode(i)
  end do
 
  lwork(circuit(1:k-1)) = .false.
 
  do i = 1, nedge
    if ( lwork(i) ) then
      if ( maxstack <= nstack ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRAPH_ARC_EULER_CIRC_CAND(): Fatal error!'
        write ( *, '(a)' ) '  Stacksize exceeded!'
        stop 1
      end if
      nstack = nstack + 1
      stack(nstack) = i
      ncan(k) = ncan(k) + 1
    end if
  end do
 
  return
end
subroutine graph_arc_euler_circ_next ( nedge, inode, jnode, circuit, stack, &
  maxstack, ncan, more )

!*****************************************************************************80
!
!! graph_arc_euler_circ_next() returns the next Euler circuit for a graph.
!
!  Discussion:
!
!    The routine produces all the Euler circuits of a graph, one at a time.
!
!    An Euler circuit of a graph is a path starting at some node, 
!    using all the edges of the graph exactly once, and returning
!    to the starting node.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 August 2000
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer NEDGE, the number of edges in the graph.
!
!    Input, integer INODE(NEDGE), JNODE(NEDGE), the edge array 
!    of the graph.  The I-th edge extends from node INODE(I) to JNODE(I).
!
!    Output, integer CIRCUIT(NEDGE).  If MORE = TRUE on output, 
!    then CIRCUIT contains the edges, in order, that constitute this circuit.
!
!    Workspace, integer STACK(MAXSTACK).  
!
!    Input, integer MAXSTACK, the dimension of STACK.
!
!    Workspace, integer NCAN(NEDGE), the number of candidates for each position.
!
!    Input/output, logical MORE.
!    On first call, set MORE to .FALSE, and do not alter it after.
!    On return, MORE is TRUE if another circuit has been returned in
!    IARRAY, and FALSE if there are no more circuits.
!
  implicit none

  integer nedge
  integer maxstack

  integer circuit(nedge)
  integer inode(nedge)
  integer, save :: indx = 0
  integer iwork(nedge)
  integer jnode(nedge)
  integer, save :: k = 0
  logical more
  integer ncan(nedge)
  integer, save :: nstack = 0
  integer stack(maxstack)

  if ( .not. more ) then
    indx = 0
    k = 0
    more = .true.
    nstack = 0
  end if
 
  do
 
    call i4vec_backtrack ( nedge, circuit, indx, k, nstack, stack, maxstack, &
      ncan )
 
    if ( indx == 1 ) then

      exit

    else if ( indx == 2 ) then

      call graph_arc_euler_circ_cand ( nedge, inode, jnode, circuit, k, &
        nstack, stack, maxstack, ncan, iwork )

    else

      more = .false.
      exit

    end if

  end do
 
  return
end
subroutine graph_arc_example_diamond ( inode, jnode, maxedge, nedge, &
  x, y, z )

!*****************************************************************************80
!
!! graph_arc_example_diamond() returns the graph of a "diamond" 3D shape.
!
!  Example:
!
!        1
!      /| |\
!     / | | \
!    2--3-4--5--(2)
!    |  | |  |
!    6--7-8--9--(6)
!     \ | | /
!      \| |/
!       10
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Output, integer INODE(MAXEDGE), JNODE(MAXEDGE), the NEDGE 
!    edges of the graph.  The I-th edge connects nodes INODE(I) and
!    JNODE(I).
!
!    Input, integer MAXEDGE, the maximum number of edges allocated
!    in the EDGE array.  MAXEDGE should be at least 20.
!
!    Output, integer NEDGE, the number of edges, which is 20.
!
!    Output, integer NNODE, the number of nodes, which is 10.
!
!    Output, real ( kind = rk ) X(10), Y(10), Z(10), the 
!    locations for the nodes.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxedge
  integer, parameter :: nnode = 10

  integer inode(maxedge)
  integer jnode(maxedge)
  integer nedge
  real ( kind = rk ) x(nnode)
  real ( kind = rk ) y(nnode)
  real ( kind = rk ) z(nnode)

  nedge = 20

  if ( maxedge < nedge ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRAPH_ARC_EXAMPLE_DIAMOND(): Fatal error!'
    write ( *, '(a,i8)' ) '  Increase MAXEDGE to at least ', nedge
    stop 1
  end if

  inode(1) = 1
  jnode(1) = 2

  inode(2) = 1
  jnode(2) = 3

  inode(3) = 1
  jnode(3) = 4

  inode(4) = 1
  jnode(4) = 5

  inode(5) = 2
  jnode(5) = 6

  inode(6) = 3
  jnode(6) = 7

  inode(7) = 4
  jnode(7) = 8

  inode(8) = 5
  jnode(8) = 9

  inode(9) = 6
  jnode(9) = 10

  inode(10) = 7
  jnode(10) = 10

  inode(11) = 8
  jnode(11) = 10

  inode(12) = 9
  jnode(12) = 10

  inode(13) = 2
  jnode(13) = 3

  inode(14) = 3
  jnode(14) = 4

  inode(15) = 4
  jnode(15) = 5

  inode(16) = 5
  jnode(16) = 2

  inode(17) = 6
  jnode(17) = 7

  inode(18) = 7
  jnode(18) = 8

  inode(19) = 8
  jnode(19) = 9

  inode(20) = 9
  jnode(20) = 6

  x(1) =  0.0D+00
  y(1) =  0.0D+00
  z(1) =  2.0D+00

  x(2) =  0.5D+00
  y(2) = -0.5D+00
  z(2) =  1.0D+00

  x(3) =  0.5D+00
  y(3) =  0.5D+00
  z(3) =  1.0D+00

  x(4) = -0.5D+00
  y(4) =  0.5D+00
  z(4) =  1.0D+00

  x(5) = -0.5D+00
  y(5) = -0.5D+00
  z(5) =  1.0D+00

  x(6) =  0.5D+00
  y(6) = -0.5D+00
  z(6) = -1.0D+00

  x(7) =  0.5D+00
  y(7) =  0.5D+00
  z(7) = -1.0D+00

  x(8) = -0.5D+00
  y(8) =  0.5D+00
  z(8) = -1.0D+00

  x(9) = -0.5D+00
  y(9) = -0.5D+00
  z(9) = -1.0D+00

  x(10) =  0.0D+00
  y(10) =  0.0D+00
  z(10) = -2.0D+00

  return
end
subroutine graph_arc_face ( face, face_count, face_order, iface, jface, &
  inode, jnode, maxface, maxorder, nedge, nface, nnode )

!*****************************************************************************80
!
!! graph_arc_face() constructs a set of faces for a plane graph.
!
!  Discussion:
!
!    Warning: This is an experimental code.
!
!    The reason this routine was written was to handle the problem of
!    converting certain forms of 3D graphics data from a point and line
!    representation to a list of faces.  While at first glance, this
!    seemed an easy task, it turned out to be one of those problems
!    that becomes harder the longer it is considered.  Particularly
!    vexing was the idea that it might be possible to do this reconstruction
!    without using any of the geometric data supplied with the connectivity
!    data.
!
!    The guiding idea was that a face ought to be a "short" cycle of
!    the graph, and that every edge ought to appear in two separate faces.
!    The resulting method should work for a connected graph which is
!    planar, or merely orientable.  A planar graph will result from a
!    reasonable "triangulation" (meaning decomposition into arbitrary
!    polygons) of the surface of a 3D object that has no holes.
!
!    This algorithm will also handle the case where the graph is not planar,
!    but results from the triangulation of a more complicated 3D object,
!    such as one that includes holes.  Even a Klein bottle, which is
!    a manifold, but not orientable, can be handled, although it may not
!    be possible then to assign a consistent orientation to the faces.
!
!    By the way, this problem is MUCH easier if we can assume that all
!    the faces use the same number of edges, such as in a triangular
!    decomposition.  This algorithm makes no such assumption.
!
!    If the graph is planar, then the decomposition into
!    faces allows us to define the dual graph.  The dual graph H of the
!    planar graph G comprises:
!
!    * nodes V(I), each of which corresponds to a face F(I) of G;
!
!    * edges (V(I), V(J)).  V(I) and V(J) share an edge in H if and only
!      if the faces F(I) and F(J) share an edge in G.  (Thus G and H
!      have the same number of edges).
!
!    In the terminology of this routine, the dual graph has NFACE nodes,
!    NEDGE edges, and the corresponding edge arrays are simply IFACE and
!    JFACE.
!
!  Formula:
!
!    If the graph is actually planar, we can regard it as the flattened
!    triangulation of a connected solid shape, and so we can apply Euler's
!    formula:
!
!      Faces + Vertices = Edges + 2
!
!    This means that we can predict beforehand that the number of faces
!    produced by this routine will be
!
!      NFACE = NEDGE + 2 - NNODE.
!
!  Notes:
!
!    The faces produced by this routine may actually overlap.  Without
!    geometric data, this is surely a possibility, since a graph may
!    have more than one embedding.  For instance, consider the following
!    two embeddings of the same graph:
!
!      A-----B       A-----B
!      |     |       |     |
!      |  E  |       D-----C
!      | / \ |        \   /
!      |/   \|         \ /
!      D-----C          E
!
!    This routine will report the two faces (A,B,C,D) and (C,D,E),
!    although in the first embedding one face seems to be part of
!    another.  This is not as bad as it might seem, since the routine
!    prefers the "smaller" face (A,B,C,D) over (A,B,C,E,D).
!
!
!    A second problem is best illustrated with a simple example.
!    Suppose we have a thin triangular rod, and that we have triangulated
!    the surface of this rod, so that the cross section of the rod
!    is a triangular graph, and the sides are made up of, say, squares.
!    Then this routine will report all the "internal" triangles as
!    faces.  It will still find the "true" faces on the sides, but
!    since it is possible to go around the diameter of the object
!    in a very few steps, the algorithm produces faces we would not
!    expect.
!
!  Restrictions:
!
!    The algorithm will fail if the graph cannot be regarded, at least
!    locally, as the triangulation of a smooth surface.  Smoothness
!    problems will not occur if the graph is planar, or results from
!    the triangulation of a 3D object, which might include holes.
!
!    The graph should be connected.
!
!    There should be no nodes of degree 1.
!
!  Method:
!
!    We have no geometric data from which to deduce physical positions
!    of the nodes.  We are only given that the graph is planar, so that
!    there is at least one embedding of the graph in the plane.
!
!    Our data structure for the method will use arrays called IFACE and JFACE.
!    For each edge I, IFACE(I) and JFACE(I) will eventually hold the
!    indices of the two faces that the edge is part of.  We begin
!    the algorithm by setting all entries of IFACE and JFACE to 0.
!
!    The second step is to find one cycle in the graph, of the shortest
!    length possible.  This cycle constitutes our first face.  We update
!    the appropriate entries of IFACE or JFACE, marking each edge as having
!    been used once.
!
!    The third step is to add one more face to our collection of faces.
!    The new face will be adjacent to the current collection of faces,
!    but will include at least one completely unused edge, if possible.
!
!    To guarantee this, we consider every edge that is part of our
!    collection of faces, and that has only been used once.  We look
!    at the endpoints of each of these edges., and
!
!      We search for an adjacent edge that has not been used.
!      If we find such an edge, then the first two edges of our next face
!      are the edge that was already part of the set of faces, and the
!      unused edge.
!
!      If we cannot find such an edge, then we repeat the search, starting
!      with an edge in the face set that has been used once.  But now
!      when we search for adjacent edges, we will consider using one that
!      has already been used once.
!
!    We then search for a path that will return us to the initial
!    node of the first edge.  Using a breadth-first search, we expect
!    to find the shortest path back to that node, and we assume that
!    this represents a face.  Again, we update the IFACE and JFACE arrays.
!
!    We repeat the third step until there are no more edges in the
!    collection of faces that have not been used twice.  Assuming the
!    graph is connected, this means that every face has been found.
!
!  Improvements:
!
!    It shouldn't be hard to modify the code to handle graphs that are
!    not connected.
!
!    If the edge arrays INODE and JNODE were sorted and indexed, some
!    operations could be done more efficiently.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer FACE(MAXORDER,MAXFACE), contains the list of 
!    edges which make up each face.  Face I is made up of the edges
!    FACE(1,I) through FACE(FACE_ORDER(I),I).
!
!    Output, integer FACE_COUNT(NEDGE).  For each edge I, 
!    FACE_COUNT(I) is the number of faces to which the edge belongs.  This 
!    value should be 0, 1 or 2.
!
!    Output, integer IFACE(NEDGE), JFACE(NEDGE).  IFACE(I) and 
!    JFACE(I) are the two faces to which edge I belongs.  Either or both may 
!    be zero if the algorithm fails.
!
!    Input, integer INODE(NEDGE), JNODE(NEDGE), the edge list for 
!    the graph.  The I-th edge connects nodes INODE(I) and JNODE(I).
!
!    Input, integer MAXFACE, the maximum number of faces for which
!    storage has been set aside in FACE and FACE_ORDER.
!
!    Input, integer MAXORDER, the maximum number of edges for 
!    which storage has been set aside in FACE.
!
!    Input, integer NEDGE, the number of edges.
!
!    Output, integer NFACE, the number of faces found by the
!    algorithm.
!
!    Input, integer NNODE, the number of nodes.
!
!    Output, integer FACE_ORDER(MAXFACE).  The number of edges used
!    in constructing face I is stored in FACE_ORDER(I).
!
  implicit none

  logical, parameter :: debug = .false.

  integer maxface
  integer maxorder
  integer nedge
  integer nnode

  integer face(maxorder,maxface)
  integer face_count(nedge)
  integer face_order(maxface)
  integer faceval
  integer i
  integer iedge
  integer iface(nedge)
  integer inode(nedge)
  integer j
  integer jface(nedge)
  integer jnode(nedge)
  integer k
  integer length
  integer nface
  integer nface_old
  integer nodes(3)
  integer nstart
!
!  Initialization.  No arc belongs to any face.
!
  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRAPH_ARC_FACE - Debug:'
    write ( *, '(a)' ) '  Initialization'
  end if

  nface = 0
  face_count(1:nedge) = 0
  iface(1:nedge) = 0
  jface(1:nedge) = 0
  face_order(1:maxface) = 0
  face(1:maxorder,1:maxface) = 0
!
!  We start here.  We may also jump back here if we have used up all the
!  connected parts of a graph, and need to jump to a new piece.
!
!  Find one new face of minimal length.
!
! 5 continue
!
  nface_old = nface

  do length = 3, nnode

    do iedge = 1, nedge

      nodes(1) = inode(iedge)
      nodes(2) = jnode(iedge)
      nstart = 2

      call graph_arc_face_next ( face, face_count, face_order, iface, jface, &  
        inode, jnode, maxface, maxorder, nedge, nface, nnode, nodes, nstart )

      if ( nface_old < nface ) then
        go to 10
      end if

    end do

  end do

  if ( nface == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRAPH_ARC_FACE - Note.'
    write ( *, '(a)' ) '  Could not find any starting face.'
  end if

  go to 60
!
!  Find an edge that is in one face, but not two.
!
10    continue

  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRAPH_ARC_FACE - Debug:'
    write ( *, '(a,i8)' ) '  Found starting face #:', nface
    write ( *, '(a,i8)' ) '  Order is ', face_order(nface)
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Vertices:'
    write ( *, '(a)' ) ' '
    do i = 1, face_order(nface)
      write ( *, '(6i8)' ) face(i,nface)
    end do
  end if

  iedge = 0
!
!  Look for an edge with FACE_COUNT of 1.
!
20    continue

  iedge = iedge + 1

  if ( face_count(iedge) == 1 ) then
    go to 30
  else if ( iedge < nedge ) then
    go to 20
  else
    if ( debug ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GRAPH_ARC_FACE - Debug:'
      write ( *, '(a)' ) '  No more nearby edges to try.'
    end if
!
!  Here, I'd like to be able to jump back and scrounge for other
!  islands of edges, but something's not right.
!
!       go to 5

    go to 60

  end if
!
!  The face will start with the two nodes of edge IEDGE.
!
30    continue

  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRAPH_ARC_FACE - Debug:'
    write ( *, '(a)' ) '  Found a starting edge:'
    write ( *, '(i8)' ) inode(iedge)
    write ( *, '(i8)' ) jnode(iedge)
  end if

  nodes(1) = inode(iedge)
  nodes(2) = jnode(iedge)
!
!  Look for an edge incident to JNODE(IEDGE).  This new edge should have
!  been used just FACEVAL times already.  (FACEVAL is preferably 0, but if
!  we can't find any at that rate, we'll try FACEVAL = 1).
!
  faceval = 0

40    continue

  do i = 1, nedge

    if ( face_count(i) == faceval ) then

      if ( inode(i) == nodes(2) .and. jnode(i) /= nodes(1) ) then
        nodes(3) = jnode(i)
        go to 50
      else if ( jnode(i) == nodes(2) .and. inode(i) /= nodes(1) ) then
        nodes(3) = inode(i)
        go to 50
      end if

    end if

  end do
!
!  If we "fell through" with FACEVAL = 0, then try the search again
!  with FACEVAL = 1.
!
  if ( faceval == 0 ) then

    faceval = 1
    go to 40
!
!  If we fell through with FACEVAL = 1, then we couldn't find any
!  way to use this edge.  Mark it as though it were used, and move on.
!
  else

    if ( debug ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GRAPH_ARC_FACE - Debug:'
      write ( *, '(a)' ) '  Failure.'
      write ( *, '(a,i8)' ) '  Cannot hook up to edge IEDGE:', iedge
      write ( *, '(2i8)' ) nodes(1), nodes(2)
    end if

    face_count(iedge) = 2
    go to 20

  end if
!
!  Now call FACENEXT to search for the shortest cycle that begins
!  NODES(1), NODES(2), NODES(3), and which involves only edges that
!  have been used once or less.
!
50    continue

  nface_old = nface
  nstart = 3

  call graph_arc_face_next ( face, face_count, face_order, iface, jface, &
    inode, jnode, maxface, maxorder, nedge, nface, nnode, nodes, nstart )

  if ( nface_old < nface ) then

    if ( debug ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GRAPH_ARC_FACE - Debug'
      write ( *, '(a,i8)' ) '  NFACE_OLD = ', nface_old
      write ( *, '(a,i8)' ) '  NFACE = ', nface
      write ( *, '(a,i8)' ) '  Order is ', face_order(nface)
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Vertices:'
      write ( *, '(a)' ) ' '
      do i = 1, face_order(nface)
        write ( *, '(6i8)' ) face(i,nface)
      end do
      write ( *, '(a)' ) '  Trying the big loop again.'
    end if

    go to 10

  end if
!
!  The algorithm has failed.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GRAPH_ARC_FACE - Error!'
  write ( *, '(a)' ) '  The algorithm has failed.'
  write ( *, '(a)' ) '  Only some of the faces were found.'
!
!  Cleanup
!
60    continue

  do i = 1, nface
    face_order(i) = min ( face_order(i), maxorder )
  end do

  do i = 1, nface
    do j = 1, face_order(i)
      k = face(j,i)
      if ( k < 0 ) then
        face(j,i) = jnode(-k)
      else
        face(j,i) = inode(k)
      end if
    end do
  end do

  return
end
subroutine graph_arc_face_next ( face, face_count, face_order, iface, jface, &
  inode, jnode, maxface, maxorder, nedge, nface, nnode, nodes, nstart )

!*****************************************************************************80
!
!! graph_arc_face_next() completes the next face, given a few starting nodes.
!
!  Discussion:
!
!    This is a utility routine, called by GRAPH_ARC_FACE, which
!    constructs all the faces of a graph.  GRAPH_ARC_FACE finds the first
!    two or three nodes of a face, and then calls this routine, which
!    attempts to complete the face by using a breadth-first search
!    from the final given node of the face.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer FACE(MAXORDER,MAXFACE), contains the 
!    list of edges which make up each face.  Face I is made up of the edges
!    FACE(1,I) through FACE(FACE_ORDER(I),I).  If a new face is found, this
!    array is updated.
!
!    Input/output, integer FACE_COUNT(NEDGE).  For each edge I, 
!    FACE_COUNT(I) is the number of faces to which the edge belongs.  This value
!    will be 0, 1 or 2.  If a new face is found, this data is updated.
!
!    Input/output, integer FACE_ORDER(MAXFACE).  The number of 
!    edges used in constructing face I is stored in FACE_ORDER(I).
!
!    Input/output, integer IFACE(NEDGE), JFACE(NEDGE).  IFACE(I) 
!    and JFACE(I) are the two faces to which edge I belongs.  Either or both 
!    may be zero if the algorithm fails.  If a new face is found, this data 
!    is updated.
!
!    Input, integer INODE(NEDGE), JNODE(NEDGE), the edge list for 
!    the graph.  The I-th edge connects nodes INODE(I) and JNODE(I).
!
!    Input, integer MAXFACE, the maximum number of faces for which 
!    storage has been set aside in FACE and FACE_ORDER.
!
!    Input, integer MAXORDER, the maximum number of edges for 
!    which storage has been set aside in FACE.
!
!    Input, integer NEDGE, the number of edges.
!
!    Input/output, integer NFACE.  NFACE is the number of faces 
!    found so far.  This value will be updated by this routine if a new face 
!    is found.
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, integer NODES(NSTART), the first few nodes in the 
!    partial cycle.
!
!    Input, integer NSTART, the number of nodes in the partial
!    cycle.
!
!  Workspace:
!
!    Workspace, integer DAD(NNODE), used during the breadth first search
!    of the graph, to point backwards from each node to its predecessor
!    in a path.
!
!    Workspace, integer INDEX(NNODE), used during the breadth first search
!    to label nodes that have been visited.
!
  implicit none

  integer maxface
  integer maxorder
  integer nedge
  integer nnode
  integer nstart

  integer dad(nnode)
  integer face(maxorder,maxface)
  integer face_count(nedge)
  integer face_order(maxface)
  integer i
  integer iedge2
  integer iface(nedge)
  integer index(nnode)
  integer inode(nedge)
  integer istart1
  integer istart2
  integer itemp
  integer jface(nedge)
  integer jnode(nedge)
  integer kedge
  integer kk
  integer nadd
  integer nface
  integer nodei
  integer nodej
  integer npath
  integer nodes(nstart)
!
!  Initialization.
!
  index(1:nnode) = 0
  dad(1:nnode) = 0

  istart1 = nodes(1)
  istart2 = nodes(2)

  do i = 1, nstart

    npath = i

    if ( i == 1 ) then
      dad(nodes(i)) = -1
    else
      dad(nodes(i)) = nodes(i-1)
    end if

    index(nodes(i)) = i

  end do
!
!  From the nodes with INDEX = NPATH, consider all neighbors.
!
10    continue

  npath = npath + 1
  nadd = 0

  do iedge2 = 1, nedge

    if ( index(inode(iedge2)) == npath-1 .and. index(jnode(iedge2)) == 0 ) then

      nodei = inode(iedge2)
      nodej = jnode(iedge2)

    else if ( index(jnode(iedge2)) == npath-1 .and. &
      index(inode(iedge2)) == 0 ) then

      nodei = jnode(iedge2)
      nodej = inode(iedge2)

    else if ( index(inode(iedge2)) == npath-1 .and. &
      jnode(iedge2) == istart1 ) then

      nodei = inode(iedge2)
      nodej = jnode(iedge2)

    else if ( index(jnode(iedge2)) == npath-1 .and. &
      inode(iedge2) == istart1 ) then

      nodei = jnode(iedge2)
      nodej = inode(iedge2)

    else

      nodei = 0
      nodej = 0

    end if

    if ( nodei /= 0 .and. nodej /= istart1 ) then

      nadd = nadd + 1
      index(nodej) = npath
      dad(nodej) = nodei
!
!  Success if the marked node is the starting point (except when
!  using the edge (ISTART2,ISTART1)).
!
    else if ( nodej == istart1 .and. nodei == istart2 ) then

    else if ( nodej == istart1 .and. nodei /= istart2 ) then

      nface = nface + 1

20    continue
!
!  Find the edge KK which joins NODEI and NODEJ.
!
      do kk = 1, nedge

        if ( ( inode(kk) == nodei .and. jnode(kk) == nodej ) .or. &
             ( jnode(kk) == nodei .and. inode(kk) == nodej ) ) then

          face_count(kk) = face_count(kk) + 1
          itemp = face_count(kk)

          if ( itemp == 1 ) then
            iface(kk) = nface
          else if ( itemp == 2 ) then
            jface(kk) = nface
          end if

          if ( inode(kk) == nodei ) then
            kedge = kk
          else
            kedge = -kk
          end if

          exit

        end if

      end do

      nodej = nodei
!
!  Add the edge to the face-to-edge list.
!
      if ( nface <= maxface ) then

        if ( face_order(nface) < maxorder ) then
          face_order(nface) = face_order(nface) + 1
        end if

        if ( face_order(nface) <= maxorder ) then
          face(face_order(nface),nface) = kedge
        end if

      end if

      if ( nodej /= istart1 ) then
        nodei = dad(nodej)
        go to 20
      end if

      return

    end if

  end do
!
!  If we were able to proceed another step, and we haven't exceeded
!  our limit, then go back and take another step.
!
  if ( 0 < nadd .and. npath <= nnode ) then
    go to 10
  end if

  return
end
subroutine graph_arc_is_eulerian ( nnode, nedge, inode, jnode, degree, result )

!*****************************************************************************80
!
!! graph_arc_is_eulerian() determines if a graph is Eulerian from its edge list.
!
!  Definition:
!
!    A graph is Eulerian if there exists a circuit through the graph
!    which uses every edge once.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, integer NEDGE, the number of edges.
!
!    Input, integer INODE(NEDGE), JNODE(NEDGE), the pairs of nodes
!    that form the edges.
!
!    Output, integer DEGREE(NNODE), the degree of each node, that 
!    is, the number of edges that include the node.
!
!    Output, integer RESULT.
!    0, the graph is not Eulerian.
!    1, the graph is Eulerian, but the starting and ending nodes are different.
!    2, the graph is Eulerian, and there is a closed Euler circuit.
!
  implicit none

  integer nedge
  integer nnode

  integer degree(nnode)
  integer i
  integer inode(nedge)
  integer jnode(nedge)
  integer nodd
  integer result

  call graph_arc_degree ( nnode, nedge, inode, jnode, degree )

  nodd = 0

  do i = 1, nnode

    if ( mod ( degree(i), 2 ) == 1 ) then
      nodd = nodd + 1
    end if

  end do

  if ( nodd == 0 ) then
    result = 2
  else if ( nodd == 2 ) then
    result = 1
  else
    result = 0
  end if

  return
end
subroutine graph_arc_is_tree ( nedge, inode, jnode, nnode, result )

!*****************************************************************************80
!
!! graph_arc_is_tree() determines whether a graph is a tree.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 April 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer INODE(NEDGE), JNODE(NEDGE).  INODE(I) and 
!    JNODE(I) are the start and end nodes of the I-th edge of the graph G.
!
!    Input, integer NEDGE, the number of edges in the graph G.
!
!    Input, integer NNODE, the number of nodes.
!
!    Output, integer RESULT.
!    0, the graph is not a tree.
!    1, the graph is a tree.
!
  implicit none

  integer nedge
  integer nnode

  integer adj(nnode,nnode)
  integer inode(nedge)
  integer jnode(nedge)
  integer result

  call graph_arc_to_graph_adj ( nedge, inode, jnode, adj, nnode )

  call graph_adj_is_tree ( adj, nnode, result )

  return
end
subroutine graph_arc_match ( nnode, nedge, inode, jnode, type, match )

!*****************************************************************************80
!
!! graph_arc_match() finds a maximum matching in a bipartite graph.
!
!  Discussion:
!
!    The nodes of the graph are assumed to be divided into two groups,
!    and it is desired to determine as matching that is as large as possible.
!    The matching is a set of pairs ( NODE(I), NODE(J) ) with the properties:
!
!    * NODE(I) is in group 1 and NODE(J) is in group 2;
!    * there is an edge between NODE(I) and NODE(J);
!    * NODE(I) and NODE(J) are not used in any other pairing in the matching.
!
!    The user inputs the edge list that defines the graph, and a set of
!    labels that classify the nodes as being in one group or the other.
!    It is not necessary that the graph actually be bipartite; edges between
!    nodes in the same group are allowed, but they will not affect the
!    outcome in any way.
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
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, integer NEDGE, the number of edges.
!
!    Input, integer INODE(NEDGE), JNODE(NEDGE), the end nodes of
!    the edges.
!
!    Input, integer TYPE(NNODE), labels the two types of nodes 
!    in the graph.  Normally, TYPE(I) would be 0 or 1, but any two distinct
!    values will do.
!
!    Output, integer MATCH(NNODE), the matching node for each
!    node, or 0 if no match was assigned.
!
  implicit none

  integer nedge
  integer nnode

  integer capflo(2,2*nedge+2*nnode)
  integer i
  integer icut(nnode+2)
  integer iendpt(2,2*nedge+2*nnode)
  integer in
  integer inode(nedge)
  integer isink
  integer isorce
  integer j
  integer jn
  integer jnode(nedge)
  integer match(nnode)
  integer nedge2
  integer nnode2
  integer node_flow(nnode+2)
  integer type(nnode)

  match(1:nnode) = 0
!
!  Create a network from the graph, with two extra nodes.
!
  isorce = nnode + 1
  isink = nnode + 2
  nnode2 = nnode + 2

  j = 0

  do i = 1, nedge

    in = inode(i)
    jn = jnode(i)

    if ( type(in) /= type(jn) ) then

      j = j + 1
      iendpt(1,j) = inode(i)
      iendpt(2,j) = jnode(i)
      capflo(1,j) = 1
      capflo(2,j) = 0

      j = j + 1
      iendpt(1,j) = jnode(i)
      iendpt(2,j) = inode(i)
      capflo(1,j) = 1
      capflo(2,j) = 0

    end if

  end do
!
!  Nodes of type 1 are connected to the source, 
!  and nodes of type 2 are connected to the sink.
!
  do i = 1, nnode

    if ( type(i) == type(1) ) then
      j = j + 1
      iendpt(1,j) = isorce
      iendpt(2,j) = i
      capflo(1,j) = 1
      capflo(2,j) = 0
      j = j + 1
      iendpt(1,j) = i
      iendpt(2,j) = isorce
      capflo(1,j) = 1
      capflo(2,j) = 0
    else
      j = j + 1
      iendpt(1,j) = i
      iendpt(2,j) = isink
      capflo(1,j) = 1
      capflo(2,j) = 0
      j = j + 1
      iendpt(1,j) = isink
      iendpt(2,j) = i
      capflo(1,j) = 1
      capflo(2,j) = 0
    end if
  end do
!
!  Determine the maximum flow on the network.
!
!  Then a pair of nodes connected by an edge that has a network flow of 1
!  are part of the maximal matching.
!
  nedge2 = j

  call network_flow_max ( nnode2, nedge2, iendpt, capflo, isorce, isink, &
    icut, node_flow )

  do i = 1, nedge2
    if ( iendpt(1,i) <= nnode .and. &
         iendpt(2,i) <= nnode .and. &
         0 < capflo(1,i) .and. &
         capflo(2,i) == 1 ) then
      in = iendpt(1,i)
      jn = iendpt(2,i)
      match(in) = jn
      match(jn) = in
    end if
  end do

  return
end
subroutine graph_arc_min_path ( nnode, nedge, inode, jnode, arcost, &
  istart, last, num_path, ispath, xlen )

!*****************************************************************************80
!
!! graph_arc_min_path() finds the shortest path between two nodes.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 September 1999
!
!  Reference:
!
!    Hang Tong Lau,
!    Combinatorial Heuristic Algorithms in FORTRAN,
!    Springer Verlag, 1986.
!
!  Parameters:
!
!    Input, integer NNODE, the number of nodes in the graph.
!
!    Input, integer NEDGE, the number of edges in the graph.
!
!    Input, integer INODE(NEDGE), JNODE(NEDGE), the edges of the
!    graph, describe by pairs of nodes.
!
!    Input, real ( kind = rk ) ARCOST(NEDGE), the distance or cost of each edge.
!
!    Input, integer ISTART, LAST, are the two nodes between which a
!    shortest path is desired.
!
!    Output, integer NUM_PATH, the number of nodes in the shortest
!    path.  NUM_PATH is zero if no path could be found.
!
!    Output, integer ISPATH(NNODE), lists the nodes in the
!    shortest path, from ISPATH(1) to ISPATH(NUM_PATH).
!
!    Output, real ( kind = rk ) XLEN, the length of the shortest path 
!    from ISTART to LAST.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nedge
  integer nnode

  real ( kind = rk ) arcost(nedge)
  real ( kind = rk ) d
  integer i
  integer ic
  integer ient
  logical ifin
  integer ij
  integer inode(nedge)
  integer ispath(nnode)
  integer istart
  logical iwork1(nnode)
  integer iwork2(nnode)
  integer iwork3(nedge)
  integer j
  integer jnode(nedge)
  integer k
  integer l
  integer last
  integer num_path
  real ( kind = rk ) wk4(nnode)
  real ( kind = rk ) xlen

  wk4(1:nnode) = huge ( wk4(1) )
  iwork1(1:nnode) = .true.
  iwork2(1:nnode) = 0

  wk4(istart) = 0.0D+00
  i = istart
  iwork1(istart) = .false.
  xlen = 0
!
!  For each forward arc originating at node I calculate
!  the length of the path to node I.
!
10    continue

  ic = 0

  do k = 1, nedge

    if ( inode(k) == i ) then
      ic = ic + 1
      iwork3(ic) = k
      ispath(ic) = jnode(k)
    end if

    if ( jnode(k) == i ) then
      ic = ic + 1
      iwork3(ic) = k
      ispath(ic) = inode(k)
    end if

  end do

  if ( 0 < ic ) then

    do l = 1, ic
      k = iwork3(l)
      j = ispath(l)
      if ( iwork1(j) ) then
        d = wk4(i) + arcost(k)
        if ( d < wk4(j) ) then
          wk4(j) = d
          iwork2(j) = k
        end if
      end if
    end do

  end if
!
!  Find the minimum potential.
!
  d = huge ( d )
  ient = 0
  ifin = .false.

  do i = 1, nnode

    if ( iwork1(i) ) then
      ifin = .true.
      if ( wk4(i) < d ) then
        d = wk4(i)
        ient = i
      end if
    end if

  end do
!
!  Include the node in the current path.
!
  if ( d < huge ( d ) ) then
    iwork1(ient) = .false.
    if ( ient /= last ) then
      i = ient
      go to 10
    end if
  else
    if ( ifin ) then
      num_path = 0
      return
    end if
  end if

  ij = last
  num_path = 1
  ispath(1) = last

  do

    k = iwork2(ij)

    if ( inode(k) == ij ) then
      ij = jnode(k)
    else
      ij = inode(k)
    end if

    num_path = num_path + 1
    ispath(num_path) = ij

    if ( ij == istart ) then
      exit
    end if

  end do

  l = num_path / 2
  j = num_path

  do i = 1, l
    call i4_swap ( ispath(i), ispath(j) )
    j = j - 1
  end do

  xlen = wk4(last)

  return
end
subroutine graph_arc_min_span_tree ( nnode, nedge, inode, jnode, cost, &
  itree, jtree, tree_cost )

!*****************************************************************************80
!
!! graph_arc_min_span_tree() finds a minimum spanning tree of a graph.
!
!  Discussion:
!
!    The input graph is represented by a list of edges.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 July 2000
!
!  Reference:
!
!    Hang Tong Lau,
!    Combinatorial Heuristic Algorithms in FORTRAN,
!    Springer Verlag, 1986.
!
!  Parameters:
!
!    Input, integer NNODE, the number of nodes in the graph.
!
!    Input, integer NEDGE, the number of edges in the graph.
!
!    Input, integer INODE(NEDGE), JNODE(NEDGE), the start and end 
!    nodes of the edges.
!
!    Input, real ( kind = rk ) COST(NEDGE), the cost or length of each edge.
!
!    Output, integer ITREE(NNODE-1), JTREE(NNODE-1), the pairs 
!    of nodes that form the edges of the spanning tree.
!
!    Output, real ( kind = rk ) TREE_COST, the total cost or length 
!    of the spanning tree.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nedge
  integer nnode

  integer best
  real ( kind = rk ) cost(nedge)
  real ( kind = rk ) d
  logical free(nnode)
  integer i
  integer ic
  integer ij
  integer inode(nedge)
  integer itr
  integer itree(nnode-1)
  integer iwork1(nnode)
  integer iwork2(nnode)
  integer iwork4(nedge)
  integer iwork5(nedge)
  integer j
  integer jnode(nedge)
  integer jtree(nnode-1)
  integer jj
  integer k
  integer kk
  integer l
  real ( kind = rk ) tree_cost
  real ( kind = rk ) wk6(nnode)

  wk6(1:nnode) = huge ( wk6(1) )
  free(1:nnode) = .true.
  iwork1(1:nnode) = 0
  iwork2(1:nnode) = 0
  itree(1:nnode-1) = 0
  jtree(1:nnode-1) = 0
!
!  Find the first non-zero arc.
!
  do ij = 1, nedge
    if ( inode(ij) /= 0 ) then
      i = inode(ij)
      exit
    end if
  end do

  wk6(i) = 0.0D+00
  free(i) = .false.

  tree_cost = 0.0D+00

  do jj = 1, nnode - 1

    wk6(1:nnode) = huge ( wk6(1) )

    do i = 1, nnode
!
!  For each forward arc originating at node I
!  calculate the length of the path to node I.
!
      if ( .not. free(i) ) then

        ic = 0

        do k = 1, nedge

          if ( inode(k) == i ) then
            ic = ic + 1
            iwork5(ic) = k
            iwork4(ic) = jnode(k)
          end if

          if ( jnode(k) == i ) then
            ic = ic + 1
            iwork5(ic) = k
            iwork4(ic) = inode(k)
          end if

        end do

        if ( 0 < ic ) then

          do l = 1, ic

            k = iwork5(l)
            j = iwork4(l)

            if ( free(j) ) then

              d = tree_cost + cost(k)

              if ( d < wk6(j) ) then
                wk6(j) = d
                iwork1(j) = i
                iwork2(j) = k
              end if

            end if

          end do

        end if

      end if

    end do
!
!  Identify the free node of minimum potential.
!
    d = huge ( d )
    best = 0

    do i = 1, nnode

      if ( free(i) ) then
        if ( wk6(i) < d ) then
          d = wk6(i)
          best = i
          itr = iwork1(i)
          kk = iwork2(i)
        end if
      end if

    end do
!
!  Add that node to the tree.
!
    if ( d < huge ( d ) ) then
      free(best) = .false.
      tree_cost = tree_cost + cost(kk)
      itree(jj) = itr
      jtree(jj) = best
    end if

  end do

  return
end
subroutine graph_arc_ncolor_print ( nedge, inode, jnode, nnode, color, title )

!*****************************************************************************80
!
!! graph_arc_ncolor_print() prints out a node-colored graph from an edge list.
!
!  Discussion:
!
!    The printout is arranged to emphasize the colors of the neighboring nodes.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NEDGE, the number of edges.
!
!    Input, integer INODE(NEDGE), JNODE(NEDGE), the beginning 
!    and end nodes of the edges.
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, integer COLOR(NNODE), the color of each node.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer nedge
  integer nnode

  integer color(nnode)
  integer i
  integer in
  integer inode(nedge)
  integer jn
  integer jnode(nedge)
  character ( len = * ) title

  if ( len_trim ( title ) /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Edge  Node 1  Node 2     Color 1 Color 2'
  write ( *, '(a)' ) ' '

  do i = 1, nedge
    in = inode(i)
    jn = jnode(i)
    write ( *, '(i8,2x,i8,2x,i8,2x,i8,2x,i8)' ) i, in, jn, color(in), color(jn)
  end do

  return
end
subroutine graph_arc_node_count ( nedge, inode, jnode, mnode, nnode )

!*****************************************************************************80
!
!! graph_arc_node_count() counts the number of nodes in a graph.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    24 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NEDGE, the number of edges.
!
!    Input, integer INODE(NEDGE), JNODE(NEDGE).  INODE(I) and 
!    JNODE(I) are the start and end nodes of the I-th edge.
!
!    Output, integer MNODE, the maximum node index.
!
!    Output, integer NNODE, the number of distinct nodes.
!
  implicit none

  integer nedge

  integer iedge
  integer inode(nedge)
  integer jnode(nedge)
  integer knode(2*nedge)
  integer mnode
  integer nnode

  mnode = max ( maxval ( inode(1:nedge) ), maxval ( jnode(1:nedge) ) )
!
!  Copy all the node labels into KNODE,
!  sort KNODE,
!  count the unique entries.  
!
!  That's NNODE.
!
  knode(1:nedge) = inode(1:nedge)

  do iedge = 1, nedge
    knode(nedge+iedge) = jnode(iedge)
  end do

  call i4vec_sort_heap_a ( 2*nedge, knode )

  call i4vec_uniq ( 2*nedge, knode, nnode )

  return
end
subroutine graph_arc_print ( nedge, inode, jnode, title )

!*****************************************************************************80
!
!! graph_arc_print() prints out a graph from an edge list.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NEDGE, the number of edges.
!
!    Input, integer INODE(NEDGE), JNODE(NEDGE), the beginning 
!    and end nodes of the edges.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer nedge

  integer i
  integer inode(nedge)
  integer jnode(nedge)
  character ( len = * ) title

  if ( len_trim ( title ) /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  do i = 1, nedge
    write ( *, '(i8,4x,2i8)' ) i, inode(i), jnode(i)
  end do

  return
end
subroutine graph_arc_to_ps ( file_name, inode, jnode, nedge, nnode, x, y )

!*****************************************************************************80
!
!! graph_arc_to_ps() writes graph information to a PostScript file.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the output file.
!
!    Input, integer INODE(NEDGE), JNODE(NEDGE), the edge array.
!    The I-th edge connects nodes INODE(I) and JNODE(I).
!
!    Input, integer NEDGE, the number of edges.
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, real ( kind = rk ) X(NNODE), Y(NNODE), the X and Y components
!    of points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nedge
  integer nnode

  real ( kind = rk ) alpha
  real ( kind = rk ) blue
  character ( len = 8 ) date
  character ( len = * ) file_name
  real ( kind = rk ) green
  integer inode(nedge)
  integer ios
  integer iunit
  integer jnode(nedge)
  integer margin
  integer pagexmax
  integer pagexmin
  integer pageymax
  integer pageymin
  integer plotxmax
  integer plotxmin
  integer plotxmin2
  integer plotymax
  integer plotymin
  integer plotymin2
  real ( kind = rk ) red
  real ( kind = rk ) x(nnode)
  real ( kind = rk ) xmax
  real ( kind = rk ) xmin
  real ( kind = rk ) y(nnode)
  real ( kind = rk ) ymax
  real ( kind = rk ) ymin
!
!  Bounding box.
!
  xmin = minval ( x(1:nnode) )
  xmax = maxval ( x(1:nnode) )
  ymin = minval ( y(1:nnode) )
  ymax = maxval ( y(1:nnode) )

  if ( xmin == xmax ) then
    xmin = x(1) - 0.5D+00
    xmax = x(1) + 0.5D+00
  end if

  if ( ymin == ymax ) then
    ymin = y(1) - 0.5D+00
    ymax = y(1) + 0.5D+00
  end if
!
!  Compute the scale factor.
!
  pagexmax = 612
  pagexmin = 0
  pageymax = 792
  pageymin = 0

  margin = 36

  plotxmax = pagexmax - margin
  plotxmin = pagexmin + margin
  plotymax = pageymax - margin
  plotymin = pageymin + margin

  alpha = min ( ( plotxmax - plotxmin ) / ( xmax - xmin ), &
                ( plotymax - plotymin ) / ( ymax - ymin ) )
!
!  Adjust PLOTXMIN and PLOTYMIN to center the image.
!
  plotxmin2 = int ( 0.5D+00 * ( plotxmin + plotxmax - alpha * ( xmax - xmin ) ) )
  plotymin2 = int ( 0.5D+00 * ( plotymin + plotymax - alpha * ( ymax - ymin ) ) )

  call get_unit ( iunit )

  open ( unit = iunit, file = file_name, status = 'replace', iostat = ios )

  if ( ios /= 0 ) then
    return
  end if
!
!  Write the prolog.
!
  write ( iunit, '(a)' ) '%!PS-Adobe-3.0'
  write ( iunit, '(a)' ) '%%Document-Fonts: Times-Roman'
  write ( iunit, '(a,a)' ) '%%Title: ' , trim ( file_name )
  write ( iunit, '(a)' ) '%%Creator: GRAFPACK(graph_arc_to_ps)'
  call date_and_time ( date )
  write ( iunit, '(a)' ) '%%CreationDate: ' // trim ( date )
  write ( iunit, '(a,4i5)' ) '%%BoundingBox', plotxmin, plotymin, plotxmax, &
    plotymax
  write ( iunit, '(a)' ) '%%LanguageLevel: 2'
  write ( iunit, '(a)' ) '%%EndComments'
  write ( iunit, '(a)' ) '%%BeginProlog'
  write ( iunit, '(a)' ) '%%EndProlog'
!
!  Set the line color.
!
  red = 0.0D+00
  green = 0.0D+00
  blue = 0.0D+00

  write ( iunit, '(3f7.4,a)' ) red, green, blue, ' setrgbcolor'
!
!  Draw lines.
!
  call edges_to_ps ( plotxmin2, plotymin2, alpha, iunit, inode, jnode, &
    nedge, nnode, x, y, xmin, ymin )
!
!  Set the fill color.
!
  red = 0.1
  green = 0.1
  blue = 0.7

  write ( iunit, '(3f7.4,a)' ) red, green, blue, ' setrgbcolor'
!
!  Draw points.
!
  call nodes_to_ps ( plotxmin2, plotymin2, alpha, iunit, nnode, x, y, &
    xmin, ymin )

  write ( iunit, '(a)' ) 'showpage'
!
!  Write the epilog.
!
  write ( iunit, '(a)' ) 'grestore'
  write ( iunit, '(a)' ) '%%Trailer'
  write ( iunit, '(a,i2)' ) '%%Pages: 1'

  close ( iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GRAPH_ARC_TO_PS'
  write ( *, '(a)' ) '  The data was written to the file: ' &
    // trim ( file_name )

  return
end
subroutine graph_arc_span_forest ( nnode, nedge, inode, jnode, ncomp, &
  component )

!*****************************************************************************80
!
!! graph_arc_span_forest() determines a graph's connectivity and spanning forest.
!
!  Discussion:
!
!    A (connected) component of a graph is a maximal subgraph which
!    is connected.
!
!    A tree is a connected graph containing no cycles.
!
!    A spanning tree of a connected graph is a subgraph which is a 
!    maximal tree.
!
!    A forest is a collection of trees, no two of which share a node.
!
!    A spanning forest of a possibly unconnected graph is a collection
!    containing a single spanning tree for each component of the graph.
!
!    The input graph may be connected or unconnected.
!
!    If the input graph is connected, this routine simply returns a
!    spanning tree for the graph.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    29 October 1999
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, integer NEDGE, the number of edges in the graph.
!
!    Input/output, integer INODE(NEDGE), JNODE(NEDGE), the edge 
!    list of the graph.  On output, this array has been rearranged.  Edges
!    belonging to the spanning tree of component 1 are first, followed
!    by edges belonging to the other spanning trees, followed last by
!    edges that were not used in any spanning tree.
!
!    Output, integer NCOMP, the number of connected components 
!    of the graph.
!
!    Input/output, integer IENDPT(2,NEDGE), the edge array of the 
!    graph.  IENDPT(1,I) and IENDPT(2,I) are the two nodes that make up edge I.
!  
!    On input, IENDPT describes the graph.
!
!    On output, the input entries of IENDPT have been reordered, so that
!    edges belonging to the spanning forest come first, followed by those
!    edges which are not part of the spanning forest.
!
!    Output, integer NCOMP, the number of connected components of 
!    the graph.
!
!    Output, integer IARRAY(NNODE).  IARRAY(I) is the component to 
!    which node I belongs, and will take on values between 1 and NCOMP.
!
  implicit none

  integer nedge
  integer nnode

  integer component(nnode)
  integer i
  integer inode(nedge)
  integer inode2(nedge)
  integer j
  integer jnode(nedge)
  integer jnode2(nedge)
  integer left
  integer ncomp
  integer nstack
  integer num
  integer prev
  integer r
  integer right
  integer stack_node(nnode)
  integer stack_prev(nnode)
  integer v
  integer x_num(nnode)

  left = 0
  right = nedge + 1
  inode2(1:nedge) = 0
  jnode2(1:nedge) = 0
!
!  Part A:
!
  component(1:nnode) = 0
  x_num(1:nnode) = 0
  ncomp = 0
  v = 1
  num = 0

  nstack = 0
!
!  Part B:
!  Scan next V.
!
10    continue

  if ( nnode < v ) then
    inode(1:nedge) = inode2(1:nedge)
    jnode(1:nedge) = jnode2(1:nedge)
    return
  end if

  if ( component(v) /= 0 ) then
    v = v + 1
    go to 10
  end if
!
!  Begin the NCOMP-th component at V.
!
  ncomp = ncomp + 1
  num = num + 1

  component(v) = ncomp
  x_num(v) = num

  nstack = nstack + 1
  stack_node(nstack) = v
  stack_prev(nstack) = 0
!
!  Part C:
!  Is component NCOMP finished?
!
  do

    if ( nstack <= 0 ) then
      v = v + 1
      go to 10
    end if

    j = stack_node(nstack)
    prev = stack_prev(nstack)
    nstack = nstack - 1
!
!  Examine each vertex R that is adjacent to node J.
!
    do i = 1, nedge

      if ( inode(i) == j ) then
        r = jnode(i)
      else if ( jnode(i) == j ) then
        r = inode(i)
      else
        r = 0
      end if

      if ( r /= 0 ) then

        if ( component(r) == 0 ) then

          num = num + 1
          component(r) = ncomp
          x_num(r) = num

          nstack = nstack + 1
          stack_node(nstack) = r
          stack_prev(nstack) = j

          left = left + 1
          inode2(left) = j
          jnode2(left) = r

        else

          if ( r == prev .or. x_num(j) < x_num(r) ) then

          else

            right = right - 1
            inode2(right) = j
            jnode2(right) = r

          end if

        end if

      end if

    end do

  end do

  return
end
subroutine graph_arc_span_tree ( nedge, inode, jnode, nnode, dad )

!*****************************************************************************80
!
!! graph_arc_span_tree() constructs the spanning tree of a graph.
!
!  Discussion:
!
!    If the graph is connected, then exactly one node will have no
!    parent, and a DAD value of -1.
!
!    If the graph is not connected, but divided into NCOMP components, then 
!    NCOMP nodes will have a DAD value of -1.
!
!    If the graph is connected, then once the tree is computed, the 
!    addition to the tree of any edge not included in the tree will 
!    form a cycle.  Since there are NNODE-1 edges in the tree, this will 
!    normally mean that there are NEDGE-(NNODE-1) "fundamental" cycles 
!    that can be generated in this way.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    24 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NEDGE, the number of edges.
!
!    Input, integer INODE(NEDGE), JNODE(NEDGE), the edge array 
!    of the graph.  The I-th edge joins nodes INODE(I) and JNODE(I).
!
!    Input, integer NNODE, the number of nodes.
!
!    Output, integer DAD(NNODE), the "father" array.  If node I is 
!    the root of the tree spanning a given component of the graph, then 
!    DAD(I) = -1.  Otherwise, DAD(I) is the index of another node J in
!    the same component, such that the edge (I,J) is the first step
!    in the path along the tree from node I to the root of its component.
! 
  implicit none

  integer nedge
  integer nnode

  integer dad(nnode)
  integer i
  integer iedge
  integer inode(nedge)
  integer jnode(nedge)
  integer nodei
  integer nodej
  integer nstacki
  integer nstackj
  integer stack(nnode)
!
!  Initialize.
!
  nstacki = 0
  nstackj = 0

  dad(1:nnode) = 0
  stack(1:nnode) = 0
!
!  Start at an unvisited node.
!
  do

    i = 0

    do

      i = i + 1

      if ( nnode < i ) then
        return
      end if

      if ( dad(i) == 0 ) then
        exit
      end if

    end do

    nodei = i
    dad(nodei) = - 1

    nstacki = 1
    stack(nstacki) = nodei
!
!  Search for unvisited neighbors of the last set of nodes.
!
    do

      do i = 1, nstacki

        nodei = stack(i)

        do iedge = 1, nedge

          if ( inode(iedge) == nodei ) then
            nodej = jnode(iedge)
          else if ( jnode(iedge) == nodei ) then
            nodej = inode(iedge)
          else
            nodej = 0
          end if
!
!  Store unvisited neighbors in STACK.
!
          if ( nodej /= 0 ) then

            if ( dad(nodej) == 0 ) then
              dad(nodej) = nodei
              nstackj = nstackj + 1
              stack(nstacki+nstackj) = nodej
            end if

          end if
 
        end do

      end do
!
!  If you picked up new neighbors on this pass, then we need to
!  search for THEIR neighbors.
!
      if ( nstackj <= 0 ) then
        exit
      end if

      stack(1:nstackj) = stack(nstacki+1:nstacki+nstackj)
      nstacki = nstackj
      nstackj = 0

    end do

  end do

  return
end
subroutine graph_arc_to_digraph_arc ( iarc, jarc, inode, jnode, maxarc, narc, &
  nedge )

!*****************************************************************************80
!
!! graph_arc_to_digraph_arc() converts an undirected to a directed graph.
!
!  Discussion:
!
!    The intent is that every edge (I,J) of the undirected graph will 
!    become two directed edges or "arcs" (I,J) and (J,I) of the directed
!    graph.  An "arc" (I,J) is a path FROM I TO J, and does not allow
!    passage back from J to I.
!
!    An edge (I,I) results in a single arc (I,I).
!
!    If the input data already includes edges (I,J) and (J,I), then 
!    the code will catch this fact, and will produce two arcs, not four.
!
!    As part of the processing, edges (I,J) in the input array are 
!    reordered if necessary so that I <= J.  Then the edge array is 
!    sorted, and duplicates are removed.  Only then are the arcs 
!    generated.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IARC(MAXARC), JARC(MAXARC), the arcs of a 
!    directed graph, with the property that every edge (I,J) in the undirected 
!    graph corresponds to a pair of arcs (I,J) and (J,I) in the directed 
!    graph, with the exception that an edge (I,I) corresponds to a single 
!    arc (I,I).  The I-th arc goes from IARC(I) to JARC(I).
!
!    Input, integer INODE(NEDGE), JNODE(NEDGE), the edge array for
!    an undirected graph.  The I-th edge connects nodes INODE(I) and JNODE(I).
!
!    Input, integer MAXARC, the maximum number of arcs for which 
!    storage has been set aside.  MAXARC = 2*NEDGE is always enough, but less
!    space may be required if there are many duplicate edges, or 
!    edges of the form (I,I).
!
!    Output, integer NARC, the actual number of arcs constructed 
!    for the directed graph.
!
!    Input, integer NEDGE, the number of edges in the undirected 
!    graph.
!
  implicit none

  integer maxarc
  integer nedge

  integer i
  integer iarc(maxarc)
  integer inode(nedge)
  integer jarc(maxarc)
  integer jnode(nedge)
  integer narc
  integer nuniq
!
!  Copy the edge array into the initial part of the arc array.
!
  narc = nedge

  iarc(1:narc) = inode(1:narc)
  jarc(1:narc) = jnode(1:narc)
!
!  Sort the edge array as though it were undirected.
!
  call graph_arc_edge_sort ( narc, iarc, jarc )
!
!  Eliminate duplicates.
!
  call i4vec2_uniq ( narc, iarc, jarc, nuniq )
!
!  Generate the extra arcs.
!
  narc = nuniq

  do i = 1, nuniq

    if ( iarc(i) /= jarc(i) ) then

      narc = narc + 1

      if ( narc <= maxarc ) then
        iarc(narc) = jarc(i)
        jarc(narc) = iarc(i)
      end if

    end if
 
  end do
!
!  Now sort the digraph edge array.
!
  call digraph_arc_edge_sort ( narc, iarc, jarc )

  return
end
subroutine graph_arc_to_graph_adj ( nedge, inode, jnode, adj, nnode )

!*****************************************************************************80
!
!! graph_arc_to_graph_adj() converts an arc list graph to an adjacency graph.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NEDGE, the number of edges.
!
!    Input, integer INODE(NEDGE), JNODE(NEDGE), the edge array for 
!    an undirected graph.  The I-th edge connects nodes INODE(I) and JNODE(I).
!
!    Output, integer ADJ(NNODE,NNODE), the adjacency information.
!
!    Input, integer NNODE, the number of nodes.
!
  implicit none

  integer nnode
  integer nedge

  integer adj(nnode,nnode)
  integer i
  integer inode(nedge)
  integer j
  integer jnode(nedge)
  integer k

  adj(1:nnode,1:nnode) = 0

  do k = 1, nedge
    i = inode(k)
    j = jnode(k)
    adj(i,j) = 1
    adj(j,i) = 1
  end do

  return
end
subroutine graph_arc_to_graph_star ( nnode, nedge, inode, jnode, arcfir, &
  fwdarc )

!*****************************************************************************80
!
!! graph_arc_to_graph_star() sets the forward star form of an undirected graph.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, integer NEDGE, the number of edges.
!
!    Input, integer INODE(NEDGE), JNODE(NEDGE); the I-th edge of 
!    the graph extends from node INODE(I) to JNODE(I).
!
!    Output, integer ARCFIR(NNODE+1); ARCFIR(I) is the number of 
!    the first edge starting at node I in the forward star representation of 
!    the graph.
!
!    Output, integer FWDARC(2*NEDGE); FWDARC(I) is the ending node 
!    of the I-th edge in the forward star representation of the graph.
!
  implicit none

  integer nedge
  integer nnode

  integer arcfir(nnode+1)
  integer fwdarc(2*nedge)
  integer i
  integer inode(nedge)
  integer j
  integer jnode(nedge)
  integer k
!
!  Set up the forward star representation of the graph.
!
  k = 0

  do i = 1, nnode

    arcfir(i) = k + 1

    do j = 1, nedge

      if ( inode(j) == i ) then
        k = k + 1
        fwdarc(k) = jnode(j)
      end if

      if ( jnode(j) == i ) then
        k = k + 1
        fwdarc(k) = inode(j)
      end if

    end do

  end do

  arcfir(nnode+1) = k + 1

  return
end
subroutine graph_arc_weight_print ( nedge, inode, jnode, wnode, title )

!*****************************************************************************80
!
!! graph_arc_weight_print() prints out a weighted graph from an edge list.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer NEDGE, the number of edges.
!
!    integer INODE(NEDGE), JNODE(NEDGE), the beginning 
!    and end nodes of the edges.
!
!    real ( kind = rk ) WNODE(NEDGE), the weights of the edges.
!
!    character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nedge

  integer i
  integer inode(nedge)
  integer jnode(nedge)
  character ( len = * ) title
  real ( kind = rk ) wnode(nedge)

  if ( len_trim ( title ) /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  do i = 1, nedge
    write ( *, '(i8,4x,2i8,g14.6)' ) i, inode(i), jnode(i), wnode(i)
  end do

  return
end
subroutine graph_dist_print ( dist, nnode, title )

!*****************************************************************************80
!
!! graph_dist_print() prints a distance matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) DIST(NNODE,NNODE), the distance matrix.  
!    DIST(I,J) is the distance from node I to node J.
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nnode

  real ( kind = rk ) dist(nnode,nnode)
  integer ihi
  integer ilo
  integer jhi
  integer jlo
  integer ncol
  integer nrow
  character ( len = * ) title

  if ( len_trim ( title ) /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  ilo = 1
  ihi = nnode
  jlo = 1
  jhi = nnode
  ncol = nnode
  nrow = nnode

  call r8mat_print ( dist, ihi, ilo, jhi, jlo, ncol, nrow )

  return
end
subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! i4_swap() switches two integer values.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer I, J.  On output, the values of I and
!    J have been interchanged.
!
  implicit none

  integer i
  integer j
  integer k

  k = i
  i = j
  j = k

  return
end
subroutine i4vec_backtrack ( n, x, indx, k, nstack, stack, maxstack, ncan )

!*****************************************************************************80
!
!! i4vec_backtrack() supervises a backtrack search for an integer vector.
!
!  Discussion:
!
!    The routine tries to construct an integer vector one index at a time,
!    using possible candidates as supplied by the user.
!
!    At any time, the partially constructed vector may be discovered to be
!    unsatisfactory, but the routine records information about where the
!    last arbitrary choice was made, so that the search can be
!    carried out efficiently, rather than starting out all over again.
!
!    First, call the routine with INDX = 0 so it can initialize itself.
!
!    Now, on each return from the routine, if INDX is:
!      1, you've just been handed a complete candidate vector;
!         Admire it, analyze it, do what you like.
!      2, please determine suitable candidates for position X(K).
!         Return the number of candidates in NCAN(K), adding each
!         candidate to the end of STACK, and increasing NSTACK.
!      3, you're done.  Stop calling the routine;
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    24 July 2000
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer N, the number of positions to be filled in the 
!    vector.
!
!    Input/output, integer X(N), the partial or complete candidate 
!    vector.
!
!    Input/output, integer INDX, a communication flag.
!    On input,
!      0 to start a search.
!    On output:
!      1, a complete output vector has been determined and returned in X(1:N);
!      2, candidates are needed for position X(K);
!      3, no more possible vectors exist.
!
!    Output, integer K, if INDX=2, the current vector index being 
!    considered.
!
!    Input/output, integer NSTACK, the current length of the stack.
!
!    Input, integer STACK(MAXSTACK), a list of all current 
!    candidates for all positions 1 through K.
!
!    Input, integer MAXSTACK, the maximum length of the stack.
!
!    Input/output, integer NCAN(N), lists the current number of 
!    candidates for positions 1 through K.
!
  implicit none

  integer n
  integer maxstack

  integer indx
  integer k
  integer ncan(n)
  integer nstack
  integer stack(maxstack)
  integer x(n)
!
!  If this is the first call, request a candidate for position 1.
!
  if ( indx == 0 ) then
    k = 1
    nstack = 0
    indx = 2
    return
  end if
!
!  Examine the stack.
!
  do
!
!  If there are candidates for position K, take the first available
!  one off the stack, and increment K.
!
!  This may cause K to reach the desired value of N, in which case
!  we need to signal the user that a complete set of candidates
!  is being returned.
!
    if ( 0 < ncan(k) ) then

      x(k) = stack(nstack)
      nstack = nstack - 1

      ncan(k) = ncan(k) - 1

      if ( k /= n ) then
        k = k + 1
        indx = 2
      else
        indx = 1
      end if

      exit
!
!  If there are no candidates for position K, then decrement K.
!  If K is still positive, repeat the examination of the stack.
!
    else

      k = k - 1

      if ( k <= 0 ) then
        indx = 3
        exit
      end if

    end if

  end do

  return
end
subroutine i4vec_heap_d ( n, a )

!*****************************************************************************80
!
!! i4vec_heap_d() reorders an array of integers into an descending heap.
!
!  Definition:
!
!    A descending heap is an array A with the property that, for every index J,
!    A(J) >= A(2*J) and A(J) >= A(2*J+1), (as long as the indices
!    2*J and 2*J+1 are legal).
!
!  Diagram:
!
!                  A(1)
!                /      \
!            A(2)         A(3)
!          /     \        /  \
!      A(4)       A(5)  A(6) A(7)
!      /  \       /   \
!    A(8) A(9) A(10) A(11)
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer N, the size of the input array.
!
!    Input/output, integer A(N).
!    On input, an unsorted array.
!    On output, the array has been reordered into a heap.
!
  implicit none

  integer n

  integer a(n)
  integer i
  integer ifree
  integer key
  integer m
!
!  Only nodes N/2 down to 1 can be "parent" nodes.
!
  do i = n/2, 1, -1
!
!  Copy the value out of the parent node.
!  Position IFREE is now "open".
!
    key = a(i)
    ifree = i

10  continue
!
!  Positions 2*IFREE and 2*IFREE + 1 are the descendants of position
!  IFREE.  (One or both may not exist because they exceed N.)
!
    m = 2 * ifree
!
!  Does the first position exist?
!
    if ( m <= n ) then
!
!  Does the second position exist?
!
      if ( m + 1 <= n ) then
!
!  If both positions exist, take the larger of the two values,
!  and update M if necessary.
!
        if ( a(m) < a(m+1) ) then
          m = m + 1
        end if

      end if
!
!  If the large descendant is larger than KEY, move it up,
!  and update IFREE, the location of the free position, and
!  consider the descendants of THIS position.
!
      if ( key < a(m) ) then
        a(ifree) = a(m)
        ifree = m
        go to 10
      end if

    end if
!
!  Once there is no more shifting to do, the value KEY
!  moves into the free spot IFREE.
!
    a(ifree) = key

  end do

  return
end
subroutine i4vec_indicator ( n, a )

!*****************************************************************************80
!
!! i4vec_indicator() sets an integer vector to the indicator vector.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements of A.
!
!    Output, integer A(N), the array to be initialized.
!
  implicit none

  integer n

  integer a(n)
  integer i

  do i = 1, n
    a(i) = i
  end do

  return
end
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! i4vec_print() prints an integer vector.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, integer A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer n

  integer a(n)
  integer i
  character ( len = * ) title

  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i8,i10)' ) i, a(i)
  end do

  return
end
subroutine i4vec_sort_heap_a ( n, a )

!*****************************************************************************80
!
!! i4vec_sort_heap_a() ascending sorts an integer array using heap sort.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input/output, integer A(N).
!    On input, the array to be sorted;
!    On output, the array has been sorted.
!
  implicit none

  integer n

  integer a(n)
  integer n1

  if ( n <= 1 ) then
    return
  end if
!
!  1: Put A into descending heap form.
!
  call i4vec_heap_d ( n, a )
!
!  2: Sort A.
!
!  The largest object in the heap is in A(1).
!  Move it to position A(N).
!
  call i4_swap ( a(1), a(n) )
!
!  Consider the diminished heap of size N1.
!
  do n1 = n - 1, 2, -1
!
!  Restore the heap structure of A(1) through A(N1).
!
    call i4vec_heap_d ( n1, a )
!
!  Take the largest object from A(1) and move it to A(N1).
!
    call i4_swap ( a(1), a(n1) )

  end do

  return
end
subroutine i4vec_uniq ( n, a, nuniq )

!*****************************************************************************80
!
!! i4vec_uniq() finds the number of unique elements in a sorted integer array.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements in A.
!
!    Input/output, integer A(N).  On input, the sorted
!    integer array.  On output, the unique elements in A.
!
!    Output, integer NUNIQ, the number of unique elements in A.
!
  implicit none

  integer n

  integer a(n)
  integer itest
  integer nuniq

  nuniq = 0

  if ( n <= 0 ) then
    return
  end if

  nuniq = 1
  
  do itest = 2, n

    if ( a(itest) /= a(nuniq) ) then
      nuniq = nuniq + 1
      a(nuniq) = a(itest)
    end if
 
  end do

  return
end
subroutine i4vec2_uniq ( n, a1, a2, nuniq )

!*****************************************************************************80
!
!! i4vec2_uniq() keeps the unique elements in a array of pairs of integers.
!
!  Discussion:
!
!    Item I is stored as the pair A1(I), A2(I).
!
!    The items must have been sorted, or at least it must be the
!    case that equal items are stored in adjacent vector locations.
!
!    If the items were not sorted, then this routine will only
!    replace a string of equal values by a single representative.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of items.
!
!    Input/output, integer A1(N), A2(N).
!    On input, the array of N items.
!    On output, an array of NUNIQ unique items.
!
!    Output, integer NUNIQ, the number of unique items.
!
  implicit none

  integer n

  integer a1(n)
  integer a2(n)
  integer itest
  integer nuniq

  nuniq = 0

  if ( n <= 0 ) then
    return
  end if

  nuniq = 1

  do itest = 2, n

    if ( a1(itest) /= a1(nuniq) .or. a2(itest) /= a2(nuniq) ) then

      nuniq = nuniq + 1

      a1(nuniq) = a1(itest)
      a2(nuniq) = a2(itest)

    end if

  end do

  return
end
subroutine network_flow_max ( nnode, nedge, iendpt, icpflo, isorce, isink, &
  icut, node_flow )

!*****************************************************************************80
!
!! network_flow_max() finds the maximal flow and a minimal cut in a network.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 July 2000
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, integer NEDGE, the number of edges.
!
!    Input/output, integer IENDPT(2,NEDGE), the edges of the
!    network, defined as pairs of nodes.  Each edge should be listed TWICE,
!    the second time in reverse order.  On output, the edges have
!    been reordered, and so the columns of IENDPT have been rearranged.
!
!    Input/output, integer ICPFLO(2,NEDGE).  Capacities and flows.
!    On input, ICPFLO(1,I) is the capacity of edge I.  On output,
!    ICPFLO(2,I) is the flow on edge I and ICPFLO(1,I) has
!    been rearranged to match the reordering of IENDPT.
!
!    Input, integer ISORCE, the designated source node.
!
!    Input, integer ISINK, the designated sink node.
!
!    Output, integer ICUT(NNODE).  ICUT(I) = 1 if node I is in the
!    minimal cut set, otherwise 0.
!
!    Output, integer NODE_FLOW(NNODE).  NODE_FLOW(I) is the value 
!    of the flow through node I.
!
  implicit none

  integer nedge
  integer nnode

  integer i
  integer iarray(nnode)
  integer icpflo(2,nedge)
  integer icut(nnode)
  integer idel
  integer ien1
  integer ien2
  integer iendpt(2,nedge)
  integer indx
  integer ip
  integer iparm
  integer iq
  integer ir
  integer iread
  integer irite
  integer is
  integer isink
  integer isorce
  integer isort
  integer it
  integer itemp
  integer iwork(nnode,2)
  integer j
  integer kz
  integer lst
  integer m
  integer node_flow(nnode)

  iarray(1:nnode) = 0
  idel = 0
 
  do i = 1, nedge

    icpflo(2,i) = 0
    ip = iendpt(1,i)

    if ( ip == isorce ) then
      idel = idel + icpflo(1,i)
    end if

    iarray(ip) = iarray(ip) + 1

  end do
 
  node_flow(isorce) = idel
  is = 1
 
  do i = 1, nnode
    it = iarray(i)
    iarray(i) = is
    iwork(i,1) = is
    is = is + it
  end do
 
  isort = 0
  ien1 = 0
  ien2 = 0
 
10    continue
 
  indx = 0
 
50    continue
 
  call sort_heap_external ( nedge, indx, ien1, ien2, is )

  if ( indx < 0 ) then
 
    is = iendpt(1,ien1) - iendpt(1,ien2)

    if ( is == 0 ) then
      is = iendpt(2,ien1) - iendpt(2,ien2)
    end if

  else if ( 0 < indx ) then
 
    do ir = 1, 2
      itemp           = iendpt(ir,ien1)
      iendpt(ir,ien1) = iendpt(ir,ien2)
      iendpt(ir,ien2) = itemp

      itemp           = icpflo(ir,ien1)
      icpflo(ir,ien1) = icpflo(ir,ien2)
      icpflo(ir,ien2) = itemp
    end do
 
  else
 
    if ( 0 < isort ) then
      return
    end if
 
    do i = 1, nedge
      iq = iendpt(2,i)
      iendpt(1,i) = iwork(iq,1)
      iwork(iq,1) = iwork(iq,1) + 1
    end do
 
    go to 100
 
  end if
 
  go to 50
 
80    continue
 
  iendpt(1,iendpt(1,ien1)) = ien2
  iendpt(1,iendpt(1,ien2)) = ien1
 
  do ir = 1, 2
    itemp           = iendpt(ir,ien1)
    iendpt(ir,ien1) = iendpt(ir,ien2)
    iendpt(ir,ien2) = itemp

    itemp           = icpflo(ir,ien1)
    icpflo(ir,ien1) = icpflo(ir,ien2)
    icpflo(ir,ien2) = itemp
  end do
 
  if ( indx < 0 ) then
    go to 270
  end if

  if ( indx == 0 ) then
    go to 170
  end if

  go to 50
 
100   continue
 
  indx = 0
 
  do i = 1, nnode

    if ( i /= isorce ) then
      node_flow(i) = 0
    end if

    iwork(i,2) = nedge + 1

    if ( i < nnode ) then
      iwork(i,2) = iarray(i+1)
    end if

    icut(i) = 0

  end do
 
  iread = 0
  irite = 1
  iwork(1,1) = isorce
  icut(isorce) = - 1
 
120   continue
 
  iread = iread + 1
 
  if ( iread <= irite ) then
 
    ip = iwork(iread,1)
    lst = iwork(ip,2) - 1
    i = iarray(ip) - 1
 
130     continue
 
    i = i + 1

    if ( lst < i ) then
      go to 120
    end if

    iq = iendpt(2,i)
    idel = icpflo(1,i) - icpflo(2,i)

    if ( icut(iq) /= 0 .or. idel == 0 ) then
      go to 130
    end if
 
    if ( iq /= isink ) then
      irite = irite + 1
      iwork(irite,1) = iq
    end if
 
    icut(iq) = - 1
    go to 130
 
  end if
 
  if ( icut(isink) == 0 ) then
 
    icut(1:nnode) = - icut(1:nnode)
 
    do i = 1, nedge
      ip = iendpt(2,iendpt(1,i))
      if ( icpflo(2,i) < 0 ) then
        node_flow(ip) = node_flow(ip) - icpflo(2,i)
      end if
      iendpt(1,i) = ip
    end do
 
    node_flow(isorce) = node_flow(isink)
    isort = 1
    go to 10
 
  end if
 
  icut(isink) = 1
 
160   continue
 
  iread = iread - 1

  if ( iread == 0 ) then
    go to 180
  end if

  ip = iwork(iread,1)
  ien1 = iarray(ip) - 1
  ien2 = iwork(ip,2) - 1
 
170   continue
 
  if ( ien1 /= ien2 ) then
 
    iq = iendpt(2,ien2)
 
    if ( icut(iq) <= 0 .or. icpflo(1,ien2) == icpflo(2,ien2) ) then
      ien2 = ien2 - 1
      go to 170
    end if
 
    iendpt(2,ien2) = - iq
    icpflo(1,ien2) = icpflo(1,ien2) - icpflo(2,ien2)
    icpflo(2,ien2) = 0
    ien1 = ien1 + 1

    if ( ien1 < ien2 ) then
      go to 80
    end if
 
  end if
 
  if ( iarray(ip) <= ien1 ) then
    icut(ip) = ien1
  end if

  go to 160
 
180   continue
 
  kz = 0
 
  do ir = 1, irite
    if ( 0 < icut(iwork(ir,1)) ) then
      kz = kz + 1
      iwork(kz,1) = iwork(ir,1)
    end if
  end do
 
  indx = - 1
  m = 1
 
200   continue
 
  ip = iwork(m,1)

  if ( 0 < node_flow(ip) ) then
    go to 250
  end if
 
210   continue
 
  m = m + 1

  if ( m <= kz ) then
    go to 200
  end if

  iparm = 0
 
220   continue
 
  m = m - 1
 
  if ( m == 1 ) then
 
    do i = 1, nedge
 
      iq = - iendpt(2,i)
 
      if ( 0 <= iq ) then

        iendpt(2,i) = iq
        j = iendpt(1,i)
        icpflo(1,i) = icpflo(1,i) - icpflo(2,j)

        idel = icpflo(2,i) - icpflo(2,j)
        icpflo(2,i) = idel
        icpflo(2,j) = - idel

      end if
 
    end do
 
    go to 100
 
  end if
 
  ip = iwork(m,1)

  if ( node_flow(ip) < 0 ) then
    go to 220
  end if
 
  if ( node_flow(ip) == 0 ) then

    lst = nedge + 1

    if ( ip < nnode ) then
      lst = iarray(ip+1)
    end if

    i = iwork(ip,2)
    iwork(ip,2) = lst
 
240     continue
 
    if ( i == lst ) then
      go to 220
    end if

    j = iendpt(1,i)
    idel = icpflo(2,j)
    icpflo(2,j) = 0
    icpflo(1,j) = icpflo(1,j) - idel
    icpflo(2,i) = icpflo(2,i) - idel
    i = i + 1
    go to 240
 
  end if
 
  if ( icut(ip) < iarray(ip) ) then
    go to 300
  end if
 
250   continue
 
  i = icut(ip) + 1
 
260   continue

  i = i - 1

  if ( i < iarray(ip) ) then
    go to 290
  end if

  iq = - iendpt(2,i)

  if ( node_flow(iq) < 0 ) then
    go to 260
  end if
 
  idel = icpflo(1,i) - icpflo(2,i)

  if ( node_flow(ip) < idel ) then
    idel = node_flow(ip)
  end if

  icpflo(2,i) = icpflo(2,i) + idel
  node_flow(ip) = node_flow(ip) - idel
  node_flow(iq) = node_flow(iq) + idel
  iparm = 1
  ien1 = iendpt(1,i)
  ien2 = iwork(iq,2) - 1

  if ( ien1 < ien2 ) then
    go to 80
  end if

  if ( ien1 /= ien2 ) then
    go to 280
  end if
 
270   continue
 
  iwork(iq,2) = ien2
 
280   continue
 
  if ( 0 < node_flow(ip) ) then
    go to 260
  end if

  if ( icpflo(1,i) == icpflo(2,i) ) then
    i = i - 1
  end if
 
290   continue
 
  icut(ip) = i

  if ( iparm /= 0 ) then
    go to 210
  end if
 
300   continue
 
  i = iwork(ip,2)
 
310   continue
 
  j = iendpt(1,i)
  idel = icpflo(2,j)

  if ( node_flow(ip) < idel ) then
    idel = node_flow(ip)
  end if

  icpflo(2,j) = icpflo(2,j) - idel
  node_flow(ip) = node_flow(ip) - idel
  iq = iendpt(2,i)
  node_flow(iq) = node_flow(iq) + idel
  i = i + 1

  if ( 0 < node_flow(ip) ) then
    go to 310
  end if

  node_flow(ip) = - 1
  go to 220
 
end
subroutine nodes_to_ps ( plotxmin2, plotymin2, alpha, iunit, nnode, x, y, &
  xmin, ymin )

!*****************************************************************************80
!
!! nodes_to_ps() writes subplot nodes to a PostScript file.
!
!  Discussion:
!
!    A small filled circle is placed at each node.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer PLOTXMIN2, PLOTYMIN2, the Postscript point
!    corresponding to the physical point XMIN, YMIN.
!
!    Input, real ( kind = rk ) ALPHA, the physical-to-Postscript scale factor.
!
!    Input, integer IUNIT, the output FORTRAN unit.
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, real ( kind = rk ) X(NNODE), Y(NNODE), the coordinates of points.
!
!    Input, real ( kind = rk ) XMIN, YMIN, the coordinates of the physical
!    origin.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nnode

  real ( kind = rk ) alpha
  integer i
  integer iunit
  integer plotxmin2
  integer plotymin2
  integer px1
  integer py1
  integer, parameter :: rad = 10
  real ( kind = rk ) x(nnode)
  real ( kind = rk ) xmin
  real ( kind = rk ) y(nnode)
  real ( kind = rk ) ymin
!
!  Draw points.
!
  do i = 1, nnode

    px1 = plotxmin2 + nint ( alpha * ( x(i) - xmin ) )
    py1 = plotymin2 + nint ( alpha * ( y(i) - ymin ) )

    write ( iunit, '(3i4,a)' ) px1, py1, rad, ' 0 360 arc closepath fill'

  end do

  return
end
subroutine poly ( n, iarray, ix0, iopt, ival )

!*****************************************************************************80
!
!! poly() performs operations on polynomials in power or factorial form.
!
!  Definition:
!
!    The power sum form of a polynomial is
!
!      P(X) = A1+A2*X+A3*X^2+...+(AN+1)*X^N
!
!    The Taylor expansion at C has the form
!
!      P(X) = A1+A2*(X-C)+A3*(X-C)^2+...+(AN+1)*(X-C)^N
!
!    The factorial form of a polynomial is
!
!      P(X) = A1+A2*X+A3*(X)*(X-1)+A4*(X)*(X-1)*(X-2)+...
!        +(AN+1)*(X)*(X-1)*...*(X-N+1)
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 April 1999
!
!  Parameters:
!
!    Input, integer N, the number of coefficients in the polynomial
!    (in other words, the polynomial degree + 1)
!
!    Input, integer IOPT, a flag describing which algorithm is to
!    be carried out:
!
!    -3: Reverse Stirling.  Input the coefficients of
!    the polynomial in factorial form, output them in
!    power sum form.
!
!    -2: Stirling.  Input the coefficients in power sum
!    form, output them in factorial form.
!
!    -1: Evaluate a polynomial which has been input
!    in factorial form.
!
!    0:  Evaluate a polynomial input in power sum form.
!
!    1 or more:  Given the coefficients of a polynomial in
!    power sum form, compute the first IOPT coefficients of
!    the polynomial in Taylor expansion form.
!
!    Input, integer IX0, for IOPT = -1, 0, or positive, the value
!    X of the argument at which the polynomial is to be evaluated, or the
!    Taylor expansion is to be carried out.
!
!    Output, integer IVAL, for IOPT = -1 or 0, the value of the
!    polynomial at the point IX0.
!
!    Input, integer IARRAY(N).  Contains the coefficients of the
!    polynomial.  Depending on the option chosen, these coefficients may
!    be overwritten by those of a different form of the polynomial.
!
  implicit none

  integer n

  integer i
  integer iarray(n)
  integer ieps
  integer iopt
  integer ival
  integer iw
  integer ix0
  integer iz
  integer m
  integer n1

  n1 = min ( n, iopt )
  n1 = max ( 1, n1 )
 
  if ( iopt < -1 ) then
    n1 = n
  end if
 
  ieps = mod ( max ( -iopt, 0 ), 2 )
 
  iw = -n * ieps

  if ( -2 < iopt ) then
    iw = iw + ix0
  end if
 
  do m = 1, n1
 
    ival = 0
    iz = iw
 
    do i = m, n
      iz = iz + ieps
      ival = iarray(n+m-i) + iz * ival
      if ( iopt /= 0 .and. iopt /= -1 ) then
        iarray(n+m-i) = ival
      end if
    end do
 
    if ( iopt < 0 ) then
      iw = iw + 1
    end if
 
  end do
 
  return
end
subroutine r8mat_print ( a, ihi, ilo, jhi, jlo, ncol, nrow )

!*****************************************************************************80
!
!! r8mat_print() prints out a portion of a dense matrix.
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
!    Input, real ( kind = rk ) A(NROW,NCOL), an NROW by NCOL matrix to be printed.
!
!    Input, integer IHI, ILO, the first and last rows to print.
!
!    Input, integer JHI, JLO, the first and last columns to print.
!
!    Input, integer NCOL, NROW, the number of rows and columns
!    in the matrix A.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: incx = 5

  integer ncol

  real ( kind = rk ) a(nrow,ncol)
  character ctemp(incx)*14
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
  integer nrow

  write ( *, '(a)' ) ' '

  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, ncol )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(''Columns '',5a14)' ) ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, nrow )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == 0.0D+00 ) then
          ctemp(j2) = '    0.0'
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ctemp(1:inc)

    end do

  end do

  return
end
subroutine sort_heap_external ( n, indx, i, j, isgn )

!*****************************************************************************80
!
!! sort_heap_external() externally sorts a list of items into ascending order.
!
!  Discussion:
!
!    The actual list of data is not passed to the routine.  Hence this
!    routine may be used to sort integers, real ( kind = rk )s, numbers, names,
!    dates, shoe sizes, and so on.  After each call, the routine asks
!    the user to compare or interchange two items, until a special
!    return value signals that the sorting is completed.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 February 2004
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer N, the number of items to be sorted.
!
!    Input/output, integer INDX, the main communication signal.
!
!    The user must set INDX to 0 before the first call.
!    Thereafter, the user should not change the value of INDX until
!    the sorting is done.
!
!    On return, if INDX is
!
!      greater than 0,
!      * interchange items I and J;
!      * call again.
!
!      less than 0,
!      * compare items I and J;
!      * set ISGN = -1 if I < J, ISGN = +1 if J < I;
!      * call again.
!
!      equal to 0, the sorting is done.
!
!    Output, integer I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ISGN, results of comparison of elements I 
!    and J.
!    (Used only when the previous call returned INDX less than 0).
!    ISGN <= 0 means I is less than or equal to J;
!    0 <= ISGN means I is greater than or equal to J.
!
  implicit none

  integer i
  integer, save :: i_save = 0
  integer indx
  integer isgn
  integer j
  integer, save :: j_save = 0
  integer, save :: k = 0
  integer, save :: k1 = 0
  integer n
  integer, save :: n1 = 0
!
!  INDX = 0: This is the first call.
!
  if ( indx == 0 ) then

    i_save = 0
    j_save = 0
    k = n / 2
    k1 = k
    n1 = n
!
!  INDX < 0: The user is returning the results of a comparison.
!
  else if ( indx < 0 ) then

    if ( indx == -2 ) then

      if ( isgn < 0 ) then
        i_save = i_save + 1
      end if

      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return

    end if

    if ( 0 < isgn ) then
      indx = 2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then

      if ( n1 == 1 ) then
        i_save = 0
        j_save = 0
        indx = 0
      else
        i_save = n1
        n1 = n1 - 1
        j_save = 1
        indx = 1
      end if

      i = i_save
      j = j_save
      return

    end if

    k = k - 1
    k1 = k
!
!  0 < INDX, the user was asked to make an interchange.
!
  else if ( indx == 1 ) then

    k1 = k

  end if

  do

    i_save = 2 * k1

    if ( i_save == n1 ) then
      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return
    else if ( i_save <= n1 ) then
      j_save = i_save + 1
      indx = -2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then
      exit
    end if

    k = k - 1
    k1 = k

  end do

  if ( n1 == 1 ) then
    i_save = 0
    j_save = 0
    indx = 0
    i = i_save
    j = j_save
  else
    i_save = n1
    n1 = n1 - 1
    j_save = 1
    indx = 1
    i = i_save
    j = j_save
  end if

  return
end
subroutine span_forest ( nnode, nedge, iendpt, k, component )

!*****************************************************************************80
!
!! span_forest() determines a graph's connectivity and spanning forest.
!
!  Discussion:
!
!    The input graph may be connected or unconnected.
!
!    If the input graph is connected, this routine simply returns a
!    spanning tree for the graph.
!
!    Definition: A (connected) component of a graph is a maximal subgraph
!    which is connected.
!
!    Definition: A tree is a connected graph containing no cycles.
!
!    Definition: A spanning tree of a connected graph is a subgraph which 
!    is a maximal tree.
!
!    Definition: A forest is a collection of trees, no two of which share 
!    a node.
!
!    Definition: A spanning forest of a possibly unconnected graph 
!    is a collection containing a single spanning tree for each component 
!    of the graph.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 May 1999
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, integer NEDGE, the number of edges in graph.
!
!    Input/output, integer IENDPT(2,NEDGE), the edge array of 
!    the graph.  IENDPT(1,I) and IENDPT(2,I) are the two nodes that make up 
!    edge I.
!  
!    On input, IENDPT describes the graph.
!
!    On output, the input entries of IENDPT have been reordered, so that
!    edges belonging to the spanning forest come first, followed by those
!    edges which are not part of the spanning forest.
!
!    Output, integer K, the number of connected components of the 
!    graph.
!
!    Output, integer COMPONENT(NNODE), the component to which each 
!    node belongs.
!
  implicit none

  integer nedge
  integer nnode

  integer component(nnode)
  integer i
  integer iendpt(2,nedge)
  integer ip
  integer iq
  integer ir
  integer iret
  integer is
  integer itemp
  integer k
  integer l
  integer l0
  integer l1
  integer loc
  integer m
  integer m0
  integer m1
  integer mm

  mm = 1 + max ( nnode, nedge )
 
  do i = 1, nnode
    component(i) = -i
  end do
 
  do m = 1, nedge
    do l = 1, 2
      ip = iendpt(l,m)
      iendpt(l,m) = component(ip)
      component(ip) = - l * mm - m
    end do
  end do
 
  k = 0
  loc = 0
 
10 continue
 
  do i = 1, nnode
 
    iq = component(i)
 
    if ( iq <= 0 ) then
 
      k = k + 1
      component(i) = k
 
      if ( iq + i < 0 ) then
        ip = i
        is = - iq
        iret = 31
        l1 = - iq / mm
        m1 = - iq - l1 * mm
        go to 110
      end if
 
    end if
 
  end do
 
  do m = 1, nedge
 
    do
 
      ir = - iendpt(1,m)
 
      if ( ir < 0 ) then
        exit
      end if

      itemp = iendpt(2,m)
      iendpt(2,m) = iendpt(2,ir)
      iendpt(2,ir) = itemp

      iendpt(1,m) = iendpt(1,ir)
      iendpt(1,ir) = component(iendpt(2,ir))

    end do
 
  end do
 
  component(iendpt(2,1:loc)) = component(iendpt(1,1:loc))
 
  return
 
90    continue
 
  if ( ir /= 0 ) then
    loc = loc + 1
    component(ip) = iendpt(1,ir) + iendpt(2,ir) - ip
    iendpt(1,ir) = - loc
    iendpt(2,ir) = ip
  end if
 
  ip = m

  if ( m <= 0 ) then
    go to 10
  end if

  is = - component(ip)
 
100   continue
 
  l = is / mm
  m = is - l * mm

  if ( l == 0 ) then
    go to 90
  end if

  l1 = 3 - l
  m1 = m
 
110   continue
 
  iq = iendpt(l1,m1)
 
  if ( 0 < iq ) then
    if ( iq <= mm ) then
      if ( 0 <= component(iq) ) then
        ir = m1
      end if
    end if
  end if
 
  if ( 0 <= iq ) then
    is = abs ( iendpt(l,m) )
    iendpt(l,m) = ip
    go to 100
  end if
 
  if ( -mm <= iq ) then
 
    iq = - iq
    iendpt(l1,m1) = 0
 
    if ( iret == 31 ) then

      l0 = l1
      m0 = m1
      ir = 0
      iret = 43

    else
 
      iendpt(l0,m0) = iq
      l0 = l1
      m0 = m1
      is = abs ( iendpt(l,m) )
      iendpt(l,m) = ip

    end if

    go to 100
 
  end if
 
  iendpt(l1,m1) = - iq
  l1 = - iq / mm
  m1 = - iq - l1 * mm
  go to 110
 
end

