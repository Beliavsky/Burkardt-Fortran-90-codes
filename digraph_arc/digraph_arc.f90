subroutine digraph_adj_print ( adj, nnode, title )

!*****************************************************************************80
!
!! digraph_adj_print() prints out an adjacency matrix for a digraph.
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
!    Input, integer ADJ(NNODE,NNODE), the adjacency matrix of a 
!    digraph.  ADJ(I,J) is 1 if there is a direct connection FROM node I TO 
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

    jhi = min ( nnode, 80 )

    do j = 1, jhi

      if ( 0 <= adj(i,j) .and. adj(i,j) <= 9 ) then
        string(j:j) = char ( 48 + adj(i,j) )
      else
        string(j:j) = '*'
      end if

    end do

    write ( *, '(i2,2x,a)' ) i, string(1:jhi)

  end do

  return
end
subroutine digraph_arc_degree ( nnode, nedge, inode, jnode, indegree, &
  outdegree )

!*****************************************************************************80
!
!! digraph_arc_degree() determines the degree of the nodes of a digraph.
!
!  Discussion:
!
!    Definition: The degree of a node is the number of edges that 
!    include the node.
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
!    integer NNODE, the number of nodes.
!
!    integer NEDGE, the number of edges.
!
!    integer INODE(NEDGE), JNODE(NEDGE), pairs of nodes forming edges.
!
!  Output:
!
!    integer INDEGREE(NNODE), OUTDEGREE(NNODE), the
!    indegree and outdegree of each node, that is, the number of edges that end 
!    with the node, and the number of edges that begin with it.
!
  implicit none

  integer nedge
  integer nnode

  integer i
  integer indegree(nnode)
  integer inode(nedge)
  integer jnode(nedge)
  integer n
  integer outdegree(nnode)

  indegree(1:nnode) = 0
  outdegree(1:nnode) = 0

  do i = 1, nedge

    n = inode(i)
    if ( 1 <= n .and. n <= nnode ) then
      outdegree(n) = outdegree(n) + 1
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DIGRAPH_ARC_DEGREE(): Fatal error!'
      write ( *, '(a,i8)' ) '  Out-of-range node value = ', n
      stop 1
    end if

    n = jnode(i)
    if ( 1 <= n .and. n <= nnode ) then
      indegree(n) = indegree(n) + 1
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DIGRAPH_ARC_DEGREE(): Fatal error!'
      write ( *, '(a,i8)' ) '  Out-of-range node value = ', n
      stop 1
    end if

  end do

  return
end
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
subroutine digraph_arc_euler_circ_cand ( nedge, inode, jnode, circuit, k, &
  nstack, stack, maxstack, ncan, iwork, lwork )

!*****************************************************************************80
!
!! digraph_arc_euler_circ_cand(): candidates for K-th edge of an Euler circuit.
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
!    17 August 2000
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
!    Input, integer NEDGE, the number of edges in the digraph.
!
!    Input, integer INODE(NEDGE), JNODE(NEDGE), the edge array of 
!    the digraph.  The I-th edge extends from node INODE(I) to JNODE(I).
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
!    Input, integer STACK(MAXSTACK), as yet unused candidates for 
!    positions 1 to K-1.
!
!    Input, integer MAXSTACK, the dimension of STACK.
!
!    Workspace, integer IWORK(NEDGE).
!
!    Workspace, logical LWORK(NEDGE).
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
    lwork(i) = it == inode(i)
  end do
 
  lwork(circuit(1:k-1)) = .false.
  
  do i = 1, nedge
    if ( lwork(i) ) then
      if ( maxstack <= nstack ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DIGRAPH_ARC_EULER_CIRC_CAND(): Fatal error!'
        write ( *, '(a)' ) '  Stack size exceeded.'
        stop 1
      end if
      nstack = nstack + 1
      stack(nstack) = i
      ncan(k) = ncan(k) + 1
    end if
  end do
 
  return
end
subroutine digraph_arc_euler_circ_next ( nedge, inode, jnode, circuit, stack, &
  maxstack, ncan, more )

!*****************************************************************************80
!
!! digraph_arc_euler_circ_next() returns the next Euler circuit for a digraph.
!
!  Discussion:
!
!    The routine produces all the Euler circuits of a digraph, one at a time.
!
!    Definition: An Euler circuit of a digraph is a path starting at some node, 
!    using all the edges of the digraph exactly once, and returning
!    to the starting node.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 August 2000
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
!    Input, integer NEDGE, the number of edges in the digraph.
!
!    Input, integer INODE(NEDGE), JNODE(NEDGE), the edge array 
!    of the digraph.  The I-th edge extends from node INODE(I) to JNODE(I).
!
!    Output, integer CIRCUIT(NEDGE).  If MORE = TRUE on output, 
!    then IARRAY contains the edges, in order, that constitute this circuit.
!
!    Workspace, integer STACK(MAXSTACK).  
!
!    Input, integer MAXSTACK, the dimension of STACK.
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
  logical lwork(nedge)
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

      call digraph_arc_euler_circ_cand ( nedge, inode, jnode, circuit, k, &
        nstack, stack, maxstack, ncan, iwork, lwork )

    else

      more = .false.
      exit

    end if

  end do
 
  return
end
subroutine digraph_arc_example_cycler ( inode, jnode )

!*****************************************************************************80
!
!! digraph_arc_example_cycler() sets arc list information for the cycler digraph.
!
!  Diagram:
!  
!           A
!           |
!           V
!    9--><--7---<--3--><---4
!    |            /|      /
!    V           A |     /
!    |          /  |    /
!    5----<----1   V   A
!    |        /    |  /
!    V       A     | /
!    |      /      |/
!    2-->---8---<--6
!     \------>----/
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
!  Output:
!
!    integer INODE(MAXEDGE), JNODE(MAXEDGE), the arc list
!    for the digraph.
!
  implicit none

  integer, parameter :: nedge = 16

  integer inode(nedge)
  integer jnode(nedge)

  inode(1) = 1
  jnode(1) = 3

  inode(2) = 1
  jnode(2) = 5

  inode(3) = 2
  jnode(3) = 6

  inode(4) = 2
  jnode(4) = 8

  inode(5) = 3
  jnode(5) = 4

  inode(6) = 3
  jnode(6) = 6

  inode(7) = 3
  jnode(7) = 7

  inode(8) = 4
  jnode(8) = 3

  inode(9) = 5
  jnode(9) = 2

  inode(10) = 6
  jnode(10) = 4

  inode(11) = 6
  jnode(11) = 8

  inode(12) = 7
  jnode(12) = 7

  inode(13) = 7
  jnode(13) = 9

  inode(14) = 8
  jnode(14) = 1

  inode(15) = 9
  jnode(15) = 5

  inode(16) = 9
  jnode(16) = 7

  return
end
subroutine digraph_arc_is_eulerian ( nnode, nedge, inode, jnode, indegree, &
  outdegree, result )

!*****************************************************************************80
!
!! digraph_arc_is_eulerian() determines if a digraph is Eulerian.
!
!  Discussion:
!
!    A digraph is Eulerian if there exists a circuit through the graph
!    which uses every edge once.
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
!    Input, integer NNODE, the number of nodes.
!
!    Input, integer NEDGE, the number of edges.
!
!    Input, integer INODE(NEDGE), JNODE(NEDGE), the pairs of nodes
!    that form the edges.
!
!    Output, integer INDEGREE(NNODE), OUTDEGREE(NODE), the
!    indegree and outdegree of each node, that is, the number of edges that 
!    end with the node, and that begin the node.
!
!    Output, integer RESULT.
!    0, the digraph is not Eulerian.
!    1, the digraph is Eulerian, but the starting and ending nodes differ.
!    2, the digraph is Eulerian, and there is a closed Euler circuit.
!
  implicit none

  integer nedge
  integer nnode

  integer i
  integer indegree(nnode)
  integer inode(nedge)
  integer jnode(nedge)
  integer n_minus
  integer n_plus
  integer outdegree(nnode)
  integer result

  call digraph_arc_degree ( nnode, nedge, inode, jnode, indegree, outdegree )

  n_plus = 0
  n_minus = 0

  do i = 1, nnode

    if ( indegree(i) == outdegree(i) ) then

    else if ( n_plus == 0 .and. indegree(i) == outdegree(i) + 1 ) then
      n_plus = 1
    else if ( n_minus == 0 .and. indegree(i) == outdegree(i) - 1 ) then
      n_minus = 1
    else
      result = 0
      return
    end if

  end do

  if ( n_plus == 0 .and. n_minus == 0 ) then
    result = 2
  else if ( n_plus == 1 .and. n_minus == 1 ) then
    result = 1
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIGRAPH_ARC_IS_EULERIAN(): Fatal error!'
    write ( *, '(a)' ) '  The algorithm failed.'
    stop 1
  end if

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
subroutine digraph_arc_to_digraph_adj ( nedge, inode, jnode, adj, nnode )

!*****************************************************************************80
!
!! digraph_arc_to_digraph_adj() converts arc list digraph to an adjacency digraph.
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
!    Input, integer NEDGE, the number of edges.
!
!    Input, integer INODE(NEDGE), JNODE(NEDGE), the edge array.
!    The I-th edge connects nodes INODE(I) and JNODE(I).
!
!    Output, integer ADJ(NNODE,NNODE), the adjacency information.
!
  implicit none

  integer nedge
  integer nnode

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
  end do

  return
end
subroutine digraph_arc_to_digraph_star ( nnode, nedge, inode, jnode, arcfir, &
  fwdarc )

!*****************************************************************************80
!
!! digraph_arc_to_digraph_star() sets forward star representation of a digraph.
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
!    Input, integer INODE(NEDGE), JNODE(NEDGE); the I-th edge
!    extends from node INODE(I) to JNODE(I).
!
!    Output, integer ARCFIR(NNODE+1); ARCFIR(I) is the number of
!    the first edge starting at node I in the forward star representation.
!
!    Output, integer FWDARC(NEDGE); FWDARC(I) is the ending node of
!    the I-th edge in the forward star representation.
!
  implicit none

  integer nedge
  integer nnode

  integer arcfir(nnode+1)
  integer fwdarc(nedge)
  integer i
  integer inode(nedge)
  integer j
  integer jnode(nedge)
  integer k
!
!  Set up the forward star representation.
!
  k = 0

  do i = 1, nnode

    arcfir(i) = k + 1

    do j = 1, nedge

      if ( inode(j) == i ) then
        k = k + 1
        fwdarc(k) = jnode(j)
      end if

    end do

  end do

  arcfir(nnode+1) = k + 1

  return
end
subroutine digraph_arc_weight_print ( nedge, inode, jnode, wnode, title )

!*****************************************************************************80
!
!! digraph_arc_weight_print() prints out a weighted digraph from an edge list.
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
!    Input, integer INODE(NEDGE), JNODE(NEDGE), the beginning and
!    end nodes of the edges.
!
!    Input, real ( kind = rk ) WNODE(NEDGE), the weights of the edges.
!
!    Input, character ( len = * ) TITLE, a title.
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
subroutine i4vec2_print ( n, a, b, title )

!*****************************************************************************80
!
!! i4vec2_print() prints a pair of integer vectors.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, integer A(N), B(N), the vectors to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer n

  integer a(n)
  integer b(n)
  integer i
  character ( len = * ) title

  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i8,2i10)' ) i, a(i), b(i)
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

