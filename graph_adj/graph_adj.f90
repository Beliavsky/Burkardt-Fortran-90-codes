subroutine dge_check ( m, n, ierror )

!*****************************************************************************80
!
!! dge_check() checks the dimensions of a general matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
!
!    Output, integer IERROR, reports whether any errors 
!    were detected.
!    IERROR is set to 0 before the checks are made, and then:
!    IERROR = IERROR + 2 if M is illegal;
!    IERROR = IERROR + 4 if N is illegal.
!
  implicit none

  integer ierror
  integer m
  integer n

  ierror = 0

  if ( m < 1 ) then
    ierror = ierror + 2
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) 'DGE_CHECK - Illegal M = ', m
  end if

  if ( n < 1 ) then
    ierror = ierror + 4
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) 'DGE_CHECK - Illegal N = ', n
  end if

  return
end
subroutine dge_det ( n, a, ipivot, det )

!*****************************************************************************80
!
!! dge_det() computes the determinant of a matrix factored by DGE_FA or DGE_TRF.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real ( kind = rk ) A(N,N), the LU factors computed 
!    by DGE_FA or DGE_TRF.
!
!    Input, integer IPIVOT(N), as computed by DGE_FA or DGE_TRF.
!
!    Output, real ( kind = rk ) DET, the determinant of the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) det
  integer i
  integer ierror
  integer ipivot(n)
!
!  Check the dimensions.
!
  call dge_check ( n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DGE_DET(): Fatal error!'
    write ( *, '(a)' ) '  Illegal dimensions.'
    return
  end if

  det = 1.0D+00

  do i = 1, n
    det = det * a(i,i)
  end do

  do i = 1, n
    if ( ipivot(i) /= i ) then
      det = - det
    end if
  end do

  return
end
subroutine dge_fa ( n, a, ipivot, info )

!*****************************************************************************80
!
!! dge_fa() factors a general matrix.
!
!  Discussion:
!
!    DGE_FA is a simplified version of the LINPACK routine DGEFA.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input/output, real ( kind = rk ) A(N,N), the matrix to be factored.
!    On output, A contains an upper triangular matrix and the multipliers
!    which were used to obtain it.  The factorization can be written
!    A = L * U, where L is a product of permutation and unit lower
!    triangular matrices and U is upper triangular.
!
!    Output, integer IPIVOT(N), a vector of pivot indices.
!
!    Output, integer INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,n)
  integer i
  integer ierror
  integer info
  integer ipivot(n)
  integer j
  integer k
  integer l
!
!  Check the dimensions.
!
  call dge_check ( n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DGE_FA(): Fatal error!'
    write ( *, '(a)' ) '  Illegal dimensions.'
    return
  end if

  info = 0

  do k = 1, n-1
!
!  Find L, the index of the pivot row.
!
    l = k
    do i = k+1, n
      if ( abs ( a(l,k) ) < abs ( a(i,k) ) ) then
        l = i
      end if
    end do

    ipivot(k) = l
!
!  If the pivot index is zero, the algorithm has failed.
!
    if ( a(l,k) == 0.0D+00 ) then
      info = k
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DGE_FA(): Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      return
    end if
!
!  Interchange rows L and K if necessary.
!
    if ( l /= k ) then
      call r8_swap ( a(l,k), a(k,k) )
    end if
!
!  Normalize the values that lie below the pivot entry A(K,K).
!
    a(k+1:n,k) = -a(k+1:n,k) / a(k,k)
!
!  Row elimination with column indexing.
!
    do j = k + 1, n

      if ( l /= k ) then
        call r8_swap ( a(l,j), a(k,j) )
      end if

      a(k+1:n,j) = a(k+1:n,j) + a(k+1:n,k) * a(k,j)

    end do

  end do

  ipivot(n) = n

  if ( a(n,n) == 0.0D+00 ) then
    info = n
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DGE_FA(): Fatal error!'
    write ( *, '(a,i8)' ) '  Zero pivot on step ', info
  end if

  return
end
subroutine graph_adj_breadth_first ( adj, nnode, dad, deep, order )

!*****************************************************************************80
!
!! graph_adj_breadth_first() carries out a breadth-first traversal of a graph.
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
!    Alan Gibbons,
!    Algorithmic Graph Theory,
!    Cambridge University Press, 1985,
!    ISBN 0-521-28881-9.
!
!  Parameters:
!
!    Input, integer adj(nnode,nnode), is the adjacency information.
!    ADJ(I,J) is nonzero if there is an edge from node
!    I to node J, and 0 otherwise.
!
!    Input, integer NNODE, the number of nodes.
!
!    Output, integer DAD(NNODE), DAD(I) is the node from which
!    node I is visited.  Node 1 is the first node in the search,
!    and has no predecessor, so DAD(1) is zero.  If there is
!    more than one connected component, then there
!    will be other nodes with DAD equal to zero.
!
!    Output, integer DEEP(NNODE), records the "depth" of the node.
!    The first node, node 1, has depth 1.  All the nodes that
!    can be reached in one step from node 1 have depth 2.  All
!    nodes that can be reached in one step from any of those nodes
!    have depth 3.  If there is more than one connected component,
!    then the depth of nodes in the second component will begin
!    one greater than the greatest depth of the first component,
!    and so on.
!
!    Output, integer ORDER(NNODE).  ORDER(I) is the step at which
!    node I is visited in the search.
!
  implicit none

  integer nnode

  integer adj(nnode,nnode)
  integer dad(nnode)
  integer deep(nnode)
  integer i
  integer order(nnode)
  integer iput
  integer queue(nnode)
  integer itake
  integer j
  integer jdeep
  integer k
  integer nudeep

  deep(1:nnode) = 0
  order(1:nnode) = 0
  dad(1:nnode) = 0
 
  k = 0
  i = 1
  iput = 1
  itake = 1
  nudeep = iput
  queue(iput) = i
  jdeep = 1
  deep(i) = jdeep
  k = k + 1
  order(i) = k
  dad(i) = 0
!
!  Find all sons of this father.
!  Store all sons in the son stack.
!
10    continue
 
  do j = 1, nnode
 
    if ( ( adj(i,j) /= 0 .or. adj(j,i) /= 0 ) .and. order(j) == 0 ) then

      iput = iput + 1

      if ( nnode < iput ) then
        iput = 1
      end if

      queue(iput) = j
      k = k + 1
      dad(j) = i
      order(j) = k
      deep(j) = jdeep + 1

    end if
 
  end do
!
!  Are there more fathers whose sons are to be searched for?
!
  if ( iput /= itake ) then
 
    if ( itake == nudeep ) then
      jdeep = jdeep + 1
      nudeep = iput
    end if
 
    i = queue(itake)
    itake = itake + 1

    if ( nnode < itake ) then
      itake = 1
    end if

    go to 10
!
!  No more fathers, no more sons.  Is there an unvisited component?
!
  else
 
    do i = 1, nnode
 
      if ( order(i) == 0 ) then
        itake = 1
        iput = 1
        queue(iput) = i
        jdeep = jdeep + 1
        nudeep = 1
        k = k + 1
        order(i) = k
        deep(i) = jdeep
        dad(i) = 0
        go to 10
      end if
 
    end do
 
  end if
 
  return
end
subroutine graph_adj_bipartite_random ( nnode1, nnode2, nedge, adj )

!*****************************************************************************80
!
!! graph_adj_bipartite_random() generates a random bipartite graph.
!
!  Definition:
!
!    A bipartite graph has the property that its nodes may be divided
!    into two groups, NODE1 and NODE2, with the property that the only
!    edges in the graph are between a node in NODE1 and a node in NODE2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    28 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NNODE1, NNODE2, the number of nodes in the 
!    first and second groups of nodes.
!
!    Output, integer NEDGE, the number of edges.
!
!    Output, integer ADJ(NNODE1+NNODE2,NNODE1+NNODE2), the adjacency matrix.
!    ADJ(I,J) is nonzero if there is an edge from node I to node J.  ADJ(I,I) 
!    will always be 0.
!
  implicit none

  integer nnode1
  integer nnode2
  integer nedge

  integer adj(nnode1+nnode2,nnode1+nnode2)
  integer i
  integer i4_uniform_ab
  integer j
  integer k
  integer nnode

  if ( nnode1 <= 0  ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRAPH_ADJ_BIPARTITE_RANDOM(): Fatal error!'
    write ( *, '(a,i8)' ) '  NNODE1 = ', nnode1
    write ( *, '(a)' ) '  but NNODE1 must be at least 1.'
    stop 1
  end if

  if ( nnode2 <= 0  ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRAPH_ADJ_BIPARTITE_RANDOM(): Fatal error!'
    write ( *, '(a,i8)' ) '  NNODE2 = ', nnode2
    write ( *, '(a)' ) '  but NNODE2 must be at least 1.'
    stop 1
  end if

  nnode = nnode1 + nnode2
  nedge = 0

  adj(1:nnode,1:nnode) = 0
!
!  For each node in the NODE1 group, 
!  consider a edge to each node in the NODE2 group.
!
  do i = 1, nnode1
    do j = nnode1 + 1, nnode1 + nnode2

      k = i4_uniform_ab ( 0, 1 )

      adj(i,j) = k
      adj(j,i) = k
      nedge = nedge + k

    end do
  end do
!
!  Now perform a random permutation of the rows and columns.
!
  call i4mat_perm_random ( nnode, adj )

  return
end
subroutine graph_adj_block ( adj, nnode, dad, order, stack, nblock ) 

!*****************************************************************************80
!
!! graph_adj_block(): blocks of an undirected graph from its adjacency list.
!
!  Definition:
!
!    A component of a graph is a connected subset of the graph.  If a node
!    is in the component, then all nodes to which it is connected are also
!    in the component.
!
!    An articulation point of a component of a graph is a node whose
!    removal causes the component to no longer be connected.
!
!    A component with no articulation points is called a block.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 April 1999
!
!  Reference:
!
!    Alan Gibbons,
!    Algorithmic Graph Theory,
!    Cambridge University Press, 1985,
!    ISBN 0-521-28881-9.
!
!  Parameters:
!
!    Input/output, integer adj(nnode,nnode).
!    On input, ADJ is the adjacency matrix.  ADJ(I,J) is
!    positive if there is an edge from node I to node J, and 0 otherwise.
!    On output, each positive entry of ADJ has been replaced
!    by the number of the block that the corresponding edge belongs to.
!
!    Input, integer NNODE, the number of nodes.
!
!    Output, integer DAD(NNODE), DAD(I) is the node from which
!    node I is visited.  Node 1 is the first node in the search,
!    and has no predecessor, so DAD(1) is zero.  If there is
!    more than one connected component in the graph, then there
!    will be other nodes with DAD equal to zero.
!
!    Output, integer ORDER(NNODE).  ORDER(I) records the order
!    in which the node was visited during the depth-first search.
!    The first node, node 1, has ORDER(1) = 1.
!    Note, however, that any node which is an articulation point
!    will have the value of ORDER negated.
!
!    Workspace, integer STACK(NNODE).
!
!    Output, integer NBLOCK, the number of blocks.
!
  implicit none

  integer nnode

  integer adj(nnode,nnode)
  integer dad(nnode)
  integer i
  integer idir
  integer ii
  integer inode(nnode)
  integer order(nnode)
  integer iroot
  integer j
  integer jedge
  integer jj
  integer jnode(nnode)
  integer k
  integer l
  integer label(nnode)
  integer lstack
  integer nblock
  integer stack(nnode)

  dad(1:nnode) = 0
  inode(1:nnode) = 0
  order(1:nnode) = 0
  stack(1:nnode) = 0
  jnode(1:nnode) = 0
  label(1:nnode) = 0
 
  nblock = 0
  k = 0
  i = 1
  lstack = 0
  jedge = 0
!
!  Find all descendants of the parent node in this connected component
!  of the graph.
!
10    continue
 
  iroot = i
  k = k + 1
  order(i) = k
  label(i) = k
  lstack = lstack + 1
  stack(lstack) = i
  idir = + 1
 
30    continue
 
  j = 0
!
!  Check the next neighbor.
!
40    continue
 
  j = j + 1

  if ( nnode < j ) then
    go to 50
  end if
 
  if ( adj(i,j) /= 0 .or. adj(j,i) /= 0 ) then
 
    if ( 0 < adj(i,j) .or. 0 < adj(j,i) ) then
      jedge = jedge + 1
      inode(jedge) = i
      jnode(jedge) = j
    end if
 
    if ( order(j) == 0 ) then
 
      dad(j) = i
      lstack = lstack + 1
      stack(lstack) = j
      idir = + 1
      k = k + 1
      i = j
      order(i) = k
      label(i) = k
      go to 30
 
    else
 
      if ( idir == +1 ) then
        label(i) = min ( label(i), abs ( order(j) ) )
      else
        label(i) = min ( label(i), label(j) )
      end if
 
    end if
 
  end if
 
  go to 40
!
!  Searched all directions from current node.  Back up one node,
!  or, if stack is exhausted, look for a node we haven't visited,
!  which therefore belongs to a new connected component.
!
50    continue
 
  lstack = lstack - 1
  idir = -1
 
  if ( 0 < lstack ) then
 
    j = i
    i = stack(lstack)
 
    if ( abs ( order(i) ) <= label(j) ) then
 
      if ( 0 < order(i) ) then
 
        if ( i /= iroot ) then
          order(i) = - order(i)
        else
          iroot = 0
        end if
 
      end if
 
      nblock = nblock + 1

      do

        ii = inode(jedge)
        jj = jnode(jedge)
        jedge = jedge - 1
        adj(ii,jj) = - nblock
        adj(jj,ii) = - nblock

        if ( ii == i .and. jj == j ) then
          exit
        end if

      end do
 
    end if
 
    go to 40
 
  else
 
    lstack = 0
 
    do l = 1, nnode
      if ( order(l) == 0 ) then
        i = l
        go to 10
      end if
    end do
 
  end if
!
!  Restore the positive sign of the adjacency matrix.
!
  adj(1:nnode,1:nnode) = abs ( adj(1:nnode,1:nnode) )
 
  return
end
subroutine graph_adj_transitive_closure ( adj, nnode, c )

!*****************************************************************************80
!
!! graph_adj_transitive_closure() generates the transitive closure of a graph.
!
!  Discussion:
!
!    The method is due to Stephen Warshall.
!
!    The transitive closure of a graph is a function REACH(I,J) so that
!
!      REACH(I,J) = 0 if node J cannot be reached from node I;
!                   1 if node J can be reached from node I.
!
!    This is an extension of the idea of adjacency.  REACH(I,J)=1 if
!    node J is adjacent to node I, or if node J is adjacent to a node
!    that is adjacent to node I, etc.
!
!    Note that if a graph is (node) connected, then its transitive closure
!    is the matrix that is 1 everywhere.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 March 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Sedgewick,
!    Algorithms,
!    Addison Wesley, 1983, page 425.
!
!  Input:
!
!    integer adj(nnode,nnode): the adjacency for the graph.
!
!    integer NNODE, the number of nodes.
!
!  Output:
!
!    integer c(nnode,nnode): the adjacency for the transitive closure.
!
  implicit none

  integer nnode

  integer adj(nnode,nnode)
  integer c(nnode,nnode)
  integer i
  integer j
  integer k

  c(1:nnode,:nnode) = adj(1:nnode,1:nnode)

  do i = 1, nnode
    c(i,i) = 1
  end do

  do i = 1, nnode
    do j = 1, nnode
      if ( c(j,i) /= 0 .or. c(i,j) /= 0 ) then
        do k = 1, nnode
          if ( c(i,k) /= 0 .or. c(k,i) /= 0 ) then
            c(j,k) = 1
            c(k,j) = 1
          end if
        end do
      end if
    end do
  end do

  return
end
subroutine graph_adj_color_cand ( adj, nnode, ncolor, color, k, nstack, &
  stack, maxstack, ncan )

!*****************************************************************************80
!
!! graph_adj_color_cand(): possible colors for a node during a graph coloring.
!
!  Discussion:
!
!    This routine is given a partial coloring of the graph.  
!    The total coloring of the graph must be done in such a way that no 
!    two nodes joined by an edge have the same color.
!
!    This routine must be used in conjunction with I4VEC_BACKTRACK.
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
!    Input, integer adj(nnode,nnode), the adjacency matrix.  ADJ(I,J)
!    is nonzero if there is an edge from node I to node J.
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, integer NCOLOR, the number of colors available.
!
!    Input, integer COLOR(NNODE).  COLOR(I) is the color of node I.
!
!    Input, integer K, node whose possible colors are to be found.
!
!    Input/output, integer NSTACK, current length of stack.
!
!    Workspace, integer STACK(MAXSTACK), candidates for the colors of the nodes.
!
!    Input, integer MAXSTACK, dimension of STACK.
!
!    input, integer NCAN(NNODE), the number of candidates for 
!    each position.
!
  implicit none

  integer nnode
  integer maxstack

  integer adj(nnode,nnode)
  integer color(nnode)
  integer i
  logical lwork(nnode)
  integer k
  integer ncan(nnode)
  integer nstack
  integer ncolor
  integer stack(maxstack)

  ncan(k) = 0

  if ( k <= 1 ) then

    stack(1) = 1
    nstack = 1
    ncan(k) = 1

  else
 
    lwork(1:ncolor) = .true.
 
    do i = 1, k-1
      if ( adj(i,k) /= 0 .or. adj(k,i) /= 0 ) then
        lwork(color(i)) = .false.
      end if
    end do
  
    do i = 1, ncolor
 
      if ( lwork(i) ) then
        nstack = nstack + 1
        stack(nstack) = i
        ncan(k) = ncan(k) + 1
      end if
 
    end do
 
  end if
 
  return
end
subroutine graph_adj_color_next ( adj, nnode, ncolor, color, stack, &
  maxstack, ncan, more )

!*****************************************************************************80
!
!! graph_adj_color_next() returns possible colorings of a graph, one at a time.
!
!  Definition:
!
!    A coloring of a graph using NCOLOR colors is an assignment to each
!    node of a label between 1 and NCOLOR, in such a way that no two
!    neighboring nodes have the same label.
!
!  Method:
!
!    This routine uses the backtracking method to produce the colorings.
!    Routine GRAPH_ADJ_COLOR_CAND produces candidates for a partial solution, 
!    and routine I4VEC_BACKTRACK assembles the total solution.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 July 1998
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
!    Input, integer adj(nnode,nnode), the adjacency matrix.  
!    ADJ(I,J) is nonzero if there is an edge between node I to node J.
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, integer NCOLOR, the number of colors available.
!
!    Output, integer COLOR(NNODE).  On return with MORE = TRUE, 
!    COLOR(I) is the color of node I.
!
!    Workspace, integer STACK(MAXSTACK), candidates for the colors of nodes
!    1 through K-1.
!
!    Input, integer MAXSTACK, dimension of STACK.
!
!    Workspace, integer NCAN(NNODE), the number of candidates for each position.
!
!    Input/output, logical MORE.
!    On first call, set MORE to .FALSE, and do not alter it after.
!    On return, MORE is TRUE if another coloring has been returned in
!    IARRAY, and FALSE if there are no more colorings.
!
  implicit none

  integer nnode
  integer maxstack

  integer adj(nnode,nnode)
  integer color(nnode)
  integer, save :: indx = 0
  integer, save :: k = 0
  logical more
  integer ncan(nnode)
  integer, save :: nstack = 0
  integer ncolor
  integer stack(maxstack)
!
!  First call.
!
  if ( .not. more ) then
    indx = 0
    k = 0
    more = .true.
    nstack = 0
  end if
 
  do
 
    call i4vec_backtrack ( nnode, color, indx, k, nstack, stack, maxstack, &
      ncan )
 
    if ( indx == 1 ) then

      exit

    else if ( indx == 2 ) then

      call graph_adj_color_cand ( adj, nnode, ncolor, color, k, nstack, &
        stack, maxstack, ncan )

    else

      more = .false.
      exit

    end if

  end do
 
  return
end
subroutine graph_adj_complement ( adj, nnode, c )

!*****************************************************************************80
!
!! graph_adj_complement(): the adjacency matrix of the complement of a graph.
!
!  Definition:
!
!    The complement of a graph G is a graph H with the property that
!    nodes u and v are connected in H if and only if they are not
!    connected in G.  
!
!    However, edges from a node to itself are not allowed.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer adj(nnode,nnode): the adjacency for the graph.
!
!    Input, integer NNODE, the number of nodes.
!
!  Output:
!
!    integer c(nnode,nnode): the adjacency for the complement.
!
  implicit none

  integer nnode

  integer adj(nnode,nnode)
  integer c(nnode,nnode)
  integer i
  integer j

  c(1:nnode,1:nnode) = adj(1:nnode,1:nnode)
!
!  Force the adjacency graph to be symmetric.
!
  call graph_adj_symmetrize ( c, nnode )
!
!  Compute the complement.
!
  do i = 1, nnode
    do j = 1, nnode

      if ( i == j ) then
        c(i,j) = 0
      else if ( c(i,j) /= 0 ) then
        c(i,j) = 0
      else if ( c(i,j) == 0 ) then
        c(i,j) = 1
      end if

    end do
  end do
 
  return
end
subroutine graph_adj_connect_random ( nnode, nedge, adj )

!*****************************************************************************80
!
!! graph_adj_connect_random() generates a random connected graph.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    28 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ADJ(NNODE,NNODE), the adjacency matrix.  
!    ADJ(I,J) is nonzero if there is an edge from node I to node J.  
!    ADJ(I,I) will always be 0.
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, integer NEDGE, the number of edges, which must be
!    between NNODE-1 and (NNODE*(NNODE-1))/2.  
!
  implicit none

  integer nnode
  integer nedge

  integer adj(nnode,nnode)
  integer code(nnode-2)
  integer i
  integer inode(nnode-1)
  integer iwork(nedge)
  integer j
  integer jnode(nnode-1)
  integer k
  integer l
  integer maxedge
  integer nchoice
  integer nchoose
!
!  Check.
!
  if ( nnode <= 0  ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRAPH_ADJ_CONNECT_RANDOM(): Fatal error!'
    write ( *, '(a,i8)' ) '  NNODE = ', nnode
    write ( *, '(a)' ) '  but NNODE must be at least 1.'
    stop 1
  end if

  maxedge = ( nnode * ( nnode - 1 ) ) / 2

  if ( nedge < nnode-1 .or. maxedge < nedge ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRAPH_ADJ_CONNECT_RANDOM(): Fatal error!'
    write ( *, '(a,i8)' ) '  NEDGE = ', nedge
    write ( *, '(a)' ) '  but NEDGE must be at least 0, and '
    write ( *, '(a,i8)' ) '  no more than ', maxedge
    stop 1
  end if
!
!  Initialize the adjacency matrix.
!
  adj(1:nnode,1:nnode) = 0
!
!  Pick a random tree.
!
  call tree_arc_random ( nnode, code, inode, jnode )
!
!  Convert information to adjacency form.
!
  call graph_arc_to_graph_adj ( nnode - 1, inode, jnode, adj, nnode )
!
!  Now we have NEDGE - ( NNODE - 1 ) more edges to add.
!
  nchoice = ( nnode * ( nnode - 1 ) ) / 2 - ( nnode - 1 )
  nchoose = nedge - ( nnode - 1 )

  call ksub_random ( nchoice, nchoose, iwork )

  k = 0
  l = 1
  do i = 1, nnode
    do j = i + 1, nnode
      if ( adj(i,j) /= 0 .or. adj(j,i) /= 0 ) then
        k = k + 1

        if ( l <= nchoose ) then
          if ( iwork(l) == k ) then
            adj(i,j) = 1
            adj(j,i) = 1
            l = l + 1
          end if
        end if

      end if
    end do
  end do

  return
end
subroutine graph_adj_cycle ( adj, nnode, dad, order, maxstack, stack )

!*****************************************************************************80
!
!! graph_adj_cycle() searches for cycles in a graph.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 February 1999
!
!  Parameters:
!
!    Input/output, integer adj(nnode,nnode), the adjacency matrix
!    for the graph.  ADJ(I,J) is 0 if there is no edge from node I to node J.
!    On input, ADJ(I,J) should be 1 if there is an edge from node I to node J.
!    On output, ADJ(I,J) will be one of the following values:
!      -1 if the edge from node I to node J is part of at least one cycle;
!      -2 if the edge from node I to node J is part of the depth first
!      search trees.
!
!    Input, integer NNODE, the number of nodes.
!
!    Output, integer DAD(NNODE), the father array for the depth
!    first search trees.  DAD(I) = 0 means that node I is the root of 
!    one of the trees.  DAD(I) = J means that the search descended
!    from node J to node I.
!
!    Output, integer ORDER(NNODE), the order in which the nodes
!    were traversed, from 1 to NNODE.
!
!    Input, integer MAXSTACK, the amount of stack space available.
!    The absolute maximum needed would be 2*(NNODE-1) though less
!    is likely.
!
!    Workspace, integer STACK(MAXSTACK).
!
  implicit none

  integer maxstack
  integer nnode

  integer adj(nnode,nnode)
  integer dad(nnode)
  integer daddy
  integer i
  integer j
  integer nstack
  integer order(nnode)
  integer rank
  integer stack(maxstack)

  dad(1:nnode) = 0
  order(1:nnode) = 0

  rank = 0

  do i = 1, nnode

    if ( order(i) == 0 ) then

      daddy = i
      nstack = 0
!
!  Visit node DAD.
!
10    continue

      rank = rank + 1
      order(daddy) = rank
      j = 0
!
!  Consider visiting node J from node DAD.
!
20    continue

      j = j + 1
!
!  If J is a reasonable value, adjacent to DAD, and unvisited,
!  then put DAD into the stack, make J the new value of DAD,
!  and go to 10.
!
      if ( j <= nnode ) then

        if ( 0 < adj(daddy,j) .or. 0 < adj(j,daddy) ) then

          if ( order(j) == 0 ) then

            adj(daddy,j) = - 2
            adj(j,daddy) = - 2

            if ( nstack+2 <= maxstack ) then
              dad(j) = daddy
              stack(nstack+1) = daddy
              stack(nstack+2) = j
              nstack = nstack + 2
              daddy = j
              go to 10
            else
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'GRAPH_CYCLE(): Fatal error!'
              write ( *, '(a)' ) '  Out of stack space.'
              stop 1
            end if
!
!  An adjacent node has already been visited.  This constitutes a cycle.
!
          else

            adj(j,daddy) = - 1
            adj(daddy,j) = - 1
            go to 20

          end if
!
!  If J is not suitable for a visit, get the next value of J.
!
        else

          go to 20

        end if
!
!  If no more neighbors to consider, back up one node.
!
      else if ( 2 <= nstack ) then

        daddy = stack(nstack-1)
        j = stack(nstack)
        nstack = nstack - 2
        go to 20
!
!  If no more nodes to consider in this tree, bail out.
!
      else

        nstack = 0

      end if

    end if

  end do

  return
end
subroutine graph_adj_degree ( adj, nnode, degree )

!*****************************************************************************80
!
!! graph_adj_degree() computes the degree of each node.
!
!  Discussion:
!
!    The degree of a node is the number of edges that are incident on it.
!    The sum of the degrees of the nodes is twice the number of edges.
!
!    The generalized case, where ADJ(I,J) can be greater than 1, indicating
!    the existence of 2 or more distinct edges between nodes I and J,
!    will be properly handled by this routine.  
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
!    integer adj(nnode,nnode): the adjacency for the graph.
!
!    integer NNODE: the number of nodes.
!
!  Output:
!
!    integer DEGREE(NNODE), the degree of the nodes.
!
  implicit none

  integer nnode

  integer adj(nnode,nnode)
  integer degree(nnode)
  integer i
  integer j

  degree(1:nnode) = 0

  do i = 1, nnode
    do j = 1, nnode
      degree(i) = degree(i) + adj(i,j)
    end do
  end do

  return
end
subroutine graph_adj_degree_max ( adj, nnode, degree_max )

!*****************************************************************************80
!
!! graph_adj_degree_max() computes the maximum node degree.
!
!  Discussion:
!
!    The maximum node degree of a graph is the maximum value of the
!    degree of the nodes of the graph.
!
!    If two graphs are isomorphic, they must have the same maximum node degree.
!
!    If two graphs have different maximum node degrees, they cannot
!    be isomorphic.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer adj(nnode,nnode), the adjacency information.
!
!    integer NNODE, the number of nodes.
!
!  Output:
!
!    integer DEGREE_MAX, the maximum node degree of the graph.
!
  implicit none

  integer nnode

  integer adj(nnode,nnode)
  integer degree
  integer degree_max
  integer i
  integer j

  degree_max = 0

  do i = 1, nnode
    degree = 0
    do j = 1, nnode
      if ( adj(i,j) /= 0 ) then
        degree = degree + adj(i,j)
      end if
    end do
    degree_max = max ( degree_max, degree )
  end do

  return
end
subroutine graph_adj_degree_sequence ( adj, nnode, seq )

!*****************************************************************************80
!
!! graph_adj_degree_sequence() computes the degree sequence of a graph.
!
!  Discussion:
!
!    The degree sequence of a graph is constructed by computing the
!    degree of each node, and then ordering these values in decreasing order.
!
!    If two graphs are isomorphic, they must have the same degree sequence.
!
!    If two graphs have different degree sequences, they cannot be isomorphic.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer adj(nnode,nnode), the adjacency information.
!    ADJ(I,J) is 1 if there is an edge from node I to node J.
!
!    Input, integer NNODE, the number of nodes.
!
!    Output, integer SEQ(NNODE), the degree sequence of the graph.
!
  implicit none

  integer nnode

  integer adj(nnode,nnode)
  integer i
  integer j
  integer seq(nnode)

  seq(1:nnode) = 0

  do i = 1, nnode
    do j = 1, nnode
      seq(i) = seq(i) + adj(i,j)
    end do
  end do

  call i4vec_sort_heap_d ( nnode, seq )

  return
end
subroutine graph_adj_depth_first ( adj, nnode, dad, order )

!*****************************************************************************80
!
!! graph_adj_depth_first() does a depth first search of a graph.
!
!  Discussion:
!
!    The routine returns:
!
!    * a list of the order in which the nodes were visited,
!    * a list of the parents of each node in the search tree,
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 July 2000
!
!  Reference:
!
!    Robert Sedgewick,
!    Algorithms,
!    Addison Wesley, 1983, page 382.
!
!  Input:
!
!    integer adj(nnode,nnode), the adjacency matrix for 
!    the graph.  ADJ(I,J) is 0 if node J is not adjacent to node I, and nonzero
!    otherwise.
!
!    integer NNODE, the number of nodes.
!
!  Output:
!
!    integer DAD(NNODE), the father array for the depth
!    first search trees.  DAD(I) = 0 means that node I is the root of 
!    one of the trees.  DAD(I) = J means that the search descended
!    from node J to node I.
!
!    integer ORDER(NNODE), the order in which the nodes
!    were traversed, from 1 to NNODE.
!
  implicit none

  integer nnode

  integer adj(nnode,nnode)
  integer dad(nnode)
  integer daddy
  integer i
  integer j
  integer maxstack
  integer nstack
  integer order(nnode)
  integer rank
  integer stack(2*(nnode-1))

  dad(1:nnode) = 0
  maxstack = 2 * ( nnode - 1 )
  order(1:nnode) = 0

  rank = 0

  do i = 1, nnode

    if ( order(i) == 0 ) then

      daddy = i
      nstack = 0
!
!  Visit node DAD.
!
10    continue

      rank = rank + 1
      order(daddy) = rank
      j = 0
!
!  Consider visiting node J from node DAD.
!
20    continue

      j = j + 1
!
!  If J is a reasonable value, adjacent to DAD, and unvisited,
!  then put DAD into the stack, make J the new value of DAD,
!  and go to 10.
!
      if ( j <= nnode ) then

        if ( adj(daddy,j) /= 0 .and. order(j) == 0 ) then

          if ( nstack+2 <= maxstack ) then
            dad(j) = daddy
            stack(nstack+1) = daddy
            stack(nstack+2) = j
            nstack = nstack + 2
            daddy = j
            go to 10
          else
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'graph_adj_depth_first(): Fatal error!'
            write ( *, '(a)' ) '  Out of stack space.'
            stop 1
          end if
!
!  If J is not suitable for a visit, get the next value of J.
!
        else

          go to 20

        end if
!
!  If no more neighbors to consider, back up one node.
!
      else if ( 2 <= nstack ) then

        daddy = stack(nstack-1)
        j = stack(nstack)
        nstack = nstack - 2
        go to 20
!
!  If no more nodes to consider in this tree, bail out.
!
      else

        nstack = 0

      end if

    end if

  end do

  return
end
subroutine graph_adj_depth_first_2 ( adj, nnode, dad, order )

!*****************************************************************************80
!
!! graph_adj_depth_first_2() does a depth-first search of a graph.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 April 1999
!
!  Reference:
!
!    Alan Gibbons,
!    Algorithmic Graph Theory,
!    Cambridge University Press, 1985,
!    ISBN 0-521-28881-9.
!
!  Parameters:
!
!    Input, integer adj(nnode,nnode), is the adjacency matrix of
!    the graph.  ADJ(I,J) is nonzero if there is an edge from node
!    I to node J.
!
!    Input, integer NNODE, the number of nodes.
!
!    Output, integer DAD(NNODE), DAD(I) is the node from which
!    node I is visited.  Node 1 is the first node in the search,
!    and has no predecessor, so DAD(1) is zero.  If there is
!    more than one connected component in the graph, then there
!    will be other nodes with DAD equal to zero.
!
!    Output, integer ORDER(NNODE).  ORDER(I) is the step at which
!    node I is visited in the search.
!
  implicit none

  integer nnode

  integer adj(nnode,nnode)
  integer dad(nnode)
  integer i
  integer order(nnode)
  integer j
  integer kount
  integer l
  integer lstack
  integer stack(nnode)

  order(1:nnode) = 0
  dad(1:nnode) = 0
  stack(1:nnode) = 0
 
  kount = 0
  i = 1
  lstack = 0
!
!  Find all descendents of the parent node I in this connected component
!  of the graph.
!
10    continue
 
  kount = kount + 1
  dad(i) = 0
  order(i) = kount
  lstack = lstack + 1
  stack(lstack) = i
!
!  Check to see if each node, J, is a "descendant" of node I.
!
30    continue
 
  j = 0
!
!  Check next neighbor, J.
!
40    continue
 
  j = j + 1

  if ( j <= nnode ) then
 
    if ( adj(i,j) /= 0 .and. order(j) == 0 ) then

      lstack = lstack + 1
      stack(lstack) = j
      dad(j) = i
      kount = kount + 1
      order(j) = kount
      i = j

      if ( kount == nnode ) then
        return
      end if

      go to 30

    end if

    go to 40

  end if
!
!  Searched all directions from current node.  Back up one node.
!
  lstack = lstack - 1
 
  if ( 0 < lstack ) then
    j = i
    i = stack(lstack)
    go to 40
  end if
!
!  The stack is exhausted.  It's time to look for another connected
!  component.
!
  lstack = 0
 
  do l = 1, nnode
    if ( order(l) == 0 ) then
      i = l
      go to 10
    end if
  end do
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_adj_depth_first_2(): Fatal error!'

  stop 1
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
subroutine graph_adj_eigen ( adj, nnode, neigen, eigen )

!*****************************************************************************80
!
!! graph_adj_eigen() computes eigenvalues of a graph from its adjacency matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer adj(nnode,nnode), the adjacency information.
!    ADJ(I,J) is 1 if there is an edge from node I to node J.
!
!    Input, integer NNODE, the number of nodes.
!
!    Output, integer NEIGEN, the number of eigenvalues computed.
!    Normally, this would be equal to NNODE, unless the algorithm failed.
!
!    Output, real ( kind = rk ) EIGEN(NNODE), contains the computed eigenvalues.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nnode

  real ( kind = rk ) a(nnode,nnode)
  integer adj(nnode,nnode)
  real ( kind = rk ) e(nnode)
  real ( kind = rk ) e2(nnode)
  real ( kind = rk ) eigen(nnode)
  integer ierr
  integer neigen

  a(1:nnode,1:nnode) = real ( adj(1:nnode,1:nnode), kind = rk )

  call tred1 ( nnode, nnode, a, eigen, e, e2 )

  call tqlrat ( nnode, eigen, e2, ierr )

  if ( ierr == 0 ) then
    neigen = nnode
  else
    neigen = ierr - 1
  end if

  return
end
subroutine graph_adj_example_bush ( adj )

!*****************************************************************************80
!
!! graph_adj_example_bush() sets up the adjacency information for the bush graph.
!
!  Diagram:
!
!        6   3
!        |   |
!    1---4---5---2
!        |
!        7
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Output:
!
!    integer ADJ(7,7): the adjacency information for the graph.  
!
  implicit none

  integer, parameter :: nnode = 7

  integer adj(nnode,nnode)

  adj(1:nnode,1:nnode) = 0

  adj(1,4) = 1

  adj(2,5) = 1

  adj(3,5) = 1

  adj(4,1) = 1
  adj(4,5) = 1
  adj(4,6) = 1
  adj(4,7) = 1

  adj(5,2) = 1
  adj(5,3) = 1
  adj(5,4) = 1

  adj(6,4) = 1

  adj(7,4) = 1

  return
end
subroutine graph_adj_example_cube ( adj )

!*****************************************************************************80
!
!! graph_adj_example_cube() sets up the adjacency information for the cube graph.
!
!  Diagram:
!
!      4-----7
!     /|    /|
!    8-----3 |
!    | |   | |
!    | 5---|-2
!    |/    |/
!    1-----6
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ADJ(8,8), the adjacency information for 
!    the graph.  ADJ(I,J) is 1 if nodes I and J are adjacent and 0 otherwise.
!
  implicit none

  integer, parameter :: nnode = 8

  integer adj(nnode,nnode)

  adj(1:nnode,1:nnode) = 0

  adj(1,5) = 1
  adj(1,6) = 1
  adj(1,8) = 1

  adj(2,5) = 1
  adj(2,6) = 1
  adj(2,7) = 1

  adj(3,6) = 1
  adj(3,7) = 1
  adj(3,8) = 1

  adj(4,5) = 1
  adj(4,7) = 1
  adj(4,8) = 1

  adj(5,1) = 1
  adj(5,2) = 1
  adj(5,4) = 1

  adj(6,1) = 1
  adj(6,2) = 1
  adj(6,3) = 1

  adj(7,2) = 1
  adj(7,3) = 1
  adj(7,4) = 1

  adj(8,1) = 1
  adj(8,3) = 1
  adj(8,4) = 1

  return
end
subroutine graph_adj_example_dodecahedron ( adj )

!*****************************************************************************80
!
!! graph_adj_example_dodecahedron(): adjacency info for the dodecahedron graph.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ADJ(NNODE,NNODE), the adjacency information for 
!    the graph.  ADJ(I,J) is 1 if nodes I and J are adjacent and 0 otherwise.
!
  implicit none

  integer, parameter :: nnode = 20

  integer adj(nnode,nnode)

  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1
  adj(1,5) = 1
  adj(1,6) = 1

  adj(2,1) = 1
  adj(2,3) = 1
  adj(2,8) = 1

  adj(3,2) = 1
  adj(3,4) = 1
  adj(3,10) = 1

  adj(4,3) = 1
  adj(4,5) = 1
  adj(4,12) = 1

  adj(5,1) = 1
  adj(5,4) = 1
  adj(5,14) = 1

  adj(6,1) = 1
  adj(6,7) = 1
  adj(6,15) = 1

  adj(7,6) = 1
  adj(7,8) = 1
  adj(7,17) = 1

  adj(8,7) = 1
  adj(8,9) = 1
  adj(8,2) = 1

  adj(9,8) = 1
  adj(9,10) = 1
  adj(9,16) = 1

  adj(10,3) = 1
  adj(10,9) = 1
  adj(10,11) = 1

  adj(11,10) = 1
  adj(11,12) = 1
  adj(11,20) = 1

  adj(12,4) = 1
  adj(12,11) = 1
  adj(12,13) = 1

  adj(13,12) = 1
  adj(13,14) = 1
  adj(13,19) = 1

  adj(14,13) = 1
  adj(14,15) = 1
  adj(14,5) = 1

  adj(15,6) = 1
  adj(15,14) = 1
  adj(15,18) = 1

  adj(16,9) = 1
  adj(16,17) = 1
  adj(16,20) = 1

  adj(17,16) = 1
  adj(17,18) = 1
  adj(17,7) = 1

  adj(18,15) = 1
  adj(18,17) = 1
  adj(18,19) = 1

  adj(19,13) = 1
  adj(19,18) = 1
  adj(19,20) = 1

  adj(20,11) = 1
  adj(20,16) = 1
  adj(20,19) = 1

  return
end
subroutine graph_adj_example_octo ( example, adj )

!*****************************************************************************80
!
!! graph_adj_example_octo() sets up an 8 node example graph.
!
!  Diagram:
!
!      1---2
!     /|   |\
!    8-+---+-3
!    | |   | |
!    7-+---+-4
!     \|   |/
!      6---5
!
!     Graph "A"
!
!    There are 7 graphs to choose from.  They are all on 8 nodes.  The first
!    5 have degree 3 at every node.  Graphs 6 and 7 have degree 5 at every
!    node.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    28 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer EXAMPLE, should be between 1 and 7, and 
!    indicates which example graph to pick.
!
!    Output, integer ADJ(8,8), the adjacency information for 
!    the graph.  ADJ(I,J) is 1 if nodes I and J are adjacent and 0 otherwise.  
!
  implicit none

  integer, parameter :: nnode = 8

  integer adj(nnode,nnode)
  integer example
  integer i
  integer i4_uniform_ab
  integer j
  integer nsave

  if ( example <= 0 ) then
    nsave = i4_uniform_ab ( 1, 7 )
  else
    example = mod ( example - 1, 7 ) + 1
    nsave = example
  end if

  adj(1:nnode,1:nnode) = 0

  do i = 1, nnode
    j = i + 1
    if ( nnode < j ) then
      j = j - nnode
    end if

    adj(i,j) = 1
    adj(j,i) = 1

  end do

  if ( nsave == 1 ) then
    
    adj(1,6) = 1
    adj(6,1) = 1
    adj(2,5) = 1
    adj(5,2) = 1
    adj(3,8) = 1
    adj(8,3) = 1
    adj(4,7) = 1
    adj(7,4) = 1

  else if ( nsave == 2 ) then

    adj(1,6) = 1
    adj(6,1) = 1
    adj(2,8) = 1
    adj(8,2) = 1
    adj(3,5) = 1
    adj(5,3) = 1
    adj(4,7) = 1
    adj(7,4) = 1

  else if ( nsave == 3 ) then

    adj(1,5) = 1
    adj(5,1) = 1
    adj(2,6) = 1
    adj(6,2) = 1
    adj(3,7) = 1
    adj(7,3) = 1
    adj(4,8) = 1
    adj(8,4) = 1

  else if ( nsave == 4 ) then

    adj(1,3) = 1
    adj(3,1) = 1
    adj(2,4) = 1
    adj(4,2) = 1
    adj(5,7) = 1
    adj(7,5) = 1
    adj(6,8) = 1
    adj(8,6) = 1

  else if ( nsave == 5 ) then

    adj(1,4) = 1
    adj(4,1) = 1
    adj(2,6) = 1
    adj(6,2) = 1
    adj(3,8) = 1
    adj(8,3) = 1
    adj(5,7) = 1
    adj(7,5) = 1

  else if ( nsave == 6 ) then

    adj(1,4) = 1
    adj(1,5) = 1
    adj(1,6) = 1

    adj(2,5) = 1
    adj(2,6) = 1
    adj(2,7) = 1

    adj(3,6) = 1
    adj(3,7) = 1
    adj(3,8) = 1

    adj(4,7) = 1
    adj(4,8) = 1
    adj(4,1) = 1

    adj(5,8) = 1
    adj(5,1) = 1
    adj(5,2) = 1

    adj(6,1) = 1
    adj(6,2) = 1
    adj(6,3) = 1

    adj(7,2) = 1
    adj(7,3) = 1
    adj(7,4) = 1

    adj(8,3) = 1
    adj(8,4) = 1
    adj(8,5) = 1

  else if ( nsave == 7 ) then

    adj(1,3) = 1
    adj(1,5) = 1
    adj(1,7) = 1

    adj(2,4) = 1
    adj(2,6) = 1
    adj(2,8) = 1

    adj(3,5) = 1
    adj(3,7) = 1
    adj(3,1) = 1

    adj(4,6) = 1
    adj(4,8) = 1
    adj(4,2) = 1

    adj(5,7) = 1
    adj(5,1) = 1
    adj(5,3) = 1

    adj(6,8) = 1
    adj(6,2) = 1
    adj(6,4) = 1

    adj(7,1) = 1
    adj(7,3) = 1
    adj(7,5) = 1

    adj(8,2) = 1
    adj(8,4) = 1
    adj(8,6) = 1

  end if
!
!  Now permute the graph.
!
  call i4mat_perm_random ( nnode, adj )

  return
end
subroutine graph_adj_example_twig ( adj )

!*****************************************************************************80
!
!! graph_adj_example_twig() sets up the adjacency information for the twig graph.
!
!  Diagram:
!
!    1---2---3
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ADJ(NNODE,NNODE), the adjacency information for
!    the graph.  ADJ(I,J) is 1 if nodes I and J are adjacent.  
!
  implicit none

  integer, parameter :: nnode = 3

  integer adj(nnode,nnode)

  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1

  adj(2,1) = 1
  adj(2,3) = 1

  adj(3,2) = 1

  return
end
subroutine graph_adj_ham_cand ( adj, nnode, circuit, k, nstack, &
  stack, maxstack, ncan )

!*****************************************************************************80
!
!! graph_adj_ham_cand(): candidates for the next node in a Hamiltonian circuit.
!
!  Discussion:
!
!    This routine is used in conjunction with I4VEC_BACKTRACK.  
!
!  Definition:
!
!    A Hamiltonian circuit of a graph is a path that starts at a given node, 
!    visits every node exactly once, and returns to the starting node.
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
!    Input, integer adj(nnode,nnode).  ADJ(I,J) = 1 if there is
!    an edge from node I to node J, 0 otherwise.
!
!    Input, integer NNODE, the number of nodes in the graph.
!
!    Input, integer CIRCUIT(NNODE), the nodes forming the circuit.
!
!    Input, integer K, index of the next node to be determined
!    for the circuit.
!
!    Input/output, integer NSTACK, the current length of stack.
!
!    Input, integer STACK(MAXSTACK), candidates for steps 1...K-1.
!
!    Input, integer MAXSTACK, the dimension of STACK.
!
!    Workspace, integer NCAN(NNODE), the number of candidates for
!    positions in the circuit.
!
  implicit none

  integer nnode
  integer maxstack

  integer adj(nnode,nnode)
  integer circuit(nnode)
  integer i
  integer iwork(nnode)
  integer k
  integer ncan(nnode)
  integer nstack
  integer stack(maxstack)

  ncan(k) = 0

  if ( k == 1 ) then
    stack(1) = 1
    nstack = 1
    ncan(k) = 1
    return
  end if
 
  iwork(1:nnode) = adj(circuit(k-1),1:nnode)
 
  iwork(circuit(1:k-1)) = 0
  
  if ( k /= nnode ) then
 
    do i = 1, nnode
      if ( iwork(i) == 1 ) then
        if ( maxstack <= nstack ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'GRAPH_ADJ_HAM_CAND(): Fatal error!'
          write ( *, '(a)' ) '  Stack size exceeded.'
          stop 1
        end if
        nstack = nstack + 1
        stack(nstack) = i
        ncan(k) = ncan(k) + 1
      end if
    end do
 
    return
 
 else if ( k == nnode ) then
 
    do i = 1, nnode
 
      if ( iwork(i) == 1 ) then
 
        if ( circuit(2) < i .or. adj(i,1) == 0 ) then

        else
          if ( maxstack <= nstack ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'GRAPH_ADJ_HAM_CAND(): Fatal error!'
            write ( *, '(a)' ) '  Stack size exceeded.'
            stop 1
          end if
          nstack = nstack + 1
          stack(nstack) = i
          ncan(k) = ncan(k) + 1
        end if

        return
 
      end if
 
    end do

  end if
 
  return
end
subroutine graph_adj_ham_next ( adj, nnode, circuit, stack, maxstack, &
  ncan, more )

!*****************************************************************************80
!
!! graph_adj_ham_next() returns the next Hamilton circuit for a graph.
!
!  Discussion:
!
!    The routine produces all the Hamilton circuits of a graph, one at a time.
!
!    A Hamiltonian circuit of a graph is a path that starts at a given
!    node, visits every node exactly once, and returns to the starting node.
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
!    Input, integer adj(nnode,nnode).  ADJ(I,J) = 1 if there is
!    an edge from node I to node J, 0 otherwise.
!
!    Input, integer NNODE, the number of nodes in graph.
!
!    Input, integer CIRCUIT(NNODE).  CIRCUIT(I) is the I-th node
!    in the circuit. 
!
!    Input, integer K, the index of the next node to be determined
!    for circuit.
!
!    Input/output, integer NSTACK, the current length of stack.
!
!    Input, integer STACK(MAXSTACK).  Candidates for steps 1...K-1.
!
!    Input, integer MAXSTACK, dimension of STACK.
!
!    Workspace, integer NCAN(NNODE), the number of candidates for each
!    position in the circuit.
!
!    Input/output, logical MORE.
!    On first call, set MORE to .FALSE, and do not alter it after.
!    On return, MORE is TRUE if another circuit has been returned in
!    IARRAY, and FALSE if there are no more circuits.
!
  implicit none

  integer nnode
  integer maxstack

  integer adj(nnode,nnode)
  integer circuit(nnode)
  integer, save :: indx = 0
  integer, save :: k = 0
  logical more
  integer ncan(nnode)
  integer, save :: nstack = 0
  integer stack(maxstack)

  if ( .not. more ) then
    indx = 0
    k = 0
    more = .true.
    nstack = 0
  end if
 
  do
 
    call i4vec_backtrack ( nnode, circuit, indx, k, nstack, stack, maxstack, &
      ncan )
 
    if ( indx == 1 ) then

      exit

    else if ( indx == 2 ) then

      call graph_adj_ham_cand ( adj, nnode, circuit, k, nstack, &
        stack, maxstack, ncan )

    else

      more = .false.
      exit

    end if

  end do
 
  return
end
subroutine graph_adj_ham_next_brute ( adj, nnode, circuit, iset )

!*****************************************************************************80
!
!! graph_adj_ham_next_brute() finds the next Hamiltonian circuit in a graph.
!
!  Discussion:
!
!    This is a brute force algorithm, and not suitable for large problems.
!    It is really only useful as a demonstration, and as a check for
!    the backtracking algorithm.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer adj(nnode,nnode), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link between nodes I and J.
!
!    Input, integer NNODE, the number of nodes.
!
!    Input/output, integer CIRCUIT(NNODE).
!
!    On input, if ISET = 0, then CIRCUIT is not presumed to contain any 
!    information.  If ISET is nonzero, then CIRCUIT contains the circuit 
!    computed on the previous call.
!
!    On output, CIRCUIT contains the circuit computed by this call.
!
!    Input/output, integer ISET.
!    On input, 0 means this is the first call for this graph.  
!    Any other value means this is a repeated call for more circuits.
!
!    On output, a 0 value means that no more circuits could be computed.
!    Otherwise, ISET is incremented by one on each call.
!
  implicit none

  integer nnode

  integer adj(nnode,nnode)
  integer circuit(nnode)
  integer i
  integer ipos
  integer iset
!
!  If ISET is 0, this is a starting call, and we set CIRCUIT
!  to the lexically first circuit to check.
!
!  Otherwise, we set CIRCUIT to the next permutation.
!
  if ( iset == 0 ) then
    ipos = 0
    circuit(1:nnode) = 0
  else
    ipos = nnode - 1
  end if
 
10 continue
 
  call perm_inc ( circuit, ipos, nnode )

  if ( ipos <= 0 .or. circuit(1) /= 1 ) then
    iset = 0
    circuit(1:nnode) = 0
    return
  end if
!
!  Check whether the entries of CIRCUIT actually form a circuit.
!  If we find a break in the circuit, store that location in IPOS
!  and move on to try the next permutation.
!
  do i = 1, nnode - 1
    ipos = i
    if ( adj(circuit(i),circuit(i+1)) == 0 ) then
      go to 10
    end if
  end do
!
!  If the circuit connects all the nodes, we only have to check whether
!  the last node connects back to the first one.
!
!  To cut down the pairs of equivalent circuits created by going one
!  way or the other over the same set of nodes, we also require that,
!  for 2 < NNODE, the last node be numbered higher than the second one.
!
  if ( adj(circuit(nnode),circuit(1)) == 0 ) then
    go to 10
  end if

  if ( 2 < nnode ) then
    if ( circuit(nnode) < circuit(2) ) then
      go to 10
    end if
  end if
 
  iset = iset + 1

  return
end
subroutine graph_adj_is_bipartite ( adj, nnode, result )

!*****************************************************************************80
!
!! graph_adj_is_bipartite() determines if a graph is bipartite.
!
!  Definition:
!
!    A graph is bipartite if its nodes can be divided into two subsets
!    in such a way that every edge joins a node from each subset.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    25 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer adj(nnode,nnode), the adjacency for the graph.  
!
!    integer NNODE, the number of nodes.
!
!  Output:
!
!    integer RESULT.
!    0, the graph is not bipartite.
!    1, the graph is bipartite.
!
  implicit none

  integer nnode

  integer adj(nnode,nnode)
  integer i
  integer j
  integer k
  integer khi
  integer klo
  integer lhi
  integer list(nnode)
  integer oldset
  integer result
  integer set
  integer subset(nnode)

  result = 1
  subset(1:nnode) = -1
!
!  Node 1 is put in subset 1.
!
  set = 1
  list(1) = 1
  subset(1) = set

  klo = 1
  khi = 1
!
!  Working from the set of nodes found on the previous step, look
!  for all in and out neighbors.  
!
10    continue

  oldset = set
  set = 1 - set

  lhi = khi
!
!  Consider each node I in the previously found set.
!
  do k = klo, khi

    i = list(k)
!
!  Look at all in and out neighbors J.
!
    do j = 1, nnode

      if ( adj(i,j) /= 0 .or. adj(j,i) /= 0 ) then
!
!  If the node is not in any subset, put it in the other one.
!
        if ( subset(j) == -1 ) then

          lhi = lhi + 1
          list(lhi) = j
          subset(j) = set
!
!  But if the node is in the same subset, bipartiteness has failed.
!
        else if ( subset(j) == oldset ) then
          result = 0
          return
        end if

      end if

    end do

  end do
!
!  Assuming we found more nodes, on this sweep, then ...
!
  if ( khi < lhi ) then
    klo = khi + 1
    khi = lhi 
    go to 10
  end if
!
!  Assuming we found no new nodes on this sweep, see if there are any
!  nodes we have missed.  These will be completely isolated from all the
!  nodes we have found so far.
!
  do i = 1, nnode

    if ( subset(i) == -1 ) then
      klo = khi + 1
      khi = klo
      subset(i) = set
      list(klo) = i
      go to 10
    end if

  end do

  result = 1

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
subroutine graph_adj_is_eulerian ( adj, nnode, result )

!*****************************************************************************80
!
!! graph_adj_is_eulerian() determines if a graph is Eulerian.
!
!  Discussion:
!
!    A graph is path-Eulerian if there exists a path through the graph
!    which uses every edge once.
!
!    A graph is circuit-Eulerian if there exists a path through the graph
!    which uses every edge once, and which starts and ends on the same node.
!
!    Note that it is NOT necessary for the path or circuit to pass through
!    every node; simply that all the edges can be used exactly once to
!    make a connected path.  This means an Eulerian graph can have isolated
!    nodes, for instance.
!
!    A graph is path-Eulerian if and only if it is edge connected, and all 
!    but two nodes are of even degree.
!
!    A graph is circuit-Eulerian if and only if it is edge connected and
!    all nodes are of even degree.
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
!  Input:
!
!    integer adj(nnode,nnode), the adjacency matrix for the  graph.
!
!    integer NNODE, the number of nodes.
!
!  Output:
!
!    integer RESULT.
!    0, the graph is not Eulerian.
!    1, the graph is path-Eulerian.
!    2, the graph is circuit-Eulerian.
!
  implicit none

  integer nnode

  integer adj(nnode,nnode)
  integer degree
  integer i
  logical is_connected
  integer j
  integer nodd
  integer result
!
!  First check that the graph is edgewise connected.
!
  call graph_adj_is_edge_connected ( adj, nnode, is_connected )

  if ( .not. is_connected ) then
    result = 0
    return
  end if
!
!  Now look at node degree.
!
  nodd = 0

  do i = 1, nnode

    degree = 0

    do j = 1, nnode
      if ( adj(i,j) /= 0 ) then
        if ( i == j ) then
          degree = degree + 2
        else
          degree = degree + 1
        end if
      end if
    end do

    if ( mod ( degree, 2 ) == 1 ) then
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
subroutine graph_adj_edges_random ( nnode, nedge, adj )

!*****************************************************************************80
!
!! graph_adj_edges_random() generates a random graph on NNODE nodes with NEDGE edges.
!
!  Discussion:
!
!    Each edge will show up twice in the adjacency matrix.
!    The number of edges must be between 0 and (NNODE*(NNODE-1))/2.
!    ADJ(I,J) is nonzero if there is an edge from node I to node J.
!    ADJ(I,I) will always be 0.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 September 2006
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
!  Output:
!
!    integer ADJ(NNODE,NNODE), the adjacency matrix.  
!
  implicit none

  integer nnode
  integer nedge

  integer adj(nnode,nnode)
  integer i
  integer iwork(nedge)
  integer j
  integer k
  integer l
  integer maxedge
!
!  Check.
!
  if ( nnode <= 0  ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'graph_adj_edges_random(): Fatal error!'
    write ( *, '(a,i8)' ) '  NNODE = ', nnode
    write ( *, '(a)' ) '  but NNODE must be at least 1.'
    stop 1
  end if

  maxedge = ( nnode * ( nnode - 1 ) ) / 2

  if ( nedge < 0 .or. maxedge < nedge ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'graph_adj_edges_random(): Fatal error!'
    write ( *, '(a,i8)' ) '  NEDGE = ', nedge
    write ( *, '(a)' ) '  but NEDGE must be at least 0, and '
    write ( *, '(a,i8)' ) '  no more than ', maxedge
    stop 1
  end if
!
!  Initialize the adjacency matrix.
!
  adj(1:nnode,1:nnode) = 0
!
!  Pick IWORK, a random subset of MAXEDGE, of size NEDGE.
!
  call ksub_random ( maxedge, nedge, iwork )
!
!  The usable spots in the superdiagonal are numbered as follows:
!
!  * 1  2   3  ...  n-1
!  * * n+1 n+2 ... 2n-3
!  ...
!  * *  *   *  ... (n*(n-1))/2
!  * *  *   *  ...   * 
!
  k = 0
  l = 1
  do i = 1, nnode - 1
    do j = i + 1, nnode

      k = k + 1

      if ( l <= nedge ) then

        if ( k == iwork(l) ) then
          adj(i,j) = 1
          adj(j,i) = 1
          l = l + 1
        end if

      end if

    end do
  end do

  return
end
subroutine graph_adj_random ( nnode, prob, adj )

!*****************************************************************************80
!
!! graph_adj_random() generates a random graph on NNODE nodes.
!
!  Discussion:
!
!    The user specifies the probability P that an edge will be generated
!    between any pair of nodes.
!
!    ADJ(I,J) is nonzero if there is an edge from node I to node J.  
!    ADJ(I,I) will always be 0.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer NNODE, the number of nodes.
!
!    real ( kind = rk ) PROB, the probability that an edge will
!    be generated between any given pair of nodes.
!
!  Output:
!
!    integer ADJ(NNODE,NNODE), the adjacency matrix.  
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nnode

  integer adj(nnode,nnode)
  integer i
  integer j
  real ( kind = rk ) p(nnode,nnode)
  real ( kind = rk ) prob

  if ( prob < 0.0D+00 .or. 1.0D+00 < prob ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'graph_adj_random(): Fatal error!'
    write ( *, '(a)' ) '  Input probability PROB is not between 0 and 1.'
    stop 1
  end if

  call random_number ( harvest = p(1:nnode,1:nnode) )

  do j = 1, nnode
    do i = j + 1, nnode
      p(i,j) = p(j,i)
    end do
  end do

  p = ( p + transpose ( p ) ) / 2.0 

  adj(1:nnode,1:nnode) = 0
  where ( p <= prob )
    adj = 1
  end where

  do i = 1, nnode
    adj(i,i) = 0
  end do

  return
end
subroutine graph_adj_reduce ( adj, nnode )

!*****************************************************************************80
!
!! graph_adj_reduce() generates a transitive reduction of a graph.
!
!  Discussion:
!
!    This routine is given an adjacency matrix B, which might be a
!    transitive closure of a graph G.
!
!    The transitive closure graph is generated from a graph G by the 
!    following procedure:
!
!      B(I,J) = 0 if node J cannot be reached from node I in graph G;
!               1 if node J can be reached from node I in graph G.
!
!    The purpose of this routine is to try to find the original, sparser
!    graph G which generated the given transitive closure graph.  Such a
!    graph G is known as a transitive reduction..  In general,
!    there is no unique solution.  In particular, any graph is a transitive
!    reduction of itself.  
!
!    Hence, the real task is to drop as many redundant edges as possible
!    from the given graph, arriving at a graph from which no more edges 
!    may be removed.
!
!  Method:
!
!    One way of explaining the algorithm is based on the adjacency matrix:
!
!    * Zero out the diagonals of the adjacency matrix.
!
!    * Consider row 1.  Any other row that can "reach" row 1 doesn't
!      need a 1 if row 1 has it.  So "subtract" all the 1's in row 1
!      from such rows.  We are done with row 1 and column 1.
!
!    * Repeat for the other rows.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ADJ(NNODE,NNODE).
!    On input, the adjacency matrix of the transitive closure graph H.
!    On output, the adjacency matrix of a transitive reduction graph G 
!    of the graph H.
!
!    Input, integer NNODE, the number of nodes.
!
  implicit none

  integer nnode

  integer adj(nnode,nnode)
  integer i
  integer j
  integer k
!
!  First discard those useless self-edges.
!
  do i = 1, nnode
    adj(i,i) = 0
  end do
!
!  If you can get from J to I and I to K, you don't need an
!  edge from J to K.
!
  do i = 1, nnode
    do j = 1, nnode
      if ( adj(j,i) /= 0 .or. adj(i,j) /= 0 ) then
        do k = 1, nnode
          if ( adj(i,k) /= 0 .or. adj(k,i) /= 0 ) then
            adj(j,k) = 0
            adj(k,j) = 0
          end if
        end do
      end if
    end do
  end do

  return
end
subroutine graph_adj_span_tree ( adj, nnode, inode, jnode )

!*****************************************************************************80
!
!! graph_adj_span_tree() finds a spanning tree of a graph.
!
!  Discussion:
!
!    If the graph is connected, NNODE-1 edges comprise the spanning tree.  
!
!    If the graph is not connected, but divided into NCOMP components, then 
!    NNODE-NCOMP edges will comprise the spanning "forest", and the other 
!    edges will be zero.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer adj(nnode,nnode), the adjacency matrix for 
!    the graph.  ADJ(I,J) is 0 if there is no edge from node I to node J.
!
!    Input, integer NNODE, the number of nodes.
!
!    Output, integer INODE(NNODE-1), JNODE(NNODE-1), the edge list 
!    for the spanning tree or forest.
!
  implicit none

  integer nnode

  integer adj(nnode,nnode)
  integer i
  integer inode(nnode-1)
  integer j
  integer jnode(nnode-1)
  integer label(nnode)
  integer level
  integer nedge
  integer nfound
  integer nlabel

  label(1:nnode) = 0

  inode(1:nnode-1) = 0
  jnode(1:nnode-1) = 0

  level = 0
  nedge = 0
  nlabel = 0
!
!  Find an unvisited node.
!
  do

    i = 0
 
    do

      i = i + 1

      if ( label(i) == 0 ) then
        exit
      end if

    end do
  
    label(i) = level + 1
    nlabel = nlabel + 1
!
!  Search for all nodes reachable from the node.
!
    do

      level = level + 1
      nfound = 0

      do i = 1, nnode

        if ( label(i) == level ) then

          do j = 1, nnode

            if ( label(j) == 0 ) then

              if ( adj(i,j) /= 0 .or. adj(j,i) /= 0 ) then
                label(j) = level + 1
                nlabel = nlabel + 1
                nfound = nfound + 1
                nedge = nedge + 1
                inode(nedge) = i
                jnode(nedge) = j
              end if

            end if

          end do

        end if

      end do

      if ( nfound <= 0 ) then
        exit
      end if

    end do
!
!  If we have labeled all nodes, exit.
!
    if ( nnode <= nlabel ) then
      exit
    end if

  end do

  return
end
subroutine graph_adj_span_tree_enum ( adj, nnode, tree_num )

!*****************************************************************************80
!
!! graph_adj_span_tree_enum() enumerates the spanning trees of a graph.
!
!  Discussion:
!
!    If ADJ is the adjacency matrix of the graph, let A be the matrix
!      A = DEG - ADJ
!    where DEG is the diagonal matrix with DEG(I,I) = degree of node I.
!    Then the number of spanning trees of the graph is equal to the
!    determinant of any cofactor of A.  A cofactor of A is obtained by
!    deleting a row and column of A.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer adj(nnode,nnode), the adjacency matrix for 
!    the graph.  ADJ(I,J) is 0 if there is no edge from node I to node J.
!
!    Input, integer NNODE, the number of nodes.
!
!    Output, integer TREE_NUM, the number of spanning trees.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nnode

  real ( kind = rk ) a(nnode,nnode)
  integer adj(nnode,nnode)
  integer degree(nnode)
  real ( kind = rk ) det
  integer i
  integer info
  integer ipivot(nnode)
  integer tree_num
!
!  Construct the matrix.
!
  call graph_adj_degree ( adj, nnode, degree )

  a(1:nnode,1:nnode) = - real ( adj(1:nnode,1:nnode), kind = rk )

  do i = 1, nnode
    a(i,i) = a(i,i) + real ( degree(i), kind = rk )
  end do
!
!  Factor the NNODE-1 order matrix.
!
  call dge_fa ( nnode-1, a(1:nnode-1,1:nnode-1), ipivot, info )

  if ( info /= 0 ) then
    tree_num = 0
    return
  end if
!
!  Get the determinant.
!
  call dge_det ( nnode-1, a(1:nnode-1,1:nnode-1), ipivot, det )

  tree_num = nint ( det )

  return
end
subroutine graph_adj_symmetrize ( adj, nnode )

!*****************************************************************************80
!
!! graph_adj_symmetrize() symmetrizes an adjacency matrix.
!
!  Discussion:
!
!    For a graph, if there is an edge from I to J, there is an edge from
!    J to I.  Therefore, the adjacency matrix should be symmetric.  
!    This routine enforces that condition.  If either ADJ(I,J) or ADJ(J,I)
!    is nonzero, the output adjacency matrix will have both entries nonzero.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer adj(nnode,nnode).  On output, the 
!    adjacency information has been symmetrized.
!
!    Input, integer NNODE, the number of nodes.
!
  implicit none

  integer nnode

  integer adj(nnode,nnode)
  integer i
  integer j
!
!  While not perfect, this method does not assume that 1 is the only
!  legal nonzero value in ADJ.
!
  do i = 1, nnode
    do j = i + 1, nnode
      if ( adj(i,j) /= 0 ) then
        adj(j,i) = adj(i,j)
      else if ( adj(j,i) /= 0 ) then
        adj(i,j) = adj(j,i)
      end if
    end do
  end do
 
  return
end
subroutine graph_adj_to_graph_arc ( adj, nnode, maxedge, nedge, inode, &
  jnode )

!*****************************************************************************80
!
!! graph_adj_to_graph_arc() converts an adjacency graph to an arc list graph.
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
!    Input, integer adj(nnode,nnode), the adjacency matrix for the 
!    graph.  ADJ(I,J) is nonzero if there is an edge from node I to node J.
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, integer MAXEDGE, the maximum number of edges.
!
!    Output, integer NEDGE, the number of edges.
!
!    Output, integer INODE(MAXEDGE), JNODE(MAXEDGE), the arc list
!    of the graph.
!
  implicit none

  integer maxedge
  integer nnode

  integer adj(nnode,nnode)
  integer i
  integer inode(maxedge)
  integer j
  integer jnode(maxedge)
  integer nedge

  nedge = 0

  inode(1:maxedge) = 0
  jnode(1:maxedge) = 0

  do j = 1, nnode
    do i = j, nnode
      if ( adj(i,j) /= 0 .or. adj(j,i) /= 0 ) then
        nedge = nedge + 1
        if ( nedge <= maxedge ) then
          inode(nedge) = i
          jnode(nedge) = j
        else
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'GRAPH_ADJ_TO_GRAPH_ARC(): Fatal error!'
          write ( *, '(a)' ) '  MAXEDGE exceeded.'
          stop 1
        end if
      end if
    end do
  end do

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
function i4_uniform_ab ( a, b )

!*****************************************************************************80
!
!! i4_uniform_ab() returns a scaled pseudorandom I4 between A and B.
!
!  Discussion:
!
!    An I4 is an integer value.
!
!    The pseudorandom number will be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Input:
!
!    integer A, B, the limits of the interval.
!
!  Output:
!
!    integer I4_UNIFORM_AB, a number between A and B.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer a
  integer b
  integer i4_uniform_ab
  real ( kind = rk ) r

  call random_number ( harvest = r )
  i4_uniform_ab = a + int ( ( b + 1 - a ) * r )

  return
end
subroutine i4mat_perm_random ( n, a )

!*****************************************************************************80
!
!! i4mat_perm_random() selects a random permutation of an I4MAT.
!
!  Discussion:
!
!    The matrix is assumed to be square.  A single permutation is
!    applied to both rows and columns.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    28 March 2005
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
!    Input, integer N, the number of rows and columns in the array.
!
!    Input/output, integer A(N,N), the array to be permuted.
!
  implicit none

  integer n

  integer a(n,n)
  integer i
  integer i2
  integer i4_uniform_ab
  integer j
!
!  Permute the rows and columns together.
!
  do i = 1, n

    i2 = i4_uniform_ab ( i, n )

    do j = 1, n
      call i4_swap ( a(i2,j), a(i,j) )
    end do

    do j = 1, n
      call i4_swap ( a(j,i2), a(j,i) )
    end do

  end do

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
subroutine i4vec_heap_a ( n, a )

!*****************************************************************************80
!
!! i4vec_heap_a() reorders an array of integers into an ascending heap.
!
!  Definition:
!
!    An ascending heap is an array A with the property that, for every index J,
!    A(J) <= A(2*J) and A(J) <= A(2*J+1), (as long as the indices
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
!  If both positions exist, take the smaller of the two values,
!  and update M if necessary.
!
        if ( a(m+1) < a(m) ) then
          m = m + 1
        end if

      end if
!
!  If the small descendant is smaller than KEY, move it up,
!  and update IFREE, the location of the free position, and
!  consider the descendants of THIS position.
!
      if ( a(m) < key ) then
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
subroutine i4vec_sort_heap_d ( n, a )

!*****************************************************************************80
!
!! i4vec_sort_heap_d() descending sorts an integer array using heap sort.
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
!  1: Put A into ascending heap form.
!
  call i4vec_heap_a ( n, a )
!
!  2: Sort A.
!
!  The smallest object in the heap is in A(1).
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
    call i4vec_heap_a ( n1, a )
!
!  Take the smallest object from A(1) and move it to A(N1).
!
    call i4_swap ( a(1), a(n1) )

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
subroutine ksub_random ( n, k, iarray )

!*****************************************************************************80
!
!! ksub_random() selects a random subset of size K from a set of size N.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    28 March 2005
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
!    Input, integer N, the size of the set from which subsets are
!    drawn.
!
!    Input, integer K, number of elements in desired subsets.  
!    K must be between 0 and N.
!
!    Output, integer IARRAY(K).  IARRAY(I) is the I-th element of
!    the output set.  The elements of IARRAY are in order.
!
  implicit none

  integer k

  integer i
  integer i4_uniform_ab
  integer iarray(k)
  integer ids
  integer ihi
  integer ip
  integer ir
  integer is
  integer ix
  integer l
  integer ll
  integer m
  integer m0
  integer n

  if ( k < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUB_RANDOM(): Fatal error!'
    write ( *, '(a,i8)' ) '  K = ', k
    write ( *, '(a)' ) '  but 0 <= K is required!'
    stop 1
  else if ( n < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUB_RANDOM(): Fatal error!'
    write ( *, '(a,i8)' ) '  N = ', n
    write ( *, '(a,i8)' ) '  K = ', k
    write ( *, '(a)' ) '  K <= N is required!'
    stop 1
  end if

  if ( k == 0 ) then
    return
  end if

  do i = 1, k
    iarray(i) = ( ( i - 1 ) * n ) / k
  end do

  do i = 1, k

    do

      ix = i4_uniform_ab ( 1, n )

      l = 1 + ( ix * k - 1 ) / n

      if ( iarray(l) < ix ) then
        exit
      end if

    end do

    iarray(l) = iarray(l) + 1

  end do

  ip = 0
  is = k

  do i = 1, k

    m = iarray(i)
    iarray(i) = 0

    if ( m /= ( (i-1) * n ) / k ) then
      ip = ip + 1
      iarray(ip) = m
    end if

  end do

  ihi = ip

  do i = 1, ihi
    ip = ihi + 1 - i
    l = 1 + ( iarray(ip) * k - 1 ) / n
    ids = iarray(ip) - ( ( l - 1 ) * n ) / k
    iarray(ip) = 0
    iarray(is) = l
    is = is - ids
  end do

  do ll = 1, k

    l = k + 1 - ll

    if ( iarray(l) /= 0 ) then
      ir = l
      m0 = 1 + ( ( iarray(l) - 1 ) * n ) / k
      m = ( iarray(l) * n ) / k - m0 + 1
    end if

    ix = i4_uniform_ab ( m0, m0+m-1 )

    i = l + 1

    do while ( i <= ir )

      if ( ix < iarray(i) ) then
        exit
      end if

      ix = ix + 1
      iarray(i-1) = iarray(i)
      i = i + 1

    end do

    iarray(i-1) = ix
    m = m - 1

  end do

  return
end
subroutine perm_inc ( iperm, ipos, n )

!*****************************************************************************80
!
!! perm_inc() "increments" a permutation to get the "next" one.
!
!  Discussion:
!
!    The routine is given IPERM, a permutation of the numbers from 1 to N,
!    and a position IPOS between 1 and N.
!
!    It returns the next permutation in the dictionary order which
!    comes after all permutations beginning IPERM(1) through IPERM(IPOS).
!
!  Example:
!
!             PERM              IPOS
!
!    Input    123456789         7
!    Output   123456798         7
!
!    Input    123456789         9
!    Output   213456789         0
!
!    Input    134826795         3
!    Output   134925678         3
!
!    Input    134826795         0
!    Output   123456789         0
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
!  Parameters:
!
!    Input/output, integer IPERM(N).
!    On input, the current permutation.
!    On output, the "incremented" permutation.
!
!    Input/output, integer IPOS.
!    On input, IPOS is the location of the end of the string of
!    "digits" in IPERM that form the test string.  That is, the
!    new permutation to be computed must be the very next one,
!    in dictionary order, which succeeds all strings whose first
!    IPOS digits agree with the input value of IPERM.
!
!    On output, IPOS is the position of the last digit of the output
!    value of IPERM which agrees with the input value of IPERM.
!
!    Input, integer N, is the number of entries in IPERM.
!
  implicit none

  integer n

  integer ipcopy
  integer iperm(n)
  integer ipos
  integer j
  integer k
  integer new

  if ( ipos == 0 ) then
    ipos = n
    call i4vec_indicator ( n, iperm )
    return
  end if
 
  ipcopy = ipos

10    continue
!
!  To get the next permutation, we need to increment the IPOS+1 "digit".
!
!  We do this by finding, if possible, a digit in positions IPOS+2
!  through N that is just larger than the current value IPOS+1 digit.
!  If we find such a digit, it becomes the IPOS+1 digit, and the
!  remaining values are sorted into increasing order.
!
  new = 0
  do j = ipcopy+2, n
    if ( new == 0 ) then
      if ( iperm(ipcopy+1) < iperm(j) ) then
        new = j
      end if
    else
      if ( iperm(ipcopy+1) < iperm(j) .and. iperm(j) < iperm(new) ) then
        new = j
      end if
    end if
  end do
!
!  There is a next candidate that agrees with IPERM through entry I.
!  Swap entries IPOS+1 and NEW, and sort the entries (IPOS+2,...,N).
!
!  The output value of IPOS equals the input value.
!
  if ( new /= 0 ) then

    call i4_swap ( iperm(new), iperm(ipcopy+1) )
 
    do j = ipcopy+2, n
 
      do k = j+1, n
        if ( iperm(k) < iperm(j) ) then
          call i4_swap ( iperm(j), iperm(k) )
        end if
      end do
 
    end do
    return
  end if
!
!  There is no next candidate that agrees with IPERM through entry 
!  IPOS.  Can we decrease IPOS and try for a next candidate that way?
!
  if ( 0 < ipcopy ) then
    ipcopy = ipcopy - 1
    go to 10
  end if
!
!  IPOS is now zero.  There is no successor to the current permutation,
!  so we start again at the first permutation.
!
  ipos = 0
  call i4vec_indicator ( n, iperm )
 
  return
end
subroutine pruefer_to_tree_arc ( nnode, iarray, inode, jnode )

!*****************************************************************************80
!
!! pruefer_to_tree_arc() is given a Pruefer code, and computes the tree.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 October 1999
!
!  Reference:
!
!    Dennis Stanton, Dennis White,
!    Constructive Combinatorics,
!    Springer Verlag, New York, 1986.
!
!  Parameters:
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, integer IARRAY(NNODE-2), the Pruefer code of the tree.
!
!    Output, integer INODE(NNODE-1), JNODE(NNODE-1), the edge
!    array of the tree.  The I-th edge joins nodes INODE(I) and JNODE(I).
!
  implicit none

  integer nnode

  integer i
  integer iarray(nnode-2)
  integer ii
  integer inode(nnode-1)
  integer iwork(nnode)
  integer j
  integer jnode(nnode-1)
!
!  Initialize IWORK(I) to count the number of neighbors of node I.
!  The Pruefer code uses each node one less time than its total
!  number of neighbors.
!
  iwork(1:nnode) = 1
 
  do i = 1, nnode-2
    iwork(iarray(i)) = iwork(iarray(i)) + 1
  end do
!
!  Now process each entry in the Pruefer code.
!
  do i = 1, nnode-2
 
    ii = 0
    do j = 1, nnode
      if ( iwork(j) == 1 ) then
        ii = j
      end if
    end do
 
    inode(i) = ii
    jnode(i) = iarray(i)
    iwork(ii) = 0
    iwork(iarray(i)) = iwork(iarray(i)) - 1
 
  end do
 
  inode(nnode-1) = iarray(nnode-2)
 
  if ( iarray(nnode-2) /= 1 ) then
    jnode(nnode-1) = 1
  else
    jnode(nnode-1) = 2
  end if
 
  return
end
function pythag ( a, b )

!*****************************************************************************80
!
!! pythag() computes SQRT ( A^2 + B^2 ) carefully.
!
!  Discussion:
!
!    The formula
!
!      PYTHAG = sqrt ( A^2 + B^2 )
!
!    is reasonably accurate, but the formula can actually fail if
!    for example, A^2 is larger than the machine overflow.  The
!    formula can lose most of its accuracy if the sum of the squares
!    is very large or very small.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 March 2000
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow, 
!    Y Ikebe, V Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, the two legs of a right triangle.
!
!    Output, real ( kind = rk ) PYTHAG, the length of the hypotenuse.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) p
  real ( kind = rk ) pythag
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) t
  real ( kind = rk ) u

  p = max ( abs ( a ), abs ( b ) )

  if ( p /= 0.0D+00 ) then

    r = ( min ( abs ( a ), abs ( b ) ) / p )**2

   10   continue

    t = 4.0D+00 + r

    if ( t /= 4.0D+00 ) then
      s = r / t
      u = 1.0D+00 + 2.0D+00 * s
      p = u * p
      r = ( s / u )**2 * r
      go to 10
    end if

  end if

  pythag = p

  return
end
subroutine r8_swap ( x, y )

!*****************************************************************************80
!
!! r8_swap() swaps two double precision values.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = rk ) X, Y.  On output, the values of X and
!    Y have been interchanged.
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
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! r8vec_print() prints an R8VEC.
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
!    Input, real ( kind = rk ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  integer i
  character ( len = * ) title

  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i8,g14.6)' ) i, a(i)
  end do

  return
end
subroutine tqlrat ( n, d, e2, ierr )

!*****************************************************************************80
!
!! tqlrat() computes all eigenvalues of a real symmetric tridiagonal matrix.
!
!  Discussion:
!
!    This subroutine finds the eigenvalues of a symmetric
!    tridiagonal matrix by the rational QL method.
!
!  Reference:
!
!    Christian Reinsch,
!    Algorithm 464, TQLRAT,
!    Communications of the ACM,
!    Volume 16, page 689, 1973.
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Y Ikebe, V Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real ( kind = rk ) D(N).  On input, D contains the diagonal
!    elements of the matrix.  On output, D contains the eigenvalues in ascending
!    order.  If an error exit was made, then the eigenvalues are correct
!    in positions 1 through IERR-1, but may not be the smallest eigenvalues.
!
!    Input/output, real ( kind = rk ) E2(N), contains in positions 2 through N
!    the squares of the subdiagonal elements of the matrix.  E2(1) is
!    arbitrary.  On output, E2 has been overwritten by workspace
!    information.
!
!    Output, integer IERR, error flag.
!    0, for no error,
!    J, if the J-th eigenvalue could not be determined after 30 iterations.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) d(n)
  real ( kind = rk ) e2(n)
  real ( kind = rk ) f
  real ( kind = rk ) g
  real ( kind = rk ) h
  integer i
  integer ierr
  integer ii
  integer j
  integer l
  integer l1
  integer m
  integer mml
  real ( kind = rk ) p
  real ( kind = rk ) pythag
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) t

  ierr = 0

  if ( n == 1 ) then
    return
  end if

  do i = 2, n
    e2(i-1) = e2(i)
  end do

  f = 0.0D+00
  t = 0.0D+00
  e2(n) = 0.0D+00

  do l = 1, n

     j = 0
     h = abs ( d(l) ) + sqrt ( e2(l) )

     if ( t <= h ) then

       t = h
       b = abs ( t ) * epsilon ( b )
       c = b * b

     end if
!
!  Look for small squared sub-diagonal element.
!
     do m = l, n
       if ( e2(m) <= c ) then
         exit
       end if
     end do

     if ( m == l ) then
       go to 210
     end if

130  continue

     if ( 30 <= j ) then
       ierr = l
       return
     end if

     j = j + 1
!
!  Form shift.
!
     l1 = l + 1
     s = sqrt ( e2(l) )
     g = d(l)
     p = ( d(l1) - g ) / ( 2.0D+00 * s )
     r = pythag ( p, 1.0D+00 )
     d(l) = s / ( p + sign ( r, p ) )
     h = g - d(l)

     do i = l1, n
       d(i) = d(i) - h
     end do

     f = f + h
!
!  Rational QL transformation.
!
     g = d(m)
     if ( g == 0.0D+00 ) g = b
     h = g
     s = 0.0D+00
     mml = m - l

     do ii = 1, mml
       i = m - ii
       p = g * h
       r = p + e2(i)
       e2(i+1) = s * r
       s = e2(i) / r
       d(i+1) = h + s * (h + d(i))
       g = d(i) - e2(i) / g
       if ( g == 0.0D+00 ) then
         g = b
       end if
       h = g * p / r
     end do

     e2(l) = s * g
     d(l) = h
!
!  Guard against underflow in convergence test.
!
     if ( h == 0.0D+00 ) go to 210
     if ( abs ( e2(l) ) <= abs ( c / h ) ) go to 210
     e2(l) = h * e2(l)
     if ( e2(l) /= 0.0D+00 ) go to 130

210  continue

     p = d(l) + f
!
!  Order the eigenvalues.
!
     do ii = 2, l
       i = l + 2 - ii
       if ( d(i-1) <= p ) then
         go to 270
       end if
       d(i) = d(i-1)
     end do

     i = 1
270  continue
     d(i) = p

  end do

  return
end
subroutine tred1 ( nm, n, a, d, e, e2 )

!*****************************************************************************80
!
!! tred1() transforms a real symmetric matrix to tridiagonal form.
!
!  Discussion:
!
!    The routine reduces a real symmetric matrix to a symmetric
!    tridiagonal matrix using orthogonal similarity transformations.
!
!  Reference:
!
!    Martin, Reinsch, James Wilkinson,
!    TRED1,
!    Numerische Mathematik,
!    Volume 11, pages 181-195, 1968.
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow, 
!    Y Ikebe, V Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the array A.
!    NM must be at least N.
!
!    Input, integer N, the order of the matrix A.
!
!    Input/output, real ( kind = rk ) A(NM,N), on input, contains the real
!    symmetric matrix.  Only the lower triangle of the matrix need be supplied.
!    On output, A contains information about the orthogonal transformations
!    used in the reduction in its strict lower triangle.
!    The full upper triangle of A is unaltered.
!
!    Output, real ( kind = rk ) D(N), contains the diagonal elements of the
!    tridiagonal matrix.
!
!    Output, real ( kind = rk ) E(N), contains the subdiagonal elements of the
!    tridiagonal matrix in its last N-1 positions.  E(1) is set to zero.
!
!    Output, real ( kind = rk ) E2(N), contains the squares of the corresponding
!    elements of E.  E2 may coincide with E if the squares are not needed.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nm

  real ( kind = rk ) a(nm,n)
  real ( kind = rk ) d(n)
  real ( kind = rk ) e(n)
  real ( kind = rk ) e2(n)
  real ( kind = rk ) f
  real ( kind = rk ) g
  real ( kind = rk ) h
  integer i
  integer ii
  integer j
  integer k
  integer l
  real ( kind = rk ) scale

  d(1:n) = a(n,1:n)

  do i = 1, n
    a(n,i) = a(i,i)
  end do

  do ii = 1, n

    i = n + 1 - ii
    l = i - 1
    h = 0.0D+00
    scale = 0.0D+00

    if ( l < 1 ) go to 130
!
!  Scale row.
!
    do k = 1, l
      scale = scale + abs ( d(k) )
    end do

    if ( scale /= 0.0D+00 ) go to 140

    do j = 1, l
      d(j) = a(l,j)
      a(l,j) = a(i,j)
      a(i,j) = 0.0D+00
    end do

130 continue

    e(i) = 0.0D+00
    e2(i) = 0.0D+00
    go to 300

140 continue

    do k = 1, l
      d(k) = d(k) / scale
      h = h + d(k) * d(k)
    end do

    e2(i) = scale * scale * h
    f = d(l)
    g = - sign ( sqrt ( h ), f )
    e(i) = scale * g
    h = h - f * g
    d(l) = f - g

    if ( l == 1 ) go to 285
!
!  Form A*U.
!
    e(1:l) = 0.0D+00

    do j = 1, l

      f = d(j)
      g = e(j) + a(j,j) * f

      do k = j+1, l
        g = g + a(k,j) * d(k)
        e(k) = e(k) + a(k,j) * f
      end do

      e(j) = g

    end do
!
!  Form P.
!
    f = 0.0D+00

    do j = 1, l
      e(j) = e(j) / h
      f = f + e(j) * d(j)
    end do

    h = f / ( h + h )
!
!  Form Q.
!
    e(1:l) = e(1:l) - h * d(1:l)
!
!  Form reduced A.
!
    do j = 1, l

      f = d(j)
      g = e(j)

      do k = j, l
        a(k,j) = a(k,j) - f * e(k) - g * d(k)
      end do

    end do

  285 continue

    do j = 1, l
      f = d(j)
      d(j) = a(l,j)
      a(l,j) = a(i,j)
      a(i,j) = f * scale
    end do

300 continue

  end do

  return
end
subroutine tree_arc_random ( nnode, code, inode, jnode )

!*****************************************************************************80
!
!! tree_arc_random() selects a random labeled tree and its Pruefer code.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    28 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NNODE, the number of nodes.
!
!    Output, integer CODE(NNODE-2), the Pruefer code for the 
!    labeled tree.
!
!    Output, integer INODE(NNODE-1), JNODE(NNODE-1), the edge 
!    array for the tree.
!
  implicit none

  integer nnode

  integer code(nnode-2)
  integer inode(nnode-1)
  integer jnode(nnode-1)

  if ( nnode <= 0  ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TREE_ARC_RANDOM(): Fatal error!'
    write ( *, '(a,i8)' ) '  NNODE = ', nnode
    write ( *, '(a)' ) '  but NNODE must be at least 1.'
    stop 1
  end if

  if ( nnode <= 2 ) then
    return
  end if

  call vec_random ( nnode-2, nnode, code )
 
  code(1:nnode-2) = code(1:nnode-2) + 1
 
  call pruefer_to_tree_arc ( nnode, code, inode, jnode )
 
  return
end
subroutine vec_random ( n, base, iarray )

!*****************************************************************************80
!
!! vec_random() selects a random N-vector of integers modulo a given base.
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
!  Parameters:
!
!    Input, integer N, the size of the vector to be generated.
!
!    Input, integer BASE, the base to be used.
!
!    Output, integer IARRAY(N), a list of N random values between
!    0 and IBASE-1.
!
  implicit none

  integer n

  integer base
  integer i
  integer i4_uniform_ab
  integer iarray(n)
  integer ival

  do i = 1, n
    ival = i4_uniform_ab ( 0, base - 1 )
    iarray(i) = ival
  end do
 
  return
end
