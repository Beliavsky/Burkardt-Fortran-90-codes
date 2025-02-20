subroutine balanc ( nm, n, a, low, igh, scale )

!*****************************************************************************80
!
!! balanc() balances a real matrix before eigenvalue calculations.
!
!  Discussion:
!
!    This subroutine balances a real matrix and isolates eigenvalues
!    whenever possible.
!
!    Suppose that the principal submatrix in rows LOW through IGH
!    has been balanced, that P(J) denotes the index interchanged
!    with J during the permutation step, and that the elements
!    of the diagonal matrix used are denoted by D(I,J).  Then
!
!      SCALE(J) = P(J),    J = 1,...,LOW-1,
!               = D(J,J),  J = LOW,...,IGH,
!               = P(J)     J = IGH+1,...,N.
!
!    The order in which the interchanges are made is N to IGH+1,
!    then 1 to LOW-1.
!
!    Note that 1 is returned for LOW if IGH is zero formally.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 December 2008
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, 
!    Ikebe, Klema, Moler.
!    This version by John Burkardt.
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
!    Input, integer NM, the leading dimension of A, which must
!    be at least N.
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real ( kind = rk ) A(NM,N), the N by N matrix.  On output,
!    the matrix has been balanced.
!
!    Output, integer LOW, IGH, indicate that A(I,J) is equal to 
!    zero if
!    (1) I is greater than J and
!    (2) J=1,...,LOW-1 or I=IGH+1,...,N.
!
!    Output, real ( kind = rk ) SCALE(N), contains information determining the
!    permutations and scaling factors used.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nm
  integer n

  real ( kind = rk ) a(nm,n)
  real ( kind = rk ) b2
  real ( kind = rk ) c
  real ( kind = rk ) f
  real ( kind = rk ) g
  integer i
  integer iexc
  integer igh
  integer j
  integer k
  integer l
  integer low
  integer m
  logical noconv
  real ( kind = rk ) r
  real ( kind = rk ), parameter :: radix = 16.0D+00
  real ( kind = rk ) s
  real ( kind = rk ) scale(n)

  iexc = 0
  j = 0
  m = 0

  b2 = radix * radix
  k = 1
  l = n
  go to 100

20 continue

  scale(m) = j

  if ( j /= m ) then

    do i = 1, l
      call r8_swap ( a(i,j), a(i,m) )
    end do

    do i = k, n
      call r8_swap ( a(j,i), a(m,i) )
    end do

  end if

  if ( iexc == 2 ) then
    go to 130
  end if
!
!  Search for rows isolating an eigenvalue and push them down.
!
  if ( l == 1 ) then
    low = k
    igh = l
    return
  end if

  l = l - 1

100 continue

  do j = l, 1, -1

     do i = 1, l
       if ( i /= j ) then
         if ( a(j,i) /= 0.0D+00 ) then
           go to 120
         end if
       end if
     end do

     m = l
     iexc = 1
     go to 20

120  continue

  end do

  go to 140
!
!  Search for columns isolating an eigenvalue and push them left.
!
130 continue

  k = k + 1

140 continue

  do j = k, l

    do i = k, l
      if ( i /= j ) then
        if ( a(i,j) /= 0.0D+00 ) then
          go to 170
        end if
      end if
    end do

    m = k
    iexc = 2
    go to 20

170 continue

  end do
!
!  Balance the submatrix in rows K to L.
!
  scale(k:l) = 1.0D+00
!
!  Iterative loop for norm reduction.
!
  noconv = .true.

  do while ( noconv )

    noconv = .false.

    do i = k, l

      c = 0.0D+00
      r = 0.0D+00

      do j = k, l
        if ( j /= i ) then
          c = c + abs ( a(j,i) )
          r = r + abs ( a(i,j) )
        end if
      end do
!
!  Guard against zero C or R due to underflow.
!
      if ( c /= 0.0D+00 .and. r /= 0.0D+00 ) then

        g = r / radix
        f = 1.0D+00
        s = c + r

        do while ( c < g )
          f = f * radix
          c = c * b2
        end do

        g = r * radix

        do while ( g <= c )
          f = f / radix
          c = c / b2
        end do
!
!  Balance.
!
        if ( ( c + r ) / f < 0.95D+00 * s ) then

          g = 1.0D+00 / f
          scale(i) = scale(i) * f
          noconv = .true.

          a(i,k:n) = a(i,k:n) * g
          a(1:l,i) = a(1:l,i) * f

        end if

      end if

    end do

  end do

  low = k
  igh = l

  return
end
subroutine digraph_adjacency_transitive_closure ( adj, nnode )

!*****************************************************************************80
!
!! digraph_adjacency_transitive_closure() generates the transitive closure of a digraph.
!
!  Discussion:
!
!    The method is due to Stephen Warshall.
!
!  Definition:
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
!  Reference:
!
!    Robert Sedgewick,
!    Algorithms,
!    Addison Wesley, 1983, page 425.
!
!  Parameters:
!
!    Input/output, integer ADJ(NNODE,NNODE).
!    On input, ADJ is the adjacency matrix.  ADJ(I,J)
!    is nonzero if there is an edge from node I to node J.
!    On output, ADJ is the transitive closure matrix.
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
!  You can "reach" a node from itself.
!
  do i = 1, nnode
    adj(i,i) = 1
  end do

  do i = 1, nnode
    do j = 1, nnode
      if ( adj(j,i) /= 0 ) then
        do k = 1, nnode
          if ( adj(i,k) /= 0 ) then
            adj(j,k) = 1
          end if
        end do
      end if
    end do
  end do

  return
end
subroutine digraph_adjacency_complement ( adj, nnode, c )

!*****************************************************************************80
!
!! digraph_adjacency_complement(): the adjacency matrix of the complement of a digraph.
!
!  Definition:
!
!    The complement of a digraph G is a digraph H with the property that
!    nodes u and v are connected in H if and only if they are not
!    connected in G.  
!
!    Edges from a node to itself are not allowed.
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
!  Input:
!
!    integer adj(nnode,nnode): the adjacency for the digraph G.
!
!    Input, integer NNODE, the number of nodes.
!
!  Output:
!
!    integer c(nnode,nnode): the adjacency for the complement of G.
!
  implicit none

  integer nnode

  integer adj(nnode,nnode)
  integer c(nnode,nnode)
  integer i
  integer j

  c(1:nnode,1:nnode) = adj(1:nnode,1:nnode)
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
subroutine digraph_adjacency_components ( adj, nnode, ncomp, comp, dad, order )

!*****************************************************************************80
!
!! digraph_adjacency_components() finds the strongly connected components of a digraph.
!
!  Discussion:
!
!    A digraph is a directed graph.
!
!    A strongly connected component of a directed graph is the largest
!    set of nodes such that there is a directed path from any node to 
!    any other node in the same component.
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
!    K Thulasiraman, M Swamy,
!    Graph Theory and Algorithms,
!    John Wiley, New York, 1992.
!
!  Parameters:
!
!    Input, integer ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Input, integer NNODE, the number of nodes.
!
!    Output, integer NCOMP, the number of strongly connected 
!    components.
!
!    Output, integer COMP(NNODE), lists the connected component 
!    to which each node belongs.
!
!    Output, integer DAD(NNODE), the father array for the depth 
!    first search trees.  DAD(I) = 0 means that node I is the root of 
!    one of the trees.  DAD(I) = J means that the search descended
!    from node J to node I.
!
!    Output, integer ORDER(NNODE), the order in which the nodes 
!    were traversed, from 1 to NNODE.
!
  implicit none

  integer nnode

  integer adj(nnode,nnode)
  integer comp(nnode)
  integer dad(nnode)
  integer iorder
  integer lowlink(nnode)
  integer mark(nnode)
  integer ncomp
  integer nstack
  integer order(nnode)
  integer point(nnode)
  integer stack(nnode)
  integer v
  integer w
  integer x
!
!  Initialization.
!
  comp(1:nnode) = 0
  dad(1:nnode) = 0
  order(1:nnode) = 0
  lowlink(1:nnode) = 0
  mark(1:nnode) = 0
  point(1:nnode) = 0

  iorder = 0
  nstack = 0
  ncomp = 0
!
!  Select any node V not stored in the stack, that is, with MARK(V) = 0.
!
  do

    v = 0

    do

      v = v + 1

      if ( nnode < v ) then
        adj(1:nnode,1:nnode) = abs ( adj(1:nnode,1:nnode) )
        return
      end if

      if ( mark(v) /= 1 ) then
        exit
      end if

    end do

    iorder = iorder + 1

    order(v) = iorder
    lowlink(v) = iorder
    mark(v) = 1
 
    nstack = nstack + 1
    stack(nstack) = v
    point(v) = 1

30  continue
!
!  Consider each node W.
!
    do w = 1, nnode
!
!  Is there an edge (V,W) and has it not been examined yet?
!
      if ( 0 < adj(v,w) ) then

        adj(v,w) = - adj(v,w)
!
!  Is the node on the other end of the edge undiscovered yet?
!
        if ( mark(w) == 0 ) then

          iorder = iorder + 1
          order(w) = iorder
          lowlink(w) = iorder
          dad(w) = v
          mark(w) = 1

          nstack = nstack + 1
          stack(nstack) = w
          point(w) = 1

          v = w

        else if ( mark(w) == 1 ) then

          if ( order(w) < order(v) .and. point(w) == 1 ) then
            lowlink(v) = min ( lowlink(v), order(w) )
          end if

        end if

        go to 30

      end if

    end do

    if ( lowlink(v) == order(v) ) then

      ncomp = ncomp + 1

      do

        if ( nstack <= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'DIGRAPH_adjacency_COMPONENTS(): Fatal error!'
          write ( *, '(a)' ) '  Illegal stack reference.'
          stop 1
        end if

        x = stack(nstack)
        nstack = nstack - 1

        point(x) = 0
        comp(x) = ncomp

        if ( x == v ) then
          exit
        end if

      end do

    end if

    if ( dad(v) /= 0 ) then
      lowlink(dad(v)) = min ( lowlink(dad(v)), lowlink(v) )
      v = dad(v)
      go to 30
    end if

  end do

  return
end
subroutine digraph_adjacency_cycle ( adj, nnode, adj2, dad, order )

!*****************************************************************************80
!
!! digraph_adjacency_cycle() searches for cycles in a digraph.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 July 2000
!
!  Parameters:
!
!    Input, integer ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Input, integer NNODE, the number of nodes.
!
!    Output, integer ADJ2(NNODE,NNODE), will be one of the following
!    values depending on the role of the edge from node I to node J:
!       0, no edge,
!       1, neither in a search tree, nor needed to disconnect a cycle;
!      -1, completes a cycle,
!      -2, part of a search tree.
!
!    Output, integer DAD(NNODE), the father array for the depth
!    first search trees.  DAD(I) = 0 means that node I is the root of 
!    one of the trees.  DAD(I) = J means that the search descended
!    from node J to node I.
!
!    Output, integer ORDER(NNODE), the order in which the nodes
!    were traversed, from 1 to NNODE.
!
  implicit none

  integer nnode

  integer adj(nnode,nnode)
  integer adj2(nnode,nnode)
  integer dad(nnode)
  integer daddy
  integer i
  integer j
  integer jj
  integer maxstack
  integer nstack
  integer order(nnode)
  integer rank
  integer stack(2*(nnode-1))
!
!  Initialization.
!
  adj2(1:nnode,1:nnode) = adj(1:nnode,1:nnode)
  dad(1:nnode) = 0
  maxstack = 2 * ( nnode - 1 )
  order(1:nnode) = 0

  rank = 0

  do i = 1, nnode

    if ( order(i) == 0 ) then

      daddy = i
      nstack = 0
!
!  Visit the unvisited node DAD.
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
!  If 
!    J is a reasonable value, 
!    J is adjacent to DAD, and 
!    J is unvisited,
!  then 
!    put DAD into the stack, 
!    make J the new value of DAD, and
!    examine J's neighbors.
!
      if ( j <= nnode ) then

        if ( 0 < adj2(daddy,j) ) then

          if ( order(j) == 0 ) then

            adj2(daddy,j) = -2

            if ( nstack+2 <= maxstack ) then
              dad(j) = daddy
              stack(nstack+1) = daddy
              stack(nstack+2) = j
              nstack = nstack + 2
              daddy = j
              go to 10
            else
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'DIGRAPH_adjacency_CYCLE(): Fatal error!'
              write ( *, '(a)' ) '  Out of stack space.'
              stop 1
            end if
!
!  Adjacent node J has already been visited.  If J is actually
!  in the current stack, then we have a cycle.
!
          else

            if ( j == daddy ) then

              adj2(daddy,j) = - 1

            else

              do jj = 1, nstack-1, 2
                if ( stack(jj) == j ) then
                  adj2(daddy,j) = - 1
                end if
              end do

            end if

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
subroutine digraph_adjacency_degree ( adj, nnode, indegree, outdegree )

!*****************************************************************************80
!
!! digraph_adjacency_degree() computes the indegree and outdegree of each node.
!
!  Discussion:
!
!    The indegree of a node is the number of directed edges that 
!    end at the node.  
!
!    The outdegree of a node is the number of directed edges that
!    begin at the node.
!
!    The sum of the indegrees and outdegrees of all the nodes is twice 
!    the number of edges.
!
!    The generalized case, where ADJ(I,J) can be greater than 1, indicating
!    the existence of 2 or more distinct edges from node I to node J,
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
!  Parameters:
!
!    Input, integer ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Input, integer NNODE, the number of nodes.
!
!    Output, integer INDEGREE(NNODE), OUTDEGREE(NNODE), 
!    the indegree and outdegree of the nodes.
!
  implicit none

  integer nnode

  integer adj(nnode,nnode)
  integer i
  integer indegree(nnode)
  integer j
  integer outdegree(nnode)

  indegree(1:nnode) = 0
  outdegree(1:nnode) = 0

  do i = 1, nnode
    do j = 1, nnode
      if ( adj(i,j) /= 0 ) then
        outdegree(i) = outdegree(i) + adj(i,j)
        indegree(j) = indegree(j) + adj(i,j)
      end if
    end do
  end do

  return
end
subroutine digraph_adjacency_degree_max ( adj, nnode, indegree_max, &
  outdegree_max, degree_max )

!*****************************************************************************80
!
!! digraph_adjacency_degree_max() computes the maximum degrees of a digraph.
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
!    Input, integer ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Input, integer NNODE, the number of nodes.
!
!    Output, integer INDEGREE_MAX, OUTDEGREE_MAX, the maximum 
!    indegree and outdegree, considered independently, which may occur at 
!    different nodes.
!
!    Output, integer DEGREE_MAX, the maximum value of the sum at 
!    each node of the indegree and outdegree.
!
  implicit none

  integer nnode

  integer adj(nnode,nnode)
  integer degree
  integer degree_max
  integer i
  integer indegree
  integer indegree_max
  integer outdegree
  integer outdegree_max

  degree_max = 0
  indegree_max = 0
  outdegree_max = 0

  do i = 1, nnode

    indegree = sum ( adj(1:nnode,i) )
    outdegree = sum ( adj(i,1:nnode) )

    degree = indegree + outdegree

    indegree_max = max ( indegree_max, indegree )
    outdegree_max = max ( outdegree_max, outdegree )
    degree_max = max ( degree_max, degree )
     
  end do

  return
end
subroutine digraph_adjacency_degree_sequence ( adj, nnode, in_seq, out_seq )

!*****************************************************************************80
!
!! digraph_adjacency_degree_sequence() computes the directed degree sequence.
!
!  Discussion:
!
!    The directed degree sequence of a graph is the sequence of indegrees
!    and the sequence of outdegrees, arranged to correspond to nodes of
!    successively decreasing total degree.  For nodes of equal degree, those
!    of higher outdegree take precedence. 
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
!    Input, integer ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Input, integer NNODE, the number of nodes.
!
!    Output, integer IN_SEQ(NNODE), OUT_SEQ(NNODE),
!    the degree sequence of the digraph.
!
  implicit none

  integer nnode

  integer adj(nnode,nnode)
  integer in_seq(nnode)
  integer out_seq(nnode)

  call digraph_adjacency_degree ( adj, nnode, in_seq, out_seq )

  call i4vec2_sort_d ( nnode, out_seq, in_seq )

  return
end
subroutine digraph_adjacency_edge_count ( adj, nnode, nedge )

!*****************************************************************************80
!
!! digraph_adjacency_edge_count() counts the number of edges in a digraph.
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
!    Input, integer ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is an edge from node I to node J.
!
!    Input, integer NNODE, the number of nodes.
!
!    Output, integer NEDGE, the number of edges in the digraph.
!
  implicit none

  integer nnode

  integer adj(nnode,nnode)
  integer nedge

  nedge = sum ( adj(1:nnode,1:nnode) )

  return
end
subroutine digraph_adjacency_eigen ( adj, nnode, neigen, eigenr, eigeni )

!*****************************************************************************80
!
!! digraph_adjacency_eigen(): eigenvalues of a digraph from its adjacency matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is an edge from node I to node J.
!
!    Input, integer NNODE, the number of nodes.
!
!    Output, integer NEIGEN, the number of eigenvalues computed.
!    Normally, this would be equal to NNODE, unless the algorithm failed.
!
!    Output, real ( kind = rk ) EIGENR(NNODE), EIGENI(NNODE), contains the real
!    and imaginary parts of the computed eigenvalues.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nnode

  real ( kind = rk ) a(nnode,nnode)
  integer adj(nnode,nnode)
  real ( kind = rk ) eigeni(nnode)
  real ( kind = rk ) eigenr(nnode)
  integer i
  integer igh
  integer ind(nnode)
  integer info
  integer low
  integer neigen
  real ( kind = rk ) scale(nnode)

  a(1:nnode,1:nnode) = real ( adj(1:nnode,1:nnode), kind = rk )

  call balanc ( nnode, nnode, a, low, igh, scale )

  call elmhes ( nnode, nnode, low, igh, a, ind )

  call hqr ( nnode, nnode, low, igh, a, eigenr, eigeni, info )

  if ( info == 0 ) then
    neigen = nnode
  else
    neigen = nnode - info
    do i = 1, neigen
      eigenr(i) = eigenr(i+info)
      eigeni(i) = eigeni(i+info)
    end do
  end if

  return
end
subroutine digraph_adjacency_example_cycler ( adj )

!*****************************************************************************80
!
!! digraph_adjacency_example_cycler() sets adjacency information for the cycler digraph.
!
!  Diagram:
!  
!           A
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
!    22 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ADJ(9,9), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
  implicit none

  integer adj(9,9)
  integer nnode

  nnode = 9

  adj(1:nnode,1:nnode) = 0

  adj(1,3) = 1
  adj(1,5) = 1

  adj(2,6) = 1
  adj(2,8) = 1

  adj(3,4) = 1
  adj(3,6) = 1
  adj(3,7) = 1

  adj(4,3) = 1

  adj(5,2) = 1

  adj(6,4) = 1
  adj(6,8) = 1

  adj(7,7) = 1
  adj(7,9) = 1

  adj(8,1) = 1

  adj(9,5) = 1
  adj(9,7) = 1

  return
end
subroutine digraph_adjacency_example_octo ( example, adj )

!*****************************************************************************80
!
!! digraph_adjacency_example_octo() sets up an 8 node example digraph.
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
!    There are 12 digraphs to choose from, all on 8 nodes.  There are 7
!    underlying graphs.  The first 5 underlying graphs have degree 3 at 
!    every node.  Graphs 6 and 7 have degree 5 at every node.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer EXAMPLE, should be between 1 and 12, and 
!    indicates which example graph to pick.
!
!    Output, integer ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if nodes I and J are adjacent and 0 otherwise.
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
    nsave = i4_uniform_ab ( 1, 12 )
  else
    example = mod ( example - 1, 12 ) + 1
    nsave = example
  end if

  adj(1:nnode,1:nnode) = 0

  do i = 1, nnode
    j = i + 1
    if ( nnode < j ) then
      j = j - nnode
    end if

    adj(i,j) = 1

  end do
!
!  Underlying graph 1.
!
  if ( nsave == 1 ) then

      adj(1,6) = 1
      adj(2,5) = 1
      adj(3,8) = 1
      adj(4,7) = 1

  else if ( nsave == 2 ) then

      adj(1,6) = 1
      adj(5,2) = 1
      adj(3,8) = 1
      adj(7,4) = 1
!
!  Underlying graph 2.
!  Digraphs 3 and 4 have different indegree/outdegree sequences.
!
  else if ( nsave == 3 ) then

    adj(1,6) = 1
    adj(6,1) = 1
    adj(2,8) = 1
    adj(8,2) = 1
    adj(3,5) = 1
    adj(5,3) = 1
    adj(4,7) = 1
    adj(7,4) = 1

  else if ( nsave == 4 ) then

    adj(1,6) = 1
    adj(2,8) = 1
    adj(3,5) = 1
    adj(4,7) = 1
!
!  Underlying graph 3
!  Digraphs 5 and 6 have the same indegree/outdegree sequences.
!
  else if ( nsave == 5 ) then

    adj(1,5) = 1
    adj(2,6) = 1
    adj(3,7) = 1
    adj(4,8) = 1

  else if ( nsave == 6 ) then

    adj(1:nnode,1:nnode) = 0

    adj(1,8) = 1
    adj(1,5) = 1
    adj(2,1) = 1
    adj(2,3) = 1
    adj(3,4) = 1
    adj(3,7) = 1
    adj(4,5) = 1
    adj(4,8) = 1
    adj(5,6) = 1
    adj(6,2) = 1
    adj(7,6) = 1
    adj(8,7) = 1
!
!  Underlying graph 4
!
  else if ( nsave == 7 ) then

    adj(3,1) = 1
    adj(4,2) = 1
    adj(5,7) = 1
    adj(6,8) = 1

  else if ( nsave == 8 ) then

    adj(3,1) = 1
    adj(4,2) = 1
    adj(5,7) = 1
    adj(8,6) = 1
!
!  Underlying graph 5
!
  else if ( nsave == 9 ) then

    adj(1,4) = 1
    adj(2,6) = 1
    adj(8,3) = 1

    adj(5,7) = 1
    adj(7,5) = 1

  else if ( nsave == 10 ) then

    adj(1,4) = 1
    adj(2,6) = 1
    adj(3,8) = 1

    adj(5,7) = 1
    adj(7,5) = 1
!
!  Underlying graph 6
!
  else if ( nsave == 11 ) then

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

    adj(5,8) = 1
!
!  Underlying graph 7
!
  else if ( nsave == 12 ) then

    adj(1,3) = 1
    adj(1,5) = 1
    adj(1,7) = 1

    adj(2,4) = 1
    adj(2,6) = 1
    adj(2,8) = 1

    adj(3,5) = 1
    adj(3,7) = 1

    adj(4,6) = 1
    adj(4,8) = 1

    adj(5,7) = 1

    adj(6,8) = 1

  end if
!
!  Now permute the graph.
!
  call i4mat_perm_random ( nnode, adj )

  return
end
subroutine digraph_adjacency_example_sixty ( adj )

!*****************************************************************************80
!
!! digraph_adjacency_example_sixty() sets the adjacency matrix for the sixty digraph.
!
!  Discussion:
!
!    The nodes of the digraph are divisors of 60.  There is a link from I to
!    J if divisor I can be divided by divisor J.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Output:
!
!    integer ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
  implicit none

  integer, parameter :: nnode = 12

  integer adj(nnode,nnode)
  integer d(nnode)
  integer i
  integer j

  d(1:12) = (/ 60, 30, 20, 15, 12, 10, 6, 5, 4, 3, 2, 1 /)

  do i = 1, nnode
    do j = 1, nnode
      if ( i == j ) then
        adj(i,j) = 0
      else if ( mod ( d(i), d(j) ) == 0 ) then
        adj(i,j) = 1
      else
        adj(i,j) = 0
      end if
    end do
  end do

  return
end
subroutine digraph_adjacency_ham_cand ( adj, nnode, circuit, k, nstack, &
  stack, maxstack, ncan )

!*****************************************************************************80
!
!! digraph_adjacency_ham_cand(): candidates for the next node in a Hamiltonian circuit.
!
!  Discussion:
!
!    This routine is used in conjunction with I4VEC_BACKTRACK.  
!
!    A Hamiltonian circuit of a digraph is a path that starts at a given node, 
!    visits every node exactly once, and returns to the starting node.
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
!    Input, integer ADJ(NNODE,NNODE).  ADJ(I,J) = 1 if there is
!    an edge from node I to node J, 0 otherwise.
!
!    Input, integer NNODE, the number of nodes in the digraph.
!
!    Input, integer CIRCUIT(NNODE), contains in CIRCUIT(1:K-1) 
!    the partial candidate circuit being constructed.
!
!    Input, integer K, the index of the next node to be determined
!    for the circuit.
!
!    Input/output, integer NSTACK, the current length of the stack.
!
!    Input, integer STACK(MAXSTACK), candidates for positions 
!    1...K-1.
!
!    Input, integer MAXSTACK, the dimension of STACK.
!
!    Input/output, integer NCAN(NNODE), the number of candidates 
!    for each position.  On input, contains values for steps 1 to K-1.  On 
!    output, the value for position K has been determined.
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
          write ( *, '(a)' ) 'DIGRAPH_adjacency_HAM_CAND(): Fatal error!'
          write ( *, '(a)' ) '  MAXSTACK is too small.'
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
 
        if ( adj(i,1) /= 0 ) then
          if ( maxstack <= nstack ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'DIGRAPH_adjacency_HAM_CAND(): Fatal error!'
            write ( *, '(a)' ) '  MAXSTACK is too small.'
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
subroutine digraph_adjacency_ham_next ( adj, nnode, circuit, stack, &
  maxstack, ncan, more )

!*****************************************************************************80
!
!! digraph_adjacency_ham_next() returns the next Hamilton circuit for a digraph.
!
!  Discussion:
!
!    The routine produces all the Hamilton circuits of a digraph, one at a time.
!
!    A Hamiltonian circuit of a digraph is a path that starts at a given
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
!    Input, integer ADJ(NNODE,NNODE), the adjacency matrix of the 
!    digraph.  ADJ(I,J) = 1 if there is an edge from node I to node J, 
!    0 otherwise.
!
!    Input, integer NNODE, the number of nodes in the digraph.
!
!    Input/output, integer CIRCUIT(NNODE).  On the first call to 
!    this routine, the contents of CIRCUIT are irrelevant.  On return, CIRCUIT 
!    contains a list of the nodes that form a cirucit.  On each subsequent 
!    call, the input value of CIRCUIT is used to construct the next solution,
!    so the user should not alter the contents of CIRCUIT during a computation.
!
!    Workspace, integer STACK(MAXSTACK).  Candidates for the positions in
!    the circuit.
!
!    Input, integer MAXSTACK, the dimension of STACK.
!
!    Workspace, integer NCAN(NNODE), a count of the number of candidates for 
!    each step.
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

      call digraph_adjacency_ham_cand ( adj, nnode, circuit, k, nstack, &
        stack, maxstack, ncan )

    else

      more = .false.
      exit

    end if

  end do
 
  return
end
subroutine digraph_adjacency_ham_next_brute ( adj, nnode, circuit, iset )

!*****************************************************************************80
!
!! digraph_adjacency_ham_next_brute() finds the next Hamiltonian circuit in a digraph.
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
!    01 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
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
 
  do
 
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
    ipos = 0
    do i = 1, nnode - 1
      if ( adj(circuit(i),circuit(i+1)) == 0 ) then
        ipos = i
        exit
      end if
    end do

    if ( ipos /= 0 ) then
      cycle
    end if
!
!  If the circuit connects all the nodes, we only have to check whether
!  the last node connects back to the first one.
!
    if ( adj(circuit(nnode),circuit(1)) /= 0 ) then
      exit
    end if

    ipos = nnode - 1

  end do

  iset = iset + 1

  return
end
subroutine digraph_adjacency_ham_path_next_brute ( adj, nnode, path, iset )

!*****************************************************************************80
!
!! digraph_adjacency_ham_path_next_brute(): next path in digraph that visits all nodes.
!
!  Discussion:
!
!    The path is not required to be a circuit.  That is, there is no requirement
!    that there be an edge from the last node visited back to the first one.
!
!    This is a brute force algorithm, and not suitable for large problems.
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
!    Input, integer ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Input, integer NNODE, the number of nodes.
!
!    Input/output, integer PATH(NNODE).
!
!    On input, if ISET = 0, then PATH is not presumed to contain any
!    information.  If ISET is nonzero, then PATH contains the
!    path computed on the previous call.
!
!    On output, PATH contains the path computed by this call.
!
!    Input/output, integer ISET.
!
!    On input, a 0 value means this is the first call for this
!    graph.  Any other value means this is a repeated call for more paths.
!
!    On output, a 0 value means that no more paths could be computed.
!    Otherwise, ISET is incremented by one on each call.
!
  implicit none

  integer nnode

  integer adj(nnode,nnode)
  integer i
  integer ipos
  integer iset
  integer path(nnode)
!
!  If ISET is 0, this is a starting call, and we set PATH
!  to the lexically first path to check.
!
!  Otherwise, we set PATH to the next permutation.
!
  if ( iset == 0 ) then
    ipos = 0
    path(1:nnode) = 0
  else
    ipos = nnode - 1
  end if
 
  do
 
    call perm_inc ( path, ipos, nnode )
 
    if ( ipos == 0 ) then
      iset = 0
      path(1:nnode) = 0
      return
    end if
!
!  Check whether the entries of PATH actually form a path.
!
    ipos = 0
    do i = 1, nnode - 1
      if ( adj(path(i),path(i+1)) == 0 ) then
        ipos = i
        exit
      end if
    end do

    if ( ipos == 0 ) then
      exit
    end if

  end do 

  iset = iset + 1
 
  return
end
subroutine digraph_adjacency_is_edge_connected ( adj, nnode, result )

!*****************************************************************************80
!
!! digraph_adjacency_is_edge_connected() determines if a digraph is edgewise connected.
!
!  Discussion:
!
!    A digraph is edgewise connected if from any edge it is possible to reach
!    any other edge.  An edgewise connected digraph may include isolated nodes.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Input, integer NNODE, the number of nodes.
!
!    Output, integer RESULT.
!    0, the digraph is not edgewise connected.
!    1, the digraph is edgewise connected.
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
!  FOUND(I) is 1 if edge I has been reached.
!  LIST(I) contains a list of the nodes as they are reached.
!
  list(1:nnode) = 0
  found(1:nnode) = 0
!
!  Find an edge.
!
  ilo = 1
  ihi = 0

  do i = 1, nnode
    do j = 1, nnode

      if ( 0 < adj(i,j) ) then

        adj(i,j) = - adj(i,j)
        ihi = ihi + 1
        list(ihi) = i
        found(i) = 1

        if ( i /= j ) then
          ihi = ihi + 1
          list(ihi) = j
          found(j) = 1
        end if

        exit

      end if

    end do

    if ( 0 < ihi ) then
      exit
    end if

  end do
!
!  A digraph with NO edges is edgewise connected!
!
  if ( ihi == 0 ) then
    result = 1
    return
  end if
!
!  From the batch of edge nodes found last time, LIST(ILO:IHI),
!  look for unfound neighbors, and store their indices in LIST(JLO:JHI).
!
  do

    jlo = ihi + 1
    jhi = ihi

    do ii = ilo, ihi

      i = list(ii)

      do j = 1, nnode

        if ( 0 < adj(i,j) ) then

          adj(i,j) = - adj(i,j)

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
!  If any edges were unvisited, then the digraph is not edgewise connected.
!
  result = 1

  do i = 1, nnode
    do j = 1, nnode
      if ( 0 < adj(i,j) ) then
        result = 0
      end if
    end do
  end do
!
!  Restore the positive sign of ADJ.
!
  adj(1:nnode,1:nnode) = abs ( adj(1:nnode,1:nnode) )

  return
end
subroutine digraph_adjacency_is_eulerian ( adj, nnode, result )

!*****************************************************************************80
!
!! digraph_adjacency_is_eulerian() determines if a digraph is Eulerian.
!
!  Discussion:
!
!    A digraph is path-Eulerian if there exists a path through the digraph
!    which uses every edge once.
!
!    A digraph is circuit-Eulerian if there exists a path through the digraph
!    which uses every edge once, and which starts and ends on the same node.
!
!    Note that it is NOT necessary for the path or circuit to pass through
!    every node; simply that all the edges can be used exactly once to
!    make a connected path.  This means an Eulerian digraph can have isolated
!    nodes, for instance.
!
!    A digraph is path-Eulerian if and only if it is edge-connected, and 
!    for all but two nodes, the indegree and outdegree are equal, and
!    for those two nodes, the indegree and outdegree, if different, differ
!    by 1.
!
!    A digraph is circuit-Eulerian if and only if it is edge connected and
!    for every node the indegree equals the outdegree.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    28 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer adj(nnode,nnode), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Input, integer NNODE, the number of nodes.
!
!    Output, integer RESULT.
!    0, the digraph is not Eulerian.
!    1, the digraph is path-Eulerian.
!    2, the digraph is circuit-Eulerian.
!
  implicit none

  integer nnode

  integer adj(nnode,nnode)
  integer i
  integer indegree(nnode)
  integer ndiff
  integer outdegree(nnode)
  integer result
!
!  First check that the digraph is edgewise connected.
!
  call digraph_adjacency_is_edge_connected ( adj, nnode, result )

  if ( result == 0 ) then
    return
  end if
!
!  Now look at node degree.
!
  call digraph_adjacency_degree ( adj, nnode, indegree, outdegree )

  ndiff = 0

  do i = 1, nnode

    if ( indegree(i) /= outdegree(i) ) then

      ndiff = ndiff + 1

      if ( 2 < ndiff ) then
        result = 0
        return
      end if

      if ( 1 < abs ( indegree(i) - outdegree(i) ) ) then
        result = 0
        return
      end if

    end if

  end do

  if ( ndiff == 0 ) then
    result = 2
  else
    result = 1
  end if

  return
end
subroutine digraph_adjacency_is_strongly_connected ( adj, nnode, result )

!*****************************************************************************80
!
!! digraph_adjacency_is_strongly_connected(): is a digraph strongly connected?
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer adj(nnode,nnode), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Input, integer NNODE, the number of nodes.
!
!    Output, integer RESULT,
!    0, the digraph is not strongly connected;
!    1, the digraph is strongly connected.
!
  implicit none

  integer nnode

  integer adj(nnode,nnode)
  integer dad(nnode)
  integer iorder
  integer lowlink(nnode)
  integer mark(nnode)
  integer ncomp
  integer nstack
  integer order(nnode)
  integer point(nnode)
  integer result
  integer stack(nnode)
  integer v
  integer w
  integer x
!
!  Initialization.
!
  dad(1:nnode) = 0
  order(1:nnode) = 0
  lowlink(1:nnode) = 0
  mark(1:nnode) = 0
  point(1:nnode) = 0

  iorder = 0
  nstack = 0
  ncomp = 0
!
!  Select any node V not stored in the stack, that is, with MARK(V) = 0.
!
  do

    v = 0

    do

      v = v + 1

      if ( nnode < v ) then

        adj(1:nnode,1:nnode) = abs ( adj(1:nnode,1:nnode) )

        if ( 1 < ncomp ) then
          result = 0
        else
          result = 1
        end if

        return
      end if

      if ( mark(v) /= 1 ) then
        exit
      end if

    end do

    iorder = iorder + 1

    order(v) = iorder
    lowlink(v) = iorder
    mark(v) = 1

    nstack = nstack + 1
    stack(nstack) = v
    point(v) = 1

30  continue
!
!  Consider each node W.
!
    do w = 1, nnode
!
!  Is there an edge (V,W) and has it not been examined yet?
!
      if ( 0 < adj(v,w) ) then

        adj(v,w) = - adj(v,w)
!
!  Is the node on the other end of the edge undiscovered yet?
!
        if ( mark(w) == 0 ) then

          iorder = iorder + 1
          order(w) = iorder
          lowlink(w) = iorder
          dad(w) = v
          mark(w) = 1

          nstack = nstack + 1
          stack(nstack) = w
          point(w) = 1

          v = w

        else if ( mark(w) == 1 ) then

          if ( order(w) < order(v) .and. point(w) == 1 ) then
            lowlink(v) = min ( lowlink(v), order(w) )
          end if

        end if

        go to 30

      end if

    end do

    if ( lowlink(v) == order(v) ) then

      ncomp = ncomp + 1

      do

        if ( nstack <= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'digraph_adjacency_is_strongly_connected(): Fatal error!'
          write ( *, '(a)' ) '  Illegal stack reference.'
          stop 1
        end if

        x = stack(nstack)
        nstack = nstack - 1

        point(x) = 0

        if ( x == v ) then
          exit
        end if

      end do

    end if

    if ( dad(v) /= 0 ) then
      lowlink(dad(v)) = min ( lowlink(dad(v)), lowlink(v) )
      v = dad(v)
      go to 30
    end if

  end do

  return
end
subroutine digraph_adjacency_is_tournament ( adj, nnode, result )

!*****************************************************************************80
!
!! digraph_adjacency_is_tournament() determines if a digraph is a tournament.
!
!  Discussion:
!
!    A digraph is a tournament if every pair of distinct nodes is connected by
!    exactly one directed edge.
!
!    We interpret the directed edge "i --> j" as "i beat j in the tournament".
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer adj(nnode,nnode), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    integer NNODE, the number of nodes.
!
!  Output:
!
!    integer RESULT.
!    0, the digraph is not a tournament.
!    1, the digraph is a tournament.
!
  implicit none

  integer nnode

  integer adj(nnode,nnode)
  integer i
  integer j
  integer result

  result = 0
!
!  No self links.
!
  do i = 1, nnode
    if ( adj(i,i) /= 0 ) then
      return
    end if
  end do
!
!  Distinct I and J must have exactly one connection.
!
  do i = 1, nnode
    do j = i + 1, nnode
      if ( .not. ( adj(i,j) == 0 .and. adj(j,i) == 1 ) .and. &
           .not. ( adj(i,j) == 1 .and. adj(j,i) == 0 ) ) then
        return
      end if
    end do
  end do

  result = 1
 
  return
end
subroutine digraph_adjacency_is_transitive ( adj, nnode, result )

!*****************************************************************************80
!
!! digraph_adjacency_is_transitive() determines if a digraph is transitive.
!
!  Discussion:
!
!    A digraph is transitive if whenever there's a long way between two
!    nodes, there's an immediate way.  Formally:
!
!      ADJ(I,J) and ADJ(J,K) nonzero imply ADJ(I,K) nonzero.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer adj(nnode,nnode), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    integer NNODE, the number of nodes.
!
!  Output:
!
!    integer RESULT.
!    0, the digraph is not transitive.
!    1, the digraph is transitive.
!
  implicit none

  integer nnode

  integer adj(nnode,nnode)
  integer i
  integer j
  integer k
  integer result

  result = 0

  do i = 1, nnode
    do j = 1, nnode
      if ( adj(i,j) /= 0 ) then
        do k = 1, nnode
          if ( adj(j,k) /= 0 ) then
            if ( adj(i,k) == 0 ) then
              return
            end if
          end if
        end do
      end if
    end do
  end do

  result = 1

  return
end
subroutine digraph_adjacency_is_weakly_connected ( adj, nnode, result )

!*****************************************************************************80
!
!! digraph_adjacency_is_weakly_connected() determines if a digraph is weakly connected.
!
!  Discussion:
!
!    A digraph is weakly connected if the underlying graph is node connected.
!    In other words, if a graph is constructed from the digraph by replacing
!    every directed edge by an undirected edge, and the it is possible to
!    travel from any node to any other node, then the digraph is weakly
!    connected.
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
!    Input, integer adj(nnode,nnode), the adjacency matrix for 
!    the digraph.  ADJ(I,J) is nonzero if there is an edge from node I 
!    to node J.
!
!    Input, integer NNODE, the number of nodes.
!
!    Output, integer RESULT.
!    0, the digraph is not weakly connected.
!    1, the digraph is weakly connected.
!
  implicit none

  integer nnode

  integer adj(nnode,nnode)
  integer result

  call graph_adjacency_is_node_connected ( adj, nnode, result )

  return
end
subroutine digraph_adjacency_print ( adj, nnode, title )

!*****************************************************************************80
!
!! digraph_adjacency_print() prints out an adjacency matrix for a digraph.
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
subroutine digraph_adjacency_random ( nnode, nedge, adj )

!*****************************************************************************80
!
!! digraph_adjacency_random() generates a random digraph.
!
!  Discussion:
!
!    A digraph is a directed graph.
!
!    The user specifies the number of nodes and edges in the digraph.
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
!    Input, integer NEDGE, the number of edges, which must be 
!    between 0 and NNODE*(NNODE-1).
!
!    Output, integer adj(nnode,nnode), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
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

  if ( nnode <= 0  ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'digraph_adjacency_random(): Fatal error!'
    write ( *, '(a,i8)' ) '  NNODE = ', nnode
    write ( *, '(a)' ) '  but NNODE must be at least 1.'
    stop 1
  end if

  maxedge = nnode * ( nnode - 1 )

  if ( nedge < 0 .or. maxedge < nedge ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'digraph_random(): Fatal error!'
    write ( *, '(a,i8)' ) '  NEDGE = ', nedge
    write ( *, '(a)' ) '  but NEDGE must be at least 0, and '
    write ( *, '(a,i8)' ) '  no more than ', maxedge
    stop 1
  end if

  adj(1:nnode,1:nnode) = 0
!
!  Pick a random NEDGE subset of NNODE*(NNODE-1).
!
  call ksub_random ( maxedge, nedge, iwork )
!
!  The usable spots in the matrix are numbered as follows:
!
!   *    1    2   3  ...      n-2        n-1
!   n    *   n+1 n+2 ...     2n-1      2(n-1)
!  2n-1  2n   *  ... ... ........  ..........
!  .... ...  ... ... ...     *     (n-1)(n-1)
!  .... ...  ... ... ...   n(n-1)       *
!
  k = 0
  l = 1
  do i = 1, nnode
    do j = 1, nnode

      if ( i /= j ) then

        k = k + 1
        if ( l <= nedge ) then

          if ( k == iwork(l) ) then
            adj(i,j) = 1
            l = l + 1
          end if

        end if

      end if

    end do
  end do

  return
end
subroutine digraph_adjacency_reduce ( adj, nnode )

!*****************************************************************************80
!
!! digraph_adjacency_reduce() generates a transitive reduction of a digraph.
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
!    15 April 1999
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
!  If you can get from J to I and I to K, you don't need a direct
!  edge from J to K.
!
  do i = 1, nnode
    do j = 1, nnode
      if ( adj(j,i) /= 0 ) then
        do k = 1, nnode
          if ( adj(i,k) /= 0 ) then
            adj(j,k) = 0
          end if
        end do
      end if
    end do
  end do

  return
end
subroutine digraph_adjacency_to_digraph_arc ( adj, nnode, maxedge, nedge, &
  inode, jnode )

!*****************************************************************************80
!
!! digraph_adjacency_to_digraph_arc() converts digraph from adjacency to arc list form.
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
!    Input, integer adj(nnode,nnode), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, integer MAXEDGE, the maximum number of edges.
!
!    Output, integer NEDGE, the number of edges.
!
!    Output, integer INODE(MAXEDGE), JNODE(MAXEDGE), the arc list 
!    of the digraph.
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
    do i = 1, nnode
      if ( adj(i,j) /= 0 ) then
        nedge = nedge + 1
        if ( nedge <= maxedge ) then
          inode(nedge) = i
          jnode(nedge) = j
        else
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'digraph_adjacency_to_digraph_arc(): Fatal error!'
          write ( *, '(a)' ) '  MAXEDGE exceeded.'
          stop 1
        end if
      end if
    end do
  end do

  return
end
subroutine digraph_adjacency_to_incidence ( adj, nnode, maxarc, narc, inc )

!*****************************************************************************80
!
!! digraph_adjacency_to_incidence() converts adjacency to incidence matrix.
!
!  Discussion:
!
!    INC(node,arc) = 0 if NODE is not the beginning or end of ARC, or
!                       if ARC is a loop;
!                     1 if NODE is the beginning of ARC;
!                    -1 if NODE is the end of ARC.
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 July 2000
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
!    integer MAXARC, the maximum number of arcs.
!
!  Output:
!
!    integer NARC, the number of arcs.
!
!    integer INC(NNODE,MAXARC), the incidence matrix.
!
  implicit none

  integer maxarc
  integer nnode

  integer adj(nnode,nnode)
  integer i
  integer inc(nnode,maxarc)
  integer j
  integer narc

  narc = 0

  do j = 1, nnode
    do i = 1, nnode

      if ( i == j ) then

      else if ( adj(i,j) /= 0 ) then
        narc = narc + 1
        if ( narc <= maxarc ) then
          inc(i,narc) = 1
          inc(j,narc) = -1
        else
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'digraph_adjacency_to_incidence(): Fatal error!'
          write ( *, '(a)' ) '  MAXARC exceeded.'
          stop 1
        end if
      end if
    end do
  end do

  return
end
subroutine digraph_adjacency_top_sort ( adj, nnode, dad, visit, node_list )

!*****************************************************************************80
!
!! digraph_adjacency_top_sort(): reverse topological sort of a directed acyclic graph.
!
!  Discussion:
!
!    The routine performs a depth first search of the DAG and returns:
!
!    * a list of the order in which the nodes were visited;
!    * a list of the parents of each node in the search trees;
!    * a list of the nodes, in a reverse topological order.
!
!    In a reverse topological sorting of the nodes of a directed
!    acyclic graph, nodes are listed "lowest" first.  That is,
!    if node A precedes node B in the list, then there may or may
!    not be an edge or indirect path from B to A, but there
!    is neither an edge or indirect path from A to B.
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
!    Robert Sedgewick,
!    Algorithms,
!    Addison Wesley, 1983, page 426.
!
!  Parameters:
!
!    Input, integer adj(nnode,nnode), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Input, integer NNODE, the number of nodes.
!
!    Output, integer DAD(NNODE), the father array for the depth
!    first search trees.  DAD(I) = 0 means that node I is the root of 
!    one of the trees.  DAD(I) = J means that the search descended
!    from node J to node I.
!
!    Output, integer VISIT(NNODE), the order in which the nodes
!    were visited, from 1 to NNODE.
!
!    Output, integer NODE_LIST(NNODE), a list of the nodes, in
!    reverse topological order.
!
  implicit none

  integer nnode

  integer adj(nnode,nnode)
  integer dad(nnode)
  integer daddy
  integer i
  integer j
  integer maxstack
  integer nsort
  integer nstack
  integer node_list(nnode)
  integer rank
  integer stack(2*(nnode-1))
  integer visit(nnode)

  dad(1:nnode) = 0
  maxstack = 2 * ( nnode - 1 )
  visit(1:nnode) = 0
  node_list(1:nnode) = 0

  rank = 0
  nsort = 0

  do i = 1, nnode
!
!  Find the next unused node and begin a new search tree.
!
    if ( visit(i) == 0 ) then

      daddy = i
      dad(daddy) = 0
      nstack = 0
!
!  Visit node DAD.
!
10    continue

      rank = rank + 1
      visit(daddy) = rank
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

        if ( adj(daddy,j) /= 0 .and. visit(j) == 0 ) then

          if ( nstack+2 <= maxstack ) then
            dad(j) = daddy
            stack(nstack+1) = daddy
            stack(nstack+2) = j
            nstack = nstack + 2
            daddy = j
            go to 10
          else
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'DIGRAPH_adjacency_TOP_SORT(): Fatal error!'
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

        nsort = nsort + 1
        node_list(nsort) = daddy

        daddy = stack(nstack-1)
        j = stack(nstack)
        nstack = nstack - 2
        go to 20
!
!  If no more nodes to consider in this tree, bail out.
!
      else

        nsort = nsort + 1
        node_list(nsort) = daddy

        nstack = 0

      end if

    end if

  end do

  return
end
subroutine digraph_adjacency_tournament_random ( nnode, adj )

!*****************************************************************************80
!
!! digraph_adjacency_tournament_random() generates a random tournament digraph.
!
!  Discussion:
!
!    Definition: A tournament is a directed graph in which every pair 
!    of nodes are joined by exactly one directed edge.
!
!    The user specifies the number of nodes in the digraph.  The number of
!    edges will be (NNODE*(NNODE-1))/2.
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
!    Output, integer adj(nnode,nnode), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
  implicit none

  integer nnode

  integer adj(nnode,nnode)
  integer i
  integer i4_uniform_ab
  integer j
  integer k

  if ( nnode <= 0  ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIGRAPH_adjacency_TOURNAMENT_RANDOM(): Fatal error!'
    write ( *, '(a,i8)' ) '  NNODE = ', nnode
    write ( *, '(a)' ) '  but NNODE must be at least 1.'
    stop 1
  end if

  adj(1:nnode,1:nnode) = 0

  do i = 1, nnode
    do j = i + 1, nnode

      k = i4_uniform_ab ( 1, 2 )

      if ( k == 1 ) then
        adj(i,j) = 1
      else
        adj(j,i) = 1
      end if

    end do
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
subroutine digraph_incidence_print ( nnode, narc, inc, title )

!*****************************************************************************80
!
!! digraph_incidence_print() prints the incidence matrix of a digraph.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, integer NARC, the number of arcs.
!
!    Input, integer INC(NNODE,NARC), the incidence matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer narc
  integer nnode

  integer i
  integer inc(nnode,narc)
  character ( len = * ) title

  if ( len_trim ( title ) /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  do i = 1, nnode
    write ( *, '(20i3)' ) inc(i,1:narc)
  end do

  return
end
subroutine elmhes ( nm, n, low, igh, a, ind )

!*****************************************************************************80
!
!! elmhes() transforms a real general matrix to upper Hessenberg form.
!
!  Discussion:
!
!    Given a real general matrix, this subroutine reduces a submatrix
!    situated in rows and columns LOW through IGH to upper Hessenberg
!    form by stabilized elementary similarity transformations.
!
!  Reference:
!
!    Martin, James Wilkinson,
!    ELMHES,
!    Numerische Mathematik,
!    Volume 12, pages 349-368, 1968.
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
!    Input, integer N, the order of the matrix.
!
!    Input, integer LOW, IGH, are determined by the balancing
!    routine BALANC.  If BALANC has not been used, set LOW = 1, IGH = N.
!
!    Input/output, real ( kind = rk ) A(NM,N).  On input, the matrix to be
!    reduced.  On output, the Hessenberg matrix.  The multipliers
!    which were used in the reduction are stored in the
!    remaining triangle under the Hessenberg matrix.
!
!    Output, integer IND(N), contains information on the rows and
!    columns interchanged in the reduction.  Only elements LOW through IGH are
!    used.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer igh
  integer n
  integer nm

  real ( kind = rk ) a(nm,n)
  integer i
  integer ind(igh)
  integer j
  integer la
  integer low
  integer m
  integer mm1
  real ( kind = rk ) x
  real ( kind = rk ) y

  la = igh - 1

  do m = low + 1, la

    mm1 = m - 1
    x = 0.0D+00
    i = m

    do j = m, igh
      if ( abs ( x ) < abs ( a(j,mm1) ) ) then
        x = a(j,mm1)
        i = j
      end if
    end do

    ind(m) = i
!
!  Interchange rows and columns of the matrix.
!
    if ( i /= m ) then

      do j = mm1, n
        call r8_swap ( a(i,j), a(m,j) )
      end do

      do j = 1, igh
        call r8_swap ( a(j,i), a(j,m) )
      end do

    end if

    if ( x /= 0.0D+00 ) then

      do i = m+1, igh

        y = a(i,mm1)

        if ( y /= 0.0D+00 ) then

          y = y / x
          a(i,mm1) = y

          do j = m, n
            a(i,j) = a(i,j) - y * a(m,j)
          end do

          do j = 1, igh
            a(j,m) = a(j,m) + y * a(j,i)
          end do

        end if

      end do

    end if

  end do

  return
end
subroutine graph_adjacency_is_node_connected ( adj, nnode, result )

!*****************************************************************************80
!
!! graph_adjacency_is_node_connected() determines if a graph is nodewise connected.
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
subroutine hqr ( nm, n, low, igh, h, wr, wi, ierr )

!*****************************************************************************80
!
!! hqr() computes all eigenvalues of a real upper Hessenberg matrix.
!
!  Discussion:
!
!    This subroutine finds the eigenvalues of a real
!    upper Hessenberg matrix by the QR method.
!
!  Reference:
!
!    Martin, Peters, James Wilkinson,
!    HQR,
!    Numerische Mathematik,
!    Volume 14, pages 219-231, 1970.
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
!    Input, integer NM, the leading dimension of H, which must
!    be at least N.
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer LOW, IGH, two integers determined by the
!    routine BALANC.  If BALANC is not used, set LOW=1, IGH=N.
!
!    Input/output, real ( kind = rk ) H(NM,N), the N by N upper Hessenberg
!    matrix.  Information about the transformations used in the reduction to
!    Hessenberg form by ELMHES or ORTHES, if performed, is stored
!    in the remaining triangle under the Hessenberg matrix.
!    On output, the information in H has been destroyed.
!
!    Output, real ( kind = rk ) WR(N), WI(N), the real and imaginary parts of the
!    eigenvalues.  The eigenvalues are unordered, except that complex
!    conjugate pairs of values appear consecutively, with the eigenvalue
!    having positive imaginary part listed first.  If an error exit
!    occurred, then the eigenvalues should be correct for indices
!    IERR+1 through N.
!
!    Output, integer IERR, error flag.
!    0, no error.
!    J, the limit of 30*N iterations was reached while searching for
!       the J-th eigenvalue.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nm

  integer en
  integer enm2
  real ( kind = rk ) h(nm,n)
  integer i
  integer ierr
  integer igh
  integer itn
  integer its
  integer j
  integer k
  integer l
  integer ll
  integer low
  integer m
  integer mm
  integer na
  real ( kind = rk ) norm
  logical notlas
  real ( kind = rk ) p
  real ( kind = rk ) q
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) t
  real ( kind = rk ) tst1
  real ( kind = rk ) tst2
  real ( kind = rk ) w
  real ( kind = rk ) wi(n)
  real ( kind = rk ) wr(n)
  real ( kind = rk ) x
  real ( kind = rk ) y
  real ( kind = rk ) zz

  ierr = 0
  norm = 0.0D+00
  k = 1
!
!  Store roots isolated by BALANC and compute matrix norm.
!
  do i = 1, n

    do j = k, n
      norm = norm + abs ( h(i,j) )
    end do

    k = i
    if (i < low .or. igh < i ) then
      wr(i) = h(i,i)
      wi(i) = 0.0D+00
    end if

  end do

  en = igh
  t = 0.0D+00
  itn = 30 * n
!
!  Search for next eigenvalues.
!
60 continue

  if ( en < low ) then
    return
  end if

  its = 0
  na = en - 1
  enm2 = na - 1
!
!  Look for a single small sub-diagonal element.
!
70 continue

  do ll = low, en
    l = en + low - ll
    if ( l == low ) then
      exit
    end if
    s = abs ( h(l-1,l-1) ) + abs ( h(l,l) )
    if ( s == 0.0D+00 ) then
      s = norm
    end if
    tst1 = s
    tst2 = tst1 + abs ( h(l,l-1))
    if ( tst2 == tst1 ) then
      exit
    end if
  end do
!
!  Form shift.
!
  x = h(en,en)
  if ( l == en ) go to 270
  y = h(na,na)
  w = h(en,na) * h(na,en)
  if ( l == na ) go to 280

  if ( itn == 0 ) then
    ierr = en
    return
  end if
!
!  Form an exceptional shift.
!
  if ( its == 10 .or. its == 20 ) then

    t = t + x

    do i = low, en
      h(i,i) = h(i,i) - x
    end do

    s = abs ( h(en,na) ) + abs ( h(na,enm2) )
    x = 0.75D+00 * s
    y = x
    w = -0.4375D+00 * s * s

  end if

  its = its + 1
  itn = itn - 1
!
!  Look for two consecutive small sub-diagonal elements.
!
  do mm = l, enm2

    m = enm2 + l - mm
    zz = h(m,m)
    r = x - zz
    s = y - zz
    p = ( r * s - w ) / h(m+1,m) + h(m,m+1)
    q = h(m+1,m+1) - zz - r - s
    r = h(m+2,m+1)
    s = abs ( p ) + abs ( q ) + abs ( r )
    p = p / s
    q = q / s
    r = r / s

    if ( m == l ) then
      exit
    end if

    tst1 = abs ( p ) * ( abs ( h(m-1,m-1) ) + abs ( zz ) + abs ( h(m+1,m+1) ) )
    tst2 = tst1 + abs ( h(m,m-1) ) * ( abs ( q ) + abs ( r ) )

    if ( tst2 == tst1 ) then
      exit
    end if

  end do

  do i = m+2, en
    h(i,i-2) = 0.0D+00
    if ( i /= m+2 ) then
      h(i,i-3) = 0.0D+00
    end if
  end do
!
!  Double QR step involving rows l to EN and columns M to EN.
!
  do k = m, na

     notlas = k /= na

     if ( k == m ) go to 170

     p = h(k,k-1)
     q = h(k+1,k-1)
     if ( notlas ) then
       r = h(k+2,k-1)
     else
       r = 0.0D+00
     end if
     x = abs ( p ) + abs ( q ) + abs ( r )
     if ( x == 0.0D+00 ) go to 260
     p = p / x
     q = q / x
     r = r / x

170  continue

     s = sign ( sqrt ( p**2 + q**2 + r**2 ), p )

     if ( k /= m ) then
       h(k,k-1) = - s * x
     else if ( l /= m ) then
       h(k,k-1) = - h(k,k-1)
     end if

     p = p + s
     x = p / s
     y = q / s
     zz = r / s
     q = q / p
     r = r / p
     if ( notlas ) go to 225
!
!  Row modification.
!
     do j = k, n
       p = h(k,j) + q * h(k+1,j)
       h(k,j) = h(k,j) - p * x
       h(k+1,j) = h(k+1,j) - p * y
     end do

     j = min ( en, k+3 )
!
!  Column modification.
!
     do i = 1, j
       p = x * h(i,k) + y * h(i,k+1)
       h(i,k) = h(i,k) - p
       h(i,k+1) = h(i,k+1) - p * q
     end do

     go to 255

225  continue
!
!  Row modification.
!
     do j = k, n
       p = h(k,j) + q * h(k+1,j) + r * h(k+2,j)
       h(k,j) = h(k,j) - p * x
       h(k+1,j) = h(k+1,j) - p * y
       h(k+2,j) = h(k+2,j) - p * zz
     end do

     j = min ( en, k+3 )
!
!  Column modification.
!
     do i = 1, j
       p = x * h(i,k) + y * h(i,k+1) + zz * h(i,k+2)
       h(i,k) = h(i,k) - p
       h(i,k+1) = h(i,k+1) - p * q
       h(i,k+2) = h(i,k+2) - p * r
     end do

255 continue

260 continue

  end do

  go to 70
!
!  One root found.
!
270 continue

  wr(en) = x + t
  wi(en) = 0.0D+00
  en = na
  go to 60
!
!  Two roots found.
!
280 continue

  p = ( y - x ) / 2.0D+00
  q = p * p + w
  zz = sqrt ( abs ( q ) )
  x = x + t
!
!  Real root, or complex pair.
!
  if ( 0.0D+00 <= q ) then

    zz = p + sign ( zz, p )
    wr(na) = x + zz
    if ( zz == 0.0D+00 ) then
      wr(en) = wr(na)
    else
      wr(en) = x - w / zz
    end if
    wi(na) = 0.0D+00
    wi(en) = 0.0D+00

  else

    wr(na) = x + p
    wr(en) = x + p
    wi(na) = zz
    wi(en) = -zz

  end if

  en = enm2
  go to 60
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
subroutine i4mat_perm ( matrix, n, p )

!*****************************************************************************80
!
!! i4mat_perm() permutes the rows and columns of a square I4MAT.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 July 2000
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    This version by John Burkardt.
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
!    Input/output, integer MATRIX(N,N).
!    On input, the matrix to be permuted.
!    On output, the permuted matrix.
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer P(N), the permutation.  P(I) is the new number
!    of row and column I.
!
  implicit none

  integer n

  integer i
  integer i1
  integer is
  integer it
  integer j
  integer j1
  integer j2
  integer k
  integer lc
  integer matrix(n,n)
  integer nc
  integer p(n)

  call perm_cycle ( n, p, is, nc, 1 )

  do i = 1, n

    i1 = - p(i)

    if ( 0 < i1 ) then

      lc = 0

      do

        i1 = p(i1)
        lc = lc + 1

        if ( i1 <= 0 ) then
          exit
        end if

      end do

      i1 = i

      do j = 1, n

        if ( p(j) <= 0 ) then

          j2 = j
          k = lc

          do

            j1 = j2
            it = matrix(i1,j1)

            do

              i1 = abs ( p(i1) )
              j1 = abs ( p(j1) )

              call i4_swap ( matrix(i1,j1), it )

              if ( j1 /= j2 ) then
                cycle
              end if

              k = k - 1

              if ( i1 == i ) then
                exit
              end if

            end do

            j2 = abs ( p(j2) )

            if ( k == 0 ) then
              exit
            end if

          end do

        end if

      end do

    end if

  end do
!
!  Restore the positive signs of the data.
!
  p(1:n) = abs ( p(1:n) )

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
!    This version by John Burkardt.
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
subroutine i4vec2_compare ( n, i4vec, jvec, i, j, isgn )

!*****************************************************************************80
!
!! i4vec2_compare() compares pairs of integers stored in two vectors.
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
!    Input, integer N, the number of data items.
!
!    Input, integer I4VEC(N), JVEC(N), contain the two components 
!    of each item.
!
!    Input, integer I, J, the items to be compared.
!
!    Output, integer ISGN, the results of the comparison:
!    -1, item I is less than item J,
!     0, item I is equal to item J,
!    +1, item I is greater than item J.
!
  implicit none

  integer n

  integer i
  integer isgn
  integer i4vec(n)
  integer j
  integer jvec(n)

  isgn = 0

  if ( i4vec(i) < i4vec(j) ) then
    isgn = -1
  else if ( i4vec(i) == i4vec(j) ) then
    if ( jvec(i) < jvec(j) ) then
      isgn = -1
    else if ( jvec(i) < jvec(j) ) then
      isgn = 0
    else if ( jvec(j) < jvec(i) ) then
      isgn = +1
    end if
  else if ( i4vec(j) < i4vec(i) ) then
    isgn = +1
  end if

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
subroutine i4vec2_sort_d ( n, a1, a2 )

!*****************************************************************************80
!
!! i4vec2_sort_d() descending sorts a vector of pairs of integers.
!
!  Discussion:
!
!    Each item to be sorted is a pair of integers (I,J), with the I
!    and J values stored in separate vectors A1 and A2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of items of data.
!
!    Input/output, integer A1(N), A2(N), the data to be sorted..
!
  implicit none

  integer n

  integer a1(n)
  integer a2(n)
  integer i
  integer indx
  integer isgn
  integer j
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      call i4_swap ( a1(i), a1(j) )
      call i4_swap ( a2(i), a2(j) )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4vec2_compare ( n, a1, a2, i, j, isgn )
      isgn = - isgn

    else if ( indx == 0 ) then

      exit

    end if

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
!    This version by John Burkardt.
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
subroutine perm_cycle ( n, isig, isgn, ncycle, iopt )

!*****************************************************************************80
!
!! perm_cycle() analyzes a permutation.
!
!  Discussion:
!
!    The routine will count cycles, find the sign of a permutation,
!    and tag a permutation.
!
!  Example:
!
!    Input:
!
!      N = 9
!      IOPT = 1
!      ISIG = 2, 3, 9, 6, 7, 8, 5, 4, 1
!
!    Output:
!
!      NCYCLE = 3
!      ISGN = +1
!      ISIG = -2, 3, 9, -6, -7, 8, 5, 4, 1
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
!    This version by John Burkardt.
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
!    Input, integer N, the number of objects being permuted.
!
!    Input/output, integer ISIG(N).  On input, ISIG describes a
!    permutation, in the sense that entry I is to be moved to ISIG(I).
!    If IOPT = 0, then ISIG will not be changed by this routine.
!    If IOPT = 1, then on output, ISIG will be "tagged".  That is,
!    one element of every cycle in ISIG will be negated.  In this way,
!    a user can traverse a cycle by starting at any entry I1 of ISIG
!    which is negative, moving to I2 = ABS(ISIG(I1)), then to
!    ISIG(I2), and so on, until returning to I1.
!
!    Output, integer ISGN, the "sign" of the permutation, which is
!    +1 if the permutation is even, -1 if odd.  Every permutation
!    may be produced by a certain number of pairwise switches.
!    If the number of switches is even, the permutation itself is
!    called even.
!
!    Output, integer NCYCLE, the number of cycles in the
!    permutation.
!
!    Input, integer IOPT, requests tagging.
!    0, the permutation will not be tagged.
!    1, the permutation will be tagged.
!
  implicit none

  integer n

  integer i
  integer i1
  integer i2
  integer iopt
  integer is
  integer isgn
  integer isig(n)
  integer ncycle

  is = 1
  ncycle = n

  do i = 1, n

    i1 = isig(i)

    do while ( i < i1 )
      ncycle = ncycle - 1
      i2 = isig(i1)
      isig(i1) = - i2
      i1 = i2
    end do

    if ( iopt /= 0 ) then
      is = - isign ( 1, isig(i) )
    end if

    isig(i) = isign ( isig(i), is )

  end do

  isgn = 1 - 2 * mod ( n-ncycle, 2 )

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
subroutine perm_inv ( n, isig )

!*****************************************************************************80
!
!! perm_inv() inverts a permutation.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 July 2000
!
!  Parameters:
!
!    Input, integer N, the number of objects being permuted.
!
!    Input/output, integer ISIG(N).
!
!    On input, ISIG describes a permutation.
!
!    ISIG is used to represent a permutation by the convention that
!    the permutation maps the letter I to ISIG(I).  Thus, if ISIG
!    contains the values (4, 1, 3, 2), then the permutation
!    represented permutes 1 to 4, 2 to 1, 3 to 3, and 4 to 2.
!
!    On output, ISIG describes the inverse permutation
!
  implicit none

  integer n

  integer i
  integer i0
  integer i1
  integer i2
  integer is
  integer isig(n)

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_INV(): Fatal error!'
    write ( *, '(a,i8)' ) '  Input value of N = ', n
    stop 1
  end if

  is = 1

  do i = 1, n

    i1 = isig(i)

    do while ( i < i1 )
      i2 = isig(i1)
      isig(i1) = - i2
      i1 = i2
    end do

    is = - isign ( 1, isig(i) )
    isig(i) = isign ( isig(i), is )

  end do

  do i = 1, n

    i1 = - isig(i)

    if ( 0 <= i1 ) then

      i0 = i

      do

        i2 = isig(i1)
        isig(i1) = i0

        if ( i2 < 0 ) then
          exit
        end if

        i0 = i1
        i1 = i2

      end do

    end if

  end do

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
subroutine r8vec2_print ( n, a1, a2, title )

!*****************************************************************************80
!
!! r8vec2_print() prints a pair of R8VEC's.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 June 2000
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
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a1(n)
  real ( kind = rk ) a2(n)
  integer i
  character ( len = * ) title

  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i8,2g14.6)' ) i, a1(i), a2(i)
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
!    This version by John Burkardt.
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

