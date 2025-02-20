subroutine catalan ( n, c )

!*****************************************************************************80
!
!! catalan() computes the Catalan numbers, from C(0) to C(N).
!
!  First values:
!
!     C(0)     1
!     C(1)     1
!     C(2)     2
!     C(3)     5
!     C(4)    14
!     C(5)    42
!     C(6)   132
!     C(7)   429
!     C(8)  1430
!     C(9)  4862
!    C(10) 16796
!
!  Formula:
!
!    C(N) = (2*N)! / ( (N+1) * (N!) * (N!) )
!         = 1 / (N+1) * COMB ( 2N, N )
!         = 1 / (2N+1) * COMB ( 2N+1, N+1).
!
!  Recursion:
!
!    C(N) = 2 * (2*N-1) * C(N-1) / (N+1)
!    C(N) = SUM ( I = 1 to N-1 ) C(I) * C(N-I)
!
!  Comments:
!
!    The Catalan number C(N) counts:
!
!    1) the number of binary trees on N vertices;
!    2) the number of ordered trees on N+1 vertices;
!    3) the number of full binary trees on 2N+1 vertices;
!    4) the number of well formed sequences of 2N parentheses;
!    5) number of ways 2N ballots can be counted, in order,
!       with N positive and N negative, so that the running sum
!       is never negative;
!    6) the number of standard tableaus in a 2 by N rectangular Ferrers diagram;
!    7) the number of monotone functions from [1..N} to [1..N} which
!       satisfy f(i) <= i for all i,
!    8) the number of ways to triangulate a polygon with N+2 vertices.
!
!  Example:
!
!    N = 3
!
!    ()()()
!    ()(())
!    (()())
!    (())()
!    ((()))
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dennis Stanton, Dennis White,
!    Constructive Combinatorics,
!    Springer Verlag, New York, 1986.
!
!  Parameters:
!
!    Input, integer N, the number of Catalan numbers desired.
!
!    Output, integer C(0:N), the Catalan numbers from C(0) to C(N).
!
  implicit none

  integer n

  integer i
  integer c(0:n)

  c(0) = 1
!
!  The extra parentheses ensure that the integer division is
!  done AFTER the integer multiplication.
!
  do i = 1, n
    c(i) = ( c(i-1) * 2 * ( 2 * i - 1 ) ) / ( i + 1 )
  end do

  return
end
subroutine catalan_values ( n_data, n, c )

!*****************************************************************************80
!
!! CATALAN_VALUES returns some values of the Catalan numbers.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Binomial[2*n,n] / ( n + 1 )
!
!  First values:
!
!     C(0)     1
!     C(1)     1
!     C(2)     2
!     C(3)     5
!     C(4)    14
!     C(5)    42
!     C(6)   132
!     C(7)   429
!     C(8)  1430
!     C(9)  4862
!    C(10) 16796
!
!  Formula:
!
!    C(N) = (2*N)! / ( (N+1) * (N!) * (N!) )
!         = 1 / (N+1) * COMB ( 2N, N )
!         = 1 / (2N+1) * COMB ( 2N+1, N+1).
!
!  Recursion:
!
!    C(N) = 2 * (2*N-1) * C(N-1) / (N+1)
!    C(N) = sum ( 1 <= I <= N-1 ) C(I) * C(N-I)
!
!  Discussion:
!
!    The Catalan number C(N) counts:
!
!    1) the number of binary trees on N vertices;
!    2) the number of ordered trees on N+1 vertices;
!    3) the number of full binary trees on 2N+1 vertices;
!    4) the number of well formed sequences of 2N parentheses;
!    5) the number of ways 2N ballots can be counted, in order,
!       with N positive and N negative, so that the running sum
!       is never negative;
!    6) the number of standard tableaus in a 2 by N rectangular Ferrers diagram;
!    7) the number of monotone functions from [1..N} to [1..N} which
!       satisfy f(i) <= i for all i;
!    8) the number of ways to triangulate a polygon with N+2 vertices.
!
!  Example:
!
!    N = 3
!
!    ()()()
!    ()(())
!    (()())
!    (())()
!    ((()))
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer N, the order of the Catalan number.
!
!    Output, integer C, the value of the Catalan number.
!
  implicit none

  integer, parameter :: n_max = 11

  integer c
  integer, save, dimension ( n_max ) :: c_vec = (/ &
    1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862, 16796 /)
  integer n
  integer n_data
  integer, save, dimension ( n_max ) :: n_vec = (/ &
     0,  1,  2,  3,  4, 5,  6,  7,  8,  9,  10 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    c = 0
  else
    n = n_vec(n_data)
    c = c_vec(n_data)
  end if

  return
end
subroutine cbt_traverse ( depth )

!*****************************************************************************80
!
!! CBT_TRAVERSE traverses a complete binary tree of given depth.
!
!  Discussion:
!
!    There will be 2^DEPTH terminal nodes of the complete binary tree.
!
!    This function traverses the tree, and prints out a binary code of 0's
!    and 1's each time it encounters a terminal node.  This results in a
!    printout of the binary digits from 0 to 2^DEPTH - 1.
!
!    The function is intended as a framework to be used to traverse a binary
!    tree.  Thus, in practice, a user would insert some action when a terminal
!    node is encountered.
!
!    Another use would occur when a combinatorial search is being made, for
!    example in a knapsack problem.  Each binary string then represents which
!    objects are to be included in the knapsack.  In that case, the traversal
!    could be speeded up by noticing cases where a nonterminal node has been
!    reached, but the knapsack is already full, in which case the only solution
!    uses none of the succeeding items, or overfull, in which case no solutions
!    exist that include this initial path segment.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 December 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DEPTH, the depth of the tree.
!
  implicit none

  integer depth

  integer b(depth)
  integer direction
  integer :: DOWNLEFT = 1
  integer k
  integer p
  integer :: UP = 3
  integer :: UPDOWNRIGHT = 2

  if ( depth < 1 ) then
    return
  end if

  b(1:depth) = 0
  p = 0
  direction = DOWNLEFT
  k = 0

  do
!
!  Try going in direction DOWNLEFT.
!
    if ( direction == DOWNLEFT ) then
      p = p + 1
      b(p) = 0
      if ( p < depth ) then
        write ( *, '(2x,a1,2x,4x,2x,10i1)' ) ' ', b(1:p)
      else
        write ( *, '(2x,a1,2x,i4,2x,10i1)' ) '(', k, b(1:depth)
        k = k + 1
        direction = UPDOWNRIGHT
      end if
    end if
!
!  Try going in direction UPDOWNRIGHT.
!
    if ( direction == UPDOWNRIGHT ) then
      b(p) = + 1
      if ( p < depth ) then
        write ( *, '(2x,a1,2x,4x,2x,10i1)' ) ' ', b(1:p)
        direction = DOWNLEFT
      else
        write ( *, '(2x,a1,2x,i4,2x,10i1)' ) ')', k, b(1:depth)
        k = k + 1
        direction = UP
      end if
    end if
!
!  Try going in direction UP.
!
    if ( direction == UP ) then
      p = p - 1
      if ( 1 <= p ) then
        write ( *, '(2x,a1,2x,4x,2x,10i1)' ) ' ', b(1:p)
        if ( b(p) == 0 ) then
          direction = UPDOWNRIGHT
        end if
      else
        exit
      end if
    end if

  end do

  return
end
subroutine graph_adj_edge_count ( adj, nnode, nedge )

!*****************************************************************************80
!
!! GRAPH_ADJ_EDGE_COUNT counts the number of edges in a graph.
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
!    ADJ(I,J) is 1 if there is an edge from node I to node J.
!
!    Input, integer NNODE, the number of nodes.
!
!    Output, integer NEDGE, the number of edges in the graph.
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
subroutine graph_adj_is_node_connected ( adj, nnode, result )

!*****************************************************************************80
!
!! GRAPH_ADJ_IS_NODE_CONNECTED determines if a graph is nodewise connected.
!
!  Discussion:
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
!    24 June 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ADJ(NNODE,NNODE), the adjacency matrix for the 
!    graph.  ADJ(I,J) is nonzero if there is an edge from node I to node J.
!
!    Input, integer NNODE, the number of nodes.
!
!    Output, integer RESULT.
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
!  If no neighbors were found, exit.
!
    if ( jhi < jlo ) then
      exit
    end if
!
!  If neighbors were found, then go back and find THEIR neighbors.
!
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
!! GRAPH_ADJ_IS_TREE determines whether a graph is a tree.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    24 June 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ADJ(NNODE,NNODE), the adjacency matrix for the 
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
subroutine graph_arc_degree ( nnode, nedge, inode, jnode, degree )

!*****************************************************************************80
!
!! GRAPH_ARC_DEGREE determines the degree of the nodes of a graph.
!
!  Discussion:
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
!  Parameters:
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, integer NEDGE, the number of edges.
!
!    Input, integer INODE(NEDGE), JNODE(NEDGE), the pairs of nodes
!    that form the edges.
!
!    Output, integer DEGREE(NNODE), the degree of each node, 
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
      write ( *, '(a)' ) 'GRAPH_ARC_DEGREE - Fatal error!'
      write ( *, '(a,i8)' ) '  Out-of-range node value = ', n
      stop
    end if

    n = jnode(i)
    if ( 1 <= n .and. n <= nnode ) then
      degree(n) = degree(n) + 1
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GRAPH_ARC_DEGREE - Fatal error!'
      write ( *, '(a,i8)' ) '  Out-of-range node value = ', n
      stop
    end if

  end do

  return
end
subroutine graph_arc_is_tree ( nedge, inode, jnode, nnode, result )

!*****************************************************************************80
!
!! GRAPH_ARC_IS_TREE determines whether a graph is a tree.
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
  integer nnode2
  integer result

  call graph_arc_to_graph_adj ( nedge, inode, jnode, adj, nnode2 )

  call graph_adj_is_tree ( adj, nnode, result )

  return
end
subroutine graph_arc_node_count ( nedge, inode, jnode, nnode )

!*****************************************************************************80
!
!! GRAPH_ARC_NODE_COUNT counts the number of nodes in a graph.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 July 2013
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
!    Output, integer NNODE, the number of distinct nodes.
!
  implicit none

  integer nedge

  integer iedge
  integer inode(nedge)
  integer jnode(nedge)
  integer knode(2*nedge)
  integer nnode
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

  call i4vec_sorted_unique_count ( 2*nedge, knode, nnode )

  return
end
subroutine graph_arc_node_max ( nedge, inode, jnode, node_max )

!*****************************************************************************80
!
!! GRAPH_ARC_NODE_MAX determines the maximum node label.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 July 2013
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
!    Output, integer NODE_MAX, the maximum node index.
!
  implicit none

  integer nedge

  integer inode(nedge)
  integer jnode(nedge)
  integer node_max

  node_max = max ( maxval ( inode(1:nedge) ), maxval ( jnode(1:nedge) ) )

  return
end
subroutine graph_arc_print ( nedge, inode, jnode, title )

!*****************************************************************************80
!
!! GRAPH_ARC_PRINT prints out a graph from an edge list.
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
!! GRAPH_ARC_TO_GRAPH_ADJ converts an arc list graph to an adjacency graph.
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
!    Output, integer NNODE, the number of nodes.
!
  implicit none

  integer nedge
  integer nnode

  integer adj(nedge+1,nedge+1)
  integer i
  integer inode(nedge)
  integer j
  integer jnode(nedge)
  integer k
!
!  Determine the number of nodes.
!
  call graph_arc_node_count ( nedge, inode, jnode, nnode )

  adj(1:nnode,1:nnode) = 0

  do k = 1, nedge
    i = inode(k)
    j = jnode(k)
    adj(i,j) = 1
    adj(j,i) = 1
  end do

  return
end
function i4_uniform_ab ( a, b, seed )

!*****************************************************************************80
!
!! I4_UNIFORM_AB returns a scaled pseudorandom I4 between A and B.
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
!  Parameters:
!
!    Input, integer A, B, the limits of the interval.
!
!    Input/output, integer SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer I4_UNIFORM_AB, a number between A and B.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer a
  integer b
  integer, parameter :: i4_huge = 2147483647
  integer i4_uniform_ab
  integer k
  real ( kind = rk ) r
  integer seed
  integer value

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_UNIFORM_AB - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r = real ( seed, kind = rk ) * 4.656612875D-10
!
!  Scale R to lie between A-0.5 and B+0.5.
!
  r = ( 1.0D+00 - r ) * ( real ( min ( a, b ), kind = rk ) - 0.5D+00 ) & 
    +             r   * ( real ( max ( a, b ), kind = rk ) + 0.5D+00 )
!
!  Use rounding to convert R to an integer between A and B.
!
  value = nint ( r )

  value = max ( value, min ( a, b ) )
  value = min ( value, max ( a, b ) )

  i4_uniform_ab = value

  return
end
subroutine i4mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! I4MAT_PRINT prints an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 June 2003
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
!    Input, integer A(M,N), the matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer m
  integer n

  integer a(m,n)
  integer ihi
  integer ilo
  integer jhi
  integer jlo
  character ( len = * ) title

  ilo = 1
  ihi = m
  jlo = 1
  jhi = n

  call i4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

  return
end
subroutine i4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! I4MAT_PRINT_SOME prints some of an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4's.
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
!    Input, integer A(M,N), an M by N matrix to be printed.
!
!    Input, integer ILO, JLO, the first row and column to print.
!
!    Input, integer IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: incx = 10
  integer m
  integer n

  integer a(m,n)
  character ( len = 8 )  ctemp(incx)
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
      write ( ctemp(j2), '(i8)' ) j
    end do

    write ( *, '(''  Col '',10a8)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        write ( ctemp(j2), '(i8)' ) a(i,j)

      end do

      write ( *, '(i5,a,10a8)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine i4vec_heap_d ( n, a )

!*****************************************************************************80
!
!! I4VEC_HEAP_D reorders an I4VEC into an descending heap.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    A descending heap is an array A with the property that, for every index J,
!    A(J) >= A(2*J) and A(J) >= A(2*J+1), (as long as the indices
!    2*J and 2*J+1 are legal).
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
!    Combinatorial Algorithms for Computers and Calculators,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
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
  do i = n / 2, 1, -1
!
!  Copy the value out of the parent node.
!  Position IFREE is now "open".
!
    key = a(i)
    ifree = i

    do
!
!  Positions 2*IFREE and 2*IFREE + 1 are the descendants of position
!  IFREE.  (One or both may not exist because they exceed N.)
!
      m = 2 * ifree
!
!  Does the first position exist?
!
      if ( n < m ) then
        exit
      end if
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
      if ( a(m) <= key ) then
        exit
      end if

      a(ifree) = a(m)
      ifree = m

    end do
!
!  Once there is no more shifting to do, KEY moves into the free spot IFREE.
!
    a(ifree) = key

  end do

  return
end
subroutine i4vec_indicator ( n, a )

!*****************************************************************************80
!
!! I4VEC_INDICATOR sets an I4VEC to the indicator vector.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 May 2007
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
!! I4VEC_PRINT prints an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 May 2010
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
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer n

  integer a(n)
  integer i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,a,2x,i12)' ) i, ':', a(i)
  end do

  return
end
subroutine i4vec_sort_heap_a ( n, a )

!*****************************************************************************80
!
!! I4VEC_SORT_HEAP_A ascending sorts an I4VEC using heap sort.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
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
  integer t

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
  t    = a(1)
  a(1) = a(n)
  a(n) = t
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
    t    = a(1)
    a(1) = a(n1)
    a(n1) = t

  end do

  return
end
subroutine i4vec_sorted_unique_count ( n, a, unique_num )

!*****************************************************************************80
!
!! I4VEC_SORTED_UNIQUE_COUNT counts the unique elements in a sorted I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    Because the array is sorted, this algorithm is O(N).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements of A.
!
!    Input, integer A(N), the sorted array to examine.
!
!    Output, integer UNIQUE_NUM, the number of unique elements.
!
  implicit none

  integer n

  integer a(n)
  integer i
  integer unique_num

  if ( n < 1 ) then
    unique_num = 0
    return
  end if

  unique_num = 1

  do i = 2, n

    if ( a(i-1) /= a(i) ) then
      unique_num = unique_num + 1
    end if

  end do

  return
end
subroutine pruefer_to_tree_arc ( nnode, iarray, inode, jnode )

!*****************************************************************************80
!
!! PRUEFER_TO_TREE_ARC is given a Pruefer code, and computes the tree.
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
 
  do i = 1, nnode - 2
    iwork(iarray(i)) = iwork(iarray(i)) + 1
  end do
!
!  Now process each entry in the Pruefer code.
!
  inode(1:nnode-1) = -1
  jnode(1:nnode-1) = -1

  do i = 1, nnode - 2
 
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
subroutine pruefer_to_tree_2 ( nnode, iarray, itree )

!*****************************************************************************80
!
!! PRUEFER_TO_TREE_2 produces the edge list of a tree from its Pruefer code.
!
!  Discussion:
!
!    One can thus exhibit all trees on N nodes, produce
!    one at random, find the M-th one on the list, etc, by
!    manipulating the Pruefer codes.
!
!    For every labeled tree on N nodes, there is a unique N-2 tuple
!    of integers A1 through AN-2, with each A between 1 and N.  There
!    are N^(N-2) such sequences, and each one is associated with exactly
!    one tree.
!
!    This routine apparently assumes that the Pruefer code is
!    generated by taking the LOWEST labeled terminal node each time.
!    This is not consistent with PRUEFER_TO_TREE and TREE_TO_PRUEFER.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 June 2013
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis. Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer NNODE, number of nodes in desired tree.
!
!    Input, integer IARRAY(NNODE).  IARRAY(I), I = 1, NNODE-2 
!    is the Pruefer code for the tree.
!
!    Output, integer ITREE(NNODE); the I-th edge of the tree
!    joins nodes I and ITREE(I).
!
  implicit none

  integer nnode

  integer i
  integer ir
  integer itree(nnode)
  integer iarray(nnode)
  integer j
  integer k
  integer kp
  integer l

  itree(1:nnode) = 0
 
  do i = nnode - 2, 1, -1
 
    l = iarray(i)
 
    if ( itree(l) == 0 ) then
      iarray(i) = - l
      itree(l) = - 1
    end if
 
  end do
 
  iarray(nnode-1) = nnode
!
!  Find next index K so that ITREE(K) is 0.
!
  k = 1
  j = 0
 
  do while ( itree(k) /= 0 )
    k = k + 1
  end do
 
  kp = k
 
  do
 
    j = j + 1
    ir = abs ( iarray(j) )
    itree(kp) = ir
 
    if ( j == nnode - 1 ) then
      exit
    end if
 
    if ( 0 < iarray(j) ) then
      do while ( itree(k) /= 0 )
        k = k + 1
      end do
      kp = k
      cycle
    end if
 
    if ( k < ir ) then
      itree(ir) = 0
      do while ( itree(k) /= 0 )
        k = k + 1
      end do
      kp = k
      cycle
    end if
 
    kp = ir

  end do
!
!  Restore the signs of IARRAY.
!
  iarray(1:nnode-2) = abs ( iarray(1:nnode-2) )
 
  return
end
function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = rk ) value.
!
!    For now, the input quantity SEED is an integer variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = rk ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: i4_huge = 2147483647
  integer k
  real ( kind = rk ) r8_uniform_01
  integer seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r8_uniform_01 = real ( seed, kind = rk ) * 4.656612875D-10

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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine tree_arc_center ( nnode, inode, jnode, center, eccent, parity )

!*****************************************************************************80
!
!! TREE_ARC_CENTER computes the center, eccentricity, and parity of a tree.
!
!  Discussion:
!
!    A tree is an undirected graph of N nodes, which uses N-1 edges,
!    and is connected.  
!
!    A graph with N-1 edges is not guaranteed to be a tree, and so this
!    routine must first check that condition before proceeding.
!
!    The edge distance between two nodes I and J is the minimum number of
!    edges that must be traversed in a path from I and J.
!
!    The eccentricity of a node I is the maximum edge distance between
!    node I and the other nodes J in the graph.
!
!    The radius of a graph is the minimum eccentricity over all nodes
!    in the graph.
!
!    The diameter of a graph is the maximum eccentricity over all nodes
!    in the graph.
!
!    The center of a graph is the set of nodes whose eccentricity is 
!    equal to the radius, that is, the set of nodes of minimum eccentricity.
!
!    For a tree, the center is either a single node, or a pair of
!    neighbor nodes.
!
!    The parity of the tree is 1 if the center is a single node, or 2 if
!    the center is 2 nodes.
!
!    The center of a tree can be found by removing all "leaves", that is,
!    nodes of degree 1.  This step is repeated until only 1 or 2 nodes
!    are left.
!
!    Thanks to Alexander Sax for pointing out that a previous version of the
!    code was failing when the tree had an odd parity, that is, a single
!    center node, 15 April 2013.
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
!    Input, integer NNODE, the number of nodes.
!
!    Input, integer INODE(NNODE-1), JNODE(NNODE-1), the edges of
!    the tree.  Edge I connects nodes INODE(I) and JNODE(I).
!
!    Output, integer CENTER(2).  CENTER(1) is the index of the
!    first node in the center.  CENTER(2) is 0 if there is only one node
!    in the center, or else the index of the second node.
!
!    Output, integer ECCENT, the eccentricity of the nodes in 
!    the center, and the radius of the the tree.
!
!    Output, integer PARITY, the parity of the tree, which is
!    normally 1 or 2.
!
  implicit none

  integer nnode

  integer center(2)
  integer degree(nnode)
  integer eccent
  integer i
  integer iedge
  integer ileaf
  integer inode(nnode-1)
  integer j
  integer jnode(nnode-1)
  integer list(nnode)
  integer nedge
  integer nleaf
  integer nnode2
  integer parity
  integer result

  eccent = 0
  center(1) = 0
  center(2) = 0
  parity = 0

  if ( nnode <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TREE_ARC_CENTER - Fatal error!'
    write ( *, '(a)' ) '  NNODE <= 0.'
    stop
  else if ( nnode == 1 ) then
    eccent = 0
    center(1) = 1
    center(2) = 0
    parity = 1
    return
  else if ( nnode == 2 ) then
    eccent = 1
    center(1) = 1
    center(2) = 2
    parity = 2
    return
  end if
!
!  Is this graph really a tree?
!
  nedge = nnode - 1
  call graph_arc_is_tree ( nedge, inode, jnode, nnode, result )

  if ( result == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TREE_ARC_CENTER - Fatal error!'
    write ( *, '(a)' ) '  This graph is NOT a tree.'
    stop
  end if
!
!  Compute the degrees.
!
  call graph_arc_degree ( nnode, nedge, inode, jnode, degree )
!
!  Defoliate the tree.
!
  nnode2 = nnode

  do

    eccent = eccent + 1
!
!  Find and mark the leaves.
!
    nleaf = 0

    do i = 1, nnode

      if ( degree(i) == 1 ) then
        nleaf = nleaf + 1
        list(nleaf) = i
      end if

    end do
!
!  Delete the leaves.
!
    do ileaf = 1, nleaf

      i = list(ileaf)

      iedge = 0
      j = 0

      do

        iedge = iedge + 1

        if ( nedge < iedge ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'TREE_ARC_CENTER - Fatal error!'
          write ( *, '(a)' ) '  Data or algorithm failure.'
          stop
        end if

        if ( inode(iedge) == i ) then
          j = jnode(iedge)
          inode(iedge) = - inode(iedge)
          jnode(iedge) = - jnode(iedge)
        else if ( jnode(iedge) == i ) then
          j = inode(iedge)
          inode(iedge) = - inode(iedge)
          jnode(iedge) = - jnode(iedge)
        end if

        if ( j /= 0 ) then
          exit
        end if

      end do

      degree(i) = -1
      nnode2 = nnode2 - 1
      degree(j) = degree(j) - 1
!
!  If the other node has degree 0, we must have just finished
!  stripping all leaves from the tree, leaving a single node.
!  Don't kill it here.  It is our odd center.
!
!     if ( degree(j) == 0 ) then
!       nnode2 = nnode2 - 1
!     end if

    end do
!
!  Find the remaining nodes.
!
    nnode2 = 0

    do i = 1, nnode

      if ( 0 <= degree(i) ) then
        nnode2 = nnode2 + 1
        list(nnode2) = i
      end if

    end do
!
!  If at least 3, more pruning is required.
!
    if ( nnode2 < 3 ) then
      exit
    end if

  end do
!
!  If only one or two nodes left, we are done.
!
  parity = nnode2

  center(1:nnode2) = list(1:nnode2)
  inode(1:nedge) = abs ( inode(1:nedge) )
  jnode(1:nedge) = abs ( jnode(1:nedge) )

  return
end
subroutine tree_arc_diam ( nnode, inode, jnode, diam, label, n1, n2 )

!*****************************************************************************80
!
!! TREE_ARC_DIAM computes the "diameter" of a tree.
!
!  Discussion:
!
!    A tree is an undirected graph of N nodes, which uses N-1 edges,
!    and is connected.  
!
!    A graph with N-1 edges is not guaranteed to be a tree, and so this
!    routine must first check that condition before proceeding.
!
!    The diameter of a graph is the length of the longest possible
!    path that never repeats an edge.
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
!    Input, integer NNODE, the number of nodes.
!
!    Input, integer INODE(NNODE-1), JNODE(NNODE-1), the edges 
!    of the tree.  Edge I connects nodes INODE(I) and JNODE(I).
!
!    Output, integer DIAM, the length of the longest path 
!    in the tree.
!
!    Output, integer LABEL(NNODE), marks the path between 
!    nodes N1 and N2.  Node I is in this path if LABEL(I) is 1.
!
!    Output, integer N1, N2, the indices of two nodes in the 
!    tree which are separated by DIAM edges.
!
  implicit none

  integer nnode

  integer degree(nnode)
  integer diam
  integer i
  integer inode(nnode-1)
  integer invals
  integer j
  integer jnode(nnode-1)
  integer k
  integer kstep
  integer label(nnode)
  integer n1
  integer n2
  integer nabe
  integer nedge
  integer result

  if ( nnode <= 0 ) then
    diam = 0
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TREE_ARC_DIAM - Fatal error!'
    write ( *, '(a)' ) '  NNODE <= 0.'
    stop
  end if

  if ( nnode == 1 ) then
    diam = 0
    n1 = 1
    n2 = 1
    return
  end if
!
!  Is this graph really a tree?
!
  nedge = nnode - 1
  call graph_arc_is_tree ( nedge, inode, jnode, nnode, result )

  if ( result == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TREE_ARC_DIAM - Fatal error!'
    write ( *, '(a)' ) '  This graph is NOT a tree.'
    stop
  end if

  call i4vec_indicator ( nnode, label )
!
!  On step K:
!
!    Identify the terminal and interior nodes.
!
!    If there are no interior nodes left, 
!
!      then there are just two nodes left at all.  The diameter is 2*K-1, 
!      and a maximal path extends between the nodes whose labels are 
!      contained in the two remaining terminal nodes.
!
!    Else
!
!      The label of each terminal node is passed to its interior neighbor.
!      If more than one label arrives, take any one.
!
!      The terminal nodes are removed.
!
  kstep = 0

  do

    kstep = kstep + 1
!
!  Compute the degree of each node.
!
    degree(1:nnode) = 0        

    do j = 1, nedge

      k = inode(j)
      if ( 0 < k ) then
        degree(k) = degree(k) + 1
      end if

      k = jnode(j)
      if ( 0 < k ) then
        degree(k) = degree(k) + 1
      end if

    end do
!
!  Count the number of interior nodes.
!
    invals = 0
    do i = 1, nnode
      if ( 1 < degree(i) ) then
        invals = invals + 1
      end if
    end do
!
!  If there are 1 or 0 interior nodes, it's time to stop.
!
    if ( invals == 1 ) then

      diam = 2 * kstep
      exit
    
    else if ( invals == 0 ) then

      diam = 2 * kstep - 1
      exit

    end if
!
!  If there are at least two interior nodes, then chop off the 
!  terminal nodes and pass their labels inward.
!
    do k = 1, nnode

      if ( degree(k) == 1 ) then

        do j = 1, nedge

          if ( inode(j) == k ) then
            nabe = jnode(j)
            label(nabe) = label(k)
            inode(j) = - inode(j)
            jnode(j) = - jnode(j)
          else if ( jnode(j) == k ) then
            nabe = inode(j)
            label(nabe) = label(k)
            inode(j) = - inode(j)
            jnode(j) = - jnode(j)
          end if

        end do

      end if

    end do

  end do
!
!  Now get the labels from two of the remaining terminal nodes.
!  The nodes represented by these labels will be a diameter apart.
!
  n1 = 0
  n2 = 0

  do i = 1, nnode
    if ( degree(i) == 1 ) then
      if ( n1 == 0 ) then
        n1 = label(i)
      else if ( n2 == 0 ) then
        n2 = label(i)
      end if
    end if
  end do
!
!  Set the labels of the interior node (if any) and nodes marked
!  N1 and N2 to 1, and all others to 0.  This will label the nodes on the path.
!
  if ( invals == 1 ) then

    do i = 1, nnode
      if ( 1 < degree(i) ) then
        label(i) = 1
      end if
    end do

  end if

  do i = 1, nnode
    if ( label(i) == n1 .or. label(i) == n2 ) then
      label(i) = 1
    else
      label(i) = 0
    end if
  end do
!
!  Clean up the arrays.
!
  do j = 1, nedge
    inode(j) = abs ( inode(j) )
    k = inode(j)
    jnode(j) = abs ( jnode(j) )
    k = jnode(j)
  end do

  return
end
subroutine tree_arc_random ( nnode, seed, code, inode, jnode )

!*****************************************************************************80
!
!! TREE_ARC_RANDOM selects a random labeled tree and its Pruefer code.
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
!    Input/output, integer SEED, a seed for the random number 
!    generator.
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
  integer seed

  if ( nnode <= 0  ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TREE_ARC_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  NNODE = ', nnode
    write ( *, '(a)' ) '  but NNODE must be at least 1.'
    stop
  end if

  if ( nnode <= 2 ) then
    return
  end if

  call vec_random ( nnode-2, nnode, seed, code )
 
  code(1:nnode-2) = code(1:nnode-2) + 1
 
  call pruefer_to_tree_arc ( nnode, code, inode, jnode )
 
  return
end
subroutine tree_arc_to_pruefer ( nnode, inode, jnode, code )

!*****************************************************************************80
!
!! TREE_ARC_TO_PRUEFER is given a labeled tree, and computes its Pruefer code.
!
!  Discussion:
!
!    A tree is an undirected graph of N nodes, which uses N-1 edges,
!    and is connected.  
!
!    A graph with N-1 edges is not guaranteed to be a tree, and so this
!    routine must first check that condition before proceeding.
!
!    The Pruefer code is a correspondence between all labeled trees of
!    N nodes, and all list of N-2 integers between 1 and N (with repetition
!    allowed).  The number of labeled trees on N nodes is therefore N^(N-2).
!
!    The Pruefer code is constructed from the tree as follows:
!
!    A terminal node on the tree is defined as a node with only one neighbor.
!
!    Consider the set of all terminal nodes on the tree.  Take the one
!    with the highest label, I.  Record the label of its neighbor, J.
!    Delete node I and the edge between node I and J.
!
!    J is the first entry in the Pruefer code for the tree.   Repeat
!    the operation a total of N-2 times to get the complete code.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 April 2013
!
!  Reference:
!
!    Dennis Stanton, Dennis White,
!    Constructive Combinatorics,
!    Springer Verlage, New York, 1986.
!
!  Parameters:
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, integer INODE(NNODE-1), JNODE(NNODE-1), the edge array 
!    of the tree.  The I-th edge joins nodes INODE(I) and JNODE(I).
!
!    Output, integer CODE(NNODE-2), the Pruefer code of the tree.
!
  implicit none

  integer nnode

  integer code(nnode-2)
  integer degree(nnode)
  integer i
  integer i2
  integer iterm
  integer inode(nnode-1)
  integer j
  integer jnode(nnode-1)
  integer jsave
  integer nedge
  integer result
!
!  Is this graph really a tree?
!
  nedge = nnode - 1
  call graph_arc_is_tree ( nedge, inode, jnode, nnode, result )

  if ( result == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TREE_ARC_TO_PRUEFER - Fatal error!'
    write ( *, '(a)' ) '  This graph is NOT a tree.'
    stop
  end if
!
!  Compute the degree of each node.
!
  nedge = nnode - 1
  call graph_arc_degree ( nnode, nedge, inode, jnode, degree )
!
!  Compute the next term of the Pruefer code.
!
  do i = 1, nnode - 2
!
!  Find the terminal node with the highest label.
!
    iterm = 0
 
    do j = 1, nnode
      if ( degree(j) == 1 ) then
        iterm = j
      end if
    end do
!
!  Find the edge that includes this node, and note the
!  index of the other node.
!
    do j = 1, nnode - 1

      jsave = j
 
      if ( inode(j) == iterm ) then
        i2 = 2
        exit
      else if ( jnode(j) == iterm ) then
        i2 = 1
        exit
      end if
 
    end do
!
!  Delete the edge from the tree.
!
    degree(inode(jsave)) = degree(inode(jsave)) - 1
    degree(jnode(jsave)) = degree(jnode(jsave)) - 1
!
!  Add the neighbor of the node to the Pruefer code.
!
    if ( i2 == 1 ) then
      code(i) = inode(jsave)
    else
      code(i) = jnode(jsave)
    end if
!
!  Negate the nodes in the edge list to mark them as used.
!
    inode(jsave) = - inode(jsave)
    jnode(jsave) = - jnode(jsave)
 
  end do
!
!  Before returning, restore the original form of the edge list.
!
  inode(1:nnode-1) = abs ( inode(1:nnode-1) )
  jnode(1:nnode-1) = abs ( jnode(1:nnode-1) )
 
  return
end
subroutine tree_enum ( nnode, ntree )

!*****************************************************************************80
!
!! TREE_ENUM enumerates the labeled trees on NNODE nodes.
!
!  Discussion:
!
!    The formula is due to Cauchy.
!
!  Example:
!
!    NNODE      NTREE
!
!    0              1
!    1              1
!    2              1
!    3              3
!    4             16
!    5            125
!    6           1296
!    7          16807
!    8         262144
!    9        4782969
!   10      100000000
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NNODE, the number of nodes in each tree.
!    NNODE must normally be at least 3, but for this routine,
!    any value of NNODE is allowed.  Values of NNODE greater than 10
!    will probably overflow.
!
!    Output, integer NTREE, the number of distinct labeled trees.
!
  implicit none

  integer nnode
  integer ntree

  if ( nnode < 0 ) then
    ntree = 0
  else if ( nnode == 0 ) then
    ntree = 1
  else if ( nnode == 1 ) then
    ntree = 1
  else if ( nnode == 2 ) then
    ntree = 1
  else
    ntree = nnode ** ( nnode - 2 )
  end if

  return
end
subroutine tree_parent_next ( nnode, iarray, code, itree, more )

!*****************************************************************************80
!
!! TREE_PARENT_NEXT generates, one at a time, all labeled trees.
!
!  Discussion:
!
!    The routine also returns the corresponding Pruefer codes.
!
!  Formula:
!
!    There are N^(N-2) labeled trees on N nodes (Cayley's formula).
!
!    The number of trees in which node I has degree D(I) is the
!    multinomial coefficient: ( N-2; D(1)-1, D(2)-1, ..., D(N)-1 ).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 August 2000
!
!  Parameters:
!
!    Input, integer NNODE, the number of nodes to be used in 
!    the trees.
!
!    Workspace, integer IARRAY(NNODE).
!
!    Output, integer CODE(NNODE).  The first NNODE-2 entries 
!    of CODE contain the Pruefer code for the given labeled tree.
!
!    Output, integer ITREE(NNODE).  The first NNODE-1 entries 
!    of ITREE describe the edges that go between the nodes.  Each pair
!    (I, ITREE(I)) represents an edge.  Thus if ITREE(5) = 3,
!    there is an edge from node 3 to node 5.
!
!    Input/output, logical MORE.  On the first call only, the
!    user is required to set MORE = .FALSE.  Then call the routine, and
!    it will return information about the first tree
!    as well as setting MORE to the value .TRUE.
!    Keep calling to get another tree until MORE is .FALSE.
!    on return, at which point there are no more trees.
!
  implicit none

  integer nnode

  integer code(nnode)
  integer, dimension ( nnode ) :: iarray
  integer itree(nnode)
  logical more

  call vec_next ( nnode-2, nnode, iarray, more )
 
  code(1:nnode-2) = iarray(1:nnode-2) + 1
 
  call pruefer_to_tree_2 ( nnode, code, itree )
 
  return
end
subroutine tree_parent_to_arc ( nnode, parent, nedge, inode, jnode )

!*****************************************************************************80
!
!! TREE_PARENT_TO_ARC converts a tree from parent to arc representation.
!
!  Discussion:
!
!    Parent representation lists the parent node of each node.  For a
!    tree of N nodes, one node has a parent of 0, representing a null link.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NNODE, the number of nodes in the tree.
!
!    Input, integer PARENT(NNODE), the parent node representation 
!    of the tree.
!
!    Output, integer NEDGE, the number of edges, normally NNODE-1.
!
!    Output, integer INODE(NEDGE), JNODE(NEDGE), pairs of nodes
!    that define the links.
!
  implicit none

  integer nnode

  integer i
  integer inode(nnode-1)
  integer jnode(nnode-1)
  integer nedge
  integer parent(nnode)

  nedge = 0

  do i = 1, nnode

    if ( parent(i) /= 0 ) then
      nedge = nedge + 1
      inode(nedge) = i
      jnode(nedge) = parent(i)
    end if

  end do

  return
end
subroutine tree_rb_enum ( n, num )

!*****************************************************************************80
!
!! TREE_RB_ENUM returns the number of rooted binary trees with N nodes.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of nodes in the rooted 
!    binary tree.  N should be odd.
!
!    Output, integer NUM, the number of rooted binary trees 
!    with N nodes.
!
  implicit none

  integer n

  integer c(0:(n-1)/2)
  integer num

  if ( n < 0 ) then

    num = 0

  else if ( n == 0 ) then

    num = 1

  else if ( mod ( n, 2 ) == 0 ) then

    num = 0

  else

    call catalan ( ( n - 1 ) / 2, c )

    num = c((n-1)/2)

  end if

  return
end
subroutine tree_rb_lex_next ( n, a, more )

!*****************************************************************************80
!
!! TREE_RB_LEX_NEXT generates rooted binary trees in lexicographic order.
!
!  Discussion:
!
!    The information definining the tree of N nodes is stored in a vector 
!    of 0's and 1's, in preorder traversal form.  Essentially, the
!    shape of the tree is traced out with a pencil that starts at the root,
!    and stops at the very last null leaf.  The first time that a (non-null) 
!    node is encountered, a 1 is added to the vector, and the left 
!    descendant of the node is visited next.  When the path returns from
!    the first descendant, the second descendant is visited.  When then path
!    returns again, the path passes back up from the node to its parent.
!    A null leaf is encountered only once, and causes a zero to be added to 
!    the vector, and the path goes back up to the parent node.  
!
!    The lexicographic order is used to order the vectors of 1's and 0's.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 August 2000
!
!  Reference:
!
!    Frank Ruskey,
!    Combinatorial Generation,
!    To appear.
!
!  Parameters:
!
!    Input, integer N, the number of nodes in the rooted binary
!    tree.  N should be odd.
!
!    Input/output, integer A(N), the preorder traversal form for
!    the previous/next rooted binary tree.
!
!    Output, logical MORE, is TRUE if the next rooted binary tree was
!    returned on this call, or FALSE if there are no more rooted binary
!    trees, and the output of the previous call was the last one.
!
  implicit none

  integer n

  integer a(n)
  integer k
  logical more
  integer p
  integer q

  if ( .not. more ) then
    a(1:n-2:2) = 1
    a(2:n-1:2) = 0
    a(n) = 0
    more = .true.
    return
  end if
!
!  Find the last 1 in A.
!
  k = n
  do while ( a(k) == 0 )
    k = k - 1
  end do
  q = n - k - 1
!
!  Find the last 0 preceding the last 1 in A.
!  If there is none, then we are done, because 11...1100..00 
!  is the final element.
!
  do 

    if ( k == 1 ) then
      more = .false.
      return
    end if

    if ( a(k) == 0 ) then
      exit
    end if

    k = k - 1

  end do
	
  p = n - k - q - 1
  a(k) = 1
  a(k+1:n-2*p+1) = 0
  a(n-2*p+2:n-2:2) = 1
  a(n-2*p+3:n-1:2) = 0
  a(n) = 0

  return
end
subroutine tree_rb_to_parent ( n, a, parent )

!*****************************************************************************80
!
!! TREE_RB_TO_PARENT converts rooted binary tree to parent node representation.
!
!  Discussion:
!
!    Parent node representation of a tree assigns to each node a "parent" node,
!    which represents the first link of the path between the node and the 
!    root node.  The root node itself is assigned a parent of 0.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of nodes in the tree.
!
!    Input, integer A(N), the preorder traversal form for the
!    rooted binary tree.
!
!    Output, integer PARENT(N), the parent node representation 
!    of the tree.
!
  implicit none

  integer n

  integer a(n)
  integer dad
  integer k
  integer node
  integer node_num
  integer parent(n)
  integer use(n)

  node = 0
  node_num = 0

  do k = 1, n

    dad = node
    node_num = node_num + 1
    node = node_num
    parent(node) = dad

    if ( a(k) == 1 ) then

      use(node) = 0

    else

      use(node) = 2

      do while ( use(node) == 2 )
        node = dad
        if ( node == 0 ) then
          exit
        end if
        use(node) = use(node) + 1
        dad = parent(node)
      end do

    end if

  end do

  return
end
subroutine tree_rb_yule ( n, seed, a )

!*****************************************************************************80
!
!! TREE_RB_YULE adds two nodes to a rooted binary tree using the Yule model.
!
!  Discussion:
!
!    The Yule model is a simulation of how an evolutionary family tree
!    develops.  We start with a root node.  The internal nodes of the tree 
!    are inactive and never change.  Each pendant or leaf node of the
!    tree represents a biological family that can spontaneously "fission",
!    developing two new distinct sub families.  In graphical terms, the node
!    becomes internal, with two new leaf nodes depending from it.
!
!    The tree is stored in inorder traversal form.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer N, the number of nodes in the input
!    tree.  On output, this number has been increased, usually by 2.
!
!    Input/output, integer SEED, a seed for the random number
!    generator.
!
!    Input/output, integer A(*), the preorder traversal form 
!    for the rooted binary tree.  The number of entries in A is N.
!
  implicit none

  integer a(*)
  integer i
  integer i4_uniform_ab
  integer ileaf
  integer j
  integer jleaf
  integer n
  integer nleaf
  integer seed

  if ( n <= 0 ) then
    n = 1
    a(1) = 0
    return
  end if
!
!  Count the expected number of leaves, which are the 0 values.
!
  nleaf = ( n + 1 ) / 2
!
!  Choose a random number between 1 and NLEAF.
!
  ileaf = i4_uniform_ab ( 1, nleaf, seed )
!
!  Locate leaf number ILEAF.
!
  j = 0
  jleaf = 0
  do i = 1, n
    if ( a(i) == 0 ) then
      jleaf = jleaf + 1
    end if
    if ( jleaf == ileaf ) then
      j = i
      exit
    end if
  end do
!
!  Replace '0' by '100'
!
  a(n+2:j+2:-1) = a(n:j:-1)
  a(j:j+1) = (/ 1, 0 /)

  n = n + 2

  return
end
subroutine tree_rooted_code ( nnode, parent, code )

!*****************************************************************************80
!
!! TREE_ROOTED_CODE returns the code of a rooted tree.
!
!  Discussion:
!
!    This code for a rooted tree depends on the node ordering, so it's actually
!    the code for a labeled rooted tree.  To eliminate the effects of node
!    labeling, one could choose as the code for a tree the maximum of all
!    the codes associated with the different possible labelings of the tree.
!    There are more effective ways of arriving at this code than simply
!    generating all possible codes and comparing them.  
!
!    For a tree with NNODES, the code is a list of 2*NNODE 0's and 1's,
!    describing a traversal of the tree starting at an imaginary node 0,
!    moving "down" to the root (a code entry of 1), and then moving
!    "down" (1) or "up" (0) as the tree is traversed in a depth first
!    manner.  The final move must be from the root up to the imaginary
!    node 0.
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
!  Parameters:
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, integer PARENT(NNODE), is the parent node of each node.
!    The node with parent 0 is the root.
!
!    Output, integer CODE(2*NNODE), the code for the tree.
!
  implicit none

  integer nnode

  integer code(2*nnode)
  integer father
  integer i
  integer k
  integer parent(nnode)
  integer son
!
!  Find the root.
!
  father = 0
  do i = 1, nnode
    if ( parent(i) == 0 ) then
      k = 1
      code(1) = 1
      father = i
      exit
    end if
  end do

  if ( father == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TREE_ROOTED_CODE - Fatal error!'
    write ( *, '(a)' ) '  Could not find the root.'
    stop
  end if

  do while ( father /= 0 ) 

    k = k + 1
    code(k) = 0

    do son = 1, nnode
      if ( parent(son) == father ) then
        code(k) = 1
        father = son
        exit
      end if
    end do

    if ( code(k) == 0 ) then
      parent(father) = - parent(father)
      father = - parent(father)
    end if

  end do

  parent(1:nnode) = - parent(1:nnode)

  return
end
subroutine tree_rooted_code_compare ( nnode, npart, code1, code2, result )

!*****************************************************************************80
!
!! TREE_ROOTED_CODE_COMPARE compares a portion of the code for two rooted trees.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, integer NPART, the number of nodes for which the code
!    has been determined.  This determines the portion of the codes to be
!    compared.  We expect 0 <= NPART <= NNODE.
!
!    Input, integer CODE1(2*NNODE), CODE2(2*NNODE), the two 
!    rooted tree codes to be compared.
!
!    Output, integer RESULT, the result of the comparison.
!    -1, CODE1 < CODE2,
!     0, CODE1 = CODE2,
!    +1, CODE1 > CODE2.
!
  implicit none

  integer nnode

  integer code1(2*nnode)
  integer code2(2*nnode)
  integer i
  integer ihi
  integer npart
  integer result

  result = 0

  if ( npart <= 0 ) then
    return
  end if

  ihi = 2 * min ( npart, nnode )

  do i = 1, ihi

    if ( code1(i) < code2(i) ) then
      result = -1
      return
    else if ( code2(i) < code1(i) ) then
      result = +1
      return
    end if

  end do

  return
end
subroutine tree_rooted_depth ( nnode, parent, depth, depth_node )

!*****************************************************************************80
!
!! TREE_ROOTED_DEPTH returns the depth of a rooted tree.
!
!  Discussion:
!
!    The depth of any node of a rooted tree is the number of edges in 
!    the shortest path from the root to the node.
!
!    The depth of the rooted tree is the maximum of the depths
!    of all the nodes.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, integer PARENT(NNODE), is the parent node of each node.
!    The node with parent 0 is the root.
!
!    Output, integer DEPTH, the depth of the tree.
!
!    Output, integer DEPTH_NODE(NNODE), the depth of each node.
!
  implicit none

  integer nnode

  integer depth
  integer depth_node(nnode)
  integer i
  integer j
  integer parent(nnode)
  integer root
!
!  Find the root.
!
  root = 0
  do i = 1, nnode
    if ( parent(i) == 0 ) then
      root = i
      exit
    end if
  end do

  if ( root == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TREE_ROOTED_DEPTH - Fatal error!'
    write ( *, '(a)' ) '  Could not find the root.'
    stop
  end if
!
!  Determine the depth of each node by moving towards the node.
!  If you reach a node whose depth is already known, stop early.
!
  depth_node(1:nnode) = 0

  do i = 1, nnode

    j = i

    do while ( j /= root )

      depth_node(i) = depth_node(i) + 1
      j = parent(j)

      if ( 0 < depth_node(j) ) then
        depth_node(i) = depth_node(i) + depth_node(j)
        exit
      end if

    end do

  end do
!
!  Determine the maximum depth.
!
  depth = maxval ( depth_node(1:nnode) )

  return
end
subroutine tree_rooted_enum ( nnode, ntree )

!*****************************************************************************80
!
!! TREE_ROOTED_ENUM counts the number of unlabeled rooted trees.
!
!  Example:
!
!    Input    Output
!
!      1         1
!      2         1
!      3         2
!      4         4
!      5         9
!      6        20
!      7        48
!      8       115
!      9       286
!     10       719
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 August 2000
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
!    Output, integer NTREE(NNODE).  NTREE(I) is the number of 
!    rooted, unlabeled trees on I nodes, for I = 1, 2, ... NNODE.
!
  implicit none

  integer nnode

  integer i
  integer id
  integer isum
  integer itd
  integer j
  integer nlast
  integer ntree(nnode)

  ntree(1) = 1
 
  do nlast = 2, nnode
 
    isum = 0
 
    do id = 1, nlast - 1
 
      i = nlast
      itd = ntree(id) * id
 
      do j = 1, nlast - 1

        i = i - id

        if ( i <= 0 ) then
          exit
        end if

        isum = isum + ntree(i) * itd

      end do
 
    end do
 
    ntree(nlast) = isum / ( nlast - 1 )
 
  end do

  return
end
subroutine tree_rooted_random ( nnode, seed, ntree, itree )

!*****************************************************************************80
!
!! TREE_ROOTED_RANDOM selects a random unlabeled rooted tree.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 June 2013
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
!    Input/output, integer SEED, a seed for the random 
!    number generator.
!
!    Output, integer NTREE(NNODE).  NTREE(I) is the number of 
!    rooted, unlabeled trees on I nodes, for I = 1, 2, ... NNODE.
!
!    Output, integer ITREE(NNODE).  (I,ITREE(I)) is the I-th edge
!    of the output tree for I = 2,NNODE.  ITREE(1)=0.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nnode

  integer i
  integer id
  integer is1
  integer is2
  integer itd
  integer itree(nnode)
  integer iz
  integer j
  integer l
  integer ll
  integer ls
  integer m
  integer ntree(nnode)
  integer nval
  real ( kind = rk ) r
  real ( kind = rk ) r8_uniform_01
  integer seed
  integer stack(2,nnode)

  if ( nnode <= 0  ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TREE_ROOTED_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  NNODE = ', nnode
    write ( *, '(a)' ) '  but NNODE must be at least 1.'
    stop
  end if
!
!  Compute a table of the number of such trees for a given number of nodes.
!
  call tree_rooted_enum ( nnode, ntree )
!
!  Now select one such tree at random.
!
  l = 0

  nval = nnode
  is1 = 0
  is2 = 0
  
  do

    do while ( 2 < nval )
 
      r = r8_uniform_01 ( seed )

      iz = int ( ( nval - 1 ) * ntree(nval) * r )

      id = 0
  
      id = id + 1
      itd = id * ntree(id)
      m = nval
      j = 0
 
      do
 
        j = j + 1
        m = m - id

        if ( m < 1 ) then
          id = id + 1
          itd = id * ntree(id)
          m = nval
          j = 0
          cycle
        end if

        iz = iz - ntree(m) * itd

        if ( iz < 0 ) then
          exit
        end if

      end do

      is1 = is1 + 1
      stack(1,is1) = j
      stack(2,is1) = id
      nval = m
 
    end do
 
    itree(is2+1) = l
    l = is2 + 1
    is2 = is2 + nval

    if ( 1 < nval ) then
      itree(is2) = is2 - 1
    end if
 
    do
 
      nval = stack(2,is1)
 
      if ( nval /= 0 ) then
        stack(2,is1) = 0
        exit
      end if
 
      j = stack(1,is1)
      is1 = is1 - 1
      m = is2 - l + 1
      ll = itree(l)
      ls = l + ( j - 1 ) * m - 1
 
      if ( j /= 1 ) then
        do i = l, ls
          itree(i+m) = itree(i) + m
          if ( mod(i-l,m) == 0 ) then
            itree(i+m) = ll
          end if
        end do
      end if
 
      is2 = ls + m
 
      if ( is2 == nnode ) then
        return
      end if

      l = ll
 
    end do

  end do
 
  return
end
subroutine vec_next ( n, ibase, iarray, more )

!*****************************************************************************80
!
!! VEC_NEXT generates all N-vectors of integers modulo a given base.
!
!  Discussion:
!
!    The items are produced one at a time.
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
!    Input, integer N, the size of the vectors to be used.
!
!    Input, integer IBASE, the base to be used.  IBASE = 2 will
!    give vectors of 0's and 1's, for instance.
!
!    Input/output, integer IARRAY(N).  On each return from VECNEX,
!    IARRAY will contain entries in the range 0 to IBASE-1.
!
!    Input/output, logical MORE.  Set this variable .FALSE. before
!    the first call.  Normally, MORE will be returned .TRUE. but
!    once all the vectors have been generated, MORE will be
!    reset .FALSE. and you should stop calling the program.
!
  implicit none

  integer n

  integer i
  integer iarray(n)
  integer ibase
  integer, save :: kount
  integer, save :: last
  logical more
  integer nn

  if ( .not. more ) then
 
    kount = 1
    last = ibase ** n
    more = .true.
    iarray(1:n) = 0
 
  else
 
    kount = kount + 1

    if ( kount == last ) then
      more = .false.
    end if

    iarray(n) = iarray(n) + 1
 
    do i = 1, n

      nn = n - i

      if ( iarray(nn+1) < ibase ) then
        return
      end if

      iarray(nn+1) = 0

      if ( nn /= 0 ) then
        iarray(nn) = iarray(nn) + 1
      end if

    end do
 
  end if
 
  return
end
subroutine vec_random ( n, base, seed, iarray )

!*****************************************************************************80
!
!! VEC_RANDOM selects a random N-vector of integers modulo a given base.
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
!    Input/output, integer SEED, a seed for the random 
!    number generator.
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
  integer seed

  do i = 1, n
    ival = i4_uniform_ab ( 0, base-1, seed )
    iarray(i) = ival
  end do
 
  return
end
