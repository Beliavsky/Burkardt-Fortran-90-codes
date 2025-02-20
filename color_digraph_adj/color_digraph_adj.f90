subroutine color_digraph_adj_color_count ( adj, nnode, mcolor, ncolor )

!*****************************************************************************80
!
!! color_digraph_adj_color_count() counts the number of colors in a color digraph.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 March 2023
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ADJ(NNODE,NNODE), the adjacency information.  
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge from node I to node J.
!
!    Input, integer NNODE, the number of nodes.
!
!    Output, integer MCOLOR, the maximum color index.
!
!    Output, integer NCOLOR, the number of colors.
!
  implicit none

  integer nnode

  integer adj(nnode,nnode)
  integer colors(nnode)
  integer i
  integer mcolor
  integer ncolor

  mcolor = 0
  do i = 1, nnode
    mcolor = max ( mcolor, adj(i,i) )
  end do

  do i = 1, nnode
    colors(i) = adj(i,i)
  end do

  call i4vec_sort_heap_a ( nnode, colors )

  call i4vec_uniq ( nnode, colors, ncolor )

  return
end
subroutine color_digraph_adj_degree ( adj, nnode, indegree, outdegree )

!*****************************************************************************80
!
!! color_digraph_adj_degree() computes the indegree and outdegree of each node.
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
!    10 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer ADJ(NNODE,NNODE), the adjacency information 
!    for graph 1.  ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is 
!    positive if there is an edge from node I to node J. 
!
!    integer NNODE, the number of nodes.
!
!  Output:
!
!    integer INDEGREE(NNODE), OUTDEGREE(NNODE), 
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
      if ( i /= j ) then
        outdegree(i) = outdegree(i) + adj(i,j)
        indegree(j) = indegree(j) + adj(i,j)
      end if
    end do
  end do

  return
end
subroutine color_digraph_adj_degree_sequence ( adj, nnode, in_seq, out_seq )

!*****************************************************************************80
!
!! color_digraph_adj_degree_sequence() computes the degree sequence of a color digraph.
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
!    04 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ADJ(NNODE,NNODE), the adjacency information.  
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge from node I to node J.
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

  call color_digraph_adj_degree ( adj, nnode, in_seq, out_seq )

  call i4vec2_sort_d ( nnode, out_seq, in_seq )

  return
end
subroutine color_digraph_adj_edge_count ( adj, nnode, nedge )

!*****************************************************************************80
!
!! color_digraph_adj_edge_count() counts the number of edges in a color digraph.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ADJ(NNODE,NNODE), the adjacency information.  
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge from node I to node J.
!
!    Input, integer NNODE, the number of nodes.
!
!    Output, integer NEDGE, the number of edges.
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

      if ( i /= j ) then
        nedge = nedge + adj(i,j)
      end if

    end do
  end do

  return
end
subroutine color_digraph_adj_example_cube ( adj )

!*****************************************************************************80
!
!! color_digraph_adj_example_cube() sets up the cube color digraph.
!
!  Diagram:
!
!
!    8B----<-----3B
!    |\          /|\
!    | A        V | |
!    |  \      /  | |
!    |  4R-->-7R  | |
!    |   |     |  | |
!    A   A     V  V A
!    |   |     |  | |
!    |   5B-<-2G  | |
!    |  /      \  | |
!    | A        A | |
!    |/          \|/
!    1G----->----6B
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Output:
!
!    integer ADJ(8,8), the adjacency information.  
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge from node I to node J.
!
  implicit none

  integer, parameter :: nnode = 8

  integer, parameter :: BLUE = 1
  integer, parameter :: GREEN = 2
  integer, parameter :: RED = 3

  integer adj(nnode,nnode)

  adj(1:nnode,1:nnode) = 0

  adj(1,1) = GREEN
  adj(1,5) = 1
  adj(1,6) = 1
  adj(1,8) = 1

  adj(2,2) = GREEN
  adj(2,5) = 1

  adj(3,3) = BLUE
  adj(3,6) = 1
  adj(3,7) = 1
  adj(3,8) = 1

  adj(4,4) = RED
  adj(4,7) = 1
  adj(4,8) = 1

  adj(5,5) = BLUE
  adj(5,4) = 1

  adj(6,6) = BLUE
  adj(6,2) = 1
  adj(6,3) = 1

  adj(7,7) = RED
  adj(7,2) = 1

  adj(8,8) = BLUE

  return
end
subroutine color_digraph_adj_example_octo ( example, adj )

!*****************************************************************************80
!
!! color_digraph_adj_example_octo() sets up an 8 node example color digraph.
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
!    Input, integer EXAMPLE, should be between 1 and 60, and 
!    indicates which example graph to pick.
!
!    Output, integer ADJ(8,8), the adjacency information.
!    ADJ(I,I) is the color of node I.
!    ADJ(I,J) is 1 if nodes I and J are adjacent and 0 otherwise.
!
  implicit none

  integer, parameter :: nnode = 8

  integer, parameter :: BLUE = 1
  integer, parameter :: GREEN = 2
  integer, parameter :: RED = 3
  integer, parameter :: YELLOW = 4

  integer adj(nnode,nnode)
  integer example
  integer i
  integer j
  integer msave
  integer nsave

  example = mod ( example - 1, 60 ) + 1
  msave = ( example - 1 ) / 12 + 1
  nsave = mod ( example - 1, 12 ) + 1

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

  if ( msave == 1 ) then

    adj(1,1) = RED
    adj(2,2) = RED
    adj(3,3) = RED
    adj(4,4) = BLUE
    adj(5,5) = BLUE
    adj(6,6) = BLUE
    adj(7,7) = GREEN
    adj(8,8) = GREEN

  else if ( msave == 2 ) then

    adj(1,1) = RED
    adj(2,2) = RED
    adj(3,3) = RED
    adj(4,4) = BLUE
    adj(5,5) = BLUE
    adj(6,6) = BLUE
    adj(7,7) = GREEN
    adj(8,8) = YELLOW

  else if ( msave == 3 ) then

    adj(1,1) = RED
    adj(2,2) = RED
    adj(3,3) = RED
    adj(4,4) = BLUE
    adj(5,5) = BLUE
    adj(6,6) = BLUE
    adj(7,7) = YELLOW
    adj(8,8) = YELLOW

  else if ( msave == 4 ) then

    adj(1,1) = RED
    adj(2,2) = RED
    adj(3,3) = RED
    adj(4,4) = BLUE
    adj(5,5) = BLUE
    adj(6,6) = GREEN
    adj(7,7) = GREEN
    adj(8,8) = GREEN

  else if ( msave == 5 ) then

    adj(1,1) = RED
    adj(2,2) = BLUE
    adj(3,3) = RED
    adj(4,4) = GREEN
    adj(5,5) = BLUE
    adj(6,6) = RED
    adj(7,7) = BLUE
    adj(8,8) = GREEN

  end if
!
!  Now permute the graph.
!
  call i4mat_perm_random ( nnode, adj )

  return
end
subroutine color_digraph_adj_print ( adj, nnode, title )

!*****************************************************************************80
!
!! color_digraph_adj_print() prints the adjacency matrix of a color digraph.
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
!    Input, integer ADJ(NNODE,NNODE), the adjacency information.  
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge from node I to node J.
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
  integer k
  character ( len = 80 ) string
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, nnode

    do j = 1, nnode

      k = ( j - 1 ) * 3 + 1
      write ( string(k:k+2), '(i3)' ) adj(i,j)

    end do

    write ( *, '(i2,2x,a)' ) i, string(1:3*nnode)

  end do

  return
end
subroutine color_digraph_adj_random ( nnode, ncolor, nedge, adj )

!*****************************************************************************80
!
!! color_digraph_adj_random() generates a random color graph.
!
!  Discussion:
!
!    Each node is assumed to have an associated color, between 1 and NCOLOR,
!    which will be chosen at random.
!
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge from node I to node J.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, integer NCOLOR, the number of colors available.  
!
!    Input, integer NEDGE, the number of edges, which must be 
!    between 0 and NNODE*(NNODE-1).
!
!    Output, integer ADJ(NNODE,NNODE), the adjacency information.  
!
  implicit none

  integer ncolor
  integer nedge
  integer nnode

  integer adj(nnode,nnode)
  integer color
  integer i
  integer i4_uniform_ab
  integer iwork(nedge)
  integer j
  integer k
  integer l
  integer maxedge
  integer perm(ncolor)
  integer subset(ncolor)

  if ( nnode <= 0  ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'color_digraph_adj_random(): Fatal error!'
    write ( *, '(a,i8)' ) '  NNODE = ', nnode
    write ( *, '(a)' ) '  but NNODE must be at least 1.'
    stop 1
  end if

  maxedge = nnode * ( nnode - 1 )

  if ( nedge < 0 .or. maxedge < nedge ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COLOR_DIGRAPH_ADJ_RANDOM(): Fatal error!'
    write ( *, '(a,i8)' ) '  NEDGE = ', nedge
    write ( *, '(a)' ) '  but NEDGE must be at least 0, and '
    write ( *, '(a,i8)' ) '  no more than ', maxedge
    stop 1
  end if
!
!  Start with no edges, no colors.
!
  adj(1:nnode,1:nnode) = 0
!
!  Choose the colors.
!
  call ksub_random ( nnode, ncolor, subset )

  call perm_random ( ncolor, perm )

  do color = 1, ncolor
    i = subset(perm(color))
    adj(i,i) = color
  end do

  do i = 1, nnode
    if ( adj(i,i) == 0 ) then
      color = i4_uniform_ab ( 1, ncolor )
      adj(i,i) = color
    end if
  end do
!
!  Pick a random NEDGE subset.
!
  call ksub_random ( maxedge, nedge, iwork )
!
!  Mark the potential edges that were chosen.
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
!    Original Fortran77 version by Albert Nijenhuis, Herbert Wilf.
!    Fortran90 version by John Burkardt.
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
subroutine perm_random ( n, iarray )

!*****************************************************************************80
!
!! perm_random() selects a random permutation of N objects.
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
!    Original Fortran77 version by Albert Nijenhuis, Herbert Wilf.
!    Fortran90 version by John Burkardt.
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
!    Input, integer N, the number of objects to be permuted.
!
!    Output, integer IARRAY(N), the random permutation.
!
  implicit none

  integer n

  integer i
  integer i4_uniform_ab
  integer iarray(n)
  integer j

  call i4vec_indicator ( n, iarray )

  do i = 1, n
    j = i4_uniform_ab ( i, n )
    call i4_swap ( iarray(j), iarray(i) )
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
!    Original Fortran77 version by Albert Nijenhuis, Herbert Wilf.
!    Fortran90 version by John Burkardt.
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

