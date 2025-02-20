subroutine color_graph_adj_color_count ( adj, nnode, mcolor, ncolor )

!*****************************************************************************80
!
!! color_graph_adj_color_count() counts the number of colors in a color graph.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ADJ(NNODE,NNODE), the adjacency information.  
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge between node I and node J.
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
subroutine color_graph_adj_color_sequence ( adj, nnode, seq )

!*****************************************************************************80
!
!! color_graph_adj_color_sequence() computes the color sequence of a color graph.
!
!  Discussion:
!
!    The color sequence of a color graph is constructed by computing the
!    color of each node, and then ordering these values in decreasing order.
!
!    If two color graphs are isomorphic, they must have the same color sequence.
!
!    If two color graphs have different color sequences, they cannot be
!    isomorphic.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ADJ(NNODE,NNODE), the adjacency information.  
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge between node I and node J.
!
!    Input, integer NNODE, the number of nodes.
!
!    Output, integer SEQ(NNODE), the color sequence.
!
  implicit none

  integer nnode

  integer adj(nnode,nnode)
  integer i
  integer seq(nnode)

  do i = 1, nnode
    seq(i) = adj(i,i)
  end do

  call i4vec_sort_heap_d ( nnode, seq )

  return
end
subroutine color_graph_adj_connect_random ( nnode, nedge, ncolor, adj )

!*****************************************************************************80
!
!! color_graph_adj_connect_random() generates a random connected color graph.
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
!    between NNODE-1 and (NNODE*(NNODE-1))/2.  
!
!    Input, integer NCOLOR, the number of colors available to 
!    choose for the nodes.  NCOLOR must be at least 1, and no more than NNODE.
!
!    Output, integer ADJ(NNODE,NNODE), the adjacency information.  
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge between node I and node J.
!
  implicit none

  integer ncolor
  integer nnode
  integer nedge

  integer adj(nnode,nnode)
  integer code(nnode-2)
  integer color
  integer i
  integer i4_uniform_ab
  integer inode(nnode-1)
  integer iwork(nedge)
  integer j
  integer jnode(nnode-1)
  integer k
  integer l
  integer maxedge
  integer nchoice
  integer nchoose
  integer perm(ncolor)
  integer subset(ncolor)
!
!  Check.
!
  if ( nnode <= 0  ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COLOR_GRAPH_ADJ_CONNECT_RANDOM(): Fatal error!'
    write ( *, '(a,i8)' ) '  NNODE = ', nnode
    write ( *, '(a)' ) '  but NNODE must be at least 1.'
    stop 1
  end if

  maxedge = ( nnode * ( nnode - 1 ) ) / 2

  if ( nedge < nnode-1 .or. maxedge < nedge ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COLOR_GRAPH_ADJ_CONNECT_RANDOM(): Fatal error!'
    write ( *, '(a,i8)' ) '  NEDGE = ', nedge
    write ( *, '(a)' ) '  but NEDGE must be at least 0, and '
    write ( *, '(a,i8)' ) '  no more than ', maxedge
    stop 1
  end if

  if ( ncolor < 1 .or. nnode < ncolor ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COLOR_GRAPH_ADJ_CONNECT_RANDOM(): Fatal error!'
    write ( *, '(a,i8)' ) '  NCOLOR = ', ncolor
    write ( *, '(a)' ) '  but NCOLOR must be at least 1, and '
    write ( *, '(a,i8)' ) '  no more than ', nnode
    stop 1
  end if
!
!  Initialize the adjacency matrix.
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
!  Pick a random tree.
!
  call tree_arc_random ( nnode, code, inode, jnode )
!
!  Convert information to adjacency form.
!
  call graph_arc_to_graph_adj ( nnode-1, inode, jnode, adj, nnode )
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
      if ( adj(i,j) /= 0 ) then
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
subroutine color_graph_adj_degree ( adj, nnode, degree )

!*****************************************************************************80
!
!! color_graph_adj_degree() computes the degree of each node.
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
!    10 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ADJ(NNODE,NNODE), the adjacency information.  
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge between node I and node J.
!
!    Input, integer NNODE, the number of nodes.
!
!    Output, integer DEGREE(NNODE), the degree of the nodes.
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
      if ( i /= j ) then
        if ( adj(i,j) /= 0 ) then
          degree(i) = degree(i) + adj(i,j)
        end if
      end if
    end do
  end do

  return
end
subroutine color_graph_adj_degree_sequence ( adj, nnode, seq )

!*****************************************************************************80
!
!! color_graph_adj_degree_sequence() computes the degree sequence of a color graph.
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
!    10 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ADJ(NNODE,NNODE), the adjacency information.  
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge between node I and node J.
!
!    Input, integer NNODE, the number of nodes.
!
!    Output, integer SEQ(NNODE), the degree sequence.
!
  implicit none

  integer nnode

  integer adj(nnode,nnode)
  integer seq(nnode)

  call color_graph_adj_degree ( adj, nnode, seq )

  call i4vec_sort_heap_d ( nnode, seq )

  return
end
subroutine color_graph_adj_edge_count ( adj, nnode, nedge )

!*****************************************************************************80
!
!! color_graph_adj_edge_count() counts the number of edges in a color graph.
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
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge between node I and node J.
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

  nedge = nedge / 2

  return
end
subroutine color_graph_adj_example_bush ( adj )

!*****************************************************************************80
!
!! color_graph_adj_example_bush() sets up the bush color graph.
!
!  Diagram:
!
!        6G  3R
!        |   |
!    1B--4G--5W--2R
!        |
!        7W
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
!  Parameters:
!
!    Output, integer ADJ(NNODE,NNODE), the adjacency information.  
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge between node I and node J.
!
  implicit none

  integer, parameter :: nnode = 7

  integer, parameter :: BLUE = 1
  integer, parameter :: GREEN = 2
  integer, parameter :: RED = 3
  integer, parameter :: WHITE = 4

  integer adj(nnode,nnode)

  adj(1:nnode,1:nnode) = 0

  adj(1,1) = BLUE
  adj(1,4) = 1

  adj(2,2) = RED
  adj(2,5) = 1

  adj(3,3) = RED
  adj(3,5) = 1

  adj(4,1) = 1
  adj(4,4) = GREEN
  adj(4,5) = 1
  adj(4,6) = 1
  adj(4,7) = 1

  adj(5,2) = 1
  adj(5,3) = 1
  adj(5,4) = 1
  adj(5,5) = WHITE

  adj(6,4) = 1
  adj(6,6) = GREEN

  adj(7,4) = 1
  adj(7,7) = WHITE

  return
end
subroutine color_graph_adj_example_cube ( adj )

!*****************************************************************************80
!
!! color_graph_adj_example_cube() sets up the cube color graph.
!
!  Diagram:
!
!      4R----7R
!     /|    /|
!    8B----3B|
!    | |   | |
!    | 5B--|-2G
!    |/    |/
!    1G----6B
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
!  Parameters:
!
!    Output, integer ADJ(NNODE,NNODE), the adjacency information.  
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge between node I and node J.
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
  adj(2,6) = 1
  adj(2,7) = 1

  adj(3,3) = BLUE
  adj(3,6) = 1
  adj(3,7) = 1
  adj(3,8) = 1

  adj(4,4) = RED
  adj(4,5) = 1
  adj(4,7) = 1
  adj(4,8) = 1

  adj(5,5) = BLUE
  adj(5,1) = 1
  adj(5,2) = 1
  adj(5,4) = 1

  adj(6,6) = BLUE
  adj(6,1) = 1
  adj(6,2) = 1
  adj(6,3) = 1

  adj(7,7) = RED
  adj(7,2) = 1
  adj(7,3) = 1
  adj(7,4) = 1

  adj(8,8) = BLUE
  adj(8,1) = 1
  adj(8,3) = 1
  adj(8,4) = 1

  return
end
subroutine color_graph_adj_example_octo ( example, adj )

!*****************************************************************************80
!
!! color_graph_adj_example_octo() sets up an 8 node example color graph.
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
!    05 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer EXAMPLE, should be between 1 and 35, and 
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

  example = mod ( example - 1, 35 ) + 1
  msave = ( ( example - 1 ) / 7 ) + 1
  nsave = mod ( example - 1, 7 ) + 1

  adj(1:nnode,1:nnode) = 0

  do i = 1, nnode
    j = i + 1
    if ( nnode < j ) then
      j = j - nnode
    end if

    adj(i,j) = 1
    adj(j,i) = 1

  end do
!
!  Underlying graph 1.
!
  if ( nsave == 1 ) then

    adj(1,6) = 1
    adj(6,1) = 1
    adj(2,5) = 1
    adj(5,2) = 1
    adj(3,8) = 1
    adj(8,3) = 1
    adj(4,7) = 1
    adj(7,4) = 1
!
!  Underlying graph 2.
!
  else if ( nsave == 2 ) then

    adj(1,6) = 1
    adj(6,1) = 1
    adj(2,8) = 1
    adj(8,2) = 1
    adj(3,5) = 1
    adj(5,3) = 1
    adj(4,7) = 1
    adj(7,4) = 1
!
!  Underlying graph 3.
!
  else if ( nsave == 3 ) then

    adj(1,5) = 1
    adj(5,1) = 1
    adj(2,6) = 1
    adj(6,2) = 1
    adj(3,7) = 1
    adj(7,3) = 1
    adj(4,8) = 1
    adj(8,4) = 1
!
!  Underlying graph 4.
!
  else if ( nsave == 4 ) then

    adj(1,3) = 1
    adj(3,1) = 1
    adj(2,4) = 1
    adj(4,2) = 1
    adj(5,7) = 1
    adj(7,5) = 1
    adj(6,8) = 1
    adj(8,6) = 1
!
!  Underlying graph 5.
!
  else if ( nsave == 5 ) then

    adj(1,4) = 1
    adj(4,1) = 1
    adj(2,6) = 1
    adj(6,2) = 1
    adj(3,8) = 1
    adj(8,3) = 1
    adj(5,7) = 1
    adj(7,5) = 1
!
!  Underlying graph 6.
!
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
!
!  Underlying graph 7.
!
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
subroutine color_graph_adj_example_twig ( adj )

!*****************************************************************************80
!
!! color_graph_adj_example_twig() sets up the twig color graph.
!
!  Diagram:
!
!    1R---2R---3B
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
!  Parameters:
!
!    Output, integer ADJ(3,3), the adjacency information.  
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge between node I and node J.
!
  implicit none

  integer, parameter :: nnode = 3

  integer, parameter :: BLUE = 1
  integer, parameter :: RED = 3

  integer adj(nnode,nnode)

  adj(1:nnode,1:nnode) = 0

  adj(1,1) = RED
  adj(1,2) = 1

  adj(2,1) = 1
  adj(2,2) = RED
  adj(2,3) = 1

  adj(3,2) = 1
  adj(3,3) = BLUE

  return
end
subroutine color_graph_adj_print ( adj, nnode, title )

!*****************************************************************************80
!
!! color_graph_adj_print() prints the adjacency matrix of a color graph.
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
!    if there is an edge between node I and node J.
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

  if ( len_trim ( title ) /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  do i = 1, nnode

    do j = 1, nnode

      k = (j-1) * 3 + 1
      write ( string(k:k+2), '(i3)' ) adj(i,j)

    end do

    write ( *, '(i2,2x,a)' ) i, string(1:3*nnode)

  end do

  return
end
subroutine color_graph_adj_random ( nnode, ncolor, nedge, adj )

!*****************************************************************************80
!
!! color_graph_adj_random() generates a random color graph.
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
!    Input, integer NCOLOR, the number of colors available to 
!    choose for the nodes.  NCOLOR must be at least 1, and no more than NNODE.
!
!    Input, integer NEDGE, the number of edges, which must be 
!    between 0 and (NNODE*(NNODE-1))/2.
!
!    Output, integer ADJ(NNODE,NNODE), the adjacency information.  
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge between node I and node J.
!
  implicit none

  integer nnode
  integer nedge

  integer adj(nnode,nnode)
  integer color
  integer i
  integer i4_uniform_ab
  integer iwork(nedge)
  integer j
  integer k
  integer l
  integer maxedge
  integer ncolor
  integer perm(ncolor)
  integer subset(ncolor)

  if ( nnode <= 0  ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COLOR_GRAPH_ADJ_RANDOM(): Fatal error!'
    write ( *, '(a,i8)' ) '  NNODE = ', nnode
    write ( *, '(a)' ) '  but NNODE must be at least 1.'
    stop 1
  end if

  maxedge = ( nnode * ( nnode - 1 ) ) / 2

  if ( nedge < 0 .or. maxedge < nedge ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COLOR_GRAPH_ADJ_RANDOM(): Fatal error!'
    write ( *, '(a,i8)' ) '  NEDGE = ', nedge
    write ( *, '(a)' ) '  but NEDGE must be at least 0, and '
    write ( *, '(a,i8)' ) '  no more than ', maxedge
    stop 1
  end if

  if ( ncolor < 1 .or. nnode < ncolor ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COLOR_GRAPH_ADJ_RANDOM(): Fatal error!'
    write ( *, '(a)' ) '  Illegal value of NCOLOR.'
    stop 1
  end if
!
!  Start out with no edges and no colors.
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
!  Pick a random NEDGE subset of (N*(N-1))/2.
!
  call ksub_random ( maxedge, nedge, iwork )
!
!  The (n*(n-1))/2 spots in the superdiagonal are numbered as follows:
!
!  * 1  2   3  ...  n-1   n
!  * * n+1 n+2 ... 2n-2  2n-1
!  ...
!  * *  *   *  ...   *   (n*(n-1))/2
!  * *  *   *  ...   *    *
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

