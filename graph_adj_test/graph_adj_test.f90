program main

!*****************************************************************************80
!
!! graph_adj_test() tests graph_adj().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 March 2023
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_adj_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  graph_adj() implements graph algorithms.'

  call graph_adj_bipartite_random_test ( )
  call graph_adj_breadth_first_test ( )
  call graph_adj_block_test ( )
  call graph_adj_color_next_test ( )
  call graph_adj_complement_test ( )
  call graph_adj_connect_random_test ( )
  call graph_adj_cycle_test ( )
  call graph_adj_degree_test ( )
  call graph_adj_degree_max_test ( )
  call graph_adj_degree_sequence_test ( )
  call graph_adj_depth_first_test ( )
  call graph_adj_depth_first_2_test ( )
  call graph_adj_edge_count_test ( )
  call graph_adj_edge_select_test ( )
  call graph_adj_edges_random_test()
  call graph_adj_eigen_test ( )
  call test036 ( )
  call test0365 ( )
  call graph_adj_ham_next_brute_test ( )
  call graph_adj_is_bipartite_test ( )
  call graph_adj_is_edge_connected_test ( )
  call graph_adj_is_node_connected_test ( )
  call graph_adj_is_tree_test ( )
  call graph_adj_random_test ( )
  call graph_adj_reduce_test ( )
  call graph_adj_span_tree_test ( )
  call graph_adj_span_tree_enum_test ( )
  call graph_adj_transitive_closure_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_adj_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine graph_adj_bipartite_random_test ( )

!*****************************************************************************80
!
!! graph_adj_bipartite_random_test() tests graph_adj_bipartite_random().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: nnode1 = 4
  integer, parameter :: nnode2 = 6
  integer, parameter :: nnode = nnode1 + nnode2

  integer adj(nnode,nnode)
  integer nedge
  integer nedge2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_adj_bipartite_random_test():'
  write ( *, '(a)' ) '  graph_adj_bipartite_random_() returns a random '
  write ( *, '(a)' ) '  bipartite graph;'
  write ( *, '(a)' ) ' '

  write ( *, '(a,i8)' ) '  Number of nodes in set 1 is ', nnode1
  write ( *, '(a,i8)' ) '  Number of nodes in set 2 is ', nnode2

  call graph_adj_bipartite_random ( nnode1, nnode2, nedge, adj )

  call graph_adj_print ( adj, nnode, '  The graph:' )
!
!  Count the edges.
!
  call graph_adj_edge_count ( adj, nnode, nedge2 )

  write ( *, '(a,i8)' ) '  Total number of edges is   ', nedge
  write ( *, '(a,i8)' ) '  Counted number of edges is ', nedge2

  return
end 
subroutine graph_adj_breadth_first_test ( )

!*****************************************************************************80
!
!! graph_adj_breadth_first_test() tests graph_adj_breadth_first().
!
!  Discussion:
!
!    This example is from page 22 of the Gibbons reference.
!
!    The correct result is
!
!    Node  Idad   Ideep Iorder
!
!     1      0       1    1
!     2      1       2    2
!     3      1       2    3
!     4      1       2    4
!     5      1       2    5
!     6      1       2    6
!     7      1       2    7
!     8      1       2    8
!     9      0       3    9
!    10      9       4   10
!    11     10       5   12
!    12     10       5   13
!    13      9       4   11
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan Gibbons,
!    Algorithmic Graph Theory,
!    Cambridge University Press, 1985
!    ISBN 0-521-28881-9
!
  implicit none

  integer, parameter :: nnode = 13

  integer i
  integer adj(nnode,nnode)
  integer dad(nnode)
  integer deep(nnode)
  integer order(nnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_adj_breadth_first_test()'
  write ( *, '(a)' ) '  graph_adj_breadth_first() sets up a breadth-first'
  write ( *, '(a)' ) '  traversal of a graph.'
  write ( *, '(a)' ) ' '
 
  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1
  adj(1,3) = 1
  adj(1,4) = 1
  adj(1,5) = 1
  adj(1,6) = 1
  adj(1,7) = 1
  adj(1,8) = 1
 
  adj(2,1) = 1
  adj(2,5) = 1
  adj(2,6) = 1
  adj(2,8) = 1
 
  adj(3,1) = 1
  adj(3,4) = 1
  adj(3,7) = 1
 
  adj(4,1) = 1
  adj(4,3) = 1
 
  adj(5,1) = 1
  adj(5,2) = 1
 
  adj(6,1) = 1
  adj(6,2) = 1
 
  adj(7,1) = 1
  adj(7,3) = 1
 
  adj(8,1) = 1
  adj(8,2) = 1
 
  adj(9,10) = 1
  adj(9,13) = 1
 
  adj(10,9) = 1
  adj(10,11) = 1
  adj(10,12) = 1
  adj(10,13) = 1
 
  adj(11,10) = 1
  adj(11,12) = 1
 
  adj(12,10) = 1
  adj(12,11) = 1
 
  adj(13,9) = 1
  adj(13,10) = 1

  call graph_adj_print ( adj, nnode, '  The graph:' )
 
  call graph_adj_breadth_first ( adj, nnode, dad, deep, order )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, dad(i), deep(i), order(i)'
  write ( *, '(a)' ) ' '
 
  do i = 1, nnode
    write ( *, '(4i8)' )  i, dad(i), deep(i), order(i)
  end do
 
  return
end
subroutine graph_adj_block_test ( )

!*****************************************************************************80
!
!! graph_adj_block_test() tests graph_adj_block().
!
!  Discussion:
!
!    The correct result is
!
!    3 blocks
!
!    Node  Idad  Iorder
!
!       1     0      -1
!       2     1       2
!       3     4       5
!       4     1      -4
!       5     4       6
!       6     2       3
!
!    Revised adjacency matrix:
!
!      0 1 0 3 3 1
!      1 0 0 0 0 1
!      0 0 0 2 0 0
!      3 0 2 0 3 0
!      3 0 0 3 0 0
!      1 1 0 0 0 0
!
!    The three blocks are defined by the edges:
!
!      1: (6,1), (2,6), (1,2)
!      2: (4,3)
!      3: (1,4), (4,5), (5,1)
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: nnode = 6

  integer adj(nnode,nnode)

  integer dad(nnode)
  integer order(nnode)
  integer nblock
  integer stack(nnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_adj_block_test():'
  write ( *, '(a)' ) '  graph_adj_block() finds the blocks in a graph.'
  write ( *, '(a)' ) ' '
 
  adj(1:nnode,1:nnode) = 0
 
  adj(1,2) = 1
  adj(1,4) = 1
  adj(1,5) = 1
  adj(1,6) = 1
 
  adj(2,1) = 1
  adj(2,6) = 1
 
  adj(3,4) = 1
 
  adj(4,1) = 1
  adj(4,3) = 1
  adj(4,5) = 1
 
  adj(5,1) = 1
  adj(5,4) = 1
 
  adj(6,1) = 1
  adj(6,2) = 1
 
  call graph_adj_block ( adj, nnode, dad, order, stack, nblock )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of blocks = ', nblock

  call i4vec2_print ( nnode, dad, order, '  I, DAD(I), ORDER(I)' )

  call graph_adj_print ( adj, nnode, '  The graph:' )
 
  return
end
subroutine graph_adj_transitive_closure_test ( )

!*****************************************************************************80
!
!! graph_adj_transitive_closure_test() tests graph_adj_transitive_closure().
!
!    1--5      2
!    | /|
!    |/ |      8--3--7
!    4  6
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
  implicit none

  integer, parameter :: nnode = 8

  integer adj(nnode,nnode)
  integer c(nnode,nnode)
  integer i
  integer j

  do i = 1, nnode
    do j = 1, nnode
      if ( i == j ) then
        adj(i,j) = 1
      else
        adj(i,j) = 0
      end if
    end do
  end do

  adj(1,4) = 1
  adj(1,5) = 1

  adj(3,7) = 1
  adj(3,8) = 1

  adj(4,1) = 1
  adj(4,5) = 1

  adj(5,1) = 1
  adj(5,4) = 1
  adj(5,6) = 1

  adj(6,5) = 1

  adj(7,3) = 1

  adj(8,3) = 1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_adj_transitive_closure_test():'
  write ( *, '(a)' ) '  graph_adj_transitive_closure() finds the transitive closure '
  write ( *, '(a)' ) '  of a graph;'

  call graph_adj_print ( adj, nnode, '  The adjacency matrix for G:' )

  call graph_adj_transitive_closure ( adj, nnode, c )

  call graph_adj_print ( c, nnode, &
    '  Adjacency matrix for the transitive closure of G:' )

  return
end
subroutine graph_adj_color_next_test ( )

!*****************************************************************************80
!
!! graph_adj_color_next_test() tests graph_adj_color_next().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: nnode = 4
  integer, parameter :: maxstack = 20

  integer adj(nnode,nnode)
  integer color(nnode)
  integer i
  integer j
  logical more
  integer ncan(nnode)
  integer :: ncolor = 3
  integer stack(maxstack)

  data ( ( adj(i,j), j = 1, nnode ), i = 1, nnode) / &
    0, 1, 0, 1, &
    1, 0, 1, 0, &
    0, 1, 0, 1, &
    1, 0, 1, 0 /

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_adj_color_next_test():'
  write ( *, '(a)' ) '  graph_adj_color_next() produces colorings of a graph'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of colors available is ', ncolor

  call graph_adj_print ( adj, nnode, '  The graph:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Possible node colorings:'
  write ( *, '(a)' ) ' '

  more = .false.

  do

    call graph_adj_color_next ( adj, nnode, ncolor, color, stack, &
      maxstack, ncan, more )

    if ( .not. more ) then
      exit
    end if

    write ( *, '(19i4)' ) color(1:nnode)

  end do

  return
end
subroutine graph_adj_complement_test ( )

!*****************************************************************************80
!
!! graph_adj_complement_test() tests graph_adj_complement().
!
!    1--5      2
!    | /|
!    |/ |      8--3--7
!    4  6
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
  implicit none

  integer, parameter :: nnode = 8

  integer adj(nnode,nnode)
  integer c(nnode,nnode)
  integer i
  integer j

  do i = 1, nnode
    do j = 1, nnode
      if ( i == j ) then
        adj(i,j) = 1
      else
        adj(i,j) = 0
      end if
    end do
  end do

  adj(1,4) = 1
  adj(1,5) = 1

  adj(3,7) = 1
  adj(3,8) = 1

  adj(4,1) = 1
  adj(4,5) = 1

  adj(5,1) = 1
  adj(5,4) = 1
  adj(5,6) = 1

  adj(6,5) = 1

  adj(7,3) = 1

  adj(8,3) = 1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_adj_complement_test():'
  write ( *, '(a)' ) '  graph_adj_complement() finds the complement '
  write ( *, '(a)' ) '  of a graph;'

  call graph_adj_print ( adj, nnode, '  The adjacency matrix for G:' )

  call graph_adj_complement ( adj, nnode, c )

  call graph_adj_print ( c, nnode, &
    '  Adjacency matrix for the complement of G:' )

  return
end
subroutine graph_adj_connect_random_test ( )

!*****************************************************************************80
!
!! graph_adj_connect_random_test() tests graph_adj_connect_random().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: nedge = 8
  integer, parameter :: nnode = 6

  integer adj(nnode,nnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_adj_connect_random_test():'
  write ( *, '(a)' ) '  graph_adj_connect_random() returns a random connected graph;'
  write ( *, '(a)' ) ' '

  write ( *, '(a,i8)' ) '  Number of nodes is ', nnode
  write ( *, '(a,i8)' ) '  Number of edges is ', nedge

  call graph_adj_connect_random ( nnode, nedge, adj )

  call graph_adj_print ( adj, nnode, '  The graph:' )

  return
end 
subroutine graph_adj_is_edge_connected_test ( )

!*****************************************************************************80
!
!! graph_adj_is_edge_connected_test() tests graph_adj_is_edge_connected().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: nedge = 8
  integer, parameter :: nnode = 6

  integer adj(nnode,nnode)
  integer result

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_adj_is_edge_connected_test():'
  write ( *, '(a)' ) '  graph_adj_is_edge_connected() reports if a'
  write ( *, '(a)' ) '  graph is edgewise connected;'
  write ( *, '(a)' ) ' '

  write ( *, '(a,i8)' ) '  Number of nodes is ', nnode
  write ( *, '(a,i8)' ) '  Number of edges is ', nedge

  call graph_adj_connect_random ( nnode, nedge, adj )

  call graph_adj_print ( adj, nnode, '  The graph:' )
!
!  Check connectedness.
!
  call graph_adj_is_edge_connected ( adj, nnode, result )

  write ( *, '(a)' ) ' '
  if ( result == 0 ) then
    write ( *, '(a)' ) '  The graph is NOT edgewise connected.'
  else
    write ( *, '(a)' ) '  The graph IS edgewise connected.'
  end if

  return
end 
subroutine graph_adj_is_node_connected_test ( )

!*****************************************************************************80
!
!! graph_adj_is_node_connected_test() tests graph_adj_is_node_connected().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: nedge = 8
  integer, parameter :: nnode = 6

  integer adj(nnode,nnode)
  integer result

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_adj_is_node_connected_test():'
  write ( *, '(a)' ) '  graph_adj_is_node_connected() reports if a'
  write ( *, '(a)' ) '  graph is node connected;'
  write ( *, '(a)' ) ' '

  write ( *, '(a,i8)' ) '  Number of nodes is ', nnode
  write ( *, '(a,i8)' ) '  Number of edges is ', nedge

  call graph_adj_connect_random ( nnode, nedge, adj )

  call graph_adj_print ( adj, nnode, '  The graph:' )

  call graph_adj_is_node_connected ( adj, nnode, result )

  if ( result == 0 ) then
    write ( *, '(a)' ) '  The graph is NOT nodewise connected.'
  else
    write ( *, '(a)' ) '  The graph IS nodewise connected.'
  end if

  return
end 
subroutine graph_adj_cycle_test ( )

!*****************************************************************************80
!
!! graph_adj_cycle_test() tests graph_adj_cycle().
!
!  Discussion:
!
!    5--2--10--1--3--6
!           |  |  | /
!           8  |  9
!           |  |  
!           4--7  
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: maxstack = 100
  integer, parameter :: nnode = 10

  integer adj(nnode,nnode)
  integer dad(nnode)
  integer i
  integer order(nnode)
  integer stack(maxstack)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_adj_cycle_test():'
  write ( *, '(a)' ) '  graph_adj_cycle() searches for cycles in a graph.'

  adj(1:nnode,1:nnode) = 0

  adj(1,3) = 1
  adj(1,7) = 1
  adj(1,10) = 1

  adj(2,5) = 1
  adj(2,10) = 1

  adj(3,1) = 1
  adj(3,6) = 1
  adj(3,9) = 1

  adj(4,7) = 1
  adj(4,8) = 1

  adj(5,2) = 1

  adj(6,3) = 1
  adj(6,9) = 1

  adj(7,1) = 1
  adj(7,4) = 1
 
  adj(8,4) = 1
  adj(8,10) = 1

  adj(9,3) = 1
  adj(9,6) = 1

  adj(10,1) = 1
  adj(10,2) = 1
  adj(10,8) = 1
   
  call graph_adj_print ( adj, nnode, '  The graph:' )

  call graph_adj_cycle ( adj, nnode, dad, order, maxstack, stack )

  call i4vec2_print ( nnode, dad, order, '  Node, Dad, Order' )

  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '  Adjacency matrix with cycles marked.'
  write ( *, '(a)' ) ' '

  do i = 1, nnode
    write ( *, '(10i3)') adj(i,1:nnode)
  end do

  return
end 
subroutine graph_adj_degree_test ( )

!*****************************************************************************80
!
!! graph_adj_degree_test() tests graph_adj_degree().
!
!
!  5--2--10--1--3--6
!         |  |  | /
!         8  |  9
!         |  |  
!         4--7  
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: nnode = 10

  integer adj(nnode,nnode)
  integer degree(nnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_adj_degree_test():'
  write ( *, '(a)' ) '  graph_adj_degree() computes the degree of the nodes;'

  adj(1:nnode,1:nnode) = 0

  adj(1,3) = 1
  adj(1,7) = 1
  adj(1,10) = 1

  adj(2,5) = 1
  adj(2,10) = 1

  adj(3,1) = 1
  adj(3,6) = 1
  adj(3,9) = 1

  adj(4,7) = 1
  adj(4,8) = 1

  adj(5,2) = 1

  adj(6,3) = 1
  adj(6,9) = 1

  adj(7,1) = 1
  adj(7,4) = 1
 
  adj(8,4) = 1
  adj(8,10) = 1

  adj(9,3) = 1
  adj(9,6) = 1

  adj(10,1) = 1
  adj(10,2) = 1
  adj(10,8) = 1
   
  call graph_adj_print ( adj, nnode, '  The graph:' )

  call graph_adj_degree ( adj, nnode, degree )

  call i4vec_print ( nnode, degree, '  Node degrees:' )

  return
end
subroutine graph_adj_degree_max_test ( )

!*****************************************************************************80
!
!! graph_adj_degree_max_test() tests graph_adj_degree_max().
!
!
!  5--2--10--1--3--6
!         |  |  | /
!         8  |  9
!         |  |  
!         4--7  
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: nnode = 10

  integer adj(nnode,nnode)
  integer degree_max

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_adj_degree_max_test():'
  write ( *, '(a)' ) '  graph_adj_degree_max() computes the maximum'
  write ( *, '(a)' ) '  degree of the nodes;'

  adj(1:nnode,1:nnode) = 0

  adj(1,3) = 1
  adj(1,7) = 1
  adj(1,10) = 1

  adj(2,5) = 1
  adj(2,10) = 1

  adj(3,1) = 1
  adj(3,6) = 1
  adj(3,9) = 1

  adj(4,7) = 1
  adj(4,8) = 1

  adj(5,2) = 1

  adj(6,3) = 1
  adj(6,9) = 1

  adj(7,1) = 1
  adj(7,4) = 1
 
  adj(8,4) = 1
  adj(8,10) = 1

  adj(9,3) = 1
  adj(9,6) = 1

  adj(10,1) = 1
  adj(10,2) = 1
  adj(10,8) = 1
   
  call graph_adj_print ( adj, nnode, '  The graph:' )

  call graph_adj_degree_max ( adj, nnode, degree_max )

  write ( *, '(a)' ) ' ' 
  write ( *, '(a,i8)' ) '  Maximum node degree is ', degree_max
  write ( *, '(a)' ) ' '

  return
end 
subroutine graph_adj_degree_sequence_test ( )

!*****************************************************************************80
!
!! graph_adj_degree_sequence_test() tests graph_adj_degree_sequence().
!
!
!  5--2--10--1--3--6
!         |  |  | /
!         8  |  9
!         |  |  
!         4--7  
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: nnode = 10

  integer adj(nnode,nnode)
  integer degree_seq(nnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_adj_degree_sequence_test():'
  write ( *, '(a)' ) '  graph_adj_degree_sequence() computes the degree sequence;'

  adj(1:nnode,1:nnode) = 0

  adj(1,3) = 1
  adj(1,7) = 1
  adj(1,10) = 1

  adj(2,5) = 1
  adj(2,10) = 1

  adj(3,1) = 1
  adj(3,6) = 1
  adj(3,9) = 1

  adj(4,7) = 1
  adj(4,8) = 1

  adj(5,2) = 1

  adj(6,3) = 1
  adj(6,9) = 1

  adj(7,1) = 1
  adj(7,4) = 1
 
  adj(8,4) = 1
  adj(8,10) = 1

  adj(9,3) = 1
  adj(9,6) = 1

  adj(10,1) = 1
  adj(10,2) = 1
  adj(10,8) = 1
   
  call graph_adj_print ( adj, nnode, '  The graph:' )

  call graph_adj_degree_sequence ( adj, nnode, degree_seq )

  call i4vec_print ( nnode, degree_seq, '  Degree sequence:' )

  return
end 
subroutine graph_adj_depth_first_test ( )

!*****************************************************************************80
!
!! graph_adj_depth_first_test() tests graph_adj_depth_first().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: nnode = 13

  integer adj(nnode,nnode)
  integer dad(nnode)
  integer order(nnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_adj_depth_first_test():'
  write ( *, '(a)' ) '  graph_adj_depth_first() does depth first search of graph.'

  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1
  adj(1,3) = 1
  adj(1,6) = 1
  adj(1,7) = 1

  adj(5,4) = 1
  adj(5,7) = 1

  adj(6,5) = 1

  adj(8,9) = 1

  adj(10,11) = 1
  adj(10,12) = 1
  adj(10,13) = 1

  adj(12,13) = 1

  call graph_adj_print ( adj, nnode, '  The graph:' )

  call graph_adj_depth_first ( adj, nnode, dad, order )

  call i4vec2_print ( nnode, dad, order, '  Node, Dad, Order' )

  return
end
subroutine graph_adj_edge_count_test ( )

!*****************************************************************************80
!
!! graph_adj_edge_count_test() tests graph_adj_edge_count().
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
  implicit none

  integer, parameter :: nnode = 7

  integer adj(nnode,nnode)
  integer nedge

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'graph_adj_edge_count_test():'
  write ( *, '(a)' ) '  graph_adj_edge_count() counts the edges in a graph.'

  call graph_adj_example_bush ( adj )

  call graph_adj_print ( adj, nnode, '  Adjacency matrix:' )

  call graph_adj_edge_count ( adj, nnode, nedge )

  write ( *, '(a)' ) ''
  write ( *, '(a,i3)' ) '  Number of edges is ', nedge

  return
end
subroutine graph_adj_edge_select_test ( )

!*****************************************************************************80
!
!! graph_adj_edge_select_test() tests graph_adj_edge_select().
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
  implicit none

  integer, parameter :: nnode = 7

  integer adj(nnode,nnode)
  integer ni
  integer nj

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'graph_adj_edge_select_test():'
  write ( *, '(a)' ) '  graph_adj_edge_select() selects an edge from'
  write ( *, '(a)' ) '  a graph defined by an adjacency matrix.'

  call graph_adj_example_bush ( adj )

  call graph_adj_print ( adj, nnode, '  Adjacency matrix for bush example' )

  call graph_adj_edge_select ( adj, nnode, ni, nj )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  An edge of this graph extends from'
  write ( *, '(a,i2,a,i2)' ) '  node ', ni, ' to node ', nj

  return
end
subroutine graph_adj_edges_random_test ( )

!*****************************************************************************80
!
!! graph_adj_edges_random_test() tests graph_adj_edges_random().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: nedge = 10
  integer, parameter :: nnode = 6

  integer adj(nnode,nnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_adj_edges_random_test():'
  write ( *, '(a)' ) '  graph_adj_edges_random() returns a random graph'
  write ( *, '(a)' ) '  with a specified number of edges.'
  write ( *, '(a)' ) ' '

  write ( *, '(a,i8)' ) '  Number of edges requested = ', nedge

  call graph_adj_edges_random ( nnode, nedge, adj )

  call graph_adj_print ( adj, nnode, '  The graph:' )

  return
end
subroutine graph_adj_eigen_test ( )

!*****************************************************************************80
!
!! graph_adj_eigen_test() tests graph_adj_eigen().
!
!
!  5--2--10--1--3--6
!         |  |  | /
!         8  |  9
!         |  |
!         4--7
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nnode = 10

  integer adj(nnode,nnode)
  real ( kind = rk ) eigen(nnode)
  integer neigen

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_adj_eigen_test():'
  write ( *, '(a)' ) '  graph_adj_eigen() computes the eigenvalues of a graph.'

  adj(1:nnode,1:nnode) = 0

  adj(1,3) = 1
  adj(1,7) = 1
  adj(1,10) = 1

  adj(2,5) = 1
  adj(2,10) = 1

  adj(3,1) = 1
  adj(3,6) = 1
  adj(3,9) = 1

  adj(4,7) = 1
  adj(4,8) = 1

  adj(5,2) = 1

  adj(6,3) = 1
  adj(6,9) = 1

  adj(7,1) = 1
  adj(7,4) = 1

  adj(8,4) = 1
  adj(8,10) = 1

  adj(9,3) = 1
  adj(9,6) = 1

  adj(10,1) = 1
  adj(10,2) = 1
  adj(10,8) = 1

  call graph_adj_print ( adj, nnode, '  The graph:' )

  call graph_adj_eigen ( adj, nnode, neigen, eigen )

  call r8vec_print ( neigen, eigen, '  The eigenvalues:' )

  if ( neigen < nnode ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Warning!  Not all eigenvalues were computed.'
  end if

  return
end
subroutine graph_adj_depth_first_2_test ( )

!*****************************************************************************80
!
!! graph_adj_depth_first_2_test() tests graph_adj_depth_first_2().
!
!  Discussion:
!
!    This example is from page 22 of
!
!    Alan Gibbons,
!    Algorithmic Graph Theory,
!    Cambridge University Press, 1985
!    ISBN 0-521-28881-9
!
!    The correct result is
!
!    Node  Idad  Iorder
!
!     1      0      1
!     2      1      2
!     3      1      6
!     4      3      7
!     5      2      3
!     6      2      4
!     7      3      8
!     8      2      5
!     9      0      9
!    10      9     10
!    11     10     11
!    12     10     12
!    13     10     13
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: nnode = 13

  integer adj(nnode,nnode)
  integer dad(nnode)
  integer order(nnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_adj_depth_first_2_test():'
  write ( *, '(a)' ) '  graph_adj_depth_first_2() sets up depth-first traversal'
  write ( *, '(a)' ) '  of a graph described by an adjacency matrix.'
  write ( *, '(a)' ) ' '
 
  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1
  adj(1,3) = 1
  adj(1,4) = 1
  adj(1,5) = 1
  adj(1,6) = 1
  adj(1,7) = 1
  adj(1,8) = 1
 
  adj(2,1) = 1
  adj(2,5) = 1
  adj(2,6) = 1
  adj(2,8) = 1
 
  adj(3,1) = 1
  adj(3,4) = 1
  adj(3,7) = 1
 
  adj(4,1) = 1
  adj(4,3) = 1
 
  adj(5,1) = 1
  adj(5,2) = 1
 
  adj(6,1) = 1
  adj(6,2) = 1
 
  adj(7,1) = 1
  adj(7,3) = 1
 
  adj(8,1) = 1
  adj(8,2) = 1
 
  adj(9,10) = 1
  adj(9,13) = 1
 
  adj(10,9) = 1
  adj(10,11) = 1
  adj(10,12) = 1
  adj(10,13) = 1
 
  adj(11,10) = 1
  adj(11,12) = 1
 
  adj(12,10) = 1
  adj(12,11) = 1
 
  adj(13,9) = 1
  adj(13,10) = 1

  call graph_adj_print ( adj, nnode, '  The graph:' )
 
  call graph_adj_depth_first_2 ( adj, nnode, dad, order )
 
  call i4vec2_print ( nnode, dad, order, '  I, DAD(I), ORDER(I)' )
 
  return
end
subroutine test036 ( )

!*****************************************************************************80
!
!! TEST036 tests GRAPH_ADJ_HAM_NEXT.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: nnode = 20
  integer, parameter :: maxstack = 100

  integer adj(nnode,nnode)
  integer circuit(nnode)
  integer i
  integer j
  logical more
  integer ncan(nnode)
  integer stack(maxstack)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST036'
  write ( *, '(a)' ) '  GRAPH_ADJ_HAM_NEXT produces Hamilton circuits;'
  write ( *, '(a)' ) ' '
 
  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1
  adj(1,8) = 1
  adj(1,20) = 1

  adj(2,1) = 1
  adj(2,3) = 1
  adj(2,15) = 1

  adj(3,2) = 1
  adj(3,7) = 1
  adj(3,4) = 1

  adj(4,3) = 1
  adj(4,5) = 1
  adj(4,14) = 1

  adj(5,4) = 1
  adj(5,6) = 1
  adj(5,12) = 1

  adj(6,10) = 1
  adj(6,7) = 1

  adj(7,3) = 1
  adj(7,6) = 1
  adj(7,8) = 1

  adj(8,1) = 1
  adj(8,7) = 1
  adj(8,9) = 1

  adj(9,8) = 1
  adj(9,10) = 1
  adj(9,19) = 1

  adj(10,6) = 1
  adj(10,9) = 1
  adj(10,11) = 1

  adj(11,10) = 1
  adj(11,12) = 1
  adj(11,18) = 1

  adj(12,5) = 1
  adj(12,11) = 1
  adj(12,13) = 1

  adj(13,12) = 1
  adj(13,14) = 1
  adj(13,17) = 1

  adj(14,4) = 1
  adj(14,13) = 1
  adj(14,15) = 1

  adj(15,2) = 1
  adj(15,14) = 1
  adj(15,16) = 1

  adj(16,15) = 1
  adj(16,17) = 1
  adj(16,20) = 1

  adj(17,13) = 1
  adj(17,16) = 1
  adj(17,18) = 1

  adj(18,11) = 1
  adj(18,17) = 1
  adj(18,19) = 1

  adj(19,9) = 1
  adj(19,18) = 1
  adj(19,20) = 1

  adj(20,1) = 1
  adj(20,16) = 1
  adj(20,19) = 1
 
  do i = 1, nnode-1
    do j = i+1, nnode
      if ( adj(i,j) == 1 ) then
        adj(j,i) = 1
      end if
    end do
  end do

  call graph_adj_print ( adj, nnode, '  The graph:' )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Circuits:'
  write ( *, '(a)' ) ' '

  i = 0

  more = .false.

  do

    call graph_adj_ham_next ( adj, nnode, circuit, stack, maxstack, &
      ncan, more )

    if ( .not. more ) then
      exit
    end if

    i = i + 1
    write ( *, '(i3,2x,20i3)' ) i, circuit(1:nnode)

  end do
 
  return
end
subroutine test0365 ( )

!*****************************************************************************80
!
!! TEST0365 tests GRAPH_ADJ_HAM_NEXT.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: nnode = 9
  integer, parameter :: maxstack = 100

  integer adj(nnode,nnode)
  integer circuit(nnode)
  integer i
  logical more
  integer ncan(nnode)
  integer stack(maxstack)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0365'
  write ( *, '(a)' ) '  GRAPH_ADJ_HAM_NEXT produces Hamilton circuits;'
  write ( *, '(a)' ) ' '

  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1
  adj(1,4) = 1
  adj(1,6) = 1
 
  adj(2,1) = 1
  adj(2,3) = 1
  adj(2,7) = 1
 
  adj(3,2) = 1
  adj(3,4) = 1
  adj(3,6) = 1
 
  adj(4,1) = 1
  adj(4,3) = 1
  adj(4,7) = 1
 
  adj(5,6) = 1
  adj(5,7) = 1
  adj(5,9) = 1
 
  adj(6,1) = 1
  adj(6,3) = 1
  adj(6,5) = 1
  adj(6,8) = 1
 
  adj(7,2) = 1
  adj(7,4) = 1
  adj(7,5) = 1
 
  adj(8,6) = 1
  adj(8,9) = 1
 
  adj(9,5) = 1
  adj(9,8) = 1

  call graph_adj_print ( adj, nnode, '  The graph:' )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Circuits:'
  write ( *, '(a)' ) ' '

  i = 0

  more = .false.

  do

    call graph_adj_ham_next ( adj, nnode, circuit, stack, maxstack, &
      ncan, more )

    if ( .not. more ) then
      exit
    end if

    i = i + 1
    write ( *, '(2x,i3,2x,20i3)' ) i, circuit(1:nnode)

  end do
 
  return
end
subroutine graph_adj_ham_next_brute_test ( )

!*****************************************************************************80
!
!! graph_adj_ham_next_brute_test() tests graph_adj_ham_next_brute().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: nnode = 9

  integer adj(nnode,nnode)
  integer circuit(nnode)
  integer i
  integer iset

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_adj_ham_next_brute_test():'
  write ( *, '(a)' ) '  graph_adj_ham_next_brute() seeks circuits'
  write ( *, '(a)' ) '  in a graph which visit every node.'
  write ( *, '(a)' ) '  A brute force algorithm is used.'

  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1
  adj(1,4) = 1
  adj(1,6) = 1
 
  adj(2,1) = 1
  adj(2,3) = 1
  adj(2,7) = 1
 
  adj(3,2) = 1
  adj(3,4) = 1
  adj(3,6) = 1
 
  adj(4,1) = 1
  adj(4,3) = 1
  adj(4,7) = 1
 
  adj(5,6) = 1
  adj(5,7) = 1
  adj(5,9) = 1
 
  adj(6,1) = 1
  adj(6,3) = 1
  adj(6,5) = 1
  adj(6,8) = 1
 
  adj(7,2) = 1
  adj(7,4) = 1
  adj(7,5) = 1
 
  adj(8,6) = 1
  adj(8,9) = 1
 
  adj(9,5) = 1
  adj(9,8) = 1

  call graph_adj_print ( adj, nnode, '  The graph:' )
 
  iset = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Circuits:'
  write ( *, '(a)' ) ' '

  i = 0
  
  do
 
    call graph_adj_ham_next_brute ( adj, nnode, circuit, iset )
 
    if ( iset == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  No more circuits were found.'
      exit
    end if

    i = i + 1
    write ( *, '(2x,i3,2x,20i3)' ) i, circuit(1:nnode)

  end do
 
  return
end
subroutine graph_adj_is_bipartite_test ( )

!*****************************************************************************80
!
!! graph_adj_is_bipartite_test() tests graph_adj_is_bipartite().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: nnode1 = 4
  integer, parameter :: nnode2 = 6
  integer, parameter :: nnode = nnode1 + nnode2

  integer adj(nnode,nnode)
  integer nedge
  integer result

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_adj_is_bipartite_test()'
  write ( *, '(a)' ) '  graph_adj_is_bipartite() reports if a graph is bipartite.'
  write ( *, '(a)' ) ' '

  call graph_adj_bipartite_random ( nnode1, nnode2, nedge, adj )

  call graph_adj_print ( adj, nnode, '  The graph:' )

  call graph_adj_is_bipartite ( adj, nnode, result )

  write ( *, '(a)' ) ' '
  if ( result == 0 ) then
    write ( *, '(a)' ) '  The graph is NOT bipartite.'
  else
    write ( *, '(a)' ) '  The graph IS bipartite.'
  end if

  return
end 
subroutine graph_adj_is_tree_test ( )

!*****************************************************************************80
!
!! graph_adj_is_tree_test() tests graph_adj_is_tree().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: nnode = 6
  integer, parameter :: nedge = nnode - 1

  integer adj(nnode,nnode)
  integer result

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_adj_is_tree_test():'
  write ( *, '(a)' ) '  GRAPH_ADJ_IS_TREE reports if a graph is a tree.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nodes is ', nnode
  write ( *, '(a,i8)' ) '  Number of edges is ', nedge

  call graph_adj_connect_random ( nnode, nedge, adj )

  call graph_adj_print ( adj, nnode, '  The graph:' )

  call graph_adj_is_tree ( adj, nnode, result )

  if ( result == 0 ) then
    write ( *, '(a)' ) '  The graph is NOT a tree.'
  else
    write ( *, '(a)' ) '  The graph IS a tree.'
  end if

  return
end 
subroutine graph_adj_random_test ( )

!*****************************************************************************80
!
!! graph_adj_random_test() tests graph_adj_random().
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nnode = 20
  integer, parameter :: test_num = 3

  integer adj(nnode,nnode)
  real ( kind = rk ) eigen(nnode)
  integer nedge
  integer nedge_max
  integer neigen
  real ( kind = rk ) prob
  real ( kind = rk ), dimension ( test_num ) :: prob_test = (/ &
    0.25D+00, 0.40D+00,  0.65D+00 /)
  real ( kind = rk ) ratio
  integer test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_adj_random_test():'
  write ( *, '(a)' ) '  graph_adj_random() returns a random graph, for which'
  write ( *, '(a)' ) '  edges are generated with a given probability.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here, we show the effect of increasing connectivity'
  write ( *, '(a)' ) '  on the singularity of the adjacency matrix.'

  do test = 1, test_num

    prob = prob_test(test)

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Probability of edge generation = ', prob

    call graph_adj_random ( nnode, prob, adj )

    nedge = sum ( adj ) / 2
    nedge_max = ( nnode * ( nnode - 1 ) ) / 2
    ratio = real ( nedge, kind = rk ) / real ( nedge_max, kind = rk )

    write ( *, '(a,i8)' ) '  Number of edges generated = ', nedge
    write ( *, '(a,i8)' ) '  Maximum number of edges = ', nedge_max
    write ( *, '(a,g14.6)' ) '  Generated / Maximum = ', ratio
      
    call graph_adj_print ( adj, nnode, '  The graph:' )

    call graph_adj_eigen ( adj, nnode, neigen, eigen )

    call r8vec_print ( neigen, eigen, '  The eigenvalues:' )

  end do

  return
end
subroutine graph_adj_reduce_test ( )

!*****************************************************************************80
!
!! graph_adj_reduce_test() tests graph_adj_reduce()
!
!    1--5      2
!    | /|
!    |/ |      8--3--7
!    4  6
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: nnode = 8

  integer adj(nnode,nnode)
  integer i
  integer j

  do i = 1, nnode
    do j = 1, nnode
      if ( i == j ) then
        adj(i,j) = 1
      else
        adj(i,j) = 0
      end if
    end do
  end do

  adj(1,4) = 1
  adj(1,5) = 1

  adj(3,7) = 1
  adj(3,8) = 1

  adj(4,1) = 1
  adj(4,5) = 1

  adj(5,1) = 1
  adj(5,4) = 1
  adj(5,6) = 1

  adj(6,5) = 1

  adj(7,3) = 1

  adj(8,3) = 1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_adj_reduce_test():'
  write ( *, '(a)' ) '  graph_adj_reduce() finds the transitive reduction'
  write ( *, '(a)' ) '  of a graph.'

  call graph_adj_print ( adj, nnode, '  The adjacency matrix for G:' )

  call graph_adj_reduce ( adj, nnode )

  call graph_adj_print ( adj, nnode, &
    '  Adjacency matrix for the transitive reduction of G:' )

  return
end
subroutine graph_adj_span_tree_test ( )

!*****************************************************************************80
!
!! graph_adj_span_tree_test() tests graph_adj_span_tree().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: nnode = 13

  integer adj(nnode,nnode)
  integer inode(nnode-1)
  integer jnode(nnode-1)

  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1
  adj(1,3) = 1
  adj(1,4) = 1
  adj(1,5) = 1
  adj(1,6) = 1
  adj(1,7) = 1
  adj(1,8) = 1
 
  adj(2,1) = 1
  adj(2,5) = 1
  adj(2,6) = 1
  adj(2,8) = 1
 
  adj(3,1) = 1
  adj(3,4) = 1
  adj(3,7) = 1
 
  adj(4,1) = 1
  adj(4,3) = 1
 
  adj(5,1) = 1
  adj(5,2) = 1
 
  adj(6,1) = 1
  adj(6,2) = 1
 
  adj(7,1) = 1
  adj(7,3) = 1
 
  adj(8,1) = 1
  adj(8,2) = 1
  adj(8,9) = 1

  adj(9,8) = 1 
  adj(9,10) = 1
  adj(9,13) = 1
 
  adj(10,9) = 1
  adj(10,11) = 1
  adj(10,12) = 1
  adj(10,13) = 1
 
  adj(11,10) = 1
  adj(11,12) = 1
 
  adj(12,10) = 1
  adj(12,11) = 1
 
  adj(13,9) = 1
  adj(13,10) = 1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_adj_span_tree_test():'
  write ( *, '(a)' ) '  graph_adj_span_tree() constructs a spanning tree of a graph.'

  call graph_adj_print ( adj, nnode, '  The graph:' )

  call graph_adj_span_tree ( adj, nnode, inode, jnode )

  call graph_arc_print ( nnode-1, inode, jnode, '  The spanning tree:' )

  return
end
subroutine graph_adj_span_tree_enum_test ( )

!*****************************************************************************80
!
!! graph_adj_span_tree_enum_test() tests graph_adj_span_tree_enum().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: nnode = 13

  integer adj(nnode,nnode)
  integer tree_num

  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1
  adj(1,3) = 1
  adj(1,4) = 1
  adj(1,5) = 1
  adj(1,6) = 1
  adj(1,7) = 1
  adj(1,8) = 1
 
  adj(2,1) = 1
  adj(2,5) = 1
  adj(2,6) = 1
  adj(2,8) = 1
 
  adj(3,1) = 1
  adj(3,4) = 1
  adj(3,7) = 1
 
  adj(4,1) = 1
  adj(4,3) = 1
 
  adj(5,1) = 1
  adj(5,2) = 1
 
  adj(6,1) = 1
  adj(6,2) = 1
 
  adj(7,1) = 1
  adj(7,3) = 1
 
  adj(8,1) = 1
  adj(8,2) = 1
  adj(8,9) = 1

  adj(9,8) = 1 
  adj(9,10) = 1
  adj(9,13) = 1
 
  adj(10,9) = 1
  adj(10,11) = 1
  adj(10,12) = 1
  adj(10,13) = 1
 
  adj(11,10) = 1
  adj(11,12) = 1
 
  adj(12,10) = 1
  adj(12,11) = 1
 
  adj(13,9) = 1
  adj(13,10) = 1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_adj_span_tree_enum_test():'
  write ( *, '(a)' ) '  graph_adj_span_tree_enum() enumerates the spanning trees'
  write ( *, '(a)' ) '  of a graph.'

  call graph_adj_print ( adj, nnode, '  The graph:' )

  call graph_adj_span_tree_enum ( adj, nnode, tree_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Total number of spanning trees is ', tree_num

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! timestamp() prints the current YMDHMS date as a time stamp.
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
!    15 August 2021
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

