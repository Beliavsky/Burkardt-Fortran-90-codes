program main

!*****************************************************************************80
!
!! graph_arc_test() tests graph_arc().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 March 2023
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_arc_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  graph_arc() implements graph algorithms.'

  call graph_arc_chromatic_test ( )
  call graph_arc_complement_test ( )
  call graph_arc_degree_test ( )
  call graph_arc_edge_con2_test ( )
  call graph_arc_edge_sort_test ( )
  call graph_arc_euler_circ_test ()
  call graph_arc_euler_circ_next_test ( )
  call graph_arc_is_eulerian_test ( )
  call graph_arc_match_test ( )
  call graph_arc_min_path_test ( )
  call graph_arc_min_span_tree_test ( )
  call graph_arc_node_count_test ( )
  call graph_arc_span_forest_test ( )
  call graph_arc_span_tree_test ( )
  call graph_arc_to_digraph_arc_test ( )
  call graph_arc_to_graph_adj_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_theory_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine graph_arc_chromatic_test ( )

!*****************************************************************************80
!
!! graph_arc_chromatic_test tests graph_arc_chromatic().
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
!    John Burkardt
!
  implicit none

  integer, parameter :: nnode = 6 
  integer, parameter :: nedge = 12
  integer, parameter :: maxstack = nnode * nedge

  integer i
  integer iarray(nnode)
  integer inode(nedge)
  integer jarray(nnode)
  integer jnode(nedge)
  integer karray(nnode)
  integer stack(2,maxstack)

  data inode / &
    1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5 /

  data jnode / &
    2, 3, 4, 5, 3, 4, 6, 5, 6, 5, 6, 6 /

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_arc_chromatic_test():'
  write ( *, '(a)' ) '  graph_arc_chromatic() finds the chromatic polynomial'
  write ( *, '(a)' ) '  of a graph.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The end point arrays:'
  write ( *, '(a)' ) ' '
  write ( *, '(19i4)' ) ( inode(i), i = 1, nedge )
  write ( *, '(19i4)' ) ( jnode(i), i = 1, nedge )
 
  call graph_arc_chromatic ( nnode, nedge, inode, jnode, iarray, jarray, &
    karray, stack, maxstack )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The chromatic polynomial:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Power sum form:'
  write ( *, '(19i4)' ) iarray(1:nnode)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Tutte or tree form:'
  write ( *, '(19i4)' ) jarray(1:nnode)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Stirling form:'
  write ( *, '(19i4)' ) karray(1:nnode)
 
  return
end
subroutine graph_arc_complement_test ( )

!*****************************************************************************80
!
!! graph_arc_complement_test() tests graph_arc_complement().
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

  integer, parameter :: maxedge = 90
  integer, parameter :: maxnode = 10

  integer inode(maxedge)
  integer inode2(maxedge)
  integer jnode(maxedge)
  integer jnode2(maxedge)
  integer nedge
  integer nedge2
  integer nnode
  real ( kind = rk ) x(maxnode)
  real ( kind = rk ) y(maxnode)
  real ( kind = rk ) z(maxnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_arc_complement_test():'
  write ( *, '(a)' ) '  graph_arc_complement() computes the complement'
  write ( *, '(a)' ) '  of a graph described by its edge array;'

  call graph_arc_example_diamond ( inode, jnode, maxedge, nedge, nnode, x, y, z )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of edges in original graph is ', nedge
  write ( *, '(a,i8)' ) '  Number of nodes is ', nnode
 
  call graph_arc_edge_sort ( nedge, inode, jnode )
 
  call graph_arc_print ( nedge, inode, jnode, '  The graph:' )
 
  call graph_arc_complement ( inode, jnode, inode2, jnode2, maxedge, nedge, &
    nedge2, nnode )
 
  write ( *, '(a,i8)' ) 'Number of edges in complement is ', nedge2

  call graph_arc_edge_sort ( nedge2, inode2, jnode2 )
 
  call graph_arc_print ( nedge, inode, jnode, '  The complement graph:' )
 
  return
end
subroutine graph_arc_degree_test ( )

!*****************************************************************************80
!
!! graph_arc_degree_test() tests graph_arc_degree().
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

  integer, parameter :: nedge = 11
  integer, parameter :: nnode = 10

  integer degree(nnode)
  integer inode(nedge)
  integer jnode(nedge)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_arc_degree_test():'
  write ( *, '(a)' ) '  graph_arc_degree() computes the degree of the nodes;'

  inode = (/ 1, 1,  1, 2,  2, 3, 3, 4, 4, 6,  8 /)
  jnode = (/ 3, 7, 10, 5, 10, 6, 9, 7, 8, 9, 10 /)

  call graph_arc_print ( nedge, inode, jnode, '  The graph:' )

  call graph_arc_degree ( nnode, nedge, inode, jnode, degree )

  call i4vec_print ( nnode, degree, '  The node degrees:' )

  return
end 
subroutine graph_arc_edge_con2_test ( )

!*****************************************************************************80
!
!! graph_arc_edge_con2_test() tests graph_arc_edge_con2().
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

  integer, parameter :: nedge = 17
  integer, parameter :: nnode = 9

  integer edge_con
  integer, dimension ( nedge ) :: inode = &
    (/ 6,2,3,6,7,1,4,7,3,4,9,6,5,4,2,9,4 /)
  integer, dimension ( nedge ) :: jnode = &
    (/ 8,5,1,3,2,8,3,5,8,1,2,1,9,8,6,7,2 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_arc_edge_con2_test():'
  write ( *, '(a)' ) '  graph_arc_edge_con2() finds graph edge connectivity.'

  call graph_arc_print ( nedge, inode, jnode, '  The arc list of the graph:' )

  call graph_arc_edge_con2 ( nnode, nedge, inode, jnode, edge_con )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The computed edge connectivity is ', edge_con

  return
end
subroutine graph_arc_edge_sort_test ( )

!*****************************************************************************80
!
!! graph_arc_edge_sort_test() tests graph_arc_edge_sort().
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

  integer, parameter :: maxedge = 90
  integer, parameter :: maxnode = 10

  integer inode(maxedge)
  integer jnode(maxedge)
  integer nedge
  integer nnode
  real ( kind = rk ) x(maxnode)
  real ( kind = rk ) y(maxnode)
  real ( kind = rk ) z(maxnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_arc_edge_sort_test():'
  write ( *, '(a)' ) '  graph_arc_edge_sort() sorts the edge array of a graph.'

  call graph_arc_example_diamond ( inode, jnode, maxedge, nedge, nnode, x, y, z )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of edges in original graph is ', nedge
  write ( *, '(a,i8)' ) '  Number of nodes is ', nnode
 
  call graph_arc_edge_sort ( nedge, inode, jnode )
 
  call graph_arc_print ( nedge, inode, jnode, '  The graph:' )

  return
end
subroutine graph_arc_euler_circ_test ( )

!*****************************************************************************80
!
!! graph_arc_euler_circ_test() tests graph_arc_euler_circ().
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
  integer, parameter :: nnode = 5

  integer circuit(nedge)
  integer, dimension ( nedge ) :: inode = (/ 1, 1, 1, 1, 2, 2, 2, 3, 3, 4 /)
  integer, dimension ( nedge ) :: jnode = (/ 2, 3, 4, 5, 3, 4, 5, 4, 5, 5 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_arc_euler_circ_test():'
  write ( *, '(a)' ) '  graph_arc_euler_circ() returns an Euler circuit'
  write ( *, '(a)' ) '  of a graph.'

  call graph_arc_print ( nedge, inode, jnode, '  The graph:' )

  call graph_arc_euler_circ ( nnode, nedge, inode, jnode, circuit )

  call i4vec_print ( nedge, circuit, '  The nodes in the Euler circuit:' )

  return
end
subroutine graph_arc_euler_circ_next_test ( )

!*****************************************************************************80
!
!! graph_arc_euler_circ_next_test() tests graph_arc_euler_circ_next().
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

  integer, parameter :: maxstack = 130
  integer, parameter :: nedge = 10
  integer, parameter :: nnode = 5

  integer circuit(nedge)
  integer i
  integer, dimension ( nedge ) :: inode = (/ 1, 1, 1, 1, 2, 2, 2, 3, 3, 4 /)
  integer, dimension ( nedge ) :: jnode = (/ 2, 3, 4, 5, 3, 4, 5, 4, 5, 5 /)
  logical more
  integer ncan(nedge)
  integer stack(maxstack)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_arc_euler_circ_next_test():'
  write ( *, '(a)' ) '  graph_arc_euler_circ_next() finds the next'
  write ( *, '(a)' ) '  Euler circuit of a graph.'
  write ( *, '(a)' ) ' '

  call graph_arc_print ( nedge, inode, jnode, '  The graph:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Circuits:'
  write ( *, '(a)' ) ' '
  i = 0
  more = .false.

  do

    call graph_arc_euler_circ_next ( nedge, inode, jnode, circuit, stack, &
      maxstack, ncan, more )

    if ( .not. more ) then
      exit
    end if

    i = i + 1
    write ( *, '(i3,2x,20i3)' ) i, circuit(1:nedge)

  end do

  return
end
subroutine graph_arc_is_eulerian_test ( )

!*****************************************************************************80
!
!! graph_arc_is_eulerian_test() tests graph_arc_is_eulerian().
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

  integer, parameter :: maxstack = 130
  integer, parameter :: nedge = 10
  integer, parameter :: nnode = 5

  integer circuit(nedge)
  integer degree(nnode)
  integer i
  integer, dimension ( nedge ) :: inode = (/ 1, 1, 1, 1, 2, 2, 2, 3, 3, 4 /)
  integer, dimension ( nedge ) :: jnode = (/ 2, 3, 4, 5, 3, 4, 5, 4, 5, 5 /)
  logical more
  integer ncan(nedge)
  integer result
  integer stack(maxstack)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_arc_is_eulerian_test()'
  write ( *, '(a)' ) '  graph_arc_is_eulerian() checks if a graph has an'
  write ( *, '(a)' ) '  Euler circuit.'

  call graph_arc_print ( nedge, inode, jnode, '  The graph:' )

  call graph_arc_is_eulerian ( nnode, nedge, inode, jnode, degree, result )

  if ( result == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The graph is NOT eulerian.'
    return
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The graph is eulerian.'
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Circuits:'
  write ( *, '(a)' ) ' '
  i = 0
  more = .false.

  do

    call graph_arc_euler_circ_next ( nedge, inode, jnode, circuit, stack, &
      maxstack, ncan, more )

    if ( .not. more ) then
      exit
    end if

    i = i + 1
    write ( *, '(i3,2x,20i3)' ) i, circuit(1:nedge)

  end do

  return
end

subroutine graph_arc_match_test()

!*****************************************************************************80
!
!  graph_arc_match_test() tests graph_arc_match().
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
  integer, parameter :: nedge = 14
  integer, parameter :: nnode = 12

  integer, dimension ( nedge ) :: inode = &
    (/ 6, 9, 3,  4, 11, 6, 4, 5,  6, 10, 3, 4, 1, 3 /)
  integer, dimension ( nedge ) :: jnode = &
    (/ 2, 7, 7, 10,  5, 8, 6, 7, 12,  2, 1, 2, 5, 5 /)
  integer, dimension ( nnode ) :: match
  integer, dimension ( nnode ) :: type = (/ &
    1, 1, 2, 1, 2, 2, 1, 2, 2, 2, 1, 1 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_arc_match_test():'
  write ( *, '(a)' ) '  graph_arc_match() finds a maximal matching in a graph.'

  call graph_arc_print ( nedge, inode, jnode, '  The edge list of the graph:' )

  call i4vec_print ( nnode, type, '  Nodes and their types:' )

  call graph_arc_match ( nnode, nedge, inode, jnode, type, match )

  call i4vec_print ( nnode, match, '  Node and matching node:' )

  return
end
subroutine graph_arc_min_path_test ( )

!*****************************************************************************80
!
!! graph_arc_min_path_test() tests graph_arc_min_path().
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

  integer, parameter :: nnode = 5
  integer, parameter :: nedge = 6

  real ( kind = rk ), save, dimension ( nedge ) :: cost = (/ &
    1.0D+00, 1.0D+00, 3.0D+00, 2.0D+00, 2.0D+00, 5.0D+00 /)
  real ( kind = rk ) dist(nnode,nnode)
  integer, save, dimension ( nedge ) :: inode = (/ 1, 1, 2, 2, 3, 3 /)
  integer istart
  integer istop
  integer, save, dimension ( nedge ) :: jnode = (/ 2, 3, 3, 5, 4, 5 /)
  integer num_path
  integer path(nnode)
  real ( kind = rk ) path_length

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_arc_min_path_test():'
  write ( *, '(a)' ) '  graph_arc_min_path() computes the shortest path from one'
  write ( *, '(a)' ) '  node to another.'
  write ( *, '(a)' ) ' '

  call graph_arc_weight_print ( nedge, inode, jnode, cost, &
    '  The weighted graph:' )

  dist(1:nnode,1:nnode) = 0.0D+00

  do istart = 1, nnode
    do istop = istart+1, nnode
      call graph_arc_min_path ( nnode, nedge, inode, jnode, cost, istart, &
        istop, num_path, path, path_length )
      dist(istart,istop) = path_length
      dist(istop,istart) = path_length
    end do
  end do

  call graph_dist_print ( dist, nnode, &
    '  The distance matrix constructed by GRAPH_ARC_MIN_PATH:' )
 
  istart = 4
  istop = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The routine actually also computes the path.'
  write ( *, '(a,i8)' ) '  For instance, starting at node ', istart
  write ( *, '(a,i8)' ) '  we compute the shortest path to node ', istop

  call graph_arc_min_path ( nnode, nedge, inode, jnode, cost, istart, &
    istop, num_path, path, path_length )

  call i4vec_print ( num_path, path, '  The path:' )

  return
end
subroutine graph_arc_min_span_tree_test ( )

!*****************************************************************************80
!
!! graph_arc_min_span_tree_test() tests graph_arc_min_span_tree().
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

  integer, parameter :: nedge = 10
  integer, parameter :: nnode = 5

  real ( kind = rk ), dimension ( nedge ) :: cost = &
    (/ 100.0, 125.0, 120.0, 110.0, 40.0, 65.0, 60.0, 45.0, 55.0, 50.0 /)
  real ( kind = rk ), dimension ( nnode-1) :: ctree
  integer, dimension ( nedge ) :: inode = (/ 1, 1, 1, 1, 2, 2, 2, 3, 3, 4 /)
  integer i
  integer itree(nnode-1)
  integer j
  integer, dimension ( nedge ) :: jnode = (/ 2, 3, 4, 5, 3, 4, 5, 4, 5, 5 /)
  integer jtree(nnode-1)
  real ( kind = rk ) tree_cost

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_arc_min_span_tree_test():'
  write ( *, '(a)' ) '  graph_arc_min_span_tree() finds a minimum length'
  write ( *, '(a)' ) '  spanning tree.'
  write ( *, '(a)' ) ' '

  call graph_arc_weight_print ( nedge, inode, jnode, cost, &
    '  The weighted graph:' )

  call graph_arc_min_span_tree ( nnode, nedge, inode, jnode, cost, &
    itree, jtree, tree_cost )

  do i = 1, nnode - 1
    ctree(i) = 0.0D+00
    do j = 1, nedge
      if ( ( inode(j) == itree(i) .and. jnode(j) == jtree(i) ) .or. &
           ( inode(j) == jtree(i) .and. jnode(j) == itree(i) ) ) then
        ctree(i) = cost(j)
        exit
      end if
    end do
  end do

  call graph_arc_weight_print ( nnode-1, itree, jtree, ctree, &
    '  The minimal spanning tree:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The length of the minimal tree is ', sum ( ctree )

  return
end
subroutine graph_arc_node_count_test ( )

!*****************************************************************************80
!
!! graph_arc_node_count_test() tests graph_arc_node_count().
!
!
!  5--2--100-1--3--0
!         |  |  | /
!        88  |  9
!         |  |  
!      (-4)--7  
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

  integer, parameter :: nedge = 11

  integer, dimension ( nedge ) :: inode = &
    (/ 1, 1,   1, 2,   2, 3, 3, -4, -4, 0,  88 /)
  integer, dimension ( nedge ) :: jnode = &
    (/ 3, 7, 100, 5, 100, 0, 9,  7, 88, 9, 100 /)
  integer mnode
  integer nnode

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_arc_node_count_test():'
  write ( *, '(a)' ) '  graph_arc_node_count() counts the nodes and'
  write ( *, '(a)' ) '  finds the highest label.'

  call graph_arc_print ( nedge, inode, jnode, '  The graph:' )

  call graph_arc_node_count ( nedge, inode, jnode, mnode, nnode )

  write ( *, '(a,i8)' ) '  Number of nodes is    ', nnode
  write ( *, '(a,i8)' ) '  Maximum node label is ', mnode

  return
end
subroutine graph_arc_span_forest_test ( )

!*****************************************************************************80
!
!! graph_arc_span_forest_test() tests graph_arc_span_forest().
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

  integer, parameter :: nnode = 14
  integer, parameter :: nedge = 10

  integer component(nnode)
  integer, save, dimension ( nedge ) :: inode = &
    (/  2,  4,  1,  7,  5,  2,  6,  2,  3,  4 /)
  integer, save, dimension ( nedge ) :: jnode = &
    (/  3,  7,  9, 11,  8,  5, 10,  8,  8, 11 /)
  integer ncomp

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_arc_span_forest_test():'
  write ( *, '(a)' ) '  graph_arc_span_forest()'
  write ( *, '(a)' ) '  computes a spanning forest for a graph'
 
  call graph_arc_print ( nedge, inode, jnode, '  The graph:' )

  call graph_arc_span_forest ( nnode, nedge, inode, jnode, ncomp, component )

  call graph_arc_print ( nedge, inode, jnode, &
    '  The reordered endpoint array:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of connected components = ', ncomp

  call i4vec_print ( nnode, component, '  Node component membership:' )

  return
end
subroutine graph_arc_span_tree_test ( )

!*****************************************************************************80
!
!! graph_arc_span_tree_test() tests graph_arc_span_tree().
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

  integer, parameter :: nedge = 18
  integer, parameter :: nnode = 13

  integer dad(nnode)
  integer inode(nedge)
  integer jnode(nedge)

  inode(1) = 1
  jnode(1) = 2
  inode(2) = 1
  jnode(2) = 3
  inode(3) = 1
  jnode(3) = 4
  inode(4) = 1
  jnode(4) = 5
  inode(5) = 1
  jnode(5) = 6
  inode(6) = 1
  jnode(6) = 7
  inode(7) = 1
  jnode(7) = 8
 
  inode(8) = 2
  jnode(8) = 5
  inode(9) = 2
  jnode(9) = 6
  inode(10) = 2
  jnode(10) = 8
 
  inode(11) = 3
  jnode(11) = 4
  inode(12) = 3
  jnode(12) = 7
 
  inode(13) = 9
  jnode(13) = 10
  inode(14) = 9
  jnode(14) = 13
 
  inode(15) = 10
  jnode(15) = 11
  inode(16) = 10
  jnode(16) = 12
  inode(17) = 10
  jnode(17) = 13
 
  inode(18) = 11
  jnode(18) = 12

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_arc_span_tree_test():'
  write ( *, '(a)' ) '  graph_arc_span_tree() constructs a spanning tree.'
  write ( *, '(a)' ) ' '

  call graph_arc_print ( nedge, inode, jnode, '  The graph:' )

  call graph_arc_span_tree ( nedge, inode, jnode, nnode, dad )

  call i4vec_print ( nnode, dad, '  Nodes and Parent Nodes:' )

  return
end
subroutine graph_arc_to_digraph_arc_test ( )

!*****************************************************************************80
!
!! graph_arc_to_digraph_arc_test() tests graph_arc_to_digraph_arc().
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
  integer, parameter :: maxarc = 2 * nedge

  integer iarc(maxarc)
  integer, dimension ( nedge ) :: inode = (/ 1, 1, 1, 2, 3, 4, 2, 4 /)
  integer jarc(maxarc)
  integer, dimension ( nedge ) :: jnode = (/ 2, 1, 4, 1, 2, 1, 3, 2 /)
  integer narc

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_arc_to_digraph_arc_test():'
  write ( *, '(a)' ) '  graph_arc_to_digraph_arc() makes a directed graph'
  write ( *, '(a)' ) '  from an undirected one.'

  call graph_arc_print ( nedge, inode, jnode, '  The graph:' )

  call graph_arc_to_digraph_arc ( iarc, jarc, inode, jnode, maxarc, narc, &
    nedge )

  call digraph_arc_print ( narc, iarc, jarc, '  The digraph:' )

  return
end
subroutine graph_arc_to_graph_adj_test ( )

!*****************************************************************************80
!
!! graph_arc_to_graph_adj_test() tests graph_arc_to_graph_adj().
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
  integer, parameter :: nnode = 5

  integer adj(nnode,nnode)
  integer, dimension ( nedge ) :: inode = (/ 1, 1, 1, 2, 3, 4, 2, 4 /)
  integer, dimension ( nedge ) :: jnode = (/ 2, 1, 4, 1, 2, 1, 3, 2 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_arc_to_graph_adj_test():'
  write ( *, '(a)' ) '  graph_arc_to_graph_adj() converts an arclist'
  write ( *, '(a)' ) '  graph to an adjacency graph.'
  write ( *, '(a)' ) ' '

  call graph_arc_print ( nedge, inode, jnode, '  The graph:' )

  call graph_arc_to_graph_adj ( nedge, inode, jnode, adj, nnode )

  call graph_adj_print ( adj, nnode, '  The graph:' )

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
