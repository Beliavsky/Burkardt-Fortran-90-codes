program main

!*****************************************************************************80
!
!! digraph_adjacency_test() tests digraph_adj().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 March 2023
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'digraph_adjacency_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test digraph_adj(), which implements digraph algorithms.'

  call digraph_adjacency_components_test ( )
  call digraph_adjacency_cycle_test ( )
  call digraph_adjacency_degree_test ( )
  call digraph_adjacency_degree_max_test ( )
  call digraph_adjacency_degree_sequence_test ( )
  call digraph_adjacency_eigen_test ( )
  call digraph_adjacency_ham_next_test ( )
  call digraph_adjacency_ham_next_brute_test ( )
  call digraph_adjacency_ham_path_next_brute_test ( )
  call digraph_adjacency_is_edge_connected_test ( )
  call digraph_adjacency_is_eulerian_test ( )
  call digraph_adjacency_is_strongly_connected_test ( )
  call digraph_adjacency_is_tournament_test ( )
  call digraph_adjacency_is_transitive_test ( )
  call digraph_adjacency_random_test ( )
  call digraph_adjacency_reduce_test ( )
  call digraph_adjacency_to_digraph_arc_test ( )
  call digraph_adjacency_to_incidence_test ( )
  call digraph_adjacency_top_sort_test ( )
  call digraph_adjacency_tournament_random_test ( )
  call digraph_adjacency_transitive_closure_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'digraph_adjacency_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine digraph_adjacency_components_test ( )

!*****************************************************************************80
!
!! digraph_adjacency_components_test() tests digraph_adjacency_components().
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
  integer comp(nnode)
  integer dad(nnode)
  integer i
  integer j
  integer ncomp
  integer order(nnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'digraph_adjacency_components_test():'
  write ( *, '(a)' ) '  digraph_adjacency_components() finds strongly connected'
  write ( *, '(a)' ) '  components of a directed graph.'

  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1
  adj(1,11) = 1

  adj(2,3) = 1
  adj(2,6) = 1

  adj(3,4) = 1
  adj(3,5) = 1

  adj(4,3) = 1

  adj(5,4) = 1

  adj(6,7) = 1
  adj(6,8) = 1

  adj(7,6) = 1

  adj(8,9) = 1
  adj(8,10) = 1

  adj(9,7) = 1

  adj(10,9) = 1

  adj(11,12) = 1
  adj(11,13) = 1

  adj(12,1) = 1

  adj(13,1) = 1
  adj(13,12) = 1

  call digraph_adjacency_print ( adj, nnode, '  The digraph' )
 
  call digraph_adjacency_components ( adj, nnode, ncomp, comp, dad, order )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of components = ', ncomp
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '  Node, Dad, Component, Order'
  write ( *, '(a)' ) ' '

  do i = 1, nnode
    write ( *, '(5i8)' ) i, dad(i), comp(i), order(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The correct components are:'
  write ( *, '(a)' ) ' (1,11,12,13), (2), (3,4,5), (6,7,8,9,10).'
!
!  Compute a reordering of the nodes.
!
  do i = 1, nnode
    order(i) = i
  end do

  do i = 2, nnode
    do j = 1, i - 1
      if ( comp(j) > comp(i) .or. &
         ( comp(j) == comp(i) .and. order(j) > order(i) ) ) then
        call i4_swap ( comp(j), comp(i) )
        call i4_swap ( order(j), order(i) )
      end if
    end do
  end do

  call i4vec2_print ( nnode, comp, order, '  I, Component(I), Node(I)' )

  call perm_inv ( nnode, order )

  call i4mat_perm ( nnode, adj, order )

  call digraph_adjacency_print ( adj, nnode, '  The graph:' )

  return
end 
subroutine digraph_adjacency_cycle_test ( )

!*****************************************************************************80
!
!! digraph_adjacency_cycle_test() tests digraph_adjacency_cycle().
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
  integer adj2(nnode,nnode)
  integer dad(nnode)
  integer i
  integer nedge
  integer order(nnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'digraph_adjacency_cycle_test():'
  write ( *, '(a)' ) '  digraph_adjacency_cycle() searches for cycles in a digraph.'

  call digraph_adjacency_example_cycler ( adj )

  call digraph_adjacency_print ( adj, nnode, '  The digraph:' )
!
!  Count the edges.
!
  call digraph_adjacency_edge_count ( adj, nnode, nedge )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of edges is ', nedge

  call digraph_adjacency_cycle ( adj, nnode, adj2, dad, order )

  call i4vec2_print ( nnode, dad, order, '  Node, Dad, Order' )

  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '  Adjacency matrix with cycles marked.'
  write ( *, '(a)' ) ' '
  
  do i = 1, nnode
    write ( *, '(10i3)' ) adj2(i,1:nnode)
  end do

  return
end
subroutine digraph_adjacency_degree_test ( )

!*****************************************************************************80
!
!! digraph_adjacency_degree_test tests digraph_adjacency_degree().
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
  integer indegree(nnode)
  integer outdegree(nnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'digraph_adjacency_degree_test():'
  write ( *, '(a)' ) '  digraph_adjacency_degree() computes the degree of the nodes;'

  call digraph_adjacency_example_cycler ( adj )

  call digraph_adjacency_print ( adj, nnode, '  The digraph:' )

  call digraph_adjacency_degree ( adj, nnode, indegree, outdegree )

  call i4vec2_print ( nnode, indegree, outdegree, '  Node, In/Outdegree' )

  return
end
subroutine digraph_adjacency_degree_max_test ( )

!*****************************************************************************80
!
!! digraph_adjacency_degree_max_test() tests digraph_adjacency_degree_max().
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
  integer degree_max
  integer indegree_max
  integer outdegree_max

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'digraph_adjacency_degree_max_test():'
  write ( *, '(a)' ) '  digraph_adjacency_degree_max_test() computes the maximum'
  write ( *, '(a)' ) '  degree of the nodes;'

  call digraph_adjacency_example_cycler ( adj )

  call digraph_adjacency_print ( adj, nnode, '  The digraph:' )

  call digraph_adjacency_degree_max ( adj, nnode, indegree_max, outdegree_max, &
    degree_max )

  write ( *, '(a)' ) ' ' 
  write ( *, '(a,i8)' ) '  Maximum  indegree is ',  indegree_max
  write ( *, '(a,i8)' ) '  Maximum outdegree is ', outdegree_max
  write ( *, '(a,i8)' ) '  Maximum    degree is ',    degree_max
  write ( *, '(a)' ) ' '

  return
end
subroutine digraph_adjacency_degree_sequence_test ( )

!*****************************************************************************80
!
!! digraph_adjacency_degree_sequence_test() tests digraph_adjacency_degree_sequence().
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
  integer indegree_seq(nnode)
  integer outdegree_seq(nnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'digraph_adjacency_degree_sequence_test():'
  write ( *, '(a)' ) '  digraph_adjacency_degree_sequence() computes the degree sequence.'

  call digraph_adjacency_example_cycler ( adj )

  call digraph_adjacency_print ( adj, nnode, '  The digraph:' )

  call digraph_adjacency_degree_sequence ( adj, nnode, indegree_seq, outdegree_seq )

  call i4vec2_print ( nnode, indegree_seq, outdegree_seq, &
    '  Node, In/Outdegree sequence' )

  return
end
subroutine digraph_adjacency_eigen_test ( )

!*****************************************************************************80
!
!! digraph_adjacency_eigen_test() tests digraph_adjacency_eigen().
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

  integer, parameter :: nnode = 9

  integer adj(nnode,nnode)
  real ( kind = rk ) eigeni(nnode)
  real ( kind = rk ) eigenr(nnode)
  integer neigen

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'digraph_adjacency_eigen_test()'
  write ( *, '(a)' ) '  digraph_adjacency_eigen() computes the eigenvalues of a digraph.'

  call digraph_adjacency_example_cycler ( adj )

  call digraph_adjacency_print ( adj, nnode, '  The digraph:' )

  call digraph_adjacency_eigen ( adj, nnode, neigen, eigenr, eigeni )

  call r8vec2_print ( neigen, eigenr, eigeni, &
    '  Real and imaginary parts of eigenvalues:' )

  if ( neigen < nnode ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Warning!  Not all eigenvalues were computed.'
  end if

  return
end
subroutine digraph_adjacency_ham_next_test ( )

!*****************************************************************************80
!
!! digraph_adjacency_ham_next_test() tests digraph_adjacency_ham_next().
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

  integer, allocatable :: adj(:,:)
  integer, allocatable :: circuit(:)
  integer i
  integer j
  logical more
  integer nnode
  integer, allocatable :: ncan(:)
  integer stack(maxstack)
  integer test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'digraph_adjacency_ham_next_test():'
  write ( *, '(a)' ) '  digraph_adjacency_ham_next() produces Hamilton circuits in'
  write ( *, '(a)' ) '  a digraph, one at a time.'
 
  do test = 1, 2

    if ( test == 1 ) then

      nnode = 20
      allocate ( adj(1:nnode,1:nnode) )
      allocate ( circuit(1:nnode) )
      allocate ( ncan(1:nnode) )

      adj(1:nnode,1:nnode) = 0

      adj(1,8) = 1
      adj(1,2) = 1
      adj(1,20) = 1
      adj(2,3) = 1
      adj(2,15) = 1
      adj(3,7) = 1
      adj(3,4) = 1
      adj(4,5) = 1
      adj(4,14) = 1
      adj(5,6) = 1
      adj(5,12) = 1
      adj(6,10) = 1
      adj(6,7) = 1
      adj(7,8) = 1
      adj(8,9) = 1
      adj(9,10) = 1
      adj(9,19) = 1
      adj(10,11) = 1
      adj(11,12) = 1
      adj(11,18) = 1
      adj(12,13) = 1
      adj(13,14) = 1
      adj(13,17) = 1
      adj(14,15) = 1
      adj(15,16) = 1
      adj(16,17) = 1
      adj(16,20) = 1
      adj(17,18) = 1
      adj(18,19) = 1
      adj(19,20) = 1
 
      do i = 1, nnode - 1
        do j = i + 1, nnode
          if ( adj(i,j) == 1 ) then
            adj(j,i) = 1
          end if
        end do
      end do

    else

      nnode = 9
      allocate ( adj(1:nnode,1:nnode) )
      allocate ( circuit(1:nnode) )
      allocate ( ncan(1:nnode) )
      adj(1:nnode,1:nnode) = 0

      adj(1,2) = 1
      adj(1,6) = 1
 
      adj(2,3) = 1
      adj(2,5) = 1
 
      adj(3,4) = 1
 
      adj(4,1) = 1
      adj(4,5) = 1
      adj(4,8) = 1
 
      adj(5,1) = 1
      adj(5,2) = 1
      adj(5,3) = 1
      adj(5,4) = 1
      adj(5,7) = 1
      adj(5,8) = 1
      adj(5,9) = 1
 
      adj(6,3) = 1
      adj(6,5) = 1
      adj(6,8) = 1
 
      adj(7,2) = 1
      adj(7,4) = 1
      adj(7,5) = 1
 
      adj(8,4) = 1
      adj(8,5) = 1
      adj(8,6) = 1
      adj(8,9) = 1
 
      adj(9,5) = 1
    end if

    call digraph_adjacency_print ( adj, nnode, '  The digraph:' )
 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Circuits:'
    write ( *, '(a)' ) ' '
    i = 0

    more = .false.

    do

      call digraph_adjacency_ham_next ( adj, nnode, circuit, stack, maxstack, &
        ncan, more )

      if ( .not. more ) then
        exit
      end if

      i = i + 1
      write ( *, '(i3,2x,20i3)' ) i, circuit(1:nnode)

    end do

    deallocate ( adj )
    deallocate ( circuit )
    deallocate ( ncan )
 
  end do

  return
end
subroutine digraph_adjacency_ham_next_brute_test ( )

!*****************************************************************************80
!
!! digraph_adjacency_ham_next_brute_test() tests digraph_adjacency_ham_next_brute().
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
  write ( *, '(a)' ) 'digraph_adjacency_ham_next_brute_test():'
  write ( *, '(a)' ) '  digraph_adjacency_ham_next_brute() seeks circuits'
  write ( *, '(a)' ) '  in a directed graph which visit every node.'
  write ( *, '(a)' ) '  A brute force algorithm is used.'

  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1
  adj(1,6) = 1
 
  adj(2,3) = 1
  adj(2,5) = 1
 
  adj(3,4) = 1
 
  adj(4,1) = 1
  adj(4,5) = 1
  adj(4,8) = 1
 
  adj(5,1) = 1
  adj(5,2) = 1
  adj(5,3) = 1
  adj(5,4) = 1
  adj(5,7) = 1
  adj(5,8) = 1
  adj(5,9) = 1
 
  adj(6,3) = 1
  adj(6,5) = 1
  adj(6,8) = 1
 
  adj(7,2) = 1
  adj(7,4) = 1
  adj(7,5) = 1
 
  adj(8,4) = 1
  adj(8,5) = 1
  adj(8,6) = 1
  adj(8,9) = 1
 
  adj(9,5) = 1

  call digraph_adjacency_print ( adj, nnode, '  The digraph:' )
 
  iset = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Circuits:'
  write ( *, '(a)' ) ' '
  i = 0

  do
 
    call digraph_adjacency_ham_next_brute ( adj, nnode, circuit, iset )
 
    if ( iset == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  No more circuits were found.'
      exit
    end if

    i = i + 1
    write ( *, '(i3,2x,20i3)' ) i, circuit(1:nnode)

  end do
 
  return
end
subroutine digraph_adjacency_ham_path_next_brute_test ( )

!*****************************************************************************80
!
!! digraph_adjacency_ham_path_next_brute_test() tests digraph_adjacency_ham_path_next_brute().
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

  integer i
  integer adj(nnode,nnode)
  integer iset
  integer j
  integer path(nnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'digraph_adjacency_ham_path_next_brute_test()'
  write ( *, '(a)' ) '  digraph_adjacency_ham_path_next_brute() seeks paths in a'
  write ( *, '(a)' ) '  digraph which visit every node once.'
  write ( *, '(a)' ) '  A brute force algorithm is used.'
!
!  Initialize the adjacency matrix to the identity.
!
  do i = 1, nnode
    do j = 1, nnode
      if ( i == j ) then
        adj(i,j) = 1
      else
        adj(i,j) = 0
      end if
    end do
  end do
!
!  Add entries for specific edges.  This is a directed graph.
!  ADJ(I, j) = 1 means there's a edge from I to J.
!
  adj(1,2) = 1
  adj(1,4) = 1
 
  adj(2,4) = 1
 
  adj(3,1) = 1
  adj(3,4) = 1
 
  adj(4,2) = 1

  call digraph_adjacency_print ( adj, nnode, '  The digraph:' )
 
  iset = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Paths:'
  write ( *, '(a)' ) ' '
  i = 0

  do
 
    call digraph_adjacency_ham_path_next_brute ( adj, nnode, path, iset )
 
    if ( iset == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  No more paths were found.'
      exit
    end if

    i = i + 1
    write ( *, '(i3,2x,20i3)' ) i, path(1:nnode)

  end do
 
  return
end
subroutine digraph_adjacency_is_edge_connected_test ( )

!*****************************************************************************80
!
!! digraph_adjacency_is_edge_connected_test() tests digraph_adjacency_is_edge_connected();
!
!  Discussion:
!
!    Here is a picture of the digraph.
!
!    1-->--2
!    |     |
!    A     A
!    |     |
!    4--<--3
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

  integer adj(nnode,nnode)
  integer result

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'digraph_adjacency_is_edge_connected_test():'
  write ( *, '(a)' ) '  digraph_adjacency_is_edge_connected() reports if a'
  write ( *, '(a)' ) '  digraph is edgewise connected;'

  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1
  adj(3,2) = 1
  adj(3,4) = 1
  adj(4,1) = 1

  call digraph_adjacency_print ( adj, nnode, '  The digraph:' )

  call digraph_adjacency_is_edge_connected ( adj, nnode, result )

  write ( *, '(a)' ) ' '
  if ( result == 0 ) then
    write ( *, '(a)' ) '  The digraph is NOT edgewise connected.'
  else
    write ( *, '(a)' ) '  The digraph IS edgewise connected.'
  end if

  return
end 
subroutine digraph_adjacency_is_eulerian_test ( )

!*****************************************************************************80
!
!! digraph_adjacency_is_eulerian_test() tests digraph_adjacency_is_eulerian();
!
!  Discussion:
!
!    Here is a picture of the digraph:
!
!    1->---2-->---3
!        A V      V
!      6<--5--<---4
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
  integer result

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'digraph_adjacency_is_eulerian_test():'
  write ( *, '(a)' ) '  digraph_adjacency_is_eulerian() reports if a digraph is Eulerian;'

  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1
  adj(2,3) = 1
  adj(3,4) = 1
  adj(4,5) = 1
  adj(5,6) = 1
  adj(6,2) = 1

  call digraph_adjacency_print ( adj, nnode, '  The digraph:' )

  call digraph_adjacency_is_eulerian ( adj, nnode, result )

  write ( *, '(a)' ) ' '
  if ( result == 0 ) then
    write ( *, '(a)' ) '  The digraph is NOT Eulerian.'
  else if ( result == 1 ) then
    write ( *, '(a)' ) '  The digraph IS path Eulerian.'
  else if ( result == 2 ) then
    write ( *, '(a)' ) '  The digraph IS circuit Eulerian.'
  end if

  return
end
subroutine digraph_adjacency_is_strongly_connected_test ( )

!*****************************************************************************80
!
!! digraph_adjacency_is_strongly_connected_test() tests digraph_adjacency_is_strongly_connected();
!
!  Discussion:
!
!    Here are pictures of the digraphs:
!
!  1)
!
!    1-->--2
!    |     |
!    A     A
!    |     |
!    4--<--3
!
!  2)
!
!    1-->--2-->--3-->--4
!    |     |     |     |
!    A     V     A     V
!    |     |     |     |
!    5--<--6     7--<--8
!
!  3)
!
!    1-->--2-->--3-->--4
!    |     |     |     |
!    A     V     A     V
!    |     |     |     |
!    5--<--6--<--7--<--8
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

  integer, allocatable, dimension (:,:) :: adj
  integer nnode
  integer result

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'digraph_adjacency_is_strongly_connected_test():'
  write ( *, '(a)' ) '  digraph_adjacency_is_strongly_connected() reports if a'
  write ( *, '(a)' ) '  digraph is strongly connected;'
!
!  Test 1
!
  nnode = 4
  allocate ( adj(1:nnode,1:nnode) )

  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1
  adj(3,2) = 1
  adj(3,4) = 1
  adj(4,1) = 1

  call digraph_adjacency_print ( adj, nnode, '  The digraph:' )

  call digraph_adjacency_is_strongly_connected ( adj, nnode, result )

  write ( *, '(a)' ) ' '
  if ( result == 0 ) then
    write ( *, '(a)' ) '  The digraph is NOT strongly connected.'
  else
    write ( *, '(a)' ) '  The digraph IS strongly connected.'
  end if

  deallocate ( adj )
!
!  Test 2
!
  nnode = 8
  allocate ( adj(1:nnode,1:nnode) )
  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1
  adj(2,3) = 1
  adj(2,6) = 1
  adj(6,5) = 1
  adj(5,1) = 1
  adj(3,4) = 1
  adj(4,8) = 1
  adj(8,7) = 1
  adj(7,3) = 1

  call digraph_adjacency_print ( adj, nnode, '  The digraph:' )

  call digraph_adjacency_is_strongly_connected ( adj, nnode, result )

  write ( *, '(a)' ) ' '
  if ( result == 0 ) then
    write ( *, '(a)' ) '  The digraph is NOT strongly connected.'
  else
    write ( *, '(a)' ) '  The digraph IS strongly connected.'
  end if

  deallocate ( adj )
!
!  Test 3
!
  nnode = 8
  allocate ( adj(1:nnode,1:nnode) )

  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1
  adj(2,3) = 1
  adj(2,6) = 1
  adj(6,5) = 1
  adj(5,1) = 1
  adj(3,4) = 1
  adj(4,8) = 1
  adj(8,7) = 1
  adj(7,3) = 1
  adj(7,6) = 1

  call digraph_adjacency_print ( adj, nnode, '  The digraph:' )

  call digraph_adjacency_is_strongly_connected ( adj, nnode, result )

  write ( *, '(a)' ) ' '
  if ( result == 0 ) then
    write ( *, '(a)' ) '  The digraph is NOT strongly connected.'
  else
    write ( *, '(a)' ) '  The digraph IS strongly connected.'
  end if

  deallocate ( adj )

  return
end 
subroutine digraph_adjacency_is_tournament_test ( )

!*****************************************************************************80
!
!! digraph_adjacency_is_tournament_test() tests digraph_adjacency_is_tournament();
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
  integer result

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'digraph_adjacency_is_tournament_test():'
  write ( *, '(a)' ) '  digraph_adjacency_is_tournament() reports if a'
  write ( *, '(a)' ) '  digraph is a tournament.'

  call digraph_adjacency_tournament_random ( nnode, adj )

  call digraph_adjacency_print ( adj, nnode, '  A random tournament digraph:' )

  call digraph_adjacency_is_tournament ( adj, nnode, result )

  write ( *, '(a)' ) ' '
  if ( result == 0 ) then
    write ( *, '(a)' ) '  The digraph is NOT a tournament.'
  else
    write ( *, '(a)' ) '  The digraph IS a tournament.'
  end if

  return
end 
subroutine digraph_adjacency_is_transitive_test ( )

!*****************************************************************************80
!
!! digraph_adjacency_is_transitive_test() tests digraph_adjacency_is_transitive();
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

  integer, parameter :: nnode = 12

  integer adj(nnode,nnode)
  integer result

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'digraph_adjacency_is_transitive_test():'
  write ( *, '(a)' ) '  digraph_adjacency_is_transitive() reports if a'
  write ( *, '(a)' ) '  digraph is transitive;'

  call digraph_adjacency_example_sixty ( adj )

  call digraph_adjacency_print ( adj, nnode, '  The digraph:' )

  call digraph_adjacency_is_transitive ( adj, nnode, result )

  write ( *, '(a)' ) ' '
  if ( result == 0 ) then
    write ( *, '(a)' ) '  The digraph is NOT transitive.'
  else
    write ( *, '(a)' ) '  The digraph IS transitive.'
  end if

  return
end 
subroutine digraph_adjacency_random_test ( )

!*****************************************************************************80
!
!! digraph_adjacency_random_test() tests digraph_adjacency_random();
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
  integer nedge

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'digraph_adjacency_random_test():'
  write ( *, '(a)' ) '  digraph_adjacency_random() returns a random digraph.'
  write ( *, '(a)' ) ' '

  nedge = 10
  write ( *, '(a,i8)' ) '  Number of edges requested = ', nedge

  call digraph_adjacency_random ( nnode, nedge, adj )

  call digraph_adjacency_print ( adj, nnode, '  The digraph:' )
!
!  Count the edges.
!
  call digraph_adjacency_edge_count ( adj, nnode, nedge )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of edges is ', nedge

  return
end
subroutine digraph_adjacency_reduce_test ( )

!*****************************************************************************80
!
!! digraph_adjacency_reduce_test() tests digraph_adjacency_reduce().
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

  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1
  adj(1,6) = 1
  adj(1,7) = 1

  adj(3,1) = 1

  adj(4,6) = 1

  adj(5,4) = 1

  adj(6,5) = 1

  adj(7,3) = 1
  adj(7,5) = 1
  adj(7,10) = 1

  adj(8,7) = 1
  adj(8,9) = 1

  adj(9,8) = 1

  adj(10,11) = 1
  adj(10,12) = 1
  adj(10,13) = 1

  adj(12,7) = 1
  adj(12,13) = 1

  adj(13,12) = 1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'digraph_adjacency_reduce_test():'
  write ( *, '(a)' ) '  digraph_adjacency_reduce() finds the transitive '
  write ( *, '(a)' ) '  reduction of a digraph.'
  write ( *, '(a)' ) ' '

  call digraph_adjacency_print ( adj, nnode, '  Adjacency matrix for G:' )

  call digraph_adjacency_reduce ( adj, nnode )

  call digraph_adjacency_print ( adj, nnode, &
    '  Adjacency matrix for the transitive reduction of G:' )

  return
end
subroutine digraph_adjacency_to_digraph_arc_test ( )

!*****************************************************************************80
!
!! digraph_adjacency_to_digraph_arc_test() tests digraph_adjacency_to_digraph_arc();
!
!    1->---2-->---3
!    |     |      |
!    A     V      V
!    |     |      |
!    6--<--5--<---4
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
  integer, parameter :: maxarc = 10

  integer adj(nnode,nnode)
  integer inode(maxarc)
  integer jnode(maxarc)
  integer narc

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'digraph_adjacency_to_digraph_arc_test():'
  write ( *, '(a)' ) '  digraph_adjacency_to_digraph_arc() converts a digraph in'
  write ( *, '(a)' ) '  adjacency form to arc list form;'

  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1
  adj(2,3) = 1
  adj(3,4) = 1
  adj(4,5) = 1
  adj(5,6) = 1
  adj(6,2) = 1

  call digraph_adjacency_print ( adj, nnode, '  The digraph in adjacency form:' )

  call digraph_adjacency_to_digraph_arc ( adj, nnode, maxarc, narc, &
    inode, jnode )

  call digraph_arc_print ( narc, inode, jnode, &
    '  The digraph in arc list form:' )

  return
end
subroutine digraph_adjacency_to_incidence_test ( )

!*****************************************************************************80
!
!! digraph_adjacency_to_incidence_test() tests digraph_adjacency_to_incidence();
!
!    1->---2-->---3
!        A V      V
!      6<--5--<---4
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
  integer, parameter :: maxarc = 10

  integer adj(nnode,nnode)
  integer inc(nnode,maxarc)
  integer narc

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'digraph_adjacency_to_incidence_test():'
  write ( *, '(a)' ) '  digraph_adjacency_to_incidence() converts a digraph in'
  write ( *, '(a)' ) '  adjacency form to incidence matrix form;'

  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1
  adj(2,3) = 1
  adj(3,4) = 1
  adj(4,5) = 1
  adj(5,6) = 1
  adj(6,2) = 1

  call digraph_adjacency_print ( adj, nnode, '  The digraph in adjacency form:' )

  call digraph_adjacency_to_incidence ( adj, nnode, maxarc, narc, inc )

  call digraph_incidence_print ( nnode, narc, inc, &
    '  The digraph in incidence form:' )

  return
end
subroutine digraph_adjacency_top_sort_test ( )

!*****************************************************************************80
!
!! digraph_adjacency_top_sort_test() tests digraph_adjacency_top_sort().
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
  integer node_list(nnode)
  integer order(nnode)
  integer visit(nnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'digraph_adjacency_top_sort_test()'
  write ( *, '(a)' ) '  digraph_adjacency_top_sort() does a topological sort'
  write ( *, '(a)' ) '  of an acyclic digraph.'

  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1
  adj(1,3) = 1
  adj(1,6) = 1

  adj(5,4) = 1

  adj(6,4) = 1
  adj(6,5) = 1

  adj(7,3) = 1
  adj(7,5) = 1
  adj(7,8) = 1

  adj(8,9) = 1

  adj(10,7) = 1
  adj(10,11) = 1
  adj(10,12) = 1
  adj(10,13) = 1

  adj(12,7) = 1
  adj(12,13) = 1

  call digraph_adjacency_print ( adj, nnode, '  The digraph:' )
 
  call digraph_adjacency_top_sort ( adj, nnode, dad, visit, node_list )

  call i4vec_print ( nnode, dad, '  Nodes and "Dads":' )

  call i4vec_print ( nnode, visit, '  Nodes and order of visit:' )

  call i4vec_print ( nnode, node_list, '  Nodes and reverse topological order:' )
!
!  Invert the listing to get a permutation.
!
  order(1:nnode) = node_list(1:nnode)

  call perm_inv ( nnode, order )
!
!  Apply reordering and print adjacency matrix.
!
  call i4mat_perm ( nnode, adj, order )

  call digraph_adjacency_print ( adj, nnode, '  The reordered digraph:' )

  return
end
subroutine digraph_adjacency_tournament_random_test ( )

!*****************************************************************************80
!
!! digraph_adjacency_tournament_random_test() tests digraph_adjacency_tournament_random();
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

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'digraph_adjacency_tournament_random_test():'
  write ( *, '(a)' ) '  digraph_adjacency_tournament_random() returns a random'
  write ( *, '(a)' ) '  tournament digraph.'

  call digraph_adjacency_tournament_random ( nnode, adj )

  call digraph_adjacency_print ( adj, nnode, '  A random tournament digraph:' )

  return
end 
subroutine digraph_adjacency_transitive_closure_test ( )

!*****************************************************************************80
!
!! digraph_adjacency_transitive_closure_test() tests digraph_adjacency_transitive_closure().
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

  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1
  adj(1,6) = 1
  adj(1,7) = 1

  adj(3,1) = 1

  adj(4,6) = 1

  adj(5,4) = 1

  adj(6,5) = 1

  adj(7,3) = 1
  adj(7,5) = 1
  adj(7,10) = 1

  adj(8,7) = 1
  adj(8,9) = 1

  adj(9,8) = 1

  adj(10,11) = 1
  adj(10,12) = 1
  adj(10,13) = 1

  adj(12,7) = 1
  adj(12,13) = 1

  adj(13,12) = 1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'digraph_adjacency_transitive_closure_test():'
  write ( *, '(a)' ) '  digraph_adjacency_transitive_closure() finds the transitive '
  write ( *, '(a)' ) '  closure of a digraph;'
  write ( *, '(a)' ) ' '

  call digraph_adjacency_print ( adj, nnode, '  Adjacency matrix for G:' )

  call digraph_adjacency_transitive_closure ( adj, nnode )

  call digraph_adjacency_print ( adj, nnode, &
    '  Adjacency matrix for H, the transitive closure of G:' )

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

