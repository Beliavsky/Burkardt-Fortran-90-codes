program main

!*****************************************************************************80
!
!! graph_theory_test() tests graph_theory().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 March 2023
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_theory_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  graph_theory() implements graph algorithms.'

  call degree_sequence_is_graphic_test ( )
  call degree_sequence_to_graph_adj_test ( )

  call face_check_test ( )

  call test060 ( )

  call grf_read_test ( )

  call network_flow_max_test ( )

  call node_relax_test ( )

  call perm_inc_test ( )

  call poly_to_tri_test ( )

  call test0695 ( )
  call test0696 ( )
  call test0697 ( )

  call sort_heap_external_test ( )

  call span_forest_test ( )
  call span_tree_next_test ( )
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
subroutine degree_sequence_is_graphic_test ( )

!*****************************************************************************80
!
!! degree_sequence_is_graphic_test() tests degree_sequence_is_graphic().
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
  integer, parameter :: ntest = 5

  integer degree_seq(nnode)
  integer i
  integer result

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'degree_sequence_is_graphic_test():'
  write ( *, '(a)' ) '  degree_sequence_is_graphic() reports whether'
  write ( *, '(a)' ) '  a given sequence can represent the degree'
  write ( *, '(a)' ) '  sequence of a graph.'
  write ( *, '(a)' ) ' '

  do i = 1, ntest

    call i4vec_uniform_ab ( nnode, 1, nnode - 1, degree_seq )

    call i4vec_sort_heap_d ( nnode, degree_seq )

    call i4vec_print ( nnode, degree_seq, '  The degree sequence:' )

    call degree_sequence_is_graphic ( nnode, degree_seq, result )

    write ( *, '(a)' ) ' '
    if ( result == 0 ) then
      write ( *, '(a)' ) '  The sequence is NOT graphic.'
    else if ( result == 1 ) then
      write ( *, '(a)' ) '  The sequence IS graphic.'
    end if

  end do

  return
end
subroutine degree_sequence_to_graph_adj_test ( )

!*****************************************************************************80
!
!! degree_sequence_to_graph_adj_test() tests degree_sequence_to_graph_adj().
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
  integer ierror
  integer, dimension ( nnode ) :: seq = (/ 5, 5, 4, 3, 3, 2 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'degree_sequence_to_graph_adj_test()'
  write ( *, '(a)' ) '  degree_sequence_to_graph_adj() is given a degree'
  write ( *, '(a)' ) '  sequence, and constructs the adjancency'
  write ( *, '(a)' ) '  matrix of a corresponding graph.'

  call i4vec_print ( nnode, seq, '  The degree sequence:' )

  call degree_sequence_to_graph_adj ( nnode, seq, adj, ierror )

  call graph_adj_print ( adj, nnode, '  The graph:' )

  return
end
subroutine face_check_test ( )

!*****************************************************************************80
!
!! face_check_test() tests face_check().
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

  integer, parameter :: max_edge = 30
  integer, parameter :: max_order = 4
  integer, parameter :: max_face = 10

  integer edge(4,max_edge)
  integer face(max_order,max_face)
  integer face_object(max_face)
  integer face_order(max_face)
  integer face_rank(max_face)
  integer face_tier(max_face)
  integer i
  integer j
  integer num_edge
  integer num_face
  integer num_object

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'face_check_test()'
  write ( *, '(a)' ) '  face_check() checks faces;'
!
!  Get the problem data.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  max_face =  ', max_face
  write ( *, '(a,i8)' ) '  max_order = ', max_order
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Get a test example'

  call face_example_pieces ( face, face_order, max_face, max_order, num_face )
!
!  List the problem data.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Face, Order, Nodes'
  write ( *, '(a)' ) ' '
  do i = 1, num_face
    write ( *, '(10i3)' ) i, face_order(i), ( face(j,i), j = 1, face_order(i) )
  end do
!
!  Check the problem data.
!
  call face_check ( edge, face, face_object, face_order, face_rank, &
    face_tier, max_edge, max_order, num_edge, num_face, num_object )

  return
end
subroutine test060 ( )

!*****************************************************************************80
!
!! TEST060 tests VLA_TO_GRAPH_ARC, GRAPH_ARC_FACE, FACE_TO_IV;
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

  integer, parameter :: maxedge = 1000
  integer, parameter :: maxface = 2000
  integer, parameter :: maxnode = 1000
  integer, parameter :: maxorder = 20

  integer face(maxorder,maxface)
  integer face_count(maxedge)
  integer face_order(maxface)
  character ( len = 80 ) :: file_in = 'fish_lines.vla'
  character ( len = 80 ) :: file_out = 'fish_faces.iv'
  integer ierror
  integer iface(maxedge)
  integer inode(maxedge)
  integer jface(maxedge)
  integer jnode(maxedge)
  integer nedge
  integer nface
  integer nnode
  real ( kind = rk ) x(maxnode)
  real ( kind = rk ) y(maxnode)
  real ( kind = rk ) z(maxnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST060'
  write ( *, '(a)' ) '  VLA_TO_GRAPH_ARC converts VLA edge data to a'
  write ( *, '(a)' ) '  graph defined by arcs;'
  write ( *, '(a)' ) '  GRAPH_ARC_FACE constructs the faces of an orientable graph.'
  write ( *, '(a)' ) '  FACE_TO_IV writes face data to an IV file.'
!
!  Get the edge array for the graph.
!
  call vla_to_graph_arc ( file_in, maxedge, maxnode, nedge, nnode, inode, &
    jnode, x, y, z, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) 'TEST060 - Fatal error!'
    write ( *, '(a)' ) '  VLA_TO_GRAPH_ARC returned an error.'
    stop 1
  end if
!
!  Sort the edge array.

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sort the edges:'

  call graph_arc_edge_sort ( nedge, inode, jnode )
!
!  Determine the faces.
!
  write ( *, '(a)' ) '  Determine the faces:'

  call graph_arc_face ( face, face_count, face_order, iface, jface, &
    inode, jnode, maxface, maxorder, nedge, nface, nnode )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of faces found was ', nface
  write ( *, '(a,i8)' ) '  Euler predicted ', nedge + 2 - nnode
!
!  Write the faces to an IV file.
!
  call face_to_iv ( file_out, face, face_order, inode, jnode, &
    nedge, maxnode, maxface, maxorder, nnode, nface, x, y, z )

  return
end
subroutine grf_read_test ( )

!*****************************************************************************80
!
!! grf_read_test() tests grf_read().
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

  integer, parameter :: maxedge = 500
  integer, parameter :: maxnode = 100

  integer adj(maxnode,maxnode)
  character ( len = 80 ) :: file_grf = 'knightstour.grf'
  character ( len = 80 ) :: file_ps = 'knightstour.eps'
  integer i
  integer inode(maxedge)
  integer jnode(maxedge)
  integer nedge
  integer nnode
  real ( kind = rk ) x(maxnode)
  real ( kind = rk ) y(maxnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'grf_read_test():'
  write ( *, '(a)' ) '  grf_read() reads a GRF file,'
  write ( *, '(a)' ) '  graph_arc_to_ps() writes a PostScript version of it.'

  call grf_read ( file_grf, inode, jnode, maxedge, maxnode, nedge, nnode, x, y )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Node, X, Y'
  write ( *, '(a)' ) ' '

  do i = 1, nnode
    write ( *, '(i8,2g14.6)' ) i, x(i), y(i)
  end do

  call graph_arc_to_graph_adj ( nedge, inode, jnode, adj, nnode )

  call graph_adj_print ( adj, nnode, '  The graph:' )
!
!  Now write out a PostScript version.
!
  call graph_arc_to_ps ( file_ps, inode, jnode, nedge, nnode, x, y )

  return
end
subroutine network_flow_max_test ( )

!*****************************************************************************80
!
!! network_flow_max_test() tests network_flow_max().
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
  integer, parameter :: nedge = 20

  integer i
  integer icut(nnode)
  integer icpflo(2,nedge)
  integer iendpt(2,nedge)
  integer :: isink = 6
  integer :: isorce = 1
  integer j
  integer node_flow(nnode)

  data ( ( iendpt(i,j), j = 1, nedge ), i = 1, 2 ) / &
    1,2, 1,3, 2,3, 2,4, 2,5, 3,4, 3,5, 4,5, 4,6, 5,6, &
    2,1, 3,1, 3,2, 4,2, 5,2, 4,3, 5,3, 5,4, 6,4, 6,5 /
 
  data ( ( icpflo(i,j), j = 1, nedge ), i = 1, 2 ) / &
    3,0,7,0,2,0,5,0,4,0,1,0,4,0,2,0,8,0,3,0, &
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'network_flow_max_test():'
  write ( *, '(a)' ) '  network_flow_max() finds the maximum flow on a network.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The source is node ', isorce
  write ( *, '(a,i8)' ) '  The sink is node   ', isink
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Endpoint array:'
  write ( *, '(a)' ) ' '
  write ( *, '(20i3)' ) ( iendpt(1,i), i = 1, nedge )
  write ( *, '(20i3)' ) ( iendpt(2,i), i = 1, nedge )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Input edge capacity array:'
  write ( *, '(a)' ) ' '
  write ( *, '(20i3)' ) ( icpflo(1,i), i = 1, nedge)
 
  call network_flow_max ( nnode, nedge, iendpt, icpflo, isorce, &
    isink, icut, node_flow )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Reordered endpoint array:'
  write ( *, '(a)' ) ' '
  write ( *, '(20i3)' ) ( iendpt(1,i), i = 1, nedge )
  write ( *, '(20i3)' ) ( iendpt(2,i), i = 1, nedge )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Output edge capacity/flow array:'
  write ( *, '(a)' ) ' '
  write ( *, '(20i3)' ) ( icpflo(1,i), i = 1, nedge )
  write ( *, '(20i3)' ) ( icpflo(2,i), i = 1, nedge )

  call i4vec_print ( nnode, icut, '  Minimal node cut vector:' )

  call i4vec_print ( nnode, node_flow, '  Nodal flow vector:' )

  return
end
subroutine node_relax_test ( )

!*****************************************************************************80
!
!! node_relax_test() tests node_relax().
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

  integer, parameter :: max_cor3 = 100
  integer, parameter :: max_face = 100
  integer, parameter :: max_order = 5

  real ( kind = rk ) cor3(3,max_cor3)
  real ( kind = rk ) cor3_new(3,max_cor3)
  integer cor3_num(max_cor3)
  integer face(max_order,max_face)
  integer face_order(max_face)
  integer j
  integer num_cor3
  integer num_face

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'node_relax_test():'
  write ( *, '(a)' ) '  node_relax() smooths a surface.'

  num_cor3 = 8

  cor3(1,1) = 0.0D+00
  cor3(2,1) = 0.0D+00
  cor3(3,1) = 0.0D+00

  cor3(1,2) = 1.0D+00
  cor3(2,2) = 0.0D+00
  cor3(3,2) = 0.0D+00

  cor3(1,3) = 1.0D+00
  cor3(2,3) = 1.0D+00
  cor3(3,3) = 0.0D+00

  cor3(1,4) = 0.0D+00
  cor3(2,4) = 1.0D+00
  cor3(3,4) = 0.0D+00

  cor3(1,5) = 0.0D+00
  cor3(2,5) = 0.0D+00
  cor3(3,5) = 1.0D+00

  cor3(1,6) = 1.0D+00
  cor3(2,6) = 0.0D+00
  cor3(3,6) = 1.0D+00

  cor3(1,7) = 1.0D+00
  cor3(2,7) = 1.0D+00
  cor3(3,7) = 1.0D+00

  cor3(1,8) = 0.0D+00
  cor3(2,8) = 1.0D+00
  cor3(3,8) = 1.0D+00

  num_face = 6

  face(1,1) = 1
  face(2,1) = 4
  face(3,1) = 3
  face(4,1) = 2

  face(1,2) = 2
  face(2,2) = 6
  face(3,2) = 7
  face(4,2) = 3

  face(1,3) = 3
  face(2,3) = 7
  face(3,3) = 8
  face(4,3) = 4

  face(1,4) = 4
  face(2,4) = 8
  face(3,4) = 5
  face(4,4) = 1

  face(1,5) = 1
  face(2,5) = 5
  face(3,5) = 6
  face(4,5) = 2

  face(1,6) = 5
  face(2,6) = 8
  face(3,6) = 7
  face(4,6) = 6

  face_order(1:num_face) = 4

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Old coordinates'
  write ( *, '(a)' ) ' '
  do j = 1, num_cor3
    write ( *, '(i4, 3g14.6)' ) j, cor3(1:3,j)
  end do

  call node_relax ( cor3, cor3_new, cor3_num, face, face_order, max_cor3, &
    max_face, max_order, num_cor3, num_face )

  write ( *, '(a)' ) ' '
  write ( *, '(a)') '  After 1 step'
  write ( *, '(a)' ) ' '

  do j = 1, num_cor3
    write ( *, '(i4, 3g14.6)' ) j, cor3_new(1:3,j)
  end do

  cor3(1:3,1:num_cor3) = cor3_new(1:3,1:num_cor3)

  call node_relax ( cor3, cor3_new, cor3_num, face, face_order, max_cor3, &
    max_face, max_order, num_cor3, num_face )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  After 2 steps'
  write ( *, '(a)' ) ' '

  do j = 1, num_cor3
    write ( *, '(i4, 3g14.6)' ) j, cor3_new(1:3,j)
  end do

  cor3(1:3,1:num_cor3) = cor3_new(1:3,1:num_cor3)

  call node_relax ( cor3, cor3_new, cor3_num, face, face_order, max_cor3, &
    max_face, max_order, num_cor3, num_face )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  After 3 steps'
  write ( *, '(a)' ) ' '

  do j = 1, num_cor3
    write ( *, '(i4, 3g14.6)' ) j, cor3_new(1:3,j)
  end do

  return
end
subroutine perm_inc_test ( )

!*****************************************************************************80
!
!! perm_inc_test() tests perm_inc().
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

  integer, parameter :: n = 4

  integer i
  integer ipos
  integer perm(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'perm_inc_test()'
  write ( *, '(a)' ) '  perm_inc() increments a permutation.'
  write ( *, '(a)' ) ' '

  i = 0
  ipos = 0

  do

    call perm_inc ( perm, ipos, n )

    if ( ipos == 0 ) then
      exit
    end if

    i = i + 1
    write ( *, '(i3,2x,4i2)' ) i, perm(1:n)

  end do

  return
end
subroutine poly_to_tri_test ( )

!*****************************************************************************80
!
!! poly_to_tri_test() tests poly_to_tri().
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

  integer, parameter :: max_face = 20
  integer, parameter :: max_vert = 5

  integer face(max_vert,max_face)
  integer i
  integer ierror
  integer j
  integer num_face
  integer num_vert(max_face)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'poly_to_tri_test():'
  write ( *, '(a)' ) '  poly_to_tri() replaces a polygonal mesh with a'
  write ( *, '(a)' ) '  triangular one.'

  num_face = 4

  num_vert(1) = 4
  face(1,1) = 1
  face(2,1) = 3
  face(3,1) = 5
  face(4,1) = 7

  num_vert(2) = 3
  face(1,2) = 2
  face(2,2) = 3
  face(3,2) = 9

  num_vert(3) = 5
  face(1,3) = 3
  face(2,3) = 7
  face(3,3) = 8
  face(4,3) = 23
  face(5,3) = 2

  num_vert(4) = 4
  face(1,4) = 4
  face(2,4) = 7
  face(3,4) = 8
  face(4,4) = 23

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of faces = ', num_face

  call i4vec_print ( num_face, num_vert, '  Faces and number of vertices:' )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Face   Vertices'
  write ( *, '(a)' ) ' '
  do i = 1, num_face
    write ( *, '(6i8)' ) i, ( face(j,i), j = 1, num_vert(i) )
  end do

  call poly_to_tri ( face, ierror, max_face, max_vert, num_face, num_vert )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The algorithm failed.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Number of faces = ', num_face

    call i4vec_print ( num_face, num_vert, '  Faces and number of vertices:' )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Face   Vertices'
    write ( *, '(a)' ) ' '
    do i = 1, num_face
      write ( *, '(6i8)' ) i, ( face(j,i), j = 1, num_vert(i) )
    end do

  end if

  return
end
subroutine test0695 ( )

!*****************************************************************************80
!
!! TEST0695 tests VLA_TO_GRAPH_ARC, SHAPE_3D_NODES_TO_PS.
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

  integer, parameter :: max_edge = 1000
  integer, parameter :: max_node = 1000

  character ( len = 80 ) :: file_in = 'fish_lines.vla'
  character ( len = 80 ) :: file_out = 'fish_nodes.ps'
  integer ierror
  integer inode(max_edge)
  integer jnode(max_edge)
  integer num_edge
  integer num_node
  real ( kind = rk ) x(max_node)
  real ( kind = rk ) y(max_node)
  real ( kind = rk ) z(max_node)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0695'
  write ( *, '(a)' ) '  VLA_TO_GRAPH_ARC reads a VLA file and converts it'
  write ( *, '(a)' ) '  to a graph defined by an arc list.'
  write ( *, '(a)' ) '  SHAPE_3D_NODES_TO_PS writes the nodes to a PostScript file.'
!
!  Get the edge array for the graph.
!
  call vla_to_graph_arc ( file_in, max_edge, max_node, num_edge, &
    num_node, inode, jnode, x, y, z, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)') '  VLA_TO_GRAPH_ARC returned an error.'
    return
  end if

  call shape_3d_nodes_to_ps ( file_out, num_node, x, y, z )

  return
end
subroutine test0696 ( )

!*****************************************************************************80
!
!! TEST0696 tests VLA_TO_GRAPH_ARC, SHAPE_3D_EDGES_TO_PS.
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

  integer, parameter :: max_edge = 1000
  integer, parameter :: max_face = 2000
  integer, parameter :: max_node = 1000
  integer, parameter :: max_order = 20

  integer face(max_order,max_face)
  integer face_count(max_edge)
  integer face_order(max_face)
  character ( len = 80 ) :: file_in = 'fish_lines.vla'
  character ( len = 80 ) :: file_out = 'fish_edges.ps'
  integer ierror
  integer iface(max_edge)
  integer inode(max_edge)
  integer jface(max_edge)
  integer jnode(max_edge)
  integer num_edge
  integer num_face
  integer num_node
  real ( kind = rk ) x(max_node)
  real ( kind = rk ) y(max_node)
  real ( kind = rk ) z(max_node)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0696'
  write ( *, '(a)' ) '  VLA_TO_GRAPH_ARC reads a VLA file and converts it'
  write ( *, '(a)' ) '  to a graph defined by an arc list.'
  write ( *, '(a)' ) '  SHAPE_3D_EDGES_TO_PS writes the edges to a PostScript file.'
!
!  Get the edge array for the graph.
!
  call vla_to_graph_arc ( file_in, max_edge, max_node, num_edge, &
    num_node, inode, jnode, x, y, z, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  VLA_TO_GRAPH_ARC returned an error.'
    return
  end if
!
!  Sort the edge array.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sort the edges:'

  call graph_arc_edge_sort ( num_edge, inode, jnode )
!
!  Determine the faces.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Determine the faces:'

  call graph_arc_face ( face, face_count, face_order, iface, jface, inode, &
    jnode, max_face, max_order, num_edge, num_face, num_node )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The faces were determined.'
  write ( *, '(a,i8)' ) '    Number of faces found was ', num_face
  write ( *, '(a,i8)' ) '    Euler predicted ', num_edge + 2 - num_node

  call shape_3d_edges_to_ps ( file_out, max_order, num_face, num_node, &
    face, face_order, x, y, z )

  return
end
subroutine test0697 ( )

!*****************************************************************************80
!
!! TEST0697 tests VLA_TO_GRAPH_ARC, SHAPE_3D_FACES_TO_PS.
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

  integer, parameter :: max_edge = 1000
  integer, parameter :: max_face = 2000
  integer, parameter :: max_node = 500
  integer, parameter :: max_order = 20

  integer face(max_order,max_face)
  integer face_count(max_edge)
  integer face_order(max_face)
  character ( len = 80 ) :: file_in = 'fish_lines.vla'
  character ( len = 80 ) :: file_out = 'fish_faces.ps'
  integer ierror
  integer iface(max_edge)
  integer inode(max_edge)
  integer jface(max_edge)
  integer jnode(max_edge)
  integer num_edge
  integer num_face
  integer num_node
  real ( kind = rk ) x(max_node)
  real ( kind = rk ) y(max_node)
  real ( kind = rk ) z(max_node)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0697'
  write ( *, '(a)' ) '  VLA_TO_GRAPH_ARC reads a VLA file and converts it'
  write ( *, '(a)' ) '  to a graph defined by an arc list.'
  write ( *, '(a)' ) '  SHAPE_3D_FACES_TO_PS writes the faces to a PostScript file.'
!
!  Get the edge array for the graph.
!
  call vla_to_graph_arc ( file_in, max_edge, max_node, num_edge, &
    num_node, inode, jnode, x, y, z, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) 'TEST0697 - Error!'
    write ( *, '(a)' ) '  VLA_TO_GRAPH_ARC returned an error.'
    return
  end if
!
!  Sort the edge array.

  call graph_arc_edge_sort ( num_edge, inode, jnode )
!
!  Determine the faces.
!
  call graph_arc_face ( face, face_count, face_order, iface, jface, inode, &
    jnode, max_face, max_order, num_edge, num_face, num_node )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of faces found was ', num_face
  write ( *, '(a,i8)' ) '  Euler predicted ', num_edge + 2 - num_node

  call shape_3d_faces_to_ps ( file_out, max_order, num_face, num_node, &
    face, face_order, x, y, z )

  return
end
subroutine sort_heap_external_test ( )

!*****************************************************************************80
!
!! sort_heap_external_test() tests sort_heap_external().
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

  integer, parameter :: n = 20

  integer a(n)
  integer i
  integer indx
  integer isgn
  integer j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'sort_heap_external_test():'
  write ( *, '(a)' ) '  sort_heap_external() sorts objects externally.'
  write ( *, '(a)' ) ' '

  indx = 0
  i = 0
  j = 0
  isgn = 0

  call i4vec_uniform_ab ( n, 1, n, a )

  call i4vec_print ( n, a, '  Before sorting:' ) 
 
  do

    call sort_heap_external ( n, indx, i, j, isgn )
 
    if ( indx < 0 ) then
      isgn = 1
      if ( a(i) <= a(j) ) then
        isgn = -1
      end if
    else if ( indx > 0 ) then
      call i4_swap ( a(i), a(j) )
    else
      exit
    end if

  end do

  call i4vec_print ( n, a, '  After sorting:' )
 
  return
end
subroutine span_forest_test ( )

!*****************************************************************************80
!
!! span_forest_test() tests span_forest().
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
  integer i
  integer iendpt(2,nedge)
  integer j
  integer k

  data ( ( iendpt(i,j), i = 1, 2 ), j = 1, nedge ) / &
    2,3, 4,7, 1,9, 7,11, 5,8, 2,5, 6,10, 2,8, 3,8, 4,11 /

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'span_forest_test():'
  write ( *, '(a)' ) '  span_forest() constructs a spanning forest for a graph.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Initial end point array:'
  write ( *, '(a)' ) ' '
  write ( *, '(19i4)' ) ( iendpt(1,j), j = 1, nedge )
  write ( *, '(19i4)' ) ( iendpt(2,j), j = 1, nedge )
 
  call span_forest ( nnode, nedge, iendpt, k, component )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Reordered endpoint array:'
  write ( *, '(a)' ) ' '
  write ( *, '(19i4)' ) ( iendpt(1,j), j = 1, nedge )
  write ( *, '(19i4)' ) ( iendpt(2,j), j = 1, nedge )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of connected components = ', k
 
  call i4vec_print ( nnode, component, '  Node, Component' )

  return
end
subroutine span_tree_next_test ( )

!*****************************************************************************80
!
!! span_tree_next_test() tests span_tree_next();
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

  integer, parameter :: nnode = 5
  integer, parameter :: nedge = 10

  integer i
  integer iarray(nnode-1)
  integer iendpt(2,nedge)
  integer j
  integer ncan(nnode-1)
  integer nspan
  integer signal

  data ( ( iendpt(i,j), i = 1, 2 ), j = 1, nedge ) / &
    1,2, 1,3, 1,4, 1,5, 2,3, 2,4, 2,5, 3,4, 3,5, 4,5 /

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'span_tree_next_test():'
  write ( *, '(a)' ) '  span_tree_next() constructs spanning trees'
  write ( *, '(a)' ) '  of a graph using a backtrack search.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Node1   Node2'
  write ( *, '(a)' ) ' '
  do i = 1, nedge
    write ( *, '(3i8)' ) iendpt(1,i), iendpt(2,i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Edges in spanning tree:'
  write ( *, '(a)' ) ' '

  nspan = 0
  signal = 0

  do

    call span_tree_next ( signal, nnode, nedge, iendpt, iarray, ncan )

    if ( signal == 0 ) then 
      exit
    end if

    nspan = nspan + 1
    write ( *, '(i4,4x,5i4)' ) nspan, iarray(1:nnode-1)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of spanning trees found was ', nspan

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

