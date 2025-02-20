program main

!*****************************************************************************80
!
!! graph_dist_test() tests graph_dist().
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
  write ( *, '(a)' ) 'graph_dist_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  graph_dist() implements graph algorithms.'

  call graph_dist_all_test ( )
  call graph_dist_check_test ( )
  call graph_dist_min_span_tree_test ( )
  call graph_dist_min_span_tree_alternate_test ( )
  call graph_dist_min_span_tree2_test ( )
  call graph_dist_min_span_tree3_test ( )
  call graph_dist_one_test ( )
  call graph_dist_pairing_greedy_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_dist_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine graph_dist_all_test ( )

!*****************************************************************************80
!
!! graph_dist_all_test() tests graph_dist_all().
!
!  The graph is:
!
!  N3 --3-- N2 --4-- N4 --5-- N5
!
!     \      |      /
!       6    2     1
!        \   |    /
!
!            N1
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

  real ( kind = rk ) dinfin
  real ( kind = rk ) dist(nnode,nnode)
  integer i
  real ( kind = rk ) path_dist(nnode,nnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_dist_all_test():'
  write ( *, '(a)' ) '  graph_dist_all() computes the distance between'
  write ( *, '(a)' ) '  all pairs of nodes.'
  write ( *, '(a)' ) ' '

  dinfin = 1000.0D+00
 
  dist(1:nnode,1:nnode) = dinfin

  do i = 1, nnode
    dist(i,i) = 0.0D+00
  end do
 
  dist(1,2) = 2.0D+00
  dist(1,3) = 6.0D+00
  dist(1,4) = 1.0D+00

  dist(2,1) = 2.0D+00
  dist(2,3) = 3.0D+00
  dist(2,4) = 4.0D+00

  dist(3,1) = 6.0D+00
  dist(3,2) = 3.0D+00

  dist(4,1) = 1.0D+00
  dist(4,2) = 4.0D+00
  dist(4,5) = 5.0D+00

  dist(5,4) = 5.0D+00
 
  call graph_dist_print ( dist, nnode, &
    '  Immediate node distance matrix:' )

  call graph_dist_all ( dist, dinfin, nnode, path_dist )
 
  call graph_dist_print ( path_dist, nnode, &
    '  Total node distance matrix:' )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Note that "infinity" is represented by ', dinfin
 
  return
end
subroutine graph_dist_check_test ( )

!*****************************************************************************80
!
!! graph_dist_check_test() tests graph_dist_check().
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

  integer, parameter :: nnode = 15

  real ( kind = rk ) a(nnode,nnode)
  integer i
  integer ierror
  integer j

  data ( ( a(i,j), j = 1, nnode ), i = 1, nnode ) / &
     0., 29., 82., 46., 68., 52., 72., 42., 51., 55., 29., 74., 23., 72., 46., &
    29.,  0., 55., 46., 42., 43., 43., 23., 23., 31., 41., 51., 11., 52., 21., &
    82., 55.,  0., 68., 46., 55., 23., 43., 41., 29., 79., 21., 64., 31., 51., &
    46., 46., 68.,  0., 82., 15., 72., 31., 62., 42., 21., 51., 51., 43., 64., &
    68., 42., 46., 82.,  0., 74., 23., 52., 21., 46., 82., 58., 46., 65., 23., &
    52., 43., 55., 15., 74.,  0., 61., 23., 55., 31., 33., 37., 51., 29., 59., &
    72., 43., 23., 72., 23., 61.,  0., 42., 23., 31., 77., 37., 51., 46., 33., &
    42., 23., 43., 31., 52., 23., 42.,  0., 33., 15., 37., 33., 33., 31., 37., &
    51., 23., 41., 62., 21., 55., 23., 33.,  0., 29., 62., 46., 29., 51., 11., &
    55., 31., 29., 42., 46., 31., 31., 15., 29.,  0., 51., 21., 41., 23., 37., &
    29., 41., 79., 21., 82., 33., 77., 37., 62., 51.,  0., 65., 42., 59., 61., &
    74., 51., 21., 51., 58., 37., 37., 33., 46., 21., 65.,  0., 61., 11., 55., &
    23., 11., 64., 51., 46., 51., 51., 33., 29., 41., 42., 61.,  0., 62., 23., &
    72., 52., 31., 43., 65., 29., 46., 31., 51., 23., 59., 11., 62.,  0., 59., &
    46., 21., 51., 64., 23., 59., 33., 37., 11., 37., 61., 55., 23., 59.,  0. /
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_dist_check_test():'
  write ( *, '(a)' ) '  graph_dist_check() checks a distance matrix.'

  call graph_dist_check ( a, nnode, ierror )

  if ( ierror == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'The distance matrix passed all tests.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) 'The distance matrix failed test ', ierror
  end if

  return
end
subroutine graph_dist_min_span_tree_test ( )

!*****************************************************************************80
!
!! graph_dist_min_span_tree_test() tests graph_dist_min_span_tree().
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

  real ( kind = rk ) dist(nnode,nnode)
  integer i
  integer itree(nnode-1)
  integer j
  integer jtree(nnode-1)
  real ( kind = rk ) wtree(nnode-1)

  data ( ( dist(i,j), i = 1, nnode ), j = 1, nnode ) / &
       0.0D+00, 100.0D+00, 125.0D+00, 120.0D+00, 110.0D+00, &
     100.0D+00,   0.0D+00,  40.0D+00,  65.0D+00,  60.0D+00, &
     125.0D+00,  40.0D+00,   0.0D+00,  45.0D+00,  55.0D+00, &
     120.0D+00,  65.0D+00,  45.0D+00,   0.0D+00,  50.0D+00, &
     110.0D+00,  60.0D+00,  55.0D+00,  50.0D+00,   0.0D+00 /

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_dist_min_span_tree_test():'
  write ( *, '(a)' ) '  graph_dist_min_span_tree() finds a minimum spanning tree.'
  write ( *, '(a)' ) ' '
 
  call graph_dist_print ( dist, nnode, '  The graph:' )

  call graph_dist_min_span_tree ( nnode, dist, itree, jtree )
 
  do i = 1, nnode-1
    wtree(i) = dist(itree(i),jtree(i))
  end do
 
  call graph_arc_weight_print ( nnode-1, itree, jtree, wtree, &
    '  The minimal spanning tree:' )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The length of the minimal tree is ', sum ( wtree )
 
  return
end
subroutine graph_dist_min_span_tree2_test ( )

!*****************************************************************************80
!
!! graph_dist_min_span_tree2_test() tests graph_dist_min_span_tree2().
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

  integer class(nnode)
  real ( kind = rk ) dist(nnode,nnode)
  integer i
  integer itree(nnode-1)
  integer j
  integer jtree(nnode-1)
  real ( kind = rk ) wtree(nnode-1)

  data ( ( dist(i,j), i = 1, nnode ), j = 1, nnode ) / &
       0.0, 100.0, 125.0, 120.0, 110.0, &
     100.0,   0.0,  40.0,  65.0,  60.0, &
     125.0,  40.0,   0.0,  45.0,  55.0, &
     120.0,  65.0,  45.0,   0.0,  50.0, &
     110.0,  60.0,  55.0,  50.0,   0.0D+00 /

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_dist_min_span_tree2_test():'
  write ( *, '(a)' ) '  graph_dist_min_span_tree2() finds a minimum spanning tree.'
  write ( *, '(a)' ) ' '
 
  call graph_dist_print ( dist, nnode, '  The graph:' )

  call graph_dist_min_span_tree2 ( nnode, dist, class, itree, jtree )
 
  do i = 1, nnode-1
    wtree(i) = dist(itree(i),jtree(i))
  end do

  call graph_arc_weight_print ( nnode-1, itree, jtree, wtree, &
    '  The minimal spanning tree:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The length of the minimal tree is ', sum ( wtree )
 
  return
end
subroutine graph_dist_min_span_tree3_test ( )

!*****************************************************************************80
!
!! graph_dist_min_span_tree3_test() tests graph_dist_min_span_tree3().
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

  real ( kind = rk ) dist(nnode,nnode)
  integer i
  integer itree(nnode-1)
  integer jtree(nnode-1)
  integer j
  real ( kind = rk ) wtree(nnode-1)

  data ( ( dist(i,j), i = 1, nnode ), j = 1, nnode ) / &
       0.0, 100.0, 125.0, 120.0, 110.0, &
     100.0,   0.0,  40.0,  65.0,  60.0, &
     125.0,  40.0,   0.0,  45.0,  55.0, &
     120.0,  65.0,  45.0,   0.0,  50.0, &
     110.0,  60.0,  55.0,  50.0,   0.0D+00 /

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_dist_min_span_tree3_test()'
  write ( *, '(a)' ) '  graph_dist_min_span_tree3() finds a minimum spanning tree.'
  write ( *, '(a)' ) ' '
 
  call graph_dist_print ( dist, nnode, '  The graph:' )

  call graph_dist_min_span_tree3 ( nnode, dist, itree, jtree )

  do i = 1, nnode-1
    wtree(i) = dist(itree(i),jtree(i))
  end do

  call graph_arc_weight_print ( nnode-1, itree, jtree, wtree, &
    '  The minimal spanning tree:' )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The length of the minimal tree is ', sum ( wtree )
 
  return
end
subroutine graph_dist_min_span_tree_alternate_test ( )

!*****************************************************************************80
!
!! graph_dist_min_span_tree_alternate_test() tests graph_dist_min_span_tree().
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

  integer, parameter :: nnode = 57

  real ( kind = rk ) dist(nnode,nnode)
  character ( len = 80 ) :: file_name = '57_city_distances.txt'
  integer i
  integer ios
  integer itree(nnode-1)
  integer iunit
  integer jtree(nnode-1)
  real ( kind = rk ) wtree(nnode-1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_dist_min_span_tree_alternate_test()'
  write ( *, '(a)' ) '  graph_dist_min_span_tree() finds a minimum '
  write ( *, '(a)' ) '  spanning tree.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read distance data for 57 cities from file.'
!
!  Read the data.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = file_name, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Problems opening the file: ' // trim ( file_name )
    write ( *, '(a)' ) '  The test was abandoned.'
    return
  end if

  do i = 1, nnode

    read ( iunit, *, iostat = ios ) dist(i,1:nnode)

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Problems reading the data.'
      write ( *, '(a)' ) '  The test was abandoned.'
      return
    end if

  end do

  close ( unit = iunit )
!
!  Compute the tree.
!
  call graph_dist_min_span_tree ( nnode, dist, itree, jtree )
 
  do i = 1, nnode-1
    wtree(i) = dist(itree(i),jtree(i))
  end do
 
  call graph_arc_weight_print ( nnode-1, itree, jtree, wtree, &
    '  The weighted tree:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The length of the minimal tree is ', sum ( wtree )

  return
end
subroutine graph_dist_one_test ( )

!*****************************************************************************80
!
!! graph_dist_one_test() tests graph_dist_one().
!
!  Discussion:
!
!    This example appears on page 15 of the reference book by Gibbons.
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

  real ( kind = rk ) dinfin
  real ( kind = rk ) dist(nnode,nnode)
  integer i
  integer idad(nnode)
  integer inode
  integer path(nnode)
  integer itemp(nnode)
  integer j
  integer length
  real ( kind = rk ) path_dist(nnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_dist_one_test():'
  write ( *, '(a)' ) '  graph_dist_one() computes the distance from one'
  write ( *, '(a)' ) '  node to all others in a graph.'
  write ( *, '(a)' ) ' '

  dinfin = 1000.0D+00
 
  do i = 1, nnode
    do j = 1, nnode
      dist(i,j) = dinfin
    end do
    dist(i,i) = 0.0D+00
  end do
 
  dist(1,2) = 1.0D+00
  dist(1,3) = 3.0D+00
 
  dist(2,1) = 2.0D+00
  dist(2,3) = 1.0D+00
  dist(2,5) = 2.0D+00
 
  dist(3,4) = 2.0D+00
  dist(3,5) = 3.0D+00
 
  dist(4,3) = 1.0D+00
 
  dist(5,1) = 1.0D+00
  dist(5,2) = 3.0D+00
  dist(5,4) = 6.0D+00

  call graph_dist_print ( dist, nnode, '  Edge Distance Matrix:' )

  inode = 5
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) 'The starting node is ', inode
  write ( *, '(a)' ) ' '

  call graph_dist_one ( dist, dinfin, path_dist, idad, inode, path, nnode )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Node    Distance   Path Idad'
  write ( *, '(a)' ) ' '
 
  do i = 1, nnode
    write ( *, '(i5,g14.6,2i5)' ) i, path_dist(i), path(i), idad(i)
  end do
 
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Note that "infinity" is represented by ', dinfin
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here are the paths for each node:'
  write ( *, '(a)' ) ' '
 
  do i = 1, nnode

    length = 1
    itemp(length) = i
 
    do while ( itemp(length) /= inode )
      length = length+1
      itemp(length) = idad(itemp(length-1))
    end do
 
    write ( *, '(5i5)' ) itemp(1:length)
 
  end do
 
  return
end
subroutine graph_dist_pairing_greedy_test ( )

!*****************************************************************************80
!
!! graph_dist_pairing_greedy_test() tests graph_dist_pairing_greedy().
!
!  Discussion:
!
!    Random data is used in setting up the problem.
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

  integer, parameter :: nnode = 15
  integer, parameter :: nnodeh = 7

  real ( kind = rk ) dist
  integer ido
  integer indx
  integer maxit
  integer nodeb(nnode)
  integer nodeb1
  integer noder(nnode)
  integer noder1
  real ( kind = rk ) tol
  real ( kind = rk ) total
  real ( kind = rk ) xb(nnode)
  real ( kind = rk ) xhi
  real ( kind = rk ) xlo
  real ( kind = rk ) xr(nnode)
  real ( kind = rk ) yb(nnode)
  real ( kind = rk ) yhi
  real ( kind = rk ) ylo
  real ( kind = rk ) yr(nnode)
!
!  IDO just tells us if this is the first or later trials.
!
  ido = 1
!
!  Set the maximum number of iterations.
!
  maxit = 10
!
!  Set the range of the X and Y coordinates.
!
  xhi = 10.0D+00
  xlo = 0.0D+00
  yhi = 5.0D+00
  ylo = 3.0D+00
!
!  Set the relative tolerance for the stepwise distance decrease.
!
  tol = 0.05D+00
!
!  Say hello.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_dist_pairing_greedy_test():'
  write ( *, '(a)' ) '  graph_dist_pairing_greedy() tries to minimize the total distance'
  write ( *, '(a)' ) '  in a pairing of black and red nodes.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Try to find a pairing of two sets of nodes'
  write ( *, '(a)' ) '  with a low discrepancy.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Relative tolerance for step decrease = ', tol
  write ( *, '(a,i8)' ) '  Maximum number of steps = ', maxit
  write ( *, '(a,g14.6,a,g14.6)' ) '  X range is ', xlo,' to ', xhi
  write ( *, '(a,g14.6,a,g14.6)' ) '  Y range is ', ylo,' to ', yhi
!
!  Make an arbitrary pairing of the nodes.
!
  do indx = 1, nnode
    nodeb(indx) = indx
    noder(indx) = indx
  end do
!
!  Make up a random set of X, Y coordinates for the nodes.
!
  call r8vec_uniform_ab ( nnode, xlo, xhi, xb )
  call r8vec_uniform_ab ( nnode, xlo, xhi, xr )
  call r8vec_uniform_ab ( nnode, ylo, yhi, yb )
  call r8vec_uniform_ab ( nnode, ylo, yhi, yr )
!
!  We will jump back here if we restart with a permuted NODER.
!
  do ido = 1, 2
 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Initial black node coordinates:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '    I   Black   X             Y'
    write ( *, '(a)' ) ' '
 
    do indx = 1, nnode
      write ( *, '(2i8,2g14.6)' ) indx, nodeb(indx), xb(indx), yb(indx)
    end do
 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Initial red node coordinates:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '    I    Red    X             Y'
    write ( *, '(a)' ) ' '
 
    do indx = 1, nnode
      write ( *, '(2i8,2g14.6)' ) indx, noder(indx), xr(indx), yr(indx)
    end do
 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Initial pairing of nodes:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '    I   Black  Red    Distance'
    write ( *, '(a)' ) ' '
 
    do indx = 1, nnode
      nodeb1 = nodeb(indx)
      noder1 = noder(indx)
      dist = sqrt ( ( xb(nodeb1) - xr(noder1) )**2 + &
                    ( yb(nodeb1) - yr(noder1) )**2 )

      write ( *, '(3i8,g14.6)' ) indx, nodeb1, noder1, dist
    end do
 
    total = 0.0D+00
    do indx = 1, nnode
      nodeb1 = nodeb(indx)
      noder1 = noder(indx)
      total = total + sqrt ( ( xb(nodeb1) - xr(noder1) )**2 &
                           + ( yb(nodeb1) - yr(noder1) )**2 )
    end do
 
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) 'Total discrepancy of initial pairing = ', total
!
!  Seek a better pairing.
!
    call graph_dist_pairing_greedy ( maxit, nodeb, noder, nnode, tol, xb, xr, yb, yr )
 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Final black node coordinates:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '    I   Black   X             Y'
    write ( *, '(a)' ) ' '
 
    do indx = 1, nnode
      write ( *, '(2i8,2g14.6)' ) indx, nodeb(indx), xb(indx), yb(indx)
    end do
 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Final red node coordinates:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '    I    Red    X             Y'
    write ( *, '(a)' ) ' '
 
    do indx = 1, nnode
      write ( *, '(2i8,2g14.6)' ) indx, noder(indx), xr(indx), yr(indx)
    end do
 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Final pairing of nodes:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '    I   Black  Red    Distance'
    write ( *, '(a)' ) ' '
 
    do indx = 1, nnode

      nodeb1 = nodeb(indx)
      noder1 = noder(indx)

      dist = sqrt ( ( xb(nodeb1) - xr(noder1) )**2 &
                  + ( yb(nodeb1) - yr(noder1) )**2 )

      write ( *, '(3i8,g14.6)') indx, nodeb1, noder1, dist

    end do
 
    total = 0.0D+00
    do indx = 1, nnode
      nodeb1 = nodeb(indx)
      noder1 = noder(indx)
      dist = sqrt ( ( xb(nodeb1) - xr(noder1) )**2 &
                  + ( yb(nodeb1) - yr(noder1) )**2 )

      total = total + dist
 
    end do
 
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Total discrepancy of final pairing = ', total
!
!  On the second try, reverse the ordering of the red nodes.
!  Any random permutation would be worth trying.
!
    if ( ido == 1 ) then
 
      do indx = 1, nnodeh
        call i4_swap ( noder(indx), noder(nnode+1-indx) )
      end do
 
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Reversing NODER!'
 
    end if

  end do
 
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

