subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! get_unit() returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
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
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IUNIT.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5 and 6).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
  implicit none

  integer i
  integer ios
  integer iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do
!
!  No free unit was found.
!
  iunit = 0

  return
end
subroutine graph_arc_weight_print ( nedge, inode, jnode, wnode, title )

!*****************************************************************************80
!
!! graph_arc_weight_print() prints out a weighted graph from an edge list.
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
!  Input:
!
!    integer NEDGE, the number of edges.
!
!    integer INODE(NEDGE), JNODE(NEDGE), the beginning 
!    and end nodes of the edges.
!
!    real ( kind = rk ) WNODE(NEDGE), the weights of the edges.
!
!    character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nedge

  integer i
  integer inode(nedge)
  integer jnode(nedge)
  character ( len = * ) title
  real ( kind = rk ) wnode(nedge)

  if ( len_trim ( title ) /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  do i = 1, nedge
    write ( *, '(i8,4x,2i8,g14.6)' ) i, inode(i), jnode(i), wnode(i)
  end do

  return
end
subroutine graph_dist_all ( dist, dinfin, nnode, path_dist )

!*****************************************************************************80
!
!! graph_dist_all() finds the distance from every node to every other node.
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
!    Input/output, real ( kind = rk ) DIST(NNODE,NNODE). 
!
!    On input, DIST(I,J) is the length of the edge FROM node I TO node J.
!    DIST(I,J) = DINFIN if there is no direct edge from I to J.
!
!    On output, DIST has been overwritten by other information.
!
!    Input, real ( kind = rk ) DINFIN, is a "large" number, larger than any 
!    entry in DIST, which is taken to be "infinity".  DIST(I,J) = DINFIN 
!    means there is no direct edge from node I to node J.  On output,
!    DIST(I,J) = DINFIN means there is no path from node I to node J.
!
!    Input, integer NNODE, the number of nodes.
!
!    Output, real ( kind = rk ) PATH_DIST(NNODE,NNODE).  This array contains the
!    lengths of the shortest paths from each node to another node.
!    PATH_DIST(I,J) is the length of the shortest path from node I
!    to node J.  If PATH_DIST(I,J) = DINFIN, there is no path from node
!    I to node J.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nnode

  real ( kind = rk ) dist(nnode,nnode)
  real ( kind = rk ) dinfin
  real ( kind = rk ) path_dist(nnode,nnode)
  integer i
  integer j
  integer k

  do k = 1, nnode
 
    do i = 1, nnode
      do j = 1, nnode
 
        path_dist(i,j) = dist(i,j)

        if ( dist(i,k) /= dinfin .and. dist(k,j) /= dinfin ) then

          path_dist(i,j) = min ( path_dist(i,j), dist(i,k) + dist(k,j) )

        end if
 
      end do
    end do
 
    dist(1:nnode,1:nnode) = path_dist(1:nnode,1:nnode)
 
  end do
 
  return
end
subroutine graph_dist_check ( dist, nnode, ierror )

!*****************************************************************************80
!
!! graph_dist_check() checks a distance matrix for consistency.
!
!  Discussion:
!
!    The checks made are:
!
!      1): DIST(I,I) = 0
!      2): DIST(I,J) > 0 for I different from J
!      3): DIST(I,J) = DIST(J,I) for I different from J.
!      4): DIST(I,J) + DIST(J,K) >= DIST(I,K).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) DIST(NNODE,NNODE), the distance matrix.  
!    DIST(I,J) is the distance FROM node I TO node J.
!
!    Input, integer NNODE, the number of nodes.
!
!    Output, integer IERROR, error flag.
!    0, no errors.
!    1, DIST(I,I) is nonzero for some I;
!    2, DIST(I,J) <= 0 for some distinct I, J
!    3, DIST(I,J) not equal to DIST(J,I) for some distinct I, J.
!    4, DIST(I,J) + DIST(J,K) < DIST(I,K) for some I, J, K.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nnode

  real ( kind = rk ) dist(nnode,nnode)
  integer i
  integer ierror
  integer j
  integer k

  ierror = 1

  do i = 1, nnode
    if ( dist(i,i) /= 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GRAPH_DIST_CHECK - Failed test #1:'
      write ( *, '(a,i8)' ) '  DIST(I,I) nonzero for I = ', i
      return
    end if
  end do

  ierror = 2
  do i = 1, nnode
    do j = 1, nnode
      if ( i /= j ) then
        if ( dist(i,j) <= 0.0D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'GRAPH_DIST_CHECK - Failed test #2:'
          write ( *, '(a,2i8)' ) '  DIST(I,J) <= 0 for I, J = ', i, j
          return
        end if
      end if
    end do
  end do

  ierror = 3
  do i = 1, nnode
    do j = 1, i - 1
      if ( dist(i,j) /= dist(j,i) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRAPH_DIST_CHECK - Failed test #3:'
        write ( *, '(a)' ) '  DIST(I,J) is not equal to DIST(J,I)'
        write ( *, '(a,2i8)' ) '  for I, J = ', i, j
        return
      end if  
    end do
  end do

  ierror = 4
  do i = 1, nnode
    do j = 1, nnode
      do k = 1, i - 1
        if ( dist(i,j) + dist(j,k) < dist(i,k) ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'GRAPH_DIST_CHECK - Failed test #4:'
          write ( *, '(a)' ) '  DIST(I,J) + DIST(J,K) < DIST(I,K)'
          write ( *, '(a,3i8)' ) '  I, J, K = ', i, j, k
          write ( *, '(a,g14.6)' ) '  DIST(I,J) = ', dist(i,j)
          write ( *, '(a,g14.6)' ) '  DIST(J,K) = ', dist(j,k)
          write ( *, '(a,g14.6)' ) '  DIST(I,K) = ', dist(i,k)
          return
        end if
      end do
    end do
  end do

  ierror = 0

  return
end
subroutine graph_dist_min_span_tree ( nnode, dist, itree, jtree )

!*****************************************************************************80
!
!! graph_dist_min_span_tree() computes a spanning tree of minimal length.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 July 2000
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
!    Input, real ( kind = rk ) DIST(NNODE,NNODE).  DIST(I,J) = distance from node I
!    to node J.
!
!    Output, integer ITREE(NNODE-1), JTREE(NNODE-1), the pairs of 
!    nodes that form the edges of the tree.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nnode

  real ( kind = rk ) d
  real ( kind = rk ) dist(nnode,nnode)
  real ( kind = rk ) dmin
  integer i
  integer imin
  integer it
  integer itree(nnode-1)
  integer j
  integer jtree(nnode-1)

  call i4vec_indicator ( nnode-1, itree )

  jtree(1:nnode-1) = -nnode
 
  do j = 1, nnode-1
!
!  Choose the node IMIN whose tree edge ( ITREE(IMIN)=IMIN, JTREE(IMIN) ) 
!  will be set.
!
    dmin = huge ( dmin )
 
    do i = 1, nnode-1
 
      it = jtree(i)
 
      if ( it < 0 ) then
 
        d = dist(-it,i)
 
        if ( d < dmin ) then
          dmin = d
          imin = i
        end if
 
      end if
 
    end do
 
    jtree(imin) = - jtree(imin)

    do i = 1, nnode-1
 
      it = jtree(i)
 
      if ( it < 0 ) then
        if ( dist(i,imin) < dist(i,-it) ) then
          jtree(i) = - imin
        end if
      end if
 
    end do
 
  end do
 
  return
end
subroutine graph_dist_min_span_tree2 ( nnode, dist, class, itree, jtree )

!*****************************************************************************80
!
!! graph_dist_min_span_tree2() computes a spanning tree of minimal length.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, real ( kind = rk ) DIST(NNODE,NNODE).  DIST(I,J) = distance from node I
!    to node J.
!
!    Output, integer CLASS(NNODE), lists the nodes in the order in
!    which they joined the tree.
!
!    Output, integer ITREE(NNODE-1), JTREE(NNODE-1), the pairs of 
!    nodes that form the edges of the tree.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nnode

  integer class(nnode)
  real ( kind = rk ) dist(nnode,nnode)
  real ( kind = rk ) dmin
  integer i
  integer ii
  integer imin
  integer itree(nnode-1)
  integer j
  integer jj
  integer jjmin
  integer jmin
  integer jtree(nnode-1)
  integer k
  integer npos
  logical smaller
  logical unset

  if ( nnode <= 1 ) then
    return
  end if
!
!  All the nodes start out in the negative class.
!
  npos = 0
  call i4vec_indicator ( nnode, class )
!
!  Find the shortest edge (I,J).
!
  unset = .true.
  dmin = 0.0D+00

  do i = 1, nnode
    do j = i + 1, nnode

      if ( unset ) then
        smaller = .true.
      else if ( dist(i,j) < dmin ) then
        smaller = .true.
      else
        smaller = .false.
      end if

      if ( smaller ) then
        imin = i
        jmin = j
        dmin = dist(i,j)
        unset = .false.
      end if

    end do
  end do
!
!  Carry nodes IMIN and JMIN into the positive class.
!
  npos = npos + 1
  call i4_swap ( class(npos), class(imin) )

  npos = npos + 1
  call i4_swap ( class(npos), class(jmin) )

  itree(1) = imin
  jtree(1) = jmin
!
!  Now, repeatedly, find the shortest edge connecting a negative
!  and positive node.  Move the negative node to the positive class and
!  repeat.
!
  do k = 2, nnode - 1

    unset = .true.
    dmin = 0.0D+00
    imin = - 99
    jmin = - 99

    do ii = 1, npos

      i = class(ii)

      do jj = npos + 1, nnode

        j = class(jj)

        if ( unset ) then
          smaller = .true.
        else if ( dist(i,j) < dmin ) then
          smaller = .true.
        else
          smaller = .false.
        end if

        if ( smaller ) then
          imin = i
          jmin = j
          jjmin = jj
          dmin = dist(i,j)
          unset = .false.
        end if

      end do

    end do

    npos = npos + 1
    call i4_swap ( class(npos), class(jjmin) )

    itree(k) = imin
    jtree(k) = jmin

  end do

  return
end
subroutine graph_dist_min_span_tree3 ( nnode, dist, inode, jnode )

!*****************************************************************************80
!
!! graph_dist_min_span_tree3() finds a minimum spanning tree.
!
!  Discussion:
!
!    The input graph is represented by a distance matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 July 2000
!
!  Reference:
!
!    Hang Tong Lau,
!    Combinatorial Heuristic Algorithms in FORTRAN,
!    Springer Verlag, 1986.
!
!  Parameters:
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, real ( kind = rk ) DIST(NNODE,NNODE), an NNODE by NNODE distance
!    matrix.  DIST(I,J) is the distance from node I to node J.  The matrix
!    should be symmetric.  If there is no arc from node I to node J,
!    set DIST(I,J) = HUGE(1.0).
!
!    Output, integer INODE(NNODE), JNODE(NNODE); entries 1 through 
!    NNODE-1 describe the edges of the spanning tree as pairs of nodes.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nnode

  real ( kind = rk ) d
  real ( kind = rk ) dist(nnode,nnode)
  integer i
  integer ient
  integer ij
  integer inode(nnode)
  integer itr
  integer iwork1(nnode)
  integer iwork2(nnode)
  integer j
  integer jj
  integer jnode(nnode)
  integer kj
  integer kk4
  real ( kind = rk ) tree_length
  real ( kind = rk ) work(nnode)

  work(1:nnode) = huge ( work(1) )
  iwork1(1:nnode) = 0
  iwork2(1:nnode) = 0
!
!  Find the first non-zero arc.
!
  do ij = 1, nnode
    do kj = 1, nnode
      if ( dist(ij,kj) < huge ( dist(1,1) ) ) then
        i = ij
        go to 10
      end if
    end do
  end do

10    continue

  work(i) = 0
  iwork1(i) = 1
  tree_length = 0.0D+00
  kk4 = nnode - 1

  do jj = 1, kk4

    work(1:nnode) = huge ( work(1) )

    do i = 1, nnode
!
!  For each forward arc originating at node I calculate
!  the length of the path to node I
!
      if ( iwork1(i) == 1 ) then
        do j = 1, nnode
          if ( dist(i,j) < huge ( dist(1,1) ) .and. iwork1(j) == 0 ) then
            d = tree_length + dist(i,j)
            if ( d < work(j) ) then
              work(j) = d
              iwork2(j) = i
            end if
          end if
        end do
      end if

    end do
!
!  Find the minimum potential.
!
    d = huge ( d )
    ient = 0

    do i = 1, nnode
      if ( iwork1(i) == 0 .and. work(i) < d ) then
        d = work(i)
        ient = i
        itr = iwork2(i)
      end if
    end do
!
!  Include the node in the current path.
!
    if ( d < huge ( d ) ) then
      iwork1(ient) = 1
      tree_length = tree_length + dist(itr,ient)
      inode(jj) = itr
      jnode(jj) = ient
    end if

  end do

  return
end
subroutine graph_dist_one ( dist, dinfin, path_dist, dad, inode, path, nnode ) 

!*****************************************************************************80
!
!! graph_dist_one() computes the distance from one node to all others in a graph.
!
!  Discussion:
!
!    This routine can handle both ordinary graphs and directed graphs.  
!
!    In an ordinary graph, a connection between two nodes is always guaranteed
!    to be "symmetric".  That is, if node I is connected to node J by
!    an edge of length D, then node J is connected to node I, and the
!    distance is again D.
!
!    In a directed graph, if node I is connect to node J by an edge of
!    length D, then nothing is known about a possible connection from 
!    node J back to node I.  In particular, it is possible that:
!
!    * there is no direct edge from node J to node I;
!    * the edge from node J to node I exists, but is a different "length"
!      than the edge from node I to node J.
!
!    The program computes:
!
!    * PATH_DIST, an array of distances from node INODE to all other nodes;
!
!    * DAD, an array which can be used to determine the path from
!      node INODE to any particular node.
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
!    Input, real ( kind = rk ) DIST(NNODE,NNODE).  DIST contains the weighted
!    adjacency information defining the graph, or directed graph.
!    The diagonal entries of DIST, that is, DIST(I,I), should be set to 0.
!    The value of the typical off diagonal element DIST(I,J) should 
!    represent the length, or weight, of the edge from node I to
!    node J.  If the graph is undirected, then DIST(I,J) should always
!    equal DIST(J,I).  For a directed graph, these quantities may differ.
!    If there is no edge from node I to node J, then it would be natural
!    to set DIST(I,J) to "infinity".  Since this is not computationally
!    possible, the user must specify a special value, called DINFIN,
!    that will be used to mark such entries.  The most natural thing
!    to do would simply be to pick DINFIN to be "very large", such
!    as DINFIN = 10,000.
!    All the entries in DIST should be non-negative.  The algorithm will
!    NOT work correctly if negative edge lengths are input.
!    Off-diagonal elements DIST(I,J) may be set to zero.  This simply
!    means that two nodes are "very close", like St Paul and Minneapolis.
!
!    Input, real ( kind = rk ) DINFIN, is a "large" number, which should be
!    larger than the length of any edge in the graph, and in fact larger
!    than the length of any reasonable path along the edges of the graph.  
!    The user should have set the DIST matrix so that DIST(I,J) = DINFIN
!    whenever there is no edge from node I to node J.  The program has to 
!    know the value of DINFIN so it can understand this information stored
!    in DIST.
!
!    Output, real ( kind = rk ) PATH_DIST(NNODE).  On output, for every value
!    of I from 1 to NNODE, PATH_DIST(I) contains the distance from node INODE 
!    to node I in the graph.  Of course, PATH_DIST(INODE) is zero.  Moreover,
!    if PATH_DIST(I) = DINFIN, then this is the program's way of reporting that
!    there is NO path from node INODE to node I.
!
!    Output, integer DAD(NNODE), information defining the shortest 
!    path from node INODE to any node I, which presumably will be of
!    total distance PATH_DIST(I).
!
!    The path from node I to node INODE, is recorded "in reverse"
!    in DAD.  The last node is INODE, of course.  The previous node
!    is DAD(INODE).  The next node is DAD(DAD(INODE)) and
!    so on, until INODE itself is reached.  
!
!    If the distance from node I to node INODE is "infinity", then
!    DAD will still record a path; it's just probably of no interest.
!
!    Input, integer INODE, the base node, from which distances to 
!    the other nodes are to be calculated.
!
!    Output, integer PATH(NNODE).  The value of PATH(I) records
!    the step on which the distance from INODE to node I was
!    determined.  There will be NNODE steps, and on each step
!    just one such distance is computed.
!
!    Input, integer NNODE, the number of nodes.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nnode

  integer dad(nnode)
  real ( kind = rk ) dist(nnode,nnode)
  real ( kind = rk ) dinfin
  real ( kind = rk ) dmin
  real ( kind = rk ) dtemp
  integer imin
  integer inode
  integer istep
  integer j
  integer path(nnode)
  real ( kind = rk ) path_dist(nnode)
!
!  Initialize the data.
!
  dad(1:nnode) = inode
  path(1:nnode) = 0
  path_dist(1:nnode) = dist(inode,1:nnode)
!
!  On step 1, we connect node INODE itself.
!
  dad(inode) = inode
  path(inode) = 1
!
!  On steps ISTEP = 2 through NNODE, we try to add just one more node.
!
!  Of all the nodes which are not yet connected to INODE (because PATH
!  is 0 for this node), choose the one whose distance is least.
!
  do istep = 2, nnode
 
    dmin = dinfin
    imin = 0
 
    do j = 1, nnode
 
      if ( path(j) == 0 ) then
        if ( path_dist(j) <= dmin ) then
          dmin = path_dist(j)
          imin = j
        end if
      end if
 
    end do
!
!  If we found no new node to add, then any remaining nodes cannot
!  be connected.
!
    if ( dmin == dinfin ) then
      return
    end if
!
!  Now add the closest node, labeled IMIN, to the list.
!
    path(imin) = istep
!
!  Update the distances of the remaining unconnected nodes.
!
    do j = 1, nnode
 
      if ( path(j) == 0 ) then
 
        dtemp = path_dist(imin) + dist(imin,j)
 
        if ( dtemp < path_dist(j) ) then
          path_dist(j) = dtemp
          dad(j) = imin
        end if
 
      end if
 
    end do
 
  end do
 
  return
end
subroutine graph_dist_pairing_greedy ( maxit, nodeb, noder, nnode, tol, xb, &
  xr, yb, yr )

!*****************************************************************************80
!
!! graph_dist_pairing_greedy() pairs two sets of nodes using the least total distance.
!
!  Discussion:
!
!    The method is iterative, and is not guaranteed to find the best
!    possible arrangement.  This is particulary true because it is a
!    "local" method, which only considers pairwise switches of the
!    red nodes that reduce the total distance.  This means that a
!    "locally minimizing" pairing might be found which is not the
!    global minimizer.
!
!    On the other hand, in the absence of a theoretical plan for how
!    to reach the global minimizer, the brute force search would
!    require that EVERY possible pairing be considered, and its total
!    distance computed.  This means that a total of NNODE!
!    graphs would have to be generated.
!
!    The approach used here, on each iterative step, looks at a
!    maximum of NNODE * (NNODE-1) graphs, which represents a
!    significantly more efficient method.
!
!    It would not be hard to extend this approach to a method which
!    considers switches of THREE red nodes at a time, though the
!    work there involve looking at NNODE * (NNODE-1) * (NNODE-2)
!    graphs, and as we increase the number of graphs we examine,
!    we begin to approach the NNODE! rate for the brute force
!    algorithm.
!
!    It also would not be hard to extend this method to a case where
!    there are three sets of nodes, arranged in triples, and again
!    the total distance is to be minimized.
!
!
!    If it is suspected that the pairing returned is only
!    a local minimizer, then the user is advised to restart the
!    calculation after randomly permuting the entries of NODER, so that
!    the routine starts from a different point in the space of graphs.
!
!    The routine is given:
!
!      an initial ordering of the black and red nodes, so that
!      ( NODEB(I), NODER(I) ) represents the I-th pair,
!
!      the X and Y coordinates of the black and red nodes,
!
!      a maximum number of iterations, and a relative distance
!      decrease requirement,
!
!    and computes:
!
!      a new ordering of the red nodes, contained in NODER, which should
!      reduce the total distance between corresponding red and black
!      nodes.
!
!
!    The approach can be applied to a variety of problems including:
!
!    1) We are given two sets of NNODE points, which we will call the
!       "red" and "black " groups, and the (X,Y) coordinates of each
!       point.  We may imagine these points as forming the two sets of
!       nodes of a bipartite graph lying in the (X,Y) plane.  We wish
!       to choose a pairing of red and black nodes which results in
!       the shortest total arc length.
!
!    2) We are given two sets of NNODE complex quantities, which we
!       believe are approximations to the same (unknown) set of
!       quantities.  We wish to arrange this data into NNODE pairs,
!       each containing a unique element from each set of data, which
!       minimizes the sum of squares of the discrepancies between the
!       pairs.  In particular, the two sets of data might be two
!       separate estimates of the complex eigenvalues of a matrix.
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
!    Input, integer MAXIT, the maximum number of iterations
!    allowed.  Each iteration considers, one at a time, each black node, and
!    seeks to switch its red neighbor for another red neighbor that
!    reduces the total distance.
!
!    Input, integer NODEB(NNODE), the "labels" of the black
!    nodes.  You probably want to just set NODEB(I) = I, for i = 1 to NNODE.  
!    The entries in NODEB will not be changed.
!
!    Input/output, integer NODER(NNODE), the "labels" of the red
!    nodes.  You probably want to just set the input value of NODER(I) = I, 
!    for i = 1 to NNODE.  The entries in NODER WILL be changed.
!
!    At all times, the values of ( NODEB(I), NODER(I) ) contain the
!    labels of the I-th pair of black and red nodes.
!
!    On output, if a better pairing of the nodes has been found,
!    this will be reflected in the newly permuted values of NODER.
!
!    Input, integer NNODE, the number of nodes in the black, 
!    and in the red sets.
!
!    Input, real ( kind = rk ) TOL.
!    TOL is the relative decrease that the user demands in the
!    total distance, after each iterative step.  If we denote
!    the distance before the iterative step as OLDTOT, and the
!    distance after the iterative step as TOTAL, then the
!    routine will try another iterative step as long as "enough"
!    progress was made on this step.  Enough progress was made
!    whenever OLDTOT - TOTAL < TOL * TOTAL
!
!    Input, real ( kind = rk ) XB(NNODE), the X coordinates of the black nodes.
!
!    Input, real ( kind = rk ) XR(NNODE), the X coordinates of the red nodes.
!
!    Input, real ( kind = rk ) YB(NNODE), the Y coordinates of the black nodes.
!
!    Input, real ( kind = rk ) YR(NNODE), the Y coordinates of the red nodes.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nnode

  real ( kind = rk ) dist1
  real ( kind = rk ) dist2
  integer indx
  integer indx1
  integer indx2
  integer it
  integer maxit
  integer nodeb(nnode)
  integer nodeb1
  integer nodeb2
  integer noder(nnode)
  integer noder1
  integer noder2
  integer nswap
  real ( kind = rk ) oldtot
  real ( kind = rk ) temp
  real ( kind = rk ) tol
  real ( kind = rk ) total
  real ( kind = rk ) xb(nnode)
  real ( kind = rk ) xr(nnode)
  real ( kind = rk ) yb(nnode)
  real ( kind = rk ) yr(nnode)
!
!  Compute the total distance of the starting pairing.
!
  total = 0.0D+00
  do indx = 1, nnode

    nodeb1 = nodeb(indx)
    noder1 = noder(indx)

    total = total + sqrt ( &
      ( xb(nodeb1) - xr(noder1) )**2 + ( yb(nodeb1) - yr(noder1) )**2 )

  end do
 
  write ( *, '(a)' ) ' '
!
!  Begin the iterations.
!
  do it = 1, maxit

    if ( total == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'graph_dist_pairing_greedy(): Early termination.'
      write ( *, '(a)' ) '  Total discrepancy is low enough.'
      return
    end if
!
!  Save the current total distance for comparison at the end of the 
!  iteration.
!
    oldtot = total
    nswap = 0
!
!  Consider each black node, by running through indices INDX1 = 1
!  through NNODE of the NODEB array.
!
    do indx1 = 1, nnode
!
!  Get the actual labels of the current INDX1-th pair of black and 
!  red nodes.
!
      nodeb1 = nodeb(indx1)
      noder1 = noder(indx1)
!
!  Now look at the black node with INDX2 = 1 through NNODE, but ignore
!  the case where INDX1 = INDX2.
!
      do indx2 = 1, nnode
!
!  Get the labels of the current INDX2-th pair of black and red nodes.
!
        nodeb2 = nodeb(indx2)
        noder2 = noder(indx2)
 
        if ( indx2 /= indx1 ) then
!
!  Compute the total distance between (NODEB1,NODER1) and 
!  (NODEB2,NODER2), and compare it to the total where we switch the 
!  red nodes.
!
          dist1 = sqrt ( ( xb(nodeb1) - xr(noder1) )**2 &
                       + ( yb(nodeb1) - yr(noder1) )**2 ) &
                + sqrt ( ( xb(nodeb2) - xr(noder2) )**2 &
                       + ( yb(nodeb2) - yr(noder2) )**2 )

          dist2 = sqrt ( ( xb(nodeb1) - xr(noder2) )**2 &
                       + ( yb(nodeb1) - yr(noder2) )**2 ) &
                + sqrt ( ( xb(nodeb2) - xr(noder1) )**2 &
                       + ( yb(nodeb2) - yr(noder1) )**2 )
! 
!  If the new arrangement is any shorter, take it, by shuffling the
!  red nodes only, and update the total distance.
!
          if ( dist2 < dist1 ) then
            call i4_swap ( noder(indx1), noder(indx2) )
            nswap = nswap + 1
          end if
 
        end if
 
      end do
 
    end do
!
!  Now that we've checked all pairs of nodes,
!  print the new total distance, and see if we may
!  continue, or should give up.
!
    total = 0.0D+00
    do indx1 = 1, nnode

      nodeb1 = nodeb(indx1)
      noder1 = noder(indx1)

      total = total + sqrt ( ( xb(nodeb1) - xr(noder1) )**2 &
                           + ( yb(nodeb1) - yr(noder1) )**2 )

    end do

    write ( *, '(a,i8)' ) '  On step ', it
    write ( *, '(a,g14.6)' ) '  discrepancy =', total
    write ( *, '(a,i8)' ) '  Swaps made was ', nswap
 
    if ( oldtot - total <= tol * oldtot ) then
 
      temp = ( oldtot - total ) / oldtot
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'graph_dist_pairing_greedy(): Warning:'
      write ( *, '(a)' ) '  The relative change in the discrepancy '
      write ( *, '(a,g14.6)' ) '  was only ', temp
      write ( *, '(a,g14.6)' ) '  which is less than the tolerance TOL =',tol
      write ( *, '(a)' ) '  Bailing out of the iteration.'
      write ( *, '(a)' ) ' '
      return
 
    end if
 
  end do
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'graph_dist_pairing_greedy(): Note:'
  write ( *, '(a)' ) '  The discrepancy has decreased by at least the'
  write ( *, '(a)' ) '  tolerance on every step.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Increasing the number of iterations might '
  write ( *, '(a)' ) '  provide further improvement at this rate.'
 
  return
end
subroutine graph_dist_print ( dist, nnode, title )

!*****************************************************************************80
!
!! graph_dist_print() prints a distance matrix.
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
!    Input, real ( kind = rk ) DIST(NNODE,NNODE), the distance matrix.  
!    DIST(I,J) is the distance from node I to node J.
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nnode

  real ( kind = rk ) dist(nnode,nnode)
  integer ihi
  integer ilo
  integer jhi
  integer jlo
  integer ncol
  integer nrow
  character ( len = * ) title

  if ( len_trim ( title ) /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  ilo = 1
  ihi = nnode
  jlo = 1
  jhi = nnode
  ncol = nnode
  nrow = nnode

  call r8mat_print ( dist, ihi, ilo, jhi, jlo, ncol, nrow )

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
subroutine r8mat_print ( a, ihi, ilo, jhi, jlo, ncol, nrow )

!*****************************************************************************80
!
!! r8mat_print() prints out a portion of a dense matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A(NROW,NCOL), an NROW by NCOL matrix to be printed.
!
!    Input, integer IHI, ILO, the first and last rows to print.
!
!    Input, integer JHI, JLO, the first and last columns to print.
!
!    Input, integer NCOL, NROW, the number of rows and columns
!    in the matrix A.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: incx = 5

  integer ncol

  real ( kind = rk ) a(nrow,ncol)
  character ctemp(incx)*14
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
  integer nrow

  write ( *, '(a)' ) ' '

  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, ncol )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(''Columns '',5a14)' ) ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, nrow )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == 0.0D+00 ) then
          ctemp(j2) = '    0.0'
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ctemp(1:inc)

    end do

  end do

  return
end
subroutine r8vec_uniform_ab ( n, a, b, r )

!*****************************************************************************80
!
!! r8vec_uniform_ab() returns a scaled pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    Each dimension ranges from A to B.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    31 May 2007
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
!    integer N, the number of entries in the vector.
!
!    real ( kind = rk ) A, B, the lower and upper limits.
!
!  Output:
!
!    real ( kind = rk ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) r(n)

  call random_number ( harvest = r(1:n) )
  r(1:n) = a + ( b - a ) * r(1:n)

  return
end

