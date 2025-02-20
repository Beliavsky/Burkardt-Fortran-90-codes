subroutine graph_arc_euler_circ ( nnode, nedge, inode, jnode, loop )

!*****************************************************************************80
!
!! graph_arc_euler_circ() finds an Euler circuit in an Eulerian graph.
!
!  Discussion:
!
!    An Euler circuit of a graph is a path that uses each edge exactly once.
!
!    A graph is Eulerian if it has an Euler circuit.
!
!    An Eulerian graph may have many circuits; this routine only finds one.
!
!  Modified:
!
!    12 October 2010
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Combinatorial Heuristic Algorithms in FORTRAN,
!    Springer Verlag, 1986,
!    LC: QA402.5 L37.
!
!  Input:
!
!    integer NNODE, the number of nodes in the graph.
!
!    integer NEDGE, the number of edges in the graph.
!
!    integer INODE(NEDGE), JNODE(NEDGE), the two
!    end nodes of each edge.
!
!  Output:
!
!    integer LOOP(NEDGE), the Euler circuit, as a series of nodes.
!
  implicit none

  integer nedge
  integer nnode

  logical copyon
  logical found
  integer i
  integer ibase
  integer iforwd
  integer inode(nedge)
  integer insert
  integer ipivot
  integer iwork1(nedge)
  integer iwork2(nedge)
  integer iwork3(nnode)
  integer iwork4(nnode)
  integer iwork5(nnode)
  integer iwork6(nnode)
  integer j
  integer jnode(nedge)
  integer k
  integer l
  integer locbas
  integer loop(nedge)
  integer nbreak
  integer ncopy
  integer numarc
  integer numnode
!
!  The number of times each node has been visited begins at 0.
!
  iwork3(1:nnode) = 0
  loop(1:nedge) = 0
  iwork1(1:nedge) = 0
  iwork2(1:nedge) = 0
!
!  Begin the Euler circuit with the edge INODE(1), JNODE(1).
!
  numarc = 1
  iwork2(1) = 1

  numnode = 1
  i = inode(1)
  iwork1(numnode) = i
  iwork3(i) = 1

  numnode = numnode + 1
  j = jnode(1)
  iwork1(numnode) = j
  iwork3(j) = 1

  ibase = j
  nbreak = 0
!
!  Look for the next arc.
!
30    continue

  do i = 2, nedge

    if ( iwork2(i) == 0 ) then

      if ( inode(i) == ibase ) then

        found = .true.
        ibase = jnode(i)

      else if ( jnode(i) == ibase ) then

        found = .true.
        ibase = inode(i)

      else

        found = .false.

      end if

      if ( found ) then
        iwork2(i) = 1
        numarc = numarc + 1
        numnode = numnode + 1
        if ( numnode <= nedge ) then
          iwork1(numnode) = ibase
        end if
        iwork3(ibase) = 1
        go to 30
      end if

    end if

  end do
!
!  A cycle has been found.
!
  if ( 0 < nbreak ) then
    numnode = numnode - 1
    iwork5(nbreak) = numnode
  end if

  if ( numarc < nedge )  then

    iwork1(numnode) = ibase
!
!  Find a node in the current Euler circuit.
!
    do i = 2, nedge

      if ( iwork2(i) == 0 ) then

        if ( iwork3(inode(i)) /= 0 ) then

          found = .true.
          j = inode(i)
          k = jnode(i)

        else if ( iwork3(jnode(i)) /= 0 ) then

          found = .true.
          j = jnode(i)
          k = inode(i)

        else

          found = .false.

        end if
!
!  Identify a path which will be added to the circuit.
!
        if ( found ) then
          nbreak = nbreak + 1
          iwork6(nbreak) = j
          ibase = k
          iwork3(k) = 1
          numnode = numnode + 1
          iwork4(nbreak) = numnode
          iwork1(numnode) = ibase
          iwork2(i) = 1
          numarc = numarc + 1
          go to 30
        end if

      end if

    end do

  end if
!
!  Form the Euler circuit.
!
  if ( nbreak == 0 ) then
    numnode = numnode - 1
    loop(1:numnode) = iwork1(1:numnode)
    return
  end if

  insert = 1
  ipivot = iwork6(insert)
  iforwd = 0

  do

    ncopy = 1
    ibase = iwork1(1)
    locbas = 1
    loop(ncopy) = ibase
!
!  A path identified before is added to the circuit.
!
80  continue

    if ( ibase == ipivot ) then

      j = iwork4(insert) + iforwd
      k = iwork5(insert) + iforwd

      do l = j, k
        ncopy = ncopy + 1
        loop(ncopy) = iwork1(l)
        iwork1(l) = 0
      end do

      ncopy = ncopy + 1
!
!  Add the intersecting node to the circuit.
!
      loop(ncopy) = ibase
      iforwd = iforwd + 1

      if ( ncopy < numnode ) then

        do

          if ( nedge <= ncopy ) then
            exit
          end if

          locbas = locbas + 1

          if ( nedge <= locbas ) then
            exit
          end if

          ibase = iwork1(locbas)

          if ( ibase /= 0 ) then
            ncopy = ncopy + 1
            loop(ncopy) = ibase
          end if

        end do

      end if

    else

      ncopy = ncopy + 1

       if ( ncopy <= numnode ) then
         locbas = locbas + 1
         ibase = iwork1(locbas)
         loop(ncopy) = ibase
         go to 80
       end if

     end if
!
!  Check if more paths are to be added to the circuit.
!
    copyon = .false.
    insert = insert + 1

    if ( insert <= nbreak ) then
      copyon = .true.
      ipivot = iwork6(insert)
    end if

    if ( .not. copyon ) then
      exit
    end if

    iwork1(1:nedge) = loop(1:nedge)

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
!  Modified:
!
!    03 January 2004
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Combinatorial Heuristic Algorithms in FORTRAN,
!    Springer Verlag, 1986,
!    ISBN: 3540171614,
!    LC: QA402.5.L37.
!
!  Input:
!
!    integer NNODE, the number of nodes.
!
!    real ( kind = rk ) DIST(NNODE,NNODE), the distance matrix.
!    DIST(I,J) is the distance from node I to node J.  The matrix
!    should be symmetric.  If there is no arc from node I to node J,
!    set DIST(I,J) = HUGE(1.0).
!
!  Output:
!
!    integer INODE(NNODE), JNODE(NNODE); entries 1
!    through NNODE-1 describe the edges of the spanning tree as pairs of nodes.
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

  work(1:nnode) = huge ( d )
  iwork1(1:nnode) = 0
  iwork2(1:nnode) = 0
!
!  Find the first non-zero arc.
!
  i = 0

  do ij = 1, nnode
    do kj = 1, nnode
      if ( dist(ij,kj) < huge ( d ) ) then
        i = ij
        exit
      end if
    end do

    if ( i /= 0 ) then
      exit
    end if

  end do

  if ( i == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'graph_dist_min_span_tree(): Fatal error!'
    write ( *, '(a)' ) '  There are no nonzero arcs in this graph.'
    stop 1
  end if

  work(i) = 0.0D+00
  iwork1(i) = 1
  tree_length = 0.0D+00
  kk4 = nnode - 1

  do jj = 1, kk4

    work(1:nnode) = huge ( d )

    do i = 1, nnode
!
!  For each forward arc originating at node I calculate
!  the length of the path to node I
!
      if ( iwork1(i) == 1 ) then
        do j = 1, nnode
          if ( dist(i,j) < huge ( 1.0D+00 ) .and. iwork1(j) == 0 ) then
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
subroutine pmatch ( n, cost, pair )

!*****************************************************************************80
!
!! pmatch() finds a minimum weight perfect matching in a graph.
!
!  Modified:
!
!    15 September 1999
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Rainer Burkard, Ulrich Derigs,
!    Assignment and Matching Problems: Solution methods with
!    FORTRAN programs,
!    Lecture Notes in Economics and Mathematical Systems,
!    Volume 184,
!    Springer, 1980,
!    ISBN: 0387102671,
!    LC: QA402.5.B86.
!
!    Hang Tong Lau,
!    Combinatorial Heuristic Algorithms in FORTRAN,
!    Springer Verlag, 1986,
!    ISBN: 3540171614,
!    LC: QA402.5.L37.
!
!  Input:
!
!    integer N, the number of nodes in the complete graph.
!    N is assumed to be even.
!
!    real ( kind = rk ) COST(N*(N-1)/2), the strict upper triangle of the
!    cost matrix, stored by rows.  In other words, the first elements of
!    COST are the costs C(1,2), C(1,3), ..., C(1,N), C(2,3), C(2,4),
!    ..., C(2,N).
!
!  Output:
!
!    integer PAIR(N), contains the minimum weight perfect
!    matching.  Node I is connected to node PAIR(I), for I = 1 to N.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) cost((n*(n-1))/2)
  real ( kind = rk ) cst
  real ( kind = rk ) cstlow
  real ( kind = rk ) cswk
  real ( kind = rk ) cwk2
  integer i
  integer ihead
  integer index
  integer isub
  integer j
  integer jwk1(n)
  integer jwk2(n)
  integer jwk3(n)
  integer jwk4(n)
  integer jwk5(n)
  integer jwk6(n)
  integer jwk7(n)
  integer jwk8(n)
  integer jwk9(n)
  integer kk1
  integer kk2
  integer kk3
  integer kk4
  integer kk5
  integer kk6
  integer ll1
  integer ll2
  integer ll3
  integer ll4
  integer ll5
  integer max
  integer min
  integer mm1
  integer mm2
  integer mm3
  integer mm4
  integer mm5
  integer nn
  integer nn2
  integer pair(n)
  real ( kind = rk ) xcst
  real ( kind = rk ) value
  real ( kind = rk ) work1(n)
  real ( kind = rk ) work2(n)
  real ( kind = rk ) work3(n)
  real ( kind = rk ) work4(n)
  real ( kind = rk ) xwk2
  real ( kind = rk ) xwk3
  real ( kind = rk ) xwork

  if ( mod ( n, 2 ) /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'pmatch(): Fatal error!'
    write ( *, '(a)' ) '  The value of N must be even.'
    write ( *, '(a,i8)' ) '  The input value was N = ', n
    stop 1
  end if
!
!  Initialization
!
  nn2 = ( n * ( n - 1 ) ) / 2

  jwk1(2) = 0

  do i = 3, n
    jwk1(i) = jwk1(i-1) + i - 2
  end do

  ihead = n + 2

  do i = 1, n
    jwk2(i) = i
    jwk3(i) = i
    jwk5(i) = i
  end do

  jwk4(1:n) = 0
  jwk6(1:n) = ihead
  jwk7(1:n) = ihead
  jwk8(1:n) = ihead
  pair(1:n) = ihead
  work1(1:n) = huge ( 1.0D+00 )
  work2(1:n) = 0.0D+00
  work3(1:n) = 0.0D+00
  work4(1:n) = huge ( 1.0D+00 )
!
!  Start procedure.
!
  do i = 1, n

    if ( pair(i) == ihead ) then

      nn = 0
      cwk2 = huge ( 1.0D+00 )

      do j = 1, n

        min = i
        max = j

        if ( i /= j ) then

          if ( j < i ) then
            max = i
            min = j
          end if

          isub = jwk1(max) + min
          xcst = cost(isub)
          cswk = cost(isub) - work2(j)

          if ( cswk <= cwk2 ) then

            if ( cswk == cwk2 ) then
              if ( nn == 0 ) then
                go to 30
              else
                cycle
              end if
            end if

            cwk2 = cswk

            nn = 0

30          continue

            if ( pair(j) == ihead ) then
              nn = j
            end if

          end if

        end if

      end do

      if ( nn /= 0 ) then
        work2(i) = cwk2
        pair(i) = nn
        pair(nn) = i
      end if

    end if

  end do
!
!  Initial labeling.
!
  nn = 0

  do i = 1, n

    if ( pair(i) == ihead ) then

      nn = nn + 1
      jwk6(i) = 0
      work4(i) = 0.0D+00
      xwk2 = work2(i)

      do j = 1, n

        min = i
        max = j

        if ( i /= j ) then

          if ( j < i ) then
            max = i
            min = j
          end if

          isub = jwk1(max) + min
          xcst = cost(isub)
          cswk = cost(isub) - xwk2 - work2(j)

          if ( cswk < work1(j) ) then
            work1(j) = cswk
            jwk4(j) = i
          end if

        end if

      end do

    end if

  end do

  if ( nn <= 1 ) then
    go to 340
  end if
!
!  Examine the labeling and prepare for the next step.
!
80 continue

  cstlow = huge ( 1.0D+00 )

  do i = 1, n

    if ( jwk2(i) == i ) then

      cst = work1(i)

      if ( jwk6(i) < ihead ) then

        cst = 0.5D+00 * ( cst + work4(i) )

        if ( cst <= cstlow ) then
          index = i
          cstlow = cst
        end if

      else

        if ( jwk7(i) < ihead ) then
          if ( jwk3(i) /= i ) then
            cst = cst + work2(i)
            if ( cst < cstlow ) then
              index = i
              cstlow = cst
            end if
          end if
        else
          if ( cst < cstlow ) then
            index = i
            cstlow = cst
          end if
        end if

      end if

    end if

  end do

  if  ( jwk7(index) < ihead ) then
    go to 190
  end if

  if  ( jwk6(index) < ihead ) then

    ll4 = jwk4(index)
    ll5 = jwk5(index)
    kk4 = index
    kk1 = kk4
    kk5 = jwk2(ll4)
    kk2 = kk5

    do

      jwk7(kk1) = kk2
      mm5 = jwk6(kk1)

      if ( mm5 == 0 ) then
        exit
      end if

      kk2 = jwk2(mm5)
      kk1 = jwk7(kk2)
      kk1 = jwk2(kk1)

    end do

    ll2 = kk1
    kk1 = kk5
    kk2 = kk4

110 continue

    if ( ihead <= jwk7(kk1) ) then

      jwk7(kk1) = kk2
      mm5 = jwk6(kk1)

      if ( mm5 == 0 ) then
        go to 280
      end if

      kk2 = jwk2(mm5)
      kk1 = jwk7(kk2)
      kk1 = jwk2(kk1)
      go to 110

    end if

120 continue

    if ( kk1 == ll2 ) then
      go to 130
    end if

    mm5 = jwk7(ll2)
    jwk7(ll2) = ihead
    ll1 = pair(mm5)
    ll2 = jwk2(ll1)
    go to 120

  end if
!
!  Growing an alternating tree, add two edges.
!
  jwk7(index) = jwk4(index)
  jwk8(index) = jwk5(index)
  ll1 = pair(index)
  ll3 = jwk2(ll1)
  work4(ll3) = cstlow
  jwk6(ll3) = pair(ll3)

  call pmatch_sub_b ( ll3, n, nn2, cost, jwk1, jwk2, jwk3, jwk4, &
    jwk5, jwk7, jwk9, work1, work2, work3, work4 )

  go to 80
!
!  Shrink a blossom.
!
130 continue

  xwork = work2(ll2) + cstlow - work4(ll2)
  work2(ll2) = 0.0D+00
  mm1 = ll2

 do

    work3(mm1) = work3(mm1) + xwork
    mm1 = jwk3(mm1)

    if ( mm1 == ll2 ) then
      exit
    end if

  end do

  mm5 = jwk3(ll2)

  if ( ll2 /= kk5 ) then
    go to 160
  end if

150 continue

  kk5 = kk4
  kk2 = jwk7(ll2)

160 continue

  jwk3(mm1) = kk2
  ll1 = pair(kk2)
  jwk6(kk2) = ll1
  xwk2 = work2(kk2) + work1(kk2) - cstlow
  mm1 = kk2

  do

    mm2 = mm1
    work3(mm2) = work3(mm2) + xwk2
    jwk2(mm2) = ll2
    mm1 = jwk3(mm2)

    if ( mm1 == kk2 ) then
      exit
    end if

  end do

  jwk5(kk2) = mm2
  work2(kk2) = xwk2
  kk1 = jwk2(ll1)
  jwk3(mm2) = kk1
  xwk2 = work2(kk1) + cstlow - work4(kk1)
  mm2 = kk1

  do

    mm1 = mm2
    work3(mm1) = work3(mm1) + xwk2
    jwk2(mm1) = ll2
    mm2 = jwk3(mm1)

    if ( mm2 == kk1 ) then
      exit
    end if

  end do

  jwk5(kk1) = mm1
  work2(kk1) = xwk2

  if  ( kk5 /= kk1 ) then
    kk2 = jwk7(kk1)
    jwk7(kk1) = jwk8(kk2)
    jwk8(kk1) = jwk7(kk2)
    go to 160
  end if

  if ( kk5 /= index ) then
    jwk7(kk5) = ll5
    jwk8(kk5) = ll4
    if ( ll2 /= index ) then
      go to 150
    end if
  else
    jwk7(index) = ll4
    jwk8(index) = ll5
  end if

  jwk3(mm1) = mm5
  kk4 = jwk3(ll2)
  jwk4(kk4) = mm5
  work4(kk4) = xwork
  jwk7(ll2) = ihead
  work4(ll2) = cstlow

  call pmatch_sub_b ( ll2, n, nn2, cost, jwk1, jwk2, jwk3, jwk4, &
    jwk5, jwk7, jwk9, work1, work2, work3, work4 )

  go to 80
!
!  Expand a T-labeled blossom.
!
190   continue

  kk4 = jwk3(index)
  kk3 = kk4
  ll4 = jwk4(kk4)
  mm2 = kk4

  do

    mm1 = mm2
    ll5 = jwk5(mm1)
    xwk2 = work2(mm1)

    do

      jwk2(mm2) = mm1
      work3(mm2)= work3(mm2) - xwk2

      if ( mm2 == ll5 ) then
        exit
      end if

      mm2 = jwk3(mm2)

    end do

    mm2 = jwk3(ll5)
    jwk3(ll5) = mm1

    if ( mm2 == ll4 ) then
      exit
    end if

  end do

  xwk2 = work4(kk4)
  work2(index) = xwk2
  jwk3(index) = ll4
  mm2 = ll4

  do

    work3(mm2) = work3(mm2) - xwk2

    if ( mm2 == index ) then
      exit
    end if

    mm2 = jwk3(mm2)

  end do

  mm1 = pair(index)
  kk1 = jwk2(mm1)
  mm2 = jwk6(kk1)
  ll2 = jwk2(mm2)

  if ( ll2 /= index ) then

    kk2 = ll2

    do

      mm5 = jwk7(kk2)
      kk1 = jwk2(mm5)

      if ( kk1 == index ) then
        exit
      end if

      kk2 = jwk6(kk1)
      kk2 = jwk2(kk2)

    end do

    jwk7(ll2) = jwk7(index)
    jwk7(index) = jwk8(kk2)
    jwk8(ll2) = jwk8(index)
    jwk8(index) = mm5
    mm3 = jwk6(ll2)
    kk3 = jwk2(mm3)
    mm4 = jwk6(kk3)
    jwk6(ll2) = ihead
    pair(ll2) = mm1
    kk1 = kk3

    do

      mm1 = jwk7(kk1)
      mm2 = jwk8(kk1)
      jwk7(kk1) = mm4
      jwk8(kk1) = mm3
      jwk6(kk1) = mm1
      pair(kk1) = mm1
      kk2 = jwk2(mm1)
      pair(kk2) = mm2
      mm3 = jwk6(kk2)
      jwk6(kk2) = mm2

      if ( kk2 == index ) then
        exit
      end if

      kk1 = jwk2(mm3)
      mm4 = jwk6(kk1)
      jwk7(kk2) = mm3
      jwk8(kk2) = mm4

    end do

  end if

  mm2 = jwk8(ll2)
  kk1 = jwk2(mm2)
  work1(kk1) = cstlow
  kk4 = 0

  if ( kk1 /= ll2 ) then

    mm1 = jwk7(kk1)
    kk3 = jwk2(mm1)
    jwk7(kk1) = jwk7(ll2)
    jwk8(kk1) = mm2

    do

      mm5 = jwk6(kk1)
      jwk6(kk1) = ihead
      kk2 = jwk2(mm5)
      mm5 = jwk7(kk2)
      jwk7(kk2) = ihead
      kk5 = jwk8(kk2)
      jwk8(kk2) = kk4
      kk4 = kk2
      work4(kk2) = cstlow
      kk1 = jwk2(mm5)
      work1(kk1) = cstlow

      if ( kk1 == ll2 ) then
        exit
      end if

    end do

    jwk7(ll2) = kk5
    jwk8(ll2) = mm5
    jwk6(ll2) = ihead

    if ( kk3 == ll2 ) then
      go to 270
    end if

  end if

  kk1 = 0
  kk2 = kk3

  do

    mm5 = jwk6(kk2)
    jwk6(kk2) = ihead
    jwk7(kk2) = ihead
    jwk8(kk2) = kk1
    kk1 = jwk2(mm5)
    mm5 = jwk7(kk1)
    jwk6(kk1) = ihead
    jwk7(kk1) = ihead
    jwk8(kk1) = kk2
    kk2 = jwk2(mm5)

    if ( kk2 == ll2 ) then
      exit
    end if

  end do

  call pmatch_sub_a ( kk1, n, nn2, cost, jwk1, jwk2, jwk3, jwk4, &
    jwk5, jwk6, jwk8, work1, work2, work3, work4 )

270  continue

  if ( kk4 == 0 ) then
    go to 80
  end if

  ll2 = kk4

  call pmatch_sub_b ( ll2, n, nn2, cost, jwk1, jwk2, jwk3, jwk4, &
    jwk5, jwk7, jwk9, work1, work2, work3, work4 )

  kk4 = jwk8(ll2)
  jwk8(ll2) = ihead
  go to 270
!
!  Augmentation of the matching.
!
!  Exchange the matching and non-matching edges along the augmenting path.
!
 280  continue

  ll2 = kk4
  mm5 = ll4

  do

    kk1 = ll2

    do

      pair(kk1) = mm5
      mm5 = jwk6(kk1)
      jwk7(kk1) = ihead

      if ( mm5 == 0 ) then
        exit
      end if

      kk2 = jwk2(mm5)
      mm1 = jwk7(kk2)
      mm5 = jwk8(kk2)
      kk1 = jwk2(mm1)
      pair(kk2) = mm1

    end do

    if ( ll2 /= kk4 ) then
      exit
    end if

    ll2 = kk5
    mm5 = ll5

  end do
!
!  Remove all labels of non-exposed base nodes.
!
  do i = 1, n

    if ( jwk2(i) == i ) then

      if ( jwk6(i) < ihead ) then

        cst = cstlow - work4(i)
        work2(i) = work2(i) + cst
        jwk6(i) = ihead

        if ( pair(i) /= ihead ) then
          work4(i) = huge ( 1.0D+00 )
        else
          jwk6(i) = 0
          work4(i) = 0.0D+00
        end if

      else

        if ( jwk7(i) < ihead ) then
          cst = work1(i) - cstlow
          work2(i) = work2(i) + cst
          jwk7(i) = ihead
          jwk8(i) = ihead
        end if

        work4(i) = huge ( 1.0D+00 )

      end if

      work1(i) = huge ( 1.0D+00 )

    end if

  end do

  nn = nn - 2
!
!  Determine the new WORK1 values.
!
  if ( nn <= 1 ) then
    go to 340
  end if

  do i = 1, n

    kk1 = jwk2(i)

    if ( jwk6(kk1) == 0 ) then

      xwk2 = work2(kk1)
      xwk3 = work3(i)

      do j = 1, n

        kk2 = jwk2(j)

        if ( kk1 /= kk2 ) then

          min = i
          max = j

          if ( i /= j ) then

            if ( j < i ) then
              max = i
              min = j
            end if

            isub = jwk1(max) + min
            xcst = cost(isub)
            cswk = cost(isub) - xwk2 - xwk3
            cswk = cswk - work2(kk2) - work3(j)

            if ( cswk < work1(kk2) ) then
              jwk4(kk2) = i
              jwk5(kk2) = j
              work1(kk2) = cswk
            end if

          end if
        end if

      end do

    end if

  end do

  go to 80
!
!  Generate the original graph by expanding all shrunken blossoms.
!
340   continue

  value = 0.0D+00

  do i = 1, n

    if ( jwk2(i) == i ) then

      if ( 0 <= jwk6(i) ) then

        kk5 = pair(i)
        kk2 = jwk2(kk5)
        kk4 = pair(kk2)
        jwk6(i) = -1
        jwk6(kk2) = -1
        min = kk4
        max = kk5

        if ( kk4 /= kk5 ) then

          if ( kk5 < kk4 ) then
            max = kk4
            min = kk5
          end if

          isub = jwk1(max) + min
          xcst = cost(isub)
          value = value + xcst
        end if

      end if

    end if

  end do

  do i = 1, n

360 continue

    ll2 = jwk2(i)

    if ( ll2 == i ) then
      cycle
    end if

    mm2 = jwk3(ll2)
    ll4 = jwk4(mm2)
    kk3 = mm2
    xwork = work4(mm2)

    do

      mm1 = mm2
      ll5 = jwk5(mm1)
      xwk2 = work2(mm1)

      do

        jwk2(mm2) = mm1
        work3(mm2) = work3(mm2) - xwk2

        if ( mm2 == ll5 ) then
          exit
        end if

        mm2 = jwk3(mm2)

      end do

      mm2 = jwk3(ll5)
      jwk3(ll5) = mm1

      if ( mm2 == ll4 ) then
        exit
      end if

    end do

    work2(ll2) = xwork
    jwk3(ll2) = ll4
    mm2 = ll4

    do

      work3(mm2) = work3(mm2) - xwork

      if ( mm2 == ll2 ) then
        exit
      end if

      mm2 = jwk3(mm2)

    end do

    mm5 = pair(ll2)
    mm1 = jwk2(mm5)
    mm1 = pair(mm1)
    kk1 = jwk2(mm1)

    if ( ll2 /= kk1 ) then

      pair(kk1) = mm5
      kk3 = jwk7(kk1)
      kk3 = jwk2(kk3)

      do

        mm3 = jwk6(kk1)
        kk2 = jwk2(mm3)
        mm1 = jwk7(kk2)
        mm2 = jwk8(kk2)
        kk1 = jwk2(mm1)
        pair(kk1) = mm2
        pair(kk2) = mm1
        min = mm1
        max = mm2

        if ( mm1 == mm2 ) then
          go to 360
        end if

        if ( mm2 < mm1 ) then
          max = mm1
          min = mm2
        end if

        isub = jwk1(max) + min
        xcst = cost(isub)
        value = value + xcst

        if ( kk1 == ll2 ) then
          exit
        end if

      end do

      if ( kk3 == ll2 ) then
        go to 360
      end if

    end if

    do

      kk5 = jwk6(kk3)
      kk2 = jwk2(kk5)
      kk6 = jwk6(kk2)
      min = kk5
      max = kk6

      if ( kk5 == kk6 ) then
        go to 360
      end if

      if ( kk6 < kk5 ) then
        max = kk5
        min = kk6
      end if

      isub = jwk1(max) + min
      xcst = cost(isub)
      value = value + xcst
      kk6 = jwk7(kk2)
      kk3 = jwk2(kk6)

      if ( kk3 == ll2 ) then
        go to 360
      end if

    end do

  end do

  return
end
subroutine pmatch_sub_a ( kk, n, nn2, cost, jwk1, jwk2, jwk3, jwk4, &
  jwk5, jwk6, jwk8, work1, work2, work3, work4 )

!*****************************************************************************80
!
!! pmatch_sub_a() is used by pmatch().
!
!  Modified:
!
!    11 September 1999
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Combinatorial Heuristic Algorithms in FORTRAN,
!    Springer Verlag, 1986,
!    ISBN: 3540171614,
!    LC: QA402.5.L37.
!
!  Input:
!
!  Output:
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nn2

  real ( kind = rk ) cost(nn2)
  real ( kind = rk ) cstwk
  real ( kind = rk ) cswk
  integer i
  integer ihead
  integer isub
  integer j
  integer jj1
  integer jj2
  integer jj3
  integer jj4
  integer jwk1(n)
  integer jwk2(n)
  integer jwk3(n)
  integer jwk4(n)
  integer jwk5(n)
  integer jwk6(n)
  integer jwk8(n)
  integer kk
  integer max
  integer min
  real ( kind = rk ) work1(n)
  real ( kind = rk ) work2(n)
  real ( kind = rk ) work3(n)
  real ( kind = rk ) work4(n)
  real ( kind = rk ) xwk2
  real ( kind = rk ) xwk3

  ihead = n + 2

  do

    jj1 = kk
    kk  = jwk8(jj1)
    jwk8(jj1) = ihead
    cstwk = huge ( 1.0D+00 )
    jj3 = 0
    jj4 = 0
    j = jj1
    xwk2 = work2(jj1)

    do

      xwk3 = work3(j)

      do i = 1, n

        jj2 = jwk2(i)

        if ( jwk6(jj2) < ihead ) then

          min = j
          max = i

          if ( j /= i ) then

            if ( i < j ) then
              max = j
              min = i
            end if

            isub = jwk1(max) + min

            cswk = cost(isub) - xwk2 - xwk3 - work2(jj2) - work3(i) + work4(jj2)

            if ( cswk < cstwk ) then
              jj3 = i
              jj4 = j
              cstwk = cswk
            end if

          end if

        end if

      end do

      j = jwk3(j)

      if ( j == jj1 ) then
        exit
      end if

    end do

    jwk4(jj1) = jj3
    jwk5(jj1) = jj4
    work1(jj1) = cstwk

    if ( kk == 0 ) then
      exit
    end if

  end do

  return
end
subroutine pmatch_sub_b ( kk, n, nn2, cost, jwk1, jwk2, jwk3, jwk4, &
  jwk5, jwk7, jwk9, work1, work2, work3, work4 )

!*****************************************************************************80
!
!! pmatch_sub_b() is used by pmatch().
!
!  Modified:
!
!    11 September 1999
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Combinatorial Heuristic Algorithms in FORTRAN,
!    Springer Verlag, 1986,
!    ISBN: 3540171614,
!    LC: QA402.5.L37.
!
!  Input:
!
!  Output:
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nn2

  real ( kind = rk ) cost(nn2)
  real ( kind = rk ) cswk
  integer i
  integer ihead
  integer ii
  integer isub
  integer jj1
  integer jj2
  integer jj3
  integer jwk1(n)
  integer jwk2(n)
  integer jwk3(n)
  integer jwk4(n)
  integer jwk5(n)
  integer jwk7(n)
  integer jwk9(n)
  integer kk
  integer max
  integer min
  real ( kind = rk ) work1(n)
  real ( kind = rk ) work2(n)
  real ( kind = rk ) work3(n)
  real ( kind = rk ) work4(n)
  real ( kind = rk ) xwk1
  real ( kind = rk ) xwk2

  ihead = n + 2
  xwk1 = work4(kk) - work2(kk)
  work1(kk) = huge ( 1.0D+00 )
  xwk2 = xwk1 - work3(kk)
  jwk7(kk) = 0
  ii = 0

  do i = 1, n

    jj3 = jwk2(i)

    if ( ihead <= jwk7(jj3) ) then

      ii = ii + 1
      jwk9(ii) = i
      min = kk
      max = i

      if ( kk /= i ) then

        if ( i < kk ) then
          max = kk
          min = i
        end if

        isub = jwk1(max) + min

        cswk = cost(isub) + xwk2 - work2(jj3) - work3(i)

        if ( cswk < work1(jj3) ) then
          jwk4(jj3) = kk
          jwk5(jj3) = i
          work1(jj3) = cswk
        end if
      end if
    end if

  end do

  jwk7(kk) = ihead
  jj1 = kk
  jj1 = jwk3(jj1)

  if ( jj1 == kk ) then
    return
  end if

  do

    xwk2 = xwk1 - work3(jj1)

    do i = 1, ii

      jj2 = jwk9(i)
      jj3 = jwk2(jj2)
      min = jj1
      max = jj2

      if ( jj1 /= jj2 ) then

        if ( jj2 < jj1 ) then
          max = jj1
          min = jj2
        end if

        isub = jwk1(max) + min

        cswk = cost(isub) + xwk2 - work2(jj3) - work3(jj2)

        if ( cswk < work1(jj3) ) then
          jwk4(jj3) = jj1
          jwk5(jj3) = jj2
          work1(jj3) = cswk
        end if

      end if

    end do

    jj1 = jwk3(jj1)

    if ( jj1 == kk ) then
      exit
    end if

  end do

  return
end
subroutine tsp ( nnode, dist, isol )

!*****************************************************************************80
!
!! tsp() is a heuristic algorithm for the traveling salesman problem.
!
!  Discussion:
!
!    Let G be a complete graph with an associated distance matrix
!    DIST(I,J) on its edges.  The traveling salesman problem is to
!    start from a node in G, visit all the other nodes exactly once,
!    and return to the starting node in such a way that the total
!    traveled distance is a minimum.
!
!    In general, it is very hard to develop efficient algorithms
!    that yield good approximate solutions to this problem.  However,
!    in the special case where the distance matrix is symmetric:
!
!      DIST(I,J) = DIST(J,I)
!
!    and satisfies the triangle inequality:
!
!      DIST(I,K) <= DIST(I,J) + DIST(J,K),
!
!    solutions close to the optimum value can be found in a relatively
!    short time.  The algorithm here will find a circuit of no worse
!    than 3/2 the optimum length, in polynomial time.
!
!  Modified:
!
!    12 October 2010
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Nicos Christofides,
!    Worst-case analysis of a new heuristic for the traveling salesman
!    problem,
!    Management Science Research Report Number 388,
!    Carnegie-Mellon University, 1976.
!
!    Hang Tong Lau,
!    Combinatorial Heuristic Algorithms in FORTRAN,
!    Springer Verlag, 1986,
!    ISBN: 3540171614,
!    LC: QA402.5.L37.
!
!  Input:
!
!    integer NNODE, the number of nodes.
!
!    real ( kind = rk ) DIST(NNODE,NNODE), the distance between
!    each pair of nodes.
!
!  Output:
!
!    integer ISOL(NNODE), contains the nodes
!    in the order in which they should be visited.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nnode

  real ( kind = rk ) dist(nnode,nnode)
  integer i
  integer ict
  integer ii
  integer inode((3*nnode)/2)
  integer isol(nnode)
  integer iwk10((3*nnode)/2)
  integer iwk11(nnode)
  integer iwk18(nnode)
  integer j
  integer jj
  integer jnode((3*nnode)/2)
  integer k
  integer match
  integer nedge
  real ( kind = rk ) wk1( ( nnode * ( nnode - 1 ) ) / 2)
!
!  From the graph described by DIST, construct a minimum length spanning
!  tree, described by INODE, JNODE.
!
  call graph_dist_min_span_tree3 ( nnode, dist, inode, jnode )
!
!  Compute the degree (in the spanning tree) of each node.
!
  iwk10(1:nnode) = 0

  do i = 1, nnode - 1
    ii = inode(i)
    jj = jnode(i)
    iwk10(ii) = iwk10(ii) + 1
    iwk10(jj) = iwk10(jj) + 1
  end do
!
!  Construct a list of the nodes of odd degree.
!
  ict = 0
  do i = 1, nnode
    if ( mod ( iwk10(i), 2 ) /= 0 ) then
      ict = ict + 1
      iwk11(ict) = i
    end if
  end do
!
!  Determine a minimum weight perfect matching for the set of odd-degree nodes.
!
  k = 0
  do i = 2, ict
    do j = 1, i - 1
      k = k + 1
      wk1(k) = dist(iwk11(j),iwk11(i))
    end do
  end do

  call pmatch ( ict, wk1, iwk18 )
!
!  Store up the edges in the perfect matching.
!
  do i = 1, ict
    iwk18(i) = iwk11(iwk18(i))
  end do

  iwk10(1:nnode) = 0

  k = nnode - 1

  do i = 1, ict

    if ( iwk10(iwk11(i)) == 0 ) then
      iwk10(iwk11(i)) = 1
      iwk10(iwk18(i)) = 1
      k = k + 1
      inode(k) = iwk11(i)
      jnode(k) = iwk18(i)
    end if

  end do
!
!  Find an Euler circuit.
!
  nedge = nnode - 1 + ( ict / 2 )
  call graph_arc_euler_circ ( nnode, nedge, inode, jnode, iwk10 )
!
!  Form the Hamiltonian circuit.
!
  isol(1) = iwk10(1)
  j = 2
  isol(j) = iwk10(j)
  k = 2

  do

    match = 1

    do while ( match == 1 )

      j = j + 1

      match = 0

      do i = 1, j - 1
        if ( iwk10(j) == iwk10(i) ) then
          match = 1
          exit
        end if
      end do

    end do

    k = k + 1
    isol(k) = iwk10(j)

    if ( nnode <= k ) then
      exit
    end if

  end do

  return
end
