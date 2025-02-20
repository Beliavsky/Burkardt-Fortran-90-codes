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
subroutine maze_diameter ( bar, degree, diam, flat, m, n, path, istart, &
  jstart, istop, jstop )

!*****************************************************************************80
!
!! maze_diameter() computes the "diameter" of a maze that has no circuits.
!
!  Discussion:
!
!    The routine also returns two cells, (ISTART,JSTART), and (ISTOP,JSTOP)
!    which are separated by a path of length DIAM.
!
!  Definition:
!
!    The diameter is the length of the longest path that never passes
!    through the same cell twice.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    24 March 2022
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer BAR(M,N+1), records the vertical "bars" 
!    in the maze.
!    -1, means "indefinite", that there is no cell of the maze on either
!        side of this position;
!     0, means "wall", that there is a cell on at least one side, and
!        a wall here;
!     1, means "open", that there are cells on both sides (or possibly
!        an opening to the exterior) and the way is open.
!
!    Output, integer DEGREE(M,N), the degree of each node.
!
!    Output, integer DIAM, the length of the longest path 
!    in the tree.
!
!    Input, integer FLAT(M+1,N), records the horizontal "flats"
!    in the maze.
!    -1, means "indefinite", that there is no cell of the maze on either
!        side of this position;
!     0, means "wall", that there is a cell on at least one side, and
!        a wall here;
!     1, means "open", that there are cells on both sides (or possibly
!        an opening to the exterior) and the way is open.
!
!    Input, integer M, N, the number of rows and columns of cells.
!
!    Output, integer PATH(M,N), marks the path between the cells
!    (ISTART,JSTART) and (ISTOP,JSTOP).  A cell (I,J) is in the path
!    if PATH(I,J) is 1.
!
!    Output, integer ISTART, JSTART, are the I and J cell
!    coordinates of the starting cell.
!
!    Output, integer ISTOP, JSTOP, are the I and J cell 
!    coordinates of the goal cell.
!
  implicit none

  integer, parameter :: OPEN = 1
  integer, parameter :: SHUT = 2

  integer m
  integer n

  integer bar(m,n+1)
  integer degree(m,n)
  integer diam
  integer flat(m+1,n)
  integer i
  integer i2
  integer invals
  integer istart
  integer istop
  integer j
  integer j2
  integer jstart
  integer jstop
  integer k
  integer kstep
  integer n1
  integer n2
  integer path(m,n)

  if ( m * n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'maze_diameter(): Fatal error!'
    write ( *, '(a)' ) '  M*N <= 0.'
    stop 1
  else if ( m * n == 1 ) then
    diam = 0
    return
  end if

  k = 0
  do j = 1, n
    do i = 1, m
      k = k + 1
      path(i,j) = k
    end do
  end do
!
!  On step KSTEP:
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
    do j = 1, n
      do i = 1, m

        degree(i,j) = 0

        if ( flat(i,j) == OPEN ) then
          degree(i,j) = degree(i,j) + 1
        end if
        if ( flat(i+1,j) == OPEN ) then
          degree(i,j) = degree(i,j) + 1
        end if
        if ( bar(i,j) == OPEN ) then
          degree(i,j) = degree(i,j) + 1
        end if
        if ( bar(i,j+1) == OPEN ) then
          degree(i,j) = degree(i,j) + 1
        end if

      end do
    end do
!
!  Count the number of interior nodes.
!
    invals = 0
    do j = 1, n
      do i = 1, m
        if ( 1 < degree(i,j) ) then
          invals = invals + 1
        end if
      end do
    end do
!
!  If there are at least two interior nodes, then chop off the
!  terminal nodes and pass their labels inward.
!
    if ( 2 <= invals ) then

      k = 0

      do j = 1, n

        do i = 1, m

          k = k + 1

          if ( degree(i,j) == 1 ) then

            if ( flat(i,j) == OPEN ) then
              i2 = i - 1
              j2 = j
              flat(i,j) = SHUT
            else if ( flat(i+1,j) == OPEN ) then
              i2 = i + 1
              j2 = j
              flat(i+1,j) = SHUT
            else if ( bar(i,j) == OPEN ) then
              i2 = i
              j2 = j - 1
              bar(i,j) = SHUT
            else if ( bar(i,j+1) == OPEN ) then
              i2 = i
              j2 = j + 1
              bar(i,j+1) = SHUT
            end if

            path(i2,j2) = path(i,j)

          end if
 
        end do

      end do
!
!  But if there are 1 or 0 interior nodes, it's time to stop.
!
    else if ( invals == 1 ) then

      diam = 2 * kstep + 2
      exit

    else if ( invals == 0 ) then

      diam = 2 * kstep + 1
      exit

    end if

  end do
!
!  Now get the labels from two of the remaining terminal nodes.
!  The nodes represented by these labels will be a diameter apart.
!
  n1 = 0
  n2 = 0

  do j = 1, n
    do i = 1, m

      if ( degree(i,j) == 1 ) then
        if ( n1 == 0 ) then
          n1 = path(i,j)
        else if ( n2 == 0 ) then
          n2 = path(i,j)
        end if
      end if

    end do
  end do
!
!  Set the labels of the interior node (if any) and nodes marked
!  N1 and N2 to 1, and all others to 0.  This will label the nodes on the path.
!
  if ( invals == 1 ) then

    do j = 1, n
      do i = 1, m
        if ( 1 < degree(i,j) ) then
          path(i,j) = 1
        end if
      end do
    end do

  end if

  do j = 1, n
    do i = 1, m

      if ( path(i,j) == n1 .or. path(i,j) == n2 ) then
        path(i,j) = 1
      else
        path(i,j) = 0
      end if

    end do
  end do
!
!  Translate N1 and N2 to row, column.
!
  jstart = ( n1 - 1 ) / m + 1
  istart = n1 - ( jstart - 1 ) * m

  jstop = ( n2 - 1 ) / m + 1
  istop = n2 - ( jstop - 1 ) * m
!
!  Clean up the DEGREE and LINKS arrays.
!
  do i = 1, m
    do j = 1, n + 1
      if ( bar(i,j) == SHUT ) then
        bar(i,j) = OPEN
      end if
    end do
  end do

  do i = 1, m + 1
    do j = 1, n
      if ( flat(i,j) == SHUT ) then
        flat(i,j) = OPEN
      end if
    end do
  end do

  do j = 1, n
    do i = 1, m

      degree(i,j) = 0

      if ( flat(i,j) == OPEN ) then
        degree(i,j) = degree(i,j) + 1
      end if
      if ( flat(i+1,j) == OPEN ) then
        degree(i,j) = degree(i,j) + 1
      end if
      if ( bar(i,j) == OPEN ) then
        degree(i,j) = degree(i,j) + 1
      end if
      if ( bar(i,j+1) == OPEN ) then
        degree(i,j) = degree(i,j) + 1
      end if

    end do
  end do

  return
end
subroutine maze_path ( bar, flat, m, n, istart, jstart, istop, jstop )

!*****************************************************************************80
!
!! maze_path() finds a path through a maze.
!
!  Warning: 
!
!    This routine has some stupid internal limits which could
!    be fixed by reprogramming.  (Use the BAR and FLAT arrays to record
!    the tentative path, for instance.)
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
!    Input/output, integer BAR(M,N+1), records the vertical 
!    "bars" in the maze, and on output, the path through open bars:
!    -1, means "indefinite", that there is no cell of the maze on either
!        side of this position;
!     0, means "wall", that there is a cell on at least one side, and
!        a wall here;
!     1, means "open", that there are cells on both sides (or possibly
!        an opening to the exterior) and the way is open.
!     2, means the path goes through this open bar.
!
!    Input/output, integer FLAT(M+1,N), records the horizontal 
!    "flats" in the maze, and on output, the path through open flats:
!    -1, means "indefinite", that there is no cell of the maze on either
!        side of this position;
!     0, means "wall", that there is a cell on at least one side, and
!        a wall here;
!     1, means "open", that there are cells on both sides (or possibly
!        an opening to the exterior) and the way is open.
!     2, means the path goes through this open flat.
!
!    Input, integer M, N, the number of rows and columns of 
!    cells.
!
!    Input, integer ISTART, JSTART, are the I and J cell 
!    coordinates of the starting cell.
!
!    Input, integer ISTOP, JSTOP, are the I and J cell coordinates 
!    of the goal cell, which will be required to be a terminal node of the tree.
!
  implicit none

  integer, parameter :: maxpath = 200
  integer, parameter :: maxstack = 500
  integer, parameter :: maxused = 500
  integer, parameter :: OPEN = 1
  integer, parameter :: PATH = 2

  integer m
  integer n

  integer bar(m,n+1)
  integer flat(m+1,n)
  integer ipath
  integer istart
  integer istop
  integer ival
  integer ival2
  integer jstart
  integer jstop
  integer jval
  integer jval2
  integer kval
  integer kval2
  integer ncan
  integer npath
  integer nstack
  integer pathlist(maxpath)
  integer stack(maxstack)
  integer used(maxused)

  if ( istart < 1 .or. m < istart ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'maze_path(): Fatal error!'
    write ( *, '(a,i8)' ) '  ISTART out of range, = ', istart
    write ( *, '(a,i8)' ) '  Must be between 1 and ', m
    stop
  else if ( jstart < 1 .or. n < jstart ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'maze_path(): Fatal error!'
    write ( *, '(a,i8)' ) '  JSTART out of range, = ', jstart
    write ( *, '(a,i8)' ) '  Must be between 1 and ', n
    stop
  else if ( istop < 1 .or. m < istop ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'maze_path(): Fatal error!'
    write ( *, '(a,i8)' ) '  ISTOP out of range, = ', istop
    write ( *, '(a,i8)' ) '  Must be between 1 and ', m
    stop
  else if ( jstop < 1 .or. n < jstop ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'maze_path(): Fatal error!'
    write ( *, '(a,i8)' ) '  JSTOP out of range, = ', jstop
    write ( *, '(a,i8)' ) '  Must be between 1 and ', n
    stop
  end if

  if ( maxused < m * n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'maze_path(): Fatal error!'
    write ( *, '(a)' ) '  M * N greater than internal limit MAXUSED.'
    stop
  end if

  used(1:m*n) = 0
  pathlist(1:m*n) = 0
!
!  Begin the path at (ISTART,JSTART).
!
  npath = 1
  ival = istart
  jval = jstart
  kval = ( jval - 1 ) * m + ival
  pathlist(npath) = kval
  used(kval) = npath

  ncan = 0
  nstack = 1
  stack(nstack) = ncan
!
!  Try to take a new step.
!
  do
!
!  Find all the accessible never-used neighbors of the current endpoint.
!  Add them to the stack, and set NCAN to the number of candidates.
!
    ncan = 0

    if ( ival /= 1 ) then

      if ( flat ( ival, jval ) == OPEN ) then

        ival2 = ival - 1
        jval2 = jval
        kval2 = ( jval2 - 1 ) * m + ival2

        if ( used(kval2) == 0 ) then
          ncan = ncan + 1
          nstack = nstack + 1
          if ( maxstack < nstack ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'maze_path(): Fatal error!'
            write ( *, '(a)' ) '  The size of the internal stack was exceeded!'
            stop 1
          end if
          stack(nstack) = kval2
        end if

      end if

    end if
  
    if ( jval /= n ) then

      if ( bar ( ival, jval+1 ) == OPEN ) then

        ival2 = ival
        jval2 = jval + 1
        kval2 = ( jval2 - 1 ) * m + ival2

        if ( used(kval2) == 0 ) then
          ncan = ncan + 1
          nstack = nstack + 1
          if ( maxstack < nstack ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'maze_path(): Fatal error!'
            write ( *, '(a)' ) '  The size of the internal stack was exceeded!'
            stop 1
          end if
          stack(nstack) = kval2
        end if

      end if

    end if

    if ( jval /= 1 ) then

      if ( bar ( ival, jval ) == OPEN ) then
        ival2 = ival
        jval2 = jval - 1
        kval2 = ( jval2 - 1 ) * m + ival2

        if ( used(kval2) == 0 ) then
          ncan = ncan + 1
          nstack = nstack + 1
          if ( maxstack < nstack ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'maze_path(): Fatal error!'
            write ( *, '(a)' ) '  The size of the internal stack was exceeded!'
            stop 1
          end if
          stack(nstack) = kval2
        end if

      end if

    end if

    if ( ival /= m ) then

      if ( flat ( ival+1, jval ) == OPEN ) then

        ival2 = ival + 1
        jval2 = jval
        kval2 = ( jval2 - 1 ) * m + ival2

        if ( used(kval2) == 0 ) then
          ncan = ncan + 1
          nstack = nstack + 1
          if ( maxstack < nstack ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'maze_path(): Fatal error!'
            write ( *, '(a)' ) '  The size of the internal stack was exceeded!'
            stop 1
          end if
          stack(nstack) = kval2
        end if

      end if

    end if
!
!  Add NCAN to the stack.
!
    nstack = nstack + 1
    if ( maxstack < nstack ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'maze_path(): Fatal error!'
      write ( *, '(a)' ) '  The size of the internal stack was exceeded!'
      stop 1
    end if

    stack(nstack) = ncan

20  continue
!
!  If NCAN=0, then...
!
    if ( ncan == 0 ) then
!
!  ...if the current cell is the starting point, we've failed.
!
      if ( npath == 1 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'maze_path(): Fatal error!'
        write ( *, '(a)' ) '  Could not find a path to the goal.'
        return
!
!  ...Else drop the current endpoint, going back to previous cell
!  on the path, pop the stack one level, (getting new value of NCAN), 
!  go to 20.
!
      else

        used(kval) = - used(kval)

        npath = npath - 1
        kval = pathlist(npath)
        ival = mod ( kval, m )
        jval = 1 + ( kval - ival ) / m 

        nstack = nstack - 1
        ncan = stack(nstack)
        go to 20

      end if
!
!  Else, take one candidate off the stack, add it to the path,
!  mark it as used, set NCAN = NCAN-1.
!
    else 

      kval = stack(nstack-1)

      npath = npath + 1

      if ( maxpath < npath ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'maze_path(): Fatal error!'
        write ( *, '(a)' ) '  NPATH exceeds internal limit MAXPATH.'
        stop 1
      end if

      pathlist(npath) = kval

      used(kval) = npath

      jval = ( kval - 1 ) / m + 1
      ival = kval - ( jval - 1 ) * m

      ncan = ncan - 1
      nstack = nstack - 1
      stack(nstack) = ncan
!
!  If the candidate is not the goal, go to 10...
!
      if ( ival /= istop .or. jval /= jstop ) then
        cycle
      end if
!
!  ...else we're done.
!
      do ipath = 1, npath-1

        kval = pathlist(ipath)
        jval = ( kval - 1 ) / m + 1
        ival = kval - ( jval - 1 ) * m

        kval2 = pathlist(ipath+1)

        if ( kval2 == kval - 1 ) then
          flat(ival,jval) = PATH
        else if ( kval2 == kval + m ) then
          bar(ival,jval+1) = PATH
        else if ( kval2 == kval - m ) then
          bar(ival,jval) = PATH
        else if ( kval2 == kval + 1 ) then
          flat(ival+1,jval) = PATH
        end if

      end do

      return

    end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'maze_path(): Fatal error!'
  write ( *, '(a)' ) '  The size of the internal stack was exceeded!'
  stop 1
end
subroutine maze_print ( bar, flat, m, n, istart, jstart, istop, jstop, title )

!*****************************************************************************80
!
!! maze_print() prints out a maze and a path.
!
!  Example:
!
!    +--+--+
!    |*****|$$
!    +**+**+**+
!    |00|*****|
!    +  +--+--+
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
!    Input, integer BAR(M,N+1), records the vertical "bars" 
!    in the maze.
!    -1, means "indefinite", that there is no cell of the maze on either
!        side of this position;
!     0, means "wall", that there is a cell on at least one side, and
!        a wall here;
!     1, means "open", that there are cells on both sides (or possibly
!        an opening to the exterior) and the way is open.
!     2, means "path", that the way is open, and the path goes this way.
!
!    Input, integer FLAT(M+1,N), records the horizontal "flats" 
!    in the maze.
!    -1, means "indefinite", that there is no cell of the maze on either
!        side of this position;
!     0, means "wall", that there is a cell on at least one side, and
!        a wall here;
!     1, means "open", that there are cells on both sides (or possibly
!        an opening to the exterior) and the way is open.
!     2, means "path", that the way is open, and the path goes this way.
!
!    Input, integer M, N, the number of rows and columns of cells.
!    Currently, the program cannot handle a maze with more than 26 columns.
!
!    Input, integer ISTART, JSTART, are the I and J cell 
!    coordinates of the starting cell.  The starting cell will be marked 
!    "00".  If no starting cell is to be specified, set ISTART = JSTART = 0.
!
!    Input, integer ISTOP, JSTOP, are the I and J cell coordinates 
!    of the goal cell.  The goal cell will be marked "$$".  If no goal cell
!    is to be specified, set ISTOP = JSTOP = 0.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: NMAX = 26
  integer, parameter :: INDEF = -1
  integer, parameter :: WALL = 0
  integer, parameter :: OPEN = 1
  integer, parameter :: PATH = 2

  integer m
  integer n

  integer bar(m,n+1)
  integer flat(m+1,n)
  integer i
  integer ilo
  integer istart
  integer istop
  integer j
  integer jstart
  integer jstop
  integer nsafe
  character ( len = 3*(NMAX+1) ) string
  character ( len = * ) title

  if ( NMAX < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'maze_print(): Warning!'
    write ( *, '(a,i8)' ) '  N may not be more than ', NMAX
    write ( *, '(a)' ) '  Only a portion of the maze will be shown.'
  end if

  if ( len_trim ( title ) /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  nsafe = min ( n, NMAX )

  do i = 1, m

    string = ' '

    ilo = 1
    do j = 1, nsafe

      if ( flat(i,j) == WALL ) then

        string(ilo:ilo+3) = '+--+'

      else if ( flat(i,j) == OPEN ) then

        string(ilo:ilo+3) = '+  +'

      else if ( flat(i,j) == PATH ) then

        string(ilo:ilo+3) = '+**+'

      else if ( flat(i,j) == INDEF ) then

      end if

      ilo = ilo + 3

    end do

    write ( *, '(a)' ) string(1:ilo)

    string = ' '

    ilo = 1
    do j = 1, nsafe+1

      if ( bar(i,j) == WALL ) then

        string(ilo:ilo) = '|'

      else if ( bar(i,j) == OPEN ) then

        string(ilo:ilo) = ' '

      else if ( bar(i,j) == PATH ) then

        string(ilo:ilo) = '*'

      else if ( bar(i,j) == INDEF ) then

      end if
!
!  Now fill in the interior of the cell.
!
      if ( i == istart .and. j == jstart ) then

        string(ilo+1:ilo+2) = '00'

      else if ( i== istop .and. j == jstop ) then

        string(ilo+1:ilo+2) = '$$'

      else if ( bar(i,j) == PATH ) then

        string(ilo+1:ilo+2) = '**'

      else if ( j <= n ) then

        if ( flat(i,j) == PATH .or. bar(i,j+1) == PATH .or. &
             flat(i+1,j) == PATH ) then

          string(ilo+1:ilo+2) = '**'

        end if

      end if

      ilo = ilo + 3

    end do

    ilo = ilo - 3

    write ( *, '(a)' ) string(1:ilo)

  end do

  string = ' '
  i = m+1
  ilo = 1
  do j = 1, nsafe

    if ( flat(i,j) == WALL ) then
      string(ilo:ilo+3) = '+--+'
    else if ( flat(i,j) == OPEN ) then
      string(ilo:ilo+3) = '+  +'
    else if ( flat(i,j) == PATH ) then
      string(ilo:ilo+3) = '+**+'
    else if ( flat(i,j) == INDEF ) then

    end if

    ilo = ilo + 3

  end do

  write ( *, '(a)' ) string(1:ilo)

  return
end
subroutine maze_random ( m, n, bar, dad, flat )

!*****************************************************************************80
!
!! maze_random() generates a random maze in a rectangular region.
!
!  Discussion:
!
!    The rectangular region is assumed to be made of a grid of M rows
!    and N columns of square cells.  The maze is to be begun in 
!    one cell, and ended in another.  The boundary of the region
!    is walled off, except possibly for entrances to the beginning
!    cell, and an exit from the ending cell.
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
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns of cells.
!
!    Output, integer BAR(M,N+1), records the vertical "bars" 
!    in the maze.
!    -1, means "indefinite", that there is no cell of the maze on either
!        side of this position;
!     0, means "wall", that there is a cell on at least one side, and
!        a wall here;
!     1, means "open", that there are cells on both sides (or possibly
!        an opening to the exterior) and the way is open.
!
!    Output, integer DAD(M,N), a rooted tree representation of
!    the maze.  The root of the tree has DAD(I,J) = 0.  All other cells
!    that are connectable to the root should have DAD(I,J) = K, where
!    K is the cell index K = ( J - 1 ) * M + I, with I and J the row
!    and column indices of the cell.  If the cell is not connectable
!    to the root, then DAD(I,J) is -1.
!
!    Output, integer FLAT(M+1,N), records the horizontal "flats" 
!    in the maze.
!    -1, means "indefinite", that there is no cell of the maze on either
!        side of this position;
!     0, means "wall", that there is a cell on at least one side, and
!        a wall here;
!     1, means "open", that there are cells on both sides (or possibly
!        an opening to the exterior) and the way is open.
!
  implicit none

  integer, parameter :: maxstack = 500
  integer, parameter :: NORTH = 1
  integer, parameter :: EAST = 2
  integer, parameter :: WEST = 3
  integer, parameter :: SOUTH = 4
  integer, parameter :: INDEF = -1
  integer, parameter :: WALL = 0
  integer, parameter :: OPEN = 1

  integer m
  integer n

  integer bar(m,n+1)
  integer dad(m,n)
  integer dir
  integer flat(m+1,n)
  integer i
  integer i4_uniform_ab
  integer ihi
  integer ival
  integer j
  integer jval
  integer nabe
  integer nbase
  integer stack(maxstack)

  if ( m < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'maze_random(): Fatal error!'
    write ( *, '(a)' ) '  M must be at least 1.'
    stop 1
  end if

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'maze_random(): Fatal error!'
    write ( *, '(a)' ) '  N must be at least 1.'
    stop 1
  end if

  if ( m == 1 .and. n == 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'maze_random(): Fatal error!'
    write ( *, '(a)' ) '  At least one of M and N must be more than 1.'
    stop 1
  end if
!
!  Initialize arrays to INDEF.
!
  bar(1:m,1:n+1) = INDEF
  flat(1:m+1,1:n) = INDEF
!
!  Set the boundaries to walls.
!
  flat(1,1:n) = WALL
  flat(m+1,1:n) = WALL
  bar(1:m,1) = WALL
  bar(1:m,n+1) = WALL
!
!  Initialize the tree pointers.
!
  dad(1:m,1:n) = -1
!
!  Pick a random starting point.
!
  ival = i4_uniform_ab ( 1, m )
  jval = i4_uniform_ab ( 1, n )

  dad(ival,jval) = 0
!
!  Count the number of neighbors of the starting cell,
!  choose randomly from the neigbors, and add it.
!
  do

    nabe = 0

    do i = 1, m
      do j = 1, n

        if ( dad(i,j) /= -1 ) then

          if ( flat(i,j) /= WALL ) then

            if ( dad(i-1,j) == -1 ) then

              if ( maxstack < nabe + 3 ) then
                write ( *, '(a)' ) ' '
                write ( *, '(a)' ) 'maze_random(): Fatal error!'
                write ( *, '(a)' ) '  Ran out of stack space.'
                stop 1
              end if

              stack(nabe+1) = i
              stack(nabe+2) = j
              stack(nabe+3) = NORTH
              nabe = nabe + 3

            end if

          end if

          if ( bar(i,j+1) /= WALL ) then

            if ( dad(i,j+1) == -1 ) then

              if ( maxstack < nabe + 3 ) then
                write ( *, '(a)' ) ' '
                write ( *, '(a)' ) 'maze_random(): Fatal error!'
                write ( *, '(a)' ) '  Ran out of stack space.'
                stop 1
              end if

              stack(nabe+1) = i
              stack(nabe+2) = j
              stack(nabe+3) = EAST
              nabe = nabe + 3

            end if

          end if

          if ( bar(i,j) /= WALL ) then

            if ( dad(i,j-1) == -1 ) then

              if ( maxstack < nabe + 3 ) then
                write ( *, '(a)' ) ' '
                write ( *, '(a)' ) 'maze_random(): Fatal error!'
                write ( *, '(a)' ) '  Ran out of stack space.'
                stop 1
              end if

              stack(nabe+1) = i
              stack(nabe+2) = j
              stack(nabe+3) = WEST
              nabe = nabe + 3

            end if

          end if

          if ( flat(i+1,j) /= WALL ) then

            if ( dad(i+1,j) == -1 ) then

              if ( maxstack < nabe + 3 ) then
                write ( *, '(a)' ) ' '
                write ( *, '(a)' ) 'maze_random(): Fatal error!'
                write ( *, '(a)' ) '  Ran out of stack space!'
                stop 1
              end if

              stack(nabe+1) = i
              stack(nabe+2) = j
              stack(nabe+3) = SOUTH
              nabe = nabe + 3

            end if

          end if

        end if

      end do
    end do
!
!  If there are accessible neighbors, randomly choose one.
!
    if ( nabe <= 0 ) then
      exit
    end if

    ihi = nabe / 3
 
    ival = i4_uniform_ab ( 1, ihi )
  
    nbase = 3 * ival - 3
    i = stack(nbase+1)
    j = stack(nbase+2)
    dir = stack(nbase+3)

    if ( dir == NORTH ) then
      flat(i,j) = OPEN
      dad(i-1,j) = ( j - 1 ) * m + i
    else if ( dir == EAST ) then
      bar(i,j+1) = OPEN
      dad(i,j+1) = ( j - 1 ) * m + i
    else if ( dir == WEST ) then
      bar(i,j) = OPEN
      dad(i,j-1) = ( j - 1 ) * m + i
    else if ( dir == SOUTH ) then
      flat(i+1,j) = OPEN
      dad(i+1,j) = ( j - 1 ) * m + i
    end if

  end do
!
!  Set all remaining INDEF's to WALLS.
!
  do i = 1, m
    do j = 1, n + 1
      if ( bar(i,j) == INDEF ) then
        bar(i,j) = WALL
      end if
    end do
  end do

  do i = 1, m + 1
    do j = 1, n
      if ( flat(i,j) == INDEF ) then
        flat(i,j) = WALL
      end if
    end do
  end do

  return
end
