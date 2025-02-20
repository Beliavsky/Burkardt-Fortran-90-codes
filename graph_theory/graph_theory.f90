subroutine balanc ( nm, n, a, low, igh, scale )

!*****************************************************************************80
!
!! balanc() balances a real matrix before eigenvalue calculations.
!
!  Discussion:
!
!    This subroutine balances a real matrix and isolates eigenvalues
!    whenever possible.
!
!    Suppose that the principal submatrix in rows LOW through IGH
!    has been balanced, that P(J) denotes the index interchanged
!    with J during the permutation step, and that the elements
!    of the diagonal matrix used are denoted by D(I,J).  Then
!
!      SCALE(J) = P(J),    J = 1,...,LOW-1,
!               = D(J,J),  J = LOW,...,IGH,
!               = P(J)     J = IGH+1,...,N.
!
!    The order in which the interchanges are made is N to IGH+1,
!    then 1 to LOW-1.
!
!    Note that 1 is returned for LOW if IGH is zero formally.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 December 2008
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, 
!    Ikebe, Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow, 
!    Y Ikebe, V Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of A, which must
!    be at least N.
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real ( kind = rk ) A(NM,N), the N by N matrix.  On output,
!    the matrix has been balanced.
!
!    Output, integer LOW, IGH, indicate that A(I,J) is equal to 
!    zero if
!    (1) I is greater than J and
!    (2) J=1,...,LOW-1 or I=IGH+1,...,N.
!
!    Output, real ( kind = rk ) SCALE(N), contains information determining the
!    permutations and scaling factors used.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nm
  integer n

  real ( kind = rk ) a(nm,n)
  real ( kind = rk ) b2
  real ( kind = rk ) c
  real ( kind = rk ) f
  real ( kind = rk ) g
  integer i
  integer iexc
  integer igh
  integer j
  integer k
  integer l
  integer low
  integer m
  logical noconv
  real ( kind = rk ) r
  real ( kind = rk ), parameter :: radix = 16.0D+00
  real ( kind = rk ) s
  real ( kind = rk ) scale(n)

  iexc = 0
  j = 0
  m = 0

  b2 = radix * radix
  k = 1
  l = n
  go to 100

20 continue

  scale(m) = j

  if ( j /= m ) then

    do i = 1, l
      call r8_swap ( a(i,j), a(i,m) )
    end do

    do i = k, n
      call r8_swap ( a(j,i), a(m,i) )
    end do

  end if

  if ( iexc == 2 ) then
    go to 130
  end if
!
!  Search for rows isolating an eigenvalue and push them down.
!
  if ( l == 1 ) then
    low = k
    igh = l
    return
  end if

  l = l - 1

100 continue

  do j = l, 1, -1

     do i = 1, l
       if ( i /= j ) then
         if ( a(j,i) /= 0.0D+00 ) then
           go to 120
         end if
       end if
     end do

     m = l
     iexc = 1
     go to 20

120  continue

  end do

  go to 140
!
!  Search for columns isolating an eigenvalue and push them left.
!
130 continue

  k = k + 1

140 continue

  do j = k, l

    do i = k, l
      if ( i /= j ) then
        if ( a(i,j) /= 0.0D+00 ) then
          go to 170
        end if
      end if
    end do

    m = k
    iexc = 2
    go to 20

170 continue

  end do
!
!  Balance the submatrix in rows K to L.
!
  scale(k:l) = 1.0D+00
!
!  Iterative loop for norm reduction.
!
  noconv = .true.

  do while ( noconv )

    noconv = .false.

    do i = k, l

      c = 0.0D+00
      r = 0.0D+00

      do j = k, l
        if ( j /= i ) then
          c = c + abs ( a(j,i) )
          r = r + abs ( a(i,j) )
        end if
      end do
!
!  Guard against zero C or R due to underflow.
!
      if ( c /= 0.0D+00 .and. r /= 0.0D+00 ) then

        g = r / radix
        f = 1.0D+00
        s = c + r

        do while ( c < g )
          f = f * radix
          c = c * b2
        end do

        g = r * radix

        do while ( g <= c )
          f = f / radix
          c = c / b2
        end do
!
!  Balance.
!
        if ( ( c + r ) / f < 0.95D+00 * s ) then

          g = 1.0D+00 / f
          scale(i) = scale(i) * f
          noconv = .true.

          a(i,k:n) = a(i,k:n) * g
          a(1:l,i) = a(1:l,i) * f

        end if

      end if

    end do

  end do

  low = k
  igh = l

  return
end
subroutine ch_cap ( c )

!*****************************************************************************80
!
!! ch_cap() capitalizes a single character.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, the character to capitalize.
!
  implicit none

  character c
  integer itemp

  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end
function ch_eqi ( c1, c2 )

!*****************************************************************************80
!
!! ch_eqi() is a case insensitive comparison of two characters for equality.
!
!  Example:
!
!    C_EQI ( 'A', 'a' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C1, C2, the characters to compare.
!
!    Output, logical CH_EQI, the result of the comparison.
!
  implicit none

  logical ch_eqi
  character c1
  character c2
  character cc1
  character cc2

  cc1 = c1
  cc2 = c2

  call ch_cap ( cc1 )
  call ch_cap ( cc2 )

  if ( cc1 == cc2 ) then
    ch_eqi = .true.
  else
    ch_eqi = .false.
  end if

  return
end
subroutine ch_to_digit ( c, digit )

!*****************************************************************************80
!
!! ch_to_digit() returns the integer value of a base 10 digit.
!
!  Example:
!
!     C   DIGIT
!    ---  -----
!    '0'    0
!    '1'    1
!    ...  ...
!    '9'    9
!    ' '    0
!    'X'   -1
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the decimal digit, '0' through '9' or blank
!    are legal.
!
!    Output, integer DIGIT, the corresponding integer value.  
!    If C was 'illegal', then DIGIT is -1.
!
  implicit none

  character c
  integer digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

    digit = ichar ( c ) - 48

  else if ( c == ' ' ) then

    digit = 0

  else

    digit = - 1

  end if

  return
end
subroutine catalan ( n, c )

!*****************************************************************************80
!
!! catalan() computes the Catalan numbers, from C(0) to C(N).
!
!  First values:
!
!     C(0)     1
!     C(1)     1
!     C(2)     2
!     C(3)     5
!     C(4)    14
!     C(5)    42
!     C(6)   132
!     C(7)   429
!     C(8)  1430
!     C(9)  4862
!    C(10) 16796
!
!  Formula:
!
!    C(N) = (2*N)! / ( (N+1) * (N!) * (N!) )
!         = 1 / (N+1) * COMB ( 2N, N )
!         = 1 / (2N+1) * COMB ( 2N+1, N+1).
!
!  Recursion:
!
!    C(N) = 2 * (2*N-1) * C(N-1) / (N+1)
!    C(N) = SUM ( I = 1 to N-1 ) C(I) * C(N-I)
!
!  Comments:
!
!    The Catalan number C(N) counts:
!
!    1) the number of binary trees on N vertices;
!    2) the number of ordered trees on N+1 vertices;
!    3) the number of full binary trees on 2N+1 vertices;
!    4) the number of well formed sequences of 2N parentheses;
!    5) number of ways 2N ballots can be counted, in order,
!       with N positive and N negative, so that the running sum
!       is never negative;
!    6) the number of standard tableaus in a 2 by N rectangular Ferrers diagram;
!    7) the number of monotone functions from [1..N} to [1..N} which
!       satisfy f(i) <= i for all i,
!    8) the number of ways to triangulate a polygon with N+2 vertices.
!
!  Example:
!
!    N = 3
!
!    ()()()
!    ()(())
!    (()())
!    (())()
!    ((()))
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dennis Stanton, Dennis White,
!    Constructive Combinatorics,
!    Springer Verlag, New York, 1986.
!
!  Parameters:
!
!    Input, integer N, the number of Catalan numbers desired.
!
!    Output, integer C(0:N), the Catalan numbers from C(0) to C(N).
!
  implicit none

  integer n

  integer i
  integer c(0:n)

  c(0) = 1
!
!  The extra parentheses ensure that the integer division is
!  done AFTER the integer multiplication.
!
  do i = 1, n
    c(i) = ( c(i-1) * 2 * ( 2 * i - 1 ) ) / ( i + 1 )
  end do

  return
end
subroutine degree_sequence_is_graphic ( nnode, seq, result )

!*****************************************************************************80
!
!! degree_sequence_is_graphic() reports whether a degree sequence represents a graph.
!
!  Discussion:
!
!    The degree sequence of a graph is constructed by computing the
!    degree of each node, and then ordering these values in decreasing order.
!
!    A sequence of NNODE nonnegative integers is said to be "graphic" if
!    there exists a graph whose degree sequence is the given sequence.
!
!    The Havel Hakimi theorem states that 
!
!      s t1 t2 ... ts d1 d2 ... dn
!
!    is graphic if and only if
!
!        t1-1 t2-1 ... ts-1 d1 d2 ... dn
!
!    is graphic (after any necessary resorting and dropping of zeroes).
!    Definitely, the one thing we cannot have is that any nonzero entry
!    is equal to or greater than the number of nonzero entries.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, integer SEQ(NNODE), the degree sequence to be tested.
!
!    Output, integer RESULT, the result.
!    0, SEQ is not graphic.
!    1, SEQ is graphic.
!
  implicit none

  integer nnode

  integer dmax
  integer i
  integer i4vec_nonzero
  integer nonzero
  integer order
  integer result
  integer seq(nnode)

  result = 0

  do i = 1, nnode
    if ( seq(i) < 0 ) then
      return
    end if
  end do
!
!  Check that the sequence is decreasing.
!
  call i4vec_order_type ( nnode, seq, order )

  if ( order == -1 .or. order == 1 .or. order == 2 ) then
    return
  end if
!
!  Now apply the Havel Hakimi theorem.
!
  do

    nonzero = i4vec_nonzero ( nnode, seq )

    if ( nonzero == 0 ) then
      result = 1
      exit
    end if

    call i4vec_sort_heap_d ( nnode, seq )

    dmax = seq(1)

    if ( nonzero <= dmax ) then
      result = 0
      exit
    end if

    seq(1) = 0
    do i = 2, dmax + 1
      seq(i) = seq(i) - 1
    end do

  end do
        
  return
end
subroutine degree_sequence_to_graph_adj ( nnode, seq, adj, ierror )

!*****************************************************************************80
!
!! degree_sequence_to_graph_adj() computes a graph with the given degree sequence.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, integer SEQ(NNODE), the degree sequence.
!
!    Output, integer ADJ(NNODE,NNODE), the adjacency information.  
!    ADJ(I,J) is nonzero if there is an edge from node I to node J.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  implicit none

  integer nnode

  integer adj(nnode,nnode)
  integer i
  integer ierror
  integer indx(nnode)
  integer nonzero
  integer s
  integer seq(nnode)
  integer seq2(nnode)

  ierror = 0

  adj(1:nnode,1:nnode) = 0

  seq2(1:nnode) = seq(1:nnode)

  do

    call i4vec_sort_heap_index_d ( nnode, seq2, indx )

    nonzero = 0
    do i = 1, nnode
      if ( seq2(i) /= 0 ) then
        nonzero = nonzero + 1
      end if
    end do

    if ( nonzero == 0 ) then
      exit
    end if

    s = seq2(indx(1))

    if ( nonzero <= s ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'degree_sequence_to_graph_adj(): Fatal error!'
      write ( *, '(a)' ) '  The degree sequence is not graphic!'
      return
    end if

    seq2(indx(1)) = 0

    do i = 2, s+1
      adj(indx(i),indx(1)) = 1
      adj(indx(1),indx(i)) = 1
      seq2(indx(i)) = seq2(indx(i)) - 1
    end do

  end do

  return
end
subroutine dge_check ( m, n, ierror )

!*****************************************************************************80
!
!! dge_check() checks the dimensions of a general matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
!
!    Output, integer IERROR, reports whether any errors 
!    were detected.
!    IERROR is set to 0 before the checks are made, and then:
!    IERROR = IERROR + 2 if M is illegal;
!    IERROR = IERROR + 4 if N is illegal.
!
  implicit none

  integer ierror
  integer m
  integer n

  ierror = 0

  if ( m < 1 ) then
    ierror = ierror + 2
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) 'DGE_CHECK - Illegal M = ', m
  end if

  if ( n < 1 ) then
    ierror = ierror + 4
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) 'DGE_CHECK - Illegal N = ', n
  end if

  return
end
subroutine dge_det ( n, a, ipivot, det )

!*****************************************************************************80
!
!! dge_det() computes the determinant of a matrix factored by DGE_FA or DGE_TRF.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real ( kind = rk ) A(N,N), the LU factors computed 
!    by DGE_FA or DGE_TRF.
!
!    Input, integer IPIVOT(N), as computed by DGE_FA or DGE_TRF.
!
!    Output, real ( kind = rk ) DET, the determinant of the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) det
  integer i
  integer ierror
  integer ipivot(n)
!
!  Check the dimensions.
!
  call dge_check ( n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DGE_DET(): Fatal error!'
    write ( *, '(a)' ) '  Illegal dimensions.'
    return
  end if

  det = 1.0D+00

  do i = 1, n
    det = det * a(i,i)
  end do

  do i = 1, n
    if ( ipivot(i) /= i ) then
      det = - det
    end if
  end do

  return
end
subroutine dge_fa ( n, a, ipivot, info )

!*****************************************************************************80
!
!! dge_fa() factors a general matrix.
!
!  Discussion:
!
!    DGE_FA is a simplified version of the LINPACK routine DGEFA.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input/output, real ( kind = rk ) A(N,N), the matrix to be factored.
!    On output, A contains an upper triangular matrix and the multipliers
!    which were used to obtain it.  The factorization can be written
!    A = L * U, where L is a product of permutation and unit lower
!    triangular matrices and U is upper triangular.
!
!    Output, integer IPIVOT(N), a vector of pivot indices.
!
!    Output, integer INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,n)
  integer i
  integer ierror
  integer info
  integer ipivot(n)
  integer j
  integer k
  integer l
!
!  Check the dimensions.
!
  call dge_check ( n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DGE_FA(): Fatal error!'
    write ( *, '(a)' ) '  Illegal dimensions.'
    return
  end if

  info = 0

  do k = 1, n-1
!
!  Find L, the index of the pivot row.
!
    l = k
    do i = k+1, n
      if ( abs ( a(l,k) ) < abs ( a(i,k) ) ) then
        l = i
      end if
    end do

    ipivot(k) = l
!
!  If the pivot index is zero, the algorithm has failed.
!
    if ( a(l,k) == 0.0D+00 ) then
      info = k
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DGE_FA(): Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      return
    end if
!
!  Interchange rows L and K if necessary.
!
    if ( l /= k ) then
      call r8_swap ( a(l,k), a(k,k) )
    end if
!
!  Normalize the values that lie below the pivot entry A(K,K).
!
    a(k+1:n,k) = -a(k+1:n,k) / a(k,k)
!
!  Row elimination with column indexing.
!
    do j = k+1, n

      if ( l /= k ) then
        call r8_swap ( a(l,j), a(k,j) )
      end if

      a(k+1:n,j) = a(k+1:n,j) + a(k+1:n,k) * a(k,j)

    end do

  end do

  ipivot(n) = n

  if ( a(n,n) == 0.0D+00 ) then
    info = n
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DGE_FA(): Fatal error!'
    write ( *, '(a,i8)' ) '  Zero pivot on step ', info
  end if

  return
end
subroutine digraph_adj_print ( adj, nnode, title )

!*****************************************************************************80
!
!! digraph_adj_print() prints out an adjacency matrix for a digraph.
!
!  Discussion:
!
!    This routine actually allows the entries of ADJ to have ANY value.
!    Values between 0 and 9 will be printed as is.  Other values will
!    be printed as '*'.
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
!    Input, integer ADJ(NNODE,NNODE), the adjacency matrix of a 
!    digraph.  ADJ(I,J) is 1 if there is a direct connection FROM node I TO 
!    node J, and is 0 otherwise.
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
  integer jhi
  character ( len = 80 ) string
  character ( len = * ) title

  if ( len_trim ( title ) /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  do i = 1, nnode

    jhi = min ( nnode, 80 )

    do j = 1, jhi

      if ( 0 <= adj(i,j) .and. adj(i,j) <= 9 ) then
        string(j:j) = char ( 48 + adj(i,j) )
      else
        string(j:j) = '*'
      end if

    end do

    write ( *, '(i2,2x,a)' ) i, string(1:jhi)

  end do

  return
end
subroutine digraph_arc_edge_sort ( nedge, inode, jnode )

!*****************************************************************************80
!
!! digraph_arc_edge_sort() sorts the edge array of a graph.
!
!  Discussion:
!
!    The edges are sorted in dictionary order.
!
!  Example:
!
!    Input:
!
!      INODE  JNODE
!
!        3      2
!        2      4
!        4      3
!        2      1
!        1      4
!
!    Output:
!
!      INODE  JNODE
!
!        1      4
!        2      1
!        2      4
!        3      2
!        4      3
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    24 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NEDGE, the number of edges.
!
!    Input/output, integer INODE(NEDGE), JNODE(NEDGE), the edge
!    array.  The I-th edge goes from node INODE(I) to node JNODE(I).
!    On output, the INODE and JNODE arrays have been sorted as described.
!
  implicit none

  integer nedge

  integer iedge
  integer indx
  integer inode(nedge)
  integer isgn
  integer jedge
  integer jnode(nedge)

  if ( nedge <= 1 ) then
    return
  end if
!
!  Sort the edges using an external heap sort.
!
  iedge = 0
  jedge = 0
  indx = 0
  isgn = 0

  do

    call sort_heap_external ( nedge, indx, iedge, jedge, isgn )
!
!  Interchange edges IEDGE and JEDGE.
!
    if ( 0 < indx ) then

      call i4_swap ( inode(iedge), inode(jedge) )
      call i4_swap ( jnode(iedge), jnode(jedge) )
!
!  Compare edges IEDGE and JEDGE.
!
    else if ( indx < 0 ) then

      if ( ( inode(iedge) < inode(jedge) ) .or. &
        ( inode(iedge) == inode(jedge) .and. &
          jnode(iedge) < jnode(jedge) ) ) then
        isgn = -1
      else
        isgn = +1
      end if

    else if ( indx == 0 ) then

      exit

    end if

  end do
 
  return
end
subroutine digraph_arc_print ( nedge, inode, jnode, title )

!*****************************************************************************80
!
!! digraph_arc_print() prints a digraph from an edge list.
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
!    Input, integer NEDGE, the number of edges.
!
!    Input, integer INODE(NEDGE), JNODE(NEDGE), the beginning and
!    end nodes of the edges.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer nedge

  integer i
  integer inode(nedge)
  integer jnode(nedge)
  character ( len = * ) title

  if ( len_trim ( title ) /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  do i = 1, nedge
    write ( *, '(i8,4x,2i8)' ) i, inode(i), jnode(i)
  end do

  return
end
subroutine digraph_dist_print ( dist, nnode, title )

!*****************************************************************************80
!
!! digraph_dist_print() prints the distance matrix defining a digraph.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 July 2000
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
subroutine digraph_inc_print ( nnode, narc, inc, title )

!*****************************************************************************80
!
!! digraph_inc_print() prints the incidence matrix of a digraph.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, integer NARC, the number of arcs.
!
!    Input, integer INC(NNODE,NARC), the incidence matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer narc
  integer nnode

  integer i
  integer inc(nnode,narc)
  character ( len = * ) title

  if ( len_trim ( title ) /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  do i = 1, nnode
    write ( *, '(20i3)' ) inc(i,1:narc)
  end do

  return
end
subroutine edge_add_nodes ( edge, max_edge, num_edge, iface, n1, n2, ierror )

!*****************************************************************************80
!
!! edge_add_nodes() adds the edge defined by two nodes to the edge list.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer EDGE(4,MAX_EDGE), edge information.
!    EDGE(1,I) is the starting node of edge I;
!    EDGE(2,I) is the ending node of edge I;
!    EDGE(3,I) is the positive face;
!    EDGE(4,I) is the negative face, if any.
!
!    Input, integer MAX_EDGE, the maximum number of edges.
!
!    Input/output, integer NUM_EDGE, the number of edges.
!
!    Input, integer IFACE, the face to which the nodes belong.
!
!    Input, integer N1, N2, two nodes which form an edge.
!
!    Output, integer IERROR, error flag, 0 = no error, 
!    nonzero = error.
!
  implicit none

  integer max_edge

  integer edge(4,max_edge)
  integer ierror
  integer iface
  integer n1
  integer n2
  integer num_edge

  if ( num_edge < max_edge ) then
    num_edge = num_edge + 1
    edge(1,num_edge) = n1
    edge(2,num_edge) = n2
    edge(3,num_edge) = iface
    edge(4,num_edge) = 0
    ierror = 0
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EDGE_ADD_NODES(): Fatal error!'
    write ( *, '(a,i8)' ) '  Exceeding MAX_EDGE = ', max_edge
    ierror = 1
  end if

  return
end
subroutine edge_bound ( edge, max_edge, num_edge )

!*****************************************************************************80
!
!! edge_bound() reports the edges which are part of the boundary.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer EDGE(4,MAX_EDGE), edge information.
!    EDGE(1,I) is the starting node of edge I;
!    EDGE(2,I) is the ending node of edge I;
!    EDGE(3,I) is the positive face;
!    EDGE(4,I) is the negative face, if any.
!
!    Input, integer MAX_EDGE, the maximum number of edges.
!
!    Input, integer NUM_EDGE, the number of edges.
!
  implicit none

  integer max_edge

  integer edge(4,max_edge)
  integer iedge
  integer num_bound
  integer num_edge

  num_bound = 0

  do iedge = 1, num_edge
    if ( ( edge(3,iedge) /= 0 .and. edge(4,iedge) == 0 ) .or. &
         ( edge(3,iedge) == 0 .and. edge(4,iedge) /= 0 ) ) then
      num_bound = num_bound + 1
    end if
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EDGE_BOUND'
  write ( *, '(a,i8)' ) '  Number of boundary edges = ', num_bound

  return
end
subroutine edge_match_face ( edge, max_edge, num_edge, facelist, n, index )

!*****************************************************************************80
!
!! edge_match_face() seeks an edge common to a face and the edge list.
!
!  Discussion:
!
!    If a common edge is found, then the information in the face node
!    list is adjusted so that the first two entries correspond to the
!    matching edge in EDGE, but in reverse order.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer EDGE(4,MAX_EDGE), edge information.
!    EDGE(1,I) is the starting node of edge I;
!    EDGE(2,I) is the ending node of edge I;
!    EDGE(3,I) is the positive face;
!    EDGE(4,I) is the negative face, if any.
!
!    Input, integer MAX_EDGE, the maximum number of edges.
!
!    Input, integer NUM_EDGE, the number of edges.
!
!    Input/output, integer FACELIST(N), the list of nodes making a
!    face.
!
!    Input, integer N, the number of nodes in the face.
!
!    Output, integer INDEX, the results of the search.
!    0, there is no edge common to the face and the EDGE array.
!    nonzero, edge INDEX is common to the face and the EDGE array.
!
  implicit none

  integer n
  integer max_edge

  integer edge(4,max_edge)
  integer facelist(n)
  integer iedge
  integer index
  integer j
  integer jp1
  integer n1
  integer n2
  integer num_edge

  index = 0

  if ( n <= 0 ) then
    return
  end if

  if ( num_edge <= 0 ) then
    return
  end if

  do j = 1, n

    if ( j == n ) then
      jp1 = 1
    else
      jp1 = j + 1
    end if

    n1 = facelist(j)
    n2 = facelist(jp1)

    do iedge = 1, num_edge

      if ( edge(1,iedge) == n2 .and. edge(2,iedge) == n1 ) then

        call i4vec_rotate ( n, 1 - j, facelist )

        index = iedge
        return

      else if ( edge(1,iedge) == n1 .and. edge(2,iedge) == n2 ) then

        call i4vec_rotate ( n, n - jp1, facelist )

        call i4vec_reverse ( n, facelist )

        index = iedge
        return

      end if

    end do
   
  end do

  return
end
subroutine edge_match_nodes ( edge, max_edge, num_edge, n1, n2, iedge )

!*****************************************************************************80
!
!! edge_match_nodes() seeks an edge of the form (N1,N2) or (N2,N1) in EDGE.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer EDGE(4,MAX_EDGE), edge information.
!    EDGE(1,I) is the starting node of edge I;
!    EDGE(2,I) is the ending node of edge I;
!    EDGE(3,I) is the positive face;
!    EDGE(4,I) is the negative face, if any.
!
!    Input, integer MAX_EDGE, the maximum number of edges.
!
!    Input, integer NUM_EDGE, the number of edges.
!
!    Input, integer N1, N2, two nodes that form an edge.
!
!    Output, integer IEDGE, the results of the search.
!    0, no matching edge was found.
!    nonzero, edge IEDGE of the EDGE array matches (N1,N2) or (N2,N1).
!
  implicit none

  integer max_edge

  integer edge(4,max_edge)
  integer i
  integer iedge
  integer n1
  integer n2
  integer num_edge

  iedge = 0
  do i = 1, num_edge

    if ( ( n1 == edge(1,i) .and. n2 == edge(2,i) ) .or. &
         ( n2 == edge(1,i) .and. n1 == edge(2,i) ) ) then
      iedge = i
      return
    end if

  end do

  return
end
subroutine edges_to_ps ( plotxmin2, plotymin2, alpha, iunit, inode, jnode, &
  nedge, nnode, x, y, xmin, ymin )

!*****************************************************************************80
!
!! edges_to_ps() writes subplot edges to a PostScript file.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer PLOTXMIN2, PLOTYMIN2, the Postscript origin.
!
!    Input, real ( kind = rk ) ALPHA, the physical-to-Postscript scale factor.
!
!    Input, integer IUNIT, the output FORTRAN unit.
!
!    Input, integer INODE(NEDGE), JNODE(NEDGE), the edge array.
!    The I-th edge connects nodes INODE(I) and JNODE(I).
!
!    Input, integer NEDGE, the number of edges.
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, real ( kind = rk ) X(NNODE), Y(NNODE), the X and Y components
!    of points.
!
!    Input, real ( kind = rk ) XMIN, YMIN, the physical origin.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nedge
  integer nnode

  real ( kind = rk ) alpha
  integer i
  integer inode(nedge)
  integer iunit
  integer jnode(nedge)
  integer node
  integer plotxmin2
  integer plotymin2
  integer px1
  integer px2
  integer py1
  integer py2
  real ( kind = rk ) x(nnode)
  real ( kind = rk ) xmin
  real ( kind = rk ) y(nnode)
  real ( kind = rk ) ymin
!
!  Draw lines.
!
  do i = 1, nedge

    node = inode(i)
    px1 = plotxmin2 + nint ( alpha * ( x(node) - xmin ) )
    py1 = plotymin2 + nint ( alpha * ( y(node) - ymin ) )

    node = jnode(i)
    px2 = plotxmin2 + nint ( alpha * ( x(node) - xmin ) )
    py2 = plotymin2 + nint ( alpha * ( y(node) - ymin ) )

    write ( iunit, '(2i4,a,2i4,a)' ) px1, py1, ' moveto ', px2, py2, &
      ' lineto stroke'

  end do

  return
end
subroutine elmhes ( nm, n, low, igh, a, ind )

!*****************************************************************************80
!
!! elmhes() transforms a real general matrix to upper Hessenberg form.
!
!  Discussion:
!
!    Given a real general matrix, this subroutine reduces a submatrix
!    situated in rows and columns LOW through IGH to upper Hessenberg
!    form by stabilized elementary similarity transformations.
!
!  Reference:
!
!    Martin, James Wilkinson,
!    ELMHES,
!    Numerische Mathematik,
!    Volume 12, pages 349-368, 1968.
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow, 
!    Y Ikebe, V Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the array A.
!    NM must be at least N.
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer LOW, IGH, are determined by the balancing
!    routine BALANC.  If BALANC has not been used, set LOW = 1, IGH = N.
!
!    Input/output, real ( kind = rk ) A(NM,N).  On input, the matrix to be
!    reduced.  On output, the Hessenberg matrix.  The multipliers
!    which were used in the reduction are stored in the
!    remaining triangle under the Hessenberg matrix.
!
!    Output, integer IND(N), contains information on the rows and
!    columns interchanged in the reduction.  Only elements LOW through IGH are
!    used.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer igh
  integer n
  integer nm

  real ( kind = rk ) a(nm,n)
  integer i
  integer ind(igh)
  integer j
  integer la
  integer low
  integer m
  integer mm1
  real ( kind = rk ) x
  real ( kind = rk ) y

  la = igh - 1

  do m = low + 1, la

    mm1 = m - 1
    x = 0.0D+00
    i = m

    do j = m, igh
      if ( abs ( x ) < abs ( a(j,mm1) ) ) then
        x = a(j,mm1)
        i = j
      end if
    end do

    ind(m) = i
!
!  Interchange rows and columns of the matrix.
!
    if ( i /= m ) then

      do j = mm1, n
        call r8_swap ( a(i,j), a(m,j) )
      end do

      do j = 1, igh
        call r8_swap ( a(j,i), a(j,m) )
      end do

    end if

    if ( x /= 0.0D+00 ) then

      do i = m+1, igh

        y = a(i,mm1)

        if ( y /= 0.0D+00 ) then

          y = y / x
          a(i,mm1) = y

          do j = m, n
            a(i,j) = a(i,j) - y * a(m,j)
          end do

          do j = 1, igh
            a(j,m) = a(j,m) + y * a(j,i)
          end do

        end if

      end do

    end if

  end do

  return
end
subroutine face_check ( edge, face, face_object, face_order, face_rank, &
  face_tier, max_edge, max_order, num_edge, num_face, num_object )

!*****************************************************************************80
!
!! face_check() checks and analyzes a set of faces.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer EDGE(4,MAX_EDGE), edge information.
!    EDGE(1,I) is the starting node of edge I;
!    EDGE(2,I) is the ending node of edge I;
!    EDGE(3,I) is the positive face;
!    EDGE(4,I) is the negative face, or 0 if the edge is used once.
!
!    Input, integer FACE(MAX_ORDER,NUM_FACE), describes the faces.
!    FACE(I,J) is the index of the I-th node in face J.  It is best
!    if the nodes of all faces are listed in counterclockwise order.
!
!    Output, integer FACE_OBJECT(NUM_FACE), describes the objects.
!    FACE_OBJECT(I) is the index of the edge-connected "object" to 
!    which face I belongs.
!
!    Input, integer FACE_ORDER(NUM_FACE), is the number of nodes
!    making up each face.
!
!    Output, integer FACE_RANK(NUM_FACE), is an ordered list of
!    faces.  FACE_RANK(1) is the index of the face in the first tier of the 
!    first object, followed by second tier faces, and so on until
!    object one is complete.  Object two follows, and so on.
!
!    Output, integer FACE_TIER(NUM_FACE).  FACE_TIER(I) is the
!    "tier" of face I in its object.  The base of the object has tier 1,
!    the neighbors of the base have tier 2, and so on.
!
!    Input, integer MAX_EDGE, the maximum number of edges.
!
!    Input, integer MAX_ORDER, is the maximum number of nodes that
!    can make up a face, required to dimension FACE.
!
!    Output, integer NUM_EDGE, the number of edges.
!
!    Input, integer NUM_FACE, the number of faces.
!
!    Output, integer NUM_OBJECT, the number of objects.
!
  implicit none

  integer max_edge
  integer max_order
  integer num_face

  integer edge(4,max_edge)
  integer face(max_order,num_face)
  integer face_object(num_face)
  integer face_order(num_face)
  integer face_rank(num_face)
  integer face_tier(num_face)
  integer i
  integer ierror
  integer j
  integer num_edge
  integer num_fix
  integer num_object
!
!  Organize the faces into layered objects.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Determine edge-connected objects.'

  call object_build ( face, face_object, face_order, face_rank, face_tier, &
    max_order, num_face, num_object )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) 'Number of objects = ', num_object
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Face, Object, Tier'
  write ( *, '(a)' ) ' '

  do i = 1, num_face
    write ( *, '(3i8)' ) i, face_object(i), face_tier(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Preferred order:'
  write ( *, '(a)' ) '  Order, Face'
  write ( *, '(a)' ) ' '
  do i = 1, num_face
    write ( *, '(2i8)' ) i, face_rank(i)
  end do
!
!  Reorder the faces by object and tier.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Reorder the faces.'

  call face_sort ( face, face_object, face_order, face_tier, max_order, &
    num_face )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Face, Label, Object, Tier'
  write ( *, '(a)' ) ' '

  do i = 1, num_face
    write ( *, '(4i8)' ) i, face_rank(i), face_object(i), face_tier(i)
  end do
!
!  Construct the edge list.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Construct the edge list.'
  write ( *, '(a)' ) '(While doing so, check for edges used more'
  write ( *, '(a)' ) 'than twice.)'

  call face_to_edge ( edge, face, face_order, ierror, max_edge, max_order, &
    num_edge, num_face )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FACE_CHECK(): Fatal error!'
    write ( *, '(a)' ) '  FACE_TO_EDGE failed.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Edge, Node1, Node2, Face1, Face2, Tier, Object'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' I, node1(i), node2(i), face1(i), face2(i)'
  write ( *, '(a)' ) ' '

  do i = 1, num_edge
    write ( *, '(10i3)' ) i, ( edge(j,i), j = 1, 4 )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Face, Order, Nodes'
  write ( *, '(a)' ) ' '
  do i = 1, num_face
    write ( *, '(10i3)' ) i, face_order(i), ( face(j,i), j = 1, face_order(i) )
  end do
!
!  Now force faces to have a consistent orientation.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Force faces to consistent orientation.'
  
  call face_flip ( edge, face, face_order, max_edge, max_order, num_edge, &
    num_face, num_fix )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Face, Order, Nodes'
  write ( *, '(a)' ) ' '
  do i = 1, num_face
    write ( *, '(10i3)' ) i, face_order(i), ( face(j,i), j = 1, face_order(i) )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'List boundary edges.'

  call edge_bound ( edge, max_edge, num_edge )

  return
end
subroutine face_example_box ( face, face_order, max_face, max_order, num_face )

!*****************************************************************************80
!
!! face_example_box() returns the faces of a simple box.
!
!  Diagram:
!
!    1---------------------------4
!    |\                         /|
!    | \                       / |
!    |  \         1           /  |
!    |   \                   /   |
!    |    2-----------------3    |
!    |    |                 |    |
!    |    |                 |    |
!    |  3 |       4         | 5  |
!    |    |                 |    |
!    |    |                 |    |
!    |    6-----------------7    |
!    |   /                   \   |
!    |  /                     \  |
!    | /          2            \ |
!    |/                         \|
!    5---------------------------8
!
!  Discussion:
!
!    This routine is used to supply some very simple data for the 
!    face checking routines.
!
!    This is "almost" a cube, except that one face is missing.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer FACE(MAX_ORDER,NUM_FACE), describes the faces.
!    FACE(I,J) is the index of the I-th node in face J.  It is best
!    if the nodes of all faces are listed in counterclockwise order.
!
!    Output, integer FACE_ORDER(NUM_FACE), is the number of nodes
!    making up each face.
!
!    Input, integer MAX_FACE, the maximum number of faces allowed.
!
!    Input, integer MAX_ORDER, is the maximum number of nodes that
!    can make up a face, required to dimension FACE.
!
!    Output, integer NUM_FACE, the number of faces.
!
  implicit none

  integer max_order
  integer max_face

  integer face(max_order,max_face)
  integer face_order(max_face)
  integer num_face

  num_face = 5

  if ( max_face < num_face ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FACE_EXAMPLE_OPEN_BOX(): Fatal error!'
    write ( *, '(a,i8)' ) '  Increase MAX_FACE to ', num_face
    stop 1
  end if

  face(1,1) = 1
  face(2,1) = 2
  face(3,1) = 3
  face(4,1) = 4

  face(1,2) = 5
  face(2,2) = 6
  face(3,2) = 7
  face(4,2) = 8

  face(1,3) = 1
  face(2,3) = 2
  face(3,3) = 6
  face(4,3) = 5

  face(1,4) = 6
  face(2,4) = 7
  face(3,4) = 3
  face(4,4) = 2

  face(1,5) = 3
  face(2,5) = 4
  face(3,5) = 8
  face(4,5) = 7

  face_order(1:num_face) = 4

  return
end
subroutine face_example_pieces ( face, face_order, max_face, max_order, &
  num_face )

!*****************************************************************************80
!
!! face_example_pieces() returns the faces of a set of three objects.
!
!  Diagram:
!
!    1---------------------------4
!    |\                         /|
!    | \                       / |       9--------10
!    |  \        7            /  |       |         |
!    |   \                   /   |       |   1     |
!    |    2-----------------3    |       |         |
!    |    |                 |    |       |         |
!    |    |                 |    |       11-------12
!    |  3 |       4         | 5  |        \       /
!    |    |                 |    |         \  6  /
!    |    |                 |    |          \   /
!    |    6-----------------7    |           \ /
!    |   /                   \   |           13
!    |  /                     \  |           / \
!    | /          8            \ |          /   \
!    |/                         \|         /  2  \
!    5---------------------------8        /       \
!                                        14-------15
!
!  Discussion:
!
!    THREE_PIECE is used to supply some very simple data for the 
!    face checking routines.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer FACE(MAX_ORDER,MAX_FACE), describes the faces.
!    FACE(I,J) is the index of the I-th node in face J.  It is best
!    if the nodes of all faces are listed in counterclockwise order.
!
!    Output, integer FACE_ORDER(MAX_FACE), is the number of nodes
!    making up each face.
!
!    Input, integer MAX_FACE, the maximum number of faces allowed.
!
!    Input, integer MAX_ORDER, is the maximum number of nodes that
!    can make up a face, required to dimension FACE.
!
!    Output, integer NUM_FACE, the number of faces.
!
  implicit none

  integer max_order
  integer max_face

  integer face(max_order,max_face)
  integer face_order(max_face)
  integer num_face

  num_face = 8

  if ( max_face < num_face ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FACE_EXAMPLE_PIECES(): Fatal error!'
    write ( *, '(a)' ) '  MAX_FACE < NUM_FACE!'
    write ( *, '(a,i8)' ) '  NUM_FACE = ', num_face
    write ( *, '(a,i8)' ) '  MAX_FACE = ', max_face
    stop 1
  end if

  face(1,1) = 9
  face(2,1) = 10
  face(3,1) = 12
  face(4,1) = 11

  face(1,2) = 14
  face(2,2) = 13
  face(3,2) = 15

  face(1,3) = 1
  face(2,3) = 2
  face(3,3) = 6
  face(4,3) = 5

  face(1,4) = 6
  face(2,4) = 7
  face(3,4) = 3
  face(4,4) = 2

  face(1,5) = 3
  face(2,5) = 4
  face(3,5) = 8
  face(4,5) = 7

  face(1,6) = 13
  face(2,6) = 12
  face(3,6) = 11

  face(1,7) = 1
  face(2,7) = 2
  face(3,7) = 3
  face(4,7) = 4

  face(1,8) = 5
  face(2,8) = 6
  face(3,8) = 7
  face(4,8) = 8

  face_order(1) = 4
  face_order(2) = 3
  face_order(3) = 4
  face_order(4) = 4
  face_order(5) = 4
  face_order(6) = 3
  face_order(7) = 4
  face_order(8) = 4

  return
end
subroutine face_flip ( edge, face, face_order, max_edge, max_order, &
  num_edge, num_face, num_fix )

!*****************************************************************************80
!
!! face_flip() flips faces to achieve a consistent orientation.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer EDGE(4,MAX_EDGE), edge information.
!    EDGE(1,I) is the starting node of edge I;
!    EDGE(2,I) is the ending node of edge I;
!    EDGE(3,I) is the positive face;
!    EDGE(4,I) is the negative face, if any.
!
!    Input, integer FACE(MAX_ORDER,NUM_FACE), describes the faces.
!    FACE(I,J) is the index of the I-th node in face J.  It is best
!    if the nodes of all faces are listed in counterclockwise order.
!
!    Input, integer FACE_ORDER(NUM_FACE), is the number of nodes
!    making up each face.
!
!    Input, integer MAX_EDGE, the maximum number of edges.
!
!    Input, integer MAX_ORDER, the maximum number of nodes that can
!    make up a face, required to dimension FACE.
!
!    Input, integer NUM_EDGE, the number of edges.
!
!    Input, integer NUM_FACE, the number of faces.
!
!    Output, integer NUM_FIX, the number of bad faces that were
!    found.
!
  implicit none

  integer max_edge
  integer max_order
  integer num_face

  integer edge(4,max_edge)
  integer f1
  integer f2
  integer face(max_order,num_face)
  integer face_order(num_face)
  integer iedge
  integer j
  integer jp1
  integer m1
  integer m2
  integer n1
  integer n2
  integer num_edge
  integer num_fix

  num_fix = 0

  do iedge = 1, num_edge

    n1 = edge(1,iedge)
    n2 = edge(2,iedge)
    f1 = edge(3,iedge)
    f2 = edge(4,iedge)
!
!  For now, just whine unless (N1,N2) is positive in F1 and negative in F2.
!
    if ( f1 /= 0 ) then

      do j = 1, face_order(f1)

        if ( j < face_order(f1) ) then
          jp1 = j + 1
        else
          jp1 = j
        end if

        m1 = face(j,f1)
        m2 = face(jp1,f1)

        if ( m1 == n1 .and. m2 == n2 ) then
          exit
        end if

        if ( m1 == n2 .and. m2 == n1 ) then
          num_fix = num_fix + 1
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'FACE_FLIP - Warning!'
          write ( *, '(a)' ) 'Bad orientation'
          write ( *, '(a,i8)' ) '  Face = ', f1
          write ( *, '(a,i8)' ) '  Side = ', j
          exit
        end if

      end do

    end if

    if ( f2 /= 0 ) then

      do j = 1, face_order(f2)

        if ( j < face_order(f2) ) then
          jp1 = j + 1
        else
          jp1 = j
        end if

        m1 = face(j,f2)
        m2 = face(jp1,f2)

        if ( m1 == n2 .and. m2 == n1 ) then
          exit
        end if

        if ( m1 == n1 .and. m2 == n2 ) then
          num_fix = num_fix + 1
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'FACE_FLIP - Warning!'
          write ( *, '(a)' ) 'Bad orientation'
          write ( *, '(a,i8)' ) '  Face = ', f2
          write ( *, '(a,i8)' ) '  Side = ', j
          exit
        end if

      end do

    end if

  end do

  if ( 0 < num_fix ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FACE_FLIP - Warning:'
    write ( *, '(a,i8)' ) '  Number of badly oriented faces = ', num_fix
  end if

  return
end
subroutine face_to_iv ( file_name, face, face_order, inode, jnode, nedge, &
  maxnode, maxface, maxorder, nnode, nface, x, y, z )

!*****************************************************************************80
!
!! face_to_iv() writes some simple graphics data to an Inventor file.
!
!  Example:
!
!     #Inventor V2.0 ascii
!
!     Separator {
!       Separator {
!         LightModel {
!           model PHONG
!         }
!         Material {
!           ambientColor  0.2 0.2 0.2
!           diffuseColor  0.8 0.8 0.8
!           emissiveColor 0.0 0.0 0.0
!           specularColor 0.0 0.0 0.0
!           shininess     0.2
!           transparency  0.0
!         }
!         Coordinate3 {
!           point [
!                8.59816       5.55317      -3.05561,
!                8.59816       2.49756      0.000000D+00,
!                ...etc...
!                2.48695       2.49756      -3.05561,
!           ]
!         }
!         IndexedFaceSet {
!           coordIndex [
!              0,    1,    2,   -1,    3,    4,    5,   -1,    7,    8,    9,
!            ...etc...
!            191,   -1,
!           ]
!         }
!       }
!     }
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the file name.
!
!    Input, integer FACE(MAX_ORDER,MAX_FACE), the nodes making
!    faces.
!
!    Input, integer FACE_ORDER(MAX_FACE), the number of nodes per
!    face.
!
!    Input, integer INODE(NEDGE), JNODE(NEDGE), node pairs for
!    edges.
!
!    Input, integer NEDGE, the number of edges.
!
!    Input, integer MAXNODE, the maximum number of nodes.
!
!    Input, integer MAXFACE, the maximum number of faces.
!
!    Input, integer MAXORDER, the maximum number of nodes per face.
!
!    Input, integer NNODE, the number of points.
!
!    Input, integer NFACE, the number of faces.
!
!    Input, real ( kind = rk ) X(MAXNODE), Y(MAXNODE), Z(MAXNODE), 
!    the coordinates of points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: OFFSET = 1

  integer maxnode
  integer maxface
  integer maxorder
  integer nedge

  integer face(maxorder,maxface)
  integer face_order(maxface)
  character ( len = * ) file_name
  integer icor3
  integer iface
  integer inode(nedge)
  integer ios
  integer itemp
  integer iunit
  integer ivert
  integer j
  integer jnode(nedge)
  integer length
  integer nnode
  integer nface
  character ( len = 200 ) text
  character ( len = 20 ) word
  real ( kind = rk ) x(maxnode)
  real ( kind = rk ) y(maxnode)
  real ( kind = rk ) z(maxnode)

  call get_unit ( iunit )

  open ( unit = iunit, file = file_name, status = 'replace', iostat = ios )

  if ( ios /= 0 ) then
    return
  end if

  write ( iunit, '(a)' ) '#Inventor V2.0 ascii'
  write ( iunit, '(a)' ) ' '
  write ( iunit, '(a)' ) 'Separator {'
  write ( iunit, '(a)' ) '  Separator {'
!
!  LightModel:
!
!    BASE_COLOR ignores light sources, and uses only diffuse color
!      and transparency.  Even without normal vector information,
!      the object will show up.  However, you won't get shadow
!      and lighting effects.
!
!    PHONG uses the Phong lighting model, accounting for light sources
!      and surface orientation.  This is the default.  I believe
!      you need accurate normal vector information in order for this
!      option to produce nice pictures.
!
!    DEPTH ignores light sources, and calculates lighting based on
!      the location of the object within the near and far planes
!      of the current camera's view volume.
!
  write ( iunit, '(a)' ) '    LightModel {'
  write ( iunit, '(a)' ) '      model PHONG'
  write ( iunit, '(a)' ) '    }'
!
!  Material
!
  write ( iunit, '(a)' ) '    Material {'
  write ( iunit, '(a)' ) '      ambientColor  0.5 0.2 0.2'
  write ( iunit, '(a)' ) '      diffuseColor  0.5 0.2 0.3'
  write ( iunit, '(a)' ) '      emissiveColor 0.5 0.0 0.0'
  write ( iunit, '(a)' ) '      specularColor 0.5 0.0 0.0'
  write ( iunit, '(a)' ) '      shininess     0.5'
  write ( iunit, '(a)' ) '      transparency  0.0'
  write ( iunit, '(a)' ) '    }'
!
!  Point coordinates.
!
  write ( iunit, '(a)' ) '    Coordinate3 {'
  write ( iunit, '(a)' ) '      point ['

  do icor3 = 1, nnode
    write ( text, '(3f12.4,'','')' ) x(icor3), y(icor3), z(icor3)
    call s_blanks_delete ( text )
    write ( iunit, '(8x,a)' ) trim ( text )
  end do

  write ( iunit, '(a)' ) '      ]'
  write ( iunit, '(a)' ) '    }'
  write ( iunit, '(a)' ) '    IndexedLineSet {'
!
!  IndexedLineSet coordIndex
!
    write ( iunit, '(a)' ) '      coordIndex ['

    do j = 1, nedge
      write ( iunit, '(8x,i8,'','',i8,'','',i8,'','')' ) &
        inode(j) - OFFSET, jnode(j)-offset, -1
    end do

    write ( iunit, '(a)' ) '      ]'

    write ( iunit, '(a)' ) '    }'
!
!  IndexedFaceSet.
!
  if ( 0 < nface ) then

    write ( iunit, '(a)' ) '    IndexedFaceSet {'
!
!  IndexedFaceSet coordIndex
!
    write ( iunit, '(a)' ) '      coordIndex ['

    text = ' '
    length = 0

    do iface = 1, nface

      do ivert = 1, face_order(iface) + 1

        if ( ivert <= face_order(iface) ) then
          itemp = face(ivert,iface) - OFFSET
        else
          itemp = 0 - OFFSET
        end if

        write ( word, '(i8,'','')' ) itemp

        call s_cat ( text, word, text )
        length = length + 1

        if ( itemp == -1 .or. 10 <= length .or. &
          ( iface == nface .and. ivert == face_order(iface) + 1 ) ) then

          call s_blanks_delete ( text )
          write ( iunit, '(8x,a)' ) trim ( text )
          text = ' '
          length = 0

        end if

      end do

    end do

    write ( iunit, '(a)' ) '      ]'

    write ( iunit, '(a)' ) '    }'

  end if
!
!  Close up the Separator node.
!
  write ( iunit, '(a)' ) '  }'
!
!  Close up the Separator node.
!
  write ( iunit, '(a)' ) '}'

  close ( unit = iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FACE_TO_IV:'
  write ( *, '(a)' ) '  The data was written to the file: ' &
    // trim ( file_name )

  return
end
subroutine face_sort ( face, face_object, face_order, face_tier, max_order, &
  num_face )

!*****************************************************************************80
!
!! face_sort() renumbers the faces in order of object and tier.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer FACE(MAX_ORDER,NUM_FACE), describes the
!    faces.  FACE(I,J) is the index of the I-th node in face J.  It is best
!    if the nodes of all faces are listed in counterclockwise order.
!
!    Input/output, integer FACE_OBJECT(NUM_FACE), describes the
!    objects.  FACE_OBJECT(I) is the index of the edge-connected "object" to 
!    which face I belongs.
!
!    Input/output, integer FACE_ORDER(NUM_FACE), is the number of
!    nodes making up each face.
!
!    Input/output, integer FACE_TIER(NUM_FACE).  FACE_TIER(I) is
!    the "tier" of face I in its object.  The base of the object has tier 1,
!    the neighbors of the base have tier 2, and so on.
!
!    Input, integer MAX_ORDER, is the maximum number of nodes that
!    can make up a face, required to dimension FACE.
!
!    Input, integer NUM_FACE, the number of faces.
!
  implicit none

  integer max_order
  integer num_face

  integer face(max_order,num_face)
  integer face_object(num_face)
  integer face_order(num_face)
  integer face_tier(num_face)
  integer i
  integer iface
  integer indx
  integer isgn
  integer jface

  iface = 0
  jface = 0
  indx = 0
  isgn = 0

  do

    call sort_heap_external ( num_face, indx, iface, jface, isgn )
!
!  Interchange faces IFACE and JFACE.
!
    if ( 0 < indx ) then

      do i = 1, max_order
        call i4_swap ( face(i,iface), face(i,jface) )
      end do

      call i4_swap ( face_object(iface), face_object(jface) )
      call i4_swap ( face_order(iface), face_order(jface) )
      call i4_swap ( face_tier(iface), face_tier(jface) )
!
!  Compare faces IFACE and JFACE.
!
    else if ( indx < 0 ) then

      if ( ( face_object(iface) < face_object(jface) ) .or. &
           ( face_object(iface) == face_object(jface) .and. &
             face_tier(iface) < face_tier(jface) ) ) then
        isgn = -1
      else
        isgn = +1
      end if

    else

      exit

    end if

  end do

  return
end
subroutine face_to_edge ( edge, face, face_order, ierror, max_edge, &
  max_order, num_edge, num_face )

!*****************************************************************************80
!
!! face_to_edge() converts face data to edge data.
!
!  Discussion:
!
!    The computation will fail if:
!
!    * More than two faces claim to share an edge (Node1,Node2).
!    * Not enough storage is set aside by MAX_EDGE.
!
!    If is expected that the edge (Node1,Node2) in Face1 is traversed in
!    the opposite sense, as (Node2,Node1), in Face2.  If this is not the
!    case, then some faces may need to be reoriented, but that will not
!    affect the computation.
!    
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer EDGE(4,MAX_EDGE), edge information.
!    EDGE(1,I) is the starting node of edge I;
!    EDGE(2,I) is the ending node of edge I;
!    EDGE(3,I) is the positive face;
!    EDGE(4,I) is the negative face, or 0 if the edge is used once.
!
!    Input, integer FACE(MAX_ORDER,NUM_FACE), describes the faces.
!    FACE(I,J) is the index of the I-th node in face J.  It is best
!    if the nodes of all faces are listed in counterclockwise order.
!
!    Input, integer FACE_ORDER(NUM_FACE), is the number of nodes
!    making up each face.
!
!    Output, integer IERROR, error flag: 0 = no error, 
!    nonzero = error.
!
!    Input, integer MAX_EDGE, the maximum number of edges.
!
!    Input, integer MAX_ORDER, the maximum number of nodes that can
!    make up a face, required to dimension FACE.
!
!    Output, integer NUM_EDGE, the number of edges.
!
!    Input, integer NUM_FACE, the number of faces.
!
  implicit none

  integer max_edge
  integer max_order
  integer num_face

  integer edge(4,max_edge)
  integer face(max_order,num_face)
  integer face_order(num_face)
  integer iedge
  integer ierror
  integer iface
  integer index
  integer j
  integer jp1
  integer n1
  integer n2
  integer num_edge
!
!  Initialize.
!
  ierror = 0

  edge(1:4,1:max_edge) = 0

  num_edge = 0
!
!  Consider face #I.
!
  do iface = 1, num_face
!
!  Seek an edge of face IFACE that already occurs in the edge list.
!  If there is one, then slide and reverse the entries in FACE(*,IFACE)
!  so that that edge occurs first, and in the opposite sense to its
!  occurrence in the edge list.
!
    call edge_match_face ( edge, max_edge, num_edge, face(1,iface), &
      face_order(iface), index )
!
!  Now, in any case, we know that the first two nodes in FACE(*,IFACE)
!  are the negative of an existing edge, or no nodes in FACE(*,IFACE)
!  occur in any existing edge.
!
    do j = 1, face_order(iface)

      n1 = face(j,iface)

      if ( j == face_order(iface) ) then
        jp1 = 1
      else
        jp1 = j + 1
      end if

      n2 = face(jp1,iface)

      call edge_match_nodes ( edge, max_edge, num_edge, n1, n2, iedge )

      if ( iedge == 0 ) then

        call edge_add_nodes ( edge, max_edge, num_edge, iface, n1, n2, ierror )

        if ( ierror /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'FACE_TO_EDGE(): Fatal error!'
          write ( *, '(a)' ) '  EDGE_ADD_NODES failed.'
          ierror = 1
          return
        end if

      else if ( edge(4,iedge) == 0 ) then

        edge(4,iedge) = iface

      else 

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FACE_TO_EDGE(): Fatal error!'
        write ( *, '(a,2i8)' ) '  Edge between nodes ', &
          edge(1,iedge), edge(2,iedge)
        write ( *, '(a)' ) '  is used at least 3 times, by faces:'
        write ( *, '(3i8)' ) edge(3,iedge), edge(4,iedge), iface
        ierror = 1
        return

      end if

    end do
  end do

  return
end
subroutine face_touch ( face, face_order, max_order, num_face, iface, jface, &
  touch )

!*****************************************************************************80
!
!! face_touch() reports whether two polygonal faces touch.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer FACE(MAX_ORDER,NUM_FACE), describes the faces.
!    FACE(I,J) is the index of the I-th node in face J.  It is best
!    if the nodes of all faces are listed in counterclockwise order.
!
!    Input, integer FACE_ORDER(NUM_FACE), is the number of nodes
!    making up each face.
!
!    Input, integer MAX_ORDER, is the maximum number of nodes that
!    can make up a face, required to dimension FACE.
!
!    Input, integer NUM_FACE, the number of faces.
!
!    Input, integer IFACE, JFACE, the faces to be checked.
!
!    Output, integer TOUCH:
!     0, the faces do not touch;
!    +1, the faces touch, both using an arc in the same direction;
!    -1, the faces touch, using an arc in opposite directions.
!
  implicit none

  integer max_order
  integer num_face

  integer face(max_order,num_face)
  integer face_order(num_face)
  integer i
  integer iface
  integer j
  integer jface
  integer m
  integer mp1
  integer mm1
  integer n
  integer np1
  integer touch

  touch = 0
!
!  Arc N1-N2 on IFACE must be matched by arc N1-N2 or N2-N1 on JFACE.
!
  do i = 1, face_order(iface)

    n = face(i,iface)
    if ( i < face_order(iface) ) then
      np1 = face(i+1,iface)
    else
      np1 = face(1,iface)
    end if

    do j = 1, face_order(jface)

      m = face(j,jface)
      if ( j < face_order(jface) ) then
        mp1 = face(j+1,jface)
      else
        mp1 = face(1,jface)
      end if

      if ( 1 < j ) then
        mm1 = face(j-1,jface)
      else
        mm1 = face(face_order(jface),jface)
      end if

      if ( n == m ) then
        if ( np1 == mp1 ) then
          touch = + 1
          return
        else if ( np1 == mm1 ) then
          touch = - 1
          return
        end if
      end if

    end do
  end do

  return
end
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
subroutine graph_adj_print ( adj, nnode, title )

!*****************************************************************************80
!
!! graph_adj_print() prints out an adjacency matrix for a graph.
!
!  Discussion:
!
!    This routine actually allows the entries of ADJ to have ANY value.
!    Values between 0 and 9 will be printed as is.  Other values will
!    be printed as '*'.
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
!    Input, integer ADJ(NNODEA,NNODE), the adjacency matrix of a 
!    graph.  ADJ(I,J) is 1 if there is a direct connection FROM node I TO 
!    node J, and is 0 otherwise.
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
  integer jhi
  character ( len = 80 ) string
  character ( len = * ) title

  if ( len_trim ( title ) /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  do i = 1, nnode

    jhi = min ( nnode, 74 )

    do j = 1, jhi

      if ( 0 <= adj(i,j) .and. adj(i,j) <= 9 ) then
        string(j:j) = char ( 48 + adj(i,j) )
      else
        string(j:j) = '*'
      end if

    end do

    write ( *, '(2x,i3,2x,a)' ) i, string(1:jhi)

  end do

  return
end
subroutine graph_arc_edge_sort ( nedge, inode, jnode )

!*****************************************************************************80
!
!! graph_arc_edge_sort() sorts the edge array of a graph.
!
!  Discussion:
!
!    The pair of nodes (I,J) representing an edge is reordered so
!    that the smaller node is listed first.
!
!    Then the edges are sorted in dictionary order.
!
!  Example:
!
!    Input:
!
!      INODE  JNODE
!
!        3      2
!        4      3
!        2      1
!        1      4
!
!    Output:
!
!      INODE  JNODE
!
!        1      2
!        1      4
!        2      3
!        3      4
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    24 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NEDGE, the number of edges.
!
!    Input/output, integer INODE(NEDGE), JNODE(NEDGE), the edge 
!    array of a graph.  The I-th edge of the graph connects nodes INODE(I) and
!    JNODE(I).
!
!    On output, the INODE and JNODE arrays have been sorted as described.
!
  implicit none

  integer nedge

  integer i
  integer iedge
  integer indx
  integer inode(nedge)
  integer isgn
  integer jedge
  integer jnode(nedge)

  if ( nedge <= 1 ) then
    return
  end if
!
!  Sort the node pairs.
!
  do i = 1, nedge
    if ( jnode(i) < inode(i) ) then
      call i4_swap ( inode(i), jnode(i) )
    end if
  end do
!
!  Sort the edges using an external heap sort.
!
  iedge = 0
  jedge = 0
  indx = 0
  isgn = 0

  do

    call sort_heap_external ( nedge, indx, iedge, jedge, isgn )
!
!  Interchange edges IEDGE and JEDGE.
!
    if ( 0 < indx ) then

      call i4_swap ( inode(iedge), inode(jedge) )
      call i4_swap ( jnode(iedge), jnode(jedge) )
!
!  Compare edges IEDGE and JEDGE.
!
    else if ( indx < 0 ) then

      if ( ( inode(iedge) < inode(jedge) ) .or. &
         ( inode(iedge) == inode(jedge) .and. &
           jnode(iedge) < jnode(jedge) ) ) then
        isgn = -1
      else
        isgn = +1
      end if

    else

      exit

    end if

  end do
 
  return
end
subroutine graph_arc_face ( face, face_count, face_order, iface, jface, &
  inode, jnode, maxface, maxorder, nedge, nface, nnode )

!*****************************************************************************80
!
!! graph_arc_face() constructs a set of faces for a plane graph.
!
!  Discussion:
!
!    Warning: This is an experimental code.
!
!    The reason this routine was written was to handle the problem of
!    converting certain forms of 3D graphics data from a point and line
!    representation to a list of faces.  While at first glance, this
!    seemed an easy task, it turned out to be one of those problems
!    that becomes harder the longer it is considered.  Particularly
!    vexing was the idea that it might be possible to do this reconstruction
!    without using any of the geometric data supplied with the connectivity
!    data.
!
!    The guiding idea was that a face ought to be a "short" cycle of
!    the graph, and that every edge ought to appear in two separate faces.
!    The resulting method should work for a connected graph which is
!    planar, or merely orientable.  A planar graph will result from a
!    reasonable "triangulation" (meaning decomposition into arbitrary
!    polygons) of the surface of a 3D object that has no holes.
!
!    This algorithm will also handle the case where the graph is not planar,
!    but results from the triangulation of a more complicated 3D object,
!    such as one that includes holes.  Even a Klein bottle, which is
!    a manifold, but not orientable, can be handled, although it may not
!    be possible then to assign a consistent orientation to the faces.
!
!    By the way, this problem is MUCH easier if we can assume that all
!    the faces use the same number of edges, such as in a triangular
!    decomposition.  This algorithm makes no such assumption.
!
!    If the graph is planar, then the decomposition into
!    faces allows us to define the dual graph.  The dual graph H of the
!    planar graph G comprises:
!
!    * nodes V(I), each of which corresponds to a face F(I) of G;
!
!    * edges (V(I), V(J)).  V(I) and V(J) share an edge in H if and only
!      if the faces F(I) and F(J) share an edge in G.  (Thus G and H
!      have the same number of edges).
!
!    In the terminology of this routine, the dual graph has NFACE nodes,
!    NEDGE edges, and the corresponding edge arrays are simply IFACE and
!    JFACE.
!
!  Formula:
!
!    If the graph is actually planar, we can regard it as the flattened
!    triangulation of a connected solid shape, and so we can apply Euler's
!    formula:
!
!      Faces + Vertices = Edges + 2
!
!    This means that we can predict beforehand that the number of faces
!    produced by this routine will be
!
!      NFACE = NEDGE + 2 - NNODE.
!
!  Notes:
!
!    The faces produced by this routine may actually overlap.  Without
!    geometric data, this is surely a possibility, since a graph may
!    have more than one embedding.  For instance, consider the following
!    two embeddings of the same graph:
!
!      A-----B       A-----B
!      |     |       |     |
!      |  E  |       D-----C
!      | / \ |        \   /
!      |/   \|         \ /
!      D-----C          E
!
!    This routine will report the two faces (A,B,C,D) and (C,D,E),
!    although in the first embedding one face seems to be part of
!    another.  This is not as bad as it might seem, since the routine
!    prefers the "smaller" face (A,B,C,D) over (A,B,C,E,D).
!
!
!    A second problem is best illustrated with a simple example.
!    Suppose we have a thin triangular rod, and that we have triangulated
!    the surface of this rod, so that the cross section of the rod
!    is a triangular graph, and the sides are made up of, say, squares.
!    Then this routine will report all the "internal" triangles as
!    faces.  It will still find the "true" faces on the sides, but
!    since it is possible to go around the diameter of the object
!    in a very few steps, the algorithm produces faces we would not
!    expect.
!
!  Restrictions:
!
!    The algorithm will fail if the graph cannot be regarded, at least
!    locally, as the triangulation of a smooth surface.  Smoothness
!    problems will not occur if the graph is planar, or results from
!    the triangulation of a 3D object, which might include holes.
!
!    The graph should be connected.
!
!    There should be no nodes of degree 1.
!
!  Method:
!
!    We have no geometric data from which to deduce physical positions
!    of the nodes.  We are only given that the graph is planar, so that
!    there is at least one embedding of the graph in the plane.
!
!    Our data structure for the method will use arrays called IFACE and JFACE.
!    For each edge I, IFACE(I) and JFACE(I) will eventually hold the
!    indices of the two faces that the edge is part of.  We begin
!    the algorithm by setting all entries of IFACE and JFACE to 0.
!
!    The second step is to find one cycle in the graph, of the shortest
!    length possible.  This cycle constitutes our first face.  We update
!    the appropriate entries of IFACE or JFACE, marking each edge as having
!    been used once.
!
!    The third step is to add one more face to our collection of faces.
!    The new face will be adjacent to the current collection of faces,
!    but will include at least one completely unused edge, if possible.
!
!    To guarantee this, we consider every edge that is part of our
!    collection of faces, and that has only been used once.  We look
!    at the endpoints of each of these edges., and
!
!      We search for an adjacent edge that has not been used.
!      If we find such an edge, then the first two edges of our next face
!      are the edge that was already part of the set of faces, and the
!      unused edge.
!
!      If we cannot find such an edge, then we repeat the search, starting
!      with an edge in the face set that has been used once.  But now
!      when we search for adjacent edges, we will consider using one that
!      has already been used once.
!
!    We then search for a path that will return us to the initial
!    node of the first edge.  Using a breadth-first search, we expect
!    to find the shortest path back to that node, and we assume that
!    this represents a face.  Again, we update the IFACE and JFACE arrays.
!
!    We repeat the third step until there are no more edges in the
!    collection of faces that have not been used twice.  Assuming the
!    graph is connected, this means that every face has been found.
!
!  Improvements:
!
!    It shouldn't be hard to modify the code to handle graphs that are
!    not connected.
!
!    If the edge arrays INODE and JNODE were sorted and indexed, some
!    operations could be done more efficiently.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer FACE(MAXORDER,MAXFACE), contains the list of 
!    edges which make up each face.  Face I is made up of the edges
!    FACE(1,I) through FACE(FACE_ORDER(I),I).
!
!    Output, integer FACE_COUNT(NEDGE).  For each edge I, 
!    FACE_COUNT(I) is the number of faces to which the edge belongs.  This 
!    value should be 0, 1 or 2.
!
!    Output, integer IFACE(NEDGE), JFACE(NEDGE).  IFACE(I) and 
!    JFACE(I) are the two faces to which edge I belongs.  Either or both may 
!    be zero if the algorithm fails.
!
!    Input, integer INODE(NEDGE), JNODE(NEDGE), the edge list for 
!    the graph.  The I-th edge connects nodes INODE(I) and JNODE(I).
!
!    Input, integer MAXFACE, the maximum number of faces for which
!    storage has been set aside in FACE and FACE_ORDER.
!
!    Input, integer MAXORDER, the maximum number of edges for 
!    which storage has been set aside in FACE.
!
!    Input, integer NEDGE, the number of edges.
!
!    Output, integer NFACE, the number of faces found by the
!    algorithm.
!
!    Input, integer NNODE, the number of nodes.
!
!    Output, integer FACE_ORDER(MAXFACE).  The number of edges used
!    in constructing face I is stored in FACE_ORDER(I).
!
  implicit none

  logical, parameter :: debug = .false.

  integer maxface
  integer maxorder
  integer nedge
  integer nnode

  integer face(maxorder,maxface)
  integer face_count(nedge)
  integer face_order(maxface)
  integer faceval
  integer i
  integer iedge
  integer iface(nedge)
  integer inode(nedge)
  integer j
  integer jface(nedge)
  integer jnode(nedge)
  integer k
  integer length
  integer nface
  integer nface_old
  integer nodes(3)
  integer nstart
!
!  Initialization.  No arc belongs to any face.
!
  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRAPH_ARC_FACE - Debug:'
    write ( *, '(a)' ) '  Initialization'
  end if

  nface = 0
  face_count(1:nedge) = 0
  iface(1:nedge) = 0
  jface(1:nedge) = 0
  face_order(1:maxface) = 0
  face(1:maxorder,1:maxface) = 0
!
!  We start here.  We may also jump back here if we have used up all the
!  connected parts of a graph, and need to jump to a new piece.
!
!  Find one new face of minimal length.
!
! 5 continue
!
  nface_old = nface

  do length = 3, nnode

    do iedge = 1, nedge

      nodes(1) = inode(iedge)
      nodes(2) = jnode(iedge)
      nstart = 2

      call graph_arc_face_next ( face, face_count, face_order, iface, jface, &  
        inode, jnode, maxface, maxorder, nedge, nface, nnode, nodes, nstart )

      if ( nface_old < nface ) then
        go to 10
      end if

    end do

  end do

  if ( nface == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRAPH_ARC_FACE - Note.'
    write ( *, '(a)' ) '  Could not find any starting face.'
  end if

  go to 60
!
!  Find an edge that is in one face, but not two.
!
10    continue

  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRAPH_ARC_FACE - Debug:'
    write ( *, '(a,i8)' ) '  Found starting face #:', nface
    write ( *, '(a,i8)' ) '  Order is ', face_order(nface)
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Vertices:'
    write ( *, '(a)' ) ' '
    do i = 1, face_order(nface)
      write ( *, '(6i8)' ) face(i,nface)
    end do
  end if

  iedge = 0
!
!  Look for an edge with FACE_COUNT of 1.
!
20    continue

  iedge = iedge + 1

  if ( face_count(iedge) == 1 ) then
    go to 30
  else if ( iedge < nedge ) then
    go to 20
  else
    if ( debug ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GRAPH_ARC_FACE - Debug:'
      write ( *, '(a)' ) '  No more nearby edges to try.'
    end if
!
!  Here, I'd like to be able to jump back and scrounge for other
!  islands of edges, but something's not right.
!
!       go to 5

    go to 60

  end if
!
!  The face will start with the two nodes of edge IEDGE.
!
30    continue

  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRAPH_ARC_FACE - Debug:'
    write ( *, '(a)' ) '  Found a starting edge:'
    write ( *, '(i8)' ) inode(iedge)
    write ( *, '(i8)' ) jnode(iedge)
  end if

  nodes(1) = inode(iedge)
  nodes(2) = jnode(iedge)
!
!  Look for an edge incident to JNODE(IEDGE).  This new edge should have
!  been used just FACEVAL times already.  (FACEVAL is preferably 0, but if
!  we can't find any at that rate, we'll try FACEVAL = 1).
!
  faceval = 0

40    continue

  do i = 1, nedge

    if ( face_count(i) == faceval ) then

      if ( inode(i) == nodes(2) .and. jnode(i) /= nodes(1) ) then
        nodes(3) = jnode(i)
        go to 50
      else if ( jnode(i) == nodes(2) .and. inode(i) /= nodes(1) ) then
        nodes(3) = inode(i)
        go to 50
      end if

    end if

  end do
!
!  If we "fell through" with FACEVAL = 0, then try the search again
!  with FACEVAL = 1.
!
  if ( faceval == 0 ) then

    faceval = 1
    go to 40
!
!  If we fell through with FACEVAL = 1, then we couldn't find any
!  way to use this edge.  Mark it as though it were used, and move on.
!
  else

    if ( debug ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GRAPH_ARC_FACE - Debug:'
      write ( *, '(a)' ) '  Failure.'
      write ( *, '(a,i8)' ) '  Cannot hook up to edge IEDGE:', iedge
      write ( *, '(2i8)' ) nodes(1), nodes(2)
    end if

    face_count(iedge) = 2
    go to 20

  end if
!
!  Now call FACENEXT to search for the shortest cycle that begins
!  NODES(1), NODES(2), NODES(3), and which involves only edges that
!  have been used once or less.
!
50    continue

  nface_old = nface
  nstart = 3

  call graph_arc_face_next ( face, face_count, face_order, iface, jface, &
    inode, jnode, maxface, maxorder, nedge, nface, nnode, nodes, nstart )

  if ( nface_old < nface ) then

    if ( debug ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GRAPH_ARC_FACE - Debug'
      write ( *, '(a,i8)' ) '  NFACE_OLD = ', nface_old
      write ( *, '(a,i8)' ) '  NFACE = ', nface
      write ( *, '(a,i8)' ) '  Order is ', face_order(nface)
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Vertices:'
      write ( *, '(a)' ) ' '
      do i = 1, face_order(nface)
        write ( *, '(6i8)' ) face(i,nface)
      end do
      write ( *, '(a)' ) '  Trying the big loop again.'
    end if

    go to 10

  end if
!
!  The algorithm has failed.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GRAPH_ARC_FACE - Error!'
  write ( *, '(a)' ) '  The algorithm has failed.'
  write ( *, '(a)' ) '  Only some of the faces were found.'
!
!  Cleanup
!
60    continue

  do i = 1, nface
    face_order(i) = min ( face_order(i), maxorder )
  end do

  do i = 1, nface
    do j = 1, face_order(i)
      k = face(j,i)
      if ( k < 0 ) then
        face(j,i) = jnode(-k)
      else
        face(j,i) = inode(k)
      end if
    end do
  end do

  return
end
subroutine graph_arc_face_next ( face, face_count, face_order, iface, jface, &
  inode, jnode, maxface, maxorder, nedge, nface, nnode, nodes, nstart )

!*****************************************************************************80
!
!! graph_arc_face_next() completes the next face, given a few starting nodes.
!
!  Discussion:
!
!    This is a utility routine, called by GRAPH_ARC_FACE, which
!    constructs all the faces of a graph.  GRAPH_ARC_FACE finds the first
!    two or three nodes of a face, and then calls this routine, which
!    attempts to complete the face by using a breadth-first search
!    from the final given node of the face.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer FACE(MAXORDER,MAXFACE), contains the 
!    list of edges which make up each face.  Face I is made up of the edges
!    FACE(1,I) through FACE(FACE_ORDER(I),I).  If a new face is found, this
!    array is updated.
!
!    Input/output, integer FACE_COUNT(NEDGE).  For each edge I, 
!    FACE_COUNT(I) is the number of faces to which the edge belongs.  This value
!    will be 0, 1 or 2.  If a new face is found, this data is updated.
!
!    Input/output, integer FACE_ORDER(MAXFACE).  The number of 
!    edges used in constructing face I is stored in FACE_ORDER(I).
!
!    Input/output, integer IFACE(NEDGE), JFACE(NEDGE).  IFACE(I) 
!    and JFACE(I) are the two faces to which edge I belongs.  Either or both 
!    may be zero if the algorithm fails.  If a new face is found, this data 
!    is updated.
!
!    Input, integer INODE(NEDGE), JNODE(NEDGE), the edge list for 
!    the graph.  The I-th edge connects nodes INODE(I) and JNODE(I).
!
!    Input, integer MAXFACE, the maximum number of faces for which 
!    storage has been set aside in FACE and FACE_ORDER.
!
!    Input, integer MAXORDER, the maximum number of edges for 
!    which storage has been set aside in FACE.
!
!    Input, integer NEDGE, the number of edges.
!
!    Input/output, integer NFACE.  NFACE is the number of faces 
!    found so far.  This value will be updated by this routine if a new face 
!    is found.
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, integer NODES(NSTART), the first few nodes in the 
!    partial cycle.
!
!    Input, integer NSTART, the number of nodes in the partial
!    cycle.
!
!  Workspace:
!
!    Workspace, integer DAD(NNODE), used during the breadth first search
!    of the graph, to point backwards from each node to its predecessor
!    in a path.
!
!    Workspace, integer INDEX(NNODE), used during the breadth first search
!    to label nodes that have been visited.
!
  implicit none

  integer maxface
  integer maxorder
  integer nedge
  integer nnode
  integer nstart

  integer dad(nnode)
  integer face(maxorder,maxface)
  integer face_count(nedge)
  integer face_order(maxface)
  integer i
  integer iedge2
  integer iface(nedge)
  integer index(nnode)
  integer inode(nedge)
  integer istart1
  integer istart2
  integer itemp
  integer jface(nedge)
  integer jnode(nedge)
  integer kedge
  integer kk
  integer nadd
  integer nface
  integer nodei
  integer nodej
  integer npath
  integer nodes(nstart)
!
!  Initialization.
!
  index(1:nnode) = 0
  dad(1:nnode) = 0

  istart1 = nodes(1)
  istart2 = nodes(2)

  do i = 1, nstart

    npath = i

    if ( i == 1 ) then
      dad(nodes(i)) = -1
    else
      dad(nodes(i)) = nodes(i-1)
    end if

    index(nodes(i)) = i

  end do
!
!  From the nodes with INDEX = NPATH, consider all neighbors.
!
10    continue

  npath = npath + 1
  nadd = 0

  do iedge2 = 1, nedge

    if ( index(inode(iedge2)) == npath-1 .and. index(jnode(iedge2)) == 0 ) then

      nodei = inode(iedge2)
      nodej = jnode(iedge2)

    else if ( index(jnode(iedge2)) == npath-1 .and. &
      index(inode(iedge2)) == 0 ) then

      nodei = jnode(iedge2)
      nodej = inode(iedge2)

    else if ( index(inode(iedge2)) == npath-1 .and. &
      jnode(iedge2) == istart1 ) then

      nodei = inode(iedge2)
      nodej = jnode(iedge2)

    else if ( index(jnode(iedge2)) == npath-1 .and. &
      inode(iedge2) == istart1 ) then

      nodei = jnode(iedge2)
      nodej = inode(iedge2)

    else

      nodei = 0
      nodej = 0

    end if

    if ( nodei /= 0 .and. nodej /= istart1 ) then

      nadd = nadd + 1
      index(nodej) = npath
      dad(nodej) = nodei
!
!  Success if the marked node is the starting point (except when
!  using the edge (ISTART2,ISTART1)).
!
    else if ( nodej == istart1 .and. nodei == istart2 ) then

    else if ( nodej == istart1 .and. nodei /= istart2 ) then

      nface = nface + 1

20    continue
!
!  Find the edge KK which joins NODEI and NODEJ.
!
      do kk = 1, nedge

        if ( ( inode(kk) == nodei .and. jnode(kk) == nodej ) .or. &
             ( jnode(kk) == nodei .and. inode(kk) == nodej ) ) then

          face_count(kk) = face_count(kk) + 1
          itemp = face_count(kk)

          if ( itemp == 1 ) then
            iface(kk) = nface
          else if ( itemp == 2 ) then
            jface(kk) = nface
          end if

          if ( inode(kk) == nodei ) then
            kedge = kk
          else
            kedge = -kk
          end if

          exit

        end if

      end do

      nodej = nodei
!
!  Add the edge to the face-to-edge list.
!
      if ( nface <= maxface ) then

        if ( face_order(nface) < maxorder ) then
          face_order(nface) = face_order(nface) + 1
        end if

        if ( face_order(nface) <= maxorder ) then
          face(face_order(nface),nface) = kedge
        end if

      end if

      if ( nodej /= istart1 ) then
        nodei = dad(nodej)
        go to 20
      end if

      return

    end if

  end do
!
!  If we were able to proceed another step, and we haven't exceeded
!  our limit, then go back and take another step.
!
  if ( 0 < nadd .and. npath <= nnode ) then
    go to 10
  end if

  return
end
subroutine graph_arc_print ( nedge, inode, jnode, title )

!*****************************************************************************80
!
!! graph_arc_print() prints out a graph from an edge list.
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
!    Input, integer NEDGE, the number of edges.
!
!    Input, integer INODE(NEDGE), JNODE(NEDGE), the beginning 
!    and end nodes of the edges.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer nedge

  integer i
  integer inode(nedge)
  integer jnode(nedge)
  character ( len = * ) title

  if ( len_trim ( title ) /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  do i = 1, nedge
    write ( *, '(i8,4x,2i8)' ) i, inode(i), jnode(i)
  end do

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
subroutine graph_arc_to_ps ( file_name, inode, jnode, nedge, nnode, x, y )

!*****************************************************************************80
!
!! graph_arc_to_ps() writes graph information to a PostScript file.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the output file.
!
!    Input, integer INODE(NEDGE), JNODE(NEDGE), the edge array.
!    The I-th edge connects nodes INODE(I) and JNODE(I).
!
!    Input, integer NEDGE, the number of edges.
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, real ( kind = rk ) X(NNODE), Y(NNODE), the X and Y components
!    of points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nedge
  integer nnode

  real ( kind = rk ) alpha
  real ( kind = rk ) blue
  character ( len = 8 ) date
  character ( len = * ) file_name
  real ( kind = rk ) green
  integer inode(nedge)
  integer ios
  integer iunit
  integer jnode(nedge)
  integer margin
  integer pagexmax
  integer pagexmin
  integer pageymax
  integer pageymin
  integer plotxmax
  integer plotxmin
  integer plotxmin2
  integer plotymax
  integer plotymin
  integer plotymin2
  real ( kind = rk ) red
  real ( kind = rk ) x(nnode)
  real ( kind = rk ) xmax
  real ( kind = rk ) xmin
  real ( kind = rk ) y(nnode)
  real ( kind = rk ) ymax
  real ( kind = rk ) ymin
!
!  Bounding box.
!
  xmin = minval ( x(1:nnode) )
  xmax = maxval ( x(1:nnode) )
  ymin = minval ( y(1:nnode) )
  ymax = maxval ( y(1:nnode) )

  if ( xmin == xmax ) then
    xmin = x(1) - 0.5D+00
    xmax = x(1) + 0.5D+00
  end if

  if ( ymin == ymax ) then
    ymin = y(1) - 0.5D+00
    ymax = y(1) + 0.5D+00
  end if
!
!  Compute the scale factor.
!
  pagexmax = 612
  pagexmin = 0
  pageymax = 792
  pageymin = 0

  margin = 36

  plotxmax = pagexmax - margin
  plotxmin = pagexmin + margin
  plotymax = pageymax - margin
  plotymin = pageymin + margin

  alpha = min ( ( plotxmax - plotxmin ) / ( xmax - xmin ), &
                ( plotymax - plotymin ) / ( ymax - ymin ) )
!
!  Adjust PLOTXMIN and PLOTYMIN to center the image.
!
  plotxmin2 = int ( 0.5D+00 * ( plotxmin + plotxmax - alpha * ( xmax - xmin ) ) )
  plotymin2 = int ( 0.5D+00 * ( plotymin + plotymax - alpha * ( ymax - ymin ) ) )

  call get_unit ( iunit )

  open ( unit = iunit, file = file_name, status = 'replace', iostat = ios )

  if ( ios /= 0 ) then
    return
  end if
!
!  Write the prolog.
!
  write ( iunit, '(a)' ) '%!PS-Adobe-3.0'
  write ( iunit, '(a)' ) '%%Document-Fonts: Times-Roman'
  write ( iunit, '(a,a)' ) '%%Title: ' , trim ( file_name )
  write ( iunit, '(a)' ) '%%Creator: GRAFPACK(graph_arc_to_ps)'
  call date_and_time ( date )
  write ( iunit, '(a)' ) '%%CreationDate: ' // trim ( date )
  write ( iunit, '(a,4i5)' ) '%%BoundingBox', plotxmin, plotymin, plotxmax, &
    plotymax
  write ( iunit, '(a)' ) '%%LanguageLevel: 2'
  write ( iunit, '(a)' ) '%%EndComments'
  write ( iunit, '(a)' ) '%%BeginProlog'
  write ( iunit, '(a)' ) '%%EndProlog'
!
!  Set the line color.
!
  red = 0.0D+00
  green = 0.0D+00
  blue = 0.0D+00

  write ( iunit, '(3f7.4,a)' ) red, green, blue, ' setrgbcolor'
!
!  Draw lines.
!
  call edges_to_ps ( plotxmin2, plotymin2, alpha, iunit, inode, jnode, &
    nedge, nnode, x, y, xmin, ymin )
!
!  Set the fill color.
!
  red = 0.1
  green = 0.1
  blue = 0.7

  write ( iunit, '(3f7.4,a)' ) red, green, blue, ' setrgbcolor'
!
!  Draw points.
!
  call nodes_to_ps ( plotxmin2, plotymin2, alpha, iunit, nnode, x, y, &
    xmin, ymin )

  write ( iunit, '(a)' ) 'showpage'
!
!  Write the epilog.
!
  write ( iunit, '(a)' ) 'grestore'
  write ( iunit, '(a)' ) '%%Trailer'
  write ( iunit, '(a,i2)' ) '%%Pages: 1'

  close ( iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GRAPH_ARC_TO_PS'
  write ( *, '(a)' ) '  The data was written to the file: ' &
    // trim ( file_name )

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
subroutine grf_read ( file_name, inode, jnode, maxedge, maxnode, nedge, nnode, &
  x, y )

!*****************************************************************************80
!
!! grf_read() reads a GRF file containing a 2D representation of a graph.
!
!  Example:
!
!    #  A graph where every node has 3 neighbors.
!    #
!    1      0.546  0.956  5      6      2    
!    2      0.144  0.650  7      3      1    
!    3      0.326  0.188  8      4      2    
!    4      0.796  0.188  9      5      3    
!    5      0.988  0.646  10     4      1    
!    6      0.552  0.814  11     12     1    
!    7      0.264  0.616  11     15     2    
!    8      0.404  0.296  15     14     3    
!    9      0.752  0.298  14     13     4    
!    10     0.846  0.624  13     12     5    
!    11     0.430  0.692  16     6      7    
!    12     0.682  0.692  17     10     6    
!    13     0.758  0.492  18     9      10   
!    14     0.566  0.358  19     8      9    
!    15     0.364  0.484  20     7      8    
!    16     0.504  0.602  11     20     17   
!    17     0.608  0.602  12     18     16   
!    18     0.634  0.510  13     19     17   
!    19     0.566  0.444  14     20     18   
!    20     0.480  0.510  15     16     19   
!
!  Discussion:
!
!    The original GRF format has been modified so that a line starting
!    with a # is considered a comment line.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the file name.
!
!    Output, integer INODE(MAXEDGE), JNODE(MAXEDGE), the edges.  
!    The I-th edge joins nodes INODE(I) and JNODE(I).
!
!    Input, integer IUNIT, the FORTRAN unit number associated with
!    the graph file, which should already have been opened by the user.
!
!    Input, integer MAXEDGE, the maximum number of edges.
!
!    Input, integer MAXNODE, the maximum number of nodes.
!
!    Output, integer NEDGE, the number of edges that were read.
!
!    Output, integer NNODE, the number of nodes that were read.
!
!    Output, real ( kind = rk ) X(MAXNODE), Y(MAXNODE), the coordinates of the
!    nodes.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: maxchr = 200

  integer maxedge
  integer maxnode

  character ( len = * ) file_name
  integer ierror
  integer inode(maxedge)
  integer ios
  integer istring
  integer iunit
  integer jnode(maxedge)
  integer lchar
  integer nbad
  integer nedge
  integer nnode
  integer nodei
  integer nodej
  integer ntext
  character ( len = maxchr ) string
  real ( kind = rk ) x(maxnode)
  real ( kind = rk ) xval
  real ( kind = rk ) y(maxnode)
  real ( kind = rk ) yval

  nbad = 0
  nedge = 0
  nnode = 0
  ntext = 0

  call get_unit ( iunit )

  open ( unit = iunit, file = file_name, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRF_READ(): Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file.'
    return
  end if
!
!  Read information about each node.
!
  do

    read ( iunit, '(a)', iostat = ios ) string

    if ( ios /= 0 ) then
      exit
    end if

    ntext = ntext + 1

    if ( len ( string ) <= 0 ) then
      cycle
    end if

    if ( string(1:1) == '#' ) then
      cycle
    end if

    istring = 1
!
!  Extract the node index, NODEI.
!
    call s_to_i4 ( string(istring:), nodei, ierror, lchar )

    if ( ierror /= 0 .or. lchar == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GRF_READ(): Fatal error!'
      write ( *, '(a)' ) '  Unreadable node index value.'
      nbad = nbad + 1
      cycle
    end if

    istring = istring + lchar

    if ( nodei < 1 .or. maxnode < nodei ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GRF_READ(): Fatal error!'
      write ( *, '(a,i8)' ) '  Illegal node index value, NODEI = ', nodei
      cycle
    end if

    if ( nodei == nnode + 1 ) then
      nnode = nnode + 1
    else if ( nnode < nodei ) then
      nnode = nodei
    end if
!
!  Extract the X, Y coordinates of the node.
!
    call s_to_r8 ( string(istring:), xval, ierror, lchar )

    if ( ierror /= 0 .or. lchar == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GRF_READ(): Fatal error!'
      write ( *, '(a)' ) '  Unreadable X coordinate for node.'
      nbad = nbad + 1
      cycle
    end if

    istring = istring + lchar

    call s_to_r8 ( string(istring:), yval, ierror, lchar )

    if ( ierror /= 0 .or. lchar == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GRF_READ(): Fatal error!'
      write ( *, '(a)' ) '  Unreadable Y coordinate for node.'
      nbad = nbad + 1
      cycle
    end if

    istring = istring + lchar

    x(nodei) = xval
    y(nodei) = yval
!
!  Read the indices of the nodes to which NODEI is connected.
!
    do

      call s_to_i4 ( string(istring:), nodej, ierror, lchar )

      if ( ierror /= 0 .and. ierror /= 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRF_READ(): Fatal error!'
        write ( *, '(a)' ) '  Unreadable node neighbor value.'
        nbad = nbad + 1
        cycle
      end if

      istring = istring + lchar

      if ( lchar <= 0 ) then
        exit
      end if

      if ( 1 <= nodej .and. nodej <= maxnode ) then

        if ( nedge < maxedge ) then
          nedge = nedge + 1
          inode(nedge) = nodei
          jnode(nedge) = nodej
        end if

      end if

      if ( maxchr < istring ) then
        exit
      end if

    end do

  end do

  close ( unit = iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GRF_READ - Input file statistics:'
  write ( *, '(a,i8)' ) '  Text lines:     ', ntext
  write ( *, '(a,i8)' ) '  Bad text lines: ', nbad
  write ( *, '(a,i8)' ) '  Nodes:          ', nnode
  write ( *, '(a,i8)' ) '  Edges:          ', nedge

  return
end
subroutine hqr ( nm, n, low, igh, h, wr, wi, ierr )

!*****************************************************************************80
!
!! hqr() computes all eigenvalues of a real upper Hessenberg matrix.
!
!  Discussion:
!
!    This subroutine finds the eigenvalues of a real
!    upper Hessenberg matrix by the QR method.
!
!  Reference:
!
!    Martin, Peters, James Wilkinson,
!    HQR,
!    Numerische Mathematik,
!    Volume 14, pages 219-231, 1970.
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Y Ikebe, V Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of H, which must
!    be at least N.
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer LOW, IGH, two integers determined by the
!    routine BALANC.  If BALANC is not used, set LOW=1, IGH=N.
!
!    Input/output, real ( kind = rk ) H(NM,N), the N by N upper Hessenberg
!    matrix.  Information about the transformations used in the reduction to
!    Hessenberg form by ELMHES or ORTHES, if performed, is stored
!    in the remaining triangle under the Hessenberg matrix.
!    On output, the information in H has been destroyed.
!
!    Output, real ( kind = rk ) WR(N), WI(N), the real and imaginary parts of the
!    eigenvalues.  The eigenvalues are unordered, except that complex
!    conjugate pairs of values appear consecutively, with the eigenvalue
!    having positive imaginary part listed first.  If an error exit
!    occurred, then the eigenvalues should be correct for indices
!    IERR+1 through N.
!
!    Output, integer IERR, error flag.
!    0, no error.
!    J, the limit of 30*N iterations was reached while searching for
!       the J-th eigenvalue.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nm

  integer en
  integer enm2
  real ( kind = rk ) h(nm,n)
  integer i
  integer ierr
  integer igh
  integer itn
  integer its
  integer j
  integer k
  integer l
  integer ll
  integer low
  integer m
  integer mm
  integer na
  real ( kind = rk ) norm
  logical notlas
  real ( kind = rk ) p
  real ( kind = rk ) q
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) t
  real ( kind = rk ) tst1
  real ( kind = rk ) tst2
  real ( kind = rk ) w
  real ( kind = rk ) wi(n)
  real ( kind = rk ) wr(n)
  real ( kind = rk ) x
  real ( kind = rk ) y
  real ( kind = rk ) zz

  ierr = 0
  norm = 0.0D+00
  k = 1
!
!  Store roots isolated by BALANC and compute matrix norm.
!
  do i = 1, n

    do j = k, n
      norm = norm + abs ( h(i,j) )
    end do

    k = i
    if (i < low .or. igh < i ) then
      wr(i) = h(i,i)
      wi(i) = 0.0D+00
    end if

  end do

  en = igh
  t = 0.0D+00
  itn = 30 * n
!
!  Search for next eigenvalues.
!
60 continue

  if ( en < low ) then
    return
  end if

  its = 0
  na = en - 1
  enm2 = na - 1
!
!  Look for a single small sub-diagonal element.
!
70 continue

  do ll = low, en
    l = en + low - ll
    if ( l == low ) then
      exit
    end if
    s = abs ( h(l-1,l-1) ) + abs ( h(l,l) )
    if ( s == 0.0D+00 ) then
      s = norm
    end if
    tst1 = s
    tst2 = tst1 + abs ( h(l,l-1))
    if ( tst2 == tst1 ) then
      exit
    end if
  end do
!
!  Form shift.
!
  x = h(en,en)
  if ( l == en ) go to 270
  y = h(na,na)
  w = h(en,na) * h(na,en)
  if ( l == na ) go to 280

  if ( itn == 0 ) then
    ierr = en
    return
  end if
!
!  Form an exceptional shift.
!
  if ( its == 10 .or. its == 20 ) then

    t = t + x

    do i = low, en
      h(i,i) = h(i,i) - x
    end do

    s = abs ( h(en,na) ) + abs ( h(na,enm2) )
    x = 0.75D+00 * s
    y = x
    w = -0.4375D+00 * s * s

  end if

  its = its + 1
  itn = itn - 1
!
!  Look for two consecutive small sub-diagonal elements.
!
  do mm = l, enm2

    m = enm2 + l - mm
    zz = h(m,m)
    r = x - zz
    s = y - zz
    p = ( r * s - w ) / h(m+1,m) + h(m,m+1)
    q = h(m+1,m+1) - zz - r - s
    r = h(m+2,m+1)
    s = abs ( p ) + abs ( q ) + abs ( r )
    p = p / s
    q = q / s
    r = r / s

    if ( m == l ) then
      exit
    end if

    tst1 = abs ( p ) * ( abs ( h(m-1,m-1) ) + abs ( zz ) + abs ( h(m+1,m+1) ) )
    tst2 = tst1 + abs ( h(m,m-1) ) * ( abs ( q ) + abs ( r ) )

    if ( tst2 == tst1 ) then
      exit
    end if

  end do

  do i = m+2, en
    h(i,i-2) = 0.0D+00
    if ( i /= m+2 ) then
      h(i,i-3) = 0.0D+00
    end if
  end do
!
!  Double QR step involving rows l to EN and columns M to EN.
!
  do k = m, na

     notlas = k /= na

     if ( k == m ) go to 170

     p = h(k,k-1)
     q = h(k+1,k-1)
     if ( notlas ) then
       r = h(k+2,k-1)
     else
       r = 0.0D+00
     end if
     x = abs ( p ) + abs ( q ) + abs ( r )
     if ( x == 0.0D+00 ) go to 260
     p = p / x
     q = q / x
     r = r / x

170  continue

     s = sign ( sqrt ( p**2 + q**2 + r**2 ), p )

     if ( k /= m ) then
       h(k,k-1) = - s * x
     else if ( l /= m ) then
       h(k,k-1) = - h(k,k-1)
     end if

     p = p + s
     x = p / s
     y = q / s
     zz = r / s
     q = q / p
     r = r / p
     if ( notlas ) go to 225
!
!  Row modification.
!
     do j = k, n
       p = h(k,j) + q * h(k+1,j)
       h(k,j) = h(k,j) - p * x
       h(k+1,j) = h(k+1,j) - p * y
     end do

     j = min ( en, k+3 )
!
!  Column modification.
!
     do i = 1, j
       p = x * h(i,k) + y * h(i,k+1)
       h(i,k) = h(i,k) - p
       h(i,k+1) = h(i,k+1) - p * q
     end do

     go to 255

225  continue
!
!  Row modification.
!
     do j = k, n
       p = h(k,j) + q * h(k+1,j) + r * h(k+2,j)
       h(k,j) = h(k,j) - p * x
       h(k+1,j) = h(k+1,j) - p * y
       h(k+2,j) = h(k+2,j) - p * zz
     end do

     j = min ( en, k+3 )
!
!  Column modification.
!
     do i = 1, j
       p = x * h(i,k) + y * h(i,k+1) + zz * h(i,k+2)
       h(i,k) = h(i,k) - p
       h(i,k+1) = h(i,k+1) - p * q
       h(i,k+2) = h(i,k+2) - p * r
     end do

255 continue

260 continue

  end do

  go to 70
!
!  One root found.
!
270 continue

  wr(en) = x + t
  wi(en) = 0.0D+00
  en = na
  go to 60
!
!  Two roots found.
!
280 continue

  p = ( y - x ) / 2.0D+00
  q = p * p + w
  zz = sqrt ( abs ( q ) )
  x = x + t
!
!  Real root, or complex pair.
!
  if ( 0.0D+00 <= q ) then

    zz = p + sign ( zz, p )
    wr(na) = x + zz
    if ( zz == 0.0D+00 ) then
      wr(en) = wr(na)
    else
      wr(en) = x - w / zz
    end if
    wi(na) = 0.0D+00
    wi(en) = 0.0D+00

  else

    wr(na) = x + p
    wr(en) = x + p
    wi(na) = zz
    wi(en) = -zz

  end if

  en = enm2
  go to 60
end
function i4_modp ( i, j )

!*****************************************************************************80
!
!! i4_modp() returns the nonnegative remainder of integer division.
!
!  Discussion:
!
!    If
!      NREM = I4_MODP ( I, J )
!      NMULT = ( I - NREM ) / J
!    then
!      I = J * NMULT + NREM
!    where NREM is always nonnegative.
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360.0) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, I4_MODP(A,360.0) is between 0 and 360, always.
!
!  Example:
!
!        I         J     MOD   I4_MODP   I4_MODP Factorization
!
!      107        50       7       7    107 =  2 *  50 + 7
!      107       -50       7       7    107 = -2 * -50 + 7
!     -107        50      -7      43   -107 = -3 *  50 + 43
!     -107       -50      -7      43   -107 =  3 * -50 + 43
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
!    Input, integer I, the number to be divided.
!
!    Input, integer J, the number that divides I.
!
!    Output, integer I4_MODP, the nonnegative remainder 
!    when I is divided by J.
!
  implicit none

  integer i
  integer j
  integer i4_modp

  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MODP(): Fatal error!'
    write ( *, '(a,i8)' ) '  I4_MODP ( I, J ) called with J = ', j
    stop 1
  end if

  i4_modp = mod ( i, j )

  if ( i4_modp < 0 ) then
    i4_modp = i4_modp + abs ( j )
  end if

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
subroutine i4vec_uniform_ab ( n, a, b, x )

!*****************************************************************************80
!
!! i4vec_uniform_ab() returns a scaled pseudorandom I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    The pseudorandom numbers should be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the dimension of the vector.
!
!    integer A, B, the limits of the interval.
!
!  Output:
!
!    integer X(N), a vector of numbers between A and B.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer a
  integer b
  integer i
  real r
  integer value
  integer x(n)

  do i = 1, n

    call random_number ( harvest = r )
!
!  Scale R to lie between A-0.5 and B+0.5.
!
    r = ( 1.0E+00 - r ) * ( real ( min ( a, b ) ) - 0.5E+00 ) &
      +             r   * ( real ( max ( a, b ) ) + 0.5E+00 )
!
!  Use rounding to convert R to an integer between A and B.
!
    value = nint ( r )

    value = max ( value, min ( a, b ) )
    value = min ( value, max ( a, b ) )

    x(i) = value

  end do

  return
end
subroutine i4col_compare ( m, n, a, i, j, isgn )

!*****************************************************************************80
!
!! i4col_compare() compares columns I and J of a integer array.
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, I = 2, J = 4
!
!      A = (
!        1  2  3  4
!        5  6  7  8
!        9 10 11 12 )
!
!    Output:
!
!      ISGN = -1
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer A(M,N), an array of N columns of vectors of
!    length M.
!
!    Input, integer I, J, the columns to be compared.
!    I and J must be between 1 and N.
!
!    Output, integer ISGN, the results of the comparison:
!    -1, column I < column J,
!     0, column I = column J,
!    +1, column I > column J.
!
  implicit none

  integer m
  integer n

  integer a(m,n)
  integer i
  integer isgn
  integer j
  integer k
!
!  Check.
!
  if ( i < 1 .or. n < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_COMPARE(): Fatal error!'
    write ( *, '(a)' ) '  Column index I is out of bounds.'
    stop 1
  end if

  if ( j < 1 .or. n < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_COMPARE(): Fatal error!'
    write ( *, '(a)' ) '  Column index J is out of bounds.'
    stop 1
  end if

  isgn = 0

  if ( i == j ) then
    return
  end if

  k = 1

  do while ( k <= m )

    if ( a(k,i) < a(k,j) ) then
      isgn = - 1
      return
    else if ( a(k,j) < a(k,i) ) then
      isgn = + 1
      return
    end if

    k = k + 1

  end do

  return
end
subroutine i4col_sort_a ( nrow, ncol, ia )

!*****************************************************************************80
!
!! i4col_sort_a() ascending sorts an I4COL.
!
!  Definition:
!
!    In lexicographic order, the statement "X < Y", applied to two real
!    vectors X and Y of length NROW, means that there is some index I, with
!    1 <= I <= NROW, with the property that
!
!      X(J) = Y(J) for J < I, and
!      X(I) < Y(I).
!
!    In other words, the first time they differ, X is smaller.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NROW, the number of rows of A, and the length
!    of a vector of data.
!
!    Input, integer NCOL, the number of columns of A.
!
!    Input/output, integer IA(NROW,NCOL).
!    On input, the array of NCOL columns of NROW-vectors.
!    On output, the columns of A have been sorted in lexicographic order.
!
  implicit none

  integer ncol
  integer nrow

  integer ia(nrow,ncol)
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

    call sort_heap_external ( ncol, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      call i4col_swap ( nrow, ncol, ia, i, j )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4col_compare ( nrow, ncol, ia, i, j, isgn )

     else

      exit

    end if

  end do

  return
end
subroutine i4col_swap ( nrow, ncol, ia, i, j )

!*****************************************************************************80
!
!! i4col_swap() swaps columns I and J of an I4COL.
!
!  Example:
!
!    Input:
!
!      NROW = 3, NCOL = 4, I = 2, J = 4
!
!      IA = (
!        1  2  3  4
!        5  6  7  8
!        9 10 11 12 )
!
!    Output:
!
!      IA = (
!        1  4  3  2
!        5  8  7  6
!        9 12 11 10 )
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NROW, NCOL, the number of rows and columns in
!    the table.
!
!    Input, integer IA(NROW,NCOL), a table of numbers, regarded as
!    NCOL columns of vectors of length NROW.
!
!    Input, integer I, J, the columns to be swapped.
!
  implicit none

  integer ncol
  integer nrow

  integer i
  integer ia(nrow,ncol)
  integer j
  integer k

  if ( 1 <= i .and. i <= ncol .and. 1 <= j .and. j <= ncol ) then

    do k = 1, nrow
      call i4_swap ( ia(k,k), ia(k,j) )
    end do

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_SWAP(): Fatal error!'
    write ( *, '(a)' ) '  I or J is out of bounds.'
    write ( *, '(a,i8)' ) '  I =    ', i
    write ( *, '(a,i8)' ) '  J =    ', j
    write ( *, '(a,i8)' ) '  NCOL = ', ncol
    stop 1

  end if

  return
end
subroutine i4col_uniq ( m, n, a, nuniq )

!*****************************************************************************80
!
!! i4col_uniq() keeps the unique elements in a sorted I4COL.
!
!  Discussion:
!
!    The array can be sorted into ascending or descending order.
!    The important point is that identical elements must be stored
!    in adjacent positions.
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
!    Input, integer M, the number of rows of A, and the length of
!    a vector of data.
!
!    Input, integer N, the number of columns of A.
!
!    Input/output, real ( kind = rk ) A(M,N).
!    On input, the sorted array of N columns of M-vectors.
!    On output, a sorted array of NUNIQ columns of M-vectors.
!
!    Output, integer NUNIQ, the number of unique columns of A.
!
  implicit none

  integer m
  integer n

  integer a(m,n)
  integer isgn
  integer itest
  integer nuniq

  nuniq = 0

  if ( n <= 0 ) then
    return
  end if

  nuniq = 1

  do itest = 2, n

    call i4col_compare ( m, n, a, itest, nuniq, isgn )

    if ( isgn /= 0 ) then
      nuniq = nuniq + 1
      a(1:m,nuniq) = a(1:m,itest)
    end if

  end do

  return
end
subroutine i4mat_perm ( matrix, n, p )

!*****************************************************************************80
!
!! i4mat_perm() permutes the rows and columns of a square I4MAT.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 July 2000
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
!    Input/output, integer MATRIX(N,N).
!    On input, the matrix to be permuted.
!    On output, the permuted matrix.
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer P(N), the permutation.  P(I) is the new number
!    of row and column I.
!
  implicit none

  integer n

  integer i
  integer i1
  integer is
  integer it
  integer j
  integer j1
  integer j2
  integer k
  integer lc
  integer matrix(n,n)
  integer nc
  integer p(n)

  call perm_cycle ( n, p, is, nc, 1 )

  do i = 1, n

    i1 = - p(i)

    if ( 0 < i1 ) then

      lc = 0

      do

        i1 = p(i1)
        lc = lc + 1

        if ( i1 <= 0 ) then
          exit
        end if

      end do

      i1 = i

      do j = 1, n

        if ( p(j) <= 0 ) then

          j2 = j
          k = lc

          do

            j1 = j2
            it = matrix(i1,j1)

            do

              i1 = abs ( p(i1) )
              j1 = abs ( p(j1) )

              call i4_swap ( matrix(i1,j1), it )

              if ( j1 /= j2 ) then
                cycle
              end if

              k = k - 1

              if ( i1 == i ) then
                exit
              end if

            end do

            j2 = abs ( p(j2) )

            if ( k == 0 ) then
              exit
            end if

          end do

        end if

      end do

    end if

  end do
!
!  Restore the positive signs of the data.
!
  p(1:n) = abs ( p(1:n) )

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
subroutine i4mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! i4mat_print() prints an I4MAT.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    08 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows in A.
!
!    Input, integer N, the number of columns in A.
!
!    Input, integer A(M,N), the matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer m
  integer n

  integer a(m,n)
  integer i
  integer j
  integer jhi
  integer jlo
  character ( len = * ) title

  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) title
  end if

  do jlo = 1, n, 10
    jhi = min ( jlo + 9, n )
    write ( *, '(a)' ) ' '
    write ( *, '(6x,10(i7))' ) ( j, j = jlo, jhi )
    write ( *, '(a)' ) ' '
    do i = 1, m
      write ( *, '(i8,10i7)' ) i, a(i,jlo:jhi)
    end do
  end do

  return
end
subroutine i4mat_row_compare ( m, n, a, row1, row2, result )

!*****************************************************************************80
!
!! i4mat_row_compare() compares two rows of an I4MAT.
!
!  Discussion:
!
!    The rows are compared in the lexicographic sense.  They are equal
!    if every entry is equal.  Otherwise, let I be the first index 
!    where they differ.  Row 1 is less or greater than row 2 as
!    the corresponding indexed values are less or greater.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, number of rows in the matrix.
!
!    Input, integer N, number of columns in the matrix.
!
!    Input, integer A(M,N), the matrix.
!
!    Input, integer ROW1, ROW2, the indices of the two rows to
!    compare.
!
!    Output, integer RESULT:
!    -1, ROW1 < ROW2,
!     0, ROW1 = ROW2,
!    +1, ROW1 > ROW2.
!
  implicit none

  integer m
  integer n

  integer a(m,n)
  integer j
  integer result
  integer row1
  integer row2

  if ( row1 < 1 .or. m < row1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4MAT_ROW_COMPARE(): Fatal error!'
    write ( *, '(a)' ) '  ROW1 index out of bounds.'
    stop 1
  end if

  if ( row2 < 1 .or. m < row2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4MAT_ROW_COMPARE(): Fatal error!'
    write ( *, '(a)' ) '  ROW2 index out of bounds.'
    stop 1
  end if

  result = 0

  do j = 1, n

    if ( a(row1,j) < a(row2,j) ) then
      result = -1 
      return
    else if ( a(row2,j) < a(row1,j) ) then
      result = + 1
      return
    end if

  end do

  return
end
subroutine i4mat_row_sort_d ( m, n, a )

!*****************************************************************************80
!
!! i4mat_row_sort_d() sorts the rows of an I4MAT into descending order.
!
!  Discussion:
!
!    Rows are compared lexicographically.  
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows in A.
!
!    Input, integer N, the number of columns in A.
!
!    Input/output, integer A(M,N).  On input, the M by N matrix
!    to be row sorted.  On output, the row-sorted matrix.
!
  implicit none

  integer m
  integer n

  integer a(m,n)
  integer indx
  integer isgn
  integer row1
  integer row2
!
!  Initialize.
!
  indx = 0
  isgn = 0
  row1 = 0
  row2 = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( m, indx, row1, row2, isgn )
!
!  Interchange the objects.
!
    if ( 0 < indx ) then

      call i4mat_row_swap ( m, n, a, row1, row2 )
!
!  Compare the objects.
!
    else if ( indx < 0 ) then

      call i4mat_row_compare ( m, n, a, row1, row2, isgn )
      isgn = - isgn

    else

      exit

    end if

  end do

  return
end
subroutine i4mat_row_swap ( m, n, a, row1, row2 )

!*****************************************************************************80
!
!! i4mat_row_swap() swaps two rows of an I4MAT.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, number of rows in the matrix.
!
!    Input, integer N, number of columns in the matrix.
!
!    Input/output, integer A(M,N), the matrix.
!
!    Input, integer ROW1, ROW2, the indices of the two rows to
!    swap.
!
  implicit none

  integer m
  integer n

  integer a(m,n)
  integer j
  integer row1
  integer row2

  if ( row1 < 1 .or. m < row1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4MAT_ROW_SWAP(): Fatal error!'
    write ( *, '(a)' ) '  ROW1 index out of bounds.'
    stop 1
  end if

  if ( row2 < 1 .or. m < row2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4MAT_ROW_SWAP(): Fatal error!'
    write ( *, '(a)' ) '  ROW2 index out of bounds.'
    stop 1
  end if

  do j = 1, n
    call i4_swap ( a(row1,j), a(row2,j) )
  end do

  return
end
subroutine i4row_compare ( m, n, a, i, j, isgn )

!*****************************************************************************80
!
!! i4row_compare() compares two rows of an I4ROW.
!
!  Example:
!
!    Input:
!
!  M = 3, N = 4, I = 2, J = 3
!
!  A = (
!    1  2  3  4
!    5  6  7  8
!    9 10 11 12 )
!
!    Output:
!
!  ISGN = -1
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer A(M,N), an array of M rows of vectors of
!    length N.
!
!    Input, integer I, J, the rows to be compared.
!    I and J must be between 1 and M.
!
!    Output, integer ISGN, the results of the comparison:
!    -1, row I < row J,
!     0, row I = row J,
!    +1, row I > row J.
!
  implicit none

  integer m
  integer n

  integer a(m,n)
  integer i
  integer isgn
  integer j
  integer k
!
!  Check that I and J are legal.
!
  if ( i < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4ROW_COMPARE(): Fatal error!'
    write ( *, '(a)' ) '  Row index I is less than 1.'
    write ( *, '(a,i8)' ) '  I = ', i
    stop 1
  else if ( m < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4ROW_COMPARE(): Fatal error!'
    write ( *, '(a)' ) '  Row index I is out of bounds.'
    write ( *, '(a,i8)' ) '  I = ', i
    write ( *, '(a,i8)' ) '  Maximum legal value is M = ', m
    stop 1
  end if

  if ( j < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4ROW_COMPARE(): Fatal error!'
    write ( *, '(a)' ) '  Row index J is less than 1.'
    write ( *, '(a,i8)' ) '  J = ', j
    stop 1
  else if ( m < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4ROW_COMPARE(): Fatal error!'
    write ( *, '(a)' ) '  Row index J is out of bounds.'
    write ( *, '(a,i8)' ) '  J = ', j
    write ( *, '(a,i8)' ) '  Maximum legal value is M = ', m
    stop 1
  end if

  isgn = 0

  if ( i == j ) then
    return
  end if

  k = 1

  do while ( k <= n )

    if ( a(i,k) < a(j,k) ) then
      isgn = - 1
      return
    else if ( a(j,k) < a(i,k) ) then
      isgn = + 1
      return
    end if

    k = k + 1

  end do

  return
end
subroutine i4row_sort_d ( m, n, a )

!*****************************************************************************80
!
!! i4row_sort_d() descending sorts the rows of an I4ROW.
!
!  Definition:
!
!    In lexicographic order, the statement "X < Y", applied to two real
!    vectors X and Y of length M, means that there is some index I, with
!    1 <= I <= M, with the property that
!
!      X(J) = Y(J) for J < I,
!    and
!      X(I) < Y(I).
!
!    In other words, the first time they differ, X is smaller.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    25 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows and columns of A.
!
!    Input/output, integer A(M,N).
!    On input, the array of M rows of N-vectors.
!    On output, the rows of A have been sorted in descending
!    lexicographic order.
!
  implicit none

  integer m
  integer n

  integer a(m,n)
  integer i
  integer indx
  integer isgn
  integer j

  if ( m <= 1 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if
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

    call sort_heap_external ( m, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      call i4row_swap ( m, n, a, i, j )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4row_compare ( m, n, a, i, j, isgn )
      isgn = - isgn

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine i4row_swap ( m, n, a, row1, row2 )

!*****************************************************************************80
!
!! i4row_swap() swaps two rows of an I4ROW.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input/output, integer A(M,N), an array of data.
!
!    Input, integer ROW1, ROW2, the two rows to swap.
!
  implicit none

  integer m
  integer n

  integer a(m,n)
  integer row1
  integer row2
  integer row(n)
!
!  Check.
!
  if ( row1 < 1 .or. m < row1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4ROW_SWAP(): Fatal error!'
    write ( *, '(a)' ) '  ROW1 is out of range.'
    stop 1
  end if

  if ( row2 < 1 .or. m < row2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4ROW_SWAP(): Fatal error!'
    write ( *, '(a)' ) '  ROW2 is out of range.'
    stop 1
  end if

  if ( row1 == row2 ) then
    return
  end if

  row(1:n) = a(row1,1:n)
  a(row1,1:n) = a(row2,1:n)
  a(row2,1:n) = row(1:n)

  return
end
function iset2_compare ( x1, y1, x2, y2 )

!*****************************************************************************80
!
!! iset2_compare() compares two I2 sets.
!
!  Discussion:
!
!    The I2 set (X1,Y1) < (X2,Y2) if
!
!      min ( X1, Y1 ) < min ( X2, Y2 ) or
!      min ( X1, Y1 ) = min ( X2, Y2 ) and max ( X1, Y1 ) < max ( X2, Y2 )
!
!    The I2 set (X1,Y1) = (X2,Y2) if
!
!      min ( X1, Y1 ) = min ( X2, Y2 ) and max ( X1, Y1 ) = max ( X2, Y2 )
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer X1, Y1, the first I2 set.
!
!    Input, integer X2, Y2, the second I2 set.
!
!    Output, character ISET2_COMPARE: '<', '>' or '=' if the first I2 set
!    is less, greater or equal to the second.
!
  implicit none

  integer a1
  integer a2
  integer b1
  integer b2
  character c
  character iset2_compare
  integer x1
  integer x2
  integer y1
  integer y2

  a1 = min ( x1, y1 )
  b1 = max ( x1, y1 )

  a2 = min ( x2, y2 )
  b2 = max ( x2, y2 )

  if ( a1 < a2 ) then
    c = '<'
  else if ( a1 > a2 ) then
    c = '>'
  else if ( b1 < b2 ) then
    c = '<'
  else if ( b1 > b2 ) then
    c = '>'
  else
    c = '='
  end if

  iset2_compare = c

  return
end
subroutine iset2_index_insert_unique ( maxn, n, x, y, indx, &
  xval, yval, ival, ierror )

!*****************************************************************************80
!
!! iset2_index_insert_unique() inserts unique I2 set value in indexed sorted list.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MAXN, the maximum size of the list.
!
!    Input/output, integer N, the size of the list.
!
!    Input/output, integer X(N), Y(N), the list of I2 sets.
!
!    Input, integer INDX(N), the sort index of the list.
!
!    Input, integer XVAL, YVAL, the value to be inserted if it is
!    not already in the list.
!
!    Output, integer IVAL, the index in INDX corresponding to the
!    value XVAL, YVAL.
!
!    Output, integer IERROR, 0 for no error, 1 if an error
!    occurred.
!
  implicit none

  integer maxn

  integer equal
  integer ierror
  integer indx(maxn)
  integer ival
  integer less
  integer more
  integer n
  integer x(maxn)
  integer xval
  integer y(maxn)
  integer yval

  ierror = 0

  if ( n <= 0 ) then

    if ( maxn <= 0 ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ISET2_INDEX_INSERT_UNIQUE(): Fatal error!'
      write ( *, '(a)' ) '  Not enough space to store new data.'
      return
    end if

    n = 1
    x(1) = min ( xval, yval )
    y(1) = max ( xval, yval )
    indx(1) = 1
    ival = 1
    return

  end if
!
!  Does ( XVAL, YVAL ) already occur in the list?
!
  call iset2_index_search ( maxn, n, x, y, indx, xval, yval, &
    less, equal, more )

  if ( equal == 0 ) then

    if ( maxn <= n ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ISET2_INDEX_INSERT_UNIQUE(): Fatal error!'
      write ( *, '(a)' ) '  Not enough space to store new data.'
      return
    end if

    x(n+1) = min ( xval, yval )
    y(n+1) = max ( xval, yval )
    ival = more
    indx(n+1:more+1:-1) = indx(n:more:-1)
    indx(more) = n + 1
    n = n + 1

  else

    ival = equal

  end if

  return
end
subroutine iset2_index_search ( maxn, n, x, y, indx, xval, yval, &
  less, equal, more )

!*****************************************************************************80
!
!! iset2_index_search() searches for an I2 set value in an indexed sorted list.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MAXN, the maximum size of the list.
!
!    Input, integer N, the size of the current list.
!
!    Input, integer X(N), Y(N), the list.
!
!    Input, integer INDX(N), the sort index of the list.
!
!    Input, integer XVAL, YVAL, the value to be sought.
!
!    Output, integer LESS, EQUAL, MORE, the indexes in INDX of the
!    list entries that are just less than, equal to, and just greater
!    than the test value.  If the test value does not occur in the list,
!    then EQUAL is zero.  If the test value is the minimum entry of the
!    list, then LESS is 0.  If the test value is the greatest entry of
!    the list, then MORE is N+1.
!
  implicit none

  integer maxn

  character c
  integer equal
  integer hi
  integer indx(maxn)
  integer less
  integer lo
  integer mid
  integer more
  integer n
  character iset2_compare
  integer x(maxn)
  integer xhi
  integer xlo
  integer xmid
  integer xval
  integer y(maxn)
  integer yhi
  integer ylo
  integer ymid
  integer yval

  if ( n <= 0 ) then
    less = 0
    equal = 0
    more = 0
    return
  end if

  lo = 1
  hi = n

  xlo = x(indx(lo))
  ylo = y(indx(lo))

  xhi = x(indx(hi))
  yhi = y(indx(hi))

  c = iset2_compare ( xval, yval, xlo, ylo )

  if ( c == '<' ) then
    less = 0
    equal = 0
    more = 1
    return
  else if ( c == '=' ) then
    less = 0
    equal = 1
    more = 2
    return
  end if

  c = iset2_compare ( xval, yval, xhi, yhi )

  if ( c == '>' ) then
    less = n
    equal = 0
    more = n + 1
    return
  else if ( c == '=' ) then
    less = n - 1
    equal = n
    more = n + 1
    return
  end if

  do

    if ( lo + 1 == hi ) then
      less = lo
      equal = 0
      more = hi
      return
    end if

    mid = ( lo + hi ) / 2
    xmid = x(indx(mid))
    ymid = y(indx(mid))

    c = iset2_compare ( xval, yval, xmid, ymid )

    if ( c == '=' ) then
      equal = mid
      less = equal - 1
      more = equal + 1
      return
    else if ( c == '<' ) then
      hi = mid
    else if ( c == '>' ) then
      lo = mid
    end if

  end do

  return
end
subroutine i4vec_backtrack ( n, x, indx, k, nstack, stack, maxstack, ncan )

!*****************************************************************************80
!
!! i4vec_backtrack() supervises a backtrack search for an integer vector.
!
!  Discussion:
!
!    The routine tries to construct an integer vector one index at a time,
!    using possible candidates as supplied by the user.
!
!    At any time, the partially constructed vector may be discovered to be
!    unsatisfactory, but the routine records information about where the
!    last arbitrary choice was made, so that the search can be
!    carried out efficiently, rather than starting out all over again.
!
!    First, call the routine with INDX = 0 so it can initialize itself.
!
!    Now, on each return from the routine, if INDX is:
!      1, you've just been handed a complete candidate vector;
!         Admire it, analyze it, do what you like.
!      2, please determine suitable candidates for position X(K).
!         Return the number of candidates in NCAN(K), adding each
!         candidate to the end of STACK, and increasing NSTACK.
!      3, you're done.  Stop calling the routine;
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    24 July 2000
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
!    Input, integer N, the number of positions to be filled in the 
!    vector.
!
!    Input/output, integer X(N), the partial or complete candidate 
!    vector.
!
!    Input/output, integer INDX, a communication flag.
!    On input,
!      0 to start a search.
!    On output:
!      1, a complete output vector has been determined and returned in X(1:N);
!      2, candidates are needed for position X(K);
!      3, no more possible vectors exist.
!
!    Output, integer K, if INDX=2, the current vector index being 
!    considered.
!
!    Input/output, integer NSTACK, the current length of the stack.
!
!    Input, integer STACK(MAXSTACK), a list of all current 
!    candidates for all positions 1 through K.
!
!    Input, integer MAXSTACK, the maximum length of the stack.
!
!    Input/output, integer NCAN(N), lists the current number of 
!    candidates for positions 1 through K.
!
  implicit none

  integer n
  integer maxstack

  integer indx
  integer k
  integer ncan(n)
  integer nstack
  integer stack(maxstack)
  integer x(n)
!
!  If this is the first call, request a candidate for position 1.
!
  if ( indx == 0 ) then
    k = 1
    nstack = 0
    indx = 2
    return
  end if
!
!  Examine the stack.
!
  do
!
!  If there are candidates for position K, take the first available
!  one off the stack, and increment K.
!
!  This may cause K to reach the desired value of N, in which case
!  we need to signal the user that a complete set of candidates
!  is being returned.
!
    if ( 0 < ncan(k) ) then

      x(k) = stack(nstack)
      nstack = nstack - 1

      ncan(k) = ncan(k) - 1

      if ( k /= n ) then
        k = k + 1
        indx = 2
      else
        indx = 1
      end if

      exit
!
!  If there are no candidates for position K, then decrement K.
!  If K is still positive, repeat the examination of the stack.
!
    else

      k = k - 1

      if ( k <= 0 ) then
        indx = 3
        exit
      end if

    end if

  end do

  return
end 
subroutine i4vec_compare ( n, a, b, isgn )

!*****************************************************************************80
!
!! i4vec_compare() compares two integer vectors.
!
!  Discussion:
!
!    The lexicographic ordering is used.
!
!  Example:
!
!    Input:
!
!      A = ( 2, 6, 2 )
!      B = ( 2, 8, 12 )
!
!    Output:
!
!      ISGN = -1
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vectors.
!
!    Input, integer A(N), B(N), the vectors to be compared.
!
!    Output, integer ISGN, the results of the comparison:
!    -1, A is lexicographically less than B,
!     0, A is equal to B,
!    +1, A is lexicographically greater than B.
!
  implicit none

  integer n

  integer a(n)
  integer b(n)
  integer isgn
  integer k

  isgn = 0

  do k = 1, n

    if ( a(k) < b(k) ) then
      isgn = - 1
      return
    else if ( b(k) < a(k) ) then
      isgn = + 1
      return
    end if

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
function i4vec_nonzero ( n, a )

!*****************************************************************************80
!
!! i4vec_nonzero() counts the nonzero entries in an integer vector
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the size of the input array.
!
!    Input, integer A(N), an array.
!
!    Output, integer I4VEC_NONZERO, the number of nonzero entries.
!
  implicit none

  integer n

  integer a(n)
  integer i
  integer i4vec_nonzero

  i4vec_nonzero = 0

  do i = 1, n
    if ( a(i) /= 0 ) then
      i4vec_nonzero = i4vec_nonzero + 1
    end if
  end do

  return
end
subroutine i4vec_order_type ( n, a, order )

!*****************************************************************************80
!
!! i4vec_order_type() determines if I4VEC is (non)strictly ascending/descending.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries of the array.
!
!    Input, integer A(N), the array to be checked.
!
!    Output, integer ORDER, order indicator:
!    -1, no discernable order;
!    0, all entries are equal;
!    1, ascending order;
!    2, strictly ascending order;
!    3, descending order;
!    4, strictly descending order.
!
  implicit none

  integer n

  integer a(n)
  integer i
  integer order
!
!  Search for the first value not equal to A(1).
!
  i = 1

  do

    i = i + 1

    if ( n < i ) then
      order = 0
      return
    end if

    if ( a(1) < a(i) ) then

      if ( i == 2 ) then
        order = 2
      else
        order = 1
      end if

      exit

    else if ( a(i) < a(1) ) then

      if ( i == 2 ) then
        order = 4
      else
        order = 3
      end if

      exit

    end if

  end do
!
!  Now we have a "direction".  Examine subsequent entries.
!
  do while ( i < n )

    i = i + 1

    if ( order == 1 ) then

      if ( a(i) < a(i-1) ) then
        order = -1
        exit
      end if

    else if ( order == 2 ) then

      if ( a(i) < a(i-1) ) then
        order = -1
        exit
      else if ( a(i) == a(i-1) ) then
        order = 1
      end if

    else if ( order == 3 ) then

      if ( a(i-1) < a(i) ) then
        order = -1
        exit
      end if

    else if ( order == 4 ) then

      if ( a(i-1) < a(i) ) then
        order = -1
        exit
      else if ( a(i) == a(i-1) ) then
        order = 3
      end if

    end if

  end do

  return
end
subroutine i4vec_perm_random ( n, a )

!*****************************************************************************80
!
!! i4vec_perm_random() selects a random permutation of an I4VEC.
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
!    Input/output, integer A(N), the vector to be permuted.
!
  implicit none

  integer n

  integer a(n)
  integer i
  integer i4_uniform_ab
  integer j

  do i = 1, n

    j = i4_uniform_ab ( i, n )

    call i4_swap ( a(i), a(j) )

  end do

  return
end
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! i4vec_print() prints an integer vector.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, integer A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer n

  integer a(n)
  integer i
  character ( len = * ) title

  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i8,i10)' ) i, a(i)
  end do

  return
end
subroutine i4vec_reverse ( n, a )

!*****************************************************************************80
!
!! i4vec_reverse() reverses the elements of an integer vector.
!
!  Example:
!
!    Input:
!
!      N = 5, A = ( 11, 12, 13, 14, 15 ).
!
!    Output:
!
!      A = ( 15, 14, 13, 12, 11 ).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input/output, integer A(N), the array to be reversed.
!
  implicit none

  integer n

  integer a(n)
  integer i

  do i = 1, n / 2
    call i4_swap ( a(i), a(n+1-i) )
  end do

  return
end
subroutine i4vec_rotate ( n, m, a )

!*****************************************************************************80
!
!! i4vec_rotate() rotates an object in place.
!
!  Example:
!
!    Input:
!
!      N = 5, M = 2
!      A = ( 1, 2, 3, 4, 5 )
!
!    Output:
!
!      A = ( 4, 5, 1, 2, 3 ).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of objects.
!
!    Input, integer M, the number of positions to the right that
!    each element should be moved.  Elements that shift pass position
!    N "wrap around" to the beginning of the array.
!
!    Input/output, integer A(N), the array to be rotated.
!
  implicit none

  integer n

  integer a(n)
  integer i4_modp
  integer iget
  integer iput
  integer istart
  integer m
  integer mcopy
  integer nset
  integer temp
!
!  Force M to be positive, between 0 and N-1.
!
  mcopy = i4_modp ( m, n )

  if ( mcopy == 0 ) then
    return
  end if

  istart = 0
  nset = 0

  do

    istart = istart + 1

    if ( n < istart ) then
      exit
    end if

    temp = a(istart)
    iget = istart
!
!  Copy the new value into the vacated entry.
!
    do

      iput = iget

      iget = iget - mcopy

      if ( iget < 1 ) then
        iget = iget + n
      end if

      if ( iget == istart ) then
        exit
      end if

      a(iput) = a(iget)
      nset = nset + 1

    end do

    a(iput) = temp
    nset = nset + 1

    if ( n <= nset ) then
      exit
    end if

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
subroutine i4vec_sort_heap_index_d ( n, a, indx )

!*****************************************************************************80
!
!! i4vec_sort_heap_index_d(): indexed heap descending sort of an I4VEC.
!
!  Discussion:
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort 
!    or index the array, or to sort or index related arrays keyed on the 
!    original array.
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      A(INDX(I)), I = 1 to N is sorted,
!
!    or explicitly, by the call
!
!      call I4VEC_PERMUTE ( N, A, INDX )
!
!    after which A(I), I = 1 to N is sorted.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input, integer A(N), an array to be index-sorted.
!
!    Output, integer INDX(N), contains the sort index.  The
!    I-th element of the sorted array is A(INDX(I)).
!
  implicit none

  integer n

  integer a(n)
  integer aval
  integer i
  integer indx(n)
  integer indxt
  integer ir
  integer j
  integer l

  call i4vec_indicator ( n, indx )

  l = n / 2 + 1
  ir = n

  do

    if ( 1 < l ) then

      l = l - 1
      indxt = indx(l)
      aval = a(indxt)

    else

      indxt = indx(ir)
      aval = a(indxt)
      indx(ir) = indx(1)
      ir = ir - 1

      if ( ir == 1 ) then
        indx(1) = indxt
        return
      end if

    end if

    i = l
    j = l + l

    do while ( j <= ir )

      if ( j < ir ) then
        if ( a(indx(j+1)) < a(indx(j)) ) then
          j = j + 1
        end if
      end if

      if ( a(indx(j)) < aval ) then
        indx(i) = indx(j)
        i = j
        j = j + j
      else
        j = ir + 1
      end if

    end do

    indx(i) = indxt

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
subroutine i4vec2_print ( n, a, b, title )

!*****************************************************************************80
!
!! i4vec2_print() prints a pair of integer vectors.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, integer A(N), B(N), the vectors to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer n

  integer a(n)
  integer b(n)
  integer i
  character ( len = * ) title

  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i8,2i10)' ) i, a(i), b(i)
  end do

  return
end
subroutine i4vec2_sort_a ( n, a1, a2 )

!*****************************************************************************80
!
!! i4vec2_sort_a() ascending sorts a vector of pairs of integers.
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

    else if ( indx == 0 ) then

      exit

    end if

  end do

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
subroutine i4vec2_uniq ( n, a1, a2, nuniq )

!*****************************************************************************80
!
!! i4vec2_uniq() keeps the unique elements in a array of pairs of integers.
!
!  Discussion:
!
!    Item I is stored as the pair A1(I), A2(I).
!
!    The items must have been sorted, or at least it must be the
!    case that equal items are stored in adjacent vector locations.
!
!    If the items were not sorted, then this routine will only
!    replace a string of equal values by a single representative.
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
!    Input, integer N, the number of items.
!
!    Input/output, integer A1(N), A2(N).
!    On input, the array of N items.
!    On output, an array of NUNIQ unique items.
!
!    Output, integer NUNIQ, the number of unique items.
!
  implicit none

  integer n

  integer a1(n)
  integer a2(n)
  integer itest
  integer nuniq

  nuniq = 0

  if ( n <= 0 ) then
    return
  end if

  nuniq = 1

  do itest = 2, n

    if ( a1(itest) /= a1(nuniq) .or. a2(itest) /= a2(nuniq) ) then

      nuniq = nuniq + 1

      a1(nuniq) = a1(itest)
      a2(nuniq) = a2(itest)

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
subroutine m_graph_adj_edge_sequence ( adj, nnode, edge_seq )

!*****************************************************************************80
!
!! m_graph_adj_edge_sequence() computes the edge sequence of a multigraph.
!
!  Discussion:
!
!    The edge sequence of a multigraph may be constructed by sorting the
!    entries of each row of the adjacency matrix in descending order, and 
!    then sorting the rows themselves in descending order.
!
!    If two multigraphs are isomorphic, they must have the same edge sequence.
!
!    If two multigraphs have different edge sequences, they cannot be
!    isomorphic.
!
!  Example:
!
!    ADJ = 
!       0 1 2 3
!       1 0 2 0
!       2 2 0 1
!       3 0 1 0
!
!    EDGE_SEQ =
!
!       3 2 1 0
!       3 1 0 0
!       2 2 1 0
!       2 1 0 0
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is an edge from node I to node J.
!
!    Input, integer NNODE, the number of nodes.
!
!    Output, integer EDGE_SEQ(NNODE,NNODE), the degree sequence of
!    the graph.
!
  implicit none

  integer nnode

  integer adj(nnode,nnode)
  integer edge_seq(nnode,nnode)
!
!  Copy the adjacency matrix.
!
  edge_seq(1:nnode,1:nnode) = adj(1:nnode,1:nnode)
!
!  Descending sort the elements of each row.
!
  call i4row_sort_d ( nnode, nnode, edge_seq )
!
!  Sort the rows of the matrix.
!
  call i4mat_row_sort_d ( nnode, nnode, edge_seq )

  return
end
subroutine network_flow_max ( nnode, nedge, iendpt, icpflo, isorce, isink, &
  icut, node_flow )

!*****************************************************************************80
!
!! network_flow_max() finds the maximal flow and a minimal cut in a network.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 July 2000
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
!    Input, integer NEDGE, the number of edges.
!
!    Input/output, integer IENDPT(2,NEDGE), the edges of the
!    network, defined as pairs of nodes.  Each edge should be listed TWICE,
!    the second time in reverse order.  On output, the edges have
!    been reordered, and so the columns of IENDPT have been rearranged.
!
!    Input/output, integer ICPFLO(2,NEDGE).  Capacities and flows.
!    On input, ICPFLO(1,I) is the capacity of edge I.  On output,
!    ICPFLO(2,I) is the flow on edge I and ICPFLO(1,I) has
!    been rearranged to match the reordering of IENDPT.
!
!    Input, integer ISORCE, the designated source node.
!
!    Input, integer ISINK, the designated sink node.
!
!    Output, integer ICUT(NNODE).  ICUT(I) = 1 if node I is in the
!    minimal cut set, otherwise 0.
!
!    Output, integer NODE_FLOW(NNODE).  NODE_FLOW(I) is the value 
!    of the flow through node I.
!
  implicit none

  integer nedge
  integer nnode

  integer i
  integer iarray(nnode)
  integer icpflo(2,nedge)
  integer icut(nnode)
  integer idel
  integer ien1
  integer ien2
  integer iendpt(2,nedge)
  integer indx
  integer ip
  integer iparm
  integer iq
  integer ir
  integer iread
  integer irite
  integer is
  integer isink
  integer isorce
  integer isort
  integer it
  integer itemp
  integer iwork(nnode,2)
  integer j
  integer kz
  integer lst
  integer m
  integer node_flow(nnode)

  iarray(1:nnode) = 0
  idel = 0
 
  do i = 1, nedge

    icpflo(2,i) = 0
    ip = iendpt(1,i)

    if ( ip == isorce ) then
      idel = idel + icpflo(1,i)
    end if

    iarray(ip) = iarray(ip) + 1

  end do
 
  node_flow(isorce) = idel
  is = 1
 
  do i = 1, nnode
    it = iarray(i)
    iarray(i) = is
    iwork(i,1) = is
    is = is + it
  end do
 
  isort = 0
  ien1 = 0
  ien2 = 0
 
10    continue
 
  indx = 0
 
50    continue
 
  call sort_heap_external ( nedge, indx, ien1, ien2, is )

  if ( indx < 0 ) then
 
    is = iendpt(1,ien1) - iendpt(1,ien2)

    if ( is == 0 ) then
      is = iendpt(2,ien1) - iendpt(2,ien2)
    end if

  else if ( 0 < indx ) then
 
    do ir = 1, 2
      itemp           = iendpt(ir,ien1)
      iendpt(ir,ien1) = iendpt(ir,ien2)
      iendpt(ir,ien2) = itemp

      itemp           = icpflo(ir,ien1)
      icpflo(ir,ien1) = icpflo(ir,ien2)
      icpflo(ir,ien2) = itemp
    end do
 
  else
 
    if ( 0 < isort ) then
      return
    end if
 
    do i = 1, nedge
      iq = iendpt(2,i)
      iendpt(1,i) = iwork(iq,1)
      iwork(iq,1) = iwork(iq,1) + 1
    end do
 
    go to 100
 
  end if
 
  go to 50
 
80    continue
 
  iendpt(1,iendpt(1,ien1)) = ien2
  iendpt(1,iendpt(1,ien2)) = ien1
 
  do ir = 1, 2
    itemp           = iendpt(ir,ien1)
    iendpt(ir,ien1) = iendpt(ir,ien2)
    iendpt(ir,ien2) = itemp

    itemp           = icpflo(ir,ien1)
    icpflo(ir,ien1) = icpflo(ir,ien2)
    icpflo(ir,ien2) = itemp
  end do
 
  if ( indx < 0 ) then
    go to 270
  end if

  if ( indx == 0 ) then
    go to 170
  end if

  go to 50
 
100   continue
 
  indx = 0
 
  do i = 1, nnode

    if ( i /= isorce ) then
      node_flow(i) = 0
    end if

    iwork(i,2) = nedge + 1

    if ( i < nnode ) then
      iwork(i,2) = iarray(i+1)
    end if

    icut(i) = 0

  end do
 
  iread = 0
  irite = 1
  iwork(1,1) = isorce
  icut(isorce) = - 1
 
120   continue
 
  iread = iread + 1
 
  if ( iread <= irite ) then
 
    ip = iwork(iread,1)
    lst = iwork(ip,2) - 1
    i = iarray(ip) - 1
 
130     continue
 
    i = i + 1

    if ( lst < i ) then
      go to 120
    end if

    iq = iendpt(2,i)
    idel = icpflo(1,i) - icpflo(2,i)

    if ( icut(iq) /= 0 .or. idel == 0 ) then
      go to 130
    end if
 
    if ( iq /= isink ) then
      irite = irite + 1
      iwork(irite,1) = iq
    end if
 
    icut(iq) = - 1
    go to 130
 
  end if
 
  if ( icut(isink) == 0 ) then
 
    icut(1:nnode) = - icut(1:nnode)
 
    do i = 1, nedge
      ip = iendpt(2,iendpt(1,i))
      if ( icpflo(2,i) < 0 ) then
        node_flow(ip) = node_flow(ip) - icpflo(2,i)
      end if
      iendpt(1,i) = ip
    end do
 
    node_flow(isorce) = node_flow(isink)
    isort = 1
    go to 10
 
  end if
 
  icut(isink) = 1
 
160   continue
 
  iread = iread - 1

  if ( iread == 0 ) then
    go to 180
  end if

  ip = iwork(iread,1)
  ien1 = iarray(ip) - 1
  ien2 = iwork(ip,2) - 1
 
170   continue
 
  if ( ien1 /= ien2 ) then
 
    iq = iendpt(2,ien2)
 
    if ( icut(iq) <= 0 .or. icpflo(1,ien2) == icpflo(2,ien2) ) then
      ien2 = ien2 - 1
      go to 170
    end if
 
    iendpt(2,ien2) = - iq
    icpflo(1,ien2) = icpflo(1,ien2) - icpflo(2,ien2)
    icpflo(2,ien2) = 0
    ien1 = ien1 + 1

    if ( ien1 < ien2 ) then
      go to 80
    end if
 
  end if
 
  if ( iarray(ip) <= ien1 ) then
    icut(ip) = ien1
  end if

  go to 160
 
180   continue
 
  kz = 0
 
  do ir = 1, irite
    if ( 0 < icut(iwork(ir,1)) ) then
      kz = kz + 1
      iwork(kz,1) = iwork(ir,1)
    end if
  end do
 
  indx = - 1
  m = 1
 
200   continue
 
  ip = iwork(m,1)

  if ( 0 < node_flow(ip) ) then
    go to 250
  end if
 
210   continue
 
  m = m + 1

  if ( m <= kz ) then
    go to 200
  end if

  iparm = 0
 
220   continue
 
  m = m - 1
 
  if ( m == 1 ) then
 
    do i = 1, nedge
 
      iq = - iendpt(2,i)
 
      if ( 0 <= iq ) then

        iendpt(2,i) = iq
        j = iendpt(1,i)
        icpflo(1,i) = icpflo(1,i) - icpflo(2,j)

        idel = icpflo(2,i) - icpflo(2,j)
        icpflo(2,i) = idel
        icpflo(2,j) = - idel

      end if
 
    end do
 
    go to 100
 
  end if
 
  ip = iwork(m,1)

  if ( node_flow(ip) < 0 ) then
    go to 220
  end if
 
  if ( node_flow(ip) == 0 ) then

    lst = nedge + 1

    if ( ip < nnode ) then
      lst = iarray(ip+1)
    end if

    i = iwork(ip,2)
    iwork(ip,2) = lst
 
240     continue
 
    if ( i == lst ) then
      go to 220
    end if

    j = iendpt(1,i)
    idel = icpflo(2,j)
    icpflo(2,j) = 0
    icpflo(1,j) = icpflo(1,j) - idel
    icpflo(2,i) = icpflo(2,i) - idel
    i = i + 1
    go to 240
 
  end if
 
  if ( icut(ip) < iarray(ip) ) then
    go to 300
  end if
 
250   continue
 
  i = icut(ip) + 1
 
260   continue

  i = i - 1

  if ( i < iarray(ip) ) then
    go to 290
  end if

  iq = - iendpt(2,i)

  if ( node_flow(iq) < 0 ) then
    go to 260
  end if
 
  idel = icpflo(1,i) - icpflo(2,i)

  if ( node_flow(ip) < idel ) then
    idel = node_flow(ip)
  end if

  icpflo(2,i) = icpflo(2,i) + idel
  node_flow(ip) = node_flow(ip) - idel
  node_flow(iq) = node_flow(iq) + idel
  iparm = 1
  ien1 = iendpt(1,i)
  ien2 = iwork(iq,2) - 1

  if ( ien1 < ien2 ) then
    go to 80
  end if

  if ( ien1 /= ien2 ) then
    go to 280
  end if
 
270   continue
 
  iwork(iq,2) = ien2
 
280   continue
 
  if ( 0 < node_flow(ip) ) then
    go to 260
  end if

  if ( icpflo(1,i) == icpflo(2,i) ) then
    i = i - 1
  end if
 
290   continue
 
  icut(ip) = i

  if ( iparm /= 0 ) then
    go to 210
  end if
 
300   continue
 
  i = iwork(ip,2)
 
310   continue
 
  j = iendpt(1,i)
  idel = icpflo(2,j)

  if ( node_flow(ip) < idel ) then
    idel = node_flow(ip)
  end if

  icpflo(2,j) = icpflo(2,j) - idel
  node_flow(ip) = node_flow(ip) - idel
  iq = iendpt(2,i)
  node_flow(iq) = node_flow(iq) + idel
  i = i + 1

  if ( 0 < node_flow(ip) ) then
    go to 310
  end if

  node_flow(ip) = - 1
  go to 220
 
end
subroutine node_order_print ( nnode, order )

!*****************************************************************************80
!
!! node_order_print() prints out a node ordering.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    29 May 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, integer ORDER(NNODE), the node ordering.  ORDER(1) is
!    the label of the node which is to be taken as the first node, and so on.
!
  implicit none

  integer nnode

  integer i
  integer ihi
  integer ilo
  integer inc
  integer order(nnode)

  inc = 15

  do ilo = 1, nnode, inc

    ihi = min ( ilo + inc - 1, nnode )

    write ( *, '(a)' ) ' '
    write ( *, '(a6,4x,15i4)' ) 'Order:', ( i, i = ilo, ihi )
    write ( *, '(a6,4x,15i4)' ) 'Label:', order(ilo:ihi)

  end do

  return
end
subroutine node_relax ( cor3, cor3_new, cor3_nabe, face, face_order, max_cor3, &
  max_face, max_order, num_cor3, num_face )

!*****************************************************************************80
!
!! node_relax() smooths a shape by an averaging operation on the node positions.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) COR3(3,MAXCOR3), the coordinates of the nodes.
!
!    Output, real ( kind = rk ) COR3_NEW(3,MAXCOR3), the new, averaged
!    coordinates of the nodes.
!
!    Workspace, integer COR3_NABE(MAXCOR3).  On output, COR3_NABE(I)
!    will contain the number of node neighbors of node I.
!
!    Input, integer FACE(MAX_ORDER,MAX_FACE), describes the faces.
!    FACE(I,J) is the index of the I-th node in face J. 
!
!    Input, integer FACE_ORDER(MAX_FACE), is the number of nodes
!    making up each face.
!
!    Input, integer MAX_FACE, the maximum number of faces.
!
!    Input, integer MAX_ORDER, is the maximum number of nodes that
!    can make up a face, required to dimension FACE.
!
!    Input, integer NUM_FACE, the number of faces.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer max_cor3
  integer max_face
  integer max_order

  real ( kind = rk ) cor3(3,max_cor3)
  real ( kind = rk ) cor3_new(3,max_cor3)
  integer cor3_nabe(max_cor3)
  integer face(max_order,max_face)
  integer face_order(max_face)
  integer icor3
  integer iface
  integer inode
  integer ivert
  integer jnode
  integer num_cor3
  integer num_face
!
!  COR3_NEW will contain the new averaged coordinates.
!
  cor3_nabe(1:num_cor3) = 0
  cor3_new(1:3,1:num_cor3) = 0.0D+00
!
!  Consider each edge.  Essentially, the edge (I,J) is a signal to
!  add the old coordinates of I to the new J coordinates, and vice versa.
!
!  Because we are using a face representation, many, perhaps all the
!  edges, will show up repeatedly, probably twice.  To keep the algorithm
!  simple, for now we will simply use an edge every time it shows up
!  in a face, which means that edges that occur in multiple faces
!  will be weighted more.
!
  do iface = 1, num_face

    inode = face(face_order(iface),iface)

    do ivert = 1, face_order(iface)
      jnode = inode
      inode = face(ivert,iface)
      cor3_nabe(inode) = cor3_nabe(inode) + 1
      cor3_nabe(jnode) = cor3_nabe(jnode) + 1
      cor3_new(1:3,jnode) = cor3_new(1:3,jnode) + cor3(1:3,inode)
      cor3_new(1:3,inode) = cor3_new(1:3,inode) + cor3(1:3,jnode)
    end do

  end do
!
!  Copy the new into the old.
!
  do icor3 = 1, num_cor3

    if ( cor3_nabe(icor3) /= 0 ) then
      cor3_new(1:3,icor3) = cor3_new(1:3,icor3) &
        / real ( cor3_nabe(icor3), kind = rk )
    end if

  end do

  return
end
subroutine nodes_to_ps ( plotxmin2, plotymin2, alpha, iunit, nnode, x, y, &
  xmin, ymin )

!*****************************************************************************80
!
!! nodes_to_ps() writes subplot nodes to a PostScript file.
!
!  Discussion:
!
!    A small filled circle is placed at each node.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer PLOTXMIN2, PLOTYMIN2, the Postscript point
!    corresponding to the physical point XMIN, YMIN.
!
!    Input, real ( kind = rk ) ALPHA, the physical-to-Postscript scale factor.
!
!    Input, integer IUNIT, the output FORTRAN unit.
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, real ( kind = rk ) X(NNODE), Y(NNODE), the coordinates of points.
!
!    Input, real ( kind = rk ) XMIN, YMIN, the coordinates of the physical
!    origin.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nnode

  real ( kind = rk ) alpha
  integer i
  integer iunit
  integer plotxmin2
  integer plotymin2
  integer px1
  integer py1
  integer, parameter :: rad = 10
  real ( kind = rk ) x(nnode)
  real ( kind = rk ) xmin
  real ( kind = rk ) y(nnode)
  real ( kind = rk ) ymin
!
!  Draw points.
!
  do i = 1, nnode

    px1 = plotxmin2 + nint ( alpha * ( x(i) - xmin ) )
    py1 = plotymin2 + nint ( alpha * ( y(i) - ymin ) )

    write ( iunit, '(3i4,a)' ) px1, py1, rad, ' 0 360 arc closepath fill'

  end do

  return
end
subroutine object_build ( face, face_object, face_order, face_rank, face_tier, &
  max_order, num_face, num_object )

!*****************************************************************************80
!
!! object_build() builds edge-connected "objects" out of polygonal faces.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer FACE(MAX_ORDER,NUM_FACE), describes the faces.
!    FACE(I,J) is the index of the I-th node in face J.  It is best
!    if the nodes of all faces are listed in counterclockwise order.
!
!    Output, integer FACE_OBJECT(NUM_FACE), describes the objects.
!    FACE_OBJECT(I) is the index of the edge-connected "object" to 
!    which face I belongs.
!
!    Input, integer FACE_ORDER(NUM_FACE), is the number of nodes
!    making up each face.
!
!    Output, integer FACE_RANK(NUM_FACE), is an ordered list of
!    faces.  FACE_RANK(1) is the index of the face in the first tier of the 
!    first object, followed by second tier faces, and so on until
!    object one is complete.  Object two follows, and so on.
!
!    Output, integer FACE_TIER(NUM_FACE).  FACE_TIER(I) is the
!    "tier" of face I in its object.  The base of the object has tier 1,
!    the neighbors of the base have tier 2, and so on.
!
!    Input, integer MAX_ORDER, is the maximum number of nodes that
!    can make up a face, required to dimension FACE.
!
!    Input, integer NUM_FACE, the number of faces.
!
!    Output, integer NUM_OBJECT, the number of objects.
!
  implicit none

  integer max_order
  integer num_face

  integer base
  integer face(max_order,num_face)
  integer face_object(num_face)
  integer face_order(num_face)
  integer face_rank(num_face)
  integer face_tier(num_face)
  integer i
  integer iface
  integer ihi
  integer ihi_next
  integer ilo
  integer ilo_next
  integer irank
  integer jface
  integer num_object
  integer tier
  integer touch
!
!  Initialization.
!
  num_object = 0

  if ( num_face <= 0 ) then
    return
  end if

  face_object(1:num_face) = 0
  face_rank(1:num_face) = 0
  face_tier(1:num_face) = 0

  irank = 0

  base = 1
!
!  Begin the next object, from face BASE.
!
10    continue

  tier = 1

  num_object = num_object + 1
  irank = irank + 1

  face_rank(irank) = base
  face_tier(base) = tier
  face_object(base) = num_object

  ilo = irank
  ihi = irank
!
!  Begin the next tier of faces, which are neighbors of faces we
!  found in the previous tier.
!
20    continue

  tier = tier + 1

  ilo_next = ihi + 1
  ihi_next = ihi

  do jface = 1, num_face

    if ( face_tier(jface) == 0 ) then

      do i = ilo, ihi

        iface = face_rank(i)

        call face_touch ( face, face_order, max_order, num_face, iface, &
          jface, touch )

        if ( touch /= 0 ) then
          ihi_next = ihi_next + 1
          irank = irank + 1
          face_rank(irank) = jface
          face_tier(jface) = tier
          face_object(jface) = num_object
          exit
        end if

      end do

    end if

  end do

  if ( ilo_next <= ihi_next ) then
    ilo = ilo_next
    ihi = ihi_next
    go to 20
  end if
!
!  No neighbors were found, so this object is complete.  
!  Search for an unused face, which will be the base of the next object.
!
  do iface = 1, num_face

    if ( face_tier(iface) == 0 ) then
      base = iface
      go to 10
    end if

  end do

  return
end
subroutine perm_cycle ( n, isig, isgn, ncycle, iopt )

!*****************************************************************************80
!
!! perm_cycle() analyzes a permutation.
!
!  Discussion:
!
!    The routine will count cycles, find the sign of a permutation,
!    and tag a permutation.
!
!  Example:
!
!    Input:
!
!      N = 9
!      IOPT = 1
!      ISIG = 2, 3, 9, 6, 7, 8, 5, 4, 1
!
!    Output:
!
!      NCYCLE = 3
!      ISGN = +1
!      ISIG = -2, 3, 9, -6, -7, 8, 5, 4, 1
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 July 2000
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
!    Input, integer N, the number of objects being permuted.
!
!    Input/output, integer ISIG(N).  On input, ISIG describes a
!    permutation, in the sense that entry I is to be moved to ISIG(I).
!    If IOPT = 0, then ISIG will not be changed by this routine.
!    If IOPT = 1, then on output, ISIG will be "tagged".  That is,
!    one element of every cycle in ISIG will be negated.  In this way,
!    a user can traverse a cycle by starting at any entry I1 of ISIG
!    which is negative, moving to I2 = ABS(ISIG(I1)), then to
!    ISIG(I2), and so on, until returning to I1.
!
!    Output, integer ISGN, the "sign" of the permutation, which is
!    +1 if the permutation is even, -1 if odd.  Every permutation
!    may be produced by a certain number of pairwise switches.
!    If the number of switches is even, the permutation itself is
!    called even.
!
!    Output, integer NCYCLE, the number of cycles in the
!    permutation.
!
!    Input, integer IOPT, requests tagging.
!    0, the permutation will not be tagged.
!    1, the permutation will be tagged.
!
  implicit none

  integer n

  integer i
  integer i1
  integer i2
  integer iopt
  integer is
  integer isgn
  integer isig(n)
  integer ncycle

  is = 1
  ncycle = n

  do i = 1, n

    i1 = isig(i)

    do while ( i < i1 )
      ncycle = ncycle - 1
      i2 = isig(i1)
      isig(i1) = - i2
      i1 = i2
    end do

    if ( iopt /= 0 ) then
      is = - isign ( 1, isig(i) )
    end if

    isig(i) = isign ( isig(i), is )

  end do

  isgn = 1 - 2 * mod ( n-ncycle, 2 )

  return
end
subroutine perm_free ( ipart, npart, ifree, nfree )

!*****************************************************************************80
!
!! perm_free() reports the number of unused items in a partial permutation.
!
!  Discussion:
!
!    It is assumed that the N objects being permuted are the integers
!    from 1 to N, and that IPART contains a "partial" permutation, that
!    is, the NPART entries of IPART represent the beginning of a
!    permutation of all N items.
!
!    The routine returns in IFREE the items that have not been used yet.
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
!    Input, integer IPART(NPART), the partial permutation, which
!    should contain, at most once, some of the integers between 1 and
!    NPART+NFREE.
!
!    Input, integer NPART, the number of entries in IPART.  NPART
!    may be 0.
!
!    Output, integer IFREE(NFREE), the integers between 1 and 
!    NPART+NFREE that were not used in IPART.
!
!    Input, integer NFREE, the number of integers that have not
!    been used in IPART.  This is simply N - NPART.  NFREE may be zero.
!
  implicit none

  integer nfree
  integer npart

  integer i
  integer ifree(nfree)
  integer ipart(npart)
  integer j
  integer k
  integer n

  n = npart + nfree

  if ( npart < 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_FREE(): Fatal error!'
    write ( *, '(a)' ) '  NPART < 0.'
    stop 1

  else if ( npart == 0 ) then

    call i4vec_indicator ( n, ifree )

  else if ( nfree < 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_FREE(): Fatal error!'
    write ( *, '(a)' ) '  NFREE < 0.'
    stop 1

  else if ( nfree == 0 ) then

    return

  else

    k = 0

    do i = 1, n

      do j = 1, npart
        if ( ipart(j) == i ) then
          go to 10
        end if
      end do

      k = k + 1

      if ( nfree < k ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PERM_FREE(): Fatal error!'
        write ( *, '(a)' ) '  The partial permutation is illegal.'
        write ( *, '(a)' ) '  It should contain, at most once, some of'
        write ( *, '(a,i8)' ) '  the integers between 1 and ', n
        stop 1
      end if

      ifree(k) = i

10    continue

    end do

  end if

  return
end
subroutine perm_inc ( iperm, ipos, n )

!*****************************************************************************80
!
!! perm_inc() "increments" a permutation to get the "next" one.
!
!  Discussion:
!
!    The routine is given IPERM, a permutation of the numbers from 1 to N,
!    and a position IPOS between 1 and N.
!
!    It returns the next permutation in the dictionary order which
!    comes after all permutations beginning IPERM(1) through IPERM(IPOS).
!
!  Example:
!
!             PERM              IPOS
!
!    Input    123456789         7
!    Output   123456798         7
!
!    Input    123456789         9
!    Output   213456789         0
!
!    Input    134826795         3
!    Output   134925678         3
!
!    Input    134826795         0
!    Output   123456789         0
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
!    Input/output, integer IPERM(N).
!    On input, the current permutation.
!    On output, the "incremented" permutation.
!
!    Input/output, integer IPOS.
!    On input, IPOS is the location of the end of the string of
!    "digits" in IPERM that form the test string.  That is, the
!    new permutation to be computed must be the very next one,
!    in dictionary order, which succeeds all strings whose first
!    IPOS digits agree with the input value of IPERM.
!
!    On output, IPOS is the position of the last digit of the output
!    value of IPERM which agrees with the input value of IPERM.
!
!    Input, integer N, is the number of entries in IPERM.
!
  implicit none

  integer n

  integer ipcopy
  integer iperm(n)
  integer ipos
  integer j
  integer k
  integer new

  if ( ipos == 0 ) then
    ipos = n
    call i4vec_indicator ( n, iperm )
    return
  end if
 
  ipcopy = ipos

10    continue
!
!  To get the next permutation, we need to increment the IPOS+1 "digit".
!
!  We do this by finding, if possible, a digit in positions IPOS+2
!  through N that is just larger than the current value IPOS+1 digit.
!  If we find such a digit, it becomes the IPOS+1 digit, and the
!  remaining values are sorted into increasing order.
!
  new = 0
  do j = ipcopy+2, n
    if ( new == 0 ) then
      if ( iperm(ipcopy+1) < iperm(j) ) then
        new = j
      end if
    else
      if ( iperm(ipcopy+1) < iperm(j) .and. iperm(j) < iperm(new) ) then
        new = j
      end if
    end if
  end do
!
!  There is a next candidate that agrees with IPERM through entry I.
!  Swap entries IPOS+1 and NEW, and sort the entries (IPOS+2,...,N).
!
!  The output value of IPOS equals the input value.
!
  if ( new /= 0 ) then

    call i4_swap ( iperm(new), iperm(ipcopy+1) )
 
    do j = ipcopy+2, n
 
      do k = j+1, n
        if ( iperm(k) < iperm(j) ) then
          call i4_swap ( iperm(j), iperm(k) )
        end if
      end do
 
    end do
    return
  end if
!
!  There is no next candidate that agrees with IPERM through entry 
!  IPOS.  Can we decrease IPOS and try for a next candidate that way?
!
  if ( 0 < ipcopy ) then
    ipcopy = ipcopy - 1
    go to 10
  end if
!
!  IPOS is now zero.  There is no successor to the current permutation,
!  so we start again at the first permutation.
!
  ipos = 0
  call i4vec_indicator ( n, iperm )
 
  return
end
subroutine perm_inv ( n, isig )

!*****************************************************************************80
!
!! perm_inv() inverts a permutation.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 July 2000
!
!  Parameters:
!
!    Input, integer N, the number of objects being permuted.
!
!    Input/output, integer ISIG(N).
!
!    On input, ISIG describes a permutation.
!
!    ISIG is used to represent a permutation by the convention that
!    the permutation maps the letter I to ISIG(I).  Thus, if ISIG
!    contains the values (4, 1, 3, 2), then the permutation
!    represented permutes 1 to 4, 2 to 1, 3 to 3, and 4 to 2.
!
!    On output, ISIG describes the inverse permutation
!
  implicit none

  integer n

  integer i
  integer i0
  integer i1
  integer i2
  integer is
  integer isig(n)

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_INV(): Fatal error!'
    write ( *, '(a,i8)' ) '  Input value of N = ', n
    stop 1
  end if

  is = 1

  do i = 1, n

    i1 = isig(i)

    do while ( i < i1 )
      i2 = isig(i1)
      isig(i1) = - i2
      i1 = i2
    end do

    is = - isign ( 1, isig(i) )
    isig(i) = isign ( isig(i), is )

  end do

  do i = 1, n

    i1 = - isig(i)

    if ( 0 <= i1 ) then

      i0 = i

      do

        i2 = isig(i1)
        isig(i1) = i0

        if ( i2 < 0 ) then
          exit
        end if

        i0 = i1
        i1 = i2

      end do

    end if

  end do

  return
end
subroutine perm_next ( n, iarray, more, even )

!*****************************************************************************80
!
!! perm_next() computes all of the permutations on N objects, one at a time.
!
!  Discussion:
!
!    If the routine is called with MORE = .TRUE., any permutation in
!    IARRAY, and EVEN = .TRUE., then the successor of the input
!    permutation will be produced, unless IARRAY is the last permutation
!    on N letters, in which case IARRAY(1) will be set to 0 on return.
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
!    Input, integer N, the number of objects being permuted.
!
!    Input/output, integer IARRAY(N).
!
!    If MORE is .TRUE., then IARRAY is assumed to contain the
!    "previous" permutation, and on IARRAY(I) is the value
!    of the I-th object under the next permutation.
!
!    Otherwise, IARRAY(I) will be set to the "first" permutation.
!
!    Input/output, logical MORE.
!
!    Set MORE to FALSE before first calling this routine.
!
!    MORE will be reset to .TRUE. and a permutation will be returned.
!
!    Each new call produces a new permutation until
!    MORE is returned .FALSE.
!
!    Output, logical EVEN.
!
!    EVEN is .TRUE. if the output permutation is even, that is,
!    involves an even number of transpositions.
!
!    EVEN is .FALSE. otherwise.
!
  implicit none

  integer n

  integer i
  integer i1
  integer ia
  integer iarray(n)
  integer id
  integer is
  integer j
  integer l
  integer m
  logical more
  logical even

  if ( .not. more ) then

    call i4vec_indicator ( n, iarray )

    more = .true.
    even = .true.

    if ( n == 1 ) then
      more = .false.
      return
    end if

    if ( iarray(n) /= 1 .or. iarray(1) /= 2 + mod ( n, 2 ) ) then
      return
    end if

    do i = 1, n - 3
      if ( iarray(i+1) /= iarray(i)+1 ) then
        return
      end if
    end do

    more = .false.

  else

    if ( n == 1 ) then
      iarray(1) = 0
      more = .false.
      return
    end if

    if ( even ) then

      ia = iarray(1)
      iarray(1) = iarray(2)
      iarray(2) = ia
      even = .false.

      if ( iarray(n) /= 1 .or. iarray(1) /= 2 + mod ( n, 2 ) ) then
        return
      end if

      do i = 1, n - 3
        if ( iarray(i+1) /= iarray(i)+1 ) then
          return
        end if
      end do

      more = .false.
      return

    else

      is = 0

      do i1 = 2, n

        ia = iarray(i1)
        i = i1-1
        id = 0

        do j = 1, i
          if ( ia < iarray(j) ) then
            id = id + 1
          end if
        end do

        is = id + is

        if ( id /= i * mod ( is, 2 ) ) then
          go to 10
        end if

      end do

      iarray(1) = 0
      more = .false.
      return

    end if

10      continue

    m = mod ( is+1, 2 ) * (n+1)

    do j = 1, i

      if ( isign(1,iarray(j)-ia) /= isign(1,iarray(j)-m) ) then
        m = iarray(j)
        l = j
      end if

    end do

    iarray(l) = ia
    iarray(i1) = m
    even = .true.

  end if

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
subroutine poly ( n, iarray, ix0, iopt, ival )

!*****************************************************************************80
!
!! poly() performs operations on polynomials in power or factorial form.
!
!  Definition:
!
!    The power sum form of a polynomial is
!
!      P(X) = A1+A2*X+A3*X^2+...+(AN+1)*X^N
!
!    The Taylor expansion at C has the form
!
!      P(X) = A1+A2*(X-C)+A3*(X-C)^2+...+(AN+1)*(X-C)^N
!
!    The factorial form of a polynomial is
!
!      P(X) = A1+A2*X+A3*(X)*(X-1)+A4*(X)*(X-1)*(X-2)+...
!        +(AN+1)*(X)*(X-1)*...*(X-N+1)
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 April 1999
!
!  Parameters:
!
!    Input, integer N, the number of coefficients in the polynomial
!    (in other words, the polynomial degree + 1)
!
!    Input, integer IOPT, a flag describing which algorithm is to
!    be carried out:
!
!    -3: Reverse Stirling.  Input the coefficients of
!    the polynomial in factorial form, output them in
!    power sum form.
!
!    -2: Stirling.  Input the coefficients in power sum
!    form, output them in factorial form.
!
!    -1: Evaluate a polynomial which has been input
!    in factorial form.
!
!    0:  Evaluate a polynomial input in power sum form.
!
!    1 or more:  Given the coefficients of a polynomial in
!    power sum form, compute the first IOPT coefficients of
!    the polynomial in Taylor expansion form.
!
!    Input, integer IX0, for IOPT = -1, 0, or positive, the value
!    X of the argument at which the polynomial is to be evaluated, or the
!    Taylor expansion is to be carried out.
!
!    Output, integer IVAL, for IOPT = -1 or 0, the value of the
!    polynomial at the point IX0.
!
!    Input, integer IARRAY(N).  Contains the coefficients of the
!    polynomial.  Depending on the option chosen, these coefficients may
!    be overwritten by those of a different form of the polynomial.
!
  implicit none

  integer n

  integer i
  integer iarray(n)
  integer ieps
  integer iopt
  integer ival
  integer iw
  integer ix0
  integer iz
  integer m
  integer n1

  n1 = min ( n, iopt )
  n1 = max ( 1, n1 )
 
  if ( iopt < -1 ) then
    n1 = n
  end if
 
  ieps = mod ( max ( -iopt, 0 ), 2 )
 
  iw = -n * ieps

  if ( -2 < iopt ) then
    iw = iw + ix0
  end if
 
  do m = 1, n1
 
    ival = 0
    iz = iw
 
    do i = m, n
      iz = iz + ieps
      ival = iarray(n+m-i) + iz * ival
      if ( iopt /= 0 .and. iopt /= -1 ) then
        iarray(n+m-i) = ival
      end if
    end do
 
    if ( iopt < 0 ) then
      iw = iw + 1
    end if
 
  end do
 
  return
end
subroutine poly_to_tri ( face, ierror, max_face, max_vert, num_face, num_vert )

!*****************************************************************************80
!
!! poly_to_tri() converts a collection of polygons into a collection of triangles.
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
!    Input/output, integer FACE(MAX_VERT,MAX_FACE), describes 
!    the faces.  FACE(I,J) is the I-th node associated with the J-th face.  
!    This array is updated on return.
!
!    Output, integer IERROR, error flag.
!    0, no error.
!    1, the algorithm failed because MAX_FACE was too small.
!    2, the algorithm failed because there were faces of order < 3.
!    3, the algorithm failed because there were faces of order > MAX_VERT.
!
!    Input, integer MAX_FACE, the maximum number of faces allowed.
!
!    Input, integer MAX_VERT, the maximum number of nodes allowed 
!    per face.
!
!    Input/output, integer NUM_FACE, the number of faces.  
!    This value is updated on return.
!
!    Input/output, integer NUM_VERT(MAX_FACE), the number of nodes
!    associated with each face.  On successful return, every entry of
!    this array will be 3.
!
  implicit none

  integer max_face
  integer max_vert

  integer face(max_vert,max_face)
  integer ierror
  integer iface
  integer iface_old
  integer ivert
  integer k
  integer num_face
  integer num_face2
  integer num_vert(max_face)

  ierror = 0
  num_face2 = 0

  do iface = 1, num_face

    if ( num_vert(iface) < 3 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POLY_TO_TRI(): Fatal error!'
      write ( *, '(a,i8)' ) '  Illegal face ', iface
      write ( *, '(a,i8)' ) '  Number of nodes is ', num_vert(iface)
      ierror = 2
      return
    else if ( max_vert < num_vert(iface) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POLY_TO_TRI(): Fatal error!'
      write ( *, '(a,i8)' ) '  Illegal face ', iface
      write ( *, '(a,i8)' ) '  Number of nodes is ', num_vert(iface)
      write ( *, '(a,i8)' ) '  MAX_VERT is ', max_vert
      ierror = 3
      return
    end if

    do ivert = 3, num_vert(iface)
      num_face2 = num_face2 + 1
    end do

  end do

  if ( max_face < num_face2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POLY_TO_TRI(): Fatal error!'
    write ( *, '(a)' ) '  MAX_FACE is too small to replace all faces'
    write ( *, '(a)' ) '  by triangles.'
    write ( *, '(a,i8)' ) '  MAX_FACE = ', max_face
    write ( *, '(a,i8)' ) '  NUM_FACE2 = ', num_face2
    ierror = 1
    return
  end if

  iface_old = num_face
  k = num_vert(iface_old)

  do iface = num_face2, 1, -1

    if ( k < 3 ) then
      iface_old = iface_old - 1
      k = num_vert(iface_old)
    end if

    num_vert(iface) = 3
    face(1,iface) = face(1,iface_old)
    do ivert = 2, 3
      face(ivert,iface) = face(k+ivert-3,iface_old)
    end do

    k = k - 1

  end do

  num_face = num_face2

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
function pythag ( a, b )

!*****************************************************************************80
!
!! pythag() computes SQRT ( A^2 + B^2 ) carefully.
!
!  Discussion:
!
!    The formula
!
!      PYTHAG = sqrt ( A^2 + B^2 )
!
!    is reasonably accurate, but the formula can actually fail if
!    for example, A^2 is larger than the machine overflow.  The
!    formula can lose most of its accuracy if the sum of the squares
!    is very large or very small.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 March 2000
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow, 
!    Y Ikebe, V Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, the two legs of a right triangle.
!
!    Output, real ( kind = rk ) PYTHAG, the length of the hypotenuse.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) p
  real ( kind = rk ) pythag
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) t
  real ( kind = rk ) u

  p = max ( abs ( a ), abs ( b ) )

  if ( p /= 0.0D+00 ) then

    r = ( min ( abs ( a ), abs ( b ) ) / p )**2

   10   continue

    t = 4.0D+00 + r

    if ( t /= 4.0D+00 ) then
      s = r / t
      u = 1.0D+00 + 2.0D+00 * s
      p = u * p
      r = ( s / u )**2 * r
      go to 10
    end if

  end if

  pythag = p

  return
end
subroutine r8_swap ( x, y )

!*****************************************************************************80
!
!! r8_swap() swaps two double precision values.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = rk ) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) x
  real ( kind = rk ) y
  real ( kind = rk ) z

  z = x
  x = y
  y = z

  return
end
subroutine r8col_find ( m, n, a, x, i4col )

!*****************************************************************************80
!
!! r8col_find() seeks a table column equal to a real vector.
!
!  Example:
!
!    Input:
!
!      M = 3,
!      N = 4,
!
!      A = (
!        1.  2.  3.  4.
!        5.  6.  7.  8.
!        9. 10. 11. 12. )
!
!      x = ( 3.,
!            7.,
!           11. )
!
!    Output:
!
!      I4COL = 3
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = rk ) A(M,N), a table of numbers, regarded as
!    N columns of vectors of length M.
!
!    Input, real ( kind = rk ) X(M), a vector to be matched with a column of A.
!
!    Output, integer I4COL, the index of the first column of A
!    which exactly matches every entry of X, or 0 if no match
!    could be found.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  integer i
  integer i4col
  integer j
  real ( kind = rk ) x(m)

  i4col = 0

  do j = 1, n

    i4col = j

    do i = 1, m
      if ( x(i) /= a(i,j) ) then
        i4col = 0
        exit
      end if
    end do

    if ( i4col /= 0 ) then
      return
    end if

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
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! r8vec_print() prints an R8VEC.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, real ( kind = rk ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  integer i
  character ( len = * ) title

  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i8,g14.6)' ) i, a(i)
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
subroutine r8vec2_print ( n, a1, a2, title )

!*****************************************************************************80
!
!! r8vec2_print() prints a pair of R8VEC's.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, real ( kind = rk ) A1(N), A2(N), the vectors to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a1(n)
  real ( kind = rk ) a2(n)
  integer i
  character ( len = * ) title

  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i8,2g14.6)' ) i, a1(i), a2(i)
  end do

  return
end
function r8vec3_compare ( x1, y1, z1, x2, y2, z2 )

!*****************************************************************************80
!
!! r8vec3_compare() compares two R8VEC's.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) X1, Y1, Z1, the first vector.
!
!    Input, real ( kind = rk ) X2, Y2, Z2, the second vector.
!
!    Output, character R8VEC3_COMPARE: '<', '>' or '=' if the first vector
!    is less, greater or equal to the second.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character c
  character r8vec3_compare
  real ( kind = rk ) x1
  real ( kind = rk ) x2
  real ( kind = rk ) y1
  real ( kind = rk ) y2
  real ( kind = rk ) z1
  real ( kind = rk ) z2

  if ( x1 < x2 ) then
    c = '<'
  else if ( x1 > x2 ) then
    c = '>'
  else if ( y1 < y2 ) then
    c = '<'
  else if ( y1 > y2 ) then
    c = '>'
  else if ( z1 < z2 ) then
    c = '<'
  else if ( z1 > z2 ) then
    c = '>'
  else
    c = '='
  end if

  r8vec3_compare = c

  return
end
subroutine r8vec3_index_insert_unique ( maxn, n, x, y, z, indx, &
  xval, yval, zval, ival, ierror )

!*****************************************************************************80
!
!! r8vec3_index_insert_unique() inserts unique value in an indexed sorted list.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MAXN, the maximum size of the list.
!
!    Input/output, integer N, the size of the list.
!
!    Input/output, real ( kind = rk ) X(N), Y(N), Z(N), the list of R3 vectors.
!
!    Input, integer INDX(N), the sort index of the list.
!
!    Input, real ( kind = rk ) XVAL, YVAL, ZVAL, the value to be inserted 
!    if it is not already in the list.
!
!    Output, integer IVAL, the index in INDX corresponding to the
!    value XVAL, YVAL, ZVAL.
!
!    Output, integer IERROR, 0 for no error, 1 if an error
!    occurred.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxn

  integer equal
  integer ierror
  integer indx(maxn)
  integer ival
  integer less
  integer more
  integer n
  real ( kind = rk ) x(maxn)
  real ( kind = rk ) xval
  real ( kind = rk ) y(maxn)
  real ( kind = rk ) yval
  real ( kind = rk ) z(maxn)
  real ( kind = rk ) zval

  ierror = 0

  if ( n <= 0 ) then

    if ( maxn <= 0 ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8VEC3_INDEX_INSERT_UNIQUE(): Fatal error!'
      write ( *, '(a)' ) '  Not enough space to store new data.'
      return
    end if

    n = 1
    x(1) = xval
    y(1) = yval
    z(1) = zval
    indx(1) = 1
    ival = 1
    return

  end if
!
!  Does ( XVAL, YVAL, ZVAL ) already occur in ( X, Y, Z)?
!
  call r8vec3_index_search ( maxn, n, x, y, z, indx, xval, yval, zval, &
    less, equal, more )

  if ( equal == 0 ) then

    if ( maxn <= n ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8VEC3_INDEX_INSERT_UNIQUE(): Fatal error!'
      write ( *, '(a)' ) '  Not enough space to store new data.'
      return
    end if

    x(n+1) = xval
    y(n+1) = yval
    z(n+1) = zval
    ival = more
    indx(n+1:more+1:-1) = indx(n:more:-1)
    indx(more) = n + 1
    n = n + 1

  else

    ival = equal

  end if

  return
end
subroutine r8vec3_index_search ( maxn, n, x, y, z, indx, xval, yval, &
  zval, less, equal, more )

!*****************************************************************************80
!
!! r8vec3_index_search() searches for an R3 value in an indexed sorted list.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MAXN, the maximum size of the list.
!
!    Input, integer N, the size of the current list.
!
!    Input, real ( kind = rk ) X(N), Y(N), Z(N), the list.
!
!    Input, integer INDX(N), the sort index of the list.
!
!    Input, real ( kind = rk ) XVAL, YVAL, ZVAL, the value to be sought.
!
!    Output, integer LESS, EQUAL, MORE, the indexes in INDX of the
!    entries of X that are just less than, equal to, and just greater
!    than XVAL.  If XVAL does not occur in X, then EQUAL is zero.
!    If XVAL is the minimum entry of X, then LESS is 0.  If XVAL
!    is the greatest entry of X, then MORE is N+1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxn

  character c
  integer equal
  integer hi
  integer indx(maxn)
  integer less
  integer lo
  integer mid
  integer more
  integer n
  character r8vec3_compare
  real ( kind = rk ) x(maxn)
  real ( kind = rk ) xhi
  real ( kind = rk ) xlo
  real ( kind = rk ) xmid
  real ( kind = rk ) xval
  real ( kind = rk ) y(maxn)
  real ( kind = rk ) yhi
  real ( kind = rk ) ylo
  real ( kind = rk ) ymid
  real ( kind = rk ) yval
  real ( kind = rk ) z(maxn)
  real ( kind = rk ) zhi
  real ( kind = rk ) zlo
  real ( kind = rk ) zmid
  real ( kind = rk ) zval

  if ( n <= 0 ) then
    less = 0
    equal = 0
    more = 0
    return
  end if

  lo = 1
  hi = n

  xlo = x(indx(lo))
  ylo = y(indx(lo))
  zlo = z(indx(lo))

  xhi = x(indx(hi))
  yhi = y(indx(hi))
  zhi = z(indx(hi))

  c = r8vec3_compare ( xval, yval, zval, xlo, ylo, zlo )

  if ( c == '<' ) then
    less = 0
    equal = 0
    more = 1
    return
  else if ( c == '=' ) then
    less = 0
    equal = 1
    more = 2
    return
  end if

  c = r8vec3_compare ( xval, yval, zval, xhi, yhi, zhi )

  if ( c == '>' ) then
    less = n
    equal = 0
    more = n + 1
    return
  else if ( c == '=' ) then
    less = n - 1
    equal = n
    more = n + 1
    return
  end if

  do

    if ( lo + 1 == hi ) then
      less = lo
      equal = 0
      more = hi
      return
    end if

    mid = ( lo + hi ) / 2
    xmid = x(indx(mid))
    ymid = y(indx(mid))
    zmid = z(indx(mid))

    c = r8vec3_compare ( xval, yval, zval, xmid, ymid, zmid )

    if ( c == '=' ) then
      equal = mid
      less = equal - 1
      more = equal + 1
      return
    else if ( c == '<' ) then
      hi = mid
    else if ( c == '>' ) then
      lo = mid
    end if

  end do

  return
end
subroutine s_blanks_delete ( s )

!*****************************************************************************80
!
!! s_blanks_delete() replaces consecutive blanks by one blank.
!
!  Discussion:
!
!    The remaining characters are left justified and right padded with blanks.
!    TAB characters are converted to spaces.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
  implicit none

  integer i
  integer j
  character newchr
  character oldchr
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

  j = 0
  newchr = ' '

  do i = 1, len ( s )

    oldchr = newchr
    newchr = s(i:i)

    if ( newchr == TAB ) then
      newchr = ' '
    end if

    s(i:i) = ' '

    if ( oldchr /= ' ' .or. newchr /= ' ' ) then
      j = j + 1
      s(j:j) = newchr
    end if

  end do

  return
end
subroutine s_cat ( s1, s2, s3 )

!*****************************************************************************80
!
!! s_cat() concatenates two strings to make a third string.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S1, the "prefix" string.
!
!    Input, character ( len = * ) S2, the "postfix" string.
!
!    Output, character ( len = * ) S3, the string made by
!    concatenating S1 and S2, ignoring any trailing blanks.
!
  implicit none

  character ( len = * ) s1
  character ( len = * ) s2
  character ( len = * ) s3

  s3 = trim ( s1 ) // trim ( s2 )

  return
end
function s_eqi ( strng1, strng2 )

!*****************************************************************************80
!
!! s_eqi() is a case insensitive comparison of two strings for equality.
!
!  Example:
!
!    S_EQI ( 'Anjana', 'ANJANA' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRNG1, STRNG2, the strings to compare.
!
!    Output, logical S_EQI, the result of the comparison.
!
  implicit none

  integer i
  integer len1
  integer len2
  integer lenc
  logical s_eqi
  character s1
  character s2
  character ( len = * ) strng1
  character ( len = * ) strng2

  len1 = len ( strng1 )
  len2 = len ( strng2 )
  lenc = min ( len1, len2 )

  s_eqi = .false.

  do i = 1, lenc

    s1 = strng1(i:i)
    s2 = strng2(i:i)
    call ch_cap ( s1 )
    call ch_cap ( s2 )

    if ( s1 /= s2 ) then
      return
    end if

  end do

  do i = lenc + 1, len1
    if ( strng1(i:i) /= ' ' ) then
      return
    end if
  end do

  do i = lenc + 1, len2
    if ( strng2(i:i) /= ' ' ) then
      return
    end if
  end do

  s_eqi = .true.

  return
end
subroutine s_to_i4 ( string, ival, ierror, last )

!*****************************************************************************80
!
!! s_to_i4() reads an I4 from a string.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRING, a string to be examined.
!
!    Output, integer IVAL, the integer value read from the string.
!    If STRING is blank, then IVAL will be returned 0.
!
!    Output, integer IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer LAST, the last character in STRING that was
!    part of the representation of IVAL.
!
  implicit none

  character c
  integer i
  integer ierror
  integer isgn
  integer istate
  integer ival
  integer last
  integer lens
  character ( len = * ) string

  ierror = 0
  istate = 0

  isgn = 1
  ival = 0

  lens = len ( string )

  i = 0

  do

    i = i + 1

    c = string(i:i)

    if ( istate == 0 ) then

      if ( c == ' ' ) then

      else if ( c == '-' ) then
        istate = 1
        isgn = -1
      else if ( c == '+' ) then
        istate = 1
        isgn = + 1
      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if

    else if ( istate == 1 ) then

      if ( c == ' ' ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if

    else if ( istate == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        ival = 10 * ival + ichar ( c ) - ichar ( '0' )
      else
        istate = 3
      end if

    end if
!
!  Continue or exit?
!
    if ( istate == 3 ) then
      ival = isgn * ival
      last = i - 1
      return
    else if ( lens <= i ) then
      if ( istate == 2 ) then
        ival = isgn * ival
        last = lens
      else
        ierror = 1
        last = 0
      end if
      return
    end if

  end do

  return
end
subroutine s_to_r8 ( s, r, ierror, lchar )

!*****************************************************************************80
!
!! s_to_r8() reads an R8 from a string.
!
!  Discussion:
!
!    This routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the real number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 spaces
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semi4colon.
!
!    with most quantities optional.
!
!  Example:
!
!    S                 R
!
!    '1'               1.0
!    '     1   '       1.0
!    '1A'              1.0
!    '12,34,56'        12.0
!    '  34 7'          34.0
!    '-1E2ABCD'        -100.0
!    '-1X2ABCD'        -1.0
!    ' 2E-1'           0.2
!    '23.45'           23.45
!    '-4.2E+2'         -420.0
!    '17d2'            1700.0
!    '-14e-2'         -0.14
!    'e2'              100.0
!    '-12.73e-9.23'   -12.73 * 10.0^(-9.23)
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string containing the
!    data to be read.  Reading will begin at position 1 and
!    terminate at the end of the string, or when no more
!    characters can be read to form a legal real.  Blanks,
!    commas, or other nonnumeric data will, in particular,
!    cause the conversion to halt.
!
!    Output, real ( kind = rk ) R, the real value that was read from the string.
!
!    Output, integer IERROR, error flag.
!
!    0, no errors occurred.
!
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer LCHAR, the number of characters read from
!    the string to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  logical ch_eqi
  character c
  integer ierror
  integer ihave
  integer isgn
  integer iterm
  integer jbot
  integer jsgn
  integer jtop
  integer lchar
  integer nchar
  integer ndig
  real ( kind = rk ) r
  real ( kind = rk ) rbot
  real ( kind = rk ) rexp
  real ( kind = rk ) rtop
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

  nchar = len_trim ( s )
  ierror = 0
  r = 0.0D+00
  lchar = - 1
  isgn = 1
  rtop = 0.0D+00
  rbot = 1.0D+00
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

  do

    lchar = lchar + 1
    c = s(lchar+1:lchar+1)
!
!  Blank or TAB character.
!
    if ( c == ' ' .or. c == TAB ) then

      if ( ihave == 2 ) then

      else if ( ihave == 6 .or. ihave == 7 ) then
        iterm = 1
      else if ( 1 < ihave ) then
        ihave = 11
      end if
!
!  Comma.
!
    else if ( c == ',' .or. c == ';' ) then

      if ( ihave /= 1 ) then
        iterm = 1
        ihave = 12
        lchar = lchar + 1
      end if
!
!  Minus sign.
!
    else if ( c == '-' ) then

      if ( ihave == 1 ) then
        ihave = 2
        isgn = - 1
      else if ( ihave == 6 ) then
        ihave = 7
        jsgn = - 1
      else
        iterm = 1
      end if
!
!  Plus sign.
!
    else if ( c == '+' ) then

      if ( ihave == 1 ) then
        ihave = 2
      else if ( ihave == 6 ) then
        ihave = 7
      else
        iterm = 1
      end if
!
!  Decimal point.
!
    else if ( c == '.' ) then

      if ( ihave < 4 ) then
        ihave = 4
      else if ( 6 <= ihave .and. ihave <= 8 ) then
        ihave = 9
      else
        iterm = 1
      end if
!
!  Exponent marker.
!
    else if ( ch_eqi ( c, 'E' ) .or. ch_eqi ( c, 'D' ) ) then

      if ( ihave < 6 ) then
        ihave = 6
      else
        iterm = 1
      end if
!
!  Digit.
!
    else if ( ihave < 11 .and. lge ( c, '0' ) .and. lle ( c, '9' ) ) then

      if ( ihave <= 2 ) then
        ihave = 3
      else if ( ihave == 4 ) then
        ihave = 5
      else if ( ihave == 6 .or. ihave == 7 ) then
        ihave = 8
      else if ( ihave == 9 ) then
        ihave = 10
      end if

      call ch_to_digit ( c, ndig )

      if ( ihave == 3 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = rk )
      else if ( ihave == 5 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = rk )
        rbot = 10.0D+00 * rbot
      else if ( ihave == 8 ) then
        jtop = 10 * jtop + ndig
      else if ( ihave == 10 ) then
        jtop = 10 * jtop + ndig
        jbot = 10 * jbot
      end if
!
!  Anything else is regarded as a terminator.
!
    else
      iterm = 1
    end if
!
!  If we haven't seen a terminator, and we haven't examined the
!  entire string, go get the next character.
!
    if ( iterm == 1 .or. nchar <= lchar+1 ) then
      exit
    end if

  end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LCHAR is equal to NCHAR.
!
  if ( iterm /= 1 .and. lchar+1 == nchar ) then
    lchar = nchar
  end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
  if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then

    ierror = ihave

    return
  end if
!
!  Number seems OK.  Form it.
!
  if ( jtop == 0 ) then
    rexp = 1.0D+00
  else

    if ( jbot == 1 ) then
      rexp = 10.0D+00**( jsgn * jtop )
    else
      rexp = jsgn * jtop
      rexp = rexp / jbot
      rexp = 10.0D+00**rexp
    end if

  end if

  r = isgn * rexp * rtop / rbot

  return
end
subroutine shape_2d_edges_to_ps ( plotxmin2, plotymin2, alpha, iunit, &
  max_order, nface, nnode, face, face_order, x, y, xmin, ymin )

!*****************************************************************************80
!
!! shape_2d_edges_to_ps() writes 2D shape edges to a PostScript file.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer PLOTXMIN2, PLOTYMIN2, the Postscript origin.
!
!    Input, real ( kind = rk ) ALPHA, the physical-to-Postscript scale factor.
!
!    Input, integer IUNIT, the output FORTRAN unit.
!
!    Input, integer MAX_ORDER, the maximum number of nodes per
!    face.
!
!    Input, integer NFACE, the number of faces.
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, integer FACE(MAX_ORDER,NFACE), the nodes making faces.
!
!    Input, integer FACE_ORDER(NFACE), the number of nodes per
!    face.
!
!    Input, real ( kind = rk ) X(NNODE), Y(NNODE), the coordinates of points.
!
!    Input, real ( kind = rk ) XMIN, YMIN, the physical origin.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer max_order
  integer nface
  integer nnode

  real ( kind = rk ) alpha
  integer face(max_order,nface)
  integer face_order(nface)
  integer iface
  integer iunit
  integer j
  integer node
  integer plotxmin2
  integer plotymin2
  integer px
  integer py
  real ( kind = rk ) x(nnode)
  real ( kind = rk ) xmin
  real ( kind = rk ) y(nnode)
  real ( kind = rk ) ymin
!
!  Draw faces and fill them.
!
  do iface = 1, nface

    write ( iunit, '(a)' ) 'newpath'

    node = face(face_order(iface),iface)
    px = plotxmin2 + nint ( alpha * ( x(node) - xmin ) )
    py = plotymin2 + nint ( alpha * ( y(node) - ymin ) )
    write ( iunit, '(2i4,a,2i4,a)' ) px, py, ' moveto '

    do j = 1, face_order(iface)

      node = face(j,iface)
      px = plotxmin2 + nint ( alpha * ( x(node) - xmin ) )
      py = plotymin2 + nint ( alpha * ( y(node) - ymin ) )
      write ( iunit, '(2i4,a,2i4,a)' ) px, py, ' lineto '

    end do

    write ( iunit, '(a)' ) 'stroke'
!   write ( iunit, '(a)' ) 'fill'

  end do

  return
end
subroutine shape_2d_faces_to_ps ( plotxmin2, plotymin2, alpha, iunit, &
  max_order, nface, nnode, face, face_order, x, y, xmin, ymin )

!*****************************************************************************80
!
!! shape_2d_faces_to_ps() writes 2D shape faces to a PostScript file.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer PLOTXMIN2, PLOTYMIN2, the Postscript origin.
!
!    Input, real ( kind = rk ) ALPHA, the physical-to-Postscript scale factor.
!
!    Input, integer IUNIT, the output FORTRAN unit.
!
!    Input, integer MAX_ORDER, the maximum number of nodes per
!    face.
!
!    Input, integer NFACE, the number of faces.
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, integer FACE(MAX_ORDER,NFACE), the nodes making faces.
!
!    Input, integer FACE_ORDER(NFACE), the number of nodes per
!    face.
!
!    Input, real ( kind = rk ) X(NNODE), Y(NNODE), the coordinates of points.
!
!    Input, real ( kind = rk ) XMIN, YMIN, the physical origin.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer max_order
  integer nface
  integer nnode

  real ( kind = rk ) alpha
  real ( kind = rk ) blue
  integer face(max_order,nface)
  real ( kind = rk ) green
  integer i
  integer iface
  integer iunit
  integer j
  integer node
  integer face_order(nface)
  integer plotxmin2
  integer plotymin2
  integer px
  integer py
  real ( kind = rk ) red
  real ( kind = rk ) x(nnode)
  real ( kind = rk ) xmin
  real ( kind = rk ) y(nnode)
  real ( kind = rk ) ymin
!
!  Draw the faces.
!
  do iface = 1, nface

    do i = 1, 2

      if ( i == 1 ) then
        red = 0.9D+00
        green = 0.9D+00
        blue = 1.0D+00
      else
        red = 0.0D+00
        green = 0.0D+00
        blue = 0.0D+00
      end if

      write ( iunit, '(3f7.4,a)' ) red, green, blue, ' setrgbcolor'

      write ( iunit, '(a)' ) 'newpath'

      node = face(face_order(iface),iface)
      px = plotxmin2 + nint ( alpha * ( x(node) - xmin ) )
      py = plotymin2 + nint ( alpha * ( y(node) - ymin ) )
      write ( iunit, '(2i4,a,2i4,a)' ) px, py, ' moveto '

      do j = 1, face_order(iface)

        node = face(j,iface)
        px = plotxmin2 + nint ( alpha * ( x(node) - xmin ) )
        py = plotymin2 + nint ( alpha * ( y(node) - ymin ) )
        write ( iunit, '(2i4,a,2i4,a)' ) px, py, ' lineto '

      end do

      if ( i == 1 ) then
        write ( iunit, '(a)' ) 'fill'
      else
        write ( iunit, '(a)' ) 'stroke'
      end if

    end do

  end do

  return
end
subroutine shape_2d_nodes_to_ps ( plotxmin2, plotymin2, alpha, iunit, &
  nnode, x, y, xmin, ymin )

!*****************************************************************************80
!
!! shape_2d_nodes_to_ps() writes 2D shape nodes to a PostScript file.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer PLOTXMIN2, PLOTYMIN2, the Postscript origin.
!
!    Input, real ( kind = rk ) ALPHA, the physical-to-Postscript scale factor.
!
!    Input, integer IUNIT, the output FORTRAN unit.
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, real ( kind = rk ) X(NNODE), Y(NNODE), the X and Y components
!    of points.
!
!    Input, real ( kind = rk ) XMIN, YMIN, the physical origin.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nnode

  real ( kind = rk ) alpha
  integer i
  integer iunit
  integer plotxmin2
  integer plotymin2
  integer px
  integer py
  real ( kind = rk ) x(nnode)
  real ( kind = rk ) xmin
  real ( kind = rk ) y(nnode)
  real ( kind = rk ) ymin
!
!  Draw the nodes.
!
  do i = 1, nnode

    px = plotxmin2 + nint ( alpha * ( x(i) - xmin ) )
    py = plotymin2 + nint ( alpha * ( y(i) - ymin ) )

    write ( iunit, '(a,2i4,a)' ) 'newpath ', px, py, &
      ' 5 0 360 arc closepath stroke'

  end do

  return
end
subroutine shape_3d_edges_to_ps ( file_name, max_order, nface, nnode, &
  face, face_order, x, y, z )

!*****************************************************************************80
!
!! shape_3d_edges_to_ps() writes 3D shape edges to a PostScript file.
!
!  Discussion:
!
!    Four views are created in one picture: XY, YZ, ZX, and XYZ.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the output file.
!
!    Input, integer MAX_ORDER, the maximum number of nodes per
!    face.
!
!    Input, integer NFACE, the number of faces.
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, integer FACE(MAX_ORDER,NFACE), the nodes making faces.
!
!    Input, integer FACE_ORDER(NFACE), the number of nodes per
!    face.
!
!    Input, real ( kind = rk ) X(NNODE), Y(NNODE), Z(NNODE), the X, Y and Z
!    components of points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer max_order
  integer nface
  integer nnode

  real ( kind = rk ) alpha
  real ( kind = rk ) blue
  character ( len = 8 ) date
  integer face(max_order,nface)
  integer face_order(nface)
  character ( len = * ) file_name
  real ( kind = rk ) green
  integer ios
  integer iunit
  integer, parameter :: margin = 36
  integer pagexmax
  integer pagexmin
  integer pageymax
  integer pageymin
  integer plotxmax
  integer plotxmin
  integer plotxmin2
  integer plotymax
  integer plotymin
  integer plotymin2
  integer px1
  integer px2
  integer px3
  integer px4
  integer px5
  integer py1
  integer py2
  integer py3
  integer py4
  integer py5
  real ( kind = rk ) red
  real ( kind = rk ) x(nnode)
  real ( kind = rk ) xmax
  real ( kind = rk ) xmin
  real ( kind = rk ) xx(nnode)
  real ( kind = rk ) y(nnode)
  real ( kind = rk ) ymax
  real ( kind = rk ) ymin
  real ( kind = rk ) yy(nnode)
  real ( kind = rk ) z(nnode)
!
!  Open the file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = file_name, status = 'replace', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SHAPE_3D_EDGES_TO_PS(): Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    return
  end if
!
!  Write the prolog.
!
  pagexmax = 612
  pagexmin = 0
  pageymax = 792
  pageymin = 0

  px1 = 0
  px2 = margin
  px3 = pagexmax / 2
  px4 = pagexmax - margin
  px5 = pagexmax

  py1 = 0
  py2 = margin
  py3 = pageymax / 2
  py4 = pageymax - margin
  py5 = pageymax

  write ( iunit, '(a)' ) '%!PS-Adobe-3.0'
  write ( iunit, '(a)' ) '%%Document-Fonts: Times-Roman'
  write ( iunit, '(a,a)' ) '%%Title: ' , trim ( file_name )
  write ( iunit, '(a)' ) '%%Creator: GRAFPACK(shape_3d_edges_to_ps)'
  call date_and_time ( date )
  write ( iunit, '(a)' ) '%%CreationDate: ' // trim ( date )
  write ( iunit, '(a)' ) '%%BoundingBox 0 0 612 794'
  write ( iunit, '(a)' ) '%%LanguageLevel: 2'
  write ( iunit, '(a)' ) '%%EndComments'
  write ( iunit, '(a)' ) '%%BeginProlog'
  write ( iunit, '(a)' ) '%%EndProlog'
!
!  Draw gray lines to separate the boxes.
!
  red = 0.5
  green = 0.5
  blue = 0.5

  write ( iunit, '(3f7.4,a)' ) red, green, blue, ' setrgbcolor'

  write ( iunit, '(2i4,a)' ) px2, py2, ' moveto '
  write ( iunit, '(2i4,a)' ) px4, py2, ' lineto '
  write ( iunit, '(2i4,a)' ) px2, py3, ' moveto '
  write ( iunit, '(2i4,a)' ) px4, py3, ' lineto '
  write ( iunit, '(2i4,a)' ) px2, py4, ' moveto '
  write ( iunit, '(2i4,a)' ) px4, py4, ' lineto '

  write ( iunit, '(2i4,a)' ) px2, py2, ' moveto '
  write ( iunit, '(2i4,a)' ) px2, py4, ' lineto '
  write ( iunit, '(2i4,a)' ) px3, py2, ' moveto '
  write ( iunit, '(2i4,a)' ) px3, py4, ' lineto '
  write ( iunit, '(2i4,a)' ) px4, py2, ' moveto '
  write ( iunit, '(2i4,a)' ) px4, py4, ' lineto '
  write ( iunit, '(a)' ) 'stroke'
!
!  Determine ALPHA, the single scale factor to be used for both
!  directions, and all four plots!
!
  xx(1:nnode) = x(1:nnode)
  yy(1:nnode) = y(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  alpha = min ( ( px3 - px2 ) / ( xmax - xmin ), &
                ( py4 - py3 ) / ( ymax - ymin ) )

  xx(1:nnode) = y(1:nnode)
  yy(1:nnode) = z(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  alpha = min ( alpha, &
                ( px4 - px3 ) / ( xmax - xmin ), &
                ( py4 - py3 ) / ( ymax - ymin ) )

  xx(1:nnode) = z(1:nnode)
  yy(1:nnode) = x(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  alpha = min ( alpha, &
                ( px3 - px2 ) / ( xmax - xmin ), &
                ( py3 - py2 ) / ( ymax - ymin ) )

  xx(1:nnode) = 0.80 * x(1:nnode) - 0.31 * y(1:nnode) + 0.50 * z(1:nnode)
  yy(1:nnode) = 0.50 * x(1:nnode) + 0.80 * y(1:nnode) - 0.31 * z(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  alpha = min ( alpha, &
                ( px4 - px3 ) / ( xmax - xmin ), &
                ( py3 - py2 ) / ( ymax - ymin ) )
!
!  Set the fill color.
!
  red = 0.9D+00
  green = 0.9D+00
  blue = 1.0D+00

  write ( iunit, '(3f7.4,a)' ) red, green, blue, ' setrgbcolor'
!
!  XY edges.
!
  plotxmin = px2
  plotxmax = px3
  plotymin = py3
  plotymax = py4

  xx(1:nnode) = x(1:nnode)
  yy(1:nnode) = y(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  plotxmin2 = int ( 0.5 * ( plotxmin + plotxmax - alpha * ( xmax - xmin ) ) )
  plotymin2 = int ( 0.5 * ( plotymin + plotymax - alpha * ( ymax - ymin ) ) )

  call shape_2d_edges_to_ps ( plotxmin2, plotymin2, alpha, iunit, &
    max_order, nface, nnode, face, face_order, xx, yy, xmin, ymin )
!
!  YZ edges.
!
  plotxmin = px3
  plotxmax = px4
  plotymin = py3
  plotymax = py4

  xx(1:nnode) = y(1:nnode)
  yy(1:nnode) = z(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  plotxmin2 = int ( 0.5 * ( plotxmin + plotxmax - alpha * ( xmax - xmin ) ) )
  plotymin2 = int ( 0.5 * ( plotymin + plotymax - alpha * ( ymax - ymin ) ) )

  call shape_2d_edges_to_ps ( plotxmin2, plotymin2, alpha, iunit, &
    max_order, nface, nnode, face, face_order, xx, yy, xmin, ymin )
!
!  ZX edges.
!
  plotxmin = px2
  plotxmax = px3
  plotymin = py2
  plotymax = py3

  xx(1:nnode) = z(1:nnode)
  yy(1:nnode) = x(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  plotxmin2 = int ( 0.5 * ( plotxmin + plotxmax - alpha * ( xmax - xmin ) ) )
  plotymin2 = int ( 0.5 * ( plotymin + plotymax - alpha * ( ymax - ymin ) ) )

  call shape_2d_edges_to_ps ( plotxmin2, plotymin2, alpha, iunit, &
    max_order, nface, nnode, face, face_order, xx, yy, xmin, ymin )
!
!  XYZ edges.
!
  plotxmin = px3
  plotxmax = px4
  plotymin = py2
  plotymax = py3

  xx(1:nnode) = 0.80 * x(1:nnode) - 0.31 * y(1:nnode) + 0.50 * z(1:nnode)
  yy(1:nnode) = 0.50 * x(1:nnode) + 0.80 * y(1:nnode) - 0.31 * z(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  plotxmin2 = int ( 0.5 * ( plotxmin + plotxmax - alpha * ( xmax - xmin ) ) )
  plotymin2 = int ( 0.5 * ( plotymin + plotymax - alpha * ( ymax - ymin ) ) )

  call shape_2d_edges_to_ps ( plotxmin2, plotymin2, alpha, iunit, &
    max_order, nface, nnode, face, face_order, xx, yy, xmin, ymin )

  write ( iunit, '(a)' ) 'showpage'
!
!  Write the epilog.
!
  write ( iunit, '(a)' ) 'grestore'
  write ( iunit, '(a)' ) '%%Trailer'
  write ( iunit, '(a,i2)' ) '%%Pages: 1'

  close ( unit = iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SHAPE_3D_EDGES_TO_PS'
  write ( *, '(a)' ) '  The data was written to the file: ' &
    // trim ( file_name )

  return
end
subroutine shape_3d_faces_to_ps ( file_name, max_order, nface, nnode, face, &
  face_order, x, y, z )

!*****************************************************************************80
!
!! shape_3d_faces_to_ps() writes 3D shape faces to a PostScript file.
!
!  Discussion:
!
!    Four views are created in one picture: XY, YZ, ZX, and XYZ.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the output file.
!
!    Input, integer MAX_ORDER, the maximum number of nodes per
!    face.
!
!    Input, integer NFACE, the number of faces.
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, integer FACE(MAX_ORDER,NFACE), the nodes making faces.
!
!    Input, integer FACE_ORDER(NFACE), the number of nodes per
!    face.
!
!    Input, real ( kind = rk ) X(NNODE), Y(NNODE), Z(NNODE), the X, Y and Z
!    components of points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer max_order
  integer nface
  integer nnode

  real ( kind = rk ) alpha
  real ( kind = rk ) blue
  character ( len = 8 ) date
  integer face(max_order,nface)
  integer face_order(nface)
  character ( len = * ) file_name
  real ( kind = rk ) green
  integer ios
  integer iunit
  integer, parameter :: margin = 36
  integer pagexmax
  integer pagexmin
  integer pageymax
  integer pageymin
  integer plotxmax
  integer plotxmin
  integer plotxmin2
  integer plotymax
  integer plotymin
  integer plotymin2
  integer px1
  integer px2
  integer px3
  integer px4
  integer px5
  integer py1
  integer py2
  integer py3
  integer py4
  integer py5
  real ( kind = rk ) red
  real ( kind = rk ) x(nnode)
  real ( kind = rk ) xmax
  real ( kind = rk ) xmin
  real ( kind = rk ) xx(nnode)
  real ( kind = rk ) y(nnode)
  real ( kind = rk ) ymax
  real ( kind = rk ) ymin
  real ( kind = rk ) yy(nnode)
  real ( kind = rk ) z(nnode)
!
!  Open the file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = file_name, status = 'replace', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SHAPE_3D_EDGES_TO_PS(): Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    return
  end if
!
!  Write the prolog.
!
  pagexmax = 612
  pagexmin = 0
  pageymax = 792
  pageymin = 0

  px1 = 0
  px2 = margin
  px3 = pagexmax / 2
  px4 = pagexmax - margin
  px5 = pagexmax

  py1 = 0
  py2 = margin
  py3 = pageymax / 2
  py4 = pageymax - margin
  py5 = pageymax

  write ( iunit, '(a)' ) '%!PS-Adobe-3.0'
  write ( iunit, '(a)' ) '%%Document-Fonts: Times-Roman'
  write ( iunit, '(a,a)' ) '%%Title: ' , trim ( file_name )
  write ( iunit, '(a)' ) '%%Creator: GRAFPACK(shape_3d_edges_to_ps)'
  call date_and_time ( date )
  write ( iunit, '(a)' ) '%%CreationDate: ' // trim ( date )
  write ( iunit, '(a)' ) '%%BoundingBox 0 0 612 794'
  write ( iunit, '(a)' ) '%%LanguageLevel: 2'
  write ( iunit, '(a)' ) '%%EndComments'
  write ( iunit, '(a)' ) '%%BeginProlog'
  write ( iunit, '(a)' ) '%%EndProlog'
!
!  Draw gray lines to separate the boxes.
!
  red = 0.5D+00
  green = 0.5D+00
  blue = 0.5D+00

  write ( iunit, '(3f7.4,a)' ) red, green, blue, ' setrgbcolor'

  write ( iunit, '(2i4,a)' ) px2, py2, ' moveto '
  write ( iunit, '(2i4,a)' ) px4, py2, ' lineto '
  write ( iunit, '(2i4,a)' ) px2, py3, ' moveto '
  write ( iunit, '(2i4,a)' ) px4, py3, ' lineto '
  write ( iunit, '(2i4,a)' ) px2, py4, ' moveto '
  write ( iunit, '(2i4,a)' ) px4, py4, ' lineto '

  write ( iunit, '(2i4,a)' ) px2, py2, ' moveto '
  write ( iunit, '(2i4,a)' ) px2, py4, ' lineto '
  write ( iunit, '(2i4,a)' ) px3, py2, ' moveto '
  write ( iunit, '(2i4,a)' ) px3, py4, ' lineto '
  write ( iunit, '(2i4,a)' ) px4, py2, ' moveto '
  write ( iunit, '(2i4,a)' ) px4, py4, ' lineto '
  write ( iunit, '(a)' ) 'stroke'
!
!  Determine ALPHA, the single scale factor to be used for both
!  directions, and all four plots!
!
  xx(1:nnode) = x(1:nnode)
  yy(1:nnode) = y(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  alpha = min ( ( px3 - px2 ) / ( xmax - xmin ), &
                ( py4 - py3 ) / ( ymax - ymin ) )

  xx(1:nnode) = y(1:nnode)
  yy(1:nnode) = z(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  alpha = min ( alpha, &
                ( px4 - px3 ) / ( xmax - xmin ), &
                ( py4 - py3 ) / ( ymax - ymin ) )

  xx(1:nnode) = z(1:nnode)
  yy(1:nnode) = x(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  alpha = min ( alpha, &
                ( px3 - px2 ) / ( xmax - xmin ), &
                ( py3 - py2 ) / ( ymax - ymin ) )

  xx(1:nnode) = 0.80 * x(1:nnode) - 0.31 * y(1:nnode) + 0.50 * z(1:nnode)
  yy(1:nnode) = 0.50 * x(1:nnode) + 0.80 * y(1:nnode) - 0.31 * z(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  alpha = min ( alpha, &
                ( px4 - px3 ) / ( xmax - xmin ), &
                ( py3 - py2 ) / ( ymax - ymin ) )
!
!  Set the fill color.
!
  red = 0.9D+00
  green = 0.9D+00
  blue = 1.0D+00

  write ( iunit, '(3f7.4,a)' ) red, green, blue, ' setrgbcolor'
!
!  XY edges.
!
  plotxmin = px2
  plotxmax = px3
  plotymin = py3
  plotymax = py4

  xx(1:nnode) = x(1:nnode)
  yy(1:nnode) = y(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  plotxmin2 = int ( 0.5 * ( plotxmin + plotxmax - alpha * ( xmax - xmin ) ) )
  plotymin2 = int ( 0.5 * ( plotymin + plotymax - alpha * ( ymax - ymin ) ) )

  call shape_2d_faces_to_ps ( plotxmin2, plotymin2, alpha, iunit, &
    max_order, nface, nnode, face, face_order, xx, yy, xmin, ymin )
!
!  YZ edges.
!
  plotxmin = px3
  plotxmax = px4
  plotymin = py3
  plotymax = py4

  xx(1:nnode) = y(1:nnode)
  yy(1:nnode) = z(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  plotxmin2 = int ( 0.5 * ( plotxmin + plotxmax - alpha * ( xmax - xmin ) ) )
  plotymin2 = int ( 0.5 * ( plotymin + plotymax - alpha * ( ymax - ymin ) ) )

  call shape_2d_faces_to_ps ( plotxmin2, plotymin2, alpha, iunit, &
    max_order, nface, nnode, face, face_order, xx, yy, xmin, ymin )
!
!  ZX edges.
!
  plotxmin = px2
  plotxmax = px3
  plotymin = py2
  plotymax = py3

  xx(1:nnode) = z(1:nnode)
  yy(1:nnode) = x(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  plotxmin2 = int ( 0.5 * ( plotxmin + plotxmax - alpha * ( xmax - xmin ) ) )
  plotymin2 = int ( 0.5 * ( plotymin + plotymax - alpha * ( ymax - ymin ) ) )

  call shape_2d_faces_to_ps ( plotxmin2, plotymin2, alpha, iunit, &
    max_order, nface, nnode, face, face_order, xx, yy, xmin, ymin )
!
!  XYZ edges.
!
  plotxmin = px3
  plotxmax = px4
  plotymin = py2
  plotymax = py3

  xx(1:nnode) = 0.80 * x(1:nnode) - 0.31 * y(1:nnode) + 0.50 * z(1:nnode)
  yy(1:nnode) = 0.50 * x(1:nnode) + 0.80 * y(1:nnode) - 0.31 * z(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  plotxmin2 = int ( 0.5 * ( plotxmin + plotxmax - alpha * ( xmax - xmin ) ) )
  plotymin2 = int ( 0.5 * ( plotymin + plotymax - alpha * ( ymax - ymin ) ) )

  call shape_2d_faces_to_ps ( plotxmin2, plotymin2, alpha, iunit, &
    max_order, nface, nnode, face, face_order, xx, yy, xmin, ymin )

  write ( iunit, '(a)' ) 'showpage'
!
!  Write the epilog.
!
  write ( iunit, '(a)' ) 'grestore'
  write ( iunit, '(a)' ) '%%Trailer'
  write ( iunit, '(a,i2)' ) '%%Pages: 1'

  close ( unit = iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SHAPE_3D_EDGES_TO_PS'
  write ( *, '(a)' ) '  The data was written to the file: ' &
    // trim ( file_name )

  return
end
subroutine shape_3d_nodes_to_ps ( file_name, nnode, x, y, z )

!*****************************************************************************80
!
!! shape_3d_nodes_to_ps() writes 3D shape nodes to a PostScript file.
!
!  Discussion:
!
!    Four views are created in one picture: XY, YZ, ZX, and XYZ.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the output file.
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, real ( kind = rk ) X(NNODE), Y(NNODE), Z(NNODE), the X, Y and Z
!    components of points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nnode

  real ( kind = rk ) alpha
  real ( kind = rk ) blue
  character ( len = 8 ) date
  character ( len = * ) file_name
  real ( kind = rk ) green
  integer ios
  integer iunit
  integer, parameter :: margin = 36
  integer pagexmax
  integer pagexmin
  integer pageymax
  integer pageymin
  integer plotxmax
  integer plotxmin
  integer plotxmin2
  integer plotymax
  integer plotymin
  integer plotymin2
  integer px1
  integer px2
  integer px3
  integer px4
  integer px5
  integer py1
  integer py2
  integer py3
  integer py4
  integer py5
  real ( kind = rk ) red
  real ( kind = rk ) x(nnode)
  real ( kind = rk ) xmax
  real ( kind = rk ) xmin
  real ( kind = rk ) xx(nnode)
  real ( kind = rk ) y(nnode)
  real ( kind = rk ) ymax
  real ( kind = rk ) ymin
  real ( kind = rk ) yy(nnode)
  real ( kind = rk ) z(nnode)
!
!  Open the file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = file_name, status = 'replace', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SHAPE_3D_NODES_TO_PS(): Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    return
  end if
!
!  Write the prolog.
!
  pagexmax = 612
  pagexmin = 0
  pageymax = 792
  pageymin = 0

  px1 = 0
  px2 = margin
  px3 = pagexmax / 2
  px4 = pagexmax - margin
  px5 = pagexmax

  py1 = 0
  py2 = margin
  py3 = pageymax / 2
  py4 = pageymax - margin
  py5 = pageymax

  write ( iunit, '(a)' ) '%!PS-Adobe-3.0'
  write ( iunit, '(a)' ) '%%Document-Fonts: Times-Roman'
  write ( iunit, '(a,a)' ) '%%Title: ' , trim ( file_name )
  write ( iunit, '(a)' ) '%%Creator: GRAFPACK(shape_3d_nodes_to_ps)'
  call date_and_time ( date )
  write ( iunit, '(a)' ) '%%CreationDate: ' // trim ( date )
  write ( iunit, '(a)' ) '%%BoundingBox 0 0 612 794'
  write ( iunit, '(a)' ) '%%LanguageLevel: 2'
  write ( iunit, '(a)' ) '%%EndComments'
  write ( iunit, '(a)' ) '%%BeginProlog'
  write ( iunit, '(a)' ) '%%EndProlog'
!
!  Draw gray lines to separate the boxes.
!
  red = 0.5D+00
  green = 0.5D+00
  blue = 0.5D+00

  write ( iunit, '(3f7.4,a)' ) red, green, blue, ' setrgbcolor'

  write ( iunit, '(2i4,a)' ) px2, py2, ' moveto '
  write ( iunit, '(2i4,a)' ) px4, py2, ' lineto '
  write ( iunit, '(2i4,a)' ) px2, py3, ' moveto '
  write ( iunit, '(2i4,a)' ) px4, py3, ' lineto '
  write ( iunit, '(2i4,a)' ) px2, py4, ' moveto '
  write ( iunit, '(2i4,a)' ) px4, py4, ' lineto '

  write ( iunit, '(2i4,a)' ) px2, py2, ' moveto '
  write ( iunit, '(2i4,a)' ) px2, py4, ' lineto '
  write ( iunit, '(2i4,a)' ) px3, py2, ' moveto '
  write ( iunit, '(2i4,a)' ) px3, py4, ' lineto '
  write ( iunit, '(2i4,a)' ) px4, py2, ' moveto '
  write ( iunit, '(2i4,a)' ) px4, py4, ' lineto '
  write ( iunit, '(a)' ) 'stroke'
!
!  Determine ALPHA, the single scale factor to be used for both
!  directions, and all four plots!
!
  xx(1:nnode) = x(1:nnode)
  yy(1:nnode) = y(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  alpha = min ( ( px3 - px2 ) / ( xmax - xmin ), &
                ( py4 - py3 ) / ( ymax - ymin ) )

  xx(1:nnode) = y(1:nnode)
  yy(1:nnode) = z(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  alpha = min ( alpha, &
                ( px4 - px3 ) / ( xmax - xmin ), &
                ( py4 - py3 ) / ( ymax - ymin ) )

  xx(1:nnode) = z(1:nnode)
  yy(1:nnode) = x(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  alpha = min ( alpha, &
                ( px3 - px2 ) / ( xmax - xmin ), &
                ( py3 - py2 ) / ( ymax - ymin ) )

  xx(1:nnode) = 0.80 * x(1:nnode) - 0.31 * y(1:nnode) + 0.50 * z(1:nnode)
  yy(1:nnode) = 0.50 * x(1:nnode) + 0.80 * y(1:nnode) - 0.31 * z(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  alpha = min ( alpha, &
                ( px4 - px3 ) / ( xmax - xmin ), &
                ( py3 - py2 ) / ( ymax - ymin ) )
!
!  Set the color.
!
  red = 0.3
  green = 0.3
  blue = 0.3

  write ( iunit, '(3f7.4,a)' ) red, green, blue, ' setrgbcolor'
!
!  XY nodes.
!
  plotxmin = px2
  plotxmax = px3
  plotymin = py3
  plotymax = py4

  xx(1:nnode) = x(1:nnode)
  yy(1:nnode) = y(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  plotxmin2 = int ( 0.5 * ( plotxmin + plotxmax - alpha * ( xmax - xmin ) ) )
  plotymin2 = int ( 0.5 * ( plotymin + plotymax - alpha * ( ymax - ymin ) ) )

  call shape_2d_nodes_to_ps ( plotxmin2, plotymin2, alpha, iunit, &
    nnode, xx, yy, xmin, ymin )
!
!  YZ edges.
!
  plotxmin = px3
  plotxmax = px4
  plotymin = py3
  plotymax = py4

  xx(1:nnode) = y(1:nnode)
  yy(1:nnode) = z(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  plotxmin2 = int ( 0.5 * ( plotxmin + plotxmax - alpha * ( xmax - xmin ) ) )
  plotymin2 = int ( 0.5 * ( plotymin + plotymax - alpha * ( ymax - ymin ) ) )

  call shape_2d_nodes_to_ps ( plotxmin2, plotymin2, alpha, iunit, &
    nnode, xx, yy, xmin, ymin )
!
!  ZX edges.
!
  plotxmin = px2
  plotxmax = px3
  plotymin = py2
  plotymax = py3

  xx(1:nnode) = z(1:nnode)
  yy(1:nnode) = x(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  plotxmin2 = int ( 0.5 * ( plotxmin + plotxmax - alpha * ( xmax - xmin ) ) )
  plotymin2 = int ( 0.5 * ( plotymin + plotymax - alpha * ( ymax - ymin ) ) )

  call shape_2d_nodes_to_ps ( plotxmin2, plotymin2, alpha, iunit, &
    nnode, xx, yy, xmin, ymin )
!
!  XYZ edges.
!
  plotxmin = px3
  plotxmax = px4
  plotymin = py2
  plotymax = py3

  xx(1:nnode) = 0.80 * x(1:nnode) - 0.31 * y(1:nnode) + 0.50 * z(1:nnode)
  yy(1:nnode) = 0.50 * x(1:nnode) + 0.80 * y(1:nnode) - 0.31 * z(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  plotxmin2 = int ( 0.5 * ( plotxmin + plotxmax - alpha * ( xmax - xmin ) ) )
  plotymin2 = int ( 0.5 * ( plotymin + plotymax - alpha * ( ymax - ymin ) ) )

  call shape_2d_nodes_to_ps ( plotxmin2, plotymin2, alpha, iunit, &
    nnode, xx, yy, xmin, ymin )

  write ( iunit, '(a)' ) 'showpage'
!
!  Write the epilog.
!
  write ( iunit, '(a)' ) 'grestore'
  write ( iunit, '(a)' ) '%%Trailer'
  write ( iunit, '(a,i2)' ) '%%Pages: 1'

  close ( unit = iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SHAPE_3D_NODES_TO_PS'
  write ( *, '(a)' ) '  The data was written to the file: ' &
    // trim ( file_name )

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
subroutine span_forest ( nnode, nedge, iendpt, k, component )

!*****************************************************************************80
!
!! span_forest() determines a graph's connectivity and spanning forest.
!
!  Discussion:
!
!    The input graph may be connected or unconnected.
!
!    If the input graph is connected, this routine simply returns a
!    spanning tree for the graph.
!
!    Definition: A (connected) component of a graph is a maximal subgraph
!    which is connected.
!
!    Definition: A tree is a connected graph containing no cycles.
!
!    Definition: A spanning tree of a connected graph is a subgraph which 
!    is a maximal tree.
!
!    Definition: A forest is a collection of trees, no two of which share 
!    a node.
!
!    Definition: A spanning forest of a possibly unconnected graph 
!    is a collection containing a single spanning tree for each component 
!    of the graph.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 May 1999
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
!    Input, integer NEDGE, the number of edges in graph.
!
!    Input/output, integer IENDPT(2,NEDGE), the edge array of 
!    the graph.  IENDPT(1,I) and IENDPT(2,I) are the two nodes that make up 
!    edge I.
!  
!    On input, IENDPT describes the graph.
!
!    On output, the input entries of IENDPT have been reordered, so that
!    edges belonging to the spanning forest come first, followed by those
!    edges which are not part of the spanning forest.
!
!    Output, integer K, the number of connected components of the 
!    graph.
!
!    Output, integer COMPONENT(NNODE), the component to which each 
!    node belongs.
!
  implicit none

  integer nedge
  integer nnode

  integer component(nnode)
  integer i
  integer iendpt(2,nedge)
  integer ip
  integer iq
  integer ir
  integer iret
  integer is
  integer itemp
  integer k
  integer l
  integer l0
  integer l1
  integer loc
  integer m
  integer m0
  integer m1
  integer mm

  mm = 1 + max ( nnode, nedge )
 
  do i = 1, nnode
    component(i) = -i
  end do
 
  do m = 1, nedge
    do l = 1, 2
      ip = iendpt(l,m)
      iendpt(l,m) = component(ip)
      component(ip) = - l * mm - m
    end do
  end do
 
  k = 0
  loc = 0
 
10 continue
 
  do i = 1, nnode
 
    iq = component(i)
 
    if ( iq <= 0 ) then
 
      k = k + 1
      component(i) = k
 
      if ( iq + i < 0 ) then
        ip = i
        is = - iq
        iret = 31
        l1 = - iq / mm
        m1 = - iq - l1 * mm
        go to 110
      end if
 
    end if
 
  end do
 
  do m = 1, nedge
 
    do
 
      ir = - iendpt(1,m)
 
      if ( ir < 0 ) then
        exit
      end if

      itemp = iendpt(2,m)
      iendpt(2,m) = iendpt(2,ir)
      iendpt(2,ir) = itemp

      iendpt(1,m) = iendpt(1,ir)
      iendpt(1,ir) = component(iendpt(2,ir))

    end do
 
  end do
 
  component(iendpt(2,1:loc)) = component(iendpt(1,1:loc))
 
  return
 
90    continue
 
  if ( ir /= 0 ) then
    loc = loc + 1
    component(ip) = iendpt(1,ir) + iendpt(2,ir) - ip
    iendpt(1,ir) = - loc
    iendpt(2,ir) = ip
  end if
 
  ip = m

  if ( m <= 0 ) then
    go to 10
  end if

  is = - component(ip)
 
100   continue
 
  l = is / mm
  m = is - l * mm

  if ( l == 0 ) then
    go to 90
  end if

  l1 = 3 - l
  m1 = m
 
110   continue
 
  iq = iendpt(l1,m1)
 
  if ( 0 < iq ) then
    if ( iq <= mm ) then
      if ( 0 <= component(iq) ) then
        ir = m1
      end if
    end if
  end if
 
  if ( 0 <= iq ) then
    is = abs ( iendpt(l,m) )
    iendpt(l,m) = ip
    go to 100
  end if
 
  if ( -mm <= iq ) then
 
    iq = - iq
    iendpt(l1,m1) = 0
 
    if ( iret == 31 ) then

      l0 = l1
      m0 = m1
      ir = 0
      iret = 43

    else
 
      iendpt(l0,m0) = iq
      l0 = l1
      m0 = m1
      is = abs ( iendpt(l,m) )
      iendpt(l,m) = ip

    end if

    go to 100
 
  end if
 
  iendpt(l1,m1) = - iq
  l1 = - iq / mm
  m1 = - iq - l1 * mm
  go to 110
 
end
subroutine span_tree_cand ( nedge, nnode, iarray, k, nstack, stack, &
  maxstack, iendpt, ierror, ncan )

!*****************************************************************************80
!
!! span_tree_cand() finds candidates for the K-th edge of a spanning tree.  
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 September 2000
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
!    Input, integer NEDGE, the number of edges.
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, integer IARRAY(NNODE).  IARRAY(I) is the I-th edge of
!    the spanning tree.
!
!    Input, integer K, index of position in IARRAY for which
!    candidates are needed.
!
!    Output, integer NSTACK, the current size of the stack.
!
!    Output, integer STACK(MAXSTACK).  List of candidates for all
!    positions.
!
!    Input, integer MAXSTACK, the maximum size of the stack.
!
!    Input, integer IENDPT(2,NEDGE).  IENDPT(1,I), IENDPT(2,I) are
!    the two nodes of edge I in graph.
!
!    Output, integer IERROR, error flag.  0 if no errors, or 1
!    if needed stack size reached available stacksize MAXSTACK.
!    You should increase the dimension of STACK and call again.
!
!    Input/output, integer NCAN(NNODE-1), the number of candidates
!    for each position.
!
  implicit none

  integer nedge
  integer nnode
  integer maxstack

  integer i
  integer iarray(nnode)
  integer iend(2,nnode)
  integer iendpt(2,nedge)
  integer ierror
  integer iwork(nnode)
  integer k
  integer ncan(nnode-1)
  integer ncomp
  integer nstack
  integer stack(maxstack)

  if ( k <= 0 ) then
    ierror = 1
    return
  end if

  ncan(k) = 0

  ierror = 0
 
  if ( k == 1 ) then
 
    nstack = nedge - nnode
 
    if ( maxstack < nstack ) then
      ierror = 1
      return
    end if
 
    call i4vec_indicator ( nstack, stack )
 
    ncan(k) = nedge - nnode
 
  else

    iend(1,1:k-1) = iendpt(1,iarray(1:k-1))
    iend(2,1:k-1) = iendpt(2,iarray(1:k-1))
 
    call span_forest ( nnode, k-1, iend, ncomp, iwork )
 
    do i = iarray(k-1)+1, nedge+k+1-nnode
 
      if ( iwork(iendpt(1,i)) /= iwork(iendpt(2,i)) ) then
 
        nstack = nstack + 1
          
        if ( maxstack < nstack ) then
          ierror = 1
          return
        end if
 
        stack(nstack) = i
        ncan(k) = ncan(k) + 1

      end if
 
    end do
 
  end if
 
  return
end
subroutine span_tree_next ( signal, nnode, nedge, iendpt, iarray, ncan )

!*****************************************************************************80
!
!! span_tree_next() uses backtracking to find spanning forests of a graph.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer SIGNAL.
!    On input, 0 means this is the first call for a new problem.
!    On output, 0 means no more solutions exist; 1 means another solution was fo
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, integer NEDGE, the number of edges.
!
!    Input, integer IENDPT(2,NEDGE), the edge array of the graph.
!
!    Output, integer IARRAY(NNODE-1).  If SIGNAL = 1, then IARRAY
!    contains the "next" spanning forest found by the routine, stored as a list
!    of edge indices.
!
!    Workspace, integer NCAN(NNODE-1), the number of candidates for each 
!    position.
!
  implicit none

  integer, parameter :: maxstack = 1000

  integer nedge
  integer nnode

  integer iarray(nnode-1)
  integer iendpt(2,nedge)
  integer ierror
  integer, save :: indx
  integer, save :: k
  integer, dimension ( nnode-1) :: ncan
  integer, save :: nstack
  integer signal
  integer, save, dimension ( maxstack ) :: stack
!
!  First call for this problem.
!
  if ( signal == 0 ) then

    iarray(1:nnode-1) = 0
    indx = 0
    k = 0
    ncan(1:nnode-1) = 0
    nstack = 0
    stack(1:maxstack) = 0

  end if
!
!  Try to extend the current partial solution.
!
  do

    call i4vec_backtrack ( nnode-1, iarray, indx, k, nstack, stack, &
      maxstack, ncan )
!
!  A full solution was found.
!
    if ( indx == 1 ) then

      signal = 1
      exit
!
!  A partial solution was found.  Seek candidates for the next entry.
!
    else if ( indx == 2 ) then

      call span_tree_cand ( nedge, nnode, iarray, k, nstack, stack, &
        maxstack, iendpt, ierror, ncan )

      if ( ierror /= 0 ) then
        signal = 0
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SPAN_TREE_NEXT(): Fatal error!'
        write ( *, '(a,i8)' ) '  Stack needs at least ', nstack
        write ( *, '(a,i8)' ) '  Available space is ', maxstack
        exit
      end if
!
!  No more found.
!
    else

      signal = 0
      exit

    end if

  end do

  return
end
subroutine tqlrat ( n, d, e2, ierr )

!*****************************************************************************80
!
!! tqlrat() computes all eigenvalues of a real symmetric tridiagonal matrix.
!
!  Discussion:
!
!    This subroutine finds the eigenvalues of a symmetric
!    tridiagonal matrix by the rational QL method.
!
!  Reference:
!
!    Christian Reinsch,
!    Algorithm 464, TQLRAT,
!    Communications of the ACM,
!    Volume 16, page 689, 1973.
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Y Ikebe, V Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real ( kind = rk ) D(N).  On input, D contains the diagonal
!    elements of the matrix.  On output, D contains the eigenvalues in ascending
!    order.  If an error exit was made, then the eigenvalues are correct
!    in positions 1 through IERR-1, but may not be the smallest eigenvalues.
!
!    Input/output, real ( kind = rk ) E2(N), contains in positions 2 through N
!    the squares of the subdiagonal elements of the matrix.  E2(1) is
!    arbitrary.  On output, E2 has been overwritten by workspace
!    information.
!
!    Output, integer IERR, error flag.
!    0, for no error,
!    J, if the J-th eigenvalue could not be determined after 30 iterations.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) d(n)
  real ( kind = rk ) e2(n)
  real ( kind = rk ) f
  real ( kind = rk ) g
  real ( kind = rk ) h
  integer i
  integer ierr
  integer ii
  integer j
  integer l
  integer l1
  integer m
  integer mml
  real ( kind = rk ) p
  real ( kind = rk ) pythag
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) t

  ierr = 0

  if ( n == 1 ) then
    return
  end if

  do i = 2, n
    e2(i-1) = e2(i)
  end do

  f = 0.0D+00
  t = 0.0D+00
  e2(n) = 0.0D+00

  do l = 1, n

     j = 0
     h = abs ( d(l) ) + sqrt ( e2(l) )

     if ( t <= h ) then

       t = h
       b = abs ( t ) * epsilon ( b )
       c = b * b

     end if
!
!  Look for small squared sub-diagonal element.
!
     do m = l, n
       if ( e2(m) <= c ) then
         exit
       end if
     end do

     if ( m == l ) then
       go to 210
     end if

130  continue

     if ( 30 <= j ) then
       ierr = l
       return
     end if

     j = j + 1
!
!  Form shift.
!
     l1 = l + 1
     s = sqrt ( e2(l) )
     g = d(l)
     p = ( d(l1) - g ) / ( 2.0D+00 * s )
     r = pythag ( p, 1.0D+00 )
     d(l) = s / ( p + sign ( r, p ) )
     h = g - d(l)

     do i = l1, n
       d(i) = d(i) - h
     end do

     f = f + h
!
!  Rational QL transformation.
!
     g = d(m)
     if ( g == 0.0D+00 ) g = b
     h = g
     s = 0.0D+00
     mml = m - l

     do ii = 1, mml
       i = m - ii
       p = g * h
       r = p + e2(i)
       e2(i+1) = s * r
       s = e2(i) / r
       d(i+1) = h + s * (h + d(i))
       g = d(i) - e2(i) / g
       if ( g == 0.0D+00 ) then
         g = b
       end if
       h = g * p / r
     end do

     e2(l) = s * g
     d(l) = h
!
!  Guard against underflow in convergence test.
!
     if ( h == 0.0D+00 ) go to 210
     if ( abs ( e2(l) ) <= abs ( c / h ) ) go to 210
     e2(l) = h * e2(l)
     if ( e2(l) /= 0.0D+00 ) go to 130

210  continue

     p = d(l) + f
!
!  Order the eigenvalues.
!
     do ii = 2, l
       i = l + 2 - ii
       if ( d(i-1) <= p ) then
         go to 270
       end if
       d(i) = d(i-1)
     end do

     i = 1
270  continue
     d(i) = p

  end do

  return
end
subroutine tred1 ( nm, n, a, d, e, e2 )

!*****************************************************************************80
!
!! tred1() transforms a real symmetric matrix to tridiagonal form.
!
!  Discussion:
!
!    The routine reduces a real symmetric matrix to a symmetric
!    tridiagonal matrix using orthogonal similarity transformations.
!
!  Reference:
!
!    Martin, Reinsch, James Wilkinson,
!    TRED1,
!    Numerische Mathematik,
!    Volume 11, pages 181-195, 1968.
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow, 
!    Y Ikebe, V Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer NM, the leading dimension of the array A.
!    NM must be at least N.
!
!    Input, integer N, the order of the matrix A.
!
!    Input/output, real ( kind = rk ) A(NM,N), on input, contains the real
!    symmetric matrix.  Only the lower triangle of the matrix need be supplied.
!    On output, A contains information about the orthogonal transformations
!    used in the reduction in its strict lower triangle.
!    The full upper triangle of A is unaltered.
!
!    Output, real ( kind = rk ) D(N), contains the diagonal elements of the
!    tridiagonal matrix.
!
!    Output, real ( kind = rk ) E(N), contains the subdiagonal elements of the
!    tridiagonal matrix in its last N-1 positions.  E(1) is set to zero.
!
!    Output, real ( kind = rk ) E2(N), contains the squares of the corresponding
!    elements of E.  E2 may coincide with E if the squares are not needed.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nm

  real ( kind = rk ) a(nm,n)
  real ( kind = rk ) d(n)
  real ( kind = rk ) e(n)
  real ( kind = rk ) e2(n)
  real ( kind = rk ) f
  real ( kind = rk ) g
  real ( kind = rk ) h
  integer i
  integer ii
  integer j
  integer k
  integer l
  real ( kind = rk ) scale

  d(1:n) = a(n,1:n)

  do i = 1, n
    a(n,i) = a(i,i)
  end do

  do ii = 1, n

    i = n + 1 - ii
    l = i - 1
    h = 0.0D+00
    scale = 0.0D+00

    if ( l < 1 ) go to 130
!
!  Scale row.
!
    do k = 1, l
      scale = scale + abs ( d(k) )
    end do

    if ( scale /= 0.0D+00 ) go to 140

    do j = 1, l
      d(j) = a(l,j)
      a(l,j) = a(i,j)
      a(i,j) = 0.0D+00
    end do

130 continue

    e(i) = 0.0D+00
    e2(i) = 0.0D+00
    go to 300

140 continue

    do k = 1, l
      d(k) = d(k) / scale
      h = h + d(k) * d(k)
    end do

    e2(i) = scale * scale * h
    f = d(l)
    g = - sign ( sqrt ( h ), f )
    e(i) = scale * g
    h = h - f * g
    d(l) = f - g

    if ( l == 1 ) go to 285
!
!  Form A*U.
!
    e(1:l) = 0.0D+00

    do j = 1, l

      f = d(j)
      g = e(j) + a(j,j) * f

      do k = j+1, l
        g = g + a(k,j) * d(k)
        e(k) = e(k) + a(k,j) * f
      end do

      e(j) = g

    end do
!
!  Form P.
!
    f = 0.0D+00

    do j = 1, l
      e(j) = e(j) / h
      f = f + e(j) * d(j)
    end do

    h = f / ( h + h )
!
!  Form Q.
!
    e(1:l) = e(1:l) - h * d(1:l)
!
!  Form reduced A.
!
    do j = 1, l

      f = d(j)
      g = e(j)

      do k = j, l
        a(k,j) = a(k,j) - f * e(k) - g * d(k)
      end do

    end do

  285 continue

    do j = 1, l
      f = d(j)
      d(j) = a(l,j)
      a(l,j) = a(i,j)
      a(i,j) = f * scale
    end do

300 continue

  end do

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
subroutine vec_next ( n, iarray, more, ibase )

!*****************************************************************************80
!
!! vec_next() generates all N-vectors of integers modulo a given base.
!
!  Discussion:
!
!    The items are produced one at a time.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 April 1999
!
!  Parameters:
!
!    Input, integer N, the size of the vectors to be used.
!
!    Output, integer IARRAY(N).  On each return from VECNEX,
!    IARRAY will contain entries in the range 0 to IBASE-1.
!
!    Input/output, logical MORE.  Set this variable .FALSE. before
!    the first call.  Normally, MORE will be returned .TRUE. but
!    once all the vectors have been generated, MORE will be
!    reset .FALSE. and you should stop calling the program.
!
!    Input, integer IBASE, the base to be used.  IBASE = 2 will
!    give vectors of 0's and 1's, for instance.
!
  implicit none

  integer n

  integer i
  integer iarray(n)
  integer ibase
  integer, save :: kount
  integer, save :: last
  logical more
  integer nn

  if ( .not. more ) then
 
    kount = 1
    last = ibase**n
    more = .true.
    iarray(1:n) = 0
 
  else
 
    kount = kount + 1

    if ( kount == last ) then
      more = .false.
    end if

    iarray(n) = iarray(n) + 1
 
    do i = 1, n

      nn = n - i

      if ( iarray(nn+1) < ibase ) then
        return
      end if

      iarray(nn+1) = 0

      if ( nn /= 0 ) then
        iarray(nn) = iarray(nn) + 1
      end if

    end do
 
  end if
 
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
subroutine vla_to_graph_arc ( file_name, maxedge, maxnode, nedge, nnode, &
  inode, jnode, x, y, z, ierror )

!*****************************************************************************80
!
!! vla_to_graph_arc() reads graphics information from a VLA file.
!
!  Discussion:
!
!    Internal comments begin with a semi4colon in column 1.
!
!    The X, Y, Z coordinates of points begin with a "P" to
!    denote the beginning of a line, and "L" to denote the
!    continuation of a line.  The fourth entry is intensity, which
!    should be between 0.0 and 1.0.
!
!    It is intended that the information read from the file can
!    either start a whole new graphics object, or simply be added
!    to a current graphics object via the '<<' command.
!
!    This is controlled by whether the input values have been zeroed
!    out or not.  This routine simply tacks on the information it
!    finds to the current graphics object.
!
!  Example:
!
!     set comment fish.vla created by IVREAD
!     set comment from data in file fish.iv
!     set comment
!     set intensity EXPLICIT
!     set parametric NON_PARAMETRIC
!     set filecontent LINES
!     set filetype NEW
!     set depthcue 0
!     set defaultdraw stellar
!     set coordsys RIGHT
!     set author IVREAD
!     set site Buhl Planetarium
!     set library_id UNKNOWN
!     ; DXF LINE entity
!     P   8.59816       5.55317      -3.05561       1.00000
!     L   8.59816       2.49756      0.000000D+00   1.00000
!     L   8.59816       2.49756      -3.05561       1.00000
!     L   8.59816       5.55317      -3.05561       1.00000
!     ; DXF LINE entity
!     P   8.59816       5.55317      0.000000D+00   1.00000
!     ...etc...
!     L   2.48695       2.49756      -3.05561       1.00000
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    13 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the file name.
!
!    Input, integer MAXEDGE, the maximum number of edges.
!
!    Input, integer MAXNODE, the maximum number of nodes.
!
!    Output, integer NEDGE, the number of edges.
!
!    Output, integer NNODE, the number of nodes.
!
!    Output, integer INODE(MAXEDGE), JNODE(MAXEDGE), node pairs 
!    of each edge.
!
!    Output, real ( kind = rk ) X(MAXNODE), Y(MAXNODE), Z(MAXNODE),
!    the coordinates of nodes.
!
!    Output, integer IERROR, 0 no error, 1 an error.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxnode
  integer maxedge

  logical done
  character ( len = * ) file_name
  integer i
  integer icor3
  integer icor3_old
  integer iedge
  integer ierror
  integer indx(maxnode)
  integer indx_edge(maxedge)
  integer inode(maxedge)
  integer ios
  integer iunit
  integer iword
  integer jcor3
  integer jcor3_old
  integer jnode(maxedge)
  integer lchar
  integer nedge
  integer num_bad
  integer num_bad_old
  integer nnode
  real ( kind = rk ) rval
  logical s_eqi
  character ( len = 255 ) text
  character ( len = 255 ) word
  character ( len = 255 ) word1
  real ( kind = rk ) x(maxnode)
  real ( kind = rk ) xval
  real ( kind = rk ) y(maxnode)
  real ( kind = rk ) yval
  real ( kind = rk ) z(maxnode)
  real ( kind = rk ) zval

  ierror = 0

  call get_unit ( iunit )

  open ( unit = iunit, file = file_name, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VLA_TO_GRAPH_ARC(): Fatal error!'
    write ( *, '(a)' ) '  Could not open the file.'
    return
  end if

  ierror = 0
  icor3 = 0
  jcor3 = 0
  num_bad = 0
  num_bad_old = 0
  nedge = 0
  nnode = 0
!
!  Read the next line.
!
  do

    read ( iunit, '(a)', iostat = ios ) text

    if ( ios /= 0 ) then
      exit
    end if

    done = .true.
    iword = 0
!
!  Read the next word.
!
    do

      call word_next_read ( text, word, done )
!
!  If no more words in this line, read a new line.
!
      if ( done ) then
        exit
      end if

      iword = iword + 1
!
!  The first word in the line tells us what's happening.
!
      if ( iword == 1 ) then
        word1 = word
      end if
!
!  If WORD1 is "SET", then we regard this line as comments.
!
      if ( s_eqi ( word1, 'set' ) ) then
!
!  If WORD1 is ";", then we regard this line as comments.
!
      else if ( word1 == ';' ) then
!
!  If WORD1 is "P", then this is the initial point on a line.
!  If WORD1 is "L", then this is a followup point on a line.
!
      else if ( s_eqi ( word1, 'P' ) .or. s_eqi ( word1, 'L' ) ) then
!
!  Read in the point coordinates.
!
        num_bad_old = num_bad

        do i = 1, 3

          call word_next_read ( text, word, done )

          if ( done ) then
            num_bad = num_bad + 1
            exit
          end if

          call s_to_r8 ( word, rval, ierror, lchar )

          if ( ierror /= 0 ) then
            num_bad = num_bad + 1
            exit
          end if

          if ( nnode <= maxnode ) then
            if ( i == 1 ) then
              xval = rval
            else if ( i == 2 ) then
              yval = rval
            else if ( i == 3 ) then
              zval = rval
            end if
          end if

        end do

        if ( num_bad_old < num_bad ) then
          exit
        end if
!
!  Assign a node index to the point.
!
        icor3_old = icor3
        jcor3_old = jcor3
!
!  ICOR3 is the index of the new value.
!  (If such a point already exists, a new one won't be added.)
!
        call r8vec3_index_insert_unique ( maxnode, nnode, x, y, z, indx, &
          xval, yval, zval, icor3, ierror )

        jcor3 = indx(icor3)

        if ( ierror /= 0 ) then
          ierror = 1
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'VLA_TO_GRAPH_ARC(): Fatal error!'
          write ( *, '(a)' ) '  R8VEC3_INDEX_INSERT_UNIQUE returned an error!'
          return
        end if
!
!  Define the line as joining JCOR3_OLD to JCOR3.
!  (If such a line already exists, a new copy won't be added.)
!
        if ( s_eqi ( word1, 'L' ) ) then

          call iset2_index_insert_unique ( maxedge, nedge, inode, jnode, &
            indx_edge, jcor3_old, jcor3, iedge, ierror )

        end if

        exit
!
!  If the first word is unrecognized, then skip the whole line.
!
      else

        num_bad = num_bad + 1
        exit

      end if

    end do

  end do

  close ( unit = iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'VLA_TO_GRAPH_ARC - Note:'
  write ( *, '(a)' ) '  The graph was read properly.'
  write ( *, '(a,i8)' ) '  Number of nodes = ', nnode
  write ( *, '(a,i8)' ) '  Number of edges = ', nedge

  return
end
subroutine word_next_read ( line, word, done )

!*****************************************************************************80
!
!! word_next_read() "reads" words from a string, one at a time.
!
!  Discussion:
!
!    The following characters are considered to be a single word,
!    whether surrounded by spaces or not:
!
!      " ( ) { } [ ]
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) LINE, a string, presumably containing words
!    separated by spaces.
!
!    Output, character ( len = * ) WORD.
!
!    If DONE is FALSE, then WORD contains the "next" word read from LINE.
!    If DONE is TRUE, then WORD is blank, because there was no more to read.
!
!    Input/output, logical DONE.
!
!    On input with a fresh value of LINE, set DONE to TRUE.
!
!    On output, the routine sets DONE:
!      FALSE if another word was read from LINE,
!      TRUE if no more words could be read (LINE is exhausted).
!
  implicit none

  logical done
  integer ilo
  integer, save :: lenc = 0
  character ( len = * ) line
  integer, save :: next = 1
  character TAB
  character ( len = * ) word

  TAB = char ( 9 )
!
!  An input value of DONE = TRUE signals a new line of text to examine.
!
  if ( done ) then

    next = 1
    done = .false.
    lenc = len_trim ( line )

    if ( lenc <= 0 ) then
      done = .true.
      word = ' '
      return
    end if

  end if
!
!  Beginning at index NEXT, search LINE for the next nonblank,
!  which signals the beginning of a word.
!
  ilo = next

10    continue
!
!  ...LINE(NEXT:) is blank.  Return with WORD = ' ' and DONE = TRUE.
!
  if ( lenc < ilo ) then
    word = ' '
    done = .true.
    next = lenc + 1
    return
  end if
!
!  If the current character is blank, skip to the next one.
!
  if ( line(ilo:ilo) == ' ' .or. line(ilo:ilo) == TAB ) then
    ilo = ilo + 1
    go to 10
  end if
!
!  ILO is the index of the next nonblank character in the string.
!
!  If this initial nonblank is a special character,
!  then that's the whole word as far as we're concerned,
!  so return immediately.
!
  if ( line(ilo:ilo) == '"' .or. &
       line(ilo:ilo) == '(' .or. &
       line(ilo:ilo) == ')' .or. &
       line(ilo:ilo) == '{' .or. &
       line(ilo:ilo) == '}' .or. &
       line(ilo:ilo) == '[' .or. &
       line(ilo:ilo) == ']' ) then

    word = line(ilo:ilo)
    next = ilo + 1
    return

  end if
!
!  Now search for the last contiguous character that is not a
!  blank, TAB, or special character.
!
  next = ilo + 1

20    continue

  if ( lenc < next ) then
    word = line(ilo:next-1)
    return
  end if

  if ( line(next:next) /= ' ' .and. &
       line(next:next) /= TAB .and. &
       line(next:next) /= '"' .and. &
       line(next:next) /= '(' .and. &
       line(next:next) /= ')' .and. &
       line(next:next) /= '{' .and. &
       line(next:next) /= '}' .and. &
       line(next:next) /= '[' .and. &
       line(next:next) /= ']' ) then

    next = next + 1
    go to 20

  end if

  word = line(ilo:next-1)

  return
end
