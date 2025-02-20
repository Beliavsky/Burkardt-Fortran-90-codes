subroutine comp_next_grlex ( kc, xc )

!*****************************************************************************80
!
!! comp_next_grlex() returns the next composition in grlex order.
!
!  Discussion:
!
!    Example:
!
!    KC = 3
!
!    #   XC(1  XC(2) XC(3)  Degree
!      +------------------------
!    1 |  0     0     0        0
!      |
!    2 |  0     0     1        1
!    3 |  0     1     0        1
!    4 |  1     0     0        1
!      |
!    5 |  0     0     2        2
!    6 |  0     1     1        2
!    7 |  0     2     0        2
!    8 |  1     0     1        2
!    9 |  1     1     0        2
!   10 |  2     0     0        2
!      |
!   11 |  0     0     3        3
!   12 |  0     1     2        3
!   13 |  0     2     1        3
!   14 |  0     3     0        3
!   15 |  1     0     2        3
!   16 |  1     1     1        3
!   17 |  1     2     0        3
!   18 |  2     0     1        3
!   19 |  2     1     0        3
!   20 |  3     0     0        3
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 December 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer KC, the number of parts of the composition.
!    1 <= KC.
!
!    Input/output, integer XC(KC), the current composition.
!    Each entry of XC must be nonnegative.
!    On return, XC has been replaced by the next composition in the
!    grlex order.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer kc

  integer i
  integer im1
  integer j
  integer t
  integer xc(kc)
!
!  Ensure that 1 <= KC.
!
  if ( kc < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COMP_NEXT_GRLEX - Fatal error!'
    write ( *, '(a)' ) '  KC < 1'
    stop 1
  end if
!
!  Ensure that 0 <= XC(I).
!
  do i = 1, kc
    if ( xc(i) < 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'COMP_NEXT_GRLEX - Fatal error!'
      write ( *, '(a)' ) '  XC(I) < 0'
      stop 1
    end if
  end do
!
!  Find I, the index of the rightmost nonzero entry of X.
!
  i = 0
  do j = kc, 1, -1
    if ( 0 < xc(j) ) then
      i = j
      exit
    end if
  end do
!
!  set T = X(I)
!  set XC(I) to zero,
!  increase XC(I-1) by 1,
!  increment XC(KC) by T-1.
!
  if ( i == 0 ) then
    xc(kc) = 1
    return
  else if ( i == 1 ) then
    t = xc(1) + 1
    im1 = kc
  else if ( 1 < i ) then
    t = xc(i)
    im1 = i - 1
  end if

  xc(i) = 0
  xc(im1) = xc(im1) + 1
  xc(kc) = xc(kc) + t - 1

  return
end
subroutine comp_random ( n, k, a )

!*****************************************************************************80
!
!! COMP_RANDOM selects a random composition of the integer N into K parts.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 June 2024
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    This version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Second Edition,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer N, the integer to be decomposed.
!
!    Input, integer K, the number of parts in the composition.
!
!    Output, integer A(K), the parts of the composition.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer k

  integer a(k)
  integer i
  integer l
  integer m
  integer n

  call ksub_random ( n + k - 1, k - 1, a )

  a(k) = n + k
  l = 0

  do i = 1, k
    m = a(i)
    a(i) = a(i) - l - 1
    l = m
  end do

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
!    01 June 2024
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
  integer, parameter :: i4_huge = 2147483647
  integer i4_uniform_ab
  real r
  integer value

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

  i4_uniform_ab = value

  return
end
subroutine i4mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! I4MAT_TRANSPOSE_PRINT prints an I4MAT, transposed.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer A(M,N), an M by N matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer a(m,n)
  character ( len = * ) title

  call i4mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine i4mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! I4MAT_TRANSPOSE_PRINT_SOME prints some of the transpose of an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer A(M,N), an M by N matrix to be printed.
!
!    Input, integer ILO, JLO, the first row and column to print.
!
!    Input, integer IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: incx = 10
  integer m
  integer n

  integer a(m,n)
  character ( len = 8 ) ctemp(incx)
  integer i
  integer i2
  integer i2hi
  integer i2lo
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)'
    return
  end if

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(i8)' ) i
    end do

    write ( *, '(''  Row '',10a8)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc

        i = i2lo - 1 + i2

        write ( ctemp(i2), '(i8)' ) a(i,j)

      end do

      write ( *, '(i5,a,10a8)' ) j, ':', ( ctemp(i), i = 1, inc )

    end do

  end do

  return
end
subroutine ksub_random ( n, k, a )

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
!    01 June 2024
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    This version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Second Edition,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer N, the size of the set from which subsets
!    are drawn.
!
!    Input, integer K, number of elements in desired subsets.
!    K must be between 0 and N.
!
!    Output, integer A(K).  A(I) is the I-th element of the
!    output set.  The elements of A are in order.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer k

  integer a(k)
  integer i
  integer i4_uniform_ab
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
    write ( *, '(a)' ) 'KSUB_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  K = ', k
    write ( *, '(a)' ) '  but 0 <= K is required!'
    stop
  else if ( n < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUB_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  N = ', n
    write ( *, '(a,i8)' ) '  K = ', k
    write ( *, '(a)' ) '  K <= N is required!'
    stop
  end if

  if ( k == 0 ) then
    return
  end if

  do i = 1, k
    a(i) = ( ( i - 1 ) * n ) / k
  end do

  do i = 1, k

    do

      ix = i4_uniform_ab ( 1, n )

      l = 1 + ( ix * k - 1 ) / n

      if ( a(l) < ix ) then
        exit
      end if

    end do

    a(l) = a(l) + 1

  end do

  ip = 0
  is = k

  do i = 1, k

    m = a(i)
    a(i) = 0

    if ( m /= ( ( i - 1 ) * n ) / k ) then
      ip = ip + 1
      a(ip) = m
    end if

  end do

  ihi = ip

  do i = 1, ihi
    ip = ihi + 1 - i
    l = 1 + ( a(ip) * k - 1 ) / n
    ids = a(ip) - ( ( l - 1 ) * n ) / k
    a(ip) = 0
    a(is) = l
    is = is - ids
  end do

  do ll = 1, k

    l = k + 1 - ll

    if ( a(l) /= 0 ) then
      ir = l
      m0 = 1 + ( ( a(l) - 1 ) * n ) / k
      m = ( a(l) * n ) / k - m0 + 1
    end if

    ix = i4_uniform_ab ( m0, m0 + m - 1 )

    i = l + 1

    do while ( i <= ir )

      if ( ix < a(i) ) then
        exit
      end if

      ix = ix + 1
      a(i-1) = a(i)
      i = i + 1

    end do

    a(i-1) = ix
    m = m - 1

  end do

  return
end
subroutine r8mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = rk ) A(M,N), an M by N matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  character ( len = * ) title

  call r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = rk ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ILO, JLO, the first row and column to print.
!
!    Input, integer IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: incx = 5
  integer m
  integer n

  real ( kind = rk ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer i
  integer i2
  integer i2hi
  integer i2lo
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)'
    return
  end if

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(i8,6x)' ) i
    end do

    write ( *, '(''  Row   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc
        i = i2lo - 1 + i2
        write ( ctemp(i2), '(g14.6)' ) a(i,j)
      end do

      write ( *, '(i5,a,5a14)' ) j, ':', ( ctemp(i), i = 1, inc )

    end do

  end do

  return
end
subroutine simplex_grid_index_all ( m, n, ng, grid )

!*****************************************************************************80
!
!! SIMPLEX_GRID_INDEX_ALL returns all the simplex grid indices.
!
!  Discussion:
!
!    The number of grid indices can be determined by calling 
!      ng = simplex_grid_size ( m, n )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of subintervals.
!
!    Input, integer NG, the number of values in the grid.
!
!    Output, integer GRID(M+1,NG), the current, and then the next,
!    grid index.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer ng

  integer g(m+1)
  integer grid(m+1,ng)
  integer i
  integer k
  integer n

  do i = 1, m
    g(i) = 0
  end do
  g(m+1) = n

  k = 1
  grid(1:m+1,k) = g(1:m+1)

  do while ( k < ng )
    call comp_next_grlex ( m + 1, g )
    k = k + 1
    grid(1:m+1,k) = g(1:m+1)
  end do

  return
end
subroutine simplex_grid_index_next ( m, n, g )

!*****************************************************************************80
!
!! SIMPLEX_GRID_INDEX_NEXT returns the next simplex grid index.
!
!  Discussion:
!
!    The vector G has dimension M+1.  The first M entries may be regarded
!    as grid coordinates.  These coordinates must have a sum between 0 and N.
!    The M+1 entry contains the remainder, that is N minus the sum of the
!    first M coordinates.
!
!    Each time the function is called, it is given a current grid index, and
!    computes the next one.  The very first index is all zero except for a 
!    final value of N, and the very last index has all zero except for an'
!    intial value of N.
!
!    For example, here are the coordinates in order for M = 3, N = 3:
!
!     0  0 0 0 3
!     1  0 0 1 2
!     2  0 0 2 1
!     3  0 0 3 0
!     4  0 1 0 2
!     5  0 1 1 1
!     6  0 1 2 0
!     7  0 2 0 1
!     8  0 2 1 0
!     9  0 3 0 0
!    10  1 0 0 2
!    11  1 0 1 1
!    12  1 0 2 0
!    13  1 1 0 1
!    14  1 1 1 0
!    15  1 2 0 0
!    16  2 0 0 1
!    17  2 0 1 0
!    18  2 1 0 0
!    19  3 0 0 0
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of subintervals.
!
!    Input/output, integer G(M+1), the current, and then the next,
!    grid index.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m

  integer g(m+1)
  integer n

  call comp_next_grlex ( m + 1, g )

  return
end
subroutine simplex_grid_index_sample ( m, n, g )

!*****************************************************************************80
!
!! simplex_grid_index_sample() returns a random simplex grid index.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 June 2024
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of subintervals in
!    each dimension.
!
!    Output, integer G(M+1), a randomly selected index in the 
!    simplex grid.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m

  integer g(m+1)
  integer n

  call comp_random ( n, m + 1, g )

  return
end
subroutine simplex_grid_index_to_point ( m, n, ng, g, v, x )

!*****************************************************************************80
!
!! simplex_grid_index_to_point() returns  points corresponding to simplex indices.
!
!  Discussion:
!
!    The M-dimensional simplex is defined by M+1 vertices.
!
!    Given a regular grid that uses N subintervals along the edge between
!    each pair of vertices, a simplex grid index G is a set of M+1 values
!    each between 0 and N, and summing to N. 
!
!    This function determines the coordinates X of the point corresponding
!    to the index G.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of subintervals.
!
!    Input, integer NG, the number of grid indices to be converted.
!
!    Input, integer G(M+1,NG), the grid indices of 1 
!    or more points.
!
!    Input, real ( kind = rk ) V(M,M+1), the coordinates of the vertices 
!    of the simplex.
!
!    Output, real ( kind = rk ) X(M,NG), the coordinates of one or more points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n
  integer ng

  integer g(m+1,ng)
  integer i
  integer j
  integer k
  real ( kind = rk ) v(m,m+1)
  real ( kind = rk ) x(m,ng)
  
  do j = 1, ng
    do i = 1, m
      x(i,j) = 0.0D+00
      do k = 1, m + 1
        x(i,j) = x(i,j) + v(i,k) * real ( g(k,j), kind = rk )
      end do
      x(i,j) = x(i,j) / real ( n, kind = rk )
    end do
  end do

  return
end
subroutine simplex_grid_size ( m, n, ng )

!*****************************************************************************80
!
!! SIMPLEX_GRID_SIZE counts the grid points inside a simplex.
!
!  Discussion:
!
!    The size of a grid with parameters M, N is C(M+N,N).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of subintervals.
!
!    Output, integer NG, the number of grid points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer m
  integer n
  integer ng

  ng = 1

  do i = 1, m
    ng = ( ng * ( n + i ) ) / i
  end do

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
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
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
