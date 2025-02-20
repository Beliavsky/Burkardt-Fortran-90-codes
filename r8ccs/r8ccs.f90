subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! get_unit() returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IUNIT, the free unit number.
!
  implicit none

  integer i
  integer ios
  integer iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
function i4_log_10 ( i )

!*****************************************************************************80
!
!! I4_LOG_10 returns the integer part of the logarithm base 10 of an I4.
!
!  Example:
!
!        I  I4_LOG_10
!    -----  --------
!        0    0
!        1    0
!        2    0
!        9    0
!       10    1
!       11    1
!       99    1
!      100    2
!      101    2
!      999    2
!     1000    3
!     1001    3
!     9999    3
!    10000    4
!
!  Discussion:
!
!    I4_LOG_10 ( I ) + 1 is the number of decimal digits in I.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    08 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, the number whose logarithm base 10 
!    is desired.
!
!    Output, integer I4_LOG_10, the integer part of the 
!    logarithm base 10 of the absolute value of X.
!
  implicit none

  integer i
  integer i_abs
  integer i4_log_10
  integer ten_pow

  if ( i == 0 ) then

    i4_log_10 = 0

  else

    i4_log_10 = 0
    ten_pow = 10

    i_abs = abs ( i )

    do while ( ten_pow <= i_abs )
      i4_log_10 = i4_log_10 + 1
      ten_pow = ten_pow * 10
    end do

  end if

  return
end
function i4_uniform_ab ( a, b, seed )

!*****************************************************************************80
!
!! I4_UNIFORM_AB returns a scaled pseudorandom I4 between A and B.
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
!  Parameters:
!
!    Input, integer A, B, the limits of the interval.
!
!    Input/output, integer SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer I4_UNIFORM_AB, a number between A and B.
!
  implicit none

  integer a
  integer b
  integer, parameter :: i4_huge = 2147483647
  integer i4_uniform_ab
  integer k
  real r
  integer seed
  integer value

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_UNIFORM_AB - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r = real ( seed ) * 4.656612875E-10
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
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 May 2010
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
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer n

  integer a(n)
  integer i
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,a,2x,i12)' ) i, ':', a(i)
  end do

  return
end
subroutine i4vec_search_binary_a ( n, a, b, indx )

!*****************************************************************************80
!
!! I4VEC_SEARCH_BINARY_A searches an ascending sorted I4VEC for a value.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    Binary search is used.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Algorithm 1.9,
!    Combinatorial Algorithms,
!    CRC Press, 1998, page 26.
!
!  Parameters:
!
!    Input, integer N, the number of elements in the vector.
!
!    Input, integer A(N), the array to be searched.  A must
!    be sorted in ascending order.
!
!    Input, integer B, the value to be searched for.
!
!    Output, integer INDX, the result of the search.
!    -1, B does not occur in A.
!    I, A(I) = B.
!
  implicit none

  integer n

  integer a(n)
  integer b
  integer high
  integer indx
  integer low
  integer mid

  indx = - 1

  low = 1
  high = n

  do while ( low <= high )

    mid = ( low + high ) / 2

    if ( a(mid) == b ) then
      indx = mid
      exit
    else if ( a(mid) < b ) then
      low = mid + 1
    else if ( b < a(mid) ) then
      high = mid - 1
    end if

  end do

  return
end
subroutine r8ccs_dif2 ( m, n, nz_num, col, row, a )

!*****************************************************************************80
!
!! r8ccs_DIF2 sets the second difference as an r8ccs matrix.
!
!  Discussion:
!
!    The r8ccs format is the double precision sparse compressed column
!    format.  Associated with this format, we have an M by N matrix
!    with NZ_NUM nonzero entries.  We construct the column pointer
!    vector COL of length N+1, such that entries of column J will be
!    stored in positions COL(J) through COL(J+1)-1.  This indexing
!    refers to both the ROW and A vectors, which store the row indices
!    and the values of the nonzero entries.  The entries of the
!    ROW vector corresponding to each column are assumed to be
!    ascending sorted.
!
!    The r8ccs format is equivalent to the MATLAB "sparse" format,
!    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 October 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows of the matrix.
!
!    Input, integer N, the number of columns of the matrix.
!
!    Input, integer NZ_NUM, the number of nonzero entries.
!
!    Output, integer COL(N+1), indicate where each column's data
!    begins.
!
!    Output, integer ROW(NZ_NUM), the row indices.
!
!    Output, real ( kind = rk ) A(NZ_NUM), the nonzero entries.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nz_num

  real ( kind = rk ) a(nz_num)
  integer col(n+1)
  integer i
  integer j
  integer k
  integer m
  integer row(nz_num)
!
!  Column pointers
!
  col(1) = 1
  col(2) = 3
  do j = 3, n
    col(j) = col(j-1) + 3
  end do
  col(n+1) = col(n) + 2
!
!  Row indices
!
  k = 0
  k = k + 1
  row(k) = 1
  k = k + 1
  row(k) = 2
  do j = 2, n - 1
    do i = j - 1, j + 1
      k = k + 1
      row(k) = i
    end do
  end do
  k = k + 1
  row(k) = m - 1
  k = k + 1
  row(k) = m
!
!  Values
!
  k = 0

  j = 1
  i = 1
  k = k + 1
  a(k) = 2.0D+00
  i = 2
  k = k + 1
  a(k) = -1.0D+00

  do j = 2, n - 1
    i = j - 1
    k = k + 1
    a(k) = -1.0D+00
    i = j
    k = k + 1
    a(k) =  2.0D+00
    i = j + 1
    k = k + 1
    a(k) = -1.0D+00
  end do

  j = n
  i = m - 1
  k = k + 1
  a(k) = -1.0D+00
  i = m
  k = k + 1
  a(k) = 2.0D+00

  return
end
subroutine r8ccs_get ( m, n, nz_num, col, row, a, i, j, aij )

!*****************************************************************************80
!
!! r8ccs_GET gets a value of an r8ccs matrix.
!
!  Discussion:
!
!    It is legal to request entries of the matrix for which no storage
!    was set aside.  In that case, a zero value will be returned.
!
!    The r8ccs format is the double precision sparse compressed column
!    format.  Associated with this format, we have an M by N matrix
!    with NZ_NUM nonzero entries.  We construct the column pointer
!    vector COL of length N+1, such that entries of column J will be
!    stored in positions COL(J) through COL(J+1)-1.  This indexing
!    refers to both the ROW and A vectors, which store the row indices
!    and the values of the nonzero entries.  The entries of the
!    ROW vector corresponding to each column are assumed to be
!    ascending sorted.
!
!    The r8ccs format is equivalent to the MATLAB "sparse" format,
!    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    October 1992
!
!  Parameters:
!
!    Input, integer M, the number of rows of the matrix.
!
!    Input, integer N, the number of columns of the matrix.
!
!    Input, integer NZ_NUM, the number of nonzero entries.
!
!    Input, integer COL(N+1), indicate where each column's data
!    begins.
!
!    Input, integer ROW(NZ_NUM), the row indices.
!
!    Input, real ( kind = rk ) A(NZ_NUM), the nonzero entries.
!
!    Input, integer I, J, the indices of the value to retrieve.
!
!    Output, real ( kind = rk ) AIJ, the value of A(I,J).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nz_num

  real ( kind = rk ) a(nz_num)
  real ( kind = rk ) aij
  integer col(n+1)
  integer i
  integer j
  integer k
  integer m
  integer row(nz_num)
!
!  Seek sparse index K corresponding to full index (I,J).
!
  call r8ccs_ijk ( m, n, nz_num, col, row, i, j, k )
!
!  If no K was found, then be merciful, and simply return 0.
!
  if ( k == -1 ) then
    aij = 0.0D+00
  else
    aij = a(k)
  end if

  return
end
subroutine r8ccs_ijk ( m, n, nz_num, col, row, i, j, k )

!*****************************************************************************80
!
!! r8ccs_IJK seeks the sparse index K of (I,J), the full index of an r8ccs matrix.
!
!  Discussion:
!
!    The r8ccs format is the double precision sparse compressed column
!    format.  Associated with this format, we have an M by N matrix
!    with NZ_NUM nonzero entries.  We construct the column pointer
!    vector COL of length N+1, such that entries of column J will be
!    stored in positions COL(J) through COL(J+1)-1.  This indexing
!    refers to both the ROW and A vectors, which store the row indices
!    and the values of the nonzero entries.  The entries of the
!    ROW vector corresponding to each column are assumed to be
!    ascending sorted.
!
!    The r8ccs format is equivalent to the MATLAB "sparse" format,
!    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    October 1992
!
!  Parameters:
!
!    Input, integer M, the number of rows of the matrix.
!
!    Input, integer N, the number of columns of the matrix.
!
!    Input, integer NZ_NUM, the number of nonzero entries.
!
!    Input, integer COL(N+1), indicate where each column's 
!    data begins.
!
!    Input, integer ROW(NZ_NUM), the row indices.
!
!    Input, integer I, J, the indices of the value to retrieve.
!
!    Output, integer K, the index of the sparse matrix in which 
!    entry (I,J) is stored, or -1 if no such entry exists.
!
  implicit none

  integer n
  integer nz_num

  integer col(n+1)
  integer i
  integer j
  integer k
  integer k1
  integer k2
  integer m
  integer row(nz_num)
!
!  Determine the part of ROW containing row indices of entries
!  in column J.
!
  k1 = col(j)
  k2 = col(j+1)-1
!
!  Seek the location K for which ROW(K) = I.
!  
  call i4vec_search_binary_a ( k2+1-k1, row(k1:k2), i, k )

  if ( k /= -1 ) then
    k = k + k1 - 1
  end if

  return
end
subroutine r8ccs_inc ( m, n, nz_num, col, row, a, i, j, aij )

!*****************************************************************************80
!
!! r8ccs_INC increments a value of an r8ccs matrix.
!
!  Discussion:
!
!    The r8ccs format is the double precision sparse compressed column
!    format.  Associated with this format, we have an M by N matrix
!    with NZ_NUM nonzero entries.  We construct the column pointer
!    vector COL of length N+1, such that entries of column J will be
!    stored in positions COL(J) through COL(J+1)-1.  This indexing
!    refers to both the ROW and A vectors, which store the row indices
!    and the values of the nonzero entries.  The entries of the
!    ROW vector corresponding to each column are assumed to be
!    ascending sorted.
!
!    The r8ccs format is equivalent to the MATLAB "sparse" format,
!    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    October 1992
!
!  Parameters:
!
!    Input, integer M, the number of rows of the matrix.
!
!    Input, integer N, the number of columns of the matrix.
!
!    Input, integer NZ_NUM, the number of nonzero entries.
!
!    Input, integer COL(N+1), indicate where each column's 
!    data begins.
!
!    Input, integer ROW(NZ_NUM), the row indices.
!
!    Input/output, real ( kind = rk ) A(NZ_NUM), the nonzero entries.
!    On output, entry (I,J) has been incremented.
!
!    Input, integer I, J, the indices of the value to retrieve.
!
!    Input, real ( kind = rk ) AIJ, the value to be added to A(I,J).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nz_num

  real ( kind = rk ) a(nz_num)
  real ( kind = rk ) aij
  integer col(n+1)
  integer i
  integer j
  integer k
  integer m
  integer row(nz_num)
!
!  Seek sparse index K corresponding to full index (I,J).
!
  call r8ccs_ijk ( m, n, nz_num, col, row, i, j, k )
!
!  If no K was found, we fail.
!
  if ( k == -1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'r8ccs_INC - Fatal error!'
    write ( *, '(a)' ) '  r8ccs_IJK could not find the entry.'
    write ( *, '(a,i8)' ) '  Row I = ', i
    write ( *, '(a,i8)' ) '  Col J = ', j
    stop 1
  end if

  a(k) = a(k) + aij

  return
end
subroutine r8ccs_indicator ( m, n, nz_num, col, row, a )

!*****************************************************************************80
!
!! r8ccs_INDICATOR sets up an r8ccs indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The r8ccs format is the double precision sparse compressed column
!    format.  Associated with this format, we have an M by N matrix
!    with NZ_NUM nonzero entries.  We construct the column pointer
!    vector COL of length N+1, such that entries of column J will be
!    stored in positions COL(J) through COL(J+1)-1.  This indexing
!    refers to both the ROW and A vectors, which store the row indices
!    and the values of the nonzero entries.  The entries of the
!    ROW vector corresponding to each column are assumed to be
!    ascending sorted.
!
!    The r8ccs format is equivalent to the MATLAB "sparse" format,
!    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    October 1992
!
!  Parameters:
!
!    Input, integer M, the number of rows of the matrix.
!
!    Input, integer N, the number of columns of the matrix.
!
!    Input, integer NZ_NUM, the number of nonzero elements in A.
!
!    Input, integer COL(N+1), points to the first element of 
!    each column.
!
!    Input, integer ROW(NZ_NUM), contains the row indices 
!    of the elements.
!
!    Output, real ( kind = rk ) A(NZ_NUM), the r8ccs matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nz_num

  real ( kind = rk ) a(nz_num)
  integer col(n+1)
  integer fac
  integer i
  integer i4_log_10
  integer j
  integer k
  integer m
  integer row(nz_num)

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  do j = 1, n
    do k = col(j), col(j+1) - 1
      i = row(k)
      a(k) = real ( fac * i + j, kind = rk )
    end do
  end do

  return
end
subroutine r8ccs_kij ( m, n, nz_num, col, row, k, i, j )

!*****************************************************************************80
!
!! r8ccs_KIJ seeks the full index (I,J) of K, the sparse index of an r8ccs matrix.
!
!  Discussion:
!
!    The r8ccs format is the double precision sparse compressed column
!    format.  Associated with this format, we have an M by N matrix
!    with NZ_NUM nonzero entries.  We construct the column pointer
!    vector COL of length N+1, such that entries of column J will be
!    stored in positions COL(J) through COL(J+1)-1.  This indexing
!    refers to both the ROW and A vectors, which store the row indices
!    and the values of the nonzero entries.  The entries of the
!    ROW vector corresponding to each column are assumed to be
!    ascending sorted.
!
!    The r8ccs format is equivalent to the MATLAB "sparse" format,
!    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    October 1992
!
!  Parameters:
!
!    Input, integer M, the number of rows of the matrix.
!
!    Input, integer N, the number of columns of the matrix.
!
!    Input, integer NZ_NUM, the number of nonzero entries.
!
!    Input, integer COL(N+1), indicate where each column's data 
!    begins.
!
!    Input, integer ROW(NZ_NUM), the row indices.
!
!    Input, integer K, the sparse index of an entry of the matrix.
!    1 <= K <= NZ_NUM.
!
!    Output, integer I, J, the full indices corresponding to the 
!    sparse index K.
!
  implicit none

  integer n
  integer nz_num

  integer col(n+1)
  integer i
  integer j
  integer jj
  integer k
  integer k1
  integer k2
  integer m
  integer row(nz_num)

  i = -1
  j = -1

  if ( k < 1 .or. nz_num < k ) then
    return
  end if
!
!  The row index is easy.
!
  i = row(k)
!
!  Determine the column by bracketing in COL.
!
  do jj = 1, n
    k1 = col(jj)
    k2 = col(jj+1)-1
    if ( k1 <= k .and. k <= k2 ) then
      j = jj
      exit
    end if
  end do

  if ( j == -1 ) then
    return
  end if

  return
end
subroutine r8ccs_mtv ( m, n, nz_num, col, row, a, x, b )

!*****************************************************************************80
!
!! r8ccs_MTV multiplies an R8VEC times an r8ccs matrix.
!
!  Discussion:
!
!    The r8ccs format is the double precision sparse compressed column
!    format.  Associated with this format, we have an M by N matrix
!    with NZ_NUM nonzero entries.  We construct the column pointer
!    vector COL of length N+1, such that entries of column J will be
!    stored in positions COL(J) through COL(J+1)-1.  This indexing
!    refers to both the ROW and A vectors, which store the row indices
!    and the values of the nonzero entries.  The entries of the
!    ROW vector corresponding to each column are assumed to be
!    ascending sorted.
!
!    The r8ccs format is equivalent to the MATLAB "sparse" format,
!    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    October 1992
!
!  Parameters:
!
!    Input, integer M, the number of rows of the matrix.
!
!    Input, integer N, the number of columns of the matrix.
!
!    Input, integer NZ_NUM, the number of nonzero elements in A.
!
!    Input, integer COL(N+1), points to the first element 
!    of each column.
!
!    Input, integer ROW(NZ_NUM), contains the row indices 
!    of the elements.
!
!    Input, real ( kind = rk ) A(NZ_NUM), the r8ccs matrix.
!
!    Input, real ( kind = rk ) X(M), the vector to be multiplied.
!
!    Output, real ( kind = rk ) B(N), the product A'*X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n
  integer nz_num

  real ( kind = rk ) a(nz_num)
  real ( kind = rk ) b(n)
  integer col(n+1)
  integer i
  integer j
  integer k
  integer row(nz_num)
  real ( kind = rk ) x(m)

  b(1:n) = 0.0D+00

  do j = 1, n
    do k = col(j), col(j+1) - 1
      i = row(k)
      b(j) = b(j) + a(k) * x(i)
    end do
  end do

  return
end
subroutine r8ccs_mv ( m, n, nz_num, col, row, a, x, b )

!*****************************************************************************80
!
!! r8ccs_MV multiplies an r8ccs matrix by an R8VEC.
!
!  Discussion:
!
!    The r8ccs format is the double precision sparse compressed column
!    format.  Associated with this format, we have an M by N matrix
!    with NZ_NUM nonzero entries.  We construct the column pointer
!    vector COL of length N+1, such that entries of column J will be
!    stored in positions COL(J) through COL(J+1)-1.  This indexing
!    refers to both the ROW and A vectors, which store the row indices
!    and the values of the nonzero entries.  The entries of the
!    ROW vector corresponding to each column are assumed to be
!    ascending sorted.
!
!    The r8ccs format is equivalent to the MATLAB "sparse" format,
!    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    October 1992
!
!  Parameters:
!
!    Input, integer M, the number of rows of the matrix.
!
!    Input, integer N, the number of columns of the matrix.
!
!    Input, integer NZ_NUM, the number of nonzero elements in A.
!
!    Input, integer COL(N+1), points to the first element of 
!    each column.
!
!    Input, integer ROW(NZ_NUM), contains the row indices of 
!    the elements.
!
!    Input, real ( kind = rk ) A(NZ_NUM), the r8ccs matrix.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied.
!
!    Output, real ( kind = rk ) B(M), the product A*X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n
  integer nz_num

  real ( kind = rk ) a(nz_num)
  real ( kind = rk ) b(m)
  integer col(n+1)
  integer i
  integer j
  integer k
  integer row(nz_num)
  real ( kind = rk ) x(n)

  b(1:m) = 0.0D+00

  do j = 1, n
    do k = col(j), col(j+1) - 1
      i = row(k)
      b(i) = b(i) + a(k) * x(j)
    end do
  end do

  return
end
subroutine r8ccs_print ( m, n, nz_num, col, row, a, title )

!*****************************************************************************80
!
!! r8ccs_PRINT prints an r8ccs matrix.
!
!  Discussion:
!
!    The r8ccs format is the double precision sparse compressed column
!    format.  Associated with this format, we have an M by N matrix
!    with NZ_NUM nonzero entries.  We construct the column pointer
!    vector COL of length N+1, such that entries of column J will be
!    stored in positions COL(J) through COL(J+1)-1.  This indexing
!    refers to both the ROW and A vectors, which store the row indices
!    and the values of the nonzero entries.  The entries of the
!    ROW vector corresponding to each column are assumed to be
!    ascending sorted.
!
!    The r8ccs format is equivalent to the MATLAB "sparse" format,
!    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    October 1992
!
!  Parameters:
!
!    Input, integer M, the number of rows of the matrix.
!
!    Input, integer N, the number of columns of the matrix.
!
!    Input, integer NZ_NUM, the number of nonzero elements in A.
!
!    Input, integer COL(N+1), points to the first element of 
!    each column.
!
!    Input, integer ROW(NZ_NUM), contains the row indices of 
!    the elements.
!
!    Input, real ( kind = rk ) A(NZ_NUM), the r8ccs matrix.
!
!    Input, character ( len = * ) TITLE, a title.
! 
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nz_num

  real ( kind = rk ) a(nz_num)
  integer col(n+1)
  integer m
  integer row(nz_num)
  character ( len = * ) title

  call r8ccs_print_some ( m, n, nz_num, col, row, a, 1, 1, n, n, title )

  return
end
subroutine r8ccs_print_some ( m, n, nz_num, col, row, a, ilo, jlo, &
  ihi, jhi, title )

!*****************************************************************************80
!
!! r8ccs_PRINT_SOME prints some of an r8ccs matrix.
!
!  Discussion:
!
!    The r8ccs format is the double precision sparse compressed column
!    format.  Associated with this format, we have an M by N matrix
!    with NZ_NUM nonzero entries.  We construct the column pointer
!    vector COL of length N+1, such that entries of column J will be
!    stored in positions COL(J) through COL(J+1)-1.  This indexing
!    refers to both the ROW and A vectors, which store the row indices
!    and the values of the nonzero entries.  The entries of the
!    ROW vector corresponding to each column are assumed to be
!    ascending sorted.
!
!    The r8ccs format is equivalent to the MATLAB "sparse" format,
!    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    October 1992
!
!  Parameters:
!
!    Input, integer M, the number of rows of the matrix.
!
!    Input, integer N, the number of columns of the matrix.
!
!    Input, integer NZ_NUM, the number of nonzero elements in A.
!
!    Input, integer COL(N+1), points to the first element of 
!    each column.
!
!    Input, integer ROW(NZ_NUM), contains the row indices of 
!    the elements.
!
!    Input, real ( kind = rk ) A(NZ_NUM), the r8ccs matrix.
!
!    Input, integer ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: incx = 5
  integer n
  integer nz_num

  real ( kind = rk ) a(nz_num)
  real ( kind = rk ) aij
  integer col(n+1)
  character ( len = 14 ) ctemp(incx)
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
  integer k
  integer m
  integer row(nz_num)
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(''  Col:  '',5a14)' ) ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
!  1) Assume everything is zero.
!
      aij = 0.0D+00
      do j2 = 1, inc
        write ( ctemp(j2), '(f8.0,6x)' ) aij
      end do
!
!  2) Now consider each column J in J2LO to J2HI, 
!     and look at every nonzero, and check if it occurs in row I.
!
      do j = j2lo, j2hi
        do k = col(j), col(j+1)-1
          if ( row(k) == i ) then
            j2 = j - j2lo + 1
            write ( ctemp(j2), '(g14.6)' ) a(k)
          end if
        end do
      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8ccs_random ( m, n, nz_num, col, row, a, seed )

!*****************************************************************************80
!
!! r8ccs_RANDOM randomizes an r8ccs matrix.
!
!  Discussion:
!
!    The r8ccs format is the double precision sparse compressed column
!    format.  Associated with this format, we have an M by N matrix
!    with NZ_NUM nonzero entries.  We construct the column pointer
!    vector COL of length N+1, such that entries of column J will be
!    stored in positions COL(J) through COL(J+1)-1.  This indexing
!    refers to both the ROW and A vectors, which store the row indices
!    and the values of the nonzero entries.  The entries of the
!    ROW vector corresponding to each column are assumed to be
!    ascending sorted.
!
!    The r8ccs format is equivalent to the MATLAB "sparse" format,
!    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 October 2015
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    October 1992
!
!  Parameters:
!
!    Input, integer M, the number of rows of the matrix.
!
!    Input, integer N, the number of columns of the matrix.
!
!    Input, integer NZ_NUM, the number of nonzero elements in A.
!
!    Input, integer COL(N+1), points to the first element of each 
!    column.
!
!    Input, integer ROW(NZ_NUM), contains the row indices of the 
!    elements.
!
!    Output, real ( kind = rk ) A(NZ_NUM), the r8ccs matrix.
!
!    Input/output, integer SEED, a seed for the random number 
!    generator.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nz_num

  real ( kind = rk ) a(nz_num)
  integer col(n+1)
  integer m
  integer row(nz_num)
  integer seed

  call r8vec_uniform_01 ( nz_num, seed, a )

  return
end
subroutine r8ccs_read ( col_file, row_file, a_file, m, n, nz_num, col, row, a )

!*****************************************************************************80
!
!! r8ccs_READ reads an r8ccs matrix from three files.
!
!  Discussion:
!
!    This routine needs the values of M, N, and NZ_NUM, which can be 
!    determined by a call to r8ccs_READ_SIZE.
!
!    The r8ccs format is the double precision sparse compressed column
!    format.  Associated with this format, we have an M by N matrix
!    with NZ_NUM nonzero entries.  We construct the column pointer
!    vector COL of length N+1, such that entries of column J will be
!    stored in positions COL(J) through COL(J+1)-1.  This indexing
!    refers to both the ROW and A vectors, which store the row indices
!    and the values of the nonzero entries.  The entries of the
!    ROW vector corresponding to each column are assumed to be
!    ascending sorted.
!
!    The r8ccs format is equivalent to the MATLAB "sparse" format,
!    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    October 1992
!
!  Parameters:
!
!    Input, character ( len = * ) COL_FILE, ROW_FILE, A_FILE, the names of the 
!    files containing the column pointers, row indices, and matrix entries.
!
!    Input, integer M, N, the number of rows and columns in the 
!    matrix.
!
!    Input, integer NZ_NUM, the number of nonzero elements in 
!    the matrix.
!
!    Output, integer COL(N+1), the column pointers.
!
!    Output, integer ROW(NZ_NUM), the row indices.
!
!    Output, real ( kind = rk ) A(NZ_NUM), the nonzero elements 
!    of the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nz_num

  real ( kind = rk ) a(nz_num)
  character ( len = * ) a_file
  integer col(n+1)
  character ( len = * ) col_file
  integer input_unit
  integer ios
  integer k
  integer m
  integer row(nz_num)
  character ( len = * ) row_file

  call get_unit ( input_unit )
!
!  Read the column information.
!
  open ( unit = input_unit, file = col_file, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'r8ccs_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file "' &
      // trim ( col_file ) // '".'
    stop 1
  end if

  do k = 1, n + 1

    read ( input_unit, *, iostat = ios ) col(k)

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'r8ccs_READ - Fatal error!'
      write ( *, '(a,i8,a)' ) '  I/O error while reading record ', k, &
        ' of "' // trim ( col_file ) // '".'
      stop 1
    end if

  end do

  close ( unit = input_unit )
!
!  Read the row information.
!
  open ( unit = input_unit, file = row_file, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'r8ccs_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file "' &
      // trim ( row_file ) // '".'
    stop 1
  end if

  do k = 1, nz_num

    read ( input_unit, *, iostat = ios ) row(k)

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'r8ccs_READ - Fatal error!'
      write ( *, '(a,i8,a)' ) '  I/O error while reading record ', k, &
        ' of "' // trim ( row_file ) // '".'
      stop 1
    end if

  end do

  close ( unit = input_unit )
!
!  Read the value information.
!
  open ( unit = input_unit, file = a_file, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'r8ccs_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file "' &
      // trim ( a_file ) // '".'
    stop 1
  end if

  do k = 1, nz_num

    read ( input_unit, *, iostat = ios ) a(k)

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'r8ccs_READ - Fatal error!'
      write ( *, '(a,i8,a)' ) '  I/O error while reading record ', k, &
        ' of "' // trim ( a_file ) // '".'
      stop 1
    end if

  end do

  close ( unit = input_unit )

  return
end
subroutine r8ccs_read_size ( col_file, row_file, m, n, nz_num, base )

!*****************************************************************************80
!
!! r8ccs_READ_SIZE reads the sizes of an r8ccs sparse matrix from a file.
!
!  Discussion:
!
!    The value of M is "guessed" to be the largest value that occurs in
!    the ROW file.  However, if a row index of 0 is encountered, then
!    the value of M is incremented by 1.
!
!    The value of N is the number of records in the COL file minus 1.
!
!    The value of NZ_NUM is simply the number of records in the ROW file.
!
!    The value of BASE is 0 or 1, depending on whether the program
!    "guesses" that the row and column indices are 0-based or 1-based.
!    Although the first entry of the COL array might be used as evidence,
!    this program makes its determination based on whether it encounters
!    a 0 index in the ROW file.
!
!    The r8ccs format is the double precision sparse compressed column
!    format.  Associated with this format, we have an M by N matrix
!    with NZ_NUM nonzero entries.  We construct the column pointer
!    vector COL of length N+1, such that entries of column J will be
!    stored in positions COL(J) through COL(J+1)-1.  This indexing
!    refers to both the ROW and A vectors, which store the row indices
!    and the values of the nonzero entries.  The entries of the
!    ROW vector corresponding to each column are assumed to be
!    ascending sorted.
!
!    The r8ccs format is equivalent to the MATLAB "sparse" format,
!    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    28 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    October 1992
!
!  Parameters:
!
!    Input, character ( len = * ) COL_FILE, ROW_FILE, the names of the 
!    column and row files that describe the structure of the matrix.
!
!    Output, integer M, N, the inferred number of rows and columns 
!    in the sparse matrix.
!
!    Output, integer NZ_NUM, the number of nonzero entries in the
!    sparse matrix.
!
!    Output, integer BASE, is 0 if the row indexing is believed
!    to be 0-based, and 1 if the row-index is believed to be
!    1-based.  In uncertain cases, BASE = 1 is the default.
!
  implicit none

  integer base
  integer col
  character ( len = * ) col_file
  integer input_unit
  integer ios
  integer m
  integer n
  integer nz_num
  integer row
  character ( len = * ) row_file
!
!  Default values.
!
  m = -1
  n = -1
  nz_num = -1
  base = -1
!
!  Check the COL file first.
!
  call get_unit ( input_unit )

  open ( unit = input_unit, file = col_file, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'r8ccs_READ_SIZE - Fatal error!'
    write ( *, '(a)' ) '  Could not open  the input file "' &
      // trim ( col_file ) // '".'
    stop 1
  end if

  n = -1

  do

    read ( input_unit, *, iostat = ios ) col

    if ( ios /= 0 ) then
      exit
    end if

    n = n + 1

  end do

  close ( unit = input_unit )
!
!  Check the ROW file.
!
  call get_unit ( input_unit )

  open ( unit = input_unit, file = row_file, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'r8ccs_READ_SIZE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file "' &
      // trim ( row_file ) // '".'
    stop 1
  end if

  base = 1
  m = 0
  nz_num = 0
  
  do

    read ( input_unit, *, iostat = ios ) row

    if ( ios /= 0 ) then
      exit
    end if

    nz_num = nz_num + 1
    m = max ( m, row )
    if ( row == 0 ) then
      base = 0
    end if

  end do

  close ( unit = input_unit )

  return
end
subroutine r8ccs_set ( m, n, nz_num, col, row, a, i, j, aij )

!*****************************************************************************80
!
!! r8ccs_SET sets a value of an r8ccs matrix.
!
!  Discussion:
!
!    The r8ccs format is the double precision sparse compressed column
!    format.  Associated with this format, we have an M by N matrix
!    with NZ_NUM nonzero entries.  We construct the column pointer
!    vector COL of length N+1, such that entries of column J will be
!    stored in positions COL(J) through COL(J+1)-1.  This indexing
!    refers to both the ROW and A vectors, which store the row indices
!    and the values of the nonzero entries.  The entries of the
!    ROW vector corresponding to each column are assumed to be
!    ascending sorted.
!
!    The r8ccs format is equivalent to the MATLAB "sparse" format,
!    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    October 1992
!
!  Parameters:
!
!    Input, integer M, the number of rows of the matrix.
!
!    Input, integer N, the number of columns of the matrix.
!
!    Input, integer NZ_NUM, the number of nonzero entries.
!
!    Input, integer COL(N+1), indicate where each column's data 
!    begins.
!
!    Input, integer ROW(NZ_NUM), the row indices.
!
!    Input/output, real ( kind = rk ) A(NZ_NUM), the nonzero entries.
!    On output, the entry of A corresponding to (I,J) has been reset.
!
!    Input, integer I, J, the indices of the value to retrieve.
!
!    Input, real ( kind = rk ) AIJ, the new value of A(I,J).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nz_num

  real ( kind = rk ) a(nz_num)
  real ( kind = rk ) aij
  integer col(n+1)
  integer i
  integer j
  integer k
  integer m
  integer row(nz_num)
!
!  Seek sparse index K corresponding to full index (I,J).
!
  call r8ccs_ijk ( m, n, nz_num, col, row, i, j, k )
!
!  If no K was found, we fail.
!
  if ( k == -1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'r8ccs_SET - Fatal error!'
    write ( *, '(a)' ) '  r8ccs_IJK could not find the entry.'
    write ( *, '(a,i8)' ) '  Row I = ', i
    write ( *, '(a,i8)' ) '  Col J = ', j
    stop 1
  end if

  a(k) = aij

  return
end
subroutine r8ccs_to_r8ge ( m, n, nz_num, col, row, a, b )

!*****************************************************************************80
!
!! r8ccs_TO_R8GE converts an r8ccs matrix to an R8GE matrix.
!
!  Discussion:
!
!    The r8ccs format is the double precision sparse compressed column
!    format.  Associated with this format, we have an M by N matrix
!    with NZ_NUM nonzero entries.  We construct the column pointer
!    vector COL of length N+1, such that entries of column J will be
!    stored in positions COL(J) through COL(J+1)-1.  This indexing
!    refers to both the ROW and A vectors, which store the row indices
!    and the values of the nonzero entries.  The entries of the
!    ROW vector corresponding to each column are assumed to be
!    ascending sorted.
!
!    The r8ccs format is equivalent to the MATLAB "sparse" format,
!    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    October 1992
!
!  Parameters:
!
!    Input, integer M, the number of rows of the matrix.
!
!    Input, integer N, the number of columns of the matrix.
!
!    Input, integer NZ_NUM, the number of nonzero elements in A.
!
!    Input, integer COL(N+1), points to the first element of 
!    each column.
!
!    Input, integer ROW(NZ_NUM), contains the row indices 
!    of the elements.
!
!    Input, real ( kind = rk ) A(NZ_NUM), the r8ccs matrix.
!
!    Output, real ( kind = rk ) B(M,N), the R8GE matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n
  integer nz_num

  real ( kind = rk ) a(nz_num)
  real ( kind = rk ) b(m,n)
  integer col(n+1)
  integer j
  integer k
  integer row(nz_num)

  b(1:m,1:n) = 0.0D+00

  do j = 1, n
    do k = col(j), col(j+1)-1
      b(row(k),j) = a(k)
    end do
  end do

  return
end
subroutine r8ccs_write ( col_file, row_file, a_file, m, n, nz_num, col, row, a )

!*****************************************************************************80
!
!! r8ccs_WRITE writes an r8ccs matrix to three files.
!
!  Discussion:
!
!    The r8ccs format is the double precision sparse compressed column
!    format.  Associated with this format, we have an M by N matrix
!    with NZ_NUM nonzero entries.  We construct the column pointer
!    vector COL of length N+1, such that entries of column J will be
!    stored in positions COL(J) through COL(J+1)-1.  This indexing
!    refers to both the ROW and A vectors, which store the row indices
!    and the values of the nonzero entries.  The entries of the
!    ROW vector corresponding to each column are assumed to be
!    ascending sorted.
!
!    The r8ccs format is equivalent to the MATLAB "sparse" format,
!    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    October 1992
!
!  Parameters:
!
!    Input, character ( len = * ) COL_FILE, ROW_FILE, A_FILE, the names of the 
!    files containing the column pointers, row entries, and matrix entries.
!
!    Input, integer M, N, the number of rows and columns 
!    in the matrix.
!
!    Input, integer NZ_NUM, the number of nonzero elements 
!    in the matrix.
!
!    Input, integer COL(N+1), the column pointers.
!
!    Input, integer ROW(NZ_NUM), the row indices.
!
!    Input, real ( kind = rk ) A(NZ_NUM), the nonzero elements 
!    of the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nz_num

  real ( kind = rk ) a(nz_num)
  character ( len = * ) a_file
  integer col(n+1)
  character ( len = * ) col_file
  integer ios
  integer k
  integer m
  integer output_unit
  integer row(nz_num)
  character ( len = * ) row_file

  call get_unit ( output_unit )
!
!  Write the column information.
!
  open ( unit = output_unit, file = col_file, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'r8ccs_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file "' &
      // trim ( col_file ) // '".'
    stop 1
  end if

  do k = 1, n + 1

    write ( output_unit, '(i8)' ) col(k)

  end do

  close ( unit = output_unit )
!
!  Write the row information.
!
  open ( unit = output_unit, file = row_file, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'r8ccs_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file "' &
      // trim ( row_file ) // '".'
    stop 1
  end if

  do k = 1, nz_num

    write ( output_unit, '(i8)' ) row(k)

  end do

  close ( unit = output_unit )
!
!  Write the value information.
!
  open ( unit = output_unit, file = a_file, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'r8ccs_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file "' &
      // trim ( a_file ) // '".'
    stop 1
  end if

  do k = 1, nz_num

    write ( output_unit, '(g14.6)' ) a(k)

  end do

  close ( unit = output_unit )

  return
end
subroutine r8ccs_zeros ( m, n, nz_num, col, row, a )

!*****************************************************************************80
!
!! r8ccs_ZEROS zeroes an r8ccs matrix.
!
!  Discussion:
!
!    The r8ccs format is the double precision sparse compressed column
!    format.  Associated with this format, we have an M by N matrix
!    with NZ_NUM nonzero entries.  We construct the column pointer
!    vector COL of length N+1, such that entries of column J will be
!    stored in positions COL(J) through COL(J+1)-1.  This indexing
!    refers to both the ROW and A vectors, which store the row indices
!    and the values of the nonzero entries.  The entries of the
!    ROW vector corresponding to each column are assumed to be
!    ascending sorted.
!
!    The r8ccs format is equivalent to the MATLAB "sparse" format,
!    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 January 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    October 1992
!
!  Parameters:
!
!    Input, integer M, the number of rows of the matrix.
!
!    Input, integer N, the number of columns of the matrix.
!
!    Input, integer NZ_NUM, the number of nonzero elements in A.
!
!    Input, integer COL(N+1), points to the first element of each 
!    column.
!
!    Input, integer ROW(NZ_NUM), contains the row indices of the 
!    elements.
!
!    Output, real ( kind = rk ) A(NZ_NUM), the r8ccs matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nz_num

  real ( kind = rk ) a(nz_num)
  integer col(n+1)
  integer m
  integer row(nz_num)

  a(1:nz_num) = 0.0D+00

  return
end
subroutine r8ge_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8GE_PRINT prints an R8GE matrix.
!
!  Discussion:
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!    R8GE storage is used by LINPACK and LAPACK.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 May 2000
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
!    Input, real ( kind = rk ) A(M,N), the R8GE matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  character ( len = * ) title

  call r8ge_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8ge_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8GE_PRINT_SOME prints some of an R8GE matrix.
!
!  Discussion:
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!    R8GE storage is used by LINPACK and LAPACK.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 March 2001
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
!    Input, real ( kind = rk ) A(M,N), the R8GE matrix.
!
!    Input, integer ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
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
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(''  Col:  '',5a14)' ) ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        write ( ctemp(j2), '(g14.6)' ) a(i,j)

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8vec_indicator1 ( n, a )

!*****************************************************************************80
!
!! R8VEC_INDICATOR1 sets an R8VEC to the indicator1 vector.
!
!  Discussion:
!
!    A(1:N) = (/ 1 : N /)
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements of A.
!
!    Output, real ( kind = rk ) A(N), the array to be initialized.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  integer i

  do i = 1, n
    a(i) = real ( i, kind = rk )
  end do

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
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
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  integer i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i8,g14.6)' ) i, a(i)
  end do

  return
end
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of real ( kind = rk ) values.
!
!    For now, the input quantity SEED is an integer variable.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer N, the number of entries 
!    in the vector.
!
!    Input/output, integer SEED, the "seed" value, 
!    which should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = rk ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer k
  integer seed
  real ( kind = rk ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + 2147483647
    end if

    r(i) = real ( seed, kind = rk ) * 4.656612875D-10

  end do

  return
end

