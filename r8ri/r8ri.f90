function i4_log_10 ( i )

!*****************************************************************************80
!
!! i4_log_10() returns the integer part of the logarithm base 10 of an I4.
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
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an I4VEC.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    28 November 2000
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
!    Input, character ( len = * ) TITLE, a title first.
!    TITLE may be blank.
!
  implicit none

  integer n

  integer a(n)
  integer big
  integer i
  character ( len = * ) title

  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  big = maxval ( abs ( a(1:n) ) )

  write ( *, '(a)' ) ' '
  if ( big < 1000 ) then
    do i = 1, n
      write ( *, '(i8,1x,i4)' ) i, a(i)
    end do
  else if ( big < 1000000 ) then
    do i = 1, n
      write ( *, '(i8,1x,i7)' ) i, a(i)
    end do
  else
    do i = 1, n
      write ( *, '(i8,i11)' ) i, a(i)
    end do
  end if

  return
end
function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = rk ) value.
!
!    For now, the input quantity SEED is an integer variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
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
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
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
!    Input/output, integer SEED, the "seed" value, 
!    which should NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = rk ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer k
  real ( kind = rk ) r8_uniform_01
  integer seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = rk ) * 4.656612875D-10

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
subroutine r8ge_to_r8ri ( n, a, nz, ija, a_r8ri )

!*****************************************************************************80
!
!! R8GE_TO_R8RI converts an R8GE matrix to R8RI form.
!
!  Discussion:
!
!    A R8GE matrix is in general storage.
!
!    An R8RI matrix is in row indexed sparse storage form, using an index
!    array IJA and a value array A.  The first N entries of A store the
!    diagonal elements in order.  The first N entries of IJA store the index
!    of the first off-diagonal element of the corresponding row; if there is
!    no off-diagonal element in that row, it is one greater than the index
!    in A of the most recently stored element in the previous row.
!    Location 1 of IJA is always equal to N+2; location N+1 of IJA is one
!    greater than the index in A of the last off-diagonal element of the
!    last row.  Location N+1 of A is not used.  Entries in A with index
!    N+2 or greater contain the off-diagonal values, ordered by row, and
!    then by column.  Entries in IJA with index N+2 or greater contain the
!    column number of the corresponding element in A.
!
!  Example:
!
!    A:
!      3 0 1 0 0
!      0 4 0 0 0
!      0 7 5 9 0
!      0 0 0 0 2
!      0 0 0 6 8
!
!    NZ = 11
!
!    IJA:
!      7  8  8 10 11 12  3  2  4  5  4
!
!    A:
!      3  4  5  0  8  *  1  7  9  2  6
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 January 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
!    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
!    Third Edition,
!    Cambridge University Press, 2007,
!    ISBN13: 978-0-521-88068-8,
!    LC: QA297.N866.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real ( kind = rk ) A(N,N), the matrix stored in GE 
!    or "general" format.
!
!    Input, integer NZ, the size required for the RI
!    or "row indexed" sparse storage.
!
!    Output, integer IJA(NZ), the index vector.
!
!    Output, real ( kind = rk ) A_R8RI(NZ), the value vector.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nz

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) a_r8ri(nz)
  integer i
  integer ija(nz)
  integer im
  integer j
  integer k

  do k = 1, n
    i = k
    j = k
    a_r8ri(k) = a(i,j)
  end do

  k = n + 1
  a_r8ri(k) = 0.0D+00

  ija(1:n+1) = 0
  im = 1

  do i = 1, n
    do j = 1, n

      if ( i == j ) then
        cycle
      end if

      if ( a(i,j) == 0.0D+00 ) then
        cycle
      end if

      k = k + 1

      if ( ija(i) == 0 ) then
        ija(im:i) = k
        im = i + 1
      end if

      ija(k) = j
      a_r8ri(k) = a(i,j)

    end do

  end do

  ija(n+1) = k + 1

  return
end
subroutine r8ge_to_r8ri_size ( n, a, nz )

!*****************************************************************************80
!
!! R8GE_TO_R8RI_SIZE determines the size of an R8RI matrix.
!
!  Discussion:
!
!    An R8RI matrix is in row indexed sparse storage form, using an index
!    array IJA and a value array A.  The first N entries of A store the
!    diagonal elements in order.  The first N entries of IJA store the index
!    of the first off-diagonal element of the corresponding row; if there is
!    no off-diagonal element in that row, it is one greater than the index
!    in A of the most recently stored element in the previous row.
!    Location 1 of IJA is always equal to N+2; location N+1 of IJA is one
!    greater than the index in A of the last off-diagonal element of the
!    last row.  Location N+1 of A is not used.  Entries in A with index
!    N+2 or greater contain the off-diagonal values, ordered by row, and
!    then by column.  Entries in IJA with index N+2 or greater contain the
!    column number of the corresponding element in A.
!
!  Example:
!
!    A:
!      3 0 1 0 0
!      0 4 0 0 0
!      0 7 5 9 0
!      0 0 0 0 2
!      0 0 0 6 8
!
!    NZ = 11
!
!    IJA:
!      7  8  8 10 11 12  3  2  4  5  4
!
!    A:
!      3  4  5  0  8  *  1  7  9  2  6
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 January 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
!    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
!    Third Edition,
!    Cambridge University Press, 2007,
!    ISBN13: 978-0-521-88068-8,
!    LC: QA297.N866.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real ( kind = rk ) A(N,N), the matrix stored in GE 
!    or "general" format.
!
!    Output, integer NZ, the size required for the RI
!    or "row indexed" sparse storage.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,n)
  integer i
  integer j
  integer nz

  nz = n + 1

  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        cycle
      end if
      if ( a(i,j) /= 0.0D+00 ) then
        nz = nz + 1
      end if
    end do
  end do

  return
end
subroutine r8ri_dif2 ( n, nz, ija, a )

!*****************************************************************************80
!
!! R8RI_DIF2 stores the second difference matrix in R8RI format.
!
!  Discussion:
!
!    An R8RI matrix is in row indexed sparse storage form, using an index
!    array IJA and a value array A.  The first N entries of A store the
!    diagonal elements in order.  The first N entries of IJA store the index
!    of the first off-diagonal element of the corresponding row; if there is
!    no off-diagonal element in that row, it is one greater than the index
!    in A of the most recently stored element in the previous row.
!    Location 1 of IJA is always equal to N+2; location N+1 of IJA is one
!    greater than the index in A of the last off-diagonal element of the
!    last row.  Location N+1 of A is not used.  Entries in A with index
!    N+2 or greater contain the off-diagonal values, ordered by row, and
!    then by column.  Entries in IJA with index N+2 or greater contain the
!    column number of the corresponding element in A.
!
!  Example:
!
!    A:
!      3 0 1 0 0
!      0 4 0 0 0
!      0 7 5 9 0
!      0 0 0 0 2
!      0 0 0 6 8
!
!    NZ = 11
!
!    IJA:
!      7  8  8 10 11 12  3  2  4  5  4
!
!    A:
!      3  4  5  0  8  *  1  7  9  2  6
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 July 2016
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
!    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
!    Third Edition,
!    Cambridge University Press, 2007,
!    ISBN13: 978-0-521-88068-8,
!    LC: QA297.N866.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer NZ, the size required for the RI
!    or "row indexed" sparse storage.  NZ = 3*N-1.
!
!    Output, integer IJA(NZ), the index vector.
!
!    Output, real ( kind = rk ) A(NZ), the value vector.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nz

  real ( kind = rk ) a(nz)
  integer i
  integer ija(nz)
  integer k
!
!  Diagonal elements of A.
!
  do i = 1, n
    a(i) = 2.0D+00
  end do
!
!  First N entries of IJA store first offdiagonal of each row.
!
  k = n + 2

  do i = 1, n
    ija(i) = k
    if ( i == 1 .or. i == n ) then
      k = k + 1
    else
      k = k + 2
    end if
  end do
!
!  IJA(N+1) stores one beyond last element of A.
!
  ija(n+1) = k
  a(n+1) = 0.0D+00
!
!  IJA(N+2), A(N+2) and beyond store column and value.
!
  k = n + 1

  do i = 1, n

    if ( i == 1 ) then
      k = k + 1
      ija(k) = i + 1
      a(k) = - 1.0D+00
    else if ( i < n ) then
      k = k + 1
      ija(k) = i - 1
      a(k) = - 1.0D+00
      k = k + 1
      ija(k) = i + 1
      a(k) = - 1.0D+00
    else if ( i == n ) then
      k = k + 1
      ija(k) = i - 1
      a(k) = - 1.0D+00
    end if

  end do

  return
end
subroutine r8ri_indicator ( n, nz, ija, a )

!*****************************************************************************80
!
!! R8RI_INDICATOR returns the R8RI indicator matrix for given sparsity.
!
!  Discussion:
!
!    An R8RI matrix is in row indexed sparse storage form, using an index
!    array IJA and a value array A.  The first N entries of A store the
!    diagonal elements in order.  The first N entries of IJA store the index
!    of the first off-diagonal element of the corresponding row; if there is
!    no off-diagonal element in that row, it is one greater than the index
!    in A of the most recently stored element in the previous row.
!    Location 1 of IJA is always equal to N+2; location N+1 of IJA is one
!    greater than the index in A of the last off-diagonal element of the
!    last row.  Location N+1 of A is not used.  Entries in A with index
!    N+2 or greater contain the off-diagonal values, ordered by row, and
!    then by column.  Entries in IJA with index N+2 or greater contain the
!    column number of the corresponding element in A.
!
!  Example:
!
!    A:
!      3 0 1 0 0
!      0 4 0 0 0
!      0 7 5 9 0
!      0 0 0 0 2
!      0 0 0 6 8
!
!    NZ = 11
!
!    IJA:
!      7  8  8 10 11 12  3  2  4  5  4
!
!    A:
!      3  4  5  0  8  *  1  7  9  2  6
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 July 2016
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
!    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
!    Third Edition,
!    Cambridge University Press, 2007,
!    ISBN13: 978-0-521-88068-8,
!    LC: QA297.N866.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer NZ, the size required for the RI
!    or "row indexed" sparse storage.  NZ = 3*N-1.
!
!    Input, integer IJA(NZ), the index vector.
!
!    Output, real ( kind = rk ) A(NZ), the value vector.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nz

  real ( kind = rk ) a(nz)
  integer fac
  integer i
  integer i4_log_10
  integer ija(nz)
  integer j
  integer k

  fac = 10 ** ( i4_log_10 ( n ) + 1 )
!
!  Diagonal elements of A.
!
  do i = 1, n
    a(i) = real ( fac * i + i, kind = rk )
  end do

  do i = 1, n
    do k = ija(i), ija(i+1) - 1
      j = ija(k)
      a(k) = real ( fac * i + j, kind = rk )
    end do
  end do

  return
end
subroutine r8ri_mtv ( n, nz, ija, a, x, b )

!*****************************************************************************80
!
!! R8RI_MTV multiplies the transpose of an R8RI matrix times a vector.
!
!  Discussion:
!
!    An R8RI matrix is in row indexed sparse storage form, using an index
!    array IJA and a value array A.  The first N entries of A store the
!    diagonal elements in order.  The first N entries of IJA store the index
!    of the first off-diagonal element of the corresponding row; if there is
!    no off-diagonal element in that row, it is one greater than the index
!    in A of the most recently stored element in the previous row.
!    Location 1 of IJA is always equal to N+2; location N+1 of IJA is one
!    greater than the index in A of the last off-diagonal element of the
!    last row.  Location N+1 of A is not used.  Entries in A with index
!    N+2 or greater contain the off-diagonal values, ordered by row, and
!    then by column.  Entries in IJA with index N+2 or greater contain the
!    column number of the corresponding element in A.
!
!  Example:
!
!    A:
!      3 0 1 0 0
!      0 4 0 0 0
!      0 7 5 9 0
!      0 0 0 0 2
!      0 0 0 6 8
!
!    NZ = 11
!
!    IJA:
!      7  8  8 10 11 12  3  2  4  5  4
!
!    A:
!      3  4  5  0  8  *  1  7  9  2  6
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
!    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
!    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
!    Third Edition,
!    Cambridge University Press, 2007,
!    ISBN13: 978-0-521-88068-8,
!    LC: QA297.N866.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer NZ, the size required for the RI
!    or "row indexed" sparse storage.
!
!    Input, integer IJA(NZ), the index vector.
!
!    Input, real ( kind = rk ) A(NZ), the value vector.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied.
!
!    Output, real ( kind = rk ) B(N), the product A'*X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nz

  real ( kind = rk ) a(nz)
  real ( kind = rk ) b(n)
  integer i
  integer ija(nz)
  integer j
  integer k
  real ( kind = rk ) x(n)

  if ( ija(1) /= n + 2 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'R8RI_MTV - Fatal error!'
    write ( *, '(a)' ) '  The values IJA(1) and N are inconsistent.'
    stop 1
  end if

  do i = 1, n
    b(i) = a(i) * x(i)
  end do

  do i = 1, n
    do k = ija(i), ija(i+1) - 1
      j = ija(k)
      b(j) = b(j) + a(k) * x(i)
    end do
  end do

  return
end
subroutine r8ri_mv ( n, nz, ija, a, x, b )

!*****************************************************************************80
!
!! R8RI_MV multiplies an R8RI matrix times a vector.
!
!  Discussion:
!
!    An R8RI matrix is in row indexed sparse storage form, using an index
!    array IJA and a value array A.  The first N entries of A store the
!    diagonal elements in order.  The first N entries of IJA store the index
!    of the first off-diagonal element of the corresponding row; if there is
!    no off-diagonal element in that row, it is one greater than the index
!    in A of the most recently stored element in the previous row.
!    Location 1 of IJA is always equal to N+2; location N+1 of IJA is one
!    greater than the index in A of the last off-diagonal element of the
!    last row.  Location N+1 of A is not used.  Entries in A with index
!    N+2 or greater contain the off-diagonal values, ordered by row, and
!    then by column.  Entries in IJA with index N+2 or greater contain the
!    column number of the corresponding element in A.
!
!  Example:
!
!    A:
!      3 0 1 0 0
!      0 4 0 0 0
!      0 7 5 9 0
!      0 0 0 0 2
!      0 0 0 6 8
!
!    NZ = 11
!
!    IJA:
!      7  8  8 10 11 12  3  2  4  5  4
!
!    A:
!      3  4  5  0  8  *  1  7  9  2  6
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
!    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
!    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
!    Third Edition,
!    Cambridge University Press, 2007,
!    ISBN13: 978-0-521-88068-8,
!    LC: QA297.N866.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer NZ, the size required for the RI
!    or "row indexed" sparse storage.
!
!    Input, integer IJA(NZ), the index vector.
!
!    Input, real ( kind = rk ) A(NZ), the value vector.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied.
!
!    Output, real ( kind = rk ) B(N), the product A*X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nz

  real ( kind = rk ) a(nz)
  real ( kind = rk ) b(n)
  integer i
  integer ija(nz)
  integer k
  real ( kind = rk ) x(n)

  if ( ija(1) /= n + 2 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'R8RI_MV - Fatal error!'
    write ( *, '(a)' ) '  The values IJA(1) and N are inconsistent.'
    stop 1
  end if

  do i = 1, n
    b(i) = a(i) * x(i)
    do k = ija(i), ija(i+1) - 1
      b(i) = b(i) + a(k) * x(ija(k))
    end do
  end do

  return
end
subroutine r8ri_print ( n, nz, ija, a, title )

!*****************************************************************************80
!
!! R8RI_PRINT prints an R8RI matrix.
!
!  Discussion:
!
!    An R8RI matrix is in row indexed sparse storage form, using an index
!    array IJA and a value array A.  The first N entries of A store the
!    diagonal elements in order.  The first N entries of IJA store the index
!    of the first off-diagonal element of the corresponding row; if there is
!    no off-diagonal element in that row, it is one greater than the index
!    in A of the most recently stored element in the previous row.
!    Location 1 of IJA is always equal to N+2; location N+1 of IJA is one
!    greater than the index in A of the last off-diagonal element of the
!    last row.  Location N+1 of A is not used.  Entries in A with index
!    N+2 or greater contain the off-diagonal values, ordered by row, and
!    then by column.  Entries in IJA with index N+2 or greater contain the
!    column number of the corresponding element in A.
!
!  Example:
!
!    A:
!      3 0 1 0 0
!      0 4 0 0 0
!      0 7 5 9 0
!      0 0 0 0 2
!      0 0 0 6 8
!
!    NZ = 11
!
!    IJA:
!      7  8  8 10 11 12  3  2  4  5  4
!
!    A:
!      3  4  5  0  8  *  1  7  9  2  6
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 July 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer NZ, the size required for the RI
!    or "row indexed" sparse storage.
!
!    Input, integer IJA(NZ), the index vector.
!
!    Input, real ( kind = rk ) A(NZ), the value vector.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nz

  real ( kind = rk ) a(nz)
  integer ija(nz)
  integer n
  character ( len = * ) title

  call r8ri_print_some ( n, nz, ija, a, 1, 1, n, n, title )

  return
end
subroutine r8ri_print_some ( n, nz, ija, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8RI_PRINT_SOME prints some of an R8RI matrix.
!
!  Discussion:
!
!    An R8RI matrix is in row indexed sparse storage form, using an index
!    array IJA and a value array A.  The first N entries of A store the
!    diagonal elements in order.  The first N entries of IJA store the index
!    of the first off-diagonal element of the corresponding row; if there is
!    no off-diagonal element in that row, it is one greater than the index
!    in A of the most recently stored element in the previous row.
!    Location 1 of IJA is always equal to N+2; location N+1 of IJA is one
!    greater than the index in A of the last off-diagonal element of the
!    last row.  Location N+1 of A is not used.  Entries in A with index
!    N+2 or greater contain the off-diagonal values, ordered by row, and
!    then by column.  Entries in IJA with index N+2 or greater contain the
!    column number of the corresponding element in A.
!
!  Example:
!
!    A:
!      3 0 1 0 0
!      0 4 0 0 0
!      0 7 5 9 0
!      0 0 0 0 2
!      0 0 0 6 8
!
!    NZ = 11
!
!    IJA:
!      7  8  8 10 11 12  3  2  4  5  4
!
!    A:
!      3  4  5  0  8  *  1  7  9  2  6
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 July 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer NZ, the size required for the RI
!    or "row indexed" sparse storage.
!
!    Input, integer IJA(NZ), the index vector.
!
!    Input, real ( kind = rk ) A(NZ), the value vector.
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
  integer nz

  real ( kind = rk ) a(nz)
  real ( kind = rk ) aij
  character ( len = 14 ) ctemp(incx)
  integer i
  integer i2hi
  integer i2lo
  integer ihi
  integer ija(nz)
  integer ilo
  integer inc
  integer j
  integer j2
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  integer k
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
    i2hi = min ( ihi, n )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, J2LO <= J <= J2HI.
!
!  1) Assume everything is zero.
!
      aij = 0.0D+00
      do j2 = 1, inc
        write ( ctemp(j2), '(f8.0,6x)' ) aij
      end do
!
!  2) Insert the diagonal entry, if appropriate.
!
      if ( j2lo <= i .and. i <= j2hi ) then
        j2 = 1 + ( i - j2lo )
        aij = a(i)
        write ( ctemp(j2), '(g14.6)' ) aij
      end if
!
!  3) Now examine all the offdiagonal entries.
!
      do k = ija(i), ija(i+1) - 1
        j = ija(k)
        if ( j2lo <= j .and. j <= j2hi ) then 
          j2 = 1 + ( j - j2lo )
          aij = a(k)
          write ( ctemp(j2), '(g14.6)' ) aij
        end if
      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8ri_random ( n, nz, ija, seed, a )

!*****************************************************************************80
!
!! R8RI_RANDOM randomizes an R8RI matrix for given sparsity.
!
!  Discussion:
!
!    An R8RI matrix is in row indexed sparse storage form, using an index
!    array IJA and a value array A.  The first N entries of A store the
!    diagonal elements in order.  The first N entries of IJA store the index
!    of the first off-diagonal element of the corresponding row; if there is
!    no off-diagonal element in that row, it is one greater than the index
!    in A of the most recently stored element in the previous row.
!    Location 1 of IJA is always equal to N+2; location N+1 of IJA is one
!    greater than the index in A of the last off-diagonal element of the
!    last row.  Location N+1 of A is not used.  Entries in A with index
!    N+2 or greater contain the off-diagonal values, ordered by row, and
!    then by column.  Entries in IJA with index N+2 or greater contain the
!    column number of the corresponding element in A.
!
!  Example:
!
!    A:
!      3 0 1 0 0
!      0 4 0 0 0
!      0 7 5 9 0
!      0 0 0 0 2
!      0 0 0 6 8
!
!    NZ = 11
!
!    IJA:
!      7  8  8 10 11 12  3  2  4  5  4
!
!    A:
!      3  4  5  0  8  *  1  7  9  2  6
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 July 2016
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
!    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
!    Third Edition,
!    Cambridge University Press, 2007,
!    ISBN13: 978-0-521-88068-8,
!    LC: QA297.N866.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer NZ, the size required for the RI
!    or "row indexed" sparse storage.  NZ = 3*N-1.
!
!    Input, integer IJA(NZ), the index vector.
!
!    Input/output, integer SEED, a seed for the random number
!    generator.
!
!    Output, real ( kind = rk ) A(NZ), the value vector.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nz

  real ( kind = rk ) a(nz)
  integer i
  integer ija(nz)
  integer j
  integer k
  real ( kind = rk ) r8_uniform_01
  integer seed
!
!  Diagonal elements of A.
!
  do i = 1, n
    a(i) = r8_uniform_01 ( seed )
  end do

  do i = 1, n
    do k = ija(i), ija(i+1) - 1
      j = ija(k)
      a(k) = r8_uniform_01 ( seed )
    end do
  end do

  return
end
subroutine r8ri_to_r8ge ( n, nz, ija, a, a_r8ge )

!*****************************************************************************80
!
!! R8RI_TO_R8GE converts an R8RI matrix to R8GE form.
!
!  Discussion:
!
!    An R8RI matrix is in row indexed sparse storage form, using an index
!    array IJA and a value array A.  The first N entries of A store the
!    diagonal elements in order.  The first N entries of IJA store the index
!    of the first off-diagonal element of the corresponding row; if there is
!    no off-diagonal element in that row, it is one greater than the index
!    in A of the most recently stored element in the previous row.
!    Location 1 of IJA is always equal to N+2; location N+1 of IJA is one
!    greater than the index in A of the last off-diagonal element of the
!    last row.  Location N+1 of A is not used.  Entries in A with index
!    N+2 or greater contain the off-diagonal values, ordered by row, and
!    then by column.  Entries in IJA with index N+2 or greater contain the
!    column number of the corresponding element in A.
!
!    A R8GE matrix is in general storage.
!
!  Example:
!
!    A:
!      3 0 1 0 0
!      0 4 0 0 0
!      0 7 5 9 0
!      0 0 0 0 2
!      0 0 0 6 8
!
!    NZ = 11
!
!    IJA:
!      7  8  8 10 11 12  3  2  4  5  4
!
!    A:
!      3  4  5  0  8  *  1  7  9  2  6
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 January 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
!    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
!    Third Edition,
!    Cambridge University Press, 2007,
!    ISBN13: 978-0-521-88068-8,
!    LC: QA297.N866.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer NZ, the size required for the RI
!    or "row indexed" sparse storage.
!
!    Input, integer IJA(NZ), the index vector.
!
!    Input, real ( kind = rk ) A(NZ), the value vector.
!
!    Output, real ( kind = rk ) A(N,N), the matrix stored in GE 
!    or "general" format.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nz

  real ( kind = rk ) a(nz)
  real ( kind = rk ) a_r8ge(n,n)
  integer i
  integer ija(nz)
  integer j
  integer k

  a_r8ge(1:n,1:n) = 0.0D+00

  do k = 1, n
    i = k
    j = k
    a_r8ge(i,j) = a(k)
  end do

  do i = 1, n
    do k = ija(i), ija(i+1) - 1
      j = ija(k)
      a_r8ge(i,j) = a(k)
    end do
  end do

  return
end
subroutine r8ri_zeros ( n, nz, ija, a )

!*****************************************************************************80
!
!! R8RI_ZEROS zeroes an R8RI matrix.
!
!  Discussion:
!
!    An R8RI matrix is in row indexed sparse storage form, using an index
!    array IJA and a value array A.  The first N entries of A store the
!    diagonal elements in order.  The first N entries of IJA store the index
!    of the first off-diagonal element of the corresponding row; if there is
!    no off-diagonal element in that row, it is one greater than the index
!    in A of the most recently stored element in the previous row.
!    Location 1 of IJA is always equal to N+2; location N+1 of IJA is one
!    greater than the index in A of the last off-diagonal element of the
!    last row.  Location N+1 of A is not used.  Entries in A with index
!    N+2 or greater contain the off-diagonal values, ordered by row, and
!    then by column.  Entries in IJA with index N+2 or greater contain the
!    column number of the corresponding element in A.
!
!  Example:
!
!    A:
!      3 0 1 0 0
!      0 4 0 0 0
!      0 7 5 9 0
!      0 0 0 0 2
!      0 0 0 6 8
!
!    NZ = 11
!
!    IJA:
!      7  8  8 10 11 12  3  2  4  5  4
!
!    A:
!      3  4  5  0  8  *  1  7  9  2  6
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
!    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
!    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
!    Third Edition,
!    Cambridge University Press, 2007,
!    ISBN13: 978-0-521-88068-8,
!    LC: QA297.N866.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer NZ, the size required for the RI
!    or "row indexed" sparse storage.
!
!    Input, integer IJA(NZ), the index vector.
!
!    Output, real ( kind = rk ) A(NZ), the value vector.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nz

  real ( kind = rk ) a(nz)
  integer n
  integer ija(nz)

  a(1:nz) = 0.0D+00

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

