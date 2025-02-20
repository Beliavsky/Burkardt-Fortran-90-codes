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
subroutine r8ncf_dif2 ( m, n, nz_num, rowcol, a )

!*****************************************************************************80
!
!! R8NCF_DIF2 sets up an R8NCF second difference matrix.
!
!  Discussion:
!
!    The R8NCF storage format stores NZ_NUM, the number of nonzeros, 
!    a real array containing the nonzero values, a 2 by NZ_NUM integer 
!    array storing the row and column of each nonzero entry.
!
!    The R8NCF format is used by NSPCG.  NSPCG requires that the information
!    for the diagonal entries of the matrix must come first.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 July 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns in
!    the matrix.
!
!    Input, integer NZ_NUM, the number of nonzero entries.
!
!    Input, integer ROWCOL(2,NZ_NUM), the coordinates of 
!    the nonzero entries.
!
!    Output, real ( kind = rk ) A(NZ_NUM), the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nz_num

  real ( kind = rk ) a(nz_num)
  integer i
  integer j
  integer k
  integer m
  integer n
  integer rowcol(2,nz_num)

  do k = 1, nz_num
    i = rowcol(1,k)
    j = rowcol(2,k)
    if ( j == i - 1 ) then
      a(k) = -1.0D+00
    else if ( j == i ) then
      a(k) = 2.0D+00
    else if ( j == i + 1 ) then
      a(k) = -1.0D+00
    else
      a(k) = 0.0D+00
    end if

  end do

  return
end
subroutine r8ncf_dif2_nz_num ( m, n, nz_num )

!*****************************************************************************80
!
!! R8NCF_DIF2_NZ_NUM counts nonzeros in an R8NCF second difference matrix.
!
!  Discussion:
!
!    The R8NCF storage format stores NZ_NUM, the number of nonzeros, 
!    a real array containing the nonzero values, a 2 by NZ_NUM integer 
!    array storing the row and column of each nonzero entry.
!
!    The R8NCF format is used by NSPCG.  NSPCG requires that the information
!    for the diagonal entries of the matrix must come first.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 July 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns in
!    the matrix.
!
!    Output, integer NZ_NUM, the number of nonzero entries.
!
  implicit none

  integer nz_num

  integer m
  integer n

  if ( m < n ) then
    nz_num = 3 * m - 1
  else if ( m == n ) then
    nz_num = 3 * n - 2
  else
    nz_num = 3 * n - 1
  end if

  return
end
subroutine r8ncf_dif2_rowcol ( m, n, nz_num, rowcol )

!*****************************************************************************80
!
!! R8NCF_DIF2_ROWCOL sets indexing for an R8NCF second difference matrix.
!
!  Discussion:
!
!    The R8NCF storage format stores NZ_NUM, the number of nonzeros, 
!    a real array containing the nonzero values, a 2 by NZ_NUM integer 
!    array storing the row and column of each nonzero entry.
!
!    The R8NCF format is used by NSPCG.  NSPCG requires that the information
!    for the diagonal entries of the matrix must come first.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 July 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns in
!    the matrix.
!
!    Input, integer NZ_NUM, the number of nonzero entries.
!
!    Output, integer ROWCOL(2,NZ_NUM), the coordinates of 
!    the nonzero entries.
!
  implicit none

  integer nz_num

  integer i
  integer j
  integer k
  integer m
  integer n
  integer rowcol(2,nz_num)

  k = 0

  do i = 1, m

    j = i - 1
    if ( 1 <= j .and. j <= n ) then
      k = k + 1
      rowcol(1,k) = i
      rowcol(2,k) = j
    end if

    j = i
    if ( j <= n ) then
      k = k + 1
      rowcol(1,k) = i
      rowcol(2,k) = j
    end if

    j = i + 1
    if ( j <= n ) then
      k = k + 1
      rowcol(1,k) = i
      rowcol(2,k) = j
    end if

  end do

  return
end
subroutine r8ncf_indicator ( m, n, nz_num, rowcol, a )

!*****************************************************************************80
!
!! R8NCF_INDICATOR sets up an R8NCF indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R8NCF storage format stores NZ_NUM, the number of nonzeros, 
!    a real array containing the nonzero values, a 2 by NZ_NUM integer 
!    array storing the row and column of each nonzero entry.
!
!    The R8NCF format is used by NSPCG.  NSPCG requires that the information
!    for the diagonal entries of the matrix must come first.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns in
!    the matrix.
!
!    Input, integer NZ_NUM, the number of nonzero entries.
!
!    Input, integer ROWCOL(2,NZ_NUM), the coordinates of 
!    the nonzero entries.
!
!    Output, real ( kind = rk ) A(NZ_NUM), the indicator matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nz_num

  real ( kind = rk ) a(nz_num)
  integer fac
  integer i
  integer i4_log_10
  integer j
  integer k
  integer m
  integer n
  integer rowcol(2,nz_num)

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  do k = 1, nz_num
    i = rowcol(1,k)
    j = rowcol(2,k)
    a(k) = real ( fac * i + j, kind = rk )
  end do

  return
end
subroutine r8ncf_mtv ( m, n, nz_num, rowcol, a, x, b )

!*****************************************************************************80
!
!! R8NCF_MTV multiplies an R8VEC times an R8NCF matrix.
!
!  Discussion:
!
!    The R8NCF storage format stores NZ_NUM, the number of nonzeros, 
!    a real array containing the nonzero values, a 2 by NZ_NUM integer 
!    array storing the row and column of each nonzero entry.
!
!    The R8NCF format is used by NSPCG.  NSPCG requires that the information
!    for the diagonal entries of the matrix must come first.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 July 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the order of the matrix.
!
!    Input, integer NZ_NUM, the number of nonzero elements in
!    the matrix.
!
!    Input, integer ROWCOL(2,NZ_NUM), the row and column
!    indices of the nonzero elements.
!
!    Input, real ( kind = rk ) A(NZ_NUM), the nonzero elements of the matrix.
!
!    Input, real ( kind = rk ) X(M), the vector to be multiplied by A'.
!
!    Output, real ( kind = rk ) B(N), the product A' * x.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n
  integer nz_num

  real ( kind = rk ) a(nz_num)
  real ( kind = rk ) b(n)
  integer i
  integer j
  integer k
  integer rowcol(2,nz_num)
  real ( kind = rk ) x(m)

  b(1:n) = 0.0D+00

  do k = 1, nz_num
    i = rowcol(1,k)
    j = rowcol(2,k)
    b(j) = b(j) + a(k) * x(i)
  end do

  return
end
subroutine r8ncf_mv ( m, n, nz_num, rowcol, a, x, b )

!*****************************************************************************80
!
!! R8NCF_MV multiplies an R8NCF matrix by an R8VEC.
!
!  Discussion:
!
!    The R8NCF storage format stores NZ_NUM, the number of nonzeros, 
!    a real array containing the nonzero values, a 2 by NZ_NUM integer 
!    array storing the row and column of each nonzero entry.
!
!    The R8NCF format is used by NSPCG.  NSPCG requires that the information
!    for the diagonal entries of the matrix must come first.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 July 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the order of the matrix.
!
!    Input, integer NZ_NUM, the number of nonzero elements in 
!    the matrix.
!
!    Input, integer ROWCOL(2,NZ_NUM), the row and column 
!    indices of the nonzero elements.
!
!    Input, real ( kind = rk ) A(NZ_NUM), the nonzero elements of the matrix.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = rk ) B(M), the product A * x.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n
  integer nz_num

  real ( kind = rk ) a(nz_num)
  real ( kind = rk ) b(m)
  integer i
  integer j
  integer k
  integer rowcol(2,nz_num)
  real ( kind = rk ) x(n)

  b(1:m) = 0.0D+00

  do k = 1, nz_num
    i = rowcol(1,k)
    j = rowcol(2,k)
    b(i) = b(i) + a(k) * x(j)
  end do

  return
end
subroutine r8ncf_print ( m, n, nz_num, rowcol, a, title )

!*****************************************************************************80
!
!! R8NCF_PRINT prints an R8NCF matrix.
!
!  Discussion:
!
!    The R8NCF storage format stores NZ_NUM, the number of nonzeros, 
!    a real array containing the nonzero values, a 2 by NZ_NUM integer 
!    array storing the row and column of each nonzero entry.
!
!    The R8NCF format is used by NSPCG.  NSPCG requires that the information
!    for the diagonal entries of the matrix must come first.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns of 
!    the matrix.
!
!    Input, integer NZ_NUM, the number of nonzero elements in 
!    the matrix.
!
!    Input, integer ROWCOL(2,NZ_NUM), the row and column indices
!    of the nonzero elements.
!
!    Input, real ( kind = rk ) A(NZ_NUM), the nonzero elements of the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
! 
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nz_num

  real ( kind = rk ) a(nz_num)
  integer m
  integer n
  integer rowcol(2,nz_num)
  character ( len = * ) title

  call r8ncf_print_some ( m, n, nz_num, rowcol, a, 1, 1, m, n, title )

  return
end
subroutine r8ncf_print_some ( m, n, nz_num, rowcol, a, ilo, jlo, &
  ihi, jhi, title )

!*****************************************************************************80
!
!! R8NCF_PRINT_SOME prints some of an R8NCF matrix.
!
!  Discussion:
!
!    The R8NCF storage format stores NZ_NUM, the number of nonzeros, 
!    a real array containing the nonzero values, a 2 by NZ_NUM integer 
!    array storing the row and column of each nonzero entry.
!
!    The R8NCF format is used by NSPCG.  NSPCG requires that the information
!    for the diagonal entries of the matrix must come first.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns of 
!    the matrix.
!
!    Input, integer NZ_NUM, the number of nonzero elements 
!    in the matrix.
!
!    Input, integer ROWCOL(2,NZ_NUM), the row and column indices
!    of the nonzero elements.
!
!    Input, real ( kind = rk ) A(NZ_NUM), the nonzero elements of the matrix.
!
!    Input, integer ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: incx = 5
  integer nz_num

  real ( kind = rk ) a(nz_num)
  real ( kind = rk ) aij
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
  integer n
  logical nonzero
  integer rowcol(2,nz_num)
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
      nonzero = .false.

      aij = 0.0D+00
      do j2 = 1, inc
        write ( ctemp(j2), '(f8.0,6x)' ) aij
      end do

      do k = 1, nz_num

        if ( &
          i == rowcol(1,k) .and. &
          j2lo <= rowcol(2,k) .and. &
          rowcol(2,k) <= j2hi ) then 

          j2 = rowcol(2,k) - j2lo + 1
          aij = a(k)

          if ( aij == 0.0D+00 ) then
            cycle
          end if

          nonzero = .true.

          write ( ctemp(j2), '(g14.6)' ) aij

        end if

      end do

      if ( nonzero ) then
        write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )
      end if

    end do

  end do

  return
end
subroutine r8ncf_random ( m, n, nz_num, rowcol, seed, a )

!*****************************************************************************80
!
!! R8NCF_RANDOM randomizes an R8NCF matrix.
!
!  Discussion:
!
!    The R8NCF storage format stores NZ_NUM, the number of nonzeros, 
!    a real array containing the nonzero values, a 2 by NZ_NUM integer 
!    array storing the row and column of each nonzero entry.
!
!    The R8NCF format is used by NSPCG.  NSPCG requires that the information
!    for the diagonal entries of the matrix must come first.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 July 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns in
!    the matrix.
!
!    Input, integer NZ_NUM, the number of nonzero entries.
!
!    Input, integer ROWCOL(2,NZ_NUM), the coordinates of 
!    the nonzero entries.
!
!    Input/output, integer SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = rk ) A(NZ_NUM), the indicator matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nz_num

  real ( kind = rk ) a(nz_num)
  integer i
  integer j
  integer k
  integer m
  integer n
  real ( kind = rk ) r8_uniform_01
  integer rowcol(2,nz_num)
  integer seed

  do k = 1, nz_num
    i = rowcol(1,k)
    j = rowcol(2,k)
    a(k) = r8_uniform_01 ( seed )
  end do

  return
end
subroutine r8ncf_to_r8ge ( m, n, nz_num, rowcol, a, a_r8ge )

!*****************************************************************************80
!
!! R8NCF_TO_R8GE converts an R8NCF matrix to R8GE format.
!
!  Discussion:
!
!    The R8NCF storage format stores NZ_NUM, the number of nonzeros, 
!    a real array containing the nonzero values, a 2 by NZ_NUM integer 
!    array storing the row and column of each nonzero entry.
!
!    The R8NCF format is used by NSPCG.  NSPCG requires that the information
!    for the diagonal entries of the matrix must come first.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 July 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns in
!    the matrix.
!
!    Input, integer NZ_NUM, the number of nonzero entries.
!
!    Input, integer ROWCOL(2,NZ_NUM), the coordinates of 
!    the nonzero entries.
!
!    Input, real ( kind = rk ) A(NZ_NUM), the matrix.
!
!    Output, real ( kind = rk ) A_R8GE(M,N), the R8GE matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n
  integer nz_num

  real ( kind = rk ) a(nz_num)
  real ( kind = rk ) a_r8ge(m,n)
  integer i
  integer j
  integer k
  integer rowcol(2,nz_num)

  a_r8ge(1:m,1:n) = 0.0D+00

  do k = 1, nz_num
    i = rowcol(1,k)
    j = rowcol(2,k)
    a_r8ge(i,j) = a_r8ge(i,j) + a(k)
  end do

  return
end
subroutine r8ncf_zeros ( m, n, nz_num, rowcol, a )

!*****************************************************************************80
!
!! R8NCF_ZEROS zeroes an R8NCF matrix.
!
!  Discussion:
!
!    The R8NCF storage format stores NZ_NUM, the number of nonzeros, 
!    a real array containing the nonzero values, a 2 by NZ_NUM integer 
!    array storing the row and column of each nonzero entry.
!
!    The R8NCF format is used by NSPCG.  NSPCG requires that the information
!    for the diagonal entries of the matrix must come first.
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
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns of 
!    the matrix.
!
!    Input, integer NZ_NUM, the number of nonzero elements 
!    in the matrix.
!
!    Input, integer ROWCOL(2,NZ_NUM), the row and column indices
!    of the nonzero elements.
!
!    Input, real ( kind = rk ) A(NZ_NUM), the nonzero elements of the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nz_num

  real ( kind = rk ) a(nz_num)
  integer m
  integer n
  integer rowcol(2,nz_num)

  a(1:nz_num) = 0.0D+00

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

