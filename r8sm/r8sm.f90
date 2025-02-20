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
subroutine r8ge_fa ( n, a, pivot, info )

!*****************************************************************************80
!
!! R8GE_FA performs a LINPACK style PLU factorization of an R8GE matrix.
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
!    R8GE_FA is a simplified version of the LINPACK routine SGEFA.
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
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
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
!    Output, integer PIVOT(N), a vector of pivot indices.
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
  integer info
  integer pivot(n)
  integer j
  integer k
  integer l
  real ( kind = rk ) t

  info = 0

  do k = 1, n - 1
!
!  Find L, the index of the pivot row.
!
    l = k
    do i = k + 1, n
      if ( abs ( a(l,k) ) < abs ( a(i,k) ) ) then
        l = i
      end if
    end do

    pivot(k) = l
!
!  If the pivot index is zero, the algorithm has failed.
!
    if ( a(l,k) == 0.0D+00 ) then
      info = k
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8GE_FA - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      stop 1
    end if
!
!  Interchange rows L and K if necessary.
!
    if ( l /= k ) then
      t      = a(l,k)
      a(l,k) = a(k,k)
      a(k,k) = t
    end if
!
!  Normalize the values that lie below the pivot entry A(K,K).
!
    a(k+1:n,k) = -a(k+1:n,k) / a(k,k)
!
!  Row elimination with column indexing.
!
    do j = k + 1, n

      if ( l /= k ) then
        t      = a(l,j)
        a(l,j) = a(k,j)
        a(k,j) = t
      end if

      a(k+1:n,j) = a(k+1:n,j) + a(k+1:n,k) * a(k,j)

    end do

  end do

  pivot(n) = n

  if ( a(n,n) == 0.0D+00 ) then
    info = n
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_FA - Fatal error!'
    write ( *, '(a,i8)' ) '  Zero pivot on step ', info
    stop 1
  end if

  return
end
subroutine r8ge_ml ( n, a_lu, pivot, x, b, job )

!*****************************************************************************80
!
!! R8GE_ML computes A * x or A' * x, using R8GE_FA factors.
!
!  Discussion:
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!    R8GE storage is used by LINPACK and LAPACK.
!
!    It is assumed that R8GE_FA has overwritten the original matrix
!    information by LU factors.  R8GE_ML is able to reconstruct the
!    original matrix from the LU factor data.
!
!    R8GE_ML allows the user to check that the solution of a linear
!    system is correct, without having to save an unfactored copy
!    of the matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    31 October 1998
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
!    Input, real ( kind = rk ) A_LU(N,N), the LU factors from R8GE_FA.
!
!    Input, integer PIVOT(N), the pivot vector computed by R8GE_FA.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied.
!
!    Output, real ( kind = rk ) B(N), the result of the multiplication.
!
!    Input, integer JOB, specifies the operation to be done:
!    JOB = 0, compute A * x.
!    JOB nonzero, compute A' * X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a_lu(n,n)
  real ( kind = rk ) b(n)
  integer i
  integer pivot(n)
  integer j
  integer job
  integer k
  real ( kind = rk ) t
  real ( kind = rk ) x(n)

  b(1:n) = x(1:n)

  if ( job == 0 ) then
!
!  Y = U * X.
!
    do j = 1, n
      b(1:j-1) = b(1:j-1) + a_lu(1:j-1,j) * b(j)
      b(j) = a_lu(j,j) * b(j)
    end do
!
!  B = PL * Y = PL * U * X = A * x.
!
    do j = n-1, 1, -1

      b(j+1:n) = b(j+1:n) - a_lu(j+1:n,j) * b(j)

      k = pivot(j)

      if ( k /= j ) then
        t    = b(k)
        b(k) = b(j)
        b(j) = t
      end if

    end do

  else
!
!  Y = (PL)' * X:
!
    do j = 1, n - 1

      k = pivot(j)

      if ( k /= j ) then
        t    = b(k)
        b(k) = b(j)
        b(j) = t
      end if

      b(j) = b(j) - sum ( b(j+1:n) * a_lu(j+1:n,j) )

    end do
!
!  B = U' * Y = ( PL * U )' * X = A' * X.
!
    do i = n, 1, -1
      b(i+1:n) = b(i+1:n) + b(i) * a_lu(i,i+1:n)
      b(i) = b(i) * a_lu(i,i)
    end do

  end if

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
subroutine r8ge_random ( m, n, seed, a )

!*****************************************************************************80
!
!! R8GE_RANDOM randomizes an R8GE matrix.
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
!    11 August 2004
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
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns in
!    the array.
!
!    Input/output, integer SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = rk ) A(M,N), the array.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  integer i
  integer, parameter :: i4_huge = 2147483647
  integer j
  integer k
  integer seed

  do j = 1, n

    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge
      end if

      a(i,j) = real ( seed, kind = rk ) * 4.656612875D-10

    end do
  end do

  return
end
subroutine r8ge_sl ( n, a_lu, pivot, b, job )

!*****************************************************************************80
!
!! R8GE_SL solves a system factored by R8GE_FA.
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
!    R8GE_SL is a simplified version of the LINPACK routine SGESL.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 March 1999
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
!    Input, real ( kind = rk ) A_LU(N,N), the LU factors from R8GE_FA.
!
!    Input, integer PIVOT(N), the pivot vector from R8GE_FA.
!
!    Input/output, real ( kind = rk ) B(N).
!    On input, the right hand side vector.
!    On output, the solution vector.
!
!    Input, integer JOB, specifies the operation.
!    0, solve A * x = b.
!    nonzero, solve A' * x = b.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a_lu(n,n)
  real ( kind = rk ) b(n)
  integer pivot(n)
  integer job
  integer k
  integer l
  real ( kind = rk ) t
!
!  Solve A * x = b.
!
  if ( job == 0 ) then
!
!  Solve PL * Y = B.
!
    do k = 1, n - 1

      l = pivot(k)

      if ( l /= k ) then
        t    = b(l)
        b(l) = b(k)
        b(k) = t
      end if

      b(k+1:n) = b(k+1:n) + a_lu(k+1:n,k) * b(k)

    end do
!
!  Solve U * X = Y.
!
    do k = n, 1, -1
      b(k) = b(k) / a_lu(k,k)
      b(1:k-1) = b(1:k-1) - a_lu(1:k-1,k) * b(k)
    end do
!
!  Solve A' * X = B.
!
  else
!
!  Solve U' * Y = B.
!
    do k = 1, n
      b(k) = ( b(k) - sum ( b(1:k-1) * a_lu(1:k-1,k) ) ) / a_lu(k,k)
    end do
!
!  Solve ( PL )' * X = Y.
!
    do k = n - 1, 1, -1

      b(k) = b(k) + sum ( b(k+1:n) * a_lu(k+1:n,k) )

      l = pivot(k)

      if ( l /= k ) then
        t    = b(l)
        b(l) = b(k)
        b(k) = t
      end if

    end do

  end if

  return
end
subroutine r8sm_indicator ( m, n, a, u, v )

!*****************************************************************************80
!
!! R8SM_INDICATOR returns the indicator matrix as an R8SM matrix.
!
!  Discussion:
!
!    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
!    which is defined by an M by N matrix A, an M vector U, and
!    an N vector V, by B = A - U * V'
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 May 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns 
!    of the matrix.
!
!    Output, real ( kind = rk ) A(M,N), the R8SM matrix.
!
!    Output, real ( kind = rk ) U(M), V(N), the R8SM vectors.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  integer fac
  integer i
  integer i4_log_10
  real ( kind = rk ) u(m)
  real ( kind = rk ) v(n)

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  u(1:m) = - 1.0D+00

  do i = 1, n
    v(i) = real ( i, kind = rk )
  end do

  do i = 1, m
    a(i,1:n) = real ( fac * i, kind = rk )
  end do

  return
end
subroutine r8sm_ml ( n, a_lu, u, v, pivot, x, b, job )

!*****************************************************************************80
!
!! R8SM_ML multiplies a factored square R8SM matrix by an R8VEC.
!
!  Discussion:
!
!    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
!    which is defined by an M by N matrix A, an M vector U, and
!    an N vector V, by B = A - U * V'
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 January 2004
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
!    Input, real ( kind = rk ) A_LU(N,N), the LU factors from R8GE_FA.
!
!    Input, real ( kind = rk ) U(N), V(N), the Sherman Morrison vectors.
!
!    Input, integer PIVOT(N), the pivot vector computed by R8GE_FA.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied.
!
!    Output, real ( kind = rk ) B(N), the result of the multiplication.
!
!    Input, integer JOB, specifies the operation to be done:
!    JOB = 0, compute (A-u*v') * x.
!    JOB nonzero, compute (A-u*v')' * x.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a_lu(n,n)
  real ( kind = rk ) b(n)
  integer pivot(n)
  integer job
  real ( kind = rk ) u(n)
  real ( kind = rk ) v(n)
  real ( kind = rk ) x(n)

  call r8ge_ml ( n, a_lu, pivot, x, b, job )

  if ( job == 0 ) then

    b(1:n) = b(1:n) - u(1:n) * dot_product ( v, x )

  else

    b(1:n) = b(1:n) - v(1:n) * dot_product ( u, x )

  end if

  return
end
subroutine r8sm_mtv ( m, n, a, u, v, x, b )

!*****************************************************************************80
!
!! R8SM_MTV multiplies an R8VEC by an R8SM matrix.
!
!  Discussion:
!
!    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
!    which is defined by an M by N matrix A, an M vector U, and
!    an N vector V, by B = A - U * V'
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns 
!    of the matrix.
!
!    Input, real ( kind = rk ) A(M,N), the R8SM matrix.
!
!    Input, real ( kind = rk ) U(M), V(N), the R8SM vectors.
!
!    Input, real ( kind = rk ) X(M), the vector to be multiplied.
!
!    Output, real ( kind = rk ) B(N), the product (A-u*v')' * X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) b(n)
  real ( kind = rk ) u(m)
  real ( kind = rk ) v(n)
  real ( kind = rk ) x(m)

  b(1:n) = matmul ( transpose ( a(1:m,1:n) ), x(1:m) ) &
    - v(1:n) * sum ( u(1:m) * x(1:m) )

  return
end
subroutine r8sm_mv ( m, n, a, u, v, x, b )

!*****************************************************************************80
!
!! R8SM_MV multiplies an R8SM matrix by an R8VEC.
!
!  Discussion:
!
!    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
!    which is defined by an M by N matrix A, an M vector U, and
!    an N vector V, by B = A - U * V'
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 May 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns 
!    of the matrix.
!
!    Input, real ( kind = rk ) A(M,N), the R8SM matrix.
!
!    Input, real ( kind = rk ) U(M), V(N), the R8SM vectors U and V.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied by (A-u*v').
!
!    Output, real ( kind = rk ) B(M), the product (A-u*v') * x.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) b(m)
  real ( kind = rk ) u(m)
  real ( kind = rk ) v(n)
  real ( kind = rk ) x(n)

  b(1:m) = matmul ( a(1:m,1:n), x(1:n) ) - u(1:m) * dot_product ( v, x )

  return
end
subroutine r8sm_print ( m, n, a, u, v, title )

!*****************************************************************************80
!
!! R8SM_PRINT prints an R8SM matrix.
!
!  Discussion:
!
!    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
!    which is defined by an M by N matrix A, an M vector U, and
!    an N vector V, by B = A - U * V'
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 January 2004
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
!    Input, real ( kind = rk ) A(M,N), the R8SM matrix.
!
!    Input, real ( kind = rk ) U(M), V(N), the R8SM vectors.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  character ( len = * ) title
  real ( kind = rk ) u(m)
  real ( kind = rk ) v(n)

  call r8sm_print_some ( m, n, a, u, v, 1, 1, m, n, title )

  return
end
subroutine r8sm_print_some ( m, n, a, u, v, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8SM_PRINT_SOME prints some of an R8SM matrix.
!
!  Discussion:
!
!    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
!    which is defined by an M by N matrix A, an M vector U, and
!    an N vector V, by B = A - U * V'
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns 
!    of the matrix.
!
!    Input, real ( kind = rk ) A(M,N), the R8SM matrix.
!
!    Input, real ( kind = rk ) U(M), V(N), the R8SM vectors.
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
  real ( kind = rk ) u(n)
  real ( kind = rk ) v(n)
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

    write ( *, '(a,5a14)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
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

        aij = a(i,j) - u(i) * v(j)

        write ( ctemp(j2), '(g14.6)' ) aij

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8sm_random ( m, n, seed, a, u, v )

!*****************************************************************************80
!
!! R8SM_RANDOM randomizes an R8SM matrix.
!
!  Discussion:
!
!    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
!    which is defined by an M by N matrix A, an M vector U, and
!    an N vector V, by B = A - U * V'
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns 
!    of the matrix.
!
!    Input/output, integer SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = rk ) A(M,N), the R8SM matrix.
!
!    Output, real ( kind = rk ) U(M), V(N), the R8SM vectors.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  integer seed
  real ( kind = rk ) u(m)
  real ( kind = rk ) v(n)

  call r8ge_random ( m, n, seed, a )
  
  call r8vec_uniform_01 ( m, seed, u )

  call r8vec_uniform_01 ( n, seed, v )

  return
end
subroutine r8sm_sl ( n, a_lu, u, v, b, ierror, pivot )

!*****************************************************************************80
!
!! R8SM_SL solves B*x=b, where B is a square R8SM matrix.
!
!  Discussion:
!
!    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
!    which is defined by an M by N matrix A, an M vector U, and
!    an N vector V, by B = A - U * V'
!
!    It is assumed that A has been decomposed into its LU factors
!    by R8GE_FA.  The Sherman Morrison formula allows
!    us to solve linear systems involving (A-u*v') by solving linear
!    systems involving A and adjusting the results.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 October 2015
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash
!    Numerical Methods and Software,
!    Prentice Hall, 1989
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real ( kind = rk ) A_LU(N,N), the LU factors from R8GE_FA.
!
!    Input, real ( kind = rk ) U(N), V(N), the R8SM vectors U and V.
!
!    Input/output, real ( kind = rk ) B(N).
!    On input, the right hand side vector.
!    On output, the solution vector.
!
!    Output, integer IERROR, an error flag.
!    0, no error occurred.  The solution was successfully computed.
!    1, an error occurred.  1 - v' * Inverse(A) * u = 0.
!    The solution was not computed.
!
!    Input, integer PIVOT(N), the pivot vector produced by R8GE_FA.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a_lu(n,n)
  real ( kind = rk ) alpha
  real ( kind = rk ) b(n)
  real ( kind = rk ) beta
  integer ierror
  integer pivot(n)
  integer job_local
  real ( kind = rk ) u(n)
  real ( kind = rk ) v(n)
  real ( kind = rk ) w(n)

  ierror = 0
!
!  Solve A' * w = v.
!
  w(1:n) = v(1:n)

  job_local = 1
  call r8ge_sl ( n, a_lu, pivot, w, job_local )
!
!  Set beta = w' * b.
!
  beta = dot_product ( w, b )
!
!  Solve A * b = b.
!
  job_local = 0
  call r8ge_sl ( n, a_lu, pivot, b, job_local )
!
!  Solve A * w = u.
!
  w(1:n) = u(1:n)

  job_local = 0
  call r8ge_sl ( n, a_lu, pivot, w, job_local )
!
!  Set alpha = 1 / ( 1 - v' * w ).
!
  alpha = 1.0D+00 - dot_product ( v, w )

  if ( alpha == 0.0D+00 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8SM_SL - Fatal error!'
    write ( *, '(a)' ) '  The divisor ALPHA is zero.'
    stop 1
  end if

  alpha = 1.0D+00 / alpha
!
!  Set b = b + alpha * beta * w.
!
  b(1:n) = b(1:n) + alpha * beta * w(1:n)

  return
end
subroutine r8sm_slt ( n, a_lu, u, v, b, ierror, pivot )

!*****************************************************************************80
!
!! R8SM_SLT solves B'*x=b, where B is an R8SM matrix.
!
!  Discussion:
!
!    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
!    which is defined by an M by N matrix A, an M vector U, and
!    an N vector V, by B = A - U * V'
!
!    It is assumed that A has been decomposed into its LU factors
!    by R8GE_FA.  The Sherman Morrison formula allows
!    us to solve linear systems involving (A-u*v') by solving linear
!    systems involving A and adjusting the results.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash
!    Numerical Methods and Software,
!    Prentice Hall, 1989
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real ( kind = rk ) A_LU(N,N), the LU factors from R8GE_FA.
!
!    Input, real ( kind = rk ) U(N), V(N), the R8SM vectors U and V.
!
!    Input/output, real ( kind = rk ) B(N).
!    On input, the right hand side vector.
!    On output, the solution vector.
!
!    Output, integer IERROR, an error flag.
!    0, no error occurred.  The solution was successfully computed.
!    1, an error occurred.  1 - v' * Inverse(A) * u = 0.
!    The solution was not computed.
!
!    Input, integer PIVOT(N), the pivot vector produced by R8GE_FA.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a_lu(n,n)
  real ( kind = rk ) alpha
  real ( kind = rk ) b(n)
  real ( kind = rk ) beta
  integer ierror
  integer pivot(n)
  integer job_local
  real ( kind = rk ) u(n)
  real ( kind = rk ) v(n)
  real ( kind = rk ) w(n)

  ierror = 0
!
!  Solve A * w = u.
!
  w(1:n) = u(1:n)

  job_local = 0
  call r8ge_sl ( n, a_lu, pivot, w, job_local )
!
!  Set beta = w' * b.
!
  beta = dot_product ( w, b )
!
!  Solve A' * b = b.
!
  job_local = 1
  call r8ge_sl ( n, a_lu, pivot, b, job_local )
!
!  Solve A' * w = v.
!
  w(1:n) = v(1:n)

  job_local = 1
  call r8ge_sl ( n, a_lu, pivot, w, job_local )
!
!  Set alpha = 1 / ( 1 - u' * w ).
!
  alpha = 1.0D+00 - dot_product ( u, w )

  if ( alpha == 0.0D+00 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8SM_SLT - Fatal error!'
    write ( *, '(a)' ) '  The divisor ALPHA is zero.'
    stop 1
  end if

  alpha = 1.0D+00 / alpha
!
!  Set b = b + alpha * beta * w.
!
  b(1:n) = b(1:n) + alpha * beta * w(1:n)

  return
end
subroutine r8sm_to_r8ge ( m, n, a, u, v, b )

!*****************************************************************************80
!
!! R8SM_TO_R8GE copies an R8SM matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
!    which is defined by an M by N matrix A, an M vector U, and
!    an N vector V, by B = A - U * V'
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
!    17 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns 
!    of the matrix.
!
!    Input, real ( kind = rk ) A(M,N), the R8SM matrix.
!
!    Input, real ( kind = rk ) U(M), V(N), the R8SM vectors.
!
!    Output, real ( kind = rk ) B(M,N), the R8GE matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) b(m,n)
  integer i
  real ( kind = rk ) u(m)
  real ( kind = rk ) v(n)

  do i = 1, m
    b(i,1:n) = a(i,1:n) - u(i) * v(1:n)
  end do

  return
end
subroutine r8sm_zeros ( m, n, a, u, v )

!*****************************************************************************80
!
!! R8SM_ZEROS zeroes an R8SM matrix.
!
!  Discussion:
!
!    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
!    which is defined by an M by N matrix A, an M vector U, and
!    an N vector V, by B = A - U * V'
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
!    Input, integer M, N, the number of rows and columns 
!    of the matrix.
!
!    Output, real ( kind = rk ) A(M,N), the R8SM matrix.
!
!    Output, real ( kind = rk ) U(M), V(N), the R8SM vectors.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) u(m)
  real ( kind = rk ) v(n)

  a(1:m,1:n) = 0.0D+00
  u(1:m) = 0.0D+00
  v(1:n) = 0.0D+00

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
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 August 2014
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
!    Input, integer N, the number of entries in the vector.
!
!    Input/output, integer SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = rk ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer, parameter :: i4_huge = 2147483647
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
      seed = seed + i4_huge
    end if

    r(i) = real ( seed, kind = rk ) * 4.656612875D-10

  end do

  return
end
subroutine r8vec2_print_some ( n, x1, x2, max_print, title )

!*****************************************************************************80
!
!! R8VEC2_PRINT_SOME prints "some" of a pair of R8VEC's.
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vectors, is no more than MAX_PRINT, then
!    the entire vectors are printed, one entry of each per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries of the vectors.
!
!    Input, real ( kind = rk ) X1(N), X2(N), the vector to be printed.
!
!    Input, integer MAX_PRINT, the maximum number of lines 
!    to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer max_print
  character ( len = * ) title
  real ( kind = rk ) x1(n)
  real ( kind = rk ) x2(n)

  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  if ( n <= max_print ) then

    do i = 1, n
      write ( *, '(i8,2x,g14.6,2x,g14.6)' ) i, x1(i), x2(i)
    end do

  else if ( 3 <= max_print ) then

    do i = 1, max_print - 2
      write ( *, '(i8,2x,g14.6,2x,g14.6)' ) i, x1(i), x2(i)
    end do
    write ( *, '(a)' ) '......  ..............  ..............'
    i = n
    write ( *, '(i8,2x,g14.6,2x,g14.6)' ) i, x1(i), x2(i)

  else

    do i = 1, max_print - 1
      write ( *, '(i8,2x,g14.6,2x,g14.6)' ) i, x1(i), x2(i)
    end do
    i = max_print
    write ( *, '(i8,2x,g14.6,2x,g14.6,2x,a)' ) i, x1(i), x2(i), &
      '...more entries...'

  end if

  return
end

