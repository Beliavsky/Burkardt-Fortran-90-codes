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
subroutine r8sto_dif2 ( n, a )

!*****************************************************************************80
!
!! R8STO_DIF2 sets the second difference as an R8STO matrix.
!
!  Discussion:
!
!    The R8STO storage format is used for a symmetric Toeplitz matrix.
!    It stores the N elements of the first row.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 September 2015
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
!    Output, real ( kind = rk ) A(N), the R8STO matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)

  a(1:n) = 0.0D+00

  a(1) = 2.0D+00
  a(2) = -1.0D+00
  
  return
end
subroutine r8sto_indicator ( n, a )

!*****************************************************************************80
!
!! R8STO_INDICATOR sets up an R8STO indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R8STO storage format is used for a symmetric Toeplitz matrix.
!    It stores the N elements of the first row.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 January 2004
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
!    Output, real ( kind = rk ) A(N), the R8STO matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  integer fac
  integer i
  integer i4_log_10
  integer j
  integer k

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  i = 1
  k = 0
  do j = 1, n
    k = k + 1
    a(k) = real ( fac * i + j, kind = rk )
  end do
  
  return
end
subroutine r8sto_inverse ( n, a, b )

!*****************************************************************************80
!
!! R8STO_INVERSE computes the inverse of an R8STO matrix.
!
!  Discussion:
!
!    The R8STO storage format is used for a symmetric Toeplitz matrix.
!    It stores the N elements of the first row.
!
!    For this routine, the matrix is also required to be positive definite.
!
!    The original implementation of the algorithm assumed that the
!    diagonal element was 1.  The algorithm has been modified so that
!    this is no longer necessary.
!
!    The inverse matrix is NOT guaranteed to be a Toeplitz matrix.  
!    It is guaranteed to be symmetric and persymmetric.
!    The inverse matrix is returned in general storage, that is,
!    as an "SGE" matrix.
!
!  Example:
!
!    To compute the inverse of
!
!     1.0 0.5 0.2
!     0.5 1.0 0.5
!     0.2 0.5 1.0
!
!    we input:
!
!      N = 3
!      A = (/ 1.0, 0.5, 0.2 /)
!
!    with output:
!
!      B(1:3,1:3) = (/ 75, -40,   5,
!                     -40,  96, -40,
!                       5, -40,  75 /) / 56
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    29 May 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Gene Golub, Charles Van Loan,
!    Section 4.7.3, "Computing the Inverse",
!    Matrix Computations,
!    Third Edition,
!    Johns Hopkins, 1996.
!
!  Parameters:
!
!    Input, integer N, the order of the system.
!
!    Input, real ( kind = rk ) A(N), the R8STO matrix.
!
!    Output, real ( kind = rk ) B(N,N), the inverse of the matrix, 
!    in R8GE format.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) a2(n-1)
  real ( kind = rk ) b(n,n)
  integer i
  integer j
  real ( kind = rk ) v(n)

  a2(1:n-1) = a(2:n) / a(1)

  call r8sto_yw_sl ( n - 1, a2, v )
!
!  Compute the N-th entry of V.
!
  v(n) = 1.0D+00 / ( 1.0D+00 + sum ( a2(1:n-1) * v(1:n-1) ) )
!
!  Reverse and scale entries 1 through N-1.
!
  v(1:n-1) = v(n-1:1:-1)

  v(1:n-1) = v(n) * v(1:n-1)
!
!  Set the boundaries of B.
!
  b(1,1:n) = v(n:1:-1)
  b(n,1:n) = v(1:n)
  b(2:n-1,1) = v(n-1:2:-1)
  b(2:n-1,n) = v(2:n-1)
!
!  Fill the interior.
!
  do i = 2, 1 + ( ( n - 1 ) / 2 )
    do j = i, n - i + 1
      b(i,j) = b(i-1,j-1) + ( v(n+1-j) * v(n+1-i) - v(i-1) * v(j-1) ) / v(n)
      b(j,i) = b(i,j)
      b(n+1-i,n+1-j) = b(i,j)
      b(n+1-j,n+1-i) = b(i,j)
    end do
  end do
!
!  Scale B.
!
  b(1:n,1:n) = b(1:n,1:n) / a(1)

  return
end
subroutine r8sto_mv ( n, a, x, b )

!*****************************************************************************80
!
!! R8STO_MV multiplies an R8STO matrix by an R8VEC.
!
!  Discussion:
!
!    The R8STO storage format is used for a symmetric Toeplitz matrix.
!    It stores the N elements of the first row.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real ( kind = rk ) A(N), the R8STO matrix.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = rk ) B(N), the product A * x.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) b(n)
  integer i
  real ( kind = rk ) x(n)

  do i = 1, n
    b(i) = sum ( a(i:2:-1) * x(1:i-1) ) + sum ( a(1:n+1-i) * x(i:n) )
  end do

  return
end
subroutine r8sto_print ( n, a, title )

!*****************************************************************************80
!
!! R8STO_PRINT prints an R8STO matrix.
!
!  Discussion:
!
!    The R8STO storage format is used for a symmetric Toeplitz matrix.
!    It stores the N elements of the first row.
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
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real ( kind = rk ) A(N), the R8STO matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  character ( len = * ) title

  call r8sto_print_some ( n, a, 1, 1, n, n, title )

  return
end
subroutine r8sto_print_some ( n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8STO_PRINT_SOME prints some of an R8STO matrix.
!
!  Discussion:
!
!    The R8STO storage format is used for a symmetric Toeplitz matrix.
!    It stores the N elements of the first row.
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
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, real ( kind = rk ) A(N), the R8STO matrix.
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

  real ( kind = rk ) a(n)
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
    i2hi = min ( ihi, n )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( i <= j ) then
          aij = a(1+j-i)
        else
          aij = a(1+i-j)
        end if

        write ( ctemp(j2), '(g14.6)' ) aij

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8sto_random ( n, seed, a )

!*****************************************************************************80
!
!! R8STO_RANDOM randomizes an R8STO matrix.
!
!  Discussion:
!
!    The R8STO storage format is used for a symmetric Toeplitz matrix.
!    It stores the N elements of the first row.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 May 2003
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
!    Input/output, integer SEED, a seed for the random number
!    generator.
!
!    Output, real ( kind = rk ) A(N), the R8STO matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  integer seed

  call r8vec_uniform_01 ( n, seed, a )

  return
end
subroutine r8sto_sl ( n, a, b, x )

!*****************************************************************************80
!
!! R8STO_SL solves an R8STO system.
!
!  Discussion:
!
!    The R8STO storage format is used for a symmetric Toeplitz matrix.
!    It stores the N elements of the first row.
!
!    For this routine, the matrix is also required to be positive definite.
!
!    This implementation of the algorithm assumes that the diagonal element
!    is 1.
!
!    The real symmetric Toeplitz matrix can be described by N numbers, which,
!    for convenience, we will label A(0:N-1).
!
!    Note that there is a typographical error in the presentation
!    of this algorithm in the reference, and another in the presentation
!    of a sample problem.  Both involve sign errors.  A minor error
!    makes the algorithm incorrect for the case N = 1.
!
!  Example:
!
!    To solve
!
!     1.0 0.5 0.2    x1    4.0
!     0.5 1.0 0.5 *  x2 = -1.0
!     0.2 0.5 1.0    x3    3.0
!
!    we input:
!
!      N = 3
!      A(0:N-1) = (/ 1.0, 0.5, 0.2 /)
!      B(1:3) = (/ 4.0, -1.0, 3.0 /)
!
!    with output:
!
!      X(1:3) = (/ 355, -376, 285 /) / 56
!             = (/ 6.339, -6.714, 5.089 /)
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Gene Golub, Charles Van Loan,
!    Section 4.7.3, "The General Right Hand Side Problem",
!    Matrix Computations,
!    Third Edition,
!    Johns Hopkins, 1996.
!
!  Parameters:
!
!    Input, integer N, the order of the system.
!
!    Input, real ( kind = rk ) A(N), the R8STO matrix.
!
!    Input, real ( kind = rk ) B(N), the right hand side of the linear system.
!
!    Output, real ( kind = rk ) X(N), the solution of the linear system.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) b(n)
  real ( kind = rk ) beta
  integer k
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)

  k = 0

  beta = 1.0D+00
  x(k+1) = b(k+1) / beta

  if ( k < n-1 ) then
    y(k+1) = -a(k+2) / beta
  end if

  do k = 1, n-1

    beta = ( 1.0D+00 - y(k) * y(k) ) * beta

    x(k+1) = ( b(k+1) - sum ( a(2:k+1) * x(k:1:-1) ) ) / beta

    x(1:k) = x(1:k) + x(k+1) * y(k:1:-1)

    if ( k < n - 1 ) then
      y(k+1) = ( -a(k+2) - sum ( a(2:k+1) * y(k:1:-1) ) ) / beta
      y(1:k) = y(1:k) + y(k+1) * y(k:1:-1)
    end if

  end do

  return
end
subroutine r8sto_to_r8ge ( n, a, b )

!*****************************************************************************80
!
!! R8STO_TO_R8GE copies an R8STO matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R8STO storage format is used for a symmetric Toeplitz matrix.
!    It stores the N elements of the first row.
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
!    22 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real ( kind = rk ) A(N), the R8STO matrix.
!
!    Output, real ( kind = rk ) B(N,N), the R8GE matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) b(n,n)
  integer i

  do i = 1, n
    b(i,1:i-1) = a(i:2:-1)
    b(i,i:n) = a(1:n-i+1)
  end do

  return
end
subroutine r8sto_yw_sl ( n, b, x )

!*****************************************************************************80
!
!! R8STO_YW_SL solves the Yule-Walker equations for an R8STO matrix.
!
!  Discussion:
!
!    The R8STO storage format is used for a symmetric Toeplitz matrix.
!    It stores the N elements of the first row.
!
!    The matrix is also required to be positive definite.
!
!    This implementation of the algorithm assumes that the diagonal element
!    is 1.
!
!    The real symmetric Toeplitz matrix can be described by N numbers, which,
!    for convenience, we will label B(0:N-1).  We assume there is one more
!    number, B(N).  If we let A be the symmetric Toeplitz matrix whose first
!    row is B(0:N-1), then the Yule-Walker equations are:
!
!      A * X = -B(1:N)
!
!  Example:
!
!    To solve
!
!     1.0 0.5 0.2    x1   0.5
!     0.5 1.0 0.5 *  x2 = 0.2
!     0.2 0.5 1.0    x3   0.1
!
!    we input:
!
!      N = 3
!      B(1:3) = (/ 0.5, 0.2, 0.1 /)
!
!    with output:
!
!      X(1:3) = (/ -75, 12, -5 /) / 140
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Gene Golub, Charles Van Loan,
!    Section 4.7.2, "Solving the Yule-Walker Equations",
!    Matrix Computations,
!    Third Edition,
!    Johns Hopkins, 1996.
!
!  Parameters:
!
!    Input, integer N, the order of the system.
!
!    Input, real ( kind = rk ) B(N), defines the linear system.  The first
!    entry of A is a 1, followed by B(1) through B(N-1).  The right hand
!    side of the system is -B(1:N).
!
!    Output, real ( kind = rk ) X(N), the solution of the linear system.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) alpha
  real ( kind = rk ) b(n)
  real ( kind = rk ) beta
  integer i
  real ( kind = rk ) x(n)

  x(1) = -b(1)
  beta = 1.0D+00
  alpha = -b(1)

  do i = 1, n-1
    beta = ( 1.0D+00 - alpha * alpha ) * beta
    alpha = - ( b(i+1) + sum ( b(i:1:-1) * x(1:i) ) ) / beta
    x(1:i) = x(1:i) + alpha * x(i:1:-1)
    x(i+1) = alpha
  end do

  return
end
subroutine r8sto_zeros ( n, a )

!*****************************************************************************80
!
!! R8STO_ZEROS zeroes an R8STO matrix.
!
!  Discussion:
!
!    The R8STO storage format is used for a symmetric Toeplitz matrix.
!    It stores the N elements of the first row.
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
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Output, real ( kind = rk ) A(N), the R8STO matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)

  a(1:n) = 0.0D+00

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

