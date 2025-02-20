function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! r8_uniform_01() returns a unit pseudorandom R8.
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
subroutine r8ge_dif2 ( m, n, a )

!*****************************************************************************80
!
!! R8GE_DIF2 returns the DIF2 matrix in R8GE format.
!
!  Example:
!
!    N = 5
!
!    2 -1  .  .  .
!   -1  2 -1  .  .
!    . -1  2 -1  .
!    .  . -1  2 -1
!    .  .  . -1  2
!
!  Properties:
!
!    A is banded, with bandwidth 3.
!
!    A is tridiagonal.
!
!    Because A is tridiagonal, it has property A (bipartite).
!
!    A is a special case of the TRIS or tridiagonal scalar matrix.
!
!    A is integral, therefore det ( A ) is integral, and 
!    det ( A ) * inverse ( A ) is integral.
!
!    A is Toeplitz: constant along diagonals.
!
!    A is symmetric: A' = A.
!
!    Because A is symmetric, it is normal.
!
!    Because A is normal, it is diagonalizable.
!
!    A is persymmetric: A(I,J) = A(N+1-J,N+1-I).
!
!    A is positive definite.
!
!    A is an M matrix.
!
!    A is weakly diagonally dominant, but not strictly diagonally dominant.
!
!    A has an LU factorization A = L * U, without pivoting.
!
!      The matrix L is lower bidiagonal with subdiagonal elements:
!
!        L(I+1,I) = -I/(I+1)
!
!      The matrix U is upper bidiagonal, with diagonal elements
!
!        U(I,I) = (I+1)/I
!
!      and superdiagonal elements which are all -1.
!
!    A has a Cholesky factorization A = L * L', with L lower bidiagonal.
!
!      L(I,I) =    sqrt ( (I+1) / I )
!      L(I,I-1) = -sqrt ( (I-1) / I )
!
!    The eigenvalues are
!
!      LAMBDA(I) = 2 + 2 * COS(I*PI/(N+1))
!                = 4 SIN^2(I*PI/(2*N+2))
!
!    The corresponding eigenvector X(I) has entries
!
!       X(I)(J) = sqrt(2/(N+1)) * sin ( I*J*PI/(N+1) ).
!
!    Simple linear systems:
!
!      x = (1,1,1,...,1,1),   A*x=(1,0,0,...,0,1)
!
!      x = (1,2,3,...,n-1,n), A*x=(0,0,0,...,0,n+1)
!
!    det ( A ) = N + 1.
!
!    The value of the determinant can be seen by induction,
!    and expanding the determinant across the first row:
!
!      det ( A(N) ) = 2 * det ( A(N-1) ) - (-1) * (-1) * det ( A(N-2) )
!                = 2 * N - (N-1)
!                = N + 1
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Gregory, David Karney,
!    A Collection of Matrices for Testing Computational Algorithms,
!    Wiley, 1969,
!    ISBN: 0882756494,
!    LC: QA263.68
!
!    Morris Newman, John Todd,
!    Example A8,
!    The evaluation of matrix inversion programs,
!    Journal of the Society for Industrial and Applied Mathematics,
!    Volume 6, Number 4, pages 466-476, 1958.
!
!    John Todd,
!    Basic Numerical Mathematics,
!    Volume 2: Numerical Algebra,
!    Birkhauser, 1980,
!    ISBN: 0817608117,
!    LC: QA297.T58.
!
!    Joan Westlake,
!    A Handbook of Numerical Matrix Inversion and Solution of 
!    Linear Equations,
!    John Wiley, 1968,
!    ISBN13: 978-0471936756,
!    LC: QA263.W47.
!
!  Parameters:
!
!    Input, integer M, N, the order of the matrix.
!
!    Output, real ( kind = rk ) A(M,N), the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  integer i
  integer j

  do j = 1, n
    do i = 1, m

      if ( j == i - 1 ) then
        a(i,j) = -1.0D+00
      else if ( j == i ) then
        a(i,j) = 2.0D+00
      else if ( j == i + 1 ) then
        a(i,j) = -1.0D+00
      else
        a(i,j) = 0.0D+00
      end if

    end do
  end do

  return
end
subroutine r8ge_mm ( n1, n2, n3, a, b, c )

!*****************************************************************************80
!
!! R8GE_MM multiplies two R8GE matrices.
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
!    12 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, N3, the order of the matrices.
!    N must be positive.
!
!    Input, real ( kind = rk ) A(N1,N2), B(N2,N3), the factors.
!
!    Output, real ( kind = rk ) C(N1,N3), the product.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n1
  integer n2
  integer n3

  real ( kind = rk ) a(n1,n2)
  real ( kind = rk ) b(n2,n3)
  real ( kind = rk ) c(n1,n3)

  c(1:n1,1:n3) = matmul ( a(1:n1,1:n2), b(1:n2,1:n3) )

  return
end
subroutine r8ge_mtv ( m, n, a, x, b )

!*****************************************************************************80
!
!! R8GE_MTV multiplies an R8VEC by an R8GE matrix.
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
!    Input, real ( kind = rk ) A(M,N), the R8GE matrix.
!
!    Input, real ( kind = rk ) X(M), the vector to be multiplied by A.
!
!    Output, real ( kind = rk ) B(N), the product A' * x.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) b(n)
  real ( kind = rk ) x(m)

  b(1:n) = matmul ( transpose ( a(1:m,1:n) ), x(1:m) )

  return
end
subroutine r8ge_mv ( m, n, a, x, b )

!*****************************************************************************80
!
!! R8GE_MV multiplies an R8GE matrix by an R8VEC.
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
!    Input, real ( kind = rk ) A(M,N), the R8GE matrix.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = rk ) B(M), the product A * x.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) b(m)
  real ( kind = rk ) x(n)

  b(1:m) = matmul ( a(1:m,1:n), x(1:n) )

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
subroutine r8ge_to_r8vec ( m, n, a, x )

!*****************************************************************************80
!
!! R8GE_TO_R8VEC copies an R8GE matrix to an R8VEC.
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
!    In C++ and FORTRAN, this routine is not really needed.  In MATLAB,
!    a data item carries its dimensionality implicitly, and so cannot be
!    regarded sometimes as a vector and sometimes as an array.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns in 
!    the array.
!
!    Input, real ( kind = rk ) A(M,N), the array to be copied.
!
!    Output, real ( kind = rk ) X(M*N), the vector.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  integer i
  integer j
  integer k
  real ( kind = rk ) x(m*n)

  k = 0
  do j = 1, n
    do i = 1, m
      k = k + 1
      x(k) = a(i,j)
    end do
  end do

  return
end
subroutine r8ge_np_det ( n, a_lu, det )

!*****************************************************************************80
!
!! R8GE_NP_DET computes the determinant of a matrix factored by R8GE_NP_FA.
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
!    25 March 2004
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
!    Input, real ( kind = rk ) A_LU(N,N), the LU factors from R8GE_NP_FA.
!
!    Output, real ( kind = rk ) DET, the determinant of the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a_lu(n,n)
  real ( kind = rk ) det
  integer i

  det = 1.0D+00

  do i = 1, n
    det = det * a_lu(i,i)
  end do

  return
end
subroutine r8ge_np_fa ( n, a, info )

!*****************************************************************************80
!
!! R8GE_NP_FA factors an R8GE matrix by nonpivoting Gaussian elimination.
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
!    R8GE_NP_FA is a version of the LINPACK routine SGEFA, but uses no
!    pivoting.  It will fail if the matrix is singular, or if any zero
!    pivot is encountered.
!
!    If R8GE_NP_FA successfully factors the matrix, R8GE_NP_SL may be called
!    to solve linear systems involving the matrix.
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
!    Input/output, real ( kind = rk ) A(N,N).
!    On input, A contains the matrix to be factored.
!    On output, A contains information about the factorization,
!    which must be passed unchanged to R8GE_NP_SL for solutions.
!
!    Output, integer INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,n)
  integer info
  integer j
  integer k

  info = 0

  do k = 1, n - 1

    if ( a(k,k) == 0.0D+00 ) then
      info = k
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8GE_NP_FA - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      stop 1
    end if

    a(k+1:n,k) = - a(k+1:n,k) / a(k,k)
    do j = k + 1, n
      a(k+1:n,j) = a(k+1:n,j) + a(k+1:n,k) * a(k,j)
    end do

  end do

  if ( a(n,n) == 0.0D+00 ) then
    info = n
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_NP_FA - Fatal error!'
    write ( *, '(a,i8)' ) '  Zero pivot on step ', info
    stop 1
  end if

  return
end
subroutine r8ge_np_inverse ( n, a )

!*****************************************************************************80
!
!! R8GE_NP_INVERSE computes the inverse of a matrix factored by R8GE_NP_FA.
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
!    12 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix A.
!
!    Input/output, real ( kind = rk ) A(N,N).
!    On input, the factor information computed by R8GE_NP_FA.
!    On output, the inverse matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,n)
  integer j
  integer k
  real ( kind = rk ) temp
  real ( kind = rk ) work(n)
!
!  Compute Inverse(U).
!
  do k = 1, n

    a(k,k) = 1.0D+00 / a(k,k)
    a(1:k-1,k) = - a(1:k-1,k) * a(k,k)

    do j = k + 1, n

      temp = a(k,j)
      a(k,j) = 0.0D+00
      a(1:k,j) = a(1:k,j) + temp * a(1:k,k)

    end do

  end do
!
!  Form Inverse(U) * Inverse(L).
!
  do k = n - 1, 1, -1

    work(k+1:n) = a(k+1:n,k)
    a(k+1:n,k) = 0.0D+00

    do j = k + 1, n
      a(1:n,k) = a(1:n,k) + a(1:n,j) * work(j)
    end do

  end do

  return
end
subroutine r8ge_np_ml ( n, a_lu, x, b, job )

!*****************************************************************************80
!
!! R8GE_NP_ML computes A * x or x * A, for a matrix factored by R8GE_NP_FA.
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
!    The matrix A is assumed to have been factored by R8GE_NP_FA.
!
!    R8GE_NP_ML allows the user to check that the solution of a linear
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
!    Input, real ( kind = rk ) A_LU(N,N), the LU factors from R8GE_NP_FA.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied.
!
!    Output, real ( kind = rk ) B(N), the result of the multiplication.
!
!    Input, integer JOB, determines the multiplication to
!    be carried out:
!    JOB = 0, compute A * x.
!    JOB nonzero, compute A' * X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a_lu(n,n)
  real ( kind = rk ) b(n)
  integer i
  integer j
  integer job
  real ( kind = rk ) x(n)

  b(1:n) = x(1:n)

  if ( job == 0 ) then
!
!  Compute U * X = Y:
!
    do i = 1, n
      b(i) = sum ( a_lu(i,i:n) * b(i:n) )
    end do
!
!  Compute L * Y = B:
!
    do j = n - 1, 1, -1
      b(j+1:n) = b(j+1:n) - a_lu(j+1:n,j) * b(j)
    end do

  else
!
!  Compute L' * X = Y:
!
    do j = 1, n - 1
      b(j) = b(j) - sum ( b(j+1:n) * a_lu(j+1:n,j) )
    end do
!
!  Compute U' * Y = B:
!
    do j = n, 1, -1
      b(j) = sum ( b(1:j) * a_lu(1:j,j) )
    end do

  end if

  return
end
subroutine r8ge_np_sl ( n, a_lu, b, job )

!*****************************************************************************80
!
!! R8GE_NP_SL solves a system factored by R8GE_NP_FA.
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
!    16 October 1998
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
!    Input, real ( kind = rk ) A_LU(N,N), the LU factors from R8GE_NP_FA.
!
!    Input/output, real ( kind = rk ) B(N).
!
!    On input, B contains the right hand side vector B.
!    On output, B contains the solution X.
!
!    Input, integer JOB.
!    If JOB is zero, the routine will solve A * x = b.
!    If JOB is nonzero, the routine will solve A' * x = b.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a_lu(n,n)
  real ( kind = rk ) b(n)
  integer job
  integer k
!
!  Solve A * x = b.
!
  if ( job == 0 ) then

    do k = 1, n - 1
      b(k+1:n) = b(k+1:n) + a_lu(k+1:n,k) * b(k)
    end do

    do k = n, 1, -1
      b(k) = b(k) / a_lu(k,k)
      b(1:k-1) = b(1:k-1) - a_lu(1:k-1,k) * b(k)
    end do
!
!  Solve A' * X = B.
!
  else

    do k = 1, n
      b(k) = ( b(k) - sum ( b(1:k-1) * a_lu(1:k-1,k) ) ) / a_lu(k,k)
    end do

    do k = n - 1, 1, -1
      b(k) = b(k) + sum ( b(k+1:n) * a_lu(k+1:n,k) )
    end do

  end if

  return
end
subroutine r8ge_np_trf ( m, n, a, info )

!*****************************************************************************80
!
!! R8GE_NP_TRF computes the LU factorization of an R8GE matrix.
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
!    R8GE_NP_TRF is a nonpivoting version of R8GE_TRF, and will fail if
!    a zero element is encountered along the diagonal.
!
!    The factorization has the form
!      A = L * U
!    where L is lower triangular with unit diagonal elements (lower
!    trapezoidal if N < M), and U is upper triangular (upper trapezoidal
!    if M < N).
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
!    Input, integer M, the number of rows of the matrix A.  
!    0 <= M.
!
!    Input, integer N, the number of columns of the matrix A.  
!    0 <= N.
!
!    Input/output, real ( kind = rk ) A(M,N).
!    On entry, the M by N matrix to be factored.
!    On exit, the factors L and U from the factorization
!    A = L*U; the unit diagonal elements of L are not stored.
!
!    Output, integer INFO.
!    = 0: successful exit
!    < 0: if INFO = -K, the K-th argument had an illegal value
!    > 0: if INFO = K, U(K,K) is exactly zero. The factorization
!         has been completed, but the factor U is exactly
!         singular, and division by zero will occur if it is used
!         to solve a system of equations.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  integer ii
  integer info
  integer j
!
!  Test the input parameters.
!
  info = 0

  if ( m < 0 ) then
    info = - 1
    return
  else if ( n < 0 ) then
    info = - 2
    return
  end if

  if ( m == 0 .or. n == 0 ) then
    return
  end if

  do j = 1, min ( m, n )
!
!  Compute elements J+1:M of the J-th column.
!
    if ( a(j,j) /= 0.0D+00 ) then
      a(j+1:m,j) = a(j+1:m,j) / a(j,j)
    else if ( info == 0 ) then
      info = j
    end if
!
!  Update the trailing submatrix.
!
    if ( j < min ( m, n ) ) then

      do ii = j + 1, m
        a(ii,j+1:n) = a(ii,j+1:n) - a(ii,j) * a(j,j+1:n)
      end do

    end if

  end do

  return
end
subroutine r8ge_np_trm ( m, n, a, x, b, job )

!*****************************************************************************80
!
!! R8GE_NP_TRM computes A * x or A' * x, for a matrix factored by R8GE_NP_TRF.
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
!    The matrix A is assumed to have been factored by R8GE_NP_TRF.
!
!    R8GE_NP_TRM allows the user to check that the solution of a linear
!    system is correct, without having to save an unfactored copy
!    of the matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    24 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
!    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
!    Sven Hammarling, Alan McKenney, Danny Sorensen,
!    LAPACK User's Guide,
!    Second Edition,
!    SIAM, 1995.
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns in 
!    the matrix.  M and N must be positive.
!
!    Input, real ( kind = rk ) A(M,N), the M by N matrix factors computed
!    by R8GE_NP_TRF.
!
!    Input, real ( kind = rk ) X(*), the vector to be multiplied.
!    If JOB is 0, X must have dimension N.
!    If JOB is nonzero, X must have dimension M.
!
!    Output, real ( kind = rk ) B(*), the result of the multiplication.
!    If JOB is 0, B must have dimension M.
!    If JOB is nonzero, B must have dimension N.
!
!    Input, integer JOB, determines the multiplication to
!    be carried out:
!    JOB = 0, compute A * x.
!    JOB nonzero, compute A' * X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) b(*)
  integer i
  integer job
  real ( kind = rk ) x(*)

  if ( job == 0 ) then

    b(1:m) = 0.0D+00
!
!  Compute U * X = Y:
!
    do i = 1, m
      b(i) = sum ( a(i,i:n) * x(i:n) )
    end do
!
!  Compute L * Y = B:
!
    do i = min ( m, n + 1 ), 2, -1
      b(i) = b(i) + sum ( a(i,1:i-1) * b(1:i-1) )
    end do

  else

    b(1:n) = 0.0D+00
!
!  Compute L' * X = Y:
!
    do i = 1, min ( m, n )
      b(i) = x(i) + sum ( a(i+1:m,i) * x(i+1:m) )
    end do
!
!  Compute U' * Y = B:
!
    do i = min ( m, n ), 1, -1
      b(i) = sum ( a(1:i,i) * b(1:i) )
    end do

  end if

  return
end
subroutine r8ge_np_trs ( n, nrhs, trans, a, b, info )

!*****************************************************************************80
!
!! R8GE_NP_TRS solves a system of linear equations factored by R8GE_NP_TRF.
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
!    R8GE_NP_TRS is a nonpivoting version of R8GE_TRS.
!
!    R8GE_TRS solves a system of linear equations
!      A * x = b  or  A' * X = B
!    with a general N by N matrix A using the LU factorization computed
!    by R8GE_NP_TRF.
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
!  Reference:
!
!    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
!    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
!    Sven Hammarling, Alan McKenney, Danny Sorensen,
!    LAPACK User's Guide,
!    Second Edition,
!    SIAM, 1995.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix A.  0 <= N.
!
!    Input, integer NRHS, the number of right hand sides.  
!    0 <= NRHS.
!
!    Input, character TRANS, pecifies the form of the system of equations:
!    'N':  A * x = b  (No transpose)
!    'T':  A'* X = B  (Transpose)
!    'C':  A'* X = B  (Conjugate transpose = Transpose)
!
!    Input, real ( kind = rk ) A(N,N), the factors L and U from the factorization
!    A = L*U as computed by R8GE_NP_TRF.
!
!    Input/output, real ( kind = rk ) B(N,NRHS).
!    On entry, the right hand side matrix B.
!    On exit, the solution matrix X.
!
!    Output, integer INFO
!    = 0:  successful exit
!    < 0:  if INFO = -I, the I-th argument had an illegal value.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nrhs

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) b(n,nrhs)
  integer info
  integer j
  integer k
  character trans

  info = 0

  if ( trans /= 'n' .and. trans /= 'N' .and. &
       trans /= 't' .and. trans /= 'T' .and. &
       trans /= 'c' .and. trans /= 'C' ) then
    info = - 1
    return
  else if ( n < 0 ) then
    info = - 2
    return
  else if ( nrhs < 0 ) then
    info = - 3
    return
  end if

  if ( n == 0 .or. nrhs == 0 ) then
    return
  end if

  if ( trans == 'n' .or. trans == 'N' ) then
!
!  Solve L * x = b, overwriting b with x.
!
    do k = 1, nrhs
      do j = 1, n - 1
        b(j+1:n,k) = b(j+1:n,k) - a(j+1:n,j) * b(j,k)
      end do
    end do
!
!  Solve U * x = b, overwriting b with x.
!
    do k = 1, nrhs
      do j = n, 1, -1
        b(j,k) = b(j,k) / a(j,j)
        b(1:j-1,k) = b(1:j-1,k) - a(1:j-1,j) * b(j,k)
      end do
    end do

  else
!
!  Solve U' * x = b, overwriting b with x.
!
    do k = 1, nrhs
      do j = 1, n
        b(j,k) = b(j,k) / a(j,j)
        b(j+1:n,k) = b(j+1:n,k) - a(j,j+1:n) * b(j,k)
      end do
    end do
!
!  Solve L' * x = b, overwriting b with x.
!
    do k = 1, nrhs
      do j = n, 2, -1
        b(1:j-1,k) = b(1:j-1,k) - a(j,1:j-1) * b(j,k)
      end do
    end do

  end if

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
subroutine r8vec_print_some ( n, a, max_print, title )

!*****************************************************************************80
!
!! R8VEC_PRINT_SOME prints "some" of an R8VEC.
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vector, is no more than MAX_PRINT, then
!    the entire vector is printed, one entry per line.
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
!    16 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries of the vector.
!
!    Input, real ( kind = rk ) A(N), the vector to be printed.
!
!    Input, integer MAX_PRINT, the maximum number of lines 
!    to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  integer i
  integer max_print
  character ( len = * ) title

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

    if ( all ( a(1:n) == aint ( a(1:n) ) ) ) then
      do i = 1, n
        write ( *, '(i8,2x,i8)' ) i, int ( a(i) )
      end do
    else if ( all ( abs ( a(1:n) ) < 1000000.0D+00 ) ) then
      do i = 1, n
        write ( *, '(i8,2x,f14.6)' ) i, a(i)
      end do
    else
      do i = 1, n
        write ( *, '(i8,2x,g14.6)' ) i, a(i)
      end do
    end if

  else if ( 3 <= max_print ) then

    if ( all ( a(1:max_print-2) == aint ( a(1:max_print-2) ) ) ) then
      do i = 1, max_print - 2
        write ( *, '(i8,2x,i8)' ) i, int ( a(i) )
      end do
    else if ( all ( abs ( a(1:max_print-2) ) < 1000000.0D+00 ) ) then
      do i = 1, max_print - 2
        write ( *, '(i8,2x,f14.6)' ) i, a(i)
      end do
    else
      do i = 1, max_print - 2
        write ( *, '(i8,2x,g14.6)' ) i, a(i)
      end do
    end if

    write ( *, '(a)' ) '......  ..............'
    i = n

    if ( a(i) == real ( int ( a(i) ), kind = rk ) ) then
      write ( *, '(i8,2x,i8)' ) i, int ( a(i) )
    else if (  abs ( a(i) ) < 1000000.0D+00 ) then
      write ( *, '(i8,2x,f14.6)' ) i, a(i)
    else
      write ( *, '(i8,2x,g14.6)' ) i, a(i)
    end if

  else

    if ( all ( a(1:max_print-1) == aint ( a(1:max_print-1) ) ) ) then
      do i = 1, max_print - 1
        write ( *, '(i8,2x,i8)' ) i, int ( a(i) )
      end do
    else if ( all ( abs ( a(1:max_print-1) ) < 1000000.0D+00 ) ) then
      do i = 1, max_print - 1
        write ( *, '(i8,2x,f14.6)' ) i, a(i)
      end do
    else
      do i = 1, max_print - 1
        write ( *, '(i8,2x,g14.6)' ) i, a(i)
      end do
    end if

    i = max_print

    if ( a(i) == aint ( a(i) ) ) then
      write ( *, '(i8,2x,i8,a)' ) i, int ( a(i) ), '...more entries...'
    else if (  abs ( a(i) ) < 1000000.0D+00 ) then
      write ( *, '(i8,2x,f14.6,a)' ) i, a(i), '...more entries...'
    else
      write ( *, '(i8,2x,g14.6,a)' ) i, a(i), '...more entries...'
    end if

  end if

  return
end
subroutine r8vec_to_r8ge ( m, n, x, a )

!*****************************************************************************80
!
!! R8VEC_TO_R8GE copies an R8VEC into an R8GE matrix.
!
!  Discussion:
!
!    In C++ and FORTRAN, this routine is not really needed.  In MATLAB,
!    a data item carries its dimensionality implicitly, and so cannot be
!    regarded sometimes as a vector and sometimes as an array.
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
!    16 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns 
!    in the array.
!
!    Input, real ( kind = rk ) X(M*N), the vector to be copied into the array.
!
!    Output, real ( kind = rk ) A(M,N), the array.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  integer i
  integer j
  integer k
  real ( kind = rk ) x(m*n)

  k = 0
  do j = 1, n
    do i = 1, m
      k = k + 1
      a(i,j) = x(k)
    end do
  end do

  return
end
subroutine r8vec2_print ( n, a1, a2, title )

!*****************************************************************************80
!
!! R8VEC2_PRINT prints an R8VEC2.
!
!  Discussion:
!
!    An R8VEC2 is a dataset consisting of N pairs of R8's, stored
!    as two separate vectors A1 and A2.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 December 2004
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
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a1(n)
  real ( kind = rk ) a2(n)
  integer i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i4,2x,g14.6,2x,g14.6)' ) i, a1(i), a2(i)
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
subroutine r8vm_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8VM_PRINT prints an R8VM matrix.
!
!  Discussion:
!
!    The R8VM storage format is used for an M by N Vandermonde matrix.
!    An M by N Vandermonde matrix is defined by the values in its second
!    row, which will be written here as X(1:N).  The matrix has a first 
!    row of 1's, a second row equal to X(1:N), a third row whose entries
!    are the squares of the X values, up to the M-th row whose entries
!    are the (M-1)th powers of the X values.  The matrix can be stored
!    compactly by listing just the values X(1:N).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    29 September 2003
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
!    Input, real ( kind = rk ) A(N), the R8VM matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  integer m
  character ( len = * ) title

  call r8vm_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8vm_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8VM_PRINT_SOME prints some of an R8VM matrix.
!
!  Discussion:
!
!    The R8VM storage format is used for an M by N Vandermonde matrix.
!    An M by N Vandermonde matrix is defined by the values in its second
!    row, which will be written here as X(1:N).  The matrix has a first 
!    row of 1's, a second row equal to X(1:N), a third row whose entries
!    are the squares of the X values, up to the M-th row whose entries
!    are the (M-1)th powers of the X values.  The matrix can be stored
!    compactly by listing just the values X(1:N).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    29 September 2003
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
!    Input, real ( kind = rk ) A(N), the R8VM matrix.
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
  integer m
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

        if ( i == 1 ) then
          aij = 1.0D+00
        else
          aij = a(j)**(i-1)
        end if

        write ( ctemp(j2), '(g14.6)' ) aij

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine sort_heap_external ( n, indx, i, j, isgn )

!*****************************************************************************80
!
!! SORT_HEAP_EXTERNAL externally sorts a list of items into linear order.
!
!  Discussion:
!
!    The actual list of data is not passed to the routine.  Hence this
!    routine may be used to sort integers, reals, numbers, names,
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
!    12 November 2000
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
!      * set ISGN = -1 if I precedes J, ISGN = +1 if J precedes I;
!      * call again.
!
!      equal to 0, the sorting is done.
!
!    Output, integer I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ISGN, results of comparison of elements 
!    I and J.  (Used only when the previous call returned INDX less than 0).
!    ISGN <= 0 means I precedes J;
!    ISGN => 0 means J precedes I.
!
  implicit none

  integer i
  integer indx
  integer isgn
  integer j
  integer, save :: k = 0
  integer, save :: k1 = 0
  integer n
  integer, save :: n1 = 0
!
!  INDX = 0: This is the first call.
!
  if ( indx == 0 ) then

    n1 = n
    k = n / 2
    k1 = k
!
!  INDX < 0: The user is returning the results of a comparison.
!
  else if ( indx < 0 ) then

    if ( indx == -2 ) then

      if ( isgn < 0 ) then
        i = i + 1
      end if

      j = k1
      k1 = i
      indx = - 1
      return

    end if

    if ( 0 < isgn ) then
      indx = 2
      return
    end if

    if ( k <= 1 ) then

      if ( n1 == 1 ) then
        indx = 0
      else
        i = n1
        n1 = n1 - 1
        j = 1
        indx = 1
      end if

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

    i = 2 * k1

    if ( i == n1 ) then
      j = k1
      k1 = i
      indx = - 1
      return
    else if ( i <= n1 ) then
      j = i + 1
      indx = - 2
      return
    end if

    if ( k <= 1 ) then
      exit
    end if

    k = k - 1
    k1 = k

  end do

  if ( n1 == 1 ) then
    indx = 0
  else
    i = n1
    n1 = n1 - 1
    j = 1
    indx = 1
  end if

  return
end

