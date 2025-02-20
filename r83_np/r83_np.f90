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
subroutine r83_dif2 ( m, n, a )

!*****************************************************************************80
!
!! R83_DIF2 returns the DIF2 matrix in R83 format.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:min(M+1,N)).
!    The diagonal in entries (2,1:min(M,N)).
!    The subdiagonal in (3,min(M-1,N)).
!    R8GE A(I,J) = R83 A(I-J+2,J).
!
!  Example:
!
!    An R83 matrix of order 3x5 would be stored:
!
!       *  A12 A23 A34  *
!      A11 A22 A33  *   *
!      A21 A32  *   *   *
!
!    An R83 matrix of order 5x5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!    An R83 matrix of order 5x3 would be stored:
!
!       *  A12 A23
!      A11 A22 A33
!      A21 A32 A43
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
!    05 September 2015
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
!    Output, real ( kind = rk ) A(3,N), the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(3,n)
  integer i
  integer j

  a(1:3,1:n) = 0.0D+00

  do j = 1, n
    do i = max ( 1, j - 1 ), min ( m, j + 1 )
      if ( i - j + 2 == 1 ) then
        a(i-j+2,j) = -1.0D+00
      else if ( i - j + 2 == 2 ) then
        a(i-j+2,j) = +2.0D+00
      else if ( i - j + 2 == 3 ) then
        a(i-j+2,j) = -1.0D+00
      end if
    end do
  end do
  
  return
end
subroutine r83_mtv ( m, n, a, x, b )

!*****************************************************************************80
!
!! R83_MTV computes A'*x=b, where A is an R83 matrix.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:min(M+1,N)).
!    The diagonal in entries (2,1:min(M,N)).
!    The subdiagonal in (3,min(M-1,N)).
!    R8GE A(I,J) = R83 A(I-J+2,J).
!
!  Example:
!
!    An R83 matrix of order 3x5 would be stored:
!
!       *  A12 A23 A34  *
!      A11 A22 A33  *   *
!      A21 A32  *   *   *
!
!    An R83 matrix of order 5x5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!    An R83 matrix of order 5x3 would be stored:
!
!       *  A12 A23
!      A11 A22 A33
!      A21 A32 A43
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 September 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the order of the linear system.
!
!    Input, real ( kind = rk ) A(3,N), the MxN R83 matrix.
!
!    Input, real ( kind = rk ) X(M), the vector to be multiplied by A'.
!
!    Output, real ( kind = rk ) B(N), the product A' * x.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(3,n)
  real ( kind = rk ) b(n)
  integer i
  integer j
  real ( kind = rk ) x(m)

  b(1:n) = 0.0D+00
!
!  Find each nonzero A(I,J), multiply by X(I), add to B(J).
!
  do j = 1, n
    do i = max ( 1, j - 1 ), min ( m, j + 1 )
      b(j) = b(j) + x(i) * a(i-j+2,j)
    end do
  end do

  return
end
subroutine r83_mv ( m, n, a, x, b )

!*****************************************************************************80
!
!! R83_MV multiplies an R83 matrix times an R8VEC.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:min(M+1,N)).
!    The diagonal in entries (2,1:min(M,N)).
!    The subdiagonal in (3,min(M-1,N)).
!    R8GE A(I,J) = R83 A(I-J+2,J).
!
!  Example:
!
!    An R83 matrix of order 3x5 would be stored:
!
!       *  A12 A23 A34  *
!      A11 A22 A33  *   *
!      A21 A32  *   *   *
!
!    An R83 matrix of order 5x5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!    An R83 matrix of order 5x3 would be stored:
!
!       *  A12 A23
!      A11 A22 A33
!      A21 A32 A43
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 September 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = rk ) A(3,N), the R83 matrix.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = rk ) B(M), the product A * x.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(3,n)
  real ( kind = rk ) b(m)
  integer i
  integer j
  real ( kind = rk ) x(n)

  b(1:m) = 0.0D+00

  do j = 1, n
    do i = max ( 1, j - 1 ), min ( m, j + 1 )
      b(i) = b(i) + a(i-j+2,j) * x(j)
    end do
  end do

  return
end
subroutine r83_print ( m, n, a, title )

!*****************************************************************************80
!
!! R83_PRINT prints an R83 matrix.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:min(M+1,N)).
!    The diagonal in entries (2,1:min(M,N)).
!    The subdiagonal in (3,min(M-1,N)).
!    R8GE A(I,J) = R83 A(I-J+2,J).
!
!  Example:
!
!    An R83 matrix of order 3x5 would be stored:
!
!       *  A12 A23 A34  *
!      A11 A22 A33  *   *
!      A21 A32  *   *   *
!
!    An R83 matrix of order 5x5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!    An R83 matrix of order 5x3 would be stored:
!
!       *  A12 A23
!      A11 A22 A33
!      A21 A32 A43
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    31 August 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the order of the matrix.
!    N must be positive.
!
!    Input, real ( kind = rk ) A(3,N), the R83 matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(3,n)
  integer m
  character ( len = * ) title

  call r83_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r83_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R83_PRINT_SOME prints some of an R83 matrix.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:min(M+1,N)).
!    The diagonal in entries (2,1:min(M,N)).
!    The subdiagonal in (3,min(M-1,N)).
!    R8GE A(I,J) = R83 A(I-J+2,J).
!
!  Example:
!
!    An R83 matrix of order 3x5 would be stored:
!
!       *  A12 A23 A34  *
!      A11 A22 A33  *   *
!      A21 A32  *   *   *
!
!    An R83 matrix of order 5x5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!    An R83 matrix of order 5x3 would be stored:
!
!       *  A12 A23
!      A11 A22 A33
!      A21 A32 A43
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 September 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the order of the matrix.
!    N must be positive.
!
!    Input, real ( kind = rk ) A(3,N), the R83 matrix.
!
!    Input, integer ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column, to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: incx = 5
  integer n

  real ( kind = rk ) a(3,n)
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

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if
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
    i2lo = max ( i2lo, j2lo - 1 )
    i2hi = min ( ihi, m )
    i2hi = min ( i2hi, j2hi + 1 )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( i - j + 2 < 1 .or. 3 < i - j + 2 ) then
          ctemp(j2) = '              '
        else
          write ( ctemp(j2), '(g14.6)' ) a(i-j+2,j)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r83_random ( m, n, seed, a )

!*****************************************************************************80
!
!! R83_RANDOM randomizes an R83 matrix.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:min(M+1,N)).
!    The diagonal in entries (2,1:min(M,N)).
!    The subdiagonal in (3,min(M-1,N)).
!    R8GE A(I,J) = R83 A(I-J+2,J).
!
!  Example:
!
!    An R83 matrix of order 3x5 would be stored:
!
!       *  A12 A23 A34  *
!      A11 A22 A33  *   *
!      A21 A32  *   *   *
!
!    An R83 matrix of order 5x5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!    An R83 matrix of order 5x3 would be stored:
!
!       *  A12 A23
!      A11 A22 A33
!      A21 A32 A43
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 September 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the linear system.
!
!    Input/output, integer SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = rk ) A(3,N), the R83 matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(3,n)
  integer i
  integer j
  integer m
  real ( kind = rk ) r8_uniform_01
  integer seed

  a(1:3,1:n) = 0.0D+00
 
  do j = 1, n
    do i = max ( 1, j - 1 ), min ( m, j + 1 )
      a(i-j+2,j) = r8_uniform_01 ( seed )
    end do
  end do

  return
end
subroutine r83_np_det ( n, a_lu, det )

!*****************************************************************************80
!
!! R83_NP_DET returns the determinant of an R83 system factored by R83_NP_FA.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:min(M+1,N)).
!    The diagonal in entries (2,1:min(M,N)).
!    The subdiagonal in (3,min(M-1,N)).
!
!  Example:
!
!    An R83 matrix of order 3x5 would be stored:
!
!       *  A12 A23 A34  *
!      A11 A22 A33  *   *
!      A21 A32  *   *   *
!
!    An R83 matrix of order 5x5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!    An R83 matrix of order 5x3 would be stored:
!
!       *  A12 A23
!      A11 A22 A33
!      A21 A32 A43
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
!    N must be at least 2.
!
!    Input, real ( kind = rk ) A_LU(3,N), the LU factors computed by R83_NP_FA.
!
!    Output, real ( kind = rk ) DET, the determinant of the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a_lu(3,n)
  real ( kind = rk ) det

  det = product ( a_lu(2,1:n) )

  return
end
subroutine r83_np_fa ( n, a, info )

!*****************************************************************************80
!
!! R83_NP_FA factors an R83 matrix without pivoting.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:min(M+1,N)).
!    The diagonal in entries (2,1:min(M,N)).
!    The subdiagonal in (3,min(M-1,N)).
!
!    Because this routine does not use pivoting, it can fail even when
!    the matrix is not singular, and it is liable to make larger
!    errors.
!
!    R83_NP_FA and R83_NP_SL may be preferable to the corresponding
!    LINPACK routine SGTSL for tridiagonal systems, which factors and solves
!    in one step, and does not save the factorization.
!
!  Example:
!
!    An R83 matrix of order 3x5 would be stored:
!
!       *  A12 A23 A34  *
!      A11 A22 A33  *   *
!      A21 A32  *   *   *
!
!    An R83 matrix of order 5x5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!    An R83 matrix of order 5x3 would be stored:
!
!       *  A12 A23
!      A11 A22 A33
!      A21 A32 A43
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 2.
!
!    Input/output, real ( kind = rk ) A(3,N).
!    On input, the tridiagonal matrix.  On output, factorization information.
!
!    Output, integer INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(3,n)
  integer i
  integer info

  info = 0

  do i = 1, n - 1

    if ( a(2,i) == 0.0D+00 ) then
      info = i
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R83_NP_FA - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      stop 1
    end if
!
!  Store the multiplier in L.
!
    a(3,i) = a(3,i) / a(2,i)
!
!  Modify the diagonal entry in the next column.
!
    a(2,i+1) = a(2,i+1) - a(3,i) * a(1,i+1)

  end do

  if ( a(2,n) == 0.0D+00 ) then
    info = n
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R83_NP_FA - Fatal error!'
    write ( *, '(a,i8)' ) '  Zero pivot on step ', info
    stop 1
  end if

  return
end
subroutine r83_np_fs ( n, a, b, x )

!*****************************************************************************80
!
!! R83_NP_FS factors and solves an R83 system.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:min(M+1,N)).
!    The diagonal in entries (2,1:min(M,N)).
!    The subdiagonal in (3,min(M-1,N)).
!
!    This algorithm requires that each diagonal entry be nonzero.
!    It does not use pivoting, and so can fail on systems that
!    are actually nonsingular.
!
!  Example:
!
!    An R83 matrix of order 3x5 would be stored:
!
!       *  A12 A23 A34  *
!      A11 A22 A33  *   *
!      A21 A32  *   *   *
!
!    An R83 matrix of order 5x5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!    An R83 matrix of order 5x3 would be stored:
!
!       *  A12 A23
!      A11 A22 A33
!      A21 A32 A43
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the linear system.
!
!    Input/output, real ( kind = rk ) A(3,N).
!    On input, the tridiagonal matrix.
!    On output, the data in these vectors has been overwritten
!    by factorization information.
!
!    Input, real ( kind = rk ) B(N), the right hand side of the linear system.
!
!    Output, real ( kind = rk ) X(N), the solution of the linear system.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(3,n)
  real ( kind = rk ) b(n)
  integer i
  real ( kind = rk ) x(n)
  real ( kind = rk ) xmult
!
!  The diagonal entries can't be zero.
!
  do i = 1, n
    if ( a(2,i) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R83_NP_FS - Fatal error!'
      write ( *, '(a,i8,a)' ) '  A(2,', i, ') = 0.'
      stop 1
    end if
  end do

  x(1:n) = b(1:n)

  do i = 2, n
    xmult = a(3,i-1) / a(2,i-1)
    a(2,i) = a(2,i) - xmult * a(1,i)
    x(i)   = x(i)   - xmult * x(i-1)
  end do

  x(n) = x(n) / a(2,n)
  do i = n - 1, 1, -1
    x(i) = ( x(i) - a(1,i+1) * x(i+1) ) / a(2,i)
  end do

  return
end
subroutine r83_np_fss ( n, a, nb, b, x )

!*****************************************************************************80
!
!! R83_NP_FSS factors and solves multiple R83 systems.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:min(M+1,N)).
!    The diagonal in entries (2,1:min(M,N)).
!    The subdiagonal in (3,min(M-1,N)).
!
!    This algorithm requires that each diagonal entry be nonzero.
!    It does not use pivoting, and so can fail on systems that
!    are actually nonsingular.
!
!  Example:
!
!    An R83 matrix of order 3x5 would be stored:
!
!       *  A12 A23 A34  *
!      A11 A22 A33  *   *
!      A21 A32  *   *   *
!
!    An R83 matrix of order 5x5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!    An R83 matrix of order 5x3 would be stored:
!
!       *  A12 A23
!      A11 A22 A33
!      A21 A32 A43
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the linear system.
!
!    Input/output, real ( kind = rk ) A(3,N).
!    On input, the tridiagonal matrix.
!    On output, the data in these vectors has been overwritten
!    by factorization information.
!
!    Input, integer NB, the number of right hand sides.
!
!    Input, real ( kind = rk ) B(N,NB), the right hand side of the linear system.
!
!    Output, real ( kind = rk ) X(N,NB), the solution of the linear system.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nb

  real ( kind = rk ) a(3,n)
  real ( kind = rk ) b(n,nb)
  integer i
  real ( kind = rk ) x(n,nb)
  real ( kind = rk ) xmult
!
!  The diagonal entries can't be zero.
!
  do i = 1, n
    if ( a(2,i) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R83_NP_FSS - Fatal error!'
      write ( *, '(a,i8,a)' ) '  A(2,', i, ') = 0.'
      stop 1
    end if
  end do

  x(1:n,1:nb) = b(1:n,1:nb)

  do i = 2, n
    xmult = a(3,i-1) / a(2,i-1)
    a(2,i) = a(2,i) - xmult * a(1,i)
    x(i,1:nb)   = x(i,1:nb)   - xmult * x(i-1,1:nb)
  end do

  x(n,1:nb) = x(n,1:nb) / a(2,n)
  do i = n - 1, 1, -1
    x(i,1:nb) = ( x(i,1:nb) - a(1,i+1) * x(i+1,1:nb) ) / a(2,i)
  end do

  return
end
subroutine r83_np_ml ( n, a_lu, x, b, job )

!*****************************************************************************80
!
!! R83_NP_ML computes A * x or x * A, where A has been factored by R83_NP_FA.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:min(M+1,N)).
!    The diagonal in entries (2,1:min(M,N)).
!    The subdiagonal in (3,min(M-1,N)).
!
!  Example:
!
!    An R83 matrix of order 3x5 would be stored:
!
!       *  A12 A23 A34  *
!      A11 A22 A33  *   *
!      A21 A32  *   *   *
!
!    An R83 matrix of order 5x5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!    An R83 matrix of order 5x3 would be stored:
!
!       *  A12 A23
!      A11 A22 A33
!      A21 A32 A43
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 2.
!
!    Input, real ( kind = rk ) A_LU(3,N), the LU factors from R83_FA.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = rk ) B(N), the product A*x or A'*x.
!
!    Input, integer JOB, specifies the product to find.
!    0, compute A * x.
!    nonzero, compute A' * x.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a_lu(3,n)
  real ( kind = rk ) b(n)
  integer i
  integer job
  real ( kind = rk ) x(n)

  b(1:n) = x(1:n)

  if ( job == 0 ) then
!
!  Compute X := U * X
!
    do i = 1, n

      b(i) = a_lu(2,i) * b(i)

      if ( i < n ) then
        b(i) = b(i) + a_lu(1,i+1) * b(i+1)
      end if

    end do
!
!  Compute X: = L * X.
!
    do i = n, 2, -1
      b(i) = b(i) + a_lu(3,i-1) * b(i-1)
    end do

  else
!
!  Compute X: = L' * X.
!
    do i = 1, n - 1
      b(i) = b(i) + a_lu(3,i) * b(i+1)
    end do
!
!  Compute X: = U' * X.
!
    do i = n, 2, -1
      b(i) = a_lu(2,i) * b(i)
      b(i) = b(i) + a_lu(1,i) * b(i-1)
    end do
    b(1) = a_lu(2,1) * b(1)

  end if

  return
end
subroutine r83_np_sl ( n, a_lu, b, job )

!*****************************************************************************80
!
!! R83_NP_SL solves an R83 system factored by R83_NP_FA.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:min(M+1,N)).
!    The diagonal in entries (2,1:min(M,N)).
!    The subdiagonal in (3,min(M-1,N)).
!
!  Example:
!
!    An R83 matrix of order 3x5 would be stored:
!
!       *  A12 A23 A34  *
!      A11 A22 A33  *   *
!      A21 A32  *   *   *
!
!    An R83 matrix of order 5x5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!    An R83 matrix of order 5x3 would be stored:
!
!       *  A12 A23
!      A11 A22 A33
!      A21 A32 A43
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 2.
!
!    Input, real ( kind = rk ) A_LU(3,N), the LU factors from R83_NP_FA.
!
!    Input/output, real ( kind = rk ) B(N).
!    On input, B contains the right hand side of the linear system.
!    On output, B contains the solution of the linear system.
!
!    Input, integer JOB, specifies the system to solve.
!    0, solve A * x = b.
!    nonzero, solve A' * x = b.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a_lu(3,n)
  real ( kind = rk ) b(n)
  integer i
  integer job

  if ( job == 0 ) then
!
!  Solve L * Y = B.
!
    do i = 2, n
      b(i) = b(i) - a_lu(3,i-1) * b(i-1)
    end do
!
!  Solve U * X = Y.
!
    do i = n, 1, -1
      b(i) = b(i) / a_lu(2,i)
      if ( 1 < i ) then
        b(i-1) = b(i-1) - a_lu(1,i) * b(i)
      end if
    end do

  else
!
!  Solve U' * Y = B
!
    do i = 1, n
      b(i) = b(i) / a_lu(2,i)
      if ( i < n ) then
        b(i+1) = b(i+1) - a_lu(1,i+1) * b(i)
      end if
    end do
!
!  Solve L' * X = Y.
!
    do i = n - 1, 1, -1
      b(i) = b(i) - a_lu(3,i) * b(i+1)
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
subroutine r8vec_indicator0 ( n, a )

!*****************************************************************************80
!
!! R8VEC_INDICATOR0 sets an R8VEC to the indicator vector (0,1,2,...).
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
!    27 September 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements of A.
!
!    Output, real ( kind = rk ) A(N), the array.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  integer i

  do i = 1, n
    a(i) = real ( i - 1, kind = rk )
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
subroutine r8vec2_print_some ( n, x1, x2, max_print, title )

!*****************************************************************************80
!
!! R8VEC2_PRINT_SOME prints "some" of an R8VEC2.
!
!  Discussion:
!
!    An R8VEC2 is a dataset consisting of N pairs of R8's, stored
!    as two separate vectors A1 and A2.
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
!    10 September 2009
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
      write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) i, x1(i), x2(i)
    end do

  else if ( 3 <= max_print ) then

    do i = 1, max_print - 2
      write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) i, x1(i), x2(i)
    end do
    write ( *, '(a)' ) '  ......  ..............  ..............'
    i = n
    write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) i, x1(i), x2(i)

  else

    do i = 1, max_print - 1
      write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) i, x1(i), x2(i)
    end do
    i = max_print
    write ( *, '(2x,i8,2x,g14.6,2x,g14.6,2x,a)' ) i, x1(i), x2(i), &
      '...more entries...'

  end if

  return
end

