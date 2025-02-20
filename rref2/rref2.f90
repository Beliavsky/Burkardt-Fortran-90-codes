function is_rref ( m, n, A )

!*****************************************************************************80
!
!! is_rref() determines if a matrix is in reduced row echelon format.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 December 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    A(m,n): a matrix.
!
!  Output:
!
!    is_ref: True if A is in reduced row echelon form.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk8 ) A(m,n)
  integer c
  integer i
  logical is_rref
  integer j
  integer p
  integer r

  c = 0

  do r = 1, m
!
!  Increment the first legal column for next pivot.
!
    c = c + 1
!
!  Search for pivot p in this row.
!  If none, set p = n + 1.
!
    p = n + 1
    do j = 1, n
      if ( A(r,j) /= 0.0 ) then
        p = j
        exit
      end if
    end do
!
!  If p == n + 1, go to next row.
!
    if ( p == n + 1 ) then
      cycle
    end if
!
!  If p is too early, fail.
!
    if ( p < c ) then
      is_rref = .false.
      return
    end if
!
!  Accept p as new c.
!
    c = p
!
!  If A(r,c) is not 1, fail
!
    if ( A(r,c) /= 1.0 ) then
      is_rref = .false.
      return
    end if
!
!  If A(r,c) is not the only nonzero in column c, fail
!
    do i = 1, m
      if ( i /= r ) then
        if ( A(i,c) /= 0.0D+00 ) then
          is_rref = .false.
          return
        end if
      end if
    end do

  end do

  is_rref = .true.

  return
end
subroutine rref_compute ( m, n, A, A_RREF, a_cols )

!*****************************************************************************80
!
!! rref_compute() computes the reduced row echelon form of a matrix.
!
!  Discussion:
!
!    A rectangular matrix is in row reduced echelon form if:
!
!    * The leading nonzero entry in each row has the value 1.
!
!    * All entries are zero above and below the leading nonzero entry 
!      in each row.
!
!    * The leading nonzero in each row occurs in a column to
!      the right of the leading nonzero in the previous row.
!
!    * Rows which are entirely zero occur last.
!
!  Example:
!
!    M = 4, N = 7
!
!    Matrix A:
!
!     1    3    0    2    6    3    1
!    -2   -6    0   -2   -8    3    1
!     3    9    0    0    6    6    2
!    -1   -3    0    1    0    9    3
!
!    RREF(A):
!
!     1    3    0    0    2    0    0
!     0    0    0    1    2    0    0
!     0    0    0    0    0    1   1/3
!     0    0    0    0    0    0    0
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 December 2024
!
!  Author:
!
!    Original Python version by ChatGPT.
!    This version by John Burkardt.
!
!  Reference:
!
!    Charles Cullen,
!    An Introduction to Numerical Linear Algebra,
!    PWS Publishing Company, 1994,
!    ISBN: 978-0534936903,
!    LC: QA185.D37.C85.
!
!  Input:
!
!    integer m, n: the number of rows and columns in the matrix A.
!
!    real A(M,N), the matrix to be analyzed. 
!
!  Output:
!
!    real A_RREF(M,N), the reduced row echelon form of the matrix.
!
!    integer a_cols(n): the columns for which a pivot entry was found.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk8 ) A(m,n)
  integer a_cols(n)
  real ( kind = rk8 ) A_RREF(m,n)
  integer c
  integer i
  integer j
  integer p
  real ( kind = rk8 ) pval
  integer r
  integer rank
  real ( kind = rk8 ) temp
  real ( kind = rk8 ) tol
!
!  Work on a copy of the original matrix.
!
  A_RREF(1:m,1:n) = A(1:m,1:n)
!
!  Initialize rank to 0.
!
  rank = 0
!
!  Initialize indices of independent columns to 0.
!
  a_cols(1:n) = 0
!
!  Set a tolerance for the size of a pivot value.
!
  tol = sqrt ( epsilon ( tol ) )
!
!  Start the pivot search in row R = 1.
!
  r = 1
!
!  Seek a pivot for column C.
!
  do c = 1, n
!
!  Exit if we have run out of rows to examine.
!
    if ( m < r ) then
      exit
    end if
!
!  Find row index P of maximum element in subvector A(R:M,C).
!
    pval = -1.0
    p = -1
    do i = r, m
      if ( pval <= abs ( A_RREF(i,c) ) ) then
        p = i
        pval = abs ( A_RREF(i,c) )
      end if
    end do
!
!  A pivot was not found.  A(r:m,c) was all tiny.
!
    if ( pval <= tol ) then
      A_RREF(p,c) = 0.0
      cycle
    end if
!
!  A pivot was found.
!
    rank = rank + 1
    a_cols(rank) = c
!
!  Swap rows R and P.
!
    do j = c, n
      temp        = A_RREF(r,j)
      A_RREF(r,j) = A_RREF(p,j)
      A_RREF(p,j) = temp
    end do
!
!  Normalize row R so that A(r,c) = 1.
!
    A_RREF(r,c:n) = A_RREF(r,c:n) / A_RREF(r,c)
!
!  Eliminate nonzeros in column C using multiples of row R, except for A(R,C).
!
    do i = 1, m
      if ( i /= r ) then
        A_RREF(i,c:n) = A_RREF(i,c:n) - A_RREF(i,c) * A_RREF(r,c:n)
      end if
    end do
!
!  Move to the next row.
!
    r = r + 1

  end do

  return
end
subroutine rref_determinant ( n, A, a_det )

!*****************************************************************************80
!
!! rref_determinant() uses reduced row echelon form (RREF) to compute a determinant.
!
!  Discussion:
!
!    The procedure will fail if A is not square.
!
!    This is simply a demonstration of how the RREF can be used to compute
!    the determinant.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 December 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer n: the number of rows and columns in the matrix A.
!
!    real A(N,N), a square matrix.
!
!  Output:
!
!    real A_DET: the determinant.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  real ( kind = rk8 ) A(n,n)
  integer a_cols(n)
  real ( kind = rk8 ) a_det
  real ( kind = rk8 ) A_RREF(n,n)
  integer i
!
!  Do an RREF on A.
!
  call rref_compute ( n, n, A, A_RREF, a_cols )
!
!  Compute product of diagonal entries.
!
  a_det = 1.0
  do i = 1, n
    a_det = a_det * A_RREF(i,i)
  end do

  return
end
subroutine rref_inverse ( n, A, A_INV )

!*****************************************************************************80
!
!! rref_inverse() uses reduced row echelon form (RREF) to compute an inverse.
!
!  Discussion:
!
!    The procedure will fail if A is not square, or not invertible.
!
!    This is simply a demonstration of how RREF can be used to compute
!    the inverse.  But:
!    * there are better ways to compute the inverse, including
!      B = inv ( A )
!    * the inverse matrix is usually not the appropriate tool for solving
!      linear systems.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 December 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer n: the number of rows and columns in the matrix A.
!
!    real A(N,N), a square, invertible matrix.
!
!  Output:
!
!    real A_INV(N,N), the inverse of A, computed using the rref_compute().
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  real ( kind = rk8 ) A(n,n)
  integer a_cols(n)
  real ( kind = rk8 ) A_INV(n,n)
  real ( kind = rk8 ) A_RREF(n,n)
  real ( kind = rk8 ) AI(n,2*n)
  real ( kind = rk8 ) AI_RREF(n,2*n)
  integer ai_cols(2*n)
  integer i
  integer j
!
!  First do an RREF on A alone.
!  The second argument is the independent columns found in the input matrix.
!
  call rref_compute ( n, n, A, A_RREF, a_cols )

  do i = 1, n
    if ( a_cols(i) == 0 ) then
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'rref_inverse(): Warning!'
      write ( *, '(a)' ) '  The input matrix seems to be singular.'
      write ( *, '(a)' ) '  The inverse could not be computed.'
      stop ( 0 )
    end if
  end do
!
!  If A has N independent columns, append the identity and repeat the RREF.
!
  AI(1:n,1:n) = A(1:n,1:n)
  do i = 1, n
    do j = n + 1, 2 * n
      AI(i,j) = 0.0
    end do
    AI(i,i+n) = 1.0
  end do

  call rref_compute ( n, n + n, AI, AI_RREF, ai_cols )
!
!  The last N columns of AI_RREF are the inverse.
!
  A_INV = AI_RREF(1:n,n+1:2*n)

  return
end
subroutine rref_rank ( m, n, A, a_rank )

!*****************************************************************************80
!
!! rref_rank() returns the rank of a matrix, using the reduced row echelon form.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 December 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer m, n: the number of rows and columns in the matrix A.
!
!    real A(M,N), the matrix to be analyzed.
!
!  Output:
!
!    integer A_RANK, the estimated rank of A.
!    0 <= A_RANK <= min ( M, N ).
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk8 ) A(m,n)
  integer a_cols(n)
  integer a_rank
  real ( kind = rk8 ) A_RREF(m,n)
  integer i

  call rref_compute ( m, n, A, A_RREF, a_cols )

  a_rank = 0
  do i = 1, n
    if ( a_cols(i) /= 0 ) then
      a_rank = a_rank + 1
    end if
  end do

  return
end
subroutine rref_solve ( m, n, A, nb, b, x )

!*****************************************************************************80
!
!! rref_solve() uses reduced row echelon form (RREF) to solve a linear system.
!
!  Discussion:
!
!    A linear system A*x=b is given, where
!    A is an MxN1 matrix, possibly singular.
!    b is an Mx1 right hand side, or MxN2 collection of right hand sides.
!    x is the desired N1x1 solution, or N1xN2 collection of solutions.
!
!    The right hand sides are assumed to be consistent, that is,
!    each right hand side is assumed to be a linear combination of 
!    columns of A.  If this is not so, the procedure will still produce
!    a result x, but it will not satisfy the equation A*x=b. 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 December 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real A(M,N), a matrix.
!
!    real B(M,NB), one or more right hand sides consistent with A.
!
!  Output:
!
!    real X(N,NB), one or more solution vectors.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer m
  integer n
  integer nb

  real ( kind = rk8 ) A(m,n)
  real ( kind = rk8 ) AI(m,n+nb)
  real ( kind = rk8 ) AI_RREF(m,n+nb)
  integer ai_cols(n+nb)
  real ( kind = rk8 ) b(m,nb)
  real ( kind = rk8 ) x(n,nb)

  AI(1:m,1:n) = A(1:m,1:n)
  AI(1:m,n+1:n+nb) = b(1:m,1:nb)

  call rref_compute ( m, n + nb, AI, AI_RREF, ai_cols )

  x(1:n,1:nb) = AI_RREF(1:n,n+1:n+nb)

  return
end

