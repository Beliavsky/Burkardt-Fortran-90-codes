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

  do i = 1, n-1

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
    do i = 1, n-1
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
    do i = n-1, 1, -1
      b(i) = b(i) - a_lu(3,i) * b(i+1)
    end do

  end if

  return
end
subroutine r83p_det ( n, a_lu, work4, det )

!*****************************************************************************80
!
!! R83P_DET computes the determinant of a matrix factored by R83P_FA.
!
!  Discussion:
!
!    The R83P storage format stores a periodic tridiagonal matrix as 
!    a 3 by N array, in which each row corresponds to a diagonal, and 
!    column locations are preserved.  The matrix value 
!    A(1,N) is stored as the array entry A(3,N), and the matrix value
!    A(N,1) is stored as the array entry A(1,1).
!
!  Example:
!
!    Here is how an R83P matrix of order 5 would be stored:
!
!      A51 A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54 A15
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 3.
!
!    Input, real ( kind = rk ) A_LU(3,N), LU factors from R83P_FA.
!
!    Input, real ( kind = rk ) WORK4, factorization information from R83P_FA.
!
!    Output, real ( kind = rk ) DET, the determinant of the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a_lu(3,n)
  real ( kind = rk ) det
  real ( kind = rk ) work4

  det = product ( a_lu(2,1:n-1) ) * work4

  return
end
subroutine r83p_fa ( n, a, info, work2, work3, work4 )

!*****************************************************************************80
!
!! R83P_FA factors an R83P matrix.
!
!  Discussion:
!
!    The R83P storage format stores a periodic tridiagonal matrix as 
!    a 3 by N array, in which each row corresponds to a diagonal, and 
!    column locations are preserved.  The matrix value 
!    A(1,N) is stored as the array entry A(3,N), and the matrix value
!    A(N,1) is stored as the array entry A(1,1).
!
!    Once the matrix has been factored by R83P_FA, R83P_SL may be called
!    to solve linear systems involving the matrix.
!
!    The logical matrix has a form which is suggested by this diagram:
!
!      D1 U1          L1
!      L2 D2 U2
!         L3 D3 U3
!            L4 D4 U4
!               L5 D5 U5
!      U6          L6 D6
!
!    The algorithm treats the matrix as a border banded matrix:
!
!      ( A1  A2 )
!      ( A3  A4 )
!
!    where:
!
!      D1 U1          | L1
!      L2 D2 U2       |  0
!         L3 D3  U3    |  0
!            L4 D4 U4 |  0
!               L5 D5 | U5
!      ---------------+---
!      U6  0  0  0 L6 | D6
!
!  Example:
!
!    Here is how an R83P matrix of order 5 would be stored:
!
!      A51 A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54 A15
!
!  Method:
!
!    The algorithm rewrites the system as:
!
!         X1 + inverse(A1) A2 X2 = inverse(A1) B1
!
!      A3 X1 +             A4 X2 = B2
!
!    The first equation can be "solved" for X1 in terms of X2:
!
!         X1 = - inverse(A1) A2 X2 + inverse(A1) B1
!
!    allowing us to rewrite the second equation for X2 explicitly:
!
!      ( A4 - A3 inverse(A1) A2 ) X2 = B2 - A3 inverse(A1) B1
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 3.
!
!    Input/output, real ( kind = rk ) A(3,N).
!    On input, the periodic tridiagonal matrix.  
!    On output, the arrays have been modified to hold information
!    defining the border-banded factorization of submatrices A1
!    and A3.
!
!    Output, integer INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
!    Output, real ( kind = rk ) WORK2(N-1), WORK3(N-1), WORK4, 
!    factorization information.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(3,n)
  integer info
  integer job
  real ( kind = rk ) work2(n-1)
  real ( kind = rk ) work3(n-1)
  real ( kind = rk ) work4
!
!  Compute inverse(A1):
!
  call r83_np_fa ( n - 1, a, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R83P_FA - Fatal error!'
    write ( *, '(a,i8)' ) '  R83_NP_FA returned INFO = ', info
    write ( *, '(a)' ) '  Factoring failed for column INFO.'
    write ( *, '(a)' ) '  The tridiagonal matrix A1 is singular.'
    write ( *, '(a)' ) '  This algorithm cannot continue!'
    stop 1
  end if
!
!  WORK2 := inverse(A1) * A2.
!
  work2(1) = a(3,n)
  work2(2:n-2) = 0.0D+00
  work2(n-1) = a(1,n)

  job = 0
  call r83_np_sl ( n - 1, a, work2, job )
!
!  WORK3 := inverse ( A1' ) * A3'.
!
  work3(1) = a(1,1)
  work3(2:n-2) = 0.0D+00
  work3(n-1) = a(3,n-1)

  job = 1
  call r83_np_sl ( n - 1, a, work3, job )
!
!  A4 := ( A4 - A3 * inverse(A1) * A2 )
!
  work4 = a(2,n) - a(1,1) * work2(1) - a(3,n-1) * work2(n-1)

  if ( work4 == 0.0D+00 ) then
    info = n
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R83P_FA - Fatal error!'
    write ( *, '(a)' ) '  The factored A4 submatrix is zero.'
    write ( *, '(a)' ) '  This algorithm cannot continue!'
    stop 1
  end if

  return
end
subroutine r83p_indicator ( n, a )

!*****************************************************************************80
!
!! R83P_INDICATOR sets up an R83P indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R83P storage format stores a periodic tridiagonal matrix as 
!    a 3 by N array, in which each row corresponds to a diagonal, and 
!    column locations are preserved.  The matrix value 
!    A(1,N) is stored as the array entry A(3,N), and the matrix value
!    A(N,1) is stored as the array entry A(1,1).
!
!  Example:
!
!    Here is how an R83P matrix of order 5 would be stored:
!
!      A51 A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54 A15
!
!    Here are the values as stored in an indicator matrix:
!
!      51 12 23 34 45
!      11 22 33 44 55
!      21 32 43 54 15
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 January 2004
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
!    Output, real ( kind = rk ) A(3,N), the R83P indicator matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(3,n)
  integer fac
  integer i
  integer i4_log_10
  integer j

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  i = n
  j = 1
  a(1,j) = real ( fac * i + j, kind = rk )
  do j = 2, n
    i = j - 1
    a(1,j) = real ( fac * i + j, kind = rk )
  end do

  do j = 1, n
    i = j
    a(2,j) = real ( fac * i + j, kind = rk )
  end do

  do j = 1, n-1
    i = j + 1
    a(3,j) = real ( fac * i + j, kind = rk )
  end do
  i = 1
  j = n
  a(3,j) = real ( fac * i + j, kind = rk )

  return
end
subroutine r83p_ml ( n, a_lu, x, b, job )

!*****************************************************************************80
!
!! R83P_ML computes A * x or x * A, where A has been factored by R83P_FA.
!
!  Discussion:
!
!    The R83P storage format stores a periodic tridiagonal matrix as 
!    a 3 by N array, in which each row corresponds to a diagonal, and 
!    column locations are preserved.  The matrix value 
!    A(1,N) is stored as the array entry A(3,N), and the matrix value
!    A(N,1) is stored as the array entry A(1,1).
!
!  Example:
!
!    Here is how an R83P matrix of order 5 would be stored:
!
!      A51 A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54 A15
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
!    N must be at least 3.
!
!    Input, real ( kind = rk ) A_LU(3,N), the LU factors ffrom R83P_FA.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied by the matrix.
!
!    Output, real ( kind = rk ) B(N), the result of the multiplication.
!
!    Input, integer JOB, indicates what product should be computed.
!    0, compute A * x.
!    nonzero, compute A' * x.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a_lu(3,n)
  real ( kind = rk ) b(n)
  integer job
  real ( kind = rk ) x(n)
!
!  Multiply A(1:N-1,1:N-1) and X(1:N-1).
!
  call r83_np_ml ( n - 1, a_lu, x, b, job )
!
!  Add terms from the border.
!
  if ( job == 0 ) then
    b(1) = b(1) + a_lu(3,n) * x(n)
    b(n-1) = b(n-1) + a_lu(1,n) * x(n)
    b(n) = a_lu(1,1) * x(1) + a_lu(3,n-1) * x(n-1) + a_lu(2,n) * x(n)
  else
    b(1) = b(1) + a_lu(1,1) * x(n)
    b(n-1) = b(n-1) + a_lu(3,n-1) * x(n)
    b(n) = a_lu(3,n) * x(1) + a_lu(1,n) * x(n-1) + a_lu(2,n) * x(n)
  end if

  return
end
subroutine r83p_mtv ( n, a, x, b )

!*****************************************************************************80
!
!! R83P_MTV multiplies an R8VEC by an R83P matrix.
!
!  Discussion:
!
!    The R83P storage format stores a periodic tridiagonal matrix as 
!    a 3 by N array, in which each row corresponds to a diagonal, and 
!    column locations are preserved.  The matrix value 
!    A(1,N) is stored as the array entry A(3,N), and the matrix value
!    A(N,1) is stored as the array entry A(1,1).
!
!  Example:
!
!    Here is how an R83P matrix of order 5 would be stored:
!
!      A51 A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54 A15
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 3.
!
!    Input, real ( kind = rk ) A(3,N), the R83P matrix.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = rk ) B(N), the product X * A.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(3,n)
  real ( kind = rk ) b(n)
  integer i
  real ( kind = rk ) x(n)

  b(1) =   a(1,1)   * x(n)   + a(2,1) * x(1) + a(3,1)   * x(2)

  do i = 2, n-1
    b(i) = a(1,i)   * x(i-1) + a(2,i) * x(i) + a(3,i)   * x(i+1)
  end do

  b(n) =   a(1,n)   * x(n-1) + a(2,n) * x(n) + a(3,n)   * x(1)

  return
end
subroutine r83p_mv ( n, a, x, b )

!*****************************************************************************80
!
!! R83P_MV multiplies an R83P matrix by an R8VEC.
!
!  Discussion:
!
!    The R83P storage format stores a periodic tridiagonal matrix as 
!    a 3 by N array, in which each row corresponds to a diagonal, and 
!    column locations are preserved.  The matrix value 
!    A(1,N) is stored as the array entry A(3,N), and the matrix value
!    A(N,1) is stored as the array entry A(1,1).
!
!  Example:
!
!    Here is how an R83P matrix of order 5 would be stored:
!
!      A51 A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54 A15
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 3.
!
!    Input, real ( kind = rk ) A(3,N), the R83P matrix.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = rk ) B(N), the product A * x.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(3,n)
  real ( kind = rk ) b(n)
  integer i
  real ( kind = rk ) x(n)

  b(1) =   a(3,n)   * x(n)   + a(2,1) * x(1) + a(1,2)   * x(2)

  do i = 2, n-1
    b(i) = a(3,i-1) * x(i-1) + a(2,i) * x(i) + a(1,i+1) * x(i+1)
  end do

  b(n) =   a(3,n-1) * x(n-1) + a(2,n) * x(n) + a(1,1)   * x(1)

  return
end
subroutine r83p_print ( n, a, title )

!*****************************************************************************80
!
!! R83P_PRINT prints an R83P matrix.
!
!  Discussion:
!
!    The R83P storage format stores a periodic tridiagonal matrix as 
!    a 3 by N array, in which each row corresponds to a diagonal, and 
!    column locations are preserved.  The matrix value 
!    A(1,N) is stored as the array entry A(3,N), and the matrix value
!    A(N,1) is stored as the array entry A(1,1).
!
!  Example:
!
!    Here is how an R83P matrix of order 5 would be stored:
!
!      A51 A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54 A15
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
!    Input, real ( kind = rk ) A(3,N), the R83P matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(3,n)
  character ( len = * ) title

  call r83p_print_some ( n, a, 1, 1, n, n, title )

  return
end
subroutine r83p_print_some ( n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R83P_PRINT_SOME prints some of an R83P matrix.
!
!  Discussion:
!
!    The R83P storage format stores a periodic tridiagonal matrix as 
!    a 3 by N array, in which each row corresponds to a diagonal, and 
!    column locations are preserved.  The matrix value 
!    A(1,N) is stored as the array entry A(3,N), and the matrix value
!    A(N,1) is stored as the array entry A(1,1).
!
!  Example:
!
!    Here is how an R83P matrix of order 5 would be stored:
!
!      A51 A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54 A15
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 January 2004
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
!    Input, real ( kind = rk ) A(3,N), the R83P matrix.
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
!  Determine the column range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )

    if ( 1 < i2lo .or. j2hi < n ) then
      i2lo = max ( i2lo, j2lo - 1 )
    end if

    i2hi = min ( ihi, n )

    if ( i2hi < n .or. 1 < j2lo ) then
      i2hi = min ( i2hi, j2hi + 1 )
    end if

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( i == n .and. j == 1 ) then
          write ( ctemp(j2), '(g14.6)' ) a(1,j)
        else if ( i == 1 .and. j == n ) then
          write ( ctemp(j2), '(g14.6)' ) a(3,j)
        else if ( 1 < i - j .or. 1 < j - i ) then
          ctemp(j2) = '              '
        else if ( j == i + 1 ) then
          write ( ctemp(j2), '(g14.6)' ) a(1,j)
        else if ( j == i ) then
          write ( ctemp(j2), '(g14.6)' ) a(2,j)
        else if ( j == i - 1 ) then
          write ( ctemp(j2), '(g14.6)' ) a(3,j)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r83p_random ( n, seed, a )

!*****************************************************************************80
!
!! R83P_RANDOM randomizes an R83P matrix.
!
!  Discussion:
!
!    The R83P storage format stores a periodic tridiagonal matrix as 
!    a 3 by N array, in which each row corresponds to a diagonal, and 
!    column locations are preserved.  The matrix value 
!    A(1,N) is stored as the array entry A(3,N), and the matrix value
!    A(N,1) is stored as the array entry A(1,1).
!
!  Example:
!
!    Here is how an R83P matrix of order 5 would be stored:
!
!      A51 A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54 A15
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 May 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 3.
!
!    Input/output, integer SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = rk ) A(3,N), the R83P matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(3,n)
  integer seed

  call r8mat_uniform_01 ( 3, n, seed, a )

  return
end
subroutine r83p_sl ( n, a_lu, b, x, job, work2, work3, work4 )

!*****************************************************************************80
!
!! R83P_SL solves an R83P system.
!
!  Discussion:
!
!    The R83P storage format stores a periodic tridiagonal matrix as 
!    a 3 by N array, in which each row corresponds to a diagonal, and 
!    column locations are preserved.  The matrix value 
!    A(1,N) is stored as the array entry A(3,N), and the matrix value
!    A(N,1) is stored as the array entry A(1,1).
!
!    The linear system must have been factored by R83P_FA.
!
!  Example:
!
!    Here is how an R83P matrix of order 5 would be stored:
!
!      A51 A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54 A15
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 3.
!
!    Input, real ( kind = rk ) A_LU(3,N), the LU factors from R83P_FA.
!
!    Input, real ( kind = rk ) B(N), the right hand side of the linear system.
!
!    Output, real ( kind = rk ) X(N), the solution to the linear system.
!
!    Input, integer JOB, specifies the system to solve.
!    0, solve A * x = b.
!    nonzero, solve A' * x = b.
!
!    Input, real ( kind = rk ) WORK2(N-1), WORK3(N-1), WORK4, 
!    factor data from R83P_FA.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a_lu(3,n)
  real ( kind = rk ) b(n)
  integer job
  real ( kind = rk ) work2(n-1)
  real ( kind = rk ) work3(n-1)
  real ( kind = rk ) work4
  real ( kind = rk ) x(n)

  x(1:n) = b(1:n)

  if ( job == 0 ) then
!
!  Solve A1 * X1 = B1.
!
    call r83_np_sl ( n - 1, a_lu, x, job )
!
!  X2 = B2 - A3 * X1
!
    x(n) = x(n) - a_lu(1,1) * x(1) - a_lu(3,n-1) * x(n-1)
!
!  Solve A4 * X2 = X2
!
    x(n) = x(n) / work4
!
!  X1 := X1 - inverse ( A1 ) * A2 * X2.
!
    x(1:n-1) = x(1:n-1) - work2(1:n-1) * x(n)

  else
!
!  Solve A1' * X1 = B1.
!
    call r83_np_sl ( n - 1, a_lu, x, job )
!
!  X2 := X2 - A2' * B1
!
    x(n) = x(n) - a_lu(3,n) * x(1) - a_lu(1,n) * x(n-1)
!
!  Solve A4 * X2 = X2.
!
    x(n) = x(n) / work4
!
!  X1 := X1 - transpose ( inverse ( A1 ) * A3 ) * X2.
!
    x(1:n-1) = x(1:n-1) - work3(1:n-1) * x(n)

  end if

  return
end
subroutine r83p_to_r8ge ( n, a, b )

!*****************************************************************************80
!
!! R83P_TO_R8GE copies an R83P matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R83P storage format stores a periodic tridiagonal matrix as 
!    a 3 by N array, in which each row corresponds to a diagonal, and 
!    column locations are preserved.  The matrix value 
!    A(1,N) is stored as the array entry A(3,N), and the matrix value
!    A(N,1) is stored as the array entry A(1,1).
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!  Example:
!
!    Here is how an R83P matrix of order 5 would be stored:
!
!      A51 A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54 A15
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 3.
!
!    Input, real ( kind = rk ) A(3,N), the R83P matrix.
!
!    Output, real ( kind = rk ) B(N,N), the R8GE matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(3,n)
  real ( kind = rk ) b(n,n)
  integer i
  integer j

  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        b(i,j) = a(2,j)
      else if ( j == i-1 ) then
        b(i,j) = a(3,j)
      else if ( j == i+1 ) then
        b(i,j) = a(1,j)
      else if ( i == 1 .and. j == n ) then
        b(i,j) = a(3,j)
      else if ( i == n .and. j == 1 ) then
        b(i,j) = a(1,j)
      else
        b(i,j) = 0.0D+00
      end if
    end do
  end do

  return
end
subroutine r83p_zeros ( n, a )

!*****************************************************************************80
!
!! R83P_ZEROS zeroes an R83P matrix.
!
!  Discussion:
!
!    The R83P storage format stores a periodic tridiagonal matrix as 
!    a 3 by N array, in which each row corresponds to a diagonal, and 
!    column locations are preserved.  The matrix value 
!    A(1,N) is stored as the array entry A(3,N), and the matrix value
!    A(N,1) is stored as the array entry A(1,1).
!
!  Example:
!
!    Here is how an R83P matrix of order 5 would be stored:
!
!      A51 A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54 A15
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
!    N must be at least 3.
!
!    Output, real ( kind = rk ) A(3,N), the R83P matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(3,n)

  a(1:3,1:n) = 0.0D+00

  return
end
subroutine r8ge_det ( n, a_lu, pivot, det )

!*****************************************************************************80
!
!! R8GE_DET: determinant of a matrix factored by R8GE_FA or R8GE_TRF.
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
!    29 March 2003
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
!    Input, real ( kind = rk ) A_LU(N,N), the LU factors from R8GE_FA 
!    or R8GE_TRF.
!
!    Input, integer PIVOT(N), as computed by R8GE_FA or R8GE_TRF.
!
!    Output, real ( kind = rk ) DET, the determinant of the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a_lu(n,n)
  real ( kind = rk ) det
  integer i
  integer pivot(n)

  det = 1.0D+00

  do i = 1, n
    det = det * a_lu(i,i)
    if ( pivot(i) /= i ) then
      det = - det
    end if
  end do

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
subroutine r8mat_uniform_01 ( m, n, seed, r )

!*****************************************************************************80
!
!! R8MAT_UNIFORM_01 fills an R8MAT with unit pseudorandom numbers.
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
!    Output, real ( kind = rk ) R(M,N), the array of pseudorandom values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer i
  integer, parameter :: i4_huge = 2147483647
  integer j
  integer k
  integer seed
  real ( kind = rk ) r(m,n)

  do j = 1, n

    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge
      end if

      r(i,j) = real ( seed, kind = rk ) * 4.656612875D-10

    end do
  end do

  return
end
subroutine r8vec_indicator1 ( n, a )

!*****************************************************************************80
!
!! R8VEC_INDICATOR1 sets an R8VEC to the indicator vector (1,2,3,...).
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
    a(i) = real ( i, kind = rk )
  end do

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
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
!    22 August 2000
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
    write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
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

