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
!  Input:
!
!    integer I, the number whose logarithm base 10 
!    is desired.
!
!  Output:
!
!    integer I4_LOG_10, the integer part of the 
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
function r8_uniform_01 ( )

!*****************************************************************************80
!
!! r8_uniform_01() returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = rk ) value.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 September 2021
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
!  Output:
!
!    real ( kind = rk ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r
  real ( kind = rk ) r8_uniform_01

  call random_number ( harvest = r )

  r8_uniform_01 = r

  return
end
subroutine r8blt_check ( n, ml, ierror )

!*****************************************************************************80
!
!! r8blt_check() checks the dimensions of a banded lower triangular matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the order of the matrix.
!    N must be positive.
!
!    integer ML, the lower bandwidth of the matrix.
!    ML must be at least 0, and no greater than N-1.
!
!  Output:
!
!    integer IERROR, reports whether any errors were detected.
!    IERROR is set to 0 before the checks are made, and then:
!    IERROR = IERROR + 1 if ML is illegal;
!    IERROR = IERROR + 2 if N is illegal.
!
  implicit none

  integer ierror
  integer ml
  integer n

  ierror = 0

  if ( ml < 0 .or. n - 1 < ml ) then
    ierror = ierror + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) 'r8blt_check(): Illegal ML = ', ml
  end if

  if ( n <= 0 ) then
    ierror = ierror + 2
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) 'r8blt_check(): Illegal N = ', n
    return
  end if

  return
end
subroutine r8blt_det ( n, ml, a, det )

!*****************************************************************************80
!
!! r8blt_det() computes the determinant of an R8BLT matrix.
!
!  Discussion:
!
!    The R8BLT storage format is used for a banded lower triangular matrix.
!    The matrix is assumed to be zero below the ML-th subdiagonal.
!    The matrix is stored in an ML+1 by N array, in which the diagonal
!    appears in the first row, followed by successive subdiagonals.
!    Columns are preserved.
!
!  Example:
!
!    N = 5, ML = 2
!
!    A11   0   0   0   0
!    A21 A22   0   0   0
!    A31 A32 A33   0   0
!      0 A42 A43 A44   0
!      0   0 A53 A54 A55
!                --- ---
!                    ---
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the order of the matrix.
!
!    integer ML, the lower bandwidth.
!
!    real ( kind = rk ) A(ML+1,N), the R8BLT matrix.
!
!  Output:
!
!    real ( kind = rk ) DET, the determinant of A.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer n

  real ( kind = rk ) a(ml+1,n)
  real ( kind = rk ) det

  det = product ( a(1,1:n) )

  return
end
subroutine r8blt_indicator ( n, ml, a )

!*****************************************************************************80
!
!! r8blt_indicator() sets up an R8BLT indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R8BLT storage format is used for a banded lower triangular matrix.
!    The matrix is assumed to be zero below the ML-th subdiagonal.
!    The matrix is stored in an ML+1 by N array, in which the diagonal
!    appears in the first row, followed by successive subdiagonals.
!    Columns are preserved.
!
!  Example:
!
!    N = 5, ML = 2
!
!    A11   0   0   0   0
!    A21 A22   0   0   0
!    A31 A32 A33   0   0
!      0 A42 A43 A44   0
!      0   0 A53 A54 A55
!                --- ---
!                    ---
!
!    The indicator matrix is stored as:
!
!      11  22  33  44  55
!      21  32  43  54   0
!      31  42  53   0   0
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of columns of the matrix.
!
!    integer ML, the lower bandwidth.
!
!  Output:
!
!    real ( kind = rk ) A(ML+1,N), the R8BLT matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer n

  real ( kind = rk ) a(ml+1,n)
  integer fac
  integer i
  integer i4_log_10
  integer j

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  do i = 1, n
    do j = max ( 1, i - ml ), i
      a(i-j+1,j) = real ( fac * i + j, kind = rk )
    end do
  end do

  do i = n+1, n+ml
    do j = i-ml, n
      a(i-j+1,j) = 0.0D+00
    end do
  end do

  return
end
subroutine r8blt_mtv ( n, ml, a, x, b )

!*****************************************************************************80
!
!! r8blt_mtv() multiplies an R8VEC by an R8BLT matrix.
!
!  Discussion:
!
!    The R8BLT storage format is used for a banded lower triangular matrix.
!    The matrix is assumed to be zero below the ML-th subdiagonal.
!    The matrix is stored in an ML+1 by N array, in which the diagonal
!    appears in the first row, followed by successive subdiagonals.
!    Columns are preserved.
!
!  Example:
!
!    N = 5, ML = 2
!
!    A11   0   0   0   0
!    A21 A22   0   0   0
!    A31 A32 A33   0   0
!      0 A42 A43 A44   0
!      0   0 A53 A54 A55
!                --- ---
!                    ---
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the order of the matrix.
!
!    integer ML, the lower bandwidth.
!
!    real ( kind = rk ) A(ML+1,N), the R8BLT matrix.
!
!    real ( kind = rk ) X(N), the vector to be multiplied by A.
!
!  Output:
!
!    real ( kind = rk ) B(N), the product X*A.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer n

  real ( kind = rk ) a(ml+1,n)
  real ( kind = rk ) b(n)
  integer i
  integer j
  integer jhi
  integer jlo
  real ( kind = rk ) x(n)

  b(1:n) = 0.0D+00

  do i = 1, n
    jlo = max ( 1, i - ml )
    jhi = i
    do j = jlo, jhi
      b(j) = b(j) + x(i) * a(i-j+1,j)
    end do
  end do

  return
end
subroutine r8blt_mv ( n, ml, a, x, b )

!*****************************************************************************80
!
!! r8blt_mv() multiplies an R8BLT matrix by an R8VEC.
!
!  Discussion:
!
!    The R8BLT storage format is used for a banded lower triangular matrix.
!    The matrix is assumed to be zero below the ML-th subdiagonal.
!    The matrix is stored in an ML+1 by N array, in which the diagonal
!    appears in the first row, followed by successive subdiagonals.
!    Columns are preserved.
!
!  Example:
!
!    N = 5, ML = 2
!
!    A11   0   0   0   0
!    A21 A22   0   0   0
!    A31 A32 A33   0   0
!      0 A42 A43 A44   0
!      0   0 A53 A54 A55
!                --- ---
!                    ---
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the order of the matrix.
!
!    integer ML, the lower bandwidth.
!
!    real ( kind = rk ) A(ML+1,N), the R8BLT matrix.
!
!    real ( kind = rk ) X(N), the vector to be multiplied by A.
!
!  Output:
!
!    real ( kind = rk ) B(N), the product A * x.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer n

  real ( kind = rk ) a(ml+1,n)
  real ( kind = rk ) b(n)
  integer i
  integer j
  integer jhi
  integer jlo
  real ( kind = rk ) x(n)

  do i = 1, n
    b(i) = 0.0D+00
    jlo = max ( 1, i - ml )
    jhi = i
    do j = jlo, jhi
      b(i) = b(i) + a(i-j+1,j) * x(j)
    end do
  end do

  return
end
subroutine r8blt_print ( n, ml, a, title )

!*****************************************************************************80
!
!! r8blt_print() prints an R8BLT matrix.
!
!  Discussion:
!
!    The R8BLT storage format is used for a banded lower triangular matrix.
!    The matrix is assumed to be zero below the ML-th subdiagonal.
!    The matrix is stored in an ML+1 by N array, in which the diagonal
!    appears in the first row, followed by successive subdiagonals.
!    Columns are preserved.
!
!  Example:
!
!    N = 5, ML = 2
!
!    A11   0   0   0   0
!    A21 A22   0   0   0
!    A31 A32 A33   0   0
!      0 A42 A43 A44   0
!      0   0 A53 A54 A55
!                --- ---
!                    ---
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the order of the matrix.
!
!    integer ML, the lower bandwidth.
!
!    real ( kind = rk ) A(ML+1,N), the R8BLT matrix.
!
!    character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer n

  real ( kind = rk ) a(ml+1,n)
  character ( len = * ) title

  call r8blt_print_some ( n, ml, a, 1, 1, n, n, title )

  return
end
subroutine r8blt_print_some ( n, ml, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! r8blt_print_some() prints some of an R8BLT matrix.
!
!  Discussion:
!
!    The R8BLT storage format is used for a banded lower triangular matrix.
!    The matrix is assumed to be zero below the ML-th subdiagonal.
!    The matrix is stored in an ML+1 by N array, in which the diagonal
!    appears in the first row, followed by successive subdiagonals.
!    Columns are preserved.
!
!  Example:
!
!    N = 5, ML = 2
!
!    A11   0   0   0   0
!    A21 A22   0   0   0
!    A31 A32 A33   0   0
!      0 A42 A43 A44   0
!      0   0 A53 A54 A55
!                --- ---
!                    ---
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the order of the matrix.
!
!    integer ML, the lower bandwidth.
!
!    real ( kind = rk ) A(ML+1,N), the R8BLT matrix.
!
!    integer ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: incx = 5
  integer ml
  integer n

  real ( kind = rk ) a(ml+1,n)
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
    i2lo = max ( i2lo, j2lo )
    i2hi = min ( ihi, n )
    i2hi = min ( i2hi, j2hi + ml )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( j <= i .and. i <= j + ml ) then
          aij = a(i-j+1,j)
        else
          aij = 0.0D+00
        end if

        if ( ml < i - j .or. 0 < j - i ) then
          ctemp(j2) = '              '
        else
          write ( ctemp(j2), '(g14.6)' ) aij
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8blt_random ( n, ml, a )

!*****************************************************************************80
!
!! r8blt_random() randomizes an R8BLT matrix.
!
!  Discussion:
!
!    The R8BLT storage format is used for a banded lower triangular matrix.
!    The matrix is assumed to be zero below the ML-th subdiagonal.
!    The matrix is stored in an ML+1 by N array, in which the diagonal
!    appears in the first row, followed by successive subdiagonals.
!    Columns are preserved.
!
!  Example:
!
!    N = 5, ML = 2
!
!    A11   0   0   0   0
!    A21 A22   0   0   0
!    A31 A32 A33   0   0
!      0 A42 A43 A44   0
!      0   0 A53 A54 A55
!                --- ---
!                    ---
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
!  Input:
!
!    integer N, the number of columns of the matrix.
!
!    integer ML, the lower bandwidth.
!
!  Output:
!
!    real ( kind = rk ) A(ML+1,N), the R8BLT matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer n

  real ( kind = rk ) a(ml+1,n)
  real ( kind = rk ) r8_uniform_01
  integer i
  integer j

  do i = 1, n
    do j = max ( 1, i - ml ), i
      a(i-j+1,j) = r8_uniform_01 ( )
    end do
  end do
!
!  The junk entries can be thought of as corresponding to
!  elements of a phantom portion of the matrix.
!
  do i = n + 1, n + ml
    do j = i-ml, n
      a(i-j+1,j) = 0.0D+00
    end do
  end do

  return
end
subroutine r8blt_sl ( n, ml, a, b )

!*****************************************************************************80
!
!! r8blt_sl() solves an R8BLT system A*x=b.
!
!  Discussion:
!
!    The R8BLT storage format is used for a banded lower triangular matrix.
!    The matrix is assumed to be zero below the ML-th subdiagonal.
!    The matrix is stored in an ML+1 by N array, in which the diagonal
!    appears in the first row, followed by successive subdiagonals.
!    Columns are preserved.
!
!    No factorization of the lower triangular matrix is required.
!
!  Example:
!
!    N = 5, ML = 2
!
!    A11   0   0   0   0
!    A21 A22   0   0   0
!    A31 A32 A33   0   0
!      0 A42 A43 A44   0
!      0   0 A53 A54 A55
!                --- ---
!                    ---
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 October 2015
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the order of the matrix.
!
!    integer ML, the lower bandwidth.
!
!    real ( kind = rk ) A(ML+1,N), the R8BLT matrix.
!
!    real ( kind = rk ) B(N): the right hand side.
!
!  Output:
!
!    real ( kind = rk ) B(N): the solution vector.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer n

  real ( kind = rk ) a(ml+1,n)
  real ( kind = rk ) b(n)
  integer i
  integer ihi
  integer j

  do j = 1, n
    b(j) = b(j) / a(1,j)
    ihi = min ( j + ml, n )
    do i = j + 1, ihi
      b(i) = b(i) - a(i-j+1,j) * b(j)
    end do
  end do

  return
end
subroutine r8blt_slt ( n, ml, a, b )

!*****************************************************************************80
!
!! R8BLT_SLT solves an R8BLT system A'*x=b.
!
!  Discussion:
!
!    The R8BLT storage format is used for a banded lower triangular matrix.
!    The matrix is assumed to be zero below the ML-th subdiagonal.
!    The matrix is stored in an ML+1 by N array, in which the diagonal
!    appears in the first row, followed by successive subdiagonals.
!    Columns are preserved.
!
!    No factorization of the lower triangular matrix is required.
!
!  Example:
!
!    N = 5, ML = 2
!
!    A11   0   0   0   0
!    A21 A22   0   0   0
!    A31 A32 A33   0   0
!      0 A42 A43 A44   0
!      0   0 A53 A54 A55
!                --- ---
!                    ---
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the order of the matrix.
!
!    integer ML, the lower bandwidth.
!
!    real ( kind = rk ) A(ML+1,N), the R8BLT matrix.
!
!    real ( kind = rk ) B(N): the right hand side.
!
!  Output:
!
!    real ( kind = rk ) B(N): the solution.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer n

  real ( kind = rk ) a(ml+1,n)
  real ( kind = rk ) b(n)
  integer i
  integer ilo
  integer j

  do j = n, 1, -1
    b(j) = b(j) / a(1,j)
    ilo = max ( j - ml, 1 )
    do i = ilo, j-1
      b(i) = b(i) - a(j-i+1,i) * b(j)
    end do
  end do

  return
end
subroutine r8blt_to_r8ge ( n, ml, a, b )

!*****************************************************************************80
!
!! R8BLT_TO_R8GE copies an R8BLT matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R8BLT storage format is used for a banded lower triangular matrix.
!    The matrix is assumed to be zero below the ML-th subdiagonal.
!    The matrix is stored in an ML+1 by N array, in which the diagonal
!    appears in the first row, followed by successive subdiagonals.
!    Columns are preserved.
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
!    N = 5, ML = 2
!
!    A11   0   0   0   0
!    A21 A22   0   0   0
!    A31 A32 A33   0   0
!      0 A42 A43 A44   0
!      0   0 A53 A54 A55
!                --- ---
!                    ---
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    28 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the order of the matrices.
!    N must be positive.
!
!    integer ML, the lower bandwidth of A.
!    ML must be nonnegative, and no greater than N-1.
!
!    real ( kind = rk ) A(ML+1,N), the R8BLT matrix.
!
!  Output:
!
!    real ( kind = rk ) B(N,N), the R8GE matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer n

  real ( kind = rk ) a(ml+1,n)
  real ( kind = rk ) b(n,n)
  integer i
  integer j

  do i = 1, n
    do j = 1, n
      if ( j <= i .and. i <= j + ml ) then
        b(i,j) = a(i-j+1,j)
      else
        b(i,j) = 0.0D+00
      end if
    end do
  end do

  return
end
subroutine r8blt_zeros ( n, ml, a )

!*****************************************************************************80
!
!! R8BLT_ZEROS zeroes an R8BLT matrix.
!
!  Discussion:
!
!    The R8BLT storage format is used for a banded lower triangular matrix.
!    The matrix is assumed to be zero below the ML-th subdiagonal.
!    The matrix is stored in an ML+1 by N array, in which the diagonal
!    appears in the first row, followed by successive subdiagonals.
!    Columns are preserved.
!
!  Example:
!
!    N = 5, ML = 2
!
!    A11   0   0   0   0
!    A21 A22   0   0   0
!    A31 A32 A33   0   0
!      0 A42 A43 A44   0
!      0   0 A53 A54 A55
!                --- ---
!                    ---
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
!  Input:
!
!    integer N, the number of columns of the matrix.
!
!    integer ML, the lower bandwidth.
!
!  Output:
!
!    real ( kind = rk ) A(ML+1,N), the R8BLT matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ml
  integer n

  real ( kind = rk ) a(ml+1,n)

  a(1:ml+1,1:n) = 0.0D+00

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
!  Input:
!
!    integer M, the number of rows of the matrix.
!    M must be positive.
!
!    integer N, the number of columns of the matrix.
!    N must be positive.
!
!    real ( kind = rk ) A(M,N), the R8GE matrix.
!
!    character ( len = * ) TITLE, a title.
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
!  Input:
!
!    integer M, the number of rows of the matrix.
!    M must be positive.
!
!    integer N, the number of columns of the matrix.
!    N must be positive.
!
!    real ( kind = rk ) A(M,N), the R8GE matrix.
!
!    integer ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    character ( len = * ) TITLE, a title.
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
!  Input:
!
!    integer N, the number of elements of A.
!
!  Output:
!
!    real ( kind = rk ) A(N), the indicator values.
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
!  Input:
!
!    integer N, the number of components of the vector.
!
!    real ( kind = rk ) A(N), the vector to be printed.
!
!    character ( len = * ) TITLE, a title.
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
subroutine r8vec_uniform_01 ( n, r )

!*****************************************************************************80
!
!! r8vec_uniform_01() returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of real ( kind = rk ) values.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
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
!  Input:
!
!    integer N, the number of entries in the vector.
!
!  Output:
!
!    real ( kind = rk ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) r(n)

  call random_number ( harvest = r(1:n) )

  return
end

