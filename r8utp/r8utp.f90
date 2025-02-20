function i4_log_10 ( i )

!*****************************************************************************80
!
!! i4_log_10() returns the integer part of the logarithm base 10 of an I4.
!
!  Discussion:
!
!    I4_LOG_10 ( I ) + 1 is the number of decimal digits in I.
!
!    An I4 is an integer value.
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
subroutine r8ge_print ( m, n, a, title )

!*****************************************************************************80
!
!! r8ge_print() prints an R8GE matrix.
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
!! r8ge_print_some() prints some of an R8GE matrix.
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
subroutine r8ge_to_r8utp ( m, n, a_ge, a_utp )

!*****************************************************************************80
!
!! r8ge_to_r8utp() copies an R8GE matrix to an R8UTP matrix.
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
!    The R8UTP storage format is appropriate for an upper triangular
!    matrix.  Only the upper triangle of the matrix is stored,
!    by successive partial columns, in an array which contains
!    (A11,A12,A22,A13,A23,A33,A14,...,AMN).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 August 2022
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, N, the number of rows and columns of 
!    the matrix.
!
!    real ( kind = rk ) A_GE(N,N), the R8GE matrix.
!
!  Output:
!
!    real ( kind = rk ) A_UTP(*), the R8UTP matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a_ge(m,n)
  real ( kind = rk ) a_utp(*)
  integer i
  integer j
  integer k

  k = 0
  do j = 1, n
    do i = 1, min ( j, m )
      k = k + 1
      a_utp(k) = a_ge(i,j)
    end do
  end do

  return
end
subroutine r8utp_det ( m, n, a, det )

!*****************************************************************************80
!
!! r8utp_det() computes the determinant of an R8UTP matrix.
!
!  Discussion:
!
!    The R8UTP storage format is appropriate for an upper triangular
!    matrix.  Only the upper triangle of the matrix is stored,
!    by successive partial columns, in an array which contains
!    (A11,A12,A22,A13,A23,A33,A14,...,AMN).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 August 2022
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, N, the order of the matrix.
!
!    real ( kind = rk ) A(*), the R8UTP matrix.
!
!  Output:
!
!    real ( kind = rk ) DET, the determinant of the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(*)
  real ( kind = rk ) det
  integer j
  integer k

  if ( m /= n ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'r8utp_det(): Fatal error!'
    write ( *, '(a)' ) '  m and n are not equal.'
    stop ( 1 )
  end if

  det = 1.0D+00
  k = 1
  do j = 1, n
    det = det * a(k)
    k = k + j + 1
  end do

  return
end
subroutine r8utp_indicator ( m, n, a )

!*****************************************************************************80
!
!! r8utp_indicator() sets up a R8UTP indicator matrix.
!
!  Discussion:
!
!    The R8UTP storage format is appropriate for an upper triangular
!    matrix.  Only the upper triangle of the matrix is stored,
!    by successive partial columns, in an array which contains
!    (A11,A12,A22,A13,A23,A33,A14,...,AMN).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 August 2022
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, N, the number of rows and columns of the matrix.
!
!  Output:
!
!    real ( kind = rk ) A(*), the R8UTP matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(*)
  integer fac
  integer i
  integer i4_log_10
  integer j
  integer k
  integer m
  integer n

  fac = 10**( i4_log_10 ( n ) + 1 )

  k = 0
  do j = 1, n
    do i = 1, min ( j, m )
      k = k + 1
      a(k) = fac * i + j
    end do
  end do

  return
end
subroutine r8utp_print ( m, n, a, title )

!*****************************************************************************80
!
!! r8utp_print() prints an R8UTP matrix.
!
!  Discussion:
!
!    The R8UTP storage format is appropriate for an upper triangular
!    matrix.  Only the upper triangle of the matrix is stored,
!    by successive partial columns, in an array which contains
!    (A11,A12,A22,A13,A23,A33,A14,...,AMN).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    13 August 2022
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, N, the order of the matrix.
!
!    real ( kind = rk ) A(*), the matrix.
!
!    character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(*)
  character ( len = * ) title

  call r8utp_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8utp_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! r8utp_print_some() prints some of an R8UTP matrix.
!
!  Discussion:
!
!    The R8UTP storage format is appropriate for an upper triangular
!    matrix.  Only the upper triangle of the matrix is stored,
!    by successive partial columns, in an array which contains
!    (A11,A12,A22,A13,A23,A33,A14,...,AMN).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 August 2022
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, N, the order of the matrix.
!
!    real ( kind = rk ) A(*), the matrix.
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

  real ( kind = rk ) a(*)
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
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( i <= j ) then
          aij = a(i+(j*(j-1))/2)
        else
          aij = 0.0D+00
        end if

        write ( ctemp(j2), '(g14.6)' ) aij

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8utp_random ( m, n, a )

!*****************************************************************************80
!
!! r8utp_random() randomizes an R8UTP matrix.
!
!  Discussion:
!
!    The R8UTP storage format is appropriate for an upper triangular
!    matrix.  Only the upper triangle of the matrix is stored,
!    by successive partial columns, in an array which contains
!    (A11,A12,A22,A13,A23,A33,A14,...,AMN).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 August 2022
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, N, the number of rows and columns of the matrix.
!    M and N must be positive.
!
!  Output:
!
!    real ( kind = rk ) A(*), the R8UTP matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(*)
  integer m
  integer mn
  integer n
  integer r8utp_size

  mn = r8utp_size ( m, n )

  call random_number ( a(1:mn ) )

  return
end
function r8utp_size ( m, n )

!*****************************************************************************80
!
!! r8utp_size() returns the size of an M x N R8UTP matrix.
!
!  Discussion:
!
!    The R8UTP storage format is appropriate for an upper triangular
!    matrix.  Only the upper triangle of the matrix is stored,
!    by successive partial columns, in an array which contains
!    (A11,A12,A22,A13,A23,A33,A14,...,AMN).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 August 2022
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, N: the number of rows and columns.
!
!  Output:
!
!    integer r8utp_size: the length of the array needed to store the matrix.
!
  implicit none

  integer m
  integer n
  integer r8utp_size
  integer value

  if ( n < m ) then
    value = ( n * ( n + 1 ) ) / 2
  else if ( m == n ) then
    value = ( m * ( m + 1 ) ) / 2
  else
    value = ( m * ( m + 1 ) / 2 ) + ( n - m ) * m
  end if

  r8utp_size = value

  return
end
subroutine r8utp_to_r8ge ( m, n, a_utp, a_ge )

!*****************************************************************************80
!
!! r8utp_to_r8ge() copies an R8UTP matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R8UTP storage format is appropriate for an upper triangular
!    matrix.  Only the upper triangle of the matrix is stored,
!    by successive partial columns, in an array which contains
!    (A11,A12,A22,A13,A23,A33,A14,...,AMN).
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
!    15 August 2022
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, N, the order of the matrix.
!
!    real ( kind = rk ) A_UTP(*), the R8UTP matrix.
!
!  Output:
!
!    real ( kind = rk ) A_GE(M,N), the R8GE matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a_ge(m,n)
  real ( kind = rk ) a_utp(*)
  integer i
  integer j
  integer k

  a_ge(1:m,1:n) = 0.0D+00

  k = 0
  do j = 1, n
    do i = 1, min ( j, m )
      k = k + 1
      a_ge(i,j) = a_utp(k)
    end do
  end do

  return
end
subroutine r8utp_zeros ( m, n, a )

!*****************************************************************************80
!
!! r8utp_zeros() zeroes an R8UTP matrix.
!
!  Discussion:
!
!    The R8UTP storage format is appropriate for an upper triangular
!    matrix.  Only the upper triangle of the matrix is stored,
!    by successive partial columns, in an array which contains
!    (A11,A12,A22,A13,A23,A33,A14,...,AMN).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 August 2022
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, N, the order of the matrix.
!
!  Output:
!
!    real ( kind = rk ) A(*), the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(*)
  integer m
  integer mn
  integer n
  integer r8utp_size

  mn = r8utp_size ( m, n )

  a(1:mn) = 0.0D+00

  return
end

