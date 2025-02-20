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
subroutine r8pp_det ( n, a_lu, det )

!*****************************************************************************80
!
!! R8PP_DET computes the determinant of an R8PP matrix factored by R8PP_FA.
!
!  Discussion:
!
!    The R8PP storage format is used for a symmetric positive
!    definite matrix.  Only the upper triangle of the matrix is stored,
!    by successive partial columns, in an array of length (N*(N+1))/2,
!    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
!
!    R8PP storage is used by LINPACK and LAPACK.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    13 June 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real ( kind = rk ) A_LU((N*(N+1))/2), the LU factors from R8PO_FA.
!
!    Output, real ( kind = rk ) DET, the determinant of A.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a_lu((n*(n+1))/2)
  real ( kind = rk ) det
  integer i
  integer k

  det = 1.0D+00

  k = 0
  do i = 1, n
    k = k + i
    det = det * a_lu(k)
  end do

  det = det * det

  return
end
subroutine r8pp_dif2 ( n, a )

!*****************************************************************************80
!
!! R8PP_DIF2 sets up an R8PP second difference matrix.
!
!  Discussion:
!
!    The R8PP storage format is appropriate for a symmetric positive
!    definite matrix.  Only the upper triangle of the matrix is stored,
!    by successive partial columns, in an array of length (N*(N+1))/2,
!    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    13 June 2016
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
!    Output, real ( kind = rk ) A((N*(N+1))/2), the R8PP matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a((n*(n+1))/2)
  integer i
  integer j
  integer k

  k = 0
  do j = 1, n
    do i = 1, j - 2
      k = k + 1
      a(k) = 0.0D+00
    end do
    if ( 1 < j ) then
      k = k + 1
      a(k) = - 1.0D+00
    end if
    k = k + 1
    a(k) = 2.0D+00
  end do

  return
end
subroutine r8pp_fa ( n, a, info )

!*****************************************************************************80
!
!! R8PP_FA factors an R8PP matrix.
!
!  Discussion:
!
!    The R8PP storage format is appropriate for a symmetric positive
!    definite matrix.  Only the upper triangle of the matrix is stored,
!    by successive partial columns, in an array of length (N*(N+1))/2,
!    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
!
!    R8PP storage is used by LINPACK and LAPACK.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 November 2008
!
!  Author:
!
!    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, (Society for Industrial and Applied Mathematics),
!    3600 University City Science Center,
!    Philadelphia, PA, 19104-2688.
!    ISBN 0-89871-172-X
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real ( kind = rk ) A((N*(N+1))/2).  On input, an R8PP matrix.
!    On output, an upper triangular matrix R, stored in packed form, 
!    so that A = R'*R.
!
!    Output, integer INFO, error flag.
!    0, for normal return.
!    K, if the leading minor of order K is not positive definite.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a((n*(n+1))/2)
  integer i
  integer info
  integer j
  integer jj
  integer k
  integer kj
  integer kk
  real ( kind = rk ) s
  real ( kind = rk ) t

  info = 0
  jj = 0

  do j = 1, n

    s = 0.0D+00
    kj = jj
    kk = 0

    do k = 1, j - 1

      kj = kj + 1
      t = a(kj)
      do i = 1, k - 1
        t = t - a(kk+i) * a(jj+i)
      end do
      kk = kk + k
      t = t / a(kk)
      a(kj) = t
      s = s + t * t

    end do

    jj = jj + j
    s = a(jj) - s

    if ( s <= 0.0D+00 ) then
      info = j
      return
    end if

    a(jj) = sqrt ( s )

  end do

  return
end
subroutine r8pp_indicator ( n, a )

!*****************************************************************************80
!
!! R8PP_INDICATOR sets up an R8PP indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R8PP storage format is appropriate for a symmetric positive
!    definite matrix.  Only the upper triangle of the matrix is stored,
!    by successive partial columns, in an array of length (N*(N+1))/2,
!    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
!
!    R8PP storage is used by LINPACK and LAPACK.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 February 2004
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
!    Output, real ( kind = rk ) A((N*(N+1))/2), the R8PP matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a((n*(n+1))/2)
  integer fac
  integer i
  integer i4_log_10
  integer j
  integer k

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  k = 0
  do j = 1, n
    do i = 1, j
      k = k + 1
      a(k) = real ( fac * i + j, kind = rk )
    end do
  end do

  return
end
subroutine r8pp_mv ( n, a, x, b )

!*****************************************************************************80
!
!! R8PP_MV multiplies an R8PP matrix by an R8VEC.
!
!  Discussion:
!
!    The R8PP storage format is appropriate for a symmetric positive
!    definite matrix.  Only the upper triangle of the matrix is stored,
!    by successive partial columns, in an array of length (N*(N+1))/2,
!    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
!
!    R8PP storage is used by LINPACK and LAPACK.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real ( kind = rk ) A((N*(N+1))/2), the R8PP matrix.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = rk ) B(N), the product A * x.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a((n*(n+1))/2)
  real ( kind = rk ) b(n)
  integer i
  integer j
  integer k
  real ( kind = rk ) x(n)

  do i = 1, n
    b(i) = 0.0D+00
    do j = 1, i - 1
      k = j + ( i * ( i - 1 ) ) / 2
      b(i) = b(i) + a(k) * x(j)
    end do
    do j = i, n
      k = i + ( j * ( j - 1 ) ) / 2
      b(i) = b(i) + a(k) * x(j)
    end do
  end do

  return
end
subroutine r8pp_print ( n, a, title )

!*****************************************************************************80
!
!! R8PP_PRINT prints an R8PP matrix.
!
!  Discussion:
!
!    The R8PP storage format is appropriate for a symmetric positive
!    definite matrix.  Only the upper triangle of the matrix is stored,
!    by successive partial columns, in an array of length (N*(N+1))/2,
!    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
!
!    R8PP storage is used by LINPACK and LAPACK.
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
!    Input, real ( kind = rk ) A((N*(N+1))/2), the R8PP matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a((n*(n+1))/2)
  character ( len = * ) title

  call r8pp_print_some ( n, a, 1, 1, n, n, title )

  return
end
subroutine r8pp_print_some ( n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8PP_PRINT_SOME prints some of an R8PP matrix.
!
!  Discussion:
!
!    The R8PP storage format is appropriate for a symmetric positive
!    definite matrix.  Only the upper triangle of the matrix is stored,
!    by successive partial columns, in an array of length (N*(N+1))/2,
!    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
!
!    R8PP storage is used by LINPACK and LAPACK.
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
!    Input, real ( kind = rk ) A((N*(N+1))/2), the R8PP matrix.
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

  real ( kind = rk ) a((n*(n+1))/2)
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
          aij = a(i+(j*(j-1))/2)
        else
          aij = a(j+(i*(i-1))/2)
        end if

        write ( ctemp(j2), '(g14.6)' ) aij

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8pp_random ( n, seed, a )

!*****************************************************************************80
!
!! R8PP_RANDOM randomizes an R8PP matrix.
!
!  Discussion:
!
!    The R8PP storage format is appropriate for a symmetric positive
!    definite matrix.  Only the upper triangle of the matrix is stored,
!    by successive partial columns, in an array of length (N*(N+1))/2,
!    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
!
!    R8PP storage is used by LINPACK and LAPACK.
!
!    The matrix is computed by setting a "random" upper triangular
!    Cholesky factor R, and then computing A = R'*R.
!    The randomness is limited by the fact that all the entries of
!    R will be between 0 and 1.  A truly random R is only required
!    to have positive entries on the diagonal.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 October 1998
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
!    Input/output, integer SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = rk ) A((N*(N+1))/2), the R8PP matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a((n*(n+1))/2)
  real ( kind = rk ) r8_uniform_01
  integer i
  integer ii
  integer ij
  integer ik
  integer j
  integer k
  integer kj
  integer seed

  a(1:(n*(n+1))/2) = 0.0D+00

  do i = n, 1, -1
!
!  Set row I of R.
!
    do j = i, n
      ij = i + ( j * ( j - 1 ) ) / 2
      a(ij) = r8_uniform_01 ( seed )
    end do
!
!  Consider element J of row I, last to first.
!
    do j = n, i, -1
!
!  Add multiples of row I to lower elements of column J.
!
      ij = i + ( j * ( j - 1 ) ) / 2

      do k = i + 1, j
        kj = k + ( j * ( j - 1 ) ) / 2
        ik = i + ( k * ( k - 1 ) ) / 2
        a(kj) = a(kj) + a(ik) * a(ij)
      end do
!
!  Reset element J.
!
      ii = i + ( i * ( i - 1 ) ) / 2
      a(ij) = a(ii) * a(ij)

    end do
  end do

  return
end
subroutine r8pp_sl ( n, a_lu, b )

!*****************************************************************************80
!
!! R8PP_SL solves an R8PP system factored by R8PP_FA.
!
!  Discussion:
!
!    The R8PP storage format is appropriate for a symmetric positive
!    definite matrix.  Only the upper triangle of the matrix is stored,
!    by successive partial columns, in an array of length (N*(N+1))/2,
!    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
!
!    R8PP storage is used by LINPACK and LAPACK.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 November 2008
!
!  Author:
!
!    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, (Society for Industrial and Applied Mathematics),
!    3600 University City Science Center,
!    Philadelphia, PA, 19104-2688.
!    ISBN 0-89871-172-X
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real ( kind = rk ) A_LU((N*(N+1))/2), the LU factors from R8PP_FA.
!
!    Input/output, real ( kind = rk ) B(N).  On input, the right hand side.
!    On output, the solution.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a_lu((n*(n+1))/2)
  real ( kind = rk ) b(n)
  integer i
  integer k
  integer kk
  real ( kind = rk ) t

  kk = 0

  do k = 1, n
    t = 0.0D+00
    do i = 1, k - 1
      t = t + a_lu(kk+i) * b(i)
    end do
    kk = kk + k
    b(k) = ( b(k) - t ) / a_lu(kk)
  end do

  do k = n, 1, -1
    b(k) = b(k) / a_lu(kk)
    kk = kk - k
    t = -b(k)
    do i = 1, k - 1
      b(i) = b(i) + t * a_lu(kk+i)
    end do
  end do

  return
end
subroutine r8pp_to_r8ge ( n, a, b )

!*****************************************************************************80
!
!! R8PP_TO_R8GE copies an R8PP matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R8PP storage format is appropriate for a symmetric positive
!    definite matrix.  Only the upper triangle of the matrix is stored,
!    by successive partial columns, in an array of length (N*(N+1))/2,
!    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
!
!    R8PP storage is used by LINPACK and LAPACK.
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
!    13 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real ( kind = rk ) A((N*(N+1))/2), the R8PP matrix.
!
!    Output, real ( kind = rk ) B(N,N), the R8GE matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a((n*(n+1))/2)
  real ( kind = rk ) b(n,n)
  integer i
  integer j

  do i = 1, n
    do j = 1, n
      if ( i <= j ) then
        b(i,j) = a(i+(j*(j-1))/2)
      else
        b(i,j) = a(j+(i*(i-1))/2)
      end if
    end do
  end do

  return
end
subroutine r8pp_zeros ( n, a )

!*****************************************************************************80
!
!! R8PP_ZEROS zeroes an R8PP matrix.
!
!  Discussion:
!
!    The R8PP storage format is appropriate for a symmetric positive
!    definite matrix.  Only the upper triangle of the matrix is stored,
!    by successive partial columns, in an array of length (N*(N+1))/2,
!    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
!
!    R8PP storage is used by LINPACK and LAPACK.
!
!    The matrix is computed by setting a "random" upper triangular
!    Cholesky factor R, and then computing A = R'*R.
!    The randomness is limited by the fact that all the entries of
!    R will be between 0 and 1.  A truly random R is only required
!    to have positive entries on the diagonal.
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
!    Output, real ( kind = rk ) A((N*(N+1))/2), the R8PP matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a((n*(n+1))/2)

  a(1:(n*(n+1))/2) = 0.0D+00

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

