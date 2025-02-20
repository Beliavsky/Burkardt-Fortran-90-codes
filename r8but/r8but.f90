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
subroutine r8but_det ( n, mu, a, det )

!*****************************************************************************80
!
!! R8BUT_DET computes the determinant of an R8BUT matrix.
!
!  Discussion:
!
!    The R8BUT storage format is used for a banded upper triangular matrix.
!    The matrix is assumed to be zero above the MU-th superdiagonal.
!    The matrix is stored in an MU+1 by N array.
!    Columns are preserved.
!
!    The diagonal is stored in row MU+1 of the array.
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
!
!  Example:
!
!    N = 5, MU = 2
!
!    A11 A12 A13   0   0
!      0 A22 A23 A24   0
!      0   0 A33 A34 A35
!      0   0   0 A44 A45
!      0   0   0   0 A55
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
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer MU, the upper bandwidth.
!
!    Input, real ( kind = rk ) A(MU+1,N), the R8BUT matrix.
!
!    Output, real ( kind = rk ) DET, the determinant of A.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer mu
  integer n

  real ( kind = rk ) a(mu+1,n)
  real ( kind = rk ) det

  det = product ( a(mu+1,1:n) )

  return
end
subroutine r8but_indicator ( n, mu, a )

!*****************************************************************************80
!
!! R8BUT_INDICATOR sets up an R8BUT indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R8BUT storage format is used for a banded upper triangular matrix.
!    The matrix is assumed to be zero above the MU-th superdiagonal.
!    The matrix is stored in an MU+1 by N array.
!    Columns are preserved.
!
!    The diagonal is stored in row MU+1 of the array.
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
!
!  Example:
!
!    N = 5, MU = 2
!
!    A11 A12 A13   0   0
!      0 A22 A23 A24   0
!      0   0 A33 A34 A35
!      0   0   0 A44 A45
!      0   0   0   0 A55
!                --- ---
!                    ---
!
!    The indicator matrix is stored as:
!
!       0   0  13  24  35
!       0  12  23  34  45
!      11  22  33  44  55
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
!  Parameters:
!
!    Input, integer N, the number of columns of the matrix.
!
!    Input, integer MU, the upper bandwidth.
!
!    Output, real ( kind = rk ) A(MU+1,N), the R8BUT matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer mu
  integer n

  real ( kind = rk ) a(mu+1,n)
  integer fac
  integer i
  integer i4_log_10
  integer j

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  do i = 1, n
    do j = i, min ( n, i + mu )
      a(i-j+mu+1,j) = real ( fac * i + j, kind = rk )
    end do
  end do

  do i = 1, mu
    do j = 1, mu+1-i
      a(i,j) = 0.0D+00
    end do
  end do

  return
end
subroutine r8but_mtv ( n, mu, a, x, b )

!*****************************************************************************80
!
!! R8BUT_MTV multiplies an R8VECr by an R8BUT matrix.
!
!  Discussion:
!
!    The R8BUT storage format is used for a banded upper triangular matrix.
!    The matrix is assumed to be zero above the MU-th superdiagonal.
!    The matrix is stored in an MU+1 by N array.
!    Columns are preserved.
!
!    The diagonal is stored in row MU+1 of the array.
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
!
!  Example:
!
!    N = 5, MU = 2
!
!    A11 A12 A13   0   0
!      0 A22 A23 A24   0
!      0   0 A33 A34 A35
!      0   0   0 A44 A45
!      0   0   0   0 A55
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
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer MU, the upper bandwidth.
!
!    Input, real ( kind = rk ) A(MU+1,N), the R8BUT matrix.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = rk ) B(N), the product X*A.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer mu
  integer n

  real ( kind = rk ) a(mu+1,n)
  real ( kind = rk ) b(n)
  integer i
  integer ilo
  real ( kind = rk ) x(n)

  do i = 1, n
    ilo = max ( 1, i - mu )
    b(i) = sum ( x(ilo:i) * a(ilo-i+mu+1:mu+1,i) )
  end do

  return
end
subroutine r8but_mv ( n, mu, a, x, b )

!*****************************************************************************80
!
!! R8BUT_MV multiplies an R8BUT matrix by an R8VEC.
!
!  Discussion:
!
!    The R8BUT storage format is used for a banded upper triangular matrix.
!    The matrix is assumed to be zero above the MU-th superdiagonal.
!    The matrix is stored in an MU+1 by N array.
!    Columns are preserved.
!
!    The diagonal is stored in row MU+1 of the array.
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
!
!  Example:
!
!    N = 5, MU = 2
!
!    A11 A12 A13   0   0
!      0 A22 A23 A24   0
!      0   0 A33 A34 A35
!      0   0   0 A44 A45
!      0   0   0   0 A55
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
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer MU, the upper bandwidth.
!
!    Input, real ( kind = rk ) A(MU+1,N), the R8BUT matrix.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = rk ) B(N), the product A * x.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer mu
  integer n

  real ( kind = rk ) a(mu+1,n)
  real ( kind = rk ) b(n)
  integer i
  integer j
  real ( kind = rk ) x(n)

  b(1:n) = 0.0D+00

  do i = 1, n
    do j = i, min ( n, i + mu )
      b(i) = b(i) + a(i-j+mu+1,j) * x(j)
    end do
  end do

  return
end
subroutine r8but_print ( n, mu, a, title )

!*****************************************************************************80
!
!! R8BUT_PRINT prints an R8BUT matrix.
!
!  Discussion:
!
!    The R8BUT storage format is used for a banded upper triangular matrix.
!    The matrix is assumed to be zero above the MU-th superdiagonal.
!    The matrix is stored in an MU+1 by N array.
!    Columns are preserved.
!
!    The diagonal is stored in row MU+1 of the array.
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
!
!  Example:
!
!    N = 5, MU = 2
!
!    A11 A12 A13   0   0
!      0 A22 A23 A24   0
!      0   0 A33 A34 A35
!      0   0   0 A44 A45
!      0   0   0   0 A55
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
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer MU, the upper bandwidth.
!
!    Input, real ( kind = rk ) A(MU+1,N), the R8BUT matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer mu
  integer n

  real ( kind = rk ) a(mu+1,n)
  character ( len = * ) title

  call r8but_print_some ( n, mu, a, 1, 1, n, n, title )

  return
end
subroutine r8but_print_some ( n, mu, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8BUT_PRINT_SOME prints some of an R8BUT matrix.
!
!  Discussion:
!
!    The R8BUT storage format is used for a banded upper triangular matrix.
!    The matrix is assumed to be zero above the MU-th superdiagonal.
!    The matrix is stored in an MU+1 by N array.
!    Columns are preserved.
!
!    The diagonal is stored in row MU+1 of the array.
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
!
!  Example:
!
!    N = 5, MU = 2
!
!    A11 A12 A13   0   0
!      0 A22 A23 A24   0
!      0   0 A33 A34 A35
!      0   0   0 A44 A45
!      0   0   0   0 A55
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
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer MU, the upper bandwidth.
!
!    Input, real ( kind = rk ) A(MU+1,N), the R8BUT matrix.
!
!    Input, integer ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: incx = 5
  integer mu
  integer n

  real ( kind = rk ) a(mu+1,n)
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
    i2hi = min ( i2hi, j2hi + mu )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( i <= j .and. j <= i + mu ) then
          aij = a(i-j+mu+1,j)
          write ( ctemp(j2), '(g14.6)' ) aij
        else
          ctemp(j2) = '              '
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8but_random ( n, mu, seed, a )

!*****************************************************************************80
!
!! R8BUT_RANDOM randomizes an R8BUT matrix.
!
!  Discussion:
!
!    The R8BUT storage format is used for a banded upper triangular matrix.
!    The matrix is assumed to be zero above the MU-th superdiagonal.
!    The matrix is stored in an MU+1 by N array.
!    Columns are preserved.
!
!    The diagonal is stored in row MU+1 of the array.
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
!
!  Example:
!
!    N = 5, MU = 2
!
!    A11 A12 A13   0   0
!      0 A22 A23 A24   0
!      0   0 A33 A34 A35
!      0   0   0 A44 A45
!      0   0   0   0 A55
!                --- ---
!                    ---
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    08 October 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of columns of the matrix.
!
!    Input, integer MU, the upper bandwidth.
!
!    Input/output, integer SEED, a seed for the random number 
!    generator.
!
!    Output, real ( kind = rk ) A(MU+1,N), the R8BUT matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer mu
  integer n

  real ( kind = rk ) a(mu+1,n)
  integer i
  integer j
  integer k
  real ( kind = rk ) r8_uniform_01
  integer seed

  do j = 1, n
    do i = max ( 1, j - mu ), j
      k = i - j + mu + 1
      a(k,j) = r8_uniform_01 ( seed )
    end do
  end do

  return
end
subroutine r8but_sl ( n, mu, a, b )

!*****************************************************************************80
!
!! R8BUT_SL solves A*x=b, where A is an R8BUT matrix.
!
!  Discussion:
!
!    The R8BUT storage format is used for a banded upper triangular matrix.
!    The matrix is assumed to be zero above the MU-th superdiagonal.
!    The matrix is stored in an MU+1 by N array.
!    Columns are preserved.
!
!    The diagonal is stored in row MU+1 of the array.
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
!
!  Example:
!
!    N = 5, MU = 2
!
!    A11 A12 A13   0   0
!      0 A22 A23 A24   0
!      0   0 A33 A34 A35
!      0   0   0 A44 A45
!      0   0   0   0 A55
!                --- ---
!                    ---
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    08 October 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer MU, the upper bandwidth.
!
!    Input, real ( kind = rk ) A(MU+1,N), the R8BUT matrix.
!
!    Input/output, real ( kind = rk ) B(N).
!    On input, the right hand side.
!    On output, the solution vector.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer mu
  integer n

  real ( kind = rk ) a(mu+1,n)
  real ( kind = rk ) b(n)
  integer j
  integer jlo

  do j = n, 1, -1
    b(j) = b(j) / a(j-j+mu+1,j)
    jlo = max ( 1, j - mu )
    b(jlo:j-1) = b(jlo:j-1) - a(jlo-j+mu+1:j-1-j+mu+1,j) * b(j)
  end do

  return
end
subroutine r8but_slt ( n, mu, a, b )

!*****************************************************************************80
!
!! R8BUT_SLT solves A'*x=b, where A is an R8BUT matrix.
!
!  Discussion:
!
!    The R8BUT storage format is used for a banded upper triangular matrix.
!    The matrix is assumed to be zero above the MU-th superdiagonal.
!    The matrix is stored in an MU+1 by N array.
!    Columns are preserved.
!
!    The diagonal is stored in row MU+1 of the array.
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
!
!  Example:
!
!    N = 5, MU = 2
!
!    A11 A12 A13   0   0
!      0 A22 A23 A24   0
!      0   0 A33 A34 A35
!      0   0   0 A44 A45
!      0   0   0   0 A55
!                --- ---
!                    ---
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    08 October 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer MU, the upper bandwidth.
!
!    Input, real ( kind = rk ) A(MU+1,N), the R8BUT matrix.
!
!    Input/output, real ( kind = rk ) B(N).
!    On input, the right hand side.
!    On output, the solution vector.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer mu
  integer n

  real ( kind = rk ) a(mu+1,n)
  real ( kind = rk ) b(n)
  integer i
  integer ihi
  integer j

  do j = 1, n
    b(j) = b(j) / a(j-j+mu+1,j)
    ihi = min ( n, j + mu )
    do i = j + 1, ihi
      b(i) = b(i) - a(j-i+mu+1,i) * b(j)
    end do
  end do

  return
end
subroutine r8but_to_r8ge ( n, mu, a, b )

!*****************************************************************************80
!
!! R8BUT_TO_R8GE copies an R8BUT matrix to an R8GE matrix.
!
!  Discussion:
!
!    The R8BUT storage format is used for a banded upper triangular matrix.
!    The matrix is assumed to be zero above the MU-th superdiagonal.
!    The matrix is stored in an MU+1 by N array.
!    Columns are preserved.
!
!    The diagonal is stored in row MU+1 of the array.
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
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
!    N = 5, MU = 2
!
!    A11 A12 A13   0   0
!      0 A22 A23 A24   0
!      0   0 A33 A34 A35
!      0   0   0 A44 A45
!      0   0   0   0 A55
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
!  Parameters:
!
!    Input, integer N, the order of the matrices.
!    N must be positive.
!
!    Input, integer MU, the upper bandwidth of A.
!    MU must be nonnegative, and no greater than N-1.
!
!    Input, real ( kind = rk ) A(MU+1,N), the R8BUT matrix.
!
!    Output, real ( kind = rk ) B(N,N), the R8GE matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer mu
  integer n

  real ( kind = rk ) a(mu+1,n)
  real ( kind = rk ) b(n,n)
  integer i
  integer j

  do i = 1, n
    do j = 1, n
      if ( i <= j .and. j <= i+mu ) then
        b(i,j) = a(i-j+mu+1,j)
      else
        b(i,j) = 0.0D+00
      end if
    end do
  end do

  return
end
subroutine r8but_zeros ( n, mu, a )

!*****************************************************************************80
!
!! R8BUT_ZEROS zeroes an R8BUT matrix.
!
!  Discussion:
!
!    The R8BUT storage format is used for a banded upper triangular matrix.
!    The matrix is assumed to be zero above the MU-th superdiagonal.
!    The matrix is stored in an MU+1 by N array.
!    Columns are preserved.
!
!    The diagonal is stored in row MU+1 of the array.
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
!
!  Example:
!
!    N = 5, MU = 2
!
!    A11 A12 A13   0   0
!      0 A22 A23 A24   0
!      0   0 A33 A34 A35
!      0   0   0 A44 A45
!      0   0   0   0 A55
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
!  Parameters:
!
!    Input, integer N, the number of columns of the matrix.
!
!    Input, integer MU, the upper bandwidth.
!
!    Output, real ( kind = rk ) A(MU+1,N), the R8BUT matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer mu
  integer n

  real ( kind = rk ) a(mu+1,n)

  a(1:mu+1,1:n) = 0.0D+00

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

