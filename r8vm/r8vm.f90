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
function r8_factorial ( n )

!*****************************************************************************80
!
!! R8_FACTORIAL computes the factorial of N.
!
!  Discussion:
!
!    factorial ( N ) = product ( 1 <= I <= N ) I
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the argument of the factorial function.
!    If N is less than 1, the function value is returned as 1.
!
!    Output, real ( kind = rk ) R8_FACTORIAL, the factorial of N.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r8_factorial
  integer i
  integer n
  real ( kind = rk ) value

  value = 1.0D+00

  do i = 1, n
    value = value * real ( i, kind = rk )
  end do

  r8_factorial = value

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
subroutine r8ge_to_r8vm ( m, n, a_ge, a_vm )

!*****************************************************************************80
!
!! R8GE_TO_R8VM copies an R8GE matrix to an R8VM matrix.
!
!  Discussion:
!
!    The assumption is made that am MxN Vandermonde matrix has been stored in
!    general format (R8GE), and that a copy is to be created in Vandermonde 
!    (R8VM) format.  This means that only the second row of the GE matrix is 
!    examined.
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
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
!    18 August 2015
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
!    Input, real ( kind = rk ) A_GE(N), the R8GE matrix.
!
!    Output, real ( kind = rk ) A_VM(N), the R8VM matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a_ge(m,n)
  real ( kind = rk ) a_vm(n)

  a_vm(1:n) = a_ge(2,1:n)

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
subroutine r8vm_det ( n, a, det )

!*****************************************************************************80
!
!! R8VM_DET computes the determinant of an R8VM matrix.
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
!    20 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of rows and columns of 
!    the matrix.
!
!    Input, real ( kind = rk ) A(N), the R8VM matrix.
!
!    Output, real ( kind = rk ) DET, the determinant of the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) det
  integer i
  integer j

  det = 1.0D+00
  do j = 1, n
    do i = j + 1, n
      det = det * ( a(i) - a(j) )
    end do
  end do

  return
end
subroutine r8vm_indicator ( m, n, a )

!*****************************************************************************80
!
!! R8VM_INDICATOR sets an R8VM indicator matrix.
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
!    The parameter M is not actually needed by this routine.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 August 2015
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
!    Output, real ( kind = rk ) A(N), the R8VM matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  integer j
  integer m

  do j = 1, n
    a(j) = real ( j, kind = rk )
  end do

  return
end
subroutine r8vm_indicator_det ( n, det )

!*****************************************************************************80
!
!! R8VM_INDICATOR_DET returns the determinant of an R8VM indicator matrix.
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
!    25 August 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of rows and columns of 
!    the matrix.
!
!    Output, real ( kind = rk ) DET, the determinant.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) det
  integer i
  integer n
  real ( kind = rk ) r8_factorial

  det = 1.0D+00
  do i = 1, n
    det = det * r8_factorial ( i - 1 )
  end do

  return
end
subroutine r8vm_mtv ( m, n, a, x, b )

!*****************************************************************************80
!
!! R8VM_MTV computes A'*x=b, where A is an R8VM matrix.
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
!    Input, real ( kind = rk ) X(M), the vector to be multiplied by A.
!
!    Output, real ( kind = rk ) B(N), the product A' * x.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) b(n)
  integer i
  integer j
  real ( kind = rk ) x(m)

  do j = 1, n
    b(j) = 0.0D+00
    do i = 1, m
      if ( i == 1 ) then
        b(j) = b(j) + x(i)
      else
        b(j) = b(j) + a(j)**(i-1) * x(i)
      end if
    end do
  end do

  return
end
subroutine r8vm_mv ( m, n, a, x, b )

!*****************************************************************************80
!
!! R8VM_MV computes A*x=b, where A is an R8VM matrix.
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
!    Input, real ( kind = rk ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = rk ) B(M), the product A * x.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) b(m)
  integer i
  integer j
  real ( kind = rk ) x(n)

  do i = 1, m
    b(i) = 0.0D+00
    do j = 1, n
      if ( i == 1 ) then
        b(i) = b(i) + x(j)
      else
        b(i) = b(i) + a(j) ** ( i - 1 ) * x(j)
      end if
    end do
  end do

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
subroutine r8vm_random ( m, n, seed, a )

!*****************************************************************************80
!
!! R8VM_RANDOM randomizes an R8VM matrix.
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
!    The parameter M is not actually needed by this routine.
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
!    Input/output, integer SEED, a seed for the random number 
!    generator.
!
!    Output, real ( kind = rk ) A(N), the R8VM matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  integer m
  integer seed

  call r8vec_uniform_01 ( n, seed, a )

  return
end
subroutine r8vm_sl ( n, a, b, x, info )

!*****************************************************************************80
!
!! R8VM_SL solves the linear system A*x=b where A is an R8VM matrix.
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
!    Vandermonde systems are very close to singularity.  The singularity
!    gets worse as N increases, and as any pair of values defining
!    the matrix get close.  Even a system as small as N = 10 will
!    involve the 9th power of the defining values.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 August 2015
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Gene Golub, Charles Van Loan,
!    Matrix Computations,
!    Third Edition,
!    Johns Hopkins, 1996.
!
!  Parameters:
!
!    Input, integer N, the number of rows and columns of 
!    the matrix.
!
!    Input, real ( kind = rk ) A(N), the R8VM matrix.
!
!    Input, real ( kind = rk ) B(N), the right hand side.
!
!    Output, real ( kind = rk ) X(N), the solution of the linear system.
!
!    Output, integer INFO.
!    0, no error.
!    nonzero, the matrix is singular.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) b(n)
  integer i
  integer info
  integer j
  real ( kind = rk ) x(n)
!
!  Check for explicit singularity.
!
  info = 0

  do j = 1, n - 1
    do i = j + 1, n
      if ( a(i) == a(j) ) then
        info = 1
        return
      end if
    end do
  end do

  x(1:n) = b(1:n)

  do j = 1, n - 1
    do i = n, j + 1, -1
      x(i) = x(i) - a(j) * x(i-1)
    end do
  end do

  do j = n - 1, 1, -1

    do i = j + 1, n
      x(i) = x(i) / ( a(i) - a(i-j) )
    end do

    do i = j, n - 1
      x(i) = x(i) - x(i+1)
    end do

  end do

  return
end
subroutine r8vm_slt ( n, a, b, x, info )

!*****************************************************************************80
!
!! R8VM_SLT solves the linear system A'*x=b where A is an R8VM matrix.
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
!    Vandermonde systems are very close to singularity.  The singularity
!    gets worse as N increases, and as any pair of values defining
!    the matrix get close.  Even a system as small as N = 10 will
!    involve the 9th power of the defining values.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 August 2015
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Gene Golub, Charles Van Loan,
!    Matrix Computations,
!    Third Edition,
!    Johns Hopkins, 1996.
!
!  Parameters:
!
!    Input, integer N, the number of rows and columns of 
!    the matrix.
!
!    Input, real ( kind = rk ) A(N), the R8VM matrix.
!
!    Input, real ( kind = rk ) B(N), the right hand side.
!
!    Output, real ( kind = rk ) X(N), the solution of the linear system.
!
!    Output, integer INFO.
!    0, no error.
!    nonzero, the matrix is singular.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) b(n)
  integer i
  integer info
  integer j
  real ( kind = rk ) x(n)
!
!  Check for explicit singularity.
!
  info = 0

  do j = 1, n - 1
    do i = j + 1, n
      if ( a(i) == a(j) ) then
        info = 1
        return
      end if
    end do
  end do

  x(1:n) = b(1:n)

  do j = 1, n - 1
    do i = n, j + 1, -1
      x(i) = ( x(i) - x(i-1) ) / ( a(i) - a(i-j) )
    end do
  end do

  do j = n - 1, 1, -1
    do i = j, n - 1
      x(i) = x(i) - x(i+1) * a(j)
    end do
  end do

  return
end
subroutine r8vm_to_r8ge ( m, n, a, b )

!*****************************************************************************80
!
!! R8VM_TO_R8GE copies an R8VM matrix to an R8GE matrix.
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
!    23 January 2004
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
!    Output, real ( kind = rk ) B(M,N), the R8GE matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) b(m,n)
  integer i
  integer j

  do i = 1, m
    do j = 1, n
      if ( i == 1 ) then
        b(i,j) = 1.0D+00
      else
        b(i,j) = b(i-1,j) * a(j)
      end if
    end do
  end do

  return
end
subroutine r8vm_zeros ( m, n, a )

!*****************************************************************************80
!
!! R8VM_ZEROS zeroes an R8VM matrix.
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
!    The parameter M is not actually needed by this routine.
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
!    Input, integer M, N, the number of rows and columns of 
!    the matrix.
!
!    Output, real ( kind = rk ) A(N), the R8VM matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  integer m

  a(1:n) = 0.0D+00

  return
end

