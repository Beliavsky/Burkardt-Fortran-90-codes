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
subroutine r85_dif2 ( n, a )

!*****************************************************************************80
!
!! R85_DIF2 sets up an R85 second difference matrix.
!
!  Discussion:
!
!    The R85 storage format represents a pentadiagonal matrix as a 5 
!    by N array, in which each row corresponds to a diagonal, and 
!    column locations are preserved.  Thus, the original matrix is 
!    "collapsed" vertically into the array.
!
!  Example:
!
!    Here is how an R85 matrix of order 6 would be stored:
!
!       *   *  A13 A24 A35 A46
!       *  A12 A23 A34 A45 A56
!      A11 A22 A33 A44 A55 A66
!      A21 A32 A43 A54 A65  *
!      A31 A42 A53 A64  *   *
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 July 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Output, real ( kind = rk ) A(5,N), the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(5,n)
  integer i
  integer j

  a(1:5,1:n) = 0.0D+00

  do j = 2, n
    i = j - 1
    a(2,j) = - 1.0D+00
  end do

  do j = 1, n
    i = j
    a(3,j) = 2.0D+00
  end do

  do j = 1, n - 1
    i = j + 1
    a(4,j) = - 1.0D+00
  end do

  return
end
subroutine r85_indicator ( n, a )

!*****************************************************************************80
!
!! R85_INDICATOR sets up an R85 indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R85 storage format represents a pentadiagonal matrix as a 5 
!    by N array, in which each row corresponds to a diagonal, and 
!    column locations are preserved.  Thus, the original matrix is 
!    "collapsed" vertically into the array.
!
!  Example:
!
!    Here is how an R85 matrix of order 6 would be stored:
!
!       *   *  A13 A24 A35 A46
!       *  A12 A23 A34 A45 A56
!      A11 A22 A33 A44 A55 A66
!      A21 A32 A43 A54 A65  *
!      A31 A42 A53 A64  *   *
!
!    Here are the values as stored in an indicator matrix:
!
!      00 00 13 24 35 46
!      00 12 23 34 45 56
!      11 22 33 44 55 66
!      21 32 43 54 65 00
!      31 42 53 64 00 00
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
!
!    Output, real ( kind = rk ) A(5,N), the indicator matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(5,n)
  integer fac
  integer i
  integer i4_log_10
  integer j

  a(1:5,1:n) = 0.0D+00

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  do j = 3, n
    i = j - 2
    a(1,j) = real ( fac * i + j, kind = rk )
  end do

  do j = 2, n
    i = j - 1
    a(2,j) = real ( fac * i + j, kind = rk )
  end do

  do j = 1, n
    i = j
    a(3,j) = real ( fac * i + j, kind = rk )
  end do

  do j = 1, n - 1
    i = j + 1
    a(4,j) = real ( fac * i + j, kind = rk )
  end do

  do j = 1, n - 2
    i = j + 2
    a(5,j) = real ( fac * i + j, kind = rk )
  end do

  return
end
subroutine r85_np_fs ( n, a, b, x )

!*****************************************************************************80
!
!! R85_NP_FS factors and solves an R85 linear system.
!
!  Discussion:
!
!    The R85 storage format represents a pentadiagonal matrix as a 5 
!    by N array, in which each row corresponds to a diagonal, and 
!    column locations are preserved.  Thus, the original matrix is 
!    "collapsed" vertically into the array.
!
!    The factorization algorithm requires that each diagonal entry be nonzero.
!
!    No pivoting is performed, and therefore the algorithm may fail
!    in simple cases where the matrix is not singular.
!
!  Example:
!
!    Here is how an R85 matrix of order 6 would be stored:
!
!       *   *  A13 A24 A35 A46
!       *  A12 A23 A34 A45 A56
!      A11 A22 A33 A44 A55 A66
!      A21 A32 A43 A54 A65  *
!      A31 A42 A53 A64  *   *
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 September 2003
!
!  Author:
!
!    Original FORTRAN77 version by Ward Cheney, David Kincaid.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Ward Cheney, David Kincaid,
!    Numerical Mathematics and Computing,
!    Brooks-Cole Publishing, 2004,
!    ISBN: 0534201121.
!
!  Parameters:
!
!    Input, integer N, the order of the linear system.
!
!    Input/output, real ( kind = rk ) A(5,N),
!    On input, the pentadiagonal matrix.
!    On output, the data has been overwritten by factorization information.
!
!    Input/output, real ( kind = rk ) B(N).
!    On input, B contains the right hand side of the linear system.
!    On output, B has been overwritten by factorization information.
!
!    Output, real ( kind = rk ) X(N), the solution of the linear system.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(5,n)
  real ( kind = rk ) b(n)
  integer i
  real ( kind = rk ) x(n)
  real ( kind = rk ) xmult

  do i = 1, n
    if ( a(3,i) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R85_NP_FS - Fatal error!'
      write ( *, '(a,i8,a)' ) '  A(3,', i, ') = 0.'
      stop 1
    end if
  end do

  do i = 2, n - 1

    xmult = a(2,i) / a(3,i-1)
    a(3,i) = a(3,i) - xmult * a(4,i-1)
    a(4,i) = a(4,i) - xmult * a(5,i-1)

    b(i) = b(i) - xmult * b(i-1)

    xmult = a(1,i+1) / a(3,i-1)
    a(2,i+1) = a(2,i+1) - xmult * a(4,i-1)
    a(3,i+1) = a(3,i+1) - xmult * a(5,i-1)

    b(i+1) = b(i+1) - xmult * b(i-1)

  end do

  xmult = a(2,n) / a(3,n-1)
  a(3,n) = a(3,n) - xmult * a(4,n-1)

  x(n) = ( b(n) - xmult * b(n-1) ) / a(3,n)
  x(n-1) = ( b(n-1) - a(4,n-1) * x(n) ) / a(3,n-1)

  do i = n - 2, 1, -1
    x(i) = ( b(i) - a(4,i) * x(i+1) - a(5,i) * x(i+2) ) / a(3,i)
  end do

  return
end
subroutine r85_mtv ( n, a, x, b )

!*****************************************************************************80
!
!! R85_MTV multiplies an R8VEC by an R85 matrix.
!
!  Discussion:
!
!    The R85 storage format represents a pentadiagonal matrix as a 5 
!    by N array, in which each row corresponds to a diagonal, and 
!    column locations are preserved.  Thus, the original matrix is 
!    "collapsed" vertically into the array.
!
!  Example:
!
!    Here is how an R85 matrix of order 6 would be stored:
!
!       *   *  A13 A24 A35 A46
!       *  A12 A23 A34 A45 A56
!      A11 A22 A33 A44 A55 A66
!      A21 A32 A43 A54 A65  *
!      A31 A42 A53 A64  *   *
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the linear system.
!
!    Input, real ( kind = rk ) A(5,N), the matrix.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied by A'.
!
!    Output, real ( kind = rk ) B(N), the product A' * x.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(5,n)
  real ( kind = rk ) b(n)
  real ( kind = rk ) x(n)

  b(1:n)   =            a(3,1:n)   * x(1:n)
  b(2:n)   = b(2:n)   + a(4,1:n-1) * x(1:n-1)
  b(3:n)   = b(3:n)   + a(5,1:n-2) * x(1:n-2)
  b(1:n-1) = b(1:n-1) + a(2,2:n)   * x(2:n)
  b(1:n-2) = b(1:n-2) + a(1,3:n)   * x(3:n)

  return
end
subroutine r85_mv ( n, a, x, b )

!*****************************************************************************80
!
!! R85_MV multiplies an R85 matrix by an R8VEC.
!
!  Discussion:
!
!    The R85 storage format represents a pentadiagonal matrix as a 5 
!    by N array, in which each row corresponds to a diagonal, and 
!    column locations are preserved.  Thus, the original matrix is 
!    "collapsed" vertically into the array.
!
!  Example:
!
!    Here is how an R85 matrix of order 6 would be stored:
!
!       *   *  A13 A24 A35 A46
!       *  A12 A23 A34 A45 A56
!      A11 A22 A33 A44 A55 A66
!      A21 A32 A43 A54 A65  *
!      A31 A42 A53 A64  *   *
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the linear system.
!
!    Input, real ( kind = rk ) A(5,N), the matrix.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = rk ) B(N), the product A * x.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(5,n)
  real ( kind = rk ) b(n)
  real ( kind = rk ) x(n)

  b(1:n)   =            a(3,1:n)   * x(1:n)

  b(3:n)   = b(3:n)   + a(1,3:n)   * x(1:n-2)
  b(2:n)   = b(2:n)   + a(2,2:n)   * x(1:n-1)
  b(1:n-1) = b(1:n-1) + a(4,1:n-1) * x(2:n)
  b(1:n-2) = b(1:n-2) + a(5,1:n-2) * x(3:n)

  return
end
subroutine r85_print ( n, a, title )

!*****************************************************************************80
!
!! R85_PRINT prints an R85 matrix.
!
!  Discussion:
!
!    The R85 storage format represents a pentadiagonal matrix as a 5 
!    by N array, in which each row corresponds to a diagonal, and 
!    column locations are preserved.  Thus, the original matrix is 
!    "collapsed" vertically into the array.
!
!  Example:
!
!    Here is how an R85 matrix of order 6 would be stored:
!
!       *   *  A13 A24 A35 A46
!       *  A12 A23 A34 A45 A56
!      A11 A22 A33 A44 A55 A66
!      A21 A32 A43 A54 A65  *
!      A31 A42 A53 A64  *   *
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real ( kind = rk ) A(5,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(5,n)
  character ( len = * ) title

  call r85_print_some ( n, a, 1, 1, n, n, title )

  return
end
subroutine r85_print_some ( n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R85_PRINT_SOME prints some of an R85 matrix.
!
!  Discussion:
!
!    The R85 storage format represents a pentadiagonal matrix as a 5 
!    by N array, in which each row corresponds to a diagonal, and 
!    column locations are preserved.  Thus, the original matrix is 
!    "collapsed" vertically into the array.
!
!  Example:
!
!    Here is how an R85 matrix of order 6 would be stored:
!
!       *   *  A13 A24 A35 A46
!       *  A12 A23 A34 A45 A56
!      A11 A22 A33 A44 A55 A66
!      A21 A32 A43 A54 A65  *
!      A31 A42 A53 A64  *   *
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
!
!    Input, real ( kind = rk ) A(5,N), the matrix.
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

  real ( kind = rk ) a(5,n)
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
    i2lo = max ( i2lo, j2lo - 2 )

    i2hi = min ( ihi, n )
    i2hi = min ( i2hi, j2hi + 2 )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( 2 < i - j .or. 2 < j - i ) then
          ctemp(j2) = '              '
        else if ( j == i + 2 ) then
          write ( ctemp(j2), '(g14.6)' ) a(1,j)
        else if ( j == i + 1 ) then
          write ( ctemp(j2), '(g14.6)' ) a(2,j)
        else if ( j == i ) then
          write ( ctemp(j2), '(g14.6)' ) a(3,j)
        else if ( j == i - 1 ) then
          write ( ctemp(j2), '(g14.6)' ) a(4,j)
        else if ( j == i - 2 ) then
          write ( ctemp(j2), '(g14.6)' ) a(5,j)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r85_random ( n, seed, a )

!*****************************************************************************80
!
!! R85_RANDOM randomizes an R85 matrix.
!
!  Discussion:
!
!    The R85 storage format represents a pentadiagonal matrix as a 5 
!    by N array, in which each row corresponds to a diagonal, and 
!    column locations are preserved.  Thus, the original matrix is 
!    "collapsed" vertically into the array.
!
!  Example:
!
!    Here is how an R85 matrix of order 6 would be stored:
!
!       *   *  A13 A24 A35 A46
!       *  A12 A23 A34 A45 A56
!      A11 A22 A33 A44 A55 A66
!      A21 A32 A43 A54 A65  *
!      A31 A42 A53 A64  *   *
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
!    Input, integer N, the order of the linear system.
!
!    Input/output, integer SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = rk ) A(5,N), the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(5,n)
  real ( kind = rk ) r8_uniform_01
  integer j
  integer seed

  a(1:5,1:n) = 0.0D+00

  do j = 3, n
    a(1,j) = r8_uniform_01 ( seed )
  end do

  do j = 2, n
    a(2,j) = r8_uniform_01 ( seed )
  end do

  do j = 1, n
    a(3,j) = r8_uniform_01 ( seed )
  end do

  do j = 1, n - 1
    a(4,j) = r8_uniform_01 ( seed )
  end do

  do j = 1, n - 2
    a(5,j) = r8_uniform_01 ( seed )
  end do

  return
end
subroutine r85_to_r8ge ( n, a, b )

!*****************************************************************************80
!
!! R85_TO_R8GE copies an R85 matrix into an R8GE matrix.
!
!  Discussion:
!
!    The R85 storage format represents a pentadiagonal matrix as a 5 
!    by N array, in which each row corresponds to a diagonal, and 
!    column locations are preserved.  Thus, the original matrix is 
!    "collapsed" vertically into the array.
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
!    Here is how an R85 matrix of order 6 would be stored:
!
!       *   *  A13 A24 A35 A46
!       *  A12 A23 A34 A45 A56
!      A11 A22 A33 A44 A55 A66
!      A21 A32 A43 A54 A65  *
!      A31 A42 A53 A64  *   *
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real ( kind = rk ) A(5,N), the R85 matrix.
!
!    Output, real ( kind = rk ) A(N,N), the R8GE matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(5,n)
  real ( kind = rk ) b(n,n)
  integer i
  integer j

  b(1:n,1:n) = 0.0D+00

  do i = 1, n
    do j = 1, n

      if ( j == i - 2 ) then
        b(i,j) = a(1,i)
      else if ( j == i - 1 ) then
        b(i,j) = a(2,i)
      else if ( i == j ) then
        b(i,j) = a(3,i)
      else if ( j == i + 1 ) then
        b(i,j) = a(4,i)
      else if ( j == i + 2 ) then
        b(i,j) = a(5,i)
      end if

    end do
  end do

  return
end
subroutine r85_zeros ( n, a )

!*****************************************************************************80
!
!! R85_ZEROS zeroes an R85 matrix.
!
!  Discussion:
!
!    The R85 storage format represents a pentadiagonal matrix as a 5 
!    by N array, in which each row corresponds to a diagonal, and 
!    column locations are preserved.  Thus, the original matrix is 
!    "collapsed" vertically into the array.
!
!  Example:
!
!    Here is how an R85 matrix of order 6 would be stored:
!
!       *   *  A13 A24 A35 A46
!       *  A12 A23 A34 A45 A56
!      A11 A22 A33 A44 A55 A66
!      A21 A32 A43 A54 A65  *
!      A31 A42 A53 A64  *   *
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
!    Input, integer N, the order of the linear system.
!
!    Output, real ( kind = rk ) A(5,N), the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(5,n)

  a(1:5,1:n) = 0.0D+00

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

