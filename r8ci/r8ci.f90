function c8_le_l2 ( x, y )

!*****************************************************************************80
!
!! c8_le_l2() := X <= Y for the L2 norm on C8 values.
!
!  Discussion:
!
!    The L2 norm can be defined here as:
!
!      C8_NORM2(X) = sqrt ( ( real (X) )^2 + ( imag (X) )^2 )
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
!  Parameters:
!
!    Input, complex ( kind = ck ) X, Y, the values to be compared.
!
!    Output, logical C8_LE_L2, is TRUE if X <= Y.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  complex ( kind = ck ) x
  complex ( kind = ck ) y
  logical c8_le_l2

  if ( ( real ( x, kind = rk ) )**2 + ( imag ( x ) )**2 <= &
       ( real ( y, kind = rk ) )**2 + ( imag ( y ) )**2 ) then
    c8_le_l2 = .true.
  else
    c8_le_l2 = .false.
  end if

  return
end
subroutine c8vec_print ( n, a, title )

!*****************************************************************************80
!
!! C8VEC_PRINT prints a C8VEC.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    08 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, complex ( kind = ck ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  integer n

  complex ( kind = ck ) a(n)
  integer i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i8,2g14.6)' ) i, a(i)
  end do

  return
end
subroutine c8vec_sort_a_l2 ( n, x )

!*****************************************************************************80
!
!! C8VEC_SORT_A_L2 ascending sorts a C8VEC by L2 norm.
!
!  Discussion:
!
!    The L2 norm of A+Bi is sqrt ( A**2 + B**2 ).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, length of input array.
!
!    Input/output, complex ( kind = ck ) X(N).
!    On input, an unsorted array.
!    On output, X has been sorted.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  integer n

  logical c8_le_l2
  integer i
  integer indx
  integer isgn
  integer j
  complex ( kind = ck ) t
  complex ( kind = ck ) x(n)

  i = 0
  indx = 0
  isgn = 0
  j = 0

  do

    call sort_heap_external ( n, indx, i, j, isgn )

    if ( 0 < indx ) then

      t    = x(i)
      x(i) = x(j)
      x(j) = t

    else if ( indx < 0 ) then

      if ( c8_le_l2 ( x(i), x(j) ) ) then
        isgn = - 1
      else
        isgn = + 1
      end if

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine c8vec_unity ( n, a )

!*****************************************************************************80
!
!! C8VEC_UNITY returns the N roots of unity as a C8VEC.
!
!  Discussion:
!
!    X(1:N) = exp ( 2 * PI * (0:N-1) / N )
!
!    X(1:N)**N = ( (1,0), (1,0), ..., (1,0) ).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements of A.
!
!    Output, complex ( kind = ck ) A(N), the N roots of unity.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  complex ( kind = ck ) a(n)
  integer i
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) theta

  do i = 1, n
    theta = r8_pi * real ( 2 * ( i - 1 ), kind = rk ) / real ( n, kind = rk )
    a(i) = cmplx ( cos ( theta ), sin ( theta ), kind = ck )
  end do

  return
end
function i4_log_10 ( i )

!*****************************************************************************80
!
!! I4_LOG_10 returns the integer part of the logarithm base 10 of an I4.
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
subroutine r8ci_det ( n, a, det )

!*****************************************************************************80
!
!! R8CI_DET returns the determinant of an R8CI matrix.
!
!  Discussion:
!
!    The R8CI storage format is used for a real N by N circulant matrix.
!    An N by N circulant matrix A has the property that the entries on
!    row I appear again on row I+1, shifted one position to the right,
!    with the final entry of row I appearing as the first of row I+1.
!    The R8CI format simply records the first row of the matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 June 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real ( kind = rk ) A(N), the R8CI matrix.
!
!    Output, real ( kind = rk ) DET, the determinant.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) det
  integer i
  complex ( kind = ck ) lambda(n)
  complex ( kind = ck ) w(n)

  call c8vec_unity ( n, w )

  lambda(1:n) = cmplx ( a(n), 0.0D+00, kind = ck )
  do i = n - 1, 1, -1
    lambda(1:n) = lambda(1:n) * w(1:n) + cmplx ( a(i), 0.0D+00, kind = ck )
  end do

  det = real ( product ( lambda(1:n) ), kind = rk )

  return
end
subroutine r8ci_dif2 ( n, a )

!*****************************************************************************80
!
!! R8CI_DIF2 sets up an R8CI second difference matrix.
!
!  Discussion:
!
!    This is actually a periodic second difference matrix.
!
!    The R8CI storage format is used for a real N by N circulant matrix.
!    An N by N circulant matrix A has the property that the entries on
!    row I appear again on row I+1, shifted one position to the right,
!    with the final entry of row I appearing as the first of row I+1.
!    The R8CI format simply records the first row of the matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 June 2016
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
!    Output, real ( kind = rk ) A(N), the R8CI matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)

  a(1:n) = 0.0D+00

  a(1) = 2.0D+00
  a(2) = -1.0D+00
  a(n) = -1.0D+00

  return
end
subroutine r8ci_eval ( n, a, lambda )

!*****************************************************************************80
!
!! R8CI_EVAL returns the eigenvalues of an R8CI matrix.
!
!  Discussion:
!
!    The R8CI storage format is used for a real N by N circulant matrix.
!    An N by N circulant matrix A has the property that the entries on
!    row I appear again on row I+1, shifted one position to the right,
!    with the final entry of row I appearing as the first of row I+1.
!    The R8CI format simply records the first row of the matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis,
!    Circulant Matrices,
!    Second Edition,
!    Chelsea, 1994,
!    ISBN: 0828403384,
!    LC: QA188.D37.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real ( kind = rk ) A(N), the R8CI matrix.
!
!    Output, complex ( kind = ck ) LAMBDA(N), the eigenvalues.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  integer i
  complex ( kind = ck ) lambda(n)
  complex ( kind = ck ) w(n)

  call c8vec_unity ( n, w )

  lambda(1:n) = cmplx ( a(n), 0.0D+00, kind = ck )
  do i = n - 1, 1, -1
    lambda(1:n) = lambda(1:n) * w(1:n) + cmplx ( a(i), 0.0D+00, kind = ck )
  end do

  call c8vec_sort_a_l2 ( n, lambda )

  return
end
subroutine r8ci_indicator ( n, a )

!*****************************************************************************80
!
!! R8CI_INDICATOR sets up an R8CI indicator matrix.
!
!  Discussion:
!
!    The "indicator matrix" simply has a value like I*10+J at every
!    entry of a dense matrix, or at every entry of a compressed storage
!    matrix for which storage is allocated. 
!
!    The R8CI storage format is used for a real N by N circulant matrix.
!    An N by N circulant matrix A has the property that the entries on
!    row I appear again on row I+1, shifted one position to the right,
!    with the final entry of row I appearing as the first of row I+1.
!    The R8CI format simply records the first row of the matrix.
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
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Output, real ( kind = rk ) A(N), the R8CI matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  integer fac
  integer i
  integer i4_log_10
  integer j

  fac = 10 ** ( i4_log_10 ( n ) + 1 )

  i = 1

  do j = 1, n
    a(j) = real ( fac * i + j, kind = rk )
  end do

  return
end
subroutine r8ci_mtv ( n, a, x, b )

!*****************************************************************************80
!
!! R8CI_MTV multiplies an R8VEC by an R8CI matrix.
!
!  Discussion:
!
!    The R8CI storage format is used for a real N by N circulant matrix.
!    An N by N circulant matrix A has the property that the entries on
!    row I appear again on row I+1, shifted one position to the right,
!    with the final entry of row I appearing as the first of row I+1.
!    The R8CI format simply records the first row of the matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real ( kind = rk ) A(N), the R8CI matrix.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = rk ) B(N), the product A' * X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) b(n)
  integer i
  real ( kind = rk ) x(n)

  do i = 1, n
    b(i) = sum ( a(i:1:-1) * x(1:i) ) &
         + sum ( a(n:i+1:-1) * x(i+1:n) )
  end do

  return
end
subroutine r8ci_mv ( n, a, x, b )

!*****************************************************************************80
!
!! R8CI_MV multiplies an R8CI matrix by an R8VEC.
!
!  Discussion:
!
!    The R8CI storage format is used for a real N by N circulant matrix.
!    An N by N circulant matrix A has the property that the entries on
!    row I appear again on row I+1, shifted one position to the right,
!    with the final entry of row I appearing as the first of row I+1.
!    The R8CI format simply records the first row of the matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real ( kind = rk ) A(N), the R8CI matrix.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = rk ) B(N), the product A * x.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) b(n)
  integer i
  real ( kind = rk ) x(n)

  do i = 1, n
    b(i) = sum ( a(n+2-i:n) * x(1:i-1) ) &
         + sum ( a(1:n+1-i) * x(i:n) )
  end do

  return
end
subroutine r8ci_print ( n, a, title )

!*****************************************************************************80
!
!! R8CI_PRINT prints an R8CI matrix.
!
!  Discussion:
!
!    The R8CI storage format is used for a real N by N circulant matrix.
!    An N by N circulant matrix A has the property that the entries on
!    row I appear again on row I+1, shifted one position to the right,
!    with the final entry of row I appearing as the first of row I+1.
!    The R8CI format simply records the first row of the matrix.
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
!    Input, real ( kind = rk ) A(N), the R8CI matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  character ( len = * ) title

  call r8ci_print_some ( n, a, 1, 1, n, n, title )

  return
end
subroutine r8ci_print_some ( n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8CI_PRINT_SOME prints some of an R8CI matrix.
!
!  Discussion:
!
!    The R8CI storage format is used for a real N by N circulant matrix.
!    An N by N circulant matrix A has the property that the entries on
!    row I appear again on row I+1, shifted one position to the right,
!    with the final entry of row I appearing as the first of row I+1.
!    The R8CI format simply records the first row of the matrix.
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
!    Input, real ( kind = rk ) A(N), the R8CI matrix.
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
          aij = a(j+1-i)
        else
          aij = a(n+j+1-i)
        end if

        write ( ctemp(j2), '(g14.6)' ) aij

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8ci_random ( n, seed, a )

!*****************************************************************************80
!
!! R8CI_RANDOM randomizes an R8CI matrix.
!
!  Discussion:
!
!    The R8CI storage format is used for a real N by N circulant matrix.
!    An N by N circulant matrix A has the property that the entries on
!    row I appear again on row I+1, shifted one position to the right,
!    with the final entry of row I appearing as the first of row I+1.
!    The R8CI format simply records the first row of the matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 May 2003
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
!    Input/output, integer SEED, a seed for the random number 
!    generator.
!
!    Output, real ( kind = rk ) A(N), the R8CI matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  integer seed

  call r8vec_uniform_01 ( n, seed, a )

  return
end
subroutine r8ci_sl ( n, a, b, x, job )

!*****************************************************************************80
!
!! R8CI_SL solves an R8CI system.
!
!  Discussion:
!
!    The R8CI storage format is used for a real N by N circulant matrix.
!    An N by N circulant matrix A has the property that the entries on
!    row I appear again on row I+1, shifted one position to the right,
!    with the final entry of row I appearing as the first of row I+1.
!    The R8CI format simply records the first row of the matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 September 1999
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real ( kind = rk ) A(N), the R8CI matrix.
!
!    Input, real ( kind = rk ) B(N), the right hand side.
!
!    Output, real ( kind = rk ) X(N), the solution of the linear system.
!
!    Input, integer JOB, specifies the system to solve.
!    0, solve A * x = b.
!    nonzero, solve A' * x = b.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) b(n)
  integer i
  integer job
  integer nsub
  real ( kind = rk ) r1
  real ( kind = rk ) r2
  real ( kind = rk ) r3
  real ( kind = rk ) r5
  real ( kind = rk ) r6
  real ( kind = rk ) work(2*n-2)
  real ( kind = rk ) x(n)

  if ( job == 0 ) then
!
!  Solve the system with the principal minor of order 1.
!
    r1 = a(1)
    x(1) = b(1) / r1

    r2 = 0.0D+00
!
!  Recurrent process for solving the system.
!
    do nsub = 2, n
!
!  Compute multiples of the first and last columns of
!  the inverse of the principal minor of order N.
!
      r5 = a(n+2-nsub)
      r6 = a(nsub)

      if ( 2 < nsub ) then

        work(nsub-1) = r2

        do i = 1, nsub - 2
          r5 = r5 + a(n+1-i) * work(nsub-i)
          r6 = r6 + a(i+1) * work(n-1+i)
        end do

      end if

      r2 = - r5 / r1
      r3 = - r6 / r1
      r1 = r1 + r5 * r3

      if ( 2 < nsub ) then

        r6 = work(n)
        work(n-1+nsub-1) = 0.0D+00
        do i = 2, nsub - 1
          r5 = work(n-1+i)
          work(n-1+i) = work(i) * r3 + r6
          work(i) = work(i) + r6 * r2
          r6 = r5
        end do

      end if

      work(n) = r3
!
!  Compute the solution of the system with the principal minor of order NSUB.
!
      r5 = 0.0D+00
      do i = 1, nsub - 1
        r5 = r5 + a(n+1-i) * x(nsub-i)
      end do

      r6 = ( b(nsub) - r5 ) / r1
      x(1:nsub-1) = x(1:nsub-1) + work(n:n+nsub-2) * r6
      x(nsub) = r6

    end do

  else
!
!  Solve the system with the principal minor of order 1.
!
    r1 = a(1)
    x(1) = b(1) / r1

    r2 = 0.0D+00
!
!  Recurrent process for solving the system.
!
    do nsub = 2, n
!
!  Compute multiples of the first and last columns of
!  the inverse of the principal minor of order N.
!
      r5 = a(nsub)
      r6 = a(n+2-nsub)

      if ( 2 < nsub ) then

        work(nsub-1) = r2

        do i = 1, nsub - 2
          r5 = r5 + a(i+1) * work(nsub-i)
          r6 = r6 + a(n+1-i) * work(n-1+i)
        end do

      end if

      r2 = - r5 / r1
      r3 = - r6 / r1
      r1 = r1 + r5 * r3

      if ( 2 < nsub ) then

        r6 = work(n)
        work(n-1+nsub-1) = 0.0D+00
        do i = 2, nsub - 1
          r5 = work(n-1+i)
          work(n-1+i) = work(i) * r3 + r6
          work(i) = work(i) + r6 * r2
          r6 = r5
        end do

      end if

      work(n) = r3
!
!  Compute the solution of the system with the principal minor of order NSUB.
!
      r5 = 0.0D+00
      do i = 1, nsub - 1
        r5 = r5 + a(i+1) * x(nsub-i)
      end do

      r6 = ( b(nsub) - r5 ) / r1
      do i = 1, nsub - 1
        x(i) = x(i) + work(n-1+i) * r6
      end do

      x(nsub) = r6

    end do

  end if

  return
end
subroutine r8ci_to_r8ge ( n, a, b )

!*****************************************************************************80
!
!! R8CI_TO_R8GE copies an R8CI matrix into an R8GE matrix.
!
!  Discussion:
!
!    The R8CI storage format is used for a real N by N circulant matrix.
!    An N by N circulant matrix A has the property that the entries on
!    row I appear again on row I+1, shifted one position to the right,
!    with the final entry of row I appearing as the first of row I+1.
!    The R8CI format simply records the first row of the matrix.
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
!    17 June 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real ( kind = rk ) A(N), the R8CI matrix.
!
!    Output, real ( kind = rk ) B(N,N), the R8GE matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) b(n,n)
  integer i

  do i = 1, n
    b(i,1:i-1) = a(n+2-i:n)
    b(i,i:n) = a(1:n+1-i)
  end do

  return
end
subroutine r8ci_zeros ( n, a )

!*****************************************************************************80
!
!! R8CI_ZEROS zeroes an R8CI matrix.
!
!  Discussion:
!
!    The R8CI storage format is used for a real N by N circulant matrix.
!    An N by N circulant matrix A has the property that the entries on
!    row I appear again on row I+1, shifted one position to the right,
!    with the final entry of row I appearing as the first of row I+1.
!    The R8CI format simply records the first row of the matrix.
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
!    Output, real ( kind = rk ) A(N), the R8CI matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)

  a(1:n) = 0.0D+00

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

