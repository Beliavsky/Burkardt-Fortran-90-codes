subroutine bivand1 ( n, alpha, beta, a )

!*****************************************************************************80
!
!! bivand1() returns a bidimensional Vandermonde1 matrix.
!
!  Discussion:
!
!    N = 3, ALPHA = ( 1, 2, 3 ), BETA = ( 10, 20, 30 )
!
!    (x,y)   | (1,10)  (2,10)  (3,10)  (1,20)  (2,20)  (1,30)
!    --------+-----------------------------------------------
!    1       |     1       1       1       1       1       1  
!    x       |     1       2       3       1       2       1
!       y    |    10      10      10      20      20      30
!    x^2     |     1       4       9       1       4       1
!    x  y    |    10      20      30      20      40      30
!    x^2y^2  |   100     100     100     400     400     900
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the data vectors.
!
!    Input, real ( kind = rk ) ALPHA(N), BETA(N), the values that define A.
!
!    Output, real ( kind = rk ) A(((N+1)*N)/2,((N+1)*N)/2), the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(((n+1)*n)/2,((n+1)*n)/2)
  real ( kind = rk ) alpha(n)
  real ( kind = rk ) beta(n)
  integer e
  integer e1
  integer e2
  integer ii
  integer j1
  integer j2
  integer jj
  integer n2

  n2 = ( n * ( n + 1 ) ) / 2

  e1 = 0
  e2 = 0
  e = 0

  do ii = 1, n2

    j1 = 1
    j2 = 1

    do jj = 1, n2

      if ( ii == 1 ) then
        a(ii,jj) = 1.0D+00
      else
        a(ii,jj) = alpha(j1) ** e1 * beta(j2) ** e2
      end if

      if ( j1 + j2 < n + 1 ) then
        j1 = j1 + 1
      else
        j1 = 1
        j2 = j2 + 1
      end if

    end do

    if ( e2 < e ) then
      e1 = e1 - 1
      e2 = e2 + 1
    else
      e = e + 1
      e1 = e
      e2 = 0
    end if

  end do

  return
end
subroutine bivand2 ( n, alpha, beta, a )

!*****************************************************************************80
!
!! BIVAND2 returns a bidimensional Vandermonde1 matrix.
!
!  Discussion:
!
!    N = 3, ALPHA = ( 1, 2, 3 ), BETA = ( 10, 20, 30 )
!
!    (x,y)   | (1,10) (2,10) (3,10) (1,20) (2,20) (3,20) (1,30) (2,30) (3,30)
!    --------+---------------------------------------------------------------
!    1       |     1      1      1      1      1      1      1      1      1  
!    x       |     1      2      3      1      2      1      1      2      3
!    x^2     |     1      4      9      1      4      1      1      4      9
!       y    |    10     10     10     20     20     20     30     30     30
!    x  y    |    10     20     30     20     40     60     30     60     90
!    x^2y    |    10     40     90     20     80    180     30    120    270
!       y^2  |   100    100    100    400    400    400    900    900    900
!    x  y^2  |   100    200    300    400    800   1200    900   1800   2700
!    x^2y^2  |   100    400    900    400   1600   3600    900   3600   8100
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the data vectors.
!
!    Input, real ( kind = rk ) ALPHA(N), BETA(N), the values that define A.
!
!    Output, real ( kind = rk ) A(N^2,N^2), the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n*n,n*n)
  real ( kind = rk ) alpha(n)
  real ( kind = rk ) beta(n)
  integer i
  integer ix
  integer iy
  integer j
  integer jx
  integer jy

  i = 0
  do iy = 1, n
    do ix = 1, n
      i = i + 1
      j = 0
      do jy = 1, n
        do jx = 1, n
          j = j + 1
          a(i,j) = alpha(jx) ** ( ix - 1 ) * beta(jy) ** ( iy - 1 )
        end do
      end do
    end do
  end do

  return
end
subroutine dvand ( n, alpha, b, x )

!*****************************************************************************80
!
!! DVAND solves a Vandermonde system A' * x = b.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Ake Bjorck, Victor Pereyra,
!    Solution of Vandermonde Systems of Equations,
!    Mathematics of Computation,
!    Volume 24, Number 112, October 1970, pages 893-903.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real ( kind = rk ) ALPHA(N), the parameters that define the matrix.
!    The values should be distinct.
!
!    Input, real ( kind = rk ) B(N), the right hand side of the linear system.
!
!    Output, real ( kind = rk ) X(N), the solution of the linear system.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) alpha(n)
  real ( kind = rk ) b(n)
  integer j
  integer k
  real ( kind = rk ) x(n)

  x(1:n) = b(1:n)

  do k = 1, n - 1
    do j = n, k + 1, -1
      x(j) = ( x(j) - x(j-1) ) / ( alpha(j) - alpha(j-k) )
    end do
  end do

  do k = n - 1, 1, -1
    do j = k, n - 1
      x(j) = x(j) - alpha(k) * x(j+1)
    end do
  end do

  return
end
subroutine dvandprg ( n, alpha, b, x, c, m )

!*****************************************************************************80
!
!! DVANDPRG solves a Vandermonde system A' * x = f progressively.
!
!  Discussion:
!
!    This function receives the solution to the system of equations A' * x = f
!    where A is a Vandermonde matrix for alpha(0) through alpha(n-1),
!    and new values alpha(n) and f(n).  It updates the solution.
!
!    To solve a system of Nbig equations, this function may be called 
!    repeatedly, with N = 1, 2, ..., Nbig.  Each time, a solution to the 
!    current subsystem is returned.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 April 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Ake Bjorck, Victor Pereyra,
!    Solution of Vandermonde Systems of Equations,
!    Mathematics of Computation,
!    Volume 24, Number 112, October 1970, pages 893-903.
!
!  Parameters:
!
!    Input, integer N, the new order of the matrix, which is 1 
!    larger than on the previous call.  For the first call, N must be 1.
!
!    Input, real ( kind = rk ) ALPHA(N), the parameters that define the matrix.
!    The values should be distinct.  The value ALPHA(N) has just been
!    added to the system.
!
!    Input, real ( kind = rk ) B(N), the right hand side of the linear system.
!
!    Input/output, real ( kind = rk ) X(N).  On input, the first N-1 entries 
!    contain the solution of the N-1xN-1 linear system.  On output, the 
!    solution to the NxN linear system.
!
!    Input/output, real ( kind = rk ) C(N), M(N).  On input, the first N-1 
!    entries contain factorization data for the N-1xN-1 linear system.  On 
!    output, factorization data for the NxN linear system.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) alpha(n)
  real ( kind = rk ) b(n)
  real ( kind = rk ) c(n)
  real ( kind = rk ) cn
  integer j
  real ( kind = rk ) m(n)
  real ( kind = rk ) x(n)

  c(n) = b(n)
  do j = n - 1, 1, -1
    c(j) = ( c(j+1) - c(j) ) / ( alpha(n) - alpha(j) )
  end do

  if ( n == 1 ) then
    m(n) = 1.0D+00
  else
    m(n) = 0.0D+00
  end if

  cn = c(1)
  x(n) = c(1)

  do j = n - 1, 1, -1
    m(j+1) = m(j+1) - alpha(n-1) * m(j)
    x(n-j) = x(n-j) + m(j+1) * cn
  end do

  return
end
subroutine pvand ( n, alpha, b, x )

!*****************************************************************************80
!
!! PVAND solves a Vandermonde system A * x = b.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Ake Bjorck, Victor Pereyra,
!    Solution of Vandermonde Systems of Equations,
!    Mathematics of Computation,
!    Volume 24, Number 112, October 1970, pages 893-903.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real ( kind = rk ) ALPHA(N), the parameters that define the matrix.
!    The values should be distinct.
!
!    Input, real ( kind = rk ) B(N), the right hand side of the linear system.
!
!    Output, real ( kind = rk ) X(N), the solution of the linear system.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) alpha(n)
  real ( kind = rk ) b(n)
  integer j
  integer k
  real ( kind = rk ) x(n)

  x(1:n) = b(1:n)

  do k = 1, n - 1
    do j = n, k + 1, -1
      x(j) = x(j) - alpha(k) * x(j-1)
    end do
  end do

  do k = n - 1, 1, -1
    do j = k + 1, n
      x(j) = x(j) / ( alpha(j) - alpha(j-k) )
    end do
    do j = k, n - 1
      x(j) = x(j) - x(j+1)
    end do
  end do

  return
end
subroutine pvandprg ( n, alpha, b, x, d, u )

!*****************************************************************************80
!
!! PVANDPRG solves a Vandermonde system A * x = f progressively.
!
!  Discussion:
!
!    This function receives the solution to the system of equations A * x = f
!    where A is a Vandermonde matrix for alpha(0) through alpha(n-1),
!    and new values alpha(n) and f(n).  It updates the solution.
!
!    To solve a system of Nbig equations, this function may be called 
!    repeatedly, with N = 1, 2, ..., Nbig.  Each time, a solution to the 
!    current subsystem is returned.
!
!    Note that the reference, which lists an Algol version of this algorithm, 
!    omits a minus sign, writing
!      u[j] := u[j] x delta;
!    where
!      u[j] := - u[j] x delta;
!    is actually necessary.  
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 April 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Ake Bjorck, Victor Pereyra,
!    Solution of Vandermonde Systems of Equations,
!    Mathematics of Computation,
!    Volume 24, Number 112, October 1970, pages 893-903.
!
!  Parameters:
!
!    Input, integer N, the new order of the matrix, which is 1 
!    larger than on the previous call.  For the first call, N must be 1.
!
!    Input, real ( kind = rk ) ALPHA(N), the parameters that define the matrix.
!    The values should be distinct.  The value ALPHA(N) has just been
!    added to the system.
!
!    Input, real ( kind = rk ) B(N), the right hand side of the linear system.
!
!    Input/output, real ( kind = rk ) X(N); on input, the solution of the 
!    N-1xN-1 linear system.  On output, the solution of the NxN linear system.
!
!    Input/output, real ( kind = rk ) D(N), U(N); on input, factorization data 
!    for the N-1xN-1 linear system.  On output, factorization data for the
!    NxN linear system.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) alpha(n)
  real ( kind = rk ) b(n)
  real ( kind = rk ) d(n)
  real ( kind = rk ) delta
  real ( kind = rk ) dn
  integer j
  real ( kind = rk ) u(n)
  real ( kind = rk ) x(n)

  d(n) = b(n)
  do j = n - 1, 1, -1
    d(j) = d(j+1) - alpha(n-j) * d(j)
  end do

  dn = d(1)
  u(n) = 1.0D+00

  do j = 1, n - 1
    delta = alpha(n) - alpha(j)
    u(j) = - u(j) * delta
    u(n) = u(n) * delta
    x(j) = x(j) + dn / u(j)
  end do

  x(n) = dn / u(n)

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
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
!    12 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows in A.
!
!    Input, integer N, the number of columns in A.
!
!    Input, real ( kind = rk ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  character ( len = * ) title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_PRINT_SOME prints some of an R8MAT.
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
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = rk ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ILO, JLO, the first row and column to print.
!
!    Input, integer IHI, JHI, the last row and column to print.
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

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)'
    return
  end if

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == real ( int ( a(i,j) ), kind = rk ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

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
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 8 ) ampm
  integer d
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  integer values(8)
  integer y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine vand1 ( n, x, a )

!*****************************************************************************80
!
!! VAND1 returns the Vandermonde1 matrix A with 1's on the first row.
!
!  Formula:
!
!    A(I,J) = X(J)^(I-1)
!
!  Example:
!
!    N = 5, X = ( 2, 3, 4, 5, 6 )
!
!    1  1   1   1   1
!    2  3   4   5   6
!    4  9  16  25  36
!    8 27  64 125  216
!   16 81 256 625 1296
!
!  Properties:
!
!    A is generally not symmetric: A' /= A.
!
!    A is nonsingular if, and only if, the X values are distinct.
!
!    det ( A ) = product ( 1 <= I <= N ) ( 1 <= J < I ) ( X(I) - X(J) ).
!             = product ( 1 <= J <= N ) X(J)
!             * product ( 1 <= I < J ) ( X(J) - X(I) ).
!
!    A is generally ill-conditioned.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Gregory, David Karney,
!    A Collection of Matrices for Testing Computational Algorithms,
!    Wiley, 1969, page 27,
!    LC: QA263.G68.
!
!    Nicholas Higham,
!    Stability analysis of algorithms for solving confluent
!    Vandermonde-like systems,
!    SIAM Journal on Matrix Analysis and Applications,
!    Volume 11, 1990, pages 23-41.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix desired.
!
!    Input, real ( kind = rk ) X(N), the values that define A.
!
!    Output, real ( kind = rk ) A(N,N), the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,n)
  integer i
  integer j
  real ( kind = rk ) x(n)

  do i = 1, n

    do j = 1, n

      if ( i == 1 .and. x(j) == 0.0D+00 ) then
        a(i,j) = 1.0D+00
      else
        a(i,j) = x(j) ** ( i - 1 )
      end if

    end do
  end do

  return
end
