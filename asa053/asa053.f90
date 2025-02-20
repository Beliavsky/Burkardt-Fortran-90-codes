subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! r8mat_print() prints an R8MAT.
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
!    27 August 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, the number of rows in A.
!
!    integer N, the number of columns in A.
!
!    real ( kind = rk ) A(M,N), the matrix.
!
!    character ( len = * ) TITLE, a title.
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
!! r8mat_print_some() prints some of an R8MAT.
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
!    27 August 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, N, the number of rows and columns.
!
!    real ( kind = rk ) A(M,N), an M by N matrix to be printed.
!
!    integer ILO, JLO, the first row and column to print.
!
!    integer IHI, JHI, the last row and column to print.
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
subroutine r8pp_print ( n, a, title )

!*****************************************************************************80
!
!! r8pp_print() prints an R8PP matrix.
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
!    27 August 2021
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
!    real ( kind = rk ) A((N*(N+1))/2), the R8PP matrix.
!
!    character ( len = * ) TITLE, a title.
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
!! r8pp_print_some() prints some of an R8PP matrix.
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
!    27 August 2021
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
!    real ( kind = rk ) A((N*(N+1))/2), the R8PP matrix.
!
!    integer ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    character ( len = * ) TITLE, a title.
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
subroutine r8utp_print ( n, a, title )

!*****************************************************************************80
!
!! r8utp_print() prints an R8UTP matrix.
!
!  Discussion:
!
!    The R8UTP storage format is appropriate for an upper triangular
!    matrix.  Only the upper triangle of the matrix is stored,
!    by successive partial columns, in an array of length (N*(N+1))/2,
!    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 August 2021
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
!    real ( kind = rk ) A((N*(N+1))/2), the matrix.
!
!    character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a((n*(n+1))/2)
  character ( len = * ) title

  call r8utp_print_some ( n, a, 1, 1, n, n, title )

  return
end
subroutine r8utp_print_some ( n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! r8utp_print_some() prints some of an R8UTP matrix.
!
!  Discussion:
!
!    The R8UTP storage format is appropriate for an upper triangular
!    matrix.  Only the upper triangle of the matrix is stored,
!    by successive partial columns, in an array of length (N*(N+1))/2,
!    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 August 2021
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
!    real ( kind = rk ) A((N*(N+1))/2), the matrix.
!
!    integer ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    character ( len = * ) TITLE, a title.
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
          aij = 0.0D+00
        end if

        write ( ctemp(j2), '(g14.6)' ) aij

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine rnorm ( u1, u2 )

!*****************************************************************************80
!
!! rnorm() returns two independent standard random normal deviates.
!
!  Discussion:
!
!    This routine sets U1 and U2 to two independent standardized 
!    random normal deviates.   This is a version of the 
!    method given in Knuth.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 August 2021
!
!  Author:
!
!    Original FORTRAN77 version by William Smith, Ronald Hocking.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Donald Knuth,
!    The Art of Computer Programming,
!    Volume 2, Seminumerical Algorithms,
!    Third Edition,
!    Addison Wesley, 1997,
!    ISBN: 0201896842,
!    LC: QA76.6.K64.
!
!  Output:
!
!    real ( kind = rk ) U1, U2, two standard random normal deviates.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) s
  real ( kind = rk ) u1
  real ( kind = rk ) u2
  real ( kind = rk ) x
  real ( kind = rk ) y

  do

    call random_number ( harvest = x )
    call random_number ( harvest = y )
    x = 2.0D+00 * x - 1.0D+00
    y = 2.0D+00 * y - 1.0D+00
    s = x * x + y * y

    if ( s <= 1.0D+00 ) then
      s = sqrt ( - 2.0D+00 * log ( s ) / s )
      u1 = x * s
      u2 = y * s
      exit
    end if

  end do

  return
end
subroutine wshrt ( d, n, np, sa )

!*****************************************************************************80
!
!! wshrt() returns a random Wishart variate.
!
!  Discussion:
!
!    This routine is a Wishart variate generator.  
!
!    SA is an upper-triangular matrix of size NP * NP,
!    written in linear form, column ordered, whose elements have a 
!    Wishart(N, SIGMA) distribution.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 August 2021
!
!  Author:
!
!    Original FORTRAN77 version by William Smith, Ronald Hocking.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Smith, Ronald Hocking,
!    Algorithm AS 53, Wishart Variate Generator,
!    Applied Statistics,
!    Volume 21, Number 3, pages 341-345, 1972.
!
!  Input:
!
!    real ( kind = rk ) D(NP*(NP+1)/2), the upper triangular array that
!    represents the Cholesky factor of the correlation matrix SIGMA.
!    D is stored in column-major form.
!
!    integer N, the number of degrees of freedom.
!    1 <= N <= NP.
!
!    integer NP, the size of variables.
!
!  Output:
!
!    real ( kind = rk ) SA(NP*(NP+1)/2), a sample from the 
!    Wishart distribution.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer np

  real ( kind = rk ) c
  real ( kind = rk ) d((np*(np+1))/2)
  real ( kind = rk ) df
  integer i
  integer ii
  integer ip
  integer j
  integer k
  integer n
  integer nnp
  integer nq
  integer nr
  integer ns
  real ( kind = rk ) rn
  real ( kind = rk ) sa((np*(np+1))/2)
  real ( kind = rk ) sb((np*(np+1))/2)
  real ( kind = rk ) u1
  real ( kind = rk ) u2

  k = 1
  nnp = ( np * ( np + 1 ) ) / 2
!
!  Load SB with independent normal (0, 1) variates.
!
  do while ( k <= nnp )

    call rnorm ( u1, u2 )

    sb(k) = u1
    k = k + 1

    if ( k <= nnp ) then
      sb(k) = u2
      k = k + 1
    end if

  end do
!
!  Load diagonal elements with square root of chi-square variates.
!
  ns = 0

  do i = 1, np
!
!  The original text read "DF = N - I + 1".
!  This should read "DF = NP - I + 1".
!
    df = real ( np - i + 1, kind = rk )
    ns = ns + i
    u1 = 2.0D+00 / ( 9.0D+00 * df )
    u2 = 1.0D+00 - u1
    u1 = sqrt ( u1 )
!
!  Wilson-Hilferty formula for approximating chi-square variates:
!  The original code did not take the absolute value!
!
    sb(ns) = sqrt ( df * abs ( ( u2 + sb(ns) * u1 ) ** 3 ) )

  end do

  rn = real ( n, kind = rk )
  nr = 1

  do i = 1, np
    nr = nr + i - 1
    do j = i, np
      ip = nr
      nq = ( j * ( j - 1 ) ) / 2 + i - 1
      c = 0.0D+00
      do k = i, j
        ip = ip + k - 1
        nq = nq + 1
        c = c + sb(ip) * d(nq)
      end do
      sa(ip) = c
    end do
  end do

  do i = 1, np
    ii = np - i + 1
    nq = nnp - np
    do j = 1, i
      ip = ( ii * ( ii - 1 ) ) / 2
      c = 0.0D+00
      do k = i, np
        ip = ip + 1
        nq = nq + 1
        c = c + sa(ip) * sa(nq)
      end do
      sa(nq) = c / rn
      nq = nq - 2 * np + i + j - 1
    end do
  end do

  return
end
 
