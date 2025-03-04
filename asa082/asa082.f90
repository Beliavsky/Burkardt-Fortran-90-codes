subroutine detq ( a, n, d, ifault )

!*****************************************************************************80
!
!! detq() computes the determinant of an orthogonal matrix.
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
!    Original FORTRAN77 version by J C Gower.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    J C Gower,
!    Algorithm AS 82:
!    The determinant of an orthogonal matrix,
!    Applied Statistics,
!    Volume 24, Number 1, 1975, page 150-153.
!
!  Input:
!
!    real ( kind = rk ) A(N,N), the orthogonal matrix stored by rows or columns.
!
!    integer N, the order of the matrix.
!
!  Output:
!
!    real ( kind = rk ) D, the determinant of A.
!
!    integer IFAULT, 
!    0, no error occurred.
!    1, an error was detected.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) a2(n*n)
  real ( kind = rk ) d
  integer i
  integer ifault
  integer j
  integer k
  integer p
  integer q
  integer r
  integer s
  real ( kind = rk ) tol
  real ( kind = rk ) x
  real ( kind = rk ) y

  ifault = 0
  tol = 0.0001D+00

  if ( n <= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'detq - Fatal error!'
    write ( *, '(a)' ) '  n <= 0'
    ifault = 1
    return
  end if

  k = 0
  do i = 1, n
    do j = 1, n
      k = k + 1
      a2(k) = a(i,j)
    end do
  end do

  d = 1.0D+00
  r = 1

  do k = 2, n + 1

    q = r
    x = a2(r)
    y = sign ( 1.0D+00, x )
    d = d * y
    y = - 1.0D+00 / ( x + y )
    x = abs ( x ) - 1.0D+00

    if ( tol < abs ( x ) ) then

      if ( 0.0D+00 < x ) then
        write ( *, '(a)' ) ''
        write ( *, '(a)' ) 'detq - Fatal error!'
        write ( *, '(a)' ) '  x < 0.0'
        write ( *, '(a,g14.6)' ) '  x = ', x
        ifault = 1
        return
      end if

      if ( k == n + 1 ) then
        write ( *, '(a)' ) ''
        write ( *, '(a)' ) 'detq - Fatal error!'
        write ( *, '(a)' ) '  k == n + 1'
        ifault = 1
        return
      end if

      do i = k, n
        q = q + n
        x = a2(q) * y
        p = r
        s = q
        do j = k, n
          p = p + 1
          s = s + 1
          a2(s) = a2(s) + x * a2(p)
        end do

      end do

    end if

    r = r + n + 1

  end do

  return
end
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

