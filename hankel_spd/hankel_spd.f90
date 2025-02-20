subroutine hankel_spd_cholesky_lower ( n, lii, liim1, l )

!*****************************************************************************80
!
!! hankel_spd_cholesky_lower() returns L such that L*L' is Hankel SPD.
!
!  Discussion:
!
!    In other words, H = L * L' is a symmetric positive definite matrix
!    with the property that H is constant along antidiagonals, so that
!
!      H(I+J) = h(k-1), for 1 <= I, J <= N, 1 <= K <= 2*N-1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 January 2017
!
!  Author:
!
!    S Al-Homidan, M Alshahrani.
!    FORTRAN90 implementation by John Burkardt.
!
!  Reference:
!
!    S Al-Homidan, M Alshahrani,
!    Positive Definite Hankel Matrices Using Cholesky Factorization,
!    Computational Methods in Applied Mathematics,
!    Volume 9, Number 3, pages 221-225, 2009.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real ( kind = rk ) LII(N), values to be used in L(I,I), 
!    for 1 <= I <= N.
!
!    Input, real ( kind = rk ) LIIM1(N-1), values to be used in L(I+1,I) 
!    for 1 <= I <= N-1.
!
!    Output, real ( kind = rk ) L(N,N), the lower Cholesky factor.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) alpha
  real ( kind = rk ) beta
  integer i
  integer j
  real ( kind = rk ) l(n,n)
  real ( kind = rk ) lii(n)
  real ( kind = rk ) liim1(n-1)
  integer q
  integer r
  integer s
  integer t

  l(1:n,1:n) = 0.0D+00

  do i = 1, n
    l(i,i) = lii(i)
  end do

  do i = 1, n - 1
    l(i+1,i) = liim1(i)
  end do

  do i = 3, n
    do j = 1, i - 2

      if ( mod ( i + j, 2 ) == 0 ) then
        q = ( i + j ) / 2
        r = q
      else
        q = ( i + j - 1 ) / 2
        r = q + 1
      end if

      alpha = 0.0D+00
      do s = 1, q
        alpha = alpha + l(q,s) * l(r,s)
      end do

      beta = 0.0D+00
      do t = 1, j - 1
        beta = beta + l(i,t) * l(j,t)
      end do

      l(i,j) = ( alpha - beta ) / l(j,j)

    end do
  end do

  return
end
subroutine r8mat_cholesky_factor ( n, a, c, flag )

!*****************************************************************************80
!
!! R8MAT_CHOLESKY_FACTOR computes the Cholesky factor of a symmetric matrix.
!
!  Discussion:
!
!    The matrix must be symmetric and positive semidefinite.
!
!    For a positive semidefinite symmetric matrix A, the Cholesky factorization
!    is a lower triangular matrix L such that:
!
!      A = L * L'
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of rows and columns of
!    the matrix A.
!
!    Input, real ( kind = rk ) A(N,N), the N by N matrix.
!
!    Output, real ( kind = rk ) C(N,N), the N by N lower triangular
!    Cholesky factor.
!
!    Output, integer FLAG:
!    0, no error occurred.
!    1, the matrix is not positive definite.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) c(n,n)
  integer flag
  integer i
  integer j
  real ( kind = rk ) sum2

  flag = 0

  c(1:n,1:n) = a(1:n,1:n)

  do j = 1, n

    c(1:j-1,j) = 0.0D+00

    do i = j, n

      sum2 = c(j,i) - dot_product ( c(j,1:j-1), c(i,1:j-1) )

      if ( i == j ) then
        if ( sum2 <= 0.0D+00 ) then
          flag = 1
          return
        else
          c(i,j) = sqrt ( sum2 )
        end if
      else
        if ( c(j,j) /= 0.0D+00 ) then
          c(i,j) = sum2 / c(j,j)
        else
          c(i,j) = 0.0D+00
        end if
      end if

    end do

  end do

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
