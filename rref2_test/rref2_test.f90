program main

!*****************************************************************************80
!
!! rref2_test() tests rref2().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 December 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'rref2_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test rref(), which analyzes matrices using the'
  write ( *, '(a)' ) '  reduced row echelon form (RREF)'

  call is_rref_test ( )
  call rref_compute_test ( )
  call rref_columns_test ( )
  call rref_determinant_test ( )
  call rref_inverse_test ( )
  call rref_rank_test ( )
  call rref_solve_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'rref2_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop
end
subroutine is_rref_test ( )

!*****************************************************************************80
!
!! is_rref_test() tests is_rref().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 December 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: m = 4
  integer, parameter :: n = 5

  logical is_rref
!
!  Zero rows must come last.
!
  real ( kind = rk8 ), dimension ( m, n ), save :: A0 = reshape ( (/ &
    1, 0, 0, 0, &
    0, 0, 0, 0, &
    0, 1, 0, 0, &
    9, 0, 0, 1, &
    4, 8, 0, 0  &
  /), (/ m, n /) )
!
!  First nonzero must be to right of first nonzero in previous row.
!
  real ( kind = rk8 ), dimension ( m, n ), save :: A1 = reshape ( (/ &
    1, 0, 0, 0, &
    0, 0, 0, 0, &
    0, 0, 1, 0, &
    9, 1, 0, 1, &
    4, 0, 8, 0  &
  /), (/ m, n /) )
!
!  First nonzero must be a 1.
!
  real ( kind = rk8 ), dimension ( m, n ), save :: A2 = reshape ( (/ &
    1, 0, 0, 0, &
    0, 1, 0, 0, &
    0, 0, 3, 0, &
    9, 2, 0, 0, &
    4, 8, 0, 0  &
  /), (/ m, n /) )
!
!  First nonzero must only nonzero in its column.
!
  real ( kind = rk8 ), dimension ( m, n ), save :: A3 = reshape ( (/ &
    1, 0, 0, 0, &
    0, 1, 0, 0, &
    3, 0, 1, 0, &
    9, 2, 0, 0, &
    4, 8, 0, 0  &
  /), (/ m, n /) )
!
!  RREF example.
!
  real ( kind = rk8 ), dimension ( m, n ), save :: A4 = reshape ( (/ &
    1, 0, 0, 0, &
    0, 1, 0, 0, &
    3, 2, 0, 0, &
    0, 0, 1, 0, &
    4, 8, 0, 0  &
  /), (/ m, n /) )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'is_rref_test():'
  write ( *, '(a)' ) '  is_rref() reports if a matrix is in reduced row echelon format.'

  call r8mat_print ( m, n, A0, '  Matrix A0:' )
  write ( *, '(a,l1)' ) '  is_rref(A0) = ', is_rref ( m, n, A0 )

  call r8mat_print ( m, n, A1, '  Matrix A1:' )
  write ( *, '(a,l1)' ) '  is_rref(A1) = ', is_rref ( m, n, A1 )

  call r8mat_print ( m, n, A2, '  Matrix A2:' )
  write ( *, '(a,l1)' ) '  is_rref(A2) = ', is_rref ( m, n, A2 )

  call r8mat_print ( m, n, A3, '  Matrix A3:' )
  write ( *, '(a,l1)' ) '  is_rref(A3) = ', is_rref ( m, n, A3 )

  call r8mat_print ( m, n, A4, '  Matrix A4:' )
  write ( *, '(a,l1)' ) '  is_rref(A4) = ', is_rref ( m, n, A4 )

  return
end
subroutine rref_columns_test ( )

!*****************************************************************************80
!
!! rref_columns_test() tests rref_columns().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 December 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: m = 7
  integer, parameter :: n = 4

  real ( kind = rk8 ), dimension ( m, n ), save :: A = reshape ( (/ &
    1,  2,  3,  4,  5,  6,  7, &
    2,  4,  6,  8, 10, 12, 14, &
    3,  9,  0,  0,  6,  6,  2, &
    1,  3,  0,  2,  6,  3,  1  &
  /), (/ m, n /) )
  integer a_cols(n)
  real ( kind = rk8 ), allocatable :: A_COLUMNS(:,:)
  real ( kind = rk8 ) A_RREF(m,n)
  integer j
  integer j2
  integer n2

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'rref_columns_test():'
  write ( *, '(a)' ) '  rref_columns() uses the reduced row echelon form (RREF)'
  write ( *, '(a)' ) '  of a matrix to find the linearly independent columns.'

  call r8mat_print ( m, n, A, '  Matrix A:' )

  call rref_compute ( m, n, A, A_RREF, a_cols )

  n2 = 0
  do j = 1, n
    if ( a_cols(j) /= 0 ) then
      n2 = n2 + 1
    end if
  end do

  write ( *, '(a)' ) ''
  write ( *, '(a,i2)' ) '  Number of independent columns is ', n2

  call i4vec_print ( n2, a_cols, '  Independent column indices' )

  allocate ( A_COLUMNS(m,n2) )
  j2 = 0
  do j = 1, n
    if ( a_cols(j) /= 0 ) then
      j2 = j2 + 1
      A_COLUMNS(1:m,j2) = A(1:m,a_cols(j))
    end if
  end do

  call r8mat_print ( m, n2, A_COLUMNS, '  Independent columns:' )

  deallocate ( A_COLUMNS )

  return
end
subroutine rref_compute_test ( )

!*****************************************************************************80
!
!! rref_compute_test() tests rref_compute().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 December 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: m = 4
  integer, parameter :: n = 7

  real ( kind = rk8 ), dimension ( m, n ), save :: A = reshape ( (/ &
    1, -2,  3, -1, &
    3, -6,  9, -3, &
    0,  0,  0,  0, &
    2, -2,  0,  1, &
    6, -8,  6,  0, &
    3,  3,  6,  9, &
    1,  1,  2,  3  &
  /), (/ m, n /) )
  integer a_cols(n)
  real ( kind = rk8 ) A_RREF(m,n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'rref_compute_test():'
  write ( *, '(a)' ) '  rref_compute() is a user-written code to compute the'
  write ( *, '(a)' ) '  reduced row echelon form (RREF) of a matrix.'

  call r8mat_print ( m, n, A, '  Matrix A:' )

  call rref_compute ( m, n, A, A_RREF, a_cols )

  call r8mat_print ( m, n, A_RREF, '  rref_compute(A):' )

  call i4vec_print ( n, a_cols, '  Column indices' )

  return
end
subroutine rref_determinant_test ( )

!*****************************************************************************80
!
!! rref_determinant_test() tests rref_determinant().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 December 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: n = 4

  real ( kind = rk8 ), dimension ( n, n ), save :: A = reshape ( (/ &
    5,  7,  6,  5, &
    7, 10,  8,  7, &
    6,  8, 10,  9, &
    5,  7,  9, 10  &
  /), (/ n, n /) )
  real ( kind = rk8 ) a_det

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'rref_determinant_test():'
  write ( *, '(a)' ) '  rref_determinant() uses the reduced row echelon form '
  write ( *, '(a)' ) '  of a square matrix to compute the determinant.'

  call r8mat_print ( n, n, A, '  matrix A:' )

  call rref_determinant ( n, A, a_det )
  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Estimated determinant of A = ', a_det

  return
end
subroutine rref_inverse_test ( )

!*****************************************************************************80
!
!! rref_inverse_test() tests rref_inverse().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 December 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: n = 4

  real ( kind = rk8 ), dimension ( n, n ), save :: A = reshape ( (/ &
    5,  7,  6,  5, &
    7, 10,  8,  7, &
    6,  8, 10,  9, &
    5,  7,  9, 10  &
  /), (/ n, n /) )
  real ( kind = rk8 ) A_INV(n,n)
  real ( kind = rk8 ) P(n,n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'rref_inverse_test():'
  write ( *, '(a)' ) '  rref_inverse() uses the reduced row echelon form '
  write ( *, '(a)' ) '  of a square matrix to compute its inverse.'

  call r8mat_print ( n, n, A, '  matrix A:' )

  call rref_inverse ( n, A, A_INV )

  call r8mat_print ( n, n, A_INV, '  Estimated inverse A_inv:' )

  P = matmul ( A_INV, A )

  call r8mat_print ( n, n, P, '  Product A_inv * A:' )

  return
end
subroutine rref_rank_test ( )

!*****************************************************************************80
!
!! rref_rank_test() tests rref_rank().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 December 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: m = 4
  integer, parameter :: n = 7

  real ( kind = rk8 ), dimension ( m, n ), save :: A = reshape ( (/ &
    1, -2,  3, -1, &
    3, -6,  9, -3, &
    0,  0,  0,  0, &
    2, -2,  0,  1, &
    6, -8,  6,  0, &
    3,  3,  6,  9, &
    1,  1,  2,  3  &
  /), (/ m, n /) )
  integer a_rank

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'rref_rank_test():'
  write ( *, '(a)' ) '  rref_rank() uses the reduced row echelon form '
  write ( *, '(a)' ) '  of a matrix to estimate its rank.'

  call r8mat_print ( m, n, A, '  matrix A:' )

  call rref_rank ( m, n, A, a_rank )
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  A has rank ', a_rank

  return
end
subroutine rref_solve_test ( )

!*****************************************************************************80
!
!! rref_solve_test() tests rref_solve().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 December 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: m = 4
  integer, parameter :: n = 4
  integer, parameter :: nb = 1

  real ( kind = rk8 ), dimension ( m, n ), save :: A = reshape ( (/ &
    5,  7,  6,  5, &
    7, 10,  8,  7, &
    6,  8, 10,  9, &
    5,  7,  9, 10  &
  /), (/ m, n /) )
  real ( kind = rk8 ) b(m,nb)
  real ( kind = rk8 ) b2(m,nb)
  real ( kind = rk8 ), dimension ( n, nb ) :: x = reshape ( (/ &
    1, 2, 3, 4 /), (/ n, nb /) )
  real ( kind = rk8 ) x2(n,nb)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'rref_solve_test():'
  write ( *, '(a)' ) '  rref_solve() uses the reduced row echelon form '
  write ( *, '(a)' ) '  of a square matrix to solve a linear system.'

  call r8mat_print ( m, n, A, '  matrix A:' )

  b = matmul ( A, x )
  call r8mat_print ( m, nb, b,  '  Right hand side b:' )

  call rref_solve ( m, n, A, nb, b, x2 )

  call r8mat_print ( n, nb, x2,  '  Estimated solution:' )

  b2 = matmul ( A, x2 )
  call r8mat_print ( m, nb, b2,  '  Product A * x:' )

  return
end
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! i4vec_print() prints an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 May 2010
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of components of the vector.
!
!    integer A(N), the vector to be printed.
!
!    character ( len = * ) TITLE, a title.
!
  implicit none

  integer n

  integer a(n)
  integer i
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,a,2x,i12)' ) i, ':', a(i)
  end do

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! r8mat_print() prints a real matrix.
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
!  Input:
!
!    integer M, the number of rows in A.
!
!    integer N, the number of columns in A.
!
!    real ( kind = rk8 ) A(M,N), the matrix.
!
!    character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk8 ) a(m,n)
  character ( len = * ) title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! r8mat_print_some prints some of a real matrix.
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
!  Input:
!
!    integer M, N, the number of rows and columns.
!
!    real ( kind = rk8 ) A(M,N), an M by N matrix to be printed.
!
!    integer ILO, JLO, the first row and column to print.
!
!    integer IHI, JHI, the last row and column to print.
!
!    character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: incx = 5
  integer m
  integer n

  real ( kind = rk8 ) a(m,n)
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

        write ( ctemp(j2), '(g14.6)' ) a(i,j)

      end do

      write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! timestamp() prints the current YMDHMS date as a time stamp.
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
!    15 August 2021
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

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end

