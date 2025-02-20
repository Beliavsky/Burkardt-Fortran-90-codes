program main

!*****************************************************************************80
!
!! eigs_test() tests eigs().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 June 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ck8 = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ), allocatable :: A(:,:)
  complex ( kind = ck8 ), allocatable :: lambda(:)
  integer n

  call timestamp ( )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'eigs_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test eigs():'

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  The Clement2 matrix is symmetric, and has'
  write ( *, '(a)' ) '  integer eigenvalues.'
  n = 5
  allocate ( A(1:n,1:n) )
  call clement2_matrix ( n, A )
  call r8mat_print ( n, n, A, '  The matrix A:' )
  allocate ( lambda(1:n) )
  call eigs ( n, A, lambda )
  call c8vec_print ( n, lambda, '  The eigenvalues LAMBDA' )
  deallocate ( A )
  deallocate ( lambda )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  The Helmert matrix is nonsymmetric, orthogonal, and has'
  write ( *, '(a)' ) '  complex eigenvalues of norm 1.'
  n = 5
  allocate ( A(1:n,1:n) )
  call helmert_matrix ( n, A )
  call r8mat_print ( n, n, A, '  The matrix A:' )
  allocate ( lambda(1:n) )
  call eigs ( n, A, lambda )
  call c8vec_print ( n, lambda, '  The eigenvalues LAMBDA' )
  deallocate ( A )
  deallocate ( lambda )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'eigs_test():'
  write ( *, '(a)' ) '  Normal end of execution.'

  call timestamp ( )

  stop 0
end
subroutine clement2_matrix ( n, a )

!*****************************************************************************80
!
!! clement2_matrix() returns the Clement2 matrix.
!
!  Formula:
!
!    if ( J = I+1 )
!      A(I,J) = sqrt(I*(N-I))
!    else if ( I = J+1 )
!      A(I,J) = sqrt(J*(N-J))
!    else
!      A(I,J) = 0
!
!  Example:
!
!    N = 5
!
!       .    sqrt(4)    .       .       .
!    sqrt(4)    .    sqrt(6)    .       .
!       .    sqrt(6)    .    sqrt(6)    .
!       .       .    sqrt(6)    .    sqrt(4)
!       .       .       .    sqrt(4)    .
!
!  Properties:
!
!    A is tridiagonal.
!
!    Because A is tridiagonal, it has property A (bipartite).
!
!    A is symmetric: A' = A.
!
!    Because A is symmetric, it is normal.
!
!    Because A is normal, it is diagonalizable.
!
!    A is persymmetric: A(I,J) = A(N+1-J,N+1-I).
!
!    The diagonal of A is zero.
!
!    A is singular if N is odd.
!
!    About 64 percent of the entries of the inverse of A are zero.
!
!    The eigenvalues are plus and minus the numbers
!
!      N-1, N-3, N-5, ..., (1 or 0).
!
!    If N is even,
!
!      det ( A ) = (-1)**(N/2) * (N-1) * (N+1)**(N/2)
!
!    and if N is odd,
!
!      det ( A ) = 0
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    P A Clement,
!    A class of triple-diagonal matrices for test purposes,
!    SIAM Review,
!    Volume 1, 1959, pages 50-52.
!
!  Input:
!
!    integer N, the order of A.
!
!  Output:
!
!    real A(N,N), the matrix.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  real ( kind = rk8 ) a(n,n)
  integer i
  integer j

  do i = 1, n
    do j = 1, n

      if ( j == i + 1 ) then
        a(i,j) = sqrt ( real ( i * ( n - i ), kind = rk8 ) )
      else if ( i == j + 1 ) then
        a(i,j) = sqrt ( real ( j * ( n - j ), kind = rk8 ) )
      else
        a(i,j) = 0.0D+00
      end if

    end do
  end do

  return
end
subroutine helmert_matrix ( n, a )

!*****************************************************************************80
!
!! helmert_matrix() returns the HELMERT matrix.
!
!  Discussion:
!
!    A matrix is a (standard) Helmert matrix if it is orthogonal,
!    and the elements which are above the diagonal and below the
!    first row are zero.
!
!    If the elements of the first row are all strictly positive,
!    then the matrix is a strictly Helmertian matrix.
!
!    It is possible to require in addition that all elements below
!    the diagonal be strictly positive, but in the reference, this
!    condition is discarded as being cumbersome and not useful.
!
!    A Helmert matrix can be regarded as a change of basis matrix
!    between a pair of orthonormal coordinate bases.  The first row
!    gives the coordinates of the first new basis vector in the old
!    basis.  Each later row describes combinations of (an increasingly
!    extensive set of) old basis vectors that constitute the new
!    basis vectors.
!
!    Helmert matrices have important applications in statistics.
!
!  Formula:
!
!    If I = 1 then
!      A(I,J) = 1 / sqrt ( N )
!    else if J < I then
!      A(I,J) = 1 / sqrt ( I * ( I - 1 ) )
!    else if J = I then
!      A(I,J) = - sqrt (I-1) / sqrt ( I )
!    else
!      A(I,J) = 0
!
!  Example:
!
!    N = 5
!
!    0.4472    0.4472    0.4472    0.4472    0.4472
!    0.7071   -0.7071         0         0         0
!    0.4082    0.4082   -0.8165         0         0
!    0.2887    0.2887    0.2887   -0.8660         0
!    0.2236    0.2236    0.2236    0.2236   -0.8944
!
!  Properties:
!
!    A is generally not symmetric: A' /= A.
!
!    A is orthogonal: A' * A = A * A' = I.
!
!    Because A is orthogonal, it is normal: A' * A = A * A'.
!
!    det ( A ) = (-1)^(N+1)
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 June 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    HO Lancaster,
!    The Helmert Matrices,
!    American Mathematical Monthly,
!    Volume 72, 1965, pages 4-12.
!
!  Input:
!
!    integer N, the order of the matrix.
!
!  Output:
!
!    real A(N,N), the matrix.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  real ( kind = rk8 ) a(n,n)
  integer i
  integer j
!
!  A begins with the first row, diagonal, and lower triangle set to 1.
!
  do i = 1, n
    do j = 1, n

      if ( i == 1 ) then
        a(i,j) = 1.0D+00 / sqrt ( real ( n, kind = rk8 ) )
      else if ( j < i ) then
        a(i,j) = 1.0D+00 / sqrt ( real ( i * ( i - 1 ), kind = rk8 ) )
      else if ( i == j ) then
        a(i,j) = - sqrt ( real ( i - 1, kind = rk8 ) ) &
                 / sqrt ( real ( i,     kind = rk8 ) )
      else
        a(i,j) = 0.0D+00
      end if

    end do
  end do

  return
end
subroutine c8vec_print ( n, a, title )

!*****************************************************************************80
!
!! c8vec_print() prints a C8VEC, with an optional title.
!
!  Discussion:
!
!    A C8VEC is a vector of C8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of components of the vector.
!
!    complex A(N), the vector to be printed.
!
!    character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: ck8 = kind ( ( 1.0D+00, 1.0D+00 ) )

  integer n

  complex ( kind = ck8 ) a(n)
  integer i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,2x,2g14.6)' ) i, a(i)
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
!    real A(M,N), the matrix.
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
!    real A(M,N), an M by N matrix to be printed.
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
