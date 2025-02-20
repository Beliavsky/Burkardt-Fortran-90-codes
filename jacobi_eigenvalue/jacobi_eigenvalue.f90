subroutine jacobi_eigenvalue ( n, a, it_max, v, d, it_num, rot_num )

!*****************************************************************************80
!
!! JACOBI_EIGENVALUE carries out the Jacobi eigenvalue iteration.
!
!  Discussion:
!
!    This function computes the eigenvalues and eigenvectors of a
!    real symmetric matrix, using Rutishauser's modfications of the classical
!    Jacobi rotation method with threshold pivoting. 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 September 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Gene Golub, Charles VanLoan,
!    Matrix Computations,
!    Third Edition,
!    Johns Hopkins, 1996,
!    ISBN: 0-8018-4513-X,
!    LC: QA188.G65.
!
!  Input:
!
!    integer N, the order of the matrix.
!
!    real ( kind = rk ) A(N,N), the matrix, which must be square, real,
!    and symmetric.
!
!    integer IT_MAX, the maximum number of iterations.
!
!  Output:
!
!    real ( kind = rk ) V(N,N), the matrix of eigenvectors.
!
!    real ( kind = rk ) D(N), the eigenvalues, in descending order.
!
!    integer IT_NUM, the total number of iterations.
!
!    integer ROT_NUM, the total number of rotations.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) bw(n)
  real ( kind = rk ) c
  real ( kind = rk ) d(n)
  real ( kind = rk ) g
  real ( kind = rk ) gapq
  real ( kind = rk ) h
  integer i
  integer it_max
  integer it_num
  integer j
  integer k
  integer l
  integer m
  integer p
  integer q
  integer rot_num
  real ( kind = rk ) s
  real ( kind = rk ) t
  real ( kind = rk ) tau
  real ( kind = rk ) term
  real ( kind = rk ) termp
  real ( kind = rk ) termq
  real ( kind = rk ) theta
  real ( kind = rk ) thresh
  real ( kind = rk ) v(n,n)
  real ( kind = rk ) w(n)
  real ( kind = rk ) zw(n)

  do j = 1, n
    do i = 1, n
      v(i,j) = 0.0D+00
    end do
    v(j,j) = 1.0D+00
  end do

  do i = 1, n
    d(i) = a(i,i)
  end do

  bw(1:n) = d(1:n)
  zw(1:n) = 0.0D+00
  it_num = 0
  rot_num = 0

  do while ( it_num < it_max )

    it_num = it_num + 1
!
!  The convergence threshold is based on the size of the elements in
!  the strict upper triangle of the matrix.
!
    thresh = 0.0D+00
    do j = 1, n
      do i = 1, j - 1
        thresh = thresh + a(i,j) ** 2
      end do
    end do

    thresh = sqrt ( thresh ) / real ( 4 * n, kind = rk )

    if ( thresh == 0.0D+00 ) then
      exit 
    end if

    do p = 1, n
      do q = p + 1, n

        gapq = 10.0D+00 * abs ( a(p,q) )
        termp = gapq + abs ( d(p) )
        termq = gapq + abs ( d(q) )
!
!  Annihilate tiny offdiagonal elements.
!
        if ( 4 < it_num .and. &
             termp == abs ( d(p) ) .and. &
             termq == abs ( d(q) ) ) then

          a(p,q) = 0.0D+00
!
!  Otherwise, apply a rotation.
!
        else if ( thresh <= abs ( a(p,q) ) ) then

          h = d(q) - d(p)
          term = abs ( h ) + gapq

          if ( term == abs ( h ) ) then
            t = a(p,q) / h
          else
            theta = 0.5D+00 * h / a(p,q)
            t = 1.0D+00 / ( abs ( theta ) + sqrt ( 1.0D+00 + theta * theta ) )
            if ( theta < 0.0D+00 ) then 
              t = - t
            end if
          end if

          c = 1.0D+00 / sqrt ( 1.0D+00 + t * t )
          s = t * c
          tau = s / ( 1.0D+00 + c )
          h = t * a(p,q)
!
!  Accumulate corrections to diagonal elements.
!
          zw(p) = zw(p) - h                  
          zw(q) = zw(q) + h
          d(p) = d(p) - h
          d(q) = d(q) + h

          a(p,q) = 0.0D+00
!
!  Rotate, using information from the upper triangle of A only.
!
          do j = 1, p - 1
            g = a(j,p)
            h = a(j,q)
            a(j,p) = g - s * ( h + g * tau )
            a(j,q) = h + s * ( g - h * tau )
          end do

          do j = p + 1, q - 1
            g = a(p,j)
            h = a(j,q)
            a(p,j) = g - s * ( h + g * tau )
            a(j,q) = h + s * ( g - h * tau )
          end do

          do j = q + 1, n
            g = a(p,j)
            h = a(q,j)
            a(p,j) = g - s * ( h + g * tau )
            a(q,j) = h + s * ( g - h * tau )
          end do
!
!  Accumulate information in the eigenvector matrix.
!
          do j = 1, n
            g = v(j,p)
            h = v(j,q)
            v(j,p) = g - s * ( h + g * tau )
            v(j,q) = h + s * ( g - h * tau )
          end do

          rot_num = rot_num + 1

        end if

      end do
    end do

    bw(1:n) = bw(1:n) + zw(1:n)
    d(1:n) = bw(1:n)
    zw(1:n) = 0.0D+00

  end do
!
!  Restore upper triangle of input matrix.
!
  do j = 1, n
    do i = 1, j - 1
      a(i,j) = a(j,i)
    end do
  end do
!
!  Ascending sort the eigenvalues and eigenvectors.
!
  do k = 1, n - 1

    m = k

    do l = k + 1, n
      if ( d(l) < d(m) ) then
        m = l
      end if
    end do

    if ( m /= k ) then

      t    = d(m)
      d(m) = d(k)
      d(k) = t

      w(1:n)   = v(1:n,m)
      v(1:n,m) = v(1:n,k)
      v(1:n,k) = w(1:n)

    end if

  end do

  return
end
subroutine r8mat_diag_get_vector ( n, a, v )

!*****************************************************************************80
!
!! R8MAT_DIAG_GET_VECTOR gets the value of the diagonal of an R8MAT.
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
!    22 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of rows and columns of
!    the matrix.
!
!    real ( kind = rk ) A(N,N), the N by N matrix.
!
!  Output:
!
!    real ( kind = rk ) V(N), the diagonal entries
!    of the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,n)
  integer i
  real ( kind = rk ) v(n)

  do i = 1, n
    v(i) = a(i,i)
  end do

  return
end
subroutine r8mat_identity ( n, a )

!*****************************************************************************80
!
!! R8MAT_IDENTITY stores the identity matrix in an R8MAT.
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
!    24 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the order of A.
!
!  Output:
!
!    real ( kind = rk ) A(N,N), the N by N identity matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,n)
  integer i

  a(1:n,1:n) = 0.0D+00

  do i = 1, n
    a(i,i) = 1.0D+00
  end do

  return
end
subroutine r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )

!*****************************************************************************80
!
!! R8MAT_IS_EIGEN_RIGHT determines the error in a (right) eigensystem.
!
!  Discussion:
!
!    An R8MAT is a matrix of real ( kind = rk ) values.
!
!    This routine computes the Frobenius norm of
!
!      A * X - X * LAMBDA
!
!    where
!
!      A is an N by N matrix,
!      X is an N by K matrix (each of K columns is an eigenvector)
!      LAMBDA is a K by K diagonal matrix of eigenvalues.
!
!    This routine assumes that A, X and LAMBDA are all real!
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the order of the matrix.
!
!    integer K, the number of eigenvectors.
!    K is usually 1 or N.
!
!    real ( kind = rk ) A(N,N), the matrix.
!
!    real ( kind = rk ) X(N,K), the K eigenvectors.
!
!    real ( kind = rk ) LAMBDA(K), the K eigenvalues.
!
!  Output:
!
!    real ( kind = rk ) ERROR_FROBENIUS, the Frobenius norm
!    of the difference matrix A * X - X * LAMBDA, which would be exactly zero
!    if X and LAMBDA were exact eigenvectors and eigenvalues of A.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer k
  integer n

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) c(n,k)
  real ( kind = rk ) error_frobenius
  integer j
  real ( kind = rk ) lambda(k)
  real ( kind = rk ) r8mat_norm_fro
  real ( kind = rk ) x(n,k)

  c(1:n,1:k) = matmul ( a(1:n,1:n), x(1:n,1:k) )

  do j = 1, k
    c(1:n,j) = c(1:n,j) - lambda(j) * x(1:n,j)
  end do

  error_frobenius = r8mat_norm_fro ( n, k, c )

  return
end
function r8mat_norm_fro ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_NORM_FRO returns the Frobenius norm of an M by N R8MAT.
!
!  Discussion:
!
!    An R8MAT is a matrix of real ( kind = rk ) values.
!
!    The Frobenius norm is defined as
!
!      R8MAT_NORM_FRO = sqrt (
!        sum ( 1 <= I <= M ) Sum ( 1 <= J <= N ) A(I,J)^2 )
!
!    The matrix Frobenius-norm is not derived from a vector norm, but
!    is compatible with the vector L2 norm, so that:
!
!      r8vec_norm_l2 ( A*x ) <= r8mat_norm_fro ( A ) * r8vec_norm_l2 ( x ).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, N, the order of the matrix.
!
!    real ( kind = rk ) A(M,N), the matrix.
!
!  Output:
!
!    real ( kind = rk ) R8MAT_NORM_FRO, the Frobenius norm of A.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) r8mat_norm_fro

  r8mat_norm_fro = sqrt ( sum ( a(1:m,1:n)**2 ) )

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
!  Input:
!
!    integer N, the number of components of the vector.
!
!    real ( kind = rk ) A(N), the vector to be printed.
!
!    character ( len = * ) TITLE, a title.
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
