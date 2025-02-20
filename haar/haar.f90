subroutine haar_1d ( n, x )

!*****************************************************************************80
!
!! HAAR_1D computes the Haar transform of a vector.
!
!  Discussion:
!
!    For the classical Haar transform, N should be a power of 2.
!    However, this is not required here.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 March 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of the vector.
!
!    Input/output, real ( kind = rk ) X(N), on input, the vector to be 
!    transformed.  On output, the transformed vector.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer k
  real ( kind = rk ) s
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)

  s = sqrt ( 2.0D+00 )
!
!  Initialize
!
  y(1:n) = 0.0D+00
!
!  Determine K, the largest power of 2 such that K <= N.
!
  k = 1
  do while ( k * 2 <= n )
    k = k * 2
  end do

  do while ( 1 < k )
  
    k = k / 2

    do i = 1, k
      y(i)   = ( x(2*i-1) + x(2*i) ) / s
      y(i+k) = ( x(2*i-1) - x(2*i) ) / s
    end do

    x(1:2*k) = y(1:2*k)

  end do

  return
end
subroutine haar_1d_inverse ( n, x )

!*****************************************************************************80
!
!! HAAR_1D_INVERSE computes the inverse Haar transform of a vector.
!
!  Discussion:
!
!    For the classical Haar transform, N should be a power of 2.
!    However, this is not required here.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 March 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of the vector.
!
!    Input/output, real ( kind = rk ) X(N), on input, the vector to be
!    transformed.  On output, the transformed vector.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer k
  real ( kind = rk ) s
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)

  s = sqrt ( 2.0D+00 )
!
!  Initialize.
!
  y(1:n) = 0.0D+00

  k = 1

  do while ( k * 2 <= n )

    do i = 1, k
      y(2*i-1) = ( x(i) + x(i+k) ) / s
      y(2*i)   = ( x(i) - x(i+k) ) / s
    end do

    x(1:2*k) = y(1:2*k)
    k = k * 2

  end do

  return
end
subroutine haar_2d ( m, n, u )

!*****************************************************************************80
!
!! HAAR_2D computes the Haar transform of an array.
!
!  Discussion:
!
!    For the classical Haar transform, M and N should be a power of 2.
!    However, this is not required here.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 March 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the dimensions of the array.
!
!    Input/output, real ( kind = rk ) U(M,N), the array to be transformed.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer k
  real ( kind = rk ) s
  real ( kind = rk ) u(m,n)
  real ( kind = rk ) v(m,n)

  s = sqrt ( 2.0D+00 )

  v(1:m,1:n) = u(1:m,1:n)
!
!  Determine K, the largest power of 2 such that K <= M.
!
  k = 1
  do while ( k * 2 <= m )
    k = k * 2
  end do
!
!  Transform all columns.
!
  do while ( 1 < k )
  
    k = k / 2

    v(  1:  k,1:n) = ( u(1:2*k-1:2,1:n) + u(2:2*k:2,1:n) ) / s
    v(k+1:k+k,1:n) = ( u(1:2*k-1:2,1:n) - u(2:2*k:2,1:n) ) / s

    u(1:2*k,1:n) = v(1:2*k,1:n)

  end do
!
!  Determine K, the largest power of 2 such that K <= N.
!
  k = 1
  do while ( k * 2 <= n )
    k = k * 2
  end do
!
!  Transform all rows.
!
  do while ( 1 < k )
  
    k = k / 2

    v(1:m,  1:  k) = ( u(1:m,1:2*k-1:2) + u(1:m,2:2*k:2) ) / s
    v(1:m,k+1:k+k) = ( u(1:m,1:2*k-1:2) - u(1:m,2:2*k:2) ) / s

    u(1:m,1:2*k) = v(1:m,1:2*k)

  end do

  return
end
subroutine haar_2d_inverse ( m, n, u )

!*****************************************************************************80
!
!! HAAR_2D_INVERSE inverts the Haar transform of an array.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 March 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the dimensions of the array.
!
!    Input/output, real ( kind = rk ) U(M,N), the array to be transformed.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer k
  real ( kind = rk ) s
  real ( kind = rk ) u(m,n)
  real ( kind = rk ) v(m,n)

  s = sqrt ( 2.0D+00 )

  v(1:m,1:n) = u(1:m,1:n)
!
!  Inverse transform of all rows.
!
  k = 1

  do while ( k * 2 <= n )

    v(1:m,1:2*k-1:2) = ( u(1:m,1:k) + u(1:m,1+k:k+k) ) / s
    v(1:m,2:2*k:2)   = ( u(1:m,1:k) - u(1:m,1+k:k+k) ) / s

    u(1:m,1:2*k) = v(1:m,1:2*k)
    k = k * 2

  end do
!
!  Inverse transform of all columns.
!
  k = 1

  do while ( k * 2 <= m )

    v(1:2*k-1:2,1:n) = ( u(1:k,1:n) + u(1+k:k+k,1:n) ) / s
    v(2:2*k:2,1:n)   = ( u(1:k,1:n) - u(1+k:k+k,1:n) ) / s

    u(1:2*k,1:n) = v(1:2*k,1:n)
    k = k * 2

  end do

  return
end
function r8mat_diff_frobenius ( m, n, a1, a2 )

!*****************************************************************************80
!
!! R8MAT_DIFF_FROBENIUS returns the Frobenius norm of an R8MAT difference.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    The Frobenius norm is defined as
!
!      R8MAT_DIFF_FROBENIUS = sqrt (
!        sum ( 1 <= I <= M ) sum ( 1 <= j <= N ) A(I,J) * A(I,J) )
!
!    The matrix Frobenius norm is not derived from a vector norm, but
!    is compatible with the vector L2 norm, so that:
!
!      r8vec_norm_l2 ( A * x ) <= 
!        r8mat_diff_frobenius ( A ) * r8vec_norm_l2 ( x ).
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
!  Parameters:
!
!    Input, integer M, the number of rows.
!
!    Input, integer N, the number of columns.
!
!    Input, real ( kind = rk ) A1(M,N), A2(M,N), the matrices for whose 
!    difference the Frobenius norm is desired.
!
!    Output, real ( kind = rk ) R8MAT_DIFF_FROBENIUSE, the Frobenius 
!    norm of A1 - A2.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a1(m,n)
  real ( kind = rk ) a2(m,n)
  real ( kind = rk ) r8mat_diff_frobenius

  r8mat_diff_frobenius = sqrt ( sum ( ( a1(1:m,1:n) - a2(1:m,1:n) )**2 ) )

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
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
!    An R8MAT is an array of R8 values.
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
function r8vec_diff_norm ( n, a, b )

!*****************************************************************************80
!
!! R8VEC_DIFF_NORM returns the L2 norm of the difference of R8VEC's.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The vector L2 norm is defined as:
!
!      R8VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 April 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in A.
!
!    Input, real ( kind = rk ) A(N), B(N), the vectors
!
!    Output, real ( kind = rk ) R8VEC_DIFF_NORM, the L2 norm of A - B.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) b(n)
  real ( kind = rk ) r8vec_diff_norm

  r8vec_diff_norm = sqrt ( sum ( ( a(1:n) - b(1:n) )**2 ) )

  return
end
subroutine r8vec_linspace ( n, a_first, a_last, a )

!*****************************************************************************80
!
!! R8VEC_LINSPACE returns a vector of linearly spaced values.
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
!    14 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real ( kind = rk ) A_FIRST, A_LAST, the first and last 
!    entries of A.
!
!    Output, real ( kind = rk ) A(N), the vector.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) a_first
  real ( kind = rk ) a_last
  integer i

  if ( n == 1 ) then
    a(1) = ( a_first + a_last ) / 2.0D+00
  else
  
    do i = 1, n
      a(i) = ( real ( n - i,     kind = rk ) * a_first &
             + real (     i - 1, kind = rk ) * a_last ) &
             / real ( n     - 1, kind = rk )
    end do

  end if

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

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
