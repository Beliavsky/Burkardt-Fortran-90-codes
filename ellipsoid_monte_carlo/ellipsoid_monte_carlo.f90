subroutine ellipsoid_sample ( m, n, a, v, r, x )

!*****************************************************************************80
!
!! ellipsoid_sample() samples points uniformly from an ellipsoid.
!
!  Discussion:
!
!    The points X in the ellipsoid are described by a M by M
!    positive definite symmetric matrix A, a "center" V, and 
!    a "radius" R, such that
!
!      (X-V)' * A * (X-V) <= R * R
!
!    The algorithm computes the Cholesky factorization of A:
!
!      A = U' * U.
!
!    A set of uniformly random points Y is generated, satisfying:
!
!      Y' * Y <= R * R.
!
!    The appropriate points in the ellipsoid are found by solving
!
!      U * X = Y
!      X = X + V
!
!    Thanks to Dr Karl-Heinz Keil for pointing out that the original
!    coding was actually correct only if A was replaced by its inverse.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 August 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Reuven Rubinstein,
!    Monte Carlo Optimization, Simulation, and Sensitivity
!    of Queueing Networks,
!    Krieger, 1992,
!    ISBN: 0894647644,
!    LC: QA298.R79.
!
!  Parameters:
!
!    Input, integer M, the dimension of the space.
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) A(M,M), the matrix that describes
!    the ellipsoid.
!
!    Input, real ( kind = rk ) V(M), the "center" of the ellipsoid.
!
!    Input, real ( kind = rk ) R, the "radius" of the ellipsoid.
!
!    Output, real ( kind = rk ) X(M,N), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,m)
  integer i
  integer info
  integer j
  real ( kind = rk ) r
  real ( kind = rk ) u(m,m)
  real ( kind = rk ) v(m)
  real ( kind = rk ) x(m,n)
!
!  Get the Cholesky factor U.
!
  u(1:m,1:m) = a(1:m,1:m)

  call r8po_fa ( m, u, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ELLIPSOID_SAMPLE - Fatal error!'
    write ( *, '(a)' ) '  R8PO_FA reports that the matrix A '
    write ( *, '(a)' ) '  is not positive definite symmetric.'
    stop 1
  end if
!
!  Get the points Y that satisfy Y' * Y <= 1.
!
  call uniform_in_sphere01_map ( m, n, x )
!
!  Get the points Y that satisfy Y' * Y <= R * R.
!
  x(1:m,1:n) = r * x(1:m,1:n)
!
!  Solve U * X = Y.
!
  do j = 1, n
    call r8po_sl ( m, u, x(1:m,j) )
  end do
!
!  X = X + V.
!
  do i = 1, m
    x(i,1:n) = x(i,1:n) + v(i)
  end do

  return
end
function ellipsoid_volume ( m, a, v, r )

!*****************************************************************************80
!
!! ELLIPSOID_VOLUME returns the volume of an ellipsoid.
!
!  Discussion:
!
!    The points X in the ellipsoid are described by an M by M
!    positive definite symmetric matrix A, an M-dimensional point V,
!    and a "radius" R, such that
!      (X-V)' * A * (X-V) <= R * R
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 August 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the spatial dimension.
!
!    Input, real ( kind = rk ) A(M,M), the matrix that describes
!    the ellipsoid.  A must be symmetric and positive definite.
!
!    Input, real ( kind = rk ) V(M), the "center" of the ellipse.
!    The value of V is not actually needed by this function.
!
!    Input, real ( kind = rk ) R, the "radius" of the ellipse.
!
!    Output, real ( kind = rk ) ELLIPSOID_VOLUME, the volume of the ellipsoid.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m

  real ( kind = rk ) a(m,m)
  real ( kind = rk ) ellipsoid_volume
  real ( kind = rk ) hypersphere_unit_volume
  integer i
  integer info
  real ( kind = rk ) r
  real ( kind = rk ) sqrt_det
  real ( kind = rk ) u(m,m)
  real ( kind = rk ) v(m)

  call r8_fake_use ( v(1) )

  u(1:m,1:m) = a(1:m,1:m)

  call r8po_fa ( m, u, info )
  
  sqrt_det = 1.0D+00
  do i = 1, m
    sqrt_det = sqrt_det * u(i,i)
  end do

  ellipsoid_volume = r ** m * hypersphere_unit_volume ( m ) / sqrt_det

  return
end
function hypersphere_unit_volume ( m )

!*****************************************************************************80
!
!! HYPERSPHERE_UNIT_VOLUME: volume of a unit hypersphere in M dimensions.
!
!  Discussion:
!
!    The unit hypersphere in M dimensions satisfies:
!
!      sum ( 1 <= I <= M ) X(I) * X(I) = 1
!
!    Results for the first few values of DIM_NUM are:
!
!     M  Volume
!
!     1    2
!     2    1        * PI
!     3  ( 4 /   3) * PI
!     4  ( 1 /   2) * PI^2
!     5  ( 8 /  15) * PI^2
!     6  ( 1 /   6) * PI^3
!     7  (16 / 105) * PI^3
!     8  ( 1 /  24) * PI^4
!     9  (32 / 945) * PI^4
!    10  ( 1 / 120) * PI^5
!
!    For the unit sphere, Volume(M) = 2 * PI * Volume(M-2)/ M
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    13 August 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the spatial dimension.
!
!    Output, real ( kind = rk ) HYPERSPHERE_UNIT_VOLUME, the volume.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) hypersphere_unit_volume
  integer i
  integer m
  integer m2
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) value

  if ( mod ( m, 2 ) == 0 ) then
    m2 = m / 2
    value = r8_pi ** m2
    do i = 1, m2
      value = value / real ( i, kind = rk )
    end do
  else
    m2 = ( m - 1 ) / 2
    value = r8_pi ** m2 * 2.0D+00 ** m
    do i = m2 + 1, 2 * m2 + 1
      value = value / real ( i, kind = rk )
    end do
  end if

  hypersphere_unit_volume = value

  return
end
subroutine monomial_value ( m, n, e, x, value )

!*****************************************************************************80
!
!! MONOMIAL_VALUE evaluates a monomial.
!
!  Discussion:
!
!    This routine evaluates a monomial of the form
!
!      product ( 1 <= i <= m ) x(i)^e(i)
!
!    where the exponents are nonnegative integers.  Note that
!    if the combination 0^0 is encountered, it should be treated
!    as 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of points at which the
!    monomial is to be evaluated.
!
!    Input, integer E(M), the exponents.
!
!    Input, real ( kind = rk ) X(M,N), the point coordinates.
!
!    Output, real ( kind = rk ) VALUE(N), the value of the monomial.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer e(m)
  integer i
  real ( kind = rk ) value(n)
  real ( kind = rk ) x(m,n)

  value(1:n) = 1.0D+00

  do i = 1, m
    if ( 0 /= e(i) ) then
      value(1:n) = value(1:n) * x(i,1:n) ** e(i)
    end if
  end do

  return
end
subroutine r8_fake_use ( x )

!*****************************************************************************80
!
!! r8_fake_use pretends to use a variable.
!
!  Discussion:
!
!    Some compilers will issue a warning if a variable is unused.
!    Sometimes there's a good reason to include a variable in a program,
!    but not to use it.  Calling this function with that variable as
!    the argument will shut the compiler up.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 April 2020
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) X, the variable to be "used".
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) x

  if ( x /= x ) then
    write ( *, '(a)' ) '  r8_fake_use: variable is NAN.'
  end if

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
subroutine r8mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
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
!    14 June 2004
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
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  character ( len = * ) title

  call r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
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
  integer i2
  integer i2hi
  integer i2lo
  integer ihi
  integer ilo
  integer inc
  integer j
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

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(i8,6x)' ) i
    end do

    write ( *, '(''  Row   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc
        i = i2lo - 1 + i2
        write ( ctemp(i2), '(g14.6)' ) a(i,j)
      end do

      write ( *, '(i5,a,5a14)' ) j, ':', ( ctemp(i), i = 1, inc )

    end do

  end do

  return
end
subroutine r8po_fa ( n, a, info )

!*****************************************************************************80
!
!! R8PO_FA factors an R8PO matrix.
!
!  Discussion:
!
!    The R8PO storage format is used for a symmetric positive definite
!    matrix and its inverse.  (The Cholesky factor of a R8PO matrix is an
!    upper triangular matrix, so it will be in DGE storage format.)
!
!    Only the diagonal and upper triangle of the square array are used.
!    This same storage scheme is used when the matrix is factored by
!    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
!    is set to zero.
!
!    R8PO storage is used by LINPACK and LAPACK.
!
!    The positive definite symmetric matrix A has a Cholesky factorization
!    of the form:
!
!      A = R' * R
!
!    where R is an upper triangular matrix with positive elements on
!    its diagonal.  This routine overwrites the matrix A with its
!    factor R.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 March 2003
!
!  Author:
!
!    Original FORTRAN77 version by Jack Dongarra, Jim Bunch,
!    Cleve Moler, Pete Stewart.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real ( kind = rk ) A(N,N).
!    On input, the matrix in R8PO storage.
!    On output, the Cholesky factor R in DGE storage.
!
!    Output, integer INFO, error flag.
!    0, normal return.
!    K, error condition.  The principal minor of order K is not
!    positive definite, and the factorization was not completed.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,n)
  integer i
  integer info
  integer j
  integer k
  real ( kind = rk ) s

  do j = 1, n

    do k = 1, j - 1
      a(k,j) = ( a(k,j) - sum ( a(1:k-1,k) * a(1:k-1,j) ) ) / a(k,k)
    end do

    s = a(j,j) - sum ( a(1:j-1,j)**2 )

    if ( s <= 0.0D+00 ) then
      info = j
      return
    end if

    a(j,j) = sqrt ( s )

  end do

  info = 0
!
!  Since the Cholesky factor is stored in DGE format, be sure to
!  zero out the lower triangle.
!
  do i = 1, n
    do j = 1, i - 1
      a(i,j) = 0.0D+00
    end do
  end do

  return
end
subroutine r8po_sl ( n, a_lu, b )

!*****************************************************************************80
!
!! R8PO_SL solves an R8PO system factored by R8PO_FA.
!
!  Discussion:
!
!    The R8PO storage format is used for a symmetric positive definite
!    matrix and its inverse.  (The Cholesky factor of a R8PO matrix is an
!    upper triangular matrix, so it will be in DGE storage format.)
!
!    Only the diagonal and upper triangle of the square array are used.
!    This same storage scheme is used when the matrix is factored by
!    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
!    is set to zero.
!
!    R8PO storage is used by LINPACK and LAPACK.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 March 1999
!
!  Author:
!
!    Original FORTRAN77 version by Jack Dongarra, Jim Bunch,
!    Cleve Moler, Pete Stewart.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real ( kind = rk ) A_LU(N,N), the Cholesky factor from R8PO_FA.
!
!    Input/output, real ( kind = rk ) B(N).
!    On input, the right hand side.
!    On output, the solution vector.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a_lu(n,n)
  real ( kind = rk ) b(n)
  integer k
!
!  Solve R' * y = b.
!
  do k = 1, n
    b(k) = ( b(k) - sum ( b(1:k-1) * a_lu(1:k-1,k) ) ) / a_lu(k,k)
  end do
!
!  Solve R * x = y.
!
  do k = n, 1, -1
    b(k) = b(k) / a_lu(k,k)
    b(1:k-1) = b(1:k-1) - a_lu(1:k-1,k) * b(k)
  end do

  return
end
subroutine r8vec_normal_01 ( n, x )

!*****************************************************************************80
!
!! R8VEC_NORMAL_01 returns a unit pseudonormal R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!    This routine can generate a vector of values on one call.  It
!    has the feature that it should provide the same results
!    in the same order no matter how we break up the task.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of values desired.  
!
!    Output, real ( kind = rk ) X(N), a sample of the standard normal PDF.
!
!  Local:
!
!    Local, real ( kind = rk ) R(N+1), is used to store some uniform
!    random values.  Its dimension is N+1, but really it is only needed
!    to be the smallest even number greater than or equal to N.
!
!    Local, integer X_LO_INDEX, X_HI_INDEX, records the range of entries of
!    X that we need to compute.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer m
  real ( kind = rk ) r(n+1)
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) x(n)
  integer x_hi_index
  integer x_lo_index
!
!  Record the range of X we need to fill in.
!
  x_lo_index = 1
  x_hi_index = n
!
!  Maybe we don't need any more values.
!
  if ( x_hi_index - x_lo_index + 1 == 1 ) then

    call random_number ( harvest = r(1:2) )

    x(x_hi_index) = &
             sqrt ( -2.0D+00 * log ( r(1) ) ) * cos ( 2.0D+00 * r8_pi * r(2) )
!
!  If we require an even number of values, that's easy.
!
  else if ( mod ( x_hi_index - x_lo_index + 1, 2 ) == 0 ) then

    m = ( x_hi_index - x_lo_index + 1 ) / 2

    call random_number ( harvest = r(1:2*m) )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * cos ( 2.0D+00 * r8_pi * r(2:2*m:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * sin ( 2.0D+00 * r8_pi * r(2:2*m:2) )
!
!  If we require an odd number of values, we generate an even number,
!  and handle the last pair specially, storing one in X(N), and
!  saving the other for later.
!
  else

    x_hi_index = x_hi_index - 1

    m = ( x_hi_index - x_lo_index + 1 ) / 2 + 1

    call random_number ( harvest = r(1:2*m) )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * cos ( 2.0D+00 * r8_pi * r(2:2*m-2:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * sin ( 2.0D+00 * r8_pi * r(2:2*m-2:2) )

    x(n) = sqrt ( -2.0D+00 * log ( r(2*m-1) ) ) &
      * cos ( 2.0D+00 * r8_pi * r(2*m) )

  end if

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
subroutine uniform_in_sphere01_map ( m, n, x )

!*****************************************************************************80
!
!! UNIFORM_IN_SPHERE01_MAP maps uniform points into the unit sphere.
!
!  Discussion:
!
!    The sphere has center 0 and radius 1.
!
!    This routine is valid for any spatial dimension DIM_NUM.
!
!    We first generate a point ON the sphere, and then distribute it
!    IN the sphere.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Russell Cheng,
!    Random Variate Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998, pages 168.
!
!    Reuven Rubinstein,
!    Monte Carlo Optimization, Simulation, and Sensitivity
!    of Queueing Networks,
!    Krieger, 1992,
!    ISBN: 0894647644,
!    LC: QA298.R79.
!
!  Parameters:
!
!    Input, integer M, the dimension of the space.
!
!    Input, integer N, the number of points.
!
!    Output, real ( kind = rk ) X(M,N), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) exponent
  integer j
  real ( kind = rk ) norm
  real ( kind = rk ) r
  real ( kind = rk ) x(m,n)

  exponent = 1.0D+00 / real ( m, kind = rk )

  do j = 1, n
!
!  Fill a vector with normally distributed values.
!
    call r8vec_normal_01 ( m, x(1:m,j) )
!
!  Compute the length of the vector.
!
    norm = sqrt ( sum ( x(1:m,j)**2 ) )
!
!  Normalize the vector.
!
    x(1:m,j) = x(1:m,j) / norm
!
!  Now compute a value to map the point ON the sphere INTO the sphere.
!
    call random_number ( harvest = r )

    x(1:m,j) = r ** exponent * x(1:m,j)

  end do

  return
end
