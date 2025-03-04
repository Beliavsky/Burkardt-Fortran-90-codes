subroutine cholesky_upper_error ( n, a, c, error_frobenius )

!*****************************************************************************80
!
!! cholesky_upper_error() determines the error in an upper Cholesky factor.
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
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real ( kind = rk ) A(N,N), the matrix.
!
!    Input, real ( kind = rk ) C(N,N), the upper triangular Cholesky factor.
!
!    Output, real ( kind = rk ) ERROR_FROBENIUS, the Frobenius norm
!    of the difference matrix A - C' * C.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) c(n,n)
  real ( kind = rk ) ctc(n,n)
  real ( kind = rk ) error_frobenius

  ctc(1:n,1:n) = matmul ( transpose ( c(1:n,1:n) ), c(1:n,1:n) )

  error_frobenius = sqrt ( sum ( ( a(1:n,1:n) - ctc(1:n,1:n) )**2 ) )

  return
end
subroutine eigen_error ( n, k, a, x, lambda, error_frobenius )

!*****************************************************************************80
!
!! EIGEN_ERROR determines the error in a (right) eigensystem.
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
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer K, the number of eigenvectors.
!    K is usually 1 or N.
!
!    Input, real ( kind = rk ) A(N,N), the matrix.
!
!    Input, real ( kind = rk ) X(N,K), the K eigenvectors.
!
!    Input, real ( kind = rk ) LAMBDA(K), the K eigenvalues.
!
!    Output, real ( kind = rk ) ERROR_FROBENIUS, the Frobenius norm
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
  real ( kind = rk ) x(n,k)

  c(1:n,1:k) = matmul ( a(1:n,1:n), x(1:n,1:k) )

  do j = 1, k
    c(1:n,j) = c(1:n,j) - lambda(j) * x(1:n,j)
  end do

  error_frobenius = sqrt ( sum ( c(1:n,1:k)**2 ) )

  return
end
subroutine inverse_error ( n, a, b, error_frobenius )

!*****************************************************************************80
!
!! INVERSE_ERROR determines the error in an inverse matrix.
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
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real ( kind = rk ) A(N,N), the matrix.
!
!    Input, real ( kind = rk ) B(N,N), the inverse.
!
!    Output, real ( kind = rk ) ERROR_FROBENIUS, the Frobenius norm
!    of (A*B-I) + (B*A-I).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) b(n,n)
  real ( kind = rk ) c(n,n)
  real ( kind = rk ) error_ab
  real ( kind = rk ) error_ba
  real ( kind = rk ) error_frobenius
  integer j

  c(1:n,1:n) = matmul ( a(1:n,1:n), b(1:n,1:n) )

  do j = 1, n
    c(j,j) = c(j,j) - 1.0D+00
  end do

  error_ab = sqrt ( sum ( c(1:n,1:n)**2 ) )

  c(1:n,1:n) = matmul ( b(1:n,1:n), a(1:n,1:n) )

  do j = 1, n
    c(j,j) = c(j,j) - 1.0D+00
  end do

  error_ba = sqrt ( sum ( c(1:n,1:n)**2 ) )

  error_frobenius = error_ab + error_ba

  return
end
subroutine l1dd_apply ( n, h, u, lu )

!*****************************************************************************80
!
!! L1DD_APPLY applies the 1D DD Laplacian to a vector.
!
!  Discussion:
!
!    The N grid points are assumed to be evenly spaced by H.
!
!    For N = 5, the discrete Laplacian with Dirichlet boundary conditions
!    at both ends of [0,6] is applied to a vector of 7 values, with a spacing
!    of H = 6/(N+1) = 1 at the points X:
!
!      0  1  2  3  4  5  6
!
!    and has the matrix form L:
!
!       2 -1  0  0  0
!      -1  2 -1  0  0
!       0 -1  2 -1  0
!       0  0 -1  2 -1
!       0  0  0 -1  2
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 October 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!    N must be at least 3.
!
!    Input, real ( kind = rk ) H, the spacing between points.
!
!    Input, real ( kind = rk ) U(N), the value at each point.
!
!    Output, real ( kind = rk ) LU(N), the Laplacian evaluated at each point.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) h
  integer i
  real ( kind = rk ) lu(n)
  real ( kind = rk ) u(n)

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'L1DD_APPLY - Fatal error!'
    write ( *, '(a)' ) '  N < 3.'
    stop 1
  end if

  i = 1
  lu(i) = ( 2.0D+00 * u(i) - u(i+1) ) / h / h
  do i = 2, n - 1
    lu(i) = ( - u(i-1) + 2.0D+00 * u(i) - u(i+1) ) / h / h
  end do
  i = n
  lu(i) = ( - u(i-1) + 2.0D+00 * u(i) ) / h / h

  return
end
subroutine l1dd_cholesky ( n, h, c )

!*****************************************************************************80
!
!! L1DD_CHOLESKY computes the Cholesky factor of the 1D DD Laplacian.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 October 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!    N must be at least 3.
!
!    Input, real ( kind = rk ) H, the spacing between points.
!
!    Output, real ( kind = rk ) C(N,N), the Cholesky factor.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) c(n,n)
  real ( kind = rk ) h
  integer i
  integer j

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'L1DD_CHOLESKY - Fatal error!'
    write ( *, '(a)' ) '  N < 3.'
    stop 1
  end if

  do j = 1, n
    do i = 1, n
      c(i,j) = 0.0D+00
    end do
  end do

  do i = 1, n
    c(i,i) = sqrt ( real ( i + 1, kind = rk ) ) &
           / sqrt ( real ( i, kind = rk ) )
  end do

  do i = 1, n - 1
    c(i,i+1) = - sqrt ( real ( i, kind = rk ) ) &
               / sqrt ( real ( i + 1, kind = rk ) )
  end do

  c(1:n,1:n) = c(1:n,1:n) / h

  return
end
subroutine l1dd_eigen ( n, h, v, lambda )

!*****************************************************************************80
!
!! L1DD_EIGEN returns eigeninformation for the 1D DD Laplacian.
!
!  Discussion:
!
!    The grid points are assumed to be evenly spaced by H.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 October 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) H, the spacing between points.
!
!    Output, real ( kind = rk ) V(N,N), the eigenvectors.
!
!    Output, real ( kind = rk ) LAMBDA(N), the eigenvalues.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) h
  integer i
  real ( kind = rk ) i_r8
  integer j
  real ( kind = rk ) j_r8
  real ( kind = rk ) lambda(n)
  real ( kind = rk ) n_r8
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) theta
  real ( kind = rk ) v(n,n)

  n_r8 = real ( n, kind = rk )

  do j = 1, n
    j_r8 = real ( j, kind = rk )
    theta = 0.5D+00 * pi * j_r8 / ( n_r8 + 1.0D+00 )
    lambda(j) = 4.0D+00 * ( sin ( theta ) / h ) ** 2
    do i = 1, n
      i_r8 = real ( i, kind = rk )
      theta = pi * i_r8 * j_r8 / ( n_r8 + 1.0D+00 )
      v(i,j) = sqrt ( 2.0D+00 / ( n_r8 + 1.0D+00 ) ) * sin ( theta )
    end do
  end do

  return
end
subroutine l1dd ( n, h, l )

!*****************************************************************************80
!
!! L1DD stores the 1D DD Laplacian as a full matrix.
!
!  Discussion:
!
!    The N grid points are assumed to be evenly spaced by H.
!
!    For N = 5, the discrete Laplacian with Dirichlet boundary conditions
!    at both ends of [0,6] has the matrix form L:
!
!       2 -1  0  0  0
!      -1  2 -1  0  0
!       0 -1  2 -1  0
!       0  0 -1  2 -1
!       0  0  0 -1  2
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 October 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!    N must be at least 3.
!
!    Input, real ( kind = rk ) H, the spacing between points.
!
!    Output, real ( kind = rk ) L(N,N), the Laplacian matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) h
  integer i
  integer j
  real ( kind = rk ) l(n,n)

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'L1DD - Fatal error!'
    write ( *, '(a)' ) '  N < 3.'
    stop 1
  end if

  do j = 1, n
    do i = 1, n
      l(i,j) = 0.0D+00
    end do
  end do

  i = 1
  l(i,i) =  2.0D+00 / h / h
  l(i,i+1) = -1.0D+00 / h / h

  do i = 2, n - 1
    l(i,i-1) = -1.0D+00 / h / h
    l(i,i) =    2.0D+00 / h / h
    l(i,i+1) = -1.0D+00 / h / h
  end do

  i = n
  l(i,i-1) = -1.0D+00 / h / h
  l(i,i) =    2.0D+00 / h / h

  return
end
subroutine l1dd_inverse ( n, h, l )

!*****************************************************************************80
!
!! L1DD_INVERSE stores the inverse of the 1D DD Laplacian.
!
!  Discussion:
!
!    The N grid points are assumed to be evenly spaced by H.
!
!    For N = 5, the discrete Laplacian with Dirichlet boundary conditions
!    at both ends of [0,6] has the matrix form L:
!
!       2 -1  0  0  0
!      -1  2 -1  0  0
!       0 -1  2 -1  0
!       0  0 -1  2 -1
!       0  0  0 -1  2
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 October 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!    N must be at least 3.
!
!    Input, real ( kind = rk ) H, the spacing between points.
!
!    Output, real ( kind = rk ) L(N,N), the inverse of the Laplacian matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) h
  integer i
  integer j
  real ( kind = rk ) l(n,n)

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'L1DD_INVERSE - Fatal error!'
    write ( *, '(a)' ) '  N < 3.'
    stop 1
  end if

  do j = 1, n
    do i = 1, n
      l(i,j) = real ( min ( i, j ) * ( n + 1 - max ( i, j ) ), kind = rk ) &
             * h * h / real ( n + 1, kind = rk )
    end do
  end do

  return
end
subroutine l1dd_lu ( n, h, l, u )

!*****************************************************************************80
!
!! L1DD_LU computes the LU factors of the 1D DD Laplacian.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 November 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!    N must be at least 3.
!
!    Input, real ( kind = rk ) H, the spacing between points.
!
!    Output, real ( kind = rk ) L(N,N), U(N,N), the LU factors.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) h
  integer i
  real ( kind = rk ) i_r8
  real ( kind = rk ) l(n,n)
  real ( kind = rk ) u(n,n)

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'L1DD_LU - Fatal error!'
    write ( *, '(a)' ) '  N < 3.'
    stop 1
  end if

  l(1:n,1:n) = 0.0D+00

  do i = 1, n
    l(i,i) = 1.0D+00
  end do

  do i = 2, n
    i_r8 = real ( i, kind = rk )
    l(i,i-1) = - ( i_r8 - 1.0D+00 ) / i_r8
  end do

  u(1:n,1:n) = 0.0D+00

  do i = 1, n
    i_r8 = real ( i, kind = rk )
    u(i,i) = ( i_r8 + 1.0D+00 ) / i_r8
  end do

  do i = 1, n - 1
    u(i,i+1) = - 1.0D+00
  end do

  u(1:n,1:n) = u(1:n,1:n) / h / h

  return
end
subroutine l1dn_apply ( n, h, u, lu )

!*****************************************************************************80
!
!! L1DN_APPLY applies the 1D DN Laplacian to a vector.
!
!  Discussion:
!
!    The N grid points are assumed to be evenly spaced by H.
!
!    For N = 5, the discrete Laplacian with left Dirichlet and right
!    Neumann condition on [0,6] has the matrix form L:
!
!       2 -1  0  0  0
!      -1  2 -1  0  0
!       0 -1  2 -1  0
!       0  0 -1  2 -1
!       0  0  0 -1  1
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 October 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!    N must be at least 3.
!
!    Input, real ( kind = rk ) H, the spacing between points.
!
!    Input, real ( kind = rk ) U(N), the value at each point.
!
!    Output, real ( kind = rk ) LU(N), the Laplacian evaluated at each point.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) h
  integer i
  real ( kind = rk ) lu(n)
  real ( kind = rk ) u(n)

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'L1DN_APPLY - Fatal error!'
    write ( *, '(a)' ) '  N < 3.'
    stop 1
  end if

  i = 1
  lu(i) = ( 2.0D+00 * u(i) - u(i+1) ) / h / h
  do i = 2, n - 1
    lu(i) = ( - u(i-1) + 2.0D+00 * u(i) - u(i+1) ) / h / h
  end do
  i = n
  lu(i) = ( - u(i-1) + u(i) ) / h / h

  return
end
subroutine l1dn_cholesky ( n, h, c )

!*****************************************************************************80
!
!! L1DN_CHOLESKY computes the Cholesky factor of the 1D DN Laplacian.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 October 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!    N must be at least 3.
!
!    Input, real ( kind = rk ) H, the spacing between points.
!
!    Output, real ( kind = rk ) C(N,N), the Cholesky factor.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) c(n,n)
  real ( kind = rk ) h
  integer i
  real ( kind = rk ) i_r8
  integer j

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'L1DN_CHOLESKY - Fatal error!'
    write ( *, '(a)' ) '  N < 3.'
    stop 1
  end if

  do j = 1, n
    do i = 1, n
      c(i,j) = 0.0D+00
    end do
  end do

  do i = 1, n - 1
    i_r8 = real ( i, kind = rk )
    c(i,i)   =   sqrt ( i_r8 + 1.0D+00 ) / sqrt ( i_r8 )
    c(i,i+1) = - sqrt ( i_r8 ) / sqrt ( i_r8 + 1.0D+00 )
  end do
  c(n,n) = 1.0D+00 / sqrt ( real ( n, kind = rk ) )

  c(1:n,1:n) = c(1:n,1:n) / h

  return
end
subroutine l1dn_eigen ( n, h, v, lambda )

!*****************************************************************************80
!
!! L1DN_EIGEN returns eigeninformation for the 1D DN Laplacian.
!
!  Discussion:
!
!    The grid points are assumed to be evenly spaced by H.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 October 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) H, the spacing between points.
!
!    Output, real ( kind = rk ) V(N,N), the eigenvectors.
!
!    Output, real ( kind = rk ) LAMBDA(N), the eigenvalues.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) h
  integer i
  real ( kind = rk ) i_r8
  integer j
  real ( kind = rk ) j_r8
  real ( kind = rk ) lambda(n)
  real ( kind = rk ) n_r8
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) theta
  real ( kind = rk ) v(n,n)

  n_r8 = real ( n, kind = rk )

  do j = 1, n
    j_r8 = real ( j, kind = rk )
    theta = pi * ( j_r8 - 0.5D+00 ) / ( 2.0D+00 * n_r8 + 1.0D+00 )
    lambda(j) = ( 2.0D+00 * sin ( theta ) / h ) ** 2
    do i = 1, n
      i_r8 = real ( i, kind = rk )
      theta = pi * i_r8 * ( 2.0D+00 * j_r8 - 1.0D+00 ) / &
        ( 2.0D+00 * n_r8 + 1.0D+00 )
      v(i,j) = sqrt ( 2.0D+00 / ( n_r8 + 0.5D+00 ) ) * sin ( theta )
    end do
  end do

  return
end
subroutine l1dn ( n, h, l )

!*****************************************************************************80
!
!! L1DN stores the 1D DN Laplacian as a full matrix.
!
!  Discussion:
!
!    The N grid points are assumed to be evenly spaced by H.
!
!    For N = 5, the discrete Laplacian with left Dirichlet and right
!    Neumann condition on [0,6] has the matrix form L:
!
!       2 -1  0  0  0
!      -1  2 -1  0  0
!       0 -1  2 -1  0
!       0  0 -1  2 -1
!       0  0  0 -1  1
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 October 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!    N must be at least 3.
!
!    Input, real ( kind = rk ) H, the spacing between points.
!
!    Output, real ( kind = rk ) L(N,N), the Laplacian matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) h
  integer i
  integer j
  real ( kind = rk ) l(n,n)

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'L1DN - Fatal error!'
    write ( *, '(a)' ) '  N < 3.'
    stop 1
  end if

  do j = 1, n
    do i = 1, n
      l(i,j) = 0.0D+00
    end do
  end do

  i = 1
  l(i,i)   =  2.0D+00 / h / h
  l(i,i+1) = -1.0D+00 / h / h

  do i = 2, n - 1
    l(i,i-1) = -1.0D+00 / h / h
    l(i,i) =    2.0D+00 / h / h
    l(i,i+1) = -1.0D+00 / h / h
  end do

  i = n
  l(i,i-1) = -1.0D+00 / h / h
  l(i,i) =    1.0D+00 / h / h

  return
end
subroutine l1dn_inverse ( n, h, l )

!*****************************************************************************80
!
!! L1DN_INVERSE stores the inverse of the 1D DN Laplacian.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 October 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!    N must be at least 3.
!
!    Input, real ( kind = rk ) H, the spacing between points.
!
!    Output, real ( kind = rk ) L(N,N), the inverse of the Laplacian matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) h
  integer i
  integer j
  real ( kind = rk ) l(n,n)

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'L1DN_INVERSE - Fatal error!'
    write ( *, '(a)' ) '  N < 3.'
    stop 1
  end if

  do j = 1, n
    do i = 1, n
      l(i,j) = real ( min ( i, j ), kind = rk ) * h * h
    end do
  end do

  return
end
subroutine l1dn_lu ( n, h, l, u )

!*****************************************************************************80
!
!! L1DN_LU computes the LU factors of the 1D DN Laplacian.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 November 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!    N must be at least 3.
!
!    Input, real ( kind = rk ) H, the spacing between points.
!
!    Output, real ( kind = rk ) L(N,N), U(N,N), the LU factors.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) h
  integer i
  real ( kind = rk ) i_r8
  real ( kind = rk ) l(n,n)
  real ( kind = rk ) u(n,n)

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'L1DN_LU - Fatal error!'
    write ( *, '(a)' ) '  N < 3.'
    stop 1
  end if

  l(1:n,1:n) = 0.0D+00

  do i = 1, n
    l(i,i) = 1.0D+00
  end do

  do i = 2, n
    i_r8 = real ( i, kind = rk )
    l(i,i-1) = - ( i_r8 - 1.0D+00 ) / i_r8
  end do

  u(1:n,1:n) = 0.0D+00

  do i = 1, n - 1
    i_r8 = real ( i, kind = rk )
    u(i,i) = ( i_r8 + 1.0D+00 ) / i_r8
  end do
  i = n
  i_r8 = real ( i, kind = rk )
  u(i,i) = 1.0D+00 / i_r8

  do i = 1, n - 1
    u(i,i+1) = - 1.0D+00
  end do

  u(1:n,1:n) = u(1:n,1:n) / h / h

  return
end
subroutine l1nd_apply ( n, h, u, lu )

!*****************************************************************************80
!
!! L1ND_APPLY applies the 1D ND Laplacian to a vector.
!
!  Discussion:
!
!    The N grid points are assumed to be evenly spaced by H.
!
!    For N = 5, the discrete Laplacian with left Neumann and right Dirichlet
!    boundary conditions on [0,6] has the matrix form L:
!
!       1 -1  0  0  0
!      -1  2 -1  0  0
!       0 -1  2 -1  0
!       0  0 -1  2 -1
!       0  0  0 -1  2
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 October 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!    N must be at least 3.
!
!    Input, real ( kind = rk ) H, the spacing between points.
!
!    Input, real ( kind = rk ) U(N), the value at each point.
!
!    Output, real ( kind = rk ) LU(N), the Laplacian evaluated at each point.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) h
  integer i
  real ( kind = rk ) lu(n)
  real ( kind = rk ) u(n)

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'L1ND_APPLY - Fatal error!'
    write ( *, '(a)' ) '  N < 3.'
    stop 1
  end if

  i = 1
  lu(i) = ( u(i) - u(i+1) ) / h / h
  do i = 2, n - 1
    lu(i) = ( - u(i-1) + 2.0D+00 * u(i) - u(i+1) ) / h / h
  end do
  i = n
  lu(i) = ( - u(i-1) + 2.0D+00 * u(i) ) / h / h

  return
end
subroutine l1nd_cholesky ( n, h, c )

!*****************************************************************************80
!
!! L1ND_CHOLESKY computes the Cholesky factor of the 1D ND Laplacian.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 October 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!    N must be at least 3.
!
!    Input, real ( kind = rk ) H, the spacing between points.
!
!    Output, real ( kind = rk ) C(N,N), the Cholesky factor.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) c(n,n)
  real ( kind = rk ) h
  integer i
  integer j

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'L1DN_CHOLESKY - Fatal error!'
    write ( *, '(a)' ) '  N < 3.'
    stop 1
  end if

  do j = 1, n
    do i = 1, n
      c(i,j) = 0.0D+00
    end do
  end do

  do i = 1, n
    c(i,i) = 1.0D+00
  end do

  do i = 1, n - 1
    c(i,i+1) = - 1.0D+00
  end do

  c(1:n,1:n) = c(1:n,1:n) / h

  return
end
subroutine l1nd_eigen ( n, h, v, lambda )

!*****************************************************************************80
!
!! L1ND_EIGEN returns eigeninformation for the 1D ND Laplacian.
!
!  Discussion:
!
!    The grid points are assumed to be evenly spaced by H.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 October 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) H, the spacing between points.
!
!    Output, real ( kind = rk ) V(N,N), the eigenvectors.
!
!    Output, real ( kind = rk ) LAMBDA(N), the eigenvalues.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) h
  integer i
  real ( kind = rk ) i_r8
  integer j
  real ( kind = rk ) j_r8
  real ( kind = rk ) lambda(n)
  real ( kind = rk ) n_r8
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) theta
  real ( kind = rk ) v(n,n)

  n_r8 = real ( n, kind = rk )

  do j = 1, n
    j_r8 = real ( j, kind = rk )
    theta = pi * ( j_r8 - 0.5D+00 ) / ( 2.0D+00 * n_r8 + 1.0D+00 )
    lambda(j) = 4.0D+00 * ( sin ( theta ) / h ) ** 2
    do i = 1, n
      i_r8 = real ( i, kind = rk )
      theta = pi * ( i_r8 - 0.5D+00 ) * ( 2.0D+00 * j_r8 - 1.0D+00 ) / &
        ( 2.0D+00 * n_r8 + 1.0D+00 )
      v(i,j) = sqrt ( 2.0D+00 / ( n_r8 + 0.5D+00 ) ) * cos ( theta )
    end do
  end do

  return
end
subroutine l1nd ( n, h, l )

!*****************************************************************************80
!
!! L1ND stores the 1D ND Laplacian as a full matrix.
!
!  Discussion:
!
!    The N grid points are assumed to be evenly spaced by H.
!
!    For N = 5, the discrete Laplacian with left Neumann and right Dirichlet
!    boundary conditions on [0,6] has the matrix form L:
!
!       1 -1  0  0  0
!      -1  2 -1  0  0
!       0 -1  2 -1  0
!       0  0 -1  2 -1
!       0  0  0 -1  2
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 October 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!    N must be at least 3.
!
!    Input, real ( kind = rk ) H, the spacing between points.
!
!    Output, real ( kind = rk ) L(N,N), the Laplacian matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) h
  integer i
  integer j
  real ( kind = rk ) l(n,n)

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'L1ND - Fatal error!'
    write ( *, '(a)' ) '  N < 3.'
    stop 1
  end if

  do j = 1, n
    do i = 1, n
      l(i,j) = 0.0D+00
    end do
  end do

  i = 1
  l(i,i)   =  1.0D+00 / h / h
  l(i,i+1) = -1.0D+00 / h / h

  do i = 2, n - 1
    l(i,i-1) = -1.0D+00 / h / h
    l(i,i) =    2.0D+00 / h / h
    l(i,i+1) = -1.0D+00 / h / h
  end do

  i = n
  l(i,i-1) = -1.0D+00 / h / h
  l(i,i) =    2.0D+00 / h / h

  return
end
subroutine l1nd_inverse ( n, h, l )

!*****************************************************************************80
!
!! L1ND_INVERSE stores the inverse of the 1D ND Laplacian.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 October 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!    N must be at least 3.
!
!    Input, real ( kind = rk ) H, the spacing between points.
!
!    Output, real ( kind = rk ) L(N,N), the inverse of the Laplacian matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) h
  integer i
  integer j
  real ( kind = rk ) l(n,n)

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'L1ND_INVERSE - Fatal error!'
    write ( *, '(a)' ) '  N < 3.'
    stop 1
  end if

  do j = 1, n
    do i = 1, n
      l(i,j) = real ( n + 1 - max ( i, j ), kind = rk ) * h * h
    end do
  end do

  return
end
subroutine l1nd_lu ( n, h, l, u )

!*****************************************************************************80
!
!! L1ND_LU computes the LU factors of the 1D ND Laplacian.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 November 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!    N must be at least 3.
!
!    Input, real ( kind = rk ) H, the spacing between points.
!
!    Output, real ( kind = rk ) L(N,N), U(N,N), the LU factors.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) h
  integer i
  real ( kind = rk ) l(n,n)
  real ( kind = rk ) u(n,n)

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'L1ND_LU - Fatal error!'
    write ( *, '(a)' ) '  N < 3.'
    stop 1
  end if

  l(1:n,1:n) = 0.0D+00

  do i = 1, n
    l(i,i) = 1.0D+00
  end do

  do i = 2, n
    l(i,i-1) = - 1.0D+00;
  end do

  u(1:n,1:n) = 0.0D+00

  do i = 1, n
    u(i,i) = 1.0D+00
  end do

  do i = 1, n - 1
    u(i,i+1) = - 1.0D+00
  end do

  u(1:n,1:n) = u(1:n,1:n) / h / h

  return
end
subroutine l1nn_apply ( n, h, u, lu )

!*****************************************************************************80
!
!! L1NN_APPLY applies the 1D NN Laplacian to a vector.
!
!  Discussion:
!
!    The N grid points are assumed to be evenly spaced by H.
!
!    For N = 5, the discrete Laplacian with left Neumann and right Neumann
!    boundary conditions on [0,6] has the matrix form L:
!
!       1 -1  0  0  0
!      -1  2 -1  0  0
!       0 -1  2 -1  0
!       0  0 -1  2 -1
!       0  0  0 -1  1
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 October 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!    N must be at least 3.
!
!    Input, real ( kind = rk ) H, the spacing between points.
!
!    Input, real ( kind = rk ) U(N), the value at each point.
!
!    Output, real ( kind = rk ) LU(N), the Laplacian evaluated at each point.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) h
  integer i
  real ( kind = rk ) lu(n)
  real ( kind = rk ) u(n)

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'L1NN_APPLY - Fatal error!'
    write ( *, '(a)' ) '  N < 3.'
    stop 1
  end if

  i = 1
  lu(i) = ( u(i) - u(i+1) ) / h / h
  do i = 2, n - 1
    lu(i) = ( - u(i-1) + 2.0D+00 * u(i) - u(i+1) ) / h / h
  end do
  i = n
  lu(i) = ( - u(i-1) +  u(i) ) / h / h

  return
end
subroutine l1nn_cholesky ( n, h, c )

!*****************************************************************************80
!
!! L1NN_CHOLESKY computes the Cholesky factor of the 1D NN Laplacian.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 October 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!    N must be at least 3.
!
!    Input, real ( kind = rk ) H, the spacing between points.
!
!    Output, real ( kind = rk ) C(N,N), the Cholesky factor.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) c(n,n)
  real ( kind = rk ) h
  integer i
  integer j

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'L1NN_CHOLESKY - Fatal error!'
    write ( *, '(a)' ) '  N < 3.'
    stop 1
  end if

  do j = 1, n
    do i = 1, n
      c(i,j) = 0.0D+00
    end do
  end do

  do i = 1, n - 1
    c(i,i)   = + 1.0D+00
    c(i,i+1) = - 1.0D+00
  end do

  c(1:n,1:n) = c(1:n,1:n) / h

  return
end
subroutine l1nn_eigen ( n, h, v, lambda )

!*****************************************************************************80
!
!! L1NN_EIGEN returns eigeninformation for the 1D NN Laplacian.
!
!  Discussion:
!
!    The grid points are assumed to be evenly spaced by H.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 October 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) H, the spacing between points.
!
!    Output, real ( kind = rk ) V(N,N), the eigenvectors.
!
!    Output, real ( kind = rk ) LAMBDA(N), the eigenvalues.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) h
  integer i
  real ( kind = rk ) i_r8
  integer j
  real ( kind = rk ) j_r8
  real ( kind = rk ) lambda(n)
  real ( kind = rk ) n_r8
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) theta
  real ( kind = rk ) v(n,n)

  n_r8 = real ( n, kind = rk )

  do j = 1, n
    j_r8 = real ( j, kind = rk )
    theta = pi * ( j_r8 - 1.0D+00 ) / ( 2.0D+00 * n_r8 )
    lambda(j) = 4.0D+00 * ( sin ( theta ) / h ) ** 2
    if ( j == 1 ) then
      do i = 1, n
        v(i,j) = sqrt ( n_r8 )
      end do
    else
      do i = 1, n
        i_r8 = real ( i, kind = rk )
        theta = pi * ( i_r8 - 0.5D+00 ) * ( j_r8 - 1.0D+00 ) / n_r8
        v(i,j) = sqrt ( 2.0D+00 / n_r8 ) * cos ( theta )
      end do
    end if
  end do

  return
end
subroutine l1nn ( n, h, l )

!*****************************************************************************80
!
!! L1NN stores the 1D NN Laplacian as a full matrix.
!
!  Discussion:
!
!    The N grid points are assumed to be evenly spaced by H.
!
!    For N = 5, the discrete Laplacian with Neumann boundary conditions
!    at both ends of [0,6] has the matrix form L:
!
!       1 -1  0  0  0
!      -1  2 -1  0  0
!       0 -1  2 -1  0
!       0  0 -1  2 -1
!       0  0  0 -1  1
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 October 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!    N must be at least 3.
!
!    Input, real ( kind = rk ) H, the spacing between points.
!
!    Output, real ( kind = rk ) L(N,N), the Laplacian matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) h
  integer i
  integer j
  real ( kind = rk ) l(n,n)

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'L1NN - Fatal error!'
    write ( *, '(a)' ) '  N < 3.'
    stop 1
  end if

  do j = 1, n
    do i = 1, n
      l(i,j) = 0.0D+00
    end do
  end do

  i = 1
  l(i,i) =  1.0D+00 / h / h
  l(i,i+1) = -1.0D+00 / h / h

  do i = 2, n - 1
    l(i,i-1) = -1.0D+00 / h / h
    l(i,i) =    2.0D+00 / h / h
    l(i,i+1) = -1.0D+00 / h / h
  end do

  i = n
  l(i,i-1) = -1.0D+00 / h / h
  l(i,i) =    1.0D+00 / h / h

  return
end
subroutine l1nn_lu ( n, h, l, u )

!*****************************************************************************80
!
!! L1NN_LU computes the LU factors of the 1D NN Laplacian.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 November 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!    N must be at least 3.
!
!    Input, real ( kind = rk ) H, the spacing between points.
!
!    Output, real ( kind = rk ) L(N,N), U(N,N), the LU factors.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) h
  integer i
  real ( kind = rk ) l(n,n)
  real ( kind = rk ) u(n,n)

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'L1NN_LU - Fatal error!'
    write ( *, '(a)' ) '  N < 3.'
    stop 1
  end if

  l(1:n,1:n) = 0.0D+00

  do i = 1, n
    l(i,i) = 1.0D+00
  end do

  do i = 2, n
    l(i,i-1) = - 1.0D+00;
  end do

  u(1:n,1:n) = 0.0D+00

  do i = 1, n - 1
    u(i,i) = 1.0D+00
  end do
  u(n,n) = 0.0D+00

  do i = 1, n - 1
    u(i,i+1) = - 1.0D+00
  end do

  u(1:n,1:n) = u(1:n,1:n) / h / h

  return
end
subroutine l1pp_apply ( n, h, u, lu )

!*****************************************************************************80
!
!! L1PP_APPLY applies the 1D PP Laplacian to a vector.
!
!  Discussion:
!
!    The N grid points are assumed to be evenly spaced by H.
!
!    For N = 5, the discrete Laplacian with periodic boundary conditions
!    on [0,6] has the matrix form L:
!
!       2 -1  0  0 -1
!      -1  2 -1  0  0
!       0 -1  2 -1  0
!       0  0 -1  2 -1
!      -1  0  0 -1  2
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 October 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!    N must be at least 3.
!
!    Input, real ( kind = rk ) H, the spacing between points.
!
!    Input, real ( kind = rk ) U(N), the value at each point.
!
!    Output, real ( kind = rk ) LU(N), the Laplacian evaluated at each point.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) h
  integer i
  real ( kind = rk ) lu(n)
  real ( kind = rk ) u(n)

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'L1PP_APPLY - Fatal error!'
    write ( *, '(a)' ) '  N < 3.'
    stop 1
  end if

  i = 1
  lu(i) = ( - u(n) + 2.0D+00 * u(i) - u(i+1) ) / h / h
  do i = 2, n - 1
    lu(i) = ( - u(i-1) + 2.0D+00 * u(i) - u(i+1) ) / h / h
  end do
  i = n
  lu(i) = ( - u(i-1) + 2.0D+00 * u(i) - u(1) ) / h / h

  return
end
subroutine l1pp_cholesky ( n, h, c )

!*****************************************************************************80
!
!! L1PP_CHOLESKY computes the Cholesky factor of the 1D PP Laplacian.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 October 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!    N must be at least 3.
!
!    Input, real ( kind = rk ) H, the spacing between points.
!
!    Output, real ( kind = rk ) C(N,N), the Cholesky factor.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) c(n,n)
  real ( kind = rk ) h
  integer i
  real ( kind = rk ) i_r8
  integer j

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'L1PP_CHOLESKY - Fatal error!'
    write ( *, '(a)' ) '  N < 3.'
    stop 1
  end if

  do j = 1, n
    do i = 1, n
      c(i,j) = 0.0D+00
    end do
  end do

  do i = 1, n - 1
    i_r8 = real ( i, kind = rk )
    c(i,i) = sqrt ( i_r8 + 1.0D+00 ) &
           / sqrt ( i_r8 )
  end do

  do i = 1, n - 2
    i_r8 = real ( i, kind = rk )
    c(i,i+1) = - i_r8 / ( i_r8 + 1.0D+00 ) &
             * sqrt ( i_r8 + 1.0D+00 )  &
             / sqrt ( i_r8 )
  end do

  do i = 1, n - 2
    i_r8 = real ( i, kind = rk )
    c(i,n) = - 1.0D+00 / ( i_r8 + 1.0D+00 ) &
             * sqrt ( i_r8 + 1.0D+00 )  &
             / sqrt ( i_r8 )
  end do

  i = n - 1
  i_r8 = real ( i, kind = rk )
  c(i,n) = - real ( n, kind = rk ) / ( i_r8 + 1.0D+00 ) &
           * sqrt ( i_r8 + 1.0D+00 )  &
           / sqrt ( i_r8 )

  c(1:n,1:n) = c(1:n,1:n) / h

  return
end
subroutine l1pp_eigen ( n, h, v, lambda )

!*****************************************************************************80
!
!! L1PP_EIGEN returns eigeninformation for the 1D PP Laplacian.
!
!  Discussion:
!
!    The grid points are assumed to be evenly spaced by H.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 October 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) H, the spacing between points.
!
!    Output, real ( kind = rk ) V(N,N), the eigenvectors.
!
!    Output, real ( kind = rk ) LAMBDA(N), the eigenvalues.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) h
  integer i
  real ( kind = rk ) i_r8
  integer j
  real ( kind = rk ) j_r8
  real ( kind = rk ) lambda(n)
  real ( kind = rk ) n_r8
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) s
  real ( kind = rk ) theta
  real ( kind = rk ) v(n,n)

  n_r8 = real ( n, kind = rk )

  do j = 1, n

    j_r8 = real ( j, kind = rk )
    if ( mod ( j, 2 ) == 1 ) then
      theta = pi * ( j_r8 - 1.0D+00 ) / ( 2.0D+00 * n_r8 )
    else
      theta = pi *   j_r8             / ( 2.0D+00 * n_r8 )
    end if
    lambda(j) = 4.0D+00 * ( sin ( theta ) / h ) ** 2

    if ( mod ( j, 2 ) == 1 ) then
      if ( j == 1 ) then
        do i = 1, n
          v(i,j) = 1.0D+00 / sqrt ( n_r8 )
        end do
      else
        do i = 1, n
          i_r8 = real ( i, kind = rk )
          theta = pi * ( i_r8 - 0.5D+00 ) * ( j_r8 - 1.0D+00 ) /  n_r8
          v(i,j) = sqrt ( 2.0D+00 / n_r8 ) * cos ( theta )
        end do
      end if
    else
      if ( j == n ) then
        s = - 1.0D+00 / sqrt ( n_r8 )
        do i = 1, n
          v(i,j) = s
          s = - s
        end do
      else
        do i = 1, n
          i_r8 = real ( i, kind = rk )
          theta = pi * ( i_r8 - 0.5D+00 ) * j_r8 / n_r8
          v(i,j) = sqrt ( 2.0D+00 / n_r8 ) * sin ( theta )
        end do
      end if
    end if

  end do

  return
end
subroutine l1pp ( n, h, l )

!*****************************************************************************80
!
!! L1PP stores the 1D PP Laplacian as a full matrix.
!
!  Discussion:
!
!    The N grid points are assumed to be evenly spaced by H.
!
!    For N = 5, the discrete Laplacian with periodic boundary conditions
!    has the matrix form L:
!
!       2 -1  0  0 -1
!      -1  2 -1  0  0
!       0 -1  2 -1  0
!       0  0 -1  2 -1
!      -1  0  0 -1  2
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 October 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!    N must be at least 3.
!
!    Input, real ( kind = rk ) H, the spacing between points.
!
!    Output, real ( kind = rk ) L(N,N), the Laplacian matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) h
  integer i
  integer j
  real ( kind = rk ) l(n,n)

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'L1PP - Fatal error!'
    write ( *, '(a)' ) '  N < 3.'
    stop 1
  end if

  do j = 1, n
    do i = 1, n
      l(i,j) = 0.0D+00
    end do
  end do

  i = 1
  l(i,i)   =  2.0D+00 / h / h
  l(i,i+1) = -1.0D+00 / h / h
  l(i,n)   = -1.0D+00 / h / h

  do i = 2, n - 1
    l(i,i-1) = -1.0D+00 / h / h
    l(i,i) =    2.0D+00 / h / h
    l(i,i+1) = -1.0D+00 / h / h
  end do

  i = n
  l(i,1) =   -1.0D+00 / h / h
  l(i,i-1) = -1.0D+00 / h / h
  l(i,i) =    2.0D+00 / h / h

  return
end
subroutine l1pp_lu ( n, h, l, u )

!*****************************************************************************80
!
!! L1PP_LU computes the LU factors of the 1D PP Laplacian.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 November 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!    N must be at least 3.
!
!    Input, real ( kind = rk ) H, the spacing between points.
!
!    Output, real ( kind = rk ) L(N,N), U(N,N), the LU factors.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) h
  integer i
  real ( kind = rk ) i_r8
  real ( kind = rk ) l(n,n)
  real ( kind = rk ) u(n,n)

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'L1PP_LU - Fatal error!'
    write ( *, '(a)' ) '  N < 3.'
    stop 1
  end if

  l(1:n,1:n) = 0.0D+00

  do i = 1, n
    l(i,i) = 1.0D+00
  end do

  do i = 2, n - 1
    i_r8 = real ( i, kind = rk )
    l(i,i-1) = - ( i_r8 - 1.0D+00 ) / i_r8
    l(n,i-1) =          - 1.0D+00   / i_r8
  end do
  l(n,n-1) = -1.0D+00

  u(1:n,1:n) = 0.0D+00

  do i = 1, n - 2
    i_r8 = real ( i, kind = rk )
    u(i,i) =   ( i_r8 + 1.0D+00 ) / i_r8
    u(i,i+1) = - 1.0D+00
    u(i,n) =   - 1.0D+00 / i_r8
  end do

  i = n - 1
  i_r8 = real ( i, kind = rk )
  u(i,i) =   ( i_r8 + 1.0D+00 ) / i_r8
  u(i,i+1) = - ( i_r8 + 1.0D+00 ) / i_r8

  i = n
  u(i,i) = 0.0D+00

  u(1:n,1:n) = u(1:n,1:n) / h / h

  return
end
subroutine laplacian3_interval_uneven_vector ( bc, n, x, u, lu )

!*****************************************************************************80
!
!! laplacian3_interval_uneven_vector() applies the 1D Discrete Laplacian to a vector.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 September 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer bc(2): defines the boundary conditions.
!    0: Dirichlet
!    1: Neumann
!    2: Robin
!
!    integer n: the number of entries in the vector.
!
!    real x(n): the coordinates of the points.
!
!    real u(n): the function value at each point.
!
!  Output:
!
!    real lu(n): the Laplacian at the interior points.
!    The boundary condition at the endpoints.  
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) alpha
  integer bc(2)
  real ( kind = rk ) beta
  real ( kind = rk ) dxl
  real ( kind = rk ) dxr
  real ( kind = rk ) gamma
  integer i
  real ( kind = rk ) lu(n)
  real ( kind = rk ) u(n)
  real ( kind = rk ) x(n)

  if ( bc(1) == 0 ) then
    lu(1) = u(1)
  else if ( bc(1) == 1 ) then
    lu(1) = ( u(2) - u(1) ) / ( x(2) - x(1) )
  else if ( bc(1) == 2 ) then
    dxl = x(n) - x(n-1)
    dxr = x(2) - x(1)
    alpha =  2.0 *   dxr         / dxl / ( dxl + dxr ) / dxr
    beta = - 2.0 * ( dxr + dxl ) / dxl / ( dxl + dxr ) / dxr
    gamma =  2.0 *         dxl   / dxl / ( dxl + dxr ) / dxr
    lu(1) = alpha * u(n-1) + beta * u(1) + gamma * u(2)
  end if

  do i = 2, n - 1

    dxl = x(i)   - x(i-1)
    dxr = x(i+1) - x(i)

    alpha =  2.0 *   dxr         / dxl / ( dxl + dxr ) / dxr
    beta = - 2.0 * ( dxr + dxl ) / dxl / ( dxl + dxr ) / dxr
    gamma =  2.0 *         dxl   / dxl / ( dxl + dxr ) / dxr

    lu(i) = alpha * u(i-1) + beta * u(i) + gamma * u(i+1)

  end do

  if ( bc(2) == 0 ) then
    lu(n) = u(n)
  else if ( bc(2) == 1 ) then
    lu(n) = ( u(n) - u(n-1) ) / ( x(n) - x(n-1) )
  else if ( bc(2) == 2 ) then
    dxl = x(n) - x(n-1)
    dxr = x(2) - x(1)
    alpha =  2.0 *   dxr         / dxl / ( dxl + dxr ) / dxr
    beta = - 2.0 * ( dxr + dxl ) / dxl / ( dxl + dxr ) / dxr
    gamma =  2.0 *         dxl   / dxl / ( dxl + dxr ) / dxr
    lu(n) = alpha * u(n-1) + beta * u(n) + gamma * u(2)
  end if

  return
end
subroutine lu_error ( n, a, l, u, error_frobenius )

!*****************************************************************************80
!
!! LU_ERROR determines the error in an LU factorization.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 November 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real ( kind = rk ) A(N,N), the matrix.
!
!    Input, real ( kind = rk ) L(N,N), U(N,N), the LU factors.
!
!    Output, real ( kind = rk ) ERROR_FROBENIUS, the Frobenius norm
!    of the difference matrix A - L * U.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) error_frobenius
  real ( kind = rk ) l(n,n)
  real ( kind = rk ) lu(n,n)
  real ( kind = rk ) u(n,n)

  lu(1:n,1:n) = matmul ( l(1:n,1:n), u(1:n,1:n) )

  error_frobenius = sqrt ( sum ( ( a(1:n,1:n) - lu(1:n,1:n) )**2 ) )

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
!  Parameters:
!
!    None
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

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
