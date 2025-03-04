subroutine differ_backward ( h, o, p, c, x )

!*****************************************************************************80
!
!! differ_backward() computes backward difference coefficients.
!
!  Discussion:
!
!    We determine coefficients C to approximate the derivative at X0
!    of order O and precision P, using equally spaced backward
!    differences, so that 
!
!      d^o f(x)/dx^o = sum ( 0 <= i <= o+p-1 ) c(i) f(x-ih) + O(h^(p))
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 November 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) H, the spacing.  0 < H.
!
!    Input, integer O, the order of the derivative to be 
!    approximated.  1 <= O.
!
!    Input, integer P, the order of the error, as a power of H.
!
!    Output, real ( kind = rk ) C(O+P), the coefficients.
!
!    Output, real ( kind = rk ) X(O+P), the evaluation points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer o
  integer p

  real ( kind = rk ) b(o+p)
  real ( kind = rk ) c(o+p)
  real ( kind = rk ) h
  integer i
  integer info
  integer job
  integer n
  real ( kind = rk ) r8_factorial
  real ( kind = rk ) x(o+p)

  n = o + p

  do i = 1, n
    x(i) = real ( i - n, kind = rk ) * h
  end do

  b(1:o+p) = 0.0D+00
  b(o+1) = 1.0D+00

  job = 0
  call r8vm_sl ( n, x, b, c, job, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIFFER_BACKWARD - Fatal error!'
    write ( *, '(a)' ) '  Vandermonde linear system is singular.'
    stop 1
  end if

  c(1:n) = c(1:n) * r8_factorial ( o )

  return
end
subroutine differ_central ( h, o, p, c, x )

!*****************************************************************************80
!
!! DIFFER_CENTRAL computes central difference coefficients.
!
!  Discussion:
!
!    We determine coefficients C to approximate the derivative at X0
!    of order O and precision P, using equally spaced central
!    differences, so that 
!
!      d^o f(x)/dx^o = sum ( 0 <= i <= o+p-1 ) c(i) f(x+(2*i-o-p+1)*h/2) 
!        + O(h^(p))
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 November 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) H, the spacing.  0 < H.
!
!    Input, integer O, the order of the derivative to 
!    be approximated.  1 <= O.
!
!    Input, integer P, the order of the error, as a power of H.
!
!    Output, real ( kind = rk ) C(O+P), the coefficients.
!
!    Output, real ( kind = rk ) X(O+P), the evaluation points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer o
  integer p

  real ( kind = rk ) b(o+p)
  real ( kind = rk ) c(o+p)
  real ( kind = rk ) h
  integer i
  integer info
  integer job
  integer n
  real ( kind = rk ) r8_factorial
  real ( kind = rk ) x(o+p)

  n = o + p

  do i = 1, n
    x(i) = real ( - n - 1 + 2 * i, kind = rk ) * h / 2.0D+00
  end do

  b(1:o+p) = 0.0D+00
  b(o+1) = 1.0D+00

  job = 0
  call r8vm_sl ( n, x, b, c, job, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIFFER_CENTRAL - Fatal error!'
    write ( *, '(a)' ) '  Vandermonde linear system is singular.'
    stop 1
  end if

  c(1:n) = c(1:n) * r8_factorial ( o )

  return
end
subroutine differ_forward ( h, o, p, c, x )

!*****************************************************************************80
!
!! DIFFER_FORWARD computes forward difference coefficients.
!
!  Discussion:
!
!    We determine coefficients C to approximate the derivative at X0
!    of order O and precision P, using equally spaced forward
!    differences, so that 
!
!      d^o f(x)/dx^o = sum ( 0 <= i <= o+p-1 ) c(i) f(x+ih) + O(h^(p))
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 November 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real H, the spacing.  0 < H.
!
!    Input, integer O, the order of the derivative to be approximated.
!    1 <= O.
!
!    Input, integer P, the order of the error, as a power of H.
!
!    Output, real C(O+P), the coefficients.
!
!    Output, real X(O+P), the evaluation points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer o
  integer p

  real ( kind = rk ) b(o+p)
  real ( kind = rk ) c(o+p)
  real ( kind = rk ) h
  integer i
  integer info
  integer job
  integer n
  real ( kind = rk ) r8_factorial
  real ( kind = rk ) x(o+p)

  n = o + p

  do i = 1, n
    x(i) = real ( i - 1, kind = rk ) * h
  end do

  b(1:o+p) = 0.0D+00
  b(o+1) = 1.0D+00

  job = 0
  call r8vm_sl ( n, x, b, c, job, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIFFER_FORWARD - Fatal error!'
    write ( *, '(a)' ) '  Vandermonde linear system is singular.'
    stop 1
  end if

  c(1:n) = c(1:n) * r8_factorial ( o )

  return
end
subroutine differ_inverse ( n, stencil, a )

!*****************************************************************************80
!
!! DIFFER_INVERSE returns the inverse of the DIFFER matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real ( kind = rk ) STENCIL(N), the values that define A.
!    The entries of STENCIL must be distinct.
!    No entry of STENCIL may be 0.
!
!    Output, real ( kind = rk ) A(N,N), the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,n)
  integer i
  integer indx
  integer j
  integer k
  real ( kind = rk ) stencil(n)

  do j = 1, n
    do i = 1, n
      if ( j == 1 ) then
        a(i,j) = 1.0D+00
      else
        a(i,j) = 0.0D+00
      end if
    end do
  end do

  do i = 1, n

    indx = 0

    do k = 1, n

      if ( k /= i ) then

        indx = indx + 1

        do j = indx + 1, 1, -1

          a(i,j) = - stencil(k) * a(i,j) / ( stencil(i) - stencil(k) )

          if ( 1 < j ) then
            a(i,j) = a(i,j) + a(i,j-1) / ( stencil(i) - stencil(k) )
          end if

        end do

      end if

    end do

  end do

  do i = 1, n
    a(i,1:n) = a(i,1:n) / stencil(i)
  end do

  return
end
subroutine differ_matrix ( n, stencil, a )

!*****************************************************************************80
!
!! DIFFER_MATRIX computes the stencil matrix from the stencil vector.
!
!  Discussion:
!
!    If N = 4, and STENCIL = ( -3, -2, -1, 1 ), then A will be
!
!    -3  -2  -1  1
!     9   4   1  1
!   -27  -8  -1  1
!    81  16   1  1
!
!    This matrix is a generalized form of a Vandermonde matrix A2:
!
!     1   1   1  1
!    -3  -2  -1  1
!     9   4   1  1
!   -27  -8  -1  1    
!
!    and if A * x = b, the A2 * x2 = b, where x2(i) = x(i) * stencil(i)
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
!    Input, integer N, the number of stencil points.
!
!    Input, real ( kind = rk ) STENCIL(N), the stencil vector.
!    The entries in this vector must be distinct.
!    No entry of STENCIL may be 0.
!
!    Output, real ( kind = rk ) A(N,N), the stencil matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  
  real ( kind = rk ) a(n,n)
  integer i
  real ( kind = rk ) stencil(n)

  a(1,1:n) = stencil(1:n)
  do i = 2, n
    a(i,1:n) = a(i-1,1:n) * stencil(1:n)
  end do

  return
end
subroutine differ_solve ( n, stencil, order, c )

!*****************************************************************************80
!
!! DIFFER_SOLVE solves for finite difference coefficients.
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
!    Input, integer N, the number of stencil points.
!
!    Input, real ( kind = rk ) STENCIL(N), the stencil vector.
!    The entries in this vector must be distinct.
!    No entry of STENCIL may be 0.
!
!    Input, integer ORDER, the order of the derivative to
!    be approximated.  1 <= ORDER <= N.
!
!    Output, real ( kind = rk ) C(N), the coefficients to be used
!    to multiply U(STENCIL(I))-U(0), so that the sum forms an
!    approximation to the derivative of order ORDER, with error 
!    of order H^(N+1-ORDER).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  
  real ( kind = rk ) a(n,n)
  real ( kind = rk ) b(n)
  real ( kind = rk ) c(n)
  integer info
  integer order
  real ( kind = rk ) stencil(n)

  call differ_matrix ( n, stencil, a )

  b(1:n) = 0.0D+00
  b(order) = 1.0D+00
!
!  Solve A * C = B.
!
  c(1:n) = b(1:n)
  call r8mat_fs ( n, a, c, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIFFER_SOLVE - Fatal error!'
    write ( *, '(a)' ) '  The DIFFER system is singular.'
    stop
  end if

  return
end
subroutine differ_stencil ( x0, o, p, x, c )

!*****************************************************************************80
!
!! DIFFER_STENCIL computes finite difference coefficients.
!
!  Discussion:
!
!    We determine coefficients C to approximate the derivative at X0
!    of order O and precision P, using finite differences, so that 
!
!      d^o f(x)/dx^o (x0) = sum ( 0 <= i <= o+p-1 ) c(i) f(x(i)) 
!        + O(h^(p))
!
!    where H is the maximum spacing between X0 and any X(I).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 November 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) X0, the point where the derivative is to 
!    be approximated.
!
!    Input, integer O, the order of the derivative to be 
!    approximated.  1 <= O.
!
!    Input, integer P, the order of the error, as a power of H.
!
!    Input, real ( kind = rk ) X(O+P), the evaluation points.
!
!    Output, real ( kind = rk ) C(O+P), the coefficients.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer o
  integer p

  real ( kind = rk ) b(o+p)
  real ( kind = rk ) c(o+p)
  real ( kind = rk ) dx(o+p)
  integer info
  integer job
  integer n
  real ( kind = rk ) r8_factorial
  real ( kind = rk ) x(o+p)
  real ( kind = rk ) x0

  n = o + p

  dx(1:n) = x(1:n) - x0

  b(1:o+p) = 0.0D+00
  b(o+1) = 1.0D+00

  job = 0
  call r8vm_sl ( n, dx, b, c, job, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIFFER_STENCIL - Fatal error!'
    write ( *, '(a)' ) '  Vandermonde linear system is singular.'
    stop 1
  end if

  c(1:n) = c(1:n) * r8_factorial ( o )

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
function r8_factorial ( n )

!*****************************************************************************80
!
!! R8_FACTORIAL computes the factorial of N.
!
!  Discussion:
!
!    factorial ( N ) = product ( 1 <= I <= N ) I
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the argument of the factorial function.
!    If N is less than 1, the function value is returned as 1.
!
!    Output, real ( kind = rk ) R8_FACTORIAL, the factorial of N.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r8_factorial
  integer i
  integer n

  r8_factorial = 1.0D+00

  do i = 1, n
    r8_factorial = r8_factorial * real ( i, kind = rk )
  end do

  return
end
subroutine r8mat_fs ( n, a, b, info )

!*****************************************************************************80
!
!! R8MAT_FS factors and solves a system with one right hand side.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    This routine differs from R8MAT_FSS in two ways:
!    * only one right hand side is allowed;
!    * the input matrix A is not modified.
!
!    This routine uses partial pivoting, but no pivot vector is required.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 January 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input/output, real ( kind = rk ) A(N,N).
!    On input, A is the coefficient matrix of the linear system.
!    On output, A is in unit upper triangular form, and
!    represents the U factor of an LU factorization of the
!    original coefficient matrix.
!
!    Input/output, real ( kind = rk ) B(N).
!    On input, the right hand side of the linear system.
!    On output, the solution of the linear systems.
!
!    Output, integer INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) a2(n,n)
  real ( kind = rk ) b(n)
  integer i
  integer info
  integer ipiv
  integer jcol
  real ( kind = rk ) piv
  real ( kind = rk ) row(n)
  real ( kind = rk ) t
  real ( kind = rk ) temp

  a2(1:n,1:n) = a(1:n,1:n)

  info = 0

  do jcol = 1, n
!
!  Find the maximum element in column I.
!
    piv = abs ( a2(jcol,jcol) )
    ipiv = jcol
    do i = jcol + 1, n
      if ( piv < abs ( a2(i,jcol) ) ) then
        piv = abs ( a2(i,jcol) )
        ipiv = i
      end if
    end do

    if ( piv == 0.0D+00 ) then
      info = jcol
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8MAT_FS - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      stop 1
    end if
!
!  Switch rows JCOL and IPIV, and B.
!
    if ( jcol /= ipiv ) then

      row(1:n) = a2(jcol,1:n)
      a2(jcol,1:n) = a2(ipiv,1:n)
      a2(ipiv,1:n) = row(1:n)

      t       = b(jcol)
      b(jcol) = b(ipiv)
      b(ipiv) = t

    end if
!
!  Scale the pivot row.
!
    a2(jcol,jcol+1:n) = a2(jcol,jcol+1:n) / a2(jcol,jcol)
    b(jcol) = b(jcol) / a2(jcol,jcol)
    a2(jcol,jcol) = 1.0D+00
!
!  Use the pivot row to eliminate lower entries in that column.
!
    do i = jcol + 1, n
      if ( a2(i,jcol) /= 0.0D+00 ) then
        temp = - a2(i,jcol)
        a2(i,jcol) = 0.0D+00
        a2(i,jcol+1:n) = a2(i,jcol+1:n) + temp * a2(jcol,jcol+1:n)
        b(i) = b(i) + temp * b(jcol)
      end if
    end do

  end do
!
!  Back solve.
!
  do jcol = n, 2, -1
    b(1:jcol-1) = b(1:jcol-1) - a2(1:jcol-1,jcol) * b(jcol)
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
subroutine r8vec2_print ( n, a1, a2, title )

!*****************************************************************************80
!
!! R8VEC2_PRINT prints an R8VEC2.
!
!  Discussion:
!
!    An R8VEC2 is a dataset consisting of N pairs of R8's, stored
!    as two separate vectors A1 and A2.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, real ( kind = rk ) A1(N), A2(N), the vectors to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a1(n)
  real ( kind = rk ) a2(n)
  integer i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i4,2x,g14.6,2x,g14.6)' ) i, a1(i), a2(i)
  end do

  return
end
subroutine r8vm_sl ( n, a, b, x, job, info )

!*****************************************************************************80
!
!! R8VM_SL solves an R8VM linear system.
!
!  Discussion:
!
!    The R8VM storage format is used for an M by N Vandermonde matrix.
!    An M by N Vandermonde matrix is defined by the values in its second
!    row, which will be written here as X(1:N).  The matrix has a first 
!    row of 1's, a second row equal to X(1:N), a third row whose entries
!    are the squares of the X values, up to the M-th row whose entries
!    are the (M-1)th powers of the X values.  The matrix can be stored
!    compactly by listing just the values X(1:N).
!
!    Vandermonde systems are very close to singularity.  The singularity
!    gets worse as N increases, and as any pair of values defining
!    the matrix get close.  Even a system as small as N = 10 will
!    involve the 9th power of the defining values.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    29 September 2003
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Gene Golub, Charles Van Loan,
!    Matrix Computations,
!    Third Edition,
!    Johns Hopkins, 1996.
!
!  Parameters:
!
!    Input, integer N, the number of rows and columns of 
!    the matrix.
!
!    Input, real ( kind = rk ) A(N), the R8VM matrix.
!
!    Input, real ( kind = rk ) B(N), the right hand side.
!
!    Output, real ( kind = rk ) X(N), the solution of the linear system.
!
!    Input, integer JOB, specifies the system to solve.
!    0, solve A * x = b.
!    nonzero, solve A' * x = b.
!
!    Output, integer INFO.
!    0, no error.
!    nonzero, at least two of the values in A are equal.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) b(n)
  integer i
  integer info
  integer j
  integer job
  real ( kind = rk ) x(n)
!
!  Check for explicit singularity.
!
  info = 0

  do j = 1, n - 1
    do i = j + 1, n
      if ( a(i) == a(j) ) then
        info = 1
        return
      end if
    end do
  end do

  x(1:n) = b(1:n)

  if ( job == 0 ) then

    do j = 1, n - 1
      do i = n, j + 1, -1
        x(i) = x(i) - a(j) * x(i-1)
      end do
    end do

    do j = n - 1, 1, -1

      do i = j + 1, n
        x(i) = x(i) / ( a(i) - a(i-j) )
      end do

      do i = j, n - 1
        x(i) = x(i) - x(i+1)
      end do

    end do

  else

    do j = 1, n - 1
      do i = n, j + 1, -1
        x(i) = ( x(i) - x(i-1) ) / ( a(i) - a(i-j) )
      end do
    end do

    do j = n - 1, 1, -1
      do i = j, n - 1
        x(i) = x(i) - x(i+1) * a(j)
      end do
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
