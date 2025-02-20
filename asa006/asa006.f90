subroutine cholesky ( a, n, nn, u, nullty, ifault )

!*****************************************************************************80
!
!! cholesky() computes the Cholesky factorization of an SPD matrix.
!
!  Discussion:
!
!    For a symmetric positive definite matrix A, the Cholesky factor U
!    is an upper triangular matrix such that A = U' * U.
!
!    This routine was originally named "CHOL", but that conflicted with
!    a built in MATLAB routine name.
!
!    The missing initialization "II = 0" has been added to the code.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 August 2021
!
!  Author:
!
!    Original FORTRAN77 version by Michael Healy.
!    Modifications by AJ Miller.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Michael Healy,
!    Algorithm AS 6:
!    Triangular decomposition of a symmetric matrix,
!    Applied Statistics,
!    Volume 17, Number 2, 1968, pages 195-197.
!
!  Input:
!
!    real ( kind = rk ) A((N*(N+1))/2), a symmetric positive definite matrix
!    stored by rows in lower triangular form as a one dimensional array,
!    in the sequence
!    A(1,1),
!    A(2,1), A(2,2),
!    A(3,1), A(3,2), A(3,3), and so on.
!
!    integer N, the order of A.
!
!    integer NN, the dimension of the array used to store A,
!    which should be at least (N*(N+1))/2.
!
!  Output:
!
!    real ( kind = rk ) U((N*(N+1))/2), an upper triangular matrix,
!    stored by columns, which is the Cholesky factor of A.  The program is
!    written in such a way that A and U can share storage.
!
!    integer NULLTY, the rank deficiency of A.  If NULLTY
!    is zero, the matrix is judged to have full rank.
!
!    integer IFAULT, an error indicator.
!    0, no error was detected;
!    1, if N < 1;
!    2, if A is not positive semi-definite.
!    3, NN < (N*(N+1))/2.
!
!  Local:
!
!    real ( kind = rk ) ETA, should be set equal to the smallest positive
!    value such that 1.0 + ETA is calculated as being greater than 1.0 in the
!    accuracy being used.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nn

  real ( kind = rk ) a(nn)
  real ( kind = rk ), parameter :: eta = 1.0D-09
  integer i
  integer icol
  integer ifault
  integer ii
  integer irow
  integer j
  integer k
  integer kk
  integer l
  integer m
  integer nullty
  real ( kind = rk ) u(nn)
  real ( kind = rk ) w
  real ( kind = rk ) x

  ifault = 0
  nullty = 0

  if ( n <= 0 ) then
    ifault = 1
    return
  end if

  if ( nn < ( n * ( n + 1 ) ) / 2 ) then
    ifault = 3
    return
  end if

  j = 1
  k = 0
  ii = 0
!
!  Factorize column by column, ICOL = column number.
!
  do icol = 1, n

    ii = ii + icol
    x = eta * eta * a(ii)
    l = 0
    kk = 0
!
!  IROW = row number within column ICOL.
!
    do irow = 1, icol

      kk = kk + irow
      k = k + 1
      w = a(k)
      m = j

      do i = 1, irow - 1
        l = l + 1
        w = w - u(l) * u(m)
        m = m + 1
      end do

      l = l + 1

      if ( irow == icol ) then
        exit
      end if

      if ( u(l) /= 0.0D+00 ) then

        u(k) = w / u(l)

      else

        u(k) = 0.0D+00

        if ( abs ( x * a(k) ) < w * w ) then
          ifault = 2
          return
        end if

      end if

    end do
!
!  End of row, estimate relative accuracy of diagonal element.
!
    if ( abs ( w ) <= abs ( eta * a(k) ) ) then

      u(k) = 0.0D+00
      nullty = nullty + 1

    else

      if ( w < 0.0D+00 ) then
        ifault = 2
        return
      end if

      u(k) = sqrt ( w )

    end if

    j = j + icol

  end do

  return
end
subroutine subchl ( a, b, n, u, nullty, ifault, ndim, det )

!*****************************************************************************80
!
!! subchl() computes the Cholesky factorization of a (submatrix of an) SPD matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 August 2021
!
!  Author:
!
!    Original FORTRAN77 version by Michael Healy, PR Freeman
!    FORTRAN90 version by  John Burkardt
!
!  Reference:
!
!    PR Freeman,
!    Remark AS R44:
!    A Remark on AS 6 and AS7: Triangular decomposition of a symmetric matrix
!    and Inversion of a positive semi-definite symmetric matrix,
!    Applied Statistics,
!    Volume 31, Number 3, 1982, pages 336-339.
!
!    Michael Healy,
!    Algorithm AS 6:
!    Triangular decomposition of a symmetric matrix,
!    Applied Statistics,
!    Volume 17, Number 2, 1968, pages 195-197.
!
!  Input:
!
!    real ( kind = rk ) A((M*(M+1))/2), a symmetric positive definite matrix
!    stored by rows in lower triangular form as a one dimensional array,
!    in the sequence
!    A(1,1),
!    A(2,1), A(2,2),
!    A(3,1), A(3,2), A(3,3), and so on.
!    In the simplest case, M, the order of A, is equal to N.
!
!    integer B(N), indicates the order in which the
!    rows and columns of A are to be used.  In the simplest case,
!    B = (1,2,3...,N).
!
!    integer N, the order of the matrix, that is,
!    the matrix formed by using B to select N rows and columns of A.
!
!    integer NDIM, the dimension of A and U, which might
!    be presumed to be (N*(N+1))/2.
!
!  Output:
!
!    real ( kind = rk ) U((N*(N+1))/2), an upper triangular matrix,
!    stored by columns, which is the Cholesky factor of A.  The program is
!    written in such a way that A and U can share storage.
!
!    integer NULLTY, the rank deficiency of A.
!    If NULLTY is zero, the matrix is judged to have full rank.
!
!    integer IFAULT, an error indicator.
!    0, no error was detected;
!    1, if N < 1;
!    2, if A is not positive semi-definite.
!
!    real ( kind = rk ) DET, the determinant of the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer ndim

  real ( kind = rk ) a(ndim)
  integer b(n)
  real ( kind = rk ) det
  real ( kind = rk ), parameter :: eta = 1.0D-09
  integer i
  integer icol
  integer ifault
  integer ii
  integer ij
  integer irow
  integer j
  integer jj
  integer k
  integer kk
  integer l
  integer m
  integer nullty
  real ( kind = rk ) u(ndim)
  real ( kind = rk ) w
  real ( kind = rk ) x

  ifault = 0
  nullty = 0
  det = 1.0D+00

  if ( n <= 0 ) then
    ifault = 1
    return
  end if

  ifault = 2
  j = 1
  k = 0

  do icol = 1, n

    ij = ( b(icol) * ( b(icol) - 1 ) ) / 2
    ii = ij + b(icol)
    x = eta * eta * a(ii)
    l = 0

    do irow = 1, icol

      kk = ( b(irow) * ( b(irow) + 1 ) ) / 2
      k = k + 1
      jj = ij + b(irow)
      w = a(jj)
      m = j

      do i = 1, irow - 1
        l = l + 1
        w = w - u(l) * u(m)
        m = m + 1
      end do

      l = l + 1

      if ( irow == icol ) then
        exit
      end if

      if ( u(l) /= 0.0D+00 ) then

        u(k) = w / u(l)

      else

        if ( abs ( x * a(kk) ) < w * w ) then
          ifault = 2
          return
        end if

        u(k) = 0.0D+00

      end if

    end do

    if ( abs ( eta * a(kk) ) <= abs ( w ) ) then

      if ( w < 0.0D+00 ) then
        ifault = 2
        return
      end if

      u(k) = sqrt ( w )

    else

      u(k) = 0.0D+00
      nullty = nullty + 1

    end if

    j = j + icol
    det = det * u(k) * u(k)

  end do

  return
end

