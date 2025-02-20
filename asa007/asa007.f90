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
!    27 August 2021
!
!  Author:
!
!    Original FORTRAN77 version by Michael Healy.
!    Modifications by AJ Miller.
!    FORTRAN90 version by John Burkardt.
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
!    integer NN, the dimension of A, (N*(N+1))/2.
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
!    3, if NN < (N*(N+1))/2.
!
!  Local:
!
!    real ( kind = rk ) ETA, should be set equal to the smallest positive
!    value such that 1.0 + ETA is calculated as being greater than 1.0 in the
!    accuracy being used.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

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
  integer n
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
subroutine syminv ( a, n, c, w, nullty, ifault )

!*****************************************************************************80
!
!! syminv() computes the inverse of a symmetric matrix.
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
!    Original FORTRAN77 version by Michael Healy
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Michael Healy,
!    Algorithm AS 7:
!    Inversion of a Positive Semi-Definite Symmetric Matrix,
!    Applied Statistics,
!    Volume 17, Number 2, 1968, pages 198-199.
!
!  Input:
!
!    real ( kind = rk ) A((N*(N+1))/2), a symmetric positive definite 
!    matrix stored by rows in lower triangular form as a one dimensional 
!    array, in the sequence
!    A(1,1),
!    A(2,1), A(2,2),
!    A(3,1), A(3,2), A(3,3), and so on.
!
!    integer N, the order of A.
!
!  Output:
!
!    real ( kind = rk ) C((N*(N+1))/2), the inverse of A, or generalized
!    inverse if A is singular, stored using the same storage scheme employed
!    for A.  The program is written in such a way that A and U can share 
!    storage.
!
!    integer NULLTY, the rank deficiency of A.  If NULLTY 
!    is zero, the matrix is judged to have full rank.
!
!    integer IFAULT, error indicator.
!    0, no error detected.
!    1, N < 1.
!    2, A is not positive semi-definite.
!
!  Workspace:
!
!    real ( kind = rk ) W(N).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a((n*(n+1))/2)
  real ( kind = rk ) c((n*(n+1))/2)
  integer i
  integer icol
  integer ifault
  integer irow
  integer j
  integer jcol
  integer k
  integer l
  integer mdiag
  integer ndiag
  integer nn
  integer nrow
  integer nullty
  real ( kind = rk ) w(n)
  real ( kind = rk ) x

  ifault = 0

  if ( n <= 0 ) then
    ifault = 1
    return
  end if

  nrow = n
!
!  Compute the Cholesky factorization of A.
!  The result is stored in C.
!
  nn = ( n * ( n + 1 ) ) / 2

  call cholesky ( a, n, nn, c, nullty, ifault )

  if ( ifault /= 0 ) then
    return
  end if
!
!  Invert C and form the product (Cinv)' * Cinv, where Cinv is the inverse
!  of C, row by row starting with the last row.
!  IROW = the row number,
!  NDIAG = location of last element in the row.
!
  irow = nrow
  ndiag = nn

  do
!
!  Special case, zero diagonal element.
!
    if ( c(ndiag) == 0.0D+00 ) then

      l = ndiag
      do j = irow, nrow
        c(l) = 0.0D+00
        l = l + j
      end do

    else

      l = ndiag
      do i = irow, nrow
        w(i) = c(l)
        l = l + i
      end do

      icol = nrow
      jcol = nn
      mdiag = nn

      do

        l = jcol

        if ( icol == irow ) then
          x = 1.0D+00 / w(irow)
        else
          x = 0.0D+00
        end if

        k = nrow

        do while ( irow < k )

          x = x - w(k) * c(l)
          k = k - 1
          l = l - 1

          if ( mdiag < l ) then
            l = l - k + 1
          end if

        end do

        c(l) = x / w(irow)

        if ( icol <= irow ) then
          exit
        end if

        mdiag = mdiag - icol
        icol = icol - 1
        jcol = jcol - 1

      end do

    end if

    ndiag = ndiag - irow
    irow = irow - 1

    if ( irow <= 0 ) then
      exit
    end if

  end do

  return
end

