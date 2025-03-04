subroutine newton_rc ( ido, n, x, fx )

!*****************************************************************************80
!
!! newton_rc() solves a small system of nonlinear equations.
!
!  Discussion:
!
!    NEWTON_RC uses Newton's method, with the jacobian approximated via
!    finite differences.
!
!    NEWTON_RC uses "reverse communication".  That is, NEWTON_RC does not call 
!    any user subroutines.  Instead, the user repeatedly calls NEWTON_RC, 
!    and the user and NEWTON_RC communicate via values of the parameter IDO.
!
!    To begin using NEWTON_RC, set IDO=0, set X to an approximate root,
!    set FX(1) to a small positive tolerance, and call NEWTON_RC.
!
!    NEWTON_RC will return with IDO=2.  Evaluate the residual at the
!    point X that is returned by NEWTON_RC, store that value in FX,
!    and call NEWTON_RC back with IDO=2.
!
!    NEWTON_RC will return with IDO=1.  Again, evaluate the residual
!    at the point X returned by NEWTON_RC, store that value in FX<
!    and call NEWTON_RC back with IDO=1.
!
!    This process will be repeated until NEWTON_RC returns with
!    IDO=0, meaning that the value X is a good approximation to
!    the root, or IDO=-1, in which case the algorithm failed.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer IDO, indicates the user's action.
!    * 0, this is a new problem, with N, X, and FX supplied.
!    * 1 or 2, FX evaluated at the X returned by previous call.
!
!    integer N, the number of equations.
!
!    real ( kind = rk ) X(N), the current iterate.
!
!    real ( kind = rk ) FX(N): the function value at the current iterate.
!
!  Output:
!
!    integer IDO.
!    * 0, the output value X is a good approximation to the solution.
!    * 1, please evaluate F(X).  (Jacobian is being approximated).
!    * 2, please evaluate F(X).
!    * -1, the algorithm failed.
!
!   real ( kind = rk ) X(N), the current point being considered.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: maxn = 20

  integer n

  real ( kind = rk ) delx(n)
  real ( kind = rk ), save :: delxj
  real ( kind = rk ), save :: eps
  real ( kind = rk ) fx(n)
  integer, save :: fxcall
  real ( kind = rk ), save :: fprime(maxn*maxn)
  real ( kind = rk ), save :: fxnrm
  real ( kind = rk ), save :: fxnrm_previous
  real ( kind = rk ), save :: fxnrm0
  real ( kind = rk ), save :: fxold(maxn)
  integer ido
  integer ihi
  integer ilo
  integer info
  integer, save :: j
  integer, parameter :: maxcall = 500
  integer, save :: ncall
  logical :: verbose = .false.
  real ( kind = rk ) x(n)
  real ( kind = rk ) xnrm
  real ( kind = rk ), save :: xold

  if ( maxn < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'newton_rc - Fatal error!'
    write ( *, '(a,i6,a,i6)' ) '  Input value n = ', n, ' exceeds limit maxn = ', maxn
    ido = -1
    return
  end if

  ncall = ncall + 1

  if ( maxcall < ncall ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'newton_rc - Fatal error!'
    write ( *, '(a,i6,a)' ) '  newton_rc called ', ncall, ' times.'
    write ( *, '(a,i6)' ) '  Internal limit is ', maxcall
    ido = -1
    return
  end if
!
!  IDO=0, user is starting a new problem.
!
  if ( ido == 0 ) then

    eps = sqrt ( epsilon ( eps ) )
    fxcall = 0
    fxnrm = sqrt ( sum ( fx(1:n)**2 ) )
    fxnrm0 = fxnrm
    ncall = 0
    ido = 2
!
!  IDO=1, user is returning F(X+delX(J)) for a jacobian estimation
!
  else if ( ido == 1 ) then

    ilo = 1 + ( j - 1 ) * n
    ihi = n + ( j - 1 ) * n
    fprime(ilo:ihi) = ( fx(1:n) - fxold(1:n) ) / delxj

    x(j) = xold

    if ( j < n ) then

      j = j + 1
      delxj = eps * ( abs ( x(j) ) + 1.0D+00 )
      xold = x(j)
      x(j) = x(j) + delxj

    else

      delx(1:n) = - fx(1:n)

      call r8mat_fs ( n, fprime, delx, info )

      if ( info /= 0 )then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'newton_rc - Fatal error!'
        write ( *, '(a)' ) '  Linear system solution failed.'
        ido = -1
        return
      end if

      x(1:n) = x(1:n) + delx(1:n)

      ido = 2

    end if
!
!  IDO = 2, user is returning F(X) for a convergence evaluation.
!
  else if ( ido == 2 ) then

    fxcall = fxcall + 1

    fxnrm_previous = fxnrm

    xnrm = sqrt ( sum ( x(1:n)**2 ) )
    fxnrm = sqrt ( sum ( fx(1:n)**2 ) )

    if ( verbose ) then
      write ( *, '(a,g14.6,a,g14.6)' ) '  ||X|| = ', xnrm, '  ||FX|| = ', fxnrm
    end if

    if ( fxnrm <= eps * ( fxnrm0 + 1.0 ) ) then
      ido = 0
      return
    end if

    if ( 15 < fxcall ) then
      if ( 0.95D+00 * fxnrm_previous < fxnrm ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'newton_rc - Warning.'
        write ( *, '(a,g14.6)' ) '  Previous ||FX|| = ', fxnrm_previous
        write ( *, '(a,g14.6)' ) '  Latest   ||FX|| = ', fxnrm
        write ( *, '(a,i6,a)' ) '  Convergence has slowed down after ', fxcall, ' steps'
        ido = -1
        return
      end if
    end if

    ido = 1
    fxold(1:n) = fx(1:n)
    j = 1
    delxj = eps * ( abs ( x(j) ) + 1.0D+00 )
    xold = x(j)
    x(j) = x(j) + delxj
!
!  Unexpected value of IDO.
!
  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'newton_rc - Fatal error.'
    write ( *, '(a,i6)' ) '  Unexpected value of IDO = ', ido
    ido = -1

  end if

  return
end
subroutine r8mat_fs ( n, a, b, info )

!*****************************************************************************80
!
!! r8mat_fs() factors and solves a system with one right hand side.
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
!  Input:
!
!    integer N, the order of the matrix.
!    N must be positive.
!
!    real ( kind = rk ) A(N,N), the coefficient matrix of the linear system.
!
!    real ( kind = rk ) B(N), the right hand side of the linear system.
!
!  Output:
!
!    real ( kind = rk ) A(N,N), a unit upper triangular matrix,
!    the U factor of an LU factorization of the
!    original coefficient matrix.
!
!    real ( kind = rk ) B(N), the solution of the linear systems.
!
!    integer INFO, singularity flag.
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
      call r8mat_print ( n, n, a, '  The Jacobian matrix:' )
      info = -1
      return
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
!! r8mat_print_some() prints some of a real matrix.
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

        write ( ctemp(j2), '(g14.6)' ) a(i,j)

      end do

      write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end

