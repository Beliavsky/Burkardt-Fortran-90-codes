program main

!*****************************************************************************80
!
!! slap_test() tests slap().
!
!  Modified:
!
!    06 September 2022
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'slap_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test slap().'

  call slap_dgmres_test ( )
  call slap_quick_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'slap_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  call timestamp ( )

  return
end
subroutine slap_dgmres_test ( )

!*****************************************************************************80
!
!! slap_dgmres_test() tests slap().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 January 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: ligw = 20
  integer, parameter :: maxl = 20
  integer, parameter :: n = 100

  integer, parameter :: lrgw = 1 + n * ( maxl + 6 ) + maxl * ( maxl + 3 )

  real ( kind = rk ), allocatable, dimension ( : ) :: a
  real ( kind = rk ) b(n)
  real ( kind = rk ) err
  integer, allocatable, dimension ( : ) :: ia
  integer ierr
  integer igwk(ligw)
  integer isym
  integer iter
  integer itmax
  integer itol
  integer iunit
  integer iwork(1)
  integer, allocatable, dimension ( : ) :: ja
  external matvec_triad
  integer mode
  external msolve_identity
  integer nelt
  real ( kind = rk ) rgwk(lrgw)
  real ( kind = rk ) rwork(1)
  real ( kind = rk ) sb(n)
  real ( kind = rk ) sx(n)
  real ( kind = rk ) tol
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'slap_dgmres_test():'
  write ( *, '(a)' ) '  Test the slap library.'

  isym = 0
  itol = 0
  tol = 0.00001D+00
  itmax = 500
  iunit = 0
  sb(1:n) = 1.0D+00
  sx(1:n) = 1.0D+00

  igwk(1) = maxl
  igwk(2) = maxl
  igwk(3) = 0
  igwk(4) = 0
  igwk(5) = 60
!
!  Determine the amount of storage needed for the matrix.
!
  mode = 0
  call matset_triad ( n, mode, nelt, ia, ja, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Problem dimension N = ', n
  write ( *, '(a,i6)' ) '  Matrix storage NELT = ', nelt

  allocate ( ia(1:nelt) )
  allocate ( ja(1:nelt) )
  allocate (  a(1:nelt) )
!
!  Assign the values for the matrix.
!
  mode = 1
  call matset_triad ( n, mode, nelt, ia, ja, a )
!
!  Set the solution X.
!
  call r8vec_indicator ( n, x )

  call r8vec_print_some ( n, x, 10, '  Part of the exact solution X:' )
!
!  Set the right hand side B!
!
  call matvec_triad ( n, x, b, nelt, ia, ja, a, isym )

  call r8vec_print_some ( n, b, 10, '  Part of the right hand side B:' )
!
!  Zero out X, because it will be used as an initial guess.
!
  x(1:n) = 0.0D+00
!
!  Call DGMRES to solve the system..
!
  call dgmres ( n, b, x, nelt, ia, ja, a, isym, matvec_triad, &
    msolve_identity, itol, tol, itmax, iter, err, ierr, iunit, sb, &
    sx, rgwk, lrgw, igwk, ligw, rwork, iwork )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of iterations:  ', iter
  write ( *, '(a,g14.6)' ) '  Convergence measure is ', rgwk(1)
  write ( *, '(a,g14.6)' ) '  Error estimate ', err
  write ( *, '(a,i6)' ) '  Error code is ', ierr

  call r8vec_print_some ( n, x, 10, '  Part of the computed solution X:' )
!
!  Free memory.
!
  deallocate ( ia )
  deallocate ( ja )
  deallocate ( a )

  return
end
subroutine matset_triad ( n, mode, nelt, ia, ja, a )

!*****************************************************************************80
!
!! MATSET_TRIAD sets the data structure for a sparse matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements in the vectors.
!
!    Input, integer MODE.
!    0, setup mode.  Count the number of nonzero elements.
!    1, assignment mode.  NELT is set.  Set the matrix entries.
!
!    Input/output, integer NELT, the number of nonzero entries.
!    On call with MODE = 0, this is an output quantity.
!    On call with MODE = 1, this is an input quantity.
!
!    Output, integer IA(NELT), JA(NELT), real ( kind = rk ) A(NELT), the data
!    structure storing the sparse matrix.  These values are only
!    set on a call with MODE = 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nelt

  real ( kind = rk ) a(nelt)
  integer i
  integer ia(nelt)
  integer ja(nelt)
  integer mode
  integer k

  k = 0
  do i = 1, n

    if ( 1 < i ) then
      k = k + 1
      if ( mode == 1 ) then
        ia(k) = i
        ja(k) = i-1
        a(k) = -1.0D+00
      end if
    end if

    k = k + 1
    if ( mode == 1 ) then
      ia(k) = i
      ja(k) = i
      a(k) = 2.0D+00
    end if

    if ( i < n ) then
      k = k + 1
      if ( mode == 1 ) then
        ia(k) = i
        ja(k) = i+1
        a(k) = -1.0D+00
      end if
    end if

  end do

  if ( mode == 0 ) then
    nelt = k
  end if

  return
end
subroutine matvec_triad ( n, x, y, nelt, ia, ja, a, isym )

!*****************************************************************************80
!
!! MATVEC_TRIAD computes A*X for a matrix A stored in SLAP Triad form.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements in the vectors.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = rk ) Y(N), the product A * X.
!
!    Input, integer NELT, the number of nonzero entries in A.
!
!    Input, integer IA(NELT), JA(NELT), real ( kind = rk ) A(NELT), the data
!    structure storing the sparse matrix.
!
!    Input, integer ISYM, is 0 if all nonzero entries of the matrix
!    are stored, and 1 if only the diagonal and upper or lower triangle
!    are stored.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nelt

  real ( kind = rk ) a(nelt)
  integer ia(nelt)
  integer isym
  integer ja(nelt)
  integer k
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)

  y(1:n) = 0.0D+00

  do k = 1, nelt
    y(ia(k)) = y(ia(k)) + a(k) * x(ja(k))
  end do

  return
end
subroutine msolve_identity ( n, r, z, nelt, ia, ja, a, isym, rwork, iwork )

!*****************************************************************************80
!
!! MSOLVE_IDENTITY applies the identity matrix preconditioner.
!
!  Discussion:
!
!    Most SLAP solver routines require a preconditioner routine
!    that can solve M * Z = R.  If no preconditioning is required,
!    then you can simply behave as though the preconditioning matrix
!    M was the identity matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements in the vectors.
!
!    Input, real ( kind = rk ) R(N), the right hand side.
!
!    Output, real ( kind = rk ) Z(N), the solution of M * Z = R.
!
!    Input, integer NELT, the number of nonzero entries in A.
!
!    Input, integer IA(NELT), JA(NELT), real ( kind = rk ) A(NELT), 
!    the data structure storing the sparse matrix.
!
!    Input, integer ISYM, is 0 if all nonzero entries of the matrix
!    are stored, and 1 if only the diagonal and upper or lower triangle
!    are stored.
!
!    Input, real ( kind = rk ) RWORK(*), a real array that
!    can be used to pass information to the preconditioner.
!
!    Input, integer IWORK(*), an integer array that
!    can be used to pass information to the preconditioner.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nelt

  real ( kind = rk ) a(nelt)
  integer ia(nelt)
  integer isym
  integer iwork(*)
  integer ja(nelt)
  real ( kind = rk ) r(n)
  real ( kind = rk ) rwork(*)
  real ( kind = rk ) z(n)

  z(1:n) = r(1:n)

  return
end
subroutine r8vec_indicator ( n, a )

!*****************************************************************************80
!
!! R8VEC_INDICATOR sets a real vector to the indicator vector.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements of A.
!
!    Output, real ( kind = rk ) A(N), the array to be initialized.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  integer i

  do i = 1, n
    a(i) = real ( i, kind = rk )
  end do

  return
end
subroutine r8vec_print_some ( n, a, max_print, title )

!*****************************************************************************80
!
!! R8VEC_PRINT_SOME prints "some" of a real vector.
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vector, is no more than MAX_PRINT, then
!    the entire vector is printed, one entry per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries of the vector.
!
!    Input, real ( kind = rk ) A(N), the vector to be printed.
!
!    Input, integer MAX_PRINT, the maximum number of lines to print.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  integer i
  integer max_print
  character ( len = * ) title

  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  if ( n <= max_print ) then

    if ( all ( a(1:n) == aint ( a(1:n) ) ) ) then
      do i = 1, n
        write ( *, '(i6,2x,i6)' ) i, int ( a(i) )
      end do
    else if ( all ( abs ( a(1:n) ) < 1000000.0D+00 ) ) then
      do i = 1, n
        write ( *, '(i6,2x,f14.6)' ) i, a(i)
      end do
    else
      do i = 1, n
        write ( *, '(i6,2x,g14.6)' ) i, a(i)
      end do
    end if

  else if ( 3 <= max_print ) then

    if ( all ( a(1:max_print-2) == aint ( a(1:max_print-2) ) ) ) then
      do i = 1, max_print-2
        write ( *, '(i6,2x,i6)' ) i, int ( a(i) )
      end do
    else if ( all ( abs ( a(1:max_print-2) ) < 1000000.0D+00 ) ) then
      do i = 1, max_print-2
        write ( *, '(i6,2x,f14.6)' ) i, a(i)
      end do
    else
      do i = 1, max_print-2
        write ( *, '(i6,2x,g14.6)' ) i, a(i)
      end do
    end if

    write ( *, '(a)' ) '......  ..............'
    i = n

    if ( a(i) == aint ( a(i) ) ) then
      write ( *, '(i6,2x,i6)' ) i, int ( a(i) )
    else if (  abs ( a(i) ) < 1000000.0D+00 ) then
      write ( *, '(i6,2x,f14.6)' ) i, a(i)
    else
      write ( *, '(i6,2x,g14.6)' ) i, a(i)
    end if

  else

    if ( all ( a(1:max_print-1) == aint ( a(1:max_print-1) ) ) ) then
      do i = 1, max_print-1
        write ( *, '(i6,2x,i6)' ) i, int ( a(i) )
      end do
    else if ( all ( abs ( a(1:max_print-1) ) < 1000000.0D+00 ) ) then
      do i = 1, max_print-1
        write ( *, '(i6,2x,f14.6)' ) i, a(i)
      end do
    else
      do i = 1, max_print-1
        write ( *, '(i6,2x,g14.6)' ) i, a(i)
      end do
    end if

    i = max_print

    if ( a(i) == aint ( a(i) ) ) then
      write ( *, '(i6,2x,i6,a)' ) i, int ( a(i) ), '...more entries...'
    else if (  abs ( a(i) ) < 1000000.0D+00 ) then
      write ( *, '(i6,2x,f14.6,a)' ) i, a(i), '...more entries...'
    else
      write ( *, '(i6,2x,g14.6,a)' ) i, a(i), '...more entries...'
    end if

  end if

  return
end
subroutine slap_quick_test ( )

!*****************************************************************************80
!
!! slap_quick_test() tests slap().
!
!  Discussion:
!
!    This is a SLATEC Quick Checks program to test the *SLAP*
!    Version 2.0 package.  It generates a "random" matrix (See
!    DRMGEN) and then runs all the various methods with all the
!    various preconditoners and all the various stop tests.
!
!  Modified:
!
!    08 August 2006
!
!  Author:
!
!    Mark Seager
!    Lawrence Livermore National Laboratory
!
!  Local Parameters:
!
!    Local, integer KPRINT, determines the amount of output.
!    0  Quick checks - No printing.
!       Driver       - Short pass or fail message printed.
!    1  Quick checks - No message printed for passed tests,
!                      short message printed for failed tests.
!       Driver       - Short pass or fail message printed.
!    2  Quick checks - Print short message for passed tests,
!                      fuller information for failed tests.
!       Driver       - Pass or fail message printed.
!    3  Quick checks - Print complete quick check results.
!       Driver       - Pass or fail message printed.
!    4  Quick checks - Print complete quick check results.
!                      Prints matricies, etc.  Very verbose.
!       Driver       - Pass or fail message printed.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: maxiw = 50000
  integer, parameter :: maxn = 441
  integer, parameter :: maxrw = 50000
  integer, parameter :: mxnelt = 50000

  real ( kind = rk ) a(mxnelt)
  real ( kind = rk ) dens
  real ( kind = rk ) err
  real ( kind = rk ) f(maxn)
  real ( kind = rk ) factor
  integer i1mach
  integer ia(mxnelt)
  integer ierr
  integer iout
  integer istdi
  integer istdo
  integer isym
  integer itdi
  integer iter
  integer itmax
  integer itol
  integer itolgm
  integer iunit
  integer iwork(maxiw)
  integer ja(mxnelt)
  integer k
  integer kase
  integer, parameter :: kprint = 2
  integer leniw
  integer lenw
  integer lun
  character ( len = 72 ) mesg
  integer n
  integer nelt
  integer neltmx
  integer nfail
  integer nmax
  integer nsave
  real ( kind = rk ) rwork(maxrw)
  real ( kind = rk ) soln
  real ( kind = rk ) tol
  real ( kind = rk ) xiter(maxn)
!
!  Supplying the exact solution in common block SOLBLK allows slap,
!  when ITOL = 11, to compare its computed answer to the exact solution.
!
  common /solblk/ soln(maxn)

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'slap_quick_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the slap library.'
!
!  Set up the error routines.
!
  istdi = i1mach(1)
  istdo = i1mach(2)
  nfail = 0

  call xsetun ( lun )

  if ( kprint <= 1 ) then
    call xsetf ( 0 )
  else
    call xsetf ( 1 )
  end if

  call xermax ( 1000 )
!
!  Set the maximum problem sizes.
!
  neltmx = mxnelt
  nmax   = maxn
  leniw  = maxiw
  lenw   = maxrw
!
!  Set some input data.
!
  n      = nmax
  itmax  = n
  iout   = kprint
  factor = 1.2D+00

  if ( iout < 3 ) then
    iunit = 0
  else
    iunit = istdo
  end if
!
!  Set the error tolerance to depend on the machine epsilon.
!
  tol = max ( 1.0D+03 * epsilon ( tol ), 1.0D-06 )
!
!  Test routines using various convergence criteria.
!
  do kase = 3, 3

    if ( kase  ==  1 .or. kase  ==  2 ) then
      itol = kase
    else if ( kase  ==  3 ) then
      itol = 11
    end if
!
!  Test routines using nonsymmetric (ISYM=0) and symmetric
!  storage (ISYM=1).  For ISYM=0 a really non-symmetric matrix
!  is generated.  The amount of non-symmetry is controlled by
!  user.
!
     do isym = 0, 1
!
!  Set up a random matrix.
!
        call drmgen ( neltmx, factor, ierr, n, nelt, &
          isym, ia, ja, a, f, soln, rwork, iwork, iwork(n+1) )

        if ( ierr /= 0 ) then
           mesg = 'slapQC -- Fatal error (i1) generating '// &
                '*random* matrix.'
           call xerrwv( mesg,len(mesg),ierr,2,1,ierr,0, &
                0,0.0,0.0 )
        end if

        if ( isym == 0 ) then
           dens = real ( nelt, kind = rk ) / real ( n*n, kind = rk )
        else
           dens = real ( 2*nelt, kind = rk ) / real ( n*n, kind = rk )
        end if

        if ( 2 <= iout ) then
          write(istdo,1020) n, nelt, dens
          write(istdo,1030) tol
        end if
!
!  Convert to the SLAP-Column format and
!  write out matrix in SLAP-Column format, if desired.
!
        call ds2y( n, nelt, ia, ja, a, isym )

        if ( 4 <= iout ) then
           write(istdo,1040) (k,ia(k),ja(k),a(k),k=1,nelt)
           call dcpplt( n, nelt, ia, ja, a, isym, istdo )
        end if
!
!  BEGIN SLAP QUICK TESTS
!
!  DSJAC.
!
        if ( 3 <= iout ) then
          write ( *, '(a)' ) ' '
          write ( *, '(2x,a6,a, i2,a,i1 )' ) &
            'DSJAC :  ITOL = ', itol, '  ISYM = ', isym
        end if

        xiter(1:n) = 0.0D+00

        call dsjac(n, f, xiter, nelt, ia, ja, a, isym, &
             itol, tol, 2*itmax, iter, err, ierr, iunit, &
             rwork, lenw, iwork, leniw )

        call duterr( 'dsjac ',ierr,iout,nfail,istdo,iter,err )
!
!  DSGS.
!
        if ( 3 <= iout ) then
          write ( *, '(a)' ) ' '
          write ( *, '(2x,a6,a, i2,a,i1 )' ) &
            'DSOS  :  ITOL = ', itol, '  ISYM = ', isym
        end if

        xiter(1:n) = 0.0D+00

        call dsgs(n, f, xiter, nelt, ia, ja, a, isym, &
             itol, tol, itmax, iter, err, ierr, iunit, &
             rwork, lenw, iwork, leniw )

        call duterr( 'dsgs  ',ierr,iout,nfail,istdo,iter,err )
!
!  DSILUR.
!
        if ( 3 <= iout ) then
          write ( *, '(a)' ) ' '
          write ( *, '(2x,a6,a, i2,a,i1 )' ) &
            'DSILUR:  ITOL = ', itol, '  ISYM = ', isym
        end if

        xiter(1:n) = 0.0D+00

        call dsilur(n, f, xiter, nelt, ia, ja, a, isym, &
             itol, tol, itmax, iter, err, ierr, iunit, &
             rwork, lenw, iwork, leniw )

        call duterr( 'dsilur',ierr,iout,nfail,istdo,iter,err )
!
!  DSDCG.
!
        if ( isym == 1 ) then

           if ( 3 <= iout ) then
              write ( *, '(a)' ) ' '
              write ( *, '(2x,a6,a, i2,a,i1 )' ) &
                'DSDCG :  ITOL = ', itol, '  ISYM = ', isym
           end if

           xiter(1:n) = 0.0D+00

           call dsdcg(n, f, xiter, nelt, ia, ja, a, isym, &
                itol, tol, itmax, iter, err, ierr, iunit, &
                rwork, lenw, iwork, leniw )

           call duterr( 'dsdcg ',ierr,iout,nfail,istdo,iter,err )
        end if
!
!  DSICCG.
!
        if ( isym == 1 ) then

           if ( 3 <= iout ) then
             write ( *, '(a)' ) ' '
             write ( *, '(2x,a6,a, i2,a,i1 )' ) &
               'DSICG :  ITOL = ', itol, '  ISYM = ', isym
           end if

           xiter(1:n) = 0.0D+00

           call dsiccg(n, f, xiter, nelt, ia, ja, a, isym, &
                itol, tol, itmax, iter, err, ierr, iunit, rwork, &
                lenw, iwork, leniw )

           call duterr( 'dsiccg',ierr,iout,nfail,istdo,iter,err )

        end if
!
!  DSDCGN.
!
        if ( 3 <= iout ) then
          write ( *, '(a)' ) ' '
          write ( *, '(2x,a6,a, i2,a,i1 )' ) &
            'DSDCGN:  ITOL = ', itol, '  ISYM = ', isym
        end if

        xiter(1:n) = 0.0D+00

        call dsdcgn(n, f, xiter, nelt, ia, ja, a, isym, itol, &
             tol, itmax, iter, err, ierr, iunit, rwork, lenw, &
             iwork, leniw )

        call duterr( 'dsdcgn',ierr,iout,nfail,istdo,iter,err )
!
!  DSLUCN.
!
        if ( 3 <= iout ) then
          write ( *, '(a)' ) ' '
          write ( *, '(2x,a6,a, i2,a,i1 )' ) &
            'DSLUCN:  ITOL = ', itol, '  ISYM = ', isym
        end if

        xiter(1:n) = 0.0D+00

        call dslucn(n, f, xiter, nelt, ia, ja, a, isym, itol, &
             tol, itmax, iter, err, ierr, iunit, rwork, lenw, &
             iwork, leniw )

        call duterr( 'dslucn',ierr,iout,nfail,istdo,iter,err )
!
!  DSDBCG.
!
        if ( 3 <= iout ) then
          write ( *, '(a)' ) ' '
          write ( *, '(2x,a6,a, i2,a,i1 )' ) &
            'dsdbcg:  itol = ', itol, '  isym = ', isym
        end if

        xiter(1:n) = 0.0d+00

        call dsdbcg(n, f, xiter, nelt, ia, ja, a, isym, itol, &
             tol, itmax, iter, err, ierr, iunit, rwork, lenw, &
             iwork, leniw )

        call duterr( 'dsdbcg',ierr,iout,nfail,istdo,iter,err )
!
!  DSLUBC.
!
        if ( 3 <= iout ) then
          write ( *, '(a)' ) ' '
          write ( *, '(2x,a6,a, i2,a,i1 )' ) &
            'dslubc:  itol = ', itol, '  isym = ', isym
        end if

        xiter(1:n) = 0.0d+00

        call dslubc(n, f, xiter, nelt, ia, ja, a, isym, &
             itol, tol, itmax, iter, err, ierr, iunit, &
             rwork, lenw, iwork, leniw )

        call duterr( 'dslubc',ierr,iout,nfail,istdo,iter,err )
!
!  DSDCGS.
!
        if ( 3 <= iout ) then
          write ( *, '(a)' ) ' '
          write ( *, '(2x,a6,a, i2,a,i1 )' ) &
            'dsdcos:  itol = ', itol, '  isym = ', isym
        end if

        xiter(1:n) = 0.0d+00

        call dsdcgs(n, f, xiter, nelt, ia, ja, a, isym, itol, &
             tol, itmax, iter, err, ierr, iunit, rwork, lenw, &
             iwork, leniw )

        call duterr( 'dsdcgs',ierr,iout,nfail,istdo,iter,err )
!
!  DSLUCS.
!
        if ( 3 <= iout ) then
          write ( *, '(a)' ) ' '
          write ( *, '(2x,a6,a, i2,a,i1 )' ) &
            'dslucs:  itol = ', itol, '  isym = ', isym
        end if

        xiter(1:n) = 0.0d+00

        call dslucs(n, f, xiter, nelt, ia, ja, a, isym, &
             itol, tol, itmax, iter, err, ierr, iunit, &
             rwork, lenw, iwork, leniw )

        call duterr( 'dslucs',ierr,iout,nfail,istdo,iter,err )
!
!  DSDOMN.
!
!VD$ NOVECTOR

        do nsave = 0, 3

           if ( 3 <= iout ) then
             write ( *, '(a)' ) ' '
             write ( *, '(2x,a6,a, i2,a,i1,a,i2 )' ) &
               'dsdomn:  itol = ', itol, '  isym = ', isym, &
               '  nsave = ', nsave
           end if

           xiter(1:n) = 0.0d+00

           call dsdomn(n, f, xiter, nelt, ia, ja, a, &
                isym, nsave, itol, tol, itmax, iter, err, ierr, &
                iunit, rwork, lenw, iwork, leniw )

           call duterr( 'dsdomn',ierr,iout,nfail,istdo,iter,err )

        end do
!
!  DSLUOM
!
!VD$ NOVECTOR

        do nsave = 0, 3

           if ( 3 <= iout ) then
             write ( *, '(a)' ) ' '
             write ( *, '(2x,a6,a, i2,a,i1,a,i2 )' ) &
               'dsluom:  itol = ', itol, '  isym = ', isym, &
               '  nsave = ', nsave
           end if

           xiter(1:n) = 0.0d+00

           call dsluom(n, f, xiter, nelt, ia, ja, a, &
                isym, nsave, itol, tol, itmax, iter, err, ierr, &
                iunit, rwork, lenw, iwork, leniw )

           call duterr( 'dsluom',ierr,iout,nfail,istdo,iter,err )

        end do
!
!  DSDGMR
!
!VD$ NOVECTOR

        do nsave = 5, 12

           if ( 3 <= iout ) then
             write ( *, '(a)' ) ' '
             write ( *, '(2x,a6,a, i2,a,i1,a,i2 )' ) &
               'dsdgmr:  itol = ', itol, '  isym = ', isym, &
               '  nsave = ', nsave
           end if

           xiter(1:n) = 0.0d+00
           itolgm = 0

           call dsdgmr(n, f, xiter, nelt, ia, ja, a, &
                isym, nsave, itolgm, tol, itmax, iter, err, ierr, &
                iunit, rwork, lenw, iwork, leniw )

           call duterr( 'dsdgmr',ierr,iout,nfail,istdo,iter,err )

        end do
!
!  DSLUGM
!
!VD$ NOVECTOR

        do nsave = 5, 12

           if ( 3 <= iout ) then
             write ( *, '(a)' ) ' '
             write ( *, '(2x,a6,a, i2,a,i1,a,i2 )' ) &
               'dslugm:  itol = ', itol, '  isym = ', isym, &
               '  nsave = ', nsave
           end if

           xiter(1:n) = 0.0d+00

           call dslugm(n, f, xiter, nelt, ia, ja, a, &
                isym, nsave, itol, tol, itmax, iter, err, ierr, &
                iunit, rwork, lenw, iwork, leniw )

           call duterr( 'dslugm',ierr,iout,nfail,istdo,iter,err )

        end do

     end do

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '

  if ( NFAIL == 0 ) then
    write ( *, '(a)' ) '*******************************************************'
    write ( *, '(a)' ) '**** All SLAP Double Precision Quick Checks Passed ****'
    write ( *, '(a)' ) '****                 No Errors                     ****'
    write ( *, '(a)' ) '*******************************************************'
  else
     write(istdo,1060) nfail
  end if

  return

 1020 FORMAT(/'                * RANDOM Matrix of size',I5,'*' &
       /'                ', &
       'Number of non-zeros & Density = ', I5,E16.7)
 1030 FORMAT('                Error tolerance = ',E16.7)
 1040 FORMAT(/'  ***** SLAP Column Matrix *****'/ &
          ' Indx   ia   ja     a'/(1X,I4,1X,I4,1X,I4,1X,E16.7))

 1060 FORMAT(// &
       '************************************************'/ &
       '**     ===>',I3,' Failures detected <===      **'/ &
       '**     SLAP Double Precision Quick Checks     **'/ &
       '** Set KPRINT = 3 for DEBUG information and   **'/ &
       '** rerun the tests to determine the problem.  **'/ &
       '************************************************')
end
subroutine DUTERR ( METHOD, IERR, IOUT, NFAIL, ISTDO, ITER, ERR )

!*****************************************************************************80
!
!! DUTERR outputs error messages for the SLAP quick checks.
!
!  Modified:
!
!    08 August 2006
!
!  Author:
!
!    Mark Seager,
!    Lawrence Livermore National Laboratory
!
!  Parameters:
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) err
  integer ierr
  integer iout
  integer istdo
  integer iter
  character*6 method
  integer nfail

  if ( ierr /= 0 ) then
    nfail = nfail+1
  end if

  if ( iout == 1 .and. ierr /= 0 ) then
     write(istdo,1000) method
  end if

  if ( iout == 2 ) then
     if ( ierr == 0 ) then
        write(istdo,1010) method
     else
        write(istdo,1020) method,ierr,iter,err
     end if
  end if

  if ( 3 <= iout ) then
     if ( ierr == 0 ) then
        write(istdo,1030) method,ierr,iter,err
     else
        write(istdo,1020) method,ierr,iter,err
     end if
  end if

  return
 1000 FORMAT( 1X,A6,' : **** FAILURE ****')
 1010 FORMAT( 1X,A6,' : **** PASSED  ****')
 1020 FORMAT(' **************** WARNING ***********************'/ &
         ' **** ',A6,' Quick Test FAILED: IERR = ',I5,' ****'/ &
         ' **************** WARNING ***********************'/ &
         ' Iteration Count = ',I3,' Stop Test = ',E12.6)
 1030 FORMAT(' ***************** PASSED ***********************'/ &
         ' **** ',A6,' Quick Test PASSED: IERR = ',I5,' ****'/ &
         ' ***************** PASSED ***********************'/ &
         ' Iteration Count = ',I3,' Stop Test = ',E12.6)
end
subroutine drmgen ( neltmx, factor, ierr, n, nelt, isym, &
  ia, ja, a, f, soln, dsum, itmp, idiag )

!*****************************************************************************80
!
!! DRMGEN generates a random matrix for SLAP quick checks.
!
!  Discussion:
!
!    The matrix is generated by choosing a random number of
!    entries for each column and then chosing negative random
!    numbers for each off diagionals.   The diagionals elements
!    are chosen to be positive and large enough so the matrix
!    is slightly diagionally dominant.  The lower triangle of
!    the matrix is generated and if isym == 0 (all matrix elements
!    stored) the upper triangle elements are chosen so that they
!    are FACTOR times the coresponding lower triangular element.
!
!  Modified:
!
!    08 August 2006
!
!  Author:
!
!    Mark Seager,
!    Lawrence Livermore National Laboratory
!
!  Parameters:
!
!    NELTMX :IN       Integer.
!    Maximum number of non-zeros that can be created by this
!    routine for storage in the IA, JA, A arrays,  see below.
!
!    FACTOR :IN       Double Precision.
!    Non-zeros in the upper triangle are set to FACTOR times
!    the coresponding entry in the lower triangle when a non-
!    symmetric matrix is requested (See ISYM, below).
!
!    IERR   :OUT      Integer.
!    Return error flag.
!    0 => everything went OK.
!    1 => Ran out of space trying to create matrix.
!    Set NELTMX to something larger and retry.
!
!    N      :IN       Integer.
!    Size of the linear system to generate (number of unknowns).
!
!    NELT   :OUT      Integer.
!    Number of non-zeros stored in the IA, JA, A arrays, see below.
!
!    ISYM   :IN       Integer.
!    Flag to indicate the type of matrix to generate:
!    0 => Non-Symmetric Matrix (See FACTOR, above).
!    1 => Symmetric Matrix.
!
!    IA     :OUT      Integer IA(NELTMX).
!    Stores the row indicies for the non-zeros.
!
!    JA     :OUT      Integer JA(NELTMX).
!    Stores the column indicies for the non-zeros.
!
!    A      :OUT      Double Precision A(NELTMX).
!    Stores the values of the non-zeros.
!
!    F      :OUT      Double Precision F(N).
!    The right hand side of the linear system.  Obtained by mult-
!    iplying the matrix time SOLN, see below.
!
!    SOLN   :OUT      Double Precision SOLN(N).
!    The true solution to the linear system.  Each component is
!    chosen at random (0.0<SOLN(I)<1.0, I=1,N)
!
!    DSUM   :WORK     Double Precision DSUM(N).
!
!    ITMP   :WORK     Integer ITMP(N).
!
!    IDIAG  :WORK     Integer IDIAG(N).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer neltmx

  real ( kind = rk ) a(neltmx)
  real ( kind = rk ) dsum(n)
  real ( kind = rk ) f(n)
  real ( kind = rk ) factor
  integer i
  integer ia(neltmx)
  integer icol
  integer idiag(n)
  integer ierr
  integer inum
  integer irow
  integer iseed
  integer isym
  integer itmp(n)
  integer ja(neltmx)
  integer k
  integer nelt
  integer nl
  real rand
  real rn
  real ( kind = rk ) soln(n)
!
!  Start by setting the random number generator seed.
!  This is done for reproducability in debugging.  Remove
!  the seed setting call for production testing.
!
  ierr = 0

  idiag(1:n) = 0
  dsum(1:n) = -1.0d+00
!
!  Set the matrix elements.
!  Loop over the columns.
!
  nelt = 0

!VD$ NOCONCUR

  do icol = 1, n

     nl = n+1-icol
!
!  To keep things sparse divide by two, three or four or ...
!
     call random_number ( harvest = rn )

     inum = ( int ( rn * nl ) + 1)/3

     call dmpl ( nl, inum, itmp )
!
!  Set up this column (and row, if non-symmetric structure).
!
!VD$ NOVECTOR
!VD$ NOCONCUR

     do irow = 1, inum

        nelt = nelt + 1

        if ( nelt > neltmx ) then
           ierr = 1
           return
        end if

        ia(nelt) = n+1-itmp(irow)
        ja(nelt) = icol

        if ( ia(nelt) == icol ) then

           idiag(icol) = nelt
        else

           call random_number ( harvest = rn )
           a(nelt) = - rn
           dsum(icol) = dsum(icol) + a(nelt)
           if ( isym == 0 ) then
!
!  Copy this element into upper triangle.
!
              nelt = nelt + 1

              if ( nelt > neltmx ) then
                 ierr = 1
                 return
              end if

              ia(nelt) = icol
              ja(nelt) = ia(nelt-1)
              a(nelt)  = a(nelt-1)*factor
              dsum(ja(nelt)) = dsum(ja(nelt)) + a(nelt)
           else
              dsum(ia(nelt)) = dsum(ia(nelt)) + a(nelt)
           end if

        end if
     end do

     if ( IDIAG(ICOL) == 0 ) then
!
!  Add a diagional to the column.
!
        nelt = nelt + 1
        if ( nelt > neltmx ) then
           ierr = 1
           return
        end if

        idiag(icol) = nelt
        a(nelt) = 0.0d0
        ia(nelt) = icol
        ja(nelt) = icol

     end if

  end do
!
!  Clean up the diagionals.
!
!VD$ NODEPCHK
!LLL. OPTION ASSERT (NOHAZARD)
!DIR$ IVDEP

  do i = 1, n
     a(idiag(i)) = -1.0001D+00 * dsum(i)
  end do
!
!  Set a random solution and determine the right-hand side.
!
!VD$ NOVECTOR
!VD$ NOCONCUR

  call random_number ( harvest = soln(1:n) )

  f(1:n) = 0.0d+00
!
!VD$ NOVECTOR
!VD$ NOCONCUR

  do k = 1, nelt
     f(ia(k)) = f(ia(k)) + a(k)*soln(ja(k))
     if ( isym /= 0 .and. ia(k) /= ja(k) ) then
        f(ja(k)) = f(ja(k)) + a(k)*soln(ia(k))
     end if
  end do

  return
end
subroutine dmpl ( n, m, indx )

!*****************************************************************************80
!
!! DMPL picks M distinct integers at random between 1 and N.
!
!  Modified:
!
!    08 August 2006
!
!  Author:
!
!    Mark Seager
!
!  Parameters:
!
!    Input, integer N, the upper limit on the values.
!
!    Input, integer M, the number of values to pick.
!
!    Output, integer INDX(M), distinct integers between 1 and N.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m

  logical found
  integer i
  integer id
  integer indx(m)
  integer j
  integer n
  real rn
!
!  Check the input.
!
  if ( n * m < 0 .or. n < m ) then
    return
  end if
!
!  Set the indices.
!
  call random_number ( harvest = rn )

  indx(1) = int ( rn * real ( n, kind = rk ) ) + 1

!VD$ NOCONCUR

  do i = 2, m

    do

      call random_number ( harvest = rn )
      id = int ( rn * real ( n, kind = rk ) ) + 1
!
!  Check to see if id has already been chosen.
!
!VD$ NOVECTOR
!VD$ NOCONCUR

      found = .true.

      do j = 1, i-1

        if ( id == indx(j) ) then
          found = .false.
          exit
        end if

      end do

      if ( found ) then
        indx(i) = id
        exit
      end if

    end do

  end do

  return
end
