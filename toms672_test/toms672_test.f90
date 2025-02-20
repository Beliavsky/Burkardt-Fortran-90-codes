program main

!*****************************************************************************80
!
!! toms672_test() tests toms672().
!
!  Discussion:
!
!    The following templates demonstrate the use of the procedure EXTEND
!    to generate sequences of extended quadrature rules for various weight
!    functions given the definitions of the orthogonal polynomial 3-term
!    recurrence relations.
!
!    ICASE selects the required initial sequence as follows:
!    * 1, 3 point Gauss-Legendre in [-1,1]
!    * 2, 2 point Gauss-Lobatto in [-1,1]
!    * 3, 6 point Radau in [-1,1]
!    * 4, 2 point Gauss-Laguerre in [0,+oo)
!    * 5, 3 point Gauss-Hermite in (-oo,+oo)
!    * 6, 3 point Gauss-Jacobi in [0,1]
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 March 2009
!
!  Author:
!
!    Original FORTRAN77 version by Thomas Patterson.
!    FORTRAN90 version by John Burkardt.
!
  implicit none

  integer icase
  integer nseq

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'toms672_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test toms672().'
!
!  Get data from terminal.
!
!  Select demonstration.
!
  do

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Choose ICASE, the initial rule to extend:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  1, 3 point Gauss-Legendre in [-1,1];'
    write ( *, '(a)' ) '  2, 2 point Gauss-Lobatto in [-1,1];'
    write ( *, '(a)' ) '  3, 6 point Radau in [-1,1];'
    write ( *, '(a)' ) '  4, 2 point Gauss-Laguerre in [0,+oo);'
    write ( *, '(a)' ) '  5, 3 point Gauss-Hermite in (-oo,+oo);'
    write ( *, '(a)' ) '  6, 3 point Gauss-Jacobi in [0,1].'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter ICASE, between 1 and 6, or -1 to stop:'
    read ( *, * ) icase
    write ( *, '(a,i8)' ) '  ICASE = ', icase

    if ( icase < 1 .or. 6 < icase ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  ICASE out of bounds.'
      exit
    end if
!
!  Select number of iterative extensions to be performed.
!
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter NSEQ, the number of nested rules to compute:'

    read ( *, * ) nseq
    write ( *, '(a,i8)' ) '  NSEQ = ', nseq

    if ( nseq < 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  NSEQ out of bounds.'
      stop 1
    end if

    if ( icase == 1 ) then
      call test01 ( nseq )
    else if ( icase == 2 ) then
      call test02 ( nseq )
    else if ( icase == 3 ) then
      call test03 ( nseq )
    else if ( icase == 4 ) then
      call test04 ( nseq )
    else if ( icase == 5 ) then
      call test05 ( nseq )
    else if ( icase == 6 ) then
      call test06 ( nseq )
    end if

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS672_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( nseq )

!*****************************************************************************80
!
!! TEST01: extension of 3 point Gauss-Legendre rule.
!
!  Discussion:
!
!    Generate a 3 point Gauss initially from 0 point rule,
!    No pre-assigned nodes.  Symmetry exploited.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 March 2009
!
!  Author:
!
!    Original FORTRAN77 version by Thomas Patterson.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer NSEQ, the number of nested sequences to compute.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: lda = 513
  integer, parameter :: ldb = 2 * lda + 1
  integer, parameter :: ntest = 4

  real ( kind = rk ) err(lda)
  real ( kind = rk ) ext(0:lda)
  real ( kind = rk ) h0
  integer i
  integer ideg
  integer, parameter :: idigit = 8
  integer ierr
  integer iflag
  integer iwork(lda)
  integer k
  integer m
  integer m0
  integer n
  integer, parameter :: nexp = 1021
  integer nodes
  integer nseq
  real ( kind = rk ) pnodes(lda)
  real ( kind = rk ) qi(lda)
  real ( kind = rk ) qr(lda)
  external recura
  logical start
  logical symmet
  real ( kind = rk ) t(0:lda)
  real ( kind = rk ) test(0:ntest)
  real ( kind = rk ) worka(lda,lda)
  real ( kind = rk ) workb(ldb,3)
  real ( kind = rk ) wt(lda)

  n = 0
  m = 3
  m0 = 0
  t(0) = 1.0D+00
  symmet = .true.
  start = .false.

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Extension of a 3 point Gauss-Legendre rule.'
  write ( *, '(a,i8)' ) ' '
  write ( *, '(a,i8)' ) '  N =      ', n
  write ( *, '(a,i8)' ) '  M =      ', m
  write ( *, '(a,i8)' ) '  M0 =     ', m0
  write ( *, '(a,l1)' ) '  SYMMET = ', symmet
  write ( *, '(a,l1)' ) '  START =  ', start

  do i = 1, nseq
!
!  Calculate the extension.
!
    h0 = 2.0D+00
    ideg = n + 2 * m - 1

    call extend ( n, m, m0, t, recura, symmet, start, pnodes, h0, &
      nexp, idigit, wt, nodes, qr, qi, err, ext, iwork, worka, lda, &
      workb, ldb, iflag )
!
!  Tests.
!
    if ( iflag == 0 .or. iflag == 6 ) then
      do k = 0, min ( ntest, ideg / 2 )
        call check ( n, pnodes, wt, k, h0, recura, test(k), ierr )
      end do
    end if
!
!  Print results.
!
    call results ( err, ext, i, ideg, iflag, iwork, lda, m, m0, &
      n, nodes, ntest, pnodes, qi, qr, symmet, t, test, wt )
!
!  On next iteration, add N+1 nodes using the pre-assigned nodes defined by
!  the T polynomial generated in the previous cycle.
!
    m = n + 1
    start = .false.

  end do

  return
end
subroutine test02 ( nseq )

!*****************************************************************************80
!
!! TEST02: extension of a 2 point Lobatto rule.
!
!  Discussion:
!
!    Add one node, pre-assign -1.0 and 1.0.  Symmetry exploited.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 March 2009
!
!  Author:
!
!    Original FORTRAN77 version by Thomas Patterson.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer NSEQ, the number of nested sequences to compute.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: lda = 257
  integer, parameter :: ldb = 2 * lda + 1
  integer, parameter :: ntest = 4

  real ( kind = rk ) err(lda)
  real ( kind = rk ) ext(0:lda)
  real ( kind = rk ) h0
  integer i
  integer ideg
  integer, parameter :: idigit = 8
  integer ierr
  integer iflag
  integer iwork(lda)
  integer k
  integer m
  integer m0
  integer n
  integer, parameter :: nexp = 38
  integer nodes
  integer nseq
  real ( kind = rk ) pnodes(lda)
  real ( kind = rk ) qi(lda)
  real ( kind = rk ) qr(lda)
  external recura
  logical start
  logical symmet
  real ( kind = rk ) t(0:lda)
  real ( kind = rk ) test(0:ntest)
  real ( kind = rk ) worka(lda,lda)
  real ( kind = rk ) workb(ldb,3)
  real ( kind = rk ) wt(lda)

  n = 2
  m = 1
  m0 = 0
  pnodes(1) = 1.0D+00
  pnodes(2) = - 1.0D+00
  symmet = .true.
  start = .true.

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Extension of a Lobatto 2 point rule.'
  write ( *, '(a,i8)' ) ' '
  write ( *, '(a,i8)' ) '  N =      ', n
  write ( *, '(a,i8)' ) '  M =      ', m
  write ( *, '(a,i8)' ) '  M0 =     ', m0
  write ( *, '(a,l1)' ) '  SYMMET = ', symmet
  write ( *, '(a,l1)' ) '  START =  ', start

  do i = 1, nseq

    h0 = 2.0D+00
    ideg = n + 2 * m - 1

    call extend ( n, m, m0, t, recura, symmet, start, pnodes, h0, &
      nexp, idigit, wt, nodes, qr, qi, err, ext, iwork, worka, lda, &
      workb, ldb, iflag )
!
!  Tests.
!
    if ( iflag == 0 .or. iflag == 6 ) then
      do k  =  0, min ( ntest, ideg/2 )
        call check ( n, pnodes, wt, k, h0, recura, test(k), ierr )
      end do
    end if
!
!  Print results.
!
    call results ( err, ext, i, ideg, iflag, iwork, lda, m, m0, &
      n, nodes, ntest, pnodes, qi, qr, symmet, t, test, wt )
!
!  On next iteration, add N-1 nodes using the pre-assigned nodes defined
!  by the T polynomial generated in the previous cycle.
!
    m = n - 1
    start = .false.

  end do

  return
end
subroutine test03 ( nseq )

!*****************************************************************************80
!
!! TEST03: extension of 6-point Radau rule.
!
!  Discussion:
!
!    Add five nodes.  Pre-assign -1.0.  No symmetry.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 March 2009
!
!  Author:
!
!    Original FORTRAN77 version by Thomas Patterson.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer NSEQ, the number of nested sequences to compute.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer lda
  parameter (lda = 257)
  integer ldb
  parameter ( ldb = 2 * lda + 1 )
  integer ntest
  parameter ( ntest = 4 )

  real ( kind = rk ) err(lda)
  real ( kind = rk ) ext(0:lda)
  real ( kind = rk ) h0
  integer i
  integer ideg
  integer, parameter :: idigit = 8
  integer ierr
  integer iflag
  integer iwork(lda)
  integer k
  integer m
  integer m0
  integer n
  integer, parameter :: nexp = 38
  integer nodes
  integer nseq
  real ( kind = rk ) pnodes(lda)
  real ( kind = rk ) qi(lda)
  real ( kind = rk ) qr(lda)
  external recura
  logical start
  logical symmet
  real ( kind = rk ) t(0:lda)
  real ( kind = rk ) test(0:ntest)
  real ( kind = rk ) worka(lda,lda)
  real ( kind = rk ) workb(ldb,3)
  real ( kind = rk ) wt(lda)

  n = 1
  m = 5
  m0 = 0
  pnodes(1) = - 1.0D+00
  symmet = .false.
  start = .true.

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Extension of a 6 point Radau rule.'
  write ( *, '(a,i8)' ) ' '
  write ( *, '(a,i8)' ) '  N =      ', n
  write ( *, '(a,i8)' ) '  M =      ', m
  write ( *, '(a,i8)' ) '  M0 =     ', m0
  write ( *, '(a,l1)' ) '  SYMMET = ', symmet
  write ( *, '(a,l1)' ) '  START =  ', start

  do i = 1, nseq

    h0 = 2.0D+00
    ideg = n + 2 * m - 1

    call extend ( n, m, m0, t, recura, symmet, start, pnodes, h0, &
      nexp, idigit, wt, nodes, qr, qi, err, ext, iwork, worka, lda, &
      workb, ldb, iflag )
!
!  Tests.
!
    if ( iflag == 0 .or. iflag == 6 ) then
      do k = 0, min ( ntest, ideg / 2 )
        call check ( n, pnodes, wt, k, h0, recura, test(k), ierr )
      end do
    end if
!
!  Print results.
!
    call results ( err, ext, i, ideg, iflag, iwork, lda, m, m0, &
      n, nodes, ntest, pnodes, qi, qr, symmet, t, test, wt )
!
!  Add N+1 nodes using the pre-assigned nodes defined by the T polynomial
!  generated in the previous cycle.
!
    m = n + 1
    start = .false.

  end do

  return
end
subroutine test04 ( nseq )

!*****************************************************************************80
!
!! TEST04: extension of 2 point Gauss-Laguerre rule.
!
!  Discussion:
!
!    Generate a 2 point rule initially from a 0 point rule.
!    No pre-assigned nodes are used.
!
!    There is no symmetry.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 March 2009
!
!  Author:
!
!    Original FORTRAN77 version by Thomas Patterson.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer NSEQ, the number of nested sequences to compute.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer lda
  parameter (lda = 257)
  integer ldb
  parameter ( ldb = 2 * lda + 1 )
  integer ntest
  parameter ( ntest = 4 )

  real ( kind = rk ) err(lda)
  real ( kind = rk ) ext(0:lda)
  real ( kind = rk ) h0
  integer i
  integer ideg
  integer, parameter :: idigit = 8
  integer ierr
  integer iflag
  integer iwork(lda)
  integer k
  integer m
  integer m0
  integer n
  integer, parameter :: nexp = 38
  integer nodes
  integer nseq
  real ( kind = rk ) pnodes(lda)
  real ( kind = rk ) qi(lda)
  real ( kind = rk ) qr(lda)
  external recurb
  logical start
  logical symmet
  real ( kind = rk ) t(0:lda)
  real ( kind = rk ) test(0:ntest)
  real ( kind = rk ) worka(lda,lda)
  real ( kind = rk ) workb(ldb,3)
  real ( kind = rk ) wt(lda)

  n = 0
  m = 2
  m0 = 0
  t(0) = 1.0D+00
  symmet = .false.
  start = .false.

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  Extension of 2 point Gauss-Laguerre rule.'
  write ( *, '(a,i8)' ) ' '
  write ( *, '(a,i8)' ) '  N =      ', n
  write ( *, '(a,i8)' ) '  M =      ', m
  write ( *, '(a,i8)' ) '  M0 =     ', m0
  write ( *, '(a,l1)' ) '  SYMMET = ', symmet
  write ( *, '(a,l1)' ) '  START =  ', start

  do i = 1, nseq
!
!  Calculate extension.
!
    h0 = 1.0D+00
    ideg = n + 2 * m - 1

    call extend ( n, m, m0, t, recurb, symmet, start, pnodes, h0, &
      nexp, idigit, wt, nodes, qr, qi, err, ext, iwork, worka, lda, &
      workb, ldb, iflag )
!
!  Tests.
!
    if ( iflag == 0 .or. iflag == 6 ) then
      do k = 0, min ( ntest, ideg / 2 )
        call check ( n, pnodes, wt, k, h0, recurb, test(k), ierr )
      end do
    end if
!
!  Print results.
!
    call results ( err, ext, i, ideg, iflag, iwork, lda, m, m0, &
      n, nodes, ntest, pnodes, qi, qr, symmet, t, test, wt )
!
!  On the next iteration, add N+1 nodes using the pre-assigned nodes defined
!  by the T polynomial generated in the previous cycle.
!
    m = n + 1
    start = .false.

  end do

  return
end
subroutine test05 ( nseq )

!*****************************************************************************80
!
!! TEST05: extension of 3-point Gauss-Hermite rule.
!
!  Discussion:
!
!    Generate a 3-point rule initially from a zero point rule.
!    That is, there are no pre-assigned nodes.
!
!    Symmetry is exploited.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 March 2009
!
!  Author:
!
!    Original FORTRAN77 version by Thomas Patterson.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer NSEQ, the number of nested sequences to compute.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer lda
  parameter (lda = 257)
  integer ldb
  parameter ( ldb = 2 * lda + 1 )
  integer ntest
  parameter ( ntest = 4 )

  real ( kind = rk ) err(lda)
  real ( kind = rk ) ext(0:lda)
  real ( kind = rk ) h0
  integer i
  integer ideg
  integer, parameter :: idigit = 8
  integer ierr
  integer iflag
  integer iwork(lda)
  integer k
  integer m
  integer m0
  integer n
  integer, parameter :: nexp = 38
  integer nodes
  integer nseq
  real ( kind = rk ) pnodes(lda)
  real ( kind = rk ) qi(lda)
  real ( kind = rk ) qr(lda)
  external recurc
  logical start
  logical symmet
  real ( kind = rk ) t(0:lda)
  real ( kind = rk ) test(0:ntest)
  real ( kind = rk ) worka(lda,lda)
  real ( kind = rk ) workb(ldb,3)
  real ( kind = rk ) wt(lda)

  m = 3
  n = 0
  m0 = 0
  t(0) = 1.0D+00
  symmet = .true.
  start = .false.

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  Extension of a 3 point Gauss-Hermite rule.'
  write ( *, '(a,i8)' ) ' '
  write ( *, '(a,i8)' ) '  N =      ', n
  write ( *, '(a,i8)' ) '  M =      ', m
  write ( *, '(a,i8)' ) '  M0 =     ', m0
  write ( *, '(a,l1)' ) '  SYMMET = ', symmet
  write ( *, '(a,l1)' ) '  START =  ', start

  do i = 1, nseq
!
!  Calculate extension.
!  Zero moment integral = sqrt(pi)
!
    h0 = 2.0D+00 * sqrt ( atan ( 1.0D+00 ) )
    ideg = n + 2 * m - 1

    call extend ( n, m, m0, t, recurc, symmet, start, pnodes, h0, &
      nexp, idigit, wt, nodes, qr, qi, err, ext, iwork, worka, lda, &
      workb, ldb, iflag )
!
!  Tests.
!
    if ( iflag == 0 .or. iflag == 6 ) then
      do k = 0, min ( ntest, ideg / 2 )
        call check ( n, pnodes, wt, k, h0, recurc, test(k), ierr )
      end do
    end if
!
!  Print results.
!
    call results ( err, ext, i, ideg, iflag, iwork, lda, m, m0, &
      n, nodes, ntest, pnodes, qi, qr, symmet, t, test, wt )
!
!  On the next iteration, add N+1 nodes using the pre-assigned nodes defined
!  by the T polynomial enerated in the previous cycle.
!
    m = n + 1
    start = .false.

  end do

  return
end
subroutine test06 ( nseq )

!*****************************************************************************80
!
!! TEST06: extension of 3-point Gauss-Jacobi rule for weight sqrt(x) in [0,1].
!
!  Discussion:
!
!    Generate 3-point rule initially from zero point rule.
!
!    No pre-assigned nodes, no symmetry.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 March 2009
!
!  Author:
!
!    Original FORTRAN77 version by Thomas Patterson.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer NSEQ, the number of nested sequences to compute.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer lda
  parameter (lda = 257)
  integer ldb
  parameter ( ldb = 2 * lda + 1 )
  integer ntest
  parameter ( ntest = 4 )

  real ( kind = rk ) err(lda)
  real ( kind = rk ) ext(0:lda)
  real ( kind = rk ) h0
  integer i
  integer ideg
  integer, parameter :: idigit = 8
  integer ierr
  integer iflag
  integer iwork(lda)
  integer k
  integer m
  integer m0
  integer n
  integer, parameter :: nexp = 38
  integer nodes
  integer nseq
  real ( kind = rk ) pnodes(lda)
  real ( kind = rk ) qi(lda)
  real ( kind = rk ) qr(lda)
  external recurd
  logical start
  logical symmet
  real ( kind = rk ) t(0:lda)
  real ( kind = rk ) test(0:ntest)
  real ( kind = rk ) worka(lda,lda)
  real ( kind = rk ) workb(ldb,3)
  real ( kind = rk ) wt(lda)

  m = 3
  n = 0
  m0 = 0
  t(0) = 1.0D+00
  symmet = .false.
  start = .false.

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  Extension of a 3 point Gauss-Jacobi rule.'
  write ( *, '(a,i8)' ) ' '
  write ( *, '(a,i8)' ) '  N =      ', n
  write ( *, '(a,i8)' ) '  M =      ', m
  write ( *, '(a,i8)' ) '  M0 =     ', m0
  write ( *, '(a,l1)' ) '  SYMMET = ', symmet
  write ( *, '(a,l1)' ) '  START =  ', start

  do i = 1, nseq
!
!  Calculate extension.
!
    h0 = 2.0D+00 / 3.0D+00
    ideg = n + 2 * m - 1

    call extend ( n, m, m0, t, recurd, symmet, start, pnodes, h0, &
      nexp, idigit, wt, nodes, qr, qi, err, ext, iwork, worka, lda, &
      workb, ldb, iflag )
!
!  Tests.
!
    if ( iflag == 0 .or. iflag == 6 ) then
      do k = 0, min ( ntest, ideg / 2 )
       call check ( n, pnodes, wt, k, h0, recurd, test(k), ierr )
      end do
    end if
!
!  Print results.
!
    call results ( err, ext, i, ideg, iflag, iwork, lda, m, m0, &
      n, nodes, ntest, pnodes, qi, qr, symmet, t, test, wt )
!
!  On the next iteration, add N+1 nodes using the pre-assigned nodes defined
!  by the T polynomial generated in the previous cycle.
!
    m = n + 1
    start = .false.

  end do

  return
end
subroutine recura ( k, c, d, e )

!*****************************************************************************80
!
!! RECURA is the recurrence used for Gauss, Lobatto and Radau.
!
!  Discussion:
!
!    This is an example of a user supplied subroutine to define the
!    orthogonal polynomials.
!
!    RECURA ( K, C, D, E ) gives the coefficients C, D and E such that:
!
!      P(K+1,X) = ( C * X + D ) * P(K,X) + E * P(K-1,X)
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 March 2009
!
!  Author:
!
!    Original FORTRAN77 version by Thomas Patterson.
!    FORTRAN90 version by John Burkardt.
!
!  Input:
!
!    integer K, the index.
!
!  Output:
!
!    real ( kind = rk ) C, D, E, the recurrence relation parameters.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) c
  real ( kind = rk ) d
  real ( kind = rk ) e
  real ( kind = rk ) f
  integer k

  f  =  real ( k + 1, kind = rk )
  c  =  real ( 2 * k + 1, kind = rk ) / f
  d  =  0.0D+00
  e  =  - real ( k, kind = rk ) / f

  return
end
subroutine recurb ( k, c, d, e )

!*****************************************************************************80
!
!! RECURB is the recurrence used for Laguerre rules.
!
!  Discussion:
!
!    This is an example of a user supplied subroutine to define the
!    orthogonal polynomials.
!
!    RECURB ( K, C, D, E ) gives the coefficients C, D and E such that:
!
!      P(K+1,X) = ( C * X + D ) * P(K,X) + E * P(K-1,X)
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 March 2009
!
!  Author:
!
!    Original FORTRAN77 version by Thomas Patterson.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer K, the index.
!
!    Output, real ( kind = rk ) C, D, E, the recurrence relation parameters.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) c
  real ( kind = rk ) d
  real ( kind = rk ) e
  real ( kind = rk ) f
  integer k

  f = real ( k + 1, kind = rk )
  c = - 1.0D+00 / f
  d = real ( 2 * k + 1, kind = rk ) / f
  e = - real ( k, kind = rk ) / f

  return
end
subroutine recurc ( k, c, d, e )

!*****************************************************************************80
!
!! RECURC is the recurrence used for Hermite.
!
!  Discussion:
!
!    This is an example of a user supplied subroutine to define the
!    orthogonal polynomials.
!
!    RECURC ( K, C, D, E ) gives the coefficients C, D and E such that:
!
!      P(K+1,X) = ( C * X + D ) * P(K,X) + E * P(K-1,X)
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 March 2009
!
!  Author:
!
!    Original FORTRAN77 version by Thomas Patterson.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer K, the index.
!
!    Output, real ( kind = rk ) C, D, E, the recurrence relation parameters.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) c
  real ( kind = rk ) d
  real ( kind = rk ) e
  integer k

  c = 1.0D+00
  d = 0.0D+00
  e = - real ( k, kind = rk ) / 2.0D+00

  return
end
subroutine recurd ( k, c, d, e )

!*****************************************************************************80
!
!! RECURD is the recurrence used for Jacobi.
!
!  Discussion:
!
!    Jacobi polynomials in [0,1].
!
!    The weight function is (1-x)^(p-q) * x^(q-1).
!
!    This case for weight sqrt(x), p = 3/2 and q = p.
!
!    This is an example of a user supplied subroutine to define the
!    orthogonal polynomials.
!
!    RECURD ( K, C, D, E ) gives the coefficients C, D and E such that:
!
!      P(K+1,X) = ( C * X + D ) * P(K,X) + E * P(K-1,X)
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 March 2009
!
!  Author:
!
!    Original FORTRAN77 version by Thomas Patterson.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer K, the index.
!
!    Output, real ( kind = rk ) C, D, E, the recurrence relation parameters.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a1
  real ( kind = rk ) a2
  real ( kind = rk ) a3
  real ( kind = rk ) a4
  real ( kind = rk ) b
  real ( kind = rk ) bp1
  real ( kind = rk ) c
  real ( kind = rk ) d
  real ( kind = rk ) e
  real ( kind = rk ) f2k
  real ( kind = rk ) fk
  integer k
  real ( kind = rk ) p
  real ( kind = rk ) q
  real ( kind = rk ) x3

  p = 1.5D+00
  q = p
  fk = real ( k, kind = rk )
  f2k = 2.0D+00 * fk
  b = f2k + p
  bp1 = b + 1.0D+00
  x3 = ( b - 2.0D+00 ) * ( b - 1.0D+00 ) * b
  a1 = ( b - 1.0D+00 ) * bp1 * x3
  a2 = - x3 * ( f2k * ( fk + p ) + q * ( p - 1.0D+00 ) )
  a3 = x3 * bp1 * ( b - 1.0D+00 )
  a4 = fk * ( fk + q - 1.0D+00 ) * ( fk + p - 1.0D+00 ) * ( fk + p - q ) * bp1
  c = a3 / a1
  d = a2 / a1
  e = - a4 / a1

  return
end
subroutine results ( err, ext, i, ideg, iflag, iwork, lda, m, m0, &
  n, nodes, ntest, pnodes, qi, qr, symmet, t, test, wt )

!*****************************************************************************80
!
!! RESULTS prints the results.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 March 2009
!
!  Author:
!
!    Original FORTRAN77 version by Thomas Patterson.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = rk ) ERR(M), a measure of the relative error in the
!    nodes.  This may be inspected if the convergence error flag has been raised
!    (IFLAG = 3) to decide if the nodes in question are acceptable.  (ERR(*)
!    actually gives the mean last correction to the quadratic factor in the
!    generalized Bairstow root finder (see BAIR).
!
!    Input, real ( kind = rk ) EXT(0:M), the coefficients of the polynomial whose
!    roots are the  extended nodes (QRNODES(*),QINODES(*)) and expressed as:
!      EXT = SUM (I = 0 to M) EXT(I)*P(I,X).
!
!    Input, integer I, the current stage of the set of nested
!    quadrature rule.
!
!    Input, integer IDEG, the expected polynomial accuracy of 
!    the quadrature rule.
!
!    Input, integer IFLAG, error flag.
!    * 0, No error detected;
!    * 1, The linear system of equations defining the polynomial
!      whose roots are the extended nodes became singular or
!      very  ill-conditioned.   (FATAL).
!    * 2, The linear system of equations used to generate the
!      polynomial T when START is TRUE became singular
!      or very ill-conditioned. (FATAL).
!    * 3, Poor convergence has been detected in the calculation
!      of the roots of EXT corresponding to the new
!      nodes or all nodes have not been found (M not equal
!      to NODES). See also ERR(*).
!    * 4, Possible imaginary nodes detected.
!    * 5, Value of N and M incompatible for SYMMET = TRUE.
!      Both cannot be odd. (FATAL)
!    * 6, Test of new quadrature rule has failed.
!
!    Input, integer IWORK(max(M,N)), convergence flags.  Elements 
!    1 to NODES give information on the convergence of the roots of the 
!    polynomial EXT corresponding to each extended node.
!    * 0, Convergence of I th root satisfactory;
!    * 1, Convergence of I th root unsatisfactory.
!
!    Input, integer LDA, the leading dimension of WORKA in the
!    calling program.
!
!    Input, integer M, the number of nodes to be optimally added.
!
!    Input, integer M0, the lower limit to the expansion of T.
!
!    Input, integer N, the number of nodes in the current 
!    quadrature rule.
!
!    Input, integer NODES, the number of extended nodes found.
!    Normally this equals M, but NODES will be less than M in cases where the
!    computation was terminated because an error condition was encountered.
!
!    Input, integer NTEST, the number of tests made of the 
!    quadrature rule.
!
!    Input, real ( kind = rk ) PNODES(N), the nodes of the extended quadrature 
!    rule.
!
!    Input, real ( kind = rk ) QINODE(M), the imaginary parts of the extended
!    nodes.
!
!    Input, real ( kind = rk ) QRNODE(M), the real parts of the extended nodes.
!
!    Input, logical SYMMET,
!    * FALSE, if no advantage is to be taken of symmetry, if any,
!      about x = 0 in the interval of integration and the
!      orthogonality  weight function. Note that if symmetry in
!      fact does exist setting this parameter to zero will still
!      produce correct results - only efficiency is effected.
!    * TRUE, if the extended rule computations should
!      exploit symmetry about x = 0 in the interval of
!      integration and the orthogonality  weight function.
!      This reduces the size of the system of linear equations
!      determining EXT by a factor of about 2 (see WORKA). If
!      symmetry does not in fact exist erroneous results will be
!      produced.
!
!    Input, real ( kind = rk ) T(0:M); the coefficients TI of the orthogonal
!    expansion whose roots are the nodes of the extended quadrature rule
!    (that is, the pre-assigned nodes plus the extended nodes) and expressed as:
!      SUM (I = M to N+M) (TI/HI)*P(I,X)
!    T(I-M) holds the value of TI.
!
!    Input, real ( kind = rk ) TEST(NTEST), the results of tests on the
!    quadrature rule.
!
!    Input, real ( kind = rk ) WT(N), the quadrature weights for
!    the extended rule associated with the nodes in PNODES.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer lda
  integer ntest

  real ( kind = rk ) err(lda)
  real ( kind = rk ) ext(0:lda)
  integer i
  integer ideg
  integer iflag
  integer iwork(lda)
  integer j
  integer k
  integer m
  integer m0
  integer n
  integer nodes
  integer num
  real ( kind = rk ) pnodes(lda)
  real ( kind = rk ) qi(lda)
  real ( kind = rk ) qr(lda)
  logical symmet
  real ( kind = rk ) t(0:lda)
  real ( kind = rk ) test(0:ntest)
  real ( kind = rk ) wt(lda)
!
!  Display results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Iteration', i
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Coefficients of expansion whose roots are the new nodes:'
  write ( *, '(a)' ) ' '
  do j = 0, m
    write ( *, '(d25.16,a,i3,a)' )  ext(j), ' * P(', j, ',X)'
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  New nodes'
  write ( *, '(a)' )  '                     real             Imaginary' // &
    ' Flag       Err'
  write ( *, '(a)' ) ' '
  do k = 1, nodes
    write ( *, '(2d25.16,i5,d10.1)') qr(k), qi(k), iwork(k), err(k)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  New full extended expansion'
  write ( *, '(a)' ) ' '
  do j = m0, n
    write ( *, '(d25.16,a,i3,a)')  t(j-m0), ' * P(',j, ',X)/HI'
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a,i2,a,i3,a,i1,a,i3)' ) &
    '  Complete extended rule: STEP = ', i, &
    '  POINTS = ', n, &
    '  IFLAG = ', iflag, &
    '  Nodes added = ', nodes

  if ( iflag /= 0 .and. iflag /= 6 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  IFLAG = ', iflag
    write ( *, '(a)' ) '    Computation terminated prematurely.'
    return
  end if
!
!  Print rule (positive nodes only if symmetry present).
!
  if ( symmet ) then
    num = ( n + 1 ) / 2
  else
    num = n
  end if
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  No.                     Node                   Weight'
  write ( *, '(a)' ) ' '
  do j = 1, num
    write ( *, '(i5,2d25.16)' ) j, pnodes(j), wt(j)
  end do
!
!  Display test results.
!
  write ( *, '(a)' ) ' '
  do k = 0, min ( ntest, ideg / 2 )
    write ( *, '(a,i2,a,d25.16)' ) '  Test(', k, ') = ', test(k)
  end do

  if ( iflag == 6 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  IFLAG = 6,'
    write ( *, '(a)' ) '    The rule test is unsatisfactory.'
  end if

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! timestamp prints the current YMDHMS date as a time stamp.
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

