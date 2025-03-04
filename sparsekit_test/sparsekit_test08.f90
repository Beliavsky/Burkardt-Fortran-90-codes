program main

!*****************************************************************************80
!
!! SPARSEKIT_PRB08 runs the Zlatev test suite.
!
!  Discussion:
!
!    Three matrices are generated, written to separate files
!    in the Harwell-Boeing format.
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSEKIT_PRB08'
  write ( *, '(a)' ) '  FORTRAN90 version'

  call test01
  call test02
  call test03

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSEKIT_PRB08'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01

!*****************************************************************************80
!
!! TEST01 runs the first test.
!
  implicit none

  integer, parameter :: nmax = 1000

  integer, parameter :: nzmax = 20 * nmax

  real ( kind = 8 ) a(nzmax)
  real ( kind = 8 ) alpha
  character ( len = 2 ) guesol
  integer ia(nzmax)
  integer ic
  integer ierr
  integer ifmt
  integer indx
  integer ios
  integer iout
  integer iwk(nmax)
  integer ja(nzmax)
  integer job
  character ( len = 3 ) key
  integer m
  integer n
  integer nn
  integer nz
  real ( kind = 8 ) rhs(1)
  character ( len = 72 ) title
  character ( len = 8 ) type

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Write Zlatev example to file.'
!
!  Call MATRF2 to set up the matrix in COO (coordinate) format.
!
  m = 100
  n = m
  ic = n / 2
  indx = 10
  alpha = 5.0D+00
  nn = nzmax

  call matrf2 ( m, n, ic, indx, alpha, nn, nz, a, ia, ja, ierr )
!
!  Convert the matrix from COO to CSR format.
!
  job = 1
  call coocsr_inplace ( n, nz, job, a, ia, ja, iwk )
!
!  Write the matrix to a file using Harwell-Boeing format.
!
  title = 'First matrix from Zlatev examples               '
  type = 'RUA'
  key = ' ZLATEV1'
  guesol = 'NN'
  ifmt = 3
  job = 2
  iout = 7

  open ( unit = iout, file = 'zlatev1_hb.txt', status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST01 - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    return
  end if

  call prtmt ( n, n, a, ia, ja, rhs, guesol, title, type, key, &
    ifmt, job, iout )

  close ( unit = iout )

  return
end
subroutine test02

!*****************************************************************************80
!
!! TEST02 runs the second test, with DCN matrices.
!
  implicit none

  integer, parameter :: nmax = 1000
  integer, parameter :: nzmax = 20 * nmax

  real ( kind = 8 ) a(nzmax)
  character ( len = 2 ) guesol
  integer ia(nzmax)
  integer ic
  integer ierr
  integer ifmt
  integer ios
  integer iout
  integer iwk(nmax)
  integer ja(nzmax)
  integer job
  character ( len = 3 ) key
  integer n
  integer ne
  integer nn
  real ( kind = 8 ) rhs(1)
  character ( len = 72 ) title
  character ( len = 8 ) type 

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Write DCN example to file.'
!
!  Call DCN to set up the matrix in COO format.
!
  n = 200
  nn = nzmax
  ic = 20
  call dcn ( a, ia, ja, n, ne, ic, nn, ierr )
!
!  Convert from COO to CSR format.
!
  job = 1

  call coocsr_inplace ( n, ne, job, a, ia, ja, iwk )
!
!  Write the matrix to file using Harwell Boeing format.
!
  title = 'second matrix from Zlatev examples               '
  guesol = 'NN'
  type = 'RUA'
  key = ' ZLATEV2'
  ifmt = 3
  job = 2
  iout = 7

  open ( unit = iout, file = 'zlatev2_hb.txt', status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST02 - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    return
  end if

  call prtmt ( n, n, a, ia, ja, rhs, guesol, title, type, key, &
    ifmt, job, iout )

  close ( unit = iout )

  return
end
subroutine test03

!*****************************************************************************80
!
!! TEST03 runs the test with ECN matrices.
!
  implicit none

  integer, parameter :: nmax = 1000
  integer, parameter :: nzmax = 20 * nmax

  real ( kind = 8 ) a(nzmax)
  character ( len = 2 ) guesol
  integer ia(nzmax)
  integer ic
  integer ierr
  integer ifmt
  integer ios
  integer iout
  integer iwk(nmax)
  integer ja(nzmax)
  integer job
  character ( len = 3 ) key
  integer n
  integer ne
  integer nn
  real ( kind = 8 ) rhs(1)
  character ( len = 72 ) title
  character ( len = 8 ) type 

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Write ECN example to file.'
!
!  Call ECN to set up the matrix in COO format.
!
  n = 200
  ic = 20
  nn = nzmax
  call ecn ( n, ic, ne, ia, ja, a, nn, ierr )
!
!  Convert the matrix from COO to CSR format.
!
  job = 1
  call coocsr_inplace ( n, ne, job, a, ia, ja, iwk )
!
!  Store matrix to a file using Harwell Boeing format.
!
  title = 'Third matrix from Zlatev examples               '
  guesol = 'NN'
  type = 'RUA'
  key = ' ZLATEV3'
  ifmt = 3
  job = 2
  iout = 7

  open ( unit = iout, file = 'zlatev3_hb.txt', status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST03 - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    return
  end if

  call prtmt ( n, n, a, ia, ja, rhs, guesol, title, type, key, &
    ifmt, job, iout )

  close ( unit = iout )

  return
end
function afun ( x, y, z )

!*****************************************************************************80
!
!! AFUN
!
  implicit none

  real ( kind = 8 ) afun
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  afun = -1.0D+00

  return
end
function bfun ( x, y, z )

!*****************************************************************************80
!
!! BFUN
!
  implicit none

  real ( kind = 8 ) bfun
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  bfun = -1.0D+00

  return
end
function cfun ( x, y, z )

!*****************************************************************************80
!
!! CFUN
!
  implicit none

  real ( kind = 8 ) cfun
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  cfun = -1.0D+00

  return
end
function dfun ( x, y, z )

!*****************************************************************************80
!
!! DFUN
!
  implicit none

  real ( kind = 8 ) dfun
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  dfun = 0.0D+00

  return
end
function efun ( x, y, z )

!*****************************************************************************80
!
!! EFUN
!
  implicit none

  real ( kind = 8 ) efun
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  efun = 0.0D+00

  return
end
function ffun ( x, y, z )

!*****************************************************************************80
!
!! FFUN
!
  implicit none

  real ( kind = 8 ) ffun
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  ffun = 0.0D+00

  return
end
function gfun ( x, y, z )

!*****************************************************************************80
!
!! GFUN
!
  implicit none

  real ( kind = 8 ) gfun
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  gfun = 0.0D+00

  return
end
subroutine afunbl ( nfree, x, y, z, coeff )

!*****************************************************************************80
!
!! AFUNBL
!
  implicit none

  real ( kind = 8 ) coeff(100)
  integer i
  integer j
  integer nfree
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  do j=1, nfree
    do i=1, nfree
      coeff((j-1)*nfree+i) = 0.0D+00
    end do
    coeff((j-1)*nfree+j) = -1.0D+00
  end do

  return
end
subroutine bfunbl ( nfree, x, y, z, coeff )

!*****************************************************************************80
!
!! BFUNBL
!
  implicit none

  real ( kind = 8 ) coeff(100)
  integer i
  integer j
  integer nfree
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  do j=1, nfree
    do i=1, nfree
      coeff((j-1)*nfree+i) = 0.0D+00
    end do
    coeff((j-1)*nfree+j) = -1.0D+00
  end do

  return
end
subroutine cfunbl ( nfree, x, y, z, coeff )

!*****************************************************************************80
!
!! CFUNBL
!
  implicit none

  real ( kind = 8 ) coeff(100)
  integer i
  integer j
  integer nfree
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  do j=1, nfree
    do i=1, nfree
      coeff((j-1)*nfree+i) = 0.0D+00
    end do
    coeff((j-1)*nfree+j) = -1.0D+00
  end do

  return
end
subroutine dfunbl ( nfree, x, y, z, coeff )

!*****************************************************************************80
!
!! DFUNBL 
!
  implicit none

  real ( kind = 8 ) coeff(100)
  integer i
  integer j
  integer nfree
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  do j=1, nfree
    do i=1, nfree
      coeff((j-1)*nfree+i) = 0.0D+00
    end do
  end do

  return
end
subroutine efunbl ( nfree, x, y, z, coeff )

!*****************************************************************************80
!
!! EFUNBL
!
  implicit none

  real ( kind = 8 ) coeff(100)
  integer i
  integer j
  integer nfree
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  do j=1, nfree
    do i=1, nfree
      coeff((j-1)*nfree+i) = 0.0D+00
    end do
  end do

  return
end
subroutine ffunbl ( nfree, x, y, z, coeff )

!*****************************************************************************80
!
!! FFUNBL
!
  implicit none

  real ( kind = 8 ) coeff(100)
  integer i
  integer j
  integer nfree
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  do j=1, nfree
    do i=1, nfree
      coeff((j-1)*nfree+i) = 0.0D+00
    end do
  end do

  return
end
subroutine gfunbl ( nfree, x, y, z, coeff )

!*****************************************************************************80
!
!! GFUNBL
!
  implicit none

  real ( kind = 8 ) coeff(100)
  integer i
  integer j
  integer nfree
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  do j=1, nfree
    do i=1, nfree
      coeff((j-1)*nfree+i) = 0.0D+00
    end do
  end do

  return
end
