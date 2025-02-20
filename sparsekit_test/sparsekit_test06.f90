program main

!*****************************************************************************80
!
!! SPARSEKIT_PRB06 demonstrates how to read a Harwell-Boeing sparse matrix file.
!
  implicit none

  integer, parameter :: nmax = 500
  integer, parameter :: nzmax = 7000

  real ( kind = 8 ) a(nzmax)
  real ( kind = 8 ) a1(nzmax)
  character ( len = 80 ) filnam
  character ( len = 2 ) guesol
  integer ia(nmax+1)
  integer ia1(nmax+1)
  integer ierr
  integer iin
  integer ios
  integer iout
  integer ja(nzmax)
  integer ja1(nzmax)
  integer job
  character ( len = 8 ) key
  integer ncol
  integer nnz
  integer nrhs
  integer nrow
  real ( kind = 8 ) rhs(1)
  character ( len = 72 ) title
  character ( len = 3 ) type
  logical valued

  iout = 6
  job = 2
  nrhs = 0

  call timestamp ( )

  write ( *, * ) ' '
  write ( *, * ) 'SPARSEKIT_PRB06'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, * ) '  This program demonstrates the use of the SPARSKIT'
  write ( *, * ) '  routines READMT and DINFO1 to read and report on'
  write ( *, * ) '  a sparse matrix stored in a file in the format'
  write ( *, * ) '  used by the Harwell-Boeing Sparse Matrix Collection'
  write ( *, * ) '  or "HBSMC".'
  write ( *, * ) ' '

  filnam = 'saylor_hb.txt'
  iin = 20

  open ( unit = iin, file = filnam, status = 'old', iostat = ios )

  if ( ios < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Error!'
    write ( *, '(a)' ) '  Unable to open file.'
    stop
  end if

  call readmt ( nmax, nzmax, job, iin, a, ja, ia, rhs, nrhs, &
    guesol, nrow, ncol, nnz, title, key, type, ierr )

  close ( unit = iin )
!
!  If not readable, return.
!
  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Error!'
    write ( *, '(a)' ) '  Unable to read matrix.'
    write ( *, '(a,i6)' ) '  READMT returned IERR = ', ierr
    stop
  end if

  valued = ( 2 <= job )

  call dinfo1 ( ncol, iout, a, ja, ia, valued, title, key, type, a1, ja1, ia1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSEKIT_PRB06'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
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
