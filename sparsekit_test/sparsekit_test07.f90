program main

!*****************************************************************************80
!
!! MAIN generates a Markov chain matrix to test eigenvalue routines.
!
!  The matrix A produced by this routine can also be used to generate
!  a related singular matrix, I-A.
!  The matrix models  simple random walk on a triangular grid.
!  see additional comments in subroutine.
!
!   will create a matrix in the HARWELL/BOEING format and put it in
!   the file markov.mat
!
  implicit none

  integer, parameter :: nmax = 5000
  integer, parameter :: nzmax= 4 * nmax

  real ( kind = 8 ) a(nzmax)
  integer ia(nmax+1)
  integer ifmt
  integer ios
  integer iout
  integer ja(nzmax)
  integer job
  character ( len = 8 ) key
  integer m
  integer n
  character ( len = 72 ) title
  character ( len = 3 ) type
  real ( kind = 8 ) rhs(1)

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSEKIT_PRB07'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Generate a Markov chain matrix to test the'
  write ( *, '(a)' ) '  eigenvalue routine.'

  open ( unit = 11, file = 'markov.mat', status = 'replace', iostat = ios )

  if ( ios < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPARSEKIT_PRB07 - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file.'
    stop
  end if
!
!  Set grid size - will not accept too large grids.
!
  m = 5

  if ( 2 * nmax < m * ( m + 1 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPARSEKIT_PRB07 - Fatal error!'
    write ( *, '(a)' ) '  M is too large - unable to produce matrix.'
    stop
  end if
!
!  Call the matrix generator.
!
  call markgen ( m, n, a, ja, ia )
!
!  Store result in file.
!
  title = ' Test matrix from SPARSKIT - Markov chain model           '
  key = 'randwk01'
  type = 'rua'
  iout = 11
  job = 2
  ifmt = 10

  call prtmt ( n, n, a, ja, ia, rhs, 'NN', title, key, type, ifmt, job, iout )

  close ( unit = iout )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix has been stored in Harwell/Boeing format'
  write ( *, '(a)' ) '  in the file "markov.mat".'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSEKIT_PRB07'
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
