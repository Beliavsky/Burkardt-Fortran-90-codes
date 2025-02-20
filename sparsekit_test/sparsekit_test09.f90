program main

!*****************************************************************************80
!
!! MAIN generates 5 point and 7-point matrices in Harwell-Boeing format. 
!
!  Creates a file containing a harwell-boeing matrix.
!
!  nz = 1 will create a 2-D problem
!
  implicit none

  integer, parameter :: nxmax = 50
  integer, parameter :: nmx = nxmax * nxmax

  real ( kind = 8 ) a(7*nmx)
  character ( len = 2 ) guesol
  integer ia(nmx)
  integer iau(nmx)
  integer ifmt
  integer iout
  integer ja(7*nmx)
  integer job
  character ( len = 8 ) key
  character ( len = 50 ) matfil
  integer n
  integer nx
  integer ny
  integer nz
  real ( kind = 8 ) rhs(1)
  real ( kind = 8 ) stencil(7)
  character ( len = 72 ) title
  character ( len = 3 ) type

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSEKIT_PRB09:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This program demonstrates the use of GEN57PT'
  write ( *, '(a)' ) '  to generate a sparse matrix derived from a 5 or'
  write ( *, '(a)' ) '  7 point stencil on an NX by NY by NZ grid.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix is then stored in Harwell-Boeing format'
  write ( *, '(a)' ) '  in a file, using routine PRTMT.'
!
!  Set data defining the matrix.
!
  nx = 10
  ny = 10
  nz = 1
!
!  Call GEN57PT to generate the matrix.
!
  call gen57pt(nx,ny,nz,a,ja,ia,iau,stencil)
!
!  Set parameters required for the Harwell-Boeing format.
!
  n = nx * ny * nz
  guesol = 'NN'
  title = ' 5-POINT TEST MATRIX FROM SPARSKIT                    '
  type = 'RUA'
  key = 'SC5POINT'
  ifmt = 104
  job = 2
!
!  Write matrix to file.
!
  iout = 7
  matfil = 'test.mat'

  open ( unit = IOUT, file = MATFIL, STATUS = 'replace' )

  call prtmt ( n, n, a, ja, ia, rhs, guesol,title,key,type,ifmt,job,iout)

  close ( unit = iout )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSEKIT_PRB09'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
function afun (x,y,z)

!*****************************************************************************80
!
  implicit none

  REAL afun, x,y, z

  afun = -1.0

  return
end
function bfun (x,y,z)

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) bfun, x,y, z

  bfun = -1.0

  return
end
function cfun (x,y,z)

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) cfun, x,y, z

  cfun = -1.0

  return
end
function dfun (x,y,z)

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) dfun, gamma, x,y, z
  data gamma /100.0/

  dfun = 10.0

  return
end
function efun (x,y,z)

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) efun, gamma, x,y, z

  data gamma /100.0/

  efun = 0.0

  return
end
function ffun (x,y,z)

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) ffun, x,y, z

  ffun = 0.0

  return
end
function gfun (x,y,z)

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) gfun, x,y, z

  gfun = 0.0

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
