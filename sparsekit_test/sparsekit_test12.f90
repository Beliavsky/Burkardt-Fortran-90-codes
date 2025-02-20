program main

!*****************************************************************************80
!
!! MAIN is a test program for exponential propagator using Arnoldi approach.
!
!  Discussion:
!
!    This main program is a very simple test using diagonal matrices
!    (Krylov subspace methods are blind to the structure of the matrix
!    except for symmetry). This provides a good way of testing the
!    accuracy of the method as well as the error estimates.
!
!  Modified:
!
!    02 July 2005
!
  implicit none

  integer, parameter :: nmax = 150
  integer, parameter :: ih0 = 60
  integer, parameter :: ndmx = 20

  integer, parameter :: nzmax = 7 * nmax

  real ( kind = 8 ) a(nzmax)
  real ( kind = 8 ) :: a0 = 0.0
  real ( kind = 8 ) :: b0 = 1.0
  real ( kind = 8 ) ddot
  real ( kind = 8 ) eps
  real ( kind = 8 ) h
  integer ioff(10)
  integer j
  integer k
  integer m
  integer n
  integer ndiag
  real ( kind = 8 ) t
  real ( kind = 8 ) tn
  real ( kind = 8 ) u(ih0*nmax)
  real ( kind = 8 ) w(nmax)
  real ( kind = 8 ) w1(nmax)
  real ( kind = 8 ) x(nmax)
  real ( kind = 8 ) y(nmax)

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSEKIT_PRB12'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test EXPPRO, which computes the matrix exponential.'
!
!  Set dimension of matrix.
!
  n = 100
!
!  Define matrix.
!  A is a single diagonal matrix (ndiag = 1 and ioff(1) = 0 )
!
  ndiag = 1
  ioff(1) = 0
!
!  entries in the diagonal are uniformly distributed.
!
  h = 1.0 / real ( n + 1, kind = 8 )
  do j = 1, n
    a(j) = real ( j + 1, kind = 8 ) / real ( n + 1, kind = 8 )
  end do
!
!  Set a value for TN
!
  tn = 2.0

  eps = 0.0001

  m = 5
  write ( *, '(a,i6)' ) '  Dimension of Krylov subspace M = ', m
!
!  Define initial conditions: chosen so that solution = (1,1,1,1..1)^T
!
  do j = 1, n
    w(j) = exp ( a(j) * tn )
    w1(j) = w(j)
  end do

  call expprod ( n, m, eps, tn, u, w, x, y, a, ioff, ndiag )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  First 10 components of final answer:'
  WRITE(*,*)(w(k),k=1,10)
  write ( *, '(a)' ) ' '

  do k = 1, n
    w1(k) = exp ( -a(k) * tn ) * w1(k)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  First 10 components of exact solution '
  WRITE(*,*)(w1(k),k=1,10)
!
!  Compute actual 2-norm of error.
!
  t = 0.0
  do k = 1, n
    t = t + ( w1(k) - w(k) )**2
  end do
  t = sqrt ( t / ddot ( n, w, 1, w, 1 ) )

  write ( *, '(a)' ) ' '
  write ( *, * ) '  RMS error (approx-exact)=', t

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSEKIT_PRB12'
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
