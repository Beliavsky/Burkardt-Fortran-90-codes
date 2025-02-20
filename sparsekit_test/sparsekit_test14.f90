program main

!*****************************************************************************80
!
!! SPARSEKIT_test14 demonstrates CSR to NCF conversion.
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSEKIT_test14'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This test demonstrates the conversion of a'
  write ( *, '(a)' ) '  sparse matrix from CSR to NCF formats.'

  call test01

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSEKIT_test14'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01

!*****************************************************************************80
!
!! TEST01 sets up a CSR matrix and converts it to NCF form.
!
!  Discussion:
!
!    The matrix is
!
!    11  12   0  14   0
!    21  22  23   0  25
!     0  32  33  34   0
!    41   0  43  44  45
!     0  52   0  54  55
!
!    The CSR format is
!
!    A = ( 11, 12, 14, 21, 22, 23, 25, 32, 33, 34, 41, 43, 44, 45, 52, 54, 55 )
!    
!    IA = ( 1,          4,              8,         11,             15,
!          18 )
!
!    JA = ( 1,  2,  4,  1,  2,  3,  5,  2,  3,  4,  1,  3,  4,  5,  2, 4,  5 )
!
!    The NCF format is
!
!    COEF =   ( 11, 22, 33, 44, 55, 12, 14, 21, 23, 25, 32, 34, 41, 43, 45, 52, 54 )
!    JCOEF' = ( 1,  2,  3,  4,  5,  1,  1,  2,  2,  2,  3,  3,  4,  4,  4,  5,  5  )
!             (  1,  2,  3,  4,  5,  2,  4,  1,  3,  5,  2,  4,  1,  3,  5,  2,  4 )
!
  implicit none

  integer, parameter :: n = 5
  integer, parameter :: nzmax = 17

  real ( kind = 8 ) a(nzmax)
  real ( kind = 8 ) coef(nzmax)
  integer i
  integer ia(n+1)
  integer ierr
  integer ja(nzmax)
  integer jcoef(nzmax,2)
  integer nz

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Set up a CSR matrix, and call'
  write ( *, '(a)' ) '  CSRNCF to convert it to NCF format.'

  a(1:nzmax) = (/ &
    11.0D+00, 12.0D+00,           14.0D+00,           &
    21.0D+00, 22.0D+00, 23.0D+00,           25.0D+00, &
              32.0D+00, 33.0D+00, 34.0D+00,           &
    41.0D+00,           43.0D+00, 44.0D+00, 45.0D+00, &
              52.0D+00,           54.0D+00, 55.0D+00 /)

  ia(1:n+1) = (/ 1, 4, 8, 11, 15, 18 /)

  ja(1:nzmax) = (/ &
    1, 2,    4,    &
    1, 2, 3,    5, &
       2, 3, 4,    &
    1,    3, 4, 5, & 
       2, 4, 5 /)

  call csrncf ( n, a, ja, ia, nzmax, nz, coef, jcoef, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST01 - Fatal error!'
    write ( *, '(a)' ) '  CSRNCF returned IERR = ', ierr
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       I       J        A(I,J)'
  write ( *, '(a)' ) ' '

  do i = 1, nz
    write ( *, '(2x,i6,2x,i6,2x,g14.6)' ) jcoef(i,1), jcoef(i,2), coef(i)
  end do

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
