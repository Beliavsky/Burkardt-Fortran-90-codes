program main

!*****************************************************************************80
!
!! EISPACK_TEST() tests EISPACK().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 February 2023
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  logical matz

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'EISPACK_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) '  Test the EISPACK library.'

  call balanc_test ( )

  call bandv_test ( )

  call bisect_test ( )

  call bqr_test ( )

  call cbal_test ( )

  matz = .false.
  call cg_lr_test ( matz )
!
!  Do I believe the eigenvectors from CG_LR?
!
  matz = .true.
  call cg_lr_test ( matz )

  matz = .false.
  call cg_qr_test ( matz )
  matz = .true.
  call cg_qr_test ( matz )

  matz = .false.
  call ch_test ( matz )
  matz = .true.
  call ch_test ( matz )

  matz = .false.
  call ch3_test ( matz )
  matz = .true.
  call ch3_test ( matz )
!
!  This CINVIT_TEST fails.
!
  call cinvit_test ( )

  call imtqlv_test ( )
!
!  Eigenvectors from INVIT disagree with MATLAB.
!
  call invit_test ( )

  call minfit_test ( )

  matz = .false.
  call rg_elm_test ( matz )
  matz = .true.
  call rg_elm_test ( matz )

  matz = .false.
  call rg_ort_test ( matz )
  matz = .true.
  call rg_ort_test ( matz )

  matz = .false.
  call rgg_test ( matz )
  matz = .true.
  call rgg_test ( matz )

  matz = .false.
  call rs_test ( matz )
  matz = .true.
  call rs_test ( matz )

  matz = .false.
  call rsb_test ( matz )
  matz = .true.
  call rsb_test ( matz )

  matz = .false.
  call rsg_test ( matz )
  matz = .true.
  call rsg_test ( matz )

  matz = .false.
  call rsgab_test ( matz )
  matz = .true.
  call rsgab_test ( matz )

  matz = .false.
  call rsgba_test ( matz )
  matz = .true.
  call rsgba_test ( matz )

  call rsm_test ( )

  matz = .false.
  call rsp_test ( matz )
  matz = .true.
  call rsp_test ( matz )

  matz = .false.
  call rspp_test ( matz )
  matz = .true.
  call rspp_test ( matz )

  matz = .false.
  call rst_test ( matz )
  matz = .true.
  call rst_test ( matz )

  matz = .false.
  call rt_test ( matz )
  matz = .true.
  call rt_test ( matz )

  call sturm_sequence_test ( )

  call svd_test ( )

  call tql1_test ( )

  call tridib_test ( )

  call tsturm_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'EISPACK_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop
end
subroutine balanc_test ( )

!*****************************************************************************80
!
!! BALANC_TEST tests BALANC.
!
!  Discussion:
!
!    BALANC balances a matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 February 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ), dimension ( n, n ) :: a = reshape ( (/ &
    110.0D+00,  0.0D+00, 310.0D+00,  0.0D+00, 510.0D+00, &
     12.0D+00, 22.0D+00,  32.0D+00,  0.0D+00,   0.0D+00, &
     13.0D+00,  0.0D+00,  33.0D+00, 43.0D+00,  53.0D+00, &
      0.0D+00,  0.0D+00,   0.0D+00, 44.0D+00,   0.0D+00, &
     15.0D+00,  0.0D+00,  35.0D+00,  0.0D+00,  55.0D+00 /), (/ n, n /) )
  integer igh
  integer low
  real ( kind = rk ) scale(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'BALANC_TEST'
  write ( *, '(a)' ) '  BALANC balances a real general matrix.'
  write ( *, '(a,i8)' ) '  Matrix order = ', n

  call r8mat_print ( n, n, a, '  The matrix A:' )

  call balanc ( n, a, low, igh, scale )
  
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  LOW = ', low
  write ( *, '(a,i4)' ) '  IGH = ', igh

  call r8vec_print ( n, scale, '  Scaling vector SCALE:' )

  call r8mat_print ( n, n, a, '  The balanced matrix A:' )

  return
end
subroutine bandv_test ( )

!*****************************************************************************80
!
!! BANDV_TEST tests RSB.
!
!  Discussion:
!
!    BANDV computes eigenvectors of a symmetric banded matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 February 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5
  integer, parameter :: mb = 2

  real ( kind = rk ) a(n,mb)
  real ( kind = rk ) a2(n,n)
  real ( kind = rk ) e21
  integer i
  integer ierr
  integer j
  integer m
  logical matz
  real ( kind = rk ) r(n,n)
  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n,n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'BANDV_TEST'
  write ( *, '(a)' ) '  BANDV computes the eigenvectors'
  write ( *, '(a)' ) '  of a real symmetric band matrix.'
  write ( *, '(a,i8)' ) '  Matrix order = ', n
  write ( *, '(a,i8)' ) '  Half bandwidth + 1 = ', mb
!
!  A contains the band matrix in banded storage.
!
  a(1:n,1:mb) = 0.0D+00
  a(1:n,mb) = 2.0D+00
  a(2:n,1) = -1.0D+00
!
!  A2 contains the band matrix in full storage.
!
  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        a2(i,j) = 2.0D+00
      else if ( abs ( i - j ) == 1 ) then
        a2(i,j) = - 1.0D+00
      else
        a2(i,j) = 0.0D+00
      end if
    end do
  end do

  call r8mat_print ( n, n, a2, '  The matrix A:' )
!
!  Get eivenvalues from RSB.
!
  matz = .false.
  call rsb ( n, mb, a, w, matz, x, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'BANDV_TEST - Warning!'
    write ( *, '(a,i8)' ) '  RSB error return flag IERR = ', ierr
    return
  end if

  call r8vec_print ( n, w, '  The eigenvalues Lambda:' )
!
!  Now get eivenvectors from BANDV.
!
!  RSB messed up A, so we have to restore it.
!
  a(1:n,1:mb) = 0.0D+00
  a(1:n,mb) = 2.0D+00
  a(2:n,1) = -1.0D+00
  e21 = 0.0D+00
  m = n
  call bandv ( n, mb, a, e21, m, w, x, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'BANDV_TEST - Warning!'
    write ( *, '(a,i8)' ) '  BANDV error return flag IERR = ', ierr
    return
  end if
  
  call r8mat_print ( n, n, x, '  The eigenvector matrix X:' )

  r = matmul ( a2, x )

  do j = 1, n
    r(1:n,j) = r(1:n,j) - w(j) * x(1:n,j)
  end do

  call r8mat_print ( n, n, r, '  The residual (A-Lambda*I)*X:' )

  return
end
subroutine bisect_test ( )

!*****************************************************************************80
!
!! BISECT_TEST tests BISECT.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 February 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) d(n)
  real ( kind = rk ) e(n)
  real ( kind = rk ) e2(n)
  real ( kind = rk ) eps1
  integer i
  integer ierr
  integer ind(n)
  integer j
  real ( kind = rk ) lb
  integer m
  integer mm
  real ( kind = rk ) ub
  real ( kind = rk ) w(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'BISECT_TEST'
  write ( *, '(a)' ) '  BISECT computes some eigenvalues'
  write ( *, '(a)' ) '  of a real symmetric tridiagonal matrix.'
  write ( *, '(a,i8)' ) '  Matrix order = ', n

  eps1 = 100.0D+00 * epsilon ( eps1 )
!
!  Here is where the matrix is defined.
!
  d(1:n) = 2.0D+00

  e(1) = 0.0D+00
  e(2:n) = -1.0D+00
  e2(1:n) = e(1:n) ** 2
  lb = -1.0D+00
  ub = +4.0D+00
  mm = 5
!
!  We only set up and store the matrix A this way in order to make it easy
!  to compute the residual.
!
  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        a(i,j) = 2.0D+00
      else if ( abs ( i - j ) == 1 ) then
        a(i,j) = - 1.0D+00
      else
        a(i,j) = 0.0D+00
      end if
    end do
  end do

  call r8mat_print ( n, n, a, '  The matrix A:' )

  call bisect ( n, eps1, d, e, e2, lb, ub, mm, m, w, ind, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'BISECT_TEST - Warning!'
    write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
    return
  end if

  call r8vec_print ( m, w, '  The eigenvalues Lambda:' )

  return
end
subroutine bqr_test ( )

!*****************************************************************************80
!
!! BQR_TEST tests BQR.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 February 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: mb = 2
  integer, parameter :: n = 5

  real ( kind = rk ) a(n,mb)
  integer i
  integer ierr
  integer j
  integer n2
  integer nm
  real ( kind = rk ) r
  real ( kind = rk ) shift
  real ( kind = rk ) t

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'BQR_TEST'
  write ( *, '(a)' ) '  BQR computes some eigenvalues'
  write ( *, '(a)' ) '  of a real symmetric band matrix.'
  write ( *, '(a,i8)' ) '  Matrix order = ', n
  write ( *, '(a,i8)' ) '  Half bandwidth+1 = ', mb
!
!  Here is where the matrix is defined.
!
  do i = 1, n
    do j = 1, mb
      if ( j == 1 ) then
        if ( i == 1 ) then
          a(i,j) = 0.0D+00
        else
          a(i,j) = -1.0D+00
        end if
      else if ( j == 2 ) then
        a(i,j) = 2.0D+00
      end if
    end do
  end do

  call r8mat_print ( n, mb, a, '  The compressed matrix A:' )

  nm = n
  n2 = n
  t = 0.0D+00
  r = 0.0D+00

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Eigenvalues:'
  write ( *, '(a)' ) ''

  do i = 1, n

    shift = t
    call bqr ( nm, n2, mb, a, t, r, ierr )

    if ( ierr /= 0 ) then
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'BQR_TEST - Warning!'
      write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
      return
    end if

    write ( *, '(2x,i2,2x,g14.6)' ) i, t

    n2 = n2 - 1

  end do

  return
end
subroutine cbal_test ( )

!*****************************************************************************80
!
!! CBAL_TEST tests CBAL.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 March 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ), dimension ( n, n ) :: ar = reshape ( (/ &
   110.0D+00,  0.0D+00, 310.0D+00,  0.0D+00, 510.0D+00, &
    12.0D+00, 22.0D+00,  32.0D+00,  0.0D+00,   0.0D+00, &
    13.0D+00,  0.0D+00,  33.0D+00, 43.0D+00,  53.0D+00, &
     0.0D+00,  0.0D+00,   0.0D+00, 44.0D+00,   0.0D+00, &
    15.0D+00,  0.0D+00,  35.0D+00,  0.0D+00,  55.0D+00 /), (/ n, n /) )
  real ( kind = rk ), dimension ( n, n ) :: ai = reshape ( (/ &
   110.5D+00,  0.0D+00, 310.5D+00,  0.0D+00, 510.5D+00, &
    12.5D+00, 22.5D+00,  32.5D+00,  0.0D+00,   0.0D+00, &
    13.5D+00,  0.0D+00,  33.5D+00, 43.5D+00,  53.5D+00, &
     0.0D+00,  0.0D+00,   0.0D+00, 44.5D+00,   0.0D+00, &
    15.5D+00,  0.0D+00,  35.5D+00,  0.0D+00,  55.5D+00 /), (/ n, n /) )
  integer igh
  integer low
  real ( kind = rk ) scale(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'CBAL_TEST'
  write ( *, '(a)' ) '  CBAL balances a complex general matrix.'
  write ( *, '(a,i8)' ) '  Matrix order = ', n

  call r8mat_print ( n, n, ar, '  The matrix Ar:' )
  call r8mat_print ( n, n, ai, '  The matrix Ar:' )

  call cbal ( n, ar, ai, low, igh, scale )
  
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  LOW = ', low
  write ( *, '(a,i4)' ) '  IGH = ', igh

  call r8vec_print ( n, scale, '  Scaling vector SCALE:' )

  call r8mat_print ( n, n, ar, '  The balanced matrix AR:' )
  call r8mat_print ( n, n, ai, '  The balanced matrix Ai:' )

  return
end
subroutine cg_lr_test ( matz )

!*****************************************************************************80
!
!! CG_LR_TEST tests CG_LR.
!
!  Discussion:
!
!    CG_LR is for the eigenvalues of a complex general matrix.
!
!    note that the eigenvalues of such a matrix are in general complex.
!    however, we will use the same example we used before, namely
!    a hermitian matrix, so the eigenvalues will in fact be real.
!
!    (3     1     0     0+2i)
!    (1     3     0-2i  0   )
!    (0     0+2i  1     1   )
!    (0-2i  0     1     1   )
!
!    The eigenvalues are 2+2*sqrt(2), 2-2*sqrt(2), 4 and 0
!
!    The eigenvector matrix is
!
!    (  1+sqrt(2),  1,                -1,          1)
!    (  1+sqrt(2),  1,                 1,         -1)
!    (     i,       -(1+sqrt(2))*i,    i,          i)
!    (    -i,        (1+sqrt(2))*i,    i,          i)
!
!    Note that the actual eigenvector matrix from EISPACK could
!    be scaled by a real value, or by i, and the columns may
!    appear in any order.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 March
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, logical MATZ.
!    FALSE, eigenvalues only.
!    TRUE, eigenvalues and eigenvectors.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 4

  real ( kind = rk ) ai(n,n)
  real ( kind = rk ) ar(n,n)
  integer i
  integer ierr
  integer j
  logical matz
  real ( kind = rk ) wi(n)
  real ( kind = rk ) wr(n)
  real ( kind = rk ) xi(n,n)
  real ( kind = rk ) xr(n,n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'CG_LR_TEST'
  write ( *, '(a)' ) '  CG_LR computes the eigenvalues and eigenvectors of '
  write ( *, '(a)' ) '  a complex general matrix,'
  write ( *, '(a)' ) '  using elementary transformations.'
  write ( *, '(a,i8)' ) '  Matrix order = ', n

  ar = reshape ( (/ &
     3.0D+00,  1.0D+00,  0.0D+00,  0.0D+00, &
     1.0D+00,  3.0D+00,  0.0D+00,  0.0D+00, &
     0.0D+00,  0.0D+00,  1.0D+00,  1.0D+00, &
     0.0D+00,  0.0D+00,  1.0D+00,  1.0D+00 /), (/ n, n /) )
  ai = reshape ( (/ &
     0.0D+00,  0.0D+00,  0.0D+00, -2.0D+00, &
     0.0D+00,  0.0D+00,  2.0D+00,  0.0D+00, &
     0.0D+00, -2.0D+00,  0.0D+00,  0.0D+00, &
     2.0D+00,  0.0D+00,  0.0D+00,  0.0D+00 /), (/ n, n /) )
!
!  matz = .false. for eigenvalues only,
!  matz = .true. for eigenvalues and eigenvectors.
!
  call r8mat_print ( n, n, ar, 'WHAT THE FUCK AR:' )
  call r8mat_print ( n, n, ai, 'WHAT THE FUCK AI:' )
  call cg_lr ( n, ar, ai, wr, wi, matz, xr, xi, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'CG_LR_TEST - Warning!'
    write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
    return
  end if

  call r8vec2_print ( n, wr, wi, '  Real and imaginary parts of eigenvalues:' )

  if ( matz ) then
    do i = 1, n
      write ( *, '(a)' ) ''
      write ( *, '(a,i8)' ) '  Eigenvector ', i
      write ( *, '(a)' ) ''
      do j = 1, n
        write ( *, '(2g14.6)' ) xr(i,j), xi(i,j)
      end do
    end do
  end if

  return
end
subroutine cg_qr_test ( matz )

!*****************************************************************************80
!
!! CG_QR_TEST tests CG_AR.
!
!  Discussion:
!
!    CG_QR is for the eigenvalues of a complex general matrix.
!
!    The eigenvalues of such a matrix are in general complex.
!    We will use the same example we used before, namely
!    a hermitian matrix, so the eigenvalues will in fact be real.
!
!    (3     1     0     0+2i)
!    (1     3     0-2i  0   )
!    (0     0+2i  1     1   )
!    (0-2i  0     1     1   )
!
!    The eigenvalues are 2+2*sqrt(2), 2-2*sqrt(2), 4 and 0
!
!    The eigenvector matrix is
!
!    (  1+sqrt(2),  1,                -1,          1)
!    (  1+sqrt(2),  1,                 1,         -1)
!    (     i,       -(1+sqrt(2))*i,    i,          i)
!    (    -i,        (1+sqrt(2))*i,    i,          i)
!
!    Note that the actual eigenvector matrix from EISPACK could
!    be scaled by a real value, or by i, and the columns may
!    appear in any order.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 March 2018
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, logical MATZ.
!    FALSE, eigenvalues only.
!    TRUE, eigenvalues and eigenvectors.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 4

  real ( kind = rk ) ar(n,n)
  real ( kind = rk ) ai(n,n)
  integer i
  integer ierr
  integer j
  logical matz
  real ( kind = rk ) wi(n)
  real ( kind = rk ) wr(n)
  real ( kind = rk ) xi(n,n)
  real ( kind = rk ) xr(n,n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'CG_QR_TEST'
  write ( *, '(a)' ) '  CG_QR computes the eigenvalues and eigenvectors of '
  write ( *, '(a)' ) '  a complex general matrix,'
  write ( *, '(a)' ) '  using unitary transformations.'
  write ( *, '(a,i8)' ) '  Matrix order = ', n

  ar = reshape ( (/ &
     3.0D+00,  1.0D+00,  0.0D+00,  0.0D+00, &
     1.0D+00,  3.0D+00,  0.0D+00,  0.0D+00, &
     0.0D+00,  0.0D+00,  1.0D+00,  1.0D+00, &
     0.0D+00,  0.0D+00,  1.0D+00,  1.0D+00 /), (/ n, n /) )
  ai = reshape ( (/ &
     0.0D+00,  0.0D+00,  0.0D+00, -2.0D+00, &
     0.0D+00,  0.0D+00,  2.0D+00,  0.0D+00, &
     0.0D+00, -2.0D+00,  0.0D+00,  0.0D+00, &
     2.0D+00,  0.0D+00,  0.0D+00,  0.0D+00 /), (/ n, n /) )
!
!  matz = .false. for eigenvalues only,
!  matz = .true. for eigenvalues and eigenvectors.
!
  call cg_qr ( n, ar, ai, wr, wi, matz, xr, xi, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'CG_QR_TEST - Warning!'
    write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
    return
  end if

  call r8vec2_print ( n, wr, wi, '  Real and imaginary parts of eigenvalues:' )

  if ( matz ) then
    do i = 1, n
      write ( *, '(a)' ) ''
      write ( *, '(a,i8)' ) '  Eigenvector ', i
      write ( *, '(a)' ) ''
      do j = 1, n
        write ( *, '(2g14.6)' ) xr(i,j), xi(i,j)
      end do
    end do
  end if

  return
end
subroutine ch_test ( matz )

!*****************************************************************************80
!
!! CH_TEST tests CH.
!
!  Discussion:
!
!    CH is for the eigenvalues of a complex hermitian matrix.
!
!    Eigenvalues and eigenvectors of a complex hermitian matrix
!
!    Note that the eigenvalues (though not the eigenvectors) of
!    a hermitian matrix are real.
!
!    (3     1     0     0+2i)
!    (1     3     0-2i  0   )
!    (0     0+2i  1     1   )
!    (0-2i  0     1     1   )
!
!    The eigenvalues are 2+2*sqrt(2), 2-2*sqrt(2), 4 and 0
!
!    The eigenvector matrix is
!
!    (  1+sqrt(2),  1,                -1,          1)
!    (  1+sqrt(2),  1,                 1,         -1)
!    (     i,       -(1+sqrt(2))*i,    i,          i)
!    (    -i,        (1+sqrt(2))*i,    i,          i)
!
!    Note that the actual eigenvector matrix from EISPACK could
!    be scaled by a real value, or by i, and the columns may
!    appear in any order.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 February 2018
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, logical MATZ.
!    FALSE, eigenvalues only.
!    TRUE, eigenvalues and eigenvectors.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 4

  real ( kind = rk ) ar(n,n)
  real ( kind = rk ) ai(n,n)
  integer i
  integer ierr
  integer j
  logical matz
  real ( kind = rk ) w(n)
  real ( kind = rk ) xr(n,n)
  real ( kind = rk ) xi(n,n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'CH_TEST'
  write ( *, '(a)' ) '  CH computes the eigenvalues and eigenvectors of'
  write ( *, '(a)' ) '  a complex hermitian matrix.'
  write ( *, '(a,i8)' ) '  Matrix order = ', n

  ar = reshape ( (/ &
     3.0D+00,  1.0D+00,  0.0D+00,  0.0D+00, &
     1.0D+00,  3.0D+00,  0.0D+00,  0.0D+00, &
     0.0D+00,  0.0D+00,  1.0D+00,  1.0D+00, &
     0.0D+00,  0.0D+00,  1.0D+00,  1.0D+00 /), (/ n, n /) )
  ai = reshape ( (/ &
     0.0D+00,  0.0D+00,  0.0D+00, -2.0D+00, &
     0.0D+00,  0.0D+00,  2.0D+00,  0.0D+00, &
     0.0D+00, -2.0D+00,  0.0D+00,  0.0D+00, &
     2.0D+00,  0.0D+00,  0.0D+00,  0.0D+00 /), (/ n, n /) )

  call ch ( n, ar, ai, w, matz, xr, xi, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'CH_TEST - Warning!'
    write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
    return
  end if

  call r8vec_print ( n, w, '  The eigenvalues Lambda:' )

  if ( matz ) then
    do i = 1, n
      write ( *, '(a)' ) ''
      write ( *, '(a,i8)' ) '  Eigenvector ', i
      write ( *, '(a)' ) ''
      do j = 1, n
        write ( *, '(2g14.6)' ) xr(i,j), xi(i,j)
      end do
    end do
  end if

  return
end
subroutine ch3_test ( matz )

!*****************************************************************************80
!
!! CH3_TEST tests CH3.
!
!  Discussion:
!
!    CH3 is for the eigenvalues of a complex hermitian matrix.
!
!    Eigenvalues and eigenvectors of a complex hermitian matrix
!
!    Note that the eigenvalues (though not the eigenvectors) of
!    a hermitian matrix are real.
!
!    (3     1     0     0+2i)
!    (1     3     0-2i  0   )
!    (0     0+2i  1     1   )
!    (0-2i  0     1     1   )
!
!    The eigenvalues are 2+2*sqrt(2), 2-2*sqrt(2), 4 and 0
!
!    The eigenvector matrix is
!
!    (  1+sqrt(2),  1,                -1,          1)
!    (  1+sqrt(2),  1,                 1,         -1)
!    (     i,       -(1+sqrt(2))*i,    i,          i)
!    (    -i,        (1+sqrt(2))*i,    i,          i)
!
!    Note that the actual eigenvector matrix from EISPACK could
!    be scaled by a real value, or by i, and the columns may
!    appear in any order.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 March 2018
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, logical MATZ.
!    FALSE, eigenvalues only.
!    TRUE, eigenvalues and eigenvectors.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 4

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) d(n)
  integer i
  integer ierr
  integer j
  logical matz
  real ( kind = rk ) xr(n,n)
  real ( kind = rk ) xi(n,n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'CH3_TEST'
  write ( *, '(a)' ) '  CH3 computes the eigenvalues and eigenvectors of'
  write ( *, '(a)' ) '  a complex hermitian matrix.'
  write ( *, '(a,i8)' ) '  Matrix order = ', n

  a = reshape ( (/ &
     3.0D+00,  1.0D+00,  0.0D+00,  0.0D+00, &
     0.0D+00,  3.0D+00,  0.0D+00,  0.0D+00, &
     0.0D+00, -2.0D+00,  1.0D+00,  1.0D+00, &
     2.0D+00,  0.0D+00,  0.0D+00,  1.0D+00 /), (/ n, n /) )

  call r8mat_print ( n, n, a, '  Compressed matrix A:' )

  call ch3 ( n, a, d, matz, xr, xi, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'CH3_TEST - Warning!'
    write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
    return
  end if

  call r8vec_print ( n, d, '  The eigenvalues Lambda:' )

  if ( matz ) then
    do i = 1, n
      write ( *, '(a)' ) ''
      write ( *, '(a,i8)' ) '  Eigenvector ', i
      write ( *, '(a)' ) ''
      do j = 1, n
        write ( *, '(2g14.6)' ) xr(i,j), xi(i,j)
      end do
    end do
  end if

  return
end
subroutine cinvit_test ( )

!*****************************************************************************80
!
!! CINVIT_TEST tests CINVIT.
!
!  Discussion:
!
!    We seek the eigenvectors of HESS4, a complex uppper 
!    Hessenberg matrix.
!
!    Matrix A:
!
!      4+8i  7+6i  7+10i 7+10i
!      9+9i  8+1i  8+10i 2 +5i
!      0     8+3i  7+ 2i 7 +8i
!      0     0     4+10i 0 +1i
!
!    Eigenvalues:
!
!       3.324431041502838 - 2.742026572531628i
!       0.568541187348097 + 6.826204344246118i
!      -5.153228803481162 - 8.729936381660266i
!      20.260256574630240 +16.645758609945791i
!
!    Eigenvectors:
!
!     -0.3301 - 0.2223i   1.0000 + 0.0000i   0.3355 - 0.0680i   0.9522 + 0.2507i
!      1.0000 + 0.0000i   0.5032 - 0.8242i  -0.7682 + 0.0105i   1.0000 + 0.0000i
!      0.2574 + 0.3091i  -0.2150 + 0.2755i   1.0000 + 0.0000i   0.5015 - 0.1722i
!     -0.8426 + 0.1978i  -0.2382 + 0.5972i  -0.9727 - 0.1040i   0.2186 + 0.0448i
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 March 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 4

  real ( kind = rk ) ar(n,n)
  real ( kind = rk ), dimension ( n, n ) :: ar_save = reshape ( (/ &
    4.0D+00,  9.0D+00,  0.0D+00,  0.0D+00, &
    7.0D+00,  8.0D+00,  8.0D+00,  0.0D+00, &
    7.0D+00,  8.0D+00,  7.0D+00,  4.0D+00, &
    7.0D+00,  2.0D+00,  7.0D+00,  0.0D+00 /), (/ 4, 4 /) )
  real ( kind = rk ) ai(n,n)
  real ( kind = rk ), dimension ( n, n ) :: ai_save = reshape ( (/ &
    8.0D+00,  9.0D+00,  0.0D+00,  0.0D+00, &
    6.0D+00,  1.0D+00,  3.0D+00,  0.0D+00, &
   10.0D+00, 10.0D+00,  2.0D+00, 10.0D+00, &
   10.0D+00,  5.0D+00,  8.0D+00,  1.0D+00 /), (/ 4, 4 /) )
  real ( kind = rk ) fv1(n)
  real ( kind = rk ) fv2(n)
  real ( kind = rk ) fv3(n)
  integer i
  integer ierr
  integer is1
  integer is2
  integer j
  integer m
  integer mm
  real ( kind = rk ) s
  logical select(n)
  real ( kind = rk ), dimension ( n ) :: wi = (/ &
   +16.645758609945791D+00, &
    -8.729936381660266D+00, &
     6.826204344246118D+00, &
    -2.742026572531628D+00 /)
  real ( kind = rk ), dimension ( n ) :: wr = (/ &
   20.260256574630240D+00, &
   -5.153228803481162D+00, &
    3.324431041502838D+00, &
    0.568541187348097D+00 /)
  real ( kind = rk ) zi(n,n)
  real ( kind = rk ) zr(n,n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'CINVIT_TEST'
  write ( *, '(a)' ) '  CINVIT computes the eigenvectors of '
  write ( *, '(a)' ) '  a complex Hessenberg matrix.'
  write ( *, '(a,i8)' ) '  Matrix order = ', n

  ar(1:n,1:n) = ar_save(1:n,1:n)
  ai(1:n,1:n) = ai_save(1:n,1:n)

  call cbal ( n, ar, ai, is1, is2, fv1 )

  call corth ( n, is1, is2, ar, ai, fv2, fv3 )

  call comqr ( n, is1, is2, ar, ai, wr, wi, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'CINVIT_TEST - Warning!'
    write ( *, '(a)' ) '  COMQR returned flag IERR = ', ierr
    return
  end if
!
!  Now prepare to get the eigenvectors.
!
  ar(1:n,1:n) = ar_save(1:n,1:n)
  ai(1:n,1:n) = ai_save(1:n,1:n)

  call r8mat_print ( n, n, ar, '  Matrix Real Part Ar:' )
  call r8mat_print ( n, n, ai, '  Matrix Imag Part Ai:' )

  call r8vec2_print ( n, wr, wi, '  Real and imaginary parts of eigenvalues:' )

  do i = 1, n
    select(i) = .true.
  end do

  mm = n

  call cinvit ( n, ar, ai, wr, wi, select, mm, m, zr, zi, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'CINVIT_TEST - Warning!'
    write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
    return
  end if
!
!  Reverse transformations on eigenvectors.
!
  call cortb ( n, is1, is2, ar, ai, fv2, fv3, m, zr, zi )

  call cbabk2 ( n, is1, is2, fv1, m, zr, zi )

  do i = 1, n
    s = 0.0D+00
    do j = 1, mm
      s = max ( s, abs ( zr(i,j) ) );
      s = max ( s, abs ( zi(i,j) ) );
    end do

    do j = 1, mm
      zr(i,j) = zr(i,j) / s
      zi(i,j) = zi(i,j) / s
    end do
    write ( *, '(a)' ) ''
    write ( *, '(a,i8)' ) '  Eigenvector ', i
    write ( *, '(a)' ) ''
    do j = 1, mm
      write ( *, '(2g14.6)' ) zr(i,j), zi(i,j)
    end do
  end do

  return
end
subroutine imtqlv_test ( )

!*****************************************************************************80
!
!! IMTQLV_TEST tests IMTQLV.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 February 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) d(n)
  real ( kind = rk ) e(n)
  real ( kind = rk ) e2(n)
  integer i
  integer ierr
  integer ind(n)
  integer j
  real ( kind = rk ) w(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'IMTQLV_TEST'
  write ( *, '(a)' ) '  IMTQLV computes the eigenvalues of a real symmetric'
  write ( *, '(a)' ) '  tridiagonal matrix.'
  write ( *, '(a,i8)' ) '  Matrix order = ', n
!
!  Here is where the matrix is defined.
!
  d(1:n) = 2.0D+00

  e(1) = 0.0D+00
  e(2:n) = -1.0D+00

  e2(1:n) = e(1:n) ** 2
!
!  We only set up and store the matrix A this way in order to make it easy
!  to compute the residual.
!
  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        a(i,j) = 2.0D+00
      else if ( abs ( i - j ) == 1 ) then
        a(i,j) = - 1.0D+00
      else
        a(i,j) = 0.0D+00
      end if
    end do
  end do

  call r8mat_print ( n, n, a, '  The matrix A:' )

  call imtqlv ( n, d, e, e2, w, ind, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'IMTQLV_TEST - Warning!'
    write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
    return
  end if

  call r8vec_print ( n, w, '  The eigenvalues Lambda:' )

  return
end
subroutine invit_test ( )

!*****************************************************************************80
!
!! INVIT_TEST tests INVIT.
!
!  Discussion:
!
!    We seek the eigenvectors of HESS5, an uppper Hessenberg matrix.
!
!    Matrix A:
!
!     9     4     1     3     2
!     4     3     1     7     1
!     0     3     1     2     4
!     0     0     5     5     1
!     0     0     0     1     2
!
!    Eigenvalues:
!
!      1.795071645585215 + 0.000000000000000i
!     -0.484650565840867 + 3.050399870879445i
!     -0.484650565840867 - 3.050399870879445i
!      7.232089690415871 + 0.000000000000000i
!     11.942139795680633 + 0.000000000000000i
!
!    Eigenvectors:
!
! -0.4048+0.0000i -0.2788-0.1981i -0.2788+0.1981i  1.0000+0.0000i 1.0000+0.0000i
!  1.0000+0.0000i  1.0000+0.0000i  1.0000+0.0000i  0.0372+0.0000i 0.5780+0.0000i
!  0.0565+0.0000i -0.0712-0.9695i -0.0712+0.9695i -0.2064+0.0000i 0.1887+0.0000i
!  0.1687+0.0000i -0.3560+0.6933i -0.3560-0.6933i -0.5057+0.0000i 0.1379+0.0000i
! -0.8231+0.0000i  0.1938-0.0411i  0.1938+0.0411i -0.0966+0.0000i 0.0139+0.0000i
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 March 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a(n,n)
  real ( kind = rk ), dimension ( 5, 5 ), save :: a_save = reshape ( (/ &
     9.0D+00, 4.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, &
     4.0D+00, 3.0D+00, 3.0D+00, 0.0D+00, 0.0D+00, &
     1.0D+00, 1.0D+00, 1.0D+00, 5.0D+00, 0.0D+00, &
     3.0D+00, 7.0D+00, 2.0D+00, 5.0D+00, 1.0D+00, &
     2.0D+00, 1.0D+00, 4.0D+00, 1.0D+00, 2.0D+00 /), (/ 5, 5 /) )
  real ( kind = rk ) fv1(n)
  integer i
  integer ierr
  integer is1
  integer is2
  integer j
  integer k
  integer m
  integer mm
  real ( kind = rk ) ort(n)
  real ( kind = rk ) s
  real ( kind = rk ) sum1
  real ( kind = rk ) sum2
  logical select(n)
  real ( kind = rk ) wi(n)
  real ( kind = rk ) wr(n)
  real ( kind = rk ) x(n,n)
  real ( kind = rk ) xi(n,n)
  real ( kind = rk ) xr(n,n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'INVIT_TEST'
  write ( *, '(a)' ) '  INVIT computes the eigenvectors of '
  write ( *, '(a)' ) '  a real Hessenberg matrix.'
  write ( *, '(a,i8)' ) '  Matrix order = ', n
!
!  Get the eigenvalues.
!
  a(1:n,1:n) = a_save(1:n,1:n)

  call balanc ( n, a, is1, is2, fv1 )

  call orthes ( n, is1, is2, a, ort )

  call hqr ( n, is1, is2, a, wr, wi, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'INVIT_TEST - Warning!'
    write ( *, '(a)' ) '  HQR returned flag IERR = ', ierr
    return
  end if
!
!  Now prepare to get the eigenvectors.
!
  a(1:n,1:n) = a_save(1:n,1:n)

  call r8mat_print ( n, n, a, '  Matrix A:' )

  call r8vec2_print ( n, wr, wi, '  Real and imaginary parts of eigenvalues:' )

  do i = 1, n
    select(i) = .true.
  end do

  mm = n

  call invit ( n, a, wr, wi, select, mm, m, x, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'INVIT_TEST - Warning!'
    write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
    return
  end if
!
!  Reverse transformations on eigenvectors.
!
  call ortbak ( n, is1, is2, a, ort, m, x )

  call balbak ( n, is1, is2, fv1, m, x )
!
!  Expand the compressed eigenvector information.
!
  do j = 1, n

    do i = 1, n
      if ( wi(j) == 0.0D+00 ) then
        xr(i,j) = x(i,j)
        xi(i,j) = 0.0D+00
      else if ( 0.0D+00 < wi(j) ) then
        xr(i,j) = x(i,j)
        xi(i,j) = x(i,j+1)
      else if ( wi(j) < 0.0D+00 ) then
        xr(i,j) = x(i,j-1)
        xi(i,j) = -x(i,j)
      end if
    end do
!
!  Normalize the eigenvectors.
!
    s = 0.0D+00
    do i = 1, n
      s = max ( s, abs ( xr(i,j) ) )
      s = max ( s, abs ( xi(i,j) ) )
    end do

    xr(1:n,j) = xr(1:n,j) / s
    xi(1:n,j) = xi(1:n,j) / s

    call r8vec2_print ( n, xr(1:n,j), xi(1:n,j), '  Eigenvector:'  )

  end do
!
!  Check.
!  First, restore the original values of A.
!
  a(1:n,1:n) = a_save(1:n,1:n)

  do k = 1, n

    write ( *, '(a)' ) ''
    write ( *, '(a,i8)' ) '  Residuals (A*x-Lambda*x) for eigenvalue ', k
    write ( *, '(a)' ) ''

    do i = 1, n
      sum1 = 0.0D+00
      sum2 = 0.0D+00
      do j = 1, n
        sum1 = sum1 + a(i,j) * xr(j,k)
        sum2 = sum2 + a(i,j) * xi(j,k)
      end do
      sum1 = sum1 - wr(k) * xr(i,k) + wi(k) * xi(i,k)
      sum2 = sum2 - wi(k) * xr(i,k) - wr(k) * xi(i,k)
      write ( *, '(2g14.6)' ) sum1, sum2
    end do

  end do

  return
end
subroutine minfit_test ( )

!*****************************************************************************80
!
!! MINFIT_TEST tests MINFIT.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: nb = 1
  integer, parameter :: n = 2

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) acopy(m,n)
  real ( kind = rk ) b(m,nb)
  integer i
  integer ierr
  integer j
  real ( kind = rk ) r(m)
  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'MINFIT_TEST'
  write ( *, '(a)' ) '  MINFIT solves an overdetermined linear system'
  write ( *, '(a)' ) '  using least squares methods.'
  write ( *, '(a,i8)' ) '  Matrix rows = ', m
  write ( *, '(a,i8)' ) '  Matrix columns = ', n

  a(1,1) =   1.00D+00
  a(2,1) =   2.05D+00
  a(3,1) =   3.06D+00
  a(4,1) = - 1.02D+00
  a(5,1) =   4.08D+00

  a(1,2) =   1.00D+00
  a(2,2) = - 1.00D+00
  a(3,2) =   1.00D+00
  a(4,2) =   2.00D+00
  a(5,2) = - 1.00D+00

  acopy(1:m,1:n) = a(1:m,1:n)

  b(1,1) = 1.98D+00
  b(2,1) = 0.95D+00
  b(3,1) = 3.98D+00
  b(4,1) = 0.92D+00
  b(5,1) = 2.90D+00

  r(1:m) = - b(1:m,1)

  call r8mat_print ( m, n, a, '  The matrix A:' )

  call r8mat_print ( m, nb, b, '  The right hand side B:' )

  call minfit ( m, m, n, a, w, nb, b, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'MINFIT_TEST - Warning!'
    write ( *, '(a,i8)' ) '  MINFIT error code IERR = ', ierr
    return
  end if

  call r8vec_print ( n, w, '  The singular values:' )
!
!  B now contains U' * B.
!  We need to divide by the singular values, and multiply by V.
!
  b(1:n,1) = b(1:n,1) / w(1:n)

  do i = 1, n
    x(i) = 0.0D+00
    do j = 1, n
      x(i) = x(i) + a(i,j) * b(j,1)
    end do
  end do

  call r8vec_print ( n, x, '  The least squares solution X:' )

  do i = 1, m
    do j = 1, n
      r(i) = r(i) + acopy(i,j) * x(j)
    end do
  end do

  call r8vec_print ( m, r, '  The residual A * X - B:' )

  return
end
subroutine rg_elm_test ( matz )

!*****************************************************************************80
!
!! RG_ELM_TEST tests RG_ELM.
!
!  Discussion:
!
!    RG is for the eigenvalues of a general real matrix.
!
!    The matrix A is nonsymmetric.  The eigenvalues may therefore be
!    complex numbers.
!
!    ( 33  16  72)
!    (-24 -10 -57)
!    ( -8  -4 -17)
!
!    The eigenvalues of A are (1,2,3)
!
!    The eigenvectors of A are
!
!    (-1.0000 -1,0000  1.0000 )
!    ( 0.8000  0.8125 -0.7500 )
!    ( 0.2667  0.2500 -0.2500 )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 March 2018
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, logical MATZ.
!    FALSE, eigenvalues only.
!    TRUE, eigenvalues and eigenvectors.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 3

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) :: a_save(n,n) = reshape ( (/ &
    33.0D+00, -24.0D+00,  -8.0D+00, &
    16.0D+00, -10.0D+00,  -4.0D+00, &
    72.0D+00, -57.0D+00, -17.0D+00 /), (/ n, n /) )
  integer i
  integer ierr
  integer j
  integer k
  logical matz
  real ( kind = rk ) s
  real ( kind = rk ) sum1
  real ( kind = rk ) sum2
  real ( kind = rk ) wi(n)
  real ( kind = rk ) wr(n)
  real ( kind = rk ) x(n,n)
  real ( kind = rk ) xi(n,n)
  real ( kind = rk ) xr(n,n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'RG_ELM_TEST'
  write ( *, '(a)' ) '  RG_ELM computes the eigenvalues and eigenvectors of'
  write ( *, '(a)' ) '  a real general matrix, using elementary transformations.'
  write ( *, '(a,i8)' ) '  Matrix order = ', n
!
!  Set the values of the matrix.
!
  a(1:n,1:n) = a_save(1:n,1:n)

  call r8mat_print ( n, n, a, '  The matrix A:' )

  call rg_elm ( n, a, wr, wi, matz, x, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'RG_ELM_TEST - Warning!'
    write ( *, '(a,i8)' ) '  The error return flag was IERR = ', ierr
    return
  end if

  call r8vec2_print ( n, wr, wi, '  Real and imaginary parts of eigenvalues:' )

  if ( matz ) then
!
!  Expand the compressed eigenvector information.
!
    do j = 1, n

      do i = 1, n
        if ( wi(j) == 0.0D+00 ) then
          xr(i,j) = x(i,j)
          xi(i,j) = 0.0D+00
        else if ( 0.0D+00 < wi(j) ) then
          xr(i,j) = x(i,j)
          xi(i,j) = x(i,j+1)
        else if ( wi(j) < 0.0D+00 ) then
          xr(i,j) = x(i,j-1)
          xi(i,j) = -x(i,j)
        end if
      end do
!
!  Normalize the eigenvectors.
!
      s = 0.0D+00
      do i = 1, n
        s = max ( s, abs ( xr(i,j) ) )
        s = max ( s, abs ( xi(i,j) ) )
      end do

      xr(1:n,j) = xr(1:n,j) / s
      xi(1:n,j) = xi(1:n,j) / s

      call r8vec2_print ( n, xr(1:n,j), xi(1:n,j), '  Eigenvector:'  )

    end do
!
!  Check.
!  First, restore the original values of A.
!
    a(1:n,1:n) = a_save(1:n,1:n)

    do k = 1, n

      write ( *, '(a)' ) ''
      write ( *, '(a,i8)' ) '  Residuals (A*x-Lambda*x) for eigenvalue ', k
      write ( *, '(a)' ) ''

      do i = 1, n
        sum1 = 0.0D+00
        sum2 = 0.0D+00
        do j = 1, n
          sum1 = sum1 + a(i,j) * xr(j,k)
          sum2 = sum2 + a(i,j) * xi(j,k)
        end do
        sum1 = sum1 - wr(k) * xr(i,k) + wi(k) * xi(i,k)
        sum2 = sum2 - wi(k) * xr(i,k) - wr(k) * xi(i,k)
        write ( *, '(2g14.6)' ) sum1, sum2
      end do

    end do

  end if

  return
end
subroutine rg_ort_test ( matz )

!*****************************************************************************80
!
!! RG_ORT_TEST tests RG_ORT.
!
!  Discussion:
!
!    RG_ORT is for the eigenvalues of a general real matrix.
!
!    The matrix A is nonsymmetric.  The eigenvalues may therefore be
!    complex numbers.
!
!    ( 33  16  72)
!    (-24 -10 -57)
!    ( -8  -4 -17)
!
!    The eigenvalues of A are (1,2,3)
!
!    The eigenvectors of A are
!
!    (-1.0000 -1,0000  1.0000 )
!    ( 0.8000  0.8125 -0.7500 )
!    ( 0.2667  0.2500 -0.2500 )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 March 2018
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, logical MATZ.
!    FALSE, eigenvalues only.
!    TRUE, eigenvalues and eigenvectors.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 3

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) :: a_save(n,n) = reshape ( (/ &
    33.0D+00, -24.0D+00,  -8.0D+00, &
    16.0D+00, -10.0D+00,  -4.0D+00, &
    72.0D+00, -57.0D+00, -17.0D+00 /), (/ n, n /) )
  integer i
  integer ierr
  integer j
  integer k
  logical matz
  real ( kind = rk ) s
  real ( kind = rk ) sum1
  real ( kind = rk ) sum2
  real ( kind = rk ) wi(n)
  real ( kind = rk ) wr(n)
  real ( kind = rk ) x(n,n)
  real ( kind = rk ) xi(n,n)
  real ( kind = rk ) xr(n,n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'RG_ORT_TEST'
  write ( *, '(a)' ) '  RG_ORT computes the eigenvalues and eigenvectors of'
  write ( *, '(a)' ) '  a real general matrix, using orthogonal transformations.'
  write ( *, '(a,i8)' ) '  Matrix order = ', n
!
!  Set the matrix.
!
  a(1:n,1:n) = a_save(1:n,1:n)

  call r8mat_print ( n, n, a, '  The matrix A:' )

  call rg_ort ( n, a, wr, wi, matz, x, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'RG_ORT_TEST - Warning!'
    write ( *, '(a,i8)' ) '  The error return flag was IERR = ', ierr
    return
  end if

  call r8vec2_print ( n, wr, wi, '  Real and imaginary parts of eigenvalues:' )

  if ( matz ) then
!
!  Expand the compressed eigenvector information.
!
    do j = 1, n

      do i = 1, n
        if ( wi(j) == 0.0D+00 ) then
          xr(i,j) = x(i,j)
          xi(i,j) = 0.0D+00
        else if ( 0.0D+00 < wi(j) ) then
          xr(i,j) = x(i,j)
          xi(i,j) = x(i,j+1)
        else if ( wi(j) < 0.0D+00 ) then
          xr(i,j) = x(i,j-1)
          xi(i,j) = -x(i,j)
        end if
      end do
!
!  Normalize the eigenvectors.
!
      s = 0.0D+00
      do i = 1, n
        s = max ( s, abs ( xr(i,j) ) )
        s = max ( s, abs ( xi(i,j) ) )
      end do

      xr(1:n,j) = xr(1:n,j) / s
      xi(1:n,j) = xi(1:n,j) / s

      call r8vec2_print ( n, xr(1:n,j), xi(1:n,j), '  Eigenvector:'  )

    end do
!
!  Check.
!  First, restore the original values of A.
!
    a(1:n,1:n) = a_save(1:n,1:n)

    do k = 1, n

      write ( *, '(a)' ) ''
      write ( *, '(a,i8)' ) '  Residuals (A*x-Lambda*x) for eigenvalue ', k
      write ( *, '(a)' ) ''

      do i = 1, n
        sum1 = 0.0D+00
        sum2 = 0.0D+00
        do j = 1, n
          sum1 = sum1 + a(i,j) * xr(j,k)
          sum2 = sum2 + a(i,j) * xi(j,k)
        end do
        sum1 = sum1 - wr(k) * xr(i,k) + wi(k) * xi(i,k)
        sum2 = sum2 - wi(k) * xr(i,k) - wr(k) * xi(i,k)
        write ( *, '(2g14.6)' ) sum1, sum2
      end do

    end do

  end if

  return
end
subroutine rgg_test ( matz )

!*****************************************************************************80
!
!! RGG_TEST tests RGG.
!
!  Discussion:
!
!    RGG is for a real generalized general eigenvalue problem.
!
!    A generalized eigenvalue problem.  Given matrices A and B, find
!    N numbers LAMBDA, and for each LAMBDA a vector X, so that
!
!      A*x = lambda*B*x
!
!    The matrix A is
!
!    ( -7 7  6  6)
!    (-10 8 10  8)
!    ( -8 3 10 11)
!    ( -4 0  4 12)
!
!    The matrix B is
!
!    (2 1 0 0)
!    (1 2 1 0)
!    (0 1 2 1)
!    (0 0 1 2)
!
!    The correct eigenvalues LAMBDA are
!
!    (1,2,3,4)
!
!    The correct eigenvectors X are
!
!    (4 3 2 1)
!    (3 3 2 1)
!    (2 2 2 1)
!    (1 1 1 1)
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 February 2018
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, logical MATZ.
!    FALSE, eigenvalues only.
!    TRUE, eigenvalues and eigenvectors.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 4

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) acopy(n,n)
  real ( kind = rk ) alfi(n)
  real ( kind = rk ) alfr(n)
  real ( kind = rk ) b(n,n)
  real ( kind = rk ) bcopy(n,n)
  real ( kind = rk ) beta(n)
  integer i
  integer ierr
  integer j
  integer k
  logical matz
  real ( kind = rk ) sum3
  real ( kind = rk ) sum1
  real ( kind = rk ) sum2
  real ( kind = rk ) x(n,n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'RGG_TEST:'
  write ( *, '(a)' ) '  RGG for real generalized problem.'
  write ( *, '(a)' ) '  Find scalars LAMBDA and vectors X so that'
  write ( *, '(a)' ) '    A*x = LAMBDA * B * x'
  write ( *, '(a,i8)' ) '  Matrix order = ', n
!
!  Set the values in the A matrix.
!
  a(1,1) = -7.0D+00
  a(1,2) = 7.0D+00
  a(1,3) = 6.0D+00
  a(1,4) = 6.0D+00

  a(2,1) = -10.0D+00
  a(2,2) = 8.0D+00
  a(2,3) = 10.0D+00
  a(2,4) = 8.0D+00

  a(3,1) = -8.0D+00
  a(3,2) = 3.0D+00
  a(3,3) = 10.0D+00
  a(3,4) = 11.0D+00

  a(4,1) = -4.0D+00
  a(4,2) = 0.0D+00
  a(4,3) = 4.0D+00
  a(4,4) = 12.0D+00
!
!  Save a copy of A.
!
  acopy(1:n,1:n) = a(1:n,1:n)
!
!  Set the values in the B matrix.
!
  b(1,1) = 2.0D+00
  b(1,2) = 1.0D+00
  b(1,3) = 0.0D+00
  b(1,4) = 0.0D+00

  b(2,1) = 1.0D+00
  b(2,2) = 2.0D+00
  b(2,3) = 1.0D+00
  b(2,4) = 0.0D+00

  b(3,1) = 0.0D+00
  b(3,2) = 1.0D+00
  b(3,3) = 2.0D+00
  b(3,4) = 1.0D+00

  b(4,1) = 0.0D+00
  b(4,2) = 0.0D+00
  b(4,3) = 1.0D+00
  b(4,4) = 2.0D+00
!
!  Save a copy of B.
!
  bcopy(1:n,1:n) = b(1:n,1:n)

  call r8mat_print ( n, n, a, '  The matrix A:' )

  call r8mat_print ( n, n, b, '  The matrix B:' )

  call rgg ( n, a, b, alfr, alfi, beta, matz, x, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'RGG_TEST - Warning!'
    write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
    return
  end if

  alfr(1:n) = alfr(1:n) / beta(1:n)
  alfi(1:n) = alfi(1:n) / beta(1:n)

  call r8vec2_print ( n, alfr, alfi, '  Real and imaginary parts of eigenvalues:' )

  if ( matz ) then
    do i = 1, n
      write ( *, '(a)' ) ''
      write ( *, '(a,i8)' ) '  Eigenvector ', i
      write ( *, '(a)' ) ''
      do j = 1, n
        write ( *, '(g14.6)' ) x(i,j)
      end do
    end do
  end if
!
!  Check.
!  First, restore the original values of A and B.
!
  if ( matz ) then

    a(1:n,1:n) = acopy(1:n,1:n)
    b(1:n,1:n) = bcopy(1:n,1:n)

    do k = 1, n
      write ( *, '(a)' ) ''
      write ( *, '(a,i8)' ) &
        '  Residuals (A*x-(Alfr+Alfi*I)*B*x) for eigenvalue ', k
      write ( *, '(a)' ) ''

      if ( alfi(k) == 0.0D+00 ) then

        do i = 1, n

          sum3 = dot_product ( a(i,1:n), x(1:n,k) )

          do j = 1, n
            sum3 = sum3 - alfr(k) * b(i,j) * x(j,k)
          end do

          write ( *, '(g14.6)' ) sum3
        end do

      else if ( 0.0D+00 < alfi(k) ) then

        do i = 1, n

          sum1 = 0.0D+00
          sum2 = 0.0D+00
          do j = 1, n
            sum1 = sum1 + a(i,j) * x(j,k)
            sum2 = sum2 + a(i,j) * x(j,k+1)
          end do

          do j = 1, n
            sum1 = sum1 - alfr(k) * b(i,j) * x(j,k) &
                        + alfi(k) * b(i,j) * x(j,k+1)

            sum2 = sum2 - alfi(k) * b(i,j) * x(j,k) &
                        - alfr(k) * b(i,j) * x(j,k+1)

          end do

          write ( *, '(2g14.6)' ) sum1, sum2
        end do

      else if ( alfi(k) < 0.0D+00 ) then

        do i = 1, n

          sum1 = 0.0D+00
          sum2 = 0.0D+00
          do j = 1, n
            sum1 = sum1 + a(i,j) * x(j,k-1)
            sum2 = sum2 - a(i,j) * x(j,k)
          end do

          do j = 1, n

            sum1 = sum1 - alfr(k) * b(i,j) * x(j,k-1) &
                        - alfi(k) * b(i,j) * x(j,k)

            sum2 = sum2 - alfi(k) * b(i,j) * x(j,k-1) &
                        + alfr(k) * b(i,j) * x(j,k)
          end do

          write ( *, '(2g14.6)' ) sum1, sum2
        end do

      end if

    end do

  end if

  return
end
subroutine rs_test ( matz )

!*****************************************************************************80
!
!! RS_TEST tests RS.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 February 2018
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, logical MATZ.
!    FALSE, eigenvalues only.
!    TRUE, eigenvalues and eigenvectors.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 4

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) a2(n,n)
  integer ierr
  integer j
  logical matz
  real ( kind = rk ) r(n,n)
  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n,n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'RS_TEST'
  write ( *, '(a)' ) '  RS computes the eigenvalues and eigenvectors'
  write ( *, '(a)' ) '  of a real symmetric matrix.'
  write ( *, '(a,i8)' ) '  Matrix order = ', n
!
!  Set the values in the matrix.
!
  a(1,1) = 5.0D+00
  a(1,2) = 4.0D+00
  a(1,3) = 1.0D+00
  a(1,4) = 1.0D+00

  a(2,1) = 4.0D+00
  a(2,2) = 5.0D+00
  a(2,3) = 1.0D+00
  a(2,4) = 1.0D+00

  a(3,1) = 1.0D+00
  a(3,2) = 1.0D+00
  a(3,3) = 4.0D+00
  a(3,4) = 2.0D+00

  a(4,1) = 1.0D+00
  a(4,2) = 1.0D+00
  a(4,3) = 2.0D+00
  a(4,4) = 4.0D+00
!
!  Save a copy of the matrix.
!
  a2(1:n,1:n) = a(1:n,1:n)

  call r8mat_print ( n, n, a, '  The matrix A:' )

  call rs ( n, a, w, matz, x, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'RS_TEST - Warning!'
    write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
    return
  end if

  call r8vec_print ( n, w, '  The eigenvalues Lambda:' )

  if ( matz ) then

    call r8mat_print ( n, n, x, '  The eigenvector matrix:' )

    r(1:n,1:n) = matmul ( a2(1:n,1:n), x(1:n,1:n) )

    do j = 1, n
      r(1:n,j) = r(1:n,j) - w(j) * x(1:n,j)
    end do

    call r8mat_print ( n, n, r, '  The residual (A-Lambda*I)*X:' )

  end if

  return
end
subroutine rsb_test ( matz )

!*****************************************************************************80
!
!! RSB_TEST tests RSB.
!
!  Discussion:
!
!    RSB solves the eigenvalue problem for a symmetric banded matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 February 2018
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, logical MATZ.
!    FALSE, eigenvalues only.
!    TRUE, eigenvalues and eigenvectors.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5
  integer, parameter :: mb = 2

  real ( kind = rk ) a(n,mb)
  real ( kind = rk ) a2(n,n)
  integer i
  integer ierr
  integer j
  logical matz
  real ( kind = rk ) r(n,n)
  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n,n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'RSB_TEST'
  write ( *, '(a)' ) '  RSB computes the eigenvalues and eigenvectors'
  write ( *, '(a)' ) '  of a real symmetric band matrix.'
  write ( *, '(a,i8)' ) '  Matrix order = ', n
!
!  A contains the band matrix in banded storage.
!
  a(1:n,1:mb) = 0.0D+00
  a(1:n,mb) = 2.0D+00
  a(2:n,1) = -1.0D+00
!
!  A2 contains the band matrix in full storage.
!
  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        a2(i,j) = 2.0D+00
      else if ( abs ( i - j ) == 1 ) then
        a2(i,j) = - 1.0D+00
      else
        a2(i,j) = 0.0D+00
      end if
    end do
  end do

  call r8mat_print ( n, n, a2, '  The matrix A:' )

  call rsb ( n, mb, a, w, matz, x, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'RSB_TEST - Warning!'
    write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
    return
  end if

  call r8vec_print ( n, w, '  The eigenvalues Lambda:' )

  if ( matz ) then

    call r8mat_print ( n, n, x, '  The eigenvector matrix X:' )

    r = matmul ( a2, x )

    do j = 1, n
      r(1:n,j) = r(1:n,j) - w(j) * x(1:n,j)
    end do

    call r8mat_print ( n, n, r, '  The residual (A-Lambda*I)*X:' )

  end if

  return
end
subroutine rsg_test ( matz )

!*****************************************************************************80
!
!! RSG_TEST tests RSG.
!
!  Discussion:
!
!    RSG is for a real generalized eigenvalue problem of the form
!
!      A*x = lambda*B*x
!
!    with A symmetric and B positive definite symmetric.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 February 2018
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, logical MATZ.
!    FALSE, eigenvalues only.
!    TRUE, eigenvalues and eigenvectors.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 4

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) a2(n,n)
  real ( kind = rk ) b(n,n)
  real ( kind = rk ) b2(n,n)
  integer i
  integer ierr
  integer j
  integer k
  logical matz
  real ( kind = rk ) sum3
  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n,n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'RSG_TEST:'
  write ( *, '(a)' ) '  RSG for real symmetric generalized problem.'
  write ( *, '(a)' ) '  Find scalars LAMBDA and vectors X so that'
  write ( *, '(a)' ) '    A*x = LAMBDA * B * x'
  write ( *, '(a,i8)' ) '  Matrix order = ', n

  do i = 1, n
    do j = 1, n
      a(i,j) = abs ( i - j )
    end do
  end do

  a2(1:n,1:n) = a(1:n,1:n)

  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        b(i,j) = 2.0D+00
      else if ( abs ( i - j ) == 1 ) then
        b(i,j) = - 1.0D+00
      else
        b(i,j) = 0.0D+00
      end if
    end do
  end do

  b2(1:n,1:n) = b(1:n,1:n)

  call r8mat_print ( n, n, a, '  The matrix A:' )

  call r8mat_print ( n, n, b, '  The matrix B:' )

  call rsg ( n, a, b, w, matz, x, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'RSG_TEST - Warning!'
    write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
    return
  end if

  call r8vec_print ( n, w, '  The eigenvalues Lambda:' )

  if ( matz ) then

    call r8mat_print ( n, n, x, '  The eigenvector matrix X:' )

    a(1:n,1:n) = a2(1:n,1:n)
    b(1:n,1:n) = b2(1:n,1:n)

    do k = 1, n

      write ( *, '(a)' ) ''
      write ( *, '(a,i8)' ) &
        '  Residuals (A*x-(w*I)*B*x) for eigenvalue ', k
      write ( *, '(a)' ) ''

        do i = 1, n

          sum3 = 0.0D+00
          do j = 1, n
            sum3 = sum3 + a(i,j) * x(j,k)
          end do

          do j = 1, n
            sum3 = sum3 - w(k) * b(i,j) * x(j,k)
          end do

          write ( *, '(g14.6)' ) sum3
        end do

    end do

  end if

  return
end
subroutine rsgab_test ( matz )

!*****************************************************************************80
!
!! RSGAB_TEST tests RSGAB.
!
!  Discussion:
!
!    RSGAB is for a real generalized eigenvalue problem of the form
!
!      A*B*x = lambda*x
!
!    with A symmetric and B positive definite symmetric.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 February 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, logical MATZ.
!    FALSE, eigenvalues only.
!    TRUE, eigenvalues and eigenvectors.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 4

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) a2(n,n)
  real ( kind = rk ) b(n,n)
  real ( kind = rk ) b2(n,n)
  integer i
  integer ierr
  integer j
  logical matz
  real ( kind = rk ) r(n,n)
  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n,n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'RSGAB_TEST:'
  write ( *, '(a)' ) '  RSGAB for real symmetric generalized problem.'
  write ( *, '(a)' ) '  Find scalars LAMBDA and vectors X so that'
  write ( *, '(a)' ) '    A*B*X = LAMBDA * X'
  write ( *, '(a,i8)' ) '  Matrix order = ', n

  do i = 1, n
    do j = 1, n
      a(i,j) = abs ( i - j )
    end do
  end do

  a2(1:n,1:n) = a(1:n,1:n)

  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        b(i,j) = 2.0D+00
      else if ( abs ( i - j ) == 1 ) then
        b(i,j) = - 1.0D+00
      else
        b(i,j) = 0.0D+00
      end if
    end do
  end do

  b2(1:n,1:n) = b(1:n,1:n)

  call r8mat_print ( n, n, a, '  The matrix A:' )

  call r8mat_print ( n, n, b, '  The matrix B:' )

  call rsgab ( n, a, b, w, matz, x, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'RSGAB_TEST - Warning!'
    write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
    return
  end if

  call r8vec_print ( n, w, '  The eigenvalues Lambda:' )

  if ( matz ) then

    call r8mat_print ( n, n, x, '  The eigenvector matrix X:' )

    r = matmul ( b2, x )

    r = matmul ( a2, r )

    do j = 1, n
      r(1:n,j) = r(1:n,j) - w(j) * x(1:n,j)
    end do

    call r8mat_print ( n, n, r, '  The residual matrix (A*B-Lambda*I)*X:' )

  end if

  return
end
subroutine rsgba_test ( matz )

!*****************************************************************************80
!
!! RSGBA_TEST tests RSGBA.
!
!  Discussion:
!
!    RSGBA is for a real generalized eigenvalue problem of the form
!
!      B * A * x = lambda * x
!
!    with A symmetric and B positive definite symmetric.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 February 2018
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, logical MATZ.
!    FALSE, eigenvalues only.
!    TRUE, eigenvalues and eigenvectors.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 4

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) a2(n,n)
  real ( kind = rk ) b(n,n)
  real ( kind = rk ) b2(n,n)
  integer i
  integer ierr
  integer j
  logical matz
  real ( kind = rk ) r(n,n)
  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n,n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'RSGBA_TEST:'
  write ( *, '(a)' ) '  RSGBA for real symmetric generalized problem.'
  write ( *, '(a)' ) '  Find scalars LAMBDA and vectors X so that'
  write ( *, '(a)' ) '    B*A*X = LAMBDA * X'
  write ( *, '(a,i8)' ) '  Matrix order = ', n

  do i = 1, n
    do j = 1, n
      a(i,j) = abs ( i - j )
    end do
  end do

  a2(1:n,1:n) = a(1:n,1:n)

  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        b(i,j) = 2.0D+00
      else if ( abs ( i - j ) == 1 ) then
        b(i,j) = - 1.0D+00
      else
        b(i,j) = 0.0D+00
      end if
    end do
  end do

  b2(1:n,1:n) = b(1:n,1:n)

  call r8mat_print ( n, n, a, '  The matrix A:' )

  call r8mat_print ( n, n, b, '  The matrix B:' )

  call rsgba ( n, a, b, w, matz, x, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'RSGBA_TEST - Warning!'
    write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
    return
  end if

  call r8vec_print ( n, w, '  The eigenvalues Lambda:' )

  if ( matz ) then

    call r8mat_print ( n, n, x, '  The eigenvector matrix X:' )

    r = matmul ( a2, x )

    r = matmul ( b2, r )

    do j = 1, n
      r(1:n,j) = r(1:n,j) - w(j) * x(1:n,j)
    end do

    call r8mat_print ( n, n, r, '  The residual matrix (B*A-Lambda*I)*X:' )

  end if

  return
end
subroutine rsm_test ( )

!*****************************************************************************80
!
!! RSM_TEST tests RSM.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 February 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 4
  integer, parameter :: m = n

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) a2(n,n)
  integer ierr
  integer j
  real ( kind = rk ) r(n,m)
  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n,m)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'RSM_TEST'
  write ( *, '(a)' ) '  RSM computes some eigenvalues and eigenvectors'
  write ( *, '(a)' ) '  of a real symmetric matrix.'
  write ( *, '(a,i8)' ) '  Matrix order = ', n
  write ( *, '(a,i8)' ) '  Number of eigenvectors desired = ', m

  a(1,1) = 5.0D+00
  a(1,2) = 4.0D+00
  a(1,3) = 1.0D+00
  a(1,4) = 1.0D+00

  a(2,1) = 4.0D+00
  a(2,2) = 5.0D+00
  a(2,3) = 1.0D+00
  a(2,4) = 1.0D+00

  a(3,1) = 1.0D+00
  a(3,2) = 1.0D+00
  a(3,3) = 4.0D+00
  a(3,4) = 2.0D+00

  a(4,1) = 1.0D+00
  a(4,2) = 1.0D+00
  a(4,3) = 2.0D+00
  a(4,4) = 4.0D+00

  a2(1:n,1:n) = a(1:n,1:n)

  call r8mat_print ( n, n, a, '  The matrix A:' )

  call rsm ( n, a, w, m, x, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'RSM_TEST - Warning!'
    write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
    return
  end if

  call r8vec_print ( n, w, '  The eigenvalues Lambda:' )

  if ( 0 < m ) then

    call r8mat_print ( n, m, x, '  The eigenvector matrix X:' )

    r = matmul ( a2, x )

    do j = 1, m
      r(1:n,j) = r(1:n,j) - w(j) * x(1:n,j)
    end do

    call r8mat_print ( n, m, r, '  The residual (A-Lambda*I)*X:' )

  end if

  return
end
subroutine rsp_test ( matz )

!*****************************************************************************80
!
!! RSP_TEST tests RSP.
!
!  Discussion:
!
!    RSP is for the eigenvalues of a real symmetric packed matrix.
!
!    A is symmetric.  Because of this, we know that the eigenvalues
!    of A must be real (rather than complex) numbers.
!
!    The entries of A are
!
!    (5 4 1 1)
!    (4 5 1 1)
!    (1 1 4 2)
!    (1 1 2 4)
!
!    The eigenvalues of A are (10, 5, 2, 1)
!
!    One set of eigenvectors of A is:
!
!    ( 2 -1  0 -1)
!    ( 2 -1  0  1)
!    ( 1  2 -1  0)
!    ( 1  2  1  0)
!
!    However, this set is not orthonormal, and EISPACK will compute
!    a different set of values.
!
!    Note that the I-th eigenvector corresponding to the I-th eigenvalue
!    consists of the I-th column of the above matrix of eigenvectors.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 February 2018
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, logical MATZ.
!    FALSE, eigenvalues only.
!    TRUE, eigenvalues and eigenvectors.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 4
  integer, parameter :: nv = ( n * ( n + 1 ) ) / 2

  real ( kind = rk ) a(nv)
  real ( kind = rk ) a2(n,n)
  integer ierr
  integer j
  logical matz
  real ( kind = rk ) r(n,n)
  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n,n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'RSP_TEST'
  write ( *, '(a)' ) '  RSP computes the eigenvalues and eigenvectors'
  write ( *, '(a)' ) '  of a real symmetric packed matrix.'
  write ( *, '(a,i8)' ) '  Matrix order = ', n
!
!  Set the values in the matrix.
!
  a(1) = 5.0D+00

  a(2) = 4.0D+00
  a(3) = 5.0D+00

  a(4) = 1.0D+00
  a(5) = 1.0D+00
  a(6) = 4.0D+00

  a(7) = 1.0D+00
  a(8) = 1.0D+00
  a(9) = 2.0D+00
  a(10) = 4.0D+00

  a2(1,1) = 5.0D+00
  a2(1,2) = 4.0D+00
  a2(1,3) = 1.0D+00
  a2(1,4) = 1.0D+00

  a2(2,1) = 4.0D+00
  a2(2,2) = 5.0D+00
  a2(2,3) = 1.0D+00
  a2(2,4) = 1.0D+00

  a2(3,1) = 1.0D+00
  a2(3,2) = 1.0D+00
  a2(3,3) = 4.0D+00
  a2(3,4) = 2.0D+00

  a2(4,1) = 1.0D+00
  a2(4,2) = 1.0D+00
  a2(4,3) = 2.0D+00
  a2(4,4) = 4.0D+00

  call r8mat_print ( n, n, a2, '  The matrix A:' )

  call rsp ( n, nv, a, w, matz, x, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'RSP_TEST - Warning!'
    write ( *, '(a,i8)' ) '  The error return flag was IERR = ', ierr
    return
  end if

  call r8vec_print ( n, w, '  The eigenvalues Lambda:' )

  if ( matz ) then

    call r8mat_print ( n, n, x, '  The eigenvector matrix X:' )

    r = matmul ( a2, x )

    do j = 1, n
      r(1:n,j) = r(1:n,j) - w(j) * x(1:n,j)
    end do

    call r8mat_print ( n, n, r, '  The residual matrix (A-Lambda*I)*X:' )

  end if

  return
end
subroutine rspp_test ( matz )

!*****************************************************************************80
!
!! RSPP_TEST tests RSPP.
!
!  Discussion:
!
!    RSPP is for some eigenvalues of a real symmetric packed matrix.
!
!    A is symmetric.  Because of this, we know that the eigenvalues
!    of A must be real (rather than complex) numbers.
!
!
!    The entries of A are
!
!    (5 4 1 1)
!    (4 5 1 1)
!    (1 1 4 2)
!    (1 1 2 4)
!
!    The eigenvalues of A are (10, 5, 2, 1)
!
!    One set of eigenvectors of A is:
!
!    ( 2 -1  0 -1)
!    ( 2 -1  0  1)
!    ( 1  2 -1  0)
!    ( 1  2  1  0)
!
!    However, this set is not orthonormal, and EISPACK will compute
!    a different set of values.
!
!    Note that the I-th eigenvector corresponding to the I-th eigenvalue
!    consists of the I-th column of the above matrix of eigenvectors.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 February 2018
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, logical MATZ.
!    FALSE, eigenvalues only.
!    TRUE, eigenvalues and eigenvectors.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 4
  integer, parameter :: m = n
  integer, parameter :: nv = ( n * ( n + 1 ) ) / 2

  real ( kind = rk ) a(nv)
  real ( kind = rk ) a2(n,n)
  integer ierr
  integer j
  logical matz
  real ( kind = rk ) r(n,m)
  logical type
  real ( kind = rk ) w(m)
  real ( kind = rk ) x(n,m)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'RSPP_TEST'
  write ( *, '(a)' ) '  RSPP finds some eigenvalues and eigenvectors of'
  write ( *, '(a)' ) '  a real symmetric packed matrix.'
  write ( *, '(a,i8)' ) '  Matrix order = ', n

!
!  Set the values in the matrix.
!
  a(1) = 5.0D+00

  a(2) = 4.0D+00
  a(3) = 5.0D+00

  a(4) = 1.0D+00
  a(5) = 1.0D+00
  a(6) = 4.0D+00

  a(7) = 1.0D+00
  a(8) = 1.0D+00
  a(9) = 2.0D+00
  a(10) = 4.0D+00

  a2(1,1) = 5.0D+00
  a2(1,2) = 4.0D+00
  a2(1,3) = 1.0D+00
  a2(1,4) = 1.0D+00

  a2(2,1) = 4.0D+00
  a2(2,2) = 5.0D+00
  a2(2,3) = 1.0D+00
  a2(2,4) = 1.0D+00

  a2(3,1) = 1.0D+00
  a2(3,2) = 1.0D+00
  a2(3,3) = 4.0D+00
  a2(3,4) = 2.0D+00

  a2(4,1) = 1.0D+00
  a2(4,2) = 1.0D+00
  a2(4,3) = 2.0D+00
  a2(4,4) = 4.0D+00

  call r8mat_print ( n, n, a2, '  The matrix A:' )
!
!  TYPE = TRUE to find smallest eigenvalues, FALSE for largest.
!
  type = .true.

  call rspp ( n, nv, a, w, matz, x, m, type, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'RSPP_TEST - Warning!'
    write ( *, '(a,i8)' ) '  The error return flag was IERR = ', ierr
    return
  end if

  call r8vec_print ( m, w, '  The eigenvalues Lambda:' )

  if ( matz ) then

    call r8mat_print ( n, m, x, '  The eigenvector matrix X:' )

    r = matmul ( a2, x )

    do j = 1, m
      r(1:n,j) = r(1:n,j) - w(j) * x(1:n,j)
    end do

    call r8mat_print ( n, m, r, '  The residual matrix (A-Lambda*I)*X:' )

  end if

  return
end
subroutine rst_test ( matz )

!*****************************************************************************80
!
!! RST_TEST tests RST.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 February 2018
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, logical MATZ.
!    FALSE, eigenvalues only.
!    TRUE, eigenvalues and eigenvectors.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) e(n)
  integer i
  integer ierr
  integer j
  logical matz
  real ( kind = rk ) r(n,n)
  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n,n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'RST_TEST'
  write ( *, '(a)' ) '  RST computes the eigenvalues and eigenvectors'
  write ( *, '(a)' ) '  of a real symmetric tridiagonal matrix.'
  write ( *, '(a,i8)' ) '  Matrix order = ', n
!
!  Here is where the matrix is defined.
!
  w(1:n) = 2.0D+00

  e(1) = 0.0D+00
  e(2:n) = -1.0D+00
!
!  We only set up and store the matrix A this way in order to make it easy
!  to compute the residual.
!
  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        a(i,j) = 2.0D+00
      else if ( abs ( i - j ) == 1 ) then
        a(i,j) = - 1.0D+00
      else
        a(i,j) = 0.0D+00
      end if
    end do
  end do

  call r8mat_print ( n, n, a, '  The matrix A:' )

  call rst ( n, w, e, matz, x, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'RST_TEST - Warning!'
    write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
    return
  end if

  call r8vec_print ( n, w, '  The eigenvalues Lambda:' )

  if ( matz ) then

    call r8mat_print ( n, n, x, '  The eigenvector matrix X:' )

    r = matmul ( a, x )

    do j = 1, n
      r(1:n,j) = r(1:n,j) - w(j) * x(1:n,j)
    end do

    call r8mat_print ( n, n, r, '  The residual matrix (A-Lambda*I)*X:' )

  end if

  return
end
subroutine rt_test ( matz )

!*****************************************************************************80
!
!! RT_TEST tests RT.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 February 2018
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, logical MATZ.
!    FALSE, eigenvalues only.
!    TRUE, eigenvalues and eigenvectors.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a(n,3)
  real ( kind = rk ) a2(n,n)
  integer i
  integer ierr
  integer j
  logical matz
  real ( kind = rk ) r(n,n)
  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n,n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'RT_TEST'
  write ( *, '(a)' ) '  RT computes the eigenvalues and eigenvectors'
  write ( *, '(a)' ) '  of a real sign-symmetric tridiagonal matrix.'
  write ( *, '(a,i8)' ) '  Matrix order = ', n
!
!  Here is where the matrix is defined.
!
  a(2:n,  1) = - 1.0D+00
  a(1:n,  2) =   2.0D+00
  a(1:n-1,3) = - 1.0D+00
!
!  We only set up and store the matrix A this way in order to make it easy
!  to compute the residual.
!
  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        a2(i,j) = 2.0D+00
      else if ( abs ( i - j ) == 1 ) then
        a2(i,j) = - 1.0D+00
      else
        a2(i,j) = 0.0D+00
      end if
    end do
  end do

  call r8mat_print ( n, n, a2, '  The matrix A:' )

  call rt ( n, a, w, matz, x, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'RT_TEST - Warning!'
    write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
    return
  end if

  call r8vec_print ( n, w, '  The eigenvalues Lambda:' )

  if ( matz ) then

    call r8mat_print ( n, n, x, '  The eigenvector matrix X:' )

    r = matmul ( a2, x )

    do j = 1, n
      r(1:n,j) = r(1:n,j) - w(j) * x(1:n,j)
    end do

    call r8mat_print ( n, n, r, '  The residual matrix (A-Lambda*I)*X:' )

  end if

  return
end
subroutine sturm_sequence_test ( )

!*****************************************************************************80
!
!! STURM_SEQUENCE_TEST tests STURM_SEQUENCE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 February 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) d(n)
  real ( kind = rk ) e(n)
  real ( kind = rk ) e2(n)
  integer i
  integer k
  integer p
  integer q
  integer sturm_sequence
  real ( kind = rk ) x1

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'STURM_SEQUENCE_TEST'
  write ( *, '(a)' ) '  STURM_SEQUENCE considers a (P,Q) submatrix of a'
  write ( *, '(a)' ) '  symmetric tridiagonal submatrix.'
  write ( *, '(a)' ) '  It counts the number of eigenvalues less than X1.'
  write ( *, '(a,i8)' ) '  Matrix order = ', n

  d(1:n) = 2.0D+00
  e(1) = 0.0D+00
  e(2:n) = -1.0D+00
  e2(1:n) = e(1:n) ** 2

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '     P     Q      X1      Count'
  write ( *, '(a)' ) ''

  p = 1
  q = n
  x1 = 4.5D+00

  do i = 1, 6
    k = sturm_sequence ( d, e, e2, n, p, q, x1 )
    write ( *, '(2x,i4,2x,i4,2x,g14.6,2x,i4)' ) p, q, x1, k
    x1 = x1 - 1.0D+00
  end do

  return
end
subroutine svd_test ( )

!*****************************************************************************80
!
!! SVD_TEST tests SVD.
!
!  Discussion:
!
!    In our special example, the matrix is square and symmetric.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 4
  integer, parameter :: m = n

  real ( kind = rk ) a(m,n)
  integer ierr
  integer j
  logical matu
  logical matv
  real ( kind = rk ) r(m,n)
  real ( kind = rk ) u(m,n)
  real ( kind = rk ) v(n,n)
  real ( kind = rk ) w(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'SVD_TEST'
  write ( *, '(a)' ) '  SVD computes the singular value decomposition'
  write ( *, '(a)' ) '  of a real general matrix.'
  write ( *, '(a,i8)' ) '  Matrix order = ', n
!
!  Set the values of the matrix.
!
  a(1,1) = 0.9900D+00
  a(1,2) = 0.0020D+00
  a(1,3) = 0.0060D+00
  a(1,4) = 0.0020D+00

  a(2,1) = 0.0020D+00
  a(2,2) = 0.9900D+00
  a(2,3) = 0.0020D+00
  a(2,4) = 0.0060D+00

  a(3,1) = 0.0060D+00
  a(3,2) = 0.0020D+00
  a(3,3) = 0.9900D+00
  a(3,4) = 0.0020D+00

  a(4,1) = 0.0020D+00
  a(4,2) = 0.0060D+00
  a(4,3) = 0.0020D+00
  a(4,4) = 0.9900D+00

  call r8mat_print ( m, n, a, '  The matrix A:' )

  matu = .true.
  matv = .true.

  call svd ( m, n, a, w, matu, u, matv, v, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'SVD_TEST - Warning!'
    write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
    return
  end if

  call r8vec_print ( n, w, '  The singular values S' )

  call r8mat_print ( m, n, u, '  The U matrix:' )

  call r8mat_print ( n, n, v, '  The V matrix:' )

  do j = 1, n
    v(1:n,j) = w(j) * v(1:n,j)
  end do

  r(1:m,1:n) = matmul ( u(1:m,1:n), transpose ( v(1:n,1:n) ) )

  call r8mat_print ( m, n, r, '  The product U * S * Transpose(V):' )

  return
end
subroutine tql1_test ( )

!*****************************************************************************80
!
!! TQL1_TEST tests TQL1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) d(n)
  real ( kind = rk ) e(n)
  integer i
  integer ierr
  integer j

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TQL1_TEST'
  write ( *, '(a)' ) '  TQL1 computes the eigenvalues'
  write ( *, '(a)' ) '  of a real symmetric tridiagonal matrix.'
  write ( *, '(a,i8)' ) '  Matrix order = ', n
!
!  Here is where the matrix is defined.
!
  d(1:n) = 2.0D+00

  e(1) = 0.0D+00
  e(2:n) = -1.0D+00
!
!  We only set up and store the matrix A this way in order to make it easy
!  to compute the residual.
!
  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        a(i,j) = d(i)
      else if ( i == j - 1 ) then
        a(i,j) = e(j)
      else if ( i == j + 1 ) then
        a(i,j) = e(i)
      else
        a(i,j) = 0.0D+00
      end if
    end do
  end do

  call r8mat_print ( n, n, a, '  The matrix A:' )

  call tql1 ( n, d, e, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'TQL1_TEST - Warning!'
    write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
    return
  end if

  call r8vec_print ( n, d, '  The eigenvalues Lambda:' )

  return
end
subroutine tridib_test ( )

!*****************************************************************************80
!
!! TRIDIB_TEST tests TRIDIB.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 March 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) d(n)
  real ( kind = rk ) e(n)
  real ( kind = rk ) e2(n)
  real ( kind = rk ) eps1
  integer i
  integer ierr
  integer ind(n)
  integer j
  real ( kind = rk ) lb
  integer m
  integer m11
  real ( kind = rk ) ub
  real ( kind = rk ) w(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TRIDIB_TEST'
  write ( *, '(a)' ) '  TRIDIB computes some eigenvalues'
  write ( *, '(a)' ) '  of a real symmetric tridiagonal matrix.'
  write ( *, '(a,i8)' ) '  Matrix order = ', n

  eps1 = 100.0D+00 * epsilon ( eps1 )
!
!  Here is where the matrix is defined.
!
  d(1:n) = 2.0D+00

  e(1) = 0.0D+00
  e(2:n) = -1.0D+00
  e2(1:n) = e(1:n) ** 2
  m11 = 1
  m = 5
!
!  We only set up and store the matrix A this way in order to make it easy
!  to compute the residual.
!
  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        a(i,j) = 2.0D+00
      else if ( abs ( i - j ) == 1 ) then
        a(i,j) = - 1.0D+00
      else
        a(i,j) = 0.0D+00
      end if
    end do
  end do

  call r8mat_print ( n, n, a, '  The matrix A:' )

  call tridib ( n, eps1, d, e, e2, lb, ub, m11, m, w, ind, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'TRIDIB_TEST - Warning!'
    write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
    return
  end if

  call r8vec_print ( m, w, '  The eigenvalues Lambda:' )

  return
end
subroutine tsturm_test ( )

!*****************************************************************************80
!
!! TSTURM_TEST tests TSTURM.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 March 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) d(n)
  real ( kind = rk ) e(n)
  real ( kind = rk ) e2(n)
  real ( kind = rk ) eps1
  integer i
  integer ierr
  integer j
  real ( kind = rk ) lb
  integer m
  integer mm
  real ( kind = rk ) r(n,n)
  real ( kind = rk ) ub
  real ( kind = rk ) w(n)
  real ( kind = rk ) z(n,n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TSTURM_TEST'
  write ( *, '(a)' ) '  TSTURM computes some eigenvalues and eigenvectors'
  write ( *, '(a)' ) '  of a real symmetric tridiagonal matrix.'
  write ( *, '(a,i8)' ) '  Matrix order = ', n
!
!  The effect of EPS1 seems COUNTERINTUITIVE.
!
  eps1 = 10.0D+00 * epsilon ( eps1 )
!
!  Here is where the matrix is defined.
!
  d(1:n) = 2.0D+00

  e(1) = 0.0D+00
  e(2:n) = -1.0D+00
  e2(1:n) = e(1:n) ** 2
  lb = -1.0D+00
  ub = +4.0D+00
  mm = 5
!
!  We only set up and store the matrix A this way in order to make it easy
!  to compute the residual.
!
  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        a(i,j) = 2.0D+00
      else if ( abs ( i - j ) == 1 ) then
        a(i,j) = - 1.0D+00
      else
        a(i,j) = 0.0D+00
      end if
    end do
  end do

  call r8mat_print ( n, n, a, '  The matrix A:' )

  call tred1 ( n, a, d, e, e2 )
  call tsturm ( n, eps1, d, e, e2, lb, ub, mm, m, w, z, ierr )
  call trbak1 ( n, a, e, m, z )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'TSTURM_TEST - Warning!'
    write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
    return
  end if

  call r8vec_print ( m, w, '  The eigenvalues Lambda:' )

  call r8mat_print ( n, m, z, '  The eigenvector matrix Z:' )
!
!  Restore matrix A.
!
  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        a(i,j) = 2.0D+00
      else if ( abs ( i - j ) == 1 ) then
        a(i,j) = - 1.0D+00
      else
        a(i,j) = 0.0D+00
      end if
    end do
  end do

  r = matmul ( a, z )

  do j = 1, m
    r(1:n,j) = r(1:n,j) - w(j) * z(1:n,j)
  end do

  call r8mat_print ( n, m, r, '  The residual matrix (A-Lambda*I)*X:' )

  return
end

