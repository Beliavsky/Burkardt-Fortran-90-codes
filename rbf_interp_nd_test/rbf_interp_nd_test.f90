program main

!*****************************************************************************80
!
!! RBF_INTERP_ND_TEST() tests RBF_INTERP_ND().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 October 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RBF_INTERP_ND_TEST:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the RBF_INTERP_ND library.'
  write ( *, '(a)' ) '  The R8LIB library is also needed.'

  call rbf_interp_nd_test01 ( )
  call rbf_interp_nd_test02 ( )
  call rbf_interp_nd_test03 ( )
  call rbf_interp_nd_test04 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RBF_INTERP_ND_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine rbf_interp_nd_test01 ( )

!*****************************************************************************80
!
!! RBF_INTERP_ND_TEST01 tests RBF_WEIGHTS and RBF_INTERP with PHI1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 June 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 2
  integer, parameter :: n1d = 5

  integer, parameter :: nd = n1d**m

  real ( kind = rk ) a
  real ( kind = rk ) app_error
  real ( kind = rk ) b
  real ( kind = rk ) fd(nd)
  real ( kind = rk ), allocatable :: fe(:)
  real ( kind = rk ), allocatable :: fi(:)
  integer i
  real ( kind = rk ) int_error
  integer ni
  external phi1
  real ( kind = rk ) r0
  real ( kind = rk ) r8vec_norm_affine
  real ( kind = rk ) w(nd)
  real ( kind = rk ) x1d(n1d)
  real ( kind = rk ) xd(m,nd)
  real ( kind = rk ), allocatable :: xi(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RBF_INTERP_ND_TEST01:'
  write ( *, '(a)' ) '  RBF_WEIGHT computes weights for RBF interpolation.'
  write ( *, '(a)' ) '  RBF_INTERP evaluates the RBF interpolant.'
  write ( *, '(a)' ) '  Use the multiquadratic basis function PHI1(R).'

  a = 0.0D+00
  b = 2.0D+00

  call r8vec_linspace ( n1d, a, b, x1d )

  do i = 1, m
    call r8vec_direct_product ( i, n1d, x1d, m, nd, xd )
  end do

  call r8mat_transpose_print ( m, nd, xd, '  The product points:' )

  r0 = ( b - a ) / real ( n1d, kind = rk )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Scale factor R0 = ', r0

  fd(1:nd) = xd(1,1:nd) * xd(2,1:nd) * exp ( - xd(1,1:nd) * xd(2,1:nd) )

  call r8vec_print ( nd, fd, '  Function data:' )

  call rbf_weight ( m, nd, xd, r0, phi1, fd, w )

  call r8vec_print ( nd, w, '  Weight vector:' )
!
!  #1: Interpolation test.  Does interpolant match function at interpolation points?
!
  ni = nd
  allocate ( xi(1:m,1:ni) )
  allocate ( fi(1:ni) )
  xi(1:m,1:ni) = xd(1:m,1:ni)

  call rbf_interp_nd ( m, nd, xd, r0, phi1, w, ni, xi, fi )

  int_error = r8vec_norm_affine ( nd, fd, fi ) / real ( nd, kind = rk )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  L2 interpolation error averaged per interpolant node = ', int_error

  deallocate ( fi )
  deallocate ( xi )
!
!  #2: Approximation test.  Estimate the integral (f-interp(f))^2.
!
  ni = 1000
  allocate ( xi(1:m,1:ni) )
  allocate ( fi(1:ni) )
  allocate ( fe(1:ni) )
  call r8mat_uniform_ab ( m, ni, a, b, xi )
  call rbf_interp_nd ( m, nd, xd, r0, phi1, w, ni, xi, fi )

  fe(1:ni) = xi(1,1:ni) * xi(2,1:ni) * exp ( - xi(1,1:ni) * xi(2,1:ni) )
  app_error = ( b - a ) ** m * r8vec_norm_affine ( ni, fi, fe ) / real ( ni, kind = rk )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  L2 approximation error averaged per 1000 samples = ', app_error

  deallocate ( fe )
  deallocate ( fi )
  deallocate ( xi )

  return
end
subroutine rbf_interp_nd_test02 ( )

!*****************************************************************************80
!
!! RBF_INTERP_ND_TEST02 tests RBF_WEIGHTS and RBF_INTERP with PHI2.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 June 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 2
  integer, parameter :: n1d = 5

  integer, parameter :: nd = n1d**m

  real ( kind = rk ) a
  real ( kind = rk ) app_error
  real ( kind = rk ) b
  real ( kind = rk ) fd(nd)
  real ( kind = rk ), allocatable :: fe(:)
  real ( kind = rk ), allocatable :: fi(:)
  integer i
  real ( kind = rk ) int_error
  integer ni
  external phi2
  real ( kind = rk ) r0
  real ( kind = rk ) r8vec_norm_affine
  real ( kind = rk ) w(nd)
  real ( kind = rk ) x1d(n1d)
  real ( kind = rk ) xd(m,nd)
  real ( kind = rk ), allocatable :: xi(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RBF_INTERP_ND_TEST02:'
  write ( *, '(a)' ) '  RBF_WEIGHT computes weights for RBF interpolation.'
  write ( *, '(a)' ) '  RBF_INTERP evaluates the RBF interpolant.'
  write ( *, '(a)' ) '  Use the inverse multiquadratic basis function PHI2(R).'

  a = 0.0D+00
  b = 2.0D+00

  call r8vec_linspace ( n1d, a, b, x1d )

  do i = 1, m
    call r8vec_direct_product ( i, n1d, x1d, m, nd, xd )
  end do

  call r8mat_transpose_print ( m, nd, xd, '  The product points:' )

  r0 = ( b - a ) / real ( n1d, kind = rk )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Scale factor R0 = ', r0

  fd(1:nd) = xd(1,1:nd) * xd(2,1:nd) * exp ( - xd(1,1:nd) * xd(2,1:nd) )

  call r8vec_print ( nd, fd, '  Function data:' )

  call rbf_weight ( m, nd, xd, r0, phi2, fd, w )

  call r8vec_print ( nd, w, '  Weight vector:' )
!
!  #1: Interpolation test.  Does interpolant match function at interpolation points?
!
  ni = nd
  allocate ( xi(1:m,1:ni) )
  allocate ( fi(1:ni) )
  xi(1:m,1:ni) = xd(1:m,1:ni)

  call rbf_interp_nd ( m, nd, xd, r0, phi2, w, ni, xi, fi )

  int_error = r8vec_norm_affine ( nd, fd, fi ) / real ( nd, kind = rk )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  L2 interpolation error averaged per interpolant node = ', int_error

  deallocate ( fi )
  deallocate ( xi )
!
!  #2: Approximation test.  Estimate the integral (f-interp(f))^2.
!
  ni = 1000
  allocate ( xi(1:m,1:ni) )
  allocate ( fi(1:ni) )
  allocate ( fe(1:ni) )
  call r8mat_uniform_ab ( m, ni, a, b, xi )
  call rbf_interp_nd ( m, nd, xd, r0, phi2, w, ni, xi, fi )

  fe(1:ni) = xi(1,1:ni) * xi(2,1:ni) * exp ( - xi(1,1:ni) * xi(2,1:ni) )
  app_error = ( b - a ) ** m * r8vec_norm_affine ( ni, fi, fe ) / real ( ni, kind = rk )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  L2 approximation error averaged per 1000 samples = ', app_error

  deallocate ( fe )
  deallocate ( fi )
  deallocate ( xi )

  return
end
subroutine rbf_interp_nd_test03 ( )

!*****************************************************************************80
!
!! RBF_INTERP_ND_TEST03 tests RBF_WEIGHTS and RBF_INTERP with PHI3.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 June 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 2
  integer, parameter :: n1d = 5

  integer, parameter :: nd = n1d**m

  real ( kind = rk ) a
  real ( kind = rk ) app_error
  real ( kind = rk ) b
  real ( kind = rk ) fd(nd)
  real ( kind = rk ), allocatable :: fe(:)
  real ( kind = rk ), allocatable :: fi(:)
  integer i
  real ( kind = rk ) int_error
  integer ni
  external phi3
  real ( kind = rk ) r0
  real ( kind = rk ) r8vec_norm_affine
  real ( kind = rk ) w(nd)
  real ( kind = rk ) x1d(n1d)
  real ( kind = rk ) xd(m,nd)
  real ( kind = rk ), allocatable :: xi(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RBF_INTERP_ND_TEST03:'
  write ( *, '(a)' ) '  RBF_WEIGHT computes weights for RBF interpolation.'
  write ( *, '(a)' ) '  RBF_INTERP evaluates the RBF interpolant.'
  write ( *, '(a)' ) '  Use the thin-plate spline basis function PHI3(R).'

  a = 0.0D+00
  b = 2.0D+00

  call r8vec_linspace ( n1d, a, b, x1d )

  do i = 1, m
    call r8vec_direct_product ( i, n1d, x1d, m, nd, xd )
  end do

  call r8mat_transpose_print ( m, nd, xd, '  The product points:' )

  r0 = ( b - a ) / real ( n1d, kind = rk )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Scale factor R0 = ', r0

  fd(1:nd) = xd(1,1:nd) * xd(2,1:nd) * exp ( - xd(1,1:nd) * xd(2,1:nd) )

  call r8vec_print ( nd, fd, '  Function data:' )

  call rbf_weight ( m, nd, xd, r0, phi3, fd, w )

  call r8vec_print ( nd, w, '  Weight vector:' )
!
!  #1: Interpolation test.  Does interpolant match function at interpolation points?
!
  ni = nd
  allocate ( xi(1:m,1:ni) )
  allocate ( fi(1:ni) )
  xi(1:m,1:ni) = xd(1:m,1:ni)

  call rbf_interp_nd ( m, nd, xd, r0, phi3, w, ni, xi, fi )

  int_error = r8vec_norm_affine ( nd, fd, fi ) / real ( nd, kind = rk )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  L2 interpolation error averaged per interpolant node = ', int_error

  deallocate ( fi )
  deallocate ( xi )
!
!  #2: Approximation test.  Estimate the integral (f-interp(f))^2.
!
  ni = 1000
  allocate ( xi(1:m,1:ni) )
  allocate ( fi(1:ni) )
  allocate ( fe(1:ni) )
  call r8mat_uniform_ab ( m, ni, a, b, xi )
  call rbf_interp_nd ( m, nd, xd, r0, phi3, w, ni, xi, fi )

  fe(1:ni) = xi(1,1:ni) * xi(2,1:ni) * exp ( - xi(1,1:ni) * xi(2,1:ni) )
  app_error = ( b - a ) ** m * r8vec_norm_affine ( ni, fi, fe ) / real ( ni, kind = rk )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  L2 approximation error averaged per 1000 samples = ', app_error

  deallocate ( fe )
  deallocate ( fi )
  deallocate ( xi )

  return
end
subroutine rbf_interp_nd_test04 ( )

!*****************************************************************************80
!
!! RBF_INTERP_ND_TEST04 tests RBF_WEIGHTS and RBF_INTERP with PHI4.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 June 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 2
  integer, parameter :: n1d = 5

  integer, parameter :: nd = n1d**m

  real ( kind = rk ) a
  real ( kind = rk ) app_error
  real ( kind = rk ) b
  real ( kind = rk ) fd(nd)
  real ( kind = rk ), allocatable :: fe(:)
  real ( kind = rk ), allocatable :: fi(:)
  integer i
  real ( kind = rk ) int_error
  integer ni
  external phi4
  real ( kind = rk ) r0
  real ( kind = rk ) r8vec_norm_affine
  real ( kind = rk ) w(nd)
  real ( kind = rk ) x1d(n1d)
  real ( kind = rk ) xd(m,nd)
  real ( kind = rk ), allocatable :: xi(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RBF_INTERP_ND_TEST04:'
  write ( *, '(a)' ) '  RBF_WEIGHT computes weights for RBF interpolation.'
  write ( *, '(a)' ) '  RBF_INTERP evaluates the RBF interpolant.'
  write ( *, '(a)' ) '  Use the gaussian basis function PHI4(R).'

  a = 0.0D+00
  b = 2.0D+00

  call r8vec_linspace ( n1d, a, b, x1d )

  do i = 1, m
    call r8vec_direct_product ( i, n1d, x1d, m, nd, xd )
  end do

  call r8mat_transpose_print ( m, nd, xd, '  The product points:' )

  r0 = ( b - a ) / real ( n1d, kind = rk )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Scale factor R0 = ', r0

  fd(1:nd) = xd(1,1:nd) * xd(2,1:nd) * exp ( - xd(1,1:nd) * xd(2,1:nd) )

  call r8vec_print ( nd, fd, '  Function data:' )

  call rbf_weight ( m, nd, xd, r0, phi4, fd, w )

  call r8vec_print ( nd, w, '  Weight vector:' )
!
!  #1: Interpolation test.  Does interpolant match function at interpolation points?
!
  ni = nd
  allocate ( xi(1:m,1:ni) )
  allocate ( fi(1:ni) )
  xi(1:m,1:ni) = xd(1:m,1:ni)

  call rbf_interp_nd ( m, nd, xd, r0, phi4, w, ni, xi, fi )

  int_error = r8vec_norm_affine ( nd, fd, fi ) / real ( nd, kind = rk )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  L2 interpolation error averaged per interpolant node = ', int_error

  deallocate ( fi )
  deallocate ( xi )
!
!  #2: Approximation test.  Estimate the integral (f-interp(f))^2.
!
  ni = 1000
  allocate ( xi(1:m,1:ni) )
  allocate ( fi(1:ni) )
  allocate ( fe(1:ni) )
  call r8mat_uniform_ab ( m, ni, a, b, xi )
  call rbf_interp_nd ( m, nd, xd, r0, phi4, w, ni, xi, fi )

  fe(1:ni) = xi(1,1:ni) * xi(2,1:ni) * exp ( - xi(1,1:ni) * xi(2,1:ni) )
  app_error = ( b - a ) ** m * r8vec_norm_affine ( ni, fi, fe ) / real ( ni, kind = rk )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  L2 approximation error averaged per 1000 samples = ', app_error

  deallocate ( fe )
  deallocate ( fi )
  deallocate ( xi )

  return
end

