program main

!*****************************************************************************80
!
!! chebyshev_interp_1d_test() tests chebyshev_interp_1d().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer prob
  integer prob_num

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'CHEBYSHEV_INTERP_1D_TEST:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the CHEBYSHEV_INTERP_1D library.'
  write ( *, '(a)' ) '  The QR_SOLVE and R8LIB libraries are needed.'
  write ( *, '(a)' ) '  The test needs the TEST_INTERP library as well.'

  call p00_prob_num ( prob_num )
  do prob = 1, prob_num
    call test01 ( prob )
  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'CHEBYSHEV_INTERP_1D_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine test01 ( prob )

!*****************************************************************************80
!
!! TEST01 tests CHEBYSHEV_VALUE_1D.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) int_error
  integer nd
  integer ni
  integer prob
  real ( kind = rk ) r8vec_norm_affine
  real ( kind = rk ), allocatable :: xd(:)
  real ( kind = rk ), allocatable :: xi(:)
  real ( kind = rk ), allocatable :: xy(:,:)
  real ( kind = rk ), allocatable :: yd(:)
  real ( kind = rk ), allocatable :: yi(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'CHEBYSHEV_INTERP_1D_TEST01:'
  write ( *, '(a,i6)' ) '  Interpolate data from TEST_INTERP problem #', prob

  call p00_data_num ( prob, nd )
  write ( *, '(a,i6)' ) '  Number of data points = ', nd

  allocate ( xy(1:2,1:nd) )

  call p00_data ( prob, 2, nd, xy )
  
  call r8mat_transpose_print ( 2, nd, xy, '  Data array:' )

  allocate ( xd(1:nd) )
  allocate ( yd(1:nd) )
  
  xd(1:nd) = xy(1,1:nd)
  yd(1:nd) = xy(2,1:nd)
!
!  #1:  Does interpolant match function at interpolation points?
!
  ni = nd

  allocate ( xi(1:ni) )
  allocate ( yi(1:ni) )

  xi(1:ni) = xd(1:nd)
  call chebyshev_interp_1d ( nd, xd, yd, ni, xi, yi )

  int_error = r8vec_norm_affine ( ni, yi, yd ) / real ( ni, kind = rk )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) &
    '  L2 interpolation error averaged per interpolant node = ', int_error

  deallocate ( xd )
  deallocate ( xi )
  deallocate ( xy )
  deallocate ( yd )
  deallocate ( yi )

  return
end
 
