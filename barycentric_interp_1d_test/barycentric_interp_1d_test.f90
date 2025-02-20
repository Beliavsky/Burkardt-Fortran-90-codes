program main

!*****************************************************************************80
!
!! barycentric_interp_1d_test() tests barycentric_interp_1d().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nd_test_num = 6

  integer i
  integer nd
  integer, dimension ( nd_test_num ) :: nd_test = (/ &
    4, 8, 16, 32, 64, 1000 /)
  integer prob
  integer prob_num

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'barycentric_interp_1d_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test barycentric_interp_1d().'
  write ( *, '(a)' ) '  The R8LIB library is needed.'
  write ( *, '(a)' ) '  The tests need the TEST_INTERP_1D library.'

  call p00_prob_num ( prob_num )

  do prob = 1, prob_num
    do i = 1, nd_test_num
      nd = nd_test(i)
      call lagcheby1_interp_1d_test ( prob, nd )
    end do
  end do

  do prob = 1, prob_num
    do i = 1, nd_test_num
      nd = nd_test(i)
      call lagcheby2_interp_1d_test ( prob, nd )
    end do
  end do

  do prob = 1, prob_num
    do i = 1, nd_test_num
      nd = nd_test(i)
      call lageven_interp_1d_test ( prob, nd )
    end do
  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'barycentric_interp_1d_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine lagcheby1_interp_1d_test ( prob, nd )

!*****************************************************************************80
!
!! lagcheby1_interp_1d_test() tests lagcheby1_interp_1d().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer PROB, the problem index.
!
!    integer ND, the number of data points to use.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) int_error
  integer nd
  integer ni
  integer prob
  real ( kind = rk ) r8vec_norm_affine
  real ( kind = rk ), allocatable :: xd(:)
  real ( kind = rk ), allocatable :: xi(:)
  real ( kind = rk ), allocatable :: yd(:)
  real ( kind = rk ), allocatable :: yi(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'lagcheby1_interp_1d_test():'
  write ( *, '(a)' ) '  lagcheb1_interp_1d() uses Chebyshev Type 1 spacing for data points.'
  write ( *, '(a,i6)' ) '  Interpolate data from TEST_INTERP_1D problem #', prob
  write ( *, '(a,i6)' ) '  Number of data points = ', nd
!
!  Define the data.
!
  a =  0.0D+00
  b = +1.0D+00
  allocate ( xd(1:nd) )
  allocate ( yd(1:nd) )
  call r8vec_cheby1space ( nd, a, b, xd )
  call p00_f ( prob, nd, xd, yd )

  if ( nd < 10 ) then
    call r8vec2_print ( nd, xd, yd, '  Data array:' )
  end if
!
!  #1:  Does the interpolant match the function at the interpolation points?
!
  ni = nd
  allocate ( xi(1:ni) )
  allocate ( yi(1:ni) )
  xi(1:ni) = xd(1:ni)
  call lagcheby1_interp_1d ( nd, xd, yd, ni, xi, yi )

  int_error = r8vec_norm_affine ( ni, yi, yd ) / real ( ni, kind = rk )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) &
    '  L2 interpolation error averaged per interpolant node = ', int_error

  deallocate ( xd )
  deallocate ( xi )
  deallocate ( yd )
  deallocate ( yi )

  return
end
subroutine lagcheby2_interp_1d_test ( prob, nd )

!*****************************************************************************80
!
!! lagcheby2_interp_1d_test() tests lagcheby2_interp_1d().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer PROB, the problem index.
!
!    integer ND, the number of data points to use.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) int_error
  integer nd
  integer ni
  integer prob
  real ( kind = rk ) r8vec_norm_affine
  real ( kind = rk ), allocatable :: xd(:)
  real ( kind = rk ), allocatable :: xi(:)
  real ( kind = rk ), allocatable :: yd(:)
  real ( kind = rk ), allocatable :: yi(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'lagcheby2_interp_1d_test():'
  write ( *, '(a)' ) '  lagcheby2_interp_1d() uses Chebyshev Type 2 spacing for data points.'
  write ( *, '(a,i6)' ) '  Interpolate data from TEST_INTERP_1D problem #', prob
  write ( *, '(a,i6)' ) '  Number of data points = ', nd
!
!  Define the data.
!
  a =  0.0D+00
  b = +1.0D+00
  allocate ( xd(1:nd) )
  allocate ( yd(1:nd) )
  call r8vec_cheby2space ( nd, a, b, xd )
  call p00_f ( prob, nd, xd, yd )

  if ( nd < 10 ) then
    call r8vec2_print ( nd, xd, yd, '  Data array:' )
  end if
!
!  #1:  Does the interpolant match the function at the interpolation points?
!
  ni = nd
  allocate ( xi(1:ni) )
  allocate ( yi(1:ni) )
  xi(1:ni) = xd(1:ni)
  call lagcheby2_interp_1d ( nd, xd, yd, ni, xi, yi )

  int_error = r8vec_norm_affine ( ni, yi, yd ) / real ( ni, kind = rk )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) &
    '  L2 interpolation error averaged per interpolant node = ', int_error

  deallocate ( xd )
  deallocate ( xi )
  deallocate ( yd )
  deallocate ( yi )

  return
end
subroutine lageven_interp_1d_test ( prob, nd )

!*****************************************************************************80
!
!! lageven_interp_1d_test() tests lageven_interp_1d().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer PROB, the problem index.
!
!    integer ND, the number of data points to use.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) int_error
  integer nd
  integer ni
  integer prob
  real ( kind = rk ) r8vec_norm_affine
  real ( kind = rk ), allocatable :: xd(:)
  real ( kind = rk ), allocatable :: xi(:)
  real ( kind = rk ), allocatable :: yd(:)
  real ( kind = rk ), allocatable :: yi(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'lageven_interp_1d_test():'
  write ( *, '(a)' ) '  lageven_interp_1d_test() uses even spacing for data points.'
  write ( *, '(a,i6)' ) '  Interpolate data from TEST_INTERP_1D problem #', prob
  write ( *, '(a,i6)' ) '  Number of data points = ', nd
!
!  Define the data.
!
  a =  0.0D+00
  b = +1.0D+00
  allocate ( xd(1:nd) )
  allocate ( yd(1:nd) )
  call r8vec_midspace ( nd, a, b, xd )
  call p00_f ( prob, nd, xd, yd )

  if ( nd < 10 ) then
    call r8vec2_print ( nd, xd, yd, '  Data array:' )
  end if
!
!  #1:  Does the interpolant match the function at the interpolation points?
!
  ni = nd
  allocate ( xi(1:ni) )
  allocate ( yi(1:ni) )
  xi(1:ni) = xd(1:ni)
  call lageven_interp_1d ( nd, xd, yd, ni, xi, yi )

  int_error = r8vec_norm_affine ( ni, yi, yd ) / real ( ni, kind = rk )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) &
    '  L2 interpolation error averaged per interpolant node = ', int_error

  deallocate ( xd )
  deallocate ( xi )
  deallocate ( yd )
  deallocate ( yi )

  return
end
 
