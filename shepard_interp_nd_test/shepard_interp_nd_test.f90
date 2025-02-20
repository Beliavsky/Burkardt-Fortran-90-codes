program main

!*****************************************************************************80
!
!! shepard_interp_nd_test() tests shepard_interp_nd().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 June 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: p_test_num = 4

  integer j
  integer m
  integer n1d
  integer nd
  real ( kind = rk ) p
  real ( kind = rk ), dimension ( p_test_num ) :: p_test = (/ &
    1.0D+00, 2.0D+00, 4.0D+00, 8.0D+00 /)
  integer prob
  integer prob_num

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'shepard_interp_nd_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test shepard_interp_nd().'
  write ( *, '(a)' ) '  The R8LIB library is needed.'
  write ( *, '(a)' ) '  This test also needs the TEST_INTERP_ND library.'
!
!  Look at Shepard interpolant on an irregular grid.
!
  nd = 25

  call p00_prob_num ( prob_num )

  do prob = 1, prob_num

    do m = 2, 5, 3

      do j = 1, p_test_num
        p = p_test(j)
        call test01 ( prob, p, m, nd )
      end do

    end do
  end do
!
!  Look at Shepard interpolant on a regular N1D^M grid.
!
  n1d = 5

  call p00_prob_num ( prob_num )

  do prob = 1, prob_num

    do m = 2, 5, 3

      do j = 1, p_test_num
        p = p_test(j)
        call test02 ( prob, p, m, n1d )
      end do

    end do
  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'shepard_interp_nd_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine test01 ( prob, p, m, nd )

!*****************************************************************************80
!
!! TEST01() tests SHEPARD_INTERP() on an irregular grid.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 June 2024
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer PROB, the problem number.
!
!    Input, real ( kind = rk ) P, the power used in the distance weighting.
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer ND, the number of data points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) app_error
  real ( kind = rk ), allocatable :: c(:)
  real ( kind = rk ) int_error
  integer m
  integer nd
  integer ni
  real ( kind = rk ) p
  integer prob
  real ( kind = rk ) r8vec_norm_affine
  real ( kind = rk ), allocatable :: w(:)
  real ( kind = rk ), allocatable :: xd(:,:)
  real ( kind = rk ), allocatable :: xi(:,:)
  real ( kind = rk ), allocatable :: zd(:)
  real ( kind = rk ), allocatable :: ze(:)
  real ( kind = rk ), allocatable :: zi(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST01():'
  write ( *, '(a,i4)' ) '  Interpolate data from TEST_INTERP_ND problem #', prob
  write ( *, '(a,g14.6)' ) '  using Shepard interpolation with P = ', p
  write ( *, '(a,i4)' ) '  spatial dimension M = ', m
  write ( *, '(a,i4,a)' ) '  and an irregular grid of ND = ', nd, ' data points.'
!
!  Set problem parameters:
!
  allocate ( c(1:m) )
  call r8vec_uniform_01 ( m, c )
  allocate ( w(1:m) )
  call r8vec_uniform_01 ( m, w )

  allocate ( xd(1:m,1:nd) )
  call r8mat_uniform_01 ( m, nd, xd )

  allocate ( zd(1:nd) )
  call p00_f ( prob, m, c, w, nd, xd, zd )
!
!  #1:  Does interpolant match function at interpolation points?
!
  ni = nd
  allocate ( xi(1:m,1:ni) )
  xi(1:m,1:ni) = xd(1:m,1:ni)
  allocate ( zi(1:ni) )
  call shepard_interp_nd ( m, nd, xd, zd, p, ni, xi, zi )

  int_error = r8vec_norm_affine ( ni, zi, zd ) / real ( ni, kind = rk )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) &
    '  L2 interpolation error averaged per interpolant node = ', int_error

  deallocate ( xi )
  deallocate ( zi )
!
!  #2: Approximation test.  Estimate the integral (f-interp(f))^2.
!
  ni = 1000
  ni = 50
  allocate ( xi(1:m,1:ni) )
  call r8mat_uniform_01 ( m, ni, xi )

  allocate ( zi(1:ni) )
  call shepard_interp_nd ( m, nd, xd, zd, p, ni, xi, zi )

  allocate ( ze(1:ni) )
  call p00_f ( prob, m, c, w, ni, xi, ze )

  app_error = r8vec_norm_affine ( ni, zi, ze ) / real ( ni, kind = rk )

  write ( *, '(a,g14.6)' ) &
    '  L2 approximation error averaged per 1000 samples =     ', app_error

  deallocate ( c )
  deallocate ( w )
  deallocate ( xd )
  deallocate ( xi )
  deallocate ( zd )
  deallocate ( ze )
  deallocate ( zi )

  return
end
subroutine test02 ( prob, p, m, n1d )

!*****************************************************************************80
!
!! TEST02() tests SHEPARD_INTERP_ND() on a regular N1D^M grid.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer PROB, the problem number.
!
!    Input, real ( kind = rk ) P, the power used in the distance weighting.
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N1D, the number of points in 1D.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) app_error
  real ( kind = rk ) b
  real ( kind = rk ), allocatable :: c(:)
  integer i
  real ( kind = rk ) int_error
  integer m
  integer n1d
  integer nd
  integer ni
  real ( kind = rk ) p
  integer prob
  real ( kind = rk ) r8vec_norm_affine
  real ( kind = rk ), allocatable :: w(:)
  real ( kind = rk ), allocatable :: x1d(:)
  real ( kind = rk ), allocatable :: xd(:,:)
  real ( kind = rk ), allocatable :: xi(:,:)
  real ( kind = rk ), allocatable :: zd(:)
  real ( kind = rk ), allocatable :: ze(:)
  real ( kind = rk ), allocatable :: zi(:)
!
!  Set problem parameters:
!
  allocate ( c(1:m) )
  call r8vec_uniform_01 ( m, c )
  allocate ( w(1:m) )
  call r8vec_uniform_01 ( m, w )

  nd = n1d ** m

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST02():'
  write ( *, '(a,i4)' ) '  Interpolate data from TEST_INTERP_ND problem #', prob
  write ( *, '(a,g14.6)' ) '  using Shepard interpolation with P = ', p
  write ( *, '(a,i4)' ) '  spatial dimension M = ', m
  write ( *, '(a,i6,a)' ) '  and a regular grid of N1D^M = ', nd, ' data points.'

  a = 0.0D+00
  b = 1.0D+00
  allocate ( x1d(1:n1d) )

  call r8vec_linspace ( n1d, a, b, x1d )

  allocate ( xd(1:m,1:nd) )
  do i = 1, m
    call r8vec_direct_product ( i, n1d, x1d, m, nd, xd )
  end do

  allocate ( zd(1:nd) )
  call p00_f ( prob, m, c, w, nd, xd, zd )
!
!  #1:  Does interpolant match function at interpolation points?
!
  ni = nd
  allocate ( xi(1:m,1:nd) )
  xi(1:m,1:nd) = xd(1:m,1:nd)
  allocate ( zi(1:ni) )
  call shepard_interp_nd ( m, nd, xd, zd, p, ni, xi, zi )

  int_error = r8vec_norm_affine ( ni, zi, zd ) / real ( ni, kind = rk )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) &
    '  L2 interpolation error averaged per interpolant node = ', int_error

  deallocate ( xi )
  deallocate ( zi )
!
!  #2: Approximation test.  Estimate the integral (f-interp(f))^2.
!
  ni = 1000
  allocate ( xi(1:m,1:ni) )
  call r8mat_uniform_01 ( m, ni, xi )

  allocate ( zi(1:ni) )
  call shepard_interp_nd ( m, nd, xd, zd, p, ni, xi, zi )

  allocate ( ze(1:ni) )
  call p00_f ( prob, m, c, w, ni, xi, ze )

  app_error = r8vec_norm_affine ( ni, zi, ze ) / real ( ni, kind = rk )

  write ( *, '(a,g14.6)' ) &
    '  L2 approximation error averaged per 1000 samples =     ', app_error

  deallocate ( c )
  deallocate ( w )
  deallocate ( xd )
  deallocate ( xi )
  deallocate ( zd )
  deallocate ( ze )
  deallocate ( zi )

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! timestamp() prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 8 ) ampm
  integer d
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  integer values(8)
  integer y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end

