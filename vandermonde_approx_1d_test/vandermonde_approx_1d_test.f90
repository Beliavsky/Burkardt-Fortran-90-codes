program main

!*****************************************************************************80
!
!! vandermonde_approx_1d_test() tests vandermonde_approx_1d().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 October 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m_test_num = 8

  integer j
  integer m
  integer, dimension(m_test_num) :: m_test = (/ &
    0, 1, 2, 3, 4, 5, 9, 12 /)
  integer prob
  integer prob_num

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'VANDERMONDE_APPROX_1D_TEST():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test VANDERMONDE_APPROX_1D().'
  write ( *, '(a)' ) '  The R8LIB library is needed.'
  write ( *, '(a)' ) '  The QR_SOLVE library is needed.'
  write ( *, '(a)' ) '  The test needs the CONDITION library.'
  write ( *, '(a)' ) '  The test needs the TEST_INTERP libary.'

  call p00_prob_num ( prob_num )
  do prob = 1, prob_num
    do j = 1, m_test_num
      m = m_test(j)
      call test01 ( prob, m )
    end do
  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'vandermonde_approx_1d_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine test01 ( prob, m )

!*****************************************************************************80
!
!! TEST01 tests VANDERMONDE_APPROX_1D_MATRIX.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer PROB, the problem number.
!
!    Input, integer M, the polynomial degree.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: a(:,:)
  real ( kind = rk ) app_error
  real ( kind = rk ), allocatable :: c(:)
  logical, parameter :: debug = .false.
  real ( kind = rk ) ld
  real ( kind = rk ) li
  integer m
  integer nd
  integer ni
  integer prob
  real ( kind = rk ) r8vec_norm_affine
  real ( kind = rk ), allocatable :: xd(:)
  real ( kind = rk ), allocatable :: xi(:)
  real ( kind = rk ) xmax
  real ( kind = rk ) xmin
  real ( kind = rk ), allocatable :: xy(:,:)
  real ( kind = rk ), allocatable :: yd(:)
  real ( kind = rk ), allocatable :: yi(:)
  real ( kind = rk ) ymax
  real ( kind = rk ) ymin

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a,i4)' ) '  Approximate data from TEST_INTERP problem #', prob

  call p00_data_num ( prob, nd )
  write ( *, '(a,i4)' ) '  Number of data points = ', nd

  allocate ( xy(1:2,1:nd) )
  call p00_data ( prob, 2, nd, xy )
  
  if ( debug ) then
    call r8mat_transpose_print ( 2, nd, xy, '  Data array:' )
  end if

  allocate ( xd(1:nd) )
  allocate ( yd(1:nd) )
  xd = xy(1,1:nd)
  yd = xy(2,1:nd)
!
!  Compute the Vandermonde matrix.
!
  write ( *, '(a,i4)' ) '  Using polynomial approximant of degree ', m

  allocate ( a(1:nd,0:m) )
  call vandermonde_approx_1d_matrix ( nd, m, xd, a )
!
!  Solve linear system.
!
  allocate ( c(0:m) )
  call qr_solve ( nd, m + 1, a, yd, c )
!
!  #1:  Does approximant match function at data points?
!
  ni = nd
  allocate ( xi(1:ni) )
  xi(1:ni) = xd(1:ni)

  allocate ( yi(1:ni) )
  call r8poly_values_horner ( m, c, ni, xi, yi )

  app_error = r8vec_norm_affine ( ni, yi, yd ) / real ( ni, kind = rk )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  L2 data approximation error = ', app_error

  deallocate ( xi )
  deallocate ( yi )
!
!  #2: Compare estimated curve length to piecewise linear (minimal) curve length.
!  Assume data is sorted, and normalize X and Y dimensions by (XMAX-XMIN) and
!  (YMAX-YMIN).
!
  xmin = minval ( xd(1:nd) )
  xmax = maxval ( xd(1:nd) )
  ymin = minval ( yd(1:nd) )
  ymax = maxval ( yd(1:nd) )

  ni = 501
  allocate ( xi(1:ni) )
  call r8vec_linspace ( ni, xmin, xmax, xi )

  allocate ( yi(1:ni) )
  call r8poly_values_horner ( m, c, ni, xi, yi )

  ld = sum ( sqrt ( ( ( xd(2:nd) - xd(1:nd-1) ) / ( xmax - xmin ) ) ** 2 &
                  + ( ( yd(2:nd) - yd(1:nd-1) ) / ( ymax - ymin ) ) ** 2 ) )

  li = sum ( sqrt ( ( ( xi(2:ni) - xi(1:ni-1) ) / ( xmax - xmin ) ) ** 2 &
                  + ( ( yi(2:ni) - yi(1:ni-1) ) / ( ymax - ymin ) ) ** 2 ) )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Normalized length of piecewise linear interpolant = ', ld
  write ( *, '(a,g14.6)' ) '  Normalized length of polynomial approximant       = ', li

  deallocate ( a )
  deallocate ( c )
  deallocate ( xd )
  deallocate ( xi )
  deallocate ( xy )
  deallocate ( yd )
  deallocate ( yi )

  return
end
subroutine r8poly_values_horner ( m, c, n, x, p )

!*****************************************************************************80
!
!! r8poly_values_horner evaluates a polynomial using Horner's method.
!
!  Discussion:
!
!    The polynomial
!
!      p(x) = c0 + c1 * x + c2 * x^2 + ... + cm * x^m
!
!    is to be evaluated at the vector of values X.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, the degree.
!
!    real ( kind = rk ) C(0:M), the polynomial coefficients.
!    C(I) is the coefficient of X^I.
!
!    integer N, the number of evaluation points.
!
!    real ( kind = rk ) X(N), the evaluation points.
!
!  Output:
!
!    real ( kind = rk ) P(N), the polynomial values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) c(0:m)
  integer i
  real ( kind = rk ) p(n)
  real ( kind = rk ) x(n)

  p(1:n) = c(m)
  do i = m - 1, 0, -1
    p(1:n) = p(1:n) * x(1:n) + c(i)
  end do

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
