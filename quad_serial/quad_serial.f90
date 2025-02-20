function quad_serial ( f, a, b, n )

!*****************************************************************************80
!
!! quad_serial() estimates an integral using quadrature.
!
!  Discussion:
!
!    This code is a candidate for parallelization.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 October 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    function f ( x ): the function whose integral is to be estimated.
!
!    real a, b: the limits of integration.
!
!    integer n: the number of sample points.
!
!  Output:
!
!    real quad_serial: an estimate for the integral of f(x) over the interval.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  external f
  real ( kind = rk8 ) f
  integer i
  integer n
  real ( kind = rk8 ) q
  real ( kind = rk8 ) quad_serial
  real ( kind = rk8 ) x

  q = 0.0D+00
  do i = 1, n
    x = ( ( n - i ) * a + ( i - 1 ) * b ) / ( n - 1 )
    q = q + f ( x )
  end do

  q = ( b - a ) * q / dble ( n )
 
  quad_serial = q

  return
end

