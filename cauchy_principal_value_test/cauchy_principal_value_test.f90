program main

!*****************************************************************************80
!
!! cauchy_principal_value_test() tests cauchy_principal_value().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'cauchy_principal_value_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test cauchy_principal_value().'

  call cauchy_principal_value_test01 ( )
  call cauchy_principal_value_test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'cauchy_principal_value_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine cauchy_principal_value_test01 ( )

!*****************************************************************************80
!
!! cauchy_principal_value_test01(): CPV of Integral ( -1 <= t <= 1 ) exp(t) / t dt
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) cpv
  real ( kind = rk ) exact
  real ( kind = rk ), external :: f01
  integer n
  real ( kind = rk ) value

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'cauchy_principal_value_test01():'
  write ( *, '(a)' ) '  cpv() computes cauchy principal value of'
  write ( *, '(a)' ) '  Integral ( -1 <= t <= 1 ) exp(t) / t dt'

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '   N           Estimate             Error'
  write ( *, '(a)' ) ''

  exact = 2.11450175075D+00
  a = -1.0D+00
  b = +1.0D+00
  do n = 2, 8, 2
    value = cpv ( f01, a, b, n )
    write ( *, '(2x,i2,2x,g24.16,2x,g14.6)' ) n, value, abs ( value - exact )
  end do

  return
end
function f01 ( t )

!*****************************************************************************80
!
!! f01() evaluates the integrand of Integral ( -1 <= t <= 1 ) exp(t) / t dt
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) T, the argument.
!
!  Output:
!
!    real ( kind = rk ) F01, the value of the integrand.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) f01
  real ( kind = rk ) t

  f01 = exp ( t )

  return
end
subroutine cauchy_principal_value_test02 ( )

!*****************************************************************************80
!
!! cauchy_principal_value_test02() is another test.
!
!  Discussion:
!
!    We seek
!      cauchy_principal_value ( Integral ( 1-delta <= t <= 1+delta ) 1/(1-t)^3 dt )
!    which we must rewrite as
!      cauchy_principal_value ( Integral ( 1-delta <= t <= 1+delta ) 1/(1+t+t^2) 1/(1-t) dt )
!    so that our "integrand" is 1/(1+t+t^2).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) cpv
  real ( kind = rk ) delta
  real ( kind = rk ) exact
  real ( kind = rk ), external :: f02
  integer k
  integer n
  real ( kind = rk ) r1
  real ( kind = rk ) r2
  real ( kind = rk ) r3
  real ( kind = rk ) r4
  real ( kind = rk ) value

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'cauchy_principal_value_test02():'
  write ( *, '(a)' ) '  cpv() computes cauchy_principal_value of'
  write ( *, '(a)' ) '  Integral ( 1-delta <= t <= 1+delta ) 1/(1-t)^3 dt )'
  write ( *, '(a)' ) '  Try this for delta = 1, 1/2, 1/4.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) &
    '   N          Estimate                  Exact                  Error'
  delta = 1.0D+00
  do k = 1, 3
    write ( *, '(a)' ) ''
    r1 = (   delta + 1.5D+00 ) ** 2 + 0.75D+00
    r2 = ( - delta + 1.5D+00 ) ** 2 + 0.75D+00
    r3 = atan ( sqrt ( 0.75D+00 ) / (   delta + 1.5D+00 ) )
    r4 = atan ( sqrt ( 0.75D+00 ) / ( - delta + 1.5D+00 ) )
    exact = - log ( r1 / r2 ) / 6.0D+00 + ( r3 - r4 ) / sqrt ( 3.0D+00 )
    do n = 2, 8, 2
      a = 1.0D+00 - delta
      b = 1.0D+00 + delta
      value = cpv ( f02, a, b, n )
      write ( *, '(2x,i2,g24.16,2x,g24.16,2x,g14.6)' ) &
        n, value, exact, abs ( exact - value )
    end do
    delta = delta / 2.0D+00
  end do

  return
end
function f02 ( t )

!*****************************************************************************80
!
!! f02(): integrand of Integral ( 1-delta <= t <= 1+delta ) 1/(1-t^3) dt
!
!  Discussion:
!
!    1/(1-t^3) = 1/(1+t+t^2) * 1/(1-t)
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) T, the argument.
!
!  Output:
!
!    real ( kind = rk ) F02, the value of the integrand.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) f02
  real ( kind = rk ) t

  f02 = 1.0D+00 / ( 1.0D+00 + t + t * t )

  return
end
