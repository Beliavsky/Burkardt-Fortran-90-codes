program main

!*****************************************************************************80
!
!! alpert_rule_test() tests alpert_rule().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'alpert_rule_test():'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) '  Test alpert_rule().'

  call monte_carlo_regular_test ( )
  call monte_carlo_log_test ( )
  call monte_carlo_power_test ( )

  call trapezoid_regular_test ( )
  call trapezoid_log_test ( )
  call trapezoid_power_test ( )

  call alpert_regular_test ( )
  call alpert_log_test ( )
  call alpert_power_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'alpert_rule_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  return
end
subroutine alpert_log_test ( )

!*****************************************************************************80
!
!! alpert_log_test() tests the Alpert rule on the log integrand.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+0 )

  integer a_l
  integer a_r
  real ( kind = rk ), allocatable :: f1(:)
  real ( kind = rk ), allocatable :: f2(:)
  real ( kind = rk ), allocatable :: f3(:)
  real ( kind = rk ) h
  integer j_l
  integer j_r
  integer n
  integer nlog
  integer num_l
  integer order_l
  integer order_r
  integer rule
  real ( kind = rk ) s1
  real ( kind = rk ) s2
  real ( kind = rk ) s3
  real ( kind = rk ) v1
  real ( kind = rk ) v2
  real ( kind = rk ), allocatable :: w_l(:)
  real ( kind = rk ), allocatable :: w_r(:)
  real ( kind = rk ), allocatable :: x_l(:)
  real ( kind = rk ), allocatable :: x_r(:)
  real ( kind = rk ), allocatable :: x1(:)
  real ( kind = rk ), allocatable :: x2(:)
  real ( kind = rk ), allocatable :: x3(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'alpert_log_test():'
  write ( *, '(a)' ) '  Test the Alpert rule on the log integrand.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Rule  Order   J   A        N     N+2J               H        Estimate           Error'
  write ( *, '(a)' ) ''

  call integral_log ( v2 )

  call num_log ( num_l )
!
!  For the righthand interval, use the regular rule of the same index.
!
  do rule = 1, num_l

    call a_log ( rule, a_l )
    call j_log ( rule, j_l )
    call order_log ( rule, order_l )
    allocate ( x_l(1:j_l) )
    allocate ( w_l(1:j_l) )
    call rule_log ( rule, j_l, x_l, w_l )

    allocate ( f1(1:j_l) )
    allocate ( x1(1:j_l) )

    call a_regular ( rule, a_r )
    call j_regular ( rule, j_r )
    call order_regular ( rule, order_r )
    allocate ( x_r(1:j_r) )
    allocate ( w_r(1:j_r) )
    call rule_regular ( rule, j_r, x_r, w_r )

    allocate ( f3(1:j_r) )
    allocate ( x3(1:j_r) )

    n = 8

    do nlog = 4, 7

      n = n * 2
      h = 1.0D+00 / real ( n + a_l + a_r - 1, kind = rk )

      x1(1:j_l) = h * x_l(1:j_l)
      call integrand_log ( j_l, x1, f1 )
      s1 = dot_product ( w_l(1:j_l), f1(1:j_l) )

      allocate ( x2(1:n) )
      call r8vec_linspace ( n, a_l * h, ( a_l + n - 1 ) * h, x2 )
      allocate ( f2(1:n) )
      call integrand_log ( n, x2, f2 )
      s2 = sum ( f2(1:n) )
      deallocate ( f2 )
      deallocate ( x2 )

      x3(1:j_r) = 1.0D+00 - h * x_r(1:j_r)
      call integrand_log ( j_r, x3, f3 )
      s3 = dot_product ( w_r(1:j_r), f3(1:j_r) )

      v1 = h * ( s1 + s2 + s3 )

      write ( *, &
        '(4x,i2,3x,i4,2x,i2,2x,i2,2x,i7,2x,i7,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
        rule, order_l, j_l, a_l, n, n + j_l + j_r, h, v1, abs ( v1 - v2 )

    end do

    write ( *, '(a)' ) ''

    deallocate ( f1 )
    deallocate ( f3 )
    deallocate ( w_l )
    deallocate ( w_r )
    deallocate ( x_l )
    deallocate ( x_r )
    deallocate ( x1 )
    deallocate ( x3 )

  end do

  write ( *, '(a)' ) ''
  write ( *, '(51x,a,g14.6)' ) 'Exact:', v2

  return
end
subroutine alpert_power_test ( )

!*****************************************************************************80
!
!! alpert_power_test() tests the Alpert rule on the power integrand.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+0 )

  integer a_p
  integer a_r
  real ( kind = rk ), allocatable :: f1(:)
  real ( kind = rk ), allocatable :: f2(:)
  real ( kind = rk ), allocatable :: f3(:)
  real ( kind = rk ) h
  integer j_p
  integer j_r
  integer n
  integer nlog
  integer num_p
  real ( kind = rk ) order_p
  integer order_r
  integer rule
  real ( kind = rk ) s1
  real ( kind = rk ) s2
  real ( kind = rk ) s3
  real ( kind = rk ) v1
  real ( kind = rk ) v2
  real ( kind = rk ), allocatable :: w_p(:)
  real ( kind = rk ), allocatable :: w_r(:)
  real ( kind = rk ), allocatable :: x_p(:)
  real ( kind = rk ), allocatable :: x_r(:)
  real ( kind = rk ), allocatable :: x1(:)
  real ( kind = rk ), allocatable :: x2(:)
  real ( kind = rk ), allocatable :: x3(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'alpert_power_test():'
  write ( *, '(a)' ) '  Test the Alpert rule on the power integrand.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Rule  Order   J   A        N     N+2J               H        Estimate           Error'
  write ( *, '(a)' ) ''

  call integral_power ( v2 )

  call num_power ( num_p )
!
!  For the righthand interval, use the regular rule of the same index.
!
  do rule = 1, num_p

    call a_power ( rule, a_p )
    call j_power ( rule, j_p )
    call order_power ( rule, order_p )
    allocate ( x_p(1:j_p) )
    allocate ( w_p(1:j_p) )
    call rule_power ( rule, j_p, x_p, w_p )

    allocate ( f1(1:j_p) )
    allocate ( x1(1:j_p) )

    call a_regular ( rule, a_r )
    call j_regular ( rule, j_r )
    call order_regular ( rule, order_r )
    allocate ( x_r(1:j_r) )
    allocate ( w_r(1:j_r) )
    call rule_regular ( rule, j_r, x_r, w_r )

    allocate ( f3(1:j_r) )
    allocate ( x3(1:j_r) )

    n = 8

    do nlog = 4, 6

      n = n * 2
      h = 1.0D+00 / real ( n + a_p + a_r - 1, kind = rk )

      x1(1:j_p) = h * x_p(1:j_p)
      call integrand_power ( j_p, x1, f1 )
      s1 = dot_product ( w_p(1:j_p), f1(1:j_p) )

      allocate ( x2(1:n) )
      call r8vec_linspace ( n, a_p * h, ( a_p + n - 1 ) * h, x2 )
      allocate ( f2(1:n) )
      call integrand_power ( n, x2, f2 )
      s2 = sum ( f2(1:n) )
      deallocate ( f2 )
      deallocate ( x2 )

      x3(1:j_r) = 1.0D+00 - h * x_r(1:j_r)
      call integrand_power ( j_r, x3, f3 )
      s3 = dot_product ( w_r(1:j_r), f3(1:j_r) )

      v1 = h * ( s1 + s2 + s3 )

      write ( *, &
        '(4x,i2,3x,f4.1,2x,i2,2x,i2,2x,i7,2x,i7,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
        rule, order_p, j_p, a_p, n, n + j_p + j_r, h, v1, abs ( v1 - v2 )

    end do

    write ( *, '(a)' ) ''

    deallocate ( f1 )
    deallocate ( f3 )
    deallocate ( w_p )
    deallocate ( w_r )
    deallocate ( x_p )
    deallocate ( x_r )
    deallocate ( x1 )
    deallocate ( x3 )

  end do

  write ( *, '(a)' ) ''
  write ( *, '(51x,a,g14.6)' ) 'Exact:', v2

  return
end
subroutine alpert_regular_test ( )

!*****************************************************************************80
!
!! alpert_regular_test() tests the Alpert rule on the regular integrand.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+0 )

  integer a
  real ( kind = rk ), allocatable :: f1(:)
  real ( kind = rk ), allocatable :: f2(:)
  real ( kind = rk ), allocatable :: f3(:)
  real ( kind = rk ) h
  integer j
  integer n
  integer nlog
  integer num
  integer order
  integer rule
  real ( kind = rk ) s1
  real ( kind = rk ) s2
  real ( kind = rk ) s3
  real ( kind = rk ) v1
  real ( kind = rk ) v2
  real ( kind = rk ), allocatable :: w(:)
  real ( kind = rk ), allocatable :: x(:)
  real ( kind = rk ), allocatable :: x1(:)
  real ( kind = rk ), allocatable :: x2(:)
  real ( kind = rk ), allocatable :: x3(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'alpert_regular_test():'
  write ( *, '(a)' ) '  Test the Alpert rule on the regular integrand.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Rule  Order   J   A        N     N+2J               H        Estimate           Error'
  write ( *, '(a)' ) ''

  call integral_regular ( v2 )

  call num_regular ( num )

  do rule = 1, num

    call a_regular ( rule, a )
    call j_regular ( rule, j )
    call order_regular ( rule, order )
    allocate ( x(1:j) )
    allocate ( w(1:j) )
    call rule_regular ( rule, j, x, w )

    allocate ( f1(1:j) )
    allocate ( x1(1:j) )
    allocate ( f3(1:j) )
    allocate ( x3(1:j) )

    n = 8

    do nlog = 4, 6

      n = n * 2
      h = 1.0D+00 / real ( n + 2 * a - 1, kind = rk )

      x1(1:j) = h * x(1:j)
      call integrand_regular ( j, x1, f1 ) 
      s1 = dot_product ( w(1:j), f1(1:j) )

      allocate ( x2(1:n) )
      call r8vec_linspace ( n, a * h, ( a + n - 1 ) * h, x2 )
      allocate ( f2(1:n) )
      call integrand_regular ( n, x2, f2 )
      s2 = sum ( f2(1:n) )

      x3(1:j) = 1.0D+00 - h * x(1:j)
      call integrand_regular ( j, x3, f3 )
      s3 = dot_product ( w(1:j), f3(1:j) )

      v1 = h * ( s1 + s2 + s3 )

      write ( *, &
        '(4x,i2,3x,i4,2x,i2,2x,i2,2x,i7,2x,i7,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
        rule, order, j, a, n, n + 2 * j, h, v1, abs ( v1 - v2 )

      deallocate ( f2 )
      deallocate ( x2 )

    end do

    write ( *, '(a)' ) ''

    deallocate ( f1 )
    deallocate ( f3 )
    deallocate ( w )
    deallocate ( x )
    deallocate ( x1 )
    deallocate ( x3 )

  end do

  write ( *, '(a)' ) ''
  write ( *, '(51x,a,g14.6)' ) 'Exact:', v2

  return
end
subroutine monte_carlo_log_test ( )

!*****************************************************************************80
!
!! monte_carlo_log_test() tests the Monte Carlo rule on the log singular integrand.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+0 )

  real ( kind = rk ), allocatable :: f(:)
  real ( kind = rk ) h
  integer n
  integer nlog
  integer seed
  real ( kind = rk ) v1
  real ( kind = rk ) v2
  real ( kind = rk ), allocatable :: x(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'monte_carlo_log_test():'
  write ( *, '(a)' ) '  Test the Monte Carlo rule on the log singular integrand.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '          N        Estimate           Error'
  write ( *, '(a)' ) ''

  call integral_log ( v2 )

  seed = 123456789

  n = 17

  do nlog = 5, 20

    n = ( n - 1 ) * 2 + 1
    h = 1.0D+00 / real ( n, kind = rk )
    allocate ( x(1:n) )
    call random_number ( harvest = x(1:n) )
    allocate ( f(1:n) )
    call integrand_log ( n, x, f )
    v1 = h * sum ( f )
    write ( *, '(2x,i9,2x,g14.6,2x,g14.6)' ) n, v1, abs ( v1 - v2 )

    deallocate ( f )
    deallocate ( x )

  end do

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '      Exact: ', v2

  return
end
subroutine monte_carlo_power_test ( )

!*****************************************************************************80
!
!! monte_carlo_power_test() tests the Monte Carlo rule on the power singular integrand.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+0 )

  real ( kind = rk ), allocatable :: f(:)
  real ( kind = rk ) h
  integer n
  integer nlog
  integer seed
  real ( kind = rk ) v1
  real ( kind = rk ) v2
  real ( kind = rk ), allocatable :: x(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'monte_carlo_power_test():'
  write ( *, '(a)' ) '  Test the Monte Carlo rule on the power singular integrand.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '          N        Estimate           Error'
  write ( *, '(a)' ) ''

  call integral_power ( v2 )

  seed = 123456789

  n = 17

  do nlog = 5, 20

    n = ( n - 1 ) * 2 + 1
    h = 1.0D+00 / real ( n, kind = rk )
    allocate ( x(1:n) )
    call random_number ( harvest = x(1:n) )
    allocate ( f(1:n) )
    call integrand_power ( n, x, f )
    v1 = h * sum ( f )
    write ( *, '(2x,i9,2x,g14.6,2x,g14.6)' ) n, v1, abs ( v1 - v2 )

    deallocate ( f )
    deallocate ( x )

  end do

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '      Exact: ', v2

  return
end
subroutine monte_carlo_regular_test ( )

!*****************************************************************************80
!
!! monte_carlo_regular_test() tests the Monte Carlo rule on the regular integrand.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+0 )

  real ( kind = rk ), allocatable :: f(:)
  real ( kind = rk ) h
  integer n
  integer nlog
  integer seed
  real ( kind = rk ) v1
  real ( kind = rk ) v2
  real ( kind = rk ), allocatable :: x(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'monte_carlo_regular_test():'
  write ( *, '(a)' ) '  Test the Monte Carlo rule on the regular integrand.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '          N        Estimate           Error'
  write ( *, '(a)' ) ''

  call integral_regular ( v2 )

  seed = 123456789

  n = 17

  do nlog = 5, 20

    n = ( n - 1 ) * 2 + 1
    h = 1.0D+00 / real ( n, kind = rk )
    allocate ( x(1:n) )
    call random_number ( harvest = x(1:n) )
    allocate ( f(1:n) )
    call integrand_regular ( n, x, f )
    v1 = h * sum ( f )
    write ( *, '(2x,i9,2x,g14.6,2x,g14.6)' ) n, v1, abs ( v1 - v2 )

    deallocate ( f )
    deallocate ( x )

  end do

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '      Exact: ', v2

  return
end
subroutine trapezoid_log_test ( )

!*****************************************************************************80
!
!! trapezoid_log_test() tests the trapezoid rule on the log-singular integrand.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+0 )

  real ( kind = rk ), allocatable :: f(:)
  real ( kind = rk ) h
  integer n
  integer nlog
  real ( kind = rk ) v1
  real ( kind = rk ) v2
  real ( kind = rk ), allocatable :: x(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'trapezoid_log_test():'
  write ( *, '(a)' ) '  Test the trapezoidal rule on the log-singular integrand.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '        N        Estimate           Error'
  write ( *, '(a)' ) ''

  call integral_log ( v2 )

  n = 17

  do nlog = 5, 12

    n = ( n - 1 ) * 2 + 1
    h = 1.0D+00 / real ( n - 1, kind = rk )
    allocate ( x(1:n) )
    call r8vec_linspace ( n, 0.0D+00, 1.0D+00, x )
    x(1) = 0.5 * ( x(1) + x(2) )
    allocate ( f(1:n) )
    call integrand_log ( n, x, f )
    v1 = h * ( sum ( f(1:n) ) - 0.5D+00 * ( f(1) + f(n) ) )
    write ( *, '(2x,i7,2x,g14.6,2x,g14.6)' ) n, v1, abs ( v1 - v2 )

    deallocate ( f )
    deallocate ( x )

  end do

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '    Exact: ', v2

  return
end
subroutine trapezoid_power_test ( )

!*****************************************************************************80
!
!! trapezoid_power_test() tests the trapezoid rule on the power-singular integrand.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+0 )

  real ( kind = rk ), allocatable :: f(:)
  real ( kind = rk ) h
  integer n
  integer nlog
  real ( kind = rk ) v1
  real ( kind = rk ) v2
  real ( kind = rk ), allocatable :: x(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'trapezoid_power_test():'
  write ( *, '(a)' ) '  Test the trapezoidal rule on the power-singular integrand.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '        N        Estimate           Error'
  write ( *, '(a)' ) ''

  call integral_power ( v2 )

  n = 17

  do nlog = 5, 12

    n = ( n - 1 ) * 2 + 1
    h = 1.0D+00 / real ( n - 1, kind = rk )
    allocate ( x(1:n) )
    call r8vec_linspace ( n, 0.0D+00, 1.0D+00, x )
    x(1) = 0.5D+00 * ( x(1) + x(2) )
    allocate ( f(1:n) )
    call integrand_power ( n, x, f )
    v1 = h * ( sum ( f(1:n) ) - 0.5D+00 * ( f(1) + f(n) ) )
    write ( *, '(2x,i7,2x,g14.6,2x,g14.6)' ) n, v1, abs ( v1 - v2 )

    deallocate ( f )
    deallocate ( x )

  end do

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '    Exact: ', v2

  return
end
subroutine trapezoid_regular_test ( )

!*****************************************************************************80
!
!! trapezoid_regular_test() tests the trapezoid rule on the regular integrand.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+0 )

  real ( kind = rk ), allocatable :: f(:)
  real ( kind = rk ) h
  integer n
  integer nlog
  real ( kind = rk ) v1
  real ( kind = rk ) v2
  real ( kind = rk ), allocatable :: x(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'trapezoid_regular_test():'
  write ( *, '(a)' ) '  Test the trapezoidal rule on the regular integrand.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '        N        Estimate           Error'
  write ( *, '(a)' ) ''

  call integral_regular ( v2 )

  n = 17

  do nlog = 5, 12

    n = ( n - 1 ) * 2 + 1
    h = 1.0D+00 / real ( n - 1, kind = rk )
    allocate ( x(1:n) )
    call r8vec_linspace ( n, 0.0D+00, 1.0D+00, x )
    allocate ( f(1:n) )
    call integrand_regular ( n, x, f )
    v1 = h * ( sum ( f(1:n) ) - 0.5D+00 * ( f(1) + f(n) ) )
    write ( *, '(2x,i7,2x,g14.6,2x,g14.6)' ) n, v1, abs ( v1 - v2 )

    deallocate ( f )
    deallocate ( x )

  end do

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '    Exact: ', v2

  return
end

