program main

!*****************************************************************************80
!
!! gegenbauer_cc_test() tests gegenbauer_cc().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    08 January 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'gegenbauer_cc_test():'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) '  Test gegenbauer_cc().'

  call chebyshev_even1_test ( )
  call chebyshev_even2_test ( )
  call gegenbauer_cc1_test ( )
  call gegenbauer_cc2_test ( )
  call i4_uniform_ab_test ( )
  call r8_mop_test ( )
  call r8vec_print_test ( )
  call r8vec2_print_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'GEGENBAUER_CC_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
function f ( x )

!*****************************************************************************80
!
!! F is the function to be integrated.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, the argument.
!
!    Output, real ( kind = rk ) F, the value.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) f
  real ( kind = rk ) x

  a = 2.0D+00
  f = cos ( a * x )

  return
end
subroutine chebyshev_even1_test ( )

!*****************************************************************************80
!
!! CHEBYSHEV_EVEN1_TEST tests CHEBYSHEV_EVEN1.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    08 January 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) a2(0:3)
  real ( kind = rk ) :: a2_exact(0:3) = (/ &
    0.4477815660D+00, &
   -0.7056685603D+00, &
    0.0680357987D+00, &
   -0.0048097159D+00 /)
  real ( kind = rk ), external :: f
  real ( kind = rk ) lambda
  integer s
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'CHEBYSHEV_EVEN1_TEST:'
  write ( *, '(a)' ) '  CHEBYSHEV_EVEN1 computes the even Chebyshev coefficients'
  write ( *, '(a)' ) '  of a function F, using the extreme points of Tn(x).'

  lambda = 0.75D+00
  a = 2.0D+00
  n = 6

  call chebyshev_even1 ( n, f, a2 )

  s = ( n / 2 )
  call r8vec2_print ( s + 1, a2, a2_exact, '  Computed and Exact Coefficients:' )

  return
end
subroutine chebyshev_even2_test ( )

!*****************************************************************************80
!
!! CHEBYSHEV_EVEN2_TEST tests CHEBYSHEV_EVEN2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    08 January 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b2(0:3)
  real ( kind = rk ), external :: f
  real ( kind = rk ) lambda
  integer n
  integer s

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'CHEBYSHEV_EVEN2_TEST:'
  write ( *, '(a)' ) '  CHEBYSHEV_EVEN2 computes the even Chebyshev coefficients'
  write ( *, '(a)' ) '  of a function F, using the zeros of Tn(x).'

  lambda = 0.75D+00
  a = 2.0D+00
  n = 6

  call chebyshev_even2 ( n, f, b2 )

  s = ( n / 2 )
  call r8vec_print ( s + 1, b2, '  Computed Coefficients:' )

  return
end
subroutine gegenbauer_cc1_test ( )

!*****************************************************************************80
!
!! GEGENBAUER_CC1_TEST tests GEGENBAUER_CC1.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) besselj
  real ( kind = rk ) exact
  real ( kind = rk ), external :: f
  real ( kind = rk ) lambda
  integer n
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) value

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'GEGENBAUER_CC1_TEST:'
  write ( *, '(a)' ) '  GEGENBAUER_CC1 estimates the Gegenbauer integral of'
  write ( *, '(a)' ) '  a function f(x) using a Clenshaw-Curtis type approach'
  write ( *, '(a)' ) '  based on the extreme points of Tn(x).'

  lambda = 0.75D+00
  a = 2.0D+00
  n = 6

  call gegenbauer_cc1 ( n, lambda, f, value )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Value = ', value
  exact = gamma ( lambda + 0.5D+00 ) * sqrt ( r8_pi ) * besselj ( lambda, a ) &
    / ( 0.5D+00 * a ) ** lambda
  write ( *, '(a,g14.6)' ) '  Exact = ', exact

  return
end
subroutine gegenbauer_cc2_test ( )

!*****************************************************************************80
!
!! GEGENBAUER_CC2_TEST tests GEGENBAUER_CC2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) besselj
  real ( kind = rk ) exact
  real ( kind = rk ), external :: f
  real ( kind = rk ) lambda
  integer n
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) value

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'GEGENBAUER_CC2_TEST:'
  write ( *, '(a)' ) '  GEGENBAUER_CC2 estimates the Gegenbauer integral of'
  write ( *, '(a)' ) '  a function f(x) using a Clenshaw-Curtis type approach'
  write ( *, '(a)' ) '  based on the zeros of Tn(x).'

  lambda = 0.75D+00
  a = 2.0D+00
  n = 6

  call gegenbauer_cc2 ( n, lambda, f, value )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Value = ', value
  exact = gamma ( lambda + 0.5D+00 ) * sqrt ( r8_pi ) * besselj ( lambda, a ) &
    / ( 0.5D+00 * a ) ** lambda
  write ( *, '(a,g14.6)' ) '  Exact = ', exact

  return
end
subroutine i4_uniform_ab_test ( )

!*****************************************************************************80
!
!! I4_UNIFORM_AB_TEST tests I4_UNIFORM_AB.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 October 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: a = -100
  integer, parameter :: b = 200
  integer i
  integer i4_uniform_ab
  integer j
  integer seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'I4_UNIFORM_AB_TEST'
  write ( *, '(a)' ) '  I4_UNIFORM_AB computes pseudorandom values '
  write ( *, '(a)' ) '  in an interval [A,B].'

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The lower endpoint A = ', a
  write ( *, '(a,i12)' ) '  The upper endpoint B = ', b
  write ( *, '(a,i12)' ) '  The initial seed is ', seed
  write ( *, '(a)' ) ' '

  do i = 1, 20

    j = i4_uniform_ab ( a, b, seed )

    write ( *, '(2x,i8,2x,i8)' ) i, j

  end do

  return
end
subroutine r8_mop_test ( )

!*****************************************************************************80
!
!! R8_MOP_TEST tests R8_MOP.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 December 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i4
  integer i4_max
  integer i4_min
  integer i4_uniform_ab
  real ( kind = rk ) r8
  real ( kind = rk ) r8_mop
  integer seed
  integer test

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8_MOP_TEST'
  write ( *, '(a)' ) '  R8_MOP evaluates (-1.0)^I4 as an R8.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '    I4  R8_MOP(I4)'
  write ( *, '(a)' ) ''

  i4_min = -100
  i4_max = +100
  seed = 123456789

  do test = 1, 10
    i4 = i4_uniform_ab ( i4_min, i4_max, seed )
    r8 = r8_mop ( i4 )
    write ( *, '(2x,i4,2x,f4.1)' ) i4, r8
  end do

  return
end
subroutine r8vec_print_test ( )

!*****************************************************************************80
!
!! R8VEC_PRINT_TEST tests R8VEC_PRINT.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 August 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 4

  real ( kind = rk ), dimension ( n ) :: a = (/ &
    123.456D+00, 0.000005D+00, -1.0D+06, 3.14159265D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8VEC_PRINT_TEST'
  write ( *, '(a)' ) '  R8VEC_PRINT prints an R8VEC.'

  call r8vec_print ( n, a, '  The R8VEC:' )

  return
end
subroutine r8vec2_print_test ( )

!*****************************************************************************80
!
!! R8VEC2_PRINT_TEST tests R8VEC2_PRINT.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 January 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) :: a(n) = (/ &
    1.0D+00, 2.0D+00, 3.0D+00, 4.0D+00, 5.0D+00 /)
  real ( kind = rk ) b(n)
  real ( kind = rk ) c(n)
  integer i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8VEC2_PRINT_TEST'
  write ( *, '(a)' ) '  R8VEC2_PRINT prints a pair of R8VEC''s.'

  do i = 1, n
    b(i) = a(i) ** 2
    c(i) = sqrt ( a(i) )
  end do

  call r8vec_print ( n, a, '  Squares and square roots:' )

  return
end
