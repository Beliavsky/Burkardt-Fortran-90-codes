program main

!*****************************************************************************80
!
!! circle_rule_test() tests circle_rule().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer nt

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'circle_rule_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test circle_rule().'

  nt = 8
  call test01 ( nt )

  nt = 32
  call test01 ( nt )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'circle_rule_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine test01 ( nt )

!*****************************************************************************80
!
!! test01() tests circle_rule().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nt

  integer e(2)
  integer e1
  integer e2
  real ( kind = rk ) exact
  integer i
  real ( kind = rk ) q
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) t(nt)
  real ( kind = rk ) w(nt)
  real ( kind = rk ) x
  real ( kind = rk ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test01():'
  write ( *, '(a)' ) '  circle_rule() computes a rule Q(f) for the unit circle'
  write ( *, '(a)' ) '  using NT equally spaced angles.'
  write ( *, '(a)' ) '  Estimate integrals I(f) where f = x^e(1) * y^e(2)'
  write ( *, '(a,i4,a)' ) '  using ', nt, ' points.'
!
!  Compute the quadrature rule.
!
  call circle_rule ( nt, w, t )
!
!  Apply it to integrands.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  E(1)  E(2)    I(f)            Q(f)' 
  write ( *, '(a)' ) ' '
!
!  Specify a monomial.
!
  do e1 = 0, 6, 2

    e(1) = e1

    do e2 = e1, 6, 2

      e(2) = e2

      q = 0.0D+00
      do i = 1, nt
        x = cos ( t(i) )
        y = sin ( t(i) )
        q = q + w(i) * x ** e(1) * y ** e(2)
      end do

      q = 2.0D+00 * r8_pi * q

      call circle01_monomial_integral ( e, exact )

      write ( *, '(3x,i2,3x,i2,2x,g14.6,2x,g14.6)' ) e(1), e(2), exact, q

    end do

  end do

  return
end

