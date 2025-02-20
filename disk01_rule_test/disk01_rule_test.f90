program main

!*****************************************************************************80
!
!! DISK01_RULE_TEST() tests DISK01_RULE().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 April 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'DISK01_RULE_TEST:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the DISK01_RULE library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'DISK01_RULE_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests DISK01_RULE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 March 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nr = 4
  integer, parameter :: nt = 8

  real ( kind = rk ) area
  integer e(2)
  integer e1
  integer e2
  real ( kind = rk ) exact
  integer i
  integer j
  real ( kind = rk ) q
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) r(nr)
  real ( kind = rk ) s
  real ( kind = rk ) t(nt)
  real ( kind = rk ) w(nr)
  real ( kind = rk ) x
  real ( kind = rk ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  DISK01_RULE can compute a rule Q(f) for the unit disk'
  write ( *, '(a)' ) '  using NT equally spaced angles and NR radial distances.'
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  NT = ', nt
  write ( *, '(a,i4)' ) '  NR = ', nr
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Estimate integrals I(f) where f = x^e(1) * y^e(2).'
!
!  Compute the quadrature rule.
!
  call disk01_rule ( nr, nt, w, r, t )
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

      s = 0.0D+00
      do j = 1, nt
        do i = 1, nr
          x = r(i) * cos ( t(j) )
          y = r(i) * sin ( t(j) )
          s = s + w(i) * x ** e(1) * y ** e(2)
        end do
      end do

      area = r8_pi

      q = r8_pi * s

      call disk01_monomial_integral ( e, exact )

      write ( *, '(3x,i2,3x,i2,2x,g14.6,2x,g14.6)' ) e(1), e(2), exact, q

    end do

  end do

  return
end

