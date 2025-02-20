program main

!*****************************************************************************80
!
!! triangle_twb_rule_test() tests triangle_twb_rule().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 April 2019
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer degree_max

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'triangle_twb_rule_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test triangle_twb_rule().'

  degree_max = 5
  call triangle_unit_quad_test ( degree_max )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'triangle_twb_rule_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine triangle_unit_quad_test ( degree_max )

!*****************************************************************************80
!
!! triangle_unit_quad_test() tests rules for the unit triangle.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 April 2019
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, int DEGREE_MAX, the maximum total degree of the monomials to check.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer degree
  integer degree_max
  integer ex
  integer ey
  integer n
  real ( kind = rk ) q
  integer strength
  real ( kind = rk ), allocatable :: v(:)
  real ( kind = rk ), allocatable :: w(:)
  real ( kind = rk ), allocatable :: x(:)
  real ( kind = rk ), allocatable :: y(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'triangle_unit_quad_test():'
  write ( *, '(a)' ) '  Approximate monomial integrals in triangle with TWB rules.'

  degree = 0
  ex = 0
  ey = degree

  do

    write ( *, '(a)' ) ''
    write ( *, '(a,i2,a,i2)' ) '  Monomial: x^', ex, ' y^', ey

    do strength = 1, 25

      call twb_rule_n ( strength, n )

      if ( n < 1 ) then
        cycle
      end if

      allocate ( w(1:n) )
      allocate ( x(1:n) )
      allocate ( y(1:n) )
      allocate ( v(1:n) )

      call twb_rule_w ( strength, w )
      call twb_rule_x ( strength, x )
      call twb_rule_y ( strength, y )
      call monomial_value_2d ( n, ex, ey, x, y, v )
      q = dot_product ( w, v )
      write ( *, '(2x,i6,2x,i6,2x,g14.6)' ) strength, n, q

      deallocate ( v )
      deallocate ( w )
      deallocate ( x )
      deallocate ( y )

    end do

    call triangle_unit_monomial ( ex, ey, q )
    write ( *, '(a,g14.6)' ) '   Exact          ', q

    if ( ex < degree ) then
      ex = ex + 1
      ey = ey - 1
    else if ( degree < degree_max ) then
      degree = degree + 1
      ex = 0
      ey = degree
    else
      exit
    end if

  end do

  return
end
