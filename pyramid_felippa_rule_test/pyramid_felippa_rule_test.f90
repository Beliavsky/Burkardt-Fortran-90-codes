program main

!*****************************************************************************80
!
!! pyramid_felippa_rule_test() tests pyramid_felippa_rule().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 August 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer degree_max

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'pyramid_felippa_rule_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test PYRAMID_FELIPPA_RULE().'

  degree_max = 4
  call pyramid_unit_monomial_test ( degree_max )

  degree_max = 5
  call pyramid_unit_quad_test ( degree_max )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PYRAMID_FELIPPA_RULE_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine pyramid_unit_monomial_test ( degree_max )

!*****************************************************************************80
!
!! PYRAMID_UNIT_MONOMIAL_TEST tests PYRAMID_UNIT_MONOMIAL.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DEGREE_MAX, the maximum total degree of the
!    monomials to check.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer alpha
  integer beta
  integer degree_max
  integer expon(3)
  integer gamma
  real ( kind = rk ) pyramid_unit_volume
  real ( kind = rk ) value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PYRAMID_UNIT_MONOMIAL_TEST'
  write ( *, '(a)' ) '  For the unit pyramid,'
  write ( *, '(a)' ) '  PYRAMID_UNIT_MONOMIAL returns the exact value of the'
  write ( *, '(a)' ) '  integral of X^ALPHA Y^BETA Z^GAMMA'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Volume = ', pyramid_unit_volume ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     ALPHA      BETA     GAMMA      INTEGRAL'
  write ( *, '(a)' ) ' '

  do alpha = 0, degree_max
    expon(1) = alpha
    do beta = 0, degree_max - alpha
      expon(2) = beta
      do gamma = 0, degree_max - alpha - beta
        expon(3) = gamma
        call pyramid_unit_monomial ( expon, value )
        write ( *, '(2x,i8,2x,i8,2x,i8,2x,g14.6)' ) expon(1:3), value
      end do
    end do
  end do

  return
end
subroutine pyramid_unit_quad_test ( degree_max )

!*****************************************************************************80
!
!! PYRAMID_UNIT_QUAD_TEST tests the rules for the unit pyramid.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DEGREE_MAX, the maximum total degree of the
!    monomials to check.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  integer degree_max
  integer expon(dim_num)
  integer h
  logical more
  integer order
  real ( kind = rk ) quad
  integer t
  real ( kind = rk ) pyramid_unit_volume
  real ( kind = rk ), allocatable :: v(:)
  real ( kind = rk ), allocatable :: w(:)
  real ( kind = rk ), allocatable :: xyz(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PYRAMID_UNIT_QUAD_TEST'
  write ( *, '(a)' ) '  For the unit pyramid,'
  write ( *, '(a)' ) '  we approximate monomial integrals with:'
  write ( *, '(a)' ) '  PYRAMID_UNIT_O01,'
  write ( *, '(a)' ) '  PYRAMID_UNIT_O05,'
  write ( *, '(a)' ) '  PYRAMID_UNIT_O06,'
  write ( *, '(a)' ) '  PYRAMID_UNIT_O08,'
  write ( *, '(a)' ) '  PYRAMID_UNIT_O08b,'
  write ( *, '(a)' ) '  PYRAMID_UNIT_O09,'
  write ( *, '(a)' ) '  PYRAMID_UNIT_O13,'
  write ( *, '(a)' ) '  PYRAMID_UNIT_O18,'
  write ( *, '(a)' ) '  PYRAMID_UNIT_O27,'
  write ( *, '(a)' ) '  PYRAMID_UNIT_O48.'

  more = .false.

  do

    call subcomp_next ( degree_max, dim_num, expon, more, h, t )

    if ( mod ( expon(1), 2 ) == 1 .or. &
         mod ( expon(2), 2 ) == 1 ) then
      cycle
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a,2x,i2,2x,i2,2x,i2)' ) &
      '  Monomial exponents: ', expon(1:dim_num)
    write ( *, '(a)' ) ' '

    order = 1
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( xyz(1:dim_num,1:order) )
    call pyramid_unit_o01 ( w, xyz )
    call monomial_value ( dim_num, order, expon, xyz, v )
    quad = pyramid_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,g14.6)' ) order, quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( xyz )

    order = 5
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( xyz(1:dim_num,1:order) )
    call pyramid_unit_o05 ( w, xyz )
    call monomial_value ( dim_num, order, expon, xyz, v )
    quad = pyramid_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,g14.6)' ) order, quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( xyz )

    order = 6
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( xyz(1:dim_num,1:order) )
    call pyramid_unit_o06 ( w, xyz )
    call monomial_value ( dim_num, order, expon, xyz, v )
    quad = pyramid_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,g14.6)' ) order, quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( xyz )

    order = 8
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( xyz(1:dim_num,1:order) )
    call pyramid_unit_o08 ( w, xyz )
    call monomial_value ( dim_num, order, expon, xyz, v )
    quad = pyramid_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,g14.6)' ) order, quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( xyz )

    order = 8
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( xyz(1:dim_num,1:order) )
    call pyramid_unit_o08b ( w, xyz )
    call monomial_value ( dim_num, order, expon, xyz, v )
    quad = pyramid_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,g14.6)' ) order, quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( xyz )

    order = 9
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( xyz(1:dim_num,1:order) )
    call pyramid_unit_o09 ( w, xyz )
    call monomial_value ( dim_num, order, expon, xyz, v )
    quad = pyramid_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,g14.6)' ) order, quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( xyz )

    order = 13
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( xyz(1:dim_num,1:order) )
    call pyramid_unit_o13 ( w, xyz )
    call monomial_value ( dim_num, order, expon, xyz, v )
    quad = pyramid_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,g14.6)' ) order, quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( xyz )

    order = 18
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( xyz(1:dim_num,1:order) )
    call pyramid_unit_o18 ( w, xyz )
    call monomial_value ( dim_num, order, expon, xyz, v )
    quad = pyramid_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,g14.6)' ) order, quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( xyz )

    order = 27
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( xyz(1:dim_num,1:order) )
    call pyramid_unit_o27 ( w, xyz )
    call monomial_value ( dim_num, order, expon, xyz, v )
    quad = pyramid_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,g14.6)' ) order, quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( xyz )

    order = 48
    allocate ( v(1:order) )
    allocate ( w(1:order) )
    allocate ( xyz(1:dim_num,1:order) )
    call pyramid_unit_o48 ( w, xyz )
    call monomial_value ( dim_num, order, expon, xyz, v )
    quad = pyramid_unit_volume ( ) * dot_product ( w(1:order), v(1:order) )
    write ( *, '(2x,i6,2x,g14.6)' ) order, quad
    deallocate ( v )
    deallocate ( w )
    deallocate ( xyz )

    write ( *, '(a)' ) ' '
    call pyramid_unit_monomial ( expon, quad )
    write ( *, '(2x,a,2x,g14.6)' ) ' Exact', quad

    if ( .not. more ) then
      exit
    end if

  end do

  return
end

