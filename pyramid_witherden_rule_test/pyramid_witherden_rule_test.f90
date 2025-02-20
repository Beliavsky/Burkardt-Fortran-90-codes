program main

!*****************************************************************************80
!
!! pyramid_witherden_rule_test() tests pyramid_witherden_rule().
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license.
!
!  Modified:
!
!    13 July 2023
!
!  Author:
!
!    John Burkardt.
!
  implicit none

  integer p
  integer p_hi
  integer p_lo

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'pyramid_witherden_rule_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test pyramid_witherden_rule().'

  p = 5
  call pyramid_witherden_rule_test01 ( p )

  p = 5
  call pyramid_witherden_rule_test02 ( p )

  p_lo = 0
  p_hi = 10
  call pyramid_witherden_rule_test03 ( p_lo, p_hi )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'pyramid_witherden_rule_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  return
end
subroutine pyramid_witherden_rule_test01 ( p )

!*****************************************************************************80
!
!! pyramid_witherden_rule_test01() computes a quadrature rule of precision P.
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license.
!
!  Modified:
!
!    26 May 2023
!
!  Author:
!
!    John Burkardt.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer n
  integer p
  real ( kind = rk ), allocatable :: w(:)
  real ( kind = rk ) w_sum
  real ( kind = rk ), allocatable :: x(:)
  real ( kind = rk ), allocatable :: y(:)
  real ( kind = rk ), allocatable :: z(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'pyramid_witherden_rule_test01():'
  write ( *, '(a)' ) '  Quadrature rule for the unit pyramid,'
  write ( *, '(a,i2)' ) '  Precision p = ', p
!
!  Retrieve the rule.
!
  call rule_order ( p, n )
  write ( *, '(a)' ) ''
  write ( *, '(a,i3)' ) '  Number of nodes N = ', n

  allocate ( x(1:n) )
  allocate ( y(1:n) )
  allocate ( z(1:n) )
  allocate ( w(1:n) )

  call pyramid_witherden_rule ( p, n, x, y, z, w )
!
!  Print the rule.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '     I      W          X           Y           Z'
  write ( *, '(a)' ) ''
  do i = 1, n
    write ( *, '(2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      i, w(i), x(i), y(i), z(i)
  end do
!
!  Verify weights sum to 1.
!
  w_sum = sum ( w )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Weight Sum ', w_sum
!
!  Free memory.
!
  deallocate ( w )
  deallocate ( x )
  deallocate ( y )
  deallocate ( z )

  return
end
subroutine pyramid_witherden_rule_test02 ( p )

!*****************************************************************************80
!
!! pyramid_witherden_rule_test02() tests a quadrature rule of precision P.
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license.
!
!  Modified:
!
!    16 May 2023
!
!  Author:
!
!    John Burkardt.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer degree
  integer dim_num
  integer expon(3)
  real ( kind = rk ) exact
  integer h
  real ( kind = rk ) max_error
  logical more
  integer n
  integer p
  real ( kind = rk ) pyramid_unit_volume
  real ( kind = rk ) quad
  real ( kind = rk ) quad_error
  integer t
  real ( kind = rk ), allocatable :: v(:)
  real ( kind = rk ), allocatable :: w(:)
  real ( kind = rk ), allocatable :: x(:)
  real ( kind = rk ), allocatable :: xyz(:,:)
  real ( kind = rk ), allocatable :: y(:)
  real ( kind = rk ), allocatable :: z(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'pyramid_witherden_rule_test02():'
  write ( *, '(a)' ) '  Test the precision of a quadrature rule for the unit pyramid,'

  dim_num = 3
!
!  Retrieve the rule.
!
  call rule_order ( p, n )
  write ( *, '(a)' ) ''
  write ( *, '(a,i3)' ) '  Number of nodes N = ', n

  allocate ( x(1:n) )
  allocate ( y(1:n) )
  allocate ( z(1:n) )
  allocate ( w(1:n) )

  call pyramid_witherden_rule ( p, n, x, y, z, w )
!
!  Pack the x, y, z vectors as columns of an array.
!
  allocate ( xyz(1:n,1:3) )
  xyz = reshape ( (/ x, y, z /), (/ n, 3 /) )

  write ( *, '(a)' ) ''
  write ( *, '(a,i2)' ) '  Stated precision of rule    = ', p 
  write ( *, '(a,i3)' ) '  Number of quadrature points = ', n
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Degree  Maximum error'
  write ( *, '(a)' ) ''

  allocate ( v(1:n) )

  do degree = 0, p + 2

    expon = (/ 0, 0, 0 /)
    more = .false.
    h = 0
    t = 0
    max_error = 0.0D+00

    do while ( .true. )

      call comp_next ( degree, dim_num, expon, more, h, t )

      call monomial_value ( dim_num, n, expon, xyz, v )

      quad = pyramid_unit_volume ( ) * dot_product ( w, v )
 
      call pyramid_unit_monomial_integral ( expon, exact )

      quad_error = abs ( quad - exact )
      max_error = max ( max_error, quad_error )

      if ( .not. more ) then
        exit
      end if

    end do

    write ( *, '(2x,i2,2x,g24.16)' ) degree, max_error

  end do
!
!  Free memory.
!
  deallocate ( v )
  deallocate ( w )
  deallocate ( x )
  deallocate ( xyz )
  deallocate ( y )
  deallocate ( z )

  return
end
subroutine pyramid_witherden_rule_test03 ( p_lo, p_hi )

!*****************************************************************************80
!
!! pyramid_witherden_rule_test03() tests absolute and relative precision.
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license.
!
!  Modified:
!
!    13 July 2023
!
!  Author:
!
!    John Burkardt.
!
!  Input:
!
!    integer p_lo, p_hi: the lowest and highest rules to check.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer degree
  integer dim_num
  integer expon(3)
  real ( kind = rk ) exact
  integer h
  real ( kind = rk ) max_abs
  real ( kind = rk ) max_rel
  logical more
  integer n
  integer p
  integer p_hi
  integer p_lo
  real ( kind = rk ) pyramid_unit_volume
  real ( kind = rk ) quad
  real ( kind = rk ) quad_error
  integer t
  real ( kind = rk ), allocatable :: v(:)
  real ( kind = rk ), allocatable :: w(:)
  real ( kind = rk ), allocatable :: x(:)
  real ( kind = rk ), allocatable :: xyz(:,:)
  real ( kind = rk ), allocatable :: y(:)
  real ( kind = rk ), allocatable :: z(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'pyramid_witherden_rule_test03():'
  write ( *, '(a)' ) '  Test the precision of quadrature rules for the unit pyramid,'
  write ( *, '(a,i2,a,i2)' ) '  Check rules of precision p =', p_lo, ' through ', p_hi
  write ( *, '(a)' ) '  for error in approximating integrals of monomials.'

  dim_num = 3

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '              maximum                   maximum'
  write ( *, '(a)' ) '   p          absolute                  relative'
  write ( *, '(a)' ) '              error                     error'
  write ( *, '(a)' ) ''

  do p = p_lo, p_hi

    call rule_order ( p, n )

    allocate ( x(1:n) )
    allocate ( y(1:n) )
    allocate ( z(1:n) )
    allocate ( w(1:n) )

    call pyramid_witherden_rule ( p, n, x, y, z, w )
!
!  Pack the x, y, z vectors as columns of an array.
!
    allocate ( xyz(1:n,1:3) )
    xyz = reshape ( (/ x, y, z /), (/ n, 3 /) )

    allocate ( v(1:n) )

    max_abs = 0.0D+00
    max_rel = 0.0D+00

    do degree = 0, p

      expon = (/ 0, 0, 0 /)
      more = .false.
      h = 0
      t = 0
 
      do while ( .true. )

        call comp_next ( degree, dim_num, expon, more, h, t )

        call monomial_value ( dim_num, n, expon, xyz, v )

        quad = pyramid_unit_volume ( ) * dot_product ( w, v )
 
        call pyramid_unit_monomial_integral ( expon, exact )

        quad_error = abs ( quad - exact )
        max_abs = max ( max_abs, quad_error )

        if ( exact /= 0.0D+00 ) then
          max_rel = max ( max_rel, quad_error / abs ( exact ) )
        end if 

        if ( .not. more ) then
          exit
        end if

      end do

    end do

    write ( *, '(2x,i2,2x,g24.16,2x,g24.16)' ) p, max_abs, max_rel

    deallocate ( v )
    deallocate ( w )
    deallocate ( x )
    deallocate ( xyz )
    deallocate ( y )
    deallocate ( z )

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
