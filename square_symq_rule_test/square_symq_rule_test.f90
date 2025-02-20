program main

!*****************************************************************************80
!
!! square_symq_rule_test() tests square_symq_rule().
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license.
!
!  Modified:
!
!    16 July 2023
!
!  Author:
!
!    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
!    This version by John Burkardt.
!
!  Reference:
!
!    Hong Xiao, Zydrunas Gimbutas,
!    A numerical algorithm for the construction of efficient quadrature
!    rules in two and higher dimensions,
!    Computers and Mathematics with Applications,
!    Volume 59, 2010, pages 663-676.
!
  implicit none

  integer p
  integer p_hi
  integer p_lo

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'square_symq_rule_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test square_symq_rule().'

  p = 5
  call test01 ( p )

  p = 5
  call test02 ( p )

  p_lo = 0
  p_hi = 20
  call test03 ( p_lo, p_hi )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'square_symq_rule_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine test01 ( p )

!*****************************************************************************80
!
!! test01() prints a quadrature rule of given precision.
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license.
!
!  Modified:
!
!    16 July 2023
!
!  Author:
!
!    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
!    This version by John Burkardt.
!
!  Reference:
!
!    Hong Xiao, Zydrunas Gimbutas,
!    A numerical algorithm for the construction of efficient quadrature
!    rules in two and higher dimensions,
!    Computers and Mathematics with Applications,
!    Volume 59, 2010, pages 663-676.
!
!  Input:
!
!    integer p: the precision of the rule.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) area
  integer j
  integer n
  integer p
  real ( kind = rk ), allocatable :: w(:)
  real ( kind = rk ) w_sum
  real ( kind = rk ), allocatable :: x(:)
  real ( kind = rk ), allocatable :: y(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'test01():'
  write ( *, '(a)' ) '  Symmetric quadrature rule for a square.'
  write ( *, '(a,i4)' )'  Precision = ', p

  area = 4.0D+00
!
!  Retrieve and print a symmetric quadrature rule.
!
  call rule_order ( p, n )

  allocate ( w(1:n) )
  allocate ( x(1:n) )
  allocate ( y(1:n) )
  call square_symq ( p, n, x, y, w )

  write ( *, '(a)' ) ''
  write ( *, '(a,i6)' ) '  Number of nodes N = ', n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '     J  W       X       Y'
  write ( *, '(a)' ) ''
  do j = 1, n
    write ( *, '(2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      j, w(j), x(j), y(j)
  end do

  w_sum = sum ( w(1:n) )

  write ( *, '(a,2x,g14.6)' ) '  Sum of weight:', w_sum

  deallocate ( w )
  deallocate ( x )
  deallocate ( y )

  return
end
subroutine test02 ( p )

!*****************************************************************************80
!
!! test02() tests a rule of precision P.
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license.
!
!  Modified:
!
!    16 July 2023
!
!  Author:
!
!    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
!    This version by John Burkardt.
!
!  Reference:
!
!    Hong Xiao, Zydrunas Gimbutas,
!    A numerical algorithm for the construction of efficient quadrature
!    rules in two and higher dimensions,
!    Computers and Mathematics with Applications,
!    Volume 59, 2010, pages 663-676.
!
!  Input:
!
!    integer p: the precision of the rule.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer degree
  integer dim_num
  integer expon(2)
  real ( kind = rk ) exact
  integer h
  real ( kind = rk ) max_error
  logical more
  integer n
  integer p
  real ( kind = rk ) quadrilateral_unit_area
  real ( kind = rk ) quad
  real ( kind = rk ) quad_error
  integer t
  real ( kind = rk ), allocatable :: v(:)
  real ( kind = rk ), allocatable :: w(:)
  real ( kind = rk ), allocatable :: x(:)
  real ( kind = rk ), allocatable :: xy(:,:)
  real ( kind = rk ), allocatable :: y(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'test02():'
  write ( *, '(a)' ) '  Get a quadrature rule for the symmetric square.'
  write ( *, '(a)' ) '  Test its accuracy.'
  write ( *, '(a,i4)' ) '  Polynomial precision = ', p

  dim_num = 2
!
!  Retrieve a symmetric quadrature rule.
!
  call rule_order ( p, n )
  allocate ( x(1:n) )
  allocate ( y(1:n) )
  allocate ( w(1:n) )

  call square_symq ( p, n, x, y, w )
!
!  Pack the x, y vectors as columns of an array.
!
  allocate ( xy(1:n,1:2) )
  xy = reshape ( (/ x, y /), (/ n, 2 /) )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Degree  Maximum error'
  write ( *, '(a)' ) ''

  allocate ( v(1:n) )

  do degree = 0, p + 2

    expon = (/ 0, 0 /)
    more = .false.
    h = 0
    t = 0
    max_error = 0.0D+00

    do while ( .true. )

      call comp_next ( degree, dim_num, expon, more, h, t )

      call monomial_value ( dim_num, n, expon, xy, v )

      quad = quadrilateral_unit_area ( ) * dot_product ( w, v )

      call quadrilateral_unit_monomial_integral ( expon, exact )

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
  deallocate ( xy )
  deallocate ( y )

  return
end
subroutine test03 ( p_lo, p_hi )

!*****************************************************************************80
!
!! test03() tests absolute and relative precision.
!
!  Licensing:
!
!    This code is distributed under the GNU GPL license.
!
!  Modified:
!
!    16 July 2023
!
!  Author:
!
!    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
!    This version by John Burkardt.
!
!  Reference:
!
!    Hong Xiao, Zydrunas Gimbutas,
!    A numerical algorithm for the construction of efficient quadrature
!    rules in two and higher dimensions,
!    Computers and Mathematics with Applications,
!    Volume 59, 2010, pages 663-676.
!
!  Input:
!
!    integer p_lo, p_hi: the lowest and highest rules to check.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer degree
  integer dim_num
  integer expon(2)
  real ( kind = rk ) exact
  integer h
  real ( kind = rk ) max_abs
  real ( kind = rk ) max_rel
  logical more
  integer n
  integer p
  integer p_hi
  integer p_lo
  real ( kind = rk ) quadrilateral_unit_area
  real ( kind = rk ) quad
  real ( kind = rk ) quad_error
  integer t
  real ( kind = rk ), allocatable :: v(:)
  real ( kind = rk ), allocatable :: w(:)
  real ( kind = rk ), allocatable :: x(:)
  real ( kind = rk ), allocatable :: xy(:,:)
  real ( kind = rk ), allocatable :: y(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'test02():'
  write ( *, '(a)' ) '  Test the precision of quadrature rules for the unit quadrilateral,'
  write ( *, '(a,i2,a,i2)' ) '  Check rules of precision p =', p_lo, ' through ', p_hi
  write ( *, '(a)' ) '  for error in approximating integrals of monomials.'

  dim_num = 2

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '              maximum                   maximum'
  write ( *, '(a)' ) '   p          absolute                  relative'
  write ( *, '(a)' ) '              error                     error'
  write ( *, '(a)' ) ''

  do p = p_lo, p_hi

    call rule_order ( p, n )
    allocate ( x(1:n) )
    allocate ( y(1:n) )
    allocate ( w(1:n) )

    call square_symq ( p, n, x, y, w )
!
!  Pack the x, y vectors as columns of an array.
!
    allocate ( xy(1:n,1:2) )
    xy = reshape ( (/ x, y /), (/ n, 2 /) )

    allocate ( v(1:n) )

    max_abs = 0.0D+00
    max_rel = 0.0D+00

    do degree = 0, p

      expon = (/ 0, 0 /)
      more = .false.
      h = 0
      t = 0

      do while ( .true. )

        call comp_next ( degree, dim_num, expon, more, h, t )

        call monomial_value ( dim_num, n, expon, xy, v )

        quad = quadrilateral_unit_area ( ) * dot_product ( w, v )

        call quadrilateral_unit_monomial_integral ( expon, exact )

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
    deallocate ( xy )
    deallocate ( y )

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
!    18 May 2013
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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end

