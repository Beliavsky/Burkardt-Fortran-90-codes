program main

!*****************************************************************************80
!
!! annulus_rule_test() tests annulus_rule().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) center(2)
  real ( kind = rk ) r1
  real ( kind = rk ) r2

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'annulus_rule_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test annulus_rule().'

  call annulus_area_test ( )
  call annulus_rule_compute_test ( )

  center = (/ 0.0D+00, 0.0D+00 /)
  r1 = 0.0D+00
  r2 = 1.0D+00
  call annulus_rule_monomial_test ( center, r1, r2 )

  center = (/ 0.0D+00, 0.0D+00 /)
  r1 = 0.5D+00
  r2 = 1.0D+00
  call annulus_rule_monomial_test ( center, r1, r2 )

  center = (/ 1.0D+00, 0.0D+00 /)
  r1 = 0.0D+00
  r2 = 1.0D+00
  call annulus_rule_monomial_test ( center, r1, r2 )

  center = (/ 0.0D+00, 0.0D+00 /)
  r1 = 0.0D+00
  r2 = 1.0D+00
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'annulus_rule_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine annulus_area_test ( )

!*****************************************************************************80
!
!! annulus_area_test() test annulus_area().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) annulus_area
  real ( kind = rk ) area
  real ( kind = rk ) center(2)
  real ( kind = rk ) dat(4)
  integer i
  real ( kind = rk ) r1
  real ( kind = rk ) r2

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'annulus_area_test():'
  write ( *, '(a)' ) '  annulus_area() computes the area of an annulus with'
  write ( *, '(a)' ) '  center = (CX,CY), inner radius R1 and outer radius R2.'
  
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  (   CX        CY     )    R1         R2         Area'
  write ( *, '(a)' ) ''

  do i = 1, 10
    call random_number ( harvest = dat(1:4) )
    center(1) = 10.0 * dat(1) - 5.0
    center(2) = 10.0 * dat(2) - 5.0
    r1 = dat(3)
    r2 = r1 + dat(4)
    area = annulus_area ( center, r1, r2 )
    write ( *, '(2x,a,f9.6,a,f9.6,a,2x,f9.6,2x,f9.6,2x,f9.6)' ) &
    '(', center(1), ',', center(2), ')', r1, r2, area
  end do

  return
end
subroutine annulus_rule_compute_test ( )

!*****************************************************************************80
!
!! annulus_rule_compute_test() tests annulus_rule_compute().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) center(2)
  integer n
  integer nr
  integer nt
  real ( kind = rk ) r1
  real ( kind = rk ) r2
  real ( kind = rk ), allocatable :: w(:)
  real ( kind = rk ), allocatable :: x(:)
  real ( kind = rk ), allocatable :: y(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'annulus_rule_compute_test():'
  write ( *, '(a)' ) '  Test annulus_rule_compute().'

  center = (/ 0.0D+00, 0.0D+00 /)
  r1 = 0.5D+00
  r2 = 1.0D+00
  nr = 3
  nt = 12
  n = nt * nr

  allocate ( w(1:n) )
  allocate ( x(1:n) )
  allocate ( y(1:n) )

  call annulus_rule_compute ( center, r1, r2, nr, nt, w, x, y )

  call r8vec3_print ( n, w, x, y, '  W, X, Y for annulus quadrature:' )

  deallocate ( w )
  deallocate ( x )
  deallocate ( y )

  return
end
subroutine annulus_rule_monomial_test ( center, r1, r2 )

!*****************************************************************************80
!
!! annulus_rule_monomial_test() estimates monomial integrals using quadrature
!
!  Discussion:
!
!    If CENTER=(0,0) and R1 = 0 and R2 = 1, then we can compare exact values.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 August 2021
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) CENTER(2), the coordinates of the center.
!
!    Input, real ( kind = rk ) R1, R2, the inner and outer radii of the annulus.
!    0.0 <= R1 <= R2.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) center(2)
  integer e(2)
  integer :: e_test(2,7) = reshape ( (/ &
    0, 0, &
    2, 0, &
    0, 2, &
    4, 0, &
    2, 2, &
    0, 4, &
    6, 0 /), (/ 2, 7 /) )
  integer j
  integer n
  integer nr
  integer nt
  real ( kind = rk ) r1
  real ( kind = rk ) r2
  real ( kind = rk ) result
  real ( kind = rk ), allocatable :: value(:)
  real ( kind = rk ), allocatable :: w(:)
  real ( kind = rk ), allocatable :: x(:)
  real ( kind = rk ), allocatable :: xy(:,:)
  real ( kind = rk ), allocatable :: y(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'annulus_rule_monomial_test():'
  write ( *, '(a)' ) '  annulus_rule_compute() supplies a quadrature rule for'
  write ( *, '(a,f8.4,a,f8.4,a,f8.4,a,f8.4)' ) &
    '  the annulus centered at (', center(1), ',', center(2), &
    ') with R1 = ', r1, ' R2 = ', r2

  write ( *, '(a)' ) ''
  write ( *, '(a)', advance = 'no' ) &
    '    NR    NT       1              X^2             Y^2        '
  write ( *, '(a)' ) '     X^4             X^2Y^2          Y^4             X^6'
  write ( *, '(a)' ) ''

  nr = 4

  do while ( nr <= 64 )

    nt = 4 * nr

    n = nr * nt

    allocate ( value(1:n) )
    allocate ( w(1:n) )
    allocate ( x(1:n) )
    allocate ( xy(1:2,1:n) )
    allocate ( y(1:n) )

    call annulus_rule_compute ( center, r1, r2, nr, nt, w, x, y )

    xy(1,1:n) = x(1:n)
    xy(2,1:n) = y(1:n)

    write ( *, '(2x,i4,2x,i4)', advance = 'no' ) nr, nt
    do j = 1, 7
      e(1:2) = e_test(1:2,j)
      call monomial_value ( 2, n, e, xy, value )
      result = dot_product ( w(1:n), value(1:n) )
      write ( *, '(2x,g14.6)', advance = 'no' ) result
    end do
    write ( *, '(a)' ) ''

    deallocate ( value )
    deallocate ( w )
    deallocate ( x )
    deallocate ( xy )
    deallocate ( y )

    nr = 2 * nr

  end do

  if ( &
    center(1) == 0.0D+00 .and. &
    center(2) == 0.0D+00 .and. &
    r1 == 0.0D+00 .and. &
    r2 == 1.0D+00 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)', advance = 'no' ) '     Exact  '
    do j = 1, 7
      e(1:2) = e_test(1:2,j)
      call disk01_monomial_integral ( e, result )
      write ( *, '(2x,g14.6)', advance = 'no' ) result
    end do
    write ( *, '(a)' ) ''

  end if

  return
end

