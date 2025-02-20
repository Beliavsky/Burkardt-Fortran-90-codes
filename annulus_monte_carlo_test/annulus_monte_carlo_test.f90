program main

!*****************************************************************************80
!
!! annulus_monte_carlo_test() tests annulus_monte_carlo().
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
  write ( *, '(a)' ) 'annulus_monte_carlo_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test annulus_monte_carlo().'

  call annulus_area_test ( )

  center = (/ 0.0_rk, 0.0_rk /)
  r1 = 0.0_rk
  r2 = 1.0_rk
  call annulus_sample_test ( center, r1, r2 )

  center = (/ 0.0_rk, 0.0_rk /)
  r1 = 0.5_rk
  r2 = 1.0_rk
  call annulus_sample_test ( center, r1, r2 )

  center = (/ 1.0_rk, 0.0_rk /)
  r1 = 0.0_rk
  r2 = 1.0_rk
  call annulus_sample_test ( center, r1, r2 )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'annulus_monte_carlo_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop ( 0 )
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
!    05 July 2018
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
    call random_number ( harvest = dat )
    center(1) = 10.0_rk * dat(1) - 5.0_rk
    center(2) = 10.0_rk * dat(2) - 5.0_rk
    r1 = dat(3)
    r2 = r1 + dat(4)
    area = annulus_area ( r1, r2 )
    write ( *, '(2x,a1,f9.6,a1,f9.6,a1,2x,f9.6,2x,f9.6,2x,f9.6)' ) &
      '(', center(1), ',', center(2), ',', r1, r2, area
  end do

  return
end
subroutine annulus_sample_test ( center, r1, r2 )

!*****************************************************************************80
!
!! annulus_sample_test() tests annulus_sample().
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
!    05 July 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) annulus_area
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
  real ( kind = rk ) r1
  real ( kind = rk ) r2
  real ( kind = rk ) result
  real ( kind = rk ), allocatable :: value(:)
  real ( kind = rk ), allocatable :: x(:,:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'annulus_sample_test():'
  write ( *, '(a)' ) '  annulus_sample() samples an annulus uniformly.'
  write ( *, '(a)' ) '  Use it to estimate integrals in the annulus'
  write ( *, '(a,f8.4,a,f8.4,a,f8.4,a,f8.4)' ) &
    '  centered at (', center(1), ',', center(2), &
    ') with R1 = ', r1, ' R2 = ', r2

  write ( *, '(a)' ) ''
  write ( *, '(a)', advance = 'no' ) &
    '         N        1              X^2             Y^2        '
  write ( *, '(a)' ) '     X^4             X^2Y^2           Y^4             X^6'
  write ( *, '(a)' ) ''

  n = 1

  do while ( n <= 65536 )

    allocate ( x(1:2,1:n) )
    allocate ( value(1:n) )

    call annulus_sample ( center, r1, r2, n, x )

    write ( *, '(2x,i8)', advance = 'no' ) n
    do j = 1, 7
      e(1:2) = e_test(1:2,j)
      call monomial_value ( 2, n, e, x, value )
      result = annulus_area ( r1, r2 ) &
        * sum ( value(1:n) ) / real ( n, kind = rk )
      write ( *, '(2x,g14.6)', advance = 'no' ) result
    end do
    write ( *, '(a)' ) ''

    deallocate ( value )
    deallocate ( x )

    n = 2 * n

  end do

  if ( &
    center(1) == 0.0_rk .and. &
    center(2) == 0.0_rk .and. &
    r1 == 0.0_rk .and. &
    r2 == 1.0_rk ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)', advance = 'no' ) '     Exact'
    do j = 1, 7
      e(1:2) = e_test(1:2,j)
      call disk01_monomial_integral ( e, result )
      write ( *, '(2x,g14.6)', advance = 'no' ) result
    end do
    write ( *, '(a)' ) ''

  end if

  return
end
 
