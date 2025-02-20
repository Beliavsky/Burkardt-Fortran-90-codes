program main

!*****************************************************************************80
!
!! sphere_triangle_monte_carlo_test() tests sphere_triangle_monte_carlo().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 October 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'sphere_triangle_monte_carlo_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test sphere_triangle_monte_carlo().'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'sphere_triangle_monte_carlo_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 uses SPHERE_TRIANGLE_SAMPLE_01 with an increasing number of points.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 3

  real ( kind = rk ) area
  integer e(m)
  integer :: e_test(m,7) = reshape ( (/ &
    0, 0, 0, &
    2, 0, 0, &
    0, 2, 0, &
    0, 0, 2, &
    4, 0, 0, &
    2, 2, 0, &
    0, 0, 4 /), (/ m, 7 /) )
  integer j
  integer k
  integer n
  real ( kind = rk ), parameter :: pi = 3.1415926535897932384626434D+00
  real ( kind = rk ) result(7)
  integer seed
  real ( kind = rk ) shrink
  real ( kind = rk ) v1(m)
  real ( kind = rk ) v2(m)
  real ( kind = rk ) v3(m)
  real ( kind = rk ) wc(m)
  real ( kind = rk ) w1(m)
  real ( kind = rk ) w2(m)
  real ( kind = rk ) w3(m)
  real ( kind = rk ), allocatable :: value(:)
  real ( kind = rk ), allocatable :: x(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Estimate monomial integrals over a sphere triangle'
  write ( *, '(a)' ) '  using the Monte Carlo method.'

  seed = 123456789
!
!  Choose three points at random to define a spherical triangle.
!
  call sphere01_sample ( 1, seed, w1 )
  call sphere01_sample ( 1, seed, w2 )
  call sphere01_sample ( 1, seed, w3 )

  wc(1:m) = ( w1(1:m) + w2(1:m) + w3(1:m) ) / 3.0D+00
  call r8vec_normalize ( m, wc )
!
!  Shrink triangle by factor F.
!
  shrink = 2.0D+00

  do k = 1, 3

    shrink = shrink / 2.0D+00

    v1(1:m) = wc(1:m) + shrink * ( w1(1:m) - wc(1:m) )
    v2(1:m) = wc(1:m) + shrink * ( w2(1:m) - wc(1:m) )
    v3(1:m) = wc(1:m) + shrink * ( w3(1:m) - wc(1:m) )

    call r8vec_normalize ( m, v1 )
    call r8vec_normalize ( m, v2 )
    call r8vec_normalize ( m, v3 )

    call sphere01_triangle_vertices_to_area ( v1, v2, v3, area )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Vertices of random spherical triangle'
    write ( *, '(a,g14.6)' ) '  with shrink factor = ', shrink
    write ( *, '(a,g14.6)' ) '  and area = ', area
    write ( *, '(a)' ) ' '
    write ( *, '(a,3g14.6)' ) '  V1:', v1(1:m)
    write ( *, '(a,3g14.6)' ) '  V2:', v2(1:m)
    write ( *, '(a,3g14.6)' ) '  V3:', v3(1:m)
!
!  Estimate integrals.
!
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '         N        1              X^2             Y^2' // &
      '             Z^2             X^4           X^2Y^2           Z^4'
    write ( *, '(a)' ) ' '

    n = 1

    do while ( n <= 4 * 65536 )

      allocate ( value(1:n) )
      allocate ( x(1:m,1:n) )

      call sphere01_triangle_sample ( n, v1, v2, v3, seed, x )

      do j = 1, 7

        e(1:m) = e_test(1:m,j)

        call monomial_value ( m, n, e, x, value )

        result(j) = area * sum ( value(1:n) ) / real ( n, kind = rk )

      end do

      write ( *, '(2x,i8,7(2x,g14.6))' ) n, result(1:7)

      deallocate ( value )
      deallocate ( x )

      n = 2 * n

    end do

  end do

  return
end
