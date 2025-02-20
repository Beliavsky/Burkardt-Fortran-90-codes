program main

!*****************************************************************************80
!
!! PWL_INTERP_2D_SCATTERED_TEST tests the PWL_INTERP_2D_SCATTERED_TEST library.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 October 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer prob
  integer prob_num

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PWL_INTERP_2D_SCATTERED_TEST:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test PWL_INTERP_2D_SCATTERED.'
  write ( *, '(a)' ) '  The R8LIB library is needed.'
  write ( *, '(a)' ) '  This test also needs the TEST_INTERP_2D library.'

  call test01 ( )
  call test02 ( )
!
!  Numerical tests.
!
  call f00_num ( prob_num )

  do prob = 1, prob_num
    call test03 ( prob )
  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PWL_INTERP_2D_SCATTERED_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests R8TRIS2.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2
  integer, parameter :: node_num = 9
  integer, parameter :: element_order = 3

  real ( kind = rk ), dimension (dim_num,node_num) :: node_xy = reshape ( (/ &
       0.0D+00, 0.0D+00, &
       0.0D+00, 1.0D+00, &
       0.2D+00, 0.5D+00, &
       0.3D+00, 0.6D+00, &
       0.4D+00, 0.5D+00, &
       0.6D+00, 0.4D+00, &
       0.6D+00, 0.5D+00, &
       1.0D+00, 0.0D+00, &
       1.0D+00, 1.0D+00 /), (/ dim_num, node_num /) )
  integer element_num
  integer triangle(element_order,2*node_num)
  integer element_neighbor(3,2*node_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  R8TRIS2 computes the Delaunay triangulation of'
  write ( *, '(a)' ) '  a set of nodes in 2D.'
!
!  Set up the Delaunay triangulation.
!
  call r8tris2 ( node_num, node_xy, element_num, triangle, element_neighbor )

  call triangulation_order3_print ( node_num, element_num, node_xy, &
    triangle, element_neighbor )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests PWL_INTERP_2D_SCATTERED_VALUE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2
  integer, parameter :: ni = 25
  integer, parameter :: node_num = 9
  integer, parameter :: element_order = 3

  integer element_neighbor(3,2*node_num)
  integer element_num
  integer i
  integer j
  integer k
  real ( kind = rk ), dimension (dim_num,node_num) :: node_xy = reshape ( (/ &
       0.0D+00, 0.0D+00, &
       0.0D+00, 1.0D+00, &
       0.2D+00, 0.5D+00, &
       0.3D+00, 0.6D+00, &
       0.4D+00, 0.5D+00, &
       0.6D+00, 0.4D+00, &
       0.6D+00, 0.5D+00, &
       1.0D+00, 0.0D+00, &
       1.0D+00, 1.0D+00 /), (/ dim_num, node_num /) )
  integer triangle(element_order,2*node_num)
  real ( kind = rk ) x
  real ( kind = rk ) xyi(2,ni)
  real ( kind = rk ) y
  real ( kind = rk ) zd(node_num)
  real ( kind = rk ) zi(ni)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  PWL_INTERP_2D_SCATTERED_VALUE evaluates a'
  write ( *, '(a)' ) '  piecewise linear interpolant to scattered data.'
!
!  Set up the Delaunay triangulation.
!
  call r8tris2 ( node_num, node_xy, element_num, triangle, element_neighbor )
!
!  Define the Z data.
!
  do i = 1, node_num
    x = node_xy(1,i)
    y = node_xy(2,i)
    zd(i) = x + 2.0D+00 * y
  end do
!
!  Define the interpolation points.
!
  k = 0
  do i = 1, 5
    do j = 1, 5
      k = k + 1
      xyi(1,k) = ( i - 1 ) / 4.0D+00
      xyi(2,k) = ( j - 1 ) / 4.0D+00
    end do
  end do
!
!  Evaluate the interpolant.
!
  call pwl_interp_2d_scattered_value ( node_num, node_xy, zd, element_num, &
    triangle, element_neighbor, ni, xyi, zi )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     K      Xi(K)       Yi(K)       Zi(K)       Z(X,Y)'
  write ( *, '(a)' ) ' '
  do k = 1, ni
    write ( *, '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4,2x,f10.4)' ) &
      k, xyi(1,k), xyi(2,k), zi(k), xyi(1,k) + 2.0 * xyi(2,k)
  end do

  return
end
subroutine test03 ( prob )

!*****************************************************************************80
!
!! TEST03 tests PWL_INTERP_2D_SCATTERED_VALUE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 October 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: ni = 25

  integer, allocatable :: element_neighbor(:,:)
  integer element_num
  integer g
  integer i
  integer j
  integer k
  integer nd
  integer prob
  real ( kind = rk ) rms
  integer, allocatable :: triangle(:,:)
  real ( kind = rk ), allocatable :: xd(:)
  real ( kind = rk ) xi(ni)
  real ( kind = rk ), allocatable :: xyd(:,:)
  real ( kind = rk ) xyi(2,ni)
  real ( kind = rk ), allocatable :: yd(:)
  real ( kind = rk ) yi(ni)
  real ( kind = rk ), allocatable :: zd(:)
  real ( kind = rk ) ze(ni)
  real ( kind = rk ) zi(ni)

  g = 2
  call g00_size ( g, nd )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  PWL_INTERP_2D_SCATTERED_VALUE evaluates a'
  write ( *, '(a)' ) '  piecewise linear interpolant to scattered data.'
  write ( *, '(a,i2)' ) '  Here, we use grid number ', g
  write ( *, '(a,i4,a)' ) '  with ', nd, ' scattered points in the unit square'
  write ( *, '(a,i2)' ) '  on problem ', prob
!
!  Get the data points and evaluate the function there.
!
  allocate ( xd(1:nd) )
  allocate ( yd(1:nd) )

  call g00_xy ( g, nd, xd, yd )

  allocate ( zd(1:nd) )

  call f00_f0 ( prob, nd, xd, yd, zd )

  allocate ( xyd(2,1:nd) )

  xyd(1,1:nd) = xd(1:nd)
  xyd(2,1:nd) = yd(1:nd)
!
!  Set up the Delaunay triangulation.
!
  allocate ( element_neighbor(3,2*nd) )
  allocate ( triangle(3,2*nd) )

  call r8tris2 ( nd, xyd, element_num, triangle, element_neighbor )
!
!  Define the interpolation points.
!
  k = 0
  do i = 1, 5
    do j = 1, 5
      k = k + 1
      xyi(1,k) = ( 2 * i - 1 ) / 10.0D+00
      xyi(2,k) = ( 2 * j - 1 ) / 10.0D+00
    end do
  end do

  xi(1:ni) = xyi(1,1:ni)
  yi(1:ni) = xyi(2,1:ni)

  call f00_f0 ( prob, ni, xi, yi, ze )
!
!  Evaluate the interpolant.
!
  call pwl_interp_2d_scattered_value ( nd, xyd, zd, element_num, &
    triangle, element_neighbor, ni, xyi, zi )

  rms = 0.0D+00
  do k = 1, ni
    rms = rms + ( zi(k) - ze(k) )**2
  end do
  rms = sqrt ( rms / real ( ni ) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,e10.2)' ) '  RMS error is ', rms

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     K      Xi(K)       Yi(K)       Zi(K)       Z(X,Y)'
  write ( *, '(a)' ) ' '

  do k = 1, ni
    write ( *, '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4,2x,f10.4)' ) &
      k, xyi(1,k), xyi(2,k), zi(k), ze(k)
  end do

  deallocate ( element_neighbor )
  deallocate ( triangle )
  deallocate ( xd )
  deallocate ( xyd )
  deallocate ( yd )
  deallocate ( zd )

  return
end
