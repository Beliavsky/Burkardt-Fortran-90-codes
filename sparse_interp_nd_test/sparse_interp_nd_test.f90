program main

!*****************************************************************************80
!
!! SPARSE_INTERP_ND_TEST() tests SPARSE_INTERP_ND().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 October 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer sparse_max

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSE_INTERP_ND_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) '  Test the SPARSE_INTERP_ND library.'
  write ( *, '(a)' ) '  The R8LIB library is also required.'

  m = 1
  sparse_max = 9
  call test01 ( m, sparse_max )

  m = 2
  sparse_max = 9
  call test01 ( m, sparse_max )

  m = 3
  sparse_max = 9
  call test01 ( m, sparse_max )

  m = 4
  sparse_max = 7
  call test01 ( m, sparse_max )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSE_INTERP_ND_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( m, sparse_max )

!*****************************************************************************80
!
!! TEST01: sequence of sparse interpolants to an M-dimensional function.
!
!  Discussion:
!
!    We have functions that can generate a Lagrange interpolant to data
!    in M dimensions, with specified order or level in each dimension.
!
!    We use the Lagrange function as the inner evaluator for a sparse
!    grid procedure. 
!
!    The procedure computes sparse interpolants of levels 0 to SPARSE_MAX
!    to a given function.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Local, integer M, the spatial dimension.
!
!    Input, integer SPARSE_MAX, the maximum sparse grid level to try.
!
!  Local:
!
!    Local, real ( kind = rk ) A(M), B(M), the upper and lower variable limits 
!    in each dimension.
!
!    Local, real ( kind = rk ) APP_ERROR, the averaged Euclidean norm of the 
!    difference between the sparse interpolant and the exact function at 
!    the interpolation points.
!
!    Local, integer C(0:L_MAX), the sparse grid coefficient vector.
!    Results at level L are weighted by C(L).
!
!    Local, integer IND(M), the 1D indices defining a Lagrange grid.
!    Each entry is a 1d "level" that specifies the order of a 
!    Clenshaw Curtis 1D grid.
!
!    Local, integer L, the current Lagrange grid level.
!
!    Local, integer L_MAX, the current sparse grid level.
!
!    Local, integer MORE, used to control the enumeration of all the
!    Lagrange grids at a current grid level.
!
!    Local, integer ND, the number of points used in a Lagrange grid.
!
!    Local, integer ND_TOTAL, the total number of points used in all the
!    Lagrange interpolants for a given sparse interpolant points that occur
!    in more than one Lagrange grid are counted multiple times.
!
!    Local, integer NI, the number of interpolant evaluation points.
!
!    Local, integer SPARSE_MIN, the minimum sparse grid level to try.
!
!    Local, real ( kind = rk ) XD(M,ND), the data points for a Lagrange grid.
!
!    Local, real ( kind = rk ) XI(M,NI), the interpolant evaluation points.
!
!    Local, real ( kind = rk ) ZD(ND), the data values for a Lagrange grid.
!
!    Local, real ( kind = rk ) ZE(NI), the exact function values at XI.
!
!    Local, real ( kind = rk ) ZI(NI), the sparse interpolant values at XI.
!
!    Local, real ( kind = rk ) ZPI(NI), one set of Lagrange interpolant values at XI.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: a(:)
  real ( kind = rk ) app_error
  real ( kind = rk ), allocatable :: b(:)
  integer, allocatable :: c(:)
  integer h
  integer, allocatable :: ind(:)
  integer l
  integer l_max
  integer l_min
  integer m
  logical more
  integer nd
  integer nd_total
  integer ni
  real ( kind = rk ) r8vec_norm_affine
  integer sparse_max
  integer sparse_min
  integer t
  integer, allocatable :: w(:)
  real ( kind = rk ), allocatable :: xd(:,:)
  real ( kind = rk ), allocatable :: xi(:,:)
  real ( kind = rk ), allocatable :: zd(:)
  real ( kind = rk ), allocatable :: ze(:)
  real ( kind = rk ), allocatable :: zi(:)
  real ( kind = rk ), allocatable :: zpi(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  Sparse interpolation for a function f(x) of M-dimensional argument.'
  write ( *, '(a)' ) '  Use a sequence of sparse grids of levels 0 through SPARSE_MAX.'
  write ( *, '(a)' ) '  Invoke a general Lagrange interpolant function to do this.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Compare the exact function and the interpolants at a grid of points.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  The "order" is the sum of the orders of all the product grids'
  write ( *, '(a)' ) '  used to make a particular sparse grid.'
!
!  User input.
!
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Spatial dimension M = ', m
  write ( *, '(a,i4)' ) '  Maximum sparse grid level = ', sparse_max

  allocate ( ind(1:m) )
!
!  Define the region.
!
  allocate ( a(1:m) )
  allocate ( b(1:m) )

  a(1:m) = 0.0D+00
  b(1:m) = 1.0D+00
!
!  Define the interpolation evaluation information.
!
  ni = 100

  allocate ( xi(m,ni) )
  call r8mat_uniform_abvec ( m, ni, a, b, xi )

  write ( *, '(a,i6)' ) '  Number of interpolation points is NI = ', ni

  allocate ( ze(1:ni) )
  call f_sinr ( m, ni, xi, ze )
!
!  Compute a sequence of sparse grid interpolants of increasing level.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '   L     Order    ApproxError'
  write ( *, '(a)' ) ''

  allocate ( zpi(1:ni) )
  allocate ( zi(1:ni) )

  sparse_min = 0

  do l_max = sparse_min, sparse_max

    allocate ( c(0:l_max) )
    allocate ( w(0:l_max) )
    call smolyak_coefficients ( l_max, m, c, w )

    zi(1:ni) = 0.0D+00
    nd_total = 0

    l_min = max ( l_max + 1 - m, 0 )

    do l = l_min, l_max

      more = .false.

      do
!
!  Define the next product grid.
!
        call comp_next ( l, m, ind, more, h, t )
!
!  Count the grid, find its points, evaluate the data there.
!
        call lagrange_interp_nd_size2 ( m, ind, nd )
        allocate ( xd(1:m,1:nd) )
        call lagrange_interp_nd_grid2 ( m, ind, a, b, nd, xd )
        allocate ( zd(1:nd) )
        call f_sinr ( m, nd, xd, zd )
!
!  Use the grid to evaluate the interpolant.
!
        call lagrange_interp_nd_value2 ( m, ind, a, b, nd, zd, ni, xi, zpi )
!
!  Weighted the interpolant values and add to the sparse grid interpolant.
!
        nd_total = nd_total + nd
        zi(1:ni) = zi(1:ni) + c(l) * zpi(1:ni)

        deallocate ( xd )
        deallocate ( zd )

        if ( .not. more ) then
          exit
        end if

      end do

    end do
!
!  Compare sparse interpolant and exact function at interpolation points.
!
    app_error = r8vec_norm_affine ( ni, zi, ze ) / real ( ni, kind = rk )

    write ( *, '(2x,i2,2x,i8,2x,e8.2)' ) l, nd_total, app_error

    deallocate ( c )
    deallocate ( w )

  end do

  deallocate ( a )
  deallocate ( b )
  deallocate ( ind )
  deallocate ( xi )
  deallocate ( ze )
  deallocate ( zi )
  deallocate ( zpi )

  return
end
subroutine f_sinr ( m, n, x, z )

!*****************************************************************************80
!
!! F_SINR is a scalar function of an M-dimensional argument, to be interpolated.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) X(M,N), the points.
!
!    Output, real ( kind = rk ) Z(N), the value of the function at each point.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) r(n)
  real ( kind = rk ) x(m,n)
  real ( kind = rk ) z(n)
      
  r(1:n) = sqrt ( sum ( x(1:m,1:n)**2, dim = 1 ) )
  z(1:n) = sin ( r(1:n) )

  return
end
