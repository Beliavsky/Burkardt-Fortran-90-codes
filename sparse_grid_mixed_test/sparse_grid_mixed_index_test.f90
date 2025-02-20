program main

!****************************************************************************80
!
!! sparse_grid_mixed_index_test() tests sparse_grid_mixed_index().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 December 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Fabio Nobile, Raul Tempone, Clayton Webster,
!    A Sparse Grid Stochastic Collocation Method for Partial Differential
!    Equations with Random Input Data,
!    SIAM Journal on Numerical Analysis,
!    Volume 46, Number 5, 2008, pages 2309-2345.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r8_epsilon
  real ( kind = rk ) tol

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSE_GRID_MIXED_INDEX_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version'

  tol = r8_epsilon ( )
  call sparse_grid_mixed_index_tests ( tol )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSE_GRID_MIXED_INDEX_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine sparse_grid_mixed_index_tests ( tol )

!****************************************************************************80
!
!! SPARSE_GRID_MIXED_INDEX_TESTS calls SPARSE_GRID_MIXED_INDEX_TEST with various arguments.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) TOL, a tolerance for point equality.
!    A value of sqrt ( eps ) is reasonable, and will allow the code to
!    consolidate points which are equal, or very nearly so.  A value of
!    -1.0, on the other hand, will force the code to use every point, regardless
!    of duplication.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable, dimension ( : ) :: alpha
  real ( kind = rk ), allocatable, dimension ( : ) :: beta
  integer dim_num
  integer level_max_max
  integer level_max_min
  integer, allocatable, dimension ( : ) :: order_1d
  integer order_nd
  integer, allocatable, dimension ( : ) :: rule
  real ( kind = rk ) tol

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSE_GRID_MIXED_INDEX_TESTS'
  write ( *, '(a)' ) '  Call SPARSE_GRID_MIXED_INDEX_TEST with various arguments.'
  write ( *, '(a,g14.6)' ) '  All tests will use a point equality tolerance of ', tol

  dim_num = 2
  level_max_min = 0
  level_max_max = 2
  allocate ( alpha(1:dim_num) )
  allocate ( beta(1:dim_num) )
  allocate ( rule(1:dim_num) )
  rule = (/ 1, 1 /)
  alpha = (/ 0.0D+00, 0.0D+00 /)
  beta = (/ 0.0D+00, 0.0D+00 /)
  call sparse_grid_mixed_index_test ( dim_num, level_max_min, level_max_max, rule, &
    alpha, beta, tol )
  deallocate ( alpha )
  deallocate ( beta )
  deallocate ( rule )

  dim_num = 2
  level_max_min = 0
  level_max_max = 2
  allocate ( alpha(1:dim_num) )
  allocate ( beta(1:dim_num) )
  allocate ( rule(1:dim_num) )
  rule = (/ 1, 3 /)
  alpha = (/ 0.0D+00, 0.0D+00 /)
  beta = (/ 0.0D+00, 0.0D+00 /)
  call sparse_grid_mixed_index_test ( dim_num, level_max_min, level_max_max, rule, &
    alpha, beta, tol )
  deallocate ( alpha )
  deallocate ( beta )
  deallocate ( rule )

  dim_num = 2
  level_max_min = 0
  level_max_max = 2
  allocate ( alpha(1:dim_num) )
  allocate ( beta(1:dim_num) )
  allocate ( rule(1:dim_num) )
  rule = (/ 1, 4 /)
  alpha = (/ 0.0D+00, 0.0D+00 /)
  beta = (/ 0.0D+00, 0.0D+00 /)
  call sparse_grid_mixed_index_test ( dim_num, level_max_min, level_max_max, rule, &
    alpha, beta, tol )
  deallocate ( alpha )
  deallocate ( beta )
  deallocate ( rule )

  dim_num = 2
  level_max_min = 0
  level_max_max = 2
  allocate ( alpha(1:dim_num) )
  allocate ( beta(1:dim_num) )
  allocate ( rule(1:dim_num) )
  rule = (/ 1, 7 /)
  alpha = (/ 0.0D+00, 0.0D+00 /)
  beta = (/ 0.0D+00, 0.0D+00 /)
  call sparse_grid_mixed_index_test ( dim_num, level_max_min, level_max_max, rule, &
    alpha, beta, tol )
  deallocate ( alpha )
  deallocate ( beta )
  deallocate ( rule )

  dim_num = 2
  level_max_min = 0
  level_max_max = 2
  allocate ( alpha(1:dim_num) )
  allocate ( beta(1:dim_num) )
  allocate ( rule(1:dim_num) )
  rule = (/ 1, 8 /)
  alpha = (/ 0.0D+00, 1.5D+00 /)
  beta = (/ 0.0D+00, 0.0D+00 /)
  call sparse_grid_mixed_index_test ( dim_num, level_max_min, level_max_max, rule, &
    alpha, beta, tol )
  deallocate ( alpha )
  deallocate ( beta )
  deallocate ( rule )

  dim_num = 2
  level_max_min = 0
  level_max_max = 2
  allocate ( alpha(1:dim_num) )
  allocate ( beta(1:dim_num) )
  allocate ( rule(1:dim_num) )
  rule = (/ 2, 9 /)
  alpha = (/ 0.0D+00, 0.5D+00 /)
  beta = (/ 0.0D+00, 1.5D+00 /)
  call sparse_grid_mixed_index_test ( dim_num, level_max_min, level_max_max, rule, &
    alpha, beta, tol )
  deallocate ( alpha )
  deallocate ( beta )
  deallocate ( rule )

  dim_num = 2
  level_max_min = 0
  level_max_max = 2
  allocate ( alpha(1:dim_num) )
  allocate ( beta(1:dim_num) )
  allocate ( rule(1:dim_num) )
  rule = (/ 6, 4 /)
  alpha = (/ 2.0D+00, 0.0D+00 /)
  beta = (/ 0.0D+00, 0.0D+00 /)
  call sparse_grid_mixed_index_test ( dim_num, level_max_min, level_max_max, rule, &
    alpha, beta, tol )
  deallocate ( alpha )
  deallocate ( beta )
  deallocate ( rule )

  dim_num = 3
  level_max_min = 0
  level_max_max = 2
  allocate ( alpha(1:dim_num) )
  allocate ( beta(1:dim_num) )
  allocate ( rule(1:dim_num) )
  rule = (/ 1, 2, 5 /)
  alpha = (/ 0.0D+00, 0.0D+00, 0.0D+00 /)
  beta = (/ 0.0D+00, 0.0D+00, 0.0D+00 /)
  call sparse_grid_mixed_index_test ( dim_num, level_max_min, level_max_max, rule, &
    alpha, beta, tol )
  deallocate ( alpha )
  deallocate ( beta )
  deallocate ( rule )
!
!  Dimension 2, Level 3, Rule 13
!
  dim_num = 2
  level_max_min = 3
  level_max_max = 3
  allocate ( alpha(1:dim_num) )
  allocate ( beta(1:dim_num) )
  allocate ( rule(1:dim_num) )
  rule = (/ 13, 13 /)
  alpha = (/ 0.0D+00, 0.0D+00 /)
  beta = (/ 0.0D+00, 0.0D+00 /)
  call sparse_grid_mixed_index_test ( dim_num, level_max_min, level_max_max, rule, &
    alpha, beta, tol )
  deallocate ( alpha )
  deallocate ( beta )
  deallocate ( rule )
!
!  Dimension 2, Level 3, Rule 16
!
  dim_num = 2
  level_max_min = 3
  level_max_max = 3
  allocate ( alpha(1:dim_num) )
  allocate ( beta(1:dim_num) )
  allocate ( rule(1:dim_num) )
  rule = (/ 16, 16 /)
  alpha = (/ 0.0D+00, 0.0D+00 /)
  beta = (/ 0.0D+00, 0.0D+00 /)
  call sparse_grid_mixed_index_test ( dim_num, level_max_min, level_max_max, rule, &
    alpha, beta, tol )
  deallocate ( alpha )
  deallocate ( beta )
  deallocate ( rule )
!
!  Dimension 2, Level 3, Rule 17
!
  dim_num = 2
  level_max_min = 3
  level_max_max = 3
  allocate ( alpha(1:dim_num) )
  allocate ( beta(1:dim_num) )
  allocate ( rule(1:dim_num) )
  rule = (/ 17, 17 /)
  alpha = (/ 0.0D+00, 0.0D+00 /)
  beta = (/ 0.0D+00, 0.0D+00 /)
  call sparse_grid_mixed_index_test ( dim_num, level_max_min, level_max_max, rule, &
    alpha, beta, tol )
  deallocate ( alpha )
  deallocate ( beta )
  deallocate ( rule )

  return
end
subroutine sparse_grid_mixed_index_test ( dim_num, level_max_min, &
  level_max_max, rule, alpha, beta, tol )

!****************************************************************************80
!
!! SPARSE_GRID_MIXED_INDEX_TEST tests SPARSE_GRID_MIXED_INDEX.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer LEVEL_MAX_MIN, LEVEL_MAX_MAX, the minimum and
!    maximum values of LEVEL_MAX.
!
!    Input, integer RULE(DIM_NUM), the rule in each dimension.
!     1, "CC",  Clenshaw Curtis, Closed Fully Nested.
!     2, "F2",  Fejer Type 2, Open Fully Nested.
!     3, "GP",  Gauss Patterson, Open Fully Nested.
!     4, "GL",  Gauss Legendre, Open Weakly Nested.
!     5, "GH",  Gauss Hermite, Open Weakly Nested.
!     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested.
!     7, "LG",  Gauss Laguerre, Open Non Nested.
!     8, "GLG", Generalized Gauss Laguerre, Open Non Nested.
!     9, "GJ",  Gauss Jacobi, Open Non Nested.
!    10, "GW",  Golub Welsch, (presumed) Open Non Nested.
!    11, "CC_SE", Clenshaw Curtis Slow Exponential, Closed Fully Nested.
!    12, "F2_SE", Fejer Type 2 Slow Exponential, Closed Fully Nested.
!    13, "GP_SE", Gauss Patterson Slow Exponential, Closed Fully Nested.
!    14, "CC_ME", Clenshaw Curtis Moderate Exponential, Closed Fully Nested.
!    15, "F2_ME", Fejer Type 2 Moderate Exponential, Closed Fully Nested.
!    16, "GP_ME", Gauss Patterson Moderate Exponential, Closed Fully Nested.
!    17, "CCN", Clenshaw Curtis Nested, Linear, Closed Fully Nested rule.
!
!    Input, real ( kind = rk ) ALPHA(DIM_NUM), BETA(DIM_NUM), parameters used for
!    Generalized Gauss Hermite, Generalized Gauss Laguerre, and Gauss Jacobi rules.
!
!    Input, real ( kind = rk ) TOL, a tolerance for point equality.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num

  real ( kind = rk ) alpha(dim_num)
  real ( kind = rk ) beta(dim_num)
  integer dim
  integer level_max
  integer level_max_max
  integer level_max_min
  integer point
  integer point_num
  integer point_total_num
  integer rule(dim_num)
  integer, allocatable, dimension ( :, : ) :: sparse_index
  integer, allocatable, dimension ( :, : ) :: sparse_order
  integer, allocatable, dimension ( : ) :: sparse_unique_index
  real ( kind = rk ) tol

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSE_GRID_MIXED_INDEX_TEST'
  write ( *, '(a)' ) '  SPARSE_GRID_MIXED_INDEX returns index and order vectors that'
  write ( *, '(a)' ) '  identify each point in a multidimensional sparse grid with mixed factors.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Each sparse grid is of spatial dimension DIM_NUM,'
  write ( *, '(a)' ) '  and is made up of product grids of levels up to LEVEL_MAX.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' Dimension      Rule     Alpha          Beta'
  write ( *, '(a)' ) ' '

  do dim = 1, dim_num
    if ( rule(dim) ==  1 ) then
      write ( *, '(2x,i8,2x,i8)' ) &
        dim, rule(dim)
    else if ( rule(dim) ==  2 ) then
      write ( *, '(2x,i8,2x,i8)' ) &
        dim, rule(dim)
    else if ( rule(dim) ==  3 ) then
      write ( *, '(2x,i8,2x,i8)' ) &
        dim, rule(dim)
    else if ( rule(dim) ==  4 ) then
      write ( *, '(2x,i8,2x,i8)' ) &
        dim, rule(dim)
    else if ( rule(dim) ==  5 ) then
      write ( *, '(2x,i8,2x,i8)' ) &
        dim, rule(dim)
    else if ( rule(dim) ==  6 ) then
      write ( *, '(2x,i8,2x,i8,2x,g14.6)' ) &
        dim, rule(dim), alpha(dim)
    else if ( rule(dim) ==  7 ) then
      write ( *, '(2x,i8,2x,i8)' ) &
        dim, rule(dim)
    else if ( rule(dim) ==  8 ) then
      write ( *, '(2x,i8,2x,i8,2x,g14.6)' ) &
        dim, rule(dim), alpha(dim)
    else if ( rule(dim) ==  9 ) then
      write ( *, '(2x,i8,2x,i8,2x,g14.6,2x,g14.6)' ) &
        dim, rule(dim), alpha(dim), beta(dim)
    else if ( rule(dim) ==  10 ) then
      write ( *, '(2x,i8,2x,i8)' ) &
        dim, rule(dim)
    else if ( rule(dim) ==  11 ) then
      write ( *, '(2x,i8,2x,i8)' ) &
        dim, rule(dim)
    else if ( rule(dim) ==  12 ) then
      write ( *, '(2x,i8,2x,i8)' ) &
        dim, rule(dim)
    else if ( rule(dim) ==  13 ) then
      write ( *, '(2x,i8,2x,i8)' ) &
        dim, rule(dim)
    else if ( rule(dim) ==  14 ) then
      write ( *, '(2x,i8,2x,i8)' ) &
        dim, rule(dim)
    else if ( rule(dim) ==  15 ) then
      write ( *, '(2x,i8,2x,i8)' ) &
        dim, rule(dim)
    else if ( rule(dim) ==  16 ) then
      write ( *, '(2x,i8,2x,i8)' ) &
        dim, rule(dim)
    else if ( rule(dim) ==  17 ) then
      write ( *, '(2x,i8,2x,i8)' ) &
        dim, rule(dim)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPARSE_GRID_MIXED_INDEX_TEST - Fatal error!'
      write ( *, '(a,i8)' ) '  Unexpected value of RULE = ', rule(dim)
      stop
    end if
  end do

  do level_max = level_max_min, level_max_max

    call sparse_grid_mixed_size_total ( dim_num, level_max, rule, &
      point_total_num )

    call sparse_grid_mixed_size ( dim_num, level_max, rule, alpha, beta, &
      tol, point_num )

    allocate ( sparse_unique_index(1:point_total_num) )

    call sparse_grid_mixed_unique_index ( dim_num, level_max, rule, alpha, &
      beta, tol, point_num, point_total_num, sparse_unique_index )

    allocate ( sparse_index(1:dim_num,1:point_num) )
    allocate ( sparse_order(1:dim_num,1:point_num) )

    call sparse_grid_mixed_index ( dim_num, level_max, rule, point_num, &
      point_total_num, sparse_unique_index, sparse_order, sparse_index )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  For LEVEL_MAX = ', level_max
    write ( *, '(a)' ) ' '
    do point = 1, point_num
      write ( *, '(2x,i4,2x,10(2x,i3, a2, i3 ))' ) &
        point, ( sparse_index(dim,point), ' /', &
        sparse_order(dim,point), dim = 1, dim_num )
    end do

    deallocate ( sparse_index )
    deallocate ( sparse_order )
    deallocate ( sparse_unique_index )

  end do

  return
end

