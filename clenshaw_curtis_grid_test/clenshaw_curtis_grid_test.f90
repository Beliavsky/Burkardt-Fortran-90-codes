program main

!*****************************************************************************80
!
!! clenshaw_curtis_grid_test() tests clenshaw_curtis_grid().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'clenshaw_curtis_grid_test()'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  test clenshaw_curtis_grid().'

  call test005 ( )
  call test01 ( )
  call test015 ( )
  call test02 ( )
  call test025 ( )
  call test03 ( )
  call test035 ( )
  call test04 ( )
  call test045 ( )
  call test05 ( )
  call test06 ( )
  call test07 ( )
  call test08 ( )
  call test09 ( )
  call test10 ( )
  call test11 ( )
  call test12 ( )
  call test13 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'clenshaw_curtis_grid_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test005 ( )

!*****************************************************************************80
!
!! test005() calls cc_grid() for the 1D problem.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 1

  real ( kind = rk ), allocatable, dimension ( :, : ) :: grid_point
  integer i
  integer j
  integer, dimension ( dim_num ) :: order_1d
  integer order_nd

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test005():'
  write ( *, '(a)' ) '  cc_grid() returns a grid of Clenshaw-Curtis points.'
  write ( *, '(a)' ) '  Here, we simply call for grids in the 1D case'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension of grid = ', dim_num

  do i = 1, 10

    write ( *, '(a)' ) ' '

    order_1d(1) = i
    order_nd = product ( order_1d(1:dim_num) )

    allocate ( grid_point(dim_num,order_nd) )

    call cc_grid ( dim_num, order_1d, order_nd, grid_point )

    do j = 1, order_nd
      write ( *, '(2x,i8,2x,f14.6)' ) j, grid_point(1:dim_num,j)
    end do

    deallocate ( grid_point )

  end do

  return
end
subroutine test01 ( )

!*****************************************************************************80
!
!! test01() simply calls cc_grid() once.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ), allocatable, dimension ( :, : ) :: grid_point
  integer j
  integer, dimension ( dim_num ) :: order_1d  = (/ 3, 4, 2 /)
  integer order_nd
  integer q

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test01:'
  write ( *, '(a)' ) '  CC_GRID returns a grid of Clenshaw-Curtis points.'
  write ( *, '(a)' ) '  Here, we simply call for a specific grid.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension of grid = ', dim_num

  order_nd = product ( order_1d(1:dim_num) )

  allocate ( grid_point(dim_num,order_nd) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Total number of points in the grid = ', order_nd
  write ( *, '(a)' ) ' '

  call cc_grid ( dim_num, order_1d, order_nd, grid_point )

  j = 1
  q = sum ( order_1d(1:dim_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I         Q          Grid orders:'
  write ( *, '(a)' ) ' '

  write ( *, '(2x,i8,2x,i8,4x,i8,2x,i8,2x,i8)' ) j, q, order_1d(1:dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Grid points:'
  write ( *, '(a)' ) ' '

  do j = 1, order_nd
    write ( *, '(2x,i8,2x,3f14.6)' ) j, grid_point(1:dim_num,j)
  end do

  deallocate ( grid_point )

  return
end
subroutine test015 ( )

!*****************************************************************************80
!
!! test015() tests cc_grid_index().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  integer, allocatable, dimension ( :, : ) :: grid_index
  integer j
  integer, dimension ( dim_num ) :: order_1d  = (/ 3, 4, 2 /)
  integer order_nd
  integer q

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test015:'
  write ( *, '(a)' ) &
    '  CC_GRID_INDEX returns an indexed grid of Clenshaw-Curtis points.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension of grid = ', dim_num

  order_nd = product ( order_1d(1:dim_num) )

  allocate ( grid_index(dim_num,order_nd) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Total number of points in the grid = ', order_nd
  write ( *, '(a)' ) ' '

  call cc_grid_index ( dim_num, order_1d, order_nd, grid_index )

  j = 1
  q = sum ( order_1d(1:dim_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I         Q          Grid orders:'
  write ( *, '(a)' ) ' '

  write ( *, '(2x,i8,2x,i8,4x,i4,2x,i4,2x,i4)' ) j, q, order_1d(1:dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Grid indexed points:'
  write ( *, '(a)' ) ' '

  do j = 1, order_nd
    write ( *, '(2x,i8,2x,8x,4x,i4,2x,i4,2x,i4)' ) j, grid_index(1:dim_num,j)
  end do

  deallocate ( grid_index )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! test02() tests cc_grids_minmax() for all points on 2D grids for Q = 3 to 5.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  integer grid_num
  integer, allocatable, dimension ( :, : ) :: grid_order
  real ( kind = rk ), allocatable, dimension ( :, : ) :: grid_point
  integer j
  integer point_num
  integer q
  integer, parameter :: q_max = 5
  integer, parameter :: q_min = 3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test02:'
  write ( *, '(a)' ) '  CC_GRIDS_MINMAX returns all Clenshaw Curtis grids'
  write ( *, '(a)' ) '  whose Q value satisfies Q_MIN <= Q <= Q_MAX.'
  write ( *, '(a)' ) '  Here, Q is the sum of the orders of the 1D rules, and'
  write ( *, '(a,i8)' ) '  Q_MIN = ', q_min
  write ( *, '(a,i8)' ) '  Q_MAX = ', q_max
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension of grids = ', dim_num

  call cc_grids_minmax_size ( dim_num, q_min, q_max, grid_num, &
    point_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of grids = ', grid_num
  write ( *, '(a,i8)' ) '  Number of points in the grids = ', point_num
!
!  Allocate the space.
!
  allocate ( grid_order(1:dim_num,1:grid_num) )
  allocate ( grid_point(1:dim_num,1:point_num) )
!
!  Compute the orders and points.
!
  call cc_grids_minmax ( dim_num, q_min, q_max, grid_num, point_num, &
    grid_order, grid_point )
!
!  Now we're done.  Print the merged grid data.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I         Q          Grid orders:'
  write ( *, '(a)' ) ' '
  do j = 1, grid_num
    q = sum ( grid_order(1:dim_num,j) )
    write ( *, '(2x,i8,2x,i8,4x,i8,2x,i8,2x,i8)' ) j, q, grid_order(1:dim_num,j)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Grid points:'
  write ( *, '(a)' ) ' '
  do j = 1, point_num
    write ( *, '(2x,i8,2x,3f14.6)' ) j, grid_point(1:dim_num,j)
  end do

  deallocate ( grid_order )
  deallocate ( grid_point )

  return
end
subroutine test025 ( )

!*****************************************************************************80
!
!! test025() tests cc_levels_minmax().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2
  integer, parameter :: test_num = 3

  integer grid_num
  integer, allocatable, dimension ( :, : ) :: grid_level
  integer, allocatable, dimension ( :, : ) :: grid_order
  real ( kind = rk ), allocatable, dimension ( :, : ) :: grid_point
  integer j
  integer level
  integer level_max
  integer, dimension ( test_num ) :: level_max_test = (/ 2, 3, 3 /)
  integer level_min
  integer, dimension ( test_num ) :: level_min_test = (/ 2, 0, 3 /)
  integer point_num
  integer test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test025:'
  write ( *, '(a)' ) '  CC_LEVELS_MINMAX returns all Clenshaw Curtis grids'
  write ( *, '(a)' ) '  whose level value satisfies'
  write ( *, '(a)' ) '    LEVEL_MIN <= LEVEL <= LEVEL_MAX.'
  write ( *, '(a)' ) '  Here, LEVEL is the sum of the levels of the 1D rules,'
  write ( *, '(a)' ) '  and the order of the rule is 2**LEVEL + 1.'

  do test = 1, test_num

    level_min = level_min_test(test)
    level_max = level_max_test(test)

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  LEVEL_MIN = ', level_min
    write ( *, '(a,i8)' ) '  LEVEL_MAX = ', level_max
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Spatial dimension of grids = ', dim_num

    call cc_levels_minmax_size ( dim_num, level_min, level_max, grid_num, &
      point_num )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Number of grids = ', grid_num
    write ( *, '(a,i8)' ) '  Number of points in the grids = ', point_num
!
!  Allocate the space.
!
    allocate ( grid_level(1:dim_num,1:grid_num) )
    allocate ( grid_order(1:dim_num,1:grid_num) )
    allocate ( grid_point(1:dim_num,1:point_num) )
!
!  Compute the orders and points.
!
    call cc_levels_minmax ( dim_num, level_min, level_max, grid_num, &
      point_num, grid_level, grid_order, grid_point )
!
!  Now we're done.  Print the merged grid data.
!
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) &
      '      Grid     Level           Grid Levels         Grid orders:'
    write ( *, '(a)' ) &
      '      ----     -----          ------------        ------------'
    write ( *, '(a)' ) ' '
    do j = 1, grid_num
      level = sum ( grid_level(1:dim_num,j) )
      write ( *, '(2x,i8,2x,i8,4x,i8,2x,i8,2x,i8,2x,i8,2x,i8)' ) &
        j, level, grid_level(1:dim_num,j), grid_order(1:dim_num,j)
    end do

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Grid points:'
    write ( *, '(a)' ) ' '
    do j = 1, point_num
      write ( *, '(2x,i8,2x,3f14.6)' ) j, grid_point(1:dim_num,j)
    end do

    deallocate ( grid_level )
    deallocate ( grid_order )
    deallocate ( grid_point )

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! test03() calls cc_grids_constrained() to collect constrained grids.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ), dimension ( dim_num ) :: alpha = (/ 2.0D+00, 3.0D+00 /)
  integer dim
  integer grid_num
  integer, allocatable, dimension ( :, : ) :: grid_order
  real ( kind = rk ), allocatable, dimension ( :, : ) :: grid_point
  integer j
  integer, dimension ( dim_num ) :: order_max
  integer, dimension ( dim_num ) :: order_min
  integer point_num
  real ( kind = rk ) q
  real ( kind = rk ), parameter :: q_max = 13.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test03:'
  write ( *, '(a)' ) '  CC_GRIDS_CONSTRAINED returns all Clenshaw Curtis grids'
  write ( *, '(a)' ) '  satisfying a set of constraints.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  ORDER(I), the order of the 1D rule in dimension I,'
  write ( *, '(a)' ) '  is constrained by '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    ORDER_MIN(I) <= ORDER(I) <= ORDER_MAX(I)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We also define the total weighted order Q'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Q = ALPHA(1) * ORDER(1) + ... + ALPHA(N) * ORDER(N)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  and further constrain our grids to satisfy'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f14.6)' ) '    Q <= Q_MAX = ', q_max
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension of grids = ', dim_num

  order_min(1:dim_num) = 1
  order_max(1:dim_num) = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' Dimension Order_min Order_max     Alpha'
  write ( *, '(a)' ) ' '

  do dim = 1, dim_num
    write ( *, '(2x,i8,2x,i8,2x,i8,2x,f14.6)' ) &
      dim, order_min(dim), order_max(dim), alpha(dim)
  end do

  call cc_grids_constrained_size ( dim_num, q_max, alpha, &
    order_min, order_max, grid_num, point_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of grids = ', grid_num
  write ( *, '(a,i8)' ) '  Number of points in the grids = ', point_num
!
!  Allocate the space.
!
  allocate ( grid_point(1:dim_num,1:point_num) )
  allocate ( grid_order(1:dim_num,1:grid_num) )

  call cc_grids_constrained ( dim_num, q_max, alpha, &
    order_min, order_max, grid_num, point_num, grid_order, grid_point )
!
!  Now we're done.  Print the merged grid data.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I         Q          Grid orders:'
  write ( *, '(a)' ) ' '
  do j = 1, grid_num

    q = dot_product ( alpha(1:dim_num), grid_order(1:dim_num,j) )

    write ( *, '(2x,i8,2x,f14.6,4x,i8,2x,i8,2x,i8)' ) &
      j, q, grid_order(1:dim_num,j)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Grid points:'
  write ( *, '(a)' ) ' '
  do j = 1, point_num
    write ( *, '(2x,i8,2x,3f14.6)' ) j, grid_point(1:dim_num,j)
  end do

  deallocate ( grid_point )
  deallocate ( grid_order )

  return
end
subroutine test035 ( )

!*****************************************************************************80
!
!! test035() calls cc_levels_constrained() to collect constrained grids.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 September 202106 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ), dimension ( dim_num ) :: alpha = (/ 2.0D+00, 3.0D+00 /)
  integer dim
  integer, allocatable, dimension ( :, : ) :: grid_level
  integer grid_num
  real ( kind = rk ), allocatable, dimension ( :, : ) :: grid_point
  integer j
  integer, dimension ( dim_num ) :: level_max
  integer, dimension ( dim_num ) :: level_min
  integer point_num
  real ( kind = rk ) q
  real ( kind = rk ), parameter :: q_max = 13.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test035:'
  write ( *, '(a)' ) '  CC_LEVELS_CONSTRAINED returns all Clenshaw Curtis grids'
  write ( *, '(a)' ) '  satisfying a set of constraints.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The constraint on the levels of the 1D Clenshaw Curtis'
  write ( *, '(a)' ) '  rule in spatial dimension I is:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    LEVEL_MIN(I) <= LEVEL(I) <= LEVEL_MAX(I) '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The constraint on the levels making up a rule is:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Sum ( 1 <= I <= DIM_NUM ) ALPHA(I) * LEVEL(I) <= Q_MAX.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  where Q_MAX = ', q_max
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The relationship of level to order is roughly '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    ORDER = 2^LEVEL+1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension of grids = ', dim_num

  level_min(1:dim_num) = 1
  level_max(1:dim_num) = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' Dimension Level_min Level_max        Alpha'
  write ( *, '(a)' ) ' '

  do dim = 1, dim_num
    write ( *, '(2x,i8,2x,i8,2x,i8,2x,f14.6)' ) &
      dim, level_min(dim), level_max(dim), alpha(dim)
  end do

  call cc_levels_constrained_size ( dim_num, q_max, alpha, &
    level_min, level_max, grid_num, point_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of grids = ', grid_num
  write ( *, '(a,i8)' ) '  Number of points in the grids = ', point_num
!
!  Allocate the space.
!
  allocate ( grid_point(1:dim_num,1:point_num) )
  allocate ( grid_level(1:dim_num,1:grid_num) )

  call cc_levels_constrained ( dim_num, q_max, alpha, &
    level_min, level_max, grid_num, point_num, grid_level, grid_point )
!
!  Now we're done.  Print the merged grid data.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I         Q               Grid level:'
  write ( *, '(a)' ) ' '
  do j = 1, grid_num

    q = dot_product ( alpha(1:dim_num), grid_level(1:dim_num,j) )

    write ( *, '(2x,i8,2x,f14.6,4x,i8,2x,i8,2x,i8)' ) &
      j, q, grid_level(1:dim_num,j)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Grid points:'
  write ( *, '(a)' ) ' '
  do j = 1, point_num
    write ( *, '(2x,i8,2x,3f14.6)' ) j, grid_point(1:dim_num,j)
  end do

  deallocate ( grid_point )
  deallocate ( grid_level )

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! test04() tests clenshaw_curtis_compute().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 September 202106 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: order_max = 16

  integer i
  integer order
  real ( kind = rk ) w(order_max)
  real ( kind = rk ) x(order_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test04'
  write ( *, '(a)' ) '  CLENSHAW_CURTIS_COMPUTE computes'
  write ( *, '(a)' ) '  a Clenshaw-Curtis quadrature rule over [-1,1]'
  write ( *, '(a)' ) '  of given order.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Order       W               X'
  write ( *, '(a)' ) ' '

  do order = 1, 10

    call clenshaw_curtis_compute ( order, x, w )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,i8)' ) order

    do i = 1, order
      write ( *, '(10x,2x,g14.6,2x,g14.6)' ) w(i), x(i)
    end do

  end do

  return
end
subroutine test045 ( )

!*****************************************************************************80
!
!! test045() tests cc_abscissa() and cc_weight().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: order = 10

  real ( kind = rk ) cc_abscissa
  real ( kind = rk ) cc_weight
  integer i
  real ( kind = rk ) w
  real ( kind = rk ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test045'
  write ( *, '(a)' ) '  To compute a single Clenshaw Curtis weight or abscissa,'
  write ( *, '(a)' ) '  CC_ABSCISSA computes one abscissa,'
  write ( *, '(a)' ) '  CC_WEIGHT computes one weight.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We use these routines wastefully,'
  write ( *, '(a)' ) '  to compute the order 10 rule one value at a time.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Order       W               X'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,i8)' ) order

  do i = 1, order

    x = cc_abscissa ( order, i )
    w = cc_weight ( order, i )

    write ( *, '(10x,2x,g14.6,2x,g14.6)' ) w, x

  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! test05() tests clenshaw_curtis_compute().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: order_max = 16

  real ( kind = rk ), external :: f1
  real ( kind = rk ), external :: f2
  real ( kind = rk ), external :: f3
  integer i
  integer order
  real ( kind = rk ), dimension ( 3 ) :: result
  real ( kind = rk ) weight(order_max)
  real ( kind = rk ) xtab(order_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test05'
  write ( *, '(a)' ) '  CLENSHAW_CURTIS_COMPUTE computes a Clenshaw-Curtis rule;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The integration interval is [-1,1].'
  write ( *, '(a)' ) '  Quadrature order will vary.'
  write ( *, '(a)' ) '  Integrand will vary.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      Order     F1              F2              F3'
  write ( *, '(a)' ) ' '

  do order = 1, order_max

    call clenshaw_curtis_compute ( order, xtab, weight )
 
    result(1) = 0.0D+00
    do i = 1, order
      result(1) = result(1) + weight(i) * f1 ( xtab(i) )
    end do
 
    result(2) = 0.0D+00
    do i = 1, order
      result(2) = result(2) + weight(i) * f2 ( xtab(i) )
    end do

    result(3) = 0.0D+00
    do i = 1, order
      result(3) = result(3) + weight(i) * f3 ( xtab(i) )
    end do

    write ( *, '(2x,i6,2x,f14.8,2x,f14.8,2x,f14.8)' ) order, result(1:3)

  end do
 
  write ( *, '(a)' ) ' '

  result(1) = 46.0D+00 * sinh ( 1.0D+00 ) / 25.0D+00 - 2.0D+00 * sin ( 1.0D+00 )
  result(2) = 1.5822329637296729331D+00
  result(3) = ( sqrt ( 2.0D+00 ) + 3.0D+00 * sqrt ( 6.0D+00 ) ) / 6.0D+00

  write ( *, '(2x,a,2x,f14.8,2x,f14.8,2x,f14.8)' ) 'Exact ', result(1:3)

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! test06() tests clenshaw_curtis_set().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: order_max = 16

  real ( kind = rk ), external :: f1
  real ( kind = rk ), external :: f2
  real ( kind = rk ), external :: f3
  integer i
  integer order
  real ( kind = rk ), dimension ( 3 ) :: result
  real ( kind = rk ) weight(order_max)
  real ( kind = rk ) xtab(order_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test06'
  write ( *, '(a)' ) '  CLENSHAW_CURTIS_SET sets up a Clenshaw-Curtis rule;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The integration interval is [-1,1].'
  write ( *, '(a)' ) '  Quadrature order will vary.'
  write ( *, '(a)' ) '  Integrand will vary.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      Order     F1              F2              F3'
  write ( *, '(a)' ) ' '

  do order = 1, order_max

    call clenshaw_curtis_set ( order, xtab, weight )
 
    result(1) = 0.0D+00
    do i = 1, order
      result(1) = result(1) + weight(i) * f1 ( xtab(i) )
    end do
 
    result(2) = 0.0D+00
    do i = 1, order
      result(2) = result(2) + weight(i) * f2 ( xtab(i) )
    end do

    result(3) = 0.0D+00
    do i = 1, order
      result(3) = result(3) + weight(i) * f3 ( xtab(i) )
    end do

    write ( *, '(2x,i6,2x,f14.8,2x,f14.8,2x,f14.8)' ) order, result(1:3)

  end do
 
  write ( *, '(a)' ) ' '

  result(1) = 46.0D+00 * sinh ( 1.0D+00 ) / 25.0D+00 - 2.0D+00 * sin ( 1.0D+00 )
  result(2) = 1.5822329637296729331D+00
  result(3) = ( sqrt ( 2.0D+00 ) + 3.0D+00 * sqrt ( 6.0D+00 ) ) / 6.0D+00

  write ( *, '(2x,a6,2x,f14.8,2x,f14.8,2x,f14.8)' ) 'Exact ', result(1:3)

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! test07() tests r8table_write0().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 1
  integer, parameter :: order = 9

  character ( len = 100 ) :: r_file = 'cc_d1_o9_r.txt'
  real ( kind = rk ) r(dim_num,2)
  character ( len = 100 ) :: w_file = 'cc_d1_o9_w.txt'
  real ( kind = rk ) w(order)
  character ( len = 100 ) :: x_file = 'cc_d1_o9_x.txt'
  real ( kind = rk ) x(dim_num,order)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test07'
  write ( *, '(a)' ) '  r8table_write0() writes a Clenshaw-Curtis '
  write ( *, '(a)' ) '  quadrature rule to a file.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension = ', dim_num
  write ( *, '(a,i8)' ) '  Computing the rule of order = ', order

  call clenshaw_curtis_compute ( order, x, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Write abscissas to file "' // &
    trim ( x_file ) // '".'

  call r8table_write0 ( x_file, dim_num, order, x )

  write ( *, '(a)' ) '  Write weights to file "' // &
    trim ( w_file ) // '".'

  call r8table_write0 ( w_file, 1, order, w )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Write range to file "' // &
    trim ( r_file ) // '".'

  r(1:dim_num,1) = -1.0D+00
  r(1:dim_num,2) = +1.0D+00
  call r8table_write0 ( r_file, dim_num, 2, r )

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! test08() tests clenshaw_curtis_compute_nd().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num
  integer order
  integer, allocatable, dimension ( : ) :: order_1d
  integer order_nd
  real ( kind = rk ), allocatable, dimension ( :, : ) :: point
  real ( kind = rk ), allocatable, dimension ( : ) :: weight
  real ( kind = rk ) weight_sum

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test08'
  write ( *, '(a)' ) '  CLENSHAW_CURTIS_COMPUTE_ND computes'
  write ( *, '(a)' ) '  a multidimensional Clenshaw-Curtis quadrature rule'
  write ( *, '(a)' ) '  over the hypercube [-1,1]^ND of given'
  write ( *, '(a)' ) '  (possibly different) orders in each dimension.'

  dim_num = 2

  allocate ( order_1d(1:dim_num) )

  order_1d(1:dim_num) = 5

  order_nd = product ( order_1d(1:dim_num) )

  allocate ( point(1:dim_num,1:order_nd) )
  allocate ( weight(1:order_nd) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this example, we use the SAME ORDER'
  write ( *, '(a)' ) '  in all dimensions.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension DIM_NUM = ', dim_num
  write ( *, '(a,i8,i8)' ) '  1D orders = ', order_1d(1:dim_num)
  write ( *, '(a,i8)' ) '  Number of points = ', order_nd

  call clenshaw_curtis_compute_nd ( dim_num, order_1d, point, weight )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      Weight           X(1)           X(2)'
  write ( *, '(a)' ) '  --------------  --------------  --------------'
  write ( *, '(a)' ) ' '

  do order = 1, order_nd

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      weight(order), point(1:dim_num,order)

  end do

  weight_sum = sum ( weight(1:order_nd) )

  write ( *, '(a)' ) ' '
  write ( *, '(2x,g14.6)' ) weight_sum

  deallocate ( order_1d )
  deallocate ( point )
  deallocate ( weight )

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! test09() tests clenshaw_curtis_compute_nd().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num
  integer order
  integer, allocatable, dimension ( : ) :: order_1d
  integer order_nd
  real ( kind = rk ), allocatable, dimension ( :, : ) :: point
  real ( kind = rk ), allocatable, dimension ( : ) :: weight
  real ( kind = rk ) weight_sum

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test09'
  write ( *, '(a)' ) '  CLENSHAW_CURTIS_COMPUTE_ND computes'
  write ( *, '(a)' ) '  a multidimensional Clenshaw-Curtis quadrature rule'
  write ( *, '(a)' ) '  over the hypercube [-1,1]^ND of given'
  write ( *, '(a)' ) '  (possibly different) orders in each dimension.'

  dim_num = 3

  allocate ( order_1d(1:dim_num) )

  order_1d(1) = 2
  order_1d(2) = 4
  order_1d(3) = 3

  order_nd = product ( order_1d(1:dim_num) )

  allocate ( point(1:dim_num,1:order_nd) )
  allocate ( weight(1:order_nd) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this example, we use DIFFERENT ORDERS'
  write ( *, '(a)' ) '  in each dimension.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension DIM_NUM = ', dim_num
  write ( *, '(a,i8,i8,i8)' ) '  1D orders = ', order_1d(1:dim_num)
  write ( *, '(a,i8)' ) '  Number of points = ', order_nd

  call clenshaw_curtis_compute_nd ( dim_num, order_1d, point, weight )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      Weight           X(1)           X(2)           X(3)'
  write ( *, '(a)' ) '  --------------  --------------  --------------  --------------'
  write ( *, '(a)' ) ' '

  do order = 1, order_nd

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      weight(order), point(1:dim_num,order)

  end do

  weight_sum = sum ( weight(1:order_nd) )

  write ( *, '(a)' ) ' '
  write ( *, '(2x,g14.6)' ) weight_sum

  deallocate ( order_1d )
  deallocate ( point )
  deallocate ( weight )

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! test10() uses r8table_write0() to write out a multidimensional CC rule.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num
  integer, allocatable, dimension ( : ) :: order_1d
  integer order_nd
  real ( kind = rk ), allocatable, dimension ( :, : ) :: point
  real ( kind = rk ), allocatable, dimension ( :, : ) :: r
  character ( len = 100 ) :: r_file = 'cc_d4_o3x3x3x3_r.txt'
  character ( len = 100 ) :: w_file = 'cc_d4_o3x3x3x3_w.txt'
  real ( kind = rk ), allocatable, dimension ( : ) :: weight
  character ( len = 100 ) :: x_file = 'cc_d4_o3x3x3x3_x.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test10'
  write ( *, '(a)' ) '  Use r8table_write() to write out a multidimensional rule.'
  write ( *, '(a)' ) '  CLENSHAW_CURTIS_COMPUTE_ND computes'
  write ( *, '(a)' ) '  a multidimensional Clenshaw-Curtis quadrature rule'
  write ( *, '(a)' ) '  over the hypercube [-1,1]^ND of given'
  write ( *, '(a)' ) '  (possibly different) orders in each dimension.'

  dim_num = 4
 
  allocate ( order_1d(1:dim_num) )

  order_1d(1:dim_num) = 3

  order_nd = product ( order_1d(1:dim_num) )

  allocate ( point(1:dim_num,1:order_nd) )
  allocate ( r(1:dim_num,1:2) )
  allocate ( weight(1:order_nd) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this example, we use the SAME ORDER'
  write ( *, '(a)' ) '  in all dimensions.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension DIM_NUM = ', dim_num
  write ( *, '(a,i8,i8,i8,i8)' ) '  1D orders = ', order_1d(1:dim_num)
  write ( *, '(a,i8)' ) '  Number of points = ', order_nd

  call clenshaw_curtis_compute_nd ( dim_num, order_1d, point, weight )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Write abscissas to file "' // &
    trim ( x_file ) // '".'

  call r8table_write0 ( x_file, dim_num, order_nd, point )

  write ( *, '(a)' ) '  Write weights to file "' // &
    trim ( w_file ) // '".'

  call r8table_write0 ( w_file, 1, order_nd, weight )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Write range to file "' // &
    trim ( r_file ) // '".'

  r(1:dim_num,1) = -1.0D+00
  r(1:dim_num,2) = +1.0D+00

  call r8table_write0 ( r_file, dim_num, 2, r )

  deallocate ( order_1d )
  deallocate ( point )
  deallocate ( r )
  deallocate ( weight )

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! test11() tests cc_abscissa_level_1d().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer base
  integer i
  integer order
  integer, allocatable, dimension ( : ) :: test_level
  integer test_num
  integer, allocatable, dimension ( : ) :: test_val

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test11'
  write ( *, '(a)' ) '  CC_ABSCISSA_LEVEL_1D reports the level on which'
  write ( *, '(a)' ) '  a Clenshaw Curtis abscissa of given index will first'
  write ( *, '(a)' ) '  be generated, assuming a series of grids that grow'
  write ( *, '(a)' ) '  in order as 2**LEVEL+1.'

  base = 5
  order = 2**base + 1
  test_num = 2**base + 1

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Base B = ', base
  write ( *, '(a,i8)' ) '  ORDER 2^B+1 = ', order

  allocate ( test_val(1:test_num) )
  allocate ( test_level(1:test_num) )

  do i = 1, test_num
    test_val(i) = i - 1
  end do

  call cc_abscissa_level_1d ( base, test_num, test_val, test_level )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I  Level(I)'
  write ( *, '(a)' ) ' '

  do i = 1, test_num
    write ( *, '(2x,i8,2x,i8)' ) test_val(i), test_level(i)
  end do

  deallocate ( test_level )
  deallocate ( test_val )

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! test12() tests cc_abscissa_level_1d().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 September 202106 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer base
  integer i
  integer i4_modp
  integer order
  integer, allocatable, dimension ( : ) :: test_level
  integer test_num
  integer, allocatable, dimension ( : ) :: test_val

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test12'
  write ( *, '(a)' ) '  CC_ABSCISSA_LEVEL_1D can also be called for values'
  write ( *, '(a)' ) '  outside the standard range of 0 through 2**LEVEL_MAX.'
  write ( *, '(a)' ) '  In that case, a MOD operation is applied first,'
  write ( *, '(a)' ) '  to make a sensible result.'

  base = 5
  order = 2 ** base + 1
  test_num = 20

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Base B = ', base
  write ( *, '(a,i8)' ) '  ORDER = 2^B+1 = ', order

  allocate ( test_val(1:test_num) )
  allocate ( test_level(1:test_num) )

  call i4vec_uniform_ab ( test_num, -20, 100, test_val )

  call cc_abscissa_level_1d ( base, test_num, test_val, test_level )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I  Mod(I,O)  Level(I)'
  write ( *, '(a)' ) ' '

  do i = 1, test_num
    write ( *, '(2x,i8,2x,i8,2x,i8)' ) &
      test_val(i), i4_modp ( test_val(i), order ), test_level(i)
  end do

  deallocate ( test_level )
  deallocate ( test_val )

  return
end
subroutine test13 ( )

!*****************************************************************************80
!
!! test13() tests cc_abscissa_level_nd().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer base
  integer dim_num
  integer i
  integer j
  integer k
  integer order
  integer, allocatable, dimension ( : ) :: test_level
  integer test_num
  integer, allocatable, dimension ( :, : ) :: test_val

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test13'
  write ( *, '(a)' ) '  CC_ABSCISSA_LEVEL_ND reports the level on which'
  write ( *, '(a)' ) '  a Clenshaw Curtis abscissa of given index will first'
  write ( *, '(a)' ) '  be generated, assuming a series of grids that grow'
  write ( *, '(a)' ) '  in order as 2**LEVEL+1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This routine is applied for multidimensional cases.'

  base = 3
  order = 2**base + 1
  dim_num = 2
  test_num = order * order

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Base B = ', base
  write ( *, '(a,i8)' ) '  ORDER 2**B+1 = ', order
  write ( *, '(a,i8)' ) '  DIM_NUM = ', dim_num

  allocate ( test_val(1:dim_num,1:test_num) )
  allocate ( test_level(1:test_num) )

  k = 0
  do i = 0, order - 1
    do j = 0, order-1
      k = k + 1
      test_val(1:2,k) = (/ i, j /)
    end do
  end do

  call cc_abscissa_level_nd ( base, dim_num, test_num, test_val, test_level )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I         J  Level(I,J)'
  write ( *, '(a)' ) ' '

  do i = 1, test_num
    write ( *, '(2x,i8,2x,i8,2x,i8)' ) test_val(1:dim_num,i), test_level(i)
  end do

  deallocate ( test_level )
  deallocate ( test_val )

  return
end
function f1 ( x )

!*****************************************************************************80
!
!! f1() evaluates F1(X) = 23 * cosh ( x ) / 25 - cos ( x ).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Charles Clenshaw, Alan Curtis,
!    A Method for Numerical Integration on an Automatic Computer,
!    Numerische Mathematik,
!    Volume 2, Number 1, December 1960, pages 197-205.
!
!  Input:
!
!    real ( kind = rk ) X, the argument.
!
!  Output:
!
!    real ( kind = rk ) F1, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) f1
  real ( kind = rk ) x

  f1 = 23.0D+00 * cosh ( x ) / 25.0D+00 - cos ( x )

  return
end
function f2 ( x )

!*****************************************************************************80
!
!! f2() evaluates F2(X) = 1 / ( x^4 + x^2 + 0.9 ).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Charles Clenshaw, Alan Curtis,
!    A Method for Numerical Integration on an Automatic Computer,
!    Numerische Mathematik,
!    Volume 2, Number 1, December 1960, pages 197-205.
!
!  Input:
!
!    real ( kind = rk ) X, the argument.
!
!  Output:
!
!    real ( kind = rk ) F2, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) f2
  real ( kind = rk ) x

  f2 = 1.0D+00 / ( x**4 + x**2 + 0.9D+00 )

  return
end
function f3 ( x )

!*****************************************************************************80
!
!! f3() evaluates F3(X) = sqrt ( abs ( x + 1/2 ) ).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Charles Clenshaw, Alan Curtis,
!    A Method for Numerical Integration on an Automatic Computer,
!    Numerische Mathematik,
!    Volume 2, Number 1, December 1960, pages 197-205.
!
!  Input:
!
!    real ( kind = rk ) X, the argument.
!
!  Output:
!
!    real ( kind = rk ) F3, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) f3
  real ( kind = rk ) x

  f3 = sqrt ( abs ( x + 0.5D+00 ) )

  return
end
 
