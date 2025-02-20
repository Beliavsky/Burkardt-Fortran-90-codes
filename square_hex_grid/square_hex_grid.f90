subroutine box_print_2d ( box )

!*****************************************************************************80
!
!! box_print_2d() prints information about a coordinate box in 2D.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) BOX(2,2), the coordinates of the lower left
!    and upper right corners of the box.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) box(2,2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Coordinate box:'
  write ( *, '(a,f8.4,a,f8.4)' ) '  ', box(1,1), ' <= X <= ', box(1,2)
  write ( *, '(a,f8.4,a,f8.4)' ) '  ', box(2,1), ' <= Y <= ', box(2,2)

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IUNIT, the free unit number.
!
  implicit none

  integer i
  integer ios
  integer iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
subroutine hex_grid_01_approximate_h ( h_goal, nodes_per_layer, h )

!*****************************************************************************80
!
!! HEX_GRID_01_APPROXIMATE_H seeks a unit square hex grid with spacing H.
!
!  Discussion:
!
!    The parameter NODES_PER_LAYER controls the number of nodes and the
!    grid spacing, but in a somewhat obscure way.  This routine experiments
!    with various values until it is convinced it has the value 
!    of NODES_PER_LAYER that produces a grid spacing that is no
!    no greater than H.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) H_GOAL, the desired grid spacing.
!
!    Output, integer NODES_PER_LAYER, the number of nodes per layer
!    which produces a mesh with grid spacing H_GOAL or less.
!
!    Output, real ( kind = rk ) H, the actual grid spacing.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) h_goal
  real ( kind = rk ) h
  integer nodes_per_layer
  integer nodes_per_layer2

  if ( h_goal <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HEX_GRID_01_APPROXIMATE_H - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Illegal input value of H_GOAL = ', h_goal
    stop
  end if

  nodes_per_layer = 1 + int ( 1.0D+00 / h_goal )
!
!  Check whether roundoff means we could use one less node per layer.
!
  if ( 2 < nodes_per_layer ) then

    nodes_per_layer2 = nodes_per_layer - 1
    h = 1.0D+00 / real ( nodes_per_layer2 - 1, kind = rk )

    if ( h <= h_goal ) then
      nodes_per_layer = nodes_per_layer2
      return
    end if

  end if

  h = 1.0D+00 / real ( nodes_per_layer - 1, kind = rk )

  return
end
subroutine hex_grid_01_approximate_n ( n_goal, nodes_per_layer, n )

!*****************************************************************************80
!
!! HEX_GRID_01_APPROXIMATE_N seeks a unit square hex grid of about N nodes.
!
!  Discussion:
!
!    The parameter NODES_PER_LAYER controls the number of nodes, but
!    in a somewhat obscure way.  This routine experiments with various
!    values until it is convinced it has the value of NODES_PER_LAYER
!    that comes as close as possible to producing N nodes.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N_GOAL, the desired number of nodes.
!
!    Output, integer NODES_PER_LAYER, the number of nodes per layer
!    which produces a mesh with about N_GOAL nodes.
!
!    Output, integer N, the number of nodes in the mesh.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer n_goal
  integer n_high
  integer n_low
  integer nodes_per_layer
  integer nodes_per_layer_high
  integer nodes_per_layer_low

  if ( n_goal <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HEX_GRID_01_APPROXIMATE_N - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal input value of N_GOAL = ', n_goal
    stop
  end if

  nodes_per_layer_low = 0
  n_low = 0

  nodes_per_layer = nint ( sqrt ( real ( n_goal, kind = rk ) ) )

  nodes_per_layer_high = n_goal 
  n_high = n_goal**2

  do 

    call hex_grid_01_n ( nodes_per_layer, n )

    if ( n == n_goal ) then
      exit
    end if

    if ( n < n_goal ) then
      nodes_per_layer_low = nodes_per_layer
      n_low = n
    else
      nodes_per_layer_high = nodes_per_layer
      n_high = n
    end if

    if ( nodes_per_layer_low + 1 == nodes_per_layer_high ) then
      if ( n - n_low <= n_high - n ) then
        nodes_per_layer = nodes_per_layer_high
        n = n_high
      else
        nodes_per_layer = nodes_per_layer_low
        n = n_low
      end if
      exit
    end if

    nodes_per_layer = ( nodes_per_layer_low + nodes_per_layer_high ) / 2

  end do

  return
end
subroutine hex_grid_01_h ( nodes_per_layer, hx, hy )

!*****************************************************************************80
!
!! HEX_GRID_01_H computes the unit square hex grid spacings.
!
!  Discussion:
!
!    This routine determines the values of HX and HY from
!    the fundamental hexagonal grid parameter NODES_PER_LAYER.
!
!    A hexagonal grid is defined in the unit square [0,1] x [0,1].
!
!    All nodes of the grid lie on one of LAYERS horizontal lines.
!    The first of these lines is the X axis, and each successive
!    line is HY units higher.
!
!    On all the odd numbered lines, there are NODES_PER_LAYER points, 
!    equally spaced from 0 to 1, with a spacing of HX.  
!
!    On the even numbered lines, there are NODES_PER_LAYER-1 points,
!    whose values are the midpoints of successive intervals on
!    an odd numbered line.  (The grid is staggered).
!
!    HY = HX * sqrt ( 3 ) / 2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODES_PER_LAYER, the number of grid points
!    on the first horizontal layer of points.
!
!    Output, real ( kind = rk ) HX, the spacing between grid points
!    on a horizontal line.
!
!    Output, real ( kind = rk ) HY, the spacing between horizontal lines.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) hx
  real ( kind = rk ) hy
  integer nodes_per_layer

  if ( nodes_per_layer < 1 ) then

    hx = 0.0D+00
    hy = 0.0D+00

  else if ( nodes_per_layer == 1 ) then

    hx = 1.0D+00
    hy = 1.0D+00

  else

    hx = 1.0D+00 / real ( nodes_per_layer - 1, kind = rk )
    hy = hx * sqrt ( 3.0D+00 ) / 2.0D+00

  end if

  return
end
subroutine hex_grid_01_layers ( nodes_per_layer, layers )

!*****************************************************************************80
!
!! HEX_GRID_01_LAYERS computes the unit square hex grid column width.
!
!  Discussion:
!
!    This routine determines the value of LAYERS, the number of
!    layers, from the fundamental hexagonal grid parameter NODES_PER_LAYER.
!
!    A hexagonal grid is defined in the unit square [0,1] x [0,1].
!
!    All nodes of the grid lie on one of LAYERS horizontal lines.
!    The first of these lines is the X axis, and each successive
!    line is HY units higher.
!
!    On all the odd numbered lines, there are NODES_PER_LAYER points, 
!    equally spaced from 0 to 1, with a spacing of HX.  
!
!    On the even numbered lines, there are NODES_PER_LAYER-1 points,
!    whose values are the midpoints of successive intervals on
!    an odd numbered line.  (The grid is staggered).
!
!    HY = HX * sqrt ( 3 ) / 2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODES_PER_LAYER, the number of grid points 
!    on the first horizontal layer of points.
!
!    Output, integer LAYERS, the number of horizontal layers.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) hx
  real ( kind = rk ) hy
  integer layers
  integer nodes_per_layer

  if ( nodes_per_layer < 1 ) then

    layers = 0
    
  else if ( nodes_per_layer == 1 ) then

    layers = 1

  else

    hx = 1.0D+00 / real ( nodes_per_layer - 1, kind = rk )
    hy = sqrt ( 3.0D+00 ) * hx / 2.0D+00
    layers = 1 + int ( 1.0D+00 / hy )

  end if

  return
end
subroutine hex_grid_01_n ( nodes_per_layer, n )

!*****************************************************************************80
!
!! HEX_GRID_01_N computes the number of unit square hex grid points.
!
!  Discussion:
!
!    This routine determines the value of N, the number of
!    hex grid points, from the fundamental hexagonal grid 
!    parameter NODES_PER_LAYER.
!
!    A hexagonal grid is defined in the unit square [0,1] x [0,1].
!
!    All nodes of the grid lie on one of LAYERS horizontal lines.
!    The first of these lines is the X axis, and each successive
!    line is HY units higher.
!
!    On all the odd numbered lines, there are NODES_PER_LAYER points, 
!    equally spaced from 0 to 1, with a spacing of HX.  
!
!    On the even numbered lines, there are NODES_PER_LAYER-1 points,
!    whose values are the midpoints of successive intervals on
!    an odd numbered line.  (The grid is staggered).
!
!    HY = HX * sqrt ( 3 ) / 2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODES_PER_LAYER, the number of grid points 
!    on the first horizontal layer of points.
!
!    Output, integer N, the number of hex grid points.
!
  implicit none

  integer n
  integer layers
  integer nodes_per_layer

  if ( nodes_per_layer < 1 ) then

    n = 0

  else if ( nodes_per_layer == 1 ) then

    n = 1

  else

    call hex_grid_01_layers ( nodes_per_layer, layers )

    n = nodes_per_layer       * ( ( layers + 1 ) / 2 ) + &
      ( nodes_per_layer - 1 ) * ( ( layers     ) / 2 )

  end if

  return
end
subroutine hex_grid_01_points ( nodes_per_layer, layers, n, p )

!*****************************************************************************80
!
!! HEX_GRID_01_POINTS returns unit square hex grid points.
!
!  Discussion:
!
!    This routine determines the coordinates of the elements of
!    a hexagonal grid in the unit square.
!
!    A hexagonal grid is defined in the unit square [0,1] x [0,1].
!
!    All nodes of the grid lie on one of LAYERS horizontal lines.
!    The first of these lines is the X axis, and each successive
!    line is HY units higher.
!
!    On all the odd numbered lines, there are NODES_PER_LAYER points, 
!    equally spaced from 0 to 1, with a spacing of HX.  
!
!    On the even numbered lines, there are NODES_PER_LAYER-1 points,
!    whose values are the midpoints of successive intervals on
!    an odd numbered line.  (The grid is staggered).
!
!    HY = HX * sqrt ( 3 ) / 2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODES_PER_LAYER, the number of grid points on 
!    the first horizontal layer of points.
!
!    Input, integer LAYERS, the number of horizontal layers.
!
!    Input, integer N, the total number of hex grid points.
!
!    Output, real (kind = rk ) P(2,N), the coordinates of the
!    mesh points, listed one horizontal layer at a time.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) hx
  real ( kind = rk ) hy
  integer i
  integer j
  integer jmod
  integer k
  integer layers
  integer nodes_per_layer
  real ( kind = rk ) p(2,n)
  real ( kind = rk ) x
  real ( kind = rk ) y

  if ( nodes_per_layer < 1 ) then
    return
  end if

  if ( nodes_per_layer == 1 ) then
    p(1:2,1) = 0.5D+00
    return
  end if

  call hex_grid_01_h ( nodes_per_layer, hx, hy )

  k = 0

  do j = 1, layers

    y = hy * real ( j - 1, kind = rk )

    jmod = mod ( j, 2 )

    if ( jmod == 1 ) then

      do i = 1, nodes_per_layer
        x = real ( i - 1, kind = rk ) / real ( nodes_per_layer - 1, kind = rk )
        k = k + 1
        if ( k <= n ) then
          p(1,k) = x
          p(2,k) = y
        end if
      end do

    else

      do i = 1, nodes_per_layer-1
        x = real ( 2 * i - 1, kind = rk ) &
          / real ( 2 * nodes_per_layer - 2, kind = rk )
        k = k + 1
        if ( k <= n ) then
          p(1,k) = x
          p(2,k) = y
        end if
      end do

    end if

  end do

  return
end
subroutine hex_grid_approximate_h ( box, h_goal, nodes_per_layer, h )

!*****************************************************************************80
!
!! HEX_GRID_APPROXIMATE_H seeks a hex grid with spacing H.
!
!  Discussion:
!
!    The parameter NODES_PER_LAYER controls the number of nodes and the
!    grid spacing, but in a somewhat obscure way.  This routine experiments
!    with various values until it is convinced it has the value 
!    of NODES_PER_LAYER that produces a grid spacing that is no
!    no greater than H.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) BOX(2,2), the lower and upper corners
!    of the rectangular region.
!
!    Input, real ( kind = rk ) H_GOAL, the desired grid spacing.
!
!    Output, integer NODES_PER_LAYER, the number of nodes per layer
!    which produces a mesh with grid spacing H_GOAL or less.
!
!    Output, real ( kind = rk ) H, the actual grid spacing.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) box(2,2)
  real ( kind = rk ) h_goal
  real ( kind = rk ) h
  integer nodes_per_layer
  integer nodes_per_layer2

  if ( h_goal <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HEX_GRID_APPROXIMATE_H - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Illegal input value of H_GOAL = ', h_goal
    stop
  end if

  nodes_per_layer = 1 + int ( ( box(1,2) - box(1,1) ) / h_goal )
!
!  Check whether roundoff means we could use one less node per layer.
!
  if ( 2 < nodes_per_layer ) then

    nodes_per_layer2 = nodes_per_layer - 1
    h = ( box(1,2) - box(1,1) ) / real ( nodes_per_layer2 - 1, kind = rk )

    if ( h <= h_goal ) then
      nodes_per_layer = nodes_per_layer2
      return
    end if

  end if

  h = ( box(1,2) - box(1,1) ) / real ( nodes_per_layer - 1, kind = rk )

  return
end
subroutine hex_grid_approximate_n ( box, n_goal, nodes_per_layer, n )

!*****************************************************************************80
!
!! HEX_GRID_APPROXIMATE_N seeks a hex grid of about N nodes.
!
!  Discussion:
!
!    The parameter NODES_PER_LAYER controls the number of nodes, but
!    in a somewhat obscure way.  This routine experiments with various
!    values until it is convinced it has the value of NODES_PER_LAYER
!    that comes as close as possible to producing N nodes.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) BOX(2,2), the lower and upper corners
!    of the rectangular region.
!
!    Input, integer N_GOAL, the desired number of nodes.
!
!    Output, integer NODES_PER_LAYER, the number of nodes per layer
!    which produces a mesh with about N_GOAL nodes.
!
!    Output, integer N, the number of nodes in the mesh.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) box(2,2)
  integer n
  integer n_goal
  integer n_high
  integer n_low
  integer nodes_per_layer
  integer nodes_per_layer_high
  integer nodes_per_layer_low

  if ( n_goal <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HEX_GRID_APPROXIMATE_N - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal input value of N_GOAL = ', n_goal
    stop
  end if

  nodes_per_layer_low = 0
  n_low = 0

  nodes_per_layer = nint ( sqrt ( real ( n_goal, kind = rk ) ) )

  nodes_per_layer_high = n_goal  
  n_high = n_goal**2

  do 

    call hex_grid_n ( nodes_per_layer, box, n )

    if ( n == n_goal ) then
      exit
    end if

    if ( n < n_goal ) then
      nodes_per_layer_low = nodes_per_layer
      n_low = n
    else
      nodes_per_layer_high = nodes_per_layer
      n_high = n
    end if

    if ( nodes_per_layer_low + 1 == nodes_per_layer_high ) then
      if ( n - n_low <= n_high - n ) then
        nodes_per_layer = nodes_per_layer_high
        n = n_high
      else
        nodes_per_layer = nodes_per_layer_low
        n = n_low
      end if
      exit
    end if

    nodes_per_layer = ( nodes_per_layer_low + nodes_per_layer_high ) / 2

  end do

  return
end
subroutine hex_grid_h ( nodes_per_layer, box, hx, hy )

!*****************************************************************************80
!
!! HEX_GRID_H computes the coordinate box hex grid spacings.
!
!  Discussion:
!
!    This routine determines the values of HX and HY from
!    the fundamental hexagonal grid parameter NODES_PER_LAYER.
!
!    A hexagonal grid is defined in the coordinate box [A,B] x [C,D].
!
!    All nodes of the grid lie on one of LAYERS horizontal lines.
!    The first of these lines is from (A,C) to (B,C), and each 
!    successive line is HY units higher.
!
!    On all the odd numbered lines, there are NODES_PER_LAYER points, 
!    equally spaced from A to B, with a spacing of HX.  
!
!    On the even numbered lines, there are NODES_PER_LAYER-1 points,
!    whose values are the midpoints of successive intervals on
!    an odd numbered line.  (The grid is staggered).
!
!    HY = HX * sqrt ( 3 ) / 2.
!
!  Modified:
!
!    02 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODES_PER_LAYER, the number of grid points 
!    on the first horizontal layer of points.
!
!    Input, real ( kind = rk ) BOX(2,2), the values of A, B, C and D
!    that define the coordinate box.
!
!    Output, real ( kind = rk ) HX, the spacing between grid points
!    on a horizontal line.
!
!    Output, real ( kind = rk ) HY, the spacing between horizontal lines.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) box(2,2)
  real ( kind = rk ) hx
  real ( kind = rk ) hy
  integer nodes_per_layer

  if ( nodes_per_layer < 1 ) then

    hx = 0.0D+00
    hy = 0.0D+00

  else if ( nodes_per_layer == 1 ) then

    hx = box(1,2) - box(1,1)
    hy = box(2,2) - box(2,1)

  else

    hx = ( box(1,2) - box(1,1) ) / real ( nodes_per_layer - 1, kind = rk )
    hy = hx * sqrt ( 3.0D+00 ) / 2.0D+00

  end if

  return
end
subroutine hex_grid_layers ( nodes_per_layer, box, layers )

!*****************************************************************************80
!
!! HEX_GRID_LAYERS computes the coordinate box hex grid column width.
!
!  Discussion:
!
!    This routine determines the value of LAYERS, the number of
!    layers, from the fundamental hexagonal grid parameter NODES_PER_LAYER.
!
!    A hexagonal grid is defined in a coordinate box [A,B] x [C,D].
!
!    All nodes of the grid lie on one of LAYERS horizontal lines.
!    The first of these lines is from (A,C) to (B,C), and each 
!    successive line is HY units higher.
!
!    On all the odd numbered lines, there are NODES_PER_LAYER points, 
!    equally spaced from A to B, with a spacing of HX.  
!
!    On the even numbered lines, there are NODES_PER_LAYER-1 points,
!    whose values are the midpoints of successive intervals on
!    an odd numbered line.  (The grid is staggered).
!
!    HY = HX * sqrt ( 3 ) / 2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODES_PER_LAYER, the number of grid points 
!    on the first horizontal layer of points.
!
!    Input, real ( kind = rk ) BOX(2,2), the values of A, B, C and D
!    that define the coordinate box.
!
!    Output, integer LAYERS, the number of horizontal layers.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) box(2,2)
  real ( kind = rk ) hx
  real ( kind = rk ) hy
  integer layers
  integer nodes_per_layer

  if ( nodes_per_layer < 1 ) then

    layers = 0
    
  else if ( nodes_per_layer == 1 ) then

    layers = 1

  else

    hx = ( box(1,2) - box(1,1) ) / real ( nodes_per_layer - 1, kind = rk )
    hy = sqrt ( 3.0D+00 ) * hx / 2.0D+00
    layers = 1 + int ( ( box(2,2) - box(2,1) ) / hy )

  end if

  return
end
subroutine hex_grid_n ( nodes_per_layer, box, n )

!*****************************************************************************80
!
!! HEX_GRID_N computes the number of coordinate box hex grid points.
!
!  Discussion:
!
!    This routine determines the value of N, the number of
!    hex grid points, from the fundamental hexagonal grid 
!    parameter NODES_PER_LAYER.
!
!    A hexagonal grid is defined in the coordinate box [A,B] x [C,D].
!
!    All nodes of the grid lie on one of LAYERS horizontal lines.
!    The first of these lines is from (A,C) to (B,C), and each 
!    successive line is HY units higher.
!
!    On all the odd numbered lines, there are NODES_PER_LAYER points, 
!    equally spaced from A to B, with a spacing of HX.  
!
!    On the even numbered lines, there are NODES_PER_LAYER-1 points,
!    whose values are the midpoints of successive intervals on
!    an odd numbered line.  (The grid is staggered).
!
!    HY = HX * sqrt ( 3 ) / 2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODES_PER_LAYER, the number of grid points 
!    on the first horizontal layer of points.
!
!    Input, real ( kind = rk ) BOX(2,2), the values of A, B, C and D
!    that define the coordinate box.
!
!    Output, integer N, the number of hex grid points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) box(2,2)
  integer n
  integer layers
  integer nodes_per_layer

  if ( nodes_per_layer < 1 ) then

    n = 0

  else if ( nodes_per_layer == 1 ) then

    n = 1

  else

    call hex_grid_layers ( nodes_per_layer, box, layers )

    n = nodes_per_layer       * ( ( layers + 1 ) / 2 ) + &
      ( nodes_per_layer - 1 ) * ( ( layers     ) / 2 )

  end if

  return
end
subroutine hex_grid_plot ( m, n, box, p, header )

!*****************************************************************************80
!
!! HEX_GRID_PLOT creates GNUPLOT graphics files of a hex grid.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 August 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the spatial dimension, which should be 2.
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) BOX(M,2), two corners of the box enclosing
!    the grid.
!
!    Input, real ( kind = rk  P(M,N), the grid points.
!
!    Input, character ( len = * ) HEADER, a header for the file names.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  character ( len = 255 ) boundary_filename
  integer boundary_unit
  real ( kind = rk ) box(m,2)
  character ( len = 255 ) command_filename
  integer command_unit
  character ( len = 255 ) data_filename
  integer data_unit
  character ( len = * ) header
  integer i
  real ( kind = rk ) p(m,n)
  character ( len = 255 ) plot_filename
!
!  Create boundary data file.
!
  call get_unit ( boundary_unit )
  boundary_filename = trim ( header ) // '_boundary.txt'
  open ( unit = boundary_unit, file = boundary_filename, status = 'replace' )
  write ( boundary_unit, '(2x,g14.6,2x,g14.6)' ) box(1,1), box(2,1)
  write ( boundary_unit, '(2x,g14.6,2x,g14.6)' ) box(1,2), box(2,1)
  write ( boundary_unit, '(2x,g14.6,2x,g14.6)' ) box(1,2), box(2,2)
  write ( boundary_unit, '(2x,g14.6,2x,g14.6)' ) box(1,1), box(2,2)
  write ( boundary_unit, '(2x,g14.6,2x,g14.6)' ) box(1,1), box(2,1)
  close ( unit = boundary_unit )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Created boundary file "' // trim ( boundary_filename ) // '".'
!
!  Create point data file.
!
  call get_unit ( data_unit )
  data_filename = trim ( header ) // '_data.txt'
  open ( unit = data_unit, file = data_filename, status = 'replace' )
  do i = 1, n
    write ( data_unit, '(2x,g14.6,2x,g14.6)' ) p(1:2,i)
  end do
  close ( unit = data_unit )
  write ( *, '(a)' ) '  Created data file "' // trim ( data_filename ) // '".'
!
!  Create graphics command file.
!
  call get_unit ( command_unit )
  command_filename = trim ( header ) // '_commands.txt'
  open ( unit = command_unit, file = command_filename, status = 'replace' )
  write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  plot_filename = trim ( header ) // '.png'
  write ( command_unit, '(a)' ) 'set output "' // trim ( plot_filename ) // '"'
  write ( command_unit, '(a)' ) 'set xlabel "<--- X --->"'
  write ( command_unit, '(a)' ) 'set ylabel "<--- Y --->"'
  write ( command_unit, '(a)' ) 'set title "Disk Grid"'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set key off'
  write ( command_unit, '(a)' ) 'set timestamp'
  write ( command_unit, '(a)' ) 'set size ratio -1'
  write ( command_unit, '(a)' ) 'set style data lines'
  write ( command_unit, '(a)' ) 'plot "' // trim ( data_filename ) // &
    '" using 1:2 with points lt 3 pt 3,\'
  write ( command_unit, '(a)' ) '    "' // trim ( boundary_filename ) // &
    '" using 1:2 lw 3 linecolor rgb "black"'
  write ( command_unit, '(a)' ) 'quit'
  close ( unit = command_unit )

  write ( *, '(a)' ) &
    '  Created command file "' // trim ( command_filename ) // '".'

  return
end
subroutine hex_grid_points ( nodes_per_layer, layers, n, box, p )

!*****************************************************************************80
!
!! HEX_GRID_POINTS returns coordinate box hex grid points.
!
!  Discussion:
!
!    This routine determines the coordinates of the elements of
!    a hexagonal grid in the unit square.
!
!    A hexagonal grid is defined in the coordinate box [A,B] x [C,D].
!
!    All nodes of the grid lie on one of LAYERS horizontal lines.
!    The first of these lines is from (A,C) to (B,C), and each 
!    successive line is HY units higher.
!
!    On all the odd numbered lines, there are NODES_PER_LAYER points, 
!    equally spaced from A to B, with a spacing of HX.  
!
!    On the even numbered lines, there are NODES_PER_LAYER-1 points,
!    whose values are the midpoints of successive intervals on
!    an odd numbered line.  (The grid is staggered).
!
!    HY = HX * sqrt ( 3 ) / 2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODES_PER_LAYER, the number of grid points on 
!    the first horizontal layer of points.
!
!    Input, integer LAYERS, the number of horizontal layers.
!
!    Input, integer N, the total number of hex grid points.
!
!    Input, real ( kind = rk ) BOX(2,2), the values of A, B, C and D
!    that define the coordinate box.
!
!    Output, real (kind = rk ) P(2,N), the coordinates of the
!    mesh points, listed one horizontal layer at a time.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) box(2,2)
  real ( kind = rk ) hx
  real ( kind = rk ) hy
  integer i
  integer j
  integer jmod
  integer k
  integer layers
  integer nodes_per_layer
  real ( kind = rk ) p(2,n)
  real ( kind = rk ) x
  real ( kind = rk ) y

  if ( nodes_per_layer < 1 ) then
    return
  end if

  if ( nodes_per_layer == 1 ) then
    do i = 1, 2
      p(i,1) = ( box(i,1) + box(i,2) ) / 2.0D+00
    end do
    return
  end if

  call hex_grid_h ( nodes_per_layer, box, hx, hy )

  k = 0

  do j = 1, layers

    y = box(2,1) + hy * real ( j - 1, kind = rk )

    jmod = mod ( j, 2 )

    if ( jmod == 1 ) then

      do i = 1, nodes_per_layer
        x = box(1,1) + ( box(1,2) - box(1,1) ) &
          * real ( i - 1, kind = rk ) &
          / real ( nodes_per_layer - 1, kind = rk )
        k = k + 1
        if ( k <= n ) then
          p(1,k) = x
          p(2,k) = y
        end if
      end do

    else

      do i = 1, nodes_per_layer-1
        x = box(1,1) + ( box(1,2) - box(1,1) ) &
          * real ( 2 * i - 1, kind = rk ) &
          / real ( 2 * nodes_per_layer - 2, kind = rk )
        k = k + 1
        if ( k <= n ) then
          p(1,k) = x
          p(2,k) = y
        end if
      end do

    end if

  end do

  return
end
subroutine r8mat_write ( output_filename, m, n, table )

!*****************************************************************************80
!
!! R8MAT_WRITE writes an R8MAT file.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the output file name.
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) TABLE(M,N), the data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer j
  character ( len = * ) output_filename
  integer output_status
  integer output_unit
  character ( len = 30 ) string
  real ( kind = rk ) table(m,n)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_filename ) // '" on unit ', output_unit
    output_unit = -1
    stop 1
  end if
!
!  Create a format string.
!
!  For less precision in the output file, try:
!
!                                            '(', m, 'g', 14, '.', 6, ')'
!
  if ( 0 < m .and. 0 < n ) then

    write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 24, '.', 16, ')'
!
!  Write the data.
!
    do j = 1, n
      write ( output_unit, string ) table(1:m,j)
    end do

  end if
!
!  Close the file.
!
  close ( unit = output_unit )

  return
end
subroutine s_blank_delete ( s )

!*****************************************************************************80
!
!! S_BLANK_DELETE removes blanks from a string, left justifying the remainder.
!
!  Discussion:
!
!    All TAB characters are also removed.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
  implicit none

  character c
  integer get
  integer put
  integer nchar
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

  put = 0
  nchar = len_trim ( s )

  do get = 1, nchar

    c = s(get:get)

    if ( c /= ' ' .and. c /= TAB ) then
      put = put + 1
      s(put:put) = c
    end if

  end do

  s(put+1:nchar) = ' '

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
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
