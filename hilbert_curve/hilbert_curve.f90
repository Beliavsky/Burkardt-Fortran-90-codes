subroutine d2xy ( m, d, x, y )

!*****************************************************************************80
!
!! d2xy() converts a 1D Hilbert coordinate to a 2D Cartesian coordinate.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 December 2015
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, the index of the Hilbert curve.
!    The number of cells is N=2^M.
!    0 < M.
!
!    integer D, the Hilbert coordinate of the cell.
!    0 <= D < N * N.
!
!  Output:
!
!    integer X, Y, the Cartesian coordinates of the cell.
!    0 <= X, Y < N.
!
  implicit none

  integer d
  integer m
  integer n
  integer rx
  integer ry
  integer s
  integer t
  integer x
  integer y

  n = 2 ** m

  x = 0
  y = 0
  t = d
  s = 1

  do while ( s < n )

    rx = mod ( t / 2, 2 )
    if ( rx == 0 ) then
      ry = mod ( t, 2 )
    else
      ry = mod ( ieor ( t, rx ), 2 )
    end if
    call rot ( s, x, y, rx, ry )
    x = x + s * rx
    y = y + s * ry
    t = t / 4

    s = s * 2

  end do

  return
end
subroutine rot ( n, x, y, rx, ry ) 

!*****************************************************************************80
!
!! rot() rotates and flips a quadrant appropriately.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 December 2015
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the length of a side of the square.  
!    N must be a power of 2.
!
!    integer X, Y, the coordinates of a point.
!
!    integer RX, RY, values of 0 or 1 which indicate
!    whether reflection or flips should be carried out.
!
!  Output:
!
!    integer X, Y, the coordinates of the modified point.
!
  implicit none

  integer n
  integer rx
  integer ry
  integer t
  integer x
  integer y

  if ( ry == 0 ) then
!
!  Reflect.
!
    if ( rx == 1 ) then
      x = n - 1 - x
      y = n - 1 - y
    end if
!
!  Flip.
!
    t = x
    x = y
    y = t

  end if

  return
end
subroutine xy2d ( m, x, y, d )

!*****************************************************************************80
!
!! xy2d() converts a 2D Cartesian coordinate to a 1D Hilbert coordinate.
!
!  Discussion:
!
!    It is assumed that a square has been divided into an NxN array of cells,
!    where N is a power of 2.
!
!    Cell (0,0) is in the lower left corner, and (N-1,N-1) in the upper 
!    right corner.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 January 2016
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, the index of the Hilbert curve.
!    The number of cells is N=2^M.
!    0 < M.
!
!    integer X, Y, the Cartesian coordinates of a cell.
!    0 <= X, Y < N.
!
!  Output:
!
!    integer D, the Hilbert coordinate of the cell.
!    0 <= D < N * N.
!
  implicit none

  integer d
  integer m
  integer n
  integer rx
  integer ry
  integer s
  integer x
  integer xcopy
  integer y
  integer ycopy

  xcopy = x
  ycopy = y

  d = 0
  n = 2 ** m

  s = n / 2

  do while ( 0 < s )

    if ( iand ( xcopy, s ) > 0 ) then
      rx = 1
    else
      rx = 0
    end if

    if ( iand ( ycopy, s ) > 0 ) then
      ry = 1
    else
      ry = 0
    end if

    d = d + s * s * ( ieor ( 3 * rx, ry ) )

    call rot ( s, xcopy, ycopy, rx, ry )

    s = s / 2

  end do

  return
end
