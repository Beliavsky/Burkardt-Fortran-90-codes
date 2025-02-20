program main

!*****************************************************************************80
!
!! maze_test() tests maze().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    28 February 2023
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'maze_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  maze() implements some maze algorithms.'

  call maze_diameter_test ( )
  call maze_path_test ( )
  call maze_print_test ( )
  call maze_random_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'maze_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
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
!    15 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

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

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine maze_diameter_test ( )

!*****************************************************************************80
!
!! maze_diameter_test() tests maze_diameter().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 8
  integer, parameter :: n = 10

  integer bar(m,n+1)
  integer dad(m,n)
  integer degree(m,n)
  integer diam
  integer flat(m+1,n)
  integer i
  integer istart
  integer istop
  integer j
  integer jstart
  integer jstop
  integer path(m,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'maze_diameter_test():'
  write ( *, '(a)' ) '  maze_diameter(): find two far apart cells;'
!
!  Print out the cell numbers for the maze.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Cell numbers for the maze:'
  write ( *, '(a)' ) ' '
  do i = 1, m
    write ( *, '(20i3)' ) ( (j-1)*m+i, j = 1, n )
  end do
!
!  Get a random maze and print it.
!
  call maze_random ( m, n, bar, dad, flat )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A random maze:'
  write ( *, '(a,i8)' ) '    Number of rows =   ', m
  write ( *, '(a,i8)' ) '    Number of columns = ', n

  istart = 0
  jstart = 0

  istop = 0
  jstop = 0

  call maze_print ( bar, flat, m, n, istart, jstart, istop, jstop, &
    '  The maze:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Rooted tree representation:'
  write ( *, '(a)' ) '  (0 is the root.  All other cells print the'
  write ( *, '(a)' ) '  cell number of their parent on the tree.)'
  write ( *, '(a)' ) ' '
  do i = 1, m
    write ( *, '(20i3)' ) dad(i,1:n)
  end do
!
!  Get start and end points that are far apart and print the maze.
!
  call maze_diameter ( bar, degree, diam, flat, m, n, path, istart, jstart, &
    istop, jstop )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Random maze with far apart ends:'
  write ( *, '(a,i8)' ) '    Diameter = ', diam
  write ( *, '(a,2i8)' ) '    Starting cell = ', istart, jstart
  write ( *, '(a,2i8)' ) '    Stopping cell = ', istop, jstop

  call maze_print ( bar, flat, m, n, istart, jstart, istop, jstop, &
    '  The maze:' )

  return
end
subroutine maze_path_test ( )

!*****************************************************************************80
!
!! maze_path_test() tests maze_path().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 8
  integer, parameter :: n = 10

  integer bar(m,n+1)
  integer dad(m,n)
  integer degree(m,n)
  integer diam
  integer flat(m+1,n)
  integer i
  integer istart
  integer istop
  integer j
  integer jstart
  integer jstop
  integer path(m,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'maze_path_test():'
  write ( *, '(a)' ) '  maze_path() generates a path from one cell to'
  write ( *, '(a)' ) '  another, in a connected maze.'
!
!  Print out the cell numbers for the maze.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Cell numbers for the maze:'
  write ( *, '(a)' ) ' '
  do i = 1, m
    write ( *, '(20i3)' ) ( (j-1)*m+i, j = 1, n )
  end do
!
!  Get a random maze and print it.
!
  call maze_random ( m, n, bar, dad, flat )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A random maze:'
  write ( *, '(a,i8)' ) '    Number of rows =   ', m
  write ( *, '(a,i8)' ) '    Number of columns = ', n

  istart = 0
  jstart = 0

  istop = 0
  jstop = 0

  call maze_print ( bar, flat, m, n, istart, jstart, istop, jstop, &
    '  The maze:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Rooted tree representation:'
  write ( *, '(a)' ) '  (0 is the root.  All other cells print the'
  write ( *, '(a)' ) '  cell number of their parent on the tree.)'
  write ( *, '(a)' ) ' '
  do i = 1, m
    write ( *, '(20i3)' ) dad(i,1:n)
  end do
!
!  Get start and end points that are far apart and print the maze.
!
  call maze_diameter ( bar, degree, diam, flat, m, n, path, istart, jstart, &
    istop, jstop )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Random maze with far apart ends:'
  write ( *, '(a,i8)' ) '    Diameter = ', diam
  write ( *, '(a,2i8)' ) '    Starting cell = ', istart, jstart
  write ( *, '(a,2i8)' ) '    Stopping cell = ', istop, jstop

  call maze_print ( bar, flat, m, n, istart, jstart, istop, jstop, &
    '  The maze:' )
!
!  Find a path and print it.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Random maze with path from start to stop:'

  call maze_path ( bar, flat, m, n, istart, jstart, istop, jstop )

  call maze_print ( bar, flat, m, n, istart, jstart, istop, jstop, &
    '  The maze' )

  return
end
subroutine maze_print_test ( )

!*****************************************************************************80
!
!! maze_print_test() tests maze_print().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 2
  integer, parameter :: n = 3
 
  integer, parameter :: INDEF = -1
  integer, parameter :: WALL = 0
  integer, parameter :: OPEN = 1

  integer bar(m,n+1)
  integer flat(m+1,n)
  integer istart
  integer istop
  integer jstart
  integer jstop

  bar(1:m,1:n+1) = WALL
  flat(1:m+1,1:n) = WALL

  bar(1,2) = OPEN
  bar(1,4) = INDEF
  bar(2,3) = OPEN

  flat(1,3) = INDEF
  flat(2,1) = OPEN
  flat(2,2) = OPEN
  flat(2,3) = OPEN
  flat(3,1) = OPEN

  istart = 2
  jstart = 1

  istop = 1
  jstop = 3
!
!  Now mark the path.
!
  flat(2,1) = 2
  bar(1,2) = 2
  flat(2,2) = 2
  bar(2,3) = 2
  flat(2,3) = 2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'maze_print_test():'
  write ( *, '(a)' ) '  maze_print() prints a maze with path marked.'
  write ( *, '(a)' ) ' '

  call maze_print ( bar, flat, m, n, istart, jstart, istop, jstop, &
    '  The maze:' )

  return
end
subroutine maze_random_test ( )

!*****************************************************************************80
!
!! maze_random_test() tests maze_random().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 8
  integer, parameter :: n = 10

  integer bar(m,n+1)
  integer dad(m,n)
  integer flat(m+1,n)
  integer i
  integer istart
  integer istop
  integer j
  integer jstart
  integer jstop

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'maze_random_test():'
  write ( *, '(a)' ) '  maze_random() generate a random maze.'
!
!  Print out the cell numbers for the maze.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Cell numbers for the maze:'
  write ( *, '(a)' ) ' '
  do i = 1, m
    write ( *, '(20i3)' ) ( (j-1)*m+i, j = 1, n )
  end do
!
!  Get a random maze and print it.
!
  call maze_random ( m, n, bar, dad, flat )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A random maze:'
  write ( *, '(a,i8)' ) '    Number of rows =   ', m
  write ( *, '(a,i8)' ) '    Number of columns = ', n

  istart = 0
  jstart = 0

  istop = 0
  jstop = 0

  call maze_print ( bar, flat, m, n, istart, jstart, istop, jstop, &
    '  The maze:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Rooted tree representation:'
  write ( *, '(a)' ) '  (0 is the root.  All other cells print the'
  write ( *, '(a)' ) '  cell number of their parent on the tree.)'
  write ( *, '(a)' ) ' '
  do i = 1, m
    write ( *, '(20i3)' ) dad(i,1:n)
  end do

  return
end

