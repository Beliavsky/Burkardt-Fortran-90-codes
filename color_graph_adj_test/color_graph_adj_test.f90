program main

!*****************************************************************************80
!
!! color_graph_adj_test() tests color_graph_adj().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 March 2023
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'color_graph_adj_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  color_graph_adj() implements graph algorithms.'

  call color_graph_adj_connect_random_test ( )
  call color_graph_adj_random_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'color_graph_adj_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine color_graph_adj_connect_random_test ( )

!*****************************************************************************80
!
!! color_graph_adj_connect_random_test() tests color_graph_adj_connect_random().
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

  integer, parameter :: nnode = 6

  integer adj(nnode,nnode)
  integer mcolor
  integer ncolor
  integer nedge
  integer result

  ncolor = 3
  nedge = 8

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'color_graph_adj_connect_random_test():'
  write ( *, '(a)' ) '  color_graph_adj_connect_random_test() returns a random '
  write ( *, '(a)' ) '  connected color graph;'
  write ( *, '(a)' ) ' '

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Random object is to have:'
  write ( *, '(a,i8)' ) '  Number of colors = ', ncolor
  write ( *, '(a,i8)' ) '  Number of nodes =  ', nnode
  write ( *, '(a,i8)' ) '  Number of edges =  ', nedge

  call color_graph_adj_connect_random ( nnode, nedge, ncolor, adj )

  call color_graph_adj_print ( adj, nnode, '  The graph:' )
!
!  Check connectedness.
!
  call graph_adj_is_edge_connected ( adj, nnode, result )

  write ( *, '(a)' ) ' '
  if ( result == 0 ) then
    write ( *, '(a)' ) '  The graph is NOT edgewise connected.'
  else
    write ( *, '(a)' ) '  The graph IS edgewise connected.'
  end if

  call graph_adj_is_node_connected ( adj, nnode, result )

  if ( result == 0 ) then
    write ( *, '(a)' ) '  The graph is NOT nodewise connected.'
  else
    write ( *, '(a)' ) '  The graph IS nodewise connected.'
  end if
!
!  Count the edges.
!
  call color_graph_adj_edge_count ( adj, nnode, nedge )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of edges is ', nedge
!
!  Count the colors.
!
  call color_graph_adj_color_count ( adj, nnode, mcolor, ncolor )
 
  write ( *, '(a,i8)' ) '  Number of colors is    ', ncolor
  write ( *, '(a,i8)' ) '  Maximum color index is ', mcolor

  return
end 
subroutine color_graph_adj_random_test ( )

!*****************************************************************************80
!
!! color_graph_adj_random_test() tests color_graph_adj_random().
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

  integer, parameter :: nnode = 6

  integer adj(nnode,nnode)
  integer mcolor
  integer ncolor
  integer nedge
 
  ncolor = 3
  nedge = 7

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'color_graph_adj_random_test():'
  write ( *, '(a)' ) '  color_graph_adj_random_test() returns a random color digraph.'
  write ( *, '(a)' ) ' '

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Random object is to have:'
  write ( *, '(a,i8)' ) '  Number of colors = ', ncolor
  write ( *, '(a,i8)' ) '  Number of nodes = ', nnode
  write ( *, '(a,i8)' ) '  Number of edges = ', nedge

  call color_graph_adj_random ( nnode, ncolor, nedge, adj )

  call color_graph_adj_print ( adj, nnode, '  The color graph:' )
!
!  Count the edges.
!
  call color_graph_adj_edge_count ( adj, nnode, nedge )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of edges is ', nedge
!
!  Count the colors.
!
  call color_graph_adj_color_count ( adj, nnode, mcolor, ncolor )
 
  write ( *, '(a,i8)' ) '  Number of colors is    ', ncolor
  write ( *, '(a,i8)' ) '  Maximum color index is ', mcolor

  return
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
