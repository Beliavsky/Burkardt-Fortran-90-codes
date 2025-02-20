program main

!*****************************************************************************80
!
!! digraph_arc_test() tests digraph_arc().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 March 2023
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'digraph_arc_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  digraph_arc() implements digraph algorithms.'

  call digraph_arc_degree_test ( )
  call digraph_arc_euler_circ_next_test ( )
  call digraph_arc_is_eulerian_test ( )
  call digraph_arc_to_digraph_adj_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'digraph_arc_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine digraph_arc_degree_test ( )

!*****************************************************************************80
!
!! digraph_arc_degree_test() tests digraph_arc_degree().
!
!  5--2--10--1--3--6
!         |  |  | /
!         8  |  9
!         |  |  
!         4--7  
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

  integer, parameter :: nedge = 11
  integer, parameter :: nnode = 10

  integer indegree(nnode)
  integer inode(nedge)
  integer jnode(nedge)
  integer outdegree(nnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'digraph_arc_degree_test():'
  write ( *, '(a)' ) '  digraph_arc_degree() computes the degree of the nodes;'

  inode = (/ 1, 1,  1, 2,  2, 3, 3, 4, 4, 6,  8 /)
  jnode = (/ 3, 7, 10, 5, 10, 6, 9, 7, 8, 9, 10 /)

  call digraph_arc_print ( nedge, inode, jnode, '  The graph:' )

  call digraph_arc_degree ( nnode, nedge, inode, jnode, indegree, outdegree )

  call i4vec2_print ( nnode, indegree, outdegree, '  Node, Indegree, Outdegree' )

  return
end 
subroutine digraph_arc_euler_circ_next_test ( )

!*****************************************************************************80
!
!! digraph_arc_euler_circ_next_test() tests digraph_arc_euler_circ_next().
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

  integer, parameter :: maxstack = 130
  integer, parameter :: nedge = 10
  integer, parameter :: nnode = 5

  integer circuit(nedge)
  integer i
  integer, dimension ( nedge ) :: inode = (/ 1, 3, 1, 5, 2, 4, 2, 4, 3, 5 /)
  integer, dimension ( nedge ) :: jnode = (/ 2, 1, 4, 1, 3, 2, 5, 3, 5, 4 /)
  logical more
  integer ncan(nedge)
  integer stack(maxstack)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'digraph_arc_euler_circ_next_test():'
  write ( *, '(a)' ) '  digraph_arc_euler_circ_next() finds the next'
  write ( *, '(a)' ) '  Euler circuit of a graph.'
  write ( *, '(a)' ) ' '

  call digraph_arc_print ( nedge, inode, jnode, '  The digraph:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Circuits:'
  write ( *, '(a)' ) ' '
  i = 0

  more = .false.

  do

    call digraph_arc_euler_circ_next ( nedge, inode, jnode, circuit, stack, &
      maxstack, ncan, more )

    if ( .not. more ) then
      exit
    end if

    i = i + 1
    write ( *, '(i3,2x,20i3)' ) i, circuit(1:nedge)

  end do

  return
end
subroutine digraph_arc_is_eulerian_test ( )

!*****************************************************************************80
!
!! digraph_arc_is_eulerian_test() tests digraph_arc_is_eulerian().
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

  integer, parameter :: maxstack = 130
  integer, parameter :: nedge = 10
  integer, parameter :: nnode = 5

  integer indegree(nnode)
  integer, dimension ( nedge ) :: inode = (/ 1, 3, 1, 5, 2, 4, 2, 4, 3, 5 /)
  integer, dimension ( nedge ) :: jnode = (/ 2, 1, 4, 1, 3, 2, 5, 3, 5, 4 /)
  integer outdegree(nnode)
  integer result

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'digraph_arc_is_eulerian_test()'
  write ( *, '(a)' ) '  digraph_arc_is_eulerian() checks if a graph'
  write ( *, '(a)' ) '  has (at least one) Euler circuit.'
  write ( *, '(a)' ) ' '

  call digraph_arc_print ( nedge, inode, jnode, '  The digraph:' )

  call digraph_arc_is_eulerian ( nnode, nedge, inode, jnode, indegree, &
    outdegree, result )

  if ( result == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The digraph is NOT eulerian.'
    return
  else if ( result == 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The digraph has an eulerian path,'
    write ( *, '(a)' ) '  but not an eulerian circuit.'
  else if ( result == 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The digraph has an eulerian circuit.'
  end if

  return
end
subroutine digraph_arc_to_digraph_adj_test ( )

!*****************************************************************************80
!
!! digraph_arc_to_digraph_adj_test() tests digraph_arc_to_digraph_adj().
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

  integer, parameter :: nedge = 16
  integer, parameter :: nnode = 20

  integer adj(nnode,nnode)
  integer inode(nedge)
  integer jnode(nedge)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'digraph_arc_to_digraph_adj_test():'
  write ( *, '(a)' ) '  digraph_arc_to_digraph_adj() converts an arclist'
  write ( *, '(a)' ) '  digraph to an adjacency digraph.'
  write ( *, '(a)' ) ' '

  call digraph_arc_example_cycler ( inode, jnode )

  call digraph_arc_print ( nedge, inode, jnode, '  The graph:' )

  call digraph_arc_to_digraph_adj ( nedge, inode, jnode, adj, nnode )

  call digraph_adj_print ( adj, nnode, '  The digraph:' )

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
