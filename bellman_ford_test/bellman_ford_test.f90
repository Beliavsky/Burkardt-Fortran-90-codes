program main

!*****************************************************************************80
!
!! bellman_ford_test() tests bellman_ford().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'bellman_ford_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test bellman_ford().'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'bellman_ford_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! test01() runs a simple test.
!
!  Discussion:
!
!    The correct distances are { 0, -6, -2, 3, 0, 0 }.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: e_num = 10
  integer, parameter :: v_num = 6

  integer, dimension (2,e_num) :: e = reshape ( (/ &
    2, 1, &
    5, 2, &
    2, 3, &
    3, 5, &
    5, 1, &
    3, 6, &
    6, 1, &
    4, 3, &
    6, 4, &
    4, 1 /), (/ 2, e_num /) )
  real ( kind = rk ), dimension ( e_num ) :: e_weight = (/ &
    -3.0, &
     6.0, &
    -4.0, &
    -1.0, &
     4.0, &
    -2.0, &
     2.0, &
     8.0, &
    -3.0, &
     3.0 /)
  integer predecessor(v_num)
  integer source
  real ( kind = rk ) v_weight(v_num)
  
  source = 1

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'test01()::'
  write ( *, '(a)' ) '  Bellman-Ford shortest path algorithm.'

  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Number of vertices = ', v_num
  write ( *, '(a,i4)' ) '  Number of edges = ', e_num
  write ( *, '(a,i4)' ) '  The reference vertex is ', source

  call i4mat_transpose_print ( 2, e_num, e, '  The edge array:' )
  call r8vec_print ( e_num, e_weight, '  The edge weights:' )

  call bellman_ford ( v_num, e_num, source, e, e_weight, v_weight, &
    predecessor )

  call r8vec_print ( v_num, v_weight, "  The shortest distances:" )

  call i4vec_print ( v_num, predecessor, &
    '  The vertex predecessor parents for the shortest paths:' )

  return
end
