program main

!*****************************************************************************80
!
!! floyd_test() tests floyd().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 March 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'floyd_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  floyd() implements the Floyd algorithm for the'
  write ( *, '(a)' ) '  distance between all pairs of nodes in a graph.'

  call floyd_test01 ( )
  call floyd_test02 ( )
  call floyd_test03 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'floyd_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine floyd_test01 ( )

!*****************************************************************************80
!
!! floyd_test01 tests floyd().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 March 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 6

  real ( kind = rk ), dimension ( n, n ) :: A = reshape ( (/ &
     0.0D+00, -1.0D+00, -1.0D+00, -1.0D+00, -1.0D+00, -1.0D+00, &
     2.0D+00,  0.0D+00, -1.0D+00, -1.0D+00, -1.0D+00,  5.0D+00, &
     5.0D+00,  7.0D+00,  0.0D+00, -1.0D+00,  2.0D+00, -1.0D+00, &
    -1.0D+00,  1.0D+00,  4.0D+00,  0.0D+00, -1.0D+00,  2.0D+00, &
    -1.0D+00, -1.0D+00, -1.0D+00,  3.0D+00,  0.0D+00,  4.0D+00, &
    -1.0D+00,  8.0D+00, -1.0D+00, -1.0D+00,  3.0D+00,  0.0D+00  &
    /), (/ n, n /) )
  real ( kind = rk ) B(n,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'floyd_test01():'
  write ( *, '(a)' ) '  floyd() uses Floyd''s algorithm to find the'
  write ( *, '(a)' ) '  shortest distance between all pairs of nodes'
  write ( *, '(a)' ) '  in a directed graph, starting from the initial array'
  write ( *, '(a)' ) '  of direct node-to-node distances.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In the initial direct distance array, if'
  write ( *, '(a)' ) '    A(I,J) = -1,'
  write ( *, '(a)' ) '  this indicates there is NO directed link from'
  write ( *, '(a)' ) '  node I to node J.  In that case, the value of'
  write ( *, '(a)' ) '  of A(I,J) is essentially "infinity".'

  call r8mat_print ( n, n, A, '  Initial direct distance array:' )

  where ( A(1:n,1:n) == - 1.0D+00 )
    A = huge ( 1.0D+00 )
  end where

  call floyd ( n, A, B )

  where ( B(1:n,1:n) == huge ( 1.0D+00 ) )
    B = - 1.0D+00
  end where

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In the final shortest distance array, if'
  write ( *, '(a)' ) '    B(I,J) = -1,'
  write ( *, '(a)' ) '  this indicates there is NO directed path from'
  write ( *, '(a)' ) '  node I to node J.'

  call r8mat_print ( n, n, B, '  Final shortest distance array:' )

  return
end
subroutine floyd_test02 ( )

!*****************************************************************************80
!
!! floyd_test02 applies Floyd's algorithm to problems of increasing size.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 March 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  real ( kind = rk ) wtime

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'floyd_test02():'
  write ( *, '(a)' ) '  Test floyd() on matrices of increasing size.'
  write ( *, '(a)' ) '  The work is roughly N^3.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N   Time (seconds)  Time/N^3'
  write ( *, '(a)' ) ' '

  n = 1
  do while ( n <= 1024 )
    call floyd_test02_sub ( n, wtime )
    write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) &
      n, wtime, 1000000.0D+00 * wtime / real ( n**3, kind = rk )
    n = n * 2
  end do

  return
end
subroutine floyd_test02_sub ( n, wtime )

!*****************************************************************************80
!
!! floyd_test02_sub tests floyd().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 March 2022
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the size of the matrix.
!
!  Output:
!
!    real ( kind = rk ) WTIME, the CPU  time required.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) A(n,n)
  real ( kind = rk ) B(n,n)
  integer i
  integer j
  real ( kind = rk ) time1
  real ( kind = rk ) time2
  real ( kind = rk ) wtime

  call random_number ( A )
  do j = 1, n
    do i = 1, n
      if ( A(i,j) < 0.75D+00 ) then
        A(i,j) = huge ( 1.0D+00 )
      end if
    end do
  end do

  call cpu_time ( time1 )

  call floyd ( n, A, B )

  call cpu_time ( time2 )

  wtime = time2 - time1

  return
end
subroutine floyd_test03 ( )

!*****************************************************************************80
!
!! floyd_test03 uses Floyd's method for a triangulation.
!
!  Discussion:
!
!     8  41--42--43--44  45--46--47--48
!     |   | \ | \ | \ |   | \ | \ | \ |
!     7  33--34--35--36  37--38--39--40
!     |   | \ |                   | \ |
!     6  29--30                  31--32
!     |   | \ |                   | \ |
!     5  25--26                  27--28
!     |   | \ |                   | \ |
!     4  21--22                  23--24
!     |   | \ |                   | \ |
!     3  17--18                  19--20
!     |   | \ |                   | \ |
!     2   9--10--11--12--13--14--15--16
!     |   | \ | \ | \ | \ | \ | \ | \ |
!     1   1---2---3---4---5---6---7---8
!     |    
!     +---1---2---3---4---5---6---7---8
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 March 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: element_num = 46
  integer, parameter :: node_num = 48
  
  real ( kind = rk ) A(node_num,node_num)
  real ( kind = rk ) B(node_num,node_num)
  integer element
  integer, dimension ( 3, element_num ) :: element_node = &
    reshape ( (/ &
     1,  2,  9, &
     2, 10,  9, &
     2,  3, 10, &
     3, 11, 10, &
     3,  4, 11, &
     4, 12, 11, &
     4,  5, 12, &
     5, 13, 12, &
     5,  6, 13, &
     6, 14, 13, &
     6,  7, 14, &
     7, 15, 14, &
     7,  8, 15, &
     8, 16, 15, &
     9, 10, 17, &
    10, 18, 17, &
    15, 16, 19, &
    16, 20, 19, &
    17, 18, 21, &
    18, 22, 21, &
    19, 20, 23, &
    20, 24, 23, &
    21, 22, 25, &
    22, 26, 25, &
    23, 24, 27, &
    24, 28, 27, &
    25, 26, 29, &
    26, 30, 29, &
    27, 28, 31, &
    28, 32, 31, &
    29, 30, 33, &
    30, 34, 33, &
    31, 32, 39, &
    32, 40, 39, &
    33, 34, 41, &
    34, 42, 41, &
    34, 35, 42, &
    35, 43, 42, &
    35, 36, 43, &
    36, 44, 43, &
    37, 38, 45, &
    38, 46, 45, &
    38, 39, 46, &
    39, 47, 46, &
    39, 40, 47, &
    40, 48, 47 /), (/ 3, element_num /) )
  integer i
  integer n1
  integer n2
  real ( kind = rk ) r8vec_diff_norm
  real ( kind = rk ), dimension ( 2, node_num ) :: xy = reshape ( (/ &
    1.0, 1.0, &
    2.0, 1.0, &
    3.0, 1.0, &
    4.0, 1.0, &
    5.0, 1.0, &
    6.0, 1.0, &
    7.0, 1.0, &
    8.0, 1.0, &
    1.0, 2.0, &
    2.0, 2.0, &
    3.0, 2.0, &
    4.0, 2.0, &
    5.0, 2.0, &
    6.0, 2.0, &
    7.0, 2.0, &
    8.0, 2.0, &
    1.0, 3.0, & 
    2.0, 3.0, &
    7.0, 3.0, &
    8.0, 3.0, &
    1.0, 4.0, &
    2.0, 4.0, &
    7.0, 4.0, &
    8.0, 4.0, &
    1.0, 5.0, &
    2.0, 5.0, &
    7.0, 5.0, &
    8.0, 5.0, &
    1.0, 6.0, &
    2.0, 6.0, &
    7.0, 6.0, &
    8.0, 6.0, &
    1.0, 7.0, &
    2.0, 7.0, &
    3.0, 7.0, &
    4.0, 7.0, &
    5.0, 7.0, &
    6.0, 7.0, &
    7.0, 7.0, &
    8.0, 7.0, &
    1.0, 8.0, & 
    2.0, 8.0, &
    3.0, 8.0, &
    4.0, 8.0, &
    5.0, 8.0, &
    6.0, 8.0, &
    7.0, 8.0, &
    8.0, 8.0 /), (/ 2, node_num /) )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'floyd_test03():'
  write ( *, '(a)' ) '  floyd() analyzes node-to-node distances for'
  write ( *, '(a)' ) '  a triangulation.'
!
!  Must initialize distances to -1!
!
  A(1:node_num,1:node_num) = -1.0D+00
!
!  Diagonals are 0.
!
  do i = 1, node_num
    A(i,i) = 0.0D+00
  end do
!
!  Initialize A to the one-step distance.
!
  do element = 1, element_num
    n1 = element_node(3,element)
    do i = 1, 3
      n2 = element_node(i,element)
      A(n1,n2) = r8vec_diff_norm ( 2, xy(1:2,n1), xy(1:2,n2) )
      A(n2,n1) = A(n1,n2)
      n1 = n2
    end do
  end do
!
!  Reset -1 values to HUGE.
!
  where ( A(1:node_num,1:node_num) == - 1.0D+00 )
    A = huge ( 1.0D+00 )
  end where
!
!  Call floyd().
!
  call floyd ( node_num, A, B )
!
!  For printing, replace HUGE by -1.'
!
  where ( B == huge ( 1.0D+00 ) )
    B = - 1.0D+00
  end where

  call r8mat_print ( node_num, node_num, B, '  Distance matrix' )

  return
end

