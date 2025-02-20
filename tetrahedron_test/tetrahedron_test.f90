program main

!*****************************************************************************80
!
!! tetrahedron_test() tests tetrahedron().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 May 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'tetrahedron_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test tetrahedron().'

  call tetrahedron_barycentric_test ( )
  call tetrahedron_centroid_test ( )
  call tetrahedron_contains_point_test ( )
  call tetrahedron_circumsphere_test ( )
  call tetrahedron_edge_length_test ( )
  call tetrahedron_insphere_test ( )
  call tetrahedron_lattice_layer_point_next_test ( )
  call tetrahedron_lattice_point_next_test ( )
  call tetrahedron_quality1_test ( )
  call tetrahedron_quality2_test ( )
  call tetrahedron_quality3_test ( )
  call tetrahedron_quality4_test ( )
  call tetrahedron_rhombic_shape_test ( )
  call tetrahedron_sample_test ( )
  call tetrahedron_shape_test ( )
  call tetrahedron_solid_angles_test ( )
  call tetrahedron_volume_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'tetrahedron_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine tetrahedron_barycentric_test ( )

!*****************************************************************************80
!
!! tetrahedron_barycentric_test() tests tetrahedron_barycentric().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer, parameter :: test_num = 10

  real ( kind = rk ) p(dim_num,test_num)
  real ( kind = rk ), dimension(dim_num,4) :: t = reshape ( (/ &
     1.0D+00, 4.0D+00, 3.0D+00, &
     2.0D+00, 4.0D+00, 3.0D+00, &
     1.0D+00, 6.0D+00, 3.0D+00, &
     1.0D+00, 4.0D+00, 4.0D+00 /), (/ dim_num, 4 /) )
  integer test
  real ( kind = rk ) xsi(dim_num+1)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'tetrahedron_barycentric_test():'
  write ( *, '(a)' ) '  tetrahedron_barycentric() converts Cartesian to'
  write ( *, '(a)' ) '  barycentric coordinates.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  We are computing the barycentric coordinates just to'
  write ( *, '(a)' ) '  verify that the points are inside the tetrahedron.'

  call r8mat_transpose_print ( dim_num, 4, t, '  Tetrahedron vertices' )

  call tetrahedron_sample ( t, test_num, p )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '      P                           Barycentric:'
  write ( *, '(a)' ) ''

  do test = 1, test_num
    call tetrahedron_barycentric ( t, p(1:3,test), xsi )
    write ( *, '(2x,3f8.4,4x,4f8.4)' ) p(1:dim_num,test), xsi(1:dim_num+1)
  end do

  return
end
subroutine tetrahedron_centroid_test ( )

!*****************************************************************************80
!
!! tetrahedron_centroid_test() tests tetrahedron_centroid();
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) centroid(dim_num)
  real ( kind = rk ), dimension (dim_num,4) :: tetra = reshape ( (/&
     0.000000D+00,  0.942809D+00, -0.333333D+00, &
    -0.816496D+00, -0.816496D+00, -0.333333D+00, &
     0.816496D+00, -0.816496D+00, -0.333333D+00, &
     0.000000D+00,  0.000000D+00,  1.000000D+00 /), (/ dim_num, 4 /) )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'tetrahedron_centroid_test():'
  write ( *, '(a)' ) '  tetrahedron_centroid() computes the centroid'
  write ( *, '(a)' ) '  of a tetrahedron.'

  call r8mat_transpose_print ( dim_num, 4, tetra, '  Tetrahedron vertices:' )

  call tetrahedron_centroid ( tetra, centroid )

  call r8vec_print ( dim_num, centroid, '  Centroid:' )

  return
end
subroutine tetrahedron_contains_point_test ( )

!*****************************************************************************80
!
!! tetrahedron_contains_point_test() tests tetrahedron_contains_point();
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer, parameter :: test_num = 3

  real ( kind = rk ) c(4)
  real ( kind = rk ), dimension(4,test_num) :: c_test = reshape ( (/ &
     0.0D+00, 0.1D+00,  0.2D+00, 0.7D+00, &
    -1.3D+00, 2.0D+00,  0.2D+00, 0.1D+00, &
     0.8D+00, 0.6D+00, -0.5D+00, 0.1D+00 /), (/ 4, test_num /) )
  logical inside
  real ( kind = rk ) p(dim_num)
  integer test
  real ( kind = rk ), dimension(dim_num,4) :: tetra = reshape ( (/&
     0.000000D+00,  0.942809D+00, -0.333333D+00, &
    -0.816496D+00, -0.816496D+00, -0.333333D+00, &
     0.816496D+00, -0.816496D+00, -0.333333D+00, &
     0.000000D+00,  0.000000D+00,  1.000000D+00 /), (/ dim_num, 4 /) )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'tetrahedron_contains_point_test():'
  write ( *, '(a)' ) '  tetrahedron_contains_point() finds if a point '
  write ( *, '(a)' ) '  is inside of a tetrahedron.'

  call r8mat_transpose_print ( dim_num, 4, tetra, '  Tetrahedron vertices:' )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  P     Inside_Tetra?'
  write ( *, '(a)' ) ''

  do test = 1, test_num

    c(1:4) = c_test(1:4,test)

    p(1:dim_num) = matmul ( tetra(1:dim_num,1:4), c(1:4) )

    call tetrahedron_contains_point ( tetra, p, inside )

    write ( *, '(2x,3g14.6,2x,l1)' ) p(1:dim_num), inside

  end do

  return
end
subroutine tetrahedron_circumsphere_test ( )

!*****************************************************************************80
!
!! tetrahedron_circumsphere_test() tests tetrahedron_circumsphere();
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) r
  real ( kind = rk ), dimension(dim_num,4) :: tetra = reshape ( (/&
     0.577350269189626D+00,  0.0D+00, 0.0D+00, &
    -0.288675134594813D+00,  0.5D+00, 0.0D+00, &
    -0.288675134594813D+00, -0.5D+00, 0.0D+00, &
     0.0D+00,                0.0D+00, 0.816496580927726D+00 /), &
    (/ dim_num, 4 /) )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'tetrahedron_circumsphere_test():'
  write ( *, '(a)' ) '  tetrahedron_circumsphere() computes the circumsphere'
  write ( *, '(a)' ) '  of a tetrahedron.'

  call r8mat_transpose_print ( dim_num, 4, tetra, '  Tetrahedron vertices:' )

  call tetrahedron_circumsphere ( tetra, r, pc )

  call r8vec_print ( dim_num, pc, '  Circumsphere center:' )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Circumsphere radius is ', r
 
  return
end
subroutine tetrahedron_edge_length_test ( )

!*****************************************************************************80
!
!! tetrahedron_edge_length_test() tests tetrahedron_edge_length();
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) edge_length(6)
  real ( kind = rk ), dimension(dim_num,4) :: tetra = reshape ( (/&
     0.577350269189626D+00,  0.0D+00, 0.0D+00, &
    -0.288675134594813D+00,  0.5D+00, 0.0D+00, &
    -0.288675134594813D+00, -0.5D+00, 0.0D+00, &
     0.0D+00,                0.0D+00, 0.816496580927726D+00 /), &
     (/ dim_num, 4 /) )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'tetrahedron_edge_length_test():'
  write ( *, '(a)' ) '  tetrahedron_edge_length() computes the edge lengths'
  write ( *, '(a)' ) '  of a tetrahedron.'

  call r8mat_transpose_print ( dim_num, 4, tetra, '  Tetrahedron vertices:' )

  call tetrahedron_edge_length ( tetra, edge_length )

  call r8vec_print ( 6, edge_length, '  Edge lengths:' )
 
  return
end
subroutine tetrahedron_insphere_test ( )

!*****************************************************************************80
!
!! tetrahedron_insphere_test() tests tetrahedron_insphere();
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    08 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) r
  real ( kind = rk ), dimension(dim_num,4) :: tetra = reshape ( (/&
     0.577350269189626D+00,  0.0D+00, 0.0D+00, &
    -0.288675134594813D+00,  0.5D+00, 0.0D+00, &
    -0.288675134594813D+00, -0.5D+00, 0.0D+00, &
     0.0D+00,                0.0D+00, 0.816496580927726D+00 /), &
        (/ dim_num, 4 /) )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'tetrahedron_insphere_test():'
  write ( *, '(a)' ) '  tetrahedron_insphere() computes the insphere'
  write ( *, '(a)' ) '  of a tetrahedron.'

  call r8mat_transpose_print ( dim_num, 4, tetra, '  Tetrahedron vertices:' )

  call tetrahedron_insphere ( tetra, r, pc )

  call r8vec_print ( dim_num, pc, '  Insphere center:' )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Insphere radius is ', r
 
  return
end
subroutine tetrahedron_lattice_layer_point_next_test ( )

!*****************************************************************************80
!
!! tetrahedron_lattice_layer_point_next_test() tests tetrahedron_lattice_layer_point_next().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 July 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 3

  integer c(n+1)
  integer i
  integer layer
  logical more
  integer v(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'tetrahedron_lattice_layer_point_next_test():'
  write ( *, '(a)' ) '  tetrahedron_lattice_layer_point_next() returns the next'
  write ( *, '(a)' ) '  point in a tetrahedron lattice layer defined by:'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '    C(4) - 1 < X(1)/C(1) + X(2)/C(2) +X(3)/C(3) <= C(4).'

  c(1) = 2
  c(2) = 3
  c(3) = 4
  v(1:n) = 0

  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  N = ', n
  write ( *, '(a)', ADVANCE = 'NO' ) '  C =       '
  do i = 1, n
    write ( *, '(2x,i4)', ADVANCE = 'NO' ) c(i)
  end do
  write ( *, '(a)', ADVANCE = 'YES' ) 

  do layer = 0, 2

    write ( *, '(a)' ) ''
    write ( *, '(a,i4)' ) '  Layer ', layer
    write ( *, '(a)' ) ''

    c(4) = layer
    more = .false.
    i = 0

    do
      call tetrahedron_lattice_layer_point_next ( c, v, more )
      if ( .not. more ) then
        write ( *, '(a)' ) '  No more.'
        exit
      end if
      i = i + 1
      write ( *, '(2x,i4,6x,10(2x,i4))' ) i, v(1:n)

    end do

  end do

  return
end
subroutine tetrahedron_lattice_point_next_test ( )

!*****************************************************************************80
!
!! tetrahedron_lattice_point_next_test() tests tetrahedron_lattice_point_next().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 July 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 3

  integer c(n+1)
  integer i
  logical more
  integer v(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'tetrahedron_lattice_point_next():'
  write ( *, '(a)' ) '  tetrahedron_lattice_point_next() returns the next lattice'
  write ( *, '(a)' ) '  point in a tetrahedron defined by:'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '    0 <= X(1)/C(1) + X(2)/C(2) + X(3)/C(3) <= C(4).'

  do i = 1, n + 1
    c(i) = n + 2 - i
  end do
  v(1:n) = 0
  more = .false.

  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  N = ', n
  write ( *, '(a)', ADVANCE = 'NO' ) '  C =       '
  do i = 1, n + 1
    write ( *, '(2x,i4)', ADVANCE = 'NO' ) c(i)
  end do
  write ( *, '(a)', ADVANCE = 'YES' ) 
  write ( *, '(a)' ) ''

  i = 0

  do
    call tetrahedron_lattice_point_next ( c, v, more )
    if ( .not. more ) then
      write ( *, '(a)' ) '  No more.'
      exit
    end if
    i = i + 1
    write ( *, '(2x,i4,6x,10(2x,i4))' ) i, v(1:n)
  end do

  return
end
subroutine tetrahedron_quality1_test ( )

!*****************************************************************************80
!
!! tetrahedron_quality1_test tests tetrahedron_quality1();
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer, parameter :: test_num = 2

  real ( kind = rk ) quality
  real ( kind = rk ), dimension(dim_num,4) :: tetra
  real ( kind = rk ), dimension(dim_num,4,test_num) :: tetra_test = reshape ( (/&
     0.577350269189626D+00,  0.0D+00, 0.0D+00, &
    -0.288675134594813D+00,  0.5D+00, 0.0D+00, &
    -0.288675134594813D+00, -0.5D+00, 0.0D+00, &
     0.0D+00,                0.0D+00, 0.816496580927726D+00, &
     0.577350269189626D+00,  0.0D+00, 0.0D+00, &
    -0.288675134594813D+00,  0.5D+00, 0.0D+00, &
    -0.288675134594813D+00, -0.5D+00, 0.0D+00, &
     0.0D+00,                0.0D+00, 0.408248290463863D+00 /), &
        (/ dim_num, 4, test_num /) )
  integer test

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'tetrahedron_quality1_test():'
  write ( *, '(a)' ) '  tetrahedron_quality1() computes quality measure #1'
  write ( *, '(a)' ) '  of a tetrahedron.'

  do test = 1, test_num

    tetra(1:dim_num,1:4) = tetra_test(1:dim_num,1:4,test)

    call r8mat_transpose_print ( dim_num, 4, tetra, '  Tetrahedron vertices:' )

    call tetrahedron_quality1 ( tetra, quality )

    write ( *, '(a,g14.6)' ) '  Tetrahedron quality is ', quality

  end do
 
  return
end
subroutine tetrahedron_quality2_test ( )

!*****************************************************************************80
!
!! tetrahedron_quality2_test() tests tetrahedron_quality2();
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer, parameter :: test_num = 2

  real ( kind = rk ) quality2
  real ( kind = rk ), dimension(dim_num,4) :: tetra
  real ( kind = rk ), dimension(dim_num,4,test_num) :: tetra_test = reshape ( (/&
     0.577350269189626D+00,  0.0D+00, 0.0D+00, &
    -0.288675134594813D+00,  0.5D+00, 0.0D+00, &
    -0.288675134594813D+00, -0.5D+00, 0.0D+00, &
     0.0D+00,                0.0D+00, 0.816496580927726D+00, &
     0.577350269189626D+00,  0.0D+00, 0.0D+00, &
    -0.288675134594813D+00,  0.5D+00, 0.0D+00, &
    -0.288675134594813D+00, -0.5D+00, 0.0D+00, &
     0.0D+00,                0.0D+00, 0.408248290463863D+00 /), &
        (/ dim_num, 4, test_num /) )
  integer test

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'tetrahedron_quality2_test():'
  write ( *, '(a)' ) '  tetrahedron_quality2() computes quality measure #2'
  write ( *, '(a)' ) '  of a tetrahedron.'

  do test = 1, test_num

    tetra(1:dim_num,1:4) = tetra_test(1:dim_num,1:4,test)

    call r8mat_transpose_print ( dim_num, 4, tetra, '  Tetrahedron vertices:' )

    call tetrahedron_quality2 ( tetra, quality2 )

    write ( *, '(a,g14.6)' ) '  Tetrahedron quality is ', quality2

  end do

  return
end
subroutine tetrahedron_quality3_test ( )

!*****************************************************************************80
!
!! tetrahedron_quality3_test tests tetrahedron_quality3();
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer, parameter :: test_num = 2

  real ( kind = rk ) quality3
  real ( kind = rk ), dimension(dim_num,4) :: tetra
  real ( kind = rk ), dimension(dim_num,4,test_num) :: tetra_test = reshape ( (/&
     0.577350269189626D+00,  0.0D+00, 0.0D+00, &
    -0.288675134594813D+00,  0.5D+00, 0.0D+00, &
    -0.288675134594813D+00, -0.5D+00, 0.0D+00, &
     0.0D+00,                0.0D+00, 0.816496580927726D+00, &
     0.577350269189626D+00,  0.0D+00, 0.0D+00, &
    -0.288675134594813D+00,  0.5D+00, 0.0D+00, &
    -0.288675134594813D+00, -0.5D+00, 0.0D+00, &
     0.0D+00,                0.0D+00, 0.408248290463863D+00 /), &
        (/ dim_num, 4, test_num /) )
  integer test

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'tetrahedron_quality3_test():'
  write ( *, '(a)' ) '  tetrahedron_quality3() computes quality measure #3'
  write ( *, '(a)' ) '  of a tetrahedron.'

  do test = 1, test_num

    tetra(1:dim_num,1:4) = tetra_test(1:dim_num,1:4,test)

    call r8mat_transpose_print ( dim_num, 4, tetra, '  Tetrahedron vertices:' )

    call tetrahedron_quality3 ( tetra, quality3 )

    write ( *, '(a,g14.6)' ) '  Tetrahedron quality is ', quality3

  end do

  return
end
subroutine tetrahedron_quality4_test ( )

!*****************************************************************************80
!
!! tetrahedron_quality4_test() tests tetrahedron_quality4();
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer, parameter :: test_num = 2

  real ( kind = rk ) quality4
  real ( kind = rk ), dimension(dim_num,4) :: tetra
  real ( kind = rk ), dimension(dim_num,4,test_num) :: tetra_test = reshape ( (/&
     0.577350269189626D+00,  0.0D+00, 0.0D+00, &
    -0.288675134594813D+00,  0.5D+00, 0.0D+00, &
    -0.288675134594813D+00, -0.5D+00, 0.0D+00, &
     0.0D+00,                0.0D+00, 0.816496580927726D+00, &
     0.577350269189626D+00,  0.0D+00, 0.0D+00, &
    -0.288675134594813D+00,  0.5D+00, 0.0D+00, &
    -0.288675134594813D+00, -0.5D+00, 0.0D+00, &
     0.0D+00,                0.0D+00, 0.408248290463863D+00 /), &
        (/ dim_num, 4, test_num /) )
  integer test

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'tetrahedron_quality4_test():'
  write ( *, '(a)' ) '  tetrahedron_quality4() computes quality measure #4'
  write ( *, '(a)' ) '  of a tetrahedron.'

  do test = 1, test_num

    tetra(1:dim_num,1:4) = tetra_test(1:dim_num,1:4,test)

    call r8mat_transpose_print ( dim_num, 4, tetra, '  Tetrahedron vertices:' )

    call tetrahedron_quality4 ( tetra, quality4 )

    write ( *, '(a,g14.6)' ) '  Tetrahedron quality is ', quality4

  end do

  return
end
subroutine tetrahedron_rhombic_shape_test ( )

!*****************************************************************************80
!
!! tetrahedron_rhombic_shape_test() tests tetrahedron_rhombic_shape().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer edge_num
  integer face_num
  integer, allocatable, dimension ( : ) ::  face_order
  integer face_order_max
  integer, allocatable, dimension ( :, : ) ::  face_point
  integer point_num
  real ( kind = rk ), allocatable, dimension ( :, : ) :: point_coord

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'tetrahedron_rhombic_shape_test():'
  write ( *, '(a)' ) '  tetrahedron_rhombic_size() returns dimension information;'
  write ( *, '(a)' ) '  tetrahedron_rhombic_shape() returns face and order information.'
  write ( *, '(a)' ) '  shape_print() prints this information.'
!
!  Get the data sizes.
!
  call tetrahedron_rhombic_size ( point_num, edge_num, face_num, &
    face_order_max )

  write ( *, '(a)' ) ''
  write ( *, '(a,i8)' ) '  Number of vertices = ', point_num
  write ( *, '(a,i8)' ) '  Number of edges =    ', edge_num
  write ( *, '(a,i8)' ) '  Number of faces =    ', face_num
  write ( *, '(a,i8)' ) '  Maximum face order = ', face_order_max
!
!  Make room for the data.
!
  allocate ( face_order(1:face_num) )
  allocate ( face_point(1:face_order_max,1:face_num) )
  allocate ( point_coord(1:3,1:point_num) )
!
!  Get the data.
!
  call tetrahedron_rhombic_shape ( point_num, face_num, face_order_max, &
    point_coord, face_order, face_point )
!
!  Print the data.
!
  call shape_print ( point_num, face_num, face_order_max, &
    point_coord, face_order, face_point )

  deallocate ( face_order )
  deallocate ( face_point )
  deallocate ( point_coord )

  return
end
subroutine tetrahedron_sample_test ( )

!*****************************************************************************80
!
!! tetrahedron_sample_test() tests tetrahedron_sample().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer, parameter :: test_num = 10

  real ( kind = rk ) p(dim_num,test_num)
  real ( kind = rk ), dimension(dim_num,4) :: t = reshape ( (/ &
     1.0D+00, 4.0D+00, 3.0D+00, &
     2.0D+00, 4.0D+00, 3.0D+00, &
     1.0D+00, 6.0D+00, 3.0D+00, &
     1.0D+00, 4.0D+00, 4.0D+00 /), (/ dim_num, 4 /) )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'tetrahedron_sample_test():'
  write ( *, '(a)' ) '  tetrahedron_sample() samples a tetrahedron.'
  write ( *, '(a)' ) '  barycentric coordinates.'

  call r8mat_transpose_print ( dim_num, 4, t, '  Tetrahedron vertices' )

  call tetrahedron_sample ( t, test_num, p )

  call r8mat_transpose_print ( dim_num, test_num, p, '  Sample points:' )

  return
end
subroutine tetrahedron_shape_test ( )

!*****************************************************************************80
!
!! tetrahedron_shape_test() tests tetrahedron_shape().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) area
  integer edge_num
  integer face
  integer face_num
  integer, allocatable, dimension ( : ) ::  face_order
  integer face_order_max
  integer, allocatable, dimension ( :, : ) ::  face_point
  integer i
  integer j
  integer k
  real ( kind = rk ) normal(dim_num)
  integer point
  integer point_num
  real ( kind = rk ), allocatable, dimension ( :, : ) :: point_coord
  real ( kind = rk ), allocatable, dimension ( :, : ) :: v
  real ( kind = rk ) vave(dim_num)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'tetrahedron_shape_test():'
  write ( *, '(a)' ) '  tetrahedron_size() returns dimension information'
  write ( *, '(a)' ) '  of a tetrahedron.'
  write ( *, '(a)' ) '  tetrahedron_shape() returns face and order info'
  write ( *, '(a)' ) '  of a tetrahedron.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  We will use this information to compute the'
  write ( *, '(a)' ) '  areas and centers of each face.'

  call tetrahedron_size ( point_num, edge_num, face_num, face_order_max )

  write ( *, '(a)' ) ''
  write ( *, '(a,i8)' ) '  Number of vertices = ', point_num
  write ( *, '(a,i8)' ) '  Number of edges =    ', edge_num
  write ( *, '(a,i8)' ) '  Number of faces =    ', face_num
  write ( *, '(a,i8)' ) '  Maximum face order = ', face_order_max

  allocate ( face_order(1:face_num) )
  allocate ( face_point(1:face_order_max,1:face_num) )
  allocate ( point_coord(1:dim_num,1:point_num) )
  allocate ( v(1:dim_num,1:face_order_max) )

  call tetrahedron_shape ( point_num, face_num, face_order_max, &
    point_coord, face_order, face_point )

  call shape_print ( point_num, face_num, face_order_max, &
    point_coord, face_order, face_point )
!
!  Compute the area of each face.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Face  Order  Area'
  write ( *, '(a)' ) ''

  do face = 1, face_num

    do j = 1, face_order(face)
      point = face_point(j,face)
      v(1:dim_num,j) = point_coord(1:dim_num,point)
    end do

    call polygon_area_3d ( face_order(face), v, area, normal )

    write ( *, '(2x,i8,i7,f8.4)' ) face, face_order(face), area

  end do
!
!  Find the center of each face.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Face  Center'
  write ( *, '(a)' ) ''

  do i = 1, face_num

    vave(1:dim_num) = 0.0D+00

    do j = 1, face_order(i)
      k = face_point(j,i)
      vave(1:dim_num) = vave(1:dim_num) + point_coord(1:dim_num,k)
    end do

    vave(1:dim_num) = vave(1:dim_num) / real ( face_order(i), kind = rk )

    write ( *, '(2x,i8,3f8.4)' ) i, vave(1:dim_num)

  end do

  deallocate ( face_order )
  deallocate ( face_point )
  deallocate ( point_coord )
  deallocate ( v )

  return
end
subroutine tetrahedron_solid_angles_test ( )

!*****************************************************************************80
!
!! tetrahedron_solid_angles_test() tests tetrahedron_solid_angles();
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    28 May 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) angle(4)
  real ( kind = rk ), dimension(3,4) :: t1 = reshape ( (/&
     0.000000D+00,  0.942809D+00, -0.333333D+00, &
    -0.816496D+00, -0.816496D+00, -0.333333D+00, &
     0.816496D+00, -0.816496D+00, -0.333333D+00, &
     0.000000D+00,  0.000000D+00,  1.000000D+00 /), (/ 3, 4 /) )
  real ( kind = rk ), dimension(3,4) :: t2 = reshape ( (/&
     0.000000D+00,  0.000000D+00,  0.000000D+00, &
     1.000000D+00,  0.000000D+00,  0.000000D+00, &
     0.000000D+00,  1.000000D+00,  0.000000D+00, &
     0.000000D+00,  0.000000D+00,  1.000000D+00 /), (/ 3, 4 /) )
  real ( kind = rk ), dimension(3,4) :: t3 = reshape ( (/&
     0.000000D+00,  0.000000D+00,  0.000000D+00, &
     1.000000D+00,  0.000000D+00,  0.000000D+00, &
     0.000000D+00,  2.000000D+00,  0.000000D+00, &
     0.000000D+00,  0.000000D+00,  4.000000D+00 /), (/ 3, 4 /) )
  real ( kind = rk ), dimension(3,4) :: t4 = reshape ( (/&
     0.000000D+00,  0.000000D+00,  0.000000D+00, &
     1.000000D+00,  0.000000D+00,  0.000000D+00, &
     0.000000D+00,  1.000000D+00,  0.000000D+00, &
     1.000000D+00,  1.000000D+00,  1.000000D+00 /), (/ 3, 4 /) )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'tetrahedron_solid_angles_test():'
  write ( *, '(a)' ) '  tetrahedron_solid_angles() computes the solid angles'
  write ( *, '(a)' ) '  associated with the vertices of a tetrahedron.'

  call r8mat_transpose_print ( 3, 4, t1, '  Tetrahedron #1' )
  call tetrahedron_solid_angles ( t1, angle )
  call r8vec_print ( 4, angle, '  Solid angles for tetrahedron #1:' )
  
  call r8mat_transpose_print ( 3, 4, t2, '  Tetrahedron #2' )
  call tetrahedron_solid_angles ( t2, angle )
  call r8vec_print ( 4, angle, '  Solid angles for tetrahedron #2:' )

  call r8mat_transpose_print ( 3, 4, t3, '  Tetrahedron #3' )
  call tetrahedron_solid_angles ( t3, angle )
  call r8vec_print ( 4, angle, '  Solid angles for tetrahedron #3:' )

  call r8mat_transpose_print ( 3, 4, t4, '  Tetrahedron #4' )
  call tetrahedron_solid_angles ( t4, angle )
  call r8vec_print ( 4, angle, '  Solid angles for tetrahedron #4:' )

  return
end
subroutine tetrahedron_volume_test ( )

!*****************************************************************************80
!
!! tetrahedron_volume_test() tests tetrahedron_volume();
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ), dimension(dim_num,4) :: tetra = reshape ( (/&
     0.000000D+00,  0.942809D+00, -0.333333D+00, &
    -0.816496D+00, -0.816496D+00, -0.333333D+00, &
     0.816496D+00, -0.816496D+00, -0.333333D+00, &
     0.000000D+00,  0.000000D+00,  1.000000D+00 /), (/ dim_num, 4 /) )
  real ( kind = rk ) volume

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'tetrahedron_volume_test():'
  write ( *, '(a)' ) '  tetrahedron_volume() computes the volume of a tetrahedron.'

  call r8mat_transpose_print ( dim_num, 4, tetra, '  Tetrahedron vertices' )

  call tetrahedron_volume ( tetra, volume )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Volume = ', volume

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

