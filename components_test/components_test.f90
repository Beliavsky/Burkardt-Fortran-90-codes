program main

!*****************************************************************************80
!
!! components_test() tests components().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 July 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'components_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test components().'

  call components_1d_test ( )
  call components_2d_test ( )
  call components_3d_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'components_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine components_1d_test ( )

!*****************************************************************************80
!
!! components_1d_test() tests components_1d().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 July 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) , parameter :: n = 28

  integer ( kind = 4 ) :: a(n) = (/ &
    0, 0, 1, 2, 4, 0, 0, 4, 0, 0, &
    0, 8, 9, 9, 1, 2, 3, 0, 0, 5, &
    0, 1, 6, 0, 0, 0, 4, 0 /)
  integer ( kind = 4 ) c(n)
  integer ( kind = 4 ) component_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'components_1d_test():'
  write ( *, '(a)' ) '  components_1d() finds and labels connected'
  write ( *, '(a)' ) '  components in a 1D integer vector.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A:'
  write ( *, '(a)' ) ' '
  write ( *, '(4x,28i1)' ) a(1:n)

  call components_1d ( n, a, component_num, c )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of components = ', component_num
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  C:'
  write ( *, '(a)' ) ' '
  write ( *, '(4x,28i1)' ) c(1:n)

  return
end
subroutine components_2d_test ( )

!*****************************************************************************80
!
!! components_2d_test() tests components_2d().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 July 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 9
  integer ( kind = 4 ) , parameter :: n = 17

  integer ( kind = 4 ) :: a(m,n) = reshape ( (/ &
    0, 0, 0, 0, 0, 0, 0, 0, 0, &
    0, 0, 1, 0, 0, 1, 0, 0, 0, &
    0, 1, 1, 0, 1, 1, 1, 0, 0, &
    0, 1, 1, 1, 1, 1, 1, 0, 0, &
    0, 0, 1, 1, 1, 0, 0, 0, 0, &
    0, 0, 1, 1, 1, 0, 0, 0, 0, &
    0, 1, 1, 1, 0, 1, 0, 1, 0, &
    0, 1, 1, 0, 0, 1, 0, 1, 0, &
    0, 0, 1, 0, 0, 0, 0, 1, 0, &
    0, 0, 0, 0, 1, 0, 1, 1, 0, &
    0, 1, 0, 1, 1, 0, 1, 0, 0, &
    0, 1, 1, 1, 1, 1, 0, 0, 0, &
    0, 0, 1, 1, 0, 1, 0, 1, 0, &
    0, 0, 1, 1, 0, 1, 0, 1, 0, &
    0, 1, 1, 0, 1, 0, 1, 1, 0, &
    0, 1, 0, 0, 1, 0, 1, 1, 0, &
    0, 0, 0, 0, 0, 0, 0, 0, 0 /), &
    (/ m, n /) )
  integer ( kind = 4 ) c(m,n)
  integer ( kind = 4 ) component_num
  integer ( kind = 4 ) i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'components_2d_test():'
  write ( *, '(a)' ) '  components_2d() finds and labels connected'
  write ( *, '(a)' ) '  components in a 2D integer array.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A:'
  write ( *, '(a)' ) ' '
  do i = 1, m
    write ( *, '(4x,17i1)' ) a(i,1:n)
  end do

  call components_2d ( m, n, a, component_num, c )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of components = ', component_num
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  C:'
  write ( *, '(a)' ) ' '
  do i = 1, m
    write ( *, '(4x,17i1)' ) c(i,1:n)
  end do

  return
end
subroutine components_3d_test ( )

!*****************************************************************************80
!
!! components_3d_test() tests components_3d().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 July 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: l = 64
  integer ( kind = 4 ), parameter :: m = 64
  integer ( kind = 4 ) , parameter :: n = 26

  integer ( kind = 4 ) a(l,m,n)
  integer ( kind = 4 ) c(l,m,n)
  integer ( kind = 4 ) component_num
  character ( len = 80 ) filename
  integer ( kind = 4 ) i
  integer ( kind = 4 ), allocatable :: indices(:,:)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) n1
  integer ( kind = 4 ), allocatable :: s(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'components_3d_test():'
  write ( *, '(a)' ) '  components_3d() finds and labels connected'
  write ( *, '(a)' ) '  components in a 3D integer block.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i2,a,i2,a,i2)' ) &
    '  A is a 3D block of order ', l, ' * ',  m, ' * ', n

  a(1:l,1:m,1:n) = 0
!
!  Retrieve the indices of nonzero data in A by reading a file.
!
  filename = 'indices.txt'

  call i4mat_header_read ( filename, m1, n1 )

  allocate ( indices(m1,n1) )
  call i4mat_data_read ( filename, m1, n1, indices )

  do j1 = 1, n1
    i = indices(1,j1)
    j = indices(2,j1)
    k = indices(3,j1)
    a(i,j,k) = 1
  end do

  deallocate ( indices )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nonzero A values is ', sum ( a )
!
!  Determine the components.
!
  call components_3d ( l, m, n, a, component_num, c )

  allocate ( s(1:component_num) )

  s(1:component_num) = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of components = ', component_num

  do i = 1, l
    do j = 1, m
      do k = 1, n
        if ( c(i,j,k) /= 0 ) then
          s(c(i,j,k)) = s(c(i,j,k)) + 1
        end if
      end do
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Component  Size'
  write ( *, '(a)' ) ' '
  do i = 1, component_num
    write ( *, '(2x,i4,2x,i8)' ) i, s(i)
  end do
  write ( *, '(a)' ) '------  --------'
  write ( *, '(a6,2x,i8)' ) ' Total', sum ( s(1:component_num) )

  deallocate ( s )

  return
end
