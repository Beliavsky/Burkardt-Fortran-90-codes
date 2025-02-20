program main

!*****************************************************************************80
!
!! hypersphere_test() tests hypersphere().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 December 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'hypersphere_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test hypersphere().'

  call hypersphere_test01 ( )
  call hypersphere_test02 ( )
  call hypersphere_test03 ( )
  call hypersphere_test04 ( )
  call hypersphere_test05 ( )
  call hypersphere_test06 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'hypersphere_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine hypersphere_test01 ( )

!*****************************************************************************80
!
!! hypersphere_test01() tests the coordinate conversion routines.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 May 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: c(:)
  real ( kind = rk ) err
  integer m
  integer n
  real ( kind = rk ), allocatable :: r(:)
  real ( kind = rk ) r8mat_norm_fro_affine
  integer test
  real ( kind = rk ), allocatable :: theta(:,:)
  real ( kind = rk ), allocatable :: x(:,:)
  real ( kind = rk ), allocatable :: x2(:,:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'hypersphere_test01():'
  write ( *, '(a)' ) '  Test the coordinate conversion routines:'
  write ( *, '(a)' ) '  CARTESIAN_TO_HYPERSPHERE: X       -> R,Theta'
  write ( *, '(a)' ) '  HYPERSPHERE_TO_CARTESIAN: R,Theta -> X.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Pick a random X, and compute X2 by converting X'
  write ( *, '(a)' ) '  to hypersphere and back.  Consider norm of difference.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  M    || X - X2 ||'

  n = 1

  allocate ( r(n) )

  do m = 1, 5

    write ( *, '(a)' ) ''

    allocate ( x(1:m,1:n) )
    allocate ( x2(1:m,1:n) )
    allocate ( c(1:m) )
    allocate ( theta(m-1,n) )

    do test = 1, 5
      call random_number ( harvest = x(1:m,1:n) )
      call random_number ( harvest = c(1:m) )
      call cartesian_to_hypersphere ( m, n, c, x, r, theta )
      call hypersphere_to_cartesian ( m, n, c, r, theta, x2 )
      err = r8mat_norm_fro_affine ( m, n, x, x2 )
      write ( *, '(2x,i2,2x,g14.6)' ) m, err
    end do

    deallocate ( c )
    deallocate ( theta )
    deallocate ( x )
    deallocate ( x2 )

  end do

  deallocate ( r )

  return
end
subroutine hypersphere_test02 ( )

!*****************************************************************************80
!
!! hypersphere_test02() tests HYPERSPHERE_01_SURFACE_UNIFORM.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 May 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n
  integer test
  real ( kind = rk ), allocatable :: x(:,:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'hypersphere_test02():'
  write ( *, '(a)' ) '  HYPERSPHERE_01_SURFACE_UNIFORM samples uniformly from the'
  write ( *, '(a)' ) '  surface of the unit hypersphere'

  n = 1
  do m = 1, 5
    allocate ( x(1:m,1:n) )
    do test = 1, 3
      call hypersphere_01_surface_uniform ( m, n, x )
      call r8vec_transpose_print ( m, x, '  Random hypersphere point:' )
    end do
    deallocate ( x )
  end do

  return
end
subroutine hypersphere_test03 ( )

!*****************************************************************************80
!
!! hypersphere_test03() tests HYPERSPHERE_01_AREA.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 May 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) area
  real ( kind = rk ) area2
  real ( kind = rk ) hypersphere_01_area
  integer m
  integer n_data

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'hypersphere_test03():'
  write ( *, '(a)' ) '  HYPERSPHERE_01_AREA evaluates the area of the unit'
  write ( *, '(a)' ) '  hypersphere in M dimensions.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '       M      Exact       Computed'
  write ( *, '(a)' ) '              Area        Area'
  write ( *, '(a)' ) ''

  n_data = 0

  do

    call hypersphere_01_area_values ( n_data, m, area )

    if ( n_data == 0 ) then
      exit
    end if

    area2 = hypersphere_01_area ( m )

    write ( *, '(2x,i6,2x,f10.4,2x,f10.4)' ) m, area, area2

  end do

  return
end
subroutine hypersphere_test04 ( )

!*****************************************************************************80
!
!! hypersphere_test04() tests HYPERSPHERE_01_VOLUME.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 December 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) hypersphere_01_volume
  integer m
  integer n_data
  real ( kind = rk ) volume
  real ( kind = rk ) volume2

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'hypersphere_test04():'
  write ( *, '(a)' ) '  HYPERSPHERE_01_VOLUME evaluates the area of the unit'
  write ( *, '(a)' ) '  hypersphere in M dimensions.'
  write ( *, '(a)' ) '  HYPERSPHERE_01_VOLUME_VALUES returns some test values.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '       M      Exact       Computed'
  write ( *, '(a)' ) '              Volume      Volume'
  write ( *, '(a)' ) ''

  n_data = 0

  do

    call hypersphere_01_volume_values ( n_data, m, volume )

    if ( n_data == 0 ) then
      exit
    end if

    volume2 = hypersphere_01_volume ( m )

    write ( *, '(2x,i6,2x,f10.4,2x,f10.4)' ) m, volume, volume2

  end do

  return
end
subroutine hypersphere_test05 ( )

!*****************************************************************************80
!
!! hypersphere_test05() tests HYPERSPHERE_AREA, HYPERSPHERE_VOLUME.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 May 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) area
  real ( kind = rk ) hypersphere_area
  real ( kind = rk ) hypersphere_volume
  integer m
  real ( kind = rk ) r
  real ( kind = rk ) volume

  r = 1.5D+00

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'hypersphere_test05():'
  write ( *, '(a)' ) '  For a hypersphere in M dimensions:'
  write ( *, '(a)' ) '  HYPERSPHERE_AREA computes the area'
  write ( *, '(a)' ) '  HYPERSPHERE_VOLUME computes the volume.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Notice that both quantities eventually decrease!'
  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  We use a radius of R = ', r
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '    M        Area          Volume    Area / Volume '
  write ( *, '(a)' ) ''

  do m = 1, 20
    area = hypersphere_area ( m, r )
    volume = hypersphere_volume ( m, r )
    write ( *, '(2x,i3,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      m, area, volume, area / volume
  end do

  return
end
subroutine hypersphere_test06 ( )

!*****************************************************************************80
!
!! hypersphere_test06() tests the stereographic mapping.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 May 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) err
  integer m
  integer n
  real ( kind = rk ) r8mat_norm_fro_affine
  integer test
  real ( kind = rk ), allocatable :: x1(:,:)
  real ( kind = rk ), allocatable :: x2(:,:)
  real ( kind = rk ), allocatable :: x3(:,:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'hypersphere_test06()'
  write ( *, '(a)' ) '  Test the stereographic mapping:'
  write ( *, '(a)' ) '  HYPERSPHERE_STEREOGRAPH maps hypersphere points to the plane.'
  write ( *, '(a)' ) '  HYPERSPHERE_STEREOGRAPH_INVERSE inverts the mapping.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Pick a random X1 on the hypersphere.'
  write ( *, '(a)' ) '  Map it to a point X2 on the plane.'
  write ( *, '(a)' ) '  Map it back to a point X3 on the hypersphere.'
  write ( *, '(a)' ) '  Consider norm of difference.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  M    || X1 - X3 ||'

  n = 1
  do m = 2, 5 
    write ( *, '(a)' ) ''
    allocate ( x1(m,n) )
    allocate ( x2(m-1,n) )
    allocate ( x3(m,n) )
    do test = 1, 5
      call hypersphere_01_surface_uniform ( m, n, x1 )
      call hypersphere_stereograph ( m, n, x1, x2 )
      call hypersphere_stereograph_inverse ( m, n, x2, x3 )
      err = r8mat_norm_fro_affine ( m, n, x1, x3 )
      write ( *, '(2x,i2,2x,g14.6)' ) m, err
    end do
    deallocate ( x1 )
    deallocate ( x2 )
    deallocate ( x3 )
  end do

  return
end

