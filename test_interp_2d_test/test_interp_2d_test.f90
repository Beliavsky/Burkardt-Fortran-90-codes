program main

!*****************************************************************************80
!
!! TEST_INTERP_2D_TEST tests the TEST_INTERP_2D library.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 January 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_INTERP_2D_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TEST_INTERP_2D library.'
  write ( *, '(a)' ) '  The test requires access to the R8LIB library.'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_INTERP_2D_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 simply prints the title of each grid and function.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 January 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer f_num
  integer fi
  character ( len = 50 ) ft
  integer g_num
  integer gi
  character ( len = 50 ) gt

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  For each grid and function, print the title.'

  call g00_num ( g_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  GRIDS:'
  write ( *, '(a)' ) '  Index  Title'
  write ( *, '(a)' ) ' '

  do gi = 1, g_num

    call g00_title ( gi, gt )

    write ( *, '(2x,i2,2x,a)' ) gi, gt

  end do

  call f00_num ( f_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FUNCTIONS:'
  write ( *, '(a)' ) '  Index  Title'
  write ( *, '(a)' ) ' '

  do fi = 1, f_num

    call f00_title ( fi, ft )

    write ( *, '(2x,i2,2x,a)' ) fi, ft

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 samples each function using each grid.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    28 January 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: f(:)
  real ( kind = rk ) f_ave
  real ( kind = rk ) f_max
  real ( kind = rk ) f_min
  integer f_num
  integer fi
  character ( len = 50 ) ft
  integer g_num
  integer gi
  integer gn
  character ( len = 50 ) gt
  real ( kind = rk ), allocatable :: gx(:)
  real ( kind = rk ), allocatable :: gy(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Sample each function over each grid.'

  call g00_num ( g_num )
  call f00_num ( f_num )

  do fi = 1, f_num

    call f00_title ( fi, ft )
    write ( *, '(a)' ) ' '
    write ( *, '(2x,i2,2x,a)' ) fi, trim ( ft )
    write ( *, '(a)' ) '        Grid Title                     Min(F)' // &
      '          Ave(F)           Max(F)'
    write ( *, '(a)' ) ' '

    do gi = 1, g_num

      call g00_title ( gi, gt )
      call g00_size ( gi, gn )

      allocate ( gx(1:gn) )
      allocate ( gy(1:gn) )
      allocate ( f(1:gn) )

      call g00_xy ( gi, gn, gx, gy )

      call f00_f0 ( fi, gn, gx, gy, f )

      f_max = maxval ( f(1:gn) )
      f_min = minval ( f(1:gn) )
      f_ave = sum ( f(1:gn) ) / real ( gn, kind = rk )

      write ( *, '(2x,i4,2x,a25,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
        gi, gt, f_min, f_ave, f_max

      deallocate ( f )
      deallocate ( gx )
      deallocate ( gy )

    end do

  end do

  return
end
