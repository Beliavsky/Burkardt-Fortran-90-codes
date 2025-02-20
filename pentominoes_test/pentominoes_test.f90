program main

!******************************************************************************/
!
!! MAIN is the main program for PENTOMINOES_TEST.
!
!  Discussion:
!
!    PENTOMINOES_TEST tests the PENTOMINOES library.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 April 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'PENTOMINOES_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the PENTOMINOES library.'

  call pentomino_matrix_test ( )
  call pentomino_plot_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'PENTOMINOES_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine pentomino_matrix_test ( )

!*****************************************************************************80
!
!! PENTOMINO_MATRIX_TEST tests PENTOMINO_MATRIX.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 April 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  character ( len = 1 ) name
  integer ( kind = 4 ) p(5*5)
  integer ( kind = 4 ) p_m
  integer ( kind = 4 ) p_n
  character ( len = 1 ) :: pentominoes(12) = (/ &
    'F', 'I', 'L', 'N', 'P', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z' /)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'PENTOMINO_MATRIX_TEST'
  write ( *, '(a)' ) '  PENTOMINO_MATRIX returns a 0/1 matrix representing a pentomino.'

  do k = 1, 12
    name = pentominoes(k)
    call pentomino_matrix ( name, p_m, p_n, p )
    write ( *, '(a)' ) ''
    write ( *, '(a,i2,a,i2,a)' ) &
      '  // trim ( name ) // pentomino (', p_m, ',', p_n, ')'
    write ( *, '(a)' ) ''
    do i = 1, p_m
      write ( *, '(a)', advance = 'no' ) '    '
      do j = 1, p_n
        write ( *, '(i1)', advance = 'no' ) p(i+j*p_m)
      end do
      write ( *, '(a)' )
    end do
  end do

  return
end
subroutine pentomino_plot_test ( )

!*****************************************************************************80
!
!! PENTOMINO_PLOT_TEST tests PENTOMINO_PLOT.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 April 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) k
  character ( len = 1 ) name
  integer ( kind = 4 ) p(5*5)
  integer ( kind = 4 ) p_m
  integer ( kind = 4 ) p_n
  character ( len = 1 ) :: pentominoes(12) = (/ &
    'F', 'I', 'L', 'N', 'P', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z' /)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'PENTOMINO_PLOT_TEST'
  write ( *, '(a)' ) '  PENTOMINO_PLOT plots a pentomino.'

  do k = 1, 12
    name = pentominoes(k)
    call pentomino_matrix ( name, p_m, p_n, p )
    call pentomino_plot ( p_m, p_n, p, name )
  end do

  return
end
