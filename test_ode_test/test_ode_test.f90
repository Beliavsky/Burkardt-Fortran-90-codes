program main

!*****************************************************************************80
!
!! test_ode_test() tests test_ode().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 October 2023
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test_ode_test()'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test test_ode().'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test05 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test_ode_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 simply lists the problems with titles and sizes.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer neqn
  integer test
  integer test_num
  character ( len = 80 ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  List the problem titles and sizes.'
!
!  Find out how many test problems are available.
!
  call p00_test_num ( test_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) '  There are ', test_num, ' test problems.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test  Size  Title'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    call p00_title ( test, title )
    call p00_neqn ( test, neqn )
    write ( *, '(2x,i4,2x,i4,2x,a)' ) test, neqn, trim ( title )
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 solves most of the problems using an Euler method.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: step_num = 500
  integer test
  integer test_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Solve problems using an Euler method.'
  write ( *, '(a,i8)' ) '  The number of steps taken is ', step_num

  call p00_test_num ( test_num )

  write ( *, '(a,i8)' ) '  The number of tests available is ', test_num
!
!  Solve each problem.
!
  do test = 1, test_num

    if ( test == 32 .or. test == 36 .or. test == 37 ) then

    else
      call euler_test ( test, step_num )
    end if

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 solves most of the problems using a Runge-Kutta method.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer order
  integer, parameter :: step_num = 500
  integer test
  integer test_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Solve problems using a Runge-Kutta method.'
  write ( *, '(a,i8)' ) '  The number of steps taken is ', step_num

  call p00_test_num ( test_num )

  write ( *, '(a,i8)' ) '  The number of tests available is ', test_num
!
!  Solve each problem.
!
  order = 3

  do test = 1, test_num

    call rk_test ( test, step_num, order )

  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 compares the Jacobian to a finite difference estimate.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 February 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) diff
  real ( kind = rk ) dy
  real ( kind = rk ) e
  real ( kind = rk ), allocatable :: f1(:)
  real ( kind = rk ), allocatable :: f2(:)
  real ( kind = rk ), allocatable :: jac1(:,:)
  real ( kind = rk ), allocatable :: jac2(:,:)
  integer j
  integer neqn
  real ( kind = rk ) r8_sign
  real ( kind = rk ) r8_uniform_ab
  real ( kind = rk ) r8mat_norm_fro_affine
  integer seed
  real ( kind = rk ) t_start
  real ( kind = rk ) t_stop
  real ( kind = rk ) t1
  integer test
  integer test_num
  real ( kind = rk ), allocatable :: y_start(:)
  real ( kind = rk ), allocatable :: y_stop(:)
  real ( kind = rk ), allocatable :: y1(:)
  real ( kind = rk ), allocatable :: y2(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  At a random time T in the time interval,'
  write ( *, '(a)' ) '  and a random vector Y, compare the jacobian dF/dY'
  write ( *, '(a)' ) '  and a finite difference estimate.'
  call p00_test_num ( test_num )

  write ( *, '(a,i8)' ) '  The number of tests available is ', test_num
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test    Difference'
  write ( *, '(a)' ) ' '

  e = sqrt ( epsilon ( dy ) )
!
!  Solve each problem.
!
  do test = 1, test_num

    call p00_neqn ( test, neqn )

    allocate ( f1(1:neqn) )
    allocate ( f2(1:neqn) )
    allocate ( jac1(1:neqn,1:neqn) )
    allocate ( jac2(1:neqn,1:neqn) )
    allocate ( y_start(1:neqn) )
    allocate ( y_stop(1:neqn) )
    allocate ( y1(1:neqn) )
    allocate ( y2(1:neqn) )

    call p00_start ( test, neqn, t_start, y_start )
    call p00_stop ( test, neqn, t_stop, y_stop )

    seed = 123456789
    t1 = r8_uniform_ab ( t_start, t_stop, seed )
    call r8vec_uniform_abvec ( neqn, y_start, y_stop, seed, y1 )

    call p00_jac ( test, neqn, t1, y1, jac1 )
    jac2(1:neqn,1:neqn) = 0.0D+00

    call p00_fun ( test, neqn, t1, y1, f1 )

    do j = 1, neqn
      y2(1:neqn) = y1(1:neqn)
      dy = r8_sign ( y2(j) ) * e * ( abs ( y2(j) ) + 1.0D+00 )
      y2(j) = y2(j) + dy
      call p00_fun ( test, neqn, t1, y2, f2 )
      jac2(1:neqn,j) = ( f2(1:neqn) - f1(1:neqn) ) / dy
    end do

    diff = r8mat_norm_fro_affine ( neqn, neqn, jac1, jac2 )

    write ( *, '(2x,i2,2x,g14.6)' ) test, diff

    deallocate ( f1 )
    deallocate ( f2 )
    deallocate ( jac1 )
    deallocate ( jac2 )
    deallocate ( y_start )
    deallocate ( y_stop )
    deallocate ( y1 )
    deallocate ( y2 )

  end do

  return
end


