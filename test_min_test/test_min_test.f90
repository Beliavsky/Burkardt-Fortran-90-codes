program main

!*****************************************************************************80
!
!! MAIN is the main program for test_min_test.
!
!  Discussion:
!
!    test_min_test tests the TEST_MIN library.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 February 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  call timestamp (  )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test_min_test'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TEST_MIN library.'

  call p00_title_test ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call p00_fmin_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test_min_test'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine p00_title_test ( )

!*****************************************************************************80
!
!! p00_title_test prints the title of each problem.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 February 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer problem_num
  integer problem
  character ( len = 50 ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'p00_title_test'
  write ( *, '(a)' ) '  For each problem, print the title.'
!
!  Get the number of problems.
!
  call p00_problem_num ( problem_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Problem Title'
  write ( *, '(a)' ) ' '

  do problem = 1, problem_num

    call p00_title ( problem, title )

    write ( *, '(2x,i8,2x,a)' ) problem, trim ( title )

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 evaluates the objective function at each starting point.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 August 2019
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) f_sol
  real ( kind = rk ) f_start
  real ( kind = rk ) p00_sol
  integer problem_num
  integer problem
  character ( len = 50 ) title
  real ( kind = rk ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  For each problem, evaluate the function'
  write ( *, '(a)' ) '  at the starting point and the solution.'
!
!  Get the number of problems.
!
  call p00_problem_num ( problem_num )

  do problem = 1, problem_num

    call p00_title ( problem, title )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Problem ', problem
    write ( *, '(2x,a)' ) trim ( title )
    write ( *, '(a)' ) ' '
 
    call p00_start ( problem, x )
    call p00_f ( problem, x, f_start )
    write ( *, '(4x,a,g16.8)' ) 'F(X_START)=', f_start

    x = p00_sol ( problem )
    call p00_f ( problem, x, f_sol )
    write ( *, '(4x,a,g16.8)' ) 'F(X_SOL)=  ', f_sol

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 compares the exact and approximate first derivatives.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 February 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) f1
  real ( kind = rk ) f1_dif
  integer problem_num
  integer problem
  character ( len = 50 ) title
  real ( kind = rk ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  For each problem, compare the exact and'
  write ( *, '(a)' ) '  approximate gradients at the starting point.'
!
!  Get the number of problems.
!
  call p00_problem_num ( problem_num )

  do problem = 1, problem_num

    call p00_title ( problem, title )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Problem ', problem
    write ( *, '(2x,a)' ) trim ( title )

    call p00_start ( problem, x )

    call p00_f1 ( problem, x, f1 )

    call p00_f1_dif ( problem, x, f1_dif )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,a)' ) 'X'
    write ( *, '(4x,5g16.8)' ) x
    write ( *, '(2x,a)' ) 'F''(X) (exact)'
    write ( *, '(4x,5g16.8)' ) f1
    write ( *, '(2x,a)' ) 'F''(X) (difference)'
    write ( *, '(4x,5g16.8)' ) f1_dif

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 compares the exact and approximate second derivatives.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 February 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) f2
  real ( kind = rk ) f2_dif
  integer problem_num
  integer problem
  character ( len = 50 ) title
  real ( kind = rk ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  For each problem, compare the exact and'
  write ( *, '(a)' ) '  approximate second derivatives at the starting point.'
!
!  Get the number of problems.
!
  call p00_problem_num ( problem_num )

  do problem = 1, problem_num

    call p00_title ( problem, title )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Problem ', problem
    write ( *, '(2x,a)' ) trim ( title )

    call p00_start ( problem, x )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  X:'
    write ( *, '(4x,5g16.8)' ) x

    call p00_f2 ( problem, x, f2 )

    write ( *, '(a)' ) '  F"(X) (exact):'
    write ( *, '(4x,6g13.5)' ) f2

    call p00_f2_dif ( problem, x, f2_dif )

    write ( *, '(a)' ) '  F"(X) (difference):'
    write ( *, '(4x,6g13.5)' ) f2_dif

  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 carries out a simple bisection method.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 February 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) fa
  real ( kind = rk ) fb
  real ( kind = rk ) fc
  real ( kind = rk ) fd
  real ( kind = rk ) fe
  integer i
  integer, parameter :: max_step = 10
  integer problem_num
  integer problem
  character ( len = 50 ) title
  real ( kind = rk ) xa
  real ( kind = rk ) xb
  real ( kind = rk ) xc
  real ( kind = rk ) xd
  real ( kind = rk ) xe

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  For each problem, take a few steps of '
  write ( *, '(a)' ) '  the bisection method.'
!
!  Get the number of problems.
!
  call p00_problem_num ( problem_num )

  do problem = 1, problem_num

    call p00_title ( problem, title )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Problem ', problem
    write ( *, '(2x,a)' ) trim ( title )

    call p00_interval ( problem, xa, xc )
    xb = 0.5D+00 * ( xa + xc )
    call p00_f ( problem, xa, fa )
    call p00_f ( problem, xc, fc )
    call p00_f ( problem, xb, fb )

    i = 0
    write ( *, '(a)' ) ' '
    write ( *, '(i6)' ) i
    write ( *, '(a,3g16.10)' ) '  X:', xa, xb, xc
    write ( *, '(a,3g16.10)' ) '  F:', fa, fb, fc

    do i = 1, max_step

      xd = 0.5D+00 * ( xa + xb )
      call p00_f ( problem, xd, fd )

      xe = 0.5D+00 * ( xb + xc )
      call p00_f ( problem, xe, fe )

      if ( fd <= fb ) then
        xc = xb
        fc = fb
        xb = xd
        fb = fd
      else if ( fe <= fb ) then
        xa = xb
        fa = fb
        xb = xe
        fb = fe
      else
        xa = xd
        fa = fd
        xc = xe
        fc = fe
      end if

      write ( *, '(i6)' ) i
      write ( *, '(a,3g16.10)' ) '  X:', xa, xb, xc
      write ( *, '(a,3g16.10)' ) '  F:', fa, fb, fc

    end do

  end do

  return
end
subroutine p00_fmin_test ( )

!*****************************************************************************80
!
!! p00_fmin_test tests p00_fmin.
!
!  Discussion:
!
!    p00_fmin carries out a version of Brent's derivative-free minimizer.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 August 2019
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) fa
  real ( kind = rk ) fb
  real ( kind = rk ) fx
  real ( kind = rk ) p00_fmin
  integer problem_num
  integer problem
  character ( len = 50 ) title
  real ( kind = rk ), parameter :: tol = 0.00000001D+00
  real ( kind = rk ) x
  real ( kind = rk ) xa
  real ( kind = rk ) xb

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'p00_fmin_test'
  write ( *, '(a)' ) '  For each problem, use Brent''s method.'
!
!  Get the number of problems.
!
  call p00_problem_num ( problem_num )

  do problem = 1, problem_num

    call p00_title ( problem, title )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Problem ', problem
    write ( *, '(2x,a)' ) trim ( title )

    call p00_interval ( problem, xa, xb )

    call p00_f ( problem, xa, fa )
    call p00_f ( problem, xb, fb )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Initial interval [A,B]:'
    write ( *, '(a)' ) ' '
    write ( *, '(a,g16.10,14x,g16.10)' ) '   A,       B:', xa,     xb
    write ( *, '(a,g16.10,14x,g16.10)' ) '  FA,      FB:', fa,     fb

    x = p00_fmin ( xa, xb, problem, tol )

    call p00_f ( problem, xa, fa )
    call p00_f ( problem, xb, fb )
    call p00_f ( problem, x, fx )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Final interval [A,X*,B]:'
    write ( *, '(a)' ) ' '
    write ( *, '(a,g16.10,g16.10,g16.10)' ) '   A,  X*,  B:', xa, x,  xb
    write ( *, '(a,g16.10,g16.10,g16.10)' ) '  FA, FX*, FB:', fa, fx, fb

  end do

  return
end
