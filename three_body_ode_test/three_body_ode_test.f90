program main

!*****************************************************************************80
!
!! MAIN is the main program for three_body_ode.
!
!  Discussion:
!
!    three_body_ode uses RKF45 as an integrator for the simple version
!    of the three-body problem.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 April 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'three_body_ode_test'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Simulate the behavior of three bodies which are'
  write ( *, '(a)' ) '  constrained to lie in a plane, moving under the'
  write ( *, '(a)' ) '  influence of gravity.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Use RKF45 for the ODE integrator.'

  call three_body_ode_run ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'three_body_ode_test'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine three_body_ode_run ( )

!*****************************************************************************80
!
!! three_body_ode_run runs the simple three body ODE system.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 April 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: neqn = 12
  integer, parameter :: step_num = 630

  real ( kind = rk ) abserr
  integer flag
  real ( kind = rk ) m0
  real ( kind = rk ) m1
  real ( kind = rk ) m2
  real ( kind = rk ) relerr
  external three_body_deriv
  integer step
  real ( kind = rk ) t
  character ( len = 80 ) :: t_filename = 'three_body_ode_t.txt'
  real ( kind = rk ) t_out
  real ( kind = rk ) t_start
  real ( kind = rk ) t_stop
  real ( kind = rk ) ts(0:step_num)
  real ( kind = rk ) y(neqn)
  character ( len = 80 ) :: y_filename = 'three_body_ode_y.txt'
  real ( kind = rk ) yp(neqn)
  real ( kind = rk ) ys(neqn,0:step_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'three_body_ode_run'
  write ( *, '(a)' ) '  Simulate the planar three-body problem as an ODE system'
  write ( *, '(a)' ) '  using RKF45 for the ODE integration.'

  m0 = 5.0D+00
  m1 = 3.0D+00
  m2 = 4.0D+00

  abserr = 1.0D-10
  relerr = 1.0D-10

  flag = 1

  t_start = 0.0D+00
  t_stop = 63.0D+00

  t = 0.0D+00
  t_out = 0.0D+00

  y(1:neqn) = (/ 1.0D+00, -1.0D+00,  0.0D+00,  0.0D+00, &
                 1.0D+00,  3.0D+00,  0.0D+00,  0.0D+00, &
                -2.0D+00, -1.0D+00,  0.0D+00,  0.0D+00 /)

  call three_body_deriv ( t, y, yp )

  ys(1:neqn,0) = y(1:neqn)
  ts(0) = t

  do step = 1, step_num

    t = ( real ( step_num - step + 1, kind = rk ) * t_start &
        + real (            step - 1, kind = rk ) * t_stop ) &
        / real ( step_num,            kind = rk )

    t_out = ( real ( step_num - step, kind = rk ) * t_start &
            + real (            step, kind = rk ) * t_stop ) &
            / real ( step_num,        kind = rk )

    call r8_rkf45 ( three_body_deriv, neqn, y, yp, t, t_out, relerr, abserr, flag )

    if ( abs ( flag ) /= 2 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'three_body_ode_run - Warning!'
      write ( *, '(a,i4,a,g14.6)' ) '  Output value of FLAG = ', flag, &
        ' at T_OUT = ', t_out
    end if

    ys(1:neqn,step) = y(1:neqn)
    ts(step) = t_out

  end do

  call r8mat_write ( t_filename, 1, step_num + 1, ts )
  call r8mat_write ( y_filename, neqn, step_num + 1, ys )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'three_body_ode_run:'
  write ( *, '(a)' ) '  Time data written to "' // trim ( t_filename ) // '".'
  write ( *, '(a)' ) '  Solution data written to "' // trim ( y_filename ) // '".'

  return
end
