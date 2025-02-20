program main

!*****************************************************************************80
!
!! ODE_TEST() tests ODE().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ODE_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the ODE library.'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ODE_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests ODE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: neqn = 2

  real ( kind = rk ) abserr
  external f01
  integer i
  integer iflag
  integer iwork(5)
  real ( kind = rk ) :: pi = 3.141592653589793D+00
  real ( kind = rk ) relerr
  integer, parameter :: step_num = 12
  real ( kind = rk ) t
  real ( kind = rk ) tout
  real ( kind = rk ) work(100+21*neqn)
  real ( kind = rk ) y(neqn)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  ODE solves a system of ordinary differential'
  write ( *, '(a)' ) '  equations.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      T           Y(1)         Y(2)'
  write ( *, '(a)' ) ' '

  abserr = 0.00001D+00
  relerr = 0.00001D+00

  iflag = 1

  t = 0.0D+00
  y(1) = 1.0D+00
  y(2) = 0.0D+00

  write ( *, '(2x,f8.4,2x,2g14.6)' ) t, y(1:neqn)

  do i = 1, step_num

    tout = real ( i, kind = rk ) * 2.0D+00 * pi / real ( step_num, kind = rk )

    call ode ( f01, neqn, y, t, tout, relerr, abserr, iflag, work, iwork )

    if ( iflag /= 2 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01 - Fatal error!'
      write ( *, '(a,i8)' ) '  ODE returned IFLAG = ', iflag
      exit
    end if

    write ( *, '(2x,f8.4,2x,2g14.6)' ) t, y(1:neqn)

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests ODE by integrating in the NEGATIVE time direction.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: neqn = 2

  real ( kind = rk ) abserr
  external f01
  integer i
  integer iflag
  integer iwork(5)
  real ( kind = rk ) :: pi = 3.141592653589793D+00
  real ( kind = rk ) relerr
  integer, parameter :: step_num = 12
  real ( kind = rk ) t
  real ( kind = rk ) tout
  real ( kind = rk ) work(100+21*neqn)
  real ( kind = rk ) y(neqn)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  ODE solves a system of ordinary differential'
  write ( *, '(a)' ) '  equations.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this example, we integrate in the negative'
  write ( *, '(a)' ) '  time direction.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      T           Y(1)         Y(2)'
  write ( *, '(a)' ) ' '

  abserr = 0.00001D+00
  relerr = 0.00001D+00

  iflag = 1

  t = 0.0D+00
  y(1) = 1.0D+00
  y(2) = 0.0D+00

  write ( *, '(2x,f8.4,2x,2g14.6)' ) t, y(1:neqn)

  do i = 1, step_num

    tout = - real ( i, kind = rk ) * 2.0D+00 * pi / real ( step_num, kind = rk )

    call ode ( f01, neqn, y, t, tout, relerr, abserr, iflag, work, iwork )

    if ( iflag /= 2 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02 - Fatal error!'
      write ( *, '(a,i8)' ) '  ODE returned IFLAG = ', iflag
      exit
    end if

    write ( *, '(2x,f8.4,2x,2g14.6)' ) t, y(1:neqn)

  end do

  return
end
subroutine f01 ( t, y, yp )

!*****************************************************************************80
!
!! F01 supplies the right hand side of the ODE for problem 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) T, the time.
!
!    Input, real ( kind = rk ) Y(), the dependent variable.
!
!    Output, real ( kind = rk ) YP(), the value of the derivative.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) t
  real ( kind = rk ) y(2)
  real ( kind = rk ) yp(2)

  yp(1) =   y(2)
  yp(2) = - y(1)

  return
end
