program main

!*****************************************************************************80
!
!! pce_ode_hermite_test() tests pce_ode_hermite().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 March 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( );
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'pce_ode_hermite_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test pce_ode_hermite().'

  call pce_ode_hermite_test01 ( )
  call pce_ode_hermite_test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'pce_ode_hermite_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine pce_ode_hermite_test01 ( )

!*****************************************************************************80
!
!! pce_ode_hermite_test01() runs a test problem with pce_ode_hermite().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 March 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: np = 4
  integer, parameter :: nt = 200

  real ( kind = rk ) alpha_mu
  real ( kind = rk ) alpha_sigma
  integer i
  real ( kind = rk ) t(0:nt)
  real ( kind = rk ) tf
  real ( kind = rk ) ti
  real ( kind = rk ) u(0:nt,0:np)
  real ( kind = rk ) uex(0:nt)
  real ( kind = rk ) ui

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'pce_ode_hermite_test01():'
  write ( *, '(a)' ) '  pce_ode_hermite() to compute a polynomial chaos expansion'
  write ( *, '(a)' ) '  for the ODE:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    u'' = - alpha * u,'
  write ( *, '(a)' ) '    u(0) = 1.'

  ti = 0.0D+00
  tf = 2.0D+00
  ui = 1.0D+00
  alpha_mu = 0.0D+00
  alpha_sigma = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Initial time         TI = ', ti
  write ( *, '(a,g14.6)' ) '  Final time           TF = ', tf
  write ( *, '(a,i6)' ) '  Number of time steps NT = ', nt
  write ( *, '(a,g14.6)' ) '  Initial condition    UI = ', ui
  write ( *, '(a,i6)' ) '  Expansion degree     NP = ', np
  write ( *, '(a,g14.6)' ) '  E(ALPHA)       ALPHA_MU = ', alpha_mu
  write ( *, '(a,g14.6)' ) '  STD(ALPHA)  ALPHA_SIGMA = ', alpha_sigma

  call pce_ode_hermite ( ti, tf, nt, ui, np, alpha_mu, alpha_sigma, t, u )
!
!  Evaluate the exact expected value function.
!
  uex(0:nt) = ui * exp ( t(0:nt)**2 / 2.0D+00 )
!
!  Compare the first computed component against the exact expected value.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' i  T(i)  E(U(T(i)))    U(T(i),0)'
  write ( *, '(a)' ) ' '
  do i = 0, nt, 10
    write ( *, '(2x,i4,2x,f6.3,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      i, t(i), uex(i), u(i,0), abs ( uex(i) - u(i,0) )
  end do

  return
end
subroutine pce_ode_hermite_test02 ( )

!*****************************************************************************80
!
!! pce_ode_hermite_test02() looks at convergence behavior for a fixed time.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 March 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nt = 2000

  real ( kind = rk ) alpha_mu
  real ( kind = rk ) alpha_sigma
  real ( kind = rk ) ep(0:5)
  integer np
  real ( kind = rk ) t(0:nt)
  real ( kind = rk ) tf
  real ( kind = rk ) ti
  real ( kind = rk ), allocatable :: u(:,:)
  real ( kind = rk ) uexf
  real ( kind = rk ) ui

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'pce_ode_hermite_test02():'
  write ( *, '(a)' ) '  Examine convergence behavior as the PCE degree increases:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    u'' = - alpha * u,'
  write ( *, '(a)' ) '    u(0) = 1.'

  ti = 0.0D+00
  tf = 2.0D+00
  ui = 1.0D+00
  alpha_mu = 0.0D+00
  alpha_sigma = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Initial time         TI = ', ti
  write ( *, '(a,g14.6)' ) '  Final time           TF = ', tf
  write ( *, '(a,i6)' ) '  Number of time steps NT = ', nt
  write ( *, '(a,g14.6)' ) '  Initial condition    UI = ', ui
  write ( *, '(a,g14.6)' ) '  E(ALPHA)       ALPHA_MU = ', alpha_mu
  write ( *, '(a,g14.6)' ) '  STD(ALPHA)  ALPHA_SIGMA = ', alpha_sigma

  uexf = ui * exp ( tf**2 / 2.0D+00 )

  do np = 0, 5

    allocate ( u(0:nt,0:np) )

    call pce_ode_hermite ( ti, tf, nt, ui, np, alpha_mu, alpha_sigma, t, u )

    ep(np) = abs ( uexf - u(nt,0) )

    deallocate ( u )

  end do
!
!  Print error in expected value as a function of the PCE degree.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    NP     Error(NP)     Log(Error(NP))'
  write ( *, '(a)' ) ' '
  do np = 0, 5
    write ( *, '(2x,i4,2x,g14.6,2x,g14.6)' ) np, ep(np), log ( ep(np) )
  end do

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
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

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

