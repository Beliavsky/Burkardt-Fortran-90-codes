subroutine flame_exact ( n, t, y )

!*****************************************************************************80
!
!! flame_exact() evaluates the exact solution for the flame ODE.
!
!  Discussion:
!
!    Evaluating the exact solution requires the Lambert W function.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 April 2021
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Cleve Moler,
!    Cleve's Corner: Stiff Differential Equations,
!    MATLAB News and Notes,
!    May 2003, pages 12-13.
!
!  Input:
!
!    integer n: the number of times.
!
!    real ( kind = rk ) t(n): the times.
!
!  Output:
!
!    real ( kind = rk ) y(n), the exact solution values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  interface
    subroutine flame_parameters ( t0_in, tstop_in, y0_in, t0_out, &
      tstop_out, y0_out )
      integer, parameter :: rk = kind ( 1.0D+00 )
      real ( kind = rk ), optional :: t0_in
      real ( kind = rk ), optional :: t0_out
      real ( kind = rk ), optional :: tstop_in
      real ( kind = rk ), optional :: tstop_out
      real ( kind = rk ), optional :: y0_in
      real ( kind = rk ), optional :: y0_out
    end subroutine
  end interface

  integer n

  real ( kind = rk ) a
  integer i
  integer l
  real ( kind = rk ) lambert_w
  integer nb
  real ( kind = rk ) t(n)
  real ( kind = rk ) t0
! real ( kind = rk ) t0_out
  real ( kind = rk ) tstop
! real ( kind = rk ) tstop_out
  real ( kind = rk ) x
  real ( kind = rk ) y(n)
  real ( kind = rk ) y0
! real ( kind = rk ) y0_out

  call flame_parameters ( t0_out = t0, y0_out = y0, tstop_out = tstop )

  a = ( 1.0D+00 - y0 ) / y0

  nb = 0
  l = 0

  do i = 1, n
    x = a * exp ( a - ( t(i) - t0 ) )
    y(i) = 1.0D+00 / ( lambert_w ( x, nb, l ) + 1.0D+00 )
  end do

  return
end
subroutine flame_parameters ( t0_in, y0_in, tstop_in, t0_out, y0_out, tstop_out )

!*****************************************************************************80
!
!! flame_parameters() returns parameters for flame_ode().
!
!  Discussion:
!
!    If input values are specified, this resets the default parameters.
!    Otherwise, the output will be the current defaults.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 April 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real t0_in: the initial time.
!
!    real y0_in: the initial condition at time T0.
!
!    real tstop_in: the final time.
!
!  Output:
!
!    real t0_out: the initial time.
!
!    real y0_out: the initial condition at time T0.
!
!    real tstop_out: the final time.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), save     :: t0_default = 0.0D+00
  real ( kind = rk ), optional :: t0_in
  real ( kind = rk ), optional :: t0_out
  real ( kind = rk ), save     :: tstop_default = 200.0D+00
  real ( kind = rk ), optional :: tstop_in
  real ( kind = rk ), optional :: tstop_out
  real ( kind = rk ), save     :: y0_default = 0.01D+00
  real ( kind = rk ), optional :: y0_in
  real ( kind = rk ), optional :: y0_out
!
!  New values, if supplied on input, overwrite the current values.
!
  if ( present ( t0_in ) ) then
    t0_default = t0_in;
  end if

  if ( present ( tstop_in ) ) then
    tstop_default = tstop_in;
  end if

  if ( present ( y0_in ) ) then
    y0_default = y0_in;
  end if
!
!  The current values are copied to the output.
!
  if ( present ( t0_out ) ) then
    t0_out = t0_default;
  end if

  if ( present ( tstop_out ) ) then
    tstop_out = tstop_default;
  end if

  if ( present ( y0_out ) ) then
    y0_out = y0_default;
  end if

  return
end
function lambert_w ( x, nb, l )

!*****************************************************************************80
!
!! lambert_w() approximates the W function.
!
!  Discussion:
!
!    The call will fail if the input value X is out of range.
!    The range requirement for the upper branch is:
!      -exp(-1) <= X.
!    The range requirement for the lower branch is:
!      -exp(-1) < X < 0.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 June 2014
!
!  Author:
!
!    Original FORTRAN77 version by Andrew Barry, S. J. Barry, 
!    Patricia Culligan-Hensley.
!    This version by John Burkardt.
!
!  Reference:
!
!    Andrew Barry, S. J. Barry, Patricia Culligan-Hensley,
!    Algorithm 743: WAPR - A Fortran routine for calculating real 
!    values of the W-function,
!    ACM Transactions on Mathematical Software,
!    Volume 21, Number 2, June 1995, pages 172-181.
!
!  Input:
!
!    real ( kind = rk ) X, the argument.
!
!    integer NB, indicates the desired branch.
!    * 0, the upper branch;
!    * nonzero, the lower branch.
!
!    integer L, indicates the interpretation of X.
!    * 1, X is actually the offset from -(exp-1), so compute W(X-exp(-1)).
!    * not 1, X is the argument; compute W(X);
!
!  Output:
!
!    real ( kind = rk ) lambert_w: the approximate value of W(X).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) an2
  real ( kind = rk ) an3
  real ( kind = rk ) an4
  real ( kind = rk ) an5
  real ( kind = rk ) an6
  real ( kind = rk ) c13
  real ( kind = rk ) c23
  real ( kind = rk ) d12
  real ( kind = rk ) delx
  real ( kind = rk ) em
  real ( kind = rk ) em2
  real ( kind = rk ) em9
  real ( kind = rk ) eta
  integer i
  integer init
  integer l
  real ( kind = rk ) lambert_w
  integer nb
  integer nbits
  integer nerror
  integer niter
  real ( kind = rk ) reta
  real ( kind = rk ) s2
  real ( kind = rk ) s21
  real ( kind = rk ) s22
  real ( kind = rk ) s23
  real ( kind = rk ) t
  real ( kind = rk ) tb
  real ( kind = rk ) tb2
  real ( kind = rk ) temp
  real ( kind = rk ) temp2
  real ( kind = rk ) ts
  real ( kind = rk ) x
  real ( kind = rk ) x0
  real ( kind = rk ) x1
  real ( kind = rk ) xx
  real ( kind = rk ) zl
  real ( kind = rk ) zn

  save an3
  save an4
  save an5
  save an6
  save c13
  save c23
  save d12
  save em
  save em2
  save em9
  save init
  save nbits
  save niter
  save s2
  save s21
  save s22
  save s23
  save tb
  save tb2
  save x0
  save x1

  data init / 0 /
  data niter / 1 /

  lambert_w = 0.0D+00
  nerror = 0

  if ( init == 0 ) then

    init = 1

    nbits = 52

    if ( 56 <= nbits ) then
      niter = 2
    end if
!
!  Various mathematical constants.
!
    em = -exp ( -1.0D+00 )
    em9 = -exp ( -9.0D+00 )
    c13 = 1.0D+00 / 3.0D+00
    c23 = 2.0D+00 * c13
    em2 = 2.0D+00 / em
    d12 = -em2
    tb = 0.5D+00 ** nbits
    tb2 = sqrt ( tb )
    x0 = tb ** ( 1.0D+00 / 6.0D+00 ) * 0.5D+00
    x1 = ( 1.0D+00 - 17.0D+00 * tb ** ( 2.0D+00 / 7.0D+00 ) ) * em
    an3 = 8.0D+00 / 3.0D+00
    an4 = 135.0D+00 / 83.0D+00
    an5 = 166.0D+00 / 39.0D+00
    an6 = 3167.0D+00 / 3549.0D+00
    s2 = sqrt ( 2.0D+00 )
    s21 = 2.0D+00 * s2 - 3.0D+00
    s22 = 4.0D+00 - 3.0D+00 * s2
    s23 = s2 - 2.0D+00

  end if

  if ( l == 1 ) then

    delx = x

    if ( delx < 0.0D+00 ) then
      nerror = 1
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'lambert_w(): Fatal error!'
      write ( *, '(a)' ) '  The offset X is negative.'
      write ( *, '(a)' ) '  It must be nonnegative.'
      stop 1
    end if

    xx = x + em

  else

    if ( x < em ) then
      nerror = 1
      return
    else if ( x == em ) then
      lambert_w = -1.0D+00
      return
    end if

    xx = x
    delx = xx - em

  end if

  if ( nb == 0 ) then
!
!  Calculations for Wp.
!
    if ( abs ( xx ) <= x0 ) then
      lambert_w = xx / ( 1.0D+00 + xx / ( 1.0D+00 + xx &
        / ( 2.0D+00 + xx / ( 0.6D+00 + 0.34D+00 * xx ))))
      return
    else if ( xx <= x1 ) then
      reta = sqrt ( d12 * delx )
      lambert_w = reta / ( 1.0D+00 + reta / ( 3.0D+00 + reta / ( reta &
        / ( an4 + reta / ( reta * an6 + an5 ) ) + an3 ) ) ) &
        - 1.0D+00
      return
    else if ( xx <= 20.0D+00 ) then
      reta = s2 * sqrt ( 1.0D+00 - xx / em )
      an2 = 4.612634277343749D+00 * sqrt ( sqrt ( reta + &
        1.09556884765625D+00 ))
      lambert_w = reta / ( 1.0D+00 + reta / ( 3.0D+00 + ( s21 * an2 &
        + s22 ) * reta / ( s23 * ( an2 + reta )))) - 1.0D+00
    else
      zl = log ( xx )
      lambert_w = log ( xx / log ( xx &
        / zl ** exp ( -1.124491989777808D+00 / &
        ( 0.4225028202459761D+00 + zl ))))
    end if
!
!  Calculations for Wm.
!
  else

    if ( 0.0D+00 <= xx ) then
      nerror = 1
      return
    else if ( xx <= x1 ) then
      reta = sqrt ( d12 * delx )
      lambert_w = reta / ( reta / ( 3.0D+00 + reta / ( reta / ( an4 &
        + reta / ( reta * an6 - an5 ) ) - an3 ) ) - 1.0D+00 ) - 1.0D+00
      return
    else if ( xx <= em9 ) then
      zl = log ( -xx )
      t = -1.0D+00 - zl
      ts = sqrt ( t )
      lambert_w = zl - ( 2.0D+00 * ts ) / ( s2 + ( c13 - t &
        / ( 270.0D+00 + ts * 127.0471381349219D+00 )) * ts )
    else
      zl = log ( -xx )
      eta = 2.0D+00 - em2 * xx
      lambert_w = log ( xx / log ( -xx / ( ( 1.0D+00 &
        - 0.5043921323068457D+00 * ( zl + 1.0D+00 ) ) &
        * ( sqrt ( eta ) + eta / 3.0D+00 ) + 1.0D+00 )))
    end if

  end if

  do i = 1, niter
    zn = log ( xx / lambert_w ) - lambert_w
    temp = 1.0D+00 + lambert_w
    temp2 = temp + c23 * zn
    temp2 = 2.0D+00 * temp * temp2
    lambert_w = lambert_w * ( 1.0D+00 + ( zn / temp ) * ( temp2 - zn ) &
      / ( temp2 - 2.0D+00 * zn ) )
  end do

  return
end
