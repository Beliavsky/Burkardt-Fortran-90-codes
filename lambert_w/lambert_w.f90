function lambert_w ( x, nb, l )

!*****************************************************************************80
!
!! lambert_w() approximates the Lambert W function.
!
!  Discussion:
!
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
!    20 June 2023
!
!  Author:
!
!    Original FORTRAN77 version by Andrew Barry, S. J. Barry, 
!    Patricia Culligan-Hensley.
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Andrew Barry, S. J. Barry, Patricia Culligan-Hensley,
!    Algorithm 743: lambert_w - A Fortran routine for calculating real 
!    values of the W-function,
!    ACM Transactions on Mathematical Software,
!    Volume 21, Number 2, June 1995, pages 172-181.
!
!  Input:
!
!    real ( kind = rk ) x: the argument.
!
!    integer nb: indicates the desired branch.
!    * 0, the upper branch;
!    * nonzero, the lower branch.
!
!    integer l: indicates the interpretation of X.
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
  integer l
  real ( kind = rk ) lambert_w
  integer nb
  integer nbits
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

  lambert_w = - huge ( 1.0D+00 )
!
!  Various mathematical constants.
!
  em = -exp ( -1.0D+00 )
  em9 = -exp ( -9.0D+00 )
  c13 = 1.0D+00 / 3.0D+00
  c23 = 2.0D+00 * c13
  em2 = 2.0D+00 / em
  d12 = -em2
  nbits = 52
  niter = 1
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

  if ( l == 1 ) then

    delx = x
!
!  The offset X must be nonnegative.
!
    if ( delx < 0.0D+00 ) then
      return
    end if

    xx = x + em

  else

    if ( x < em ) then
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
