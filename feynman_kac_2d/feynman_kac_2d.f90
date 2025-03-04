program main

!*****************************************************************************80
!
!! feynman_kac_2d() uses the Feynman-Kac method on a 2D PDE.
!
!  Discussion:
!
!    This program is derived from section 2.5, exercise 2.2 of Petersen 
!    and Arbenz.
!
!    The problem is to determine the solution U(X,Y) of the following 
!    partial differential equation:
!
!      (1/2) Laplacian U - V(X,Y) * U = 0,
!
!    inside the elliptic domain D:
! 
!      D = { (X,Y) | (X/A)^2+(Y/B)^2 <= 1 }
!   
!    with the boundary condition U(boundary(D)) = 1.
!
!    V(X,Y) is the potential function:
!
!      V = 2 * ( (X/A^2)^2 + (Y/B^2)^2 ) + 1/A^2 + 1/B^2.
!
!    The analytic solution of this problem is already known:
!
!      U(X,Y) = exp ( (X/A)^2 + (Y/B)^2 - 1 ).
!
!    Our method is via the Feynman-Kac Formula.
!
!    The idea is to start from any (x,y) in D, and
!    compute (x+Wx(t),y+Wy(t)) where 2D Brownian motion
!    (Wx,Wy) is updated each step by sqrt(h)*(z1,z2),
!    each z1,z2 are independent approximately Gaussian 
!    random variables with zero mean and variance 1. 
!
!    Each (x1(t),x2(t) ) is advanced until (x1,x2 ) exits the domain D.  
!
!    Upon its first exit from D, the sample path (x1,x2) is stopped and a 
!    new sample path at (x,y) is started until N such paths are completed.
! 
!    The Feynman-Kac formula gives the solution here as
!
!      U(X,Y) = (1/N) sum(1 <= I <= N) Y(tau_i),
!
!    where
!
!      Y(tau) = exp( -int(s=0..tau) v(x1(s),x2(s)) ds),
!
!    and tau = first exit time for path (x1,x2). 
!
!    The integration procedure is a second order weak accurate method:
!
!      X(t+h)  = [ x1(t) + sqrt ( h ) * z1 ]
!                [ x2(t) + sqrt ( h ) * z2 ]
!
!    Here Z1 and Z2 are approximately normal univariate Gaussians. 
!
!    An Euler predictor approximates Y at the end of the step
!
!      Y_e     = (1 - h*v(X(t)) * Y(t), 
!
!    A trapezoidal rule completes the step:
!
!      Y(t+h)  = Y(t) - (h/2)*[v(X(t+h))*Y_e + v(X(t))*Y(t)].
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 May 2012
!
!  Author:
!
!    Original C 3D version by Wesley Petersen.
!    FORTRAN90 2D version by John Burkardt.
!
!  Reference:
!
!    Peter Arbenz, Wesley Petersen,
!    Introduction to Parallel Computing:
!    A Practical Guide with Examples in C,
!    Oxford, 2004,
!    ISBN: 0-19-851577-4,
!    LC: QA76.59.P47.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) :: a = 2.0D+00
  real ( kind = rk ) :: b = 1.0D+00
  real ( kind = rk ) chk
  integer :: dim = 2
  real ( kind = rk ) dx
  real ( kind = rk ) dy
  real ( kind = rk ) err
  real ( kind = rk ) :: h = 0.0001D+00
  integer i
  integer i4_ceiling
  integer j
  integer k
  integer :: n = 10000
  integer n_inside
  integer ni
  integer nj
  real ( kind = rk ) rth
  integer steps
  integer steps_ave
  real ( kind = rk ) us
  real ( kind = rk ) ut
  real ( kind = rk ) vh
  real ( kind = rk ) vs
  real ( kind = rk ) x
  real ( kind = rk ) x1
  real ( kind = rk ) x2
  real ( kind = rk ) y
  real ( kind = rk ) w
  real ( kind = rk ) w_exact
  real ( kind = rk ) we
  real ( kind = rk ) wt

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEYNMAN_KAC_2D():'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Program parameters:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The calculation takes place inside a 2D ellipse.'
  write ( *, '(a)' ) '  A rectangular grid of points will be defined.'
  write ( *, '(a)' ) '  The solution will be estimated for those grid points'
  write ( *, '(a)' ) '  that lie inside the ellipse.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) &
    '  Each solution will be estimated by computing ', n, &
    ' trajectories'
  write ( *, '(a)' ) '  from the point to the boundary.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    (X/A)^2 + (Y/B)^2 = 1'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The ellipse parameters A and B are:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '    A = ', a
  write ( *, '(a,g14.6)' ) '    B = ', b
  write ( *, '(a,g14.6)' ) '  Stepsize H = ', h
!
!  RTH is the scaled stepsize.
!
  rth = sqrt ( real ( dim, kind = rk ) * h )
!
!  Choose the spacing so we have about 10 points in the shorter direction.
!
  if ( a < b ) then
    ni = 11
    nj = 1 + i4_ceiling ( b / a ) * ( ni - 1 )
  else
    nj = 11
    ni = 1 + i4_ceiling ( a / b ) * ( nj - 1 )
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4,a)' ) '  X coordinate marked by ', ni, ' points'
  write ( *, '(a,i4,a)' ) '  Y coordinate marked by ', nj, ' points'

  err = 0.0D+00
  n_inside = 0
!
!  Loop over the grid points.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '     X            Y             W Approx' // &
    '      W Exact          Error      Ave Steps'
  write ( *, '(a)' ) ' '

  do j = 1, nj

    x = ( real ( nj - j,     kind = rk ) * ( - a ) &
        + real (      j - 1, kind = rk ) *     a ) &
        / real ( nj     - 1, kind = rk )

    do i = 1, ni

      y = ( real ( ni - i,     kind = rk ) * ( - b ) &
          + real (      i - 1, kind = rk ) *     b ) & 
          / real ( ni     - 1, kind = rk )

      chk = ( x / a )**2 + ( y / b )**2

      if ( 1.0D+00 < chk ) then
        w_exact = 1.0D+00
        wt = 1.0D+00
        steps_ave = 0
        write ( *, &
          '(5(2x,g12.4),2x,i8)' ) &
          x, y, wt, w_exact, abs ( w_exact - wt ), steps_ave
        cycle
      end if

      n_inside = n_inside + 1
!
!  Compute the exact solution at this point (x,y,z).
!
      w_exact = exp ( ( x / a )**2 + ( y / b )**2 - 1.0D+00 )
!
!  Now try to estimate the solution at this point.
!
      wt = 0.0D+00
      steps = 0

      do k = 1, n

        x1 = x
        x2 = y
! 
!  W = exp(-int(s=0..t) v(X)ds) 
!
        w = 1.0D+00
!
!  CHK is < 1.0 while the point is inside the ellipsoid.
!
        chk = 0.0D+00

        do while ( chk < 1.0D+00 )
!
!  Determine DX, DY.
!
          call random_number ( harvest = ut )
!
!  How is 1/2 or 1/3 chosen here?
!
          if ( ut < 1.0D+00 / 2.0D+00 ) then
            call random_number ( harvest = us )
            us = us - 0.5D+00
            if ( us < 0.0D+00 ) then
              dx = - rth
            else
              dx = + rth
            end if
          else
            dx = 0.0D+00
          end if

          call random_number ( harvest = ut )
          if ( ut < 1.0D+00 / 2.0D+00 ) then
            call random_number ( harvest = us )
            us = us - 0.5D+00
            if ( us < 0.0D+00 ) then
              dy = - rth
            else
              dy = + rth
            end if
          else
            dy = 0.0D+00
          end if

          call potential ( a, b, x1, x2, vs )
!
!  Move to the new point.
!
          x1 = x1 + dx
          x2 = x2 + dy

          steps = steps + 1

          call potential ( a, b, x1, x2, vh )

          we = ( 1.0D+00 - h * vs ) * w
          w = w - 0.5D+00 * h * ( vh * we + vs * w ) 

          chk = ( x1 / a )**2 + ( x2 / b )**2

        end do

        wt = wt + w

      end do
!
!  WT is the average of the sum of the different trials.
!
      wt = wt / real ( n, kind = rk )
      steps_ave = steps / n
!
!  Add error in WT to the running L2 error in the solution.
!
      err = err + ( w_exact - wt )**2

      write ( *, '(5(2x,g12.4),2x,i8)' ) &
        x, y, wt, w_exact, abs ( w_exact - wt ), steps_ave

    end do

  end do
!
!  Compute the RMS error for all the points.
!
  err = sqrt ( err / real ( n_inside, kind = rk ) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  RMS absolute error in solution = ', err
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEYNMAN_KAC_2D:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '

  call timestamp ( )

  stop
end
function i4_ceiling ( r )

!*****************************************************************************80
!
!! I4_CEILING rounds an R8 "up" (towards +oo) to the next I4.
!
!  Example:
!
!    R     Value
!
!   -1.1  -1
!   -1.0  -1
!   -0.9   0
!    0.0   0
!    5.0   5
!    5.1   6
!    5.9   6
!    6.0   6
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) R, the value to be rounded up.
!
!    Output, integer I4_CEILING, the rounded value.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i4_ceiling
  real ( kind = rk ) r
  integer value

  value = int ( r )
  if ( real ( value, kind = rk ) < r ) then
    value = value + 1
  end if

  i4_ceiling = value

  return
end
subroutine potential ( a, b, x, y, v )

!*****************************************************************************80
!
!! POTENTIAL evaluates the potential function V(X,Y).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    08 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, the parameters that define the ellipse.
!
!    Input, real ( kind = rk ) X, Y, the coordinates of the point.
!
!    Output, real ( kind = rk ) V, the value of the potential function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) v
  real ( kind = rk ) x
  real ( kind = rk ) y

  v = 2.0D+00 * ( ( x / a**2 )**2 + ( y / b**2 )**2 ) &
    + 1.0D+00 / a**2 + 1.0D+00 / b**2

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
