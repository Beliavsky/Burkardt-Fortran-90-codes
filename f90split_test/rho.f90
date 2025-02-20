subroutine rho ( abserr, fxname, iprint, maxstp, maxtab, &
  relerr, xroot, xtry1, xtry2 )

!*****************************************************************************80
!
!! rho() includes an INTERFACE statement.
!
!  Discussion:
!
!    The method uses the idea of considering the related function
!
!      H(X) = 1 / F(X)
!
!    The iteration begins with two estimates for the root supplied by
!    the user.
!
!    From the most recent approximation to the root, X(K), the next
!    approximation X(K+1) is determined by:
!
!      X(K+1) = X(K) + H(X(K-R),...,X(K-1)) / H(X(K-R),...,X(K-1),X(K))
!
!    where K-R = 1 until the maximal order NTAB is reached.
!
!    Generally, the next iterate X(K+1) is the zero of a rational function
!    which passes through the previous data points.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    FM Larkin,
!    Root Finding by Divided Differences,
!    Numerische Mathematik,
!    Volume 37, pages 93-104, 1981.
!
!  Parameters:
!
!    Input, real ( kind = rk ) ABSERR, a positive absolute error tolerance.
!    If an estimate X for the root is found with ABS ( F(X) ) <= ABSERR,
!    the iteration is stopped.
!
!    Input, external FXNAME, the name of the function routine which evaluates
!    F(X).  The form of FXNAME must be similar to the following function which
!    has F(X) = ( X - 1 ) * ( X + 1 ).
!
!    function parab ( x )
!
!      real ( kind = rk ) parab
!      real ( kind = rk ) x
!
!      parab = ( x - 1.0D+00 ) * ( x + 1.0D+00 )
!
!      return
!    end
!
!    Input, integer IPRINT, a switch controlling printed output:
!    0, only print error messages.
!    nonzero, also print a table of the iterative process.
!
!    Input, integer MAXSTP, the limit on how many iterations
!    may be tried.
!
!    Input, integer MAXTAB, the limit on how high an order can be
!    used in the divided difference table.  MAXTAB must be at least 2, and
!    probably should not be too large.  Perhaps a value of 5 or 6 is reasonable,
!    20 is too large.
!
!    Input, real ( kind = rk ) RELERR, a tolerance on the size of the change
!    in the root estimates.  If a step is taken, and the change in the root
!    estimate is less than RELERR, the iteration will stop.
!
!    Output, real ( kind = rk ) XROOT, the point which the program has
!    produced as an approximate root.
!    Either ABS ( F(XROOT) ) <= ABSERR, or the maximum number of steps was
!    reached, or the current estimate of the root could not be significantly
!    improved.
!
!    Input, real ( kind = rk ) XTRY1, XTRY2, two initial approximations to
!    the root, supplied by the user, which must be distinct.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, intent(in) :: maxtab

  real ( kind = rk ), intent(in) :: abserr
  real ( kind = rk ) diftab(maxtab)
  real ( kind = rk ) froot
  real ( kind = rk ) ftemp1
  real ( kind = rk ) ftemp2
!  real ( kind = rk ), external :: fxname
  integer, intent(in) :: iprint
  integer istep
  integer, intent(in) :: maxstp
  integer ntab
  real ( kind = rk ), intent(in) :: relerr
  real ( kind = rk ) xdelt
  real ( kind = rk ) xold
  real ( kind = rk ), intent(out) :: xroot
  real ( kind = rk ) xtab(maxtab)
  real ( kind = rk ), intent(in) :: xtry1, xtry2
  real ( kind = rk ) yval
!
  interface
    function fxname(x)
      implicit none
      real ( kind = rk ) :: fxname
      real ( kind = rk ), intent(in) :: x
    end function
  end interface
!
!  Make sure XTRY1 and XTRY2 are not equal.
!
  if ( xtry1 == xtry2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'rho - Fatal error!'
    write ( *, '(a)' ) '  XTRY1 = XTRY2 on input.'
    stop 1
  end if
!
!  Make sure MAXTAB is at least 2.
!
  if ( maxtab < 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'rho - Fatal error!'
    write ( *, '(a)' ) '  MAXTAB < 2 on input!'
    stop 1
  end if

  xtab(1) = xtry1
  xtab(2) = xtry2
  ftemp1 = fxname ( xtry1 )
  ftemp2 = fxname ( xtry2 )

  if ( abs ( ftemp2 ) < abs ( ftemp1 ) ) then
    xtab(1) = xtry2
    xtab(2) = xtry1
    call r8_swap ( ftemp1, ftemp2 )
  end if
!
!  Initialize the number of steps.
!
  istep = 0
!
!  Initialize the number of data points.
!
  ntab = 2

  if ( 0 < iprint ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '   Step  NTAB    XROOT        F(XROOT)      XDELT'
    write ( *, '(a)' ) ' '
  end if
!
!  Initialize the divided difference table data.
!
  diftab(1) = 1.0D+00 / ftemp1
  diftab(2) = 1.0D+00 / ftemp2

  call data_to_dif ( ntab, xtab, diftab, diftab )
!
!  Initialize values used in the iteration.
!
  xroot = xtry1
  froot = ftemp1
  xdelt = xtry1 - xtry2
!
!  Does the starting data already satisfy the function norm
!  error tolerance ABSERR, or the interval norm error tolerance
!  RELERR?
!
  do

    if ( 0 < iprint ) then
      write ( *, '(3x,i4,4x,i2, 3g14.6)' ) istep, ntab, xroot, froot, xdelt
    end if

    if ( abs ( froot ) <= abserr ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'rho - Absolute convergence,'
      write ( *, '(a)' ) '  The function value meets the error tolerance.'
      exit
    end if

    if ( abs ( xdelt ) <= relerr ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'rho - Relative convergence.'
      write ( *, '(a)' ) '  The stepsize meets the error tolerance.'
      exit
    end if
!
!  Check the number of steps taken.
!
    if ( maxstp <= istep ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'rho - Nonconvergence!'
      write ( *, '(a)' ) '  The maximum number of steps was taken.'
      exit
    end if
!
!  Generate the next point, XVAL.
!
    xold = xroot
    istep = istep + 1
!
!  Algorithm breakdown: The divisor DIFTAB(NTAB) is zero.
!
    if ( diftab(ntab) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'rho - Fatal error!'
      write ( *, '(a,i6)' ) '  Algorithm using differences of order ', ntab
      write ( *, '(a)' ) '  A zero-divisor was computed.'
      write ( *, '(a)' ) '  The algorithm has broken down.'
      write ( *, '(a)' ) '  Examine the results.  They may be useful.'
      write ( *, '(a)' ) '  Perhaps a lower value of MAXTAB would help.'
      stop 1
    end if

    xroot = xtab(ntab) + ( diftab(ntab-1) / diftab(ntab) )
    xdelt = xroot - xold
    froot = fxname ( xroot )

    if ( abs ( froot ) <= abserr ) then
      cycle
    end if

    yval = 1.0D+00 / froot
!
!  If we are now using MAXTAB points, we have to remove an old
!  one before adding the new one.
!
    if ( maxtab <= ntab ) then
      ntab = ntab - 1
    end if

    call dif_append ( ntab, xtab, diftab, xroot, yval, ntab, xtab, diftab )

  end do

  return
end
