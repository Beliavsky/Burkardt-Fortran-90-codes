program main

!*****************************************************************************80
!
!! toms577_test() tests toms577().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 May 2018
!
!  Author:
!
!    Original FORTRAN77 version by Bille Carlson, Elaine Notis.
!    This FORTRAN90 version by John Burkardt.
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TOMS577_TEST():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  TOMS577() evaluates Carlson''s elliptic functions.'

  call rc_test ( )
  call rc_test2 ( )
  call rd_test ( )
  call rf_test ( )
  call rj_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TOMS577_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine rc_test ( )

!*****************************************************************************80
!
!! RC_TEST tests RC.
!
!  Discussion:
!
!    This driver tests the real ( kind = rk ) function subroutine for the
!    integral RC(X,Y).  The first six sets of values of X and Y are
!    extreme points of the region of valid arguments defined by the
!    machine-dependent constants LOLIM and UPLIM.  The values of LOLIM,
!    UPLIM, X, Y, and ERRTOL (see comments in subroutine) may be used on
!    most machines but provide a severe test of robustness only on the
!    ibm 360/370 series.  The seventh set tests the failure exit.  The
!    next three sets are check values: RC(0,0.25) = RC(0.0625,0.125) = PI
!    and RC(2.25,2) = LN(2).  The remaining sets show the dependence on X
!    when Y = 1.  Fixing Y entails no loss here because RC is homogeneous.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 May 2018
!
!  Author:
!
!    Original FORTRAN77 version by Bille Carlson, Elaine Notis.
!    This FORTRAN90 version by John Burkardt.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) eliptc
  real ( kind = rk ) errtol
  integer i
  integer ierr
  real ( kind = rk ) rc
  real ( kind = rk ) x(43)
  real ( kind = rk ) y(43)

  save x
  save y

  data x / &
   1.51D-78, &
   3.01D-78, &
   0.00D+00, &
   0.99D+75, &
   0.00D+00, &
   0.99D+75, &
   0.00D+00, &
   0.00D+00, &
   6.25D-02, &
   2.25D+00, &
   0.01D+00, &
   0.02D+00, &
   0.05D+00, &
   0.10D+00, &
   0.20D+00, &
   0.40D+00, &
   0.60D+00, &
   0.80D+00, &
   1.00D+00, &
   1.20D+00, &
   1.50D+00, &
   2.00D+00, &
   3.00D+00, &
   4.00D+00, &
   5.00D+00, &
   1.00D+01, &
   2.00D+01, &
   5.00D+01, &
   1.00D+02, &
   1.00D+03, &
   1.00D+04, &
   1.00D+05, &
   1.00D+06, &
   1.00D+07, &
   1.00D+08, &
   1.00D+09, &
   1.00D+10, &
   1.00D+12, &
   1.00D+15, &
   1.00D+20, &
   1.00D+30, &
   1.00D+40, &
   1.00D+50 /

  data y / &
   1.51D-78, &
   0.55D-78, &
   3.01D-78, &
   0.55D-78, &
   0.99D+75, &
   0.99D+75, &
   2.00D-78, &
   2.50D-01, &
   1.25D-01, &
   2.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00 /

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'RC_TEST'
  write ( *, '(a)' ) &
    '  RC evaluates the elementary integral RC(X,Y)'

  write ( *, '(a)' ) ''
  write ( *, '(14x,a,26x,a,25x,a)' ) 'X', 'Y', 'RC(X,Y)'
  write ( *, '(a)' ) ''

  errtol = 1.0D-3

  do i = 1, 43
    eliptc = rc ( x(i), y(i), errtol, ierr )
    if ( ierr == 0 ) then
      write ( *, '(2x,3d27.16)' ) x(i), y(i), eliptc
    else
      write ( *, '(2x,2d27.16,a)' ) x(i), y(i), '  ***Error***'
    end if
  end do

  return
end
subroutine rc_test2 ( )

!*****************************************************************************80
!
!! RC_TEST2 checks RC by examining special values.
!
!  Discussion:
!
!    This driver compares values of (LOG X)/(X-1) and ARCTAN(X)
!    calculated on one hand from the subroutine RC and on the other
!    from library LOG and ARCTAN routines.  to avoid over/underflows
!    for extreme values of X, we write (LOG X)/(X-1) = RC(Y,X/Y)/SQRT(Y),
!    where Y = (1+X)/2, and ARCTAN(X) = SQRT(X)*RC(Z,Z+X), where Z = 1/X.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 May 2018
!
!  Author:
!
!    Original FORTRAN77 version by Bille Carlson, Elaine Notis.
!    This FORTRAN90 version by John Burkardt.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) errtol
  integer i
  real ( kind = rk ) ibmarc
  real ( kind = rk ) ibmlog
  integer ierr
  integer ipower
  integer j
  integer m
  real ( kind = rk ) myarc
  real ( kind = rk ) mylog
  real ( kind = rk ) rc
  real ( kind = rk ) v
  real ( kind = rk ) w
  real ( kind = rk ) x
  real ( kind = rk ) x_vec(13)
  real ( kind = rk ) y
  real ( kind = rk ) z

  save x_vec

  data x_vec / &
   1.0D-75, &
   1.0D-15, &
   1.0D-03, &
   1.0D-01, &
   2.0D-01, &
   5.0D-01, &
   1.0D+00, &
   2.0D+00, &
   5.0D+00, &
   1.0D+01, &
   1.0D+03, &
   1.0D+15, &
   1.0D+75 /

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'RC_TEST2'
  write ( *, '(a)' ) '  Compare LOG(X)/(X-1) and ARCTAN(X) with'
  write ( *, '(a)' ) '  values based on RC.'

  write ( *, '(a)' ) ''
  write ( *, '(5x,a,15x,a,19x,a)' ) 'X', 'From LOG', 'From RC'
  write ( *, '(a)' ) ''

  errtol = 1.0D-3

  do j = 1, 10
    x = 0.2D0 * dble ( j )
    y = ( 1.0D0 + x ) / 2.0D0
    v = x / y
    mylog = rc ( y, v, errtol, ierr ) / sqrt ( y )
    if ( j == 5 ) then
      write ( *, '(d9.1,5x,a,d27.16)' ) &
    x, '**** ZERO DIVIDE *****', mylog
    else
      ibmlog = log ( x ) / ( x - 1.0D0 )
      write ( *, '(d9.1,5x,2d27.16)' ) x, ibmlog, mylog
    end if
  end do

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Extreme values of X'
  write ( *, '(a)' ) ''
  write ( *, '(5x,a,15x,a,19x,a)' ) 'X', 'From LOG', 'From RC'
  write ( *, '(a)' ) ''

  do i = 1, 16
    ipower = - 75 + 10 * ( i - 1 )
    x = 10.D0 ** ipower
    y = ( 1.0D0 + x ) / 2.D0
    v = x / y
    mylog = rc ( y, v, errtol, ierr ) / sqrt ( y )
    ibmlog = log ( x ) / ( x - 1.0D0 )
    write ( *, '(d9.1,2d27.16)' ) x, ibmlog, mylog
  end do

  write ( *, '(a)' ) ''
  write ( *, '(5xa,14x,a,17x,a)' ) 'X','From ARCTAN', 'From RC'
  write ( *, '(a)' ) ''

  do m = 1, 13
    x = x_vec(m)
    z = 1.0D0 / x
    w = z + x
    myarc = sqrt ( x ) * rc ( z, w, errtol, ierr )
    ibmarc = atan ( x )
    write ( *, '(d9.1,2d27.16)' ) x, ibmarc, myarc
  end do

  return
end
subroutine rd_test ( )

!*****************************************************************************80
!
!! RD_TEST tests RD.
!
!  Discussion:
!
!    This driver tests the real ( kind = rk ) function subroutine for the
!    integral RD(X,Y,Z), which is symmetric in X and Y.  The first
!    twelve sets of values of X, Y, Z are extreme points of the region of
!    valid arguments defined by the machine-dependent constants LOLIM
!    and UPLIM.  The values of LOLIM, UPLIM, X, Y, Z, and ERRTOL (see
!    comments in subroutine) may be used on most machines but provide a
!    severe test of robustness only on the ibm 360/370 series.  The
!    thirteenth set tests the failure exit.  The fourteenth set is a
!    check value: RD(0,2,1) = 3B = 3(PI)/4A, where A and B are the
!    lemniscate constants.  The remaining sets show the dependence
!    on Z when Y = 1 (no loss of generality because of homogeneity)
!    and X = 0.5 (midway between the complete case X = 0 and the
!    degenerate case X = Y).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 May 2018
!
!  Author:
!
!    Original FORTRAN77 version by Bille Carlson, Elaine Notis.
!    This FORTRAN90 version by John Burkardt.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) eliptc
  real ( kind = rk ) errtol
  integer i
  integer ierr
  real ( kind = rk ) rd
  real ( kind = rk ) x(27)
  real ( kind = rk ) y(27)
  real ( kind = rk ) z(27)

  save x
  save y
  save z

  data x / &
   0.00D+00, &
   0.55D-78, &
   0.00D+00, &
   0.55D-78, &
   0.00D+00, &
   0.55D-78, &
   0.00D+00, &
   0.55D-78, &
   3.01D-51, &
   3.01D-51, &
   0.99D+48, &
   0.99D+48, &
   0.00D+00, &
   0.00D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00 /

  data y / &
   6.01D-51, &
   6.01D-51, &
   6.01D-51, &
   6.01D-51, &
   0.99D+48, &
   0.99D+48, &
   0.99D+48, &
   0.99D+48, &
   3.01D-51, &
   3.01D-51, &
   0.99D+48, &
   0.99D+48, &
   3.01D-51, &
   2.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00 /

  data z / &
   6.01D-51, &
   6.01D-51, &
   0.99D+48, &
   0.99D+48, &
   6.01D-51, &
   6.01D-51, &
   0.99D+48, &
   0.99D+48, &
   6.01D-51, &
   0.99D+48, &
   6.01D-51, &
   0.99D+48, &
   1.00D+00, &
   1.00D+00, &
   1.00D-10, &
   1.00D-05, &
   1.00D-02, &
   1.00D-01, &
   2.00D-01, &
   5.00D-01, &
   1.00D+00, &
   2.00D+00, &
   5.00D+00, &
   1.00D+01, &
   1.00D+02, &
   1.00D+05, &
   1.00D+10 /

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'RD_TEST'
  write ( *, '(a)' ) '  RD evaluates the Carlson elliptic integral'
  write ( *, '(a)' ) '  of the second kind, RD(X,Y,Z)'
  write ( *, '(a)' ) ''
  write ( *, '(15x,a,26x,a,26x,a,25x,a)' ) &
    'X', 'Y', 'Z', 'RD(X,Y,Z)'
  write ( *, '(a)' ) ''

  errtol = 1.0d-3

  do i = 1, 27
    eliptc = rd ( x(i), y(i), z(i), errtol, ierr )
    if (ierr == 0 ) then
      write ( *, '(4d27.16)' ) x(i), y(i), z(i), eliptc
    else
      write ( *, '(3d27.16,a)' ) x(i), y(i), z(i), '  ***Error***'
    end if
  end do

  return
end
subroutine rf_test ( )

!*****************************************************************************80
!
!! RF_TEST tests RF.
!
!  Discussion:
!
!    This driver tests the real ( kind = rk ) function subroutine for the
!    integral RF(X,Y,Z), which is symmetric in X, Y, Z.  The first nine
!    sets of values of X, Y, Z are extreme points of the region of valid
!    arguments defined by the machine-dependent constants LOLIM and
!    UPLIM.  The values of LOLIM, UPLIM, X, Y, Z, and ERRTOL (see
!    comments in subroutine) may be used on most machines but provide a
!    severe test of robustness only on the ibm 360/370 series.  The
!    tenth set tests the failure exit.  The eleventh set is a check
!    value: RF(0,1,2) = A, where A is the first lemniscate constant.
!    The remaining sets show the dependence on Z when Y = 1 (no loss of
!    generality because of homogeneity) and X = 0.5 (midway between the
!    complete case X = 0 and the degenerate case X = Y).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 May 2018
!
!  Author:
!
!    Original FORTRAN77 version by Bille Carlson, Elaine Notis.
!    This FORTRAN90 version by John Burkardt.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) eliptc
  real ( kind = rk ) errtol
  integer i
  integer ierr
  real ( kind = rk ) rf
  real ( kind = rk ) x(55)
  real ( kind = rk ) y(55)
  real ( kind = rk ) z(55)

  save x
  save y
  save z

  data x / &
   1.51D-78, &
   1.51D-78, &
   0.00D+00, &
   0.00D+00, &
   0.00D+00, &
   0.99D+75, &
   0.55D-78, &
   0.55D-78, &
   0.55D-78, &
   0.00D+00, &
   0.00D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00 /

  data y / &
   1.51D-78, &
   1.51D-78, &
   3.01D-78, &
   3.01D-78, &
   0.99D+75, &
   0.99D+75, &
   3.01D-78, &
   3.01D-78, &
   0.99D+75, &
   2.00D-78, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00 /

  data z / &
   1.51D-78, &
   0.99D+75, &
   3.01D-78, &
   0.99D+75, &
   0.99D+75, &
   0.99D+75, &
   3.01D-78, &
   0.99D+75, &
   0.99D+75, &
   1.00D+00, &
   2.00D+00, &
   1.00D+00, &
   1.10D+00, &
   1.20D+00, &
   1.30D+00, &
   1.40D+00, &
   1.50D+00, &
   1.60D+00, &
   1.70D+00, &
   1.80D+00, &
   1.90D+00, &
   2.00D+00, &
   2.20D+00, &
   2.40D+00, &
   2.60D+00, &
   2.80D+00, &
   3.00D+00, &
   3.50D+00, &
   4.00D+00, &
   4.50D+00, &
   5.00D+00, &
   6.00D+00, &
   7.00D+00, &
   8.00D+00, &
   9.00D+00, &
   1.00D+01, &
   2.00D+01, &
   3.00D+01, &
   4.00D+01, &
   5.00D+01, &
   1.00D+02, &
   2.00D+02, &
   5.00D+02, &
   1.00D+03, &
   1.00D+04, &
   1.00D+05, &
   1.00D+06, &
   1.00D+08, &
   1.00D+10, &
   1.00D+12, &
   1.00D+15, &
   1.00D+20, &
   1.00D+30, &
   1.00D+40, &
   1.00D+50 /

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'RF_TEST'
  write ( *, '(a)' ) '  RF evaluates the Carlson elliptic integral'
  write ( *, '(a)' ) '  of the first kind, RF(X,Y,Z)'
  write ( *, '(a)' ) ''
  write ( *, '(15x,a,26x,a,26x,a,25x,a)' ) &
    'X', 'Y', 'Z', 'RF(X,Y,Z)'
  write ( *, '(a)' ) ''

  errtol = 1.0D-3

  do i = 1, 55
    eliptc = rf ( x(i), y(i), z(i), errtol, ierr )
    if (ierr == 0 ) then
      write ( *, '(4d27.16)' ) x(i), y(i), z(i), eliptc
    else
      write ( *, '(3d27.16,a)' ) x(i), y(i), z(i), '  ***Error***'
    end if
  end do

  return
end
subroutine rj_test ( )

!*****************************************************************************80
!
!! RJ_TEST tests RJ.
!
!  Discussion:
!
!    This driver tests the real ( kind = rk ) function subroutine for the
!    integral Rj(X,Y,Z,P), which is symmetric in X, Y, Z.  The first
!    twenty sets of values of X, Y, Z, P are extreme points of the region
!    of valid arguments defined by the machine-dependent constants
!    LOLIM and UPLIM.  The values of LOLIM, UPLIM, X, Y, Z, P, and
!    ERRTOL (see comments in subroutine) may be used on most machines
!    but provide a severe test of robustness only on the ibm 360/370
!    series.  The twenty-first set tests the failure exit.  The twenty-
!    second set is a check value:
!      RJ(2,3,4,5) = 0.1429757966715675383323308.
!    The remaining sets show the dependence on Z and P
!    when Y = 1 (no loss of generality because of homogeneity) and
!    X = 0.5 (midway between the complete case x = 0 and the degenerate
!    case X = Y).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 May 2018
!
!  Author:
!
!    Original FORTRAN77 version by Bille Carlson, Elaine Notis.
!    This FORTRAN90 version by John Burkardt.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) eliptc
  real ( kind = rk ) errtol
  integer i
  integer ierr
  real ( kind = rk ) p(42)
  real ( kind = rk ) rj
  real ( kind = rk ) x(42)
  real ( kind = rk ) y(42)
  real ( kind = rk ) z(42)

  save p
  save x
  save y
  save z

  data p / &
   2.01D-26, &
   2.01D-26, &
   2.01D-26, &
   2.01D-26, &
   2.01D-26, &
   2.01D-26, &
   2.01D-26, &
   2.01D-26, &
   2.01D-26, &
   2.01D-26, &
   2.99D+24, &
   2.99D+24, &
   2.99D+24, &
   2.99D+24, &
   2.99D+24, &
   2.99D+24, &
   2.99D+24, &
   2.99D+24, &
   2.99D+24, &
   2.99D+24, &
   1.00D+00, &
   5.00D+00, &
   0.25D+00, &
   0.75D+00, &
   1.00D+00, &
   2.00D+00, &
   0.25D+00, &
   0.75D+00, &
   1.50D+00, &
   4.00D+00, &
   0.25D+00, &
   0.75D+00, &
   3.00D+00, &
   1.00D+01, &
   0.25D+00, &
   0.75D+00, &
   5.00D+00, &
   2.00D+01, &
   0.25D+00, &
   0.75D+00, &
   5.00D+01, &
   2.00D+02 /

  data x / &
   1.01D-26, &
   1.01D-26, &
   0.00D+00, &
   0.00D+00, &
   0.00D+00, &
   2.99D+24, &
   0.55D-78, &
   0.55D-78, &
   0.55D-78, &
   2.01D-26, &
   1.01D-26, &
   1.01D-26, &
   0.00D+00, &
   0.00D+00, &
   0.00D+00, &
   2.99D+24, &
   0.55D-78, &
   0.55D-78, &
   0.55D-78, &
   2.01D-26, &
   0.00D+00, &
   2.00D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00, &
   0.50D+00 /

  data y / &
   1.01D-26, &
   1.01D-26, &
   2.01D-26, &
   2.01D-26, &
   2.99D+24, &
   2.99D+24, &
   2.01D-26, &
   2.01D-26, &
   2.99D+24, &
   2.01D-26, &
   1.01D-26, &
   1.01D-26, &
   2.01D-26, &
   2.01D-26, &
   2.99D+24, &
   2.99D+24, &
   2.01D-26, &
   2.01D-26, &
   2.99D+24, &
   2.01D-26, &
   1.90D-26, &
   3.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00 /

  data z / &
   1.01D-26, &
   2.99D+24, &
   2.01D-26, &
   2.99D+24, &
   2.99D+24, &
   2.99D+24, &
   2.01D-26, &
   2.99D+24, &
   2.99D+24, &
   2.01D-26, &
   1.01D-26, &
   2.99D+24, &
   2.01D-26, &
   2.99D+24, &
   2.99D+24, &
   2.99D+24, &
   2.01D-26, &
   2.99D+24, &
   2.99D+24, &
   2.01D-26, &
   1.90D-26, &
   4.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   1.00D+00, &
   2.00D+00, &
   2.00D+00, &
   2.00D+00, &
   2.00D+00, &
   5.00D+00, &
   5.00D+00, &
   5.00D+00, &
   5.00D+00, &
   1.00D+01, &
   1.00D+01, &
   1.00D+01, &
   1.00D+01, &
   1.00D+02, &
   1.00D+02, &
   1.00D+02, &
   1.00D+02 /

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'RJ_TEST'
  write ( *, '(a)' ) '  RJ evaluates the Carlson elliptic integral'
  write ( *, '(a)' ) '  of the third kind, RJ(X,Y,Z,P)'
  write ( *, '(a)' ) ''
  write ( *, '(15x,a,26x,a,26x,a,26x,a,25x,a)' ) &
    'X', 'Y', 'Z', 'P', 'RJ(X,Y,Z,P)'
  write ( *, '(a)' ) ''

  errtol = 1.0D-3

  do i = 1, 42
    eliptc = rj ( x(i), y(i), z(i), p(i), errtol, ierr )
    if (ierr == 0 ) then
      write ( *, '(5d27.16)' ) x(i), y(i), z(i), p(i), eliptc
    else
      write ( *, '(4d27.16,a)' ) x(i), y(i), z(i), p(i), '  ***Error***'
    end if
  end do

  return
end

