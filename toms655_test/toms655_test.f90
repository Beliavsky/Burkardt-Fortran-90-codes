program main

!*****************************************************************************80
!
!! toms655_test() tests toms655().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 February 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) alpha
  real ( kind = rk ) b
  real ( kind = rk ) beta
  integer kinda
  integer nt

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'toms655_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test TOMS655().'

  call test01 ( )
  call test02 ( )
  call ceiqfs_test ( )
  call ceiqf_test ( )
  call cliqfs_test ( )
  call cliqf_test ( )
  call cegqf_test ( )
  call cegqfs_test ( )
  call cgqfs_test ( )
!
!  Compute 15 points of an example of each rule, with default A, B.
!
  do kinda = 1, 9
    nt = 15
    if ( kinda == 8 ) then
      alpha = 1.0D+00
      beta = - alpha - 2 * nt - 2
    else
      alpha = 0.0D+00
      beta = 0.0D+00
    end if
    call test10 ( nt, kinda, alpha, beta )
  end do
!
!  Compute 15 points of an example of each rule using nondefault A, B.
!
  do kinda = 1, 9

    nt = 15

    if ( kinda == 1 ) then
      alpha = 0.0D+00
      beta = 0.0D+00
      a = 0.0D+00
      b = 1.0D+00
    else if ( kinda == 2 ) then
      alpha = 0.0D+00
      beta = 0.0D+00
      a = 0.0D+00
      b = 1.0D+00
    else if ( kinda == 3 ) then
      alpha = 1.0D+00
      beta = 0.0D+00
      a = 0.0D+00
      b = 1.0D+00
    else if ( kinda == 4 ) then
      alpha = 1.5D+00
      beta = 0.5D+00
      a = 0.0D+00
      b = 1.0D+00
    else if ( kinda == 5 ) then
      alpha = 1.0D+00
      beta = 0.0D+00
      a = 1.0D+00
      b = 1.0D+00
    else if ( kinda == 6 ) then
      alpha = 1.0D+00
      beta = 0.0D+00
      a = 0.0D+00
      b = 0.5D+00
    else if ( kinda == 7 ) then
      alpha = 1.0D+00
      beta = 0.0D+00
      a = 0.0D+00
      b = 1.0D+00
    else if ( kinda == 8 ) then
      alpha = 1.0D+00
      beta = - alpha - 2 * nt - 2
      a = 0.0D+00
      b = 1.0D+00
    else if ( kinda == 9 ) then
      alpha = 0.0D+00
      beta = 0.0D+00
      a = 0.0D+00
      b = 1.0D+00
    end if

    call test11 ( nt, kinda, alpha, beta, a, b )

  end do

  call wm_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS655_TEST()'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests CIQFS.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 January 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) alpha
  real ( kind = rk ) beta
  integer i
  integer key
  integer kinda
  integer lu
  integer, allocatable :: mlt(:)
  integer, allocatable :: ndx(:)
  integer nt
  integer nwts
  real ( kind = rk ), parameter :: pi = 3.14159265358979323846264338327950D+00
  real ( kind = rk ), allocatable :: t(:)
  real ( kind = rk ), allocatable :: wts(:)

  write ( *, '(a)' ) '  ----------------------------------------'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Test CIQFS.'
!
!  Number of knots.
!
  nt = 5
!
!  Set the knots in the default interval [-1,+1].
!
  allocate ( t(1:nt) )

  do i = 1, nt
    t(i) = cos ( real ( 2 * i - 1, kind = rk ) * pi / real ( 2 * nt, kind = rk ) )
  end do
!
!  Set the knot multiplicities.
!
  allocate ( mlt(1:nt) )
  mlt(1:nt) = 2
!
!  Set the size of the weights array.
!
  nwts = sum ( mlt(1:nt) )
!
!  Because KEY = 1, NDX will be set up for us.
!
  allocate ( ndx(1:nt) )
!
!  KEY = 1 indicates that the WTS array should hold the weights
!  in the usual order.
!
  key = 1
!
!  Request Legendre weight function.
!
  kinda = 1
!
!  ALPHA, BETA not used in Legendre weight function but set anyway.
!
  alpha = 0.0D+00
  beta  = 0.0D+00
!
!  LU controls printing.
!  A positive value requests that we compute and print weights, and
!  conduct a moments check.
!
  lu = 6
!
!  This call returns the WTS array.
!
  allocate ( wts(1:nwts) )

  call ciqfs ( nt, t, mlt, nwts, ndx, key, kinda, alpha, beta, lu, wts )

  deallocate ( mlt )
  deallocate ( ndx )
  deallocate ( t )
  deallocate ( wts )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests CIQFS.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 February 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) alpha
  real ( kind = rk ) b
  real ( kind = rk ) beta
  integer key
  integer kinda
  integer lu
  integer, allocatable :: mlt(:)
  integer, allocatable :: ndx(:)
  integer nt
  integer nwts
  real ( kind = rk ), allocatable :: t(:)
  real ( kind = rk ), allocatable :: wts(:)

  write ( *, '(a)' ) '  ----------------------------------------'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Test CIQF, CIQFS, CGQF and CGQFS'
  write ( *, '(a)' ) '  with all classical weight functions.'
!
!  Try all weight functions.
!
  do kinda = 1, 9
!
!  Number of knots.
!
    nt = 5
!
!  Set parameters ALPHA and BETA.
!
    alpha = 0.5D+00
    if ( kinda /= 8 ) then
      beta  = 2.0D+00
    else
      beta = - 16.0D+00
    end if
!
!  Set A and B.
!
    a = - 0.5D+00
    b = 2.0D+00
    lu = 6
!
!  Have CGQF compute the knots and weights.
!
    allocate ( t(1:nt) )
    allocate ( wts(1:nt) )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Knots and weights of Gauss quadrature formula'
    write ( *, '(a)' ) '  computed by CGQF.'

    call cgqf ( nt, kinda, alpha, beta, a, b, lu, t, wts )
!
!  Now compute the weights for the same knots by CIQF.
!
!  Set the knot multiplicities.
!
    allocate ( mlt(1:nt) )
    mlt(1:nt) = 2
!
!  Set the size of the weights array.
!
    nwts = sum ( mlt(1:nt) )
!
!  We need to deallocate and reallocate WTS because it is now of
!  dimension NWTS rather than NT.
!
    deallocate ( wts )
    allocate ( wts(1:nwts) )
!
!  Because KEY = 1, NDX will be set up for us.
!
    allocate ( ndx(1:nt) )
!
!  KEY = 1 indicates that the WTS array should hold the weights
!  in the usual order.
!
    key = 1
!
!  LU controls printing.
!  A positive value requests that we compute and print weights, and
!  conduct a moments check.
!
    lu = 6

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Weights of Gauss quadrature formula computed from the'
    write ( *, '(a)' ) '  knots by CIQF.'

    call ciqf ( nt, t, mlt, nwts, ndx, key, kinda, alpha, beta, a, b, lu, wts )

    deallocate ( mlt )
    deallocate ( ndx )
    deallocate ( t )
    deallocate ( wts )

  end do

  return
end
subroutine ceiqfs_test ( )

!*****************************************************************************80
!
!! CEIQFS_TEST tests CEIQFS.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 January 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) alpha
  real ( kind = rk ) beta
  real ( kind = rk ), external :: f
  integer i
  integer kinda
  integer, allocatable :: mlt(:)
  integer nt
  real ( kind = rk ), parameter :: pi = 3.14159265358979323846264338327950D+00
  real ( kind = rk ) qfsum
  real ( kind = rk ) qfsx
  real ( kind = rk ), allocatable :: t(:)

  write ( *, '(a)' ) '  ----------------------------------------'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CEIQFS_TEST'
  write ( *, '(a)' ) '  Test CEIQFS.'
!
!  Number of knots.
!
  nt = 5
!
!  Set the knots in the default interval [-1,+1].
!
  allocate ( t(1:nt) )

  do i = 1, nt
    t(i) = cos ( real ( 2 * i - 1, kind = rk ) * pi / real ( 2 * nt, kind = rk ) )
  end do
!
!  Set the knot multiplicities.
!
  allocate ( mlt(1:nt) )
  mlt(1:nt) = 2
!
!  Set KINDA to the Legendre weight function.
!
  kinda = 1
!
!  ALPHA, BETA not used in Legendre weight function but set anyway.
!
  alpha = 0.0D+00
  beta  = 0.0D+00
!
!  Call CEIQFS to set up the quadrature formula and evaluate it on F.
!
  call ceiqfs ( nt, t, mlt, kinda, alpha, beta, f, qfsum )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integral of sin(x) on -1, 1 by Fejer type rule'
  write ( *, '(a,i4,a,i4)' ) '  with ', nt, ' points of multiplicity 2.'
  write ( *, '(a,g24.16)' ) '  Quadrature formula:', qfsum

  qfsx = cos ( - 1.0D+00 ) - cos ( 1.0D+00 )
  write ( *, '(a,g24.16)' ) '  Exact value       :', qfsx
  write ( *, '(a,g24.16)' ) '  Error             :', abs ( qfsum - qfsx )

  deallocate ( mlt )
  deallocate ( t )

  return
end
subroutine ceiqf_test ( )

!*****************************************************************************80
!
!! CEIQF_TEST tests CEIQF.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 January 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) alpha
  real ( kind = rk ) b
  real ( kind = rk ) beta
  real ( kind = rk ), external :: f
  integer i
  integer kinda
  integer, allocatable :: mlt(:)
  integer nt
  real ( kind = rk ), parameter :: pi = 3.14159265358979323846264338327950D+00
  real ( kind = rk ) qfsum
  real ( kind = rk ) qfsx
  real ( kind = rk ), allocatable :: t(:)

  write ( *, '(a)' ) '  ----------------------------------------'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CEIQF_TEST'
  write ( *, '(a)' ) '  Test CEIQF.'
!
!  Number of knots.
!
  nt = 5
!
!  Set the knots in the default interval [-1,+1].
!
  allocate ( t(1:nt) )

  do i = 1, nt
    t(i) = cos ( real ( 2 * i - 1, kind = rk ) * pi / real ( 2 * nt, kind = rk ) )
  end do
!
!  Set the knot multiplicities.
!
  allocate ( mlt(1:nt) )
  mlt(1:nt) = 2
!
!  Set KINDA to the Legendre weight function.
!
  kinda = 1
!
!  ALPHA, BETA not used in Legendre weight function but set anyway.
!
  alpha = 0.0D+00
  beta  = 0.0D+00
!
!  Set nonstandard interval A, B.
!
  a = -0.5D+00
  b = 2.0D+00
!
!  Shift knots from [-1,1] to [A,B].
!
  do i = 1, nt
    t(i) = ( ( b - a ) * t(i) + ( a + b ) ) / 2.0D+00
  end do
!
!  Call CEIQF to set up the quadrature formula and evaluate it on F.
!
  call ceiqf ( nt, t, mlt, kinda, alpha, beta, a, b, f, qfsum )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,g14.6)' ) '  Integral of sin(x) from ', a, ' to ', b
  write ( *, '(a,i4,a)' ) '  by Fejer type rule with ', nt, ' points'
  write ( *, '(a)' ) '  of multiplicity 2.'
  write ( *, '(a,g24.16)' ) '  Quadrature formula:', qfsum

  qfsx = cos ( a ) - cos ( b )
  write ( *, '(a,g24.16)' ) '  Exact value       :', qfsx
  write ( *, '(a,g24.16)' ) '  Error             :', abs ( qfsum - qfsx )

  deallocate ( mlt )
  deallocate ( t )

  return
end
subroutine cliqfs_test ( )

!*****************************************************************************80
!
!! CLIQFS_TEST tests CLIQFS.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 January 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) alpha
  real ( kind = rk ) beta
  integer i
  integer kinda
  integer lu
  integer nt
  real ( kind = rk ), parameter :: pi = 3.14159265358979323846264338327950D+00
  real ( kind = rk ), allocatable :: t(:)
  real ( kind = rk ), allocatable :: wts(:)

  write ( *, '(a)' ) '  ----------------------------------------'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CLIQFS_TEST'
  write ( *, '(a)' ) '  Test CLIQFS.'
!
!  Number of knots.
!
  nt = 5
!
!  Set the knots in the default interval [-1,+1].
!
  allocate ( t(1:nt) )

  do i = 1, nt
    t(i) = cos ( real ( 2 * i - 1, kind = rk ) * pi / real ( 2 * nt, kind = rk ) )
  end do
!
!  Request Legendre weight function.
!
  kinda = 1
!
!  ALPHA, BETA not used in Legendre weight function but set anyway.
!
  alpha = 0.0D+00
  beta  = 0.0D+00
!
!  LU controls printing.
!  A positive value requests that we compute and print weights, and
!  conduct a moments check.
!
  lu = 6
!
!  This call returns the WTS array.
!
  allocate ( wts(1:nt) )

  call cliqfs ( nt, t, kinda, alpha, beta, lu, wts )

  deallocate ( t )
  deallocate ( wts )

  return
end
subroutine cliqf_test ( )

!*****************************************************************************80
!
!! CLIQF_TEST tests CLIQF and EIQFS.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 January 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) alpha
  real ( kind = rk ) b
  real ( kind = rk ) beta
  real ( kind = rk ), external :: f
  integer i
  integer kinda
  integer lu
  integer nt
  real ( kind = rk ), parameter :: pi = 3.14159265358979323846264338327950D+00
  real ( kind = rk ) qfsum
  real ( kind = rk ) qfsx
  real ( kind = rk ), allocatable :: t(:)
  real ( kind = rk ), allocatable :: wts(:)

  write ( *, '(a)' ) '  ----------------------------------------'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CLIQF_TEST'
  write ( *, '(a)' ) '  Test CLIQF and EIQFS.'
!
!  Number of knots.
!
  nt = 5
!
!  Set the knots in the default interval [-1,+1].
!
  allocate ( t(1:nt) )

  do i = 1, nt
    t(i) = cos ( real ( 2 * i - 1, kind = rk ) * pi / real ( 2 * nt, kind = rk ) )
  end do
!
!  Set KINDA to the Legendre weight function.
!
  kinda = 1
!
!  ALPHA, BETA not used in Legendre weight function but set anyway.
!
  alpha = 0.0D+00
  beta  = 0.0D+00
!
!  Set nonstandard interval A, B.
!
  a = -0.5D+00
  b = 2.0D+00
!
!  Shift knots from [-1,1] to [A,B].
!
  do i = 1, nt
    t(i) = ( ( b - a ) * t(i) + ( a + b ) ) / 2.0D+00
  end do
!
!  LU controls printout.
!
  lu = 6
!
!  Allocate space for WTS.
!
  allocate ( wts(1:nt) )
!
!  Call CLIQF to set up the quadrature formula.
!
  call cliqf ( nt, t, kinda, alpha, beta, a, b, lu, wts )
!
!  Call EIQFS to evaluate the quadrature formula.
!
  call eiqfs ( nt, t, wts, f, qfsum )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,g14.6)' ) '  Integral of sin(x) from ', a, ' to ', b
  write ( *, '(a,i4,a)' ) '  by Fejer type rule with ', nt, ' points'
  write ( *, '(a)' ) '  of multiplicity 1.'
  write ( *, '(a,g24.16)' ) '  Quadrature formula:', qfsum

  qfsx = cos ( a ) - cos ( b )
  write ( *, '(a,g24.16)' ) '  Exact value       :', qfsx
  write ( *, '(a,g24.16)' ) '  Error             :', abs ( qfsum - qfsx )

  deallocate ( t )
  deallocate ( wts )

  return
end
subroutine cegqf_test ( )

!*****************************************************************************80
!
!! CEGQF_TEST tests CEGQF.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 January 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) alpha
  real ( kind = rk ) b
  real ( kind = rk ) beta
  real ( kind = rk ), external :: f
  integer kinda
  integer nt
  real ( kind = rk ) qfsum
  real ( kind = rk ) qfsx

  write ( *, '(a)' ) '  ----------------------------------------'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CEGQF_TEST'
  write ( *, '(a)' ) '  Test CEGQF.'
!
!  Number of knots.
!
  nt = 12
!
!  Request exponential weight function.
!
  kinda = 7
!
!  Set ALPHA and BETA.
!
  alpha = 1.0D+00
  beta  = 0.0D+00
!
!  Set interval [A,B].
!
  a = -0.5D+00
  b = 2.0D+00
!
!  Call CEGQF to compute and evaluate the Gauss quadrature formula.
!
  call cegqf ( nt, kinda, alpha, beta, a, b, f, qfsum )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,g14.6)' ) '  Integral of x*sin(x) from ', a, ' to ', b
  write ( *, '(a,i4,a)' ) '  by Gauss-exponential rule with ', nt, ' points'
  write ( *, '(a,g24.16)' ) '  Quadrature formula:', qfsum

  qfsx = ( b - a ) * 0.5D+00 * ( cos ( a ) - cos ( b ) ) &
    + sin ( b ) + sin ( a ) - 2.0D+00 * sin ( ( a + b ) / 2.0D+00 )

  write ( *, '(a,g24.16)' ) '  Exact value       :', qfsx
  write ( *, '(a,g24.16)' ) '  Error             :', abs ( qfsum - qfsx )

  return
end
subroutine cegqfs_test ( )

!*****************************************************************************80
!
!! CEGQFS_TEST tests CEGQFS.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 January 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) alpha
  real ( kind = rk ) beta
  real ( kind = rk ), external :: f
  integer kinda
  integer nt
  real ( kind = rk ) qfsum
  real ( kind = rk ) qfsx

  write ( *, '(a)' ) '  ----------------------------------------'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CEGQFS_TEST'
  write ( *, '(a)' ) '  Test CEGQFS.'
!
!  Number of knots.
!
  nt = 12
!
!  Request exponential weight function.
!
  kinda = 7
!
!  Set ALPHA and BETA.
!
  alpha = 1.0D+00
  beta  = 0.0D+00
!
!  Call CEGQFS to compute and evaluate the Gauss quadrature formula.
!
  call cegqfs ( nt, kinda, alpha, beta, f, qfsum )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integral of x*sin(x) from -1 to +1'
  write ( *, '(a,i4,a)' ) '  by Gauss-exponential rule with ', nt, ' points'
  write ( *, '(a,g24.16)' ) '  Quadrature formula:', qfsum

  qfsx = cos ( -1.0D+00 ) - cos ( +1.0D+00 )

  write ( *, '(a,g24.16)' ) '  Exact value       :', qfsx
  write ( *, '(a,g24.16)' ) '  Error             :', abs ( qfsum - qfsx )

  return
end
subroutine cgqfs_test ( )

!*****************************************************************************80
!
!! CGQFS_TEST calls CGQFS to compute and print generalized Gauss-Hermite rules.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) alpha
  real ( kind = rk ) beta
  integer io
  integer kinda
  integer nt
  real ( kind = rk ), allocatable :: t(:)
  real ( kind = rk ), allocatable :: wts(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CGQFS_TEST'
  write ( *, '(a)' ) '  Call CGQFS to compute generalized Hermite rules.'

  nt = 15
  kinda = 6
  alpha = 1.0D+00
  beta = 0.0D+00
  io = - 6
  allocate ( t(nt) )
  allocate ( wts(nt) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  NT = ', nt
  write ( *, '(a,g14.6)' ) '  ALPHA = ', alpha

  call cgqfs ( nt, kinda, alpha, beta, io, t, wts )

  deallocate ( t )
  deallocate ( wts )

  return
end
subroutine test10 ( nt, kinda, alpha, beta )

!*****************************************************************************80
!
!! TEST10 calls CDGQF to compute a quadrature formula.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 January 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nt

  real ( kind = rk ) alpha
  real ( kind = rk ) beta
  integer i
  integer kinda
  real ( kind = rk ) t(nt)
  real ( kind = rk ) wts(nt)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  Call CDGQF to compute a quadrature formula.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  KINDA = ', kinda
  write ( *, '(a,g14.6)' ) '  ALPHA = ', alpha
  write ( *, '(a,g14.6)' ) '  BETA  = ', beta

  call cdgqf ( nt, kinda, alpha, beta, t, wts )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' Index     Abscissas                 Weights'
  write ( *, '(a)' ) ' '
  do i = 1, nt
    write ( *, '(2x,i4,2x,g24.16,2x,g24.16)' ) i, t(i), wts(i)
  end do

  return
end
subroutine test11 ( nt, kinda, alpha, beta, a, b )

!*****************************************************************************80
!
!! TEST11 calls CGQF to compute a quadrature formula.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 February 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nt

  real ( kind = rk ) a
  real ( kind = rk ) alpha
  real ( kind = rk ) b
  real ( kind = rk ) beta
  integer i
  integer kinda
  integer lu
  real ( kind = rk ) t(nt)
  real ( kind = rk ) wts(nt)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  Call CGQF to compute a quadrature formula'
  write ( *, '(a)' ) '  with nondefault values of parameters A, B.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  KINDA = ', kinda
  write ( *, '(a,g14.6)' ) '  ALPHA = ', alpha
  write ( *, '(a,g14.6)' ) '  BETA  = ', beta
  write ( *, '(a,g14.6)' ) '  A =     ', a
  write ( *, '(a,g14.6)' ) '  B =     ', b

  lu = 0
  call cgqf ( nt, kinda, alpha, beta, a, b, lu, t, wts )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' Index     Abscissas                 Weights'
  write ( *, '(a)' ) ' '
  do i = 1, nt
    write ( *, '(2x,i4,2x,g24.16,2x,g24.16)' ) i, t(i), wts(i)
  end do

  return
end
subroutine wm_test ( )

!*****************************************************************************80
!
!! WM_TEST calls WM_TESTER with various parameter values.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 November 2015
!
!  Author:
!
!    John Burkardt.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) alpha
  real ( kind = rk ) beta
  integer kinda
  integer m

  m = 5
  kinda = 1
  alpha = 0.0D+00
  beta = 0.0D+00
  call wm_tester ( m, kinda, alpha, beta )

  m = 5
  kinda = 2
  alpha = 0.0D+00
  beta = 0.0D+00
  call wm_tester ( m, kinda, alpha, beta )

  m = 5
  kinda = 3
  alpha = 0.5D+00
  beta = 0.0D+00
  call wm_tester ( m, kinda, alpha, beta )

  m = 5
  kinda =  4
  alpha = 0.25D+00
  beta = 0.75D+00
  call wm_tester ( m, kinda, alpha, beta )

  m = 5
  kinda = 5
  alpha = 2.0D+00
  beta = 0.0D+00
  call wm_tester ( m, kinda, alpha, beta )

  m = 5
  kinda = 6
  alpha = 1.0D+00
  beta = 0.0D+00
  call wm_tester ( m, kinda, alpha, beta )

  m = 5
  kinda = 7
  alpha = 2.0D+00
  beta = 0.0D+00
  call wm_tester ( m, kinda, alpha, beta )

  m = 5
  kinda = 8
  alpha = -0.5D+00
  beta = -6.0D+00
  call wm_tester ( m, kinda, alpha, beta )

  m = 5
  kinda = 9
  alpha = 0.0D+00
  beta = 0.0D+00
  call wm_tester ( m, kinda, alpha, beta )

  return
end
subroutine wm_tester ( m, kinda, alpha, beta )

!*****************************************************************************80
!
!! WM_TESTER tests WM.
!
!  Discussion:
!
!    Moment(K) = Integral ( A <= X <= B ) X^(K-1) * W(X) dx
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 November 2015
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer M, the number of moments to evaluate.
!
!    Input, integer KINDA, the rule.
!    1, Legendre,             (a,b)       1.0
!    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
!    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
!    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
!    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
!    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
!    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
!    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
!    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
!
!    Input, real ( kind = rk ) ALPHA, the value of Alpha, if needed.
!
!    Input, real BETA ( kind = rk ), the value of Beta, if needed.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m

  real ( kind = rk ) alpha
  real ( kind = rk ) beta
  integer i
  integer kinda
  real ( kind = rk ) w(m)

  call wm ( m, kinda, alpha, beta, w )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'WM_TESTER:'
  write ( *, '(a,i1)' ) '  WM_TEST computes moments for rule ', kinda
  write ( *, '(a,g14.6,a,g14.6)' ) '  with ALPHA = ', alpha, ' BETA = ', beta
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Order          Moment'
  write ( *, '(a)' ) ''
  do i = 1, m
    write ( *, '(5x,i2,2x,g14.6)' ) i - 1, w(i)
  end do

  return
end
function f ( x, i )

!*****************************************************************************80
!
!! F returns values of the integrand or its derivatives.
!
!  Discussion:
!
!    This function is an example of an integrand function.
!
!    The package can generate quadrature formulas that use derivative
!    information as well as function values.  Therefore, this routine is
!    set up to provide derivatives of any order as well as the function
!    value.  In an actual application, the highest derivative needed
!    is of order one less than the highest knot multiplicity.
!
!    In other words, in the usual case where knots are not repeated,
!    this routine only needs to return function values, not any derivatives.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 January 2010
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, the evaluation point.
!
!    Input, integer I, the order of the derivative of F to
!    be evaluated.
!
!    Output, real ( kind = rk ) F, the value of the I-th derivative of F at X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) f
  integer i
  integer l
  real ( kind = rk ) x

  l = mod ( i, 4 )

  if ( l == 0 ) then
    f = sin ( x )
  else if ( l == 1 ) then
    f = cos ( x )
  else if ( l == 2 ) then
    f = - sin ( x )
  else if ( l == 3 ) then
    f = - cos ( x )
  end if

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

