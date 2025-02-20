program main

!*****************************************************************************80
!
!! MAIN is the main program for ELLIPSE.
!
!  Discussion:
!
!    This driver computes the interpolation of the Franke function
!    on the ellipse E((C1,C2),ALPHA,BETA) = E((0.5,0.5),0.5,0.5)  
!    at the first family of Padua points. 
!
!    The ellipse has the equation:
!
!      ( ( X - C1 ) / ALPHA )^2 + ( ( Y - C2 ) / BETA )^2 = 1
!
!    The degree of interpolation DEG = 60 and the number of target 
!    points is NTG = NTG1 ^ 2 - 2 * NTG1 + 2, NTG1 = 100.  
!
!    The maps from the reference square [-1,1]^2 to the current domain 
!    are SIGMA1 and SIGMA2 with inverses ISIGM1 and ISIGM2.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!  
!  Modified:
!
!    12 February 2014
!
!  Author:
!
!    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi, 
!    Marco Vianello.
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Marco Caliari, Stefano de Marchi, Marco Vianello,
!    Algorithm 886:
!    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
!    ACM Transactions on Mathematical Software,
!    Volume 35, Number 3, October 2008, Article 21, 11 pages.
!
!  Local Parameters:
!
!    Local, integer DEGMAX, the maximum degree of interpolation.
!
!    Local, integer NPDMAX, the maximum number of Padua points
!    = (DEGMAX + 1) * (DEGMAX + 2) / 2.
!
!    Local, integer NTG1MX, the maximum value of the parameter
!    determining the number of target points.
!
!    Local, integer NTGMAX, the maximum number of target points,
!    dependent on NTG1MX.
!
!    Local, integer DEG, the degree of interpolation.
!
!    Local, integer NTG1, the parameter determining the number of 
!    target points.
!
!    Local, integer NPD, the number of Padua points 
!    = (DEG + 1) * (DEG + 2) / 2.
!
!    Local, integer NTG, the number of target points, dependent on NTG1.
!
!    Local, real ( kind = rk ) PD1(NPDMAX), the first coordinates of 
!    the Padua points.
!
!    Local, real ( kind = rk ) PD2(NPDMAX), the second coordinates of the 
!    Padua points.
!
!    Local, real ( kind = rk ) WPD(NPDMAX), the weights.
!
!    Local, real ( kind = rk ) FPD(NPDMAX), the function at the Padua points.
!
!    Workspace, real ( kind = rk ) RAUX1(DEGMAX+1)*(DEGMAX+2)).
!
!    Workspace, real ( kind = rk ) RAUX2(DEGMAX+1)*(DEGMAX+2)).
!
!    Local, real ( kind = rk ) C0(0:DEGMAX+1,0:DEGMAX+1), the coefficient matrix.
!
!    Local, real ( kind = rk ) TG1(NTGMAX), the first coordinates of the 
!    target points.
!
!    Local, real ( kind = rk ) TG2(NTGMAX), the second coordinates of the 
!    target points.
!
!    Local, real ( kind = rk ) INTFTG(NTGMAX), the values of the 
!    interpolated function.
!
!    Local, real ( kind = rk ) MAXERR, the maximum norm of the error at target 
!    points.
!
!    Local, real ( kind = rk ) ESTERR, the estimated error.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: degmax = 60
  integer, parameter :: ntg1mx = 100

  integer, parameter :: npdmax = ( degmax + 1 ) * ( degmax + 2 ) / 2
  integer, parameter :: ntgmax = ntg1mx ** 2 - 2 * ntg1mx + 2

  real ( kind = rk ) alpha
  real ( kind = rk ) beta
  real ( kind = rk ) c0(0:degmax+1,0:degmax+1)
  real ( kind = rk ) c1
  real ( kind = rk ) c2
  integer deg
  real ( kind = rk ) esterr
  integer family
  character ( len = 255 ) filename
  real ( kind = rk ) fmax
  real ( kind = rk ) fmin
  real ( kind = rk ) fpd(npdmax)
  real ( kind = rk ) franke
  real ( kind = rk ) fxy
  integer i
  real ( kind = rk ) intftg(ntgmax)
  real ( kind = rk ) isigm1
  real ( kind = rk ) isigm2
  real ( kind = rk ) ixy
  real ( kind = rk ) maxdev
  real ( kind = rk ) maxerr
  real ( kind = rk ) mean
  integer npd
  integer ntg
  integer ntg1
  real ( kind = rk ) pd1(npdmax)
  real ( kind = rk ) pd2(npdmax)
  real ( kind = rk ) pd2val
  real ( kind = rk ) r8_huge
  real ( kind = rk ) raux1((degmax+1)*(degmax+2))
  real ( kind = rk ) raux2((degmax+1)*(degmax+2))
  real ( kind = rk ) sigma1
  real ( kind = rk ) sigma2
  real ( kind = rk ) tg1(ntgmax)
  real ( kind = rk ) tg2(ntgmax)
  real ( kind = rk ) wpd(npdmax)
  real ( kind = rk ) x
  real ( kind = rk ) y

  alpha = 0.5D+00
  beta = 0.5D+00
  c1 = 0.5D+00
  c2 = 0.5D+00
  family = 1
  deg = 60
  ntg1 = 100

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ELLIPSE():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Interpolation of the Franke function'
  write ( *, '(a)' ) '  on the disk with center = (0.5,0.5) and radius = 0.5'
  write ( *, '(a,i6)' ) '  of degree = ', deg

  if ( degmax < deg ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ELLIPSE - Fatal error!'
    write ( *, '(a)' ) '  DEGMAX < DEG.'
    write ( *, '(a,i6)' ) '  DEG =    ', deg
    write ( *, '(a,i6)' ) '  DEGMAX = ', degmax
    stop 1
  end if
!     
!  Build the first family of Padua points in the square [-1,1]^2.
!   
  call pdpts ( deg, pd1, pd2, wpd, npd )
!     
!  Compute the Franke function at Padua points mapped to the region.
!  
  do i = 1, npd
    x = sigma1 ( pd1(i), pd2(i), c1, c2, alpha, beta )
    y = sigma2 ( pd1(i), pd2(i), c1, c2, alpha, beta )
    fpd(i) = franke ( x, y )
  end do
!
!  Write X, Y, F(X,Y) to a file.
!
  filename = 'ellipse_fpd.txt'
  open ( unit = 10, file = filename, status = 'replace' )
  do i = 1, npd
    x = sigma1 ( pd1(i), pd2(i), c1, c2, alpha, beta )
    y = sigma2 ( pd1(i), pd2(i), c1, c2, alpha, beta )
    write ( 10, '(3g14.6)' ) x, y, fpd(i)
  end do
  close ( unit = 10 )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Wrote F(x,y) at Padua points in "' &
    // trim ( filename ) // '".'
!     
!  Compute the matrix C0 of the coefficients in the bivariate
!  orthonormal Chebyshev basis.
!    
  call padua2 ( deg, degmax, npd, wpd, fpd, raux1, raux2, c0, esterr )
!     
!  Evaluate the target points in the region.
!     
  call target ( c1, c2, alpha, beta, ntg1, ntgmax, tg1, tg2, ntg )
!     
!  Evaluate the interpolant at the target points.
!    
  do i = 1, ntg
    x = isigm1 ( tg1(i), tg2(i), c1, c2, alpha, beta )
    y = isigm2 ( tg1(i), tg2(i), c1, c2, alpha, beta )
    intftg(i) = pd2val ( deg, degmax, c0, x, y )
  end do
!
!  Write the function value at target points to a file.
!
  filename = 'ellipse_ftg.txt'
  open ( unit = 10, file = filename, status = 'replace' )
  do i = 1, ntg
    write ( 10, '(3g14.6)' ) tg1(i), tg2(i), franke ( tg1(i), tg2(i) )
  end do
  close ( unit = 10 )
  write ( *, '(a)' ) '  Wrote F(x,y) at target points in "' &
    // trim ( filename ) // '".'
!
!  Write the interpolated function value at target points to a file.
!
  filename = 'ellipse_itg.txt'
  open ( unit = 10, file = filename, status = 'replace' )
  do i = 1, ntg
    write ( 10, '(3g14.6)' ) tg1(i), tg2(i), intftg(i)
  end do
  close ( unit = 10 )
  write ( *, '(a)' ) '  Wrote I(F)(x,y) at target points in "' &
    // trim ( filename ) // '".'
!
!  Compute the error relative to the max deviation from the mean.
!     
  maxerr = 0.0D+00
  mean = 0.0D+00
  fmax = - r8_huge ( )
  fmin = + r8_huge ( )

  do i = 1, ntg
    fxy = franke ( tg1(i), tg2(i) )
    ixy = intftg(i)
    maxerr = max ( maxerr, abs ( fxy - ixy ) )
    mean = mean + fxy
    fmax = max ( fmax, fxy )
    fmin = min ( fmin, fxy )
  end do
 
  if ( fmax == fmin ) then
    maxdev = 1.0D+00
  else
    mean = mean / real ( ntg, kind = rk )
    maxdev = max ( fmax - mean, mean - fmin )
  end if
!
!  Print error ratios.
!
  write ( *, '(a)' ) ''
  write ( *, '(a,e10.4)' ) '  Estimated error:  ', esterr / maxdev
  write ( *, '(a,e10.4)' ) '  Actual error:     ', maxerr / maxdev
  write ( *, '(a,e10.4)' ) '  Expected error:   ', 0.1769D-09
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ELLIPSE():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
function sigma1 ( t1, t2, c1, c2, alpha, beta )

!*****************************************************************************80
!
!! SIGMA1 maps first coordinate from square to ellipse.
!
!  Discussion:
!
!    This function returns the first component of the map 
!    from the square [-1,1]^2 to the ellipse E((C1,C2),ALPHA,BETA).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!  
!  Modified:
!
!    12 February 2014
!
!  Author:
!
!    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi, 
!    Marco Vianello.
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Marco Caliari, Stefano de Marchi, Marco Vianello,
!    Algorithm 886:
!    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
!    ACM Transactions on Mathematical Software,
!    Volume 35, Number 3, October 2008, Article 21, 11 pages.
!
!  Parameters:
!
!    Input, real ( kind = rk ) T1, T2, the coordinates of a point in the square.
!
!    Input, real ( kind = rk ) C1, C2, ALPHA, BETA, the center and scale
!    parameters of the ellipse.
!
!    Output, real ( kind = rk ) SIGMA1, the X coordinate of the corresponding
!    point in the ellipse.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) alpha
  real ( kind = rk ) beta
  real ( kind = rk ) c1
  real ( kind = rk ) c2
  real ( kind = rk ) phi
  real ( kind = rk ) sigma1
  real ( kind = rk ) t1
  real ( kind = rk ) t2

  sigma1 = c1 - alpha * t2 * sin ( phi ( t1 ) )

  return
end
function isigm1 ( sigma1, sigma2, c1, c2, alpha, beta )

!*****************************************************************************80
!
!! ISIGM1 maps the first coordinate from the ellipse to the square.
!
!  Discussion:
!
!    This function returns the first component of the map 
!    from the ellipse E((C1,C2),ALPHA,BETA) to the square [-1,1]^2. 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!  
!  Modified:
!
!    09 February 2014
!
!  Author:
!
!    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi, 
!    Marco Vianello.
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Marco Caliari, Stefano de Marchi, Marco Vianello,
!    Algorithm 886:
!    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
!    ACM Transactions on Mathematical Software,
!    Volume 35, Number 3, October 2008, Article 21, 11 pages.
!
!  Parameters:
!
!    Input, real ( kind = rk ) SIGMA1, SIGMA2, the coordinates of a point 
!    in the ellipse.
!
!    Input, real ( kind = rk ) C1, C2, ALPHA, BETA, the center and scale
!    parameters of the ellipse.
!
!    Output, real ( kind = rk ) ISIGM1, the X coordinate of the corresponding
!    point in the square.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) alpha
  real ( kind = rk ) beta
  real ( kind = rk ) c1
  real ( kind = rk ) c2
  real ( kind = rk ) iphi
  real ( kind = rk ) isigm1
  real ( kind = rk ) sigma1
  real ( kind = rk ) sigma2

  if ( sigma2 == c2 ) then
    isigm1 = 1.0D+00
  else
    isigm1 = iphi ( atan ( beta * ( c1 - sigma1 ) / &
      ( alpha * ( sigma2 - c2 ) ) ) )
  end if   

  return
end
function sigma2 ( t1, t2, c1, c2, alpha, beta )

!*****************************************************************************80
!
!! SIGMA2 maps the second coordinate from square to ellipse.
!
!  Discussion:
!
!    This function returns the second component of the map 
!    from the square [-1,1]^2 to the ellipse E((C1,C2),ALPHA,BETA).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!  
!  Modified:
!
!    12 February 2014
!
!  Author:
!
!    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi, 
!    Marco Vianello.
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Marco Caliari, Stefano de Marchi, Marco Vianello,
!    Algorithm 886:
!    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
!    ACM Transactions on Mathematical Software,
!    Volume 35, Number 3, October 2008, Article 21, 11 pages.
!
!  Parameters:
!
!    Input, real ( kind = rk ) T1, T2, the coordinates of a point in the square.
!
!    Input, real ( kind = rk ) C1, C2, ALPHA, BETA, the center and scale
!    parameters of the ellipse.
!
!    Output, real ( kind = rk ) SIGMA2, the Y coordinate of the corresponding
!    point in the ellipse.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) alpha
  real ( kind = rk ) beta
  real ( kind = rk ) c1
  real ( kind = rk ) c2
  real ( kind = rk ) phi
  real ( kind = rk ) sigma2
  real ( kind = rk ) t1
  real ( kind = rk ) t2

  sigma2 = c2 + beta * t2 * cos ( phi ( t1 ) )

  return
end
function isigm2 ( sigma1, sigma2, c1, c2, alpha, beta )

!*****************************************************************************80
!
!! ISIGM2 maps second coordinate from ellipse to the square.
!
!  Discussion:
!
!    This function returns the second component of the map 
!    from the ellipse E((C1,C2),ALPHA,BETA) to the square [-1,1]^2. 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!  
!  Modified:
!
!    09 February 2014
!
!  Author:
!
!    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi, 
!    Marco Vianello.
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Marco Caliari, Stefano de Marchi, Marco Vianello,
!    Algorithm 886:
!    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
!    ACM Transactions on Mathematical Software,
!    Volume 35, Number 3, October 2008, Article 21, 11 pages.
!
!  Parameters:
!
!    Input, real ( kind = rk ) SIGMA1, SIGMA2, the coordinates of a point 
!    in the ellipse.
!
!    Input, real ( kind = rk ) C1, C2, ALPHA, BETA, the center and scale
!    parameters of the ellipse.
!
!    Output, real ( kind = rk ) ISIGM2, the Y coordinate of the corresponding
!    point in the square.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) alpha
  real ( kind = rk ) beta
  real ( kind = rk ) c1
  real ( kind = rk ) c2
  real ( kind = rk ) isigm2
  real ( kind = rk ) sigma1
  real ( kind = rk ) sigma2

  if ( sigma2 == c2 ) then
    isigm2 = ( c1 - sigma1 ) / alpha
  else
    isigm2 = sqrt ( beta ** 2 * ( c1 - sigma1 ) ** 2 + &
      alpha ** 2 * ( c2 - sigma2 ) ** 2 ) / ( alpha * beta ) * &
      sign ( 1.0D+00, sigma2 - c2 )
  end if   

  return
end
function phi ( x )

!*****************************************************************************80
!
!! PHI maps from [-1,+1] to [-pi/2,+pi/2].
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!  
!  Modified:
!
!    12 February 2014
!
!  Author:
!
!    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi, 
!    Marco Vianello.
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Marco Caliari, Stefano de Marchi, Marco Vianello,
!    Algorithm 886:
!    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
!    ACM Transactions on Mathematical Software,
!    Volume 35, Number 3, October 2008, Article 21, 11 pages.
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, a point in [-1,+1];
!
!    Output, real ( kind = rk ) PHI, a corresponding point in [-pi/2,+pi/2].
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) phi
  real ( kind = rk ), parameter :: pi = 3.1415926535897931D+00
  real ( kind = rk ) x

  phi = pi * x / 2.0D+00

  return
end
function iphi ( x )

!*****************************************************************************80
!
!! IPHI maps from [-pi/2,+pi/2] to [-1,+1].
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!  
!  Modified:
!
!    12 February 2014
!
!  Author:
!
!    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi, 
!    Marco Vianello.
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Marco Caliari, Stefano de Marchi, Marco Vianello,
!    Algorithm 886:
!    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
!    ACM Transactions on Mathematical Software,
!    Volume 35, Number 3, October 2008, Article 21, 11 pages.
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, a point in [-pi/2,+pi/2].
!
!    Output, real ( kind = rk ) IPHI, a corresponding point in [-1,+1].
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) iphi
  real ( kind = rk ), parameter :: pi = 3.1415926535897931D+00
  real ( kind = rk ) x

  iphi = 2.0D+00 * x / pi

  return
end
subroutine target ( c1, c2, alpha, beta, ntg1, ntgmax, tg1, tg2, ntg )

!*****************************************************************************80
!
!! TARGET returns the target points on the ellipse.
!
!  Discussion:
!
!    Target points on the ellipse E((C1,C2),ALPHA,BETA).
!    The number of target points is NTG = NTG1^2 - 2 * NTG1 + 2.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!  
!  Modified:
!
!    12 February 2014
!
!  Author:
!
!    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi, 
!    Marco Vianello.
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Marco Caliari, Stefano de Marchi, Marco Vianello,
!    Algorithm 886:
!    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
!    ACM Transactions on Mathematical Software,
!    Volume 35, Number 3, October 2008, Article 21, 11 pages.
!
!  Parameters:
!
!    Input, real ( kind = rk ) C1, C2, ALPHA, BETA, the center and scale
!    parameters of the ellipse.
!
!    Input, integer NTG1, a parameter determining the number 
!    of target points.  2 <= NTG1.
!
!    Input, integer NTGMAX, the maximum number of target points.
!
!    Output, real ( kind = rk ) TG1(NTG), TG2(NTG), the X and Y coordinates
!    of the target points.
!
!    Output, integer NTG, the number of target points computed.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ntgmax

  real ( kind = rk ) alpha
  real ( kind = rk ) beta
  real ( kind = rk ) c1
  real ( kind = rk ) c2
  integer i
  integer j
  integer ntg
  integer ntg1
  real ( kind = rk ) tg1(ntgmax)
  real ( kind = rk ) tg2(ntgmax)

  if ( ntg1 < 2 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'TARGET - Fatal error!'
    write ( *, '(a)' ) '  NTG1 < 2'
    write ( *, '(a,i4)' ) '  NTG1 = ', ntg1
    stop 1
  end if

  if ( ntgmax < ntg1 * ntg1 - 2 * ntg1 + 2 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'TARGET - Fatal error!'
    write ( *, '(a)' ) '  NTGMAX < NTG1 * NTG1 - 2 * NTG1 + 2.'
    write ( *, '(a,i4)' ) '  NTG1 = ', ntg1
    write ( *, '(a,i4)' ) '  NTGMAX = ', ntgmax
    stop 1
  end if      

  i = 1
  j = 1
  ntg = 0
  ntg = ntg + 1

  tg1(ntg) = alpha * &
    ( - 1.0D+00 + real ( i - 1, kind = rk ) * 2.0D+00 &
    / real ( ntg1 - 1, kind = rk ) ) + c1

  tg2(ntg) =  beta * &
    ( - 1.0D+00 + real ( j - 1, kind = rk ) * 2.0D+00 &
    / real ( ntg1 - 1, kind = rk ) ) * &
    sqrt ( 1.0D+00 - ( - 1.0D+00 + real ( i - 1, kind = rk ) * 2.0D+00 &
    / real ( ntg1 - 1, kind = rk ) ) ** 2 ) &
    + c2

  do i = 2, ntg1 - 1
    do j = 1, ntg1
      ntg = ntg + 1

      tg1(ntg) = alpha * &
        ( - 1.0D+00 + real ( i - 1, kind = rk ) * 2.0D+00 &
        / real ( ntg1 - 1, kind = rk ) ) + c1

      tg2(ntg) =  beta * &
        ( - 1.0D+00 + real ( j - 1, kind = rk ) * 2.0D+00 &
        / real ( ntg1 - 1, kind = rk ) ) * &
        sqrt ( 1.0D+00 - ( - 1.0D+00 + real ( i - 1, kind = rk ) &
        * 2.0D+00 / real ( ntg1 - 1, kind = rk ) ) ** 2 ) &
        + c2

    end do
  end do

  i = ntg1
  j = 1
  ntg = ntg + 1

  tg1(ntg) = alpha * &
    ( -1.0D+00 + real ( i - 1, kind = rk ) * 2.0D+00 &
    / real ( ntg1 - 1, kind = rk ) ) + c1

  tg2(ntg) = beta * &
    ( -1.0D+00 + real ( j - 1, kind = rk ) * 2.0D+00 &
    / real ( ntg1 - 1, kind = rk ) ) * &
    sqrt ( 1.0D+00 - ( - 1.0D+00 + real ( i - 1, kind = rk ) * 2.0D+00 &
    / real ( ntg1 - 1, kind = rk ) ) ** 2 ) &
    + c2

  return
end
