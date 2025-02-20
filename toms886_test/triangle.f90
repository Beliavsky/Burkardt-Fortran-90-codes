program main

!*****************************************************************************80
!
!! MAIN is the main program for TRIANGLE.
!
!  Discussion:
!
!    This driver computes the interpolation of the Franke function
!    on the triangle T(U,V,W) with vertices U=(U1,U2)=(0,0), 
!    V=(V1,V2)=(1,0) and W=(W1,W2)=(0,1) (unit triangle) 
!    at the first family of Padua points. 
!
!    The degree of interpolation is DEG = 60 and the number of target 
!    points is NTG = NTG1 ** 2 - NTG1 + 1, NTG1 = 100.
!
!    The maps from the reference square [-1,1]^2 to the triangle
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
!  Parameters:
!
!    Local, integer DEGMAX, the maximum degree of interpolation.
!
!    Local, integer NPDMAX, the maximum number of Padua points
!    = (DEGMAX + 1) * (DEGMAX + 2) / 2.
!
!    Local, integer NTG1MX, the maximum value of the parameter determining 
!    the number of target points.
!
!    Local, integer NTGMAX, the maximum number of target points,
!    dependent on NTG1MX.
!
!    Local, integer DEG, the degree of interpolation.
!
!    Local, integer NTG1, the parameter determining the number 
!   of target points.
!
!    Local, integer NPD, the number of Padua points = (DEG + 1) * (DEG + 2) / 2.
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
  integer, parameter :: ntgmax = ntg1mx ** 2 - ntg1mx + 1

  real ( kind = rk ) c0(0:degmax+1,0:degmax+1)
  integer deg
  real ( kind = rk ) esterr
  integer family
  character * ( 255 ) filename
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
  real ( kind = rk ) u1
  real ( kind = rk ) u2
  real ( kind = rk ) v1
  real ( kind = rk ) v2
  real ( kind = rk ) w1
  real ( kind = rk ) w2
  real ( kind = rk ) wpd(npdmax)
  real ( kind = rk ) x
  real ( kind = rk ) y

  u1 = 0.0D+00
  u2 = 0.0D+00
  v1 = 1.0D+00
  v2 = 0.0D+00
  w1 = 0.0D+00
  w2 = 1.0D+00
  family = 1
  deg = 60
  ntg1 = 100

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIANGLE:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Interpolation of the Franke function'
  write ( *, '(a)' ) '  on the unit triangle T((0,0),(1,0),(0,1))'
  write ( *, '(a,i6)' ) '  at degree = ', deg

  if ( degmax < deg ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGLE - Fatal error!'
    write ( *, '(a)' ) '  DEGMAX < DEG.'
    write ( *, '(a,i6)' ) '  DEG =    ', deg
    write ( *, '(a,i6)' ) '  DEGMAX = ', degmax
    stop 1
  end if
!
!  Build the first family of Padua points in the square [-1,1]^2
!     
  call pdpts ( deg, pd1, pd2, wpd, npd )
!     
!  Compute the Franke function at Padua points mapped to T(U,V,W).
!   
  do i = 1, npd
    x = sigma1 ( pd1(i), pd2(i), u1, u2, v1, v2, w1, w2 )
    y = sigma2 ( pd1(i), pd2(i), u1, u2, v1, v2, w1, w2 )
    fpd(i) = franke ( x, y )
  end do
!
!  Write X, Y, F(X,Y) to a file.
!
  filename = 'triangle_fpd.txt'
  open ( unit = 10, file = filename, status = 'replace' )
  do i = 1, npd
    x = sigma1 ( pd1(i), pd2(i), u1, u2, v1, v2, w1, w2 )
    y = sigma2 ( pd1(i), pd2(i), u1, u2, v1, v2, w1, w2 )
    write ( 10, '(3g14.6)' ) x, y, fpd(i)
  end do
  close ( unit = 10 )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Wrote F(x,y) at Padua points in "' &
    // trim ( filename ) // '".'
!     
!  Compute the matrix C0 of the coefficients in the bivariate
!  orthonormal Chebyshev basis
!     
  call padua2 ( deg, degmax, npd, wpd, fpd, raux1, raux2, c0, esterr )
!     
!  Build the set of target points on T(U,V,W)
!     
  call target ( u1, u2, v1, v2, w1, w2, ntg1, ntgmax, tg1, tg2, ntg )
!     
!  Evaluate the interpolant at the target points.
!    
  do i = 1, ntg
    x = isigm1 ( tg1(i), tg2(i), u1, u2, v1, v2, w1, w2 )
    y = isigm2 ( tg1(i), tg2(i), u1, u2, v1, v2, w1, w2 )
    intftg(i) = pd2val ( deg, degmax, c0, x, y )
  end do
!
!  Write the function value at target points to a file.
!
  filename = 'triangle_ftg.txt'
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
  filename = 'triangle_itg.txt'
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
  write ( *, '(a,e10.4)' ) '  Expected error:   ', 0.1226D-09
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIANGLE:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
function sigma1 ( t1, t2, u1, u2, v1, v2, w1, w2 )

!*****************************************************************************80
!
!! SIGMA1 maps first coordinate from square to triangle.
!
!  Discussion:
!
!    This function returns the first component of the map
!    from the square [-1,1]^2 to the triangle T(U,V,W). 
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
!    Input, real ( kind = rk ) U1, U2, V1, V2, W1, W2, the coordinates of the 
!    vertices of the triangle.
!
!    Output, real ( kind = rk ) SIGMA1, the X coordinate of the corresponding
!    point in the triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) sigma1
  real ( kind = rk ) t1
  real ( kind = rk ) t2
  real ( kind = rk ) u1
  real ( kind = rk ) u2
  real ( kind = rk ) v1
  real ( kind = rk ) v2
  real ( kind = rk ) w1
  real ( kind = rk ) w2

  sigma1 = ( v1 - u1 ) * ( 1.0D+00 + t1 ) &
    * ( 1.0D+00 - t2 ) / 4.0D+00 &
    + ( w1 - u1 ) * ( 1.0D+00 + t2 ) / 2.0D+00 + u1

  return
end
function isigm1 ( sigma1, sigma2, u1, u2, v1, v2, w1, w2 )

!*****************************************************************************80
!
!! ISIGM1 maps first coordinate from triangle to the square.
!
!  Discussion:
!
!    This functions returns the first component of the map
!    from the the triangle T(U,V,W) to the square [-1,1]^2.
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
!    Input, real ( kind = rk ) SIGMA1, SIGMA2, the coordinates of a point 
!    in the triangle.
!
!    Input, real ( kind = rk ) U1, U2, V1, V2, W1, W2, the coordinates of the 
!    vertices of the triangle.
!
!    Output, real ( kind = rk ) ISIGM1, the X coordinate of the corresponding
!    point in the square.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) isigm1
  real ( kind = rk ) rho1
  real ( kind = rk ) rho2
  real ( kind = rk ) sigma1
  real ( kind = rk ) sigma2
  real ( kind = rk ) u1
  real ( kind = rk ) u2
  real ( kind = rk ) v1
  real ( kind = rk ) v2
  real ( kind = rk ) w1
  real ( kind = rk ) w2

  rho1 = ( sigma1 * ( w2 - u2 ) - sigma2 * ( w1 - u1 ) &
    + ( w1 - u1 ) * u2 - ( w2 - u2 ) * u1 ) / &
    ( ( v1 - u1 ) * ( w2 - u2 ) - ( v2 - u2 ) * ( w1 - u1 ) )

  rho2 = ( sigma1 * ( v2 - u2 ) - sigma2 * ( v1 - u1 ) &
    + ( v1 - u1 ) * u2 - ( v2 - u2 ) * u1 ) / &
    ( ( w1 - u1 ) * ( v2 - u2 ) - ( w2 - u2 ) * ( v1 - u1 ) )

  if ( rho2 == 1.0D+00 ) then
    isigm1 = 0.0D+00
  else
    isigm1 = 2.0D+00 * rho1 / ( 1.0D+00 - rho2 ) - 1.0D+00
  end if  
 
  return
end
function sigma2 ( t1, t2, u1, u2, v1, v2, w1, w2 )

!*****************************************************************************80
!
!! SIGMA2 maps the second coordinate from square to triangle.
!
!  Discussion:
!
!    This functions returns the second component of the map
!    from the square [-1,1]^2 to the triangle T(U,V,W).
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
!    Input, real ( kind = rk ) U1, U2, V1, V2, W1, W2, the coordinates of the 
!    vertices of the triangle.
!
!    Output, real ( kind = rk ) SIGMA2, the Y coordinate of the corresponding
!    point in the triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) sigma2
  real ( kind = rk ) t1
  real ( kind = rk ) t2
  real ( kind = rk ) u1
  real ( kind = rk ) u2
  real ( kind = rk ) v1
  real ( kind = rk ) v2
  real ( kind = rk ) w1
  real ( kind = rk ) w2

  sigma2 = ( v2 - u2 ) * ( 1.0D+00 + t1 ) &
    * ( 1.0D+00 - t2 ) / 4.0D+00 + ( w2 - u2 ) &
    * ( 1.0D+00 + t2 ) / 2.0D+00 + u2

  return
end
function isigm2 ( sigma1, sigma2, u1, u2, v1, v2, w1, w2 )

!*****************************************************************************80
!
!! ISIGM2 maps second coordinate from triangle to the square.
!
!  Discussion:
!
!    This functions returns the second component of the map
!    from the the triangle T(U,V,W) to the square [-1,1]^2.
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
!    Input, real ( kind = rk ) SIGMA1, SIGMA2, the coordinates of a point 
!    in the ellipse.
!
!    Input, real ( kind = rk ) U1, U2, V1, V2, W1, W2, the coordinates of the 
!    vertices of the triangle.
!
!    Output, real ( kind = rk ) ISIGM2, the Y coordinate of the corresponding
!    point in the triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) isigm2
  real ( kind = rk ) rho2
  real ( kind = rk ) sigma1
  real ( kind = rk ) sigma2
  real ( kind = rk ) u1
  real ( kind = rk ) u2
  real ( kind = rk ) v1
  real ( kind = rk ) v2
  real ( kind = rk ) w1
  real ( kind = rk ) w2

  rho2 = ( sigma1 * ( v2 - u2 ) - &
    sigma2 * ( v1 - u1) + ( v1 - u1 ) * u2 - ( v2 - u2 ) * u1 ) / &
    ( ( w1 - u1 ) * ( v2 - u2 ) - ( w2 - u2 ) * ( v1 - u1 ) )

  if ( rho2 == 1.0D+00 ) then
    isigm2 = 1.0D+00
  else
    isigm2 = 2.0D+00 * rho2 - 1.0D+00
  end if 
  
  return
end
subroutine target ( u1, u2, v1, v2, w1, w2, ntg1, ntgmax, tg1, tg2, ntg )

!*****************************************************************************80
!
!! TARGET returns the target points on the triangle.
!
!  Discussion:
!
!    Target points on the triangle T(U,V,W).
!    The number of target points is NTG = NTG1^2 - NGT1 + 1.
!
!  Licensing:
!
!    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi, 
!    Marco Vianello.
!    This FORTRAN90 version by John Burkardt.
!  
!  Modified:
!
!    12 February 2014
!
!  Author:
!
!    Marco Caliari, Stefano De Marchi, Marco Vianello
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
!    Input, real ( kind = rk ) U1, U2, V1, V2, W1, W2, the coordinates of the 
!    vertices of the triangle.
!
!    Input, integer NTG1, a parameter determining the number 
!    of target points
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

  integer i
  integer j
  integer ntg
  integer ntg1
  real ( kind = rk ) tg1(ntgmax)
  real ( kind = rk ) tg2(ntgmax)
  real ( kind = rk ) u1
  real ( kind = rk ) u2
  real ( kind = rk ) v1
  real ( kind = rk ) v2
  real ( kind = rk ) w1
  real ( kind = rk ) w2

  if ( ntg1 < 2 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'TARGET - Fatal error!'
    write ( *, '(a)' ) '  NTG1 < 2'
    write ( *, '(a,i4)' ) '  NTG1 = ', ntg1
    stop 1
  end if  
   
  if ( ntgmax < ntg1 ** 2 - ntg1 + 1 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'TARGET - Fatal error!'
    write ( *, '(a)' ) '  NTGMAX < NTG1 * NTG1 - NTG1 + 1.'
    write ( *, '(a,i4)' ) '  NTG1 = ', ntg1
    write ( *, '(a,i4)' ) '  NTGMAX = ', ntgmax
    stop 1
  end if    

  ntg = 0
  do i = 1, ntg1 - 1
    do j = 1, ntg1

      ntg = ntg + 1

      tg1(ntg) = ( v1 - u1 ) * real ( i - 1, kind = rk ) / real ( ntg1 - 1, kind = rk ) &
        + ( w1 - u1 ) * ( real ( j - 1, kind = rk ) / real ( ntg1 - 1, kind = rk ) ) &
        * ( 1.0D+00 - real ( i - 1, kind = rk ) / real ( ntg1 - 1, kind = rk ) ) + u1

      tg2(ntg) = ( v2 - u2 ) * real ( i - 1, kind = rk ) / real ( ntg1 - 1, kind = rk ) &
        + ( w2 - u2 ) * ( real ( j - 1, kind = rk ) / real ( ntg1 - 1, kind = rk ) ) &
        * ( 1.0D+00 - real ( i - 1, kind = rk ) / real ( ntg1 - 1, kind = rk ) ) + u2

    end do
  end do

  i = ntg1
  j = 1
  ntg = ntg + 1

  tg1(ntg) = ( v1 - u1 ) * real ( i - 1, kind = rk ) / real ( ntg1 - 1, kind = rk ) &
    + ( w1 - u1 ) * ( real ( j - 1, kind = rk ) / real ( ntg1 - 1, kind = rk ) ) &
    * ( 1.0D+00 - real ( i - 1, kind = rk ) / real ( ntg1 - 1, kind = rk ) ) + u1

  tg2(ntg) = ( v2 - u2 ) * real ( i - 1, kind = rk ) / real ( ntg1 - 1, kind = rk ) &
    + ( w2 - u2 ) * ( real ( j - 1, kind = rk ) / real ( ntg1 - 1, kind = rk ) ) &
    * ( 1.0D+00 - real ( i - 1, kind = rk ) / real ( ntg1 - 1, kind = rk ) ) + u2

  return
end
