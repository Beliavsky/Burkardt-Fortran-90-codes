program main

!*****************************************************************************80
!
!! MAIN is the main program for RECTANGLE.
!
!  Discussion:
!
!    This driver computes the interpolation of the Franke function
!    on the rectangle R(A,B) = [A1,B1] x [A2,B2] with A=(A1,A2)=(0,0) 
!    and B=(B1,B2)=(1,1) (unit square) at the FAMILY = 1 of Padua points. 
!
!    The degree of interpolation is DEG = 60 and the number of target 
!    points is NTG = NTG1^2, NTG1 = 100. 
!
!    The maps from the reference square [-1,1]^2 to the rectangle
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
!    Local, integer DEG, the degree of interpolation,
!
!    Local, integer NTG1, the parameter determining the number 
!    of target points.
!
!    Local, integer FAMILY, specifies the desired family of Padua points.
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
  integer, parameter :: ntgmax = ntg1mx ** 2

  real ( kind = rk ) a1
  real ( kind = rk ) a2
  real ( kind = rk ) b1
  real ( kind = rk ) b2
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
  real ( kind = rk ) wpd(npdmax)
  real ( kind = rk ) x
  real ( kind = rk ) y

  a1 = 0.0D+00
  a2 = 0.0D+00
  b1 = 1.0D+00
  b2 = 1.0D+00
  family = 1
  deg = 60
  ntg1 = 100
 
  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RECTANGLE:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Interpolation of the Franke function'
  write ( *, '(a)' ) '  on the unit square [0,1] x [0,1]'
  write ( *, '(a,i6)' ) '  of degree = ', deg

  if ( degmax < deg ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'RECTANGLE - Fatal error!'
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
!  Compute the Franke function at Padua points mapped to the region.
!
  do i = 1, npd    
    x = sigma1 ( pd1(i), pd2(i), a1, a2, b1, b2, family, deg )
    y = sigma2 ( pd1(i), pd2(i), a1, a2, b1, b2, family, deg )
    fpd(i) = franke ( x, y )
  end do
!
!  Write X, Y, F(X,Y) to a file.
!
  filename = 'rectangle_fpd.txt'
  open ( unit = 10, file = filename, status = 'replace' )
  do i = 1, npd
    x = sigma1 ( pd1(i), pd2(i), a1, a2, b1, b2, family, deg )
    y = sigma2 ( pd1(i), pd2(i), a1, a2, b1, b2, family, deg )
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
!  Evaluate the target points in the region.
!     
  call target ( a1, b1, a2, b2, ntg1, ntgmax, tg1, tg2, ntg )
!     
!  Evaluate the interpolant at the target points.
! 
  do i = 1, ntg    
    x = isigm1 ( tg1(i), tg2(i), a1, a2, b1, b2, family, deg )
    y = isigm2 ( tg1(i), tg2(i), a1, a2, b1, b2, family, deg )
    intftg(i) = pd2val ( deg, degmax, c0, x, y )
  end do
!
!  Write the function value at target points to a file.
!
  filename = 'rectangle_ftg.txt'
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
  filename = 'rectangle_itg.txt'
  open ( unit = 10, file = filename, status = 'replace' )
  do i = 1, ntg
    write ( 10, '(3g14.6)' ) tg1(i), tg2(i), intftg(i)
  end do
  close ( unit = 10 )
  write ( *, '(a)' ) '  Wrote I(F)(x,y) at target points in "' &
    // trim ( filename ) // '".'
!  
!  Compute the error relative to the max deviation from the mean
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
  write ( *, '(a,e10.4)' ) '  Expected error:   ', 0.2468D-10
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RECTANGLE:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
function sigma1 ( t1, t2, a1, a2, b1, b2, family, deg )

!*****************************************************************************80
!
!! SIGMA1 maps first coordinate from square to the rectangle.
!
!  Discussion:
!
!    This function returns the first component of the map 
!    from the square [-1,1]^2 to the rectangle [A1,B1] x [A2,B2]. 
!    FAMILY and DEG select the rotation in order to get 
!    the corresponding FAMILY of Padua points at degree DEG.
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
!    Input, real ( kind = rk ) A1, B1, A2, B2, the coordinates of the extreme
!    corners of the rectangle.
!
!    Input, integer FAMILY, DEG, select the family of Padua points at 
!    degree DEG.
!
!    Output, real ( kind = rk ) SIGMA1, the X coordinate of the corresponding
!    point in the rectangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a1
  real ( kind = rk ) a2
  real ( kind = rk ) b1
  real ( kind = rk ) b2
  integer deg
  integer family
  real ( kind = rk ), parameter :: pi = 3.1415926535897931D+00
  real ( kind = rk ) sigma1
  real ( kind = rk ) t1
  real ( kind = rk ) t2
  real ( kind = rk ) theta

  theta = real ( 2 * mod ( deg, 2 ) - 1, kind = rk ) &
    * real ( family - 1, kind = rk ) * pi / 2.0D+00
  sigma1 = t1 * cos ( theta ) - t2 * sin ( theta )
  sigma1 = ( ( b1 - a1 ) * sigma1 + ( b1 + a1 ) ) / 2.0D+00

  return
end
function isigm1 ( sigma1, sigma2, a1, a2, b1, b2, family, deg )

!*****************************************************************************80
!
!! ISIGM1 maps first coordinate from the rectangle to the square.
!
!  Discussion:
!
!    This function returns the first component of the map 
!    from the rectangle [A1,B1] x [A2,B2] to the square [-1,1]^2. 
!    FAMILY and DEG select the rotation in order to get 
!    the corresponding FAMILY of Padua points at degree DEG.
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
!    in the rectangle.
!
!    Input, real ( kind = rk ) A1, B1, A2, B2, the coordinates of the extreme
!    corners of the rectangle.
!
!    Input, integer FAMILY, DEG, select the family of Padua points at 
!    degree DEG.
!
!    Output, real ( kind = rk ) ISIGM1, the X coordinate of the corresponding
!    point in the square.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a1
  real ( kind = rk ) a2
  real ( kind = rk ) b1
  real ( kind = rk ) b2
  integer deg
  integer family
  real ( kind = rk ) isigm1
  real ( kind = rk ) isigm2
  real ( kind = rk ), parameter :: pi = 3.1415926535897931D+00
  real ( kind = rk ) sigma1
  real ( kind = rk ) sigma2
  real ( kind = rk ) theta

  theta = real ( 2 * mod ( deg, 2 ) - 1, kind = rk ) &
    * real ( family - 1, kind = rk ) * pi / 2.0D+00
  isigm1 = ( 2.0D+00 * sigma1 - ( b1 + a1 ) ) / ( b1 - a1 )
  isigm2 = ( 2.0D+00 * sigma2 - ( b2 + a2 ) ) / ( b2 - a2 )
  isigm1 = isigm1 * cos ( theta ) + isigm2 * sin ( theta )

  return
end
function sigma2 ( t1, t2, a1, a2, b1, b2, family, deg )

!*****************************************************************************80
!
!! SIGMA2 maps second coordinate from square to the rectangle.
!
!  Discussion:
!
!    This function returns the second component of the map 
!    from the square [-1,1]^2 to the rectangle [A1,B1] x [A2,B2]. 
!    FAMILY and DEG select the rotation in order to get 
!    the corresponding FAMILY of Padua points at degree DEG.
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
!    Input, real ( kind = rk ) A1, B1, A2, B2, the coordinates of the extreme
!    corners of the rectangle.
!
!    Input, integer FAMILY, DEG, select the family of Padua points at 
!    degree DEG.
!
!    Output, real ( kind = rk ) SIGMA2, the Y coordinate of the corresponding
!    point in the rectangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a1
  real ( kind = rk ) a2
  real ( kind = rk ) b1
  real ( kind = rk ) b2
  integer deg
  integer family
  real ( kind = rk ), parameter :: pi = 3.1415926535897931D+00
  real ( kind = rk ) sigma2
  real ( kind = rk ) t1
  real ( kind = rk ) t2
  real ( kind = rk ) theta

  theta = real ( 2 * mod ( deg, 2 ) - 1, kind = rk ) &
    * real ( family - 1, kind = rk ) * pi / 2.0D+00
  sigma2 = t1 * sin ( theta ) + t2 * cos ( theta )
  sigma2 = ( ( b2 - a2 ) * sigma2 + ( b2 + a2 ) ) / 2.0D+00

  return
end
function isigm2 ( sigma1, sigma2, a1, a2, b1, b2, family, deg )

!*****************************************************************************80
!
!! ISIGM2 maps the second coordinate from the rectangle to the square.
!
!  Discussion:
!
!    This function returns the second component of the map 
!    from the rectangle [A1,B1] x [A2,B2] to the square [-1,1]^2. 
!
!    FAMILY and DEG select the rotation in order to get 
!    the corresponding FAMILY of Padua points at degree DEG.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!  
!  Modified:
!
!    11 February 2014
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
!    Input, real ( kind = rk ) A1, B1, A2, B2, the coordinates of the extreme
!    corners of the rectangle.
!
!    Input, integer FAMILY, DEG, select the family of Padua points at 
!    degree DEG.
!
!    Output, real ( kind = rk ) ISIGM2, the Y coordinate of the corresponding
!    point in the rectangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a1
  real ( kind = rk ) a2
  real ( kind = rk ) b1
  real ( kind = rk ) b2
  integer deg
  integer family
  real ( kind = rk ) isigm1
  real ( kind = rk ) isigm2
  real ( kind = rk ), parameter :: pi = 3.1415926535897931D+00
  real ( kind = rk ) sigma1
  real ( kind = rk ) sigma2
  real ( kind = rk ) theta

  theta = real ( 2 * mod ( deg, 2 ) - 1, kind = rk ) &
    * real ( family - 1, kind = rk ) * pi / 2.0D+00
  isigm1 = ( 2.0D+00 * sigma1 - ( b1 + a1 ) ) / ( b1 - a1 )
  isigm2 = ( 2.0D+00 * sigma2 - ( b2 + a2 ) ) / ( b2 - a2 )
  isigm2 = - isigm1 * sin ( theta ) + isigm2 * cos ( theta )

  return
end
subroutine target ( a1, b1, a2, b2, ntg1, ntgmax, tg1, tg2, ntg )

!*****************************************************************************80
!
!! TARGET returns the target points on the rectangle.
!
!  Discussion:
!
!    Target points (uniform grid) on the rectangle [A1,B1] x [A2,B2].
!    The number of target points is NTG = NTG1^2.
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
!    Input, real ( kind = rk ) A1, B1, A2, B2, the coordinates of the extreme
!    corners of the rectangle.
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

  real ( kind = rk ) a1
  real ( kind = rk ) b1
  real ( kind = rk ) a2
  real ( kind = rk ) b2
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

  if ( ntgmax < ntg1 ** 2 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'TARGET - Fatal error!'
    write ( *, '(a)' ) '  NTGMAX < NTG1 * NTG1.'
    write ( *, '(a,i4)' ) '  NTG1 = ', ntg1
    write ( *, '(a,i4)' ) '  NTGMAX = ', ntgmax
    stop 1
  end if
  
  ntg = 0
  do i = 1, ntg1
    do j = 1, ntg1
      ntg = ntg + 1
      tg1(ntg) = a1 + real ( j - 1, kind = rk ) * ( b1 - a1 ) &
        / real ( ntg1 - 1, kind = rk ) 
      tg2(ntg) = a2 + real ( i - 1, kind = rk ) * ( b2 - a2 ) &
        / real ( ntg1 - 1, kind = rk )
    end do
  end do

  return
end
