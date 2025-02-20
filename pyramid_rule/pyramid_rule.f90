program main

!*****************************************************************************80
!
!! pyramid_rule() computes a quadrature rule for a pyramid.
!
!  Discussion:
!
!    This program computes a quadrature rule for a pyramid
!    and writes it to a file.
!
!    The user specifies:
!    * the LEGENDRE_ORDER (number of points in the X and Y dimensions)
!    * the JACOBI_ORDER (number of points in the Z dimension)
!    * FILENAME, the root name of the output files.
!
!    The integration region is:
!
!      - ( 1 - Z ) <= X <= 1 - Z
!      - ( 1 - Z ) <= Y <= 1 - Z
!                0 <= Z <= 1.
!
!    When Z is zero, the integration region is a square lying in the (X,Y) 
!    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
!    radius of the square diminishes, and when Z reaches 1, the square has 
!    contracted to the single point (0,0,1).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 July 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer arg_num
  character ( len = 255 ) filename
  integer iarg
  integer iargc
  integer ierror
  integer jacobi_order
  integer last
  integer legendre_order
  character ( len = 255 ) string

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'pyramid_rule():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Compute a quadrature rule for approximating'
  write ( *, '(a)' ) '  the integral of a function over a pyramid.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The user specifies:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  LEGENDRE_ORDER, the order of the Legendre rule for X and Y.'
  write ( *, '(a)' ) '  JACOBI_ORDER, the order of the Jacobi rule for Z,'
  write ( *, '(a)' ) '  FILENAME, the prefix of the three output files:'
  write ( *, '(a)' ) '    filename_w.txt - the weight file'
  write ( *, '(a)' ) '    filename_x.txt - the abscissa file.'
  write ( *, '(a)' ) '    filename_r.txt - the region file.'
!
!  Get the number of command line arguments.
!
  arg_num = iargc ( )
!
!  Get the Legendre order.
!
  if ( 1 <= arg_num ) then
  
    iarg = 1
    call getarg ( iarg, string )
    call s_to_i4 ( string, legendre_order, ierror, last )
    
  else
  
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter the Legendre rule order:'
    read ( *, * ) legendre_order
    
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The requested Legendre order of the rule is ', &
    legendre_order
!
!  Get the Jacobi order.
!
  if ( 2 <= arg_num ) then
  
    iarg = 2
    call getarg ( iarg, string )
    call s_to_i4 ( string, jacobi_order, ierror, last )
    
  else
  
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter the Jacobi rule order:'
    read ( *, * ) jacobi_order
    
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The requested Jacobi order of the rule is ', &
    jacobi_order
!
!  Get the output option or quadrature file root name:
!
  if ( 3 <= arg_num ) then

    iarg = 3
    call getarg ( iarg, filename )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter the "root name" of the quadrature files.'

    read ( *, '(a)' ) filename

  end if

  call pyramid_handle ( legendre_order, jacobi_order, filename )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PYRAMID_RULE:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 January 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IUNIT, the free unit number.
!
  implicit none

  integer i
  integer ios
  integer iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
subroutine jacobi_compute ( order, alpha, beta, xtab, weight )

!*****************************************************************************80
!
!! JACOBI_COMPUTE computes a Gauss-Jacobi quadrature rule.
!
!  Discussion:
!
!    The weight function is w(x) = (1-X)^ALPHA * (1+X)^BETA.
!
!    The integral to approximate:
!
!      Integral ( -1 <= X <= 1 ) (1-X)^ALPHA * (1+X)^BETA * F(X) dX
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
!
!    Thanks to Xu Xiang of Fudan University for pointing out that
!    an earlier implementation of this routine was incorrect!
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 May 2007
!
!  Author:
!
!    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input, integer ORDER, the order of the quadrature rule
!    to be computed.
!
!    Input, real ( kind = rk ) ALPHA, BETA, the exponents of (1-X) and
!    (1+X) in the quadrature rule.  For simple Gauss-Legendre quadrature,
!    set ALPHA = BETA = 0.0.  -1.0 < ALPHA and -1.0 < BETA are required.
!
!    Output, real ( kind = rk ) XTAB(ORDER), the abscissas.
!
!    Output, real ( kind = rk ) WEIGHT(ORDER), the weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer order

  real ( kind = rk ) alpha
  real ( kind = rk ) an
  real ( kind = rk ) b(order)
  real ( kind = rk ) beta
  real ( kind = rk ) bn
  real ( kind = rk ) c(order)
  real ( kind = rk ) cc
  real ( kind = rk ) delta
  real ( kind = rk ) dp2
  integer i
  real ( kind = rk ) p1
  real ( kind = rk ) r1
  real ( kind = rk ) r2
  real ( kind = rk ) r3
  real ( kind = rk ) weight(order)
  real ( kind = rk ) x
  real ( kind = rk ) xtab(order)
!
!  Check ALPHA and BETA.
!
  if ( alpha <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'JACOBI_COMPUTE - Fatal error!'
    write ( *, '(a)' ) '  -1.0 < ALPHA is required.'
    stop 1
  end if

  if ( beta <= -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'JACOBI_COMPUTE - Fatal error!'
    write ( *, '(a)' ) '  -1.0 < BETA is required.'
    stop 1
  end if
!
!  Set the recursion coefficients.
!
  do i = 1, order

    if ( alpha + beta == 0.0D+00 .or. beta - alpha == 0.0D+00 ) then

      b(i) = 0.0D+00

    else

      b(i) = ( alpha + beta ) * ( beta - alpha ) / &
            ( ( alpha + beta + real ( 2 * i, kind = rk ) ) &
            * ( alpha + beta + real ( 2 * i - 2, kind = rk ) ) )

    end if

    if ( i == 1 ) then

      c(i) = 0.0D+00

    else

      c(i) = 4.0D+00 * real ( i - 1, kind = rk ) &
            * ( alpha + real ( i - 1, kind = rk ) ) &
            * ( beta + real ( i - 1, kind = rk ) ) &
            * ( alpha + beta + real ( i - 1, kind = rk ) ) / &
            ( ( alpha + beta + real ( 2 * i - 1, kind = rk ) ) &
            * ( alpha + beta + real ( 2 * i - 2, kind = rk ) )**2 &
            * ( alpha + beta + real ( 2 * i - 3, kind = rk ) ) )

    end if

  end do

  delta = gamma ( alpha        + 1.0D+00 ) &
        * gamma (         beta + 1.0D+00 ) &
        / gamma ( alpha + beta + 2.0D+00 )

  cc = delta * 2.0D+00**( alpha + beta + 1.0D+00 ) * product ( c(2:order) )

  do i = 1, order

    if ( i == 1 ) then

      an = alpha / real ( order, kind = rk )
      bn = beta / real ( order, kind = rk )

      r1 = ( 1.0D+00 + alpha ) &
        * ( 2.78D+00 / ( 4.0D+00 + real ( order**2, kind = rk ) ) &
        + 0.768D+00 * an / real ( order, kind = rk ) )

      r2 = 1.0D+00 + 1.48D+00 * an + 0.96D+00 * bn &
        + 0.452D+00 * an**2 + 0.83D+00 * an * bn

      x = ( r2 - r1 ) / r2

    else if ( i == 2 ) then

      r1 = ( 4.1D+00 + alpha ) / &
        ( ( 1.0D+00 + alpha ) * ( 1.0D+00 + 0.156D+00 * alpha ) )

      r2 = 1.0D+00 + 0.06D+00 * ( real ( order, kind = rk ) - 8.0D+00 ) * &
        ( 1.0D+00 + 0.12D+00 * alpha ) / real ( order, kind = rk )

      r3 = 1.0D+00 + 0.012D+00 * beta * &
        ( 1.0D+00 + 0.25D+00 * abs ( alpha ) ) / real ( order, kind = rk )

      x = x - r1 * r2 * r3 * ( 1.0D+00 - x )

    else if ( i == 3 ) then

      r1 = ( 1.67D+00 + 0.28D+00 * alpha ) / ( 1.0D+00 + 0.37D+00 * alpha )

      r2 = 1.0D+00 + 0.22D+00 * ( real ( order, kind = rk ) - 8.0D+00 ) &
        / real ( order, kind = rk )

      r3 = 1.0D+00 + 8.0D+00 * beta / &
        ( ( 6.28D+00 + beta ) * real ( order**2, kind = rk ) )

      x = x - r1 * r2 * r3 * ( xtab(1) - x )

    else if ( i < order - 1 ) then

      x = 3.0D+00 * xtab(i-1) - 3.0D+00 * xtab(i-2) + xtab(i-3)

    else if ( i == order - 1 ) then

      r1 = ( 1.0D+00 + 0.235D+00 * beta ) / ( 0.766D+00 + 0.119D+00 * beta )

      r2 = 1.0D+00 / ( 1.0D+00 + 0.639D+00 &
        * ( real ( order, kind = rk ) - 4.0D+00 ) &
        / ( 1.0D+00 + 0.71D+00 * ( real ( order, kind = rk ) - 4.0D+00 ) ) )

      r3 = 1.0D+00 / ( 1.0D+00 + 20.0D+00 * alpha / ( ( 7.5D+00 + alpha ) * &
        real ( order**2, kind = rk ) ) )

      x = x + r1 * r2 * r3 * ( x - xtab(i-2) )

    else if ( i == order ) then

      r1 = ( 1.0D+00 + 0.37D+00 * beta ) / ( 1.67D+00 + 0.28D+00 * beta )

      r2 = 1.0D+00 / &
        ( 1.0D+00 + 0.22D+00 * ( real ( order, kind = rk ) - 8.0D+00 ) &
        / real ( order, kind = rk ) )

      r3 = 1.0D+00 / ( 1.0D+00 + 8.0D+00 * alpha / &
        ( ( 6.28D+00 + alpha ) * real ( order**2, kind = rk ) ) )

      x = x + r1 * r2 * r3 * ( x - xtab(i-2) )

    end if

    call jacobi_root ( x, order, alpha, beta, dp2, p1, b, c )

    xtab(i) = x
    weight(i) = cc / ( dp2 * p1 )

  end do
!
!  Reverse the order of the data.
!
  xtab(1:order) = xtab(order:1:-1)
  weight(1:order) = weight(order:1:-1)

  return
end
subroutine jacobi_recur ( p2, dp2, p1, x, order, alpha, beta, b, c )

!*****************************************************************************80
!
!! JACOBI_RECUR finds the value and derivative of a Jacobi polynomial.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Output, real ( kind = rk ) P2, the value of J(ORDER)(X).
!
!    Output, real ( kind = rk ) DP2, the value of J'(ORDER)(X).
!
!    Output, real ( kind = rk ) P1, the value of J(ORDER-1)(X).
!
!    Input, real ( kind = rk ) X, the point at which polynomials are evaluated.
!
!    Input, integer ORDER, the order of the polynomial 
!    to be computed.
!
!    Input, real ( kind = rk ) ALPHA, BETA, the exponents of (1-X) and
!    (1+X) in the quadrature rule.
!
!    Input, real ( kind = rk ) B(ORDER), C(ORDER), the recursion
!    coefficients.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer order

  real ( kind = rk ) alpha
  real ( kind = rk ) b(order)
  real ( kind = rk ) beta
  real ( kind = rk ) c(order)
  real ( kind = rk ) dp0
  real ( kind = rk ) dp1
  real ( kind = rk ) dp2
  integer i
  real ( kind = rk ) p0
  real ( kind = rk ) p1
  real ( kind = rk ) p2
  real ( kind = rk ) x

  p1 = 1.0D+00
  dp1 = 0.0D+00

  p2 = x + ( alpha - beta ) / ( alpha + beta + 2.0D+00 )
  dp2 = 1.0D+00

  do i = 2, order

    p0 = p1
    dp0 = dp1

    p1 = p2
    dp1 = dp2

    p2 = ( x - b(i) ) * p1 - c(i) * p0
    dp2 = ( x - b(i) ) * dp1 + p1 - c(i) * dp0

  end do

  return
end
subroutine jacobi_root ( x, order, alpha, beta, dp2, p1, b, c )

!*****************************************************************************80
!
!! JACOBI_ROOT improves an approximate root of a Jacobi polynomial.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input/output, real ( kind = rk ) X, the approximate root, which
!    should be improved on output.
!
!    Input, integer ORDER, the order of the polynomial 
!    to be computed.
!
!    Input, real ( kind = rk ) ALPHA, BETA, the exponents of (1-X) and
!    (1+X) in the quadrature rule.
!
!    Output, real ( kind = rk ) DP2, the value of J'(ORDER)(X).
!
!    Output, real ( kind = rk ) P1, the value of J(ORDER-1)(X).
!
!    Input, real ( kind = rk ) B(ORDER), C(ORDER), the recursion coefficients.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer order

  real ( kind = rk ) alpha
  real ( kind = rk ) b(order)
  real ( kind = rk ) beta
  real ( kind = rk ) c(order)
  real ( kind = rk ) d
  real ( kind = rk ) dp2
  real ( kind = rk ) eps
  real ( kind = rk ) p1
  real ( kind = rk ) p2
  real ( kind = rk ) r8_epsilon
  integer step
  integer, parameter :: step_max = 10
  real ( kind = rk ) x

  eps = r8_epsilon ( )

  do step = 1, step_max

    call jacobi_recur ( p2, dp2, p1, x, order, alpha, beta, b, c )

    d = p2 / dp2
    x = x - d

    if ( abs ( d ) <= eps * ( abs ( x ) + 1.0D+00 ) ) then
      return
    end if

  end do

  return
end
subroutine legendre_compute ( order, xtab, weight )

!*****************************************************************************80
!
!! LEGENDRE_COMPUTE computes a Gauss-Legendre quadrature rule.
!
!  Discussion:
!
!    The integration interval is [ -1, 1 ].
!
!    The weight function is w(x) = 1.0.
!
!    The integral to approximate:
!
!      Integral ( -1 <= X <= 1 ) F(X) dX
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ORDER, the order of the rule.
!    ORDER must be greater than 0.
!
!    Output, real ( kind = rk ) XTAB(ORDER), the abscissas.
!
!    Output, real ( kind = rk ) WEIGHT(ORDER), the weights.
!    The weights are positive, symmetric, and should sum to 2.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer order

  real ( kind = rk ) d1
  real ( kind = rk ) d2pn
  real ( kind = rk ) d3pn
  real ( kind = rk ) d4pn
  real ( kind = rk ) dp
  real ( kind = rk ) dpn
  real ( kind = rk ) e1
  real ( kind = rk ) fx
  real ( kind = rk ) h
  integer i
  integer iback
  integer k
  integer m
  integer mp1mi
  integer ncopy
  integer nmove
  real ( kind = rk ) p
  real ( kind = rk ) :: pi = 3.141592653589793D+00
  real ( kind = rk ) pk
  real ( kind = rk ) pkm1
  real ( kind = rk ) pkp1
  real ( kind = rk ) t
  real ( kind = rk ) u
  real ( kind = rk ) v
  real ( kind = rk ) x0
  real ( kind = rk ) xtab(order)
  real ( kind = rk ) xtemp
  real ( kind = rk ) weight(order)

  if ( order < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_COMPUTE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of ORDER = ', order
    stop 1
  end if

  e1 = real ( order * ( order + 1 ), kind = rk )

  m = ( order + 1 ) / 2

  do i = 1, m

    mp1mi = m + 1 - i

    t = real ( 4 * i - 1, kind = rk ) * pi &
      / real ( 4 * order + 2, kind = rk )

    x0 = cos ( t ) * ( 1.0D+00 - ( 1.0D+00 - 1.0D+00 &
      / real ( order, kind = rk ) ) &
      / real ( 8 * order * order, kind = rk ) )

    pkm1 = 1.0D+00
    pk = x0

    do k = 2, order
      pkp1 = 2.0D+00 * x0 * pk - pkm1 - ( x0 * pk - pkm1 ) &
        / real ( k, kind = rk )
      pkm1 = pk
      pk = pkp1
    end do

    d1 = real ( order, kind = rk ) * ( pkm1 - x0 * pk )

    dpn = d1 / ( 1.0D+00 - x0 * x0 )

    d2pn = ( 2.0D+00 * x0 * dpn - e1 * pk ) / ( 1.0D+00 - x0 * x0 )

    d3pn = ( 4.0D+00 * x0 * d2pn + ( 2.0D+00 - e1 ) * dpn ) &
      / ( 1.0D+00 - x0 * x0 )

    d4pn = ( 6.0D+00 * x0 * d3pn + ( 6.0D+00 - e1 ) * d2pn ) / &
      ( 1.0D+00 - x0 * x0 )

    u = pk / dpn
    v = d2pn / dpn
!
!  Initial approximation H:
!
    h = - u * ( 1.0D+00 + 0.5D+00 * u * ( v + u * ( v * v - d3pn / &
      ( 3.0D+00 * dpn ) ) ) )
!
!  Refine H using one step of Newton's method:
!
    p = pk + h * ( dpn + 0.5D+00 * h * ( d2pn + h / 3.0D+00 &
      * ( d3pn + 0.25D+00 * h * d4pn ) ) )

    dp = dpn + h * ( d2pn + 0.5D+00 * h * ( d3pn + h * d4pn / 3.0D+00 ) )

    h = h - p / dp

    xtemp = x0 + h

    xtab(mp1mi) = xtemp

    fx = d1 - h * e1 * ( pk + 0.5D+00 * h * ( dpn + h / 3.0D+00 &
      * ( d2pn + 0.25D+00 * h * ( d3pn + 0.2D+00 * h * d4pn ) ) ) )

    weight(mp1mi) = 2.0D+00 * ( 1.0D+00 - xtemp * xtemp ) / ( fx * fx )

  end do

  if ( mod ( order, 2 ) == 1 ) then
    xtab(1) = 0.0D+00
  end if
!
!  Shift the data up.
!
  nmove = ( order + 1 ) / 2
  ncopy = order - nmove

  do i = 1, nmove
    iback = order + 1 - i
    xtab(iback) = xtab(iback-ncopy)
    weight(iback) = weight(iback-ncopy)
  end do
!
!  Reflect values for the negative abscissas.
!
  do i = 1, order - nmove
    xtab(i) = - xtab(order+1-i)
    weight(i) = weight(order+1-i)
  end do

  return
end
subroutine pyramid_handle ( legendre_order, jacobi_order, filename )

!*****************************************************************************80
!
!! PYRAMID_HANDLE computes the requested pyramid rule and outputs it.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 July 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LEGENDRE_ORDER, JACOBI_ORDER, the orders 
!    of the component Legendre and Jacobi rules.
!
!    Input, character ( len = * ) FILENAME, the rootname for the files,  
!    write files 'file_w.txt' and 'file_x.txt', and 'file_r.txt', weights,
!    abscissas, and region.
! 
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  character ( len = * )  filename
  character ( len = 255 ) filename_r
  character ( len = 255 ) filename_w
  character ( len = 255 ) filename_x
  integer i
  integer j
  real ( kind = rk ) jacobi_alpha
  real ( kind = rk ) jacobi_beta
  integer jacobi_order
  real ( kind = rk ), allocatable, dimension ( : ) :: jacobi_w
  real ( kind = rk ), allocatable, dimension ( : ) :: jacobi_x
  integer k
  integer l
  integer legendre_order
  real ( kind = rk ), allocatable, dimension ( : ) :: legendre_w
  real ( kind = rk ), allocatable, dimension ( : ) :: legendre_x
  integer pyramid_order
  real ( kind = rk ) pyramid_r(dim_num,5)
  real ( kind = rk ), allocatable, dimension ( : ) :: pyramid_w
  real ( kind = rk ), allocatable, dimension ( :, : ) :: pyramid_x
  real ( kind = rk ) volume
  real ( kind = rk ) wi
  real ( kind = rk ) wj
  real ( kind = rk ) wk
  real ( kind = rk ) xi
  real ( kind = rk ) xj
  real ( kind = rk ) xk
!
!  Compute the factor rules.
!
  allocate ( legendre_w(legendre_order) )
  allocate ( legendre_x(legendre_order) )

  call legendre_compute ( legendre_order, legendre_x, legendre_w )

  allocate ( jacobi_w(jacobi_order) )
  allocate ( jacobi_x(jacobi_order) )

  jacobi_alpha = 2.0D+00
  jacobi_beta = 0.0D+00

  call jacobi_compute ( jacobi_order, jacobi_alpha, jacobi_beta, jacobi_x, &
    jacobi_w )
!
!  Compute the pyramid rule.
!
  pyramid_order = legendre_order * legendre_order * jacobi_order

  allocate ( pyramid_w(pyramid_order) )
  allocate ( pyramid_x(1:dim_num,pyramid_order) )

  volume = 4.0D+00 / 3.0D+00

  l = 0
  do k = 1, jacobi_order
    xk = ( jacobi_x(k) + 1.0D+00 ) / 2.0D+00
    wk = jacobi_w(k) / 2.0D+00
    do j = 1, legendre_order
      xj = legendre_x(j)
      wj = legendre_w(j)
      do i = 1, legendre_order
        xi = legendre_x(i)
        wi = legendre_w(i)
        l = l + 1
        pyramid_w(l) = wi * wj * wk / 4.0D+00 / volume
        pyramid_x(1:dim_num,l) = (/ xi * ( 1.0D+00 - xk ), &
          xj * ( 1.0D+00 - xk ), xk /)
      end do
    end do
  end do

  deallocate ( jacobi_w )
  deallocate ( jacobi_x )
  deallocate ( legendre_w )
  deallocate ( legendre_x )

  pyramid_r(1:dim_num,1) = (/ -1.0D+00, -1.0D+00, 0.0D+00 /)
  pyramid_r(1:dim_num,2) = (/ +1.0D+00, -1.0D+00, 0.0D+00 /)
  pyramid_r(1:dim_num,3) = (/ -1.0D+00, +1.0D+00, 0.0D+00 /)
  pyramid_r(1:dim_num,4) = (/ +1.0D+00, +1.0D+00, 0.0D+00 /)
  pyramid_r(1:dim_num,5) = (/  0.0D+00,  0.0D+00, 1.0D+00 /)
!
!  Write the rule to files.
!
  filename_w = trim ( filename ) // '_w.txt'
  filename_x = trim ( filename ) // '_x.txt'
  filename_r = trim ( filename ) // '_r.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Creating quadrature files.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  "Root" file name is   "' // trim ( filename ) // '".'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Weight file will be   "' // trim ( filename_w ) // '".'
  write ( *, '(a)' ) '  Abscissa file will be "' // trim ( filename_x ) // '".'
  write ( *, '(a)' ) '  Region file will be   "' // trim ( filename_r ) // '".'

  call r8mat_write ( filename_w, 1,       pyramid_order, pyramid_w )
  call r8mat_write ( filename_x, dim_num, pyramid_order, pyramid_x )
  call r8mat_write ( filename_r, dim_num, 5,             pyramid_r )

  deallocate ( pyramid_w )
  deallocate ( pyramid_x )

  return
end
function r8_epsilon ( )

!*****************************************************************************80
!
!! R8_EPSILON returns the R8 roundoff unit.
!
!  Discussion:
!
!    The roundoff unit is a number R which is a power of 2 with the
!    property that, to the precision of the computer's arithmetic,
!      1 < 1 + R
!    but
!      1 = ( 1 + R / 2 )
!
!    FORTRAN90 provides the superior library routine
!
!      EPSILON ( X )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) R8_EPSILON, the round-off unit.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r8_epsilon

  r8_epsilon = 2.220446049250313D-016

  return
end
subroutine r8mat_write ( output_filename, m, n, table )

!*****************************************************************************80
!
!! R8MAT_WRITE writes an R8MAT file.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the output file name.
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) TABLE(M,N), the table data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer j
  character ( len = * ) output_filename
  integer output_status
  integer output_unit
  character ( len = 30 ) string
  real ( kind = rk ) table(m,n)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_filename ) // '" on unit ', output_unit
    output_unit = -1
    stop 1
  end if
!
!  Create a format string.
!
!  For less precision in the output file, try:
!
!                                            '(', m, 'g', 14, '.', 6, ')'
!
  if ( 0 < m .and. 0 < n ) then

    write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 24, '.', 16, ')'
!
!  Write the data.
!
    do j = 1, n
      write ( output_unit, string ) table(1:m,j)
    end do

  end if
!
!  Close the file.
!
  close ( unit = output_unit )

  return
end
subroutine s_to_i4 ( s, ival, ierror, length )

!*****************************************************************************80
!
!! S_TO_I4 reads an I4 from a string.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 January 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string to be examined.
!
!    Output, integer IVAL, the integer value read from the string.
!    If the string is blank, then IVAL will be returned 0.
!
!    Output, integer IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer LENGTH, the number of characters of S 
!    used to make IVAL.
!
  implicit none

  character c
  integer i
  integer ierror
  integer isgn
  integer istate
  integer ival
  integer length
  character ( len = * ) s

  ierror = 0
  istate = 0
  isgn = 1
  ival = 0

  do i = 1, len_trim ( s )

    c = s(i:i)
!
!  Haven't read anything.
!
    if ( istate == 0 ) then

      if ( c == ' ' ) then

      else if ( c == '-' ) then
        istate = 1
        isgn = -1
      else if ( c == '+' ) then
        istate = 1
        isgn = + 1
      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read the sign, expecting digits.
!
    else if ( istate == 1 ) then

      if ( c == ' ' ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read at least one digit, expecting more.
!
    else if ( istate == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        ival = 10 * ival + ichar ( c ) - ichar ( '0' )
      else
        ival = isgn * ival
        length = i - 1
        return
      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( istate == 2 ) then
    ival = isgn * ival
    length = len_trim ( s )
  else
    ierror = 1
    length = 0
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

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
