subroutine diffusivity_1d_xk ( dc0, m, omega, n, x, dc )

!*****************************************************************************80
!
!! diffusivity_1d_xk() evaluates a 1D stochastic diffusivity function.
!
!  Discussion:
!
!    The 1D diffusion equation has the form
!
!      - d/dx ( DC(X) Del U(X) ) = F(X)
!
!    where DC(X) is a function called the diffusivity.
!
!    In the stochastic version of the problem, the diffusivity function
!    includes the influence of stochastic parameters:
!
!      - d/dx ( DC(X;OMEGA) d/dx U(X) ) = F(X).
!
!    In this function, the domain is assumed to be the unit interval [0.1].
!
!
!    For DC0 = 1 and F(X) = 0, with boundary conditions U(0:OMEGA) = 0,
!    U(1;OMEGA) = 1, the exact solution is
!
!    If OMEGA ~= 0:
!
!      U(X;OMEGA) = log ( 1 + OMEGA * X ) / log ( 1 + OMEGA )
!
!    If OMEGA = 0:
!
!      U(X;OMEGA) = X
!
!    In the numerical experiments described in the paper, OMEGA was taken
!    to be a random variable with a Beta, or Uniform, or Gaussian or 
!    Poisson or Binomial distribution.
!
!    For the Gaussian and Poisson distributions, the positivity requirement 
!    could not be guaranteed, and the experiments were simply made with a 
!    "small" variance of 0.1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 July 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu, George Karniadakis,
!    Modeling uncertainty in steady state diffusion problems via
!    generalized polynomial chaos,
!    Computer Methods in Applied Mechanics and Engineering,
!    Volume 191, 2002, pages 4927-4948.
!
!  Parameters:
!
!    Input, real ( kind = rk8 ) DC0, the constant term in the expansion of the 
!    diffusion coefficient.
!
!    Input, integer M, the number of stochastic parameters.
!
!    Input, real ( kind = rk8 ) OMEGA(M), the stochastic parameters.  
!
!    Input, integer N, the number of evaluation points.
!
!    Input, real ( kind = rk8 ) X(N), the point where the diffusion coefficient 
!    is to be evaluated.
!
!    Output, real ( kind = rk8 ) DC(N), the value of the diffusion coefficient 
!    at X.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk8 ) dc(n)
  real ( kind = rk8 ) dc0
  integer k
  real ( kind = rk8 ) omega(m)
  real ( kind = rk8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk8 ) w
  real ( kind = rk8 ) x(n)

  k = 0
  w = 1.0D+00

  dc(1:n) = 0.0D+00

  do while ( k < m )

    if ( k < m ) then
      k = k + 1
      dc(1:n) = dc(1:n) + omega(k) * sin ( w * pi * x(1:n) )
    end if

    if ( k < m ) then
      k = k + 1
      dc(1:n) = dc(1:n) + omega(k) * cos ( w * pi * x(1:n) )    
    end if

    w = w + 1.0D+00

  end do

  dc(1:n) = exp ( - 0.125D+00 ) * dc(1:n)

  dc(1:n) = dc0 + exp ( dc(1:n) )

  return
end
subroutine diffusivity_2d_bnt ( dc0, omega, n, x, y, dc )

!*****************************************************************************80
!
!! diffusivity_2d_bnt() evaluates a 2D stochastic diffusivity function.
!
!  Discussion:
!
!    The 2D diffusion equation has the form
!
!      - Del ( DC(X,Y) Del U(X,Y) ) = F(X,Y)
!
!    where DC(X,Y) is a function called the diffusivity.
!
!    In the stochastic version of the problem, the diffusivity function
!    includes the influence of stochastic parameters:
!
!      - Del ( DC(X,Y;OMEGA) Del U(X,Y;OMEGA) ) = F(X,Y).
!
!    In this function, the domain is the rectangle [-1.5,0]x[-0.4,0.8].
!
!    The four stochastic parameters OMEGA(1:4) are assumed to be independent
!    identically distributed random variables with mean value zero and 
!    variance 1.  The distribution is typically taken to be Gaussian or
!    uniform.
!
!    A collocation approach to this problem would then use the roots of
!    Hermite or Legendre polynomials.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 July 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Ivo Babuska, Fabio Nobile, Raul Tempone,
!    A stochastic collocation method for elliptic partial differential equations
!    with random input data,
!    SIAM Journal on Numerical Analysis,
!    Volume 45, Number 3, 2007, pages 1005-1034.
!
!  Parameters:
!
!    Input, real ( kind = rk8 ) DC0, the constant term in the expansion of the 
!    diffusion coefficient.  Take DC0 = 10.
!
!    Input, real ( kind = rk8 ) OMEGA(4), the stochastic parameters.
!
!    Input, integer N, the number of evaluation points.
!
!    Input, real ( kind = rk8 ) X(N), Y(N), the points where the diffusion 
!    coefficient is to be evaluated.
!
!    Output, real ( kind = rk8 ) DC(N), the value of the diffusion coefficient 
!    at (X,Y).
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  real ( kind = rk8 ) arg(n)
  real ( kind = rk8 ) dc(n)
  real ( kind = rk8 ) dc0
  real ( kind = rk8 ) omega(4)
  real ( kind = rk8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk8 ) x(n)
  real ( kind = rk8 ) y(n)

  arg(1:n) = omega(1) * cos ( pi * x(1:n) ) &
           + omega(2) * sin ( pi * x(1:n) ) &
           + omega(3) * cos ( pi * y(1:n) ) &
           + omega(4) * sin ( pi * y(1:n) )

  arg(1:n)= exp ( - 0.125D+00 ) * arg(1:n)

  dc(1:n) = dc0 + exp ( arg(1:n) )

  return
end
subroutine diffusivity_2d_elman ( a, cl, dc0, m_1d, omega, n1, n2, x, y, dc )

!*****************************************************************************80
!
!! diffusivity_2d_elman() evaluates a 2D stochastic diffusivity function.
!
!  Discussion:
!
!    The 2D diffusion equation has the form
!
!      - Del ( DC(X,Y) Del U(X,Y) ) = F(X,Y)
!
!    where DC(X,Y) is a function called the diffusivity.
!
!    In the stochastic version of the problem, the diffusivity function
!    includes the influence of stochastic parameters:
!
!      - Del ( DC(X,Y;OMEGA) Del U(X,Y;OMEGA) ) = F(X,Y).
!
!    In this function, the domain is assumed to be the square [-A,+A]x[-A,+A].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 July 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Howard Elman, Darran Furnaval,
!    Solving the stochastic steady-state diffusion problem using multigrid,
!    IMA Journal on Numerical Analysis,
!    Volume 27, Number 4, 2007, pages 675-688.
!
!    Roger Ghanem, Pol Spanos,
!    Stochastic Finite Elements: A Spectral Approach,
!    Revised Edition,
!    Dover, 2003,
!    ISBN: 0486428184,
!    LC: TA347.F5.G56.
!
!  Parameters:
!
!    Input, real ( kind = rk8 ) A, the "radius" of the square region.  The region
!    is assumed to be [-A,+A]x[-A,+A].
!    0 < A.
!
!    Input, real ( kind = rk8 ) CL, the correlation length.
!    0 < CL.
!
!    Input, real ( kind = rk8 ) DC0, the constant term in the expansion of the 
!    diffusion coefficient.  Take DC0 = 10.
!
!    Input, integer M_1D, the first and second dimensions of the
!    stochastic parameter array.
!
!    Input, real ( kind = rk8 ) OMEGA(M_1D,M_1D), the stochastic parameters.
!
!    Input, integer N1, N2, the dimensions of the X and Y arrays.
!
!    Input, real ( kind = rk8 ) X(N1,N2), Y(N1,N2), the points where the
!    diffusion coefficient is to be evaluated.
!
!    Output, real ( kind = rk8 ) DC(N1,N2), the value of the diffusion 
!    coefficient at X.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer m_1d
  integer n1
  integer n2

  real ( kind = rk8 ) a
  real ( kind = rk8 ) c_1dx(m_1d,n1,n2)
  real ( kind = rk8 ) c_1dy(m_1d,n1,n2)
  real ( kind = rk8 ) cl
  real ( kind = rk8 ) dc(n1,n2)
  real ( kind = rk8 ) dc0
  integer i
  integer j
  integer k
  real ( kind = rk8 ) lambda_1d(m_1d)
  integer m
  real ( kind = rk8 ) omega(m_1d,m_1d)
  real ( kind = rk8 ) theta_1d(m_1d)
  real ( kind = rk8 ) x(n1,n2)
  real ( kind = rk8 ) y(n1,n2)

  m = m_1d * m_1d
!
!  Compute THETA.
!
  call theta_solve ( a, cl, m_1d, theta_1d )
!
!  Compute LAMBDA_1D.
!
  lambda_1d(1:m_1d) = 2.0D+00 * cl &
    / ( 1.0D+00 + cl * cl * theta_1d(1:m_1d) ** 2 )
!
!  Compute C_1DX(1:M1D) and C_1DY(1:M1D) at (X,Y).
!
  c_1dx(1:m_1d,1:n1,1:n2) = 0.0D+00
  c_1dy(1:m_1d,1:n1,1:n2) = 0.0D+00

  k = 0

  do

    if ( m_1d <= k ) then
      exit
    end if

    k = k + 1

    c_1dx(k,1:n1,1:n2) = cos ( theta_1d(k) * a * x(1:n1,1:n2) ) &
      / sqrt ( a + sin ( 2.0D+00 * theta_1d(k) * a ) &
      / ( 2.0D+00 * theta_1d(k) ) )

    c_1dy(k,1:n1,1:n2) = cos ( theta_1d(k) * a * y(1:n1,1:n2) ) &
      / sqrt ( a + sin ( 2.0D+00 * theta_1d(k) * a ) &
      / ( 2.0D+00 * theta_1d(k) ) )

    if ( m_1d <= k ) then
      exit
    end if

    k = k + 1

    c_1dx(k,1:n1,1:n2) = sin ( theta_1d(k) * a * x(1:n1,1:n2) ) &
      / sqrt ( a - sin ( 2.0D+00 * theta_1d(k) * a ) &
      / ( 2.0D+00 * theta_1d(k) ) )

    c_1dy(k,1:n1,1:n2) = sin ( theta_1d(k) * a * y(1:n1,1:n2) ) &
      / sqrt ( a - sin ( 2.0D+00 * theta_1d(k) * a ) &
      / ( 2.0D+00 * theta_1d(k) ) )

  end do
!
!  Evaluate the diffusion coefficient DC at (X,Y).
!
  dc(1:n1,1:n2) = dc0
  do j = 1, m_1d
    do i = 1, m_1d
      dc(1:n1,1:n2) = dc(1:n1,1:n2) + sqrt ( lambda_1d(i) * lambda_1d(j) ) &
        * c_1dx(i,1:n1,1:n2) * c_1dy(j,1:n1,1:n2) * omega(i,j)
    end do
  end do

  return
end
subroutine diffusivity_2d_ntw ( cl, dc0, m, omega, n, x, y, dc )

!*****************************************************************************80
!
!! diffusivity_2d_ntw() evaluates a 2D stochastic diffusivity function.
!
!  Discussion:
!
!    The 2D diffusion equation has the form
!
!      - Del ( DC(X,Y) Del U(X,Y) ) = F(X,Y)
!
!    where DC(X,Y) is a function called the diffusivity.
!
!    In the stochastic version of the problem, the diffusivity function
!    includes the influence of stochastic parameters:
!
!      - Del ( DC(X,Y;OMEGA) Del U(X,Y;OMEGA) ) = F(X,Y).
!
!    In this function, the domain is the rectangle [0,D]x[0,D] where D = 1.
!
!    Note that in this problem the diffusivity has a one-dimensional
!    spatial dependence on X, but not on Y!
!
!    The random variables OMEGA are independent, have zero mean and unit
!    variance, and are uniformly distributed in [-sqrt(3),+sqrt(3)].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 July 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Xiang Ma, Nicholas Zabaras,
!    An adaptive hierarchical sparse grid collocation algorithm for the solution
!    of stochastic differential equations,
!    Journal of Computational Physics,
!    Volume 228, pages 3084-3113, 2009.
!
!    Fabio Nobile, Raul Tempone, Clayton Webster,
!    A Sparse Grid Stochastic Collocation Method for Partial Differential
!    Equations with Random Input Data,
!    SIAM Journal on Numerical Analysis,
!    Volume 46, Number 5, 2008, pages 2309-2345.
!
!  Parameters:
!
!    Input, real ( kind = rk8 ) CL, the desired physical correlation length for 
!    the coefficient.
!
!    Input, real ( kind = rk8 ) DC0, the constant term in the expansion of the 
!    diffusion coefficient.  Take DC0 = 0.5.
!
!    Input, integer M, the number of terms in the expansion.
!
!    Input, real ( kind = rk8 ) OMEGA(M), the stochastic parameters.
!
!    Input, integer N, the number of evaluation points.
!
!    Input, real ( kind = rk8 ) X(N), Y(N), the points where the diffusion 
!    coefficient is to be evaluated.
!
!    Output, real ( kind = rk8 ) DC(N), the value of the diffusion coefficient
!    at (X,Y).
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk8 ) cl
  real ( kind = rk8 ) d
  real ( kind = rk8 ) dc(n)
  real ( kind = rk8 ) dc_arg(n)
  real ( kind = rk8 ) dc0
  integer i
  real ( kind = rk8 ) ihalf_r8
  real ( kind = rk8 ) l
  real ( kind = rk8 ) lp
  real ( kind = rk8 ) omega(m)
  real ( kind = rk8 ) phi(n)
  real ( kind = rk8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk8 ) x(n)
  real ( kind = rk8 ) y(n)
  real ( kind = rk8 ) zeta
  real ( kind = rk8 ) zeta_arg

  d = 1.0D+00
  lp = max ( d, 2.0D+00 * cl )
  l = cl / lp

  dc_arg(1:n) = 1.0D+00 + omega(1) * sqrt ( sqrt ( pi ) * l / 2.0D+00 )

  do i = 2, m

    ihalf_r8 = real ( i / 2, kind = rk8 )
    zeta_arg = - ( ihalf_r8 * pi * l ) ** 2 / 8.0D+00
    zeta = sqrt ( sqrt ( pi ) * l ) * exp ( zeta_arg )

    if ( mod ( i, 2 ) == 0 ) then
      phi(1:n) = sin ( ihalf_r8 * pi * x(1:n) / lp )
    else
      phi(1:n) = cos ( ihalf_r8 * pi * x(1:n) / lp )
    end if

    dc_arg(1:n) = dc_arg(1:n) + zeta * phi(1:n) * omega(i)

  end do

  dc(1:n) = dc0 + exp ( dc_arg(1:n) )

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! get_unit() returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 October 2008
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
function r8mat_max ( m, n, a )

!*****************************************************************************80
!
!! r8mat_max() returns the maximum entry of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows in A.
!
!    Input, integer N, the number of columns in A.
!
!    Input, real ( kind = rk8 ) A(M,N), the M by N matrix.
!
!    Output, real ( kind = rk8 ) R8MAT_MAX, the maximum entry of A.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk8 ) a(m,n)
  real ( kind = rk8 ) r8mat_max

  r8mat_max = maxval ( a(1:m,1:n) )

  return
end
subroutine r8vec_linspace ( n, a, b, x )

!*****************************************************************************80
!
!! r8vec_linspace() creates a vector of linearly spaced values.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    4 points evenly spaced between 0 and 12 will yield 0, 4, 8, 12.
!
!    In other words, the interval is divided into N-1 even subintervals,
!    and the endpoints of intervals are used as the points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real ( kind = rk8 ) A_FIRST, A_LAST, the first and last entries.
!
!    Output, real ( kind = rk8 ) X(N), a vector of linearly spaced data.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  integer i
  real ( kind = rk8 ) x(n)

  if ( n == 1 ) then

    x(1) = ( a + b ) / 2.0D+00

  else

    do i = 1, n
      x(i) = ( real ( n - i,     kind = rk8 ) * a   &
             + real (     i - 1, kind = rk8 ) * b ) &
             / real ( n     - 1, kind = rk8 )
    end do

  end if

  return
end
subroutine r8vec_max ( n, a, amax )

!*****************************************************************************80
!
!! r8vec_max() returns the maximum value in an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input, real ( kind = rk8 ) A(N), the array.
!
!    Output, real ( kind = rk8 ) AMAX, the value of the largest entry.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  real ( kind = rk8 ) a(n)
  real ( kind = rk8 ) amax

  amax = maxval ( a(1:n) )

  return
end
subroutine r8vec_mesh_2d ( nx, ny, xvec, yvec, xmat, ymat )

!*****************************************************************************80
!
!! r8vec_mesh_2d() creates a 2D mesh from X and Y vectors.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    NX = 2
!    XVEC = ( 1, 2, 3 )
!    NY = 3
!    YVEC = ( 4, 5 )
!
!    XMAT = (
!      1, 2, 3
!      1, 2, 3 )
!
!    YMAT = (
!      4, 4, 4
!      5, 5, 5 ) 
!    
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 July 2013
!
!  Parameters:
!
!    Input, integer NX, NY, the number of X and Y values.
!
!    Input, real ( kind = rk8 ) XVEC(NX), YVEC(NY), the X and Y coordinate
!    values.
!
!    Output, real ( kind = rk8 ) XMAT(NX,NY), YMAT(NX,NY), the coordinate
!    values of points on an NX by NY mesh.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer nx
  integer ny

  integer j
  real ( kind = rk8 ) xmat(nx,ny)
  real ( kind = rk8 ) xvec(nx)
  real ( kind = rk8 ) ymat(nx,ny)
  real ( kind = rk8 ) yvec(ny)

  do j = 1, ny
    xmat(1:nx,j) = xvec(1:nx)
  end do

  do j = 1, ny
    ymat(1:nx,j) = yvec(j)
  end do

 return
end
subroutine r8vec_normal_01 ( n, x )

!*****************************************************************************80
!
!! r8vec_normal_01() returns a unit pseudonormal R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!    This routine can generate a vector of values on one call.  It
!    has the feature that it should provide the same results
!    in the same order no matter how we break up the task.
!
!    The Box-Muller method is used, which is efficient, but
!    generates an even number of values each time.  On any call
!    to this routine, an even number of new values are generated.
!    Depending on the situation, one value may be left over.
!    
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of values desired.
!
!    Output, real ( kind = rk8 ) X(N), a sample of the standard normal PDF.
!
!  Local:
!
!    real ( kind = rk8 ) R(N+1), is used to store some uniform
!    random values.  Its dimension is N+1, but really it is only needed
!    to be the smallest even number greater than or equal to N.
!
!    integer X_LO_INDEX, X_HI_INDEX, records the range of entries of
!    X that we need to compute.  This starts off as 1:N, but is adjusted
!    if we have a saved value that can be immediately stored in X(1),
!    and so on.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  integer m
  real ( kind = rk8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk8 ) r(n+1)
  real ( kind = rk8 ) x(n)
  integer x_hi_index
  integer x_lo_index
!
!  Record the range of X we need to fill in.
!
  x_lo_index = 1
  x_hi_index = n

  if ( x_hi_index - x_lo_index + 1 == 1 ) then

    call random_number ( harvest = r(1:2) )

    x(x_hi_index) = &
             sqrt ( -2.0D+00 * log ( r(1) ) ) * cos ( 2.0D+00 * pi * r(2) )
!
!  If we require an even number of values, that's easy.
!
  else if ( mod ( x_hi_index - x_lo_index + 1, 2 ) == 0 ) then

    m = ( x_hi_index - x_lo_index + 1 ) / 2

    call random_number ( harvest = r(1:2*m) )
 
    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * cos ( 2.0D+00 * pi * r(2:2*m:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * sin ( 2.0D+00 * pi * r(2:2*m:2) )
!
!  If we require an odd number of values, we generate an even number,
!  and handle the last pair specially, storing one in X(N), and
!  saving the other for later.
!
  else

    x_hi_index = x_hi_index - 1

    m = ( x_hi_index - x_lo_index + 1 ) / 2 + 1

    call random_number ( harvest = r(1:2*m) )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * cos ( 2.0D+00 * pi * r(2:2*m-2:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * sin ( 2.0D+00 * pi * r(2:2*m-2:2) )

    x(n) = sqrt ( -2.0D+00 * log ( r(2*m-1) ) ) &
      * cos ( 2.0D+00 * pi * r(2*m) )

  end if

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! r8vec_print() prints an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, real ( kind = rk8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  real ( kind = rk8 ) a(n)
  integer i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
  end do

  return
end
subroutine theta_solve ( a, cl, m, theta )

!*****************************************************************************80
!
!! theta_solve() solves a pair of transcendental equations.
!
!  Discussion:
!
!    The vector THETA returned by this function is needed in order to define
!    the terms in a Karhunen-Loeve expansion of a diffusion coefficient.
!
!    The two equations are:
!
!      1/CL - THETA * TAN ( A * THETA ) = 0
!      THETA - 1/CL * TAN ( A * THETA ) = 0
!
!    A and CL are taken to be positive.  Over each open interval 
!
!      ( n - 1/2 pi, n + 1/2 pi ) / A, for N = 0, 1, ...
!
!    the function TAN ( A * THETA ) monotonically rises from -oo to +00; 
!    therefore, it can be shown that there is one root of each equation 
!    in every interval of this form.  Moreover, because of the positivity
!    of A and CL, we can restrict our search to the interval 
!
!      [ n pi, n + 1/2 pi ) / A, for N = 0, 1, ...
!
!    This function computes K such roots, starting in the first interval,
!    finding those two roots, moving to the next interval, and so on, until
!    the requested number of roots have been found.  Odd index roots will
!    correspond to the first equation, and even index roots to the second.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 July 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Howard Elman, Darran Furnival,
!    Solving the Stochastic Steady-State Diffusion Problem Using Multigrid,
!    University of Maryland Department of Computer Science,
!    Technical Report TR-4786.
!
!  Parameters:
!
!    Input, real ( kind = rk8 ) A, the "radius" of the domain, D = (-A,A)x(-A,A).
!    0 < A.
!
!    Input, real ( kind = rk8 ) CL, the correlation length.
!    0 < CL.
!
!    Input, integer M, the number of values to compute.
!
!    Output, real ( kind = rk8 ) THETA(M), the values of Theta.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer m

  real ( kind = rk8 ) a
  real ( kind = rk8 ) bmatol
  real ( kind = rk8 ) cl
  real ( kind = rk8 ) eps
  real ( kind = rk8 ) fa
  real ( kind = rk8 ) fb
  real ( kind = rk8 ) fc
  real ( kind = rk8 ) ftol
  integer k
  real ( kind = rk8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk8 ) theta(m)
  real ( kind = rk8 ) xa
  real ( kind = rk8 ) xa_init
  real ( kind = rk8 ) xb
  real ( kind = rk8 ) xb_init
  real ( kind = rk8 ) xc

  k = 0
  theta(1:m) = 0.0D+00
!
!  [ XA_INIT, XB_INIT] = [ n * pi, n+1/2 pi ] / a, n = 0, 1, 2, ...
!
  xa_init = 0.0D+00
  xb_init = ( pi / 2.0D+00 ) / a
  eps =  epsilon ( eps )

  do
!
!  Seek root of equation 1 in interval.
!
    if ( m <= k ) then
      exit
    end if

    k = k + 1
    xa = xa_init
    fa = 1.0D+00 / cl - xa * tan ( a * xa )
    ftol = eps * ( abs ( fa ) + 1.0D+00 )
    xb = xb_init
    fb = - fa
    fc = fa
    bmatol = 100.0D+00 * eps * ( abs ( xa ) + abs ( xb ) )

    do while ( bmatol < xb - xa )

      xc = ( xa + xb ) / 2.0D+00
      fc = 1.0D+00 / cl - xc * tan ( a * xc )

      if ( abs ( fc ) <= ftol ) then
        exit
      else if ( 0.0D+00 < fc ) then
        xa = xc
      else
        xb = xc
      end if

    end do

    theta(k) = xc
!
!  Seek root of equation 2 in interval.
!
    if ( m <= k ) then
      exit
    end if

    k = k + 1
!
!  In the first interval, we need to skip the zero root of equation 2.
!
    if ( k == 2 ) then

      k = k - 1

    else

      xa = xa_init
      fa = xa - tan ( a * xa ) / cl
      ftol = eps * ( abs ( fa ) + 1.0D+00 )
      xb = xb_init
      fb = - fa

      do while ( bmatol < xb - xa )

        xc = ( xa + xb ) / 2.0D+00
        fc = xc - tan ( a * xc ) / cl

        if ( abs ( fc ) <= ftol ) then
          exit
        else if ( 0.0D+00 < fc ) then
          xa = xc
        else
          xb = xc
        end if

      end do

      theta(k) = xc

    end if
!
!  Advance the interval.
!
    xa_init = xa_init + pi / a
    xb_init = xb_init + pi / a

  end do

  return
end

