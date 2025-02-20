subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
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
!    This code is distributed under the MIT license.
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

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer ios
  integer iunit
  logical ( kind = 4 ) lopen

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
subroutine grid_2d ( x_num, x_lo, x_hi, y_num, y_lo, y_hi, x, y )

!*****************************************************************************80
!
!! GRID_2D returns a regular 2D grid.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 January 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer X_NUM, the number of X values to use.
!
!    Input, real ( kind = rk ) X_LO, X_HI, the range of X values.
!
!    Input, integer Y_NUM, the number of Y values to use.
!
!    Input, real ( kind = rk ) Y_LO, Y_HI, the range of Y values.
!
!    Output, real ( kind = rk ) X(X_NUM,Y_NUM), Y(X_NUM,Y_NUM), 
!    the coordinates of the grid.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer x_num
  integer y_num

  integer i
  integer j
  real ( kind = rk ) x(x_num,y_num)
  real ( kind = rk ) x_hi
  real ( kind = rk ) x_lo
  real ( kind = rk ) xi
  real ( kind = rk ) y(x_num,y_num)
  real ( kind = rk ) y_hi
  real ( kind = rk ) y_lo
  real ( kind = rk ) yj

  if ( x_num == 1 ) then
    x(1:x_num,1:y_num) = ( x_lo + x_hi ) / 2.0D+00
  else
    do i = 1, x_num
      xi = ( real ( x_num - i,     kind = rk ) * x_lo   &
           + real (         i - 1, kind = rk ) * x_hi ) &
           / real ( x_num     - 1, kind = rk )
      x(i,1:y_num) = xi
    end do
  end if

  if ( y_num == 1 ) then
    y(1:x_num,1:y_num) = ( y_lo + y_hi ) / 2.0D+00
  else
    do j = 1, y_num
      yj = ( real ( y_num - j,     kind = rk ) * y_lo   &
           + real (         j - 1, kind = rk ) * y_hi ) &
           / real ( y_num     - 1, kind = rk )
      y(1:x_num,j) = yj
    end do
  end if

  return
end
subroutine r8_fake_use ( x )

!*****************************************************************************80
!
!! r8_fake_use() pretends to use an R8 variable.
!
!  Discussion:
!
!    Some compilers will issue a warning if a variable is unused.
!    Sometimes there's a good reason to include a variable in a program,
!    but not to use it.  Calling this function with that variable as
!    the argument will shut the compiler up.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 April 2020
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) X, the variable to be "used".
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) x

  if ( x /= x ) then
    write ( *, '(a)' ) '  r8_fake_use: variable is NAN.'
  end if

  return
end
function r8vec_norm_l2 ( n, a )

!*****************************************************************************80
!
!! R8VEC_NORM_L2 returns the L2 norm of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The vector L2 norm is defined as:
!
!      R8VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in A.
!
!    Input, real ( kind = rk ) A(N), the vector whose L2 norm is desired.
!
!    Output, real ( kind = rk ) R8VEC_NORM_L2, the L2 norm of A.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) r8vec_norm_l2

  r8vec_norm_l2 = sqrt ( sum ( a(1:n)**2 ) )

  return
end
subroutine r8vec_uniform_ab ( n, a, b, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_AB returns a scaled pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    Each dimension ranges from A to B.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real ( kind = rk ) A, B, the lower and upper limits.
!
!    Input/output, integer SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = rk ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a
  real ( kind = rk ) b
  integer i
  integer, parameter :: i4_huge = 2147483647
  integer k
  integer seed
  real ( kind = rk ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_AB - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    r(i) = a + ( b - a ) * real ( seed, kind = rk ) * 4.656612875D-10

  end do

  return
end
subroutine resid_stokes1 ( n, x, y, ur, vr, pr )

!*****************************************************************************80
!
!! RESID_STOKES1 returns residuals of the exact Stokes solution #1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 January 2015
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Junping Wang, Yanqiu Wang, Xiu Ye,
!    A robust numerical method for Stokes equations based on divergence-free
!    H(div) finite element methods,
!    SIAM Journal on Scientific Computing,
!    Volume 31, Number 4, 2009, pages 2784-2802.
!
!  Parameters:
!
!    Input, integer N, the number of evaluation points.
!
!    Input, real ( kind = rk ) X(N), Y(N), the coordinates of the points.
!
!    Output, real ( kind = rk ) UR(N), VR(N), PR(N), the residuals in the U,
!    V and P equations.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) f(n)
  real ( kind = rk ) g(n)
  real ( kind = rk ) h(n)
  integer i
  real ( kind = rk ) p
  real ( kind = rk ) pr(n)
  real ( kind = rk ) px
  real ( kind = rk ) py
  real ( kind = rk ) u
  real ( kind = rk ) ur(n)
  real ( kind = rk ) ux
  real ( kind = rk ) uxx
  real ( kind = rk ) uy
  real ( kind = rk ) uyy
  real ( kind = rk ) v
  real ( kind = rk ) vr(n)
  real ( kind = rk ) vx
  real ( kind = rk ) vxx
  real ( kind = rk ) vy
  real ( kind = rk ) vyy
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)
!
!  Get the right hand sides.
!
  call rhs_stokes1 ( n, x, y, f, g, h )
!
!  Form the functions and derivatives.
!
  do i = 1, n

    u = - 2.0D+00 &
          * x(i) ** 2 * ( x(i) - 1.0D+00 ) ** 2  &
          * y(i) * ( y(i) - 1.0D+00 ) * ( 2.0D+00 * y(i) - 1.0D+00 )

    ux = - 2.0D+00  &
          * ( 4.0D+00 * x(i) ** 3 - 6.0D+00 * x(i) ** 2  &
          + 2.0D+00 * x(i) ) &
          * y(i) * ( y(i) - 1.0D+00 ) * ( 2.0D+00 * y(i) - 1.0D+00 )

    uxx = - 2.0D+00  &
          * ( 12.0D+00 * x(i) ** 2 - 12.0D+00 * x(i) + 2.0D+00 ) &
          * ( 2.0D+00 * y(i) ** 3 - 3.0D+00 * y(i) **2 + y(i) )

    uy = - 2.0D+00  &
          * x(i) ** 2 * ( x(i) - 1.0D+00 ) ** 2  &
          * ( 6.0D+00 * y(i) ** 2 - 3.0D+00 * y(i) + 1.0D+00 )

    uyy = - 2.0D+00  &
          * ( x(i) ** 4 - 2.0D+00 * x(i) ** 3 + x(i) ** 2 ) &
          * ( 12.0D+00 * y(i) - 6.0D+00 )

    v =   2.0D+00  &
          * x(i) * ( x(i) - 1.0D+00 ) * ( 2.0D+00 * x(i) - 1.0D+00 ) &
          * y(i) ** 2 * ( y(i) - 1.0D+00 ) ** 2 

    vx =   2.0D+00  &
          * ( 6.0D+00 * x(i) ** 2 - 6.0D+00 * x(i) + 1.0D+00 ) &
          * y(i) ** 2 * ( y(i) - 1.0D+00 ) ** 2 

    vxx =   2.0D+00  &
          * ( 12.0D+00 * x(i) - 6.0D+00 ) &
          * y(i) ** 2 * ( y(i) - 1.0D+00 ) ** 2 

    vy =   2.0D+00  &
          * x(i) * ( x(i) - 1.0D+00 ) * ( 2.0D+00 * x(i) - 1.0D+00 ) &
          * ( 4.0D+00 * y(i) ** 3 - 6.0D+00 * y(i) ** 2  &
          + 2.0D+00 * y(i) )

    vyy =   2.0D+00  &
          * x(i) * ( x(i) - 1.0D+00 ) * ( 2.0D+00 * x(i) - 1.0D+00 ) &
          * ( 12.0D+00 * y(i) ** 2 - 12.0D+00 * y(i) + 2.0D+00 )

    p = 0.0D+00
    px = 0.0D+00
    py = 0.0D+00

    ur(i) = px - ( uxx + uyy ) - f(i)
    vr(i) = py - ( vxx + vyy ) - g(i)
    pr(i) = ux + vy - h(i)

  end do

  return
end
subroutine resid_stokes2 ( n, x, y, ur, vr, pr )

!*****************************************************************************80
!
!! RESID_STOKES2 returns residuals of the exact Stokes solution #2.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 January 2015
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Junping Wang, Yanqiu Wang, Xiu Ye,
!    A robust numerical method for Stokes equations based on divergence-free
!    H(div) finite element methods,
!    SIAM Journal on Scientific Computing,
!    Volume 31, Number 4, 2009, pages 2784-2802.
!
!  Parameters:
!
!    Input, integer N, the number of evaluation points.
!
!    Input, real ( kind = rk ) X(N), Y(N), the coordinates of the points.
!
!    Output, real ( kind = rk ) UR(N), VR(N), PR(N), the residuals in the U,
!    V and P equations.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) f(n)
  real ( kind = rk ) g(n)
  real ( kind = rk ) h(n)
  integer i
  real ( kind = rk ) p
  real ( kind = rk ) pr(n)
  real ( kind = rk ) px
  real ( kind = rk ) py
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) u
  real ( kind = rk ) ur(n)
  real ( kind = rk ) ux
  real ( kind = rk ) uxx
  real ( kind = rk ) uy
  real ( kind = rk ) uyy
  real ( kind = rk ) v
  real ( kind = rk ) vr(n)
  real ( kind = rk ) vx
  real ( kind = rk ) vxx
  real ( kind = rk ) vy
  real ( kind = rk ) vyy
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)
!
!  Get the right hand sides.
!
  call rhs_stokes2 ( n, x, y, f, g, h )

  do i = 1, n

    u =   2.0D+00 &
          * sin ( 2.0D+00 * r8_pi * x(i) ) &
          * cos ( 2.0D+00 * r8_pi * y(i) )

    ux =   4.0D+00 * r8_pi &
          * cos ( 2.0D+00 * r8_pi * x(i) ) &
          * cos ( 2.0D+00 * r8_pi * y(i) )

    uxx = - 8.0D+00 * r8_pi ** 2 &
          * sin ( 2.0D+00 * r8_pi * x(i) ) &
          * cos ( 2.0D+00 * r8_pi * y(i) )

    uy = - 4.0D+00 * r8_pi &
          * sin ( 2.0D+00 * r8_pi * x(i) ) &
          * sin ( 2.0D+00 * r8_pi * y(i) )

    uyy = - 8.0D+00 * r8_pi ** 2 &
          * sin ( 2.0D+00 * r8_pi * x(i) ) &
          * cos ( 2.0D+00 * r8_pi * y(i) )

    v = - 2.0D+00 &
          * cos ( 2.0D+00 * r8_pi * x(i) ) &
          * sin ( 2.0D+00 * r8_pi * y(i) )

    vx =   4.0D+00 * r8_pi &
          * sin ( 2.0D+00 * r8_pi * x(i) ) &
          * sin ( 2.0D+00 * r8_pi * y(i) )

    vxx =   8.0D+00 * r8_pi ** 2 &
          * cos ( 2.0D+00 * r8_pi * x(i) ) &
          * sin ( 2.0D+00 * r8_pi * y(i) )

    vy = - 4.0D+00 * r8_pi &
          * cos ( 2.0D+00 * r8_pi * x(i) ) &
          * cos ( 2.0D+00 * r8_pi * y(i) )

    vyy =   8.0D+00 * r8_pi ** 2 &
          * cos ( 2.0D+00 * r8_pi * x(i) ) &
          * sin ( 2.0D+00 * r8_pi * y(i) )

    p = x(i) ** 2 + y(i) ** 2

    px = 2.0D+00 * x(i)
    py = 2.0D+00 * y(i)

    ur(i) = px - ( uxx + uyy ) - f(i)
    vr(i) = py - ( vxx + vyy ) - g(i)
    pr(i) = ux + uy - h(i)

  end do

  return
end
subroutine resid_stokes3 ( n, x, y, ur, vr, pr )

!*****************************************************************************80
!
!! RESID_STOKES3 returns residuals of the exact Stokes solution #3.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 February 2015
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Howard Elman, Alison Ramage, David Silvester,
!    Finite Elements and Fast Iterative Solvers with
!    Applications in Incompressible Fluid Dynamics,
!    Oxford, 2005,
!    ISBN: 978-0198528678,
!    LC: QA911.E39.
!
!  Parameters:
!
!    Input, integer N, the number of evaluation points.
!
!    Input, real ( kind = rk ) X(N), Y(N), the coordinates of the points.
!
!    Output, real ( kind = rk ) UR(N), VR(N), PR(N), the residuals in the U,
!    V and P equations.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) f(n)
  real ( kind = rk ) g(n)
  real ( kind = rk ) h(n)
  integer i
  real ( kind = rk ) p
  real ( kind = rk ) pr(n)
  real ( kind = rk ) px
  real ( kind = rk ) py
  real ( kind = rk ) u
  real ( kind = rk ) ur(n)
  real ( kind = rk ) ux
  real ( kind = rk ) uxx
  real ( kind = rk ) uy
  real ( kind = rk ) uyy
  real ( kind = rk ) v
  real ( kind = rk ) vr(n)
  real ( kind = rk ) vx
  real ( kind = rk ) vxx
  real ( kind = rk ) vy
  real ( kind = rk ) vyy
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)
!
!  Get the right hand sides.
!
  call rhs_stokes3 ( n, x, y, f, g, h )
!
!  Form the functions and derivatives.
!
  do i = 1, n

    u =   20.0D+00 * x(i) * y(i) ** 3
    ux = 20.0D+00 * y(i) ** 3
    uxx = 0.0D+00
    uy = 60.0D+00 * x(i) * y(i) ** 2
    uyy = 120.0D+00 * x(i) * y(i)

    v = 5.0D+00 * ( x(i) ** 4  - y(i) ** 4 )
    vx = 20.0D+00 * x(i) ** 3
    vxx = 60.0D+00 * x(i) ** 2
    vy = - 20.0D+00 * y(i) ** 3
    vyy = - 60.0D+00 * y(i) ** 2

    p =   60.0D+00 * x(i) ** 2 * y(i) - 20.0D+00 * y(i) ** 3 + 10.0D+00
    px = 120.0D+00 * x(i) * y(i)
    py =  60.0D+00 * x(i) ** 2 - 60.0D+00 * y(i) ** 2

    ur(i) = px - ( uxx + uyy ) - f(i)
    vr(i) = py - ( vxx + vyy ) - g(i)
    pr(i) = ux + vy - h(i)

  end do

  return
end
subroutine rhs_stokes1 ( n, x, y, f, g, h )

!*****************************************************************************80
!
!! RHS_STOKES1 returns the right hand sides of the exact Stokes solution #1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 January 2015
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Junping Wang, Yanqiu Wang, Xiu Ye,
!    A robust numerical method for Stokes equations based on divergence-free
!    H(div) finite element methods,
!    SIAM Journal on Scientific Computing,
!    Volume 31, Number 4, 2009, pages 2784-2802.
!
!  Parameters:
!
!    Input, integer N, the number of evaluation points.
!
!    Input, real ( kind = rk ) X(N), Y(N), the coordinates of the points.
!
!    Output, real ( kind = rk ) F(N), G(N), H(N), the right hand sides in the U,
!    V and P equations.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) f(n)
  real ( kind = rk ) g(n)
  real ( kind = rk ) h(n)
  integer i
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)

  do i = 1, n

    f(i) = + 2.0D+00 &
          * ( 12.0D+00 * x(i) ** 2 - 12.0D+00 * x(i) + 2.0D+00 ) &
          * ( 2.0D+00 * y(i) ** 3 - 3.0D+00 * y(i) **2 + y(i) ) &
          + 2.0D+00 &
          * ( x(i) ** 4 - 2.0D+00 * x(i) ** 3 + x(i) ** 2 ) &
          * ( 12.0D+00 * y(i) - 6.0D+00 )

    g(i) = - 2.0D+00 &
          * ( 12.0D+00 * x(i) - 6.0D+00 ) &
          * ( y(i) ** 4 - 2.0D+00 * y(i) ** 3 + y(i) ** 2 ) &
          - 2.0D+00 &
          * ( 2.0D+00 * x(i) ** 3 - 3.0D+00 * x(i) ** 2 + x(i) ) &
          * ( 12.0D+00 * y(i) ** 2 - 12.0D+00 * y(i) + 2.0D+00 )

    h(i) = 0.0D+00

  end do

  return
end
subroutine rhs_stokes2 ( n, x, y, f, g, h )

!*****************************************************************************80
!
!! RHS_STOKES2 returns the right hand sides of the exact Stokes solution #2.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 January 2015
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Junping Wang, Yanqiu Wang, Xiu Ye,
!    A robust numerical method for Stokes equations based on divergence-free
!    H(div) finite element methods,
!    SIAM Journal on Scientific Computing,
!    Volume 31, Number 4, 2009, pages 2784-2802.
!
!  Parameters:
!
!    Input, integer N, the number of evaluation points.
!
!    Input, real ( kind = rk ) X(N), Y(N), the coordinates of the points.
!
!    Output, real ( kind = rk ) F(N), G(N), H(N), the right hand sides in the U,
!    V and P equations.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) f(n)
  real ( kind = rk ) g(n)
  real ( kind = rk ) h(n)
  integer i
  real ( kind = rk ) p
  real ( kind = rk ) px
  real ( kind = rk ) py
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) u
  real ( kind = rk ) ux
  real ( kind = rk ) uxx
  real ( kind = rk ) uy
  real ( kind = rk ) uyy
  real ( kind = rk ) v
  real ( kind = rk ) vx
  real ( kind = rk ) vxx
  real ( kind = rk ) vy
  real ( kind = rk ) vyy
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)

  do i = 1, n

    u =   2.0D+00 &
          * sin ( 2.0D+00 * r8_pi * x(i) ) &
          * cos ( 2.0D+00 * r8_pi * y(i) )

    ux =   4.0D+00 * r8_pi &
          * cos ( 2.0D+00 * r8_pi * x(i) ) &
          * cos ( 2.0D+00 * r8_pi * y(i) )

    uxx = - 8.0D+00 * r8_pi ** 2 &
          * sin ( 2.0D+00 * r8_pi * x(i) ) &
          * cos ( 2.0D+00 * r8_pi * y(i) )

    uy = - 4.0D+00 * r8_pi &
          * sin ( 2.0D+00 * r8_pi * x(i) ) &
          * sin ( 2.0D+00 * r8_pi * y(i) )

    uyy = - 8.0D+00 * r8_pi ** 2 &
          * sin ( 2.0D+00 * r8_pi * x(i) ) &
          * cos ( 2.0D+00 * r8_pi * y(i) )

    v = - 2.0D+00 &
          * cos ( 2.0D+00 * r8_pi * x(i) ) &
          * sin ( 2.0D+00 * r8_pi * y(i) )

    vx =   4.0D+00 * r8_pi &
          * sin ( 2.0D+00 * r8_pi * x(i) ) &
          * sin ( 2.0D+00 * r8_pi * y(i) )

    vxx =   8.0D+00 * r8_pi ** 2 &
          * cos ( 2.0D+00 * r8_pi * x(i) ) &
          * sin ( 2.0D+00 * r8_pi * y(i) )

    vy = - 4.0D+00 * r8_pi &
          * cos ( 2.0D+00 * r8_pi * x(i) ) &
          * cos ( 2.0D+00 * r8_pi * y(i) )

    vyy =   8.0D+00 * r8_pi ** 2 &
          * cos ( 2.0D+00 * r8_pi * x(i) ) &
          * sin ( 2.0D+00 * r8_pi * y(i) )

    p = x(i) ** 2 + y(i) ** 2

    px = 2.0D+00 * x(i)
    py = 2.0D+00 * y(i)

    f(i) = px - ( uxx + uyy )
    g(i) = py - ( vxx + vyy )
    h(i) = ux + vy

  end do

  return
end
subroutine rhs_stokes3 ( n, x, y, f, g, h )

!*****************************************************************************80
!
!! RHS_STOKES3 returns the right hand sides of the exact Stokes solution #3.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 February 2015
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Howard Elman, Alison Ramage, David Silvester,
!    Finite Elements and Fast Iterative Solvers with
!    Applications in Incompressible Fluid Dynamics,
!    Oxford, 2005,
!    ISBN: 978-0198528678,
!    LC: QA911.E39.
!
!  Parameters:
!
!    Input, integer N, the number of evaluation points.
!
!    Input, real ( kind = rk ) X(N), Y(N), the coordinates of the points.
!
!    Output, real ( kind = rk ) F(N), G(N), H(N), the right hand sides in the U,
!    V and P equations.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) f(n)
  real ( kind = rk ) g(n)
  real ( kind = rk ) h(n)
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)

  call r8_fake_use ( x(1) )
  call r8_fake_use ( y(1) )

  f(1:n) = 0.0D+00
  g(1:n) = 0.0D+00
  h(1:n) = 0.0D+00

  return
end
subroutine stokes_gnuplot ( header, n, x, y, u, v, s )

!*****************************************************************************80
!
!! STOKES_GNUPLOT writes the Stokes vector field to files for GNUPLOT.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 January 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) HEADER, a header to be used to name the files.
!
!    Input, integer N, the number of evaluation points.
!
!    Input, real ( kind = rk ) X(N), Y(N), the coordinates of the 
!    evaluation points.
!
!    Input, real ( kind = rk ) U(N), V(N), the velocity components.
!
!    Input, real ( kind = rk ) S, a scale factor for the velocity vectors.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  character ( len = 255 ) command_filename
  integer command_unit
  character ( len = 255 ) data_filename
  integer data_unit
  character ( len = * ) header
  integer i
  character ( len = 255 ) plot_filename
  real ( kind = rk ) s
  real ( kind = rk ) u(n)
  real ( kind = rk ) v(n)
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)
!
!  Write the data file.
!
  data_filename = trim ( header ) // '_data.txt'

  call get_unit ( data_unit )

  open ( unit = data_unit, file = data_filename, status = 'replace' )

  do i = 1, n
    write ( data_unit, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      x(i), y(i), s * u(i), s * v(i)
  end do

  close ( unit = data_unit )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Data written to "' // trim ( data_filename ) // '".'
!
!  Write the command file.
!
  command_filename = trim ( header ) // '_commands.txt'
  call get_unit ( command_unit )

  plot_filename = trim ( header ) // '.png'

  open ( unit = command_unit, file = command_filename, status = 'replace' )

  write ( command_unit, '(a)' ) '#  ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) 'set output "' // trim ( plot_filename ) // '"'
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '#  Add titles and labels.'
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set xlabel "<--- X --->"'
  write ( command_unit, '(a)' ) 'set ylabel "<--- Y --->"'
  write ( command_unit, '(a)' ) 'set title "Stokes flow"'
  write ( command_unit, '(a)' ) 'unset key'
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '#  Add grid lines.'
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set size ratio -1'
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '#  Timestamp the plot.'
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set timestamp'
  write ( command_unit, '(a)' ) 'plot "' // trim ( data_filename ) // &
    '" using 1:2:3:4 with vectors \'
  write ( command_unit, '(a)' ) '  head filled lt 2 linecolor rgb "blue"'
  write ( command_unit, '(a)' ) 'quit'

  close ( unit = command_unit )

  write ( *, '(a)' ) '  Commands written to "' // &
    trim ( command_filename ) // '".'

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
!  Parameters:
!
!    None
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

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
subroutine uvp_stokes1 (  n, x, y, u, v, p )

!*****************************************************************************80
!
!! UVP_STOKES1 evaluates the exact Stokes solution #1.
!
!  Discussion:
!
!    The solution is defined over the unit square [0,1]x[0,1].
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 January 2015
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Junping Wang, Yanqiu Wang, Xiu Ye,
!    A robust numerical method for Stokes equations based on divergence-free
!    H(div) finite element methods,
!    SIAM Journal on Scientific Computing,
!    Volume 31, Number 4, 2009, pages 2784-2802.
!
!  Parameters:
!
!    Input, integer N, the number of evaluation points.
!
!    Input, real ( kind = rk ) X(N), Y(N), the coordinates of the points.
!
!    Output, real ( kind = rk ) U(N), V(N), P(N), the velocity components and
!    pressure at each of the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  real ( kind = rk ) p(n)
  real ( kind = rk ) u(n)
  real ( kind = rk ) v(n)
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)

  do i = 1, n

    u(i) = - 2.0D+00 &
          * x(i) ** 2 * ( x(i) - 1.0D+00 ) ** 2 &
          * y(i) * ( y(i) - 1.0D+00 ) * ( 2.0D+00 * y(i) - 1.0D+00 )

    v(i) =   2.0D+00 &
          * x(i) * ( x(i) - 1.0D+00 ) * ( 2.0D+00 * x(i) - 1.0D+00 ) &
          * y(i) ** 2 * ( y(i) - 1.0D+00 ) ** 2 

    p(i) = 0.0D+00

  end do

  return
end
subroutine uvp_stokes2 (  n, x, y, u, v, p )

!*****************************************************************************80
!
!! UVP_STOKES2 evaluates the exact Stokes solution #2.
!
!  Discussion:
!
!    The solution is defined over the unit square [0,1]x[0,1].
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 January 2015
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Junping Wang, Yanqiu Wang, Xiu Ye,
!    A robust numerical method for Stokes equations based on divergence-free
!    H(div) finite element methods,
!    SIAM Journal on Scientific Computing,
!    Volume 31, Number 4, 2009, pages 2784-2802.
!
!  Parameters:
!
!    Input, integer N, the number of evaluation points.
!
!    Input, real ( kind = rk ) X(N), Y(N), the coordinates of the points.
!
!    Output, real ( kind = rk ) U(N), V(N), P(N), the velocity components and
!    pressure at each of the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  real ( kind = rk ) p(n)
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) u(n)
  real ( kind = rk ) v(n)
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)

  do i = 1, n

    u(i) =   2.0D+00 &
          * sin ( 2.0D+00 * r8_pi * x(i) ) &
          * cos ( 2.0D+00 * r8_pi * y(i) )


    v(i) = - 2.0D+00 &
          * cos ( 2.0D+00 * r8_pi * x(i) ) &
          * sin ( 2.0D+00 * r8_pi * y(i) )

    p(i) = x(i) ** 2 + y(i) ** 2

  end do

  return
end
subroutine uvp_stokes3 (  n, x, y, u, v, p )

!*****************************************************************************80
!
!! UVP_STOKES3 evaluates the exact Stokes solution #3.
!
!  Discussion:
!
!    The solution is defined over the unit square [-1,+1]x[-1,+1].
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 February 2015
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Howard Elman, Alison Ramage, David Silvester,
!    Finite Elements and Fast Iterative Solvers with
!    Applications in Incompressible Fluid Dynamics,
!    Oxford, 2005,
!    ISBN: 978-0198528678,
!    LC: QA911.E39.
!
!  Parameters:
!
!    Input, integer N, the number of evaluation points.
!
!    Input, real ( kind = rk ) X(N), Y(N), the coordinates of the points.
!
!    Output, real ( kind = rk ) U(N), V(N), P(N), the velocity components and
!    pressure at each of the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  real ( kind = rk ) p(n)
  real ( kind = rk ) u(n)
  real ( kind = rk ) v(n)
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)

  do i = 1, n

    u(i) =   20.0D+00 * x(i) * y(i) ** 3
    v(i) =    5.0D+00 * ( x(i) ** 4  - y(i) ** 4 )
    p(i) =   60.0D+00 * x(i) ** 2 * y(i) - 2.0D+00 * y(i) ** 3 + 10.0D+00

  end do

  return
end
