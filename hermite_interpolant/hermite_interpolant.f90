subroutine dif_deriv ( nd, xd, yd, ndp, xdp, ydp )

!*****************************************************************************80
!
!! dif_deriv() computes the derivative of a polynomial in divided difference form.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 June 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, integer ND, the size of the input table.
!
!    Input, real ( kind = rk ) XD(ND), the abscissas for the divided
!    difference table.
!
!    Input, real ( kind = rk ) YD(ND), the divided difference table.
!
!    Output, integer NDP, the size of the output table, 
!    which is ND-1.
!
!    Input, real ( kind = rk ) XDP(NDP), the abscissas for the divided
!    difference table for the derivative.
!
!    Output, real ( kind = rk ) YDP(NDP), the divided difference 
!    table for the derivative.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nd

  integer i
  integer ndp
  real ( kind = rk ) xd(nd)
  real ( kind = rk ) xd_temp(nd)
  real ( kind = rk ) xdp(nd-1)
  real ( kind = rk ) yd(nd)
  real ( kind = rk ) yd_temp(nd)
  real ( kind = rk ) ydp(nd-1)
!
!  Using a temporary copy of the difference table, shift the abscissas to zero.
!
  xd_temp(1:nd) = xd(1:nd)
  yd_temp(1:nd) = yd(1:nd)

  call dif_shift_zero ( nd, xd_temp, yd_temp )
!
!  Now construct the derivative.
!
  ndp = nd - 1

  xdp(1:ndp) = 0.0D+00

  do i = 1, ndp
    ydp(i) = real ( i, kind = rk ) * yd_temp(i+1)
  end do

  return
end
subroutine dif_shift_x ( nd, xd, yd, xv )

!*****************************************************************************80
!
!! DIF_SHIFT_X replaces one abscissa of a divided difference table.
!
!  Discussion:
!
!    The routine shifts the representation of a divided difference polynomial by
!    dropping the last X value in XD, and adding a new value XV to the
!    beginning of the XD array, suitably modifying the coefficients stored
!    in YD.
!
!    The representation of the polynomial is changed, but the polynomial itself
!    should be identical.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 June 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, integer ND, the number of divided difference 
!    coefficients, and the number of entries in XD.
!
!    Input/output, real ( kind = rk ) XD(ND), the X values used in 
!    the representation of the divided difference polynomial.
!    After a call to this routine, the last entry of XD has been dropped,
!    the other entries have shifted up one index, and XV has been inserted
!    at the beginning of the array.
!
!    Input/output, real ( kind = rk ) YD(ND), the divided difference
!    coefficients corresponding to the XD array.  On output, this 
!    array has been adjusted.
!
!    Input, real ( kind = rk ) XV, a new X value which is to be used in 
!    the representation of the polynomial.  On output, XD(1) equals 
!    XV and the representation of the polynomial has been suitably changed.
!    Note that XV does not have to be distinct from any of the original XD
!    values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nd

  integer i
  real ( kind = rk ) xd(nd)
  real ( kind = rk ) xv
  real ( kind = rk ) yd(nd)
!
!  Recompute the divided difference coefficients.
!
  do i = nd - 1, 1, -1
    yd(i) = yd(i) + ( xv - xd(i) ) * yd(i+1)
  end do
!
!  Shift the XD values up one position, and insert XV at the beginning.
!
  xd(2:nd) = xd(1:nd-1)

  xd(1) = xv

  return
end
subroutine dif_shift_zero ( nd, xd, yd )

!*****************************************************************************80
!
!! DIF_SHIFT_ZERO shifts a divided difference table so all abscissas are zero.
!
!  Discussion:
!
!    When the abscissas are changed, the coefficients naturally
!    must also be changed.
!
!    The resulting pair (XD, YD) still represents the
!    same polynomial, but the entries in YD are now the
!    standard polynomial coefficients.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 June 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, integer ND, the length of the XD and YD arrays.
!
!    Input/output, real ( kind = rk ) XD(ND), the X values that 
!    correspond to the divided difference table.  On output, XD 
!    contains only zeroes.
!
!    Input/output, real ( kind = rk ) YD(ND), the divided difference table
!    for the polynomial.  On output, YD is also the coefficient array for 
!    the standard representation of the polynomial.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nd

  integer i
  real ( kind = rk ) xd(nd)
  real ( kind = rk ) xv
  real ( kind = rk ) yd(nd)

  xv = 0.0D+00

  do i = 1, nd
    call dif_shift_x ( nd, xd, yd, xv )
  end do

  return
end
subroutine dif_to_r8poly ( nd, xd, yd, c )

!*****************************************************************************80
!
!! DIF_TO_R8POLY converts a divided difference table to a standard polynomial.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, integer ND, the number of coefficients, and abscissas.
!
!    Input, real ( kind = rk ) XD(ND), the X values used in the divided
!    difference representation of the polynomial.
!
!    Input, real ( kind = rk ) YD(ND) the divided difference table.
!
!    Output, real ( kind = rk ) C(ND), the polyomial coefficients.
!    C(1) is the constant term, and C(ND) is the coefficient of X^(ND-1).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nd

  real ( kind = rk ) c(nd)
  integer i
  integer j
  real ( kind = rk ) xd(nd)
  real ( kind = rk ) yd(nd)

  c(1:nd) = yd(1:nd)
!
!  Recompute the divided difference coefficients.
!
  do j = 1, nd - 1
    do i = 1, nd - j
      c(nd-i) = c(nd-i) - xd(nd-i-j+1) * c(nd-i+1)
    end do
  end do

  return
end
subroutine dif_vals ( nd, xd, yd, nv, xv, yv )

!*****************************************************************************80
!
!! DIF_VALS evaluates a divided difference polynomial at a set of points.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, integer ND, the order of the difference table.
!
!    Input, real ( kind = rk ) XD(ND), the X values of the difference table.
!
!    Input, real ( kind = rk ) YD(ND), the divided differences.
!
!    Input, integer NV, the number of evaluation points.
!
!    Input, real ( kind = rk ) XV(NV), the evaluation points.
!
!    Output, real ( kind = rk ) YV(NV), the value of the divided difference
!    polynomial at the evaluation points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nd
  integer nv

  integer i
  real ( kind = rk ) xd(nd)
  real ( kind = rk ) xv(nv)
  real ( kind = rk ) yd(nd)
  real ( kind = rk ) yv(nv)

  yv(1:nv) = yd(nd)
  do i = 1, nd - 1
    yv(1:nv) = yd(nd-i) + ( xv(1:nv) - xd(nd-i) ) * yv(1:nv)
  end do

  return
end
subroutine hermite_basis_0 ( n, x, i, xv, value )

!*****************************************************************************80
!
!! HERMITE_BASIS_0 evaluates a zero-order Hermite interpolation basis function.
!
!  Discussion:
!
!    Given ND points XD, with values YD and derivative values YPD, the
!    Hermite interpolant can be written as:
!
!      H(X) = sum ( 1 <= I <= ND ) YD(I)  * H0(I;X)
!           + sum ( 1 <= I <= ND ) YPD(I) * H1(I;X)
!
!    where H0(I;X) is the I-th zero order Hermite interpolation basis function,
!    and H1(I;X) is the I-th first order Hermite interpolation basis function.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of abscissas.
!
!    Input, real ( kind = rk ) X(N), the abscissas.
!
!    Input, integer I, the index of the first-order basis function.
!
!    Input, real ( kind = rk ) XV, the evaluation point.
!
!    Output, real ( kind = rk ) VALUE, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) factor(n)
  integer i
  integer j
  real ( kind = rk ) li
  real ( kind = rk ) lp
  real ( kind = rk ) lpp
  real ( kind = rk ) value
  real ( kind = rk ) x(n)
  real ( kind = rk ) xv

  if ( i < 1 .or. n < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HERMITE_BASIS_0 - Fatal error!'
    write ( *, '(a)' ) '  I < 1 or N < I.'
    stop 1
  end if
!
!  L(X) = product ( X - X(1:N) )
!
!  L'(X(I)).
!
  factor(1:n) = x(i) - x(1:n)
  factor(i) = 1.0D+00

  lp = product ( factor(1:n) )
!
!  LI(X) = L(X) / ( X - X(I) ) / L'(X(I))
!
  factor(1:n) = xv - x(1:n)
  factor(i) = 1.0D+00

  li = product ( factor(1:n) ) / lp
!
!  L''(X(I)).
!
  lpp = 0.0D+00
  factor(1:n) = x(i) - x(1:n)
  factor(i) = 1.0D+00

  do j = 1, n
    if ( j /= i ) then
      factor(j) = 1.0D+00
      lpp = lpp + 2.0D+00 * product ( factor )
      factor(j) = x(i) - x(j)
    end if
  end do

  value = ( 1.0D+00 - ( xv - x(i) ) * lpp / lp ) * li * li

  return
end
subroutine hermite_basis_1 ( n, x, i, xv, value )

!*****************************************************************************80
!
!! HERMITE_BASIS_1 evaluates a first-order Hermite interpolation basis function.
!
!  Discussion:
!
!    Given ND points XD, with values YD and derivative values YPD, the
!    Hermite interpolant can be written as:
!
!      H(X) = sum ( 1 <= I <= ND ) YD(I)  * H0(I;X)
!           + sum ( 1 <= I <= ND ) YPD(I) * H1(I;X)
!
!    where H0(I;X) is the I-th zero order Hermite interpolation basis function,
!    and H1(I;X) is the I-th first order Hermite interpolation basis function.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of abscissas.
!
!    Input, real ( kind = rk ) X(N), the abscissas.
!
!    Input, integer I, the index of the first-order basis function.
!
!    Input, real ( kind = rk ) XV, the evaluation point.
!
!    Output, real ( kind = rk ) VALUE, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) bot
  real ( kind = rk ) factor(n)
  integer i
  real ( kind = rk ) top
  real ( kind = rk ) value
  real ( kind = rk ) x(n)
  real ( kind = rk ) xv

  if ( i < 1 .or. n < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HERMITE_BASIS_1 - Fatal error!'
    write ( *, '(a)' ) '  I < 1 or N < I.'
    stop 1
  end if

  factor(1:n) = xv - x(1:n)
  factor(i) = 1.0D+00
  top = product ( factor(1:n) )

  factor(1:n) = x(i) - x(1:n)
  factor(i) = 1.0D+00
  bot = product ( factor(1:n) )

  value = ( xv - x(i) ) * ( top / bot ) * ( top / bot )

  return
end
subroutine hermite_demo ( n, x, y, yp )

!*****************************************************************************80
!
!! HERMITE_DEMO computes and prints Hermite interpolant information for data.
!
!  Discussion:
!
!    Given a set of Hermite data, this routine calls HDATA_TO_DIF to determine
!    and print the divided difference table, and then DIF_TO_R8POLY to 
!    determine and print the coefficients of the polynomial in standard form.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data points.
!
!    Input, real ( kind = rk ) X(N), the abscissas.
!
!    Input, real ( kind = rk ) Y(N), YP(N), the function and derivative
!    values at the abscissas.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) cd(0:2*n-1)
  integer i
  integer nd
  integer ndp
  integer nv
  real ( kind = rk ) x(n)
  real ( kind = rk ) xd(2*n)
  real ( kind = rk ) xdp(2*n-1)
  real ( kind = rk ) xv(n)
  real ( kind = rk ) y(n)
  real ( kind = rk ) yd(2*n)
  real ( kind = rk ) ydp(2*n-1)
  real ( kind = rk ) yp(n)
  real ( kind = rk ) yv(n)
  real ( kind = rk ) yvp(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_DEMO'
  write ( *, '(a)' ) '  Compute coefficients CD of the Hermite polynomial'
  write ( *, '(a)' ) '  interpolant to given data (x,y,yp).'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data:'
  write ( *, '(a)' ) '              X           Y           Y'''
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4)' ) i, x(i), y(i), yp(i)
  end do

  call hermite_interpolant ( n, x, y, yp, xd, yd, xdp, ydp )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Difference table for interpolant:'
  write ( *, '(a)' ) '              XD          YD'
  write ( *, '(a)' ) ' '

  nd = 2 * n

  do i = 1, nd
    write ( *, '(2x,i4,2x,f10.4,2x,f10.4)' ) i, xd(i), yd(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Difference table for derivative:'
  write ( *, '(a)' ) '              XDP         YDP'
  write ( *, '(a)' ) ' '

  ndp = 2 * n - 1

  do i = 1, ndp
    write ( *, '(2x,i4,2x,f10.4,2x,f10.4)' ) i, xdp(i), ydp(i)
  end do

  call dif_to_r8poly ( nd, xd, yd, cd )

  call r8poly_print ( nd - 1, cd, '  Hermite interpolating polynomial:' )
!
!  Verify interpolation claim!
!
  nv = n
  xv(1:nv) = x(1:n)

  call hermite_interpolant_value ( nd, xd, yd, xdp, ydp, nv, xv, yv, yvp )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data Versus Interpolant:'
  write ( *, '(a)' ) &
    '              X           Y           H           YP          HP'
  write ( *, '(a)' ) ' '
  do i = 1, nv
    write ( *, '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4,2x,f10.4,2x,f10.4)' ) &
      i, xv(i), y(i), yv(i), yp(i), yvp(i)
  end do

  return
end
subroutine hermite_interpolant ( n, x, y, yp, xd, yd, xdp, ydp )

!*****************************************************************************80
!
!! HERMITE_INTERPOLANT sets up a divided difference table from Hermite data.
!
!  Discussion:
!
!    The polynomial represented by the divided difference table can be
!    evaluated by calling HERMITE_INTERPOLANT_VALUE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, integer N, the number of items of data 
!    ( X(I), Y(I), YP(I) ).
!
!    Input, real ( kind = rk ) X(N), the abscissas.
!    These values must be distinct.
!
!    Input, real ( kind = rk ) Y(N), YP(N), the function and 
!    derivative values.
!
!    Output, real ( kind = rk ) XD(2*N), YD(2*N), the divided difference table
!    for the interpolant value.
!
!    Output, real ( kind = rk ) XDP(2*N-1), YDP(2*N-1), the divided difference 
!    table for the interpolant derivative.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer j
  integer nd
  integer ndp
  real ( kind = rk ) x(n)
  real ( kind = rk ) xd(2*n)
  real ( kind = rk ) xdp(2*n-1)
  real ( kind = rk ) y(n)
  real ( kind = rk ) yd(2*n)
  real ( kind = rk ) ydp(2*n-1)
  real ( kind = rk ) yp(n)
!
!  Copy the data:
!
  nd = 2 * n
  xd(1:nd-1:2) = x(1:n)
  xd(2:nd  :2) = x(1:n)
!
!  Carry out the first step of differencing.
!
  yd(1) = y(1)
  yd(3:nd-1:2) = ( y(2:n) - y(1:n-1) ) / ( x(2:n) - x(1:n-1) )
  yd(2:nd  :2) = yp(1:n)  
!
!  Carry out the remaining steps in the usual way.
!
  do i = 3, nd
    do j = nd, i, -1

      yd(j) = ( yd(j) - yd(j-1) ) / ( xd(j) - xd(j+1-i) )

    end do

  end do
!
!  Compute the difference table for the derivative.
!
  call dif_deriv ( nd, xd, yd, ndp, xdp, ydp )

  return
end
subroutine hermite_interpolant_rule ( n, a, b, x, w )

!*****************************************************************************80
!
!! HERMITE_INTERPOLANT_RULE: quadrature rule for a Hermite interpolant.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 June 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of abscissas.
!
!    Input, real ( kind = rk ) A, B, the integration limits.
!
!    Input, real ( kind = rk ) X(N), the abscissas.
!
!    Output, real ( kind = rk ) W(2*N), the quadrature coefficients, given as
!    pairs for function and derivative values at each abscissa.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a
  real ( kind = rk ) a_value
  real ( kind = rk ) b
  real ( kind = rk ) b_value
  real ( kind = rk ) c(0:2*n-1)
  integer i
  integer k
  integer nd
  real ( kind = rk ) w(2*n)
  real ( kind = rk ) x(n)
  real ( kind = rk ) xd(2*n)
  real ( kind = rk ) xdp(2*n-1)
  real ( kind = rk ) y(n)
  real ( kind = rk ) yd(2*n)
  real ( kind = rk ) ydp(2*n-1)
  real ( kind = rk ) yp(n)

  nd = 2 * n

  k = 0

  do i = 1, n

    k = k + 1
    y(1:n) = 0.0D+00
    y(i) = 1.0D+00
    yp(1:n) = 0.0D+00
    call hermite_interpolant ( n, x, y, yp, xd, yd, xdp, ydp )
    call dif_to_r8poly ( nd, xd, yd, c )
    call r8poly_ant_val ( n, c, a, a_value )
    call r8poly_ant_val ( n, c, b, b_value )
    w(k) = b_value - a_value

    k = k + 1
    y(1:n) = 0.0D+00
    yp(1:n) = 0.0D+00
    yp(i) = 1.0D+00
    call hermite_interpolant ( n, x, y, yp, xd, yd, xdp, ydp )
    call dif_to_r8poly ( nd, xd, yd, c )
    call r8poly_ant_val ( n, c, a, a_value )
    call r8poly_ant_val ( n, c, b, b_value )
    w(k) = b_value - a_value

  end do

  return
end
subroutine hermite_interpolant_value ( nd, xd, yd, xdp, ydp, nv, xv, yv, yvp )

!*****************************************************************************80
!
!! HERMITE_INTERPOLANT_VALUE evaluates the Hermite interpolant polynomial.
!
!  Discussion:
!
!    In fact, this function will evaluate an arbitrary polynomial that is
!    represented by a difference table.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, integer ND, the order of the difference table.
!
!    Input, real ( kind = rk ) XD(ND), YD(ND), the difference table for the
!    interpolant value.
!
!    Input, real ( kind = rk ) XDP(ND-1), YDP(ND-1), the difference table for
!    the interpolant derivative.
!
!    Input, integer NV, the number of evaluation points.
!
!    Input, real ( kind = rk ) XV(NV), the evaluation points.
!
!    Output, real ( kind = rk ) YV(NV), YVP(NV), the value of the interpolant and
!    its derivative at the evaluation points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nd
  integer nv

  integer i
  integer ndp
  real ( kind = rk ) xd(nd)
  real ( kind = rk ) xdp(nd-1)
  real ( kind = rk ) xv(nv)
  real ( kind = rk ) yd(nd)
  real ( kind = rk ) ydp(nd-1)
  real ( kind = rk ) yv(nv)
  real ( kind = rk ) yvp(nv)

  yv(1:nv) = yd(nd)
  do i = nd - 1, 1, -1
    yv(1:nv) = yd(i) + ( xv(1:nv) - xd(i) ) * yv(1:nv)
  end do

  ndp = nd - 1

  yvp(1:nv) = ydp(ndp)
  do i = ndp - 1, 1, -1
    yvp(1:nv) = ydp(i) + ( xv(1:nv) - xdp(i) ) * yvp(1:nv)
  end do

  return
end
subroutine r8poly_ant_val ( n, c, xv, yv )

!*****************************************************************************80
!
!! R8POLY_ANT_VAL evaluates the antiderivative of a polynomial in standard form.
!
!  Discussion:
!
!    The constant term of the antiderivative is taken to be zero.
!
!    Thus, if 
!      p(x) = c(1)     + c(2) * x   + c(3) * x^2   + ... + c(n) * x^(n-1)
!    then q(x), the antiderivative, is taken to be:
!      q(x) = c(1) * x + c(2) * x/2 + c(3) * x^3/3 + ... + c(n) * x^n/n
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the polynomial.
!
!    Input, real ( kind = rk ) C(N), the polynomial coefficients.
!    C(1) is the constant term, and C(N) is the coefficient of X^(N-1).
!
!    Input, real ( kind = rk ) XV, the evaluation point.
!
!    Output, real ( kind = rk ) YV, the value of the antiderivative of 
!    the polynomial at XV.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) c(n)
  integer i
  real ( kind = rk ) xv
  real ( kind = rk ) yv

  yv = 0.0D+00

  do i = n, 1, -1
    yv = ( yv + c(i) / real ( i, kind = rk ) ) * xv
  end do

  return
end
subroutine r8poly_degree ( na, a, degree )

!*****************************************************************************80
!
!! R8POLY_DEGREE returns the degree of a polynomial.
!
!  Discussion:
!
!    The degree of a polynomial is the index of the highest power
!    of X with a nonzero coefficient.
!
!    The degree of a constant polynomial is 0.  The degree of the
!    zero polynomial is debatable, but this routine returns the
!    degree as 0.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NA, the dimension of A.
!
!    Input, real ( kind = rk ) A(0:NA), the coefficients of the polynomials.
!
!    Output, integer DEGREE, the degree of A.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer na

  real ( kind = rk ) a(0:na)
  integer degree

  degree = na

  do while ( 0 < degree )

    if ( a(degree) /= 0.0D+00 ) then
      return
    end if

    degree = degree - 1

  end do

  return
end
subroutine r8poly_print ( n, a, title )

!*****************************************************************************80
!
!! R8POLY_PRINT prints out a polynomial.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of A.
!
!    Input, real ( kind = rk ) A(0:N), the polynomial coefficients.
!    A(0) is the constant term and
!    A(N) is the coefficient of X^N.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(0:n)
  integer i
  real ( kind = rk ) mag
  integer n2
  character plus_minus
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  call r8poly_degree ( n, a, n2 )

  if ( n2 <= 0 ) then
    write ( *, '( ''  p(x) = 0'' )' )
    return
  end if

  if ( a(n2) < 0.0D+00 ) then
    plus_minus = '-'
  else
    plus_minus = ' '
  end if

  mag = abs ( a(n2) )

  if ( 2 <= n2 ) then
    write ( *, '( ''  p(x) = '', a1, g14.6, '' * x ^ '', i3 )' ) &
      plus_minus, mag, n2
  else if ( n2 == 1 ) then
    write ( *, '( ''  p(x) = '', a1, g14.6, '' * x'' )' ) &
      plus_minus, mag
  else if ( n2 == 0 ) then
    write ( *, '( ''  p(x) = '', a1, g14.6 )' ) plus_minus, mag
  end if

  do i = n2 - 1, 0, -1

    if ( a(i) < 0.0D+00 ) then
      plus_minus = '-'
    else
      plus_minus = '+'
    end if

    mag = abs ( a(i) )

    if ( mag /= 0.0D+00 ) then

      if ( 2 <= i ) then
        write ( *, ' ( ''         '', a1, g14.6, '' * x ^ '', i3 )' ) &
          plus_minus, mag, i
      else if ( i == 1 ) then
        write ( *, ' ( ''         '', a1, g14.6, '' * x'' )' ) plus_minus, mag
      else if ( i == 0 ) then
        write ( *, ' ( ''         '', a1, g14.6 )' ) plus_minus, mag
      end if
    end if

  end do

  return
end
subroutine r8vec_chebyshev ( n, a_first, a_last, a )

!*****************************************************************************80
!
!! R8VEC_CHEBYSHEV creates a vector of Chebyshev spaced values.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 June 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real ( kind = rk ) A_FIRST, A_LAST, the first and last entries.
!
!    Output, real ( kind = rk ) A(N), a vector of Chebyshev spaced data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) a_first
  real ( kind = rk ) a_last
  real ( kind = rk ) c
  integer i
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) theta

  if ( n == 1 ) then

    a(1) = ( a_first + a_last ) / 2.0D+00

  else

    do i = 1, n

      theta = real ( n - i, kind = rk ) * pi &
            / real ( n - 1, kind = rk )

      c = cos ( theta )

      a(i) = ( ( 1.0D+00 - c ) * a_first  &
             + ( 1.0D+00 + c ) * a_last ) &
             /   2.0D+00

    end do

  end if

  return
end
subroutine r8vec_linspace ( n, a_first, a_last, a )

!*****************************************************************************80
!
!! R8VEC_LINSPACE creates a vector of linearly spaced values.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
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
!    Input, real ( kind = rk ) A_FIRST, A_LAST, the first and last entries.
!
!    Output, real ( kind = rk ) A(N), a vector of linearly spaced data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) a_first
  real ( kind = rk ) a_last
  integer i

  if ( n == 1 ) then

    a(1) = ( a_first + a_last ) / 2.0D+00

  else

    do i = 1, n
      a(i) = ( real ( n - i,     kind = rk ) * a_first &
             + real (     i - 1, kind = rk ) * a_last ) &
             / real ( n     - 1, kind = rk )
    end do

  end if

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
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
!    Input, real ( kind = rk ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
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
!    06 August 2005
!
!  Author:
!
!    John Burkardt
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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
