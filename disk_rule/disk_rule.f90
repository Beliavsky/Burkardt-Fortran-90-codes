subroutine disk_rule ( nr, nt, xc, yc, rc, w, x, y )

!*****************************************************************************80
!
!! disk_rule() computes a quadrature rule for a general disk.
!
!  Discussion:
!
!    The general disk is the region:
!
!      ( x - xc ) ^ 2 + ( y - yc ) ^ 2 <= rc ^ 2.
!
!    The integral I(f) is then approximated by
!
!      S(f) = sum ( 1 <= i <= NT * NR ) W(i) * F ( X(i), Y(i) ).
!
!      Area = pi * RC ^ 2
!
!      Q(f) = Area * S(f)
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 April 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NR, the number of points in the radial rule.
!
!    Input, integer NT, the number of angles to use.
!
!    Input, real ( kind = rk ) XC, YC, the center of the disk.
!
!    Input, real ( kind = rk ) RC, the radius of the disk.
!
!    Output, real ( kind = rk ) W(NR*NT), the weights for the rule.
!
!    Output, real ( kind = rk ) X(NR*NT), Y(NR*NT), the points for the rule.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nr
  integer nt

  integer j
  real ( kind = rk ) r01(nr)
  real ( kind = rk ) rc
  real ( kind = rk ) t01(nt)
  real ( kind = rk ) w01(nr)
  real ( kind = rk ) w(nr,nt)
  real ( kind = rk ) x(nr,nt)
  real ( kind = rk ) xc
  real ( kind = rk ) y(nr,nt)
  real ( kind = rk ) yc

  call disk01_rule ( nr, nt, w01, r01, t01 )
!
!  Recompute the rule for the general circle in terms of X, Y.
!
  do j = 1, nt
    w(1:nr,j) = w01(1:nr)
    x(1:nr,j) = xc + rc * r01(1:nr) * cos ( t01(j) )
    y(1:nr,j) = yc + rc * r01(1:nr) * sin ( t01(j) )
  end do

  return
end
subroutine disk01_monomial_integral ( e, integral )

!*****************************************************************************80
!
!! DISK01_MONOMIAL_INTEGRAL returns monomial integrals in the unit disk in 2D.
!
!  Discussion:
!
!    The integration region is 
!
!      X^2 + Y^2 <= 1.
!
!    The monomial is F(X,Y) = X^E(1) * Y^E(2).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 January 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer E(2), the exponents of X and Y in the 
!    monomial.  Each exponent must be nonnegative.
!
!    Output, real ( kind = rk ) INTEGRAL, the integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) arg
  integer e(2)
  integer i
  real ( kind = rk ) integral
  real ( kind = rk ), parameter :: r = 1.0D+00
  integer s

  if ( any ( e(1:2) < 0 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DISK01_MONOMIAL_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  All exponents must be nonnegative.'
    stop 1
  end if

  if ( any ( mod ( e(1:2), 2 ) == 1 ) ) then

    integral = 0.0D+00

  else

    integral = 2.0D+00

    do i = 1, 2
      arg = 0.5D+00 * real ( e(i) + 1, kind = rk )
      integral = integral * gamma ( arg )
    end do

    arg = 0.5D+00 * ( real ( sum ( e(1:2) + 1 ), kind = rk ) )
    integral = integral / gamma ( arg )

  end if
!
!  The surface integral is now adjusted to give the volume integral.
!
  s = sum ( e(1:2) ) + 2

  integral = integral * r ** s / real ( s, kind = rk )

  return
end
subroutine disk01_rule ( nr, nt, w, r, t )

!*****************************************************************************80
!
!! DISK01_RULE computes a quadrature rule for the unit disk.
!
!  Discussion:
!
!    The unit disk is the region:
!
!      x * x + y * y <= 1.
!
!    The integral I(f) is then approximated by
!
!      Q(f) = pi * sum ( 1 <= j <= NT ) sum ( 1 <= i <= NR ) 
!        W(i) * F ( R(i) * cos(T(j)), R(i) * sin(T(j)) ).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 March 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NR, the number of points in the radial rule.
!
!    Input, integer NT, the number of angles to use.
!
!    Output, real ( kind = rk ) W(NR), the weights for the rule.
!
!    Output, real ( kind = rk ) R(NR), T(NT), the (R,Theta) points for the rule.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nr
  integer nt

  integer ir
  integer it
  real ( kind = rk ) r(nr)
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) t(nt)
  real ( kind = rk ) w(nr)
  real ( kind = rk ) wr(nr)
  real ( kind = rk ) xr(nr)
!
!  Request a Legendre rule for [-1,+1].
!
  call legendre_ek_compute ( nr, xr, wr )
!
!  Shift the rule to [0,1].
!
  xr(1:nr) = ( xr(1:nr) + 1.0D+00 ) / 2.0D+00
  wr(1:nr) = wr(1:nr) / 2.0D+00
!
!  Compute the disk rule.
!
  do it = 1, nt
    t(it) = 2.0D+00 * r8_pi * real ( it - 1, kind = rk ) / real ( nt, kind = rk )
  end do

  do ir = 1, nr
    w(ir) = wr(ir) / real ( nt, kind = rk )
    r(ir) = sqrt ( xr(ir) )
  end do

  return
end
subroutine imtqlx ( n, d, e, z )

!*****************************************************************************80
!
!! IMTQLX diagonalizes a symmetric tridiagonal matrix.
!
!  Discussion:
!
!    This routine is a slightly modified version of the EISPACK routine to
!    perform the implicit QL algorithm on a symmetric tridiagonal matrix.
!
!    The authors thank the authors of EISPACK for permission to use this
!    routine.
!
!    It has been modified to produce the product Q' * Z, where Z is an input
!    vector and Q is the orthogonal matrix diagonalizing the input matrix.
!    The changes consist (essentially) of applying the orthogonal
!    transformations directly to Z as they are generated.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 December 2009
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!    Roger Martin, James Wilkinson,
!    The Implicit QL Algorithm,
!    Numerische Mathematik,
!    Volume 12, Number 5, December 1968, pages 377-383.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real ( kind = rk ) D(N), the diagonal entries of the matrix.
!    On output, the information in D has been overwritten.
!
!    Input/output, real ( kind = rk ) E(N), the subdiagonal entries of the
!    matrix, in entries E(1) through E(N-1).  On output, the information in
!    E has been overwritten.
!
!    Input/output, real ( kind = rk ) Z(N).  On input, a vector.  On output,
!    the value of Q' * Z, where Q is the matrix that diagonalizes the
!    input symmetric tridiagonal matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) d(n)
  real ( kind = rk ) e(n)
  real ( kind = rk ) f
  real ( kind = rk ) g
  integer i
  integer ii
  integer, parameter :: itn = 30
  integer j
  integer k
  integer l
  integer m
  integer mml
  real ( kind = rk ) p
  real ( kind = rk ) prec
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) z(n)

  prec = epsilon ( prec )

  if ( n == 1 ) then
    return
  end if

  e(n) = 0.0D+00

  do l = 1, n

    j = 0

    do

      do m = l, n

        if ( m == n ) then
          exit
        end if

        if ( abs ( e(m) ) <= prec * ( abs ( d(m) ) + abs ( d(m+1) ) ) ) then
          exit
        end if

      end do

      p = d(l)

      if ( m == l ) then
        exit
      end if

      if ( itn <= j ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IMTQLX - Fatal error!'
        write ( *, '(a)' ) '  Iteration limit exceeded.'
        write ( *, '(a,i8)' ) '  J = ', j
        write ( *, '(a,i8)' ) '  L = ', l
        write ( *, '(a,i8)' ) '  M = ', m
        write ( *, '(a,i8)' ) '  N = ', n
        stop
      end if

      j = j + 1
      g = ( d(l+1) - p ) / ( 2.0D+00 * e(l) )
      r =  sqrt ( g * g + 1.0D+00 )
      g = d(m) - p + e(l) / ( g + sign ( r, g ) )
      s = 1.0D+00
      c = 1.0D+00
      p = 0.0D+00
      mml = m - l

      do ii = 1, mml

        i = m - ii
        f = s * e(i)
        b = c * e(i)

        if ( abs ( g ) <= abs ( f ) ) then
          c = g / f
          r =  sqrt ( c * c + 1.0D+00 )
          e(i+1) = f * r
          s = 1.0D+00 / r
          c = c * s
        else
          s = f / g
          r =  sqrt ( s * s + 1.0D+00 )
          e(i+1) = g * r
          c = 1.0D+00 / r
          s = s * c
        end if

        g = d(i+1) - p
        r = ( d(i) - g ) * s + 2.0D+00 * c * b
        p = s * r
        d(i+1) = g + p
        g = c * r - b
        f = z(i+1)
        z(i+1) = s * z(i) + c * f
        z(i) = c * z(i) - s * f

      end do

      d(l) = d(l) - p
      e(l) = g
      e(m) = 0.0D+00

    end do

  end do
!
!  Sorting.
!
  do ii = 2, n

    i = ii - 1
    k = i
    p = d(i)

    do j = ii, n
      if ( d(j) < p ) then
        k = j
        p = d(j)
      end if
    end do

    if ( k /= i ) then
      d(k) = d(i)
      d(i) = p
      p = z(i)
      z(i) = z(k)
      z(k) = p
    end if

  end do

  return
end
subroutine legendre_ek_compute ( n, x, w )

!*****************************************************************************80
!
!! LEGENDRE_EK_COMPUTE: Legendre quadrature rule by the Elhay-Kautsky method.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 April 2011
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
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
!    Input, integer N, the order.
!
!    Output, real ( kind = rk ) X(N), the abscissas.
!
!    Output, real ( kind = rk ) W(N), the weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) bj(n)
  integer i
  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n)
  real ( kind = rk ) zemu
!
!  Define the zero-th moment.
!
  zemu = 2.0D+00
!
!  Define the Jacobi matrix.
!
  do i = 1, n
    bj(i) = real ( i * i, kind = rk ) / real ( 4 * i * i - 1, kind = rk )
  end do
  bj(1:n) = sqrt ( bj(1:n) )

  x(1:n) = 0.0D+00

  w(1) = sqrt ( zemu )
  w(2:n) = 0.0D+00
!
!  Diagonalize the Jacobi matrix.
!
  call imtqlx ( n, x, bj, w )

  w(1:n) = w(1:n)**2

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
