function annulus_area ( center, r1, r2 )

!*****************************************************************************80
!
!! annulus_area() computes the area of a circular annulus in 2D.
!
!  Discussion:
!
!    A circular annulus with center (XC,YC), inner radius R1 and
!    outer radius R2, is the set of points (X,Y) so that
!
!      R1^2 <= (X-XC)^2 + (Y-YC)^2 <= R2^2
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 August 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) CENTER(2), the coordinates of the center.
!    This data is actually not necessary for area calculations.
!
!    real ( kind = rk ) R1, R2, the inner and outer radii of the annulus.
!    0.0 <= R1 <= R2.
!
!  Output:
!
!    real ( kind = rk ) ANNULUS_AREA, the area of the annulus.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) annulus_area
  real ( kind = rk ) center(2)
  real ( kind = rk ) r1
  real ( kind = rk ) r2
  real ( kind = rk ), parameter :: r8_pi = 3.1415926535897932384626434D+00
  real ( kind = rk ) value

  if ( r1 < 0.0D+00 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'ANNULUS_AREA - Fatal error!'
    write ( *, '(a)' ) '  Inner radius R1 < 0.0.'
    stop 1
  end if

  if ( r2 < r1 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'ANNULUS_AREA - Fatal error!'
    write ( *, '(a)' ) '  Outer radius R1 < R1 = inner radius.'
    stop 1
  end if

  value = r8_pi * ( r2 + r1 ) * ( r2 - r1 )

  annulus_area = value

  return
end
subroutine annulus_rule_compute ( center, r1, r2, nr, nt, w, x, y )

!*****************************************************************************80
!
!! annulus_rule_compute() computes a quadrature rule for an annulus.
!
!  Discussion:
!
!    The integration region is points (X,Y) such that
!
!      R1^2 <= ( X - CENTER(1) )^2 + ( Y - CENTER(2) )^2 <= R2^2
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 August 2021
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Peirce,
!    Numerical Integration Over the Planar Annulus,
!    Journal of the Society for Industrial and Applied Mathematics,
!    Volume 5, Issue 2, June 1957, pages 66-73.
!
!  Input:
!
!    real ( kind = rk ) CENTER(2), the coordinates of the center.
!
!    real ( kind = rk ) R1, R2, the inner and outer radii of the annulus.
!    0.0 <= R1 <= R2.
!
!    integer NR, the number of points in the radial rule.
!
!    integer NT, the number of angles to use.
!    The value NT=4*NR is recommended.
!
!  Output:
!
!    real ( kind = rk ) W(NR*NT), the weights for the rule.
!
!    real ( kind = rk ) X(NR*NT), Y(NR*NT), the points for the rule.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nr
  integer nt

  real ( kind = rk ) a
  real ( kind = rk ) annulus_area
  real ( kind = rk ) area
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) center(2)
  real ( kind = rk ) d
  integer i
  integer j
  integer k
  real ( kind = rk ) r1
  real ( kind = rk ) r2
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) ra(nr)
  real ( kind = rk ) rw(nr)
  real ( kind = rk ) t
  real ( kind = rk ) tw
  real ( kind = rk ) w(nr*nt)
  real ( kind = rk ) x(nr*nt)
  real ( kind = rk ) y(nr*nt)
!
!  Request a Legendre rule for [-1,+1].
!
  call legendre_ek_compute ( nr, ra, rw )
!
!  Adjust the rule from [-1,+1] to [r1^2,r2^2].
!
  a = -1.0D+00
  b = +1.0D+00
  c = r1 * r1
  d = r2 * r2
  call rule_adjust ( a, b, c, d, nr, ra, rw )
!
!  Convert from R^2 to R.
!
  ra(1:nr) = sqrt ( ra(1:nr) )
  rw(1:nr) = rw(1:nr) / ( r2 + r1 ) / ( r2 - r1 )
!
!  Set the angular weight.
!
  tw = 1.0D+00 / real ( nt, kind = rk )
!
!  Get area of annulus.
!
  area = annulus_area ( center, r1, r2 )
!
!  Form the abscissa and weight vectors.
!
  k = 0
  do i = 0, nt - 1
    t = 2.0D+00 * r8_pi * real ( i, kind = rk ) / real ( nt, kind = rk )
    do j = 1, nr
      k = k + 1
      x(k) = center(1) + ra(j) * cos ( t )
      y(k) = center(2) + ra(j) * sin ( t )
      w(k) = area * tw * rw(j)
    end do
  end do

  return
end
subroutine disk01_monomial_integral ( e, integral )

!*****************************************************************************80
!
!! disk01_monomial_integral() returns monomial integrals in the unit disk in 2D.
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
!    26 August 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer E(2), the exponents of X and Y in the 
!    monomial.  Each exponent must be nonnegative.
!
!  Output:
!
!    real ( kind = rk ) INTEGRAL, the integral.
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
subroutine imtqlx ( n, d, e, z )

!*****************************************************************************80
!
!! imtqlx() diagonalizes a symmetric tridiagonal matrix.
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
!    26 August 2021
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
!  Input:
!
!    integer N, the order of the matrix.
!
!    real ( kind = rk ) D(N), the diagonal entries of the matrix.
!
!    real ( kind = rk ) E(N), the subdiagonal entries of the
!    matrix, in entries E(1) through E(N-1).
!
!    real ( kind = rk ) Z(N), a vector.  
!
!  Output:
!
!    real ( kind = rk ) D(N), overwritten.
!
!    real ( kind = rk ) E(N), overwritten.
!
!    real ( kind = rk ) Z(N), the value of Q' * Z, where Q is the matrix 
!    that diagonalizes the input symmetric tridiagonal matrix.
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
!! legendre_ek_compute(): Legendre quadrature rule by the Elhay-Kautsky method.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 August 2021
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
!  Input:
!
!    integer N, the order.
!
!  Output:
!
!    real ( kind = rk ) X(N), the abscissas.
!
!    real ( kind = rk ) W(N), the weights.
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

  w(1:n) = w(1:n) ** 2

  return
end
subroutine monomial_value ( m, n, e, x, v )

!*****************************************************************************80
!
!! monomial_value() evaluates a monomial.
!
!  Discussion:
!
!    This routine evaluates a monomial of the form
!
!      product ( 1 <= i <= m ) x(i)^e(i)
!
!    The combination 0.0^0, if encountered, is treated as 1.0.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 August 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, the spatial dimension.
!
!    integer N, the number of evaluation points.
!
!    integer E(M), the exponents.
!
!    real ( kind = rk ) X(M,N), the point coordinates.
!
!  Output:
!
!    real ( kind = rk ) V(N), the monomial values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer e(m)
  integer i
  real ( kind = rk ) v(n)
  real ( kind = rk ) x(m,n)

  v(1:n) = 1.0D+00

  do i = 1, m
    if ( 0 /= e(i) ) then
      v(1:n) = v(1:n) * x(i,1:n) ** e(i)
    end if
  end do

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! r8vec_print() prints a real vector.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 August 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of components of the vector.
!
!    real ( kind = rk ) A(N), the vector to be printed.
!
!    character ( len = * ) TITLE, a title.
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
subroutine r8vec3_print ( n, a1, a2, a3, title )

!*****************************************************************************80
!
!! r8vec3_print() prints a trio of real vectors.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 August 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of components of the vector.
!
!    real ( kind = rk ) A1(N), A2(N), A3(N), the vectors to be printed.
!
!    character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a1(n)
  real ( kind = rk ) a2(n)
  real ( kind = rk ) a3(n)
  integer i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(i8,3g14.6)' ) i, a1(i), a2(i), a3(i)
  end do

  return
end
subroutine rule_adjust ( a, b, c, d, order, x, w )

!*****************************************************************************80
!
!! rule_adjust() maps a quadrature rule from [A,B] to [C,D].
!
!  Discussion:
!
!    Most quadrature rules are defined on a special interval, like
!    [-1,1] or [0,1].  To integrate over an interval, the abscissas
!    and weights must be adjusted.  This can be done on the fly,
!    or by calling this routine.
!
!    If the weight function W(X) is not 1, then the W vector will
!    require further adjustment by the user.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 August 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) A, B, the endpoints of the definition interval.
!
!    real ( kind = rk ) C, D, the endpoints of the integration interval.
!
!    integer ORDER, the number of abscissas and weights.
!
!    real ( kind = rk ) X(ORDER), W(ORDER), the abscissas and weights.
!
!  Output:
!
!    real ( kind = rk ) X(ORDER), W(ORDER), the adjusted abscissas and weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer order

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) d
  real ( kind = rk ) w(order)
  real ( kind = rk ) x(order)

  x(1:order) = ( ( b - x(1:order)     ) * c   &
               + (     x(1:order) - a ) * d ) &
               / ( b              - a )

  w(1:order) = ( ( d - c ) / ( b - a ) ) * w(1:order)

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! timestamp() prints the current YMDHMS date as a time stamp.
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
!    26 August 2021
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
 
