subroutine comp_next ( n, k, a, more, h, t )

!*****************************************************************************80
!
!! comp_next() computes the compositions of the integer N into K parts.
!
!  Discussion:
!
!    A composition of the integer N into K parts is an ordered sequence
!    of K nonnegative integers which sum to N.  The compositions (1,2,1)
!    and (1,1,2) are considered to be distinct.
!
!    The routine computes one composition on each call until there are no more.
!    For instance, one composition of 6 into 3 parts is
!    3+2+1, another would be 6+0+0.
!
!    On the first call to this routine, set MORE = FALSE.  The routine
!    will compute the first element in the sequence of compositions, and
!    return it, as well as setting MORE = TRUE.  If more compositions
!    are desired, call again, and again.  Each time, the routine will
!    return with a new composition.
!
!    However, when the LAST composition in the sequence is computed
!    and returned, the routine will reset MORE to FALSE, signaling that
!    the end of the sequence has been reached.
!
!    This routine originally used a SAVE statement to maintain the
!    variables H and T.  I have decided that it is safer
!    to pass these variables as arguments, even though the user should
!    never alter them.  This allows this routine to safely shuffle
!    between several ongoing calculations.
!
!
!    There are 28 compositions of 6 into three parts.  This routine will
!    produce those compositions in the following order:
!
!     I         A
!     -     ---------
!     1     6   0   0
!     2     5   1   0
!     3     4   2   0
!     4     3   3   0
!     5     2   4   0
!     6     1   5   0
!     7     0   6   0
!     8     5   0   1
!     9     4   1   1
!    10     3   2   1
!    11     2   3   1
!    12     1   4   1
!    13     0   5   1
!    14     4   0   2
!    15     3   1   2
!    16     2   2   2
!    17     1   3   2
!    18     0   4   2
!    19     3   0   3
!    20     2   1   3
!    21     1   2   3
!    22     0   3   3
!    23     2   0   4
!    24     1   1   4
!    25     0   2   4
!    26     1   0   5
!    27     0   1   5
!    28     0   0   6
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 July 2008
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Second Edition,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer N, the integer whose compositions are desired.
!
!    Input, integer K, the number of parts in the composition.
!
!    Input/output, integer A(K), the parts of the composition.
!
!    Input/output, logical MORE, set by the user to start the computation,
!    and by the routine to terminate it.
!
!    Input/output, integer  H, T, two internal parameters needed
!    for the computation.  The user should allocate space for these in the
!    calling program, include them in the calling sequence, but never alter
!    them!
!
  implicit none

  integer k

  integer a(k)
  integer h
  logical more
  integer n
  integer t
!
!  The first computation.
!
  if ( .not. more ) then

    t = n
    h = 0
    a(1) = n
    a(2:k) = 0
!
!  The next computation.
!
  else
!
!  If the first entry A(1) is positive, then set H to zero,
!  so that when we increment H, it points to A(1); we will decrement A(1) by 1
!  and increment A(2).
!
    if ( 1 < t ) then
      h = 0
    end if
!
!  Otherwise, A(1) is 0.  Then by H + 1 is the entry we incremented last time.
!  Set H = H + 1, zero A(H), adding all but one of its value to A(1),
!  and incrementing A(H+1) by 1.
!
    h = h + 1
    t = a(h)
    a(h) = 0
    a(1) = t - 1
    a(h+1) = a(h+1) + 1

  end if
!
!  This is the last element of the sequence if all the
!  items are in the last slot.
!
  more = ( a(k) /= n )

  return
end
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
subroutine hermite_compute ( n, x, w )

!*****************************************************************************80
!
!! HERMITE_COMPUTE computes a Gauss-Hermite quadrature rule.
!
!  Discussion:
!
!    The code uses an algorithm by Elhay and Kautsky.
!
!    The abscissas are the zeros of the N-th order Hermite polynomial.
!
!    The integral:
!
!      integral ( -oo < x < +oo ) exp ( - x * x ) * f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 April 2011
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
!    Input, integer N, the number of abscissas.
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
  real ( kind = rk ) r8_gamma
  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n)
  real ( kind = rk ) zemu
!
!  Define the zero-th moment.
!
  zemu = r8_gamma ( 1.0D+00 / 2.0D+00 )
!
!  Define the Jacobi matrix.
!
  do i = 1, n
    bj(i) = real ( i, kind = rk ) / 2.0D+00
  end do
  bj(1:n) = sqrt ( bj(1:n) )

  x(1:n) = 0.0D+00

  w(1) = sqrt ( zemu )
  w(2:n) = 0.0D+00
!
!  Diagonalize the Jacobi matrix.
!
  call imtqlx ( n, x, bj, w )
!
!  If N is odd, force the center value to be zero.
!
  if ( mod ( n, 2 ) == 1 ) then
    x((n+1)/2) = 0.0D+00
  end if

  w(1:n) = w(1:n)**2

  return
end
subroutine hermite_compute_points ( n, x )

!*****************************************************************************80
!
!! HERMITE_COMPUTE_POINTS computes points of a Hermite quadrature rule.
!
!  Discussion:
!
!    This is simply a convenient interface.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 August 2009
!
!  Author:
!
!    John Burkardt
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
!    Input, integer N, the order.
!    1 <= N.
!
!    Output, real ( kind = rk ) X(N), the abscissas.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n)

  call hermite_compute ( n, x, w )

  return
end
subroutine hermite_compute_weights ( n, w )

!*****************************************************************************80
!
!! HERMITE_COMPUTE_WEIGHTS computes weights of a Hermite quadrature rule.
!
!  Discussion:
!
!    This is simply a convenient interface.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 August 2009
!
!  Author:
!
!    John Burkardt
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
!    Input, integer N, the order.
!    1 <= ORDER.
!
!    Output, real ( kind = rk ) W(N), the weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n)

  call hermite_compute ( n, x, w )

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
        stop 1
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
subroutine level_to_order_hermite ( dim_num, level, order )

!*****************************************************************************80
!
!! LEVEL_TO_ORDER_HERMITE: default growth for Hermite sparse grids.
!
!  Discussion:
!
!    This function uses:
!
!    * exponential growth rates for fully nested rules,
!      ( "CC", "F2", "GP");
!
!    * linear growth rates for most rules:
!      ( "GL", "GH", "GGH", "LG", "GLG", "GJ", "GW" ).
!
!    * slow exponential growth alternative for fully nested rules:
!      ("CC_SE", "F2_SE", "GP_SE").
!
!    * moderate exponential growth alternative for fully nested rules:
!      ("CC_ME", "F2_ME", "GP_ME").
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer LEVEL(DIM_NUM), the 1D levels.
!
!    Output, integer ORDER(DIM_NUM), the 1D orders
!    (number of points).
!
  implicit none

  integer dim_num

  integer dim
  integer level(dim_num)
  integer order(dim_num)

  do dim = 1, dim_num

    order(dim) = 2 ** ( level(dim) + 1 ) - 1

  end do

  return
end
subroutine monomial_integral_hermite ( dim_num, expon, value )

!*****************************************************************************80
!
!! MONOMIAL_INTEGRAL_HERMITE integrates a Hermite monomial.
!
!  Discussion:
!
!    H(d,n) = Integral ( -oo < x < +oo )
!      x1^n1 * x2^n2...*xd^nd * exp(-x1^2-x2^2...-xd^2 ) dx
!
!    H(d,n) is 0 if any n(i) odd.
!
!    H(d,n) = product ( 1 <= i <= d )
!      ( (n(i)-1)!! * sqrt(pi) / 2^(n(i)/2) for all n(i) even.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the dimension of the integral.
!
!    Input, integer EXPON(DIM_NUM), the order of the integral.
!    0 <= EXPON(1:DIM_NUM).
!
!    Output, real ( kind = rk ) VALUE, the value of the integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num

  integer dim
  integer expon(dim_num)
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) r8_factorial2
  real ( kind = rk ) r8_huge
  real ( kind = rk ) value

  if ( any ( expon(1:dim_num) < 0 ) ) then

    value = - r8_huge ( )

  else if ( any ( mod ( expon(1:dim_num), 2 ) == 1 ) ) then

    value = 0.0D+00

  else

    value = 1.0D+00
    do dim = 1, dim_num
      value = value * r8_factorial2 ( expon(dim) - 1 ) * sqrt ( pi ) &
        / 2.0D+00**( expon(dim) / 2 )
    end do

  end if

  return
end
subroutine monomial_quadrature_hermite ( dim_num, expon, point_num, weight, &
  x, exact, quad, quad_error )

!*****************************************************************************80
!
!! MONOMIAL_QUADRATURE_HERMITE applies a quadrature rule to a monomial.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer EXPON(DIM_NUM), the exponents.
!
!    Input, integer POINT_NUM, the number of points in the rule.
!
!    Input, real ( kind = rk ) WEIGHT(POINT_NUM), the quadrature weights.
!
!    Input, real ( kind = rk ) X(DIM_NUM,POINT_NUM), the quadrature points.
!
!    Output, real ( kind = rk ) EXACT, the exact integral.
!
!    Output, real ( kind = rk ) QUAD, the estimated integral.
!
!    Output, real ( kind = rk ) QUAD_ERROR, the relative quadrature error.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num

  real ( kind = rk ) exact
  integer expon(dim_num)
  integer point_num
  real ( kind = rk ) quad
  real ( kind = rk ) quad_error
  real ( kind = rk ) value(point_num)
  real ( kind = rk ) weight(point_num)
  real ( kind = rk ) x(dim_num,point_num)
!
!  Get the exact value of the integral of the unscaled monomial.
!
  call monomial_integral_hermite ( dim_num, expon, exact )
!
!  Evaluate the monomial at the quadrature points.
!
  call monomial_value ( dim_num, point_num, x, expon, value )
!
!  Compute the weighted sum.
!
  quad = dot_product ( weight, value )
!
!  If exact value is nonzero, use it to scale the data.
!
  if ( exact == 0.0D+00 ) then
    quad_error = abs ( quad )
  else
    quad_error = abs ( ( quad - exact ) / exact )
  end if

  return
end
subroutine monomial_value ( dim_num, point_num, x, expon, value )

!*****************************************************************************80
!
!! MONOMIAL_VALUE evaluates a monomial.
!
!  Discussion:
!
!    This routine evaluates a monomial of the form
!
!      product ( 1 <= dim <= dim_num ) x(dim)^expon(dim)
!
!    where the exponents are nonnegative integers.  Note that
!    if the combination 0^0 is encountered, it should be treated
!    as 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 October 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer POINT_NUM, the number of points at which the
!    monomial is to be evaluated.
!
!    Input, real ( kind = rk ) X(DIM_NUM,POINT_NUM), the point coordinates.
!
!    Input, integer EXPON(DIM_NUM), the exponents.
!
!    Output, real ( kind = rk ) VALUE(POINT_NUM), the value of the monomial.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num
  integer point_num

  integer dim
  integer expon(dim_num)
  integer point
  real ( kind = rk ) value(point_num)
  real ( kind = rk ) x(dim_num,point_num)

  value(1:point_num) = 1.0D+00

  do dim = 1, dim_num
    do point = 1, point_num
      if ( x(dim,point) /= 0.0D+00 ) then
        value(point) = value(point) * x(dim,point)**expon(dim)
      else if ( expon(dim) == 0 ) then
        value(point) = value(point)
      else
        value(point) = 0.0D+00
      end if
    end do
  end do

  return
end
subroutine product_hermite_weight ( dim_num, order_1d, order_nd, weight_nd )

!*****************************************************************************80
!
!! PRODUCT_HERMITE_WEIGHT computes the weights of a Hermite product rule.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer ORDER_1D(DIM_NUM), the order of the 1D rules.
!
!    Input, integer ORDER_ND, the order of the product rule.
!
!    Output, real ( kind = rk ) WEIGHT_ND(ORDER_ND), the product rule weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim
  integer dim_num
  integer order_nd

  integer order_1d(dim_num)
  real ( kind = rk ), allocatable, dimension ( : ) :: weight_1d
  real ( kind = rk ) weight_nd(order_nd)

  weight_nd(1:order_nd) = 1.0D+00

  do dim = 1, dim_num

    allocate ( weight_1d(1:order_1d(dim) ) )

    call hermite_compute_weights ( order_1d(dim), weight_1d )

    call r8vec_direct_product2 ( dim, order_1d(dim), weight_1d, &
      dim_num, order_nd, weight_nd )

    deallocate ( weight_1d )

  end do

  return
end
function r8_choose ( n, k )

!*****************************************************************************80
!
!! R8_CHOOSE computes the binomial coefficient C(N,K) as an R8.
!
!  Discussion:
!
!    The value is calculated in such a way as to avoid overflow and
!    roundoff.  The calculation is done in R8 arithmetic.
!
!    The formula used is:
!
!      C(N,K) = N! / ( K! * (N-K)! )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    ML Wolfson, HV Wright,
!    Algorithm 160:
!    Combinatorial of M Things Taken N at a Time,
!    Communications of the ACM,
!    Volume 6, Number 4, April 1963, page 161.
!
!  Parameters:
!
!    Input, integer N, K, are the values of N and K.
!
!    Output, real ( kind = rk ) R8_CHOOSE, the number of combinations of N
!    things taken K at a time.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer k
  integer mn
  integer mx
  integer n
  real ( kind = rk ) r8_choose
  real ( kind = rk ) value

  mn = min ( k, n - k )

  if ( mn < 0 ) then

    value = 0.0D+00

  else if ( mn == 0 ) then

    value = 1.0D+00

  else

    mx = max ( k, n - k )
    value = real ( mx + 1, kind = rk )

    do i = 2, mn
      value = ( value * real ( mx + i, kind = rk ) ) / real ( i, kind = rk )
    end do

  end if

  r8_choose = value

  return
end
function r8_factorial2 ( n )

!*****************************************************************************80
!
!! R8_FACTORIAL2 computes the double factorial function.
!
!  Formula:
!
!    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
!                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
!
!  Example:
!
!     N    Factorial2(N)
!
!     0     1
!     1     1
!     2     2
!     3     3
!     4     8
!     5    15
!     6    48
!     7   105
!     8   384
!     9   945
!    10  3840
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the argument of the double factorial
!    function.  If N is less than 1, R8_FACTORIAL2 is returned as 1.0.
!
!    Output, real ( kind = rk ) R8_FACTORIAL2, the value of the double
!    factorial of N.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  real ( kind = rk ) r8_factorial2
  real ( kind = rk ) r8_n

  if ( n < 1 ) then
    r8_factorial2 = 1.0D+00
    return
  end if

  r8_n = real ( n, kind = rk )
  r8_factorial2 = 1.0D+00

  do while ( 1.0D+00 < r8_n )
    r8_factorial2 = r8_factorial2 * r8_n
    r8_n = r8_n - 2.0D+00
  end do

  return
end
function r8_gamma ( x )

!*****************************************************************************80
!
!! R8_GAMMA evaluates Gamma(X) for a real argument.
!
!  Discussion:
!
!    This routine calculates the gamma function for a real argument X.
!
!    Computation is based on an algorithm outlined in reference 1.
!    The program uses rational functions that approximate the gamma
!    function to at least 20 significant decimal digits.  Coefficients
!    for the approximation over the interval (1,2) are unpublished.
!    Those for the approximation for 12 <= X are from reference 2.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody,
!    An Overview of Software Development for Special Functions,
!    in Numerical Analysis Dundee, 1975,
!    edited by GA Watson,
!    Lecture Notes in Mathematics 506,
!    Springer, 1976.
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thatcher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, the argument of the function.
!
!    Output, real ( kind = rk ) R8_GAMMA, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

!
!  Coefficients for minimax approximation over (12, INF).
!
  real ( kind = rk ), dimension ( 7 ) :: c = (/ &
   -1.910444077728D-03, &
    8.4171387781295D-04, &
   -5.952379913043012D-04, &
    7.93650793500350248D-04, &
   -2.777777777777681622553D-03, &
    8.333333333333333331554247D-02, &
    5.7083835261D-03 /)
  real ( kind = rk ), parameter :: eps = 2.22D-16
  real ( kind = rk ) fact
  real ( kind = rk ), parameter :: half = 0.5D+00
  integer i
  integer n
  real ( kind = rk ), parameter :: one = 1.0D+00
  real ( kind = rk ), dimension ( 8 ) :: p = (/ &
    -1.71618513886549492533811D+00, &
     2.47656508055759199108314D+01, &
    -3.79804256470945635097577D+02, &
     6.29331155312818442661052D+02, &
     8.66966202790413211295064D+02, &
    -3.14512729688483675254357D+04, &
    -3.61444134186911729807069D+04, &
     6.64561438202405440627855D+04 /)
  logical parity
  real ( kind = rk ), parameter :: pi = 3.1415926535897932384626434D+00
  real ( kind = rk ), dimension ( 8 ) :: q = (/ &
    -3.08402300119738975254353D+01, &
     3.15350626979604161529144D+02, &
    -1.01515636749021914166146D+03, &
    -3.10777167157231109440444D+03, &
     2.25381184209801510330112D+04, &
     4.75584627752788110767815D+03, &
    -1.34659959864969306392456D+05, &
    -1.15132259675553483497211D+05 /)
  real ( kind = rk ) r8_gamma
  real ( kind = rk ) res
  real ( kind = rk ), parameter :: sqrtpi = 0.9189385332046727417803297D+00
  real ( kind = rk ) sum
  real ( kind = rk ), parameter :: twelve = 12.0D+00
  real ( kind = rk ), parameter :: two = 2.0D+00
  real ( kind = rk ) x
  real ( kind = rk ), parameter :: xbig = 171.624D+00
  real ( kind = rk ) xden
  real ( kind = rk ), parameter :: xinf = 1.0D+30
  real ( kind = rk ), parameter :: xminin = 2.23D-308
  real ( kind = rk ) xnum
  real ( kind = rk ) y
  real ( kind = rk ) y1
  real ( kind = rk ) ysq
  real ( kind = rk ) z
  real ( kind = rk ), parameter :: zero = 0.0D+00

  parity = .false.
  fact = one
  n = 0
  y = x
!
!  Argument is negative.
!
  if ( y <= zero ) then

    y = - x
    y1 = aint ( y )
    res = y - y1

    if ( res /= zero ) then

      if ( y1 /= aint ( y1 * half ) * two ) then
        parity = .true.
      end if

      fact = - pi / sin ( pi * res )
      y = y + one

    else

      res = xinf
      r8_gamma = res
      return

    end if

  end if
!
!  Argument is positive.
!
  if ( y < eps ) then
!
!  Argument < EPS.
!
    if ( xminin <= y ) then
      res = one / y
    else
      res = xinf
      r8_gamma = res
      return
    end if

  else if ( y < twelve ) then

    y1 = y
!
!  0.0 < argument < 1.0.
!
    if ( y < one ) then

      z = y
      y = y + one
!
!  1.0 < argument < 12.0.
!  Reduce argument if necessary.
!
    else

      n = int ( y ) - 1
      y = y - real ( n, kind = rk )
      z = y - one

    end if
!
!  Evaluate approximation for 1.0 < argument < 2.0.
!
    xnum = zero
    xden = one
    do i = 1, 8
      xnum = ( xnum + p(i) ) * z
      xden = xden * z + q(i)
    end do

    res = xnum / xden + one
!
!  Adjust result for case  0.0 < argument < 1.0.
!
    if ( y1 < y ) then

      res = res / y1
!
!  Adjust result for case 2.0 < argument < 12.0.
!
    else if ( y < y1 ) then

      do i = 1, n
        res = res * y
        y = y + one
      end do

    end if

  else
!
!  Evaluate for 12.0 <= argument.
!
    if ( y <= xbig ) then

      ysq = y * y
      sum = c(7)
      do i = 1, 6
        sum = sum / ysq + c(i)
      end do
      sum = sum / y - y + sqrtpi
      sum = sum + ( y - half ) * log ( y )
      res = exp ( sum )

    else

      res = xinf
      r8_gamma = res
      return

    end if

  end if
!
!  Final adjustments and return.
!
  if ( parity ) then
    res = - res
  end if

  if ( fact /= one ) then
    res = fact / res
  end if

  r8_gamma = res

  return
end
function r8_huge ( )

!*****************************************************************************80
!
!! R8_HUGE returns a very large R8.
!
!  Discussion:
!
!    The value returned by this function is NOT required to be the
!    maximum representable R8.  This value varies from machine to machine,
!    from compiler to compiler, and may cause problems when being printed.
!    We simply want a "very large" but non-infinite number.
!
!    FORTRAN90 provides a built-in routine HUGE ( X ) that
!    can return the maximum representable number of the same datatype
!    as X, if that is what is really desired.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 October 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) R8_HUGE, a "huge" value.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r8_huge

  r8_huge = 1.0D+30

  return
end
function r8_mop ( i )

!*****************************************************************************80
!
!! R8_MOP returns the I-th power of -1 as an R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = rk ) value.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, the power of -1.
!
!    Output, real ( kind = rk ) R8_MOP, the I-th power of -1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  real ( kind = rk ) r8_mop

  if ( mod ( i, 2 ) == 0 ) then
    r8_mop = + 1.0D+00
  else
    r8_mop = - 1.0D+00
  end if

  return
end
subroutine r8col_compare ( m, n, a, i, j, value )

!*****************************************************************************80
!
!! R8COL_COMPARE compares columns in an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, I = 2, J = 4
!
!      A = (
!        1.  2.  3.  4.
!        5.  6.  7.  8.
!        9. 10. 11. 12. )
!
!    Output:
!
!      VALUE = -1
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = rk ) A(M,N), the M by N array.
!
!    Input, integer I, J, the columns to be compared.
!    I and J must be between 1 and N.
!
!    Output, integer VALUE, the results of the comparison:
!    -1, column I < column J,
!     0, column I = column J,
!    +1, column J < column I.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  integer i
  integer j
  integer k
  integer value
!
!  Check.
!
  if ( i < 1 .or. n < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8COL_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Column index I is out of bounds.'
    write ( *, '(a,i8)' ) '  I = ', i
    stop 1
  end if

  if ( j < 1 .or. n < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8COL_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Column index J is out of bounds.'
    write ( *, '(a,i8)' ) '  J = ', j
    stop 1
  end if

  value = 0

  if ( i == j ) then
    return
  end if

  k = 1

  do while ( k <= m )

    if ( a(k,i) < a(k,j) ) then
      value = -1
      return
    else if ( a(k,j) < a(k,i) ) then
      value = +1
      return
    end if

    k = k + 1

  end do

  return
end
subroutine r8col_sort_heap_a ( m, n, a )

!*****************************************************************************80
!
!! R8COL_SORT_HEAP_A ascending heapsorts an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!    In lexicographic order, the statement "X < Y", applied to two real
!    vectors X and Y of length M, means that there is some index I, with
!    1 <= I <= M, with the property that
!
!      X(J) = Y(J) for J < I,
!    and
!      X(I) < Y(I).
!
!    In other words, the first time they differ, X is smaller.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input/output, real ( kind = rk ) A(M,N).
!    On input, the array of N columns of M-vectors.
!    On output, the columns of A have been sorted in lexicographic order.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  integer i
  integer indx
  integer isgn
  integer j

  if ( m <= 0 ) then
    return
  end if

  if ( n <= 1 ) then
    return
  end if
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      call r8col_swap ( m, n, a, i, j )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call r8col_compare ( m, n, a, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine r8col_sort_heap_index_a ( m, n, a, indx )

!*****************************************************************************80
!
!! R8COL_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    A(*,J1) < A(*,J2) if the first nonzero entry of A(*,J1)-A(*,J2)
!    is negative.
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      A(*,INDX(*)) is sorted,
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows in each column of A.
!
!    Input, integer N, the number of columns in A.
!
!    Input, real ( kind = rk ) A(M,N), the array.
!
!    Output, integer INDX(N), the sort index.  The I-th element
!    of the sorted array is column INDX(I).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) column(m)
  integer i
  integer indx(n)
  integer indxt
  integer ir
  integer isgn
  integer j
  integer l

  if ( n < 1 ) then
    return
  end if

  do i = 1, n
    indx(i) = i
  end do

  if ( n == 1 ) then
    return
  end if

  l = ( n / 2 ) + 1
  ir = n

  do

    if ( 1 < l ) then

      l = l - 1
      indxt = indx(l)
      column(1:m) = a(1:m,indxt)

    else

      indxt = indx(ir)
      column(1:m) = a(1:m,indxt)
      indx(ir) = indx(1)
      ir = ir - 1

      if ( ir == 1 ) then
        indx(1) = indxt
        exit
      end if

    end if

    i = l
    j = l + l

    do while ( j <= ir )

      if ( j < ir ) then

        call r8vec_compare ( m, a(1:m,indx(j)), a(1:m,indx(j+1)), isgn )

        if ( isgn < 0 ) then
          j = j + 1
        end if

      end if

      call r8vec_compare ( m, column, a(1:m,indx(j)), isgn )

      if ( isgn < 0 ) then
        indx(i) = indx(j)
        i = j
        j = j + j
      else
        j = ir + 1
      end if

    end do

    indx(i) = indxt

  end do

  return
end
subroutine r8col_sorted_unique_count ( m, n, a, tol, unique_num )

!*****************************************************************************80
!
!! R8COL_SORTED_UNIQUE_COUNT counts unique elements in a sorted R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!    The columns of the array may be ascending or descending sorted.
!
!    If the tolerance is large enough, then the concept of uniqueness
!    can become ambiguous.  If we have a tolerance of 1.5, then in the
!    list ( 1, 2, 3, 4, 5, 6, 7, 8, 9 ) is it fair to say we have only
!    one unique entry?  That would be because 1 may be regarded as unique,
!    and then 2 is too close to 1 to be unique, and 3 is too close to 2 to
!    be unique and so on.
!
!    This seems wrongheaded.  So I prefer the idea that an item is not
!    unique under a tolerance only if it is close to something that IS unique.
!    Thus, the unique items are guaranteed to cover the space if we include
!    a disk of radius TOL around each one.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = rk ) A(M,N), a sorted array, containing
!    N columns of data.
!
!    Input, real ( kind = rk ) TOL, a tolerance for equality.
!
!    Output, integer UNIQUE_NUM, the number of unique columns.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  integer j1
  integer j2
  real ( kind = rk ) tol
  integer unique_num

  unique_num = 0

  if ( n <= 0 ) then
    return
  end if

  unique_num = 1
  j1 = 1

  do j2 = 2, n

    if ( tol < maxval ( abs ( a(1:m,j1) - a(1:m,j2) ) ) ) then
      unique_num = unique_num + 1
      j1 = j2
    end if

  end do

  return
end
subroutine r8col_swap ( m, n, a, j1, j2 )

!*****************************************************************************80
!
!! R8COL_SWAP swaps columns I and J of an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, J1 = 2, J2 = 4
!
!      A = (
!        1.  2.  3.  4.
!        5.  6.  7.  8.
!        9. 10. 11. 12. )
!
!    Output:
!
!      A = (
!        1.  4.  3.  2.
!        5.  8.  7.  6.
!        9. 12. 11. 10. )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input/output, real ( kind = rk ) A(M,N), the M by N array.
!
!    Input, integer J1, J2, the columns to be swapped.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) col(m)
  integer j1
  integer j2

  if ( j1 < 1 .or. n < j1 .or. j2 < 1 .or. n < j2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8COL_SWAP - Fatal error!'
    write ( *, '(a)' ) '  J1 or J2 is out of bounds.'
    write ( *, '(a,i8)' ) '  J1 =    ', j1
    write ( *, '(a,i8)' ) '  J2 =    ', j2
    write ( *, '(a,i8)' ) '  NCOL = ', n
    stop 1
  end if

  if ( j1 == j2 ) then
    return
  end if

  col(1:m) = a(1:m,j1)
  a(1:m,j1) = a(1:m,j2)
  a(1:m,j2) = col(1:m)

  return
end
subroutine r8col_tol_unique_count ( m, n, a, tol, unique_num )

!*****************************************************************************80
!
!! R8COL_TOL_UNIQUE_COUNT counts tolerably unique entries in an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!    If the tolerance is large enough, then the concept of uniqueness
!    can become ambiguous.  If we have a tolerance of 1.5, then in the
!    list ( 1, 2, 3, 4, 5, 6, 7, 8, 9 ) is it fair to say we have only
!    one unique entry?  That would be because 1 may be regarded as unique,
!    and then 2 is too close to 1 to be unique, and 3 is too close to 2 to
!    be unique and so on.
!
!    This seems wrongheaded.  So I prefer the idea that an item is not
!    unique under a tolerance only if it is close to something that IS unique.
!    Thus, the unique items are guaranteed to cover the space if we include
!    a disk of radius TOL around each one.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 July 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows.
!
!    Input, integer N, the number of columns.
!
!    Input, real ( kind = rk ) A(M,N), the array of N columns of data.
!
!    Input, real ( kind = rk ) TOL, a nonnegative tolerance for equality.
!
!    Output, integer UNIQUE_NUM, the number of unique columns.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) diff
  integer i
  integer indx(n)
  integer j
  integer k
  real ( kind = rk ) tol
  integer undx(n)
  logical unique
  integer unique_num
!
!  Implicitly sort the array.
!
  call r8col_sort_heap_index_a ( m, n, a, indx )
!
!  Consider entry I = 1.
!  It is unique, so set the number of unique items to K.
!  Set the K-th unique item to I.
!  Set the representative of item I to the K-th unique item.
!
  i = 1
  k = 1
  undx(k) = indx(i)
!
!  Consider entry I.
!
!  If it is unique, increase the unique count K, set the
!  K-th unique item to I, and set the representative of I to K.
!
!  If it is not unique, set the representative of item I to a
!  previously determined unique item that is close to it.
!
  do i = 2, n

    unique = .true.

    do j = 1, k
      diff = maxval ( abs ( a(1:m,indx(i)) - a(1:m,undx(j)) ) )
      if ( diff <= tol ) then
        unique = .false.
        exit
      end if
    end do

    if ( unique ) then
      k = k + 1
      undx(k) = indx(i)
    end if

  end do

  unique_num = k

  return
end
subroutine r8col_undex ( x_dim, x_num, x_val, x_unique_num, tol, undx, xdnu )

!*****************************************************************************80
!
!! R8COL_UNDEX returns unique sorted indexes for an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!    The goal of this routine is to determine a vector UNDX,
!    which points, to the unique elements of X, in sorted order,
!    and a vector XDNU, which identifies, for each entry of X, the index of
!    the unique sorted element of X.
!
!    This is all done with index vectors, so that the elements of
!    X are never moved.
!
!    The first step of the algorithm requires the indexed sorting
!    of X, which creates arrays INDX and XDNI.  (If all the entries
!    of X are unique, then these arrays are the same as UNDX and XDNU.)
!
!    We then use INDX to examine the entries of X in sorted order,
!    noting the unique entries, creating the entries of XDNU and
!    UNDX as we go.
!
!    Once this process has been completed, the object X could be
!    replaced by a compressed object XU, containing the unique entries
!    of X in sorted order, using the formula
!
!      XU(*) = X(UNDX(*)).
!
!    We could then, if we wished, reconstruct the entire vector X, or
!    any element of it, by index, as follows:
!
!      X(I) = XU(XDNU(I)).
!
!    We could then replace X by the combination of XU and XDNU.
!
!    Later, when we need the I-th entry of X, we can locate it as
!    the XDNU(I)-th entry of XU.
!
!    Here is an example of a vector X, the sort and inverse sort
!    index vectors, and the unique sort and inverse unique sort vectors
!    and the compressed unique sorted vector.
!
!      I    X   Indx  Xdni      XU   Undx  Xdnu
!    ----+-----+-----+-----+--------+-----+-----+
!      1 | 11.     1     1 |    11.     1     1
!      2 | 22.     3     5 |    22.     2     2
!      3 | 11.     6     2 |    33.     4     1
!      4 | 33.     9     8 |    55.     5     3
!      5 | 55.     2     9 |                  4
!      6 | 11.     7     3 |                  1
!      7 | 22.     8     6 |                  2
!      8 | 22.     4     7 |                  2
!      9 | 11.     5     4 |                  1
!
!    INDX(2) = 3 means that sorted item(2) is X(3).
!    XDNI(2) = 5 means that X(2) is sorted item(5).
!
!    UNDX(3) = 4 means that unique sorted item(3) is at X(4).
!    XDNU(8) = 2 means that X(8) is at unique sorted item(2).
!
!    XU(XDNU(I))) = X(I).
!    XU(I)        = X(UNDX(I)).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 November 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer X_DIM, the dimension of the data values.
!    (the number of rows in the R8COL).
!
!    Input, integer X_NUM, the number of data values.
!    (the number of columns in the R8COL).
!
!    Input, real ( kind = rk ) X_VAL(X_DIM,X_NUM), the data values.
!
!    Input, integer X_UNIQUE_NUM, the number of unique values
!    in X_VAL.  This value is only required for languages in which the size of
!    UNDX must be known in advance.
!
!    Input, real ( kind = rk ) TOL, a tolerance for equality.
!
!    Output, integer UNDX(X_UNIQUE_NUM), the UNDX vector.
!
!    Output, integer XDNU(X_NUM), the XDNU vector.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer x_dim
  integer x_num
  integer x_unique_num

  integer i
  integer indx(x_num)
  integer j
  real ( kind = rk ) tol
  integer undx(x_unique_num)
  real ( kind = rk ) x_val(x_dim,x_num)
  integer xdnu(x_num)
!
!  Implicitly sort the array.
!
  call r8col_sort_heap_index_a ( x_dim, x_num, x_val, indx )
!
!  Walk through the implicitly sorted array X.
!
  i = 1

  j = 1
  undx(j) = indx(i)

  xdnu(indx(i)) = j

  do i = 2, x_num

    if ( tol < &
         maxval ( abs ( x_val(1:x_dim,indx(i)) - x_val(1:x_dim,undx(j)) ) ) &
    ) then
      j = j + 1
      undx(j) = indx(i)
    end if

    xdnu(indx(i)) = j

  end do

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 September 2004
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
!    Input, real ( kind = rk ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  character ( len = * ) title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_PRINT_SOME prints some of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = rk ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ILO, JLO, the first row and column to print.
!
!    Input, integer IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: incx = 5
  integer m
  integer n

  real ( kind = rk ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer i
  integer i2hi
  integer i2lo
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)'
    return
  end if

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == real ( int ( a(i,j) ), kind = rk ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine r8mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = rk ) A(M,N), an M by N matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  character ( len = * ) title

  call r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = rk ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ILO, JLO, the first row and column to print.
!
!    Input, integer IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: incx = 5
  integer m
  integer n

  real ( kind = rk ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer i
  integer i2
  integer i2hi
  integer i2lo
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)'
    return
  end if

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(i8,6x)' ) i
    end do

    write ( *, '(''  Row   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc
        i = i2lo - 1 + i2
        write ( ctemp(i2), '(g14.6)' ) a(i,j)
      end do

      write ( *, '(i5,a,5a14)' ) j, ':', ( ctemp(i), i = 1, inc )

    end do

  end do

  return
end
subroutine r8mat_write ( output_filename, m, n, table )

!*****************************************************************************80
!
!! R8MAT_WRITE writes an R8MAT file.
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
!  For greater precision in the output file, try:
!
!                                            '(', m, 'g', 24, '.', 16, ')'
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
subroutine r8vec_compare ( n, a1, a2, isgn )

!*****************************************************************************80
!
!! R8VEC_COMPARE compares two R8VEC's.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The lexicographic ordering is used.
!
!  Example:
!
!    Input:
!
!      A1 = ( 2.0, 6.0, 2.0 )
!      A2 = ( 2.0, 8.0, 12.0 )
!
!    Output:
!
!      ISGN = -1
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vectors.
!
!    Input, real ( kind = rk ) A1(N), A2(N), the vectors to be compared.
!
!    Output, integer ISGN, the results of the comparison:
!    -1, A1 < A2,
!     0, A1 = A2,
!    +1, A1 > A2.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a1(n)
  real ( kind = rk ) a2(n)
  integer isgn
  integer k

  isgn = 0

  k = 1

  do while ( k <= n )

    if ( a1(k) < a2(k) ) then
      isgn = -1
      return
    else if ( a2(k) < a1(k) ) then
      isgn = + 1
      return
    end if

    k = k + 1

  end do

  return
end
subroutine r8vec_direct_product2 ( factor_index, factor_order, factor_value, &
  factor_num, point_num, w )

!*****************************************************************************80
!
!! R8VEC_DIRECT_PRODUCT2 creates a direct product of R8VEC's.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    To explain what is going on here, suppose we had to construct
!    a multidimensional quadrature rule as the product of K rules
!    for 1D quadrature.
!
!    The product rule will be represented as a list of points and weights.
!
!    The J-th item in the product rule will be associated with
!      item J1 of 1D rule 1,
!      item J2 of 1D rule 2,
!      ...,
!      item JK of 1D rule K.
!
!    In particular,
!      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
!    and
!      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)
!
!    So we can construct the quadrature rule if we can properly
!    distribute the information in the 1D quadrature rules.
!
!    This routine carries out the task involving the weights W.
!
!    Another way to do this would be to compute, one by one, the
!    set of all possible indices (J1,J2,...,JK), and then index
!    the appropriate information.  An advantage of the method shown
!    here is that you can process the K-th set of information and
!    then discard it.
!
!  Example:
!
!    Rule 1:
!      Order = 4
!      W(1:4) = ( 2, 3, 5, 7 )
!
!    Rule 2:
!      Order = 3
!      W(1:3) = ( 11, 13, 17 )
!
!    Rule 3:
!      Order = 2
!      W(1:2) = ( 19, 23 )
!
!    Product Rule:
!      Order = 24
!      W(1:24) =
!        ( 2 * 11 * 19 )
!        ( 3 * 11 * 19 )
!        ( 4 * 11 * 19 )
!        ( 7 * 11 * 19 )
!        ( 2 * 13 * 19 )
!        ( 3 * 13 * 19 )
!        ( 5 * 13 * 19 )
!        ( 7 * 13 * 19 )
!        ( 2 * 17 * 19 )
!        ( 3 * 17 * 19 )
!        ( 5 * 17 * 19 )
!        ( 7 * 17 * 19 )
!        ( 2 * 11 * 23 )
!        ( 3 * 11 * 23 )
!        ( 5 * 11 * 23 )
!        ( 7 * 11 * 23 )
!        ( 2 * 13 * 23 )
!        ( 3 * 13 * 23 )
!        ( 5 * 13 * 23 )
!        ( 7 * 13 * 23 )
!        ( 2 * 17 * 23 )
!        ( 3 * 17 * 23 )
!        ( 5 * 17 * 23 )
!        ( 7 * 17 * 23 )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer FACTOR_INDEX, the index of the factor being
!    processed.  The first factor processed must be factor 1!
!
!    Input, integer FACTOR_ORDER, the order of the factor.
!
!    Input, real ( kind = rk ) FACTOR_VALUE(FACTOR_ORDER), the factor values
!    for factor FACTOR_INDEX.
!
!    Input, integer FACTOR_NUM, the number of factors.
!
!    Input, integer POINT_NUM, the number of elements in the
!    direct product.
!
!    Input/output, real ( kind = rk ) W(POINT_NUM), the elements of the
!    direct product, which are built up gradually.
!
!  Local Parameters:
!
!    Local, integer START, the first location of a block of values
!    to set.
!
!    Local, integer CONTIG, the number of consecutive values
!    to set.
!
!    Local, integer SKIP, the distance from the current value of START
!    to the next location of a block of values to set.
!
!    Local, integer REP, the number of blocks of values to set.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer factor_num
  integer factor_order
  integer point_num

  integer, save :: contig
  integer factor_index
  real ( kind = rk ) factor_value(factor_order)
  integer j
  integer k
  integer, save :: rep
  integer, save :: skip
  integer start
  real ( kind = rk ) w(point_num)

  if ( factor_index == 1 ) then
    contig = 1
    skip = 1
    rep = point_num
    w(1:point_num) = 1.0D+00
  end if

  rep = rep / factor_order
  skip = skip * factor_order

  do j = 1, factor_order

    start = 1 + ( j - 1 ) * contig

    do k = 1, rep
      w(start:start+contig-1) = w(start:start+contig-1) * factor_value(j)
      start = start + skip
    end do

  end do

  contig = contig * factor_order

  return
end
subroutine s_blank_delete ( s )

!*****************************************************************************80
!
!! S_BLANK_DELETE removes blanks from a string, left justifying the remainder.
!
!  Discussion:
!
!    All TAB characters are also removed.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
  implicit none

  character ch
  integer get
  integer put
  character ( len = * ) s
  integer s_length
  character, parameter :: tab = achar ( 9 )

  put = 0
  s_length = len_trim ( s )

  do get = 1, s_length

    ch = s(get:get)

    if ( ch /= ' ' .and. ch /= tab ) then
      put = put + 1
      s(put:put) = ch
    end if

  end do

  s(put+1:s_length) = ' '

  return
end
subroutine sort_heap_external ( n, indx, i, j, isgn )

!*****************************************************************************80
!
!! SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
!
!  Discussion:
!
!    The actual list of data is not passed to the routine.  Hence this
!    routine may be used to sort integers, reals, numbers, names,
!    dates, shoe sizes, and so on.  After each call, the routine asks
!    the user to compare or interchange two items, until a special
!    return value signals that the sorting is completed.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 February 2004
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer N, the number of items to be sorted.
!
!    Input/output, integer INDX, the main communication signal.
!
!    The user must set INDX to 0 before the first call.
!    Thereafter, the user should not change the value of INDX until
!    the sorting is done.
!
!    On return, if INDX is
!
!      greater than 0,
!      * interchange items I and J;
!      * call again.
!
!      less than 0,
!      * compare items I and J;
!      * set ISGN = -1 if I < J, ISGN = +1 if J < I;
!      * call again.
!
!      equal to 0, the sorting is done.
!
!    Output, integer I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ISGN, results of comparison of elements
!    I and J. (Used only when the previous call returned INDX less than 0).
!    ISGN <= 0 means I is less than or equal to J;
!    0 <= ISGN means I is greater than or equal to J.
!
  implicit none

  integer i
  integer, save :: i_save = 0
  integer indx
  integer isgn
  integer j
  integer, save :: j_save = 0
  integer, save :: k = 0
  integer, save :: k1 = 0
  integer n
  integer, save :: n1 = 0
!
!  INDX = 0: This is the first call.
!
  if ( indx == 0 ) then

    i_save = 0
    j_save = 0
    k = n / 2
    k1 = k
    n1 = n
!
!  INDX < 0: The user is returning the results of a comparison.
!
  else if ( indx < 0 ) then

    if ( indx == -2 ) then

      if ( isgn < 0 ) then
        i_save = i_save + 1
      end if

      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return

    end if

    if ( 0 < isgn ) then
      indx = 2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then

      if ( n1 == 1 ) then
        i_save = 0
        j_save = 0
        indx = 0
      else
        i_save = n1
        n1 = n1 - 1
        j_save = 1
        indx = 1
      end if

      i = i_save
      j = j_save
      return

    end if

    k = k - 1
    k1 = k
!
!  0 < INDX, the user was asked to make an interchange.
!
  else if ( indx == 1 ) then

    k1 = k

  end if

  do

    i_save = 2 * k1

    if ( i_save == n1 ) then
      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return
    else if ( i_save <= n1 ) then
      j_save = i_save + 1
      indx = -2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then
      exit
    end if

    k = k - 1
    k1 = k

  end do

  if ( n1 == 1 ) then
    i_save = 0
    j_save = 0
    indx = 0
    i = i_save
    j = j_save
  else
    i_save = n1
    n1 = n1 - 1
    j_save = 1
    indx = 1
    i = i_save
    j = j_save
  end if

  return
end
subroutine sparse_grid_hermite ( dim_num, level_max, point_total_num, &
  point_num, sparse_weight, sparse_point )

!****************************************************************************80
!
!! SPARSE_GRID_HERMITE computes a Hermite sparse grid.
!
!  Discussion:
!
!    The quadrature rule is associated with a sparse grid derived from
!    a Smolyak construction using a 1D Gauss-Hermite quadrature rule.
!
!    The user specifies:
!    * the spatial dimension of the quadrature region,
!    * the level that defines the Smolyak grid.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Fabio Nobile, Raul Tempone, Clayton Webster,
!    A Sparse Grid Stochastic Collocation Method for Partial Differential
!    Equations with Random Input Data,
!    SIAM Journal on Numerical Analysis,
!    Volume 46, Number 5, 2008, pages 2309-2345.
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer LEVEL_MAX, controls the size of the
!    sparse grid.
!
!    Input, integer POINT_TOTAL_NUM, the number of points in the 
!    grid, as determined by SPARSE_GRID_HERMITE_SIZE.
!
!    Input, integer POINT_NUM, the number of points in the grid,
!    as determined by SPARSE_GRID_HERMITE_SIZE.
!
!    Output, real ( kind = rk ) SPARSE_WEIGHT(POINT_NUM), the weights.
!
!    Output, real ( kind = rk ) SPARSE_POINT(DIM_NUM,POINT_NUM), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num
  integer point_num

  integer level_max
  integer point_total_num
  integer, allocatable, dimension ( :, : ) :: sparse_index
  integer, allocatable, dimension ( :, : ) :: sparse_order
  real ( kind = rk ) sparse_point(dim_num,point_num)
  integer, allocatable, dimension ( : ) :: sparse_unique_index
  real ( kind = rk ) sparse_weight(point_num)
  real ( kind = rk ) tol

  tol = 1.0D-07

  allocate ( sparse_unique_index(1:point_total_num) )

  call sparse_grid_hermite_unique_index ( dim_num, level_max, tol, &
    point_num, point_total_num, sparse_unique_index )

  allocate ( sparse_index(1:dim_num,1:point_num) )
  allocate ( sparse_order(1:dim_num,1:point_num) )

  call sparse_grid_hermite_index ( dim_num, level_max, point_num, &
    point_total_num, sparse_unique_index, sparse_order, sparse_index )

  call sparse_grid_hermite_point ( dim_num, level_max, point_num, &
    sparse_order, sparse_index, sparse_point )

  call sparse_grid_hermite_weight ( dim_num, level_max, point_num, &
    point_total_num, sparse_unique_index, sparse_weight )

  deallocate ( sparse_index )
  deallocate ( sparse_order )
  deallocate ( sparse_unique_index )

  return
end
subroutine sparse_grid_hermite_index ( dim_num, level_max, point_num, &
  point_total_num, sparse_unique_index, sparse_order, sparse_index )

!*****************************************************************************80
!
!! SPARSE_GRID_HERMITE_INDEX indexes a Hermite sparse grid.
!
!  Discussion:
!
!    For each "unique" point in the sparse grid, we return its INDEX and ORDER.
!
!    That is, for the I-th unique point P, we determine the product grid which
!    first generated this point, and we return in SPARSE_ORDER the orders of
!    the 1D rules in that grid, and in SPARSE_INDEX the component indexes in
!    those rules that generated this specific point.
!
!    For instance, say P was first generated by a rule which was a 3D product
!    of a 9th order CC rule and a 15th order GL rule, and that to generate P,
!    we used the 7-th point of the CC rule and the 3rd point of the GL rule.
!    Then the SPARSE_ORDER information would be (9,15) and the SPARSE_INDEX
!    information would be (7,3).  This, combined with the information in RULE,
!    is enough to regenerate the value of P.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 December 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Fabio Nobile, Raul Tempone, Clayton Webster,
!    A Sparse Grid Stochastic Collocation Method for Partial Differential
!    Equations with Random Input Data,
!    SIAM Journal on Numerical Analysis,
!    Volume 46, Number 5, 2008, pages 2309-2345.
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer LEVEL_MAX, the maximum value of LEVEL.
!
!    Input, integer POINT_NUM, the number of unique points
!    in the grid.
!
!    Input, integer POINT_TOTAL_NUM, the total number of points
!    in the grid.
!
!    Input, integer SPARSE_UNIQUE_INDEX(POINT_TOTAL_NUM),
!    associates each point in the grid with its unique representative.
!
!    Output, integer SPARSE_ORDER(DIM_NUM,POINT_NUM), lists,
!    for each point, the order of the 1D rules used in the grid that
!    generated it.
!
!    Output, integer SPARSE_INDEX(DIM_NUM,POINT_NUM), lists, for
!    each point, its index in each of the 1D rules in the grid that generated
!    it.  The indices are 1-based.
!
  implicit none

  integer dim_num
  integer point_num
  integer point_total_num

  integer h
  integer level
  integer level_1d(dim_num)
  integer level_max
  integer level_min
  logical more_grids
  logical more_points
  integer order_1d(dim_num)
  integer point_count
  integer point_index(dim_num)
  integer point_unique
  integer sparse_index(dim_num,point_num)
  integer sparse_order(dim_num,point_num)
  integer sparse_unique_index(point_total_num)
  integer t
!
!  Special cases.
!
  if ( level_max < 0 ) then
    return
  end if

  if ( level_max == 0 ) then
    sparse_order(1:dim_num,1) = 1
    sparse_index(1:dim_num,1) = 1
    return
  end if
!
!  Initialize to -1 to help catch errors.
!
  sparse_order(1:dim_num,1:point_num) = -1
  sparse_index(1:dim_num,1:point_num) = -1

  point_count = 0
!
!  The outer loop generates values of LEVEL.
!
  level_min = max ( 0, level_max + 1 - dim_num )

  do level = level_min, level_max
!
!  The middle loop generates a GRID,
!  based on the next partition that adds up to LEVEL.
!
    more_grids = .false.
    h = 0
    t = 0

    do

      call comp_next ( level, dim_num, level_1d, more_grids, h, t )

      call level_to_order_hermite ( dim_num, level_1d, order_1d )
!
!  The inner loop generates a POINT of the GRID of the LEVEL.
!
      more_points = .false.

      do

        call vec_colex_next3 ( dim_num, order_1d, point_index, more_points )

        if ( .not. more_points ) then
          exit
        end if

        point_count = point_count + 1
        point_unique = sparse_unique_index(point_count)
        sparse_order(1:dim_num,point_unique) = order_1d(1:dim_num)
        sparse_index(1:dim_num,point_unique) = point_index(1:dim_num)

      end do

      if ( .not. more_grids ) then
        exit
      end if

    end do

  end do

  return
end
subroutine sparse_grid_hermite_point ( dim_num, level_max, point_num, &
  sparse_order, sparse_index, sparse_point )

!*****************************************************************************80
!
!! SPARSE_GRID_HERMITE_POINT computes the points of a Hermite sparse grid.
!
!  Discussion:
!
!    The sparse grid is the logical sum of low degree product rules.
!
!    Each product rule is the product of 1D factor rules.
!
!    The user specifies:
!    * the spatial dimension of the quadrature region,
!    * the level that defines the Smolyak grid.
!    * the quadrature rules.
!    * the number of points.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Fabio Nobile, Raul Tempone, Clayton Webster,
!    A Sparse Grid Stochastic Collocation Method for Partial Differential
!    Equations with Random Input Data,
!    SIAM Journal on Numerical Analysis,
!    Volume 46, Number 5, 2008, pages 2309-2345.
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer LEVEL_MAX, controls the size of the final
!    sparse grid.
!
!    Input, integer POINT_NUM, the number of unique points
!    in the grid.
!
!    Input, integer SPARSE_ORDER(DIM_NUM,POINT_NUM), lists, for
!    each point, the order of the 1D rules used in the grid that generated it.
!
!    Input, integer SPARSE_INDEX(DIM_NUM,POINT_NUM), lists, for
!    each point, its index in each of the 1D rules in the grid that generated
!    it.  The indices are 1-based.
!
!    Output, real ( kind = rk ) SPARSE_POINT(DIM_NUM,POINT_NUM), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num
  integer point_num

  integer dim
  integer level
  integer level_max
  integer order
  integer point
  real ( kind = rk ), allocatable, dimension ( : ) :: points
  real ( kind = rk ) r8_huge
  integer sparse_index(dim_num,point_num)
  integer sparse_order(dim_num,point_num)
  real ( kind = rk ) sparse_point(dim_num,point_num)
!
!  Compute the point coordinates.
!
  sparse_point(1:dim_num,1:point_num) = - r8_huge ( )

  do dim = 1, dim_num

    do level = 0, level_max

      call level_to_order_hermite ( 1, level, order )

      allocate ( points(1:order) )

      call hermite_compute_points ( order, points )

      do point = 1, point_num
        if ( sparse_order(dim,point) == order ) then
          sparse_point(dim,point) = points ( sparse_index(dim,point) )
        end if
      end do

      deallocate ( points )

    end do
  end do

  return
end
subroutine sparse_grid_hermite_size ( dim_num, level_max, point_num, &
  point_total_num )

!*****************************************************************************80
!
!! SPARSE_GRID_HERMITE_SIZE sizes a Hermite sparse grid, discounting duplicates.
!
!  Discussion:
!
!    The sparse grid is the logical sum of product grids with total LEVEL
!    between LEVEL_MIN and LEVEL_MAX.
!
!    Depending on the 1D rules involved, there may be many duplicate points
!    in the sparse grid.
!
!    This routine counts the unique points in the sparse grid.  It does this
!    in a straightforward way, by actually generating all the points, and
!    comparing them, with a tolerance for equality.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Fabio Nobile, Raul Tempone, Clayton Webster,
!    A Sparse Grid Stochastic Collocation Method for Partial Differential
!    Equations with Random Input Data,
!    SIAM Journal on Numerical Analysis,
!    Volume 46, Number 5, 2008, pages 2309-2345.
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer LEVEL_MAX, the maximum value of LEVEL.
!
!    Output, integer POINT_NUM, the number of unique points
!    in the grid.
!
!    Output, integer POINT_TOTAL_NUM, the number of points
!    in the grid, counting duplicates.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num

  integer dim
  integer h
  integer level
  integer level_1d(dim_num)
  integer level_max
  integer level_min
  logical more_grids
  logical more_points
  integer order
  integer order_1d(dim_num)
  integer point
  integer point_index(dim_num)
  integer point_num
  integer point_total_num
  integer point_total_num2
  real ( kind = rk ), allocatable, dimension ( : ) :: points
  real ( kind = rk ) r8_huge
  integer, allocatable, dimension ( :, : ) :: sparse_total_index
  integer, allocatable, dimension ( :, : ) :: sparse_total_order
  real ( kind = rk ), allocatable, dimension ( :, : ) :: sparse_total_point
  integer t
  real ( kind = rk ) tol

  tol = 1.0D-07
!
!  Special cases.
!
  if ( level_max < 0 ) then
    point_num = -1
    return
  end if

  if ( level_max == 0 ) then
    point_num = 1
    return
  end if
!
!  Get total number of points, including duplicates.
!
  call sparse_grid_hermite_size_total ( dim_num, level_max, point_total_num )
!
!  Generate SPARSE_TOTAL_ORDER and SPARSE_TOTAL_INDEX arrays
!  for the TOTAL set of points.
!
  allocate ( sparse_total_order(1:dim_num,1:point_total_num ) )
  allocate ( sparse_total_index(1:dim_num,1:point_total_num ) )

  point_total_num2 = 0
!
!  The outer loop generates values of LEVEL.
!
  level_min = max ( 0, level_max + 1 - dim_num )

  do level = level_min, level_max
!
!  The middle loop generates a GRID,
!  based on the next partition that adds up to LEVEL.
!
    more_grids = .false.
    h = 0
    t = 0

    do

      call comp_next ( level, dim_num, level_1d, more_grids, h, t )

      call level_to_order_hermite ( dim_num, level_1d, order_1d )
!
!  The inner loop generates a POINT of the GRID of the LEVEL.
!
      more_points = .false.

      do

        call vec_colex_next3 ( dim_num, order_1d, point_index, more_points )

        if ( .not. more_points ) then
          exit
        end if

        point_total_num2 = point_total_num2 + 1
        sparse_total_order(1:dim_num,point_total_num2) = order_1d(1:dim_num)
        sparse_total_index(1:dim_num,point_total_num2) = point_index(1:dim_num)

      end do

      if ( .not. more_grids ) then
        exit
      end if

    end do

  end do
!
!  Now compute the coordinates of the TOTAL set of points.
!
  allocate ( sparse_total_point(1:dim_num,1:point_total_num) )
  sparse_total_point(1:dim_num,1:point_total_num) = r8_huge ( )

  do dim = 1, dim_num

    do level = 0, level_max

      call level_to_order_hermite ( 1, level, order )

      allocate ( points(1:order) )

      call hermite_compute_points ( order, points )

      do point = 1, point_total_num
        if ( sparse_total_order(dim,point) == order ) then
          sparse_total_point(dim,point) = &
            points ( sparse_total_index(dim,point) )
        end if
      end do

      deallocate ( points )

    end do

  end do
!
!  Sort the columns.
!
  call r8col_sort_heap_a ( dim_num, point_total_num, sparse_total_point )
!
!  Count the unique columns.
!
  call r8col_sorted_unique_count ( dim_num, point_total_num, &
    sparse_total_point, tol, point_num )

  deallocate ( sparse_total_index )
  deallocate ( sparse_total_order )
  deallocate ( sparse_total_point )

  return
end
subroutine sparse_grid_hermite_size_total ( dim_num, level_max, &
  point_total_num )

!*****************************************************************************80
!
!! SPARSE_GRID_HERMITE_SIZE_TOTAL Hermite sparse grid size counting duplicates.
!
!  Discussion:
!
!    The sparse grid is the logical sum of product grids with total LEVEL
!    between LEVEL_MIN and LEVEL_MAX.
!
!    In some cases, the same point may occur in different product grids
!    used to form the sparse grid.
!
!    This routine counts the total number of points used to construct the sparse
!    grid; if the same point occurs several times, each occurrence is added
!    to the sum.
!
!    This computation is useful in order to be able to allocate enough
!    space for the full set of points, before they are compressed by removing
!    duplicates.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 December 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Fabio Nobile, Raul Tempone, Clayton Webster,
!    A Sparse Grid Stochastic Collocation Method for Partial Differential
!    Equations with Random Input Data,
!    SIAM Journal on Numerical Analysis,
!    Volume 46, Number 5, 2008, pages 2309-2345.
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer LEVEL_MAX, the maximum value of LEVEL.
!
!    Output, integer POINT_TOTAL_NUM, the total number of points
!    in the grid.
!
  implicit none

  integer dim_num

  integer h
  integer level
  integer level_1d(dim_num)
  integer level_max
  integer level_min
  logical more_grids
  integer order_1d(dim_num)
  integer point_total_num
  integer t
!
!  Special case.
!
  if ( level_max == 0 ) then
    point_total_num = 1
    return
  end if

  point_total_num = 0
!
!  The outer loop generates values of LEVEL.
!
  level_min = max ( 0, level_max + 1 - dim_num )

  do level = level_min, level_max
!
!  The middle loop generates a GRID,
!  based on the next partition that adds up to LEVEL.
!
    more_grids = .false.
    h = 0
    t = 0

    do

      call comp_next ( level, dim_num, level_1d, more_grids, h, t )

      call level_to_order_hermite ( dim_num, level_1d, order_1d )

      point_total_num = point_total_num + product ( order_1d(1:dim_num) )

      if ( .not. more_grids ) then
        exit
      end if

    end do

  end do

  return
end
subroutine sparse_grid_hermite_unique_index ( dim_num, level_max, tol, &
  point_num, point_total_num, sparse_unique_index )

!*****************************************************************************80
!
!! SPARSE_GRID_HERMITE_UNIQUE_INDEX maps nonunique points to unique points.
!
!  Discussion:
!
!    The sparse grid usually contains many points that occur in more
!    than one product grid.
!
!    When generating the point locations, it is easy to realize that a point
!    has already been generated.
!
!    But when it's time to compute the weights of the sparse grids, it is
!    necessary to handle situations in which weights corresponding to
!    the same point generated in multiple grids must be collected together.
!
!    This routine generates ALL the points, including their multiplicities,
!    and figures out a mapping from them to the collapsed set of unique points.
!
!    This mapping can then be used during the weight calculation so that
!    a contribution to the weight gets to the right place.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 July 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Fabio Nobile, Raul Tempone, Clayton Webster,
!    A Sparse Grid Stochastic Collocation Method for Partial Differential
!    Equations with Random Input Data,
!    SIAM Journal on Numerical Analysis,
!    Volume 46, Number 5, 2008, pages 2309-2345.
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer LEVEL_MAX, the maximum value of LEVEL.
!
!    Input, real ( kind = rk ) TOL, the tolerance for point equality.
!
!    Input, integer POINT_NUM, the number of unique points in
!    the grid.
!
!    Input, integer POINT_TOTAL_NUM, the total number of points
!    in the grid.
!
!    Output, integer SPARSE_UNIQUE_INDEX(POINT_TOTAL_NUM), lists,
!    for each (nonunique) point, the corresponding index of the same point in
!    the unique listing.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num
  integer point_total_num

  integer dim
  integer h
  integer level
  integer level_1d(dim_num)
  integer level_max
  integer level_min
  logical more_grids
  logical more_points
  integer order
  integer order_1d(dim_num)
  integer point
  integer point_index(dim_num)
  integer point_num
  integer point_total_num2
  real ( kind = rk ), allocatable, dimension ( : ) :: points
  real ( kind = rk ) r8_huge
  integer, allocatable, dimension ( :, : ) :: sparse_total_index
  integer, allocatable, dimension ( :, : ) :: sparse_total_order
  real ( kind = rk ), allocatable, dimension ( :, : ) :: sparse_total_point
  integer sparse_unique_index(point_total_num)
  integer t
  real ( kind = rk ) tol
  integer, allocatable, dimension ( : ) :: undx
!
!  Special cases.
!
  if ( level_max < 0 ) then
    return
  end if

  if ( level_max == 0 ) then
    sparse_unique_index(1) = 1
    return
  end if
!
!  Generate SPARSE_TOTAL_ORDER and SPARSE_TOTAL_INDEX arrays
!  for the TOTAL set of points.
!
  allocate ( sparse_total_order(1:dim_num,1:point_total_num ) )
  allocate ( sparse_total_index(1:dim_num,1:point_total_num ) )

  point_total_num2 = 0
!
!  The outer loop generates values of LEVEL.
!
  level_min = max ( 0, level_max + 1 - dim_num )

  do level = level_min, level_max
!
!  The middle loop generates a GRID,
!  based on the next partition that adds up to LEVEL.
!
    more_grids = .false.
    h = 0
    t = 0

    do

      call comp_next ( level, dim_num, level_1d, more_grids, h, t )

      call level_to_order_hermite ( dim_num, level_1d, order_1d )
!
!  The inner loop generates a POINT of the GRID of the LEVEL.
!
      more_points = .false.

      do

        call vec_colex_next3 ( dim_num, order_1d, point_index, more_points )

        if ( .not. more_points ) then
          exit
        end if

        point_total_num2 = point_total_num2 + 1
        sparse_total_order(1:dim_num,point_total_num2) = order_1d(1:dim_num)
        sparse_total_index(1:dim_num,point_total_num2) = point_index(1:dim_num)

      end do

      if ( .not. more_grids ) then
        exit
      end if

    end do

  end do
!
!  Now compute the coordinates of the TOTAL set of points.
!
  allocate ( sparse_total_point(1:dim_num,1:point_total_num) )
  sparse_total_point(1:dim_num,1:point_total_num) = r8_huge ( )

  do dim = 1, dim_num

    do level = 0, level_max

      call level_to_order_hermite ( 1, level, order )

      allocate ( points(1:order) )

      call hermite_compute_points ( order, points )

      do point = 1, point_total_num
        if ( sparse_total_order(dim,point) == order ) then
          sparse_total_point(dim,point) = &
            points ( sparse_total_index(dim,point) )
        end if
      end do

      deallocate ( points )

    end do

  end do
!
!  Now determine the mapping from nonunique points to unique points.
!  We can't really use the UNDX output right now.
!
  allocate ( undx(1:point_num) )

  call r8col_undex ( dim_num, point_total_num, sparse_total_point, point_num, &
    tol, undx, sparse_unique_index )

! call r8col_tol_undex ( dim_num, point_total_num, sparse_total_point, &
!   point_num, tol, undx, sparse_unique_index )

  deallocate ( sparse_total_index )
  deallocate ( sparse_total_order )
  deallocate ( sparse_total_point )
  deallocate ( undx )

  return
end
subroutine sparse_grid_hermite_weight ( dim_num, level_max, point_num, &
  point_total_num, sparse_unique_index, sparse_weight )

!*****************************************************************************80
!
!! SPARSE_GRID_HERMITE_WEIGHT computes Hermite sparse grid weights.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 December 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Fabio Nobile, Raul Tempone, Clayton Webster,
!    A Sparse Grid Stochastic Collocation Method for Partial Differential
!    Equations with Random Input Data,
!    SIAM Journal on Numerical Analysis,
!    Volume 46, Number 5, 2008, pages 2309-2345.
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer LEVEL_MAX, the maximum value of LEVEL.
!
!    Input, integer POINT_NUM, the number of unique points in
!    the grid.
!
!    Input, integer POINT_TOTAL_NUM, the total number of points
!    in the grid.
!
!    Input, integer SPARSE_UNIQUE_INDEX(POINT_TOTAL_NUM), lists,
!    for each (nonunique) point, the corresponding index of the same point in
!    the unique listing.
!
!    Output, real ( kind = rk ) SPARSE_WEIGHT(POINT_NUM), the weights
!    associated with the sparse grid points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num
  integer point_num
  integer point_total_num

  real ( kind = rk ) coeff
  real ( kind = rk ), allocatable, dimension ( : ) :: grid_weight
  integer h
  integer level
  integer level_1d(dim_num)
  integer level_max
  integer level_min
  logical more_grids
  integer order
  integer order_1d(dim_num)
  integer order_nd
  integer point_total
  integer point_unique
  real ( kind = rk ) r8_choose
  real ( kind = rk ) r8_mop
  integer sparse_unique_index(point_total_num)
  real ( kind = rk ) sparse_weight(point_num)
  integer t

  sparse_weight(1:point_num) = 0.0D+00

  point_total = 0

  level_min = max ( 0, level_max + 1 - dim_num )

  do level = level_min, level_max
!
!  The middle loop generates the next partition LEVEL_1D(1:DIM_NUM)
!  that adds up to LEVEL.
!
    more_grids = .false.
    h = 0
    t = 0

    do

      call comp_next ( level, dim_num, level_1d, more_grids, h, t )
!
!  Transform each 1D level to a corresponding 1D order.
!
      call level_to_order_hermite ( dim_num, level_1d, order_1d )
!
!  The product of the 1D orders gives us the number of points in this grid.
!
      order_nd = product ( order_1d(1:dim_num) )
!
!  Compute the weights for this grid.
!
!  The correct transfer of data from the product grid to the sparse grid
!  depends on the fact that the product rule weights are stored under colex
!  order of the points, and this is the same ordering implicitly used in
!  generating the SPARSE_UNIQUE_INDEX array.
!
      allocate ( grid_weight(1:order_nd) )

      call product_hermite_weight ( dim_num, order_1d, order_nd, grid_weight )
!
!  Compute Smolyak's binomial coefficient for this grid.
!
      coeff = r8_mop ( level_max - level ) &
        * r8_choose ( dim_num - 1, level_max - level )
!
!  Add these weights to the rule.
!
      do order = 1, order_nd

        point_total = point_total + 1

        point_unique = sparse_unique_index(point_total)

        sparse_weight(point_unique) = sparse_weight(point_unique) &
          + coeff * grid_weight(order)

      end do

      deallocate ( grid_weight )

      if ( .not. more_grids ) then
        exit
      end if

    end do

  end do

  return
end
subroutine sparse_grid_hermite_write ( dim_num, point_num, sparse_weight, &
  sparse_point, filename )

!*****************************************************************************80
!
!! SPARSE_GRID_HERMITE_WRITE writes a Hermite sparse grid rule to files.
!
!  Discussion:
!
!    The files are:
!    * the "R" file stores the region, as a DIM_NUM x 2 list.
!    * the "W" file stores the weights as a POINT_NUM list;
!    * the "X" file stores the abscissas as a DIM_NUM x POINT_NUM list;
!
!    The entries in the "R" file are the two corners of the DIM_NUM dimensional
!    rectangle that constitutes the integration region.  Coordinates that
!    should be infinite are set to 1.0E+30.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Fabio Nobile, Raul Tempone, Clayton Webster,
!    A Sparse Grid Stochastic Collocation Method for Partial Differential
!    Equations with Random Input Data,
!    SIAM Journal on Numerical Analysis,
!    Volume 46, Number 5, 2008, pages 2309-2345.
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer POINT_NUM, the number of unique points
!    in the grid.
!
!    Input, real ( kind = rk ) SPARSE_WEIGHT(POINT_NUM), the weights.
!
!    Input, real ( kind = rk ) SPARSE_POINT(DIM_NUM,POINT_NUM), the points.
!
!    Input, character ( len = * ) FILENAME, the main part of the file name.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num
  integer point_num

  integer dim
  character ( len = *  ) filename
  character ( len = 80 ) filename_r
  character ( len = 80 ) filename_w
  character ( len = 80 ) filename_x
  real ( kind = rk ) r8_huge
  real ( kind = rk ) sparse_point(dim_num,point_num)
  real ( kind = rk ) sparse_region(dim_num,2)
  real ( kind = rk ) sparse_weight(point_num)

  do dim = 1, dim_num
    sparse_region(dim,1) = - r8_huge ( )
    sparse_region(dim,2) = + r8_huge ( )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSE_GRID_HERMITE_WRITE:'

  filename_r = trim ( filename ) // '_r.txt'
  call r8mat_write ( filename_r, dim_num, 2, sparse_region )
  write ( *, '(a)' ) '  Wrote the R file = "' // trim ( filename_r ) // '".'

  filename_w = trim ( filename ) // '_w.txt'
  call r8mat_write ( filename_w, 1, point_num, sparse_weight )
  write ( *, '(a)' ) '  Wrote the W file = "' // trim ( filename_w ) // '".'

  filename_x = trim ( filename ) // '_x.txt'
  call r8mat_write ( filename_x, dim_num, point_num, sparse_point )
  write ( *, '(a)' ) '  Wrote the X file = "' // trim ( filename_x ) // '".'

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
subroutine vec_colex_next3 ( dim_num, base, a, more )

!*****************************************************************************80
!
!! VEC_COLEX_NEXT3 generates vectors in colex order.
!
!  Discussion:
!
!    The vectors are produced in colexical order, starting with
!
!    (1,        1,        ...,1),
!    (2,        1,        ...,1),
!     ...
!    (BASE(1),  1,        ...,1)
!
!    (1,        2,        ...,1)
!    (2,        2,        ...,1)
!    ...
!    (BASE(1),  2,        ...,0)
!
!    (1,        3,        ...,1)
!    (2,        3,        ...,1)
!    ...
!    (BASE(1),  BASE(2),  ...,BASE(DIM_NUM)).
!
!  Example:
!
!    DIM_NUM = 2,
!    BASE = ( 3, 3 )
!
!    1   1
!    2   1
!    3   1
!    1   2
!    2   2
!    3   2
!    1   3
!    2   3
!    3   3
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 August 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer BASE(DIM_NUM), the bases to be used in each
!    dimension.  In dimension I, entries will range from 1 to BASE(I).
!
!    Input/output, integer A(DIM_NUM).
!
!    Input/output, logical MORE.  Set this variable FALSE before
!    the first call.  On return, MORE is TRUE if another vector has
!    been computed.  If MORE is returned FALSE, ignore the output
!    vector and stop calling the routine.
!
  implicit none

  integer dim_num

  integer a(dim_num)
  integer base(dim_num)
  integer i
  logical more

  if ( .not. more ) then

    a(1:dim_num) = 1
    more = .true.

  else

    do i = 1, dim_num

      a(i) = a(i) + 1

      if ( a(i) <= base(i) ) then
        return
      end if

      a(i) = 1

    end do

    more = .false.

  end if

  return
end
