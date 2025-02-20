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
!    Input/output, integer H, T, two internal parameters needed
!    for the computation.  The user should allocate space for these in the
!    calling program, include them in the calling sequence, but never alter
!    them!
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

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
subroutine monomial_value ( d, n, e, x, v )

!*****************************************************************************80
!
!! monomial_value() evaluates a monomial.
!
!  Discussion:
!
!    This routine evaluates a monomial of the form
!
!      product ( 1 <= i <= d ) x(i)^e(i)
!
!    The combination 0.0^0, if encountered, is treated as 1.0.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 May 2023
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer D, the spatial dimension.
!
!    integer N, the number of evaluation points.
!
!    integer E(D), the exponents.
!
!    real ( kind = rk ) X(N,D), the point coordinates.
!
!  Output:
!
!    real ( kind = rk ) V(N), the monomial values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer d
  integer n

  integer e(d)
  integer j
  real ( kind = rk ) v(n)
  real ( kind = rk ) x(n,d)

  v(1:n) = 1.0D+00

  do j = 1, d
    if ( 0 /= e(j) ) then
      v(1:n) = v(1:n) * x(1:n,j) ** e(j)
    end if
  end do

  return
end
subroutine prism_jaskowiec_rule ( p, n, x, y, z, w )

!*****************************************************************************80
!
!! prism_jaskowiec_rule() returns a prism quadrature rule of given precision.
!
!  Discussion:
!
!    A unit triangular prism has a triangular base in the (x,y) plane,
!    projected vertically one unit in the z direction.  
!
!    The vertices are 
!      (1,0,0), (0,1,0), (0,0,0), 
!      (1,0,1), (0,1,1), (0,0,1), 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 May 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order symmetric cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 122, Number 1, pages 148-171, 24 August 2020.
!
!  Input:
!
!    int p: the precision, 0 <= p <= 20.
!
!    int n: the order of the rule.
!
!  Output:
!
!    double x(n), y(n), z(n): the coordinates of quadrature points.
!
!    double w(n): the quadrature weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer p
  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)
  real ( kind = rk ) z(n)

  if ( p < 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'prism_jaskowiec_rule(): Fatal error!'
    write ( *, '(a)' ) '  Input p < 0.'
    stop ( 1 )
  end if

  if ( 20 < p ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'prism_jaskowiec_rule(): Fatal error!'
    write ( *, '(a)' ) '  Input 20 < p.'
    stop ( 1 )
  end if

  if ( p == 0 ) then
    call rule00 ( n, x, y, z, w )
  else if ( p == 1 ) then
    call rule01 ( n, x, y, z, w )
  else if ( p == 2 ) then
    call rule02 ( n, x, y, z, w )
  else if ( p == 3 ) then
    call rule03 ( n, x, y, z, w )
  else if ( p == 4 ) then
    call rule04 ( n, x, y, z, w )
  else if ( p == 5 ) then
    call rule05 ( n, x, y, z, w )
  else if ( p == 6 ) then
    call rule06 ( n, x, y, z, w )
  else if ( p == 7 ) then
    call rule07 ( n, x, y, z, w )
  else if ( p == 8 ) then
    call rule08 ( n, x, y, z, w )
  else if ( p == 9 ) then
    call rule09 ( n, x, y, z, w )
  else if ( p == 10 ) then
    call rule10 ( n, x, y, z, w )
  else if ( p == 11 ) then
    call rule11 ( n, x, y, z, w )
  else if ( p == 12 ) then
    call rule12 ( n, x, y, z, w )
  else if ( p == 13 ) then
    call rule13 ( n, x, y, z, w )
  else if ( p == 14 ) then
    call rule14 ( n, x, y, z, w )
  else if ( p == 15 ) then
    call rule15 ( n, x, y, z, w )
  else if ( p == 16 ) then
    call rule16 ( n, x, y, z, w )
  else if ( p == 17 ) then
    call rule17 ( n, x, y, z, w )
  else if ( p == 18 ) then
    call rule18 ( n, x, y, z, w )
  else if ( p == 19 ) then
    call rule19 ( n, x, y, z, w )
  else if ( p == 20 ) then
    call rule20 ( n, x, y, z, w )
  end if

  return
end
subroutine prism_unit_monomial_integral ( expon, value )

!*****************************************************************************80
!
!! prism_unit_monomial_integral(): monomial integral in a unit prism.
!
!  Discussion:
!
!    This routine returns the integral of
!
!      product ( 1 <= I <= 3 ) X(I)^EXPON(I)
!
!    over the unit triangular prism.
!
!    A unit triangular prism has a triangular base in the (x,y) plane,
!    projected vertically one unit in the z direction.  
!
!    The vertices are 
!      (1,0,0), (0,1,0), (0,0,0), 
!      (1,0,1), (0,1,1), (0,0,1), 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 May 2023
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer EXPON(3), the exponents.
!
!  Output:
!
!    real ( kind = rk ) VALUE, the integral of the monomial.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer expon(3)
  real ( kind = rk ) r8_factorial
  real ( kind = rk ) value

  value = r8_factorial ( expon(1) ) &
        * r8_factorial ( expon(2) ) &
        / ( expon(3) + 1 ) &
        / r8_factorial ( expon(1) + expon(2) + 2 )

  return
end
function prism_unit_volume ( )

!*****************************************************************************80
!
!! prism_unit_volume() returns the volume of a unit prism.
!
!  Discussion:
!
!    A unit triangular prism has a triangular base in the (x,y) plane,
!    projected vertically one unit in the z direction.  
!
!    The vertices are 
!      (1,0,0), (0,1,0), (0,0,0), 
!      (1,0,1), (0,1,1), (0,0,1), 
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 May 2023
!
!  Author:
!
!    John Burkardt
!
!  Output:
!
!    real value: the volume.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) prism_unit_volume

  prism_unit_volume = 0.5D+00

  return
end
function r8_factorial ( n )

!*****************************************************************************80
!
!! r8_factorial() computes the factorial of N.
!
!  Discussion:
!
!    factorial ( N ) = product ( 1 <= I <= N ) I
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N: the argument of the factorial function.
!    If N is less than 1, the function value is returned as 1.
!
!  Output:
!
!    real ( kind = rk ) R8_FACTORIAL: the factorial of N.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r8_factorial
  integer i
  integer n
  real ( kind = rk ) value

  value = 1.0D+00

  do i = 1, n
    value = value * real ( i, kind = rk )
  end do

  r8_factorial = value

  return
end
subroutine rule_order ( p, order )

!*****************************************************************************80
!
!! rule_order() returns the order of a prism quadrature rule of given precision.
!
!  Discussion:
!
!    The "order" of a quadrature rule is the number of points involved.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 April 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order symmetric cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 122, Number 1, pages 148-171, 24 August 2020.
!
!  Input:
!
!    integer p: the precision, 0 <= p <= 20.
!
!  Output:
!
!    integer order: the order of the rule.
!
  implicit none

  integer order
  integer, dimension ( 0 : 20 ) :: order_vec = (/ &
      1, &
      1,   5,   8,  11,  16,  28,  35,  46,  59,  84, &
     99, 136, 162, 194, 238, 287, 338, 396, 420, 518 /)
  integer p

  if ( p < 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'rule_order(): Fatal error.'
    write ( *, '(a)' ) '  Input p < 0.'
  end if

  if ( 20 < p ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'rule_order(): Fatal error.'
    write ( *, '(a)' ) '  Input 20 < p.'
  end if

  order = order_vec(p)

  return
end
subroutine rule00 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule00() returns the prism rule of precision 0.
!
!  Discussion:
!
!    We suppose we are given a triangular prism P with vertices
!    (1,0,0), (0,1,0), (0,0,0),
!    (1,0,1), (0,1,1), (0,0,1),
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over P is 
!    approximated by Q as follows:
!
!    Q = volume(P) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 May 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order symmetric cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 122, Number 1, pages 148-171, 24 August 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real x(n), y(n), z(n): the coordinates of quadrature points.
!
!    real w(n): the quadrature weights, which sum to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: n_save = 1

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
          0.3333333333333333D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.3333333333333333D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.5000000000000000D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          1.0000000000000000D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule01 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule01() returns the prism rule of precision 1.
!
!  Discussion:
!
!    We suppose we are given a triangular prism P with vertices
!    (1,0,0), (0,1,0), (0,0,0),
!    (1,0,1), (0,1,1), (0,0,1),
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over P is 
!    approximated by Q as follows:
!
!    Q = volume(P) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 May 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order symmetric cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 122, Number 1, pages 148-171, 24 August 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real x(n), y(n), z(n): the coordinates of quadrature points.
!
!    real w(n): the quadrature weights, which sum to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: n_save = 1

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
          0.3333333333333333D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.3333333333333333D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.5000000000000000D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          1.0000000000000000D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule02 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule02() returns the prism rule of precision 2.
!
!  Discussion:
!
!    We suppose we are given a triangular prism P with vertices
!    (1,0,0), (0,1,0), (0,0,0),
!    (1,0,1), (0,1,1), (0,0,1),
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over P is 
!    approximated by Q as follows:
!
!    Q = volume(P) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 May 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order symmetric cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 122, Number 1, pages 148-171, 24 August 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real x(n), y(n), z(n): the coordinates of quadrature points.
!
!    real w(n): the quadrature weights, which sum to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: n_save = 5

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.1046424703769979D+00, &
          0.1046424703769979D+00, &
          0.7907150592460042D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.1046424703769979D+00, &
          0.7907150592460042D+00, &
          0.1046424703769979D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.0784174655667168D+00, &
          0.9215825344332832D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          0.2344355869392759D+00, &
          0.2344355869392759D+00, &
          0.1770429420404827D+00, &
          0.1770429420404827D+00, &
          0.1770429420404827D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule03 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule03() returns the prism rule of precision 3.
!
!  Discussion:
!
!    We suppose we are given a triangular prism P with vertices
!    (1,0,0), (0,1,0), (0,0,0),
!    (1,0,1), (0,1,1), (0,0,1),
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over P is 
!    approximated by Q as follows:
!
!    Q = volume(P) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 May 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order symmetric cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 122, Number 1, pages 148-171, 24 August 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real x(n), y(n), z(n): the coordinates of quadrature points.
!
!    real w(n): the quadrature weights, which sum to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: n_save = 8

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.4890588576053607D+00, &
          0.4890588576053607D+00, &
          0.0218822847892786D+00, &
          0.0778317780273062D+00, &
          0.0778317780273062D+00, &
          0.8443364439453875D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.4890588576053607D+00, &
          0.0218822847892786D+00, &
          0.4890588576053607D+00, &
          0.0778317780273062D+00, &
          0.8443364439453875D+00, &
          0.0778317780273062D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.0192797910055989D+00, &
          0.9807202089944012D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          0.1803034341765672D+00, &
          0.1803034341765672D+00, &
          0.1134313729984015D+00, &
          0.1134313729984015D+00, &
          0.1134313729984015D+00, &
          0.0996996708838870D+00, &
          0.0996996708838870D+00, &
          0.0996996708838870D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule04 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule04() returns the prism rule of precision 4.
!
!  Discussion:
!
!    We suppose we are given a triangular prism P with vertices
!    (1,0,0), (0,1,0), (0,0,0),
!    (1,0,1), (0,1,1), (0,0,1),
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over P is 
!    approximated by Q as follows:
!
!    Q = volume(P) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 May 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order symmetric cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 122, Number 1, pages 148-171, 24 August 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real x(n), y(n), z(n): the coordinates of quadrature points.
!
!    real w(n): the quadrature weights, which sum to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: n_save = 11

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.4686558098619952D+00, &
          0.4686558098619952D+00, &
          0.0626883802760096D+00, &
          0.1007404057989106D+00, &
          0.1007404057989106D+00, &
          0.7985191884021787D+00, &
          0.1007404057989106D+00, &
          0.1007404057989106D+00, &
          0.7985191884021787D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.4686558098619952D+00, &
          0.0626883802760096D+00, &
          0.4686558098619952D+00, &
          0.1007404057989106D+00, &
          0.7985191884021787D+00, &
          0.1007404057989106D+00, &
          0.1007404057989106D+00, &
          0.7985191884021787D+00, &
          0.1007404057989106D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.0665690129954826D+00, &
          0.9334309870045174D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.8378199118411298D+00, &
          0.8378199118411298D+00, &
          0.8378199118411298D+00, &
          0.1621800881588701D+00, &
          0.1621800881588701D+00, &
          0.1621800881588701D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          0.1079119748155355D+00, &
          0.1079119748155355D+00, &
          0.1364146126054776D+00, &
          0.1364146126054776D+00, &
          0.1364146126054776D+00, &
          0.0624887020920827D+00, &
          0.0624887020920827D+00, &
          0.0624887020920827D+00, &
          0.0624887020920827D+00, &
          0.0624887020920827D+00, &
          0.0624887020920827D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule05 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule05() returns the prism rule of precision 5.
!
!  Discussion:
!
!    We suppose we are given a triangular prism P with vertices
!    (1,0,0), (0,1,0), (0,0,0),
!    (1,0,1), (0,1,1), (0,0,1),
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over P is 
!    approximated by Q as follows:
!
!    Q = volume(P) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 May 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order symmetric cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 122, Number 1, pages 148-171, 24 August 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real x(n), y(n), z(n): the coordinates of quadrature points.
!
!    real w(n): the quadrature weights, which sum to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: n_save = 16

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
          0.3333333333333333D+00, &
          0.0517646178271648D+00, &
          0.0517646178271648D+00, &
          0.8964707643456705D+00, &
          0.1663967696311171D+00, &
          0.1663967696311171D+00, &
          0.6672064607377658D+00, &
          0.1663967696311171D+00, &
          0.1663967696311171D+00, &
          0.6672064607377658D+00, &
          0.4976649895838920D+00, &
          0.4976649895838920D+00, &
          0.0046700208322159D+00, &
          0.4976649895838920D+00, &
          0.4976649895838920D+00, &
          0.0046700208322159D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.3333333333333333D+00, &
          0.0517646178271648D+00, &
          0.8964707643456705D+00, &
          0.0517646178271648D+00, &
          0.1663967696311171D+00, &
          0.6672064607377658D+00, &
          0.1663967696311171D+00, &
          0.1663967696311171D+00, &
          0.6672064607377658D+00, &
          0.1663967696311171D+00, &
          0.4976649895838920D+00, &
          0.0046700208322159D+00, &
          0.4976649895838920D+00, &
          0.4976649895838920D+00, &
          0.0046700208322159D+00, &
          0.4976649895838920D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9035817431942219D+00, &
          0.9035817431942219D+00, &
          0.9035817431942219D+00, &
          0.0964182568057780D+00, &
          0.0964182568057780D+00, &
          0.0964182568057780D+00, &
          0.3013691627751696D+00, &
          0.3013691627751696D+00, &
          0.3013691627751696D+00, &
          0.6986308372248304D+00, &
          0.6986308372248304D+00, &
          0.6986308372248304D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          0.2071428343483058D+00, &
          0.0380755890309976D+00, &
          0.0380755890309976D+00, &
          0.0380755890309976D+00, &
          0.0763742613922619D+00, &
          0.0763742613922619D+00, &
          0.0763742613922619D+00, &
          0.0763742613922619D+00, &
          0.0763742613922619D+00, &
          0.0763742613922619D+00, &
          0.0367308050341883D+00, &
          0.0367308050341883D+00, &
          0.0367308050341883D+00, &
          0.0367308050341883D+00, &
          0.0367308050341883D+00, &
          0.0367308050341883D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule06 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule06() returns the prism rule of precision 6.
!
!  Discussion:
!
!    We suppose we are given a triangular prism P with vertices
!    (1,0,0), (0,1,0), (0,0,0),
!    (1,0,1), (0,1,1), (0,0,1),
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over P is 
!    approximated by Q as follows:
!
!    Q = volume(P) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 May 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order symmetric cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 122, Number 1, pages 148-171, 24 August 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real x(n), y(n), z(n): the coordinates of quadrature points.
!
!    real w(n): the quadrature weights, which sum to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: n_save = 28

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.0175042446586512D+00, &
          0.0175042446586512D+00, &
          0.9649915106826975D+00, &
          0.1617417813899514D+00, &
          0.1617417813899514D+00, &
          0.6765164372200974D+00, &
          0.4656535513495914D+00, &
          0.4656535513495914D+00, &
          0.0686928973008173D+00, &
          0.4656535513495914D+00, &
          0.4656535513495914D+00, &
          0.0686928973008173D+00, &
          0.0345948698524570D+00, &
          0.2025039451729335D+00, &
          0.0345948698524570D+00, &
          0.7629011849746096D+00, &
          0.2025039451729335D+00, &
          0.7629011849746096D+00, &
          0.0345948698524570D+00, &
          0.2025039451729335D+00, &
          0.0345948698524570D+00, &
          0.7629011849746096D+00, &
          0.2025039451729335D+00, &
          0.7629011849746096D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.0175042446586512D+00, &
          0.9649915106826975D+00, &
          0.0175042446586512D+00, &
          0.1617417813899514D+00, &
          0.6765164372200974D+00, &
          0.1617417813899514D+00, &
          0.4656535513495914D+00, &
          0.0686928973008173D+00, &
          0.4656535513495914D+00, &
          0.4656535513495914D+00, &
          0.0686928973008173D+00, &
          0.4656535513495914D+00, &
          0.2025039451729335D+00, &
          0.0345948698524570D+00, &
          0.7629011849746096D+00, &
          0.0345948698524570D+00, &
          0.7629011849746096D+00, &
          0.2025039451729335D+00, &
          0.2025039451729335D+00, &
          0.0345948698524570D+00, &
          0.7629011849746096D+00, &
          0.0345948698524570D+00, &
          0.7629011849746096D+00, &
          0.2025039451729335D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.9925192962207799D+00, &
          0.0074807037792201D+00, &
          0.2480010523941030D+00, &
          0.7519989476058970D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.7405713004233492D+00, &
          0.7405713004233492D+00, &
          0.7405713004233492D+00, &
          0.2594286995766508D+00, &
          0.2594286995766508D+00, &
          0.2594286995766508D+00, &
          0.0952682547954233D+00, &
          0.0952682547954233D+00, &
          0.0952682547954233D+00, &
          0.0952682547954233D+00, &
          0.0952682547954233D+00, &
          0.0952682547954233D+00, &
          0.9047317452045767D+00, &
          0.9047317452045767D+00, &
          0.9047317452045767D+00, &
          0.9047317452045767D+00, &
          0.9047317452045767D+00, &
          0.9047317452045767D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          0.0361446293820950D+00, &
          0.0361446293820950D+00, &
          0.0554469020242208D+00, &
          0.0554469020242208D+00, &
          0.0116429635765844D+00, &
          0.0116429635765844D+00, &
          0.0116429635765844D+00, &
          0.0768806419265571D+00, &
          0.0768806419265571D+00, &
          0.0768806419265571D+00, &
          0.0496268551496232D+00, &
          0.0496268551496232D+00, &
          0.0496268551496232D+00, &
          0.0496268551496232D+00, &
          0.0496268551496232D+00, &
          0.0496268551496232D+00, &
          0.0211237491483504D+00, &
          0.0211237491483504D+00, &
          0.0211237491483504D+00, &
          0.0211237491483504D+00, &
          0.0211237491483504D+00, &
          0.0211237491483504D+00, &
          0.0211237491483504D+00, &
          0.0211237491483504D+00, &
          0.0211237491483504D+00, &
          0.0211237491483504D+00, &
          0.0211237491483504D+00, &
          0.0211237491483504D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule07 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule07() returns the prism rule of precision 7.
!
!  Discussion:
!
!    We suppose we are given a triangular prism P with vertices
!    (1,0,0), (0,1,0), (0,0,0),
!    (1,0,1), (0,1,1), (0,0,1),
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over P is 
!    approximated by Q as follows:
!
!    Q = volume(P) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 May 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order symmetric cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 122, Number 1, pages 148-171, 24 August 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real x(n), y(n), z(n): the coordinates of quadrature points.
!
!    real w(n): the quadrature weights, which sum to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: n_save = 35

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.0083375954660215D+00, &
          0.0083375954660215D+00, &
          0.9833248090679571D+00, &
          0.4815219753291366D+00, &
          0.4815219753291366D+00, &
          0.0369560493417268D+00, &
          0.4815219753291366D+00, &
          0.4815219753291366D+00, &
          0.0369560493417268D+00, &
          0.0954832483714894D+00, &
          0.0954832483714894D+00, &
          0.8090335032570213D+00, &
          0.0954832483714894D+00, &
          0.0954832483714894D+00, &
          0.8090335032570213D+00, &
          0.7429966820728956D+00, &
          0.0121491315983783D+00, &
          0.7429966820728956D+00, &
          0.2448541863287261D+00, &
          0.0121491315983783D+00, &
          0.2448541863287261D+00, &
          0.1529845984247976D+00, &
          0.3051562164322261D+00, &
          0.1529845984247976D+00, &
          0.5418591851429763D+00, &
          0.3051562164322261D+00, &
          0.5418591851429763D+00, &
          0.1529845984247976D+00, &
          0.3051562164322261D+00, &
          0.1529845984247976D+00, &
          0.5418591851429763D+00, &
          0.3051562164322261D+00, &
          0.5418591851429763D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.0083375954660215D+00, &
          0.9833248090679571D+00, &
          0.0083375954660215D+00, &
          0.4815219753291366D+00, &
          0.0369560493417268D+00, &
          0.4815219753291366D+00, &
          0.4815219753291366D+00, &
          0.0369560493417268D+00, &
          0.4815219753291366D+00, &
          0.0954832483714894D+00, &
          0.8090335032570213D+00, &
          0.0954832483714894D+00, &
          0.0954832483714894D+00, &
          0.8090335032570213D+00, &
          0.0954832483714894D+00, &
          0.0121491315983783D+00, &
          0.7429966820728956D+00, &
          0.2448541863287261D+00, &
          0.7429966820728956D+00, &
          0.2448541863287261D+00, &
          0.0121491315983783D+00, &
          0.3051562164322261D+00, &
          0.1529845984247976D+00, &
          0.5418591851429763D+00, &
          0.1529845984247976D+00, &
          0.5418591851429763D+00, &
          0.3051562164322261D+00, &
          0.3051562164322261D+00, &
          0.1529845984247976D+00, &
          0.5418591851429763D+00, &
          0.1529845984247976D+00, &
          0.5418591851429763D+00, &
          0.3051562164322261D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.9901140347954458D+00, &
          0.0098859652045542D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.0793063228967370D+00, &
          0.0793063228967370D+00, &
          0.0793063228967370D+00, &
          0.9206936771032630D+00, &
          0.9206936771032630D+00, &
          0.9206936771032630D+00, &
          0.1020754547065085D+00, &
          0.1020754547065085D+00, &
          0.1020754547065085D+00, &
          0.8979245452934915D+00, &
          0.8979245452934915D+00, &
          0.8979245452934915D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.7019672802660549D+00, &
          0.7019672802660549D+00, &
          0.7019672802660549D+00, &
          0.7019672802660549D+00, &
          0.7019672802660549D+00, &
          0.7019672802660549D+00, &
          0.2980327197339451D+00, &
          0.2980327197339451D+00, &
          0.2980327197339451D+00, &
          0.2980327197339451D+00, &
          0.2980327197339451D+00, &
          0.2980327197339451D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          0.0284456811126893D+00, &
          0.0284456811126893D+00, &
          0.0061182421730950D+00, &
          0.0061182421730950D+00, &
          0.0061182421730950D+00, &
          0.0215508640862265D+00, &
          0.0215508640862265D+00, &
          0.0215508640862265D+00, &
          0.0215508640862265D+00, &
          0.0215508640862265D+00, &
          0.0215508640862265D+00, &
          0.0291785249020985D+00, &
          0.0291785249020985D+00, &
          0.0291785249020985D+00, &
          0.0291785249020985D+00, &
          0.0291785249020985D+00, &
          0.0291785249020985D+00, &
          0.0255148563351493D+00, &
          0.0255148563351493D+00, &
          0.0255148563351493D+00, &
          0.0255148563351493D+00, &
          0.0255148563351493D+00, &
          0.0255148563351493D+00, &
          0.0389407032762076D+00, &
          0.0389407032762076D+00, &
          0.0389407032762076D+00, &
          0.0389407032762076D+00, &
          0.0389407032762076D+00, &
          0.0389407032762076D+00, &
          0.0389407032762076D+00, &
          0.0389407032762076D+00, &
          0.0389407032762076D+00, &
          0.0389407032762076D+00, &
          0.0389407032762076D+00, &
          0.0389407032762076D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule08 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule08() returns the prism rule of precision 8.
!
!  Discussion:
!
!    We suppose we are given a triangular prism P with vertices
!    (1,0,0), (0,1,0), (0,0,0),
!    (1,0,1), (0,1,1), (0,0,1),
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over P is 
!    approximated by Q as follows:
!
!    Q = volume(P) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 May 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order symmetric cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 122, Number 1, pages 148-171, 24 August 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real x(n), y(n), z(n): the coordinates of quadrature points.
!
!    real w(n): the quadrature weights, which sum to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: n_save = 46

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.4600889628137106D+00, &
          0.4600889628137106D+00, &
          0.0798220743725788D+00, &
          0.0534123509369071D+00, &
          0.0534123509369071D+00, &
          0.8931752981261857D+00, &
          0.0472387858397694D+00, &
          0.0472387858397694D+00, &
          0.9055224283204613D+00, &
          0.0472387858397694D+00, &
          0.0472387858397694D+00, &
          0.9055224283204613D+00, &
          0.1740616079243704D+00, &
          0.1740616079243704D+00, &
          0.6518767841512593D+00, &
          0.1740616079243704D+00, &
          0.1740616079243704D+00, &
          0.6518767841512593D+00, &
          0.1597492639425890D+00, &
          0.1597492639425890D+00, &
          0.6805014721148220D+00, &
          0.1597492639425890D+00, &
          0.1597492639425890D+00, &
          0.6805014721148220D+00, &
          0.4585690687909513D+00, &
          0.4585690687909513D+00, &
          0.0828618624180973D+00, &
          0.4585690687909513D+00, &
          0.4585690687909513D+00, &
          0.0828618624180973D+00, &
          0.0085881275077590D+00, &
          0.7285980718010000D+00, &
          0.0085881275077590D+00, &
          0.2628138006912410D+00, &
          0.7285980718010000D+00, &
          0.2628138006912410D+00, &
          0.0085881275077590D+00, &
          0.7285980718010000D+00, &
          0.0085881275077590D+00, &
          0.2628138006912410D+00, &
          0.7285980718010000D+00, &
          0.2628138006912410D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.4600889628137106D+00, &
          0.0798220743725788D+00, &
          0.4600889628137106D+00, &
          0.0534123509369071D+00, &
          0.8931752981261857D+00, &
          0.0534123509369071D+00, &
          0.0472387858397694D+00, &
          0.9055224283204613D+00, &
          0.0472387858397694D+00, &
          0.0472387858397694D+00, &
          0.9055224283204613D+00, &
          0.0472387858397694D+00, &
          0.1740616079243704D+00, &
          0.6518767841512593D+00, &
          0.1740616079243704D+00, &
          0.1740616079243704D+00, &
          0.6518767841512593D+00, &
          0.1740616079243704D+00, &
          0.1597492639425890D+00, &
          0.6805014721148220D+00, &
          0.1597492639425890D+00, &
          0.1597492639425890D+00, &
          0.6805014721148220D+00, &
          0.1597492639425890D+00, &
          0.4585690687909513D+00, &
          0.0828618624180973D+00, &
          0.4585690687909513D+00, &
          0.4585690687909513D+00, &
          0.0828618624180973D+00, &
          0.4585690687909513D+00, &
          0.7285980718010000D+00, &
          0.0085881275077590D+00, &
          0.2628138006912410D+00, &
          0.0085881275077590D+00, &
          0.2628138006912410D+00, &
          0.7285980718010000D+00, &
          0.7285980718010000D+00, &
          0.0085881275077590D+00, &
          0.2628138006912410D+00, &
          0.0085881275077590D+00, &
          0.2628138006912410D+00, &
          0.7285980718010000D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.6122169330529953D+00, &
          0.3877830669470047D+00, &
          0.1590872792145671D+00, &
          0.8409127207854329D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9144077215293400D+00, &
          0.9144077215293400D+00, &
          0.9144077215293400D+00, &
          0.0855922784706599D+00, &
          0.0855922784706599D+00, &
          0.0855922784706599D+00, &
          0.7030421146772952D+00, &
          0.7030421146772952D+00, &
          0.7030421146772952D+00, &
          0.2969578853227048D+00, &
          0.2969578853227048D+00, &
          0.2969578853227048D+00, &
          0.9829325832703272D+00, &
          0.9829325832703272D+00, &
          0.9829325832703272D+00, &
          0.0170674167296728D+00, &
          0.0170674167296728D+00, &
          0.0170674167296728D+00, &
          0.9380979550001272D+00, &
          0.9380979550001272D+00, &
          0.9380979550001272D+00, &
          0.0619020449998729D+00, &
          0.0619020449998729D+00, &
          0.0619020449998729D+00, &
          0.7886751345948129D+00, &
          0.7886751345948129D+00, &
          0.7886751345948129D+00, &
          0.7886751345948129D+00, &
          0.7886751345948129D+00, &
          0.7886751345948129D+00, &
          0.2113248654051871D+00, &
          0.2113248654051871D+00, &
          0.2113248654051871D+00, &
          0.2113248654051871D+00, &
          0.2113248654051871D+00, &
          0.2113248654051871D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          0.0207349366428555D+00, &
          0.0207349366428555D+00, &
          0.0510894271121815D+00, &
          0.0510894271121815D+00, &
          0.0526997449065072D+00, &
          0.0526997449065072D+00, &
          0.0526997449065072D+00, &
          0.0181522487841497D+00, &
          0.0181522487841497D+00, &
          0.0181522487841497D+00, &
          0.0071453577535386D+00, &
          0.0071453577535386D+00, &
          0.0071453577535386D+00, &
          0.0071453577535386D+00, &
          0.0071453577535386D+00, &
          0.0071453577535386D+00, &
          0.0406292626589893D+00, &
          0.0406292626589893D+00, &
          0.0406292626589893D+00, &
          0.0406292626589893D+00, &
          0.0406292626589893D+00, &
          0.0406292626589893D+00, &
          0.0111430273484653D+00, &
          0.0111430273484653D+00, &
          0.0111430273484653D+00, &
          0.0111430273484653D+00, &
          0.0111430273484653D+00, &
          0.0111430273484653D+00, &
          0.0210865092935865D+00, &
          0.0210865092935865D+00, &
          0.0210865092935865D+00, &
          0.0210865092935865D+00, &
          0.0210865092935865D+00, &
          0.0210865092935865D+00, &
          0.0136475290908731D+00, &
          0.0136475290908731D+00, &
          0.0136475290908731D+00, &
          0.0136475290908731D+00, &
          0.0136475290908731D+00, &
          0.0136475290908731D+00, &
          0.0136475290908731D+00, &
          0.0136475290908731D+00, &
          0.0136475290908731D+00, &
          0.0136475290908731D+00, &
          0.0136475290908731D+00, &
          0.0136475290908731D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule09 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule09() returns the prism rule of precision 9.
!
!  Discussion:
!
!    We suppose we are given a triangular prism P with vertices
!    (1,0,0), (0,1,0), (0,0,0),
!    (1,0,1), (0,1,1), (0,0,1),
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over P is 
!    approximated by Q as follows:
!
!    Q = volume(P) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 May 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order symmetric cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 122, Number 1, pages 148-171, 24 August 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real x(n), y(n), z(n): the coordinates of quadrature points.
!
!    real w(n): the quadrature weights, which sum to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: n_save = 59

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.4938630638969568D+00, &
          0.4938630638969568D+00, &
          0.0122738722060863D+00, &
          0.2362540169543293D+00, &
          0.2362540169543293D+00, &
          0.5274919660913414D+00, &
          0.2362540169543293D+00, &
          0.2362540169543293D+00, &
          0.5274919660913414D+00, &
          0.0722383394163824D+00, &
          0.0722383394163824D+00, &
          0.8555233211672353D+00, &
          0.0722383394163824D+00, &
          0.0722383394163824D+00, &
          0.8555233211672353D+00, &
          0.4383137607101617D+00, &
          0.4383137607101617D+00, &
          0.1233724785796766D+00, &
          0.4383137607101617D+00, &
          0.4383137607101617D+00, &
          0.1233724785796766D+00, &
          0.0364340164940779D+00, &
          0.0364340164940779D+00, &
          0.9271319670118442D+00, &
          0.0364340164940779D+00, &
          0.0364340164940779D+00, &
          0.9271319670118442D+00, &
          0.4828779929693860D+00, &
          0.4828779929693860D+00, &
          0.0342440140612279D+00, &
          0.4828779929693860D+00, &
          0.4828779929693860D+00, &
          0.0342440140612279D+00, &
          0.1628698857202373D+00, &
          0.1628698857202373D+00, &
          0.6742602285595254D+00, &
          0.1628698857202373D+00, &
          0.1628698857202373D+00, &
          0.6742602285595254D+00, &
          0.8213377527237301D+00, &
          0.1626087609745086D+00, &
          0.8213377527237301D+00, &
          0.0160534863017613D+00, &
          0.1626087609745086D+00, &
          0.0160534863017613D+00, &
          0.0424944495063928D+00, &
          0.2561926710584905D+00, &
          0.0424944495063928D+00, &
          0.7013128794351167D+00, &
          0.2561926710584905D+00, &
          0.7013128794351167D+00, &
          0.0424944495063928D+00, &
          0.2561926710584905D+00, &
          0.0424944495063928D+00, &
          0.7013128794351167D+00, &
          0.2561926710584905D+00, &
          0.7013128794351167D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.4938630638969568D+00, &
          0.0122738722060863D+00, &
          0.4938630638969568D+00, &
          0.2362540169543293D+00, &
          0.5274919660913414D+00, &
          0.2362540169543293D+00, &
          0.2362540169543293D+00, &
          0.5274919660913414D+00, &
          0.2362540169543293D+00, &
          0.0722383394163824D+00, &
          0.8555233211672353D+00, &
          0.0722383394163824D+00, &
          0.0722383394163824D+00, &
          0.8555233211672353D+00, &
          0.0722383394163824D+00, &
          0.4383137607101617D+00, &
          0.1233724785796766D+00, &
          0.4383137607101617D+00, &
          0.4383137607101617D+00, &
          0.1233724785796766D+00, &
          0.4383137607101617D+00, &
          0.0364340164940779D+00, &
          0.9271319670118442D+00, &
          0.0364340164940779D+00, &
          0.0364340164940779D+00, &
          0.9271319670118442D+00, &
          0.0364340164940779D+00, &
          0.4828779929693860D+00, &
          0.0342440140612279D+00, &
          0.4828779929693860D+00, &
          0.4828779929693860D+00, &
          0.0342440140612279D+00, &
          0.4828779929693860D+00, &
          0.1628698857202373D+00, &
          0.6742602285595254D+00, &
          0.1628698857202373D+00, &
          0.1628698857202373D+00, &
          0.6742602285595254D+00, &
          0.1628698857202373D+00, &
          0.1626087609745086D+00, &
          0.8213377527237301D+00, &
          0.0160534863017613D+00, &
          0.8213377527237301D+00, &
          0.0160534863017613D+00, &
          0.1626087609745086D+00, &
          0.2561926710584905D+00, &
          0.0424944495063928D+00, &
          0.7013128794351167D+00, &
          0.0424944495063928D+00, &
          0.7013128794351167D+00, &
          0.2561926710584905D+00, &
          0.2561926710584905D+00, &
          0.0424944495063928D+00, &
          0.7013128794351167D+00, &
          0.0424944495063928D+00, &
          0.7013128794351167D+00, &
          0.2561926710584905D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.3624907904961116D+00, &
          0.6375092095038883D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9356193788144789D+00, &
          0.9356193788144789D+00, &
          0.9356193788144789D+00, &
          0.0643806211855212D+00, &
          0.0643806211855212D+00, &
          0.0643806211855212D+00, &
          0.9775968701180802D+00, &
          0.9775968701180802D+00, &
          0.9775968701180802D+00, &
          0.0224031298819198D+00, &
          0.0224031298819198D+00, &
          0.0224031298819198D+00, &
          0.7182058908565040D+00, &
          0.7182058908565040D+00, &
          0.7182058908565040D+00, &
          0.2817941091434960D+00, &
          0.2817941091434960D+00, &
          0.2817941091434960D+00, &
          0.7515028180283210D+00, &
          0.7515028180283210D+00, &
          0.7515028180283210D+00, &
          0.2484971819716791D+00, &
          0.2484971819716791D+00, &
          0.2484971819716791D+00, &
          0.0091271384992068D+00, &
          0.0091271384992068D+00, &
          0.0091271384992068D+00, &
          0.9908728615007932D+00, &
          0.9908728615007932D+00, &
          0.9908728615007932D+00, &
          0.4125665631901645D+00, &
          0.4125665631901645D+00, &
          0.4125665631901645D+00, &
          0.5874334368098355D+00, &
          0.5874334368098355D+00, &
          0.5874334368098355D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.8418079279883068D+00, &
          0.8418079279883068D+00, &
          0.8418079279883068D+00, &
          0.8418079279883068D+00, &
          0.8418079279883068D+00, &
          0.8418079279883068D+00, &
          0.1581920720116932D+00, &
          0.1581920720116932D+00, &
          0.1581920720116932D+00, &
          0.1581920720116932D+00, &
          0.1581920720116932D+00, &
          0.1581920720116932D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          0.0342834712909237D+00, &
          0.0342834712909237D+00, &
          0.0163927977044536D+00, &
          0.0163927977044536D+00, &
          0.0163927977044536D+00, &
          0.0234717806987032D+00, &
          0.0234717806987032D+00, &
          0.0234717806987032D+00, &
          0.0234717806987032D+00, &
          0.0234717806987032D+00, &
          0.0234717806987032D+00, &
          0.0057618481582789D+00, &
          0.0057618481582789D+00, &
          0.0057618481582789D+00, &
          0.0057618481582789D+00, &
          0.0057618481582789D+00, &
          0.0057618481582789D+00, &
          0.0349247930577802D+00, &
          0.0349247930577802D+00, &
          0.0349247930577802D+00, &
          0.0349247930577802D+00, &
          0.0349247930577802D+00, &
          0.0349247930577802D+00, &
          0.0073738996232123D+00, &
          0.0073738996232123D+00, &
          0.0073738996232123D+00, &
          0.0073738996232123D+00, &
          0.0073738996232123D+00, &
          0.0073738996232123D+00, &
          0.0057224608464489D+00, &
          0.0057224608464489D+00, &
          0.0057224608464489D+00, &
          0.0057224608464489D+00, &
          0.0057224608464489D+00, &
          0.0057224608464489D+00, &
          0.0243387430144453D+00, &
          0.0243387430144453D+00, &
          0.0243387430144453D+00, &
          0.0243387430144453D+00, &
          0.0243387430144453D+00, &
          0.0243387430144453D+00, &
          0.0094129636631352D+00, &
          0.0094129636631352D+00, &
          0.0094129636631352D+00, &
          0.0094129636631352D+00, &
          0.0094129636631352D+00, &
          0.0094129636631352D+00, &
          0.0180179774943973D+00, &
          0.0180179774943973D+00, &
          0.0180179774943973D+00, &
          0.0180179774943973D+00, &
          0.0180179774943973D+00, &
          0.0180179774943973D+00, &
          0.0180179774943973D+00, &
          0.0180179774943973D+00, &
          0.0180179774943973D+00, &
          0.0180179774943973D+00, &
          0.0180179774943973D+00, &
          0.0180179774943973D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule10 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule10() returns the prism rule of precision 10.
!
!  Discussion:
!
!    We suppose we are given a triangular prism P with vertices
!    (1,0,0), (0,1,0), (0,0,0),
!    (1,0,1), (0,1,1), (0,0,1),
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over P is 
!    approximated by Q as follows:
!
!    Q = volume(P) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 May 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order symmetric cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 122, Number 1, pages 148-171, 24 August 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real x(n), y(n), z(n): the coordinates of quadrature points.
!
!    real w(n): the quadrature weights, which sum to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: n_save = 84

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.4300104731739727D+00, &
          0.4300104731739727D+00, &
          0.1399790536520547D+00, &
          0.1095696683547513D+00, &
          0.1095696683547513D+00, &
          0.7808606632904974D+00, &
          0.0158148366313868D+00, &
          0.0158148366313868D+00, &
          0.9683703267372263D+00, &
          0.4356018236697312D+00, &
          0.4356018236697312D+00, &
          0.1287963526605377D+00, &
          0.4356018236697312D+00, &
          0.4356018236697312D+00, &
          0.1287963526605377D+00, &
          0.1485976165307029D+00, &
          0.1485976165307029D+00, &
          0.7028047669385943D+00, &
          0.1485976165307029D+00, &
          0.1485976165307029D+00, &
          0.7028047669385943D+00, &
          0.0519664835822357D+00, &
          0.0519664835822357D+00, &
          0.8960670328355285D+00, &
          0.0519664835822357D+00, &
          0.0519664835822357D+00, &
          0.8960670328355285D+00, &
          0.4990219694442680D+00, &
          0.4990219694442680D+00, &
          0.0019560611114641D+00, &
          0.4990219694442680D+00, &
          0.4990219694442680D+00, &
          0.0019560611114641D+00, &
          0.0415525016027702D+00, &
          0.0415525016027702D+00, &
          0.9168949967944596D+00, &
          0.0415525016027702D+00, &
          0.0415525016027702D+00, &
          0.9168949967944596D+00, &
          0.2343945196820784D+00, &
          0.2343945196820784D+00, &
          0.5312109606358433D+00, &
          0.2343945196820784D+00, &
          0.2343945196820784D+00, &
          0.5312109606358433D+00, &
          0.0528115168465621D+00, &
          0.2716521744885937D+00, &
          0.0528115168465621D+00, &
          0.6755363086648442D+00, &
          0.2716521744885937D+00, &
          0.6755363086648442D+00, &
          0.0528115168465621D+00, &
          0.2716521744885937D+00, &
          0.0528115168465621D+00, &
          0.6755363086648442D+00, &
          0.2716521744885937D+00, &
          0.6755363086648442D+00, &
          0.1675738337212976D+00, &
          0.0101912520986929D+00, &
          0.1675738337212976D+00, &
          0.8222349141800095D+00, &
          0.0101912520986929D+00, &
          0.8222349141800095D+00, &
          0.1675738337212976D+00, &
          0.0101912520986929D+00, &
          0.1675738337212976D+00, &
          0.8222349141800095D+00, &
          0.0101912520986929D+00, &
          0.8222349141800095D+00, &
          0.3291869417398026D+00, &
          0.0509081627669518D+00, &
          0.3291869417398026D+00, &
          0.6199048954932456D+00, &
          0.0509081627669518D+00, &
          0.6199048954932456D+00, &
          0.3291869417398026D+00, &
          0.0509081627669518D+00, &
          0.3291869417398026D+00, &
          0.6199048954932456D+00, &
          0.0509081627669518D+00, &
          0.6199048954932456D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.4300104731739727D+00, &
          0.1399790536520547D+00, &
          0.4300104731739727D+00, &
          0.1095696683547513D+00, &
          0.7808606632904974D+00, &
          0.1095696683547513D+00, &
          0.0158148366313868D+00, &
          0.9683703267372263D+00, &
          0.0158148366313868D+00, &
          0.4356018236697312D+00, &
          0.1287963526605377D+00, &
          0.4356018236697312D+00, &
          0.4356018236697312D+00, &
          0.1287963526605377D+00, &
          0.4356018236697312D+00, &
          0.1485976165307029D+00, &
          0.7028047669385943D+00, &
          0.1485976165307029D+00, &
          0.1485976165307029D+00, &
          0.7028047669385943D+00, &
          0.1485976165307029D+00, &
          0.0519664835822357D+00, &
          0.8960670328355285D+00, &
          0.0519664835822357D+00, &
          0.0519664835822357D+00, &
          0.8960670328355285D+00, &
          0.0519664835822357D+00, &
          0.4990219694442680D+00, &
          0.0019560611114641D+00, &
          0.4990219694442680D+00, &
          0.4990219694442680D+00, &
          0.0019560611114641D+00, &
          0.4990219694442680D+00, &
          0.0415525016027702D+00, &
          0.9168949967944596D+00, &
          0.0415525016027702D+00, &
          0.0415525016027702D+00, &
          0.9168949967944596D+00, &
          0.0415525016027702D+00, &
          0.2343945196820784D+00, &
          0.5312109606358433D+00, &
          0.2343945196820784D+00, &
          0.2343945196820784D+00, &
          0.5312109606358433D+00, &
          0.2343945196820784D+00, &
          0.2716521744885937D+00, &
          0.0528115168465621D+00, &
          0.6755363086648442D+00, &
          0.0528115168465621D+00, &
          0.6755363086648442D+00, &
          0.2716521744885937D+00, &
          0.2716521744885937D+00, &
          0.0528115168465621D+00, &
          0.6755363086648442D+00, &
          0.0528115168465621D+00, &
          0.6755363086648442D+00, &
          0.2716521744885937D+00, &
          0.0101912520986929D+00, &
          0.1675738337212976D+00, &
          0.8222349141800095D+00, &
          0.1675738337212976D+00, &
          0.8222349141800095D+00, &
          0.0101912520986929D+00, &
          0.0101912520986929D+00, &
          0.1675738337212976D+00, &
          0.8222349141800095D+00, &
          0.1675738337212976D+00, &
          0.8222349141800095D+00, &
          0.0101912520986929D+00, &
          0.0509081627669518D+00, &
          0.3291869417398026D+00, &
          0.6199048954932456D+00, &
          0.3291869417398026D+00, &
          0.6199048954932456D+00, &
          0.0509081627669518D+00, &
          0.0509081627669518D+00, &
          0.3291869417398026D+00, &
          0.6199048954932456D+00, &
          0.3291869417398026D+00, &
          0.6199048954932456D+00, &
          0.0509081627669518D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.5000000000000000D+00, &
          0.9825433408084738D+00, &
          0.0174566591915262D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.8914517952143866D+00, &
          0.8914517952143866D+00, &
          0.8914517952143866D+00, &
          0.1085482047856134D+00, &
          0.1085482047856134D+00, &
          0.1085482047856134D+00, &
          0.8836628766030707D+00, &
          0.8836628766030707D+00, &
          0.8836628766030707D+00, &
          0.1163371233969292D+00, &
          0.1163371233969292D+00, &
          0.1163371233969292D+00, &
          0.7446456974623019D+00, &
          0.7446456974623019D+00, &
          0.7446456974623019D+00, &
          0.2553543025376981D+00, &
          0.2553543025376981D+00, &
          0.2553543025376981D+00, &
          0.8458358146384862D+00, &
          0.8458358146384862D+00, &
          0.8458358146384862D+00, &
          0.1541641853615138D+00, &
          0.1541641853615138D+00, &
          0.1541641853615138D+00, &
          0.9455005421560339D+00, &
          0.9455005421560339D+00, &
          0.9455005421560339D+00, &
          0.0544994578439661D+00, &
          0.0544994578439661D+00, &
          0.0544994578439661D+00, &
          0.7325671383641421D+00, &
          0.7325671383641421D+00, &
          0.7325671383641421D+00, &
          0.2674328616358579D+00, &
          0.2674328616358579D+00, &
          0.2674328616358579D+00, &
          0.9783107845287782D+00, &
          0.9783107845287782D+00, &
          0.9783107845287782D+00, &
          0.9783107845287782D+00, &
          0.9783107845287782D+00, &
          0.9783107845287782D+00, &
          0.0216892154712218D+00, &
          0.0216892154712218D+00, &
          0.0216892154712218D+00, &
          0.0216892154712218D+00, &
          0.0216892154712218D+00, &
          0.0216892154712218D+00, &
          0.7695898729425880D+00, &
          0.7695898729425880D+00, &
          0.7695898729425880D+00, &
          0.7695898729425880D+00, &
          0.7695898729425880D+00, &
          0.7695898729425880D+00, &
          0.2304101270574120D+00, &
          0.2304101270574120D+00, &
          0.2304101270574120D+00, &
          0.2304101270574120D+00, &
          0.2304101270574120D+00, &
          0.2304101270574120D+00, &
          0.3678576424324543D+00, &
          0.3678576424324543D+00, &
          0.3678576424324543D+00, &
          0.3678576424324543D+00, &
          0.3678576424324543D+00, &
          0.3678576424324543D+00, &
          0.6321423575675457D+00, &
          0.6321423575675457D+00, &
          0.6321423575675457D+00, &
          0.6321423575675457D+00, &
          0.6321423575675457D+00, &
          0.6321423575675457D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          0.0321960882383680D+00, &
          0.0132395953233874D+00, &
          0.0132395953233874D+00, &
          0.0166307172641394D+00, &
          0.0166307172641394D+00, &
          0.0166307172641394D+00, &
          0.0223711737567938D+00, &
          0.0223711737567938D+00, &
          0.0223711737567938D+00, &
          0.0029651442680284D+00, &
          0.0029651442680284D+00, &
          0.0029651442680284D+00, &
          0.0204122071267436D+00, &
          0.0204122071267436D+00, &
          0.0204122071267436D+00, &
          0.0204122071267436D+00, &
          0.0204122071267436D+00, &
          0.0204122071267436D+00, &
          0.0141443496976585D+00, &
          0.0141443496976585D+00, &
          0.0141443496976585D+00, &
          0.0141443496976585D+00, &
          0.0141443496976585D+00, &
          0.0141443496976585D+00, &
          0.0035457088242710D+00, &
          0.0035457088242710D+00, &
          0.0035457088242710D+00, &
          0.0035457088242710D+00, &
          0.0035457088242710D+00, &
          0.0035457088242710D+00, &
          0.0060580804829928D+00, &
          0.0060580804829928D+00, &
          0.0060580804829928D+00, &
          0.0060580804829928D+00, &
          0.0060580804829928D+00, &
          0.0060580804829928D+00, &
          0.0035558880399586D+00, &
          0.0035558880399586D+00, &
          0.0035558880399586D+00, &
          0.0035558880399586D+00, &
          0.0035558880399586D+00, &
          0.0035558880399586D+00, &
          0.0303810246652175D+00, &
          0.0303810246652175D+00, &
          0.0303810246652175D+00, &
          0.0303810246652175D+00, &
          0.0303810246652175D+00, &
          0.0303810246652175D+00, &
          0.0062118567570921D+00, &
          0.0062118567570921D+00, &
          0.0062118567570921D+00, &
          0.0062118567570921D+00, &
          0.0062118567570921D+00, &
          0.0062118567570921D+00, &
          0.0062118567570921D+00, &
          0.0062118567570921D+00, &
          0.0062118567570921D+00, &
          0.0062118567570921D+00, &
          0.0062118567570921D+00, &
          0.0062118567570921D+00, &
          0.0066184191661461D+00, &
          0.0066184191661461D+00, &
          0.0066184191661461D+00, &
          0.0066184191661461D+00, &
          0.0066184191661461D+00, &
          0.0066184191661461D+00, &
          0.0066184191661461D+00, &
          0.0066184191661461D+00, &
          0.0066184191661461D+00, &
          0.0066184191661461D+00, &
          0.0066184191661461D+00, &
          0.0066184191661461D+00, &
          0.0160730625956718D+00, &
          0.0160730625956718D+00, &
          0.0160730625956718D+00, &
          0.0160730625956718D+00, &
          0.0160730625956718D+00, &
          0.0160730625956718D+00, &
          0.0160730625956718D+00, &
          0.0160730625956718D+00, &
          0.0160730625956718D+00, &
          0.0160730625956718D+00, &
          0.0160730625956718D+00, &
          0.0160730625956718D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule11 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule11() returns the prism rule of precision 11.
!
!  Discussion:
!
!    We suppose we are given a triangular prism P with vertices
!    (1,0,0), (0,1,0), (0,0,0),
!    (1,0,1), (0,1,1), (0,0,1),
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over P is 
!    approximated by Q as follows:
!
!    Q = volume(P) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 May 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order symmetric cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 122, Number 1, pages 148-171, 24 August 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real x(n), y(n), z(n): the coordinates of quadrature points.
!
!    real w(n): the quadrature weights, which sum to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: n_save = 99

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.4269221881685347D+00, &
          0.4269221881685347D+00, &
          0.1461556236629306D+00, &
          0.4396327659106838D+00, &
          0.4396327659106838D+00, &
          0.1207344681786323D+00, &
          0.4396327659106838D+00, &
          0.4396327659106838D+00, &
          0.1207344681786323D+00, &
          0.2009128761138035D+00, &
          0.2009128761138035D+00, &
          0.5981742477723929D+00, &
          0.2009128761138035D+00, &
          0.2009128761138035D+00, &
          0.5981742477723929D+00, &
          0.0215850980431769D+00, &
          0.0215850980431769D+00, &
          0.9568298039136461D+00, &
          0.0215850980431769D+00, &
          0.0215850980431769D+00, &
          0.9568298039136461D+00, &
          0.4941546402080231D+00, &
          0.4941546402080231D+00, &
          0.0116907195839538D+00, &
          0.4941546402080231D+00, &
          0.4941546402080231D+00, &
          0.0116907195839538D+00, &
          0.0599602772654436D+00, &
          0.0599602772654436D+00, &
          0.8800794454691128D+00, &
          0.0599602772654436D+00, &
          0.0599602772654436D+00, &
          0.8800794454691128D+00, &
          0.2328102236807494D+00, &
          0.2328102236807494D+00, &
          0.5343795526385012D+00, &
          0.2328102236807494D+00, &
          0.2328102236807494D+00, &
          0.5343795526385012D+00, &
          0.1213201391549184D+00, &
          0.1213201391549184D+00, &
          0.7573597216901632D+00, &
          0.1213201391549184D+00, &
          0.1213201391549184D+00, &
          0.7573597216901632D+00, &
          0.4987543911906001D+00, &
          0.4987543911906001D+00, &
          0.0024912176187997D+00, &
          0.4987543911906001D+00, &
          0.4987543911906001D+00, &
          0.0024912176187997D+00, &
          0.1046140524481813D+00, &
          0.0165184963342511D+00, &
          0.1046140524481813D+00, &
          0.8788674512175675D+00, &
          0.0165184963342511D+00, &
          0.8788674512175675D+00, &
          0.0719632569755848D+00, &
          0.2979854667459965D+00, &
          0.0719632569755848D+00, &
          0.6300512762784186D+00, &
          0.2979854667459965D+00, &
          0.6300512762784186D+00, &
          0.0719632569755848D+00, &
          0.2979854667459965D+00, &
          0.0719632569755848D+00, &
          0.6300512762784186D+00, &
          0.2979854667459965D+00, &
          0.6300512762784186D+00, &
          0.1998026706474004D+00, &
          0.0172094825510263D+00, &
          0.1998026706474004D+00, &
          0.7829878468015733D+00, &
          0.0172094825510263D+00, &
          0.7829878468015733D+00, &
          0.1998026706474004D+00, &
          0.0172094825510263D+00, &
          0.1998026706474004D+00, &
          0.7829878468015733D+00, &
          0.0172094825510263D+00, &
          0.7829878468015733D+00, &
          0.3197359768880742D+00, &
          0.0462852386525127D+00, &
          0.3197359768880742D+00, &
          0.6339787844594131D+00, &
          0.0462852386525127D+00, &
          0.6339787844594131D+00, &
          0.3197359768880742D+00, &
          0.0462852386525127D+00, &
          0.3197359768880742D+00, &
          0.6339787844594131D+00, &
          0.0462852386525127D+00, &
          0.6339787844594131D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.4269221881685347D+00, &
          0.1461556236629306D+00, &
          0.4269221881685347D+00, &
          0.4396327659106838D+00, &
          0.1207344681786323D+00, &
          0.4396327659106838D+00, &
          0.4396327659106838D+00, &
          0.1207344681786323D+00, &
          0.4396327659106838D+00, &
          0.2009128761138035D+00, &
          0.5981742477723929D+00, &
          0.2009128761138035D+00, &
          0.2009128761138035D+00, &
          0.5981742477723929D+00, &
          0.2009128761138035D+00, &
          0.0215850980431769D+00, &
          0.9568298039136461D+00, &
          0.0215850980431769D+00, &
          0.0215850980431769D+00, &
          0.9568298039136461D+00, &
          0.0215850980431769D+00, &
          0.4941546402080231D+00, &
          0.0116907195839538D+00, &
          0.4941546402080231D+00, &
          0.4941546402080231D+00, &
          0.0116907195839538D+00, &
          0.4941546402080231D+00, &
          0.0599602772654436D+00, &
          0.8800794454691128D+00, &
          0.0599602772654436D+00, &
          0.0599602772654436D+00, &
          0.8800794454691128D+00, &
          0.0599602772654436D+00, &
          0.2328102236807494D+00, &
          0.5343795526385012D+00, &
          0.2328102236807494D+00, &
          0.2328102236807494D+00, &
          0.5343795526385012D+00, &
          0.2328102236807494D+00, &
          0.1213201391549184D+00, &
          0.7573597216901632D+00, &
          0.1213201391549184D+00, &
          0.1213201391549184D+00, &
          0.7573597216901632D+00, &
          0.1213201391549184D+00, &
          0.4987543911906001D+00, &
          0.0024912176187997D+00, &
          0.4987543911906001D+00, &
          0.4987543911906001D+00, &
          0.0024912176187997D+00, &
          0.4987543911906001D+00, &
          0.0165184963342511D+00, &
          0.1046140524481813D+00, &
          0.8788674512175675D+00, &
          0.1046140524481813D+00, &
          0.8788674512175675D+00, &
          0.0165184963342511D+00, &
          0.2979854667459965D+00, &
          0.0719632569755848D+00, &
          0.6300512762784186D+00, &
          0.0719632569755848D+00, &
          0.6300512762784186D+00, &
          0.2979854667459965D+00, &
          0.2979854667459965D+00, &
          0.0719632569755848D+00, &
          0.6300512762784186D+00, &
          0.0719632569755848D+00, &
          0.6300512762784186D+00, &
          0.2979854667459965D+00, &
          0.0172094825510263D+00, &
          0.1998026706474004D+00, &
          0.7829878468015733D+00, &
          0.1998026706474004D+00, &
          0.7829878468015733D+00, &
          0.0172094825510263D+00, &
          0.0172094825510263D+00, &
          0.1998026706474004D+00, &
          0.7829878468015733D+00, &
          0.1998026706474004D+00, &
          0.7829878468015733D+00, &
          0.0172094825510263D+00, &
          0.0462852386525127D+00, &
          0.3197359768880742D+00, &
          0.6339787844594131D+00, &
          0.3197359768880742D+00, &
          0.6339787844594131D+00, &
          0.0462852386525127D+00, &
          0.0462852386525127D+00, &
          0.3197359768880742D+00, &
          0.6339787844594131D+00, &
          0.3197359768880742D+00, &
          0.6339787844594131D+00, &
          0.0462852386525127D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.2225225095750167D+00, &
          0.7774774904249833D+00, &
          0.9852161050240122D+00, &
          0.0147838949759878D+00, &
          0.4588747274117720D+00, &
          0.5411252725882280D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.8379989630425995D+00, &
          0.8379989630425995D+00, &
          0.8379989630425995D+00, &
          0.1620010369574005D+00, &
          0.1620010369574005D+00, &
          0.1620010369574005D+00, &
          0.9180267718000730D+00, &
          0.9180267718000730D+00, &
          0.9180267718000730D+00, &
          0.0819732281999270D+00, &
          0.0819732281999270D+00, &
          0.0819732281999270D+00, &
          0.7930952762440336D+00, &
          0.7930952762440336D+00, &
          0.7930952762440336D+00, &
          0.2069047237559664D+00, &
          0.2069047237559664D+00, &
          0.2069047237559664D+00, &
          0.9873500827745887D+00, &
          0.9873500827745887D+00, &
          0.9873500827745887D+00, &
          0.0126499172254113D+00, &
          0.0126499172254113D+00, &
          0.0126499172254113D+00, &
          0.9725136338749176D+00, &
          0.9725136338749176D+00, &
          0.9725136338749176D+00, &
          0.0274863661250824D+00, &
          0.0274863661250824D+00, &
          0.0274863661250824D+00, &
          0.6317612421960503D+00, &
          0.6317612421960503D+00, &
          0.6317612421960503D+00, &
          0.3682387578039497D+00, &
          0.3682387578039497D+00, &
          0.3682387578039497D+00, &
          0.7219895510260312D+00, &
          0.7219895510260312D+00, &
          0.7219895510260312D+00, &
          0.2780104489739688D+00, &
          0.2780104489739688D+00, &
          0.2780104489739688D+00, &
          0.7875076710403349D+00, &
          0.7875076710403349D+00, &
          0.7875076710403349D+00, &
          0.2124923289596651D+00, &
          0.2124923289596651D+00, &
          0.2124923289596651D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9685666190309952D+00, &
          0.9685666190309952D+00, &
          0.9685666190309952D+00, &
          0.9685666190309952D+00, &
          0.9685666190309952D+00, &
          0.9685666190309952D+00, &
          0.0314333809690049D+00, &
          0.0314333809690049D+00, &
          0.0314333809690049D+00, &
          0.0314333809690049D+00, &
          0.0314333809690049D+00, &
          0.0314333809690049D+00, &
          0.8651322718570146D+00, &
          0.8651322718570146D+00, &
          0.8651322718570146D+00, &
          0.8651322718570146D+00, &
          0.8651322718570146D+00, &
          0.8651322718570146D+00, &
          0.1348677281429854D+00, &
          0.1348677281429854D+00, &
          0.1348677281429854D+00, &
          0.1348677281429854D+00, &
          0.1348677281429854D+00, &
          0.1348677281429854D+00, &
          0.3827016872275200D+00, &
          0.3827016872275200D+00, &
          0.3827016872275200D+00, &
          0.3827016872275200D+00, &
          0.3827016872275200D+00, &
          0.3827016872275200D+00, &
          0.6172983127724800D+00, &
          0.6172983127724800D+00, &
          0.6172983127724800D+00, &
          0.6172983127724800D+00, &
          0.6172983127724800D+00, &
          0.6172983127724800D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          0.0216801440673047D+00, &
          0.0216801440673047D+00, &
          0.0104576079858723D+00, &
          0.0104576079858723D+00, &
          0.0034934625077134D+00, &
          0.0034934625077134D+00, &
          0.0161061041817141D+00, &
          0.0161061041817141D+00, &
          0.0161061041817141D+00, &
          0.0220331028546643D+00, &
          0.0220331028546643D+00, &
          0.0220331028546643D+00, &
          0.0220331028546643D+00, &
          0.0220331028546643D+00, &
          0.0220331028546643D+00, &
          0.0127441869618096D+00, &
          0.0127441869618096D+00, &
          0.0127441869618096D+00, &
          0.0127441869618096D+00, &
          0.0127441869618096D+00, &
          0.0127441869618096D+00, &
          0.0028425261018310D+00, &
          0.0028425261018310D+00, &
          0.0028425261018310D+00, &
          0.0028425261018310D+00, &
          0.0028425261018310D+00, &
          0.0028425261018310D+00, &
          0.0015098577228106D+00, &
          0.0015098577228106D+00, &
          0.0015098577228106D+00, &
          0.0015098577228106D+00, &
          0.0015098577228106D+00, &
          0.0015098577228106D+00, &
          0.0036504270647864D+00, &
          0.0036504270647864D+00, &
          0.0036504270647864D+00, &
          0.0036504270647864D+00, &
          0.0036504270647864D+00, &
          0.0036504270647864D+00, &
          0.0210029067169834D+00, &
          0.0210029067169834D+00, &
          0.0210029067169834D+00, &
          0.0210029067169834D+00, &
          0.0210029067169834D+00, &
          0.0210029067169834D+00, &
          0.0187619585366729D+00, &
          0.0187619585366729D+00, &
          0.0187619585366729D+00, &
          0.0187619585366729D+00, &
          0.0187619585366729D+00, &
          0.0187619585366729D+00, &
          0.0051948191291646D+00, &
          0.0051948191291646D+00, &
          0.0051948191291646D+00, &
          0.0051948191291646D+00, &
          0.0051948191291646D+00, &
          0.0051948191291646D+00, &
          0.0070505879387705D+00, &
          0.0070505879387705D+00, &
          0.0070505879387705D+00, &
          0.0070505879387705D+00, &
          0.0070505879387705D+00, &
          0.0070505879387705D+00, &
          0.0059012699482948D+00, &
          0.0059012699482948D+00, &
          0.0059012699482948D+00, &
          0.0059012699482948D+00, &
          0.0059012699482948D+00, &
          0.0059012699482948D+00, &
          0.0059012699482948D+00, &
          0.0059012699482948D+00, &
          0.0059012699482948D+00, &
          0.0059012699482948D+00, &
          0.0059012699482948D+00, &
          0.0059012699482948D+00, &
          0.0062126654655306D+00, &
          0.0062126654655306D+00, &
          0.0062126654655306D+00, &
          0.0062126654655306D+00, &
          0.0062126654655306D+00, &
          0.0062126654655306D+00, &
          0.0062126654655306D+00, &
          0.0062126654655306D+00, &
          0.0062126654655306D+00, &
          0.0062126654655306D+00, &
          0.0062126654655306D+00, &
          0.0062126654655306D+00, &
          0.0138591496001844D+00, &
          0.0138591496001844D+00, &
          0.0138591496001844D+00, &
          0.0138591496001844D+00, &
          0.0138591496001844D+00, &
          0.0138591496001844D+00, &
          0.0138591496001844D+00, &
          0.0138591496001844D+00, &
          0.0138591496001844D+00, &
          0.0138591496001844D+00, &
          0.0138591496001844D+00, &
          0.0138591496001844D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule12 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule12() returns the prism rule of precision 12.
!
!  Discussion:
!
!    We suppose we are given a triangular prism P with vertices
!    (1,0,0), (0,1,0), (0,0,0),
!    (1,0,1), (0,1,1), (0,0,1),
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over P is 
!    approximated by Q as follows:
!
!    Q = volume(P) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 May 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order symmetric cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 122, Number 1, pages 148-171, 24 August 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real x(n), y(n), z(n): the coordinates of quadrature points.
!
!    real w(n): the quadrature weights, which sum to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: n_save = 136

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.0485337200788840D+00, &
          0.0485337200788840D+00, &
          0.9029325598422321D+00, &
          0.0003340541016131D+00, &
          0.0003340541016131D+00, &
          0.9993318917967738D+00, &
          0.4849074297077183D+00, &
          0.4849074297077183D+00, &
          0.0301851405845634D+00, &
          0.4849074297077183D+00, &
          0.4849074297077183D+00, &
          0.0301851405845634D+00, &
          0.1497566246841529D+00, &
          0.1497566246841529D+00, &
          0.7004867506316943D+00, &
          0.1497566246841529D+00, &
          0.1497566246841529D+00, &
          0.7004867506316943D+00, &
          0.2344536517846724D+00, &
          0.2344536517846724D+00, &
          0.5310926964306552D+00, &
          0.2344536517846724D+00, &
          0.2344536517846724D+00, &
          0.5310926964306552D+00, &
          0.4409030162469282D+00, &
          0.4409030162469282D+00, &
          0.1181939675061436D+00, &
          0.4409030162469282D+00, &
          0.4409030162469282D+00, &
          0.1181939675061436D+00, &
          0.4870412467653055D+00, &
          0.4870412467653055D+00, &
          0.0259175064693890D+00, &
          0.4870412467653055D+00, &
          0.4870412467653055D+00, &
          0.0259175064693890D+00, &
          0.0248638193002129D+00, &
          0.0248638193002129D+00, &
          0.9502723613995742D+00, &
          0.0248638193002129D+00, &
          0.0248638193002129D+00, &
          0.9502723613995742D+00, &
          0.1118542147928236D+00, &
          0.1118542147928236D+00, &
          0.7762915704143528D+00, &
          0.1118542147928236D+00, &
          0.1118542147928236D+00, &
          0.7762915704143528D+00, &
          0.1207867185816364D+00, &
          0.7097236881695401D+00, &
          0.1207867185816364D+00, &
          0.1694895932488235D+00, &
          0.7097236881695401D+00, &
          0.1694895932488235D+00, &
          0.6481099336610571D+00, &
          0.3419306029008594D+00, &
          0.6481099336610571D+00, &
          0.0099594634380836D+00, &
          0.3419306029008594D+00, &
          0.0099594634380836D+00, &
          0.0200903650277176D+00, &
          0.1335404714654308D+00, &
          0.0200903650277176D+00, &
          0.8463691635068515D+00, &
          0.1335404714654308D+00, &
          0.8463691635068515D+00, &
          0.0200903650277176D+00, &
          0.1335404714654308D+00, &
          0.0200903650277176D+00, &
          0.8463691635068515D+00, &
          0.1335404714654308D+00, &
          0.8463691635068515D+00, &
          0.3214917379706315D+00, &
          0.1648190492804087D+00, &
          0.3214917379706315D+00, &
          0.5136892127489598D+00, &
          0.1648190492804087D+00, &
          0.5136892127489598D+00, &
          0.3214917379706315D+00, &
          0.1648190492804087D+00, &
          0.3214917379706315D+00, &
          0.5136892127489598D+00, &
          0.1648190492804087D+00, &
          0.5136892127489598D+00, &
          0.0689956505491457D+00, &
          0.6636404691861656D+00, &
          0.0689956505491457D+00, &
          0.2673638802646887D+00, &
          0.6636404691861656D+00, &
          0.2673638802646887D+00, &
          0.0689956505491457D+00, &
          0.6636404691861656D+00, &
          0.0689956505491457D+00, &
          0.2673638802646887D+00, &
          0.6636404691861656D+00, &
          0.2673638802646887D+00, &
          0.0199626366433414D+00, &
          0.2615529990296567D+00, &
          0.0199626366433414D+00, &
          0.7184843643270019D+00, &
          0.2615529990296567D+00, &
          0.7184843643270019D+00, &
          0.0199626366433414D+00, &
          0.2615529990296567D+00, &
          0.0199626366433414D+00, &
          0.7184843643270019D+00, &
          0.2615529990296567D+00, &
          0.7184843643270019D+00, &
          0.0900297723891655D+00, &
          0.0155218563345725D+00, &
          0.0900297723891655D+00, &
          0.8944483712762620D+00, &
          0.0155218563345725D+00, &
          0.8944483712762620D+00, &
          0.0900297723891655D+00, &
          0.0155218563345725D+00, &
          0.0900297723891655D+00, &
          0.8944483712762620D+00, &
          0.0155218563345725D+00, &
          0.8944483712762620D+00, &
          0.1172214513692467D+00, &
          0.5641016399337317D+00, &
          0.1172214513692467D+00, &
          0.3186769086970216D+00, &
          0.5641016399337317D+00, &
          0.3186769086970216D+00, &
          0.1172214513692467D+00, &
          0.5641016399337317D+00, &
          0.1172214513692467D+00, &
          0.3186769086970216D+00, &
          0.5641016399337317D+00, &
          0.3186769086970216D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.0485337200788840D+00, &
          0.9029325598422321D+00, &
          0.0485337200788840D+00, &
          0.0003340541016131D+00, &
          0.9993318917967738D+00, &
          0.0003340541016131D+00, &
          0.4849074297077183D+00, &
          0.0301851405845634D+00, &
          0.4849074297077183D+00, &
          0.4849074297077183D+00, &
          0.0301851405845634D+00, &
          0.4849074297077183D+00, &
          0.1497566246841529D+00, &
          0.7004867506316943D+00, &
          0.1497566246841529D+00, &
          0.1497566246841529D+00, &
          0.7004867506316943D+00, &
          0.1497566246841529D+00, &
          0.2344536517846724D+00, &
          0.5310926964306552D+00, &
          0.2344536517846724D+00, &
          0.2344536517846724D+00, &
          0.5310926964306552D+00, &
          0.2344536517846724D+00, &
          0.4409030162469282D+00, &
          0.1181939675061436D+00, &
          0.4409030162469282D+00, &
          0.4409030162469282D+00, &
          0.1181939675061436D+00, &
          0.4409030162469282D+00, &
          0.4870412467653055D+00, &
          0.0259175064693890D+00, &
          0.4870412467653055D+00, &
          0.4870412467653055D+00, &
          0.0259175064693890D+00, &
          0.4870412467653055D+00, &
          0.0248638193002129D+00, &
          0.9502723613995742D+00, &
          0.0248638193002129D+00, &
          0.0248638193002129D+00, &
          0.9502723613995742D+00, &
          0.0248638193002129D+00, &
          0.1118542147928236D+00, &
          0.7762915704143528D+00, &
          0.1118542147928236D+00, &
          0.1118542147928236D+00, &
          0.7762915704143528D+00, &
          0.1118542147928236D+00, &
          0.7097236881695401D+00, &
          0.1207867185816364D+00, &
          0.1694895932488235D+00, &
          0.1207867185816364D+00, &
          0.1694895932488235D+00, &
          0.7097236881695401D+00, &
          0.3419306029008594D+00, &
          0.6481099336610571D+00, &
          0.0099594634380836D+00, &
          0.6481099336610571D+00, &
          0.0099594634380836D+00, &
          0.3419306029008594D+00, &
          0.1335404714654308D+00, &
          0.0200903650277176D+00, &
          0.8463691635068515D+00, &
          0.0200903650277176D+00, &
          0.8463691635068515D+00, &
          0.1335404714654308D+00, &
          0.1335404714654308D+00, &
          0.0200903650277176D+00, &
          0.8463691635068515D+00, &
          0.0200903650277176D+00, &
          0.8463691635068515D+00, &
          0.1335404714654308D+00, &
          0.1648190492804087D+00, &
          0.3214917379706315D+00, &
          0.5136892127489598D+00, &
          0.3214917379706315D+00, &
          0.5136892127489598D+00, &
          0.1648190492804087D+00, &
          0.1648190492804087D+00, &
          0.3214917379706315D+00, &
          0.5136892127489598D+00, &
          0.3214917379706315D+00, &
          0.5136892127489598D+00, &
          0.1648190492804087D+00, &
          0.6636404691861656D+00, &
          0.0689956505491457D+00, &
          0.2673638802646887D+00, &
          0.0689956505491457D+00, &
          0.2673638802646887D+00, &
          0.6636404691861656D+00, &
          0.6636404691861656D+00, &
          0.0689956505491457D+00, &
          0.2673638802646887D+00, &
          0.0689956505491457D+00, &
          0.2673638802646887D+00, &
          0.6636404691861656D+00, &
          0.2615529990296567D+00, &
          0.0199626366433414D+00, &
          0.7184843643270019D+00, &
          0.0199626366433414D+00, &
          0.7184843643270019D+00, &
          0.2615529990296567D+00, &
          0.2615529990296567D+00, &
          0.0199626366433414D+00, &
          0.7184843643270019D+00, &
          0.0199626366433414D+00, &
          0.7184843643270019D+00, &
          0.2615529990296567D+00, &
          0.0155218563345725D+00, &
          0.0900297723891655D+00, &
          0.8944483712762620D+00, &
          0.0900297723891655D+00, &
          0.8944483712762620D+00, &
          0.0155218563345725D+00, &
          0.0155218563345725D+00, &
          0.0900297723891655D+00, &
          0.8944483712762620D+00, &
          0.0900297723891655D+00, &
          0.8944483712762620D+00, &
          0.0155218563345725D+00, &
          0.5641016399337317D+00, &
          0.1172214513692467D+00, &
          0.3186769086970216D+00, &
          0.1172214513692467D+00, &
          0.3186769086970216D+00, &
          0.5641016399337317D+00, &
          0.5641016399337317D+00, &
          0.1172214513692467D+00, &
          0.3186769086970216D+00, &
          0.1172214513692467D+00, &
          0.3186769086970216D+00, &
          0.5641016399337317D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.9624802640887508D+00, &
          0.0375197359112492D+00, &
          0.6587282524060957D+00, &
          0.3412717475939043D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.7793516203235177D+00, &
          0.7793516203235177D+00, &
          0.7793516203235177D+00, &
          0.2206483796764823D+00, &
          0.2206483796764823D+00, &
          0.2206483796764823D+00, &
          0.0004872322680031D+00, &
          0.0004872322680031D+00, &
          0.0004872322680031D+00, &
          0.9995127677319969D+00, &
          0.9995127677319969D+00, &
          0.9995127677319969D+00, &
          0.6654124101220354D+00, &
          0.6654124101220354D+00, &
          0.6654124101220354D+00, &
          0.3345875898779646D+00, &
          0.3345875898779646D+00, &
          0.3345875898779646D+00, &
          0.4636840317697777D+00, &
          0.4636840317697777D+00, &
          0.4636840317697777D+00, &
          0.5363159682302223D+00, &
          0.5363159682302223D+00, &
          0.5363159682302223D+00, &
          0.9686710476965903D+00, &
          0.9686710476965903D+00, &
          0.9686710476965903D+00, &
          0.0313289523034098D+00, &
          0.0313289523034098D+00, &
          0.0313289523034098D+00, &
          0.8407183935395043D+00, &
          0.8407183935395043D+00, &
          0.8407183935395043D+00, &
          0.1592816064604957D+00, &
          0.1592816064604957D+00, &
          0.1592816064604957D+00, &
          0.8773254399543938D+00, &
          0.8773254399543938D+00, &
          0.8773254399543938D+00, &
          0.1226745600456062D+00, &
          0.1226745600456062D+00, &
          0.1226745600456062D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.6944689221677809D+00, &
          0.6944689221677809D+00, &
          0.6944689221677809D+00, &
          0.6944689221677809D+00, &
          0.6944689221677809D+00, &
          0.6944689221677809D+00, &
          0.3055310778322191D+00, &
          0.3055310778322191D+00, &
          0.3055310778322191D+00, &
          0.3055310778322191D+00, &
          0.3055310778322191D+00, &
          0.3055310778322191D+00, &
          0.8559006310012134D+00, &
          0.8559006310012134D+00, &
          0.8559006310012134D+00, &
          0.8559006310012134D+00, &
          0.8559006310012134D+00, &
          0.8559006310012134D+00, &
          0.1440993689987866D+00, &
          0.1440993689987866D+00, &
          0.1440993689987866D+00, &
          0.1440993689987866D+00, &
          0.1440993689987866D+00, &
          0.1440993689987866D+00, &
          0.2888455393681256D+00, &
          0.2888455393681256D+00, &
          0.2888455393681256D+00, &
          0.2888455393681256D+00, &
          0.2888455393681256D+00, &
          0.2888455393681256D+00, &
          0.7111544606318744D+00, &
          0.7111544606318744D+00, &
          0.7111544606318744D+00, &
          0.7111544606318744D+00, &
          0.7111544606318744D+00, &
          0.7111544606318744D+00, &
          0.9221446748533075D+00, &
          0.9221446748533075D+00, &
          0.9221446748533075D+00, &
          0.9221446748533075D+00, &
          0.9221446748533075D+00, &
          0.9221446748533075D+00, &
          0.0778553251466925D+00, &
          0.0778553251466925D+00, &
          0.0778553251466925D+00, &
          0.0778553251466925D+00, &
          0.0778553251466925D+00, &
          0.0778553251466925D+00, &
          0.9762657414250897D+00, &
          0.9762657414250897D+00, &
          0.9762657414250897D+00, &
          0.9762657414250897D+00, &
          0.9762657414250897D+00, &
          0.9762657414250897D+00, &
          0.0237342585749103D+00, &
          0.0237342585749103D+00, &
          0.0237342585749103D+00, &
          0.0237342585749103D+00, &
          0.0237342585749103D+00, &
          0.0237342585749103D+00, &
          0.0278439909861404D+00, &
          0.0278439909861404D+00, &
          0.0278439909861404D+00, &
          0.0278439909861404D+00, &
          0.0278439909861404D+00, &
          0.0278439909861404D+00, &
          0.9721560090138596D+00, &
          0.9721560090138596D+00, &
          0.9721560090138596D+00, &
          0.9721560090138596D+00, &
          0.9721560090138596D+00, &
          0.9721560090138596D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          0.0112254009526498D+00, &
          0.0112254009526498D+00, &
          0.0225097870163436D+00, &
          0.0225097870163436D+00, &
          0.0054894240839370D+00, &
          0.0054894240839370D+00, &
          0.0054894240839370D+00, &
          0.0004964406537954D+00, &
          0.0004964406537954D+00, &
          0.0004964406537954D+00, &
          0.0102090400109948D+00, &
          0.0102090400109948D+00, &
          0.0102090400109948D+00, &
          0.0102090400109948D+00, &
          0.0102090400109948D+00, &
          0.0102090400109948D+00, &
          0.0021006495508094D+00, &
          0.0021006495508094D+00, &
          0.0021006495508094D+00, &
          0.0021006495508094D+00, &
          0.0021006495508094D+00, &
          0.0021006495508094D+00, &
          0.0152094206670507D+00, &
          0.0152094206670507D+00, &
          0.0152094206670507D+00, &
          0.0152094206670507D+00, &
          0.0152094206670507D+00, &
          0.0152094206670507D+00, &
          0.0136582892001806D+00, &
          0.0136582892001806D+00, &
          0.0136582892001806D+00, &
          0.0136582892001806D+00, &
          0.0136582892001806D+00, &
          0.0136582892001806D+00, &
          0.0027887393220490D+00, &
          0.0027887393220490D+00, &
          0.0027887393220490D+00, &
          0.0027887393220490D+00, &
          0.0027887393220490D+00, &
          0.0027887393220490D+00, &
          0.0022520950855930D+00, &
          0.0022520950855930D+00, &
          0.0022520950855930D+00, &
          0.0022520950855930D+00, &
          0.0022520950855930D+00, &
          0.0022520950855930D+00, &
          0.0099174482589678D+00, &
          0.0099174482589678D+00, &
          0.0099174482589678D+00, &
          0.0099174482589678D+00, &
          0.0099174482589678D+00, &
          0.0099174482589678D+00, &
          0.0095245941038326D+00, &
          0.0095245941038326D+00, &
          0.0095245941038326D+00, &
          0.0095245941038326D+00, &
          0.0095245941038326D+00, &
          0.0095245941038326D+00, &
          0.0066881354290901D+00, &
          0.0066881354290901D+00, &
          0.0066881354290901D+00, &
          0.0066881354290901D+00, &
          0.0066881354290901D+00, &
          0.0066881354290901D+00, &
          0.0056673224713314D+00, &
          0.0056673224713314D+00, &
          0.0056673224713314D+00, &
          0.0056673224713314D+00, &
          0.0056673224713314D+00, &
          0.0056673224713314D+00, &
          0.0056673224713314D+00, &
          0.0056673224713314D+00, &
          0.0056673224713314D+00, &
          0.0056673224713314D+00, &
          0.0056673224713314D+00, &
          0.0056673224713314D+00, &
          0.0125866454585098D+00, &
          0.0125866454585098D+00, &
          0.0125866454585098D+00, &
          0.0125866454585098D+00, &
          0.0125866454585098D+00, &
          0.0125866454585098D+00, &
          0.0125866454585098D+00, &
          0.0125866454585098D+00, &
          0.0125866454585098D+00, &
          0.0125866454585098D+00, &
          0.0125866454585098D+00, &
          0.0125866454585098D+00, &
          0.0114370118000321D+00, &
          0.0114370118000321D+00, &
          0.0114370118000321D+00, &
          0.0114370118000321D+00, &
          0.0114370118000321D+00, &
          0.0114370118000321D+00, &
          0.0114370118000321D+00, &
          0.0114370118000321D+00, &
          0.0114370118000321D+00, &
          0.0114370118000321D+00, &
          0.0114370118000321D+00, &
          0.0114370118000321D+00, &
          0.0044955975561441D+00, &
          0.0044955975561441D+00, &
          0.0044955975561441D+00, &
          0.0044955975561441D+00, &
          0.0044955975561441D+00, &
          0.0044955975561441D+00, &
          0.0044955975561441D+00, &
          0.0044955975561441D+00, &
          0.0044955975561441D+00, &
          0.0044955975561441D+00, &
          0.0044955975561441D+00, &
          0.0044955975561441D+00, &
          0.0011897740209102D+00, &
          0.0011897740209102D+00, &
          0.0011897740209102D+00, &
          0.0011897740209102D+00, &
          0.0011897740209102D+00, &
          0.0011897740209102D+00, &
          0.0011897740209102D+00, &
          0.0011897740209102D+00, &
          0.0011897740209102D+00, &
          0.0011897740209102D+00, &
          0.0011897740209102D+00, &
          0.0011897740209102D+00, &
          0.0046637786995230D+00, &
          0.0046637786995230D+00, &
          0.0046637786995230D+00, &
          0.0046637786995230D+00, &
          0.0046637786995230D+00, &
          0.0046637786995230D+00, &
          0.0046637786995230D+00, &
          0.0046637786995230D+00, &
          0.0046637786995230D+00, &
          0.0046637786995230D+00, &
          0.0046637786995230D+00, &
          0.0046637786995230D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule13 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule13() returns the prism rule of precision 13.
!
!  Discussion:
!
!    We suppose we are given a triangular prism P with vertices
!    (1,0,0), (0,1,0), (0,0,0),
!    (1,0,1), (0,1,1), (0,0,1),
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over P is 
!    approximated by Q as follows:
!
!    Q = volume(P) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 May 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order symmetric cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 122, Number 1, pages 148-171, 24 August 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real x(n), y(n), z(n): the coordinates of quadrature points.
!
!    real w(n): the quadrature weights, which sum to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: n_save = 162

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
          0.4898439392870481D+00, &
          0.4898439392870481D+00, &
          0.0203121214259039D+00, &
          0.4898439392870481D+00, &
          0.4898439392870481D+00, &
          0.0203121214259039D+00, &
          0.4652402431593082D+00, &
          0.4652402431593082D+00, &
          0.0695195136813836D+00, &
          0.4652402431593082D+00, &
          0.4652402431593082D+00, &
          0.0695195136813836D+00, &
          0.1189161740974999D+00, &
          0.1189161740974999D+00, &
          0.7621676518050002D+00, &
          0.1189161740974999D+00, &
          0.1189161740974999D+00, &
          0.7621676518050002D+00, &
          0.4065656818638698D+00, &
          0.4065656818638698D+00, &
          0.1868686362722604D+00, &
          0.4065656818638698D+00, &
          0.4065656818638698D+00, &
          0.1868686362722604D+00, &
          0.3731243598486834D+00, &
          0.3731243598486834D+00, &
          0.2537512803026333D+00, &
          0.3731243598486834D+00, &
          0.3731243598486834D+00, &
          0.2537512803026333D+00, &
          0.0222257711836755D+00, &
          0.0222257711836755D+00, &
          0.9555484576326490D+00, &
          0.0222257711836755D+00, &
          0.0222257711836755D+00, &
          0.9555484576326490D+00, &
          0.4252710267490136D+00, &
          0.4252710267490136D+00, &
          0.1494579465019728D+00, &
          0.4252710267490136D+00, &
          0.4252710267490136D+00, &
          0.1494579465019728D+00, &
          0.2896929731593648D+00, &
          0.2896929731593648D+00, &
          0.4206140536812703D+00, &
          0.2896929731593648D+00, &
          0.2896929731593648D+00, &
          0.4206140536812703D+00, &
          0.1072868932263340D+00, &
          0.1072868932263340D+00, &
          0.7854262135473320D+00, &
          0.1072868932263340D+00, &
          0.1072868932263340D+00, &
          0.7854262135473320D+00, &
          0.2080125480772502D+00, &
          0.2080125480772502D+00, &
          0.5839749038454997D+00, &
          0.2080125480772502D+00, &
          0.2080125480772502D+00, &
          0.5839749038454997D+00, &
          0.2230030950549982D+00, &
          0.2230030950549982D+00, &
          0.5539938098900037D+00, &
          0.2230030950549982D+00, &
          0.2230030950549982D+00, &
          0.5539938098900037D+00, &
          0.0261743916084927D+00, &
          0.0261743916084927D+00, &
          0.9476512167830147D+00, &
          0.0261743916084927D+00, &
          0.0261743916084927D+00, &
          0.9476512167830147D+00, &
          0.0069141701592508D+00, &
          0.2742659465495736D+00, &
          0.0069141701592508D+00, &
          0.7188198832911755D+00, &
          0.2742659465495736D+00, &
          0.7188198832911755D+00, &
          0.1210652398684976D+00, &
          0.0234237758449634D+00, &
          0.1210652398684976D+00, &
          0.8555109842865389D+00, &
          0.0234237758449634D+00, &
          0.8555109842865389D+00, &
          0.0935360177546358D+00, &
          0.5415308571572887D+00, &
          0.0935360177546358D+00, &
          0.3649331250880756D+00, &
          0.5415308571572887D+00, &
          0.3649331250880756D+00, &
          0.0190809537819619D+00, &
          0.4425451813256520D+00, &
          0.0190809537819619D+00, &
          0.5383738648923861D+00, &
          0.4425451813256520D+00, &
          0.5383738648923861D+00, &
          0.0190809537819619D+00, &
          0.4425451813256520D+00, &
          0.0190809537819619D+00, &
          0.5383738648923861D+00, &
          0.4425451813256520D+00, &
          0.5383738648923861D+00, &
          0.2923396969545124D+00, &
          0.0181776522602859D+00, &
          0.2923396969545124D+00, &
          0.6894826507852017D+00, &
          0.0181776522602859D+00, &
          0.6894826507852017D+00, &
          0.2923396969545124D+00, &
          0.0181776522602859D+00, &
          0.2923396969545124D+00, &
          0.6894826507852017D+00, &
          0.0181776522602859D+00, &
          0.6894826507852017D+00, &
          0.0212309381671597D+00, &
          0.8615902483468726D+00, &
          0.0212309381671597D+00, &
          0.1171788134859677D+00, &
          0.8615902483468726D+00, &
          0.1171788134859677D+00, &
          0.0212309381671597D+00, &
          0.8615902483468726D+00, &
          0.0212309381671597D+00, &
          0.1171788134859677D+00, &
          0.8615902483468726D+00, &
          0.1171788134859677D+00, &
          0.2541452627735283D+00, &
          0.0752550719401233D+00, &
          0.2541452627735283D+00, &
          0.6705996652863484D+00, &
          0.0752550719401233D+00, &
          0.6705996652863484D+00, &
          0.2541452627735283D+00, &
          0.0752550719401233D+00, &
          0.2541452627735283D+00, &
          0.6705996652863484D+00, &
          0.0752550719401233D+00, &
          0.6705996652863484D+00, &
          0.1414034499801619D+00, &
          0.0248530932165772D+00, &
          0.1414034499801619D+00, &
          0.8337434568032609D+00, &
          0.0248530932165772D+00, &
          0.8337434568032609D+00, &
          0.1414034499801619D+00, &
          0.0248530932165772D+00, &
          0.1414034499801619D+00, &
          0.8337434568032609D+00, &
          0.0248530932165772D+00, &
          0.8337434568032609D+00, &
          0.5920519581168684D+00, &
          0.1079442302776815D+00, &
          0.5920519581168684D+00, &
          0.3000038116054501D+00, &
          0.1079442302776815D+00, &
          0.3000038116054501D+00, &
          0.5920519581168684D+00, &
          0.1079442302776815D+00, &
          0.5920519581168684D+00, &
          0.3000038116054501D+00, &
          0.1079442302776815D+00, &
          0.3000038116054501D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.4898439392870481D+00, &
          0.0203121214259039D+00, &
          0.4898439392870481D+00, &
          0.4898439392870481D+00, &
          0.0203121214259039D+00, &
          0.4898439392870481D+00, &
          0.4652402431593082D+00, &
          0.0695195136813836D+00, &
          0.4652402431593082D+00, &
          0.4652402431593082D+00, &
          0.0695195136813836D+00, &
          0.4652402431593082D+00, &
          0.1189161740974999D+00, &
          0.7621676518050002D+00, &
          0.1189161740974999D+00, &
          0.1189161740974999D+00, &
          0.7621676518050002D+00, &
          0.1189161740974999D+00, &
          0.4065656818638698D+00, &
          0.1868686362722604D+00, &
          0.4065656818638698D+00, &
          0.4065656818638698D+00, &
          0.1868686362722604D+00, &
          0.4065656818638698D+00, &
          0.3731243598486834D+00, &
          0.2537512803026333D+00, &
          0.3731243598486834D+00, &
          0.3731243598486834D+00, &
          0.2537512803026333D+00, &
          0.3731243598486834D+00, &
          0.0222257711836755D+00, &
          0.9555484576326490D+00, &
          0.0222257711836755D+00, &
          0.0222257711836755D+00, &
          0.9555484576326490D+00, &
          0.0222257711836755D+00, &
          0.4252710267490136D+00, &
          0.1494579465019728D+00, &
          0.4252710267490136D+00, &
          0.4252710267490136D+00, &
          0.1494579465019728D+00, &
          0.4252710267490136D+00, &
          0.2896929731593648D+00, &
          0.4206140536812703D+00, &
          0.2896929731593648D+00, &
          0.2896929731593648D+00, &
          0.4206140536812703D+00, &
          0.2896929731593648D+00, &
          0.1072868932263340D+00, &
          0.7854262135473320D+00, &
          0.1072868932263340D+00, &
          0.1072868932263340D+00, &
          0.7854262135473320D+00, &
          0.1072868932263340D+00, &
          0.2080125480772502D+00, &
          0.5839749038454997D+00, &
          0.2080125480772502D+00, &
          0.2080125480772502D+00, &
          0.5839749038454997D+00, &
          0.2080125480772502D+00, &
          0.2230030950549982D+00, &
          0.5539938098900037D+00, &
          0.2230030950549982D+00, &
          0.2230030950549982D+00, &
          0.5539938098900037D+00, &
          0.2230030950549982D+00, &
          0.0261743916084927D+00, &
          0.9476512167830147D+00, &
          0.0261743916084927D+00, &
          0.0261743916084927D+00, &
          0.9476512167830147D+00, &
          0.0261743916084927D+00, &
          0.2742659465495736D+00, &
          0.0069141701592508D+00, &
          0.7188198832911755D+00, &
          0.0069141701592508D+00, &
          0.7188198832911755D+00, &
          0.2742659465495736D+00, &
          0.0234237758449634D+00, &
          0.1210652398684976D+00, &
          0.8555109842865389D+00, &
          0.1210652398684976D+00, &
          0.8555109842865389D+00, &
          0.0234237758449634D+00, &
          0.5415308571572887D+00, &
          0.0935360177546358D+00, &
          0.3649331250880756D+00, &
          0.0935360177546358D+00, &
          0.3649331250880756D+00, &
          0.5415308571572887D+00, &
          0.4425451813256520D+00, &
          0.0190809537819619D+00, &
          0.5383738648923861D+00, &
          0.0190809537819619D+00, &
          0.5383738648923861D+00, &
          0.4425451813256520D+00, &
          0.4425451813256520D+00, &
          0.0190809537819619D+00, &
          0.5383738648923861D+00, &
          0.0190809537819619D+00, &
          0.5383738648923861D+00, &
          0.4425451813256520D+00, &
          0.0181776522602859D+00, &
          0.2923396969545124D+00, &
          0.6894826507852017D+00, &
          0.2923396969545124D+00, &
          0.6894826507852017D+00, &
          0.0181776522602859D+00, &
          0.0181776522602859D+00, &
          0.2923396969545124D+00, &
          0.6894826507852017D+00, &
          0.2923396969545124D+00, &
          0.6894826507852017D+00, &
          0.0181776522602859D+00, &
          0.8615902483468726D+00, &
          0.0212309381671597D+00, &
          0.1171788134859677D+00, &
          0.0212309381671597D+00, &
          0.1171788134859677D+00, &
          0.8615902483468726D+00, &
          0.8615902483468726D+00, &
          0.0212309381671597D+00, &
          0.1171788134859677D+00, &
          0.0212309381671597D+00, &
          0.1171788134859677D+00, &
          0.8615902483468726D+00, &
          0.0752550719401233D+00, &
          0.2541452627735283D+00, &
          0.6705996652863484D+00, &
          0.2541452627735283D+00, &
          0.6705996652863484D+00, &
          0.0752550719401233D+00, &
          0.0752550719401233D+00, &
          0.2541452627735283D+00, &
          0.6705996652863484D+00, &
          0.2541452627735283D+00, &
          0.6705996652863484D+00, &
          0.0752550719401233D+00, &
          0.0248530932165772D+00, &
          0.1414034499801619D+00, &
          0.8337434568032609D+00, &
          0.1414034499801619D+00, &
          0.8337434568032609D+00, &
          0.0248530932165772D+00, &
          0.0248530932165772D+00, &
          0.1414034499801619D+00, &
          0.8337434568032609D+00, &
          0.1414034499801619D+00, &
          0.8337434568032609D+00, &
          0.0248530932165772D+00, &
          0.1079442302776815D+00, &
          0.5920519581168684D+00, &
          0.3000038116054501D+00, &
          0.5920519581168684D+00, &
          0.3000038116054501D+00, &
          0.1079442302776815D+00, &
          0.1079442302776815D+00, &
          0.5920519581168684D+00, &
          0.3000038116054501D+00, &
          0.5920519581168684D+00, &
          0.3000038116054501D+00, &
          0.1079442302776815D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.9628749136752521D+00, &
          0.9628749136752521D+00, &
          0.9628749136752521D+00, &
          0.0371250863247479D+00, &
          0.0371250863247479D+00, &
          0.0371250863247479D+00, &
          0.1596095438845206D+00, &
          0.1596095438845206D+00, &
          0.1596095438845206D+00, &
          0.8403904561154794D+00, &
          0.8403904561154794D+00, &
          0.8403904561154794D+00, &
          0.9056527042389668D+00, &
          0.9056527042389668D+00, &
          0.9056527042389668D+00, &
          0.0943472957610332D+00, &
          0.0943472957610332D+00, &
          0.0943472957610332D+00, &
          0.5096252745686303D+00, &
          0.5096252745686303D+00, &
          0.5096252745686303D+00, &
          0.4903747254313698D+00, &
          0.4903747254313698D+00, &
          0.4903747254313698D+00, &
          0.9622105990610692D+00, &
          0.9622105990610692D+00, &
          0.9622105990610692D+00, &
          0.0377894009389309D+00, &
          0.0377894009389309D+00, &
          0.0377894009389309D+00, &
          0.6451848165106220D+00, &
          0.6451848165106220D+00, &
          0.6451848165106220D+00, &
          0.3548151834893781D+00, &
          0.3548151834893781D+00, &
          0.3548151834893781D+00, &
          0.8358323923694781D+00, &
          0.8358323923694781D+00, &
          0.8358323923694781D+00, &
          0.1641676076305219D+00, &
          0.1641676076305219D+00, &
          0.1641676076305219D+00, &
          0.6970082735594627D+00, &
          0.6970082735594627D+00, &
          0.6970082735594627D+00, &
          0.3029917264405373D+00, &
          0.3029917264405373D+00, &
          0.3029917264405373D+00, &
          0.5806209783551364D+00, &
          0.5806209783551364D+00, &
          0.5806209783551364D+00, &
          0.4193790216448636D+00, &
          0.4193790216448636D+00, &
          0.4193790216448636D+00, &
          0.6310517191712730D+00, &
          0.6310517191712730D+00, &
          0.6310517191712730D+00, &
          0.3689482808287270D+00, &
          0.3689482808287270D+00, &
          0.3689482808287270D+00, &
          0.8758646517044397D+00, &
          0.8758646517044397D+00, &
          0.8758646517044397D+00, &
          0.1241353482955603D+00, &
          0.1241353482955603D+00, &
          0.1241353482955603D+00, &
          0.9265272967069804D+00, &
          0.9265272967069804D+00, &
          0.9265272967069804D+00, &
          0.0734727032930196D+00, &
          0.0734727032930196D+00, &
          0.0734727032930196D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.6831055612129073D+00, &
          0.6831055612129073D+00, &
          0.6831055612129073D+00, &
          0.6831055612129073D+00, &
          0.6831055612129073D+00, &
          0.6831055612129073D+00, &
          0.3168944387870927D+00, &
          0.3168944387870927D+00, &
          0.3168944387870927D+00, &
          0.3168944387870927D+00, &
          0.3168944387870927D+00, &
          0.3168944387870927D+00, &
          0.8946991133344783D+00, &
          0.8946991133344783D+00, &
          0.8946991133344783D+00, &
          0.8946991133344783D+00, &
          0.8946991133344783D+00, &
          0.8946991133344783D+00, &
          0.1053008866655217D+00, &
          0.1053008866655217D+00, &
          0.1053008866655217D+00, &
          0.1053008866655217D+00, &
          0.1053008866655217D+00, &
          0.1053008866655217D+00, &
          0.7783161730088111D+00, &
          0.7783161730088111D+00, &
          0.7783161730088111D+00, &
          0.7783161730088111D+00, &
          0.7783161730088111D+00, &
          0.7783161730088111D+00, &
          0.2216838269911889D+00, &
          0.2216838269911889D+00, &
          0.2216838269911889D+00, &
          0.2216838269911889D+00, &
          0.2216838269911889D+00, &
          0.2216838269911889D+00, &
          0.2838393063545257D+00, &
          0.2838393063545257D+00, &
          0.2838393063545257D+00, &
          0.2838393063545257D+00, &
          0.2838393063545257D+00, &
          0.2838393063545257D+00, &
          0.7161606936454743D+00, &
          0.7161606936454743D+00, &
          0.7161606936454743D+00, &
          0.7161606936454743D+00, &
          0.7161606936454743D+00, &
          0.7161606936454743D+00, &
          0.9877489005025584D+00, &
          0.9877489005025584D+00, &
          0.9877489005025584D+00, &
          0.9877489005025584D+00, &
          0.9877489005025584D+00, &
          0.9877489005025584D+00, &
          0.0122510994974416D+00, &
          0.0122510994974416D+00, &
          0.0122510994974416D+00, &
          0.0122510994974416D+00, &
          0.0122510994974416D+00, &
          0.0122510994974416D+00, &
          0.9801022671238835D+00, &
          0.9801022671238835D+00, &
          0.9801022671238835D+00, &
          0.9801022671238835D+00, &
          0.9801022671238835D+00, &
          0.9801022671238835D+00, &
          0.0198977328761166D+00, &
          0.0198977328761166D+00, &
          0.0198977328761166D+00, &
          0.0198977328761166D+00, &
          0.0198977328761166D+00, &
          0.0198977328761166D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          0.0026290597067936D+00, &
          0.0026290597067936D+00, &
          0.0026290597067936D+00, &
          0.0026290597067936D+00, &
          0.0026290597067936D+00, &
          0.0026290597067936D+00, &
          0.0059734237493126D+00, &
          0.0059734237493126D+00, &
          0.0059734237493126D+00, &
          0.0059734237493126D+00, &
          0.0059734237493126D+00, &
          0.0059734237493126D+00, &
          0.0080935494241557D+00, &
          0.0080935494241557D+00, &
          0.0080935494241557D+00, &
          0.0080935494241557D+00, &
          0.0080935494241557D+00, &
          0.0080935494241557D+00, &
          0.0068191162590357D+00, &
          0.0068191162590357D+00, &
          0.0068191162590357D+00, &
          0.0068191162590357D+00, &
          0.0068191162590357D+00, &
          0.0068191162590357D+00, &
          0.0049334513687730D+00, &
          0.0049334513687730D+00, &
          0.0049334513687730D+00, &
          0.0049334513687730D+00, &
          0.0049334513687730D+00, &
          0.0049334513687730D+00, &
          0.0020147080132884D+00, &
          0.0020147080132884D+00, &
          0.0020147080132884D+00, &
          0.0020147080132884D+00, &
          0.0020147080132884D+00, &
          0.0020147080132884D+00, &
          0.0123357707370753D+00, &
          0.0123357707370753D+00, &
          0.0123357707370753D+00, &
          0.0123357707370753D+00, &
          0.0123357707370753D+00, &
          0.0123357707370753D+00, &
          0.0113647943376711D+00, &
          0.0113647943376711D+00, &
          0.0113647943376711D+00, &
          0.0113647943376711D+00, &
          0.0113647943376711D+00, &
          0.0113647943376711D+00, &
          0.0066668649606825D+00, &
          0.0066668649606825D+00, &
          0.0066668649606825D+00, &
          0.0066668649606825D+00, &
          0.0066668649606825D+00, &
          0.0066668649606825D+00, &
          0.0119452201363976D+00, &
          0.0119452201363976D+00, &
          0.0119452201363976D+00, &
          0.0119452201363976D+00, &
          0.0119452201363976D+00, &
          0.0119452201363976D+00, &
          0.0106878809406112D+00, &
          0.0106878809406112D+00, &
          0.0106878809406112D+00, &
          0.0106878809406112D+00, &
          0.0106878809406112D+00, &
          0.0106878809406112D+00, &
          0.0016663443274632D+00, &
          0.0016663443274632D+00, &
          0.0016663443274632D+00, &
          0.0016663443274632D+00, &
          0.0016663443274632D+00, &
          0.0016663443274632D+00, &
          0.0045120516170652D+00, &
          0.0045120516170652D+00, &
          0.0045120516170652D+00, &
          0.0045120516170652D+00, &
          0.0045120516170652D+00, &
          0.0045120516170652D+00, &
          0.0034887068394433D+00, &
          0.0034887068394433D+00, &
          0.0034887068394433D+00, &
          0.0034887068394433D+00, &
          0.0034887068394433D+00, &
          0.0034887068394433D+00, &
          0.0111803506944158D+00, &
          0.0111803506944158D+00, &
          0.0111803506944158D+00, &
          0.0111803506944158D+00, &
          0.0111803506944158D+00, &
          0.0111803506944158D+00, &
          0.0042052602324505D+00, &
          0.0042052602324505D+00, &
          0.0042052602324505D+00, &
          0.0042052602324505D+00, &
          0.0042052602324505D+00, &
          0.0042052602324505D+00, &
          0.0042052602324505D+00, &
          0.0042052602324505D+00, &
          0.0042052602324505D+00, &
          0.0042052602324505D+00, &
          0.0042052602324505D+00, &
          0.0042052602324505D+00, &
          0.0043767917497427D+00, &
          0.0043767917497427D+00, &
          0.0043767917497427D+00, &
          0.0043767917497427D+00, &
          0.0043767917497427D+00, &
          0.0043767917497427D+00, &
          0.0043767917497427D+00, &
          0.0043767917497427D+00, &
          0.0043767917497427D+00, &
          0.0043767917497427D+00, &
          0.0043767917497427D+00, &
          0.0043767917497427D+00, &
          0.0044067550419370D+00, &
          0.0044067550419370D+00, &
          0.0044067550419370D+00, &
          0.0044067550419370D+00, &
          0.0044067550419370D+00, &
          0.0044067550419370D+00, &
          0.0044067550419370D+00, &
          0.0044067550419370D+00, &
          0.0044067550419370D+00, &
          0.0044067550419370D+00, &
          0.0044067550419370D+00, &
          0.0044067550419370D+00, &
          0.0120350326772681D+00, &
          0.0120350326772681D+00, &
          0.0120350326772681D+00, &
          0.0120350326772681D+00, &
          0.0120350326772681D+00, &
          0.0120350326772681D+00, &
          0.0120350326772681D+00, &
          0.0120350326772681D+00, &
          0.0120350326772681D+00, &
          0.0120350326772681D+00, &
          0.0120350326772681D+00, &
          0.0120350326772681D+00, &
          0.0013986921884649D+00, &
          0.0013986921884649D+00, &
          0.0013986921884649D+00, &
          0.0013986921884649D+00, &
          0.0013986921884649D+00, &
          0.0013986921884649D+00, &
          0.0013986921884649D+00, &
          0.0013986921884649D+00, &
          0.0013986921884649D+00, &
          0.0013986921884649D+00, &
          0.0013986921884649D+00, &
          0.0013986921884649D+00, &
          0.0047551548873781D+00, &
          0.0047551548873781D+00, &
          0.0047551548873781D+00, &
          0.0047551548873781D+00, &
          0.0047551548873781D+00, &
          0.0047551548873781D+00, &
          0.0047551548873781D+00, &
          0.0047551548873781D+00, &
          0.0047551548873781D+00, &
          0.0047551548873781D+00, &
          0.0047551548873781D+00, &
          0.0047551548873781D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule14 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule14() returns the prism rule of precision 14.
!
!  Discussion:
!
!    We suppose we are given a triangular prism P with vertices
!    (1,0,0), (0,1,0), (0,0,0),
!    (1,0,1), (0,1,1), (0,0,1),
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over P is 
!    approximated by Q as follows:
!
!    Q = volume(P) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 May 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order symmetric cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 122, Number 1, pages 148-171, 24 August 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real x(n), y(n), z(n): the coordinates of quadrature points.
!
!    real w(n): the quadrature weights, which sum to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: n_save = 194

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.1531409494078462D+00, &
          0.1531409494078462D+00, &
          0.6937181011843075D+00, &
          0.1094642266754876D+00, &
          0.1094642266754876D+00, &
          0.7810715466490248D+00, &
          0.4852139749601931D+00, &
          0.4852139749601931D+00, &
          0.0295720500796137D+00, &
          0.4852139749601931D+00, &
          0.4852139749601931D+00, &
          0.0295720500796137D+00, &
          0.4811857851081537D+00, &
          0.4811857851081537D+00, &
          0.0376284297836926D+00, &
          0.4811857851081537D+00, &
          0.4811857851081537D+00, &
          0.0376284297836926D+00, &
          0.0854252549169795D+00, &
          0.0854252549169795D+00, &
          0.8291494901660410D+00, &
          0.0854252549169795D+00, &
          0.0854252549169795D+00, &
          0.8291494901660410D+00, &
          0.4949498254464672D+00, &
          0.4949498254464672D+00, &
          0.0101003491070657D+00, &
          0.4949498254464672D+00, &
          0.4949498254464672D+00, &
          0.0101003491070657D+00, &
          0.2006722689888837D+00, &
          0.2006722689888837D+00, &
          0.5986554620222326D+00, &
          0.2006722689888837D+00, &
          0.2006722689888837D+00, &
          0.5986554620222326D+00, &
          0.0134631402389198D+00, &
          0.0134631402389198D+00, &
          0.9730737195221603D+00, &
          0.0134631402389198D+00, &
          0.0134631402389198D+00, &
          0.9730737195221603D+00, &
          0.3888073656098609D+00, &
          0.3888073656098609D+00, &
          0.2223852687802782D+00, &
          0.3888073656098609D+00, &
          0.3888073656098609D+00, &
          0.2223852687802782D+00, &
          0.2589597694248802D+00, &
          0.2589597694248802D+00, &
          0.4820804611502396D+00, &
          0.2589597694248802D+00, &
          0.2589597694248802D+00, &
          0.4820804611502396D+00, &
          0.0656674639158517D+00, &
          0.0656674639158517D+00, &
          0.8686650721682967D+00, &
          0.0656674639158517D+00, &
          0.0656674639158517D+00, &
          0.8686650721682967D+00, &
          0.4462504387663027D+00, &
          0.4462504387663027D+00, &
          0.1074991224673947D+00, &
          0.4462504387663027D+00, &
          0.4462504387663027D+00, &
          0.1074991224673947D+00, &
          0.1155506883468244D+00, &
          0.1155506883468244D+00, &
          0.7688986233063512D+00, &
          0.1155506883468244D+00, &
          0.1155506883468244D+00, &
          0.7688986233063512D+00, &
          0.0257939798379175D+00, &
          0.0257939798379175D+00, &
          0.9484120403241649D+00, &
          0.0257939798379175D+00, &
          0.0257939798379175D+00, &
          0.9484120403241649D+00, &
          0.0231313365079955D+00, &
          0.3737411598567404D+00, &
          0.0231313365079955D+00, &
          0.6031275036352640D+00, &
          0.3737411598567404D+00, &
          0.6031275036352640D+00, &
          0.0780837656724579D+00, &
          0.0087167837437324D+00, &
          0.0780837656724579D+00, &
          0.9131994505838098D+00, &
          0.0087167837437324D+00, &
          0.9131994505838098D+00, &
          0.1932090397525293D+00, &
          0.4408189862889988D+00, &
          0.1932090397525293D+00, &
          0.3659719739584719D+00, &
          0.4408189862889988D+00, &
          0.3659719739584719D+00, &
          0.0110695610243579D+00, &
          0.3495104935335414D+00, &
          0.0110695610243579D+00, &
          0.6394199454421007D+00, &
          0.3495104935335414D+00, &
          0.6394199454421007D+00, &
          0.0110695610243579D+00, &
          0.3495104935335414D+00, &
          0.0110695610243579D+00, &
          0.6394199454421007D+00, &
          0.3495104935335414D+00, &
          0.6394199454421007D+00, &
          0.2813786819004687D+00, &
          0.0725240717490409D+00, &
          0.2813786819004687D+00, &
          0.6460972463504904D+00, &
          0.0725240717490409D+00, &
          0.6460972463504904D+00, &
          0.2813786819004687D+00, &
          0.0725240717490409D+00, &
          0.2813786819004687D+00, &
          0.6460972463504904D+00, &
          0.0725240717490409D+00, &
          0.6460972463504904D+00, &
          0.0136067459730556D+00, &
          0.7694543610610796D+00, &
          0.0136067459730556D+00, &
          0.2169388929658648D+00, &
          0.7694543610610796D+00, &
          0.2169388929658648D+00, &
          0.0136067459730556D+00, &
          0.7694543610610796D+00, &
          0.0136067459730556D+00, &
          0.2169388929658648D+00, &
          0.7694543610610796D+00, &
          0.2169388929658648D+00, &
          0.3160642637515643D+00, &
          0.0980612921350026D+00, &
          0.3160642637515643D+00, &
          0.5858744441134331D+00, &
          0.0980612921350026D+00, &
          0.5858744441134331D+00, &
          0.3160642637515643D+00, &
          0.0980612921350026D+00, &
          0.3160642637515643D+00, &
          0.5858744441134331D+00, &
          0.0980612921350026D+00, &
          0.5858744441134331D+00, &
          0.2010219519914552D+00, &
          0.0101533805746944D+00, &
          0.2010219519914552D+00, &
          0.7888246674338504D+00, &
          0.0101533805746944D+00, &
          0.7888246674338504D+00, &
          0.2010219519914552D+00, &
          0.0101533805746944D+00, &
          0.2010219519914552D+00, &
          0.7888246674338504D+00, &
          0.0101533805746944D+00, &
          0.7888246674338504D+00, &
          0.5404862370418329D+00, &
          0.1675040220811526D+00, &
          0.5404862370418329D+00, &
          0.2920097408770145D+00, &
          0.1675040220811526D+00, &
          0.2920097408770145D+00, &
          0.5404862370418329D+00, &
          0.1675040220811526D+00, &
          0.5404862370418329D+00, &
          0.2920097408770145D+00, &
          0.1675040220811526D+00, &
          0.2920097408770145D+00, &
          0.1016331118955575D+00, &
          0.0168250918958854D+00, &
          0.1016331118955575D+00, &
          0.8815417962085570D+00, &
          0.0168250918958854D+00, &
          0.8815417962085570D+00, &
          0.1016331118955575D+00, &
          0.0168250918958854D+00, &
          0.1016331118955575D+00, &
          0.8815417962085570D+00, &
          0.0168250918958854D+00, &
          0.8815417962085570D+00, &
          0.7553391402298588D+00, &
          0.0552850680433709D+00, &
          0.7553391402298588D+00, &
          0.1893757917267703D+00, &
          0.0552850680433709D+00, &
          0.1893757917267703D+00, &
          0.7553391402298588D+00, &
          0.0552850680433709D+00, &
          0.7553391402298588D+00, &
          0.1893757917267703D+00, &
          0.0552850680433709D+00, &
          0.1893757917267703D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.1531409494078462D+00, &
          0.6937181011843075D+00, &
          0.1531409494078462D+00, &
          0.1094642266754876D+00, &
          0.7810715466490248D+00, &
          0.1094642266754876D+00, &
          0.4852139749601931D+00, &
          0.0295720500796137D+00, &
          0.4852139749601931D+00, &
          0.4852139749601931D+00, &
          0.0295720500796137D+00, &
          0.4852139749601931D+00, &
          0.4811857851081537D+00, &
          0.0376284297836926D+00, &
          0.4811857851081537D+00, &
          0.4811857851081537D+00, &
          0.0376284297836926D+00, &
          0.4811857851081537D+00, &
          0.0854252549169795D+00, &
          0.8291494901660410D+00, &
          0.0854252549169795D+00, &
          0.0854252549169795D+00, &
          0.8291494901660410D+00, &
          0.0854252549169795D+00, &
          0.4949498254464672D+00, &
          0.0101003491070657D+00, &
          0.4949498254464672D+00, &
          0.4949498254464672D+00, &
          0.0101003491070657D+00, &
          0.4949498254464672D+00, &
          0.2006722689888837D+00, &
          0.5986554620222326D+00, &
          0.2006722689888837D+00, &
          0.2006722689888837D+00, &
          0.5986554620222326D+00, &
          0.2006722689888837D+00, &
          0.0134631402389198D+00, &
          0.9730737195221603D+00, &
          0.0134631402389198D+00, &
          0.0134631402389198D+00, &
          0.9730737195221603D+00, &
          0.0134631402389198D+00, &
          0.3888073656098609D+00, &
          0.2223852687802782D+00, &
          0.3888073656098609D+00, &
          0.3888073656098609D+00, &
          0.2223852687802782D+00, &
          0.3888073656098609D+00, &
          0.2589597694248802D+00, &
          0.4820804611502396D+00, &
          0.2589597694248802D+00, &
          0.2589597694248802D+00, &
          0.4820804611502396D+00, &
          0.2589597694248802D+00, &
          0.0656674639158517D+00, &
          0.8686650721682967D+00, &
          0.0656674639158517D+00, &
          0.0656674639158517D+00, &
          0.8686650721682967D+00, &
          0.0656674639158517D+00, &
          0.4462504387663027D+00, &
          0.1074991224673947D+00, &
          0.4462504387663027D+00, &
          0.4462504387663027D+00, &
          0.1074991224673947D+00, &
          0.4462504387663027D+00, &
          0.1155506883468244D+00, &
          0.7688986233063512D+00, &
          0.1155506883468244D+00, &
          0.1155506883468244D+00, &
          0.7688986233063512D+00, &
          0.1155506883468244D+00, &
          0.0257939798379175D+00, &
          0.9484120403241649D+00, &
          0.0257939798379175D+00, &
          0.0257939798379175D+00, &
          0.9484120403241649D+00, &
          0.0257939798379175D+00, &
          0.3737411598567404D+00, &
          0.0231313365079955D+00, &
          0.6031275036352640D+00, &
          0.0231313365079955D+00, &
          0.6031275036352640D+00, &
          0.3737411598567404D+00, &
          0.0087167837437324D+00, &
          0.0780837656724579D+00, &
          0.9131994505838098D+00, &
          0.0780837656724579D+00, &
          0.9131994505838098D+00, &
          0.0087167837437324D+00, &
          0.4408189862889988D+00, &
          0.1932090397525293D+00, &
          0.3659719739584719D+00, &
          0.1932090397525293D+00, &
          0.3659719739584719D+00, &
          0.4408189862889988D+00, &
          0.3495104935335414D+00, &
          0.0110695610243579D+00, &
          0.6394199454421007D+00, &
          0.0110695610243579D+00, &
          0.6394199454421007D+00, &
          0.3495104935335414D+00, &
          0.3495104935335414D+00, &
          0.0110695610243579D+00, &
          0.6394199454421007D+00, &
          0.0110695610243579D+00, &
          0.6394199454421007D+00, &
          0.3495104935335414D+00, &
          0.0725240717490409D+00, &
          0.2813786819004687D+00, &
          0.6460972463504904D+00, &
          0.2813786819004687D+00, &
          0.6460972463504904D+00, &
          0.0725240717490409D+00, &
          0.0725240717490409D+00, &
          0.2813786819004687D+00, &
          0.6460972463504904D+00, &
          0.2813786819004687D+00, &
          0.6460972463504904D+00, &
          0.0725240717490409D+00, &
          0.7694543610610796D+00, &
          0.0136067459730556D+00, &
          0.2169388929658648D+00, &
          0.0136067459730556D+00, &
          0.2169388929658648D+00, &
          0.7694543610610796D+00, &
          0.7694543610610796D+00, &
          0.0136067459730556D+00, &
          0.2169388929658648D+00, &
          0.0136067459730556D+00, &
          0.2169388929658648D+00, &
          0.7694543610610796D+00, &
          0.0980612921350026D+00, &
          0.3160642637515643D+00, &
          0.5858744441134331D+00, &
          0.3160642637515643D+00, &
          0.5858744441134331D+00, &
          0.0980612921350026D+00, &
          0.0980612921350026D+00, &
          0.3160642637515643D+00, &
          0.5858744441134331D+00, &
          0.3160642637515643D+00, &
          0.5858744441134331D+00, &
          0.0980612921350026D+00, &
          0.0101533805746944D+00, &
          0.2010219519914552D+00, &
          0.7888246674338504D+00, &
          0.2010219519914552D+00, &
          0.7888246674338504D+00, &
          0.0101533805746944D+00, &
          0.0101533805746944D+00, &
          0.2010219519914552D+00, &
          0.7888246674338504D+00, &
          0.2010219519914552D+00, &
          0.7888246674338504D+00, &
          0.0101533805746944D+00, &
          0.1675040220811526D+00, &
          0.5404862370418329D+00, &
          0.2920097408770145D+00, &
          0.5404862370418329D+00, &
          0.2920097408770145D+00, &
          0.1675040220811526D+00, &
          0.1675040220811526D+00, &
          0.5404862370418329D+00, &
          0.2920097408770145D+00, &
          0.5404862370418329D+00, &
          0.2920097408770145D+00, &
          0.1675040220811526D+00, &
          0.0168250918958854D+00, &
          0.1016331118955575D+00, &
          0.8815417962085570D+00, &
          0.1016331118955575D+00, &
          0.8815417962085570D+00, &
          0.0168250918958854D+00, &
          0.0168250918958854D+00, &
          0.1016331118955575D+00, &
          0.8815417962085570D+00, &
          0.1016331118955575D+00, &
          0.8815417962085570D+00, &
          0.0168250918958854D+00, &
          0.0552850680433709D+00, &
          0.7553391402298588D+00, &
          0.1893757917267703D+00, &
          0.7553391402298588D+00, &
          0.1893757917267703D+00, &
          0.0552850680433709D+00, &
          0.0552850680433709D+00, &
          0.7553391402298588D+00, &
          0.1893757917267703D+00, &
          0.7553391402298588D+00, &
          0.1893757917267703D+00, &
          0.0552850680433709D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.2284245901167349D+00, &
          0.7715754098832651D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9749706927792001D+00, &
          0.9749706927792001D+00, &
          0.9749706927792001D+00, &
          0.0250293072207999D+00, &
          0.0250293072207999D+00, &
          0.0250293072207999D+00, &
          0.3291597905587964D+00, &
          0.3291597905587964D+00, &
          0.3291597905587964D+00, &
          0.6708402094412036D+00, &
          0.6708402094412036D+00, &
          0.6708402094412036D+00, &
          0.9996712623901741D+00, &
          0.9996712623901741D+00, &
          0.9996712623901741D+00, &
          0.0003287376098259D+00, &
          0.0003287376098259D+00, &
          0.0003287376098259D+00, &
          0.7021882096439400D+00, &
          0.7021882096439400D+00, &
          0.7021882096439400D+00, &
          0.2978117903560600D+00, &
          0.2978117903560600D+00, &
          0.2978117903560600D+00, &
          0.8352161845608261D+00, &
          0.8352161845608261D+00, &
          0.8352161845608261D+00, &
          0.1647838154391739D+00, &
          0.1647838154391739D+00, &
          0.1647838154391739D+00, &
          0.7201109571464277D+00, &
          0.7201109571464277D+00, &
          0.7201109571464277D+00, &
          0.2798890428535724D+00, &
          0.2798890428535724D+00, &
          0.2798890428535724D+00, &
          0.9332471062622827D+00, &
          0.9332471062622827D+00, &
          0.9332471062622827D+00, &
          0.0667528937377173D+00, &
          0.0667528937377173D+00, &
          0.0667528937377173D+00, &
          0.3097333917149071D+00, &
          0.3097333917149071D+00, &
          0.3097333917149071D+00, &
          0.6902666082850929D+00, &
          0.6902666082850929D+00, &
          0.6902666082850929D+00, &
          0.6590250213320482D+00, &
          0.6590250213320482D+00, &
          0.6590250213320482D+00, &
          0.3409749786679518D+00, &
          0.3409749786679518D+00, &
          0.3409749786679518D+00, &
          0.8368689390766315D+00, &
          0.8368689390766315D+00, &
          0.8368689390766315D+00, &
          0.1631310609233684D+00, &
          0.1631310609233684D+00, &
          0.1631310609233684D+00, &
          0.9263853262636454D+00, &
          0.9263853262636454D+00, &
          0.9263853262636454D+00, &
          0.0736146737363546D+00, &
          0.0736146737363546D+00, &
          0.0736146737363546D+00, &
          0.9586297434779142D+00, &
          0.9586297434779142D+00, &
          0.9586297434779142D+00, &
          0.0413702565220858D+00, &
          0.0413702565220858D+00, &
          0.0413702565220858D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.8444056700436123D+00, &
          0.8444056700436123D+00, &
          0.8444056700436123D+00, &
          0.8444056700436123D+00, &
          0.8444056700436123D+00, &
          0.8444056700436123D+00, &
          0.1555943299563877D+00, &
          0.1555943299563877D+00, &
          0.1555943299563877D+00, &
          0.1555943299563877D+00, &
          0.1555943299563877D+00, &
          0.1555943299563877D+00, &
          0.9409140558660263D+00, &
          0.9409140558660263D+00, &
          0.9409140558660263D+00, &
          0.9409140558660263D+00, &
          0.9409140558660263D+00, &
          0.9409140558660263D+00, &
          0.0590859441339738D+00, &
          0.0590859441339738D+00, &
          0.0590859441339738D+00, &
          0.0590859441339738D+00, &
          0.0590859441339738D+00, &
          0.0590859441339738D+00, &
          0.6041114898852048D+00, &
          0.6041114898852048D+00, &
          0.6041114898852048D+00, &
          0.6041114898852048D+00, &
          0.6041114898852048D+00, &
          0.6041114898852048D+00, &
          0.3958885101147953D+00, &
          0.3958885101147953D+00, &
          0.3958885101147953D+00, &
          0.3958885101147953D+00, &
          0.3958885101147953D+00, &
          0.3958885101147953D+00, &
          0.3519734130158997D+00, &
          0.3519734130158997D+00, &
          0.3519734130158997D+00, &
          0.3519734130158997D+00, &
          0.3519734130158997D+00, &
          0.3519734130158997D+00, &
          0.6480265869841003D+00, &
          0.6480265869841003D+00, &
          0.6480265869841003D+00, &
          0.6480265869841003D+00, &
          0.6480265869841003D+00, &
          0.6480265869841003D+00, &
          0.9783845500663495D+00, &
          0.9783845500663495D+00, &
          0.9783845500663495D+00, &
          0.9783845500663495D+00, &
          0.9783845500663495D+00, &
          0.9783845500663495D+00, &
          0.0216154499336504D+00, &
          0.0216154499336504D+00, &
          0.0216154499336504D+00, &
          0.0216154499336504D+00, &
          0.0216154499336504D+00, &
          0.0216154499336504D+00, &
          0.9921911717290701D+00, &
          0.9921911717290701D+00, &
          0.9921911717290701D+00, &
          0.9921911717290701D+00, &
          0.9921911717290701D+00, &
          0.9921911717290701D+00, &
          0.0078088282709299D+00, &
          0.0078088282709299D+00, &
          0.0078088282709299D+00, &
          0.0078088282709299D+00, &
          0.0078088282709299D+00, &
          0.0078088282709299D+00, &
          0.8588015965488265D+00, &
          0.8588015965488265D+00, &
          0.8588015965488265D+00, &
          0.8588015965488265D+00, &
          0.8588015965488265D+00, &
          0.8588015965488265D+00, &
          0.1411984034511735D+00, &
          0.1411984034511735D+00, &
          0.1411984034511735D+00, &
          0.1411984034511735D+00, &
          0.1411984034511735D+00, &
          0.1411984034511735D+00, &
          0.2510713639088549D+00, &
          0.2510713639088549D+00, &
          0.2510713639088549D+00, &
          0.2510713639088549D+00, &
          0.2510713639088549D+00, &
          0.2510713639088549D+00, &
          0.7489286360911451D+00, &
          0.7489286360911451D+00, &
          0.7489286360911451D+00, &
          0.7489286360911451D+00, &
          0.7489286360911451D+00, &
          0.7489286360911451D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          0.0117723371203579D+00, &
          0.0117723371203579D+00, &
          0.0135361913058678D+00, &
          0.0135361913058678D+00, &
          0.0135361913058678D+00, &
          0.0028636244590190D+00, &
          0.0028636244590190D+00, &
          0.0028636244590190D+00, &
          0.0028317101953144D+00, &
          0.0028317101953144D+00, &
          0.0028317101953144D+00, &
          0.0028317101953144D+00, &
          0.0028317101953144D+00, &
          0.0028317101953144D+00, &
          0.0043546495593520D+00, &
          0.0043546495593520D+00, &
          0.0043546495593520D+00, &
          0.0043546495593520D+00, &
          0.0043546495593520D+00, &
          0.0043546495593520D+00, &
          0.0009145278212111D+00, &
          0.0009145278212111D+00, &
          0.0009145278212111D+00, &
          0.0009145278212111D+00, &
          0.0009145278212111D+00, &
          0.0009145278212111D+00, &
          0.0016406658653518D+00, &
          0.0016406658653518D+00, &
          0.0016406658653518D+00, &
          0.0016406658653518D+00, &
          0.0016406658653518D+00, &
          0.0016406658653518D+00, &
          0.0121732957893404D+00, &
          0.0121732957893404D+00, &
          0.0121732957893404D+00, &
          0.0121732957893404D+00, &
          0.0121732957893404D+00, &
          0.0121732957893404D+00, &
          0.0010621280231715D+00, &
          0.0010621280231715D+00, &
          0.0010621280231715D+00, &
          0.0010621280231715D+00, &
          0.0010621280231715D+00, &
          0.0010621280231715D+00, &
          0.0081462190313536D+00, &
          0.0081462190313536D+00, &
          0.0081462190313536D+00, &
          0.0081462190313536D+00, &
          0.0081462190313536D+00, &
          0.0081462190313536D+00, &
          0.0118819825086521D+00, &
          0.0118819825086521D+00, &
          0.0118819825086521D+00, &
          0.0118819825086521D+00, &
          0.0118819825086521D+00, &
          0.0118819825086521D+00, &
          0.0045106060735437D+00, &
          0.0045106060735437D+00, &
          0.0045106060735437D+00, &
          0.0045106060735437D+00, &
          0.0045106060735437D+00, &
          0.0045106060735437D+00, &
          0.0120105850348785D+00, &
          0.0120105850348785D+00, &
          0.0120105850348785D+00, &
          0.0120105850348785D+00, &
          0.0120105850348785D+00, &
          0.0120105850348785D+00, &
          0.0047664383021567D+00, &
          0.0047664383021567D+00, &
          0.0047664383021567D+00, &
          0.0047664383021567D+00, &
          0.0047664383021567D+00, &
          0.0047664383021567D+00, &
          0.0009777024948606D+00, &
          0.0009777024948606D+00, &
          0.0009777024948606D+00, &
          0.0009777024948606D+00, &
          0.0009777024948606D+00, &
          0.0009777024948606D+00, &
          0.0051118807026059D+00, &
          0.0051118807026059D+00, &
          0.0051118807026059D+00, &
          0.0051118807026059D+00, &
          0.0051118807026059D+00, &
          0.0051118807026059D+00, &
          0.0021717621070222D+00, &
          0.0021717621070222D+00, &
          0.0021717621070222D+00, &
          0.0021717621070222D+00, &
          0.0021717621070222D+00, &
          0.0021717621070222D+00, &
          0.0104022745738753D+00, &
          0.0104022745738753D+00, &
          0.0104022745738753D+00, &
          0.0104022745738753D+00, &
          0.0104022745738753D+00, &
          0.0104022745738753D+00, &
          0.0035903626892995D+00, &
          0.0035903626892995D+00, &
          0.0035903626892995D+00, &
          0.0035903626892995D+00, &
          0.0035903626892995D+00, &
          0.0035903626892995D+00, &
          0.0035903626892995D+00, &
          0.0035903626892995D+00, &
          0.0035903626892995D+00, &
          0.0035903626892995D+00, &
          0.0035903626892995D+00, &
          0.0035903626892995D+00, &
          0.0052108482781140D+00, &
          0.0052108482781140D+00, &
          0.0052108482781140D+00, &
          0.0052108482781140D+00, &
          0.0052108482781140D+00, &
          0.0052108482781140D+00, &
          0.0052108482781140D+00, &
          0.0052108482781140D+00, &
          0.0052108482781140D+00, &
          0.0052108482781140D+00, &
          0.0052108482781140D+00, &
          0.0052108482781140D+00, &
          0.0029024210927214D+00, &
          0.0029024210927214D+00, &
          0.0029024210927214D+00, &
          0.0029024210927214D+00, &
          0.0029024210927214D+00, &
          0.0029024210927214D+00, &
          0.0029024210927214D+00, &
          0.0029024210927214D+00, &
          0.0029024210927214D+00, &
          0.0029024210927214D+00, &
          0.0029024210927214D+00, &
          0.0029024210927214D+00, &
          0.0109370973987031D+00, &
          0.0109370973987031D+00, &
          0.0109370973987031D+00, &
          0.0109370973987031D+00, &
          0.0109370973987031D+00, &
          0.0109370973987031D+00, &
          0.0109370973987031D+00, &
          0.0109370973987031D+00, &
          0.0109370973987031D+00, &
          0.0109370973987031D+00, &
          0.0109370973987031D+00, &
          0.0109370973987031D+00, &
          0.0010706231542004D+00, &
          0.0010706231542004D+00, &
          0.0010706231542004D+00, &
          0.0010706231542004D+00, &
          0.0010706231542004D+00, &
          0.0010706231542004D+00, &
          0.0010706231542004D+00, &
          0.0010706231542004D+00, &
          0.0010706231542004D+00, &
          0.0010706231542004D+00, &
          0.0010706231542004D+00, &
          0.0010706231542004D+00, &
          0.0023660542319623D+00, &
          0.0023660542319623D+00, &
          0.0023660542319623D+00, &
          0.0023660542319623D+00, &
          0.0023660542319623D+00, &
          0.0023660542319623D+00, &
          0.0023660542319623D+00, &
          0.0023660542319623D+00, &
          0.0023660542319623D+00, &
          0.0023660542319623D+00, &
          0.0023660542319623D+00, &
          0.0023660542319623D+00, &
          0.0026196720900109D+00, &
          0.0026196720900109D+00, &
          0.0026196720900109D+00, &
          0.0026196720900109D+00, &
          0.0026196720900109D+00, &
          0.0026196720900109D+00, &
          0.0026196720900109D+00, &
          0.0026196720900109D+00, &
          0.0026196720900109D+00, &
          0.0026196720900109D+00, &
          0.0026196720900109D+00, &
          0.0026196720900109D+00, &
          0.0070960302290289D+00, &
          0.0070960302290289D+00, &
          0.0070960302290289D+00, &
          0.0070960302290289D+00, &
          0.0070960302290289D+00, &
          0.0070960302290289D+00, &
          0.0070960302290289D+00, &
          0.0070960302290289D+00, &
          0.0070960302290289D+00, &
          0.0070960302290289D+00, &
          0.0070960302290289D+00, &
          0.0070960302290289D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule15 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule15() returns the prism rule of precision 15.
!
!  Discussion:
!
!    We suppose we are given a triangular prism P with vertices
!    (1,0,0), (0,1,0), (0,0,0),
!    (1,0,1), (0,1,1), (0,0,1),
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over P is 
!    approximated by Q as follows:
!
!    Q = volume(P) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 May 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order symmetric cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 122, Number 1, pages 148-171, 24 August 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real x(n), y(n), z(n): the coordinates of quadrature points.
!
!    real w(n): the quadrature weights, which sum to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: n_save = 238

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.2180179116942069D+00, &
          0.2180179116942069D+00, &
          0.5639641766115862D+00, &
          0.1819996803717089D+00, &
          0.1819996803717089D+00, &
          0.6360006392565821D+00, &
          0.0174014011631261D+00, &
          0.0174014011631261D+00, &
          0.9651971976737478D+00, &
          0.0542711501859966D+00, &
          0.0542711501859966D+00, &
          0.8914576996280067D+00, &
          0.3702896314191779D+00, &
          0.3702896314191779D+00, &
          0.2594207371616443D+00, &
          0.3702896314191779D+00, &
          0.3702896314191779D+00, &
          0.2594207371616443D+00, &
          0.0683792324219036D+00, &
          0.0683792324219036D+00, &
          0.8632415351561927D+00, &
          0.0683792324219036D+00, &
          0.0683792324219036D+00, &
          0.8632415351561927D+00, &
          0.4041029865296590D+00, &
          0.4041029865296590D+00, &
          0.1917940269406819D+00, &
          0.4041029865296590D+00, &
          0.4041029865296590D+00, &
          0.1917940269406819D+00, &
          0.4924009437532922D+00, &
          0.4924009437532922D+00, &
          0.0151981124934156D+00, &
          0.4924009437532922D+00, &
          0.4924009437532922D+00, &
          0.0151981124934156D+00, &
          0.4679592294726652D+00, &
          0.4679592294726652D+00, &
          0.0640815410546695D+00, &
          0.4679592294726652D+00, &
          0.4679592294726652D+00, &
          0.0640815410546695D+00, &
          0.2538613877625208D+00, &
          0.2538613877625208D+00, &
          0.4922772244749584D+00, &
          0.2538613877625208D+00, &
          0.2538613877625208D+00, &
          0.4922772244749584D+00, &
          0.1490233678882002D+00, &
          0.1490233678882002D+00, &
          0.7019532642235996D+00, &
          0.1490233678882002D+00, &
          0.1490233678882002D+00, &
          0.7019532642235996D+00, &
          0.4943016520035221D+00, &
          0.4943016520035221D+00, &
          0.0113966959929558D+00, &
          0.4943016520035221D+00, &
          0.4943016520035221D+00, &
          0.0113966959929558D+00, &
          0.2380227290434666D+00, &
          0.0191912097570085D+00, &
          0.2380227290434666D+00, &
          0.7427860611995248D+00, &
          0.0191912097570085D+00, &
          0.7427860611995248D+00, &
          0.2658819381343187D+00, &
          0.0947024276247834D+00, &
          0.2658819381343187D+00, &
          0.6394156342408979D+00, &
          0.0947024276247834D+00, &
          0.6394156342408979D+00, &
          0.2658819381343187D+00, &
          0.0947024276247834D+00, &
          0.2658819381343187D+00, &
          0.6394156342408979D+00, &
          0.0947024276247834D+00, &
          0.6394156342408979D+00, &
          0.0282488031371648D+00, &
          0.0067976462534296D+00, &
          0.0282488031371648D+00, &
          0.9649535506094056D+00, &
          0.0067976462534296D+00, &
          0.9649535506094056D+00, &
          0.0282488031371648D+00, &
          0.0067976462534296D+00, &
          0.0282488031371648D+00, &
          0.9649535506094056D+00, &
          0.0067976462534296D+00, &
          0.9649535506094056D+00, &
          0.6410939576871652D+00, &
          0.0151121945323486D+00, &
          0.6410939576871652D+00, &
          0.3437938477804863D+00, &
          0.0151121945323486D+00, &
          0.3437938477804863D+00, &
          0.6410939576871652D+00, &
          0.0151121945323486D+00, &
          0.6410939576871652D+00, &
          0.3437938477804863D+00, &
          0.0151121945323486D+00, &
          0.3437938477804863D+00, &
          0.1570345607673937D+00, &
          0.0792048851268470D+00, &
          0.1570345607673937D+00, &
          0.7637605541057593D+00, &
          0.0792048851268470D+00, &
          0.7637605541057593D+00, &
          0.1570345607673937D+00, &
          0.0792048851268470D+00, &
          0.1570345607673937D+00, &
          0.7637605541057593D+00, &
          0.0792048851268470D+00, &
          0.7637605541057593D+00, &
          0.2618681572601245D+00, &
          0.1978262294209248D+00, &
          0.2618681572601245D+00, &
          0.5403056133189508D+00, &
          0.1978262294209248D+00, &
          0.5403056133189508D+00, &
          0.2618681572601245D+00, &
          0.1978262294209248D+00, &
          0.2618681572601245D+00, &
          0.5403056133189508D+00, &
          0.1978262294209248D+00, &
          0.5403056133189508D+00, &
          0.0065324588746805D+00, &
          0.0486805622526683D+00, &
          0.0065324588746805D+00, &
          0.9447869788726512D+00, &
          0.0486805622526683D+00, &
          0.9447869788726512D+00, &
          0.0065324588746805D+00, &
          0.0486805622526683D+00, &
          0.0065324588746805D+00, &
          0.9447869788726512D+00, &
          0.0486805622526683D+00, &
          0.9447869788726512D+00, &
          0.3702624031957258D+00, &
          0.0912977133028943D+00, &
          0.3702624031957258D+00, &
          0.5384398835013798D+00, &
          0.0912977133028943D+00, &
          0.5384398835013798D+00, &
          0.3702624031957258D+00, &
          0.0912977133028943D+00, &
          0.3702624031957258D+00, &
          0.5384398835013798D+00, &
          0.0912977133028943D+00, &
          0.5384398835013798D+00, &
          0.1377079768509805D+00, &
          0.0373008369533812D+00, &
          0.1377079768509805D+00, &
          0.8249911861956383D+00, &
          0.0373008369533812D+00, &
          0.8249911861956383D+00, &
          0.1377079768509805D+00, &
          0.0373008369533812D+00, &
          0.1377079768509805D+00, &
          0.8249911861956383D+00, &
          0.0373008369533812D+00, &
          0.8249911861956383D+00, &
          0.4784950101267365D+00, &
          0.1635654688089934D+00, &
          0.4784950101267365D+00, &
          0.3579395210642701D+00, &
          0.1635654688089934D+00, &
          0.3579395210642701D+00, &
          0.4784950101267365D+00, &
          0.1635654688089934D+00, &
          0.4784950101267365D+00, &
          0.3579395210642701D+00, &
          0.1635654688089934D+00, &
          0.3579395210642701D+00, &
          0.3278240979071815D+00, &
          0.0249615816629217D+00, &
          0.3278240979071815D+00, &
          0.6472143204298968D+00, &
          0.0249615816629217D+00, &
          0.6472143204298968D+00, &
          0.3278240979071815D+00, &
          0.0249615816629217D+00, &
          0.3278240979071815D+00, &
          0.6472143204298968D+00, &
          0.0249615816629217D+00, &
          0.6472143204298968D+00, &
          0.3475735063324100D+00, &
          0.0871302288391398D+00, &
          0.3475735063324100D+00, &
          0.5652962648284502D+00, &
          0.0871302288391398D+00, &
          0.5652962648284502D+00, &
          0.3475735063324100D+00, &
          0.0871302288391398D+00, &
          0.3475735063324100D+00, &
          0.5652962648284502D+00, &
          0.0871302288391398D+00, &
          0.5652962648284502D+00, &
          0.1096327891311607D+00, &
          0.0128395677390919D+00, &
          0.1096327891311607D+00, &
          0.8775276431297473D+00, &
          0.0128395677390919D+00, &
          0.8775276431297473D+00, &
          0.1096327891311607D+00, &
          0.0128395677390919D+00, &
          0.1096327891311607D+00, &
          0.8775276431297473D+00, &
          0.0128395677390919D+00, &
          0.8775276431297473D+00, &
          0.2078784801536998D+00, &
          0.0849045675529205D+00, &
          0.2078784801536998D+00, &
          0.7072169522933797D+00, &
          0.0849045675529205D+00, &
          0.7072169522933797D+00, &
          0.2078784801536998D+00, &
          0.0849045675529205D+00, &
          0.2078784801536998D+00, &
          0.7072169522933797D+00, &
          0.0849045675529205D+00, &
          0.7072169522933797D+00, &
          0.1953330193360589D+00, &
          0.0143045787175747D+00, &
          0.1953330193360589D+00, &
          0.7903624019463664D+00, &
          0.0143045787175747D+00, &
          0.7903624019463664D+00, &
          0.1953330193360589D+00, &
          0.0143045787175747D+00, &
          0.1953330193360589D+00, &
          0.7903624019463664D+00, &
          0.0143045787175747D+00, &
          0.7903624019463664D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.2180179116942069D+00, &
          0.5639641766115862D+00, &
          0.2180179116942069D+00, &
          0.1819996803717089D+00, &
          0.6360006392565821D+00, &
          0.1819996803717089D+00, &
          0.0174014011631261D+00, &
          0.9651971976737478D+00, &
          0.0174014011631261D+00, &
          0.0542711501859966D+00, &
          0.8914576996280067D+00, &
          0.0542711501859966D+00, &
          0.3702896314191779D+00, &
          0.2594207371616443D+00, &
          0.3702896314191779D+00, &
          0.3702896314191779D+00, &
          0.2594207371616443D+00, &
          0.3702896314191779D+00, &
          0.0683792324219036D+00, &
          0.8632415351561927D+00, &
          0.0683792324219036D+00, &
          0.0683792324219036D+00, &
          0.8632415351561927D+00, &
          0.0683792324219036D+00, &
          0.4041029865296590D+00, &
          0.1917940269406819D+00, &
          0.4041029865296590D+00, &
          0.4041029865296590D+00, &
          0.1917940269406819D+00, &
          0.4041029865296590D+00, &
          0.4924009437532922D+00, &
          0.0151981124934156D+00, &
          0.4924009437532922D+00, &
          0.4924009437532922D+00, &
          0.0151981124934156D+00, &
          0.4924009437532922D+00, &
          0.4679592294726652D+00, &
          0.0640815410546695D+00, &
          0.4679592294726652D+00, &
          0.4679592294726652D+00, &
          0.0640815410546695D+00, &
          0.4679592294726652D+00, &
          0.2538613877625208D+00, &
          0.4922772244749584D+00, &
          0.2538613877625208D+00, &
          0.2538613877625208D+00, &
          0.4922772244749584D+00, &
          0.2538613877625208D+00, &
          0.1490233678882002D+00, &
          0.7019532642235996D+00, &
          0.1490233678882002D+00, &
          0.1490233678882002D+00, &
          0.7019532642235996D+00, &
          0.1490233678882002D+00, &
          0.4943016520035221D+00, &
          0.0113966959929558D+00, &
          0.4943016520035221D+00, &
          0.4943016520035221D+00, &
          0.0113966959929558D+00, &
          0.4943016520035221D+00, &
          0.0191912097570085D+00, &
          0.2380227290434666D+00, &
          0.7427860611995248D+00, &
          0.2380227290434666D+00, &
          0.7427860611995248D+00, &
          0.0191912097570085D+00, &
          0.0947024276247834D+00, &
          0.2658819381343187D+00, &
          0.6394156342408979D+00, &
          0.2658819381343187D+00, &
          0.6394156342408979D+00, &
          0.0947024276247834D+00, &
          0.0947024276247834D+00, &
          0.2658819381343187D+00, &
          0.6394156342408979D+00, &
          0.2658819381343187D+00, &
          0.6394156342408979D+00, &
          0.0947024276247834D+00, &
          0.0067976462534296D+00, &
          0.0282488031371648D+00, &
          0.9649535506094056D+00, &
          0.0282488031371648D+00, &
          0.9649535506094056D+00, &
          0.0067976462534296D+00, &
          0.0067976462534296D+00, &
          0.0282488031371648D+00, &
          0.9649535506094056D+00, &
          0.0282488031371648D+00, &
          0.9649535506094056D+00, &
          0.0067976462534296D+00, &
          0.0151121945323486D+00, &
          0.6410939576871652D+00, &
          0.3437938477804863D+00, &
          0.6410939576871652D+00, &
          0.3437938477804863D+00, &
          0.0151121945323486D+00, &
          0.0151121945323486D+00, &
          0.6410939576871652D+00, &
          0.3437938477804863D+00, &
          0.6410939576871652D+00, &
          0.3437938477804863D+00, &
          0.0151121945323486D+00, &
          0.0792048851268470D+00, &
          0.1570345607673937D+00, &
          0.7637605541057593D+00, &
          0.1570345607673937D+00, &
          0.7637605541057593D+00, &
          0.0792048851268470D+00, &
          0.0792048851268470D+00, &
          0.1570345607673937D+00, &
          0.7637605541057593D+00, &
          0.1570345607673937D+00, &
          0.7637605541057593D+00, &
          0.0792048851268470D+00, &
          0.1978262294209248D+00, &
          0.2618681572601245D+00, &
          0.5403056133189508D+00, &
          0.2618681572601245D+00, &
          0.5403056133189508D+00, &
          0.1978262294209248D+00, &
          0.1978262294209248D+00, &
          0.2618681572601245D+00, &
          0.5403056133189508D+00, &
          0.2618681572601245D+00, &
          0.5403056133189508D+00, &
          0.1978262294209248D+00, &
          0.0486805622526683D+00, &
          0.0065324588746805D+00, &
          0.9447869788726512D+00, &
          0.0065324588746805D+00, &
          0.9447869788726512D+00, &
          0.0486805622526683D+00, &
          0.0486805622526683D+00, &
          0.0065324588746805D+00, &
          0.9447869788726512D+00, &
          0.0065324588746805D+00, &
          0.9447869788726512D+00, &
          0.0486805622526683D+00, &
          0.0912977133028943D+00, &
          0.3702624031957258D+00, &
          0.5384398835013798D+00, &
          0.3702624031957258D+00, &
          0.5384398835013798D+00, &
          0.0912977133028943D+00, &
          0.0912977133028943D+00, &
          0.3702624031957258D+00, &
          0.5384398835013798D+00, &
          0.3702624031957258D+00, &
          0.5384398835013798D+00, &
          0.0912977133028943D+00, &
          0.0373008369533812D+00, &
          0.1377079768509805D+00, &
          0.8249911861956383D+00, &
          0.1377079768509805D+00, &
          0.8249911861956383D+00, &
          0.0373008369533812D+00, &
          0.0373008369533812D+00, &
          0.1377079768509805D+00, &
          0.8249911861956383D+00, &
          0.1377079768509805D+00, &
          0.8249911861956383D+00, &
          0.0373008369533812D+00, &
          0.1635654688089934D+00, &
          0.4784950101267365D+00, &
          0.3579395210642701D+00, &
          0.4784950101267365D+00, &
          0.3579395210642701D+00, &
          0.1635654688089934D+00, &
          0.1635654688089934D+00, &
          0.4784950101267365D+00, &
          0.3579395210642701D+00, &
          0.4784950101267365D+00, &
          0.3579395210642701D+00, &
          0.1635654688089934D+00, &
          0.0249615816629217D+00, &
          0.3278240979071815D+00, &
          0.6472143204298968D+00, &
          0.3278240979071815D+00, &
          0.6472143204298968D+00, &
          0.0249615816629217D+00, &
          0.0249615816629217D+00, &
          0.3278240979071815D+00, &
          0.6472143204298968D+00, &
          0.3278240979071815D+00, &
          0.6472143204298968D+00, &
          0.0249615816629217D+00, &
          0.0871302288391398D+00, &
          0.3475735063324100D+00, &
          0.5652962648284502D+00, &
          0.3475735063324100D+00, &
          0.5652962648284502D+00, &
          0.0871302288391398D+00, &
          0.0871302288391398D+00, &
          0.3475735063324100D+00, &
          0.5652962648284502D+00, &
          0.3475735063324100D+00, &
          0.5652962648284502D+00, &
          0.0871302288391398D+00, &
          0.0128395677390919D+00, &
          0.1096327891311607D+00, &
          0.8775276431297473D+00, &
          0.1096327891311607D+00, &
          0.8775276431297473D+00, &
          0.0128395677390919D+00, &
          0.0128395677390919D+00, &
          0.1096327891311607D+00, &
          0.8775276431297473D+00, &
          0.1096327891311607D+00, &
          0.8775276431297473D+00, &
          0.0128395677390919D+00, &
          0.0849045675529205D+00, &
          0.2078784801536998D+00, &
          0.7072169522933797D+00, &
          0.2078784801536998D+00, &
          0.7072169522933797D+00, &
          0.0849045675529205D+00, &
          0.0849045675529205D+00, &
          0.2078784801536998D+00, &
          0.7072169522933797D+00, &
          0.2078784801536998D+00, &
          0.7072169522933797D+00, &
          0.0849045675529205D+00, &
          0.0143045787175747D+00, &
          0.1953330193360589D+00, &
          0.7903624019463664D+00, &
          0.1953330193360589D+00, &
          0.7903624019463664D+00, &
          0.0143045787175747D+00, &
          0.0143045787175747D+00, &
          0.1953330193360589D+00, &
          0.7903624019463664D+00, &
          0.1953330193360589D+00, &
          0.7903624019463664D+00, &
          0.0143045787175747D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.8401109114804578D+00, &
          0.1598890885195422D+00, &
          0.9524833096339549D+00, &
          0.0475166903660451D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5508863349912149D+00, &
          0.5508863349912149D+00, &
          0.5508863349912149D+00, &
          0.4491136650087851D+00, &
          0.4491136650087851D+00, &
          0.4491136650087851D+00, &
          0.8664406024867690D+00, &
          0.8664406024867690D+00, &
          0.8664406024867690D+00, &
          0.1335593975132310D+00, &
          0.1335593975132310D+00, &
          0.1335593975132310D+00, &
          0.3291089745865182D+00, &
          0.3291089745865182D+00, &
          0.3291089745865182D+00, &
          0.6708910254134818D+00, &
          0.6708910254134818D+00, &
          0.6708910254134818D+00, &
          0.4287046442912852D+00, &
          0.4287046442912852D+00, &
          0.4287046442912852D+00, &
          0.5712953557087147D+00, &
          0.5712953557087147D+00, &
          0.5712953557087147D+00, &
          0.8306396134186287D+00, &
          0.8306396134186287D+00, &
          0.8306396134186287D+00, &
          0.1693603865813713D+00, &
          0.1693603865813713D+00, &
          0.1693603865813713D+00, &
          0.9826991514876333D+00, &
          0.9826991514876333D+00, &
          0.9826991514876333D+00, &
          0.0173008485123666D+00, &
          0.0173008485123666D+00, &
          0.0173008485123666D+00, &
          0.9603593395589844D+00, &
          0.9603593395589844D+00, &
          0.9603593395589844D+00, &
          0.0396406604410156D+00, &
          0.0396406604410156D+00, &
          0.0396406604410156D+00, &
          0.9264398795821308D+00, &
          0.9264398795821308D+00, &
          0.9264398795821308D+00, &
          0.0735601204178692D+00, &
          0.0735601204178692D+00, &
          0.0735601204178692D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9025925504145782D+00, &
          0.9025925504145782D+00, &
          0.9025925504145782D+00, &
          0.9025925504145782D+00, &
          0.9025925504145782D+00, &
          0.9025925504145782D+00, &
          0.0974074495854217D+00, &
          0.0974074495854217D+00, &
          0.0974074495854217D+00, &
          0.0974074495854217D+00, &
          0.0974074495854217D+00, &
          0.0974074495854217D+00, &
          0.7775063050268207D+00, &
          0.7775063050268207D+00, &
          0.7775063050268207D+00, &
          0.7775063050268207D+00, &
          0.7775063050268207D+00, &
          0.7775063050268207D+00, &
          0.2224936949731793D+00, &
          0.2224936949731793D+00, &
          0.2224936949731793D+00, &
          0.2224936949731793D+00, &
          0.2224936949731793D+00, &
          0.2224936949731793D+00, &
          0.7482800508722226D+00, &
          0.7482800508722226D+00, &
          0.7482800508722226D+00, &
          0.7482800508722226D+00, &
          0.7482800508722226D+00, &
          0.7482800508722226D+00, &
          0.2517199491277774D+00, &
          0.2517199491277774D+00, &
          0.2517199491277774D+00, &
          0.2517199491277774D+00, &
          0.2517199491277774D+00, &
          0.2517199491277774D+00, &
          0.6403701196964353D+00, &
          0.6403701196964353D+00, &
          0.6403701196964353D+00, &
          0.6403701196964353D+00, &
          0.6403701196964353D+00, &
          0.6403701196964353D+00, &
          0.3596298803035647D+00, &
          0.3596298803035647D+00, &
          0.3596298803035647D+00, &
          0.3596298803035647D+00, &
          0.3596298803035647D+00, &
          0.3596298803035647D+00, &
          0.7769630399535563D+00, &
          0.7769630399535563D+00, &
          0.7769630399535563D+00, &
          0.7769630399535563D+00, &
          0.7769630399535563D+00, &
          0.7769630399535563D+00, &
          0.2230369600464437D+00, &
          0.2230369600464437D+00, &
          0.2230369600464437D+00, &
          0.2230369600464437D+00, &
          0.2230369600464437D+00, &
          0.2230369600464437D+00, &
          0.9565498865251021D+00, &
          0.9565498865251021D+00, &
          0.9565498865251021D+00, &
          0.9565498865251021D+00, &
          0.9565498865251021D+00, &
          0.9565498865251021D+00, &
          0.0434501134748980D+00, &
          0.0434501134748980D+00, &
          0.0434501134748980D+00, &
          0.0434501134748980D+00, &
          0.0434501134748980D+00, &
          0.0434501134748980D+00, &
          0.9928793309605464D+00, &
          0.9928793309605464D+00, &
          0.9928793309605464D+00, &
          0.9928793309605464D+00, &
          0.9928793309605464D+00, &
          0.9928793309605464D+00, &
          0.0071206690394536D+00, &
          0.0071206690394536D+00, &
          0.0071206690394536D+00, &
          0.0071206690394536D+00, &
          0.0071206690394536D+00, &
          0.0071206690394536D+00, &
          0.9880631245383409D+00, &
          0.9880631245383409D+00, &
          0.9880631245383409D+00, &
          0.9880631245383409D+00, &
          0.9880631245383409D+00, &
          0.9880631245383409D+00, &
          0.0119368754616592D+00, &
          0.0119368754616592D+00, &
          0.0119368754616592D+00, &
          0.0119368754616592D+00, &
          0.0119368754616592D+00, &
          0.0119368754616592D+00, &
          0.9154311485206492D+00, &
          0.9154311485206492D+00, &
          0.9154311485206492D+00, &
          0.9154311485206492D+00, &
          0.9154311485206492D+00, &
          0.9154311485206492D+00, &
          0.0845688514793507D+00, &
          0.0845688514793507D+00, &
          0.0845688514793507D+00, &
          0.0845688514793507D+00, &
          0.0845688514793507D+00, &
          0.0845688514793507D+00, &
          0.9705250000766562D+00, &
          0.9705250000766562D+00, &
          0.9705250000766562D+00, &
          0.9705250000766562D+00, &
          0.9705250000766562D+00, &
          0.9705250000766562D+00, &
          0.0294749999233438D+00, &
          0.0294749999233438D+00, &
          0.0294749999233438D+00, &
          0.0294749999233438D+00, &
          0.0294749999233438D+00, &
          0.0294749999233438D+00, &
          0.3853017195325944D+00, &
          0.3853017195325944D+00, &
          0.3853017195325944D+00, &
          0.3853017195325944D+00, &
          0.3853017195325944D+00, &
          0.3853017195325944D+00, &
          0.6146982804674056D+00, &
          0.6146982804674056D+00, &
          0.6146982804674056D+00, &
          0.6146982804674056D+00, &
          0.6146982804674056D+00, &
          0.6146982804674056D+00, &
          0.6900924303976952D+00, &
          0.6900924303976952D+00, &
          0.6900924303976952D+00, &
          0.6900924303976952D+00, &
          0.6900924303976952D+00, &
          0.6900924303976952D+00, &
          0.3099075696023049D+00, &
          0.3099075696023049D+00, &
          0.3099075696023049D+00, &
          0.3099075696023049D+00, &
          0.3099075696023049D+00, &
          0.3099075696023049D+00, &
          0.7816502945768592D+00, &
          0.7816502945768592D+00, &
          0.7816502945768592D+00, &
          0.7816502945768592D+00, &
          0.7816502945768592D+00, &
          0.7816502945768592D+00, &
          0.2183497054231408D+00, &
          0.2183497054231408D+00, &
          0.2183497054231408D+00, &
          0.2183497054231408D+00, &
          0.2183497054231408D+00, &
          0.2183497054231408D+00, &
          0.8883851866879513D+00, &
          0.8883851866879513D+00, &
          0.8883851866879513D+00, &
          0.8883851866879513D+00, &
          0.8883851866879513D+00, &
          0.8883851866879513D+00, &
          0.1116148133120487D+00, &
          0.1116148133120487D+00, &
          0.1116148133120487D+00, &
          0.1116148133120487D+00, &
          0.1116148133120487D+00, &
          0.1116148133120487D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          0.0109734885358618D+00, &
          0.0109734885358618D+00, &
          0.0032613958307349D+00, &
          0.0032613958307349D+00, &
          0.0093941405464380D+00, &
          0.0093941405464380D+00, &
          0.0093941405464380D+00, &
          0.0074721925384407D+00, &
          0.0074721925384407D+00, &
          0.0074721925384407D+00, &
          0.0011251486610045D+00, &
          0.0011251486610045D+00, &
          0.0011251486610045D+00, &
          0.0040278942594456D+00, &
          0.0040278942594456D+00, &
          0.0040278942594456D+00, &
          0.0046987405274126D+00, &
          0.0046987405274126D+00, &
          0.0046987405274126D+00, &
          0.0046987405274126D+00, &
          0.0046987405274126D+00, &
          0.0046987405274126D+00, &
          0.0044231068731859D+00, &
          0.0044231068731859D+00, &
          0.0044231068731859D+00, &
          0.0044231068731859D+00, &
          0.0044231068731859D+00, &
          0.0044231068731859D+00, &
          0.0119087119167916D+00, &
          0.0119087119167916D+00, &
          0.0119087119167916D+00, &
          0.0119087119167916D+00, &
          0.0119087119167916D+00, &
          0.0119087119167916D+00, &
          0.0034190583552652D+00, &
          0.0034190583552652D+00, &
          0.0034190583552652D+00, &
          0.0034190583552652D+00, &
          0.0034190583552652D+00, &
          0.0034190583552652D+00, &
          0.0076078695830692D+00, &
          0.0076078695830692D+00, &
          0.0076078695830692D+00, &
          0.0076078695830692D+00, &
          0.0076078695830692D+00, &
          0.0076078695830692D+00, &
          0.0034877218143443D+00, &
          0.0034877218143443D+00, &
          0.0034877218143443D+00, &
          0.0034877218143443D+00, &
          0.0034877218143443D+00, &
          0.0034877218143443D+00, &
          0.0041265711237682D+00, &
          0.0041265711237682D+00, &
          0.0041265711237682D+00, &
          0.0041265711237682D+00, &
          0.0041265711237682D+00, &
          0.0041265711237682D+00, &
          0.0020552290439392D+00, &
          0.0020552290439392D+00, &
          0.0020552290439392D+00, &
          0.0020552290439392D+00, &
          0.0020552290439392D+00, &
          0.0020552290439392D+00, &
          0.0053334369811175D+00, &
          0.0053334369811175D+00, &
          0.0053334369811175D+00, &
          0.0053334369811175D+00, &
          0.0053334369811175D+00, &
          0.0053334369811175D+00, &
          0.0045658362940629D+00, &
          0.0045658362940629D+00, &
          0.0045658362940629D+00, &
          0.0045658362940629D+00, &
          0.0045658362940629D+00, &
          0.0045658362940629D+00, &
          0.0045658362940629D+00, &
          0.0045658362940629D+00, &
          0.0045658362940629D+00, &
          0.0045658362940629D+00, &
          0.0045658362940629D+00, &
          0.0045658362940629D+00, &
          0.0005459132755231D+00, &
          0.0005459132755231D+00, &
          0.0005459132755231D+00, &
          0.0005459132755231D+00, &
          0.0005459132755231D+00, &
          0.0005459132755231D+00, &
          0.0005459132755231D+00, &
          0.0005459132755231D+00, &
          0.0005459132755231D+00, &
          0.0005459132755231D+00, &
          0.0005459132755231D+00, &
          0.0005459132755231D+00, &
          0.0043191948291063D+00, &
          0.0043191948291063D+00, &
          0.0043191948291063D+00, &
          0.0043191948291063D+00, &
          0.0043191948291063D+00, &
          0.0043191948291063D+00, &
          0.0043191948291063D+00, &
          0.0043191948291063D+00, &
          0.0043191948291063D+00, &
          0.0043191948291063D+00, &
          0.0043191948291063D+00, &
          0.0043191948291063D+00, &
          0.0053774143502532D+00, &
          0.0053774143502532D+00, &
          0.0053774143502532D+00, &
          0.0053774143502532D+00, &
          0.0053774143502532D+00, &
          0.0053774143502532D+00, &
          0.0053774143502532D+00, &
          0.0053774143502532D+00, &
          0.0053774143502532D+00, &
          0.0053774143502532D+00, &
          0.0053774143502532D+00, &
          0.0053774143502532D+00, &
          0.0071616945998409D+00, &
          0.0071616945998409D+00, &
          0.0071616945998409D+00, &
          0.0071616945998409D+00, &
          0.0071616945998409D+00, &
          0.0071616945998409D+00, &
          0.0071616945998409D+00, &
          0.0071616945998409D+00, &
          0.0071616945998409D+00, &
          0.0071616945998409D+00, &
          0.0071616945998409D+00, &
          0.0071616945998409D+00, &
          0.0005425741614460D+00, &
          0.0005425741614460D+00, &
          0.0005425741614460D+00, &
          0.0005425741614460D+00, &
          0.0005425741614460D+00, &
          0.0005425741614460D+00, &
          0.0005425741614460D+00, &
          0.0005425741614460D+00, &
          0.0005425741614460D+00, &
          0.0005425741614460D+00, &
          0.0005425741614460D+00, &
          0.0005425741614460D+00, &
          0.0014376691849096D+00, &
          0.0014376691849096D+00, &
          0.0014376691849096D+00, &
          0.0014376691849096D+00, &
          0.0014376691849096D+00, &
          0.0014376691849096D+00, &
          0.0014376691849096D+00, &
          0.0014376691849096D+00, &
          0.0014376691849096D+00, &
          0.0014376691849096D+00, &
          0.0014376691849096D+00, &
          0.0014376691849096D+00, &
          0.0011160390387250D+00, &
          0.0011160390387250D+00, &
          0.0011160390387250D+00, &
          0.0011160390387250D+00, &
          0.0011160390387250D+00, &
          0.0011160390387250D+00, &
          0.0011160390387250D+00, &
          0.0011160390387250D+00, &
          0.0011160390387250D+00, &
          0.0011160390387250D+00, &
          0.0011160390387250D+00, &
          0.0011160390387250D+00, &
          0.0048559812857500D+00, &
          0.0048559812857500D+00, &
          0.0048559812857500D+00, &
          0.0048559812857500D+00, &
          0.0048559812857500D+00, &
          0.0048559812857500D+00, &
          0.0048559812857500D+00, &
          0.0048559812857500D+00, &
          0.0048559812857500D+00, &
          0.0048559812857500D+00, &
          0.0048559812857500D+00, &
          0.0048559812857500D+00, &
          0.0017971531168555D+00, &
          0.0017971531168555D+00, &
          0.0017971531168555D+00, &
          0.0017971531168555D+00, &
          0.0017971531168555D+00, &
          0.0017971531168555D+00, &
          0.0017971531168555D+00, &
          0.0017971531168555D+00, &
          0.0017971531168555D+00, &
          0.0017971531168555D+00, &
          0.0017971531168555D+00, &
          0.0017971531168555D+00, &
          0.0096050585986694D+00, &
          0.0096050585986694D+00, &
          0.0096050585986694D+00, &
          0.0096050585986694D+00, &
          0.0096050585986694D+00, &
          0.0096050585986694D+00, &
          0.0096050585986694D+00, &
          0.0096050585986694D+00, &
          0.0096050585986694D+00, &
          0.0096050585986694D+00, &
          0.0096050585986694D+00, &
          0.0096050585986694D+00, &
          0.0027021277312842D+00, &
          0.0027021277312842D+00, &
          0.0027021277312842D+00, &
          0.0027021277312842D+00, &
          0.0027021277312842D+00, &
          0.0027021277312842D+00, &
          0.0027021277312842D+00, &
          0.0027021277312842D+00, &
          0.0027021277312842D+00, &
          0.0027021277312842D+00, &
          0.0027021277312842D+00, &
          0.0027021277312842D+00, &
          0.0055044460347859D+00, &
          0.0055044460347859D+00, &
          0.0055044460347859D+00, &
          0.0055044460347859D+00, &
          0.0055044460347859D+00, &
          0.0055044460347859D+00, &
          0.0055044460347859D+00, &
          0.0055044460347859D+00, &
          0.0055044460347859D+00, &
          0.0055044460347859D+00, &
          0.0055044460347859D+00, &
          0.0055044460347859D+00, &
          0.0023946829935762D+00, &
          0.0023946829935762D+00, &
          0.0023946829935762D+00, &
          0.0023946829935762D+00, &
          0.0023946829935762D+00, &
          0.0023946829935762D+00, &
          0.0023946829935762D+00, &
          0.0023946829935762D+00, &
          0.0023946829935762D+00, &
          0.0023946829935762D+00, &
          0.0023946829935762D+00, &
          0.0023946829935762D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule16 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule16() returns the prism rule of precision 16.
!
!  Discussion:
!
!    We suppose we are given a triangular prism P with vertices
!    (1,0,0), (0,1,0), (0,0,0),
!    (1,0,1), (0,1,1), (0,0,1),
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over P is 
!    approximated by Q as follows:
!
!    Q = volume(P) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 May 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order symmetric cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 122, Number 1, pages 148-171, 24 August 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real x(n), y(n), z(n): the coordinates of quadrature points.
!
!    real w(n): the quadrature weights, which sum to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: n_save = 287

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.0089540840312146D+00, &
          0.0089540840312146D+00, &
          0.9820918319375707D+00, &
          0.0553908100440216D+00, &
          0.0553908100440216D+00, &
          0.8892183799119567D+00, &
          0.3916317627335951D+00, &
          0.3916317627335951D+00, &
          0.2167364745328097D+00, &
          0.0561000778858347D+00, &
          0.0561000778858347D+00, &
          0.8877998442283307D+00, &
          0.0561000778858347D+00, &
          0.0561000778858347D+00, &
          0.8877998442283307D+00, &
          0.2669967593271997D+00, &
          0.2669967593271997D+00, &
          0.4660064813456005D+00, &
          0.2669967593271997D+00, &
          0.2669967593271997D+00, &
          0.4660064813456005D+00, &
          0.1619375672046858D+00, &
          0.1619375672046858D+00, &
          0.6761248655906285D+00, &
          0.1619375672046858D+00, &
          0.1619375672046858D+00, &
          0.6761248655906285D+00, &
          0.4259891438745869D+00, &
          0.4259891438745869D+00, &
          0.1480217122508263D+00, &
          0.4259891438745869D+00, &
          0.4259891438745869D+00, &
          0.1480217122508263D+00, &
          0.4943682456985981D+00, &
          0.4943682456985981D+00, &
          0.0112635086028039D+00, &
          0.4943682456985981D+00, &
          0.4943682456985981D+00, &
          0.0112635086028039D+00, &
          0.4745817664118026D+00, &
          0.4745817664118026D+00, &
          0.0508364671763949D+00, &
          0.4745817664118026D+00, &
          0.4745817664118026D+00, &
          0.0508364671763949D+00, &
          0.2785739189962955D+00, &
          0.2785739189962955D+00, &
          0.4428521620074090D+00, &
          0.2785739189962955D+00, &
          0.2785739189962955D+00, &
          0.4428521620074090D+00, &
          0.1564426776295505D+00, &
          0.1564426776295505D+00, &
          0.6871146447408989D+00, &
          0.1564426776295505D+00, &
          0.1564426776295505D+00, &
          0.6871146447408989D+00, &
          0.4933597406415689D+00, &
          0.4933597406415689D+00, &
          0.0132805187168622D+00, &
          0.4933597406415689D+00, &
          0.4933597406415689D+00, &
          0.0132805187168622D+00, &
          0.5254267651000835D+00, &
          0.0993149454903146D+00, &
          0.5254267651000835D+00, &
          0.3752582894096019D+00, &
          0.0993149454903146D+00, &
          0.3752582894096019D+00, &
          0.3482743570815052D+00, &
          0.0212715561472286D+00, &
          0.3482743570815052D+00, &
          0.6304540867712661D+00, &
          0.0212715561472286D+00, &
          0.6304540867712661D+00, &
          0.1654756461617898D+00, &
          0.0860850979983652D+00, &
          0.1654756461617898D+00, &
          0.7484392558398449D+00, &
          0.0860850979983652D+00, &
          0.7484392558398449D+00, &
          0.1639763763638563D+00, &
          0.0551722290870520D+00, &
          0.1639763763638563D+00, &
          0.7808513945490917D+00, &
          0.0551722290870520D+00, &
          0.7808513945490917D+00, &
          0.1639763763638563D+00, &
          0.0551722290870520D+00, &
          0.1639763763638563D+00, &
          0.7808513945490917D+00, &
          0.0551722290870520D+00, &
          0.7808513945490917D+00, &
          0.0370382522819345D+00, &
          0.0105036665914126D+00, &
          0.0370382522819345D+00, &
          0.9524580811266530D+00, &
          0.0105036665914126D+00, &
          0.9524580811266530D+00, &
          0.0370382522819345D+00, &
          0.0105036665914126D+00, &
          0.0370382522819345D+00, &
          0.9524580811266530D+00, &
          0.0105036665914126D+00, &
          0.9524580811266530D+00, &
          0.6350080118967869D+00, &
          0.0055639594221052D+00, &
          0.6350080118967869D+00, &
          0.3594280286811078D+00, &
          0.0055639594221052D+00, &
          0.3594280286811078D+00, &
          0.6350080118967869D+00, &
          0.0055639594221052D+00, &
          0.6350080118967869D+00, &
          0.3594280286811078D+00, &
          0.0055639594221052D+00, &
          0.3594280286811078D+00, &
          0.1297894012222233D+00, &
          0.0724868197207572D+00, &
          0.1297894012222233D+00, &
          0.7977237790570195D+00, &
          0.0724868197207572D+00, &
          0.7977237790570195D+00, &
          0.1297894012222233D+00, &
          0.0724868197207572D+00, &
          0.1297894012222233D+00, &
          0.7977237790570195D+00, &
          0.0724868197207572D+00, &
          0.7977237790570195D+00, &
          0.3059321030072774D+00, &
          0.0621201882612809D+00, &
          0.3059321030072774D+00, &
          0.6319477087314417D+00, &
          0.0621201882612809D+00, &
          0.6319477087314417D+00, &
          0.3059321030072774D+00, &
          0.0621201882612809D+00, &
          0.3059321030072774D+00, &
          0.6319477087314417D+00, &
          0.0621201882612809D+00, &
          0.6319477087314417D+00, &
          0.1677554820715698D+00, &
          0.2612303617350463D+00, &
          0.1677554820715698D+00, &
          0.5710141561933839D+00, &
          0.2612303617350463D+00, &
          0.5710141561933839D+00, &
          0.1677554820715698D+00, &
          0.2612303617350463D+00, &
          0.1677554820715698D+00, &
          0.5710141561933839D+00, &
          0.2612303617350463D+00, &
          0.5710141561933839D+00, &
          0.0015296219296330D+00, &
          0.0337686922990924D+00, &
          0.0015296219296330D+00, &
          0.9647016857712746D+00, &
          0.0337686922990924D+00, &
          0.9647016857712746D+00, &
          0.0015296219296330D+00, &
          0.0337686922990924D+00, &
          0.0015296219296330D+00, &
          0.9647016857712746D+00, &
          0.0337686922990924D+00, &
          0.9647016857712746D+00, &
          0.3532796606132861D+00, &
          0.0949242221588698D+00, &
          0.3532796606132861D+00, &
          0.5517961172278441D+00, &
          0.0949242221588698D+00, &
          0.5517961172278441D+00, &
          0.3532796606132861D+00, &
          0.0949242221588698D+00, &
          0.3532796606132861D+00, &
          0.5517961172278441D+00, &
          0.0949242221588698D+00, &
          0.5517961172278441D+00, &
          0.0959262325920420D+00, &
          0.0248817170289358D+00, &
          0.0959262325920420D+00, &
          0.8791920503790223D+00, &
          0.0248817170289358D+00, &
          0.8791920503790223D+00, &
          0.0959262325920420D+00, &
          0.0248817170289358D+00, &
          0.0959262325920420D+00, &
          0.8791920503790223D+00, &
          0.0248817170289358D+00, &
          0.8791920503790223D+00, &
          0.4963305263335406D+00, &
          0.1757029526928022D+00, &
          0.4963305263335406D+00, &
          0.3279665209736572D+00, &
          0.1757029526928022D+00, &
          0.3279665209736572D+00, &
          0.4963305263335406D+00, &
          0.1757029526928022D+00, &
          0.4963305263335406D+00, &
          0.3279665209736572D+00, &
          0.1757029526928022D+00, &
          0.3279665209736572D+00, &
          0.2738501035919119D+00, &
          0.0140786122092260D+00, &
          0.2738501035919119D+00, &
          0.7120712841988621D+00, &
          0.0140786122092260D+00, &
          0.7120712841988621D+00, &
          0.2738501035919119D+00, &
          0.0140786122092260D+00, &
          0.2738501035919119D+00, &
          0.7120712841988621D+00, &
          0.0140786122092260D+00, &
          0.7120712841988621D+00, &
          0.3828844998895772D+00, &
          0.0565171623200343D+00, &
          0.3828844998895772D+00, &
          0.5605983377903885D+00, &
          0.0565171623200343D+00, &
          0.5605983377903885D+00, &
          0.3828844998895772D+00, &
          0.0565171623200343D+00, &
          0.3828844998895772D+00, &
          0.5605983377903885D+00, &
          0.0565171623200343D+00, &
          0.5605983377903885D+00, &
          0.0899365323953004D+00, &
          0.0192412570201629D+00, &
          0.0899365323953004D+00, &
          0.8908222105845367D+00, &
          0.0192412570201629D+00, &
          0.8908222105845367D+00, &
          0.0899365323953004D+00, &
          0.0192412570201629D+00, &
          0.0899365323953004D+00, &
          0.8908222105845367D+00, &
          0.0192412570201629D+00, &
          0.8908222105845367D+00, &
          0.2516074150213084D+00, &
          0.0611370184542401D+00, &
          0.2516074150213084D+00, &
          0.6872555665244514D+00, &
          0.0611370184542401D+00, &
          0.6872555665244514D+00, &
          0.2516074150213084D+00, &
          0.0611370184542401D+00, &
          0.2516074150213084D+00, &
          0.6872555665244514D+00, &
          0.0611370184542401D+00, &
          0.6872555665244514D+00, &
          0.1524244609819593D+00, &
          0.0050182699561948D+00, &
          0.1524244609819593D+00, &
          0.8425572690618459D+00, &
          0.0050182699561948D+00, &
          0.8425572690618459D+00, &
          0.1524244609819593D+00, &
          0.0050182699561948D+00, &
          0.1524244609819593D+00, &
          0.8425572690618459D+00, &
          0.0050182699561948D+00, &
          0.8425572690618459D+00, &
          0.2475011786624472D+00, &
          0.1365910574355268D+00, &
          0.2475011786624472D+00, &
          0.6159077639020261D+00, &
          0.1365910574355268D+00, &
          0.6159077639020261D+00, &
          0.2475011786624472D+00, &
          0.1365910574355268D+00, &
          0.2475011786624472D+00, &
          0.6159077639020261D+00, &
          0.1365910574355268D+00, &
          0.6159077639020261D+00, &
          0.7591852986956938D+00, &
          0.0173543718109579D+00, &
          0.7591852986956938D+00, &
          0.2234603294933483D+00, &
          0.0173543718109579D+00, &
          0.2234603294933483D+00, &
          0.7591852986956938D+00, &
          0.0173543718109579D+00, &
          0.7591852986956938D+00, &
          0.2234603294933483D+00, &
          0.0173543718109579D+00, &
          0.2234603294933483D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.0089540840312146D+00, &
          0.9820918319375707D+00, &
          0.0089540840312146D+00, &
          0.0553908100440216D+00, &
          0.8892183799119567D+00, &
          0.0553908100440216D+00, &
          0.3916317627335951D+00, &
          0.2167364745328097D+00, &
          0.3916317627335951D+00, &
          0.0561000778858347D+00, &
          0.8877998442283307D+00, &
          0.0561000778858347D+00, &
          0.0561000778858347D+00, &
          0.8877998442283307D+00, &
          0.0561000778858347D+00, &
          0.2669967593271997D+00, &
          0.4660064813456005D+00, &
          0.2669967593271997D+00, &
          0.2669967593271997D+00, &
          0.4660064813456005D+00, &
          0.2669967593271997D+00, &
          0.1619375672046858D+00, &
          0.6761248655906285D+00, &
          0.1619375672046858D+00, &
          0.1619375672046858D+00, &
          0.6761248655906285D+00, &
          0.1619375672046858D+00, &
          0.4259891438745869D+00, &
          0.1480217122508263D+00, &
          0.4259891438745869D+00, &
          0.4259891438745869D+00, &
          0.1480217122508263D+00, &
          0.4259891438745869D+00, &
          0.4943682456985981D+00, &
          0.0112635086028039D+00, &
          0.4943682456985981D+00, &
          0.4943682456985981D+00, &
          0.0112635086028039D+00, &
          0.4943682456985981D+00, &
          0.4745817664118026D+00, &
          0.0508364671763949D+00, &
          0.4745817664118026D+00, &
          0.4745817664118026D+00, &
          0.0508364671763949D+00, &
          0.4745817664118026D+00, &
          0.2785739189962955D+00, &
          0.4428521620074090D+00, &
          0.2785739189962955D+00, &
          0.2785739189962955D+00, &
          0.4428521620074090D+00, &
          0.2785739189962955D+00, &
          0.1564426776295505D+00, &
          0.6871146447408989D+00, &
          0.1564426776295505D+00, &
          0.1564426776295505D+00, &
          0.6871146447408989D+00, &
          0.1564426776295505D+00, &
          0.4933597406415689D+00, &
          0.0132805187168622D+00, &
          0.4933597406415689D+00, &
          0.4933597406415689D+00, &
          0.0132805187168622D+00, &
          0.4933597406415689D+00, &
          0.0993149454903146D+00, &
          0.5254267651000835D+00, &
          0.3752582894096019D+00, &
          0.5254267651000835D+00, &
          0.3752582894096019D+00, &
          0.0993149454903146D+00, &
          0.0212715561472286D+00, &
          0.3482743570815052D+00, &
          0.6304540867712661D+00, &
          0.3482743570815052D+00, &
          0.6304540867712661D+00, &
          0.0212715561472286D+00, &
          0.0860850979983652D+00, &
          0.1654756461617898D+00, &
          0.7484392558398449D+00, &
          0.1654756461617898D+00, &
          0.7484392558398449D+00, &
          0.0860850979983652D+00, &
          0.0551722290870520D+00, &
          0.1639763763638563D+00, &
          0.7808513945490917D+00, &
          0.1639763763638563D+00, &
          0.7808513945490917D+00, &
          0.0551722290870520D+00, &
          0.0551722290870520D+00, &
          0.1639763763638563D+00, &
          0.7808513945490917D+00, &
          0.1639763763638563D+00, &
          0.7808513945490917D+00, &
          0.0551722290870520D+00, &
          0.0105036665914126D+00, &
          0.0370382522819345D+00, &
          0.9524580811266530D+00, &
          0.0370382522819345D+00, &
          0.9524580811266530D+00, &
          0.0105036665914126D+00, &
          0.0105036665914126D+00, &
          0.0370382522819345D+00, &
          0.9524580811266530D+00, &
          0.0370382522819345D+00, &
          0.9524580811266530D+00, &
          0.0105036665914126D+00, &
          0.0055639594221052D+00, &
          0.6350080118967869D+00, &
          0.3594280286811078D+00, &
          0.6350080118967869D+00, &
          0.3594280286811078D+00, &
          0.0055639594221052D+00, &
          0.0055639594221052D+00, &
          0.6350080118967869D+00, &
          0.3594280286811078D+00, &
          0.6350080118967869D+00, &
          0.3594280286811078D+00, &
          0.0055639594221052D+00, &
          0.0724868197207572D+00, &
          0.1297894012222233D+00, &
          0.7977237790570195D+00, &
          0.1297894012222233D+00, &
          0.7977237790570195D+00, &
          0.0724868197207572D+00, &
          0.0724868197207572D+00, &
          0.1297894012222233D+00, &
          0.7977237790570195D+00, &
          0.1297894012222233D+00, &
          0.7977237790570195D+00, &
          0.0724868197207572D+00, &
          0.0621201882612809D+00, &
          0.3059321030072774D+00, &
          0.6319477087314417D+00, &
          0.3059321030072774D+00, &
          0.6319477087314417D+00, &
          0.0621201882612809D+00, &
          0.0621201882612809D+00, &
          0.3059321030072774D+00, &
          0.6319477087314417D+00, &
          0.3059321030072774D+00, &
          0.6319477087314417D+00, &
          0.0621201882612809D+00, &
          0.2612303617350463D+00, &
          0.1677554820715698D+00, &
          0.5710141561933839D+00, &
          0.1677554820715698D+00, &
          0.5710141561933839D+00, &
          0.2612303617350463D+00, &
          0.2612303617350463D+00, &
          0.1677554820715698D+00, &
          0.5710141561933839D+00, &
          0.1677554820715698D+00, &
          0.5710141561933839D+00, &
          0.2612303617350463D+00, &
          0.0337686922990924D+00, &
          0.0015296219296330D+00, &
          0.9647016857712746D+00, &
          0.0015296219296330D+00, &
          0.9647016857712746D+00, &
          0.0337686922990924D+00, &
          0.0337686922990924D+00, &
          0.0015296219296330D+00, &
          0.9647016857712746D+00, &
          0.0015296219296330D+00, &
          0.9647016857712746D+00, &
          0.0337686922990924D+00, &
          0.0949242221588698D+00, &
          0.3532796606132861D+00, &
          0.5517961172278441D+00, &
          0.3532796606132861D+00, &
          0.5517961172278441D+00, &
          0.0949242221588698D+00, &
          0.0949242221588698D+00, &
          0.3532796606132861D+00, &
          0.5517961172278441D+00, &
          0.3532796606132861D+00, &
          0.5517961172278441D+00, &
          0.0949242221588698D+00, &
          0.0248817170289358D+00, &
          0.0959262325920420D+00, &
          0.8791920503790223D+00, &
          0.0959262325920420D+00, &
          0.8791920503790223D+00, &
          0.0248817170289358D+00, &
          0.0248817170289358D+00, &
          0.0959262325920420D+00, &
          0.8791920503790223D+00, &
          0.0959262325920420D+00, &
          0.8791920503790223D+00, &
          0.0248817170289358D+00, &
          0.1757029526928022D+00, &
          0.4963305263335406D+00, &
          0.3279665209736572D+00, &
          0.4963305263335406D+00, &
          0.3279665209736572D+00, &
          0.1757029526928022D+00, &
          0.1757029526928022D+00, &
          0.4963305263335406D+00, &
          0.3279665209736572D+00, &
          0.4963305263335406D+00, &
          0.3279665209736572D+00, &
          0.1757029526928022D+00, &
          0.0140786122092260D+00, &
          0.2738501035919119D+00, &
          0.7120712841988621D+00, &
          0.2738501035919119D+00, &
          0.7120712841988621D+00, &
          0.0140786122092260D+00, &
          0.0140786122092260D+00, &
          0.2738501035919119D+00, &
          0.7120712841988621D+00, &
          0.2738501035919119D+00, &
          0.7120712841988621D+00, &
          0.0140786122092260D+00, &
          0.0565171623200343D+00, &
          0.3828844998895772D+00, &
          0.5605983377903885D+00, &
          0.3828844998895772D+00, &
          0.5605983377903885D+00, &
          0.0565171623200343D+00, &
          0.0565171623200343D+00, &
          0.3828844998895772D+00, &
          0.5605983377903885D+00, &
          0.3828844998895772D+00, &
          0.5605983377903885D+00, &
          0.0565171623200343D+00, &
          0.0192412570201629D+00, &
          0.0899365323953004D+00, &
          0.8908222105845367D+00, &
          0.0899365323953004D+00, &
          0.8908222105845367D+00, &
          0.0192412570201629D+00, &
          0.0192412570201629D+00, &
          0.0899365323953004D+00, &
          0.8908222105845367D+00, &
          0.0899365323953004D+00, &
          0.8908222105845367D+00, &
          0.0192412570201629D+00, &
          0.0611370184542401D+00, &
          0.2516074150213084D+00, &
          0.6872555665244514D+00, &
          0.2516074150213084D+00, &
          0.6872555665244514D+00, &
          0.0611370184542401D+00, &
          0.0611370184542401D+00, &
          0.2516074150213084D+00, &
          0.6872555665244514D+00, &
          0.2516074150213084D+00, &
          0.6872555665244514D+00, &
          0.0611370184542401D+00, &
          0.0050182699561948D+00, &
          0.1524244609819593D+00, &
          0.8425572690618459D+00, &
          0.1524244609819593D+00, &
          0.8425572690618459D+00, &
          0.0050182699561948D+00, &
          0.0050182699561948D+00, &
          0.1524244609819593D+00, &
          0.8425572690618459D+00, &
          0.1524244609819593D+00, &
          0.8425572690618459D+00, &
          0.0050182699561948D+00, &
          0.1365910574355268D+00, &
          0.2475011786624472D+00, &
          0.6159077639020261D+00, &
          0.2475011786624472D+00, &
          0.6159077639020261D+00, &
          0.1365910574355268D+00, &
          0.1365910574355268D+00, &
          0.2475011786624472D+00, &
          0.6159077639020261D+00, &
          0.2475011786624472D+00, &
          0.6159077639020261D+00, &
          0.1365910574355268D+00, &
          0.0173543718109579D+00, &
          0.7591852986956938D+00, &
          0.2234603294933483D+00, &
          0.7591852986956938D+00, &
          0.2234603294933483D+00, &
          0.0173543718109579D+00, &
          0.0173543718109579D+00, &
          0.7591852986956938D+00, &
          0.2234603294933483D+00, &
          0.7591852986956938D+00, &
          0.2234603294933483D+00, &
          0.0173543718109579D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.8403117926958672D+00, &
          0.1596882073041329D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.8842725336194912D+00, &
          0.8842725336194912D+00, &
          0.8842725336194912D+00, &
          0.1157274663805087D+00, &
          0.1157274663805087D+00, &
          0.1157274663805087D+00, &
          0.6593583359504520D+00, &
          0.6593583359504520D+00, &
          0.6593583359504520D+00, &
          0.3406416640495479D+00, &
          0.3406416640495479D+00, &
          0.3406416640495479D+00, &
          0.8712606670369063D+00, &
          0.8712606670369063D+00, &
          0.8712606670369063D+00, &
          0.1287393329630936D+00, &
          0.1287393329630936D+00, &
          0.1287393329630936D+00, &
          0.2304229580101690D+00, &
          0.2304229580101690D+00, &
          0.2304229580101690D+00, &
          0.7695770419898310D+00, &
          0.7695770419898310D+00, &
          0.7695770419898310D+00, &
          0.6455038513606328D+00, &
          0.6455038513606328D+00, &
          0.6455038513606328D+00, &
          0.3544961486393671D+00, &
          0.3544961486393671D+00, &
          0.3544961486393671D+00, &
          0.8851421910096543D+00, &
          0.8851421910096543D+00, &
          0.8851421910096543D+00, &
          0.1148578089903457D+00, &
          0.1148578089903457D+00, &
          0.1148578089903457D+00, &
          0.9818829756657210D+00, &
          0.9818829756657210D+00, &
          0.9818829756657210D+00, &
          0.0181170243342790D+00, &
          0.0181170243342790D+00, &
          0.0181170243342790D+00, &
          0.9850118189004804D+00, &
          0.9850118189004804D+00, &
          0.9850118189004804D+00, &
          0.0149881810995195D+00, &
          0.0149881810995195D+00, &
          0.0149881810995195D+00, &
          0.9701315317214787D+00, &
          0.9701315317214787D+00, &
          0.9701315317214787D+00, &
          0.0298684682785212D+00, &
          0.0298684682785212D+00, &
          0.0298684682785212D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9387595239709540D+00, &
          0.9387595239709540D+00, &
          0.9387595239709540D+00, &
          0.9387595239709540D+00, &
          0.9387595239709540D+00, &
          0.9387595239709540D+00, &
          0.0612404760290460D+00, &
          0.0612404760290460D+00, &
          0.0612404760290460D+00, &
          0.0612404760290460D+00, &
          0.0612404760290460D+00, &
          0.0612404760290460D+00, &
          0.7620680922209941D+00, &
          0.7620680922209941D+00, &
          0.7620680922209941D+00, &
          0.7620680922209941D+00, &
          0.7620680922209941D+00, &
          0.7620680922209941D+00, &
          0.2379319077790059D+00, &
          0.2379319077790059D+00, &
          0.2379319077790059D+00, &
          0.2379319077790059D+00, &
          0.2379319077790059D+00, &
          0.2379319077790059D+00, &
          0.8273450356745347D+00, &
          0.8273450356745347D+00, &
          0.8273450356745347D+00, &
          0.8273450356745347D+00, &
          0.8273450356745347D+00, &
          0.8273450356745347D+00, &
          0.1726549643254653D+00, &
          0.1726549643254653D+00, &
          0.1726549643254653D+00, &
          0.1726549643254653D+00, &
          0.1726549643254653D+00, &
          0.1726549643254653D+00, &
          0.7401886228316812D+00, &
          0.7401886228316812D+00, &
          0.7401886228316812D+00, &
          0.7401886228316812D+00, &
          0.7401886228316812D+00, &
          0.7401886228316812D+00, &
          0.2598113771683188D+00, &
          0.2598113771683188D+00, &
          0.2598113771683188D+00, &
          0.2598113771683188D+00, &
          0.2598113771683188D+00, &
          0.2598113771683188D+00, &
          0.9210357305129615D+00, &
          0.9210357305129615D+00, &
          0.9210357305129615D+00, &
          0.9210357305129615D+00, &
          0.9210357305129615D+00, &
          0.9210357305129615D+00, &
          0.0789642694870386D+00, &
          0.0789642694870386D+00, &
          0.0789642694870386D+00, &
          0.0789642694870386D+00, &
          0.0789642694870386D+00, &
          0.0789642694870386D+00, &
          0.8233326926574609D+00, &
          0.8233326926574609D+00, &
          0.8233326926574609D+00, &
          0.8233326926574609D+00, &
          0.8233326926574609D+00, &
          0.8233326926574609D+00, &
          0.1766673073425390D+00, &
          0.1766673073425390D+00, &
          0.1766673073425390D+00, &
          0.1766673073425390D+00, &
          0.1766673073425390D+00, &
          0.1766673073425390D+00, &
          0.9478854137232893D+00, &
          0.9478854137232893D+00, &
          0.9478854137232893D+00, &
          0.9478854137232893D+00, &
          0.9478854137232893D+00, &
          0.9478854137232893D+00, &
          0.0521145862767107D+00, &
          0.0521145862767107D+00, &
          0.0521145862767107D+00, &
          0.0521145862767107D+00, &
          0.0521145862767107D+00, &
          0.0521145862767107D+00, &
          0.9916391512107949D+00, &
          0.9916391512107949D+00, &
          0.9916391512107949D+00, &
          0.9916391512107949D+00, &
          0.9916391512107949D+00, &
          0.9916391512107949D+00, &
          0.0083608487892050D+00, &
          0.0083608487892050D+00, &
          0.0083608487892050D+00, &
          0.0083608487892050D+00, &
          0.0083608487892050D+00, &
          0.0083608487892050D+00, &
          0.9898679390727328D+00, &
          0.9898679390727328D+00, &
          0.9898679390727328D+00, &
          0.9898679390727328D+00, &
          0.9898679390727328D+00, &
          0.9898679390727328D+00, &
          0.0101320609272672D+00, &
          0.0101320609272672D+00, &
          0.0101320609272672D+00, &
          0.0101320609272672D+00, &
          0.0101320609272672D+00, &
          0.0101320609272672D+00, &
          0.9305712292308971D+00, &
          0.9305712292308971D+00, &
          0.9305712292308971D+00, &
          0.9305712292308971D+00, &
          0.9305712292308971D+00, &
          0.9305712292308971D+00, &
          0.0694287707691029D+00, &
          0.0694287707691029D+00, &
          0.0694287707691029D+00, &
          0.0694287707691029D+00, &
          0.0694287707691029D+00, &
          0.0694287707691029D+00, &
          0.9758011913142881D+00, &
          0.9758011913142881D+00, &
          0.9758011913142881D+00, &
          0.9758011913142881D+00, &
          0.9758011913142881D+00, &
          0.9758011913142881D+00, &
          0.0241988086857118D+00, &
          0.0241988086857118D+00, &
          0.0241988086857118D+00, &
          0.0241988086857118D+00, &
          0.0241988086857118D+00, &
          0.0241988086857118D+00, &
          0.3036894659036268D+00, &
          0.3036894659036268D+00, &
          0.3036894659036268D+00, &
          0.3036894659036268D+00, &
          0.3036894659036268D+00, &
          0.3036894659036268D+00, &
          0.6963105340963731D+00, &
          0.6963105340963731D+00, &
          0.6963105340963731D+00, &
          0.6963105340963731D+00, &
          0.6963105340963731D+00, &
          0.6963105340963731D+00, &
          0.6138684175664644D+00, &
          0.6138684175664644D+00, &
          0.6138684175664644D+00, &
          0.6138684175664644D+00, &
          0.6138684175664644D+00, &
          0.6138684175664644D+00, &
          0.3861315824335356D+00, &
          0.3861315824335356D+00, &
          0.3861315824335356D+00, &
          0.3861315824335356D+00, &
          0.3861315824335356D+00, &
          0.3861315824335356D+00, &
          0.7885638236910046D+00, &
          0.7885638236910046D+00, &
          0.7885638236910046D+00, &
          0.7885638236910046D+00, &
          0.7885638236910046D+00, &
          0.7885638236910046D+00, &
          0.2114361763089954D+00, &
          0.2114361763089954D+00, &
          0.2114361763089954D+00, &
          0.2114361763089954D+00, &
          0.2114361763089954D+00, &
          0.2114361763089954D+00, &
          0.8514618541948253D+00, &
          0.8514618541948253D+00, &
          0.8514618541948253D+00, &
          0.8514618541948253D+00, &
          0.8514618541948253D+00, &
          0.8514618541948253D+00, &
          0.1485381458051747D+00, &
          0.1485381458051747D+00, &
          0.1485381458051747D+00, &
          0.1485381458051747D+00, &
          0.1485381458051747D+00, &
          0.1485381458051747D+00, &
          0.3722497735842279D+00, &
          0.3722497735842279D+00, &
          0.3722497735842279D+00, &
          0.3722497735842279D+00, &
          0.3722497735842279D+00, &
          0.3722497735842279D+00, &
          0.6277502264157721D+00, &
          0.6277502264157721D+00, &
          0.6277502264157721D+00, &
          0.6277502264157721D+00, &
          0.6277502264157721D+00, &
          0.6277502264157721D+00, &
          0.6329858368666343D+00, &
          0.6329858368666343D+00, &
          0.6329858368666343D+00, &
          0.6329858368666343D+00, &
          0.6329858368666343D+00, &
          0.6329858368666343D+00, &
          0.3670141631333657D+00, &
          0.3670141631333657D+00, &
          0.3670141631333657D+00, &
          0.3670141631333657D+00, &
          0.3670141631333657D+00, &
          0.3670141631333657D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          0.0126082800938750D+00, &
          0.0126082800938750D+00, &
          0.0005804416586274D+00, &
          0.0005804416586274D+00, &
          0.0005804416586274D+00, &
          0.0014515275277035D+00, &
          0.0014515275277035D+00, &
          0.0014515275277035D+00, &
          0.0110392519457591D+00, &
          0.0110392519457591D+00, &
          0.0110392519457591D+00, &
          0.0023896607884456D+00, &
          0.0023896607884456D+00, &
          0.0023896607884456D+00, &
          0.0023896607884456D+00, &
          0.0023896607884456D+00, &
          0.0023896607884456D+00, &
          0.0111852959227564D+00, &
          0.0111852959227564D+00, &
          0.0111852959227564D+00, &
          0.0111852959227564D+00, &
          0.0111852959227564D+00, &
          0.0111852959227564D+00, &
          0.0044527106261489D+00, &
          0.0044527106261489D+00, &
          0.0044527106261489D+00, &
          0.0044527106261489D+00, &
          0.0044527106261489D+00, &
          0.0044527106261489D+00, &
          0.0099310429205233D+00, &
          0.0099310429205233D+00, &
          0.0099310429205233D+00, &
          0.0099310429205233D+00, &
          0.0099310429205233D+00, &
          0.0099310429205233D+00, &
          0.0027983315098658D+00, &
          0.0027983315098658D+00, &
          0.0027983315098658D+00, &
          0.0027983315098658D+00, &
          0.0027983315098658D+00, &
          0.0027983315098658D+00, &
          0.0045950791440008D+00, &
          0.0045950791440008D+00, &
          0.0045950791440008D+00, &
          0.0045950791440008D+00, &
          0.0045950791440008D+00, &
          0.0045950791440008D+00, &
          0.0026497419827358D+00, &
          0.0026497419827358D+00, &
          0.0026497419827358D+00, &
          0.0026497419827358D+00, &
          0.0026497419827358D+00, &
          0.0026497419827358D+00, &
          0.0023791749450044D+00, &
          0.0023791749450044D+00, &
          0.0023791749450044D+00, &
          0.0023791749450044D+00, &
          0.0023791749450044D+00, &
          0.0023791749450044D+00, &
          0.0012902115379941D+00, &
          0.0012902115379941D+00, &
          0.0012902115379941D+00, &
          0.0012902115379941D+00, &
          0.0012902115379941D+00, &
          0.0012902115379941D+00, &
          0.0071749704457678D+00, &
          0.0071749704457678D+00, &
          0.0071749704457678D+00, &
          0.0071749704457678D+00, &
          0.0071749704457678D+00, &
          0.0071749704457678D+00, &
          0.0039041969129065D+00, &
          0.0039041969129065D+00, &
          0.0039041969129065D+00, &
          0.0039041969129065D+00, &
          0.0039041969129065D+00, &
          0.0039041969129065D+00, &
          0.0058968479213494D+00, &
          0.0058968479213494D+00, &
          0.0058968479213494D+00, &
          0.0058968479213494D+00, &
          0.0058968479213494D+00, &
          0.0058968479213494D+00, &
          0.0024349402323106D+00, &
          0.0024349402323106D+00, &
          0.0024349402323106D+00, &
          0.0024349402323106D+00, &
          0.0024349402323106D+00, &
          0.0024349402323106D+00, &
          0.0024349402323106D+00, &
          0.0024349402323106D+00, &
          0.0024349402323106D+00, &
          0.0024349402323106D+00, &
          0.0024349402323106D+00, &
          0.0024349402323106D+00, &
          0.0007867926281754D+00, &
          0.0007867926281754D+00, &
          0.0007867926281754D+00, &
          0.0007867926281754D+00, &
          0.0007867926281754D+00, &
          0.0007867926281754D+00, &
          0.0007867926281754D+00, &
          0.0007867926281754D+00, &
          0.0007867926281754D+00, &
          0.0007867926281754D+00, &
          0.0007867926281754D+00, &
          0.0007867926281754D+00, &
          0.0019733458950933D+00, &
          0.0019733458950933D+00, &
          0.0019733458950933D+00, &
          0.0019733458950933D+00, &
          0.0019733458950933D+00, &
          0.0019733458950933D+00, &
          0.0019733458950933D+00, &
          0.0019733458950933D+00, &
          0.0019733458950933D+00, &
          0.0019733458950933D+00, &
          0.0019733458950933D+00, &
          0.0019733458950933D+00, &
          0.0039028044764867D+00, &
          0.0039028044764867D+00, &
          0.0039028044764867D+00, &
          0.0039028044764867D+00, &
          0.0039028044764867D+00, &
          0.0039028044764867D+00, &
          0.0039028044764867D+00, &
          0.0039028044764867D+00, &
          0.0039028044764867D+00, &
          0.0039028044764867D+00, &
          0.0039028044764867D+00, &
          0.0039028044764867D+00, &
          0.0028651948489933D+00, &
          0.0028651948489933D+00, &
          0.0028651948489933D+00, &
          0.0028651948489933D+00, &
          0.0028651948489933D+00, &
          0.0028651948489933D+00, &
          0.0028651948489933D+00, &
          0.0028651948489933D+00, &
          0.0028651948489933D+00, &
          0.0028651948489933D+00, &
          0.0028651948489933D+00, &
          0.0028651948489933D+00, &
          0.0042658629737291D+00, &
          0.0042658629737291D+00, &
          0.0042658629737291D+00, &
          0.0042658629737291D+00, &
          0.0042658629737291D+00, &
          0.0042658629737291D+00, &
          0.0042658629737291D+00, &
          0.0042658629737291D+00, &
          0.0042658629737291D+00, &
          0.0042658629737291D+00, &
          0.0042658629737291D+00, &
          0.0042658629737291D+00, &
          0.0002763846307136D+00, &
          0.0002763846307136D+00, &
          0.0002763846307136D+00, &
          0.0002763846307136D+00, &
          0.0002763846307136D+00, &
          0.0002763846307136D+00, &
          0.0002763846307136D+00, &
          0.0002763846307136D+00, &
          0.0002763846307136D+00, &
          0.0002763846307136D+00, &
          0.0002763846307136D+00, &
          0.0002763846307136D+00, &
          0.0015947204945141D+00, &
          0.0015947204945141D+00, &
          0.0015947204945141D+00, &
          0.0015947204945141D+00, &
          0.0015947204945141D+00, &
          0.0015947204945141D+00, &
          0.0015947204945141D+00, &
          0.0015947204945141D+00, &
          0.0015947204945141D+00, &
          0.0015947204945141D+00, &
          0.0015947204945141D+00, &
          0.0015947204945141D+00, &
          0.0005967018601525D+00, &
          0.0005967018601525D+00, &
          0.0005967018601525D+00, &
          0.0005967018601525D+00, &
          0.0005967018601525D+00, &
          0.0005967018601525D+00, &
          0.0005967018601525D+00, &
          0.0005967018601525D+00, &
          0.0005967018601525D+00, &
          0.0005967018601525D+00, &
          0.0005967018601525D+00, &
          0.0005967018601525D+00, &
          0.0050890108573681D+00, &
          0.0050890108573681D+00, &
          0.0050890108573681D+00, &
          0.0050890108573681D+00, &
          0.0050890108573681D+00, &
          0.0050890108573681D+00, &
          0.0050890108573681D+00, &
          0.0050890108573681D+00, &
          0.0050890108573681D+00, &
          0.0050890108573681D+00, &
          0.0050890108573681D+00, &
          0.0050890108573681D+00, &
          0.0010232527285653D+00, &
          0.0010232527285653D+00, &
          0.0010232527285653D+00, &
          0.0010232527285653D+00, &
          0.0010232527285653D+00, &
          0.0010232527285653D+00, &
          0.0010232527285653D+00, &
          0.0010232527285653D+00, &
          0.0010232527285653D+00, &
          0.0010232527285653D+00, &
          0.0010232527285653D+00, &
          0.0010232527285653D+00, &
          0.0048336118379361D+00, &
          0.0048336118379361D+00, &
          0.0048336118379361D+00, &
          0.0048336118379361D+00, &
          0.0048336118379361D+00, &
          0.0048336118379361D+00, &
          0.0048336118379361D+00, &
          0.0048336118379361D+00, &
          0.0048336118379361D+00, &
          0.0048336118379361D+00, &
          0.0048336118379361D+00, &
          0.0048336118379361D+00, &
          0.0022565704494986D+00, &
          0.0022565704494986D+00, &
          0.0022565704494986D+00, &
          0.0022565704494986D+00, &
          0.0022565704494986D+00, &
          0.0022565704494986D+00, &
          0.0022565704494986D+00, &
          0.0022565704494986D+00, &
          0.0022565704494986D+00, &
          0.0022565704494986D+00, &
          0.0022565704494986D+00, &
          0.0022565704494986D+00, &
          0.0043942174418955D+00, &
          0.0043942174418955D+00, &
          0.0043942174418955D+00, &
          0.0043942174418955D+00, &
          0.0043942174418955D+00, &
          0.0043942174418955D+00, &
          0.0043942174418955D+00, &
          0.0043942174418955D+00, &
          0.0043942174418955D+00, &
          0.0043942174418955D+00, &
          0.0043942174418955D+00, &
          0.0043942174418955D+00, &
          0.0014412804208529D+00, &
          0.0014412804208529D+00, &
          0.0014412804208529D+00, &
          0.0014412804208529D+00, &
          0.0014412804208529D+00, &
          0.0014412804208529D+00, &
          0.0014412804208529D+00, &
          0.0014412804208529D+00, &
          0.0014412804208529D+00, &
          0.0014412804208529D+00, &
          0.0014412804208529D+00, &
          0.0014412804208529D+00, &
          0.0076389038664824D+00, &
          0.0076389038664824D+00, &
          0.0076389038664824D+00, &
          0.0076389038664824D+00, &
          0.0076389038664824D+00, &
          0.0076389038664824D+00, &
          0.0076389038664824D+00, &
          0.0076389038664824D+00, &
          0.0076389038664824D+00, &
          0.0076389038664824D+00, &
          0.0076389038664824D+00, &
          0.0076389038664824D+00, &
          0.0032669200631480D+00, &
          0.0032669200631480D+00, &
          0.0032669200631480D+00, &
          0.0032669200631480D+00, &
          0.0032669200631480D+00, &
          0.0032669200631480D+00, &
          0.0032669200631480D+00, &
          0.0032669200631480D+00, &
          0.0032669200631480D+00, &
          0.0032669200631480D+00, &
          0.0032669200631480D+00, &
          0.0032669200631480D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule17 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule17() returns the prism rule of precision 17.
!
!  Discussion:
!
!    We suppose we are given a triangular prism P with vertices
!    (1,0,0), (0,1,0), (0,0,0),
!    (1,0,1), (0,1,1), (0,0,1),
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over P is 
!    approximated by Q as follows:
!
!    Q = volume(P) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 May 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order symmetric cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 122, Number 1, pages 148-171, 24 August 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real x(n), y(n), z(n): the coordinates of quadrature points.
!
!    real w(n): the quadrature weights, which sum to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: n_save = 338

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.4665556153120432D+00, &
          0.4665556153120432D+00, &
          0.0668887693759136D+00, &
          0.4665556153120432D+00, &
          0.4665556153120432D+00, &
          0.0668887693759136D+00, &
          0.4923809674098066D+00, &
          0.4923809674098066D+00, &
          0.0152380651803867D+00, &
          0.4923809674098066D+00, &
          0.4923809674098066D+00, &
          0.0152380651803867D+00, &
          0.4466963397682847D+00, &
          0.4466963397682847D+00, &
          0.1066073204634306D+00, &
          0.4466963397682847D+00, &
          0.4466963397682847D+00, &
          0.1066073204634306D+00, &
          0.2058899494758731D+00, &
          0.2058899494758731D+00, &
          0.5882201010482538D+00, &
          0.2058899494758731D+00, &
          0.2058899494758731D+00, &
          0.5882201010482538D+00, &
          0.2305779360323416D+00, &
          0.2305779360323416D+00, &
          0.5388441279353168D+00, &
          0.2305779360323416D+00, &
          0.2305779360323416D+00, &
          0.5388441279353168D+00, &
          0.0070726502342596D+00, &
          0.0070726502342596D+00, &
          0.9858546995314807D+00, &
          0.0070726502342596D+00, &
          0.0070726502342596D+00, &
          0.9858546995314807D+00, &
          0.0524838545504923D+00, &
          0.0524838545504923D+00, &
          0.8950322908990154D+00, &
          0.0524838545504923D+00, &
          0.0524838545504923D+00, &
          0.8950322908990154D+00, &
          0.2780743897936196D+00, &
          0.2780743897936196D+00, &
          0.4438512204127608D+00, &
          0.2780743897936196D+00, &
          0.2780743897936196D+00, &
          0.4438512204127608D+00, &
          0.3893995813144771D+00, &
          0.3893995813144771D+00, &
          0.2212008373710459D+00, &
          0.3893995813144771D+00, &
          0.3893995813144771D+00, &
          0.2212008373710459D+00, &
          0.1132999955833483D+00, &
          0.1132999955833483D+00, &
          0.7734000088333033D+00, &
          0.1132999955833483D+00, &
          0.1132999955833483D+00, &
          0.7734000088333033D+00, &
          0.4687774749283292D+00, &
          0.4687774749283292D+00, &
          0.0624450501433417D+00, &
          0.4687774749283292D+00, &
          0.4687774749283292D+00, &
          0.0624450501433417D+00, &
          0.4942337889964747D+00, &
          0.4942337889964747D+00, &
          0.0115324220070506D+00, &
          0.4942337889964747D+00, &
          0.4942337889964747D+00, &
          0.0115324220070506D+00, &
          0.2906999098818323D+00, &
          0.2906999098818323D+00, &
          0.4186001802363355D+00, &
          0.2906999098818323D+00, &
          0.2906999098818323D+00, &
          0.4186001802363355D+00, &
          0.0937693039918497D+00, &
          0.0937693039918497D+00, &
          0.8124613920163005D+00, &
          0.0937693039918497D+00, &
          0.0937693039918497D+00, &
          0.8124613920163005D+00, &
          0.1632812980104028D+00, &
          0.1632812980104028D+00, &
          0.6734374039791945D+00, &
          0.1632812980104028D+00, &
          0.1632812980104028D+00, &
          0.6734374039791945D+00, &
          0.1774626418569745D+00, &
          0.1774626418569745D+00, &
          0.6450747162860511D+00, &
          0.1774626418569745D+00, &
          0.1774626418569745D+00, &
          0.6450747162860511D+00, &
          0.3930786561208929D+00, &
          0.3930786561208929D+00, &
          0.2138426877582142D+00, &
          0.3930786561208929D+00, &
          0.3930786561208929D+00, &
          0.2138426877582142D+00, &
          0.2711153947647855D+00, &
          0.2711153947647855D+00, &
          0.4577692104704291D+00, &
          0.2711153947647855D+00, &
          0.2711153947647855D+00, &
          0.4577692104704291D+00, &
          0.0243062210015648D+00, &
          0.0243062210015648D+00, &
          0.9513875579968704D+00, &
          0.0243062210015648D+00, &
          0.0243062210015648D+00, &
          0.9513875579968704D+00, &
          0.1273859447314594D+00, &
          0.2765687474993688D+00, &
          0.1273859447314594D+00, &
          0.5960453077691718D+00, &
          0.2765687474993688D+00, &
          0.5960453077691718D+00, &
          0.2450510833431427D+00, &
          0.0852468501828055D+00, &
          0.2450510833431427D+00, &
          0.6697020664740517D+00, &
          0.0852468501828055D+00, &
          0.6697020664740517D+00, &
          0.2450510833431427D+00, &
          0.0852468501828055D+00, &
          0.2450510833431427D+00, &
          0.6697020664740517D+00, &
          0.0852468501828055D+00, &
          0.6697020664740517D+00, &
          0.3192911229101512D+00, &
          0.0276449852112552D+00, &
          0.3192911229101512D+00, &
          0.6530638918785936D+00, &
          0.0276449852112552D+00, &
          0.6530638918785936D+00, &
          0.3192911229101512D+00, &
          0.0276449852112552D+00, &
          0.3192911229101512D+00, &
          0.6530638918785936D+00, &
          0.0276449852112552D+00, &
          0.6530638918785936D+00, &
          0.4577741357714497D+00, &
          0.1385428807092584D+00, &
          0.4577741357714497D+00, &
          0.4036829835192919D+00, &
          0.1385428807092584D+00, &
          0.4036829835192919D+00, &
          0.4577741357714497D+00, &
          0.1385428807092584D+00, &
          0.4577741357714497D+00, &
          0.4036829835192919D+00, &
          0.1385428807092584D+00, &
          0.4036829835192919D+00, &
          0.2421105196372634D+00, &
          0.0124568391026713D+00, &
          0.2421105196372634D+00, &
          0.7454326412600653D+00, &
          0.0124568391026713D+00, &
          0.7454326412600653D+00, &
          0.2421105196372634D+00, &
          0.0124568391026713D+00, &
          0.2421105196372634D+00, &
          0.7454326412600653D+00, &
          0.0124568391026713D+00, &
          0.7454326412600653D+00, &
          0.7269945541994681D+00, &
          0.0888093024099239D+00, &
          0.7269945541994681D+00, &
          0.1841961433906080D+00, &
          0.0888093024099239D+00, &
          0.1841961433906080D+00, &
          0.7269945541994681D+00, &
          0.0888093024099239D+00, &
          0.7269945541994681D+00, &
          0.1841961433906080D+00, &
          0.0888093024099239D+00, &
          0.1841961433906080D+00, &
          0.4870517086755036D+00, &
          0.1725779658472295D+00, &
          0.4870517086755036D+00, &
          0.3403703254772670D+00, &
          0.1725779658472295D+00, &
          0.3403703254772670D+00, &
          0.4870517086755036D+00, &
          0.1725779658472295D+00, &
          0.4870517086755036D+00, &
          0.3403703254772670D+00, &
          0.1725779658472295D+00, &
          0.3403703254772670D+00, &
          0.1582123188291182D+00, &
          0.3108975273704328D+00, &
          0.1582123188291182D+00, &
          0.5308901538004490D+00, &
          0.3108975273704328D+00, &
          0.5308901538004490D+00, &
          0.1582123188291182D+00, &
          0.3108975273704328D+00, &
          0.1582123188291182D+00, &
          0.5308901538004490D+00, &
          0.3108975273704328D+00, &
          0.5308901538004490D+00, &
          0.5973833093262283D+00, &
          0.0432967856355362D+00, &
          0.5973833093262283D+00, &
          0.3593199050382355D+00, &
          0.0432967856355362D+00, &
          0.3593199050382355D+00, &
          0.5973833093262283D+00, &
          0.0432967856355362D+00, &
          0.5973833093262283D+00, &
          0.3593199050382355D+00, &
          0.0432967856355362D+00, &
          0.3593199050382355D+00, &
          0.1111929938785087D+00, &
          0.0368260832055792D+00, &
          0.1111929938785087D+00, &
          0.8519809229159122D+00, &
          0.0368260832055792D+00, &
          0.8519809229159122D+00, &
          0.1111929938785087D+00, &
          0.0368260832055792D+00, &
          0.1111929938785087D+00, &
          0.8519809229159122D+00, &
          0.0368260832055792D+00, &
          0.8519809229159122D+00, &
          0.0387865006290567D+00, &
          0.0070454102530681D+00, &
          0.0387865006290567D+00, &
          0.9541680891178752D+00, &
          0.0070454102530681D+00, &
          0.9541680891178752D+00, &
          0.0387865006290567D+00, &
          0.0070454102530681D+00, &
          0.0387865006290567D+00, &
          0.9541680891178752D+00, &
          0.0070454102530681D+00, &
          0.9541680891178752D+00, &
          0.1417441449927015D+00, &
          0.0066561350058640D+00, &
          0.1417441449927015D+00, &
          0.8515997200014345D+00, &
          0.0066561350058640D+00, &
          0.8515997200014345D+00, &
          0.1417441449927015D+00, &
          0.0066561350058640D+00, &
          0.1417441449927015D+00, &
          0.8515997200014345D+00, &
          0.0066561350058640D+00, &
          0.8515997200014345D+00, &
          0.1343472140817544D+00, &
          0.0269814183822196D+00, &
          0.1343472140817544D+00, &
          0.8386713675360260D+00, &
          0.0269814183822196D+00, &
          0.8386713675360260D+00, &
          0.1343472140817544D+00, &
          0.0269814183822196D+00, &
          0.1343472140817544D+00, &
          0.8386713675360260D+00, &
          0.0269814183822196D+00, &
          0.8386713675360260D+00, &
          0.3195031824862086D+00, &
          0.0813440561372205D+00, &
          0.3195031824862086D+00, &
          0.5991527613765709D+00, &
          0.0813440561372205D+00, &
          0.5991527613765709D+00, &
          0.3195031824862086D+00, &
          0.0813440561372205D+00, &
          0.3195031824862086D+00, &
          0.5991527613765709D+00, &
          0.0813440561372205D+00, &
          0.5991527613765709D+00, &
          0.0826308438588340D+00, &
          0.0102635465039764D+00, &
          0.0826308438588340D+00, &
          0.9071056096371896D+00, &
          0.0102635465039764D+00, &
          0.9071056096371896D+00, &
          0.0826308438588340D+00, &
          0.0102635465039764D+00, &
          0.0826308438588340D+00, &
          0.9071056096371896D+00, &
          0.0102635465039764D+00, &
          0.9071056096371896D+00, &
          0.3827807429336583D+00, &
          0.0067490134150570D+00, &
          0.3827807429336583D+00, &
          0.6104702436512848D+00, &
          0.0067490134150570D+00, &
          0.6104702436512848D+00, &
          0.3827807429336583D+00, &
          0.0067490134150570D+00, &
          0.3827807429336583D+00, &
          0.6104702436512848D+00, &
          0.0067490134150570D+00, &
          0.6104702436512848D+00, &
          0.2446543644502634D+00, &
          0.0056501992086932D+00, &
          0.2446543644502634D+00, &
          0.7496954363410434D+00, &
          0.0056501992086932D+00, &
          0.7496954363410434D+00, &
          0.2446543644502634D+00, &
          0.0056501992086932D+00, &
          0.2446543644502634D+00, &
          0.7496954363410434D+00, &
          0.0056501992086932D+00, &
          0.7496954363410434D+00, &
          0.0263707136798594D+00, &
          0.6264598409479015D+00, &
          0.0263707136798594D+00, &
          0.3471694453722390D+00, &
          0.6264598409479015D+00, &
          0.3471694453722390D+00, &
          0.0263707136798594D+00, &
          0.6264598409479015D+00, &
          0.0263707136798594D+00, &
          0.3471694453722390D+00, &
          0.6264598409479015D+00, &
          0.3471694453722390D+00, &
          0.0434074542015924D+00, &
          0.1947680793017931D+00, &
          0.0434074542015924D+00, &
          0.7618244664966145D+00, &
          0.1947680793017931D+00, &
          0.7618244664966145D+00, &
          0.0434074542015924D+00, &
          0.1947680793017931D+00, &
          0.0434074542015924D+00, &
          0.7618244664966145D+00, &
          0.1947680793017931D+00, &
          0.7618244664966145D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.4665556153120432D+00, &
          0.0668887693759136D+00, &
          0.4665556153120432D+00, &
          0.4665556153120432D+00, &
          0.0668887693759136D+00, &
          0.4665556153120432D+00, &
          0.4923809674098066D+00, &
          0.0152380651803867D+00, &
          0.4923809674098066D+00, &
          0.4923809674098066D+00, &
          0.0152380651803867D+00, &
          0.4923809674098066D+00, &
          0.4466963397682847D+00, &
          0.1066073204634306D+00, &
          0.4466963397682847D+00, &
          0.4466963397682847D+00, &
          0.1066073204634306D+00, &
          0.4466963397682847D+00, &
          0.2058899494758731D+00, &
          0.5882201010482538D+00, &
          0.2058899494758731D+00, &
          0.2058899494758731D+00, &
          0.5882201010482538D+00, &
          0.2058899494758731D+00, &
          0.2305779360323416D+00, &
          0.5388441279353168D+00, &
          0.2305779360323416D+00, &
          0.2305779360323416D+00, &
          0.5388441279353168D+00, &
          0.2305779360323416D+00, &
          0.0070726502342596D+00, &
          0.9858546995314807D+00, &
          0.0070726502342596D+00, &
          0.0070726502342596D+00, &
          0.9858546995314807D+00, &
          0.0070726502342596D+00, &
          0.0524838545504923D+00, &
          0.8950322908990154D+00, &
          0.0524838545504923D+00, &
          0.0524838545504923D+00, &
          0.8950322908990154D+00, &
          0.0524838545504923D+00, &
          0.2780743897936196D+00, &
          0.4438512204127608D+00, &
          0.2780743897936196D+00, &
          0.2780743897936196D+00, &
          0.4438512204127608D+00, &
          0.2780743897936196D+00, &
          0.3893995813144771D+00, &
          0.2212008373710459D+00, &
          0.3893995813144771D+00, &
          0.3893995813144771D+00, &
          0.2212008373710459D+00, &
          0.3893995813144771D+00, &
          0.1132999955833483D+00, &
          0.7734000088333033D+00, &
          0.1132999955833483D+00, &
          0.1132999955833483D+00, &
          0.7734000088333033D+00, &
          0.1132999955833483D+00, &
          0.4687774749283292D+00, &
          0.0624450501433417D+00, &
          0.4687774749283292D+00, &
          0.4687774749283292D+00, &
          0.0624450501433417D+00, &
          0.4687774749283292D+00, &
          0.4942337889964747D+00, &
          0.0115324220070506D+00, &
          0.4942337889964747D+00, &
          0.4942337889964747D+00, &
          0.0115324220070506D+00, &
          0.4942337889964747D+00, &
          0.2906999098818323D+00, &
          0.4186001802363355D+00, &
          0.2906999098818323D+00, &
          0.2906999098818323D+00, &
          0.4186001802363355D+00, &
          0.2906999098818323D+00, &
          0.0937693039918497D+00, &
          0.8124613920163005D+00, &
          0.0937693039918497D+00, &
          0.0937693039918497D+00, &
          0.8124613920163005D+00, &
          0.0937693039918497D+00, &
          0.1632812980104028D+00, &
          0.6734374039791945D+00, &
          0.1632812980104028D+00, &
          0.1632812980104028D+00, &
          0.6734374039791945D+00, &
          0.1632812980104028D+00, &
          0.1774626418569745D+00, &
          0.6450747162860511D+00, &
          0.1774626418569745D+00, &
          0.1774626418569745D+00, &
          0.6450747162860511D+00, &
          0.1774626418569745D+00, &
          0.3930786561208929D+00, &
          0.2138426877582142D+00, &
          0.3930786561208929D+00, &
          0.3930786561208929D+00, &
          0.2138426877582142D+00, &
          0.3930786561208929D+00, &
          0.2711153947647855D+00, &
          0.4577692104704291D+00, &
          0.2711153947647855D+00, &
          0.2711153947647855D+00, &
          0.4577692104704291D+00, &
          0.2711153947647855D+00, &
          0.0243062210015648D+00, &
          0.9513875579968704D+00, &
          0.0243062210015648D+00, &
          0.0243062210015648D+00, &
          0.9513875579968704D+00, &
          0.0243062210015648D+00, &
          0.2765687474993688D+00, &
          0.1273859447314594D+00, &
          0.5960453077691718D+00, &
          0.1273859447314594D+00, &
          0.5960453077691718D+00, &
          0.2765687474993688D+00, &
          0.0852468501828055D+00, &
          0.2450510833431427D+00, &
          0.6697020664740517D+00, &
          0.2450510833431427D+00, &
          0.6697020664740517D+00, &
          0.0852468501828055D+00, &
          0.0852468501828055D+00, &
          0.2450510833431427D+00, &
          0.6697020664740517D+00, &
          0.2450510833431427D+00, &
          0.6697020664740517D+00, &
          0.0852468501828055D+00, &
          0.0276449852112552D+00, &
          0.3192911229101512D+00, &
          0.6530638918785936D+00, &
          0.3192911229101512D+00, &
          0.6530638918785936D+00, &
          0.0276449852112552D+00, &
          0.0276449852112552D+00, &
          0.3192911229101512D+00, &
          0.6530638918785936D+00, &
          0.3192911229101512D+00, &
          0.6530638918785936D+00, &
          0.0276449852112552D+00, &
          0.1385428807092584D+00, &
          0.4577741357714497D+00, &
          0.4036829835192919D+00, &
          0.4577741357714497D+00, &
          0.4036829835192919D+00, &
          0.1385428807092584D+00, &
          0.1385428807092584D+00, &
          0.4577741357714497D+00, &
          0.4036829835192919D+00, &
          0.4577741357714497D+00, &
          0.4036829835192919D+00, &
          0.1385428807092584D+00, &
          0.0124568391026713D+00, &
          0.2421105196372634D+00, &
          0.7454326412600653D+00, &
          0.2421105196372634D+00, &
          0.7454326412600653D+00, &
          0.0124568391026713D+00, &
          0.0124568391026713D+00, &
          0.2421105196372634D+00, &
          0.7454326412600653D+00, &
          0.2421105196372634D+00, &
          0.7454326412600653D+00, &
          0.0124568391026713D+00, &
          0.0888093024099239D+00, &
          0.7269945541994681D+00, &
          0.1841961433906080D+00, &
          0.7269945541994681D+00, &
          0.1841961433906080D+00, &
          0.0888093024099239D+00, &
          0.0888093024099239D+00, &
          0.7269945541994681D+00, &
          0.1841961433906080D+00, &
          0.7269945541994681D+00, &
          0.1841961433906080D+00, &
          0.0888093024099239D+00, &
          0.1725779658472295D+00, &
          0.4870517086755036D+00, &
          0.3403703254772670D+00, &
          0.4870517086755036D+00, &
          0.3403703254772670D+00, &
          0.1725779658472295D+00, &
          0.1725779658472295D+00, &
          0.4870517086755036D+00, &
          0.3403703254772670D+00, &
          0.4870517086755036D+00, &
          0.3403703254772670D+00, &
          0.1725779658472295D+00, &
          0.3108975273704328D+00, &
          0.1582123188291182D+00, &
          0.5308901538004490D+00, &
          0.1582123188291182D+00, &
          0.5308901538004490D+00, &
          0.3108975273704328D+00, &
          0.3108975273704328D+00, &
          0.1582123188291182D+00, &
          0.5308901538004490D+00, &
          0.1582123188291182D+00, &
          0.5308901538004490D+00, &
          0.3108975273704328D+00, &
          0.0432967856355362D+00, &
          0.5973833093262283D+00, &
          0.3593199050382355D+00, &
          0.5973833093262283D+00, &
          0.3593199050382355D+00, &
          0.0432967856355362D+00, &
          0.0432967856355362D+00, &
          0.5973833093262283D+00, &
          0.3593199050382355D+00, &
          0.5973833093262283D+00, &
          0.3593199050382355D+00, &
          0.0432967856355362D+00, &
          0.0368260832055792D+00, &
          0.1111929938785087D+00, &
          0.8519809229159122D+00, &
          0.1111929938785087D+00, &
          0.8519809229159122D+00, &
          0.0368260832055792D+00, &
          0.0368260832055792D+00, &
          0.1111929938785087D+00, &
          0.8519809229159122D+00, &
          0.1111929938785087D+00, &
          0.8519809229159122D+00, &
          0.0368260832055792D+00, &
          0.0070454102530681D+00, &
          0.0387865006290567D+00, &
          0.9541680891178752D+00, &
          0.0387865006290567D+00, &
          0.9541680891178752D+00, &
          0.0070454102530681D+00, &
          0.0070454102530681D+00, &
          0.0387865006290567D+00, &
          0.9541680891178752D+00, &
          0.0387865006290567D+00, &
          0.9541680891178752D+00, &
          0.0070454102530681D+00, &
          0.0066561350058640D+00, &
          0.1417441449927015D+00, &
          0.8515997200014345D+00, &
          0.1417441449927015D+00, &
          0.8515997200014345D+00, &
          0.0066561350058640D+00, &
          0.0066561350058640D+00, &
          0.1417441449927015D+00, &
          0.8515997200014345D+00, &
          0.1417441449927015D+00, &
          0.8515997200014345D+00, &
          0.0066561350058640D+00, &
          0.0269814183822196D+00, &
          0.1343472140817544D+00, &
          0.8386713675360260D+00, &
          0.1343472140817544D+00, &
          0.8386713675360260D+00, &
          0.0269814183822196D+00, &
          0.0269814183822196D+00, &
          0.1343472140817544D+00, &
          0.8386713675360260D+00, &
          0.1343472140817544D+00, &
          0.8386713675360260D+00, &
          0.0269814183822196D+00, &
          0.0813440561372205D+00, &
          0.3195031824862086D+00, &
          0.5991527613765709D+00, &
          0.3195031824862086D+00, &
          0.5991527613765709D+00, &
          0.0813440561372205D+00, &
          0.0813440561372205D+00, &
          0.3195031824862086D+00, &
          0.5991527613765709D+00, &
          0.3195031824862086D+00, &
          0.5991527613765709D+00, &
          0.0813440561372205D+00, &
          0.0102635465039764D+00, &
          0.0826308438588340D+00, &
          0.9071056096371896D+00, &
          0.0826308438588340D+00, &
          0.9071056096371896D+00, &
          0.0102635465039764D+00, &
          0.0102635465039764D+00, &
          0.0826308438588340D+00, &
          0.9071056096371896D+00, &
          0.0826308438588340D+00, &
          0.9071056096371896D+00, &
          0.0102635465039764D+00, &
          0.0067490134150570D+00, &
          0.3827807429336583D+00, &
          0.6104702436512848D+00, &
          0.3827807429336583D+00, &
          0.6104702436512848D+00, &
          0.0067490134150570D+00, &
          0.0067490134150570D+00, &
          0.3827807429336583D+00, &
          0.6104702436512848D+00, &
          0.3827807429336583D+00, &
          0.6104702436512848D+00, &
          0.0067490134150570D+00, &
          0.0056501992086932D+00, &
          0.2446543644502634D+00, &
          0.7496954363410434D+00, &
          0.2446543644502634D+00, &
          0.7496954363410434D+00, &
          0.0056501992086932D+00, &
          0.0056501992086932D+00, &
          0.2446543644502634D+00, &
          0.7496954363410434D+00, &
          0.2446543644502634D+00, &
          0.7496954363410434D+00, &
          0.0056501992086932D+00, &
          0.6264598409479015D+00, &
          0.0263707136798594D+00, &
          0.3471694453722390D+00, &
          0.0263707136798594D+00, &
          0.3471694453722390D+00, &
          0.6264598409479015D+00, &
          0.6264598409479015D+00, &
          0.0263707136798594D+00, &
          0.3471694453722390D+00, &
          0.0263707136798594D+00, &
          0.3471694453722390D+00, &
          0.6264598409479015D+00, &
          0.1947680793017931D+00, &
          0.0434074542015924D+00, &
          0.7618244664966145D+00, &
          0.0434074542015924D+00, &
          0.7618244664966145D+00, &
          0.1947680793017931D+00, &
          0.1947680793017931D+00, &
          0.0434074542015924D+00, &
          0.7618244664966145D+00, &
          0.0434074542015924D+00, &
          0.7618244664966145D+00, &
          0.1947680793017931D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.8666490488529237D+00, &
          0.1333509511470762D+00, &
          0.8703524297116281D+00, &
          0.8703524297116281D+00, &
          0.8703524297116281D+00, &
          0.1296475702883719D+00, &
          0.1296475702883719D+00, &
          0.1296475702883719D+00, &
          0.5962892344246855D+00, &
          0.5962892344246855D+00, &
          0.5962892344246855D+00, &
          0.4037107655753145D+00, &
          0.4037107655753145D+00, &
          0.4037107655753145D+00, &
          0.9862457308604132D+00, &
          0.9862457308604132D+00, &
          0.9862457308604132D+00, &
          0.0137542691395868D+00, &
          0.0137542691395868D+00, &
          0.0137542691395868D+00, &
          0.6193405024389006D+00, &
          0.6193405024389006D+00, &
          0.6193405024389006D+00, &
          0.3806594975610994D+00, &
          0.3806594975610994D+00, &
          0.3806594975610994D+00, &
          0.6506114926978297D+00, &
          0.6506114926978297D+00, &
          0.6506114926978297D+00, &
          0.3493885073021704D+00, &
          0.3493885073021704D+00, &
          0.3493885073021704D+00, &
          0.8258627018814777D+00, &
          0.8258627018814777D+00, &
          0.8258627018814777D+00, &
          0.1741372981185223D+00, &
          0.1741372981185223D+00, &
          0.1741372981185223D+00, &
          0.2477808238551293D+00, &
          0.2477808238551293D+00, &
          0.2477808238551293D+00, &
          0.7522191761448707D+00, &
          0.7522191761448707D+00, &
          0.7522191761448707D+00, &
          0.2921358009316097D+00, &
          0.2921358009316097D+00, &
          0.2921358009316097D+00, &
          0.7078641990683903D+00, &
          0.7078641990683903D+00, &
          0.7078641990683903D+00, &
          0.0981873258283363D+00, &
          0.0981873258283363D+00, &
          0.0981873258283363D+00, &
          0.9018126741716637D+00, &
          0.9018126741716637D+00, &
          0.9018126741716637D+00, &
          0.7125198639059245D+00, &
          0.7125198639059245D+00, &
          0.7125198639059245D+00, &
          0.2874801360940755D+00, &
          0.2874801360940755D+00, &
          0.2874801360940755D+00, &
          0.3516395946331653D+00, &
          0.3516395946331653D+00, &
          0.3516395946331653D+00, &
          0.6483604053668347D+00, &
          0.6483604053668347D+00, &
          0.6483604053668347D+00, &
          0.9459742145088597D+00, &
          0.9459742145088597D+00, &
          0.9459742145088597D+00, &
          0.0540257854911404D+00, &
          0.0540257854911404D+00, &
          0.0540257854911404D+00, &
          0.5658612299497008D+00, &
          0.5658612299497008D+00, &
          0.5658612299497008D+00, &
          0.4341387700502993D+00, &
          0.4341387700502993D+00, &
          0.4341387700502993D+00, &
          0.9273724754855923D+00, &
          0.9273724754855923D+00, &
          0.9273724754855923D+00, &
          0.0726275245144077D+00, &
          0.0726275245144077D+00, &
          0.0726275245144077D+00, &
          0.9989623565279326D+00, &
          0.9989623565279326D+00, &
          0.9989623565279326D+00, &
          0.0010376434720673D+00, &
          0.0010376434720673D+00, &
          0.0010376434720673D+00, &
          0.1604307070751930D+00, &
          0.1604307070751930D+00, &
          0.1604307070751930D+00, &
          0.8395692929248070D+00, &
          0.8395692929248070D+00, &
          0.8395692929248070D+00, &
          0.6345291354152611D+00, &
          0.6345291354152611D+00, &
          0.6345291354152611D+00, &
          0.3654708645847389D+00, &
          0.3654708645847389D+00, &
          0.3654708645847389D+00, &
          0.9811993443803730D+00, &
          0.9811993443803730D+00, &
          0.9811993443803730D+00, &
          0.0188006556196270D+00, &
          0.0188006556196270D+00, &
          0.0188006556196270D+00, &
          0.9712793580501632D+00, &
          0.9712793580501632D+00, &
          0.9712793580501632D+00, &
          0.0287206419498368D+00, &
          0.0287206419498368D+00, &
          0.0287206419498368D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9553022983267381D+00, &
          0.9553022983267381D+00, &
          0.9553022983267381D+00, &
          0.9553022983267381D+00, &
          0.9553022983267381D+00, &
          0.9553022983267381D+00, &
          0.0446977016732619D+00, &
          0.0446977016732619D+00, &
          0.0446977016732619D+00, &
          0.0446977016732619D+00, &
          0.0446977016732619D+00, &
          0.0446977016732619D+00, &
          0.5287768623754378D+00, &
          0.5287768623754378D+00, &
          0.5287768623754378D+00, &
          0.5287768623754378D+00, &
          0.5287768623754378D+00, &
          0.5287768623754378D+00, &
          0.4712231376245622D+00, &
          0.4712231376245622D+00, &
          0.4712231376245622D+00, &
          0.4712231376245622D+00, &
          0.4712231376245622D+00, &
          0.4712231376245622D+00, &
          0.4538216112941401D+00, &
          0.4538216112941401D+00, &
          0.4538216112941401D+00, &
          0.4538216112941401D+00, &
          0.4538216112941401D+00, &
          0.4538216112941401D+00, &
          0.5461783887058599D+00, &
          0.5461783887058599D+00, &
          0.5461783887058599D+00, &
          0.5461783887058599D+00, &
          0.5461783887058599D+00, &
          0.5461783887058599D+00, &
          0.3624591754389214D+00, &
          0.3624591754389214D+00, &
          0.3624591754389214D+00, &
          0.3624591754389214D+00, &
          0.3624591754389214D+00, &
          0.3624591754389214D+00, &
          0.6375408245610786D+00, &
          0.6375408245610786D+00, &
          0.6375408245610786D+00, &
          0.6375408245610786D+00, &
          0.6375408245610786D+00, &
          0.6375408245610786D+00, &
          0.6176517353683332D+00, &
          0.6176517353683332D+00, &
          0.6176517353683332D+00, &
          0.6176517353683332D+00, &
          0.6176517353683332D+00, &
          0.6176517353683332D+00, &
          0.3823482646316668D+00, &
          0.3823482646316668D+00, &
          0.3823482646316668D+00, &
          0.3823482646316668D+00, &
          0.3823482646316668D+00, &
          0.3823482646316668D+00, &
          0.2130333812420883D+00, &
          0.2130333812420883D+00, &
          0.2130333812420883D+00, &
          0.2130333812420883D+00, &
          0.2130333812420883D+00, &
          0.2130333812420883D+00, &
          0.7869666187579117D+00, &
          0.7869666187579117D+00, &
          0.7869666187579117D+00, &
          0.7869666187579117D+00, &
          0.7869666187579117D+00, &
          0.7869666187579117D+00, &
          0.9281160730015593D+00, &
          0.9281160730015593D+00, &
          0.9281160730015593D+00, &
          0.9281160730015593D+00, &
          0.9281160730015593D+00, &
          0.9281160730015593D+00, &
          0.0718839269984407D+00, &
          0.0718839269984407D+00, &
          0.0718839269984407D+00, &
          0.0718839269984407D+00, &
          0.0718839269984407D+00, &
          0.0718839269984407D+00, &
          0.8760469217983575D+00, &
          0.8760469217983575D+00, &
          0.8760469217983575D+00, &
          0.8760469217983575D+00, &
          0.8760469217983575D+00, &
          0.8760469217983575D+00, &
          0.1239530782016424D+00, &
          0.1239530782016424D+00, &
          0.1239530782016424D+00, &
          0.1239530782016424D+00, &
          0.1239530782016424D+00, &
          0.1239530782016424D+00, &
          0.5105770631364416D+00, &
          0.5105770631364416D+00, &
          0.5105770631364416D+00, &
          0.5105770631364416D+00, &
          0.5105770631364416D+00, &
          0.5105770631364416D+00, &
          0.4894229368635584D+00, &
          0.4894229368635584D+00, &
          0.4894229368635584D+00, &
          0.4894229368635584D+00, &
          0.4894229368635584D+00, &
          0.4894229368635584D+00, &
          0.6124742043718856D+00, &
          0.6124742043718856D+00, &
          0.6124742043718856D+00, &
          0.6124742043718856D+00, &
          0.6124742043718856D+00, &
          0.6124742043718856D+00, &
          0.3875257956281143D+00, &
          0.3875257956281143D+00, &
          0.3875257956281143D+00, &
          0.3875257956281143D+00, &
          0.3875257956281143D+00, &
          0.3875257956281143D+00, &
          0.6912097488491545D+00, &
          0.6912097488491545D+00, &
          0.6912097488491545D+00, &
          0.6912097488491545D+00, &
          0.6912097488491545D+00, &
          0.6912097488491545D+00, &
          0.3087902511508455D+00, &
          0.3087902511508455D+00, &
          0.3087902511508455D+00, &
          0.3087902511508455D+00, &
          0.3087902511508455D+00, &
          0.3087902511508455D+00, &
          0.9870924037489929D+00, &
          0.9870924037489929D+00, &
          0.9870924037489929D+00, &
          0.9870924037489929D+00, &
          0.9870924037489929D+00, &
          0.9870924037489929D+00, &
          0.0129075962510071D+00, &
          0.0129075962510071D+00, &
          0.0129075962510071D+00, &
          0.0129075962510071D+00, &
          0.0129075962510071D+00, &
          0.0129075962510071D+00, &
          0.7190879938703654D+00, &
          0.7190879938703654D+00, &
          0.7190879938703654D+00, &
          0.7190879938703654D+00, &
          0.7190879938703654D+00, &
          0.7190879938703654D+00, &
          0.2809120061296346D+00, &
          0.2809120061296346D+00, &
          0.2809120061296346D+00, &
          0.2809120061296346D+00, &
          0.2809120061296346D+00, &
          0.2809120061296346D+00, &
          0.8842315513280089D+00, &
          0.8842315513280089D+00, &
          0.8842315513280089D+00, &
          0.8842315513280089D+00, &
          0.8842315513280089D+00, &
          0.8842315513280089D+00, &
          0.1157684486719911D+00, &
          0.1157684486719911D+00, &
          0.1157684486719911D+00, &
          0.1157684486719911D+00, &
          0.1157684486719911D+00, &
          0.1157684486719911D+00, &
          0.7732864296384829D+00, &
          0.7732864296384829D+00, &
          0.7732864296384829D+00, &
          0.7732864296384829D+00, &
          0.7732864296384829D+00, &
          0.7732864296384829D+00, &
          0.2267135703615171D+00, &
          0.2267135703615171D+00, &
          0.2267135703615171D+00, &
          0.2267135703615171D+00, &
          0.2267135703615171D+00, &
          0.2267135703615171D+00, &
          0.9276833503490310D+00, &
          0.9276833503490310D+00, &
          0.9276833503490310D+00, &
          0.9276833503490310D+00, &
          0.9276833503490310D+00, &
          0.9276833503490310D+00, &
          0.0723166496509690D+00, &
          0.0723166496509690D+00, &
          0.0723166496509690D+00, &
          0.0723166496509690D+00, &
          0.0723166496509690D+00, &
          0.0723166496509690D+00, &
          0.9877761641032086D+00, &
          0.9877761641032086D+00, &
          0.9877761641032086D+00, &
          0.9877761641032086D+00, &
          0.9877761641032086D+00, &
          0.9877761641032086D+00, &
          0.0122238358967914D+00, &
          0.0122238358967914D+00, &
          0.0122238358967914D+00, &
          0.0122238358967914D+00, &
          0.0122238358967914D+00, &
          0.0122238358967914D+00, &
          0.8232871808501311D+00, &
          0.8232871808501311D+00, &
          0.8232871808501311D+00, &
          0.8232871808501311D+00, &
          0.8232871808501311D+00, &
          0.8232871808501311D+00, &
          0.1767128191498689D+00, &
          0.1767128191498689D+00, &
          0.1767128191498689D+00, &
          0.1767128191498689D+00, &
          0.1767128191498689D+00, &
          0.1767128191498689D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          0.0075216525929737D+00, &
          0.0075216525929737D+00, &
          0.0039523401420815D+00, &
          0.0039523401420815D+00, &
          0.0039523401420815D+00, &
          0.0039523401420815D+00, &
          0.0039523401420815D+00, &
          0.0039523401420815D+00, &
          0.0026587971023890D+00, &
          0.0026587971023890D+00, &
          0.0026587971023890D+00, &
          0.0026587971023890D+00, &
          0.0026587971023890D+00, &
          0.0026587971023890D+00, &
          0.0023032616801324D+00, &
          0.0023032616801324D+00, &
          0.0023032616801324D+00, &
          0.0023032616801324D+00, &
          0.0023032616801324D+00, &
          0.0023032616801324D+00, &
          0.0043006096205643D+00, &
          0.0043006096205643D+00, &
          0.0043006096205643D+00, &
          0.0043006096205643D+00, &
          0.0043006096205643D+00, &
          0.0043006096205643D+00, &
          0.0030221060153018D+00, &
          0.0030221060153018D+00, &
          0.0030221060153018D+00, &
          0.0030221060153018D+00, &
          0.0030221060153018D+00, &
          0.0030221060153018D+00, &
          0.0003076949950084D+00, &
          0.0003076949950084D+00, &
          0.0003076949950084D+00, &
          0.0003076949950084D+00, &
          0.0003076949950084D+00, &
          0.0003076949950084D+00, &
          0.0030734403191510D+00, &
          0.0030734403191510D+00, &
          0.0030734403191510D+00, &
          0.0030734403191510D+00, &
          0.0030734403191510D+00, &
          0.0030734403191510D+00, &
          0.0049055295500793D+00, &
          0.0049055295500793D+00, &
          0.0049055295500793D+00, &
          0.0049055295500793D+00, &
          0.0049055295500793D+00, &
          0.0049055295500793D+00, &
          0.0027066696484092D+00, &
          0.0027066696484092D+00, &
          0.0027066696484092D+00, &
          0.0027066696484092D+00, &
          0.0027066696484092D+00, &
          0.0027066696484092D+00, &
          0.0034577749276135D+00, &
          0.0034577749276135D+00, &
          0.0034577749276135D+00, &
          0.0034577749276135D+00, &
          0.0034577749276135D+00, &
          0.0034577749276135D+00, &
          0.0051705174541971D+00, &
          0.0051705174541971D+00, &
          0.0051705174541971D+00, &
          0.0051705174541971D+00, &
          0.0051705174541971D+00, &
          0.0051705174541971D+00, &
          0.0014022583115306D+00, &
          0.0014022583115306D+00, &
          0.0014022583115306D+00, &
          0.0014022583115306D+00, &
          0.0014022583115306D+00, &
          0.0014022583115306D+00, &
          0.0042401026043811D+00, &
          0.0042401026043811D+00, &
          0.0042401026043811D+00, &
          0.0042401026043811D+00, &
          0.0042401026043811D+00, &
          0.0042401026043811D+00, &
          0.0032841487411702D+00, &
          0.0032841487411702D+00, &
          0.0032841487411702D+00, &
          0.0032841487411702D+00, &
          0.0032841487411702D+00, &
          0.0032841487411702D+00, &
          0.0011342259853124D+00, &
          0.0011342259853124D+00, &
          0.0011342259853124D+00, &
          0.0011342259853124D+00, &
          0.0011342259853124D+00, &
          0.0011342259853124D+00, &
          0.0074708548163608D+00, &
          0.0074708548163608D+00, &
          0.0074708548163608D+00, &
          0.0074708548163608D+00, &
          0.0074708548163608D+00, &
          0.0074708548163608D+00, &
          0.0034340384613521D+00, &
          0.0034340384613521D+00, &
          0.0034340384613521D+00, &
          0.0034340384613521D+00, &
          0.0034340384613521D+00, &
          0.0034340384613521D+00, &
          0.0031462609694126D+00, &
          0.0031462609694126D+00, &
          0.0031462609694126D+00, &
          0.0031462609694126D+00, &
          0.0031462609694126D+00, &
          0.0031462609694126D+00, &
          0.0005752613942102D+00, &
          0.0005752613942102D+00, &
          0.0005752613942102D+00, &
          0.0005752613942102D+00, &
          0.0005752613942102D+00, &
          0.0005752613942102D+00, &
          0.0065552795530324D+00, &
          0.0065552795530324D+00, &
          0.0065552795530324D+00, &
          0.0065552795530324D+00, &
          0.0065552795530324D+00, &
          0.0065552795530324D+00, &
          0.0029454095619872D+00, &
          0.0029454095619872D+00, &
          0.0029454095619872D+00, &
          0.0029454095619872D+00, &
          0.0029454095619872D+00, &
          0.0029454095619872D+00, &
          0.0029454095619872D+00, &
          0.0029454095619872D+00, &
          0.0029454095619872D+00, &
          0.0029454095619872D+00, &
          0.0029454095619872D+00, &
          0.0029454095619872D+00, &
          0.0024530630415632D+00, &
          0.0024530630415632D+00, &
          0.0024530630415632D+00, &
          0.0024530630415632D+00, &
          0.0024530630415632D+00, &
          0.0024530630415632D+00, &
          0.0024530630415632D+00, &
          0.0024530630415632D+00, &
          0.0024530630415632D+00, &
          0.0024530630415632D+00, &
          0.0024530630415632D+00, &
          0.0024530630415632D+00, &
          0.0028703966404990D+00, &
          0.0028703966404990D+00, &
          0.0028703966404990D+00, &
          0.0028703966404990D+00, &
          0.0028703966404990D+00, &
          0.0028703966404990D+00, &
          0.0028703966404990D+00, &
          0.0028703966404990D+00, &
          0.0028703966404990D+00, &
          0.0028703966404990D+00, &
          0.0028703966404990D+00, &
          0.0028703966404990D+00, &
          0.0018268277900660D+00, &
          0.0018268277900660D+00, &
          0.0018268277900660D+00, &
          0.0018268277900660D+00, &
          0.0018268277900660D+00, &
          0.0018268277900660D+00, &
          0.0018268277900660D+00, &
          0.0018268277900660D+00, &
          0.0018268277900660D+00, &
          0.0018268277900660D+00, &
          0.0018268277900660D+00, &
          0.0018268277900660D+00, &
          0.0046842877207978D+00, &
          0.0046842877207978D+00, &
          0.0046842877207978D+00, &
          0.0046842877207978D+00, &
          0.0046842877207978D+00, &
          0.0046842877207978D+00, &
          0.0046842877207978D+00, &
          0.0046842877207978D+00, &
          0.0046842877207978D+00, &
          0.0046842877207978D+00, &
          0.0046842877207978D+00, &
          0.0046842877207978D+00, &
          0.0061987832367774D+00, &
          0.0061987832367774D+00, &
          0.0061987832367774D+00, &
          0.0061987832367774D+00, &
          0.0061987832367774D+00, &
          0.0061987832367774D+00, &
          0.0061987832367774D+00, &
          0.0061987832367774D+00, &
          0.0061987832367774D+00, &
          0.0061987832367774D+00, &
          0.0061987832367774D+00, &
          0.0061987832367774D+00, &
          0.0041978630490507D+00, &
          0.0041978630490507D+00, &
          0.0041978630490507D+00, &
          0.0041978630490507D+00, &
          0.0041978630490507D+00, &
          0.0041978630490507D+00, &
          0.0041978630490507D+00, &
          0.0041978630490507D+00, &
          0.0041978630490507D+00, &
          0.0041978630490507D+00, &
          0.0041978630490507D+00, &
          0.0041978630490507D+00, &
          0.0025528142014630D+00, &
          0.0025528142014630D+00, &
          0.0025528142014630D+00, &
          0.0025528142014630D+00, &
          0.0025528142014630D+00, &
          0.0025528142014630D+00, &
          0.0025528142014630D+00, &
          0.0025528142014630D+00, &
          0.0025528142014630D+00, &
          0.0025528142014630D+00, &
          0.0025528142014630D+00, &
          0.0025528142014630D+00, &
          0.0020901151709172D+00, &
          0.0020901151709172D+00, &
          0.0020901151709172D+00, &
          0.0020901151709172D+00, &
          0.0020901151709172D+00, &
          0.0020901151709172D+00, &
          0.0020901151709172D+00, &
          0.0020901151709172D+00, &
          0.0020901151709172D+00, &
          0.0020901151709172D+00, &
          0.0020901151709172D+00, &
          0.0020901151709172D+00, &
          0.0007568671403396D+00, &
          0.0007568671403396D+00, &
          0.0007568671403396D+00, &
          0.0007568671403396D+00, &
          0.0007568671403396D+00, &
          0.0007568671403396D+00, &
          0.0007568671403396D+00, &
          0.0007568671403396D+00, &
          0.0007568671403396D+00, &
          0.0007568671403396D+00, &
          0.0007568671403396D+00, &
          0.0007568671403396D+00, &
          0.0013000657288395D+00, &
          0.0013000657288395D+00, &
          0.0013000657288395D+00, &
          0.0013000657288395D+00, &
          0.0013000657288395D+00, &
          0.0013000657288395D+00, &
          0.0013000657288395D+00, &
          0.0013000657288395D+00, &
          0.0013000657288395D+00, &
          0.0013000657288395D+00, &
          0.0013000657288395D+00, &
          0.0013000657288395D+00, &
          0.0008479873802968D+00, &
          0.0008479873802968D+00, &
          0.0008479873802968D+00, &
          0.0008479873802968D+00, &
          0.0008479873802968D+00, &
          0.0008479873802968D+00, &
          0.0008479873802968D+00, &
          0.0008479873802968D+00, &
          0.0008479873802968D+00, &
          0.0008479873802968D+00, &
          0.0008479873802968D+00, &
          0.0008479873802968D+00, &
          0.0064046130401019D+00, &
          0.0064046130401019D+00, &
          0.0064046130401019D+00, &
          0.0064046130401019D+00, &
          0.0064046130401019D+00, &
          0.0064046130401019D+00, &
          0.0064046130401019D+00, &
          0.0064046130401019D+00, &
          0.0064046130401019D+00, &
          0.0064046130401019D+00, &
          0.0064046130401019D+00, &
          0.0064046130401019D+00, &
          0.0011323520941966D+00, &
          0.0011323520941966D+00, &
          0.0011323520941966D+00, &
          0.0011323520941966D+00, &
          0.0011323520941966D+00, &
          0.0011323520941966D+00, &
          0.0011323520941966D+00, &
          0.0011323520941966D+00, &
          0.0011323520941966D+00, &
          0.0011323520941966D+00, &
          0.0011323520941966D+00, &
          0.0011323520941966D+00, &
          0.0020644102531444D+00, &
          0.0020644102531444D+00, &
          0.0020644102531444D+00, &
          0.0020644102531444D+00, &
          0.0020644102531444D+00, &
          0.0020644102531444D+00, &
          0.0020644102531444D+00, &
          0.0020644102531444D+00, &
          0.0020644102531444D+00, &
          0.0020644102531444D+00, &
          0.0020644102531444D+00, &
          0.0020644102531444D+00, &
          0.0010481934583324D+00, &
          0.0010481934583324D+00, &
          0.0010481934583324D+00, &
          0.0010481934583324D+00, &
          0.0010481934583324D+00, &
          0.0010481934583324D+00, &
          0.0010481934583324D+00, &
          0.0010481934583324D+00, &
          0.0010481934583324D+00, &
          0.0010481934583324D+00, &
          0.0010481934583324D+00, &
          0.0010481934583324D+00, &
          0.0009390285812953D+00, &
          0.0009390285812953D+00, &
          0.0009390285812953D+00, &
          0.0009390285812953D+00, &
          0.0009390285812953D+00, &
          0.0009390285812953D+00, &
          0.0009390285812953D+00, &
          0.0009390285812953D+00, &
          0.0009390285812953D+00, &
          0.0009390285812953D+00, &
          0.0009390285812953D+00, &
          0.0009390285812953D+00, &
          0.0042160603323249D+00, &
          0.0042160603323249D+00, &
          0.0042160603323249D+00, &
          0.0042160603323249D+00, &
          0.0042160603323249D+00, &
          0.0042160603323249D+00, &
          0.0042160603323249D+00, &
          0.0042160603323249D+00, &
          0.0042160603323249D+00, &
          0.0042160603323249D+00, &
          0.0042160603323249D+00, &
          0.0042160603323249D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule18 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule18() returns the prism rule of precision 18.
!
!  Discussion:
!
!    We suppose we are given a triangular prism P with vertices
!    (1,0,0), (0,1,0), (0,0,0),
!    (1,0,1), (0,1,1), (0,0,1),
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over P is 
!    approximated by Q as follows:
!
!    Q = volume(P) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 May 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order symmetric cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 122, Number 1, pages 148-171, 24 August 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real x(n), y(n), z(n): the coordinates of quadrature points.
!
!    real w(n): the quadrature weights, which sum to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: n_save = 396

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
          0.4728631477307805D+00, &
          0.4728631477307805D+00, &
          0.0542737045384390D+00, &
          0.4728631477307805D+00, &
          0.4728631477307805D+00, &
          0.0542737045384390D+00, &
          0.1090318158369727D+00, &
          0.1090318158369727D+00, &
          0.7819363683260546D+00, &
          0.1090318158369727D+00, &
          0.1090318158369727D+00, &
          0.7819363683260546D+00, &
          0.2141954002982068D+00, &
          0.2141954002982068D+00, &
          0.5716091994035865D+00, &
          0.2141954002982068D+00, &
          0.2141954002982068D+00, &
          0.5716091994035865D+00, &
          0.0084150887533426D+00, &
          0.0084150887533426D+00, &
          0.9831698224933149D+00, &
          0.0084150887533426D+00, &
          0.0084150887533426D+00, &
          0.9831698224933149D+00, &
          0.0494156440294096D+00, &
          0.0494156440294096D+00, &
          0.9011687119411808D+00, &
          0.0494156440294096D+00, &
          0.0494156440294096D+00, &
          0.9011687119411808D+00, &
          0.4574258170355346D+00, &
          0.4574258170355346D+00, &
          0.0851483659289308D+00, &
          0.4574258170355346D+00, &
          0.4574258170355346D+00, &
          0.0851483659289308D+00, &
          0.2945570490820424D+00, &
          0.2945570490820424D+00, &
          0.4108859018359152D+00, &
          0.2945570490820424D+00, &
          0.2945570490820424D+00, &
          0.4108859018359152D+00, &
          0.4148583469418042D+00, &
          0.4148583469418042D+00, &
          0.1702833061163916D+00, &
          0.4148583469418042D+00, &
          0.4148583469418042D+00, &
          0.1702833061163916D+00, &
          0.1801567120271234D+00, &
          0.1801567120271234D+00, &
          0.6396865759457532D+00, &
          0.1801567120271234D+00, &
          0.1801567120271234D+00, &
          0.6396865759457532D+00, &
          0.4092592655375933D+00, &
          0.4092592655375933D+00, &
          0.1814814689248134D+00, &
          0.4092592655375933D+00, &
          0.4092592655375933D+00, &
          0.1814814689248134D+00, &
          0.4106876572727113D+00, &
          0.4106876572727113D+00, &
          0.1786246854545774D+00, &
          0.4106876572727113D+00, &
          0.4106876572727113D+00, &
          0.1786246854545774D+00, &
          0.4942616172602019D+00, &
          0.4942616172602019D+00, &
          0.0114767654795962D+00, &
          0.4942616172602019D+00, &
          0.4942616172602019D+00, &
          0.0114767654795962D+00, &
          0.1455821883151533D+00, &
          0.1455821883151533D+00, &
          0.7088356233696933D+00, &
          0.1455821883151533D+00, &
          0.1455821883151533D+00, &
          0.7088356233696933D+00, &
          0.0775928696462990D+00, &
          0.0775928696462990D+00, &
          0.8448142607074021D+00, &
          0.0775928696462990D+00, &
          0.0775928696462990D+00, &
          0.8448142607074021D+00, &
          0.1688466480848478D+00, &
          0.1688466480848478D+00, &
          0.6623067038303044D+00, &
          0.1688466480848478D+00, &
          0.1688466480848478D+00, &
          0.6623067038303044D+00, &
          0.2737821968608845D+00, &
          0.2737821968608845D+00, &
          0.4524356062782311D+00, &
          0.2737821968608845D+00, &
          0.2737821968608845D+00, &
          0.4524356062782311D+00, &
          0.3193703110432156D+00, &
          0.3193703110432156D+00, &
          0.3612593779135688D+00, &
          0.3193703110432156D+00, &
          0.3193703110432156D+00, &
          0.3612593779135688D+00, &
          0.2789278090499041D+00, &
          0.2789278090499041D+00, &
          0.4421443819001918D+00, &
          0.2789278090499041D+00, &
          0.2789278090499041D+00, &
          0.4421443819001918D+00, &
          0.0147343308678136D+00, &
          0.0147343308678136D+00, &
          0.9705313382643729D+00, &
          0.0147343308678136D+00, &
          0.0147343308678136D+00, &
          0.9705313382643729D+00, &
          0.2803931989784251D+00, &
          0.1730708483947652D+00, &
          0.2803931989784251D+00, &
          0.5465359526268097D+00, &
          0.1730708483947652D+00, &
          0.5465359526268097D+00, &
          0.3067402740370590D+00, &
          0.0472583316199373D+00, &
          0.3067402740370590D+00, &
          0.6460013943430037D+00, &
          0.0472583316199373D+00, &
          0.6460013943430037D+00, &
          0.3067402740370590D+00, &
          0.0472583316199373D+00, &
          0.3067402740370590D+00, &
          0.6460013943430037D+00, &
          0.0472583316199373D+00, &
          0.6460013943430037D+00, &
          0.5587472044782793D+00, &
          0.0637486493547561D+00, &
          0.5587472044782793D+00, &
          0.3775041461669645D+00, &
          0.0637486493547561D+00, &
          0.3775041461669645D+00, &
          0.5587472044782793D+00, &
          0.0637486493547561D+00, &
          0.5587472044782793D+00, &
          0.3775041461669645D+00, &
          0.0637486493547561D+00, &
          0.3775041461669645D+00, &
          0.1295467219204834D+00, &
          0.2834002542935580D+00, &
          0.1295467219204834D+00, &
          0.5870530237859586D+00, &
          0.2834002542935580D+00, &
          0.5870530237859586D+00, &
          0.1295467219204834D+00, &
          0.2834002542935580D+00, &
          0.1295467219204834D+00, &
          0.5870530237859586D+00, &
          0.2834002542935580D+00, &
          0.5870530237859586D+00, &
          0.2943165499688927D+00, &
          0.0323606160924521D+00, &
          0.2943165499688927D+00, &
          0.6733228339386552D+00, &
          0.0323606160924521D+00, &
          0.6733228339386552D+00, &
          0.2943165499688927D+00, &
          0.0323606160924521D+00, &
          0.2943165499688927D+00, &
          0.6733228339386552D+00, &
          0.0323606160924521D+00, &
          0.6733228339386552D+00, &
          0.7045618883179275D+00, &
          0.0682727115900655D+00, &
          0.7045618883179275D+00, &
          0.2271654000920070D+00, &
          0.0682727115900655D+00, &
          0.2271654000920070D+00, &
          0.7045618883179275D+00, &
          0.0682727115900655D+00, &
          0.7045618883179275D+00, &
          0.2271654000920070D+00, &
          0.0682727115900655D+00, &
          0.2271654000920070D+00, &
          0.5336092397180819D+00, &
          0.0823227847933820D+00, &
          0.5336092397180819D+00, &
          0.3840679754885361D+00, &
          0.0823227847933820D+00, &
          0.3840679754885361D+00, &
          0.5336092397180819D+00, &
          0.0823227847933820D+00, &
          0.5336092397180819D+00, &
          0.3840679754885361D+00, &
          0.0823227847933820D+00, &
          0.3840679754885361D+00, &
          0.1659343236267790D+00, &
          0.2958566157742775D+00, &
          0.1659343236267790D+00, &
          0.5382090605989435D+00, &
          0.2958566157742775D+00, &
          0.5382090605989435D+00, &
          0.1659343236267790D+00, &
          0.2958566157742775D+00, &
          0.1659343236267790D+00, &
          0.5382090605989435D+00, &
          0.2958566157742775D+00, &
          0.5382090605989435D+00, &
          0.5873912217089957D+00, &
          0.0091538552254334D+00, &
          0.5873912217089957D+00, &
          0.4034549230655710D+00, &
          0.0091538552254334D+00, &
          0.4034549230655710D+00, &
          0.5873912217089957D+00, &
          0.0091538552254334D+00, &
          0.5873912217089957D+00, &
          0.4034549230655710D+00, &
          0.0091538552254334D+00, &
          0.4034549230655710D+00, &
          0.1216404492275597D+00, &
          0.0492455683577680D+00, &
          0.1216404492275597D+00, &
          0.8291139824146723D+00, &
          0.0492455683577680D+00, &
          0.8291139824146723D+00, &
          0.1216404492275597D+00, &
          0.0492455683577680D+00, &
          0.1216404492275597D+00, &
          0.8291139824146723D+00, &
          0.0492455683577680D+00, &
          0.8291139824146723D+00, &
          0.0420809533226036D+00, &
          0.0077537345196640D+00, &
          0.0420809533226036D+00, &
          0.9501653121577324D+00, &
          0.0077537345196640D+00, &
          0.9501653121577324D+00, &
          0.0420809533226036D+00, &
          0.0077537345196640D+00, &
          0.0420809533226036D+00, &
          0.9501653121577324D+00, &
          0.0077537345196640D+00, &
          0.9501653121577324D+00, &
          0.4516351625773616D+00, &
          0.0228219169717125D+00, &
          0.4516351625773616D+00, &
          0.5255429204509259D+00, &
          0.0228219169717125D+00, &
          0.5255429204509259D+00, &
          0.4516351625773616D+00, &
          0.0228219169717125D+00, &
          0.4516351625773616D+00, &
          0.5255429204509259D+00, &
          0.0228219169717125D+00, &
          0.5255429204509259D+00, &
          0.1252580641579170D+00, &
          0.0073963274047276D+00, &
          0.1252580641579170D+00, &
          0.8673456084373554D+00, &
          0.0073963274047276D+00, &
          0.8673456084373554D+00, &
          0.1252580641579170D+00, &
          0.0073963274047276D+00, &
          0.1252580641579170D+00, &
          0.8673456084373554D+00, &
          0.0073963274047276D+00, &
          0.8673456084373554D+00, &
          0.1765722195894258D+00, &
          0.0653002115725378D+00, &
          0.1765722195894258D+00, &
          0.7581275688380364D+00, &
          0.0653002115725378D+00, &
          0.7581275688380364D+00, &
          0.1765722195894258D+00, &
          0.0653002115725378D+00, &
          0.1765722195894258D+00, &
          0.7581275688380364D+00, &
          0.0653002115725378D+00, &
          0.7581275688380364D+00, &
          0.2404267015692479D+00, &
          0.0085879015506819D+00, &
          0.2404267015692479D+00, &
          0.7509853968800702D+00, &
          0.0085879015506819D+00, &
          0.7509853968800702D+00, &
          0.2404267015692479D+00, &
          0.0085879015506819D+00, &
          0.2404267015692479D+00, &
          0.7509853968800702D+00, &
          0.0085879015506819D+00, &
          0.7509853968800702D+00, &
          0.0644058552263157D+00, &
          0.0122678034366560D+00, &
          0.0644058552263157D+00, &
          0.9233263413370284D+00, &
          0.0122678034366560D+00, &
          0.9233263413370284D+00, &
          0.0644058552263157D+00, &
          0.0122678034366560D+00, &
          0.0644058552263157D+00, &
          0.9233263413370284D+00, &
          0.0122678034366560D+00, &
          0.9233263413370284D+00, &
          0.2681069831840066D+00, &
          0.0044793265613672D+00, &
          0.2681069831840066D+00, &
          0.7274136902546262D+00, &
          0.0044793265613672D+00, &
          0.7274136902546262D+00, &
          0.2681069831840066D+00, &
          0.0044793265613672D+00, &
          0.2681069831840066D+00, &
          0.7274136902546262D+00, &
          0.0044793265613672D+00, &
          0.7274136902546262D+00, &
          0.1727205252244992D+00, &
          0.0117334535244670D+00, &
          0.1727205252244992D+00, &
          0.8155460212510338D+00, &
          0.0117334535244670D+00, &
          0.8155460212510338D+00, &
          0.1727205252244992D+00, &
          0.0117334535244670D+00, &
          0.1727205252244992D+00, &
          0.8155460212510338D+00, &
          0.0117334535244670D+00, &
          0.8155460212510338D+00, &
          0.1118329953359960D+00, &
          0.5446775610503505D+00, &
          0.1118329953359960D+00, &
          0.3434894436136536D+00, &
          0.5446775610503505D+00, &
          0.3434894436136536D+00, &
          0.1118329953359960D+00, &
          0.5446775610503505D+00, &
          0.1118329953359960D+00, &
          0.3434894436136536D+00, &
          0.5446775610503505D+00, &
          0.3434894436136536D+00, &
          0.5899536365475773D+00, &
          0.0080175874177025D+00, &
          0.5899536365475773D+00, &
          0.4020287760347202D+00, &
          0.0080175874177025D+00, &
          0.4020287760347202D+00, &
          0.5899536365475773D+00, &
          0.0080175874177025D+00, &
          0.5899536365475773D+00, &
          0.4020287760347202D+00, &
          0.0080175874177025D+00, &
          0.4020287760347202D+00, &
          0.0961540597861764D+00, &
          0.2067498834662211D+00, &
          0.0961540597861764D+00, &
          0.6970960567476024D+00, &
          0.2067498834662211D+00, &
          0.6970960567476024D+00, &
          0.0961540597861764D+00, &
          0.2067498834662211D+00, &
          0.0961540597861764D+00, &
          0.6970960567476024D+00, &
          0.2067498834662211D+00, &
          0.6970960567476024D+00, &
          0.3068968355307726D+00, &
          0.0117746553239110D+00, &
          0.3068968355307726D+00, &
          0.6813285091453164D+00, &
          0.0117746553239110D+00, &
          0.6813285091453164D+00, &
          0.3068968355307726D+00, &
          0.0117746553239110D+00, &
          0.3068968355307726D+00, &
          0.6813285091453164D+00, &
          0.0117746553239110D+00, &
          0.6813285091453164D+00, &
          0.0798486447852760D+00, &
          0.0174063356195718D+00, &
          0.0798486447852760D+00, &
          0.9027450195951523D+00, &
          0.0174063356195718D+00, &
          0.9027450195951523D+00, &
          0.0798486447852760D+00, &
          0.0174063356195718D+00, &
          0.0798486447852760D+00, &
          0.9027450195951523D+00, &
          0.0174063356195718D+00, &
          0.9027450195951523D+00, &
          0.0360777790627768D+00, &
          0.1618562004775664D+00, &
          0.0360777790627768D+00, &
          0.8020660204596567D+00, &
          0.1618562004775664D+00, &
          0.8020660204596567D+00, &
          0.0360777790627768D+00, &
          0.1618562004775664D+00, &
          0.0360777790627768D+00, &
          0.8020660204596567D+00, &
          0.1618562004775664D+00, &
          0.8020660204596567D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.4728631477307805D+00, &
          0.0542737045384390D+00, &
          0.4728631477307805D+00, &
          0.4728631477307805D+00, &
          0.0542737045384390D+00, &
          0.4728631477307805D+00, &
          0.1090318158369727D+00, &
          0.7819363683260546D+00, &
          0.1090318158369727D+00, &
          0.1090318158369727D+00, &
          0.7819363683260546D+00, &
          0.1090318158369727D+00, &
          0.2141954002982068D+00, &
          0.5716091994035865D+00, &
          0.2141954002982068D+00, &
          0.2141954002982068D+00, &
          0.5716091994035865D+00, &
          0.2141954002982068D+00, &
          0.0084150887533426D+00, &
          0.9831698224933149D+00, &
          0.0084150887533426D+00, &
          0.0084150887533426D+00, &
          0.9831698224933149D+00, &
          0.0084150887533426D+00, &
          0.0494156440294096D+00, &
          0.9011687119411808D+00, &
          0.0494156440294096D+00, &
          0.0494156440294096D+00, &
          0.9011687119411808D+00, &
          0.0494156440294096D+00, &
          0.4574258170355346D+00, &
          0.0851483659289308D+00, &
          0.4574258170355346D+00, &
          0.4574258170355346D+00, &
          0.0851483659289308D+00, &
          0.4574258170355346D+00, &
          0.2945570490820424D+00, &
          0.4108859018359152D+00, &
          0.2945570490820424D+00, &
          0.2945570490820424D+00, &
          0.4108859018359152D+00, &
          0.2945570490820424D+00, &
          0.4148583469418042D+00, &
          0.1702833061163916D+00, &
          0.4148583469418042D+00, &
          0.4148583469418042D+00, &
          0.1702833061163916D+00, &
          0.4148583469418042D+00, &
          0.1801567120271234D+00, &
          0.6396865759457532D+00, &
          0.1801567120271234D+00, &
          0.1801567120271234D+00, &
          0.6396865759457532D+00, &
          0.1801567120271234D+00, &
          0.4092592655375933D+00, &
          0.1814814689248134D+00, &
          0.4092592655375933D+00, &
          0.4092592655375933D+00, &
          0.1814814689248134D+00, &
          0.4092592655375933D+00, &
          0.4106876572727113D+00, &
          0.1786246854545774D+00, &
          0.4106876572727113D+00, &
          0.4106876572727113D+00, &
          0.1786246854545774D+00, &
          0.4106876572727113D+00, &
          0.4942616172602019D+00, &
          0.0114767654795962D+00, &
          0.4942616172602019D+00, &
          0.4942616172602019D+00, &
          0.0114767654795962D+00, &
          0.4942616172602019D+00, &
          0.1455821883151533D+00, &
          0.7088356233696933D+00, &
          0.1455821883151533D+00, &
          0.1455821883151533D+00, &
          0.7088356233696933D+00, &
          0.1455821883151533D+00, &
          0.0775928696462990D+00, &
          0.8448142607074021D+00, &
          0.0775928696462990D+00, &
          0.0775928696462990D+00, &
          0.8448142607074021D+00, &
          0.0775928696462990D+00, &
          0.1688466480848478D+00, &
          0.6623067038303044D+00, &
          0.1688466480848478D+00, &
          0.1688466480848478D+00, &
          0.6623067038303044D+00, &
          0.1688466480848478D+00, &
          0.2737821968608845D+00, &
          0.4524356062782311D+00, &
          0.2737821968608845D+00, &
          0.2737821968608845D+00, &
          0.4524356062782311D+00, &
          0.2737821968608845D+00, &
          0.3193703110432156D+00, &
          0.3612593779135688D+00, &
          0.3193703110432156D+00, &
          0.3193703110432156D+00, &
          0.3612593779135688D+00, &
          0.3193703110432156D+00, &
          0.2789278090499041D+00, &
          0.4421443819001918D+00, &
          0.2789278090499041D+00, &
          0.2789278090499041D+00, &
          0.4421443819001918D+00, &
          0.2789278090499041D+00, &
          0.0147343308678136D+00, &
          0.9705313382643729D+00, &
          0.0147343308678136D+00, &
          0.0147343308678136D+00, &
          0.9705313382643729D+00, &
          0.0147343308678136D+00, &
          0.1730708483947652D+00, &
          0.2803931989784251D+00, &
          0.5465359526268097D+00, &
          0.2803931989784251D+00, &
          0.5465359526268097D+00, &
          0.1730708483947652D+00, &
          0.0472583316199373D+00, &
          0.3067402740370590D+00, &
          0.6460013943430037D+00, &
          0.3067402740370590D+00, &
          0.6460013943430037D+00, &
          0.0472583316199373D+00, &
          0.0472583316199373D+00, &
          0.3067402740370590D+00, &
          0.6460013943430037D+00, &
          0.3067402740370590D+00, &
          0.6460013943430037D+00, &
          0.0472583316199373D+00, &
          0.0637486493547561D+00, &
          0.5587472044782793D+00, &
          0.3775041461669645D+00, &
          0.5587472044782793D+00, &
          0.3775041461669645D+00, &
          0.0637486493547561D+00, &
          0.0637486493547561D+00, &
          0.5587472044782793D+00, &
          0.3775041461669645D+00, &
          0.5587472044782793D+00, &
          0.3775041461669645D+00, &
          0.0637486493547561D+00, &
          0.2834002542935580D+00, &
          0.1295467219204834D+00, &
          0.5870530237859586D+00, &
          0.1295467219204834D+00, &
          0.5870530237859586D+00, &
          0.2834002542935580D+00, &
          0.2834002542935580D+00, &
          0.1295467219204834D+00, &
          0.5870530237859586D+00, &
          0.1295467219204834D+00, &
          0.5870530237859586D+00, &
          0.2834002542935580D+00, &
          0.0323606160924521D+00, &
          0.2943165499688927D+00, &
          0.6733228339386552D+00, &
          0.2943165499688927D+00, &
          0.6733228339386552D+00, &
          0.0323606160924521D+00, &
          0.0323606160924521D+00, &
          0.2943165499688927D+00, &
          0.6733228339386552D+00, &
          0.2943165499688927D+00, &
          0.6733228339386552D+00, &
          0.0323606160924521D+00, &
          0.0682727115900655D+00, &
          0.7045618883179275D+00, &
          0.2271654000920070D+00, &
          0.7045618883179275D+00, &
          0.2271654000920070D+00, &
          0.0682727115900655D+00, &
          0.0682727115900655D+00, &
          0.7045618883179275D+00, &
          0.2271654000920070D+00, &
          0.7045618883179275D+00, &
          0.2271654000920070D+00, &
          0.0682727115900655D+00, &
          0.0823227847933820D+00, &
          0.5336092397180819D+00, &
          0.3840679754885361D+00, &
          0.5336092397180819D+00, &
          0.3840679754885361D+00, &
          0.0823227847933820D+00, &
          0.0823227847933820D+00, &
          0.5336092397180819D+00, &
          0.3840679754885361D+00, &
          0.5336092397180819D+00, &
          0.3840679754885361D+00, &
          0.0823227847933820D+00, &
          0.2958566157742775D+00, &
          0.1659343236267790D+00, &
          0.5382090605989435D+00, &
          0.1659343236267790D+00, &
          0.5382090605989435D+00, &
          0.2958566157742775D+00, &
          0.2958566157742775D+00, &
          0.1659343236267790D+00, &
          0.5382090605989435D+00, &
          0.1659343236267790D+00, &
          0.5382090605989435D+00, &
          0.2958566157742775D+00, &
          0.0091538552254334D+00, &
          0.5873912217089957D+00, &
          0.4034549230655710D+00, &
          0.5873912217089957D+00, &
          0.4034549230655710D+00, &
          0.0091538552254334D+00, &
          0.0091538552254334D+00, &
          0.5873912217089957D+00, &
          0.4034549230655710D+00, &
          0.5873912217089957D+00, &
          0.4034549230655710D+00, &
          0.0091538552254334D+00, &
          0.0492455683577680D+00, &
          0.1216404492275597D+00, &
          0.8291139824146723D+00, &
          0.1216404492275597D+00, &
          0.8291139824146723D+00, &
          0.0492455683577680D+00, &
          0.0492455683577680D+00, &
          0.1216404492275597D+00, &
          0.8291139824146723D+00, &
          0.1216404492275597D+00, &
          0.8291139824146723D+00, &
          0.0492455683577680D+00, &
          0.0077537345196640D+00, &
          0.0420809533226036D+00, &
          0.9501653121577324D+00, &
          0.0420809533226036D+00, &
          0.9501653121577324D+00, &
          0.0077537345196640D+00, &
          0.0077537345196640D+00, &
          0.0420809533226036D+00, &
          0.9501653121577324D+00, &
          0.0420809533226036D+00, &
          0.9501653121577324D+00, &
          0.0077537345196640D+00, &
          0.0228219169717125D+00, &
          0.4516351625773616D+00, &
          0.5255429204509259D+00, &
          0.4516351625773616D+00, &
          0.5255429204509259D+00, &
          0.0228219169717125D+00, &
          0.0228219169717125D+00, &
          0.4516351625773616D+00, &
          0.5255429204509259D+00, &
          0.4516351625773616D+00, &
          0.5255429204509259D+00, &
          0.0228219169717125D+00, &
          0.0073963274047276D+00, &
          0.1252580641579170D+00, &
          0.8673456084373554D+00, &
          0.1252580641579170D+00, &
          0.8673456084373554D+00, &
          0.0073963274047276D+00, &
          0.0073963274047276D+00, &
          0.1252580641579170D+00, &
          0.8673456084373554D+00, &
          0.1252580641579170D+00, &
          0.8673456084373554D+00, &
          0.0073963274047276D+00, &
          0.0653002115725378D+00, &
          0.1765722195894258D+00, &
          0.7581275688380364D+00, &
          0.1765722195894258D+00, &
          0.7581275688380364D+00, &
          0.0653002115725378D+00, &
          0.0653002115725378D+00, &
          0.1765722195894258D+00, &
          0.7581275688380364D+00, &
          0.1765722195894258D+00, &
          0.7581275688380364D+00, &
          0.0653002115725378D+00, &
          0.0085879015506819D+00, &
          0.2404267015692479D+00, &
          0.7509853968800702D+00, &
          0.2404267015692479D+00, &
          0.7509853968800702D+00, &
          0.0085879015506819D+00, &
          0.0085879015506819D+00, &
          0.2404267015692479D+00, &
          0.7509853968800702D+00, &
          0.2404267015692479D+00, &
          0.7509853968800702D+00, &
          0.0085879015506819D+00, &
          0.0122678034366560D+00, &
          0.0644058552263157D+00, &
          0.9233263413370284D+00, &
          0.0644058552263157D+00, &
          0.9233263413370284D+00, &
          0.0122678034366560D+00, &
          0.0122678034366560D+00, &
          0.0644058552263157D+00, &
          0.9233263413370284D+00, &
          0.0644058552263157D+00, &
          0.9233263413370284D+00, &
          0.0122678034366560D+00, &
          0.0044793265613672D+00, &
          0.2681069831840066D+00, &
          0.7274136902546262D+00, &
          0.2681069831840066D+00, &
          0.7274136902546262D+00, &
          0.0044793265613672D+00, &
          0.0044793265613672D+00, &
          0.2681069831840066D+00, &
          0.7274136902546262D+00, &
          0.2681069831840066D+00, &
          0.7274136902546262D+00, &
          0.0044793265613672D+00, &
          0.0117334535244670D+00, &
          0.1727205252244992D+00, &
          0.8155460212510338D+00, &
          0.1727205252244992D+00, &
          0.8155460212510338D+00, &
          0.0117334535244670D+00, &
          0.0117334535244670D+00, &
          0.1727205252244992D+00, &
          0.8155460212510338D+00, &
          0.1727205252244992D+00, &
          0.8155460212510338D+00, &
          0.0117334535244670D+00, &
          0.5446775610503505D+00, &
          0.1118329953359960D+00, &
          0.3434894436136536D+00, &
          0.1118329953359960D+00, &
          0.3434894436136536D+00, &
          0.5446775610503505D+00, &
          0.5446775610503505D+00, &
          0.1118329953359960D+00, &
          0.3434894436136536D+00, &
          0.1118329953359960D+00, &
          0.3434894436136536D+00, &
          0.5446775610503505D+00, &
          0.0080175874177025D+00, &
          0.5899536365475773D+00, &
          0.4020287760347202D+00, &
          0.5899536365475773D+00, &
          0.4020287760347202D+00, &
          0.0080175874177025D+00, &
          0.0080175874177025D+00, &
          0.5899536365475773D+00, &
          0.4020287760347202D+00, &
          0.5899536365475773D+00, &
          0.4020287760347202D+00, &
          0.0080175874177025D+00, &
          0.2067498834662211D+00, &
          0.0961540597861764D+00, &
          0.6970960567476024D+00, &
          0.0961540597861764D+00, &
          0.6970960567476024D+00, &
          0.2067498834662211D+00, &
          0.2067498834662211D+00, &
          0.0961540597861764D+00, &
          0.6970960567476024D+00, &
          0.0961540597861764D+00, &
          0.6970960567476024D+00, &
          0.2067498834662211D+00, &
          0.0117746553239110D+00, &
          0.3068968355307726D+00, &
          0.6813285091453164D+00, &
          0.3068968355307726D+00, &
          0.6813285091453164D+00, &
          0.0117746553239110D+00, &
          0.0117746553239110D+00, &
          0.3068968355307726D+00, &
          0.6813285091453164D+00, &
          0.3068968355307726D+00, &
          0.6813285091453164D+00, &
          0.0117746553239110D+00, &
          0.0174063356195718D+00, &
          0.0798486447852760D+00, &
          0.9027450195951523D+00, &
          0.0798486447852760D+00, &
          0.9027450195951523D+00, &
          0.0174063356195718D+00, &
          0.0174063356195718D+00, &
          0.0798486447852760D+00, &
          0.9027450195951523D+00, &
          0.0798486447852760D+00, &
          0.9027450195951523D+00, &
          0.0174063356195718D+00, &
          0.1618562004775664D+00, &
          0.0360777790627768D+00, &
          0.8020660204596567D+00, &
          0.0360777790627768D+00, &
          0.8020660204596567D+00, &
          0.1618562004775664D+00, &
          0.1618562004775664D+00, &
          0.0360777790627768D+00, &
          0.8020660204596567D+00, &
          0.0360777790627768D+00, &
          0.8020660204596567D+00, &
          0.1618562004775664D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.9534387216382612D+00, &
          0.9534387216382612D+00, &
          0.9534387216382612D+00, &
          0.0465612783617388D+00, &
          0.0465612783617388D+00, &
          0.0465612783617388D+00, &
          0.7595655027404185D+00, &
          0.7595655027404185D+00, &
          0.7595655027404185D+00, &
          0.2404344972595815D+00, &
          0.2404344972595815D+00, &
          0.2404344972595815D+00, &
          0.8487633516520783D+00, &
          0.8487633516520783D+00, &
          0.8487633516520783D+00, &
          0.1512366483479217D+00, &
          0.1512366483479217D+00, &
          0.1512366483479217D+00, &
          0.7352481798898231D+00, &
          0.7352481798898231D+00, &
          0.7352481798898231D+00, &
          0.2647518201101769D+00, &
          0.2647518201101769D+00, &
          0.2647518201101769D+00, &
          0.3114476852954629D+00, &
          0.3114476852954629D+00, &
          0.3114476852954629D+00, &
          0.6885523147045371D+00, &
          0.6885523147045371D+00, &
          0.6885523147045371D+00, &
          0.3645354584203186D+00, &
          0.3645354584203186D+00, &
          0.3645354584203186D+00, &
          0.6354645415796814D+00, &
          0.6354645415796814D+00, &
          0.6354645415796814D+00, &
          0.4121101005024969D+00, &
          0.4121101005024969D+00, &
          0.4121101005024969D+00, &
          0.5878898994975031D+00, &
          0.5878898994975031D+00, &
          0.5878898994975031D+00, &
          0.3548440062961510D+00, &
          0.3548440062961510D+00, &
          0.3548440062961510D+00, &
          0.6451559937038489D+00, &
          0.6451559937038489D+00, &
          0.6451559937038489D+00, &
          0.3501843143023031D+00, &
          0.3501843143023031D+00, &
          0.3501843143023031D+00, &
          0.6498156856976969D+00, &
          0.6498156856976969D+00, &
          0.6498156856976969D+00, &
          0.0656042965456428D+00, &
          0.0656042965456428D+00, &
          0.0656042965456428D+00, &
          0.9343957034543572D+00, &
          0.9343957034543572D+00, &
          0.9343957034543572D+00, &
          0.1448541902653249D+00, &
          0.1448541902653249D+00, &
          0.1448541902653249D+00, &
          0.8551458097346751D+00, &
          0.8551458097346751D+00, &
          0.8551458097346751D+00, &
          0.9806685335400167D+00, &
          0.9806685335400167D+00, &
          0.9806685335400167D+00, &
          0.0193314664599833D+00, &
          0.0193314664599833D+00, &
          0.0193314664599833D+00, &
          0.6290398587069839D+00, &
          0.6290398587069839D+00, &
          0.6290398587069839D+00, &
          0.3709601412930161D+00, &
          0.3709601412930161D+00, &
          0.3709601412930161D+00, &
          0.9287263768670988D+00, &
          0.9287263768670988D+00, &
          0.9287263768670988D+00, &
          0.0712736231329013D+00, &
          0.0712736231329013D+00, &
          0.0712736231329013D+00, &
          0.9721416113157295D+00, &
          0.9721416113157295D+00, &
          0.9721416113157295D+00, &
          0.0278583886842705D+00, &
          0.0278583886842705D+00, &
          0.0278583886842705D+00, &
          0.2369147251638602D+00, &
          0.2369147251638602D+00, &
          0.2369147251638602D+00, &
          0.7630852748361399D+00, &
          0.7630852748361399D+00, &
          0.7630852748361399D+00, &
          0.8919378093663139D+00, &
          0.8919378093663139D+00, &
          0.8919378093663139D+00, &
          0.1080621906336861D+00, &
          0.1080621906336861D+00, &
          0.1080621906336861D+00, &
          0.9842060169219227D+00, &
          0.9842060169219227D+00, &
          0.9842060169219227D+00, &
          0.0157939830780774D+00, &
          0.0157939830780774D+00, &
          0.0157939830780774D+00, &
          0.9527746747185892D+00, &
          0.9527746747185892D+00, &
          0.9527746747185892D+00, &
          0.0472253252814108D+00, &
          0.0472253252814108D+00, &
          0.0472253252814108D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9343231379314123D+00, &
          0.9343231379314123D+00, &
          0.9343231379314123D+00, &
          0.9343231379314123D+00, &
          0.9343231379314123D+00, &
          0.9343231379314123D+00, &
          0.0656768620685876D+00, &
          0.0656768620685876D+00, &
          0.0656768620685876D+00, &
          0.0656768620685876D+00, &
          0.0656768620685876D+00, &
          0.0656768620685876D+00, &
          0.5513722678894569D+00, &
          0.5513722678894569D+00, &
          0.5513722678894569D+00, &
          0.5513722678894569D+00, &
          0.5513722678894569D+00, &
          0.5513722678894569D+00, &
          0.4486277321105431D+00, &
          0.4486277321105431D+00, &
          0.4486277321105431D+00, &
          0.4486277321105431D+00, &
          0.4486277321105431D+00, &
          0.4486277321105431D+00, &
          0.2743627864784612D+00, &
          0.2743627864784612D+00, &
          0.2743627864784612D+00, &
          0.2743627864784612D+00, &
          0.2743627864784612D+00, &
          0.2743627864784612D+00, &
          0.7256372135215388D+00, &
          0.7256372135215388D+00, &
          0.7256372135215388D+00, &
          0.7256372135215388D+00, &
          0.7256372135215388D+00, &
          0.7256372135215388D+00, &
          0.2817575164863711D+00, &
          0.2817575164863711D+00, &
          0.2817575164863711D+00, &
          0.2817575164863711D+00, &
          0.2817575164863711D+00, &
          0.2817575164863711D+00, &
          0.7182424835136290D+00, &
          0.7182424835136290D+00, &
          0.7182424835136290D+00, &
          0.7182424835136290D+00, &
          0.7182424835136290D+00, &
          0.7182424835136290D+00, &
          0.5891160926744499D+00, &
          0.5891160926744499D+00, &
          0.5891160926744499D+00, &
          0.5891160926744499D+00, &
          0.5891160926744499D+00, &
          0.5891160926744499D+00, &
          0.4108839073255501D+00, &
          0.4108839073255501D+00, &
          0.4108839073255501D+00, &
          0.4108839073255501D+00, &
          0.4108839073255501D+00, &
          0.4108839073255501D+00, &
          0.1714235980915265D+00, &
          0.1714235980915265D+00, &
          0.1714235980915265D+00, &
          0.1714235980915265D+00, &
          0.1714235980915265D+00, &
          0.1714235980915265D+00, &
          0.8285764019084735D+00, &
          0.8285764019084735D+00, &
          0.8285764019084735D+00, &
          0.8285764019084735D+00, &
          0.8285764019084735D+00, &
          0.8285764019084735D+00, &
          0.9394301860037320D+00, &
          0.9394301860037320D+00, &
          0.9394301860037320D+00, &
          0.9394301860037320D+00, &
          0.9394301860037320D+00, &
          0.9394301860037320D+00, &
          0.0605698139962680D+00, &
          0.0605698139962680D+00, &
          0.0605698139962680D+00, &
          0.0605698139962680D+00, &
          0.0605698139962680D+00, &
          0.0605698139962680D+00, &
          0.8744417684033916D+00, &
          0.8744417684033916D+00, &
          0.8744417684033916D+00, &
          0.8744417684033916D+00, &
          0.8744417684033916D+00, &
          0.8744417684033916D+00, &
          0.1255582315966083D+00, &
          0.1255582315966083D+00, &
          0.1255582315966083D+00, &
          0.1255582315966083D+00, &
          0.1255582315966083D+00, &
          0.1255582315966083D+00, &
          0.4576004985221995D+00, &
          0.4576004985221995D+00, &
          0.4576004985221995D+00, &
          0.4576004985221995D+00, &
          0.4576004985221995D+00, &
          0.4576004985221995D+00, &
          0.5423995014778005D+00, &
          0.5423995014778005D+00, &
          0.5423995014778005D+00, &
          0.5423995014778005D+00, &
          0.5423995014778005D+00, &
          0.5423995014778005D+00, &
          0.5469802767947162D+00, &
          0.5469802767947162D+00, &
          0.5469802767947162D+00, &
          0.5469802767947162D+00, &
          0.5469802767947162D+00, &
          0.5469802767947162D+00, &
          0.4530197232052838D+00, &
          0.4530197232052838D+00, &
          0.4530197232052838D+00, &
          0.4530197232052838D+00, &
          0.4530197232052838D+00, &
          0.4530197232052838D+00, &
          0.7227042014132130D+00, &
          0.7227042014132130D+00, &
          0.7227042014132130D+00, &
          0.7227042014132130D+00, &
          0.7227042014132130D+00, &
          0.7227042014132130D+00, &
          0.2772957985867870D+00, &
          0.2772957985867870D+00, &
          0.2772957985867870D+00, &
          0.2772957985867870D+00, &
          0.2772957985867870D+00, &
          0.2772957985867870D+00, &
          0.6586462869796975D+00, &
          0.6586462869796975D+00, &
          0.6586462869796975D+00, &
          0.6586462869796975D+00, &
          0.6586462869796975D+00, &
          0.6586462869796975D+00, &
          0.3413537130203025D+00, &
          0.3413537130203025D+00, &
          0.3413537130203025D+00, &
          0.3413537130203025D+00, &
          0.3413537130203025D+00, &
          0.3413537130203025D+00, &
          0.9865056469526946D+00, &
          0.9865056469526946D+00, &
          0.9865056469526946D+00, &
          0.9865056469526946D+00, &
          0.9865056469526946D+00, &
          0.9865056469526946D+00, &
          0.0134943530473054D+00, &
          0.0134943530473054D+00, &
          0.0134943530473054D+00, &
          0.0134943530473054D+00, &
          0.0134943530473054D+00, &
          0.0134943530473054D+00, &
          0.5827285974639408D+00, &
          0.5827285974639408D+00, &
          0.5827285974639408D+00, &
          0.5827285974639408D+00, &
          0.5827285974639408D+00, &
          0.5827285974639408D+00, &
          0.4172714025360593D+00, &
          0.4172714025360593D+00, &
          0.4172714025360593D+00, &
          0.4172714025360593D+00, &
          0.4172714025360593D+00, &
          0.4172714025360593D+00, &
          0.8529096418585074D+00, &
          0.8529096418585074D+00, &
          0.8529096418585074D+00, &
          0.8529096418585074D+00, &
          0.8529096418585074D+00, &
          0.8529096418585074D+00, &
          0.1470903581414927D+00, &
          0.1470903581414927D+00, &
          0.1470903581414927D+00, &
          0.1470903581414927D+00, &
          0.1470903581414927D+00, &
          0.1470903581414927D+00, &
          0.8283291277845266D+00, &
          0.8283291277845266D+00, &
          0.8283291277845266D+00, &
          0.8283291277845266D+00, &
          0.8283291277845266D+00, &
          0.8283291277845266D+00, &
          0.1716708722154733D+00, &
          0.1716708722154733D+00, &
          0.1716708722154733D+00, &
          0.1716708722154733D+00, &
          0.1716708722154733D+00, &
          0.1716708722154733D+00, &
          0.9370173839969559D+00, &
          0.9370173839969559D+00, &
          0.9370173839969559D+00, &
          0.9370173839969559D+00, &
          0.9370173839969559D+00, &
          0.9370173839969559D+00, &
          0.0629826160030441D+00, &
          0.0629826160030441D+00, &
          0.0629826160030441D+00, &
          0.0629826160030441D+00, &
          0.0629826160030441D+00, &
          0.0629826160030441D+00, &
          0.9939758307965130D+00, &
          0.9939758307965130D+00, &
          0.9939758307965130D+00, &
          0.9939758307965130D+00, &
          0.9939758307965130D+00, &
          0.9939758307965130D+00, &
          0.0060241692034870D+00, &
          0.0060241692034870D+00, &
          0.0060241692034870D+00, &
          0.0060241692034870D+00, &
          0.0060241692034870D+00, &
          0.0060241692034870D+00, &
          0.5786001946994722D+00, &
          0.5786001946994722D+00, &
          0.5786001946994722D+00, &
          0.5786001946994722D+00, &
          0.5786001946994722D+00, &
          0.5786001946994722D+00, &
          0.4213998053005278D+00, &
          0.4213998053005278D+00, &
          0.4213998053005278D+00, &
          0.4213998053005278D+00, &
          0.4213998053005278D+00, &
          0.4213998053005278D+00, &
          0.8755215159278712D+00, &
          0.8755215159278712D+00, &
          0.8755215159278712D+00, &
          0.8755215159278712D+00, &
          0.8755215159278712D+00, &
          0.8755215159278712D+00, &
          0.1244784840721288D+00, &
          0.1244784840721288D+00, &
          0.1244784840721288D+00, &
          0.1244784840721288D+00, &
          0.1244784840721288D+00, &
          0.1244784840721288D+00, &
          0.9893261490512200D+00, &
          0.9893261490512200D+00, &
          0.9893261490512200D+00, &
          0.9893261490512200D+00, &
          0.9893261490512200D+00, &
          0.9893261490512200D+00, &
          0.0106738509487800D+00, &
          0.0106738509487800D+00, &
          0.0106738509487800D+00, &
          0.0106738509487800D+00, &
          0.0106738509487800D+00, &
          0.0106738509487800D+00, &
          0.9893821585030006D+00, &
          0.9893821585030006D+00, &
          0.9893821585030006D+00, &
          0.9893821585030006D+00, &
          0.9893821585030006D+00, &
          0.9893821585030006D+00, &
          0.0106178414969994D+00, &
          0.0106178414969994D+00, &
          0.0106178414969994D+00, &
          0.0106178414969994D+00, &
          0.0106178414969994D+00, &
          0.0106178414969994D+00, &
          0.7966282776345323D+00, &
          0.7966282776345323D+00, &
          0.7966282776345323D+00, &
          0.7966282776345323D+00, &
          0.7966282776345323D+00, &
          0.7966282776345323D+00, &
          0.2033717223654677D+00, &
          0.2033717223654677D+00, &
          0.2033717223654677D+00, &
          0.2033717223654677D+00, &
          0.2033717223654677D+00, &
          0.2033717223654677D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          0.0022576928480960D+00, &
          0.0022576928480960D+00, &
          0.0022576928480960D+00, &
          0.0022576928480960D+00, &
          0.0022576928480960D+00, &
          0.0022576928480960D+00, &
          0.0035220009936119D+00, &
          0.0035220009936119D+00, &
          0.0035220009936119D+00, &
          0.0035220009936119D+00, &
          0.0035220009936119D+00, &
          0.0035220009936119D+00, &
          0.0044617234977870D+00, &
          0.0044617234977870D+00, &
          0.0044617234977870D+00, &
          0.0044617234977870D+00, &
          0.0044617234977870D+00, &
          0.0044617234977870D+00, &
          0.0003707999359392D+00, &
          0.0003707999359392D+00, &
          0.0003707999359392D+00, &
          0.0003707999359392D+00, &
          0.0003707999359392D+00, &
          0.0003707999359392D+00, &
          0.0024236611427218D+00, &
          0.0024236611427218D+00, &
          0.0024236611427218D+00, &
          0.0024236611427218D+00, &
          0.0024236611427218D+00, &
          0.0024236611427218D+00, &
          0.0028149967000162D+00, &
          0.0028149967000162D+00, &
          0.0028149967000162D+00, &
          0.0028149967000162D+00, &
          0.0028149967000162D+00, &
          0.0028149967000162D+00, &
          0.0053376421669576D+00, &
          0.0053376421669576D+00, &
          0.0053376421669576D+00, &
          0.0053376421669576D+00, &
          0.0053376421669576D+00, &
          0.0053376421669576D+00, &
          0.0085742666644647D+00, &
          0.0085742666644647D+00, &
          0.0085742666644647D+00, &
          0.0085742666644647D+00, &
          0.0085742666644647D+00, &
          0.0085742666644647D+00, &
          0.0022657258714867D+00, &
          0.0022657258714867D+00, &
          0.0022657258714867D+00, &
          0.0022657258714867D+00, &
          0.0022657258714867D+00, &
          0.0022657258714867D+00, &
          0.0029551579597074D+00, &
          0.0029551579597074D+00, &
          0.0029551579597074D+00, &
          0.0029551579597074D+00, &
          0.0029551579597074D+00, &
          0.0029551579597074D+00, &
          0.0039194664968460D+00, &
          0.0039194664968460D+00, &
          0.0039194664968460D+00, &
          0.0039194664968460D+00, &
          0.0039194664968460D+00, &
          0.0039194664968460D+00, &
          0.0006030674708196D+00, &
          0.0006030674708196D+00, &
          0.0006030674708196D+00, &
          0.0006030674708196D+00, &
          0.0006030674708196D+00, &
          0.0006030674708196D+00, &
          0.0045779809945426D+00, &
          0.0045779809945426D+00, &
          0.0045779809945426D+00, &
          0.0045779809945426D+00, &
          0.0045779809945426D+00, &
          0.0045779809945426D+00, &
          0.0024306636423047D+00, &
          0.0024306636423047D+00, &
          0.0024306636423047D+00, &
          0.0024306636423047D+00, &
          0.0024306636423047D+00, &
          0.0024306636423047D+00, &
          0.0021599699258395D+00, &
          0.0021599699258395D+00, &
          0.0021599699258395D+00, &
          0.0021599699258395D+00, &
          0.0021599699258395D+00, &
          0.0021599699258395D+00, &
          0.0079010510892335D+00, &
          0.0079010510892335D+00, &
          0.0079010510892335D+00, &
          0.0079010510892335D+00, &
          0.0079010510892335D+00, &
          0.0079010510892335D+00, &
          0.0021738453915315D+00, &
          0.0021738453915315D+00, &
          0.0021738453915315D+00, &
          0.0021738453915315D+00, &
          0.0021738453915315D+00, &
          0.0021738453915315D+00, &
          0.0022486463151748D+00, &
          0.0022486463151748D+00, &
          0.0022486463151748D+00, &
          0.0022486463151748D+00, &
          0.0022486463151748D+00, &
          0.0022486463151748D+00, &
          0.0003442548195872D+00, &
          0.0003442548195872D+00, &
          0.0003442548195872D+00, &
          0.0003442548195872D+00, &
          0.0003442548195872D+00, &
          0.0003442548195872D+00, &
          0.0080699355816128D+00, &
          0.0080699355816128D+00, &
          0.0080699355816128D+00, &
          0.0080699355816128D+00, &
          0.0080699355816128D+00, &
          0.0080699355816128D+00, &
          0.0023877231635901D+00, &
          0.0023877231635901D+00, &
          0.0023877231635901D+00, &
          0.0023877231635901D+00, &
          0.0023877231635901D+00, &
          0.0023877231635901D+00, &
          0.0023877231635901D+00, &
          0.0023877231635901D+00, &
          0.0023877231635901D+00, &
          0.0023877231635901D+00, &
          0.0023877231635901D+00, &
          0.0023877231635901D+00, &
          0.0034305282210409D+00, &
          0.0034305282210409D+00, &
          0.0034305282210409D+00, &
          0.0034305282210409D+00, &
          0.0034305282210409D+00, &
          0.0034305282210409D+00, &
          0.0034305282210409D+00, &
          0.0034305282210409D+00, &
          0.0034305282210409D+00, &
          0.0034305282210409D+00, &
          0.0034305282210409D+00, &
          0.0034305282210409D+00, &
          0.0060295453672165D+00, &
          0.0060295453672165D+00, &
          0.0060295453672165D+00, &
          0.0060295453672165D+00, &
          0.0060295453672165D+00, &
          0.0060295453672165D+00, &
          0.0060295453672165D+00, &
          0.0060295453672165D+00, &
          0.0060295453672165D+00, &
          0.0060295453672165D+00, &
          0.0060295453672165D+00, &
          0.0060295453672165D+00, &
          0.0032768238791581D+00, &
          0.0032768238791581D+00, &
          0.0032768238791581D+00, &
          0.0032768238791581D+00, &
          0.0032768238791581D+00, &
          0.0032768238791581D+00, &
          0.0032768238791581D+00, &
          0.0032768238791581D+00, &
          0.0032768238791581D+00, &
          0.0032768238791581D+00, &
          0.0032768238791581D+00, &
          0.0032768238791581D+00, &
          0.0036635492105407D+00, &
          0.0036635492105407D+00, &
          0.0036635492105407D+00, &
          0.0036635492105407D+00, &
          0.0036635492105407D+00, &
          0.0036635492105407D+00, &
          0.0036635492105407D+00, &
          0.0036635492105407D+00, &
          0.0036635492105407D+00, &
          0.0036635492105407D+00, &
          0.0036635492105407D+00, &
          0.0036635492105407D+00, &
          0.0044244882628077D+00, &
          0.0044244882628077D+00, &
          0.0044244882628077D+00, &
          0.0044244882628077D+00, &
          0.0044244882628077D+00, &
          0.0044244882628077D+00, &
          0.0044244882628077D+00, &
          0.0044244882628077D+00, &
          0.0044244882628077D+00, &
          0.0044244882628077D+00, &
          0.0044244882628077D+00, &
          0.0044244882628077D+00, &
          0.0027600635232146D+00, &
          0.0027600635232146D+00, &
          0.0027600635232146D+00, &
          0.0027600635232146D+00, &
          0.0027600635232146D+00, &
          0.0027600635232146D+00, &
          0.0027600635232146D+00, &
          0.0027600635232146D+00, &
          0.0027600635232146D+00, &
          0.0027600635232146D+00, &
          0.0027600635232146D+00, &
          0.0027600635232146D+00, &
          0.0013921098894439D+00, &
          0.0013921098894439D+00, &
          0.0013921098894439D+00, &
          0.0013921098894439D+00, &
          0.0013921098894439D+00, &
          0.0013921098894439D+00, &
          0.0013921098894439D+00, &
          0.0013921098894439D+00, &
          0.0013921098894439D+00, &
          0.0013921098894439D+00, &
          0.0013921098894439D+00, &
          0.0013921098894439D+00, &
          0.0021384358726425D+00, &
          0.0021384358726425D+00, &
          0.0021384358726425D+00, &
          0.0021384358726425D+00, &
          0.0021384358726425D+00, &
          0.0021384358726425D+00, &
          0.0021384358726425D+00, &
          0.0021384358726425D+00, &
          0.0021384358726425D+00, &
          0.0021384358726425D+00, &
          0.0021384358726425D+00, &
          0.0021384358726425D+00, &
          0.0005166067073489D+00, &
          0.0005166067073489D+00, &
          0.0005166067073489D+00, &
          0.0005166067073489D+00, &
          0.0005166067073489D+00, &
          0.0005166067073489D+00, &
          0.0005166067073489D+00, &
          0.0005166067073489D+00, &
          0.0005166067073489D+00, &
          0.0005166067073489D+00, &
          0.0005166067073489D+00, &
          0.0005166067073489D+00, &
          0.0019539531411333D+00, &
          0.0019539531411333D+00, &
          0.0019539531411333D+00, &
          0.0019539531411333D+00, &
          0.0019539531411333D+00, &
          0.0019539531411333D+00, &
          0.0019539531411333D+00, &
          0.0019539531411333D+00, &
          0.0019539531411333D+00, &
          0.0019539531411333D+00, &
          0.0019539531411333D+00, &
          0.0019539531411333D+00, &
          0.0012550085983638D+00, &
          0.0012550085983638D+00, &
          0.0012550085983638D+00, &
          0.0012550085983638D+00, &
          0.0012550085983638D+00, &
          0.0012550085983638D+00, &
          0.0012550085983638D+00, &
          0.0012550085983638D+00, &
          0.0012550085983638D+00, &
          0.0012550085983638D+00, &
          0.0012550085983638D+00, &
          0.0012550085983638D+00, &
          0.0010423341141818D+00, &
          0.0010423341141818D+00, &
          0.0010423341141818D+00, &
          0.0010423341141818D+00, &
          0.0010423341141818D+00, &
          0.0010423341141818D+00, &
          0.0010423341141818D+00, &
          0.0010423341141818D+00, &
          0.0010423341141818D+00, &
          0.0010423341141818D+00, &
          0.0010423341141818D+00, &
          0.0010423341141818D+00, &
          0.0012793420418862D+00, &
          0.0012793420418862D+00, &
          0.0012793420418862D+00, &
          0.0012793420418862D+00, &
          0.0012793420418862D+00, &
          0.0012793420418862D+00, &
          0.0012793420418862D+00, &
          0.0012793420418862D+00, &
          0.0012793420418862D+00, &
          0.0012793420418862D+00, &
          0.0012793420418862D+00, &
          0.0012793420418862D+00, &
          0.0010877722341374D+00, &
          0.0010877722341374D+00, &
          0.0010877722341374D+00, &
          0.0010877722341374D+00, &
          0.0010877722341374D+00, &
          0.0010877722341374D+00, &
          0.0010877722341374D+00, &
          0.0010877722341374D+00, &
          0.0010877722341374D+00, &
          0.0010877722341374D+00, &
          0.0010877722341374D+00, &
          0.0010877722341374D+00, &
          0.0008506443220554D+00, &
          0.0008506443220554D+00, &
          0.0008506443220554D+00, &
          0.0008506443220554D+00, &
          0.0008506443220554D+00, &
          0.0008506443220554D+00, &
          0.0008506443220554D+00, &
          0.0008506443220554D+00, &
          0.0008506443220554D+00, &
          0.0008506443220554D+00, &
          0.0008506443220554D+00, &
          0.0008506443220554D+00, &
          0.0010717384362760D+00, &
          0.0010717384362760D+00, &
          0.0010717384362760D+00, &
          0.0010717384362760D+00, &
          0.0010717384362760D+00, &
          0.0010717384362760D+00, &
          0.0010717384362760D+00, &
          0.0010717384362760D+00, &
          0.0010717384362760D+00, &
          0.0010717384362760D+00, &
          0.0010717384362760D+00, &
          0.0010717384362760D+00, &
          0.0012223898383723D+00, &
          0.0012223898383723D+00, &
          0.0012223898383723D+00, &
          0.0012223898383723D+00, &
          0.0012223898383723D+00, &
          0.0012223898383723D+00, &
          0.0012223898383723D+00, &
          0.0012223898383723D+00, &
          0.0012223898383723D+00, &
          0.0012223898383723D+00, &
          0.0012223898383723D+00, &
          0.0012223898383723D+00, &
          0.0014264821649004D+00, &
          0.0014264821649004D+00, &
          0.0014264821649004D+00, &
          0.0014264821649004D+00, &
          0.0014264821649004D+00, &
          0.0014264821649004D+00, &
          0.0014264821649004D+00, &
          0.0014264821649004D+00, &
          0.0014264821649004D+00, &
          0.0014264821649004D+00, &
          0.0014264821649004D+00, &
          0.0014264821649004D+00, &
          0.0038341778997734D+00, &
          0.0038341778997734D+00, &
          0.0038341778997734D+00, &
          0.0038341778997734D+00, &
          0.0038341778997734D+00, &
          0.0038341778997734D+00, &
          0.0038341778997734D+00, &
          0.0038341778997734D+00, &
          0.0038341778997734D+00, &
          0.0038341778997734D+00, &
          0.0038341778997734D+00, &
          0.0038341778997734D+00, &
          0.0004758021594541D+00, &
          0.0004758021594541D+00, &
          0.0004758021594541D+00, &
          0.0004758021594541D+00, &
          0.0004758021594541D+00, &
          0.0004758021594541D+00, &
          0.0004758021594541D+00, &
          0.0004758021594541D+00, &
          0.0004758021594541D+00, &
          0.0004758021594541D+00, &
          0.0004758021594541D+00, &
          0.0004758021594541D+00, &
          0.0003711298716255D+00, &
          0.0003711298716255D+00, &
          0.0003711298716255D+00, &
          0.0003711298716255D+00, &
          0.0003711298716255D+00, &
          0.0003711298716255D+00, &
          0.0003711298716255D+00, &
          0.0003711298716255D+00, &
          0.0003711298716255D+00, &
          0.0003711298716255D+00, &
          0.0003711298716255D+00, &
          0.0003711298716255D+00, &
          0.0027364096600290D+00, &
          0.0027364096600290D+00, &
          0.0027364096600290D+00, &
          0.0027364096600290D+00, &
          0.0027364096600290D+00, &
          0.0027364096600290D+00, &
          0.0027364096600290D+00, &
          0.0027364096600290D+00, &
          0.0027364096600290D+00, &
          0.0027364096600290D+00, &
          0.0027364096600290D+00, &
          0.0027364096600290D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule19 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule19() returns the prism rule of precision 19.
!
!  Discussion:
!
!    We suppose we are given a triangular prism P with vertices
!    (1,0,0), (0,1,0), (0,0,0),
!    (1,0,1), (0,1,1), (0,0,1),
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over P is 
!    approximated by Q as follows:
!
!    Q = volume(P) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 May 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order symmetric cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 122, Number 1, pages 148-171, 24 August 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real x(n), y(n), z(n): the coordinates of quadrature points.
!
!    real w(n): the quadrature weights, which sum to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: n_save = 420

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.4935398626809867D+00, &
          0.4935398626809867D+00, &
          0.0129202746380266D+00, &
          0.2767659367008514D+00, &
          0.2767659367008514D+00, &
          0.4464681265982973D+00, &
          0.0545777948502843D+00, &
          0.0545777948502843D+00, &
          0.8908444102994315D+00, &
          0.0139969735218886D+00, &
          0.0139969735218886D+00, &
          0.9720060529562228D+00, &
          0.1418740786751734D+00, &
          0.1418740786751734D+00, &
          0.7162518426496531D+00, &
          0.1418740786751734D+00, &
          0.1418740786751734D+00, &
          0.7162518426496531D+00, &
          0.0083955990101041D+00, &
          0.0083955990101041D+00, &
          0.9832088019797917D+00, &
          0.0083955990101041D+00, &
          0.0083955990101041D+00, &
          0.9832088019797917D+00, &
          0.1035973359771996D+00, &
          0.1035973359771996D+00, &
          0.7928053280456009D+00, &
          0.1035973359771996D+00, &
          0.1035973359771996D+00, &
          0.7928053280456009D+00, &
          0.3956080996828434D+00, &
          0.3956080996828434D+00, &
          0.2087838006343133D+00, &
          0.3956080996828434D+00, &
          0.3956080996828434D+00, &
          0.2087838006343133D+00, &
          0.0140019311977436D+00, &
          0.0140019311977436D+00, &
          0.9719961376045129D+00, &
          0.0140019311977436D+00, &
          0.0140019311977436D+00, &
          0.9719961376045129D+00, &
          0.4801897764603656D+00, &
          0.4801897764603656D+00, &
          0.0396204470792688D+00, &
          0.4801897764603656D+00, &
          0.4801897764603656D+00, &
          0.0396204470792688D+00, &
          0.1623484102426076D+00, &
          0.1623484102426076D+00, &
          0.6753031795147848D+00, &
          0.1623484102426076D+00, &
          0.1623484102426076D+00, &
          0.6753031795147848D+00, &
          0.4972236736522475D+00, &
          0.4972236736522475D+00, &
          0.0055526526955051D+00, &
          0.4972236736522475D+00, &
          0.4972236736522475D+00, &
          0.0055526526955051D+00, &
          0.1990216655966096D+00, &
          0.1990216655966096D+00, &
          0.6019566688067807D+00, &
          0.1990216655966096D+00, &
          0.1990216655966096D+00, &
          0.6019566688067807D+00, &
          0.1202889101503504D+00, &
          0.1202889101503504D+00, &
          0.7594221796992991D+00, &
          0.1202889101503504D+00, &
          0.1202889101503504D+00, &
          0.7594221796992991D+00, &
          0.0532191526098831D+00, &
          0.0532191526098831D+00, &
          0.8935616947802337D+00, &
          0.0532191526098831D+00, &
          0.0532191526098831D+00, &
          0.8935616947802337D+00, &
          0.4012796612843944D+00, &
          0.4012796612843944D+00, &
          0.1974406774312114D+00, &
          0.4012796612843944D+00, &
          0.4012796612843944D+00, &
          0.1974406774312114D+00, &
          0.1089883669604236D+00, &
          0.1089883669604236D+00, &
          0.7820232660791528D+00, &
          0.1089883669604236D+00, &
          0.1089883669604236D+00, &
          0.7820232660791528D+00, &
          0.4472748897677901D+00, &
          0.4472748897677901D+00, &
          0.1054502204644198D+00, &
          0.4472748897677901D+00, &
          0.4472748897677901D+00, &
          0.1054502204644198D+00, &
          0.0692335448571578D+00, &
          0.0692335448571578D+00, &
          0.8615329102856844D+00, &
          0.0692335448571578D+00, &
          0.0692335448571578D+00, &
          0.8615329102856844D+00, &
          0.2475145234973103D+00, &
          0.2475145234973103D+00, &
          0.5049709530053793D+00, &
          0.2475145234973103D+00, &
          0.2475145234973103D+00, &
          0.5049709530053793D+00, &
          0.4533377451558121D+00, &
          0.4533377451558121D+00, &
          0.0933245096883758D+00, &
          0.4533377451558121D+00, &
          0.4533377451558121D+00, &
          0.0933245096883758D+00, &
          0.0180837411268107D+00, &
          0.0180837411268107D+00, &
          0.9638325177463786D+00, &
          0.0180837411268107D+00, &
          0.0180837411268107D+00, &
          0.9638325177463786D+00, &
          0.1297897131077726D+00, &
          0.0069174965695053D+00, &
          0.1297897131077726D+00, &
          0.8632927903227221D+00, &
          0.0069174965695053D+00, &
          0.8632927903227221D+00, &
          0.3779418721990778D+00, &
          0.0310061742359910D+00, &
          0.3779418721990778D+00, &
          0.5910519535649312D+00, &
          0.0310061742359910D+00, &
          0.5910519535649312D+00, &
          0.0644631311955231D+00, &
          0.2203251202566757D+00, &
          0.0644631311955231D+00, &
          0.7152117485478011D+00, &
          0.2203251202566757D+00, &
          0.7152117485478011D+00, &
          0.1625899734623090D+00, &
          0.2824754920851798D+00, &
          0.1625899734623090D+00, &
          0.5549345344525113D+00, &
          0.2824754920851798D+00, &
          0.5549345344525113D+00, &
          0.1625899734623090D+00, &
          0.2824754920851798D+00, &
          0.1625899734623090D+00, &
          0.5549345344525113D+00, &
          0.2824754920851798D+00, &
          0.5549345344525113D+00, &
          0.0081184148094398D+00, &
          0.0547883511283449D+00, &
          0.0081184148094398D+00, &
          0.9370932340622153D+00, &
          0.0547883511283449D+00, &
          0.9370932340622153D+00, &
          0.0081184148094398D+00, &
          0.0547883511283449D+00, &
          0.0081184148094398D+00, &
          0.9370932340622153D+00, &
          0.0547883511283449D+00, &
          0.9370932340622153D+00, &
          0.6582273496659641D+00, &
          0.1178005364736766D+00, &
          0.6582273496659641D+00, &
          0.2239721138603593D+00, &
          0.1178005364736766D+00, &
          0.2239721138603593D+00, &
          0.6582273496659641D+00, &
          0.1178005364736766D+00, &
          0.6582273496659641D+00, &
          0.2239721138603593D+00, &
          0.1178005364736766D+00, &
          0.2239721138603593D+00, &
          0.2876874815978630D+00, &
          0.0355268541454468D+00, &
          0.2876874815978630D+00, &
          0.6767856642566902D+00, &
          0.0355268541454468D+00, &
          0.6767856642566902D+00, &
          0.2876874815978630D+00, &
          0.0355268541454468D+00, &
          0.2876874815978630D+00, &
          0.6767856642566902D+00, &
          0.0355268541454468D+00, &
          0.6767856642566902D+00, &
          0.2839484400308330D+00, &
          0.2073917255967352D+00, &
          0.2839484400308330D+00, &
          0.5086598343724318D+00, &
          0.2073917255967352D+00, &
          0.5086598343724318D+00, &
          0.2839484400308330D+00, &
          0.2073917255967352D+00, &
          0.2839484400308330D+00, &
          0.5086598343724318D+00, &
          0.2073917255967352D+00, &
          0.5086598343724318D+00, &
          0.0282457003834953D+00, &
          0.0862230026956853D+00, &
          0.0282457003834953D+00, &
          0.8855312969208194D+00, &
          0.0862230026956853D+00, &
          0.8855312969208194D+00, &
          0.0282457003834953D+00, &
          0.0862230026956853D+00, &
          0.0282457003834953D+00, &
          0.8855312969208194D+00, &
          0.0862230026956853D+00, &
          0.8855312969208194D+00, &
          0.3915043117238656D+00, &
          0.0056196516036861D+00, &
          0.3915043117238656D+00, &
          0.6028760366724483D+00, &
          0.0056196516036861D+00, &
          0.6028760366724483D+00, &
          0.3915043117238656D+00, &
          0.0056196516036861D+00, &
          0.3915043117238656D+00, &
          0.6028760366724483D+00, &
          0.0056196516036861D+00, &
          0.6028760366724483D+00, &
          0.0737667613724345D+00, &
          0.3389265718024346D+00, &
          0.0737667613724345D+00, &
          0.5873066668251310D+00, &
          0.3389265718024346D+00, &
          0.5873066668251310D+00, &
          0.0737667613724345D+00, &
          0.3389265718024346D+00, &
          0.0737667613724345D+00, &
          0.5873066668251310D+00, &
          0.3389265718024346D+00, &
          0.5873066668251310D+00, &
          0.1582419524500696D+00, &
          0.0080205188978798D+00, &
          0.1582419524500696D+00, &
          0.8337375286520506D+00, &
          0.0080205188978798D+00, &
          0.8337375286520506D+00, &
          0.1582419524500696D+00, &
          0.0080205188978798D+00, &
          0.1582419524500696D+00, &
          0.8337375286520506D+00, &
          0.0080205188978798D+00, &
          0.8337375286520506D+00, &
          0.2270465393157224D+00, &
          0.0761183063986999D+00, &
          0.2270465393157224D+00, &
          0.6968351542855776D+00, &
          0.0761183063986999D+00, &
          0.6968351542855776D+00, &
          0.2270465393157224D+00, &
          0.0761183063986999D+00, &
          0.2270465393157224D+00, &
          0.6968351542855776D+00, &
          0.0761183063986999D+00, &
          0.6968351542855776D+00, &
          0.0107402320272877D+00, &
          0.0674860196923538D+00, &
          0.0107402320272877D+00, &
          0.9217737482803585D+00, &
          0.0674860196923538D+00, &
          0.9217737482803585D+00, &
          0.0107402320272877D+00, &
          0.0674860196923538D+00, &
          0.0107402320272877D+00, &
          0.9217737482803585D+00, &
          0.0674860196923538D+00, &
          0.9217737482803585D+00, &
          0.0450629754431290D+00, &
          0.1495289831435912D+00, &
          0.0450629754431290D+00, &
          0.8054080414132798D+00, &
          0.1495289831435912D+00, &
          0.8054080414132798D+00, &
          0.0450629754431290D+00, &
          0.1495289831435912D+00, &
          0.0450629754431290D+00, &
          0.8054080414132798D+00, &
          0.1495289831435912D+00, &
          0.8054080414132798D+00, &
          0.3899263834124326D+00, &
          0.0536599901206231D+00, &
          0.3899263834124326D+00, &
          0.5564136264669443D+00, &
          0.0536599901206231D+00, &
          0.5564136264669443D+00, &
          0.3899263834124326D+00, &
          0.0536599901206231D+00, &
          0.3899263834124326D+00, &
          0.5564136264669443D+00, &
          0.0536599901206231D+00, &
          0.5564136264669443D+00, &
          0.0142616821078597D+00, &
          0.4039638109964547D+00, &
          0.0142616821078597D+00, &
          0.5817745068956857D+00, &
          0.4039638109964547D+00, &
          0.5817745068956857D+00, &
          0.0142616821078597D+00, &
          0.4039638109964547D+00, &
          0.0142616821078597D+00, &
          0.5817745068956857D+00, &
          0.4039638109964547D+00, &
          0.5817745068956857D+00, &
          0.3314062155351160D+00, &
          0.0088330831166114D+00, &
          0.3314062155351160D+00, &
          0.6597607013482726D+00, &
          0.0088330831166114D+00, &
          0.6597607013482726D+00, &
          0.3314062155351160D+00, &
          0.0088330831166114D+00, &
          0.3314062155351160D+00, &
          0.6597607013482726D+00, &
          0.0088330831166114D+00, &
          0.6597607013482726D+00, &
          0.1223644156202090D+00, &
          0.0103923793593637D+00, &
          0.1223644156202090D+00, &
          0.8672432050204273D+00, &
          0.0103923793593637D+00, &
          0.8672432050204273D+00, &
          0.1223644156202090D+00, &
          0.0103923793593637D+00, &
          0.1223644156202090D+00, &
          0.8672432050204273D+00, &
          0.0103923793593637D+00, &
          0.8672432050204273D+00, &
          0.2929652829983945D+00, &
          0.1324448044285512D+00, &
          0.2929652829983945D+00, &
          0.5745899125730543D+00, &
          0.1324448044285512D+00, &
          0.5745899125730543D+00, &
          0.2929652829983945D+00, &
          0.1324448044285512D+00, &
          0.2929652829983945D+00, &
          0.5745899125730543D+00, &
          0.1324448044285512D+00, &
          0.5745899125730543D+00, &
          0.2534748953911949D+00, &
          0.0535147091599154D+00, &
          0.2534748953911949D+00, &
          0.6930103954488898D+00, &
          0.0535147091599154D+00, &
          0.6930103954488898D+00, &
          0.2534748953911949D+00, &
          0.0535147091599154D+00, &
          0.2534748953911949D+00, &
          0.6930103954488898D+00, &
          0.0535147091599154D+00, &
          0.6930103954488898D+00, &
          0.0440785872085001D+00, &
          0.1461188761758480D+00, &
          0.0440785872085001D+00, &
          0.8098025366156519D+00, &
          0.1461188761758480D+00, &
          0.8098025366156519D+00, &
          0.0440785872085001D+00, &
          0.1461188761758480D+00, &
          0.0440785872085001D+00, &
          0.8098025366156519D+00, &
          0.1461188761758480D+00, &
          0.8098025366156519D+00, &
          0.5031457633505478D+00, &
          0.1264604222462178D+00, &
          0.5031457633505478D+00, &
          0.3703938144032344D+00, &
          0.1264604222462178D+00, &
          0.3703938144032344D+00, &
          0.5031457633505478D+00, &
          0.1264604222462178D+00, &
          0.5031457633505478D+00, &
          0.3703938144032344D+00, &
          0.1264604222462178D+00, &
          0.3703938144032344D+00, &
          0.0087499979841172D+00, &
          0.2433623618774898D+00, &
          0.0087499979841172D+00, &
          0.7478876401383930D+00, &
          0.2433623618774898D+00, &
          0.7478876401383930D+00, &
          0.0087499979841172D+00, &
          0.2433623618774898D+00, &
          0.0087499979841172D+00, &
          0.7478876401383930D+00, &
          0.2433623618774898D+00, &
          0.7478876401383930D+00, &
          0.4433279475253553D+00, &
          0.2569581578007861D+00, &
          0.4433279475253553D+00, &
          0.2997138946738586D+00, &
          0.2569581578007861D+00, &
          0.2997138946738586D+00, &
          0.4433279475253553D+00, &
          0.2569581578007861D+00, &
          0.4433279475253553D+00, &
          0.2997138946738586D+00, &
          0.2569581578007861D+00, &
          0.2997138946738586D+00, &
          0.2485412372956145D+00, &
          0.0112771540337380D+00, &
          0.2485412372956145D+00, &
          0.7401816086706475D+00, &
          0.0112771540337380D+00, &
          0.7401816086706475D+00, &
          0.2485412372956145D+00, &
          0.0112771540337380D+00, &
          0.2485412372956145D+00, &
          0.7401816086706475D+00, &
          0.0112771540337380D+00, &
          0.7401816086706475D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.4935398626809867D+00, &
          0.0129202746380266D+00, &
          0.4935398626809867D+00, &
          0.2767659367008514D+00, &
          0.4464681265982973D+00, &
          0.2767659367008514D+00, &
          0.0545777948502843D+00, &
          0.8908444102994315D+00, &
          0.0545777948502843D+00, &
          0.0139969735218886D+00, &
          0.9720060529562228D+00, &
          0.0139969735218886D+00, &
          0.1418740786751734D+00, &
          0.7162518426496531D+00, &
          0.1418740786751734D+00, &
          0.1418740786751734D+00, &
          0.7162518426496531D+00, &
          0.1418740786751734D+00, &
          0.0083955990101041D+00, &
          0.9832088019797917D+00, &
          0.0083955990101041D+00, &
          0.0083955990101041D+00, &
          0.9832088019797917D+00, &
          0.0083955990101041D+00, &
          0.1035973359771996D+00, &
          0.7928053280456009D+00, &
          0.1035973359771996D+00, &
          0.1035973359771996D+00, &
          0.7928053280456009D+00, &
          0.1035973359771996D+00, &
          0.3956080996828434D+00, &
          0.2087838006343133D+00, &
          0.3956080996828434D+00, &
          0.3956080996828434D+00, &
          0.2087838006343133D+00, &
          0.3956080996828434D+00, &
          0.0140019311977436D+00, &
          0.9719961376045129D+00, &
          0.0140019311977436D+00, &
          0.0140019311977436D+00, &
          0.9719961376045129D+00, &
          0.0140019311977436D+00, &
          0.4801897764603656D+00, &
          0.0396204470792688D+00, &
          0.4801897764603656D+00, &
          0.4801897764603656D+00, &
          0.0396204470792688D+00, &
          0.4801897764603656D+00, &
          0.1623484102426076D+00, &
          0.6753031795147848D+00, &
          0.1623484102426076D+00, &
          0.1623484102426076D+00, &
          0.6753031795147848D+00, &
          0.1623484102426076D+00, &
          0.4972236736522475D+00, &
          0.0055526526955051D+00, &
          0.4972236736522475D+00, &
          0.4972236736522475D+00, &
          0.0055526526955051D+00, &
          0.4972236736522475D+00, &
          0.1990216655966096D+00, &
          0.6019566688067807D+00, &
          0.1990216655966096D+00, &
          0.1990216655966096D+00, &
          0.6019566688067807D+00, &
          0.1990216655966096D+00, &
          0.1202889101503504D+00, &
          0.7594221796992991D+00, &
          0.1202889101503504D+00, &
          0.1202889101503504D+00, &
          0.7594221796992991D+00, &
          0.1202889101503504D+00, &
          0.0532191526098831D+00, &
          0.8935616947802337D+00, &
          0.0532191526098831D+00, &
          0.0532191526098831D+00, &
          0.8935616947802337D+00, &
          0.0532191526098831D+00, &
          0.4012796612843944D+00, &
          0.1974406774312114D+00, &
          0.4012796612843944D+00, &
          0.4012796612843944D+00, &
          0.1974406774312114D+00, &
          0.4012796612843944D+00, &
          0.1089883669604236D+00, &
          0.7820232660791528D+00, &
          0.1089883669604236D+00, &
          0.1089883669604236D+00, &
          0.7820232660791528D+00, &
          0.1089883669604236D+00, &
          0.4472748897677901D+00, &
          0.1054502204644198D+00, &
          0.4472748897677901D+00, &
          0.4472748897677901D+00, &
          0.1054502204644198D+00, &
          0.4472748897677901D+00, &
          0.0692335448571578D+00, &
          0.8615329102856844D+00, &
          0.0692335448571578D+00, &
          0.0692335448571578D+00, &
          0.8615329102856844D+00, &
          0.0692335448571578D+00, &
          0.2475145234973103D+00, &
          0.5049709530053793D+00, &
          0.2475145234973103D+00, &
          0.2475145234973103D+00, &
          0.5049709530053793D+00, &
          0.2475145234973103D+00, &
          0.4533377451558121D+00, &
          0.0933245096883758D+00, &
          0.4533377451558121D+00, &
          0.4533377451558121D+00, &
          0.0933245096883758D+00, &
          0.4533377451558121D+00, &
          0.0180837411268107D+00, &
          0.9638325177463786D+00, &
          0.0180837411268107D+00, &
          0.0180837411268107D+00, &
          0.9638325177463786D+00, &
          0.0180837411268107D+00, &
          0.0069174965695053D+00, &
          0.1297897131077726D+00, &
          0.8632927903227221D+00, &
          0.1297897131077726D+00, &
          0.8632927903227221D+00, &
          0.0069174965695053D+00, &
          0.0310061742359910D+00, &
          0.3779418721990778D+00, &
          0.5910519535649312D+00, &
          0.3779418721990778D+00, &
          0.5910519535649312D+00, &
          0.0310061742359910D+00, &
          0.2203251202566757D+00, &
          0.0644631311955231D+00, &
          0.7152117485478011D+00, &
          0.0644631311955231D+00, &
          0.7152117485478011D+00, &
          0.2203251202566757D+00, &
          0.2824754920851798D+00, &
          0.1625899734623090D+00, &
          0.5549345344525113D+00, &
          0.1625899734623090D+00, &
          0.5549345344525113D+00, &
          0.2824754920851798D+00, &
          0.2824754920851798D+00, &
          0.1625899734623090D+00, &
          0.5549345344525113D+00, &
          0.1625899734623090D+00, &
          0.5549345344525113D+00, &
          0.2824754920851798D+00, &
          0.0547883511283449D+00, &
          0.0081184148094398D+00, &
          0.9370932340622153D+00, &
          0.0081184148094398D+00, &
          0.9370932340622153D+00, &
          0.0547883511283449D+00, &
          0.0547883511283449D+00, &
          0.0081184148094398D+00, &
          0.9370932340622153D+00, &
          0.0081184148094398D+00, &
          0.9370932340622153D+00, &
          0.0547883511283449D+00, &
          0.1178005364736766D+00, &
          0.6582273496659641D+00, &
          0.2239721138603593D+00, &
          0.6582273496659641D+00, &
          0.2239721138603593D+00, &
          0.1178005364736766D+00, &
          0.1178005364736766D+00, &
          0.6582273496659641D+00, &
          0.2239721138603593D+00, &
          0.6582273496659641D+00, &
          0.2239721138603593D+00, &
          0.1178005364736766D+00, &
          0.0355268541454468D+00, &
          0.2876874815978630D+00, &
          0.6767856642566902D+00, &
          0.2876874815978630D+00, &
          0.6767856642566902D+00, &
          0.0355268541454468D+00, &
          0.0355268541454468D+00, &
          0.2876874815978630D+00, &
          0.6767856642566902D+00, &
          0.2876874815978630D+00, &
          0.6767856642566902D+00, &
          0.0355268541454468D+00, &
          0.2073917255967352D+00, &
          0.2839484400308330D+00, &
          0.5086598343724318D+00, &
          0.2839484400308330D+00, &
          0.5086598343724318D+00, &
          0.2073917255967352D+00, &
          0.2073917255967352D+00, &
          0.2839484400308330D+00, &
          0.5086598343724318D+00, &
          0.2839484400308330D+00, &
          0.5086598343724318D+00, &
          0.2073917255967352D+00, &
          0.0862230026956853D+00, &
          0.0282457003834953D+00, &
          0.8855312969208194D+00, &
          0.0282457003834953D+00, &
          0.8855312969208194D+00, &
          0.0862230026956853D+00, &
          0.0862230026956853D+00, &
          0.0282457003834953D+00, &
          0.8855312969208194D+00, &
          0.0282457003834953D+00, &
          0.8855312969208194D+00, &
          0.0862230026956853D+00, &
          0.0056196516036861D+00, &
          0.3915043117238656D+00, &
          0.6028760366724483D+00, &
          0.3915043117238656D+00, &
          0.6028760366724483D+00, &
          0.0056196516036861D+00, &
          0.0056196516036861D+00, &
          0.3915043117238656D+00, &
          0.6028760366724483D+00, &
          0.3915043117238656D+00, &
          0.6028760366724483D+00, &
          0.0056196516036861D+00, &
          0.3389265718024346D+00, &
          0.0737667613724345D+00, &
          0.5873066668251310D+00, &
          0.0737667613724345D+00, &
          0.5873066668251310D+00, &
          0.3389265718024346D+00, &
          0.3389265718024346D+00, &
          0.0737667613724345D+00, &
          0.5873066668251310D+00, &
          0.0737667613724345D+00, &
          0.5873066668251310D+00, &
          0.3389265718024346D+00, &
          0.0080205188978798D+00, &
          0.1582419524500696D+00, &
          0.8337375286520506D+00, &
          0.1582419524500696D+00, &
          0.8337375286520506D+00, &
          0.0080205188978798D+00, &
          0.0080205188978798D+00, &
          0.1582419524500696D+00, &
          0.8337375286520506D+00, &
          0.1582419524500696D+00, &
          0.8337375286520506D+00, &
          0.0080205188978798D+00, &
          0.0761183063986999D+00, &
          0.2270465393157224D+00, &
          0.6968351542855776D+00, &
          0.2270465393157224D+00, &
          0.6968351542855776D+00, &
          0.0761183063986999D+00, &
          0.0761183063986999D+00, &
          0.2270465393157224D+00, &
          0.6968351542855776D+00, &
          0.2270465393157224D+00, &
          0.6968351542855776D+00, &
          0.0761183063986999D+00, &
          0.0674860196923538D+00, &
          0.0107402320272877D+00, &
          0.9217737482803585D+00, &
          0.0107402320272877D+00, &
          0.9217737482803585D+00, &
          0.0674860196923538D+00, &
          0.0674860196923538D+00, &
          0.0107402320272877D+00, &
          0.9217737482803585D+00, &
          0.0107402320272877D+00, &
          0.9217737482803585D+00, &
          0.0674860196923538D+00, &
          0.1495289831435912D+00, &
          0.0450629754431290D+00, &
          0.8054080414132798D+00, &
          0.0450629754431290D+00, &
          0.8054080414132798D+00, &
          0.1495289831435912D+00, &
          0.1495289831435912D+00, &
          0.0450629754431290D+00, &
          0.8054080414132798D+00, &
          0.0450629754431290D+00, &
          0.8054080414132798D+00, &
          0.1495289831435912D+00, &
          0.0536599901206231D+00, &
          0.3899263834124326D+00, &
          0.5564136264669443D+00, &
          0.3899263834124326D+00, &
          0.5564136264669443D+00, &
          0.0536599901206231D+00, &
          0.0536599901206231D+00, &
          0.3899263834124326D+00, &
          0.5564136264669443D+00, &
          0.3899263834124326D+00, &
          0.5564136264669443D+00, &
          0.0536599901206231D+00, &
          0.4039638109964547D+00, &
          0.0142616821078597D+00, &
          0.5817745068956857D+00, &
          0.0142616821078597D+00, &
          0.5817745068956857D+00, &
          0.4039638109964547D+00, &
          0.4039638109964547D+00, &
          0.0142616821078597D+00, &
          0.5817745068956857D+00, &
          0.0142616821078597D+00, &
          0.5817745068956857D+00, &
          0.4039638109964547D+00, &
          0.0088330831166114D+00, &
          0.3314062155351160D+00, &
          0.6597607013482726D+00, &
          0.3314062155351160D+00, &
          0.6597607013482726D+00, &
          0.0088330831166114D+00, &
          0.0088330831166114D+00, &
          0.3314062155351160D+00, &
          0.6597607013482726D+00, &
          0.3314062155351160D+00, &
          0.6597607013482726D+00, &
          0.0088330831166114D+00, &
          0.0103923793593637D+00, &
          0.1223644156202090D+00, &
          0.8672432050204273D+00, &
          0.1223644156202090D+00, &
          0.8672432050204273D+00, &
          0.0103923793593637D+00, &
          0.0103923793593637D+00, &
          0.1223644156202090D+00, &
          0.8672432050204273D+00, &
          0.1223644156202090D+00, &
          0.8672432050204273D+00, &
          0.0103923793593637D+00, &
          0.1324448044285512D+00, &
          0.2929652829983945D+00, &
          0.5745899125730543D+00, &
          0.2929652829983945D+00, &
          0.5745899125730543D+00, &
          0.1324448044285512D+00, &
          0.1324448044285512D+00, &
          0.2929652829983945D+00, &
          0.5745899125730543D+00, &
          0.2929652829983945D+00, &
          0.5745899125730543D+00, &
          0.1324448044285512D+00, &
          0.0535147091599154D+00, &
          0.2534748953911949D+00, &
          0.6930103954488898D+00, &
          0.2534748953911949D+00, &
          0.6930103954488898D+00, &
          0.0535147091599154D+00, &
          0.0535147091599154D+00, &
          0.2534748953911949D+00, &
          0.6930103954488898D+00, &
          0.2534748953911949D+00, &
          0.6930103954488898D+00, &
          0.0535147091599154D+00, &
          0.1461188761758480D+00, &
          0.0440785872085001D+00, &
          0.8098025366156519D+00, &
          0.0440785872085001D+00, &
          0.8098025366156519D+00, &
          0.1461188761758480D+00, &
          0.1461188761758480D+00, &
          0.0440785872085001D+00, &
          0.8098025366156519D+00, &
          0.0440785872085001D+00, &
          0.8098025366156519D+00, &
          0.1461188761758480D+00, &
          0.1264604222462178D+00, &
          0.5031457633505478D+00, &
          0.3703938144032344D+00, &
          0.5031457633505478D+00, &
          0.3703938144032344D+00, &
          0.1264604222462178D+00, &
          0.1264604222462178D+00, &
          0.5031457633505478D+00, &
          0.3703938144032344D+00, &
          0.5031457633505478D+00, &
          0.3703938144032344D+00, &
          0.1264604222462178D+00, &
          0.2433623618774898D+00, &
          0.0087499979841172D+00, &
          0.7478876401383930D+00, &
          0.0087499979841172D+00, &
          0.7478876401383930D+00, &
          0.2433623618774898D+00, &
          0.2433623618774898D+00, &
          0.0087499979841172D+00, &
          0.7478876401383930D+00, &
          0.0087499979841172D+00, &
          0.7478876401383930D+00, &
          0.2433623618774898D+00, &
          0.2569581578007861D+00, &
          0.4433279475253553D+00, &
          0.2997138946738586D+00, &
          0.4433279475253553D+00, &
          0.2997138946738586D+00, &
          0.2569581578007861D+00, &
          0.2569581578007861D+00, &
          0.4433279475253553D+00, &
          0.2997138946738586D+00, &
          0.4433279475253553D+00, &
          0.2997138946738586D+00, &
          0.2569581578007861D+00, &
          0.0112771540337380D+00, &
          0.2485412372956145D+00, &
          0.7401816086706475D+00, &
          0.2485412372956145D+00, &
          0.7401816086706475D+00, &
          0.0112771540337380D+00, &
          0.0112771540337380D+00, &
          0.2485412372956145D+00, &
          0.7401816086706475D+00, &
          0.2485412372956145D+00, &
          0.7401816086706475D+00, &
          0.0112771540337380D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.7622415410729753D+00, &
          0.2377584589270247D+00, &
          0.5372666694077471D+00, &
          0.4627333305922529D+00, &
          0.9783048646074941D+00, &
          0.0216951353925058D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9776297496930976D+00, &
          0.9776297496930976D+00, &
          0.9776297496930976D+00, &
          0.0223702503069024D+00, &
          0.0223702503069024D+00, &
          0.0223702503069024D+00, &
          0.8365993771201422D+00, &
          0.8365993771201422D+00, &
          0.8365993771201422D+00, &
          0.1634006228798578D+00, &
          0.1634006228798578D+00, &
          0.1634006228798578D+00, &
          0.3997141071880904D+00, &
          0.3997141071880904D+00, &
          0.3997141071880904D+00, &
          0.6002858928119097D+00, &
          0.6002858928119097D+00, &
          0.6002858928119097D+00, &
          0.6540149239138839D+00, &
          0.6540149239138839D+00, &
          0.6540149239138839D+00, &
          0.3459850760861161D+00, &
          0.3459850760861161D+00, &
          0.3459850760861161D+00, &
          0.6177200943808741D+00, &
          0.6177200943808741D+00, &
          0.6177200943808741D+00, &
          0.3822799056191259D+00, &
          0.3822799056191259D+00, &
          0.3822799056191259D+00, &
          0.7414020374054957D+00, &
          0.7414020374054957D+00, &
          0.7414020374054957D+00, &
          0.2585979625945043D+00, &
          0.2585979625945043D+00, &
          0.2585979625945043D+00, &
          0.3998937993818196D+00, &
          0.3998937993818196D+00, &
          0.3998937993818196D+00, &
          0.6001062006181804D+00, &
          0.6001062006181804D+00, &
          0.6001062006181804D+00, &
          0.8605804354127558D+00, &
          0.8605804354127558D+00, &
          0.8605804354127558D+00, &
          0.1394195645872443D+00, &
          0.1394195645872443D+00, &
          0.1394195645872443D+00, &
          0.8884925585076502D+00, &
          0.8884925585076502D+00, &
          0.8884925585076502D+00, &
          0.1115074414923499D+00, &
          0.1115074414923499D+00, &
          0.1115074414923499D+00, &
          0.1501005238707534D+00, &
          0.1501005238707534D+00, &
          0.1501005238707534D+00, &
          0.8498994761292467D+00, &
          0.8498994761292467D+00, &
          0.8498994761292467D+00, &
          0.8161162302018612D+00, &
          0.8161162302018612D+00, &
          0.8161162302018612D+00, &
          0.1838837697981388D+00, &
          0.1838837697981388D+00, &
          0.1838837697981388D+00, &
          0.9492984016563156D+00, &
          0.9492984016563156D+00, &
          0.9492984016563156D+00, &
          0.0507015983436844D+00, &
          0.0507015983436844D+00, &
          0.0507015983436844D+00, &
          0.8207059136960952D+00, &
          0.8207059136960952D+00, &
          0.8207059136960952D+00, &
          0.1792940863039048D+00, &
          0.1792940863039048D+00, &
          0.1792940863039048D+00, &
          0.5709032950041624D+00, &
          0.5709032950041624D+00, &
          0.5709032950041624D+00, &
          0.4290967049958376D+00, &
          0.4290967049958376D+00, &
          0.4290967049958376D+00, &
          0.9796258147828216D+00, &
          0.9796258147828216D+00, &
          0.9796258147828216D+00, &
          0.0203741852171783D+00, &
          0.0203741852171783D+00, &
          0.0203741852171783D+00, &
          0.9955700370275176D+00, &
          0.9955700370275176D+00, &
          0.9955700370275176D+00, &
          0.0044299629724824D+00, &
          0.0044299629724824D+00, &
          0.0044299629724824D+00, &
          0.9940911680624349D+00, &
          0.9940911680624349D+00, &
          0.9940911680624349D+00, &
          0.0059088319375651D+00, &
          0.0059088319375651D+00, &
          0.0059088319375651D+00, &
          0.9768533807402040D+00, &
          0.9768533807402040D+00, &
          0.9768533807402040D+00, &
          0.0231466192597961D+00, &
          0.0231466192597961D+00, &
          0.0231466192597961D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5784175741511786D+00, &
          0.5784175741511786D+00, &
          0.5784175741511786D+00, &
          0.5784175741511786D+00, &
          0.5784175741511786D+00, &
          0.5784175741511786D+00, &
          0.4215824258488214D+00, &
          0.4215824258488214D+00, &
          0.4215824258488214D+00, &
          0.4215824258488214D+00, &
          0.4215824258488214D+00, &
          0.4215824258488214D+00, &
          0.3016244371993055D+00, &
          0.3016244371993055D+00, &
          0.3016244371993055D+00, &
          0.3016244371993055D+00, &
          0.3016244371993055D+00, &
          0.3016244371993055D+00, &
          0.6983755628006945D+00, &
          0.6983755628006945D+00, &
          0.6983755628006945D+00, &
          0.6983755628006945D+00, &
          0.6983755628006945D+00, &
          0.6983755628006945D+00, &
          0.2557424176244907D+00, &
          0.2557424176244907D+00, &
          0.2557424176244907D+00, &
          0.2557424176244907D+00, &
          0.2557424176244907D+00, &
          0.2557424176244907D+00, &
          0.7442575823755093D+00, &
          0.7442575823755093D+00, &
          0.7442575823755093D+00, &
          0.7442575823755093D+00, &
          0.7442575823755093D+00, &
          0.7442575823755093D+00, &
          0.7933947373630322D+00, &
          0.7933947373630322D+00, &
          0.7933947373630322D+00, &
          0.7933947373630322D+00, &
          0.7933947373630322D+00, &
          0.7933947373630322D+00, &
          0.2066052626369678D+00, &
          0.2066052626369678D+00, &
          0.2066052626369678D+00, &
          0.2066052626369678D+00, &
          0.2066052626369678D+00, &
          0.2066052626369678D+00, &
          0.7415891076871005D+00, &
          0.7415891076871005D+00, &
          0.7415891076871005D+00, &
          0.7415891076871005D+00, &
          0.7415891076871005D+00, &
          0.7415891076871005D+00, &
          0.2584108923128995D+00, &
          0.2584108923128995D+00, &
          0.2584108923128995D+00, &
          0.2584108923128995D+00, &
          0.2584108923128995D+00, &
          0.2584108923128995D+00, &
          0.3990327902043654D+00, &
          0.3990327902043654D+00, &
          0.3990327902043654D+00, &
          0.3990327902043654D+00, &
          0.3990327902043654D+00, &
          0.3990327902043654D+00, &
          0.6009672097956347D+00, &
          0.6009672097956347D+00, &
          0.6009672097956347D+00, &
          0.6009672097956347D+00, &
          0.6009672097956347D+00, &
          0.6009672097956347D+00, &
          0.6658797023892886D+00, &
          0.6658797023892886D+00, &
          0.6658797023892886D+00, &
          0.6658797023892886D+00, &
          0.6658797023892886D+00, &
          0.6658797023892886D+00, &
          0.3341202976107114D+00, &
          0.3341202976107114D+00, &
          0.3341202976107114D+00, &
          0.3341202976107114D+00, &
          0.3341202976107114D+00, &
          0.3341202976107114D+00, &
          0.3477027664621174D+00, &
          0.3477027664621174D+00, &
          0.3477027664621174D+00, &
          0.3477027664621174D+00, &
          0.3477027664621174D+00, &
          0.3477027664621174D+00, &
          0.6522972335378826D+00, &
          0.6522972335378826D+00, &
          0.6522972335378826D+00, &
          0.6522972335378826D+00, &
          0.6522972335378826D+00, &
          0.6522972335378826D+00, &
          0.8121138484281050D+00, &
          0.8121138484281050D+00, &
          0.8121138484281050D+00, &
          0.8121138484281050D+00, &
          0.8121138484281050D+00, &
          0.8121138484281050D+00, &
          0.1878861515718951D+00, &
          0.1878861515718951D+00, &
          0.1878861515718951D+00, &
          0.1878861515718951D+00, &
          0.1878861515718951D+00, &
          0.1878861515718951D+00, &
          0.8889931649237686D+00, &
          0.8889931649237686D+00, &
          0.8889931649237686D+00, &
          0.8889931649237686D+00, &
          0.8889931649237686D+00, &
          0.8889931649237686D+00, &
          0.1110068350762314D+00, &
          0.1110068350762314D+00, &
          0.1110068350762314D+00, &
          0.1110068350762314D+00, &
          0.1110068350762314D+00, &
          0.1110068350762314D+00, &
          0.9152677185839518D+00, &
          0.9152677185839518D+00, &
          0.9152677185839518D+00, &
          0.9152677185839518D+00, &
          0.9152677185839518D+00, &
          0.9152677185839518D+00, &
          0.0847322814160482D+00, &
          0.0847322814160482D+00, &
          0.0847322814160482D+00, &
          0.0847322814160482D+00, &
          0.0847322814160482D+00, &
          0.0847322814160482D+00, &
          0.6968471564249525D+00, &
          0.6968471564249525D+00, &
          0.6968471564249525D+00, &
          0.6968471564249525D+00, &
          0.6968471564249525D+00, &
          0.6968471564249525D+00, &
          0.3031528435750474D+00, &
          0.3031528435750474D+00, &
          0.3031528435750474D+00, &
          0.3031528435750474D+00, &
          0.3031528435750474D+00, &
          0.3031528435750474D+00, &
          0.9181070186825648D+00, &
          0.9181070186825648D+00, &
          0.9181070186825648D+00, &
          0.9181070186825648D+00, &
          0.9181070186825648D+00, &
          0.9181070186825648D+00, &
          0.0818929813174352D+00, &
          0.0818929813174352D+00, &
          0.0818929813174352D+00, &
          0.0818929813174352D+00, &
          0.0818929813174352D+00, &
          0.0818929813174352D+00, &
          0.9772924850792593D+00, &
          0.9772924850792593D+00, &
          0.9772924850792593D+00, &
          0.9772924850792593D+00, &
          0.9772924850792593D+00, &
          0.9772924850792593D+00, &
          0.0227075149207407D+00, &
          0.0227075149207407D+00, &
          0.0227075149207407D+00, &
          0.0227075149207407D+00, &
          0.0227075149207407D+00, &
          0.0227075149207407D+00, &
          0.8548397574539421D+00, &
          0.8548397574539421D+00, &
          0.8548397574539421D+00, &
          0.8548397574539421D+00, &
          0.8548397574539421D+00, &
          0.8548397574539421D+00, &
          0.1451602425460579D+00, &
          0.1451602425460579D+00, &
          0.1451602425460579D+00, &
          0.1451602425460579D+00, &
          0.1451602425460579D+00, &
          0.1451602425460579D+00, &
          0.9950140347594980D+00, &
          0.9950140347594980D+00, &
          0.9950140347594980D+00, &
          0.9950140347594980D+00, &
          0.9950140347594980D+00, &
          0.9950140347594980D+00, &
          0.0049859652405020D+00, &
          0.0049859652405020D+00, &
          0.0049859652405020D+00, &
          0.0049859652405020D+00, &
          0.0049859652405020D+00, &
          0.0049859652405020D+00, &
          0.9589371765386895D+00, &
          0.9589371765386895D+00, &
          0.9589371765386895D+00, &
          0.9589371765386895D+00, &
          0.9589371765386895D+00, &
          0.9589371765386895D+00, &
          0.0410628234613105D+00, &
          0.0410628234613105D+00, &
          0.0410628234613105D+00, &
          0.0410628234613105D+00, &
          0.0410628234613105D+00, &
          0.0410628234613105D+00, &
          0.9922849190190410D+00, &
          0.9922849190190410D+00, &
          0.9922849190190410D+00, &
          0.9922849190190410D+00, &
          0.9922849190190410D+00, &
          0.9922849190190410D+00, &
          0.0077150809809589D+00, &
          0.0077150809809589D+00, &
          0.0077150809809589D+00, &
          0.0077150809809589D+00, &
          0.0077150809809589D+00, &
          0.0077150809809589D+00, &
          0.9284954503109000D+00, &
          0.9284954503109000D+00, &
          0.9284954503109000D+00, &
          0.9284954503109000D+00, &
          0.9284954503109000D+00, &
          0.9284954503109000D+00, &
          0.0715045496891001D+00, &
          0.0715045496891001D+00, &
          0.0715045496891001D+00, &
          0.0715045496891001D+00, &
          0.0715045496891001D+00, &
          0.0715045496891001D+00, &
          0.8243592634966279D+00, &
          0.8243592634966279D+00, &
          0.8243592634966279D+00, &
          0.8243592634966279D+00, &
          0.8243592634966279D+00, &
          0.8243592634966279D+00, &
          0.1756407365033721D+00, &
          0.1756407365033721D+00, &
          0.1756407365033721D+00, &
          0.1756407365033721D+00, &
          0.1756407365033721D+00, &
          0.1756407365033721D+00, &
          0.9514300035294792D+00, &
          0.9514300035294792D+00, &
          0.9514300035294792D+00, &
          0.9514300035294792D+00, &
          0.9514300035294792D+00, &
          0.9514300035294792D+00, &
          0.0485699964705208D+00, &
          0.0485699964705208D+00, &
          0.0485699964705208D+00, &
          0.0485699964705208D+00, &
          0.0485699964705208D+00, &
          0.0485699964705208D+00, &
          0.8768250142401974D+00, &
          0.8768250142401974D+00, &
          0.8768250142401974D+00, &
          0.8768250142401974D+00, &
          0.8768250142401974D+00, &
          0.8768250142401974D+00, &
          0.1231749857598025D+00, &
          0.1231749857598025D+00, &
          0.1231749857598025D+00, &
          0.1231749857598025D+00, &
          0.1231749857598025D+00, &
          0.1231749857598025D+00, &
          0.6209054808190849D+00, &
          0.6209054808190849D+00, &
          0.6209054808190849D+00, &
          0.6209054808190849D+00, &
          0.6209054808190849D+00, &
          0.6209054808190849D+00, &
          0.3790945191809150D+00, &
          0.3790945191809150D+00, &
          0.3790945191809150D+00, &
          0.3790945191809150D+00, &
          0.3790945191809150D+00, &
          0.3790945191809150D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          0.0067665109724298D+00, &
          0.0067665109724298D+00, &
          0.0026669976869157D+00, &
          0.0026669976869157D+00, &
          0.0021982966836122D+00, &
          0.0021982966836122D+00, &
          0.0018035576610240D+00, &
          0.0018035576610240D+00, &
          0.0018035576610240D+00, &
          0.0063544305083746D+00, &
          0.0063544305083746D+00, &
          0.0063544305083746D+00, &
          0.0016368347217208D+00, &
          0.0016368347217208D+00, &
          0.0016368347217208D+00, &
          0.0004893654007837D+00, &
          0.0004893654007837D+00, &
          0.0004893654007837D+00, &
          0.0018464185130403D+00, &
          0.0018464185130403D+00, &
          0.0018464185130403D+00, &
          0.0018464185130403D+00, &
          0.0018464185130403D+00, &
          0.0018464185130403D+00, &
          0.0002773402798746D+00, &
          0.0002773402798746D+00, &
          0.0002773402798746D+00, &
          0.0002773402798746D+00, &
          0.0002773402798746D+00, &
          0.0002773402798746D+00, &
          0.0028236344372154D+00, &
          0.0028236344372154D+00, &
          0.0028236344372154D+00, &
          0.0028236344372154D+00, &
          0.0028236344372154D+00, &
          0.0028236344372154D+00, &
          0.0077260759989331D+00, &
          0.0077260759989331D+00, &
          0.0077260759989331D+00, &
          0.0077260759989331D+00, &
          0.0077260759989331D+00, &
          0.0077260759989331D+00, &
          0.0002193983773241D+00, &
          0.0002193983773241D+00, &
          0.0002193983773241D+00, &
          0.0002193983773241D+00, &
          0.0002193983773241D+00, &
          0.0002193983773241D+00, &
          0.0043666503668884D+00, &
          0.0043666503668884D+00, &
          0.0043666503668884D+00, &
          0.0043666503668884D+00, &
          0.0043666503668884D+00, &
          0.0043666503668884D+00, &
          0.0040292939228197D+00, &
          0.0040292939228197D+00, &
          0.0040292939228197D+00, &
          0.0040292939228197D+00, &
          0.0040292939228197D+00, &
          0.0040292939228197D+00, &
          0.0011307686923958D+00, &
          0.0011307686923958D+00, &
          0.0011307686923958D+00, &
          0.0011307686923958D+00, &
          0.0011307686923958D+00, &
          0.0011307686923958D+00, &
          0.0042878636066862D+00, &
          0.0042878636066862D+00, &
          0.0042878636066862D+00, &
          0.0042878636066862D+00, &
          0.0042878636066862D+00, &
          0.0042878636066862D+00, &
          0.0015922560724046D+00, &
          0.0015922560724046D+00, &
          0.0015922560724046D+00, &
          0.0015922560724046D+00, &
          0.0015922560724046D+00, &
          0.0015922560724046D+00, &
          0.0020526689955830D+00, &
          0.0020526689955830D+00, &
          0.0020526689955830D+00, &
          0.0020526689955830D+00, &
          0.0020526689955830D+00, &
          0.0020526689955830D+00, &
          0.0038947953255298D+00, &
          0.0038947953255298D+00, &
          0.0038947953255298D+00, &
          0.0038947953255298D+00, &
          0.0038947953255298D+00, &
          0.0038947953255298D+00, &
          0.0018034020896591D+00, &
          0.0018034020896591D+00, &
          0.0018034020896591D+00, &
          0.0018034020896591D+00, &
          0.0018034020896591D+00, &
          0.0018034020896591D+00, &
          0.0049066738140960D+00, &
          0.0049066738140960D+00, &
          0.0049066738140960D+00, &
          0.0049066738140960D+00, &
          0.0049066738140960D+00, &
          0.0049066738140960D+00, &
          0.0008385245625735D+00, &
          0.0008385245625735D+00, &
          0.0008385245625735D+00, &
          0.0008385245625735D+00, &
          0.0008385245625735D+00, &
          0.0008385245625735D+00, &
          0.0013444227016777D+00, &
          0.0013444227016777D+00, &
          0.0013444227016777D+00, &
          0.0013444227016777D+00, &
          0.0013444227016777D+00, &
          0.0013444227016777D+00, &
          0.0012123512024300D+00, &
          0.0012123512024300D+00, &
          0.0012123512024300D+00, &
          0.0012123512024300D+00, &
          0.0012123512024300D+00, &
          0.0012123512024300D+00, &
          0.0002704967330299D+00, &
          0.0002704967330299D+00, &
          0.0002704967330299D+00, &
          0.0002704967330299D+00, &
          0.0002704967330299D+00, &
          0.0002704967330299D+00, &
          0.0011944957381629D+00, &
          0.0011944957381629D+00, &
          0.0011944957381629D+00, &
          0.0011944957381629D+00, &
          0.0011944957381629D+00, &
          0.0011944957381629D+00, &
          0.0031031422250087D+00, &
          0.0031031422250087D+00, &
          0.0031031422250087D+00, &
          0.0031031422250087D+00, &
          0.0031031422250087D+00, &
          0.0031031422250087D+00, &
          0.0051316606545523D+00, &
          0.0051316606545523D+00, &
          0.0051316606545523D+00, &
          0.0051316606545523D+00, &
          0.0051316606545523D+00, &
          0.0051316606545523D+00, &
          0.0048030743388154D+00, &
          0.0048030743388154D+00, &
          0.0048030743388154D+00, &
          0.0048030743388154D+00, &
          0.0048030743388154D+00, &
          0.0048030743388154D+00, &
          0.0048030743388154D+00, &
          0.0048030743388154D+00, &
          0.0048030743388154D+00, &
          0.0048030743388154D+00, &
          0.0048030743388154D+00, &
          0.0048030743388154D+00, &
          0.0008150376242850D+00, &
          0.0008150376242850D+00, &
          0.0008150376242850D+00, &
          0.0008150376242850D+00, &
          0.0008150376242850D+00, &
          0.0008150376242850D+00, &
          0.0008150376242850D+00, &
          0.0008150376242850D+00, &
          0.0008150376242850D+00, &
          0.0008150376242850D+00, &
          0.0008150376242850D+00, &
          0.0008150376242850D+00, &
          0.0044984898722334D+00, &
          0.0044984898722334D+00, &
          0.0044984898722334D+00, &
          0.0044984898722334D+00, &
          0.0044984898722334D+00, &
          0.0044984898722334D+00, &
          0.0044984898722334D+00, &
          0.0044984898722334D+00, &
          0.0044984898722334D+00, &
          0.0044984898722334D+00, &
          0.0044984898722334D+00, &
          0.0044984898722334D+00, &
          0.0027317773513874D+00, &
          0.0027317773513874D+00, &
          0.0027317773513874D+00, &
          0.0027317773513874D+00, &
          0.0027317773513874D+00, &
          0.0027317773513874D+00, &
          0.0027317773513874D+00, &
          0.0027317773513874D+00, &
          0.0027317773513874D+00, &
          0.0027317773513874D+00, &
          0.0027317773513874D+00, &
          0.0027317773513874D+00, &
          0.0039179954921652D+00, &
          0.0039179954921652D+00, &
          0.0039179954921652D+00, &
          0.0039179954921652D+00, &
          0.0039179954921652D+00, &
          0.0039179954921652D+00, &
          0.0039179954921652D+00, &
          0.0039179954921652D+00, &
          0.0039179954921652D+00, &
          0.0039179954921652D+00, &
          0.0039179954921652D+00, &
          0.0039179954921652D+00, &
          0.0008946420147168D+00, &
          0.0008946420147168D+00, &
          0.0008946420147168D+00, &
          0.0008946420147168D+00, &
          0.0008946420147168D+00, &
          0.0008946420147168D+00, &
          0.0008946420147168D+00, &
          0.0008946420147168D+00, &
          0.0008946420147168D+00, &
          0.0008946420147168D+00, &
          0.0008946420147168D+00, &
          0.0008946420147168D+00, &
          0.0012726195819812D+00, &
          0.0012726195819812D+00, &
          0.0012726195819812D+00, &
          0.0012726195819812D+00, &
          0.0012726195819812D+00, &
          0.0012726195819812D+00, &
          0.0012726195819812D+00, &
          0.0012726195819812D+00, &
          0.0012726195819812D+00, &
          0.0012726195819812D+00, &
          0.0012726195819812D+00, &
          0.0012726195819812D+00, &
          0.0045391175168163D+00, &
          0.0045391175168163D+00, &
          0.0045391175168163D+00, &
          0.0045391175168163D+00, &
          0.0045391175168163D+00, &
          0.0045391175168163D+00, &
          0.0045391175168163D+00, &
          0.0045391175168163D+00, &
          0.0045391175168163D+00, &
          0.0045391175168163D+00, &
          0.0045391175168163D+00, &
          0.0045391175168163D+00, &
          0.0012280771399416D+00, &
          0.0012280771399416D+00, &
          0.0012280771399416D+00, &
          0.0012280771399416D+00, &
          0.0012280771399416D+00, &
          0.0012280771399416D+00, &
          0.0012280771399416D+00, &
          0.0012280771399416D+00, &
          0.0012280771399416D+00, &
          0.0012280771399416D+00, &
          0.0012280771399416D+00, &
          0.0012280771399416D+00, &
          0.0024813526864083D+00, &
          0.0024813526864083D+00, &
          0.0024813526864083D+00, &
          0.0024813526864083D+00, &
          0.0024813526864083D+00, &
          0.0024813526864083D+00, &
          0.0024813526864083D+00, &
          0.0024813526864083D+00, &
          0.0024813526864083D+00, &
          0.0024813526864083D+00, &
          0.0024813526864083D+00, &
          0.0024813526864083D+00, &
          0.0006742435916066D+00, &
          0.0006742435916066D+00, &
          0.0006742435916066D+00, &
          0.0006742435916066D+00, &
          0.0006742435916066D+00, &
          0.0006742435916066D+00, &
          0.0006742435916066D+00, &
          0.0006742435916066D+00, &
          0.0006742435916066D+00, &
          0.0006742435916066D+00, &
          0.0006742435916066D+00, &
          0.0006742435916066D+00, &
          0.0030856887684142D+00, &
          0.0030856887684142D+00, &
          0.0030856887684142D+00, &
          0.0030856887684142D+00, &
          0.0030856887684142D+00, &
          0.0030856887684142D+00, &
          0.0030856887684142D+00, &
          0.0030856887684142D+00, &
          0.0030856887684142D+00, &
          0.0030856887684142D+00, &
          0.0030856887684142D+00, &
          0.0030856887684142D+00, &
          0.0028277448828289D+00, &
          0.0028277448828289D+00, &
          0.0028277448828289D+00, &
          0.0028277448828289D+00, &
          0.0028277448828289D+00, &
          0.0028277448828289D+00, &
          0.0028277448828289D+00, &
          0.0028277448828289D+00, &
          0.0028277448828289D+00, &
          0.0028277448828289D+00, &
          0.0028277448828289D+00, &
          0.0028277448828289D+00, &
          0.0007880156766982D+00, &
          0.0007880156766982D+00, &
          0.0007880156766982D+00, &
          0.0007880156766982D+00, &
          0.0007880156766982D+00, &
          0.0007880156766982D+00, &
          0.0007880156766982D+00, &
          0.0007880156766982D+00, &
          0.0007880156766982D+00, &
          0.0007880156766982D+00, &
          0.0007880156766982D+00, &
          0.0007880156766982D+00, &
          0.0009734356783438D+00, &
          0.0009734356783438D+00, &
          0.0009734356783438D+00, &
          0.0009734356783438D+00, &
          0.0009734356783438D+00, &
          0.0009734356783438D+00, &
          0.0009734356783438D+00, &
          0.0009734356783438D+00, &
          0.0009734356783438D+00, &
          0.0009734356783438D+00, &
          0.0009734356783438D+00, &
          0.0009734356783438D+00, &
          0.0002344473397349D+00, &
          0.0002344473397349D+00, &
          0.0002344473397349D+00, &
          0.0002344473397349D+00, &
          0.0002344473397349D+00, &
          0.0002344473397349D+00, &
          0.0002344473397349D+00, &
          0.0002344473397349D+00, &
          0.0002344473397349D+00, &
          0.0002344473397349D+00, &
          0.0002344473397349D+00, &
          0.0002344473397349D+00, &
          0.0027549281349627D+00, &
          0.0027549281349627D+00, &
          0.0027549281349627D+00, &
          0.0027549281349627D+00, &
          0.0027549281349627D+00, &
          0.0027549281349627D+00, &
          0.0027549281349627D+00, &
          0.0027549281349627D+00, &
          0.0027549281349627D+00, &
          0.0027549281349627D+00, &
          0.0027549281349627D+00, &
          0.0027549281349627D+00, &
          0.0008530955496580D+00, &
          0.0008530955496580D+00, &
          0.0008530955496580D+00, &
          0.0008530955496580D+00, &
          0.0008530955496580D+00, &
          0.0008530955496580D+00, &
          0.0008530955496580D+00, &
          0.0008530955496580D+00, &
          0.0008530955496580D+00, &
          0.0008530955496580D+00, &
          0.0008530955496580D+00, &
          0.0008530955496580D+00, &
          0.0017207910666049D+00, &
          0.0017207910666049D+00, &
          0.0017207910666049D+00, &
          0.0017207910666049D+00, &
          0.0017207910666049D+00, &
          0.0017207910666049D+00, &
          0.0017207910666049D+00, &
          0.0017207910666049D+00, &
          0.0017207910666049D+00, &
          0.0017207910666049D+00, &
          0.0017207910666049D+00, &
          0.0017207910666049D+00, &
          0.0052380601608220D+00, &
          0.0052380601608220D+00, &
          0.0052380601608220D+00, &
          0.0052380601608220D+00, &
          0.0052380601608220D+00, &
          0.0052380601608220D+00, &
          0.0052380601608220D+00, &
          0.0052380601608220D+00, &
          0.0052380601608220D+00, &
          0.0052380601608220D+00, &
          0.0052380601608220D+00, &
          0.0052380601608220D+00, &
          0.0007780035531901D+00, &
          0.0007780035531901D+00, &
          0.0007780035531901D+00, &
          0.0007780035531901D+00, &
          0.0007780035531901D+00, &
          0.0007780035531901D+00, &
          0.0007780035531901D+00, &
          0.0007780035531901D+00, &
          0.0007780035531901D+00, &
          0.0007780035531901D+00, &
          0.0007780035531901D+00, &
          0.0007780035531901D+00, &
          0.0027997199341340D+00, &
          0.0027997199341340D+00, &
          0.0027997199341340D+00, &
          0.0027997199341340D+00, &
          0.0027997199341340D+00, &
          0.0027997199341340D+00, &
          0.0027997199341340D+00, &
          0.0027997199341340D+00, &
          0.0027997199341340D+00, &
          0.0027997199341340D+00, &
          0.0027997199341340D+00, &
          0.0027997199341340D+00, &
          0.0018871269258399D+00, &
          0.0018871269258399D+00, &
          0.0018871269258399D+00, &
          0.0018871269258399D+00, &
          0.0018871269258399D+00, &
          0.0018871269258399D+00, &
          0.0018871269258399D+00, &
          0.0018871269258399D+00, &
          0.0018871269258399D+00, &
          0.0018871269258399D+00, &
          0.0018871269258399D+00, &
          0.0018871269258399D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule20 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule20() returns the prism rule of precision 20.
!
!  Discussion:
!
!    We suppose we are given a triangular prism P with vertices
!    (1,0,0), (0,1,0), (0,0,0),
!    (1,0,1), (0,1,1), (0,0,1),
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over P is 
!    approximated by Q as follows:
!
!    Q = volume(P) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 May 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order symmetric cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 122, Number 1, pages 148-171, 24 August 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real x(n), y(n), z(n): the coordinates of quadrature points.
!
!    real w(n): the quadrature weights, which sum to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: n_save = 518

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.0144110873698983D+00, &
          0.0144110873698983D+00, &
          0.9711778252602034D+00, &
          0.0594974134225245D+00, &
          0.0594974134225245D+00, &
          0.8810051731549510D+00, &
          0.4716249531173324D+00, &
          0.4716249531173324D+00, &
          0.0567500937653351D+00, &
          0.4716249531173324D+00, &
          0.4716249531173324D+00, &
          0.0567500937653351D+00, &
          0.4084460043805324D+00, &
          0.4084460043805324D+00, &
          0.1831079912389352D+00, &
          0.4084460043805324D+00, &
          0.4084460043805324D+00, &
          0.1831079912389352D+00, &
          0.2544280410411358D+00, &
          0.2544280410411358D+00, &
          0.4911439179177284D+00, &
          0.2544280410411358D+00, &
          0.2544280410411358D+00, &
          0.4911439179177284D+00, &
          0.4774877095488475D+00, &
          0.4774877095488475D+00, &
          0.0450245809023050D+00, &
          0.4774877095488475D+00, &
          0.4774877095488475D+00, &
          0.0450245809023050D+00, &
          0.0080292636682187D+00, &
          0.0080292636682187D+00, &
          0.9839414726635627D+00, &
          0.0080292636682187D+00, &
          0.0080292636682187D+00, &
          0.9839414726635627D+00, &
          0.4968613511179721D+00, &
          0.4968613511179721D+00, &
          0.0062772977640557D+00, &
          0.4968613511179721D+00, &
          0.4968613511179721D+00, &
          0.0062772977640557D+00, &
          0.4952615769053836D+00, &
          0.4952615769053836D+00, &
          0.0094768461892327D+00, &
          0.4952615769053836D+00, &
          0.4952615769053836D+00, &
          0.0094768461892327D+00, &
          0.0353578097171312D+00, &
          0.0353578097171312D+00, &
          0.9292843805657376D+00, &
          0.0353578097171312D+00, &
          0.0353578097171312D+00, &
          0.9292843805657376D+00, &
          0.2045970282680664D+00, &
          0.2045970282680664D+00, &
          0.5908059434638671D+00, &
          0.2045970282680664D+00, &
          0.2045970282680664D+00, &
          0.5908059434638671D+00, &
          0.1406131665045265D+00, &
          0.1406131665045265D+00, &
          0.7187736669909469D+00, &
          0.1406131665045265D+00, &
          0.1406131665045265D+00, &
          0.7187736669909469D+00, &
          0.0982561053377771D+00, &
          0.0982561053377771D+00, &
          0.8034877893244458D+00, &
          0.0982561053377771D+00, &
          0.0982561053377771D+00, &
          0.8034877893244458D+00, &
          0.0505951664588140D+00, &
          0.0505951664588140D+00, &
          0.8988096670823720D+00, &
          0.0505951664588140D+00, &
          0.0505951664588140D+00, &
          0.8988096670823720D+00, &
          0.4494558038329708D+00, &
          0.4494558038329708D+00, &
          0.1010883923340583D+00, &
          0.4494558038329708D+00, &
          0.4494558038329708D+00, &
          0.1010883923340583D+00, &
          0.2188296334498393D+00, &
          0.2188296334498393D+00, &
          0.5623407331003214D+00, &
          0.2188296334498393D+00, &
          0.2188296334498393D+00, &
          0.5623407331003214D+00, &
          0.3839280551265914D+00, &
          0.3839280551265914D+00, &
          0.2321438897468171D+00, &
          0.3839280551265914D+00, &
          0.3839280551265914D+00, &
          0.2321438897468171D+00, &
          0.0599706591707595D+00, &
          0.0599706591707595D+00, &
          0.8800586816584810D+00, &
          0.0599706591707595D+00, &
          0.0599706591707595D+00, &
          0.8800586816584810D+00, &
          0.2735494152095392D+00, &
          0.2735494152095392D+00, &
          0.4529011695809215D+00, &
          0.2735494152095392D+00, &
          0.2735494152095392D+00, &
          0.4529011695809215D+00, &
          0.4907393281824849D+00, &
          0.4907393281824849D+00, &
          0.0185213436350301D+00, &
          0.4907393281824849D+00, &
          0.4907393281824849D+00, &
          0.0185213436350301D+00, &
          0.0114930108302544D+00, &
          0.0114930108302544D+00, &
          0.9770139783394911D+00, &
          0.0114930108302544D+00, &
          0.0114930108302544D+00, &
          0.9770139783394911D+00, &
          0.4048678509173179D+00, &
          0.4048678509173179D+00, &
          0.1902642981653643D+00, &
          0.4048678509173179D+00, &
          0.4048678509173179D+00, &
          0.1902642981653643D+00, &
          0.2917111538437177D+00, &
          0.2917111538437177D+00, &
          0.4165776923125646D+00, &
          0.2917111538437177D+00, &
          0.2917111538437177D+00, &
          0.4165776923125646D+00, &
          0.4459340787080981D+00, &
          0.4459340787080981D+00, &
          0.1081318425838039D+00, &
          0.4459340787080981D+00, &
          0.4459340787080981D+00, &
          0.1081318425838039D+00, &
          0.2107382124681030D+00, &
          0.0313712403793080D+00, &
          0.2107382124681030D+00, &
          0.7578905471525890D+00, &
          0.0313712403793080D+00, &
          0.7578905471525890D+00, &
          0.2031438036613264D+00, &
          0.1114400380412596D+00, &
          0.2031438036613264D+00, &
          0.6854161582974140D+00, &
          0.1114400380412596D+00, &
          0.6854161582974140D+00, &
          0.3734168270486775D+00, &
          0.0362901207122587D+00, &
          0.3734168270486775D+00, &
          0.5902930522390638D+00, &
          0.0362901207122587D+00, &
          0.5902930522390638D+00, &
          0.5331515253263366D+00, &
          0.1335223252058013D+00, &
          0.5331515253263366D+00, &
          0.3333261494678621D+00, &
          0.1335223252058013D+00, &
          0.3333261494678621D+00, &
          0.5331515253263366D+00, &
          0.1335223252058013D+00, &
          0.5331515253263366D+00, &
          0.3333261494678621D+00, &
          0.1335223252058013D+00, &
          0.3333261494678621D+00, &
          0.0031052202964266D+00, &
          0.0650661770452139D+00, &
          0.0031052202964266D+00, &
          0.9318286026583595D+00, &
          0.0650661770452139D+00, &
          0.9318286026583595D+00, &
          0.0031052202964266D+00, &
          0.0650661770452139D+00, &
          0.0031052202964266D+00, &
          0.9318286026583595D+00, &
          0.0650661770452139D+00, &
          0.9318286026583595D+00, &
          0.7387825606686489D+00, &
          0.0932063677657537D+00, &
          0.7387825606686489D+00, &
          0.1680110715655975D+00, &
          0.0932063677657537D+00, &
          0.1680110715655975D+00, &
          0.7387825606686489D+00, &
          0.0932063677657537D+00, &
          0.7387825606686489D+00, &
          0.1680110715655975D+00, &
          0.0932063677657537D+00, &
          0.1680110715655975D+00, &
          0.2402209181076414D+00, &
          0.0342816904321169D+00, &
          0.2402209181076414D+00, &
          0.7254973914602417D+00, &
          0.0342816904321169D+00, &
          0.7254973914602417D+00, &
          0.2402209181076414D+00, &
          0.0342816904321169D+00, &
          0.2402209181076414D+00, &
          0.7254973914602417D+00, &
          0.0342816904321169D+00, &
          0.7254973914602417D+00, &
          0.0816129838618181D+00, &
          0.3087274419444517D+00, &
          0.0816129838618181D+00, &
          0.6096595741937303D+00, &
          0.3087274419444517D+00, &
          0.6096595741937303D+00, &
          0.0816129838618181D+00, &
          0.3087274419444517D+00, &
          0.0816129838618181D+00, &
          0.6096595741937303D+00, &
          0.3087274419444517D+00, &
          0.6096595741937303D+00, &
          0.1728159430622389D+00, &
          0.3228866908185781D+00, &
          0.1728159430622389D+00, &
          0.5042973661191831D+00, &
          0.3228866908185781D+00, &
          0.5042973661191831D+00, &
          0.1728159430622389D+00, &
          0.3228866908185781D+00, &
          0.1728159430622389D+00, &
          0.5042973661191831D+00, &
          0.3228866908185781D+00, &
          0.5042973661191831D+00, &
          0.2247127800382396D+00, &
          0.1327427633967232D+00, &
          0.2247127800382396D+00, &
          0.6425444565650372D+00, &
          0.1327427633967232D+00, &
          0.6425444565650372D+00, &
          0.2247127800382396D+00, &
          0.1327427633967232D+00, &
          0.2247127800382396D+00, &
          0.6425444565650372D+00, &
          0.1327427633967232D+00, &
          0.6425444565650372D+00, &
          0.0285010351208627D+00, &
          0.1058138309132070D+00, &
          0.0285010351208627D+00, &
          0.8656851339659303D+00, &
          0.1058138309132070D+00, &
          0.8656851339659303D+00, &
          0.0285010351208627D+00, &
          0.1058138309132070D+00, &
          0.0285010351208627D+00, &
          0.8656851339659303D+00, &
          0.1058138309132070D+00, &
          0.8656851339659303D+00, &
          0.3896100790636586D+00, &
          0.0109712296494740D+00, &
          0.3896100790636586D+00, &
          0.5994186912868675D+00, &
          0.0109712296494740D+00, &
          0.5994186912868675D+00, &
          0.3896100790636586D+00, &
          0.0109712296494740D+00, &
          0.3896100790636586D+00, &
          0.5994186912868675D+00, &
          0.0109712296494740D+00, &
          0.5994186912868675D+00, &
          0.2493335646225524D+00, &
          0.3397868885786173D+00, &
          0.2493335646225524D+00, &
          0.4108795467988303D+00, &
          0.3397868885786173D+00, &
          0.4108795467988303D+00, &
          0.2493335646225524D+00, &
          0.3397868885786173D+00, &
          0.2493335646225524D+00, &
          0.4108795467988303D+00, &
          0.3397868885786173D+00, &
          0.4108795467988303D+00, &
          0.0070420236000393D+00, &
          0.3177207979690921D+00, &
          0.0070420236000393D+00, &
          0.6752371784308686D+00, &
          0.3177207979690921D+00, &
          0.6752371784308686D+00, &
          0.0070420236000393D+00, &
          0.3177207979690921D+00, &
          0.0070420236000393D+00, &
          0.6752371784308686D+00, &
          0.3177207979690921D+00, &
          0.6752371784308686D+00, &
          0.1229469656661441D+00, &
          0.0156675786078851D+00, &
          0.1229469656661441D+00, &
          0.8613854557259707D+00, &
          0.0156675786078851D+00, &
          0.8613854557259707D+00, &
          0.1229469656661441D+00, &
          0.0156675786078851D+00, &
          0.1229469656661441D+00, &
          0.8613854557259707D+00, &
          0.0156675786078851D+00, &
          0.8613854557259707D+00, &
          0.2259160262165568D+00, &
          0.0521518846441873D+00, &
          0.2259160262165568D+00, &
          0.7219320891392559D+00, &
          0.0521518846441873D+00, &
          0.7219320891392559D+00, &
          0.2259160262165568D+00, &
          0.0521518846441873D+00, &
          0.2259160262165568D+00, &
          0.7219320891392559D+00, &
          0.0521518846441873D+00, &
          0.7219320891392559D+00, &
          0.0094714933266370D+00, &
          0.0496441415264223D+00, &
          0.0094714933266370D+00, &
          0.9408843651469406D+00, &
          0.0496441415264223D+00, &
          0.9408843651469406D+00, &
          0.0094714933266370D+00, &
          0.0496441415264223D+00, &
          0.0094714933266370D+00, &
          0.9408843651469406D+00, &
          0.0496441415264223D+00, &
          0.9408843651469406D+00, &
          0.0524311424127234D+00, &
          0.1452721772236918D+00, &
          0.0524311424127234D+00, &
          0.8022966803635848D+00, &
          0.1452721772236918D+00, &
          0.8022966803635848D+00, &
          0.0524311424127234D+00, &
          0.1452721772236918D+00, &
          0.0524311424127234D+00, &
          0.8022966803635848D+00, &
          0.1452721772236918D+00, &
          0.8022966803635848D+00, &
          0.0088722025555106D+00, &
          0.3816302433543118D+00, &
          0.0088722025555106D+00, &
          0.6094975540901776D+00, &
          0.3816302433543118D+00, &
          0.6094975540901776D+00, &
          0.0088722025555106D+00, &
          0.3816302433543118D+00, &
          0.0088722025555106D+00, &
          0.6094975540901776D+00, &
          0.3816302433543118D+00, &
          0.6094975540901776D+00, &
          0.2592518129005248D+00, &
          0.5551397100853696D+00, &
          0.2592518129005248D+00, &
          0.1856084770141057D+00, &
          0.5551397100853696D+00, &
          0.1856084770141057D+00, &
          0.2592518129005248D+00, &
          0.5551397100853696D+00, &
          0.2592518129005248D+00, &
          0.1856084770141057D+00, &
          0.5551397100853696D+00, &
          0.1856084770141057D+00, &
          0.3527097683072426D+00, &
          0.0454609714340920D+00, &
          0.3527097683072426D+00, &
          0.6018292602586655D+00, &
          0.0454609714340920D+00, &
          0.6018292602586655D+00, &
          0.3527097683072426D+00, &
          0.0454609714340920D+00, &
          0.3527097683072426D+00, &
          0.6018292602586655D+00, &
          0.0454609714340920D+00, &
          0.6018292602586655D+00, &
          0.1450709676542589D+00, &
          0.0118170254385161D+00, &
          0.1450709676542589D+00, &
          0.8431120069072251D+00, &
          0.0118170254385161D+00, &
          0.8431120069072251D+00, &
          0.1450709676542589D+00, &
          0.0118170254385161D+00, &
          0.1450709676542589D+00, &
          0.8431120069072251D+00, &
          0.0118170254385161D+00, &
          0.8431120069072251D+00, &
          0.1874143309040665D+00, &
          0.1371981166173774D+00, &
          0.1874143309040665D+00, &
          0.6753875524785560D+00, &
          0.1371981166173774D+00, &
          0.6753875524785560D+00, &
          0.1874143309040665D+00, &
          0.1371981166173774D+00, &
          0.1874143309040665D+00, &
          0.6753875524785560D+00, &
          0.1371981166173774D+00, &
          0.6753875524785560D+00, &
          0.0681401608659605D+00, &
          0.1498311852045099D+00, &
          0.0681401608659605D+00, &
          0.7820286539295297D+00, &
          0.1498311852045099D+00, &
          0.7820286539295297D+00, &
          0.0681401608659605D+00, &
          0.1498311852045099D+00, &
          0.0681401608659605D+00, &
          0.7820286539295297D+00, &
          0.1498311852045099D+00, &
          0.7820286539295297D+00, &
          0.3340585216441213D+00, &
          0.0549696286134690D+00, &
          0.3340585216441213D+00, &
          0.6109718497424097D+00, &
          0.0549696286134690D+00, &
          0.6109718497424097D+00, &
          0.3340585216441213D+00, &
          0.0549696286134690D+00, &
          0.3340585216441213D+00, &
          0.6109718497424097D+00, &
          0.0549696286134690D+00, &
          0.6109718497424097D+00, &
          0.0678515797670052D+00, &
          0.1205147270278029D+00, &
          0.0678515797670052D+00, &
          0.8116336932051920D+00, &
          0.1205147270278029D+00, &
          0.8116336932051920D+00, &
          0.0678515797670052D+00, &
          0.1205147270278029D+00, &
          0.0678515797670052D+00, &
          0.8116336932051920D+00, &
          0.1205147270278029D+00, &
          0.8116336932051920D+00, &
          0.5023136055537912D+00, &
          0.1076311573233394D+00, &
          0.5023136055537912D+00, &
          0.3900552371228693D+00, &
          0.1076311573233394D+00, &
          0.3900552371228693D+00, &
          0.5023136055537912D+00, &
          0.1076311573233394D+00, &
          0.5023136055537912D+00, &
          0.3900552371228693D+00, &
          0.1076311573233394D+00, &
          0.3900552371228693D+00, &
          0.0078193997434242D+00, &
          0.2556542053755758D+00, &
          0.0078193997434242D+00, &
          0.7365263948810000D+00, &
          0.2556542053755758D+00, &
          0.7365263948810000D+00, &
          0.0078193997434242D+00, &
          0.2556542053755758D+00, &
          0.0078193997434242D+00, &
          0.7365263948810000D+00, &
          0.2556542053755758D+00, &
          0.7365263948810000D+00, &
          0.6609230621965876D+00, &
          0.0983956958507364D+00, &
          0.6609230621965876D+00, &
          0.2406812419526760D+00, &
          0.0983956958507364D+00, &
          0.2406812419526760D+00, &
          0.6609230621965876D+00, &
          0.0983956958507364D+00, &
          0.6609230621965876D+00, &
          0.2406812419526760D+00, &
          0.0983956958507364D+00, &
          0.2406812419526760D+00, &
          0.1784736485903538D+00, &
          0.0017965399943124D+00, &
          0.1784736485903538D+00, &
          0.8197298114153339D+00, &
          0.0017965399943124D+00, &
          0.8197298114153339D+00, &
          0.1784736485903538D+00, &
          0.0017965399943124D+00, &
          0.1784736485903538D+00, &
          0.8197298114153339D+00, &
          0.0017965399943124D+00, &
          0.8197298114153339D+00, &
          0.2532227359651285D+00, &
          0.0131217590858908D+00, &
          0.2532227359651285D+00, &
          0.7336555049489807D+00, &
          0.0131217590858908D+00, &
          0.7336555049489807D+00, &
          0.2532227359651285D+00, &
          0.0131217590858908D+00, &
          0.2532227359651285D+00, &
          0.7336555049489807D+00, &
          0.0131217590858908D+00, &
          0.7336555049489807D+00, &
          0.0144350536713311D+00, &
          0.0686738759058901D+00, &
          0.0144350536713311D+00, &
          0.9168910704227788D+00, &
          0.0686738759058901D+00, &
          0.9168910704227788D+00, &
          0.0144350536713311D+00, &
          0.0686738759058901D+00, &
          0.0144350536713311D+00, &
          0.9168910704227788D+00, &
          0.0686738759058901D+00, &
          0.9168910704227788D+00, &
          0.2794904874495157D+00, &
          0.1274003485132245D+00, &
          0.2794904874495157D+00, &
          0.5931091640372598D+00, &
          0.1274003485132245D+00, &
          0.5931091640372598D+00, &
          0.2794904874495157D+00, &
          0.1274003485132245D+00, &
          0.2794904874495157D+00, &
          0.5931091640372598D+00, &
          0.1274003485132245D+00, &
          0.5931091640372598D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.3333333333333333D+00, &
          0.3333333333333333D+00, &
          0.0144110873698983D+00, &
          0.9711778252602034D+00, &
          0.0144110873698983D+00, &
          0.0594974134225245D+00, &
          0.8810051731549510D+00, &
          0.0594974134225245D+00, &
          0.4716249531173324D+00, &
          0.0567500937653351D+00, &
          0.4716249531173324D+00, &
          0.4716249531173324D+00, &
          0.0567500937653351D+00, &
          0.4716249531173324D+00, &
          0.4084460043805324D+00, &
          0.1831079912389352D+00, &
          0.4084460043805324D+00, &
          0.4084460043805324D+00, &
          0.1831079912389352D+00, &
          0.4084460043805324D+00, &
          0.2544280410411358D+00, &
          0.4911439179177284D+00, &
          0.2544280410411358D+00, &
          0.2544280410411358D+00, &
          0.4911439179177284D+00, &
          0.2544280410411358D+00, &
          0.4774877095488475D+00, &
          0.0450245809023050D+00, &
          0.4774877095488475D+00, &
          0.4774877095488475D+00, &
          0.0450245809023050D+00, &
          0.4774877095488475D+00, &
          0.0080292636682187D+00, &
          0.9839414726635627D+00, &
          0.0080292636682187D+00, &
          0.0080292636682187D+00, &
          0.9839414726635627D+00, &
          0.0080292636682187D+00, &
          0.4968613511179721D+00, &
          0.0062772977640557D+00, &
          0.4968613511179721D+00, &
          0.4968613511179721D+00, &
          0.0062772977640557D+00, &
          0.4968613511179721D+00, &
          0.4952615769053836D+00, &
          0.0094768461892327D+00, &
          0.4952615769053836D+00, &
          0.4952615769053836D+00, &
          0.0094768461892327D+00, &
          0.4952615769053836D+00, &
          0.0353578097171312D+00, &
          0.9292843805657376D+00, &
          0.0353578097171312D+00, &
          0.0353578097171312D+00, &
          0.9292843805657376D+00, &
          0.0353578097171312D+00, &
          0.2045970282680664D+00, &
          0.5908059434638671D+00, &
          0.2045970282680664D+00, &
          0.2045970282680664D+00, &
          0.5908059434638671D+00, &
          0.2045970282680664D+00, &
          0.1406131665045265D+00, &
          0.7187736669909469D+00, &
          0.1406131665045265D+00, &
          0.1406131665045265D+00, &
          0.7187736669909469D+00, &
          0.1406131665045265D+00, &
          0.0982561053377771D+00, &
          0.8034877893244458D+00, &
          0.0982561053377771D+00, &
          0.0982561053377771D+00, &
          0.8034877893244458D+00, &
          0.0982561053377771D+00, &
          0.0505951664588140D+00, &
          0.8988096670823720D+00, &
          0.0505951664588140D+00, &
          0.0505951664588140D+00, &
          0.8988096670823720D+00, &
          0.0505951664588140D+00, &
          0.4494558038329708D+00, &
          0.1010883923340583D+00, &
          0.4494558038329708D+00, &
          0.4494558038329708D+00, &
          0.1010883923340583D+00, &
          0.4494558038329708D+00, &
          0.2188296334498393D+00, &
          0.5623407331003214D+00, &
          0.2188296334498393D+00, &
          0.2188296334498393D+00, &
          0.5623407331003214D+00, &
          0.2188296334498393D+00, &
          0.3839280551265914D+00, &
          0.2321438897468171D+00, &
          0.3839280551265914D+00, &
          0.3839280551265914D+00, &
          0.2321438897468171D+00, &
          0.3839280551265914D+00, &
          0.0599706591707595D+00, &
          0.8800586816584810D+00, &
          0.0599706591707595D+00, &
          0.0599706591707595D+00, &
          0.8800586816584810D+00, &
          0.0599706591707595D+00, &
          0.2735494152095392D+00, &
          0.4529011695809215D+00, &
          0.2735494152095392D+00, &
          0.2735494152095392D+00, &
          0.4529011695809215D+00, &
          0.2735494152095392D+00, &
          0.4907393281824849D+00, &
          0.0185213436350301D+00, &
          0.4907393281824849D+00, &
          0.4907393281824849D+00, &
          0.0185213436350301D+00, &
          0.4907393281824849D+00, &
          0.0114930108302544D+00, &
          0.9770139783394911D+00, &
          0.0114930108302544D+00, &
          0.0114930108302544D+00, &
          0.9770139783394911D+00, &
          0.0114930108302544D+00, &
          0.4048678509173179D+00, &
          0.1902642981653643D+00, &
          0.4048678509173179D+00, &
          0.4048678509173179D+00, &
          0.1902642981653643D+00, &
          0.4048678509173179D+00, &
          0.2917111538437177D+00, &
          0.4165776923125646D+00, &
          0.2917111538437177D+00, &
          0.2917111538437177D+00, &
          0.4165776923125646D+00, &
          0.2917111538437177D+00, &
          0.4459340787080981D+00, &
          0.1081318425838039D+00, &
          0.4459340787080981D+00, &
          0.4459340787080981D+00, &
          0.1081318425838039D+00, &
          0.4459340787080981D+00, &
          0.0313712403793080D+00, &
          0.2107382124681030D+00, &
          0.7578905471525890D+00, &
          0.2107382124681030D+00, &
          0.7578905471525890D+00, &
          0.0313712403793080D+00, &
          0.1114400380412596D+00, &
          0.2031438036613264D+00, &
          0.6854161582974140D+00, &
          0.2031438036613264D+00, &
          0.6854161582974140D+00, &
          0.1114400380412596D+00, &
          0.0362901207122587D+00, &
          0.3734168270486775D+00, &
          0.5902930522390638D+00, &
          0.3734168270486775D+00, &
          0.5902930522390638D+00, &
          0.0362901207122587D+00, &
          0.1335223252058013D+00, &
          0.5331515253263366D+00, &
          0.3333261494678621D+00, &
          0.5331515253263366D+00, &
          0.3333261494678621D+00, &
          0.1335223252058013D+00, &
          0.1335223252058013D+00, &
          0.5331515253263366D+00, &
          0.3333261494678621D+00, &
          0.5331515253263366D+00, &
          0.3333261494678621D+00, &
          0.1335223252058013D+00, &
          0.0650661770452139D+00, &
          0.0031052202964266D+00, &
          0.9318286026583595D+00, &
          0.0031052202964266D+00, &
          0.9318286026583595D+00, &
          0.0650661770452139D+00, &
          0.0650661770452139D+00, &
          0.0031052202964266D+00, &
          0.9318286026583595D+00, &
          0.0031052202964266D+00, &
          0.9318286026583595D+00, &
          0.0650661770452139D+00, &
          0.0932063677657537D+00, &
          0.7387825606686489D+00, &
          0.1680110715655975D+00, &
          0.7387825606686489D+00, &
          0.1680110715655975D+00, &
          0.0932063677657537D+00, &
          0.0932063677657537D+00, &
          0.7387825606686489D+00, &
          0.1680110715655975D+00, &
          0.7387825606686489D+00, &
          0.1680110715655975D+00, &
          0.0932063677657537D+00, &
          0.0342816904321169D+00, &
          0.2402209181076414D+00, &
          0.7254973914602417D+00, &
          0.2402209181076414D+00, &
          0.7254973914602417D+00, &
          0.0342816904321169D+00, &
          0.0342816904321169D+00, &
          0.2402209181076414D+00, &
          0.7254973914602417D+00, &
          0.2402209181076414D+00, &
          0.7254973914602417D+00, &
          0.0342816904321169D+00, &
          0.3087274419444517D+00, &
          0.0816129838618181D+00, &
          0.6096595741937303D+00, &
          0.0816129838618181D+00, &
          0.6096595741937303D+00, &
          0.3087274419444517D+00, &
          0.3087274419444517D+00, &
          0.0816129838618181D+00, &
          0.6096595741937303D+00, &
          0.0816129838618181D+00, &
          0.6096595741937303D+00, &
          0.3087274419444517D+00, &
          0.3228866908185781D+00, &
          0.1728159430622389D+00, &
          0.5042973661191831D+00, &
          0.1728159430622389D+00, &
          0.5042973661191831D+00, &
          0.3228866908185781D+00, &
          0.3228866908185781D+00, &
          0.1728159430622389D+00, &
          0.5042973661191831D+00, &
          0.1728159430622389D+00, &
          0.5042973661191831D+00, &
          0.3228866908185781D+00, &
          0.1327427633967232D+00, &
          0.2247127800382396D+00, &
          0.6425444565650372D+00, &
          0.2247127800382396D+00, &
          0.6425444565650372D+00, &
          0.1327427633967232D+00, &
          0.1327427633967232D+00, &
          0.2247127800382396D+00, &
          0.6425444565650372D+00, &
          0.2247127800382396D+00, &
          0.6425444565650372D+00, &
          0.1327427633967232D+00, &
          0.1058138309132070D+00, &
          0.0285010351208627D+00, &
          0.8656851339659303D+00, &
          0.0285010351208627D+00, &
          0.8656851339659303D+00, &
          0.1058138309132070D+00, &
          0.1058138309132070D+00, &
          0.0285010351208627D+00, &
          0.8656851339659303D+00, &
          0.0285010351208627D+00, &
          0.8656851339659303D+00, &
          0.1058138309132070D+00, &
          0.0109712296494740D+00, &
          0.3896100790636586D+00, &
          0.5994186912868675D+00, &
          0.3896100790636586D+00, &
          0.5994186912868675D+00, &
          0.0109712296494740D+00, &
          0.0109712296494740D+00, &
          0.3896100790636586D+00, &
          0.5994186912868675D+00, &
          0.3896100790636586D+00, &
          0.5994186912868675D+00, &
          0.0109712296494740D+00, &
          0.3397868885786173D+00, &
          0.2493335646225524D+00, &
          0.4108795467988303D+00, &
          0.2493335646225524D+00, &
          0.4108795467988303D+00, &
          0.3397868885786173D+00, &
          0.3397868885786173D+00, &
          0.2493335646225524D+00, &
          0.4108795467988303D+00, &
          0.2493335646225524D+00, &
          0.4108795467988303D+00, &
          0.3397868885786173D+00, &
          0.3177207979690921D+00, &
          0.0070420236000393D+00, &
          0.6752371784308686D+00, &
          0.0070420236000393D+00, &
          0.6752371784308686D+00, &
          0.3177207979690921D+00, &
          0.3177207979690921D+00, &
          0.0070420236000393D+00, &
          0.6752371784308686D+00, &
          0.0070420236000393D+00, &
          0.6752371784308686D+00, &
          0.3177207979690921D+00, &
          0.0156675786078851D+00, &
          0.1229469656661441D+00, &
          0.8613854557259707D+00, &
          0.1229469656661441D+00, &
          0.8613854557259707D+00, &
          0.0156675786078851D+00, &
          0.0156675786078851D+00, &
          0.1229469656661441D+00, &
          0.8613854557259707D+00, &
          0.1229469656661441D+00, &
          0.8613854557259707D+00, &
          0.0156675786078851D+00, &
          0.0521518846441873D+00, &
          0.2259160262165568D+00, &
          0.7219320891392559D+00, &
          0.2259160262165568D+00, &
          0.7219320891392559D+00, &
          0.0521518846441873D+00, &
          0.0521518846441873D+00, &
          0.2259160262165568D+00, &
          0.7219320891392559D+00, &
          0.2259160262165568D+00, &
          0.7219320891392559D+00, &
          0.0521518846441873D+00, &
          0.0496441415264223D+00, &
          0.0094714933266370D+00, &
          0.9408843651469406D+00, &
          0.0094714933266370D+00, &
          0.9408843651469406D+00, &
          0.0496441415264223D+00, &
          0.0496441415264223D+00, &
          0.0094714933266370D+00, &
          0.9408843651469406D+00, &
          0.0094714933266370D+00, &
          0.9408843651469406D+00, &
          0.0496441415264223D+00, &
          0.1452721772236918D+00, &
          0.0524311424127234D+00, &
          0.8022966803635848D+00, &
          0.0524311424127234D+00, &
          0.8022966803635848D+00, &
          0.1452721772236918D+00, &
          0.1452721772236918D+00, &
          0.0524311424127234D+00, &
          0.8022966803635848D+00, &
          0.0524311424127234D+00, &
          0.8022966803635848D+00, &
          0.1452721772236918D+00, &
          0.3816302433543118D+00, &
          0.0088722025555106D+00, &
          0.6094975540901776D+00, &
          0.0088722025555106D+00, &
          0.6094975540901776D+00, &
          0.3816302433543118D+00, &
          0.3816302433543118D+00, &
          0.0088722025555106D+00, &
          0.6094975540901776D+00, &
          0.0088722025555106D+00, &
          0.6094975540901776D+00, &
          0.3816302433543118D+00, &
          0.5551397100853696D+00, &
          0.2592518129005248D+00, &
          0.1856084770141057D+00, &
          0.2592518129005248D+00, &
          0.1856084770141057D+00, &
          0.5551397100853696D+00, &
          0.5551397100853696D+00, &
          0.2592518129005248D+00, &
          0.1856084770141057D+00, &
          0.2592518129005248D+00, &
          0.1856084770141057D+00, &
          0.5551397100853696D+00, &
          0.0454609714340920D+00, &
          0.3527097683072426D+00, &
          0.6018292602586655D+00, &
          0.3527097683072426D+00, &
          0.6018292602586655D+00, &
          0.0454609714340920D+00, &
          0.0454609714340920D+00, &
          0.3527097683072426D+00, &
          0.6018292602586655D+00, &
          0.3527097683072426D+00, &
          0.6018292602586655D+00, &
          0.0454609714340920D+00, &
          0.0118170254385161D+00, &
          0.1450709676542589D+00, &
          0.8431120069072251D+00, &
          0.1450709676542589D+00, &
          0.8431120069072251D+00, &
          0.0118170254385161D+00, &
          0.0118170254385161D+00, &
          0.1450709676542589D+00, &
          0.8431120069072251D+00, &
          0.1450709676542589D+00, &
          0.8431120069072251D+00, &
          0.0118170254385161D+00, &
          0.1371981166173774D+00, &
          0.1874143309040665D+00, &
          0.6753875524785560D+00, &
          0.1874143309040665D+00, &
          0.6753875524785560D+00, &
          0.1371981166173774D+00, &
          0.1371981166173774D+00, &
          0.1874143309040665D+00, &
          0.6753875524785560D+00, &
          0.1874143309040665D+00, &
          0.6753875524785560D+00, &
          0.1371981166173774D+00, &
          0.1498311852045099D+00, &
          0.0681401608659605D+00, &
          0.7820286539295297D+00, &
          0.0681401608659605D+00, &
          0.7820286539295297D+00, &
          0.1498311852045099D+00, &
          0.1498311852045099D+00, &
          0.0681401608659605D+00, &
          0.7820286539295297D+00, &
          0.0681401608659605D+00, &
          0.7820286539295297D+00, &
          0.1498311852045099D+00, &
          0.0549696286134690D+00, &
          0.3340585216441213D+00, &
          0.6109718497424097D+00, &
          0.3340585216441213D+00, &
          0.6109718497424097D+00, &
          0.0549696286134690D+00, &
          0.0549696286134690D+00, &
          0.3340585216441213D+00, &
          0.6109718497424097D+00, &
          0.3340585216441213D+00, &
          0.6109718497424097D+00, &
          0.0549696286134690D+00, &
          0.1205147270278029D+00, &
          0.0678515797670052D+00, &
          0.8116336932051920D+00, &
          0.0678515797670052D+00, &
          0.8116336932051920D+00, &
          0.1205147270278029D+00, &
          0.1205147270278029D+00, &
          0.0678515797670052D+00, &
          0.8116336932051920D+00, &
          0.0678515797670052D+00, &
          0.8116336932051920D+00, &
          0.1205147270278029D+00, &
          0.1076311573233394D+00, &
          0.5023136055537912D+00, &
          0.3900552371228693D+00, &
          0.5023136055537912D+00, &
          0.3900552371228693D+00, &
          0.1076311573233394D+00, &
          0.1076311573233394D+00, &
          0.5023136055537912D+00, &
          0.3900552371228693D+00, &
          0.5023136055537912D+00, &
          0.3900552371228693D+00, &
          0.1076311573233394D+00, &
          0.2556542053755758D+00, &
          0.0078193997434242D+00, &
          0.7365263948810000D+00, &
          0.0078193997434242D+00, &
          0.7365263948810000D+00, &
          0.2556542053755758D+00, &
          0.2556542053755758D+00, &
          0.0078193997434242D+00, &
          0.7365263948810000D+00, &
          0.0078193997434242D+00, &
          0.7365263948810000D+00, &
          0.2556542053755758D+00, &
          0.0983956958507364D+00, &
          0.6609230621965876D+00, &
          0.2406812419526760D+00, &
          0.6609230621965876D+00, &
          0.2406812419526760D+00, &
          0.0983956958507364D+00, &
          0.0983956958507364D+00, &
          0.6609230621965876D+00, &
          0.2406812419526760D+00, &
          0.6609230621965876D+00, &
          0.2406812419526760D+00, &
          0.0983956958507364D+00, &
          0.0017965399943124D+00, &
          0.1784736485903538D+00, &
          0.8197298114153339D+00, &
          0.1784736485903538D+00, &
          0.8197298114153339D+00, &
          0.0017965399943124D+00, &
          0.0017965399943124D+00, &
          0.1784736485903538D+00, &
          0.8197298114153339D+00, &
          0.1784736485903538D+00, &
          0.8197298114153339D+00, &
          0.0017965399943124D+00, &
          0.0131217590858908D+00, &
          0.2532227359651285D+00, &
          0.7336555049489807D+00, &
          0.2532227359651285D+00, &
          0.7336555049489807D+00, &
          0.0131217590858908D+00, &
          0.0131217590858908D+00, &
          0.2532227359651285D+00, &
          0.7336555049489807D+00, &
          0.2532227359651285D+00, &
          0.7336555049489807D+00, &
          0.0131217590858908D+00, &
          0.0686738759058901D+00, &
          0.0144350536713311D+00, &
          0.9168910704227788D+00, &
          0.0144350536713311D+00, &
          0.9168910704227788D+00, &
          0.0686738759058901D+00, &
          0.0686738759058901D+00, &
          0.0144350536713311D+00, &
          0.9168910704227788D+00, &
          0.0144350536713311D+00, &
          0.9168910704227788D+00, &
          0.0686738759058901D+00, &
          0.1274003485132245D+00, &
          0.2794904874495157D+00, &
          0.5931091640372598D+00, &
          0.2794904874495157D+00, &
          0.5931091640372598D+00, &
          0.1274003485132245D+00, &
          0.1274003485132245D+00, &
          0.2794904874495157D+00, &
          0.5931091640372598D+00, &
          0.2794904874495157D+00, &
          0.5931091640372598D+00, &
          0.1274003485132245D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.9210083004731364D+00, &
          0.0789916995268636D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9076441430436997D+00, &
          0.9076441430436997D+00, &
          0.9076441430436997D+00, &
          0.0923558569563003D+00, &
          0.0923558569563003D+00, &
          0.0923558569563003D+00, &
          0.7654479960865263D+00, &
          0.7654479960865263D+00, &
          0.7654479960865263D+00, &
          0.2345520039134737D+00, &
          0.2345520039134737D+00, &
          0.2345520039134737D+00, &
          0.8572582404180138D+00, &
          0.8572582404180138D+00, &
          0.8572582404180138D+00, &
          0.1427417595819861D+00, &
          0.1427417595819861D+00, &
          0.1427417595819861D+00, &
          0.6651011543921541D+00, &
          0.6651011543921541D+00, &
          0.6651011543921541D+00, &
          0.3348988456078459D+00, &
          0.3348988456078459D+00, &
          0.3348988456078459D+00, &
          0.7626983072263107D+00, &
          0.7626983072263107D+00, &
          0.7626983072263107D+00, &
          0.2373016927736893D+00, &
          0.2373016927736893D+00, &
          0.2373016927736893D+00, &
          0.5839411665251317D+00, &
          0.5839411665251317D+00, &
          0.5839411665251317D+00, &
          0.4160588334748682D+00, &
          0.4160588334748682D+00, &
          0.4160588334748682D+00, &
          0.8529900626158916D+00, &
          0.8529900626158916D+00, &
          0.8529900626158916D+00, &
          0.1470099373841084D+00, &
          0.1470099373841084D+00, &
          0.1470099373841084D+00, &
          0.6690062435950173D+00, &
          0.6690062435950173D+00, &
          0.6690062435950173D+00, &
          0.3309937564049827D+00, &
          0.3309937564049827D+00, &
          0.3309937564049827D+00, &
          0.4508779258616127D+00, &
          0.4508779258616127D+00, &
          0.4508779258616127D+00, &
          0.5491220741383873D+00, &
          0.5491220741383873D+00, &
          0.5491220741383873D+00, &
          0.8420489345379816D+00, &
          0.8420489345379816D+00, &
          0.8420489345379816D+00, &
          0.1579510654620185D+00, &
          0.1579510654620185D+00, &
          0.1579510654620185D+00, &
          0.3603208005838872D+00, &
          0.3603208005838872D+00, &
          0.3603208005838872D+00, &
          0.6396791994161128D+00, &
          0.6396791994161128D+00, &
          0.6396791994161128D+00, &
          0.7720261679832757D+00, &
          0.7720261679832757D+00, &
          0.7720261679832757D+00, &
          0.2279738320167243D+00, &
          0.2279738320167243D+00, &
          0.2279738320167243D+00, &
          0.9830894614249917D+00, &
          0.9830894614249917D+00, &
          0.9830894614249917D+00, &
          0.0169105385750084D+00, &
          0.0169105385750084D+00, &
          0.0169105385750084D+00, &
          0.7580515195220092D+00, &
          0.7580515195220092D+00, &
          0.7580515195220092D+00, &
          0.2419484804779908D+00, &
          0.2419484804779908D+00, &
          0.2419484804779908D+00, &
          0.8207076808620348D+00, &
          0.8207076808620348D+00, &
          0.8207076808620348D+00, &
          0.1792923191379651D+00, &
          0.1792923191379651D+00, &
          0.1792923191379651D+00, &
          0.9471959220764619D+00, &
          0.9471959220764619D+00, &
          0.9471959220764619D+00, &
          0.0528040779235381D+00, &
          0.0528040779235381D+00, &
          0.0528040779235381D+00, &
          0.9896592480532946D+00, &
          0.9896592480532946D+00, &
          0.9896592480532946D+00, &
          0.0103407519467054D+00, &
          0.0103407519467054D+00, &
          0.0103407519467054D+00, &
          0.9943125669152897D+00, &
          0.9943125669152897D+00, &
          0.9943125669152897D+00, &
          0.0056874330847103D+00, &
          0.0056874330847103D+00, &
          0.0056874330847103D+00, &
          0.9633861983838395D+00, &
          0.9633861983838395D+00, &
          0.9633861983838395D+00, &
          0.0366138016161605D+00, &
          0.0366138016161605D+00, &
          0.0366138016161605D+00, &
          0.9582038855682118D+00, &
          0.9582038855682118D+00, &
          0.9582038855682118D+00, &
          0.0417961144317883D+00, &
          0.0417961144317883D+00, &
          0.0417961144317883D+00, &
          0.3139113946603033D+00, &
          0.3139113946603033D+00, &
          0.3139113946603033D+00, &
          0.6860886053396967D+00, &
          0.6860886053396967D+00, &
          0.6860886053396967D+00, &
          0.5189148057278750D+00, &
          0.5189148057278750D+00, &
          0.5189148057278750D+00, &
          0.4810851942721250D+00, &
          0.4810851942721250D+00, &
          0.4810851942721250D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.8901251913253538D+00, &
          0.8901251913253538D+00, &
          0.8901251913253538D+00, &
          0.8901251913253538D+00, &
          0.8901251913253538D+00, &
          0.8901251913253538D+00, &
          0.1098748086746462D+00, &
          0.1098748086746462D+00, &
          0.1098748086746462D+00, &
          0.1098748086746462D+00, &
          0.1098748086746462D+00, &
          0.1098748086746462D+00, &
          0.3557614873681827D+00, &
          0.3557614873681827D+00, &
          0.3557614873681827D+00, &
          0.3557614873681827D+00, &
          0.3557614873681827D+00, &
          0.3557614873681827D+00, &
          0.6442385126318173D+00, &
          0.6442385126318173D+00, &
          0.6442385126318173D+00, &
          0.6442385126318173D+00, &
          0.6442385126318173D+00, &
          0.6442385126318173D+00, &
          0.3778572463523129D+00, &
          0.3778572463523129D+00, &
          0.3778572463523129D+00, &
          0.3778572463523129D+00, &
          0.3778572463523129D+00, &
          0.3778572463523129D+00, &
          0.6221427536476871D+00, &
          0.6221427536476871D+00, &
          0.6221427536476871D+00, &
          0.6221427536476871D+00, &
          0.6221427536476871D+00, &
          0.6221427536476871D+00, &
          0.7130345491338336D+00, &
          0.7130345491338336D+00, &
          0.7130345491338336D+00, &
          0.7130345491338336D+00, &
          0.7130345491338336D+00, &
          0.7130345491338336D+00, &
          0.2869654508661664D+00, &
          0.2869654508661664D+00, &
          0.2869654508661664D+00, &
          0.2869654508661664D+00, &
          0.2869654508661664D+00, &
          0.2869654508661664D+00, &
          0.6210154934017852D+00, &
          0.6210154934017852D+00, &
          0.6210154934017852D+00, &
          0.6210154934017852D+00, &
          0.6210154934017852D+00, &
          0.6210154934017852D+00, &
          0.3789845065982148D+00, &
          0.3789845065982148D+00, &
          0.3789845065982148D+00, &
          0.3789845065982148D+00, &
          0.3789845065982148D+00, &
          0.3789845065982148D+00, &
          0.3811077243266478D+00, &
          0.3811077243266478D+00, &
          0.3811077243266478D+00, &
          0.3811077243266478D+00, &
          0.3811077243266478D+00, &
          0.3811077243266478D+00, &
          0.6188922756733521D+00, &
          0.6188922756733521D+00, &
          0.6188922756733521D+00, &
          0.6188922756733521D+00, &
          0.6188922756733521D+00, &
          0.6188922756733521D+00, &
          0.7341832387945393D+00, &
          0.7341832387945393D+00, &
          0.7341832387945393D+00, &
          0.7341832387945393D+00, &
          0.7341832387945393D+00, &
          0.7341832387945393D+00, &
          0.2658167612054607D+00, &
          0.2658167612054607D+00, &
          0.2658167612054607D+00, &
          0.2658167612054607D+00, &
          0.2658167612054607D+00, &
          0.2658167612054607D+00, &
          0.4405032526890183D+00, &
          0.4405032526890183D+00, &
          0.4405032526890183D+00, &
          0.4405032526890183D+00, &
          0.4405032526890183D+00, &
          0.4405032526890183D+00, &
          0.5594967473109818D+00, &
          0.5594967473109818D+00, &
          0.5594967473109818D+00, &
          0.5594967473109818D+00, &
          0.5594967473109818D+00, &
          0.5594967473109818D+00, &
          0.7406281263771474D+00, &
          0.7406281263771474D+00, &
          0.7406281263771474D+00, &
          0.7406281263771474D+00, &
          0.7406281263771474D+00, &
          0.7406281263771474D+00, &
          0.2593718736228526D+00, &
          0.2593718736228526D+00, &
          0.2593718736228526D+00, &
          0.2593718736228526D+00, &
          0.2593718736228526D+00, &
          0.2593718736228526D+00, &
          0.4825288515428756D+00, &
          0.4825288515428756D+00, &
          0.4825288515428756D+00, &
          0.4825288515428756D+00, &
          0.4825288515428756D+00, &
          0.4825288515428756D+00, &
          0.5174711484571244D+00, &
          0.5174711484571244D+00, &
          0.5174711484571244D+00, &
          0.5174711484571244D+00, &
          0.5174711484571244D+00, &
          0.5174711484571244D+00, &
          0.6039123747541559D+00, &
          0.6039123747541559D+00, &
          0.6039123747541559D+00, &
          0.6039123747541559D+00, &
          0.6039123747541559D+00, &
          0.6039123747541559D+00, &
          0.3960876252458441D+00, &
          0.3960876252458441D+00, &
          0.3960876252458441D+00, &
          0.3960876252458441D+00, &
          0.3960876252458441D+00, &
          0.3960876252458441D+00, &
          0.8047375675549896D+00, &
          0.8047375675549896D+00, &
          0.8047375675549896D+00, &
          0.8047375675549896D+00, &
          0.8047375675549896D+00, &
          0.8047375675549896D+00, &
          0.1952624324450104D+00, &
          0.1952624324450104D+00, &
          0.1952624324450104D+00, &
          0.1952624324450104D+00, &
          0.1952624324450104D+00, &
          0.1952624324450104D+00, &
          0.9180415438599450D+00, &
          0.9180415438599450D+00, &
          0.9180415438599450D+00, &
          0.9180415438599450D+00, &
          0.9180415438599450D+00, &
          0.9180415438599450D+00, &
          0.0819584561400550D+00, &
          0.0819584561400550D+00, &
          0.0819584561400550D+00, &
          0.0819584561400550D+00, &
          0.0819584561400550D+00, &
          0.0819584561400550D+00, &
          0.8806920867999144D+00, &
          0.8806920867999144D+00, &
          0.8806920867999144D+00, &
          0.8806920867999144D+00, &
          0.8806920867999144D+00, &
          0.8806920867999144D+00, &
          0.1193079132000855D+00, &
          0.1193079132000855D+00, &
          0.1193079132000855D+00, &
          0.1193079132000855D+00, &
          0.1193079132000855D+00, &
          0.1193079132000855D+00, &
          0.7498396447851496D+00, &
          0.7498396447851496D+00, &
          0.7498396447851496D+00, &
          0.7498396447851496D+00, &
          0.7498396447851496D+00, &
          0.7498396447851496D+00, &
          0.2501603552148504D+00, &
          0.2501603552148504D+00, &
          0.2501603552148504D+00, &
          0.2501603552148504D+00, &
          0.2501603552148504D+00, &
          0.2501603552148504D+00, &
          0.9438500554714535D+00, &
          0.9438500554714535D+00, &
          0.9438500554714535D+00, &
          0.9438500554714535D+00, &
          0.9438500554714535D+00, &
          0.9438500554714535D+00, &
          0.0561499445285465D+00, &
          0.0561499445285465D+00, &
          0.0561499445285465D+00, &
          0.0561499445285465D+00, &
          0.0561499445285465D+00, &
          0.0561499445285465D+00, &
          0.9395243139937932D+00, &
          0.9395243139937932D+00, &
          0.9395243139937932D+00, &
          0.9395243139937932D+00, &
          0.9395243139937932D+00, &
          0.9395243139937932D+00, &
          0.0604756860062068D+00, &
          0.0604756860062068D+00, &
          0.0604756860062068D+00, &
          0.0604756860062068D+00, &
          0.0604756860062068D+00, &
          0.0604756860062068D+00, &
          0.8184995763288090D+00, &
          0.8184995763288090D+00, &
          0.8184995763288090D+00, &
          0.8184995763288090D+00, &
          0.8184995763288090D+00, &
          0.8184995763288090D+00, &
          0.1815004236711910D+00, &
          0.1815004236711910D+00, &
          0.1815004236711910D+00, &
          0.1815004236711910D+00, &
          0.1815004236711910D+00, &
          0.1815004236711910D+00, &
          0.9476062382889998D+00, &
          0.9476062382889998D+00, &
          0.9476062382889998D+00, &
          0.9476062382889998D+00, &
          0.9476062382889998D+00, &
          0.9476062382889998D+00, &
          0.0523937617110002D+00, &
          0.0523937617110002D+00, &
          0.0523937617110002D+00, &
          0.0523937617110002D+00, &
          0.0523937617110002D+00, &
          0.0523937617110002D+00, &
          0.9567123523504457D+00, &
          0.9567123523504457D+00, &
          0.9567123523504457D+00, &
          0.9567123523504457D+00, &
          0.9567123523504457D+00, &
          0.9567123523504457D+00, &
          0.0432876476495543D+00, &
          0.0432876476495543D+00, &
          0.0432876476495543D+00, &
          0.0432876476495543D+00, &
          0.0432876476495543D+00, &
          0.0432876476495543D+00, &
          0.9875892078343605D+00, &
          0.9875892078343605D+00, &
          0.9875892078343605D+00, &
          0.9875892078343605D+00, &
          0.9875892078343605D+00, &
          0.9875892078343605D+00, &
          0.0124107921656395D+00, &
          0.0124107921656395D+00, &
          0.0124107921656395D+00, &
          0.0124107921656395D+00, &
          0.0124107921656395D+00, &
          0.0124107921656395D+00, &
          0.9678315615098414D+00, &
          0.9678315615098414D+00, &
          0.9678315615098414D+00, &
          0.9678315615098414D+00, &
          0.9678315615098414D+00, &
          0.9678315615098414D+00, &
          0.0321684384901587D+00, &
          0.0321684384901587D+00, &
          0.0321684384901587D+00, &
          0.0321684384901587D+00, &
          0.0321684384901587D+00, &
          0.0321684384901587D+00, &
          0.8827873971246739D+00, &
          0.8827873971246739D+00, &
          0.8827873971246739D+00, &
          0.8827873971246739D+00, &
          0.8827873971246739D+00, &
          0.8827873971246739D+00, &
          0.1172126028753261D+00, &
          0.1172126028753261D+00, &
          0.1172126028753261D+00, &
          0.1172126028753261D+00, &
          0.1172126028753261D+00, &
          0.1172126028753261D+00, &
          0.7591086801871259D+00, &
          0.7591086801871259D+00, &
          0.7591086801871259D+00, &
          0.7591086801871259D+00, &
          0.7591086801871259D+00, &
          0.7591086801871259D+00, &
          0.2408913198128742D+00, &
          0.2408913198128742D+00, &
          0.2408913198128742D+00, &
          0.2408913198128742D+00, &
          0.2408913198128742D+00, &
          0.2408913198128742D+00, &
          0.8583410425577678D+00, &
          0.8583410425577678D+00, &
          0.8583410425577678D+00, &
          0.8583410425577678D+00, &
          0.8583410425577678D+00, &
          0.8583410425577678D+00, &
          0.1416589574422322D+00, &
          0.1416589574422322D+00, &
          0.1416589574422322D+00, &
          0.1416589574422322D+00, &
          0.1416589574422322D+00, &
          0.1416589574422322D+00, &
          0.8401797426849694D+00, &
          0.8401797426849694D+00, &
          0.8401797426849694D+00, &
          0.8401797426849694D+00, &
          0.8401797426849694D+00, &
          0.8401797426849694D+00, &
          0.1598202573150306D+00, &
          0.1598202573150306D+00, &
          0.1598202573150306D+00, &
          0.1598202573150306D+00, &
          0.1598202573150306D+00, &
          0.1598202573150306D+00, &
          0.6433267612797380D+00, &
          0.6433267612797380D+00, &
          0.6433267612797380D+00, &
          0.6433267612797380D+00, &
          0.6433267612797380D+00, &
          0.6433267612797380D+00, &
          0.3566732387202620D+00, &
          0.3566732387202620D+00, &
          0.3566732387202620D+00, &
          0.3566732387202620D+00, &
          0.3566732387202620D+00, &
          0.3566732387202620D+00, &
          0.9914944215982142D+00, &
          0.9914944215982142D+00, &
          0.9914944215982142D+00, &
          0.9914944215982142D+00, &
          0.9914944215982142D+00, &
          0.9914944215982142D+00, &
          0.0085055784017858D+00, &
          0.0085055784017858D+00, &
          0.0085055784017858D+00, &
          0.0085055784017858D+00, &
          0.0085055784017858D+00, &
          0.0085055784017858D+00, &
          0.9917464272284122D+00, &
          0.9917464272284122D+00, &
          0.9917464272284122D+00, &
          0.9917464272284122D+00, &
          0.9917464272284122D+00, &
          0.9917464272284122D+00, &
          0.0082535727715877D+00, &
          0.0082535727715877D+00, &
          0.0082535727715877D+00, &
          0.0082535727715877D+00, &
          0.0082535727715877D+00, &
          0.0082535727715877D+00, &
          0.9953960310331887D+00, &
          0.9953960310331887D+00, &
          0.9953960310331887D+00, &
          0.9953960310331887D+00, &
          0.9953960310331887D+00, &
          0.9953960310331887D+00, &
          0.0046039689668113D+00, &
          0.0046039689668113D+00, &
          0.0046039689668113D+00, &
          0.0046039689668113D+00, &
          0.0046039689668113D+00, &
          0.0046039689668113D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          0.0047470442329770D+00, &
          0.0047470442329770D+00, &
          0.0006640322548036D+00, &
          0.0006640322548036D+00, &
          0.0006640322548036D+00, &
          0.0009620914054916D+00, &
          0.0009620914054916D+00, &
          0.0009620914054916D+00, &
          0.0027421731795808D+00, &
          0.0027421731795808D+00, &
          0.0027421731795808D+00, &
          0.0027421731795808D+00, &
          0.0027421731795808D+00, &
          0.0027421731795808D+00, &
          0.0028143241025262D+00, &
          0.0028143241025262D+00, &
          0.0028143241025262D+00, &
          0.0028143241025262D+00, &
          0.0028143241025262D+00, &
          0.0028143241025262D+00, &
          0.0038647591741263D+00, &
          0.0038647591741263D+00, &
          0.0038647591741263D+00, &
          0.0038647591741263D+00, &
          0.0038647591741263D+00, &
          0.0038647591741263D+00, &
          0.0037745380549365D+00, &
          0.0037745380549365D+00, &
          0.0037745380549365D+00, &
          0.0037745380549365D+00, &
          0.0037745380549365D+00, &
          0.0037745380549365D+00, &
          0.0002695893315348D+00, &
          0.0002695893315348D+00, &
          0.0002695893315348D+00, &
          0.0002695893315348D+00, &
          0.0002695893315348D+00, &
          0.0002695893315348D+00, &
          0.0011175414425365D+00, &
          0.0011175414425365D+00, &
          0.0011175414425365D+00, &
          0.0011175414425365D+00, &
          0.0011175414425365D+00, &
          0.0011175414425365D+00, &
          0.0011681839673208D+00, &
          0.0011681839673208D+00, &
          0.0011681839673208D+00, &
          0.0011681839673208D+00, &
          0.0011681839673208D+00, &
          0.0011681839673208D+00, &
          0.0009039109611475D+00, &
          0.0009039109611475D+00, &
          0.0009039109611475D+00, &
          0.0009039109611475D+00, &
          0.0009039109611475D+00, &
          0.0009039109611475D+00, &
          0.0033698332015479D+00, &
          0.0033698332015479D+00, &
          0.0033698332015479D+00, &
          0.0033698332015479D+00, &
          0.0033698332015479D+00, &
          0.0033698332015479D+00, &
          0.0026397641711200D+00, &
          0.0026397641711200D+00, &
          0.0026397641711200D+00, &
          0.0026397641711200D+00, &
          0.0026397641711200D+00, &
          0.0026397641711200D+00, &
          0.0025993019789248D+00, &
          0.0025993019789248D+00, &
          0.0025993019789248D+00, &
          0.0025993019789248D+00, &
          0.0025993019789248D+00, &
          0.0025993019789248D+00, &
          0.0012822330144373D+00, &
          0.0012822330144373D+00, &
          0.0012822330144373D+00, &
          0.0012822330144373D+00, &
          0.0012822330144373D+00, &
          0.0012822330144373D+00, &
          0.0012901963397495D+00, &
          0.0012901963397495D+00, &
          0.0012901963397495D+00, &
          0.0012901963397495D+00, &
          0.0012901963397495D+00, &
          0.0012901963397495D+00, &
          0.0036340343281761D+00, &
          0.0036340343281761D+00, &
          0.0036340343281761D+00, &
          0.0036340343281761D+00, &
          0.0036340343281761D+00, &
          0.0036340343281761D+00, &
          0.0038655552411094D+00, &
          0.0038655552411094D+00, &
          0.0038655552411094D+00, &
          0.0038655552411094D+00, &
          0.0038655552411094D+00, &
          0.0038655552411094D+00, &
          0.0010549067548651D+00, &
          0.0010549067548651D+00, &
          0.0010549067548651D+00, &
          0.0010549067548651D+00, &
          0.0010549067548651D+00, &
          0.0010549067548651D+00, &
          0.0015623117217385D+00, &
          0.0015623117217385D+00, &
          0.0015623117217385D+00, &
          0.0015623117217385D+00, &
          0.0015623117217385D+00, &
          0.0015623117217385D+00, &
          0.0004729101647326D+00, &
          0.0004729101647326D+00, &
          0.0004729101647326D+00, &
          0.0004729101647326D+00, &
          0.0004729101647326D+00, &
          0.0004729101647326D+00, &
          0.0001740638131578D+00, &
          0.0001740638131578D+00, &
          0.0001740638131578D+00, &
          0.0001740638131578D+00, &
          0.0001740638131578D+00, &
          0.0001740638131578D+00, &
          0.0028959599653198D+00, &
          0.0028959599653198D+00, &
          0.0028959599653198D+00, &
          0.0028959599653198D+00, &
          0.0028959599653198D+00, &
          0.0028959599653198D+00, &
          0.0052347442788098D+00, &
          0.0052347442788098D+00, &
          0.0052347442788098D+00, &
          0.0052347442788098D+00, &
          0.0052347442788098D+00, &
          0.0052347442788098D+00, &
          0.0031766394449831D+00, &
          0.0031766394449831D+00, &
          0.0031766394449831D+00, &
          0.0031766394449831D+00, &
          0.0031766394449831D+00, &
          0.0031766394449831D+00, &
          0.0027162231233183D+00, &
          0.0027162231233183D+00, &
          0.0027162231233183D+00, &
          0.0027162231233183D+00, &
          0.0027162231233183D+00, &
          0.0027162231233183D+00, &
          0.0032818257977382D+00, &
          0.0032818257977382D+00, &
          0.0032818257977382D+00, &
          0.0032818257977382D+00, &
          0.0032818257977382D+00, &
          0.0032818257977382D+00, &
          0.0031878946984893D+00, &
          0.0031878946984893D+00, &
          0.0031878946984893D+00, &
          0.0031878946984893D+00, &
          0.0031878946984893D+00, &
          0.0031878946984893D+00, &
          0.0035903411322359D+00, &
          0.0035903411322359D+00, &
          0.0035903411322359D+00, &
          0.0035903411322359D+00, &
          0.0035903411322359D+00, &
          0.0035903411322359D+00, &
          0.0035903411322359D+00, &
          0.0035903411322359D+00, &
          0.0035903411322359D+00, &
          0.0035903411322359D+00, &
          0.0035903411322359D+00, &
          0.0035903411322359D+00, &
          0.0005184823629820D+00, &
          0.0005184823629820D+00, &
          0.0005184823629820D+00, &
          0.0005184823629820D+00, &
          0.0005184823629820D+00, &
          0.0005184823629820D+00, &
          0.0005184823629820D+00, &
          0.0005184823629820D+00, &
          0.0005184823629820D+00, &
          0.0005184823629820D+00, &
          0.0005184823629820D+00, &
          0.0005184823629820D+00, &
          0.0020440538003019D+00, &
          0.0020440538003019D+00, &
          0.0020440538003019D+00, &
          0.0020440538003019D+00, &
          0.0020440538003019D+00, &
          0.0020440538003019D+00, &
          0.0020440538003019D+00, &
          0.0020440538003019D+00, &
          0.0020440538003019D+00, &
          0.0020440538003019D+00, &
          0.0020440538003019D+00, &
          0.0020440538003019D+00, &
          0.0024543847025850D+00, &
          0.0024543847025850D+00, &
          0.0024543847025850D+00, &
          0.0024543847025850D+00, &
          0.0024543847025850D+00, &
          0.0024543847025850D+00, &
          0.0024543847025850D+00, &
          0.0024543847025850D+00, &
          0.0024543847025850D+00, &
          0.0024543847025850D+00, &
          0.0024543847025850D+00, &
          0.0024543847025850D+00, &
          0.0041194067328505D+00, &
          0.0041194067328505D+00, &
          0.0041194067328505D+00, &
          0.0041194067328505D+00, &
          0.0041194067328505D+00, &
          0.0041194067328505D+00, &
          0.0041194067328505D+00, &
          0.0041194067328505D+00, &
          0.0041194067328505D+00, &
          0.0041194067328505D+00, &
          0.0041194067328505D+00, &
          0.0041194067328505D+00, &
          0.0057629469800418D+00, &
          0.0057629469800418D+00, &
          0.0057629469800418D+00, &
          0.0057629469800418D+00, &
          0.0057629469800418D+00, &
          0.0057629469800418D+00, &
          0.0057629469800418D+00, &
          0.0057629469800418D+00, &
          0.0057629469800418D+00, &
          0.0057629469800418D+00, &
          0.0057629469800418D+00, &
          0.0057629469800418D+00, &
          0.0032546495070173D+00, &
          0.0032546495070173D+00, &
          0.0032546495070173D+00, &
          0.0032546495070173D+00, &
          0.0032546495070173D+00, &
          0.0032546495070173D+00, &
          0.0032546495070173D+00, &
          0.0032546495070173D+00, &
          0.0032546495070173D+00, &
          0.0032546495070173D+00, &
          0.0032546495070173D+00, &
          0.0032546495070173D+00, &
          0.0013769528611881D+00, &
          0.0013769528611881D+00, &
          0.0013769528611881D+00, &
          0.0013769528611881D+00, &
          0.0013769528611881D+00, &
          0.0013769528611881D+00, &
          0.0013769528611881D+00, &
          0.0013769528611881D+00, &
          0.0013769528611881D+00, &
          0.0013769528611881D+00, &
          0.0013769528611881D+00, &
          0.0013769528611881D+00, &
          0.0013901442662305D+00, &
          0.0013901442662305D+00, &
          0.0013901442662305D+00, &
          0.0013901442662305D+00, &
          0.0013901442662305D+00, &
          0.0013901442662305D+00, &
          0.0013901442662305D+00, &
          0.0013901442662305D+00, &
          0.0013901442662305D+00, &
          0.0013901442662305D+00, &
          0.0013901442662305D+00, &
          0.0013901442662305D+00, &
          0.0019734826201223D+00, &
          0.0019734826201223D+00, &
          0.0019734826201223D+00, &
          0.0019734826201223D+00, &
          0.0019734826201223D+00, &
          0.0019734826201223D+00, &
          0.0019734826201223D+00, &
          0.0019734826201223D+00, &
          0.0019734826201223D+00, &
          0.0019734826201223D+00, &
          0.0019734826201223D+00, &
          0.0019734826201223D+00, &
          0.0011120290537863D+00, &
          0.0011120290537863D+00, &
          0.0011120290537863D+00, &
          0.0011120290537863D+00, &
          0.0011120290537863D+00, &
          0.0011120290537863D+00, &
          0.0011120290537863D+00, &
          0.0011120290537863D+00, &
          0.0011120290537863D+00, &
          0.0011120290537863D+00, &
          0.0011120290537863D+00, &
          0.0011120290537863D+00, &
          0.0011961101336802D+00, &
          0.0011961101336802D+00, &
          0.0011961101336802D+00, &
          0.0011961101336802D+00, &
          0.0011961101336802D+00, &
          0.0011961101336802D+00, &
          0.0011961101336802D+00, &
          0.0011961101336802D+00, &
          0.0011961101336802D+00, &
          0.0011961101336802D+00, &
          0.0011961101336802D+00, &
          0.0011961101336802D+00, &
          0.0018889821414862D+00, &
          0.0018889821414862D+00, &
          0.0018889821414862D+00, &
          0.0018889821414862D+00, &
          0.0018889821414862D+00, &
          0.0018889821414862D+00, &
          0.0018889821414862D+00, &
          0.0018889821414862D+00, &
          0.0018889821414862D+00, &
          0.0018889821414862D+00, &
          0.0018889821414862D+00, &
          0.0018889821414862D+00, &
          0.0005659037184544D+00, &
          0.0005659037184544D+00, &
          0.0005659037184544D+00, &
          0.0005659037184544D+00, &
          0.0005659037184544D+00, &
          0.0005659037184544D+00, &
          0.0005659037184544D+00, &
          0.0005659037184544D+00, &
          0.0005659037184544D+00, &
          0.0005659037184544D+00, &
          0.0005659037184544D+00, &
          0.0005659037184544D+00, &
          0.0018321887972236D+00, &
          0.0018321887972236D+00, &
          0.0018321887972236D+00, &
          0.0018321887972236D+00, &
          0.0018321887972236D+00, &
          0.0018321887972236D+00, &
          0.0018321887972236D+00, &
          0.0018321887972236D+00, &
          0.0018321887972236D+00, &
          0.0018321887972236D+00, &
          0.0018321887972236D+00, &
          0.0018321887972236D+00, &
          0.0008460253945738D+00, &
          0.0008460253945738D+00, &
          0.0008460253945738D+00, &
          0.0008460253945738D+00, &
          0.0008460253945738D+00, &
          0.0008460253945738D+00, &
          0.0008460253945738D+00, &
          0.0008460253945738D+00, &
          0.0008460253945738D+00, &
          0.0008460253945738D+00, &
          0.0008460253945738D+00, &
          0.0008460253945738D+00, &
          0.0015975131973438D+00, &
          0.0015975131973438D+00, &
          0.0015975131973438D+00, &
          0.0015975131973438D+00, &
          0.0015975131973438D+00, &
          0.0015975131973438D+00, &
          0.0015975131973438D+00, &
          0.0015975131973438D+00, &
          0.0015975131973438D+00, &
          0.0015975131973438D+00, &
          0.0015975131973438D+00, &
          0.0015975131973438D+00, &
          0.0026781578703448D+00, &
          0.0026781578703448D+00, &
          0.0026781578703448D+00, &
          0.0026781578703448D+00, &
          0.0026781578703448D+00, &
          0.0026781578703448D+00, &
          0.0026781578703448D+00, &
          0.0026781578703448D+00, &
          0.0026781578703448D+00, &
          0.0026781578703448D+00, &
          0.0026781578703448D+00, &
          0.0026781578703448D+00, &
          0.0007373212126546D+00, &
          0.0007373212126546D+00, &
          0.0007373212126546D+00, &
          0.0007373212126546D+00, &
          0.0007373212126546D+00, &
          0.0007373212126546D+00, &
          0.0007373212126546D+00, &
          0.0007373212126546D+00, &
          0.0007373212126546D+00, &
          0.0007373212126546D+00, &
          0.0007373212126546D+00, &
          0.0007373212126546D+00, &
          0.0013900594133789D+00, &
          0.0013900594133789D+00, &
          0.0013900594133789D+00, &
          0.0013900594133789D+00, &
          0.0013900594133789D+00, &
          0.0013900594133789D+00, &
          0.0013900594133789D+00, &
          0.0013900594133789D+00, &
          0.0013900594133789D+00, &
          0.0013900594133789D+00, &
          0.0013900594133789D+00, &
          0.0013900594133789D+00, &
          0.0007964285504550D+00, &
          0.0007964285504550D+00, &
          0.0007964285504550D+00, &
          0.0007964285504550D+00, &
          0.0007964285504550D+00, &
          0.0007964285504550D+00, &
          0.0007964285504550D+00, &
          0.0007964285504550D+00, &
          0.0007964285504550D+00, &
          0.0007964285504550D+00, &
          0.0007964285504550D+00, &
          0.0007964285504550D+00, &
          0.0014601164425740D+00, &
          0.0014601164425740D+00, &
          0.0014601164425740D+00, &
          0.0014601164425740D+00, &
          0.0014601164425740D+00, &
          0.0014601164425740D+00, &
          0.0014601164425740D+00, &
          0.0014601164425740D+00, &
          0.0014601164425740D+00, &
          0.0014601164425740D+00, &
          0.0014601164425740D+00, &
          0.0014601164425740D+00, &
          0.0013483736840244D+00, &
          0.0013483736840244D+00, &
          0.0013483736840244D+00, &
          0.0013483736840244D+00, &
          0.0013483736840244D+00, &
          0.0013483736840244D+00, &
          0.0013483736840244D+00, &
          0.0013483736840244D+00, &
          0.0013483736840244D+00, &
          0.0013483736840244D+00, &
          0.0013483736840244D+00, &
          0.0013483736840244D+00, &
          0.0032860508781285D+00, &
          0.0032860508781285D+00, &
          0.0032860508781285D+00, &
          0.0032860508781285D+00, &
          0.0032860508781285D+00, &
          0.0032860508781285D+00, &
          0.0032860508781285D+00, &
          0.0032860508781285D+00, &
          0.0032860508781285D+00, &
          0.0032860508781285D+00, &
          0.0032860508781285D+00, &
          0.0032860508781285D+00, &
          0.0011247201187777D+00, &
          0.0011247201187777D+00, &
          0.0011247201187777D+00, &
          0.0011247201187777D+00, &
          0.0011247201187777D+00, &
          0.0011247201187777D+00, &
          0.0011247201187777D+00, &
          0.0011247201187777D+00, &
          0.0011247201187777D+00, &
          0.0011247201187777D+00, &
          0.0011247201187777D+00, &
          0.0011247201187777D+00, &
          0.0021271072748908D+00, &
          0.0021271072748908D+00, &
          0.0021271072748908D+00, &
          0.0021271072748908D+00, &
          0.0021271072748908D+00, &
          0.0021271072748908D+00, &
          0.0021271072748908D+00, &
          0.0021271072748908D+00, &
          0.0021271072748908D+00, &
          0.0021271072748908D+00, &
          0.0021271072748908D+00, &
          0.0021271072748908D+00, &
          0.0007135219223491D+00, &
          0.0007135219223491D+00, &
          0.0007135219223491D+00, &
          0.0007135219223491D+00, &
          0.0007135219223491D+00, &
          0.0007135219223491D+00, &
          0.0007135219223491D+00, &
          0.0007135219223491D+00, &
          0.0007135219223491D+00, &
          0.0007135219223491D+00, &
          0.0007135219223491D+00, &
          0.0007135219223491D+00, &
          0.0003719588400836D+00, &
          0.0003719588400836D+00, &
          0.0003719588400836D+00, &
          0.0003719588400836D+00, &
          0.0003719588400836D+00, &
          0.0003719588400836D+00, &
          0.0003719588400836D+00, &
          0.0003719588400836D+00, &
          0.0003719588400836D+00, &
          0.0003719588400836D+00, &
          0.0003719588400836D+00, &
          0.0003719588400836D+00, &
          0.0002151829924432D+00, &
          0.0002151829924432D+00, &
          0.0002151829924432D+00, &
          0.0002151829924432D+00, &
          0.0002151829924432D+00, &
          0.0002151829924432D+00, &
          0.0002151829924432D+00, &
          0.0002151829924432D+00, &
          0.0002151829924432D+00, &
          0.0002151829924432D+00, &
          0.0002151829924432D+00, &
          0.0002151829924432D+00, &
          0.0008163226512583D+00, &
          0.0008163226512583D+00, &
          0.0008163226512583D+00, &
          0.0008163226512583D+00, &
          0.0008163226512583D+00, &
          0.0008163226512583D+00, &
          0.0008163226512583D+00, &
          0.0008163226512583D+00, &
          0.0008163226512583D+00, &
          0.0008163226512583D+00, &
          0.0008163226512583D+00, &
          0.0008163226512583D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end

