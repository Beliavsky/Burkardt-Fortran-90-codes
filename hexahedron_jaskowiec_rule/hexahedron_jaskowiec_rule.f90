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
subroutine hexahedron_jaskowiec_rule ( p, n, x, y, z, w )

!*****************************************************************************80
!
!! hexahedron_jaskowiec_rule() returns a hexahedron quadrature rule of given precision.
!
!  Discussion:
!
!    The unit hexahedron H is bounded by the vertices
!    (0,0,0), (1,0,0), (1,1,0), (0,1,0), (0,0,1), (1,0,1), (1,1,1), (0,1,1).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 May 2023
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
!    int p: the precision, 0 <= p <= 21.
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
    write ( *, '(a)' ) 'hexahedron_jaskowiec_rule(): Fatal error!'
    write ( *, '(a)' ) '  Input p < 0.'
    stop ( 1 )
  end if

  if ( 21 < p ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'hexahedron_jaskowiec_rule(): Fatal error!'
    write ( *, '(a)' ) '  Input 21 < p.'
    stop ( 1 )
  end if

  if ( p <= 1 ) then
    call rule01 ( n, x, y, z, w )
  else if ( p <= 3 ) then
    call rule03 ( n, x, y, z, w )
  else if ( p <= 5 ) then
    call rule05 ( n, x, y, z, w )
  else if ( p <= 7 ) then
    call rule07 ( n, x, y, z, w )
  else if ( p <= 9 ) then
    call rule09 ( n, x, y, z, w )
  else if ( p <= 11 ) then
    call rule11 ( n, x, y, z, w )
  else if ( p <= 13 ) then
    call rule13 ( n, x, y, z, w )
  else if ( p <= 15 ) then
    call rule15 ( n, x, y, z, w )
  else if ( p <= 17 ) then
    call rule17 ( n, x, y, z, w )
  else if ( p <= 19 ) then
    call rule19 ( n, x, y, z, w )
  else if ( p <= 21 ) then
    call rule21 ( n, x, y, z, w )
  end if

  return
end
subroutine hexahedron_unit_monomial_integral ( expon, value )

!*****************************************************************************80
!
!! hexahedron_unit_monomial_integral(): monomial integral in a unit hexahedron.
!
!  Discussion:
!
!    This routine returns the integral of
!
!      product ( 1 <= I <= 3 ) X(I)^EXPON(I)
!
!    over the unit hexahedron.
!
!    The unit hexahedron H is bounded by the vertices
!    (0,0,0), (1,0,0), (1,1,0), (0,1,0), (0,0,1), (1,0,1), (1,1,1), (0,1,1).
! 
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 May 2023
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer EXPON(3): the exponents.
!
!  Output:
!
!    real VALUE: the integral of the monomial.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer expon(3)
  real ( kind = rk ) value

  value = 1.0D+00 / ( expon(1) + 1 ) &
        * 1.0D+00 / ( expon(2) + 1 ) &
        * 1.0D+00 / ( expon(3) + 1 )

  return
end
function hexahedron_unit_volume ( )

!*****************************************************************************80
!
!! hexahedron_unit_volume() returns the volume of a unit hexahedron.
!
!  Discussion:
!
!    The unit hexahedron has vertices 
!    (0,0,0), (1,0,0), (1,1,0), (0,1,0), (0,0,1), (1,0,1), (1,1,1), (0,1,1).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 May 2023
!
!  Author:
!
!    John Burkardt
!
!  Output:
!
!    real hexahedron_unit_volume: the volume.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) hexahedron_unit_volume

  hexahedron_unit_volume = 1.0

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
subroutine rule_order ( p, order )

!*****************************************************************************80
!
!! rule_order() returns the order of a hexahedron quadrature rule of given precision.
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
!    05 May 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order symmetric cubature rules for tetrahedra and pyramids,
!    International Journal for Numerical Methods in Engineering,
!    Volume 122, Number 1, pages 148-171, 24 August 2020.
!
!  Input:
!
!    integer p: the precision, 0 <= p <= 21.
!
!  Output:
!
!    integer order: the order of the rule.
!
  implicit none

  integer order
  integer, save, dimension(0:21) :: order_vec = (/ &
        1,  &
        1,   6,   6,  14,  14,  34,  34,  58,  58,  90, &
       90, 154, 154, 256, 256, 346, 346, 454, 454, 580, &
      580 /)
  integer p

  if ( p < 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'rule_order(): Fatal error!'
    write ( *, '(a)' ) '  Input p < 0.'
  end if

  if ( 21 < p ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'rule_order(): Fatal error!'
    write ( *, '(a)' ) '  Input 21 < p'
  end if

  order = order_vec(p)

  return
end
subroutine rule01 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule01() returns the hexahedron quadrature rule of precision 1.
!
!  Discussion:
!
!    We suppose we are given a hexahedron H with vertices
!    (0,0,0), (1,0,0), (1,1,0), (0,1,0), (0,0,1), (1,0,1), (1,1,1), (0,1,1),
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over H is 
!    approximated by Q as follows:
!
!    Q = volume(H) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 May 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order symmetric cubature rules for tetrahedra and pyramids,
!    International Journal for Numerical Methods in Engineering,
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
!    real w(n): the quadrature weights, which add to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 1

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
          0.5000000000000000D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.5000000000000000D+00 /)

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
subroutine rule03 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule03() returns the hexahedron quadrature rule of precision 3.
!
!  Discussion:
!
!    We suppose we are given a hexahedron H with vertices
!    (0,0,0), (1,0,0), (1,1,0), (0,1,0), (0,0,1), (1,0,1), (1,1,1), (0,1,1),
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over H is 
!    approximated by Q as follows:
!
!    Q = volume(H) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 May 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order symmetric cubature rules for tetrahedra and pyramids,
!    International Journal for Numerical Methods in Engineering,
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
!    real w(n): the quadrature weights, which add to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 6

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
          1.0000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.0000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.5000000000000000D+00, &
          1.0000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.0000000000000000D+00, &
          0.5000000000000000D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          1.0000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.0000000000000000D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          0.1666666666666667D+00, &
          0.1666666666666667D+00, &
          0.1666666666666667D+00, &
          0.1666666666666667D+00, &
          0.1666666666666667D+00, &
          0.1666666666666667D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule05 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule05() returns the hexahedron quadrature rule of precision 5.
!
!  Discussion:
!
!    We suppose we are given a hexahedron H with vertices
!    (0,0,0), (1,0,0), (1,1,0), (0,1,0), (0,0,1), (1,0,1), (1,1,1), (0,1,1),
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over H is 
!    approximated by Q as follows:
!
!    Q = volume(H) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 May 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order symmetric cubature rules for tetrahedra and pyramids,
!    International Journal for Numerical Methods in Engineering,
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
!    real w(n): the quadrature weights, which add to 1.
!
   implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 14

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
          0.1020887871228893D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.8979112128771107D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.8793934553196641D+00, &
          0.1206065446803359D+00, &
          0.1206065446803359D+00, &
          0.8793934553196641D+00, &
          0.8793934553196641D+00, &
          0.1206065446803359D+00, &
          0.8793934553196641D+00, &
          0.1206065446803359D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.5000000000000000D+00, &
          0.1020887871228893D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.8979112128771107D+00, &
          0.5000000000000000D+00, &
          0.8793934553196641D+00, &
          0.1206065446803359D+00, &
          0.8793934553196641D+00, &
          0.1206065446803359D+00, &
          0.8793934553196641D+00, &
          0.1206065446803359D+00, &
          0.1206065446803359D+00, &
          0.8793934553196641D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.1020887871228893D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.8979112128771107D+00, &
          0.8793934553196641D+00, &
          0.1206065446803359D+00, &
          0.8793934553196641D+00, &
          0.8793934553196641D+00, &
          0.1206065446803359D+00, &
          0.8793934553196641D+00, &
          0.1206065446803359D+00, &
          0.1206065446803359D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          0.1108033240997230D+00, &
          0.1108033240997230D+00, &
          0.1108033240997230D+00, &
          0.1108033240997230D+00, &
          0.1108033240997230D+00, &
          0.1108033240997230D+00, &
          0.0418975069252078D+00, &
          0.0418975069252078D+00, &
          0.0418975069252078D+00, &
          0.0418975069252078D+00, &
          0.0418975069252078D+00, &
          0.0418975069252078D+00, &
          0.0418975069252078D+00, &
          0.0418975069252078D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule07 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule07() returns the hexahedron quadrature rule of precision 7.
!
!  Discussion:
!
!    We suppose we are given a hexahedron H with vertices
!    (0,0,0), (1,0,0), (1,1,0), (0,1,0), (0,0,1), (1,0,1), (1,1,1), (0,1,1),
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over H is 
!    approximated by Q as follows:
!
!    Q = volume(H) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 May 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order symmetric cubature rules for tetrahedra and pyramids,
!    International Journal for Numerical Methods in Engineering,
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
!    real w(n): the quadrature weights, which add to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 34

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
          0.0305973139469912D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9694026860530087D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.8731668050064080D+00, &
          0.1268331949935920D+00, &
          0.1268331949935920D+00, &
          0.8731668050064080D+00, &
          0.8731668050064080D+00, &
          0.1268331949935920D+00, &
          0.8731668050064080D+00, &
          0.1268331949935920D+00, &
          0.7050654491660072D+00, &
          0.2949345508339928D+00, &
          0.2949345508339928D+00, &
          0.7050654491660072D+00, &
          0.7050654491660072D+00, &
          0.2949345508339928D+00, &
          0.7050654491660072D+00, &
          0.2949345508339928D+00, &
          0.9532303450614186D+00, &
          0.5000000000000000D+00, &
          0.9532303450614186D+00, &
          0.0467696549385814D+00, &
          0.5000000000000000D+00, &
          0.0467696549385814D+00, &
          0.9532303450614186D+00, &
          0.5000000000000000D+00, &
          0.9532303450614186D+00, &
          0.0467696549385814D+00, &
          0.5000000000000000D+00, &
          0.0467696549385814D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.5000000000000000D+00, &
          0.0305973139469912D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9694026860530087D+00, &
          0.5000000000000000D+00, &
          0.8731668050064080D+00, &
          0.1268331949935920D+00, &
          0.8731668050064080D+00, &
          0.1268331949935920D+00, &
          0.8731668050064080D+00, &
          0.1268331949935920D+00, &
          0.1268331949935920D+00, &
          0.8731668050064080D+00, &
          0.7050654491660072D+00, &
          0.2949345508339928D+00, &
          0.7050654491660072D+00, &
          0.2949345508339928D+00, &
          0.7050654491660072D+00, &
          0.2949345508339928D+00, &
          0.2949345508339928D+00, &
          0.7050654491660072D+00, &
          0.9532303450614186D+00, &
          0.9532303450614186D+00, &
          0.5000000000000000D+00, &
          0.9532303450614186D+00, &
          0.0467696549385814D+00, &
          0.5000000000000000D+00, &
          0.0467696549385814D+00, &
          0.9532303450614186D+00, &
          0.5000000000000000D+00, &
          0.0467696549385814D+00, &
          0.0467696549385814D+00, &
          0.5000000000000000D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.0305973139469912D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9694026860530087D+00, &
          0.8731668050064080D+00, &
          0.1268331949935920D+00, &
          0.8731668050064080D+00, &
          0.8731668050064080D+00, &
          0.1268331949935920D+00, &
          0.8731668050064080D+00, &
          0.1268331949935920D+00, &
          0.1268331949935920D+00, &
          0.7050654491660072D+00, &
          0.2949345508339928D+00, &
          0.7050654491660072D+00, &
          0.7050654491660072D+00, &
          0.2949345508339928D+00, &
          0.7050654491660072D+00, &
          0.2949345508339928D+00, &
          0.2949345508339928D+00, &
          0.5000000000000000D+00, &
          0.9532303450614186D+00, &
          0.9532303450614186D+00, &
          0.5000000000000000D+00, &
          0.9532303450614186D+00, &
          0.9532303450614186D+00, &
          0.5000000000000000D+00, &
          0.0467696549385814D+00, &
          0.0467696549385814D+00, &
          0.5000000000000000D+00, &
          0.0467696549385814D+00, &
          0.0467696549385814D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          0.0340045847498003D+00, &
          0.0340045847498003D+00, &
          0.0340045847498003D+00, &
          0.0340045847498003D+00, &
          0.0340045847498003D+00, &
          0.0340045847498003D+00, &
          0.0252967244013519D+00, &
          0.0252967244013519D+00, &
          0.0252967244013519D+00, &
          0.0252967244013519D+00, &
          0.0252967244013519D+00, &
          0.0252967244013519D+00, &
          0.0252967244013519D+00, &
          0.0252967244013519D+00, &
          0.0541706353348564D+00, &
          0.0541706353348564D+00, &
          0.0541706353348564D+00, &
          0.0541706353348564D+00, &
          0.0541706353348564D+00, &
          0.0541706353348564D+00, &
          0.0541706353348564D+00, &
          0.0541706353348564D+00, &
          0.0133528011342943D+00, &
          0.0133528011342943D+00, &
          0.0133528011342943D+00, &
          0.0133528011342943D+00, &
          0.0133528011342943D+00, &
          0.0133528011342943D+00, &
          0.0133528011342943D+00, &
          0.0133528011342943D+00, &
          0.0133528011342943D+00, &
          0.0133528011342943D+00, &
          0.0133528011342943D+00, &
          0.0133528011342943D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule09 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule09() returns the hexahedron quadrature rule of precision 9.
!
!  Discussion:
!
!    We suppose we are given a hexahedron H with vertices
!    (0,0,0), (1,0,0), (1,1,0), (0,1,0), (0,0,1), (1,0,1), (1,1,1), (0,1,1),
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over H is 
!    approximated by Q as follows:
!
!    Q = volume(H) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 May 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order symmetric cubature rules for tetrahedra and pyramids,
!    International Journal for Numerical Methods in Engineering,
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
!    real w(n): the quadrature weights, which add to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 58

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
          0.1931592652041455D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.8068407347958545D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9350498923309880D+00, &
          0.0649501076690120D+00, &
          0.0649501076690120D+00, &
          0.9350498923309880D+00, &
          0.9350498923309880D+00, &
          0.0649501076690120D+00, &
          0.9350498923309880D+00, &
          0.0649501076690120D+00, &
          0.7820554035100150D+00, &
          0.2179445964899850D+00, &
          0.2179445964899850D+00, &
          0.7820554035100150D+00, &
          0.7820554035100150D+00, &
          0.2179445964899850D+00, &
          0.7820554035100150D+00, &
          0.2179445964899850D+00, &
          0.9388435616288391D+00, &
          0.5000000000000000D+00, &
          0.9388435616288391D+00, &
          0.0611564383711609D+00, &
          0.5000000000000000D+00, &
          0.0611564383711609D+00, &
          0.9388435616288391D+00, &
          0.5000000000000000D+00, &
          0.9388435616288391D+00, &
          0.0611564383711609D+00, &
          0.5000000000000000D+00, &
          0.0611564383711609D+00, &
          0.7161339513154311D+00, &
          0.9692652109323359D+00, &
          0.7161339513154311D+00, &
          0.2838660486845689D+00, &
          0.9692652109323359D+00, &
          0.2838660486845689D+00, &
          0.2838660486845689D+00, &
          0.9692652109323359D+00, &
          0.2838660486845689D+00, &
          0.7161339513154311D+00, &
          0.9692652109323359D+00, &
          0.7161339513154311D+00, &
          0.7161339513154311D+00, &
          0.0307347890676641D+00, &
          0.7161339513154311D+00, &
          0.2838660486845689D+00, &
          0.0307347890676641D+00, &
          0.2838660486845689D+00, &
          0.2838660486845689D+00, &
          0.0307347890676641D+00, &
          0.2838660486845689D+00, &
          0.7161339513154311D+00, &
          0.0307347890676641D+00, &
          0.7161339513154311D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.5000000000000000D+00, &
          0.1931592652041455D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.8068407347958545D+00, &
          0.5000000000000000D+00, &
          0.9350498923309880D+00, &
          0.0649501076690120D+00, &
          0.9350498923309880D+00, &
          0.0649501076690120D+00, &
          0.9350498923309880D+00, &
          0.0649501076690120D+00, &
          0.0649501076690120D+00, &
          0.9350498923309880D+00, &
          0.7820554035100150D+00, &
          0.2179445964899850D+00, &
          0.7820554035100150D+00, &
          0.2179445964899850D+00, &
          0.7820554035100150D+00, &
          0.2179445964899850D+00, &
          0.2179445964899850D+00, &
          0.7820554035100150D+00, &
          0.9388435616288391D+00, &
          0.9388435616288391D+00, &
          0.5000000000000000D+00, &
          0.9388435616288391D+00, &
          0.0611564383711609D+00, &
          0.5000000000000000D+00, &
          0.0611564383711609D+00, &
          0.9388435616288391D+00, &
          0.5000000000000000D+00, &
          0.0611564383711609D+00, &
          0.0611564383711609D+00, &
          0.5000000000000000D+00, &
          0.9692652109323359D+00, &
          0.7161339513154311D+00, &
          0.7161339513154311D+00, &
          0.9692652109323359D+00, &
          0.2838660486845689D+00, &
          0.7161339513154311D+00, &
          0.9692652109323359D+00, &
          0.2838660486845689D+00, &
          0.2838660486845689D+00, &
          0.9692652109323359D+00, &
          0.7161339513154311D+00, &
          0.2838660486845689D+00, &
          0.0307347890676641D+00, &
          0.7161339513154311D+00, &
          0.7161339513154311D+00, &
          0.0307347890676641D+00, &
          0.2838660486845689D+00, &
          0.7161339513154311D+00, &
          0.0307347890676641D+00, &
          0.2838660486845689D+00, &
          0.2838660486845689D+00, &
          0.0307347890676641D+00, &
          0.7161339513154311D+00, &
          0.2838660486845689D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.1931592652041455D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.8068407347958545D+00, &
          0.9350498923309880D+00, &
          0.0649501076690120D+00, &
          0.9350498923309880D+00, &
          0.9350498923309880D+00, &
          0.0649501076690120D+00, &
          0.9350498923309880D+00, &
          0.0649501076690120D+00, &
          0.0649501076690120D+00, &
          0.7820554035100150D+00, &
          0.2179445964899850D+00, &
          0.7820554035100150D+00, &
          0.7820554035100150D+00, &
          0.2179445964899850D+00, &
          0.7820554035100150D+00, &
          0.2179445964899850D+00, &
          0.2179445964899850D+00, &
          0.5000000000000000D+00, &
          0.9388435616288391D+00, &
          0.9388435616288391D+00, &
          0.5000000000000000D+00, &
          0.9388435616288391D+00, &
          0.9388435616288391D+00, &
          0.5000000000000000D+00, &
          0.0611564383711609D+00, &
          0.0611564383711609D+00, &
          0.5000000000000000D+00, &
          0.0611564383711609D+00, &
          0.0611564383711609D+00, &
          0.7161339513154311D+00, &
          0.7161339513154311D+00, &
          0.9692652109323359D+00, &
          0.7161339513154311D+00, &
          0.7161339513154311D+00, &
          0.9692652109323359D+00, &
          0.2838660486845689D+00, &
          0.2838660486845689D+00, &
          0.9692652109323359D+00, &
          0.2838660486845689D+00, &
          0.2838660486845689D+00, &
          0.9692652109323359D+00, &
          0.7161339513154311D+00, &
          0.7161339513154311D+00, &
          0.0307347890676641D+00, &
          0.7161339513154311D+00, &
          0.7161339513154311D+00, &
          0.0307347890676641D+00, &
          0.2838660486845689D+00, &
          0.2838660486845689D+00, &
          0.0307347890676641D+00, &
          0.2838660486845689D+00, &
          0.2838660486845689D+00, &
          0.0307347890676641D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          0.0541593744687068D+00, &
          0.0541593744687068D+00, &
          0.0541593744687068D+00, &
          0.0541593744687068D+00, &
          0.0541593744687068D+00, &
          0.0541593744687068D+00, &
          0.0062685994124186D+00, &
          0.0062685994124186D+00, &
          0.0062685994124186D+00, &
          0.0062685994124186D+00, &
          0.0062685994124186D+00, &
          0.0062685994124186D+00, &
          0.0062685994124186D+00, &
          0.0062685994124186D+00, &
          0.0248574797680029D+00, &
          0.0248574797680029D+00, &
          0.0248574797680029D+00, &
          0.0248574797680029D+00, &
          0.0248574797680029D+00, &
          0.0248574797680029D+00, &
          0.0248574797680029D+00, &
          0.0248574797680029D+00, &
          0.0114737257670222D+00, &
          0.0114737257670222D+00, &
          0.0114737257670222D+00, &
          0.0114737257670222D+00, &
          0.0114737257670222D+00, &
          0.0114737257670222D+00, &
          0.0114737257670222D+00, &
          0.0114737257670222D+00, &
          0.0114737257670222D+00, &
          0.0114737257670222D+00, &
          0.0114737257670222D+00, &
          0.0114737257670222D+00, &
          0.0120146004391717D+00, &
          0.0120146004391717D+00, &
          0.0120146004391717D+00, &
          0.0120146004391717D+00, &
          0.0120146004391717D+00, &
          0.0120146004391717D+00, &
          0.0120146004391717D+00, &
          0.0120146004391717D+00, &
          0.0120146004391717D+00, &
          0.0120146004391717D+00, &
          0.0120146004391717D+00, &
          0.0120146004391717D+00, &
          0.0120146004391717D+00, &
          0.0120146004391717D+00, &
          0.0120146004391717D+00, &
          0.0120146004391717D+00, &
          0.0120146004391717D+00, &
          0.0120146004391717D+00, &
          0.0120146004391717D+00, &
          0.0120146004391717D+00, &
          0.0120146004391717D+00, &
          0.0120146004391717D+00, &
          0.0120146004391717D+00, &
          0.0120146004391717D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule11 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule11() returns the hexahedron quadrature rule of precision 11.
!
!  Discussion:
!
!    We suppose we are given a hexahedron H with vertices
!    (0,0,0), (1,0,0), (1,1,0), (0,1,0), (0,0,1), (1,0,1), (1,1,1), (0,1,1),
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over H is 
!    approximated by Q as follows:
!
!    Q = volume(H) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 May 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order symmetric cubature rules for tetrahedra and pyramids,
!    International Journal for Numerical Methods in Engineering,
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
!    real w(n): the quadrature weights, which add to 1.
!
   implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 90

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
          0.8610665194372092D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.1389334805627908D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9047441009815494D+00, &
          0.0952558990184506D+00, &
          0.0952558990184506D+00, &
          0.9047441009815494D+00, &
          0.9047441009815494D+00, &
          0.0952558990184506D+00, &
          0.9047441009815494D+00, &
          0.0952558990184506D+00, &
          0.2331729955597515D+00, &
          0.7668270044402485D+00, &
          0.7668270044402485D+00, &
          0.2331729955597515D+00, &
          0.2331729955597515D+00, &
          0.7668270044402485D+00, &
          0.2331729955597515D+00, &
          0.7668270044402485D+00, &
          0.6403862933256372D+00, &
          0.3596137066743628D+00, &
          0.3596137066743628D+00, &
          0.6403862933256372D+00, &
          0.6403862933256372D+00, &
          0.3596137066743628D+00, &
          0.6403862933256372D+00, &
          0.3596137066743628D+00, &
          0.9019667336076422D+00, &
          0.5000000000000000D+00, &
          0.9019667336076422D+00, &
          0.0980332663923578D+00, &
          0.5000000000000000D+00, &
          0.0980332663923578D+00, &
          0.9019667336076422D+00, &
          0.5000000000000000D+00, &
          0.9019667336076422D+00, &
          0.0980332663923578D+00, &
          0.5000000000000000D+00, &
          0.0980332663923578D+00, &
          0.9900497455045357D+00, &
          0.2346080844030868D+00, &
          0.9900497455045357D+00, &
          0.0099502544954643D+00, &
          0.2346080844030868D+00, &
          0.0099502544954643D+00, &
          0.0099502544954643D+00, &
          0.2346080844030868D+00, &
          0.0099502544954643D+00, &
          0.9900497455045357D+00, &
          0.2346080844030868D+00, &
          0.9900497455045357D+00, &
          0.9900497455045357D+00, &
          0.7653919155969132D+00, &
          0.9900497455045357D+00, &
          0.0099502544954643D+00, &
          0.7653919155969132D+00, &
          0.0099502544954643D+00, &
          0.0099502544954643D+00, &
          0.7653919155969132D+00, &
          0.0099502544954643D+00, &
          0.9900497455045357D+00, &
          0.7653919155969132D+00, &
          0.9900497455045357D+00, &
          0.7028429900975481D+00, &
          0.9772916094647830D+00, &
          0.7028429900975481D+00, &
          0.2971570099024518D+00, &
          0.9772916094647830D+00, &
          0.2971570099024518D+00, &
          0.2971570099024518D+00, &
          0.9772916094647830D+00, &
          0.2971570099024518D+00, &
          0.7028429900975481D+00, &
          0.9772916094647830D+00, &
          0.7028429900975481D+00, &
          0.7028429900975481D+00, &
          0.0227083905352169D+00, &
          0.7028429900975481D+00, &
          0.2971570099024518D+00, &
          0.0227083905352169D+00, &
          0.2971570099024518D+00, &
          0.2971570099024518D+00, &
          0.0227083905352169D+00, &
          0.2971570099024518D+00, &
          0.7028429900975481D+00, &
          0.0227083905352169D+00, &
          0.7028429900975481D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.5000000000000000D+00, &
          0.8610665194372092D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.1389334805627908D+00, &
          0.5000000000000000D+00, &
          0.9047441009815494D+00, &
          0.0952558990184506D+00, &
          0.9047441009815494D+00, &
          0.0952558990184506D+00, &
          0.9047441009815494D+00, &
          0.0952558990184506D+00, &
          0.0952558990184506D+00, &
          0.9047441009815494D+00, &
          0.2331729955597515D+00, &
          0.7668270044402485D+00, &
          0.2331729955597515D+00, &
          0.7668270044402485D+00, &
          0.2331729955597515D+00, &
          0.7668270044402485D+00, &
          0.7668270044402485D+00, &
          0.2331729955597515D+00, &
          0.6403862933256372D+00, &
          0.3596137066743628D+00, &
          0.6403862933256372D+00, &
          0.3596137066743628D+00, &
          0.6403862933256372D+00, &
          0.3596137066743628D+00, &
          0.3596137066743628D+00, &
          0.6403862933256372D+00, &
          0.9019667336076422D+00, &
          0.9019667336076422D+00, &
          0.5000000000000000D+00, &
          0.9019667336076422D+00, &
          0.0980332663923578D+00, &
          0.5000000000000000D+00, &
          0.0980332663923578D+00, &
          0.9019667336076422D+00, &
          0.5000000000000000D+00, &
          0.0980332663923578D+00, &
          0.0980332663923578D+00, &
          0.5000000000000000D+00, &
          0.2346080844030868D+00, &
          0.9900497455045357D+00, &
          0.9900497455045357D+00, &
          0.2346080844030868D+00, &
          0.0099502544954643D+00, &
          0.9900497455045357D+00, &
          0.2346080844030868D+00, &
          0.0099502544954643D+00, &
          0.0099502544954643D+00, &
          0.2346080844030868D+00, &
          0.9900497455045357D+00, &
          0.0099502544954643D+00, &
          0.7653919155969132D+00, &
          0.9900497455045357D+00, &
          0.9900497455045357D+00, &
          0.7653919155969132D+00, &
          0.0099502544954643D+00, &
          0.9900497455045357D+00, &
          0.7653919155969132D+00, &
          0.0099502544954643D+00, &
          0.0099502544954643D+00, &
          0.7653919155969132D+00, &
          0.9900497455045357D+00, &
          0.0099502544954643D+00, &
          0.9772916094647830D+00, &
          0.7028429900975481D+00, &
          0.7028429900975481D+00, &
          0.9772916094647830D+00, &
          0.2971570099024518D+00, &
          0.7028429900975481D+00, &
          0.9772916094647830D+00, &
          0.2971570099024518D+00, &
          0.2971570099024518D+00, &
          0.9772916094647830D+00, &
          0.7028429900975481D+00, &
          0.2971570099024518D+00, &
          0.0227083905352169D+00, &
          0.7028429900975481D+00, &
          0.7028429900975481D+00, &
          0.0227083905352169D+00, &
          0.2971570099024518D+00, &
          0.7028429900975481D+00, &
          0.0227083905352169D+00, &
          0.2971570099024518D+00, &
          0.2971570099024518D+00, &
          0.0227083905352169D+00, &
          0.7028429900975481D+00, &
          0.2971570099024518D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.8610665194372092D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.1389334805627908D+00, &
          0.9047441009815494D+00, &
          0.0952558990184506D+00, &
          0.9047441009815494D+00, &
          0.9047441009815494D+00, &
          0.0952558990184506D+00, &
          0.9047441009815494D+00, &
          0.0952558990184506D+00, &
          0.0952558990184506D+00, &
          0.2331729955597515D+00, &
          0.7668270044402485D+00, &
          0.2331729955597515D+00, &
          0.2331729955597515D+00, &
          0.7668270044402485D+00, &
          0.2331729955597515D+00, &
          0.7668270044402485D+00, &
          0.7668270044402485D+00, &
          0.6403862933256372D+00, &
          0.3596137066743628D+00, &
          0.6403862933256372D+00, &
          0.6403862933256372D+00, &
          0.3596137066743628D+00, &
          0.6403862933256372D+00, &
          0.3596137066743628D+00, &
          0.3596137066743628D+00, &
          0.5000000000000000D+00, &
          0.9019667336076422D+00, &
          0.9019667336076422D+00, &
          0.5000000000000000D+00, &
          0.9019667336076422D+00, &
          0.9019667336076422D+00, &
          0.5000000000000000D+00, &
          0.0980332663923578D+00, &
          0.0980332663923578D+00, &
          0.5000000000000000D+00, &
          0.0980332663923578D+00, &
          0.0980332663923578D+00, &
          0.9900497455045357D+00, &
          0.9900497455045357D+00, &
          0.2346080844030868D+00, &
          0.9900497455045357D+00, &
          0.9900497455045357D+00, &
          0.2346080844030868D+00, &
          0.0099502544954643D+00, &
          0.0099502544954643D+00, &
          0.2346080844030868D+00, &
          0.0099502544954643D+00, &
          0.0099502544954643D+00, &
          0.2346080844030868D+00, &
          0.9900497455045357D+00, &
          0.9900497455045357D+00, &
          0.7653919155969132D+00, &
          0.9900497455045357D+00, &
          0.9900497455045357D+00, &
          0.7653919155969132D+00, &
          0.0099502544954643D+00, &
          0.0099502544954643D+00, &
          0.7653919155969132D+00, &
          0.0099502544954643D+00, &
          0.0099502544954643D+00, &
          0.7653919155969132D+00, &
          0.7028429900975481D+00, &
          0.7028429900975481D+00, &
          0.9772916094647830D+00, &
          0.7028429900975481D+00, &
          0.7028429900975481D+00, &
          0.9772916094647830D+00, &
          0.2971570099024518D+00, &
          0.2971570099024518D+00, &
          0.9772916094647830D+00, &
          0.2971570099024518D+00, &
          0.2971570099024518D+00, &
          0.9772916094647830D+00, &
          0.7028429900975481D+00, &
          0.7028429900975481D+00, &
          0.0227083905352169D+00, &
          0.7028429900975481D+00, &
          0.7028429900975481D+00, &
          0.0227083905352169D+00, &
          0.2971570099024518D+00, &
          0.2971570099024518D+00, &
          0.0227083905352169D+00, &
          0.2971570099024518D+00, &
          0.2971570099024518D+00, &
          0.0227083905352169D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          0.0296169245240206D+00, &
          0.0296169245240206D+00, &
          0.0296169245240206D+00, &
          0.0296169245240206D+00, &
          0.0296169245240206D+00, &
          0.0296169245240206D+00, &
          0.0077511837875231D+00, &
          0.0077511837875231D+00, &
          0.0077511837875231D+00, &
          0.0077511837875231D+00, &
          0.0077511837875231D+00, &
          0.0077511837875231D+00, &
          0.0077511837875231D+00, &
          0.0077511837875231D+00, &
          0.0222210081090505D+00, &
          0.0222210081090505D+00, &
          0.0222210081090505D+00, &
          0.0222210081090505D+00, &
          0.0222210081090505D+00, &
          0.0222210081090505D+00, &
          0.0222210081090505D+00, &
          0.0222210081090505D+00, &
          0.0169887021445520D+00, &
          0.0169887021445520D+00, &
          0.0169887021445520D+00, &
          0.0169887021445520D+00, &
          0.0169887021445520D+00, &
          0.0169887021445520D+00, &
          0.0169887021445520D+00, &
          0.0169887021445520D+00, &
          0.0160272817769356D+00, &
          0.0160272817769356D+00, &
          0.0160272817769356D+00, &
          0.0160272817769356D+00, &
          0.0160272817769356D+00, &
          0.0160272817769356D+00, &
          0.0160272817769356D+00, &
          0.0160272817769356D+00, &
          0.0160272817769356D+00, &
          0.0160272817769356D+00, &
          0.0160272817769356D+00, &
          0.0160272817769356D+00, &
          0.0016188211900323D+00, &
          0.0016188211900323D+00, &
          0.0016188211900323D+00, &
          0.0016188211900323D+00, &
          0.0016188211900323D+00, &
          0.0016188211900323D+00, &
          0.0016188211900323D+00, &
          0.0016188211900323D+00, &
          0.0016188211900323D+00, &
          0.0016188211900323D+00, &
          0.0016188211900323D+00, &
          0.0016188211900323D+00, &
          0.0016188211900323D+00, &
          0.0016188211900323D+00, &
          0.0016188211900323D+00, &
          0.0016188211900323D+00, &
          0.0016188211900323D+00, &
          0.0016188211900323D+00, &
          0.0016188211900323D+00, &
          0.0016188211900323D+00, &
          0.0016188211900323D+00, &
          0.0016188211900323D+00, &
          0.0016188211900323D+00, &
          0.0016188211900323D+00, &
          0.0089763421101196D+00, &
          0.0089763421101196D+00, &
          0.0089763421101196D+00, &
          0.0089763421101196D+00, &
          0.0089763421101196D+00, &
          0.0089763421101196D+00, &
          0.0089763421101196D+00, &
          0.0089763421101196D+00, &
          0.0089763421101196D+00, &
          0.0089763421101196D+00, &
          0.0089763421101196D+00, &
          0.0089763421101196D+00, &
          0.0089763421101196D+00, &
          0.0089763421101196D+00, &
          0.0089763421101196D+00, &
          0.0089763421101196D+00, &
          0.0089763421101196D+00, &
          0.0089763421101196D+00, &
          0.0089763421101196D+00, &
          0.0089763421101196D+00, &
          0.0089763421101196D+00, &
          0.0089763421101196D+00, &
          0.0089763421101196D+00, &
          0.0089763421101196D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule13 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule13() returns the hexahedron quadrature rule of precision 13.
!
!  Discussion:
!
!    We suppose we are given a hexahedron H with vertices
!    (0,0,0), (1,0,0), (1,1,0), (0,1,0), (0,0,1), (1,0,1), (1,1,1), (0,1,1),
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over H is 
!    approximated by Q as follows:
!
!    Q = volume(H) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 May 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order symmetric cubature rules for tetrahedra and pyramids,
!    International Journal for Numerical Methods in Engineering,
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
!    real w(n): the quadrature weights, which add to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 154

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
          0.1827950053146007D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.8172049946853993D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9539270239814177D+00, &
          0.0460729760185823D+00, &
          0.0460729760185823D+00, &
          0.9539270239814177D+00, &
          0.9539270239814177D+00, &
          0.0460729760185823D+00, &
          0.9539270239814177D+00, &
          0.0460729760185823D+00, &
          0.3818973707741769D+00, &
          0.6181026292258232D+00, &
          0.6181026292258232D+00, &
          0.3818973707741769D+00, &
          0.3818973707741769D+00, &
          0.6181026292258232D+00, &
          0.3818973707741769D+00, &
          0.6181026292258232D+00, &
          0.8687873171566183D+00, &
          0.5000000000000000D+00, &
          0.8687873171566183D+00, &
          0.1312126828433817D+00, &
          0.5000000000000000D+00, &
          0.1312126828433817D+00, &
          0.8687873171566183D+00, &
          0.5000000000000000D+00, &
          0.8687873171566183D+00, &
          0.1312126828433817D+00, &
          0.5000000000000000D+00, &
          0.1312126828433817D+00, &
          0.7303580738482240D+00, &
          0.2696419261517760D+00, &
          0.7303580738482240D+00, &
          0.2696419261517760D+00, &
          0.7303580738482240D+00, &
          0.2696419261517760D+00, &
          0.7303580738482240D+00, &
          0.2696419261517760D+00, &
          0.9666726015054183D+00, &
          0.0333273984945817D+00, &
          0.0333273984945817D+00, &
          0.9666726015054183D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9666726015054183D+00, &
          0.0333273984945817D+00, &
          0.0333273984945817D+00, &
          0.9666726015054183D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.6866920757899998D+00, &
          0.8321081947965769D+00, &
          0.6866920757899998D+00, &
          0.3133079242100002D+00, &
          0.8321081947965769D+00, &
          0.3133079242100002D+00, &
          0.3133079242100002D+00, &
          0.8321081947965769D+00, &
          0.3133079242100002D+00, &
          0.6866920757899998D+00, &
          0.8321081947965769D+00, &
          0.6866920757899998D+00, &
          0.6866920757899998D+00, &
          0.1678918052034231D+00, &
          0.6866920757899998D+00, &
          0.3133079242100002D+00, &
          0.1678918052034231D+00, &
          0.3133079242100002D+00, &
          0.3133079242100002D+00, &
          0.1678918052034231D+00, &
          0.3133079242100002D+00, &
          0.6866920757899998D+00, &
          0.1678918052034231D+00, &
          0.6866920757899998D+00, &
          0.9747555407321407D+00, &
          0.6288381939284606D+00, &
          0.9747555407321407D+00, &
          0.0252444592678593D+00, &
          0.6288381939284606D+00, &
          0.0252444592678593D+00, &
          0.0252444592678593D+00, &
          0.6288381939284606D+00, &
          0.0252444592678593D+00, &
          0.9747555407321407D+00, &
          0.6288381939284606D+00, &
          0.9747555407321407D+00, &
          0.9747555407321407D+00, &
          0.3711618060715394D+00, &
          0.9747555407321407D+00, &
          0.0252444592678593D+00, &
          0.3711618060715394D+00, &
          0.0252444592678593D+00, &
          0.0252444592678593D+00, &
          0.3711618060715394D+00, &
          0.0252444592678593D+00, &
          0.9747555407321407D+00, &
          0.3711618060715394D+00, &
          0.9747555407321407D+00, &
          0.8328347413800889D+00, &
          0.9974350214009486D+00, &
          0.8328347413800889D+00, &
          0.1671652586199110D+00, &
          0.9974350214009486D+00, &
          0.1671652586199110D+00, &
          0.1671652586199110D+00, &
          0.9974350214009486D+00, &
          0.1671652586199110D+00, &
          0.8328347413800889D+00, &
          0.9974350214009486D+00, &
          0.8328347413800889D+00, &
          0.8328347413800889D+00, &
          0.0025649785990514D+00, &
          0.8328347413800889D+00, &
          0.1671652586199110D+00, &
          0.0025649785990514D+00, &
          0.1671652586199110D+00, &
          0.1671652586199110D+00, &
          0.0025649785990514D+00, &
          0.1671652586199110D+00, &
          0.8328347413800889D+00, &
          0.0025649785990514D+00, &
          0.8328347413800889D+00, &
          0.0974495447483840D+00, &
          0.2739297470225850D+00, &
          0.0974495447483840D+00, &
          0.9025504552516159D+00, &
          0.2739297470225850D+00, &
          0.9025504552516159D+00, &
          0.9025504552516159D+00, &
          0.2739297470225850D+00, &
          0.9025504552516159D+00, &
          0.0974495447483840D+00, &
          0.2739297470225850D+00, &
          0.0974495447483840D+00, &
          0.0974495447483840D+00, &
          0.7260702529774150D+00, &
          0.0974495447483840D+00, &
          0.9025504552516159D+00, &
          0.7260702529774150D+00, &
          0.9025504552516159D+00, &
          0.9025504552516159D+00, &
          0.7260702529774150D+00, &
          0.9025504552516159D+00, &
          0.0974495447483840D+00, &
          0.7260702529774150D+00, &
          0.0974495447483840D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.5000000000000000D+00, &
          0.1827950053146007D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.8172049946853993D+00, &
          0.5000000000000000D+00, &
          0.9539270239814177D+00, &
          0.0460729760185823D+00, &
          0.9539270239814177D+00, &
          0.0460729760185823D+00, &
          0.9539270239814177D+00, &
          0.0460729760185823D+00, &
          0.0460729760185823D+00, &
          0.9539270239814177D+00, &
          0.3818973707741769D+00, &
          0.6181026292258232D+00, &
          0.3818973707741769D+00, &
          0.6181026292258232D+00, &
          0.3818973707741769D+00, &
          0.6181026292258232D+00, &
          0.6181026292258232D+00, &
          0.3818973707741769D+00, &
          0.8687873171566183D+00, &
          0.8687873171566183D+00, &
          0.5000000000000000D+00, &
          0.8687873171566183D+00, &
          0.1312126828433817D+00, &
          0.5000000000000000D+00, &
          0.1312126828433817D+00, &
          0.8687873171566183D+00, &
          0.5000000000000000D+00, &
          0.1312126828433817D+00, &
          0.1312126828433817D+00, &
          0.5000000000000000D+00, &
          0.9666726015054183D+00, &
          0.0333273984945817D+00, &
          0.0333273984945817D+00, &
          0.9666726015054183D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.7303580738482240D+00, &
          0.2696419261517760D+00, &
          0.7303580738482240D+00, &
          0.2696419261517760D+00, &
          0.7303580738482240D+00, &
          0.2696419261517760D+00, &
          0.7303580738482240D+00, &
          0.2696419261517760D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9666726015054183D+00, &
          0.0333273984945817D+00, &
          0.0333273984945817D+00, &
          0.9666726015054183D+00, &
          0.8321081947965769D+00, &
          0.6866920757899998D+00, &
          0.6866920757899998D+00, &
          0.8321081947965769D+00, &
          0.3133079242100002D+00, &
          0.6866920757899998D+00, &
          0.8321081947965769D+00, &
          0.3133079242100002D+00, &
          0.3133079242100002D+00, &
          0.8321081947965769D+00, &
          0.6866920757899998D+00, &
          0.3133079242100002D+00, &
          0.1678918052034231D+00, &
          0.6866920757899998D+00, &
          0.6866920757899998D+00, &
          0.1678918052034231D+00, &
          0.3133079242100002D+00, &
          0.6866920757899998D+00, &
          0.1678918052034231D+00, &
          0.3133079242100002D+00, &
          0.3133079242100002D+00, &
          0.1678918052034231D+00, &
          0.6866920757899998D+00, &
          0.3133079242100002D+00, &
          0.6288381939284606D+00, &
          0.9747555407321407D+00, &
          0.9747555407321407D+00, &
          0.6288381939284606D+00, &
          0.0252444592678593D+00, &
          0.9747555407321407D+00, &
          0.6288381939284606D+00, &
          0.0252444592678593D+00, &
          0.0252444592678593D+00, &
          0.6288381939284606D+00, &
          0.9747555407321407D+00, &
          0.0252444592678593D+00, &
          0.3711618060715394D+00, &
          0.9747555407321407D+00, &
          0.9747555407321407D+00, &
          0.3711618060715394D+00, &
          0.0252444592678593D+00, &
          0.9747555407321407D+00, &
          0.3711618060715394D+00, &
          0.0252444592678593D+00, &
          0.0252444592678593D+00, &
          0.3711618060715394D+00, &
          0.9747555407321407D+00, &
          0.0252444592678593D+00, &
          0.9974350214009486D+00, &
          0.8328347413800889D+00, &
          0.8328347413800889D+00, &
          0.9974350214009486D+00, &
          0.1671652586199110D+00, &
          0.8328347413800889D+00, &
          0.9974350214009486D+00, &
          0.1671652586199110D+00, &
          0.1671652586199110D+00, &
          0.9974350214009486D+00, &
          0.8328347413800889D+00, &
          0.1671652586199110D+00, &
          0.0025649785990514D+00, &
          0.8328347413800889D+00, &
          0.8328347413800889D+00, &
          0.0025649785990514D+00, &
          0.1671652586199110D+00, &
          0.8328347413800889D+00, &
          0.0025649785990514D+00, &
          0.1671652586199110D+00, &
          0.1671652586199110D+00, &
          0.0025649785990514D+00, &
          0.8328347413800889D+00, &
          0.1671652586199110D+00, &
          0.2739297470225850D+00, &
          0.0974495447483840D+00, &
          0.0974495447483840D+00, &
          0.2739297470225850D+00, &
          0.9025504552516159D+00, &
          0.0974495447483840D+00, &
          0.2739297470225850D+00, &
          0.9025504552516159D+00, &
          0.9025504552516159D+00, &
          0.2739297470225850D+00, &
          0.0974495447483840D+00, &
          0.9025504552516159D+00, &
          0.7260702529774150D+00, &
          0.0974495447483840D+00, &
          0.0974495447483840D+00, &
          0.7260702529774150D+00, &
          0.9025504552516159D+00, &
          0.0974495447483840D+00, &
          0.7260702529774150D+00, &
          0.9025504552516159D+00, &
          0.9025504552516159D+00, &
          0.7260702529774150D+00, &
          0.0974495447483840D+00, &
          0.9025504552516159D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.1827950053146007D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.8172049946853993D+00, &
          0.9539270239814177D+00, &
          0.0460729760185823D+00, &
          0.9539270239814177D+00, &
          0.9539270239814177D+00, &
          0.0460729760185823D+00, &
          0.9539270239814177D+00, &
          0.0460729760185823D+00, &
          0.0460729760185823D+00, &
          0.3818973707741769D+00, &
          0.6181026292258232D+00, &
          0.3818973707741769D+00, &
          0.3818973707741769D+00, &
          0.6181026292258232D+00, &
          0.3818973707741769D+00, &
          0.6181026292258232D+00, &
          0.6181026292258232D+00, &
          0.5000000000000000D+00, &
          0.8687873171566183D+00, &
          0.8687873171566183D+00, &
          0.5000000000000000D+00, &
          0.8687873171566183D+00, &
          0.8687873171566183D+00, &
          0.5000000000000000D+00, &
          0.1312126828433817D+00, &
          0.1312126828433817D+00, &
          0.5000000000000000D+00, &
          0.1312126828433817D+00, &
          0.1312126828433817D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9666726015054183D+00, &
          0.0333273984945817D+00, &
          0.0333273984945817D+00, &
          0.9666726015054183D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9666726015054183D+00, &
          0.0333273984945817D+00, &
          0.0333273984945817D+00, &
          0.9666726015054183D+00, &
          0.7303580738482240D+00, &
          0.2696419261517760D+00, &
          0.7303580738482240D+00, &
          0.2696419261517760D+00, &
          0.7303580738482240D+00, &
          0.2696419261517760D+00, &
          0.7303580738482240D+00, &
          0.2696419261517760D+00, &
          0.6866920757899998D+00, &
          0.6866920757899998D+00, &
          0.8321081947965769D+00, &
          0.6866920757899998D+00, &
          0.6866920757899998D+00, &
          0.8321081947965769D+00, &
          0.3133079242100002D+00, &
          0.3133079242100002D+00, &
          0.8321081947965769D+00, &
          0.3133079242100002D+00, &
          0.3133079242100002D+00, &
          0.8321081947965769D+00, &
          0.6866920757899998D+00, &
          0.6866920757899998D+00, &
          0.1678918052034231D+00, &
          0.6866920757899998D+00, &
          0.6866920757899998D+00, &
          0.1678918052034231D+00, &
          0.3133079242100002D+00, &
          0.3133079242100002D+00, &
          0.1678918052034231D+00, &
          0.3133079242100002D+00, &
          0.3133079242100002D+00, &
          0.1678918052034231D+00, &
          0.9747555407321407D+00, &
          0.9747555407321407D+00, &
          0.6288381939284606D+00, &
          0.9747555407321407D+00, &
          0.9747555407321407D+00, &
          0.6288381939284606D+00, &
          0.0252444592678593D+00, &
          0.0252444592678593D+00, &
          0.6288381939284606D+00, &
          0.0252444592678593D+00, &
          0.0252444592678593D+00, &
          0.6288381939284606D+00, &
          0.9747555407321407D+00, &
          0.9747555407321407D+00, &
          0.3711618060715394D+00, &
          0.9747555407321407D+00, &
          0.9747555407321407D+00, &
          0.3711618060715394D+00, &
          0.0252444592678593D+00, &
          0.0252444592678593D+00, &
          0.3711618060715394D+00, &
          0.0252444592678593D+00, &
          0.0252444592678593D+00, &
          0.3711618060715394D+00, &
          0.8328347413800889D+00, &
          0.8328347413800889D+00, &
          0.9974350214009486D+00, &
          0.8328347413800889D+00, &
          0.8328347413800889D+00, &
          0.9974350214009486D+00, &
          0.1671652586199110D+00, &
          0.1671652586199110D+00, &
          0.9974350214009486D+00, &
          0.1671652586199110D+00, &
          0.1671652586199110D+00, &
          0.9974350214009486D+00, &
          0.8328347413800889D+00, &
          0.8328347413800889D+00, &
          0.0025649785990514D+00, &
          0.8328347413800889D+00, &
          0.8328347413800889D+00, &
          0.0025649785990514D+00, &
          0.1671652586199110D+00, &
          0.1671652586199110D+00, &
          0.0025649785990514D+00, &
          0.1671652586199110D+00, &
          0.1671652586199110D+00, &
          0.0025649785990514D+00, &
          0.0974495447483840D+00, &
          0.0974495447483840D+00, &
          0.2739297470225850D+00, &
          0.0974495447483840D+00, &
          0.0974495447483840D+00, &
          0.2739297470225850D+00, &
          0.9025504552516159D+00, &
          0.9025504552516159D+00, &
          0.2739297470225850D+00, &
          0.9025504552516159D+00, &
          0.9025504552516159D+00, &
          0.2739297470225850D+00, &
          0.0974495447483840D+00, &
          0.0974495447483840D+00, &
          0.7260702529774150D+00, &
          0.0974495447483840D+00, &
          0.0974495447483840D+00, &
          0.7260702529774150D+00, &
          0.9025504552516159D+00, &
          0.9025504552516159D+00, &
          0.7260702529774150D+00, &
          0.9025504552516159D+00, &
          0.9025504552516159D+00, &
          0.7260702529774150D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          0.0173912824874704D+00, &
          0.0173912824874704D+00, &
          0.0173912824874704D+00, &
          0.0173912824874704D+00, &
          0.0173912824874704D+00, &
          0.0173912824874704D+00, &
          0.0018386577800611D+00, &
          0.0018386577800611D+00, &
          0.0018386577800611D+00, &
          0.0018386577800611D+00, &
          0.0018386577800611D+00, &
          0.0018386577800611D+00, &
          0.0018386577800611D+00, &
          0.0018386577800611D+00, &
          0.0112910248314292D+00, &
          0.0112910248314292D+00, &
          0.0112910248314292D+00, &
          0.0112910248314292D+00, &
          0.0112910248314292D+00, &
          0.0112910248314292D+00, &
          0.0112910248314292D+00, &
          0.0112910248314292D+00, &
          0.0059757596352895D+00, &
          0.0059757596352895D+00, &
          0.0059757596352895D+00, &
          0.0059757596352895D+00, &
          0.0059757596352895D+00, &
          0.0059757596352895D+00, &
          0.0059757596352895D+00, &
          0.0059757596352895D+00, &
          0.0059757596352895D+00, &
          0.0059757596352895D+00, &
          0.0059757596352895D+00, &
          0.0059757596352895D+00, &
          0.0086912315244171D+00, &
          0.0086912315244171D+00, &
          0.0086912315244171D+00, &
          0.0086912315244171D+00, &
          0.0086912315244171D+00, &
          0.0086912315244171D+00, &
          0.0086912315244171D+00, &
          0.0086912315244171D+00, &
          0.0086912315244171D+00, &
          0.0086912315244171D+00, &
          0.0086912315244171D+00, &
          0.0086912315244171D+00, &
          0.0086912315244171D+00, &
          0.0086912315244171D+00, &
          0.0086912315244171D+00, &
          0.0086912315244171D+00, &
          0.0086912315244171D+00, &
          0.0086912315244171D+00, &
          0.0086912315244171D+00, &
          0.0086912315244171D+00, &
          0.0086912315244171D+00, &
          0.0086912315244171D+00, &
          0.0086912315244171D+00, &
          0.0086912315244171D+00, &
          0.0106145587103426D+00, &
          0.0106145587103426D+00, &
          0.0106145587103426D+00, &
          0.0106145587103426D+00, &
          0.0106145587103426D+00, &
          0.0106145587103426D+00, &
          0.0106145587103426D+00, &
          0.0106145587103426D+00, &
          0.0106145587103426D+00, &
          0.0106145587103426D+00, &
          0.0106145587103426D+00, &
          0.0106145587103426D+00, &
          0.0106145587103426D+00, &
          0.0106145587103426D+00, &
          0.0106145587103426D+00, &
          0.0106145587103426D+00, &
          0.0106145587103426D+00, &
          0.0106145587103426D+00, &
          0.0106145587103426D+00, &
          0.0106145587103426D+00, &
          0.0106145587103426D+00, &
          0.0106145587103426D+00, &
          0.0106145587103426D+00, &
          0.0106145587103426D+00, &
          0.0016903589166661D+00, &
          0.0016903589166661D+00, &
          0.0016903589166661D+00, &
          0.0016903589166661D+00, &
          0.0016903589166661D+00, &
          0.0016903589166661D+00, &
          0.0016903589166661D+00, &
          0.0016903589166661D+00, &
          0.0016903589166661D+00, &
          0.0016903589166661D+00, &
          0.0016903589166661D+00, &
          0.0016903589166661D+00, &
          0.0016903589166661D+00, &
          0.0016903589166661D+00, &
          0.0016903589166661D+00, &
          0.0016903589166661D+00, &
          0.0016903589166661D+00, &
          0.0016903589166661D+00, &
          0.0016903589166661D+00, &
          0.0016903589166661D+00, &
          0.0016903589166661D+00, &
          0.0016903589166661D+00, &
          0.0016903589166661D+00, &
          0.0016903589166661D+00, &
          0.0022842405528399D+00, &
          0.0022842405528399D+00, &
          0.0022842405528399D+00, &
          0.0022842405528399D+00, &
          0.0022842405528399D+00, &
          0.0022842405528399D+00, &
          0.0022842405528399D+00, &
          0.0022842405528399D+00, &
          0.0022842405528399D+00, &
          0.0022842405528399D+00, &
          0.0022842405528399D+00, &
          0.0022842405528399D+00, &
          0.0022842405528399D+00, &
          0.0022842405528399D+00, &
          0.0022842405528399D+00, &
          0.0022842405528399D+00, &
          0.0022842405528399D+00, &
          0.0022842405528399D+00, &
          0.0022842405528399D+00, &
          0.0022842405528399D+00, &
          0.0022842405528399D+00, &
          0.0022842405528399D+00, &
          0.0022842405528399D+00, &
          0.0022842405528399D+00, &
          0.0066740156523919D+00, &
          0.0066740156523919D+00, &
          0.0066740156523919D+00, &
          0.0066740156523919D+00, &
          0.0066740156523919D+00, &
          0.0066740156523919D+00, &
          0.0066740156523919D+00, &
          0.0066740156523919D+00, &
          0.0066740156523919D+00, &
          0.0066740156523919D+00, &
          0.0066740156523919D+00, &
          0.0066740156523919D+00, &
          0.0066740156523919D+00, &
          0.0066740156523919D+00, &
          0.0066740156523919D+00, &
          0.0066740156523919D+00, &
          0.0066740156523919D+00, &
          0.0066740156523919D+00, &
          0.0066740156523919D+00, &
          0.0066740156523919D+00, &
          0.0066740156523919D+00, &
          0.0066740156523919D+00, &
          0.0066740156523919D+00, &
          0.0066740156523919D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule15 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule15() returns the hexahedron quadrature rule of precision 15.
!
!  Discussion:
!
!    We suppose we are given a hexahedron H with vertices
!    (0,0,0), (1,0,0), (1,1,0), (0,1,0), (0,0,1), (1,0,1), (1,1,1), (0,1,1),
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over H is 
!    approximated by Q as follows:
!
!    Q = volume(H) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 May 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order symmetric cubature rules for tetrahedra and pyramids,
!    International Journal for Numerical Methods in Engineering,
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
!    real w(n): the quadrature weights, which add to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 256

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
          0.9619577574752056D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.0380422425247944D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.3913025693166751D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.6086974306833248D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.6692365934747126D+00, &
          0.3307634065252874D+00, &
          0.3307634065252874D+00, &
          0.6692365934747126D+00, &
          0.6692365934747126D+00, &
          0.3307634065252874D+00, &
          0.6692365934747126D+00, &
          0.3307634065252874D+00, &
          0.8965445655196689D+00, &
          0.1034554344803311D+00, &
          0.1034554344803311D+00, &
          0.8965445655196689D+00, &
          0.8965445655196689D+00, &
          0.1034554344803311D+00, &
          0.8965445655196689D+00, &
          0.1034554344803311D+00, &
          0.8482759513633358D+00, &
          0.5000000000000000D+00, &
          0.8482759513633358D+00, &
          0.1517240486366642D+00, &
          0.5000000000000000D+00, &
          0.1517240486366642D+00, &
          0.8482759513633358D+00, &
          0.5000000000000000D+00, &
          0.8482759513633358D+00, &
          0.1517240486366642D+00, &
          0.5000000000000000D+00, &
          0.1517240486366642D+00, &
          0.6321333796314236D+00, &
          0.3678666203685764D+00, &
          0.6321333796314236D+00, &
          0.3678666203685764D+00, &
          0.6321333796314236D+00, &
          0.3678666203685764D+00, &
          0.6321333796314236D+00, &
          0.3678666203685764D+00, &
          0.7924949743436622D+00, &
          0.2075050256563378D+00, &
          0.2075050256563378D+00, &
          0.7924949743436622D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.7924949743436622D+00, &
          0.2075050256563378D+00, &
          0.2075050256563378D+00, &
          0.7924949743436622D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.7764752956685337D+00, &
          0.2235247043314663D+00, &
          0.7764752956685337D+00, &
          0.2235247043314663D+00, &
          0.7764752956685337D+00, &
          0.2235247043314663D+00, &
          0.7764752956685337D+00, &
          0.2235247043314663D+00, &
          0.9914024363900886D+00, &
          0.0085975636099114D+00, &
          0.0085975636099114D+00, &
          0.9914024363900886D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9914024363900886D+00, &
          0.0085975636099114D+00, &
          0.0085975636099114D+00, &
          0.9914024363900886D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.6360162945430453D+00, &
          0.9390447038366581D+00, &
          0.6360162945430453D+00, &
          0.3639837054569547D+00, &
          0.9390447038366581D+00, &
          0.3639837054569547D+00, &
          0.3639837054569547D+00, &
          0.9390447038366581D+00, &
          0.3639837054569547D+00, &
          0.6360162945430453D+00, &
          0.9390447038366581D+00, &
          0.6360162945430453D+00, &
          0.6360162945430453D+00, &
          0.0609552961633420D+00, &
          0.6360162945430453D+00, &
          0.3639837054569547D+00, &
          0.0609552961633420D+00, &
          0.3639837054569547D+00, &
          0.3639837054569547D+00, &
          0.0609552961633420D+00, &
          0.3639837054569547D+00, &
          0.6360162945430453D+00, &
          0.0609552961633420D+00, &
          0.6360162945430453D+00, &
          0.8093239300427761D+00, &
          0.9699906676537717D+00, &
          0.8093239300427761D+00, &
          0.1906760699572239D+00, &
          0.9699906676537717D+00, &
          0.1906760699572239D+00, &
          0.1906760699572239D+00, &
          0.9699906676537717D+00, &
          0.1906760699572239D+00, &
          0.8093239300427761D+00, &
          0.9699906676537717D+00, &
          0.8093239300427761D+00, &
          0.8093239300427761D+00, &
          0.0300093323462282D+00, &
          0.8093239300427761D+00, &
          0.1906760699572239D+00, &
          0.0300093323462282D+00, &
          0.1906760699572239D+00, &
          0.1906760699572239D+00, &
          0.0300093323462282D+00, &
          0.1906760699572239D+00, &
          0.8093239300427761D+00, &
          0.0300093323462282D+00, &
          0.8093239300427761D+00, &
          0.9315772312765189D+00, &
          0.6521112975876302D+00, &
          0.9315772312765189D+00, &
          0.0684227687234810D+00, &
          0.6521112975876302D+00, &
          0.0684227687234810D+00, &
          0.0684227687234810D+00, &
          0.6521112975876302D+00, &
          0.0684227687234810D+00, &
          0.9315772312765189D+00, &
          0.6521112975876302D+00, &
          0.9315772312765189D+00, &
          0.9315772312765189D+00, &
          0.3478887024123698D+00, &
          0.9315772312765189D+00, &
          0.0684227687234810D+00, &
          0.3478887024123698D+00, &
          0.0684227687234810D+00, &
          0.0684227687234810D+00, &
          0.3478887024123698D+00, &
          0.0684227687234810D+00, &
          0.9315772312765189D+00, &
          0.3478887024123698D+00, &
          0.9315772312765189D+00, &
          0.7308028863234025D+00, &
          0.8685680591748091D+00, &
          0.7308028863234025D+00, &
          0.2691971136765975D+00, &
          0.8685680591748091D+00, &
          0.2691971136765975D+00, &
          0.2691971136765975D+00, &
          0.8685680591748091D+00, &
          0.2691971136765975D+00, &
          0.7308028863234025D+00, &
          0.8685680591748091D+00, &
          0.7308028863234025D+00, &
          0.7308028863234025D+00, &
          0.1314319408251909D+00, &
          0.7308028863234025D+00, &
          0.2691971136765975D+00, &
          0.1314319408251909D+00, &
          0.2691971136765975D+00, &
          0.2691971136765975D+00, &
          0.1314319408251909D+00, &
          0.2691971136765975D+00, &
          0.7308028863234025D+00, &
          0.1314319408251909D+00, &
          0.7308028863234025D+00, &
          0.9801570902708128D+00, &
          0.8835073040170505D+00, &
          0.9801570902708128D+00, &
          0.0198429097291872D+00, &
          0.8835073040170505D+00, &
          0.0198429097291872D+00, &
          0.0198429097291872D+00, &
          0.8835073040170505D+00, &
          0.0198429097291872D+00, &
          0.9801570902708128D+00, &
          0.8835073040170505D+00, &
          0.9801570902708128D+00, &
          0.9801570902708128D+00, &
          0.1164926959829496D+00, &
          0.9801570902708128D+00, &
          0.0198429097291872D+00, &
          0.1164926959829496D+00, &
          0.0198429097291872D+00, &
          0.0198429097291872D+00, &
          0.1164926959829496D+00, &
          0.0198429097291872D+00, &
          0.9801570902708128D+00, &
          0.1164926959829496D+00, &
          0.9801570902708128D+00, &
          0.6401259198958773D+00, &
          0.6401259198958773D+00, &
          0.9442969164789898D+00, &
          0.9971300972924322D+00, &
          0.9442969164789898D+00, &
          0.9971300972924322D+00, &
          0.3598740801041227D+00, &
          0.3598740801041227D+00, &
          0.9442969164789898D+00, &
          0.9971300972924322D+00, &
          0.9442969164789898D+00, &
          0.9971300972924322D+00, &
          0.6401259198958773D+00, &
          0.6401259198958773D+00, &
          0.0557030835210102D+00, &
          0.9971300972924322D+00, &
          0.0557030835210102D+00, &
          0.9971300972924322D+00, &
          0.6401259198958773D+00, &
          0.6401259198958773D+00, &
          0.9442969164789898D+00, &
          0.0028699027075679D+00, &
          0.9442969164789898D+00, &
          0.0028699027075679D+00, &
          0.3598740801041227D+00, &
          0.3598740801041227D+00, &
          0.0557030835210102D+00, &
          0.9971300972924322D+00, &
          0.0557030835210102D+00, &
          0.9971300972924322D+00, &
          0.3598740801041227D+00, &
          0.3598740801041227D+00, &
          0.9442969164789898D+00, &
          0.0028699027075679D+00, &
          0.9442969164789898D+00, &
          0.0028699027075679D+00, &
          0.6401259198958773D+00, &
          0.6401259198958773D+00, &
          0.0557030835210102D+00, &
          0.0028699027075679D+00, &
          0.0557030835210102D+00, &
          0.0028699027075679D+00, &
          0.3598740801041227D+00, &
          0.3598740801041227D+00, &
          0.0557030835210102D+00, &
          0.0028699027075679D+00, &
          0.0557030835210102D+00, &
          0.0028699027075679D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.5000000000000000D+00, &
          0.9619577574752056D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.0380422425247944D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.3913025693166751D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.6086974306833248D+00, &
          0.5000000000000000D+00, &
          0.6692365934747126D+00, &
          0.3307634065252874D+00, &
          0.6692365934747126D+00, &
          0.3307634065252874D+00, &
          0.6692365934747126D+00, &
          0.3307634065252874D+00, &
          0.3307634065252874D+00, &
          0.6692365934747126D+00, &
          0.8965445655196689D+00, &
          0.1034554344803311D+00, &
          0.8965445655196689D+00, &
          0.1034554344803311D+00, &
          0.8965445655196689D+00, &
          0.1034554344803311D+00, &
          0.1034554344803311D+00, &
          0.8965445655196689D+00, &
          0.8482759513633358D+00, &
          0.8482759513633358D+00, &
          0.5000000000000000D+00, &
          0.8482759513633358D+00, &
          0.1517240486366642D+00, &
          0.5000000000000000D+00, &
          0.1517240486366642D+00, &
          0.8482759513633358D+00, &
          0.5000000000000000D+00, &
          0.1517240486366642D+00, &
          0.1517240486366642D+00, &
          0.5000000000000000D+00, &
          0.7924949743436622D+00, &
          0.2075050256563378D+00, &
          0.2075050256563378D+00, &
          0.7924949743436622D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.6321333796314236D+00, &
          0.3678666203685764D+00, &
          0.6321333796314236D+00, &
          0.3678666203685764D+00, &
          0.6321333796314236D+00, &
          0.3678666203685764D+00, &
          0.6321333796314236D+00, &
          0.3678666203685764D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.7924949743436622D+00, &
          0.2075050256563378D+00, &
          0.2075050256563378D+00, &
          0.7924949743436622D+00, &
          0.9914024363900886D+00, &
          0.0085975636099114D+00, &
          0.0085975636099114D+00, &
          0.9914024363900886D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.7764752956685337D+00, &
          0.2235247043314663D+00, &
          0.7764752956685337D+00, &
          0.2235247043314663D+00, &
          0.7764752956685337D+00, &
          0.2235247043314663D+00, &
          0.7764752956685337D+00, &
          0.2235247043314663D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9914024363900886D+00, &
          0.0085975636099114D+00, &
          0.0085975636099114D+00, &
          0.9914024363900886D+00, &
          0.9390447038366581D+00, &
          0.6360162945430453D+00, &
          0.6360162945430453D+00, &
          0.9390447038366581D+00, &
          0.3639837054569547D+00, &
          0.6360162945430453D+00, &
          0.9390447038366581D+00, &
          0.3639837054569547D+00, &
          0.3639837054569547D+00, &
          0.9390447038366581D+00, &
          0.6360162945430453D+00, &
          0.3639837054569547D+00, &
          0.0609552961633420D+00, &
          0.6360162945430453D+00, &
          0.6360162945430453D+00, &
          0.0609552961633420D+00, &
          0.3639837054569547D+00, &
          0.6360162945430453D+00, &
          0.0609552961633420D+00, &
          0.3639837054569547D+00, &
          0.3639837054569547D+00, &
          0.0609552961633420D+00, &
          0.6360162945430453D+00, &
          0.3639837054569547D+00, &
          0.9699906676537717D+00, &
          0.8093239300427761D+00, &
          0.8093239300427761D+00, &
          0.9699906676537717D+00, &
          0.1906760699572239D+00, &
          0.8093239300427761D+00, &
          0.9699906676537717D+00, &
          0.1906760699572239D+00, &
          0.1906760699572239D+00, &
          0.9699906676537717D+00, &
          0.8093239300427761D+00, &
          0.1906760699572239D+00, &
          0.0300093323462282D+00, &
          0.8093239300427761D+00, &
          0.8093239300427761D+00, &
          0.0300093323462282D+00, &
          0.1906760699572239D+00, &
          0.8093239300427761D+00, &
          0.0300093323462282D+00, &
          0.1906760699572239D+00, &
          0.1906760699572239D+00, &
          0.0300093323462282D+00, &
          0.8093239300427761D+00, &
          0.1906760699572239D+00, &
          0.6521112975876302D+00, &
          0.9315772312765189D+00, &
          0.9315772312765189D+00, &
          0.6521112975876302D+00, &
          0.0684227687234810D+00, &
          0.9315772312765189D+00, &
          0.6521112975876302D+00, &
          0.0684227687234810D+00, &
          0.0684227687234810D+00, &
          0.6521112975876302D+00, &
          0.9315772312765189D+00, &
          0.0684227687234810D+00, &
          0.3478887024123698D+00, &
          0.9315772312765189D+00, &
          0.9315772312765189D+00, &
          0.3478887024123698D+00, &
          0.0684227687234810D+00, &
          0.9315772312765189D+00, &
          0.3478887024123698D+00, &
          0.0684227687234810D+00, &
          0.0684227687234810D+00, &
          0.3478887024123698D+00, &
          0.9315772312765189D+00, &
          0.0684227687234810D+00, &
          0.8685680591748091D+00, &
          0.7308028863234025D+00, &
          0.7308028863234025D+00, &
          0.8685680591748091D+00, &
          0.2691971136765975D+00, &
          0.7308028863234025D+00, &
          0.8685680591748091D+00, &
          0.2691971136765975D+00, &
          0.2691971136765975D+00, &
          0.8685680591748091D+00, &
          0.7308028863234025D+00, &
          0.2691971136765975D+00, &
          0.1314319408251909D+00, &
          0.7308028863234025D+00, &
          0.7308028863234025D+00, &
          0.1314319408251909D+00, &
          0.2691971136765975D+00, &
          0.7308028863234025D+00, &
          0.1314319408251909D+00, &
          0.2691971136765975D+00, &
          0.2691971136765975D+00, &
          0.1314319408251909D+00, &
          0.7308028863234025D+00, &
          0.2691971136765975D+00, &
          0.8835073040170505D+00, &
          0.9801570902708128D+00, &
          0.9801570902708128D+00, &
          0.8835073040170505D+00, &
          0.0198429097291872D+00, &
          0.9801570902708128D+00, &
          0.8835073040170505D+00, &
          0.0198429097291872D+00, &
          0.0198429097291872D+00, &
          0.8835073040170505D+00, &
          0.9801570902708128D+00, &
          0.0198429097291872D+00, &
          0.1164926959829496D+00, &
          0.9801570902708128D+00, &
          0.9801570902708128D+00, &
          0.1164926959829496D+00, &
          0.0198429097291872D+00, &
          0.9801570902708128D+00, &
          0.1164926959829496D+00, &
          0.0198429097291872D+00, &
          0.0198429097291872D+00, &
          0.1164926959829496D+00, &
          0.9801570902708128D+00, &
          0.0198429097291872D+00, &
          0.9442969164789898D+00, &
          0.9971300972924322D+00, &
          0.6401259198958773D+00, &
          0.6401259198958773D+00, &
          0.9971300972924322D+00, &
          0.9442969164789898D+00, &
          0.9442969164789898D+00, &
          0.9971300972924322D+00, &
          0.3598740801041227D+00, &
          0.3598740801041227D+00, &
          0.9971300972924322D+00, &
          0.9442969164789898D+00, &
          0.0557030835210102D+00, &
          0.9971300972924322D+00, &
          0.6401259198958773D+00, &
          0.6401259198958773D+00, &
          0.9971300972924322D+00, &
          0.0557030835210102D+00, &
          0.9442969164789898D+00, &
          0.0028699027075679D+00, &
          0.6401259198958773D+00, &
          0.6401259198958773D+00, &
          0.0028699027075679D+00, &
          0.9442969164789898D+00, &
          0.0557030835210102D+00, &
          0.9971300972924322D+00, &
          0.3598740801041227D+00, &
          0.3598740801041227D+00, &
          0.9971300972924322D+00, &
          0.0557030835210102D+00, &
          0.9442969164789898D+00, &
          0.0028699027075679D+00, &
          0.3598740801041227D+00, &
          0.3598740801041227D+00, &
          0.0028699027075679D+00, &
          0.9442969164789898D+00, &
          0.0557030835210102D+00, &
          0.0028699027075679D+00, &
          0.6401259198958773D+00, &
          0.6401259198958773D+00, &
          0.0028699027075679D+00, &
          0.0557030835210102D+00, &
          0.0557030835210102D+00, &
          0.0028699027075679D+00, &
          0.3598740801041227D+00, &
          0.3598740801041227D+00, &
          0.0028699027075679D+00, &
          0.0557030835210102D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9619577574752056D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.0380422425247944D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.3913025693166751D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.6086974306833248D+00, &
          0.6692365934747126D+00, &
          0.3307634065252874D+00, &
          0.6692365934747126D+00, &
          0.6692365934747126D+00, &
          0.3307634065252874D+00, &
          0.6692365934747126D+00, &
          0.3307634065252874D+00, &
          0.3307634065252874D+00, &
          0.8965445655196689D+00, &
          0.1034554344803311D+00, &
          0.8965445655196689D+00, &
          0.8965445655196689D+00, &
          0.1034554344803311D+00, &
          0.8965445655196689D+00, &
          0.1034554344803311D+00, &
          0.1034554344803311D+00, &
          0.5000000000000000D+00, &
          0.8482759513633358D+00, &
          0.8482759513633358D+00, &
          0.5000000000000000D+00, &
          0.8482759513633358D+00, &
          0.8482759513633358D+00, &
          0.5000000000000000D+00, &
          0.1517240486366642D+00, &
          0.1517240486366642D+00, &
          0.5000000000000000D+00, &
          0.1517240486366642D+00, &
          0.1517240486366642D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.7924949743436622D+00, &
          0.2075050256563378D+00, &
          0.2075050256563378D+00, &
          0.7924949743436622D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.7924949743436622D+00, &
          0.2075050256563378D+00, &
          0.2075050256563378D+00, &
          0.7924949743436622D+00, &
          0.6321333796314236D+00, &
          0.3678666203685764D+00, &
          0.6321333796314236D+00, &
          0.3678666203685764D+00, &
          0.6321333796314236D+00, &
          0.3678666203685764D+00, &
          0.6321333796314236D+00, &
          0.3678666203685764D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9914024363900886D+00, &
          0.0085975636099114D+00, &
          0.0085975636099114D+00, &
          0.9914024363900886D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9914024363900886D+00, &
          0.0085975636099114D+00, &
          0.0085975636099114D+00, &
          0.9914024363900886D+00, &
          0.7764752956685337D+00, &
          0.2235247043314663D+00, &
          0.7764752956685337D+00, &
          0.2235247043314663D+00, &
          0.7764752956685337D+00, &
          0.2235247043314663D+00, &
          0.7764752956685337D+00, &
          0.2235247043314663D+00, &
          0.6360162945430453D+00, &
          0.6360162945430453D+00, &
          0.9390447038366581D+00, &
          0.6360162945430453D+00, &
          0.6360162945430453D+00, &
          0.9390447038366581D+00, &
          0.3639837054569547D+00, &
          0.3639837054569547D+00, &
          0.9390447038366581D+00, &
          0.3639837054569547D+00, &
          0.3639837054569547D+00, &
          0.9390447038366581D+00, &
          0.6360162945430453D+00, &
          0.6360162945430453D+00, &
          0.0609552961633420D+00, &
          0.6360162945430453D+00, &
          0.6360162945430453D+00, &
          0.0609552961633420D+00, &
          0.3639837054569547D+00, &
          0.3639837054569547D+00, &
          0.0609552961633420D+00, &
          0.3639837054569547D+00, &
          0.3639837054569547D+00, &
          0.0609552961633420D+00, &
          0.8093239300427761D+00, &
          0.8093239300427761D+00, &
          0.9699906676537717D+00, &
          0.8093239300427761D+00, &
          0.8093239300427761D+00, &
          0.9699906676537717D+00, &
          0.1906760699572239D+00, &
          0.1906760699572239D+00, &
          0.9699906676537717D+00, &
          0.1906760699572239D+00, &
          0.1906760699572239D+00, &
          0.9699906676537717D+00, &
          0.8093239300427761D+00, &
          0.8093239300427761D+00, &
          0.0300093323462282D+00, &
          0.8093239300427761D+00, &
          0.8093239300427761D+00, &
          0.0300093323462282D+00, &
          0.1906760699572239D+00, &
          0.1906760699572239D+00, &
          0.0300093323462282D+00, &
          0.1906760699572239D+00, &
          0.1906760699572239D+00, &
          0.0300093323462282D+00, &
          0.9315772312765189D+00, &
          0.9315772312765189D+00, &
          0.6521112975876302D+00, &
          0.9315772312765189D+00, &
          0.9315772312765189D+00, &
          0.6521112975876302D+00, &
          0.0684227687234810D+00, &
          0.0684227687234810D+00, &
          0.6521112975876302D+00, &
          0.0684227687234810D+00, &
          0.0684227687234810D+00, &
          0.6521112975876302D+00, &
          0.9315772312765189D+00, &
          0.9315772312765189D+00, &
          0.3478887024123698D+00, &
          0.9315772312765189D+00, &
          0.9315772312765189D+00, &
          0.3478887024123698D+00, &
          0.0684227687234810D+00, &
          0.0684227687234810D+00, &
          0.3478887024123698D+00, &
          0.0684227687234810D+00, &
          0.0684227687234810D+00, &
          0.3478887024123698D+00, &
          0.7308028863234025D+00, &
          0.7308028863234025D+00, &
          0.8685680591748091D+00, &
          0.7308028863234025D+00, &
          0.7308028863234025D+00, &
          0.8685680591748091D+00, &
          0.2691971136765975D+00, &
          0.2691971136765975D+00, &
          0.8685680591748091D+00, &
          0.2691971136765975D+00, &
          0.2691971136765975D+00, &
          0.8685680591748091D+00, &
          0.7308028863234025D+00, &
          0.7308028863234025D+00, &
          0.1314319408251909D+00, &
          0.7308028863234025D+00, &
          0.7308028863234025D+00, &
          0.1314319408251909D+00, &
          0.2691971136765975D+00, &
          0.2691971136765975D+00, &
          0.1314319408251909D+00, &
          0.2691971136765975D+00, &
          0.2691971136765975D+00, &
          0.1314319408251909D+00, &
          0.9801570902708128D+00, &
          0.9801570902708128D+00, &
          0.8835073040170505D+00, &
          0.9801570902708128D+00, &
          0.9801570902708128D+00, &
          0.8835073040170505D+00, &
          0.0198429097291872D+00, &
          0.0198429097291872D+00, &
          0.8835073040170505D+00, &
          0.0198429097291872D+00, &
          0.0198429097291872D+00, &
          0.8835073040170505D+00, &
          0.9801570902708128D+00, &
          0.9801570902708128D+00, &
          0.1164926959829496D+00, &
          0.9801570902708128D+00, &
          0.9801570902708128D+00, &
          0.1164926959829496D+00, &
          0.0198429097291872D+00, &
          0.0198429097291872D+00, &
          0.1164926959829496D+00, &
          0.0198429097291872D+00, &
          0.0198429097291872D+00, &
          0.1164926959829496D+00, &
          0.9971300972924322D+00, &
          0.9442969164789898D+00, &
          0.9971300972924322D+00, &
          0.9442969164789898D+00, &
          0.6401259198958773D+00, &
          0.6401259198958773D+00, &
          0.9971300972924322D+00, &
          0.9442969164789898D+00, &
          0.9971300972924322D+00, &
          0.9442969164789898D+00, &
          0.3598740801041227D+00, &
          0.3598740801041227D+00, &
          0.9971300972924322D+00, &
          0.0557030835210102D+00, &
          0.9971300972924322D+00, &
          0.0557030835210102D+00, &
          0.6401259198958773D+00, &
          0.6401259198958773D+00, &
          0.0028699027075679D+00, &
          0.9442969164789898D+00, &
          0.0028699027075679D+00, &
          0.9442969164789898D+00, &
          0.6401259198958773D+00, &
          0.6401259198958773D+00, &
          0.9971300972924322D+00, &
          0.0557030835210102D+00, &
          0.9971300972924322D+00, &
          0.0557030835210102D+00, &
          0.3598740801041227D+00, &
          0.3598740801041227D+00, &
          0.0028699027075679D+00, &
          0.9442969164789898D+00, &
          0.0028699027075679D+00, &
          0.9442969164789898D+00, &
          0.3598740801041227D+00, &
          0.3598740801041227D+00, &
          0.0028699027075679D+00, &
          0.0557030835210102D+00, &
          0.0028699027075679D+00, &
          0.0557030835210102D+00, &
          0.6401259198958773D+00, &
          0.6401259198958773D+00, &
          0.0028699027075679D+00, &
          0.0557030835210102D+00, &
          0.0028699027075679D+00, &
          0.0557030835210102D+00, &
          0.3598740801041227D+00, &
          0.3598740801041227D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          0.0033403050106945D+00, &
          0.0033403050106945D+00, &
          0.0033403050106945D+00, &
          0.0033403050106945D+00, &
          0.0033403050106945D+00, &
          0.0033403050106945D+00, &
          0.0050185523378773D+00, &
          0.0050185523378773D+00, &
          0.0050185523378773D+00, &
          0.0050185523378773D+00, &
          0.0050185523378773D+00, &
          0.0050185523378773D+00, &
          0.0090827588172248D+00, &
          0.0090827588172248D+00, &
          0.0090827588172248D+00, &
          0.0090827588172248D+00, &
          0.0090827588172248D+00, &
          0.0090827588172248D+00, &
          0.0090827588172248D+00, &
          0.0090827588172248D+00, &
          0.0029364087383876D+00, &
          0.0029364087383876D+00, &
          0.0029364087383876D+00, &
          0.0029364087383876D+00, &
          0.0029364087383876D+00, &
          0.0029364087383876D+00, &
          0.0029364087383876D+00, &
          0.0029364087383876D+00, &
          0.0093182386470488D+00, &
          0.0093182386470488D+00, &
          0.0093182386470488D+00, &
          0.0093182386470488D+00, &
          0.0093182386470488D+00, &
          0.0093182386470488D+00, &
          0.0093182386470488D+00, &
          0.0093182386470488D+00, &
          0.0093182386470488D+00, &
          0.0093182386470488D+00, &
          0.0093182386470488D+00, &
          0.0093182386470488D+00, &
          0.0069907862357511D+00, &
          0.0069907862357511D+00, &
          0.0069907862357511D+00, &
          0.0069907862357511D+00, &
          0.0069907862357511D+00, &
          0.0069907862357511D+00, &
          0.0069907862357511D+00, &
          0.0069907862357511D+00, &
          0.0069907862357511D+00, &
          0.0069907862357511D+00, &
          0.0069907862357511D+00, &
          0.0069907862357511D+00, &
          0.0069907862357511D+00, &
          0.0069907862357511D+00, &
          0.0069907862357511D+00, &
          0.0069907862357511D+00, &
          0.0069907862357511D+00, &
          0.0069907862357511D+00, &
          0.0069907862357511D+00, &
          0.0069907862357511D+00, &
          0.0069907862357511D+00, &
          0.0069907862357511D+00, &
          0.0069907862357511D+00, &
          0.0069907862357511D+00, &
          0.0028410057474906D+00, &
          0.0028410057474906D+00, &
          0.0028410057474906D+00, &
          0.0028410057474906D+00, &
          0.0028410057474906D+00, &
          0.0028410057474906D+00, &
          0.0028410057474906D+00, &
          0.0028410057474906D+00, &
          0.0028410057474906D+00, &
          0.0028410057474906D+00, &
          0.0028410057474906D+00, &
          0.0028410057474906D+00, &
          0.0028410057474906D+00, &
          0.0028410057474906D+00, &
          0.0028410057474906D+00, &
          0.0028410057474906D+00, &
          0.0028410057474906D+00, &
          0.0028410057474906D+00, &
          0.0028410057474906D+00, &
          0.0028410057474906D+00, &
          0.0028410057474906D+00, &
          0.0028410057474906D+00, &
          0.0028410057474906D+00, &
          0.0028410057474906D+00, &
          0.0053095740723308D+00, &
          0.0053095740723308D+00, &
          0.0053095740723308D+00, &
          0.0053095740723308D+00, &
          0.0053095740723308D+00, &
          0.0053095740723308D+00, &
          0.0053095740723308D+00, &
          0.0053095740723308D+00, &
          0.0053095740723308D+00, &
          0.0053095740723308D+00, &
          0.0053095740723308D+00, &
          0.0053095740723308D+00, &
          0.0053095740723308D+00, &
          0.0053095740723308D+00, &
          0.0053095740723308D+00, &
          0.0053095740723308D+00, &
          0.0053095740723308D+00, &
          0.0053095740723308D+00, &
          0.0053095740723308D+00, &
          0.0053095740723308D+00, &
          0.0053095740723308D+00, &
          0.0053095740723308D+00, &
          0.0053095740723308D+00, &
          0.0053095740723308D+00, &
          0.0035708945365673D+00, &
          0.0035708945365673D+00, &
          0.0035708945365673D+00, &
          0.0035708945365673D+00, &
          0.0035708945365673D+00, &
          0.0035708945365673D+00, &
          0.0035708945365673D+00, &
          0.0035708945365673D+00, &
          0.0035708945365673D+00, &
          0.0035708945365673D+00, &
          0.0035708945365673D+00, &
          0.0035708945365673D+00, &
          0.0035708945365673D+00, &
          0.0035708945365673D+00, &
          0.0035708945365673D+00, &
          0.0035708945365673D+00, &
          0.0035708945365673D+00, &
          0.0035708945365673D+00, &
          0.0035708945365673D+00, &
          0.0035708945365673D+00, &
          0.0035708945365673D+00, &
          0.0035708945365673D+00, &
          0.0035708945365673D+00, &
          0.0035708945365673D+00, &
          0.0035058721795857D+00, &
          0.0035058721795857D+00, &
          0.0035058721795857D+00, &
          0.0035058721795857D+00, &
          0.0035058721795857D+00, &
          0.0035058721795857D+00, &
          0.0035058721795857D+00, &
          0.0035058721795857D+00, &
          0.0035058721795857D+00, &
          0.0035058721795857D+00, &
          0.0035058721795857D+00, &
          0.0035058721795857D+00, &
          0.0035058721795857D+00, &
          0.0035058721795857D+00, &
          0.0035058721795857D+00, &
          0.0035058721795857D+00, &
          0.0035058721795857D+00, &
          0.0035058721795857D+00, &
          0.0035058721795857D+00, &
          0.0035058721795857D+00, &
          0.0035058721795857D+00, &
          0.0035058721795857D+00, &
          0.0035058721795857D+00, &
          0.0035058721795857D+00, &
          0.0066651990511383D+00, &
          0.0066651990511383D+00, &
          0.0066651990511383D+00, &
          0.0066651990511383D+00, &
          0.0066651990511383D+00, &
          0.0066651990511383D+00, &
          0.0066651990511383D+00, &
          0.0066651990511383D+00, &
          0.0066651990511383D+00, &
          0.0066651990511383D+00, &
          0.0066651990511383D+00, &
          0.0066651990511383D+00, &
          0.0066651990511383D+00, &
          0.0066651990511383D+00, &
          0.0066651990511383D+00, &
          0.0066651990511383D+00, &
          0.0066651990511383D+00, &
          0.0066651990511383D+00, &
          0.0066651990511383D+00, &
          0.0066651990511383D+00, &
          0.0066651990511383D+00, &
          0.0066651990511383D+00, &
          0.0066651990511383D+00, &
          0.0066651990511383D+00, &
          0.0007781518386120D+00, &
          0.0007781518386120D+00, &
          0.0007781518386120D+00, &
          0.0007781518386120D+00, &
          0.0007781518386120D+00, &
          0.0007781518386120D+00, &
          0.0007781518386120D+00, &
          0.0007781518386120D+00, &
          0.0007781518386120D+00, &
          0.0007781518386120D+00, &
          0.0007781518386120D+00, &
          0.0007781518386120D+00, &
          0.0007781518386120D+00, &
          0.0007781518386120D+00, &
          0.0007781518386120D+00, &
          0.0007781518386120D+00, &
          0.0007781518386120D+00, &
          0.0007781518386120D+00, &
          0.0007781518386120D+00, &
          0.0007781518386120D+00, &
          0.0007781518386120D+00, &
          0.0007781518386120D+00, &
          0.0007781518386120D+00, &
          0.0007781518386120D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00, &
          0.0006249800796597D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule17 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule17() returns the hexahedron quadrature rule of precision 17.
!
!  Discussion:
!
!    We suppose we are given a hexahedron H with vertices
!    (0,0,0), (1,0,0), (1,1,0), (0,1,0), (0,0,1), (1,0,1), (1,1,1), (0,1,1),
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over H is 
!    approximated by Q as follows:
!
!    Q = volume(H) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 May 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order symmetric cubature rules for tetrahedra and pyramids,
!    International Journal for Numerical Methods in Engineering,
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
!    real w(n): the quadrature weights, which add to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 346

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
          0.7426575151151338D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.2573424848848662D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.8703394079370492D+00, &
          0.1296605920629508D+00, &
          0.1296605920629508D+00, &
          0.8703394079370492D+00, &
          0.8703394079370492D+00, &
          0.1296605920629508D+00, &
          0.8703394079370492D+00, &
          0.1296605920629508D+00, &
          0.9999997917911259D+00, &
          0.0000002082088741D+00, &
          0.0000002082088741D+00, &
          0.9999997917911259D+00, &
          0.9999997917911259D+00, &
          0.0000002082088741D+00, &
          0.9999997917911259D+00, &
          0.0000002082088741D+00, &
          0.8773134026065881D+00, &
          0.5000000000000000D+00, &
          0.8773134026065881D+00, &
          0.1226865973934119D+00, &
          0.5000000000000000D+00, &
          0.1226865973934119D+00, &
          0.8773134026065881D+00, &
          0.5000000000000000D+00, &
          0.8773134026065881D+00, &
          0.1226865973934119D+00, &
          0.5000000000000000D+00, &
          0.1226865973934119D+00, &
          0.7748409070843456D+00, &
          0.5000000000000000D+00, &
          0.7748409070843456D+00, &
          0.2251590929156544D+00, &
          0.5000000000000000D+00, &
          0.2251590929156544D+00, &
          0.7748409070843456D+00, &
          0.5000000000000000D+00, &
          0.7748409070843456D+00, &
          0.2251590929156544D+00, &
          0.5000000000000000D+00, &
          0.2251590929156544D+00, &
          0.9998342518711405D+00, &
          0.5000000000000000D+00, &
          0.9998342518711405D+00, &
          0.0001657481288594D+00, &
          0.5000000000000000D+00, &
          0.0001657481288594D+00, &
          0.9998342518711405D+00, &
          0.5000000000000000D+00, &
          0.9998342518711405D+00, &
          0.0001657481288594D+00, &
          0.5000000000000000D+00, &
          0.0001657481288594D+00, &
          0.8078645495842367D+00, &
          0.1921354504157633D+00, &
          0.8078645495842367D+00, &
          0.1921354504157633D+00, &
          0.8078645495842367D+00, &
          0.1921354504157633D+00, &
          0.8078645495842367D+00, &
          0.1921354504157633D+00, &
          0.9765209703524975D+00, &
          0.0234790296475025D+00, &
          0.0234790296475025D+00, &
          0.9765209703524975D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9765209703524975D+00, &
          0.0234790296475025D+00, &
          0.0234790296475025D+00, &
          0.9765209703524975D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.6378406365996775D+00, &
          0.3621593634003225D+00, &
          0.6378406365996775D+00, &
          0.3621593634003225D+00, &
          0.6378406365996775D+00, &
          0.3621593634003225D+00, &
          0.6378406365996775D+00, &
          0.3621593634003225D+00, &
          0.9085620340535459D+00, &
          0.0914379659464542D+00, &
          0.0914379659464542D+00, &
          0.9085620340535459D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9085620340535459D+00, &
          0.0914379659464542D+00, &
          0.0914379659464542D+00, &
          0.9085620340535459D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.6276714770186411D+00, &
          0.9848870223309694D+00, &
          0.6276714770186411D+00, &
          0.3723285229813589D+00, &
          0.9848870223309694D+00, &
          0.3723285229813589D+00, &
          0.3723285229813589D+00, &
          0.9848870223309694D+00, &
          0.3723285229813589D+00, &
          0.6276714770186411D+00, &
          0.9848870223309694D+00, &
          0.6276714770186411D+00, &
          0.6276714770186411D+00, &
          0.0151129776690307D+00, &
          0.6276714770186411D+00, &
          0.3723285229813589D+00, &
          0.0151129776690307D+00, &
          0.3723285229813589D+00, &
          0.3723285229813589D+00, &
          0.0151129776690307D+00, &
          0.3723285229813589D+00, &
          0.6276714770186411D+00, &
          0.0151129776690307D+00, &
          0.6276714770186411D+00, &
          0.9676181924490221D+00, &
          0.8733616537091691D+00, &
          0.9676181924490221D+00, &
          0.0323818075509779D+00, &
          0.8733616537091691D+00, &
          0.0323818075509779D+00, &
          0.0323818075509779D+00, &
          0.8733616537091691D+00, &
          0.0323818075509779D+00, &
          0.9676181924490221D+00, &
          0.8733616537091691D+00, &
          0.9676181924490221D+00, &
          0.9676181924490221D+00, &
          0.1266383462908310D+00, &
          0.9676181924490221D+00, &
          0.0323818075509779D+00, &
          0.1266383462908310D+00, &
          0.0323818075509779D+00, &
          0.0323818075509779D+00, &
          0.1266383462908310D+00, &
          0.0323818075509779D+00, &
          0.9676181924490221D+00, &
          0.1266383462908310D+00, &
          0.9676181924490221D+00, &
          0.6715295730694326D+00, &
          0.8093319789042792D+00, &
          0.6715295730694326D+00, &
          0.3284704269305673D+00, &
          0.8093319789042792D+00, &
          0.3284704269305673D+00, &
          0.3284704269305673D+00, &
          0.8093319789042792D+00, &
          0.3284704269305673D+00, &
          0.6715295730694326D+00, &
          0.8093319789042792D+00, &
          0.6715295730694326D+00, &
          0.6715295730694326D+00, &
          0.1906680210957207D+00, &
          0.6715295730694326D+00, &
          0.3284704269305673D+00, &
          0.1906680210957207D+00, &
          0.3284704269305673D+00, &
          0.3284704269305673D+00, &
          0.1906680210957207D+00, &
          0.3284704269305673D+00, &
          0.6715295730694326D+00, &
          0.1906680210957207D+00, &
          0.6715295730694326D+00, &
          0.6277920800575361D+00, &
          0.4229648752978111D+00, &
          0.6277920800575361D+00, &
          0.3722079199424639D+00, &
          0.4229648752978111D+00, &
          0.3722079199424639D+00, &
          0.3722079199424639D+00, &
          0.4229648752978111D+00, &
          0.3722079199424639D+00, &
          0.6277920800575361D+00, &
          0.4229648752978111D+00, &
          0.6277920800575361D+00, &
          0.6277920800575361D+00, &
          0.5770351247021889D+00, &
          0.6277920800575361D+00, &
          0.3722079199424639D+00, &
          0.5770351247021889D+00, &
          0.3722079199424639D+00, &
          0.3722079199424639D+00, &
          0.5770351247021889D+00, &
          0.3722079199424639D+00, &
          0.6277920800575361D+00, &
          0.5770351247021889D+00, &
          0.6277920800575361D+00, &
          0.7708105708234440D+00, &
          0.9478284758927406D+00, &
          0.7708105708234440D+00, &
          0.2291894291765560D+00, &
          0.9478284758927406D+00, &
          0.2291894291765560D+00, &
          0.2291894291765560D+00, &
          0.9478284758927406D+00, &
          0.2291894291765560D+00, &
          0.7708105708234440D+00, &
          0.9478284758927406D+00, &
          0.7708105708234440D+00, &
          0.7708105708234440D+00, &
          0.0521715241072594D+00, &
          0.7708105708234440D+00, &
          0.2291894291765560D+00, &
          0.0521715241072594D+00, &
          0.2291894291765560D+00, &
          0.2291894291765560D+00, &
          0.0521715241072594D+00, &
          0.2291894291765560D+00, &
          0.7708105708234440D+00, &
          0.0521715241072594D+00, &
          0.7708105708234440D+00, &
          0.9437280545720244D+00, &
          0.6207324889467394D+00, &
          0.9437280545720244D+00, &
          0.0562719454279757D+00, &
          0.6207324889467394D+00, &
          0.0562719454279757D+00, &
          0.0562719454279757D+00, &
          0.6207324889467394D+00, &
          0.0562719454279757D+00, &
          0.9437280545720244D+00, &
          0.6207324889467394D+00, &
          0.9437280545720244D+00, &
          0.9437280545720244D+00, &
          0.3792675110532606D+00, &
          0.9437280545720244D+00, &
          0.0562719454279757D+00, &
          0.3792675110532606D+00, &
          0.0562719454279757D+00, &
          0.0562719454279757D+00, &
          0.3792675110532606D+00, &
          0.0562719454279757D+00, &
          0.9437280545720244D+00, &
          0.3792675110532606D+00, &
          0.9437280545720244D+00, &
          0.9007366871015028D+00, &
          0.9007366871015028D+00, &
          0.9943923590724830D+00, &
          0.7341796540172194D+00, &
          0.9943923590724830D+00, &
          0.7341796540172194D+00, &
          0.0992633128984972D+00, &
          0.0992633128984972D+00, &
          0.9943923590724830D+00, &
          0.7341796540172194D+00, &
          0.9943923590724830D+00, &
          0.7341796540172194D+00, &
          0.9007366871015028D+00, &
          0.9007366871015028D+00, &
          0.0056076409275170D+00, &
          0.7341796540172194D+00, &
          0.0056076409275170D+00, &
          0.7341796540172194D+00, &
          0.9007366871015028D+00, &
          0.9007366871015028D+00, &
          0.9943923590724830D+00, &
          0.2658203459827806D+00, &
          0.9943923590724830D+00, &
          0.2658203459827806D+00, &
          0.0992633128984972D+00, &
          0.0992633128984972D+00, &
          0.0056076409275170D+00, &
          0.7341796540172194D+00, &
          0.0056076409275170D+00, &
          0.7341796540172194D+00, &
          0.0992633128984972D+00, &
          0.0992633128984972D+00, &
          0.9943923590724830D+00, &
          0.2658203459827806D+00, &
          0.9943923590724830D+00, &
          0.2658203459827806D+00, &
          0.9007366871015028D+00, &
          0.9007366871015028D+00, &
          0.0056076409275170D+00, &
          0.2658203459827806D+00, &
          0.0056076409275170D+00, &
          0.2658203459827806D+00, &
          0.0992633128984972D+00, &
          0.0992633128984972D+00, &
          0.0056076409275170D+00, &
          0.2658203459827806D+00, &
          0.0056076409275170D+00, &
          0.2658203459827806D+00, &
          0.8049570330025153D+00, &
          0.8049570330025153D+00, &
          0.8934115499715557D+00, &
          0.6602554101905883D+00, &
          0.8934115499715557D+00, &
          0.6602554101905883D+00, &
          0.1950429669974847D+00, &
          0.1950429669974847D+00, &
          0.8934115499715557D+00, &
          0.6602554101905883D+00, &
          0.8934115499715557D+00, &
          0.6602554101905883D+00, &
          0.8049570330025153D+00, &
          0.8049570330025153D+00, &
          0.1065884500284443D+00, &
          0.6602554101905883D+00, &
          0.1065884500284443D+00, &
          0.6602554101905883D+00, &
          0.8049570330025153D+00, &
          0.8049570330025153D+00, &
          0.8934115499715557D+00, &
          0.3397445898094117D+00, &
          0.8934115499715557D+00, &
          0.3397445898094117D+00, &
          0.1950429669974847D+00, &
          0.1950429669974847D+00, &
          0.1065884500284443D+00, &
          0.6602554101905883D+00, &
          0.1065884500284443D+00, &
          0.6602554101905883D+00, &
          0.1950429669974847D+00, &
          0.1950429669974847D+00, &
          0.8934115499715557D+00, &
          0.3397445898094117D+00, &
          0.8934115499715557D+00, &
          0.3397445898094117D+00, &
          0.8049570330025153D+00, &
          0.8049570330025153D+00, &
          0.1065884500284443D+00, &
          0.3397445898094117D+00, &
          0.1065884500284443D+00, &
          0.3397445898094117D+00, &
          0.1950429669974847D+00, &
          0.1950429669974847D+00, &
          0.1065884500284443D+00, &
          0.3397445898094117D+00, &
          0.1065884500284443D+00, &
          0.3397445898094117D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.5000000000000000D+00, &
          0.7426575151151338D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.2573424848848662D+00, &
          0.5000000000000000D+00, &
          0.8703394079370492D+00, &
          0.1296605920629508D+00, &
          0.8703394079370492D+00, &
          0.1296605920629508D+00, &
          0.8703394079370492D+00, &
          0.1296605920629508D+00, &
          0.1296605920629508D+00, &
          0.8703394079370492D+00, &
          0.9999997917911259D+00, &
          0.0000002082088741D+00, &
          0.9999997917911259D+00, &
          0.0000002082088741D+00, &
          0.9999997917911259D+00, &
          0.0000002082088741D+00, &
          0.0000002082088741D+00, &
          0.9999997917911259D+00, &
          0.8773134026065881D+00, &
          0.8773134026065881D+00, &
          0.5000000000000000D+00, &
          0.8773134026065881D+00, &
          0.1226865973934119D+00, &
          0.5000000000000000D+00, &
          0.1226865973934119D+00, &
          0.8773134026065881D+00, &
          0.5000000000000000D+00, &
          0.1226865973934119D+00, &
          0.1226865973934119D+00, &
          0.5000000000000000D+00, &
          0.7748409070843456D+00, &
          0.7748409070843456D+00, &
          0.5000000000000000D+00, &
          0.7748409070843456D+00, &
          0.2251590929156544D+00, &
          0.5000000000000000D+00, &
          0.2251590929156544D+00, &
          0.7748409070843456D+00, &
          0.5000000000000000D+00, &
          0.2251590929156544D+00, &
          0.2251590929156544D+00, &
          0.5000000000000000D+00, &
          0.9998342518711405D+00, &
          0.9998342518711405D+00, &
          0.5000000000000000D+00, &
          0.9998342518711405D+00, &
          0.0001657481288594D+00, &
          0.5000000000000000D+00, &
          0.0001657481288594D+00, &
          0.9998342518711405D+00, &
          0.5000000000000000D+00, &
          0.0001657481288594D+00, &
          0.0001657481288594D+00, &
          0.5000000000000000D+00, &
          0.9765209703524975D+00, &
          0.0234790296475025D+00, &
          0.0234790296475025D+00, &
          0.9765209703524975D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.8078645495842367D+00, &
          0.1921354504157633D+00, &
          0.8078645495842367D+00, &
          0.1921354504157633D+00, &
          0.8078645495842367D+00, &
          0.1921354504157633D+00, &
          0.8078645495842367D+00, &
          0.1921354504157633D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9765209703524975D+00, &
          0.0234790296475025D+00, &
          0.0234790296475025D+00, &
          0.9765209703524975D+00, &
          0.9085620340535459D+00, &
          0.0914379659464542D+00, &
          0.0914379659464542D+00, &
          0.9085620340535459D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.6378406365996775D+00, &
          0.3621593634003225D+00, &
          0.6378406365996775D+00, &
          0.3621593634003225D+00, &
          0.6378406365996775D+00, &
          0.3621593634003225D+00, &
          0.6378406365996775D+00, &
          0.3621593634003225D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9085620340535459D+00, &
          0.0914379659464542D+00, &
          0.0914379659464542D+00, &
          0.9085620340535459D+00, &
          0.9848870223309694D+00, &
          0.6276714770186411D+00, &
          0.6276714770186411D+00, &
          0.9848870223309694D+00, &
          0.3723285229813589D+00, &
          0.6276714770186411D+00, &
          0.9848870223309694D+00, &
          0.3723285229813589D+00, &
          0.3723285229813589D+00, &
          0.9848870223309694D+00, &
          0.6276714770186411D+00, &
          0.3723285229813589D+00, &
          0.0151129776690307D+00, &
          0.6276714770186411D+00, &
          0.6276714770186411D+00, &
          0.0151129776690307D+00, &
          0.3723285229813589D+00, &
          0.6276714770186411D+00, &
          0.0151129776690307D+00, &
          0.3723285229813589D+00, &
          0.3723285229813589D+00, &
          0.0151129776690307D+00, &
          0.6276714770186411D+00, &
          0.3723285229813589D+00, &
          0.8733616537091691D+00, &
          0.9676181924490221D+00, &
          0.9676181924490221D+00, &
          0.8733616537091691D+00, &
          0.0323818075509779D+00, &
          0.9676181924490221D+00, &
          0.8733616537091691D+00, &
          0.0323818075509779D+00, &
          0.0323818075509779D+00, &
          0.8733616537091691D+00, &
          0.9676181924490221D+00, &
          0.0323818075509779D+00, &
          0.1266383462908310D+00, &
          0.9676181924490221D+00, &
          0.9676181924490221D+00, &
          0.1266383462908310D+00, &
          0.0323818075509779D+00, &
          0.9676181924490221D+00, &
          0.1266383462908310D+00, &
          0.0323818075509779D+00, &
          0.0323818075509779D+00, &
          0.1266383462908310D+00, &
          0.9676181924490221D+00, &
          0.0323818075509779D+00, &
          0.8093319789042792D+00, &
          0.6715295730694326D+00, &
          0.6715295730694326D+00, &
          0.8093319789042792D+00, &
          0.3284704269305673D+00, &
          0.6715295730694326D+00, &
          0.8093319789042792D+00, &
          0.3284704269305673D+00, &
          0.3284704269305673D+00, &
          0.8093319789042792D+00, &
          0.6715295730694326D+00, &
          0.3284704269305673D+00, &
          0.1906680210957207D+00, &
          0.6715295730694326D+00, &
          0.6715295730694326D+00, &
          0.1906680210957207D+00, &
          0.3284704269305673D+00, &
          0.6715295730694326D+00, &
          0.1906680210957207D+00, &
          0.3284704269305673D+00, &
          0.3284704269305673D+00, &
          0.1906680210957207D+00, &
          0.6715295730694326D+00, &
          0.3284704269305673D+00, &
          0.4229648752978111D+00, &
          0.6277920800575361D+00, &
          0.6277920800575361D+00, &
          0.4229648752978111D+00, &
          0.3722079199424639D+00, &
          0.6277920800575361D+00, &
          0.4229648752978111D+00, &
          0.3722079199424639D+00, &
          0.3722079199424639D+00, &
          0.4229648752978111D+00, &
          0.6277920800575361D+00, &
          0.3722079199424639D+00, &
          0.5770351247021889D+00, &
          0.6277920800575361D+00, &
          0.6277920800575361D+00, &
          0.5770351247021889D+00, &
          0.3722079199424639D+00, &
          0.6277920800575361D+00, &
          0.5770351247021889D+00, &
          0.3722079199424639D+00, &
          0.3722079199424639D+00, &
          0.5770351247021889D+00, &
          0.6277920800575361D+00, &
          0.3722079199424639D+00, &
          0.9478284758927406D+00, &
          0.7708105708234440D+00, &
          0.7708105708234440D+00, &
          0.9478284758927406D+00, &
          0.2291894291765560D+00, &
          0.7708105708234440D+00, &
          0.9478284758927406D+00, &
          0.2291894291765560D+00, &
          0.2291894291765560D+00, &
          0.9478284758927406D+00, &
          0.7708105708234440D+00, &
          0.2291894291765560D+00, &
          0.0521715241072594D+00, &
          0.7708105708234440D+00, &
          0.7708105708234440D+00, &
          0.0521715241072594D+00, &
          0.2291894291765560D+00, &
          0.7708105708234440D+00, &
          0.0521715241072594D+00, &
          0.2291894291765560D+00, &
          0.2291894291765560D+00, &
          0.0521715241072594D+00, &
          0.7708105708234440D+00, &
          0.2291894291765560D+00, &
          0.6207324889467394D+00, &
          0.9437280545720244D+00, &
          0.9437280545720244D+00, &
          0.6207324889467394D+00, &
          0.0562719454279757D+00, &
          0.9437280545720244D+00, &
          0.6207324889467394D+00, &
          0.0562719454279757D+00, &
          0.0562719454279757D+00, &
          0.6207324889467394D+00, &
          0.9437280545720244D+00, &
          0.0562719454279757D+00, &
          0.3792675110532606D+00, &
          0.9437280545720244D+00, &
          0.9437280545720244D+00, &
          0.3792675110532606D+00, &
          0.0562719454279757D+00, &
          0.9437280545720244D+00, &
          0.3792675110532606D+00, &
          0.0562719454279757D+00, &
          0.0562719454279757D+00, &
          0.3792675110532606D+00, &
          0.9437280545720244D+00, &
          0.0562719454279757D+00, &
          0.9943923590724830D+00, &
          0.7341796540172194D+00, &
          0.9007366871015028D+00, &
          0.9007366871015028D+00, &
          0.7341796540172194D+00, &
          0.9943923590724830D+00, &
          0.9943923590724830D+00, &
          0.7341796540172194D+00, &
          0.0992633128984972D+00, &
          0.0992633128984972D+00, &
          0.7341796540172194D+00, &
          0.9943923590724830D+00, &
          0.0056076409275170D+00, &
          0.7341796540172194D+00, &
          0.9007366871015028D+00, &
          0.9007366871015028D+00, &
          0.7341796540172194D+00, &
          0.0056076409275170D+00, &
          0.9943923590724830D+00, &
          0.2658203459827806D+00, &
          0.9007366871015028D+00, &
          0.9007366871015028D+00, &
          0.2658203459827806D+00, &
          0.9943923590724830D+00, &
          0.0056076409275170D+00, &
          0.7341796540172194D+00, &
          0.0992633128984972D+00, &
          0.0992633128984972D+00, &
          0.7341796540172194D+00, &
          0.0056076409275170D+00, &
          0.9943923590724830D+00, &
          0.2658203459827806D+00, &
          0.0992633128984972D+00, &
          0.0992633128984972D+00, &
          0.2658203459827806D+00, &
          0.9943923590724830D+00, &
          0.0056076409275170D+00, &
          0.2658203459827806D+00, &
          0.9007366871015028D+00, &
          0.9007366871015028D+00, &
          0.2658203459827806D+00, &
          0.0056076409275170D+00, &
          0.0056076409275170D+00, &
          0.2658203459827806D+00, &
          0.0992633128984972D+00, &
          0.0992633128984972D+00, &
          0.2658203459827806D+00, &
          0.0056076409275170D+00, &
          0.8934115499715557D+00, &
          0.6602554101905883D+00, &
          0.8049570330025153D+00, &
          0.8049570330025153D+00, &
          0.6602554101905883D+00, &
          0.8934115499715557D+00, &
          0.8934115499715557D+00, &
          0.6602554101905883D+00, &
          0.1950429669974847D+00, &
          0.1950429669974847D+00, &
          0.6602554101905883D+00, &
          0.8934115499715557D+00, &
          0.1065884500284443D+00, &
          0.6602554101905883D+00, &
          0.8049570330025153D+00, &
          0.8049570330025153D+00, &
          0.6602554101905883D+00, &
          0.1065884500284443D+00, &
          0.8934115499715557D+00, &
          0.3397445898094117D+00, &
          0.8049570330025153D+00, &
          0.8049570330025153D+00, &
          0.3397445898094117D+00, &
          0.8934115499715557D+00, &
          0.1065884500284443D+00, &
          0.6602554101905883D+00, &
          0.1950429669974847D+00, &
          0.1950429669974847D+00, &
          0.6602554101905883D+00, &
          0.1065884500284443D+00, &
          0.8934115499715557D+00, &
          0.3397445898094117D+00, &
          0.1950429669974847D+00, &
          0.1950429669974847D+00, &
          0.3397445898094117D+00, &
          0.8934115499715557D+00, &
          0.1065884500284443D+00, &
          0.3397445898094117D+00, &
          0.8049570330025153D+00, &
          0.8049570330025153D+00, &
          0.3397445898094117D+00, &
          0.1065884500284443D+00, &
          0.1065884500284443D+00, &
          0.3397445898094117D+00, &
          0.1950429669974847D+00, &
          0.1950429669974847D+00, &
          0.3397445898094117D+00, &
          0.1065884500284443D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.7426575151151338D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.2573424848848662D+00, &
          0.8703394079370492D+00, &
          0.1296605920629508D+00, &
          0.8703394079370492D+00, &
          0.8703394079370492D+00, &
          0.1296605920629508D+00, &
          0.8703394079370492D+00, &
          0.1296605920629508D+00, &
          0.1296605920629508D+00, &
          0.9999997917911259D+00, &
          0.0000002082088741D+00, &
          0.9999997917911259D+00, &
          0.9999997917911259D+00, &
          0.0000002082088741D+00, &
          0.9999997917911259D+00, &
          0.0000002082088741D+00, &
          0.0000002082088741D+00, &
          0.5000000000000000D+00, &
          0.8773134026065881D+00, &
          0.8773134026065881D+00, &
          0.5000000000000000D+00, &
          0.8773134026065881D+00, &
          0.8773134026065881D+00, &
          0.5000000000000000D+00, &
          0.1226865973934119D+00, &
          0.1226865973934119D+00, &
          0.5000000000000000D+00, &
          0.1226865973934119D+00, &
          0.1226865973934119D+00, &
          0.5000000000000000D+00, &
          0.7748409070843456D+00, &
          0.7748409070843456D+00, &
          0.5000000000000000D+00, &
          0.7748409070843456D+00, &
          0.7748409070843456D+00, &
          0.5000000000000000D+00, &
          0.2251590929156544D+00, &
          0.2251590929156544D+00, &
          0.5000000000000000D+00, &
          0.2251590929156544D+00, &
          0.2251590929156544D+00, &
          0.5000000000000000D+00, &
          0.9998342518711405D+00, &
          0.9998342518711405D+00, &
          0.5000000000000000D+00, &
          0.9998342518711405D+00, &
          0.9998342518711405D+00, &
          0.5000000000000000D+00, &
          0.0001657481288594D+00, &
          0.0001657481288594D+00, &
          0.5000000000000000D+00, &
          0.0001657481288594D+00, &
          0.0001657481288594D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9765209703524975D+00, &
          0.0234790296475025D+00, &
          0.0234790296475025D+00, &
          0.9765209703524975D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9765209703524975D+00, &
          0.0234790296475025D+00, &
          0.0234790296475025D+00, &
          0.9765209703524975D+00, &
          0.8078645495842367D+00, &
          0.1921354504157633D+00, &
          0.8078645495842367D+00, &
          0.1921354504157633D+00, &
          0.8078645495842367D+00, &
          0.1921354504157633D+00, &
          0.8078645495842367D+00, &
          0.1921354504157633D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9085620340535459D+00, &
          0.0914379659464542D+00, &
          0.0914379659464542D+00, &
          0.9085620340535459D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9085620340535459D+00, &
          0.0914379659464542D+00, &
          0.0914379659464542D+00, &
          0.9085620340535459D+00, &
          0.6378406365996775D+00, &
          0.3621593634003225D+00, &
          0.6378406365996775D+00, &
          0.3621593634003225D+00, &
          0.6378406365996775D+00, &
          0.3621593634003225D+00, &
          0.6378406365996775D+00, &
          0.3621593634003225D+00, &
          0.6276714770186411D+00, &
          0.6276714770186411D+00, &
          0.9848870223309694D+00, &
          0.6276714770186411D+00, &
          0.6276714770186411D+00, &
          0.9848870223309694D+00, &
          0.3723285229813589D+00, &
          0.3723285229813589D+00, &
          0.9848870223309694D+00, &
          0.3723285229813589D+00, &
          0.3723285229813589D+00, &
          0.9848870223309694D+00, &
          0.6276714770186411D+00, &
          0.6276714770186411D+00, &
          0.0151129776690307D+00, &
          0.6276714770186411D+00, &
          0.6276714770186411D+00, &
          0.0151129776690307D+00, &
          0.3723285229813589D+00, &
          0.3723285229813589D+00, &
          0.0151129776690307D+00, &
          0.3723285229813589D+00, &
          0.3723285229813589D+00, &
          0.0151129776690307D+00, &
          0.9676181924490221D+00, &
          0.9676181924490221D+00, &
          0.8733616537091691D+00, &
          0.9676181924490221D+00, &
          0.9676181924490221D+00, &
          0.8733616537091691D+00, &
          0.0323818075509779D+00, &
          0.0323818075509779D+00, &
          0.8733616537091691D+00, &
          0.0323818075509779D+00, &
          0.0323818075509779D+00, &
          0.8733616537091691D+00, &
          0.9676181924490221D+00, &
          0.9676181924490221D+00, &
          0.1266383462908310D+00, &
          0.9676181924490221D+00, &
          0.9676181924490221D+00, &
          0.1266383462908310D+00, &
          0.0323818075509779D+00, &
          0.0323818075509779D+00, &
          0.1266383462908310D+00, &
          0.0323818075509779D+00, &
          0.0323818075509779D+00, &
          0.1266383462908310D+00, &
          0.6715295730694326D+00, &
          0.6715295730694326D+00, &
          0.8093319789042792D+00, &
          0.6715295730694326D+00, &
          0.6715295730694326D+00, &
          0.8093319789042792D+00, &
          0.3284704269305673D+00, &
          0.3284704269305673D+00, &
          0.8093319789042792D+00, &
          0.3284704269305673D+00, &
          0.3284704269305673D+00, &
          0.8093319789042792D+00, &
          0.6715295730694326D+00, &
          0.6715295730694326D+00, &
          0.1906680210957207D+00, &
          0.6715295730694326D+00, &
          0.6715295730694326D+00, &
          0.1906680210957207D+00, &
          0.3284704269305673D+00, &
          0.3284704269305673D+00, &
          0.1906680210957207D+00, &
          0.3284704269305673D+00, &
          0.3284704269305673D+00, &
          0.1906680210957207D+00, &
          0.6277920800575361D+00, &
          0.6277920800575361D+00, &
          0.4229648752978111D+00, &
          0.6277920800575361D+00, &
          0.6277920800575361D+00, &
          0.4229648752978111D+00, &
          0.3722079199424639D+00, &
          0.3722079199424639D+00, &
          0.4229648752978111D+00, &
          0.3722079199424639D+00, &
          0.3722079199424639D+00, &
          0.4229648752978111D+00, &
          0.6277920800575361D+00, &
          0.6277920800575361D+00, &
          0.5770351247021889D+00, &
          0.6277920800575361D+00, &
          0.6277920800575361D+00, &
          0.5770351247021889D+00, &
          0.3722079199424639D+00, &
          0.3722079199424639D+00, &
          0.5770351247021889D+00, &
          0.3722079199424639D+00, &
          0.3722079199424639D+00, &
          0.5770351247021889D+00, &
          0.7708105708234440D+00, &
          0.7708105708234440D+00, &
          0.9478284758927406D+00, &
          0.7708105708234440D+00, &
          0.7708105708234440D+00, &
          0.9478284758927406D+00, &
          0.2291894291765560D+00, &
          0.2291894291765560D+00, &
          0.9478284758927406D+00, &
          0.2291894291765560D+00, &
          0.2291894291765560D+00, &
          0.9478284758927406D+00, &
          0.7708105708234440D+00, &
          0.7708105708234440D+00, &
          0.0521715241072594D+00, &
          0.7708105708234440D+00, &
          0.7708105708234440D+00, &
          0.0521715241072594D+00, &
          0.2291894291765560D+00, &
          0.2291894291765560D+00, &
          0.0521715241072594D+00, &
          0.2291894291765560D+00, &
          0.2291894291765560D+00, &
          0.0521715241072594D+00, &
          0.9437280545720244D+00, &
          0.9437280545720244D+00, &
          0.6207324889467394D+00, &
          0.9437280545720244D+00, &
          0.9437280545720244D+00, &
          0.6207324889467394D+00, &
          0.0562719454279757D+00, &
          0.0562719454279757D+00, &
          0.6207324889467394D+00, &
          0.0562719454279757D+00, &
          0.0562719454279757D+00, &
          0.6207324889467394D+00, &
          0.9437280545720244D+00, &
          0.9437280545720244D+00, &
          0.3792675110532606D+00, &
          0.9437280545720244D+00, &
          0.9437280545720244D+00, &
          0.3792675110532606D+00, &
          0.0562719454279757D+00, &
          0.0562719454279757D+00, &
          0.3792675110532606D+00, &
          0.0562719454279757D+00, &
          0.0562719454279757D+00, &
          0.3792675110532606D+00, &
          0.7341796540172194D+00, &
          0.9943923590724830D+00, &
          0.7341796540172194D+00, &
          0.9943923590724830D+00, &
          0.9007366871015028D+00, &
          0.9007366871015028D+00, &
          0.7341796540172194D+00, &
          0.9943923590724830D+00, &
          0.7341796540172194D+00, &
          0.9943923590724830D+00, &
          0.0992633128984972D+00, &
          0.0992633128984972D+00, &
          0.7341796540172194D+00, &
          0.0056076409275170D+00, &
          0.7341796540172194D+00, &
          0.0056076409275170D+00, &
          0.9007366871015028D+00, &
          0.9007366871015028D+00, &
          0.2658203459827806D+00, &
          0.9943923590724830D+00, &
          0.2658203459827806D+00, &
          0.9943923590724830D+00, &
          0.9007366871015028D+00, &
          0.9007366871015028D+00, &
          0.7341796540172194D+00, &
          0.0056076409275170D+00, &
          0.7341796540172194D+00, &
          0.0056076409275170D+00, &
          0.0992633128984972D+00, &
          0.0992633128984972D+00, &
          0.2658203459827806D+00, &
          0.9943923590724830D+00, &
          0.2658203459827806D+00, &
          0.9943923590724830D+00, &
          0.0992633128984972D+00, &
          0.0992633128984972D+00, &
          0.2658203459827806D+00, &
          0.0056076409275170D+00, &
          0.2658203459827806D+00, &
          0.0056076409275170D+00, &
          0.9007366871015028D+00, &
          0.9007366871015028D+00, &
          0.2658203459827806D+00, &
          0.0056076409275170D+00, &
          0.2658203459827806D+00, &
          0.0056076409275170D+00, &
          0.0992633128984972D+00, &
          0.0992633128984972D+00, &
          0.6602554101905883D+00, &
          0.8934115499715557D+00, &
          0.6602554101905883D+00, &
          0.8934115499715557D+00, &
          0.8049570330025153D+00, &
          0.8049570330025153D+00, &
          0.6602554101905883D+00, &
          0.8934115499715557D+00, &
          0.6602554101905883D+00, &
          0.8934115499715557D+00, &
          0.1950429669974847D+00, &
          0.1950429669974847D+00, &
          0.6602554101905883D+00, &
          0.1065884500284443D+00, &
          0.6602554101905883D+00, &
          0.1065884500284443D+00, &
          0.8049570330025153D+00, &
          0.8049570330025153D+00, &
          0.3397445898094117D+00, &
          0.8934115499715557D+00, &
          0.3397445898094117D+00, &
          0.8934115499715557D+00, &
          0.8049570330025153D+00, &
          0.8049570330025153D+00, &
          0.6602554101905883D+00, &
          0.1065884500284443D+00, &
          0.6602554101905883D+00, &
          0.1065884500284443D+00, &
          0.1950429669974847D+00, &
          0.1950429669974847D+00, &
          0.3397445898094117D+00, &
          0.8934115499715557D+00, &
          0.3397445898094117D+00, &
          0.8934115499715557D+00, &
          0.1950429669974847D+00, &
          0.1950429669974847D+00, &
          0.3397445898094117D+00, &
          0.1065884500284443D+00, &
          0.3397445898094117D+00, &
          0.1065884500284443D+00, &
          0.8049570330025153D+00, &
          0.8049570330025153D+00, &
          0.3397445898094117D+00, &
          0.1065884500284443D+00, &
          0.3397445898094117D+00, &
          0.1065884500284443D+00, &
          0.1950429669974847D+00, &
          0.1950429669974847D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          0.0096145921289297D+00, &
          0.0096145921289297D+00, &
          0.0096145921289297D+00, &
          0.0096145921289297D+00, &
          0.0096145921289297D+00, &
          0.0096145921289297D+00, &
          0.0043938276021590D+00, &
          0.0043938276021590D+00, &
          0.0043938276021590D+00, &
          0.0043938276021590D+00, &
          0.0043938276021590D+00, &
          0.0043938276021590D+00, &
          0.0043938276021590D+00, &
          0.0043938276021590D+00, &
          0.0000496472519452D+00, &
          0.0000496472519452D+00, &
          0.0000496472519452D+00, &
          0.0000496472519452D+00, &
          0.0000496472519452D+00, &
          0.0000496472519452D+00, &
          0.0000496472519452D+00, &
          0.0000496472519452D+00, &
          0.0025062844648325D+00, &
          0.0025062844648325D+00, &
          0.0025062844648325D+00, &
          0.0025062844648325D+00, &
          0.0025062844648325D+00, &
          0.0025062844648325D+00, &
          0.0025062844648325D+00, &
          0.0025062844648325D+00, &
          0.0025062844648325D+00, &
          0.0025062844648325D+00, &
          0.0025062844648325D+00, &
          0.0025062844648325D+00, &
          0.0061332324820703D+00, &
          0.0061332324820703D+00, &
          0.0061332324820703D+00, &
          0.0061332324820703D+00, &
          0.0061332324820703D+00, &
          0.0061332324820703D+00, &
          0.0061332324820703D+00, &
          0.0061332324820703D+00, &
          0.0061332324820703D+00, &
          0.0061332324820703D+00, &
          0.0061332324820703D+00, &
          0.0061332324820703D+00, &
          0.0003052513689971D+00, &
          0.0003052513689971D+00, &
          0.0003052513689971D+00, &
          0.0003052513689971D+00, &
          0.0003052513689971D+00, &
          0.0003052513689971D+00, &
          0.0003052513689971D+00, &
          0.0003052513689971D+00, &
          0.0003052513689971D+00, &
          0.0003052513689971D+00, &
          0.0003052513689971D+00, &
          0.0003052513689971D+00, &
          0.0028254377330331D+00, &
          0.0028254377330331D+00, &
          0.0028254377330331D+00, &
          0.0028254377330331D+00, &
          0.0028254377330331D+00, &
          0.0028254377330331D+00, &
          0.0028254377330331D+00, &
          0.0028254377330331D+00, &
          0.0028254377330331D+00, &
          0.0028254377330331D+00, &
          0.0028254377330331D+00, &
          0.0028254377330331D+00, &
          0.0028254377330331D+00, &
          0.0028254377330331D+00, &
          0.0028254377330331D+00, &
          0.0028254377330331D+00, &
          0.0028254377330331D+00, &
          0.0028254377330331D+00, &
          0.0028254377330331D+00, &
          0.0028254377330331D+00, &
          0.0028254377330331D+00, &
          0.0028254377330331D+00, &
          0.0028254377330331D+00, &
          0.0028254377330331D+00, &
          0.0048177045987999D+00, &
          0.0048177045987999D+00, &
          0.0048177045987999D+00, &
          0.0048177045987999D+00, &
          0.0048177045987999D+00, &
          0.0048177045987999D+00, &
          0.0048177045987999D+00, &
          0.0048177045987999D+00, &
          0.0048177045987999D+00, &
          0.0048177045987999D+00, &
          0.0048177045987999D+00, &
          0.0048177045987999D+00, &
          0.0048177045987999D+00, &
          0.0048177045987999D+00, &
          0.0048177045987999D+00, &
          0.0048177045987999D+00, &
          0.0048177045987999D+00, &
          0.0048177045987999D+00, &
          0.0048177045987999D+00, &
          0.0048177045987999D+00, &
          0.0048177045987999D+00, &
          0.0048177045987999D+00, &
          0.0048177045987999D+00, &
          0.0048177045987999D+00, &
          0.0022831423389536D+00, &
          0.0022831423389536D+00, &
          0.0022831423389536D+00, &
          0.0022831423389536D+00, &
          0.0022831423389536D+00, &
          0.0022831423389536D+00, &
          0.0022831423389536D+00, &
          0.0022831423389536D+00, &
          0.0022831423389536D+00, &
          0.0022831423389536D+00, &
          0.0022831423389536D+00, &
          0.0022831423389536D+00, &
          0.0022831423389536D+00, &
          0.0022831423389536D+00, &
          0.0022831423389536D+00, &
          0.0022831423389536D+00, &
          0.0022831423389536D+00, &
          0.0022831423389536D+00, &
          0.0022831423389536D+00, &
          0.0022831423389536D+00, &
          0.0022831423389536D+00, &
          0.0022831423389536D+00, &
          0.0022831423389536D+00, &
          0.0022831423389536D+00, &
          0.0011774374872787D+00, &
          0.0011774374872787D+00, &
          0.0011774374872787D+00, &
          0.0011774374872787D+00, &
          0.0011774374872787D+00, &
          0.0011774374872787D+00, &
          0.0011774374872787D+00, &
          0.0011774374872787D+00, &
          0.0011774374872787D+00, &
          0.0011774374872787D+00, &
          0.0011774374872787D+00, &
          0.0011774374872787D+00, &
          0.0011774374872787D+00, &
          0.0011774374872787D+00, &
          0.0011774374872787D+00, &
          0.0011774374872787D+00, &
          0.0011774374872787D+00, &
          0.0011774374872787D+00, &
          0.0011774374872787D+00, &
          0.0011774374872787D+00, &
          0.0011774374872787D+00, &
          0.0011774374872787D+00, &
          0.0011774374872787D+00, &
          0.0011774374872787D+00, &
          0.0061840607506279D+00, &
          0.0061840607506279D+00, &
          0.0061840607506279D+00, &
          0.0061840607506279D+00, &
          0.0061840607506279D+00, &
          0.0061840607506279D+00, &
          0.0061840607506279D+00, &
          0.0061840607506279D+00, &
          0.0061840607506279D+00, &
          0.0061840607506279D+00, &
          0.0061840607506279D+00, &
          0.0061840607506279D+00, &
          0.0061840607506279D+00, &
          0.0061840607506279D+00, &
          0.0061840607506279D+00, &
          0.0061840607506279D+00, &
          0.0061840607506279D+00, &
          0.0061840607506279D+00, &
          0.0061840607506279D+00, &
          0.0061840607506279D+00, &
          0.0061840607506279D+00, &
          0.0061840607506279D+00, &
          0.0061840607506279D+00, &
          0.0061840607506279D+00, &
          0.0026746558375937D+00, &
          0.0026746558375937D+00, &
          0.0026746558375937D+00, &
          0.0026746558375937D+00, &
          0.0026746558375937D+00, &
          0.0026746558375937D+00, &
          0.0026746558375937D+00, &
          0.0026746558375937D+00, &
          0.0026746558375937D+00, &
          0.0026746558375937D+00, &
          0.0026746558375937D+00, &
          0.0026746558375937D+00, &
          0.0026746558375937D+00, &
          0.0026746558375937D+00, &
          0.0026746558375937D+00, &
          0.0026746558375937D+00, &
          0.0026746558375937D+00, &
          0.0026746558375937D+00, &
          0.0026746558375937D+00, &
          0.0026746558375937D+00, &
          0.0026746558375937D+00, &
          0.0026746558375937D+00, &
          0.0026746558375937D+00, &
          0.0026746558375937D+00, &
          0.0034191776775788D+00, &
          0.0034191776775788D+00, &
          0.0034191776775788D+00, &
          0.0034191776775788D+00, &
          0.0034191776775788D+00, &
          0.0034191776775788D+00, &
          0.0034191776775788D+00, &
          0.0034191776775788D+00, &
          0.0034191776775788D+00, &
          0.0034191776775788D+00, &
          0.0034191776775788D+00, &
          0.0034191776775788D+00, &
          0.0034191776775788D+00, &
          0.0034191776775788D+00, &
          0.0034191776775788D+00, &
          0.0034191776775788D+00, &
          0.0034191776775788D+00, &
          0.0034191776775788D+00, &
          0.0034191776775788D+00, &
          0.0034191776775788D+00, &
          0.0034191776775788D+00, &
          0.0034191776775788D+00, &
          0.0034191776775788D+00, &
          0.0034191776775788D+00, &
          0.0023346175566150D+00, &
          0.0023346175566150D+00, &
          0.0023346175566150D+00, &
          0.0023346175566150D+00, &
          0.0023346175566150D+00, &
          0.0023346175566150D+00, &
          0.0023346175566150D+00, &
          0.0023346175566150D+00, &
          0.0023346175566150D+00, &
          0.0023346175566150D+00, &
          0.0023346175566150D+00, &
          0.0023346175566150D+00, &
          0.0023346175566150D+00, &
          0.0023346175566150D+00, &
          0.0023346175566150D+00, &
          0.0023346175566150D+00, &
          0.0023346175566150D+00, &
          0.0023346175566150D+00, &
          0.0023346175566150D+00, &
          0.0023346175566150D+00, &
          0.0023346175566150D+00, &
          0.0023346175566150D+00, &
          0.0023346175566150D+00, &
          0.0023346175566150D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0010527901648871D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00, &
          0.0027438309407639D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule19 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule19() returns the hexahedron quadrature rule of precision 19.
!
!  Discussion:
!
!    We suppose we are given a hexahedron H with vertices
!    (0,0,0), (1,0,0), (1,1,0), (0,1,0), (0,0,1), (1,0,1), (1,1,1), (0,1,1),
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over H is 
!    approximated by Q as follows:
!
!    Q = volume(H) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 May 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order symmetric cubature rules for tetrahedra and pyramids,
!    International Journal for Numerical Methods in Engineering,
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
!    real w(n): the quadrature weights, which add to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 454

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
          0.6178975194487687D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.3821024805512313D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9442475200220319D+00, &
          0.0557524799779681D+00, &
          0.0557524799779681D+00, &
          0.9442475200220319D+00, &
          0.9442475200220319D+00, &
          0.0557524799779681D+00, &
          0.9442475200220319D+00, &
          0.0557524799779681D+00, &
          0.8345909234965025D+00, &
          0.1654090765034974D+00, &
          0.1654090765034974D+00, &
          0.8345909234965025D+00, &
          0.8345909234965025D+00, &
          0.1654090765034974D+00, &
          0.8345909234965025D+00, &
          0.1654090765034974D+00, &
          0.8561550198509273D+00, &
          0.1438449801490727D+00, &
          0.8561550198509273D+00, &
          0.1438449801490727D+00, &
          0.8561550198509273D+00, &
          0.1438449801490727D+00, &
          0.8561550198509273D+00, &
          0.1438449801490727D+00, &
          0.9991488896425311D+00, &
          0.0008511103574689D+00, &
          0.0008511103574689D+00, &
          0.9991488896425311D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9991488896425311D+00, &
          0.0008511103574689D+00, &
          0.0008511103574689D+00, &
          0.9991488896425311D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.6998570581371140D+00, &
          0.3001429418628860D+00, &
          0.6998570581371140D+00, &
          0.3001429418628860D+00, &
          0.6998570581371140D+00, &
          0.3001429418628860D+00, &
          0.6998570581371140D+00, &
          0.3001429418628860D+00, &
          0.5860120647237821D+00, &
          0.4139879352762179D+00, &
          0.4139879352762179D+00, &
          0.5860120647237821D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5860120647237821D+00, &
          0.4139879352762179D+00, &
          0.4139879352762179D+00, &
          0.5860120647237821D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9382902632360840D+00, &
          0.0617097367639160D+00, &
          0.9382902632360840D+00, &
          0.0617097367639160D+00, &
          0.9382902632360840D+00, &
          0.0617097367639160D+00, &
          0.9382902632360840D+00, &
          0.0617097367639160D+00, &
          0.2840454955532042D+00, &
          0.7159545044467958D+00, &
          0.7159545044467958D+00, &
          0.2840454955532042D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.2840454955532042D+00, &
          0.7159545044467958D+00, &
          0.7159545044467958D+00, &
          0.2840454955532042D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.3449134012977228D+00, &
          0.7705521820675761D+00, &
          0.3449134012977228D+00, &
          0.6550865987022771D+00, &
          0.7705521820675761D+00, &
          0.6550865987022771D+00, &
          0.6550865987022771D+00, &
          0.7705521820675761D+00, &
          0.6550865987022771D+00, &
          0.3449134012977228D+00, &
          0.7705521820675761D+00, &
          0.3449134012977228D+00, &
          0.3449134012977228D+00, &
          0.2294478179324239D+00, &
          0.3449134012977228D+00, &
          0.6550865987022771D+00, &
          0.2294478179324239D+00, &
          0.6550865987022771D+00, &
          0.6550865987022771D+00, &
          0.2294478179324239D+00, &
          0.6550865987022771D+00, &
          0.3449134012977228D+00, &
          0.2294478179324239D+00, &
          0.3449134012977228D+00, &
          0.9920926274878630D+00, &
          0.9103324974446185D+00, &
          0.9920926274878630D+00, &
          0.0079073725121370D+00, &
          0.9103324974446185D+00, &
          0.0079073725121370D+00, &
          0.0079073725121370D+00, &
          0.9103324974446185D+00, &
          0.0079073725121370D+00, &
          0.9920926274878630D+00, &
          0.9103324974446185D+00, &
          0.9920926274878630D+00, &
          0.9920926274878630D+00, &
          0.0896675025553814D+00, &
          0.9920926274878630D+00, &
          0.0079073725121370D+00, &
          0.0896675025553814D+00, &
          0.0079073725121370D+00, &
          0.0079073725121370D+00, &
          0.0896675025553814D+00, &
          0.0079073725121370D+00, &
          0.9920926274878630D+00, &
          0.0896675025553814D+00, &
          0.9920926274878630D+00, &
          0.9221936167657927D+00, &
          0.4263116836838957D+00, &
          0.9221936167657927D+00, &
          0.0778063832342073D+00, &
          0.4263116836838957D+00, &
          0.0778063832342073D+00, &
          0.0778063832342073D+00, &
          0.4263116836838957D+00, &
          0.0778063832342073D+00, &
          0.9221936167657927D+00, &
          0.4263116836838957D+00, &
          0.9221936167657927D+00, &
          0.9221936167657927D+00, &
          0.5736883163161043D+00, &
          0.9221936167657927D+00, &
          0.0778063832342073D+00, &
          0.5736883163161043D+00, &
          0.0778063832342073D+00, &
          0.0778063832342073D+00, &
          0.5736883163161043D+00, &
          0.0778063832342073D+00, &
          0.9221936167657927D+00, &
          0.5736883163161043D+00, &
          0.9221936167657927D+00, &
          0.7437996999426106D+00, &
          0.8837457372516819D+00, &
          0.7437996999426106D+00, &
          0.2562003000573894D+00, &
          0.8837457372516819D+00, &
          0.2562003000573894D+00, &
          0.2562003000573894D+00, &
          0.8837457372516819D+00, &
          0.2562003000573894D+00, &
          0.7437996999426106D+00, &
          0.8837457372516819D+00, &
          0.7437996999426106D+00, &
          0.7437996999426106D+00, &
          0.1162542627483181D+00, &
          0.7437996999426106D+00, &
          0.2562003000573894D+00, &
          0.1162542627483181D+00, &
          0.2562003000573894D+00, &
          0.2562003000573894D+00, &
          0.1162542627483181D+00, &
          0.2562003000573894D+00, &
          0.7437996999426106D+00, &
          0.1162542627483181D+00, &
          0.7437996999426106D+00, &
          0.8820659248983722D+00, &
          0.6931955182506035D+00, &
          0.8820659248983722D+00, &
          0.1179340751016278D+00, &
          0.6931955182506035D+00, &
          0.1179340751016278D+00, &
          0.1179340751016278D+00, &
          0.6931955182506035D+00, &
          0.1179340751016278D+00, &
          0.8820659248983722D+00, &
          0.6931955182506035D+00, &
          0.8820659248983722D+00, &
          0.8820659248983722D+00, &
          0.3068044817493964D+00, &
          0.8820659248983722D+00, &
          0.1179340751016278D+00, &
          0.3068044817493964D+00, &
          0.1179340751016278D+00, &
          0.1179340751016278D+00, &
          0.3068044817493964D+00, &
          0.1179340751016278D+00, &
          0.8820659248983722D+00, &
          0.3068044817493964D+00, &
          0.8820659248983722D+00, &
          0.8556268694112152D+00, &
          0.9516006474964566D+00, &
          0.8556268694112152D+00, &
          0.1443731305887847D+00, &
          0.9516006474964566D+00, &
          0.1443731305887847D+00, &
          0.1443731305887847D+00, &
          0.9516006474964566D+00, &
          0.1443731305887847D+00, &
          0.8556268694112152D+00, &
          0.9516006474964566D+00, &
          0.8556268694112152D+00, &
          0.8556268694112152D+00, &
          0.0483993525035434D+00, &
          0.8556268694112152D+00, &
          0.1443731305887847D+00, &
          0.0483993525035434D+00, &
          0.1443731305887847D+00, &
          0.1443731305887847D+00, &
          0.0483993525035434D+00, &
          0.1443731305887847D+00, &
          0.8556268694112152D+00, &
          0.0483993525035434D+00, &
          0.8556268694112152D+00, &
          0.9775744747969544D+00, &
          0.6113080322285057D+00, &
          0.9775744747969544D+00, &
          0.0224255252030456D+00, &
          0.6113080322285057D+00, &
          0.0224255252030456D+00, &
          0.0224255252030456D+00, &
          0.6113080322285057D+00, &
          0.0224255252030456D+00, &
          0.9775744747969544D+00, &
          0.6113080322285057D+00, &
          0.9775744747969544D+00, &
          0.9775744747969544D+00, &
          0.3886919677714943D+00, &
          0.9775744747969544D+00, &
          0.0224255252030456D+00, &
          0.3886919677714943D+00, &
          0.0224255252030456D+00, &
          0.0224255252030456D+00, &
          0.3886919677714943D+00, &
          0.0224255252030456D+00, &
          0.9775744747969544D+00, &
          0.3886919677714943D+00, &
          0.9775744747969544D+00, &
          0.5848459901934986D+00, &
          0.9825814337796190D+00, &
          0.5848459901934986D+00, &
          0.4151540098065015D+00, &
          0.9825814337796190D+00, &
          0.4151540098065015D+00, &
          0.4151540098065015D+00, &
          0.9825814337796190D+00, &
          0.4151540098065015D+00, &
          0.5848459901934986D+00, &
          0.9825814337796190D+00, &
          0.5848459901934986D+00, &
          0.5848459901934986D+00, &
          0.0174185662203810D+00, &
          0.5848459901934986D+00, &
          0.4151540098065015D+00, &
          0.0174185662203810D+00, &
          0.4151540098065015D+00, &
          0.4151540098065015D+00, &
          0.0174185662203810D+00, &
          0.4151540098065015D+00, &
          0.5848459901934986D+00, &
          0.0174185662203810D+00, &
          0.5848459901934986D+00, &
          0.8042667629417029D+00, &
          0.5470636340177077D+00, &
          0.8042667629417029D+00, &
          0.1957332370582971D+00, &
          0.5470636340177077D+00, &
          0.1957332370582971D+00, &
          0.1957332370582971D+00, &
          0.5470636340177077D+00, &
          0.1957332370582971D+00, &
          0.8042667629417029D+00, &
          0.5470636340177077D+00, &
          0.8042667629417029D+00, &
          0.8042667629417029D+00, &
          0.4529363659822923D+00, &
          0.8042667629417029D+00, &
          0.1957332370582971D+00, &
          0.4529363659822923D+00, &
          0.1957332370582971D+00, &
          0.1957332370582971D+00, &
          0.4529363659822923D+00, &
          0.1957332370582971D+00, &
          0.8042667629417029D+00, &
          0.4529363659822923D+00, &
          0.8042667629417029D+00, &
          0.7132758893843609D+00, &
          0.9892390565548020D+00, &
          0.7132758893843609D+00, &
          0.2867241106156391D+00, &
          0.9892390565548020D+00, &
          0.2867241106156391D+00, &
          0.2867241106156391D+00, &
          0.9892390565548020D+00, &
          0.2867241106156391D+00, &
          0.7132758893843609D+00, &
          0.9892390565548020D+00, &
          0.7132758893843609D+00, &
          0.7132758893843609D+00, &
          0.0107609434451980D+00, &
          0.7132758893843609D+00, &
          0.2867241106156391D+00, &
          0.0107609434451980D+00, &
          0.2867241106156391D+00, &
          0.2867241106156391D+00, &
          0.0107609434451980D+00, &
          0.2867241106156391D+00, &
          0.7132758893843609D+00, &
          0.0107609434451980D+00, &
          0.7132758893843609D+00, &
          0.5918460351304471D+00, &
          0.8728877113378553D+00, &
          0.5918460351304471D+00, &
          0.4081539648695529D+00, &
          0.8728877113378553D+00, &
          0.4081539648695529D+00, &
          0.4081539648695529D+00, &
          0.8728877113378553D+00, &
          0.4081539648695529D+00, &
          0.5918460351304471D+00, &
          0.8728877113378553D+00, &
          0.5918460351304471D+00, &
          0.5918460351304471D+00, &
          0.1271122886621447D+00, &
          0.5918460351304471D+00, &
          0.4081539648695529D+00, &
          0.1271122886621447D+00, &
          0.4081539648695529D+00, &
          0.4081539648695529D+00, &
          0.1271122886621447D+00, &
          0.4081539648695529D+00, &
          0.5918460351304471D+00, &
          0.1271122886621447D+00, &
          0.5918460351304471D+00, &
          0.8267665443298317D+00, &
          0.8267665443298317D+00, &
          0.9566761563025079D+00, &
          0.6444314510393020D+00, &
          0.9566761563025079D+00, &
          0.6444314510393020D+00, &
          0.1732334556701683D+00, &
          0.1732334556701683D+00, &
          0.9566761563025079D+00, &
          0.6444314510393020D+00, &
          0.9566761563025079D+00, &
          0.6444314510393020D+00, &
          0.8267665443298317D+00, &
          0.8267665443298317D+00, &
          0.0433238436974921D+00, &
          0.6444314510393020D+00, &
          0.0433238436974921D+00, &
          0.6444314510393020D+00, &
          0.8267665443298317D+00, &
          0.8267665443298317D+00, &
          0.9566761563025079D+00, &
          0.3555685489606980D+00, &
          0.9566761563025079D+00, &
          0.3555685489606980D+00, &
          0.1732334556701683D+00, &
          0.1732334556701683D+00, &
          0.0433238436974921D+00, &
          0.6444314510393020D+00, &
          0.0433238436974921D+00, &
          0.6444314510393020D+00, &
          0.1732334556701683D+00, &
          0.1732334556701683D+00, &
          0.9566761563025079D+00, &
          0.3555685489606980D+00, &
          0.9566761563025079D+00, &
          0.3555685489606980D+00, &
          0.8267665443298317D+00, &
          0.8267665443298317D+00, &
          0.0433238436974921D+00, &
          0.3555685489606980D+00, &
          0.0433238436974921D+00, &
          0.3555685489606980D+00, &
          0.1732334556701683D+00, &
          0.1732334556701683D+00, &
          0.0433238436974921D+00, &
          0.3555685489606980D+00, &
          0.0433238436974921D+00, &
          0.3555685489606980D+00, &
          0.9237817545492393D+00, &
          0.9237817545492393D+00, &
          0.9874241069713281D+00, &
          0.7767138770406730D+00, &
          0.9874241069713281D+00, &
          0.7767138770406730D+00, &
          0.0762182454507607D+00, &
          0.0762182454507607D+00, &
          0.9874241069713281D+00, &
          0.7767138770406730D+00, &
          0.9874241069713281D+00, &
          0.7767138770406730D+00, &
          0.9237817545492393D+00, &
          0.9237817545492393D+00, &
          0.0125758930286720D+00, &
          0.7767138770406730D+00, &
          0.0125758930286720D+00, &
          0.7767138770406730D+00, &
          0.9237817545492393D+00, &
          0.9237817545492393D+00, &
          0.9874241069713281D+00, &
          0.2232861229593271D+00, &
          0.9874241069713281D+00, &
          0.2232861229593271D+00, &
          0.0762182454507607D+00, &
          0.0762182454507607D+00, &
          0.0125758930286720D+00, &
          0.7767138770406730D+00, &
          0.0125758930286720D+00, &
          0.7767138770406730D+00, &
          0.0762182454507607D+00, &
          0.0762182454507607D+00, &
          0.9874241069713281D+00, &
          0.2232861229593271D+00, &
          0.9874241069713281D+00, &
          0.2232861229593271D+00, &
          0.9237817545492393D+00, &
          0.9237817545492393D+00, &
          0.0125758930286720D+00, &
          0.2232861229593271D+00, &
          0.0125758930286720D+00, &
          0.2232861229593271D+00, &
          0.0762182454507607D+00, &
          0.0762182454507607D+00, &
          0.0125758930286720D+00, &
          0.2232861229593271D+00, &
          0.0125758930286720D+00, &
          0.2232861229593271D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.5000000000000000D+00, &
          0.6178975194487687D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.3821024805512313D+00, &
          0.5000000000000000D+00, &
          0.9442475200220319D+00, &
          0.0557524799779681D+00, &
          0.9442475200220319D+00, &
          0.0557524799779681D+00, &
          0.9442475200220319D+00, &
          0.0557524799779681D+00, &
          0.0557524799779681D+00, &
          0.9442475200220319D+00, &
          0.8345909234965025D+00, &
          0.1654090765034974D+00, &
          0.8345909234965025D+00, &
          0.1654090765034974D+00, &
          0.8345909234965025D+00, &
          0.1654090765034974D+00, &
          0.1654090765034974D+00, &
          0.8345909234965025D+00, &
          0.9991488896425311D+00, &
          0.0008511103574689D+00, &
          0.0008511103574689D+00, &
          0.9991488896425311D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.8561550198509273D+00, &
          0.1438449801490727D+00, &
          0.8561550198509273D+00, &
          0.1438449801490727D+00, &
          0.8561550198509273D+00, &
          0.1438449801490727D+00, &
          0.8561550198509273D+00, &
          0.1438449801490727D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9991488896425311D+00, &
          0.0008511103574689D+00, &
          0.0008511103574689D+00, &
          0.9991488896425311D+00, &
          0.5860120647237821D+00, &
          0.4139879352762179D+00, &
          0.4139879352762179D+00, &
          0.5860120647237821D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.6998570581371140D+00, &
          0.3001429418628860D+00, &
          0.6998570581371140D+00, &
          0.3001429418628860D+00, &
          0.6998570581371140D+00, &
          0.3001429418628860D+00, &
          0.6998570581371140D+00, &
          0.3001429418628860D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5860120647237821D+00, &
          0.4139879352762179D+00, &
          0.4139879352762179D+00, &
          0.5860120647237821D+00, &
          0.2840454955532042D+00, &
          0.7159545044467958D+00, &
          0.7159545044467958D+00, &
          0.2840454955532042D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9382902632360840D+00, &
          0.0617097367639160D+00, &
          0.9382902632360840D+00, &
          0.0617097367639160D+00, &
          0.9382902632360840D+00, &
          0.0617097367639160D+00, &
          0.9382902632360840D+00, &
          0.0617097367639160D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.2840454955532042D+00, &
          0.7159545044467958D+00, &
          0.7159545044467958D+00, &
          0.2840454955532042D+00, &
          0.7705521820675761D+00, &
          0.3449134012977228D+00, &
          0.3449134012977228D+00, &
          0.7705521820675761D+00, &
          0.6550865987022771D+00, &
          0.3449134012977228D+00, &
          0.7705521820675761D+00, &
          0.6550865987022771D+00, &
          0.6550865987022771D+00, &
          0.7705521820675761D+00, &
          0.3449134012977228D+00, &
          0.6550865987022771D+00, &
          0.2294478179324239D+00, &
          0.3449134012977228D+00, &
          0.3449134012977228D+00, &
          0.2294478179324239D+00, &
          0.6550865987022771D+00, &
          0.3449134012977228D+00, &
          0.2294478179324239D+00, &
          0.6550865987022771D+00, &
          0.6550865987022771D+00, &
          0.2294478179324239D+00, &
          0.3449134012977228D+00, &
          0.6550865987022771D+00, &
          0.9103324974446185D+00, &
          0.9920926274878630D+00, &
          0.9920926274878630D+00, &
          0.9103324974446185D+00, &
          0.0079073725121370D+00, &
          0.9920926274878630D+00, &
          0.9103324974446185D+00, &
          0.0079073725121370D+00, &
          0.0079073725121370D+00, &
          0.9103324974446185D+00, &
          0.9920926274878630D+00, &
          0.0079073725121370D+00, &
          0.0896675025553814D+00, &
          0.9920926274878630D+00, &
          0.9920926274878630D+00, &
          0.0896675025553814D+00, &
          0.0079073725121370D+00, &
          0.9920926274878630D+00, &
          0.0896675025553814D+00, &
          0.0079073725121370D+00, &
          0.0079073725121370D+00, &
          0.0896675025553814D+00, &
          0.9920926274878630D+00, &
          0.0079073725121370D+00, &
          0.4263116836838957D+00, &
          0.9221936167657927D+00, &
          0.9221936167657927D+00, &
          0.4263116836838957D+00, &
          0.0778063832342073D+00, &
          0.9221936167657927D+00, &
          0.4263116836838957D+00, &
          0.0778063832342073D+00, &
          0.0778063832342073D+00, &
          0.4263116836838957D+00, &
          0.9221936167657927D+00, &
          0.0778063832342073D+00, &
          0.5736883163161043D+00, &
          0.9221936167657927D+00, &
          0.9221936167657927D+00, &
          0.5736883163161043D+00, &
          0.0778063832342073D+00, &
          0.9221936167657927D+00, &
          0.5736883163161043D+00, &
          0.0778063832342073D+00, &
          0.0778063832342073D+00, &
          0.5736883163161043D+00, &
          0.9221936167657927D+00, &
          0.0778063832342073D+00, &
          0.8837457372516819D+00, &
          0.7437996999426106D+00, &
          0.7437996999426106D+00, &
          0.8837457372516819D+00, &
          0.2562003000573894D+00, &
          0.7437996999426106D+00, &
          0.8837457372516819D+00, &
          0.2562003000573894D+00, &
          0.2562003000573894D+00, &
          0.8837457372516819D+00, &
          0.7437996999426106D+00, &
          0.2562003000573894D+00, &
          0.1162542627483181D+00, &
          0.7437996999426106D+00, &
          0.7437996999426106D+00, &
          0.1162542627483181D+00, &
          0.2562003000573894D+00, &
          0.7437996999426106D+00, &
          0.1162542627483181D+00, &
          0.2562003000573894D+00, &
          0.2562003000573894D+00, &
          0.1162542627483181D+00, &
          0.7437996999426106D+00, &
          0.2562003000573894D+00, &
          0.6931955182506035D+00, &
          0.8820659248983722D+00, &
          0.8820659248983722D+00, &
          0.6931955182506035D+00, &
          0.1179340751016278D+00, &
          0.8820659248983722D+00, &
          0.6931955182506035D+00, &
          0.1179340751016278D+00, &
          0.1179340751016278D+00, &
          0.6931955182506035D+00, &
          0.8820659248983722D+00, &
          0.1179340751016278D+00, &
          0.3068044817493964D+00, &
          0.8820659248983722D+00, &
          0.8820659248983722D+00, &
          0.3068044817493964D+00, &
          0.1179340751016278D+00, &
          0.8820659248983722D+00, &
          0.3068044817493964D+00, &
          0.1179340751016278D+00, &
          0.1179340751016278D+00, &
          0.3068044817493964D+00, &
          0.8820659248983722D+00, &
          0.1179340751016278D+00, &
          0.9516006474964566D+00, &
          0.8556268694112152D+00, &
          0.8556268694112152D+00, &
          0.9516006474964566D+00, &
          0.1443731305887847D+00, &
          0.8556268694112152D+00, &
          0.9516006474964566D+00, &
          0.1443731305887847D+00, &
          0.1443731305887847D+00, &
          0.9516006474964566D+00, &
          0.8556268694112152D+00, &
          0.1443731305887847D+00, &
          0.0483993525035434D+00, &
          0.8556268694112152D+00, &
          0.8556268694112152D+00, &
          0.0483993525035434D+00, &
          0.1443731305887847D+00, &
          0.8556268694112152D+00, &
          0.0483993525035434D+00, &
          0.1443731305887847D+00, &
          0.1443731305887847D+00, &
          0.0483993525035434D+00, &
          0.8556268694112152D+00, &
          0.1443731305887847D+00, &
          0.6113080322285057D+00, &
          0.9775744747969544D+00, &
          0.9775744747969544D+00, &
          0.6113080322285057D+00, &
          0.0224255252030456D+00, &
          0.9775744747969544D+00, &
          0.6113080322285057D+00, &
          0.0224255252030456D+00, &
          0.0224255252030456D+00, &
          0.6113080322285057D+00, &
          0.9775744747969544D+00, &
          0.0224255252030456D+00, &
          0.3886919677714943D+00, &
          0.9775744747969544D+00, &
          0.9775744747969544D+00, &
          0.3886919677714943D+00, &
          0.0224255252030456D+00, &
          0.9775744747969544D+00, &
          0.3886919677714943D+00, &
          0.0224255252030456D+00, &
          0.0224255252030456D+00, &
          0.3886919677714943D+00, &
          0.9775744747969544D+00, &
          0.0224255252030456D+00, &
          0.9825814337796190D+00, &
          0.5848459901934986D+00, &
          0.5848459901934986D+00, &
          0.9825814337796190D+00, &
          0.4151540098065015D+00, &
          0.5848459901934986D+00, &
          0.9825814337796190D+00, &
          0.4151540098065015D+00, &
          0.4151540098065015D+00, &
          0.9825814337796190D+00, &
          0.5848459901934986D+00, &
          0.4151540098065015D+00, &
          0.0174185662203810D+00, &
          0.5848459901934986D+00, &
          0.5848459901934986D+00, &
          0.0174185662203810D+00, &
          0.4151540098065015D+00, &
          0.5848459901934986D+00, &
          0.0174185662203810D+00, &
          0.4151540098065015D+00, &
          0.4151540098065015D+00, &
          0.0174185662203810D+00, &
          0.5848459901934986D+00, &
          0.4151540098065015D+00, &
          0.5470636340177077D+00, &
          0.8042667629417029D+00, &
          0.8042667629417029D+00, &
          0.5470636340177077D+00, &
          0.1957332370582971D+00, &
          0.8042667629417029D+00, &
          0.5470636340177077D+00, &
          0.1957332370582971D+00, &
          0.1957332370582971D+00, &
          0.5470636340177077D+00, &
          0.8042667629417029D+00, &
          0.1957332370582971D+00, &
          0.4529363659822923D+00, &
          0.8042667629417029D+00, &
          0.8042667629417029D+00, &
          0.4529363659822923D+00, &
          0.1957332370582971D+00, &
          0.8042667629417029D+00, &
          0.4529363659822923D+00, &
          0.1957332370582971D+00, &
          0.1957332370582971D+00, &
          0.4529363659822923D+00, &
          0.8042667629417029D+00, &
          0.1957332370582971D+00, &
          0.9892390565548020D+00, &
          0.7132758893843609D+00, &
          0.7132758893843609D+00, &
          0.9892390565548020D+00, &
          0.2867241106156391D+00, &
          0.7132758893843609D+00, &
          0.9892390565548020D+00, &
          0.2867241106156391D+00, &
          0.2867241106156391D+00, &
          0.9892390565548020D+00, &
          0.7132758893843609D+00, &
          0.2867241106156391D+00, &
          0.0107609434451980D+00, &
          0.7132758893843609D+00, &
          0.7132758893843609D+00, &
          0.0107609434451980D+00, &
          0.2867241106156391D+00, &
          0.7132758893843609D+00, &
          0.0107609434451980D+00, &
          0.2867241106156391D+00, &
          0.2867241106156391D+00, &
          0.0107609434451980D+00, &
          0.7132758893843609D+00, &
          0.2867241106156391D+00, &
          0.8728877113378553D+00, &
          0.5918460351304471D+00, &
          0.5918460351304471D+00, &
          0.8728877113378553D+00, &
          0.4081539648695529D+00, &
          0.5918460351304471D+00, &
          0.8728877113378553D+00, &
          0.4081539648695529D+00, &
          0.4081539648695529D+00, &
          0.8728877113378553D+00, &
          0.5918460351304471D+00, &
          0.4081539648695529D+00, &
          0.1271122886621447D+00, &
          0.5918460351304471D+00, &
          0.5918460351304471D+00, &
          0.1271122886621447D+00, &
          0.4081539648695529D+00, &
          0.5918460351304471D+00, &
          0.1271122886621447D+00, &
          0.4081539648695529D+00, &
          0.4081539648695529D+00, &
          0.1271122886621447D+00, &
          0.5918460351304471D+00, &
          0.4081539648695529D+00, &
          0.9566761563025079D+00, &
          0.6444314510393020D+00, &
          0.8267665443298317D+00, &
          0.8267665443298317D+00, &
          0.6444314510393020D+00, &
          0.9566761563025079D+00, &
          0.9566761563025079D+00, &
          0.6444314510393020D+00, &
          0.1732334556701683D+00, &
          0.1732334556701683D+00, &
          0.6444314510393020D+00, &
          0.9566761563025079D+00, &
          0.0433238436974921D+00, &
          0.6444314510393020D+00, &
          0.8267665443298317D+00, &
          0.8267665443298317D+00, &
          0.6444314510393020D+00, &
          0.0433238436974921D+00, &
          0.9566761563025079D+00, &
          0.3555685489606980D+00, &
          0.8267665443298317D+00, &
          0.8267665443298317D+00, &
          0.3555685489606980D+00, &
          0.9566761563025079D+00, &
          0.0433238436974921D+00, &
          0.6444314510393020D+00, &
          0.1732334556701683D+00, &
          0.1732334556701683D+00, &
          0.6444314510393020D+00, &
          0.0433238436974921D+00, &
          0.9566761563025079D+00, &
          0.3555685489606980D+00, &
          0.1732334556701683D+00, &
          0.1732334556701683D+00, &
          0.3555685489606980D+00, &
          0.9566761563025079D+00, &
          0.0433238436974921D+00, &
          0.3555685489606980D+00, &
          0.8267665443298317D+00, &
          0.8267665443298317D+00, &
          0.3555685489606980D+00, &
          0.0433238436974921D+00, &
          0.0433238436974921D+00, &
          0.3555685489606980D+00, &
          0.1732334556701683D+00, &
          0.1732334556701683D+00, &
          0.3555685489606980D+00, &
          0.0433238436974921D+00, &
          0.9874241069713281D+00, &
          0.7767138770406730D+00, &
          0.9237817545492393D+00, &
          0.9237817545492393D+00, &
          0.7767138770406730D+00, &
          0.9874241069713281D+00, &
          0.9874241069713281D+00, &
          0.7767138770406730D+00, &
          0.0762182454507607D+00, &
          0.0762182454507607D+00, &
          0.7767138770406730D+00, &
          0.9874241069713281D+00, &
          0.0125758930286720D+00, &
          0.7767138770406730D+00, &
          0.9237817545492393D+00, &
          0.9237817545492393D+00, &
          0.7767138770406730D+00, &
          0.0125758930286720D+00, &
          0.9874241069713281D+00, &
          0.2232861229593271D+00, &
          0.9237817545492393D+00, &
          0.9237817545492393D+00, &
          0.2232861229593271D+00, &
          0.9874241069713281D+00, &
          0.0125758930286720D+00, &
          0.7767138770406730D+00, &
          0.0762182454507607D+00, &
          0.0762182454507607D+00, &
          0.7767138770406730D+00, &
          0.0125758930286720D+00, &
          0.9874241069713281D+00, &
          0.2232861229593271D+00, &
          0.0762182454507607D+00, &
          0.0762182454507607D+00, &
          0.2232861229593271D+00, &
          0.9874241069713281D+00, &
          0.0125758930286720D+00, &
          0.2232861229593271D+00, &
          0.9237817545492393D+00, &
          0.9237817545492393D+00, &
          0.2232861229593271D+00, &
          0.0125758930286720D+00, &
          0.0125758930286720D+00, &
          0.2232861229593271D+00, &
          0.0762182454507607D+00, &
          0.0762182454507607D+00, &
          0.2232861229593271D+00, &
          0.0125758930286720D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.6178975194487687D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.3821024805512313D+00, &
          0.9442475200220319D+00, &
          0.0557524799779681D+00, &
          0.9442475200220319D+00, &
          0.9442475200220319D+00, &
          0.0557524799779681D+00, &
          0.9442475200220319D+00, &
          0.0557524799779681D+00, &
          0.0557524799779681D+00, &
          0.8345909234965025D+00, &
          0.1654090765034974D+00, &
          0.8345909234965025D+00, &
          0.8345909234965025D+00, &
          0.1654090765034974D+00, &
          0.8345909234965025D+00, &
          0.1654090765034974D+00, &
          0.1654090765034974D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9991488896425311D+00, &
          0.0008511103574689D+00, &
          0.0008511103574689D+00, &
          0.9991488896425311D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9991488896425311D+00, &
          0.0008511103574689D+00, &
          0.0008511103574689D+00, &
          0.9991488896425311D+00, &
          0.8561550198509273D+00, &
          0.1438449801490727D+00, &
          0.8561550198509273D+00, &
          0.1438449801490727D+00, &
          0.8561550198509273D+00, &
          0.1438449801490727D+00, &
          0.8561550198509273D+00, &
          0.1438449801490727D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5860120647237821D+00, &
          0.4139879352762179D+00, &
          0.4139879352762179D+00, &
          0.5860120647237821D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5860120647237821D+00, &
          0.4139879352762179D+00, &
          0.4139879352762179D+00, &
          0.5860120647237821D+00, &
          0.6998570581371140D+00, &
          0.3001429418628860D+00, &
          0.6998570581371140D+00, &
          0.3001429418628860D+00, &
          0.6998570581371140D+00, &
          0.3001429418628860D+00, &
          0.6998570581371140D+00, &
          0.3001429418628860D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.2840454955532042D+00, &
          0.7159545044467958D+00, &
          0.7159545044467958D+00, &
          0.2840454955532042D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.2840454955532042D+00, &
          0.7159545044467958D+00, &
          0.7159545044467958D+00, &
          0.2840454955532042D+00, &
          0.9382902632360840D+00, &
          0.0617097367639160D+00, &
          0.9382902632360840D+00, &
          0.0617097367639160D+00, &
          0.9382902632360840D+00, &
          0.0617097367639160D+00, &
          0.9382902632360840D+00, &
          0.0617097367639160D+00, &
          0.3449134012977228D+00, &
          0.3449134012977228D+00, &
          0.7705521820675761D+00, &
          0.3449134012977228D+00, &
          0.3449134012977228D+00, &
          0.7705521820675761D+00, &
          0.6550865987022771D+00, &
          0.6550865987022771D+00, &
          0.7705521820675761D+00, &
          0.6550865987022771D+00, &
          0.6550865987022771D+00, &
          0.7705521820675761D+00, &
          0.3449134012977228D+00, &
          0.3449134012977228D+00, &
          0.2294478179324239D+00, &
          0.3449134012977228D+00, &
          0.3449134012977228D+00, &
          0.2294478179324239D+00, &
          0.6550865987022771D+00, &
          0.6550865987022771D+00, &
          0.2294478179324239D+00, &
          0.6550865987022771D+00, &
          0.6550865987022771D+00, &
          0.2294478179324239D+00, &
          0.9920926274878630D+00, &
          0.9920926274878630D+00, &
          0.9103324974446185D+00, &
          0.9920926274878630D+00, &
          0.9920926274878630D+00, &
          0.9103324974446185D+00, &
          0.0079073725121370D+00, &
          0.0079073725121370D+00, &
          0.9103324974446185D+00, &
          0.0079073725121370D+00, &
          0.0079073725121370D+00, &
          0.9103324974446185D+00, &
          0.9920926274878630D+00, &
          0.9920926274878630D+00, &
          0.0896675025553814D+00, &
          0.9920926274878630D+00, &
          0.9920926274878630D+00, &
          0.0896675025553814D+00, &
          0.0079073725121370D+00, &
          0.0079073725121370D+00, &
          0.0896675025553814D+00, &
          0.0079073725121370D+00, &
          0.0079073725121370D+00, &
          0.0896675025553814D+00, &
          0.9221936167657927D+00, &
          0.9221936167657927D+00, &
          0.4263116836838957D+00, &
          0.9221936167657927D+00, &
          0.9221936167657927D+00, &
          0.4263116836838957D+00, &
          0.0778063832342073D+00, &
          0.0778063832342073D+00, &
          0.4263116836838957D+00, &
          0.0778063832342073D+00, &
          0.0778063832342073D+00, &
          0.4263116836838957D+00, &
          0.9221936167657927D+00, &
          0.9221936167657927D+00, &
          0.5736883163161043D+00, &
          0.9221936167657927D+00, &
          0.9221936167657927D+00, &
          0.5736883163161043D+00, &
          0.0778063832342073D+00, &
          0.0778063832342073D+00, &
          0.5736883163161043D+00, &
          0.0778063832342073D+00, &
          0.0778063832342073D+00, &
          0.5736883163161043D+00, &
          0.7437996999426106D+00, &
          0.7437996999426106D+00, &
          0.8837457372516819D+00, &
          0.7437996999426106D+00, &
          0.7437996999426106D+00, &
          0.8837457372516819D+00, &
          0.2562003000573894D+00, &
          0.2562003000573894D+00, &
          0.8837457372516819D+00, &
          0.2562003000573894D+00, &
          0.2562003000573894D+00, &
          0.8837457372516819D+00, &
          0.7437996999426106D+00, &
          0.7437996999426106D+00, &
          0.1162542627483181D+00, &
          0.7437996999426106D+00, &
          0.7437996999426106D+00, &
          0.1162542627483181D+00, &
          0.2562003000573894D+00, &
          0.2562003000573894D+00, &
          0.1162542627483181D+00, &
          0.2562003000573894D+00, &
          0.2562003000573894D+00, &
          0.1162542627483181D+00, &
          0.8820659248983722D+00, &
          0.8820659248983722D+00, &
          0.6931955182506035D+00, &
          0.8820659248983722D+00, &
          0.8820659248983722D+00, &
          0.6931955182506035D+00, &
          0.1179340751016278D+00, &
          0.1179340751016278D+00, &
          0.6931955182506035D+00, &
          0.1179340751016278D+00, &
          0.1179340751016278D+00, &
          0.6931955182506035D+00, &
          0.8820659248983722D+00, &
          0.8820659248983722D+00, &
          0.3068044817493964D+00, &
          0.8820659248983722D+00, &
          0.8820659248983722D+00, &
          0.3068044817493964D+00, &
          0.1179340751016278D+00, &
          0.1179340751016278D+00, &
          0.3068044817493964D+00, &
          0.1179340751016278D+00, &
          0.1179340751016278D+00, &
          0.3068044817493964D+00, &
          0.8556268694112152D+00, &
          0.8556268694112152D+00, &
          0.9516006474964566D+00, &
          0.8556268694112152D+00, &
          0.8556268694112152D+00, &
          0.9516006474964566D+00, &
          0.1443731305887847D+00, &
          0.1443731305887847D+00, &
          0.9516006474964566D+00, &
          0.1443731305887847D+00, &
          0.1443731305887847D+00, &
          0.9516006474964566D+00, &
          0.8556268694112152D+00, &
          0.8556268694112152D+00, &
          0.0483993525035434D+00, &
          0.8556268694112152D+00, &
          0.8556268694112152D+00, &
          0.0483993525035434D+00, &
          0.1443731305887847D+00, &
          0.1443731305887847D+00, &
          0.0483993525035434D+00, &
          0.1443731305887847D+00, &
          0.1443731305887847D+00, &
          0.0483993525035434D+00, &
          0.9775744747969544D+00, &
          0.9775744747969544D+00, &
          0.6113080322285057D+00, &
          0.9775744747969544D+00, &
          0.9775744747969544D+00, &
          0.6113080322285057D+00, &
          0.0224255252030456D+00, &
          0.0224255252030456D+00, &
          0.6113080322285057D+00, &
          0.0224255252030456D+00, &
          0.0224255252030456D+00, &
          0.6113080322285057D+00, &
          0.9775744747969544D+00, &
          0.9775744747969544D+00, &
          0.3886919677714943D+00, &
          0.9775744747969544D+00, &
          0.9775744747969544D+00, &
          0.3886919677714943D+00, &
          0.0224255252030456D+00, &
          0.0224255252030456D+00, &
          0.3886919677714943D+00, &
          0.0224255252030456D+00, &
          0.0224255252030456D+00, &
          0.3886919677714943D+00, &
          0.5848459901934986D+00, &
          0.5848459901934986D+00, &
          0.9825814337796190D+00, &
          0.5848459901934986D+00, &
          0.5848459901934986D+00, &
          0.9825814337796190D+00, &
          0.4151540098065015D+00, &
          0.4151540098065015D+00, &
          0.9825814337796190D+00, &
          0.4151540098065015D+00, &
          0.4151540098065015D+00, &
          0.9825814337796190D+00, &
          0.5848459901934986D+00, &
          0.5848459901934986D+00, &
          0.0174185662203810D+00, &
          0.5848459901934986D+00, &
          0.5848459901934986D+00, &
          0.0174185662203810D+00, &
          0.4151540098065015D+00, &
          0.4151540098065015D+00, &
          0.0174185662203810D+00, &
          0.4151540098065015D+00, &
          0.4151540098065015D+00, &
          0.0174185662203810D+00, &
          0.8042667629417029D+00, &
          0.8042667629417029D+00, &
          0.5470636340177077D+00, &
          0.8042667629417029D+00, &
          0.8042667629417029D+00, &
          0.5470636340177077D+00, &
          0.1957332370582971D+00, &
          0.1957332370582971D+00, &
          0.5470636340177077D+00, &
          0.1957332370582971D+00, &
          0.1957332370582971D+00, &
          0.5470636340177077D+00, &
          0.8042667629417029D+00, &
          0.8042667629417029D+00, &
          0.4529363659822923D+00, &
          0.8042667629417029D+00, &
          0.8042667629417029D+00, &
          0.4529363659822923D+00, &
          0.1957332370582971D+00, &
          0.1957332370582971D+00, &
          0.4529363659822923D+00, &
          0.1957332370582971D+00, &
          0.1957332370582971D+00, &
          0.4529363659822923D+00, &
          0.7132758893843609D+00, &
          0.7132758893843609D+00, &
          0.9892390565548020D+00, &
          0.7132758893843609D+00, &
          0.7132758893843609D+00, &
          0.9892390565548020D+00, &
          0.2867241106156391D+00, &
          0.2867241106156391D+00, &
          0.9892390565548020D+00, &
          0.2867241106156391D+00, &
          0.2867241106156391D+00, &
          0.9892390565548020D+00, &
          0.7132758893843609D+00, &
          0.7132758893843609D+00, &
          0.0107609434451980D+00, &
          0.7132758893843609D+00, &
          0.7132758893843609D+00, &
          0.0107609434451980D+00, &
          0.2867241106156391D+00, &
          0.2867241106156391D+00, &
          0.0107609434451980D+00, &
          0.2867241106156391D+00, &
          0.2867241106156391D+00, &
          0.0107609434451980D+00, &
          0.5918460351304471D+00, &
          0.5918460351304471D+00, &
          0.8728877113378553D+00, &
          0.5918460351304471D+00, &
          0.5918460351304471D+00, &
          0.8728877113378553D+00, &
          0.4081539648695529D+00, &
          0.4081539648695529D+00, &
          0.8728877113378553D+00, &
          0.4081539648695529D+00, &
          0.4081539648695529D+00, &
          0.8728877113378553D+00, &
          0.5918460351304471D+00, &
          0.5918460351304471D+00, &
          0.1271122886621447D+00, &
          0.5918460351304471D+00, &
          0.5918460351304471D+00, &
          0.1271122886621447D+00, &
          0.4081539648695529D+00, &
          0.4081539648695529D+00, &
          0.1271122886621447D+00, &
          0.4081539648695529D+00, &
          0.4081539648695529D+00, &
          0.1271122886621447D+00, &
          0.6444314510393020D+00, &
          0.9566761563025079D+00, &
          0.6444314510393020D+00, &
          0.9566761563025079D+00, &
          0.8267665443298317D+00, &
          0.8267665443298317D+00, &
          0.6444314510393020D+00, &
          0.9566761563025079D+00, &
          0.6444314510393020D+00, &
          0.9566761563025079D+00, &
          0.1732334556701683D+00, &
          0.1732334556701683D+00, &
          0.6444314510393020D+00, &
          0.0433238436974921D+00, &
          0.6444314510393020D+00, &
          0.0433238436974921D+00, &
          0.8267665443298317D+00, &
          0.8267665443298317D+00, &
          0.3555685489606980D+00, &
          0.9566761563025079D+00, &
          0.3555685489606980D+00, &
          0.9566761563025079D+00, &
          0.8267665443298317D+00, &
          0.8267665443298317D+00, &
          0.6444314510393020D+00, &
          0.0433238436974921D+00, &
          0.6444314510393020D+00, &
          0.0433238436974921D+00, &
          0.1732334556701683D+00, &
          0.1732334556701683D+00, &
          0.3555685489606980D+00, &
          0.9566761563025079D+00, &
          0.3555685489606980D+00, &
          0.9566761563025079D+00, &
          0.1732334556701683D+00, &
          0.1732334556701683D+00, &
          0.3555685489606980D+00, &
          0.0433238436974921D+00, &
          0.3555685489606980D+00, &
          0.0433238436974921D+00, &
          0.8267665443298317D+00, &
          0.8267665443298317D+00, &
          0.3555685489606980D+00, &
          0.0433238436974921D+00, &
          0.3555685489606980D+00, &
          0.0433238436974921D+00, &
          0.1732334556701683D+00, &
          0.1732334556701683D+00, &
          0.7767138770406730D+00, &
          0.9874241069713281D+00, &
          0.7767138770406730D+00, &
          0.9874241069713281D+00, &
          0.9237817545492393D+00, &
          0.9237817545492393D+00, &
          0.7767138770406730D+00, &
          0.9874241069713281D+00, &
          0.7767138770406730D+00, &
          0.9874241069713281D+00, &
          0.0762182454507607D+00, &
          0.0762182454507607D+00, &
          0.7767138770406730D+00, &
          0.0125758930286720D+00, &
          0.7767138770406730D+00, &
          0.0125758930286720D+00, &
          0.9237817545492393D+00, &
          0.9237817545492393D+00, &
          0.2232861229593271D+00, &
          0.9874241069713281D+00, &
          0.2232861229593271D+00, &
          0.9874241069713281D+00, &
          0.9237817545492393D+00, &
          0.9237817545492393D+00, &
          0.7767138770406730D+00, &
          0.0125758930286720D+00, &
          0.7767138770406730D+00, &
          0.0125758930286720D+00, &
          0.0762182454507607D+00, &
          0.0762182454507607D+00, &
          0.2232861229593271D+00, &
          0.9874241069713281D+00, &
          0.2232861229593271D+00, &
          0.9874241069713281D+00, &
          0.0762182454507607D+00, &
          0.0762182454507607D+00, &
          0.2232861229593271D+00, &
          0.0125758930286720D+00, &
          0.2232861229593271D+00, &
          0.0125758930286720D+00, &
          0.9237817545492393D+00, &
          0.9237817545492393D+00, &
          0.2232861229593271D+00, &
          0.0125758930286720D+00, &
          0.2232861229593271D+00, &
          0.0125758930286720D+00, &
          0.0762182454507607D+00, &
          0.0762182454507607D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          0.0020573710057311D+00, &
          0.0020573710057311D+00, &
          0.0020573710057311D+00, &
          0.0020573710057311D+00, &
          0.0020573710057311D+00, &
          0.0020573710057311D+00, &
          0.0009660877915728D+00, &
          0.0009660877915728D+00, &
          0.0009660877915728D+00, &
          0.0009660877915728D+00, &
          0.0009660877915728D+00, &
          0.0009660877915728D+00, &
          0.0009660877915728D+00, &
          0.0009660877915728D+00, &
          0.0016359014467992D+00, &
          0.0016359014467992D+00, &
          0.0016359014467992D+00, &
          0.0016359014467992D+00, &
          0.0016359014467992D+00, &
          0.0016359014467992D+00, &
          0.0016359014467992D+00, &
          0.0016359014467992D+00, &
          0.0008499912813105D+00, &
          0.0008499912813105D+00, &
          0.0008499912813105D+00, &
          0.0008499912813105D+00, &
          0.0008499912813105D+00, &
          0.0008499912813105D+00, &
          0.0008499912813105D+00, &
          0.0008499912813105D+00, &
          0.0008499912813105D+00, &
          0.0008499912813105D+00, &
          0.0008499912813105D+00, &
          0.0008499912813105D+00, &
          0.0008499912813105D+00, &
          0.0008499912813105D+00, &
          0.0008499912813105D+00, &
          0.0008499912813105D+00, &
          0.0008499912813105D+00, &
          0.0008499912813105D+00, &
          0.0008499912813105D+00, &
          0.0008499912813105D+00, &
          0.0008499912813105D+00, &
          0.0008499912813105D+00, &
          0.0008499912813105D+00, &
          0.0008499912813105D+00, &
          0.0035950037217364D+00, &
          0.0035950037217364D+00, &
          0.0035950037217364D+00, &
          0.0035950037217364D+00, &
          0.0035950037217364D+00, &
          0.0035950037217364D+00, &
          0.0035950037217364D+00, &
          0.0035950037217364D+00, &
          0.0035950037217364D+00, &
          0.0035950037217364D+00, &
          0.0035950037217364D+00, &
          0.0035950037217364D+00, &
          0.0035950037217364D+00, &
          0.0035950037217364D+00, &
          0.0035950037217364D+00, &
          0.0035950037217364D+00, &
          0.0035950037217364D+00, &
          0.0035950037217364D+00, &
          0.0035950037217364D+00, &
          0.0035950037217364D+00, &
          0.0035950037217364D+00, &
          0.0035950037217364D+00, &
          0.0035950037217364D+00, &
          0.0035950037217364D+00, &
          0.0038360634296226D+00, &
          0.0038360634296226D+00, &
          0.0038360634296226D+00, &
          0.0038360634296226D+00, &
          0.0038360634296226D+00, &
          0.0038360634296226D+00, &
          0.0038360634296226D+00, &
          0.0038360634296226D+00, &
          0.0038360634296226D+00, &
          0.0038360634296226D+00, &
          0.0038360634296226D+00, &
          0.0038360634296226D+00, &
          0.0038360634296226D+00, &
          0.0038360634296226D+00, &
          0.0038360634296226D+00, &
          0.0038360634296226D+00, &
          0.0038360634296226D+00, &
          0.0038360634296226D+00, &
          0.0038360634296226D+00, &
          0.0038360634296226D+00, &
          0.0038360634296226D+00, &
          0.0038360634296226D+00, &
          0.0038360634296226D+00, &
          0.0038360634296226D+00, &
          0.0050441284857989D+00, &
          0.0050441284857989D+00, &
          0.0050441284857989D+00, &
          0.0050441284857989D+00, &
          0.0050441284857989D+00, &
          0.0050441284857989D+00, &
          0.0050441284857989D+00, &
          0.0050441284857989D+00, &
          0.0050441284857989D+00, &
          0.0050441284857989D+00, &
          0.0050441284857989D+00, &
          0.0050441284857989D+00, &
          0.0050441284857989D+00, &
          0.0050441284857989D+00, &
          0.0050441284857989D+00, &
          0.0050441284857989D+00, &
          0.0050441284857989D+00, &
          0.0050441284857989D+00, &
          0.0050441284857989D+00, &
          0.0050441284857989D+00, &
          0.0050441284857989D+00, &
          0.0050441284857989D+00, &
          0.0050441284857989D+00, &
          0.0050441284857989D+00, &
          0.0001951830532879D+00, &
          0.0001951830532879D+00, &
          0.0001951830532879D+00, &
          0.0001951830532879D+00, &
          0.0001951830532879D+00, &
          0.0001951830532879D+00, &
          0.0001951830532879D+00, &
          0.0001951830532879D+00, &
          0.0001951830532879D+00, &
          0.0001951830532879D+00, &
          0.0001951830532879D+00, &
          0.0001951830532879D+00, &
          0.0001951830532879D+00, &
          0.0001951830532879D+00, &
          0.0001951830532879D+00, &
          0.0001951830532879D+00, &
          0.0001951830532879D+00, &
          0.0001951830532879D+00, &
          0.0001951830532879D+00, &
          0.0001951830532879D+00, &
          0.0001951830532879D+00, &
          0.0001951830532879D+00, &
          0.0001951830532879D+00, &
          0.0001951830532879D+00, &
          0.0017584548406612D+00, &
          0.0017584548406612D+00, &
          0.0017584548406612D+00, &
          0.0017584548406612D+00, &
          0.0017584548406612D+00, &
          0.0017584548406612D+00, &
          0.0017584548406612D+00, &
          0.0017584548406612D+00, &
          0.0017584548406612D+00, &
          0.0017584548406612D+00, &
          0.0017584548406612D+00, &
          0.0017584548406612D+00, &
          0.0017584548406612D+00, &
          0.0017584548406612D+00, &
          0.0017584548406612D+00, &
          0.0017584548406612D+00, &
          0.0017584548406612D+00, &
          0.0017584548406612D+00, &
          0.0017584548406612D+00, &
          0.0017584548406612D+00, &
          0.0017584548406612D+00, &
          0.0017584548406612D+00, &
          0.0017584548406612D+00, &
          0.0017584548406612D+00, &
          0.0038873202629322D+00, &
          0.0038873202629322D+00, &
          0.0038873202629322D+00, &
          0.0038873202629322D+00, &
          0.0038873202629322D+00, &
          0.0038873202629322D+00, &
          0.0038873202629322D+00, &
          0.0038873202629322D+00, &
          0.0038873202629322D+00, &
          0.0038873202629322D+00, &
          0.0038873202629322D+00, &
          0.0038873202629322D+00, &
          0.0038873202629322D+00, &
          0.0038873202629322D+00, &
          0.0038873202629322D+00, &
          0.0038873202629322D+00, &
          0.0038873202629322D+00, &
          0.0038873202629322D+00, &
          0.0038873202629322D+00, &
          0.0038873202629322D+00, &
          0.0038873202629322D+00, &
          0.0038873202629322D+00, &
          0.0038873202629322D+00, &
          0.0038873202629322D+00, &
          0.0018874010836757D+00, &
          0.0018874010836757D+00, &
          0.0018874010836757D+00, &
          0.0018874010836757D+00, &
          0.0018874010836757D+00, &
          0.0018874010836757D+00, &
          0.0018874010836757D+00, &
          0.0018874010836757D+00, &
          0.0018874010836757D+00, &
          0.0018874010836757D+00, &
          0.0018874010836757D+00, &
          0.0018874010836757D+00, &
          0.0018874010836757D+00, &
          0.0018874010836757D+00, &
          0.0018874010836757D+00, &
          0.0018874010836757D+00, &
          0.0018874010836757D+00, &
          0.0018874010836757D+00, &
          0.0018874010836757D+00, &
          0.0018874010836757D+00, &
          0.0018874010836757D+00, &
          0.0018874010836757D+00, &
          0.0018874010836757D+00, &
          0.0018874010836757D+00, &
          0.0016991651835215D+00, &
          0.0016991651835215D+00, &
          0.0016991651835215D+00, &
          0.0016991651835215D+00, &
          0.0016991651835215D+00, &
          0.0016991651835215D+00, &
          0.0016991651835215D+00, &
          0.0016991651835215D+00, &
          0.0016991651835215D+00, &
          0.0016991651835215D+00, &
          0.0016991651835215D+00, &
          0.0016991651835215D+00, &
          0.0016991651835215D+00, &
          0.0016991651835215D+00, &
          0.0016991651835215D+00, &
          0.0016991651835215D+00, &
          0.0016991651835215D+00, &
          0.0016991651835215D+00, &
          0.0016991651835215D+00, &
          0.0016991651835215D+00, &
          0.0016991651835215D+00, &
          0.0016991651835215D+00, &
          0.0016991651835215D+00, &
          0.0016991651835215D+00, &
          0.0007588378288378D+00, &
          0.0007588378288378D+00, &
          0.0007588378288378D+00, &
          0.0007588378288378D+00, &
          0.0007588378288378D+00, &
          0.0007588378288378D+00, &
          0.0007588378288378D+00, &
          0.0007588378288378D+00, &
          0.0007588378288378D+00, &
          0.0007588378288378D+00, &
          0.0007588378288378D+00, &
          0.0007588378288378D+00, &
          0.0007588378288378D+00, &
          0.0007588378288378D+00, &
          0.0007588378288378D+00, &
          0.0007588378288378D+00, &
          0.0007588378288378D+00, &
          0.0007588378288378D+00, &
          0.0007588378288378D+00, &
          0.0007588378288378D+00, &
          0.0007588378288378D+00, &
          0.0007588378288378D+00, &
          0.0007588378288378D+00, &
          0.0007588378288378D+00, &
          0.0014557629541032D+00, &
          0.0014557629541032D+00, &
          0.0014557629541032D+00, &
          0.0014557629541032D+00, &
          0.0014557629541032D+00, &
          0.0014557629541032D+00, &
          0.0014557629541032D+00, &
          0.0014557629541032D+00, &
          0.0014557629541032D+00, &
          0.0014557629541032D+00, &
          0.0014557629541032D+00, &
          0.0014557629541032D+00, &
          0.0014557629541032D+00, &
          0.0014557629541032D+00, &
          0.0014557629541032D+00, &
          0.0014557629541032D+00, &
          0.0014557629541032D+00, &
          0.0014557629541032D+00, &
          0.0014557629541032D+00, &
          0.0014557629541032D+00, &
          0.0014557629541032D+00, &
          0.0014557629541032D+00, &
          0.0014557629541032D+00, &
          0.0014557629541032D+00, &
          0.0040079050757571D+00, &
          0.0040079050757571D+00, &
          0.0040079050757571D+00, &
          0.0040079050757571D+00, &
          0.0040079050757571D+00, &
          0.0040079050757571D+00, &
          0.0040079050757571D+00, &
          0.0040079050757571D+00, &
          0.0040079050757571D+00, &
          0.0040079050757571D+00, &
          0.0040079050757571D+00, &
          0.0040079050757571D+00, &
          0.0040079050757571D+00, &
          0.0040079050757571D+00, &
          0.0040079050757571D+00, &
          0.0040079050757571D+00, &
          0.0040079050757571D+00, &
          0.0040079050757571D+00, &
          0.0040079050757571D+00, &
          0.0040079050757571D+00, &
          0.0040079050757571D+00, &
          0.0040079050757571D+00, &
          0.0040079050757571D+00, &
          0.0040079050757571D+00, &
          0.0013991860755843D+00, &
          0.0013991860755843D+00, &
          0.0013991860755843D+00, &
          0.0013991860755843D+00, &
          0.0013991860755843D+00, &
          0.0013991860755843D+00, &
          0.0013991860755843D+00, &
          0.0013991860755843D+00, &
          0.0013991860755843D+00, &
          0.0013991860755843D+00, &
          0.0013991860755843D+00, &
          0.0013991860755843D+00, &
          0.0013991860755843D+00, &
          0.0013991860755843D+00, &
          0.0013991860755843D+00, &
          0.0013991860755843D+00, &
          0.0013991860755843D+00, &
          0.0013991860755843D+00, &
          0.0013991860755843D+00, &
          0.0013991860755843D+00, &
          0.0013991860755843D+00, &
          0.0013991860755843D+00, &
          0.0013991860755843D+00, &
          0.0013991860755843D+00, &
          0.0040992807889267D+00, &
          0.0040992807889267D+00, &
          0.0040992807889267D+00, &
          0.0040992807889267D+00, &
          0.0040992807889267D+00, &
          0.0040992807889267D+00, &
          0.0040992807889267D+00, &
          0.0040992807889267D+00, &
          0.0040992807889267D+00, &
          0.0040992807889267D+00, &
          0.0040992807889267D+00, &
          0.0040992807889267D+00, &
          0.0040992807889267D+00, &
          0.0040992807889267D+00, &
          0.0040992807889267D+00, &
          0.0040992807889267D+00, &
          0.0040992807889267D+00, &
          0.0040992807889267D+00, &
          0.0040992807889267D+00, &
          0.0040992807889267D+00, &
          0.0040992807889267D+00, &
          0.0040992807889267D+00, &
          0.0040992807889267D+00, &
          0.0040992807889267D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0019518612573850D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00, &
          0.0009537937942919D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule21 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule21() returns the hexahedron quadrature rule of precision 21.
!
!  Discussion:
!
!    We suppose we are given a hexahedron H with vertices
!    (0,0,0), (1,0,0), (1,1,0), (0,1,0), (0,0,1), (1,0,1), (1,1,1), (0,1,1),
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over H is 
!    approximated by Q as follows:
!
!    Q = volume(H) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 May 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order symmetric cubature rules for tetrahedra and pyramids,
!    International Journal for Numerical Methods in Engineering,
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
!    real w(n): the quadrature weights, which add to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 580

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
          0.1912521437383576D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.8087478562616424D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.3543162299281049D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.6456837700718951D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9311291591147122D+00, &
          0.0688708408852879D+00, &
          0.0688708408852879D+00, &
          0.9311291591147122D+00, &
          0.9311291591147122D+00, &
          0.0688708408852879D+00, &
          0.9311291591147122D+00, &
          0.0688708408852879D+00, &
          0.7833399058047521D+00, &
          0.2166600941952480D+00, &
          0.2166600941952480D+00, &
          0.7833399058047521D+00, &
          0.7833399058047521D+00, &
          0.2166600941952480D+00, &
          0.7833399058047521D+00, &
          0.2166600941952480D+00, &
          0.8981144288438190D+00, &
          0.1018855711561810D+00, &
          0.8981144288438190D+00, &
          0.1018855711561810D+00, &
          0.8981144288438190D+00, &
          0.1018855711561810D+00, &
          0.8981144288438190D+00, &
          0.1018855711561810D+00, &
          0.9999999935557286D+00, &
          0.0000000064442714D+00, &
          0.0000000064442714D+00, &
          0.9999999935557286D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9999999935557286D+00, &
          0.0000000064442714D+00, &
          0.0000000064442714D+00, &
          0.9999999935557286D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.8467830886160004D+00, &
          0.1532169113839996D+00, &
          0.8467830886160004D+00, &
          0.1532169113839996D+00, &
          0.8467830886160004D+00, &
          0.1532169113839996D+00, &
          0.8467830886160004D+00, &
          0.1532169113839996D+00, &
          0.6731903387258514D+00, &
          0.3268096612741486D+00, &
          0.3268096612741486D+00, &
          0.6731903387258514D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.6731903387258514D+00, &
          0.3268096612741486D+00, &
          0.3268096612741486D+00, &
          0.6731903387258514D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9359181010210101D+00, &
          0.0640818989789899D+00, &
          0.9359181010210101D+00, &
          0.0640818989789899D+00, &
          0.9359181010210101D+00, &
          0.0640818989789899D+00, &
          0.9359181010210101D+00, &
          0.0640818989789899D+00, &
          0.1873892499846415D+00, &
          0.8126107500153585D+00, &
          0.8126107500153585D+00, &
          0.1873892499846415D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.1873892499846415D+00, &
          0.8126107500153585D+00, &
          0.8126107500153585D+00, &
          0.1873892499846415D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.2917727295202502D+00, &
          0.9026740712514576D+00, &
          0.2917727295202502D+00, &
          0.7082272704797498D+00, &
          0.9026740712514576D+00, &
          0.7082272704797498D+00, &
          0.7082272704797498D+00, &
          0.9026740712514576D+00, &
          0.7082272704797498D+00, &
          0.2917727295202502D+00, &
          0.9026740712514576D+00, &
          0.2917727295202502D+00, &
          0.2917727295202502D+00, &
          0.0973259287485423D+00, &
          0.2917727295202502D+00, &
          0.7082272704797498D+00, &
          0.0973259287485423D+00, &
          0.7082272704797498D+00, &
          0.7082272704797498D+00, &
          0.0973259287485423D+00, &
          0.7082272704797498D+00, &
          0.2917727295202502D+00, &
          0.0973259287485423D+00, &
          0.2917727295202502D+00, &
          0.9869259796545771D+00, &
          0.9372501994187228D+00, &
          0.9869259796545771D+00, &
          0.0130740203454229D+00, &
          0.9372501994187228D+00, &
          0.0130740203454229D+00, &
          0.0130740203454229D+00, &
          0.9372501994187228D+00, &
          0.0130740203454229D+00, &
          0.9869259796545771D+00, &
          0.9372501994187228D+00, &
          0.9869259796545771D+00, &
          0.9869259796545771D+00, &
          0.0627498005812772D+00, &
          0.9869259796545771D+00, &
          0.0130740203454229D+00, &
          0.0627498005812772D+00, &
          0.0130740203454229D+00, &
          0.0130740203454229D+00, &
          0.0627498005812772D+00, &
          0.0130740203454229D+00, &
          0.9869259796545771D+00, &
          0.0627498005812772D+00, &
          0.9869259796545771D+00, &
          0.8645117252804405D+00, &
          0.6225837716647646D+00, &
          0.8645117252804405D+00, &
          0.1354882747195595D+00, &
          0.6225837716647646D+00, &
          0.1354882747195595D+00, &
          0.1354882747195595D+00, &
          0.6225837716647646D+00, &
          0.1354882747195595D+00, &
          0.8645117252804405D+00, &
          0.6225837716647646D+00, &
          0.8645117252804405D+00, &
          0.8645117252804405D+00, &
          0.3774162283352354D+00, &
          0.8645117252804405D+00, &
          0.1354882747195595D+00, &
          0.3774162283352354D+00, &
          0.1354882747195595D+00, &
          0.1354882747195595D+00, &
          0.3774162283352354D+00, &
          0.1354882747195595D+00, &
          0.8645117252804405D+00, &
          0.3774162283352354D+00, &
          0.8645117252804405D+00, &
          0.8248058551406905D+00, &
          0.9220515632141263D+00, &
          0.8248058551406905D+00, &
          0.1751941448593096D+00, &
          0.9220515632141263D+00, &
          0.1751941448593096D+00, &
          0.1751941448593096D+00, &
          0.9220515632141263D+00, &
          0.1751941448593096D+00, &
          0.8248058551406905D+00, &
          0.9220515632141263D+00, &
          0.8248058551406905D+00, &
          0.8248058551406905D+00, &
          0.0779484367858737D+00, &
          0.8248058551406905D+00, &
          0.1751941448593096D+00, &
          0.0779484367858737D+00, &
          0.1751941448593096D+00, &
          0.1751941448593096D+00, &
          0.0779484367858737D+00, &
          0.1751941448593096D+00, &
          0.8248058551406905D+00, &
          0.0779484367858737D+00, &
          0.8248058551406905D+00, &
          0.9459853758476011D+00, &
          0.5650967942783419D+00, &
          0.9459853758476011D+00, &
          0.0540146241523989D+00, &
          0.5650967942783419D+00, &
          0.0540146241523989D+00, &
          0.0540146241523989D+00, &
          0.5650967942783419D+00, &
          0.0540146241523989D+00, &
          0.9459853758476011D+00, &
          0.5650967942783419D+00, &
          0.9459853758476011D+00, &
          0.9459853758476011D+00, &
          0.4349032057216580D+00, &
          0.9459853758476011D+00, &
          0.0540146241523989D+00, &
          0.4349032057216580D+00, &
          0.0540146241523989D+00, &
          0.0540146241523989D+00, &
          0.4349032057216580D+00, &
          0.0540146241523989D+00, &
          0.9459853758476011D+00, &
          0.4349032057216580D+00, &
          0.9459853758476011D+00, &
          0.8916635718645309D+00, &
          0.9766093536267083D+00, &
          0.8916635718645309D+00, &
          0.1083364281354691D+00, &
          0.9766093536267083D+00, &
          0.1083364281354691D+00, &
          0.1083364281354691D+00, &
          0.9766093536267083D+00, &
          0.1083364281354691D+00, &
          0.8916635718645309D+00, &
          0.9766093536267083D+00, &
          0.8916635718645309D+00, &
          0.8916635718645309D+00, &
          0.0233906463732917D+00, &
          0.8916635718645309D+00, &
          0.1083364281354691D+00, &
          0.0233906463732917D+00, &
          0.1083364281354691D+00, &
          0.1083364281354691D+00, &
          0.0233906463732917D+00, &
          0.1083364281354691D+00, &
          0.8916635718645309D+00, &
          0.0233906463732917D+00, &
          0.8916635718645309D+00, &
          0.9849705936644424D+00, &
          0.6355079664366210D+00, &
          0.9849705936644424D+00, &
          0.0150294063355575D+00, &
          0.6355079664366210D+00, &
          0.0150294063355575D+00, &
          0.0150294063355575D+00, &
          0.6355079664366210D+00, &
          0.0150294063355575D+00, &
          0.9849705936644424D+00, &
          0.6355079664366210D+00, &
          0.9849705936644424D+00, &
          0.9849705936644424D+00, &
          0.3644920335633790D+00, &
          0.9849705936644424D+00, &
          0.0150294063355575D+00, &
          0.3644920335633790D+00, &
          0.0150294063355575D+00, &
          0.0150294063355575D+00, &
          0.3644920335633790D+00, &
          0.0150294063355575D+00, &
          0.9849705936644424D+00, &
          0.3644920335633790D+00, &
          0.9849705936644424D+00, &
          0.5860169886103869D+00, &
          0.9923845018898745D+00, &
          0.5860169886103869D+00, &
          0.4139830113896131D+00, &
          0.9923845018898745D+00, &
          0.4139830113896131D+00, &
          0.4139830113896131D+00, &
          0.9923845018898745D+00, &
          0.4139830113896131D+00, &
          0.5860169886103869D+00, &
          0.9923845018898745D+00, &
          0.5860169886103869D+00, &
          0.5860169886103869D+00, &
          0.0076154981101255D+00, &
          0.5860169886103869D+00, &
          0.4139830113896131D+00, &
          0.0076154981101255D+00, &
          0.4139830113896131D+00, &
          0.4139830113896131D+00, &
          0.0076154981101255D+00, &
          0.4139830113896131D+00, &
          0.5860169886103869D+00, &
          0.0076154981101255D+00, &
          0.5860169886103869D+00, &
          0.7690492963095925D+00, &
          0.3896184326892077D+00, &
          0.7690492963095925D+00, &
          0.2309507036904074D+00, &
          0.3896184326892077D+00, &
          0.2309507036904074D+00, &
          0.2309507036904074D+00, &
          0.3896184326892077D+00, &
          0.2309507036904074D+00, &
          0.7690492963095925D+00, &
          0.3896184326892077D+00, &
          0.7690492963095925D+00, &
          0.7690492963095925D+00, &
          0.6103815673107923D+00, &
          0.7690492963095925D+00, &
          0.2309507036904074D+00, &
          0.6103815673107923D+00, &
          0.2309507036904074D+00, &
          0.2309507036904074D+00, &
          0.6103815673107923D+00, &
          0.2309507036904074D+00, &
          0.7690492963095925D+00, &
          0.6103815673107923D+00, &
          0.7690492963095925D+00, &
          0.7768424517946938D+00, &
          0.9994864489115713D+00, &
          0.7768424517946938D+00, &
          0.2231575482053062D+00, &
          0.9994864489115713D+00, &
          0.2231575482053062D+00, &
          0.2231575482053062D+00, &
          0.9994864489115713D+00, &
          0.2231575482053062D+00, &
          0.7768424517946938D+00, &
          0.9994864489115713D+00, &
          0.7768424517946938D+00, &
          0.7768424517946938D+00, &
          0.0005135510884287D+00, &
          0.7768424517946938D+00, &
          0.2231575482053062D+00, &
          0.0005135510884287D+00, &
          0.2231575482053062D+00, &
          0.2231575482053062D+00, &
          0.0005135510884287D+00, &
          0.2231575482053062D+00, &
          0.7768424517946938D+00, &
          0.0005135510884287D+00, &
          0.7768424517946938D+00, &
          0.6117790069755276D+00, &
          0.7209288728023533D+00, &
          0.6117790069755276D+00, &
          0.3882209930244724D+00, &
          0.7209288728023533D+00, &
          0.3882209930244724D+00, &
          0.3882209930244724D+00, &
          0.7209288728023533D+00, &
          0.3882209930244724D+00, &
          0.6117790069755276D+00, &
          0.7209288728023533D+00, &
          0.6117790069755276D+00, &
          0.6117790069755276D+00, &
          0.2790711271976467D+00, &
          0.6117790069755276D+00, &
          0.3882209930244724D+00, &
          0.2790711271976467D+00, &
          0.3882209930244724D+00, &
          0.3882209930244724D+00, &
          0.2790711271976467D+00, &
          0.3882209930244724D+00, &
          0.6117790069755276D+00, &
          0.2790711271976467D+00, &
          0.6117790069755276D+00, &
          0.9350312907963207D+00, &
          0.7494073604468848D+00, &
          0.9350312907963207D+00, &
          0.0649687092036792D+00, &
          0.7494073604468848D+00, &
          0.0649687092036792D+00, &
          0.0649687092036792D+00, &
          0.7494073604468848D+00, &
          0.0649687092036792D+00, &
          0.9350312907963207D+00, &
          0.7494073604468848D+00, &
          0.9350312907963207D+00, &
          0.9350312907963207D+00, &
          0.2505926395531152D+00, &
          0.9350312907963207D+00, &
          0.0649687092036792D+00, &
          0.2505926395531152D+00, &
          0.0649687092036792D+00, &
          0.0649687092036792D+00, &
          0.2505926395531152D+00, &
          0.0649687092036792D+00, &
          0.9350312907963207D+00, &
          0.2505926395531152D+00, &
          0.9350312907963207D+00, &
          0.8750719993006615D+00, &
          0.8750719993006615D+00, &
          0.9719926774632096D+00, &
          0.6986057011288138D+00, &
          0.9719926774632096D+00, &
          0.6986057011288138D+00, &
          0.1249280006993386D+00, &
          0.1249280006993386D+00, &
          0.9719926774632096D+00, &
          0.6986057011288138D+00, &
          0.9719926774632096D+00, &
          0.6986057011288138D+00, &
          0.8750719993006615D+00, &
          0.8750719993006615D+00, &
          0.0280073225367904D+00, &
          0.6986057011288138D+00, &
          0.0280073225367904D+00, &
          0.6986057011288138D+00, &
          0.8750719993006615D+00, &
          0.8750719993006615D+00, &
          0.9719926774632096D+00, &
          0.3013942988711862D+00, &
          0.9719926774632096D+00, &
          0.3013942988711862D+00, &
          0.1249280006993386D+00, &
          0.1249280006993386D+00, &
          0.0280073225367904D+00, &
          0.6986057011288138D+00, &
          0.0280073225367904D+00, &
          0.6986057011288138D+00, &
          0.1249280006993386D+00, &
          0.1249280006993386D+00, &
          0.9719926774632096D+00, &
          0.3013942988711862D+00, &
          0.9719926774632096D+00, &
          0.3013942988711862D+00, &
          0.8750719993006615D+00, &
          0.8750719993006615D+00, &
          0.0280073225367904D+00, &
          0.3013942988711862D+00, &
          0.0280073225367904D+00, &
          0.3013942988711862D+00, &
          0.1249280006993386D+00, &
          0.1249280006993386D+00, &
          0.0280073225367904D+00, &
          0.3013942988711862D+00, &
          0.0280073225367904D+00, &
          0.3013942988711862D+00, &
          0.9512929061479376D+00, &
          0.9512929061479376D+00, &
          0.9923118941018740D+00, &
          0.8153563222457214D+00, &
          0.9923118941018740D+00, &
          0.8153563222457214D+00, &
          0.0487070938520624D+00, &
          0.0487070938520624D+00, &
          0.9923118941018740D+00, &
          0.8153563222457214D+00, &
          0.9923118941018740D+00, &
          0.8153563222457214D+00, &
          0.9512929061479376D+00, &
          0.9512929061479376D+00, &
          0.0076881058981259D+00, &
          0.8153563222457214D+00, &
          0.0076881058981259D+00, &
          0.8153563222457214D+00, &
          0.9512929061479376D+00, &
          0.9512929061479376D+00, &
          0.9923118941018740D+00, &
          0.1846436777542786D+00, &
          0.9923118941018740D+00, &
          0.1846436777542786D+00, &
          0.0487070938520624D+00, &
          0.0487070938520624D+00, &
          0.0076881058981259D+00, &
          0.8153563222457214D+00, &
          0.0076881058981259D+00, &
          0.8153563222457214D+00, &
          0.0487070938520624D+00, &
          0.0487070938520624D+00, &
          0.9923118941018740D+00, &
          0.1846436777542786D+00, &
          0.9923118941018740D+00, &
          0.1846436777542786D+00, &
          0.9512929061479376D+00, &
          0.9512929061479376D+00, &
          0.0076881058981259D+00, &
          0.1846436777542786D+00, &
          0.0076881058981259D+00, &
          0.1846436777542786D+00, &
          0.0487070938520624D+00, &
          0.0487070938520624D+00, &
          0.0076881058981259D+00, &
          0.1846436777542786D+00, &
          0.0076881058981259D+00, &
          0.1846436777542786D+00, &
          0.6108897298239304D+00, &
          0.6108897298239304D+00, &
          0.5491142465864416D+00, &
          0.9360994289291793D+00, &
          0.5491142465864416D+00, &
          0.9360994289291793D+00, &
          0.3891102701760696D+00, &
          0.3891102701760696D+00, &
          0.5491142465864416D+00, &
          0.9360994289291793D+00, &
          0.5491142465864416D+00, &
          0.9360994289291793D+00, &
          0.6108897298239304D+00, &
          0.6108897298239304D+00, &
          0.4508857534135584D+00, &
          0.9360994289291793D+00, &
          0.4508857534135584D+00, &
          0.9360994289291793D+00, &
          0.6108897298239304D+00, &
          0.6108897298239304D+00, &
          0.5491142465864416D+00, &
          0.0639005710708208D+00, &
          0.5491142465864416D+00, &
          0.0639005710708208D+00, &
          0.3891102701760696D+00, &
          0.3891102701760696D+00, &
          0.4508857534135584D+00, &
          0.9360994289291793D+00, &
          0.4508857534135584D+00, &
          0.9360994289291793D+00, &
          0.3891102701760696D+00, &
          0.3891102701760696D+00, &
          0.5491142465864416D+00, &
          0.0639005710708208D+00, &
          0.5491142465864416D+00, &
          0.0639005710708208D+00, &
          0.6108897298239304D+00, &
          0.6108897298239304D+00, &
          0.4508857534135584D+00, &
          0.0639005710708208D+00, &
          0.4508857534135584D+00, &
          0.0639005710708208D+00, &
          0.3891102701760696D+00, &
          0.3891102701760696D+00, &
          0.4508857534135584D+00, &
          0.0639005710708208D+00, &
          0.4508857534135584D+00, &
          0.0639005710708208D+00, &
          0.9784211281771267D+00, &
          0.9784211281771267D+00, &
          0.6040696445748173D+00, &
          0.7479204633014319D+00, &
          0.6040696445748173D+00, &
          0.7479204633014319D+00, &
          0.0215788718228733D+00, &
          0.0215788718228733D+00, &
          0.6040696445748173D+00, &
          0.7479204633014319D+00, &
          0.6040696445748173D+00, &
          0.7479204633014319D+00, &
          0.9784211281771267D+00, &
          0.9784211281771267D+00, &
          0.3959303554251827D+00, &
          0.7479204633014319D+00, &
          0.3959303554251827D+00, &
          0.7479204633014319D+00, &
          0.9784211281771267D+00, &
          0.9784211281771267D+00, &
          0.6040696445748173D+00, &
          0.2520795366985681D+00, &
          0.6040696445748173D+00, &
          0.2520795366985681D+00, &
          0.0215788718228733D+00, &
          0.0215788718228733D+00, &
          0.3959303554251827D+00, &
          0.7479204633014319D+00, &
          0.3959303554251827D+00, &
          0.7479204633014319D+00, &
          0.0215788718228733D+00, &
          0.0215788718228733D+00, &
          0.6040696445748173D+00, &
          0.2520795366985681D+00, &
          0.6040696445748173D+00, &
          0.2520795366985681D+00, &
          0.9784211281771267D+00, &
          0.9784211281771267D+00, &
          0.3959303554251827D+00, &
          0.2520795366985681D+00, &
          0.3959303554251827D+00, &
          0.2520795366985681D+00, &
          0.0215788718228733D+00, &
          0.0215788718228733D+00, &
          0.3959303554251827D+00, &
          0.2520795366985681D+00, &
          0.3959303554251827D+00, &
          0.2520795366985681D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.5000000000000000D+00, &
          0.1912521437383576D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.8087478562616424D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.3543162299281049D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.6456837700718951D+00, &
          0.5000000000000000D+00, &
          0.9311291591147122D+00, &
          0.0688708408852879D+00, &
          0.9311291591147122D+00, &
          0.0688708408852879D+00, &
          0.9311291591147122D+00, &
          0.0688708408852879D+00, &
          0.0688708408852879D+00, &
          0.9311291591147122D+00, &
          0.7833399058047521D+00, &
          0.2166600941952480D+00, &
          0.7833399058047521D+00, &
          0.2166600941952480D+00, &
          0.7833399058047521D+00, &
          0.2166600941952480D+00, &
          0.2166600941952480D+00, &
          0.7833399058047521D+00, &
          0.9999999935557286D+00, &
          0.0000000064442714D+00, &
          0.0000000064442714D+00, &
          0.9999999935557286D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.8981144288438190D+00, &
          0.1018855711561810D+00, &
          0.8981144288438190D+00, &
          0.1018855711561810D+00, &
          0.8981144288438190D+00, &
          0.1018855711561810D+00, &
          0.8981144288438190D+00, &
          0.1018855711561810D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9999999935557286D+00, &
          0.0000000064442714D+00, &
          0.0000000064442714D+00, &
          0.9999999935557286D+00, &
          0.6731903387258514D+00, &
          0.3268096612741486D+00, &
          0.3268096612741486D+00, &
          0.6731903387258514D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.8467830886160004D+00, &
          0.1532169113839996D+00, &
          0.8467830886160004D+00, &
          0.1532169113839996D+00, &
          0.8467830886160004D+00, &
          0.1532169113839996D+00, &
          0.8467830886160004D+00, &
          0.1532169113839996D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.6731903387258514D+00, &
          0.3268096612741486D+00, &
          0.3268096612741486D+00, &
          0.6731903387258514D+00, &
          0.1873892499846415D+00, &
          0.8126107500153585D+00, &
          0.8126107500153585D+00, &
          0.1873892499846415D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9359181010210101D+00, &
          0.0640818989789899D+00, &
          0.9359181010210101D+00, &
          0.0640818989789899D+00, &
          0.9359181010210101D+00, &
          0.0640818989789899D+00, &
          0.9359181010210101D+00, &
          0.0640818989789899D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.1873892499846415D+00, &
          0.8126107500153585D+00, &
          0.8126107500153585D+00, &
          0.1873892499846415D+00, &
          0.9026740712514576D+00, &
          0.2917727295202502D+00, &
          0.2917727295202502D+00, &
          0.9026740712514576D+00, &
          0.7082272704797498D+00, &
          0.2917727295202502D+00, &
          0.9026740712514576D+00, &
          0.7082272704797498D+00, &
          0.7082272704797498D+00, &
          0.9026740712514576D+00, &
          0.2917727295202502D+00, &
          0.7082272704797498D+00, &
          0.0973259287485423D+00, &
          0.2917727295202502D+00, &
          0.2917727295202502D+00, &
          0.0973259287485423D+00, &
          0.7082272704797498D+00, &
          0.2917727295202502D+00, &
          0.0973259287485423D+00, &
          0.7082272704797498D+00, &
          0.7082272704797498D+00, &
          0.0973259287485423D+00, &
          0.2917727295202502D+00, &
          0.7082272704797498D+00, &
          0.9372501994187228D+00, &
          0.9869259796545771D+00, &
          0.9869259796545771D+00, &
          0.9372501994187228D+00, &
          0.0130740203454229D+00, &
          0.9869259796545771D+00, &
          0.9372501994187228D+00, &
          0.0130740203454229D+00, &
          0.0130740203454229D+00, &
          0.9372501994187228D+00, &
          0.9869259796545771D+00, &
          0.0130740203454229D+00, &
          0.0627498005812772D+00, &
          0.9869259796545771D+00, &
          0.9869259796545771D+00, &
          0.0627498005812772D+00, &
          0.0130740203454229D+00, &
          0.9869259796545771D+00, &
          0.0627498005812772D+00, &
          0.0130740203454229D+00, &
          0.0130740203454229D+00, &
          0.0627498005812772D+00, &
          0.9869259796545771D+00, &
          0.0130740203454229D+00, &
          0.6225837716647646D+00, &
          0.8645117252804405D+00, &
          0.8645117252804405D+00, &
          0.6225837716647646D+00, &
          0.1354882747195595D+00, &
          0.8645117252804405D+00, &
          0.6225837716647646D+00, &
          0.1354882747195595D+00, &
          0.1354882747195595D+00, &
          0.6225837716647646D+00, &
          0.8645117252804405D+00, &
          0.1354882747195595D+00, &
          0.3774162283352354D+00, &
          0.8645117252804405D+00, &
          0.8645117252804405D+00, &
          0.3774162283352354D+00, &
          0.1354882747195595D+00, &
          0.8645117252804405D+00, &
          0.3774162283352354D+00, &
          0.1354882747195595D+00, &
          0.1354882747195595D+00, &
          0.3774162283352354D+00, &
          0.8645117252804405D+00, &
          0.1354882747195595D+00, &
          0.9220515632141263D+00, &
          0.8248058551406905D+00, &
          0.8248058551406905D+00, &
          0.9220515632141263D+00, &
          0.1751941448593096D+00, &
          0.8248058551406905D+00, &
          0.9220515632141263D+00, &
          0.1751941448593096D+00, &
          0.1751941448593096D+00, &
          0.9220515632141263D+00, &
          0.8248058551406905D+00, &
          0.1751941448593096D+00, &
          0.0779484367858737D+00, &
          0.8248058551406905D+00, &
          0.8248058551406905D+00, &
          0.0779484367858737D+00, &
          0.1751941448593096D+00, &
          0.8248058551406905D+00, &
          0.0779484367858737D+00, &
          0.1751941448593096D+00, &
          0.1751941448593096D+00, &
          0.0779484367858737D+00, &
          0.8248058551406905D+00, &
          0.1751941448593096D+00, &
          0.5650967942783419D+00, &
          0.9459853758476011D+00, &
          0.9459853758476011D+00, &
          0.5650967942783419D+00, &
          0.0540146241523989D+00, &
          0.9459853758476011D+00, &
          0.5650967942783419D+00, &
          0.0540146241523989D+00, &
          0.0540146241523989D+00, &
          0.5650967942783419D+00, &
          0.9459853758476011D+00, &
          0.0540146241523989D+00, &
          0.4349032057216580D+00, &
          0.9459853758476011D+00, &
          0.9459853758476011D+00, &
          0.4349032057216580D+00, &
          0.0540146241523989D+00, &
          0.9459853758476011D+00, &
          0.4349032057216580D+00, &
          0.0540146241523989D+00, &
          0.0540146241523989D+00, &
          0.4349032057216580D+00, &
          0.9459853758476011D+00, &
          0.0540146241523989D+00, &
          0.9766093536267083D+00, &
          0.8916635718645309D+00, &
          0.8916635718645309D+00, &
          0.9766093536267083D+00, &
          0.1083364281354691D+00, &
          0.8916635718645309D+00, &
          0.9766093536267083D+00, &
          0.1083364281354691D+00, &
          0.1083364281354691D+00, &
          0.9766093536267083D+00, &
          0.8916635718645309D+00, &
          0.1083364281354691D+00, &
          0.0233906463732917D+00, &
          0.8916635718645309D+00, &
          0.8916635718645309D+00, &
          0.0233906463732917D+00, &
          0.1083364281354691D+00, &
          0.8916635718645309D+00, &
          0.0233906463732917D+00, &
          0.1083364281354691D+00, &
          0.1083364281354691D+00, &
          0.0233906463732917D+00, &
          0.8916635718645309D+00, &
          0.1083364281354691D+00, &
          0.6355079664366210D+00, &
          0.9849705936644424D+00, &
          0.9849705936644424D+00, &
          0.6355079664366210D+00, &
          0.0150294063355575D+00, &
          0.9849705936644424D+00, &
          0.6355079664366210D+00, &
          0.0150294063355575D+00, &
          0.0150294063355575D+00, &
          0.6355079664366210D+00, &
          0.9849705936644424D+00, &
          0.0150294063355575D+00, &
          0.3644920335633790D+00, &
          0.9849705936644424D+00, &
          0.9849705936644424D+00, &
          0.3644920335633790D+00, &
          0.0150294063355575D+00, &
          0.9849705936644424D+00, &
          0.3644920335633790D+00, &
          0.0150294063355575D+00, &
          0.0150294063355575D+00, &
          0.3644920335633790D+00, &
          0.9849705936644424D+00, &
          0.0150294063355575D+00, &
          0.9923845018898745D+00, &
          0.5860169886103869D+00, &
          0.5860169886103869D+00, &
          0.9923845018898745D+00, &
          0.4139830113896131D+00, &
          0.5860169886103869D+00, &
          0.9923845018898745D+00, &
          0.4139830113896131D+00, &
          0.4139830113896131D+00, &
          0.9923845018898745D+00, &
          0.5860169886103869D+00, &
          0.4139830113896131D+00, &
          0.0076154981101255D+00, &
          0.5860169886103869D+00, &
          0.5860169886103869D+00, &
          0.0076154981101255D+00, &
          0.4139830113896131D+00, &
          0.5860169886103869D+00, &
          0.0076154981101255D+00, &
          0.4139830113896131D+00, &
          0.4139830113896131D+00, &
          0.0076154981101255D+00, &
          0.5860169886103869D+00, &
          0.4139830113896131D+00, &
          0.3896184326892077D+00, &
          0.7690492963095925D+00, &
          0.7690492963095925D+00, &
          0.3896184326892077D+00, &
          0.2309507036904074D+00, &
          0.7690492963095925D+00, &
          0.3896184326892077D+00, &
          0.2309507036904074D+00, &
          0.2309507036904074D+00, &
          0.3896184326892077D+00, &
          0.7690492963095925D+00, &
          0.2309507036904074D+00, &
          0.6103815673107923D+00, &
          0.7690492963095925D+00, &
          0.7690492963095925D+00, &
          0.6103815673107923D+00, &
          0.2309507036904074D+00, &
          0.7690492963095925D+00, &
          0.6103815673107923D+00, &
          0.2309507036904074D+00, &
          0.2309507036904074D+00, &
          0.6103815673107923D+00, &
          0.7690492963095925D+00, &
          0.2309507036904074D+00, &
          0.9994864489115713D+00, &
          0.7768424517946938D+00, &
          0.7768424517946938D+00, &
          0.9994864489115713D+00, &
          0.2231575482053062D+00, &
          0.7768424517946938D+00, &
          0.9994864489115713D+00, &
          0.2231575482053062D+00, &
          0.2231575482053062D+00, &
          0.9994864489115713D+00, &
          0.7768424517946938D+00, &
          0.2231575482053062D+00, &
          0.0005135510884287D+00, &
          0.7768424517946938D+00, &
          0.7768424517946938D+00, &
          0.0005135510884287D+00, &
          0.2231575482053062D+00, &
          0.7768424517946938D+00, &
          0.0005135510884287D+00, &
          0.2231575482053062D+00, &
          0.2231575482053062D+00, &
          0.0005135510884287D+00, &
          0.7768424517946938D+00, &
          0.2231575482053062D+00, &
          0.7209288728023533D+00, &
          0.6117790069755276D+00, &
          0.6117790069755276D+00, &
          0.7209288728023533D+00, &
          0.3882209930244724D+00, &
          0.6117790069755276D+00, &
          0.7209288728023533D+00, &
          0.3882209930244724D+00, &
          0.3882209930244724D+00, &
          0.7209288728023533D+00, &
          0.6117790069755276D+00, &
          0.3882209930244724D+00, &
          0.2790711271976467D+00, &
          0.6117790069755276D+00, &
          0.6117790069755276D+00, &
          0.2790711271976467D+00, &
          0.3882209930244724D+00, &
          0.6117790069755276D+00, &
          0.2790711271976467D+00, &
          0.3882209930244724D+00, &
          0.3882209930244724D+00, &
          0.2790711271976467D+00, &
          0.6117790069755276D+00, &
          0.3882209930244724D+00, &
          0.7494073604468848D+00, &
          0.9350312907963207D+00, &
          0.9350312907963207D+00, &
          0.7494073604468848D+00, &
          0.0649687092036792D+00, &
          0.9350312907963207D+00, &
          0.7494073604468848D+00, &
          0.0649687092036792D+00, &
          0.0649687092036792D+00, &
          0.7494073604468848D+00, &
          0.9350312907963207D+00, &
          0.0649687092036792D+00, &
          0.2505926395531152D+00, &
          0.9350312907963207D+00, &
          0.9350312907963207D+00, &
          0.2505926395531152D+00, &
          0.0649687092036792D+00, &
          0.9350312907963207D+00, &
          0.2505926395531152D+00, &
          0.0649687092036792D+00, &
          0.0649687092036792D+00, &
          0.2505926395531152D+00, &
          0.9350312907963207D+00, &
          0.0649687092036792D+00, &
          0.9719926774632096D+00, &
          0.6986057011288138D+00, &
          0.8750719993006615D+00, &
          0.8750719993006615D+00, &
          0.6986057011288138D+00, &
          0.9719926774632096D+00, &
          0.9719926774632096D+00, &
          0.6986057011288138D+00, &
          0.1249280006993386D+00, &
          0.1249280006993386D+00, &
          0.6986057011288138D+00, &
          0.9719926774632096D+00, &
          0.0280073225367904D+00, &
          0.6986057011288138D+00, &
          0.8750719993006615D+00, &
          0.8750719993006615D+00, &
          0.6986057011288138D+00, &
          0.0280073225367904D+00, &
          0.9719926774632096D+00, &
          0.3013942988711862D+00, &
          0.8750719993006615D+00, &
          0.8750719993006615D+00, &
          0.3013942988711862D+00, &
          0.9719926774632096D+00, &
          0.0280073225367904D+00, &
          0.6986057011288138D+00, &
          0.1249280006993386D+00, &
          0.1249280006993386D+00, &
          0.6986057011288138D+00, &
          0.0280073225367904D+00, &
          0.9719926774632096D+00, &
          0.3013942988711862D+00, &
          0.1249280006993386D+00, &
          0.1249280006993386D+00, &
          0.3013942988711862D+00, &
          0.9719926774632096D+00, &
          0.0280073225367904D+00, &
          0.3013942988711862D+00, &
          0.8750719993006615D+00, &
          0.8750719993006615D+00, &
          0.3013942988711862D+00, &
          0.0280073225367904D+00, &
          0.0280073225367904D+00, &
          0.3013942988711862D+00, &
          0.1249280006993386D+00, &
          0.1249280006993386D+00, &
          0.3013942988711862D+00, &
          0.0280073225367904D+00, &
          0.9923118941018740D+00, &
          0.8153563222457214D+00, &
          0.9512929061479376D+00, &
          0.9512929061479376D+00, &
          0.8153563222457214D+00, &
          0.9923118941018740D+00, &
          0.9923118941018740D+00, &
          0.8153563222457214D+00, &
          0.0487070938520624D+00, &
          0.0487070938520624D+00, &
          0.8153563222457214D+00, &
          0.9923118941018740D+00, &
          0.0076881058981259D+00, &
          0.8153563222457214D+00, &
          0.9512929061479376D+00, &
          0.9512929061479376D+00, &
          0.8153563222457214D+00, &
          0.0076881058981259D+00, &
          0.9923118941018740D+00, &
          0.1846436777542786D+00, &
          0.9512929061479376D+00, &
          0.9512929061479376D+00, &
          0.1846436777542786D+00, &
          0.9923118941018740D+00, &
          0.0076881058981259D+00, &
          0.8153563222457214D+00, &
          0.0487070938520624D+00, &
          0.0487070938520624D+00, &
          0.8153563222457214D+00, &
          0.0076881058981259D+00, &
          0.9923118941018740D+00, &
          0.1846436777542786D+00, &
          0.0487070938520624D+00, &
          0.0487070938520624D+00, &
          0.1846436777542786D+00, &
          0.9923118941018740D+00, &
          0.0076881058981259D+00, &
          0.1846436777542786D+00, &
          0.9512929061479376D+00, &
          0.9512929061479376D+00, &
          0.1846436777542786D+00, &
          0.0076881058981259D+00, &
          0.0076881058981259D+00, &
          0.1846436777542786D+00, &
          0.0487070938520624D+00, &
          0.0487070938520624D+00, &
          0.1846436777542786D+00, &
          0.0076881058981259D+00, &
          0.5491142465864416D+00, &
          0.9360994289291793D+00, &
          0.6108897298239304D+00, &
          0.6108897298239304D+00, &
          0.9360994289291793D+00, &
          0.5491142465864416D+00, &
          0.5491142465864416D+00, &
          0.9360994289291793D+00, &
          0.3891102701760696D+00, &
          0.3891102701760696D+00, &
          0.9360994289291793D+00, &
          0.5491142465864416D+00, &
          0.4508857534135584D+00, &
          0.9360994289291793D+00, &
          0.6108897298239304D+00, &
          0.6108897298239304D+00, &
          0.9360994289291793D+00, &
          0.4508857534135584D+00, &
          0.5491142465864416D+00, &
          0.0639005710708208D+00, &
          0.6108897298239304D+00, &
          0.6108897298239304D+00, &
          0.0639005710708208D+00, &
          0.5491142465864416D+00, &
          0.4508857534135584D+00, &
          0.9360994289291793D+00, &
          0.3891102701760696D+00, &
          0.3891102701760696D+00, &
          0.9360994289291793D+00, &
          0.4508857534135584D+00, &
          0.5491142465864416D+00, &
          0.0639005710708208D+00, &
          0.3891102701760696D+00, &
          0.3891102701760696D+00, &
          0.0639005710708208D+00, &
          0.5491142465864416D+00, &
          0.4508857534135584D+00, &
          0.0639005710708208D+00, &
          0.6108897298239304D+00, &
          0.6108897298239304D+00, &
          0.0639005710708208D+00, &
          0.4508857534135584D+00, &
          0.4508857534135584D+00, &
          0.0639005710708208D+00, &
          0.3891102701760696D+00, &
          0.3891102701760696D+00, &
          0.0639005710708208D+00, &
          0.4508857534135584D+00, &
          0.6040696445748173D+00, &
          0.7479204633014319D+00, &
          0.9784211281771267D+00, &
          0.9784211281771267D+00, &
          0.7479204633014319D+00, &
          0.6040696445748173D+00, &
          0.6040696445748173D+00, &
          0.7479204633014319D+00, &
          0.0215788718228733D+00, &
          0.0215788718228733D+00, &
          0.7479204633014319D+00, &
          0.6040696445748173D+00, &
          0.3959303554251827D+00, &
          0.7479204633014319D+00, &
          0.9784211281771267D+00, &
          0.9784211281771267D+00, &
          0.7479204633014319D+00, &
          0.3959303554251827D+00, &
          0.6040696445748173D+00, &
          0.2520795366985681D+00, &
          0.9784211281771267D+00, &
          0.9784211281771267D+00, &
          0.2520795366985681D+00, &
          0.6040696445748173D+00, &
          0.3959303554251827D+00, &
          0.7479204633014319D+00, &
          0.0215788718228733D+00, &
          0.0215788718228733D+00, &
          0.7479204633014319D+00, &
          0.3959303554251827D+00, &
          0.6040696445748173D+00, &
          0.2520795366985681D+00, &
          0.0215788718228733D+00, &
          0.0215788718228733D+00, &
          0.2520795366985681D+00, &
          0.6040696445748173D+00, &
          0.3959303554251827D+00, &
          0.2520795366985681D+00, &
          0.9784211281771267D+00, &
          0.9784211281771267D+00, &
          0.2520795366985681D+00, &
          0.3959303554251827D+00, &
          0.3959303554251827D+00, &
          0.2520795366985681D+00, &
          0.0215788718228733D+00, &
          0.0215788718228733D+00, &
          0.2520795366985681D+00, &
          0.3959303554251827D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.1912521437383576D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.8087478562616424D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.3543162299281049D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.6456837700718951D+00, &
          0.9311291591147122D+00, &
          0.0688708408852879D+00, &
          0.9311291591147122D+00, &
          0.9311291591147122D+00, &
          0.0688708408852879D+00, &
          0.9311291591147122D+00, &
          0.0688708408852879D+00, &
          0.0688708408852879D+00, &
          0.7833399058047521D+00, &
          0.2166600941952480D+00, &
          0.7833399058047521D+00, &
          0.7833399058047521D+00, &
          0.2166600941952480D+00, &
          0.7833399058047521D+00, &
          0.2166600941952480D+00, &
          0.2166600941952480D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9999999935557286D+00, &
          0.0000000064442714D+00, &
          0.0000000064442714D+00, &
          0.9999999935557286D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.9999999935557286D+00, &
          0.0000000064442714D+00, &
          0.0000000064442714D+00, &
          0.9999999935557286D+00, &
          0.8981144288438190D+00, &
          0.1018855711561810D+00, &
          0.8981144288438190D+00, &
          0.1018855711561810D+00, &
          0.8981144288438190D+00, &
          0.1018855711561810D+00, &
          0.8981144288438190D+00, &
          0.1018855711561810D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.6731903387258514D+00, &
          0.3268096612741486D+00, &
          0.3268096612741486D+00, &
          0.6731903387258514D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.6731903387258514D+00, &
          0.3268096612741486D+00, &
          0.3268096612741486D+00, &
          0.6731903387258514D+00, &
          0.8467830886160004D+00, &
          0.1532169113839996D+00, &
          0.8467830886160004D+00, &
          0.1532169113839996D+00, &
          0.8467830886160004D+00, &
          0.1532169113839996D+00, &
          0.8467830886160004D+00, &
          0.1532169113839996D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.1873892499846415D+00, &
          0.8126107500153585D+00, &
          0.8126107500153585D+00, &
          0.1873892499846415D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.5000000000000000D+00, &
          0.1873892499846415D+00, &
          0.8126107500153585D+00, &
          0.8126107500153585D+00, &
          0.1873892499846415D+00, &
          0.9359181010210101D+00, &
          0.0640818989789899D+00, &
          0.9359181010210101D+00, &
          0.0640818989789899D+00, &
          0.9359181010210101D+00, &
          0.0640818989789899D+00, &
          0.9359181010210101D+00, &
          0.0640818989789899D+00, &
          0.2917727295202502D+00, &
          0.2917727295202502D+00, &
          0.9026740712514576D+00, &
          0.2917727295202502D+00, &
          0.2917727295202502D+00, &
          0.9026740712514576D+00, &
          0.7082272704797498D+00, &
          0.7082272704797498D+00, &
          0.9026740712514576D+00, &
          0.7082272704797498D+00, &
          0.7082272704797498D+00, &
          0.9026740712514576D+00, &
          0.2917727295202502D+00, &
          0.2917727295202502D+00, &
          0.0973259287485423D+00, &
          0.2917727295202502D+00, &
          0.2917727295202502D+00, &
          0.0973259287485423D+00, &
          0.7082272704797498D+00, &
          0.7082272704797498D+00, &
          0.0973259287485423D+00, &
          0.7082272704797498D+00, &
          0.7082272704797498D+00, &
          0.0973259287485423D+00, &
          0.9869259796545771D+00, &
          0.9869259796545771D+00, &
          0.9372501994187228D+00, &
          0.9869259796545771D+00, &
          0.9869259796545771D+00, &
          0.9372501994187228D+00, &
          0.0130740203454229D+00, &
          0.0130740203454229D+00, &
          0.9372501994187228D+00, &
          0.0130740203454229D+00, &
          0.0130740203454229D+00, &
          0.9372501994187228D+00, &
          0.9869259796545771D+00, &
          0.9869259796545771D+00, &
          0.0627498005812772D+00, &
          0.9869259796545771D+00, &
          0.9869259796545771D+00, &
          0.0627498005812772D+00, &
          0.0130740203454229D+00, &
          0.0130740203454229D+00, &
          0.0627498005812772D+00, &
          0.0130740203454229D+00, &
          0.0130740203454229D+00, &
          0.0627498005812772D+00, &
          0.8645117252804405D+00, &
          0.8645117252804405D+00, &
          0.6225837716647646D+00, &
          0.8645117252804405D+00, &
          0.8645117252804405D+00, &
          0.6225837716647646D+00, &
          0.1354882747195595D+00, &
          0.1354882747195595D+00, &
          0.6225837716647646D+00, &
          0.1354882747195595D+00, &
          0.1354882747195595D+00, &
          0.6225837716647646D+00, &
          0.8645117252804405D+00, &
          0.8645117252804405D+00, &
          0.3774162283352354D+00, &
          0.8645117252804405D+00, &
          0.8645117252804405D+00, &
          0.3774162283352354D+00, &
          0.1354882747195595D+00, &
          0.1354882747195595D+00, &
          0.3774162283352354D+00, &
          0.1354882747195595D+00, &
          0.1354882747195595D+00, &
          0.3774162283352354D+00, &
          0.8248058551406905D+00, &
          0.8248058551406905D+00, &
          0.9220515632141263D+00, &
          0.8248058551406905D+00, &
          0.8248058551406905D+00, &
          0.9220515632141263D+00, &
          0.1751941448593096D+00, &
          0.1751941448593096D+00, &
          0.9220515632141263D+00, &
          0.1751941448593096D+00, &
          0.1751941448593096D+00, &
          0.9220515632141263D+00, &
          0.8248058551406905D+00, &
          0.8248058551406905D+00, &
          0.0779484367858737D+00, &
          0.8248058551406905D+00, &
          0.8248058551406905D+00, &
          0.0779484367858737D+00, &
          0.1751941448593096D+00, &
          0.1751941448593096D+00, &
          0.0779484367858737D+00, &
          0.1751941448593096D+00, &
          0.1751941448593096D+00, &
          0.0779484367858737D+00, &
          0.9459853758476011D+00, &
          0.9459853758476011D+00, &
          0.5650967942783419D+00, &
          0.9459853758476011D+00, &
          0.9459853758476011D+00, &
          0.5650967942783419D+00, &
          0.0540146241523989D+00, &
          0.0540146241523989D+00, &
          0.5650967942783419D+00, &
          0.0540146241523989D+00, &
          0.0540146241523989D+00, &
          0.5650967942783419D+00, &
          0.9459853758476011D+00, &
          0.9459853758476011D+00, &
          0.4349032057216580D+00, &
          0.9459853758476011D+00, &
          0.9459853758476011D+00, &
          0.4349032057216580D+00, &
          0.0540146241523989D+00, &
          0.0540146241523989D+00, &
          0.4349032057216580D+00, &
          0.0540146241523989D+00, &
          0.0540146241523989D+00, &
          0.4349032057216580D+00, &
          0.8916635718645309D+00, &
          0.8916635718645309D+00, &
          0.9766093536267083D+00, &
          0.8916635718645309D+00, &
          0.8916635718645309D+00, &
          0.9766093536267083D+00, &
          0.1083364281354691D+00, &
          0.1083364281354691D+00, &
          0.9766093536267083D+00, &
          0.1083364281354691D+00, &
          0.1083364281354691D+00, &
          0.9766093536267083D+00, &
          0.8916635718645309D+00, &
          0.8916635718645309D+00, &
          0.0233906463732917D+00, &
          0.8916635718645309D+00, &
          0.8916635718645309D+00, &
          0.0233906463732917D+00, &
          0.1083364281354691D+00, &
          0.1083364281354691D+00, &
          0.0233906463732917D+00, &
          0.1083364281354691D+00, &
          0.1083364281354691D+00, &
          0.0233906463732917D+00, &
          0.9849705936644424D+00, &
          0.9849705936644424D+00, &
          0.6355079664366210D+00, &
          0.9849705936644424D+00, &
          0.9849705936644424D+00, &
          0.6355079664366210D+00, &
          0.0150294063355575D+00, &
          0.0150294063355575D+00, &
          0.6355079664366210D+00, &
          0.0150294063355575D+00, &
          0.0150294063355575D+00, &
          0.6355079664366210D+00, &
          0.9849705936644424D+00, &
          0.9849705936644424D+00, &
          0.3644920335633790D+00, &
          0.9849705936644424D+00, &
          0.9849705936644424D+00, &
          0.3644920335633790D+00, &
          0.0150294063355575D+00, &
          0.0150294063355575D+00, &
          0.3644920335633790D+00, &
          0.0150294063355575D+00, &
          0.0150294063355575D+00, &
          0.3644920335633790D+00, &
          0.5860169886103869D+00, &
          0.5860169886103869D+00, &
          0.9923845018898745D+00, &
          0.5860169886103869D+00, &
          0.5860169886103869D+00, &
          0.9923845018898745D+00, &
          0.4139830113896131D+00, &
          0.4139830113896131D+00, &
          0.9923845018898745D+00, &
          0.4139830113896131D+00, &
          0.4139830113896131D+00, &
          0.9923845018898745D+00, &
          0.5860169886103869D+00, &
          0.5860169886103869D+00, &
          0.0076154981101255D+00, &
          0.5860169886103869D+00, &
          0.5860169886103869D+00, &
          0.0076154981101255D+00, &
          0.4139830113896131D+00, &
          0.4139830113896131D+00, &
          0.0076154981101255D+00, &
          0.4139830113896131D+00, &
          0.4139830113896131D+00, &
          0.0076154981101255D+00, &
          0.7690492963095925D+00, &
          0.7690492963095925D+00, &
          0.3896184326892077D+00, &
          0.7690492963095925D+00, &
          0.7690492963095925D+00, &
          0.3896184326892077D+00, &
          0.2309507036904074D+00, &
          0.2309507036904074D+00, &
          0.3896184326892077D+00, &
          0.2309507036904074D+00, &
          0.2309507036904074D+00, &
          0.3896184326892077D+00, &
          0.7690492963095925D+00, &
          0.7690492963095925D+00, &
          0.6103815673107923D+00, &
          0.7690492963095925D+00, &
          0.7690492963095925D+00, &
          0.6103815673107923D+00, &
          0.2309507036904074D+00, &
          0.2309507036904074D+00, &
          0.6103815673107923D+00, &
          0.2309507036904074D+00, &
          0.2309507036904074D+00, &
          0.6103815673107923D+00, &
          0.7768424517946938D+00, &
          0.7768424517946938D+00, &
          0.9994864489115713D+00, &
          0.7768424517946938D+00, &
          0.7768424517946938D+00, &
          0.9994864489115713D+00, &
          0.2231575482053062D+00, &
          0.2231575482053062D+00, &
          0.9994864489115713D+00, &
          0.2231575482053062D+00, &
          0.2231575482053062D+00, &
          0.9994864489115713D+00, &
          0.7768424517946938D+00, &
          0.7768424517946938D+00, &
          0.0005135510884287D+00, &
          0.7768424517946938D+00, &
          0.7768424517946938D+00, &
          0.0005135510884287D+00, &
          0.2231575482053062D+00, &
          0.2231575482053062D+00, &
          0.0005135510884287D+00, &
          0.2231575482053062D+00, &
          0.2231575482053062D+00, &
          0.0005135510884287D+00, &
          0.6117790069755276D+00, &
          0.6117790069755276D+00, &
          0.7209288728023533D+00, &
          0.6117790069755276D+00, &
          0.6117790069755276D+00, &
          0.7209288728023533D+00, &
          0.3882209930244724D+00, &
          0.3882209930244724D+00, &
          0.7209288728023533D+00, &
          0.3882209930244724D+00, &
          0.3882209930244724D+00, &
          0.7209288728023533D+00, &
          0.6117790069755276D+00, &
          0.6117790069755276D+00, &
          0.2790711271976467D+00, &
          0.6117790069755276D+00, &
          0.6117790069755276D+00, &
          0.2790711271976467D+00, &
          0.3882209930244724D+00, &
          0.3882209930244724D+00, &
          0.2790711271976467D+00, &
          0.3882209930244724D+00, &
          0.3882209930244724D+00, &
          0.2790711271976467D+00, &
          0.9350312907963207D+00, &
          0.9350312907963207D+00, &
          0.7494073604468848D+00, &
          0.9350312907963207D+00, &
          0.9350312907963207D+00, &
          0.7494073604468848D+00, &
          0.0649687092036792D+00, &
          0.0649687092036792D+00, &
          0.7494073604468848D+00, &
          0.0649687092036792D+00, &
          0.0649687092036792D+00, &
          0.7494073604468848D+00, &
          0.9350312907963207D+00, &
          0.9350312907963207D+00, &
          0.2505926395531152D+00, &
          0.9350312907963207D+00, &
          0.9350312907963207D+00, &
          0.2505926395531152D+00, &
          0.0649687092036792D+00, &
          0.0649687092036792D+00, &
          0.2505926395531152D+00, &
          0.0649687092036792D+00, &
          0.0649687092036792D+00, &
          0.2505926395531152D+00, &
          0.6986057011288138D+00, &
          0.9719926774632096D+00, &
          0.6986057011288138D+00, &
          0.9719926774632096D+00, &
          0.8750719993006615D+00, &
          0.8750719993006615D+00, &
          0.6986057011288138D+00, &
          0.9719926774632096D+00, &
          0.6986057011288138D+00, &
          0.9719926774632096D+00, &
          0.1249280006993386D+00, &
          0.1249280006993386D+00, &
          0.6986057011288138D+00, &
          0.0280073225367904D+00, &
          0.6986057011288138D+00, &
          0.0280073225367904D+00, &
          0.8750719993006615D+00, &
          0.8750719993006615D+00, &
          0.3013942988711862D+00, &
          0.9719926774632096D+00, &
          0.3013942988711862D+00, &
          0.9719926774632096D+00, &
          0.8750719993006615D+00, &
          0.8750719993006615D+00, &
          0.6986057011288138D+00, &
          0.0280073225367904D+00, &
          0.6986057011288138D+00, &
          0.0280073225367904D+00, &
          0.1249280006993386D+00, &
          0.1249280006993386D+00, &
          0.3013942988711862D+00, &
          0.9719926774632096D+00, &
          0.3013942988711862D+00, &
          0.9719926774632096D+00, &
          0.1249280006993386D+00, &
          0.1249280006993386D+00, &
          0.3013942988711862D+00, &
          0.0280073225367904D+00, &
          0.3013942988711862D+00, &
          0.0280073225367904D+00, &
          0.8750719993006615D+00, &
          0.8750719993006615D+00, &
          0.3013942988711862D+00, &
          0.0280073225367904D+00, &
          0.3013942988711862D+00, &
          0.0280073225367904D+00, &
          0.1249280006993386D+00, &
          0.1249280006993386D+00, &
          0.8153563222457214D+00, &
          0.9923118941018740D+00, &
          0.8153563222457214D+00, &
          0.9923118941018740D+00, &
          0.9512929061479376D+00, &
          0.9512929061479376D+00, &
          0.8153563222457214D+00, &
          0.9923118941018740D+00, &
          0.8153563222457214D+00, &
          0.9923118941018740D+00, &
          0.0487070938520624D+00, &
          0.0487070938520624D+00, &
          0.8153563222457214D+00, &
          0.0076881058981259D+00, &
          0.8153563222457214D+00, &
          0.0076881058981259D+00, &
          0.9512929061479376D+00, &
          0.9512929061479376D+00, &
          0.1846436777542786D+00, &
          0.9923118941018740D+00, &
          0.1846436777542786D+00, &
          0.9923118941018740D+00, &
          0.9512929061479376D+00, &
          0.9512929061479376D+00, &
          0.8153563222457214D+00, &
          0.0076881058981259D+00, &
          0.8153563222457214D+00, &
          0.0076881058981259D+00, &
          0.0487070938520624D+00, &
          0.0487070938520624D+00, &
          0.1846436777542786D+00, &
          0.9923118941018740D+00, &
          0.1846436777542786D+00, &
          0.9923118941018740D+00, &
          0.0487070938520624D+00, &
          0.0487070938520624D+00, &
          0.1846436777542786D+00, &
          0.0076881058981259D+00, &
          0.1846436777542786D+00, &
          0.0076881058981259D+00, &
          0.9512929061479376D+00, &
          0.9512929061479376D+00, &
          0.1846436777542786D+00, &
          0.0076881058981259D+00, &
          0.1846436777542786D+00, &
          0.0076881058981259D+00, &
          0.0487070938520624D+00, &
          0.0487070938520624D+00, &
          0.9360994289291793D+00, &
          0.5491142465864416D+00, &
          0.9360994289291793D+00, &
          0.5491142465864416D+00, &
          0.6108897298239304D+00, &
          0.6108897298239304D+00, &
          0.9360994289291793D+00, &
          0.5491142465864416D+00, &
          0.9360994289291793D+00, &
          0.5491142465864416D+00, &
          0.3891102701760696D+00, &
          0.3891102701760696D+00, &
          0.9360994289291793D+00, &
          0.4508857534135584D+00, &
          0.9360994289291793D+00, &
          0.4508857534135584D+00, &
          0.6108897298239304D+00, &
          0.6108897298239304D+00, &
          0.0639005710708208D+00, &
          0.5491142465864416D+00, &
          0.0639005710708208D+00, &
          0.5491142465864416D+00, &
          0.6108897298239304D+00, &
          0.6108897298239304D+00, &
          0.9360994289291793D+00, &
          0.4508857534135584D+00, &
          0.9360994289291793D+00, &
          0.4508857534135584D+00, &
          0.3891102701760696D+00, &
          0.3891102701760696D+00, &
          0.0639005710708208D+00, &
          0.5491142465864416D+00, &
          0.0639005710708208D+00, &
          0.5491142465864416D+00, &
          0.3891102701760696D+00, &
          0.3891102701760696D+00, &
          0.0639005710708208D+00, &
          0.4508857534135584D+00, &
          0.0639005710708208D+00, &
          0.4508857534135584D+00, &
          0.6108897298239304D+00, &
          0.6108897298239304D+00, &
          0.0639005710708208D+00, &
          0.4508857534135584D+00, &
          0.0639005710708208D+00, &
          0.4508857534135584D+00, &
          0.3891102701760696D+00, &
          0.3891102701760696D+00, &
          0.7479204633014319D+00, &
          0.6040696445748173D+00, &
          0.7479204633014319D+00, &
          0.6040696445748173D+00, &
          0.9784211281771267D+00, &
          0.9784211281771267D+00, &
          0.7479204633014319D+00, &
          0.6040696445748173D+00, &
          0.7479204633014319D+00, &
          0.6040696445748173D+00, &
          0.0215788718228733D+00, &
          0.0215788718228733D+00, &
          0.7479204633014319D+00, &
          0.3959303554251827D+00, &
          0.7479204633014319D+00, &
          0.3959303554251827D+00, &
          0.9784211281771267D+00, &
          0.9784211281771267D+00, &
          0.2520795366985681D+00, &
          0.6040696445748173D+00, &
          0.2520795366985681D+00, &
          0.6040696445748173D+00, &
          0.9784211281771267D+00, &
          0.9784211281771267D+00, &
          0.7479204633014319D+00, &
          0.3959303554251827D+00, &
          0.7479204633014319D+00, &
          0.3959303554251827D+00, &
          0.0215788718228733D+00, &
          0.0215788718228733D+00, &
          0.2520795366985681D+00, &
          0.6040696445748173D+00, &
          0.2520795366985681D+00, &
          0.6040696445748173D+00, &
          0.0215788718228733D+00, &
          0.0215788718228733D+00, &
          0.2520795366985681D+00, &
          0.3959303554251827D+00, &
          0.2520795366985681D+00, &
          0.3959303554251827D+00, &
          0.9784211281771267D+00, &
          0.9784211281771267D+00, &
          0.2520795366985681D+00, &
          0.3959303554251827D+00, &
          0.2520795366985681D+00, &
          0.3959303554251827D+00, &
          0.0215788718228733D+00, &
          0.0215788718228733D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          0.0030079175616226D+00, &
          0.0030079175616226D+00, &
          0.0030079175616226D+00, &
          0.0030079175616226D+00, &
          0.0030079175616226D+00, &
          0.0030079175616226D+00, &
          0.0061575518559863D+00, &
          0.0061575518559863D+00, &
          0.0061575518559863D+00, &
          0.0061575518559863D+00, &
          0.0061575518559863D+00, &
          0.0061575518559863D+00, &
          0.0007451864938265D+00, &
          0.0007451864938265D+00, &
          0.0007451864938265D+00, &
          0.0007451864938265D+00, &
          0.0007451864938265D+00, &
          0.0007451864938265D+00, &
          0.0007451864938265D+00, &
          0.0007451864938265D+00, &
          0.0034736498758246D+00, &
          0.0034736498758246D+00, &
          0.0034736498758246D+00, &
          0.0034736498758246D+00, &
          0.0034736498758246D+00, &
          0.0034736498758246D+00, &
          0.0034736498758246D+00, &
          0.0034736498758246D+00, &
          0.0005231065457280D+00, &
          0.0005231065457280D+00, &
          0.0005231065457280D+00, &
          0.0005231065457280D+00, &
          0.0005231065457280D+00, &
          0.0005231065457280D+00, &
          0.0005231065457280D+00, &
          0.0005231065457280D+00, &
          0.0005231065457280D+00, &
          0.0005231065457280D+00, &
          0.0005231065457280D+00, &
          0.0005231065457280D+00, &
          0.0005231065457280D+00, &
          0.0005231065457280D+00, &
          0.0005231065457280D+00, &
          0.0005231065457280D+00, &
          0.0005231065457280D+00, &
          0.0005231065457280D+00, &
          0.0005231065457280D+00, &
          0.0005231065457280D+00, &
          0.0005231065457280D+00, &
          0.0005231065457280D+00, &
          0.0005231065457280D+00, &
          0.0005231065457280D+00, &
          0.0047212614224319D+00, &
          0.0047212614224319D+00, &
          0.0047212614224319D+00, &
          0.0047212614224319D+00, &
          0.0047212614224319D+00, &
          0.0047212614224319D+00, &
          0.0047212614224319D+00, &
          0.0047212614224319D+00, &
          0.0047212614224319D+00, &
          0.0047212614224319D+00, &
          0.0047212614224319D+00, &
          0.0047212614224319D+00, &
          0.0047212614224319D+00, &
          0.0047212614224319D+00, &
          0.0047212614224319D+00, &
          0.0047212614224319D+00, &
          0.0047212614224319D+00, &
          0.0047212614224319D+00, &
          0.0047212614224319D+00, &
          0.0047212614224319D+00, &
          0.0047212614224319D+00, &
          0.0047212614224319D+00, &
          0.0047212614224319D+00, &
          0.0047212614224319D+00, &
          0.0023372223737635D+00, &
          0.0023372223737635D+00, &
          0.0023372223737635D+00, &
          0.0023372223737635D+00, &
          0.0023372223737635D+00, &
          0.0023372223737635D+00, &
          0.0023372223737635D+00, &
          0.0023372223737635D+00, &
          0.0023372223737635D+00, &
          0.0023372223737635D+00, &
          0.0023372223737635D+00, &
          0.0023372223737635D+00, &
          0.0023372223737635D+00, &
          0.0023372223737635D+00, &
          0.0023372223737635D+00, &
          0.0023372223737635D+00, &
          0.0023372223737635D+00, &
          0.0023372223737635D+00, &
          0.0023372223737635D+00, &
          0.0023372223737635D+00, &
          0.0023372223737635D+00, &
          0.0023372223737635D+00, &
          0.0023372223737635D+00, &
          0.0023372223737635D+00, &
          0.0038417548985660D+00, &
          0.0038417548985660D+00, &
          0.0038417548985660D+00, &
          0.0038417548985660D+00, &
          0.0038417548985660D+00, &
          0.0038417548985660D+00, &
          0.0038417548985660D+00, &
          0.0038417548985660D+00, &
          0.0038417548985660D+00, &
          0.0038417548985660D+00, &
          0.0038417548985660D+00, &
          0.0038417548985660D+00, &
          0.0038417548985660D+00, &
          0.0038417548985660D+00, &
          0.0038417548985660D+00, &
          0.0038417548985660D+00, &
          0.0038417548985660D+00, &
          0.0038417548985660D+00, &
          0.0038417548985660D+00, &
          0.0038417548985660D+00, &
          0.0038417548985660D+00, &
          0.0038417548985660D+00, &
          0.0038417548985660D+00, &
          0.0038417548985660D+00, &
          0.0001507767719634D+00, &
          0.0001507767719634D+00, &
          0.0001507767719634D+00, &
          0.0001507767719634D+00, &
          0.0001507767719634D+00, &
          0.0001507767719634D+00, &
          0.0001507767719634D+00, &
          0.0001507767719634D+00, &
          0.0001507767719634D+00, &
          0.0001507767719634D+00, &
          0.0001507767719634D+00, &
          0.0001507767719634D+00, &
          0.0001507767719634D+00, &
          0.0001507767719634D+00, &
          0.0001507767719634D+00, &
          0.0001507767719634D+00, &
          0.0001507767719634D+00, &
          0.0001507767719634D+00, &
          0.0001507767719634D+00, &
          0.0001507767719634D+00, &
          0.0001507767719634D+00, &
          0.0001507767719634D+00, &
          0.0001507767719634D+00, &
          0.0001507767719634D+00, &
          0.0032140794095183D+00, &
          0.0032140794095183D+00, &
          0.0032140794095183D+00, &
          0.0032140794095183D+00, &
          0.0032140794095183D+00, &
          0.0032140794095183D+00, &
          0.0032140794095183D+00, &
          0.0032140794095183D+00, &
          0.0032140794095183D+00, &
          0.0032140794095183D+00, &
          0.0032140794095183D+00, &
          0.0032140794095183D+00, &
          0.0032140794095183D+00, &
          0.0032140794095183D+00, &
          0.0032140794095183D+00, &
          0.0032140794095183D+00, &
          0.0032140794095183D+00, &
          0.0032140794095183D+00, &
          0.0032140794095183D+00, &
          0.0032140794095183D+00, &
          0.0032140794095183D+00, &
          0.0032140794095183D+00, &
          0.0032140794095183D+00, &
          0.0032140794095183D+00, &
          0.0022454466226130D+00, &
          0.0022454466226130D+00, &
          0.0022454466226130D+00, &
          0.0022454466226130D+00, &
          0.0022454466226130D+00, &
          0.0022454466226130D+00, &
          0.0022454466226130D+00, &
          0.0022454466226130D+00, &
          0.0022454466226130D+00, &
          0.0022454466226130D+00, &
          0.0022454466226130D+00, &
          0.0022454466226130D+00, &
          0.0022454466226130D+00, &
          0.0022454466226130D+00, &
          0.0022454466226130D+00, &
          0.0022454466226130D+00, &
          0.0022454466226130D+00, &
          0.0022454466226130D+00, &
          0.0022454466226130D+00, &
          0.0022454466226130D+00, &
          0.0022454466226130D+00, &
          0.0022454466226130D+00, &
          0.0022454466226130D+00, &
          0.0022454466226130D+00, &
          0.0010413676731518D+00, &
          0.0010413676731518D+00, &
          0.0010413676731518D+00, &
          0.0010413676731518D+00, &
          0.0010413676731518D+00, &
          0.0010413676731518D+00, &
          0.0010413676731518D+00, &
          0.0010413676731518D+00, &
          0.0010413676731518D+00, &
          0.0010413676731518D+00, &
          0.0010413676731518D+00, &
          0.0010413676731518D+00, &
          0.0010413676731518D+00, &
          0.0010413676731518D+00, &
          0.0010413676731518D+00, &
          0.0010413676731518D+00, &
          0.0010413676731518D+00, &
          0.0010413676731518D+00, &
          0.0010413676731518D+00, &
          0.0010413676731518D+00, &
          0.0010413676731518D+00, &
          0.0010413676731518D+00, &
          0.0010413676731518D+00, &
          0.0010413676731518D+00, &
          0.0005920000965024D+00, &
          0.0005920000965024D+00, &
          0.0005920000965024D+00, &
          0.0005920000965024D+00, &
          0.0005920000965024D+00, &
          0.0005920000965024D+00, &
          0.0005920000965024D+00, &
          0.0005920000965024D+00, &
          0.0005920000965024D+00, &
          0.0005920000965024D+00, &
          0.0005920000965024D+00, &
          0.0005920000965024D+00, &
          0.0005920000965024D+00, &
          0.0005920000965024D+00, &
          0.0005920000965024D+00, &
          0.0005920000965024D+00, &
          0.0005920000965024D+00, &
          0.0005920000965024D+00, &
          0.0005920000965024D+00, &
          0.0005920000965024D+00, &
          0.0005920000965024D+00, &
          0.0005920000965024D+00, &
          0.0005920000965024D+00, &
          0.0005920000965024D+00, &
          0.0004100172337698D+00, &
          0.0004100172337698D+00, &
          0.0004100172337698D+00, &
          0.0004100172337698D+00, &
          0.0004100172337698D+00, &
          0.0004100172337698D+00, &
          0.0004100172337698D+00, &
          0.0004100172337698D+00, &
          0.0004100172337698D+00, &
          0.0004100172337698D+00, &
          0.0004100172337698D+00, &
          0.0004100172337698D+00, &
          0.0004100172337698D+00, &
          0.0004100172337698D+00, &
          0.0004100172337698D+00, &
          0.0004100172337698D+00, &
          0.0004100172337698D+00, &
          0.0004100172337698D+00, &
          0.0004100172337698D+00, &
          0.0004100172337698D+00, &
          0.0004100172337698D+00, &
          0.0004100172337698D+00, &
          0.0004100172337698D+00, &
          0.0004100172337698D+00, &
          0.0007465850926621D+00, &
          0.0007465850926621D+00, &
          0.0007465850926621D+00, &
          0.0007465850926621D+00, &
          0.0007465850926621D+00, &
          0.0007465850926621D+00, &
          0.0007465850926621D+00, &
          0.0007465850926621D+00, &
          0.0007465850926621D+00, &
          0.0007465850926621D+00, &
          0.0007465850926621D+00, &
          0.0007465850926621D+00, &
          0.0007465850926621D+00, &
          0.0007465850926621D+00, &
          0.0007465850926621D+00, &
          0.0007465850926621D+00, &
          0.0007465850926621D+00, &
          0.0007465850926621D+00, &
          0.0007465850926621D+00, &
          0.0007465850926621D+00, &
          0.0007465850926621D+00, &
          0.0007465850926621D+00, &
          0.0007465850926621D+00, &
          0.0007465850926621D+00, &
          0.0037432896839704D+00, &
          0.0037432896839704D+00, &
          0.0037432896839704D+00, &
          0.0037432896839704D+00, &
          0.0037432896839704D+00, &
          0.0037432896839704D+00, &
          0.0037432896839704D+00, &
          0.0037432896839704D+00, &
          0.0037432896839704D+00, &
          0.0037432896839704D+00, &
          0.0037432896839704D+00, &
          0.0037432896839704D+00, &
          0.0037432896839704D+00, &
          0.0037432896839704D+00, &
          0.0037432896839704D+00, &
          0.0037432896839704D+00, &
          0.0037432896839704D+00, &
          0.0037432896839704D+00, &
          0.0037432896839704D+00, &
          0.0037432896839704D+00, &
          0.0037432896839704D+00, &
          0.0037432896839704D+00, &
          0.0037432896839704D+00, &
          0.0037432896839704D+00, &
          0.0004534357504253D+00, &
          0.0004534357504253D+00, &
          0.0004534357504253D+00, &
          0.0004534357504253D+00, &
          0.0004534357504253D+00, &
          0.0004534357504253D+00, &
          0.0004534357504253D+00, &
          0.0004534357504253D+00, &
          0.0004534357504253D+00, &
          0.0004534357504253D+00, &
          0.0004534357504253D+00, &
          0.0004534357504253D+00, &
          0.0004534357504253D+00, &
          0.0004534357504253D+00, &
          0.0004534357504253D+00, &
          0.0004534357504253D+00, &
          0.0004534357504253D+00, &
          0.0004534357504253D+00, &
          0.0004534357504253D+00, &
          0.0004534357504253D+00, &
          0.0004534357504253D+00, &
          0.0004534357504253D+00, &
          0.0004534357504253D+00, &
          0.0004534357504253D+00, &
          0.0041828675781781D+00, &
          0.0041828675781781D+00, &
          0.0041828675781781D+00, &
          0.0041828675781781D+00, &
          0.0041828675781781D+00, &
          0.0041828675781781D+00, &
          0.0041828675781781D+00, &
          0.0041828675781781D+00, &
          0.0041828675781781D+00, &
          0.0041828675781781D+00, &
          0.0041828675781781D+00, &
          0.0041828675781781D+00, &
          0.0041828675781781D+00, &
          0.0041828675781781D+00, &
          0.0041828675781781D+00, &
          0.0041828675781781D+00, &
          0.0041828675781781D+00, &
          0.0041828675781781D+00, &
          0.0041828675781781D+00, &
          0.0041828675781781D+00, &
          0.0041828675781781D+00, &
          0.0041828675781781D+00, &
          0.0041828675781781D+00, &
          0.0041828675781781D+00, &
          0.0007389797381823D+00, &
          0.0007389797381823D+00, &
          0.0007389797381823D+00, &
          0.0007389797381823D+00, &
          0.0007389797381823D+00, &
          0.0007389797381823D+00, &
          0.0007389797381823D+00, &
          0.0007389797381823D+00, &
          0.0007389797381823D+00, &
          0.0007389797381823D+00, &
          0.0007389797381823D+00, &
          0.0007389797381823D+00, &
          0.0007389797381823D+00, &
          0.0007389797381823D+00, &
          0.0007389797381823D+00, &
          0.0007389797381823D+00, &
          0.0007389797381823D+00, &
          0.0007389797381823D+00, &
          0.0007389797381823D+00, &
          0.0007389797381823D+00, &
          0.0007389797381823D+00, &
          0.0007389797381823D+00, &
          0.0007389797381823D+00, &
          0.0007389797381823D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0013712825533878D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0003675472168789D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0012580194835158D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00, &
          0.0015165655616946D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
