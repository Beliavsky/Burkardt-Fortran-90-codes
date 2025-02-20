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
function r8mat_det_4d ( a )

!*****************************************************************************80
!
!! r8mat_det_4d() computes the determinant of a 4 by 4 matrix.
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
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) A(4,4): the matrix whose determinant is desired.
!
!  Output:
!
!    real ( kind = rk ) R8MAT_DET_4D: the determinant of the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(4,4)
  real ( kind = rk ) r8mat_det_4d

  r8mat_det_4d = &
         a(1,1) * ( &
             a(2,2) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
           - a(2,3) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
           + a(2,4) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) ) &
       - a(1,2) * ( &
             a(2,1) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
           - a(2,3) * ( a(3,1) * a(4,4) - a(3,4) * a(4,1) ) &
           + a(2,4) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) ) &
       + a(1,3) * ( &
             a(2,1) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
           - a(2,2) * ( a(3,1) * a(4,4) - a(3,4) * a(4,1) ) &
           + a(2,4) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) ) &
       - a(1,4) * ( &
             a(2,1) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) &
           - a(2,2) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) &
           + a(2,3) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) )

  return
end
subroutine rule_order ( p, order )

!*****************************************************************************80
!
!! rule_order() returns the order of a tetrahedron quadrature rule of given precision.
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
!    24 April 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Freddie Witherden, Peter Vincent,
!    On the identification of symmetric quadrature rules for finite element methods,
!    Computers and Mathematics with Applications,
!    Volume 69, pages 1232-1241, 2015.
!
!  Input:
!
!    integer p: the precision, 0 <= p <= 10.
!
!  Output:
!
!    integer order: the order of the rule.
!
  implicit none

  integer order
  integer, save, dimension ( 0:10 ) :: order_vec =  (/ &
    1, &
    1,   4,   8,  14,  14, &
   24,  35,  46,  59,  81 /)
  integer p

  if ( p < 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'rule_order(): Fatal error!'
    write ( *, '(a)' ) '  Input p < 0.'
    call exit ( 1 )
  end if

  if ( 10 < p ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'rule_order(): Fatal error!'
    write ( *, '(a)' ) '  Input 10 < p.'
    call exit ( 1 )
  end if

  order = order_vec(p)

  return
end
subroutine rule00 ( n, a, b, c, d, w )

!*****************************************************************************80
!
!! rule00() returns the rule of precision 0.
!
!  Discussion:
!
!    The data is given for the reference tetrahedron
!    with vertices (0,0,0), (1,0,0), (0,1,0), (0,0,1).
!
!    We suppose we are given a tetrahedron T with vertices A, B, C, D.
!    We call a rule with n points, returning barycentric coordinates
!    a, b, c, d, and weights w.  Then the integral I of f(x,y,z) over T is 
!    approximated by Q as follows:
!
!    (x,y,z) = a(1:n) * A + b(1:n) * B + c(1:n) * C + d(1:n * D
!    Q = area(T) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 April 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Freddie Witherden, Peter Vincent,
!    On the identification of symmetric quadrature rules for finite element methods,
!    Computers and Mathematics with Applications,
!    Volume 69, pages 1232-1241, 2015.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real a(n), b(n), c(n), d(n): the barycentric coordinates of quadrature points.
!
!    real w(n): the quadrature weights, which sum to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 1

  real ( kind = rk ) a(n)
  real ( kind = rk ), save, dimension ( n_save ) :: a_save = (/  &
      0.25000000000000000000D+00 /)

  real ( kind = rk ) b(n)
  real ( kind = rk ), save, dimension ( n_save ) :: b_save = (/ &
      0.25000000000000000000D+00 /)

  real ( kind = rk ) c(n)
  real ( kind = rk ), save, dimension ( n_save ) :: c_save = (/ &
      0.25000000000000000000D+00 /)

  real ( kind = rk ) d(n)
  real ( kind = rk ), save, dimension ( n_save ) :: d_save = (/ &
      0.25000000000000000000D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
      1.00000000000000000000D+00 /)

  a(1:n) = a_save(1:n)
  b(1:n) = b_save(1:n)
  c(1:n) = c_save(1:n)
  d(1:n) = d_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule01 ( n, a, b, c, d, w )

!*****************************************************************************80
!
!! rule01() returns the rule of precision 1.
!
!  Discussion:
!
!    The data is given for the reference tetrahedron
!    with vertices (0,0,0), (1,0,0), (0,1,0), (0,0,1).
!
!    We suppose we are given a tetrahedron T with vertices A, B, C, D.
!    We call a rule with n points, returning barycentric coordinates
!    a, b, c, d, and weights w.  Then the integral I of f(x,y,z) over T is 
!    approximated by Q as follows:
!
!    (x,y,z) = a(1:n) * A + b(1:n) * B + c(1:n) * C + d(1:n * D
!    Q = area(T) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 April 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Freddie Witherden, Peter Vincent,
!    On the identification of symmetric quadrature rules for finite element methods,
!    Computers and Mathematics with Applications,
!    Volume 69, pages 1232-1241, 2015.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real a(n), b(n), c(n), d(n): the barycentric coordinates of quadrature points.
!
!    real w(n): the quadrature weights, which sum to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 1

  real ( kind = rk ) a(n)
  real ( kind = rk ), save, dimension ( n_save ) :: a_save = (/  &
      0.25000000000000000000D+00 /)

  real ( kind = rk ) b(n)
  real ( kind = rk ), save, dimension ( n_save ) :: b_save = (/ &
      0.25000000000000000000D+00 /)

  real ( kind = rk ) c(n)
  real ( kind = rk ), save, dimension ( n_save ) :: c_save = (/ &
      0.25000000000000000000D+00 /)

  real ( kind = rk ) d(n)
  real ( kind = rk ), save, dimension ( n_save ) :: d_save = (/ &
      0.25000000000000000000D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
      1.00000000000000000000D+00 /)

  a(1:n) = a_save(1:n)
  b(1:n) = b_save(1:n)
  c(1:n) = c_save(1:n)
  d(1:n) = d_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule02 ( n, a, b, c, d, w )

!*****************************************************************************80
!
!! rule02() returns the rule of precision 2.
!
!  Discussion:
!
!    The data is given for the reference tetrahedron
!    with vertices (0,0,0), (1,0,0), (0,1,0), (0,0,1).
!
!    We suppose we are given a tetrahedron T with vertices A, B, C, D.
!    We call a rule with n points, returning barycentric coordinates
!    a, b, c, d, and weights w.  Then the integral I of f(x,y,z) over T is 
!    approximated by Q as follows:
!
!    (x,y,z) = a(1:n) * A + b(1:n) * B + c(1:n) * C + d(1:n * D
!    Q = area(T) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 April 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Freddie Witherden, Peter Vincent,
!    On the identification of symmetric quadrature rules for finite element methods,
!    Computers and Mathematics with Applications,
!    Volume 69, pages 1232-1241, 2015.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real a(n), b(n), c(n), d(n): the barycentric coordinates of quadrature points.
!
!    real w(n): the quadrature weights, which sum to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 4

  real ( kind = rk ) a(n)
  real ( kind = rk ), save, dimension ( n_save ) :: a_save = (/  &
      0.13819660112501053195D+00, &
      0.13819660112501053195D+00, &
      0.58541019662496851517D+00, &
      0.13819660112501053195D+00 /)

  real ( kind = rk ) b(n)
  real ( kind = rk ), save, dimension ( n_save ) :: b_save = (/ &
      0.13819660112501053195D+00, &
      0.58541019662496851517D+00, &
      0.13819660112501053195D+00, &
      0.13819660112501053195D+00 /)

  real ( kind = rk ) c(n)
  real ( kind = rk ), save, dimension ( n_save ) :: c_save = (/ &
      0.58541019662496851517D+00, &
      0.13819660112501053195D+00, &
      0.13819660112501053195D+00, &
      0.13819660112501053195D+00 /)

  real ( kind = rk ) d(n)
  real ( kind = rk ), save, dimension ( n_save ) :: d_save = (/ &
      0.13819660112501042093D+00, &
      0.13819660112501042093D+00, &
      0.13819660112501042093D+00, &
      0.58541019662496840414D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
      0.25000000000000000000D+00, &
      0.25000000000000000000D+00, &
      0.25000000000000000000D+00, &
      0.25000000000000000000D+00 /)

  a(1:n) = a_save(1:n)
  b(1:n) = b_save(1:n)
  c(1:n) = c_save(1:n)
  d(1:n) = d_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule03 ( n, a, b, c, d, w )

!*****************************************************************************80
!
!! rule03() returns the rule of precision 3.
!
!  Discussion:
!
!    The data is given for the reference tetrahedron
!    with vertices (0,0,0), (1,0,0), (0,1,0), (0,0,1).
!
!    We suppose we are given a tetrahedron T with vertices A, B, C, D.
!    We call a rule with n points, returning barycentric coordinates
!    a, b, c, d, and weights w.  Then the integral I of f(x,y,z) over T is 
!    approximated by Q as follows:
!
!    (x,y,z) = a(1:n) * A + b(1:n) * B + c(1:n) * C + d(1:n * D
!    Q = area(T) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 April 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Freddie Witherden, Peter Vincent,
!    On the identification of symmetric quadrature rules for finite element methods,
!    Computers and Mathematics with Applications,
!    Volume 69, pages 1232-1241, 2015.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real a(n), b(n), c(n), d(n): the barycentric coordinates of quadrature points.
!
!    real w(n): the quadrature weights, which sum to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 8

  real ( kind = rk ) a(n)
  real ( kind = rk ), save, dimension ( n_save ) :: a_save = (/  &
      0.32816330251638170523D+00, &
      0.32816330251638170523D+00, &
      0.01551009245085493982D+00, &
      0.32816330251638170523D+00, &
      0.10804724989842862115D+00, &
      0.10804724989842862115D+00, &
      0.67585825030471413655D+00, &
      0.10804724989842862115D+00 /)

  real ( kind = rk ) b(n)
  real ( kind = rk ), save, dimension ( n_save ) :: b_save = (/ &
      0.32816330251638170523D+00, &
      0.01551009245085493982D+00, &
      0.32816330251638170523D+00, &
      0.32816330251638170523D+00, &
      0.10804724989842862115D+00, &
      0.67585825030471413655D+00, &
      0.10804724989842862115D+00, &
      0.10804724989842862115D+00 /)

  real ( kind = rk ) c(n)
  real ( kind = rk ), save, dimension ( n_save ) :: c_save = (/ &
      0.01551009245085493982D+00, &
      0.32816330251638170523D+00, &
      0.32816330251638170523D+00, &
      0.32816330251638170523D+00, &
      0.67585825030471413655D+00, &
      0.10804724989842862115D+00, &
      0.10804724989842862115D+00, &
      0.10804724989842862115D+00 /)

  real ( kind = rk ) d(n)
  real ( kind = rk ), save, dimension ( n_save ) :: d_save = (/ &
      0.32816330251638164972D+00, &
      0.32816330251638170523D+00, &
      0.32816330251638170523D+00, &
      0.01551009245085488431D+00, &
      0.10804724989842862115D+00, &
      0.10804724989842862115D+00, &
      0.10804724989842862115D+00, &
      0.67585825030471413655D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
      0.13621784253708735246D+00, &
      0.13621784253708735246D+00, &
      0.13621784253708735246D+00, &
      0.13621784253708735246D+00, &
      0.11378215746291264754D+00, &
      0.11378215746291264754D+00, &
      0.11378215746291264754D+00, &
      0.11378215746291264754D+00 /)

  a(1:n) = a_save(1:n)
  b(1:n) = b_save(1:n)
  c(1:n) = c_save(1:n)
  d(1:n) = d_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule04 ( n, a, b, c, d, w )

!*****************************************************************************80
!
!! rule04() returns the rule of precision 4.
!
!  Discussion:
!
!    The data is given for the reference tetrahedron
!    with vertices (0,0,0), (1,0,0), (0,1,0), (0,0,1).
!
!    We suppose we are given a tetrahedron T with vertices A, B, C, D.
!    We call a rule with n points, returning barycentric coordinates
!    a, b, c, d, and weights w.  Then the integral I of f(x,y,z) over T is 
!    approximated by Q as follows:
!
!    (x,y,z) = a(1:n) * A + b(1:n) * B + c(1:n) * C + d(1:n * D
!    Q = area(T) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 April 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Freddie Witherden, Peter Vincent,
!    On the identification of symmetric quadrature rules for finite element methods,
!    Computers and Mathematics with Applications,
!    Volume 69, pages 1232-1241, 2015.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real a(n), b(n), c(n), d(n): the barycentric coordinates of quadrature points.
!
!    real w(n): the quadrature weights, which sum to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 14

  real ( kind = rk ) a(n)
  real ( kind = rk ), save, dimension ( n_save ) :: a_save = (/  &
      0.31088591926330061410D+00, &
      0.31088591926330061410D+00, &
      0.06734224221009815770D+00, &
      0.31088591926330061410D+00, &
      0.09273525031089124848D+00, &
      0.09273525031089124848D+00, &
      0.72179424906732636558D+00, &
      0.09273525031089124848D+00, &
      0.04550370412564963551D+00, &
      0.45449629587435036449D+00, &
      0.04550370412564963551D+00, &
      0.04550370412564963551D+00, &
      0.45449629587435036449D+00, &
      0.45449629587435036449D+00 /)

  real ( kind = rk ) b(n)
  real ( kind = rk ), save, dimension ( n_save ) :: b_save = (/ &
      0.31088591926330061410D+00, &
      0.06734224221009815770D+00, &
      0.31088591926330061410D+00, &
      0.31088591926330061410D+00, &
      0.09273525031089124848D+00, &
      0.72179424906732636558D+00, &
      0.09273525031089124848D+00, &
      0.09273525031089124848D+00, &
      0.45449629587435036449D+00, &
      0.04550370412564963551D+00, &
      0.04550370412564963551D+00, &
      0.45449629587435036449D+00, &
      0.04550370412564963551D+00, &
      0.45449629587435036449D+00 /)

  real ( kind = rk ) c(n)
  real ( kind = rk ), save, dimension ( n_save ) :: c_save = (/ &
      0.06734224221009815770D+00, &
      0.31088591926330061410D+00, &
      0.31088591926330061410D+00, &
      0.31088591926330061410D+00, &
      0.72179424906732636558D+00, &
      0.09273525031089124848D+00, &
      0.09273525031089124848D+00, &
      0.09273525031089124848D+00, &
      0.45449629587435036449D+00, &
      0.45449629587435036449D+00, &
      0.45449629587435036449D+00, &
      0.04550370412564963551D+00, &
      0.04550370412564963551D+00, &
      0.04550370412564963551D+00 /)

  real ( kind = rk ) d(n)
  real ( kind = rk ), save, dimension ( n_save ) :: d_save = (/ &
      0.31088591926330055859D+00, &
      0.31088591926330050308D+00, &
      0.31088591926330050308D+00, &
      0.06734224221009810218D+00, &
      0.09273525031089113746D+00, &
      0.09273525031089113746D+00, &
      0.09273525031089113746D+00, &
      0.72179424906732625455D+00, &
      0.04550370412564963551D+00, &
      0.04550370412564963551D+00, &
      0.45449629587435036449D+00, &
      0.45449629587435036449D+00, &
      0.45449629587435036449D+00, &
      0.04550370412564963551D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
      0.11268792571801586333D+00, &
      0.11268792571801586333D+00, &
      0.11268792571801586333D+00, &
      0.11268792571801586333D+00, &
      0.07349304311636195575D+00, &
      0.07349304311636195575D+00, &
      0.07349304311636195575D+00, &
      0.07349304311636195575D+00, &
      0.04254602077708146551D+00, &
      0.04254602077708146551D+00, &
      0.04254602077708146551D+00, &
      0.04254602077708146551D+00, &
      0.04254602077708146551D+00, &
      0.04254602077708146551D+00 /)

  a(1:n) = a_save(1:n)
  b(1:n) = b_save(1:n)
  c(1:n) = c_save(1:n)
  d(1:n) = d_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule05 ( n, a, b, c, d, w )

!*****************************************************************************80
!
!! rule05() returns the rule of precision 5.
!
!  Discussion:
!
!    The data is given for the reference tetrahedron
!    with vertices (0,0,0), (1,0,0), (0,1,0), (0,0,1).
!
!    We suppose we are given a tetrahedron T with vertices A, B, C, D.
!    We call a rule with n points, returning barycentric coordinates
!    a, b, c, d, and weights w.  Then the integral I of f(x,y,z) over T is 
!    approximated by Q as follows:
!
!    (x,y,z) = a(1:n) * A + b(1:n) * B + c(1:n) * C + d(1:n * D
!    Q = area(T) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 April 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Freddie Witherden, Peter Vincent,
!    On the identification of symmetric quadrature rules for finite element methods,
!    Computers and Mathematics with Applications,
!    Volume 69, pages 1232-1241, 2015.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real a(n), b(n), c(n), d(n): the barycentric coordinates of quadrature points.
!
!    real w(n): the quadrature weights, which sum to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 14

  real ( kind = rk ) a(n)
  real ( kind = rk ), save, dimension ( n_save ) :: a_save = (/  &
      0.31088591926330061410D+00, &
      0.31088591926330061410D+00, &
      0.06734224221009815770D+00, &
      0.31088591926330061410D+00, &
      0.09273525031089124848D+00, &
      0.09273525031089124848D+00, &
      0.72179424906732636558D+00, &
      0.09273525031089124848D+00, &
      0.04550370412564963551D+00, &
      0.45449629587435036449D+00, &
      0.04550370412564963551D+00, &
      0.04550370412564963551D+00, &
      0.45449629587435036449D+00, &
      0.45449629587435036449D+00 /)

  real ( kind = rk ) b(n)
  real ( kind = rk ), save, dimension ( n_save ) :: b_save = (/ &
      0.31088591926330061410D+00, &
      0.06734224221009815770D+00, &
      0.31088591926330061410D+00, &
      0.31088591926330061410D+00, &
      0.09273525031089124848D+00, &
      0.72179424906732636558D+00, &
      0.09273525031089124848D+00, &
      0.09273525031089124848D+00, &
      0.45449629587435036449D+00, &
      0.04550370412564963551D+00, &
      0.04550370412564963551D+00, &
      0.45449629587435036449D+00, &
      0.04550370412564963551D+00, &
      0.45449629587435036449D+00 /)

  real ( kind = rk ) c(n)
  real ( kind = rk ), save, dimension ( n_save ) :: c_save = (/ &
      0.06734224221009815770D+00, &
      0.31088591926330061410D+00, &
      0.31088591926330061410D+00, &
      0.31088591926330061410D+00, &
      0.72179424906732636558D+00, &
      0.09273525031089124848D+00, &
      0.09273525031089124848D+00, &
      0.09273525031089124848D+00, &
      0.45449629587435036449D+00, &
      0.45449629587435036449D+00, &
      0.45449629587435036449D+00, &
      0.04550370412564963551D+00, &
      0.04550370412564963551D+00, &
      0.04550370412564963551D+00 /)

  real ( kind = rk ) d(n)
  real ( kind = rk ), save, dimension ( n_save ) :: d_save = (/ &
      0.31088591926330055859D+00, &
      0.31088591926330050308D+00, &
      0.31088591926330050308D+00, &
      0.06734224221009810218D+00, &
      0.09273525031089113746D+00, &
      0.09273525031089113746D+00, &
      0.09273525031089113746D+00, &
      0.72179424906732625455D+00, &
      0.04550370412564963551D+00, &
      0.04550370412564963551D+00, &
      0.45449629587435036449D+00, &
      0.45449629587435036449D+00, &
      0.45449629587435036449D+00, &
      0.04550370412564963551D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
      0.11268792571801586333D+00, &
      0.11268792571801586333D+00, &
      0.11268792571801586333D+00, &
      0.11268792571801586333D+00, &
      0.07349304311636195575D+00, &
      0.07349304311636195575D+00, &
      0.07349304311636195575D+00, &
      0.07349304311636195575D+00, &
      0.04254602077708146551D+00, &
      0.04254602077708146551D+00, &
      0.04254602077708146551D+00, &
      0.04254602077708146551D+00, &
      0.04254602077708146551D+00, &
      0.04254602077708146551D+00 /)

  a(1:n) = a_save(1:n)
  b(1:n) = b_save(1:n)
  c(1:n) = c_save(1:n)
  d(1:n) = d_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule06 ( n, a, b, c, d, w )

!*****************************************************************************80
!
!! rule06() returns the rule of precision 6.
!
!  Discussion:
!
!    The data is given for the reference tetrahedron
!    with vertices (0,0,0), (1,0,0), (0,1,0), (0,0,1).
!
!    We suppose we are given a tetrahedron T with vertices A, B, C, D.
!    We call a rule with n points, returning barycentric coordinates
!    a, b, c, d, and weights w.  Then the integral I of f(x,y,z) over T is 
!    approximated by Q as follows:
!
!    (x,y,z) = a(1:n) * A + b(1:n) * B + c(1:n) * C + d(1:n * D
!    Q = area(T) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 April 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Freddie Witherden, Peter Vincent,
!    On the identification of symmetric quadrature rules for finite element methods,
!    Computers and Mathematics with Applications,
!    Volume 69, pages 1232-1241, 2015.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real a(n), b(n), c(n), d(n): the barycentric coordinates of quadrature points.
!
!    real w(n): the quadrature weights, which sum to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 24

  real ( kind = rk ) a(n)
  real ( kind = rk ), save, dimension ( n_save ) :: a_save = (/  &
      0.04067395853461136523D+00, &
      0.04067395853461136523D+00, &
      0.87797812439616595981D+00, &
      0.04067395853461136523D+00, &
      0.32233789014227554048D+00, &
      0.32233789014227554048D+00, &
      0.03298632957317348957D+00, &
      0.32233789014227554048D+00, &
      0.21460287125915200601D+00, &
      0.21460287125915200601D+00, &
      0.35619138622254387094D+00, &
      0.21460287125915200601D+00, &
      0.60300566479164918743D+00, &
      0.60300566479164918743D+00, &
      0.06366100187501749774D+00, &
      0.26967233145831581709D+00, &
      0.06366100187501749774D+00, &
      0.06366100187501749774D+00, &
      0.26967233145831581709D+00, &
      0.06366100187501749774D+00, &
      0.06366100187501749774D+00, &
      0.06366100187501749774D+00, &
      0.26967233145831581709D+00, &
      0.60300566479164918743D+00 /)

  real ( kind = rk ) b(n)
  real ( kind = rk ), save, dimension ( n_save ) :: b_save = (/ &
      0.04067395853461136523D+00, &
      0.87797812439616595981D+00, &
      0.04067395853461136523D+00, &
      0.04067395853461136523D+00, &
      0.32233789014227554048D+00, &
      0.03298632957317348957D+00, &
      0.32233789014227554048D+00, &
      0.32233789014227554048D+00, &
      0.21460287125915200601D+00, &
      0.35619138622254387094D+00, &
      0.21460287125915200601D+00, &
      0.21460287125915200601D+00, &
      0.06366100187501749774D+00, &
      0.06366100187501749774D+00, &
      0.06366100187501749774D+00, &
      0.60300566479164918743D+00, &
      0.26967233145831581709D+00, &
      0.60300566479164918743D+00, &
      0.06366100187501749774D+00, &
      0.26967233145831581709D+00, &
      0.06366100187501749774D+00, &
      0.60300566479164918743D+00, &
      0.06366100187501749774D+00, &
      0.26967233145831581709D+00 /)

  real ( kind = rk ) c(n)
  real ( kind = rk ), save, dimension ( n_save ) :: c_save = (/ &
      0.87797812439616595981D+00, &
      0.04067395853461136523D+00, &
      0.04067395853461136523D+00, &
      0.04067395853461136523D+00, &
      0.03298632957317348957D+00, &
      0.32233789014227554048D+00, &
      0.32233789014227554048D+00, &
      0.32233789014227554048D+00, &
      0.35619138622254387094D+00, &
      0.21460287125915200601D+00, &
      0.21460287125915200601D+00, &
      0.21460287125915200601D+00, &
      0.26967233145831581709D+00, &
      0.06366100187501749774D+00, &
      0.60300566479164918743D+00, &
      0.06366100187501749774D+00, &
      0.60300566479164918743D+00, &
      0.06366100187501749774D+00, &
      0.60300566479164918743D+00, &
      0.06366100187501749774D+00, &
      0.26967233145831581709D+00, &
      0.26967233145831581709D+00, &
      0.06366100187501749774D+00, &
      0.06366100187501749774D+00 /)

  real ( kind = rk ) d(n)
  real ( kind = rk ), save, dimension ( n_save ) :: d_save = (/ &
      0.04067395853461119870D+00, &
      0.04067395853461125421D+00, &
      0.04067395853461130972D+00, &
      0.87797812439616573776D+00, &
      0.32233789014227542946D+00, &
      0.32233789014227542946D+00, &
      0.32233789014227542946D+00, &
      0.03298632957317337855D+00, &
      0.21460287125915211703D+00, &
      0.21460287125915211703D+00, &
      0.21460287125915211703D+00, &
      0.35619138622254398197D+00, &
      0.06366100187501749774D+00, &
      0.26967233145831581709D+00, &
      0.26967233145831592811D+00, &
      0.06366100187501749774D+00, &
      0.06366100187501755325D+00, &
      0.26967233145831587260D+00, &
      0.06366100187501755325D+00, &
      0.60300566479164929845D+00, &
      0.60300566479164929845D+00, &
      0.06366100187501755325D+00, &
      0.60300566479164929845D+00, &
      0.06366100187501749774D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
      0.01007721105532064301D+00, &
      0.01007721105532064301D+00, &
      0.01007721105532064301D+00, &
      0.01007721105532064301D+00, &
      0.05535718154365472377D+00, &
      0.05535718154365472377D+00, &
      0.05535718154365472377D+00, &
      0.05535718154365472377D+00, &
      0.03992275025816749423D+00, &
      0.03992275025816749423D+00, &
      0.03992275025816749423D+00, &
      0.03992275025816749423D+00, &
      0.04821428571428570953D+00, &
      0.04821428571428570953D+00, &
      0.04821428571428570953D+00, &
      0.04821428571428570953D+00, &
      0.04821428571428570953D+00, &
      0.04821428571428570953D+00, &
      0.04821428571428570953D+00, &
      0.04821428571428570953D+00, &
      0.04821428571428570953D+00, &
      0.04821428571428570953D+00, &
      0.04821428571428570953D+00, &
      0.04821428571428570953D+00 /)

  a(1:n) = a_save(1:n)
  b(1:n) = b_save(1:n)
  c(1:n) = c_save(1:n)
  d(1:n) = d_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule07 ( n, a, b, c, d, w )

!*****************************************************************************80
!
!! rule07() returns the rule of precision 7.
!
!  Discussion:
!
!    The data is given for the reference tetrahedron
!    with vertices (0,0,0), (1,0,0), (0,1,0), (0,0,1).
!
!    We suppose we are given a tetrahedron T with vertices A, B, C, D.
!    We call a rule with n points, returning barycentric coordinates
!    a, b, c, d, and weights w.  Then the integral I of f(x,y,z) over T is 
!    approximated by Q as follows:
!
!    (x,y,z) = a(1:n) * A + b(1:n) * B + c(1:n) * C + d(1:n * D
!    Q = area(T) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 April 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Freddie Witherden, Peter Vincent,
!    On the identification of symmetric quadrature rules for finite element methods,
!    Computers and Mathematics with Applications,
!    Volume 69, pages 1232-1241, 2015.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real a(n), b(n), c(n), d(n): the barycentric coordinates of quadrature points.
!
!    real w(n): the quadrature weights, which sum to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 35

  real ( kind = rk ) a(n)
  real ( kind = rk ), save, dimension ( n_save ) :: a_save = (/  &
      0.25000000000000000000D+00, &
      0.31570114977820279423D+00, &
      0.31570114977820279423D+00, &
      0.05289655066539161732D+00, &
      0.31570114977820279423D+00, &
      0.05048982259839635001D+00, &
      0.44951017740160364999D+00, &
      0.05048982259839635001D+00, &
      0.05048982259839635001D+00, &
      0.44951017740160364999D+00, &
      0.44951017740160364999D+00, &
      0.57517163758699996201D+00, &
      0.57517163758699996201D+00, &
      0.18883383102600104220D+00, &
      0.04716070036099789808D+00, &
      0.18883383102600104220D+00, &
      0.18883383102600104220D+00, &
      0.04716070036099789808D+00, &
      0.18883383102600104220D+00, &
      0.18883383102600104220D+00, &
      0.18883383102600104220D+00, &
      0.04716070036099789808D+00, &
      0.57517163758699996201D+00, &
      0.81083024109854862083D+00, &
      0.81083024109854862083D+00, &
      0.02126547254148325461D+00, &
      0.14663881381848492547D+00, &
      0.02126547254148325461D+00, &
      0.02126547254148325461D+00, &
      0.14663881381848492547D+00, &
      0.02126547254148325461D+00, &
      0.02126547254148325461D+00, &
      0.02126547254148325461D+00, &
      0.14663881381848492547D+00, &
      0.81083024109854862083D+00 /)

  real ( kind = rk ) b(n)
  real ( kind = rk ), save, dimension ( n_save ) :: b_save = (/ &
      0.25000000000000000000D+00, &
      0.31570114977820279423D+00, &
      0.05289655066539161732D+00, &
      0.31570114977820279423D+00, &
      0.31570114977820279423D+00, &
      0.44951017740160364999D+00, &
      0.05048982259839635001D+00, &
      0.05048982259839635001D+00, &
      0.44951017740160364999D+00, &
      0.05048982259839635001D+00, &
      0.44951017740160364999D+00, &
      0.18883383102600104220D+00, &
      0.18883383102600104220D+00, &
      0.18883383102600104220D+00, &
      0.57517163758699996201D+00, &
      0.04716070036099789808D+00, &
      0.57517163758699996201D+00, &
      0.18883383102600104220D+00, &
      0.04716070036099789808D+00, &
      0.18883383102600104220D+00, &
      0.57517163758699996201D+00, &
      0.18883383102600104220D+00, &
      0.04716070036099789808D+00, &
      0.02126547254148325461D+00, &
      0.02126547254148325461D+00, &
      0.02126547254148325461D+00, &
      0.81083024109854862083D+00, &
      0.14663881381848492547D+00, &
      0.81083024109854862083D+00, &
      0.02126547254148325461D+00, &
      0.14663881381848492547D+00, &
      0.02126547254148325461D+00, &
      0.81083024109854862083D+00, &
      0.02126547254148325461D+00, &
      0.14663881381848492547D+00 /)

  real ( kind = rk ) c(n)
  real ( kind = rk ), save, dimension ( n_save ) :: c_save = (/ &
      0.25000000000000000000D+00, &
      0.05289655066539161732D+00, &
      0.31570114977820279423D+00, &
      0.31570114977820279423D+00, &
      0.31570114977820279423D+00, &
      0.44951017740160364999D+00, &
      0.44951017740160364999D+00, &
      0.44951017740160364999D+00, &
      0.05048982259839635001D+00, &
      0.05048982259839635001D+00, &
      0.05048982259839635001D+00, &
      0.04716070036099789808D+00, &
      0.18883383102600104220D+00, &
      0.57517163758699996201D+00, &
      0.18883383102600104220D+00, &
      0.57517163758699996201D+00, &
      0.18883383102600104220D+00, &
      0.57517163758699996201D+00, &
      0.18883383102600104220D+00, &
      0.04716070036099789808D+00, &
      0.04716070036099789808D+00, &
      0.18883383102600104220D+00, &
      0.18883383102600104220D+00, &
      0.14663881381848492547D+00, &
      0.02126547254148325461D+00, &
      0.81083024109854862083D+00, &
      0.02126547254148325461D+00, &
      0.81083024109854862083D+00, &
      0.02126547254148325461D+00, &
      0.81083024109854862083D+00, &
      0.02126547254148325461D+00, &
      0.14663881381848492547D+00, &
      0.14663881381848492547D+00, &
      0.02126547254148325461D+00, &
      0.02126547254148325461D+00 /)

  real ( kind = rk ) d(n)
  real ( kind = rk ), save, dimension ( n_save ) :: d_save = (/ &
      0.25000000000000000000D+00, &
      0.31570114977820284974D+00, &
      0.31570114977820290525D+00, &
      0.31570114977820290525D+00, &
      0.05289655066539167283D+00, &
      0.05048982259839635001D+00, &
      0.05048982259839629450D+00, &
      0.44951017740160376102D+00, &
      0.44951017740160364999D+00, &
      0.44951017740160359448D+00, &
      0.05048982259839629450D+00, &
      0.18883383102600109771D+00, &
      0.04716070036099795360D+00, &
      0.04716070036099795360D+00, &
      0.18883383102600104220D+00, &
      0.18883383102600115322D+00, &
      0.04716070036099795360D+00, &
      0.18883383102600104220D+00, &
      0.57517163758700007303D+00, &
      0.57517163758699996201D+00, &
      0.18883383102600109771D+00, &
      0.57517163758699996201D+00, &
      0.18883383102600109771D+00, &
      0.02126547254148319910D+00, &
      0.14663881381848486996D+00, &
      0.14663881381848486996D+00, &
      0.02126547254148325461D+00, &
      0.02126547254148314359D+00, &
      0.14663881381848486996D+00, &
      0.02126547254148325461D+00, &
      0.81083024109854850980D+00, &
      0.81083024109854862083D+00, &
      0.02126547254148319910D+00, &
      0.81083024109854862083D+00, &
      0.02126547254148319910D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
      0.09548528946413084584D+00, &
      0.04232958120996702794D+00, &
      0.04232958120996702794D+00, &
      0.04232958120996702794D+00, &
      0.04232958120996702794D+00, &
      0.03189692783285758004D+00, &
      0.03189692783285758004D+00, &
      0.03189692783285758004D+00, &
      0.03189692783285758004D+00, &
      0.03189692783285758004D+00, &
      0.03189692783285758004D+00, &
      0.03720713072833461976D+00, &
      0.03720713072833461976D+00, &
      0.03720713072833461976D+00, &
      0.03720713072833461976D+00, &
      0.03720713072833461976D+00, &
      0.03720713072833461976D+00, &
      0.03720713072833461976D+00, &
      0.03720713072833461976D+00, &
      0.03720713072833461976D+00, &
      0.03720713072833461976D+00, &
      0.03720713072833461976D+00, &
      0.03720713072833461976D+00, &
      0.00811077082990334201D+00, &
      0.00811077082990334201D+00, &
      0.00811077082990334201D+00, &
      0.00811077082990334201D+00, &
      0.00811077082990334201D+00, &
      0.00811077082990334201D+00, &
      0.00811077082990334201D+00, &
      0.00811077082990334201D+00, &
      0.00811077082990334201D+00, &
      0.00811077082990334201D+00, &
      0.00811077082990334201D+00, &
      0.00811077082990334201D+00 /)

  a(1:n) = a_save(1:n)
  b(1:n) = b_save(1:n)
  c(1:n) = c_save(1:n)
  d(1:n) = d_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule08 ( n, a, b, c, d, w )

!*****************************************************************************80
!
!! rule08() returns the rule of precision 8.
!
!  Discussion:
!
!    The data is given for the reference tetrahedron
!    with vertices (0,0,0), (1,0,0), (0,1,0), (0,0,1).
!
!    We suppose we are given a tetrahedron T with vertices A, B, C, D.
!    We call a rule with n points, returning barycentric coordinates
!    a, b, c, d, and weights w.  Then the integral I of f(x,y,z) over T is 
!    approximated by Q as follows:
!
!    (x,y,z) = a(1:n) * A + b(1:n) * B + c(1:n) * C + d(1:n * D
!    Q = area(T) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 April 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Freddie Witherden, Peter Vincent,
!    On the identification of symmetric quadrature rules for finite element methods,
!    Computers and Mathematics with Applications,
!    Volume 69, pages 1232-1241, 2015.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real a(n), b(n), c(n), d(n): the barycentric coordinates of quadrature points.
!
!    real w(n): the quadrature weights, which sum to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 46

  real ( kind = rk ) a(n)
  real ( kind = rk ), save, dimension ( n_save ) :: a_save = (/  &
      0.10795272496221086644D+00, &
      0.10795272496221086644D+00, &
      0.67614182511336751169D+00, &
      0.10795272496221086644D+00, &
      0.18510948778258656811D+00, &
      0.18510948778258656811D+00, &
      0.44467153665224029568D+00, &
      0.18510948778258656811D+00, &
      0.04231654368476728267D+00, &
      0.04231654368476728267D+00, &
      0.87305036894569809647D+00, &
      0.04231654368476728267D+00, &
      0.31418170912403897699D+00, &
      0.31418170912403897699D+00, &
      0.05745487262788301353D+00, &
      0.31418170912403897699D+00, &
      0.43559132858383020626D+00, &
      0.06440867141616979374D+00, &
      0.43559132858383020626D+00, &
      0.43559132858383020626D+00, &
      0.06440867141616979374D+00, &
      0.06440867141616979374D+00, &
      0.71746406342630830721D+00, &
      0.71746406342630830721D+00, &
      0.02143393012713057377D+00, &
      0.23966807631943054524D+00, &
      0.02143393012713057377D+00, &
      0.02143393012713057377D+00, &
      0.23966807631943054524D+00, &
      0.02143393012713057377D+00, &
      0.02143393012713057377D+00, &
      0.02143393012713057377D+00, &
      0.23966807631943054524D+00, &
      0.71746406342630830721D+00, &
      0.58379737830214439853D+00, &
      0.58379737830214439853D+00, &
      0.20413933387602911651D+00, &
      0.00792395394579736845D+00, &
      0.20413933387602911651D+00, &
      0.20413933387602911651D+00, &
      0.00792395394579736845D+00, &
      0.20413933387602911651D+00, &
      0.20413933387602911651D+00, &
      0.20413933387602911651D+00, &
      0.00792395394579736845D+00, &
      0.58379737830214439853D+00 /)

  real ( kind = rk ) b(n)
  real ( kind = rk ), save, dimension ( n_save ) :: b_save = (/ &
      0.10795272496221086644D+00, &
      0.67614182511336751169D+00, &
      0.10795272496221086644D+00, &
      0.10795272496221086644D+00, &
      0.18510948778258656811D+00, &
      0.44467153665224029568D+00, &
      0.18510948778258656811D+00, &
      0.18510948778258656811D+00, &
      0.04231654368476728267D+00, &
      0.87305036894569809647D+00, &
      0.04231654368476728267D+00, &
      0.04231654368476728267D+00, &
      0.31418170912403897699D+00, &
      0.05745487262788301353D+00, &
      0.31418170912403897699D+00, &
      0.31418170912403897699D+00, &
      0.06440867141616979374D+00, &
      0.43559132858383020626D+00, &
      0.43559132858383020626D+00, &
      0.06440867141616979374D+00, &
      0.43559132858383020626D+00, &
      0.06440867141616979374D+00, &
      0.02143393012713057377D+00, &
      0.02143393012713057377D+00, &
      0.02143393012713057377D+00, &
      0.71746406342630830721D+00, &
      0.23966807631943054524D+00, &
      0.71746406342630830721D+00, &
      0.02143393012713057377D+00, &
      0.23966807631943054524D+00, &
      0.02143393012713057377D+00, &
      0.71746406342630830721D+00, &
      0.02143393012713057377D+00, &
      0.23966807631943054524D+00, &
      0.20413933387602911651D+00, &
      0.20413933387602911651D+00, &
      0.20413933387602911651D+00, &
      0.58379737830214439853D+00, &
      0.00792395394579736845D+00, &
      0.58379737830214439853D+00, &
      0.20413933387602911651D+00, &
      0.00792395394579736845D+00, &
      0.20413933387602911651D+00, &
      0.58379737830214439853D+00, &
      0.20413933387602911651D+00, &
      0.00792395394579736845D+00 /)

  real ( kind = rk ) c(n)
  real ( kind = rk ), save, dimension ( n_save ) :: c_save = (/ &
      0.67614182511336751169D+00, &
      0.10795272496221086644D+00, &
      0.10795272496221086644D+00, &
      0.10795272496221086644D+00, &
      0.44467153665224029568D+00, &
      0.18510948778258656811D+00, &
      0.18510948778258656811D+00, &
      0.18510948778258656811D+00, &
      0.87305036894569809647D+00, &
      0.04231654368476728267D+00, &
      0.04231654368476728267D+00, &
      0.04231654368476728267D+00, &
      0.05745487262788301353D+00, &
      0.31418170912403897699D+00, &
      0.31418170912403897699D+00, &
      0.31418170912403897699D+00, &
      0.06440867141616979374D+00, &
      0.06440867141616979374D+00, &
      0.06440867141616979374D+00, &
      0.43559132858383020626D+00, &
      0.43559132858383020626D+00, &
      0.43559132858383020626D+00, &
      0.23966807631943054524D+00, &
      0.02143393012713057377D+00, &
      0.71746406342630830721D+00, &
      0.02143393012713057377D+00, &
      0.71746406342630830721D+00, &
      0.02143393012713057377D+00, &
      0.71746406342630830721D+00, &
      0.02143393012713057377D+00, &
      0.23966807631943054524D+00, &
      0.23966807631943054524D+00, &
      0.02143393012713057377D+00, &
      0.02143393012713057377D+00, &
      0.00792395394579736845D+00, &
      0.20413933387602911651D+00, &
      0.58379737830214439853D+00, &
      0.20413933387602911651D+00, &
      0.58379737830214439853D+00, &
      0.20413933387602911651D+00, &
      0.58379737830214439853D+00, &
      0.20413933387602911651D+00, &
      0.00792395394579736845D+00, &
      0.00792395394579736845D+00, &
      0.20413933387602911651D+00, &
      0.20413933387602911651D+00 /)

  real ( kind = rk ) d(n)
  real ( kind = rk ), save, dimension ( n_save ) :: d_save = (/ &
      0.10795272496221075542D+00, &
      0.10795272496221075542D+00, &
      0.10795272496221075542D+00, &
      0.67614182511336740067D+00, &
      0.18510948778258667913D+00, &
      0.18510948778258662362D+00, &
      0.18510948778258662362D+00, &
      0.44467153665224040671D+00, &
      0.04231654368476744921D+00, &
      0.04231654368476739370D+00, &
      0.04231654368476733818D+00, &
      0.87305036894569831851D+00, &
      0.31418170912403903250D+00, &
      0.31418170912403897699D+00, &
      0.31418170912403908801D+00, &
      0.05745487262788306904D+00, &
      0.43559132858383020626D+00, &
      0.43559132858383020626D+00, &
      0.06440867141616979374D+00, &
      0.06440867141616979374D+00, &
      0.06440867141616979374D+00, &
      0.43559132858383020626D+00, &
      0.02143393012713057377D+00, &
      0.23966807631943054524D+00, &
      0.23966807631943054524D+00, &
      0.02143393012713057377D+00, &
      0.02143393012713057377D+00, &
      0.23966807631943054524D+00, &
      0.02143393012713057377D+00, &
      0.71746406342630830721D+00, &
      0.71746406342630830721D+00, &
      0.02143393012713057377D+00, &
      0.71746406342630830721D+00, &
      0.02143393012713057377D+00, &
      0.20413933387602911651D+00, &
      0.00792395394579736845D+00, &
      0.00792395394579736845D+00, &
      0.20413933387602911651D+00, &
      0.20413933387602911651D+00, &
      0.00792395394579736845D+00, &
      0.20413933387602911651D+00, &
      0.58379737830214439853D+00, &
      0.58379737830214439853D+00, &
      0.20413933387602911651D+00, &
      0.58379737830214439853D+00, &
      0.20413933387602911651D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
      0.02642665090840883024D+00, &
      0.02642665090840883024D+00, &
      0.02642665090840883024D+00, &
      0.02642665090840883024D+00, &
      0.05203174756373853127D+00, &
      0.05203174756373853127D+00, &
      0.05203174756373853127D+00, &
      0.05203174756373853127D+00, &
      0.00752525615354019892D+00, &
      0.00752525615354019892D+00, &
      0.00752525615354019892D+00, &
      0.00752525615354019892D+00, &
      0.04176378285693489734D+00, &
      0.04176378285693489734D+00, &
      0.04176378285693489734D+00, &
      0.04176378285693489734D+00, &
      0.03628093026130882470D+00, &
      0.03628093026130882470D+00, &
      0.03628093026130882470D+00, &
      0.03628093026130882470D+00, &
      0.03628093026130882470D+00, &
      0.03628093026130882470D+00, &
      0.00715690289084443265D+00, &
      0.00715690289084443265D+00, &
      0.00715690289084443265D+00, &
      0.00715690289084443265D+00, &
      0.00715690289084443265D+00, &
      0.00715690289084443265D+00, &
      0.00715690289084443265D+00, &
      0.00715690289084443265D+00, &
      0.00715690289084443265D+00, &
      0.00715690289084443265D+00, &
      0.00715690289084443265D+00, &
      0.00715690289084443265D+00, &
      0.01545348615096033516D+00, &
      0.01545348615096033516D+00, &
      0.01545348615096033516D+00, &
      0.01545348615096033516D+00, &
      0.01545348615096033516D+00, &
      0.01545348615096033516D+00, &
      0.01545348615096033516D+00, &
      0.01545348615096033516D+00, &
      0.01545348615096033516D+00, &
      0.01545348615096033516D+00, &
      0.01545348615096033516D+00, &
      0.01545348615096033516D+00 /)

  a(1:n) = a_save(1:n)
  b(1:n) = b_save(1:n)
  c(1:n) = c_save(1:n)
  d(1:n) = d_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule09 ( n, a, b, c, d, w )

!*****************************************************************************80
!
!! rule09() returns the rule of precision 9.
!
!  Discussion:
!
!    The data is given for the reference tetrahedron
!    with vertices (0,0,0), (1,0,0), (0,1,0), (0,0,1).
!
!    We suppose we are given a tetrahedron T with vertices A, B, C, D.
!    We call a rule with n points, returning barycentric coordinates
!    a, b, c, d, and weights w.  Then the integral I of f(x,y,z) over T is 
!    approximated by Q as follows:
!
!    (x,y,z) = a(1:n) * A + b(1:n) * B + c(1:n) * C + d(1:n * D
!    Q = area(T) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 April 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Freddie Witherden, Peter Vincent,
!    On the identification of symmetric quadrature rules for finite element methods,
!    Computers and Mathematics with Applications,
!    Volume 69, pages 1232-1241, 2015.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real a(n), b(n), c(n), d(n): the barycentric coordinates of quadrature points.
!
!    real w(n): the quadrature weights, which sum to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 59

  real ( kind = rk ) a(n)
  real ( kind = rk ), save, dimension ( n_save ) :: a_save = (/  &
      0.25000000000000000000D+00, &
      0.00000000061981697552D+00, &
      0.00000000061981697552D+00, &
      0.99999999814054896241D+00, &
      0.00000000061981697552D+00, &
      0.16077453539526159743D+00, &
      0.16077453539526159743D+00, &
      0.51767639381421526323D+00, &
      0.16077453539526159743D+00, &
      0.32227652182142096926D+00, &
      0.32227652182142096926D+00, &
      0.03317043453573709222D+00, &
      0.32227652182142096926D+00, &
      0.04510891834541358447D+00, &
      0.04510891834541358447D+00, &
      0.86467324496375930210D+00, &
      0.04510891834541358447D+00, &
      0.11229654600437605216D+00, &
      0.38770345399562394784D+00, &
      0.11229654600437605216D+00, &
      0.11229654600437605216D+00, &
      0.38770345399562394784D+00, &
      0.38770345399562394784D+00, &
      0.00255457923304130974D+00, &
      0.00255457923304130974D+00, &
      0.45887144875245927667D+00, &
      0.07970252326204013693D+00, &
      0.45887144875245927667D+00, &
      0.45887144875245927667D+00, &
      0.07970252326204013693D+00, &
      0.45887144875245927667D+00, &
      0.45887144875245927667D+00, &
      0.45887144875245927667D+00, &
      0.07970252326204013693D+00, &
      0.00255457923304130974D+00, &
      0.71835032644207452712D+00, &
      0.71835032644207452712D+00, &
      0.03377587068533860482D+00, &
      0.21409793218724831876D+00, &
      0.03377587068533860482D+00, &
      0.03377587068533860482D+00, &
      0.21409793218724831876D+00, &
      0.03377587068533860482D+00, &
      0.03377587068533860482D+00, &
      0.03377587068533860482D+00, &
      0.21409793218724831876D+00, &
      0.71835032644207452712D+00, &
      0.03441591057817527943D+00, &
      0.03441591057817527943D+00, &
      0.18364136980992790127D+00, &
      0.59830134980196891803D+00, &
      0.18364136980992790127D+00, &
      0.18364136980992790127D+00, &
      0.59830134980196891803D+00, &
      0.18364136980992790127D+00, &
      0.18364136980992790127D+00, &
      0.18364136980992790127D+00, &
      0.59830134980196891803D+00, &
      0.03441591057817527943D+00 /)

  real ( kind = rk ) b(n)
  real ( kind = rk ), save, dimension ( n_save ) :: b_save = (/ &
      0.25000000000000000000D+00, &
      0.00000000061981697552D+00, &
      0.99999999814054896241D+00, &
      0.00000000061981697552D+00, &
      0.00000000061981697552D+00, &
      0.16077453539526159743D+00, &
      0.51767639381421526323D+00, &
      0.16077453539526159743D+00, &
      0.16077453539526159743D+00, &
      0.32227652182142096926D+00, &
      0.03317043453573709222D+00, &
      0.32227652182142096926D+00, &
      0.32227652182142096926D+00, &
      0.04510891834541358447D+00, &
      0.86467324496375930210D+00, &
      0.04510891834541358447D+00, &
      0.04510891834541358447D+00, &
      0.38770345399562394784D+00, &
      0.11229654600437605216D+00, &
      0.11229654600437605216D+00, &
      0.38770345399562394784D+00, &
      0.11229654600437605216D+00, &
      0.38770345399562394784D+00, &
      0.45887144875245927667D+00, &
      0.45887144875245927667D+00, &
      0.45887144875245927667D+00, &
      0.00255457923304130974D+00, &
      0.07970252326204013693D+00, &
      0.00255457923304130974D+00, &
      0.45887144875245927667D+00, &
      0.07970252326204013693D+00, &
      0.45887144875245927667D+00, &
      0.00255457923304130974D+00, &
      0.45887144875245927667D+00, &
      0.07970252326204013693D+00, &
      0.03377587068533860482D+00, &
      0.03377587068533860482D+00, &
      0.03377587068533860482D+00, &
      0.71835032644207452712D+00, &
      0.21409793218724831876D+00, &
      0.71835032644207452712D+00, &
      0.03377587068533860482D+00, &
      0.21409793218724831876D+00, &
      0.03377587068533860482D+00, &
      0.71835032644207452712D+00, &
      0.03377587068533860482D+00, &
      0.21409793218724831876D+00, &
      0.18364136980992790127D+00, &
      0.18364136980992790127D+00, &
      0.18364136980992790127D+00, &
      0.03441591057817527943D+00, &
      0.59830134980196891803D+00, &
      0.03441591057817527943D+00, &
      0.18364136980992790127D+00, &
      0.59830134980196891803D+00, &
      0.18364136980992790127D+00, &
      0.03441591057817527943D+00, &
      0.18364136980992790127D+00, &
      0.59830134980196891803D+00 /)

  real ( kind = rk ) c(n)
  real ( kind = rk ), save, dimension ( n_save ) :: c_save = (/ &
      0.25000000000000000000D+00, &
      0.99999999814054896241D+00, &
      0.00000000061981697552D+00, &
      0.00000000061981697552D+00, &
      0.00000000061981697552D+00, &
      0.51767639381421526323D+00, &
      0.16077453539526159743D+00, &
      0.16077453539526159743D+00, &
      0.16077453539526159743D+00, &
      0.03317043453573709222D+00, &
      0.32227652182142096926D+00, &
      0.32227652182142096926D+00, &
      0.32227652182142096926D+00, &
      0.86467324496375930210D+00, &
      0.04510891834541358447D+00, &
      0.04510891834541358447D+00, &
      0.04510891834541358447D+00, &
      0.38770345399562394784D+00, &
      0.38770345399562394784D+00, &
      0.38770345399562394784D+00, &
      0.11229654600437605216D+00, &
      0.11229654600437605216D+00, &
      0.11229654600437605216D+00, &
      0.07970252326204013693D+00, &
      0.45887144875245927667D+00, &
      0.00255457923304130974D+00, &
      0.45887144875245927667D+00, &
      0.00255457923304130974D+00, &
      0.45887144875245927667D+00, &
      0.00255457923304130974D+00, &
      0.45887144875245927667D+00, &
      0.07970252326204013693D+00, &
      0.07970252326204013693D+00, &
      0.45887144875245927667D+00, &
      0.45887144875245927667D+00, &
      0.21409793218724831876D+00, &
      0.03377587068533860482D+00, &
      0.71835032644207452712D+00, &
      0.03377587068533860482D+00, &
      0.71835032644207452712D+00, &
      0.03377587068533860482D+00, &
      0.71835032644207452712D+00, &
      0.03377587068533860482D+00, &
      0.21409793218724831876D+00, &
      0.21409793218724831876D+00, &
      0.03377587068533860482D+00, &
      0.03377587068533860482D+00, &
      0.59830134980196891803D+00, &
      0.18364136980992790127D+00, &
      0.03441591057817527943D+00, &
      0.18364136980992790127D+00, &
      0.03441591057817527943D+00, &
      0.18364136980992790127D+00, &
      0.03441591057817527943D+00, &
      0.18364136980992790127D+00, &
      0.59830134980196891803D+00, &
      0.59830134980196891803D+00, &
      0.18364136980992790127D+00, &
      0.18364136980992790127D+00 /)

  real ( kind = rk ) d(n)
  real ( kind = rk ), save, dimension ( n_save ) :: d_save = (/ &
      0.25000000000000000000D+00, &
      0.00000000061981708654D+00, &
      0.00000000061981708654D+00, &
      0.00000000061981708654D+00, &
      0.99999999814054907343D+00, &
      0.16077453539526143089D+00, &
      0.16077453539526148640D+00, &
      0.16077453539526154191D+00, &
      0.51767639381421504119D+00, &
      0.32227652182142096926D+00, &
      0.32227652182142096926D+00, &
      0.32227652182142096926D+00, &
      0.03317043453573709222D+00, &
      0.04510891834541341794D+00, &
      0.04510891834541347345D+00, &
      0.04510891834541352896D+00, &
      0.86467324496375908005D+00, &
      0.11229654600437599665D+00, &
      0.11229654600437605216D+00, &
      0.38770345399562383681D+00, &
      0.38770345399562389233D+00, &
      0.38770345399562394784D+00, &
      0.11229654600437610767D+00, &
      0.45887144875245938769D+00, &
      0.07970252326204024795D+00, &
      0.07970252326204019244D+00, &
      0.45887144875245938769D+00, &
      0.45887144875245933218D+00, &
      0.07970252326204024795D+00, &
      0.45887144875245933218D+00, &
      0.00255457923304136525D+00, &
      0.00255457923304136525D+00, &
      0.45887144875245938769D+00, &
      0.00255457923304136525D+00, &
      0.45887144875245938769D+00, &
      0.03377587068533854930D+00, &
      0.21409793218724826325D+00, &
      0.21409793218724826325D+00, &
      0.03377587068533860482D+00, &
      0.03377587068533860482D+00, &
      0.21409793218724826325D+00, &
      0.03377587068533860482D+00, &
      0.71835032644207452712D+00, &
      0.71835032644207452712D+00, &
      0.03377587068533854930D+00, &
      0.71835032644207452712D+00, &
      0.03377587068533854930D+00, &
      0.18364136980992795678D+00, &
      0.59830134980196891803D+00, &
      0.59830134980196880701D+00, &
      0.18364136980992790127D+00, &
      0.18364136980992784576D+00, &
      0.59830134980196891803D+00, &
      0.18364136980992790127D+00, &
      0.03441591057817522392D+00, &
      0.03441591057817516841D+00, &
      0.18364136980992784576D+00, &
      0.03441591057817527943D+00, &
      0.18364136980992790127D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
      0.05801054891248025314D+00, &
      0.00006431928175925639D+00, &
      0.00006431928175925639D+00, &
      0.00006431928175925639D+00, &
      0.00006431928175925639D+00, &
      0.02317333846242545722D+00, &
      0.02317333846242545722D+00, &
      0.02317333846242545722D+00, &
      0.02317333846242545722D+00, &
      0.02956291233542928526D+00, &
      0.02956291233542928526D+00, &
      0.02956291233542928526D+00, &
      0.02956291233542928526D+00, &
      0.00806397997961618221D+00, &
      0.00806397997961618221D+00, &
      0.00806397997961618221D+00, &
      0.00806397997961618221D+00, &
      0.03813408010370246404D+00, &
      0.03813408010370246404D+00, &
      0.03813408010370246404D+00, &
      0.03813408010370246404D+00, &
      0.03813408010370246404D+00, &
      0.03813408010370246404D+00, &
      0.00838442219829855194D+00, &
      0.00838442219829855194D+00, &
      0.00838442219829855194D+00, &
      0.00838442219829855194D+00, &
      0.00838442219829855194D+00, &
      0.00838442219829855194D+00, &
      0.00838442219829855194D+00, &
      0.00838442219829855194D+00, &
      0.00838442219829855194D+00, &
      0.00838442219829855194D+00, &
      0.00838442219829855194D+00, &
      0.00838442219829855194D+00, &
      0.01023455935274532845D+00, &
      0.01023455935274532845D+00, &
      0.01023455935274532845D+00, &
      0.01023455935274532845D+00, &
      0.01023455935274532845D+00, &
      0.01023455935274532845D+00, &
      0.01023455935274532845D+00, &
      0.01023455935274532845D+00, &
      0.01023455935274532845D+00, &
      0.01023455935274532845D+00, &
      0.01023455935274532845D+00, &
      0.01023455935274532845D+00, &
      0.02052491596798813878D+00, &
      0.02052491596798813878D+00, &
      0.02052491596798813878D+00, &
      0.02052491596798813878D+00, &
      0.02052491596798813878D+00, &
      0.02052491596798813878D+00, &
      0.02052491596798813878D+00, &
      0.02052491596798813878D+00, &
      0.02052491596798813878D+00, &
      0.02052491596798813878D+00, &
      0.02052491596798813878D+00, &
      0.02052491596798813878D+00 /)

  a(1:n) = a_save(1:n)
  b(1:n) = b_save(1:n)
  c(1:n) = c_save(1:n)
  d(1:n) = d_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule10 ( n, a, b, c, d, w )

!*****************************************************************************80
!
!! rule10() returns the rule of precision 10.
!
!  Discussion:
!
!    The data is given for the reference tetrahedron
!    with vertices (0,0,0), (1,0,0), (0,1,0), (0,0,1).
!
!    We suppose we are given a tetrahedron T with vertices A, B, C, D.
!    We call a rule with n points, returning barycentric coordinates
!    a, b, c, d, and weights w.  Then the integral I of f(x,y,z) over T is 
!    approximated by Q as follows:
!
!    (x,y,z) = a(1:n) * A + b(1:n) * B + c(1:n) * C + d(1:n * D
!    Q = area(T) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 April 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Freddie Witherden, Peter Vincent,
!    On the identification of symmetric quadrature rules for finite element methods,
!    Computers and Mathematics with Applications,
!    Volume 69, pages 1232-1241, 2015.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real a(n), b(n), c(n), d(n): the barycentric coordinates of quadrature points.
!
!    real w(n): the quadrature weights, which sum to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 81

  real ( kind = rk ) a(n)
  real ( kind = rk ), save, dimension ( n_save ) :: a_save = (/  &
      0.25000000000000000000D+00, &
      0.31225006869518867614D+00, &
      0.31225006869518867614D+00, &
      0.06324979391443408261D+00, &
      0.31225006869518867614D+00, &
      0.11430965385734614959D+00, &
      0.11430965385734614959D+00, &
      0.65707103842796155124D+00, &
      0.11430965385734614959D+00, &
      0.16548602561961106572D+00, &
      0.16548602561961106572D+00, &
      0.41043073921896550127D+00, &
      0.01365249594245798725D+00, &
      0.41043073921896550127D+00, &
      0.41043073921896550127D+00, &
      0.01365249594245798725D+00, &
      0.41043073921896550127D+00, &
      0.41043073921896550127D+00, &
      0.41043073921896550127D+00, &
      0.01365249594245798725D+00, &
      0.16548602561961106572D+00, &
      0.94298876734520487020D+00, &
      0.94298876734520487020D+00, &
      0.00613800882479076382D+00, &
      0.04473521500521365768D+00, &
      0.00613800882479076382D+00, &
      0.00613800882479076382D+00, &
      0.04473521500521365768D+00, &
      0.00613800882479076382D+00, &
      0.00613800882479076382D+00, &
      0.00613800882479076382D+00, &
      0.04473521500521365768D+00, &
      0.94298876734520487020D+00, &
      0.47719037990428037066D+00, &
      0.47719037990428037066D+00, &
      0.12105018114558940834D+00, &
      0.28070925780454081266D+00, &
      0.12105018114558940834D+00, &
      0.12105018114558940834D+00, &
      0.28070925780454081266D+00, &
      0.12105018114558940834D+00, &
      0.12105018114558940834D+00, &
      0.12105018114558940834D+00, &
      0.28070925780454081266D+00, &
      0.47719037990428037066D+00, &
      0.59425626948000698224D+00, &
      0.59425626948000698224D+00, &
      0.03277946821644267539D+00, &
      0.34018479408710766698D+00, &
      0.03277946821644267539D+00, &
      0.03277946821644267539D+00, &
      0.34018479408710766698D+00, &
      0.03277946821644267539D+00, &
      0.03277946821644267539D+00, &
      0.03277946821644267539D+00, &
      0.34018479408710766698D+00, &
      0.59425626948000698224D+00, &
      0.80117728465834436857D+00, &
      0.80117728465834436857D+00, &
      0.03248528156482305418D+00, &
      0.13385215221200952307D+00, &
      0.03248528156482305418D+00, &
      0.03248528156482305418D+00, &
      0.13385215221200952307D+00, &
      0.03248528156482305418D+00, &
      0.03248528156482305418D+00, &
      0.03248528156482305418D+00, &
      0.13385215221200952307D+00, &
      0.80117728465834436857D+00, &
      0.62807184547536598629D+00, &
      0.62807184547536598629D+00, &
      0.17497934218393901284D+00, &
      0.02196947015675593251D+00, &
      0.17497934218393901284D+00, &
      0.17497934218393901284D+00, &
      0.02196947015675593251D+00, &
      0.17497934218393901284D+00, &
      0.17497934218393901284D+00, &
      0.17497934218393901284D+00, &
      0.02196947015675593251D+00, &
      0.62807184547536598629D+00 /)

  real ( kind = rk ) b(n)
  real ( kind = rk ), save, dimension ( n_save ) :: b_save = (/ &
      0.25000000000000000000D+00, &
      0.31225006869518867614D+00, &
      0.06324979391443408261D+00, &
      0.31225006869518867614D+00, &
      0.31225006869518867614D+00, &
      0.11430965385734614959D+00, &
      0.65707103842796155124D+00, &
      0.11430965385734614959D+00, &
      0.11430965385734614959D+00, &
      0.41043073921896550127D+00, &
      0.41043073921896550127D+00, &
      0.41043073921896550127D+00, &
      0.16548602561961106572D+00, &
      0.01365249594245798725D+00, &
      0.16548602561961106572D+00, &
      0.41043073921896550127D+00, &
      0.01365249594245798725D+00, &
      0.41043073921896550127D+00, &
      0.16548602561961106572D+00, &
      0.41043073921896550127D+00, &
      0.01365249594245798725D+00, &
      0.00613800882479076382D+00, &
      0.00613800882479076382D+00, &
      0.00613800882479076382D+00, &
      0.94298876734520487020D+00, &
      0.04473521500521365768D+00, &
      0.94298876734520487020D+00, &
      0.00613800882479076382D+00, &
      0.04473521500521365768D+00, &
      0.00613800882479076382D+00, &
      0.94298876734520487020D+00, &
      0.00613800882479076382D+00, &
      0.04473521500521365768D+00, &
      0.12105018114558940834D+00, &
      0.12105018114558940834D+00, &
      0.12105018114558940834D+00, &
      0.47719037990428037066D+00, &
      0.28070925780454081266D+00, &
      0.47719037990428037066D+00, &
      0.12105018114558940834D+00, &
      0.28070925780454081266D+00, &
      0.12105018114558940834D+00, &
      0.47719037990428037066D+00, &
      0.12105018114558940834D+00, &
      0.28070925780454081266D+00, &
      0.03277946821644267539D+00, &
      0.03277946821644267539D+00, &
      0.03277946821644267539D+00, &
      0.59425626948000698224D+00, &
      0.34018479408710766698D+00, &
      0.59425626948000698224D+00, &
      0.03277946821644267539D+00, &
      0.34018479408710766698D+00, &
      0.03277946821644267539D+00, &
      0.59425626948000698224D+00, &
      0.03277946821644267539D+00, &
      0.34018479408710766698D+00, &
      0.03248528156482305418D+00, &
      0.03248528156482305418D+00, &
      0.03248528156482305418D+00, &
      0.80117728465834436857D+00, &
      0.13385215221200952307D+00, &
      0.80117728465834436857D+00, &
      0.03248528156482305418D+00, &
      0.13385215221200952307D+00, &
      0.03248528156482305418D+00, &
      0.80117728465834436857D+00, &
      0.03248528156482305418D+00, &
      0.13385215221200952307D+00, &
      0.17497934218393901284D+00, &
      0.17497934218393901284D+00, &
      0.17497934218393901284D+00, &
      0.62807184547536598629D+00, &
      0.02196947015675593251D+00, &
      0.62807184547536598629D+00, &
      0.17497934218393901284D+00, &
      0.02196947015675593251D+00, &
      0.17497934218393901284D+00, &
      0.62807184547536598629D+00, &
      0.17497934218393901284D+00, &
      0.02196947015675593251D+00 /)

  real ( kind = rk ) c(n)
  real ( kind = rk ), save, dimension ( n_save ) :: c_save = (/ &
      0.25000000000000000000D+00, &
      0.06324979391443408261D+00, &
      0.31225006869518867614D+00, &
      0.31225006869518867614D+00, &
      0.31225006869518867614D+00, &
      0.65707103842796155124D+00, &
      0.11430965385734614959D+00, &
      0.11430965385734614959D+00, &
      0.11430965385734614959D+00, &
      0.01365249594245798725D+00, &
      0.41043073921896550127D+00, &
      0.16548602561961106572D+00, &
      0.41043073921896550127D+00, &
      0.16548602561961106572D+00, &
      0.41043073921896550127D+00, &
      0.16548602561961106572D+00, &
      0.41043073921896550127D+00, &
      0.01365249594245798725D+00, &
      0.01365249594245798725D+00, &
      0.41043073921896550127D+00, &
      0.41043073921896550127D+00, &
      0.04473521500521365768D+00, &
      0.00613800882479076382D+00, &
      0.94298876734520487020D+00, &
      0.00613800882479076382D+00, &
      0.94298876734520487020D+00, &
      0.00613800882479076382D+00, &
      0.94298876734520487020D+00, &
      0.00613800882479076382D+00, &
      0.04473521500521365768D+00, &
      0.04473521500521365768D+00, &
      0.00613800882479076382D+00, &
      0.00613800882479076382D+00, &
      0.28070925780454081266D+00, &
      0.12105018114558940834D+00, &
      0.47719037990428037066D+00, &
      0.12105018114558940834D+00, &
      0.47719037990428037066D+00, &
      0.12105018114558940834D+00, &
      0.47719037990428037066D+00, &
      0.12105018114558940834D+00, &
      0.28070925780454081266D+00, &
      0.28070925780454081266D+00, &
      0.12105018114558940834D+00, &
      0.12105018114558940834D+00, &
      0.34018479408710766698D+00, &
      0.03277946821644267539D+00, &
      0.59425626948000698224D+00, &
      0.03277946821644267539D+00, &
      0.59425626948000698224D+00, &
      0.03277946821644267539D+00, &
      0.59425626948000698224D+00, &
      0.03277946821644267539D+00, &
      0.34018479408710766698D+00, &
      0.34018479408710766698D+00, &
      0.03277946821644267539D+00, &
      0.03277946821644267539D+00, &
      0.13385215221200952307D+00, &
      0.03248528156482305418D+00, &
      0.80117728465834436857D+00, &
      0.03248528156482305418D+00, &
      0.80117728465834436857D+00, &
      0.03248528156482305418D+00, &
      0.80117728465834436857D+00, &
      0.03248528156482305418D+00, &
      0.13385215221200952307D+00, &
      0.13385215221200952307D+00, &
      0.03248528156482305418D+00, &
      0.03248528156482305418D+00, &
      0.02196947015675593251D+00, &
      0.17497934218393901284D+00, &
      0.62807184547536598629D+00, &
      0.17497934218393901284D+00, &
      0.62807184547536598629D+00, &
      0.17497934218393901284D+00, &
      0.62807184547536598629D+00, &
      0.17497934218393901284D+00, &
      0.02196947015675593251D+00, &
      0.02196947015675593251D+00, &
      0.17497934218393901284D+00, &
      0.17497934218393901284D+00 /)

  real ( kind = rk ) d(n)
  real ( kind = rk ), save, dimension ( n_save ) :: d_save = (/ &
      0.25000000000000000000D+00, &
      0.31225006869518856512D+00, &
      0.31225006869518856512D+00, &
      0.31225006869518856512D+00, &
      0.06324979391443397159D+00, &
      0.11430965385734614959D+00, &
      0.11430965385734614959D+00, &
      0.11430965385734614959D+00, &
      0.65707103842796155124D+00, &
      0.41043073921896539025D+00, &
      0.01365249594245787623D+00, &
      0.01365249594245798725D+00, &
      0.41043073921896539025D+00, &
      0.41043073921896550127D+00, &
      0.01365249594245798725D+00, &
      0.41043073921896550127D+00, &
      0.16548602561961106572D+00, &
      0.16548602561961106572D+00, &
      0.41043073921896550127D+00, &
      0.16548602561961106572D+00, &
      0.41043073921896539025D+00, &
      0.00613800882479070831D+00, &
      0.04473521500521360217D+00, &
      0.04473521500521360217D+00, &
      0.00613800882479076382D+00, &
      0.00613800882479065280D+00, &
      0.04473521500521360217D+00, &
      0.00613800882479076382D+00, &
      0.94298876734520475917D+00, &
      0.94298876734520487020D+00, &
      0.00613800882479070831D+00, &
      0.94298876734520487020D+00, &
      0.00613800882479070831D+00, &
      0.12105018114558940834D+00, &
      0.28070925780454081266D+00, &
      0.28070925780454081266D+00, &
      0.12105018114558940834D+00, &
      0.12105018114558940834D+00, &
      0.28070925780454081266D+00, &
      0.12105018114558940834D+00, &
      0.47719037990428037066D+00, &
      0.47719037990428037066D+00, &
      0.12105018114558940834D+00, &
      0.47719037990428037066D+00, &
      0.12105018114558940834D+00, &
      0.03277946821644267539D+00, &
      0.34018479408710766698D+00, &
      0.34018479408710777800D+00, &
      0.03277946821644267539D+00, &
      0.03277946821644273090D+00, &
      0.34018479408710772249D+00, &
      0.03277946821644273090D+00, &
      0.59425626948000709326D+00, &
      0.59425626948000709326D+00, &
      0.03277946821644273090D+00, &
      0.59425626948000709326D+00, &
      0.03277946821644267539D+00, &
      0.03248528156482305418D+00, &
      0.13385215221200952307D+00, &
      0.13385215221200941205D+00, &
      0.03248528156482305418D+00, &
      0.03248528156482299867D+00, &
      0.13385215221200946756D+00, &
      0.03248528156482299867D+00, &
      0.80117728465834425755D+00, &
      0.80117728465834425755D+00, &
      0.03248528156482299867D+00, &
      0.80117728465834425755D+00, &
      0.03248528156482305418D+00, &
      0.17497934218393906836D+00, &
      0.02196947015675598802D+00, &
      0.02196947015675609904D+00, &
      0.17497934218393912387D+00, &
      0.17497934218393917938D+00, &
      0.02196947015675604353D+00, &
      0.17497934218393917938D+00, &
      0.62807184547536620833D+00, &
      0.62807184547536620833D+00, &
      0.17497934218393912387D+00, &
      0.62807184547536620833D+00, &
      0.17497934218393906836D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
      0.04739977355602073561D+00, &
      0.02693705999226870401D+00, &
      0.02693705999226870401D+00, &
      0.02693705999226870401D+00, &
      0.02693705999226870401D+00, &
      0.00986915971679338221D+00, &
      0.00986915971679338221D+00, &
      0.00986915971679338221D+00, &
      0.00986915971679338221D+00, &
      0.01139388122019523164D+00, &
      0.01139388122019523164D+00, &
      0.01139388122019523164D+00, &
      0.01139388122019523164D+00, &
      0.01139388122019523164D+00, &
      0.01139388122019523164D+00, &
      0.01139388122019523164D+00, &
      0.01139388122019523164D+00, &
      0.01139388122019523164D+00, &
      0.01139388122019523164D+00, &
      0.01139388122019523164D+00, &
      0.01139388122019523164D+00, &
      0.00036194434433925362D+00, &
      0.00036194434433925362D+00, &
      0.00036194434433925362D+00, &
      0.00036194434433925362D+00, &
      0.00036194434433925362D+00, &
      0.00036194434433925362D+00, &
      0.00036194434433925362D+00, &
      0.00036194434433925362D+00, &
      0.00036194434433925362D+00, &
      0.00036194434433925362D+00, &
      0.00036194434433925362D+00, &
      0.00036194434433925362D+00, &
      0.02573973198045606883D+00, &
      0.02573973198045606883D+00, &
      0.02573973198045606883D+00, &
      0.02573973198045606883D+00, &
      0.02573973198045606883D+00, &
      0.02573973198045606883D+00, &
      0.02573973198045606883D+00, &
      0.02573973198045606883D+00, &
      0.02573973198045606883D+00, &
      0.02573973198045606883D+00, &
      0.02573973198045606883D+00, &
      0.02573973198045606883D+00, &
      0.01013587167975579274D+00, &
      0.01013587167975579274D+00, &
      0.01013587167975579274D+00, &
      0.01013587167975579274D+00, &
      0.01013587167975579274D+00, &
      0.01013587167975579274D+00, &
      0.01013587167975579274D+00, &
      0.01013587167975579274D+00, &
      0.01013587167975579274D+00, &
      0.01013587167975579274D+00, &
      0.01013587167975579274D+00, &
      0.01013587167975579274D+00, &
      0.00657614727703590383D+00, &
      0.00657614727703590383D+00, &
      0.00657614727703590383D+00, &
      0.00657614727703590383D+00, &
      0.00657614727703590383D+00, &
      0.00657614727703590383D+00, &
      0.00657614727703590383D+00, &
      0.00657614727703590383D+00, &
      0.00657614727703590383D+00, &
      0.00657614727703590383D+00, &
      0.00657614727703590383D+00, &
      0.00657614727703590383D+00, &
      0.01290703579886198944D+00, &
      0.01290703579886198944D+00, &
      0.01290703579886198944D+00, &
      0.01290703579886198944D+00, &
      0.01290703579886198944D+00, &
      0.01290703579886198944D+00, &
      0.01290703579886198944D+00, &
      0.01290703579886198944D+00, &
      0.01290703579886198944D+00, &
      0.01290703579886198944D+00, &
      0.01290703579886198944D+00, &
      0.01290703579886198944D+00 /)

  a(1:n) = a_save(1:n)
  b(1:n) = b_save(1:n)
  c(1:n) = c_save(1:n)
  d(1:n) = d_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
function tetrahedron_unit_monomial_integral ( expon )

!*****************************************************************************80
!
!! tetrahedron_unit_monomial_integral() integrates a monomial over the unit tetrahedron.
!
!  Discussion:
!
!    This routine integrates a monomial of the form
!
!      product ( 1 <= dim <= dim_num ) x(dim)^expon(dim)
!
!    where the exponents are nonnegative integers.  Note that
!    if the combination 0^0 is encountered, it should be treated
!    as 1.
!
!    Integral ( over unit tetrahedron ) x^l y^m z^n dx dy = 
!    l! * m! * n! / ( m + n + 3 )!
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 May 2023
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
!    real ( kind = rk ) tetrahedron_unit_monomial_integral, the integral value.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer expon(3)
  integer i
  integer k
  real ( kind = rk ) tetrahedron_unit_monomial_integral
  real ( kind = rk ) value
!
!  The first computation ends with VALUE = 1.0;
!
  value = 1.0D+00

  k = 0

  do i = 1, expon(1)
    k = k + 1
!   value = value * real ( i, kind = rk ) / real ( k, kind = rk )
  end do

  do i = 1, expon(2)
    k = k + 1
    value = value * real ( i, kind = rk ) / real ( k, kind = rk )
  end do

  do i = 1, expon(3)
    k = k + 1
    value = value * real ( i, kind = rk ) / real ( k, kind = rk )
  end do

  k = k + 1
  value = value / real ( k, kind = rk )

  k = k + 1
  value = value / real ( k, kind = rk )

  k = k + 1
  value = value / real ( k, kind = rk )

  tetrahedron_unit_monomial_integral = value

  return
end
function tetrahedron_unit_volume ( )

!*****************************************************************************80
!
!! tetrahedron_unit_volume() computes the volume of a unit tetrahedron.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 May 2023
!
!  Author:
!
!    John Burkardt
!
!  Output:
!
!    real ( kind = rk ) tetrahedron_unit_volume: the volume of 
!    the unit tetrahedron.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) tetrahedron_unit_volume

  tetrahedron_unit_volume = 1.0D+00 / 6.0D+00

  return
end
function tetrahedron_volume ( tetra )

!*****************************************************************************80
!
!! tetrahedron_volume() computes the volume of a tetrahedron.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) tetra(4,3): the vertices of the tetrahedron.
!
!  Output:
!
!    real ( kind = rk ) tetrahedron_volume: the volume of the tetrahedron.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(4,4)
  real ( kind = rk ) r8mat_det_4d
  real ( kind = rk ) tetra(4,3)
  real ( kind = rk ) tetrahedron_volume

  a(1:4,1:3) = tetra(1:4,1:3)
  a(1:4,4) = 1.0D+00

  tetrahedron_volume = abs ( r8mat_det_4d ( a ) ) / 6.0D+00

  return
end
subroutine tetrahedron_witherden_rule ( p, n, a, b, c, d, w )

!*****************************************************************************80
!
!! tetrahedron_witherden_rule() returns a tetrahedron quadrature rule of given precision.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 May 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Freddie Witherden, Peter Vincent,
!    On the identification of symmetric quadrature rules for finite element methods,
!    Computers and Mathematics with Applications,
!    Volume 69, pages 1232-1241, 2015.
!
!  Input:
!
!    integer p: the precision, 0 <= p <= 10.
!
!    integer n: the order of the rule.
!
!  Output:
!
!    real a(n), b(n), c(n), d(n): the barycentric coordinates of quadrature points.
!
!    real w(n): the weights of quadrature points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) b(n)
  real ( kind = rk ) c(n)
  real ( kind = rk ) d(n)
  integer p
  real ( kind = rk ) w(n)

  if ( p < 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'tetrahedron_witherden_rule(): Fatal error!'
    write ( *, '(a)' ) '  Input p < 0.'
    call exit ( 1 )
  end if

  if ( 10 < p ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'tetrahedron_witherden_rule(): Fatal error!'
    write ( *, '(a)' ) '  Input 10 < p.'
    call exit ( 1 )
  end if

  if ( p == 0 ) then
    call rule00 ( n, a, b, c, d, w )
  else if ( p == 1 ) then
    call rule01 ( n, a, b, c, d, w )
  else if ( p == 2 ) then
    call rule02 ( n, a, b, c, d, w )
  else if ( p == 3 ) then
    call rule03 ( n, a, b, c, d, w )
  else if ( p == 4 ) then
    call rule04 ( n, a, b, c, d, w )
  else if ( p == 5 ) then
    call rule05 ( n, a, b, c, d, w )
  else if ( p == 6 ) then
    call rule06 ( n, a, b, c, d, w )
  else if ( p == 7 ) then
    call rule07 ( n, a, b, c, d, w )
  else if ( p == 8 ) then
    call rule08 ( n, a, b, c, d, w )
  else if ( p == 9 ) then
    call rule09 ( n, a, b, c, d, w )
  else if ( p == 10 ) then
    call rule10 ( n, a, b, c, d, w )
  end if

  return
end
