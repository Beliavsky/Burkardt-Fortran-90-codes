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
subroutine pyramid_witherden_rule ( p, n, x, y, z, w )

!*****************************************************************************80
!
!! pyramid_witherden_rule() returns a pyramid quadrature rule of given precision.
!
!  Discussion:
!
!    The unit pyramid with square base is the region
!
!      -1 <= X <= 1
!      -1 <= Y <= 1
!       0 <= Z <= 1.
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
!    Freddie Witherden, Peter Vincent,
!    On the identification of symmetric quadrature rules for finite element methods,
!    Computers and Mathematics with Applications,
!    Volume 69, pages 1232-1241, 2015.
!
!  Input:
!
!    int p: the precision, 0 <= p <= 10.
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
    write ( *, '(a)' ) 'pyramid_witherden_rule(): Fatal error!'
    write ( *, '(a)' ) '  Input p < 0.'
    stop ( 1 )
  end if

  if ( 10 < p ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'pyramid_witherden_rule(): Fatal error!'
    write ( *, '(a)' ) '  Input 10 < p.'
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
  end if

  return
end
subroutine pyramid_unit_monomial_integral ( expon, value )

!*****************************************************************************80
!
!! pyramid_unit_monomial_integral(): monomial integral in a unit pyramid.
!
!  Discussion:
!
!    This routine returns the integral of
!
!      product ( 1 <= I <= 3 ) X(I)^EXPON(I)
!
!    over the unit pyramid.
!
!    The integration region is:
!
!    - ( 1 - Z ) <= X <= 1 - Z
!    - ( 1 - Z ) <= Y <= 1 - Z
!              0 <= Z <= 1.
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
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
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
  integer i
  integer i_hi
  real ( kind = rk ) mop
  real ( kind = rk ) r8_choose
  real ( kind = rk ) value

  value = 0.0D+00

  if ( mod ( expon(1), 2 ) == 0 .and. mod ( expon(2), 2 ) == 0 ) then

    i_hi = 2 + expon(1) + expon(2)

    mop = 1.0D+00
    do i = 0, i_hi
      value = value + mop * r8_choose ( i_hi, i ) &
        / real ( i + expon(3) + 1, kind = rk )
      mop = - mop
    end do

    value = value &
          * 2.0D+00 / real ( expon(1) + 1, kind = rk ) &
          * 2.0D+00 / real ( expon(2) + 1, kind = rk )

  end if

  return
end
function pyramid_unit_volume ( )

!*****************************************************************************80
!
!! pyramid_unit_volume() returns the volume of a unit pyramid.
!
!  Discussion:
!
!    A pyramid with square base can be regarded as the upper half of a
!    3D octahedron.
!
!    The integration region:
!
!      - ( 1 - Z ) <= X <= 1 - Z
!      - ( 1 - Z ) <= Y <= 1 - Z
!                0 <= Z <= 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    31 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Output:
!
!    real pyramid_unit_volume: the volume of the pyramid.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) pyramid_unit_volume

  pyramid_unit_volume = 4.0D+00 / 3.0D+00

  return
end
function pyramid_volume ( r, h )

!*****************************************************************************80
!
!! pyramid_volume() returns the volume of a pyramid with square base in 3D.
!
!  Discussion:
!
!    A pyramid with square base can be regarded as the upper half of a
!    3D octahedron.
!
!    Z - R <= X <= R - Z
!    Z - R <= Y <= R - Z
!    0 <= Z <= H.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 July 2018
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real R, the "radius" of the pyramid, that is, half the
!    length of one of the sides of the square base.
!
!    real H, the height of the pyramid.
!
!  Output:
!
!    real pyramid_volume: the volume of the pyramid.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) h
  real ( kind = rk ) pyramid_volume
  real ( kind = rk ) r

  pyramid_volume = ( 4.0D+00 / 3.0D+00 ) * h * r * r

  return
end
function r8_choose ( n, k )

!*****************************************************************************80
!
!! r8_choose() computes the combinatorial coefficient C(N,K).
!
!  Discussion:
!
!    Real arithmetic is used, and C(N,K) is computed directly, via
!    Gamma functions, rather than recursively.
!
!    C(N,K) is the number of distinct combinations of K objects
!    chosen from a set of N distinct objects.  A combination is
!    like a set, in that order does not matter.
!
!    C(N,K) = N! / ( (N-K)! * K! )
!
!  Example:
!
!    The number of combinations of 2 things chosen from 5 is 10.
!
!    C(5,2) = ( 5 * 4 * 3 * 2 * 1 ) / ( ( 3 * 2 * 1 ) * ( 2 * 1 ) ) = 10.
!
!    The actual combinations may be represented as:
!
!      (1,2), (1,3), (1,4), (1,5), (2,3),
!      (2,4), (2,5), (3,4), (3,5), (4,5).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the value of N.
!
!    integer K, the value of K.
!
!  Output:
!
!    real ( kind = rk ) r8_choose: the value of C(N,K)
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) arg
  real ( kind = rk ) fack
  real ( kind = rk ) facn
  real ( kind = rk ) facnmk
  integer k
  integer n
  real ( kind = rk ) r8_choose
  real ( kind = rk ) r8_gamma_log
  real ( kind = rk ) value

  if ( n < 0 ) then

    value = 0.0D+00

  else if ( k == 0 ) then

    value = 1.0D+00

  else if ( k == 1 ) then

    value = real ( n, kind = rk )

  else if ( 1 < k .and. k < n - 1 ) then

    arg = real ( n + 1, kind = rk )
    facn = r8_gamma_log ( arg )

    arg = real ( k + 1, kind = rk )
    fack = r8_gamma_log ( arg )

    arg = real ( n - k + 1, kind = rk )
    facnmk = r8_gamma_log ( arg )

    value = real ( nint ( exp ( facn - fack - facnmk ) ),  kind = rk )

  else if ( k == n - 1 ) then

    value = real ( n, kind = rk )

  else if ( k == n ) then

    value = 1.0D+00

  else

    value = 0.0D+00

  end if

  r8_choose = value

  return
end
function r8_gamma_log ( x )

!*****************************************************************************80
!
!! r8_gamma_log() evaluates the logarithm of the gamma function.
!
!  Discussion:
!
!    This routine calculates the LOG(GAMMA) function for a positive real
!    argument X.  Computation is based on an algorithm outlined in
!    references 1 and 2.  The program uses rational functions that
!    theoretically approximate LOG(GAMMA) to at least 18 significant
!    decimal digits.  The approximation for X > 12 is from reference
!    3, while approximations for X < 12.0 are similar to those in
!    reference 1, but are unpublished.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 April 2013
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody, Kenneth Hillstrom,
!    Chebyshev Approximations for the Natural Logarithm of the
!    Gamma Function,
!    Mathematics of Computation,
!    Volume 21, Number 98, April 1967, pages 198-203.
!
!    Kenneth Hillstrom,
!    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
!    May 1969.
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thatcher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Input:
!
!    real ( kind = rk ) X, the argument of the function.
!
!  Output:
!
!    real ( kind = rk ) R8_GAMMA_LOG, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), dimension ( 7 ) :: c = (/ &
    -1.910444077728D-03, &
     8.4171387781295D-04, &
    -5.952379913043012D-04, &
     7.93650793500350248D-04, &
    -2.777777777777681622553D-03, &
     8.333333333333333331554247D-02, &
     5.7083835261D-03 /)
  real ( kind = rk ) corr
  real ( kind = rk ) :: d1 = -5.772156649015328605195174D-01
  real ( kind = rk ) :: d2 = 4.227843350984671393993777D-01
  real ( kind = rk ) :: d4 = 1.791759469228055000094023D+00
  real ( kind = rk ), parameter :: frtbig = 2.25D+76
  integer i
  real ( kind = rk ), dimension ( 8 ) :: p1 = (/ &
    4.945235359296727046734888D+00, &
    2.018112620856775083915565D+02, &
    2.290838373831346393026739D+03, &
    1.131967205903380828685045D+04, &
    2.855724635671635335736389D+04, &
    3.848496228443793359990269D+04, &
    2.637748787624195437963534D+04, &
    7.225813979700288197698961D+03 /)
  real ( kind = rk ), dimension ( 8 ) :: p2 = (/ &
    4.974607845568932035012064D+00, &
    5.424138599891070494101986D+02, &
    1.550693864978364947665077D+04, &
    1.847932904445632425417223D+05, &
    1.088204769468828767498470D+06, &
    3.338152967987029735917223D+06, &
    5.106661678927352456275255D+06, &
    3.074109054850539556250927D+06 /)
  real ( kind = rk ), dimension ( 8 ) :: p4 = (/ &
    1.474502166059939948905062D+04, &
    2.426813369486704502836312D+06, &
    1.214755574045093227939592D+08, &
    2.663432449630976949898078D+09, &
    2.940378956634553899906876D+10, &
    1.702665737765398868392998D+11, &
    4.926125793377430887588120D+11, &
    5.606251856223951465078242D+11 /)
  real ( kind = rk ), dimension ( 8 ) :: q1 = (/ &
    6.748212550303777196073036D+01, &
    1.113332393857199323513008D+03, &
    7.738757056935398733233834D+03, &
    2.763987074403340708898585D+04, &
    5.499310206226157329794414D+04, &
    6.161122180066002127833352D+04, &
    3.635127591501940507276287D+04, &
    8.785536302431013170870835D+03 /)
  real ( kind = rk ), dimension ( 8 ) :: q2 = (/ &
    1.830328399370592604055942D+02, &
    7.765049321445005871323047D+03, &
    1.331903827966074194402448D+05, &
    1.136705821321969608938755D+06, &
    5.267964117437946917577538D+06, &
    1.346701454311101692290052D+07, &
    1.782736530353274213975932D+07, &
    9.533095591844353613395747D+06 /)
  real ( kind = rk ), dimension ( 8 ) :: q4 = (/ &
    2.690530175870899333379843D+03, &
    6.393885654300092398984238D+05, &
    4.135599930241388052042842D+07, &
    1.120872109616147941376570D+09, &
    1.488613728678813811542398D+10, &
    1.016803586272438228077304D+11, &
    3.417476345507377132798597D+11, &
    4.463158187419713286462081D+11 /)
  real ( kind = rk ) r8_gamma_log
  real ( kind = rk ) res
  real ( kind = rk ), parameter :: sqrtpi = 0.9189385332046727417803297D+00
  real ( kind = rk ) x
  real ( kind = rk ), parameter :: xbig = 2.55D+305
  real ( kind = rk ) xden
  real ( kind = rk ), parameter :: xinf = 1.79D+308
  real ( kind = rk ) xm1
  real ( kind = rk ) xm2
  real ( kind = rk ) xm4
  real ( kind = rk ) xnum
  real ( kind = rk ) y
  real ( kind = rk ) ysq

  y = x

  if ( 0.0D+00 < y .and. y <= xbig ) then

    if ( y <= epsilon ( y ) ) then

      res = - log ( y )
!
!  EPS < X <= 1.5.
!
    else if ( y <= 1.5D+00 ) then

      if ( y < 0.6796875D+00 ) then
        corr = -log ( y )
        xm1 = y
      else
        corr = 0.0D+00
        xm1 = ( y - 0.5D+00 ) - 0.5D+00
      end if

      if ( y <= 0.5D+00 .or. 0.6796875D+00 <= y ) then

        xden = 1.0D+00
        xnum = 0.0D+00
        do i = 1, 8
          xnum = xnum * xm1 + p1(i)
          xden = xden * xm1 + q1(i)
        end do

        res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) )

      else

        xm2 = ( y - 0.5D+00 ) - 0.5D+00
        xden = 1.0D+00
        xnum = 0.0D+00
        do i = 1, 8
          xnum = xnum * xm2 + p2(i)
          xden = xden * xm2 + q2(i)
        end do

        res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) )

      end if
!
!  1.5 < X <= 4.0.
!
    else if ( y <= 4.0D+00 ) then

      xm2 = y - 2.0D+00
      xden = 1.0D+00
      xnum = 0.0D+00
      do i = 1, 8
        xnum = xnum * xm2 + p2(i)
        xden = xden * xm2 + q2(i)
      end do

      res = xm2 * ( d2 + xm2 * ( xnum / xden ) )
!
!  4.0 < X <= 12.0.
!
    else if ( y <= 12.0D+00 ) then

      xm4 = y - 4.0D+00
      xden = -1.0D+00
      xnum = 0.0D+00
      do i = 1, 8
        xnum = xnum * xm4 + p4(i)
        xden = xden * xm4 + q4(i)
      end do

      res = d4 + xm4 * ( xnum / xden )
!
!  Evaluate for 12 <= argument.
!
    else

      res = 0.0D+00

      if ( y <= frtbig ) then

        res = c(7)
        ysq = y * y

        do i = 1, 6
          res = res / ysq + c(i)
        end do

      end if

      res = res / y
      corr = log ( y )
      res = res + sqrtpi - 0.5D+00 * corr
      res = res + y * ( corr - 1.0D+00 )

    end if
!
!  Return for bad arguments.
!
  else

    res = xinf

  end if
!
!  Final adjustments and return.
!
  r8_gamma_log = res

  return
end
subroutine rule_order ( p, order )

!*****************************************************************************80
!
!! rule_order() returns the order of a pyramid quadrature rule of given precision.
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
!    25 April 2023
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
  integer, dimension ( 0:10 ) :: order_vec = (/ &
      1, &
      1,   5,   6,  10,  15, &
     24,  31,  47,  62,  83 /)
  integer p

  if ( p < 0 ) then
    write ( *, '(a)' ) 'rule_order(): Input p < 0.'
  end if

  if ( 10 < p ) then
    write ( *, '(a)' ) 'rule_order(): Input 10 < p.'
  end if

  order = order_vec(p)

  return
end
subroutine rule00 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule00() returns the pyramid quadrature rule of precision 0.
!
!  Discussion:
!
!    We suppose we are given a pyramid P with vertices
!      (-1,-1,0), (-1,+1,0), (+1,+1,0), (+1,-1,0), (0,0,1).
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over
!    P is approximated by Q as follows:
!
!    Q = volume(P) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 April 2023
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
          0.0000000000000000D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.0000000000000000D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.2500000000000000D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          1.0000000000000000D+00 /)

  x = x_save
  y = y_save
  z = z_save
  w = w_save

  return
end
subroutine rule01 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule01() returns the pyramid quadrature rule of precision 1.
!
!  Discussion:
!
!    We suppose we are given a pyramid P with vertices
!      (-1,-1,0), (-1,+1,0), (+1,+1,0), (+1,-1,0), (0,0,1).
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over
!    P is approximated by Q as follows:
!
!    Q = volume(P) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 April 2023
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
          0.0000000000000000D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.0000000000000000D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.2500000000000000D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          1.0000000000000000D+00 /)

  x = x_save
  y = y_save
  z = z_save
  w = w_save

  return
end
subroutine rule02 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule02() returns the pyramid quadrature rule of precision 2.
!
!  Discussion:
!
!    We suppose we are given a pyramid P with vertices
!      (-1,-1,0), (-1,+1,0), (+1,+1,0), (+1,-1,0), (0,0,1).
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over
!    P is approximated by Q as follows:
!
!    Q = volume(P) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 April 2023
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
          0.0000000000000000D+00, &
          0.7189210558117962D+00, &
          0.0000000000000000D+00, &
         -0.7189210558117962D+00, &
          0.0000000000000000D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.0000000000000000D+00, &
          0.0000000000000000D+00, &
          0.7189210558117962D+00, &
          0.0000000000000000D+00, &
         -0.7189210558117962D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.6082910385597788D+00, &
          0.1453364835728557D+00, &
          0.1453364835728557D+00, &
          0.1453364835728557D+00, &
          0.1453364835728557D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          0.2260773013241023D+00, &
          0.1934806746689745D+00, &
          0.1934806746689745D+00, &
          0.1934806746689745D+00, &
          0.1934806746689745D+00 /)

  x = x_save
  y = y_save
  z = z_save
  w = w_save

  return
end
subroutine rule03 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule03() returns the pyramid quadrature rule of precision 3.
!
!  Discussion:
!
!    We suppose we are given a pyramid P with vertices
!      (-1,-1,0), (-1,+1,0), (+1,+1,0), (+1,-1,0), (0,0,1).
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over
!    P is approximated by Q as follows:
!
!    Q = volume(P) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 April 2023
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
!    real x(n), y(n), z(n): the coordinates of quadrature points.
!
!    real w(n): the quadrature weights, which sum to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: n_save = 6

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
          0.0000000000000000D+00, &
          0.0000000000000000D+00, &
          0.5610836110587396D+00, &
          0.5610836110587396D+00, &
         -0.5610836110587396D+00, &
         -0.5610836110587396D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.0000000000000000D+00, &
          0.0000000000000000D+00, &
          0.5610836110587396D+00, &
         -0.5610836110587396D+00, &
          0.5610836110587396D+00, &
         -0.5610836110587396D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.5714285703860683D+00, &
          0.0000000056758500D+00, &
          0.1666666666666667D+00, &
          0.1666666666666667D+00, &
          0.1666666666666667D+00, &
          0.1666666666666667D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          0.2522058839227606D+00, &
          0.1125000060660650D+00, &
          0.1588235275027936D+00, &
          0.1588235275027936D+00, &
          0.1588235275027936D+00, &
          0.1588235275027936D+00 /)

  x = x_save
  y = y_save
  z = z_save
  w = w_save

  return
end
subroutine rule04 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule04() returns the pyramid quadrature rule of precision 4.
!
!  Discussion:
!
!    We suppose we are given a pyramid P with vertices
!      (-1,-1,0), (-1,+1,0), (+1,+1,0), (+1,-1,0), (0,0,1).
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over
!    P is approximated by Q as follows:
!
!    Q = volume(P) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 April 2023
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
!    real x(n), y(n), z(n): the coordinates of quadrature points.
!
!    real w(n): the quadrature weights, which sum to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: n_save = 10

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
          0.0000000000000000D+00, &
          0.0000000000000000D+00, &
          0.6505815563982326D+00, &
          0.0000000000000000D+00, &
         -0.6505815563982326D+00, &
          0.0000000000000000D+00, &
          0.6579669971216900D+00, &
          0.6579669971216900D+00, &
         -0.6579669971216900D+00, &
         -0.6579669971216900D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.0000000000000000D+00, &
          0.0000000000000000D+00, &
          0.0000000000000000D+00, &
          0.6505815563982326D+00, &
          0.0000000000000000D+00, &
         -0.6505815563982326D+00, &
          0.6579669971216900D+00, &
         -0.6579669971216900D+00, &
          0.6579669971216900D+00, &
         -0.6579669971216900D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.6772327888861374D+00, &
          0.1251369531087465D+00, &
          0.3223841495782137D+00, &
          0.3223841495782137D+00, &
          0.3223841495782137D+00, &
          0.3223841495782137D+00, &
          0.0392482838988154D+00, &
          0.0392482838988154D+00, &
          0.0392482838988154D+00, &
          0.0392482838988154D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          0.1137418831706419D+00, &
          0.2068834025895523D+00, &
          0.1063245878893255D+00, &
          0.1063245878893255D+00, &
          0.1063245878893255D+00, &
          0.1063245878893255D+00, &
          0.0635190906706259D+00, &
          0.0635190906706259D+00, &
          0.0635190906706259D+00, &
          0.0635190906706259D+00 /)

  x = x_save
  y = y_save
  z = z_save
  w = w_save

  return
end
subroutine rule05 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule05() returns the pyramid quadrature rule of precision 5.
!
!  Discussion:
!
!    We suppose we are given a pyramid P with vertices
!      (-1,-1,0), (-1,+1,0), (+1,+1,0), (+1,-1,0), (0,0,1).
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over
!    P is approximated by Q as follows:
!
!    Q = volume(P) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 April 2023
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
!    real x(n), y(n), z(n): the coordinates of quadrature points.
!
!    real w(n): the quadrature weights, which sum to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: n_save = 15

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
          0.0000000000000000D+00, &
          0.0000000000000000D+00, &
          0.0000000000000000D+00, &
          0.7065260315463245D+00, &
          0.0000000000000000D+00, &
         -0.7065260315463245D+00, &
          0.0000000000000000D+00, &
          0.7051171227788277D+00, &
          0.7051171227788277D+00, &
         -0.7051171227788277D+00, &
         -0.7051171227788277D+00, &
          0.4328828641035410D+00, &
          0.4328828641035410D+00, &
         -0.4328828641035410D+00, &
         -0.4328828641035410D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.0000000000000000D+00, &
          0.0000000000000000D+00, &
          0.0000000000000000D+00, &
          0.0000000000000000D+00, &
          0.7065260315463245D+00, &
          0.0000000000000000D+00, &
         -0.7065260315463245D+00, &
          0.7051171227788277D+00, &
         -0.7051171227788277D+00, &
          0.7051171227788277D+00, &
         -0.7051171227788277D+00, &
          0.4328828641035410D+00, &
         -0.4328828641035410D+00, &
          0.4328828641035410D+00, &
         -0.4328828641035410D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.7298578807825067D+00, &
          0.3004010208137690D+00, &
          0.0000000064917722D+00, &
          0.1250000000000000D+00, &
          0.1250000000000000D+00, &
          0.1250000000000000D+00, &
          0.1250000000000000D+00, &
          0.0611119070620230D+00, &
          0.0611119070620230D+00, &
          0.0611119070620230D+00, &
          0.0611119070620230D+00, &
          0.4236013371197248D+00, &
          0.4236013371197248D+00, &
          0.4236013371197248D+00, &
          0.4236013371197248D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          0.0684353699091401D+00, &
          0.1693971144927240D+00, &
          0.0587045358285745D+00, &
          0.0764412931481202D+00, &
          0.0764412931481202D+00, &
          0.0764412931481202D+00, &
          0.0764412931481202D+00, &
          0.0396709015796455D+00, &
          0.0396709015796455D+00, &
          0.0396709015796455D+00, &
          0.0396709015796455D+00, &
          0.0597535502146247D+00, &
          0.0597535502146247D+00, &
          0.0597535502146247D+00, &
          0.0597535502146247D+00 /)

  x = x_save
  y = y_save
  z = z_save
  w = w_save

  return
end
subroutine rule06 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule06() returns the pyramid quadrature rule of precision 6.
!
!  Discussion:
!
!    We suppose we are given a pyramid P with vertices
!      (-1,-1,0), (-1,+1,0), (+1,+1,0), (+1,-1,0), (0,0,1).
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over
!    P is approximated by Q as follows:
!
!    Q = volume(P) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 April 2023
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
!    real x(n), y(n), z(n): the coordinates of quadrature points.
!
!    real w(n): the quadrature weights, which sum to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: n_save = 24

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
          0.0000000000000000D+00, &
          0.0000000000000000D+00, &
          0.0000000000000000D+00, &
          0.0000000000000000D+00, &
          0.8345953511147084D+00, &
          0.0000000000000000D+00, &
         -0.8345953511147084D+00, &
          0.0000000000000000D+00, &
          0.4339254093766991D+00, &
          0.0000000000000000D+00, &
         -0.4339254093766991D+00, &
          0.0000000000000000D+00, &
          0.5656808544256755D+00, &
          0.5656808544256755D+00, &
         -0.5656808544256755D+00, &
         -0.5656808544256755D+00, &
          0.4980790917807059D+00, &
          0.4980790917807059D+00, &
         -0.4980790917807059D+00, &
         -0.4980790917807059D+00, &
          0.9508994872144825D+00, &
          0.9508994872144825D+00, &
         -0.9508994872144825D+00, &
         -0.9508994872144825D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.0000000000000000D+00, &
          0.0000000000000000D+00, &
          0.0000000000000000D+00, &
          0.0000000000000000D+00, &
          0.0000000000000000D+00, &
          0.8345953511147084D+00, &
          0.0000000000000000D+00, &
         -0.8345953511147084D+00, &
          0.0000000000000000D+00, &
          0.4339254093766991D+00, &
          0.0000000000000000D+00, &
         -0.4339254093766991D+00, &
          0.5656808544256755D+00, &
         -0.5656808544256755D+00, &
          0.5656808544256755D+00, &
         -0.5656808544256755D+00, &
          0.4980790917807059D+00, &
         -0.4980790917807059D+00, &
          0.4980790917807059D+00, &
         -0.4980790917807059D+00, &
          0.9508994872144825D+00, &
         -0.9508994872144825D+00, &
          0.9508994872144825D+00, &
         -0.9508994872144825D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.8076457976939595D+00, &
          0.0017638088528196D+00, &
          0.1382628064637306D+00, &
          0.4214239119356371D+00, &
          0.0974473410254620D+00, &
          0.0974473410254620D+00, &
          0.0974473410254620D+00, &
          0.0974473410254620D+00, &
          0.5660745906233009D+00, &
          0.5660745906233009D+00, &
          0.5660745906233009D+00, &
          0.5660745906233009D+00, &
          0.0294777308457207D+00, &
          0.0294777308457207D+00, &
          0.0294777308457207D+00, &
          0.0294777308457207D+00, &
          0.2649158632121295D+00, &
          0.2649158632121295D+00, &
          0.2649158632121295D+00, &
          0.2649158632121295D+00, &
          0.0482490706319360D+00, &
          0.0482490706319360D+00, &
          0.0482490706319360D+00, &
          0.0482490706319360D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          0.0254628936626420D+00, &
          0.0160535131751913D+00, &
          0.1195795544525238D+00, &
          0.1030606701991518D+00, &
          0.0369505879363295D+00, &
          0.0369505879363295D+00, &
          0.0369505879363295D+00, &
          0.0369505879363295D+00, &
          0.0315875794881733D+00, &
          0.0315875794881733D+00, &
          0.0315875794881733D+00, &
          0.0315875794881733D+00, &
          0.0372001293894483D+00, &
          0.0372001293894483D+00, &
          0.0372001293894483D+00, &
          0.0372001293894483D+00, &
          0.0738823846769269D+00, &
          0.0738823846769269D+00, &
          0.0738823846769269D+00, &
          0.0738823846769269D+00, &
          0.0043401606367449D+00, &
          0.0043401606367449D+00, &
          0.0043401606367449D+00, &
          0.0043401606367449D+00 /)

  x = x_save
  y = y_save
  z = z_save
  w = w_save

  return
end
subroutine rule07 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule07() returns the pyramid quadrature rule of precision 7.
!
!  Discussion:
!
!    We suppose we are given a pyramid P with vertices
!      (-1,-1,0), (-1,+1,0), (+1,+1,0), (+1,-1,0), (0,0,1).
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over
!    P is approximated by Q as follows:
!
!    Q = volume(P) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 April 2023
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
!    real x(n), y(n), z(n): the coordinates of quadrature points.
!
!    real w(n): the quadrature weights, which sum to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: n_save = 31

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
          0.0000000000000000D+00, &
          0.0000000000000000D+00, &
          0.0000000000000000D+00, &
          0.8640987597877147D+00, &
          0.0000000000000000D+00, &
         -0.8640987597877147D+00, &
          0.0000000000000000D+00, &
          0.6172133998483676D+00, &
          0.0000000000000000D+00, &
         -0.6172133998483676D+00, &
          0.0000000000000000D+00, &
          0.3541523808681161D+00, &
          0.3541523808681161D+00, &
         -0.3541523808681161D+00, &
         -0.3541523808681161D+00, &
          0.8027258501000878D+00, &
          0.8027258501000878D+00, &
         -0.8027258501000878D+00, &
         -0.8027258501000878D+00, &
          0.2541353468618572D+00, &
          0.2541353468618572D+00, &
         -0.2541353468618572D+00, &
         -0.2541353468618572D+00, &
          0.6143051077207853D+00, &
          0.6143051077207853D+00, &
         -0.6143051077207853D+00, &
         -0.6143051077207853D+00, &
          0.5248326982543755D+00, &
          0.5248326982543755D+00, &
         -0.5248326982543755D+00, &
         -0.5248326982543755D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.0000000000000000D+00, &
          0.0000000000000000D+00, &
          0.0000000000000000D+00, &
          0.0000000000000000D+00, &
          0.8640987597877147D+00, &
          0.0000000000000000D+00, &
         -0.8640987597877147D+00, &
          0.0000000000000000D+00, &
          0.6172133998483676D+00, &
          0.0000000000000000D+00, &
         -0.6172133998483676D+00, &
          0.3541523808681161D+00, &
         -0.3541523808681161D+00, &
          0.3541523808681161D+00, &
         -0.3541523808681161D+00, &
          0.8027258501000878D+00, &
         -0.8027258501000878D+00, &
          0.8027258501000878D+00, &
         -0.8027258501000878D+00, &
          0.2541353468618572D+00, &
         -0.2541353468618572D+00, &
          0.2541353468618572D+00, &
         -0.2541353468618572D+00, &
          0.6143051077207853D+00, &
         -0.6143051077207853D+00, &
          0.6143051077207853D+00, &
         -0.6143051077207853D+00, &
          0.5248326982543755D+00, &
         -0.5248326982543755D+00, &
          0.5248326982543755D+00, &
         -0.5248326982543755D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.0000454802082028D+00, &
          0.3935828767576542D+00, &
          0.8386827828168515D+00, &
          0.0666666666666667D+00, &
          0.0666666666666667D+00, &
          0.0666666666666667D+00, &
          0.0666666666666667D+00, &
          0.3333333333333334D+00, &
          0.3333333333333334D+00, &
          0.3333333333333334D+00, &
          0.3333333333333334D+00, &
          0.1292864095395495D+00, &
          0.1292864095395495D+00, &
          0.1292864095395495D+00, &
          0.1292864095395495D+00, &
          0.0801385301797781D+00, &
          0.0801385301797781D+00, &
          0.0801385301797781D+00, &
          0.0801385301797781D+00, &
          0.6055199301105999D+00, &
          0.6055199301105999D+00, &
          0.6055199301105999D+00, &
          0.6055199301105999D+00, &
          0.0000000000897583D+00, &
          0.0000000000897583D+00, &
          0.0000000000897583D+00, &
          0.0000000000897583D+00, &
          0.2905430754945976D+00, &
          0.2905430754945976D+00, &
          0.2905430754945976D+00, &
          0.2905430754945976D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          0.0250005770341970D+00, &
          0.1005301826998931D+00, &
          0.0157071744270173D+00, &
          0.0266917592930029D+00, &
          0.0266917592930029D+00, &
          0.0266917592930029D+00, &
          0.0266917592930029D+00, &
          0.0287109375000000D+00, &
          0.0287109375000000D+00, &
          0.0287109375000000D+00, &
          0.0287109375000000D+00, &
          0.0659601871499773D+00, &
          0.0659601871499773D+00, &
          0.0659601871499773D+00, &
          0.0659601871499773D+00, &
          0.0147861379562357D+00, &
          0.0147861379562357D+00, &
          0.0147861379562357D+00, &
          0.0147861379562357D+00, &
          0.0295191096942239D+00, &
          0.0295191096942239D+00, &
          0.0295191096942239D+00, &
          0.0295191096942239D+00, &
          0.0133038449696241D+00, &
          0.0133038449696241D+00, &
          0.0133038449696241D+00, &
          0.0133038449696241D+00, &
          0.0357185398966593D+00, &
          0.0357185398966593D+00, &
          0.0357185398966593D+00, &
          0.0357185398966593D+00 /)

  x = x_save
  y = y_save
  z = z_save
  w = w_save

  return
end
subroutine rule08 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule08() returns the pyramid quadrature rule of precision 8.
!
!  Discussion:
!
!    We suppose we are given a pyramid P with vertices
!      (-1,-1,0), (-1,+1,0), (+1,+1,0), (+1,-1,0), (0,0,1).
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over
!    P is approximated by Q as follows:
!
!    Q = volume(P) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 April 2023
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
!    real x(n), y(n), z(n): the coordinates of quadrature points.
!
!    real w(n): the quadrature weights, which sum to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: n_save = 47

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
          0.0000000000000000D+00, &
          0.0000000000000000D+00, &
          0.0000000000000000D+00, &
          0.7960960887742582D+00, &
          0.0000000000000000D+00, &
         -0.7960960887742582D+00, &
          0.0000000000000000D+00, &
          0.5648889587418399D+00, &
          0.0000000000000000D+00, &
         -0.5648889587418399D+00, &
          0.0000000000000000D+00, &
          0.7530226935489476D+00, &
          0.0000000000000000D+00, &
         -0.7530226935489476D+00, &
          0.0000000000000000D+00, &
          0.2465779424622870D+00, &
          0.2465779424622870D+00, &
         -0.2465779424622870D+00, &
         -0.2465779424622870D+00, &
          0.7338895325467361D+00, &
          0.7338895325467361D+00, &
         -0.7338895325467361D+00, &
         -0.7338895325467361D+00, &
          0.0000826826715128D+00, &
          0.0000826826715128D+00, &
         -0.0000826826715128D+00, &
         -0.0000826826715128D+00, &
          0.3989965570952277D+00, &
          0.3989965570952277D+00, &
         -0.3989965570952277D+00, &
         -0.3989965570952277D+00, &
          0.3653124144795900D+00, &
          0.3653124144795900D+00, &
         -0.3653124144795900D+00, &
         -0.3653124144795900D+00, &
          0.4834252440408382D+00, &
          0.4834252440408382D+00, &
         -0.4834252440408382D+00, &
         -0.4834252440408382D+00, &
          0.6314047976005609D+00, &
          0.8948561964915306D+00, &
          0.6314047976005609D+00, &
         -0.8948561964915306D+00, &
         -0.6314047976005609D+00, &
          0.8948561964915306D+00, &
         -0.6314047976005609D+00, &
         -0.8948561964915306D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.0000000000000000D+00, &
          0.0000000000000000D+00, &
          0.0000000000000000D+00, &
          0.0000000000000000D+00, &
          0.7960960887742582D+00, &
          0.0000000000000000D+00, &
         -0.7960960887742582D+00, &
          0.0000000000000000D+00, &
          0.5648889587418399D+00, &
          0.0000000000000000D+00, &
         -0.5648889587418399D+00, &
          0.0000000000000000D+00, &
          0.7530226935489476D+00, &
          0.0000000000000000D+00, &
         -0.7530226935489476D+00, &
          0.2465779424622870D+00, &
         -0.2465779424622870D+00, &
          0.2465779424622870D+00, &
         -0.2465779424622870D+00, &
          0.7338895325467361D+00, &
         -0.7338895325467361D+00, &
          0.7338895325467361D+00, &
         -0.7338895325467361D+00, &
          0.0000826826715128D+00, &
         -0.0000826826715128D+00, &
          0.0000826826715128D+00, &
         -0.0000826826715128D+00, &
          0.3989965570952277D+00, &
         -0.3989965570952277D+00, &
          0.3989965570952277D+00, &
         -0.3989965570952277D+00, &
          0.3653124144795900D+00, &
         -0.3653124144795900D+00, &
          0.3653124144795900D+00, &
         -0.3653124144795900D+00, &
          0.4834252440408382D+00, &
         -0.4834252440408382D+00, &
          0.4834252440408382D+00, &
         -0.4834252440408382D+00, &
          0.8948561964915306D+00, &
          0.6314047976005609D+00, &
         -0.8948561964915306D+00, &
          0.6314047976005609D+00, &
          0.8948561964915306D+00, &
         -0.6314047976005609D+00, &
         -0.8948561964915306D+00, &
         -0.6314047976005609D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.5848054341081260D+00, &
          0.2714503954191895D+00, &
          0.0714404302154540D+00, &
          0.0000000073236390D+00, &
          0.0000000073236390D+00, &
          0.0000000073236390D+00, &
          0.0000000073236390D+00, &
          0.4221837006010586D+00, &
          0.4221837006010586D+00, &
          0.4221837006010586D+00, &
          0.4221837006010586D+00, &
          0.1477639972965497D+00, &
          0.1477639972965497D+00, &
          0.1477639972965497D+00, &
          0.1477639972965497D+00, &
          0.6734528479254847D+00, &
          0.6734528479254847D+00, &
          0.6734528479254847D+00, &
          0.6734528479254847D+00, &
          0.2455708050362506D+00, &
          0.2455708050362506D+00, &
          0.2455708050362506D+00, &
          0.2455708050362506D+00, &
          0.8582503626505816D+00, &
          0.8582503626505816D+00, &
          0.8582503626505816D+00, &
          0.8582503626505816D+00, &
          0.0474032316194875D+00, &
          0.0474032316194875D+00, &
          0.0474032316194875D+00, &
          0.0474032316194875D+00, &
          0.4021857062205061D+00, &
          0.4021857062205061D+00, &
          0.4021857062205061D+00, &
          0.4021857062205061D+00, &
          0.1903148210864091D+00, &
          0.1903148210864091D+00, &
          0.1903148210864091D+00, &
          0.1903148210864091D+00, &
          0.0454758042237511D+00, &
          0.0454758042237511D+00, &
          0.0454758042237511D+00, &
          0.0454758042237511D+00, &
          0.0454758042237511D+00, &
          0.0454758042237511D+00, &
          0.0454758042237511D+00, &
          0.0454758042237511D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          0.0494409951557794D+00, &
          0.0873546894411460D+00, &
          0.0384605774968570D+00, &
          0.0093760594001275D+00, &
          0.0093760594001275D+00, &
          0.0093760594001275D+00, &
          0.0093760594001275D+00, &
          0.0163110192164366D+00, &
          0.0163110192164366D+00, &
          0.0163110192164366D+00, &
          0.0163110192164366D+00, &
          0.0322397220710709D+00, &
          0.0322397220710709D+00, &
          0.0322397220710709D+00, &
          0.0322397220710709D+00, &
          0.0126350996218853D+00, &
          0.0126350996218853D+00, &
          0.0126350996218853D+00, &
          0.0126350996218853D+00, &
          0.0062155450261151D+00, &
          0.0062155450261151D+00, &
          0.0062155450261151D+00, &
          0.0062155450261151D+00, &
          0.0026132528336493D+00, &
          0.0026132528336493D+00, &
          0.0026132528336493D+00, &
          0.0026132528336493D+00, &
          0.0293249681950410D+00, &
          0.0293249681950410D+00, &
          0.0293249681950410D+00, &
          0.0293249681950410D+00, &
          0.0363593842641747D+00, &
          0.0363593842641747D+00, &
          0.0363593842641747D+00, &
          0.0363593842641747D+00, &
          0.0409820609459212D+00, &
          0.0409820609459212D+00, &
          0.0409820609459212D+00, &
          0.0409820609459212D+00, &
          0.0100644114510665D+00, &
          0.0100644114510665D+00, &
          0.0100644114510665D+00, &
          0.0100644114510665D+00, &
          0.0100644114510665D+00, &
          0.0100644114510665D+00, &
          0.0100644114510665D+00, &
          0.0100644114510665D+00 /)

  x = x_save
  y = y_save
  z = z_save
  w = w_save

  return
end
subroutine rule09 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule09() returns the pyramid quadrature rule of precision 9.
!
!  Discussion:
!
!    We suppose we are given a pyramid P with vertices
!      (-1,-1,0), (-1,+1,0), (+1,+1,0), (+1,-1,0), (0,0,1).
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over
!    P is approximated by Q as follows:
!
!    Q = volume(P) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 April 2023
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
!    real x(n), y(n), z(n): the coordinates of quadrature points.
!
!    real w(n): the quadrature weights, which sum to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: n_save = 62

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
          0.0000000000000000D+00, &
          0.0000000000000000D+00, &
          0.9258200958608426D+00, &
          0.0000000000000000D+00, &
         -0.9258200958608426D+00, &
          0.0000000000000000D+00, &
          0.6760969921163639D+00, &
          0.0000000000000000D+00, &
         -0.6760969921163639D+00, &
          0.0000000000000000D+00, &
          0.6163321778509252D+00, &
          0.0000000000000000D+00, &
         -0.6163321778509252D+00, &
          0.0000000000000000D+00, &
          0.4319711692851803D+00, &
          0.0000000000000000D+00, &
         -0.4319711692851803D+00, &
          0.0000000000000000D+00, &
          0.2418018136307699D+00, &
          0.2418018136307699D+00, &
         -0.2418018136307699D+00, &
         -0.2418018136307699D+00, &
          0.8440791629031899D+00, &
          0.8440791629031899D+00, &
         -0.8440791629031899D+00, &
         -0.8440791629031899D+00, &
          0.2279669387608253D+00, &
          0.2279669387608253D+00, &
         -0.2279669387608253D+00, &
         -0.2279669387608253D+00, &
          0.5027849771350101D+00, &
          0.5027849771350101D+00, &
         -0.5027849771350101D+00, &
         -0.5027849771350101D+00, &
          0.2605010749834311D+00, &
          0.2605010749834311D+00, &
         -0.2605010749834311D+00, &
         -0.2605010749834311D+00, &
          0.0926958730867243D+00, &
          0.0926958730867243D+00, &
         -0.0926958730867243D+00, &
         -0.0926958730867243D+00, &
          0.4832161680706445D+00, &
          0.4832161680706445D+00, &
         -0.4832161680706445D+00, &
         -0.4832161680706445D+00, &
          0.5671855056082324D+00, &
          0.5671855056082324D+00, &
         -0.5671855056082324D+00, &
         -0.5671855056082324D+00, &
          0.8463949915138761D+00, &
          0.8463949915138761D+00, &
         -0.8463949915138761D+00, &
         -0.8463949915138761D+00, &
          0.8471778681177453D+00, &
          0.4641907129661463D+00, &
          0.8471778681177453D+00, &
         -0.4641907129661463D+00, &
         -0.8471778681177453D+00, &
          0.4641907129661463D+00, &
         -0.8471778681177453D+00, &
         -0.4641907129661463D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.0000000000000000D+00, &
          0.0000000000000000D+00, &
          0.0000000000000000D+00, &
          0.9258200958608426D+00, &
          0.0000000000000000D+00, &
         -0.9258200958608426D+00, &
          0.0000000000000000D+00, &
          0.6760969921163639D+00, &
          0.0000000000000000D+00, &
         -0.6760969921163639D+00, &
          0.0000000000000000D+00, &
          0.6163321778509252D+00, &
          0.0000000000000000D+00, &
         -0.6163321778509252D+00, &
          0.0000000000000000D+00, &
          0.4319711692851803D+00, &
          0.0000000000000000D+00, &
         -0.4319711692851803D+00, &
          0.2418018136307699D+00, &
         -0.2418018136307699D+00, &
          0.2418018136307699D+00, &
         -0.2418018136307699D+00, &
          0.8440791629031899D+00, &
         -0.8440791629031899D+00, &
          0.8440791629031899D+00, &
         -0.8440791629031899D+00, &
          0.2279669387608253D+00, &
         -0.2279669387608253D+00, &
          0.2279669387608253D+00, &
         -0.2279669387608253D+00, &
          0.5027849771350101D+00, &
         -0.5027849771350101D+00, &
          0.5027849771350101D+00, &
         -0.5027849771350101D+00, &
          0.2605010749834311D+00, &
         -0.2605010749834311D+00, &
          0.2605010749834311D+00, &
         -0.2605010749834311D+00, &
          0.0926958730867243D+00, &
         -0.0926958730867243D+00, &
          0.0926958730867243D+00, &
         -0.0926958730867243D+00, &
          0.4832161680706445D+00, &
         -0.4832161680706445D+00, &
          0.4832161680706445D+00, &
         -0.4832161680706445D+00, &
          0.5671855056082324D+00, &
         -0.5671855056082324D+00, &
          0.5671855056082324D+00, &
         -0.5671855056082324D+00, &
          0.8463949915138761D+00, &
         -0.8463949915138761D+00, &
          0.8463949915138761D+00, &
         -0.8463949915138761D+00, &
          0.4641907129661463D+00, &
          0.8471778681177453D+00, &
         -0.4641907129661463D+00, &
          0.8471778681177453D+00, &
          0.4641907129661463D+00, &
         -0.8471778681177453D+00, &
         -0.4641907129661463D+00, &
         -0.8471778681177453D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.0371459309315055D+00, &
          0.6942002631703261D+00, &
          0.0000000042251285D+00, &
          0.0000000042251285D+00, &
          0.0000000042251285D+00, &
          0.0000000042251285D+00, &
          0.2697317845200571D+00, &
          0.2697317845200571D+00, &
          0.2697317845200571D+00, &
          0.2697317845200571D+00, &
          0.0833333333333333D+00, &
          0.0833333333333333D+00, &
          0.0833333333333333D+00, &
          0.0833333333333333D+00, &
          0.5334178104457834D+00, &
          0.5334178104457834D+00, &
          0.5334178104457834D+00, &
          0.5334178104457834D+00, &
          0.6619880087153507D+00, &
          0.6619880087153507D+00, &
          0.6619880087153507D+00, &
          0.6619880087153507D+00, &
          0.1522047098054636D+00, &
          0.1522047098054636D+00, &
          0.1522047098054636D+00, &
          0.1522047098054636D+00, &
          0.4222655365894455D+00, &
          0.4222655365894455D+00, &
          0.4222655365894455D+00, &
          0.4222655365894455D+00, &
          0.3967643698988902D+00, &
          0.3967643698988902D+00, &
          0.3967643698988902D+00, &
          0.3967643698988902D+00, &
          0.2008207219371707D+00, &
          0.2008207219371707D+00, &
          0.2008207219371707D+00, &
          0.2008207219371707D+00, &
          0.8661453707796378D+00, &
          0.8661453707796378D+00, &
          0.8661453707796378D+00, &
          0.8661453707796378D+00, &
          0.0195040124924115D+00, &
          0.0195040124924115D+00, &
          0.0195040124924115D+00, &
          0.0195040124924115D+00, &
          0.2131686017065601D+00, &
          0.2131686017065601D+00, &
          0.2131686017065601D+00, &
          0.2131686017065601D+00, &
          0.0153477972480097D+00, &
          0.0153477972480097D+00, &
          0.0153477972480097D+00, &
          0.0153477972480097D+00, &
          0.0833333333333333D+00, &
          0.0833333333333333D+00, &
          0.0833333333333333D+00, &
          0.0833333333333333D+00, &
          0.0833333333333333D+00, &
          0.0833333333333333D+00, &
          0.0833333333333333D+00, &
          0.0833333333333333D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          0.0320478691353425D+00, &
          0.0237932209450781D+00, &
          0.0041733723029944D+00, &
          0.0041733723029944D+00, &
          0.0041733723029944D+00, &
          0.0041733723029944D+00, &
          0.0222937913596180D+00, &
          0.0222937913596180D+00, &
          0.0222937913596180D+00, &
          0.0222937913596180D+00, &
          0.0291280227779532D+00, &
          0.0291280227779532D+00, &
          0.0291280227779532D+00, &
          0.0291280227779532D+00, &
          0.0115418065483010D+00, &
          0.0115418065483010D+00, &
          0.0115418065483010D+00, &
          0.0115418065483010D+00, &
          0.0101005748273597D+00, &
          0.0101005748273597D+00, &
          0.0101005748273597D+00, &
          0.0101005748273597D+00, &
          0.0020589782014206D+00, &
          0.0020589782014206D+00, &
          0.0020589782014206D+00, &
          0.0020589782014206D+00, &
          0.0342456887657178D+00, &
          0.0342456887657178D+00, &
          0.0342456887657178D+00, &
          0.0342456887657178D+00, &
          0.0129094136418972D+00, &
          0.0129094136418972D+00, &
          0.0129094136418972D+00, &
          0.0129094136418972D+00, &
          0.0371122811478335D+00, &
          0.0371122811478335D+00, &
          0.0371122811478335D+00, &
          0.0371122811478335D+00, &
          0.0022186981787511D+00, &
          0.0022186981787511D+00, &
          0.0022186981787511D+00, &
          0.0022186981787511D+00, &
          0.0159562091361951D+00, &
          0.0159562091361951D+00, &
          0.0159562091361951D+00, &
          0.0159562091361951D+00, &
          0.0230673181273276D+00, &
          0.0230673181273276D+00, &
          0.0230673181273276D+00, &
          0.0230673181273276D+00, &
          0.0045994005726981D+00, &
          0.0045994005726981D+00, &
          0.0045994005726981D+00, &
          0.0045994005726981D+00, &
          0.0133170859459138D+00, &
          0.0133170859459138D+00, &
          0.0133170859459138D+00, &
          0.0133170859459138D+00, &
          0.0133170859459138D+00, &
          0.0133170859459138D+00, &
          0.0133170859459138D+00, &
          0.0133170859459138D+00 /)

  x = x_save
  y = y_save
  z = z_save
  w = w_save

  return
end
subroutine rule10 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule10() returns the pyramid quadrature rule of precision 10.
!
!  Discussion:
!
!    We suppose we are given a pyramid P with vertices
!      (-1,-1,0), (-1,+1,0), (+1,+1,0), (+1,-1,0), (0,0,1).
!    We call a rule with n points, returning coordinates
!    x, y, z, and weights w.  Then the integral I of f(x,y,z) over
!    P is approximated by Q as follows:
!
!    Q = volume(P) * sum ( 1 <= i <= n ) w(i) * f(x(i),y(i),z(i)) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 April 2023
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
!    real x(n), y(n), z(n): the coordinates of quadrature points.
!
!    real w(n): the quadrature weights, which sum to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: n_save = 83

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
          0.0000000000000000D+00, &
          0.0000000000000000D+00, &
          0.0000000000000000D+00, &
          0.9211301560146403D+00, &
          0.0000000000000000D+00, &
         -0.9211301560146403D+00, &
          0.0000000000000000D+00, &
          0.4546948844252537D+00, &
          0.0000000000000000D+00, &
         -0.4546948844252537D+00, &
          0.0000000000000000D+00, &
          0.3495469844967067D+00, &
          0.0000000000000000D+00, &
         -0.3495469844967067D+00, &
          0.0000000000000000D+00, &
          0.4036909359989291D+00, &
          0.4036909359989291D+00, &
         -0.4036909359989291D+00, &
         -0.4036909359989291D+00, &
          0.3714608909064742D+00, &
          0.3714608909064742D+00, &
         -0.3714608909064742D+00, &
         -0.3714608909064742D+00, &
          0.2995979966434122D+00, &
          0.2995979966434122D+00, &
         -0.2995979966434122D+00, &
         -0.2995979966434122D+00, &
          0.1582262545541785D+00, &
          0.1582262545541785D+00, &
         -0.1582262545541785D+00, &
         -0.1582262545541785D+00, &
          0.1819979349519966D+00, &
          0.1819979349519966D+00, &
         -0.1819979349519966D+00, &
         -0.1819979349519966D+00, &
          0.8455841498012413D+00, &
          0.8455841498012413D+00, &
         -0.8455841498012413D+00, &
         -0.8455841498012413D+00, &
          0.6468694059373429D+00, &
          0.6468694059373429D+00, &
         -0.6468694059373429D+00, &
         -0.6468694059373429D+00, &
          0.5875870294087208D+00, &
          0.5875870294087208D+00, &
         -0.5875870294087208D+00, &
         -0.5875870294087208D+00, &
          0.2237620599838161D+00, &
          0.2237620599838161D+00, &
         -0.2237620599838161D+00, &
         -0.2237620599838161D+00, &
          0.7455657084697210D+00, &
          0.2247174318849062D+00, &
          0.7455657084697210D+00, &
         -0.2247174318849062D+00, &
         -0.7455657084697210D+00, &
          0.2247174318849062D+00, &
         -0.7455657084697210D+00, &
         -0.2247174318849062D+00, &
          0.4056619527318771D+00, &
          0.7487267628018788D+00, &
          0.4056619527318771D+00, &
         -0.7487267628018788D+00, &
         -0.4056619527318771D+00, &
          0.7487267628018788D+00, &
         -0.4056619527318771D+00, &
         -0.7487267628018788D+00, &
          0.0707425682812096D+00, &
          0.6004569111268736D+00, &
          0.0707425682812096D+00, &
         -0.6004569111268736D+00, &
         -0.0707425682812096D+00, &
          0.6004569111268736D+00, &
         -0.0707425682812096D+00, &
         -0.6004569111268736D+00, &
          0.9422191470796681D+00, &
          0.6375631040785387D+00, &
          0.9422191470796681D+00, &
         -0.6375631040785387D+00, &
         -0.9422191470796681D+00, &
          0.6375631040785387D+00, &
         -0.9422191470796681D+00, &
         -0.6375631040785387D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
          0.0000000000000000D+00, &
          0.0000000000000000D+00, &
          0.0000000000000000D+00, &
          0.0000000000000000D+00, &
          0.9211301560146403D+00, &
          0.0000000000000000D+00, &
         -0.9211301560146403D+00, &
          0.0000000000000000D+00, &
          0.4546948844252537D+00, &
          0.0000000000000000D+00, &
         -0.4546948844252537D+00, &
          0.0000000000000000D+00, &
          0.3495469844967067D+00, &
          0.0000000000000000D+00, &
         -0.3495469844967067D+00, &
          0.4036909359989291D+00, &
         -0.4036909359989291D+00, &
          0.4036909359989291D+00, &
         -0.4036909359989291D+00, &
          0.3714608909064742D+00, &
         -0.3714608909064742D+00, &
          0.3714608909064742D+00, &
         -0.3714608909064742D+00, &
          0.2995979966434122D+00, &
         -0.2995979966434122D+00, &
          0.2995979966434122D+00, &
         -0.2995979966434122D+00, &
          0.1582262545541785D+00, &
         -0.1582262545541785D+00, &
          0.1582262545541785D+00, &
         -0.1582262545541785D+00, &
          0.1819979349519966D+00, &
         -0.1819979349519966D+00, &
          0.1819979349519966D+00, &
         -0.1819979349519966D+00, &
          0.8455841498012413D+00, &
         -0.8455841498012413D+00, &
          0.8455841498012413D+00, &
         -0.8455841498012413D+00, &
          0.6468694059373429D+00, &
         -0.6468694059373429D+00, &
          0.6468694059373429D+00, &
         -0.6468694059373429D+00, &
          0.5875870294087208D+00, &
         -0.5875870294087208D+00, &
          0.5875870294087208D+00, &
         -0.5875870294087208D+00, &
          0.2237620599838161D+00, &
         -0.2237620599838161D+00, &
          0.2237620599838161D+00, &
         -0.2237620599838161D+00, &
          0.2247174318849062D+00, &
          0.7455657084697210D+00, &
         -0.2247174318849062D+00, &
          0.7455657084697210D+00, &
          0.2247174318849062D+00, &
         -0.7455657084697210D+00, &
         -0.2247174318849062D+00, &
         -0.7455657084697210D+00, &
          0.7487267628018788D+00, &
          0.4056619527318771D+00, &
         -0.7487267628018788D+00, &
          0.4056619527318771D+00, &
          0.7487267628018788D+00, &
         -0.4056619527318771D+00, &
         -0.7487267628018788D+00, &
         -0.4056619527318771D+00, &
          0.6004569111268736D+00, &
          0.0707425682812096D+00, &
         -0.6004569111268736D+00, &
          0.0707425682812096D+00, &
          0.6004569111268736D+00, &
         -0.0707425682812096D+00, &
         -0.6004569111268736D+00, &
         -0.0707425682812096D+00, &
          0.6375631040785387D+00, &
          0.9422191470796681D+00, &
         -0.6375631040785387D+00, &
          0.9422191470796681D+00, &
          0.6375631040785387D+00, &
         -0.9422191470796681D+00, &
         -0.6375631040785387D+00, &
         -0.9422191470796681D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
          0.9167957817791272D+00, &
          0.3511368063403526D+00, &
          0.7419135679633000D+00, &
          0.0631615972799145D+00, &
          0.0631615972799145D+00, &
          0.0631615972799145D+00, &
          0.0631615972799145D+00, &
          0.1767304111617324D+00, &
          0.1767304111617324D+00, &
          0.1767304111617324D+00, &
          0.1767304111617324D+00, &
          0.6020938107079140D+00, &
          0.6020938107079140D+00, &
          0.6020938107079140D+00, &
          0.6020938107079140D+00, &
          0.5270494697681531D+00, &
          0.5270494697681531D+00, &
          0.5270494697681531D+00, &
          0.5270494697681531D+00, &
          0.3491104418985491D+00, &
          0.3491104418985491D+00, &
          0.3491104418985491D+00, &
          0.3491104418985491D+00, &
          0.0069339728754634D+00, &
          0.0069339728754634D+00, &
          0.0069339728754634D+00, &
          0.0069339728754634D+00, &
          0.7752643179331966D+00, &
          0.7752643179331966D+00, &
          0.7752643179331966D+00, &
          0.7752643179331966D+00, &
          0.5470650181463391D+00, &
          0.5470650181463391D+00, &
          0.5470650181463391D+00, &
          0.5470650181463391D+00, &
          0.0728984114756605D+00, &
          0.0728984114756605D+00, &
          0.0728984114756605D+00, &
          0.0728984114756605D+00, &
          0.2836223397977548D+00, &
          0.2836223397977548D+00, &
          0.2836223397977548D+00, &
          0.2836223397977548D+00, &
          0.0598130217927265D+00, &
          0.0598130217927265D+00, &
          0.0598130217927265D+00, &
          0.0598130217927265D+00, &
          0.0658308799806233D+00, &
          0.0658308799806233D+00, &
          0.0658308799806233D+00, &
          0.0658308799806233D+00, &
          0.0272502552746740D+00, &
          0.0272502552746740D+00, &
          0.0272502552746740D+00, &
          0.0272502552746740D+00, &
          0.0272502552746740D+00, &
          0.0272502552746740D+00, &
          0.0272502552746740D+00, &
          0.0272502552746740D+00, &
          0.1618000129270508D+00, &
          0.1618000129270508D+00, &
          0.1618000129270508D+00, &
          0.1618000129270508D+00, &
          0.1618000129270508D+00, &
          0.1618000129270508D+00, &
          0.1618000129270508D+00, &
          0.1618000129270508D+00, &
          0.3584705635890507D+00, &
          0.3584705635890507D+00, &
          0.3584705635890507D+00, &
          0.3584705635890507D+00, &
          0.3584705635890507D+00, &
          0.3584705635890507D+00, &
          0.3584705635890507D+00, &
          0.3584705635890507D+00, &
          0.0074406410768847D+00, &
          0.0074406410768847D+00, &
          0.0074406410768847D+00, &
          0.0074406410768847D+00, &
          0.0074406410768847D+00, &
          0.0074406410768847D+00, &
          0.0074406410768847D+00, &
          0.0074406410768847D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
          0.0024012267625214D+00, &
          0.0469159279097637D+00, &
          0.0098466607485292D+00, &
          0.0056773923179313D+00, &
          0.0056773923179313D+00, &
          0.0056773923179313D+00, &
          0.0056773923179313D+00, &
          0.0412212597418143D+00, &
          0.0412212597418143D+00, &
          0.0412212597418143D+00, &
          0.0412212597418143D+00, &
          0.0092087532680317D+00, &
          0.0092087532680317D+00, &
          0.0092087532680317D+00, &
          0.0092087532680317D+00, &
          0.0070683961650691D+00, &
          0.0070683961650691D+00, &
          0.0070683961650691D+00, &
          0.0070683961650691D+00, &
          0.0298306858263750D+00, &
          0.0298306858263750D+00, &
          0.0298306858263750D+00, &
          0.0298306858263750D+00, &
          0.0058699928018168D+00, &
          0.0058699928018168D+00, &
          0.0058699928018168D+00, &
          0.0058699928018168D+00, &
          0.0049543654347689D+00, &
          0.0049543654347689D+00, &
          0.0049543654347689D+00, &
          0.0049543654347689D+00, &
          0.0159309372716684D+00, &
          0.0159309372716684D+00, &
          0.0159309372716684D+00, &
          0.0159309372716684D+00, &
          0.0050915274851829D+00, &
          0.0050915274851829D+00, &
          0.0050915274851829D+00, &
          0.0050915274851829D+00, &
          0.0068457933453088D+00, &
          0.0068457933453088D+00, &
          0.0068457933453088D+00, &
          0.0068457933453088D+00, &
          0.0150380203803664D+00, &
          0.0150380203803664D+00, &
          0.0150380203803664D+00, &
          0.0150380203803664D+00, &
          0.0156151631237942D+00, &
          0.0156151631237942D+00, &
          0.0156151631237942D+00, &
          0.0156151631237942D+00, &
          0.0089093600668980D+00, &
          0.0089093600668980D+00, &
          0.0089093600668980D+00, &
          0.0089093600668980D+00, &
          0.0089093600668980D+00, &
          0.0089093600668980D+00, &
          0.0089093600668980D+00, &
          0.0089093600668980D+00, &
          0.0167917214598521D+00, &
          0.0167917214598521D+00, &
          0.0167917214598521D+00, &
          0.0167917214598521D+00, &
          0.0167917214598521D+00, &
          0.0167917214598521D+00, &
          0.0167917214598521D+00, &
          0.0167917214598521D+00, &
          0.0082057491053038D+00, &
          0.0082057491053038D+00, &
          0.0082057491053038D+00, &
          0.0082057491053038D+00, &
          0.0082057491053038D+00, &
          0.0082057491053038D+00, &
          0.0082057491053038D+00, &
          0.0082057491053038D+00, &
          0.0025215488592804D+00, &
          0.0025215488592804D+00, &
          0.0025215488592804D+00, &
          0.0025215488592804D+00, &
          0.0025215488592804D+00, &
          0.0025215488592804D+00, &
          0.0025215488592804D+00, &
          0.0025215488592804D+00 /)

  x = x_save
  y = y_save
  z = z_save
  w = w_save

  return
end

