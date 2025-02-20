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
subroutine pyramid_jaskowiec_rule ( p, n, x, y, z, w )

!*****************************************************************************80
!
!! pyramid_jaskowiec_rule() returns a pyramid quadrature rule of given precision.
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
    write ( *, '(a)' ) 'pyramid_jaskowiec_rule(): Fatal error!'
    write ( *, '(a)' ) '  Input p < 0.'
    stop ( 1 )
  end if

  if ( 20 < p ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'pyramid_jaskowiec_rule(): Fatal error!'
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
!    15 April 2023
!
!  Author:
!
!    John Burkardt
!
!  Output:
!
!    real ( kind = rk ) pyramid_unit_volume: the volume of the pyramid.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) pyramid_unit_volume
  real ( kind = rk ) volume

  volume = 4.0D+00 / 3.0D+00

  pyramid_unit_volume = volume

  return
end
subroutine pyramid_volume ( h, s, volume )

!*****************************************************************************80
!
!! pyramid_volume() computes the volume of a pyramid with square base in 3D.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 April 2023
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) h, s: the height of the pyramid, and the 
!    length of one side of the square base.
!
!  Output:
!
!    real ( kind = rk ) volume: the volume of the pyramid.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) h
  real ( kind = rk ) s
  real ( kind = rk ) volume

  volume = s * s * h / 3.0D+00

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

!****************************************************************************80
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
!    15 April 2023
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
  integer, save, dimension(0:20) :: order_save = (/ &
      1, &
      1,   5,   6,  10,  15, &
     23,  31,  47,  62,  80, &
    103, 127, 152, 184, 234, &
    285, 319, 357, 418, 489 /)
  integer p

  if ( p < 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'rule_order(): Fatal error!'
    write ( *, '(a)' ) '  Input p < 0.'
    stop ( 1 )
  end if

  if ( 20 < p ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'rule_order(): Fatal error!'
    write ( *, '(a)' ) '  Input 20 < p.'
    stop ( 1 )
  end if

  order = order_save(p)

  return
end

subroutine rule00 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule00() returns the pyramid quadrature rule of precision 0.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 April 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 121, Number 11, pages 2418-2436, 15 June 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real ( kind = rk ) x[n], y[n], z[n]: the coordinates of
!    quadrature points.
!
!    real ( kind = rk ) w[n]: the quadrature weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 1

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
      0.00000000000000000000D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
      0.00000000000000000000D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
      0.25000000000000000000D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
      1.00000000000000000000D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule01 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule01() returns the pyramid quadrature rule of precision 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 April 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 121, Number 11, pages 2418-2436, 15 June 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real ( kind = rk ) x[n], y[n], z[n]: the coordinates of
!    quadrature points.
!
!    real ( kind = rk ) w[n]: the quadrature weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 1

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
      0.00000000000000000000D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
      0.00000000000000000000D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
      0.25000000000000000000D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
      1.00000000000000000000D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule02 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule02() returns the pyramid quadrature rule of precision 2.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 April 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 121, Number 11, pages 2418-2436, 15 June 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real ( kind = rk ) x[n], y[n], z[n]: the coordinates of
!    quadrature points.
!
!    real ( kind = rk ) w[n]: the quadrature weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 5

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
      0.00000000000000000000D+00, &
      0.52699748736717488828D+00, &
     -0.52699748736717488828D+00, &
      0.52699748736717488828D+00, &
     -0.52699748736717488828D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
      0.00000000000000000000D+00, &
      0.52699748736717488828D+00, &
      0.52699748736717488828D+00, &
     -0.52699748736717488828D+00, &
     -0.52699748736717488828D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
      0.56063221253561712487D+00, &
      0.12927845700902559911D+00, &
      0.12927845700902559911D+00, &
      0.12927845700902559911D+00, &
      0.12927845700902559911D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
      0.27986667890163369199D+00, &
      0.18003333027459159088D+00, &
      0.18003333027459159088D+00, &
      0.18003333027459159088D+00, &
      0.18003333027459159088D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule03 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule03() returns the pyramid quadrature rule of precision 3.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 April 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 121, Number 11, pages 2418-2436, 15 June 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real ( kind = rk ) x[n], y[n], z[n]: the coordinates of
!    quadrature points.
!
!    real ( kind = rk ) w[n]: the quadrature weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 6

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.58459636639471157515D+00, &
     -0.58459636639471157515D+00, &
      0.58459636639471157515D+00, &
     -0.58459636639471157515D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.58459636639471157515D+00, &
      0.58459636639471157515D+00, &
     -0.58459636639471157515D+00, &
     -0.58459636639471157515D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
      0.03032132711145601317D+00, &
      0.56560718797897435728D+00, &
      0.16666666666666665741D+00, &
      0.16666666666666665741D+00, &
      0.16666666666666665741D+00, &
      0.16666666666666665741D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
      0.15345064748545933497D+00, &
      0.26133122207480513621D+00, &
      0.14630453260993386833D+00, &
      0.14630453260993386833D+00, &
      0.14630453260993386833D+00, &
      0.14630453260993386833D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule04 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule04() returns the pyramid quadrature rule of precision 4.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 April 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 121, Number 11, pages 2418-2436, 15 June 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real ( kind = rk ) x[n], y[n], z[n]: the coordinates of
!    quadrature points.
!
!    real ( kind = rk ) w[n]: the quadrature weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 10

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.65058155639823256333D+00, &
      0.00000000000000000000D+00, &
     -0.65058155639823256333D+00, &
      0.00000000000000000000D+00, &
      0.65796699712169004481D+00, &
     -0.65796699712169004481D+00, &
      0.65796699712169004481D+00, &
     -0.65796699712169004481D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.65058155639823256333D+00, &
      0.00000000000000000000D+00, &
     -0.65058155639823256333D+00, &
      0.65796699712169004481D+00, &
      0.65796699712169004481D+00, &
     -0.65796699712169004481D+00, &
     -0.65796699712169004481D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
      0.12513695310874645150D+00, &
      0.67723278888613736015D+00, &
      0.32238414957821365237D+00, &
      0.32238414957821365237D+00, &
      0.32238414957821365237D+00, &
      0.32238414957821365237D+00, &
      0.03924828389881535040D+00, &
      0.03924828389881535040D+00, &
      0.03924828389881535040D+00, &
      0.03924828389881535040D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
      0.20688340258955226214D+00, &
      0.11374188317064193310D+00, &
      0.10632458788932550031D+00, &
      0.10632458788932550031D+00, &
      0.10632458788932550031D+00, &
      0.10632458788932550031D+00, &
      0.06351909067062594394D+00, &
      0.06351909067062594394D+00, &
      0.06351909067062594394D+00, &
      0.06351909067062594394D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule05 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule05() returns the pyramid quadrature rule of precision 5.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 April 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 121, Number 11, pages 2418-2436, 15 June 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real ( kind = rk ) x[n], y[n], z[n]: the coordinates of
!    quadrature points.
!
!    real ( kind = rk ) w[n]: the quadrature weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 15

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.75344061307932941318D+00, &
      0.00000000000000000000D+00, &
     -0.75344061307932941318D+00, &
      0.00000000000000000000D+00, &
      0.41715200242575134482D+00, &
     -0.41715200242575134482D+00, &
      0.41715200242575134482D+00, &
     -0.41715200242575134482D+00, &
      0.67402251647787037037D+00, &
     -0.67402251647787037037D+00, &
      0.67402251647787037037D+00, &
     -0.67402251647787037037D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.75344061307932941318D+00, &
      0.00000000000000000000D+00, &
     -0.75344061307932941318D+00, &
      0.41715200242575134482D+00, &
      0.41715200242575134482D+00, &
     -0.41715200242575134482D+00, &
     -0.41715200242575134482D+00, &
      0.67402251647787037037D+00, &
      0.67402251647787037037D+00, &
     -0.67402251647787037037D+00, &
     -0.67402251647787037037D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
      0.73070946955479043616D+00, &
      0.00619723285819058847D+00, &
      0.26844580953431373960D+00, &
      0.12500000000000000000D+00, &
      0.12500000000000000000D+00, &
      0.12500000000000000000D+00, &
      0.12500000000000000000D+00, &
      0.42182171100285947851D+00, &
      0.42182171100285947851D+00, &
      0.42182171100285947851D+00, &
      0.42182171100285947851D+00, &
      0.06579572180745926757D+00, &
      0.06579572180745926757D+00, &
      0.06579572180745926757D+00, &
      0.06579572180745926757D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
      0.06773442693037112772D+00, &
      0.06470893518150579171D+00, &
      0.17727154901514516339D+00, &
      0.05910777216655192096D+00, &
      0.05910777216655192096D+00, &
      0.05910777216655192096D+00, &
      0.05910777216655192096D+00, &
      0.06537546219121122271D+00, &
      0.06537546219121122271D+00, &
      0.06537546219121122271D+00, &
      0.06537546219121122271D+00, &
      0.04808803786048133910D+00, &
      0.04808803786048133910D+00, &
      0.04808803786048133910D+00, &
      0.04808803786048133910D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule06 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule06() returns the pyramid quadrature rule of precision 6.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 April 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 121, Number 11, pages 2418-2436, 15 June 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real ( kind = rk ) x[n], y[n], z[n]: the coordinates of
!    quadrature points.
!
!    real ( kind = rk ) w[n]: the quadrature weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 23

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.42104595182782333929D+00, &
      0.00000000000000000000D+00, &
     -0.42104595182782333929D+00, &
      0.00000000000000000000D+00, &
      0.83584092506524387822D+00, &
      0.00000000000000000000D+00, &
     -0.83584092506524387822D+00, &
      0.00000000000000000000D+00, &
      0.51341781341302172859D+00, &
     -0.51341781341302172859D+00, &
      0.51341781341302172859D+00, &
     -0.51341781341302172859D+00, &
      0.87197953364266822529D+00, &
     -0.87197953364266822529D+00, &
      0.87197953364266822529D+00, &
     -0.87197953364266822529D+00, &
      0.47733155776773072976D+00, &
     -0.47733155776773072976D+00, &
      0.47733155776773072976D+00, &
     -0.47733155776773072976D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.42104595182782333929D+00, &
      0.00000000000000000000D+00, &
     -0.42104595182782333929D+00, &
      0.00000000000000000000D+00, &
      0.83584092506524387822D+00, &
      0.00000000000000000000D+00, &
     -0.83584092506524387822D+00, &
      0.51341781341302172859D+00, &
      0.51341781341302172859D+00, &
     -0.51341781341302172859D+00, &
     -0.51341781341302172859D+00, &
      0.87197953364266822529D+00, &
      0.87197953364266822529D+00, &
     -0.87197953364266822529D+00, &
     -0.87197953364266822529D+00, &
      0.47733155776773072976D+00, &
      0.47733155776773072976D+00, &
     -0.47733155776773072976D+00, &
     -0.47733155776773072976D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
      0.13353121706321477435D+00, &
      0.80839181878746035892D+00, &
      0.37840352066355309457D+00, &
      0.55635774022808082151D+00, &
      0.55635774022808082151D+00, &
      0.55635774022808082151D+00, &
      0.55635774022808082151D+00, &
      0.09682668434012106640D+00, &
      0.09682668434012106640D+00, &
      0.09682668434012106640D+00, &
      0.09682668434012106640D+00, &
      0.25547807503740499468D+00, &
      0.25547807503740499468D+00, &
      0.25547807503740499468D+00, &
      0.25547807503740499468D+00, &
      0.03348911098405843445D+00, &
      0.03348911098405843445D+00, &
      0.03348911098405843445D+00, &
      0.03348911098405843445D+00, &
      0.02776222122928558023D+00, &
      0.02776222122928558023D+00, &
      0.02776222122928558023D+00, &
      0.02776222122928558023D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
      0.10236994192337051102D+00, &
      0.02544552509057920395D+00, &
      0.10744358342269332007D+00, &
      0.03715744178992643615D+00, &
      0.03715744178992643615D+00, &
      0.03715744178992643615D+00, &
      0.03715744178992643615D+00, &
      0.03663269740345383857D+00, &
      0.03663269740345383857D+00, &
      0.03663269740345383857D+00, &
      0.03663269740345383857D+00, &
      0.07134885171305939411D+00, &
      0.07134885171305939411D+00, &
      0.07134885171305939411D+00, &
      0.07134885171305939411D+00, &
      0.00865946139444006419D+00, &
      0.00865946139444006419D+00, &
      0.00865946139444006419D+00, &
      0.00865946139444006419D+00, &
      0.03738678508995950389D+00, &
      0.03738678508995950389D+00, &
      0.03738678508995950389D+00, &
      0.03738678508995950389D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule07 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule07() returns the pyramid quadrature rule of precision 7.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 April 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 121, Number 11, pages 2418-2436, 15 June 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real ( kind = rk ) x[n], y[n], z[n]: the coordinates of
!    quadrature points.
!
!    real ( kind = rk ) w[n]: the quadrature weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 31

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.61721339984836764980D+00, &
      0.00000000000000000000D+00, &
     -0.61721339984836764980D+00, &
      0.00000000000000000000D+00, &
      0.86409875978771466531D+00, &
      0.00000000000000000000D+00, &
     -0.86409875978771466531D+00, &
      0.00000000000000000000D+00, &
      0.52488756030374572603D+00, &
     -0.52488756030374572603D+00, &
      0.52488756030374572603D+00, &
     -0.52488756030374572603D+00, &
      0.25419682219463812789D+00, &
     -0.25419682219463812789D+00, &
      0.25419682219463812789D+00, &
     -0.25419682219463812789D+00, &
      0.35405111881016937403D+00, &
     -0.35405111881016937403D+00, &
      0.35405111881016937403D+00, &
     -0.35405111881016937403D+00, &
      0.61427194545119712110D+00, &
     -0.61427194545119712110D+00, &
      0.61427194545119712110D+00, &
     -0.61427194545119712110D+00, &
      0.80282248626994900942D+00, &
     -0.80282248626994900942D+00, &
      0.80282248626994900942D+00, &
     -0.80282248626994900942D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.61721339984836764980D+00, &
      0.00000000000000000000D+00, &
     -0.61721339984836764980D+00, &
      0.00000000000000000000D+00, &
      0.86409875978771466531D+00, &
      0.00000000000000000000D+00, &
     -0.86409875978771466531D+00, &
      0.52488756030374572603D+00, &
      0.52488756030374572603D+00, &
     -0.52488756030374572603D+00, &
     -0.52488756030374572603D+00, &
      0.25419682219463812789D+00, &
      0.25419682219463812789D+00, &
     -0.25419682219463812789D+00, &
     -0.25419682219463812789D+00, &
      0.35405111881016937403D+00, &
      0.35405111881016937403D+00, &
     -0.35405111881016937403D+00, &
     -0.35405111881016937403D+00, &
      0.61427194545119712110D+00, &
      0.61427194545119712110D+00, &
     -0.61427194545119712110D+00, &
     -0.61427194545119712110D+00, &
      0.80282248626994900942D+00, &
      0.80282248626994900942D+00, &
     -0.80282248626994900942D+00, &
     -0.80282248626994900942D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
      0.39365048525928414413D+00, &
      0.83863414272299030561D+00, &
      0.00001985131073852604D+00, &
      0.33333333333333331483D+00, &
      0.33333333333333331483D+00, &
      0.33333333333333331483D+00, &
      0.33333333333333331483D+00, &
      0.06666666666666666574D+00, &
      0.06666666666666666574D+00, &
      0.06666666666666666574D+00, &
      0.06666666666666666574D+00, &
      0.29045491084254104752D+00, &
      0.29045491084254104752D+00, &
      0.29045491084254104752D+00, &
      0.29045491084254104752D+00, &
      0.60547835568141594731D+00, &
      0.60547835568141594731D+00, &
      0.60547835568141594731D+00, &
      0.60547835568141594731D+00, &
      0.12931884631055998169D+00, &
      0.12931884631055998169D+00, &
      0.12931884631055998169D+00, &
      0.12931884631055998169D+00, &
      0.00010086339268113571D+00, &
      0.00010086339268113571D+00, &
      0.00010086339268113571D+00, &
      0.00010086339268113571D+00, &
      0.08012951317750569014D+00, &
      0.08012951317750569014D+00, &
      0.08012951317750569014D+00, &
      0.08012951317750569014D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
      0.10051398177493840735D+00, &
      0.01571901760701541889D+00, &
      0.02499658963028165634D+00, &
      0.02871093749999999861D+00, &
      0.02871093749999999861D+00, &
      0.02871093749999999861D+00, &
      0.02871093749999999861D+00, &
      0.02669175929300291600D+00, &
      0.02669175929300291600D+00, &
      0.02669175929300291600D+00, &
      0.02669175929300291600D+00, &
      0.03572750182264943647D+00, &
      0.03572750182264943647D+00, &
      0.03572750182264943647D+00, &
      0.03572750182264943647D+00, &
      0.02951904528668866309D+00, &
      0.02951904528668866309D+00, &
      0.02951904528668866309D+00, &
      0.02951904528668866309D+00, &
      0.06594160872648228977D+00, &
      0.06594160872648228977D+00, &
      0.06594160872648228977D+00, &
      0.06594160872648228977D+00, &
      0.01333274388639105884D+00, &
      0.01333274388639105884D+00, &
      0.01333274388639105884D+00, &
      0.01333274388639105884D+00, &
      0.01476900623172677438D+00, &
      0.01476900623172677438D+00, &
      0.01476900623172677438D+00, &
      0.01476900623172677438D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule08 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule08() returns the pyramid quadrature rule of precision 8.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 April 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 121, Number 11, pages 2418-2436, 15 June 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real ( kind = rk ) x[n], y[n], z[n]: the coordinates of
!    quadrature points.
!
!    real ( kind = rk ) w[n]: the quadrature weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 47

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.25367856157821822016D+00, &
      0.00000000000000000000D+00, &
     -0.25367856157821822016D+00, &
      0.00000000000000000000D+00, &
      0.71027375776907275551D+00, &
      0.00000000000000000000D+00, &
     -0.71027375776907275551D+00, &
      0.00000000000000000000D+00, &
      0.63643362359838895337D+00, &
      0.00000000000000000000D+00, &
     -0.63643362359838895337D+00, &
      0.00000000000000000000D+00, &
      0.62337928196226433109D+00, &
      0.00000000000000000000D+00, &
     -0.62337928196226433109D+00, &
      0.00000000000000000000D+00, &
      0.51225968171005897833D+00, &
     -0.51225968171005897833D+00, &
      0.51225968171005897833D+00, &
     -0.51225968171005897833D+00, &
      0.37965901379429650708D+00, &
     -0.37965901379429650708D+00, &
      0.37965901379429650708D+00, &
     -0.37965901379429650708D+00, &
      0.69140086940529510429D+00, &
     -0.69140086940529510429D+00, &
      0.69140086940529510429D+00, &
     -0.69140086940529510429D+00, &
      0.86742190798549856368D+00, &
     -0.86742190798549856368D+00, &
      0.86742190798549856368D+00, &
     -0.86742190798549856368D+00, &
      0.31712322629106232119D+00, &
     -0.31712322629106232119D+00, &
      0.31712322629106232119D+00, &
     -0.31712322629106232119D+00, &
      0.89374797161835639603D+00, &
      0.40528326346564669258D+00, &
     -0.89374797161835639603D+00, &
      0.40528326346564669258D+00, &
      0.89374797161835639603D+00, &
     -0.40528326346564669258D+00, &
     -0.89374797161835639603D+00, &
     -0.40528326346564669258D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.25367856157821822016D+00, &
      0.00000000000000000000D+00, &
     -0.25367856157821822016D+00, &
      0.00000000000000000000D+00, &
      0.71027375776907275551D+00, &
      0.00000000000000000000D+00, &
     -0.71027375776907275551D+00, &
      0.00000000000000000000D+00, &
      0.63643362359838895337D+00, &
      0.00000000000000000000D+00, &
     -0.63643362359838895337D+00, &
      0.00000000000000000000D+00, &
      0.62337928196226433109D+00, &
      0.00000000000000000000D+00, &
     -0.62337928196226433109D+00, &
      0.51225968171005897833D+00, &
      0.51225968171005897833D+00, &
     -0.51225968171005897833D+00, &
     -0.51225968171005897833D+00, &
      0.37965901379429650708D+00, &
      0.37965901379429650708D+00, &
     -0.37965901379429650708D+00, &
     -0.37965901379429650708D+00, &
      0.69140086940529510429D+00, &
      0.69140086940529510429D+00, &
     -0.69140086940529510429D+00, &
     -0.69140086940529510429D+00, &
      0.86742190798549856368D+00, &
      0.86742190798549856368D+00, &
     -0.86742190798549856368D+00, &
     -0.86742190798549856368D+00, &
      0.31712322629106232119D+00, &
      0.31712322629106232119D+00, &
     -0.31712322629106232119D+00, &
     -0.31712322629106232119D+00, &
      0.40528326346564669258D+00, &
      0.89374797161835639603D+00, &
      0.40528326346564669258D+00, &
     -0.89374797161835639603D+00, &
     -0.40528326346564669258D+00, &
      0.89374797161835639603D+00, &
     -0.40528326346564669258D+00, &
     -0.89374797161835639603D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
      0.07395194949759914538D+00, &
      0.48064180778578036168D+00, &
      0.89787700126494018882D+00, &
      0.70399494392200201442D+00, &
      0.70399494392200201442D+00, &
      0.70399494392200201442D+00, &
      0.70399494392200201442D+00, &
      0.15999762291015326432D+00, &
      0.15999762291015326432D+00, &
      0.15999762291015326432D+00, &
      0.15999762291015326432D+00, &
      0.34740846408160180880D+00, &
      0.34740846408160180880D+00, &
      0.34740846408160180880D+00, &
      0.34740846408160180880D+00, &
      0.01127682420195143774D+00, &
      0.01127682420195143774D+00, &
      0.01127682420195143774D+00, &
      0.01127682420195143774D+00, &
      0.06351022006373874262D+00, &
      0.06351022006373874262D+00, &
      0.06351022006373874262D+00, &
      0.06351022006373874262D+00, &
      0.46284226887006968409D+00, &
      0.46284226887006968409D+00, &
      0.46284226887006968409D+00, &
      0.46284226887006968409D+00, &
      0.19177130509938980496D+00, &
      0.19177130509938980496D+00, &
      0.19177130509938980496D+00, &
      0.19177130509938980496D+00, &
      0.01630913438364359896D+00, &
      0.01630913438364359896D+00, &
      0.01630913438364359896D+00, &
      0.01630913438364359896D+00, &
      0.23681987030130632887D+00, &
      0.23681987030130632887D+00, &
      0.23681987030130632887D+00, &
      0.23681987030130632887D+00, &
      0.05005997974535449785D+00, &
      0.05005997974535449785D+00, &
      0.05005997974535449785D+00, &
      0.05005997974535449785D+00, &
      0.05005997974535449785D+00, &
      0.05005997974535449785D+00, &
      0.05005997974535449785D+00, &
      0.05005997974535449785D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
      0.05595524252285600381D+00, &
      0.06694668391641565852D+00, &
      0.00475252306983395441D+00, &
      0.01488586682102476834D+00, &
      0.01488586682102476834D+00, &
      0.01488586682102476834D+00, &
      0.01488586682102476834D+00, &
      0.02455698624881830563D+00, &
      0.02455698624881830563D+00, &
      0.02455698624881830563D+00, &
      0.02455698624881830563D+00, &
      0.01608138988371909592D+00, &
      0.01608138988371909592D+00, &
      0.01608138988371909592D+00, &
      0.01608138988371909592D+00, &
      0.01442915622061930955D+00, &
      0.01442915622061930955D+00, &
      0.01442915622061930955D+00, &
      0.01442915622061930955D+00, &
      0.02488836268558412140D+00, &
      0.02488836268558412140D+00, &
      0.02488836268558412140D+00, &
      0.02488836268558412140D+00, &
      0.03016542061786949003D+00, &
      0.03016542061786949003D+00, &
      0.03016542061786949003D+00, &
      0.03016542061786949003D+00, &
      0.01825943823062004326D+00, &
      0.01825943823062004326D+00, &
      0.01825943823062004326D+00, &
      0.01825943823062004326D+00, &
      0.00405270511408486935D+00, &
      0.00405270511408486935D+00, &
      0.00405270511408486935D+00, &
      0.00405270511408486935D+00, &
      0.05109969513593878160D+00, &
      0.05109969513593878160D+00, &
      0.05109969513593878160D+00, &
      0.05109969513593878160D+00, &
      0.00983368333222240688D+00, &
      0.00983368333222240688D+00, &
      0.00983368333222240688D+00, &
      0.00983368333222240688D+00, &
      0.00983368333222240688D+00, &
      0.00983368333222240688D+00, &
      0.00983368333222240688D+00, &
      0.00983368333222240688D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule09 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule09() returns the pyramid quadrature rule of precision 9.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 April 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 121, Number 11, pages 2418-2436, 15 June 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real ( kind = rk ) x[n], y[n], z[n]: the coordinates of
!    quadrature points.
!
!    real ( kind = rk ) w[n]: the quadrature weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 62

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.21993584526338003093D+00, &
      0.00000000000000000000D+00, &
     -0.21993584526338003093D+00, &
      0.00000000000000000000D+00, &
      0.25580155737930071469D+00, &
      0.00000000000000000000D+00, &
     -0.25580155737930071469D+00, &
      0.00000000000000000000D+00, &
      0.87392422284165005575D+00, &
      0.00000000000000000000D+00, &
     -0.87392422284165005575D+00, &
      0.00000000000000000000D+00, &
      0.62069706482035569284D+00, &
      0.00000000000000000000D+00, &
     -0.62069706482035569284D+00, &
      0.00000000000000000000D+00, &
      0.48088720239801563405D+00, &
      0.00000000000000000000D+00, &
     -0.48088720239801563405D+00, &
      0.00000000000000000000D+00, &
      0.54766337559892563913D+00, &
      0.00000000000000000000D+00, &
     -0.54766337559892563913D+00, &
      0.00000000000000000000D+00, &
      0.31716782120661740629D+00, &
     -0.31716782120661740629D+00, &
      0.31716782120661740629D+00, &
     -0.31716782120661740629D+00, &
      0.58114214633277661015D+00, &
     -0.58114214633277661015D+00, &
      0.58114214633277661015D+00, &
     -0.58114214633277661015D+00, &
      0.52308134436113074006D+00, &
     -0.52308134436113074006D+00, &
      0.52308134436113074006D+00, &
     -0.52308134436113074006D+00, &
      0.87278863855859856180D+00, &
     -0.87278863855859856180D+00, &
      0.87278863855859856180D+00, &
     -0.87278863855859856180D+00, &
      0.27447830261314615230D+00, &
     -0.27447830261314615230D+00, &
      0.27447830261314615230D+00, &
     -0.27447830261314615230D+00, &
      0.45787401929460586070D+00, &
      0.76874559686548982196D+00, &
     -0.45787401929460586070D+00, &
      0.76874559686548982196D+00, &
      0.45787401929460586070D+00, &
     -0.76874559686548982196D+00, &
     -0.45787401929460586070D+00, &
     -0.76874559686548982196D+00, &
      0.52945920866195028687D+00, &
      0.86516329503292388470D+00, &
     -0.52945920866195028687D+00, &
      0.86516329503292388470D+00, &
      0.52945920866195028687D+00, &
     -0.86516329503292388470D+00, &
     -0.52945920866195028687D+00, &
     -0.86516329503292388470D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.21993584526338003093D+00, &
      0.00000000000000000000D+00, &
     -0.21993584526338003093D+00, &
      0.00000000000000000000D+00, &
      0.25580155737930071469D+00, &
      0.00000000000000000000D+00, &
     -0.25580155737930071469D+00, &
      0.00000000000000000000D+00, &
      0.87392422284165005575D+00, &
      0.00000000000000000000D+00, &
     -0.87392422284165005575D+00, &
      0.00000000000000000000D+00, &
      0.62069706482035569284D+00, &
      0.00000000000000000000D+00, &
     -0.62069706482035569284D+00, &
      0.00000000000000000000D+00, &
      0.48088720239801563405D+00, &
      0.00000000000000000000D+00, &
     -0.48088720239801563405D+00, &
      0.00000000000000000000D+00, &
      0.54766337559892563913D+00, &
      0.00000000000000000000D+00, &
     -0.54766337559892563913D+00, &
      0.31716782120661740629D+00, &
      0.31716782120661740629D+00, &
     -0.31716782120661740629D+00, &
     -0.31716782120661740629D+00, &
      0.58114214633277661015D+00, &
      0.58114214633277661015D+00, &
     -0.58114214633277661015D+00, &
     -0.58114214633277661015D+00, &
      0.52308134436113074006D+00, &
      0.52308134436113074006D+00, &
     -0.52308134436113074006D+00, &
     -0.52308134436113074006D+00, &
      0.87278863855859856180D+00, &
      0.87278863855859856180D+00, &
     -0.87278863855859856180D+00, &
     -0.87278863855859856180D+00, &
      0.27447830261314615230D+00, &
      0.27447830261314615230D+00, &
     -0.27447830261314615230D+00, &
     -0.27447830261314615230D+00, &
      0.76874559686548982196D+00, &
      0.45787401929460586070D+00, &
      0.76874559686548982196D+00, &
     -0.45787401929460586070D+00, &
     -0.76874559686548982196D+00, &
      0.45787401929460586070D+00, &
     -0.76874559686548982196D+00, &
     -0.45787401929460586070D+00, &
      0.86516329503292388470D+00, &
      0.52945920866195028687D+00, &
      0.86516329503292388470D+00, &
     -0.52945920866195028687D+00, &
     -0.86516329503292388470D+00, &
      0.52945920866195028687D+00, &
     -0.86516329503292388470D+00, &
     -0.52945920866195028687D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
      0.90338060334285785746D+00, &
      0.56496444359597119966D+00, &
      0.74856124635194676298D+00, &
      0.74856124635194676298D+00, &
      0.74856124635194676298D+00, &
      0.74856124635194676298D+00, &
      0.13469976629555441283D+00, &
      0.13469976629555441283D+00, &
      0.13469976629555441283D+00, &
      0.13469976629555441283D+00, &
      0.05467616134251215843D+00, &
      0.05467616134251215843D+00, &
      0.05467616134251215843D+00, &
      0.05467616134251215843D+00, &
      0.18946269484050054510D+00, &
      0.18946269484050054510D+00, &
      0.18946269484050054510D+00, &
      0.18946269484050054510D+00, &
      0.02248423582708249449D+00, &
      0.02248423582708249449D+00, &
      0.02248423582708249449D+00, &
      0.02248423582708249449D+00, &
      0.41340822589276282617D+00, &
      0.41340822589276282617D+00, &
      0.41340822589276282617D+00, &
      0.41340822589276282617D+00, &
      0.55077233639120526387D+00, &
      0.55077233639120526387D+00, &
      0.55077233639120526387D+00, &
      0.55077233639120526387D+00, &
      0.31914147555912980581D+00, &
      0.31914147555912980581D+00, &
      0.31914147555912980581D+00, &
      0.31914147555912980581D+00, &
      0.08935238781868591607D+00, &
      0.08935238781868591607D+00, &
      0.08935238781868591607D+00, &
      0.08935238781868591607D+00, &
      0.06151583405729108001D+00, &
      0.06151583405729108001D+00, &
      0.06151583405729108001D+00, &
      0.06151583405729108001D+00, &
      0.32642198515966086569D+00, &
      0.32642198515966086569D+00, &
      0.32642198515966086569D+00, &
      0.32642198515966086569D+00, &
      0.17557242043675652665D+00, &
      0.17557242043675652665D+00, &
      0.17557242043675652665D+00, &
      0.17557242043675652665D+00, &
      0.17557242043675652665D+00, &
      0.17557242043675652665D+00, &
      0.17557242043675652665D+00, &
      0.17557242043675652665D+00, &
      0.01570840570495662947D+00, &
      0.01570840570495662947D+00, &
      0.01570840570495662947D+00, &
      0.01570840570495662947D+00, &
      0.01570840570495662947D+00, &
      0.01570840570495662947D+00, &
      0.01570840570495662947D+00, &
      0.01570840570495662947D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
      0.00351798245822128233D+00, &
      0.04216520928063208912D+00, &
      0.00852104325745563912D+00, &
      0.00852104325745563912D+00, &
      0.00852104325745563912D+00, &
      0.00852104325745563912D+00, &
      0.02087072229798967934D+00, &
      0.02087072229798967934D+00, &
      0.02087072229798967934D+00, &
      0.02087072229798967934D+00, &
      0.01073969585101612438D+00, &
      0.01073969585101612438D+00, &
      0.01073969585101612438D+00, &
      0.01073969585101612438D+00, &
      0.02521321637389669496D+00, &
      0.02521321637389669496D+00, &
      0.02521321637389669496D+00, &
      0.02521321637389669496D+00, &
      0.01961258547003659133D+00, &
      0.01961258547003659133D+00, &
      0.01961258547003659133D+00, &
      0.01961258547003659133D+00, &
      0.01643323197880764905D+00, &
      0.01643323197880764905D+00, &
      0.01643323197880764905D+00, &
      0.01643323197880764905D+00, &
      0.01830589291063876300D+00, &
      0.01830589291063876300D+00, &
      0.01830589291063876300D+00, &
      0.01830589291063876300D+00, &
      0.01130841811263377274D+00, &
      0.01130841811263377274D+00, &
      0.01130841811263377274D+00, &
      0.01130841811263377274D+00, &
      0.02507245299443831496D+00, &
      0.02507245299443831496D+00, &
      0.02507245299443831496D+00, &
      0.02507245299443831496D+00, &
      0.00441940990434765493D+00, &
      0.00441940990434765493D+00, &
      0.00441940990434765493D+00, &
      0.00441940990434765493D+00, &
      0.04065272607298718588D+00, &
      0.04065272607298718588D+00, &
      0.04065272607298718588D+00, &
      0.04065272607298718588D+00, &
      0.01219793270502279002D+00, &
      0.01219793270502279002D+00, &
      0.01219793270502279002D+00, &
      0.01219793270502279002D+00, &
      0.01219793270502279002D+00, &
      0.01219793270502279002D+00, &
      0.01219793270502279002D+00, &
      0.01219793270502279002D+00, &
      0.00651697071549649839D+00, &
      0.00651697071549649839D+00, &
      0.00651697071549649839D+00, &
      0.00651697071549649839D+00, &
      0.00651697071549649839D+00, &
      0.00651697071549649839D+00, &
      0.00651697071549649839D+00, &
      0.00651697071549649839D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule10 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule10() returns the pyramid quadrature rule of precision 10.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 April 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 121, Number 11, pages 2418-2436, 15 June 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real ( kind = rk ) x[n], y[n], z[n]: the coordinates of
!    quadrature points.
!
!    real ( kind = rk ) w[n]: the quadrature weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 80

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.65420170176119341043D+00, &
      0.00000000000000000000D+00, &
     -0.65420170176119341043D+00, &
      0.00000000000000000000D+00, &
      0.76633927991746164654D+00, &
      0.00000000000000000000D+00, &
     -0.76633927991746164654D+00, &
      0.00000000000000000000D+00, &
      0.80421056101815546757D+00, &
      0.00000000000000000000D+00, &
     -0.80421056101815546757D+00, &
      0.00000000000000000000D+00, &
      0.21719757883626783501D+00, &
      0.00000000000000000000D+00, &
     -0.21719757883626783501D+00, &
      0.00000000000000000000D+00, &
      0.20937347351036650345D+00, &
      0.00000000000000000000D+00, &
     -0.20937347351036650345D+00, &
      0.00000000000000000000D+00, &
      0.52830709053289259813D+00, &
      0.00000000000000000000D+00, &
     -0.52830709053289259813D+00, &
      0.00000000000000000000D+00, &
      0.34681490190497532566D+00, &
     -0.34681490190497532566D+00, &
      0.34681490190497532566D+00, &
     -0.34681490190497532566D+00, &
      0.54234774233239257946D+00, &
     -0.54234774233239257946D+00, &
      0.54234774233239257946D+00, &
     -0.54234774233239257946D+00, &
      0.42965769544556714488D+00, &
     -0.42965769544556714488D+00, &
      0.42965769544556714488D+00, &
     -0.42965769544556714488D+00, &
      0.71871414601724947779D+00, &
     -0.71871414601724947779D+00, &
      0.71871414601724947779D+00, &
     -0.71871414601724947779D+00, &
      0.76244008578748778682D+00, &
     -0.76244008578748778682D+00, &
      0.76244008578748778682D+00, &
     -0.76244008578748778682D+00, &
      0.96427350086319296718D+00, &
     -0.96427350086319296718D+00, &
      0.96427350086319296718D+00, &
     -0.96427350086319296718D+00, &
      0.27897388698743774693D+00, &
     -0.27897388698743774693D+00, &
      0.27897388698743774693D+00, &
     -0.27897388698743774693D+00, &
      0.73124529835235163588D+00, &
      0.37115067890003261564D+00, &
     -0.73124529835235163588D+00, &
      0.37115067890003261564D+00, &
      0.73124529835235163588D+00, &
     -0.37115067890003261564D+00, &
     -0.73124529835235163588D+00, &
     -0.37115067890003261564D+00, &
      0.92194807057971739361D+00, &
      0.41966522833393532510D+00, &
     -0.92194807057971739361D+00, &
      0.41966522833393532510D+00, &
      0.92194807057971739361D+00, &
     -0.41966522833393532510D+00, &
     -0.92194807057971739361D+00, &
     -0.41966522833393532510D+00, &
      0.19954261519157223681D+00, &
      0.44340097447242488027D+00, &
     -0.19954261519157223681D+00, &
      0.44340097447242488027D+00, &
      0.19954261519157223681D+00, &
     -0.44340097447242488027D+00, &
     -0.19954261519157223681D+00, &
     -0.44340097447242488027D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.65420170176119341043D+00, &
      0.00000000000000000000D+00, &
     -0.65420170176119341043D+00, &
      0.00000000000000000000D+00, &
      0.76633927991746164654D+00, &
      0.00000000000000000000D+00, &
     -0.76633927991746164654D+00, &
      0.00000000000000000000D+00, &
      0.80421056101815546757D+00, &
      0.00000000000000000000D+00, &
     -0.80421056101815546757D+00, &
      0.00000000000000000000D+00, &
      0.21719757883626783501D+00, &
      0.00000000000000000000D+00, &
     -0.21719757883626783501D+00, &
      0.00000000000000000000D+00, &
      0.20937347351036650345D+00, &
      0.00000000000000000000D+00, &
     -0.20937347351036650345D+00, &
      0.00000000000000000000D+00, &
      0.52830709053289259813D+00, &
      0.00000000000000000000D+00, &
     -0.52830709053289259813D+00, &
      0.34681490190497532566D+00, &
      0.34681490190497532566D+00, &
     -0.34681490190497532566D+00, &
     -0.34681490190497532566D+00, &
      0.54234774233239257946D+00, &
      0.54234774233239257946D+00, &
     -0.54234774233239257946D+00, &
     -0.54234774233239257946D+00, &
      0.42965769544556714488D+00, &
      0.42965769544556714488D+00, &
     -0.42965769544556714488D+00, &
     -0.42965769544556714488D+00, &
      0.71871414601724947779D+00, &
      0.71871414601724947779D+00, &
     -0.71871414601724947779D+00, &
     -0.71871414601724947779D+00, &
      0.76244008578748778682D+00, &
      0.76244008578748778682D+00, &
     -0.76244008578748778682D+00, &
     -0.76244008578748778682D+00, &
      0.96427350086319296718D+00, &
      0.96427350086319296718D+00, &
     -0.96427350086319296718D+00, &
     -0.96427350086319296718D+00, &
      0.27897388698743774693D+00, &
      0.27897388698743774693D+00, &
     -0.27897388698743774693D+00, &
     -0.27897388698743774693D+00, &
      0.37115067890003261564D+00, &
      0.73124529835235163588D+00, &
      0.37115067890003261564D+00, &
     -0.73124529835235163588D+00, &
     -0.37115067890003261564D+00, &
      0.73124529835235163588D+00, &
     -0.37115067890003261564D+00, &
     -0.73124529835235163588D+00, &
      0.41966522833393532510D+00, &
      0.92194807057971739361D+00, &
      0.41966522833393532510D+00, &
     -0.92194807057971739361D+00, &
     -0.41966522833393532510D+00, &
      0.92194807057971739361D+00, &
     -0.41966522833393532510D+00, &
     -0.92194807057971739361D+00, &
      0.44340097447242488027D+00, &
      0.19954261519157223681D+00, &
      0.44340097447242488027D+00, &
     -0.19954261519157223681D+00, &
     -0.44340097447242488027D+00, &
      0.19954261519157223681D+00, &
     -0.44340097447242488027D+00, &
     -0.19954261519157223681D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
      0.91766654896784238815D+00, &
      0.70900937404687269794D+00, &
      0.22713967401042231553D+00, &
      0.06874011069480644165D+00, &
      0.12392544997399050632D+00, &
      0.12392544997399050632D+00, &
      0.12392544997399050632D+00, &
      0.12392544997399050632D+00, &
      0.02096299791450871586D+00, &
      0.02096299791450871586D+00, &
      0.02096299791450871586D+00, &
      0.02096299791450871586D+00, &
      0.16577313546766450636D+00, &
      0.16577313546766450636D+00, &
      0.16577313546766450636D+00, &
      0.16577313546766450636D+00, &
      0.49340520933953563310D+00, &
      0.49340520933953563310D+00, &
      0.49340520933953563310D+00, &
      0.49340520933953563310D+00, &
      0.78969767134839341516D+00, &
      0.78969767134839341516D+00, &
      0.78969767134839341516D+00, &
      0.78969767134839341516D+00, &
      0.47169288904858786005D+00, &
      0.47169288904858786005D+00, &
      0.47169288904858786005D+00, &
      0.47169288904858786005D+00, &
      0.02011394567580290782D+00, &
      0.02011394567580290782D+00, &
      0.02011394567580290782D+00, &
      0.02011394567580290782D+00, &
      0.36822101998927786459D+00, &
      0.36822101998927786459D+00, &
      0.36822101998927786459D+00, &
      0.36822101998927786459D+00, &
      0.12661594867701117528D+00, &
      0.12661594867701117528D+00, &
      0.12661594867701117528D+00, &
      0.12661594867701117528D+00, &
      0.02433819982878694319D+00, &
      0.02433819982878694319D+00, &
      0.02433819982878694319D+00, &
      0.02433819982878694319D+00, &
      0.12701442502569038062D+00, &
      0.12701442502569038062D+00, &
      0.12701442502569038062D+00, &
      0.12701442502569038062D+00, &
      0.01067692573406054529D+00, &
      0.01067692573406054529D+00, &
      0.01067692573406054529D+00, &
      0.01067692573406054529D+00, &
      0.59704297313030385563D+00, &
      0.59704297313030385563D+00, &
      0.59704297313030385563D+00, &
      0.59704297313030385563D+00, &
      0.21756166315495037433D+00, &
      0.21756166315495037433D+00, &
      0.21756166315495037433D+00, &
      0.21756166315495037433D+00, &
      0.21756166315495037433D+00, &
      0.21756166315495037433D+00, &
      0.21756166315495037433D+00, &
      0.21756166315495037433D+00, &
      0.04137857207280886546D+00, &
      0.04137857207280886546D+00, &
      0.04137857207280886546D+00, &
      0.04137857207280886546D+00, &
      0.04137857207280886546D+00, &
      0.04137857207280886546D+00, &
      0.04137857207280886546D+00, &
      0.04137857207280886546D+00, &
      0.32615365723990619173D+00, &
      0.32615365723990619173D+00, &
      0.32615365723990619173D+00, &
      0.32615365723990619173D+00, &
      0.32615365723990619173D+00, &
      0.32615365723990619173D+00, &
      0.32615365723990619173D+00, &
      0.32615365723990619173D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
      0.00222688724675269151D+00, &
      0.01676862119517166794D+00, &
      0.05218524659101939772D+00, &
      0.02912684113372536118D+00, &
      0.02061250827144200243D+00, &
      0.02061250827144200243D+00, &
      0.02061250827144200243D+00, &
      0.02061250827144200243D+00, &
      0.01001318317045431922D+00, &
      0.01001318317045431922D+00, &
      0.01001318317045431922D+00, &
      0.01001318317045431922D+00, &
      0.00380589471978953851D+00, &
      0.00380589471978953851D+00, &
      0.00380589471978953851D+00, &
      0.00380589471978953851D+00, &
      0.01603356460939567643D+00, &
      0.01603356460939567643D+00, &
      0.01603356460939567643D+00, &
      0.01603356460939567643D+00, &
      0.00421057398582402181D+00, &
      0.00421057398582402181D+00, &
      0.00421057398582402181D+00, &
      0.00421057398582402181D+00, &
      0.00812973230941097316D+00, &
      0.00812973230941097316D+00, &
      0.00812973230941097316D+00, &
      0.00812973230941097316D+00, &
      0.01345646795638386767D+00, &
      0.01345646795638386767D+00, &
      0.01345646795638386767D+00, &
      0.01345646795638386767D+00, &
      0.01016376299548378027D+00, &
      0.01016376299548378027D+00, &
      0.01016376299548378027D+00, &
      0.01016376299548378027D+00, &
      0.02981645762492043245D+00, &
      0.02981645762492043245D+00, &
      0.02981645762492043245D+00, &
      0.02981645762492043245D+00, &
      0.00865712783423741410D+00, &
      0.00865712783423741410D+00, &
      0.00865712783423741410D+00, &
      0.00865712783423741410D+00, &
      0.00953199916493123987D+00, &
      0.00953199916493123987D+00, &
      0.00953199916493123987D+00, &
      0.00953199916493123987D+00, &
      0.00071395826747580517D+00, &
      0.00071395826747580517D+00, &
      0.00071395826747580517D+00, &
      0.00071395826747580517D+00, &
      0.01519293397619959039D+00, &
      0.01519293397619959039D+00, &
      0.01519293397619959039D+00, &
      0.01519293397619959039D+00, &
      0.01008533639437561233D+00, &
      0.01008533639437561233D+00, &
      0.01008533639437561233D+00, &
      0.01008533639437561233D+00, &
      0.01008533639437561233D+00, &
      0.01008533639437561233D+00, &
      0.01008533639437561233D+00, &
      0.01008533639437561233D+00, &
      0.00585933822996072967D+00, &
      0.00585933822996072967D+00, &
      0.00585933822996072967D+00, &
      0.00585933822996072967D+00, &
      0.00585933822996072967D+00, &
      0.00585933822996072967D+00, &
      0.00585933822996072967D+00, &
      0.00585933822996072967D+00, &
      0.02134779341185568530D+00, &
      0.02134779341185568530D+00, &
      0.02134779341185568530D+00, &
      0.02134779341185568530D+00, &
      0.02134779341185568530D+00, &
      0.02134779341185568530D+00, &
      0.02134779341185568530D+00, &
      0.02134779341185568530D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule11 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule11() returns the pyramid quadrature rule of precision 11.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 April 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 121, Number 11, pages 2418-2436, 15 June 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real ( kind = rk ) x[n], y[n], z[n]: the coordinates of
!    quadrature points.
!
!    real ( kind = rk ) w[n]: the quadrature weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 103

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.30959144111215058937D+00, &
      0.00000000000000000000D+00, &
     -0.30959144111215058937D+00, &
      0.00000000000000000000D+00, &
      0.76993922081029786408D+00, &
      0.00000000000000000000D+00, &
     -0.76993922081029786408D+00, &
      0.00000000000000000000D+00, &
      0.70027038117446127607D+00, &
      0.00000000000000000000D+00, &
     -0.70027038117446127607D+00, &
      0.00000000000000000000D+00, &
      0.46095035429072822586D+00, &
      0.00000000000000000000D+00, &
     -0.46095035429072822586D+00, &
      0.00000000000000000000D+00, &
      0.29928938216959510843D+00, &
      0.00000000000000000000D+00, &
     -0.29928938216959510843D+00, &
      0.00000000000000000000D+00, &
      0.59114943757155702375D+00, &
      0.00000000000000000000D+00, &
     -0.59114943757155702375D+00, &
      0.00000000000000000000D+00, &
      0.01861308147704909821D+00, &
      0.00000000000000000000D+00, &
     -0.01861308147704909821D+00, &
      0.00000000000000000000D+00, &
      0.16034290975464657314D+00, &
      0.00000000000000000000D+00, &
     -0.16034290975464657314D+00, &
      0.00000000000000000000D+00, &
      0.14450753554447465232D+00, &
     -0.14450753554447465232D+00, &
      0.14450753554447465232D+00, &
     -0.14450753554447465232D+00, &
      0.59369450979675508773D+00, &
     -0.59369450979675508773D+00, &
      0.59369450979675508773D+00, &
     -0.59369450979675508773D+00, &
      0.27490141802950085470D+00, &
     -0.27490141802950085470D+00, &
      0.27490141802950085470D+00, &
     -0.27490141802950085470D+00, &
      0.29680077247241393179D+00, &
     -0.29680077247241393179D+00, &
      0.29680077247241393179D+00, &
     -0.29680077247241393179D+00, &
      0.59849263706569033605D+00, &
     -0.59849263706569033605D+00, &
      0.59849263706569033605D+00, &
     -0.59849263706569033605D+00, &
      0.61031415307328662490D+00, &
     -0.61031415307328662490D+00, &
      0.61031415307328662490D+00, &
     -0.61031415307328662490D+00, &
      0.37345023081856359992D+00, &
     -0.37345023081856359992D+00, &
      0.37345023081856359992D+00, &
     -0.37345023081856359992D+00, &
      0.26245991463277401623D+00, &
     -0.26245991463277401623D+00, &
      0.26245991463277401623D+00, &
     -0.26245991463277401623D+00, &
      0.36860635719790568743D+00, &
     -0.36860635719790568743D+00, &
      0.36860635719790568743D+00, &
     -0.36860635719790568743D+00, &
      0.87026293442361135622D+00, &
     -0.87026293442361135622D+00, &
      0.87026293442361135622D+00, &
     -0.87026293442361135622D+00, &
      0.81025745002423155139D+00, &
     -0.81025745002423155139D+00, &
      0.81025745002423155139D+00, &
     -0.81025745002423155139D+00, &
      0.78967657467472973654D+00, &
      0.35527213835047594115D+00, &
     -0.78967657467472973654D+00, &
      0.35527213835047594115D+00, &
      0.78967657467472973654D+00, &
     -0.35527213835047594115D+00, &
     -0.78967657467472973654D+00, &
     -0.35527213835047594115D+00, &
      0.92805034944130604391D+00, &
      0.41484210234755791724D+00, &
     -0.92805034944130604391D+00, &
      0.41484210234755791724D+00, &
      0.92805034944130604391D+00, &
     -0.41484210234755791724D+00, &
     -0.92805034944130604391D+00, &
     -0.41484210234755791724D+00, &
      0.23313789083599262275D+00, &
      0.55628779191329114084D+00, &
     -0.23313789083599262275D+00, &
      0.55628779191329114084D+00, &
      0.23313789083599262275D+00, &
     -0.55628779191329114084D+00, &
     -0.23313789083599262275D+00, &
     -0.55628779191329114084D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.30959144111215058937D+00, &
      0.00000000000000000000D+00, &
     -0.30959144111215058937D+00, &
      0.00000000000000000000D+00, &
      0.76993922081029786408D+00, &
      0.00000000000000000000D+00, &
     -0.76993922081029786408D+00, &
      0.00000000000000000000D+00, &
      0.70027038117446127607D+00, &
      0.00000000000000000000D+00, &
     -0.70027038117446127607D+00, &
      0.00000000000000000000D+00, &
      0.46095035429072822586D+00, &
      0.00000000000000000000D+00, &
     -0.46095035429072822586D+00, &
      0.00000000000000000000D+00, &
      0.29928938216959510843D+00, &
      0.00000000000000000000D+00, &
     -0.29928938216959510843D+00, &
      0.00000000000000000000D+00, &
      0.59114943757155702375D+00, &
      0.00000000000000000000D+00, &
     -0.59114943757155702375D+00, &
      0.00000000000000000000D+00, &
      0.01861308147704909821D+00, &
      0.00000000000000000000D+00, &
     -0.01861308147704909821D+00, &
      0.00000000000000000000D+00, &
      0.16034290975464657314D+00, &
      0.00000000000000000000D+00, &
     -0.16034290975464657314D+00, &
      0.14450753554447465232D+00, &
      0.14450753554447465232D+00, &
     -0.14450753554447465232D+00, &
     -0.14450753554447465232D+00, &
      0.59369450979675508773D+00, &
      0.59369450979675508773D+00, &
     -0.59369450979675508773D+00, &
     -0.59369450979675508773D+00, &
      0.27490141802950085470D+00, &
      0.27490141802950085470D+00, &
     -0.27490141802950085470D+00, &
     -0.27490141802950085470D+00, &
      0.29680077247241393179D+00, &
      0.29680077247241393179D+00, &
     -0.29680077247241393179D+00, &
     -0.29680077247241393179D+00, &
      0.59849263706569033605D+00, &
      0.59849263706569033605D+00, &
     -0.59849263706569033605D+00, &
     -0.59849263706569033605D+00, &
      0.61031415307328662490D+00, &
      0.61031415307328662490D+00, &
     -0.61031415307328662490D+00, &
     -0.61031415307328662490D+00, &
      0.37345023081856359992D+00, &
      0.37345023081856359992D+00, &
     -0.37345023081856359992D+00, &
     -0.37345023081856359992D+00, &
      0.26245991463277401623D+00, &
      0.26245991463277401623D+00, &
     -0.26245991463277401623D+00, &
     -0.26245991463277401623D+00, &
      0.36860635719790568743D+00, &
      0.36860635719790568743D+00, &
     -0.36860635719790568743D+00, &
     -0.36860635719790568743D+00, &
      0.87026293442361135622D+00, &
      0.87026293442361135622D+00, &
     -0.87026293442361135622D+00, &
     -0.87026293442361135622D+00, &
      0.81025745002423155139D+00, &
      0.81025745002423155139D+00, &
     -0.81025745002423155139D+00, &
     -0.81025745002423155139D+00, &
      0.35527213835047594115D+00, &
      0.78967657467472973654D+00, &
      0.35527213835047594115D+00, &
     -0.78967657467472973654D+00, &
     -0.35527213835047594115D+00, &
      0.78967657467472973654D+00, &
     -0.35527213835047594115D+00, &
     -0.78967657467472973654D+00, &
      0.41484210234755791724D+00, &
      0.92805034944130604391D+00, &
      0.41484210234755791724D+00, &
     -0.92805034944130604391D+00, &
     -0.41484210234755791724D+00, &
      0.92805034944130604391D+00, &
     -0.41484210234755791724D+00, &
     -0.92805034944130604391D+00, &
      0.55628779191329114084D+00, &
      0.23313789083599262275D+00, &
      0.55628779191329114084D+00, &
     -0.23313789083599262275D+00, &
     -0.55628779191329114084D+00, &
      0.23313789083599262275D+00, &
     -0.55628779191329114084D+00, &
     -0.23313789083599262275D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
      0.69295639647508833203D+00, &
      0.91101523446275922691D+00, &
      0.25032929072548470995D+00, &
      0.05040798390564713710D+00, &
      0.05040798390564713710D+00, &
      0.05040798390564713710D+00, &
      0.05040798390564713710D+00, &
      0.01178348626907169128D+00, &
      0.01178348626907169128D+00, &
      0.01178348626907169128D+00, &
      0.01178348626907169128D+00, &
      0.10103940902556154957D+00, &
      0.10103940902556154957D+00, &
      0.10103940902556154957D+00, &
      0.10103940902556154957D+00, &
      0.34272757605804599068D+00, &
      0.34272757605804599068D+00, &
      0.34272757605804599068D+00, &
      0.34272757605804599068D+00, &
      0.62869898216776509692D+00, &
      0.62869898216776509692D+00, &
      0.62869898216776509692D+00, &
      0.62869898216776509692D+00, &
      0.25822848523085090156D+00, &
      0.25822848523085090156D+00, &
      0.25822848523085090156D+00, &
      0.25822848523085090156D+00, &
      0.43583449548328001555D+00, &
      0.43583449548328001555D+00, &
      0.43583449548328001555D+00, &
      0.43583449548328001555D+00, &
      0.21928319414297889334D+00, &
      0.21928319414297889334D+00, &
      0.21928319414297889334D+00, &
      0.21928319414297889334D+00, &
      0.78797433053468768360D+00, &
      0.78797433053468768360D+00, &
      0.78797433053468768360D+00, &
      0.78797433053468768360D+00, &
      0.30892958869107134401D+00, &
      0.30892958869107134401D+00, &
      0.30892958869107134401D+00, &
      0.30892958869107134401D+00, &
      0.01781185684646453826D+00, &
      0.01781185684646453826D+00, &
      0.01781185684646453826D+00, &
      0.01781185684646453826D+00, &
      0.11127033608982600521D+00, &
      0.11127033608982600521D+00, &
      0.11127033608982600521D+00, &
      0.11127033608982600521D+00, &
      0.02444716310077480956D+00, &
      0.02444716310077480956D+00, &
      0.02444716310077480956D+00, &
      0.02444716310077480956D+00, &
      0.12041875519309763742D+00, &
      0.12041875519309763742D+00, &
      0.12041875519309763742D+00, &
      0.12041875519309763742D+00, &
      0.57752532555241498091D+00, &
      0.57752532555241498091D+00, &
      0.57752532555241498091D+00, &
      0.57752532555241498091D+00, &
      0.49169801908816329616D+00, &
      0.49169801908816329616D+00, &
      0.49169801908816329616D+00, &
      0.49169801908816329616D+00, &
      0.26311325712772259955D+00, &
      0.26311325712772259955D+00, &
      0.26311325712772259955D+00, &
      0.26311325712772259955D+00, &
      0.01908681361607608359D+00, &
      0.01908681361607608359D+00, &
      0.01908681361607608359D+00, &
      0.01908681361607608359D+00, &
      0.11692619926762114202D+00, &
      0.11692619926762114202D+00, &
      0.11692619926762114202D+00, &
      0.11692619926762114202D+00, &
      0.17968095390177224457D+00, &
      0.17968095390177224457D+00, &
      0.17968095390177224457D+00, &
      0.17968095390177224457D+00, &
      0.17968095390177224457D+00, &
      0.17968095390177224457D+00, &
      0.17968095390177224457D+00, &
      0.17968095390177224457D+00, &
      0.03535559961728908934D+00, &
      0.03535559961728908934D+00, &
      0.03535559961728908934D+00, &
      0.03535559961728908934D+00, &
      0.03535559961728908934D+00, &
      0.03535559961728908934D+00, &
      0.03535559961728908934D+00, &
      0.03535559961728908934D+00, &
      0.41464402024224555898D+00, &
      0.41464402024224555898D+00, &
      0.41464402024224555898D+00, &
      0.41464402024224555898D+00, &
      0.41464402024224555898D+00, &
      0.41464402024224555898D+00, &
      0.41464402024224555898D+00, &
      0.41464402024224555898D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
      0.01286363653695535118D+00, &
      0.00259572122542773513D+00, &
      0.00817638570047552610D+00, &
      0.00511941439804289126D+00, &
      0.00511941439804289126D+00, &
      0.00511941439804289126D+00, &
      0.00511941439804289126D+00, &
      0.00684699666646234988D+00, &
      0.00684699666646234988D+00, &
      0.00684699666646234988D+00, &
      0.00684699666646234988D+00, &
      0.02021405278036775971D+00, &
      0.02021405278036775971D+00, &
      0.02021405278036775971D+00, &
      0.02021405278036775971D+00, &
      0.01240017518258633086D+00, &
      0.01240017518258633086D+00, &
      0.01240017518258633086D+00, &
      0.01240017518258633086D+00, &
      0.01050060859922616330D+00, &
      0.01050060859922616330D+00, &
      0.01050060859922616330D+00, &
      0.01050060859922616330D+00, &
      0.01207523991498536502D+00, &
      0.01207523991498536502D+00, &
      0.01207523991498536502D+00, &
      0.01207523991498536502D+00, &
      0.00794549439937168242D+00, &
      0.00794549439937168242D+00, &
      0.00794549439937168242D+00, &
      0.00794549439937168242D+00, &
      0.00897930801963979161D+00, &
      0.00897930801963979161D+00, &
      0.00897930801963979161D+00, &
      0.00897930801963979161D+00, &
      0.00432861007717266853D+00, &
      0.00432861007717266853D+00, &
      0.00432861007717266853D+00, &
      0.00432861007717266853D+00, &
      0.00959150593850333251D+00, &
      0.00959150593850333251D+00, &
      0.00959150593850333251D+00, &
      0.00959150593850333251D+00, &
      0.00845889607542490291D+00, &
      0.00845889607542490291D+00, &
      0.00845889607542490291D+00, &
      0.00845889607542490291D+00, &
      0.02093917618186586990D+00, &
      0.02093917618186586990D+00, &
      0.02093917618186586990D+00, &
      0.02093917618186586990D+00, &
      0.01007335618421271312D+00, &
      0.01007335618421271312D+00, &
      0.01007335618421271312D+00, &
      0.01007335618421271312D+00, &
      0.01474935950898293054D+00, &
      0.01474935950898293054D+00, &
      0.01474935950898293054D+00, &
      0.01474935950898293054D+00, &
      0.00405432020025909591D+00, &
      0.00405432020025909591D+00, &
      0.00405432020025909591D+00, &
      0.00405432020025909591D+00, &
      0.01768012056788626288D+00, &
      0.01768012056788626288D+00, &
      0.01768012056788626288D+00, &
      0.01768012056788626288D+00, &
      0.02349367280688881635D+00, &
      0.02349367280688881635D+00, &
      0.02349367280688881635D+00, &
      0.02349367280688881635D+00, &
      0.00300355345952720053D+00, &
      0.00300355345952720053D+00, &
      0.00300355345952720053D+00, &
      0.00300355345952720053D+00, &
      0.00433057971890565482D+00, &
      0.00433057971890565482D+00, &
      0.00433057971890565482D+00, &
      0.00433057971890565482D+00, &
      0.00788297764391489818D+00, &
      0.00788297764391489818D+00, &
      0.00788297764391489818D+00, &
      0.00788297764391489818D+00, &
      0.00788297764391489818D+00, &
      0.00788297764391489818D+00, &
      0.00788297764391489818D+00, &
      0.00788297764391489818D+00, &
      0.00509491185532656536D+00, &
      0.00509491185532656536D+00, &
      0.00509491185532656536D+00, &
      0.00509491185532656536D+00, &
      0.00509491185532656536D+00, &
      0.00509491185532656536D+00, &
      0.00509491185532656536D+00, &
      0.00509491185532656536D+00, &
      0.00667542222774531943D+00, &
      0.00667542222774531943D+00, &
      0.00667542222774531943D+00, &
      0.00667542222774531943D+00, &
      0.00667542222774531943D+00, &
      0.00667542222774531943D+00, &
      0.00667542222774531943D+00, &
      0.00667542222774531943D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule12 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule12() returns the pyramid quadrature rule of precision 12.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 April 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 121, Number 11, pages 2418-2436, 15 June 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real ( kind = rk ) x[n], y[n], z[n]: the coordinates of
!    quadrature points.
!
!    real ( kind = rk ) w[n]: the quadrature weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 127

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.61752408228797095457D+00, &
      0.00000000000000000000D+00, &
     -0.61752408228797095457D+00, &
      0.00000000000000000000D+00, &
      0.95383417813475923630D+00, &
      0.00000000000000000000D+00, &
     -0.95383417813475923630D+00, &
      0.00000000000000000000D+00, &
      0.59503005801385822071D+00, &
      0.00000000000000000000D+00, &
     -0.59503005801385822071D+00, &
      0.00000000000000000000D+00, &
      0.43539359820232903520D+00, &
      0.00000000000000000000D+00, &
     -0.43539359820232903520D+00, &
      0.00000000000000000000D+00, &
      0.13291355561208040292D+00, &
      0.00000000000000000000D+00, &
     -0.13291355561208040292D+00, &
      0.00000000000000000000D+00, &
      0.79714363482446959353D+00, &
      0.00000000000000000000D+00, &
     -0.79714363482446959353D+00, &
      0.00000000000000000000D+00, &
      0.31078295049399401462D+00, &
      0.00000000000000000000D+00, &
     -0.31078295049399401462D+00, &
      0.00000000000000000000D+00, &
      0.46896896703369211901D+00, &
      0.00000000000000000000D+00, &
     -0.46896896703369211901D+00, &
      0.00000000000000000000D+00, &
      0.06329314929065793516D+00, &
      0.00000000000000000000D+00, &
     -0.06329314929065793516D+00, &
      0.00000000000000000000D+00, &
      0.11694124473665239161D+00, &
     -0.11694124473665239161D+00, &
      0.11694124473665239161D+00, &
     -0.11694124473665239161D+00, &
      0.47710850655214198657D+00, &
     -0.47710850655214198657D+00, &
      0.47710850655214198657D+00, &
     -0.47710850655214198657D+00, &
      0.34796735400888462175D+00, &
     -0.34796735400888462175D+00, &
      0.34796735400888462175D+00, &
     -0.34796735400888462175D+00, &
      0.21748093989891267852D+00, &
     -0.21748093989891267852D+00, &
      0.21748093989891267852D+00, &
     -0.21748093989891267852D+00, &
      0.67426207541025851011D+00, &
     -0.67426207541025851011D+00, &
      0.67426207541025851011D+00, &
     -0.67426207541025851011D+00, &
      0.65415312577816908668D+00, &
     -0.65415312577816908668D+00, &
      0.65415312577816908668D+00, &
     -0.65415312577816908668D+00, &
      0.26043868995957136780D+00, &
     -0.26043868995957136780D+00, &
      0.26043868995957136780D+00, &
     -0.26043868995957136780D+00, &
      0.27442107868377230151D+00, &
     -0.27442107868377230151D+00, &
      0.27442107868377230151D+00, &
     -0.27442107868377230151D+00, &
      0.51259022737109982693D+00, &
     -0.51259022737109982693D+00, &
      0.51259022737109982693D+00, &
     -0.51259022737109982693D+00, &
      0.86878873539049084052D+00, &
     -0.86878873539049084052D+00, &
      0.86878873539049084052D+00, &
     -0.86878873539049084052D+00, &
      0.71181681800308360675D+00, &
     -0.71181681800308360675D+00, &
      0.71181681800308360675D+00, &
     -0.71181681800308360675D+00, &
      0.36754569532388259301D+00, &
     -0.36754569532388259301D+00, &
      0.36754569532388259301D+00, &
     -0.36754569532388259301D+00, &
      0.84896937482456835689D+00, &
      0.51203646317203110883D+00, &
     -0.84896937482456835689D+00, &
      0.51203646317203110883D+00, &
      0.84896937482456835689D+00, &
     -0.51203646317203110883D+00, &
     -0.84896937482456835689D+00, &
     -0.51203646317203110883D+00, &
      0.97385403135997028468D+00, &
      0.74670340940557777820D+00, &
     -0.97385403135997028468D+00, &
      0.74670340940557777820D+00, &
      0.97385403135997028468D+00, &
     -0.74670340940557777820D+00, &
     -0.97385403135997028468D+00, &
     -0.74670340940557777820D+00, &
      0.32302405518012644592D+00, &
      0.62555613253804132068D+00, &
     -0.32302405518012644592D+00, &
      0.62555613253804132068D+00, &
      0.32302405518012644592D+00, &
     -0.62555613253804132068D+00, &
     -0.32302405518012644592D+00, &
     -0.62555613253804132068D+00, &
      0.19450076345283720536D+00, &
      0.69805959364273262313D+00, &
     -0.19450076345283720536D+00, &
      0.69805959364273262313D+00, &
      0.19450076345283720536D+00, &
     -0.69805959364273262313D+00, &
     -0.19450076345283720536D+00, &
     -0.69805959364273262313D+00, &
      0.32894885387189143344D+00, &
      0.85383152896625225114D+00, &
     -0.32894885387189143344D+00, &
      0.85383152896625225114D+00, &
      0.32894885387189143344D+00, &
     -0.85383152896625225114D+00, &
     -0.32894885387189143344D+00, &
     -0.85383152896625225114D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.61752408228797095457D+00, &
      0.00000000000000000000D+00, &
     -0.61752408228797095457D+00, &
      0.00000000000000000000D+00, &
      0.95383417813475923630D+00, &
      0.00000000000000000000D+00, &
     -0.95383417813475923630D+00, &
      0.00000000000000000000D+00, &
      0.59503005801385822071D+00, &
      0.00000000000000000000D+00, &
     -0.59503005801385822071D+00, &
      0.00000000000000000000D+00, &
      0.43539359820232903520D+00, &
      0.00000000000000000000D+00, &
     -0.43539359820232903520D+00, &
      0.00000000000000000000D+00, &
      0.13291355561208040292D+00, &
      0.00000000000000000000D+00, &
     -0.13291355561208040292D+00, &
      0.00000000000000000000D+00, &
      0.79714363482446959353D+00, &
      0.00000000000000000000D+00, &
     -0.79714363482446959353D+00, &
      0.00000000000000000000D+00, &
      0.31078295049399401462D+00, &
      0.00000000000000000000D+00, &
     -0.31078295049399401462D+00, &
      0.00000000000000000000D+00, &
      0.46896896703369211901D+00, &
      0.00000000000000000000D+00, &
     -0.46896896703369211901D+00, &
      0.00000000000000000000D+00, &
      0.06329314929065793516D+00, &
      0.00000000000000000000D+00, &
     -0.06329314929065793516D+00, &
      0.11694124473665239161D+00, &
      0.11694124473665239161D+00, &
     -0.11694124473665239161D+00, &
     -0.11694124473665239161D+00, &
      0.47710850655214198657D+00, &
      0.47710850655214198657D+00, &
     -0.47710850655214198657D+00, &
     -0.47710850655214198657D+00, &
      0.34796735400888462175D+00, &
      0.34796735400888462175D+00, &
     -0.34796735400888462175D+00, &
     -0.34796735400888462175D+00, &
      0.21748093989891267852D+00, &
      0.21748093989891267852D+00, &
     -0.21748093989891267852D+00, &
     -0.21748093989891267852D+00, &
      0.67426207541025851011D+00, &
      0.67426207541025851011D+00, &
     -0.67426207541025851011D+00, &
     -0.67426207541025851011D+00, &
      0.65415312577816908668D+00, &
      0.65415312577816908668D+00, &
     -0.65415312577816908668D+00, &
     -0.65415312577816908668D+00, &
      0.26043868995957136780D+00, &
      0.26043868995957136780D+00, &
     -0.26043868995957136780D+00, &
     -0.26043868995957136780D+00, &
      0.27442107868377230151D+00, &
      0.27442107868377230151D+00, &
     -0.27442107868377230151D+00, &
     -0.27442107868377230151D+00, &
      0.51259022737109982693D+00, &
      0.51259022737109982693D+00, &
     -0.51259022737109982693D+00, &
     -0.51259022737109982693D+00, &
      0.86878873539049084052D+00, &
      0.86878873539049084052D+00, &
     -0.86878873539049084052D+00, &
     -0.86878873539049084052D+00, &
      0.71181681800308360675D+00, &
      0.71181681800308360675D+00, &
     -0.71181681800308360675D+00, &
     -0.71181681800308360675D+00, &
      0.36754569532388259301D+00, &
      0.36754569532388259301D+00, &
     -0.36754569532388259301D+00, &
     -0.36754569532388259301D+00, &
      0.51203646317203110883D+00, &
      0.84896937482456835689D+00, &
      0.51203646317203110883D+00, &
     -0.84896937482456835689D+00, &
     -0.51203646317203110883D+00, &
      0.84896937482456835689D+00, &
     -0.51203646317203110883D+00, &
     -0.84896937482456835689D+00, &
      0.74670340940557777820D+00, &
      0.97385403135997028468D+00, &
      0.74670340940557777820D+00, &
     -0.97385403135997028468D+00, &
     -0.74670340940557777820D+00, &
      0.97385403135997028468D+00, &
     -0.74670340940557777820D+00, &
     -0.97385403135997028468D+00, &
      0.62555613253804132068D+00, &
      0.32302405518012644592D+00, &
      0.62555613253804132068D+00, &
     -0.32302405518012644592D+00, &
     -0.62555613253804132068D+00, &
      0.32302405518012644592D+00, &
     -0.62555613253804132068D+00, &
     -0.32302405518012644592D+00, &
      0.69805959364273262313D+00, &
      0.19450076345283720536D+00, &
      0.69805959364273262313D+00, &
     -0.19450076345283720536D+00, &
     -0.69805959364273262313D+00, &
      0.19450076345283720536D+00, &
     -0.69805959364273262313D+00, &
     -0.19450076345283720536D+00, &
      0.85383152896625225114D+00, &
      0.32894885387189143344D+00, &
      0.85383152896625225114D+00, &
     -0.32894885387189143344D+00, &
     -0.85383152896625225114D+00, &
      0.32894885387189143344D+00, &
     -0.85383152896625225114D+00, &
     -0.32894885387189143344D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
      0.45700700424071905026D+00, &
      0.96248480104808953328D+00, &
      0.02656531684537740551D+00, &
      0.00537472082577391923D+00, &
      0.00537472082577391923D+00, &
      0.00537472082577391923D+00, &
      0.00537472082577391923D+00, &
      0.03942123626649530338D+00, &
      0.03942123626649530338D+00, &
      0.03942123626649530338D+00, &
      0.03942123626649530338D+00, &
      0.07619184421748394220D+00, &
      0.07619184421748394220D+00, &
      0.07619184421748394220D+00, &
      0.07619184421748394220D+00, &
      0.52461875533620183631D+00, &
      0.52461875533620183631D+00, &
      0.52461875533620183631D+00, &
      0.52461875533620183631D+00, &
      0.83871994031142604875D+00, &
      0.83871994031142604875D+00, &
      0.83871994031142604875D+00, &
      0.83871994031142604875D+00, &
      0.20075189777897944898D+00, &
      0.20075189777897944898D+00, &
      0.20075189777897944898D+00, &
      0.20075189777897944898D+00, &
      0.68265469105314702247D+00, &
      0.68265469105314702247D+00, &
      0.68265469105314702247D+00, &
      0.68265469105314702247D+00, &
      0.31713382467521444852D+00, &
      0.31713382467521444852D+00, &
      0.31713382467521444852D+00, &
      0.31713382467521444852D+00, &
      0.12709629034671870995D+00, &
      0.12709629034671870995D+00, &
      0.12709629034671870995D+00, &
      0.12709629034671870995D+00, &
      0.67124174294922867023D+00, &
      0.67124174294922867023D+00, &
      0.67124174294922867023D+00, &
      0.67124174294922867023D+00, &
      0.45862802718642248223D+00, &
      0.45862802718642248223D+00, &
      0.45862802718642248223D+00, &
      0.45862802718642248223D+00, &
      0.11938890613701966248D+00, &
      0.11938890613701966248D+00, &
      0.11938890613701966248D+00, &
      0.11938890613701966248D+00, &
      0.25741784866521960629D+00, &
      0.25741784866521960629D+00, &
      0.25741784866521960629D+00, &
      0.25741784866521960629D+00, &
      0.02260502558938975656D+00, &
      0.02260502558938975656D+00, &
      0.02260502558938975656D+00, &
      0.02260502558938975656D+00, &
      0.08746307755558260788D+00, &
      0.08746307755558260788D+00, &
      0.08746307755558260788D+00, &
      0.08746307755558260788D+00, &
      0.67140164127972612462D+00, &
      0.67140164127972612462D+00, &
      0.67140164127972612462D+00, &
      0.67140164127972612462D+00, &
      0.47604255048878973966D+00, &
      0.47604255048878973966D+00, &
      0.47604255048878973966D+00, &
      0.47604255048878973966D+00, &
      0.22777955815540271156D+00, &
      0.22777955815540271156D+00, &
      0.22777955815540271156D+00, &
      0.22777955815540271156D+00, &
      0.05850272176027389998D+00, &
      0.05850272176027389998D+00, &
      0.05850272176027389998D+00, &
      0.05850272176027389998D+00, &
      0.22053995731942055425D+00, &
      0.22053995731942055425D+00, &
      0.22053995731942055425D+00, &
      0.22053995731942055425D+00, &
      0.02548724290196442005D+00, &
      0.02548724290196442005D+00, &
      0.02548724290196442005D+00, &
      0.02548724290196442005D+00, &
      0.11145836537804799937D+00, &
      0.11145836537804799937D+00, &
      0.11145836537804799937D+00, &
      0.11145836537804799937D+00, &
      0.11145836537804799937D+00, &
      0.11145836537804799937D+00, &
      0.11145836537804799937D+00, &
      0.11145836537804799937D+00, &
      0.00339220824679912656D+00, &
      0.00339220824679912656D+00, &
      0.00339220824679912656D+00, &
      0.00339220824679912656D+00, &
      0.00339220824679912656D+00, &
      0.00339220824679912656D+00, &
      0.00339220824679912656D+00, &
      0.00339220824679912656D+00, &
      0.32725151644051669875D+00, &
      0.32725151644051669875D+00, &
      0.32725151644051669875D+00, &
      0.32725151644051669875D+00, &
      0.32725151644051669875D+00, &
      0.32725151644051669875D+00, &
      0.32725151644051669875D+00, &
      0.32725151644051669875D+00, &
      0.14805166861898322317D+00, &
      0.14805166861898322317D+00, &
      0.14805166861898322317D+00, &
      0.14805166861898322317D+00, &
      0.14805166861898322317D+00, &
      0.14805166861898322317D+00, &
      0.14805166861898322317D+00, &
      0.14805166861898322317D+00, &
      0.02487563023274973195D+00, &
      0.02487563023274973195D+00, &
      0.02487563023274973195D+00, &
      0.02487563023274973195D+00, &
      0.02487563023274973195D+00, &
      0.02487563023274973195D+00, &
      0.02487563023274973195D+00, &
      0.02487563023274973195D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
      0.02715854620871203245D+00, &
      0.00043531424244899530D+00, &
      0.01391222140872406220D+00, &
      0.00387288146649252939D+00, &
      0.00387288146649252939D+00, &
      0.00387288146649252939D+00, &
      0.00387288146649252939D+00, &
      0.00217967971593618394D+00, &
      0.00217967971593618394D+00, &
      0.00217967971593618394D+00, &
      0.00217967971593618394D+00, &
      0.01267133485574918096D+00, &
      0.01267133485574918096D+00, &
      0.01267133485574918096D+00, &
      0.01267133485574918096D+00, &
      0.00873464163035380693D+00, &
      0.00873464163035380693D+00, &
      0.00873464163035380693D+00, &
      0.00873464163035380693D+00, &
      0.00289146731615850605D+00, &
      0.00289146731615850605D+00, &
      0.00289146731615850605D+00, &
      0.00289146731615850605D+00, &
      0.00336686848511522736D+00, &
      0.00336686848511522736D+00, &
      0.00336686848511522736D+00, &
      0.00336686848511522736D+00, &
      0.00147603803134411875D+00, &
      0.00147603803134411875D+00, &
      0.00147603803134411875D+00, &
      0.00147603803134411875D+00, &
      0.01998444773599242219D+00, &
      0.01998444773599242219D+00, &
      0.01998444773599242219D+00, &
      0.01998444773599242219D+00, &
      0.00618927170644437228D+00, &
      0.00618927170644437228D+00, &
      0.00618927170644437228D+00, &
      0.00618927170644437228D+00, &
      0.00820078705478863301D+00, &
      0.00820078705478863301D+00, &
      0.00820078705478863301D+00, &
      0.00820078705478863301D+00, &
      0.00509651546458162342D+00, &
      0.00509651546458162342D+00, &
      0.00509651546458162342D+00, &
      0.00509651546458162342D+00, &
      0.01739298850473944280D+00, &
      0.01739298850473944280D+00, &
      0.01739298850473944280D+00, &
      0.01739298850473944280D+00, &
      0.01948797782027694020D+00, &
      0.01948797782027694020D+00, &
      0.01948797782027694020D+00, &
      0.01948797782027694020D+00, &
      0.00538896243467025700D+00, &
      0.00538896243467025700D+00, &
      0.00538896243467025700D+00, &
      0.00538896243467025700D+00, &
      0.00676911731451347500D+00, &
      0.00676911731451347500D+00, &
      0.00676911731451347500D+00, &
      0.00676911731451347500D+00, &
      0.00467587609366508332D+00, &
      0.00467587609366508332D+00, &
      0.00467587609366508332D+00, &
      0.00467587609366508332D+00, &
      0.01838172037009152757D+00, &
      0.01838172037009152757D+00, &
      0.01838172037009152757D+00, &
      0.01838172037009152757D+00, &
      0.01407964506457294943D+00, &
      0.01407964506457294943D+00, &
      0.01407964506457294943D+00, &
      0.01407964506457294943D+00, &
      0.00283370063460703292D+00, &
      0.00283370063460703292D+00, &
      0.00283370063460703292D+00, &
      0.00283370063460703292D+00, &
      0.00479137222817773486D+00, &
      0.00479137222817773486D+00, &
      0.00479137222817773486D+00, &
      0.00479137222817773486D+00, &
      0.00990680633819146332D+00, &
      0.00990680633819146332D+00, &
      0.00990680633819146332D+00, &
      0.00990680633819146332D+00, &
      0.00535574363023301383D+00, &
      0.00535574363023301383D+00, &
      0.00535574363023301383D+00, &
      0.00535574363023301383D+00, &
      0.00535574363023301383D+00, &
      0.00535574363023301383D+00, &
      0.00535574363023301383D+00, &
      0.00535574363023301383D+00, &
      0.00081681851071822213D+00, &
      0.00081681851071822213D+00, &
      0.00081681851071822213D+00, &
      0.00081681851071822213D+00, &
      0.00081681851071822213D+00, &
      0.00081681851071822213D+00, &
      0.00081681851071822213D+00, &
      0.00081681851071822213D+00, &
      0.00854175372189276694D+00, &
      0.00854175372189276694D+00, &
      0.00854175372189276694D+00, &
      0.00854175372189276694D+00, &
      0.00854175372189276694D+00, &
      0.00854175372189276694D+00, &
      0.00854175372189276694D+00, &
      0.00854175372189276694D+00, &
      0.01018769584471248754D+00, &
      0.01018769584471248754D+00, &
      0.01018769584471248754D+00, &
      0.01018769584471248754D+00, &
      0.01018769584471248754D+00, &
      0.01018769584471248754D+00, &
      0.01018769584471248754D+00, &
      0.01018769584471248754D+00, &
      0.00572367792672661812D+00, &
      0.00572367792672661812D+00, &
      0.00572367792672661812D+00, &
      0.00572367792672661812D+00, &
      0.00572367792672661812D+00, &
      0.00572367792672661812D+00, &
      0.00572367792672661812D+00, &
      0.00572367792672661812D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule13 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule13() returns the pyramid quadrature rule of precision 13.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 April 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 121, Number 11, pages 2418-2436, 15 June 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real ( kind = rk ) x[n], y[n], z[n]: the coordinates of
!    quadrature points.
!
!    real ( kind = rk ) w[n]: the quadrature weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 152

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.41223892088328983885D+00, &
      0.00000000000000000000D+00, &
     -0.41223892088328983885D+00, &
      0.00000000000000000000D+00, &
      0.94482408797958883362D+00, &
      0.00000000000000000000D+00, &
     -0.94482408797958883362D+00, &
      0.00000000000000000000D+00, &
      0.15012306550295120844D+00, &
      0.00000000000000000000D+00, &
     -0.15012306550295120844D+00, &
      0.00000000000000000000D+00, &
      0.36338283577562219273D+00, &
      0.00000000000000000000D+00, &
     -0.36338283577562219273D+00, &
      0.00000000000000000000D+00, &
      0.32897122373450704558D+00, &
      0.00000000000000000000D+00, &
     -0.32897122373450704558D+00, &
      0.00000000000000000000D+00, &
      0.55704935891608620135D+00, &
      0.00000000000000000000D+00, &
     -0.55704935891608620135D+00, &
      0.00000000000000000000D+00, &
      0.71060354437421247020D+00, &
      0.00000000000000000000D+00, &
     -0.71060354437421247020D+00, &
      0.00000000000000000000D+00, &
      0.88531748708681257121D+00, &
      0.00000000000000000000D+00, &
     -0.88531748708681257121D+00, &
      0.00000000000000000000D+00, &
      0.48884754503841615358D+00, &
     -0.48884754503841615358D+00, &
      0.48884754503841615358D+00, &
     -0.48884754503841615358D+00, &
      0.45031772456672747307D+00, &
     -0.45031772456672747307D+00, &
      0.45031772456672747307D+00, &
     -0.45031772456672747307D+00, &
      0.65872956915547864476D+00, &
     -0.65872956915547864476D+00, &
      0.65872956915547864476D+00, &
     -0.65872956915547864476D+00, &
      0.62588271911798842861D+00, &
     -0.62588271911798842861D+00, &
      0.62588271911798842861D+00, &
     -0.62588271911798842861D+00, &
      0.25178960033193620305D+00, &
     -0.25178960033193620305D+00, &
      0.25178960033193620305D+00, &
     -0.25178960033193620305D+00, &
      0.91338191280449132492D+00, &
     -0.91338191280449132492D+00, &
      0.91338191280449132492D+00, &
     -0.91338191280449132492D+00, &
      0.11117861641720541699D+00, &
     -0.11117861641720541699D+00, &
      0.11117861641720541699D+00, &
     -0.11117861641720541699D+00, &
      0.37814924987375336807D+00, &
     -0.37814924987375336807D+00, &
      0.37814924987375336807D+00, &
     -0.37814924987375336807D+00, &
      0.23146195691577758913D+00, &
     -0.23146195691577758913D+00, &
      0.23146195691577758913D+00, &
     -0.23146195691577758913D+00, &
      0.18598667282576661353D+00, &
     -0.18598667282576661353D+00, &
      0.18598667282576661353D+00, &
     -0.18598667282576661353D+00, &
      0.72941490542892961635D+00, &
     -0.72941490542892961635D+00, &
      0.72941490542892961635D+00, &
     -0.72941490542892961635D+00, &
      0.80419241777926819825D+00, &
     -0.80419241777926819825D+00, &
      0.80419241777926819825D+00, &
     -0.80419241777926819825D+00, &
      0.47036082838417470064D+00, &
     -0.47036082838417470064D+00, &
      0.47036082838417470064D+00, &
     -0.47036082838417470064D+00, &
      0.27535657520011547206D+00, &
     -0.27535657520011547206D+00, &
      0.27535657520011547206D+00, &
     -0.27535657520011547206D+00, &
      0.25840443662578793660D+00, &
     -0.25840443662578793660D+00, &
      0.25840443662578793660D+00, &
     -0.25840443662578793660D+00, &
      0.77456301686419115615D+00, &
      0.43978938971942055369D+00, &
     -0.77456301686419115615D+00, &
      0.43978938971942055369D+00, &
      0.77456301686419115615D+00, &
     -0.43978938971942055369D+00, &
     -0.77456301686419115615D+00, &
     -0.43978938971942055369D+00, &
      0.72014146551340840752D+00, &
      0.22174790443082656455D+00, &
     -0.72014146551340840752D+00, &
      0.22174790443082656455D+00, &
      0.72014146551340840752D+00, &
     -0.22174790443082656455D+00, &
     -0.72014146551340840752D+00, &
     -0.22174790443082656455D+00, &
      0.17509849556381359981D+00, &
      0.58560797514860796209D+00, &
     -0.17509849556381359981D+00, &
      0.58560797514860796209D+00, &
      0.17509849556381359981D+00, &
     -0.58560797514860796209D+00, &
     -0.17509849556381359981D+00, &
     -0.58560797514860796209D+00, &
      0.93762547318055500245D+00, &
      0.73982675325131319610D+00, &
     -0.93762547318055500245D+00, &
      0.73982675325131319610D+00, &
      0.93762547318055500245D+00, &
     -0.73982675325131319610D+00, &
     -0.93762547318055500245D+00, &
     -0.73982675325131319610D+00, &
      0.92588822421370353677D+00, &
      0.50002093145521131490D+00, &
     -0.92588822421370353677D+00, &
      0.50002093145521131490D+00, &
      0.92588822421370353677D+00, &
     -0.50002093145521131490D+00, &
     -0.92588822421370353677D+00, &
     -0.50002093145521131490D+00, &
      0.26820739210416322251D+00, &
      0.56185198689199644662D+00, &
     -0.26820739210416322251D+00, &
      0.56185198689199644662D+00, &
      0.26820739210416322251D+00, &
     -0.56185198689199644662D+00, &
     -0.26820739210416322251D+00, &
     -0.56185198689199644662D+00, &
      0.36207959439737902319D+00, &
      0.79245628070986395830D+00, &
     -0.36207959439737902319D+00, &
      0.79245628070986395830D+00, &
      0.36207959439737902319D+00, &
     -0.79245628070986395830D+00, &
     -0.36207959439737902319D+00, &
     -0.79245628070986395830D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.41223892088328983885D+00, &
      0.00000000000000000000D+00, &
     -0.41223892088328983885D+00, &
      0.00000000000000000000D+00, &
      0.94482408797958883362D+00, &
      0.00000000000000000000D+00, &
     -0.94482408797958883362D+00, &
      0.00000000000000000000D+00, &
      0.15012306550295120844D+00, &
      0.00000000000000000000D+00, &
     -0.15012306550295120844D+00, &
      0.00000000000000000000D+00, &
      0.36338283577562219273D+00, &
      0.00000000000000000000D+00, &
     -0.36338283577562219273D+00, &
      0.00000000000000000000D+00, &
      0.32897122373450704558D+00, &
      0.00000000000000000000D+00, &
     -0.32897122373450704558D+00, &
      0.00000000000000000000D+00, &
      0.55704935891608620135D+00, &
      0.00000000000000000000D+00, &
     -0.55704935891608620135D+00, &
      0.00000000000000000000D+00, &
      0.71060354437421247020D+00, &
      0.00000000000000000000D+00, &
     -0.71060354437421247020D+00, &
      0.00000000000000000000D+00, &
      0.88531748708681257121D+00, &
      0.00000000000000000000D+00, &
     -0.88531748708681257121D+00, &
      0.48884754503841615358D+00, &
      0.48884754503841615358D+00, &
     -0.48884754503841615358D+00, &
     -0.48884754503841615358D+00, &
      0.45031772456672747307D+00, &
      0.45031772456672747307D+00, &
     -0.45031772456672747307D+00, &
     -0.45031772456672747307D+00, &
      0.65872956915547864476D+00, &
      0.65872956915547864476D+00, &
     -0.65872956915547864476D+00, &
     -0.65872956915547864476D+00, &
      0.62588271911798842861D+00, &
      0.62588271911798842861D+00, &
     -0.62588271911798842861D+00, &
     -0.62588271911798842861D+00, &
      0.25178960033193620305D+00, &
      0.25178960033193620305D+00, &
     -0.25178960033193620305D+00, &
     -0.25178960033193620305D+00, &
      0.91338191280449132492D+00, &
      0.91338191280449132492D+00, &
     -0.91338191280449132492D+00, &
     -0.91338191280449132492D+00, &
      0.11117861641720541699D+00, &
      0.11117861641720541699D+00, &
     -0.11117861641720541699D+00, &
     -0.11117861641720541699D+00, &
      0.37814924987375336807D+00, &
      0.37814924987375336807D+00, &
     -0.37814924987375336807D+00, &
     -0.37814924987375336807D+00, &
      0.23146195691577758913D+00, &
      0.23146195691577758913D+00, &
     -0.23146195691577758913D+00, &
     -0.23146195691577758913D+00, &
      0.18598667282576661353D+00, &
      0.18598667282576661353D+00, &
     -0.18598667282576661353D+00, &
     -0.18598667282576661353D+00, &
      0.72941490542892961635D+00, &
      0.72941490542892961635D+00, &
     -0.72941490542892961635D+00, &
     -0.72941490542892961635D+00, &
      0.80419241777926819825D+00, &
      0.80419241777926819825D+00, &
     -0.80419241777926819825D+00, &
     -0.80419241777926819825D+00, &
      0.47036082838417470064D+00, &
      0.47036082838417470064D+00, &
     -0.47036082838417470064D+00, &
     -0.47036082838417470064D+00, &
      0.27535657520011547206D+00, &
      0.27535657520011547206D+00, &
     -0.27535657520011547206D+00, &
     -0.27535657520011547206D+00, &
      0.25840443662578793660D+00, &
      0.25840443662578793660D+00, &
     -0.25840443662578793660D+00, &
     -0.25840443662578793660D+00, &
      0.43978938971942055369D+00, &
      0.77456301686419115615D+00, &
      0.43978938971942055369D+00, &
     -0.77456301686419115615D+00, &
     -0.43978938971942055369D+00, &
      0.77456301686419115615D+00, &
     -0.43978938971942055369D+00, &
     -0.77456301686419115615D+00, &
      0.22174790443082656455D+00, &
      0.72014146551340840752D+00, &
      0.22174790443082656455D+00, &
     -0.72014146551340840752D+00, &
     -0.22174790443082656455D+00, &
      0.72014146551340840752D+00, &
     -0.22174790443082656455D+00, &
     -0.72014146551340840752D+00, &
      0.58560797514860796209D+00, &
      0.17509849556381359981D+00, &
      0.58560797514860796209D+00, &
     -0.17509849556381359981D+00, &
     -0.58560797514860796209D+00, &
      0.17509849556381359981D+00, &
     -0.58560797514860796209D+00, &
     -0.17509849556381359981D+00, &
      0.73982675325131319610D+00, &
      0.93762547318055500245D+00, &
      0.73982675325131319610D+00, &
     -0.93762547318055500245D+00, &
     -0.73982675325131319610D+00, &
      0.93762547318055500245D+00, &
     -0.73982675325131319610D+00, &
     -0.93762547318055500245D+00, &
      0.50002093145521131490D+00, &
      0.92588822421370353677D+00, &
      0.50002093145521131490D+00, &
     -0.92588822421370353677D+00, &
     -0.50002093145521131490D+00, &
      0.92588822421370353677D+00, &
     -0.50002093145521131490D+00, &
     -0.92588822421370353677D+00, &
      0.56185198689199644662D+00, &
      0.26820739210416322251D+00, &
      0.56185198689199644662D+00, &
     -0.26820739210416322251D+00, &
     -0.56185198689199644662D+00, &
      0.26820739210416322251D+00, &
     -0.56185198689199644662D+00, &
     -0.26820739210416322251D+00, &
      0.79245628070986395830D+00, &
      0.36207959439737902319D+00, &
      0.79245628070986395830D+00, &
     -0.36207959439737902319D+00, &
     -0.79245628070986395830D+00, &
      0.36207959439737902319D+00, &
     -0.79245628070986395830D+00, &
     -0.36207959439737902319D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
      0.92974005864256004106D+00, &
      0.58387036044479012631D+00, &
      0.10942458150673148309D+00, &
      0.30043598129646631456D+00, &
      0.40351765144169587929D+00, &
      0.40351765144169587929D+00, &
      0.40351765144169587929D+00, &
      0.40351765144169587929D+00, &
      0.00554139647030791119D+00, &
      0.00554139647030791119D+00, &
      0.00554139647030791119D+00, &
      0.00554139647030791119D+00, &
      0.84979294403286831372D+00, &
      0.84979294403286831372D+00, &
      0.84979294403286831372D+00, &
      0.84979294403286831372D+00, &
      0.60719504362173426504D+00, &
      0.60719504362173426504D+00, &
      0.60719504362173426504D+00, &
      0.60719504362173426504D+00, &
      0.02033792783786766978D+00, &
      0.02033792783786766978D+00, &
      0.02033792783786766978D+00, &
      0.02033792783786766978D+00, &
      0.10601199873730246526D+00, &
      0.10601199873730246526D+00, &
      0.10601199873730246526D+00, &
      0.10601199873730246526D+00, &
      0.24910419915320258788D+00, &
      0.24910419915320258788D+00, &
      0.24910419915320258788D+00, &
      0.24910419915320258788D+00, &
      0.08107956240970674855D+00, &
      0.08107956240970674855D+00, &
      0.08107956240970674855D+00, &
      0.08107956240970674855D+00, &
      0.00950815168106079051D+00, &
      0.00950815168106079051D+00, &
      0.00950815168106079051D+00, &
      0.00950815168106079051D+00, &
      0.48819360509214093646D+00, &
      0.48819360509214093646D+00, &
      0.48819360509214093646D+00, &
      0.48819360509214093646D+00, &
      0.27690831156608641805D+00, &
      0.27690831156608641805D+00, &
      0.27690831156608641805D+00, &
      0.27690831156608641805D+00, &
      0.12248093820428122835D+00, &
      0.12248093820428122835D+00, &
      0.12248093820428122835D+00, &
      0.12248093820428122835D+00, &
      0.22323057279864513824D+00, &
      0.22323057279864513824D+00, &
      0.22323057279864513824D+00, &
      0.22323057279864513824D+00, &
      0.01207165002118426068D+00, &
      0.01207165002118426068D+00, &
      0.01207165002118426068D+00, &
      0.01207165002118426068D+00, &
      0.75518724716690510679D+00, &
      0.75518724716690510679D+00, &
      0.75518724716690510679D+00, &
      0.75518724716690510679D+00, &
      0.08424825204838617965D+00, &
      0.08424825204838617965D+00, &
      0.08424825204838617965D+00, &
      0.08424825204838617965D+00, &
      0.63188577231776799081D+00, &
      0.63188577231776799081D+00, &
      0.63188577231776799081D+00, &
      0.63188577231776799081D+00, &
      0.39892977659234107879D+00, &
      0.39892977659234107879D+00, &
      0.39892977659234107879D+00, &
      0.39892977659234107879D+00, &
      0.02502635713295603415D+00, &
      0.02502635713295603415D+00, &
      0.02502635713295603415D+00, &
      0.02502635713295603415D+00, &
      0.12038000765135539738D+00, &
      0.12038000765135539738D+00, &
      0.12038000765135539738D+00, &
      0.12038000765135539738D+00, &
      0.29270359357039776871D+00, &
      0.29270359357039776871D+00, &
      0.29270359357039776871D+00, &
      0.29270359357039776871D+00, &
      0.70742748921558806785D+00, &
      0.70742748921558806785D+00, &
      0.70742748921558806785D+00, &
      0.70742748921558806785D+00, &
      0.52832746509410000169D+00, &
      0.52832746509410000169D+00, &
      0.52832746509410000169D+00, &
      0.52832746509410000169D+00, &
      0.19538409012904675577D+00, &
      0.19538409012904675577D+00, &
      0.19538409012904675577D+00, &
      0.19538409012904675577D+00, &
      0.19538409012904675577D+00, &
      0.19538409012904675577D+00, &
      0.19538409012904675577D+00, &
      0.19538409012904675577D+00, &
      0.02260035188498128386D+00, &
      0.02260035188498128386D+00, &
      0.02260035188498128386D+00, &
      0.02260035188498128386D+00, &
      0.02260035188498128386D+00, &
      0.02260035188498128386D+00, &
      0.02260035188498128386D+00, &
      0.02260035188498128386D+00, &
      0.21608520533637087802D+00, &
      0.21608520533637087802D+00, &
      0.21608520533637087802D+00, &
      0.21608520533637087802D+00, &
      0.21608520533637087802D+00, &
      0.21608520533637087802D+00, &
      0.21608520533637087802D+00, &
      0.21608520533637087802D+00, &
      0.06232437306564607427D+00, &
      0.06232437306564607427D+00, &
      0.06232437306564607427D+00, &
      0.06232437306564607427D+00, &
      0.06232437306564607427D+00, &
      0.06232437306564607427D+00, &
      0.06232437306564607427D+00, &
      0.06232437306564607427D+00, &
      0.01929314718254383776D+00, &
      0.01929314718254383776D+00, &
      0.01929314718254383776D+00, &
      0.01929314718254383776D+00, &
      0.01929314718254383776D+00, &
      0.01929314718254383776D+00, &
      0.01929314718254383776D+00, &
      0.01929314718254383776D+00, &
      0.40102466550487736452D+00, &
      0.40102466550487736452D+00, &
      0.40102466550487736452D+00, &
      0.40102466550487736452D+00, &
      0.40102466550487736452D+00, &
      0.40102466550487736452D+00, &
      0.40102466550487736452D+00, &
      0.40102466550487736452D+00, &
      0.09200529105614789482D+00, &
      0.09200529105614789482D+00, &
      0.09200529105614789482D+00, &
      0.09200529105614789482D+00, &
      0.09200529105614789482D+00, &
      0.09200529105614789482D+00, &
      0.09200529105614789482D+00, &
      0.09200529105614789482D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
      0.00120366956557754735D+00, &
      0.01890734566258640142D+00, &
      0.02412883548353303431D+00, &
      0.01430958560741040810D+00, &
      0.01462233467606358257D+00, &
      0.01462233467606358257D+00, &
      0.01462233467606358257D+00, &
      0.01462233467606358257D+00, &
      0.00131047133506746131D+00, &
      0.00131047133506746131D+00, &
      0.00131047133506746131D+00, &
      0.00131047133506746131D+00, &
      0.00131534859459711507D+00, &
      0.00131534859459711507D+00, &
      0.00131534859459711507D+00, &
      0.00131534859459711507D+00, &
      0.00604172938779333306D+00, &
      0.00604172938779333306D+00, &
      0.00604172938779333306D+00, &
      0.00604172938779333306D+00, &
      0.00828934369195859742D+00, &
      0.00828934369195859742D+00, &
      0.00828934369195859742D+00, &
      0.00828934369195859742D+00, &
      0.01293341373858935758D+00, &
      0.01293341373858935758D+00, &
      0.01293341373858935758D+00, &
      0.01293341373858935758D+00, &
      0.00563885419154383245D+00, &
      0.00563885419154383245D+00, &
      0.00563885419154383245D+00, &
      0.00563885419154383245D+00, &
      0.00450499139881360439D+00, &
      0.00450499139881360439D+00, &
      0.00450499139881360439D+00, &
      0.00450499139881360439D+00, &
      0.00363618420610894565D+00, &
      0.00363618420610894565D+00, &
      0.00363618420610894565D+00, &
      0.00363618420610894565D+00, &
      0.00388842116802182820D+00, &
      0.00388842116802182820D+00, &
      0.00388842116802182820D+00, &
      0.00388842116802182820D+00, &
      0.00402365320826360111D+00, &
      0.00402365320826360111D+00, &
      0.00402365320826360111D+00, &
      0.00402365320826360111D+00, &
      0.00955901949244871221D+00, &
      0.00955901949244871221D+00, &
      0.00955901949244871221D+00, &
      0.00955901949244871221D+00, &
      0.01938656716858040696D+00, &
      0.01938656716858040696D+00, &
      0.01938656716858040696D+00, &
      0.01938656716858040696D+00, &
      0.00096875488373819618D+00, &
      0.00096875488373819618D+00, &
      0.00096875488373819618D+00, &
      0.00096875488373819618D+00, &
      0.00502214350286199571D+00, &
      0.00502214350286199571D+00, &
      0.00502214350286199571D+00, &
      0.00502214350286199571D+00, &
      0.01448769412031481602D+00, &
      0.01448769412031481602D+00, &
      0.01448769412031481602D+00, &
      0.01448769412031481602D+00, &
      0.00459976431383042356D+00, &
      0.00459976431383042356D+00, &
      0.00459976431383042356D+00, &
      0.00459976431383042356D+00, &
      0.01107636061349879207D+00, &
      0.01107636061349879207D+00, &
      0.01107636061349879207D+00, &
      0.01107636061349879207D+00, &
      0.00469260200432727610D+00, &
      0.00469260200432727610D+00, &
      0.00469260200432727610D+00, &
      0.00469260200432727610D+00, &
      0.00285194982449880072D+00, &
      0.00285194982449880072D+00, &
      0.00285194982449880072D+00, &
      0.00285194982449880072D+00, &
      0.01323033835935080553D+00, &
      0.01323033835935080553D+00, &
      0.01323033835935080553D+00, &
      0.01323033835935080553D+00, &
      0.00111972343972882474D+00, &
      0.00111972343972882474D+00, &
      0.00111972343972882474D+00, &
      0.00111972343972882474D+00, &
      0.00918246462138119605D+00, &
      0.00918246462138119605D+00, &
      0.00918246462138119605D+00, &
      0.00918246462138119605D+00, &
      0.00456490653701050169D+00, &
      0.00456490653701050169D+00, &
      0.00456490653701050169D+00, &
      0.00456490653701050169D+00, &
      0.00456490653701050169D+00, &
      0.00456490653701050169D+00, &
      0.00456490653701050169D+00, &
      0.00456490653701050169D+00, &
      0.00547503483442879964D+00, &
      0.00547503483442879964D+00, &
      0.00547503483442879964D+00, &
      0.00547503483442879964D+00, &
      0.00547503483442879964D+00, &
      0.00547503483442879964D+00, &
      0.00547503483442879964D+00, &
      0.00547503483442879964D+00, &
      0.00992261207645742799D+00, &
      0.00992261207645742799D+00, &
      0.00992261207645742799D+00, &
      0.00992261207645742799D+00, &
      0.00992261207645742799D+00, &
      0.00992261207645742799D+00, &
      0.00992261207645742799D+00, &
      0.00992261207645742799D+00, &
      0.00110382912156457261D+00, &
      0.00110382912156457261D+00, &
      0.00110382912156457261D+00, &
      0.00110382912156457261D+00, &
      0.00110382912156457261D+00, &
      0.00110382912156457261D+00, &
      0.00110382912156457261D+00, &
      0.00110382912156457261D+00, &
      0.00236228293958670472D+00, &
      0.00236228293958670472D+00, &
      0.00236228293958670472D+00, &
      0.00236228293958670472D+00, &
      0.00236228293958670472D+00, &
      0.00236228293958670472D+00, &
      0.00236228293958670472D+00, &
      0.00236228293958670472D+00, &
      0.00630666584934871672D+00, &
      0.00630666584934871672D+00, &
      0.00630666584934871672D+00, &
      0.00630666584934871672D+00, &
      0.00630666584934871672D+00, &
      0.00630666584934871672D+00, &
      0.00630666584934871672D+00, &
      0.00630666584934871672D+00, &
      0.00675492513102409885D+00, &
      0.00675492513102409885D+00, &
      0.00675492513102409885D+00, &
      0.00675492513102409885D+00, &
      0.00675492513102409885D+00, &
      0.00675492513102409885D+00, &
      0.00675492513102409885D+00, &
      0.00675492513102409885D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule14 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule14() returns the pyramid quadrature rule of precision 14.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 April 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 121, Number 11, pages 2418-2436, 15 June 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real ( kind = rk ) x[n], y[n], z[n]: the coordinates of
!    quadrature points.
!
!    real ( kind = rk ) w[n]: the quadrature weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 184

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.42182613673485974681D+00, &
      0.00000000000000000000D+00, &
     -0.42182613673485974681D+00, &
      0.00000000000000000000D+00, &
      0.08465613435886885918D+00, &
      0.00000000000000000000D+00, &
     -0.08465613435886885918D+00, &
      0.00000000000000000000D+00, &
      0.95167131965923157377D+00, &
      0.00000000000000000000D+00, &
     -0.95167131965923157377D+00, &
      0.00000000000000000000D+00, &
      0.12212881781145627780D+00, &
      0.00000000000000000000D+00, &
     -0.12212881781145627780D+00, &
      0.00000000000000000000D+00, &
      0.34597557868730360875D+00, &
      0.00000000000000000000D+00, &
     -0.34597557868730360875D+00, &
      0.00000000000000000000D+00, &
      0.32500209212470398956D+00, &
      0.00000000000000000000D+00, &
     -0.32500209212470398956D+00, &
      0.00000000000000000000D+00, &
      0.85274656175184737616D+00, &
      0.00000000000000000000D+00, &
     -0.85274656175184737616D+00, &
      0.00000000000000000000D+00, &
      0.63827756477152641779D+00, &
      0.00000000000000000000D+00, &
     -0.63827756477152641779D+00, &
      0.00000000000000000000D+00, &
      0.67231147414585012978D+00, &
      0.00000000000000000000D+00, &
     -0.67231147414585012978D+00, &
      0.00000000000000000000D+00, &
      0.85374960595109017358D+00, &
      0.00000000000000000000D+00, &
     -0.85374960595109017358D+00, &
      0.00000000000000000000D+00, &
      0.32351828528149095821D+00, &
      0.00000000000000000000D+00, &
     -0.32351828528149095821D+00, &
      0.00000000000000000000D+00, &
      0.49697514832376565863D+00, &
      0.00000000000000000000D+00, &
     -0.49697514832376565863D+00, &
      0.00000000000000000000D+00, &
      0.50417062955844949013D+00, &
     -0.50417062955844949013D+00, &
      0.50417062955844949013D+00, &
     -0.50417062955844949013D+00, &
      0.63853202881303905425D+00, &
     -0.63853202881303905425D+00, &
      0.63853202881303905425D+00, &
     -0.63853202881303905425D+00, &
      0.76432244824860651189D+00, &
     -0.76432244824860651189D+00, &
      0.76432244824860651189D+00, &
     -0.76432244824860651189D+00, &
      0.56133306193833287789D+00, &
     -0.56133306193833287789D+00, &
      0.56133306193833287789D+00, &
     -0.56133306193833287789D+00, &
      0.28786805395155073972D+00, &
     -0.28786805395155073972D+00, &
      0.28786805395155073972D+00, &
     -0.28786805395155073972D+00, &
      0.90163083827576329110D+00, &
     -0.90163083827576329110D+00, &
      0.90163083827576329110D+00, &
     -0.90163083827576329110D+00, &
      0.09459747575878657555D+00, &
     -0.09459747575878657555D+00, &
      0.09459747575878657555D+00, &
     -0.09459747575878657555D+00, &
      0.26275744057316757774D+00, &
     -0.26275744057316757774D+00, &
      0.26275744057316757774D+00, &
     -0.26275744057316757774D+00, &
      0.38640915363814193340D+00, &
     -0.38640915363814193340D+00, &
      0.38640915363814193340D+00, &
     -0.38640915363814193340D+00, &
      0.22823005779855418118D+00, &
     -0.22823005779855418118D+00, &
      0.22823005779855418118D+00, &
     -0.22823005779855418118D+00, &
      0.70578722692835427210D+00, &
     -0.70578722692835427210D+00, &
      0.70578722692835427210D+00, &
     -0.70578722692835427210D+00, &
      0.91478755244276455105D+00, &
     -0.91478755244276455105D+00, &
      0.91478755244276455105D+00, &
     -0.91478755244276455105D+00, &
      0.50883064513682729757D+00, &
     -0.50883064513682729757D+00, &
      0.50883064513682729757D+00, &
     -0.50883064513682729757D+00, &
      0.21833905510549736495D+00, &
     -0.21833905510549736495D+00, &
      0.21833905510549736495D+00, &
     -0.21833905510549736495D+00, &
      0.20524210043072979581D+00, &
     -0.20524210043072979581D+00, &
      0.20524210043072979581D+00, &
     -0.20524210043072979581D+00, &
      0.29109013901653718603D+00, &
     -0.29109013901653718603D+00, &
      0.29109013901653718603D+00, &
     -0.29109013901653718603D+00, &
      0.74865402480464293689D+00, &
     -0.74865402480464293689D+00, &
      0.74865402480464293689D+00, &
     -0.74865402480464293689D+00, &
      0.73782464660583058080D+00, &
      0.49284956480916058963D+00, &
     -0.73782464660583058080D+00, &
      0.49284956480916058963D+00, &
      0.73782464660583058080D+00, &
     -0.49284956480916058963D+00, &
     -0.73782464660583058080D+00, &
     -0.49284956480916058963D+00, &
      0.85890606200982022589D+00, &
      0.37381968226304934655D+00, &
     -0.85890606200982022589D+00, &
      0.37381968226304934655D+00, &
      0.85890606200982022589D+00, &
     -0.37381968226304934655D+00, &
     -0.85890606200982022589D+00, &
     -0.37381968226304934655D+00, &
      0.22436150649227576404D+00, &
      0.63179494567048766207D+00, &
     -0.22436150649227576404D+00, &
      0.63179494567048766207D+00, &
      0.22436150649227576404D+00, &
     -0.63179494567048766207D+00, &
     -0.22436150649227576404D+00, &
     -0.63179494567048766207D+00, &
      0.90017577589364738966D+00, &
      0.62618678850176523465D+00, &
     -0.90017577589364738966D+00, &
      0.62618678850176523465D+00, &
      0.90017577589364738966D+00, &
     -0.62618678850176523465D+00, &
     -0.90017577589364738966D+00, &
     -0.62618678850176523465D+00, &
      0.97895018498408203911D+00, &
      0.62565783993238355265D+00, &
     -0.97895018498408203911D+00, &
      0.62565783993238355265D+00, &
      0.97895018498408203911D+00, &
     -0.62565783993238355265D+00, &
     -0.97895018498408203911D+00, &
     -0.62565783993238355265D+00, &
      0.25607583733637889756D+00, &
      0.51194810665765877467D+00, &
     -0.25607583733637889756D+00, &
      0.51194810665765877467D+00, &
      0.25607583733637889756D+00, &
     -0.51194810665765877467D+00, &
     -0.25607583733637889756D+00, &
     -0.51194810665765877467D+00, &
      0.29741762068713811784D+00, &
      0.77553958650832321986D+00, &
     -0.29741762068713811784D+00, &
      0.77553958650832321986D+00, &
      0.29741762068713811784D+00, &
     -0.77553958650832321986D+00, &
     -0.29741762068713811784D+00, &
     -0.77553958650832321986D+00, &
      0.29079659566096877077D+00, &
      0.61926152278826551711D+00, &
     -0.29079659566096877077D+00, &
      0.61926152278826551711D+00, &
      0.29079659566096877077D+00, &
     -0.61926152278826551711D+00, &
     -0.29079659566096877077D+00, &
     -0.61926152278826551711D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.42182613673485974681D+00, &
      0.00000000000000000000D+00, &
     -0.42182613673485974681D+00, &
      0.00000000000000000000D+00, &
      0.08465613435886885918D+00, &
      0.00000000000000000000D+00, &
     -0.08465613435886885918D+00, &
      0.00000000000000000000D+00, &
      0.95167131965923157377D+00, &
      0.00000000000000000000D+00, &
     -0.95167131965923157377D+00, &
      0.00000000000000000000D+00, &
      0.12212881781145627780D+00, &
      0.00000000000000000000D+00, &
     -0.12212881781145627780D+00, &
      0.00000000000000000000D+00, &
      0.34597557868730360875D+00, &
      0.00000000000000000000D+00, &
     -0.34597557868730360875D+00, &
      0.00000000000000000000D+00, &
      0.32500209212470398956D+00, &
      0.00000000000000000000D+00, &
     -0.32500209212470398956D+00, &
      0.00000000000000000000D+00, &
      0.85274656175184737616D+00, &
      0.00000000000000000000D+00, &
     -0.85274656175184737616D+00, &
      0.00000000000000000000D+00, &
      0.63827756477152641779D+00, &
      0.00000000000000000000D+00, &
     -0.63827756477152641779D+00, &
      0.00000000000000000000D+00, &
      0.67231147414585012978D+00, &
      0.00000000000000000000D+00, &
     -0.67231147414585012978D+00, &
      0.00000000000000000000D+00, &
      0.85374960595109017358D+00, &
      0.00000000000000000000D+00, &
     -0.85374960595109017358D+00, &
      0.00000000000000000000D+00, &
      0.32351828528149095821D+00, &
      0.00000000000000000000D+00, &
     -0.32351828528149095821D+00, &
      0.00000000000000000000D+00, &
      0.49697514832376565863D+00, &
      0.00000000000000000000D+00, &
     -0.49697514832376565863D+00, &
      0.50417062955844949013D+00, &
      0.50417062955844949013D+00, &
     -0.50417062955844949013D+00, &
     -0.50417062955844949013D+00, &
      0.63853202881303905425D+00, &
      0.63853202881303905425D+00, &
     -0.63853202881303905425D+00, &
     -0.63853202881303905425D+00, &
      0.76432244824860651189D+00, &
      0.76432244824860651189D+00, &
     -0.76432244824860651189D+00, &
     -0.76432244824860651189D+00, &
      0.56133306193833287789D+00, &
      0.56133306193833287789D+00, &
     -0.56133306193833287789D+00, &
     -0.56133306193833287789D+00, &
      0.28786805395155073972D+00, &
      0.28786805395155073972D+00, &
     -0.28786805395155073972D+00, &
     -0.28786805395155073972D+00, &
      0.90163083827576329110D+00, &
      0.90163083827576329110D+00, &
     -0.90163083827576329110D+00, &
     -0.90163083827576329110D+00, &
      0.09459747575878657555D+00, &
      0.09459747575878657555D+00, &
     -0.09459747575878657555D+00, &
     -0.09459747575878657555D+00, &
      0.26275744057316757774D+00, &
      0.26275744057316757774D+00, &
     -0.26275744057316757774D+00, &
     -0.26275744057316757774D+00, &
      0.38640915363814193340D+00, &
      0.38640915363814193340D+00, &
     -0.38640915363814193340D+00, &
     -0.38640915363814193340D+00, &
      0.22823005779855418118D+00, &
      0.22823005779855418118D+00, &
     -0.22823005779855418118D+00, &
     -0.22823005779855418118D+00, &
      0.70578722692835427210D+00, &
      0.70578722692835427210D+00, &
     -0.70578722692835427210D+00, &
     -0.70578722692835427210D+00, &
      0.91478755244276455105D+00, &
      0.91478755244276455105D+00, &
     -0.91478755244276455105D+00, &
     -0.91478755244276455105D+00, &
      0.50883064513682729757D+00, &
      0.50883064513682729757D+00, &
     -0.50883064513682729757D+00, &
     -0.50883064513682729757D+00, &
      0.21833905510549736495D+00, &
      0.21833905510549736495D+00, &
     -0.21833905510549736495D+00, &
     -0.21833905510549736495D+00, &
      0.20524210043072979581D+00, &
      0.20524210043072979581D+00, &
     -0.20524210043072979581D+00, &
     -0.20524210043072979581D+00, &
      0.29109013901653718603D+00, &
      0.29109013901653718603D+00, &
     -0.29109013901653718603D+00, &
     -0.29109013901653718603D+00, &
      0.74865402480464293689D+00, &
      0.74865402480464293689D+00, &
     -0.74865402480464293689D+00, &
     -0.74865402480464293689D+00, &
      0.49284956480916058963D+00, &
      0.73782464660583058080D+00, &
      0.49284956480916058963D+00, &
     -0.73782464660583058080D+00, &
     -0.49284956480916058963D+00, &
      0.73782464660583058080D+00, &
     -0.49284956480916058963D+00, &
     -0.73782464660583058080D+00, &
      0.37381968226304934655D+00, &
      0.85890606200982022589D+00, &
      0.37381968226304934655D+00, &
     -0.85890606200982022589D+00, &
     -0.37381968226304934655D+00, &
      0.85890606200982022589D+00, &
     -0.37381968226304934655D+00, &
     -0.85890606200982022589D+00, &
      0.63179494567048766207D+00, &
      0.22436150649227576404D+00, &
      0.63179494567048766207D+00, &
     -0.22436150649227576404D+00, &
     -0.63179494567048766207D+00, &
      0.22436150649227576404D+00, &
     -0.63179494567048766207D+00, &
     -0.22436150649227576404D+00, &
      0.62618678850176523465D+00, &
      0.90017577589364738966D+00, &
      0.62618678850176523465D+00, &
     -0.90017577589364738966D+00, &
     -0.62618678850176523465D+00, &
      0.90017577589364738966D+00, &
     -0.62618678850176523465D+00, &
     -0.90017577589364738966D+00, &
      0.62565783993238355265D+00, &
      0.97895018498408203911D+00, &
      0.62565783993238355265D+00, &
     -0.97895018498408203911D+00, &
     -0.62565783993238355265D+00, &
      0.97895018498408203911D+00, &
     -0.62565783993238355265D+00, &
     -0.97895018498408203911D+00, &
      0.51194810665765877467D+00, &
      0.25607583733637889756D+00, &
      0.51194810665765877467D+00, &
     -0.25607583733637889756D+00, &
     -0.51194810665765877467D+00, &
      0.25607583733637889756D+00, &
     -0.51194810665765877467D+00, &
     -0.25607583733637889756D+00, &
      0.77553958650832321986D+00, &
      0.29741762068713811784D+00, &
      0.77553958650832321986D+00, &
     -0.29741762068713811784D+00, &
     -0.77553958650832321986D+00, &
      0.29741762068713811784D+00, &
     -0.77553958650832321986D+00, &
     -0.29741762068713811784D+00, &
      0.61926152278826551711D+00, &
      0.29079659566096877077D+00, &
      0.61926152278826551711D+00, &
     -0.29079659566096877077D+00, &
     -0.61926152278826551711D+00, &
      0.29079659566096877077D+00, &
     -0.61926152278826551711D+00, &
     -0.29079659566096877077D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
      0.94550960783801574205D+00, &
      0.61382435648399280570D+00, &
      0.24540497180596210214D+00, &
      0.41832825390859845749D+00, &
      0.40175363160567706400D+00, &
      0.40175363160567706400D+00, &
      0.40175363160567706400D+00, &
      0.40175363160567706400D+00, &
      0.08958821481037450296D+00, &
      0.08958821481037450296D+00, &
      0.08958821481037450296D+00, &
      0.08958821481037450296D+00, &
      0.03128290548598175458D+00, &
      0.03128290548598175458D+00, &
      0.03128290548598175458D+00, &
      0.03128290548598175458D+00, &
      0.86562635222401107526D+00, &
      0.86562635222401107526D+00, &
      0.86562635222401107526D+00, &
      0.86562635222401107526D+00, &
      0.56154915636029745230D+00, &
      0.56154915636029745230D+00, &
      0.56154915636029745230D+00, &
      0.56154915636029745230D+00, &
      0.01398939848661482320D+00, &
      0.01398939848661482320D+00, &
      0.01398939848661482320D+00, &
      0.01398939848661482320D+00, &
      0.00165571829628404746D+00, &
      0.00165571829628404746D+00, &
      0.00165571829628404746D+00, &
      0.00165571829628404746D+00, &
      0.06493816225859559699D+00, &
      0.06493816225859559699D+00, &
      0.06493816225859559699D+00, &
      0.06493816225859559699D+00, &
      0.32265244995852410126D+00, &
      0.32265244995852410126D+00, &
      0.32265244995852410126D+00, &
      0.32265244995852410126D+00, &
      0.14284602353229125526D+00, &
      0.14284602353229125526D+00, &
      0.14284602353229125526D+00, &
      0.14284602353229125526D+00, &
      0.66509821819012426847D+00, &
      0.66509821819012426847D+00, &
      0.66509821819012426847D+00, &
      0.66509821819012426847D+00, &
      0.19278467665636450645D+00, &
      0.19278467665636450645D+00, &
      0.19278467665636450645D+00, &
      0.19278467665636450645D+00, &
      0.05324772752904030626D+00, &
      0.05324772752904030626D+00, &
      0.05324772752904030626D+00, &
      0.05324772752904030626D+00, &
      0.34999224799132661046D+00, &
      0.34999224799132661046D+00, &
      0.34999224799132661046D+00, &
      0.34999224799132661046D+00, &
      0.16096112315804175785D+00, &
      0.16096112315804175785D+00, &
      0.16096112315804175785D+00, &
      0.16096112315804175785D+00, &
      0.16039988924094228384D+00, &
      0.16039988924094228384D+00, &
      0.16039988924094228384D+00, &
      0.16039988924094228384D+00, &
      0.13732763905737063737D+00, &
      0.13732763905737063737D+00, &
      0.13732763905737063737D+00, &
      0.13732763905737063737D+00, &
      0.01032358209500776509D+00, &
      0.01032358209500776509D+00, &
      0.01032358209500776509D+00, &
      0.01032358209500776509D+00, &
      0.76779179416407739023D+00, &
      0.76779179416407739023D+00, &
      0.76779179416407739023D+00, &
      0.76779179416407739023D+00, &
      0.04943554173327800727D+00, &
      0.04943554173327800727D+00, &
      0.04943554173327800727D+00, &
      0.04943554173327800727D+00, &
      0.55602661412166731747D+00, &
      0.55602661412166731747D+00, &
      0.55602661412166731747D+00, &
      0.55602661412166731747D+00, &
      0.49627919975552486909D+00, &
      0.49627919975552486909D+00, &
      0.49627919975552486909D+00, &
      0.49627919975552486909D+00, &
      0.00200880715410161345D+00, &
      0.00200880715410161345D+00, &
      0.00200880715410161345D+00, &
      0.00200880715410161345D+00, &
      0.04956231677804342345D+00, &
      0.04956231677804342345D+00, &
      0.04956231677804342345D+00, &
      0.04956231677804342345D+00, &
      0.34635586968731740809D+00, &
      0.34635586968731740809D+00, &
      0.34635586968731740809D+00, &
      0.34635586968731740809D+00, &
      0.74542883631057244020D+00, &
      0.74542883631057244020D+00, &
      0.74542883631057244020D+00, &
      0.74542883631057244020D+00, &
      0.64153164712892607469D+00, &
      0.64153164712892607469D+00, &
      0.64153164712892607469D+00, &
      0.64153164712892607469D+00, &
      0.30622738335271937338D+00, &
      0.30622738335271937338D+00, &
      0.30622738335271937338D+00, &
      0.30622738335271937338D+00, &
      0.04936880562940795109D+00, &
      0.04936880562940795109D+00, &
      0.04936880562940795109D+00, &
      0.04936880562940795109D+00, &
      0.23699400699902964385D+00, &
      0.23699400699902964385D+00, &
      0.23699400699902964385D+00, &
      0.23699400699902964385D+00, &
      0.23699400699902964385D+00, &
      0.23699400699902964385D+00, &
      0.23699400699902964385D+00, &
      0.23699400699902964385D+00, &
      0.02776225858911670133D+00, &
      0.02776225858911670133D+00, &
      0.02776225858911670133D+00, &
      0.02776225858911670133D+00, &
      0.02776225858911670133D+00, &
      0.02776225858911670133D+00, &
      0.02776225858911670133D+00, &
      0.02776225858911670133D+00, &
      0.25928305502878512545D+00, &
      0.25928305502878512545D+00, &
      0.25928305502878512545D+00, &
      0.25928305502878512545D+00, &
      0.25928305502878512545D+00, &
      0.25928305502878512545D+00, &
      0.25928305502878512545D+00, &
      0.25928305502878512545D+00, &
      0.07806621986847779582D+00, &
      0.07806621986847779582D+00, &
      0.07806621986847779582D+00, &
      0.07806621986847779582D+00, &
      0.07806621986847779582D+00, &
      0.07806621986847779582D+00, &
      0.07806621986847779582D+00, &
      0.07806621986847779582D+00, &
      0.00443001451925785755D+00, &
      0.00443001451925785755D+00, &
      0.00443001451925785755D+00, &
      0.00443001451925785755D+00, &
      0.00443001451925785755D+00, &
      0.00443001451925785755D+00, &
      0.00443001451925785755D+00, &
      0.00443001451925785755D+00, &
      0.45428572768859826203D+00, &
      0.45428572768859826203D+00, &
      0.45428572768859826203D+00, &
      0.45428572768859826203D+00, &
      0.45428572768859826203D+00, &
      0.45428572768859826203D+00, &
      0.45428572768859826203D+00, &
      0.45428572768859826203D+00, &
      0.11761104673332581361D+00, &
      0.11761104673332581361D+00, &
      0.11761104673332581361D+00, &
      0.11761104673332581361D+00, &
      0.11761104673332581361D+00, &
      0.11761104673332581361D+00, &
      0.11761104673332581361D+00, &
      0.11761104673332581361D+00, &
      0.00783483739567485994D+00, &
      0.00783483739567485994D+00, &
      0.00783483739567485994D+00, &
      0.00783483739567485994D+00, &
      0.00783483739567485994D+00, &
      0.00783483739567485994D+00, &
      0.00783483739567485994D+00, &
      0.00783483739567485994D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
      0.00060098284109148011D+00, &
      0.01255487583666288487D+00, &
      0.02319377258130698310D+00, &
      0.01905556551345680832D+00, &
      0.01191795956205729361D+00, &
      0.01191795956205729361D+00, &
      0.01191795956205729361D+00, &
      0.01191795956205729361D+00, &
      0.00386718441233280599D+00, &
      0.00386718441233280599D+00, &
      0.00386718441233280599D+00, &
      0.00386718441233280599D+00, &
      0.00192376309517474221D+00, &
      0.00192376309517474221D+00, &
      0.00192376309517474221D+00, &
      0.00192376309517474221D+00, &
      0.00123417980449306002D+00, &
      0.00123417980449306002D+00, &
      0.00123417980449306002D+00, &
      0.00123417980449306002D+00, &
      0.00471664803911589926D+00, &
      0.00471664803911589926D+00, &
      0.00471664803911589926D+00, &
      0.00471664803911589926D+00, &
      0.00513404497338453437D+00, &
      0.00513404497338453437D+00, &
      0.00513404497338453437D+00, &
      0.00513404497338453437D+00, &
      0.00115456474640818034D+00, &
      0.00115456474640818034D+00, &
      0.00115456474640818034D+00, &
      0.00115456474640818034D+00, &
      0.01181264118995303640D+00, &
      0.01181264118995303640D+00, &
      0.01181264118995303640D+00, &
      0.01181264118995303640D+00, &
      0.00243693783463101070D+00, &
      0.00243693783463101070D+00, &
      0.00243693783463101070D+00, &
      0.00243693783463101070D+00, &
      0.00244218039395145909D+00, &
      0.00244218039395145909D+00, &
      0.00244218039395145909D+00, &
      0.00244218039395145909D+00, &
      0.00276027206518877129D+00, &
      0.00276027206518877129D+00, &
      0.00276027206518877129D+00, &
      0.00276027206518877129D+00, &
      0.01491875937488000743D+00, &
      0.01491875937488000743D+00, &
      0.01491875937488000743D+00, &
      0.01491875937488000743D+00, &
      0.00860788548052002223D+00, &
      0.00860788548052002223D+00, &
      0.00860788548052002223D+00, &
      0.00860788548052002223D+00, &
      0.00092010540262242146D+00, &
      0.00092010540262242146D+00, &
      0.00092010540262242146D+00, &
      0.00092010540262242146D+00, &
      0.00350658071335885509D+00, &
      0.00350658071335885509D+00, &
      0.00350658071335885509D+00, &
      0.00350658071335885509D+00, &
      0.01184517921997258437D+00, &
      0.01184517921997258437D+00, &
      0.01184517921997258437D+00, &
      0.01184517921997258437D+00, &
      0.01517200950119655840D+00, &
      0.01517200950119655840D+00, &
      0.01517200950119655840D+00, &
      0.01517200950119655840D+00, &
      0.00079885709646654790D+00, &
      0.00079885709646654790D+00, &
      0.00079885709646654790D+00, &
      0.00079885709646654790D+00, &
      0.00346532205672137571D+00, &
      0.00346532205672137571D+00, &
      0.00346532205672137571D+00, &
      0.00346532205672137571D+00, &
      0.00564901752891650531D+00, &
      0.00564901752891650531D+00, &
      0.00564901752891650531D+00, &
      0.00564901752891650531D+00, &
      0.00314499257654337911D+00, &
      0.00314499257654337911D+00, &
      0.00314499257654337911D+00, &
      0.00314499257654337911D+00, &
      0.01042467052306092316D+00, &
      0.01042467052306092316D+00, &
      0.01042467052306092316D+00, &
      0.01042467052306092316D+00, &
      0.00135335334419937545D+00, &
      0.00135335334419937545D+00, &
      0.00135335334419937545D+00, &
      0.00135335334419937545D+00, &
      0.00068639460491782677D+00, &
      0.00068639460491782677D+00, &
      0.00068639460491782677D+00, &
      0.00068639460491782677D+00, &
      0.00848846646315219895D+00, &
      0.00848846646315219895D+00, &
      0.00848846646315219895D+00, &
      0.00848846646315219895D+00, &
      0.00138976471399791693D+00, &
      0.00138976471399791693D+00, &
      0.00138976471399791693D+00, &
      0.00138976471399791693D+00, &
      0.00564687172166134736D+00, &
      0.00564687172166134736D+00, &
      0.00564687172166134736D+00, &
      0.00564687172166134736D+00, &
      0.01822286920356053219D+00, &
      0.01822286920356053219D+00, &
      0.01822286920356053219D+00, &
      0.01822286920356053219D+00, &
      0.00475553584446375757D+00, &
      0.00475553584446375757D+00, &
      0.00475553584446375757D+00, &
      0.00475553584446375757D+00, &
      0.00330839378673254258D+00, &
      0.00330839378673254258D+00, &
      0.00330839378673254258D+00, &
      0.00330839378673254258D+00, &
      0.00330839378673254258D+00, &
      0.00330839378673254258D+00, &
      0.00330839378673254258D+00, &
      0.00330839378673254258D+00, &
      0.00407004510316002176D+00, &
      0.00407004510316002176D+00, &
      0.00407004510316002176D+00, &
      0.00407004510316002176D+00, &
      0.00407004510316002176D+00, &
      0.00407004510316002176D+00, &
      0.00407004510316002176D+00, &
      0.00407004510316002176D+00, &
      0.00826848930450256592D+00, &
      0.00826848930450256592D+00, &
      0.00826848930450256592D+00, &
      0.00826848930450256592D+00, &
      0.00826848930450256592D+00, &
      0.00826848930450256592D+00, &
      0.00826848930450256592D+00, &
      0.00826848930450256592D+00, &
      0.00222137881102773540D+00, &
      0.00222137881102773540D+00, &
      0.00222137881102773540D+00, &
      0.00222137881102773540D+00, &
      0.00222137881102773540D+00, &
      0.00222137881102773540D+00, &
      0.00222137881102773540D+00, &
      0.00222137881102773540D+00, &
      0.00053414917220099471D+00, &
      0.00053414917220099471D+00, &
      0.00053414917220099471D+00, &
      0.00053414917220099471D+00, &
      0.00053414917220099471D+00, &
      0.00053414917220099471D+00, &
      0.00053414917220099471D+00, &
      0.00053414917220099471D+00, &
      0.00487946234737815943D+00, &
      0.00487946234737815943D+00, &
      0.00487946234737815943D+00, &
      0.00487946234737815943D+00, &
      0.00487946234737815943D+00, &
      0.00487946234737815943D+00, &
      0.00487946234737815943D+00, &
      0.00487946234737815943D+00, &
      0.00789729756886463005D+00, &
      0.00789729756886463005D+00, &
      0.00789729756886463005D+00, &
      0.00789729756886463005D+00, &
      0.00789729756886463005D+00, &
      0.00789729756886463005D+00, &
      0.00789729756886463005D+00, &
      0.00789729756886463005D+00, &
      0.00269662856611711495D+00, &
      0.00269662856611711495D+00, &
      0.00269662856611711495D+00, &
      0.00269662856611711495D+00, &
      0.00269662856611711495D+00, &
      0.00269662856611711495D+00, &
      0.00269662856611711495D+00, &
      0.00269662856611711495D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule15 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule15() returns the pyramid quadrature rule of precision 15.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 April 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 121, Number 11, pages 2418-2436, 15 June 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real ( kind = rk ) x[n], y[n], z[n]: the coordinates of
!    quadrature points.
!
!    real ( kind = rk ) w[n]: the quadrature weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 234

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.70014035742397184858D+00, &
      0.00000000000000000000D+00, &
     -0.70014035742397184858D+00, &
      0.00000000000000000000D+00, &
      0.91936314848339839578D+00, &
      0.00000000000000000000D+00, &
     -0.91936314848339839578D+00, &
      0.00000000000000000000D+00, &
      0.11013348894423610758D+00, &
      0.00000000000000000000D+00, &
     -0.11013348894423610758D+00, &
      0.00000000000000000000D+00, &
      0.44650085670154981976D+00, &
      0.00000000000000000000D+00, &
     -0.44650085670154981976D+00, &
      0.00000000000000000000D+00, &
      0.80205385758014868802D+00, &
      0.00000000000000000000D+00, &
     -0.80205385758014868802D+00, &
      0.00000000000000000000D+00, &
      0.53015124527324641868D+00, &
      0.00000000000000000000D+00, &
     -0.53015124527324641868D+00, &
      0.00000000000000000000D+00, &
      0.97322734428927526462D+00, &
      0.00000000000000000000D+00, &
     -0.97322734428927526462D+00, &
      0.00000000000000000000D+00, &
      0.50557186634064188446D+00, &
      0.00000000000000000000D+00, &
     -0.50557186634064188446D+00, &
      0.00000000000000000000D+00, &
      0.04618938260051261985D+00, &
      0.00000000000000000000D+00, &
     -0.04618938260051261985D+00, &
      0.00000000000000000000D+00, &
      0.62301364063828001960D+00, &
      0.00000000000000000000D+00, &
     -0.62301364063828001960D+00, &
      0.00000000000000000000D+00, &
      0.29509310458272342004D+00, &
      0.00000000000000000000D+00, &
     -0.29509310458272342004D+00, &
      0.00000000000000000000D+00, &
      0.47760058576112396356D+00, &
      0.00000000000000000000D+00, &
     -0.47760058576112396356D+00, &
      0.00000000000000000000D+00, &
      0.36898015805621120489D+00, &
      0.00000000000000000000D+00, &
     -0.36898015805621120489D+00, &
      0.00000000000000000000D+00, &
      0.37237829048305948199D+00, &
     -0.37237829048305948199D+00, &
      0.37237829048305948199D+00, &
     -0.37237829048305948199D+00, &
      0.53419020022883723087D+00, &
     -0.53419020022883723087D+00, &
      0.53419020022883723087D+00, &
     -0.53419020022883723087D+00, &
      0.06063808134784291065D+00, &
     -0.06063808134784291065D+00, &
      0.06063808134784291065D+00, &
     -0.06063808134784291065D+00, &
      0.47770687706587527943D+00, &
     -0.47770687706587527943D+00, &
      0.47770687706587527943D+00, &
     -0.47770687706587527943D+00, &
      0.19875431205737298379D+00, &
     -0.19875431205737298379D+00, &
      0.19875431205737298379D+00, &
     -0.19875431205737298379D+00, &
      0.75996956110375202265D+00, &
     -0.75996956110375202265D+00, &
      0.75996956110375202265D+00, &
     -0.75996956110375202265D+00, &
      0.50554584202664287762D+00, &
     -0.50554584202664287762D+00, &
      0.50554584202664287762D+00, &
     -0.50554584202664287762D+00, &
      0.22064824236265745405D+00, &
     -0.22064824236265745405D+00, &
      0.22064824236265745405D+00, &
     -0.22064824236265745405D+00, &
      0.34405034181784766023D+00, &
     -0.34405034181784766023D+00, &
      0.34405034181784766023D+00, &
     -0.34405034181784766023D+00, &
      0.19624129641007534430D+00, &
     -0.19624129641007534430D+00, &
      0.19624129641007534430D+00, &
     -0.19624129641007534430D+00, &
      0.91947901583738034237D+00, &
     -0.91947901583738034237D+00, &
      0.91947901583738034237D+00, &
     -0.91947901583738034237D+00, &
      0.26209258399888046842D+00, &
     -0.26209258399888046842D+00, &
      0.26209258399888046842D+00, &
     -0.26209258399888046842D+00, &
      0.17448879428144709047D+00, &
     -0.17448879428144709047D+00, &
      0.17448879428144709047D+00, &
     -0.17448879428144709047D+00, &
      0.13256324285803838814D+00, &
     -0.13256324285803838814D+00, &
      0.13256324285803838814D+00, &
     -0.13256324285803838814D+00, &
      0.30777107074463094794D+00, &
     -0.30777107074463094794D+00, &
      0.30777107074463094794D+00, &
     -0.30777107074463094794D+00, &
      0.45437881454568046502D+00, &
     -0.45437881454568046502D+00, &
      0.45437881454568046502D+00, &
     -0.45437881454568046502D+00, &
      0.18203893212319594008D+00, &
     -0.18203893212319594008D+00, &
      0.18203893212319594008D+00, &
     -0.18203893212319594008D+00, &
      0.18804430068288069400D+00, &
     -0.18804430068288069400D+00, &
      0.18804430068288069400D+00, &
     -0.18804430068288069400D+00, &
      0.85285866036694601977D+00, &
     -0.85285866036694601977D+00, &
      0.85285866036694601977D+00, &
     -0.85285866036694601977D+00, &
      0.19395358431436507396D+00, &
     -0.19395358431436507396D+00, &
      0.19395358431436507396D+00, &
     -0.19395358431436507396D+00, &
      0.63620569449244868121D+00, &
     -0.63620569449244868121D+00, &
      0.63620569449244868121D+00, &
     -0.63620569449244868121D+00, &
      0.71609629503236627013D+00, &
     -0.71609629503236627013D+00, &
      0.71609629503236627013D+00, &
     -0.71609629503236627013D+00, &
      0.72970588623268628492D+00, &
     -0.72970588623268628492D+00, &
      0.72970588623268628492D+00, &
     -0.72970588623268628492D+00, &
      0.85828508770542477624D+00, &
      0.57684038354992128728D+00, &
     -0.85828508770542477624D+00, &
      0.57684038354992128728D+00, &
      0.85828508770542477624D+00, &
     -0.57684038354992128728D+00, &
     -0.85828508770542477624D+00, &
     -0.57684038354992128728D+00, &
      0.76336766190340221705D+00, &
      0.21167809325271402798D+00, &
     -0.76336766190340221705D+00, &
      0.21167809325271402798D+00, &
      0.76336766190340221705D+00, &
     -0.21167809325271402798D+00, &
     -0.76336766190340221705D+00, &
     -0.21167809325271402798D+00, &
      0.21613035096350702302D+00, &
      0.58392235225835387169D+00, &
     -0.21613035096350702302D+00, &
      0.58392235225835387169D+00, &
      0.21613035096350702302D+00, &
     -0.58392235225835387169D+00, &
     -0.21613035096350702302D+00, &
     -0.58392235225835387169D+00, &
      0.95615842996744093707D+00, &
      0.64241827941687201786D+00, &
     -0.95615842996744093707D+00, &
      0.64241827941687201786D+00, &
      0.95615842996744093707D+00, &
     -0.64241827941687201786D+00, &
     -0.95615842996744093707D+00, &
     -0.64241827941687201786D+00, &
      0.87592219252594005763D+00, &
      0.37797565831285062643D+00, &
     -0.87592219252594005763D+00, &
      0.37797565831285062643D+00, &
      0.87592219252594005763D+00, &
     -0.37797565831285062643D+00, &
     -0.87592219252594005763D+00, &
     -0.37797565831285062643D+00, &
      0.24453940621889402873D+00, &
      0.46719416054032014696D+00, &
     -0.24453940621889402873D+00, &
      0.46719416054032014696D+00, &
      0.24453940621889402873D+00, &
     -0.46719416054032014696D+00, &
     -0.24453940621889402873D+00, &
     -0.46719416054032014696D+00, &
      0.37060589700094320742D+00, &
      0.54852690614804233693D+00, &
     -0.37060589700094320742D+00, &
      0.54852690614804233693D+00, &
      0.37060589700094320742D+00, &
     -0.54852690614804233693D+00, &
     -0.37060589700094320742D+00, &
     -0.54852690614804233693D+00, &
      0.61567982840573709513D+00, &
      0.00714017158789090767D+00, &
     -0.61567982840573709513D+00, &
      0.00714017158789090767D+00, &
      0.61567982840573709513D+00, &
     -0.00714017158789090767D+00, &
     -0.61567982840573709513D+00, &
     -0.00714017158789090767D+00, &
      0.71013273965512191399D+00, &
      0.27918941987191842058D+00, &
     -0.71013273965512191399D+00, &
      0.27918941987191842058D+00, &
      0.71013273965512191399D+00, &
     -0.27918941987191842058D+00, &
     -0.71013273965512191399D+00, &
     -0.27918941987191842058D+00, &
      0.82587110596774426785D+00, &
      0.35683026213317575737D+00, &
     -0.82587110596774426785D+00, &
      0.35683026213317575737D+00, &
      0.82587110596774426785D+00, &
     -0.35683026213317575737D+00, &
     -0.82587110596774426785D+00, &
     -0.35683026213317575737D+00, &
      0.47557663390649029811D+00, &
      0.69286495782397328203D+00, &
     -0.47557663390649029811D+00, &
      0.69286495782397328203D+00, &
      0.47557663390649029811D+00, &
     -0.69286495782397328203D+00, &
     -0.47557663390649029811D+00, &
     -0.69286495782397328203D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.70014035742397184858D+00, &
      0.00000000000000000000D+00, &
     -0.70014035742397184858D+00, &
      0.00000000000000000000D+00, &
      0.91936314848339839578D+00, &
      0.00000000000000000000D+00, &
     -0.91936314848339839578D+00, &
      0.00000000000000000000D+00, &
      0.11013348894423610758D+00, &
      0.00000000000000000000D+00, &
     -0.11013348894423610758D+00, &
      0.00000000000000000000D+00, &
      0.44650085670154981976D+00, &
      0.00000000000000000000D+00, &
     -0.44650085670154981976D+00, &
      0.00000000000000000000D+00, &
      0.80205385758014868802D+00, &
      0.00000000000000000000D+00, &
     -0.80205385758014868802D+00, &
      0.00000000000000000000D+00, &
      0.53015124527324641868D+00, &
      0.00000000000000000000D+00, &
     -0.53015124527324641868D+00, &
      0.00000000000000000000D+00, &
      0.97322734428927526462D+00, &
      0.00000000000000000000D+00, &
     -0.97322734428927526462D+00, &
      0.00000000000000000000D+00, &
      0.50557186634064188446D+00, &
      0.00000000000000000000D+00, &
     -0.50557186634064188446D+00, &
      0.00000000000000000000D+00, &
      0.04618938260051261985D+00, &
      0.00000000000000000000D+00, &
     -0.04618938260051261985D+00, &
      0.00000000000000000000D+00, &
      0.62301364063828001960D+00, &
      0.00000000000000000000D+00, &
     -0.62301364063828001960D+00, &
      0.00000000000000000000D+00, &
      0.29509310458272342004D+00, &
      0.00000000000000000000D+00, &
     -0.29509310458272342004D+00, &
      0.00000000000000000000D+00, &
      0.47760058576112396356D+00, &
      0.00000000000000000000D+00, &
     -0.47760058576112396356D+00, &
      0.00000000000000000000D+00, &
      0.36898015805621120489D+00, &
      0.00000000000000000000D+00, &
     -0.36898015805621120489D+00, &
      0.37237829048305948199D+00, &
      0.37237829048305948199D+00, &
     -0.37237829048305948199D+00, &
     -0.37237829048305948199D+00, &
      0.53419020022883723087D+00, &
      0.53419020022883723087D+00, &
     -0.53419020022883723087D+00, &
     -0.53419020022883723087D+00, &
      0.06063808134784291065D+00, &
      0.06063808134784291065D+00, &
     -0.06063808134784291065D+00, &
     -0.06063808134784291065D+00, &
      0.47770687706587527943D+00, &
      0.47770687706587527943D+00, &
     -0.47770687706587527943D+00, &
     -0.47770687706587527943D+00, &
      0.19875431205737298379D+00, &
      0.19875431205737298379D+00, &
     -0.19875431205737298379D+00, &
     -0.19875431205737298379D+00, &
      0.75996956110375202265D+00, &
      0.75996956110375202265D+00, &
     -0.75996956110375202265D+00, &
     -0.75996956110375202265D+00, &
      0.50554584202664287762D+00, &
      0.50554584202664287762D+00, &
     -0.50554584202664287762D+00, &
     -0.50554584202664287762D+00, &
      0.22064824236265745405D+00, &
      0.22064824236265745405D+00, &
     -0.22064824236265745405D+00, &
     -0.22064824236265745405D+00, &
      0.34405034181784766023D+00, &
      0.34405034181784766023D+00, &
     -0.34405034181784766023D+00, &
     -0.34405034181784766023D+00, &
      0.19624129641007534430D+00, &
      0.19624129641007534430D+00, &
     -0.19624129641007534430D+00, &
     -0.19624129641007534430D+00, &
      0.91947901583738034237D+00, &
      0.91947901583738034237D+00, &
     -0.91947901583738034237D+00, &
     -0.91947901583738034237D+00, &
      0.26209258399888046842D+00, &
      0.26209258399888046842D+00, &
     -0.26209258399888046842D+00, &
     -0.26209258399888046842D+00, &
      0.17448879428144709047D+00, &
      0.17448879428144709047D+00, &
     -0.17448879428144709047D+00, &
     -0.17448879428144709047D+00, &
      0.13256324285803838814D+00, &
      0.13256324285803838814D+00, &
     -0.13256324285803838814D+00, &
     -0.13256324285803838814D+00, &
      0.30777107074463094794D+00, &
      0.30777107074463094794D+00, &
     -0.30777107074463094794D+00, &
     -0.30777107074463094794D+00, &
      0.45437881454568046502D+00, &
      0.45437881454568046502D+00, &
     -0.45437881454568046502D+00, &
     -0.45437881454568046502D+00, &
      0.18203893212319594008D+00, &
      0.18203893212319594008D+00, &
     -0.18203893212319594008D+00, &
     -0.18203893212319594008D+00, &
      0.18804430068288069400D+00, &
      0.18804430068288069400D+00, &
     -0.18804430068288069400D+00, &
     -0.18804430068288069400D+00, &
      0.85285866036694601977D+00, &
      0.85285866036694601977D+00, &
     -0.85285866036694601977D+00, &
     -0.85285866036694601977D+00, &
      0.19395358431436507396D+00, &
      0.19395358431436507396D+00, &
     -0.19395358431436507396D+00, &
     -0.19395358431436507396D+00, &
      0.63620569449244868121D+00, &
      0.63620569449244868121D+00, &
     -0.63620569449244868121D+00, &
     -0.63620569449244868121D+00, &
      0.71609629503236627013D+00, &
      0.71609629503236627013D+00, &
     -0.71609629503236627013D+00, &
     -0.71609629503236627013D+00, &
      0.72970588623268628492D+00, &
      0.72970588623268628492D+00, &
     -0.72970588623268628492D+00, &
     -0.72970588623268628492D+00, &
      0.57684038354992128728D+00, &
      0.85828508770542477624D+00, &
      0.57684038354992128728D+00, &
     -0.85828508770542477624D+00, &
     -0.57684038354992128728D+00, &
      0.85828508770542477624D+00, &
     -0.57684038354992128728D+00, &
     -0.85828508770542477624D+00, &
      0.21167809325271402798D+00, &
      0.76336766190340221705D+00, &
      0.21167809325271402798D+00, &
     -0.76336766190340221705D+00, &
     -0.21167809325271402798D+00, &
      0.76336766190340221705D+00, &
     -0.21167809325271402798D+00, &
     -0.76336766190340221705D+00, &
      0.58392235225835387169D+00, &
      0.21613035096350702302D+00, &
      0.58392235225835387169D+00, &
     -0.21613035096350702302D+00, &
     -0.58392235225835387169D+00, &
      0.21613035096350702302D+00, &
     -0.58392235225835387169D+00, &
     -0.21613035096350702302D+00, &
      0.64241827941687201786D+00, &
      0.95615842996744093707D+00, &
      0.64241827941687201786D+00, &
     -0.95615842996744093707D+00, &
     -0.64241827941687201786D+00, &
      0.95615842996744093707D+00, &
     -0.64241827941687201786D+00, &
     -0.95615842996744093707D+00, &
      0.37797565831285062643D+00, &
      0.87592219252594005763D+00, &
      0.37797565831285062643D+00, &
     -0.87592219252594005763D+00, &
     -0.37797565831285062643D+00, &
      0.87592219252594005763D+00, &
     -0.37797565831285062643D+00, &
     -0.87592219252594005763D+00, &
      0.46719416054032014696D+00, &
      0.24453940621889402873D+00, &
      0.46719416054032014696D+00, &
     -0.24453940621889402873D+00, &
     -0.46719416054032014696D+00, &
      0.24453940621889402873D+00, &
     -0.46719416054032014696D+00, &
     -0.24453940621889402873D+00, &
      0.54852690614804233693D+00, &
      0.37060589700094320742D+00, &
      0.54852690614804233693D+00, &
     -0.37060589700094320742D+00, &
     -0.54852690614804233693D+00, &
      0.37060589700094320742D+00, &
     -0.54852690614804233693D+00, &
     -0.37060589700094320742D+00, &
      0.00714017158789090767D+00, &
      0.61567982840573709513D+00, &
      0.00714017158789090767D+00, &
     -0.61567982840573709513D+00, &
     -0.00714017158789090767D+00, &
      0.61567982840573709513D+00, &
     -0.00714017158789090767D+00, &
     -0.61567982840573709513D+00, &
      0.27918941987191842058D+00, &
      0.71013273965512191399D+00, &
      0.27918941987191842058D+00, &
     -0.71013273965512191399D+00, &
     -0.27918941987191842058D+00, &
      0.71013273965512191399D+00, &
     -0.27918941987191842058D+00, &
     -0.71013273965512191399D+00, &
      0.35683026213317575737D+00, &
      0.82587110596774426785D+00, &
      0.35683026213317575737D+00, &
     -0.82587110596774426785D+00, &
     -0.35683026213317575737D+00, &
      0.82587110596774426785D+00, &
     -0.35683026213317575737D+00, &
     -0.82587110596774426785D+00, &
      0.69286495782397328203D+00, &
      0.47557663390649029811D+00, &
      0.69286495782397328203D+00, &
     -0.47557663390649029811D+00, &
     -0.69286495782397328203D+00, &
      0.47557663390649029811D+00, &
     -0.69286495782397328203D+00, &
     -0.47557663390649029811D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
      0.95167840606392661851D+00, &
      0.29890092144247565331D+00, &
      0.10372522326204723642D+00, &
      0.10372522326204723642D+00, &
      0.10372522326204723642D+00, &
      0.10372522326204723642D+00, &
      0.06052669044962753070D+00, &
      0.06052669044962753070D+00, &
      0.06052669044962753070D+00, &
      0.06052669044962753070D+00, &
      0.87767072134508183900D+00, &
      0.87767072134508183900D+00, &
      0.87767072134508183900D+00, &
      0.87767072134508183900D+00, &
      0.33061431908844512995D+00, &
      0.33061431908844512995D+00, &
      0.33061431908844512995D+00, &
      0.33061431908844512995D+00, &
      0.18913057164406818500D+00, &
      0.18913057164406818500D+00, &
      0.18913057164406818500D+00, &
      0.18913057164406818500D+00, &
      0.15549793195551472880D+00, &
      0.15549793195551472880D+00, &
      0.15549793195551472880D+00, &
      0.15549793195551472880D+00, &
      0.00490511032939973824D+00, &
      0.00490511032939973824D+00, &
      0.00490511032939973824D+00, &
      0.00490511032939973824D+00, &
      0.06238803484862592841D+00, &
      0.06238803484862592841D+00, &
      0.06238803484862592841D+00, &
      0.06238803484862592841D+00, &
      0.53284146417918532013D+00, &
      0.53284146417918532013D+00, &
      0.53284146417918532013D+00, &
      0.53284146417918532013D+00, &
      0.37695593807084121218D+00, &
      0.37695593807084121218D+00, &
      0.37695593807084121218D+00, &
      0.37695593807084121218D+00, &
      0.68148251405524751245D+00, &
      0.68148251405524751245D+00, &
      0.68148251405524751245D+00, &
      0.68148251405524751245D+00, &
      0.23097535663334242684D+00, &
      0.23097535663334242684D+00, &
      0.23097535663334242684D+00, &
      0.23097535663334242684D+00, &
      0.50335666905043052743D+00, &
      0.50335666905043052743D+00, &
      0.50335666905043052743D+00, &
      0.50335666905043052743D+00, &
      0.42353675293216364039D+00, &
      0.42353675293216364039D+00, &
      0.42353675293216364039D+00, &
      0.42353675293216364039D+00, &
      0.41296679793793455993D+00, &
      0.41296679793793455993D+00, &
      0.41296679793793455993D+00, &
      0.41296679793793455993D+00, &
      0.78640579578998992538D+00, &
      0.78640579578998992538D+00, &
      0.78640579578998992538D+00, &
      0.78640579578998992538D+00, &
      0.32653775213007651956D+00, &
      0.32653775213007651956D+00, &
      0.32653775213007651956D+00, &
      0.32653775213007651956D+00, &
      0.05291423306214974864D+00, &
      0.05291423306214974864D+00, &
      0.05291423306214974864D+00, &
      0.05291423306214974864D+00, &
      0.00693099618919339969D+00, &
      0.00693099618919339969D+00, &
      0.00693099618919339969D+00, &
      0.00693099618919339969D+00, &
      0.01799129130112891647D+00, &
      0.01799129130112891647D+00, &
      0.01799129130112891647D+00, &
      0.01799129130112891647D+00, &
      0.00830477370838420340D+00, &
      0.00830477370838420340D+00, &
      0.00830477370838420340D+00, &
      0.00830477370838420340D+00, &
      0.61788739915987689333D+00, &
      0.61788739915987689333D+00, &
      0.61788739915987689333D+00, &
      0.61788739915987689333D+00, &
      0.19265498366863681445D+00, &
      0.19265498366863681445D+00, &
      0.19265498366863681445D+00, &
      0.19265498366863681445D+00, &
      0.01617257265834035757D+00, &
      0.01617257265834035757D+00, &
      0.01617257265834035757D+00, &
      0.01617257265834035757D+00, &
      0.59548810603233670591D+00, &
      0.59548810603233670591D+00, &
      0.59548810603233670591D+00, &
      0.59548810603233670591D+00, &
      0.76735204589989336466D+00, &
      0.76735204589989336466D+00, &
      0.76735204589989336466D+00, &
      0.76735204589989336466D+00, &
      0.64457148839122024864D+00, &
      0.64457148839122024864D+00, &
      0.64457148839122024864D+00, &
      0.64457148839122024864D+00, &
      0.26823389272477243805D+00, &
      0.26823389272477243805D+00, &
      0.26823389272477243805D+00, &
      0.26823389272477243805D+00, &
      0.19249663769902353172D+00, &
      0.19249663769902353172D+00, &
      0.19249663769902353172D+00, &
      0.19249663769902353172D+00, &
      0.35451205561493442930D+00, &
      0.35451205561493442930D+00, &
      0.35451205561493442930D+00, &
      0.35451205561493442930D+00, &
      0.13784397009929885702D+00, &
      0.13784397009929885702D+00, &
      0.13784397009929885702D+00, &
      0.13784397009929885702D+00, &
      0.08819642588395044946D+00, &
      0.08819642588395044946D+00, &
      0.08819642588395044946D+00, &
      0.08819642588395044946D+00, &
      0.44126341704052479686D+00, &
      0.44126341704052479686D+00, &
      0.44126341704052479686D+00, &
      0.44126341704052479686D+00, &
      0.18126706536632553046D+00, &
      0.18126706536632553046D+00, &
      0.18126706536632553046D+00, &
      0.18126706536632553046D+00, &
      0.06273138420630826328D+00, &
      0.06273138420630826328D+00, &
      0.06273138420630826328D+00, &
      0.06273138420630826328D+00, &
      0.22238800479499257201D+00, &
      0.22238800479499257201D+00, &
      0.22238800479499257201D+00, &
      0.22238800479499257201D+00, &
      0.12311329897268918909D+00, &
      0.12311329897268918909D+00, &
      0.12311329897268918909D+00, &
      0.12311329897268918909D+00, &
      0.12311329897268918909D+00, &
      0.12311329897268918909D+00, &
      0.12311329897268918909D+00, &
      0.12311329897268918909D+00, &
      0.02653658688975024313D+00, &
      0.02653658688975024313D+00, &
      0.02653658688975024313D+00, &
      0.02653658688975024313D+00, &
      0.02653658688975024313D+00, &
      0.02653658688975024313D+00, &
      0.02653658688975024313D+00, &
      0.02653658688975024313D+00, &
      0.32949702529445068500D+00, &
      0.32949702529445068500D+00, &
      0.32949702529445068500D+00, &
      0.32949702529445068500D+00, &
      0.32949702529445068500D+00, &
      0.32949702529445068500D+00, &
      0.32949702529445068500D+00, &
      0.32949702529445068500D+00, &
      0.02375108271645150898D+00, &
      0.02375108271645150898D+00, &
      0.02375108271645150898D+00, &
      0.02375108271645150898D+00, &
      0.02375108271645150898D+00, &
      0.02375108271645150898D+00, &
      0.02375108271645150898D+00, &
      0.02375108271645150898D+00, &
      0.01156156746353601689D+00, &
      0.01156156746353601689D+00, &
      0.01156156746353601689D+00, &
      0.01156156746353601689D+00, &
      0.01156156746353601689D+00, &
      0.01156156746353601689D+00, &
      0.01156156746353601689D+00, &
      0.01156156746353601689D+00, &
      0.50478124488277686943D+00, &
      0.50478124488277686943D+00, &
      0.50478124488277686943D+00, &
      0.50478124488277686943D+00, &
      0.50478124488277686943D+00, &
      0.50478124488277686943D+00, &
      0.50478124488277686943D+00, &
      0.50478124488277686943D+00, &
      0.09182493970820698737D+00, &
      0.09182493970820698737D+00, &
      0.09182493970820698737D+00, &
      0.09182493970820698737D+00, &
      0.09182493970820698737D+00, &
      0.09182493970820698737D+00, &
      0.09182493970820698737D+00, &
      0.09182493970820698737D+00, &
      0.01061641469954416848D+00, &
      0.01061641469954416848D+00, &
      0.01061641469954416848D+00, &
      0.01061641469954416848D+00, &
      0.01061641469954416848D+00, &
      0.01061641469954416848D+00, &
      0.01061641469954416848D+00, &
      0.01061641469954416848D+00, &
      0.18699248305545079774D+00, &
      0.18699248305545079774D+00, &
      0.18699248305545079774D+00, &
      0.18699248305545079774D+00, &
      0.18699248305545079774D+00, &
      0.18699248305545079774D+00, &
      0.18699248305545079774D+00, &
      0.18699248305545079774D+00, &
      0.07551964039133612916D+00, &
      0.07551964039133612916D+00, &
      0.07551964039133612916D+00, &
      0.07551964039133612916D+00, &
      0.07551964039133612916D+00, &
      0.07551964039133612916D+00, &
      0.07551964039133612916D+00, &
      0.07551964039133612916D+00, &
      0.28631001132523614672D+00, &
      0.28631001132523614672D+00, &
      0.28631001132523614672D+00, &
      0.28631001132523614672D+00, &
      0.28631001132523614672D+00, &
      0.28631001132523614672D+00, &
      0.28631001132523614672D+00, &
      0.28631001132523614672D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
      0.00042533015548600539D+00, &
      0.01525748792328103336D+00, &
      0.00528017662585705642D+00, &
      0.00528017662585705642D+00, &
      0.00528017662585705642D+00, &
      0.00528017662585705642D+00, &
      0.00244639458990403638D+00, &
      0.00244639458990403638D+00, &
      0.00244639458990403638D+00, &
      0.00244639458990403638D+00, &
      0.00097370603129723528D+00, &
      0.00097370603129723528D+00, &
      0.00097370603129723528D+00, &
      0.00097370603129723528D+00, &
      0.00777612495382315819D+00, &
      0.00777612495382315819D+00, &
      0.00777612495382315819D+00, &
      0.00777612495382315819D+00, &
      0.00267697092767919202D+00, &
      0.00267697092767919202D+00, &
      0.00267697092767919202D+00, &
      0.00267697092767919202D+00, &
      0.00619081284266621932D+00, &
      0.00619081284266621932D+00, &
      0.00619081284266621932D+00, &
      0.00619081284266621932D+00, &
      0.00065618656167288063D+00, &
      0.00065618656167288063D+00, &
      0.00065618656167288063D+00, &
      0.00065618656167288063D+00, &
      0.00656994056618269082D+00, &
      0.00656994056618269082D+00, &
      0.00656994056618269082D+00, &
      0.00656994056618269082D+00, &
      0.00267308743214515841D+00, &
      0.00267308743214515841D+00, &
      0.00267308743214515841D+00, &
      0.00267308743214515841D+00, &
      0.00175949530216048481D+00, &
      0.00175949530216048481D+00, &
      0.00175949530216048481D+00, &
      0.00175949530216048481D+00, &
      0.00331285760856183092D+00, &
      0.00331285760856183092D+00, &
      0.00331285760856183092D+00, &
      0.00331285760856183092D+00, &
      0.00751034506572093984D+00, &
      0.00751034506572093984D+00, &
      0.00751034506572093984D+00, &
      0.00751034506572093984D+00, &
      0.00893166529089345351D+00, &
      0.00893166529089345351D+00, &
      0.00893166529089345351D+00, &
      0.00893166529089345351D+00, &
      0.00556841825996671948D+00, &
      0.00556841825996671948D+00, &
      0.00556841825996671948D+00, &
      0.00556841825996671948D+00, &
      0.00264945739278132876D+00, &
      0.00264945739278132876D+00, &
      0.00264945739278132876D+00, &
      0.00264945739278132876D+00, &
      0.00197028999055473379D+00, &
      0.00197028999055473379D+00, &
      0.00197028999055473379D+00, &
      0.00197028999055473379D+00, &
      0.00513483903966016272D+00, &
      0.00513483903966016272D+00, &
      0.00513483903966016272D+00, &
      0.00513483903966016272D+00, &
      0.00742192548155929527D+00, &
      0.00742192548155929527D+00, &
      0.00742192548155929527D+00, &
      0.00742192548155929527D+00, &
      0.00164485721629859362D+00, &
      0.00164485721629859362D+00, &
      0.00164485721629859362D+00, &
      0.00164485721629859362D+00, &
      0.00553453693910099451D+00, &
      0.00553453693910099451D+00, &
      0.00553453693910099451D+00, &
      0.00553453693910099451D+00, &
      0.00325628156420546379D+00, &
      0.00325628156420546379D+00, &
      0.00325628156420546379D+00, &
      0.00325628156420546379D+00, &
      0.00151506082426114191D+00, &
      0.00151506082426114191D+00, &
      0.00151506082426114191D+00, &
      0.00151506082426114191D+00, &
      0.00591126556727756668D+00, &
      0.00591126556727756668D+00, &
      0.00591126556727756668D+00, &
      0.00591126556727756668D+00, &
      0.00083923797814166039D+00, &
      0.00083923797814166039D+00, &
      0.00083923797814166039D+00, &
      0.00083923797814166039D+00, &
      0.00402240401025226024D+00, &
      0.00402240401025226024D+00, &
      0.00402240401025226024D+00, &
      0.00402240401025226024D+00, &
      0.00223169575519892738D+00, &
      0.00223169575519892738D+00, &
      0.00223169575519892738D+00, &
      0.00223169575519892738D+00, &
      0.00615688200816193708D+00, &
      0.00615688200816193708D+00, &
      0.00615688200816193708D+00, &
      0.00615688200816193708D+00, &
      0.00864232058920194786D+00, &
      0.00864232058920194786D+00, &
      0.00864232058920194786D+00, &
      0.00864232058920194786D+00, &
      0.00745270338750553870D+00, &
      0.00745270338750553870D+00, &
      0.00745270338750553870D+00, &
      0.00745270338750553870D+00, &
      0.00502880383507017109D+00, &
      0.00502880383507017109D+00, &
      0.00502880383507017109D+00, &
      0.00502880383507017109D+00, &
      0.00879638512225545434D+00, &
      0.00879638512225545434D+00, &
      0.00879638512225545434D+00, &
      0.00879638512225545434D+00, &
      0.00174781246088704756D+00, &
      0.00174781246088704756D+00, &
      0.00174781246088704756D+00, &
      0.00174781246088704756D+00, &
      0.00918989383822925963D+00, &
      0.00918989383822925963D+00, &
      0.00918989383822925963D+00, &
      0.00918989383822925963D+00, &
      0.00620807391936788796D+00, &
      0.00620807391936788796D+00, &
      0.00620807391936788796D+00, &
      0.00620807391936788796D+00, &
      0.00568398055245949909D+00, &
      0.00568398055245949909D+00, &
      0.00568398055245949909D+00, &
      0.00568398055245949909D+00, &
      0.00183658571623211762D+00, &
      0.00183658571623211762D+00, &
      0.00183658571623211762D+00, &
      0.00183658571623211762D+00, &
      0.00221869264777755162D+00, &
      0.00221869264777755162D+00, &
      0.00221869264777755162D+00, &
      0.00221869264777755162D+00, &
      0.00221869264777755162D+00, &
      0.00221869264777755162D+00, &
      0.00221869264777755162D+00, &
      0.00221869264777755162D+00, &
      0.00260988740586402639D+00, &
      0.00260988740586402639D+00, &
      0.00260988740586402639D+00, &
      0.00260988740586402639D+00, &
      0.00260988740586402639D+00, &
      0.00260988740586402639D+00, &
      0.00260988740586402639D+00, &
      0.00260988740586402639D+00, &
      0.00616875610859360119D+00, &
      0.00616875610859360119D+00, &
      0.00616875610859360119D+00, &
      0.00616875610859360119D+00, &
      0.00616875610859360119D+00, &
      0.00616875610859360119D+00, &
      0.00616875610859360119D+00, &
      0.00616875610859360119D+00, &
      0.00122995761226799547D+00, &
      0.00122995761226799547D+00, &
      0.00122995761226799547D+00, &
      0.00122995761226799547D+00, &
      0.00122995761226799547D+00, &
      0.00122995761226799547D+00, &
      0.00122995761226799547D+00, &
      0.00122995761226799547D+00, &
      0.00195493878423380618D+00, &
      0.00195493878423380618D+00, &
      0.00195493878423380618D+00, &
      0.00195493878423380618D+00, &
      0.00195493878423380618D+00, &
      0.00195493878423380618D+00, &
      0.00195493878423380618D+00, &
      0.00195493878423380618D+00, &
      0.00346763933769407258D+00, &
      0.00346763933769407258D+00, &
      0.00346763933769407258D+00, &
      0.00346763933769407258D+00, &
      0.00346763933769407258D+00, &
      0.00346763933769407258D+00, &
      0.00346763933769407258D+00, &
      0.00346763933769407258D+00, &
      0.00677026604334400768D+00, &
      0.00677026604334400768D+00, &
      0.00677026604334400768D+00, &
      0.00677026604334400768D+00, &
      0.00677026604334400768D+00, &
      0.00677026604334400768D+00, &
      0.00677026604334400768D+00, &
      0.00677026604334400768D+00, &
      0.00167829416355688515D+00, &
      0.00167829416355688515D+00, &
      0.00167829416355688515D+00, &
      0.00167829416355688515D+00, &
      0.00167829416355688515D+00, &
      0.00167829416355688515D+00, &
      0.00167829416355688515D+00, &
      0.00167829416355688515D+00, &
      0.00713079759016688013D+00, &
      0.00713079759016688013D+00, &
      0.00713079759016688013D+00, &
      0.00713079759016688013D+00, &
      0.00713079759016688013D+00, &
      0.00713079759016688013D+00, &
      0.00713079759016688013D+00, &
      0.00713079759016688013D+00, &
      0.00467547289885053390D+00, &
      0.00467547289885053390D+00, &
      0.00467547289885053390D+00, &
      0.00467547289885053390D+00, &
      0.00467547289885053390D+00, &
      0.00467547289885053390D+00, &
      0.00467547289885053390D+00, &
      0.00467547289885053390D+00, &
      0.00253420752420773854D+00, &
      0.00253420752420773854D+00, &
      0.00253420752420773854D+00, &
      0.00253420752420773854D+00, &
      0.00253420752420773854D+00, &
      0.00253420752420773854D+00, &
      0.00253420752420773854D+00, &
      0.00253420752420773854D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule16 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule16() returns the pyramid quadrature rule of precision 16.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 April 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 121, Number 11, pages 2418-2436, 15 June 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real ( kind = rk ) x[n], y[n], z[n]: the coordinates of
!    quadrature points.
!
!    real ( kind = rk ) w[n]: the quadrature weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 285

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
      0.00000000000000000000D+00, &
      0.51870388022974789344D+00, &
      0.00000000000000000000D+00, &
     -0.51870388022974789344D+00, &
      0.00000000000000000000D+00, &
      0.91531292456766799592D+00, &
      0.00000000000000000000D+00, &
     -0.91531292456766799592D+00, &
      0.00000000000000000000D+00, &
      0.01655472656731888509D+00, &
      0.00000000000000000000D+00, &
     -0.01655472656731888509D+00, &
      0.00000000000000000000D+00, &
      0.36912135575078097727D+00, &
      0.00000000000000000000D+00, &
     -0.36912135575078097727D+00, &
      0.00000000000000000000D+00, &
      0.02060330211667108294D+00, &
      0.00000000000000000000D+00, &
     -0.02060330211667108294D+00, &
      0.00000000000000000000D+00, &
      0.87194193710104261896D+00, &
      0.00000000000000000000D+00, &
     -0.87194193710104261896D+00, &
      0.00000000000000000000D+00, &
      0.31538420787936460865D+00, &
      0.00000000000000000000D+00, &
     -0.31538420787936460865D+00, &
      0.00000000000000000000D+00, &
      0.61399555563153074278D+00, &
      0.00000000000000000000D+00, &
     -0.61399555563153074278D+00, &
      0.00000000000000000000D+00, &
      0.09765007634893953237D+00, &
      0.00000000000000000000D+00, &
     -0.09765007634893953237D+00, &
      0.00000000000000000000D+00, &
      0.73176464480784220168D+00, &
      0.00000000000000000000D+00, &
     -0.73176464480784220168D+00, &
      0.00000000000000000000D+00, &
      0.84054220943522905962D+00, &
      0.00000000000000000000D+00, &
     -0.84054220943522905962D+00, &
      0.00000000000000000000D+00, &
      0.55574609505629035677D+00, &
      0.00000000000000000000D+00, &
     -0.55574609505629035677D+00, &
      0.00000000000000000000D+00, &
      0.26438021717351822826D+00, &
      0.00000000000000000000D+00, &
     -0.26438021717351822826D+00, &
      0.00000000000000000000D+00, &
      0.13755031468645742554D+00, &
     -0.13755031468645742554D+00, &
      0.13755031468645742554D+00, &
     -0.13755031468645742554D+00, &
      0.55599929975376194413D+00, &
     -0.55599929975376194413D+00, &
      0.55599929975376194413D+00, &
     -0.55599929975376194413D+00, &
      0.68037602529077279012D+00, &
     -0.68037602529077279012D+00, &
      0.68037602529077279012D+00, &
     -0.68037602529077279012D+00, &
      0.68987074662201319786D+00, &
     -0.68987074662201319786D+00, &
      0.68987074662201319786D+00, &
     -0.68987074662201319786D+00, &
      0.46656016616381063011D+00, &
     -0.46656016616381063011D+00, &
      0.46656016616381063011D+00, &
     -0.46656016616381063011D+00, &
      0.93657427078716448676D+00, &
     -0.93657427078716448676D+00, &
      0.93657427078716448676D+00, &
     -0.93657427078716448676D+00, &
      0.10668464068510881415D+00, &
     -0.10668464068510881415D+00, &
      0.10668464068510881415D+00, &
     -0.10668464068510881415D+00, &
      0.56362866319030435758D+00, &
     -0.56362866319030435758D+00, &
      0.56362866319030435758D+00, &
     -0.56362866319030435758D+00, &
      0.08721470315267314255D+00, &
     -0.08721470315267314255D+00, &
      0.08721470315267314255D+00, &
     -0.08721470315267314255D+00, &
      0.28826482911237666373D+00, &
     -0.28826482911237666373D+00, &
      0.28826482911237666373D+00, &
     -0.28826482911237666373D+00, &
      0.23503249705404929970D+00, &
     -0.23503249705404929970D+00, &
      0.23503249705404929970D+00, &
     -0.23503249705404929970D+00, &
      0.24508582185951371946D+00, &
     -0.24508582185951371946D+00, &
      0.24508582185951371946D+00, &
     -0.24508582185951371946D+00, &
      0.50030649786215641850D+00, &
     -0.50030649786215641850D+00, &
      0.50030649786215641850D+00, &
     -0.50030649786215641850D+00, &
      0.05973130190403576345D+00, &
     -0.05973130190403576345D+00, &
      0.05973130190403576345D+00, &
     -0.05973130190403576345D+00, &
      0.36285963939710730308D+00, &
     -0.36285963939710730308D+00, &
      0.36285963939710730308D+00, &
     -0.36285963939710730308D+00, &
      0.83163011923015617288D+00, &
     -0.83163011923015617288D+00, &
      0.83163011923015617288D+00, &
     -0.83163011923015617288D+00, &
      0.59815151207893801910D+00, &
     -0.59815151207893801910D+00, &
      0.59815151207893801910D+00, &
     -0.59815151207893801910D+00, &
      0.32452468623385005708D+00, &
     -0.32452468623385005708D+00, &
      0.32452468623385005708D+00, &
     -0.32452468623385005708D+00, &
      0.18400437828694599096D+00, &
     -0.18400437828694599096D+00, &
      0.18400437828694599096D+00, &
     -0.18400437828694599096D+00, &
      0.80685585076678034699D+00, &
     -0.80685585076678034699D+00, &
      0.80685585076678034699D+00, &
     -0.80685585076678034699D+00, &
      0.39038172839426421579D+00, &
     -0.39038172839426421579D+00, &
      0.39038172839426421579D+00, &
     -0.39038172839426421579D+00, &
      0.16205656720464603482D+00, &
     -0.16205656720464603482D+00, &
      0.16205656720464603482D+00, &
     -0.16205656720464603482D+00, &
      0.79557791018541468286D+00, &
      0.52651319421754327887D+00, &
     -0.79557791018541468286D+00, &
      0.52651319421754327887D+00, &
      0.79557791018541468286D+00, &
     -0.52651319421754327887D+00, &
     -0.79557791018541468286D+00, &
     -0.52651319421754327887D+00, &
      0.82001165757451899285D+00, &
      0.58763897947885146422D+00, &
     -0.82001165757451899285D+00, &
      0.58763897947885146422D+00, &
      0.82001165757451899285D+00, &
     -0.58763897947885146422D+00, &
     -0.82001165757451899285D+00, &
     -0.58763897947885146422D+00, &
      0.48106775949054608743D+00, &
      0.14913701615330382522D+00, &
     -0.48106775949054608743D+00, &
      0.14913701615330382522D+00, &
      0.48106775949054608743D+00, &
     -0.14913701615330382522D+00, &
     -0.48106775949054608743D+00, &
     -0.14913701615330382522D+00, &
      0.91698017235976392314D+00, &
      0.63033943569284833774D+00, &
     -0.91698017235976392314D+00, &
      0.63033943569284833774D+00, &
      0.91698017235976392314D+00, &
     -0.63033943569284833774D+00, &
     -0.91698017235976392314D+00, &
     -0.63033943569284833774D+00, &
      0.95781996021629389748D+00, &
      0.33277831868591489783D+00, &
     -0.95781996021629389748D+00, &
      0.33277831868591489783D+00, &
      0.95781996021629389748D+00, &
     -0.33277831868591489783D+00, &
     -0.95781996021629389748D+00, &
     -0.33277831868591489783D+00, &
      0.21226792819096282350D+00, &
      0.42143671745401867224D+00, &
     -0.21226792819096282350D+00, &
      0.42143671745401867224D+00, &
      0.21226792819096282350D+00, &
     -0.42143671745401867224D+00, &
     -0.21226792819096282350D+00, &
     -0.42143671745401867224D+00, &
      0.24859573800356360440D+00, &
      0.62623217679932552393D+00, &
     -0.24859573800356360440D+00, &
      0.62623217679932552393D+00, &
      0.24859573800356360440D+00, &
     -0.62623217679932552393D+00, &
     -0.24859573800356360440D+00, &
     -0.62623217679932552393D+00, &
      0.55347901385730624568D+00, &
      0.18815805920102898763D+00, &
     -0.55347901385730624568D+00, &
      0.18815805920102898763D+00, &
      0.55347901385730624568D+00, &
     -0.18815805920102898763D+00, &
     -0.55347901385730624568D+00, &
     -0.18815805920102898763D+00, &
      0.82422827852456481690D+00, &
      0.36376205736175593053D+00, &
     -0.82422827852456481690D+00, &
      0.36376205736175593053D+00, &
      0.82422827852456481690D+00, &
     -0.36376205736175593053D+00, &
     -0.82422827852456481690D+00, &
     -0.36376205736175593053D+00, &
      0.14606542728759722150D+00, &
      0.44268310223761137001D+00, &
     -0.14606542728759722150D+00, &
      0.44268310223761137001D+00, &
      0.14606542728759722150D+00, &
     -0.44268310223761137001D+00, &
     -0.14606542728759722150D+00, &
     -0.44268310223761137001D+00, &
      0.38810720735022874450D+00, &
      0.62355125720878634699D+00, &
     -0.38810720735022874450D+00, &
      0.62355125720878634699D+00, &
      0.38810720735022874450D+00, &
     -0.62355125720878634699D+00, &
     -0.38810720735022874450D+00, &
     -0.62355125720878634699D+00, &
      0.95674484283290128772D+00, &
      0.79844517503668932523D+00, &
     -0.95674484283290128772D+00, &
      0.79844517503668932523D+00, &
      0.95674484283290128772D+00, &
     -0.79844517503668932523D+00, &
     -0.95674484283290128772D+00, &
     -0.79844517503668932523D+00, &
      0.60984849485865177954D+00, &
      0.31205150762818156807D+00, &
     -0.60984849485865177954D+00, &
      0.31205150762818156807D+00, &
      0.60984849485865177954D+00, &
     -0.31205150762818156807D+00, &
     -0.60984849485865177954D+00, &
     -0.31205150762818156807D+00, &
      0.24847427215400066935D+00, &
      0.13406531787546346890D+00, &
     -0.24847427215400066935D+00, &
      0.13406531787546346890D+00, &
      0.24847427215400066935D+00, &
     -0.13406531787546346890D+00, &
     -0.24847427215400066935D+00, &
     -0.13406531787546346890D+00, &
      0.69918943535356659069D+00, &
      0.17767876450874769967D+00, &
     -0.69918943535356659069D+00, &
      0.17767876450874769967D+00, &
      0.69918943535356659069D+00, &
     -0.17767876450874769967D+00, &
     -0.69918943535356659069D+00, &
     -0.17767876450874769967D+00, &
      0.36015415281659141078D+00, &
      0.43011776020837644285D+00, &
     -0.36015415281659141078D+00, &
      0.43011776020837644285D+00, &
      0.36015415281659141078D+00, &
     -0.43011776020837644285D+00, &
     -0.36015415281659141078D+00, &
     -0.43011776020837644285D+00, &
      0.73668742050512570074D+00, &
      0.30298735173674529175D+00, &
     -0.73668742050512570074D+00, &
      0.30298735173674529175D+00, &
      0.73668742050512570074D+00, &
     -0.30298735173674529175D+00, &
     -0.73668742050512570074D+00, &
     -0.30298735173674529175D+00, &
      0.31424806302114011158D+00, &
      0.02522065148266495332D+00, &
     -0.31424806302114011158D+00, &
      0.02522065148266495332D+00, &
      0.31424806302114011158D+00, &
     -0.02522065148266495332D+00, &
     -0.31424806302114011158D+00, &
     -0.02522065148266495332D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.51870388022974789344D+00, &
      0.00000000000000000000D+00, &
     -0.51870388022974789344D+00, &
      0.00000000000000000000D+00, &
      0.91531292456766799592D+00, &
      0.00000000000000000000D+00, &
     -0.91531292456766799592D+00, &
      0.00000000000000000000D+00, &
      0.01655472656731888509D+00, &
      0.00000000000000000000D+00, &
     -0.01655472656731888509D+00, &
      0.00000000000000000000D+00, &
      0.36912135575078097727D+00, &
      0.00000000000000000000D+00, &
     -0.36912135575078097727D+00, &
      0.00000000000000000000D+00, &
      0.02060330211667108294D+00, &
      0.00000000000000000000D+00, &
     -0.02060330211667108294D+00, &
      0.00000000000000000000D+00, &
      0.87194193710104261896D+00, &
      0.00000000000000000000D+00, &
     -0.87194193710104261896D+00, &
      0.00000000000000000000D+00, &
      0.31538420787936460865D+00, &
      0.00000000000000000000D+00, &
     -0.31538420787936460865D+00, &
      0.00000000000000000000D+00, &
      0.61399555563153074278D+00, &
      0.00000000000000000000D+00, &
     -0.61399555563153074278D+00, &
      0.00000000000000000000D+00, &
      0.09765007634893953237D+00, &
      0.00000000000000000000D+00, &
     -0.09765007634893953237D+00, &
      0.00000000000000000000D+00, &
      0.73176464480784220168D+00, &
      0.00000000000000000000D+00, &
     -0.73176464480784220168D+00, &
      0.00000000000000000000D+00, &
      0.84054220943522905962D+00, &
      0.00000000000000000000D+00, &
     -0.84054220943522905962D+00, &
      0.00000000000000000000D+00, &
      0.55574609505629035677D+00, &
      0.00000000000000000000D+00, &
     -0.55574609505629035677D+00, &
      0.00000000000000000000D+00, &
      0.26438021717351822826D+00, &
      0.00000000000000000000D+00, &
     -0.26438021717351822826D+00, &
      0.13755031468645742554D+00, &
      0.13755031468645742554D+00, &
     -0.13755031468645742554D+00, &
     -0.13755031468645742554D+00, &
      0.55599929975376194413D+00, &
      0.55599929975376194413D+00, &
     -0.55599929975376194413D+00, &
     -0.55599929975376194413D+00, &
      0.68037602529077279012D+00, &
      0.68037602529077279012D+00, &
     -0.68037602529077279012D+00, &
     -0.68037602529077279012D+00, &
      0.68987074662201319786D+00, &
      0.68987074662201319786D+00, &
     -0.68987074662201319786D+00, &
     -0.68987074662201319786D+00, &
      0.46656016616381063011D+00, &
      0.46656016616381063011D+00, &
     -0.46656016616381063011D+00, &
     -0.46656016616381063011D+00, &
      0.93657427078716448676D+00, &
      0.93657427078716448676D+00, &
     -0.93657427078716448676D+00, &
     -0.93657427078716448676D+00, &
      0.10668464068510881415D+00, &
      0.10668464068510881415D+00, &
     -0.10668464068510881415D+00, &
     -0.10668464068510881415D+00, &
      0.56362866319030435758D+00, &
      0.56362866319030435758D+00, &
     -0.56362866319030435758D+00, &
     -0.56362866319030435758D+00, &
      0.08721470315267314255D+00, &
      0.08721470315267314255D+00, &
     -0.08721470315267314255D+00, &
     -0.08721470315267314255D+00, &
      0.28826482911237666373D+00, &
      0.28826482911237666373D+00, &
     -0.28826482911237666373D+00, &
     -0.28826482911237666373D+00, &
      0.23503249705404929970D+00, &
      0.23503249705404929970D+00, &
     -0.23503249705404929970D+00, &
     -0.23503249705404929970D+00, &
      0.24508582185951371946D+00, &
      0.24508582185951371946D+00, &
     -0.24508582185951371946D+00, &
     -0.24508582185951371946D+00, &
      0.50030649786215641850D+00, &
      0.50030649786215641850D+00, &
     -0.50030649786215641850D+00, &
     -0.50030649786215641850D+00, &
      0.05973130190403576345D+00, &
      0.05973130190403576345D+00, &
     -0.05973130190403576345D+00, &
     -0.05973130190403576345D+00, &
      0.36285963939710730308D+00, &
      0.36285963939710730308D+00, &
     -0.36285963939710730308D+00, &
     -0.36285963939710730308D+00, &
      0.83163011923015617288D+00, &
      0.83163011923015617288D+00, &
     -0.83163011923015617288D+00, &
     -0.83163011923015617288D+00, &
      0.59815151207893801910D+00, &
      0.59815151207893801910D+00, &
     -0.59815151207893801910D+00, &
     -0.59815151207893801910D+00, &
      0.32452468623385005708D+00, &
      0.32452468623385005708D+00, &
     -0.32452468623385005708D+00, &
     -0.32452468623385005708D+00, &
      0.18400437828694599096D+00, &
      0.18400437828694599096D+00, &
     -0.18400437828694599096D+00, &
     -0.18400437828694599096D+00, &
      0.80685585076678034699D+00, &
      0.80685585076678034699D+00, &
     -0.80685585076678034699D+00, &
     -0.80685585076678034699D+00, &
      0.39038172839426421579D+00, &
      0.39038172839426421579D+00, &
     -0.39038172839426421579D+00, &
     -0.39038172839426421579D+00, &
      0.16205656720464603482D+00, &
      0.16205656720464603482D+00, &
     -0.16205656720464603482D+00, &
     -0.16205656720464603482D+00, &
      0.52651319421754327887D+00, &
      0.79557791018541468286D+00, &
      0.52651319421754327887D+00, &
     -0.79557791018541468286D+00, &
     -0.52651319421754327887D+00, &
      0.79557791018541468286D+00, &
     -0.52651319421754327887D+00, &
     -0.79557791018541468286D+00, &
      0.58763897947885146422D+00, &
      0.82001165757451899285D+00, &
      0.58763897947885146422D+00, &
     -0.82001165757451899285D+00, &
     -0.58763897947885146422D+00, &
      0.82001165757451899285D+00, &
     -0.58763897947885146422D+00, &
     -0.82001165757451899285D+00, &
      0.14913701615330382522D+00, &
      0.48106775949054608743D+00, &
      0.14913701615330382522D+00, &
     -0.48106775949054608743D+00, &
     -0.14913701615330382522D+00, &
      0.48106775949054608743D+00, &
     -0.14913701615330382522D+00, &
     -0.48106775949054608743D+00, &
      0.63033943569284833774D+00, &
      0.91698017235976392314D+00, &
      0.63033943569284833774D+00, &
     -0.91698017235976392314D+00, &
     -0.63033943569284833774D+00, &
      0.91698017235976392314D+00, &
     -0.63033943569284833774D+00, &
     -0.91698017235976392314D+00, &
      0.33277831868591489783D+00, &
      0.95781996021629389748D+00, &
      0.33277831868591489783D+00, &
     -0.95781996021629389748D+00, &
     -0.33277831868591489783D+00, &
      0.95781996021629389748D+00, &
     -0.33277831868591489783D+00, &
     -0.95781996021629389748D+00, &
      0.42143671745401867224D+00, &
      0.21226792819096282350D+00, &
      0.42143671745401867224D+00, &
     -0.21226792819096282350D+00, &
     -0.42143671745401867224D+00, &
      0.21226792819096282350D+00, &
     -0.42143671745401867224D+00, &
     -0.21226792819096282350D+00, &
      0.62623217679932552393D+00, &
      0.24859573800356360440D+00, &
      0.62623217679932552393D+00, &
     -0.24859573800356360440D+00, &
     -0.62623217679932552393D+00, &
      0.24859573800356360440D+00, &
     -0.62623217679932552393D+00, &
     -0.24859573800356360440D+00, &
      0.18815805920102898763D+00, &
      0.55347901385730624568D+00, &
      0.18815805920102898763D+00, &
     -0.55347901385730624568D+00, &
     -0.18815805920102898763D+00, &
      0.55347901385730624568D+00, &
     -0.18815805920102898763D+00, &
     -0.55347901385730624568D+00, &
      0.36376205736175593053D+00, &
      0.82422827852456481690D+00, &
      0.36376205736175593053D+00, &
     -0.82422827852456481690D+00, &
     -0.36376205736175593053D+00, &
      0.82422827852456481690D+00, &
     -0.36376205736175593053D+00, &
     -0.82422827852456481690D+00, &
      0.44268310223761137001D+00, &
      0.14606542728759722150D+00, &
      0.44268310223761137001D+00, &
     -0.14606542728759722150D+00, &
     -0.44268310223761137001D+00, &
      0.14606542728759722150D+00, &
     -0.44268310223761137001D+00, &
     -0.14606542728759722150D+00, &
      0.62355125720878634699D+00, &
      0.38810720735022874450D+00, &
      0.62355125720878634699D+00, &
     -0.38810720735022874450D+00, &
     -0.62355125720878634699D+00, &
      0.38810720735022874450D+00, &
     -0.62355125720878634699D+00, &
     -0.38810720735022874450D+00, &
      0.79844517503668932523D+00, &
      0.95674484283290128772D+00, &
      0.79844517503668932523D+00, &
     -0.95674484283290128772D+00, &
     -0.79844517503668932523D+00, &
      0.95674484283290128772D+00, &
     -0.79844517503668932523D+00, &
     -0.95674484283290128772D+00, &
      0.31205150762818156807D+00, &
      0.60984849485865177954D+00, &
      0.31205150762818156807D+00, &
     -0.60984849485865177954D+00, &
     -0.31205150762818156807D+00, &
      0.60984849485865177954D+00, &
     -0.31205150762818156807D+00, &
     -0.60984849485865177954D+00, &
      0.13406531787546346890D+00, &
      0.24847427215400066935D+00, &
      0.13406531787546346890D+00, &
     -0.24847427215400066935D+00, &
     -0.13406531787546346890D+00, &
      0.24847427215400066935D+00, &
     -0.13406531787546346890D+00, &
     -0.24847427215400066935D+00, &
      0.17767876450874769967D+00, &
      0.69918943535356659069D+00, &
      0.17767876450874769967D+00, &
     -0.69918943535356659069D+00, &
     -0.17767876450874769967D+00, &
      0.69918943535356659069D+00, &
     -0.17767876450874769967D+00, &
     -0.69918943535356659069D+00, &
      0.43011776020837644285D+00, &
      0.36015415281659141078D+00, &
      0.43011776020837644285D+00, &
     -0.36015415281659141078D+00, &
     -0.43011776020837644285D+00, &
      0.36015415281659141078D+00, &
     -0.43011776020837644285D+00, &
     -0.36015415281659141078D+00, &
      0.30298735173674529175D+00, &
      0.73668742050512570074D+00, &
      0.30298735173674529175D+00, &
     -0.73668742050512570074D+00, &
     -0.30298735173674529175D+00, &
      0.73668742050512570074D+00, &
     -0.30298735173674529175D+00, &
     -0.73668742050512570074D+00, &
      0.02522065148266495332D+00, &
      0.31424806302114011158D+00, &
      0.02522065148266495332D+00, &
     -0.31424806302114011158D+00, &
     -0.02522065148266495332D+00, &
      0.31424806302114011158D+00, &
     -0.02522065148266495332D+00, &
     -0.31424806302114011158D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
      0.95860354636208366941D+00, &
      0.14065559413657438559D+00, &
      0.14065559413657438559D+00, &
      0.14065559413657438559D+00, &
      0.14065559413657438559D+00, &
      0.04316491102619297859D+00, &
      0.04316491102619297859D+00, &
      0.04316491102619297859D+00, &
      0.04316491102619297859D+00, &
      0.39522832646013361657D+00, &
      0.39522832646013361657D+00, &
      0.39522832646013361657D+00, &
      0.39522832646013361657D+00, &
      0.29890125561593927639D+00, &
      0.29890125561593927639D+00, &
      0.29890125561593927639D+00, &
      0.29890125561593927639D+00, &
      0.06870461262955433746D+00, &
      0.06870461262955433746D+00, &
      0.06870461262955433746D+00, &
      0.06870461262955433746D+00, &
      0.11193495967860304929D+00, &
      0.11193495967860304929D+00, &
      0.11193495967860304929D+00, &
      0.11193495967860304929D+00, &
      0.56665435040432854397D+00, &
      0.56665435040432854397D+00, &
      0.56665435040432854397D+00, &
      0.56665435040432854397D+00, &
      0.15957732786354647536D+00, &
      0.15957732786354647536D+00, &
      0.15957732786354647536D+00, &
      0.15957732786354647536D+00, &
      0.89220214923004481644D+00, &
      0.89220214923004481644D+00, &
      0.89220214923004481644D+00, &
      0.89220214923004481644D+00, &
      0.24754993371512326594D+00, &
      0.24754993371512326594D+00, &
      0.24754993371512326594D+00, &
      0.24754993371512326594D+00, &
      0.01014565417563986251D+00, &
      0.01014565417563986251D+00, &
      0.01014565417563986251D+00, &
      0.01014565417563986251D+00, &
      0.42885512420774396514D+00, &
      0.42885512420774396514D+00, &
      0.42885512420774396514D+00, &
      0.42885512420774396514D+00, &
      0.71649796495360695836D+00, &
      0.71649796495360695836D+00, &
      0.71649796495360695836D+00, &
      0.71649796495360695836D+00, &
      0.17877523402579234557D+00, &
      0.17877523402579234557D+00, &
      0.17877523402579234557D+00, &
      0.17877523402579234557D+00, &
      0.29295501084965730465D+00, &
      0.29295501084965730465D+00, &
      0.29295501084965730465D+00, &
      0.29295501084965730465D+00, &
      0.27447523572941867620D+00, &
      0.27447523572941867620D+00, &
      0.27447523572941867620D+00, &
      0.27447523572941867620D+00, &
      0.15714127239405312197D+00, &
      0.15714127239405312197D+00, &
      0.15714127239405312197D+00, &
      0.15714127239405312197D+00, &
      0.25035924766443973244D+00, &
      0.25035924766443973244D+00, &
      0.25035924766443973244D+00, &
      0.25035924766443973244D+00, &
      0.02983913639912040561D+00, &
      0.02983913639912040561D+00, &
      0.02983913639912040561D+00, &
      0.02983913639912040561D+00, &
      0.70189657043113384827D+00, &
      0.70189657043113384827D+00, &
      0.70189657043113384827D+00, &
      0.70189657043113384827D+00, &
      0.05534643995880605266D+00, &
      0.05534643995880605266D+00, &
      0.05534643995880605266D+00, &
      0.05534643995880605266D+00, &
      0.58803180785949038523D+00, &
      0.58803180785949038523D+00, &
      0.58803180785949038523D+00, &
      0.58803180785949038523D+00, &
      0.35760040041294161028D+00, &
      0.35760040041294161028D+00, &
      0.35760040041294161028D+00, &
      0.35760040041294161028D+00, &
      0.63808499834946807994D+00, &
      0.63808499834946807994D+00, &
      0.63808499834946807994D+00, &
      0.63808499834946807994D+00, &
      0.11180046767985268863D+00, &
      0.11180046767985268863D+00, &
      0.11180046767985268863D+00, &
      0.11180046767985268863D+00, &
      0.45799604201643168144D+00, &
      0.45799604201643168144D+00, &
      0.45799604201643168144D+00, &
      0.45799604201643168144D+00, &
      0.81600939659279747573D+00, &
      0.81600939659279747573D+00, &
      0.81600939659279747573D+00, &
      0.81600939659279747573D+00, &
      0.05597294995869231404D+00, &
      0.05597294995869231404D+00, &
      0.05597294995869231404D+00, &
      0.05597294995869231404D+00, &
      0.12573106601315930941D+00, &
      0.12573106601315930941D+00, &
      0.12573106601315930941D+00, &
      0.12573106601315930941D+00, &
      0.10680199134354546875D+00, &
      0.10680199134354546875D+00, &
      0.10680199134354546875D+00, &
      0.10680199134354546875D+00, &
      0.64658620070888050968D+00, &
      0.64658620070888050968D+00, &
      0.64658620070888050968D+00, &
      0.64658620070888050968D+00, &
      0.25409417136336337473D+00, &
      0.25409417136336337473D+00, &
      0.25409417136336337473D+00, &
      0.25409417136336337473D+00, &
      0.05465918013108818363D+00, &
      0.05465918013108818363D+00, &
      0.05465918013108818363D+00, &
      0.05465918013108818363D+00, &
      0.45345939462310336232D+00, &
      0.45345939462310336232D+00, &
      0.45345939462310336232D+00, &
      0.45345939462310336232D+00, &
      0.78917660780490361816D+00, &
      0.78917660780490361816D+00, &
      0.78917660780490361816D+00, &
      0.78917660780490361816D+00, &
      0.18689788541072827055D+00, &
      0.18689788541072827055D+00, &
      0.18689788541072827055D+00, &
      0.18689788541072827055D+00, &
      0.18689788541072827055D+00, &
      0.18689788541072827055D+00, &
      0.18689788541072827055D+00, &
      0.18689788541072827055D+00, &
      0.01494308906905271982D+00, &
      0.01494308906905271982D+00, &
      0.01494308906905271982D+00, &
      0.01494308906905271982D+00, &
      0.01494308906905271982D+00, &
      0.01494308906905271982D+00, &
      0.01494308906905271982D+00, &
      0.01494308906905271982D+00, &
      0.39770246065633879651D+00, &
      0.39770246065633879651D+00, &
      0.39770246065633879651D+00, &
      0.39770246065633879651D+00, &
      0.39770246065633879651D+00, &
      0.39770246065633879651D+00, &
      0.39770246065633879651D+00, &
      0.39770246065633879651D+00, &
      0.06620609897561705037D+00, &
      0.06620609897561705037D+00, &
      0.06620609897561705037D+00, &
      0.06620609897561705037D+00, &
      0.06620609897561705037D+00, &
      0.06620609897561705037D+00, &
      0.06620609897561705037D+00, &
      0.06620609897561705037D+00, &
      0.01432824105432797292D+00, &
      0.01432824105432797292D+00, &
      0.01432824105432797292D+00, &
      0.01432824105432797292D+00, &
      0.01432824105432797292D+00, &
      0.01432824105432797292D+00, &
      0.01432824105432797292D+00, &
      0.01432824105432797292D+00, &
      0.55099638203681522430D+00, &
      0.55099638203681522430D+00, &
      0.55099638203681522430D+00, &
      0.55099638203681522430D+00, &
      0.55099638203681522430D+00, &
      0.55099638203681522430D+00, &
      0.55099638203681522430D+00, &
      0.55099638203681522430D+00, &
      0.28000300498314673048D+00, &
      0.28000300498314673048D+00, &
      0.28000300498314673048D+00, &
      0.28000300498314673048D+00, &
      0.28000300498314673048D+00, &
      0.28000300498314673048D+00, &
      0.28000300498314673048D+00, &
      0.28000300498314673048D+00, &
      0.24363304655846618196D+00, &
      0.24363304655846618196D+00, &
      0.24363304655846618196D+00, &
      0.24363304655846618196D+00, &
      0.24363304655846618196D+00, &
      0.24363304655846618196D+00, &
      0.24363304655846618196D+00, &
      0.24363304655846618196D+00, &
      0.07506931441576017439D+00, &
      0.07506931441576017439D+00, &
      0.07506931441576017439D+00, &
      0.07506931441576017439D+00, &
      0.07506931441576017439D+00, &
      0.07506931441576017439D+00, &
      0.07506931441576017439D+00, &
      0.07506931441576017439D+00, &
      0.05795059845476016602D+00, &
      0.05795059845476016602D+00, &
      0.05795059845476016602D+00, &
      0.05795059845476016602D+00, &
      0.05795059845476016602D+00, &
      0.05795059845476016602D+00, &
      0.05795059845476016602D+00, &
      0.05795059845476016602D+00, &
      0.35631158631326498298D+00, &
      0.35631158631326498298D+00, &
      0.35631158631326498298D+00, &
      0.35631158631326498298D+00, &
      0.35631158631326498298D+00, &
      0.35631158631326498298D+00, &
      0.35631158631326498298D+00, &
      0.35631158631326498298D+00, &
      0.00597525511512164969D+00, &
      0.00597525511512164969D+00, &
      0.00597525511512164969D+00, &
      0.00597525511512164969D+00, &
      0.00597525511512164969D+00, &
      0.00597525511512164969D+00, &
      0.00597525511512164969D+00, &
      0.00597525511512164969D+00, &
      0.01006650964023972882D+00, &
      0.01006650964023972882D+00, &
      0.01006650964023972882D+00, &
      0.01006650964023972882D+00, &
      0.01006650964023972882D+00, &
      0.01006650964023972882D+00, &
      0.01006650964023972882D+00, &
      0.01006650964023972882D+00, &
      0.47619683731341805322D+00, &
      0.47619683731341805322D+00, &
      0.47619683731341805322D+00, &
      0.47619683731341805322D+00, &
      0.47619683731341805322D+00, &
      0.47619683731341805322D+00, &
      0.47619683731341805322D+00, &
      0.47619683731341805322D+00, &
      0.05836635447459245785D+00, &
      0.05836635447459245785D+00, &
      0.05836635447459245785D+00, &
      0.05836635447459245785D+00, &
      0.05836635447459245785D+00, &
      0.05836635447459245785D+00, &
      0.05836635447459245785D+00, &
      0.05836635447459245785D+00, &
      0.16495562918454476087D+00, &
      0.16495562918454476087D+00, &
      0.16495562918454476087D+00, &
      0.16495562918454476087D+00, &
      0.16495562918454476087D+00, &
      0.16495562918454476087D+00, &
      0.16495562918454476087D+00, &
      0.16495562918454476087D+00, &
      0.15963984721353552398D+00, &
      0.15963984721353552398D+00, &
      0.15963984721353552398D+00, &
      0.15963984721353552398D+00, &
      0.15963984721353552398D+00, &
      0.15963984721353552398D+00, &
      0.15963984721353552398D+00, &
      0.15963984721353552398D+00, &
      0.01331460692424766938D+00, &
      0.01331460692424766938D+00, &
      0.01331460692424766938D+00, &
      0.01331460692424766938D+00, &
      0.01331460692424766938D+00, &
      0.01331460692424766938D+00, &
      0.01331460692424766938D+00, &
      0.01331460692424766938D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
      0.00027672120207682897D+00, &
      0.00734944257707945050D+00, &
      0.00734944257707945050D+00, &
      0.00734944257707945050D+00, &
      0.00734944257707945050D+00, &
      0.00148740437072881278D+00, &
      0.00148740437072881278D+00, &
      0.00148740437072881278D+00, &
      0.00148740437072881278D+00, &
      0.00387693638421187639D+00, &
      0.00387693638421187639D+00, &
      0.00387693638421187639D+00, &
      0.00387693638421187639D+00, &
      0.01030614824039441553D+00, &
      0.01030614824039441553D+00, &
      0.01030614824039441553D+00, &
      0.01030614824039441553D+00, &
      0.00272673198843895789D+00, &
      0.00272673198843895789D+00, &
      0.00272673198843895789D+00, &
      0.00272673198843895789D+00, &
      0.00244629213408558690D+00, &
      0.00244629213408558690D+00, &
      0.00244629213408558690D+00, &
      0.00244629213408558690D+00, &
      0.00554924399761744409D+00, &
      0.00554924399761744409D+00, &
      0.00554924399761744409D+00, &
      0.00554924399761744409D+00, &
      0.00502262417466419701D+00, &
      0.00502262417466419701D+00, &
      0.00502262417466419701D+00, &
      0.00502262417466419701D+00, &
      0.00067111170084793834D+00, &
      0.00067111170084793834D+00, &
      0.00067111170084793834D+00, &
      0.00067111170084793834D+00, &
      0.00316873518673754882D+00, &
      0.00316873518673754882D+00, &
      0.00316873518673754882D+00, &
      0.00316873518673754882D+00, &
      0.00230707268685653678D+00, &
      0.00230707268685653678D+00, &
      0.00230707268685653678D+00, &
      0.00230707268685653678D+00, &
      0.00223355434035575010D+00, &
      0.00223355434035575010D+00, &
      0.00223355434035575010D+00, &
      0.00223355434035575010D+00, &
      0.00228685693132126217D+00, &
      0.00228685693132126217D+00, &
      0.00228685693132126217D+00, &
      0.00228685693132126217D+00, &
      0.00583888699964928409D+00, &
      0.00583888699964928409D+00, &
      0.00583888699964928409D+00, &
      0.00583888699964928409D+00, &
      0.00412038756118528714D+00, &
      0.00412038756118528714D+00, &
      0.00412038756118528714D+00, &
      0.00412038756118528714D+00, &
      0.00173037948168740284D+00, &
      0.00173037948168740284D+00, &
      0.00173037948168740284D+00, &
      0.00173037948168740284D+00, &
      0.00410680906602167630D+00, &
      0.00410680906602167630D+00, &
      0.00410680906602167630D+00, &
      0.00410680906602167630D+00, &
      0.00635438948067259222D+00, &
      0.00635438948067259222D+00, &
      0.00635438948067259222D+00, &
      0.00635438948067259222D+00, &
      0.00043262101844333546D+00, &
      0.00043262101844333546D+00, &
      0.00043262101844333546D+00, &
      0.00043262101844333546D+00, &
      0.00345263292981801964D+00, &
      0.00345263292981801964D+00, &
      0.00345263292981801964D+00, &
      0.00345263292981801964D+00, &
      0.00424071455736505301D+00, &
      0.00424071455736505301D+00, &
      0.00424071455736505301D+00, &
      0.00424071455736505301D+00, &
      0.00340526631678994204D+00, &
      0.00340526631678994204D+00, &
      0.00340526631678994204D+00, &
      0.00340526631678994204D+00, &
      0.00700942929302337137D+00, &
      0.00700942929302337137D+00, &
      0.00700942929302337137D+00, &
      0.00700942929302337137D+00, &
      0.00371184509938299149D+00, &
      0.00371184509938299149D+00, &
      0.00371184509938299149D+00, &
      0.00371184509938299149D+00, &
      0.00655769989421819445D+00, &
      0.00655769989421819445D+00, &
      0.00655769989421819445D+00, &
      0.00655769989421819445D+00, &
      0.00169996297833195664D+00, &
      0.00169996297833195664D+00, &
      0.00169996297833195664D+00, &
      0.00169996297833195664D+00, &
      0.00134740004835729398D+00, &
      0.00134740004835729398D+00, &
      0.00134740004835729398D+00, &
      0.00134740004835729398D+00, &
      0.00363396950461261376D+00, &
      0.00363396950461261376D+00, &
      0.00363396950461261376D+00, &
      0.00363396950461261376D+00, &
      0.00121619462900863910D+00, &
      0.00121619462900863910D+00, &
      0.00121619462900863910D+00, &
      0.00121619462900863910D+00, &
      0.00549579708316756677D+00, &
      0.00549579708316756677D+00, &
      0.00549579708316756677D+00, &
      0.00549579708316756677D+00, &
      0.00097330273562793888D+00, &
      0.00097330273562793888D+00, &
      0.00097330273562793888D+00, &
      0.00097330273562793888D+00, &
      0.00725393003232855348D+00, &
      0.00725393003232855348D+00, &
      0.00725393003232855348D+00, &
      0.00725393003232855348D+00, &
      0.00261600067059077025D+00, &
      0.00261600067059077025D+00, &
      0.00261600067059077025D+00, &
      0.00261600067059077025D+00, &
      0.00548333692337083910D+00, &
      0.00548333692337083910D+00, &
      0.00548333692337083910D+00, &
      0.00548333692337083910D+00, &
      0.00150187064516151079D+00, &
      0.00150187064516151079D+00, &
      0.00150187064516151079D+00, &
      0.00150187064516151079D+00, &
      0.00224211276351870626D+00, &
      0.00224211276351870626D+00, &
      0.00224211276351870626D+00, &
      0.00224211276351870626D+00, &
      0.00224211276351870626D+00, &
      0.00224211276351870626D+00, &
      0.00224211276351870626D+00, &
      0.00224211276351870626D+00, &
      0.00254357012140252583D+00, &
      0.00254357012140252583D+00, &
      0.00254357012140252583D+00, &
      0.00254357012140252583D+00, &
      0.00254357012140252583D+00, &
      0.00254357012140252583D+00, &
      0.00254357012140252583D+00, &
      0.00254357012140252583D+00, &
      0.00536144481406866224D+00, &
      0.00536144481406866224D+00, &
      0.00536144481406866224D+00, &
      0.00536144481406866224D+00, &
      0.00536144481406866224D+00, &
      0.00536144481406866224D+00, &
      0.00536144481406866224D+00, &
      0.00536144481406866224D+00, &
      0.00150566238487620605D+00, &
      0.00150566238487620605D+00, &
      0.00150566238487620605D+00, &
      0.00150566238487620605D+00, &
      0.00150566238487620605D+00, &
      0.00150566238487620605D+00, &
      0.00150566238487620605D+00, &
      0.00150566238487620605D+00, &
      0.00112651998000158019D+00, &
      0.00112651998000158019D+00, &
      0.00112651998000158019D+00, &
      0.00112651998000158019D+00, &
      0.00112651998000158019D+00, &
      0.00112651998000158019D+00, &
      0.00112651998000158019D+00, &
      0.00112651998000158019D+00, &
      0.00283697215367184055D+00, &
      0.00283697215367184055D+00, &
      0.00283697215367184055D+00, &
      0.00283697215367184055D+00, &
      0.00283697215367184055D+00, &
      0.00283697215367184055D+00, &
      0.00283697215367184055D+00, &
      0.00283697215367184055D+00, &
      0.00423443304160556894D+00, &
      0.00423443304160556894D+00, &
      0.00423443304160556894D+00, &
      0.00423443304160556894D+00, &
      0.00423443304160556894D+00, &
      0.00423443304160556894D+00, &
      0.00423443304160556894D+00, &
      0.00423443304160556894D+00, &
      0.00318414785741326324D+00, &
      0.00318414785741326324D+00, &
      0.00318414785741326324D+00, &
      0.00318414785741326324D+00, &
      0.00318414785741326324D+00, &
      0.00318414785741326324D+00, &
      0.00318414785741326324D+00, &
      0.00318414785741326324D+00, &
      0.00402977998773114843D+00, &
      0.00402977998773114843D+00, &
      0.00402977998773114843D+00, &
      0.00402977998773114843D+00, &
      0.00402977998773114843D+00, &
      0.00402977998773114843D+00, &
      0.00402977998773114843D+00, &
      0.00402977998773114843D+00, &
      0.00366973681687540396D+00, &
      0.00366973681687540396D+00, &
      0.00366973681687540396D+00, &
      0.00366973681687540396D+00, &
      0.00366973681687540396D+00, &
      0.00366973681687540396D+00, &
      0.00366973681687540396D+00, &
      0.00366973681687540396D+00, &
      0.00248888728757531682D+00, &
      0.00248888728757531682D+00, &
      0.00248888728757531682D+00, &
      0.00248888728757531682D+00, &
      0.00248888728757531682D+00, &
      0.00248888728757531682D+00, &
      0.00248888728757531682D+00, &
      0.00248888728757531682D+00, &
      0.00048888922605729013D+00, &
      0.00048888922605729013D+00, &
      0.00048888922605729013D+00, &
      0.00048888922605729013D+00, &
      0.00048888922605729013D+00, &
      0.00048888922605729013D+00, &
      0.00048888922605729013D+00, &
      0.00048888922605729013D+00, &
      0.00325264515509947555D+00, &
      0.00325264515509947555D+00, &
      0.00325264515509947555D+00, &
      0.00325264515509947555D+00, &
      0.00325264515509947555D+00, &
      0.00325264515509947555D+00, &
      0.00325264515509947555D+00, &
      0.00325264515509947555D+00, &
      0.00495505909145715732D+00, &
      0.00495505909145715732D+00, &
      0.00495505909145715732D+00, &
      0.00495505909145715732D+00, &
      0.00495505909145715732D+00, &
      0.00495505909145715732D+00, &
      0.00495505909145715732D+00, &
      0.00495505909145715732D+00, &
      0.00457024218025166729D+00, &
      0.00457024218025166729D+00, &
      0.00457024218025166729D+00, &
      0.00457024218025166729D+00, &
      0.00457024218025166729D+00, &
      0.00457024218025166729D+00, &
      0.00457024218025166729D+00, &
      0.00457024218025166729D+00, &
      0.00481061764241353957D+00, &
      0.00481061764241353957D+00, &
      0.00481061764241353957D+00, &
      0.00481061764241353957D+00, &
      0.00481061764241353957D+00, &
      0.00481061764241353957D+00, &
      0.00481061764241353957D+00, &
      0.00481061764241353957D+00, &
      0.00538791794721728626D+00, &
      0.00538791794721728626D+00, &
      0.00538791794721728626D+00, &
      0.00538791794721728626D+00, &
      0.00538791794721728626D+00, &
      0.00538791794721728626D+00, &
      0.00538791794721728626D+00, &
      0.00538791794721728626D+00, &
      0.00246928056742645339D+00, &
      0.00246928056742645339D+00, &
      0.00246928056742645339D+00, &
      0.00246928056742645339D+00, &
      0.00246928056742645339D+00, &
      0.00246928056742645339D+00, &
      0.00246928056742645339D+00, &
      0.00246928056742645339D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule17 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule17() returns the pyramid quadrature rule of precision 17.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 April 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 121, Number 11, pages 2418-2436, 15 June 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real ( kind = rk ) x[n], y[n], z[n]: the coordinates of
!    quadrature points.
!
!    real ( kind = rk ) w[n]: the quadrature weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 319

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.95241408581244602072D+00, &
      0.00000000000000000000D+00, &
     -0.95241408581244602072D+00, &
      0.00000000000000000000D+00, &
      0.08275692125030435775D+00, &
      0.00000000000000000000D+00, &
     -0.08275692125030435775D+00, &
      0.00000000000000000000D+00, &
      0.26578545043127727032D+00, &
      0.00000000000000000000D+00, &
     -0.26578545043127727032D+00, &
      0.00000000000000000000D+00, &
      0.81155940028479511827D+00, &
      0.00000000000000000000D+00, &
     -0.81155940028479511827D+00, &
      0.00000000000000000000D+00, &
      0.35232799954016802424D+00, &
      0.00000000000000000000D+00, &
     -0.35232799954016802424D+00, &
      0.00000000000000000000D+00, &
      0.60303771218972423984D+00, &
      0.00000000000000000000D+00, &
     -0.60303771218972423984D+00, &
      0.00000000000000000000D+00, &
      0.78402604202159054125D+00, &
      0.00000000000000000000D+00, &
     -0.78402604202159054125D+00, &
      0.00000000000000000000D+00, &
      0.49346472198646390561D+00, &
      0.00000000000000000000D+00, &
     -0.49346472198646390561D+00, &
      0.00000000000000000000D+00, &
      0.35038751648906635294D+00, &
      0.00000000000000000000D+00, &
     -0.35038751648906635294D+00, &
      0.00000000000000000000D+00, &
      0.00336404312107426752D+00, &
      0.00000000000000000000D+00, &
     -0.00336404312107426752D+00, &
      0.00000000000000000000D+00, &
      0.43074512679067450405D+00, &
      0.00000000000000000000D+00, &
     -0.43074512679067450405D+00, &
      0.00000000000000000000D+00, &
      0.59254801282947688890D+00, &
      0.00000000000000000000D+00, &
     -0.59254801282947688890D+00, &
      0.00000000000000000000D+00, &
      0.22550686036585329552D+00, &
      0.00000000000000000000D+00, &
     -0.22550686036585329552D+00, &
      0.00000000000000000000D+00, &
      0.63975863740336180729D+00, &
     -0.63975863740336180729D+00, &
      0.63975863740336180729D+00, &
     -0.63975863740336180729D+00, &
      0.16463474879631256886D+00, &
     -0.16463474879631256886D+00, &
      0.16463474879631256886D+00, &
     -0.16463474879631256886D+00, &
      0.49303563319629611916D+00, &
     -0.49303563319629611916D+00, &
      0.49303563319629611916D+00, &
     -0.49303563319629611916D+00, &
      0.74791280396309745004D+00, &
     -0.74791280396309745004D+00, &
      0.74791280396309745004D+00, &
     -0.74791280396309745004D+00, &
      0.26251317864285966808D+00, &
     -0.26251317864285966808D+00, &
      0.26251317864285966808D+00, &
     -0.26251317864285966808D+00, &
      0.17484215589868074003D+00, &
     -0.17484215589868074003D+00, &
      0.17484215589868074003D+00, &
     -0.17484215589868074003D+00, &
      0.27602468785670569718D+00, &
     -0.27602468785670569718D+00, &
      0.27602468785670569718D+00, &
     -0.27602468785670569718D+00, &
      0.43512552408987620334D+00, &
     -0.43512552408987620334D+00, &
      0.43512552408987620334D+00, &
     -0.43512552408987620334D+00, &
      0.95307694157738065410D+00, &
     -0.95307694157738065410D+00, &
      0.95307694157738065410D+00, &
     -0.95307694157738065410D+00, &
      0.48425597905194911474D+00, &
     -0.48425597905194911474D+00, &
      0.48425597905194911474D+00, &
     -0.48425597905194911474D+00, &
      0.29381061843327571648D+00, &
     -0.29381061843327571648D+00, &
      0.29381061843327571648D+00, &
     -0.29381061843327571648D+00, &
      0.16057200503692023452D+00, &
     -0.16057200503692023452D+00, &
      0.16057200503692023452D+00, &
     -0.16057200503692023452D+00, &
      0.81014529214207364749D+00, &
     -0.81014529214207364749D+00, &
      0.81014529214207364749D+00, &
     -0.81014529214207364749D+00, &
      0.18009325754915589402D+00, &
     -0.18009325754915589402D+00, &
      0.18009325754915589402D+00, &
     -0.18009325754915589402D+00, &
      0.50136799407459931022D+00, &
     -0.50136799407459931022D+00, &
      0.50136799407459931022D+00, &
     -0.50136799407459931022D+00, &
      0.92284407855436589863D+00, &
     -0.92284407855436589863D+00, &
      0.92284407855436589863D+00, &
     -0.92284407855436589863D+00, &
      0.32077249570607596629D+00, &
     -0.32077249570607596629D+00, &
      0.32077249570607596629D+00, &
     -0.32077249570607596629D+00, &
      0.02145401504974505866D+00, &
     -0.02145401504974505866D+00, &
      0.02145401504974505866D+00, &
     -0.02145401504974505866D+00, &
      0.65293078741937671250D+00, &
     -0.65293078741937671250D+00, &
      0.65293078741937671250D+00, &
     -0.65293078741937671250D+00, &
      0.38109206736816036987D+00, &
     -0.38109206736816036987D+00, &
      0.38109206736816036987D+00, &
     -0.38109206736816036987D+00, &
      0.41413124714030785656D+00, &
     -0.41413124714030785656D+00, &
      0.41413124714030785656D+00, &
     -0.41413124714030785656D+00, &
      0.14311887945931867083D+00, &
     -0.14311887945931867083D+00, &
      0.14311887945931867083D+00, &
     -0.14311887945931867083D+00, &
      0.06816745760360433393D+00, &
     -0.06816745760360433393D+00, &
      0.06816745760360433393D+00, &
     -0.06816745760360433393D+00, &
      0.22977242422012536527D+00, &
     -0.22977242422012536527D+00, &
      0.22977242422012536527D+00, &
     -0.22977242422012536527D+00, &
      0.78211200270072955831D+00, &
     -0.78211200270072955831D+00, &
      0.78211200270072955831D+00, &
     -0.78211200270072955831D+00, &
      0.47646421490609353055D+00, &
     -0.47646421490609353055D+00, &
      0.47646421490609353055D+00, &
     -0.47646421490609353055D+00, &
      0.66550056762002873789D+00, &
     -0.66550056762002873789D+00, &
      0.66550056762002873789D+00, &
     -0.66550056762002873789D+00, &
      0.15755371600385464914D+00, &
     -0.15755371600385464914D+00, &
      0.15755371600385464914D+00, &
     -0.15755371600385464914D+00, &
      0.74760255490899130137D+00, &
      0.13357335558335847736D+00, &
     -0.74760255490899130137D+00, &
      0.13357335558335847736D+00, &
      0.74760255490899130137D+00, &
     -0.13357335558335847736D+00, &
     -0.74760255490899130137D+00, &
     -0.13357335558335847736D+00, &
      0.64024497237564748087D+00, &
      0.36978998725981265805D+00, &
     -0.64024497237564748087D+00, &
      0.36978998725981265805D+00, &
      0.64024497237564748087D+00, &
     -0.36978998725981265805D+00, &
     -0.64024497237564748087D+00, &
     -0.36978998725981265805D+00, &
      0.27725002720745439699D+00, &
      0.65880359279798383909D+00, &
     -0.27725002720745439699D+00, &
      0.65880359279798383909D+00, &
      0.27725002720745439699D+00, &
     -0.65880359279798383909D+00, &
     -0.27725002720745439699D+00, &
     -0.65880359279798383909D+00, &
      0.88435969563492811130D+00, &
      0.76901445510268306993D+00, &
     -0.88435969563492811130D+00, &
      0.76901445510268306993D+00, &
      0.88435969563492811130D+00, &
     -0.76901445510268306993D+00, &
     -0.88435969563492811130D+00, &
     -0.76901445510268306993D+00, &
      0.95432083454417460100D+00, &
      0.27565292566365701132D+00, &
     -0.95432083454417460100D+00, &
      0.27565292566365701132D+00, &
      0.95432083454417460100D+00, &
     -0.27565292566365701132D+00, &
     -0.95432083454417460100D+00, &
     -0.27565292566365701132D+00, &
      0.19610686376179234380D+00, &
      0.48815934289407264535D+00, &
     -0.19610686376179234380D+00, &
      0.48815934289407264535D+00, &
      0.19610686376179234380D+00, &
     -0.48815934289407264535D+00, &
     -0.19610686376179234380D+00, &
     -0.48815934289407264535D+00, &
      0.52539286886562275303D+00, &
      0.90365402298313357576D+00, &
     -0.52539286886562275303D+00, &
      0.90365402298313357576D+00, &
      0.52539286886562275303D+00, &
     -0.90365402298313357576D+00, &
     -0.52539286886562275303D+00, &
     -0.90365402298313357576D+00, &
      0.78710506457062645591D+00, &
      0.53056669085425134380D+00, &
     -0.78710506457062645591D+00, &
      0.53056669085425134380D+00, &
      0.78710506457062645591D+00, &
     -0.53056669085425134380D+00, &
     -0.78710506457062645591D+00, &
     -0.53056669085425134380D+00, &
      0.87305124922384047537D+00, &
      0.58249614598785404151D+00, &
     -0.87305124922384047537D+00, &
      0.58249614598785404151D+00, &
      0.87305124922384047537D+00, &
     -0.58249614598785404151D+00, &
     -0.87305124922384047537D+00, &
     -0.58249614598785404151D+00, &
      0.49539855010659772372D+00, &
      0.21910184606142177333D+00, &
     -0.49539855010659772372D+00, &
      0.21910184606142177333D+00, &
      0.49539855010659772372D+00, &
     -0.21910184606142177333D+00, &
     -0.49539855010659772372D+00, &
     -0.21910184606142177333D+00, &
      0.86957254282469464979D+00, &
      0.30557223260338794990D+00, &
     -0.86957254282469464979D+00, &
      0.30557223260338794990D+00, &
      0.86957254282469464979D+00, &
     -0.30557223260338794990D+00, &
     -0.86957254282469464979D+00, &
     -0.30557223260338794990D+00, &
      0.82556821539161207024D+00, &
      0.18201665219265533713D+00, &
     -0.82556821539161207024D+00, &
      0.18201665219265533713D+00, &
      0.82556821539161207024D+00, &
     -0.18201665219265533713D+00, &
     -0.82556821539161207024D+00, &
     -0.18201665219265533713D+00, &
      0.22633675831897015485D+00, &
      0.35538499166679399233D+00, &
     -0.22633675831897015485D+00, &
      0.35538499166679399233D+00, &
      0.22633675831897015485D+00, &
     -0.35538499166679399233D+00, &
     -0.22633675831897015485D+00, &
     -0.35538499166679399233D+00, &
      0.55136786560508410648D+00, &
      0.74868759522365924131D+00, &
     -0.55136786560508410648D+00, &
      0.74868759522365924131D+00, &
      0.55136786560508410648D+00, &
     -0.74868759522365924131D+00, &
     -0.55136786560508410648D+00, &
     -0.74868759522365924131D+00, &
      0.64168041273857601148D+00, &
      0.25896549169892302267D+00, &
     -0.64168041273857601148D+00, &
      0.25896549169892302267D+00, &
      0.64168041273857601148D+00, &
     -0.25896549169892302267D+00, &
     -0.64168041273857601148D+00, &
     -0.25896549169892302267D+00, &
      0.48623936010342988512D+00, &
      0.11286859018565884027D+00, &
     -0.48623936010342988512D+00, &
      0.11286859018565884027D+00, &
      0.48623936010342988512D+00, &
     -0.11286859018565884027D+00, &
     -0.48623936010342988512D+00, &
     -0.11286859018565884027D+00, &
      0.42728691837450794022D+00, &
      0.62836164145805561976D+00, &
     -0.42728691837450794022D+00, &
      0.62836164145805561976D+00, &
      0.42728691837450794022D+00, &
     -0.62836164145805561976D+00, &
     -0.42728691837450794022D+00, &
     -0.62836164145805561976D+00, &
      0.97870776397134839897D+00, &
      0.73410552715894472620D+00, &
     -0.97870776397134839897D+00, &
      0.73410552715894472620D+00, &
      0.97870776397134839897D+00, &
     -0.73410552715894472620D+00, &
     -0.97870776397134839897D+00, &
     -0.73410552715894472620D+00, &
      0.70585052586641572336D+00, &
      0.41358916179473048658D+00, &
     -0.70585052586641572336D+00, &
      0.41358916179473048658D+00, &
      0.70585052586641572336D+00, &
     -0.41358916179473048658D+00, &
     -0.70585052586641572336D+00, &
     -0.41358916179473048658D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.95241408581244602072D+00, &
      0.00000000000000000000D+00, &
     -0.95241408581244602072D+00, &
      0.00000000000000000000D+00, &
      0.08275692125030435775D+00, &
      0.00000000000000000000D+00, &
     -0.08275692125030435775D+00, &
      0.00000000000000000000D+00, &
      0.26578545043127727032D+00, &
      0.00000000000000000000D+00, &
     -0.26578545043127727032D+00, &
      0.00000000000000000000D+00, &
      0.81155940028479511827D+00, &
      0.00000000000000000000D+00, &
     -0.81155940028479511827D+00, &
      0.00000000000000000000D+00, &
      0.35232799954016802424D+00, &
      0.00000000000000000000D+00, &
     -0.35232799954016802424D+00, &
      0.00000000000000000000D+00, &
      0.60303771218972423984D+00, &
      0.00000000000000000000D+00, &
     -0.60303771218972423984D+00, &
      0.00000000000000000000D+00, &
      0.78402604202159054125D+00, &
      0.00000000000000000000D+00, &
     -0.78402604202159054125D+00, &
      0.00000000000000000000D+00, &
      0.49346472198646390561D+00, &
      0.00000000000000000000D+00, &
     -0.49346472198646390561D+00, &
      0.00000000000000000000D+00, &
      0.35038751648906635294D+00, &
      0.00000000000000000000D+00, &
     -0.35038751648906635294D+00, &
      0.00000000000000000000D+00, &
      0.00336404312107426752D+00, &
      0.00000000000000000000D+00, &
     -0.00336404312107426752D+00, &
      0.00000000000000000000D+00, &
      0.43074512679067450405D+00, &
      0.00000000000000000000D+00, &
     -0.43074512679067450405D+00, &
      0.00000000000000000000D+00, &
      0.59254801282947688890D+00, &
      0.00000000000000000000D+00, &
     -0.59254801282947688890D+00, &
      0.00000000000000000000D+00, &
      0.22550686036585329552D+00, &
      0.00000000000000000000D+00, &
     -0.22550686036585329552D+00, &
      0.63975863740336180729D+00, &
      0.63975863740336180729D+00, &
     -0.63975863740336180729D+00, &
     -0.63975863740336180729D+00, &
      0.16463474879631256886D+00, &
      0.16463474879631256886D+00, &
     -0.16463474879631256886D+00, &
     -0.16463474879631256886D+00, &
      0.49303563319629611916D+00, &
      0.49303563319629611916D+00, &
     -0.49303563319629611916D+00, &
     -0.49303563319629611916D+00, &
      0.74791280396309745004D+00, &
      0.74791280396309745004D+00, &
     -0.74791280396309745004D+00, &
     -0.74791280396309745004D+00, &
      0.26251317864285966808D+00, &
      0.26251317864285966808D+00, &
     -0.26251317864285966808D+00, &
     -0.26251317864285966808D+00, &
      0.17484215589868074003D+00, &
      0.17484215589868074003D+00, &
     -0.17484215589868074003D+00, &
     -0.17484215589868074003D+00, &
      0.27602468785670569718D+00, &
      0.27602468785670569718D+00, &
     -0.27602468785670569718D+00, &
     -0.27602468785670569718D+00, &
      0.43512552408987620334D+00, &
      0.43512552408987620334D+00, &
     -0.43512552408987620334D+00, &
     -0.43512552408987620334D+00, &
      0.95307694157738065410D+00, &
      0.95307694157738065410D+00, &
     -0.95307694157738065410D+00, &
     -0.95307694157738065410D+00, &
      0.48425597905194911474D+00, &
      0.48425597905194911474D+00, &
     -0.48425597905194911474D+00, &
     -0.48425597905194911474D+00, &
      0.29381061843327571648D+00, &
      0.29381061843327571648D+00, &
     -0.29381061843327571648D+00, &
     -0.29381061843327571648D+00, &
      0.16057200503692023452D+00, &
      0.16057200503692023452D+00, &
     -0.16057200503692023452D+00, &
     -0.16057200503692023452D+00, &
      0.81014529214207364749D+00, &
      0.81014529214207364749D+00, &
     -0.81014529214207364749D+00, &
     -0.81014529214207364749D+00, &
      0.18009325754915589402D+00, &
      0.18009325754915589402D+00, &
     -0.18009325754915589402D+00, &
     -0.18009325754915589402D+00, &
      0.50136799407459931022D+00, &
      0.50136799407459931022D+00, &
     -0.50136799407459931022D+00, &
     -0.50136799407459931022D+00, &
      0.92284407855436589863D+00, &
      0.92284407855436589863D+00, &
     -0.92284407855436589863D+00, &
     -0.92284407855436589863D+00, &
      0.32077249570607596629D+00, &
      0.32077249570607596629D+00, &
     -0.32077249570607596629D+00, &
     -0.32077249570607596629D+00, &
      0.02145401504974505866D+00, &
      0.02145401504974505866D+00, &
     -0.02145401504974505866D+00, &
     -0.02145401504974505866D+00, &
      0.65293078741937671250D+00, &
      0.65293078741937671250D+00, &
     -0.65293078741937671250D+00, &
     -0.65293078741937671250D+00, &
      0.38109206736816036987D+00, &
      0.38109206736816036987D+00, &
     -0.38109206736816036987D+00, &
     -0.38109206736816036987D+00, &
      0.41413124714030785656D+00, &
      0.41413124714030785656D+00, &
     -0.41413124714030785656D+00, &
     -0.41413124714030785656D+00, &
      0.14311887945931867083D+00, &
      0.14311887945931867083D+00, &
     -0.14311887945931867083D+00, &
     -0.14311887945931867083D+00, &
      0.06816745760360433393D+00, &
      0.06816745760360433393D+00, &
     -0.06816745760360433393D+00, &
     -0.06816745760360433393D+00, &
      0.22977242422012536527D+00, &
      0.22977242422012536527D+00, &
     -0.22977242422012536527D+00, &
     -0.22977242422012536527D+00, &
      0.78211200270072955831D+00, &
      0.78211200270072955831D+00, &
     -0.78211200270072955831D+00, &
     -0.78211200270072955831D+00, &
      0.47646421490609353055D+00, &
      0.47646421490609353055D+00, &
     -0.47646421490609353055D+00, &
     -0.47646421490609353055D+00, &
      0.66550056762002873789D+00, &
      0.66550056762002873789D+00, &
     -0.66550056762002873789D+00, &
     -0.66550056762002873789D+00, &
      0.15755371600385464914D+00, &
      0.15755371600385464914D+00, &
     -0.15755371600385464914D+00, &
     -0.15755371600385464914D+00, &
      0.13357335558335847736D+00, &
      0.74760255490899130137D+00, &
      0.13357335558335847736D+00, &
     -0.74760255490899130137D+00, &
     -0.13357335558335847736D+00, &
      0.74760255490899130137D+00, &
     -0.13357335558335847736D+00, &
     -0.74760255490899130137D+00, &
      0.36978998725981265805D+00, &
      0.64024497237564748087D+00, &
      0.36978998725981265805D+00, &
     -0.64024497237564748087D+00, &
     -0.36978998725981265805D+00, &
      0.64024497237564748087D+00, &
     -0.36978998725981265805D+00, &
     -0.64024497237564748087D+00, &
      0.65880359279798383909D+00, &
      0.27725002720745439699D+00, &
      0.65880359279798383909D+00, &
     -0.27725002720745439699D+00, &
     -0.65880359279798383909D+00, &
      0.27725002720745439699D+00, &
     -0.65880359279798383909D+00, &
     -0.27725002720745439699D+00, &
      0.76901445510268306993D+00, &
      0.88435969563492811130D+00, &
      0.76901445510268306993D+00, &
     -0.88435969563492811130D+00, &
     -0.76901445510268306993D+00, &
      0.88435969563492811130D+00, &
     -0.76901445510268306993D+00, &
     -0.88435969563492811130D+00, &
      0.27565292566365701132D+00, &
      0.95432083454417460100D+00, &
      0.27565292566365701132D+00, &
     -0.95432083454417460100D+00, &
     -0.27565292566365701132D+00, &
      0.95432083454417460100D+00, &
     -0.27565292566365701132D+00, &
     -0.95432083454417460100D+00, &
      0.48815934289407264535D+00, &
      0.19610686376179234380D+00, &
      0.48815934289407264535D+00, &
     -0.19610686376179234380D+00, &
     -0.48815934289407264535D+00, &
      0.19610686376179234380D+00, &
     -0.48815934289407264535D+00, &
     -0.19610686376179234380D+00, &
      0.90365402298313357576D+00, &
      0.52539286886562275303D+00, &
      0.90365402298313357576D+00, &
     -0.52539286886562275303D+00, &
     -0.90365402298313357576D+00, &
      0.52539286886562275303D+00, &
     -0.90365402298313357576D+00, &
     -0.52539286886562275303D+00, &
      0.53056669085425134380D+00, &
      0.78710506457062645591D+00, &
      0.53056669085425134380D+00, &
     -0.78710506457062645591D+00, &
     -0.53056669085425134380D+00, &
      0.78710506457062645591D+00, &
     -0.53056669085425134380D+00, &
     -0.78710506457062645591D+00, &
      0.58249614598785404151D+00, &
      0.87305124922384047537D+00, &
      0.58249614598785404151D+00, &
     -0.87305124922384047537D+00, &
     -0.58249614598785404151D+00, &
      0.87305124922384047537D+00, &
     -0.58249614598785404151D+00, &
     -0.87305124922384047537D+00, &
      0.21910184606142177333D+00, &
      0.49539855010659772372D+00, &
      0.21910184606142177333D+00, &
     -0.49539855010659772372D+00, &
     -0.21910184606142177333D+00, &
      0.49539855010659772372D+00, &
     -0.21910184606142177333D+00, &
     -0.49539855010659772372D+00, &
      0.30557223260338794990D+00, &
      0.86957254282469464979D+00, &
      0.30557223260338794990D+00, &
     -0.86957254282469464979D+00, &
     -0.30557223260338794990D+00, &
      0.86957254282469464979D+00, &
     -0.30557223260338794990D+00, &
     -0.86957254282469464979D+00, &
      0.18201665219265533713D+00, &
      0.82556821539161207024D+00, &
      0.18201665219265533713D+00, &
     -0.82556821539161207024D+00, &
     -0.18201665219265533713D+00, &
      0.82556821539161207024D+00, &
     -0.18201665219265533713D+00, &
     -0.82556821539161207024D+00, &
      0.35538499166679399233D+00, &
      0.22633675831897015485D+00, &
      0.35538499166679399233D+00, &
     -0.22633675831897015485D+00, &
     -0.35538499166679399233D+00, &
      0.22633675831897015485D+00, &
     -0.35538499166679399233D+00, &
     -0.22633675831897015485D+00, &
      0.74868759522365924131D+00, &
      0.55136786560508410648D+00, &
      0.74868759522365924131D+00, &
     -0.55136786560508410648D+00, &
     -0.74868759522365924131D+00, &
      0.55136786560508410648D+00, &
     -0.74868759522365924131D+00, &
     -0.55136786560508410648D+00, &
      0.25896549169892302267D+00, &
      0.64168041273857601148D+00, &
      0.25896549169892302267D+00, &
     -0.64168041273857601148D+00, &
     -0.25896549169892302267D+00, &
      0.64168041273857601148D+00, &
     -0.25896549169892302267D+00, &
     -0.64168041273857601148D+00, &
      0.11286859018565884027D+00, &
      0.48623936010342988512D+00, &
      0.11286859018565884027D+00, &
     -0.48623936010342988512D+00, &
     -0.11286859018565884027D+00, &
      0.48623936010342988512D+00, &
     -0.11286859018565884027D+00, &
     -0.48623936010342988512D+00, &
      0.62836164145805561976D+00, &
      0.42728691837450794022D+00, &
      0.62836164145805561976D+00, &
     -0.42728691837450794022D+00, &
     -0.62836164145805561976D+00, &
      0.42728691837450794022D+00, &
     -0.62836164145805561976D+00, &
     -0.42728691837450794022D+00, &
      0.73410552715894472620D+00, &
      0.97870776397134839897D+00, &
      0.73410552715894472620D+00, &
     -0.97870776397134839897D+00, &
     -0.73410552715894472620D+00, &
      0.97870776397134839897D+00, &
     -0.73410552715894472620D+00, &
     -0.97870776397134839897D+00, &
      0.41358916179473048658D+00, &
      0.70585052586641572336D+00, &
      0.41358916179473048658D+00, &
     -0.70585052586641572336D+00, &
     -0.41358916179473048658D+00, &
      0.70585052586641572336D+00, &
     -0.41358916179473048658D+00, &
     -0.70585052586641572336D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
      0.97220323387036922114D+00, &
      0.63870352926681406291D+00, &
      0.68020043208055336326D+00, &
      0.04194980824924493534D+00, &
      0.04194980824924493534D+00, &
      0.04194980824924493534D+00, &
      0.04194980824924493534D+00, &
      0.90141183204156405395D+00, &
      0.90141183204156405395D+00, &
      0.90141183204156405395D+00, &
      0.90141183204156405395D+00, &
      0.68432529322768964608D+00, &
      0.68432529322768964608D+00, &
      0.68432529322768964608D+00, &
      0.68432529322768964608D+00, &
      0.00000008245324867689D+00, &
      0.00000008245324867689D+00, &
      0.00000008245324867689D+00, &
      0.00000008245324867689D+00, &
      0.55417390873649197136D+00, &
      0.55417390873649197136D+00, &
      0.55417390873649197136D+00, &
      0.55417390873649197136D+00, &
      0.36258805988768089135D+00, &
      0.36258805988768089135D+00, &
      0.36258805988768089135D+00, &
      0.36258805988768089135D+00, &
      0.13875708503462674814D+00, &
      0.13875708503462674814D+00, &
      0.13875708503462674814D+00, &
      0.13875708503462674814D+00, &
      0.24260167486456202246D+00, &
      0.24260167486456202246D+00, &
      0.24260167486456202246D+00, &
      0.24260167486456202246D+00, &
      0.40130775806887936108D+00, &
      0.40130775806887936108D+00, &
      0.40130775806887936108D+00, &
      0.40130775806887936108D+00, &
      0.38152973269372325582D+00, &
      0.38152973269372325582D+00, &
      0.38152973269372325582D+00, &
      0.38152973269372325582D+00, &
      0.11901456915448617446D+00, &
      0.11901456915448617446D+00, &
      0.11901456915448617446D+00, &
      0.11901456915448617446D+00, &
      0.04928012549325393871D+00, &
      0.04928012549325393871D+00, &
      0.04928012549325393871D+00, &
      0.04928012549325393871D+00, &
      0.76206391519358240849D+00, &
      0.76206391519358240849D+00, &
      0.76206391519358240849D+00, &
      0.76206391519358240849D+00, &
      0.31153469497121605292D+00, &
      0.31153469497121605292D+00, &
      0.31153469497121605292D+00, &
      0.31153469497121605292D+00, &
      0.25650991726217964306D+00, &
      0.25650991726217964306D+00, &
      0.25650991726217964306D+00, &
      0.25650991726217964306D+00, &
      0.07449733603189459541D+00, &
      0.07449733603189459541D+00, &
      0.07449733603189459541D+00, &
      0.07449733603189459541D+00, &
      0.00130788548084566922D+00, &
      0.00130788548084566922D+00, &
      0.00130788548084566922D+00, &
      0.00130788548084566922D+00, &
      0.52858580218461959088D+00, &
      0.52858580218461959088D+00, &
      0.52858580218461959088D+00, &
      0.52858580218461959088D+00, &
      0.00955221072473505343D+00, &
      0.00955221072473505343D+00, &
      0.00955221072473505343D+00, &
      0.00955221072473505343D+00, &
      0.34142989905969745035D+00, &
      0.34142989905969745035D+00, &
      0.34142989905969745035D+00, &
      0.34142989905969745035D+00, &
      0.00027022925967801637D+00, &
      0.00027022925967801637D+00, &
      0.00027022925967801637D+00, &
      0.00027022925967801637D+00, &
      0.04692242455320395911D+00, &
      0.04692242455320395911D+00, &
      0.04692242455320395911D+00, &
      0.04692242455320395911D+00, &
      0.34449713283505384309D+00, &
      0.34449713283505384309D+00, &
      0.34449713283505384309D+00, &
      0.34449713283505384309D+00, &
      0.66811923014208407512D+00, &
      0.66811923014208407512D+00, &
      0.66811923014208407512D+00, &
      0.66811923014208407512D+00, &
      0.67391287157124502016D+00, &
      0.67391287157124502016D+00, &
      0.67391287157124502016D+00, &
      0.67391287157124502016D+00, &
      0.03437584591482227558D+00, &
      0.03437584591482227558D+00, &
      0.03437584591482227558D+00, &
      0.03437584591482227558D+00, &
      0.15542872578297051156D+00, &
      0.15542872578297051156D+00, &
      0.15542872578297051156D+00, &
      0.15542872578297051156D+00, &
      0.26971828733855934823D+00, &
      0.26971828733855934823D+00, &
      0.26971828733855934823D+00, &
      0.26971828733855934823D+00, &
      0.01181428825288092545D+00, &
      0.01181428825288092545D+00, &
      0.01181428825288092545D+00, &
      0.01181428825288092545D+00, &
      0.06366317285532444026D+00, &
      0.06366317285532444026D+00, &
      0.06366317285532444026D+00, &
      0.06366317285532444026D+00, &
      0.06762724457801919109D+00, &
      0.06762724457801919109D+00, &
      0.06762724457801919109D+00, &
      0.06762724457801919109D+00, &
      0.20149664605408518225D+00, &
      0.20149664605408518225D+00, &
      0.20149664605408518225D+00, &
      0.20149664605408518225D+00, &
      0.48586744666719888786D+00, &
      0.48586744666719888786D+00, &
      0.48586744666719888786D+00, &
      0.48586744666719888786D+00, &
      0.19071052464250978775D+00, &
      0.19071052464250978775D+00, &
      0.19071052464250978775D+00, &
      0.19071052464250978775D+00, &
      0.51037375853206956577D+00, &
      0.51037375853206956577D+00, &
      0.51037375853206956577D+00, &
      0.51037375853206956577D+00, &
      0.80351245513246505325D+00, &
      0.80351245513246505325D+00, &
      0.80351245513246505325D+00, &
      0.80351245513246505325D+00, &
      0.04121355156423638089D+00, &
      0.04121355156423638089D+00, &
      0.04121355156423638089D+00, &
      0.04121355156423638089D+00, &
      0.17760016176122497833D+00, &
      0.17760016176122497833D+00, &
      0.17760016176122497833D+00, &
      0.17760016176122497833D+00, &
      0.48439012439328182902D+00, &
      0.48439012439328182902D+00, &
      0.48439012439328182902D+00, &
      0.48439012439328182902D+00, &
      0.04796979480806326523D+00, &
      0.04796979480806326523D+00, &
      0.04796979480806326523D+00, &
      0.04796979480806326523D+00, &
      0.80108123458344093759D+00, &
      0.80108123458344093759D+00, &
      0.80108123458344093759D+00, &
      0.80108123458344093759D+00, &
      0.24101422703502115019D+00, &
      0.24101422703502115019D+00, &
      0.24101422703502115019D+00, &
      0.24101422703502115019D+00, &
      0.24101422703502115019D+00, &
      0.24101422703502115019D+00, &
      0.24101422703502115019D+00, &
      0.24101422703502115019D+00, &
      0.01069891762813818432D+00, &
      0.01069891762813818432D+00, &
      0.01069891762813818432D+00, &
      0.01069891762813818432D+00, &
      0.01069891762813818432D+00, &
      0.01069891762813818432D+00, &
      0.01069891762813818432D+00, &
      0.01069891762813818432D+00, &
      0.23988631383389671936D+00, &
      0.23988631383389671936D+00, &
      0.23988631383389671936D+00, &
      0.23988631383389671936D+00, &
      0.23988631383389671936D+00, &
      0.23988631383389671936D+00, &
      0.23988631383389671936D+00, &
      0.23988631383389671936D+00, &
      0.08514410320936129095D+00, &
      0.08514410320936129095D+00, &
      0.08514410320936129095D+00, &
      0.08514410320936129095D+00, &
      0.08514410320936129095D+00, &
      0.08514410320936129095D+00, &
      0.08514410320936129095D+00, &
      0.08514410320936129095D+00, &
      0.00847110534891616467D+00, &
      0.00847110534891616467D+00, &
      0.00847110534891616467D+00, &
      0.00847110534891616467D+00, &
      0.00847110534891616467D+00, &
      0.00847110534891616467D+00, &
      0.00847110534891616467D+00, &
      0.00847110534891616467D+00, &
      0.49562245065536153499D+00, &
      0.49562245065536153499D+00, &
      0.49562245065536153499D+00, &
      0.49562245065536153499D+00, &
      0.49562245065536153499D+00, &
      0.49562245065536153499D+00, &
      0.49562245065536153499D+00, &
      0.49562245065536153499D+00, &
      0.04642983560786508729D+00, &
      0.04642983560786508729D+00, &
      0.04642983560786508729D+00, &
      0.04642983560786508729D+00, &
      0.04642983560786508729D+00, &
      0.04642983560786508729D+00, &
      0.04642983560786508729D+00, &
      0.04642983560786508729D+00, &
      0.19656151852896430743D+00, &
      0.19656151852896430743D+00, &
      0.19656151852896430743D+00, &
      0.19656151852896430743D+00, &
      0.19656151852896430743D+00, &
      0.19656151852896430743D+00, &
      0.19656151852896430743D+00, &
      0.19656151852896430743D+00, &
      0.00687572309208568308D+00, &
      0.00687572309208568308D+00, &
      0.00687572309208568308D+00, &
      0.00687572309208568308D+00, &
      0.00687572309208568308D+00, &
      0.00687572309208568308D+00, &
      0.00687572309208568308D+00, &
      0.00687572309208568308D+00, &
      0.38670513013019941484D+00, &
      0.38670513013019941484D+00, &
      0.38670513013019941484D+00, &
      0.38670513013019941484D+00, &
      0.38670513013019941484D+00, &
      0.38670513013019941484D+00, &
      0.38670513013019941484D+00, &
      0.38670513013019941484D+00, &
      0.10896106798490250156D+00, &
      0.10896106798490250156D+00, &
      0.10896106798490250156D+00, &
      0.10896106798490250156D+00, &
      0.10896106798490250156D+00, &
      0.10896106798490250156D+00, &
      0.10896106798490250156D+00, &
      0.10896106798490250156D+00, &
      0.04416884971526349735D+00, &
      0.04416884971526349735D+00, &
      0.04416884971526349735D+00, &
      0.04416884971526349735D+00, &
      0.04416884971526349735D+00, &
      0.04416884971526349735D+00, &
      0.04416884971526349735D+00, &
      0.04416884971526349735D+00, &
      0.62489254847899400325D+00, &
      0.62489254847899400325D+00, &
      0.62489254847899400325D+00, &
      0.62489254847899400325D+00, &
      0.62489254847899400325D+00, &
      0.62489254847899400325D+00, &
      0.62489254847899400325D+00, &
      0.62489254847899400325D+00, &
      0.11727329599993376041D+00, &
      0.11727329599993376041D+00, &
      0.11727329599993376041D+00, &
      0.11727329599993376041D+00, &
      0.11727329599993376041D+00, &
      0.11727329599993376041D+00, &
      0.11727329599993376041D+00, &
      0.11727329599993376041D+00, &
      0.11958436259454885420D+00, &
      0.11958436259454885420D+00, &
      0.11958436259454885420D+00, &
      0.11958436259454885420D+00, &
      0.11958436259454885420D+00, &
      0.11958436259454885420D+00, &
      0.11958436259454885420D+00, &
      0.11958436259454885420D+00, &
      0.01562091696416078104D+00, &
      0.01562091696416078104D+00, &
      0.01562091696416078104D+00, &
      0.01562091696416078104D+00, &
      0.01562091696416078104D+00, &
      0.01562091696416078104D+00, &
      0.01562091696416078104D+00, &
      0.01562091696416078104D+00, &
      0.35132620557004257122D+00, &
      0.35132620557004257122D+00, &
      0.35132620557004257122D+00, &
      0.35132620557004257122D+00, &
      0.35132620557004257122D+00, &
      0.35132620557004257122D+00, &
      0.35132620557004257122D+00, &
      0.35132620557004257122D+00, &
      0.01323405399603627339D+00, &
      0.01323405399603627339D+00, &
      0.01323405399603627339D+00, &
      0.01323405399603627339D+00, &
      0.01323405399603627339D+00, &
      0.01323405399603627339D+00, &
      0.01323405399603627339D+00, &
      0.01323405399603627339D+00, &
      0.02986233869137579211D+00, &
      0.02986233869137579211D+00, &
      0.02986233869137579211D+00, &
      0.02986233869137579211D+00, &
      0.02986233869137579211D+00, &
      0.02986233869137579211D+00, &
      0.02986233869137579211D+00, &
      0.02986233869137579211D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
      0.00012216821078901854D+00, &
      0.00563943629581533951D+00, &
      0.00187116838043227779D+00, &
      0.00087900554061397717D+00, &
      0.00087900554061397717D+00, &
      0.00087900554061397717D+00, &
      0.00087900554061397717D+00, &
      0.00063973688323901037D+00, &
      0.00063973688323901037D+00, &
      0.00063973688323901037D+00, &
      0.00063973688323901037D+00, &
      0.00191942477958428285D+00, &
      0.00191942477958428285D+00, &
      0.00191942477958428285D+00, &
      0.00191942477958428285D+00, &
      0.00108856800183791767D+00, &
      0.00108856800183791767D+00, &
      0.00108856800183791767D+00, &
      0.00108856800183791767D+00, &
      0.00579178348482820276D+00, &
      0.00579178348482820276D+00, &
      0.00579178348482820276D+00, &
      0.00579178348482820276D+00, &
      0.00333567983247391557D+00, &
      0.00333567983247391557D+00, &
      0.00333567983247391557D+00, &
      0.00333567983247391557D+00, &
      0.00483993401414843216D+00, &
      0.00483993401414843216D+00, &
      0.00483993401414843216D+00, &
      0.00483993401414843216D+00, &
      0.01171223275540323253D+00, &
      0.01171223275540323253D+00, &
      0.01171223275540323253D+00, &
      0.01171223275540323253D+00, &
      0.00781803283741232084D+00, &
      0.00781803283741232084D+00, &
      0.00781803283741232084D+00, &
      0.00781803283741232084D+00, &
      0.00316389891636655388D+00, &
      0.00316389891636655388D+00, &
      0.00316389891636655388D+00, &
      0.00316389891636655388D+00, &
      0.00813627882273969955D+00, &
      0.00813627882273969955D+00, &
      0.00813627882273969955D+00, &
      0.00813627882273969955D+00, &
      0.00471982615752600254D+00, &
      0.00471982615752600254D+00, &
      0.00471982615752600254D+00, &
      0.00471982615752600254D+00, &
      0.00096237134938021330D+00, &
      0.00096237134938021330D+00, &
      0.00096237134938021330D+00, &
      0.00096237134938021330D+00, &
      0.00162378774013503748D+00, &
      0.00162378774013503748D+00, &
      0.00162378774013503748D+00, &
      0.00162378774013503748D+00, &
      0.00871345927599746599D+00, &
      0.00871345927599746599D+00, &
      0.00871345927599746599D+00, &
      0.00871345927599746599D+00, &
      0.00412008618466326321D+00, &
      0.00412008618466326321D+00, &
      0.00412008618466326321D+00, &
      0.00412008618466326321D+00, &
      0.00052206472289359268D+00, &
      0.00052206472289359268D+00, &
      0.00052206472289359268D+00, &
      0.00052206472289359268D+00, &
      0.00389716456343351799D+00, &
      0.00389716456343351799D+00, &
      0.00389716456343351799D+00, &
      0.00389716456343351799D+00, &
      0.00222992187534444422D+00, &
      0.00222992187534444422D+00, &
      0.00222992187534444422D+00, &
      0.00222992187534444422D+00, &
      0.00941806762258396400D+00, &
      0.00941806762258396400D+00, &
      0.00941806762258396400D+00, &
      0.00941806762258396400D+00, &
      0.00066158322303103113D+00, &
      0.00066158322303103113D+00, &
      0.00066158322303103113D+00, &
      0.00066158322303103113D+00, &
      0.00010739422189895962D+00, &
      0.00010739422189895962D+00, &
      0.00010739422189895962D+00, &
      0.00010739422189895962D+00, &
      0.00303147852588211194D+00, &
      0.00303147852588211194D+00, &
      0.00303147852588211194D+00, &
      0.00303147852588211194D+00, &
      0.00081932641066923071D+00, &
      0.00081932641066923071D+00, &
      0.00081932641066923071D+00, &
      0.00081932641066923071D+00, &
      0.00465121449070608067D+00, &
      0.00465121449070608067D+00, &
      0.00465121449070608067D+00, &
      0.00465121449070608067D+00, &
      0.00167606716811515987D+00, &
      0.00167606716811515987D+00, &
      0.00167606716811515987D+00, &
      0.00167606716811515987D+00, &
      0.00825768957465656227D+00, &
      0.00825768957465656227D+00, &
      0.00825768957465656227D+00, &
      0.00825768957465656227D+00, &
      0.00397545859544364086D+00, &
      0.00397545859544364086D+00, &
      0.00397545859544364086D+00, &
      0.00397545859544364086D+00, &
      0.00051722520547986095D+00, &
      0.00051722520547986095D+00, &
      0.00051722520547986095D+00, &
      0.00051722520547986095D+00, &
      0.00624136348038179748D+00, &
      0.00624136348038179748D+00, &
      0.00624136348038179748D+00, &
      0.00624136348038179748D+00, &
      0.00249954122416236671D+00, &
      0.00249954122416236671D+00, &
      0.00249954122416236671D+00, &
      0.00249954122416236671D+00, &
      0.00344134582405346756D+00, &
      0.00344134582405346756D+00, &
      0.00344134582405346756D+00, &
      0.00344134582405346756D+00, &
      0.00373085764189791565D+00, &
      0.00373085764189791565D+00, &
      0.00373085764189791565D+00, &
      0.00373085764189791565D+00, &
      0.00985379501780648571D+00, &
      0.00985379501780648571D+00, &
      0.00985379501780648571D+00, &
      0.00985379501780648571D+00, &
      0.00736757952028449341D+00, &
      0.00736757952028449341D+00, &
      0.00736757952028449341D+00, &
      0.00736757952028449341D+00, &
      0.00178048730560625978D+00, &
      0.00178048730560625978D+00, &
      0.00178048730560625978D+00, &
      0.00178048730560625978D+00, &
      0.00242313168519621295D+00, &
      0.00242313168519621295D+00, &
      0.00242313168519621295D+00, &
      0.00242313168519621295D+00, &
      0.00105734008948304265D+00, &
      0.00105734008948304265D+00, &
      0.00105734008948304265D+00, &
      0.00105734008948304265D+00, &
      0.00145036343491798387D+00, &
      0.00145036343491798387D+00, &
      0.00145036343491798387D+00, &
      0.00145036343491798387D+00, &
      0.00271057993574526021D+00, &
      0.00271057993574526021D+00, &
      0.00271057993574526021D+00, &
      0.00271057993574526021D+00, &
      0.00112014566884451267D+00, &
      0.00112014566884451267D+00, &
      0.00112014566884451267D+00, &
      0.00112014566884451267D+00, &
      0.00120917104359697324D+00, &
      0.00120917104359697324D+00, &
      0.00120917104359697324D+00, &
      0.00120917104359697324D+00, &
      0.00120917104359697324D+00, &
      0.00120917104359697324D+00, &
      0.00120917104359697324D+00, &
      0.00120917104359697324D+00, &
      0.00162471395073617649D+00, &
      0.00162471395073617649D+00, &
      0.00162471395073617649D+00, &
      0.00162471395073617649D+00, &
      0.00162471395073617649D+00, &
      0.00162471395073617649D+00, &
      0.00162471395073617649D+00, &
      0.00162471395073617649D+00, &
      0.00620858966914663961D+00, &
      0.00620858966914663961D+00, &
      0.00620858966914663961D+00, &
      0.00620858966914663961D+00, &
      0.00620858966914663961D+00, &
      0.00620858966914663961D+00, &
      0.00620858966914663961D+00, &
      0.00620858966914663961D+00, &
      0.00135724577802302663D+00, &
      0.00135724577802302663D+00, &
      0.00135724577802302663D+00, &
      0.00135724577802302663D+00, &
      0.00135724577802302663D+00, &
      0.00135724577802302663D+00, &
      0.00135724577802302663D+00, &
      0.00135724577802302663D+00, &
      0.00077688230182242947D+00, &
      0.00077688230182242947D+00, &
      0.00077688230182242947D+00, &
      0.00077688230182242947D+00, &
      0.00077688230182242947D+00, &
      0.00077688230182242947D+00, &
      0.00077688230182242947D+00, &
      0.00077688230182242947D+00, &
      0.00190538850038805656D+00, &
      0.00190538850038805656D+00, &
      0.00190538850038805656D+00, &
      0.00190538850038805656D+00, &
      0.00190538850038805656D+00, &
      0.00190538850038805656D+00, &
      0.00190538850038805656D+00, &
      0.00190538850038805656D+00, &
      0.00183635585567423158D+00, &
      0.00183635585567423158D+00, &
      0.00183635585567423158D+00, &
      0.00183635585567423158D+00, &
      0.00183635585567423158D+00, &
      0.00183635585567423158D+00, &
      0.00183635585567423158D+00, &
      0.00183635585567423158D+00, &
      0.00173312265591197071D+00, &
      0.00173312265591197071D+00, &
      0.00173312265591197071D+00, &
      0.00173312265591197071D+00, &
      0.00173312265591197071D+00, &
      0.00173312265591197071D+00, &
      0.00173312265591197071D+00, &
      0.00173312265591197071D+00, &
      0.00087558681440168869D+00, &
      0.00087558681440168869D+00, &
      0.00087558681440168869D+00, &
      0.00087558681440168869D+00, &
      0.00087558681440168869D+00, &
      0.00087558681440168869D+00, &
      0.00087558681440168869D+00, &
      0.00087558681440168869D+00, &
      0.00529919935281060114D+00, &
      0.00529919935281060114D+00, &
      0.00529919935281060114D+00, &
      0.00529919935281060114D+00, &
      0.00529919935281060114D+00, &
      0.00529919935281060114D+00, &
      0.00529919935281060114D+00, &
      0.00529919935281060114D+00, &
      0.00190923140339424222D+00, &
      0.00190923140339424222D+00, &
      0.00190923140339424222D+00, &
      0.00190923140339424222D+00, &
      0.00190923140339424222D+00, &
      0.00190923140339424222D+00, &
      0.00190923140339424222D+00, &
      0.00190923140339424222D+00, &
      0.00333270329648946212D+00, &
      0.00333270329648946212D+00, &
      0.00333270329648946212D+00, &
      0.00333270329648946212D+00, &
      0.00333270329648946212D+00, &
      0.00333270329648946212D+00, &
      0.00333270329648946212D+00, &
      0.00333270329648946212D+00, &
      0.00137258273942558089D+00, &
      0.00137258273942558089D+00, &
      0.00137258273942558089D+00, &
      0.00137258273942558089D+00, &
      0.00137258273942558089D+00, &
      0.00137258273942558089D+00, &
      0.00137258273942558089D+00, &
      0.00137258273942558089D+00, &
      0.00440562895172808174D+00, &
      0.00440562895172808174D+00, &
      0.00440562895172808174D+00, &
      0.00440562895172808174D+00, &
      0.00440562895172808174D+00, &
      0.00440562895172808174D+00, &
      0.00440562895172808174D+00, &
      0.00440562895172808174D+00, &
      0.00690265770456989416D+00, &
      0.00690265770456989416D+00, &
      0.00690265770456989416D+00, &
      0.00690265770456989416D+00, &
      0.00690265770456989416D+00, &
      0.00690265770456989416D+00, &
      0.00690265770456989416D+00, &
      0.00690265770456989416D+00, &
      0.00202734469527626228D+00, &
      0.00202734469527626228D+00, &
      0.00202734469527626228D+00, &
      0.00202734469527626228D+00, &
      0.00202734469527626228D+00, &
      0.00202734469527626228D+00, &
      0.00202734469527626228D+00, &
      0.00202734469527626228D+00, &
      0.00204643595334019031D+00, &
      0.00204643595334019031D+00, &
      0.00204643595334019031D+00, &
      0.00204643595334019031D+00, &
      0.00204643595334019031D+00, &
      0.00204643595334019031D+00, &
      0.00204643595334019031D+00, &
      0.00204643595334019031D+00, &
      0.00034532641662137929D+00, &
      0.00034532641662137929D+00, &
      0.00034532641662137929D+00, &
      0.00034532641662137929D+00, &
      0.00034532641662137929D+00, &
      0.00034532641662137929D+00, &
      0.00034532641662137929D+00, &
      0.00034532641662137929D+00, &
      0.00242508950332979241D+00, &
      0.00242508950332979241D+00, &
      0.00242508950332979241D+00, &
      0.00242508950332979241D+00, &
      0.00242508950332979241D+00, &
      0.00242508950332979241D+00, &
      0.00242508950332979241D+00, &
      0.00242508950332979241D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule18 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule18() returns the pyramid quadrature rule of precision 18.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 April 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 121, Number 11, pages 2418-2436, 15 June 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real ( kind = rk ) x[n], y[n], z[n]: the coordinates of
!    quadrature points.
!
!    real ( kind = rk ) w[n]: the quadrature weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 357

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
      0.00000000000000000000D+00, &
      0.73869275831950009792D+00, &
      0.00000000000000000000D+00, &
     -0.73869275831950009792D+00, &
      0.00000000000000000000D+00, &
      0.76806657302952874300D+00, &
      0.00000000000000000000D+00, &
     -0.76806657302952874300D+00, &
      0.00000000000000000000D+00, &
      0.31673944809898013641D+00, &
      0.00000000000000000000D+00, &
     -0.31673944809898013641D+00, &
      0.00000000000000000000D+00, &
      0.61877572218072185439D+00, &
      0.00000000000000000000D+00, &
     -0.61877572218072185439D+00, &
      0.00000000000000000000D+00, &
      0.37235294402955443349D+00, &
      0.00000000000000000000D+00, &
     -0.37235294402955443349D+00, &
      0.00000000000000000000D+00, &
      0.50387193570782862206D+00, &
      0.00000000000000000000D+00, &
     -0.50387193570782862206D+00, &
      0.00000000000000000000D+00, &
      0.72207836507269651527D+00, &
      0.00000000000000000000D+00, &
     -0.72207836507269651527D+00, &
      0.00000000000000000000D+00, &
      0.39403113359582681019D+00, &
      0.00000000000000000000D+00, &
     -0.39403113359582681019D+00, &
      0.00000000000000000000D+00, &
      0.11819479806124476295D+00, &
      0.00000000000000000000D+00, &
     -0.11819479806124476295D+00, &
      0.00000000000000000000D+00, &
      0.51179862044264223808D+00, &
      0.00000000000000000000D+00, &
     -0.51179862044264223808D+00, &
      0.00000000000000000000D+00, &
      0.09344056049843052492D+00, &
      0.00000000000000000000D+00, &
     -0.09344056049843052492D+00, &
      0.00000000000000000000D+00, &
      0.97786305288383967849D+00, &
      0.00000000000000000000D+00, &
     -0.97786305288383967849D+00, &
      0.00000000000000000000D+00, &
      0.26893118028825491184D+00, &
      0.00000000000000000000D+00, &
     -0.26893118028825491184D+00, &
      0.00000000000000000000D+00, &
      0.18588876717193755783D+00, &
      0.00000000000000000000D+00, &
     -0.18588876717193755783D+00, &
      0.00000000000000000000D+00, &
      0.88032062123229137818D+00, &
     -0.88032062123229137818D+00, &
      0.88032062123229137818D+00, &
     -0.88032062123229137818D+00, &
      0.68145212136761557087D+00, &
     -0.68145212136761557087D+00, &
      0.68145212136761557087D+00, &
     -0.68145212136761557087D+00, &
      0.49699036434274163065D+00, &
     -0.49699036434274163065D+00, &
      0.49699036434274163065D+00, &
     -0.49699036434274163065D+00, &
      0.30634018082228819946D+00, &
     -0.30634018082228819946D+00, &
      0.30634018082228819946D+00, &
     -0.30634018082228819946D+00, &
      0.58539966455385050725D+00, &
     -0.58539966455385050725D+00, &
      0.58539966455385050725D+00, &
     -0.58539966455385050725D+00, &
      0.77694395476943478585D+00, &
     -0.77694395476943478585D+00, &
      0.77694395476943478585D+00, &
     -0.77694395476943478585D+00, &
      0.18695368401965414829D+00, &
     -0.18695368401965414829D+00, &
      0.18695368401965414829D+00, &
     -0.18695368401965414829D+00, &
      0.64402844061291575350D+00, &
     -0.64402844061291575350D+00, &
      0.64402844061291575350D+00, &
     -0.64402844061291575350D+00, &
      0.19939684409397287479D+00, &
     -0.19939684409397287479D+00, &
      0.19939684409397287479D+00, &
     -0.19939684409397287479D+00, &
      0.49151662722959021945D+00, &
     -0.49151662722959021945D+00, &
      0.49151662722959021945D+00, &
     -0.49151662722959021945D+00, &
      0.21971883736597888626D+00, &
     -0.21971883736597888626D+00, &
      0.21971883736597888626D+00, &
     -0.21971883736597888626D+00, &
      0.38181860275292400786D+00, &
     -0.38181860275292400786D+00, &
      0.38181860275292400786D+00, &
     -0.38181860275292400786D+00, &
      0.13233574436496933768D+00, &
     -0.13233574436496933768D+00, &
      0.13233574436496933768D+00, &
     -0.13233574436496933768D+00, &
      0.26434906680031583504D+00, &
     -0.26434906680031583504D+00, &
      0.26434906680031583504D+00, &
     -0.26434906680031583504D+00, &
      0.26259016722585765136D+00, &
     -0.26259016722585765136D+00, &
      0.26259016722585765136D+00, &
     -0.26259016722585765136D+00, &
      0.63683049968859806178D+00, &
     -0.63683049968859806178D+00, &
      0.63683049968859806178D+00, &
     -0.63683049968859806178D+00, &
      0.50333230035733345087D+00, &
     -0.50333230035733345087D+00, &
      0.50333230035733345087D+00, &
     -0.50333230035733345087D+00, &
      0.43053066487698737896D+00, &
     -0.43053066487698737896D+00, &
      0.43053066487698737896D+00, &
     -0.43053066487698737896D+00, &
      0.44121177996231647489D+00, &
     -0.44121177996231647489D+00, &
      0.44121177996231647489D+00, &
     -0.44121177996231647489D+00, &
      0.12338526534910826404D+00, &
     -0.12338526534910826404D+00, &
      0.12338526534910826404D+00, &
     -0.12338526534910826404D+00, &
      0.28196466746787440805D+00, &
     -0.28196466746787440805D+00, &
      0.28196466746787440805D+00, &
     -0.28196466746787440805D+00, &
      0.12302065942509179952D+00, &
     -0.12302065942509179952D+00, &
      0.12302065942509179952D+00, &
     -0.12302065942509179952D+00, &
      0.96995340856169265376D+00, &
     -0.96995340856169265376D+00, &
      0.96995340856169265376D+00, &
     -0.96995340856169265376D+00, &
      0.15367597423015844083D+00, &
     -0.15367597423015844083D+00, &
      0.15367597423015844083D+00, &
     -0.15367597423015844083D+00, &
      0.32324916002943732130D+00, &
     -0.32324916002943732130D+00, &
      0.32324916002943732130D+00, &
     -0.32324916002943732130D+00, &
      0.07979061002341172881D+00, &
     -0.07979061002341172881D+00, &
      0.07979061002341172881D+00, &
     -0.07979061002341172881D+00, &
      0.05572751587573096521D+00, &
     -0.05572751587573096521D+00, &
      0.05572751587573096521D+00, &
     -0.05572751587573096521D+00, &
      0.80407180079417317486D+00, &
     -0.80407180079417317486D+00, &
      0.80407180079417317486D+00, &
     -0.80407180079417317486D+00, &
      0.10630794929321137066D+00, &
     -0.10630794929321137066D+00, &
      0.10630794929321137066D+00, &
     -0.10630794929321137066D+00, &
      0.14968116552008028930D+00, &
     -0.14968116552008028930D+00, &
      0.14968116552008028930D+00, &
     -0.14968116552008028930D+00, &
      0.72179802358977040999D+00, &
     -0.72179802358977040999D+00, &
      0.72179802358977040999D+00, &
     -0.72179802358977040999D+00, &
      0.81797170822487763608D+00, &
      0.64320215485058829241D+00, &
     -0.81797170822487763608D+00, &
      0.64320215485058829241D+00, &
      0.81797170822487763608D+00, &
     -0.64320215485058829241D+00, &
     -0.81797170822487763608D+00, &
     -0.64320215485058829241D+00, &
      0.45477772302529712034D+00, &
      0.16636418997995866542D+00, &
     -0.45477772302529712034D+00, &
      0.16636418997995866542D+00, &
      0.45477772302529712034D+00, &
     -0.16636418997995866542D+00, &
     -0.45477772302529712034D+00, &
     -0.16636418997995866542D+00, &
      0.22077509021471022899D+00, &
      0.64919419286301172090D+00, &
     -0.22077509021471022899D+00, &
      0.64919419286301172090D+00, &
      0.22077509021471022899D+00, &
     -0.64919419286301172090D+00, &
     -0.22077509021471022899D+00, &
     -0.64919419286301172090D+00, &
      0.85849144377065222944D+00, &
      0.26442320916257716634D+00, &
     -0.85849144377065222944D+00, &
      0.26442320916257716634D+00, &
      0.85849144377065222944D+00, &
     -0.26442320916257716634D+00, &
     -0.85849144377065222944D+00, &
     -0.26442320916257716634D+00, &
      0.19087247051624386951D+00, &
      0.45655985846318541954D+00, &
     -0.19087247051624386951D+00, &
      0.45655985846318541954D+00, &
      0.19087247051624386951D+00, &
     -0.45655985846318541954D+00, &
     -0.19087247051624386951D+00, &
     -0.45655985846318541954D+00, &
      0.37787816613860053527D+00, &
      0.66406293254753101518D+00, &
     -0.37787816613860053527D+00, &
      0.66406293254753101518D+00, &
      0.37787816613860053527D+00, &
     -0.66406293254753101518D+00, &
     -0.37787816613860053527D+00, &
     -0.66406293254753101518D+00, &
      0.89395556646395457623D+00, &
      0.17283329178347239807D+00, &
     -0.89395556646395457623D+00, &
      0.17283329178347239807D+00, &
      0.89395556646395457623D+00, &
     -0.17283329178347239807D+00, &
     -0.89395556646395457623D+00, &
     -0.17283329178347239807D+00, &
      0.40991415053935187363D+00, &
      0.05003413197772252352D+00, &
     -0.40991415053935187363D+00, &
      0.05003413197772252352D+00, &
      0.40991415053935187363D+00, &
     -0.05003413197772252352D+00, &
     -0.40991415053935187363D+00, &
     -0.05003413197772252352D+00, &
      0.93198108654420652730D+00, &
      0.10662173986061880548D+00, &
     -0.93198108654420652730D+00, &
      0.10662173986061880548D+00, &
      0.93198108654420652730D+00, &
     -0.10662173986061880548D+00, &
     -0.93198108654420652730D+00, &
     -0.10662173986061880548D+00, &
      0.28018220607206512085D+00, &
      0.84993195986801317598D+00, &
     -0.28018220607206512085D+00, &
      0.84993195986801317598D+00, &
      0.28018220607206512085D+00, &
     -0.84993195986801317598D+00, &
     -0.28018220607206512085D+00, &
     -0.84993195986801317598D+00, &
      0.57728165413613941048D+00, &
      0.46201429744731292715D+00, &
     -0.57728165413613941048D+00, &
      0.46201429744731292715D+00, &
      0.57728165413613941048D+00, &
     -0.46201429744731292715D+00, &
     -0.57728165413613941048D+00, &
     -0.46201429744731292715D+00, &
      0.84451580808527704214D+00, &
      0.61756791755060624904D+00, &
     -0.84451580808527704214D+00, &
      0.61756791755060624904D+00, &
      0.84451580808527704214D+00, &
     -0.61756791755060624904D+00, &
     -0.84451580808527704214D+00, &
     -0.61756791755060624904D+00, &
      0.52941572893811861267D+00, &
      0.23348533584682021336D+00, &
     -0.52941572893811861267D+00, &
      0.23348533584682021336D+00, &
      0.52941572893811861267D+00, &
     -0.23348533584682021336D+00, &
     -0.52941572893811861267D+00, &
     -0.23348533584682021336D+00, &
      0.92827521801894674613D+00, &
      0.67941786279817284466D+00, &
     -0.92827521801894674613D+00, &
      0.67941786279817284466D+00, &
      0.92827521801894674613D+00, &
     -0.67941786279817284466D+00, &
     -0.92827521801894674613D+00, &
     -0.67941786279817284466D+00, &
      0.96836516384095050469D+00, &
      0.46144251438466377113D+00, &
     -0.96836516384095050469D+00, &
      0.46144251438466377113D+00, &
      0.96836516384095050469D+00, &
     -0.46144251438466377113D+00, &
     -0.96836516384095050469D+00, &
     -0.46144251438466377113D+00, &
      0.16744571840030086918D+00, &
      0.31430414794159466929D+00, &
     -0.16744571840030086918D+00, &
      0.31430414794159466929D+00, &
      0.16744571840030086918D+00, &
     -0.31430414794159466929D+00, &
     -0.16744571840030086918D+00, &
     -0.31430414794159466929D+00, &
      0.49232718454496093852D+00, &
      0.73988842925403419670D+00, &
     -0.49232718454496093852D+00, &
      0.73988842925403419670D+00, &
      0.49232718454496093852D+00, &
     -0.73988842925403419670D+00, &
     -0.49232718454496093852D+00, &
     -0.73988842925403419670D+00, &
      0.85659974228270108210D+00, &
      0.54061065617200576572D+00, &
     -0.85659974228270108210D+00, &
      0.54061065617200576572D+00, &
      0.85659974228270108210D+00, &
     -0.54061065617200576572D+00, &
     -0.85659974228270108210D+00, &
     -0.54061065617200576572D+00, &
      0.26702277361260906563D+00, &
      0.57531393529153806998D+00, &
     -0.26702277361260906563D+00, &
      0.57531393529153806998D+00, &
      0.26702277361260906563D+00, &
     -0.57531393529153806998D+00, &
     -0.26702277361260906563D+00, &
     -0.57531393529153806998D+00, &
      0.41506457665342233465D+00, &
      0.71731018735152574095D+00, &
     -0.41506457665342233465D+00, &
      0.71731018735152574095D+00, &
      0.41506457665342233465D+00, &
     -0.71731018735152574095D+00, &
     -0.41506457665342233465D+00, &
     -0.71731018735152574095D+00, &
      0.94466239579161759288D+00, &
      0.81619080186595505122D+00, &
     -0.94466239579161759288D+00, &
      0.81619080186595505122D+00, &
      0.94466239579161759288D+00, &
     -0.81619080186595505122D+00, &
     -0.94466239579161759288D+00, &
     -0.81619080186595505122D+00, &
      0.71450522624081691525D+00, &
      0.31179783227055996031D+00, &
     -0.71450522624081691525D+00, &
      0.31179783227055996031D+00, &
      0.71450522624081691525D+00, &
     -0.31179783227055996031D+00, &
     -0.71450522624081691525D+00, &
     -0.31179783227055996031D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.73869275831950009792D+00, &
      0.00000000000000000000D+00, &
     -0.73869275831950009792D+00, &
      0.00000000000000000000D+00, &
      0.76806657302952874300D+00, &
      0.00000000000000000000D+00, &
     -0.76806657302952874300D+00, &
      0.00000000000000000000D+00, &
      0.31673944809898013641D+00, &
      0.00000000000000000000D+00, &
     -0.31673944809898013641D+00, &
      0.00000000000000000000D+00, &
      0.61877572218072185439D+00, &
      0.00000000000000000000D+00, &
     -0.61877572218072185439D+00, &
      0.00000000000000000000D+00, &
      0.37235294402955443349D+00, &
      0.00000000000000000000D+00, &
     -0.37235294402955443349D+00, &
      0.00000000000000000000D+00, &
      0.50387193570782862206D+00, &
      0.00000000000000000000D+00, &
     -0.50387193570782862206D+00, &
      0.00000000000000000000D+00, &
      0.72207836507269651527D+00, &
      0.00000000000000000000D+00, &
     -0.72207836507269651527D+00, &
      0.00000000000000000000D+00, &
      0.39403113359582681019D+00, &
      0.00000000000000000000D+00, &
     -0.39403113359582681019D+00, &
      0.00000000000000000000D+00, &
      0.11819479806124476295D+00, &
      0.00000000000000000000D+00, &
     -0.11819479806124476295D+00, &
      0.00000000000000000000D+00, &
      0.51179862044264223808D+00, &
      0.00000000000000000000D+00, &
     -0.51179862044264223808D+00, &
      0.00000000000000000000D+00, &
      0.09344056049843052492D+00, &
      0.00000000000000000000D+00, &
     -0.09344056049843052492D+00, &
      0.00000000000000000000D+00, &
      0.97786305288383967849D+00, &
      0.00000000000000000000D+00, &
     -0.97786305288383967849D+00, &
      0.00000000000000000000D+00, &
      0.26893118028825491184D+00, &
      0.00000000000000000000D+00, &
     -0.26893118028825491184D+00, &
      0.00000000000000000000D+00, &
      0.18588876717193755783D+00, &
      0.00000000000000000000D+00, &
     -0.18588876717193755783D+00, &
      0.88032062123229137818D+00, &
      0.88032062123229137818D+00, &
     -0.88032062123229137818D+00, &
     -0.88032062123229137818D+00, &
      0.68145212136761557087D+00, &
      0.68145212136761557087D+00, &
     -0.68145212136761557087D+00, &
     -0.68145212136761557087D+00, &
      0.49699036434274163065D+00, &
      0.49699036434274163065D+00, &
     -0.49699036434274163065D+00, &
     -0.49699036434274163065D+00, &
      0.30634018082228819946D+00, &
      0.30634018082228819946D+00, &
     -0.30634018082228819946D+00, &
     -0.30634018082228819946D+00, &
      0.58539966455385050725D+00, &
      0.58539966455385050725D+00, &
     -0.58539966455385050725D+00, &
     -0.58539966455385050725D+00, &
      0.77694395476943478585D+00, &
      0.77694395476943478585D+00, &
     -0.77694395476943478585D+00, &
     -0.77694395476943478585D+00, &
      0.18695368401965414829D+00, &
      0.18695368401965414829D+00, &
     -0.18695368401965414829D+00, &
     -0.18695368401965414829D+00, &
      0.64402844061291575350D+00, &
      0.64402844061291575350D+00, &
     -0.64402844061291575350D+00, &
     -0.64402844061291575350D+00, &
      0.19939684409397287479D+00, &
      0.19939684409397287479D+00, &
     -0.19939684409397287479D+00, &
     -0.19939684409397287479D+00, &
      0.49151662722959021945D+00, &
      0.49151662722959021945D+00, &
     -0.49151662722959021945D+00, &
     -0.49151662722959021945D+00, &
      0.21971883736597888626D+00, &
      0.21971883736597888626D+00, &
     -0.21971883736597888626D+00, &
     -0.21971883736597888626D+00, &
      0.38181860275292400786D+00, &
      0.38181860275292400786D+00, &
     -0.38181860275292400786D+00, &
     -0.38181860275292400786D+00, &
      0.13233574436496933768D+00, &
      0.13233574436496933768D+00, &
     -0.13233574436496933768D+00, &
     -0.13233574436496933768D+00, &
      0.26434906680031583504D+00, &
      0.26434906680031583504D+00, &
     -0.26434906680031583504D+00, &
     -0.26434906680031583504D+00, &
      0.26259016722585765136D+00, &
      0.26259016722585765136D+00, &
     -0.26259016722585765136D+00, &
     -0.26259016722585765136D+00, &
      0.63683049968859806178D+00, &
      0.63683049968859806178D+00, &
     -0.63683049968859806178D+00, &
     -0.63683049968859806178D+00, &
      0.50333230035733345087D+00, &
      0.50333230035733345087D+00, &
     -0.50333230035733345087D+00, &
     -0.50333230035733345087D+00, &
      0.43053066487698737896D+00, &
      0.43053066487698737896D+00, &
     -0.43053066487698737896D+00, &
     -0.43053066487698737896D+00, &
      0.44121177996231647489D+00, &
      0.44121177996231647489D+00, &
     -0.44121177996231647489D+00, &
     -0.44121177996231647489D+00, &
      0.12338526534910826404D+00, &
      0.12338526534910826404D+00, &
     -0.12338526534910826404D+00, &
     -0.12338526534910826404D+00, &
      0.28196466746787440805D+00, &
      0.28196466746787440805D+00, &
     -0.28196466746787440805D+00, &
     -0.28196466746787440805D+00, &
      0.12302065942509179952D+00, &
      0.12302065942509179952D+00, &
     -0.12302065942509179952D+00, &
     -0.12302065942509179952D+00, &
      0.96995340856169265376D+00, &
      0.96995340856169265376D+00, &
     -0.96995340856169265376D+00, &
     -0.96995340856169265376D+00, &
      0.15367597423015844083D+00, &
      0.15367597423015844083D+00, &
     -0.15367597423015844083D+00, &
     -0.15367597423015844083D+00, &
      0.32324916002943732130D+00, &
      0.32324916002943732130D+00, &
     -0.32324916002943732130D+00, &
     -0.32324916002943732130D+00, &
      0.07979061002341172881D+00, &
      0.07979061002341172881D+00, &
     -0.07979061002341172881D+00, &
     -0.07979061002341172881D+00, &
      0.05572751587573096521D+00, &
      0.05572751587573096521D+00, &
     -0.05572751587573096521D+00, &
     -0.05572751587573096521D+00, &
      0.80407180079417317486D+00, &
      0.80407180079417317486D+00, &
     -0.80407180079417317486D+00, &
     -0.80407180079417317486D+00, &
      0.10630794929321137066D+00, &
      0.10630794929321137066D+00, &
     -0.10630794929321137066D+00, &
     -0.10630794929321137066D+00, &
      0.14968116552008028930D+00, &
      0.14968116552008028930D+00, &
     -0.14968116552008028930D+00, &
     -0.14968116552008028930D+00, &
      0.72179802358977040999D+00, &
      0.72179802358977040999D+00, &
     -0.72179802358977040999D+00, &
     -0.72179802358977040999D+00, &
      0.64320215485058829241D+00, &
      0.81797170822487763608D+00, &
      0.64320215485058829241D+00, &
     -0.81797170822487763608D+00, &
     -0.64320215485058829241D+00, &
      0.81797170822487763608D+00, &
     -0.64320215485058829241D+00, &
     -0.81797170822487763608D+00, &
      0.16636418997995866542D+00, &
      0.45477772302529712034D+00, &
      0.16636418997995866542D+00, &
     -0.45477772302529712034D+00, &
     -0.16636418997995866542D+00, &
      0.45477772302529712034D+00, &
     -0.16636418997995866542D+00, &
     -0.45477772302529712034D+00, &
      0.64919419286301172090D+00, &
      0.22077509021471022899D+00, &
      0.64919419286301172090D+00, &
     -0.22077509021471022899D+00, &
     -0.64919419286301172090D+00, &
      0.22077509021471022899D+00, &
     -0.64919419286301172090D+00, &
     -0.22077509021471022899D+00, &
      0.26442320916257716634D+00, &
      0.85849144377065222944D+00, &
      0.26442320916257716634D+00, &
     -0.85849144377065222944D+00, &
     -0.26442320916257716634D+00, &
      0.85849144377065222944D+00, &
     -0.26442320916257716634D+00, &
     -0.85849144377065222944D+00, &
      0.45655985846318541954D+00, &
      0.19087247051624386951D+00, &
      0.45655985846318541954D+00, &
     -0.19087247051624386951D+00, &
     -0.45655985846318541954D+00, &
      0.19087247051624386951D+00, &
     -0.45655985846318541954D+00, &
     -0.19087247051624386951D+00, &
      0.66406293254753101518D+00, &
      0.37787816613860053527D+00, &
      0.66406293254753101518D+00, &
     -0.37787816613860053527D+00, &
     -0.66406293254753101518D+00, &
      0.37787816613860053527D+00, &
     -0.66406293254753101518D+00, &
     -0.37787816613860053527D+00, &
      0.17283329178347239807D+00, &
      0.89395556646395457623D+00, &
      0.17283329178347239807D+00, &
     -0.89395556646395457623D+00, &
     -0.17283329178347239807D+00, &
      0.89395556646395457623D+00, &
     -0.17283329178347239807D+00, &
     -0.89395556646395457623D+00, &
      0.05003413197772252352D+00, &
      0.40991415053935187363D+00, &
      0.05003413197772252352D+00, &
     -0.40991415053935187363D+00, &
     -0.05003413197772252352D+00, &
      0.40991415053935187363D+00, &
     -0.05003413197772252352D+00, &
     -0.40991415053935187363D+00, &
      0.10662173986061880548D+00, &
      0.93198108654420652730D+00, &
      0.10662173986061880548D+00, &
     -0.93198108654420652730D+00, &
     -0.10662173986061880548D+00, &
      0.93198108654420652730D+00, &
     -0.10662173986061880548D+00, &
     -0.93198108654420652730D+00, &
      0.84993195986801317598D+00, &
      0.28018220607206512085D+00, &
      0.84993195986801317598D+00, &
     -0.28018220607206512085D+00, &
     -0.84993195986801317598D+00, &
      0.28018220607206512085D+00, &
     -0.84993195986801317598D+00, &
     -0.28018220607206512085D+00, &
      0.46201429744731292715D+00, &
      0.57728165413613941048D+00, &
      0.46201429744731292715D+00, &
     -0.57728165413613941048D+00, &
     -0.46201429744731292715D+00, &
      0.57728165413613941048D+00, &
     -0.46201429744731292715D+00, &
     -0.57728165413613941048D+00, &
      0.61756791755060624904D+00, &
      0.84451580808527704214D+00, &
      0.61756791755060624904D+00, &
     -0.84451580808527704214D+00, &
     -0.61756791755060624904D+00, &
      0.84451580808527704214D+00, &
     -0.61756791755060624904D+00, &
     -0.84451580808527704214D+00, &
      0.23348533584682021336D+00, &
      0.52941572893811861267D+00, &
      0.23348533584682021336D+00, &
     -0.52941572893811861267D+00, &
     -0.23348533584682021336D+00, &
      0.52941572893811861267D+00, &
     -0.23348533584682021336D+00, &
     -0.52941572893811861267D+00, &
      0.67941786279817284466D+00, &
      0.92827521801894674613D+00, &
      0.67941786279817284466D+00, &
     -0.92827521801894674613D+00, &
     -0.67941786279817284466D+00, &
      0.92827521801894674613D+00, &
     -0.67941786279817284466D+00, &
     -0.92827521801894674613D+00, &
      0.46144251438466377113D+00, &
      0.96836516384095050469D+00, &
      0.46144251438466377113D+00, &
     -0.96836516384095050469D+00, &
     -0.46144251438466377113D+00, &
      0.96836516384095050469D+00, &
     -0.46144251438466377113D+00, &
     -0.96836516384095050469D+00, &
      0.31430414794159466929D+00, &
      0.16744571840030086918D+00, &
      0.31430414794159466929D+00, &
     -0.16744571840030086918D+00, &
     -0.31430414794159466929D+00, &
      0.16744571840030086918D+00, &
     -0.31430414794159466929D+00, &
     -0.16744571840030086918D+00, &
      0.73988842925403419670D+00, &
      0.49232718454496093852D+00, &
      0.73988842925403419670D+00, &
     -0.49232718454496093852D+00, &
     -0.73988842925403419670D+00, &
      0.49232718454496093852D+00, &
     -0.73988842925403419670D+00, &
     -0.49232718454496093852D+00, &
      0.54061065617200576572D+00, &
      0.85659974228270108210D+00, &
      0.54061065617200576572D+00, &
     -0.85659974228270108210D+00, &
     -0.54061065617200576572D+00, &
      0.85659974228270108210D+00, &
     -0.54061065617200576572D+00, &
     -0.85659974228270108210D+00, &
      0.57531393529153806998D+00, &
      0.26702277361260906563D+00, &
      0.57531393529153806998D+00, &
     -0.26702277361260906563D+00, &
     -0.57531393529153806998D+00, &
      0.26702277361260906563D+00, &
     -0.57531393529153806998D+00, &
     -0.26702277361260906563D+00, &
      0.71731018735152574095D+00, &
      0.41506457665342233465D+00, &
      0.71731018735152574095D+00, &
     -0.41506457665342233465D+00, &
     -0.71731018735152574095D+00, &
      0.41506457665342233465D+00, &
     -0.71731018735152574095D+00, &
     -0.41506457665342233465D+00, &
      0.81619080186595505122D+00, &
      0.94466239579161759288D+00, &
      0.81619080186595505122D+00, &
     -0.94466239579161759288D+00, &
     -0.81619080186595505122D+00, &
      0.94466239579161759288D+00, &
     -0.81619080186595505122D+00, &
     -0.94466239579161759288D+00, &
      0.31179783227055996031D+00, &
      0.71450522624081691525D+00, &
      0.31179783227055996031D+00, &
     -0.71450522624081691525D+00, &
     -0.31179783227055996031D+00, &
      0.71450522624081691525D+00, &
     -0.31179783227055996031D+00, &
     -0.71450522624081691525D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
      0.96891318773147550036D+00, &
      0.11218085315386774892D+00, &
      0.11218085315386774892D+00, &
      0.11218085315386774892D+00, &
      0.11218085315386774892D+00, &
      0.20952037314221555464D+00, &
      0.20952037314221555464D+00, &
      0.20952037314221555464D+00, &
      0.20952037314221555464D+00, &
      0.45043339209409832824D+00, &
      0.45043339209409832824D+00, &
      0.45043339209409832824D+00, &
      0.45043339209409832824D+00, &
      0.36713809241502448621D+00, &
      0.36713809241502448621D+00, &
      0.36713809241502448621D+00, &
      0.36713809241502448621D+00, &
      0.51496945555138884387D+00, &
      0.51496945555138884387D+00, &
      0.51496945555138884387D+00, &
      0.51496945555138884387D+00, &
      0.24658740648660407158D+00, &
      0.24658740648660407158D+00, &
      0.24658740648660407158D+00, &
      0.24658740648660407158D+00, &
      0.03347288009925680763D+00, &
      0.03347288009925680763D+00, &
      0.03347288009925680763D+00, &
      0.03347288009925680763D+00, &
      0.17604915946010571415D+00, &
      0.17604915946010571415D+00, &
      0.17604915946010571415D+00, &
      0.17604915946010571415D+00, &
      0.12627670900302459533D+00, &
      0.12627670900302459533D+00, &
      0.12627670900302459533D+00, &
      0.12627670900302459533D+00, &
      0.01146186842907397259D+00, &
      0.01146186842907397259D+00, &
      0.01146186842907397259D+00, &
      0.01146186842907397259D+00, &
      0.82865571811563620841D+00, &
      0.82865571811563620841D+00, &
      0.82865571811563620841D+00, &
      0.82865571811563620841D+00, &
      0.00411477663570613411D+00, &
      0.00411477663570613411D+00, &
      0.00411477663570613411D+00, &
      0.00411477663570613411D+00, &
      0.66431760085930258164D+00, &
      0.66431760085930258164D+00, &
      0.66431760085930258164D+00, &
      0.66431760085930258164D+00, &
      0.80714531646844267510D+00, &
      0.80714531646844267510D+00, &
      0.80714531646844267510D+00, &
      0.80714531646844267510D+00, &
      0.06759778909905436728D+00, &
      0.06759778909905436728D+00, &
      0.06759778909905436728D+00, &
      0.06759778909905436728D+00, &
      0.29403431132680074578D+00, &
      0.29403431132680074578D+00, &
      0.29403431132680074578D+00, &
      0.29403431132680074578D+00, &
      0.49665508093160370962D+00, &
      0.49665508093160370962D+00, &
      0.49665508093160370962D+00, &
      0.49665508093160370962D+00, &
      0.27911320549015106174D+00, &
      0.27911320549015106174D+00, &
      0.27911320549015106174D+00, &
      0.27911320549015106174D+00, &
      0.07876704807666434771D+00, &
      0.07876704807666434771D+00, &
      0.07876704807666434771D+00, &
      0.07876704807666434771D+00, &
      0.04825635484814582571D+00, &
      0.04825635484814582571D+00, &
      0.04825635484814582571D+00, &
      0.04825635484814582571D+00, &
      0.00881247486166711683D+00, &
      0.00881247486166711683D+00, &
      0.00881247486166711683D+00, &
      0.00881247486166711683D+00, &
      0.02005872107121197886D+00, &
      0.02005872107121197886D+00, &
      0.02005872107121197886D+00, &
      0.02005872107121197886D+00, &
      0.75639094968056375112D+00, &
      0.75639094968056375112D+00, &
      0.75639094968056375112D+00, &
      0.75639094968056375112D+00, &
      0.15155505938254110188D+00, &
      0.15155505938254110188D+00, &
      0.15155505938254110188D+00, &
      0.15155505938254110188D+00, &
      0.60842758819371955958D+00, &
      0.60842758819371955958D+00, &
      0.60842758819371955958D+00, &
      0.60842758819371955958D+00, &
      0.52569545459372502005D+00, &
      0.52569545459372502005D+00, &
      0.52569545459372502005D+00, &
      0.52569545459372502005D+00, &
      0.40281226004252063122D+00, &
      0.40281226004252063122D+00, &
      0.40281226004252063122D+00, &
      0.40281226004252063122D+00, &
      0.04148931197347632133D+00, &
      0.04148931197347632133D+00, &
      0.04148931197347632133D+00, &
      0.04148931197347632133D+00, &
      0.14661700764838450639D+00, &
      0.14661700764838450639D+00, &
      0.14661700764838450639D+00, &
      0.14661700764838450639D+00, &
      0.25661597420550030790D+00, &
      0.25661597420550030790D+00, &
      0.25661597420550030790D+00, &
      0.25661597420550030790D+00, &
      0.25552367243392443141D+00, &
      0.25552367243392443141D+00, &
      0.25552367243392443141D+00, &
      0.25552367243392443141D+00, &
      0.36981856719481731588D+00, &
      0.36981856719481731588D+00, &
      0.36981856719481731588D+00, &
      0.36981856719481731588D+00, &
      0.01554708555805913925D+00, &
      0.01554708555805913925D+00, &
      0.01554708555805913925D+00, &
      0.01554708555805913925D+00, &
      0.04951092940950102550D+00, &
      0.04951092940950102550D+00, &
      0.04951092940950102550D+00, &
      0.04951092940950102550D+00, &
      0.45729924127195364925D+00, &
      0.45729924127195364925D+00, &
      0.45729924127195364925D+00, &
      0.45729924127195364925D+00, &
      0.25865392811548787444D+00, &
      0.25865392811548787444D+00, &
      0.25865392811548787444D+00, &
      0.25865392811548787444D+00, &
      0.01471561063742727986D+00, &
      0.01471561063742727986D+00, &
      0.01471561063742727986D+00, &
      0.01471561063742727986D+00, &
      0.84302938509543134948D+00, &
      0.84302938509543134948D+00, &
      0.84302938509543134948D+00, &
      0.84302938509543134948D+00, &
      0.64221316316195331542D+00, &
      0.64221316316195331542D+00, &
      0.64221316316195331542D+00, &
      0.64221316316195331542D+00, &
      0.70815468949755844275D+00, &
      0.70815468949755844275D+00, &
      0.70815468949755844275D+00, &
      0.70815468949755844275D+00, &
      0.90756616766447684164D+00, &
      0.90756616766447684164D+00, &
      0.90756616766447684164D+00, &
      0.90756616766447684164D+00, &
      0.15801098912608321778D+00, &
      0.15801098912608321778D+00, &
      0.15801098912608321778D+00, &
      0.15801098912608321778D+00, &
      0.55864014423757180072D+00, &
      0.55864014423757180072D+00, &
      0.55864014423757180072D+00, &
      0.55864014423757180072D+00, &
      0.76553411212536370822D+00, &
      0.76553411212536370822D+00, &
      0.76553411212536370822D+00, &
      0.76553411212536370822D+00, &
      0.12038628177040819334D+00, &
      0.12038628177040819334D+00, &
      0.12038628177040819334D+00, &
      0.12038628177040819334D+00, &
      0.16540042455816997280D+00, &
      0.16540042455816997280D+00, &
      0.16540042455816997280D+00, &
      0.16540042455816997280D+00, &
      0.16540042455816997280D+00, &
      0.16540042455816997280D+00, &
      0.16540042455816997280D+00, &
      0.16540042455816997280D+00, &
      0.07191647566286844817D+00, &
      0.07191647566286844817D+00, &
      0.07191647566286844817D+00, &
      0.07191647566286844817D+00, &
      0.07191647566286844817D+00, &
      0.07191647566286844817D+00, &
      0.07191647566286844817D+00, &
      0.07191647566286844817D+00, &
      0.24369366228681560438D+00, &
      0.24369366228681560438D+00, &
      0.24369366228681560438D+00, &
      0.24369366228681560438D+00, &
      0.24369366228681560438D+00, &
      0.24369366228681560438D+00, &
      0.24369366228681560438D+00, &
      0.24369366228681560438D+00, &
      0.03997911405124788403D+00, &
      0.03997911405124788403D+00, &
      0.03997911405124788403D+00, &
      0.03997911405124788403D+00, &
      0.03997911405124788403D+00, &
      0.03997911405124788403D+00, &
      0.03997911405124788403D+00, &
      0.03997911405124788403D+00, &
      0.52403522696450499652D+00, &
      0.52403522696450499652D+00, &
      0.52403522696450499652D+00, &
      0.52403522696450499652D+00, &
      0.52403522696450499652D+00, &
      0.52403522696450499652D+00, &
      0.52403522696450499652D+00, &
      0.52403522696450499652D+00, &
      0.06332251756914623886D+00, &
      0.06332251756914623886D+00, &
      0.06332251756914623886D+00, &
      0.06332251756914623886D+00, &
      0.06332251756914623886D+00, &
      0.06332251756914623886D+00, &
      0.06332251756914623886D+00, &
      0.06332251756914623886D+00, &
      0.00627167365301165187D+00, &
      0.00627167365301165187D+00, &
      0.00627167365301165187D+00, &
      0.00627167365301165187D+00, &
      0.00627167365301165187D+00, &
      0.00627167365301165187D+00, &
      0.00627167365301165187D+00, &
      0.00627167365301165187D+00, &
      0.34612550470442166040D+00, &
      0.34612550470442166040D+00, &
      0.34612550470442166040D+00, &
      0.34612550470442166040D+00, &
      0.34612550470442166040D+00, &
      0.34612550470442166040D+00, &
      0.34612550470442166040D+00, &
      0.34612550470442166040D+00, &
      0.04934944134498662344D+00, &
      0.04934944134498662344D+00, &
      0.04934944134498662344D+00, &
      0.04934944134498662344D+00, &
      0.04934944134498662344D+00, &
      0.04934944134498662344D+00, &
      0.04934944134498662344D+00, &
      0.04934944134498662344D+00, &
      0.12050259192784877615D+00, &
      0.12050259192784877615D+00, &
      0.12050259192784877615D+00, &
      0.12050259192784877615D+00, &
      0.12050259192784877615D+00, &
      0.12050259192784877615D+00, &
      0.12050259192784877615D+00, &
      0.12050259192784877615D+00, &
      0.39428120571642771841D+00, &
      0.39428120571642771841D+00, &
      0.39428120571642771841D+00, &
      0.39428120571642771841D+00, &
      0.39428120571642771841D+00, &
      0.39428120571642771841D+00, &
      0.39428120571642771841D+00, &
      0.39428120571642771841D+00, &
      0.00838019453790634500D+00, &
      0.00838019453790634500D+00, &
      0.00838019453790634500D+00, &
      0.00838019453790634500D+00, &
      0.00838019453790634500D+00, &
      0.00838019453790634500D+00, &
      0.00838019453790634500D+00, &
      0.00838019453790634500D+00, &
      0.37991521751516699190D+00, &
      0.37991521751516699190D+00, &
      0.37991521751516699190D+00, &
      0.37991521751516699190D+00, &
      0.37991521751516699190D+00, &
      0.37991521751516699190D+00, &
      0.37991521751516699190D+00, &
      0.37991521751516699190D+00, &
      0.06658426856019065976D+00, &
      0.06658426856019065976D+00, &
      0.06658426856019065976D+00, &
      0.06658426856019065976D+00, &
      0.06658426856019065976D+00, &
      0.06658426856019065976D+00, &
      0.06658426856019065976D+00, &
      0.06658426856019065976D+00, &
      0.01360996700708025364D+00, &
      0.01360996700708025364D+00, &
      0.01360996700708025364D+00, &
      0.01360996700708025364D+00, &
      0.01360996700708025364D+00, &
      0.01360996700708025364D+00, &
      0.01360996700708025364D+00, &
      0.01360996700708025364D+00, &
      0.67377911630059594827D+00, &
      0.67377911630059594827D+00, &
      0.67377911630059594827D+00, &
      0.67377911630059594827D+00, &
      0.67377911630059594827D+00, &
      0.67377911630059594827D+00, &
      0.67377911630059594827D+00, &
      0.67377911630059594827D+00, &
      0.15122396864513820702D+00, &
      0.15122396864513820702D+00, &
      0.15122396864513820702D+00, &
      0.15122396864513820702D+00, &
      0.15122396864513820702D+00, &
      0.15122396864513820702D+00, &
      0.15122396864513820702D+00, &
      0.15122396864513820702D+00, &
      0.05960704457999711076D+00, &
      0.05960704457999711076D+00, &
      0.05960704457999711076D+00, &
      0.05960704457999711076D+00, &
      0.05960704457999711076D+00, &
      0.05960704457999711076D+00, &
      0.05960704457999711076D+00, &
      0.05960704457999711076D+00, &
      0.15557252093495238521D+00, &
      0.15557252093495238521D+00, &
      0.15557252093495238521D+00, &
      0.15557252093495238521D+00, &
      0.15557252093495238521D+00, &
      0.15557252093495238521D+00, &
      0.15557252093495238521D+00, &
      0.15557252093495238521D+00, &
      0.26467312835339185106D+00, &
      0.26467312835339185106D+00, &
      0.26467312835339185106D+00, &
      0.26467312835339185106D+00, &
      0.26467312835339185106D+00, &
      0.26467312835339185106D+00, &
      0.26467312835339185106D+00, &
      0.26467312835339185106D+00, &
      0.01193065095103094247D+00, &
      0.01193065095103094247D+00, &
      0.01193065095103094247D+00, &
      0.01193065095103094247D+00, &
      0.01193065095103094247D+00, &
      0.01193065095103094247D+00, &
      0.01193065095103094247D+00, &
      0.01193065095103094247D+00, &
      0.00729413025649284976D+00, &
      0.00729413025649284976D+00, &
      0.00729413025649284976D+00, &
      0.00729413025649284976D+00, &
      0.00729413025649284976D+00, &
      0.00729413025649284976D+00, &
      0.00729413025649284976D+00, &
      0.00729413025649284976D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
      0.00013400520132544013D+00, &
      0.00653282923653435328D+00, &
      0.00653282923653435328D+00, &
      0.00653282923653435328D+00, &
      0.00653282923653435328D+00, &
      0.00246267782986044981D+00, &
      0.00246267782986044981D+00, &
      0.00246267782986044981D+00, &
      0.00246267782986044981D+00, &
      0.00236266319877590381D+00, &
      0.00236266319877590381D+00, &
      0.00236266319877590381D+00, &
      0.00236266319877590381D+00, &
      0.00241284626019594179D+00, &
      0.00241284626019594179D+00, &
      0.00241284626019594179D+00, &
      0.00241284626019594179D+00, &
      0.00563939918988414120D+00, &
      0.00563939918988414120D+00, &
      0.00563939918988414120D+00, &
      0.00563939918988414120D+00, &
      0.00530879464537392220D+00, &
      0.00530879464537392220D+00, &
      0.00530879464537392220D+00, &
      0.00530879464537392220D+00, &
      0.00330320967760970065D+00, &
      0.00330320967760970065D+00, &
      0.00330320967760970065D+00, &
      0.00330320967760970065D+00, &
      0.00715709687293028703D+00, &
      0.00715709687293028703D+00, &
      0.00715709687293028703D+00, &
      0.00715709687293028703D+00, &
      0.00356020900079922868D+00, &
      0.00356020900079922868D+00, &
      0.00356020900079922868D+00, &
      0.00356020900079922868D+00, &
      0.00286908910615377545D+00, &
      0.00286908910615377545D+00, &
      0.00286908910615377545D+00, &
      0.00286908910615377545D+00, &
      0.00119238795542696220D+00, &
      0.00119238795542696220D+00, &
      0.00119238795542696220D+00, &
      0.00119238795542696220D+00, &
      0.00023570897730274608D+00, &
      0.00023570897730274608D+00, &
      0.00023570897730274608D+00, &
      0.00023570897730274608D+00, &
      0.00332415202652309181D+00, &
      0.00332415202652309181D+00, &
      0.00332415202652309181D+00, &
      0.00332415202652309181D+00, &
      0.00064424498387906620D+00, &
      0.00064424498387906620D+00, &
      0.00064424498387906620D+00, &
      0.00064424498387906620D+00, &
      0.00093040828397316332D+00, &
      0.00093040828397316332D+00, &
      0.00093040828397316332D+00, &
      0.00093040828397316332D+00, &
      0.00061531263855177734D+00, &
      0.00061531263855177734D+00, &
      0.00061531263855177734D+00, &
      0.00061531263855177734D+00, &
      0.00028907482286127771D+00, &
      0.00028907482286127771D+00, &
      0.00028907482286127771D+00, &
      0.00028907482286127771D+00, &
      0.00957928057788382838D+00, &
      0.00957928057788382838D+00, &
      0.00957928057788382838D+00, &
      0.00957928057788382838D+00, &
      0.00256243378291242307D+00, &
      0.00256243378291242307D+00, &
      0.00256243378291242307D+00, &
      0.00256243378291242307D+00, &
      0.00197009710848232650D+00, &
      0.00197009710848232650D+00, &
      0.00197009710848232650D+00, &
      0.00197009710848232650D+00, &
      0.00224853801623185919D+00, &
      0.00224853801623185919D+00, &
      0.00224853801623185919D+00, &
      0.00224853801623185919D+00, &
      0.00178378503667077598D+00, &
      0.00178378503667077598D+00, &
      0.00178378503667077598D+00, &
      0.00178378503667077598D+00, &
      0.00046154884778758861D+00, &
      0.00046154884778758861D+00, &
      0.00046154884778758861D+00, &
      0.00046154884778758861D+00, &
      0.00383567352584877853D+00, &
      0.00383567352584877853D+00, &
      0.00383567352584877853D+00, &
      0.00383567352584877853D+00, &
      0.00424323389268952030D+00, &
      0.00424323389268952030D+00, &
      0.00424323389268952030D+00, &
      0.00424323389268952030D+00, &
      0.00297709850925886605D+00, &
      0.00297709850925886605D+00, &
      0.00297709850925886605D+00, &
      0.00297709850925886605D+00, &
      0.00750097658421174009D+00, &
      0.00750097658421174009D+00, &
      0.00750097658421174009D+00, &
      0.00750097658421174009D+00, &
      0.00208064231905125911D+00, &
      0.00208064231905125911D+00, &
      0.00208064231905125911D+00, &
      0.00208064231905125911D+00, &
      0.00676296550091519920D+00, &
      0.00676296550091519920D+00, &
      0.00676296550091519920D+00, &
      0.00676296550091519920D+00, &
      0.00300587484767854103D+00, &
      0.00300587484767854103D+00, &
      0.00300587484767854103D+00, &
      0.00300587484767854103D+00, &
      0.00569473175276393340D+00, &
      0.00569473175276393340D+00, &
      0.00569473175276393340D+00, &
      0.00569473175276393340D+00, &
      0.00329345228074993208D+00, &
      0.00329345228074993208D+00, &
      0.00329345228074993208D+00, &
      0.00329345228074993208D+00, &
      0.00319343074693561636D+00, &
      0.00319343074693561636D+00, &
      0.00319343074693561636D+00, &
      0.00319343074693561636D+00, &
      0.00300785963957934336D+00, &
      0.00300785963957934336D+00, &
      0.00300785963957934336D+00, &
      0.00300785963957934336D+00, &
      0.00621018081426024000D+00, &
      0.00621018081426024000D+00, &
      0.00621018081426024000D+00, &
      0.00621018081426024000D+00, &
      0.00690380129393601029D+00, &
      0.00690380129393601029D+00, &
      0.00690380129393601029D+00, &
      0.00690380129393601029D+00, &
      0.00010061838637633462D+00, &
      0.00010061838637633462D+00, &
      0.00010061838637633462D+00, &
      0.00010061838637633462D+00, &
      0.00013563078282338245D+00, &
      0.00013563078282338245D+00, &
      0.00013563078282338245D+00, &
      0.00013563078282338245D+00, &
      0.00077893767219553941D+00, &
      0.00077893767219553941D+00, &
      0.00077893767219553941D+00, &
      0.00077893767219553941D+00, &
      0.00267441470336475876D+00, &
      0.00267441470336475876D+00, &
      0.00267441470336475876D+00, &
      0.00267441470336475876D+00, &
      0.00047754594130653132D+00, &
      0.00047754594130653132D+00, &
      0.00047754594130653132D+00, &
      0.00047754594130653132D+00, &
      0.00074066104011807129D+00, &
      0.00074066104011807129D+00, &
      0.00074066104011807129D+00, &
      0.00074066104011807129D+00, &
      0.00529137043945265752D+00, &
      0.00529137043945265752D+00, &
      0.00529137043945265752D+00, &
      0.00529137043945265752D+00, &
      0.00150066461105524858D+00, &
      0.00150066461105524858D+00, &
      0.00150066461105524858D+00, &
      0.00150066461105524858D+00, &
      0.00262627399936291505D+00, &
      0.00262627399936291505D+00, &
      0.00262627399936291505D+00, &
      0.00262627399936291505D+00, &
      0.00106543689917693658D+00, &
      0.00106543689917693658D+00, &
      0.00106543689917693658D+00, &
      0.00106543689917693658D+00, &
      0.00106543689917693658D+00, &
      0.00106543689917693658D+00, &
      0.00106543689917693658D+00, &
      0.00106543689917693658D+00, &
      0.00608072256149797419D+00, &
      0.00608072256149797419D+00, &
      0.00608072256149797419D+00, &
      0.00608072256149797419D+00, &
      0.00608072256149797419D+00, &
      0.00608072256149797419D+00, &
      0.00608072256149797419D+00, &
      0.00608072256149797419D+00, &
      0.00531002301058563358D+00, &
      0.00531002301058563358D+00, &
      0.00531002301058563358D+00, &
      0.00531002301058563358D+00, &
      0.00531002301058563358D+00, &
      0.00531002301058563358D+00, &
      0.00531002301058563358D+00, &
      0.00531002301058563358D+00, &
      0.00202003342483236250D+00, &
      0.00202003342483236250D+00, &
      0.00202003342483236250D+00, &
      0.00202003342483236250D+00, &
      0.00202003342483236250D+00, &
      0.00202003342483236250D+00, &
      0.00202003342483236250D+00, &
      0.00202003342483236250D+00, &
      0.00216470393083745358D+00, &
      0.00216470393083745358D+00, &
      0.00216470393083745358D+00, &
      0.00216470393083745358D+00, &
      0.00216470393083745358D+00, &
      0.00216470393083745358D+00, &
      0.00216470393083745358D+00, &
      0.00216470393083745358D+00, &
      0.00435330698915271062D+00, &
      0.00435330698915271062D+00, &
      0.00435330698915271062D+00, &
      0.00435330698915271062D+00, &
      0.00435330698915271062D+00, &
      0.00435330698915271062D+00, &
      0.00435330698915271062D+00, &
      0.00435330698915271062D+00, &
      0.00068634034778216149D+00, &
      0.00068634034778216149D+00, &
      0.00068634034778216149D+00, &
      0.00068634034778216149D+00, &
      0.00068634034778216149D+00, &
      0.00068634034778216149D+00, &
      0.00068634034778216149D+00, &
      0.00068634034778216149D+00, &
      0.00401906783115373901D+00, &
      0.00401906783115373901D+00, &
      0.00401906783115373901D+00, &
      0.00401906783115373901D+00, &
      0.00401906783115373901D+00, &
      0.00401906783115373901D+00, &
      0.00401906783115373901D+00, &
      0.00401906783115373901D+00, &
      0.00079127284739705225D+00, &
      0.00079127284739705225D+00, &
      0.00079127284739705225D+00, &
      0.00079127284739705225D+00, &
      0.00079127284739705225D+00, &
      0.00079127284739705225D+00, &
      0.00079127284739705225D+00, &
      0.00079127284739705225D+00, &
      0.00244671052896765998D+00, &
      0.00244671052896765998D+00, &
      0.00244671052896765998D+00, &
      0.00244671052896765998D+00, &
      0.00244671052896765998D+00, &
      0.00244671052896765998D+00, &
      0.00244671052896765998D+00, &
      0.00244671052896765998D+00, &
      0.00195038190205695967D+00, &
      0.00195038190205695967D+00, &
      0.00195038190205695967D+00, &
      0.00195038190205695967D+00, &
      0.00195038190205695967D+00, &
      0.00195038190205695967D+00, &
      0.00195038190205695967D+00, &
      0.00195038190205695967D+00, &
      0.00110591593882085814D+00, &
      0.00110591593882085814D+00, &
      0.00110591593882085814D+00, &
      0.00110591593882085814D+00, &
      0.00110591593882085814D+00, &
      0.00110591593882085814D+00, &
      0.00110591593882085814D+00, &
      0.00110591593882085814D+00, &
      0.00503011482374766256D+00, &
      0.00503011482374766256D+00, &
      0.00503011482374766256D+00, &
      0.00503011482374766256D+00, &
      0.00503011482374766256D+00, &
      0.00503011482374766256D+00, &
      0.00503011482374766256D+00, &
      0.00503011482374766256D+00, &
      0.00065967081570031318D+00, &
      0.00065967081570031318D+00, &
      0.00065967081570031318D+00, &
      0.00065967081570031318D+00, &
      0.00065967081570031318D+00, &
      0.00065967081570031318D+00, &
      0.00065967081570031318D+00, &
      0.00065967081570031318D+00, &
      0.00061321639781390126D+00, &
      0.00061321639781390126D+00, &
      0.00061321639781390126D+00, &
      0.00061321639781390126D+00, &
      0.00061321639781390126D+00, &
      0.00061321639781390126D+00, &
      0.00061321639781390126D+00, &
      0.00061321639781390126D+00, &
      0.00089657142446607598D+00, &
      0.00089657142446607598D+00, &
      0.00089657142446607598D+00, &
      0.00089657142446607598D+00, &
      0.00089657142446607598D+00, &
      0.00089657142446607598D+00, &
      0.00089657142446607598D+00, &
      0.00089657142446607598D+00, &
      0.00363983560477118984D+00, &
      0.00363983560477118984D+00, &
      0.00363983560477118984D+00, &
      0.00363983560477118984D+00, &
      0.00363983560477118984D+00, &
      0.00363983560477118984D+00, &
      0.00363983560477118984D+00, &
      0.00363983560477118984D+00, &
      0.00211128919028042220D+00, &
      0.00211128919028042220D+00, &
      0.00211128919028042220D+00, &
      0.00211128919028042220D+00, &
      0.00211128919028042220D+00, &
      0.00211128919028042220D+00, &
      0.00211128919028042220D+00, &
      0.00211128919028042220D+00, &
      0.00570062224543081177D+00, &
      0.00570062224543081177D+00, &
      0.00570062224543081177D+00, &
      0.00570062224543081177D+00, &
      0.00570062224543081177D+00, &
      0.00570062224543081177D+00, &
      0.00570062224543081177D+00, &
      0.00570062224543081177D+00, &
      0.00195532776601042716D+00, &
      0.00195532776601042716D+00, &
      0.00195532776601042716D+00, &
      0.00195532776601042716D+00, &
      0.00195532776601042716D+00, &
      0.00195532776601042716D+00, &
      0.00195532776601042716D+00, &
      0.00195532776601042716D+00, &
      0.00062384953506176009D+00, &
      0.00062384953506176009D+00, &
      0.00062384953506176009D+00, &
      0.00062384953506176009D+00, &
      0.00062384953506176009D+00, &
      0.00062384953506176009D+00, &
      0.00062384953506176009D+00, &
      0.00062384953506176009D+00, &
      0.00151792165402074965D+00, &
      0.00151792165402074965D+00, &
      0.00151792165402074965D+00, &
      0.00151792165402074965D+00, &
      0.00151792165402074965D+00, &
      0.00151792165402074965D+00, &
      0.00151792165402074965D+00, &
      0.00151792165402074965D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule19 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule19() returns the pyramid quadrature rule of precision 19.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 April 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 121, Number 11, pages 2418-2436, 15 June 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real ( kind = rk ) x[n], y[n], z[n]: the coordinates of
!    quadrature points.
!
!    real ( kind = rk ) w[n]: the quadrature weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 418

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.37365721634528792361D+00, &
      0.00000000000000000000D+00, &
     -0.37365721634528792361D+00, &
      0.00000000000000000000D+00, &
      0.43752827030428453892D+00, &
      0.00000000000000000000D+00, &
     -0.43752827030428453892D+00, &
      0.00000000000000000000D+00, &
      0.55782756293242685697D+00, &
      0.00000000000000000000D+00, &
     -0.55782756293242685697D+00, &
      0.00000000000000000000D+00, &
      0.27275879229149846417D+00, &
      0.00000000000000000000D+00, &
     -0.27275879229149846417D+00, &
      0.00000000000000000000D+00, &
      0.28872411753825932834D+00, &
      0.00000000000000000000D+00, &
     -0.28872411753825932834D+00, &
      0.00000000000000000000D+00, &
      0.11277465623250763904D+00, &
      0.00000000000000000000D+00, &
     -0.11277465623250763904D+00, &
      0.00000000000000000000D+00, &
      0.65470230415224206499D+00, &
      0.00000000000000000000D+00, &
     -0.65470230415224206499D+00, &
      0.00000000000000000000D+00, &
      0.22295101282662088682D+00, &
      0.00000000000000000000D+00, &
     -0.22295101282662088682D+00, &
      0.00000000000000000000D+00, &
      0.35176689121153298379D+00, &
      0.00000000000000000000D+00, &
     -0.35176689121153298379D+00, &
      0.00000000000000000000D+00, &
      0.20697735413665671600D+00, &
      0.00000000000000000000D+00, &
     -0.20697735413665671600D+00, &
      0.00000000000000000000D+00, &
      0.93111189027726670719D+00, &
      0.00000000000000000000D+00, &
     -0.93111189027726670719D+00, &
      0.00000000000000000000D+00, &
      0.84303358445951348532D+00, &
      0.00000000000000000000D+00, &
     -0.84303358445951348532D+00, &
      0.00000000000000000000D+00, &
      0.07303948013428374075D+00, &
      0.00000000000000000000D+00, &
     -0.07303948013428374075D+00, &
      0.00000000000000000000D+00, &
      0.67807351338628540915D+00, &
      0.00000000000000000000D+00, &
     -0.67807351338628540915D+00, &
      0.00000000000000000000D+00, &
      0.88198482628664209848D+00, &
      0.00000000000000000000D+00, &
     -0.88198482628664209848D+00, &
      0.00000000000000000000D+00, &
      0.38587811686103651310D+00, &
      0.00000000000000000000D+00, &
     -0.38587811686103651310D+00, &
      0.00000000000000000000D+00, &
      0.13708884737518356456D+00, &
      0.00000000000000000000D+00, &
     -0.13708884737518356456D+00, &
      0.00000000000000000000D+00, &
      0.52558521108268729805D+00, &
     -0.52558521108268729805D+00, &
      0.52558521108268729805D+00, &
     -0.52558521108268729805D+00, &
      0.45301642570589800707D+00, &
     -0.45301642570589800707D+00, &
      0.45301642570589800707D+00, &
     -0.45301642570589800707D+00, &
      0.67125116320175959306D+00, &
     -0.67125116320175959306D+00, &
      0.67125116320175959306D+00, &
     -0.67125116320175959306D+00, &
      0.36951418649168821240D+00, &
     -0.36951418649168821240D+00, &
      0.36951418649168821240D+00, &
     -0.36951418649168821240D+00, &
      0.28838223435586252119D+00, &
     -0.28838223435586252119D+00, &
      0.28838223435586252119D+00, &
     -0.28838223435586252119D+00, &
      0.82857238389465803774D+00, &
     -0.82857238389465803774D+00, &
      0.82857238389465803774D+00, &
     -0.82857238389465803774D+00, &
      0.67246070085447195996D+00, &
     -0.67246070085447195996D+00, &
      0.67246070085447195996D+00, &
     -0.67246070085447195996D+00, &
      0.48736689338962557727D+00, &
     -0.48736689338962557727D+00, &
      0.48736689338962557727D+00, &
     -0.48736689338962557727D+00, &
      0.95036644315079843448D+00, &
     -0.95036644315079843448D+00, &
      0.95036644315079843448D+00, &
     -0.95036644315079843448D+00, &
      0.09132353058612172059D+00, &
     -0.09132353058612172059D+00, &
      0.09132353058612172059D+00, &
     -0.09132353058612172059D+00, &
      0.15520695551416591185D+00, &
     -0.15520695551416591185D+00, &
      0.15520695551416591185D+00, &
     -0.15520695551416591185D+00, &
      0.05052806321356113906D+00, &
     -0.05052806321356113906D+00, &
      0.05052806321356113906D+00, &
     -0.05052806321356113906D+00, &
      0.39781392804166187949D+00, &
     -0.39781392804166187949D+00, &
      0.39781392804166187949D+00, &
     -0.39781392804166187949D+00, &
      0.37607913140641990868D+00, &
     -0.37607913140641990868D+00, &
      0.37607913140641990868D+00, &
     -0.37607913140641990868D+00, &
      0.16436382679828262510D+00, &
     -0.16436382679828262510D+00, &
      0.16436382679828262510D+00, &
     -0.16436382679828262510D+00, &
      0.34688925556861166521D+00, &
     -0.34688925556861166521D+00, &
      0.34688925556861166521D+00, &
     -0.34688925556861166521D+00, &
      0.11913325793719439782D+00, &
     -0.11913325793719439782D+00, &
      0.11913325793719439782D+00, &
     -0.11913325793719439782D+00, &
      0.74545838141881970440D+00, &
     -0.74545838141881970440D+00, &
      0.74545838141881970440D+00, &
     -0.74545838141881970440D+00, &
      0.76104952097303479874D+00, &
     -0.76104952097303479874D+00, &
      0.76104952097303479874D+00, &
     -0.76104952097303479874D+00, &
      0.24208120588771878112D+00, &
     -0.24208120588771878112D+00, &
      0.24208120588771878112D+00, &
     -0.24208120588771878112D+00, &
      0.84457983815668913330D+00, &
     -0.84457983815668913330D+00, &
      0.84457983815668913330D+00, &
     -0.84457983815668913330D+00, &
      0.16456465781662479864D+00, &
     -0.16456465781662479864D+00, &
      0.16456465781662479864D+00, &
     -0.16456465781662479864D+00, &
      0.27478423472388602278D+00, &
     -0.27478423472388602278D+00, &
      0.27478423472388602278D+00, &
     -0.27478423472388602278D+00, &
      0.58499194505600049521D+00, &
     -0.58499194505600049521D+00, &
      0.58499194505600049521D+00, &
     -0.58499194505600049521D+00, &
      0.14897381246033866709D+00, &
     -0.14897381246033866709D+00, &
      0.14897381246033866709D+00, &
     -0.14897381246033866709D+00, &
      0.43892903976173064384D+00, &
     -0.43892903976173064384D+00, &
      0.43892903976173064384D+00, &
     -0.43892903976173064384D+00, &
      0.87981971200037256686D+00, &
     -0.87981971200037256686D+00, &
      0.87981971200037256686D+00, &
     -0.87981971200037256686D+00, &
      0.11685331930792351718D+00, &
     -0.11685331930792351718D+00, &
      0.11685331930792351718D+00, &
     -0.11685331930792351718D+00, &
      0.34287164858560836844D+00, &
     -0.34287164858560836844D+00, &
      0.34287164858560836844D+00, &
     -0.34287164858560836844D+00, &
      0.86433279841860011228D+00, &
      0.64222367435367522237D+00, &
     -0.86433279841860011228D+00, &
      0.64222367435367522237D+00, &
      0.86433279841860011228D+00, &
     -0.64222367435367522237D+00, &
     -0.86433279841860011228D+00, &
     -0.64222367435367522237D+00, &
      0.57840479596659710726D+00, &
      0.15700371479211361336D+00, &
     -0.57840479596659710726D+00, &
      0.15700371479211361336D+00, &
      0.57840479596659710726D+00, &
     -0.15700371479211361336D+00, &
     -0.57840479596659710726D+00, &
     -0.15700371479211361336D+00, &
      0.31087991149352928177D+00, &
      0.53916009809927645247D+00, &
     -0.31087991149352928177D+00, &
      0.53916009809927645247D+00, &
      0.31087991149352928177D+00, &
     -0.53916009809927645247D+00, &
     -0.31087991149352928177D+00, &
     -0.53916009809927645247D+00, &
      0.20472565121250529963D+00, &
      0.56505069862686030380D+00, &
     -0.20472565121250529963D+00, &
      0.56505069862686030380D+00, &
      0.20472565121250529963D+00, &
     -0.56505069862686030380D+00, &
     -0.20472565121250529963D+00, &
     -0.56505069862686030380D+00, &
      0.15477148314547456431D+00, &
      0.80087300859132681818D+00, &
     -0.15477148314547456431D+00, &
      0.80087300859132681818D+00, &
      0.15477148314547456431D+00, &
     -0.80087300859132681818D+00, &
     -0.15477148314547456431D+00, &
     -0.80087300859132681818D+00, &
      0.91122224850183075606D+00, &
      0.36196927102445491942D+00, &
     -0.91122224850183075606D+00, &
      0.36196927102445491942D+00, &
      0.91122224850183075606D+00, &
     -0.36196927102445491942D+00, &
     -0.91122224850183075606D+00, &
     -0.36196927102445491942D+00, &
      0.23763934442358639054D+00, &
      0.46054474775257370212D+00, &
     -0.23763934442358639054D+00, &
      0.46054474775257370212D+00, &
      0.23763934442358639054D+00, &
     -0.46054474775257370212D+00, &
     -0.23763934442358639054D+00, &
     -0.46054474775257370212D+00, &
      0.93026675012099768747D+00, &
      0.32159807642681254025D+00, &
     -0.93026675012099768747D+00, &
      0.32159807642681254025D+00, &
      0.93026675012099768747D+00, &
     -0.32159807642681254025D+00, &
     -0.93026675012099768747D+00, &
     -0.32159807642681254025D+00, &
      0.27330670991417554960D+00, &
      0.75417578802823215245D+00, &
     -0.27330670991417554960D+00, &
      0.75417578802823215245D+00, &
      0.27330670991417554960D+00, &
     -0.75417578802823215245D+00, &
     -0.27330670991417554960D+00, &
     -0.75417578802823215245D+00, &
      0.61977908311233775862D+00, &
      0.40094999523846291956D+00, &
     -0.61977908311233775862D+00, &
      0.40094999523846291956D+00, &
      0.61977908311233775862D+00, &
     -0.40094999523846291956D+00, &
     -0.61977908311233775862D+00, &
     -0.40094999523846291956D+00, &
      0.96731797310728928618D+00, &
      0.71750306782387462956D+00, &
     -0.96731797310728928618D+00, &
      0.71750306782387462956D+00, &
      0.96731797310728928618D+00, &
     -0.71750306782387462956D+00, &
     -0.96731797310728928618D+00, &
     -0.71750306782387462956D+00, &
      0.28107596241726462427D+00, &
      0.14469640753388121612D+00, &
     -0.28107596241726462427D+00, &
      0.14469640753388121612D+00, &
      0.28107596241726462427D+00, &
     -0.14469640753388121612D+00, &
     -0.28107596241726462427D+00, &
     -0.14469640753388121612D+00, &
      0.95723494789717944453D+00, &
      0.80544313431837344375D+00, &
     -0.95723494789717944453D+00, &
      0.80544313431837344375D+00, &
      0.95723494789717944453D+00, &
     -0.80544313431837344375D+00, &
     -0.95723494789717944453D+00, &
     -0.80544313431837344375D+00, &
      0.93415043892659577196D+00, &
      0.56426321835446102693D+00, &
     -0.93415043892659577196D+00, &
      0.56426321835446102693D+00, &
      0.93415043892659577196D+00, &
     -0.56426321835446102693D+00, &
     -0.93415043892659577196D+00, &
     -0.56426321835446102693D+00, &
      0.22124438393817494330D+00, &
      0.27637233377156195102D+00, &
     -0.22124438393817494330D+00, &
      0.27637233377156195102D+00, &
      0.22124438393817494330D+00, &
     -0.27637233377156195102D+00, &
     -0.22124438393817494330D+00, &
     -0.27637233377156195102D+00, &
      0.51440860154860901243D+00, &
      0.67110969363741279636D+00, &
     -0.51440860154860901243D+00, &
      0.67110969363741279636D+00, &
      0.51440860154860901243D+00, &
     -0.67110969363741279636D+00, &
     -0.51440860154860901243D+00, &
     -0.67110969363741279636D+00, &
      0.84669681763194848401D+00, &
      0.64763095243614277052D+00, &
     -0.84669681763194848401D+00, &
      0.64763095243614277052D+00, &
      0.84669681763194848401D+00, &
     -0.64763095243614277052D+00, &
     -0.84669681763194848401D+00, &
     -0.64763095243614277052D+00, &
      0.58372340157348212575D+00, &
      0.75878491584810681125D+00, &
     -0.58372340157348212575D+00, &
      0.75878491584810681125D+00, &
      0.58372340157348212575D+00, &
     -0.75878491584810681125D+00, &
     -0.58372340157348212575D+00, &
     -0.75878491584810681125D+00, &
      0.23125484278203506383D+00, &
      0.75795786639659157302D+00, &
     -0.23125484278203506383D+00, &
      0.75795786639659157302D+00, &
      0.23125484278203506383D+00, &
     -0.75795786639659157302D+00, &
     -0.23125484278203506383D+00, &
     -0.75795786639659157302D+00, &
      0.65351981646876278198D+00, &
      0.21186490812432900999D+00, &
     -0.65351981646876278198D+00, &
      0.21186490812432900999D+00, &
      0.65351981646876278198D+00, &
     -0.21186490812432900999D+00, &
     -0.65351981646876278198D+00, &
     -0.21186490812432900999D+00, &
      0.69534501524896896729D+00, &
      0.18605131637198926708D+00, &
     -0.69534501524896896729D+00, &
      0.18605131637198926708D+00, &
      0.69534501524896896729D+00, &
     -0.18605131637198926708D+00, &
     -0.69534501524896896729D+00, &
     -0.18605131637198926708D+00, &
      0.08127234831798690884D+00, &
      0.52410330742380673019D+00, &
     -0.08127234831798690884D+00, &
      0.52410330742380673019D+00, &
      0.08127234831798690884D+00, &
     -0.52410330742380673019D+00, &
     -0.08127234831798690884D+00, &
     -0.52410330742380673019D+00, &
      0.67316244188134488624D+00, &
      0.47074593361034250405D+00, &
     -0.67316244188134488624D+00, &
      0.47074593361034250405D+00, &
      0.67316244188134488624D+00, &
     -0.47074593361034250405D+00, &
     -0.67316244188134488624D+00, &
     -0.47074593361034250405D+00, &
      0.82977152170712553669D+00, &
      0.38456345752643983360D+00, &
     -0.82977152170712553669D+00, &
      0.38456345752643983360D+00, &
      0.82977152170712553669D+00, &
     -0.38456345752643983360D+00, &
     -0.82977152170712553669D+00, &
     -0.38456345752643983360D+00, &
      0.85308428733291663537D+00, &
      0.48551299054356161777D+00, &
     -0.85308428733291663537D+00, &
      0.48551299054356161777D+00, &
      0.85308428733291663537D+00, &
     -0.48551299054356161777D+00, &
     -0.85308428733291663537D+00, &
     -0.48551299054356161777D+00, &
      0.44670731850691264286D+00, &
      0.17757095991022572856D+00, &
     -0.44670731850691264286D+00, &
      0.17757095991022572856D+00, &
      0.44670731850691264286D+00, &
     -0.17757095991022572856D+00, &
     -0.44670731850691264286D+00, &
     -0.17757095991022572856D+00, &
      0.67009672120301455589D+00, &
      0.46199941611318601220D+00, &
     -0.67009672120301455589D+00, &
      0.46199941611318601220D+00, &
      0.67009672120301455589D+00, &
     -0.46199941611318601220D+00, &
     -0.67009672120301455589D+00, &
     -0.46199941611318601220D+00, &
      0.97673360371074724462D+00, &
      0.26790019276269744219D+00, &
     -0.97673360371074724462D+00, &
      0.26790019276269744219D+00, &
      0.97673360371074724462D+00, &
     -0.26790019276269744219D+00, &
     -0.97673360371074724462D+00, &
     -0.26790019276269744219D+00, &
      0.31887900072125263673D+00, &
      0.45381764631367643714D+00, &
     -0.31887900072125263673D+00, &
      0.45381764631367643714D+00, &
      0.31887900072125263673D+00, &
     -0.45381764631367643714D+00, &
     -0.31887900072125263673D+00, &
     -0.45381764631367643714D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.37365721634528792361D+00, &
      0.00000000000000000000D+00, &
     -0.37365721634528792361D+00, &
      0.00000000000000000000D+00, &
      0.43752827030428453892D+00, &
      0.00000000000000000000D+00, &
     -0.43752827030428453892D+00, &
      0.00000000000000000000D+00, &
      0.55782756293242685697D+00, &
      0.00000000000000000000D+00, &
     -0.55782756293242685697D+00, &
      0.00000000000000000000D+00, &
      0.27275879229149846417D+00, &
      0.00000000000000000000D+00, &
     -0.27275879229149846417D+00, &
      0.00000000000000000000D+00, &
      0.28872411753825932834D+00, &
      0.00000000000000000000D+00, &
     -0.28872411753825932834D+00, &
      0.00000000000000000000D+00, &
      0.11277465623250763904D+00, &
      0.00000000000000000000D+00, &
     -0.11277465623250763904D+00, &
      0.00000000000000000000D+00, &
      0.65470230415224206499D+00, &
      0.00000000000000000000D+00, &
     -0.65470230415224206499D+00, &
      0.00000000000000000000D+00, &
      0.22295101282662088682D+00, &
      0.00000000000000000000D+00, &
     -0.22295101282662088682D+00, &
      0.00000000000000000000D+00, &
      0.35176689121153298379D+00, &
      0.00000000000000000000D+00, &
     -0.35176689121153298379D+00, &
      0.00000000000000000000D+00, &
      0.20697735413665671600D+00, &
      0.00000000000000000000D+00, &
     -0.20697735413665671600D+00, &
      0.00000000000000000000D+00, &
      0.93111189027726670719D+00, &
      0.00000000000000000000D+00, &
     -0.93111189027726670719D+00, &
      0.00000000000000000000D+00, &
      0.84303358445951348532D+00, &
      0.00000000000000000000D+00, &
     -0.84303358445951348532D+00, &
      0.00000000000000000000D+00, &
      0.07303948013428374075D+00, &
      0.00000000000000000000D+00, &
     -0.07303948013428374075D+00, &
      0.00000000000000000000D+00, &
      0.67807351338628540915D+00, &
      0.00000000000000000000D+00, &
     -0.67807351338628540915D+00, &
      0.00000000000000000000D+00, &
      0.88198482628664209848D+00, &
      0.00000000000000000000D+00, &
     -0.88198482628664209848D+00, &
      0.00000000000000000000D+00, &
      0.38587811686103651310D+00, &
      0.00000000000000000000D+00, &
     -0.38587811686103651310D+00, &
      0.00000000000000000000D+00, &
      0.13708884737518356456D+00, &
      0.00000000000000000000D+00, &
     -0.13708884737518356456D+00, &
      0.52558521108268729805D+00, &
      0.52558521108268729805D+00, &
     -0.52558521108268729805D+00, &
     -0.52558521108268729805D+00, &
      0.45301642570589800707D+00, &
      0.45301642570589800707D+00, &
     -0.45301642570589800707D+00, &
     -0.45301642570589800707D+00, &
      0.67125116320175959306D+00, &
      0.67125116320175959306D+00, &
     -0.67125116320175959306D+00, &
     -0.67125116320175959306D+00, &
      0.36951418649168821240D+00, &
      0.36951418649168821240D+00, &
     -0.36951418649168821240D+00, &
     -0.36951418649168821240D+00, &
      0.28838223435586252119D+00, &
      0.28838223435586252119D+00, &
     -0.28838223435586252119D+00, &
     -0.28838223435586252119D+00, &
      0.82857238389465803774D+00, &
      0.82857238389465803774D+00, &
     -0.82857238389465803774D+00, &
     -0.82857238389465803774D+00, &
      0.67246070085447195996D+00, &
      0.67246070085447195996D+00, &
     -0.67246070085447195996D+00, &
     -0.67246070085447195996D+00, &
      0.48736689338962557727D+00, &
      0.48736689338962557727D+00, &
     -0.48736689338962557727D+00, &
     -0.48736689338962557727D+00, &
      0.95036644315079843448D+00, &
      0.95036644315079843448D+00, &
     -0.95036644315079843448D+00, &
     -0.95036644315079843448D+00, &
      0.09132353058612172059D+00, &
      0.09132353058612172059D+00, &
     -0.09132353058612172059D+00, &
     -0.09132353058612172059D+00, &
      0.15520695551416591185D+00, &
      0.15520695551416591185D+00, &
     -0.15520695551416591185D+00, &
     -0.15520695551416591185D+00, &
      0.05052806321356113906D+00, &
      0.05052806321356113906D+00, &
     -0.05052806321356113906D+00, &
     -0.05052806321356113906D+00, &
      0.39781392804166187949D+00, &
      0.39781392804166187949D+00, &
     -0.39781392804166187949D+00, &
     -0.39781392804166187949D+00, &
      0.37607913140641990868D+00, &
      0.37607913140641990868D+00, &
     -0.37607913140641990868D+00, &
     -0.37607913140641990868D+00, &
      0.16436382679828262510D+00, &
      0.16436382679828262510D+00, &
     -0.16436382679828262510D+00, &
     -0.16436382679828262510D+00, &
      0.34688925556861166521D+00, &
      0.34688925556861166521D+00, &
     -0.34688925556861166521D+00, &
     -0.34688925556861166521D+00, &
      0.11913325793719439782D+00, &
      0.11913325793719439782D+00, &
     -0.11913325793719439782D+00, &
     -0.11913325793719439782D+00, &
      0.74545838141881970440D+00, &
      0.74545838141881970440D+00, &
     -0.74545838141881970440D+00, &
     -0.74545838141881970440D+00, &
      0.76104952097303479874D+00, &
      0.76104952097303479874D+00, &
     -0.76104952097303479874D+00, &
     -0.76104952097303479874D+00, &
      0.24208120588771878112D+00, &
      0.24208120588771878112D+00, &
     -0.24208120588771878112D+00, &
     -0.24208120588771878112D+00, &
      0.84457983815668913330D+00, &
      0.84457983815668913330D+00, &
     -0.84457983815668913330D+00, &
     -0.84457983815668913330D+00, &
      0.16456465781662479864D+00, &
      0.16456465781662479864D+00, &
     -0.16456465781662479864D+00, &
     -0.16456465781662479864D+00, &
      0.27478423472388602278D+00, &
      0.27478423472388602278D+00, &
     -0.27478423472388602278D+00, &
     -0.27478423472388602278D+00, &
      0.58499194505600049521D+00, &
      0.58499194505600049521D+00, &
     -0.58499194505600049521D+00, &
     -0.58499194505600049521D+00, &
      0.14897381246033866709D+00, &
      0.14897381246033866709D+00, &
     -0.14897381246033866709D+00, &
     -0.14897381246033866709D+00, &
      0.43892903976173064384D+00, &
      0.43892903976173064384D+00, &
     -0.43892903976173064384D+00, &
     -0.43892903976173064384D+00, &
      0.87981971200037256686D+00, &
      0.87981971200037256686D+00, &
     -0.87981971200037256686D+00, &
     -0.87981971200037256686D+00, &
      0.11685331930792351718D+00, &
      0.11685331930792351718D+00, &
     -0.11685331930792351718D+00, &
     -0.11685331930792351718D+00, &
      0.34287164858560836844D+00, &
      0.34287164858560836844D+00, &
     -0.34287164858560836844D+00, &
     -0.34287164858560836844D+00, &
      0.64222367435367522237D+00, &
      0.86433279841860011228D+00, &
      0.64222367435367522237D+00, &
     -0.86433279841860011228D+00, &
     -0.64222367435367522237D+00, &
      0.86433279841860011228D+00, &
     -0.64222367435367522237D+00, &
     -0.86433279841860011228D+00, &
      0.15700371479211361336D+00, &
      0.57840479596659710726D+00, &
      0.15700371479211361336D+00, &
     -0.57840479596659710726D+00, &
     -0.15700371479211361336D+00, &
      0.57840479596659710726D+00, &
     -0.15700371479211361336D+00, &
     -0.57840479596659710726D+00, &
      0.53916009809927645247D+00, &
      0.31087991149352928177D+00, &
      0.53916009809927645247D+00, &
     -0.31087991149352928177D+00, &
     -0.53916009809927645247D+00, &
      0.31087991149352928177D+00, &
     -0.53916009809927645247D+00, &
     -0.31087991149352928177D+00, &
      0.56505069862686030380D+00, &
      0.20472565121250529963D+00, &
      0.56505069862686030380D+00, &
     -0.20472565121250529963D+00, &
     -0.56505069862686030380D+00, &
      0.20472565121250529963D+00, &
     -0.56505069862686030380D+00, &
     -0.20472565121250529963D+00, &
      0.80087300859132681818D+00, &
      0.15477148314547456431D+00, &
      0.80087300859132681818D+00, &
     -0.15477148314547456431D+00, &
     -0.80087300859132681818D+00, &
      0.15477148314547456431D+00, &
     -0.80087300859132681818D+00, &
     -0.15477148314547456431D+00, &
      0.36196927102445491942D+00, &
      0.91122224850183075606D+00, &
      0.36196927102445491942D+00, &
     -0.91122224850183075606D+00, &
     -0.36196927102445491942D+00, &
      0.91122224850183075606D+00, &
     -0.36196927102445491942D+00, &
     -0.91122224850183075606D+00, &
      0.46054474775257370212D+00, &
      0.23763934442358639054D+00, &
      0.46054474775257370212D+00, &
     -0.23763934442358639054D+00, &
     -0.46054474775257370212D+00, &
      0.23763934442358639054D+00, &
     -0.46054474775257370212D+00, &
     -0.23763934442358639054D+00, &
      0.32159807642681254025D+00, &
      0.93026675012099768747D+00, &
      0.32159807642681254025D+00, &
     -0.93026675012099768747D+00, &
     -0.32159807642681254025D+00, &
      0.93026675012099768747D+00, &
     -0.32159807642681254025D+00, &
     -0.93026675012099768747D+00, &
      0.75417578802823215245D+00, &
      0.27330670991417554960D+00, &
      0.75417578802823215245D+00, &
     -0.27330670991417554960D+00, &
     -0.75417578802823215245D+00, &
      0.27330670991417554960D+00, &
     -0.75417578802823215245D+00, &
     -0.27330670991417554960D+00, &
      0.40094999523846291956D+00, &
      0.61977908311233775862D+00, &
      0.40094999523846291956D+00, &
     -0.61977908311233775862D+00, &
     -0.40094999523846291956D+00, &
      0.61977908311233775862D+00, &
     -0.40094999523846291956D+00, &
     -0.61977908311233775862D+00, &
      0.71750306782387462956D+00, &
      0.96731797310728928618D+00, &
      0.71750306782387462956D+00, &
     -0.96731797310728928618D+00, &
     -0.71750306782387462956D+00, &
      0.96731797310728928618D+00, &
     -0.71750306782387462956D+00, &
     -0.96731797310728928618D+00, &
      0.14469640753388121612D+00, &
      0.28107596241726462427D+00, &
      0.14469640753388121612D+00, &
     -0.28107596241726462427D+00, &
     -0.14469640753388121612D+00, &
      0.28107596241726462427D+00, &
     -0.14469640753388121612D+00, &
     -0.28107596241726462427D+00, &
      0.80544313431837344375D+00, &
      0.95723494789717944453D+00, &
      0.80544313431837344375D+00, &
     -0.95723494789717944453D+00, &
     -0.80544313431837344375D+00, &
      0.95723494789717944453D+00, &
     -0.80544313431837344375D+00, &
     -0.95723494789717944453D+00, &
      0.56426321835446102693D+00, &
      0.93415043892659577196D+00, &
      0.56426321835446102693D+00, &
     -0.93415043892659577196D+00, &
     -0.56426321835446102693D+00, &
      0.93415043892659577196D+00, &
     -0.56426321835446102693D+00, &
     -0.93415043892659577196D+00, &
      0.27637233377156195102D+00, &
      0.22124438393817494330D+00, &
      0.27637233377156195102D+00, &
     -0.22124438393817494330D+00, &
     -0.27637233377156195102D+00, &
      0.22124438393817494330D+00, &
     -0.27637233377156195102D+00, &
     -0.22124438393817494330D+00, &
      0.67110969363741279636D+00, &
      0.51440860154860901243D+00, &
      0.67110969363741279636D+00, &
     -0.51440860154860901243D+00, &
     -0.67110969363741279636D+00, &
      0.51440860154860901243D+00, &
     -0.67110969363741279636D+00, &
     -0.51440860154860901243D+00, &
      0.64763095243614277052D+00, &
      0.84669681763194848401D+00, &
      0.64763095243614277052D+00, &
     -0.84669681763194848401D+00, &
     -0.64763095243614277052D+00, &
      0.84669681763194848401D+00, &
     -0.64763095243614277052D+00, &
     -0.84669681763194848401D+00, &
      0.75878491584810681125D+00, &
      0.58372340157348212575D+00, &
      0.75878491584810681125D+00, &
     -0.58372340157348212575D+00, &
     -0.75878491584810681125D+00, &
      0.58372340157348212575D+00, &
     -0.75878491584810681125D+00, &
     -0.58372340157348212575D+00, &
      0.75795786639659157302D+00, &
      0.23125484278203506383D+00, &
      0.75795786639659157302D+00, &
     -0.23125484278203506383D+00, &
     -0.75795786639659157302D+00, &
      0.23125484278203506383D+00, &
     -0.75795786639659157302D+00, &
     -0.23125484278203506383D+00, &
      0.21186490812432900999D+00, &
      0.65351981646876278198D+00, &
      0.21186490812432900999D+00, &
     -0.65351981646876278198D+00, &
     -0.21186490812432900999D+00, &
      0.65351981646876278198D+00, &
     -0.21186490812432900999D+00, &
     -0.65351981646876278198D+00, &
      0.18605131637198926708D+00, &
      0.69534501524896896729D+00, &
      0.18605131637198926708D+00, &
     -0.69534501524896896729D+00, &
     -0.18605131637198926708D+00, &
      0.69534501524896896729D+00, &
     -0.18605131637198926708D+00, &
     -0.69534501524896896729D+00, &
      0.52410330742380673019D+00, &
      0.08127234831798690884D+00, &
      0.52410330742380673019D+00, &
     -0.08127234831798690884D+00, &
     -0.52410330742380673019D+00, &
      0.08127234831798690884D+00, &
     -0.52410330742380673019D+00, &
     -0.08127234831798690884D+00, &
      0.47074593361034250405D+00, &
      0.67316244188134488624D+00, &
      0.47074593361034250405D+00, &
     -0.67316244188134488624D+00, &
     -0.47074593361034250405D+00, &
      0.67316244188134488624D+00, &
     -0.47074593361034250405D+00, &
     -0.67316244188134488624D+00, &
      0.38456345752643983360D+00, &
      0.82977152170712553669D+00, &
      0.38456345752643983360D+00, &
     -0.82977152170712553669D+00, &
     -0.38456345752643983360D+00, &
      0.82977152170712553669D+00, &
     -0.38456345752643983360D+00, &
     -0.82977152170712553669D+00, &
      0.48551299054356161777D+00, &
      0.85308428733291663537D+00, &
      0.48551299054356161777D+00, &
     -0.85308428733291663537D+00, &
     -0.48551299054356161777D+00, &
      0.85308428733291663537D+00, &
     -0.48551299054356161777D+00, &
     -0.85308428733291663537D+00, &
      0.17757095991022572856D+00, &
      0.44670731850691264286D+00, &
      0.17757095991022572856D+00, &
     -0.44670731850691264286D+00, &
     -0.17757095991022572856D+00, &
      0.44670731850691264286D+00, &
     -0.17757095991022572856D+00, &
     -0.44670731850691264286D+00, &
      0.46199941611318601220D+00, &
      0.67009672120301455589D+00, &
      0.46199941611318601220D+00, &
     -0.67009672120301455589D+00, &
     -0.46199941611318601220D+00, &
      0.67009672120301455589D+00, &
     -0.46199941611318601220D+00, &
     -0.67009672120301455589D+00, &
      0.26790019276269744219D+00, &
      0.97673360371074724462D+00, &
      0.26790019276269744219D+00, &
     -0.97673360371074724462D+00, &
     -0.26790019276269744219D+00, &
      0.97673360371074724462D+00, &
     -0.26790019276269744219D+00, &
     -0.97673360371074724462D+00, &
      0.45381764631367643714D+00, &
      0.31887900072125263673D+00, &
      0.45381764631367643714D+00, &
     -0.31887900072125263673D+00, &
     -0.45381764631367643714D+00, &
      0.31887900072125263673D+00, &
     -0.45381764631367643714D+00, &
     -0.31887900072125263673D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
      0.96739133523703502160D+00, &
      0.54366532464638728239D+00, &
      0.41887560622759456574D+00, &
      0.41887560622759456574D+00, &
      0.41887560622759456574D+00, &
      0.41887560622759456574D+00, &
      0.29774894445685839983D+00, &
      0.29774894445685839983D+00, &
      0.29774894445685839983D+00, &
      0.29774894445685839983D+00, &
      0.43632558904359264318D+00, &
      0.43632558904359264318D+00, &
      0.43632558904359264318D+00, &
      0.43632558904359264318D+00, &
      0.13123475865230219140D+00, &
      0.13123475865230219140D+00, &
      0.13123475865230219140D+00, &
      0.13123475865230219140D+00, &
      0.56882169546815719574D+00, &
      0.56882169546815719574D+00, &
      0.56882169546815719574D+00, &
      0.56882169546815719574D+00, &
      0.66230573631886058283D+00, &
      0.66230573631886058283D+00, &
      0.66230573631886058283D+00, &
      0.66230573631886058283D+00, &
      0.21963789124598576130D+00, &
      0.21963789124598576130D+00, &
      0.21963789124598576130D+00, &
      0.21963789124598576130D+00, &
      0.22732696332834559372D+00, &
      0.22732696332834559372D+00, &
      0.22732696332834559372D+00, &
      0.22732696332834559372D+00, &
      0.00541751512625303269D+00, &
      0.00541751512625303269D+00, &
      0.00541751512625303269D+00, &
      0.00541751512625303269D+00, &
      0.79294906136483345183D+00, &
      0.79294906136483345183D+00, &
      0.79294906136483345183D+00, &
      0.79294906136483345183D+00, &
      0.04180662980331771583D+00, &
      0.04180662980331771583D+00, &
      0.04180662980331771583D+00, &
      0.04180662980331771583D+00, &
      0.13615518064579532065D+00, &
      0.13615518064579532065D+00, &
      0.13615518064579532065D+00, &
      0.13615518064579532065D+00, &
      0.91875475602806166986D+00, &
      0.91875475602806166986D+00, &
      0.91875475602806166986D+00, &
      0.91875475602806166986D+00, &
      0.30191905621780101843D+00, &
      0.30191905621780101843D+00, &
      0.30191905621780101843D+00, &
      0.30191905621780101843D+00, &
      0.00938749984908901144D+00, &
      0.00938749984908901144D+00, &
      0.00938749984908901144D+00, &
      0.00938749984908901144D+00, &
      0.60658598457462375997D+00, &
      0.60658598457462375997D+00, &
      0.60658598457462375997D+00, &
      0.60658598457462375997D+00, &
      0.75353482202105526166D+00, &
      0.75353482202105526166D+00, &
      0.75353482202105526166D+00, &
      0.75353482202105526166D+00, &
      0.44086853690608979184D+00, &
      0.44086853690608979184D+00, &
      0.44086853690608979184D+00, &
      0.44086853690608979184D+00, &
      0.41889775081556313019D+00, &
      0.41889775081556313019D+00, &
      0.41889775081556313019D+00, &
      0.41889775081556313019D+00, &
      0.01251810432139164979D+00, &
      0.01251810432139164979D+00, &
      0.01251810432139164979D+00, &
      0.01251810432139164979D+00, &
      0.62460148383811031625D+00, &
      0.62460148383811031625D+00, &
      0.62460148383811031625D+00, &
      0.62460148383811031625D+00, &
      0.05439973273857341696D+00, &
      0.05439973273857341696D+00, &
      0.05439973273857341696D+00, &
      0.05439973273857341696D+00, &
      0.16191862193166678408D+00, &
      0.16191862193166678408D+00, &
      0.16191862193166678408D+00, &
      0.16191862193166678408D+00, &
      0.29300331383428240839D+00, &
      0.29300331383428240839D+00, &
      0.29300331383428240839D+00, &
      0.29300331383428240839D+00, &
      0.26342658893870185555D+00, &
      0.26342658893870185555D+00, &
      0.26342658893870185555D+00, &
      0.26342658893870185555D+00, &
      0.01127614473867753070D+00, &
      0.01127614473867753070D+00, &
      0.01127614473867753070D+00, &
      0.01127614473867753070D+00, &
      0.35807984614999605055D+00, &
      0.35807984614999605055D+00, &
      0.35807984614999605055D+00, &
      0.35807984614999605055D+00, &
      0.47141886716708586436D+00, &
      0.47141886716708586436D+00, &
      0.47141886716708586436D+00, &
      0.47141886716708586436D+00, &
      0.85446407112095401626D+00, &
      0.85446407112095401626D+00, &
      0.85446407112095401626D+00, &
      0.85446407112095401626D+00, &
      0.12930311691456014556D+00, &
      0.12930311691456014556D+00, &
      0.12930311691456014556D+00, &
      0.12930311691456014556D+00, &
      0.07899233073306564934D+00, &
      0.07899233073306564934D+00, &
      0.07899233073306564934D+00, &
      0.07899233073306564934D+00, &
      0.01672886065733771022D+00, &
      0.01672886065733771022D+00, &
      0.01672886065733771022D+00, &
      0.01672886065733771022D+00, &
      0.38850504372376110096D+00, &
      0.38850504372376110096D+00, &
      0.38850504372376110096D+00, &
      0.38850504372376110096D+00, &
      0.84754739757493258168D+00, &
      0.84754739757493258168D+00, &
      0.84754739757493258168D+00, &
      0.84754739757493258168D+00, &
      0.16263602800181631292D+00, &
      0.16263602800181631292D+00, &
      0.16263602800181631292D+00, &
      0.16263602800181631292D+00, &
      0.07862268330374026781D+00, &
      0.07862268330374026781D+00, &
      0.07862268330374026781D+00, &
      0.07862268330374026781D+00, &
      0.32028064376194109730D+00, &
      0.32028064376194109730D+00, &
      0.32028064376194109730D+00, &
      0.32028064376194109730D+00, &
      0.01382133423647079022D+00, &
      0.01382133423647079022D+00, &
      0.01382133423647079022D+00, &
      0.01382133423647079022D+00, &
      0.76877522504357953537D+00, &
      0.76877522504357953537D+00, &
      0.76877522504357953537D+00, &
      0.76877522504357953537D+00, &
      0.53014024087305255950D+00, &
      0.53014024087305255950D+00, &
      0.53014024087305255950D+00, &
      0.53014024087305255950D+00, &
      0.30089238244412036538D+00, &
      0.30089238244412036538D+00, &
      0.30089238244412036538D+00, &
      0.30089238244412036538D+00, &
      0.61673192700785928189D+00, &
      0.61673192700785928189D+00, &
      0.61673192700785928189D+00, &
      0.61673192700785928189D+00, &
      0.01376377254881575897D+00, &
      0.01376377254881575897D+00, &
      0.01376377254881575897D+00, &
      0.01376377254881575897D+00, &
      0.07210305012229555055D+00, &
      0.07210305012229555055D+00, &
      0.07210305012229555055D+00, &
      0.07210305012229555055D+00, &
      0.06238118616203952582D+00, &
      0.06238118616203952582D+00, &
      0.06238118616203952582D+00, &
      0.06238118616203952582D+00, &
      0.58557953491658876199D+00, &
      0.58557953491658876199D+00, &
      0.58557953491658876199D+00, &
      0.58557953491658876199D+00, &
      0.11791920647798490029D+00, &
      0.11791920647798490029D+00, &
      0.11791920647798490029D+00, &
      0.11791920647798490029D+00, &
      0.11791920647798490029D+00, &
      0.11791920647798490029D+00, &
      0.11791920647798490029D+00, &
      0.11791920647798490029D+00, &
      0.11142523412928505289D+00, &
      0.11142523412928505289D+00, &
      0.11142523412928505289D+00, &
      0.11142523412928505289D+00, &
      0.11142523412928505289D+00, &
      0.11142523412928505289D+00, &
      0.11142523412928505289D+00, &
      0.11142523412928505289D+00, &
      0.27502672812187933804D+00, &
      0.27502672812187933804D+00, &
      0.27502672812187933804D+00, &
      0.27502672812187933804D+00, &
      0.27502672812187933804D+00, &
      0.27502672812187933804D+00, &
      0.27502672812187933804D+00, &
      0.27502672812187933804D+00, &
      0.34728469726559213493D+00, &
      0.34728469726559213493D+00, &
      0.34728469726559213493D+00, &
      0.34728469726559213493D+00, &
      0.34728469726559213493D+00, &
      0.34728469726559213493D+00, &
      0.34728469726559213493D+00, &
      0.34728469726559213493D+00, &
      0.07549551730236130076D+00, &
      0.07549551730236130076D+00, &
      0.07549551730236130076D+00, &
      0.07549551730236130076D+00, &
      0.07549551730236130076D+00, &
      0.07549551730236130076D+00, &
      0.07549551730236130076D+00, &
      0.07549551730236130076D+00, &
      0.03464296962031117311D+00, &
      0.03464296962031117311D+00, &
      0.03464296962031117311D+00, &
      0.03464296962031117311D+00, &
      0.03464296962031117311D+00, &
      0.03464296962031117311D+00, &
      0.03464296962031117311D+00, &
      0.03464296962031117311D+00, &
      0.20441529779235670383D+00, &
      0.20441529779235670383D+00, &
      0.20441529779235670383D+00, &
      0.20441529779235670383D+00, &
      0.20441529779235670383D+00, &
      0.20441529779235670383D+00, &
      0.20441529779235670383D+00, &
      0.20441529779235670383D+00, &
      0.06972695561067981940D+00, &
      0.06972695561067981940D+00, &
      0.06972695561067981940D+00, &
      0.06972695561067981940D+00, &
      0.06972695561067981940D+00, &
      0.06972695561067981940D+00, &
      0.06972695561067981940D+00, &
      0.06972695561067981940D+00, &
      0.21998980465815576313D+00, &
      0.21998980465815576313D+00, &
      0.21998980465815576313D+00, &
      0.21998980465815576313D+00, &
      0.21998980465815576313D+00, &
      0.21998980465815576313D+00, &
      0.21998980465815576313D+00, &
      0.21998980465815576313D+00, &
      0.36462254947959177320D+00, &
      0.36462254947959177320D+00, &
      0.36462254947959177320D+00, &
      0.36462254947959177320D+00, &
      0.36462254947959177320D+00, &
      0.36462254947959177320D+00, &
      0.36462254947959177320D+00, &
      0.36462254947959177320D+00, &
      0.00458633496122378467D+00, &
      0.00458633496122378467D+00, &
      0.00458633496122378467D+00, &
      0.00458633496122378467D+00, &
      0.00458633496122378467D+00, &
      0.00458633496122378467D+00, &
      0.00458633496122378467D+00, &
      0.00458633496122378467D+00, &
      0.66642300488713057671D+00, &
      0.66642300488713057671D+00, &
      0.66642300488713057671D+00, &
      0.66642300488713057671D+00, &
      0.66642300488713057671D+00, &
      0.66642300488713057671D+00, &
      0.66642300488713057671D+00, &
      0.66642300488713057671D+00, &
      0.03975302308577006311D+00, &
      0.03975302308577006311D+00, &
      0.03975302308577006311D+00, &
      0.03975302308577006311D+00, &
      0.03975302308577006311D+00, &
      0.03975302308577006311D+00, &
      0.03975302308577006311D+00, &
      0.03975302308577006311D+00, &
      0.03073598825835915579D+00, &
      0.03073598825835915579D+00, &
      0.03073598825835915579D+00, &
      0.03073598825835915579D+00, &
      0.03073598825835915579D+00, &
      0.03073598825835915579D+00, &
      0.03073598825835915579D+00, &
      0.03073598825835915579D+00, &
      0.72077351789673316240D+00, &
      0.72077351789673316240D+00, &
      0.72077351789673316240D+00, &
      0.72077351789673316240D+00, &
      0.72077351789673316240D+00, &
      0.72077351789673316240D+00, &
      0.72077351789673316240D+00, &
      0.72077351789673316240D+00, &
      0.14036094696120685055D+00, &
      0.14036094696120685055D+00, &
      0.14036094696120685055D+00, &
      0.14036094696120685055D+00, &
      0.14036094696120685055D+00, &
      0.14036094696120685055D+00, &
      0.14036094696120685055D+00, &
      0.14036094696120685055D+00, &
      0.05264191965632582237D+00, &
      0.05264191965632582237D+00, &
      0.05264191965632582237D+00, &
      0.05264191965632582237D+00, &
      0.05264191965632582237D+00, &
      0.05264191965632582237D+00, &
      0.05264191965632582237D+00, &
      0.05264191965632582237D+00, &
      0.22702332489421667150D+00, &
      0.22702332489421667150D+00, &
      0.22702332489421667150D+00, &
      0.22702332489421667150D+00, &
      0.22702332489421667150D+00, &
      0.22702332489421667150D+00, &
      0.22702332489421667150D+00, &
      0.22702332489421667150D+00, &
      0.02881326163838642679D+00, &
      0.02881326163838642679D+00, &
      0.02881326163838642679D+00, &
      0.02881326163838642679D+00, &
      0.02881326163838642679D+00, &
      0.02881326163838642679D+00, &
      0.02881326163838642679D+00, &
      0.02881326163838642679D+00, &
      0.00532138936818496097D+00, &
      0.00532138936818496097D+00, &
      0.00532138936818496097D+00, &
      0.00532138936818496097D+00, &
      0.00532138936818496097D+00, &
      0.00532138936818496097D+00, &
      0.00532138936818496097D+00, &
      0.00532138936818496097D+00, &
      0.15917093130673495849D+00, &
      0.15917093130673495849D+00, &
      0.15917093130673495849D+00, &
      0.15917093130673495849D+00, &
      0.15917093130673495849D+00, &
      0.15917093130673495849D+00, &
      0.15917093130673495849D+00, &
      0.15917093130673495849D+00, &
      0.04286780869626788393D+00, &
      0.04286780869626788393D+00, &
      0.04286780869626788393D+00, &
      0.04286780869626788393D+00, &
      0.04286780869626788393D+00, &
      0.04286780869626788393D+00, &
      0.04286780869626788393D+00, &
      0.04286780869626788393D+00, &
      0.22819215603891684907D+00, &
      0.22819215603891684907D+00, &
      0.22819215603891684907D+00, &
      0.22819215603891684907D+00, &
      0.22819215603891684907D+00, &
      0.22819215603891684907D+00, &
      0.22819215603891684907D+00, &
      0.22819215603891684907D+00, &
      0.11497584756253050042D+00, &
      0.11497584756253050042D+00, &
      0.11497584756253050042D+00, &
      0.11497584756253050042D+00, &
      0.11497584756253050042D+00, &
      0.11497584756253050042D+00, &
      0.11497584756253050042D+00, &
      0.11497584756253050042D+00, &
      0.00654660357363916808D+00, &
      0.00654660357363916808D+00, &
      0.00654660357363916808D+00, &
      0.00654660357363916808D+00, &
      0.00654660357363916808D+00, &
      0.00654660357363916808D+00, &
      0.00654660357363916808D+00, &
      0.00654660357363916808D+00, &
      0.48896580613114126734D+00, &
      0.48896580613114126734D+00, &
      0.48896580613114126734D+00, &
      0.48896580613114126734D+00, &
      0.48896580613114126734D+00, &
      0.48896580613114126734D+00, &
      0.48896580613114126734D+00, &
      0.48896580613114126734D+00, &
      0.05843691849737928795D+00, &
      0.05843691849737928795D+00, &
      0.05843691849737928795D+00, &
      0.05843691849737928795D+00, &
      0.05843691849737928795D+00, &
      0.05843691849737928795D+00, &
      0.05843691849737928795D+00, &
      0.05843691849737928795D+00, &
      0.00610498541327420680D+00, &
      0.00610498541327420680D+00, &
      0.00610498541327420680D+00, &
      0.00610498541327420680D+00, &
      0.00610498541327420680D+00, &
      0.00610498541327420680D+00, &
      0.00610498541327420680D+00, &
      0.00610498541327420680D+00, &
      0.53408664015576723383D+00, &
      0.53408664015576723383D+00, &
      0.53408664015576723383D+00, &
      0.53408664015576723383D+00, &
      0.53408664015576723383D+00, &
      0.53408664015576723383D+00, &
      0.53408664015576723383D+00, &
      0.53408664015576723383D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
      0.00012870022397472401D+00, &
      0.00539076125430030910D+00, &
      0.00666687105780559754D+00, &
      0.00666687105780559754D+00, &
      0.00666687105780559754D+00, &
      0.00666687105780559754D+00, &
      0.00619801121379879964D+00, &
      0.00619801121379879964D+00, &
      0.00619801121379879964D+00, &
      0.00619801121379879964D+00, &
      0.00119690394217259041D+00, &
      0.00119690394217259041D+00, &
      0.00119690394217259041D+00, &
      0.00119690394217259041D+00, &
      0.00900035197624313339D+00, &
      0.00900035197624313339D+00, &
      0.00900035197624313339D+00, &
      0.00900035197624313339D+00, &
      0.00416463710538401977D+00, &
      0.00416463710538401977D+00, &
      0.00416463710538401977D+00, &
      0.00416463710538401977D+00, &
      0.00182449720872727593D+00, &
      0.00182449720872727593D+00, &
      0.00182449720872727593D+00, &
      0.00182449720872727593D+00, &
      0.00452371868333050885D+00, &
      0.00452371868333050885D+00, &
      0.00452371868333050885D+00, &
      0.00452371868333050885D+00, &
      0.00816111606481412406D+00, &
      0.00816111606481412406D+00, &
      0.00816111606481412406D+00, &
      0.00816111606481412406D+00, &
      0.00164832682863707940D+00, &
      0.00164832682863707940D+00, &
      0.00164832682863707940D+00, &
      0.00164832682863707940D+00, &
      0.00055714923095003237D+00, &
      0.00055714923095003237D+00, &
      0.00055714923095003237D+00, &
      0.00055714923095003237D+00, &
      0.00111821525219482697D+00, &
      0.00111821525219482697D+00, &
      0.00111821525219482697D+00, &
      0.00111821525219482697D+00, &
      0.00188581015876115280D+00, &
      0.00188581015876115280D+00, &
      0.00188581015876115280D+00, &
      0.00188581015876115280D+00, &
      0.00028281012262636759D+00, &
      0.00028281012262636759D+00, &
      0.00028281012262636759D+00, &
      0.00028281012262636759D+00, &
      0.00181586354615248828D+00, &
      0.00181586354615248828D+00, &
      0.00181586354615248828D+00, &
      0.00181586354615248828D+00, &
      0.00124708244822032614D+00, &
      0.00124708244822032614D+00, &
      0.00124708244822032614D+00, &
      0.00124708244822032614D+00, &
      0.00107416255302727809D+00, &
      0.00107416255302727809D+00, &
      0.00107416255302727809D+00, &
      0.00107416255302727809D+00, &
      0.00257977077737831396D+00, &
      0.00257977077737831396D+00, &
      0.00257977077737831396D+00, &
      0.00257977077737831396D+00, &
      0.00091701544972603221D+00, &
      0.00091701544972603221D+00, &
      0.00091701544972603221D+00, &
      0.00091701544972603221D+00, &
      0.00306877872117356605D+00, &
      0.00306877872117356605D+00, &
      0.00306877872117356605D+00, &
      0.00306877872117356605D+00, &
      0.00187105817191728823D+00, &
      0.00187105817191728823D+00, &
      0.00187105817191728823D+00, &
      0.00187105817191728823D+00, &
      0.00015369093908062739D+00, &
      0.00015369093908062739D+00, &
      0.00015369093908062739D+00, &
      0.00015369093908062739D+00, &
      0.00384529073211040779D+00, &
      0.00384529073211040779D+00, &
      0.00384529073211040779D+00, &
      0.00384529073211040779D+00, &
      0.00022431161642440681D+00, &
      0.00022431161642440681D+00, &
      0.00022431161642440681D+00, &
      0.00022431161642440681D+00, &
      0.00084020876261339610D+00, &
      0.00084020876261339610D+00, &
      0.00084020876261339610D+00, &
      0.00084020876261339610D+00, &
      0.00304195767454043140D+00, &
      0.00304195767454043140D+00, &
      0.00304195767454043140D+00, &
      0.00304195767454043140D+00, &
      0.00022948417669879317D+00, &
      0.00022948417669879317D+00, &
      0.00022948417669879317D+00, &
      0.00022948417669879317D+00, &
      0.00418242059186963843D+00, &
      0.00418242059186963843D+00, &
      0.00418242059186963843D+00, &
      0.00418242059186963843D+00, &
      0.00632709108150545339D+00, &
      0.00632709108150545339D+00, &
      0.00632709108150545339D+00, &
      0.00632709108150545339D+00, &
      0.00075928340059512709D+00, &
      0.00075928340059512709D+00, &
      0.00075928340059512709D+00, &
      0.00075928340059512709D+00, &
      0.00459683666227672476D+00, &
      0.00459683666227672476D+00, &
      0.00459683666227672476D+00, &
      0.00459683666227672476D+00, &
      0.00276636427423579932D+00, &
      0.00276636427423579932D+00, &
      0.00276636427423579932D+00, &
      0.00276636427423579932D+00, &
      0.00214485080151055087D+00, &
      0.00214485080151055087D+00, &
      0.00214485080151055087D+00, &
      0.00214485080151055087D+00, &
      0.00528314448618283795D+00, &
      0.00528314448618283795D+00, &
      0.00528314448618283795D+00, &
      0.00528314448618283795D+00, &
      0.00048810523226141780D+00, &
      0.00048810523226141780D+00, &
      0.00048810523226141780D+00, &
      0.00048810523226141780D+00, &
      0.00212069545834849948D+00, &
      0.00212069545834849948D+00, &
      0.00212069545834849948D+00, &
      0.00212069545834849948D+00, &
      0.00161506019948201812D+00, &
      0.00161506019948201812D+00, &
      0.00161506019948201812D+00, &
      0.00161506019948201812D+00, &
      0.00694225709294844437D+00, &
      0.00694225709294844437D+00, &
      0.00694225709294844437D+00, &
      0.00694225709294844437D+00, &
      0.00100268302008121369D+00, &
      0.00100268302008121369D+00, &
      0.00100268302008121369D+00, &
      0.00100268302008121369D+00, &
      0.00092077531674916973D+00, &
      0.00092077531674916973D+00, &
      0.00092077531674916973D+00, &
      0.00092077531674916973D+00, &
      0.00355616947744542573D+00, &
      0.00355616947744542573D+00, &
      0.00355616947744542573D+00, &
      0.00355616947744542573D+00, &
      0.00218617871329652470D+00, &
      0.00218617871329652470D+00, &
      0.00218617871329652470D+00, &
      0.00218617871329652470D+00, &
      0.00269165166203209818D+00, &
      0.00269165166203209818D+00, &
      0.00269165166203209818D+00, &
      0.00269165166203209818D+00, &
      0.00301674010475732127D+00, &
      0.00301674010475732127D+00, &
      0.00301674010475732127D+00, &
      0.00301674010475732127D+00, &
      0.00074666452415505266D+00, &
      0.00074666452415505266D+00, &
      0.00074666452415505266D+00, &
      0.00074666452415505266D+00, &
      0.00299748233690233252D+00, &
      0.00299748233690233252D+00, &
      0.00299748233690233252D+00, &
      0.00299748233690233252D+00, &
      0.00144784702606324633D+00, &
      0.00144784702606324633D+00, &
      0.00144784702606324633D+00, &
      0.00144784702606324633D+00, &
      0.00107280962935336075D+00, &
      0.00107280962935336075D+00, &
      0.00107280962935336075D+00, &
      0.00107280962935336075D+00, &
      0.00107280962935336075D+00, &
      0.00107280962935336075D+00, &
      0.00107280962935336075D+00, &
      0.00107280962935336075D+00, &
      0.00483903734837801830D+00, &
      0.00483903734837801830D+00, &
      0.00483903734837801830D+00, &
      0.00483903734837801830D+00, &
      0.00483903734837801830D+00, &
      0.00483903734837801830D+00, &
      0.00483903734837801830D+00, &
      0.00483903734837801830D+00, &
      0.00274709905982003805D+00, &
      0.00274709905982003805D+00, &
      0.00274709905982003805D+00, &
      0.00274709905982003805D+00, &
      0.00274709905982003805D+00, &
      0.00274709905982003805D+00, &
      0.00274709905982003805D+00, &
      0.00274709905982003805D+00, &
      0.00433501688902940795D+00, &
      0.00433501688902940795D+00, &
      0.00433501688902940795D+00, &
      0.00433501688902940795D+00, &
      0.00433501688902940795D+00, &
      0.00433501688902940795D+00, &
      0.00433501688902940795D+00, &
      0.00433501688902940795D+00, &
      0.00271810015984954700D+00, &
      0.00271810015984954700D+00, &
      0.00271810015984954700D+00, &
      0.00271810015984954700D+00, &
      0.00271810015984954700D+00, &
      0.00271810015984954700D+00, &
      0.00271810015984954700D+00, &
      0.00271810015984954700D+00, &
      0.00091685473986105473D+00, &
      0.00091685473986105473D+00, &
      0.00091685473986105473D+00, &
      0.00091685473986105473D+00, &
      0.00091685473986105473D+00, &
      0.00091685473986105473D+00, &
      0.00091685473986105473D+00, &
      0.00091685473986105473D+00, &
      0.00680895101968847279D+00, &
      0.00680895101968847279D+00, &
      0.00680895101968847279D+00, &
      0.00680895101968847279D+00, &
      0.00680895101968847279D+00, &
      0.00680895101968847279D+00, &
      0.00680895101968847279D+00, &
      0.00680895101968847279D+00, &
      0.00040138568966109102D+00, &
      0.00040138568966109102D+00, &
      0.00040138568966109102D+00, &
      0.00040138568966109102D+00, &
      0.00040138568966109102D+00, &
      0.00040138568966109102D+00, &
      0.00040138568966109102D+00, &
      0.00040138568966109102D+00, &
      0.00215500289545228485D+00, &
      0.00215500289545228485D+00, &
      0.00215500289545228485D+00, &
      0.00215500289545228485D+00, &
      0.00215500289545228485D+00, &
      0.00215500289545228485D+00, &
      0.00215500289545228485D+00, &
      0.00215500289545228485D+00, &
      0.00161154257302305293D+00, &
      0.00161154257302305293D+00, &
      0.00161154257302305293D+00, &
      0.00161154257302305293D+00, &
      0.00161154257302305293D+00, &
      0.00161154257302305293D+00, &
      0.00161154257302305293D+00, &
      0.00161154257302305293D+00, &
      0.00027434552667444364D+00, &
      0.00027434552667444364D+00, &
      0.00027434552667444364D+00, &
      0.00027434552667444364D+00, &
      0.00027434552667444364D+00, &
      0.00027434552667444364D+00, &
      0.00027434552667444364D+00, &
      0.00027434552667444364D+00, &
      0.00230008297506218524D+00, &
      0.00230008297506218524D+00, &
      0.00230008297506218524D+00, &
      0.00230008297506218524D+00, &
      0.00230008297506218524D+00, &
      0.00230008297506218524D+00, &
      0.00230008297506218524D+00, &
      0.00230008297506218524D+00, &
      0.00029458357590517713D+00, &
      0.00029458357590517713D+00, &
      0.00029458357590517713D+00, &
      0.00029458357590517713D+00, &
      0.00029458357590517713D+00, &
      0.00029458357590517713D+00, &
      0.00029458357590517713D+00, &
      0.00029458357590517713D+00, &
      0.00060073122167725204D+00, &
      0.00060073122167725204D+00, &
      0.00060073122167725204D+00, &
      0.00060073122167725204D+00, &
      0.00060073122167725204D+00, &
      0.00060073122167725204D+00, &
      0.00060073122167725204D+00, &
      0.00060073122167725204D+00, &
      0.00026201750375773298D+00, &
      0.00026201750375773298D+00, &
      0.00026201750375773298D+00, &
      0.00026201750375773298D+00, &
      0.00026201750375773298D+00, &
      0.00026201750375773298D+00, &
      0.00026201750375773298D+00, &
      0.00026201750375773298D+00, &
      0.00419283249870963659D+00, &
      0.00419283249870963659D+00, &
      0.00419283249870963659D+00, &
      0.00419283249870963659D+00, &
      0.00419283249870963659D+00, &
      0.00419283249870963659D+00, &
      0.00419283249870963659D+00, &
      0.00419283249870963659D+00, &
      0.00173629430868813472D+00, &
      0.00173629430868813472D+00, &
      0.00173629430868813472D+00, &
      0.00173629430868813472D+00, &
      0.00173629430868813472D+00, &
      0.00173629430868813472D+00, &
      0.00173629430868813472D+00, &
      0.00173629430868813472D+00, &
      0.00094391016694345316D+00, &
      0.00094391016694345316D+00, &
      0.00094391016694345316D+00, &
      0.00094391016694345316D+00, &
      0.00094391016694345316D+00, &
      0.00094391016694345316D+00, &
      0.00094391016694345316D+00, &
      0.00094391016694345316D+00, &
      0.00217551360850374046D+00, &
      0.00217551360850374046D+00, &
      0.00217551360850374046D+00, &
      0.00217551360850374046D+00, &
      0.00217551360850374046D+00, &
      0.00217551360850374046D+00, &
      0.00217551360850374046D+00, &
      0.00217551360850374046D+00, &
      0.00128339469018335490D+00, &
      0.00128339469018335490D+00, &
      0.00128339469018335490D+00, &
      0.00128339469018335490D+00, &
      0.00128339469018335490D+00, &
      0.00128339469018335490D+00, &
      0.00128339469018335490D+00, &
      0.00128339469018335490D+00, &
      0.00272615049960569979D+00, &
      0.00272615049960569979D+00, &
      0.00272615049960569979D+00, &
      0.00272615049960569979D+00, &
      0.00272615049960569979D+00, &
      0.00272615049960569979D+00, &
      0.00272615049960569979D+00, &
      0.00272615049960569979D+00, &
      0.00318533960385062025D+00, &
      0.00318533960385062025D+00, &
      0.00318533960385062025D+00, &
      0.00318533960385062025D+00, &
      0.00318533960385062025D+00, &
      0.00318533960385062025D+00, &
      0.00318533960385062025D+00, &
      0.00318533960385062025D+00, &
      0.00266309286664456963D+00, &
      0.00266309286664456963D+00, &
      0.00266309286664456963D+00, &
      0.00266309286664456963D+00, &
      0.00266309286664456963D+00, &
      0.00266309286664456963D+00, &
      0.00266309286664456963D+00, &
      0.00266309286664456963D+00, &
      0.00232666846528824418D+00, &
      0.00232666846528824418D+00, &
      0.00232666846528824418D+00, &
      0.00232666846528824418D+00, &
      0.00232666846528824418D+00, &
      0.00232666846528824418D+00, &
      0.00232666846528824418D+00, &
      0.00232666846528824418D+00, &
      0.00097346311519703153D+00, &
      0.00097346311519703153D+00, &
      0.00097346311519703153D+00, &
      0.00097346311519703153D+00, &
      0.00097346311519703153D+00, &
      0.00097346311519703153D+00, &
      0.00097346311519703153D+00, &
      0.00097346311519703153D+00, &
      0.00348345685139948022D+00, &
      0.00348345685139948022D+00, &
      0.00348345685139948022D+00, &
      0.00348345685139948022D+00, &
      0.00348345685139948022D+00, &
      0.00348345685139948022D+00, &
      0.00348345685139948022D+00, &
      0.00348345685139948022D+00, &
      0.00406195421696304172D+00, &
      0.00406195421696304172D+00, &
      0.00406195421696304172D+00, &
      0.00406195421696304172D+00, &
      0.00406195421696304172D+00, &
      0.00406195421696304172D+00, &
      0.00406195421696304172D+00, &
      0.00406195421696304172D+00, &
      0.00033245593568842263D+00, &
      0.00033245593568842263D+00, &
      0.00033245593568842263D+00, &
      0.00033245593568842263D+00, &
      0.00033245593568842263D+00, &
      0.00033245593568842263D+00, &
      0.00033245593568842263D+00, &
      0.00033245593568842263D+00, &
      0.00092328174175319150D+00, &
      0.00092328174175319150D+00, &
      0.00092328174175319150D+00, &
      0.00092328174175319150D+00, &
      0.00092328174175319150D+00, &
      0.00092328174175319150D+00, &
      0.00092328174175319150D+00, &
      0.00092328174175319150D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
subroutine rule20 ( n, x, y, z, w )

!*****************************************************************************80
!
!! rule20() returns the pyramid quadrature rule of precision 20.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 April 2023
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jan Jaskowiec, Natarajan Sukumar,
!    High order cubature rules for tetrahedra and pyramids,
!    International Journal of Numerical Methods in Engineering,
!    Volume 121, Number 11, pages 2418-2436, 15 June 2020.
!
!  Input:
!
!    integer n: the number of quadrature points.
!
!  Output:
!
!    real ( kind = rk ) x[n], y[n], z[n]: the coordinates of
!    quadrature points.
!
!    real ( kind = rk ) w[n]: the quadrature weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00)

  integer n
  integer, parameter :: n_save = 489

  real ( kind = rk ) x(n)
  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &
      0.00000000000000000000D+00, &
      0.48169354344312520499D+00, &
      0.00000000000000000000D+00, &
     -0.48169354344312520499D+00, &
      0.00000000000000000000D+00, &
      0.69858038869870031640D+00, &
      0.00000000000000000000D+00, &
     -0.69858038869870031640D+00, &
      0.00000000000000000000D+00, &
      0.53468657951770970360D+00, &
      0.00000000000000000000D+00, &
     -0.53468657951770970360D+00, &
      0.00000000000000000000D+00, &
      0.42228380878442672852D+00, &
      0.00000000000000000000D+00, &
     -0.42228380878442672852D+00, &
      0.00000000000000000000D+00, &
      0.30680451888779119995D+00, &
      0.00000000000000000000D+00, &
     -0.30680451888779119995D+00, &
      0.00000000000000000000D+00, &
      0.22877498497689821577D+00, &
      0.00000000000000000000D+00, &
     -0.22877498497689821577D+00, &
      0.00000000000000000000D+00, &
      0.55389687412807175892D+00, &
      0.00000000000000000000D+00, &
     -0.55389687412807175892D+00, &
      0.00000000000000000000D+00, &
      0.38489228042331330437D+00, &
      0.00000000000000000000D+00, &
     -0.38489228042331330437D+00, &
      0.00000000000000000000D+00, &
      0.25711056873844140291D+00, &
      0.00000000000000000000D+00, &
     -0.25711056873844140291D+00, &
      0.00000000000000000000D+00, &
      0.94954208334711565076D+00, &
      0.00000000000000000000D+00, &
     -0.94954208334711565076D+00, &
      0.00000000000000000000D+00, &
      0.09635090648066924057D+00, &
      0.00000000000000000000D+00, &
     -0.09635090648066924057D+00, &
      0.00000000000000000000D+00, &
      0.91044931261818451418D+00, &
      0.00000000000000000000D+00, &
     -0.91044931261818451418D+00, &
      0.00000000000000000000D+00, &
      0.26830824357428106897D+00, &
      0.00000000000000000000D+00, &
     -0.26830824357428106897D+00, &
      0.00000000000000000000D+00, &
      0.82744656686299722370D+00, &
      0.00000000000000000000D+00, &
     -0.82744656686299722370D+00, &
      0.00000000000000000000D+00, &
      0.03789270083108932374D+00, &
      0.00000000000000000000D+00, &
     -0.03789270083108932374D+00, &
      0.00000000000000000000D+00, &
      0.57352052408158127328D+00, &
      0.00000000000000000000D+00, &
     -0.57352052408158127328D+00, &
      0.00000000000000000000D+00, &
      0.74231157227003041754D+00, &
      0.00000000000000000000D+00, &
     -0.74231157227003041754D+00, &
      0.00000000000000000000D+00, &
      0.28469096413200573048D+00, &
      0.00000000000000000000D+00, &
     -0.28469096413200573048D+00, &
      0.00000000000000000000D+00, &
      0.11281627950690893691D+00, &
      0.00000000000000000000D+00, &
     -0.11281627950690893691D+00, &
      0.00000000000000000000D+00, &
      0.45356762399670425001D+00, &
     -0.45356762399670425001D+00, &
      0.45356762399670425001D+00, &
     -0.45356762399670425001D+00, &
      0.33361103852611595499D+00, &
     -0.33361103852611595499D+00, &
      0.33361103852611595499D+00, &
     -0.33361103852611595499D+00, &
      0.37733375278222114346D+00, &
     -0.37733375278222114346D+00, &
      0.37733375278222114346D+00, &
     -0.37733375278222114346D+00, &
      0.30611174916927796907D+00, &
     -0.30611174916927796907D+00, &
      0.30611174916927796907D+00, &
     -0.30611174916927796907D+00, &
      0.29044600765245415230D+00, &
     -0.29044600765245415230D+00, &
      0.29044600765245415230D+00, &
     -0.29044600765245415230D+00, &
      0.59141068840923671779D+00, &
     -0.59141068840923671779D+00, &
      0.59141068840923671779D+00, &
     -0.59141068840923671779D+00, &
      0.76605496491352953470D+00, &
     -0.76605496491352953470D+00, &
      0.76605496491352953470D+00, &
     -0.76605496491352953470D+00, &
      0.44962169263179840861D+00, &
     -0.44962169263179840861D+00, &
      0.44962169263179840861D+00, &
     -0.44962169263179840861D+00, &
      0.61256910754587590162D+00, &
     -0.61256910754587590162D+00, &
      0.61256910754587590162D+00, &
     -0.61256910754587590162D+00, &
      0.96517775637880776074D+00, &
     -0.96517775637880776074D+00, &
      0.96517775637880776074D+00, &
     -0.96517775637880776074D+00, &
      0.18509443085238544424D+00, &
     -0.18509443085238544424D+00, &
      0.18509443085238544424D+00, &
     -0.18509443085238544424D+00, &
      0.13053399140728452754D+00, &
     -0.13053399140728452754D+00, &
      0.13053399140728452754D+00, &
     -0.13053399140728452754D+00, &
      0.08713330755302489683D+00, &
     -0.08713330755302489683D+00, &
      0.08713330755302489683D+00, &
     -0.08713330755302489683D+00, &
      0.11626579551616780805D+00, &
     -0.11626579551616780805D+00, &
      0.11626579551616780805D+00, &
     -0.11626579551616780805D+00, &
      0.03353132384174680597D+00, &
     -0.03353132384174680597D+00, &
      0.03353132384174680597D+00, &
     -0.03353132384174680597D+00, &
      0.63007044649757071308D+00, &
     -0.63007044649757071308D+00, &
      0.63007044649757071308D+00, &
     -0.63007044649757071308D+00, &
      0.10081553511278160129D+00, &
     -0.10081553511278160129D+00, &
      0.10081553511278160129D+00, &
     -0.10081553511278160129D+00, &
      0.09712721271803037570D+00, &
     -0.09712721271803037570D+00, &
      0.09712721271803037570D+00, &
     -0.09712721271803037570D+00, &
      0.65036744612024988133D+00, &
     -0.65036744612024988133D+00, &
      0.65036744612024988133D+00, &
     -0.65036744612024988133D+00, &
      0.78816705553909693904D+00, &
     -0.78816705553909693904D+00, &
      0.78816705553909693904D+00, &
     -0.78816705553909693904D+00, &
      0.46326102535176988395D+00, &
     -0.46326102535176988395D+00, &
      0.46326102535176988395D+00, &
     -0.46326102535176988395D+00, &
      0.87426225625453524160D+00, &
     -0.87426225625453524160D+00, &
      0.87426225625453524160D+00, &
     -0.87426225625453524160D+00, &
      0.09693365579654888986D+00, &
     -0.09693365579654888986D+00, &
      0.09693365579654888986D+00, &
     -0.09693365579654888986D+00, &
      0.52651234369915456135D+00, &
     -0.52651234369915456135D+00, &
      0.52651234369915456135D+00, &
     -0.52651234369915456135D+00, &
      0.21058018350559978837D+00, &
     -0.21058018350559978837D+00, &
      0.21058018350559978837D+00, &
     -0.21058018350559978837D+00, &
      0.48341786018142590686D+00, &
     -0.48341786018142590686D+00, &
      0.48341786018142590686D+00, &
     -0.48341786018142590686D+00, &
      0.27232398053992901144D+00, &
     -0.27232398053992901144D+00, &
      0.27232398053992901144D+00, &
     -0.27232398053992901144D+00, &
      0.88748390701146573356D+00, &
     -0.88748390701146573356D+00, &
      0.88748390701146573356D+00, &
     -0.88748390701146573356D+00, &
      0.34004340177081920915D+00, &
     -0.34004340177081920915D+00, &
      0.34004340177081920915D+00, &
     -0.34004340177081920915D+00, &
      0.24657134521158877161D+00, &
     -0.24657134521158877161D+00, &
      0.24657134521158877161D+00, &
     -0.24657134521158877161D+00, &
      0.72032611852462202773D+00, &
     -0.72032611852462202773D+00, &
      0.72032611852462202773D+00, &
     -0.72032611852462202773D+00, &
      0.86934672471843543740D+00, &
      0.66053776743762293577D+00, &
     -0.86934672471843543740D+00, &
      0.66053776743762293577D+00, &
      0.86934672471843543740D+00, &
     -0.66053776743762293577D+00, &
     -0.86934672471843543740D+00, &
     -0.66053776743762293577D+00, &
      0.51657003644145271792D+00, &
      0.12681144322431647797D+00, &
     -0.51657003644145271792D+00, &
      0.12681144322431647797D+00, &
      0.51657003644145271792D+00, &
     -0.12681144322431647797D+00, &
     -0.51657003644145271792D+00, &
     -0.12681144322431647797D+00, &
      0.21419492835914627493D+00, &
      0.53514205206517340141D+00, &
     -0.21419492835914627493D+00, &
      0.53514205206517340141D+00, &
      0.21419492835914627493D+00, &
     -0.53514205206517340141D+00, &
     -0.21419492835914627493D+00, &
     -0.53514205206517340141D+00, &
      0.38606155250360069120D+00, &
      0.14309874624364349316D+00, &
     -0.38606155250360069120D+00, &
      0.14309874624364349316D+00, &
      0.38606155250360069120D+00, &
     -0.14309874624364349316D+00, &
     -0.38606155250360069120D+00, &
     -0.14309874624364349316D+00, &
      0.31482652231165053625D+00, &
      0.57887349054423364869D+00, &
     -0.31482652231165053625D+00, &
      0.57887349054423364869D+00, &
      0.31482652231165053625D+00, &
     -0.57887349054423364869D+00, &
     -0.31482652231165053625D+00, &
     -0.57887349054423364869D+00, &
      0.13337762200744199270D+00, &
      0.86857348087766383937D+00, &
     -0.13337762200744199270D+00, &
      0.86857348087766383937D+00, &
      0.13337762200744199270D+00, &
     -0.86857348087766383937D+00, &
     -0.13337762200744199270D+00, &
     -0.86857348087766383937D+00, &
      0.89334514265657827270D+00, &
      0.29700939611864535239D+00, &
     -0.89334514265657827270D+00, &
      0.29700939611864535239D+00, &
      0.89334514265657827270D+00, &
     -0.29700939611864535239D+00, &
     -0.89334514265657827270D+00, &
     -0.29700939611864535239D+00, &
      0.64207321281802853807D+00, &
      0.21859043909362887992D+00, &
     -0.64207321281802853807D+00, &
      0.21859043909362887992D+00, &
      0.64207321281802853807D+00, &
     -0.21859043909362887992D+00, &
     -0.64207321281802853807D+00, &
     -0.21859043909362887992D+00, &
      0.90155669259026749440D+00, &
      0.28504501057958070431D+00, &
     -0.90155669259026749440D+00, &
      0.28504501057958070431D+00, &
      0.90155669259026749440D+00, &
     -0.28504501057958070431D+00, &
     -0.90155669259026749440D+00, &
     -0.28504501057958070431D+00, &
      0.17950781738619853156D+00, &
      0.67966060973500286302D+00, &
     -0.17950781738619853156D+00, &
      0.67966060973500286302D+00, &
      0.17950781738619853156D+00, &
     -0.67966060973500286302D+00, &
     -0.17950781738619853156D+00, &
     -0.67966060973500286302D+00, &
      0.52877164545115751260D+00, &
      0.37525379145011561466D+00, &
     -0.52877164545115751260D+00, &
      0.37525379145011561466D+00, &
      0.52877164545115751260D+00, &
     -0.37525379145011561466D+00, &
     -0.52877164545115751260D+00, &
     -0.37525379145011561466D+00, &
      0.85684776838109344421D+00, &
      0.70614261226539476457D+00, &
     -0.85684776838109344421D+00, &
      0.70614261226539476457D+00, &
      0.85684776838109344421D+00, &
     -0.70614261226539476457D+00, &
     -0.85684776838109344421D+00, &
     -0.70614261226539476457D+00, &
      0.23744649710478388238D+00, &
      0.09039264484495425356D+00, &
     -0.23744649710478388238D+00, &
      0.09039264484495425356D+00, &
      0.23744649710478388238D+00, &
     -0.09039264484495425356D+00, &
     -0.23744649710478388238D+00, &
     -0.09039264484495425356D+00, &
      0.94065987701065645332D+00, &
      0.74375119609914763785D+00, &
     -0.94065987701065645332D+00, &
      0.74375119609914763785D+00, &
      0.94065987701065645332D+00, &
     -0.74375119609914763785D+00, &
     -0.94065987701065645332D+00, &
     -0.74375119609914763785D+00, &
      0.96544985579484265958D+00, &
      0.55816486672448317741D+00, &
     -0.96544985579484265958D+00, &
      0.55816486672448317741D+00, &
      0.96544985579484265958D+00, &
     -0.55816486672448317741D+00, &
     -0.96544985579484265958D+00, &
     -0.55816486672448317741D+00, &
      0.11859977159394578805D+00, &
      0.19730553454891092136D+00, &
     -0.11859977159394578805D+00, &
      0.19730553454891092136D+00, &
      0.11859977159394578805D+00, &
     -0.19730553454891092136D+00, &
     -0.11859977159394578805D+00, &
     -0.19730553454891092136D+00, &
      0.43836828746731176798D+00, &
      0.68010534456895455069D+00, &
     -0.43836828746731176798D+00, &
      0.68010534456895455069D+00, &
      0.43836828746731176798D+00, &
     -0.68010534456895455069D+00, &
     -0.43836828746731176798D+00, &
     -0.68010534456895455069D+00, &
      0.85015015838954666183D+00, &
      0.57767859879985306026D+00, &
     -0.85015015838954666183D+00, &
      0.57767859879985306026D+00, &
      0.85015015838954666183D+00, &
     -0.57767859879985306026D+00, &
     -0.85015015838954666183D+00, &
     -0.57767859879985306026D+00, &
      0.40977883537749759668D+00, &
      0.59926491678603310831D+00, &
     -0.40977883537749759668D+00, &
      0.59926491678603310831D+00, &
      0.40977883537749759668D+00, &
     -0.59926491678603310831D+00, &
     -0.40977883537749759668D+00, &
     -0.59926491678603310831D+00, &
      0.62540250575452105419D+00, &
      0.77938113963987465382D+00, &
     -0.62540250575452105419D+00, &
      0.77938113963987465382D+00, &
      0.62540250575452105419D+00, &
     -0.77938113963987465382D+00, &
     -0.62540250575452105419D+00, &
     -0.77938113963987465382D+00, &
      0.16587947711557049502D+00, &
      0.65137881658599783297D+00, &
     -0.16587947711557049502D+00, &
      0.65137881658599783297D+00, &
      0.16587947711557049502D+00, &
     -0.65137881658599783297D+00, &
     -0.16587947711557049502D+00, &
     -0.65137881658599783297D+00, &
      0.50625591443667894431D+00, &
      0.24708268775667227568D+00, &
     -0.50625591443667894431D+00, &
      0.24708268775667227568D+00, &
      0.50625591443667894431D+00, &
     -0.24708268775667227568D+00, &
     -0.50625591443667894431D+00, &
     -0.24708268775667227568D+00, &
      0.74292936328073289065D+00, &
      0.16823593457058672040D+00, &
     -0.74292936328073289065D+00, &
      0.16823593457058672040D+00, &
      0.74292936328073289065D+00, &
     -0.16823593457058672040D+00, &
     -0.74292936328073289065D+00, &
     -0.16823593457058672040D+00, &
      0.01867609649526593210D+00, &
      0.25362045997818172260D+00, &
     -0.01867609649526593210D+00, &
      0.25362045997818172260D+00, &
      0.01867609649526593210D+00, &
     -0.25362045997818172260D+00, &
     -0.01867609649526593210D+00, &
     -0.25362045997818172260D+00, &
      0.26086005279281720970D+00, &
      0.40946055588941798753D+00, &
     -0.26086005279281720970D+00, &
      0.40946055588941798753D+00, &
      0.26086005279281720970D+00, &
     -0.40946055588941798753D+00, &
     -0.26086005279281720970D+00, &
     -0.40946055588941798753D+00, &
      0.79388083129835063101D+00, &
      0.46918818998244915530D+00, &
     -0.79388083129835063101D+00, &
      0.46918818998244915530D+00, &
      0.79388083129835063101D+00, &
     -0.46918818998244915530D+00, &
     -0.79388083129835063101D+00, &
     -0.46918818998244915530D+00, &
      0.79301165544060303603D+00, &
      0.29864577765030031475D+00, &
     -0.79301165544060303603D+00, &
      0.29864577765030031475D+00, &
      0.79301165544060303603D+00, &
     -0.29864577765030031475D+00, &
     -0.79301165544060303603D+00, &
     -0.29864577765030031475D+00, &
      0.44929996069833155747D+00, &
      0.16967965727402492537D+00, &
     -0.44929996069833155747D+00, &
      0.16967965727402492537D+00, &
      0.44929996069833155747D+00, &
     -0.16967965727402492537D+00, &
     -0.44929996069833155747D+00, &
     -0.16967965727402492537D+00, &
      0.76873831458793862037D+00, &
      0.39677450550481596636D+00, &
     -0.76873831458793862037D+00, &
      0.39677450550481596636D+00, &
      0.76873831458793862037D+00, &
     -0.39677450550481596636D+00, &
     -0.76873831458793862037D+00, &
     -0.39677450550481596636D+00, &
      0.97859975977571478367D+00, &
      0.19094705775922465874D+00, &
     -0.97859975977571478367D+00, &
      0.19094705775922465874D+00, &
      0.97859975977571478367D+00, &
     -0.19094705775922465874D+00, &
     -0.97859975977571478367D+00, &
     -0.19094705775922465874D+00, &
      0.21358909276843840441D+00, &
      0.36041615296155676829D+00, &
     -0.21358909276843840441D+00, &
      0.36041615296155676829D+00, &
      0.21358909276843840441D+00, &
     -0.36041615296155676829D+00, &
     -0.21358909276843840441D+00, &
     -0.36041615296155676829D+00, &
      0.79625257738750132575D+00, &
      0.15628526976368214974D+00, &
     -0.79625257738750132575D+00, &
      0.15628526976368214974D+00, &
      0.79625257738750132575D+00, &
     -0.15628526976368214974D+00, &
     -0.79625257738750132575D+00, &
     -0.15628526976368214974D+00, &
      0.91939856562697519493D+00, &
      0.41486074510619452838D+00, &
     -0.91939856562697519493D+00, &
      0.41486074510619452838D+00, &
      0.91939856562697519493D+00, &
     -0.41486074510619452838D+00, &
     -0.91939856562697519493D+00, &
     -0.41486074510619452838D+00, &
      0.51507330251609151350D+00, &
      0.68962930070823635909D+00, &
     -0.51507330251609151350D+00, &
      0.68962930070823635909D+00, &
      0.51507330251609151350D+00, &
     -0.68962930070823635909D+00, &
     -0.51507330251609151350D+00, &
     -0.68962930070823635909D+00, &
      0.96361213891168562284D+00, &
      0.85957679427925348659D+00, &
     -0.96361213891168562284D+00, &
      0.85957679427925348659D+00, &
      0.96361213891168562284D+00, &
     -0.85957679427925348659D+00, &
     -0.96361213891168562284D+00, &
     -0.85957679427925348659D+00, &
      0.72140711490936071382D+00, &
      0.49209543700208141503D+00, &
     -0.72140711490936071382D+00, &
      0.49209543700208141503D+00, &
      0.72140711490936071382D+00, &
     -0.49209543700208141503D+00, &
     -0.72140711490936071382D+00, &
     -0.49209543700208141503D+00 /)

  real ( kind = rk ) y(n)
  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &
      0.00000000000000000000D+00, &
      0.00000000000000000000D+00, &
      0.48169354344312520499D+00, &
      0.00000000000000000000D+00, &
     -0.48169354344312520499D+00, &
      0.00000000000000000000D+00, &
      0.69858038869870031640D+00, &
      0.00000000000000000000D+00, &
     -0.69858038869870031640D+00, &
      0.00000000000000000000D+00, &
      0.53468657951770970360D+00, &
      0.00000000000000000000D+00, &
     -0.53468657951770970360D+00, &
      0.00000000000000000000D+00, &
      0.42228380878442672852D+00, &
      0.00000000000000000000D+00, &
     -0.42228380878442672852D+00, &
      0.00000000000000000000D+00, &
      0.30680451888779119995D+00, &
      0.00000000000000000000D+00, &
     -0.30680451888779119995D+00, &
      0.00000000000000000000D+00, &
      0.22877498497689821577D+00, &
      0.00000000000000000000D+00, &
     -0.22877498497689821577D+00, &
      0.00000000000000000000D+00, &
      0.55389687412807175892D+00, &
      0.00000000000000000000D+00, &
     -0.55389687412807175892D+00, &
      0.00000000000000000000D+00, &
      0.38489228042331330437D+00, &
      0.00000000000000000000D+00, &
     -0.38489228042331330437D+00, &
      0.00000000000000000000D+00, &
      0.25711056873844140291D+00, &
      0.00000000000000000000D+00, &
     -0.25711056873844140291D+00, &
      0.00000000000000000000D+00, &
      0.94954208334711565076D+00, &
      0.00000000000000000000D+00, &
     -0.94954208334711565076D+00, &
      0.00000000000000000000D+00, &
      0.09635090648066924057D+00, &
      0.00000000000000000000D+00, &
     -0.09635090648066924057D+00, &
      0.00000000000000000000D+00, &
      0.91044931261818451418D+00, &
      0.00000000000000000000D+00, &
     -0.91044931261818451418D+00, &
      0.00000000000000000000D+00, &
      0.26830824357428106897D+00, &
      0.00000000000000000000D+00, &
     -0.26830824357428106897D+00, &
      0.00000000000000000000D+00, &
      0.82744656686299722370D+00, &
      0.00000000000000000000D+00, &
     -0.82744656686299722370D+00, &
      0.00000000000000000000D+00, &
      0.03789270083108932374D+00, &
      0.00000000000000000000D+00, &
     -0.03789270083108932374D+00, &
      0.00000000000000000000D+00, &
      0.57352052408158127328D+00, &
      0.00000000000000000000D+00, &
     -0.57352052408158127328D+00, &
      0.00000000000000000000D+00, &
      0.74231157227003041754D+00, &
      0.00000000000000000000D+00, &
     -0.74231157227003041754D+00, &
      0.00000000000000000000D+00, &
      0.28469096413200573048D+00, &
      0.00000000000000000000D+00, &
     -0.28469096413200573048D+00, &
      0.00000000000000000000D+00, &
      0.11281627950690893691D+00, &
      0.00000000000000000000D+00, &
     -0.11281627950690893691D+00, &
      0.45356762399670425001D+00, &
      0.45356762399670425001D+00, &
     -0.45356762399670425001D+00, &
     -0.45356762399670425001D+00, &
      0.33361103852611595499D+00, &
      0.33361103852611595499D+00, &
     -0.33361103852611595499D+00, &
     -0.33361103852611595499D+00, &
      0.37733375278222114346D+00, &
      0.37733375278222114346D+00, &
     -0.37733375278222114346D+00, &
     -0.37733375278222114346D+00, &
      0.30611174916927796907D+00, &
      0.30611174916927796907D+00, &
     -0.30611174916927796907D+00, &
     -0.30611174916927796907D+00, &
      0.29044600765245415230D+00, &
      0.29044600765245415230D+00, &
     -0.29044600765245415230D+00, &
     -0.29044600765245415230D+00, &
      0.59141068840923671779D+00, &
      0.59141068840923671779D+00, &
     -0.59141068840923671779D+00, &
     -0.59141068840923671779D+00, &
      0.76605496491352953470D+00, &
      0.76605496491352953470D+00, &
     -0.76605496491352953470D+00, &
     -0.76605496491352953470D+00, &
      0.44962169263179840861D+00, &
      0.44962169263179840861D+00, &
     -0.44962169263179840861D+00, &
     -0.44962169263179840861D+00, &
      0.61256910754587590162D+00, &
      0.61256910754587590162D+00, &
     -0.61256910754587590162D+00, &
     -0.61256910754587590162D+00, &
      0.96517775637880776074D+00, &
      0.96517775637880776074D+00, &
     -0.96517775637880776074D+00, &
     -0.96517775637880776074D+00, &
      0.18509443085238544424D+00, &
      0.18509443085238544424D+00, &
     -0.18509443085238544424D+00, &
     -0.18509443085238544424D+00, &
      0.13053399140728452754D+00, &
      0.13053399140728452754D+00, &
     -0.13053399140728452754D+00, &
     -0.13053399140728452754D+00, &
      0.08713330755302489683D+00, &
      0.08713330755302489683D+00, &
     -0.08713330755302489683D+00, &
     -0.08713330755302489683D+00, &
      0.11626579551616780805D+00, &
      0.11626579551616780805D+00, &
     -0.11626579551616780805D+00, &
     -0.11626579551616780805D+00, &
      0.03353132384174680597D+00, &
      0.03353132384174680597D+00, &
     -0.03353132384174680597D+00, &
     -0.03353132384174680597D+00, &
      0.63007044649757071308D+00, &
      0.63007044649757071308D+00, &
     -0.63007044649757071308D+00, &
     -0.63007044649757071308D+00, &
      0.10081553511278160129D+00, &
      0.10081553511278160129D+00, &
     -0.10081553511278160129D+00, &
     -0.10081553511278160129D+00, &
      0.09712721271803037570D+00, &
      0.09712721271803037570D+00, &
     -0.09712721271803037570D+00, &
     -0.09712721271803037570D+00, &
      0.65036744612024988133D+00, &
      0.65036744612024988133D+00, &
     -0.65036744612024988133D+00, &
     -0.65036744612024988133D+00, &
      0.78816705553909693904D+00, &
      0.78816705553909693904D+00, &
     -0.78816705553909693904D+00, &
     -0.78816705553909693904D+00, &
      0.46326102535176988395D+00, &
      0.46326102535176988395D+00, &
     -0.46326102535176988395D+00, &
     -0.46326102535176988395D+00, &
      0.87426225625453524160D+00, &
      0.87426225625453524160D+00, &
     -0.87426225625453524160D+00, &
     -0.87426225625453524160D+00, &
      0.09693365579654888986D+00, &
      0.09693365579654888986D+00, &
     -0.09693365579654888986D+00, &
     -0.09693365579654888986D+00, &
      0.52651234369915456135D+00, &
      0.52651234369915456135D+00, &
     -0.52651234369915456135D+00, &
     -0.52651234369915456135D+00, &
      0.21058018350559978837D+00, &
      0.21058018350559978837D+00, &
     -0.21058018350559978837D+00, &
     -0.21058018350559978837D+00, &
      0.48341786018142590686D+00, &
      0.48341786018142590686D+00, &
     -0.48341786018142590686D+00, &
     -0.48341786018142590686D+00, &
      0.27232398053992901144D+00, &
      0.27232398053992901144D+00, &
     -0.27232398053992901144D+00, &
     -0.27232398053992901144D+00, &
      0.88748390701146573356D+00, &
      0.88748390701146573356D+00, &
     -0.88748390701146573356D+00, &
     -0.88748390701146573356D+00, &
      0.34004340177081920915D+00, &
      0.34004340177081920915D+00, &
     -0.34004340177081920915D+00, &
     -0.34004340177081920915D+00, &
      0.24657134521158877161D+00, &
      0.24657134521158877161D+00, &
     -0.24657134521158877161D+00, &
     -0.24657134521158877161D+00, &
      0.72032611852462202773D+00, &
      0.72032611852462202773D+00, &
     -0.72032611852462202773D+00, &
     -0.72032611852462202773D+00, &
      0.66053776743762293577D+00, &
      0.86934672471843543740D+00, &
      0.66053776743762293577D+00, &
     -0.86934672471843543740D+00, &
     -0.66053776743762293577D+00, &
      0.86934672471843543740D+00, &
     -0.66053776743762293577D+00, &
     -0.86934672471843543740D+00, &
      0.12681144322431647797D+00, &
      0.51657003644145271792D+00, &
      0.12681144322431647797D+00, &
     -0.51657003644145271792D+00, &
     -0.12681144322431647797D+00, &
      0.51657003644145271792D+00, &
     -0.12681144322431647797D+00, &
     -0.51657003644145271792D+00, &
      0.53514205206517340141D+00, &
      0.21419492835914627493D+00, &
      0.53514205206517340141D+00, &
     -0.21419492835914627493D+00, &
     -0.53514205206517340141D+00, &
      0.21419492835914627493D+00, &
     -0.53514205206517340141D+00, &
     -0.21419492835914627493D+00, &
      0.14309874624364349316D+00, &
      0.38606155250360069120D+00, &
      0.14309874624364349316D+00, &
     -0.38606155250360069120D+00, &
     -0.14309874624364349316D+00, &
      0.38606155250360069120D+00, &
     -0.14309874624364349316D+00, &
     -0.38606155250360069120D+00, &
      0.57887349054423364869D+00, &
      0.31482652231165053625D+00, &
      0.57887349054423364869D+00, &
     -0.31482652231165053625D+00, &
     -0.57887349054423364869D+00, &
      0.31482652231165053625D+00, &
     -0.57887349054423364869D+00, &
     -0.31482652231165053625D+00, &
      0.86857348087766383937D+00, &
      0.13337762200744199270D+00, &
      0.86857348087766383937D+00, &
     -0.13337762200744199270D+00, &
     -0.86857348087766383937D+00, &
      0.13337762200744199270D+00, &
     -0.86857348087766383937D+00, &
     -0.13337762200744199270D+00, &
      0.29700939611864535239D+00, &
      0.89334514265657827270D+00, &
      0.29700939611864535239D+00, &
     -0.89334514265657827270D+00, &
     -0.29700939611864535239D+00, &
      0.89334514265657827270D+00, &
     -0.29700939611864535239D+00, &
     -0.89334514265657827270D+00, &
      0.21859043909362887992D+00, &
      0.64207321281802853807D+00, &
      0.21859043909362887992D+00, &
     -0.64207321281802853807D+00, &
     -0.21859043909362887992D+00, &
      0.64207321281802853807D+00, &
     -0.21859043909362887992D+00, &
     -0.64207321281802853807D+00, &
      0.28504501057958070431D+00, &
      0.90155669259026749440D+00, &
      0.28504501057958070431D+00, &
     -0.90155669259026749440D+00, &
     -0.28504501057958070431D+00, &
      0.90155669259026749440D+00, &
     -0.28504501057958070431D+00, &
     -0.90155669259026749440D+00, &
      0.67966060973500286302D+00, &
      0.17950781738619853156D+00, &
      0.67966060973500286302D+00, &
     -0.17950781738619853156D+00, &
     -0.67966060973500286302D+00, &
      0.17950781738619853156D+00, &
     -0.67966060973500286302D+00, &
     -0.17950781738619853156D+00, &
      0.37525379145011561466D+00, &
      0.52877164545115751260D+00, &
      0.37525379145011561466D+00, &
     -0.52877164545115751260D+00, &
     -0.37525379145011561466D+00, &
      0.52877164545115751260D+00, &
     -0.37525379145011561466D+00, &
     -0.52877164545115751260D+00, &
      0.70614261226539476457D+00, &
      0.85684776838109344421D+00, &
      0.70614261226539476457D+00, &
     -0.85684776838109344421D+00, &
     -0.70614261226539476457D+00, &
      0.85684776838109344421D+00, &
     -0.70614261226539476457D+00, &
     -0.85684776838109344421D+00, &
      0.09039264484495425356D+00, &
      0.23744649710478388238D+00, &
      0.09039264484495425356D+00, &
     -0.23744649710478388238D+00, &
     -0.09039264484495425356D+00, &
      0.23744649710478388238D+00, &
     -0.09039264484495425356D+00, &
     -0.23744649710478388238D+00, &
      0.74375119609914763785D+00, &
      0.94065987701065645332D+00, &
      0.74375119609914763785D+00, &
     -0.94065987701065645332D+00, &
     -0.74375119609914763785D+00, &
      0.94065987701065645332D+00, &
     -0.74375119609914763785D+00, &
     -0.94065987701065645332D+00, &
      0.55816486672448317741D+00, &
      0.96544985579484265958D+00, &
      0.55816486672448317741D+00, &
     -0.96544985579484265958D+00, &
     -0.55816486672448317741D+00, &
      0.96544985579484265958D+00, &
     -0.55816486672448317741D+00, &
     -0.96544985579484265958D+00, &
      0.19730553454891092136D+00, &
      0.11859977159394578805D+00, &
      0.19730553454891092136D+00, &
     -0.11859977159394578805D+00, &
     -0.19730553454891092136D+00, &
      0.11859977159394578805D+00, &
     -0.19730553454891092136D+00, &
     -0.11859977159394578805D+00, &
      0.68010534456895455069D+00, &
      0.43836828746731176798D+00, &
      0.68010534456895455069D+00, &
     -0.43836828746731176798D+00, &
     -0.68010534456895455069D+00, &
      0.43836828746731176798D+00, &
     -0.68010534456895455069D+00, &
     -0.43836828746731176798D+00, &
      0.57767859879985306026D+00, &
      0.85015015838954666183D+00, &
      0.57767859879985306026D+00, &
     -0.85015015838954666183D+00, &
     -0.57767859879985306026D+00, &
      0.85015015838954666183D+00, &
     -0.57767859879985306026D+00, &
     -0.85015015838954666183D+00, &
      0.59926491678603310831D+00, &
      0.40977883537749759668D+00, &
      0.59926491678603310831D+00, &
     -0.40977883537749759668D+00, &
     -0.59926491678603310831D+00, &
      0.40977883537749759668D+00, &
     -0.59926491678603310831D+00, &
     -0.40977883537749759668D+00, &
      0.77938113963987465382D+00, &
      0.62540250575452105419D+00, &
      0.77938113963987465382D+00, &
     -0.62540250575452105419D+00, &
     -0.77938113963987465382D+00, &
      0.62540250575452105419D+00, &
     -0.77938113963987465382D+00, &
     -0.62540250575452105419D+00, &
      0.65137881658599783297D+00, &
      0.16587947711557049502D+00, &
      0.65137881658599783297D+00, &
     -0.16587947711557049502D+00, &
     -0.65137881658599783297D+00, &
      0.16587947711557049502D+00, &
     -0.65137881658599783297D+00, &
     -0.16587947711557049502D+00, &
      0.24708268775667227568D+00, &
      0.50625591443667894431D+00, &
      0.24708268775667227568D+00, &
     -0.50625591443667894431D+00, &
     -0.24708268775667227568D+00, &
      0.50625591443667894431D+00, &
     -0.24708268775667227568D+00, &
     -0.50625591443667894431D+00, &
      0.16823593457058672040D+00, &
      0.74292936328073289065D+00, &
      0.16823593457058672040D+00, &
     -0.74292936328073289065D+00, &
     -0.16823593457058672040D+00, &
      0.74292936328073289065D+00, &
     -0.16823593457058672040D+00, &
     -0.74292936328073289065D+00, &
      0.25362045997818172260D+00, &
      0.01867609649526593210D+00, &
      0.25362045997818172260D+00, &
     -0.01867609649526593210D+00, &
     -0.25362045997818172260D+00, &
      0.01867609649526593210D+00, &
     -0.25362045997818172260D+00, &
     -0.01867609649526593210D+00, &
      0.40946055588941798753D+00, &
      0.26086005279281720970D+00, &
      0.40946055588941798753D+00, &
     -0.26086005279281720970D+00, &
     -0.40946055588941798753D+00, &
      0.26086005279281720970D+00, &
     -0.40946055588941798753D+00, &
     -0.26086005279281720970D+00, &
      0.46918818998244915530D+00, &
      0.79388083129835063101D+00, &
      0.46918818998244915530D+00, &
     -0.79388083129835063101D+00, &
     -0.46918818998244915530D+00, &
      0.79388083129835063101D+00, &
     -0.46918818998244915530D+00, &
     -0.79388083129835063101D+00, &
      0.29864577765030031475D+00, &
      0.79301165544060303603D+00, &
      0.29864577765030031475D+00, &
     -0.79301165544060303603D+00, &
     -0.29864577765030031475D+00, &
      0.79301165544060303603D+00, &
     -0.29864577765030031475D+00, &
     -0.79301165544060303603D+00, &
      0.16967965727402492537D+00, &
      0.44929996069833155747D+00, &
      0.16967965727402492537D+00, &
     -0.44929996069833155747D+00, &
     -0.16967965727402492537D+00, &
      0.44929996069833155747D+00, &
     -0.16967965727402492537D+00, &
     -0.44929996069833155747D+00, &
      0.39677450550481596636D+00, &
      0.76873831458793862037D+00, &
      0.39677450550481596636D+00, &
     -0.76873831458793862037D+00, &
     -0.39677450550481596636D+00, &
      0.76873831458793862037D+00, &
     -0.39677450550481596636D+00, &
     -0.76873831458793862037D+00, &
      0.19094705775922465874D+00, &
      0.97859975977571478367D+00, &
      0.19094705775922465874D+00, &
     -0.97859975977571478367D+00, &
     -0.19094705775922465874D+00, &
      0.97859975977571478367D+00, &
     -0.19094705775922465874D+00, &
     -0.97859975977571478367D+00, &
      0.36041615296155676829D+00, &
      0.21358909276843840441D+00, &
      0.36041615296155676829D+00, &
     -0.21358909276843840441D+00, &
     -0.36041615296155676829D+00, &
      0.21358909276843840441D+00, &
     -0.36041615296155676829D+00, &
     -0.21358909276843840441D+00, &
      0.15628526976368214974D+00, &
      0.79625257738750132575D+00, &
      0.15628526976368214974D+00, &
     -0.79625257738750132575D+00, &
     -0.15628526976368214974D+00, &
      0.79625257738750132575D+00, &
     -0.15628526976368214974D+00, &
     -0.79625257738750132575D+00, &
      0.41486074510619452838D+00, &
      0.91939856562697519493D+00, &
      0.41486074510619452838D+00, &
     -0.91939856562697519493D+00, &
     -0.41486074510619452838D+00, &
      0.91939856562697519493D+00, &
     -0.41486074510619452838D+00, &
     -0.91939856562697519493D+00, &
      0.68962930070823635909D+00, &
      0.51507330251609151350D+00, &
      0.68962930070823635909D+00, &
     -0.51507330251609151350D+00, &
     -0.68962930070823635909D+00, &
      0.51507330251609151350D+00, &
     -0.68962930070823635909D+00, &
     -0.51507330251609151350D+00, &
      0.85957679427925348659D+00, &
      0.96361213891168562284D+00, &
      0.85957679427925348659D+00, &
     -0.96361213891168562284D+00, &
     -0.85957679427925348659D+00, &
      0.96361213891168562284D+00, &
     -0.85957679427925348659D+00, &
     -0.96361213891168562284D+00, &
      0.49209543700208141503D+00, &
      0.72140711490936071382D+00, &
      0.49209543700208141503D+00, &
     -0.72140711490936071382D+00, &
     -0.49209543700208141503D+00, &
      0.72140711490936071382D+00, &
     -0.49209543700208141503D+00, &
     -0.72140711490936071382D+00 /)

  real ( kind = rk ) z(n)
  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &
      0.60412613281560589851D+00, &
      0.12832541213701875726D+00, &
      0.12832541213701875726D+00, &
      0.12832541213701875726D+00, &
      0.12832541213701875726D+00, &
      0.23636452793876766565D+00, &
      0.23636452793876766565D+00, &
      0.23636452793876766565D+00, &
      0.23636452793876766565D+00, &
      0.36948143968736257836D+00, &
      0.36948143968736257836D+00, &
      0.36948143968736257836D+00, &
      0.36948143968736257836D+00, &
      0.57413522292238361455D+00, &
      0.57413522292238361455D+00, &
      0.57413522292238361455D+00, &
      0.57413522292238361455D+00, &
      0.57097027691958568418D+00, &
      0.57097027691958568418D+00, &
      0.57097027691958568418D+00, &
      0.57097027691958568418D+00, &
      0.44532726334384131750D+00, &
      0.44532726334384131750D+00, &
      0.44532726334384131750D+00, &
      0.44532726334384131750D+00, &
      0.18507199264697138386D+00, &
      0.18507199264697138386D+00, &
      0.18507199264697138386D+00, &
      0.18507199264697138386D+00, &
      0.28281898722495468768D+00, &
      0.28281898722495468768D+00, &
      0.28281898722495468768D+00, &
      0.28281898722495468768D+00, &
      0.00745381355051214466D+00, &
      0.00745381355051214466D+00, &
      0.00745381355051214466D+00, &
      0.00745381355051214466D+00, &
      0.00273060888853154026D+00, &
      0.00273060888853154026D+00, &
      0.00273060888853154026D+00, &
      0.00273060888853154026D+00, &
      0.89474683073191085825D+00, &
      0.89474683073191085825D+00, &
      0.89474683073191085825D+00, &
      0.89474683073191085825D+00, &
      0.05489260378459437373D+00, &
      0.05489260378459437373D+00, &
      0.05489260378459437373D+00, &
      0.05489260378459437373D+00, &
      0.17890547617113325418D+00, &
      0.17890547617113325418D+00, &
      0.17890547617113325418D+00, &
      0.17890547617113325418D+00, &
      0.14438727929533640149D+00, &
      0.14438727929533640149D+00, &
      0.14438727929533640149D+00, &
      0.14438727929533640149D+00, &
      0.96005329029751285130D+00, &
      0.96005329029751285130D+00, &
      0.96005329029751285130D+00, &
      0.96005329029751285130D+00, &
      0.42262373956574911249D+00, &
      0.42262373956574911249D+00, &
      0.42262373956574911249D+00, &
      0.42262373956574911249D+00, &
      0.00937341178432237430D+00, &
      0.00937341178432237430D+00, &
      0.00937341178432237430D+00, &
      0.00937341178432237430D+00, &
      0.70511644808626805503D+00, &
      0.70511644808626805503D+00, &
      0.70511644808626805503D+00, &
      0.70511644808626805503D+00, &
      0.80976037637526798729D+00, &
      0.80976037637526798729D+00, &
      0.80976037637526798729D+00, &
      0.80976037637526798729D+00, &
      0.52107002257796464217D+00, &
      0.52107002257796464217D+00, &
      0.52107002257796464217D+00, &
      0.52107002257796464217D+00, &
      0.47622191262040014514D+00, &
      0.47622191262040014514D+00, &
      0.47622191262040014514D+00, &
      0.47622191262040014514D+00, &
      0.04553198949867700435D+00, &
      0.04553198949867700435D+00, &
      0.04553198949867700435D+00, &
      0.04553198949867700435D+00, &
      0.68236516778452971366D+00, &
      0.68236516778452971366D+00, &
      0.68236516778452971366D+00, &
      0.68236516778452971366D+00, &
      0.36766518307523687881D+00, &
      0.36766518307523687881D+00, &
      0.36766518307523687881D+00, &
      0.36766518307523687881D+00, &
      0.03140155892114703667D+00, &
      0.03140155892114703667D+00, &
      0.03140155892114703667D+00, &
      0.03140155892114703667D+00, &
      0.20854763179877680579D+00, &
      0.20854763179877680579D+00, &
      0.20854763179877680579D+00, &
      0.20854763179877680579D+00, &
      0.32095943995186076991D+00, &
      0.32095943995186076991D+00, &
      0.32095943995186076991D+00, &
      0.32095943995186076991D+00, &
      0.35590541054356356065D+00, &
      0.35590541054356356065D+00, &
      0.35590541054356356065D+00, &
      0.35590541054356356065D+00, &
      0.01963481046583391565D+00, &
      0.01963481046583391565D+00, &
      0.01963481046583391565D+00, &
      0.01963481046583391565D+00, &
      0.58110848325376152079D+00, &
      0.58110848325376152079D+00, &
      0.58110848325376152079D+00, &
      0.58110848325376152079D+00, &
      0.22030183551151610866D+00, &
      0.22030183551151610866D+00, &
      0.22030183551151610866D+00, &
      0.22030183551151610866D+00, &
      0.70743374610149956094D+00, &
      0.70743374610149956094D+00, &
      0.70743374610149956094D+00, &
      0.70743374610149956094D+00, &
      0.49754618392934674143D+00, &
      0.49754618392934674143D+00, &
      0.49754618392934674143D+00, &
      0.49754618392934674143D+00, &
      0.90197994879909670907D+00, &
      0.90197994879909670907D+00, &
      0.90197994879909670907D+00, &
      0.90197994879909670907D+00, &
      0.15830174812382535876D+00, &
      0.15830174812382535876D+00, &
      0.15830174812382535876D+00, &
      0.15830174812382535876D+00, &
      0.34703489296133355202D+00, &
      0.34703489296133355202D+00, &
      0.34703489296133355202D+00, &
      0.34703489296133355202D+00, &
      0.87727590777505048969D+00, &
      0.87727590777505048969D+00, &
      0.87727590777505048969D+00, &
      0.87727590777505048969D+00, &
      0.25098088788741812483D+00, &
      0.25098088788741812483D+00, &
      0.25098088788741812483D+00, &
      0.25098088788741812483D+00, &
      0.12287487061124809096D+00, &
      0.12287487061124809096D+00, &
      0.12287487061124809096D+00, &
      0.12287487061124809096D+00, &
      0.21569557215246823456D+00, &
      0.21569557215246823456D+00, &
      0.21569557215246823456D+00, &
      0.21569557215246823456D+00, &
      0.04316688220319717106D+00, &
      0.04316688220319717106D+00, &
      0.04316688220319717106D+00, &
      0.04316688220319717106D+00, &
      0.10724275430905566564D+00, &
      0.10724275430905566564D+00, &
      0.10724275430905566564D+00, &
      0.10724275430905566564D+00, &
      0.32964277824701476716D+00, &
      0.32964277824701476716D+00, &
      0.32964277824701476716D+00, &
      0.32964277824701476716D+00, &
      0.71600680775916891729D+00, &
      0.71600680775916891729D+00, &
      0.71600680775916891729D+00, &
      0.71600680775916891729D+00, &
      0.41973631072686923282D+00, &
      0.41973631072686923282D+00, &
      0.41973631072686923282D+00, &
      0.41973631072686923282D+00, &
      0.10771970940234577852D+00, &
      0.10771970940234577852D+00, &
      0.10771970940234577852D+00, &
      0.10771970940234577852D+00, &
      0.09137076099343852120D+00, &
      0.09137076099343852120D+00, &
      0.09137076099343852120D+00, &
      0.09137076099343852120D+00, &
      0.57433224555240058873D+00, &
      0.57433224555240058873D+00, &
      0.57433224555240058873D+00, &
      0.57433224555240058873D+00, &
      0.29652789568067272619D+00, &
      0.29652789568067272619D+00, &
      0.29652789568067272619D+00, &
      0.29652789568067272619D+00, &
      0.06861567128774997970D+00, &
      0.06861567128774997970D+00, &
      0.06861567128774997970D+00, &
      0.06861567128774997970D+00, &
      0.12112687504129497629D+00, &
      0.12112687504129497629D+00, &
      0.12112687504129497629D+00, &
      0.12112687504129497629D+00, &
      0.12112687504129497629D+00, &
      0.12112687504129497629D+00, &
      0.12112687504129497629D+00, &
      0.12112687504129497629D+00, &
      0.06232922095164156878D+00, &
      0.06232922095164156878D+00, &
      0.06232922095164156878D+00, &
      0.06232922095164156878D+00, &
      0.06232922095164156878D+00, &
      0.06232922095164156878D+00, &
      0.06232922095164156878D+00, &
      0.06232922095164156878D+00, &
      0.28878375602826145130D+00, &
      0.28878375602826145130D+00, &
      0.28878375602826145130D+00, &
      0.28878375602826145130D+00, &
      0.28878375602826145130D+00, &
      0.28878375602826145130D+00, &
      0.28878375602826145130D+00, &
      0.28878375602826145130D+00, &
      0.41981996327982368244D+00, &
      0.41981996327982368244D+00, &
      0.41981996327982368244D+00, &
      0.41981996327982368244D+00, &
      0.41981996327982368244D+00, &
      0.41981996327982368244D+00, &
      0.41981996327982368244D+00, &
      0.41981996327982368244D+00, &
      0.36719716335126328932D+00, &
      0.36719716335126328932D+00, &
      0.36719716335126328932D+00, &
      0.36719716335126328932D+00, &
      0.36719716335126328932D+00, &
      0.36719716335126328932D+00, &
      0.36719716335126328932D+00, &
      0.36719716335126328932D+00, &
      0.03539982573087990109D+00, &
      0.03539982573087990109D+00, &
      0.03539982573087990109D+00, &
      0.03539982573087990109D+00, &
      0.03539982573087990109D+00, &
      0.03539982573087990109D+00, &
      0.03539982573087990109D+00, &
      0.03539982573087990109D+00, &
      0.00689409549706292753D+00, &
      0.00689409549706292753D+00, &
      0.00689409549706292753D+00, &
      0.00689409549706292753D+00, &
      0.00689409549706292753D+00, &
      0.00689409549706292753D+00, &
      0.00689409549706292753D+00, &
      0.00689409549706292753D+00, &
      0.18034097822548403323D+00, &
      0.18034097822548403323D+00, &
      0.18034097822548403323D+00, &
      0.18034097822548403323D+00, &
      0.18034097822548403323D+00, &
      0.18034097822548403323D+00, &
      0.18034097822548403323D+00, &
      0.18034097822548403323D+00, &
      0.08701177937660883877D+00, &
      0.08701177937660883877D+00, &
      0.08701177937660883877D+00, &
      0.08701177937660883877D+00, &
      0.08701177937660883877D+00, &
      0.08701177937660883877D+00, &
      0.08701177937660883877D+00, &
      0.08701177937660883877D+00, &
      0.30184034189885056154D+00, &
      0.30184034189885056154D+00, &
      0.30184034189885056154D+00, &
      0.30184034189885056154D+00, &
      0.30184034189885056154D+00, &
      0.30184034189885056154D+00, &
      0.30184034189885056154D+00, &
      0.30184034189885056154D+00, &
      0.46375316476153072287D+00, &
      0.46375316476153072287D+00, &
      0.46375316476153072287D+00, &
      0.46375316476153072287D+00, &
      0.46375316476153072287D+00, &
      0.46375316476153072287D+00, &
      0.46375316476153072287D+00, &
      0.46375316476153072287D+00, &
      0.01259447070049873885D+00, &
      0.01259447070049873885D+00, &
      0.01259447070049873885D+00, &
      0.01259447070049873885D+00, &
      0.01259447070049873885D+00, &
      0.01259447070049873885D+00, &
      0.01259447070049873885D+00, &
      0.01259447070049873885D+00, &
      0.68464554842762070930D+00, &
      0.68464554842762070930D+00, &
      0.68464554842762070930D+00, &
      0.68464554842762070930D+00, &
      0.68464554842762070930D+00, &
      0.68464554842762070930D+00, &
      0.68464554842762070930D+00, &
      0.68464554842762070930D+00, &
      0.04669618994654996941D+00, &
      0.04669618994654996941D+00, &
      0.04669618994654996941D+00, &
      0.04669618994654996941D+00, &
      0.04669618994654996941D+00, &
      0.04669618994654996941D+00, &
      0.04669618994654996941D+00, &
      0.04669618994654996941D+00, &
      0.00997909979720432855D+00, &
      0.00997909979720432855D+00, &
      0.00997909979720432855D+00, &
      0.00997909979720432855D+00, &
      0.00997909979720432855D+00, &
      0.00997909979720432855D+00, &
      0.00997909979720432855D+00, &
      0.00997909979720432855D+00, &
      0.79170760938063011736D+00, &
      0.79170760938063011736D+00, &
      0.79170760938063011736D+00, &
      0.79170760938063011736D+00, &
      0.79170760938063011736D+00, &
      0.79170760938063011736D+00, &
      0.79170760938063011736D+00, &
      0.79170760938063011736D+00, &
      0.23239746807120234551D+00, &
      0.23239746807120234551D+00, &
      0.23239746807120234551D+00, &
      0.23239746807120234551D+00, &
      0.23239746807120234551D+00, &
      0.23239746807120234551D+00, &
      0.23239746807120234551D+00, &
      0.23239746807120234551D+00, &
      0.06212941698545271230D+00, &
      0.06212941698545271230D+00, &
      0.06212941698545271230D+00, &
      0.06212941698545271230D+00, &
      0.06212941698545271230D+00, &
      0.06212941698545271230D+00, &
      0.06212941698545271230D+00, &
      0.06212941698545271230D+00, &
      0.10250007471008426574D+00, &
      0.10250007471008426574D+00, &
      0.10250007471008426574D+00, &
      0.10250007471008426574D+00, &
      0.10250007471008426574D+00, &
      0.10250007471008426574D+00, &
      0.10250007471008426574D+00, &
      0.10250007471008426574D+00, &
      0.19638105364011448906D+00, &
      0.19638105364011448906D+00, &
      0.19638105364011448906D+00, &
      0.19638105364011448906D+00, &
      0.19638105364011448906D+00, &
      0.19638105364011448906D+00, &
      0.19638105364011448906D+00, &
      0.19638105364011448906D+00, &
      0.03807225396126873163D+00, &
      0.03807225396126873163D+00, &
      0.03807225396126873163D+00, &
      0.03807225396126873163D+00, &
      0.03807225396126873163D+00, &
      0.03807225396126873163D+00, &
      0.03807225396126873163D+00, &
      0.03807225396126873163D+00, &
      0.00969993829969578725D+00, &
      0.00969993829969578725D+00, &
      0.00969993829969578725D+00, &
      0.00969993829969578725D+00, &
      0.00969993829969578725D+00, &
      0.00969993829969578725D+00, &
      0.00969993829969578725D+00, &
      0.00969993829969578725D+00, &
      0.08712730173221996943D+00, &
      0.08712730173221996943D+00, &
      0.08712730173221996943D+00, &
      0.08712730173221996943D+00, &
      0.08712730173221996943D+00, &
      0.08712730173221996943D+00, &
      0.08712730173221996943D+00, &
      0.08712730173221996943D+00, &
      0.04169093747605617795D+00, &
      0.04169093747605617795D+00, &
      0.04169093747605617795D+00, &
      0.04169093747605617795D+00, &
      0.04169093747605617795D+00, &
      0.04169093747605617795D+00, &
      0.04169093747605617795D+00, &
      0.04169093747605617795D+00, &
      0.18982222362689099571D+00, &
      0.18982222362689099571D+00, &
      0.18982222362689099571D+00, &
      0.18982222362689099571D+00, &
      0.18982222362689099571D+00, &
      0.18982222362689099571D+00, &
      0.18982222362689099571D+00, &
      0.18982222362689099571D+00, &
      0.12762173979455287975D+00, &
      0.12762173979455287975D+00, &
      0.12762173979455287975D+00, &
      0.12762173979455287975D+00, &
      0.12762173979455287975D+00, &
      0.12762173979455287975D+00, &
      0.12762173979455287975D+00, &
      0.12762173979455287975D+00, &
      0.19550632719707150553D+00, &
      0.19550632719707150553D+00, &
      0.19550632719707150553D+00, &
      0.19550632719707150553D+00, &
      0.19550632719707150553D+00, &
      0.19550632719707150553D+00, &
      0.19550632719707150553D+00, &
      0.19550632719707150553D+00, &
      0.49356773754307026181D+00, &
      0.49356773754307026181D+00, &
      0.49356773754307026181D+00, &
      0.49356773754307026181D+00, &
      0.49356773754307026181D+00, &
      0.49356773754307026181D+00, &
      0.49356773754307026181D+00, &
      0.49356773754307026181D+00, &
      0.03472481178210238412D+00, &
      0.03472481178210238412D+00, &
      0.03472481178210238412D+00, &
      0.03472481178210238412D+00, &
      0.03472481178210238412D+00, &
      0.03472481178210238412D+00, &
      0.03472481178210238412D+00, &
      0.03472481178210238412D+00, &
      0.01518094278770998860D+00, &
      0.01518094278770998860D+00, &
      0.01518094278770998860D+00, &
      0.01518094278770998860D+00, &
      0.01518094278770998860D+00, &
      0.01518094278770998860D+00, &
      0.01518094278770998860D+00, &
      0.01518094278770998860D+00, &
      0.61542031522620999073D+00, &
      0.61542031522620999073D+00, &
      0.61542031522620999073D+00, &
      0.61542031522620999073D+00, &
      0.61542031522620999073D+00, &
      0.61542031522620999073D+00, &
      0.61542031522620999073D+00, &
      0.61542031522620999073D+00, &
      0.12575204752531654595D+00, &
      0.12575204752531654595D+00, &
      0.12575204752531654595D+00, &
      0.12575204752531654595D+00, &
      0.12575204752531654595D+00, &
      0.12575204752531654595D+00, &
      0.12575204752531654595D+00, &
      0.12575204752531654595D+00, &
      0.04457888016467764086D+00, &
      0.04457888016467764086D+00, &
      0.04457888016467764086D+00, &
      0.04457888016467764086D+00, &
      0.04457888016467764086D+00, &
      0.04457888016467764086D+00, &
      0.04457888016467764086D+00, &
      0.04457888016467764086D+00, &
      0.30170781950633751567D+00, &
      0.30170781950633751567D+00, &
      0.30170781950633751567D+00, &
      0.30170781950633751567D+00, &
      0.30170781950633751567D+00, &
      0.30170781950633751567D+00, &
      0.30170781950633751567D+00, &
      0.30170781950633751567D+00, &
      0.00577128438055363874D+00, &
      0.00577128438055363874D+00, &
      0.00577128438055363874D+00, &
      0.00577128438055363874D+00, &
      0.00577128438055363874D+00, &
      0.00577128438055363874D+00, &
      0.00577128438055363874D+00, &
      0.00577128438055363874D+00, &
      0.00239530919218184334D+00, &
      0.00239530919218184334D+00, &
      0.00239530919218184334D+00, &
      0.00239530919218184334D+00, &
      0.00239530919218184334D+00, &
      0.00239530919218184334D+00, &
      0.00239530919218184334D+00, &
      0.00239530919218184334D+00 /)

  real ( kind = rk ) w(n)
  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &
      0.00537430973147919565D+00, &
      0.00567173429377825552D+00, &
      0.00567173429377825552D+00, &
      0.00567173429377825552D+00, &
      0.00567173429377825552D+00, &
      0.00324886088172405109D+00, &
      0.00324886088172405109D+00, &
      0.00324886088172405109D+00, &
      0.00324886088172405109D+00, &
      0.00345209452668273666D+00, &
      0.00345209452668273666D+00, &
      0.00345209452668273666D+00, &
      0.00345209452668273666D+00, &
      0.00066660334447695205D+00, &
      0.00066660334447695205D+00, &
      0.00066660334447695205D+00, &
      0.00066660334447695205D+00, &
      0.00399579246038160064D+00, &
      0.00399579246038160064D+00, &
      0.00399579246038160064D+00, &
      0.00399579246038160064D+00, &
      0.00499277895953007846D+00, &
      0.00499277895953007846D+00, &
      0.00499277895953007846D+00, &
      0.00499277895953007846D+00, &
      0.00246962297642099688D+00, &
      0.00246962297642099688D+00, &
      0.00246962297642099688D+00, &
      0.00246962297642099688D+00, &
      0.00706615575413061577D+00, &
      0.00706615575413061577D+00, &
      0.00706615575413061577D+00, &
      0.00706615575413061577D+00, &
      0.00191698003537197552D+00, &
      0.00191698003537197552D+00, &
      0.00191698003537197552D+00, &
      0.00191698003537197552D+00, &
      0.00023875415581936437D+00, &
      0.00023875415581936437D+00, &
      0.00023875415581936437D+00, &
      0.00023875415581936437D+00, &
      0.00016802396701300472D+00, &
      0.00016802396701300472D+00, &
      0.00016802396701300472D+00, &
      0.00016802396701300472D+00, &
      0.00099443273732332323D+00, &
      0.00099443273732332323D+00, &
      0.00099443273732332323D+00, &
      0.00099443273732332323D+00, &
      0.00345494167713882722D+00, &
      0.00345494167713882722D+00, &
      0.00345494167713882722D+00, &
      0.00345494167713882722D+00, &
      0.00137778065777168987D+00, &
      0.00137778065777168987D+00, &
      0.00137778065777168987D+00, &
      0.00137778065777168987D+00, &
      0.00006135678341459227D+00, &
      0.00006135678341459227D+00, &
      0.00006135678341459227D+00, &
      0.00006135678341459227D+00, &
      0.00099615960761120759D+00, &
      0.00099615960761120759D+00, &
      0.00099615960761120759D+00, &
      0.00099615960761120759D+00, &
      0.00173468884468413759D+00, &
      0.00173468884468413759D+00, &
      0.00173468884468413759D+00, &
      0.00173468884468413759D+00, &
      0.00064089052901790020D+00, &
      0.00064089052901790020D+00, &
      0.00064089052901790020D+00, &
      0.00064089052901790020D+00, &
      0.00174770174992903656D+00, &
      0.00174770174992903656D+00, &
      0.00174770174992903656D+00, &
      0.00174770174992903656D+00, &
      0.00053561830982834806D+00, &
      0.00053561830982834806D+00, &
      0.00053561830982834806D+00, &
      0.00053561830982834806D+00, &
      0.00374412815566966747D+00, &
      0.00374412815566966747D+00, &
      0.00374412815566966747D+00, &
      0.00374412815566966747D+00, &
      0.00364849854393893731D+00, &
      0.00364849854393893731D+00, &
      0.00364849854393893731D+00, &
      0.00364849854393893731D+00, &
      0.00027221046304926988D+00, &
      0.00027221046304926988D+00, &
      0.00027221046304926988D+00, &
      0.00027221046304926988D+00, &
      0.00403805862793428852D+00, &
      0.00403805862793428852D+00, &
      0.00403805862793428852D+00, &
      0.00403805862793428852D+00, &
      0.00284945501315263380D+00, &
      0.00284945501315263380D+00, &
      0.00284945501315263380D+00, &
      0.00284945501315263380D+00, &
      0.00049463185151189392D+00, &
      0.00049463185151189392D+00, &
      0.00049463185151189392D+00, &
      0.00049463185151189392D+00, &
      0.00343575024351558140D+00, &
      0.00343575024351558140D+00, &
      0.00343575024351558140D+00, &
      0.00343575024351558140D+00, &
      0.00076547936426975662D+00, &
      0.00076547936426975662D+00, &
      0.00076547936426975662D+00, &
      0.00076547936426975662D+00, &
      0.00009270297025443059D+00, &
      0.00009270297025443059D+00, &
      0.00009270297025443059D+00, &
      0.00009270297025443059D+00, &
      0.00444044740059815115D+00, &
      0.00444044740059815115D+00, &
      0.00444044740059815115D+00, &
      0.00444044740059815115D+00, &
      0.00499528418027490092D+00, &
      0.00499528418027490092D+00, &
      0.00499528418027490092D+00, &
      0.00499528418027490092D+00, &
      0.00232760564344115621D+00, &
      0.00232760564344115621D+00, &
      0.00232760564344115621D+00, &
      0.00232760564344115621D+00, &
      0.00305846111057695401D+00, &
      0.00305846111057695401D+00, &
      0.00305846111057695401D+00, &
      0.00305846111057695401D+00, &
      0.00027167033342833963D+00, &
      0.00027167033342833963D+00, &
      0.00027167033342833963D+00, &
      0.00027167033342833963D+00, &
      0.00359204041455288758D+00, &
      0.00359204041455288758D+00, &
      0.00359204041455288758D+00, &
      0.00359204041455288758D+00, &
      0.00450623232462683300D+00, &
      0.00450623232462683300D+00, &
      0.00450623232462683300D+00, &
      0.00450623232462683300D+00, &
      0.00026620595222711166D+00, &
      0.00026620595222711166D+00, &
      0.00026620595222711166D+00, &
      0.00026620595222711166D+00, &
      0.00172242225189580926D+00, &
      0.00172242225189580926D+00, &
      0.00172242225189580926D+00, &
      0.00172242225189580926D+00, &
      0.00152990216308249529D+00, &
      0.00152990216308249529D+00, &
      0.00152990216308249529D+00, &
      0.00152990216308249529D+00, &
      0.00433172345565203346D+00, &
      0.00433172345565203346D+00, &
      0.00433172345565203346D+00, &
      0.00433172345565203346D+00, &
      0.00086984955905943873D+00, &
      0.00086984955905943873D+00, &
      0.00086984955905943873D+00, &
      0.00086984955905943873D+00, &
      0.00316655761594148268D+00, &
      0.00316655761594148268D+00, &
      0.00316655761594148268D+00, &
      0.00316655761594148268D+00, &
      0.00145516775846500309D+00, &
      0.00145516775846500309D+00, &
      0.00145516775846500309D+00, &
      0.00145516775846500309D+00, &
      0.00113761603493167554D+00, &
      0.00113761603493167554D+00, &
      0.00113761603493167554D+00, &
      0.00113761603493167554D+00, &
      0.00188490438846629360D+00, &
      0.00188490438846629360D+00, &
      0.00188490438846629360D+00, &
      0.00188490438846629360D+00, &
      0.00604980306252401318D+00, &
      0.00604980306252401318D+00, &
      0.00604980306252401318D+00, &
      0.00604980306252401318D+00, &
      0.00027850531603908424D+00, &
      0.00027850531603908424D+00, &
      0.00027850531603908424D+00, &
      0.00027850531603908424D+00, &
      0.00164429091495763209D+00, &
      0.00164429091495763209D+00, &
      0.00164429091495763209D+00, &
      0.00164429091495763209D+00, &
      0.00438466987463943282D+00, &
      0.00438466987463943282D+00, &
      0.00438466987463943282D+00, &
      0.00438466987463943282D+00, &
      0.00250821736032557372D+00, &
      0.00250821736032557372D+00, &
      0.00250821736032557372D+00, &
      0.00250821736032557372D+00, &
      0.00064031929971637594D+00, &
      0.00064031929971637594D+00, &
      0.00064031929971637594D+00, &
      0.00064031929971637594D+00, &
      0.00064031929971637594D+00, &
      0.00064031929971637594D+00, &
      0.00064031929971637594D+00, &
      0.00064031929971637594D+00, &
      0.00274360981913367684D+00, &
      0.00274360981913367684D+00, &
      0.00274360981913367684D+00, &
      0.00274360981913367684D+00, &
      0.00274360981913367684D+00, &
      0.00274360981913367684D+00, &
      0.00274360981913367684D+00, &
      0.00274360981913367684D+00, &
      0.00442513215733170758D+00, &
      0.00442513215733170758D+00, &
      0.00442513215733170758D+00, &
      0.00442513215733170758D+00, &
      0.00442513215733170758D+00, &
      0.00442513215733170758D+00, &
      0.00442513215733170758D+00, &
      0.00442513215733170758D+00, &
      0.00379353354761960502D+00, &
      0.00379353354761960502D+00, &
      0.00379353354761960502D+00, &
      0.00379353354761960502D+00, &
      0.00379353354761960502D+00, &
      0.00379353354761960502D+00, &
      0.00379353354761960502D+00, &
      0.00379353354761960502D+00, &
      0.00259873188089142833D+00, &
      0.00259873188089142833D+00, &
      0.00259873188089142833D+00, &
      0.00259873188089142833D+00, &
      0.00259873188089142833D+00, &
      0.00259873188089142833D+00, &
      0.00259873188089142833D+00, &
      0.00259873188089142833D+00, &
      0.00110056207140428312D+00, &
      0.00110056207140428312D+00, &
      0.00110056207140428312D+00, &
      0.00110056207140428312D+00, &
      0.00110056207140428312D+00, &
      0.00110056207140428312D+00, &
      0.00110056207140428312D+00, &
      0.00110056207140428312D+00, &
      0.00075594512745462764D+00, &
      0.00075594512745462764D+00, &
      0.00075594512745462764D+00, &
      0.00075594512745462764D+00, &
      0.00075594512745462764D+00, &
      0.00075594512745462764D+00, &
      0.00075594512745462764D+00, &
      0.00075594512745462764D+00, &
      0.00416443295717523060D+00, &
      0.00416443295717523060D+00, &
      0.00416443295717523060D+00, &
      0.00416443295717523060D+00, &
      0.00416443295717523060D+00, &
      0.00416443295717523060D+00, &
      0.00416443295717523060D+00, &
      0.00416443295717523060D+00, &
      0.00079864720063192577D+00, &
      0.00079864720063192577D+00, &
      0.00079864720063192577D+00, &
      0.00079864720063192577D+00, &
      0.00079864720063192577D+00, &
      0.00079864720063192577D+00, &
      0.00079864720063192577D+00, &
      0.00079864720063192577D+00, &
      0.00127477283170923661D+00, &
      0.00127477283170923661D+00, &
      0.00127477283170923661D+00, &
      0.00127477283170923661D+00, &
      0.00127477283170923661D+00, &
      0.00127477283170923661D+00, &
      0.00127477283170923661D+00, &
      0.00127477283170923661D+00, &
      0.00077420720151948649D+00, &
      0.00077420720151948649D+00, &
      0.00077420720151948649D+00, &
      0.00077420720151948649D+00, &
      0.00077420720151948649D+00, &
      0.00077420720151948649D+00, &
      0.00077420720151948649D+00, &
      0.00077420720151948649D+00, &
      0.00107986786877537822D+00, &
      0.00107986786877537822D+00, &
      0.00107986786877537822D+00, &
      0.00107986786877537822D+00, &
      0.00107986786877537822D+00, &
      0.00107986786877537822D+00, &
      0.00107986786877537822D+00, &
      0.00107986786877537822D+00, &
      0.00155093578240393421D+00, &
      0.00155093578240393421D+00, &
      0.00155093578240393421D+00, &
      0.00155093578240393421D+00, &
      0.00155093578240393421D+00, &
      0.00155093578240393421D+00, &
      0.00155093578240393421D+00, &
      0.00155093578240393421D+00, &
      0.00049937862088166734D+00, &
      0.00049937862088166734D+00, &
      0.00049937862088166734D+00, &
      0.00049937862088166734D+00, &
      0.00049937862088166734D+00, &
      0.00049937862088166734D+00, &
      0.00049937862088166734D+00, &
      0.00049937862088166734D+00, &
      0.00042623664206464421D+00, &
      0.00042623664206464421D+00, &
      0.00042623664206464421D+00, &
      0.00042623664206464421D+00, &
      0.00042623664206464421D+00, &
      0.00042623664206464421D+00, &
      0.00042623664206464421D+00, &
      0.00042623664206464421D+00, &
      0.00052445285235233835D+00, &
      0.00052445285235233835D+00, &
      0.00052445285235233835D+00, &
      0.00052445285235233835D+00, &
      0.00052445285235233835D+00, &
      0.00052445285235233835D+00, &
      0.00052445285235233835D+00, &
      0.00052445285235233835D+00, &
      0.00295483873714207741D+00, &
      0.00295483873714207741D+00, &
      0.00295483873714207741D+00, &
      0.00295483873714207741D+00, &
      0.00295483873714207741D+00, &
      0.00295483873714207741D+00, &
      0.00295483873714207741D+00, &
      0.00295483873714207741D+00, &
      0.00154064844826157059D+00, &
      0.00154064844826157059D+00, &
      0.00154064844826157059D+00, &
      0.00154064844826157059D+00, &
      0.00154064844826157059D+00, &
      0.00154064844826157059D+00, &
      0.00154064844826157059D+00, &
      0.00154064844826157059D+00, &
      0.00484986123625866239D+00, &
      0.00484986123625866239D+00, &
      0.00484986123625866239D+00, &
      0.00484986123625866239D+00, &
      0.00484986123625866239D+00, &
      0.00484986123625866239D+00, &
      0.00484986123625866239D+00, &
      0.00484986123625866239D+00, &
      0.00074953722109848120D+00, &
      0.00074953722109848120D+00, &
      0.00074953722109848120D+00, &
      0.00074953722109848120D+00, &
      0.00074953722109848120D+00, &
      0.00074953722109848120D+00, &
      0.00074953722109848120D+00, &
      0.00074953722109848120D+00, &
      0.00167732903590775981D+00, &
      0.00167732903590775981D+00, &
      0.00167732903590775981D+00, &
      0.00167732903590775981D+00, &
      0.00167732903590775981D+00, &
      0.00167732903590775981D+00, &
      0.00167732903590775981D+00, &
      0.00167732903590775981D+00, &
      0.00224204091471395155D+00, &
      0.00224204091471395155D+00, &
      0.00224204091471395155D+00, &
      0.00224204091471395155D+00, &
      0.00224204091471395155D+00, &
      0.00224204091471395155D+00, &
      0.00224204091471395155D+00, &
      0.00224204091471395155D+00, &
      0.00305423895198708934D+00, &
      0.00305423895198708934D+00, &
      0.00305423895198708934D+00, &
      0.00305423895198708934D+00, &
      0.00305423895198708934D+00, &
      0.00305423895198708934D+00, &
      0.00305423895198708934D+00, &
      0.00305423895198708934D+00, &
      0.00236370915044872495D+00, &
      0.00236370915044872495D+00, &
      0.00236370915044872495D+00, &
      0.00236370915044872495D+00, &
      0.00236370915044872495D+00, &
      0.00236370915044872495D+00, &
      0.00236370915044872495D+00, &
      0.00236370915044872495D+00, &
      0.00412477376791669549D+00, &
      0.00412477376791669549D+00, &
      0.00412477376791669549D+00, &
      0.00412477376791669549D+00, &
      0.00412477376791669549D+00, &
      0.00412477376791669549D+00, &
      0.00412477376791669549D+00, &
      0.00412477376791669549D+00, &
      0.00252281437387082835D+00, &
      0.00252281437387082835D+00, &
      0.00252281437387082835D+00, &
      0.00252281437387082835D+00, &
      0.00252281437387082835D+00, &
      0.00252281437387082835D+00, &
      0.00252281437387082835D+00, &
      0.00252281437387082835D+00, &
      0.00114336403720023422D+00, &
      0.00114336403720023422D+00, &
      0.00114336403720023422D+00, &
      0.00114336403720023422D+00, &
      0.00114336403720023422D+00, &
      0.00114336403720023422D+00, &
      0.00114336403720023422D+00, &
      0.00114336403720023422D+00, &
      0.00278500419257705868D+00, &
      0.00278500419257705868D+00, &
      0.00278500419257705868D+00, &
      0.00278500419257705868D+00, &
      0.00278500419257705868D+00, &
      0.00278500419257705868D+00, &
      0.00278500419257705868D+00, &
      0.00278500419257705868D+00, &
      0.00222262578256617873D+00, &
      0.00222262578256617873D+00, &
      0.00222262578256617873D+00, &
      0.00222262578256617873D+00, &
      0.00222262578256617873D+00, &
      0.00222262578256617873D+00, &
      0.00222262578256617873D+00, &
      0.00222262578256617873D+00, &
      0.00026099460442016807D+00, &
      0.00026099460442016807D+00, &
      0.00026099460442016807D+00, &
      0.00026099460442016807D+00, &
      0.00026099460442016807D+00, &
      0.00026099460442016807D+00, &
      0.00026099460442016807D+00, &
      0.00026099460442016807D+00, &
      0.00125217421068957169D+00, &
      0.00125217421068957169D+00, &
      0.00125217421068957169D+00, &
      0.00125217421068957169D+00, &
      0.00125217421068957169D+00, &
      0.00125217421068957169D+00, &
      0.00125217421068957169D+00, &
      0.00125217421068957169D+00, &
      0.00127967033864278819D+00, &
      0.00127967033864278819D+00, &
      0.00127967033864278819D+00, &
      0.00127967033864278819D+00, &
      0.00127967033864278819D+00, &
      0.00127967033864278819D+00, &
      0.00127967033864278819D+00, &
      0.00127967033864278819D+00, &
      0.00079904956567299025D+00, &
      0.00079904956567299025D+00, &
      0.00079904956567299025D+00, &
      0.00079904956567299025D+00, &
      0.00079904956567299025D+00, &
      0.00079904956567299025D+00, &
      0.00079904956567299025D+00, &
      0.00079904956567299025D+00, &
      0.00078122198597418405D+00, &
      0.00078122198597418405D+00, &
      0.00078122198597418405D+00, &
      0.00078122198597418405D+00, &
      0.00078122198597418405D+00, &
      0.00078122198597418405D+00, &
      0.00078122198597418405D+00, &
      0.00078122198597418405D+00, &
      0.00021725880518099310D+00, &
      0.00021725880518099310D+00, &
      0.00021725880518099310D+00, &
      0.00021725880518099310D+00, &
      0.00021725880518099310D+00, &
      0.00021725880518099310D+00, &
      0.00021725880518099310D+00, &
      0.00021725880518099310D+00, &
      0.00075955776144381827D+00, &
      0.00075955776144381827D+00, &
      0.00075955776144381827D+00, &
      0.00075955776144381827D+00, &
      0.00075955776144381827D+00, &
      0.00075955776144381827D+00, &
      0.00075955776144381827D+00, &
      0.00075955776144381827D+00 /)

  x(1:n) = x_save(1:n)
  y(1:n) = y_save(1:n)
  z(1:n) = z_save(1:n)
  w(1:n) = w_save(1:n)

  return
end
