function binom ( n, k )

!*****************************************************************************80
!
!! binom() computes the binomial coefficient.
!
!  Discussion:
!
!    This is ACM algorithm 160 translated to Fortran.
!
!    It calculates the number of combinations of N things taken K at a time.
!
!  Modified:
!
!    28 March 2016
!
!  Author:
!
!    Bill Buckles, Matthew Lybanon
!
!  Reference:
!
!    Bill Buckles, Matthew Lybanon,
!    Algorithm 515: Generation of a Vector from the Lexicographical Index,
!    ACM Transactions on Mathematical Software,
!    Volume 3, Number 2, June 1977, pages 180-182.
!
!  Parameters:
!
!    Input, integer N, K, the parameters for the binomial 
!    coefficient.
!
!    Output, integer BINOM, the binomial coefficient.
!
  implicit none

  integer binom
  integer i
  integer k
  integer k1
  integer n
  integer p
  integer r

  k1 = k
  p = n - k1

  if ( k1 < p ) then
    p = k1
    k1 = n - p
  end if

  if ( p == 0 ) then
    r = 1
  else
    r = k1 + 1
  end if

  do i = 2, p
    r = ( r * ( k1 + i ) ) / i
  end do

  binom = r

  return
end
subroutine comb ( n, p, l, c )

!*****************************************************************************80
!
!! COMB selects a subset of order P from a set of order N.
!
!  Discussion:
!
!    This subroutine finds the combination set of N things taken
!    P at a time for a given lexicographic index.
!
!  Modified:
!
!    27 March 2016
!
!  Author:
!
!    Bill Buckles, Matthew Lybanon
!
!  Reference:
!
!    Bill Buckles, Matthew Lybanon,
!    Algorithm 515: Generation of a Vector from the Lexicographical Index,
!    ACM Transactions on Mathematical Software,
!    Volume 3, Number 2, June 1977, pages 180-182.
!
!  Parameters:
!
!    Input, integer N, the number of things in the set.
!
!    Input, integer P, the number of things in each combination.
!    0 < P < N.
!
!    Input, integer L, the lexicographic index of the 
!    desired combination.  1 <= L <= choose(N,P).
!
!    Output, integer C(P), the combination set.
!
  implicit none

  integer p

  integer binom
  integer c(p)
  integer i
  integer k
  integer l
  integer n
  integer p1
  integer r
!
!  Special case: P = 1
!
  if ( p == 1 ) then
    c(1) = l
    return
  end if
!
!  Initialize lower bound index.
!
  k = 0
!
!  Select elements in ascending order.
!
  p1 = p - 1
  c(1) = 0

  do i = 1, p1
!
!  Update lower bound as the previously selected element.
!
    if ( 1 < i ) then
      c(i) = c(i-1)
    end if
!
!  Check validity of each entry.
!
    do

      c(i) = c(i) + 1
      r = binom ( n - c(i), p - i )
      k = k + r

      if ( l <= k ) then
        exit
      end if

    end do

    k = k - r

  end do

  c(p) = c(p1) + l - k

  return
end
function i4_choose_check ( n, k )

!*****************************************************************************80
!
!! I4_CHOOSE_CHECK reports whether the binomial coefficient can be computed.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 March 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, K, the binomial parameters.
!
!    Output, logical I4_CHOOSE_CHECK is:
!    TRUE, if C(N,K) < maximum integer.
!    FALSE, otherwise.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  logical i4_choose_check
  real ( kind = 8 ) choose_nk_log
  integer, parameter :: i4_huge = 2147483647
  real ( kind = 8 ) i4_huge_log
  integer n
  integer k
  real ( kind = 8 ) r8_gamma_log

  i4_huge_log = log ( dble ( i4_huge ) )

  choose_nk_log = &
      r8_gamma_log ( real ( n + 1, kind = 8 ) ) &
    - r8_gamma_log ( real ( k + 1, kind = 8 ) ) &
    - r8_gamma_log ( real ( n - k + 1, kind = 8 ) )

  i4_choose_check = ( choose_nk_log < i4_huge_log )
        
  return
end
function i4_uniform_ab ( a, b, seed )

!*****************************************************************************80
!
!! I4_UNIFORM_AB returns a scaled pseudorandom I4 between A and B.
!
!  Discussion:
!
!    An I4 is an integer value.
!
!    The pseudorandom number will be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 October 2012
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
!    Input, integer A, B, the limits of the interval.
!
!    Input/output, integer SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer I4_UNIFORM_AB, a number between A and B.
!
  implicit none

  integer a
  integer b
  integer, parameter :: i4_huge = 2147483647
  integer i4_uniform_ab
  integer k
  real r
  integer seed
  integer value

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_UNIFORM_AB - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r = real ( seed ) * 4.656612875E-10
!
!  Scale R to lie between A-0.5 and B+0.5.
!
  r = ( 1.0E+00 - r ) * ( real ( min ( a, b ) ) - 0.5E+00 ) & 
    +             r   * ( real ( max ( a, b ) ) + 0.5E+00 )
!
!  Use rounding to convert R to an integer between A and B.
!
  value = nint ( r )

  value = max ( value, min ( a, b ) )
  value = min ( value, max ( a, b ) )

  i4_uniform_ab = value

  return
end
function r8_gamma_log ( x )

!*****************************************************************************80
!
!! R8_GAMMA_LOG evaluates the logarithm of the gamma function.
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
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) R8_GAMMA_LOG, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = 8 ), dimension ( 7 ) :: c = (/ &
    -1.910444077728D-03, &
     8.4171387781295D-04, &
    -5.952379913043012D-04, &
     7.93650793500350248D-04, &
    -2.777777777777681622553D-03, &
     8.333333333333333331554247D-02, &
     5.7083835261D-03 /)
  real ( kind = 8 ) corr
  real ( kind = 8 ) :: d1 = -5.772156649015328605195174D-01
  real ( kind = 8 ) :: d2 = 4.227843350984671393993777D-01
  real ( kind = 8 ) :: d4 = 1.791759469228055000094023D+00
  real ( kind = 8 ), parameter :: frtbig = 2.25D+76
  integer i
  real ( kind = 8 ), dimension ( 8 ) :: p1 = (/ &
    4.945235359296727046734888D+00, &
    2.018112620856775083915565D+02, &
    2.290838373831346393026739D+03, &
    1.131967205903380828685045D+04, &
    2.855724635671635335736389D+04, &
    3.848496228443793359990269D+04, &
    2.637748787624195437963534D+04, &
    7.225813979700288197698961D+03 /)
  real ( kind = 8 ), dimension ( 8 ) :: p2 = (/ &
    4.974607845568932035012064D+00, &
    5.424138599891070494101986D+02, &
    1.550693864978364947665077D+04, &
    1.847932904445632425417223D+05, &
    1.088204769468828767498470D+06, &
    3.338152967987029735917223D+06, &
    5.106661678927352456275255D+06, &
    3.074109054850539556250927D+06 /)
  real ( kind = 8 ), dimension ( 8 ) :: p4 = (/ &
    1.474502166059939948905062D+04, &
    2.426813369486704502836312D+06, &
    1.214755574045093227939592D+08, &
    2.663432449630976949898078D+09, &
    2.940378956634553899906876D+10, &
    1.702665737765398868392998D+11, &
    4.926125793377430887588120D+11, &
    5.606251856223951465078242D+11 /)
  real ( kind = 8 ), dimension ( 8 ) :: q1 = (/ &
    6.748212550303777196073036D+01, &
    1.113332393857199323513008D+03, &
    7.738757056935398733233834D+03, &
    2.763987074403340708898585D+04, &
    5.499310206226157329794414D+04, &
    6.161122180066002127833352D+04, &
    3.635127591501940507276287D+04, &
    8.785536302431013170870835D+03 /)
  real ( kind = 8 ), dimension ( 8 ) :: q2 = (/ &
    1.830328399370592604055942D+02, &
    7.765049321445005871323047D+03, &
    1.331903827966074194402448D+05, &
    1.136705821321969608938755D+06, &
    5.267964117437946917577538D+06, &
    1.346701454311101692290052D+07, &
    1.782736530353274213975932D+07, &
    9.533095591844353613395747D+06 /)
  real ( kind = 8 ), dimension ( 8 ) :: q4 = (/ &
    2.690530175870899333379843D+03, &
    6.393885654300092398984238D+05, &
    4.135599930241388052042842D+07, &
    1.120872109616147941376570D+09, &
    1.488613728678813811542398D+10, &
    1.016803586272438228077304D+11, &
    3.417476345507377132798597D+11, &
    4.463158187419713286462081D+11 /)
  real ( kind = 8 ) r8_gamma_log
  real ( kind = 8 ) res
  real ( kind = 8 ), parameter :: sqrtpi = 0.9189385332046727417803297D+00
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter :: xbig = 2.55D+305
  real ( kind = 8 ) xden
  real ( kind = 8 ), parameter :: xinf = 1.79D+308
  real ( kind = 8 ) xm1
  real ( kind = 8 ) xm2
  real ( kind = 8 ) xm4
  real ( kind = 8 ) xnum
  real ( kind = 8 ) y
  real ( kind = 8 ) ysq

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
