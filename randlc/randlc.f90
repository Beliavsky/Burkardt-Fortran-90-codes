function randlc ( x )

!*****************************************************************************80
!
!! randlc() returns a uniform pseudorandom value.
!
!  Discussion:
!
!    The number returned is in the range (0, 1).  
!
!    The algorithm uses the linear congruential generator:
!
!      X(K+1) = A * X(K)  mod 2^46
!
!    This scheme generates 2^44 numbers before repeating.  
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Bailey, Eric Barszcz, John Barton, D Browning, Robert Carter, 
!    Leonardo Dagum, Rod Fatoohi,
!    Samuel Fineberg, Paul Frederickson, Thomas Lasinski, Robert Schreiber, 
!    Horst Simon, V Venkatakrishnan, Sisira Weeratunga,
!    The NAS Parallel Benchmarks,
!    RNR Technical Report RNR-94-007,
!    March 1994.
!
!    Donald Knuth,
!    The Art of Computer Programming,
!    Volume 2, Seminumerical Algorithms,
!    Third Edition,
!    Addison Wesley, 1997,
!    ISBN: 0201896842,
!    LC: QA76.6.K64.
!
!  Parameters:
!
!    Input/output, real ( kind = rk8 ) X, the seed.  X should be an 
!    odd integer such that 1 <= X <= 2^46.
!
!    Output, real ( kind = rk8 ) RANDLC, the next pseudorandom value.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ), parameter :: a = 1220703125.0D+00
  real ( kind = rk8 ), save :: a1
  real ( kind = rk8 ), save :: a2
  integer i
  integer, save :: ks = 0
  real ( kind = rk8 ), save :: r23
  real ( kind = rk8 ), save :: r46
  real ( kind = rk8 ) randlc
  real ( kind = rk8 ) t1
  real ( kind = rk8 ) t2
  real ( kind = rk8 ), save :: t23
  real ( kind = rk8 ) t3
  real ( kind = rk8 ) t4
  real ( kind = rk8 ), save :: t46
  real ( kind = rk8 ) x
  real ( kind = rk8 ) x1
  real ( kind = rk8 ) x2
  real ( kind = rk8 ) z
!
!  If this is the first call, compute 
!
!    R23 = 2 ^ -23, 
!    R46 = 2 ^ -46,
!    T23 = 2 ^ 23, 
!    T46 = 2 ^ 46.  
!
!  These are computed in loops, rather than by merely using the power operator, 
!  in order to insure that the results are exact on all systems.  
!
  if ( ks == 0 ) then

    r23 = 1.0D+00
    r46 = 1.0D+00
    t23 = 1.0D+00
    t46 = 1.0D+00

    do i = 1, 23
      r23 = 0.5D+00 * r23
      t23 = 2.0D+00 * t23
    end do

    do i = 1, 46
      r46 = 0.50D+00 * r46
      t46 = 2.0D+00 * t46
    end do
!
!  Break A into two parts such that A = 2^23 * A1 + A2.
!
    t1 = r23 * a
    a1 = real ( int ( t1 ), kind = rk8 )
    a2 = a - t23 * a1

    ks = 1

  end if
!
!  Deal with a 0 input value of X.
!
  if ( x == 0.0D+00 ) then
    x = 314159265.0D+00
  end if
!
!  Deal somewhat arbitrarily with negative input X.
!
  if ( x < 0.0D+00 ) then
    x = - x
  end if
!
!  Break X into two parts X1 and X2 such that:
!
!    X = 2^23 * X1 + X2, 
!
!  then compute
!
!    Z = A1 * X2 + A2 * X1  (mod 2^23)
!    X = 2^23 * Z + A2 * X2  (mod 2^46).
!
  t1 = r23 * x
  x1 = real ( int ( t1 ), kind = rk8 )
  x2 = x - t23 * x1

  t1 = a1 * x2 + a2 * x1
  t2 = real ( int ( r23 * t1 ), kind = rk8 )
  z = t1 - t23 * t2

  t3 = t23 * z + a2 * x2
  t4 = real ( int ( r46 * t3 ), kind = rk8 )
  x = t3 - t46 * t4

  randlc = r46 * x

  return
end
function randlc_jump ( x, k )

!*****************************************************************************80
!
!! randlc_jump() returns the K-th element of a uniform pseudorandom sequence.
!
!  Discussion:
!
!    The sequence uses the linear congruential generator:
!
!      X(K+1) = A * X(K)  mod 2^46
!
!    The K-th element, which can be represented as
!
!      X(K) = A^K * X(0)  mod 2^46
!
!    is computed directly using the binary algorithm for exponentiation.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Bailey, Eric Barszcz, John Barton, D Browning, Robert Carter, 
!    Leonardo Dagum, Rod Fatoohi,
!    Samuel Fineberg, Paul Frederickson, Thomas Lasinski, Robert Schreiber, 
!    Horst Simon, V Venkatakrishnan, Sisira Weeratunga,
!    The NAS Parallel Benchmarks,
!    RNR Technical Report RNR-94-007,
!    March 1994.
!
!    Donald Knuth,
!    The Art of Computer Programming,
!    Volume 2, Seminumerical Algorithms,
!    Third Edition,
!    Addison Wesley, 1997,
!    ISBN: 0201896842,
!    LC: QA76.6.K64.
!
!  Parameters:
!
!    Input, real ( kind = rk8 ) X, the initial seed (with index 0).  
!
!    Input, integer K, the index of the desired value.
!
!    Output, real ( kind = rk8 ) RANDLC_JUMP, the K-th value in the sequence.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ), parameter :: a = 1220703125.0D+00
  real ( kind = rk8 ), save :: a1
  real ( kind = rk8 ), save :: a2
  real ( kind = rk8 ) b
  real ( kind = rk8 ) b1
  real ( kind = rk8 ) b2
  integer i
  integer j
  integer k
  integer k2
  integer, save :: ks = 0
  integer m
  real ( kind = rk8 ), save :: r23
  real ( kind = rk8 ), save :: r46
  real ( kind = rk8 ) randlc_jump
  real ( kind = rk8 ) t1
  real ( kind = rk8 ) t2
  real ( kind = rk8 ), save :: t23
  real ( kind = rk8 ) t3
  real ( kind = rk8 ) t4
  real ( kind = rk8 ), save :: t46
  integer twom
  real ( kind = rk8 ) x
  real ( kind = rk8 ) x1
  real ( kind = rk8 ) x2
  real ( kind = rk8 ) xk
  real ( kind = rk8 ) z
!
!  If this is the first call, compute 
!
!    R23 = 2 ^ -23, 
!    R46 = 2 ^ -46,
!    T23 = 2 ^ 23, 
!    T46 = 2 ^ 46.  
!
!  These are computed in loops, rather than by merely using the power operator, 
!  in order to insure that the results are exact on all systems.  
!
  if ( ks == 0 ) then

    r23 = 1.0D+00
    r46 = 1.0D+00
    t23 = 1.0D+00
    t46 = 1.0D+00

    do i = 1, 23
      r23 = 0.5D+00 * r23
      t23 = 2.0D+00 * t23
    end do

    do i = 1, 46
      r46 = 0.50D+00 * r46
      t46 = 2.0D+00 * t46
    end do
!
!  Break A into two parts such that A = 2^23 * A1 + A2.
!
    t1 = r23 * a
    a1 = real ( int ( t1 ), kind = rk8 )
    a2 = a - t23 * a1

    ks = 1

  end if

  if ( k < 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RANDLC_JUMP - Fatal error!'
    write ( *, '(a)' ) '  K < 0.'
    stop 1

  else if ( k == 0 ) then

    xk = x
!
!  Find M so that K < 2^M.
!
  else

    k2 = k
    xk = x

    m = 1
    twom = 2
    do while ( twom <= k )
      twom = twom * 2
      m = m + 1
    end do

    b = a
    b1 = a1
    b2 = a2

    do i = 1, m

      j = k2 / 2
!
!  Replace X by A * X, if appropriate.
!
      if ( 2 * j /= k2 ) then

        t1 = r23 * xk
        x1 = real ( int ( t1 ), kind = rk8 )
        x2 = xk - t23 * x1

        t1 = b1 * x2 + b2 * x1
        t2 = real ( int ( r23 * t1 ), kind = rk8 )
        z = t1 - t23 * t2

        t3 = t23 * z + b2 * x2
        t4 = real ( int ( r46 * t3 ), kind = rk8 )
        xk = t3 - t46 * t4

      end if
!
!  Replace A by A * A mod 2^46.
!
      t1 = r23 * b
      x1 = real ( int ( t1 ), kind = rk8 )
      x2 = b - t23 * x1

      t1 = b1 * x2 + b2 * x1
      t2 = real ( int ( r23 * t1 ), kind = rk8 )
      z = t1 - t23 * t2

      t3 = t23 * z + b2 * x2
      t4 = real ( int ( r46 * t3 ), kind = rk8 )
      b = t3 - t46 * t4
!
!  Update A1, A2.
!
      t1 = r23 * b
      b1 = real ( int ( t1 ), kind = rk8 )
      b2 = b - t23 * b1

      k2 = j

    end do

  end if

  randlc_jump = r46 * xk

  return
end

