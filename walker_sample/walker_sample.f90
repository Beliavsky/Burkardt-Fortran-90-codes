function i4_choose ( n, k )

!*****************************************************************************80
!
!! i4_choose() computes the binomial coefficient C(N,K) as an I4.
!
!  Discussion:
!
!    The value is calculated in such a way as to avoid overflow and
!    roundoff.  The calculation is done in integer arithmetic.
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
!    30 October 2014
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
!    Output, integer I4_CHOOSE, the number of combinations of N
!    things taken K at a time.
!
  implicit none

  integer i
  integer i4_choose
  integer k
  integer mn
  integer mx
  integer n
  integer value

  mn = min ( k, n - k )
  mx = max ( k, n - k )

  if ( mn < 0 ) then

    value = 0

  else if ( mn == 0 ) then

    value = 1

  else

    value = mx + 1

    do i = 2, mn
      value = ( value * ( mx + i ) ) / i
    end do

  end if

  i4_choose = value

  return
end
subroutine i4_choose_test ( )

!*****************************************************************************80
!
!! I4_CHOOSE_TEST tests I4_CHOOSE.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 October 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer cnk
  integer i4_choose
  integer k
  integer n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'I4_CHOOSE_TEST'
  write ( *, '(a)' ) '  I4_CHOOSE evaluates C(N,K).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N         K       CNK'
 
  do n = 0, 4
    write ( *, '(a)' ) ' '
    do k = 0, n
      cnk = i4_choose ( n, k )
      write ( *, '(2x,i8,2x,i8,2x,i8)' ) n, k, cnk
    end do
  end do
 
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

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer a
  integer b
  integer, parameter :: i4_huge = 2147483647
  integer i4_uniform_ab
  integer k
  real ( kind = rk ) r
  integer seed
  integer value

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'i4_uniform_ab(): Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r = real ( seed, kind = rk ) * 4.656612875D-10
!
!  Scale R to lie between A-0.5 and B+0.5.
!
  r = ( 1.0D+00 - r ) * ( real ( min ( a, b ), kind = rk ) - 0.5D+00 ) & 
    +             r   * ( real ( max ( a, b ), kind = rk ) + 0.5D+00 )
!
!  Use rounding to convert R to an integer between A and B.
!
  value = nint ( r )

  value = max ( value, min ( a, b ) )
  value = min ( value, max ( a, b ) )

  i4_uniform_ab = value

  return
end
subroutine i4_uniform_ab_test ( )

!*****************************************************************************80
!
!! I4_UNIFORM_AB_TEST tests I4_UNIFORM_AB.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 October 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: a = -100
  integer, parameter :: b = 200
  integer i
  integer i4_uniform_ab
  integer j
  integer seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'I4_UNIFORM_AB_TEST'
  write ( *, '(a)' ) '  I4_UNIFORM_AB computes pseudorandom values '
  write ( *, '(a)' ) '  in an interval [A,B].'

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The lower endpoint A = ', a
  write ( *, '(a,i12)' ) '  The upper endpoint B = ', b
  write ( *, '(a,i12)' ) '  The initial seed is ', seed
  write ( *, '(a)' ) ' '

  do i = 1, 20

    j = i4_uniform_ab ( a, b, seed )

    write ( *, '(2x,i8,2x,i8)' ) i, j

  end do

  return
end
subroutine normalize ( n, x )

!*****************************************************************************80
!
!! NORMALIZE scales a vector X so its entries sum to 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 February 2016
!
!  Author:
!
!    Original C version by Warren Smith.
!    This FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer N, indicates the size of X.
!
!    Input/output, real ( kind = rk ) X(0:N+1), the vector to be normalized.
!    Entries X(1) through X(N) will sum to 1 on output.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) s
  real ( kind = rk ) x(0:n+1)

  s = sum ( abs ( x(1:n) ) )
  x(1:n) = x(1:n) / s

  return
end
subroutine normalize_test ( )

!*****************************************************************************80
!
!! NORMALIZE_TEST tests NORMALIZE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  real ( kind = rk ), allocatable :: x(:)
  real ( kind = rk ) x_norm

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'normalize_test():'
  write ( *, '(a)' ) '  normalize() normalizes entries 1 through N of a vector'
  write ( *, '(a)' ) '  of length N+2.'

  n = 5
  allocate ( x(0:n+2) )
  call random_number ( harvest = x(1:n+2) )
  call r8vec_print ( n + 2, x, '  Initial X:' )

  x_norm = sum ( abs ( x(1:n) ) )
  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Initial L1 norm of X(1:N) = ', x_norm

  call normalize ( n, x )

  call r8vec_print ( n + 2, x, '  Normalized X:' )

  x_norm = sum ( abs ( x(1:n) ) )
  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Final L1 norm of X(1:N) = ', x_norm
!
!  Free memory.
!
  deallocate ( x )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'NORMALIZE_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'

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
subroutine r8vec_print_test ( )

!*****************************************************************************80
!
!! R8VEC_PRINT_TEST tests R8VEC_PRINT.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 August 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 4

  real ( kind = rk ), dimension ( n ) :: a = (/ &
    123.456D+00, 0.000005D+00, -1.0D+06, 3.14159265D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8vec_print_test():'
  write ( *, '(a)' ) '  R8VEC_PRINT prints an R8VEC.'

  call r8vec_print ( n, a, '  The R8VEC:' )

  return
end
subroutine random_permutation ( n, x, seed )

!*****************************************************************************80
!
!! RANDOM_PERMUTATION applies a random permutation to an array.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 February 2016
!
!  Author:
!
!    Original C version by Warren Smith.
!    This FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer N, indicates the size of X.
!
!    Input/output, real ( kind = rk ) X(0:N+1).  On output, entries 
!    X(1) through X(N) have been randomly permuted.
!
!    Input/output, integer SEED, a seed for the random 
!    number generator.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer i4_uniform_ab
  integer j
  integer seed
  real ( kind = rk ) t
  real ( kind = rk ) x(0:n+1)

  do i = 1, n - 1

    j = i4_uniform_ab ( i, n, seed )

    t    = x(i)
    x(i) = x(j)
    x(j) = t 
     
  end do

  return
end
subroutine random_permutation_test ( )

!*****************************************************************************80
!
!! RANDOM_PERMUTATION_TEST tests RANDOM_PERMUTATION.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer n
  integer seed
  real ( kind = rk ), allocatable :: x(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'RANDOM_PERMUTATION_TEST'
  write ( *, '(a)' ) '  RANDOM_PERMUTATION randomly permutes entries 1 through'
  write ( *, '(a)' ) '  N of a vector X[0:N+1].'

  n = 5
  allocate ( x(0:n+1) )
  do i = 0, n + 1
    x(i) = real ( i, kind = rk )
  end do
  seed = 123456789
  call r8vec_print ( n + 2, x, '  Initial X:' )
  call random_permutation ( n, x, seed )
  call r8vec_print ( n + 2, x, '  Permuted X:' )
!
!  Free memory.
!
  deallocate ( x )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'RANDOM_PERMUTATION_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'

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
!  Parameters:
!
!    None
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
subroutine walker_build ( n, x, y, a )

!*****************************************************************************80
!
!! WALKER_BUILD sets up the data for a Walker sampler.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 February 2016
!
!  Author:
!
!    Original C version by Warren Smith.
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Warren Smith,
!    How to sample from a probability distribution,
!    April 18, 2002.
!
!  Parameters:
!
!    Input, integer N, indicates the size of X.
!
!    Input, real ( kind = rk ) X(0:N+1), contains in X(1) through X(N) the
!    probabilities of outcomes 1 through N.  
!
!    Output, real ( kind = rk ) Y(0:N+1), the Walker threshold vector.
!
!    Output, integer A(0:N+1), the Walker index vector.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer a(0:n+1)
  integer b(0:n+1)
  integer i
  integer j
  integer k
  real ( kind = rk ) x(0:n+1)
  real ( kind = rk ) y(0:n+1)
!
!  Initialize A.
!
  do i = 0, n + 1
    a(i) = i
  end do
!
!  Initialize B to the "stay here" value, and set sentinel values at the ends.
!
  do i = 0, n + 1
    b(i) = i
  end do
!
!  Copy Y from X.
!  Scale the probability vector and set sentinel values at the ends.
!
  y(0) = 0.0D+00
  y(1:n) = x(1:n) * real ( n, kind = rk )
  y(n+1) = 2.0D+00

  i = 0
  j = n + 1

  do
!
!  Find i so Y(B(i)) needs more.
!
    do  
      i = i + 1
      if ( 1.0D+00 <= y(b(i)) ) then
        exit
      end if
    end do
!
!  Find j so Y(B(j)) wants less.
!
    do
      j = j - 1
      if ( y(b(j)) < 1.0D+00 ) then
        exit
      end if
    end do

    if ( j <= i ) then
      exit
    end if
!
!  Swap B(i) and B(j).
!
    k    = b(i)
    b(i) = b(j)
    b(j) = k

  end do

  i = j
  j = j + 1

  do while ( 0 < i )
!
!  Find J such that Y(B(j)) needs more.
!
    do while ( y(b(j)) <= 1.0 )
      j = j + 1
    end do
!
!  Meanwhile, Y(B(i)) wants less.
!
    if ( n < j ) then
      exit
    end if
!
!  B(i) will donate to B(j) to fix up.
!
    y(b(j)) = y(b(j)) - ( 1.0D+00 - y(b(i)) )     
    a(b(i)) = b(j)             
! 
!  Y(B(j)) now wants less so readjust ordering.
!
    if ( y(b(j)) < 1.0D+00 ) then

      k    = b(i)
      b(i) = b(j)
      b(j) = k
      j = j + 1

    else

      i = i - 1

    end if

  end do

  return
end
subroutine walker_build_test ( )

!*****************************************************************************80
!
!! WALKER_BUILD_TEST tests WALKER_BUILD.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, allocatable :: a(:)
  integer i
  integer i4_choose
  integer n
  real ( kind = rk ), allocatable :: x(:)
  real ( kind = rk ), allocatable :: y(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'WALKER_BUILD_TEST'
  write ( *, '(a)' ) '  WALKER_BUILD builds the Walker sampler data vectors Y'
  write ( *, '(a)' ) '  and A, given a probability vector X.'

  n = 5
  allocate ( x(0:n+1) )

  do i = 1, n
    x(i) = real ( i4_choose ( n - 1, i - 1 ), kind = rk ) &
         / real ( 2 ** ( n - 1 ), kind = rk )
  end do

  call r8vec_print ( n + 2, x, &
    '  Binomial PDF (ignore first and last entries):' )

  allocate ( y(0:n+1) )
  allocate ( a(0:n+1) )

  call walker_build ( n, x, y, a )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '   I    A[I]    Y[i] (ignore first and last entries)'
  write ( *, '(a)' ) ''
  do i = 0, n + 1
    write ( *, '(2x,i2,2x,i2,2x,g14.6)' ) i, a(i), y(i)
  end do
!
!  Free memory.
!
  deallocate ( a )
  deallocate ( x )
  deallocate ( y )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'WALKER_BUILD_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'

  return
end
subroutine walker_sampler ( n, y, a, seed, i )

!*****************************************************************************80
!
!! WALKER_SAMPLER returns a random sample i=1..N with prob X(i).
!
!  Discussion:
!
!    Implementation of algorithm for sampling from a discrete
!    probability N-vector X(1), X(2), ..., X(N).  (N>=1.)
!    Runs on O(1) worst case time per sample,
!    and uses one integer and one real N-element array for storage.
!    Preprocessing consumes O(N) time and temporarily uses one 
!    additional integer array (B(0..N+1)) for bookkeeping. 
!    X(0) and X(N+1) are also used as sentinels in the Build() algorithm.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 February 2016
!
!  Author:
!
!    Original C version by Warren Smith.
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Warren Smith,
!    How to sample from a probability distribution,
!    April 18, 2002.
!
!  Parameters:
!
!    Input, integer N, indicates the size of X.
!
!    Input, real ( kind = rk ) Y(0:N+1), the Walker threshold vector.
!
!    Input, integer A(0:N+1), the Walker index vector.
!
!    Input/output, integer SEED, a seed for the random 
!    number generator.
!
!    Output, integer I, a sample value between 1 and N,
!    selected according to the probability vector X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer a(0:n+1)
  integer n
  integer i 
  integer i4_uniform_ab
  real ( kind = rk ) r
  integer seed
  real ( kind = rk ) y(0:n+1)
! 
!  Let i = random uniform integer from {1,2,...N}  
!
  i = i4_uniform_ab ( 1, n, seed ) 

  call random_number ( harvest = r )

  if ( y(i) < r ) then
    i = a(i)
  end if

  return
end
subroutine walker_sampler_test ( )

!*****************************************************************************80
!
!! WALKER_SAMPLER_TEST tests WALKER_SAMPLER.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 February 2016
!
!  Author:
!
!    Original C version by Warren Smith.
!    This FORTRAN90 version by John Burkardt.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, allocatable :: a(:)
  integer, allocatable :: count(:)
  real ( kind = rk ) expval
  integer i
  integer j
  integer n
  real ( kind = rk ) p
  real ( kind = rk ) s
  integer seed
  real ( kind = rk ) t
  real ( kind = rk ), allocatable :: x(:)
  real ( kind = rk ), allocatable :: y(:)

  seed = 123456789
  n = 10
  p = 2.0D+00

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'WALKER_SAMPLER_TEST:'
  write ( *, '(a)' ) '  WALKER_SAMPLER creates Walker sample vectors Y and A'
  write ( *, '(a)' ) '  for efficiently sampling a discrete probability vector.'
  write ( *, '(a)' ) '  Test the Walker sampler with a Zipf-type probability.'
!
!  Generate a standard Zipf probability vector for cases 1,...,N,
!  with parameter P.
!
  allocate ( x(0:n+1) )
  call zipf_probability ( n, p, x )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Zipf probabilities'
  write ( *, '(a,i4)' ) '  for N = ', n
  write ( *, '(a,g14.6)' ) '  and parameter P = ', p
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '     I     X(I)'
  write ( *, '(a)' ) ''
  do i = 1, n
    write ( *, '(2x,i4,2x,g14.6)' ) i, x(i)
  end do
!
!  For better testing, randomly scramble the probabilities.
!
  seed = 123456789
  call random_permutation ( n, x, seed )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Randomly permuted X:'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '     I     X(I)'
  write ( *, '(a)' ) ''
  do i = 1, n
    write ( *, '(2x,i4,2x,g14.6)' ) i, x(i)
  end do
!
!  Build the Walker sampler.
!
  allocate ( y(0:n+1) )
  allocate ( a(0:n+1) )

  call walker_build ( n, x, y, a )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Built the sampler'
  write ( *, '(a)' ) '  i Y(i) A(i):'
  write ( *, '(a)' ) ''

  do i = 1, n
    write ( *, '(2x,i4,2x,g14.6,2x,i4)' ) i, y(i), a(i)
  end do
!
!  Prepare to count the frequency of each outcome.
!
  allocate ( count(0:n+1) )
  count(0:n+1) = 0
!
!  Call the sampler many times.
!
  do i = 1, 100000
    call walker_sampler ( n, y, a, seed, j )
    count(j) = count(j) + 1
  end do
!
!  Compare normalized sample frequencies to the original probabilities in X.
!
  s = 0.0D+00
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  100000 samples:'
  write ( *, '(a)' ) '  prob   #samples:'
  write ( *, '(a)' ) ''

  do i = 1, n
    write ( *, '(2x,g14.6,2x,i6)' ) x(i), count(i)
    expval = x(i) * 100000
    t = expval - count(i)
    s = s + t * t / expval
  end do
  s = s / real ( n, kind = rk )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6,a)' ) '  sumvar = ', s, ' (should be about 1)'

  return
end
subroutine walker_verify ( n, x, y, a, v )

!*****************************************************************************80
!
!! WALKER_VERIFY verifies a Walker Sampler structure.
!
!  Discussion:
!
!    This test applies the sampling algorithms to a Zipfian distribution.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 February 2016
!
!  Author:
!
!    Original C version by Warren Smith.
!    This FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer N, indicates the size of X.
!
!    Input, real ( kind = rk ) X(0:N+1), contains in X(1) through X(N) the
!    probabilities of outcomes 1 through N.
!
!    Input, real ( kind = rk ) Y(0:N+1), the Walker threshold vector.
!
!    Input, integer A(0:N+1), the Walker index vector.
!
!    Output, real ( kind = rk ) V, the verification sum, which
!    should be close to zero.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer a(0:n+1)
  integer i
  real ( kind = rk ) v
  real ( kind = rk ) x(0:n+1)
  real ( kind = rk ) y(0:n+1)
  real ( kind = rk ) z(0:n+1)
!
!  Reverse the scaling.
!
  z(0:n+1) = y(0:n+1) / real ( n, kind = rk )
!
!  Add back the adjustments.
!  (Don't try to vectorize this statement!)
!
  do i = 1, n 
    z(a(i)) = z(a(i)) + ( 1.0D+00 - y(i) ) / real ( n, kind = rk )
  end do
!
!  Check for discrepancies between Z and X.
!
  v = sum ( abs ( z(1:n) - x(1:n) ) )

  return
end
subroutine walker_verify_test ( )

!*****************************************************************************80
!
!! WALKER_VERIFY_TEST tests WALKER_VERIFY.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, allocatable :: a(:)
  integer i
  integer n
  real ( kind = rk ) v
  real ( kind = rk ), allocatable :: x(:)
  real ( kind = rk ), allocatable :: y(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'WALKER_VERIFY_TEST'
  write ( *, '(a)' ) '  WALKER_VERIFY verifies the Walker sampler data'
  write ( *, '(a)' ) '  vectors Y and A,for a given probability vector X.'

  n = 9
  allocate ( x(0:n+1) )

  x(0) = 0.0D+00
  do i = 1, n
    x(i) = log ( 1.0D+00 + 1.0D+00 / real ( i, kind = rk ) ) &
      / log ( real ( n + 1, kind = rk ) )
  end do
  x(n+1) = 0.0D+00

  call r8vec_print ( n + 2, x, &
    '  Benford PDF (ignore first and last entries):' )

  allocate ( y(0:n+1) )
  allocate ( a(0:n+1) )

  call walker_build ( n, x, y, a )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '   I    A(I)    Y(i) (ignore first and last entries)'
  write ( *, '(a)' ) ''
  do i = 0, n + 1
    write ( *, '(2x,i2,2x,i2,g14.6)' ) i, a(i), y(i)
  end do

  call walker_verify ( n, x, y, a, v )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  The verification sum = ', v
  write ( *, '(a)' ) '  It should be very close to zero.'
!
!  Free memory.
!
  deallocate ( a )
  deallocate ( x )
  deallocate ( y )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'WALKER_VERIFY_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'

  return
end
subroutine zipf_probability ( n, p, x )

!*****************************************************************************80
!
!! ZIPF_PROBABILITY sets up a Zipf probability vector.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 February 2016
!
!  Author:
!
!    Original C version by Warren Smith.
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    George Zipf,
!    The Psychobiology of Language,
!    1935.
!
!  Parameters:
!
!    Input, integer N, indicates the size of X.
!
!    Input, real ( kind = rk ) P, the Zipf parameter.
!    1.0 < P.
!
!    Output, real ( kind = rk ) X(0:N+1), contains in X(1) through X(N) the
!    probabilities of outcomes 1 through N.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  real ( kind = rk ) p
  real ( kind = rk ) x(0:n+1)

  x(0) = 0.0D+00
  do i = 1, n
    x(i) = real ( i, kind = rk ) ** ( - p )
  end do
  x(n+1) = 0.0D+00

  call normalize ( n, x )

  return
end
subroutine zipf_probability_test ( )

!*****************************************************************************80
!
!! ZIPF_PROBABILITY_TEST tests ZIPF_PROBABILITY.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  real ( kind = rk ) p
  real ( kind = rk ), allocatable :: x(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'ZIPF_PROBABILITY_TEST'
  write ( *, '(a)' ) '  ZIPF_PROBABILITY sets up a probablity vector X of N+2'
  write ( *, '(a)' ) '  elements containing in X[1:N] the probabilities of'
  write ( *, '(a)' ) '  outcomes 1 through Nin a Zipf distribution with'
  write ( *, '(a)' ) '  parameter P.'

  n = 5
  p = 1.0D+00
  allocate ( x(0:n+1) )
  call zipf_probability ( n, p, x )
  call r8vec_print ( n + 2, x, '  X for N = 5, P = 1.0' )
  deallocate ( x )

  n = 5
  p = 2.0D+00
  allocate ( x(0:n+1) )
  call zipf_probability ( n, p, x )
  call r8vec_print ( n + 2, x, '  X for N = 5, P = 2.0' )
  deallocate ( x )

  n = 10
  p = 2.0D+00
  allocate ( x(0:n+1) )
  call zipf_probability ( n, p, x )
  call r8vec_print ( n + 2, x, '  X for N = 10, P = 2.0' )
  deallocate ( x )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'ZIPF_PROBABILITY_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'

  return
end
