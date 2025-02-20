function i4_factorial ( n )

!*****************************************************************************80
!
!! i4_factorial() computes the factorial of N.
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
!    26 June 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the argument of the factorial function.
!    If N is less than 1, the function value is returned as 1.
!    0 <= N <= 13 is required.
!
!    Output, integer I4_FACTORIAL, the factorial of N.
!
  implicit none

  integer i
  integer i4_factorial
  integer n

  i4_factorial = 1

  if ( 13 < n ) then
    i4_factorial = - 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_FACTORIAL - Fatal error!'
    write ( *, '(a)' ) '  I4_FACTORIAL(N) cannot be computed as an integer'
    write ( *, '(a)' ) '  for 13 < N.'
    write ( *, '(a,i8)' ) '  Input value N = ', n
    stop
  end if

  do i = 1, n
    i4_factorial = i4_factorial * i
  end do

  return
end
function i4_modp ( i, j )

!*****************************************************************************80
!
!! i4_modp() returns the nonnegative remainder of I4 division.
!
!  Discussion:
!
!    If
!      NREM = I4_MODP ( I, J )
!      NMULT = ( I - NREM ) / J
!    then
!      I = J * NMULT + NREM
!    where NREM is always nonnegative.
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
!
!    An I4 is an integer value.
!
!  Example:
!
!        I     J     MOD I4_MODP    Factorization
!
!      107    50       7       7    107 =  2 *  50 + 7
!      107   -50       7       7    107 = -2 * -50 + 7
!     -107    50      -7      43   -107 = -3 *  50 + 43
!     -107   -50      -7      43   -107 =  3 * -50 + 43
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, the number to be divided.
!
!    Input, integer J, the number that divides I.
!
!    Output, integer I4_MODP, the nonnegative remainder when I is
!    divided by J.
!
  implicit none

  integer i
  integer i4_modp
  integer j
  integer value

  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MODP - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal divisor J = ', j
    stop
  end if

  value = mod ( i, j )

  if ( value < 0 ) then
    value = value + abs ( j )
  end if

  i4_modp = value

  return
end
function i4_uniform_ab ( a, b )

!*****************************************************************************80
!
!! i4_uniform_ab() returns a scaled pseudorandom I4 between A and B.
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
!  Input:
!
!    integer A, B, the limits of the interval.
!
!  Output:
!
!    integer I4_UNIFORM_AB, a number between A and B.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer a
  integer b
  integer i4_uniform_ab
  real ( kind = rk ) r

  call random_number ( harvest = r )
  i4_uniform_ab = a + int ( ( b + 1 - a ) * r )

  return
end
function i4_wrap ( ival, ilo, ihi )

!*****************************************************************************80
!
!! i4_wrap() forces an I4 to lie between given limits by wrapping.
!
!  Discussion:
!
!    An I4 is an integer value.
!
!    There appears to be a bug in the GFORTRAN compiler which can lead to
!    erroneous results when the first argument of I4_WRAP is an expression.
!    In particular:
!
!    do i = 1, 3
!      if ( test ) then
!        i4 = i4_wrap ( i + 1, 1, 3 )
!      end if
!    end do
!
!    was, when I = 3, returning I4 = 3.  So I had to replace this with
!
!    do i = 1, 3
!      if ( test ) then
!        i4 = i + 1
!        i4 = i4_wrap ( i4, 1, 3 )
!      end if
!    end do
!
!  Example:
!
!    ILO = 4, IHI = 8
!
!    I  Value
!
!    -2     8
!    -1     4
!     0     5
!     1     6
!     2     7
!     3     8
!     4     4
!     5     5
!     6     6
!     7     7
!     8     8
!     9     4
!    10     5
!    11     6
!    12     7
!    13     8
!    14     4
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer IVAL, a value.
!
!    Input, integer ILO, IHI, the desired bounds.
!
!    Output, integer I4_WRAP, a "wrapped" version of the value.
!
  implicit none

  integer i4_modp
  integer i4_wrap
  integer ihi
  integer ilo
  integer ival
  integer jhi
  integer jlo
  integer value
  integer wide

  jlo = min ( ilo, ihi )
  jhi = max ( ilo, ihi )

  wide = jhi - jlo + 1

  if ( wide == 1 ) then
    value = jlo
  else
    value = jlo + i4_modp ( ival - jlo, wide )
  end if

  i4_wrap = value

  return
end
subroutine i4vec_indicator ( n, a )

!*****************************************************************************80
!
!! i4vec_indicator() sets an I4VEC to the indicator vector.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements of A.
!
!    Output, integer A(N), the array to be initialized.
!
  implicit none

  integer n

  integer a(n)
  integer i

  do i = 1, n
    a(i) = i
  end do

  return
end
subroutine i4vec_reverse ( n, a )

!*****************************************************************************80
!
!! i4vec_reverse() reverses the elements of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    In FORTRAN90, call I4VEC_REVERSE is equivalent to:
!
!      A(1:N) = A(N:1:-1)
!
!  Example:
!
!    Input:
!
!      N = 5,
!      A = ( 11, 12, 13, 14, 15 ).
!
!    Output:
!
!      A = ( 15, 14, 13, 12, 11 ).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input/output, integer A(N), the array to be reversed.
!
  implicit none

  integer n

  integer a(n)

  a(1:n) = a(n:1:-1)

  return
end
subroutine perm_check ( n, p )

!*****************************************************************************80
!
!! perm_check() checks a representation of a permutation.
!
!  Discussion:
!
!    The routine is given N and P, a vector of length N.
!    P is a legal represention of a permutation of the integers from
!    1 to N if and only if every integer from 1 to N occurs
!    as a value of P(I) for some I between 1 and N.
!
!    If an error is observed in the input, this routine prints a message
!    and stops.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of values being permuted.
!    N must be positive.
!
!    Input, integer P(N), the array to check.
!
  implicit none

  integer n

  integer i
  integer ifind
  integer iseek
  integer p(n)

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_CHECK - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    write ( *, '(a,i12)' ) '  N = ', n
    stop
  end if

  do i = 1, n
    if ( p(i) < 1 .or. n < p(i) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PERM_CHECK - Fatal error!'
      write ( *, '(a)' ) '  P(I) < 1 or N < P(I).'
      write ( *, '(a,i12)' ) '  N = ', n
      write ( *, '(a,i12)' ) '  I = ', i
      write ( *, '(a,i12)' ) '  P(I) = ', p(i)
      stop
    end if
  end do

  do iseek = 1, n

    ifind = -1

    do i = 1, n
      if ( p(i) == iseek ) then
        ifind = i
        exit
      end if
    end do

    if ( ifind == -1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PERM_CHECK - Fatal error!'
      write ( *, '(a)' ) '  Every I from 1 to N must occur in P.'
      write ( *, '(a)' ) '  Could not find occurrence of ', iseek
      stop
    end if

  end do

  return
end
subroutine perm_enum ( n, nperm )

!*****************************************************************************80
!
!! perm_enum() enumerates the permutations on N digits.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of values being permuted.
!    N must be nonnegative.
!
!    Output, integer NPERM, the number of distinct elements.
!
  implicit none

  integer i4_factorial
  integer n
  integer nperm

  nperm = i4_factorial ( n )

  return
end
subroutine perm_inverse ( n, p, pinv  )

!*****************************************************************************80
!
!! perm_inverse() computes the inverse of a permutation.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer N, the number of values being permuted.
!    N must be positive.
!
!    Input, integer P(N), describes the permutation.
!    P(I) is the item which is permuted into the I-th place
!    by the permutation.
!
!    Output, integer PINV(N), the inverse permutation.
!
  implicit none

  integer n

  integer i
  integer p(n)
  integer pinv(n)
!
!  Check.
!
  call perm_check ( n, p )

  do i = 1, n
    pinv(p(i)) = i
  end do

  return
end
function perm_is_unicycle ( n, p )

!*****************************************************************************80
!
!! perm_is_unicycle() is TRUE if a permutation is a unicycle.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 June 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of objects in the permutation.
!
!    Input, integer P(N), the permutation.
!
!    Output, logical PERM_IS_UNICYCLE, is TRUE if the permutation is a unicycle.
!
  implicit none

  integer n

  integer i
  integer j
  integer p(n)
  logical perm_is_unicycle

  perm_is_unicycle = .false.

  call perm_check ( n, p )
!
!  From 1, you must be able to take N-1 steps without reaching 1...
!
  i = 1
  do j = 1, n - 1
    i = p(i)
    if ( i == 1 ) then
      return
    end if
  end do
!
!  ...and the N-th step must reach 1.
!
  i = p(i)
  if ( i == 1 ) then
    perm_is_unicycle = .true.
  end if

  return
end
subroutine perm_lex_next ( n, p, rank )

!*****************************************************************************80
!
!! perm_lex_next() computes the lexicographic permutation successor.
!
!  Example:
!
!    RANK  Permutation
!
!       0  1 2 3 4
!       1  1 2 4 3
!       2  1 3 2 4
!       3  1 3 4 2
!       4  1 4 2 3
!       5  1 4 3 2
!       6  2 1 3 4
!       ...
!      23  4 3 2 1
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer N, the number of values being permuted.
!    N must be positive.
!
!    Input/output, integer P(N), describes the permutation.
!    P(I) is the item which is permuted into the I-th place
!    by the permutation.
!
!    Input/output, integer RANK, the rank.
!    If RANK = -1 on input, then the routine understands that this is
!    the first call, and that the user wishes the routine to supply
!    the first element in the ordering, which has RANK = 0.
!    In general, the input value of RANK is increased by 1 for output,
!    unless the very last element of the ordering was input, in which
!    case the output value of RANK is -1.
!
  implicit none

  integer n

  integer i
  integer j
  integer p(n)
  integer rank
  integer temp
!
!  Return the first element.
!
  if ( rank == -1 ) then
    call i4vec_indicator ( n, p )
    rank = 0
    return
  end if
!
!  Check.
!
  call perm_check ( n, p )
!
!  Seek I, the highest index for which the next element is bigger.
!
  i = n - 1

  do

    if ( i <= 0 ) then
      exit
    end if

    if ( p(i) <= p(i+1) ) then
      exit
    end if

    i = i - 1

  end do
!
!  If no I could be found, then we have reach the final permutation,
!  N, N-1, ..., 2, 1.  Time to start over again.
!
  if ( i == 0 ) then
    call i4vec_indicator ( n, p )
    rank = -1
  else
!
!  Otherwise, look for the the highest index after I whose element
!  is bigger than I's.  We know that I+1 is one such value, so the
!  loop will never fail.
!
    j = n
    do while ( p(j) < p(i) )
      j = j - 1
    end do
!
!  Interchange elements I and J.
!
    temp = p(i)
    p(i) = p(j)
    p(j) = temp
!
!  Reverse the elements from I+1 to N.
!
    call i4vec_reverse ( n - i, p(i+1:n) )
    rank = rank + 1

  end if

  return
end
subroutine perm_lex_rank ( n, p, rank )

!*****************************************************************************80
!
!! perm_lex_rank() computes the lexicographic rank of a permutation.
!
!  Discussion:
!
!    The original code altered the input permutation.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer N, the number of values being permuted.
!    N must be positive.
!
!    Input, integer P(N), describes the permutation.
!    P(I) is the item which is permuted into the I-th place
!    by the permutation.
!
!    Output, integer RANK, the rank of the permutation.
!
  implicit none

  integer n

  integer i
  integer i4_factorial
  integer j
  integer p(n)
  integer pcopy(n)
  integer rank
!
!  Check.
!
  call perm_check ( n, p )

  rank = 0
  pcopy(1:n) = p(1:n)

  do j = 1, n

    rank = rank + ( pcopy(j) - 1 ) * i4_factorial ( n - j )

    do i = j + 1, n
      if ( pcopy(j) < pcopy(i) ) then
        pcopy(i) = pcopy(i) - 1
      end if
    end do

  end do

  return
end
subroutine perm_lex_unrank ( n, rank, p )

!*****************************************************************************80
!
!! perm_lex_unrank() computes the permutation of given lexicographic rank.
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
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer N, the number of values being permuted.
!    N must be positive.
!
!    Input, integer RANK, the rank of the permutation.
!
!    Output, integer P(N), describes the permutation.
!
  implicit none

  integer n

  integer d
  integer i
  integer i4_factorial
  integer j
  integer nperm
  integer p(n)
  integer rank
  integer rank_copy
!
!  Check.
!
  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_LEX_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  Input N is illegal.'
    stop
  end if

  call perm_enum ( n, nperm )

  if ( rank < 0 .or. nperm < rank ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_LEX_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  The input rank is illegal.'
    stop
  end if

  rank_copy = rank

  p(n) = 1

  do j = 1, n - 1

    d = mod ( rank_copy, i4_factorial ( j + 1 ) ) / i4_factorial ( j )
    rank_copy = rank_copy - d * i4_factorial ( j )
    p(n-j) = d + 1

    do i = n - j + 1, n

      if ( d < p(i) ) then
        p(i) = p(i) + 1
      end if

    end do

  end do

  return
end
subroutine perm_print ( n, p, title )

!*****************************************************************************80
!
!! perm_print() prints a permutation.
!
!  Example:
!
!    Input:
!
!      P = 7 2 4 1 5 3 6
!
!    Printed output:
!
!      "This is the permutation:"
!
!      1 2 3 4 5 6 7
!      7 2 4 1 5 3 6
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of objects permuted.
!
!    Input, integer P(N), the permutation, in standard index form.
!
!    Input, character ( len = * ) TITLE, a title.
!    If no title is supplied, then only the permutation is printed.
!
  implicit none

  integer n

  integer i
  integer ihi
  integer ilo
  integer, parameter :: inc = 20
  integer p(n)
  character ( len = * ) title

  if ( len_trim ( title ) /= 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )

    do ilo = 1, n, inc
      ihi = min ( n, ilo + inc - 1 )
      write ( *, '(a)' ) ' '
      write ( *, '(2x,20i4)' ) ( i, i = ilo, ihi )
      write ( *, '(2x,20i4)' ) p(ilo:ihi)
    end do

  else

    do ilo = 1, n, inc
      ihi = min ( n, ilo + inc - 1 )
      write ( *, '(2x,20i4)' ) p(ilo:ihi)
    end do

  end if

  return
end
subroutine perm_random ( n, p )

!*****************************************************************************80
!
!! perm_random() selects a random permutation of N objects.
!
!  Discussion:
!
!    The routine assumes the objects are labeled 1, 2, ... N.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 May 2002
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    This version by John Burkardt.
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
!  Input:
!
!    integer N, the number of objects to be permuted.
!
!  Output:
!
!    integer P(N), a permutation of ( 1, 2, ..., N ),
!    in standard index form.
!
  implicit none

  integer n

  integer i
  integer i4_uniform_ab
  integer j
  integer p(n)
  integer t

  call i4vec_indicator ( n, p )

  do i = 1, n - 1

    j = i4_uniform_ab ( i, n )

    t    = p(i)
    p(i) = p(j)
    p(j) = t

  end do

  return
end
subroutine unicycle_check ( n, p )

!*****************************************************************************80
!
!! unicycle_check() checks that a vector represents a unicycle.
!
!  Discussion:
!
!    A unicycle is a permutation with a single cycle.  This might be called
!    a cyclic permutation, except that that name is used with at least two
!    other meanings.  So, to be clear, a unicycle is a permutation of N
!    objects in which each object is returned to itself precisely after 
!    N applications of the permutation.
!
!    This routine verifies that each of the integers from 1
!    to N occurs among the N entries of the permutation.
!
!    Any permutation of the integers 1 to N describes a unicycle.
!    The permutation ( 3, 4, 2, 1 ) indicates that the unicycle
!    sends 3 to 4, 4 to 2, 2 to 1 and 1 to 3.  This is the sequential
!    description of a unicycle.
!
!    The standard sequence "rotates" the permutation so that it begins
!    with 1.  The above sequence becomes a standard sequence when
!    written as ( 1, 3, 4, 2 ).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 June 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries.
!
!    Input, integer P(N), the unicycle sequence vector
!
  implicit none

  integer n

  integer ierror
  integer ifind
  integer iseek
  integer p(n)

  ierror = 0

  do iseek = 1, n

    ierror = iseek

    do ifind = 1, n
      if ( p(ifind) == iseek ) then
        ierror = 0
        exit
      end if
    end do

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'UNICYCLE_CHECK - Fatal error!'
      write ( *, '(a)' ) '  The input array does not represent'
      write ( *, '(a)' ) '  a unicycle.  In particular, the'
      write ( *, '(a,i8)' ) '  array is missing the value ', ierror
      stop
    end if

  end do

  return
end
subroutine unicycle_enum ( n, num )

!*****************************************************************************80
!
!! unicycle_enum() enumerates the unicycles.
!
!  Discussion:
!
!    Each standard sequence corresponds to a unique unicycle.  Since the
!    first element of a standard sequence is always 1, the number of standard
!    sequences, and hence the number of unicycles, is (n-1)!.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 June 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the unicyle.
!
!    Output, integer NUM, the number of unicycles.
!
  implicit none

  integer i4_factorial
  integer n
  integer num

  num = i4_factorial ( n - 1 )

  return
end
subroutine unicycle_index ( n, u, u_index )

!*****************************************************************************80
!
!! unicycle_index() returns the index form of a unicycle.
!
!  Example:
!
!    N = 4
!
!    U       = 1 3 4 2
!    U_INDEX = 3 1 4 2
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 June 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the unicycles.
!
!    Input, integer U(N), the unicycle sequence vector.
!
!    Output, integer U_INDEX(N), the unicycle index vector.
!
  implicit none

  integer n

  integer i
  integer i4_wrap
  integer ip1
  integer u(n)
  integer u_index(n)

  do i = 1, n
    ip1 = i4_wrap ( i + 1, 1, n )
    u_index(u(i)) = u(ip1)
  end do
  
  return
end
subroutine unicycle_index_print ( n, u_index, title )

!*****************************************************************************80
!
!! unicycle_index_print() prints a unicycle given in index form.
!
!  Example:
!
!    Input:
!
!      U_INDEX = 7 1 4 5 2 3 6
!
!    Printed output:
!
!      1 2 3 4 5 6 7
!      7 1 4 5 2 3 6
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 June 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the unicycle.
!
!    Input, integer U_INDEX(N), the unicycle index vector.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer n

  integer i
  integer ihi
  integer ilo
  integer, parameter :: inc = 20
  integer u_index(n)
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do ilo = 1, n, inc
    ihi = min ( n, ilo + inc - 1 )
    write ( *, '(a)' ) ' '
    write ( *, '(2x,20i4)' ) ( i, i = ilo, ihi )
    write ( *, '(2x,20i4)' ) u_index(ilo:ihi)
  end do

  return
end
subroutine unicycle_index_to_sequence ( n, u_index, u )

!*****************************************************************************80
!
!! unicycle_index_to_sequence() converts a unicycle from index to sequence form.
!
!  Example:
!
!    N = 4
!
!    U_INDEX = 3 1 4 2
!    U       = 1 3 4 2
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 June 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the unicycles.
!
!    Output, integer U_INDEX(N), the unicycle index vector.
!
!    Input, integer U(N), the unicycle sequence vector.
!
  implicit none

  integer n

  integer i
  integer j
  integer u(n)
  integer u_index(n)

  u(1) = 1
  i = 1

  do j = 2, n

    i = u_index(i)
    u(j) = i

    if ( i == 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'UNICYCLE_INDEX_TO_SEQUENCE - Fatal error!'
      write ( *, '(a)' ) '  The index vector does not represent a unicycle.'
      write ( *, '(a,i4,a,i4,a)' ) '  On step ', j, ' u_index(', i, ') = 1.'
      stop
    end if

  end do
  
  return
end
subroutine unicycle_inverse ( n, u, u_inverse )

!*****************************************************************************80
!
!! unicycle_inverse() returns the inverse of a unicycle.
!
!  Example:
!
!    N = 4
!
!    U         = 1 3 4 2
!    U_INVERSE = 1 2 4 3
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 June 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the unicycles.
!
!    Input, integer U(N), the unicycle sequence vector.
!
!    Output, integer U_INVERSE(N), the inverse unicycle.
!
  implicit none

  integer n

  integer u(n)
  integer u_inverse(n)

  u_inverse(1) = 1
  u_inverse(2:n) = u(n:2:-1)

  return
end
subroutine unicycle_next ( n, u, rank )

!*****************************************************************************80
!
!! unicycle_next() generates unicycles in lexical order, one at a time.
!
!  Example:
!
!    N = 4
!
!    1   1 2 3 4
!    2   1 2 4 3
!    3   1 3 2 4
!    4   1 3 4 2
!    5   1 4 2 3
!    6   1 4 3 2
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 June 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the unicycles.
!
!    Input/output, integer U(N); on first call with MORE = FALSE,
!    this value is not used.  Otherwise, the input value is the previous
!    unicycle.  The output value is the next unicycle.
!
!    Input/output, integer RANK, the rank.
!    If RANK = -1 on input, then the routine understands that this is
!    the first call, and that the user wishes the routine to supply
!    the first element in the ordering, which has RANK = 0.
!    In general, the input value of RANK is increased by 1 for output,
!    unless the very last element of the ordering was input, in which
!    case the output value of RANK is -1.
!
  implicit none

  integer n

  integer p(n-1)
  integer rank
  integer u(n)

  if ( rank == -1 ) then
    u(1) = 1
  else
    p(1:n-1) = u(2:n) - 1
  end if

  call perm_lex_next ( n - 1, p, rank )
 
  u(2:n) = p(1:n-1) + 1

  return
end
subroutine unicycle_print ( n, u, title )

!*****************************************************************************80
!
!! unicycle_print() prints a unicycle given in sequence form.
!
!  Example:
!
!    Input:
!
!      U = 7 1 4 5 2 3 6
!
!    Printed output:
!
!      7 1 4 5 2 3 6
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 June 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the unicycle.
!
!    Input, integer U(N), the unicycle sequence vector.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer n

  integer ihi
  integer ilo
  integer, parameter :: inc = 20
  integer u(n)
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
    write ( *, '(a)' ) ' '
  end if

  do ilo = 1, n, inc
    ihi = min ( n, ilo + inc - 1 )
    write ( *, '(2x,20i4)' ) u(ilo:ihi)
  end do

  return
end
subroutine unicycle_random ( n, u )

!*****************************************************************************80
!
!! unicycle_random() selects a random unicycle of N objects.
!
!  Discussion:
!
!    The routine assumes the objects are labeled 1, 2, ... N.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 June 2012
!
!  Author:
!
!    John Burkardt.
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
!  Input:
!
!    integer N, the number of objects to be permuted.
!
!  Output:
!
!    integer U(N), a unicycle in sequence form.
!
  implicit none

  integer n

  integer i
  integer i4_uniform_ab
  integer j
  integer u(n)
  integer t

  call i4vec_indicator ( n, u )

  do i = 2, n
    j = i4_uniform_ab ( i, n )
    t = u(i)
    u(i) = u(j)
    u(j) = t
  end do

  return
end
subroutine unicycle_rank ( n, u, rank )

!*****************************************************************************80
!
!! unicycle_rank() computes the rank of a unicycle.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 June 2012
!
!  Author:
!
!    John Burkardt.
!
!  Parameters:
!
!    Input, integer N, the order of the unicycle.
!
!    Input, integer U(N), a unicycle in sequence form.
!
!    Output, integer RANK, the rank of the unicycle.
!
  implicit none

  integer n

  integer p(n-1)
  integer rank
  integer u(n)

  p(1:n-1) = u(2:n) - 1

  call perm_lex_rank ( n - 1, p, rank )

  return
end
subroutine unicycle_unrank ( n, rank, u )

!*****************************************************************************80
!
!! unicycle_unrank() "unranks" a unicycle.
!
!  Discussion:
!
!    That is, given a rank, it computes the corresponding unicycle.
!
!    The value of the rank should be between 0 and (N-1)!-1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 June 2012
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Dennis Stanton, Dennis White,
!    Constructive Combinatorics,
!    Springer, 1986,
!    ISBN: 0387963472,
!    LC: QA164.S79.
!
!  Parameters:
!
!    Input, integer N, the number of elements in the set.
!
!    Input, integer RANK, the desired rank of the permutation.
!
!    Output, integer U(N), the unicycle.
!
  implicit none

  integer n

  integer p(n-1)
  integer rank
  integer u(n)

  call perm_lex_unrank ( n - 1, rank, p )

  u(1) = 1
  u(2:n) = p(1:n-1) + 1

  return
end
