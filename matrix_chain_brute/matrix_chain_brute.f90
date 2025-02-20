function catalan_number ( n )

!*****************************************************************************80
!
!! catalan_number() computes the Nth Catalan number.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 April 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the index of the Catalan number.
!
!  Output:
!
!    integer CATALAN_NUMBER: the value of the Catalan number.
!
  implicit none

  integer c
  integer catalan_number
  integer i
  integer n

  if ( n < 0 ) then
    catalan_number = 0
    return
  end if

  c = 1
!
!  The extra parentheses ensure that the integer division is
!  done AFTER the integer multiplication.
!
  do i = 1, n
    c = ( c * 2 * ( 2 * i - 1 ) ) / ( i + 1 )
  end do

  catalan_number = c

  return
end
function i4_factorial ( n )

!*****************************************************************************80
!
!! i4_factorial() computes the factorial of an I4.
!
!  Discussion:
!
!    Factorial ( N ) = Product ( 1 <= I <= N ) I
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 January 1999
!
!  Input:
!
!    integer N, the argument of the factorial function.
!    If N is less than 1, I4_FACTORIAL is returned as 1.
!
!  Output:
!
!    integer I4_FACTORIAL, the factorial of N.
!
  implicit none

  integer i
  integer i4_factorial
  integer n

  i4_factorial = 1

  do i = 1, n
    i4_factorial = i4_factorial * i
  end do

  return
end
subroutine matrix_chain_brute ( n_mats, dims, cost, p )

!*****************************************************************************80
!
!! matrix_chain_brute() finds the lowest cost to form a multiple matrix product.
!
!  Discussion:
!
!    This code represents a brute force approach.
!
!    An "efficient" brute force approach would only go through every
!    possible parenthesization of the multiplication, and a cost of
!    catalan(n_mats-1).  But it's not clear to me how to generate a sequence
!    of parentheses and properly interpret them as a multiplication ordering.
!
!    Instead, we rely on the fact that, to multiply N matrices, there are
!    N-1 multiplications to carry out, in any order.  This means that,
!    if we do some careful bookkeeping, we have N-1 choices for the first
!    multiplication, N-2 for the second, and a total of factorial(N-1)    
!    distinct cases to consider.  Some of these cases collapse to the
!    same paretheses representation, but we don't care.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 June 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer n_mats: the number of matrices in the product.
!
!    integer dims(n_mats+1): matrix dimension information.  Matrix A(i)
!    has dimensions dims(i) x dims(i+1).  All entries must be positive.
!
!  Output:
!
!    integer cost: the minimal cost, in terms of scalar multiplications,
!    for the optimal ordering of the matrix multiplications.
!
!    integer p(n_mats-1): defines the order of the multiplications.
!
  implicit none

  integer n_mats

  integer cost
  integer dims(n_mats+1)
  integer i
  integer n_dims
  integer n_mults
  integer p(n_mats-1)
  integer pivot_sequence_num
  integer this_cost
  integer this_p(n_mats-1)

  n_dims = n_mats + 1
  n_mults = n_mats - 1
!
!  Deal with stupidity.
!
  if ( n_mats == 1 ) then
    cost = 0
    return
  end if

  if ( any ( dims <= 0 ) ) then
    cost = 0
    p(1:n_mults) = 1
    return
  end if
!
!  Initialize the output.
!
  cost = huge ( 1 )
  do i = 1, n_mults
    p(i) = n_mults + 1 - i
  end do
!
!  Prepare to loop over all orderings.
!
  this_p(1:n_mults) = p(1:n_mults)
  call pivot_sequence_enum ( n_mults, pivot_sequence_num )

  do i = 1, pivot_sequence_num

    call pivot_sequence_successor ( n_mults, this_p )
    call pivot_sequence_to_matrix_chain_cost ( n_mats, this_p, dims, this_cost )

    if ( this_cost < cost ) then
      cost = this_cost
      p(1:n_mults) = this_p(1:n_mults)
    end if

  end do

  return
end
subroutine pivot_sequence_check ( n, t, check )

!*****************************************************************************80
!
!! pivot_sequence_check() checks a pivot sequence.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 June 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer n: the number of pivot steps.
!
!    integer t(n): a pivot sequence.
!
!  Output:
!
!    logical CHECK, error flag.
!    true, T is a pivot sequence.
!    false, T is not a legal pivot sequence.
!
  implicit none

  integer n

  logical check
  integer i
  integer t(n)
  logical verbose

  verbose = .true.
  check = .true.

  do i = 1, n

    if ( t(i) <= 0 ) then
      if ( verbose ) then
        write ( *, '(a)' ) ''
        write ( *, '(a)' ) 'pivot_sequence_check(): Fatal error!'
        write ( *, '(a,i2,a,i8,a)' ) '  t(', i, ') = ', t(i), ' <= 0'
      end if
      check = .false.
      return
    else if ( n + 1 - i < t(i) ) then
      if ( verbose ) then
        write ( *, '(a)' ) ''
        write ( *, '(a)' ) 'pivot_seq_check - Fatal error!'
        write ( *, '(a,i2,a,i2,a,i8)' ) &
          '  n + 1 - i = ', n + 1 - i, ' < t(', i, ') = ', t(i)
      end if
      check = .false.
      return
    end if

  end do

  return
end
subroutine pivot_sequence_enum ( n, pivot_sequence_num )

!*****************************************************************************80
!
!! pivot_sequence_enum() enumerates pivot sequences.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 June 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer n: the number of pivot steps.
!
!  Output:
!
!    integer pivot_sequence_num: the number of pivot sequences of n steps.
!
  implicit none

  integer i4_factorial
  integer n
  integer pivot_sequence_num

  pivot_sequence_num = i4_factorial ( n )

  return
end
subroutine pivot_sequence_successor ( n, t )

!*****************************************************************************80
!
!! pivot_sequence_successor() computes the pivot sequence successor.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 June 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer n: the number of pivot steps.
!
!    integer t(n): the previous pivot sequence.
!    To initiate the routine, call with t=linspace(n,1,n).
!
!  Output:
!
!    integer t(n): the lexical successor of the input.
!
  implicit none

  integer n

  logical check
  integer i
  integer last
  integer t(n)
!
!  Check.
!
  call pivot_sequence_check ( n, t, check )

  if ( .not. check ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'pivot_sequence_successor(): Fatal error!'
    write ( *, '(a)' ) '  The input array is illegal.'
    stop
  end if
!
!  Find last entry that is less than its maximum value.
!
  last = 0
  do i = n, 1, -1
    if ( t(i) < n + 1 - i ) then
      last = i
      exit
    end if
  end do

  if ( last == 0 ) then
    do i = 1, n
      t(i) = 1
    end do
  else
    t(last) = t(last) + 1
    t(last+1:n) = 1
  end if

  return
end
subroutine pivot_sequence_to_matrix_chain_cost ( n_mats, p, dims, cost )

!*****************************************************************************80
!
!! pivot_sequence_to_matrix_chain_cost() evaluates a particular matrix chain product.
!
!  Discussion:
!
!    There are n_mats matrices to multiply.
!    There are n_mats+1 dimensions to consider.
!    There are n_mats-1 matrix multiplications to carry out.
!
!    The cost, in terms of multiplications of pairs of real numbers,
!    of a multiplying a single pair of matrices A and B,
!    of dimensions LxM and MxN, is L*M*N.
!
!    The cost of multiplying a sequence of N matrices A*B*C*...*Z
!    will in general vary, depending on the order in which the N-1
!    multiplications are carried out.
!
!    This function assumes that a particular multiplication ordering 
!    has been specified in the pivot sequence p(), and returns the
!    corresponding cost.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 June 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer n_mats: the number of matrices to multiply.
!
!    integer p(n_mats-1): the pivot sequence.
!    On the i-th step of the multiplication procedure, there are
!    n-i multiplications to choose from, and we choose the p(i) one.
!    It must be the case that, for each i, 1 <= p(i) <= n - i
!
!    integer dims(n_mats+1): the matrix dimensions.
!    For 1 <= i < n, the i-th matrix has dimension dim(i) x dim(i+1).
!
!  Output:
!
!    integer cost: the matrix multiplication cost.
!
  implicit none

  integer n_mats

  integer cost
  integer d1
  integer d2
  integer d3
  integer dims(n_mats+1)
  integer dims2(n_mats+1)
  integer i
  integer j
  integer n_dims
  integer n_dims2
  integer n_mults
  integer p(n_mats-1)

  n_dims = n_mats + 1
  n_mults = n_mats - 1

  cost = 0
  n_dims2 = n_dims
  dims2(1:n_dims) = dims(1:n_dims)
!
!  Carry out n_mults - 1 multiplications.
!  On step I, we carry out multiplication J.
!  Multiplication J involves 
!    dim(j)dim(j+1) * dim(j+1)dim(j+2) => dim(j)dim(j+2)
!  at an additional cost of dim(j) * dim(j+1) * dim(j+2),
!  while removing dim(j+1) from the list of dimensions.
!
  do i = 1, n_mults
    j = p(i)
    d1 = dims2(j)
    d2 = dims2(j+1)
    d3 = dims2(j+2)
    cost = cost + d1 * d2 * d3
    dims2(j+1:n_dims2-1) = dims2(j+2:n_dims2)
    n_dims2 = n_dims2 - 1
  end do

  return
end
