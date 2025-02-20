subroutine bvec_add ( n, bvec1, bvec2, bvec3 )

!*****************************************************************************80
!
!! bvec_add() adds two (signed) binary vectors.
!
!  Discussion:
!
!    A BVEC is a vector of binary digits representing an integer.
!
!    BVEC(1) is 0 for positive values and 1 for negative values, which
!    are stored in 2's complement form.
!
!    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
!    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
!    so that printing the digits in order gives the binary form of the number.
!
!  Example:
!
!    N = 5
!
!      BVEC1       +   BVEC2       =   BVEC3
!
!    ( 0 0 0 0 1 ) + ( 0 0 0 1 1 ) = ( 0 0 1 0 0 )
!
!              1   +           3   =           4
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the length of the vectors.
!
!    integer BVEC1(N), BVEC2(N), the vectors to be added.
!
!  Output:
!
!    integer BVEC3(N), the sum of the two input vectors.
!
  implicit none

  integer n

  integer, parameter :: base = 2
  integer bvec1(n)
  integer bvec2(n)
  integer bvec3(n)
  integer i
  logical overflow

  overflow = .false.

  bvec3(1:n) = bvec1(1:n) + bvec2(1:n)

  do i = n, 1, -1

    do while ( base <= bvec3(i) )

      bvec3(i) = bvec3(i) - base

      if ( 1 < i ) then
        bvec3(i-1) = bvec3(i-1) + 1
      else
        overflow = .true.
      end if

    end do

  end do

  return
end
subroutine bvec_and ( n, bvec1, bvec2, bvec3 )

!*****************************************************************************80
!
!! bvec_and() computes the AND of two binary vectors.
!
!  Discussion:
!
!    A BVEC is a vector of binary digits representing an integer.
!
!    BVEC(1) is 0 for positive values and 1 for negative values, which
!    are stored in 2's complement form.
!
!    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
!    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
!    so that printing the digits in order gives the binary form of the number.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the length of the vectors.
!
!    integer BVEC1(N), BVEC2(N), the binary vectors.
!
!    integer BVEC3(N), the AND of the two vectors.
!
  implicit none

  integer n

  integer bvec1(n)
  integer bvec2(n)
  integer bvec3(n)

  bvec3(1:n) = min ( bvec1(1:n), bvec2(1:n) )

  return
end
subroutine bvec_check ( n, bvec, ierror )

!*****************************************************************************80
!
!! bvec_check() checks a binary vector.
!
!  Discussion:
!
!    A BVEC is a vector of binary digits representing an integer.
!
!    BVEC(1) is 0 for positive values and 1 for negative values, which
!    are stored in 2's complement form.
!
!    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
!    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
!    so that printing the digits in order gives the binary form of the number.
!
!    The only check made is that the entries are all 0 or 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the length of the vectors.
!
!    integer BVEC(N), the vector to be checked.
!
!  Output:
!
!    integer IERROR, is nonzero if an error occurred.
!
  implicit none

  integer n

  integer, parameter :: base = 2
  integer bvec(n)
  integer i
  integer ierror

  ierror = 0

  do i = 1, n
    if ( bvec(i) < 0 .or. base <= bvec(i) ) then
      ierror = i
      return
    end if
  end do

  return
end
subroutine bvec_complement2 ( n, bvec1, bvec2 )

!*****************************************************************************80
!
!! bvec_complement2() computes the two's complement of a binary vector.
!
!  Discussion:
!
!    A BVEC is a vector of binary digits representing an integer.
!
!    BVEC(1) is 0 for positive values and 1 for negative values, which
!    are stored in 2's complement form.
!
!    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
!    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
!    so that printing the digits in order gives the binary form of the number.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the length of the vectors.
!
!    integer BVEC1(N), the vector to be complemented.
!
!  Output:
!
!    integer BVEC2(N), the two's complemented vector.
!
  implicit none

  integer n

  integer, parameter :: base = 2
  integer bvec1(n)
  integer bvec2(n)
  integer bvec3(n)
  integer bvec4(n)

  bvec3(1:n) = ( base - 1 ) - bvec1(1:n)

  bvec4(1:n-1) = 0
  bvec4(n) = 1

  call bvec_add ( n, bvec3, bvec4, bvec2 )

  return
end
function bvec_enum ( n )

!*****************************************************************************80
!
!! bvec_enum() enumerates the binary vectors of length N.
!
!  Discussion:
!
!    A BVEC is a vector of binary digits representing an integer.
!
!    BVEC(1) is 0 for positive values and 1 for negative values, which
!    are stored in 2's complement form.
!
!    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
!    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
!    so that printing the digits in order gives the binary form of the number.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the length of the vectors.
!
!  Output:
!
!    integer BVEC_ENUM, the number of binary vectors.
!
  implicit none

  integer bvec_enum
  integer n
  integer value

  value = 2 ** n

  bvec_enum = value

  return
end
subroutine bvec_mul ( n, bvec1, bvec2, bvec3 )

!*****************************************************************************80
!
!! bvec_mul() computes the product of two binary vectors.
!
!  Discussion:
!
!    A BVEC is a vector of binary digits representing an integer.
!
!    BVEC(1) is 0 for positive values and 1 for negative values, which
!    are stored in 2's complement form.
!
!    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
!    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
!    so that printing the digits in order gives the binary form of the number.
!
!    Since the user may want to make calls like
!
!      call bvec_mul ( n, bvec1, bvec1, bvec3 )
!    or even
!      call bvec_mul ( n, bvec1, bvec1, bvec1 )
!
!    we need to copy the arguments, work on them, and then copy out the result.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the length of the vectors.
!
!    integer BVEC1(N), BVEC2(N), the vectors to
!    be multiplied.
!
!  Output:
!
!    integer BVEC3(N), the product of the two
!    input vectors.
!
  implicit none

  integer n

  integer, parameter :: base = 2
  integer carry
  integer bvec1(n)
  integer bvec2(n)
  integer bvec3(n)
  integer bveca(n)
  integer bvecb(n)
  integer bvecc(n)
  integer i
  integer product_sign
!
!  Copy the input.
!
  bveca(1:n) = bvec1(1:n)
  bvecb(1:n) = bvec2(1:n)
!
!  Record the sign of the product.
!  Make the factors positive.
!
  product_sign = 1

  if ( bveca(1) /= 0 ) then
    product_sign = - product_sign
    call bvec_complement2 ( n, bveca, bveca )
  end if

  if ( bvecb(1) /= 0 ) then
    product_sign = - product_sign
    call bvec_complement2 ( n, bvecb, bvecb )
  end if

  bvecc(1:n) = 0
!
!  Multiply.
!
  do i = 2, n
    bvecc(2:n+2-i) = bvecc(2:n+2-i) + bveca(n+2-i) * bvecb(i:n)
  end do
!
!  Take care of carries.
!
  do i = n, 2, -1

    carry = bvecc(i) / base
    bvecc(i) = bvecc(i) - carry * base
!
!  Unlike the case of BVEC_ADD, we do NOT allow carries into
!  the first position when multiplying.
!
    if ( 2 < i ) then
      bvecc(i-1) = bvecc(i-1) + carry
    end if

  end do
!
!  Take care of the sign of the product.
!
  if ( product_sign < 0 ) then
    call bvec_complement2 ( n, bvecc, bvecc )
  end if
!
!  Copy the output.
!
  bvec3(1:n) = bvecc(1:n)

  return
end
subroutine bvec_next ( n, bvec )

!*****************************************************************************80
!
!! bvec_next() generates the next BVEC.
!
!  Discussion:
!
!    A BVEC is a binary vector, an N vector whose entries are 0 or 1.
!
!    The vectors are produced in the order:
!
!    (0,0,...,0),
!    (0,0,...,1),
!    ...
!    (1,1,...,1)
!
!    and the "next" vector after (1,1,...,1) is (0,0,...,0).  That is,
!    we allow wrap around.
!
!  Example:
!
!    N = 3
!
!    Input      Output
!    -----      ------
!    0 0 0  =>  0 0 1
!    0 0 1  =>  0 1 0
!    0 1 0  =>  0 1 1
!    0 1 1  =>  1 0 0
!    1 0 0  =>  1 0 1
!    1 0 1  =>  1 1 0
!    1 1 0  =>  1 1 1
!    1 1 1  =>  0 0 0
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the dimension of the vectors.
!
!    integer BVEC(N), a binary vector.
!
!  Output:
!
!    integer BVEC(N): the successor to the input vector.
!
  implicit none

  integer n

  integer bvec(n)
  integer i

  do i = n, 1, -1

    if ( bvec(i) == 0 ) then
      bvec(i) = 1
      return
    end if

    bvec(i) = 0

  end do

  return
end
subroutine bvec_next_grlex ( n, bvec )

!*****************************************************************************80
!
!! bvec_next_grlex() generates the next binary vector in GRLEX order.
!
!  Discussion:
!
!    N = 3
!
!    Input      Output
!    -----      ------
!    0 0 0  =>  0 0 1
!    0 0 1  =>  0 1 0
!    0 1 0  =>  1 0 0
!    1 0 0  =>  0 1 1
!    0 1 1  =>  1 0 1
!    1 0 1  =>  1 1 0
!    1 1 0  =>  1 1 1
!    1 1 1  =>  0 0 0
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the dimension.
!
!    integer BVEC(N), the binary vector whose 
!    successor is desired.
!
!  Output:
!
!    integer BVEC(N), the successor to the input vector.
!
  implicit none

  integer n

  integer bvec(n)
  integer i
  integer o
  integer s
  integer z
!
!  Initialize locations of 0 and 1.
!
  if ( bvec(1) == 0 ) then
    z = 1
    o = 0
  else
    z = 0
    o = 1
  end if
!
!  Moving from right to left, search for a "1", preceded by a "0".
!
  do i = n, 2, -1
    if ( bvec(i) == 1 ) then
      o = i
      if ( bvec(i-1) == 0 ) then
        z = i - 1
        exit
      end if
    end if
  end do
!
!  BVEC = 0
!
  if ( o == 0 ) then
    bvec(n) = 1
!
!  01 never occurs.  So for sure, B(1) = 1.
!
  else if ( z == 0 ) then
    s = sum ( bvec(1:n) )
    if ( s == n ) then
      bvec(1:n) = 0
    else
      bvec(1:n-s-1) = 0
      bvec(n-s:n) = 1
    end if
!
!  Found the rightmost "01" string.
!  Replace it by "10".
!  Shift following 1's to the right.
!
  else
    bvec(z) = 1
    bvec(o) = 0
    s = sum ( bvec(o+1:n) )
    bvec(o+1:n-s) = 0
    bvec(n+1-s:n) = 1
  end if

  return
end
subroutine bvec_not ( n, bvec1, bvec2 )

!*****************************************************************************80
!
!! bvec_not() "negates" or takes the 1's complement of a binary vector.
!
!  Discussion:
!
!    A BVEC is a vector of binary digits representing an integer.
!
!    BVEC(1) is 0 for positive values and 1 for negative values, which
!    are stored in 2's complement form.
!
!    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
!    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
!    so that printing the digits in order gives the binary form of the number.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the length of the vectors.
!
!    integer BVEC1(N), the vector to be negated.
!
!  Output:
!
!    integer BVEC2(N), the negated vector.
!
  implicit none

  integer n

  integer, parameter :: base = 2
  integer bvec1(n)
  integer bvec2(n)

  bvec2(1:n) = ( base - 1 ) - bvec1(1:n)

  return
end
subroutine bvec_or ( n, bvec1, bvec2, bvec3 )

!*****************************************************************************80
!
!! bvec_or() computes the inclusive OR of two binary vectors.
!
!  Discussion:
!
!    A BVEC is a vector of binary digits representing an integer.
!
!    BVEC(1) is 0 for positive values and 1 for negative values, which
!    are stored in 2's complement form.
!
!    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
!    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
!    so that printing the digits in order gives the binary form of the number.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the length of the vectors.
!
!    integer BVEC1(N), BVEC2(N), the binary vectors.
!
!  Output:
!
!    integer BVEC3(N), the inclusive OR of the two vectors.
!
  implicit none

  integer n

  integer bvec1(n)
  integer bvec2(n)
  integer bvec3(n)

  bvec3(1:n) = max ( bvec1(1:n), bvec2(1:n) )

  return
end
subroutine bvec_print ( n, bvec, title )

!*****************************************************************************80
!
!! bvec_print() prints a BVEC, with an optional title.
!
!  Discussion:
!
!    A BVEC is a vector of binary digits representing an integer.
!
!    BVEC(1) is 0 for positive values and 1 for negative values, which
!    are stored in 2's complement form.
!
!    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
!    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
!    so that printing the digits in order gives the binary form of the number.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of components of the vector.
!
!    integer BVEC(N), the vector to be printed.
!
!    character ( len = * ) TITLE, a title.
!
  implicit none

  integer n

  integer bvec(n)
  integer ihi
  integer ilo
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do ilo = 1, n, 70
    ihi = min ( ilo + 70 - 1, n )
    write ( *, '(2x,70i1)' ) bvec(ilo:ihi)
  end do

  return
end
subroutine bvec_reverse ( n, bvec1, bvec2 )

!*****************************************************************************80
!
!! bvec_reverse() reverses a binary vector.
!
!  Discussion:
!
!    A BVEC is a vector of binary digits representing an integer.
!
!    BVEC(1) is 0 for positive values and 1 for negative values, which
!    are stored in 2's complement form.
!
!    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
!    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
!    so that printing the digits in order gives the binary form of the number.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the length of the vectors.
!
!    integer BVEC1(N), the vector to be reversed.
!
!  Output:
!
!    integer BVEC2(N), the reversed vector.
!
  implicit none

  integer n

  integer bvec1(n)
  integer bvec2(n)

  bvec2(1:n) = bvec1(n:1:-1)

  return
end
subroutine bvec_sub ( n, bvec1, bvec2, bvec3 )

!*****************************************************************************80
!
!! bvec_sub() subtracts two binary vectors.
!
!  Discussion:
!
!    A BVEC is a vector of binary digits representing an integer.
!
!    BVEC(1) is 0 for positive values and 1 for negative values, which
!    are stored in 2's complement form.
!
!    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
!    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
!    so that printing the digits in order gives the binary form of the number.
!
!  Example:
!
!    N = 4
!
!    BVEC1         BVEC2         BVEC3
!    -------       -------       -------
!    0 1 0 0   -   0 0 0 1   =   0 0 1 1
!          4             1             3
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the length of the vectors.
!
!    integer BVEC1(N), BVEC2(N), the vectors to
!    be subtracted.
!
!  Output:
!
!    integer BVEC3(N), the value of BVEC1 - BVEC2.
!
  implicit none

  integer n

  integer bvec1(n)
  integer bvec2(n)
  integer bvec3(n)
  integer bvec4(n)

  call bvec_complement2 ( n, bvec2, bvec4 )

  call bvec_add ( n, bvec1, bvec4, bvec3 )

  return
end
subroutine bvec_to_i4 ( n, bvec, i4 )

!*****************************************************************************80
!
!! bvec_to_i4() makes an integer from a (signed) binary vector.
!
!  Discussion:
!
!    A BVEC is a vector of binary digits representing an integer.
!
!    BVEC(1) is 0 for positive values and 1 for negative values, which
!    are stored in 2's complement form.
!
!    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
!    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
!    so that printing the digits in order gives the binary form of the number.
!
!  Example:
!
!         BVEC   binary  I
!    ----------  -----  --
!    1  2  3  4
!    ----------
!    0  0  0  1       1  1
!    0  0  1  0      10  2
!    1  1  0  0    -100 -4
!    0  1  0  0     100  4
!    1  0  0  1    -111 -9
!    1  1  1  1      -0  0
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the dimension of the vector.
!
!    integer BVEC(N), the binary representation.
!
!  Output:
!
!    integer I4, the integer.
!
  implicit none

  integer n

  integer, parameter :: base = 2
  integer bvec(n)
  integer bvec2(n)
  integer i
  integer i_sign
  integer i4

  bvec2(1:n) = bvec(1:n)

  if ( bvec2(1) == base - 1 ) then
    i_sign = -1
    bvec2(1) = 0
    call bvec_complement2 ( n - 1, bvec2(2:n), bvec2(2:n) )
  else
    i_sign = 1
  end if

  i4 = 0
  do i = 2, n
    i4 = base * i4 + bvec2(i)
  end do

  i4 = i_sign * i4

  return
end
subroutine bvec_uniform ( n, bvec )

!*****************************************************************************80
!
!! bvec_uniform() returns a pseudorandom BVEC.
!
!  Discussion:
!
!    An BVEC is a vector of binary (0/1) values, representing an integer.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 September 2021
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
!    integer N, the order of the vector.
!
!  Output:
!
!    integer BVEC(N), a pseudorandom binary vector.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer bvec(n)
  integer i
  real ( kind = rk ) r

  do i = 1, n

    call random_number ( harvest = r )

    if ( r < 0.5D+00 ) then
      bvec(i) = 0
    else
      bvec(i) = 1
    end if

  end do

  return
end
subroutine bvec_xor ( n, bvec1, bvec2, bvec3 )

!*****************************************************************************80
!
!! bvec_xor() computes the exclusive OR of two binary vectors.
!
!  Discussion:
!
!    A BVEC is a vector of binary digits representing an integer.
!
!    BVEC(1) is 0 for positive values and 1 for negative values, which
!    are stored in 2's complement form.
!
!    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
!    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
!    so that printing the digits in order gives the binary form of the number.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the length of the vectors.
!
!    integer BVEC1(N), BVEC2(N), the binary vectors
!    to be XOR'ed.
!
!  Output:
!
!    integer BVEC3(N), the exclusive OR of the two vectors.
!
  implicit none

  integer n

  integer bvec1(n)
  integer bvec2(n)
  integer bvec3(n)

  bvec3(1:n) = mod ( bvec1(1:n) + bvec2(1:n), 2 )

  return
end
function i4_bclr ( i4, pos )

!*****************************************************************************80
!
!! i4_bclr() returns a copy of an I4 in which the POS-th bit is set to 0.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Military Standard 1753,
!    FORTRAN, DoD Supplement To American National Standard X3.9-1978,
!    9 November 1978.
!
!  Input:
!
!    integer I4, the integer to be tested.
!
!    integer POS, the bit position, between 0 and 31.
!
!  Output:
!
!    integer I4_BCLR, a copy of I4, but with the POS-th bit
!    set to 0.
!
  implicit none

  integer i4
  integer i4_bclr
  integer, parameter :: i4_huge = 2147483647
  integer j
  integer k
  integer pos
  integer sub
  integer value

  value = i4

  if ( pos < 0 ) then

  else if ( pos < 31 ) then

    sub = 1

    if ( 0 <= i4 ) then
      j = i4
    else
      j = ( i4_huge + i4 ) + 1
    end if

    do k = 1, pos
      j = j / 2
      sub = sub * 2
    end do

    if ( mod ( j, 2 ) == 1 ) then
      value = i4 - sub
    end if

  else if ( pos == 31 ) then

    if ( i4 < 0 ) then
      value = ( i4_huge + i4 ) + 1
    end if

  else if ( 31 < pos ) then

    value = i4

  end if

  i4_bclr = value

  return
end
function i4_bset ( i4, pos )

!*****************************************************************************80
!
!! i4_bset() returns a copy of an I4 in which the POS-th bit is set to 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Military Standard 1753,
!    FORTRAN, DoD Supplement To American National Standard X3.9-1978,
!    9 November 1978.
!
!  Input:
!
!    integer I4, the integer to be tested.
!
!    integer POS, the bit position, between 0 and 31.
!
!  Output:
!
!    integer I4_BSET, a copy of I4, but with the POS-th bit
!    set to 1.
!
  implicit none

  integer add
  integer i4
  integer i4_bset
  integer, parameter :: i4_huge = 2147483647
  integer j
  integer k
  integer pos
  integer value

  value = i4

  if ( pos < 0 ) then

  else if ( pos < 31 ) then

    add = 1

    if ( 0 <= i4 ) then
      j = i4
    else
      j = ( i4_huge + i4 ) + 1
    end if

    do k = 1, pos
      j = j / 2
      add = add * 2
    end do

    if ( mod ( j, 2 ) == 0 ) then
      value = i4 + add
    end if

  else if ( pos == 31 ) then

    if ( 0 < i4 ) then
      value = - ( i4_huge - i4 ) - 1
    end if

  else if ( 31 < pos ) then

    value = i4

  end if

  i4_bset = value

  return
end
function i4_btest ( i4, pos )

!*****************************************************************************80
!
!! i4_btest() returns TRUE if the POS-th bit of an I4 is 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Military Standard 1753,
!    FORTRAN, DoD Supplement To American National Standard X3.9-1978,
!    9 November 1978.
!
!  Input:
!
!    integer I4, the integer to be tested.
!
!    integer POS, the bit position, between 0 and 31.
!
!  Output:
!
!    logical I4_BTEST, is TRUE if the POS-th bit of I4 is 1.
!
  implicit none

  integer i4
  logical i4_btest
  integer, parameter :: i4_huge = 2147483647
  integer j
  integer k
  integer pos

  if ( pos < 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_BTEST - Fatal error!'
    write ( *, '(a)' ) '  POS < 0.'
    stop 1

  else if ( pos < 31 ) then

    if ( 0 <= i4 ) then
      j = i4
    else
      j = ( i4_huge + i4 ) + 1
    end if

    do k = 1, pos
      j = j / 2
    end do

    if ( mod ( j, 2 ) == 0 ) then
      i4_btest = .false.
    else
      i4_btest = .true.
    end if

  else if ( pos == 31 ) then

    if ( i4 < 0 ) then
      i4_btest = .true.
    else
      i4_btest = .false.
    end if

  else if ( 31 < pos ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_BTEST - Fatal error!'
    write ( *, '(a)' ) '  31 < POS.'
    stop 1

  end if

  return
end
subroutine i4_to_bvec ( i4, n, bvec )

!*****************************************************************************80
!
!! i4_to_bvec() makes a signed binary vector from an I4.
!
!  Discussion:
!
!    A BVEC is a vector of binary digits representing an integer.
!
!    BVEC(1) is 0 for positive values and 1 for negative values, which
!    are stored in 2's complement form.
!
!    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
!    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
!    so that printing the digits in order gives the binary form of the number.
!
!    Negative values have a two's complement operation applied.
!
!    To guarantee that there will be enough space for any
!    value of I, it would be necessary to set N = 32.
!
!  Example:
!
!    I4       BVEC         binary
!    --  ----------------  ------
!     1  1  0  0  0  0  1      1
!     2  0  0  0  0  1  0     10
!     3  0  0  0  0  1  1     11
!     4  0  0  0  1  0  0    100
!     9  0  0  1  0  0  1   1001
!    -9  1  1  0  1  1  1  -1001 = 110111 (2's complement)
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer I4, an integer to be represented.
!
!    integer N, the dimension of the vector.
!
!  Output:
!
!    integer BVEC(N), the signed binary representation.
!
  implicit none

  integer n

  integer, parameter :: base = 2
  integer bvec(n)
  integer i4
  integer i4_copy
  integer j

  i4_copy = abs ( i4 )

  do j = n, 2, - 1

    bvec(j) = mod ( i4_copy, base )

    i4_copy = i4_copy / base

  end do

  bvec(1) = 0

  if ( i4 < 0 ) then
    call bvec_complement2 ( n, bvec, bvec )
  end if

  return
end
function i4_uniform_ab ( a, b )

!*****************************************************************************80
!
!! i4_uniform_ab() returns a scaled pseudorandom I4.
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
!    03 September 2021
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
  integer value

  call random_number ( harvest = r )
!
!  Scale R to lie between A-0.5 and B+0.5.
!
  r = ( 1.0E+00 - r ) * ( real ( min ( a, b ), kind = rk ) - 0.5E+00 ) & 
    +             r   * ( real ( max ( a, b ), kind = rk ) + 0.5E+00 )
!
!  Use rounding to convert R to an integer between A and B.
!
  value = nint ( r )

  value = max ( value, min ( a, b ) )
  value = min ( value, max ( a, b ) )

  i4_uniform_ab = value

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! timestamp() prints the current YMDHMS date as a time stamp.
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
!    03 September 2021
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
