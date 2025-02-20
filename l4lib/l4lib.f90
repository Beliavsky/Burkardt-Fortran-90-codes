subroutine boolean_to_string ( b, s )

!*****************************************************************************80
!
!! boolean_to_string() returns 'True' or 'False' depending on a boolean value.
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
!    boolean B: the boolean value.
!
!  Output:
!
!    string S: the corresponding string.
!
  implicit none

  logical b
  character * ( * ) s

  if ( b ) then
    s = 'True'
  else
    s = 'False'
  end if

  return
end
function i4_to_l4 ( i4 )

!*****************************************************************************80
!
!! I4_TO_L4() converts an I4 to an L4.
!
!  Discussion:
!
!    0 is FALSE, and anything else if TRUE.
!
!    An I4 is an integer value.
!    An L4 is a logical value.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I4, an integer.
!
!    Output, logical I4_TO_L4, the logical value of I4.
!
  implicit none

  integer i4
  logical i4_to_l4
  logical value

  value = ( i4 /= 0 )

  i4_to_l4 = value

  return
end
subroutine i4_to_l4vec ( i4, n, l4vec )

!*****************************************************************************80
!
!! i4_to_l4vec() converts an I4 into an L4VEC.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 November 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I4, the integer.
!
!    Input, integer N, the dimension of the vector.
!
!    Output, logical L4VEC(N), the vector.
!
  implicit none

  integer n

  logical l4vec(n)
  integer i
  integer i4
  integer i4_copy

  i4_copy = i4

  do i = n, 1, -1
    if ( mod ( i4_copy, 2 ) == 0 ) then
      l4vec(i) = .false.
    else
      l4vec(i) = .true.
    end if
    i4_copy = i4_copy / 2
  end do

  return
end
subroutine l4_to_i4 ( l4, i4 )

!*****************************************************************************80
!
!! l4_to_i4() converts an L4 to an I4.
!
!  Discussion:
!
!    0 is FALSE, and anything else if TRUE.
!
!    An I4 is an integer value.
!    An L4 is a logical value.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 November 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, logical L4, a logical value.
!
!    Output, integer I4, the integer value of L.
!
  implicit none

  integer i4
  logical l4

  if ( l4 ) then
    i4 = 1
  else
    i4 = 0
  end if

  return
end
subroutine l4_to_s ( l4, s )

!*****************************************************************************80
!
!! l4_to_s() converts an L4 to a string.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 August 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, logical L4, a logical value.
!
!    Output, character ( len = 5 ) S, the string.
!
  implicit none

  logical l4
  character ( len = 5 ) s

  if ( l4 ) then
    s = 'True '
  else
    s = 'False'
  end if

  return
end
function l4_uniform ( seed )

!*****************************************************************************80
!
!! l4_uniform() returns a pseudorandom L4.
!
!  Discussion:
!
!    An L4 is a LOGICAL value.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 December 2007
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
!    Input/output, integer SEED, the "seed" value, which should
!    NOT be 0.  On output, SEED has been updated.
!
!    Output, logical L4_UNIFORM, a pseudorandom logical value.
!
  implicit none

  integer, parameter :: i4_huge      = 2147483647
  integer, parameter :: i4_huge_half = 1073741823
  integer k
  logical l4_uniform
  integer seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'L4_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  l4_uniform = ( i4_huge_half < seed )

  return
end
function l4_xnor ( l1, l2 )

!*****************************************************************************80
!
!! l4_xnor() returns the complement exclusive OR of two L4's.
!
!  Discussion:
!
!    An L4 is a logical value.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 December 2017
!
!  Author:
!
!   John Burkardt
!
!  Parameters:
!
!    Input, logical L1, L2, two values.
!
!    Output, logical L4_XNOR, the complement exclusive OR 
!    of L1 and L2.
!
  implicit none

  logical l1
  logical l2
  logical l4_xnor
  logical value

  if ( l1 .and. l2 ) then
    value = .true.
  else if ( ( .not. l1 ) .and. ( .not. l2 ) ) then
    value = .true.
  else
    value = .false.
  end if

  l4_xnor = value

  return
end
function l4_xor ( l1, l2 )

!*****************************************************************************80
!
!! l4_xor() returns the exclusive OR of two L4's.
!
!  Discussion:
!
!    An L4 is a logical value.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 May 2015
!
!  Author:
!
!   John Burkardt
!
!  Parameters:
!
!    Input, logical L1, L2, two values.
!
!    Output, logical L4_XOR, the exclusive OR of L1 and L2.
!
  implicit none

  logical l1
  logical l2
  logical l4_xor
  logical value

  if ( l1 .and. ( .not. l2 ) ) then
    value = .true.
  else if ( ( .not. l1 ) .and. l2 ) then
    value = .true.
  else
    value = .false.
  end if

  l4_xor = value

  return
end
subroutine l4mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! l4mat_print() prints an L4MAT.
!
!  Discussion:
!
!    An L4MAT is an array of L4 values.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows in A.
!
!    Input, integer N, the number of columns in A.
!
!    Input, logical A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer m
  integer n

  logical a(m,n)
  character ( len = * ) title

  call l4mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine l4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! l4mat_print_some() prints some of an L4MAT.
!
!  Discussion:
!
!    An L4MAT is an array of L4 values.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, logical A(M,N), an M by N matrix to be printed.
!
!    Input, integer ILO, JLO, the first row and column to print.
!
!    Input, integer IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: incx = 35
  integer m
  integer n

  logical a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer i
  integer i2hi
  integer i2lo
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    if ( 100 <= j2hi ) then
      do j = j2lo, j2hi
        j2 = j + 1 - j2lo
        write ( ctemp(j2), '(1x,i1)' ) j / 100
      end do
      write ( *, '(''      '',35a2)' ) ctemp(1:inc)
    end if

    if ( 10 <= j2hi ) then
      do j = j2lo, j2hi
        j2 = j + 1 - j2lo
        write ( ctemp(j2), '(1x,i1)' ) mod ( j / 10, 10 )
      end do
      write ( *, '(''      '',35a2)' ) ctemp(1:inc)
    end if

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(1x,i1)' ) mod ( j, 10 )
    end do
    write ( *, '(''  Col '',35a2)' ) ctemp(1:inc)

    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi
      write ( *, '(i5,a1,35(1x,l1))' ) i, ':', a(i,j2lo:j2hi)
    end do

  end do

  return
end
subroutine l4mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! l4mat_transpose_print() prints an L4MAT, transposed.
!
!  Discussion:
!
!    An L4MAT is an array of L4 values.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, logical A(M,N), an M by N matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer m
  integer n

  logical a(m,n)
  character ( len = * ) title

  call l4mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine l4mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! l4mat_transpose_print_some() prints some of an L4MAT, transposed.
!
!  Discussion:
!
!    An L4MAT is an array of L4 values.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, logical A(M,N), an M by N matrix to be printed.
!
!    Input, integer ILO, JLO, the first row and column to print.
!
!    Input, integer IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: incx = 35
  integer m
  integer n

  logical a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer i
  integer i2
  integer i2hi
  integer i2lo
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    if ( 100 <= i2hi ) then
      do i = i2lo, i2hi
        i2 = i + 1 - i2lo
        write ( ctemp(i2), '(1x,i1)' ) i / 100
      end do
      write ( *, '(''      '',35a2)' ) ctemp(1:inc)
    end if

    if ( 10 <= i2hi ) then
      do i = i2lo, i2hi
        i2 = i + 1 - i2lo
        write ( ctemp(i2), '(1x,i1)' ) mod ( i / 10, 10 )
      end do
      write ( *, '(''      '',35a2)' ) ctemp(1:inc)
    end if

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(1x,i1)' ) mod ( i, 10 )
    end do
    write ( *, '(''  Row '',35a2)' ) ctemp(1:inc)

    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi
      write ( *, '(i5,a1,35(1x,l1))' ) j, ':', a(i2lo:i2hi,j)
    end do

  end do

  return
end
subroutine l4mat_uniform ( m, n, seed, lmat )

!*****************************************************************************80
!
!! l4mat_uniform() returns a pseudorandom L4MAT.
!
!  Discussion:
!
!    An L4MAT is a two dimensional array of L4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 December 2007
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
!    Input, integer M, N, the order of the matrix.
!
!    Input/output, integer SEED, the "seed" value, which should
!    NOT be 0.  On output, SEED has been updated.
!
!    Output, logical LMAT(M,N), a pseudorandom logical matrix.
!
  implicit none

  integer m
  integer n

  integer, parameter :: i4_huge      = 2147483647
  integer, parameter :: i4_huge_half = 1073741823
  integer i
  integer j
  integer k
  logical lmat(m,n)
  integer seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'L4MAT_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  do j = 1, n

    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge
      end if

      lmat(i,j) = ( i4_huge_half < seed )

    end do

  end do

  return
end
subroutine l4vec_next ( n, l4vec )

!*****************************************************************************80
!
!! l4vec_next() generates the next logical vector.
!
!  Discussion:
!
!    In the following discussion, we will let '0' stand for FALSE and
!    '1' for TRUE.
!
!    The logical vectors have the order
!
!      (0,0,...,0),
!      (0,0,...,1),
!      ...
!      (1,1,...,1)
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
!    31 May 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of the vectors.
!
!    Input/output, logical L4VEC(N), on output, the successor to the
!    input vector.
!
  implicit none

  integer n

  integer i
  logical l4vec(n)

  do i = n, 1, -1

    if ( .not. l4vec(i) ) then
      l4vec(i) = .true.
      return
    end if

    l4vec(i) = .false.

  end do

  return
end
subroutine l4vec_print ( n, a, title )

!*****************************************************************************80
!
!! l4vec_print() prints an L4VEC.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 May 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, logical A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer  n

  logical a(n)
  integer i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,a,1x,l1)' ) i, ':', a(i)
  end do

  return
end
subroutine l4vec_print_some ( n, a, i_lo, i_hi, title )

!*****************************************************************************80
!
!! l4vec_print_some() prints "some" of an L4VEC.
!
!  Discussion:
!
!    An L4VEC is a vector of logical values.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries of the vector.
!
!    Input, logical A(N), the vector to be printed.
!
!    Input, integer I_LO, I_HI, the first and last indices to 
!    print. The routine expects 1 <= I_LO <= I_HI <= N.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer n

  logical a(n)
  integer i
  integer i_hi
  integer i_lo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = max ( i_lo, 1 ), min ( i_hi, n )
    write ( *, '(2x,i8,a,1x,l1)' ) i, ':', a(i)
  end do

  return
end
subroutine l4vec_uniform ( n, seed, lvec )

!*****************************************************************************80
!
!! l4vec_uniform() returns a pseudorandom L4VEC.
!
!  Discussion:
!
!    An L4VEC is a vector of L4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 December 2007
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
!    Input, integer N, the order of the vector.
!
!    Input/output, integer SEED, the "seed" value, which should
!    NOT be 0.  On output, SEED has been updated.
!
!    Output, logical LVEC(N), a pseudorandom logical vector.
!
  implicit none

  integer n

  integer, parameter :: i4_huge      = 2147483647
  integer, parameter :: i4_huge_half = 1073741823
  integer i
  integer k
  logical lvec(n)
  integer seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'L4VEC_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    lvec(i) = ( i4_huge_half < seed )

  end do

  return
end
function s_to_l4 ( s )

!*****************************************************************************80
!
!! s_to_l4() reads a logical value from a string.
!
!  Discussion:
!
!    There are several ways of representing logical data that this routine
!    recognizes:
!
!      False   True
!      -----   ----
!
!      0       1
!      F       T
!      f       t
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 December 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be read.
!
!    Output, logical S_TO_L4, the logical value read from the string.
!
  implicit none

  integer i
  character ( len = * ) s
  integer s_length
  logical s_to_l4

  s_length = len_trim ( s )

  do i = 1, s_length
    if ( s(i:i) == '0' .or. s(i:i) == 'F' .or. s(i:i) == 'f' ) then
      s_to_l4 = .false.
      return
    else if ( s(i:i) == '1' .or. s(i:i) == 'T' .or. s(i:i) == 't' ) then
      s_to_l4 = .true.
      return
    end if
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'S_TO_L4 - Fatal error!'
  write ( *, '(a)' ) '  Input text did not contain logical data.'

  stop 1
end

