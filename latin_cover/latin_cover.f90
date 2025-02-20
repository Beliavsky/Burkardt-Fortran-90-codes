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
function i4_uniform ( a, b, seed )

!*****************************************************************************80
!
!! I4_UNIFORM returns a scaled pseudorandom I4.
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
!    12 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer A, B, the limits of the interval.
!
!    Input/output, integer SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer I4_UNIFORM, a number between A and B.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer a
  integer b
  integer, parameter :: i4_huge = 2147483647
  integer i4_uniform
  integer k
  real ( kind = rk ) r
  integer seed
  integer value

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r = real ( seed, kind = rk ) * 4.656612875E-10
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

  i4_uniform = value

  return
end
function i4_wrap ( ival, ilo, ihi )

!*****************************************************************************80
!
!! I4_WRAP forces an I4 to lie between given limits by wrapping.
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
subroutine i4block_print ( l, m, n, a, title )

!*****************************************************************************80
!
!! I4BLOCK_PRINT prints an I4BLOCK.
!
!  Discussion:
!
!    An I4BLOCK is a 3D array of I4 values.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer L, M, N, the dimensions of the block.
!
!    Input, integer A(L,M,N), the matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer l
  integer m
  integer n

  integer a(l,m,n)
  integer i
  integer j
  integer jhi
  integer jlo
  integer k
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  do k = 1, n

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  K = ', k

    do jlo = 1, m, 10
      jhi = min ( jlo + 10 - 1, m )
      write ( *, '(a)' ) ' '
      write ( *, '(8x,a2,10(2x,i6))' ) 'J:', ( j, j = jlo, jhi )
      write ( *, '(7x,a2)' ) 'I:'
      do i = 1, l
        write ( *, '(2x,i6,a1,1x,10(2x,i6))' ) i, ':', a(i,jlo:jhi,k)
      end do
    end do

  end do

  return
end
subroutine i4mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! I4MAT_PRINT prints an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4 values.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 June 2003
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
!    Input, integer A(M,N), the matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer m
  integer n

  integer a(m,n)
  integer ihi
  integer ilo
  integer jhi
  integer jlo
  character ( len = * ) title

  ilo = 1
  ihi = m
  jlo = 1
  jhi = n

  call i4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

  return
end
subroutine i4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! I4MAT_PRINT_SOME prints some of an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4 values.
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
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer A(M,N), an M by N matrix to be printed.
!
!    Input, integer ILO, JLO, the first row and column to print.
!
!    Input, integer IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: incx = 10
  integer m
  integer n

  integer a(m,n)
  character ( len = 8 )  ctemp(incx)
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

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)'
    return
  end if

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8)' ) j
    end do

    write ( *, '(''  Col '',10a8)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        write ( ctemp(j2), '(i8)' ) a(i,j)

      end do

      write ( *, '(i5,a,10a8)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine latin_cover ( n, p, a )

!*****************************************************************************80
!
!! LATIN_COVER returns a 2D Latin Square Covering.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    24 June 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!
!    Input, integer P(N), a permutation which describes the
!    first Latin square.
!
!    Output, integer A(N,N), the Latin cover.  A(I,J) = K
!    means that (I,J) is one element of the K-th Latin square.
!
  implicit none

  integer n

  integer a(n,n)
  integer i
  integer i4_wrap
  integer ik
  integer k
  integer p(n)

  call perm_check ( n, p )

  do i = 1, n
    do k = 1, n
      ik = i4_wrap ( i + k - 1, 1, n )
      a(i,p(ik)) = k
    end do
  end do

  return
end
subroutine latin_cover_2d ( n, p, a )

!*****************************************************************************80
!
!! LATIN_COVER_2D returns a 2D Latin Square Covering.
!
!  Discussion:
!
!    This procedure has a chance of being extended to M dimensions.
!
!    A basic solution is computed, and the user is permitted to permute
!    both the I and J coordinates.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    24 June 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!
!    Input, integer P(2,N), permutations to be applied
!    to the spatial dimensions.
!
!    Output, integer A(N,N), the Latin cover.  A(I,J) = K
!    means that (I,J) is one element of the K-th Latin square.
!
  implicit none

  integer n

  integer a(n,n)
  integer b(n,n)
  integer :: base = 1
  integer i
  integer i4_wrap
  integer j
  integer p(2,n)

  call perm_check ( n, p(1,1:n) )
  call perm_check ( n, p(2,1:n) )
!
!  Set up the basic solution.
!
  do i = 1, n
    do j = 1, n
      a(i,j) = i4_wrap ( i - j + base, 0 + base, n - 1 + base )
    end do
  end do
!
!  Apply permutation to dimension I.
!
  do i = 1, n
    b(p(1,i),1:n) = a(i,1:n) 
  end do
!
!  Apply permutation to dimension J.
!
  do j = 1, n
    a(1:n,p(2,j)) = b(1:n,j) 
  end do

  return
end
subroutine latin_cover_3d ( n, p, a )

!*****************************************************************************80
!
!! LATIN_COVER_3D returns a 3D Latin Square Covering.
!
!  Discussion:
!
!    A basic solution is computed, and the user is permitted to permute
!    I, J and K coordinates.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    24 June 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!
!    Input, integer P(3,N), permutations to be applied
!    to the spatial dimensions.
!
!    Output, integer A(N,N,N), the Latin cover.  A(I,J,K) = L
!    means that (I,J,K) is one element of the L-th Latin square.
!
  implicit none

  integer n

  integer a(n,n,n)
  integer b(n,n,n)
  integer i
  integer i4_wrap
  integer ik
  integer j
  integer jk
  integer k
  integer p(3,n)

  call perm_check ( n, p(1,1:n) )
  call perm_check ( n, p(2,1:n) )
  call perm_check ( n, p(3,1:n) )
!
!  Set up the basic solution.
!
  do i = 1, n
    do j = 1, n
      do k = 1, n
        ik = i4_wrap ( i + 1 - k, 1, n )
        jk = i4_wrap ( j + 1 - k, 1, n )
        b(i,j,k) = ik + ( jk - 1 ) * n
      end do
    end do
  end do
!
!  Apply permutation to dimension I.
!
  do i = 1, n
    a(p(1,i),1:n,1:n) = b(i,1:n,1:n) 
  end do
!
!  Apply permutation to dimension J.
!
  do j = 1, n
    b(1:n,p(2,j),1:n) = a(1:n,j,1:n) 
  end do
!
!  Apply permutation to dimension K.
!
  do k = 1, n
    a(1:n,1:n,p(3,k)) = b(1:n,1:n,k) 
  end do

  return
end
subroutine perm_check ( n, p )

!*****************************************************************************80
!
!! PERM_CHECK checks that a vector represents a permutation.
!
!  Discussion:
!
!    The routine verifies that each of the integers from 1
!    to N occurs among the N entries of the permutation.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries.
!
!    Input, integer P(N), the permutation, in standard index form.
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
      write ( *, '(a)' ) 'PERM_CHECK - Fatal error!'
      write ( *, '(a)' ) '  The input array does not represent'
      write ( *, '(a)' ) '  a proper permutation.  In particular, the'
      write ( *, '(a,i8)' ) '  array is missing the value ', ierror
      stop
    end if

  end do

  return
end
subroutine perm_print ( n, p, title )

!*****************************************************************************80
!
!! PERM_PRINT prints a permutation.
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
subroutine perm_uniform ( n, base, seed, p )

!*****************************************************************************80
!
!! PERM_UNIFORM selects a random permutation of N objects.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 November 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer N, the number of objects to be permuted.
!
!    Input, integer BASE, is 0 for a 0-based permutation and 1 for 
!    a 1-based permutation.
!
!    Input/output, integer SEED, a seed for the random
!    number generator.
!
!    Output, integer P(N), the permutation.  P(I) is the "new"
!    location of the object originally at I.
!
  implicit none

  integer n

  integer base
  integer i
  integer i4_uniform
  integer j
  integer k
  integer p(n)
  integer seed

  do i = 1, n
    p(i) = ( i - 1 ) + base
  end do

  do i = 1, n
    j = i4_uniform ( i, n, seed )
    k    = p(i)
    p(i) = p(j)
    p(j) = k
  end do

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
