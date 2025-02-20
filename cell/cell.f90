subroutine i4cvv_iget ( mn, a, m, roff, i, j, aij )

!*****************************************************************************80
!
!! I4CVV_IGET gets item J from row I in an I4CVV.
!
!  Discussion:
!
!    An I4CVV is a "vector of vectors" of I4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 December 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MN, the size of the cell array.
!
!    Input, integer A(MN), the cell array.
!
!    Input, integer M, the number of rows in the array.
!
!    Input, integer ROFF(M+1), the row offsets.
!
!    Input, integer I, the row of the item.
!    1 <= I <= M.
!
!    Input, integer J, the column of the item.
!    1 <= J <= NR(I).
!
!    Output, integer AIJ, the value of item A(I,J).
!
  implicit none

  integer m
  integer mn

  integer a(mn)
  integer aij
  integer i
  integer j
  integer k
  integer roff(m+1)

  k = roff(i) + j
  aij = a(k)

  return
end
subroutine i4cvv_iinc ( mn, a, m, roff, i, j, daij )

!*****************************************************************************80
!
!! I4CVV_IINC increments item J from row I in an I4CVV.
!
!  Discussion:
!
!    An I4CVV is a "vector of vectors" of I4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 December 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MN, the size of the cell array.
!
!    Input/output, integer A(MN), the cell array.
!
!    Input, integer M, the number of rows in the array.
!
!    Input, integer ROFF(M+1), the row offsets.
!
!    Input, integer I, the row of the item.
!    1 <= I <= M.
!
!    Input, integer J, the column of the item.
!    1 <= J <= NR(I).
!
!    Input, integer DAIJ, the increment to item A(I,J).
!
  implicit none

  integer m
  integer mn

  integer a(mn)
  integer daij
  integer i
  integer j
  integer k
  integer roff(m+1)

  k = roff(i) + j
  a(k) = a(k) + daij

  return
end
subroutine i4cvv_indx ( m, roff, i, j, k )

!*****************************************************************************80
!
!! I4CVV_INDX gets the index K of item J from row I in an I4CVV.
!
!  Discussion:
!
!    An I4CVV is a "vector of vectors" of I4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 December 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows in the array.
!
!    Input, integer ROFF(M+1), the row offsets.
!
!    Input, integer I, the row of the item.
!    1 <= I <= M.
!
!    Input, integer J, the column of the item.
!    1 <= J <= NR(I).
!
!    Output, integer K, the index of item J in row I.
!
  implicit none

  integer m

  integer i
  integer j
  integer k
  integer roff(m+1)

  k = roff(i) + j

  return
end
subroutine i4cvv_iset ( mn, a, m, roff, i, j, aij )

!*****************************************************************************80
!
!! I4CVV_ISET sets item J from row I in an I4CVV.
!
!  Discussion:
!
!    An I4CVV is a "vector of vectors" of I4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 December 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MN, the size of the cell array.
!
!    Input/output, integer A(MN), the cell array.
!
!    Input, integer M, the number of rows in the array.
!
!    Input, integer ROFF(M+1), the row offsets.
!
!    Input, integer I, the row of the item.
!    1 <= I <= M.
!
!    Input, integer J, the column of the item.
!    1 <= J <= NR(I).
!
!    Input, integer AIJ, the new value of item A(I,J).
!
  implicit none

  integer m
  integer mn

  integer a(mn)
  integer aij
  integer i
  integer j
  integer k
  integer roff(m+1)

  k = roff(i) + j
  a(k) = aij

  return
end
subroutine i4cvv_nget ( mn, a, m, roff, nn, in, jn, vn )

!*****************************************************************************80
!
!! I4CVV_NGET gets N items JN(*) from rows IN(*) in an I4CVV.
!
!  Discussion:
!
!    An I4CVV is a "vector of vectors" of I4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 December 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MN, the size of the cell array.
!
!    Input, integer A(MN), the cell array.
!
!    Input, integer M, the number of rows in the array.
!
!    Input, integer ROFF(M+1), the row offsets.
!
!    Input, integer NN, the number of items.
!
!    Input, integer IN(NN), the rows of the items.
!    1 <= IN(*) <= M.
!
!    Input, integer JN(NN), the columns of the items.
!    1 <= JN(*) <= NR(IN(*)).
!
!    Output, integer VN(NN), the value of items A(IN(*),JN(*)).
!
  implicit none

  integer m
  integer mn
  integer nn

  integer a(mn)
  integer i
  integer in(nn)
  integer jn(nn)
  integer k
  integer roff(m+1)
  integer vn(nn)

  do i = 1, nn
    k = roff(in(i)) + jn(i)
    vn(i) = a(k)
  end do

  return
end
subroutine i4cvv_ninc ( mn, a, m, roff, nn, in, jn, dvn )

!*****************************************************************************80
!
!! I4CVV_NINC increments items JN(*) from rows IN(*) in an I4CVV.
!
!  Discussion:
!
!    An I4CVV is a "vector of vectors" of I4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 December 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MN, the size of the cell array.
!
!    Input/output, integer A(MN), the cell array.
!
!    Input, integer M, the number of rows in the array.
!
!    Input, integer ROFF(M+1), the row offsets.
!
!    Input, integer NN, the number of items.
!
!    Input, integer IN(NN), the rows of the items.
!    1 <= IN(*) <= M.
!
!    Input, integer JN(NN), the columns of the items.
!    1 <= JN(*) <= NR(IN(*)).
!
!    Input, integer DVN(NN), the increments of 
!    items A(IN(*),JN(*)).
!
  implicit none

  integer m
  integer mn
  integer nn

  integer a(mn)
  integer dvn(nn)
  integer i
  integer in(nn)
  integer jn(nn)
  integer k
  integer roff(m+1)

  do i = 1, nn
    k = roff(in(i)) + jn(i)
    a(k) = a(k) + dvn(i)
  end do

  return
end
subroutine i4cvv_nndx ( m, roff, nn, in, jn, kn )

!*****************************************************************************80
!
!! I4CVV_NNDX gets the indices KN of N items JN(*) from rows IN(*) in an I4CVV.
!
!  Discussion:
!
!    An I4CVV is a "vector of vectors" of I4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 December 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows in the array.
!
!    Input, integer ROFF(M+1), the row offsets.
!
!    Input, integer NN, the number of items.
!
!    Input, integer IN(NN), the rows of the items.
!    1 <= IN(*) <= M.
!
!    Input, integer JN(NN), the columns of the items.
!    1 <= JN(*) <= NR(IN(*)).
!
!    Output, integer KN(NN), the indices of items A(IN(*),JN(*)).
!
  implicit none

  integer m
  integer nn

  integer i
  integer in(nn)
  integer jn(nn)
  integer kn(nn)
  integer roff(m+1)

  do i = 1, nn
    kn(i) = roff(in(i)) + jn(i)
  end do

  return
end
subroutine i4cvv_nset ( mn, a, m, roff, nn, in, jn, vn )

!*****************************************************************************80
!
!! I4CVV_NSET sets items JN(*) from rows IN(*) in an I4CVV.
!
!  Discussion:
!
!    An I4CVV is a "vector of vectors" of I4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 December 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MN, the size of the cell array.
!
!    Input/output, integer A(MN), the cell array.
!
!    Input, integer M, the number of rows in the array.
!
!    Input, integer ROFF(M+1), the row offsets.
!
!    Input, integer NN, the number of items.
!
!    Input, integer IN(NN), the rows of the items.
!    1 <= IN(*) <= M.
!
!    Input, integer JN(NN), the columns of the items.
!    1 <= JN(*) <= NR(IN(*)).
!
!    Input, integer VN(NN), the new value of items A(IN(*),JN(*)).
!
  implicit none

  integer m
  integer mn
  integer nn

  integer a(mn)
  integer i
  integer in(nn)
  integer jn(nn)
  integer k
  integer roff(m+1)
  integer vn(nn)

  do i = 1, nn
    k = roff(in(i)) + jn(i)
    a(k) = vn(i)
  end do

  return
end
subroutine i4cvv_offset ( m, nr, roff )

!*****************************************************************************80
!
!! I4CVV_OFFSET determines the row offsets of an I4CVV.
!
!  Discussion:
!
!    An I4CVV is a "vector of vectors" of I4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 December 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows in the array.
!
!    Input, integer NR(M), the row sizes.
!
!    Output, integer ROFF(M+1), the row offsets.
!
  implicit none

  integer m

  integer i
  integer roff(m+1)
  integer nr(m)

  roff(1) = 0
  do i = 1, m
    roff(i+1) = roff(i) + nr(i)
  end do

  return
end
subroutine i4cvv_print ( mn, a, m, roff, title )

!*****************************************************************************80
!
!! I4CVV_PRINT prints an I4CVV.
!
!  Discussion:
!
!    An I4CVV is a "vector of vectors" of I4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 December 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MN, the size of the cell array.
!
!    Input, integer A(MN), the cell array.
!
!    Input, integer M, the number of rows in the array.
!
!    Input, integer ROFF(M+1), the row offsets.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer m
  integer mn

  integer a(mn)
  integer i
  integer k1
  integer k2
  integer khi
  integer klo
  integer roff(m+1)
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, m

    k1 = roff(i) + 1
    k2 = roff(i+1)

    do klo = k1, k2, 10
      khi = min ( klo + 10 - 1, k2 )
      if ( klo == k1 ) then
        write ( *, '(i5,2x, 10i7)' ) i, a(klo:khi)
      else
        write ( *, '(5x,2x, 10i7)' )    a(klo:khi)
      end if
    end do

  end do

  return
end
subroutine i4cvv_rget ( mn, a, m, roff, i, ai )

!*****************************************************************************80
!
!! I4CVV_RGET gets row I from an I4CVV.
!
!  Discussion:
!
!    An I4CVV is a "vector of vectors" of I4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 December 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MN, the size of the cell array.
!
!    Input, integer A(MN), the cell array.
!
!    Input, integer M, the number of rows in the array.
!
!    Input, integer ROFF(M+1), the row offsets.
!
!    Input, integer I, the row.
!    1 <= I <= M.
!
!    Output, integer AI(NR(I)), the value of A(I,*).
!
  implicit none

  integer m
  integer mn

  integer a(mn)
  integer ai(*)
  integer i
  integer k1
  integer k2
  integer nv
  integer roff(m+1)

  k1 = roff(i) + 1
  k2 = roff(i+1)
  nv = k2 + 1 - k1
  ai(1:nv) = a(k1:k2)

  return
end
subroutine i4cvv_rinc ( mn, a, m, roff, i, dai )

!*****************************************************************************80
!
!! I4CVV_RINC increments row I in an I4CVV.
!
!  Discussion:
!
!    An I4CVV is a "vector of vectors" of I4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 December 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MN, the size of the cell array.
!
!    Input/output, integer A(MN), the cell array.
!
!    Input, integer M, the number of rows in the array.
!
!    Input, integer ROFF(M+1), the row offsets.
!
!    Input, integer I, the row.
!    1 <= I <= M.
!
!    Input, integer DAI(NR(I)), the increment for A(I,*).
!
  implicit none

  integer m
  integer mn

  integer a(mn)
  integer dai(*)
  integer i
  integer k1
  integer k2
  integer nv
  integer roff(m+1)

  k1 = roff(i) + 1
  k2 = roff(i+1)
  nv = k2 + 1 - k1
  a(k1:k2) = a(k1:k2) + dai(1:nv)

  return
end
subroutine i4cvv_rndx ( m, roff, i, k )

!*****************************************************************************80
!
!! I4CVV_RNDX gets the index K of row I in an I4CVV.
!
!  Discussion:
!
!    The index of a row is the index of the first item in the row.
!
!    An I4CVV is a "vector of vectors" of I4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 December 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows in the array.
!
!    Input, integer ROFF(M+1), the row offsets.
!
!    Input, integer I, the row of the item.
!    1 <= I <= M.
!
!    Output, integer K, the index of row I.
!
  implicit none

  integer m

  integer i
  integer k
  integer roff(m+1)

  k = roff(i) + 1

  return
end
subroutine i4cvv_rset ( mn, a, m, roff, i, ai )

!*****************************************************************************80
!
!! I4CVV_RSET sets row I from an I4CVV.
!
!  Discussion:
!
!    An I4CVV is a "vector of vectors" of I4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 December 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MN, the size of the cell array.
!
!    Input/output, integer A(MN), the cell array.
!
!    Input, integer M, the number of rows in the array.
!
!    Input, integer ROFF(M+1), the row offsets.
!
!    Input, integer I, the row.
!    1 <= I <= M.
!
!    Input, integer AI(NR(I)), the new value of A(I,*).
!
  implicit none

  integer m
  integer mn

  integer a(mn)
  integer ai(*)
  integer i
  integer k1
  integer k2
  integer nv
  integer roff(m+1)

  k1 = roff(i) + 1
  k2 = roff(i+1)
  nv = k2 + 1 - k1
  a(k1:k2) = ai(1:nv)

  return
end
subroutine i4cvv_size ( m, nr, mn )

!*****************************************************************************80
!
!! I4CVV_SIZE determines the size of an I4CVV.
!
!  Discussion:
!
!    An I4CVV is a "vector of vectors" of I4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 December 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows in the array.
!
!    Input, integer NR(M), the size of each row.
!
!    Output, integer MN, the size of the cell array.
!
  implicit none

  integer m

  integer mn
  integer nr(m)

  mn = sum ( nr(1:m) )

  return
end
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an I4VEC.
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
!    02 May 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, integer A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer n

  integer a(n)
  integer i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,a,2x,i12)' ) i, ':', a(i)
  end do

  return
end
subroutine i4vec_transpose_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_TRANSPOSE_PRINT prints an I4VEC "transposed".
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Example:
!
!    A = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 /)
!    TITLE = 'My vector:  '
!
!    My vector:
!
!        1    2    3    4    5
!        6    7    8    9   10
!       11
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, integer A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer n

  integer a(n)
  integer ihi
  integer ilo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do ilo = 1, n, 5
    ihi = min ( ilo + 5 - 1, n )
    write ( *, '(5i12)' ) a(ilo:ihi)
  end do

  return
end
subroutine r8cvv_iget ( mn, a, m, roff, i, j, aij )

!*****************************************************************************80
!
!! R8CVV_IGET gets item J from row I in an R8CVV.
!
!  Discussion:
!
!    An R8CVV is a "vector of vectors" of R8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MN, the size of the cell array.
!
!    Input, real ( kind = rk ) A(MN), the cell array.
!
!    Input, integer M, the number of rows in the array.
!
!    Input, integer ROFF(M+1), the row offsets.
!
!    Input, integer I, the row of the item.
!    1 <= I <= M.
!
!    Input, integer J, the column of the item.
!    1 <= J <= NR(I).
!
!    Output, real ( kind = rk ) AIJ, the value of item A(I,J).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer mn

  real ( kind = rk ) a(mn)
  real ( kind = rk ) aij
  integer i
  integer j
  integer k
  integer roff(m+1)

  k = roff(i) + j
  aij = a(k)

  return
end
subroutine r8cvv_iinc ( mn, a, m, roff, i, j, daij )

!*****************************************************************************80
!
!! R8CVV_IINC increments item J from row I in an R8CVV.
!
!  Discussion:
!
!    An R8CVV is a "vector of vectors" of R8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MN, the size of the cell array.
!
!    Input/output, real ( kind = rk ) A(MN), the cell array.
!
!    Input, integer M, the number of rows in the array.
!
!    Input, integer ROFF(M+1), the row offsets.
!
!    Input, integer I, the row of the item.
!    1 <= I <= M.
!
!    Input, integer J, the column of the item.
!    1 <= J <= NR(I).
!
!    Input, real ( kind = rk ) DAIJ, the increment to the value of item A(I,J).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer mn

  real ( kind = rk ) a(mn)
  real ( kind = rk ) daij
  integer i
  integer j
  integer k
  integer roff(m+1)

  k = roff(i) + j
  a(k) = a(k) + daij

  return
end
subroutine r8cvv_indx ( m, roff, i, j, k )

!*****************************************************************************80
!
!! R8CVV_INDX gets the index K of item J from row I in an R8CVV.
!
!  Discussion:
!
!    An R8CVV is a "vector of vectors" of R8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 December 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows in the array.
!
!    Input, integer ROFF(M+1), the row offsets.
!
!    Input, integer I, the row of the item.
!    1 <= I <= M.
!
!    Input, integer J, the column of the item.
!    1 <= J <= NR(I).
!
!    Output, integer K, the index of item J in row I.
!
  implicit none

  integer m

  integer i
  integer j
  integer k
  integer roff(m+1)

  k = roff(i) + j

  return
end
subroutine r8cvv_iset ( mn, a, m, roff, i, j, aij )

!*****************************************************************************80
!
!! R8CVV_ISET sets item J from row I in an R8CVV.
!
!  Discussion:
!
!    An R8CVV is a "vector of vectors" of R8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MN, the size of the cell array.
!
!    Input/output, real ( kind = rk ) A(MN), the cell array.
!
!    Input, integer M, the number of rows in the array.
!
!    Input, integer ROFF(M+1), the row offsets.
!
!    Input, integer I, the row of the item.
!    1 <= I <= M.
!
!    Input, integer J, the column of the item.
!    1 <= J <= NR(I).
!
!    Input, real ( kind = rk ) AIJ, the new value of item A(I,J).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer mn

  real ( kind = rk ) a(mn)
  real ( kind = rk ) aij
  integer i
  integer j
  integer k
  integer roff(m+1)

  k = roff(i) + j
  a(k) = aij

  return
end
subroutine r8cvv_nget ( mn, a, m, roff, nn, in, jn, vn )

!*****************************************************************************80
!
!! R8CVV_NGET gets N items JN(*) from rows IN(*) in an R8CVV.
!
!  Discussion:
!
!    An R8CVV is a "vector of vectors" of R8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MN, the size of the cell array.
!
!    Input, real ( kind = rk ) A(MN), the cell array.
!
!    Input, integer M, the number of rows in the array.
!
!    Input, integer ROFF(M+1), the row offsets.
!
!    Input, integer NN, the number of items.
!
!    Input, integer IN(NN), the rows of the items.
!    1 <= IN(*) <= M.
!
!    Input, integer JN(NN), the columns of the items.
!    1 <= JN(*) <= NR(IN(*)).
!
!    Output, real ( kind = rk ) VN(NN), the value of items A(IN(*),JN(*)).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer mn
  integer nn

  real ( kind = rk ) a(mn)
  integer i
  integer in(nn)
  integer jn(nn)
  integer k
  integer roff(m+1)
  real ( kind = rk ) vn(nn)

  do i = 1, nn
    k = roff(in(i)) + jn(i)
    vn(i) = a(k)
  end do

  return
end
subroutine r8cvv_ninc ( mn, a, m, roff, nn, in, jn, dvn )

!*****************************************************************************80
!
!! R8CVV_NINC increments items JN(*) from rows IN(*) in an R8CVV.
!
!  Discussion:
!
!    An R8CVV is a "vector of vectors" of R8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MN, the size of the cell array.
!
!    Input/output, real ( kind = rk ) A(MN), the cell array.
!
!    Input, integer M, the number of rows in the array.
!
!    Input, integer ROFF(M+1), the row offsets.
!
!    Input, integer NN, the number of items.
!
!    Input, integer IN(NN), the rows of the items.
!    1 <= IN(*) <= M.
!
!    Input, integer JN(NN), the columns of the items.
!    1 <= JN(*) <= NR(IN(*)).
!
!    Input, real ( kind = rk ) DVN(NN), the increments of items A(IN(*),JN(*)).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer mn
  integer nn

  real ( kind = rk ) a(mn)
  real ( kind = rk ) dvn(nn)
  integer i
  integer in(nn)
  integer jn(nn)
  integer k
  integer roff(m+1)

  do i = 1, nn
    k = roff(in(i)) + jn(i)
    a(k) = a(k) + dvn(i)
  end do

  return
end
subroutine r8cvv_nndx ( m, roff, nn, in, jn, kn )

!*****************************************************************************80
!
!! R8CVV_NNDX gets the indices KN of N items JN(*) from rows IN(*) in an R8CVV.
!
!  Discussion:
!
!    An R8CVV is a "vector of vectors" of R8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 December 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows in the array.
!
!    Input, integer ROFF(M+1), the row offsets.
!
!    Input, integer NN, the number of items.
!
!    Input, integer IN(NN), the rows of the items.
!    1 <= IN(*) <= M.
!
!    Input, integer JN(NN), the columns of the items.
!    1 <= JN(*) <= NR(IN(*)).
!
!    Output, integer KN(NN), the indices of items A(IN(*),JN(*)).
!
  implicit none

  integer m
  integer nn

  integer i
  integer in(nn)
  integer jn(nn)
  integer kn(nn)
  integer roff(m+1)

  do i = 1, nn
    kn(i) = roff(in(i)) + jn(i)
  end do

  return
end
subroutine r8cvv_nset ( mn, a, m, roff, nn, in, jn, vn )

!*****************************************************************************80
!
!! R8CVV_NSET sets items JN(*) from rows IN(*) in an R8CVV.
!
!  Discussion:
!
!    An R8CVV is a "vector of vectors" of R8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MN, the size of the cell array.
!
!    Input/output, real ( kind = rk ) A(MN), the cell array.
!
!    Input, integer M, the number of rows in the array.
!
!    Input, integer ROFF(M+1), the row offsets.
!
!    Input, integer NN, the number of items.
!
!    Input, integer IN(NN), the rows of the items.
!    1 <= IN(*) <= M.
!
!    Input, integer JN(NN), the columns of the items.
!    1 <= JN(*) <= NR(IN(*)).
!
!    Input, real ( kind = rk ) VN(NN), the new value of items A(IN(*),JN(*)).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer mn
  integer nn

  real ( kind = rk ) a(mn)
  integer i
  integer in(nn)
  integer jn(nn)
  integer k
  integer roff(m+1)
  real ( kind = rk ) vn(nn)

  do i = 1, nn
    k = roff(in(i)) + jn(i)
    a(k) = vn(i)
  end do

  return
end
subroutine r8cvv_offset ( m, nr, roff )

!*****************************************************************************80
!
!! R8CVV_OFFSET determines the row offsets of an R8CVV.
!
!  Discussion:
!
!    An R8CVV is a "vector of vectors" of R8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows in the array.
!
!    Input, integer NR(M), the row sizes.
!
!    Output, integer ROFF(M+1), the row offsets.
!
  implicit none

  integer m

  integer i
  integer roff(m+1)
  integer nr(m)

  roff(1) = 0
  do i = 1, m
    roff(i+1) = roff(i) + nr(i)
  end do

  return
end
subroutine r8cvv_print ( mn, a, m, roff, title )

!*****************************************************************************80
!
!! R8CVV_PRINT prints an R8CVV.
!
!  Discussion:
!
!    An R8CVV is a "vector of vectors" of R8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MN, the size of the cell array.
!
!    Input, real ( kind = rk ) A(MN), the cell array.
!
!    Input, integer M, the number of rows in the array.
!
!    Input, integer ROFF(M+1), the row offsets.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer mn

  real ( kind = rk ) a(mn)
  integer i
  integer k1
  integer k2
  integer khi
  integer klo
  integer roff(m+1)
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, m

    k1 = roff(i) + 1
    k2 = roff(i+1)

    do klo = k1, k2, 5
      khi = min ( klo + 5 - 1, k2 )
      if ( klo == k1 ) then
        write ( *, '(i5,2x, 5g14.6)' ) i, a(klo:khi)
      else
        write ( *, '(5x,2x, 5g14.6)' )    a(klo:khi)
      end if
    end do

  end do

  return
end
subroutine r8cvv_rget ( mn, a, m, roff, i, ai )

!*****************************************************************************80
!
!! R8CVV_RGET gets row I from an R8CVV.
!
!  Discussion:
!
!    An R8CVV is a "vector of vectors" of R8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MN, the size of the cell array.
!
!    Input, real ( kind = rk ) A(MN), the cell array.
!
!    Input, integer M, the number of rows in the array.
!
!    Input, integer ROFF(M+1), the row offsets.
!
!    Input, integer I, the row.
!    1 <= I <= M.
!
!    Output, real ( kind = rk ) AI(NR(I)), the value of A(I,*).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer mn

  real ( kind = rk ) a(mn)
  real ( kind = rk ) ai(*)
  integer i
  integer k1
  integer k2
  integer nv
  integer roff(m+1)

  k1 = roff(i) + 1
  k2 = roff(i+1)
  nv = k2 + 1 - k1
  ai(1:nv) = a(k1:k2)

  return
end
subroutine r8cvv_rinc ( mn, a, m, roff, i, dai )

!*****************************************************************************80
!
!! R8CVV_RINC increments row I in an R8CVV.
!
!  Discussion:
!
!    An R8CVV is a "vector of vectors" of R8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MN, the size of the cell array.
!
!    Input/output, real ( kind = rk ) A(MN), the cell array.
!
!    Input, integer M, the number of rows in the array.
!
!    Input, integer ROFF(M+1), the row offsets.
!
!    Input, integer I, the row.
!    1 <= I <= M.
!
!    Input, real ( kind = rk ) DAI(NR(I)), the increment for A(I,*).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer mn

  real ( kind = rk ) a(mn)
  real ( kind = rk ) dai(*)
  integer i
  integer k1
  integer k2
  integer nv
  integer roff(m+1)

  k1 = roff(i) + 1
  k2 = roff(i+1)
  nv = k2 + 1 - k1
  a(k1:k2) = a(k1:k2) + dai(1:nv)

  return
end
subroutine r8cvv_rndx ( m, roff, i, k )

!*****************************************************************************80
!
!! R8CVV_RNDX gets the index K of row I in an R8CVV.
!
!  Discussion:
!
!    The index of a row is the index of the first item in the row.
!
!    An R8CVV is a "vector of vectors" of R8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 December 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows in the array.
!
!    Input, integer ROFF(M+1), the row offsets.
!
!    Input, integer I, the row of the item.
!    1 <= I <= M.
!
!    Output, integer K, the index of row I.
!
  implicit none

  integer m

  integer i
  integer k
  integer roff(m+1)

  k = roff(i) + 1

  return
end
subroutine r8cvv_rset ( mn, a, m, roff, i, ai )

!*****************************************************************************80
!
!! R8CVV_RSET sets row I from an R8CVV.
!
!  Discussion:
!
!    An R8CVV is a "vector of vectors" of R8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MN, the size of the cell array.
!
!    Input/output, real ( kind = rk ) A(MN), the cell array.
!
!    Input, integer M, the number of rows in the array.
!
!    Input, integer ROFF(M+1), the row offsets.
!
!    Input, integer I, the row.
!    1 <= I <= M.
!
!    Input, real ( kind = rk ) AI(NR(I)), the new value of A(I,*).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer mn

  real ( kind = rk ) a(mn)
  real ( kind = rk ) ai(*)
  integer i
  integer k1
  integer k2
  integer nv
  integer roff(m+1)

  k1 = roff(i) + 1
  k2 = roff(i+1)
  nv = k2 + 1 - k1
  a(k1:k2) = ai(1:nv)

  return
end
subroutine r8cvv_size ( m, nr, mn )

!*****************************************************************************80
!
!! R8CVV_SIZE determines the size of an R8CVV.
!
!  Discussion:
!
!    An R8CVV is a "vector of vectors" of R8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows in the array.
!
!    Input, integer NR(M), the size of each row.
!
!    Output, integer MN, the size of the cell array.
!
  implicit none

  integer m

  integer mn
  integer nr(m)

  mn = sum ( nr(1:m) )

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
subroutine r8vec_transpose_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_TRANSPOSE_PRINT prints an R8VEC "transposed".
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Example:
!
!    A = (/ 1.0, 2.1, 3.2, 4.3, 5.4, 6.5, 7.6, 8.7, 9.8, 10.9, 11.0 /)
!    TITLE = 'My vector:  '
!
!    My vector:
!
!        1.0    2.1    3.2    4.3    5.4
!        6.5    7.6    8.7    9.8   10.9
!       11.0
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 November 2010
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
  integer ihi
  integer ilo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do ilo = 1, n, 5
    ihi = min ( ilo + 5 - 1, n )
    write ( *, '(5g14.6)' ) a(ilo:ihi)
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
!    06 August 2005
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
