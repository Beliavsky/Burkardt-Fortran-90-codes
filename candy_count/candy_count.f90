subroutine candy_count_box ( c, l, m, n, counts )

!*****************************************************************************80
!
!! candy_count_box() counts candy types in an LxMxN box.
!
!  Discussion:
!
!    We are given a box of candy containing C distinct types, and
!    asked to report how many of each type there are.
!
!    In this case, the box is an L by M by N matrix represented as A(I,J,K).  
!
!    The box entry in row I, column J, level K will store candy type C, where
!    A(i,j,k) = mod ( I + J + K - 3, C ) + 1.
!    (If we start I, J, K and C indexing at 0, this simplifies to
!    mod(I+J+K,C)! )
!
!    The effect of this numbering scheme is that the candy type is
!    constant along diagonal lines and sheets.
!
!    The task is to determine, for a given set of C, L, M and N,
!    the number of candies of each type.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 June 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer c: the number of types of candy.
!
!    integer l, m, n: the number of rows, columns, and levels in the candy box.
!
!  Output:
!
!    integer counts(c): the number of each type of candy in the box.
!
  implicit none

  integer c

  integer counts(c)
  integer counts_2d(c)
  integer ic
  integer ic2
  integer k
  integer l
  integer lf
  integer lr
  integer m
  integer n

  counts(1:c) = 0
!
!  Do the computation for the first level:
!
  call candy_count_matrix ( c, m, n, counts_2d )
!
!  L = LF * C + LR
!  Note that Fortran will automatically return a rounded down integer
!  result for (l/c).
!
  lf = ( l / c )
  lr = l - lf * c
!
!  Part of the box can be regarded as LF repetitions 
!  of a stack of sheets for which the item in the (K,1,1) position is 
!  successively 1, 2, ..., C.
!
!  Warning: Fortran mod(a,b) is perfectly capable of returning a NEGATIVE value.
!  Hence, we replace "mod(ic-k,c)" by "mod(c+ic-k,c)".
!
  do k = 1, c
    do ic = 1, c
      ic2 = mod ( c + ic - k, c ) + 1
      counts(ic) = counts(ic) + lf * counts_2d(ic2)
    end do
  end do
!
!  Now there are 0 <= LR < C more sheets, for which the item in the (K,1,1) 
!  position successively is 1, 2, ..., LR.
!
  do k = 1, lr
    do ic = 1, c
      ic2 = mod ( c + ic - k, c ) + 1
      counts(ic) = counts(ic) + counts_2d(ic2)
    end do
  end do

  return
end
subroutine candy_count_box_sum ( c, l, m, n, counts )

!*****************************************************************************80
!
!! candy_count_box_sum() counts candy types in an LxMxN box.
!
!  Discussion:
!
!    This function is a "stupid" version of candy_count_box().  Instead of
!    using a formula, it sets up the candy box, and counts items one by one.
!    It is useful as a check of the intelligent version.
!
!    We are given a box of candy containing C distinct types, and
!    asked to report how many of each type there are.
!
!    In this case, the box is an L by M by N matrix represented as A(I,J,K).  
!
!    The box entry in row I, column J, level K will store candy type C, where
!    A(i,j,k) = mod ( I + J + K - 3, C ) + 1.
!    (If we start I, J, K and C indexing at 0, this simplifies to
!    mod(I+J+K,C)! )
!
!    The effect of this numbering scheme is that the candy type is
!    constant along diagonal lines and sheets.
!
!    The task is to determine, for a given set of C, L, M and N,
!    the number of candies of each type.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 June 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer c: the number of types of candy.
!
!    integer l, m, n: the number of rows, columns, and levels in the candy box.
!
!  Output:
!
!    integer counts(c): the number of each type of candy in the box.
!
  implicit none

  integer c

  integer cijk
  integer counts(c)
  integer i
  integer j
  integer k
  integer l
  integer m
  integer n

  counts(1:c) = 0

  do i = 1, l
    do j = 1, m
      do k = 1, n
        cijk = mod ( i + j + k - 3, c ) + 1
        counts(cijk) = counts(cijk) + 1
      end do
    end do
  end do

  return
end
subroutine candy_count_matrix ( c, m, n, counts )

!*****************************************************************************80
!
!! candy_count_matrix() counts candy types in a matrix.
!
!  Discussion:
!
!    We are given a box of candy containing C distinct types, and
!    asked to report how many of each type there are.
!
!    In this case, the box is an M by N matrix represented as A(I,J).  
!    We will assume that rows and columns are indexed by I and J.
!
!    The box entry in row I, column J will store candy type C, where
!    A(i,j) = mod ( I + J - 2, C ) + 1.
!    (If we start I, J and C indexing at 0, this simplifies to
!    mod(I+J,C)! )
!
!    The effect of this numbering scheme is that the candy type is
!    constant along diagonal lines.
!
!    The task is to determine, for a given set of C, M and N,
!    the number of candies of each type.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 June 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer c: the number of types of candy.
!
!    integer m, n: the number of rows and columns in the candy box.
!
!  Output:
!
!    integer counts(c): the number of each type of candy in the box.
!
  implicit none

  integer c

  integer a
  integer b
  integer counts(c)
  integer i
  integer m
  integer n
  integer nmin

  a = mod ( m, c )
  b = mod ( n, c )

  nmin = ( ( m * n - a * b ) / c )

  do i = 1, c
    if ( i <= a + b - 1 - c ) then
      counts(i) = nmin + a + b - c
    else if ( i < min ( a, b ) ) then
      counts(i) = nmin + i
    else if ( i <= max ( a, b ) ) then
      counts(i) = nmin + min ( a, b )
    else if ( i < a + b - 1 ) then
      counts(i) = nmin + a + b - i
    else
      counts(i) = nmin
    end if
  end do

  return
end
subroutine candy_count_matrix_sum ( c, m, n, counts )

!*****************************************************************************80
!
!! candy_count_matrix_sum() counts candy types in a matrix.
!
!  Discussion:
!
!    This function is a "stupid" version of candy_count_matrix().  Instead of
!    using a formula, it sets up the candy box, and counts items one by one.
!    It is useful as a check of the intelligent version.
!
!    We are given a box of candy containing C distinct types, and
!    asked to report how many of each type there are.
!
!    In this case, the box is an M by N matrix represented as A(I,J).  
!    We will assume that rows and columns are indexed by I and J.
!
!    The box entry in row I, column J will store candy type C, where
!    A(i,j) = mod ( I + J - 2, C ) + 1.
!    (If we start I, J and C indexing at 0, this simplifies to
!    mod(I+J,C)! )
!
!    The effect of this numbering scheme is that the candy type is
!    constant along diagonal lines.
!
!    The task is to determine, for a given set of C, M and N,
!    the number of candies of each type.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 June 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer c: the number of types of candy.
!
!    integer m, n: the number of rows and columns in the candy box.
!
!  Output:
!
!    integer counts(c): the number of each type of candy in the box.
!
  implicit none

  integer c

  integer cij
  integer counts(c)
  integer i
  integer j
  integer m
  integer n

  counts(1:c) = 0

  do i = 1, m
    do j = 1, n
      cij = mod ( i + j - 2, c ) + 1
      counts(cij) = counts(cij) + 1
    end do
  end do

  return
end
subroutine candy_count_vector ( c, n, counts )

!*****************************************************************************80
!
!! candy_count_vector() counts candy types in a vector.
!
!  Discussion:
!
!    We are given a box of candy containing C distinct types, and
!    asked to report how many of each type there are.
!
!    In this case, the box is a vector, which can hold N items,
!    indexed 1 through N.
!
!    The candy types occur in sequence, that is, the first item is
!    candy type 1, item 2 is type 2, item C is type C.  Then
!    the types repeat, so item C+1 is type 1, and so on.
!
!    This version of the problem is simple to understand and solve.
!    The problem is more interesting if the candy is stored in a 
!    rectangular array, or a 3-dimensional box.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 June 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer c: the number of types of candy.
!
!    integer n: the number of columns in the candy box.
!
!  Output:
!
!    integer counts(c): the number of each type of candy in the box.
!
  implicit none

  integer c

  integer counts(c)
  integer n
  integer n_min
  integer n_rem
!
!  N = N_MIN * C + N_REM
!
  n_min = ( n / c )
!
!  Everybody gets at least N_MIN.
!
  counts(1:c) = n_min
!
!  There are N_REM candies remaining.
!  Give 1 each to types 1 through N_REM.
!
  n_rem = n - n_min * c

  counts(1:n_rem) = counts(1:n_rem) + 1

  return
end
subroutine candy_count_vector_sum ( c, n, counts )

!*****************************************************************************80
!
!! candy_count_vector_sum() counts candy types in a vector.
!
!  Discussion:
!
!    This function is a "stupid" version of candy_count_vector().  Instead of
!    using a formula, it sets up the candy box, and counts items one by one.
!    It is useful as a check of the intelligent version.
!
!    We are given a box of candy containing C distinct types, and
!    asked to report how many of each type there are.
!
!    In this case, the box is a vector, which can hold N items,
!    indexed 1 through N.
!
!    The candy types occur in sequence, that is, the first item is
!    candy type 1, item 2 is type 2, item C is type C.  Then
!    the types repeat, so item C+1 is type 1, and so on.
!
!    This version of the problem is simple to understand and solve.
!    The problem is more interesting if the candy is stored in a 
!    rectangular array, or a 3-dimensional box.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 June 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer c: the number of types of candy.
!
!    integer n: the number of columns in the candy box.
!
!  Output:
!
!    integer counts(c): the number of each type of candy in the box.
!
  implicit none

  integer c

  integer ci
  integer counts(c)
  integer i
  integer n

  counts(1:c) = 0

  do i = 1, n
    ci = mod ( i - 1, c ) + 1
    counts(ci) = counts(ci) + 1
  end do

  return
end

