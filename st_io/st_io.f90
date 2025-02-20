subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IUNIT, the free unit number.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer ios
  integer iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
subroutine i4_fake_use ( n )

!*****************************************************************************80
!
!! i4_fake_use() pretends to use a variable.
!
!  Discussion:
!
!    Some compilers will issue a warning if a variable is unused.
!    Sometimes there's a good reason to include a variable in a program,
!    but not to use it.  Calling this function with that variable as
!    the argument will shut the compiler up.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 April 2020
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the variable to be "used".
!
  implicit none

  integer n

  if ( n /= n ) then
    write ( *, '(a)' ) '  i4_fake_use: variable is NAN.'
  end if

  return
end
subroutine i4vec_dec ( n, a )

!*****************************************************************************80
!
!! I4VEC_DEC decrements an I4VEC.
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
!    18 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input/output, integer A(N), the vector to be decremented.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer a(n)

  a(1:n) = a(1:n) - 1

  return
end
subroutine i4vec_inc ( n, a )

!*****************************************************************************80
!
!! I4VEC_INC increments an I4VEC.
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
!    18 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input/output, integer A(N), the vector to be incremented.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer a(n)

  a(1:n) = a(1:n) + 1

  return
end
function i4vec_max ( n, a )

!*****************************************************************************80
!
!! I4VEC_MAX computes the maximum element of an I4VEC.
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
!    30 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input, integer A(N), the array.
!
!    Output, integer I4VEC_MAX, the value of the largest entry.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer a(n)
  integer i4vec_max

  i4vec_max = maxval ( a(1:n) )

  return
end
function i4vec_min ( n, a )

!*****************************************************************************80
!
!! I4VEC_MIN computes the minimum element of an I4VEC.
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
!    21 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input, integer A(N), the array.
!
!    Output, integer I4VEC_MIN, the value of the smallest entry.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer a(n)
  integer i4vec_min

  i4vec_min = minval ( a(1:n) )

  return
end
subroutine r8_fake_use ( x )

!*****************************************************************************80
!
!! r8_fake_use() pretends to use an R8 variable.
!
!  Discussion:
!
!    Some compilers will issue a warning if a variable is unused.
!    Sometimes there's a good reason to include a variable in a program,
!    but not to use it.  Calling this function with that variable as
!    the argument will shut the compiler up.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 April 2020
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) X, the variable to be "used".
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) x

  if ( x /= x ) then
    write ( *, '(a)' ) '  r8_fake_use: variable is NAN.'
  end if

  return
end
subroutine sort_heap_external ( n, indx, i, j, isgn )

!*****************************************************************************80
!
!! SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
!
!  Discussion:
!
!    The actual list of data is not passed to the routine.  Hence this
!    routine may be used to sort integers, reals, numbers, names,
!    dates, shoe sizes, and so on.  After each call, the routine asks
!    the user to compare or interchange two items, until a special
!    return value signals that the sorting is completed.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 February 2004
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf
!    FORTRAN90 version by John Burkardt
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
!    Input, integer N, the number of items to be sorted.
!
!    Input/output, integer INDX, the main communication signal.
!    The user must set INDX to 0 before the first call.
!    Thereafter, the user should not change the value of INDX until
!    the sorting is done.
!    On return, if INDX is
!      greater than 0,
!      * interchange items I and J;
!      * call again.
!      less than 0,
!      * compare items I and J;
!      * set ISGN = -1 if I < J, ISGN = +1 if J < I;
!      * call again.
!      equal to 0, the sorting is done.
!
!    Output, integer I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ISGN, results of comparison of elements
!    I and J. (Used only when the previous call returned INDX less than 0).
!    ISGN <= 0 means I is less than or equal to J;
!    0 <= ISGN means I is greater than or equal to J.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer, save :: i_save = 0
  integer indx
  integer isgn
  integer j
  integer, save :: j_save = 0
  integer, save :: k = 0
  integer, save :: k1 = 0
  integer n
  integer, save :: n1 = 0
!
!  INDX = 0: This is the first call.
!
  if ( indx == 0 ) then

    i_save = 0
    j_save = 0
    k = n / 2
    k1 = k
    n1 = n
!
!  INDX < 0: The user is returning the results of a comparison.
!
  else if ( indx < 0 ) then

    if ( indx == -2 ) then

      if ( isgn < 0 ) then
        i_save = i_save + 1
      end if

      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return

    end if

    if ( 0 < isgn ) then
      indx = 2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then

      if ( n1 == 1 ) then
        i_save = 0
        j_save = 0
        indx = 0
      else
        i_save = n1
        n1 = n1 - 1
        j_save = 1
        indx = 1
      end if

      i = i_save
      j = j_save
      return

    end if

    k = k - 1
    k1 = k
!
!  0 < INDX, the user was asked to make an interchange.
!
  else if ( indx == 1 ) then

    k1 = k

  end if

  do

    i_save = 2 * k1

    if ( i_save == n1 ) then
      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return
    else if ( i_save <= n1 ) then
      j_save = i_save + 1
      indx = -2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then
      exit
    end if

    k = k - 1
    k1 = k

  end do

  if ( n1 == 1 ) then
    i_save = 0
    j_save = 0
    indx = 0
    i = i_save
    j = j_save
  else
    i_save = n1
    n1 = n1 - 1
    j_save = 1
    indx = 1
    i = i_save
    j = j_save
  end if

  return
end
subroutine r8st_data_read ( input_filename, m, n, nst, ist, jst, ast )

!*****************************************************************************80
!
!! R8ST_DATA_READ reads the data of an R8ST file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 September 2018
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the file.
!
!    Input, integer M, the number of rows.
!
!    Input, integer N, the number of columns.
!
!    Input, integer NST, the number of nonzeros.
!
!    Output, integer IST(NST), JST(NST), the row and 
!    column indices.
!
!    Output, real ( kind = rk ) AST(NST), the nonzero values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nst

  real ( kind = rk ) aij
  real ( kind = rk ) ast(nst)
  integer i
  character ( len = * ) input_filename
  integer input_unit
  integer ios
  integer ist(nst)
  integer j
  integer jst(nst)
  integer k
  integer m
  integer n

  call get_unit ( input_unit )

  open ( file = input_filename, unit = input_unit, status = 'old', &
    iostat = ios )

  do k = 1, nst

    read ( input_unit, *, iostat = ios ) i, j, aij

    if ( ios /= 0 ) then
      exit
    end if

    if ( i < 1 .or. m < i .or. 1 < j .or. n < j ) then
      write ( *, '(a,i4,2x,i4)' ) 'ST_DATA_READ - Error - Bad (I,J) = ', i, j
    end if

    ist(k) = i
    jst(k) = j
    ast(k) = aij

  end do

  close ( unit = input_unit )

  return
end
subroutine r8st_header_print ( i_min, i_max, j_min, j_max, m, n, nst )

!*****************************************************************************80
!
!! R8ST_HEADER_PRINT prints the header of an R8ST file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I_MIN, I_MAX, the minimum and maximum rows.
!
!    Input, integer J_MIN, J_MAX, the minimum and maximum columns.
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer NST, the number of nonzeros.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i_max
  integer i_min
  integer j_max
  integer j_min
  integer m
  integer n
  integer nst

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sparse Triplet (ST) header information:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Minimum row index I_MIN = ', i_min
  write ( *, '(a,i8)' ) '  Maximum row index I_MAX = ', i_max
  write ( *, '(a,i8)' ) '  Minimum col index J_MIN = ', j_min
  write ( *, '(a,i8)' ) '  Maximum col index J_MAX = ', j_max
  write ( *, '(a,i8)' ) '  Number of rows        M = ', m
  write ( *, '(a,i8)' ) '  Number of columns     N = ', n
  write ( *, '(a,i8)' ) '  Number of nonzeros  NST = ', nst

  return
end
subroutine r8st_header_read ( input_filename, i_min, i_max, j_min, j_max, &
  m, n, nst )

!*****************************************************************************80
!
!! R8ST_HEADER_READ reads the header of an R8ST file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the file.
!
!    Output, integer I_MIN, I_MAX, the minimum and maximum rows.
!
!    Output, integer J_MIN, J_MAX, the minimum and maximum columns.
!
!    Output, integer M, N, the number of rows and columns.
!
!    Output, integer NST, the number of nonzeros.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) aij
  integer i
  integer i_max
  integer i_min
  integer, parameter :: i4_huge = 2147483647
  character ( len = * ) input_filename
  integer input_unit
  integer ios
  integer j
  integer j_max
  integer j_min
  integer m
  integer n
  integer nst

  call get_unit ( input_unit )

  open ( file = input_filename, unit = input_unit, status = 'old', &
    iostat = ios )

  nst = 0
  i_min = + i4_huge
  i_max = - i4_huge
  j_min = + i4_huge
  j_max = - i4_huge

  do

    read ( input_unit, *, iostat = ios ) i, j, aij

    if ( ios /= 0 ) then
      exit
    end if

    nst = nst + 1
    i_min = min ( i_min, i )
    i_max = max ( i_max, i )
    j_min = min ( j_min, j )
    j_max = max ( j_max, j )

  end do

  close ( unit = input_unit )

  m = i_max - i_min + 1
  n = j_max - j_min + 1

  return
end
subroutine r8st_print ( m, n, nst, ist, jst, ast, title )

!*****************************************************************************80
!
!! R8ST_PRINT prints a sparse matrix in R8ST format.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 September 2018
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows.
!
!    Input, integer N, the number of columns.
!
!    Input, integer NST, the number of elements.
!
!    Input, integer IST(NST), JST(NST), the rows and columns.
!
!    Input, real ( kind = rk ) AST(NST), the values.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nst

  real ( kind = rk ) ast(nst)
  integer ist(nst)
  integer jst(nst)
  integer k
  integer m
  integer n
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(2x,i4,a,i4,a)' ) m, ' rows by ', n, ' columns'
  write ( *, '(a)' ) '     #     I     J       A'
  write ( *, '(a)' ) '  ----  ----  ----  --------------'
  write ( *, '(a)' ) ' '

  do k = 1, nst
    write ( *, '(2x,i4,2x,i4,2x,i4,2x,g16.8)' ) k, ist(k), jst(k), ast(k)
  end do

  return
end
subroutine r8st_print_some ( i_min, i_max, j_min, j_max, nst, ist, jst, ast, &
  title )

!*****************************************************************************80
!
!! R8ST_PRINT_SOME prints some of a sparse matrix in R8ST format.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I_MIN, IMAX, the first and last rows to print.
!
!    Input, integer J_MIN, J_MAX, the first and last columns 
!    to print.
!
!    Input, integer NST, the number of elements.
!
!    Input, integer IST(NST), JST(NST), the ows and columns.
!
!    Input, real ( kind = rk ) AST(NST), the values.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nst

  real ( kind = rk ) ast(nst)
  integer i_max
  integer i_min
  integer ist(nst)
  integer j_max
  integer j_min
  integer jst(nst)
  integer k
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) '     #     I     J       A'
  write ( *, '(a)' ) '  ----  ----  ----  --------------'
  write ( *, '(a)' ) ' '

  do k = 1, nst
    if ( i_min <= ist(k) .and. ist(k) <= i_max .and. &
         j_min <= jst(k) .and. jst(k) <= j_max ) then
      write ( *, '(2x,i4,2x,i4,2x,i4,2x,g16.8)' ) k, ist(k), jst(k), ast(k)
    end if
  end do

  return
end
subroutine r8st_sort_a ( m, n, nst, ist, jst, ast )

!*****************************************************************************80
!
!! R8ST_SORT_A sorts the entries of an R8ST matrix by column.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 November 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows.
!
!    Input, integer N, the number of columns.
!
!    Input, integer NST, the number of nonzeros.
!
!    Input/output, integer IST(NST), JST(NST), the row and 
!    column indices.
!
!    Input/output, real ( kind = rk ) AST(NST), the nonzero values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nst

  real ( kind = rk ) aij
  real ( kind = rk ) ast(nst)
  integer cij
  integer i
  integer indx
  integer isgn
  integer ist(nst)
  integer j
  integer jst(nst)
  integer m
  integer n
  integer rij

  call i4_fake_use ( m )
  call i4_fake_use ( n )
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( nst, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      rij    = ist(i)
      ist(i) = ist(j)
      ist(j) = rij

      cij    = jst(i)
      jst(i) = jst(j)
      jst(j) = cij

      aij    = ast(i)
      ast(i) = ast(j)
      ast(j) = aij
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      if ( jst(i) == jst(j) ) then

        if ( ist(i) < ist(j) ) then
          isgn = - 1
        else if ( ist(i) == ist(j) ) then
          isgn = 0
        else if ( ist(j) < ist(i) ) then
          isgn = + 1
        end if

      else if ( jst(i) < jst(j) ) then

        isgn = - 1

      else if ( jst(j) < jst(i) ) then

        isgn = + 1

      end if

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine r8st_transpose ( m, n, nst, ist, jst, ast )

!*****************************************************************************80
!
!! R8ST_TRANSPOSE transposes an R8ST matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer M, the number of rows.
!
!    Input/output, integer N, the number of columns.
!
!    Input, integer NST, the number of nonzeros.
!
!    Input/output, integer IST(NST), JST(NST), the row and 
!    column indices.
!
!    Input, real ( kind = rk ) AST(NST), the nonzero values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nst

  real ( kind = rk ) ast(nst)
  integer ist(nst)
  integer jst(nst)
  integer k
  integer m
  integer n
  integer t

  call r8_fake_use ( ast(1) )

  t = m
  m = n
  n = t

  do k = 1, nst

    t      = ist(k)
    ist(k) = jst(k)
    jst(k) = t

  end do

  return
end
subroutine r8st_write ( output_filename, m, n, nst, ist, jst, ast )

!*****************************************************************************80
!
!! R8ST_WRITE writes an R8ST file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the name of the file.
!
!    Input, integer M, the number of rows.
!
!    Input, integer N, the number of columns.
!
!    Input, integer NST, the number of nonzeros.
!
!    Input, integer IST(NST), JST(NST), the row and 
!    column indices.
!
!    Input, real ( kind = rk ) AST(NST), the nonzero values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nst

  real ( kind = rk ) ast(nst)
  integer ist(nst)
  integer ios
  integer jst(nst)
  integer k
  integer m
  integer n
  character ( len = * ) output_filename
  integer output_unit

  call i4_fake_use ( m )
  call i4_fake_use ( n )

  call get_unit ( output_unit )

  open ( file = output_filename, unit = output_unit, status = 'replace', &
    iostat = ios )

  do k = 1, nst

    write ( output_unit, '(2x,i8,2x,i8,2x,g16.8)', iostat = ios ) &
      ist(k), jst(k), ast(k)

  end do

  close ( unit = output_unit )

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

  integer, parameter :: rk = kind ( 1.0D+00 )

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

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
