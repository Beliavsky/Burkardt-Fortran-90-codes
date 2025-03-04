program main

!*****************************************************************************80
!
!! st_to_hb() converts an ST sparse matrix file to HB format.
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
!  Usage:
!
!    st_to_hb file.st file.hb
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable, dimension ( : ) :: a
  integer arg_num
  integer base0
  integer base1
  integer, allocatable, dimension ( : ) :: col
  integer, allocatable, dimension ( : ) :: colptr
  integer i4vec_dummy(1)
  integer iarg
  integer iargc
  integer indcrd
  character ( len = 20 ) indfmt
  character ( len = 255 ) input_filename
  integer j
  integer k
  character ( len = 8 ) key
  character ( len = 3 ) mxtype
  integer ncol
  integer neltvl
  integer nnzero
  integer nrhs
  integer nrhsix
  integer nrow
  character ( len = 255 ) output_filename
  integer output_unit
  integer ptrcrd
  character ( len = 20 ) ptrfmt
  real ( kind = rk ) r8mat_dummy(1,1)
  real ( kind = rk ) r8vec_dummy(1)
  integer rhscrd
  character ( len = 20 ) rhsfmt
  character ( len = 3 )  rhstyp
  integer, allocatable, dimension ( : ) :: row
  character ( len = 72 ) title
  integer totcrd
  integer valcrd
  character ( len = 20 ) valfmt

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'st_to_hb():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Read an ST sparse matrix file,'
  write ( *, '(a)' ) &
    '  write the corresponding Harwell-Boeing sparse matrix file.'
!
!  Get the number of command line arguments.
!
  arg_num = iargc ( )
!
!  If at least one command line argument, it's the input file name.
!
  if ( 1 <= arg_num ) then

    iarg = 1
    call getarg ( iarg, input_filename )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ST_TO_HB:'
    write ( *, '(a)' ) '  Please enter the name of the input file.'

    read ( *, '(a)' ) input_filename

  end if
!
!  Second command line argument is the output file name.
!
  if ( 2 <= arg_num ) then

    iarg = 2
    call getarg ( iarg, output_filename )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ST_TO_HB:'
    write ( *, '(a)' ) '  Please enter the name of the output file.'

    read ( *, '(a)' ) output_filename

  end if
!
!  Size ST information.
!
  call st_header_read ( input_filename, nrow, ncol, nnzero )
!
!  Read ST information.
!
  allocate ( a(1:nnzero) )
  allocate ( row(1:nnzero) )
  allocate ( col(1:nnzero) )

  call st_data_read ( input_filename, nrow, ncol, nnzero, row, col, a )
!
!  Adjust ST information.
!
  base0 = 0
  base1 = 1
  call st_rebase ( base0, base1, nnzero, row )
  call st_rebase ( base0, base1, nnzero, col )
!
!  Sort (ROW,COL,A) by columns.
!
  call st_sort_a ( nrow, ncol, nnzero, row, col, a )
!
!  Determine the COLPTR array.
!
  allocate ( colptr(1:ncol+1) )

  colptr(1) = 1
  j = 2
  do k = 1, nnzero
    if ( col(k) == j ) then
      colptr(j) = k
      j = j + 1
    end if
  end do
  colptr(ncol+1) = nnzero + 1
!
!  Write HB information.
!
  call get_unit ( output_unit )

  open ( file = output_filename, unit = output_unit, status = 'replace' )

  title(1:72) = input_filename(1:72)
  key = 'Key     '

  ptrcrd = ( ncol - 1 ) / 16 + 1
  indcrd = ( nnzero - 1 ) / 16 + 1
  valcrd = ( nnzero - 1 ) / 10 + 1
  rhscrd = 0

  totcrd = ptrcrd + indcrd + valcrd + rhscrd

  mxtype = 'RUA'
  neltvl = 0
  ptrfmt = '(16I5)'
  indfmt = '(16I5)'
  valfmt = '(10F7.1)'
  rhsfmt = '(10F7.1)'
  rhstyp = 'FGX'
  nrhs = 0
  nrhsix = 0

  call hb_file_write ( output_unit, title, key, totcrd, ptrcrd, indcrd, &
    valcrd, rhscrd, mxtype, nrow, ncol, nnzero, neltvl, ptrfmt, indfmt, &
    valfmt, rhsfmt, rhstyp, nrhs, nrhsix, colptr, row, a, &
    r8mat_dummy, i4vec_dummy, i4vec_dummy, r8vec_dummy, r8mat_dummy, &
    r8mat_dummy )

  close ( unit = output_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ST_TO_HB:'
  write ( *, '(a)' ) '  Read data from ST file: "' &
    // trim ( input_filename ) // '".'
  write ( *, '(a)' ) '  Created output HB file: "' &
    // trim ( output_filename ) // '".'
!
!  Free memory.
!
  deallocate ( a )
  deallocate ( col )
  deallocate ( colptr )
  deallocate ( row )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ST_TO_HB:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine sort_heap_external ( n, indx, i, j, isgn )

!*****************************************************************************80
!
!! sort_heap_external() externally sorts a list of items into ascending order.
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
!
!    The user must set INDX to 0 before the first call.
!    Thereafter, the user should not change the value of INDX until
!    the sorting is done.
!
!    On return, if INDX is
!
!      greater than 0,
!      * interchange items I and J;
!      * call again.
!
!      less than 0,
!      * compare items I and J;
!      * set ISGN = -1 if I < J, ISGN = +1 if J < I;
!      * call again.
!
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
subroutine st_data_read ( input_filename, nrow, ncol, nnzero, row, col, a )

!*****************************************************************************80
!
!! ST_DATA_READ reads the data of an ST file.
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
!    Input, character ( len = * ) INPUT_FILENAME, the name of the ST file.
!
!    Input, integer NROW, the assumed number of rows in 
!    the matrix.
!
!    Input, integer NCOL, the assumed number of columns 
!    in the matrix.
!
!    Input, integer NNZERO, the assumed number of nonzeros 
!    in the matrix.
!
!    Output, integer ROW(NNZERO), COL(NNZERO), the 0-based row 
!    and column index of a nonzero matrix entry.
!
!    Output, real ( kind = rk ) A(NNZERO), the value of a nonzero matrix entry.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nnzero

  real ( kind = rk ) a(nnzero)
  real ( kind = rk ) aij
  integer col(nnzero)
  integer i
  character ( len = * )  input_filename
  integer input_unit
  integer ios
  integer j
  integer k
  integer ncol
  integer nrow
  integer row(nnzero)

  call i4_fake_use ( nrow )
  call i4_fake_use ( ncol )

  call get_unit ( input_unit )

  open ( file = input_filename, unit = input_unit, status = 'old', &
    iostat = ios )

  do k = 1, nnzero

    read ( input_unit, *, iostat = ios ) i, j, aij

    if ( ios /= 0 ) then
      exit
    end if

    row(k) = i
    col(k) = j
    a(k) = aij

  end do

  close ( unit = input_unit )

  return
end
subroutine st_header_read ( input_filename, nrow, ncol, nnzero )

!*****************************************************************************80
!
!! ST_HEADER_READ reads the header of an ST file.
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
!    Input, character ( len = * ) INPUT_FILENAME, the name of the ST file.
!
!    Output, integer NROW, the assumed number of rows in 
!    the matrix.
!
!    Output, integer NCOL, the assumed number of columns 
!    in the matrix.
!
!    Output, integer NNZERO, the assumed number of nonzeros 
!    in the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) aij
  integer i
  character ( len = * ) input_filename
  integer input_unit
  integer ios
  integer j
  integer ncol
  integer nnzero
  integer nrow

  call get_unit ( input_unit )

  open ( file = input_filename, unit = input_unit, status = 'old', &
    iostat = ios )

  nnzero = 0
  nrow = 0
  ncol = 0

  do

    read ( input_unit, *, iostat = ios ) i, j, aij

    if ( ios /= 0 ) then
      exit
    end if

    nnzero = nnzero + 1
    nrow = max ( nrow, i + 1 )
    ncol = max ( ncol, j + 1 )

  end do

  close ( unit = input_unit )

  return
end
subroutine st_rebase ( base1, base2, nnzero, indx )

!*****************************************************************************80
!
!! ST_REBASE changes the base of an index array.
!
!  Discussion:
!
!    Both the ROW and COL arrays are ordinarily 0-based in the ST format.
!    FORTRAN and MATLAB expect 1-based indices.
!
!    To convert ROW and COL from 0-based to 1-based form, call this routine
!    with BASE1=0, BASE2=1.
!
!    If ROW and COL from FORTRAN or MATLAB are to be converted to ST format,
!    call this routine with BASE1=1 and BASE2=0.
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
!    Input, integer BASE1, BASE2, the old and new bases.
!
!    Input, integer NNZERO, the number of nonzeros in the matrix.
!
!    Input/output, integer INDX(NNZERO), the index vector
!    to be rebased.
!
  implicit none

  integer nnzero

  integer base1
  integer base2
  integer indx(nnzero)

  indx(1:nnzero) = indx(1:nnzero) - base1 + base2

  return
end
subroutine st_sort_a ( nrow, ncol, nnzero, row, col, a )

!*****************************************************************************80
!
!! ST_SORT_A sorts the entries of an ST matrix by column.
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
!    Input, integer NROW, the number of rows in the matrix.
!
!    Input, integer NCOL, the number of columns in the matrix.
!
!    Input, integer NNZERO, the number of nonzeros in the matrix.
!
!    Input/output, integer ROW(NNZERO), COL(NNZERO), the
!    0-based row and column index of a nonzero matrix entry.
!
!    Input/output, real ( kind = rk ) A(NNZERO), the nonzero matrix entries.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nnzero

  real ( kind = rk ) a(nnzero)
  real ( kind = rk ) aij
  integer cij
  integer col(nnzero)
  integer i
  integer indx
  integer isgn
  integer j
  integer ncol
  integer nrow
  integer rij
  integer row(nnzero)

  call i4_fake_use ( nrow )
  call i4_fake_use ( ncol )
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

    call sort_heap_external ( nnzero, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      rij    = row(i)
      row(i) = row(j)
      row(j) = rij

      cij    = col(i)
      col(i) = col(j)
      col(j) = cij

      aij  = a(i)
      a(i) = a(j)
      a(j) = aij
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      if ( col(i) == col(j) ) then

        if ( row(i) < row(j) ) then
          isgn = - 1
        else if ( row(i) == row(j) ) then
          isgn = 0
        else if ( row(j) < row(i) ) then
          isgn = + 1
        end if

      else if ( col(i) < col(j) ) then

        isgn = - 1

      else if ( col(j) < col(i) ) then

        isgn = + 1

      end if

    else if ( indx == 0 ) then

      exit

    end if

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
