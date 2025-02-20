module hb_file_module

!*****************************************************************************80
!
!! hb_file_module() stores data from an HB file.
!
!  Discussion:
!
!    In particular, the routine HB_FILE_READ does not return information
!    through its argument list, but instead, stores a copy of all the
!    HB data in this module.  In this way, a data structure of any
!    size can be read, and shared with the user calling routine.
!
!    A routine that wants to access this information must include the
!    statement "use hb_file_module".  To load data into the module,
!    call HB_FILE_READ.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    Technical Report TR/PA/92/86, CERFACS,
!    October 1992.
!
!  Parameters:
!
!    Local, character ( len = 72 ) TITLE, a title for the matrix.
!
!    Local, character ( len = 8 ) KEY, an identifier for the matrix.
!
!    Local, integer TOTCRD, the total number of lines of data.
!
!    Local, integer PTRCRD, the number of input lines for pointers.
!
!    Local, integer INDCRD, the number of input lines for
!    row indices.
!
!    Local, integer VALCRD, the number of input lines for
!    numerical values.
!
!    Local, integer RHSCRD, the number of input lines for
!    right hand sides.
!
!    Local, character ( len = 3 ) MXTYPE, the matrix type.
!    First character is R for Real, C for complex, P for pattern only.
!    Second character is S for symmetric, U for unsymmetric, H for
!      Hermitian, Z for skew symmetric, R for rectangular.
!    Third character is A for assembled and E for unassembled
!      finite element matrices.
!
!    Local, integer NROW, the number of rows or variables.
!
!    Local, integer NCOL, the number of columns or elements.
!
!    Local, integer NNZERO.  In the case of assembled sparse
!    matrices, this is the number of nonzeroes.  In the case of unassembled
!    finite element matrices, in which the right hand side vectors are also
!    stored as unassembled finite element vectors, this is the total
!    number of entries in a single unassembled right hand side vector.
!
!    Local, integer NELTVL, the number of finite element matrix
!    entries, set to 0 in the case of assembled matrices.
!
!    Local, character ( len = 16 ) PTRFMT, the format for reading pointers.
!
!    Local, character ( len = 16 ) INDFMT, the format for reading indices.
!
!    Local, character ( len = 20 ) VALFMT, the format for reading values.
!
!    Local, character ( len = 20 ) RHSFMT, the format for reading values
!    of the right hand side.
!
!    Local, character ( len = 3 ) RHSTYP, the right hand side type.
!    First character is F for full storage or M for same as matrix.
!    Second character is G if starting "guess" vectors are supplied.
!    Third character is X if exact solution vectors are supplied.
!    Ignored if NRHS = 0.
!
!    Local, integer NRHS, the number of right hand sides.
!
!    Local, integer NRHSIX, the number of row indices (set to 0
!    in the case of unassembled matrices.)  Ignored if NRHS = 0.
!
!    Local, integer COLPTR(NCOL+1), COLPTR(I) points to the
!    location of the first entry of column I in the sparse matrix structure.
!
!    If MXTYPE(3:3) == 'A':
!
!      Local, integer ROWIND(NNZERO), the row index of each item.
!
!      Local, real ( kind = rk ) VALUES(NNZERO), the nonzero values of
!      the matrix.
!
!    If MXTYPE(3:3) == 'E':
!
!      Local, integer ROWIND(NELTVL), the row index of each item.
!
!      Local, real ( kind = rk ) VALUES(NELTVL), the nonzero values of
!      the matrix.
!
!    If RHSTYP(1:1) == 'F':
!
!      Local, real ( kind = rk ) RHSVAL(NROW,NRHS), contains NRHS dense
!      right hand side vectors.
!
!      Local, integer RHSPTR(*), is not used.
!
!      Local, integer RHSIND(*), is not used.
!
!      Local, real ( kind = rk ) RHSVEC(*), is not used.
!
!    If RHSTYP(1:1) = 'M' and MXTYPE(3:3) = 'A':
!
!      Local, real ( kind = rk ) RHSVAL(*,*), is not used.
!
!      Local, integer RHSPTR(NRHS+1), RHSPTR(I) points to the
!      location of the first entry of right hand side I in the sparse right hand
!      side vector.
!
!      Local, integer RHSIND(NRHSIX), indicates, for each entry of
!      RHSVEC, the corresponding row index.
!
!      Local, real ( kind = rk ) RHSVEC(NRHSIX), contains the value of the
!      right hand side entries.
!
!    If RHSTYP(1:1) = 'M' and MXTYPE(3:3) = 'E':
!
!      Local, real ( kind = rk ) RHSVAL(NNZERO,NRHS), contains NRHS unassembled
!      finite element vector right hand sides.
!
!      Local, integer RHSPTR(*), is not used.
!
!      Local, integer RHSIND(*), is not used.
!
!      Local, real ( kind = rk ) RHSVEC(*), is not used.
!
!    Local, real ( kind = rk ) GUESS(NROW,NRHS), the starting guess vectors.
!
!    Local, real ( kind = rk ) EXACT(NROW,NRHS), the exact solution vectors.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, save, allocatable, dimension ( : ) :: colptr
  real ( kind = rk ), save, allocatable, dimension ( :, : ) :: exact
  real ( kind = rk ), save, allocatable, dimension ( :, : ) :: guess
  integer, save :: indcrd
  character ( len = 16 ), save :: indfmt
  character ( len =  8 ), save :: key
  character ( len =  3 ), save :: mxtype
  integer, save :: ncol
  integer, save :: neltvl
  integer, save :: nnzero
  integer, save :: nrhs
  integer, save :: nrhsix
  integer, save :: nrow
  integer, save :: ptrcrd
  character ( len = 16 ), save :: ptrfmt
  integer, save :: rhscrd
  character ( len = 20 ), save :: rhsfmt
  integer, save, allocatable, dimension ( : ) :: rhsind
  integer, save, allocatable, dimension ( : ) :: rhsptr
  character ( len =  3 ), save :: rhstyp
  real ( kind = rk ), save, allocatable, dimension ( :, : ) :: rhsval
  real ( kind = rk ), save, allocatable, dimension ( : ) :: rhsvec
  integer, save, allocatable, dimension ( : ) :: rowind
  character ( len = 72 ), save :: title
  integer, save :: totcrd
  integer, save :: valcrd
  character ( len = 20 ), save :: valfmt
  real ( kind = rk ), save, allocatable, dimension ( : ) :: values

end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! get_unit() returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
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
subroutine hb_exact_read ( input_unit, nrow, nrhs, rhscrd, rhsfmt, rhstyp, &
  exact )

!*****************************************************************************80
!
!! HB_EXACT_READ reads the exact solution vectors in an HB file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    Technical Report TR/PA/92/86, CERFACS,
!    October 1992.
!
!  Parameters:
!
!    Input, integer INPUT_UNIT, the unit from which data is read.
!
!    Input, integer NROW, the number of rows or variables.
!
!    Input, integer NRHS, the number of right hand sides.
!
!    Input, integer RHSCRD, the number of lines in the file for
!    right hand sides.
!
!    Input, character ( len = 20 ) RHSFMT, the format for reading values
!    of the right hand side.
!
!    Input, character ( len = 3 ) RHSTYP, the right hand side type.
!    First character is F for full storage or M for same as matrix.
!    Second character is G if starting "guess" vectors are supplied.
!    Third character is X if exact solution vectors are supplied.
!    Ignored if NRHS = 0.
!
!    Output, real ( kind = rk ) EXACT(NROW,NRHS), the exact solution vectors.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nrhs
  integer nrow

  real ( kind = rk ) exact(nrow,nrhs)
  integer input_unit
  integer ios
  integer rhscrd
  character ( len = 20 ) rhsfmt
  character ( len = 3 ) rhstyp

  if ( 0 < rhscrd ) then

    if ( rhstyp(3:3) == 'X' ) then

      read ( input_unit, rhsfmt, iostat = ios ) exact(1:nrow,1:nrhs)

      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'HB_EXACT_READ - Fatal error!'
        write ( *, '(a)' ) '  Error while reading EXACT.'
        stop 1
      end if

    end if

  end if

  return
end
subroutine hb_exact_write ( output_unit, nrow, nrhs, rhscrd, rhsfmt, rhstyp, &
  exact )

!*****************************************************************************80
!
!! HB_EXACT_WRITE writes the exact solution vectors to an HB file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    Technical Report TR/PA/92/86, CERFACS,
!    October 1992.
!
!  Parameters:
!
!    Input, integer OUTPUT_UNIT, the unit to which data
!    is written.
!
!    Input, integer NROW, the number of rows or variables.
!
!    Input, integer NRHS, the number of right hand sides.
!
!    Input, integer RHSCRD, the number of lines in the file for
!    right hand sides.
!
!    Input, character ( len = 20 ) RHSFMT, the format for reading values
!    of the right hand side.
!
!    Input, character ( len = 3 ) RHSTYP, the right hand side type.
!    First character is F for full storage or M for same as matrix.
!    Second character is G if starting "guess" vectors are supplied.
!    Third character is X if exact solution vectors are supplied.
!    Ignored if NRHS = 0.
!
!    Input, real ( kind = rk ) EXACT(NROW,NRHS), the exact solution vectors.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nrhs
  integer nrow

  real ( kind = rk ) exact(nrow,nrhs)
  integer output_unit
  integer rhscrd
  character ( len = 20 ) rhsfmt
  character ( len = 3 ) rhstyp

  if ( 0 < rhscrd ) then

    if ( rhstyp(3:3) == 'X' ) then
      write ( output_unit, rhsfmt ) exact(1:nrow,1:nrhs)
    end if

  end if

  return
end
subroutine hb_file_read ( input_unit )

!*****************************************************************************80
!
!! HB_FILE_READ reads an HB file.
!
!  Discussion:
!
!    This routine reads all the information from an HB file.
!
!    Since such a file may include arrays of arbitrary size, it
!    is convenient to invoke a module, which can be shared by this
!    routine, which allocates and stores the data, and the user.
!
!    A user routine which wishes to access the data stored in the
!    module must add the statement "use hb_file_module".
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    Technical Report TR/PA/92/86, CERFACS,
!    October 1992.
!
!  Parameters:
!
!    Input, integer INPUT_UNIT, the unit from data is read.
!
  use hb_file_module

  integer input_unit
!
!  Read the header block.
!
  call hb_header_read ( input_unit, title, key, totcrd, ptrcrd, indcrd, &
    valcrd, rhscrd, mxtype, nrow, ncol, nnzero, neltvl, ptrfmt, indfmt, &
    valfmt, rhsfmt, rhstyp, nrhs, nrhsix )
!
!  Read the matrix structure.
!
  if ( 0 < ptrcrd ) then

    if ( allocated ( colptr ) ) then
      deallocate ( colptr )
    end if
    allocate ( colptr(ncol+1) )

  end if

  if ( 0 < indcrd ) then

    if ( allocated ( rowind ) ) then
      deallocate ( rowind )
    end if

    if ( mxtype(3:3) == 'A' ) then
      allocate ( rowind(nnzero) )
    else if ( mxtype(3:3) == 'E' ) then
      allocate ( rowind(neltvl) )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HB_FILE_READ - Fatal error!'
      write ( *, '(a)' ) '  Illegal value of MXTYPE(3:3)!'
      stop 1
    end if

  end if

  call hb_structure_read ( input_unit, ncol, mxtype, nnzero, neltvl, &
    ptrcrd, ptrfmt, indcrd, indfmt, colptr, rowind )
!
!  Read the matrix values.
!
  if ( 0 < valcrd ) then

    if ( allocated ( values ) ) then
      deallocate ( values )
    end if

    if ( mxtype(3:3) == 'A' ) then
      allocate ( values(1:nnzero) )
    else if ( mxtype(3:3) == 'E' ) then
      allocate ( values(1:neltvl) )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HB_FILE_READ - Fatal error!'
      write ( *, '(a)' ) '  Illegal value of MXTYPE(3:3)!'
      stop 1
    end if

    call hb_values_read ( input_unit, valcrd, mxtype, nnzero, neltvl, &
      valfmt, values )

  end if
!
!  Read the right hand sides.
!
  if ( 0 < rhscrd ) then

    if ( rhstyp(1:1) == 'F' ) then

      if ( allocated ( rhsval ) ) then
        deallocate ( rhsval )
      end if

      allocate ( rhsval(nrow,nrhs) )

    else if ( rhstyp(1:1) == 'M' .and. mxtype(3:3) == 'A' ) then

      if ( allocated ( rhsptr ) ) then
        deallocate ( rhsptr )
      end if

      allocate ( rhsptr(nrhs+1) )

      if ( allocated ( rhsind ) ) then
        deallocate ( rhsind )
      end if

      allocate ( rhsind(nrhsix) )

      if ( allocated ( rhsvec ) ) then
        deallocate ( rhsvec )
      end if

      allocate ( rhsvec(nrhsix) )

    else if ( rhstyp(1:1) == 'M' .and. mxtype(3:3) == 'E' ) then

      if ( allocated ( rhsval ) ) then
        deallocate ( rhsval )
      end if

      allocate ( rhsval(nnzero,nrhs) )

    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HB_FILE_READ - Fatal error!'
      write ( *, '(a)' ) '  Illegal combination of RHSTYP(1:1) and MXTYPE(3:3)!'
      stop 1

    end if

    call hb_rhs_read ( input_unit, nrow, nnzero, nrhs, nrhsix, &
      rhscrd, ptrfmt, indfmt, rhsfmt, mxtype, rhstyp, rhsval, &
      rhsind, rhsptr, rhsvec )
!
!  Read the starting guesses.
!
    if ( rhstyp(2:2) == 'G' ) then

      if ( allocated ( guess ) ) then
        deallocate ( guess )
      end if

      allocate ( guess(1:nrow,1:nrhs) )

      call hb_guess_read ( input_unit, nrow, nrhs, rhscrd, rhsfmt, rhstyp, &
        guess )

    end if
!
!  Read the exact solutions.
!
    if ( rhstyp(3:3) == 'X' ) then

      if ( allocated ( exact ) ) then
        deallocate ( exact )
      end if

      allocate ( exact(1:nrow,1:nrhs) )

      call hb_exact_read ( input_unit, nrow, nrhs, rhscrd, rhsfmt, rhstyp, &
        exact )

    end if

  end if

  return
end
subroutine hb_file_write ( output_unit, title, key, totcrd, ptrcrd, indcrd, &
  valcrd, rhscrd, mxtype, nrow, ncol, nnzero, neltvl, ptrfmt, indfmt, &
  valfmt, rhsfmt, rhstyp, nrhs, nrhsix, colptr, rowind, values, &
  rhsval, rhsptr, rhsind, rhsvec, guess, exact )

!*****************************************************************************80
!
!! HB_FILE_WRITE writes an HB file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    Technical Report TR/PA/92/86, CERFACS,
!    October 1992.
!
!  Parameters:
!
!    Input, integer OUTPUT_UNIT, the unit to which data
!    is written.
!
!    Input, character ( len = 72 ) TITLE, a title for the matrix.
!
!    Input, character ( len = 8 ) KEY, an identifier for the matrix.
!
!    Input, integer TOTCRD, the total number of lines of data.
!
!    Input, integer PTRCRD, the number of input lines for pointers.
!
!    Input, integer INDCRD, the number of input lines for
!    row indices.
!
!    Input, integer VALCRD, the number of input lines for
!    numerical values.
!
!    Input, integer RHSCRD, the number of input lines for right
!    hand sides.
!
!    Input, character ( len = 3 ) MXTYPE, the matrix type.
!    First character is R for Real, C for complex, P for pattern only.
!    Second character is S for symmetric, U for unsymmetric, H for
!      Hermitian, Z for skew symmetric, R for rectangular.
!    Third character is A for assembled and E for unassembled
!      finite element matrices.
!
!    Input, integer NROW, the number of rows or variables.
!
!    Input, integer NCOL, the number of columns or elements.
!
!    Input, integer NNZERO.  In the case of assembled sparse
!    matrices, this is the number of nonzeroes.  In the case of unassembled
!    finite element matrices, in which the right hand side vectors are also
!    stored as unassembled finite element vectors, this is the total
!    number of entries in a single unassembled right hand side vector.
!
!    Input, integer NELTVL, the number of finite element matrix
!    entries, set to 0 in the case of assembled matrices.
!
!    Input, character ( len = 16 ) PTRFMT, the format for reading pointers.
!
!    Input, character ( len = 16 ) INDFMT, the format for reading indices.
!
!    Input, character ( len = 20 ) VALFMT, the format for reading values.
!
!    Input, character ( len = 20 ) RHSFMT, the format for reading values
!    of the right hand side.
!
!    Input, character ( len = 3 ) RHSTYP, the right hand side type.
!    First character is F for full storage or M for same as matrix.
!    Second character is G if starting "guess" vectors are supplied.
!    Third character is X if exact solution vectors are supplied.
!    Ignored if NRHS = 0.
!
!    Input, integer NRHS, the number of right hand sides.
!
!    Input, integer NRHSIX, the number of row indices (set to 0
!    in the case of unassembled matrices.)  Ignored if NRHS = 0.
!
!    Input, integer COLPTR(NCOL+1), COLPTR(I) points to the
!    location of the first entry of column I in the sparse matrix structure.
!
!    If MXTYPE(3:3) == 'A':
!
!      Input, integer ROWIND(NNZERO), the row index of each item.
!
!      Input, real ( kind = rk ) VALUES(NNZERO), the nonzero values of
!      the matrix.
!
!    If MXTYPE(3:3) == 'E':
!
!      Input, integer ROWIND(NELTVL), the row index of each item.
!
!      Input, real ( kind = rk ) VALUES(NELTVL), the nonzero values of
!      the matrix.
!
!    If RHSTYP(1:1) == 'F':
!
!      Input, real ( kind = rk ) RHSVAL(NROW,NRHS), contains NRHS dense
!      right hand side vectors.
!
!      Input, integer RHSPTR(*), is not used.
!
!      Input, integer RHSIND(*), is not used.
!
!      Input, real ( kind = rk ) RHSVEC(*), is not used.
!
!    If RHSTYP(1:1) = 'M' and MXTYPE(3:3) = 'A':
!
!      Input, real ( kind = rk ) RHSVAL(*,*), is not used.
!
!      Input, integer RHSPTR(NRHS+1), RHSPTR(I) points to the
!      location of the first entry of right hand side I in the sparse right
!      hand side vector.
!
!      Input, integer RHSIND(NRHSIX), indicates, for each entry of
!      RHSVEC, the corresponding row index.
!
!      Input, real ( kind = rk ) RHSVEC(NRHSIX), contains the value of the
!      right hand side entries.
!
!    If RHSTYP(1:1) = 'M' and MXTYPE(3:3) = 'E':
!
!      Input, real ( kind = rk ) RHSVAL(NNZERO,NRHS), contains NRHS unassembled
!      finite element vector right hand sides.
!
!      Input, integer RHSPTR(*), is not used.
!
!      Input, integer RHSIND(*), is not used.
!
!      Input, real ( kind = rk ) RHSVEC(*), is not used.
!
!    Input, real ( kind = rk ) GUESS(NROW,NRHS), the starting guess vectors.
!
!    Input, real ( kind = rk ) EXACT(NROW,NRHS), the exact solution vectors.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ncol
  integer neltvl
  integer nnzero

  integer colptr(ncol+1)
  real ( kind = rk ) exact(*)
  real ( kind = rk ) guess(*)
  integer indcrd
  character ( len = 16 ) indfmt
  character ( len = 8 ) key
  character ( len = 3 ) mxtype
  integer nrhs
  integer nrhsix
  integer nrow
  integer output_unit
  integer ptrcrd
  character ( len = 16 ) ptrfmt
  integer rhscrd
  character ( len = 20 ) rhsfmt
  integer rhsind(*)
  integer rhsptr(*)
  character ( len = 3 ) rhstyp
  real ( kind = rk ) rhsval(*)
  real ( kind = rk ) rhsvec(*)
  integer rowind(*)
  character ( len = 72 ) title
  integer totcrd
  integer valcrd
  character ( len = 20 ) valfmt
  real ( kind = rk ) values(*)
!
!  Write the header block.
!
  call hb_header_write ( output_unit, title, key, totcrd, ptrcrd, indcrd, &
    valcrd, rhscrd, mxtype, nrow, ncol, nnzero, neltvl, ptrfmt, indfmt, &
    valfmt, rhsfmt, rhstyp, nrhs, nrhsix )
!
!  Write the matrix structure.
!
  call hb_structure_write ( output_unit, ncol, mxtype, nnzero, neltvl, &
    ptrfmt, indfmt, colptr, rowind )
!
!  Write the matrix values.
!
  call hb_values_write ( output_unit, valcrd, mxtype, nnzero, neltvl, &
    valfmt, values )
!
!  Write the right hand sides.
!
  call hb_rhs_write ( output_unit, nrow, nnzero, nrhs, nrhsix, &
    rhscrd, ptrfmt, indfmt, rhsfmt, mxtype, rhstyp, rhsval, &
    rhsind, rhsptr, rhsvec )
!
!  Write the starting guesses.
!
  call hb_guess_write ( output_unit, nrow, nrhs, rhscrd, rhsfmt, rhstyp, &
    guess )
!
!  Write the exact solutions.
!
  call hb_exact_write ( output_unit, nrow, nrhs, rhscrd, rhsfmt, rhstyp, &
    exact )

  return
end
subroutine hb_guess_read ( input_unit, nrow, nrhs, rhscrd, rhsfmt, rhstyp, &
  guess )

!*****************************************************************************80
!
!! HB_GUESS_READ reads the starting guess vectors in an HB file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    Technical Report TR/PA/92/86, CERFACS,
!    October 1992.
!
!  Parameters:
!
!    Input, integer INPUT_UNIT, the unit from which data
!    is read.
!
!    Input, integer NROW, the number of rows or variables.
!
!    Input, integer NRHS, the number of right hand sides.
!
!    Input, integer RHSCRD, the number of lines in the file for
!    right hand sides.
!
!    Input, character ( len = 20 ) RHSFMT, the format for reading values
!    of the right hand side.
!
!    Input, character ( len = 3 ) RHSTYP, the right hand side type.
!    First character is F for full storage or M for same as matrix.
!    Second character is G if starting "guess" vectors are supplied.
!    Third character is X if exact solution vectors are supplied.
!    Ignored if NRHS = 0.
!
!    Output, real ( kind = rk ) GUESS(NROW,NRHS), the starting guess vectors.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nrhs
  integer nrow

  real ( kind = rk ) guess(nrow,nrhs)
  integer input_unit
  integer ios
  integer rhscrd
  character ( len = 20 ) rhsfmt
  character ( len = 3 ) rhstyp

  if ( 0 < rhscrd ) then

    if ( rhstyp(2:2) == 'G' ) then
      read ( input_unit, rhsfmt, iostat = ios ) guess(1:nrow,1:nrhs)

      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'HB_GUESS_READ - Fatal error!'
        write ( *, '(a)' ) '  Error while reading GUESS.'
        stop 1
      end if

    end if

  end if

  return
end
subroutine hb_guess_write ( output_unit, nrow, nrhs, rhscrd, rhsfmt, rhstyp, &
  guess )

!*****************************************************************************80
!
!! HB_GUESS_WRITE writes the starting guess vectors to an HB file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    Technical Report TR/PA/92/86, CERFACS,
!    October 1992.
!
!  Parameters:
!
!    Input, integer OUTPUT_UNIT, the unit to which data is written.
!
!    Input, integer NROW, the number of rows or variables.
!
!    Input, integer NRHS, the number of right hand sides.
!
!    Input, integer RHSCRD, the number of lines in the file for
!    right hand sides.
!
!    Input, character ( len = 20 ) RHSFMT, the format for reading values
!    of the right hand side.
!
!    Input, character ( len = 3 ) RHSTYP, the right hand side type.
!    First character is F for full storage or M for same as matrix.
!    Second character is G if starting "guess" vectors are supplied.
!    Third character is X if exact solution vectors are supplied.
!    Ignored if NRHS = 0.
!
!    Input, real ( kind = rk ) GUESS(NROW,NRHS), the starting guess vectors.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nrhs
  integer nrow

  real ( kind = rk ) guess(nrow,nrhs)
  integer output_unit
  integer rhscrd
  character ( len = 20 ) rhsfmt
  character ( len = 3 ) rhstyp

  if ( 0 < rhscrd ) then

    if ( rhstyp(2:2) == 'G' ) then
      write ( output_unit, rhsfmt ) guess(1:nrow,1:nrhs)
    end if

  end if

  return
end
subroutine hb_header_print ( title, key, totcrd, ptrcrd, indcrd, valcrd, &
  rhscrd, mxtype, nrow, ncol, nnzero, neltvl, ptrfmt, indfmt, valfmt, &
  rhsfmt, rhstyp, nrhs, nrhsix )

!*****************************************************************************80
!
!! HB_HEADER_PRINT prints the header of an HB file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    Technical Report TR/PA/92/86, CERFACS,
!    October 1992.
!
!  Parameters:
!
!    Input, character ( len = 72 ) TITLE, a title for the matrix.
!
!    Input, character ( len = 8 ) KEY, an identifier for the matrix.
!
!    Input, integer TOTCRD, the total number of lines of data.
!
!    Input, integer PTRCRD, the number of input lines for pointers.
!
!    Input, integer INDCRD, the number of input lines for
!    row indices.
!
!    Input, integer VALCRD, the number of input lines for
!    numerical values.
!
!    Input, integer RHSCRD, the number of input lines for
!    right hand sides.
!
!    Input, character ( len = 3 ) MXTYPE, the matrix type.
!    First character is R for Real, C for complex, P for pattern only.
!    Second character is S for symmetric, U for unsymmetric, H for
!      Hermitian, Z for skew symmetric, R for rectangular.
!    Third character is A for assembled and E for unassembled
!      finite element matrices.
!
!    Input, integer NROW, the number of rows or variables.
!
!    Input, integer NCOL, the number of columns or elements.
!
!    Input, integer NNZERO.  In the case of assembled sparse
!    matrices, this is the number of nonzeroes.  In the case of unassembled
!    finite element matrices, in which the right hand side vectors are also
!    stored as unassembled finite element vectors, this is the total
!    number of entries in a single unassembled right hand side vector.
!
!    Input, integer NELTVL, the number of finite element matrix
!    entries, set to 0 in the case of assembled matrices.
!
!    Input, character ( len = 16 ) PTRFMT, the format for reading pointers.
!
!    Input, character ( len = 16 ) INDFMT, the format for reading indices.
!
!    Input, character ( len = 20 ) VALFMT, the format for reading values.
!
!    Input, character ( len = 20 ) RHSFMT, the format for reading values
!    of the right hand side.
!
!    Input, character ( len = 3 ) RHSTYP, the right hand side type.
!    First character is F for full storage or M for same as matrix.
!    Second character is G if starting "guess" vectors are supplied.
!    Third character is X if exact solution vectors are supplied.
!
!    Input, integer NRHS, the number of right hand sides.
!
!    Input, integer NRHSIX, the number of entries of storage for
!    right hand side values, in the case where RHSTYP(1:1) = 'M' and
!    MXTYPE(3:3) = 'A'.
!
  implicit none

  integer indcrd
  character ( len = 16 ) indfmt
  character ( len = 8 ) key
  character ( len = 3 ) mxtype
  integer ncol
  integer neltvl
  integer nnzero
  integer nrhs
  integer nrhsix
  integer nrow
  integer ptrcrd
  character ( len = 16 ) ptrfmt
  integer rhscrd
  character ( len = 20 ) rhsfmt
  character ( len = 3 ) rhstyp
  character ( len = 72 ) title
  integer totcrd
  integer valcrd
  character ( len = 20 ) valfmt

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  TOTCRD = ', totcrd
  write ( *, '(a,i8)' ) '  PTRCRD = ', ptrcrd
  write ( *, '(a,i8)' ) '  INDCRD = ', indcrd
  write ( *, '(a,i8)' ) '  VALCRD = ', valcrd
  write ( *, '(a,i8)' ) '  RHSCRD = ', rhscrd
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  KEY =    "' // trim ( key )    // '".'
  write ( *, '(a)' ) '  MXTYPE = "' // trim ( mxtype ) // '".'
  write ( *, '(a)' ) '  RHSTYP = "' // trim ( rhstyp ) // '".'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  NROW =   ', nrow
  write ( *, '(a,i8)' ) '  NCOL =   ', ncol
  write ( *, '(a,i8)' ) '  NNZERO = ', nnzero
  write ( *, '(a,i8)' ) '  NELTVL = ', neltvl
  write ( *, '(a,i8)' ) '  NRHS =   ', nrhs
  write ( *, '(a,i8)' ) '  NRHSIX = ', nrhsix
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  PTRFMT = "' // trim ( ptrfmt ) // '".'
  write ( *, '(a)' ) '  INDFMT = "' // trim ( indfmt ) // '".'
  write ( *, '(a)' ) '  VALFMT = "' // trim ( valfmt ) // '".'
  write ( *, '(a)' ) '  RHSFMT = "' // trim ( rhsfmt ) // '".'

  return
end
subroutine hb_header_read ( input_unit, title, key, totcrd, ptrcrd, indcrd, &
  valcrd, rhscrd, mxtype, nrow, ncol, nnzero, neltvl, ptrfmt, indfmt, valfmt, &
  rhsfmt, rhstyp, nrhs, nrhsix )

!*****************************************************************************80
!
!! HB_HEADER_READ reads the header of an HB file.
!
!  Discussion:
!
!    The user should already have opened the file, and positioned it
!    to the first record.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    Technical Report TR/PA/92/86, CERFACS,
!    October 1992.
!
!  Parameters:
!
!    Input, integer INPUT_UNIT, the unit from which data is read.
!
!    Output, character ( len = 72 ) TITLE, a title for the matrix.
!
!    Output, character ( len = 8 ) KEY, an identifier for the matrix.
!
!    Output, integer TOTCRD, the total number of lines of data.
!
!    Output, integer PTRCRD, the number of input lines
!    for pointers.
!
!    Output, integer INDCRD, the number of input lines for
!    row indices.
!
!    Output, integer VALCRD, the number of input lines for
!    numerical values.
!
!    Output, integer RHSCRD, the number of input lines for right
!    hand sides.
!
!    Output, character ( len = 3 ) MXTYPE, the matrix type.
!    First character is R for Real, C for complex, P for pattern only.
!    Second character is S for symmetric, U for unsymmetric, H for
!      Hermitian, Z for skew symmetric, R for rectangular.
!    Third character is A for assembled and E for unassembled
!      finite element matrices.
!
!    Output, integer NROW, the number of rows or variables.
!
!    Output, integer NCOL, the number of columns or elements.
!
!    Output, integer NNZERO.  In the case of assembled sparse
!    matrices, this is the number of nonzeroes.  In the case of unassembled
!    finite element matrices, in which the right hand side vectors are also
!    stored as unassembled finite element vectors, this is the total
!    number of entries in a single unassembled right hand side vector.
!
!    Output, integer NELTVL, the number of finite element matrix
!    entries, set to 0 in the case of assembled matrices.
!
!    Output, character ( len = 16 ) PTRFMT, the format for reading pointers.
!
!    Output, character ( len = 16 ) INDFMT, the format for reading indices.
!
!    Output, character ( len = 20 ) VALFMT, the format for reading values.
!
!    Output, character ( len = 20 ) RHSFMT, the format for reading values
!    of the right hand side.
!
!    Output, character ( len = 3 ) RHSTYP, the right hand side type.
!    First character is F for full storage or M for same as matrix.
!    Second character is G if starting "guess" vectors are supplied.
!    Third character is X if exact solution vectors are supplied.
!
!    Output, integer NRHS, the number of right hand sides.
!
!    Output, integer NRHSIX, the number of entries of storage for
!    right hand side values, in the case where RHSTYP(1:1) = 'M' and
!    MXTYPE(3:3) = 'A'.
!
  implicit none

  integer indcrd
  character ( len = 16 ) indfmt
  integer input_unit
  integer ios
  character ( len = 8 ) key
  character ( len = 3 ) mxtype
  integer ncol
  integer neltvl
  integer nnzero
  integer nrhs
  integer nrhsix
  integer nrow
  integer ptrcrd
  character ( len = 16 ) ptrfmt
  integer rhscrd
  character ( len = 20 ) rhsfmt
  character ( len = 3 ) rhstyp
  character ( len = 72 ) title
  integer totcrd
  integer valcrd
  character ( len = 20 ) valfmt
!
!  Read the header block.
!
  read ( input_unit, '(a72,a8)', iostat = ios ) &
    title,  key

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HB_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  I/O error reading header line 1.'
    stop 1
  end if

  read ( input_unit, '(5i14)', iostat = ios ) &
    totcrd, ptrcrd, indcrd, valcrd, rhscrd

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HB_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  I/O error reading header line 2.'
    stop 1
  end if

  read ( input_unit, '(a3,11x,4i14)', iostat = ios ) &
    mxtype, nrow,   ncol,   nnzero, neltvl

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HB_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  I/O error reading header line 3.'
    stop 1
  end if

  read ( input_unit, '(2a16,2a20)', iostat = ios ) &
    ptrfmt, indfmt, valfmt, rhsfmt

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HB_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  I/O error reading header line 4.'
    stop 1
  end if

  if ( 0 < rhscrd ) then

    read ( input_unit, '(a3,11x,2i14)', iostat = ios ) &
      rhstyp, nrhs, nrhsix

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HB_HEADER_READ - Fatal error!'
      write ( *, '(a)' ) '  I/O error reading header line 5.'
      stop 1
    end if

  else

    rhstyp = ' '
    nrhs = 0
    nrhsix = 0

  end if

  return
end
subroutine hb_header_write ( output_unit, title, key, totcrd, ptrcrd, indcrd, &
  valcrd, rhscrd, mxtype, nrow, ncol, nnzero, neltvl, ptrfmt, indfmt, &
  valfmt, rhsfmt, rhstyp, nrhs, nrhsix )

!*****************************************************************************80
!
!! HB_HEADER_WRITE writes the header of an HB file.
!
!  Discussion:
!
!    The user should already have opened the file, and positioned it
!    to the first record.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    Technical Report TR/PA/92/86, CERFACS,
!    October 1992.
!
!  Parameters:
!
!    Input, integer OUTPUT_UNIT, the unit to which data is written.
!
!    Input, character ( len = 72 ) TITLE, a title for the matrix.
!
!    Input, character ( len = 8 ) KEY, an identifier for the matrix.
!
!    Input, integer TOTCRD, the total number of lines of data.
!
!    Input, integer PTRCRD, the number of input lines for pointers.
!
!    Input, integer INDCRD, the number of input lines for
!    row indices.
!
!    Input, integer VALCRD, the number of input lines for
!    numerical values.
!
!    Input, integer RHSCRD, the number of input lines for right
!    hand sides.
!
!    Input, character ( len = 3 ) MXTYPE, the matrix type.
!    First character is R for Real, C for complex, P for pattern only.
!    Second character is S for symmetric, U for unsymmetric, H for
!      Hermitian, Z for skew symmetric, R for rectangular.
!    Third character is A for assembled and E for unassembled
!      finite element matrices.
!
!    Input, integer NROW, the number of rows or variables.
!
!    Input, integer NCOL, the number of columns or elements.
!
!    Input, integer NNZERO.  In the case of assembled sparse
!    matrices, this is the number of nonzeroes.  In the case of unassembled
!    finite element matrices, in which the right hand side vectors are also
!    stored as unassembled finite element vectors, this is the total
!    number of entries in a single unassembled right hand side vector.
!
!    Input, integer NELTVL, the number of finite element matrix
!    entries, set to 0 in the case of assembled matrices.
!
!    Input, character ( len = 16 ) PTRFMT, the format for reading pointers.
!
!    Input, character ( len = 16 ) INDFMT, the format for reading indices.
!
!    Input, character ( len = 20 ) VALFMT, the format for reading values.
!
!    Input, character ( len = 20 ) RHSFMT, the format for reading values
!    of the right hand side.
!
!    Input, character ( len = 3 ) RHSTYP, the right hand side type.
!    First character is F for full storage or M for same as matrix.
!    Second character is G if starting "guess" vectors are supplied.
!    Third character is X if exact solution vectors are supplied.
!    Ignored if NRHS = 0.
!
!    Input, integer NRHS, the number of right hand sides.
!
!    Input, integer NRHSIX, the number of row indices (set to 0
!    in the case of unassembled matrices.)  Ignored if NRHS = 0.
!
  implicit none

  integer indcrd
  character ( len = 16 ) indfmt
  character ( len = 8 ) key
  character ( len = 3 ) mxtype
  integer ncol
  integer neltvl
  integer nnzero
  integer nrhs
  integer nrhsix
  integer nrow
  integer output_unit
  integer ptrcrd
  character ( len = 16 ) ptrfmt
  integer rhscrd
  character ( len = 20 ) rhsfmt
  character ( len = 3 ) rhstyp
  character ( len = 72 ) title
  integer totcrd
  integer valcrd
  character ( len = 20 ) valfmt

  write ( output_unit, '(a72,a8)' )       title,  key
  write ( output_unit, '(5i14)' )         totcrd, ptrcrd, indcrd, valcrd, rhscrd
  write ( output_unit, '(a3,11x,4i14)' )  mxtype, nrow,   ncol,   nnzero, neltvl
  write ( output_unit, '(2a16,2a20)' )    ptrfmt, indfmt, valfmt, rhsfmt

  if ( 0 < rhscrd ) then
    write ( output_unit, '(a3,11x,2i14)' ) rhstyp, nrhs, nrhsix
  end if

  return
end
subroutine hb_matvec_a_mem ( nrow, ncol, nnzero, nrhs, colptr, rowind, values, &
  exact, rhsval )

!*****************************************************************************80
!
!! HB_MATVEC_A_MEM multiplies an assembled Harwell Boeing matrix times a vector.
!
!  Discussion:
!
!    In this "A_MEM" version of the routine, the matrix is assumed to be in
!    "assembled" form, and all the data is assumed to be small enough
!    to reside completely in memory; the matrix and multiplicand vectors
!    are assumed to have been read into memory before this routine is called.
!
!    It is assumed that MXTYPE(3:3) = 'A', that is, that the matrix is
!    stored in the "assembled" format.
!
!    Also, the storage used for the vectors X and the products A*X
!    corresponds to RHSTYP(1:1) = 'F', that is, the "full" storage mode
!    for vectors.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    Technical Report TR/PA/92/86, CERFACS,
!    October 1992.
!
!  Parameters:
!
!    Input, integer NROW, the number of rows or variables.
!
!    Input, integer NCOL, the number of columns or elements.
!
!    Input, integer NNZERO.  In the case of assembled sparse
!    matrices, this is the number of nonzeroes.
!
!    Input, integer NRHS, the number of right hand sides.
!
!    Input, integer COLPTR(NCOL+1), COLPTR(I) points to the
!    location of the first entry of column I in the sparse matrix structure.
!
!    Input, integer ROWIND(NNZERO), the row index of each item.
!
!    Input, real ( kind = rk ) VALUES(NNZERO), the nonzero values of the matrix.
!
!    Input, real ( kind = rk ) EXACT(NCOL,NRHS), contains NRHS dense vectors.
!
!    Output, real ( kind = rk ) RHSVAL(NROW,NRHS), the product vectors A*X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ncol
  integer nnzero
  integer nrhs
  integer nrow

  integer colptr(ncol+1)
  integer column
  real ( kind = rk ) exact(ncol,nrhs)
  integer k
  real ( kind = rk ) rhsval(nrow,nrhs)
  integer row
  integer rowind(nnzero)
  real ( kind = rk ) values(nnzero)
!
!  Zero out the result vectors.
!
  rhsval(1:nrow,1:nrhs) = 0.0D+00
!
!  For each column J of the matrix,
!
  do column = 1, ncol
!
!  For nonzero entry K
!
    do k = colptr(column), colptr(column+1)-1

      row = rowind(k)
!
!  For each right hand side vector:
!
!    B(I,1:NRHS) = B(I,1:NRHS) + A(I,J) * X(J,1:NRHS)
!
      rhsval(row,1:nrhs) = rhsval(row,1:nrhs) + values(k) * exact(column,1:nrhs)

    end do

  end do

  return
end
subroutine hb_rhs_read ( input_unit, nrow, nnzero, nrhs, nrhsix, &
  rhscrd, ptrfmt, indfmt, rhsfmt, mxtype, rhstyp, rhsval, &
  rhsind, rhsptr, rhsvec )

!*****************************************************************************80
!
!! HB_RHS_READ reads the right hand side information in an HB file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    Technical Report TR/PA/92/86, CERFACS,
!    October 1992.
!
!  Parameters:
!
!    Input, integer INPUT_UNIT, the unit from which data is read.
!
!    Input, integer NROW, the number of rows or variables.
!
!    Input, integer NNZERO.  In the case of assembled sparse
!    matrices, this is the number of nonzeroes.  In the case of unassembled
!    finite element matrices, in which the right hand side vectors are also
!    stored as unassembled finite element vectors, this is the total
!    number of entries in a single unassembled right hand side vector.
!
!    Input, integer NRHS, the number of right hand sides.
!
!    Input, integer NRHSIX, the number of entries of storage for
!    right hand side values, in the case where RHSTYP(1:1) = 'M' and
!    MXTYPE(3:3) = 'A'.
!
!    Input, integer RHSCRD, the number of lines in the file for
!    right hand sides.
!
!    Input, character ( len = 16 ) PTRFMT, the format for reading pointers.
!
!    Input, character ( len = 16 ) INDFMT, the format for reading indices.
!
!    Input, character ( len = 20 ) RHSFMT, the format for reading values
!    of the right hand side.
!
!    Input, character ( len = 3 ) MXTYPE, the matrix type.
!    First character is R for Real, C for complex, P for pattern only.
!    Second character is S for symmetric, U for unsymmetric, H for
!      Hermitian, Z for skew symmetric, R for rectangular.
!    Third character is A for assembled and E for unassembled
!      finite element matrices.
!
!    Input, character ( len = 3 ) RHSTYP, the right hand side type.
!    First character is F for full storage or M for same as matrix.
!    Second character is G if starting "guess" vectors are supplied.
!    Third character is X if exact solution vectors are supplied.
!    Ignored if NRHS = 0.
!
!    If RHSTYP(1:1) == 'F':
!
!      Output, real ( kind = rk ) RHSVAL(NROW,NRHS), contains NRHS dense
!      right hand side vectors.
!
!      Output, integer RHSPTR(*), is not used.
!
!      Output, integer RHSIND(*), is not used.
!
!      Output, integer RHSVEC(*), is not used.
!
!    If RHSTYP(1:1) = 'M' and MXTYPE(3:3) = 'A':
!
!      Output, real ( kind = rk ) RHSVAL(*,*), is not used.
!
!      Output, integer RHSPTR(NRHS+1), RHSPTR(I) points to the
!      location of the first entry of right hand side I in the sparse right hand
!      side vector.
!
!      Output, integer RHSIND(NRHSIX), indicates, for each entry of
!      RHSVEC, the corresponding row index.
!
!      Output, real ( kind = rk ) RHSVEC(NRHSIX), contains the value of the
!      right hand side entries.
!
!    If RHSTYP(1:1) = 'M' and MXTYPE(3:3) = 'E':
!
!      Output, real ( kind = rk ) RHSVAL(NNZERO,NRHS), contains NRHS unassembled
!      finite element vector right hand sides.
!
!      Output, integer RHSPTR(*), is not used.
!
!      Output, integer RHSIND(*), is not used.
!
!      Output, real ( kind = rk ) RHSVEC(*), is not used.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nnzero
  integer nrhs
  integer nrhsix
  integer nrow

  character ( len = 16 ) indfmt
  integer input_unit
  integer ios
  character ( len = 3 ) mxtype
  character ( len = 16 ) ptrfmt
  integer rhscrd
  character ( len = 20 ) rhsfmt
  integer rhsind(*)
  integer rhsptr(nrhs+1)
  character ( len = 3 ) rhstyp
  real ( kind = rk ) rhsval(*)
  real ( kind = rk ) rhsvec(*)
!
!  Read the right hand sides.
!    case F                             = "full" or "dense";
!    case not F + matrix storage is "A" = sparse pointer RHS
!    case not F + matrix storage is "E" = finite element RHS
!
  if ( 0 < rhscrd ) then
!
!  Dense right hand sides:
!
    if ( rhstyp(1:1) == 'F' ) then

      read ( input_unit, rhsfmt, iostat = ios ) rhsval(1:nrow*nrhs)

      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'HB_RHS_READ - Fatal error!'
        write ( *, '(a)' ) '  I/O error while reading RHSVAL.'
        stop 1
      end if
!
!  Sparse right-hand sides stored like the matrix.
!  Read pointer array, indices, and values.
!
    else if ( rhstyp(1:1) == 'M' ) then

      if ( mxtype(3:3) == 'A' ) then

        read ( input_unit, ptrfmt, iostat = ios ) rhsptr(1:nrhs+1)

        if ( ios /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'HB_RHS_READ - Fatal error!'
          write ( *, '(a)' ) '  I/O error while reading RHSPTR.'
          stop 1
        end if

        read ( input_unit, indfmt, iostat = ios ) rhsind(1:nrhsix)

        if ( ios /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'HB_RHS_READ - Fatal error!'
          write ( *, '(a)' ) '  I/O error while reading RHSIND.'
          stop 1
        end if

        read ( input_unit, rhsfmt, iostat = ios ) rhsvec(1:nrhsix)

        if ( ios /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'HB_RHS_READ - Fatal error!'
          write ( *, '(a)' ) '  I/O error while reading RHSVEC.'
          stop 1
        end if
!
!  Sparse right hand sides in finite element format.
!
      else if ( mxtype(3:3) == 'E' ) then

        read ( input_unit, rhsfmt, iostat = ios ) rhsval(1:nnzero*nrhs)

        if ( ios /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'HB_RHS_READ - Fatal error!'
          write ( *, '(a)' ) '  I/O error while reading RHSVAL.'
          stop 1
        end if

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'HB_RHS_READ - Fatal error!'
        write ( *, '(a)' ) '  Illegal value of MXTYPE(3:3)!'
        stop 1

      end if
!
!  0 < RHSCRD, but RHSTYP not recognized,
!
    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HB_RHS_READ - Fatal error!'
      write ( *, '(a)' ) '  Illegal value of RHSTYP(1:1)!'
      stop 1

    end if

  end if

  return
end
subroutine hb_rhs_write ( output_unit, nrow, nnzero, nrhs, nrhsix, &
  rhscrd, ptrfmt, indfmt, rhsfmt, mxtype, rhstyp, rhsval, &
  rhsind, rhsptr, rhsvec )

!*****************************************************************************80
!
!! HB_RHS_WRITE writes the right hand side information to an HB file.
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
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    Technical Report TR/PA/92/86, CERFACS,
!    October 1992.
!
!  Parameters:
!
!    Input, integer OUTPUT_UNIT, the unit to which data is written.
!
!    Input, integer NROW, the number of rows or variables.
!
!    Input, integer NNZERO.  In the case of assembled sparse
!    matrices, this is the number of nonzeroes.  In the case of unassembled
!    finite element matrices, in which the right hand side vectors are also
!    stored as unassembled finite element vectors, this is the total
!    number of entries in a single unassembled right hand side vector.
!
!    Input, integer NRHS, the number of right hand sides.
!
!    Input, integer NRHSIX, the number of entries of storage for
!    right hand side values, in the case where RHSTYP(1:1) = 'M' and
!    MXTYPE(3:3) = 'A'.
!
!    Input, integer RHSCRD, the number of lines in the file for
!    right hand sides.
!
!    Input, character ( len = 16 ) PTRFMT, the format for reading pointers.
!
!    Input, character ( len = 16 ) INDFMT, the format for reading indices.
!
!    Input, character ( len = 20 ) RHSFMT, the format for reading values
!    of the right hand side.
!
!    Input, character ( len = 3 ) MXTYPE, the matrix type.
!    First character is R for Real, C for complex, P for pattern only.
!    Second character is S for symmetric, U for unsymmetric, H for
!      Hermitian, Z for skew symmetric, R for rectangular.
!    Third character is A for assembled and E for unassembled
!      finite element matrices.
!
!    Input, character ( len = 3 ) RHSTYP, the right hand side type.
!    First character is F for full storage or M for same as matrix.
!    Second character is G if starting "guess" vectors are supplied.
!    Third character is X if exact solution vectors are supplied.
!    Ignored if NRHS = 0.
!
!    If RHSTYP(1:1) == 'F':
!
!      Input, real ( kind = rk ) RHSVAL(NROW,NRHS), contains NRHS dense
!      right hand side vectors.
!
!      Input, integer RHSPTR(*), is not used.
!
!      Input, integer RHSIND(*), is not used.
!
!      Input, real ( kind = rk ) RHSVEC(*), is not used.
!
!    If RHSTYP(1:1) = 'M' and MXTYPE(3:3) = 'A':
!
!      Input, real ( kind = rk ) RHSVAL(*,*), is not used.
!
!      Input, integer RHSPTR(NRHS+1), RHSPTR(I) points to the
!      location of the first entry of right hand side I in the sparse right hand
!      side vector.
!
!      Input, integer RHSIND(NRHSIX), indicates, for each entry of
!      RHSVEC, the corresponding row index.
!
!      Input, real ( kind = rk ) RHSVEC(NRHSIX), contains the value of the
!      right hand side entries.
!
!    If RHSTYP(1:1) = 'M' and MXTYPE(3:3) = 'E':
!
!      Input, real ( kind = rk ) RHSVAL(NNZERO,NRHS), contains NRHS unassembled
!      finite element vector right hand sides.
!
!      Input, integer RHSPTR(*), is not used.
!
!      Input, integer RHSIND(*), is not used.
!
!      Input, real ( kind = rk ) RHSVEC(*), is not used.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nnzero
  integer nrhs
  integer nrhsix
  integer nrow

  character ( len = 16 ) indfmt
  character ( len = 3 ) mxtype
  integer output_unit
  character ( len = 16 ) ptrfmt
  integer rhscrd
  character ( len = 20 ) rhsfmt
  integer rhsind(*)
  integer rhsptr(nrhs+1)
  character ( len = 3 ) rhstyp
  real ( kind = rk ) rhsval(*)
  real ( kind = rk ) rhsvec(*)
!
!  Read the right hand sides.
!    case F                             = "full" or "dense";
!    case not F + matrix storage is "A" = sparse pointer RHS
!    case not F + matrix storage is "E" = finite element RHS
!
  if ( 0 < rhscrd ) then
!
!  Dense right hand sides:
!
    if ( rhstyp(1:1) == 'F' ) then

      write ( output_unit, rhsfmt ) rhsval(1:nrow*nrhs)
!
!  Sparse right-hand sides stored like the matrix.
!  Read pointer array, indices, and values.
!
    else if ( rhstyp(1:1) == 'M' ) then

      if ( mxtype(3:3) == 'A' ) then

        write ( output_unit, ptrfmt ) rhsptr(1:nrhs+1)

        write ( output_unit, indfmt ) rhsind(1:nrhsix)

        write ( output_unit, rhsfmt ) rhsvec(1:nrhsix)
!
!  Sparse right hand sides in finite element format.
!
      else if ( mxtype(3:3) == 'E' ) then

        write ( output_unit, rhsfmt ) rhsval(1:nnzero*nrhs)

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'HB_RHS_WRITE - Fatal error!'
        write ( *, '(a)' ) &
          '  Illegal value of MXTYPE(3:3) = "' // mxtype(3:3) // '".'
        stop 1

      end if

    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HB_RHS_WRITE - Fatal error!'
      write ( *, '(a)' ) &
        '  Illegal value of RHSTYP(1:1) = "' // rhstyp(1:1) // '".'
      stop 1

    end if

  end if

  return
end
subroutine hb_structure_print ( ncol, mxtype, nnzero, neltvl, colptr, rowind )

!*****************************************************************************80
!
!! HB_STRUCTURE_PRINT prints the structure of an HB matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    Technical Report TR/PA/92/86, CERFACS,
!    October 1992.
!
!  Parameters:
!
!    Input, integer NCOL, the number of columns.
!
!    Input, character ( len = 3 ) MXTYPE, the matrix type.
!    First character is R for Real, C for complex, P for pattern only.
!    Second character is S for symmetric, U for unsymmetric, H for
!      Hermitian, Z for skew symmetric, R for rectangular.
!    Third character is A for assembled and E for unassembled
!      finite element matrices.
!
!    Input, integer NNZERO.  In the case of assembled sparse
!    matrices, this is the number of nonzeroes.  In the case of unassembled
!    finite element matrices, in which the right hand side vectors are also
!    stored as unassembled finite element vectors, this is the total
!    number of entries in a single unassembled right hand side vector.
!
!    Input, integer NELTVL, the number of finite element matrix
!    entries, set to 0 in the case of assembled matrices.
!
!    Input, integer COLPTR(NCOL+1), COLPTR(I) points to the
!    location of the first entry of column I in the sparse matrix structure.
!
!    If MXTYPE(3:3) == 'A':
!
!      Input, integer ROWIND(NNZERO), the row index of each item.
!
!    If MXTYPE(3:3) == 'F':
!
!      Input, integer ROWIND(NELTVL), the row index of each item.
!
  implicit none

  integer ncol
  integer neltvl
  integer nnzero

  integer colptr(ncol+1)
  integer j
  integer khi
  integer klo
  character ( len = 3 ) mxtype
  integer rowind(*)

  call i4_fake_use ( neltvl )
  call i4_fake_use ( nnzero )

  if ( mxtype(3:3) == 'A' ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) &
      'Column Begin   End   ----------------------------------------'
    write ( *, '(a)' ) ' '
    do j = 1, min ( ncol, 10 )

      if ( colptr(j+1)-1 < colptr(j) ) then

        write ( *, '(i8,3x,a)' ) j, 'EMPTY'

      else

        do klo = colptr(j), colptr(j+1)-1, 10
          khi = min ( klo + 9, colptr(j+1)-1 )
          if ( klo == colptr(j) ) then
            write ( *, '(3i8,3x,10i4)' ) &
              j, colptr(j), colptr(j+1)-1, rowind(klo:khi)
          else
            write ( *, '(21x,10i4)' ) rowind(klo:khi)
          end if
        end do

      end if

    end do

    write ( *, '(a)' ) &
      '                     ----------------------------------------'

  else if ( mxtype(3:3) == 'E' ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) &
      'Column Begin   End   ----------------------------------------'
    write ( *, '(a)' ) &
      '                     ----------------------------------------'

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  I haven''t thought about how to print an'
    write ( *, '(a)' ) '  unassembled matrix yet!'

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HB_STRUCTURE_PRINT - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of MXTYPE(3:3) = ' // mxtype(3:3)
    stop 1

  end if

  return
end
subroutine hb_structure_read ( input_unit, ncol, mxtype, nnzero, neltvl, &
  ptrcrd, ptrfmt, indcrd, indfmt, colptr, rowind )

!*****************************************************************************80
!
!! HB_STRUCTURE_READ reads the structure of an HB matrix.
!
!  Discussion:
!
!    The user should already have opened the file, and positioned it
!    to just after the header records.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    Technical Report TR/PA/92/86, CERFACS,
!    October 1992.
!
!  Parameters:
!
!    Input, integer INPUT_UNIT, the unit from which data is read.
!
!    Input, integer NCOL, the number of columns.
!
!    Input, character ( len = 3 ) MXTYPE, the matrix type.
!    First character is R for Real, C for complex, P for pattern only.
!    Second character is S for symmetric, U for unsymmetric, H for
!      Hermitian, Z for skew symmetric, R for rectangular.
!    Third character is A for assembled and E for unassembled
!      finite element matrices.
!
!    Input, integer NNZERO.  In the case of assembled sparse
!    matrices, this is the number of nonzeroes.  In the case of unassembled
!    finite element matrices, in which the right hand side vectors are also
!    stored as unassembled finite element vectors, this is the total
!    number of entries in a single unassembled right hand side vector.
!
!    Input, integer NELTVL, the number of finite element matrix
!    entries, set to 0 in the case of assembled matrices.
!
!    Input, integer PTRCRD, the number of input lines for pointers.
!
!    Input, character ( len = 16 ) PTRFMT, the format for reading pointers.
!
!    Input, integer INDCRD, the number of input lines for indices.
!
!    Input, character ( len = 16 ) INDFMT, the format for reading indices.
!
!    Output, integer COLPTR(NCOL+1), COLPTR(I) points to the
!    location of the first entry of column I in the sparse matrix structure.
!
!    If MXTYPE(3:3) == 'A':
!
!      Output, integer ROWIND(NNZERO), the row index of each item.
!
!    If MXTYPE(3:3) == 'F':
!
!      Output, integer ROWIND(NELTVL), the row index of each item.
!
  implicit none

  integer ncol
  integer neltvl
  integer nnzero

  integer colptr(ncol+1)
  integer indcrd
  character ( len = 16 ) indfmt
  integer input_unit
  integer ios
  character ( len = 3 ) mxtype
  integer number
  integer ptrcrd
  character ( len = 16 ) ptrfmt
  integer rowind(*)

  call i4_fake_use ( neltvl )

  if ( 0 < ptrcrd ) then

    read ( input_unit, ptrfmt, iostat = ios ) colptr(1:ncol+1)

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HB_STRUCTURE_READ - Fatal error!'
      write ( *, '(a)' ) '  I/O error reading COLPTR.'
      stop 1
    end if

  end if

  if ( 0 < indcrd ) then

    if ( mxtype(3:3) == 'A' ) then

      read ( input_unit, indfmt, iostat = ios ) rowind(1:nnzero)

      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'HB_STRUCTURE_READ - Fatal error!'
        write ( *, '(a)' ) '  I/O error reading ROWIND.'
        stop 1
      end if

    else if ( mxtype(3:3) == 'E' ) then

      number = colptr(ncol) - colptr(1)

      read ( input_unit, indfmt, iostat = ios ) rowind(1:number)

      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'HB_STRUCTURE_READ - Fatal error!'
        write ( *, '(a)' ) '  I/O error reading ROWIND.'
        stop 1
      end if

    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HB_STRUCTURE_READ - Fatal error!'
      write ( *, '(a)' ) '  Illegal value of MXTYPE(3:3).'
      stop 1

    end if

  end if

  return
end
subroutine hb_structure_write ( output_unit, ncol, mxtype, nnzero, neltvl, &
  ptrfmt, indfmt, colptr, rowind )

!*****************************************************************************80
!
!! HB_STRUCTURE_WRITE writes the structure of an HB matrix.
!
!  Discussion:
!
!    If the user is creating an HB file, then the user should
!    already have opened the file, and written the header records.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    Technical Report TR/PA/92/86, CERFACS,
!    October 1992.
!
!  Parameters:
!
!    Input, integer OUTPUT_UNIT, the unit to which data is written.
!
!    Input, integer NCOL, the number of columns.
!
!    Input, character ( len = 3 ) MXTYPE, the matrix type.
!    First character is R for Real, C for complex, P for pattern only.
!    Second character is S for symmetric, U for unsymmetric, H for
!      Hermitian, Z for skew symmetric, R for rectangular.
!    Third character is A for assembled and E for unassembled
!      finite element matrices.
!
!    Input, integer NNZERO.  In the case of assembled sparse 
!    matrices, this is the number of nonzeroes.  In the case of unassembled 
!    finite element matrices, in which the right hand side vectors are also
!    stored as unassembled finite element vectors, this is the total
!    number of entries in a single unassembled right hand side vector.
!
!    Input, integer NELTVL, the number of finite element matrix
!    entries, set to 0 in the case of assembled matrices.
!
!    Input, character ( len = 16 ) PTRFMT, the format for reading pointers.
!
!    Input, character ( len = 16 ) INDFMT, the format for reading indices.
!
!    Input, integer COLPTR(NCOL+1), COLPTR(I) points to the location of
!    the first entry of column I in the sparse matrix structure.
!
!    If MXTYPE(3:3) == 'A':
!
!      Input, integer ROWIND(NNZERO), the row index of each item.
!
!    If MXTYPE(3:3) == 'F':
!
!      Input, integer ROWIND(NELTVL), the row index of each item.
!
  implicit none

  integer ncol
  integer neltvl
  integer nnzero

  integer colptr(ncol+1)
  character ( len = 16 ) indfmt
  character ( len = 3 ) mxtype
  integer number
  integer output_unit
  character ( len = 16 ) ptrfmt
  integer rowind(*)

  call i4_fake_use ( neltvl )

  write ( output_unit, ptrfmt ) colptr(1:ncol+1)

  if ( mxtype(3:3) == 'A' ) then

    write ( output_unit, indfmt ) rowind(1:nnzero)

  else if ( mxtype(3:3) == 'E' ) then

    number = colptr(ncol) - colptr(1)

    write ( output_unit, indfmt ) rowind(1:number)

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HB_STRUCTURE_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of MXTYPE(3:3).'
    stop 1

  end if

  return
end
subroutine hb_ua_colind ( ncol, colptr, nnzero, colind )

!*****************************************************************************80
!
!! HB_UA_COLUMN_INDEX: a column index for an unsymmetric assembled matrix.
!
!  Discussion:
!
!    It is assumed that the input data corresponds to a Harwell-Boeing
!    matrix which is unsymmetric, and assembled.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    Technical Report TR/PA/92/86, CERFACS,
!    October 1992.
!
!  Parameters:
!
!    Input, integer NCOL, the number of columns.
!
!    Input, integer COLPTR(NCOL+1), COLPTR(I) points to the
!    location of the first entry of column I in the sparse matrix structure.
!
!    Input, integer NNZERO, the number of nonzeros.
!
!    Output, integer COLIND(NNZERO), the column index of each
!    matrix value.
!
  implicit none

  integer ncol
  integer nnzero

  integer colind(nnzero)
  integer colptr(ncol+1)
  integer i

  do i = 1, ncol
    colind(colptr(i):colptr(i+1)-1) = i
  end do

  return
end
subroutine hb_values_print ( ncol, colptr, mxtype, nnzero, neltvl, values )

!*****************************************************************************80
!
!! HB_VALUES_PRINT prints the values of an HB matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    Technical Report TR/PA/92/86, CERFACS,
!    October 1992.
!
!  Parameters:
!
!    Input, integer NCOL, the number of columns.
!
!    Input, integer COLPTR(NCOL+1), COLPTR(I) points to the
!    location of the first entry of column I in the sparse matrix structure.
!
!    Input, character ( len = 3 ) MXTYPE, the matrix type.
!    First character is R for Real, C for complex, P for pattern only.
!    Second character is S for symmetric, U for unsymmetric, H for
!      Hermitian, Z for skew symmetric, R for rectangular.
!    Third character is A for assembled and E for unassembled
!      finite element matrices.
!
!    Input, integer NNZERO.  In the case of assembled sparse
!    matrices, this is the number of nonzeroes.  In the case of unassembled
!    finite element matrices, in which the right hand side vectors are also
!    stored as unassembled finite element vectors, this is the total
!    number of entries in a single unassembled right hand side vector.
!
!    Input, integer NELTVL, the number of finite element matrix
!    entries, set to 0 in the case of assembled matrices.
!
!    If MXTYPE(3:3) == 'A':
!
!      Input, real ( kind = rk ) VALUES(NNZERO), the nonzero values of
!      the matrix.
!
!    If MXTYPE(3:3) == 'E':
!
!      Input, real ( kind = rk ) VALUES(NELTVL), the nonzero values of
!      the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ncol
  integer neltvl
  integer nnzero

  integer colptr(ncol+1)
  integer j
  integer khi
  integer klo
  character ( len = 3 ) mxtype
  real ( kind = rk ) values(*)

  call i4_fake_use ( neltvl )
  call i4_fake_use ( nnzero )

  if ( mxtype(3:3) == 'A' ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) &
      'Column Begin   End   ----------------------------------------'
    write ( *, '(a)' ) ' '

    do j = 1, ncol

      if ( 5 < j .and. j < ncol ) then
        cycle
      end if

      if ( j == ncol .and. 6 < ncol ) then
        write ( *, '(a)' ) '(Skipping intermediate columns...)'
      end if

      if ( colptr(j+1)-1 < colptr(j) ) then

        write ( *, '(i8)' ) j, '   EMPTY'

      else

        do klo = colptr(j), colptr(j+1)-1, 5
          khi = min ( klo + 4, colptr(j+1)-1 )
          if ( klo == colptr(j) ) then
            write ( *, '(3i5,3x,5g12.4)' ) &
              j, colptr(j), colptr(j+1)-1, values(klo:khi)
          else
            write ( *, '(18x,5g12.4)' ) values(klo:khi)
          end if
        end do

      end if

    end do

    write ( *, '(a)' ) &
      '                     ----------------------------------------'
  else if ( mxtype(3:3) == 'E' ) then


    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) &
      'Column Begin   End   ----------------------------------------'
    write ( *, '(a)' ) &
      '                     ----------------------------------------'

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I haven''t thought about how to print an'
    write ( *, '(a)' ) 'unassembled matrix yet!'

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HB_VALUES_PRINT - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of MXTYPE(3:3) = ' // mxtype(3:3)
    stop 1

  end if

  return
end
subroutine hb_values_read ( input_unit, valcrd, mxtype, nnzero, neltvl, &
  valfmt, values )

!*****************************************************************************80
!
!! HB_VALUES_READ reads the values of an HB matrix.
!
!  Discussion:
!
!    The user should already have opened the file, and positioned it
!    to just after the header and structure records.
!
!    Values are contained in an HB file if the VALCRD parameter
!    is nonzero.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    Technical Report TR/PA/92/86, CERFACS,
!    October 1992.
!
!  Parameters:
!
!    Input, integer INPUT_UNIT, the unit from which data is read.
!
!    Input, integer VALCRD, the number of input lines for
!    numerical values.
!
!    Input, character ( len = 3 ) MXTYPE, the matrix type.
!    First character is R for Real, C for complex, P for pattern only.
!    Second character is S for symmetric, U for unsymmetric, H for
!      Hermitian, Z for skew symmetric, R for rectangular.
!    Third character is A for assembled and E for unassembled
!      finite element matrices.
!
!    Input, integer NNZERO.  In the case of assembled sparse
!    matrices, this is the number of nonzeroes.  In the case of unassembled
!    finite element matrices, in which the right hand side vectors are also
!    stored as unassembled finite element vectors, this is the total
!    number of entries in a single unassembled right hand side vector.
!
!    Input, integer NELTVL, the number of finite element matrix
!    entries, set to 0 in the case of assembled matrices.
!
!    Input, character ( len = 20 ) VALFMT, the format for reading values.
!
!    If MXTYPE(3:3) == 'A':
!
!      Output, real ( kind = rk ) VALUES(NNZERO), the nonzero values of
!      the matrix.
!
!    If MXTYPE(3:3) == 'E':
!
!      Output, real ( kind = rk ) VALUES(NELTVL), the nonzero values
!      of the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer neltvl
  integer nnzero

  integer input_unit
  integer ios
  character ( len = 3 ) mxtype
  integer valcrd
  character ( len = 20 ) valfmt
  real ( kind = rk ) values(*)
!
!  Read the matrix values.
!    case "A" = assembled;
!    case "E" = unassembled finite element matrices.
!
  if ( 0 < valcrd ) then

    if ( mxtype(3:3) == 'A' ) then

      read ( input_unit, valfmt, iostat = ios ) values(1:nnzero)

      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'HB_VALUES_READ - Fatal error!'
        write ( *, '(a)' ) '  I/O error reading assembled matrix values.'
        stop 1
      end if

    else if ( mxtype(3:3) == 'E' ) then

      read ( input_unit, valfmt, iostat = ios ) values(1:neltvl)

      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'HB_VALUES_READ - Fatal error!'
        write ( *, '(a)' ) '  I/O error reading unassembled matrix values.'
        stop 1
      end if

    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HB_VALUES_READ - Fatal error!'
      write ( *, '(a)' ) '  Illegal value of MXTYPE(3:3).'
      stop 1

    end if

  end if

  return
end
subroutine hb_values_write ( output_unit, valcrd, mxtype, nnzero, neltvl, &
  valfmt, values )

!*****************************************************************************80
!
!! HB_VALUES_WRITE writes the values of an HB matrix.
!
!  Discussion:
!
!    If the user is creating an HB file, then the user should already
!    have opened the file, and written the header and structure records.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    Technical Report TR/PA/92/86, CERFACS,
!    October 1992.
!
!  Parameters:
!
!    Input, integer OUTPUT_UNIT, the unit to which data is written.
!
!    Input, integer VALCRD, the number of input lines for
!    numerical values.
!
!    Input, character ( len = 3 ) MXTYPE, the matrix type.
!    First character is R for Real, C for complex, P for pattern only.
!    Second character is S for symmetric, U for unsymmetric, H for
!      Hermitian, Z for skew symmetric, R for rectangular.
!    Third character is A for assembled and E for unassembled
!      finite element matrices.
!
!    Input, integer NNZERO.  In the case of assembled sparse
!    matrices, this is the number of nonzeroes.  In the case of unassembled
!    finite element matrices, in which the right hand side vectors are also
!    stored as unassembled finite element vectors, this is the total
!    number of entries in a single unassembled right hand side vector.
!
!    Input, integer NELTVL, the number of finite element matrix
!    entries, set to 0 in the case of assembled matrices.
!
!    Input, character ( len = 20 ) VALFMT, the format for reading values.
!
!    If MXTYPE(3:3) == 'A':
!
!      Input, real ( kind = rk ) VALUES(NNZERO), the nonzero values of
!      the matrix.
!
!    If MXTYPE(3:3) == 'E':
!
!      Input, real ( kind = rk ) VALUES(NELTVL), the nonzero values
!      of the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer neltvl
  integer nnzero

  character ( len = 3 ) mxtype
  integer output_unit
  integer valcrd
  character ( len = 20 ) valfmt
  real ( kind = rk ) values(*)

  if ( 0 < valcrd ) then

    if ( mxtype(3:3) == 'A' ) then

      write ( output_unit, valfmt ) values(1:nnzero)

    else if ( mxtype(3:3) == 'E' ) then

      write ( output_unit, valfmt ) values(1:neltvl)

    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HB_VALUES_WRITE - Fatal error!'
      write ( *, '(a)' ) '  Illegal value of MXTYPE(3:3).'
      stop 1

    end if

  end if

  return
end
subroutine hb_vecmat_a_mem ( nrow, ncol, nnzero, nrhs, colptr, rowind, values, &
  exact, rhsval )

!*****************************************************************************80
!
!! HB_VECMAT_A_MEM multiplies a vector times an assembled Harwell Boeing matrix.
!
!  Discussion:
!
!    In this "A_MEM" version of the routine, the matrix is assumed to be in
!    "assembled" form, and all the data is assumed to be small enough
!    to reside completely in memory; the matrix and multiplicand vectors
!    are assumed to have been read into memory before this routine is called.
!
!    It is assumed that MXTYPE(3:3) = 'A', that is, that the matrix is
!    stored in the "assembled" format.
!
!    Also, the storage used for the vectors X and the products A*X
!    corresponds to RHSTYP(1:1) = 'F', that is, the "full" storage mode
!    for vectors.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Iain Duff, Roger Grimes, John Lewis,
!    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
!    Technical Report TR/PA/92/86, CERFACS,
!    October 1992.
!
!  Parameters:
!
!    Input, integer NROW, the number of rows or variables.
!
!    Input, integer NCOL, the number of columns or elements.
!
!    Input, integer NNZERO.  In the case of assembled sparse
!    matrices, this is the number of nonzeroes.
!
!    Input, integer NRHS, the number of right hand sides.
!
!    Input, integer COLPTR(NCOL+1), COLPTR(I) points to the
!    location of the first entry of column I in the sparse matrix structure.
!
!    Input, integer ROWIND(NNZERO), the row index of each item.
!
!    Input, real ( kind = rk ) VALUES(NNZERO), the nonzero values of the matrix.
!
!    Input, real ( kind = rk ) EXACT(NROW,NRHS), contains NRHS dense vectors.
!
!    Output, real ( kind = rk ) RHSVAL(NCOL,NRHS), the product vectors A'*X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ncol
  integer nnzero
  integer nrhs
  integer nrow

  integer colptr(ncol+1)
  integer column
  real ( kind = rk ) exact(nrow,nrhs)
  integer k
  real ( kind = rk ) rhsval(ncol,nrhs)
  integer row
  integer rowind(nnzero)
  real ( kind = rk ) values(nnzero)
!
!  Zero out the result vectors.
!
  rhsval(1:ncol,1:nrhs) = 0.0D+00
!
!  For each column J of the matrix,
!
  do column = 1, ncol
!
!  For nonzero entry K
!
    do k = colptr(column), colptr(column+1)-1

      row = rowind(k)
!
!  For each right hand side vector:
!
!    B(J,1:NRHS) = B(J,1:NRHS) + X(I,1:NRHS) * A(I,J)
!
      rhsval(column,1:nrhs) = rhsval(column,1:nrhs) &
        + exact(row,1:nrhs) * values(k)

    end do

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
    write ( *, '(a)' ) '  i4_fake_use(): variable is NAN.'
  end if

  return
end
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 November 2000
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
    write ( *, '(2x,i8,2x,i12)' ) i, a(i)
  end do

  return
end
subroutine i4vec_print_some ( n, a, i_lo, i_hi, title )

!*****************************************************************************80
!
!! I4VEC_PRINT_SOME prints "some" of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4 values.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 October 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries of the vector.
!
!    Input, integer A(N), the vector to be printed.
!
!    Input, integer I_LO, I_HI, the first and last indices
!    to print.  The routine expects 1 <= I_LO <= I_HI <= N.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer n

  integer a(n)
  integer i
  integer i_hi
  integer i_lo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = max ( i_lo, 1 ), min ( i_hi, n )
    write ( *, '(2x,i8,2x,i8)' ) i, a(i)
  end do

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints a R8MAT.
!
!  Discussion:
!
!    A R8MAT is a two dimensional matrix of double precision real values.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 September 2004
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
!    Input, real ( kind = rk ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  character ( len = * ) title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_PRINT_SOME prints some of a R8MAT.
!
!  Discussion:
!
!    A R8MAT is a two dimensional matrix of double precision real values.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = rk ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ILO, JLO, the first row and column to print.
!
!    Input, integer IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: incx = 5
  integer m
  integer n

  real ( kind = rk ) a(m,n)
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

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == real ( int ( a(i,j) ), kind = rk ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints a R8VEC.
!
!  Discussion:
!
!    A R8VEC is an array of double precision real values.
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
    write ( *, '(2x,i8,2x,g16.8)' ) i, a(i)
  end do

  return
end
subroutine r8vec_print_some ( n, a, i_lo, i_hi, title )

!*****************************************************************************80
!
!! R8VEC_PRINT_SOME prints "some" of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8 values.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 October 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries of the vector.
!
!    Input, real ( kind = rk ) A(N), the vector to be printed.
!
!    Input, integer I_LO, I_HI, the first and last indices
!    to print.  The routine expects 1 <= I_LO <= I_HI <= N.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  integer i
  integer i_hi
  integer i_lo
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = max ( i_lo, 1 ), min ( i_hi, n )
    write ( *, '(2x,i8,2x,g14.8)' ) i, a(i)
  end do

  return
end

