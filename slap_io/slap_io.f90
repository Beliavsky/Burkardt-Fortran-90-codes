subroutine slap_file_print ( n, nelt, isym, irhs, isoln, ia, ja, a, rhs, &
  soln, title )

!*****************************************************************************80
!
!! slap_file_print() prints a slap linear system that was stored in a file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Mark Seager,
!    A SLAP for the Masses,
!    Lawrence Livermore National Laboratory,
!    Technical Report UCRL-100267, December 1988.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer NELT, the number of non-zeros stored in A.
!
!    Input, integer ISYM, a flag to indicate symmetric 
!    storage format.
!    * 0, all nonzero entries of the matrix are stored.
!    * 1, the matrix is symmetric, and only the lower triangle of the
!    matrix is stored.
!
!    Input, integer IRHS, is 1 if a right hand side vector 
!    is included.
!
!    Input, integer ISOLN, is 1 if a solution vector is included.
!
!    Input, integer IA(NELT), integer JA(NELT),
!    real ( kind = rk ) A(NELT), the slap triad matrix description.
!
!    Input, real ( kind = rk ) RHS(N), the right hand side vector.
!
!    Input, real ( kind = rk ) SOLN(N), the solution to the linear system.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nelt

  real ( kind = rk ) a(nelt)
  integer ia(nelt)
  integer irhs
  integer isoln
  integer isym
  integer ja(nelt)
  real ( kind = rk ) rhs(n)
  real ( kind = rk ) soln(n)
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  call slap_header_print ( n, nelt, isym, irhs, isoln )
!
!  Write out the matrix.
!
  call r8sp_print ( n, n, nelt, isym, ia, ja, a, '  The sparse matrix' )
!
!  Write the right hand side.
!
  if ( irhs == 1 ) then
    call slap_rhs_print ( n, rhs )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  No right hand side vector was supplied.'
  end if
!
!  Write the solution.
!
  if ( isoln == 1 ) then
    call slap_soln_print ( n, soln )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  No solution vector was supplied.'
  end if

  return
end
subroutine slap_file_read ( n_max, nelt_max, n, nelt, isym, irhs, isoln, &
  ia, ja, a, rhs, soln, iunit )

!*****************************************************************************80
!
!! slap_FILE_READ reads in a slap matrix contained in a file.
!
!  Discussion:
!
!    This routine reads in a slap Triad Format Linear System,
!    including the matrix, right hand side, and solution, if known.
!
!    The original version of this program seems to have a minor
!    logical flaw.  If the user requests the solution but not
!    the right hand side, and the file contains both, the original
!    program would not correctly read past the right hand side
!    to get to the solution.  The current version should fix
!    that flaw.
!
!    The expected format of the file is as follows.  On the first line
!    are counters and flags: N, NELT, ISYM, IRHS, ISOLN.  N, NELT
!    and ISYM are described below.  IRHS is a flag indicating if
!    the RHS was written out (1 is yes, 0 is  no).  ISOLN is a
!    flag indicating if the SOLN was written out  (1 is yes, 0 is
!    no).  The format for the first line is: 5i10.  Then comes the
!    NELT Triad's IA(I), JA(I) and A(I), I = 1, NELT.  The format
!    for these lines is   :  1X,I5,1X,I5,1X,E16.7.   Then comes
!    RHS(I), I = 1, N, if IRHS = 1.  Then comes SOLN(I), I  = 1,
!    N, if ISOLN = 1.  The format for these lines is: 1X,E16.7.
!
!
!    This routine requires that the  matrix A be stored in the
!    slap Triad format.  In this format only the non-zeros  are
!    stored.  They may appear in ANY order.  The user supplies
!    three arrays of length NELT, where NELT is the number of
!    non-zeros in the matrix: (IA(NELT), JA(NELT), A(NELT)).  For
!    each non-zero the user puts the row and column index of that
!    matrix element in the IA and JA arrays.  The value of the
!    non-zero matrix element is placed in the corresponding
!    location of the A array.   This is an extremely easy data
!    structure to generate.  On the other hand it is not too
!    efficient on vector computers for the iterative solution of
!    linear systems.  Hence,  slap changes this input data
!    structure to the slap Column format for the iteration (but
!    does not change it back).
!
!    Here is an example of the slap Triad storage format for a
!    5x5 Matrix.  Recall that the entries may appear in any order.
!
!        5x5 Matrix       slap Triad format for 5x5 matrix on left.
!                              1  2  3  4  5  6  7  8  9 10 11
!    |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21
!    |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2
!    | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1
!    | 0  0  0 44  0|
!    |51  0 53  0 55|
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Mark Seager,
!    A SLAP for the Masses,
!    Lawrence Livermore National Laboratory,
!    Technical Report UCRL-100267, December 1988.
!
!  Parameters:
!
!    Input, integer N_MAX, the maximum value of N for which storage
!    has been allocated.
!
!    Input, integer NELT_MAX, the maximum value of NELT for which
!    storage has been allocated.
!
!    Output, integer N, the order of the matrix.
!
!    Output, integer NELT, the number of non-zeros stored in A.
!
!    Output, integer ISYM, a flag to indicate symmetric storage 
!    format.
!    * 0, all nonzero entries of the matrix are stored.
!    * 1, the matrix is symmetric, and only the lower triangle of the
!    matrix is stored.
!
!    Output, integer IRHS, is 1 if a right hand side vector is 
!    included.
!
!    Output, integer ISOLN, is 1 if a solution vector is included.
!
!    Output, integer IA(NELT), integer JA(NELT),
!    real ( kind = rk ) A(NELT).  On output these arrays hold the matrix A in the
!    slap Triad format.
!
!    Output, real ( kind = rk ) RHS(N), the right hand side vector.
!
!    Output, real ( kind = rk ) SOLN(N), the solution to the linear system, 
!    if present.
!
!    Input, integer IUNIT, the FORTRAN device unit from which the
!    matrix is to be read.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n_max
  integer nelt_max

  real ( kind = rk ) a(nelt_max)
  integer i
  integer ia(nelt_max)
  integer ios
  integer irhs
  integer isoln
  integer isym
  integer iunit
  integer ja(nelt_max)
  integer n
  integer nelt
  real ( kind = rk ) rhs(n_max)
  real ( kind = rk ) soln(n_max)
!
!  Read the header line.
!
  call slap_header_read ( iunit, n, nelt, isym, irhs, isoln, ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'slap_FILE_READ - Fatal error!'
    write ( *, '(a)' ) '  Error while reading header line of slap file.'
    stop 1
  end if

  if ( nelt_max < nelt ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'slap_FILE_READ - Fatal error!'
    write ( *, '(a)' ) '  NELT_MAX < NELT.'
    write ( *, '(a,i8)' ) '  NELT_MAX = ', nelt_max
    write ( *, '(a,i8)' ) '  NELT     = ', nelt
    stop 1
  end if
!
!  Read the nonzero matrix entries in slap Triad format.
!
  do i = 1, nelt

    read ( iunit, '(1x,i5,1x,i5,1x,e16.7)', iostat = ios ) ia(i), ja(i), a(i)

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'slap_FILE_READ - Fatal error!'
      write ( *, '(a,i8)' ) '  Error while reading matrix element ', i+1
      stop 1
    end if

  end do
!
!  If a value for RHS is available in the file, read it in.
!
  if ( irhs == 1 ) then

    if ( n_max < n ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'slap_FILE_READ - Fatal error!'
      write ( *, '(a)' ) '  N_MAX < N.'
      write ( *, '(a,i8)' ) '  N_MAX = ', n_max
      write ( *, '(a,i8)' ) '  N     = ', n
      stop 1
    end if

    call slap_rhs_read ( iunit, n, rhs, ios )

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'slap_FILE_READ - Fatal error!'
      write ( *, '(a)' ) '  Error while reading RHS from slap file.'
      stop 1
    end if

  end if
!
!  If a value of SOLN is available in the file, read it.
!
  if ( isoln == 1 ) then

    if ( n_max < n ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'slap_FILE_READ - Fatal error!'
      write ( *, '(a)' ) '  N_MAX < N.'
      write ( *, '(a,i8)' ) '  N_MAX = ', n_max
      write ( *, '(a,i8)' ) '  N     = ', n
      stop 1
    end if

    call slap_soln_read ( iunit, n, soln, ios )

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'slap_FILE_READ - Fatal error!'
      write ( *, '(a)' ) '  Error while reading SOLN from slap file.'
      stop 1
    end if

  end if

  return
end
subroutine slap_file_write ( n, nelt, isym, irhs, isoln, ia, ja, a, rhs, &
  soln, iunit )

!*****************************************************************************80
!
!! slap_FILE_WRITE writes out slap Triad Format Linear System.
!
!  Discussion:
!
!    This routine writes out a slap Triad Format Linear System.
!    including the matrix, right hand side, and solution to the
!    system, if known.
!
!    The format for the output is as follows.  On  the first line
!    are counters and flags:
!
!      N, NELT, ISYM, IRHS, ISOLN.
!
!    N, NELT and ISYM are described below.  IRHS is a flag indicating if
!    the RHS was written out (1 is  yes, 0 is  no).  ISOLN  is a
!    flag indicating if the SOLN was written out  (1 is yes, 0 is
!    no).  The format for the first line is: 5i10.  Then comes the
!    NELT Triad's IA(I), JA(I) and A(I), I = 1, NELT.  The format
!    for  these lines is   :  1X,I5,1X,I5,1X,E16.7.   Then comes
!    RHS(I), I = 1, N, if IRHS = 1.  Then comes SOLN(I), I  = 1,
!    N, if ISOLN = 1.  The format for these lines is: 1X,E16.7.
!
!    This routine requires that the  matrix A be stored in the
!    slap Triad format.  In this format only the non-zeros  are
!    stored.  They may appear in ANY order.  The user supplies
!    three arrays of length NELT, where NELT is the number of
!    non-zeros in the matrix: (IA(NELT), JA(NELT), A(NELT)).  For
!    each non-zero the user puts the row and column index of that
!    matrix element in the IA and JA arrays.  The value of the
!    non-zero matrix element is placed in the corresponding
!    location of the A array.   This is an extremely easy data
!    structure to generate.  On the other hand it is not too
!    efficient on vector computers for the iterative solution of
!    linear systems.  Hence,  slap changes this input data
!    structure to the slap Column format for the iteration (but
!    does not change it back).
!
!    Here is an example of the slap Triad storage format for a
!    5x5 Matrix.  Recall that the entries may appear in any order.
!
!        5x5 Matrix       slap Triad format for 5x5 matrix on left.
!                              1  2  3  4  5  6  7  8  9 10 11
!    |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21
!    |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2
!    | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1
!    | 0  0  0 44  0|
!    |51  0 53  0 55|
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Mark Seager,
!    A SLAP for the Masses,
!    Lawrence Livermore National Laboratory,
!    Technical Report UCRL-100267, December 1988.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer NELT, the number of non-zeros stored in A.
!
!    Input, integer ISYM, indicates symmetric storage format.
!    * 0, all nonzero entries of the matrix are stored.
!    * 1, the matrix is symmetric, and only the lower triangle of
!      the matrix is stored.
!
!    Input, integer IRHS, is 1 if a right hand side vector 
!    is included.
!
!    Input, integer ISOLN, is 1 if a solution vector is included.
!
!    Input, integer IA(NELT), integer JA(NELT),
!    real ( kind = rk ) A(NELT), the slap triad matrix description.
!
!    Input, real ( kind = rk ) RHS(N), the right hand side vector.  This array is
!    accessed if JOB is set to print it out.
!
!    Input, real ( kind = rk ) SOLN(N), the solution to the linear system, 
!    if known.  This array is accessed if and only if JOB is set to print it.
!
!    Input, integer IUNIT, the FORTRAN device unit number to which
!    the matrix information is to be written.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nelt

  real ( kind = rk ) a(nelt)
  integer i
  integer ia(nelt)
  integer irhs
  integer isoln
  integer isym
  integer iunit
  integer ja(nelt)
  real ( kind = rk ) rhs(n)
  real ( kind = rk ) soln(n)

  call slap_header_write ( iunit, n, nelt, isym, irhs, isoln )
!
!  Write the matrix non-zeros in Triad format.
!
  do i = 1, nelt
    write ( iunit, '(1x,i5,1x,i5,1x,e16.7)' ) ia(i), ja(i), a(i)
  end do
!
!  Write the right hand side.
!
  if ( irhs == 1 ) then
    call slap_rhs_write ( iunit, n, rhs )
  end if
!
!  Write the solution.
!
  if ( isoln == 1 ) then
    call slap_soln_write ( iunit, n, soln )
  end if

  return
end
subroutine slap_header_print ( n, nelt, isym, irhs, isoln )

!*****************************************************************************80
!
!! slap_HEADER_PRINT prints the header line of a slap file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Mark Seager,
!    A SLAP for the Masses,
!    Lawrence Livermore National Laboratory,
!    Technical Report UCRL-100267, December 1988.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer NELT, the number of non-zeros stored in A.
!
!    Input, integer ISYM, indicates symmetric storage format.
!    * 0, all nonzero entries of the matrix are stored.
!    * 1, the matrix is symmetric, and only the lower triangle of
!      the matrix is stored.
!
!    Input, integer IRHS, is 1 if a right hand side vector 
!    is included.
!
!    Input, integer ISOLN, is 1 if a solution vector is included.
!
  implicit none

  integer isoln
  integer isym
  integer irhs
  integer n
  integer nelt

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  slap Sparse Matrix File Header:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  N,     the matrix order =                ', n
  write ( *, '(a,i8)' ) '  NELT,  the number of nonzeros stored =   ', nelt
  write ( *, '(a,i8)' ) '  ISYM,  1 if symmetric storage used =     ', isym
  write ( *, '(a,i8)' ) '  IRHS,  1 if a right hand side included = ', irhs
  write ( *, '(a,i8)' ) '  ISOLN, 1 if a solution vector included = ', isoln

  return
end
subroutine slap_header_read ( iunit, n, nelt, isym, irhs, isoln, ios )

!*****************************************************************************80
!
!! slap_HEADER_READ reads the header line from a slap file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Mark Seager,
!    A SLAP for the Masses,
!    Lawrence Livermore National Laboratory,
!    Technical Report UCRL-100267, December 1988.
!
!  Parameters:
!
!    Input, integer IUNIT, the FORTRAN device unit number to which
!    the matrix information is to be written.
!
!    Output, integer N, the order of the matrix.
!
!    Output, integer NELT, the number of non-zeros stored in A.
!
!    Output, integer ISYM, indicates symmetric storage format.
!    * 0, all nonzero entries of the matrix are stored.
!    * 1, the matrix is symmetric, and only the lower triangle of
!      the matrix is stored.
!
!    Output, integer IRHS, is 1 if a right hand side vector
!    is included.
!
!    Output, integer ISOLN, is 1 if a solution vector is included.
!
!    Output, integer IOS, the I/O status variable, which is 0 if
!    no I/O error occurred.
!
  implicit none

  integer ios
  integer isoln
  integer isym
  integer irhs
  integer iunit
  integer n
  integer nelt

  read ( iunit, *, iostat = ios ) n, nelt, isym, irhs, isoln

  return
end
subroutine slap_header_write ( iunit, n, nelt, isym, irhs, isoln )

!*****************************************************************************80
!
!! slap_HEADER_WRITE writes the header line to a slap file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Mark Seager,
!    A SLAP for the Masses,
!    Lawrence Livermore National Laboratory,
!    Technical Report UCRL-100267, December 1988.
!
!  Parameters:
!
!    Input, integer IUNIT, the FORTRAN device unit number to which
!    the matrix information is to be written.
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer NELT, the number of non-zeros stored in A.
!
!    Input, integer ISYM, indicates symmetric storage format.
!    * 0, all nonzero entries of the matrix are stored.
!    * 1, the matrix is symmetric, and only the lower triangle of
!      the matrix is stored.
!
!    Input, integer IRHS, is 1 if a right hand side is included.
!
!    Input, integer ISOLN, is 1 if a solution vector is included.
!
  implicit none

  integer isoln
  integer isym
  integer irhs
  integer iunit
  integer n
  integer nelt

  write ( iunit, '(5i10)' ) n, nelt, isym, irhs, isoln

  return
end
subroutine slap_rhs_print ( n, rhs )

!*****************************************************************************80
!
!! slap_RHS_PRINT prints the right hand side vector from a slap file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Mark Seager,
!    A SLAP for the Masses,
!    Lawrence Livermore National Laboratory,
!    Technical Report UCRL-100267, December 1988.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real ( kind = rk ) RHS(N), the right hand side vector to be written.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) rhs(n)

  call r8vec_print ( n, rhs, '  slap right hand side vector:' )

  return
end
subroutine slap_rhs_read ( iunit, n, rhs, ios )

!*****************************************************************************80
!
!! slap_RHS_READ reads the right hand side vector from a slap file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 January 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Mark Seager,
!    A SLAP for the Masses,
!    Lawrence Livermore National Laboratory,
!    Technical Report UCRL-100267, December 1988.
!
!  Parameters:
!
!    Input, integer IUNIT, the FORTRAN device unit number to which
!    the matrix information is to be written.
!
!    Input, integer N, the order of the matrix.
!
!    Output, real ( kind = rk ) RHS(N), the right hand side vector.
!
!    Output, integer IOS, the I/O status variable, which is 0 if
!    no I/O error occurred.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer ios
  integer iunit
  real ( kind = rk ) rhs(n)

  read ( iunit, *, iostat = ios ) rhs(1:n)

  return
end
subroutine slap_rhs_write ( iunit, n, rhs )

!*****************************************************************************80
!
!! slap_RHS_WRITE writes a right hand side vector to a slap file.
!
!  Modified:
!
!    31 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Mark Seager,
!    A SLAP for the Masses,
!    Lawrence Livermore National Laboratory,
!    Technical Report UCRL-100267, December 1988.
!
!  Parameters:
!
!    Input, integer IUNIT, the FORTRAN device unit number to which
!    the matrix information is to be written.
!
!    Input, integer N, the order of the matrix.
!
!    Input, real ( kind = rk ) RHS(N), the right hand side vector to be written.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer iunit
  real ( kind = rk ) rhs(n)

  write ( iunit, '(1x,e16.7)' ) rhs(1:n)

  return
end
subroutine slap_soln_print ( n, soln )

!*****************************************************************************80
!
!! slap_SOLN_PRINT prints the solution vector from a slap file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Mark Seager,
!    A SLAP for the Masses,
!    Lawrence Livermore National Laboratory,
!    Technical Report UCRL-100267, December 1988.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real ( kind = rk ) SOLN(N), the solution vector to be written.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) soln(n)

  call r8vec_print ( n, soln, '  slap solution vector:' )

  return
end
subroutine slap_soln_read ( iunit, n, soln, ios )

!*****************************************************************************80
!
!! slap_SOLN_READ reads the solution vector from a slap file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 January 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Mark Seager,
!    A SLAP for the Masses,
!    Lawrence Livermore National Laboratory,
!    Technical Report UCRL-100267, December 1988.
!
!  Parameters:
!
!    Input, integer IUNIT, the FORTRAN device unit number to which
!    the matrix information is to be written.
!
!    Input, integer N, the order of the matrix.
!
!    Output, real ( kind = rk ) SOLN(N), the solution vector.
!
!    Output, integer IOS, the I/O status variable, which is 0 if
!    no I/O error occurred.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer ios
  integer iunit
  real ( kind = rk ) soln(n)

  read ( iunit, *, iostat = ios ) soln(1:n)

  return
end
subroutine slap_soln_write ( iunit, n, soln )

!*****************************************************************************80
!
!! slap_SOLN_WRITE writes a solution vector to a slap file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Mark Seager,
!    A SLAP for the Masses,
!    Lawrence Livermore National Laboratory,
!    Technical Report UCRL-100267, December 1988.
!
!  Parameters:
!
!    Input, integer IUNIT, the FORTRAN device unit number to which
!    the matrix information is to be written.
!
!    Input, integer N, the order of the matrix.
!
!    Input, real ( kind = rk ) SOLN(N), the solution vector to be written.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer iunit
  real ( kind = rk ) soln(n)

  write ( iunit, '(2x,e16.7)' ) soln(1:n)

  return
end
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
subroutine r8sp_print ( m, n, nz_num, isym, row, col, a, title )

!*****************************************************************************80
!
!! R8SP_PRINT prints an R8SP matrix.
!
!  Discussion:
!
!    This version of R8SP_PRINT has been specifically modified to allow,
!    and correctly handle, the case in which a single matrix location
!    A(I,J) is referenced more than once by the sparse matrix structure.
!    In such cases, the routine prints out the sum of all the values.
!
!    The R8SP storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.
!
!    It is possible that a pair of indices (I,J) may occur more than
!    once.  Presumably, in this case, the intent is that the actual value
!    of A(I,J) is the sum of all such entries.  This is not a good thing
!    to do, but I seem to have come across this in MATLAB.
!
!    The R8SP format is used by CSPARSE ("sparse triplet"), slap/SLAP 
!    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns of 
!    the matrix.
!
!    Input, integer NZ_NUM, the number of nonzero elements in 
!    the matrix.
!
!    Input, integer ISYM, is 0 if the matrix is not symmetric, 
!    and 1 if the matrix is symmetric.  The symmetric case only makes sense
!    if the matrix is also square, that is, M = N.  In this case, only
!    the nonzeroes on the diagonal and in the lower triangle are stored.
!
!    Input, integer ROW(NZ_NUM), COL(NZ_NUM), the row and column 
!    indices of the nonzero elements.
!
!    Input, real ( kind = rk ) A(NZ_NUM), the nonzero elements of the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
! 
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nz_num

  real ( kind = rk ) a(nz_num)
  integer col(nz_num)
  integer isym
  integer m
  integer n
  integer row(nz_num)
  character ( len = * ) title

  call r8sp_print_some ( m, n, nz_num, isym, row, col, a, 1, 1, m, n, title )

  return
end
subroutine r8sp_print_some ( m, n, nz_num, isym, row, col, a, ilo, jlo, &
  ihi, jhi, title )

!*****************************************************************************80
!
!! R8SP_PRINT_SOME prints some of an R8SP matrix.
!
!  Discussion:
!
!    This version of R8SP_PRINT_SOME has been specifically modified to allow,
!    and correctly handle, the case in which a single matrix location
!    A(I,J) is referenced more than once by the sparse matrix structure.
!    In such cases, the routine prints out the sum of all the values.
!
!    The R8SP storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.
!
!    It is possible that a pair of indices (I,J) may occur more than
!    once.  Presumably, in this case, the intent is that the actual value
!    of A(I,J) is the sum of all such entries.  This is not a good thing
!    to do, but I seem to have come across this in MATLAB.
!
!    The R8SP format is used by CSPARSE ("sparse triplet"), slap/SLAP 
!    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns 
!    of the matrix.
!
!    Input, integer NZ_NUM, the number of nonzero elements
!    in the matrix.
!
!    Input, integer ISYM, is 0 if the matrix is not symmetric, 
!    and 1 if the matrix is symmetric.  The symmetric case only makes sense
!    if the matrix is also square, that is, M = N.  In this case, only
!    the nonzeroes on the diagonal and in the lower triangle are stored.
!
!    Input, integer ROW(NZ_NUM), COL(NZ_NUM), the row and column
!    indices of the nonzero elements.
!
!    Input, real ( kind = rk ) A(NZ_NUM), the nonzero elements of the matrix.
!
!    Input, integer ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: incx = 5
  integer nz_num

  real ( kind = rk ) a(nz_num)
  real ( kind = rk ) aij(incx)
  integer col(nz_num)
  integer i
  integer i2hi
  integer i2lo
  integer ihi
  integer ilo
  integer inc
  integer isym
  integer j
  integer j2
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  integer k
  integer m
  integer n
  integer row(nz_num)
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '
    write ( *, '(''  Col:  '',5(i7,7x))' ) ( j, j = j2lo, j2hi )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      aij(1:inc) = 0.0D+00
!
!  Is matrix entry K actually the value of A(I,J), with J2LO <= J <= J2HI?
!  Because MATLAB seems to allow for multiple (I,J,A) entries, we have
!  to sum up what we find.
! 
      do k = 1, nz_num

        if ( i == row(k) .and. &
             j2lo <= col(k) .and. &
             col(k) <= j2hi ) then 

          j2 = col(k) - j2lo + 1
          aij(j2) = aij(j2) + a(k)

        else if ( isym == 1 .and. &
                  m == n .and. &
                  i == col(k) .and. &
                  j2lo <= row(k) .and.  &
                  row(k) <= j2hi ) then

          j2 = row(k) - j2lo + 1
          aij(j2) = aij(j2) + a(k)
 
        end if

      end do

      if ( any ( aij(1:inc) /= 0.0D+00 ) ) then
        write ( *, '(i5,1x,5g14.6)' ) i, aij(1:inc)
      end if

    end do

  end do

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
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
    write ( *, '(2x,i8,g16.8)' ) i, a(i)
  end do

  return
end
subroutine r8vec_print_some ( n, a, max_print, title )

!*****************************************************************************80
!
!! R8VEC_PRINT_SOME prints "some" of an R8VEC.
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vector, is no more than MAX_PRINT, then
!    the entire vector is printed, one entry per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 December 2001
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
!    Input, integer MAX_PRINT, the maximum number of lines 
!    to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  integer i
  integer max_print
  character ( len = * ) title

  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  if ( n <= max_print ) then

    do i = 1, n
      write ( *, '(2x,i8,2x,g14.6)' ) i, a(i)
    end do

  else if ( 3 <= max_print ) then

    do i = 1, max_print-2
      write ( *, '(2x,i8,2x,g14.6)' ) i, a(i)
    end do
    write ( *, '(a)' ) '  ......  ..............'
    i = n
    write ( *, '(2x,i8,2x,g14.6)' ) i, a(i)

  else

    do i = 1, max_print - 1
      write ( *, '(2x,i8,2x,g14.6)' ) i, a(i)
    end do
    i = max_print
    write ( *, '(2x,i8,2x,g14.6,2x,a)' ) i, a(i), '...more entries...'

  end if

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
