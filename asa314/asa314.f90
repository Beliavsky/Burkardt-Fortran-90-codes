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
subroutine i4mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! I4MAT_PRINT prints an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4's.
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
!    An I4MAT is a rectangular array of I4's.
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
subroutine invmod ( mat, imat, rmod, cmod, nrow, ifault )

!*****************************************************************************80
!
!! INVMOD inverts a matrix using modulo arithmetic.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 May 2013
!
!  Author:
!
!    Original FORTRAN77 version by Roger Payne.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Roger Payne,
!    Inversion of matrices with contents subject to modulo arithmetic,
!    Applied Statistics,
!    Volume 46, Number 2, 1997, pages 295-298.
!
!  Parameters:
!
!    Input/output, integer MAT(NROW*NROW).
!    On input, the matrix to be inverted.
!    On output, the product of the input matrix and IMAT.
!
!    Output, integer IMAT(NROW*NROW), the inverse matrix.  
!    If IFAULT = -1 on output, then IMAT is only a left inverse.
!
!    Input, integer RMOD(NROW), the modulus for values in each row.
!
!    Input, integer CMOD(NROW), the modulus for values 
!    in each column.
!
!    Input, integer NROW, the order of the matrix.
!
!    Output, integer IFAULT, an error flag.
!    0, no error was detected.
!    -1, only a left inverse could be formed.
!    1, the matrix contains elements that are negative, or too large.
!    2, the matrix contains nonzero elements in mixed modulus positions.
!    3, the matrix cannot be inverted.
!
  implicit none

  integer nrow

  logical all_zero
  integer cmod(nrow)
  integer csort(nrow)
  integer i
  integer ifault
  integer imat(nrow*nrow)
  integer ir
  integer j
  integer k
  integer kir
  integer kjr
  integer mat(nrow*nrow)
  integer n
  integer rmod(nrow)
  integer rsort(nrow)
!
!  Check that elements in 'mixed-moduli' positions are all zero.
!
  n = 0
  do i = 1, nrow
    do j = 1, nrow

      n = n + 1

      if ( ( rmod(i) /= cmod(j) ) .and. ( 0 < mat(n) ) ) then
        ifault = 2
        return
      end if

      if ( ( rmod(i) < mat(n) ) .or. ( mat(n) < 0 ) ) then
        ifault = 1
        return
      end if

    end do
  end do

  n = 0
  do i = 1, nrow
    do j = 1, nrow
      n = n + 1
      imat(n) = 0
    end do
  end do
!
!  Sort rows and columns into ascending order of moduli
!
  call msort ( mat, imat, rmod, cmod, rsort, csort, nrow )
!
!  Complete initialization of inverse matrix 
!
  do n = 1, nrow * nrow, nrow + 1
    imat(n) = 1
  end do
!
!  Invert the matrix.
!
  do ir = 1, nrow

    kir = ( ir - 1 ) * nrow

    if ( mat(kir+ir) == 0 ) then
!
!  Find a row JR below IR such that K(JR,IR)>0
!
      all_zero = .true.

      do kjr = kir + nrow + ir, nrow * nrow, nrow
        if ( 0 < mat(kjr) ) then
          all_zero = .false.
          exit
        end if
      end do
!
!  Column IR contains all zeros in rows IR or below:
!  look for a row above with zeros to left of column IR 
!  and K(JR,IR)>0
!
      if ( all_zero ) then
        do kjr = ir, kir, nrow
          if ( 0 < mat(kjr) ) then
            do i = kjr - ir + 1, kjr - 1
              if ( 0 < mat(i) ) then
                ifault = 3
                return
              end if
            end do
            all_zero = .false.
            exit
          end if
        end do
      end if
!
!  Column IR contains all zeros
!
      if ( all_zero ) then
        cycle
      end if
!
!  Switch row JR with row IR
!
      kjr = kjr - ir

      do i = 1, nrow

        k = mat(kir+i)
        mat(kir+i) = mat(kjr+i)
        mat(kjr+i) = k

        k = imat(kir+i)
        imat(kir+i) = imat(kjr+i)
        imat(kjr+i) = k

      end do

    end if
!
!  Find a multiplier N such that N*MAT(IR,IR)=1 mod(P{IR})
!
    k = mat(kir+ir)
    do n = 1, rmod(ir) - 1
      if ( mod ( n * k, rmod(ir) ) == 1 ) then
        exit
      end if
    end do
!
!  Multiply row IR by N.
!
    if ( 1 < n ) then
      do i = kir + 1, ir * nrow
        mat(i) = mat(i) * n
        imat(i) = imat(i) * n
      end do
    end if
!
!  Subtract MAT(JR,IR) * row IR from each row JR
!
    do kjr = 0, nrow * nrow - 1, nrow
      n = rmod(ir) - mat(kjr+ir)
      if ( ( kjr /= kir ) .and. ( n /= 0 ) ) then
        do i = 1, nrow
          mat(kjr+i)  = mod (  mat(kjr+i) + n *  mat(kir+i), cmod(i) )
          imat(kjr+i) = mod ( imat(kjr+i) + n * imat(kir+i), cmod(i) )
        end do
      end if
    end do

  end do
!
!  Check inversion was possible - that result has
!  non-zero elements only on diagonal.
!
  ifault = 0
!
!  If we encounter a zero diagonal element, then only a left inverse
!  will be formed.
!
  do n = 1, nrow * nrow, nrow + 1
    if ( mat(n) == 0 ) then
      ifault = -1
    end if
    mat(n) = - mat(n)
  end do

  do n = 1, nrow * nrow
    if ( 0 < mat(n) ) then
      ifault = 3
      return
    end if
  end do

  do n = 1, nrow * nrow, nrow + 1
    mat(n) = - mat(n)
  end do
!
!  Unsort the rows and columns back into their original order.
!
  call musort ( mat, imat, rmod, cmod, rsort, csort, nrow )

  return
end
subroutine msort ( mat, imat, rmod, cmod, rsort, csort, nrow )

!*****************************************************************************80
!
!! MSORT sorts matrix rows and columns in ascending order of moduli.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 May 2013
!
!  Author:
!
!    Original FORTRAN77 version by Roger Payne.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Roger Payne,
!    Inversion of matrices with contents subject to modulo arithmetic,
!    Applied Statistics,
!    Volume 46, Number 2, 1997, pages 295-298.
!
!  Parameters:
!
!    Input/output, integer MAT(NROW*NROW).
!    On output, the matrix has been sorted.
!
!    Ignoreput, integer IMAT(NROW*NROW).  
!    This quantity is ignored.
!
!    Input/output, integer RMOD(NROW), the modulus for values in 
!    each row.  On output, these have been rearranged according to the sorting.
!
!    Input/output, integer CMOD(NROW), the modulus for values in 
!    each column.  On output, these have been rearranged according to the 
!    sorting.
!
!    Output, integer RSORT(NROW), the sorted row indices.
!
!    Output, integer CSORT(NROW), the sorted column indices.
!
!    Input, integer NROW, the order of the matrix.
!
  implicit none

  integer nrow

  integer cmod(nrow)
  integer csort(nrow)
  integer i
  integer imat(nrow*nrow)
  integer irc
  integer j
  integer jrc
  integer kirc
  integer kjrc
  integer mat(nrow*nrow)
  integer p
  integer rmod(nrow)
  integer rsort(nrow)
!
!  Keep compiler happy by pretending to use an entry of imat.
!
  call i4_fake_use ( imat(1) )
!
!  Initialize row and column addresses.
!
  do i = 1, nrow
    rsort(i) = i
    csort(i) = i
  end do
!
!  Sort the rows.
!
  do irc = 1, nrow
!
!  Find the next row.
!
    jrc = irc
    p = rmod(irc)

    do i = irc + 1, nrow
      if ( rmod(i) < p ) then
        p = rmod(i)
        jrc = i
      end if
    end do

    if ( irc /= jrc ) then

      i = rmod(irc)
      rmod(irc) = rmod(jrc)
      rmod(jrc) = i

      i = rsort(irc)
      rsort(irc) = rsort(jrc)
      rsort(jrc) = i
!
!  Switch the rows.
!
      kirc = ( irc - 1 ) * nrow
      kjrc = ( jrc - 1 ) * nrow

      do j = 1, nrow
        i = mat(kirc+j)
        mat(kirc+j) = mat(kjrc+j)
        mat(kjrc+j) = i
      end do

    end if

  end do
!
!  Sort the columns.
!
  do irc = 1, nrow
!
!  Find the next column.
!
    jrc = irc
    p = cmod(irc)

    do i = irc + 1, nrow
      if ( cmod(i) < p ) then
        p = cmod(i)
        jrc = i
      end if
    end do

    if ( irc /= jrc ) then

      i = cmod(irc)
      cmod(irc) = cmod(jrc)
      cmod(jrc) = i

      i = csort(irc)
      csort(irc) = csort(jrc)
      csort(jrc) = i
!
!  Switch the columns.
!
      do j = 0, nrow * nrow - 1, nrow
        i = mat(irc+j)
        mat(irc+j) = mat(jrc+j)
        mat(jrc+j) = i
      end do

    end if

  end do

  return
end
subroutine musort ( mat, imat, rmod, cmod, rsort, csort, nrow )

!*****************************************************************************80
!
!! MUSORT unsorts the inverse matrix rows and columns into the original order.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 May 2013
!
!  Author:
!
!    Original FORTRAN77 version by Roger Payne.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Roger Payne,
!    Inversion of matrices with contents subject to modulo arithmetic,
!    Applied Statistics,
!    Volume 46, Number 2, 1997, pages 295-298.
!
!  Parameters:
!
!    Input/output, integer MAT(NROW*NROW).
!    On output, the matrix has been "unsorted".
!
!    Input/output, integer IMAT(NROW*NROW).
!    On output, the matrix has been "unsorted".
!
!    Input/output, integer RMOD(NROW), the modulus for values in 
!    each row.  On output, these have been restored to their original ordering.
!
!    Input/output, integer CMOD(NROW), the modulus for values in
!    each column.  On output, these have been restored to their original 
!    ordering.
!
!    Input/output, integer RSORT(NROW), the sorted row indices.
!
!    Input/output, integer CSORT(NROW), the sorted column indices.
!
!    Input, integer NROW, the order of the matrix.
!
  implicit none

  integer nrow

  integer cmod(nrow)
  integer csort(nrow)
  integer i
  integer imat(nrow*nrow)
  integer irc
  integer j
  integer jrc
  integer kirc
  integer kjrc
  integer mat(nrow*nrow)
  integer rmod(nrow)
  integer rsort(nrow)
!
!  Sort rows of inverse (= columns of original).
!
  do irc = 1, nrow
!
!  Find next row.
!
    if ( csort(irc) /= irc ) then

      do jrc = irc + 1, nrow
        if ( csort(jrc) == irc ) then
          exit
        end if
      end do

      i = cmod(irc)
      cmod(irc) = cmod(jrc)
      cmod(jrc) = i

      i = csort(irc)
      csort(irc) = csort(jrc)
      csort(jrc) = i
!
!  Switch rows.
!
      kirc = ( irc - 1 ) * nrow
      kjrc = ( jrc - 1 ) * nrow

      do j = 1, nrow
        i = imat(kirc+j)
        imat(kirc+j) = imat(kjrc+j)
        imat(kjrc+j) = i
      end do

    end if

  end do
!
!  Sort the columns of the inverse (= rows of original).
!
  do irc = 1, nrow
!
!  Find the next column.
!
    if ( rsort(irc) /= irc ) then

      do jrc = irc + 1, nrow
        if ( rsort(jrc) == irc ) then
          exit
        end if
      end do

      i = rmod(irc)
      rmod(irc) = rmod(jrc)
      rmod(jrc) = i

      i = rsort(irc)
      rsort(irc) = rsort(jrc)
      rsort(jrc) = i
!
!  Switch the columns of IMAT.
!
      do j = 0, nrow * nrow - 1, nrow
        i = imat(irc+j)
        imat(irc+j) = imat(jrc+j)
        imat(jrc+j) = i
      end do
!
!  Switch the diagonal elements of MAT (others are zero).
!
      kirc = ( irc - 1 ) * nrow + irc
      kjrc = ( jrc - 1 ) * nrow + jrc

      i = mat(kirc)
      mat(kirc) = mat(kjrc)
      mat(kjrc) = i

    end if

  end do

  return
end

