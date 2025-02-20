subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! i4vec_print() prints an I4VEC.
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
!    04 November 2003
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
  character ( len = 8 ) ctemp(incx)
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

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
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

      write ( *, '(i5,1x,10a8)' ) i, ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine rcont2 ( nrow, ncol, nrowt, ncolt, key, matrix, ierror )

!*****************************************************************************80
!
!! rcont2() constructs a random two-way contingency table with given sums.
!
!  Discussion:
!
!    It is possible to specify row and column sum vectors which
!    correspond to no table at all.  As far as I can see, this routine does
!    not detect such a case.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 March 2009
!
!  Author:
!
!    Original FORTRAN77 version by WM Patefield.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    WM Patefield,
!    Algorithm AS 159:
!    An Efficient Method of Generating RXC Tables with
!    Given Row and Column Totals,
!    Applied Statistics,
!    Volume 30, Number 1, 1981, pages 91-97.
!
!  Parameters:
!
!    Input, integer NROW, NCOL, the number of rows and columns
!    in the table.  NROW and NCOL must each be at least 2.
!
!    Input, integer NROWT(NROW), NCOLT(NCOL), the row and column
!    sums.  Each entry must be positive.
!
!    Input/output, logical KEY, a flag that indicates whether data has
!    been initialized for this problem.  Set KEY = .FALSE. before the first
!    call.
!
!    Output, integer MATRIX(NROW,NCOL), the matrix.
!
!    Output, integer IERROR, an error flag, which is returned
!    as 0 if no error occurred.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: maxtot = 5000

  integer ncol
  integer nrow

  logical done1
  logical done2
  real ( kind = rk ), save, dimension ( maxtot+1 ) :: fact
  integer i
  integer ia
  integer iap
  integer ib
  integer ic
  integer id
  integer idp
  integer ie
  integer ierror
  integer igp
  integer ihp
  integer ii
  integer iip
  integer j
  integer jc
  integer jwork(ncol)
  logical key
  integer l
  logical lsm
  logical lsp
  integer m
  integer matrix(nrow,ncol)
  integer ncolt(ncol)
  integer nll
  integer nlm
  integer nlmp
  integer nrowt(nrow)
  integer nrowtl
  integer, save :: ntotal = 0
  real ( kind = rk ) r
  real ( kind = rk ) sumprb
  real ( kind = rk ) x
  real ( kind = rk ) y

  ierror = 0
!
!  On user's signal, set up the factorial table.
!
  if ( .not. key ) then

    key = .true.

    if ( nrow <= 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RCONT - Fatal error!'
      write ( *, '(a)' ) '  Input number of rows is less than 2.'
      ierror = 1
      return
    end if

    if ( ncol <= 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RCONT - Fatal error!'
      write ( *, '(a)' ) '  The number of columns is less than 2.'
      ierror = 2
      return
    end if

    do i = 1, nrow
      if ( nrowt(i) <= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'RCONT - Fatal error!'
        write ( *, '(a)' ) '  An entry in the row sum vector is not positive.'
        ierror = 3
        return
      end if
    end do

    do j = 1, ncol
      if ( ncolt(j) <= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'RCONT - Fatal error!'
        write ( *, '(a)' ) &
         '  An entry in the column sum vector is not positive.'
        ierror = 4
        return
      end if
    end do

    if ( sum ( ncolt(1:ncol) ) /= sum ( nrowt(1:nrow) ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RCONT - Fatal error!'
      write ( *, '(a)' ) &
        '  The row and column sum vectors do not have the same sum.'
      ierror = 6
      return
    end if

    ntotal = sum ( ncolt(1:ncol) )

    if ( maxtot < ntotal ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RCONT - Fatal error!'
      write ( *, '(a)' ) &
        '  The sum of the column sum vector entries is too large.'
      ierror = 5
      return
    end if
!
!  Calculate log-factorials.
!
    x = 0.0D+00
    fact(1) = 0.0D+00
    do i = 1, ntotal
      x = x + log ( real ( i, kind = rk ) )
      fact(i+1) = x
    end do

  end if
!
!  Construct a random matrix.
!
  jwork(1:ncol-1) = ncolt(1:ncol-1)

  jc = ntotal

  do l = 1, nrow - 1

    nrowtl = nrowt(l)
    ia = nrowtl
    ic = jc
    jc = jc - nrowtl

    do m = 1, ncol - 1

      id = jwork(m)
      ie = ic
      ic = ic - id
      ib = ie - ia
      ii = ib - id
!
!  Test for zero entries in matrix.
!
      if ( ie == 0 ) then
        ia = 0
        matrix(l,m:ncol) = 0
        exit
      end if
!
!  Generate a pseudo-random number.
!
      call random_number ( harvest = r )
!
!  Compute the conditional expected value of MATRIX(L,M).
!
      done1 = .false.

      do

        nlm = int ( &
          real ( ia * id, kind = rk ) / real ( ie , kind = rk ) + 0.5D+00 )

        iap = ia + 1
        idp = id + 1
        igp = idp - nlm
        ihp = iap - nlm
        nlmp = nlm + 1
        iip = ii + nlmp
        x = exp ( fact(iap) + fact(ib+1) + fact(ic+1) + fact(idp) - &
          fact(ie+1) - fact(nlmp) - fact(igp) - fact(ihp) - fact(iip) )

        if ( r <= x ) then
          exit
        end if

        sumprb = x
        y = x
        nll = nlm
        lsp = .false.
        lsm = .false.
!
!  Increment entry in row L, column M.
!
        do while ( .not. lsp )

          j = ( id - nlm ) * ( ia - nlm )

          if ( j == 0 ) then

            lsp = .true.

          else

            nlm = nlm + 1
            x = x * real ( j, kind = rk ) &
              / real ( nlm * ( ii + nlm ), kind = rk )
            sumprb = sumprb + x

            if ( r <= sumprb ) then
              done1 = .true.
              exit
            end if

          end if

          done2 = .false.

          do while ( .not. lsm )
!
!  Decrement the entry in row L, column M.
!
            j = nll * ( ii + nll )

            if ( j == 0 ) then
              lsm = .true.
              exit
            end if

            nll = nll - 1
            y = y * real ( j, kind = rk ) &
              / real ( ( id - nll ) * ( ia - nll ), kind = rk )
            sumprb = sumprb + y

            if ( r <= sumprb ) then
              nlm = nll
              done2 = .true.
              exit
            end if

            if ( .not. lsp ) then
              exit
            end if

          end do

          if ( done2 ) then
            exit
          end if

        end do

        if ( done1 ) then
          exit
        end if

        if ( done2 ) then
          exit
        end if

        call random_number ( harvest = r )
        r = sumprb * r

      end do

      matrix(l,m) = nlm
      ia = ia - nlm
      jwork(m) = jwork(m) - nlm

    end do

    matrix(l,ncol) = ia

  end do
!
!  Compute the last row.
!
  matrix(nrow,1:ncol-1) = jwork(1:ncol-1)
  matrix(nrow,ncol) = ib - matrix(nrow,ncol-1)

  return
end

