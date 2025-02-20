subroutine i4mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! i4mat_print() prints an I4MAT.
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
!    28 August 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, the number of rows in A.
!
!    integer N, the number of columns in A.
!
!    integer A(M,N), the matrix to be printed.
!
!    character ( len = * ) TITLE, a title.
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
!! i4mat_print_some() prints some of an I4MAT.
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
!    28 August 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, N, the number of rows and columns.
!
!    integer A(M,N), an M by N matrix to be printed.
!
!    integer ILO, JLO, the first row and column to print.
!
!    integer IHI, JHI, the last row and column to print.
!
!    character ( len = * ) TITLE, a title.
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
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! i4vec_print() prints an I4VEC.
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
!    28 August 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of components of the vector.
!
!    integer A(N), the vector to be printed.
!
!    character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer n

  integer a(n)
  integer i
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,2x,i12)' ) i, a(i)
  end do

  return
end
subroutine rcont ( nrow, ncol, nrowt, ncolt, nsubt, matrix, key, ifault )

!*****************************************************************************80
!
!! rcont() generates a random two-way table with given marginal totals.
!
!  Discussion:
!
!    Each time the program is called, another table will be randomly
!    generated.
!
!    Note that it should be the case that the sum of the row totals
!    is equal to the sum of the column totals.  However, this program
!    does not check for that condition.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 August 2021
!
!  Author:
!
!    Original FORTRAN77 version by James Boyett
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    James Boyett,
!    Algorithm AS 144:
!    Random R x C Tables with Given Row and Column Totals,
!    Applied Statistics,
!    Volume 28, Number 3, pages 329-332, 1979.
!
!  Input:
!
!    integer NROW, the number of rows in the observed
!    matrix.
!
!    integer NCOL, the number of columns in the observed
!    matrix.
!
!    integer NROWT(NROW), the row totals of the observed matrix.
!
!    integer NCOLT(NCOL), the column totals of the
!    observed matrix.
!
!    integer NSUBT(NCOL), partial column sums.  Must not be changed 
!    by the calling program.
!
!    logical KEY, should be set to FALSE by the user before the initial call.  
!
!  Output:
!
!    integer NSUBT(NCOL), updated partial column sums.
!
!    integer MATRIX(NROW,NCOL), the random matrix.
!
!    logical KEY, reset it to TRUE, and should be left
!    at that value for subsequent calls in which the same values of NROW,
!    NCOL, NROWT and NCOLT are being used.
!
!    integer IFAULT, fault indicator.
!    0, no error occured.
!    1, NROW <= 0.
!    2, NCOL <= 1.
!    3, some entry of NROWT is less than 0.
!    4, some entry of NCOLT is less than 0.
!    5, the sample size (sum of the column totals) is too large.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ncol
  integer nrow

  integer, parameter :: nvec_max = 200

  integer i
  integer ifault
  integer ii
  integer j
  integer k
  logical key
  integer limit
  integer matrix(nrow,ncol)
  integer ncolt(ncol)
  integer nnvect(nvec_max)
  integer noct
  integer nrowt(nrow)
  integer nsubt(ncol)
  integer ntemp
  integer, save :: ntotal
  integer, save, dimension ( nvec_max ) :: nvect
  real ( kind = rk ) r

  ifault = 0

  if ( .not. key ) then
!
!  Set KEY for subsequent calls.
!
    key = .true.
!
!  Check for faults and prepare for future calls.
!
    if ( nrow <= 0 ) then
      ifault = 1
      return
    end if

    if ( ncol <= 1 ) then
      ifault = 2
      return
    end if

    do i = 1, nrow
      if ( nrowt(i) <= 0 ) then
        ifault = 3
        return
      end if
    end do

    if ( ncolt(1) <= 0 ) then
      ifault = 4
      return
    end if

    nsubt(1) = ncolt(1)

    do j = 2, ncol

      if ( ncolt(j) <= 0 ) then
        ifault = 4
        return
      end if

      nsubt(j) = nsubt(j-1) + ncolt(j)

    end do

    ntotal = nsubt(ncol)

    if ( nvec_max < ntotal ) then
      ifault = 5
      return
    end if
!
!  Initialize vector to be permuted.
!
    do i = 1, ntotal
      nvect(i) = i
    end do

  end if
!
!  Initialize vector to be permuted.
!
  nnvect(1:ntotal) = nvect(1:ntotal)
!
!  Permute the vector.
!
  ntemp = ntotal
  do i = 1, ntotal
    call random_number ( harvest = r )
    noct = int ( r * real ( ntemp, kind = rk ) + 1.0D+00 )
    nvect(i) = nnvect(noct)
    nnvect(noct) = nnvect(ntemp)
    ntemp = ntemp - 1
  end do
!
!  Construct the random matrix.
!
  matrix(1:nrow,1:ncol) = 0

  ii = 1

  do i = 1, nrow

    limit = nrowt(i)

    do k = 1, limit

      do j = 1, ncol
        if ( nvect(ii) <= nsubt(j) ) then
          ii = ii + 1
          matrix(i,j) = matrix(i,j) + 1
          exit
        end if
      end do

    end do

  end do

  return
end

