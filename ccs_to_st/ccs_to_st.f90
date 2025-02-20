subroutine ccs_print ( m, n, ncc, icc, ccc, acc, title )

!*****************************************************************************80
!
!! ccs_print() prints a sparse matrix in CCS format.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 July 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows in the matrix.
!
!    Input, integer N, the number of columns in the matrix.
!
!    Input, integer NCC, the number of CCS elements.
!
!    Input, integer ICC(NCC), the CCS rows.
!
!    Input, integer CCC(N+1), the compressed CCS columns.
!
!    Input, real ( kind = rk ) ACC(NCC), the CCS values.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer ncc

  real ( kind = rk ) acc(ncc)
  integer ccc(n+1)
  integer i
  integer icc(ncc)
  integer j
  integer k
  integer m
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) '     #     I     J         A'
  write ( *, '(a)' ) '  ----  ----  ----  ----------------'
  write ( *, '(a)' ) ' '

  if ( ccc(1) == 0 ) then

    j = 0
    do k = 1, ncc
      i = icc(k)
      do while ( ccc(j+2) <= k - 1 )
        j = j + 1
      end do
      write ( *, '(2x,i4,2x,i4,2x,i4,2x,g16.8)' ) k - 1, i, j, acc(k)
    end do

  else

    j = 1
    do k = 1, ncc
      i = icc(k)
      do while ( ccc(j+1) <= k )
        j = j + 1
      end do
      write ( *, '(2x,i4,2x,i4,2x,i4,2x,g16.8)' ) k, i, j, acc(k)
    end do

  end if

  return
end
subroutine ccs_to_st ( m, n, ncc, icc, ccc, acc, nst, ist, jst, ast )

!*****************************************************************************80
!
!! ccs_to_st() converts sparse matrix information from CCS to ST format.
!
!  Discussion:
!
!    Only JST actually needs to be computed.  The other three output 
!    quantities are simply copies.  
!    
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 July 2014
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
!    Input, integer NCC, the number of CCS elements.
!
!    Input, integer ICC(NCC), the CCS rows.
!
!    Input, integer CCC(N+1), the CCS compressed columns.
!
!    Input, real ( kind = rk ) ACC(NCC), the CCS values.
!
!    Output, integer NST, the number of ST elements.
!
!    Output, integer IST(NST), JST(NST), the ST rows and columns.
!
!    Output, real ( kind = rk ) AST(NST), the ST values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer ncc

  real ( kind = rk ) ast(ncc)
  real ( kind = rk ) acc(ncc)
  integer ccc(n+1)
  integer icc(ncc)
  integer ist(ncc)
  integer j
  integer jhi
  integer jlo
  integer jst(ncc)
  integer k
  integer khi
  integer klo
  integer m
  integer nst

  nst = 0

  if ( ccc(1) == 0 ) then

    jlo = 0
    jhi = n - 1
  
    do j = jlo, jhi

      klo = ccc(j+1)
      khi = ccc(j+2) - 1

      do k = klo, khi

        nst = nst + 1
        ist(nst) = icc(k+1)
        jst(nst) = j
        ast(nst) = acc(k+1)

      end do

    end do

  else

    jlo = 1
    jhi = n
  
    do j = jlo, jhi

      klo = ccc(j)
      khi = ccc(j+1) - 1

      do k = klo, khi

        nst = nst + 1
        ist(nst) = icc(k)
        jst(nst) = j
        ast(nst) = acc(k)

      end do

    end do

  end if

  return
end
subroutine st_print ( m, n, nst, ist, jst, ast, title )

!*****************************************************************************80
!
!! st_print() prints a sparse matrix in ST format.
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
!    Input, integer M, the number of rows.
!
!    Input, integer N, the number of columns.
!
!    Input, integer NST, the number of ST elements.
!
!    Input, integer IST(NST), JST(NST), the ST rows and columns.
!
!    Input, real ( kind = rk ) AST(NST), the ST values.
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
  write ( *, '(a)' ) '     #     I     J       A'
  write ( *, '(a)' ) '  ----  ----  ----  --------------'
  write ( *, '(a)' ) ' '

  do k = 1, nst
    write ( *, '(2x,i4,2x,i4,2x,i4,2x,g16.8)' ) k, ist(k), jst(k), ast(k)
  end do

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! timestamp() prints the current YMDHMS date as a time stamp.
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

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end

