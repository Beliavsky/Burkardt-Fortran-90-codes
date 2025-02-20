program main

!*****************************************************************************80
!
!! asa205_test() tests asa205().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA205_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the ASA205 library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA205_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 examines a simple case with no repeated sum values.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 3
  integer, parameter :: n = 4

  external eval01
  integer ifault
  integer, dimension ( n ) :: colsum = (/ 2, 3, 1, 4 /)
  integer, dimension ( m ) :: rowsum = (/ 5, 3, 2 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  The tables will not have any multiplicities.'

  call i4vec_print ( m, rowsum, '  The row sums:' )

  call i4vec_print ( n, colsum, '  The column sums:' )

  call enum ( m, n, rowsum, colsum, eval01, ifault )

  if ( ifault /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  ENUM returned error flag IFAULT = ', ifault
  end if

  return
end
subroutine eval01 ( iflag, table, m, n, rowsum, colsum, prob, mult )

!*****************************************************************************80
!
!! EVAL01 is called by ENUM every time a new contingency table is determined.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer IFLAG, input flag.
!    1, this is the first call.  No table is input.
!    2, this is a call with a new table.
!    3, this is the last call.  No table is input.
!
!    Input, integer TABLE(M,N), the current contingency table.
!
!    Input, integer M, the number of rows.
!
!    Input, integer N, the number of columns.
!
!    Input, integer ROWSUM(M), the row sums.
!
!    Input, integer COLSUM(N), the column sums.
!
!    Input, real PROB, the logarithm of the hypergeometric probability
!    of this table.
!
!    Input, integer MULT, the multiplicity of this table, that is,
!    the number of different tables that still have the same set of
!    entries, but differ by a permutation of some rows and columns.
!
  implicit none

  integer m
  integer n

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer colsum(n)
  integer, save :: count1 = 0
  integer, save :: count2 = 0
  integer iflag
  integer mult
  real ( kind = rk ) prob
  real ( kind = rk ), save :: psum = 0.0D+00
  integer rowsum(m)
  integer table(m,n)
!
!  First call, no table, initialize.
!
  if ( iflag == 1 ) then

    count1 = 0
    count2 = 0
    psum = 0.0D+00

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EVAL01'
    write ( *, '(a)' ) '  Only first ten cases will be printed.'
    write ( *, '(a)' ) ' '
!
!  Call with a new table.
!
  else if ( iflag == 2 ) then

    count1 = count1 + 1
    count2 = count2 + mult
    psum = psum + mult * exp ( prob )

    if ( count1 <= 10 ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EVAL01:'
      write ( *, '(i3,i3,g14.6)' ) count1, mult, prob

      call i4mat_print ( m, n, table, '  Table' )

    end if
!
!  Last call, no table.
!
  else if ( iflag == 3 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EVAL01 summary'
    write ( *, '(a,i8)' ) '  Number of cases (ignoring multiplicity):', count1
    write ( *, '(a,i8)' ) '  Number of cases (allowing multiplicity):', count2
    write ( *, '(a,g14.6)' ) '  Probability sum = ', psum

  end if

  call i4_fake_use ( colsum(1) )
  call i4_fake_use ( rowsum(1) )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 examines a case where a sum value is repeated.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 2
  integer, parameter :: n = 3

  external eval02
  integer ifault
  integer, dimension ( n ) :: colsum = (/ 1, 2, 1 /)
  integer, dimension ( m ) :: rowsum = (/ 3, 1 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  The data will have multiplicities.'

  call i4vec_print ( m, rowsum, '  The row sums:' )

  call i4vec_print ( n, colsum, '  The column sums:' )

  call enum ( m, n, rowsum, colsum, eval02, ifault )

  if ( ifault /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  ENUM returned error flag IFAULT = ', ifault
  end if

  return
end
subroutine eval02 ( iflag, table, m, n, rowsum, colsum, prob, mult )

!*****************************************************************************80
!
!! EVAL02 is called by ENUM every time a new contingency table is determined.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer IFLAG, input flag.
!    1, this is the first call.  No table is input.
!    2, this is a call with a new table.
!    3, this is the last call.  No table is input.
!
!    Input, integer TABLE(M,N), the current contingency table.
!
!    Input, integer M, the number of rows.
!
!    Input, integer N, the number of columns.
!
!    Input, integer ROWSUM(M), the row sums.
!
!    Input, integer COLSUM(N), the column sums.
!
!    Input, real PROB, the logarithm of the hypergeometric probability
!    of this table.
!
!    Input, integer MULT, the multiplicity of this table, that is,
!    the number of different tables that still have the same set of
!    entries, but differ by a permutation of some rows and columns.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer colsum(n)
  integer, save :: count1 = 0
  integer, save :: count2 = 0
  integer iflag
  integer mult
  real ( kind = rk ) prob
  real ( kind = rk ), save :: psum = 0.0D+00
  integer rowsum(m)
  integer table(m,n)
!
!  First call, no table, initialize.
!
  if ( iflag == 1 ) then

    count1 = 0
    count2 = 0
    psum = 0.0D+00
!
!  Call with a new table.
!
  else if ( iflag == 2 ) then

    count1 = count1 + 1
    count2 = count2 + mult
    psum = psum + mult * exp ( prob )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EVAL02:'
    write ( *, '(i3,i3,g14.6)' ) count1, mult, prob

    call i4mat_print ( m, n, table, '  Table' )
!
!  Last call, no table.
!
  else if ( iflag == 3 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EVAL02 summary'
    write ( *, '(a,i8)' ) '  Number of cases (ignoring multiplicity):', count1
    write ( *, '(a,i8)' ) '  Number of cases (allowing multiplicity):', count2
    write ( *, '(a,g14.6)' ) '  Probability sum = ', psum

  end if

  call i4_fake_use ( colsum(1) )
  call i4_fake_use ( rowsum(1) )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 examines a test case from the paper, with 118489 tables.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 4
  integer, parameter :: n = 3

  external eval03
  integer ifault
  integer, dimension ( n ) :: colsum = (/ 4, 57, 59 /)
  integer, dimension ( m ) :: rowsum = (/ 3, 38, 39, 40 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Big problem test from the paper.'

  call i4vec_print ( m, rowsum, '  The row sums:' )

  call i4vec_print ( n, colsum, '  The column sums:' )

  call enum ( m, n, rowsum, colsum, eval03, ifault )

  if ( ifault /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  ENUM returned error flag IFAULT = ', ifault
  end if

  return
end
subroutine eval03 ( iflag, table, m, n, rowsum, colsum, prob, mult )

!*****************************************************************************80
!
!! EVAL03 is called by ENUM every time a new contingency table is determined.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer IFLAG, input flag.
!    1, this is the first call.  No table is input.
!    2, this is a call with a new table.
!    3, this is the last call.  No table is input.
!
!    Input, integer TABLE(M,N), the current contingency table.
!
!    Input, integer M, the number of rows.
!
!    Input, integer N, the number of columns.
!
!    Input, integer ROWSUM(M), the row sums.
!
!    Input, integer COLSUM(N), the column sums.
!
!    Input, real PROB, the logarithm of the hypergeometric probability
!    of this table.
!
!    Input, integer MULT, the multiplicity of this table, that is,
!    the number of different tables that still have the same set of
!    entries, but differ by a permutation of some rows and columns.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer colsum(n)
  integer, save :: count1 = 0
  integer, save :: count2 = 0
  integer iflag
  integer mult
  real ( kind = rk ) prob
  real ( kind = rk ), save :: psum = 0.0D+00
  integer rowsum(m)
  integer table(m,n)
!
!  First call, no table, initialize.
!
  if ( iflag == 1 ) then

    count1 = 0
    count2 = 0
    psum = 0.0D+00
!
!  Call with a new table.
!
  else if ( iflag == 2 ) then

    count1 = count1 + 1
    count2 = count2 + mult
    psum = psum + mult * exp ( prob )
!
!  Last call, no table.
!
  else if ( iflag == 3 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EVAL03 summary'
    write ( *, '(a,i8)' ) '  Number of cases (ignoring multiplicity): ', count1
    write ( *, '(a,i8)' ) '  Number of cases (allowing multiplicity): ', count2
    write ( *, '(a,g14.6)' ) '  Probability sum = ', psum
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Result from paper:'
    write ( *, '(a,i8)' ) '  Number of cases (ignoring multiplicity): ', 118489

  end if

  call i4_fake_use ( colsum(1) )
  call i4_fake_use ( rowsum(1) )
  call i4_fake_use ( table(1,1) )
  
  return
end
subroutine i4_fake_use ( n )

!*****************************************************************************80
!
!! i4_fake_use pretends to use a variable.
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

