program main

!*****************************************************************************80
!
!! ASA172_TEST tests ASA172.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 July 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA172_TEST:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the ASA172 library.'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA172_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 compares indices computed by a triple loop.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 July 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: kdim = 3

  integer i
  integer i1
  integer i2
  integer i3
  integer ifault
  integer iprod(kdim)
  integer ivec(kdim)
  integer j
  integer jsub
  integer n
  integer, dimension ( kdim ) :: nr = (/ 3, 2, 4 /)
  logical qfor
  logical qind

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  SIMDO can convert between compressed and'
  write ( *, '(a)' ) '  vector indices representing a nested loop.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here, we set QFOR = FALSE, meaning we do'
  write ( *, '(a)' ) '  NOT want to convert from FORTRAN ordering'
  write ( *, '(a)' ) '  to lexical ordering.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here, we actually carry out a triple loop'
  write ( *, '(a)' ) '  list the indices, and then compare.'

  qfor = .false.
!
!  If QFOR is FALSE, then the definition of IPROD is reversed...
!
  iprod(1) = nr(kdim)
  do i = 2, kdim
    iprod(i) = iprod(i-1) * nr(kdim+1-i)
  end do

  n = iprod(kdim)
!
!  Carry out the nested loops, and use JSUB to count each iteration.
!  In the inmost loop, print JSUB and the corresponding (I1,I2,I3) vector.
!
  jsub = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  #1: Generate JSUB by counting as we DO the loops:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  DO I1 = 1, N1'
  write ( *, '(a)' ) '    DO I2 = 1, N2'
  write ( *, '(a)' ) '      DO I3 = 1, N3'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      JSUB            I1        I2        I3'
  write ( *, '(a)' ) ' '
  do i1 = 1, nr(1)
    ivec(1) = i1
    do i2 = 1, nr(2)
      ivec(2) = i2
      do i3 = 1, nr(3)
        ivec(3) = i3
        jsub = jsub + 1
        write ( *, '(2x,i8,6x,i8,2x,i8,2x,i8)' ) jsub, i1, i2, i3
      end do
    end do
  end do
!
!  Now for each value of JSUB, retrieve the corresponding index subscript.
!  In order to use the QFOR = .FALSE. switch, I apparently have to reverse
!  the sense of the NR vector!
!
  qind = .true.

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  #2: Loop on JSUB, retrieve loop indices'
  write ( *, '(a)' ) '      QIND = TRUE J ->I(J)'
  write ( *, '(a)' ) '      QFOR = FALSE'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      JSUB            I1        I2        I3'
  write ( *, '(a)' ) ' '

  do j = 1, n
    jsub = j
    call simdo ( qind, qfor, iprod, kdim, jsub, ivec, ifault )
    write ( *, '(2x,i8,6x,i8,2x,i8,2x,i8)' ) jsub, ( ivec(i), i = 1, kdim )
  end do
!
!  Carry out the nested loops, and DO NOT compute JSUB.
!  Have SIMDO determine JSUB.
!
  qind = .false.

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  #3: For any set of loop indices, retrieve JSUB'
  write ( *, '(a)' ) '      QIND = FALSE I(J) -> J'
  write ( *, '(a)' ) '      QFOR = FALSE'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      JSUB            I1        I2        I3'
  write ( *, '(a)' ) ' '
  do i1 = 1, nr(1)
    ivec(1) = i1
    do i2 = 1, nr(2)
      ivec(2) = i2
      do i3 = 1, nr(3)
        ivec(3) = i3
        call simdo ( qind, qfor, iprod, kdim, jsub, ivec, ifault )
        write ( *, '(2x,i8,6x,i8,2x,i8,2x,i8)' ) jsub, i1, i2, i3
      end do
    end do
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 compares indices computed by a triple loop.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 July 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: kdim = 3

  integer i
  integer i1
  integer i2
  integer i3
  integer ifault
  integer iprod(kdim)
  integer ivec(kdim)
  integer j
  integer jsub
  integer n
  integer, dimension ( kdim ) :: nr = (/ 3, 2, 4 /)
  logical qfor
  logical qind

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02:'
  write ( *, '(a)' ) '  SIMDO can convert between compressed and'
  write ( *, '(a)' ) '  vector indices representing a nested loop.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here, we set QFOR = TRUE, meaning we DO'
  write ( *, '(a)' ) '  want to convert from the FORTRAN '
  write ( *, '(a)' ) '  ordering to lexical convention.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here, we actually carry out a triple loop'
  write ( *, '(a)' ) '  list the indices, and then compare.'

  qfor = .true.

  iprod(1) = nr(1)
  do i = 2, kdim
    iprod(i) = iprod(i-1) * nr(i)
  end do

  n = iprod(kdim)
!
!  Carry out the nested loops, and use JSUB to count each iteration.
!  In the inmost loop, print JSUB and the corresponding (I1,I2,I3) vector.
!
  jsub = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  #1: Generate JSUB by counting as we do the loops.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  DO I3 = 1, N3'
  write ( *, '(a)' ) '    DO I2 = 1, N2'
  write ( *, '(a)' ) '      DO I1 = 1, N1'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      JSUB            I1        I2        I3'
  write ( *, '(a)' ) ' '
  do i3 = 1, nr(3)
    ivec(3) = i3
    do i2 = 1, nr(2)
      ivec(2) = i2
      do i1 = 1, nr(1)
        ivec(1) = i1
        jsub = jsub + 1
        write ( *, '(2x,i8,6x,i8,2x,i8,2x,i8)' ) jsub, i1, i2, i3
      end do
    end do
  end do
!
!  Reverse the order, so that the loop indices are generated in lexical order.
!
  qind = .true.

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  #2: Setting QFOR false means loop indices'
  write ( *, '(a)' ) '  are generated in lexical order.'
  write ( *, '(a)' ) '      QIND = TRUE J -> I(J)'
  write ( *, '(a)' ) '      QFOR = TRUE'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      JSUB            I1        I2        I3'
  write ( *, '(a)' ) ' '

  do j = 1, n
    jsub = j
    call simdo ( qind, qfor, iprod, kdim, jsub, ivec, ifault )
    write ( *, '(2x,i8,6x,i8,2x,i8,2x,i8)' )  jsub, ( ivec(i), i = 1, kdim )
  end do
!
!  Carry out the nested loops, and DO NOT compute JSUB.
!  Have SIMDO determine JSUB.
!
  qind = .false.

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  #3: For any set of loop indices, retrieve JSUB'
  write ( *, '(a)' ) '      QIND = FALSE I(J) -> J'
  write ( *, '(a)' ) '      QFOR = TRUE'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      JSUB            I1        I2        I3'
  write ( *, '(a)' ) ' '
  do i3 = 1, nr(3)
    ivec(3) = i3
    do i2 = 1, nr(2)
      ivec(2) = i2
      do i1 = 1, nr(1)
        ivec(1) = i1
        call simdo ( qind, qfor, iprod, kdim, jsub, ivec, ifault )
        write ( *, '(2x,i8,6x,i8,2x,i8,2x,i8)' ) jsub, i1, i2, i3
      end do
    end do
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
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

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

