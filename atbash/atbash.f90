subroutine atbash_encrypt ( plain, crypt )

!*****************************************************************************80
!
!! atbash_encrypt() encrypts a string using the ATBASH code.
!
!  Discussion:
!
!    Only the alphabetic characters will be encrypted.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 August 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    character ( len = * ) PLAIN, the string to be encrypted.
!
!  Output:
!
!    character ( len = * ) CRYPT, the encrypted version of the string.
!
  implicit none

  integer big_a
  integer big_z
  character ( len = * ) crypt
  integer i
  integer j
  integer n
  character ( len = * ) plain
  integer plain_i
  integer small_a
  integer small_z
!
!  How long is the plain text?
!
  n = len ( plain )
!
!  Initialize the encrypted string to the plain text.  This way,
!  nonalphabetic characters are passed through.
!
  crypt = plain
!
!  Determine the numeric equivalents of the endpoints of the A-Z and a-z intervals.
!
  small_a = ichar ( 'a' )
  small_z = ichar ( 'z' )
  big_a = ichar ( 'A' )
  big_z = ichar ( 'Z' )
!
!  Convert plain text characters to unsigned integers,
!  and reverse the index.
!
  do i = 1, n
    plain_i = ichar ( plain(i:i) )
    if ( small_a <= plain_i .and. plain_i <= small_z ) then
      j = small_a + small_z - plain_i
      crypt(i:i) = char ( j );
    else if ( big_a <= plain_i .and. plain_i <= big_z ) then
      j = big_a + big_z - plain_i
      crypt(i:i) = char ( j )
    end if
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
!    29 August 2021
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
 
