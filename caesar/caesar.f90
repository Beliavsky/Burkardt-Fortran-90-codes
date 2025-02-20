function i4_modp ( i, j )

!*****************************************************************************80
!
!! i4_modp() returns the nonnegative remainder of I4 division.
!
!  Discussion:
!
!    If
!      NREM = I4_MODP ( I, J )
!      NMULT = ( I - NREM ) / J
!    then
!      I = J * NMULT + NREM
!    where NREM is always nonnegative.
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
!
!    An I4 is an integer value.
!
!  Example:
!
!        I     J     MOD I4_MODP    Factorization
!
!      107    50       7       7    107 =  2 *  50 + 7
!      107   -50       7       7    107 = -2 * -50 + 7
!     -107    50      -7      43   -107 = -3 *  50 + 43
!     -107   -50      -7      43   -107 =  3 * -50 + 43
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer I, the number to be divided.
!
!    integer J, the number that divides I.
!
!  Output:
!
!    integer I4_MODP, the nonnegative remainder when I is
!    divided by J.
!
  implicit none

  integer i
  integer i4_modp
  integer j
  integer value

  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MODP - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal divisor J = ', j
    stop 1
  end if

  value = mod ( i, j )

  if ( value < 0 ) then
    value = value + abs ( j )
  end if

  i4_modp = value

  return
end
subroutine s_to_caesar ( s1, k, s2 )

!*****************************************************************************80
!
!! s_to_caesar() applies a Caesar shift cipher to a string.
!
!  Discussion:
!
!    The Caesar shift cipher incremented each letter by 1, with Z going to A.
!
!    This function can apply a Caesar shift cipher to a string of characters,
!    using an arbitrary shift K, which can be positive, negative or zero.
!
!    Letters A through Z will be shifted by K, mod 26.
!    Letters a through z will be shifted by K, mod 26.
!
!    call s_to_caesar ( s1, 1, s2 ) will apply the traditional Caesar cipher.
!    call s_to_caesar ( s2, -1, s3 ) will decipher the result.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    character ( len = * ) S1, a string of characters.
!
!    integer K, the increment.
!
!  Output:
!
!    character ( len = * ) S2, the string of enciphered characters.
!
  implicit none

  integer i
  integer i1
  integer i2
  integer i4_modp
  integer iacap
  integer ialow
  integer k
  character ( len = * ) s1
  integer s1_len
  character ( len = * ) s2

  iacap = iachar ( 'A' );
  ialow = iachar ( 'a' );

  s1_len = len_trim ( s1 )
  s2 = ''

  do i = 1, s1_len

    i1 = iachar ( s1(i:i) )

    if ( iacap <= i1 .and. i1 <= iacap + 25 ) then
      i2 = i4_modp ( i1 + k - iacap, 26 ) + iacap
      s2(i:i) = achar ( i2 )
    else if ( ialow <= i1 .and. i1 <= ialow + 25 ) then
      i2 = i4_modp ( i1 + k - ialow, 26 ) + ialow
      s2(i:i) = achar ( i2 )
    else
      s2(i:i) = s1(i:i)
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
!    04 September 2021
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
 
