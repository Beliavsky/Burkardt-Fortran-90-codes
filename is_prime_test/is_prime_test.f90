program main

!*****************************************************************************80
!
!! is_prime_test() tests is_prime_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 December 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'is_prime_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test is_prime().'

  call is_prime1_test ( )
  call is_prime2_test ( )
  call is_prime3_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'is_prime_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine is_prime1_test ( )

!*****************************************************************************80
!
!! is_prime1_test() tests is_prime1().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 December 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 5 ) ch1
  character ( len = 5 ) ch2
  logical is_prime1
  integer n
  integer n_data
  logical tf1
  logical tf2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'is_prime1_test():'
  write ( *, '(a)' ) '  is_prime1(n) is true if n is a prime:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         n      is_prime(n)     is_prime1(n)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call is_prime_values ( n_data, n, tf1 )

    if ( n_data == 0 ) then
      exit
    end if

    tf2 = is_prime1 ( n )

    if ( tf1 ) then
      ch1 = 'true'
    else
      ch1 = 'false'
    end if

    if ( tf2 ) then
      ch2 = 'true'
    else
      ch2 = 'false'
    end if

    write ( *, '(2x,i10,2x,a,2x,a)' ) n, ch1, ch2

  end do

  return
end
subroutine is_prime2_test ( )

!*****************************************************************************80
!
!! is_prime2_test() tests is_prime2().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 December 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 5 ) ch1
  character ( len = 5 ) ch2
  logical is_prime2
  integer n
  integer n_data
  logical tf1
  logical tf2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'is_prime2_test():'
  write ( *, '(a)' ) '  is_prime2(n) is true if n is a prime:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         n      is_prime(n)     is_prime2(n)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call is_prime_values ( n_data, n, tf1 )

    if ( n_data == 0 ) then
      exit
    end if

    tf2 = is_prime2 ( n )

    if ( tf1 ) then
      ch1 = 'true'
    else
      ch1 = 'false'
    end if

    if ( tf2 ) then
      ch2 = 'true'
    else
      ch2 = 'false'
    end if

    write ( *, '(2x,i10,2x,a,2x,a)' ) n, ch1, ch2

  end do

  return
end
subroutine is_prime3_test ( )

!*****************************************************************************80
!
!! is_prime3_test() tests is_prime3().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 December 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 5 ) ch1
  character ( len = 5 ) ch2
  logical is_prime3
  integer n
  integer n_data
  logical tf1
  logical tf2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'is_prime3_test():'
  write ( *, '(a)' ) '  is_prime3(n) is true if n is a prime:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         n      is_prime(n)     is_prime3(n)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call is_prime_values ( n_data, n, tf1 )

    if ( n_data == 0 ) then
      exit
    end if

    tf2 = is_prime3 ( n )

    if ( tf1 ) then
      ch1 = 'true'
    else
      ch1 = 'false'
    end if

    if ( tf2 ) then
      ch2 = 'true'
    else
      ch2 = 'false'
    end if

    write ( *, '(2x,i10,2x,a,2x,a)' ) n, ch1, ch2

  end do

  return
end
subroutine is_prime_values ( n_data, n, tf )

!*****************************************************************************80
!
!! is_prime_values() returns some values of the is_prime() function.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 December 2022
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer n_data: The user sets N_DATA to 0 before the first call.
!
!  Output:
!
!    integer n_data: On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    integer n: an integer.
!
!    logical tf: true if n is a prime.
!
  implicit none

  integer, parameter :: n_max = 22

  integer n
  integer n_data
  integer, save, dimension ( n_max ) :: n_vec = (/ &
          1, &
          2, &
         12, &
          3, &
         91, &
         53, &
        437, &
        311, &
       1333, &
        719, &
      16483, &
       7919, &
     223609, &
      81799, &
     873599, &
     800573, &
    5693761, &
    7559173, &
   90166053, &
   69600977, &
    6110601, &
  145253029  /)
  logical tf
  logical, save, dimension ( n_max ) :: tf_vec = (/ &
     .false., &
     .true., &
     .false., &
     .true., & 
     .false., &
     .true., &
     .false., &
     .true., &
     .false., &
     .true., &
     .false., &
     .true., &
     .false., &
     .true., &
     .false., &
     .true., &
     .false., &
     .true., &
     .false., &
     .true., &
     .false., &
     .true. /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    tf = .false.
  else
    n = n_vec(n_data)
    tf = tf_vec(n_data)
  end if

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
!    15 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

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
