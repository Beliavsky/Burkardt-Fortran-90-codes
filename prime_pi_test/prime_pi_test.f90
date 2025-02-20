program main

!*****************************************************************************80
!
!! prime_pi_test() tests prime_pi().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 December 2022
!
!  Author:
!
!    John Burkardt
!
  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'prime_pi_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test prime_pi()'

  call prime_pi1_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'prime_pi_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  return
end
subroutine pi_values ( n_data, n, p )

!*****************************************************************************80
!
!! pi_values() returns values of the Pi function.
!
!  Discussion:
!
!    Pi[n] is the number of primes less than or equal to n.
!
!    In Mathematica, the function can be evaluated by:
!
!      PrimePi[n]
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 December 2022
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Input:
!
!    integer N_DATA.  The user sets N_DATA to 0 before the first call. 
!
!  Output:
!
!    integer N_DATA.  The routine increments N_DATA by 1,
!    and returns the corresponding data when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    integer N, the argument.
!
!    integer P, the value of the function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 21

  integer n
  integer n_data
  integer, save, dimension ( n_max ) :: n_vec = (/ &
            1, &
            2, &
            4, &
            8, &
           16, &
           32, &
           64, &
          128, &
          256, &
          512, &
         1024, &
         2048, &
         4096, &
         8192, &
        16384, &
        32768, &
        65536, &
       131072, &
       262144, &
       524288, &
      1048576 /)
  integer p
  integer, save, dimension ( n_max ) :: p_vec = (/ &
             0, &
             1, &
             2, &
             4, &
             6, &
            11, &
            18, &
            31, &
            54, &
            97, &
           172, &
           309, &
           564, &
          1028, &
          1900, &
          3512, &
          6542, &
         12251, &
         23000, &
         43390, &
         82025 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    p = 0
  else
    n = n_vec(n_data)
    p = p_vec(n_data)
  end if

  return
end
subroutine prime_pi1_test ( )

!*****************************************************************************80
!
!! prime_pi1_test() tests prime_pi1().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 December 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer n
  integer n_data
  integer pi1
  integer pi2
  integer prime_pi1

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'prime_pi1_test():'
  write ( *, '(a)' ) '  Test prime_pi1()'
  write ( *, '(a)' ) '           n          Pi(n)          prime_pi1(n)'
  write ( *, '(a)' ) ''

  n_data = 0

  do while ( .true. )

    call pi_values ( n_data, n, pi1 )

    if ( n_data == 0 ) then
      exit
    end if

    pi2 = prime_pi1 ( n )

    write ( *, '(i12,2x,i10,2x,i10)' ) n, pi1, pi2

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
