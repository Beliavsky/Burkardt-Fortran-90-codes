subroutine colored_noise_test01 ( n, q_d, alpha )

!*****************************************************************************80
!
!! colored_noise_test01() calls F_ALPHA with particular parameters.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements of the sequence 
!    to generate.
!
!    Input, real ( kind = rk ) Q_D, the variance of the sequence.
!
!    Input, real ( kind = rk ) ALPHA, the exponent of the power law.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer  n

  real ( kind = rk ) alpha
  integer i
  integer ios
  character ( len = 255 ) output_filename
  integer output_unit
  real ( kind = rk )  q_d
  real ( kind = rk ) x(n)

  write ( output_filename, '(a,f4.2,a)' ) "alpha_", alpha, '.txt'
!
!  Report parameters.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COLORED_NOISE_TEST01:'
  write ( *, '(a,i8,a)' ) '  Generating ', n, ' sample points.'
  write ( *, '(a,g14.6)' ) '  1/F^ALPHA noise has ALPHA = ', alpha
  write ( *, '(a,g14.6)' ) '  Variance is ', q_d

  call f_alpha ( n, q_d, alpha, x )
!
!  Print no more than 10 entries of the data.
!
  call r8vec_print_part ( n, x, 10, '  Noise sample:' )
!
!  Write the data to a file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COLORED_NOISE_TEST01 - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    stop
  end if

  do i = 1, n
    write ( output_unit, '(g14.6)' ) x(i)
  end do

  close ( unit = output_unit )

  write ( *, '(a)' ) '  Data written to file "' &
    // trim ( output_filename ) // '."'

  return
end
subroutine f_alpha ( n, q_d, alpha, x )

!*****************************************************************************80
!
!! F_ALPHA generates a 1/F^ALPHA noise sequence.
!
!  Discussion:
!
!    Thanks to Miro Stoyanov for pointing out that the second half of
!    the data returned by the inverse Fourier transform should be
!    discarded, 24 August 2010.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    24 August 2010
!
!  Author:
!
!    Original C version by Todd Walter.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jeremy Kasdin,
!    Discrete Simulation of Colored Noise and Stochastic Processes
!    and 1/f^a Power Law Noise Generation,
!    Proceedings of the IEEE,
!    Volume 83, Number 5, 1995, pages 802-827.
!
!  Parameters:
!
!    Input, integer N, the number of samples to generate.
!
!    Input, real ( kind = rk ) Q_D, the variance of the noise.
!
!    Input, real ( kind = rk ) ALPHA, the exponent for the noise.
!
!    Output, real ( kind = rk ) X(N), a sequence sampled with the given
!    power law.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) alpha
  real ( kind = rk ) h_a(n)
  real ( kind = rk ) h_azero
  real ( kind = rk ) h_b(n)
  real ( kind = rk ) hfa(2*n)
  integer i
  real ( kind = rk ) q_d
  real ( kind = rk ) r8_normal_01
  real ( kind = rk ) w_a(n)
  real ( kind = rk ) w_azero
  real ( kind = rk ) w_b(n)
  real ( kind = rk ) wfa(2*n)
  real ( kind = rk ) wi
  real ( kind = rk ) wr
  real ( kind = rk ) x(n)
  real ( kind = rk ) x2(2*n)
!
!  Set the deviation of the noise.
!
  q_d = sqrt ( q_d )
!
!  Generate the coefficients Hk.
!
  hfa(1) = 1.0D+00
  do i = 2, n
    hfa(i) = hfa(i-1) * ( 0.5D+00 * alpha + real ( i - 2, kind = rk ) ) &
      / ( real ( i - 1, kind = rk ) )
  end do
  hfa(n+1:2*n) = 0.0D+00
!
!  Fill Wk with white noise.
!
  do i = 1, n
    wfa(i) = q_d * r8_normal_01 ( )
  end do
  wfa(n+1:2*n) = 0.0D+00
!
!  Perform the discrete Fourier transforms of Hk and Wk.
!
  call r8vec_sftf ( 2 * n, hfa, h_azero, h_a, h_b )

  call r8vec_sftf ( 2 * n, wfa, w_azero, w_a, w_b )
!
!  Multiply the two complex vectors.
!
  w_azero = w_azero * h_azero

  do i = 1, n
    wr = w_a(i)
    wi = w_b(i)
    w_a(i) = wr * h_a(i) - wi * h_b(i)
    w_b(i) = wi * h_a(i) + wr * h_b(i)
  end do
!
!  This scaling is introduced only to match the behavior
!  of the Numerical Recipes code...
!
  w_azero = w_azero * real ( 2 * n, kind = rk )

  w_a(1:n-1) = w_a(1:n-1) * real ( n, kind = rk )
  w_b(1:n-1) = w_b(1:n-1) * real ( n, kind = rk )

  w_a(n) = w_a(n) * real ( 2 * n, kind = rk )
  w_b(n) = w_b(n) * real ( 2 * n, kind = rk )
!
!  Take the inverse Fourier transform of the result.
!
  call r8vec_sftb ( 2 * n, w_azero, w_a, w_b, x2 )
!
!  Only return the first N inverse Fourier transform values.
!
  x(1:n) = x2(1:n)

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IUNIT, the free unit number.
!
  implicit none

  integer i
  integer ios
  integer iunit
  logical lopen
 
  iunit = 0
 
  do i = 1, 99
 
    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )
 
      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if
 
  end do

  return
end
function r8_normal_01 ( )

!*****************************************************************************80
!
!! r8_normal_01() returns a unit pseudonormal R8.
!
!  Discussion:
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Output:
!
!    real ( kind = rk ) R8_NORMAL_01, a sample of the standard
!    normal PDF.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r1
  real ( kind = rk ) r2
  real ( kind = rk ) r8_normal_01
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) x

  call random_number ( harvest = r1 )
  call random_number ( harvest = r2 )

  x = sqrt ( - 2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * r8_pi * r2 )

  r8_normal_01 = x

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, real ( kind = rk ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  integer i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
  end do

  return
end
subroutine r8vec_print_part ( n, a, max_print, title )

!*****************************************************************************80
!
!! R8VEC_PRINT_PART prints "part" of an R8VEC.
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vector, is no more than MAX_PRINT, then
!    the entire vector is printed, one entry per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries of the vector.
!
!    Input, real ( kind = rk ) A(N), the vector to be printed.
!
!    Input, integer MAX_PRINT, the maximum number of lines
!    to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  integer i
  integer max_print
  character ( len = * ) title

  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  if ( n <= max_print ) then

    do i = 1, n
      write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)
    end do

  else if ( 3 <= max_print ) then

    do i = 1, max_print - 2
      write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)
    end do
    write ( *, '(a)' ) '  ......  ..............'
    i = n
    write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)

  else

    do i = 1, max_print - 1
      write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)
    end do
    i = max_print
    write ( *, '(2x,i8,a,1x,g14.6,2x,a)' ) i, ':', a(i), '...more entries...'

  end if

  return
end
subroutine r8vec_sftb ( n, azero, a, b, r )

!*****************************************************************************80
!
!! R8VEC_SFTB computes a "slow" backward Fourier transform of an R8VEC.
!
!  Discussion:
!
!    SFTB and SFTF are inverses of each other.  If we begin with data
!    AZERO, A, and B, and apply SFTB to it, and then apply SFTF to the
!    resulting R vector, we should get back the original AZERO, A and B.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data values.
!
!    Input, real ( kind = rk ) AZERO, the constant Fourier coefficient.
!
!    Input, real ( kind = rk ) A(N/2), B(N/2), the Fourier coefficients.
!
!    Output, real ( kind = rk ) R(N), the reconstructed data sequence.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n/2)
  real ( kind = rk ) azero
  real ( kind = rk ) b(n/2)
  integer i
  integer k
  real ( kind = rk ) r(n)
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) theta

  r(1:n) = azero
  do i = 1, n
    do k = 1, n / 2
      theta = real ( k * ( i - 1 ) * 2, kind = rk ) * r8_pi &
        / real ( n, kind = rk )
      r(i) = r(i) + a(k) * cos ( theta ) + b(k) * sin ( theta )
    end do
  end do

  return
end
subroutine r8vec_sftf ( n, r, azero, a, b )

!*****************************************************************************80
!
!! R8VEC_SFTF computes a "slow" forward Fourier transform of an R8VEC.
!
!  Discussion:
!
!    SFTF and SFTB are inverses of each other.  If we begin with data
!    R and apply SFTB to it, and then apply SFTB to the resulting AZERO, 
!    A, and B, we should get back the original R.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 July 2017
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data values.
!
!    Input, real ( kind = rk ) R(N), the data to be transformed.
!
!    Output, real ( kind = rk ) AZERO, = sum ( 1 <= I <= N ) R(I) / N.
!
!    Output, real ( kind = rk ) A(N/2), B(N/2), the Fourier coefficients.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(1:n/2)
  real ( kind = rk ) azero
  real ( kind = rk ) b(1:n/2)
  integer i
  integer j
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) r(n)
  real ( kind = rk ) theta

  azero = sum ( r(1:n) ) / real ( n, kind = rk )

  do i = 1, n / 2

    a(i) = 0.0D+00
    b(i) = 0.0D+00

    do j = 1, n
      theta = real ( 2 * i * ( j - 1 ), kind = rk ) * r8_pi &
        / real ( n, kind = rk )
      a(i) = a(i) + r(j) * cos ( theta )
      b(i) = b(i) + r(j) * sin ( theta )
    end do

    a(i) = a(i) / real ( n, kind = rk )
    b(i) = b(i) / real ( n, kind = rk )

    if ( mod ( n, 2 ) == 1 .or. i /= ( n / 2 ) ) then
      a(i) = 2.0D+00 * a(i)
      b(i) = 2.0D+00 * b(i)
    end if

  end do

  return
end
subroutine r8vec_sftf_test ( )

!*****************************************************************************80
!
!! R8VEC_SFTF_TEST tests R8VEC_SFTF.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 July 2017
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) azero
  real ( kind = rk ), allocatable :: a(:)
  real ( kind = rk ), allocatable :: b(:)
  integer i
  integer n
  integer nhalf
  real ( kind = rk ), allocatable :: r(:)
  real ( kind = rk ), allocatable :: r2(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8VEC_SFTF_TEST:'
  write ( *, '(a)' ) '  FORTRAN version'
  write ( *, '(a)' ) '  R8VEC_SFTF computes "slow" Fourier transform (forward)'
  write ( *, '(a)' ) '  of a vector of real data.'
  write ( *, '(a)' ) '  The original data can be recovered using R8VEC_SFTB.'

  n = 15

  allocate ( r(1:n) )
  call random_number ( harvest = r(1:n) )

  nhalf = ( n / 2 )
  allocate ( a(1:nhalf) )
  allocate ( b(1:nhalf) )
  call r8vec_sftf ( n, r, azero, a, b )
  
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Fourier coefficients:'
  write ( *, '(a)' ) ''
  write ( *, '(2x,f10.6)' ) azero
  do i = 1, nhalf
    write ( *, '(2x,f10.6,2x,f10.6)' ) a(i), b(i)
  end do

  allocate ( r2(1:n) )
  call r8vec_sftb ( n, azero, a, b, r2 )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Compare data R and recovered data R2:'
  write ( *, '(a)' ) ''
  do i = 1, n
    write ( *, '(2x,f10.6,2x,f10.6)' ) r(i), r2(i) 
  end do
!
!  Free memory.
!
  deallocate ( a )
  deallocate ( b )
  deallocate ( r )
  deallocate ( r2 )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8VEC_SFTF_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'

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
