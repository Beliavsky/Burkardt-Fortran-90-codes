program main

!*****************************************************************************80
!
!! mxm() computes a matrix-matrix product and reports the elapsed CPU time.
!
!  Discussion:
!
!    The multiplication carried out is
!
!      A(1:N1,1:N3) = B(1:N1,1:N2) * C(1:N2,1:N3)
!
!    where B and C are real double precision matrices whose entries
!    are assigned randomly.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Usage:
!
!    mxm n1 n2 n3
!
!  Parameters:
!
!    Command line argument, integer N1, N2, N3, defines the number of
!    rows and columns in the two matrices.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: a(:,:)
  integer arg_num
  real ( kind = rk ), allocatable :: b(:,:)
  real ( kind = rk ), allocatable :: c(:,:)
  real ( kind = rk ) cpu_seconds
  integer flop_count
  integer iarg
  integer iargc
  integer ierror
  integer last
  real ( kind = rk ) mflops
  integer n1
  integer n2
  integer n3
  integer seed
  character ( len = 255 ) string
  real ( kind = rk ) time_estimate

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MXM:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Compute matrix-matrix product A = B * C'
!
!  Get the number of command line arguments.
!
  arg_num = iargc ( )
!
!  Get N1.
!
  if ( 1 <= arg_num ) then
  
    iarg = 1
    call getarg ( iarg, string )
    call s_to_i4 ( string, n1, ierror, last )
    
  else
  
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter N1, the number of rows in B.'
    read ( *, * ) n1
    
  end if
!
!  Get N2.
!
  if ( 2 <= arg_num ) then
  
    iarg = 2
    call getarg ( iarg, string )
    call s_to_i4 ( string, n2, ierror, last )
    
  else
  
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter N2, the number of columns in B and rows in C.'
    read ( *, * ) n2
    
  end if
!
!  Get N3.
!
  if ( 3 <= arg_num ) then
  
    iarg = 3
    call getarg ( iarg, string )
    call s_to_i4 ( string, n3, ierror, last )
    
  else
  
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter N3, the number of columns in C.'
    read ( *, * ) n3
    
  end if
!
!  Record the amount of work.
!  Each of the N1 * N3 entries of A requires N2 multiplies and N2 adds.
!
  flop_count = 2 * n1 * n2 * n3

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a,i6)' ) '  Matrix B is ', n1, ' by ', n2
  write ( *, '(a,i6,a,i6)' ) '  Matrix C is ', n2, ' by ', n3
  write ( *, '(a,i6,a,i6)' ) '  Matrix A will be ', n1, ' by ', n3
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Number of floating point operations = ', flop_count
  time_estimate = real ( flop_count, kind = rk ) / 2.6E+09
  write ( *, '(a,g10.2,a)' ) &
    '  Estimated CPU time is ', time_estimate, ' seconds.'
!
!  Generate A, B, C.
!
  allocate ( a(1:n1,1:n3) )
  allocate ( b(1:n1,1:n2) )
  allocate ( c(1:n2,1:n3) )
!
!  Set B and C.
!
  seed = 1325
  call matgen ( n1, n2, seed, b )
  call matgen ( n2, n3, seed, c )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Method     Cpu Seconds       MegaFlopS'
  write ( *, '(a)' ) '  ------  --------------  --------------'
!
!  IJK
!
  a(1:n1,1:n3) = 0.0D+00
  call mxm_ijk ( n1, n2, n3, a, b, c, cpu_seconds )

  if ( 0.0D+00 < cpu_seconds ) then
    mflops = real ( flop_count, kind = rk ) / cpu_seconds / 1000000.0D+00
  else
    mflops = -1.0D+00
  end if

  write ( *, '(2x,a6,2x,g14.6,2x,g14.6)' ) 'IJK   ', cpu_seconds, mflops
!
!  IKJ
!
  a(1:n1,1:n3) = 0.0D+00
  call mxm_ikj ( n1, n2, n3, a, b, c, cpu_seconds )

  if ( 0.0D+00 < cpu_seconds ) then
    mflops = real ( flop_count, kind = rk ) / cpu_seconds / 1000000.0D+00
  else
    mflops = -1.0D+00
  end if

  write ( *, '(2x,a6,2x,g14.6,2x,g14.6)' ) 'IKJ   ', cpu_seconds, mflops
!
!  JIK
!
  a(1:n1,1:n3) = 0.0D+00
  call mxm_jik ( n1, n2, n3, a, b, c, cpu_seconds )

  if ( 0.0D+00 < cpu_seconds ) then
    mflops = real ( flop_count, kind = rk ) / cpu_seconds / 1000000.0D+00
  else
    mflops = -1.0D+00
  end if

  write ( *, '(2x,a6,2x,g14.6,2x,g14.6)' ) 'JIK   ', cpu_seconds, mflops
!
!  JKI
!
  a(1:n1,1:n3) = 0.0D+00
  call mxm_jki ( n1, n2, n3, a, b, c, cpu_seconds )

  if ( 0.0D+00 < cpu_seconds ) then
    mflops = real ( flop_count, kind = rk ) / cpu_seconds / 1000000.0D+00
  else
    mflops = -1.0D+00
  end if

  write ( *, '(2x,a6,2x,g14.6,2x,g14.6)' ) 'JKI   ', cpu_seconds, mflops
!
!  KIJ
!
  a(1:n1,1:n3) = 0.0D+00
  call mxm_kij ( n1, n2, n3, a, b, c, cpu_seconds )

  if ( 0.0D+00 < cpu_seconds ) then
    mflops = real ( flop_count, kind = rk ) / cpu_seconds / 1000000.0D+00
  else
    mflops = -1.0D+00
  end if

  write ( *, '(2x,a6,2x,g14.6,2x,g14.6)' ) 'KIJ   ', cpu_seconds, mflops
!
!  KJI
!
  a(1:n1,1:n3) = 0.0D+00
  call mxm_kji ( n1, n2, n3, a, b, c, cpu_seconds )

  if ( 0.0D+00 < cpu_seconds ) then
    mflops = real ( flop_count, kind = rk ) / cpu_seconds / 1000000.0D+00
  else
    mflops = -1.0D+00
  end if

  write ( *, '(2x,a6,2x,g14.6,2x,g14.6)' ) 'KJI   ', cpu_seconds, mflops
!
!  MATMUL
!
  call mxm_matmul ( n1, n2, n3, a, b, c, cpu_seconds )

  if ( 0.0D+00 < cpu_seconds ) then
    mflops = real ( flop_count, kind = rk ) / cpu_seconds / 1000000.0D+00
  else
    mflops = -1.0D+00
  end if

  write ( *, '(2x,a6,2x,g14.6,2x,g14.6)' ) 'MATMUL', cpu_seconds, mflops
!
!  Deallocate arrays.
!
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MXM:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine matgen ( m, n, seed, a )

!*****************************************************************************80
!
!! matgen() generates a random matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    08 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns 
!    of the matrix.
!
!    Input, integer SEED, a seed for the random number
!    generator.
!
!    Output, real ( kind = rk ) A(M,N), the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  integer i
  integer j
  integer seed
!
!  Set the matrix A.
!
  do j = 1, n
    do i = 1, m
      seed = mod ( ( 3125 * seed ), 65536 )
      a(i,j) = ( seed - 32768.0 ) / 16384.0
    end do
  end do

  return
end
subroutine mxm_ijk ( n1, n2, n3, a, b, c, cpu_seconds )

!*****************************************************************************80
!
!! mxm_ijk() computes A = B * C using DO I, DO J, DO K loops.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, N3, define the orders of the
!    matrices.
!
!    Output, real ( kind = rk ) A(N1,N3), the result matrix.
!
!    Input, real ( kind = rk ) B(N1,N2), C(N2,N3), the factor matrices.
!
!    Output, real ( kind = rk ) CPU_SECONDS, the elapsed CPU time.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n1
  integer n2
  integer n3

  real ( kind = rk ) a(n1,n3)
  real ( kind = rk ) b(n1,n2)
  real ( kind = rk ) c(n2,n3)
  real ( kind = rk ) cpu_seconds
  integer i
  integer j
  integer k
  real ( kind = rk ) time1
  real ( kind = rk ) time2

  call cpu_time ( time1 )

  do i = 1, n1
    do j = 1, n3
      do k = 1, n2
        a(i,j) = a(i,j) + b(i,k) * c(k,j)
      end do
    end do
  end do

  call cpu_time ( time2 )

  cpu_seconds = time2 - time1

  return
end
subroutine mxm_ikj ( n1, n2, n3, a, b, c, cpu_seconds )

!*****************************************************************************80
!
!! mxm_ikj() computes A = B * C using DO I, DO K, DO J loops.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, N3, define the orders of the
!    matrices.
!
!    Output, real ( kind = rk ) A(N1,N3), the result matrix.
!
!    Input, real ( kind = rk ) B(N1,N2), C(N2,N3), the factor matrices.
!
!    Output, real ( kind = rk ) CPU_SECONDS, the elapsed CPU time.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n1
  integer n2
  integer n3

  real ( kind = rk ) a(n1,n3)
  real ( kind = rk ) b(n1,n2)
  real ( kind = rk ) c(n2,n3)
  real ( kind = rk ) cpu_seconds
  integer i
  integer j
  integer k
  real ( kind = rk ) time1
  real ( kind = rk ) time2

  call cpu_time ( time1 )

  do i = 1, n1
    do k = 1, n2
      do j = 1, n3
        a(i,j) = a(i,j) + b(i,k) * c(k,j)
      end do
    end do
  end do

  call cpu_time ( time2 )

  cpu_seconds = time2 - time1

  return
end
subroutine mxm_jik ( n1, n2, n3, a, b, c, cpu_seconds )

!*****************************************************************************80
!
!! mxm_jik() computes A = B * C using DO J, DO I, DO K loops.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, N3, define the orders of the
!    matrices.
!
!    Output, real ( kind = rk ) A(N1,N3), the result matrix.
!
!    Input, real ( kind = rk ) B(N1,N2), C(N2,N3), the factor matrices.
!
!    Output, real ( kind = rk ) CPU_SECONDS, the elapsed CPU time.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n1
  integer n2
  integer n3

  real ( kind = rk ) a(n1,n3)
  real ( kind = rk ) b(n1,n2)
  real ( kind = rk ) c(n2,n3)
  real ( kind = rk ) cpu_seconds
  integer i
  integer j
  integer k
  real ( kind = rk ) time1
  real ( kind = rk ) time2

  call cpu_time ( time1 )

  do j = 1, n3
    do i = 1, n1
      do k = 1, n2
        a(i,j) = a(i,j) + b(i,k) * c(k,j)
      end do
    end do
  end do

  call cpu_time ( time2 )

  cpu_seconds = time2 - time1

  return
end
subroutine mxm_jki ( n1, n2, n3, a, b, c, cpu_seconds )

!*****************************************************************************80
!
!! mxm_jki() computes A = B * C using DO J, DO K, DO I loops.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, N3, define the orders of the
!    matrices.
!
!    Output, real ( kind = rk ) A(N1,N3), the result matrix.
!
!    Input, real ( kind = rk ) B(N1,N2), C(N2,N3), the factor matrices.
!
!    Output, real ( kind = rk ) CPU_SECONDS, the elapsed CPU time.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n1
  integer n2
  integer n3

  real ( kind = rk ) a(n1,n3)
  real ( kind = rk ) b(n1,n2)
  real ( kind = rk ) c(n2,n3)
  real ( kind = rk ) cpu_seconds
  integer i
  integer j
  integer k
  real ( kind = rk ) time1
  real ( kind = rk ) time2

  call cpu_time ( time1 )

  do j = 1, n3
    do k = 1, n2
      do i = 1, n1
        a(i,j) = a(i,j) + b(i,k) * c(k,j)
      end do
    end do
  end do

  call cpu_time ( time2 )

  cpu_seconds = time2 - time1

  return
end
subroutine mxm_kij ( n1, n2, n3, a, b, c, cpu_seconds )

!*****************************************************************************80
!
!! mxm_kij() computes A = B * C using DO K, DO I, DO J loops.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, N3, define the orders of the
!    matrices.
!
!    Output, real ( kind = rk ) A(N1,N3), the result matrix.
!
!    Input, real ( kind = rk ) B(N1,N2), C(N2,N3), the factor matrices.
!
!    Output, real ( kind = rk ) CPU_SECONDS, the elapsed CPU time.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n1
  integer n2
  integer n3

  real ( kind = rk ) a(n1,n3)
  real ( kind = rk ) b(n1,n2)
  real ( kind = rk ) c(n2,n3)
  real ( kind = rk ) cpu_seconds
  integer i
  integer j
  integer k
  real ( kind = rk ) time1
  real ( kind = rk ) time2

  call cpu_time ( time1 )

  do k = 1, n2
    do i = 1, n1
      do j = 1, n3
        a(i,j) = a(i,j) + b(i,k) * c(k,j)
      end do
    end do
  end do

  call cpu_time ( time2 )

  cpu_seconds = time2 - time1

  return
end
subroutine mxm_kji ( n1, n2, n3, a, b, c, cpu_seconds )

!*****************************************************************************80
!
!! mxm_kji() computes A = B * C using DO K, DO J, DO I loops.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, N3, define the orders of the
!    matrices.
!
!    Output, real ( kind = rk ) A(N1,N3), the result matrix.
!
!    Input, real ( kind = rk ) B(N1,N2), C(N2,N3), the factor matrices.
!
!    Output, real ( kind = rk ) CPU_SECONDS, the elapsed CPU time.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n1
  integer n2
  integer n3

  real ( kind = rk ) a(n1,n3)
  real ( kind = rk ) b(n1,n2)
  real ( kind = rk ) c(n2,n3)
  real ( kind = rk ) cpu_seconds
  integer i
  integer j
  integer k
  real ( kind = rk ) time1
  real ( kind = rk ) time2

  call cpu_time ( time1 )

  do k = 1, n2
    do j = 1, n3
      do i = 1, n1
        a(i,j) = a(i,j) + b(i,k) * c(k,j)
      end do
    end do
  end do

  call cpu_time ( time2 )

  cpu_seconds = time2 - time1

  return
end
subroutine mxm_matmul ( n1, n2, n3, a, b, c, cpu_seconds )

!*****************************************************************************80
!
!! mxm_matmul() computes A = B * C, using the FORTRAN90 MATMUL function.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 May 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, N3, define the orders of the
!    matrices.
!
!    Input, real ( kind = rk ) A(N1,N3), the result matrix.
!
!    Input, real ( kind = rk ) B(N1,N2), C(N2,N3), the factor matrices.
!
!    Output, real ( kind = rk ) CPU_SECONDS, the elapsed CPU time.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n1
  integer n2
  integer n3

  real ( kind = rk ) a(n1,n3)
  real ( kind = rk ) b(n1,n2)
  real ( kind = rk ) c(n2,n3)
  real ( kind = rk ) cpu_seconds
  real ( kind = rk ) time1
  real ( kind = rk ) time2

  call cpu_time ( time1 )

  a(1:n1,1:n3) = matmul ( b(1:n1,1:n2), c(1:n2,1:n3) )

  call cpu_time ( time2 )

  cpu_seconds = time2 - time1

  return
end
subroutine s_to_i4 ( s, ival, ierror, length )

!*****************************************************************************80
!
!! s_to_i4() reads an I4 from a string.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    28 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string to be examined.
!
!    Output, integer IVAL, the integer value read from the string.
!    If the string is blank, then IVAL will be returned 0.
!
!    Output, integer IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer LENGTH, the number of characters of S 
!    used to make IVAL.
!
  implicit none

  character c
  integer i
  integer ierror
  integer isgn
  integer istate
  integer ival
  integer length
  character ( len = *  ) s

  ierror = 0
  istate = 0
  isgn = 1
  ival = 0

  do i = 1, len_trim ( s )

    c = s(i:i)
!
!  Haven't read anything.
!
    if ( istate == 0 ) then

      if ( c == ' ' ) then

      else if ( c == '-' ) then
        istate = 1
        isgn = -1
      else if ( c == '+' ) then
        istate = 1
        isgn = + 1
      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read the sign, expecting digits.
!
    else if ( istate == 1 ) then

      if ( c == ' ' ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read at least one digit, expecting more.
!
    else if ( istate == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        ival = 10 * ival + ichar ( c ) - ichar ( '0' )
      else
        ival = isgn * ival
        length = i - 1
        return
      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( istate == 2 ) then
    ival = isgn * ival
    length = len_trim ( s )
  else
    ierror = 1
    length = 0
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
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
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
