program main

!*****************************************************************************80
!
!! test_nearest() compares the performance of nearest neighbor routines.
!
!  Discussion:
!
!    We are given R, a set of NR points in M dimensions.
!
!    We are given S, a set of NS points in M dimensions.
!
!    For each S(I) in S, we seek the index J of the point R(J)
!    which is nearest to S(I) over all points in R.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 December 2012
!
!  Author:
!
!    John Burkardt
!
!  Local:
!
!    integer M, the spatial dimension.
!
!    integer NR, the number of data points.
!
!    integer NS, the number of sample points.
!
!    real ( kind = rk8 ) R(M,NR), the data points.
!
!    real ( kind = rk8 ) RT(NR,M), the transposed data points.
!
!    real ( kind = rk8 ) S(M,NS), the sample points. 
!
!    real ( kind = rk8 ) ST(NS,M), the transposed sample points. 
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: m_test_num = 3
  integer, parameter :: n_test_num = 6

  integer m
  integer m_test
  integer, dimension ( m_test_num ) :: m_test_data = (/ 2, 4, 8 /)
  integer, allocatable :: nearest(:)
  integer nr
  integer, dimension ( n_test_num ) :: nr_test_data = (/ &
    1000000, 100000, 10000,  1000,    100,      10 /)
  integer ns
  integer, dimension ( n_test_num ) :: ns_test_data = (/ &
      10,    100,  1000, 10000, 100000, 1000000 /)
  real ( kind = rk8 ), allocatable :: r(:,:)
  real ( kind = rk8 ), allocatable :: s(:,:)
  integer seed
  real ( kind = rk8 ) t1
  real ( kind = rk8 ) t2
  integer test

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'test_nearest():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Consider various nearest neighbor algorithms.'

  do m_test = 1, m_test_num

    m = m_test_data(m_test);

    do test = 1, n_test_num

      nr = nr_test_data(test)
      allocate ( r(1:m,1:nr) )
      ns = ns_test_data(test)
      allocate ( s(1:m,1:ns) )
      allocate ( nearest(ns) )

      seed = 123456789
      call r8mat_uniform_01 ( m, ns, seed, s )
      call r8mat_uniform_01 ( m, nr, seed, r )

      write ( *, '(a)' ) ''
      write ( *, '(a,i8,a,i8,a,i8)' ) '  M = ', m, ' NR = ', nr, '  NS = ', ns

      call cpu_time ( t1 )
      call find_closest1 ( m, nr, r, ns, s, nearest )
      call cpu_time ( t2 )
      write ( *, '(a,g14.6,a,i8,a,i8)' ) &
        '  #1 time: ', t2 - t1, '  size = ', ns, '  i(1) = ', nearest(1)

      deallocate ( nearest )
      deallocate ( r )
      deallocate ( s )

    end do

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST_NEAREST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  return
end
subroutine find_closest1 ( m, nr, r, ns, s, nearest )

!*****************************************************************************80
!
!! FIND_CLOSEST1 finds the nearest R point to each S point.
!
!  Discussion:
!
!    We are given R, a set of NR points in M dimensions.
!
!    We are given S, a set of NS points in M dimensions.
!
!    For each S(I) in S, we seek the index J of the point R(J)
!    which is nearest to S(I) over all points in R.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 December 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer NR, the number of data points.
!
!    Input, real ( kind = rk8 ) R(M,NR), the data points.
!
!    Input, integer NS, the number of sample points.
!
!    Input, real ( kind = rk8 ) S(M,NS), the sample points.
!
!    Output, integer NEAREST(NS), the index of the nearest 
!    data point.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer m
  integer nr
  integer ns

  real ( kind = rk8 ) dist_sq
  real ( kind = rk8 ) dist_sq_min
  integer jr
  integer js
  integer nearest(ns)
  real ( kind = rk8 ) r(m,nr)
  real ( kind = rk8 ) s(m,ns)

  do js = 1, ns

    dist_sq_min = huge ( dist_sq_min )
    nearest(js) = -1

    do jr = 1, nr

      dist_sq = sum ( ( r(1:m,jr) - s(1:m,js) ) **2 )

      if ( dist_sq < dist_sq_min ) then
        dist_sq_min = dist_sq
        nearest(js) = jr
      end if

    end do

  end do

  return
end
subroutine r8mat_uniform_01 ( m, n, seed, r )

!*****************************************************************************80
!
!! R8MAT_UNIFORM_01 fills an R8MAT with unit pseudorandom numbers.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns in
!    the array.
!
!    Input/output, integer SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = rk8 ) R(M,N), the array of pseudorandom values.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer m
  integer n

  integer i
  integer, parameter :: i4_huge = 2147483647
  integer j
  integer k
  integer seed
  real ( kind = rk8 ) r(m,n)

  do j = 1, n

    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge
      end if

      r(i,j) = real ( seed, kind = rk8 ) * 4.656612875D-10

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
