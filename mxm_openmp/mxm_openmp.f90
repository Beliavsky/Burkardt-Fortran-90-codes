program main

!*****************************************************************************80
!
!! mxm_openmp() times a matrix-matrix multiply using OpenMP.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    13 October 2011
!
!  Author:
!
!    John Burkardt
!
  use omp_lib

  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: a(:,:)
  real ( kind = rk ) angle
  real ( kind = rk ), allocatable :: b(:,:)
  real ( kind = rk ), allocatable :: c(:,:)
  integer i
  integer j
  integer k
  integer n
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) s
  integer thread_num
  real ( kind = rk ) wtime

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'mxm_openmp():'
  write ( *, '(a)' ) '  FORTRAN90/OpenMP version'
  write ( *, '(a)' ) '  Compute matrix product C = A * B using OpenMP.'

  thread_num = omp_get_max_threads ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) &
    '  The number of processors available = ', omp_get_num_procs ( )
  write ( *, '(a,i8)' ) '  The number of threads available    = ', thread_num

  n = 2000
  write ( *, '(a,i8)' ) '  The matrix order N                 = ', n
  allocate ( a(1:n,1:n) )
  allocate ( b(1:n,1:n) )
  allocate ( c(1:n,1:n) )
!
!  Loop 1: Evaluate A.
!
  s = 1.0D+00 / sqrt ( real ( n, kind = rk ) )

  wtime = omp_get_wtime ( )

!$omp parallel shared ( a, b, c, n, s ) private ( angle, i, j, k )

  !$omp do
  do i = 1, n
    do j = 1, n
      angle = 2.0D+00 * pi * ( i - 1 ) * ( j - 1 ) / real ( n, kind = rk )
      a(i,j) = s * ( sin ( angle ) + cos ( angle ) ) 
    end do
  end do
  !$omp end do
!
!  Loop 2: Copy A into B.
!
  !$omp do
  do i = 1, n
    do j = 1, n
      b(i,j) = a(i,j)
    end do
  end do
  !$omp end do
!
!  Loop 3: Compute C = A * B.
!
  !$omp do
  do i = 1, n
    do j = 1, n
      c(i,j) = 0.0D+00
      do k = 1, n
        c(i,j) = c(i,j) + a(i,k) * b(k,j)
      end do
    end do
  end do
  !$omp end do

!$omp end parallel

  wtime = omp_get_wtime ( ) - wtime
  write ( *, '(a,g14.6)' ) '  Elapsed seconds = ', wtime
  write ( *, '(a,g14.6)' ) '  C(100,100)  = ', c(100,100)
!
!  Free memory.
!
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'mxm_openmp():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
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
