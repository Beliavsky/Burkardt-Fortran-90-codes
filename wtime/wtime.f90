function wtime ( )

!*****************************************************************************80
!
!! WTIME returns a reading of the wall clock time.
!
!  Discussion:
!
!    To get the elapsed wall clock time, call WTIME before and after a given
!    operation, and subtract the first reading from the second.
!
!    This function is meant to suggest the similar routines:
!
!      "omp_get_wtime ( )" in OpenMP,
!      "MPI_Wtime ( )" in MPI,
!      and "tic" and "toc" in MATLAB.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) WTIME, the wall clock reading, in seconds.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer clock_max
  integer clock_rate
  integer clock_reading
  real ( kind = rk ) wtime

  call system_clock ( clock_reading, clock_rate, clock_max )

  wtime = real ( clock_reading, kind = rk ) &
        / real ( clock_rate, kind = rk )

  return
end
