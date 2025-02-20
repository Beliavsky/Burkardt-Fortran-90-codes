program main

!*****************************************************************************80
!
!! quad_mpi() estimates an integral, with MPI for parallel processing.
!
!  Discussion:
!
!    Thanks to Chad Mitchell for pointing out that "status" must be a vector
!    with dimension MPI_STATUS_SIZE, 16 September 2011.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 September 2011
!
!  Author:
!
!    John Burkardt
!
  use mpi

  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) error
  integer error_flag
  real ( kind = rk ) exact
  external f
  real ( kind = rk ) f
  integer i
  integer master
  real ( kind = rk ) my_a
  real ( kind = rk ) my_b
  integer my_id
  integer my_n
  real ( kind = rk ) my_total
  integer p
  integer p_num
  integer n
  integer source
  integer status(MPI_STATUS_SIZE)
  integer tag
  integer target
  real ( kind = rk ) total
  real ( kind = rk ) wtime
  real ( kind = rk ) x

  a =  0.0D+00
  b = 10.0D+00
  n = 10000000
  exact = 0.49936338107645674464D+00

  master = 0

  call MPI_Init ( error_flag )

  call MPI_Comm_size ( MPI_COMM_WORLD, p_num, error_flag )

  call MPI_Comm_rank ( MPI_COMM_WORLD, my_id, error_flag )
!
!  Process 0 reads in the quadrature rule, and parcels out the
!  evaluation points among the processes.
!
  if ( my_id == 0 ) then
!
!  We want N to be the total number of evaluations.
!  If necessary, we adjust N to be divisible by the number of processors.
!
    my_n = n / ( p_num - 1 )
    n = ( p_num - 1 ) * my_n

    wtime = MPI_Wtime ( )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'quad_mpi():'
    write ( *, '(a)' ) '  FORTRAN90/MPI version'
    write ( *, '(a)' ) '  Estimate an integral of f(x) from A to B.'
    write ( *, '(a)' ) '  f(x) = 50 / (pi * ( 2500 * x * x + 1 ) )'
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  A        = ', a
    write ( *, '(a,g14.6)' ) '  B        = ', b
    write ( *, '(a,i12)' ) '  N        = ', n
    write ( *, '(a,g24.16)' ) '  Exact    = ', exact
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Use MPI to divide the computation among'
    write ( *, '(a)' ) '  multiple processes.'

  end if

  source = master

  call MPI_Bcast ( my_n, 1, MPI_INTEGER, source, MPI_COMM_WORLD, &
    error_flag )
!
!  Process 0 assigns each process a subinterval of [A,B].
!
  if ( my_id == 0 ) then

    do p = 1, p_num - 1

      my_a = ( real ( p_num - p,     kind = rk ) * a   &
             + real (         p - 1, kind = rk ) * b ) &
             / real ( p_num     - 1, kind = rk )

      target = p
      tag = 1
      call MPI_Send ( my_a, 1, MPI_DOUBLE_PRECISION, target, tag, &
        MPI_COMM_WORLD, error_flag )

      my_b = ( real ( p_num - p - 1, kind = rk ) * a   &
             + real (         p,     kind = rk ) * b ) &
             / real ( p_num     - 1, kind = rk )

      target = p
      tag = 2
      call MPI_Send ( my_b, 1, MPI_DOUBLE_PRECISION, target, tag, &
        MPI_COMM_WORLD, error_flag )

    end do

    total = 0.0D+00
    my_total = 0.0D+00
!
!  Processes receive MY_A, MY_B, and compute their part of the integral.
!
  else

    source = master
    tag = 1

    call MPI_Recv ( my_a, 1, MPI_DOUBLE_PRECISION, source, tag, &
      MPI_COMM_WORLD, status, error_flag )

    source = master
    tag = 2

    call MPI_Recv ( my_b, 1, MPI_DOUBLE_PRECISION, source, tag, &
      MPI_COMM_WORLD, status, error_flag )

    my_total = 0.0D+00
    do i = 1, my_n
      x = ( real ( my_n - i,     kind = rk ) * my_a   &
          + real (        i - 1, kind = rk ) * my_b ) &
          / real ( my_n     - 1, kind = rk )
      my_total = my_total + f ( x )
    end do

    my_total = ( my_b - my_a ) * my_total / real ( my_n, kind = rk )

    write ( *, '(a,i8,a,g14.6)' ) &
      '  Process ', my_id, ' contributes MY_TOTAL = ', my_total

  end if
!
!  Each process sends its value of MY_TOTAL to the master process, to
!  be summed in TOTAL.
!
  call MPI_Reduce ( my_total, total, 1, MPI_DOUBLE_PRECISION, &
    MPI_SUM, master, MPI_COMM_WORLD, error_flag )
!
!  Report the results.
!
  if ( my_id == master ) then

    error = abs ( total - exact )
    wtime = MPI_Wtime ( ) - wtime

    write ( *, '(a)' ) ' '
    write ( *, '(a,g24.16)' ) '  Estimate = ', total
    write ( *, '(a,g14.6)' ) '  Error    = ', error
    write ( *, '(a,g14.6)' ) '  Time     = ', wtime

  end if
!
!  Terminate MPI.
!
  call MPI_Finalize ( error_flag )
!
!  Terminate.
!
  if ( my_id == master ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'QUAD_MPI:'
    write ( *, '(a)' ) '  Normal end of execution.'
  end if

  stop
end
function f ( x )

!*****************************************************************************80
!
!! F evaluates the function.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, the argument.
!
!    Output, real ( kind = rk ) F, the function value.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) f
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) x

  f = 50.0D+00 / ( pi * ( 2500.0D+00 * x * x + 1.0D+00 ) )

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
