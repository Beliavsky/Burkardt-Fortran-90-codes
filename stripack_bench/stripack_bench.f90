program main

!*****************************************************************************80
!
!! stripack_bench() times stripack() on random data sets of increasing size.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 September 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer bisect
  real ( kind = 8 ), allocatable :: dist(:)
  integer ier
  integer, allocatable :: lend(:)
  integer, allocatable :: list(:)
  integer lnew
  integer, allocatable :: lptr(:)
  integer, allocatable :: ltri(:,:)
  integer n
  integer, allocatable :: near(:)
  integer, allocatable :: next(:)
  integer, parameter :: nrow = 6
  integer nt
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) time3
  real ( kind = 8 ), allocatable :: xyz(:,:)

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'stripack_bench():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Time STRIPACK''s Delaunay triangularion calculation'
  write ( *, '(a)' ) '  for sets of N random nodes of increasing size.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '        N        Time (seconds)'
  write ( *, '(a)' ) ''
!
!  Allocate space for the nodes.
!
  n = 12
!
!  Originally, we wanted to go to level 10, but already level 8 takes
!  30 minutes.
!
  do bisect = 0, 8

    allocate ( xyz(1:3,1:n) )
!
!  Generate N random node values on the sphere.
!
    call uniform_on_sphere01_map ( 3, n, xyz )
!
!  Start the timer.
!
    call cpu_time ( time1 )
!
!  Compute the Delaunay triangulation.
!
    allocate ( dist(1:n) )
    allocate ( lend(1:n) )
    allocate ( list(1:6*(n-2)) )
    allocate ( lptr(1:6*(n-2)) )
    allocate ( ltri(1:nrow,1:2*(n-2)) )
    allocate ( near(1:n) )
    allocate ( next(1:n) )

    call trmesh ( n, xyz(1,:), xyz(2,:), xyz(3,:), list, lptr, lend, lnew, near, &
      next, dist, ier )

    if ( ier /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'stripack_bench(): Fatal error!'
      write ( *, '(a,i8)' ) '  TRMESH returns IER = ', ier
      stop 1
    end if
!
!  Create the triangle list.
!
    call trlist ( n, list, lptr, lend, nrow, nt, ltri, ier )

    if ( ier /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'stripack_bench(): Fatal error!'
      write ( *, '(a,i8)' ) '  TRLIST returns IER = ', ier
      stop 1
    end if
!
!  Stop the timer.
!
    call cpu_time ( time2 )

    time3 = time2 - time1

    write ( *, '(2x,i8,2x,g14.6)') n, time3

    deallocate ( dist )
    deallocate ( lend )
    deallocate ( list )
    deallocate ( lptr )
    deallocate ( ltri )
    deallocate ( near )
    deallocate ( next )
    deallocate ( xyz )

    n = n + ( 4 ** bisect ) * 30

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'stripack_bench():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine r8vec_normal_01 ( n, x )

!*****************************************************************************80
!
!! r8vec_normal_01() returns a unit pseudonormal R8VEC.
!
!  Discussion:
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!    The Box-Muller method is used, which is efficient, but
!    generates an even number of values each time.  On any call
!    to this routine, an even number of new values are generated.
!    Depending on the situation, one value may be left over.
!    In that case, it is saved for the next call.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of values desired.  If N is
!    negative, then the code will flush its internal memory; in particular,
!    if there is a saved value to be used on the next call, it is
!    instead discarded.
!
!  Output:
!
!    real ( kind = rk ) X(N), a sample of the standard normal PDF.
!
!  Local:
!
!    integer MADE, records the number of values that have
!    been computed.  On input with negative N, this value overwrites
!    the return value of N, so the user can get an accounting of
!    how much work has been done.
!
!    real ( kind = rk ) R(N+1), is used to store some uniform
!    random values.  Its dimension is N+1, but really it is only needed
!    to be the smallest even number greater than or equal to N.
!
!    integer SAVED, is 0 or 1 depending on whether there is a
!    single saved value left over from the previous call.
!
!    integer X_LO_INDEX, X_HI_INDEX, records the range of entries of
!    X that we need to compute.  This starts off as 1:N, but is adjusted
!    if we have a saved value that can be immediately stored in X(1),
!    and so on.
!
!    real ( kind = rk ) Y, the value saved from the previous call, if
!    SAVED is 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer m
  integer, save :: made = 0
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) r(n+1)
  integer, save :: saved = 0
  real ( kind = rk ) x(n)
  integer x_hi_index
  integer x_lo_index
  real ( kind = rk ), save :: y = 0.0D+00
!
!  I'd like to allow the user to reset the internal data.
!  But this won't work properly if we have a saved value Y.
!  I'm making a crock option that allows the user to signal
!  explicitly that any internal memory should be flushed,
!  by passing in a negative value for N.
!
  if ( n < 0 ) then
    n = made
    made = 0
    saved = 0
    y = 0.0D+00
    return
  else if ( n == 0 ) then
    return
  end if
!
!  Record the range of X we need to fill in.
!
  x_lo_index = 1
  x_hi_index = n
!
!  Use up the old value, if we have it.
!
  if ( saved == 1 ) then
    x(1) = y
    saved = 0
    x_lo_index = 2
  end if
!
!  Maybe we don't need any more values.
!
  if ( x_hi_index - x_lo_index + 1 == 0 ) then
!
!  If we need just one new value, do that here to avoid null arrays.
!
  else if ( x_hi_index - x_lo_index + 1 == 1 ) then

    call random_number ( harvest = r(1:2) )

    x(x_hi_index) = &
             sqrt ( - 2.0D+00 * log ( r(1) ) ) * cos ( 2.0D+00 * pi * r(2) )
    y =      sqrt ( - 2.0D+00 * log ( r(1) ) ) * sin ( 2.0D+00 * pi * r(2) )

    saved = 1

    made = made + 2
!
!  If we require an even number of values, that's easy.
!
  else if ( mod ( x_hi_index - x_lo_index + 1, 2 ) == 0 ) then

    m = ( x_hi_index - x_lo_index + 1 ) / 2

    call random_number ( harvest = r(1:2*m) )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * cos ( 2.0D+00 * pi * r(2:2*m:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * sin ( 2.0D+00 * pi * r(2:2*m:2) )

    made = made + x_hi_index - x_lo_index + 1
!
!  If we require an odd number of values, we generate an even number,
!  and handle the last pair specially, storing one in X(N), and
!  saving the other for later.
!
  else

    x_hi_index = x_hi_index - 1

    m = ( x_hi_index - x_lo_index + 1 ) / 2 + 1

    call random_number ( harvest = r(1:2*m) )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * cos ( 2.0D+00 * pi * r(2:2*m-2:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * sin ( 2.0D+00 * pi * r(2:2*m-2:2) )

    x(n) = sqrt ( -2.0D+00 * log ( r(2*m-1) ) ) &
      * cos ( 2.0D+00 * pi * r(2*m) )

    y = sqrt ( -2.0D+00 * log ( r(2*m-1) ) ) &
      * sin ( 2.0D+00 * pi * r(2*m) )

    saved = 1

    made = made + x_hi_index - x_lo_index + 2

  end if

  return
end
subroutine uniform_on_sphere01_map ( dim_num, n, x )

!*****************************************************************************80
!
!! uniform_on_sphere01_map() maps uniform points onto the unit sphere.
!
!  Discussion:
!
!    The sphere has center 0 and radius 1.
!
!    This procedure is valid for any spatial dimension DIM_NUM.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 November 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Russell Cheng,
!    Random Variate Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998, pages 168.
!
!    George Marsaglia,
!    Choosing a point from the surface of a sphere,
!    Annals of Mathematical Statistics,
!    Volume 43, Number 2, April 1972, pages 645-646.
!
!    Reuven Rubinstein,
!    Monte Carlo Optimization, Simulation, and Sensitivity
!    of Queueing Networks,
!    Krieger, 1992,
!    ISBN: 0894647644,
!    LC: QA298.R79.
!
!  Input:
!
!    integer DIM_NUM, the dimension of the space.
!
!    integer N, the number of points.
!
!  Output:
!
!    real X(DIM_NUM,N), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num
  integer n

  integer j
  real ( kind = rk ) norm
  real ( kind = rk ) x(dim_num,n)
!
!  Fill a matrix with normally distributed values.
!
  call r8vec_normal_01 ( dim_num * n, x )
!
!  Normalize each column.
!
  do j = 1, n
!
!  Compute the length of the vector.
!
    norm = sqrt ( sum ( x(1:dim_num,j)**2 ) )
!
!  Normalize the vector.
!
    x(1:dim_num,j) = x(1:dim_num,j) / norm

  end do

  return
end
