subroutine bad_in_simplex01 ( dim_num, point_num, x )

!*****************************************************************************80
!
!! bad_in_simplex01() is a "bad" (nonuniform) sampling of the unit simplex.
!
!  Discussion:
!
!    The interior of the unit DIM_NUM-dimensional simplex is the set of
!    points X(1:DIM_NUM) such that each X(I) is nonnegative, and
!    sum(X(1:DIM_NUM)) <= 1.
!
!    Any point in the unit simplex CAN be chosen by this algorithm.
!
!    However, the points that are chosen tend to be clustered near
!    the centroid.
!
!    This routine is supplied as an example of "bad" sampling.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the dimension of the space.
!
!    Input, integer POINT_NUM, the number of points.
!
!    Output, real ( kind = rk ) X(DIM_NUM,POINT_NUM), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num
  integer point_num

  real ( kind = rk ) e(dim_num+1)
  real ( kind = rk ) e_sum
  integer j
  real ( kind = rk ) x(dim_num,point_num)

  do j = 1, point_num

    call random_number ( harvest = e(1:dim_num+1) )

    e_sum = sum ( e(1:dim_num+1) )

    e(1:dim_num+1) = e(1:dim_num+1) / e_sum
!
!  We may take the values E(1:DIM_NUM+1) as being the barycentric
!  coordinates of the point.
!
    x(1:dim_num,j) = e(1:dim_num)

  end do

  return
end
subroutine brownian ( dim_num, n, x )

!*****************************************************************************80
!
!! BROWNIAN creates Brownian motion points.
!
!  Discussion:
!
!    A starting point is generated at the origin.  The next point
!    is generated at a uniformly random angle and a (0,1) normally
!    distributed distance from the previous point.
!
!    It is up to the user to rescale the data, if desired.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the dimension of the space.
!
!    Input, integer N, the number of points.
!
!    Output, real ( kind = rk ) X(DIM_NUM,N), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num
  integer n

  real ( kind = rk ) direction(dim_num)
  integer j
  real ( kind = rk ) r
  real ( kind = rk ) r8_normal_01
  real ( kind = rk ) x(dim_num,n)
!
!  Initial point.
!
  x(1:dim_num,1) = 0.0D+00
!
!  Generate angles and steps.
!
  do j = 2, n

    r = r8_normal_01 ( )
    r = abs ( r )

    call direction_uniform_nd ( dim_num, direction )

    x(1:dim_num,j) = x(1:dim_num,j-1) + r * direction(1:dim_num)

  end do

  return
end
subroutine direction_uniform_nd ( dim_num, w )

!*****************************************************************************80
!
!! DIRECTION_UNIFORM_ND generates a random direction vector.
!
!  Discussion:
!
!    This is actually simply a random point on the unit sphere.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the dimension of the space.
!
!    Output, real ( kind = rk ) W(DIM_NUM), a random direction vector,
!    with unit norm.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num

  real ( kind = rk ) norm
  real ( kind = rk ) w(dim_num)
!
!  Sample the standard normal distribution.
!
  call r8vec_normal_01 ( dim_num, w )
!
!  Compute the length of the vector.
!
  norm = sqrt ( sum ( w(1:dim_num)**2 ) )
!
!  Normalize the vector.
!
  w(1:dim_num) = w(1:dim_num) / norm

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a vlue between 1 and 99 which
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
!    18 September 2005
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

  integer, parameter :: rk = kind ( 1.0D+00 )

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
subroutine grid_in_cube01 ( dim_num, n, center, r )

!*****************************************************************************80
!
!! GRID_IN_CUBE01 generates grid points in the unit hypercube.
!
!  Discussion:
!
!    N points are needed in an DIM_NUM dimensional space.
!
!    The points are to lie on a uniform grid of side N_SIDE.
!
!    Unless the N = N_SIDE**DIM_NUM for some N_SIDE, we can't use all the
!    points on a grid.  What we do is find the smallest N_SIDE
!    that's big enough, and randomly omit some points.
!
!    If N_SIDE is 4, then the choices in 1D are:
!
!    A: 0,   1/3, 2/3, 1
!    B: 1/5, 2/5, 3/5, 4/5
!    C: 0,   1/4, 2/4, 3/4
!    D: 1/4, 2/4, 3/4, 1
!    E: 1/8, 3/8, 5/8, 7/8
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer N, the number of points.
!
!    Input, integer CENTER, specifies the 1D grid centering:
!    1: first point is 0.0, last point is 1.0;
!    2: first point is 1/(N+1), last point is N/(N+1);
!    3: first point is 0, last point is (N-1)/N;
!    4: first point is 1/N, last point is 1;
!    5: first point is 1/(2*N), last point is (2*N-1)/(2*N);
!
!    Output, real ( kind = rk ) R(DIM_NUM,N), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num
  integer n

  integer center
  integer j
  integer n_grid
  integer n_side
  real ( kind = rk ) r( dim_num, n )
  integer rank
  integer rank_list(n)
  integer tuple(dim_num)
!
!  Find the dimension of the smallest grid with N points.
!
  call grid_side ( dim_num, n, n_side )
!
!  We need to select N points out of N_SIDE**DIM_NUM set.
!
  n_grid = n_side**dim_num
!
!  Generate a random subset of N items from a set of size N_GRID.
!
  call ksub_random2 ( n_grid, n, rank_list )
!
!  Must make one dummy call to TUPLE_NEXT_FAST with RANK = -1.
!
  rank = -1
  call tuple_next_fast ( n_side, dim_num, rank, tuple )
!
!  Now generate the appropriate indices, and "center" them.
!
  do j = 1, n

    rank = rank_list(j) - 1

    call tuple_next_fast ( n_side, dim_num, rank, tuple )

    if ( center == 1 ) then
      r(1:dim_num,j) = real (     tuple(1:dim_num) - 1, kind = rk ) &
                     / real (     n_side - 1,     kind = rk )
    else if ( center == 2 ) then
      r(1:dim_num,j) = real (     tuple(1:dim_num), kind = rk ) &
                     / real (     n_side + 1, kind = rk )
    else if ( center == 3 ) then
      r(1:dim_num,j) = real (     tuple(1:dim_num) - 1, kind = rk ) &
                     / real (     n_side, kind = rk )
    else if ( center == 4 ) then
      r(1:dim_num,j) = real (     tuple(1:dim_num), kind = rk ) &
                     / real (     n_side, kind = rk )
    else if ( center == 5 ) then
      r(1:dim_num,j) = real ( 2 * tuple(1:dim_num) - 1, kind = rk ) &
                     / real ( 2 * n_side, kind = rk )
    end if

  end do

  return
end
subroutine grid_side ( dim_num, n, n_side )

!*****************************************************************************80
!
!! GRID_SIDE finds the smallest grid containing at least N points.
!
!  Discussion:
!
!    Each coordinate of the grid will have N_SIDE distinct values.
!    Thus the total number of points in the grid is N_SIDE**DIM_NUM.
!    This routine seeks the smallest N_SIDE such that N <= N_SIDE**M.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer N, the number of points.
!
!    Output, integer N_SIDE, the length of one side of the
!    smallest grid in M dimensions that contains at least N points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num
  real ( kind = rk ) exponent
  integer n
  integer n_side

  if ( n <= 0 ) then
    n_side = 0
    return
  end if

  if ( dim_num <= 0 ) then
    n_side = -1
    return
  end if

  exponent = 1.0D+00 / real ( dim_num, kind = rk )

  n_side = int ( ( real ( n, kind = rk ) )**exponent )

  if ( n_side**dim_num < n ) then
    n_side = n_side + 1
  end if

  return
end
function i4_factorial ( n )

!*****************************************************************************80
!
!! I4_FACTORIAL computes the factorial N!
!
!  Discussion:
!
!    FACTORIAL ( N ) = PRODUCT ( 1 <= I <= N ) I
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the argument of the factorial function.
!    If N is less than 1, I4_FACTORIAL is returned as 1.
!
!    Output, integer I4_FACTORIAL, the factorial of N.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer i4_factorial
  integer n

  i4_factorial = 1

  do i = 1, n
    i4_factorial = i4_factorial * i
  end do

  return
end
function i4_modp ( i, j )

!*****************************************************************************80
!
!! I4_MODP returns the nonnegative remainder of integer division.
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
!  Example:
!
!        I     J     MOD  I4_MODP    Factorization
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
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, the number to be divided.
!
!    Input, integer J, the number that divides I.
!
!    Output, integer I4_MODP, the nonnegative remainder when I is
!    divided by J.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer i4_modp
  integer j

  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MODP - Fatal error!'
    write ( *, '(a,i8)' ) '  I4_MODP ( I, J ) called with J = ', j
    stop 1
  end if

  i4_modp = mod ( i, j )

  if ( i4_modp < 0 ) then
    i4_modp = i4_modp + abs ( j )
  end if

  return
end
function i4_uniform_ab ( a, b )

!*****************************************************************************80
!
!! i4_uniform_ab() returns a scaled pseudorandom I4 between A and B.
!
!  Discussion:
!
!    An I4 is an integer value.
!
!    The pseudorandom number will be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer A, B, the limits of the interval.
!
!    Output, integer I4_UNIFORM_AB, a number between A and B.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer a
  integer b
  integer i4_uniform_ab
  real r
  integer value

  call random_number ( harvest = r )
!
!  Scale R to lie between A-0.5 and B+0.5.
!
  r = ( 1.0E+00 - r ) * ( real ( min ( a, b ) ) - 0.5E+00 ) & 
    +             r   * ( real ( max ( a, b ) ) + 0.5E+00 )
!
!  Use rounding to convert R to an integer between A and B.
!
  value = nint ( r )

  value = max ( value, min ( a, b ) )
  value = min ( value, max ( a, b ) )

  i4_uniform_ab = value

  return
end
subroutine i4vec_transpose_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_TRANSPOSE_PRINT prints an I4VEC "transposed".
!
!  Example:
!
!    A = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 /)
!    TITLE = 'My vector:  '
!
!    My vector:      1    2    3    4    5
!                    6    7    8    9   10
!                   11
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, integer A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer a(n)
  integer ihi
  integer ilo
  character ( len = 11 ) string
  character ( len = * ) title
  integer title_len

  title_len = len ( title )

  write ( string, '(a,i3,a)' ) '(', title_len, 'x,5i12)'

  do ilo = 1, n, 5
    ihi = min ( ilo + 5 - 1, n )
    if ( ilo == 1 ) then
      write ( *, '(a, 5i12)' ) title, a(ilo:ihi)
    else
      write ( *, string      )        a(ilo:ihi)
    end if
  end do

  return
end
subroutine ksub_random2 ( n, k, a )

!*****************************************************************************80
!
!! KSUB_RANDOM2 selects a random subset of size K from a set of size N.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 April 2003
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    This version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer N, the size of the set.
!
!    Input, integer K, the size of the subset, between 0 and N.
!
!    Output, integer A(K), the indices of the selected elements.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer k

  integer a(k)
  integer available
  integer candidate
  integer have
  integer n
  integer need
  real ( kind = rk ) r

  if ( k < 0 .or. n < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUB_RANDOM2 - Fatal error!'
    write ( *, '(a,i8)' ) '  N = ', n
    write ( *, '(a,i8)' ) '  K = ', k
    write ( *, '(a)' ) '  but 0 <= K <= N is required!'
    stop 1
  end if

  if ( k == 0 ) then
    return
  end if

  need = k
  have = 0

  available = n
  candidate = 0

  do

    candidate = candidate + 1

    call random_number ( harvest = r )

    if ( real ( available, kind = rk ) * r <= real ( need, kind = rk ) ) then

      need = need - 1
      have = have + 1
      a(have) = candidate

      if ( need <= 0 ) then
        exit
      end if

    end if

    available = available - 1

  end do

  return
end
subroutine normal ( dim_num, n, r, mu, x )

!*****************************************************************************80
!
!! NORMAL creates normally distributed points.
!
!  Discussion:
!
!    The multivariate normal distribution for the DIM_NUM dimensional vector X
!    has the form:
!
!      pdf(X) = (2*pi*det(V))^(-DIM_NUM/2)
!        * exp(-0.5*(X-MU)'*inverse(V)*(X-MU))
!
!    where MU is the mean vector, and V is a positive definite symmetric
!    matrix called the variance-covariance matrix.
!
!    This routine requires that the user supply the upper triangular
!    Cholesky factor R, which has the property that
!
!      V = R' * R
!
!    This factorization always exists if V is actually symmetric and
!    positive definite.  This factorization can be computed by the
!    routine R8PO_FA.
!
!    The user also supplies the mean vector MU.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the dimension of the space.
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) R(DIM_NUM,DIM_NUM), the upper triangular
!    Cholesky factor of the variance-covariance matrix.
!
!    Input, real ( kind = rk ) MU(DIM_NUM), the mean vector.
!
!    Output, real ( kind = rk ) X(DIM_NUM,N), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num
  integer n

  integer j
  real ( kind = rk ) mu(dim_num)
  real ( kind = rk ) r(dim_num,dim_num)
  real ( kind = rk ) x(dim_num,n)
!
!  Get a matrix of normal data.
!
  call r8vec_normal_01 ( dim_num*n, x(1:dim_num,1:n) )
!
!  Compute R' * X.
!  We actually carry out this computation in the equivalent form X' * R.
!
  do j = 1, n
    x(1:dim_num,j) = mu(1:dim_num) &
      + matmul ( x(1:dim_num,j), r(1:dim_num,1:dim_num) )
  end do

  return
end
subroutine normal_circular ( dim_num, n, x )

!*****************************************************************************80
!
!! NORMAL_CIRCULAR creates circularly normal points.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    US Department of Commerce, 1964, page 936.
!
!  Parameters:
!
!    Input, integer DIM_NUM, the dimension of the space,
!    which must be 2.
!
!    Input, integer N, the number of points.
!
!    Output, real ( kind = rk ) X(DIM_NUM,N), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num
  integer n

  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) r(n)
  real ( kind = rk ) t(n)
  real ( kind = rk ) x(dim_num,n)
!
!  The angle varies uniformly from 0 to 2 pi.
!
  call random_number ( harvest = t(1:n) )

  t(1:n) = 2.0D+00 * pi * t(1:n)
!
!  The radius is normally distributed.
!
  call r8vec_normal_01 ( n, r )

  x(1,1:n) = r(1:n) * cos ( t(1:n) )
  x(2,1:n) = r(1:n) * sin ( t(1:n) )

  return
end
subroutine normal_multivariate ( m, n, r, mu, x )

!*****************************************************************************80
!
!! NORMAL_MULTIVARIATE samples a multivariate normal distribution.
!
!  Discussion:
!
!    The multivariate normal distribution for the DIM_NUM dimensional vector X
!    has the form:
!
!      pdf(X) = (2*pi*det(V))^(-M/2) * exp(-0.5*(X-MU)'*inverse(V)*(X-MU))
!
!    where MU is the mean vector, and V is a positive definite symmetric
!    matrix called the variance-covariance matrix.
!
!    This routine samples points associated with the M-dimensional
!    normal distribution with mean MU and covariance matrix V.
!
!    This routine requires that the user supply the upper triangular
!    Cholesky factor R of V, which has the property that
!
!      V = R' * R
!
!    This factorization always exists if V is actually symmetric and
!    positive definite.  This factorization can be computed by the
!    routine R8PO_FA.
!
!    The user also supplies the mean vector MU.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 August 2004
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
!    Wiley, 1998, pages 167-168.
!
!  Parameters:
!
!    Input, integer M, the dimension of the space.
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) R(M,M), the upper triangular
!    Cholesky factor of the variance-covariance matrix.
!
!    Input, real ( kind = rk ) MU(M), the mean vector.
!
!    Output, real ( kind = rk ) X(M,N), corresponding points associated
!    with the multivariate normal distribution.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer j
  real ( kind = rk ) mu(m)
  real ( kind = rk ) r(m,m)
  real ( kind = rk ) u(m,n)
  real ( kind = rk ) x(m,n)
!
!  Create an M by N array U of samples of the standard normal distribution.
!
  call r8vec_normal_01 ( m * n, u )
!
!  Compute X = MU + R' * U.
!  We actually carry out this computation in the equivalent form MU + U' * R.
!
  do j = 1, n
    x(1:m,j) = mu(1:m) + matmul ( u(1:m,j), r(1:m,1:m) )
  end do

  return
end
subroutine normal_simple ( dim_num, n, x )

!*****************************************************************************80
!
!! NORMAL_SIMPLE creates normally distributed points.
!
!  Discussion:
!
!    The multivariate normal distribution has the form:
!
!      f(x) = (2*pi*det(V))**(-DIM_NUM/2) * exp(-0.5*(x-mu)'*inverse(V)*(x-mu))
!
!    where mu is the mean vector, and V is a positive definite symmetric
!    matrix called the variance-covariance matrix.
!
!    This routine implements the simplest version of a multivariate
!    normal distribution.  The variance-covariance matrix is the identity,
!    and the mean vector is entirely zero.  Thus, a sample on N points
!    is simply DIM_NUM*N scalar values generated under the univariate
!    normal distribution with zero mean and unit variance.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the dimension of the space.
!
!    Input, integer N, the number of points.
!
!    Output, real ( kind = rk ) X(DIM_NUM,N), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num
  integer n

  real ( kind = rk ) x(dim_num,n)

  call r8vec_normal_01 ( dim_num*n, x(1:dim_num,1:n) )

  return
end
subroutine polygon_centroid_2d ( n, v, centroid )

!*****************************************************************************80
!
!! polygon_centroid_2d() computes the centroid of a polygon in 2D.
!
!  Formula:
!
!    Denoting the centroid coordinates by CENTROID, then
!
!      CENTROID(1) = Integral ( Polygon interior ) x dx dy / Area ( Polygon )
!      CENTROID(2) = Integral ( Polygon interior ) y dx dy / Area ( Polygon ).
!
!    Green's theorem states that
!
!      Integral ( Polygon boundary ) ( M dx + N dy ) =
!      Integral ( Polygon interior ) ( dN/dx - dM/dy ) dx dy.
!
!    Using M = 0 and N = x * x / 2, we get:
!
!      CENTROID(1) = 0.5 * Integral ( Polygon boundary ) x * x dy,
!
!    which becomes
!
!      CENTROID(1) = 1/6 Sum ( 1 <= I <= N )
!        ( X(I+1) + X(I) ) * ( X(I) * Y(I+1) - X(I+1) * Y(I))
!
!    where, when I = N, the index "I+1" is replaced by 1.
!
!    A similar calculation gives us a formula for CENTROID(2).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Gerard Bashein, Paul Detmer,
!    Centroid of a Polygon,
!    Graphics Gems IV, edited by Paul Heckbert,
!    AP Professional, 1994.
!
!  Parameters:
!
!    Input, integer N, the number of sides of the polygonal shape.
!
!    Input, real ( kind = rk ) V(2,N), the coordinates of the vertices
!    of the shape.
!
!    Output, real ( kind = rk ) CENTROID(2), the coordinates of the
!    centroid of the shape.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) area
  real ( kind = rk ) centroid(2)
  integer i
  integer ip1
  real ( kind = rk ) temp
  real ( kind = rk ) v(2,n)

  area = 0.0D+00
  centroid(1:2) = 0.0D+00

  do i = 1, n

    if ( i < n ) then
      ip1 = i + 1
    else
      ip1 = 1
    end if

    temp = ( v(1,i) * v(2,ip1) - v(1,ip1) * v(2,i) )

    area = area + temp

    centroid(1:2) = centroid(1:2) + ( v(1:2,ip1) + v(1:2,i) ) * temp

  end do

  area = area / 2.0D+00

  centroid(1:2) = centroid(1:2) / ( 6.0D+00 * area )

  return
end
function prime ( n )

!*****************************************************************************80
!
!! PRIME returns any of the first PRIME_MAX prime numbers.
!
!  Discussion:
!
!    PRIME_MAX is 1600, and the largest prime stored is 13499.
!
!    Thanks to Bart Vandewoestyne for pointing out a typo, 18 February 2005.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    US Department of Commerce, 1964, pages 870-873.
!
!    Daniel Zwillinger,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996, pages 95-98.
!
!  Parameters:
!
!    Input, integer N, the index of the desired prime number.
!    In general, is should be true that 0 <= N <= PRIME_MAX.
!    N = -1 returns PRIME_MAX, the index of the largest prime available.
!    N = 0 is legal, returning PRIME = 1.
!
!    Output, integer PRIME, the N-th prime.  If N is out of range,
!    PRIME is returned as -1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: prime_max = 1600

  integer, save :: icall = 0
  integer n
  integer, save, dimension ( prime_max ) :: npvec
  integer prime

  if ( icall == 0 ) then

    icall = 1

    npvec(1:100) = (/ &
        2,    3,    5,    7,   11,   13,   17,   19,   23,   29, &
       31,   37,   41,   43,   47,   53,   59,   61,   67,   71, &
       73,   79,   83,   89,   97,  101,  103,  107,  109,  113, &
      127,  131,  137,  139,  149,  151,  157,  163,  167,  173, &
      179,  181,  191,  193,  197,  199,  211,  223,  227,  229, &
      233,  239,  241,  251,  257,  263,  269,  271,  277,  281, &
      283,  293,  307,  311,  313,  317,  331,  337,  347,  349, &
      353,  359,  367,  373,  379,  383,  389,  397,  401,  409, &
      419,  421,  431,  433,  439,  443,  449,  457,  461,  463, &
      467,  479,  487,  491,  499,  503,  509,  521,  523,  541 /)

    npvec(101:200) = (/ &
      547,  557,  563,  569,  571,  577,  587,  593,  599,  601, &
      607,  613,  617,  619,  631,  641,  643,  647,  653,  659, &
      661,  673,  677,  683,  691,  701,  709,  719,  727,  733, &
      739,  743,  751,  757,  761,  769,  773,  787,  797,  809, &
      811,  821,  823,  827,  829,  839,  853,  857,  859,  863, &
      877,  881,  883,  887,  907,  911,  919,  929,  937,  941, &
      947,  953,  967,  971,  977,  983,  991,  997, 1009, 1013, &
     1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, &
     1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, &
     1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223 /)

    npvec(201:300) = (/ &
     1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, &
     1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, &
     1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, &
     1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, &
     1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, &
     1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, &
     1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, &
     1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, &
     1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, &
     1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987 /)

    npvec(301:400) = (/ &
     1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053, &
     2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129, &
     2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, &
     2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287, &
     2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, &
     2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, &
     2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, &
     2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, &
     2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687, &
     2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741 /)

    npvec(401:500) = (/ &
     2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, &
     2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903, &
     2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, &
     3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, &
     3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181, &
     3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257, &
     3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, &
     3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413, &
     3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, &
     3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571 /)

    npvec(501:600) = (/ &
     3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, &
     3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727, &
     3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821, &
     3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907, &
     3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989, &
     4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057, &
     4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139, &
     4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231, &
     4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297, &
     4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409 /)

    npvec(601:700) = (/ &
     4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493, &
     4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583, &
     4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657, &
     4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751, &
     4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, &
     4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937, &
     4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003, &
     5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087, &
     5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179, &
     5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279 /)

    npvec(701:800) = (/ &
     5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387, &
     5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443, &
     5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521, &
     5527, 5531, 5557, 5563, 5569, 5573, 5581, 5591, 5623, 5639, &
     5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693, &
     5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791, &
     5801, 5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857, &
     5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939, &
     5953, 5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053, &
     6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133 /)

    npvec(801:900) = (/ &
     6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221, &
     6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301, &
     6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367, &
     6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473, &
     6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571, &
     6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673, &
     6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761, &
     6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833, &
     6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917, &
     6947, 6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997 /)

    npvec(901:1000) = (/ &
     7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103, &
     7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207, &
     7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297, &
     7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411, &
     7417, 7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499, &
     7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561, &
     7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643, &
     7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723, &
     7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829, &
     7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919 /)

    npvec(1001:1100) = (/ &
     7927, 7933, 7937, 7949, 7951, 7963, 7993, 8009, 8011, 8017, &
     8039, 8053, 8059, 8069, 8081, 8087, 8089, 8093, 8101, 8111, &
     8117, 8123, 8147, 8161, 8167, 8171, 8179, 8191, 8209, 8219, &
     8221, 8231, 8233, 8237, 8243, 8263, 8269, 8273, 8287, 8291, &
     8293, 8297, 8311, 8317, 8329, 8353, 8363, 8369, 8377, 8387, &
     8389, 8419, 8423, 8429, 8431, 8443, 8447, 8461, 8467, 8501, &
     8513, 8521, 8527, 8537, 8539, 8543, 8563, 8573, 8581, 8597, &
     8599, 8609, 8623, 8627, 8629, 8641, 8647, 8663, 8669, 8677, &
     8681, 8689, 8693, 8699, 8707, 8713, 8719, 8731, 8737, 8741, &
     8747, 8753, 8761, 8779, 8783, 8803, 8807, 8819, 8821, 8831 /)

    npvec(1101:1200) = (/ &
     8837, 8839, 8849, 8861, 8863, 8867, 8887, 8893, 8923, 8929, &
     8933, 8941, 8951, 8963, 8969, 8971, 8999, 9001, 9007, 9011, &
     9013, 9029, 9041, 9043, 9049, 9059, 9067, 9091, 9103, 9109, &
     9127, 9133, 9137, 9151, 9157, 9161, 9173, 9181, 9187, 9199, &
     9203, 9209, 9221, 9227, 9239, 9241, 9257, 9277, 9281, 9283, &
     9293, 9311, 9319, 9323, 9337, 9341, 9343, 9349, 9371, 9377, &
     9391, 9397, 9403, 9413, 9419, 9421, 9431, 9433, 9437, 9439, &
     9461, 9463, 9467, 9473, 9479, 9491, 9497, 9511, 9521, 9533, &
     9539, 9547, 9551, 9587, 9601, 9613, 9619, 9623, 9629, 9631, &
     9643, 9649, 9661, 9677, 9679, 9689, 9697, 9719, 9721, 9733 /)

    npvec(1201:1300) = (/ &
     9739, 9743, 9749, 9767, 9769, 9781, 9787, 9791, 9803, 9811, &
     9817, 9829, 9833, 9839, 9851, 9857, 9859, 9871, 9883, 9887, &
     9901, 9907, 9923, 9929, 9931, 9941, 9949, 9967, 9973,10007, &
    10009,10037,10039,10061,10067,10069,10079,10091,10093,10099, &
    10103,10111,10133,10139,10141,10151,10159,10163,10169,10177, &
    10181,10193,10211,10223,10243,10247,10253,10259,10267,10271, &
    10273,10289,10301,10303,10313,10321,10331,10333,10337,10343, &
    10357,10369,10391,10399,10427,10429,10433,10453,10457,10459, &
    10463,10477,10487,10499,10501,10513,10529,10531,10559,10567, &
    10589,10597,10601,10607,10613,10627,10631,10639,10651,10657 /)

    npvec(1301:1400) = (/ &
    10663,10667,10687,10691,10709,10711,10723,10729,10733,10739, &
    10753,10771,10781,10789,10799,10831,10837,10847,10853,10859, &
    10861,10867,10883,10889,10891,10903,10909,10937,10939,10949, &
    10957,10973,10979,10987,10993,11003,11027,11047,11057,11059, &
    11069,11071,11083,11087,11093,11113,11117,11119,11131,11149, &
    11159,11161,11171,11173,11177,11197,11213,11239,11243,11251, &
    11257,11261,11273,11279,11287,11299,11311,11317,11321,11329, &
    11351,11353,11369,11383,11393,11399,11411,11423,11437,11443, &
    11447,11467,11471,11483,11489,11491,11497,11503,11519,11527, &
    11549,11551,11579,11587,11593,11597,11617,11621,11633,11657 /)

    npvec(1401:1500) = (/ &
    11677,11681,11689,11699,11701,11717,11719,11731,11743,11777, &
    11779,11783,11789,11801,11807,11813,11821,11827,11831,11833, &
    11839,11863,11867,11887,11897,11903,11909,11923,11927,11933, &
    11939,11941,11953,11959,11969,11971,11981,11987,12007,12011, &
    12037,12041,12043,12049,12071,12073,12097,12101,12107,12109, &
    12113,12119,12143,12149,12157,12161,12163,12197,12203,12211, &
    12227,12239,12241,12251,12253,12263,12269,12277,12281,12289, &
    12301,12323,12329,12343,12347,12373,12377,12379,12391,12401, &
    12409,12413,12421,12433,12437,12451,12457,12473,12479,12487, &
    12491,12497,12503,12511,12517,12527,12539,12541,12547,12553 /)

   npvec(1501:1600) = (/ &
    12569,12577,12583,12589,12601,12611,12613,12619,12637,12641, &
    12647,12653,12659,12671,12689,12697,12703,12713,12721,12739, &
    12743,12757,12763,12781,12791,12799,12809,12821,12823,12829, &
    12841,12853,12889,12893,12899,12907,12911,12917,12919,12923, &
    12941,12953,12959,12967,12973,12979,12983,13001,13003,13007, &
    13009,13033,13037,13043,13049,13063,13093,13099,13103,13109, &
    13121,13127,13147,13151,13159,13163,13171,13177,13183,13187, &
    13217,13219,13229,13241,13249,13259,13267,13291,13297,13309, &
    13313,13327,13331,13337,13339,13367,13381,13397,13399,13411, &
    13417,13421,13441,13451,13457,13463,13469,13477,13487,13499 /)

  end if

  if ( n == -1 ) then
    prime = prime_max
  else if ( n == 0 ) then
    prime = 1
  else if ( n <= prime_max ) then
    prime = npvec(n)
  else
    prime = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PRIME - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal prime index N = ', n
    write ( *, '(a,i8)' ) '  N should be between 1 and PRIME_MAX =', prime_max
    stop 1
  end if

  return
end
function r8_acos ( c )

!*****************************************************************************80
!
!! R8_ACOS computes the arc cosine function, with argument truncation.
!
!  Discussion:
!
!    If you call your system ACOS routine with an input argument that is
!    even slightly outside the range [-1.0, 1.0 ], you may get an unpleasant
!    surprise (I did).
!
!    This routine simply truncates arguments outside the range.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) C, the argument.
!
!    Output, real ( kind = rk ) R8_ACOS, an angle whose cosine is C.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) c
  real ( kind = rk ) c2
  real ( kind = rk ) r8_acos

  c2 = c
  c2 = max ( c2, -1.0D+00 )
  c2 = min ( c2, +1.0D+00 )

  r8_acos = acos ( c2 )

  return
end
function r8_normal_01 ( )

!*****************************************************************************80
!
!! R8_NORMAL_01 returns a unit pseudonormal R8.
!
!  Discussion:
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!    Because this routine uses the Box Muller method, it requires pairs
!    of uniform random values to generate a pair of normal random values.
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
!  Parameters:
!
!    Output, real ( kind = rk ) R8_NORMAL_01, a sample of the standard
!    normal PDF.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) r1
  real ( kind = rk ) r2
  real ( kind = rk ) r8_normal_01
  integer, save :: used = 0
  real ( kind = rk ) x
  real ( kind = rk ), save :: y = 0.0D+00
!
!  On odd numbered calls, generate two uniforms, create two normals,
!  return the first normal.
!
  if ( mod ( used, 2 ) == 0 ) then

    call random_number ( harvest = r1 )
    call random_number ( harvest = r2 )

    x = sqrt ( -2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * pi * r2 )
    y = sqrt ( -2.0D+00 * log ( r1 ) ) * sin ( 2.0D+00 * pi * r2 )
!
!  On odd calls, return the second normal.
!
  else

    x = y

  end if

  used = used + 1

  r8_normal_01 = x

  return
end
subroutine r8mat_normal_01 ( m, n, r )

!*****************************************************************************80
!
!! r8mat_normal_01() returns a unit pseudonormal R8MAT.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, 1969, pages 136-143.
!
!  Input:
!
!    integer M, N, the number of rows and columns
!    in the array.
!
!  Output:
!
!    real ( kind = rk ) R(M,N), the array of pseudonormal values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) r(m,n)

  call r8vec_normal_01 ( m * n, r )

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints a R8MAT.
!
!  Discussion:
!
!    A R8MAT is a two dimensional matrix of double precision real values.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows in A.
!
!    Input, integer N, the number of columns in A.
!
!    Input, real ( kind = rk ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  character ( len = * ) title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_PRINT_SOME prints some of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = rk ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ILO, JLO, the first row and column to print.
!
!    Input, integer IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: incx = 5
  integer m
  integer n

  real ( kind = rk ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer i
  integer i2hi
  integer i2lo
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == real ( int ( a(i,j) ), kind = rk ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine r8mat_write ( output_file_name, m, n, table )

!*****************************************************************************80
!
!! R8MAT_WRITE writes an R8MAT file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILE_NAME, the output file name.
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) TABLE(M,N), the table data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer j
  character ( len = * ) output_file_name
  integer output_status
  integer output_unit
  character ( len = 30 ) string
  real ( kind = rk ) table(m,n)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_file_name, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_file_name ) // '" on unit ', output_unit
    output_unit = -1
    stop 1
  end if
!
!  Create a format string.
!
!  For greater precision in the output file, try:
!
!                                            '(', m, 'g', 24, '.', 16, ')'
!
  write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 14, '.', 6, ')'
!
!  Write the data.
!
  do j = 1, n
    write ( output_unit, string ) table(1:m,j)
  end do
!
!  Close the file.
!
  close ( unit = output_unit )

  return
end
subroutine r8po_fa ( n, a, info )

!*****************************************************************************80
!
!! R8PO_FA factors an R8PO matrix.
!
!  Discussion:
!
!    The R8PO storage format is used for a symmetric positive definite
!    matrix and its inverse.  (The Cholesky factor of a R8PO matrix is an
!    upper triangular matrix, so it will be in DGE storage format.)
!
!    Only the diagonal and upper triangle of the square array are used.
!    This same storage scheme is used when the matrix is factored by
!    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
!    is set to zero.
!
!    R8PO storage is used by LINPACK and LAPACK.
!
!    The positive definite symmetric matrix A has a Cholesky factorization
!    of the form:
!
!      A = R' * R
!
!    where R is an upper triangular matrix with positive elements on
!    its diagonal.  This routine overwrites the matrix A with its
!    factor R.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 March 2003
!
!  Author:
!
!    Original FORTRAN77 version by Jack Dongarra, Jim Bunch,
!    Cleve Moler, Pete Stewart.
!    This version by John Burkardt.
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real ( kind = rk ) A(N,N).
!    On input, the matrix in R8PO storage.
!    On output, the Cholesky factor R in DGE storage.
!
!    Output, integer INFO, error flag.
!    0, normal return.
!    K, error condition.  The principal minor of order K is not
!    positive definite, and the factorization was not completed.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,n)
  integer i
  integer info
  integer j
  integer k
  real ( kind = rk ) s

  do j = 1, n

    do k = 1, j - 1
      a(k,j) = ( a(k,j) - sum ( a(1:k-1,k) * a(1:k-1,j) ) ) / a(k,k)
    end do

    s = a(j,j) - sum ( a(1:j-1,j)**2 )

    if ( s <= 0.0D+00 ) then
      info = j
      return
    end if

    a(j,j) = sqrt ( s )

  end do

  info = 0
!
!  Since the Cholesky factor is stored in DGE format, be sure to
!  zero out the lower triangle.
!
  do i = 1, n
    do j = 1, i - 1
      a(i,j) = 0.0D+00
    end do
  end do

  return
end
subroutine r8po_sl ( n, a_lu, b )

!*****************************************************************************80
!
!! R8PO_SL solves an R8PO system factored by R8PO_FA.
!
!  Discussion:
!
!    The R8PO storage format is used for a symmetric positive definite
!    matrix and its inverse.  (The Cholesky factor of a R8PO matrix is an
!    upper triangular matrix, so it will be in DGE storage format.)
!
!    Only the diagonal and upper triangle of the square array are used.
!    This same storage scheme is used when the matrix is factored by
!    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
!    is set to zero.
!
!    R8PO storage is used by LINPACK and LAPACK.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 March 1999
!
!  Author:
!
!    Original FORTRAN77 version by Jack Dongarra, Jim Bunch,
!    Cleve Moler, Pete Stewart.
!    This version by John Burkardt.
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real ( kind = rk ) A_LU(N,N), the Cholesky factor from R8PO_FA.
!
!    Input/output, real ( kind = rk ) B(N).
!    On input, the right hand side.
!    On output, the solution vector.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a_lu(n,n)
  real ( kind = rk ) b(n)
  integer k
!
!  Solve R' * y = b.
!
  do k = 1, n
    b(k) = ( b(k) - sum ( b(1:k-1) * a_lu(1:k-1,k) ) ) / a_lu(k,k)
  end do
!
!  Solve R * x = y.
!
  do k = n, 1, -1
    b(k) = b(k) / a_lu(k,k)
    b(1:k-1) = b(1:k-1) - a_lu(1:k-1,k) * b(k)
  end do

  return
end
function r8vec_norm ( n, a )

!*****************************************************************************80
!
!! R8VEC_NORM returns the L2 norm of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The vector L2 norm is defined as:
!
!      R8VEC_NORM = sqrt ( sum ( 1 <= I <= N ) A(I)**2 ).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in A.
!
!    Input, real ( kind = rk ) A(N), the vector whose L2 norm is desired.
!
!    Output, real ( kind = rk ) R8VEC_NORM, the L2 norm of A.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) r8vec_norm

  r8vec_norm = sqrt ( sum ( a(1:n)**2 ) )

  return
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
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! r8vec_print() prints a R8VEC.
!
!  Discussion:
!
!    A R8VEC is an array of double precision real values.
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
!    Input, character ( len = * ) TITLE, an optional title.
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
    write ( *, '(2x,i8,2x,g16.8)' ) i, a(i)
  end do

  return
end
subroutine scale_from_simplex01 ( dim_num, n, t, x )

!*****************************************************************************80
!
!! SCALE_FROM_SIMPLEX01 rescales data from a unit to non-unit simplex.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 August 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Reuven Rubinstein,
!    Monte Carlo Optimization, Simulation, and Sensitivity
!    of Queueing Networks,
!    Krieger, 1992,
!    ISBN: 0894647644,
!    LC: QA298.R79.
!
!  Parameters:
!
!    Input, integer DIM_NUM, the dimension of the space.
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) T(DIM_NUM,0:DIM_NUM), the coordinates of the
!    DIM_NUM+1 points that define the simplex.  T(1:DIM_NUM,0) corresponds
!    to the origin, and T(1:DIM_NUM,J) will be the image of the J-th unit
!    coordinate vector.
!
!    Input/output, real ( kind = rk ) X(DIM_NUM,N), the data to be modified.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num
  integer n

  real ( kind = rk ) a(dim_num,dim_num)
  integer j
  real ( kind = rk ) t(dim_num,0:dim_num)
  real ( kind = rk ) x(dim_num,n)

  a(1:dim_num,1:dim_num) = t(1:dim_num,1:dim_num)

  do j = 1, dim_num
    a(1:dim_num,j) = a(1:dim_num,j) - t(1:dim_num,0)
  end do

  do j = 1, n
    x(1:dim_num,j) = matmul ( a(1:dim_num,1:dim_num), x(1:dim_num,j) ) &
      + t(1:dim_num,0)
  end do

  return
end
subroutine scale_to_ball01 ( dim_num, n, x )

!*****************************************************************************80
!
!! SCALE_TO_BALL01 translates and rescales data to fit within the unit ball.
!
!  Discussion:
!
!    Completely arbitrary input data is given.
!
!    The average of the data is computed, and taken as the coordinates
!    of the center C of a sphere.  The radius R of that sphere is the
!    distance from the center to the furthest point in the data set.
!
!    Then each point is transformed to the ball of center 0 and radius
!    1 by subtracting C and dividing by R:
!
!      X(1:DIM_NUM,J) -> ( X(1:DIM_NUM,J) - C(1:DIM_NUM) ) / R
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the dimension of the space.
!
!    Input, integer N, the number of points.
!
!    Input/output, real ( kind = rk ) X(DIM_NUM,N), the data to be modified.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num
  integer n

  integer i
  integer j
  real ( kind = rk ) r(n)
  real ( kind = rk ) scale
  real ( kind = rk ) x(dim_num,n)
  real ( kind = rk ) xave(dim_num)
!
!  Determine the center.
!
  do i = 1, dim_num
    xave(i) = sum ( x(i,1:n) ) / real ( n, kind = rk )
  end do
!
!  Determine SCALE, the maximum distance of any point X from the center.
!
  do j = 1, n
    r(j) = sum ( ( x(1:dim_num,j) - xave(1:dim_num) )**2 )
  end do

  scale = sqrt ( maxval ( r(1:n) ) )
!
!  Dividing all values by SCALE will guarantee that every point is
!  inside the unit sphere, and one point at least is ON the sphere.
!
  if ( 0.0D+00 < scale ) then
    do i = 1, dim_num
      x(i,1:n) = ( x(i,1:n) - xave(i) ) / scale
    end do
  else
    x(1:dim_num,1:n) = 0.0D+00
  end if

  return
end
subroutine scale_to_block01 ( dim_num, n, x )

!*****************************************************************************80
!
!! SCALE_TO_BLOCK01 translates and rescales data to fit in the unit block.
!
!  Discussion:
!
!    The minimum and maximum coordinate values M1(I) and M2(I) are
!    determined, and the maximum of M2(I) - M1(I) is used to scale
!    all the coordinates by the same factor.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the dimension of the space.
!
!    Input, integer N, the number of points.
!
!    Input/output, real ( kind = rk ) X(DIM_NUM,N), the data to be modified.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num
  integer n

  integer i
  real ( kind = rk ) x(dim_num,n)
  real ( kind = rk ) xmax(dim_num)
  real ( kind = rk ) xmin(dim_num)
  real ( kind = rk ) xrange
  real ( kind = rk ) xrange2
!
!  Determine the extremes in each dimension.
!
  xrange = 0.0D+00
  do i = 1, dim_num
    xmin(i) = minval ( x(i,1:n) )
    xmax(i) = maxval ( x(i,1:n) )
    xrange = max ( xrange, xmax(i) - xmin(i) )
  end do
!
!  Extend all the extremes so that the range is the same in each dimension.
!
  do i = 1, dim_num
    xrange2 = xrange - ( xmax(i) - xmin(i) )
    xmax(i) = xmax(i) + 0.5D+00 * xrange2
    xmin(i) = xmin(i) - 0.5D+00 * xrange2
  end do
!
!  Now map the data to [0,1], using a single dilation factor for all dimensions.
!
  if ( 0.0D+00 == xrange ) then

    x(1:dim_num,1:n) = 0.5D+00

  else

    do i = 1, dim_num
      x(i,1:n) = ( x(i,1:n) - xmin(i) ) / xrange
    end do

  end if

  return
end
subroutine scale_to_cube01 ( dim_num, n, x )

!*****************************************************************************80
!
!! SCALE_TO_CUBE01 translates and rescales data to the unit hypercube.
!
!  Discussion:
!
!    In each coordinate dimension I, the minimum and maximum coordinate
!    values M1(I) and M2(I) are determined.
!
!    Then, in each coordinate, the points are rescaled as
!
!      X(I) -> ( X(I) - M1(I) ) / ( M2(I) - M1(I) ).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the dimension of the space.
!
!    Input, integer N, the number of points.
!
!    Input/output, real ( kind = rk ) X(DIM_NUM,N), the data to be modified.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num
  integer n

  integer i
  real ( kind = rk ) x(dim_num,n)
  real ( kind = rk ) xmax(dim_num)
  real ( kind = rk ) xmin(dim_num)

  do i = 1, dim_num
    xmin(i) = minval ( x(i,1:n) )
    xmax(i) = maxval ( x(i,1:n) )
  end do

  do i = 1, dim_num
    if ( 0.0D+00 < xmax(i) - xmin(i) ) then
      x(i,1:n) = ( x(i,1:n) - xmin(i) ) / ( xmax(i) - xmin(i) )
    else
      x(i,1:n) = 0.5D+00
    end if
  end do

  return
end
subroutine stri_angles_to_area ( r, a, b, c, area )

!*****************************************************************************80
!
!! STRI_ANGLES_TO_AREA computes the area of a spherical triangle.
!
!  Discussion:
!
!    A sphere centered at 0 in 3D satisfies the equation:
!
!      X*X + Y*Y + Z*Z = R*R
!
!    A spherical triangle is specified by three points on the surface
!    of the sphere.
!
!    The area formula is known as Girard's formula.
!
!    The area of a spherical triangle is:
!
!      AREA = ( A + B + C - PI ) * R*R
!
!    where A, B and C are the (surface) angles of the triangle.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) R, the radius of the sphere.
!
!    Input, real ( kind = rk ) A, B, C, the angles of the triangle.
!
!    Output, real ( kind = rk ) AREA, the area of the spherical triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) area
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) r
!
!  Apply Girard's formula.
!
  area = r * r * ( a + b + c - pi )

  return
end
subroutine stri_sides_to_angles ( r, as, bs, cs, a, b, c )

!*****************************************************************************80
!
!! STRI_SIDES_TO_ANGLES computes spherical triangle angles.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) R, the radius of the sphere.
!
!    Input, real ( kind = rk ) AS, BS, CS, the (geodesic) length of the
!    sides of the triangle.
!
!    Output, real ( kind = rk ) A, B, C, the spherical angles of the triangle.
!    Angle A is opposite the side of length AS, and so on.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) as
  real ( kind = rk ) asu
  real ( kind = rk ) b
  real ( kind = rk ) bs
  real ( kind = rk ) bsu
  real ( kind = rk ) c
  real ( kind = rk ) cs
  real ( kind = rk ) csu
  real ( kind = rk ) r
  real ( kind = rk ) ssu
  real ( kind = rk ) tan_a2
  real ( kind = rk ) tan_b2
  real ( kind = rk ) tan_c2

  asu = as / r
  bsu = bs / r
  csu = cs / r
  ssu = ( asu + bsu + csu ) / 2.0D+00

  tan_a2 = sqrt ( ( sin ( ssu - bsu ) * sin ( ssu - csu ) ) / &
                  ( sin ( ssu ) * sin ( ssu - asu )     ) )

  a = 2.0D+00 * atan ( tan_a2 )

  tan_b2 = sqrt ( ( sin ( ssu - asu ) * sin ( ssu - csu ) ) / &
                  ( sin ( ssu ) * sin ( ssu - bsu )     ) )

  b = 2.0D+00 * atan ( tan_b2 )

  tan_c2 = sqrt ( ( sin ( ssu - asu ) * sin ( ssu - bsu ) ) / &
                  ( sin ( ssu ) * sin ( ssu - csu )     ) )

  c = 2.0D+00 * atan ( tan_c2 )

  return
end
subroutine stri_vertices_to_sides ( r, v1, v2, v3, as, bs, cs )

!*****************************************************************************80
!
!! STRI_VERTICES_TO_SIDES_3D computes spherical triangle sides.
!
!  Discussion:
!
!    We can use the ACOS system call here, but the R8_ACOS routine
!    will automatically take care of cases where the input argument is
!    (usually slightly) out of bounds.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 August 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) R, the radius of the sphere.
!
!    Input, real ( kind = rk ) V1(3), V2(3), V3(3), the vertices of the spherical
!    triangle.
!
!    Output, real ( kind = rk ) AS, BS, CS, the (geodesic) length of the sides
!    of the triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) as
  real ( kind = rk ) bs
  real ( kind = rk ) cs
  real ( kind = rk ) r
  real ( kind = rk ) r8_acos
  real ( kind = rk ) v1(3)
  real ( kind = rk ) v2(3)
  real ( kind = rk ) v3(3)

  as = r * r8_acos ( dot_product ( v2(1:3), v3(1:3) ) / r**2 )
  bs = r * r8_acos ( dot_product ( v3(1:3), v1(1:3) ) / r**2 )
  cs = r * r8_acos ( dot_product ( v1(1:3), v2(1:3) ) / r**2 )

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
!  Parameters:
!
!    None
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
subroutine triangle_area_2d ( t, area )

!*****************************************************************************80
!
!! TRIANGLE_AREA_2D computes the area of a triangle in 2D.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) T(2,3), the triangle vertices.
!
!    Output, real ( kind = rk ) AREA, the absolute area of the triangle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2

  real ( kind = rk ) area
  real ( kind = rk ) t(dim_num,3)

  area = 0.5D+00 * abs ( &
      t(1,1) * ( t(2,2) - t(2,3) ) &
    + t(1,2) * ( t(2,3) - t(2,1) ) &
    + t(1,3) * ( t(2,1) - t(2,2) ) )

  return
end
subroutine tuple_next_fast ( m, n, rank, x )

!*****************************************************************************80
!
!! TUPLE_NEXT_FAST computes the next element of a tuple space, "fast".
!
!  Discussion:
!
!    The elements are N vectors.  Each entry is constrained to lie
!    between 1 and M.  The elements are produced one at a time.
!    The first element is
!      (1,1,...,1)
!    and the last element is
!      (M,M,...,M)
!    Intermediate elements are produced in lexicographic order.
!
!    This code was written as a possibly faster version of TUPLE_NEXT.
!
!  Example:
!
!    N = 2,
!    M = 3
!
!    INPUT        OUTPUT
!    -------      -------
!    Rank          X
!    ----          ----
!   -1            -1 -1
!
!    0             1  1
!    1             1  2
!    2             1  3
!    3             2  1
!    4             2  2
!    5             2  3
!    6             3  1
!    7             3  2
!    8             3  3
!    9             1  1
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
!  Parameters:
!
!    Input, integer M, the maximum entry in any component.
!    M must be greater than 0.
!
!    Input, integer N, the number of components.
!    N must be greater than 0.
!
!    Input, integer RANK, indicates the rank of the tuple.
!    Typically, 0 <= RANK < N**M.  Values of RANK greater than
!    N**M are legal and meaningful; they are equivalent to the
!    corresponding value mod (N**M).  If RANK < 0, this indicates
!    that this is the first call for the given values of (M,N).
!    Initialization is done, and X is set to a dummy value.
!
!    Output, integer X(N), the next tuple, or a dummy value if
!    initialization has just been done.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer, save, allocatable, dimension ( : ) :: base
  integer i
  integer m
  integer rank
  integer x(n)

  if ( rank < 0 ) then

    if ( m <= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TUPLE_NEXT_FAST - Fatal error!'
      write ( *, '(a)' ) '  The value M <= 0 is not allowed.'
      write ( *, '(a,i8)' ) '  M = ', m
      stop 1
    end if

    if ( n <= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TUPLE_NEXT_FAST - Fatal error!'
      write ( *, '(a)' ) '  The value N <= 0 is not allowed.'
      write ( *, '(a,i8)' ) '  N = ', n
      stop 1
    end if

    if ( allocated ( base ) ) then
      deallocate ( base )
    end if
    allocate ( base(1:n) )

    base(n) = 1
    do i = n-1, 1, -1
      base(i) = base(i+1) * m
    end do

    x(1:n) = -1

  else

    x(1:n) = mod ( rank / base(1:n), m ) + 1

  end if

  return
end
subroutine uniform_in_annulus ( pc, r1, r2, n, p )

!*****************************************************************************80
!
!! UNIFORM_IN_ANNULUS samples a circular annulus.
!
!  Discussion:
!
!    A circular annulus with center PC, inner radius R1 and
!    outer radius R2, is the set of points P so that
!
!      R1^2 <= (P(1)-PC(1))^2 + (P(2)-PC(2))^2 <= R2^2
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Peter Shirley,
!    Nonuniform Random Point Sets Via Warping,
!    Graphics Gems, Volume III,
!    edited by David Kirk,
!    AP Professional, 1992,
!    ISBN: 0122861663,
!    LC: T385.G6973.
!
!  Parameters:
!
!    Input, real ( kind = rk ) PC(2), the center.
!
!    Input, real ( kind = rk ) R1, R2, the inner and outer radii.
!
!    Input, integer N, the number of points to generate.
!
!    Output, real ( kind = rk ) P(2,N), sample points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2
  integer n

  real ( kind = rk ) p(dim_num,n)
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) r(n)
  real ( kind = rk ) r1
  real ( kind = rk ) r2
  real ( kind = rk ) theta(n)
  real ( kind = rk ) u(n)
  real ( kind = rk ) v(n)

  call random_number ( harvest = u(1:n) )

  theta(1:n) = u(1:n) * 2.0D+00 * r8_pi

  call random_number ( harvest = v(1:n) )

  r(1:n) = sqrt ( ( 1.0D+00 - v(1:n) ) * r1 * r1 &
           +                  v(1:n)   * r2 * r2 )

  p(1,1:n) = pc(1) + r(1:n) * cos ( theta(1:n) )
  p(2,1:n) = pc(2) + r(1:n) * sin ( theta(1:n) )

  return
end
subroutine uniform_in_annulus_accept ( pc, r1, r2, n, x )

!*****************************************************************************80
!
!! UNIFORM_IN_ANNULUS_ACCEPT accepts points in an annulus.
!
!  Discussion:
!
!    A circular annulus with center PC, inner radius R1 and
!    outer radius R2, is the set of points P so that
!
!      R1**2 <= (P(1)-PC(1))**2 + (P(2)-PC(2))**2 <= R2**2
!
!    The acceptance/rejection method is used.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) PC(2), the center.
!
!    Input, real ( kind = rk ) R1, R2, the inner and outer radii.
!
!    Input, integer N, the number of points.
!
!    Output, real ( kind = rk ) X(DIM_NUM,N), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer, parameter :: dim_num = 2

  integer j
  real ( kind = rk ) r1
  real ( kind = rk ) r2
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) u(dim_num)
  real ( kind = rk ) x(dim_num,n)

  if ( r2 <= r1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'UNIFORM_IN_ANNULUS_ACCEPT - Fatal error!'
    write ( *, '(a)' ) '  R2 <= R1.'
    return
  end if
!
!  Generate points in a square of "radius" R2.
!  Accept those points which lie inside the circle of radius R2, and outside
!  the circle of radius R1.
!
  do j = 1, n

    do

      call random_number ( harvest = u(1:dim_num) )

      u(1:dim_num) = ( 2.0D+00 * u(1:dim_num) - 1.0D+00 ) * r2

      if ( r1 * r1 <= sum ( u(1:dim_num)**2 ) .and. &
        sum ( u(1:dim_num)**2 ) <= r2 * r2 ) then
        exit
      end if

    end do

    x(1:dim_num,j) = pc(1:dim_num) + u(1:dim_num)

  end do

  return
end
subroutine uniform_in_annulus_sector ( pc, r1, r2, theta1, theta2, n, p )

!*****************************************************************************80
!
!! UNIFORM_IN_ANNULUS_SECTOR samples an annular sector in 2D.
!
!  Discussion:
!
!    An annular sector with center PC, inner radius R1 and
!    outer radius R2, and angles THETA1, THETA2, is the set of points
!    P so that
!
!      R1**2 <= (P(1)-PC(1))**2 + (P(2)-PC(2))**2 <= R2**2
!
!    and
!
!      THETA1 <= THETA ( P - PC ) <= THETA2
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Peter Shirley,
!    Nonuniform Random Point Sets Via Warping,
!    Graphics Gems, Volume III,
!    edited by David Kirk,
!    AP Professional, 1992,
!    ISBN: 0122861663,
!    LC: T385.G6973.
!
!  Parameters:
!
!    Input, real ( kind = rk ) PC(2), the center.
!
!    Input, real ( kind = rk ) R1, R2, the inner and outer radii.
!
!    Input, real ( kind = rk ) THETA1, THETA2, the angles.
!
!    Input, integer N, the number of points to generate.
!
!    Output, real ( kind = rk ) P(2,N), sample points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2
  integer n

  real ( kind = rk ) p(2,n)
  real ( kind = rk ) pc(dim_num)
  real ( kind = rk ) r(1:n)
  real ( kind = rk ) r1
  real ( kind = rk ) r2
  real ( kind = rk ) theta(1:n)
  real ( kind = rk ) theta1
  real ( kind = rk ) theta2
  real ( kind = rk ) u(1:n)
  real ( kind = rk ) v(1:n)

  call random_number ( harvest = u(1:n) )

  theta(1:n) = ( 1.0D+00 - u(1:n) ) * theta1 &
              +            u(1:n)   * theta2

  call random_number ( harvest = v(1:n) )

  r(1:n) = sqrt ( ( 1.0D+00 - v(1:n) ) * r1 * r1 &
           +                  v(1:n)   * r2 * r2 )

  p(1,1:n) = pc(1) + r(1:n) * cos ( theta(1:n) )
  p(2,1:n) = pc(2) + r(1:n) * sin ( theta(1:n) )

  return
end
subroutine uniform_in_circle01_map ( n, x )

!*****************************************************************************80
!
!! UNIFORM_IN_CIRCLE01_MAP maps uniform points into the unit circle.
!
!  Discussion:
!
!    The unit circle has center at the origin, and radius 1.
!
!    This routine is valid for spatial dimension DIM_NUM = 2.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the dimension of the space
!    (which must be 2).
!
!    Input, integer N, the number of points.
!
!    Output, real ( kind = rk ) X(DIM_NUM,N), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2
  integer n

  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) r(n)
  real ( kind = rk ) t(n)
  real ( kind = rk ) x(dim_num,n)

  call random_number ( harvest = r(1:n) )
  r(1:n) = sqrt ( r(1:n) )

  call random_number ( harvest = t(1:n) )
  t(1:n) = 2.0D+00 * pi * t(1:n)

  x(1,1:n) = r(1:n) * cos ( t(1:n) )
  x(2,1:n) = r(1:n) * sin ( t(1:n) )

  return
end
subroutine uniform_in_cube01 ( dim_num, n, x )

!*****************************************************************************80
!
!! UNIFORM_IN_CUBE01 creates uniform points in the unit hypercube.
!
!  Discussion:
!
!    The unit hypercube is defined as points whose components are between
!    0 and 1.
!
!    This routine is valid for any spatial dimension DIM_NUM.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the dimension of the space.
!
!    Input, integer N, the number of points.
!
!    Output, real ( kind = rk ) X(DIM_NUM,N), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num
  integer n

  real ( kind = rk ) x(dim_num,n)

  call random_number ( harvest = x(1:dim_num,1:n) )

  return
end
subroutine uniform_in_ellipsoid_map ( dim_num, n, a, r, x )

!*****************************************************************************80
!
!! UNIFORM_IN_ELLIPSOID_MAP maps uniform points into an ellipsoid.
!
!  Discussion:
!
!    The points X in the ellipsoid are described by a DIM_NUM by DIM_NUM
!    positive definite symmetric matrix A, and a "radius" R, such that
!
!      X' * A * X <= R * R
!
!    The algorithm computes the Cholesky factorization of A:
!
!      A = U' * U.
!
!    A set of uniformly random points Y is generated, satisfying:
!
!      Y' * Y <= R * R.
!
!    The appropriate points in the ellipsoid are found by solving
!
!      U * X = Y
!
!    Thanks to Dr Karl-Heinz Keil for pointing out that the original
!    coding was actually correct only if A was replaced by its inverse.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Reuven Rubinstein,
!    Monte Carlo Optimization, Simulation, and Sensitivity
!    of Queueing Networks,
!    Krieger, 1992,
!    ISBN: 0894647644,
!    LC: QA298.R79.
!
!  Parameters:
!
!    Input, integer DIM_NUM, the dimension of the space.
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) A(DIM_NUM,DIM_NUM), the matrix that describes
!    the ellipsoid.
!
!    Input, real ( kind = rk ) R, the right hand side of the ellipsoid equation.
!
!    Output, real ( kind = rk ) X(DIM_NUM,N), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num
  integer n

  real ( kind = rk ) a(dim_num,dim_num)
  integer info
  integer j
  real ( kind = rk ) r
  real ( kind = rk ) u(dim_num,dim_num)
  real ( kind = rk ) x(dim_num,n)
!
!  Get the Cholesky factor U.
!
  u(1:dim_num,1:dim_num) = a(1:dim_num,1:dim_num)

  call r8po_fa ( dim_num, u, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'UNIFORM_IN_ELLIPSOID_MAP - Fatal error!'
    write ( *, '(a)' ) '  R8PO_FA reports that the matrix A '
    write ( *, '(a)' ) '  is not positive definite symmetric.'
    stop 1
  end if
!
!  Get the points Y that satisfy Y' * Y <= R * R.
!
  call uniform_in_sphere01_map ( dim_num, n, x )

  x(1:dim_num,1:n) = r * x(1:dim_num,1:n)
!
!  Solve U * X = Y.
!
  do j = 1, n
    call r8po_sl ( dim_num, u, x(1:dim_num,j) )
  end do

  return
end
subroutine uniform_in_parallelogram_map ( v1, v2, v3, n, x )

!*****************************************************************************80
!
!! UNIFORM_IN_PARALLELOGRAM_MAP maps uniform points into a parallelogram.
!
!  Discussion:
!
!    The parallelogram is defined by three vertices, V1, V2 and V3.
!    The missing vertex V4 is equal to V2+V3-V1.
!
!    This routine is valid for spatial dimension DIM_NUM = 2.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Greg Turk,
!    Generating Random Points in a Triangle,
!    in Graphics Gems,
!    edited by Andrew Glassner,
!    AP Professional, 1990, pages 24-28.
!
!  Parameters:
!
!    Input, real ( kind = rk ) V1(2), V2(2), V3(2), the vertices.
!
!    Input, integer N, the number of points.
!
!    Output, real ( kind = rk ) X(2,N), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2
  integer n

  integer j
  real ( kind = rk ) r(2)
  real ( kind = rk ) v1(dim_num)
  real ( kind = rk ) v2(dim_num)
  real ( kind = rk ) v3(dim_num)
  real ( kind = rk ) x(dim_num,n)

  do j = 1, n

     call random_number ( harvest = r(1:2) )

    x(1:dim_num,j) = ( 1.0D+00 - r(1) - r(2) ) * v1(1:dim_num) &
                               + r(1)          * v2(1:dim_num) &
                                      + r(2)   * v3(1:dim_num)

  end do

  return
end
subroutine uniform_in_polygon_convex ( nv, v, n, x )

!*****************************************************************************80
!
!! uniform_in_polygon_convex() uniformly samples points in a convex polygon.
!
!  Discussion:
!
!    If the polygon is regular, or convex, or at least star-shaped,
!    this routine will work.
!
!    This routine assumes that all points between the centroid and
!    any point on the boundary lie within the polygon.
!
!    This routine is valid for spatial dimension DIM_NUM = 2.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 January 2025
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer NV, the number of vertices.
!
!    real ( kind = rk ) V(2,NV), the vertices of the polygon, listed in
!    clockwise or counterclockwise order.
!
!    integer N, the number of points to create.
!
!  Output:
!
!    real ( kind = rk ) X(2,N), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2
  integer n
  integer nv

  real ( kind = rk ) area(nv)
  real ( kind = rk ) area_percent
  real ( kind = rk ) centroid(dim_num)
  integer i
  integer ip1
  integer j
  integer k
  real ( kind = rk ) r(2)
  real ( kind = rk ) t(dim_num,3)
  real ( kind = rk ) x(dim_num,n)
  real ( kind = rk ) v(dim_num,nv)
!
!  Find the centroid.
!
  call polygon_centroid_2d ( nv, v, centroid )
!
!  Determine the areas of each triangle.
!
  do i = 1, nv

    if ( i < nv ) then
      ip1 = i + 1
    else
      ip1 = 1
    end if

    t(1:2,1) = v(1:2,i)
    t(1:2,2) = v(1:2,ip1)
    t(1:2,3) = centroid(1:2)

    call triangle_area_2d ( t, area(i) )

  end do
!
!  Normalize the areas.
!
  area(1:nv) = area(1:nv) / sum ( area(1:nv) )
!
!  Replace each area by the sum of itself and all previous ones.
!
  do i = 2, nv
    area(i) = area(i) + area(i-1)
  end do

  do j = 1, n
!
!  Choose a triangle at random, based on areas.
!
    call random_number ( harvest = area_percent )

    do k = 1, nv

      if ( area_percent <= area(k) ) then
        i = k
        exit
      end if

    end do
!
!  Now choose a point at random in the triangle.
!
    if ( i < nv ) then
      ip1 = i + 1
    else
      ip1 = 1
    end if

    call random_number ( harvest = r(1:dim_num) )

    if ( 1.0D+00 < sum ( r(1:dim_num) ) ) then
      r(1:dim_num) = 1.0D+00 - r(1:dim_num)
    end if

    x(1:dim_num,j) = ( 1.0D+00 - r(1) - r(2) ) * v(1:dim_num,i) &
                               + r(1)          * v(1:dim_num,ip1) &
                                      + r(2)   * centroid(1:dim_num)

  end do

  return
end
subroutine uniform_in_sector_map ( r1, r2, t1, t2, n, x )

!*****************************************************************************80
!
!! UNIFORM_IN_SECTOR_MAP maps uniform points into a circular sector.
!
!  Discussion:
!
!    The sector lies between circles with center at 0 and radius R1 and R2,
!    and between rays from the center at the angles T1 and T2.
!
!    This routine is valid for spatial dimension DIM_NUM = 2.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Peter Shirley,
!    Nonuniform Random Point Sets Via Warping,
!    Graphics Gems, Volume III,
!    edited by David Kirk,
!    AP Professional, 1992,
!    ISBN: 0122861663,
!    LC: T385.G6973.
!
!  Parameters:
!
!    Input, real ( kind = rk ) R1, R2, the two radii.
!
!    Input, real ( kind = rk ) T1, T2, the two angles, which should
!    be measured in radians, with T1 < T2.
!
!    Input, integer N, the number of points.
!
!    Output, real ( kind = rk ) X(2,N), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2
  integer n

  real ( kind = rk ) r(n)
  real ( kind = rk ) r1
  real ( kind = rk ) r2
  real ( kind = rk ) t(n)
  real ( kind = rk ) t1
  real ( kind = rk ) t2
  real ( kind = rk ) u(n)
  real ( kind = rk ) v(n)
  real ( kind = rk ) x(dim_num,n)

  call random_number ( harvest = u(1:n) )
  call random_number ( harvest = v(1:n) )

  t(1:n) =        ( 1.0D+00 - u(1:n) ) * t1    + u(1:n) * t2
  r(1:n) = sqrt ( ( 1.0D+00 - v(1:n) ) * r1**2 + v(1:n) * r2**2 )

  x(1,1:n) = r(1:n) * cos ( t(1:n) )
  x(2,1:n) = r(1:n) * sin ( t(1:n) )

  return
end
subroutine uniform_in_simplex01_map ( dim_num, point_num, x )

!*****************************************************************************80
!
!! UNIFORM_IN_SIMPLEX01_MAP maps uniform points into the unit simplex.
!
!  Discussion:
!
!    The interior of the unit DIM_NUM-dimensional simplex is the set of
!    points X(1:DIM_NUM) such that each X(I) is nonnegative, and
!    sum(X(1:DIM_NUM)) <= 1.
!
!    This routine is valid for any spatial dimension DIM_NUM.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Reuven Rubinstein,
!    Monte Carlo Optimization, Simulation, and Sensitivity
!    of Queueing Networks,
!    Krieger, 1992,
!    ISBN: 0894647644,
!    LC: QA298.R79.
!
!  Parameters:
!
!    Input, integer DIM_NUM, the dimension of the space.
!
!    Input, integer POINT_NUM, the number of points.
!
!    Output, real ( kind = rk ) X(DIM_NUM,POINT_NUM), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num
  integer point_num

  real ( kind = rk ) e(dim_num+1)
  integer j
  real ( kind = rk ) x(dim_num,point_num)
!
!  The construction begins by sampling DIM_NUM+1 points from the
!  exponential distribution with parameter 1.
!
  do j = 1, point_num

    call random_number ( harvest = e(1:dim_num+1) )

    e(1:dim_num+1) = -log ( e(1:dim_num+1) )

    x(1:dim_num,j) = e(1:dim_num) / sum ( e(1:dim_num+1) )

  end do

  return
end
subroutine uniform_in_sphere01_map ( dim_num, n, x )

!*****************************************************************************80
!
!! UNIFORM_IN_SPHERE01_MAP maps uniform points into the unit sphere.
!
!  Discussion:
!
!    The sphere has center 0 and radius 1.
!
!    This routine is valid for any spatial dimension DIM_NUM.
!
!    We first generate a point ON the sphere, and then distribute it
!    IN the sphere.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 August 2004
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
!    Reuven Rubinstein,
!    Monte Carlo Optimization, Simulation, and Sensitivity
!    of Queueing Networks,
!    Krieger, 1992,
!    ISBN: 0894647644,
!    LC: QA298.R79.
!
!  Parameters:
!
!    Input, integer DIM_NUM, the dimension of the space.
!
!    Input, integer N, the number of points.
!
!    Output, real ( kind = rk ) X(DIM_NUM,N), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num
  integer n

  real ( kind = rk ) exponent
  integer j
  real ( kind = rk ) norm
  real ( kind = rk ) r
  real ( kind = rk ) x(dim_num,n)

  exponent = 1.0D+00 / real ( dim_num, kind = rk )

  do j = 1, n
!
!  Fill a vector with normally distributed values.
!
    call r8vec_normal_01 ( dim_num, x(1:dim_num,j) )
!
!  Compute the length of the vector.
!
    norm = sqrt ( sum ( x(1:dim_num,j)**2 ) )
!
!  Normalize the vector.
!
    x(1:dim_num,j) = x(1:dim_num,j) / norm
!
!  Now compute a value to map the point ON the sphere INTO the sphere.
!
    call random_number ( harvest = r )

    x(1:dim_num,j) = r ** exponent * x(1:dim_num,j)

  end do

  return
end
subroutine uniform_in_tetrahedron ( v, n, x )

!*****************************************************************************80
!
!! UNIFORM_IN_TETRAHEDRON returns uniform points in a tetrahedron.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Claudio Rocchini, Paolo Cignoni,
!    Generating Random Points in a Tetrahedron,
!    Journal of Graphics Tools,
!    Volume 5, Number 5, 2000, pages 9-12.
!
!  Parameters:
!
!    Input, real ( kind = rk ) V(3,4), the vertices of the tetrahedron.
!
!    Input, integer N, the number of points.
!
!    Output, real ( kind = rk ) X(3,N), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) c(4)
  integer j
  real ( kind = rk ) t
  real ( kind = rk ) v(3,4)
  real ( kind = rk ) x(3,n)

  do j = 1, n

    call random_number ( harvest = c(1:3) )

    if ( 1.0D+00 < c(1) + c(2) ) then
      c(1) = 1.0D+00 - c(1)
      c(2) = 1.0D+00 - c(2)
    end if

    if ( 1.0D+00 < c(2) + c(3) ) then
      t = c(3)
      c(3) = 1.0D+00 - c(1) - c(2)
      c(2) = 1.0D+00 - t
    else if ( 1.0D+00 < c(1) + c(2) + c(3) ) then
       t = c(3)
       c(3) = c(1) + c(2) + c(3) - 1.0D+00
       c(1) = 1.0D+00 - c(2) - t
    end if

    c(4) = 1.0D+00 - c(1) - c(2) - c(3)
!
!  C(1:4) are the barycentric coordinates of the point.
!
    x(1:3,j) = matmul ( v(1:3,1:4), c(1:4) )

  end do

  return
end
subroutine uniform_in_triangle_map1 ( v1, v2, v3, n, x )

!*****************************************************************************80
!
!! UNIFORM_IN_TRIANGLE_MAP1 maps uniform points into a triangle.
!
!  Discussion:
!
!    The triangle is defined by three vertices.  This routine
!    uses Turk's rule #1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Greg Turk,
!    Generating Random Points in a Triangle,
!    in Graphics Gems,
!    edited by Andrew Glassner,
!    AP Professional, 1990, pages 24-28.
!
!  Parameters:
!
!    Input, real ( kind = rk ) V1(2), V2(2), V3(2), the vertices.
!
!    Input, integer N, the number of points.
!
!    Output, real ( kind = rk ) X(2,N), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2
  integer n

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  integer j
  real ( kind = rk ) r(3)
  real ( kind = rk ) x(dim_num,n)
  real ( kind = rk ) v1(dim_num)
  real ( kind = rk ) v2(dim_num)
  real ( kind = rk ) v3(dim_num)
!
!  Generate the points using Turk's rule 1.
!
  do j = 1, n

    call random_number ( harvest = r(1:2) )

    a = 1.0D+00            - sqrt ( r(2) )
    b = ( 1.0D+00 - r(1) ) * sqrt ( r(2) )
    c =             r(1)   * sqrt ( r(2) )

    x(1:dim_num,j) = a * v1(1:dim_num) &
                   + b * v2(1:dim_num) &
                   + c * v3(1:dim_num)

  end do

  return
end
subroutine uniform_in_triangle_map2 ( v1, v2, v3, n, x )

!*****************************************************************************80
!
!! UNIFORM_IN_TRIANGLE_MAP2 maps uniform points into a triangle.
!
!  Discussion:
!
!    The triangle is defined by three vertices.  This routine
!    uses Turk's rule #2.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Greg Turk,
!    Generating Random Points in a Triangle,
!    in Graphics Gems,
!    edited by Andrew Glassner,
!    AP Professional, 1990, pages 24-28.
!
!  Parameters:
!
!    Input, real ( kind = rk ) V1(2), V2(2), V3(2), the vertices.
!
!    Input, integer N, the number of points.
!
!    Output, real ( kind = rk ) X(2,N), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2
  integer n

  integer j
  real ( kind = rk ) r(3)
  real ( kind = rk ) x(dim_num,n)
  real ( kind = rk ) v1(dim_num)
  real ( kind = rk ) v2(dim_num)
  real ( kind = rk ) v3(dim_num)

  do j = 1, n

    call random_number ( harvest = r(1:2) )

    if ( 1.0D+00 < r(1) + r(2) ) then
      r(1:2) = 1.0D+00 - r(1:2)
    end if

    x(1:dim_num,j) = ( 1.0D+00 - r(1) - r(2) ) * v1(1:dim_num) &
                               + r(1)          * v2(1:dim_num) &
                                      + r(2)   * v3(1:dim_num)

  end do

  return
end
subroutine uniform_in_triangle01_map ( n, x )

!*****************************************************************************80
!
!! UNIFORM_IN_TRIANGLE01_MAP maps uniform points into the unit triangle.
!
!  Discussion:
!
!    The triangle is defined by the three vertices (1,0), (0,1) and (0,0).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!
!    Output, real ( kind = rk ) X(2,N), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2
  integer n

  integer j
  real ( kind = rk ) r(dim_num)
  real ( kind = rk ) x(dim_num,n)
!
!  Generate the points using barycentric coordinates.
!
  do j = 1, n

    call random_number ( harvest = r(1:dim_num) )

    if ( 1.0D+00 < sum ( r(1:dim_num) ) ) then
      r(1:dim_num) = 1.0D+00 - r(1:dim_num)
    end if

    x(1:dim_num,j) = r(1:dim_num)

  end do

  return
end
subroutine uniform_on_cube ( m, n, c, r, x )

!*****************************************************************************80
!
!! UNIFORM_ON_CUBE returns random points on the surface of a cube.
!
!  Discussion:
!
!    The cube is assumed to be aligned with the coordinate axes.
!
!    The cube has center C and radius R.  Any point on the surface of
!    the cube is described by
!
!      X = C + R * PM
!
!    where PM is an M-dimensional vector whose entries are between
!    -1 and +1, and for which at least one value has norm 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 April 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the spatial dimension.
!    1 <= M.
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) C(M), the coordinates of the center.
!
!    Input, real ( kind = rk ) R, the radius.
!
!    Output, real ( kind = rk ) X(M,N), the coordinates of N points, chosen
!    uniformly at random from the surface of the M-cube of center C and 
!    radius R.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) c(m)
  integer i
  integer i4_uniform_ab
  integer j
  integer k
  real ( kind = rk ) r
  real ( kind = rk ) x(m,n)
!
!  Choose random points within the cube of radius 1.
!
  call random_number ( harvest = x(1:m,1:n) )

  x(1:m,1:n) = 2.0D+00 * x(1:m,1:n) - 1.0D+00
!
!  For each point, select a coordinate at random, and set it to +1 or -1.
!
  do j = 1, n 
    i = i4_uniform_ab ( 1, m )
    k = i4_uniform_ab ( 0, 1 )
    if ( k == 0 ) then
      x(i,j) = 1.0D+00
    else
      x(i,j) = -1.0D+00
    end if
  end do
!
!  Shift by C and scale by R.
!
  do j = 1, n
    x(1:m,j) = c(1:m) + r * x(1:m,j)
  end do

  return
end
subroutine uniform_on_cube01 ( m, n, x )

!*****************************************************************************80
!
!! UNIFORM_ON_CUBE01 returns random points on the surface of the unit cube.
!
!  Discussion:
!
!    The cube is assumed to be aligned with the coordinate axes.
!
!    The cube has center at the origin and radius 1. Any point on the 
!    surface of is described by
!
!      X = PM
!
!    where PM is an M-dimensional vector whose entries are between
!    -1 and +1, and for which at least one value has norm 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 April 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the spatial dimension.
!    1 <= M.
!
!    Input, integer N, the number of points.
!
!    Output, real ( kind = rk ) X(M,N), the coordinates of N points, chosen
!    uniformly at random from the surface of the unit M-cube.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer i
  integer i4_uniform_ab
  integer j
  integer k
  real ( kind = rk ) x(m,n)
!
!  Choose random points within the cube of radius 1.
!
  call random_number ( harvest = x(1:m,1:n) )

  x(1:m,1:n) = 2.0D+00 * x(1:m,1:n) - 1.0D+00
!
!  For each point, select a coordinate at random, and set it to +1 or -1.
!
  do j = 1, n 
    i = i4_uniform_ab ( 1, m )
    k = i4_uniform_ab ( 0, 1 )
    if ( k == 0 ) then
      x(i,j) = 1.0D+00
    else
      x(i,j) = -1.0D+00
    end if
  end do

  return
end
subroutine uniform_on_ellipsoid_map ( dim_num, n, a, r, x )

!*****************************************************************************80
!
!! UNIFORM_ON_ELLIPSOID_MAP maps uniform points onto an ellipsoid.
!
!  Discussion:
!
!    The points X on the ellipsoid are described by an M by M positive
!    definite symmetric matrix A, and a "radius" R, such that
!
!      X' * A * X = R * R
!
!    The algorithm computes the Cholesky factorization of A:
!
!      A = U' * U.
!
!    A set of uniformly random points Y is generated, satisfying:
!
!      Y' * Y = R * R.
!
!    The appropriate points in the ellipsoid are found by solving
!
!      U * X = Y
!
!    Thanks to Dr Karl-Heinz Keil for pointing out that the original
!    coding was actually correct only if A was replaced by its inverse.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Reuven Rubinstein,
!    Monte Carlo Optimization, Simulation, and Sensitivity
!    of Queueing Networks,
!    Krieger, 1992,
!    ISBN: 0894647644,
!    LC: QA298.R79.
!
!  Parameters:
!
!    Input, integer DIM_NUM, the dimension of the space.
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) A(DIM_NUM,DIM_NUM), the matrix that describes
!    the ellipsoid.
!
!    Input, real ( kind = rk ) R, the right hand side of the ellipsoid equation.
!
!    Output, real ( kind = rk ) X(DIM_NUM,N), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num
  integer n

  real ( kind = rk ) a(dim_num,dim_num)
  integer info
  integer j
  real ( kind = rk ) r
  real ( kind = rk ) u(dim_num,dim_num)
  real ( kind = rk ) x(dim_num,n)
!
!  Get the factor U.
!
  u(1:dim_num,1:dim_num) = a(1:dim_num,1:dim_num)

  call r8po_fa ( dim_num, u, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'UNIFORM_ON_ELLIPSOID_MAP - Fatal error!'
    write ( *, '(a)' ) '  R8PO_FA reports that the matrix A '
    write ( *, '(a)' ) '  is not positive definite symmetric.'
    stop 1
  end if
!
!  Get the points Y that satisfy Y' * Y = R * R.
!
  call uniform_on_sphere01_map ( dim_num, n, x )

  x(1:dim_num,1:n) = r * x(1:dim_num,1:n)
!
!  Solve U * X = Y.
!
  do j = 1, n
    call r8po_sl ( dim_num, u, x(1:dim_num,j) )
  end do

  return
end
subroutine uniform_on_hemisphere01_phong ( n, m, x )

!*****************************************************************************80
!
!! uniform_on_hemisphere01_phong() maps uniform points onto the unit hemisphere.
!
!  Discussion:
!
!    The sphere has center 0 and radius 1.
!
!    The Phong density is used, with exponent M:
!
!    rho ( theta, phi; m ) = ( m + 1 ) * cos ( phi )**M / ( 2 * pi )
!
!    This routine is valid for spatial dimension DIM_NUM = 3.
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
!  Reference:
!
!    Peter Shirley,
!    Nonuniform Random Point Sets Via Warping,
!    Graphics Gems, Volume III,
!    edited by David Kirk,
!    AP Professional, 1992,
!    ISBN: 0122861663,
!    LC: T385.G6973.
!
!  Parameters:
!
!    Input, integer N, the number of points.
!
!    Input, integer M, the Phong exponent.
!
!    Output, real ( kind = rk ) X(3,N), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer n

  integer m
  real ( kind = rk ) phi(n)
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) power
  real ( kind = rk ) theta(n)
  real ( kind = rk ) x(dim_num,n)

  power = 1.0D+00 / real ( m + 1, kind = rk )
  call random_number ( harvest = phi(1:n) )
  phi(1:n) = acos ( ( 1.0D+00 - phi(1:n) )**power )

  call random_number ( harvest = theta(1:n) )
  theta(1:n) = 2.0D+00 * pi * theta(1:n)

  x(1,1:n) = cos ( theta(1:n) ) * sin ( phi(1:n) )
  x(2,1:n) = sin ( theta(1:n) ) * sin ( phi(1:n) )
  x(3,1:n) = cos ( phi(1:n) )

  return
end
subroutine uniform_on_simplex01_map ( dim_num, n, x )

!*****************************************************************************80
!
!! uniform_on_simplex01_map() maps uniform points onto the unit simplex.
!
!  Discussion:
!
!    The surface of the unit DIM_NUM-dimensional simplex is the set of points
!    X(1:DIM_NUM) such that each X(I) is nonnegative,
!    every X(I) is no greater than 1, and
!
!    ( X(I) = 0 for some I, or sum ( X(1:M) ) = 1. )
!
!    In DIM_NUM dimensions, there are DIM_NUM sides, and one main face.
!    This code picks a point uniformly with respect to "area".
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 April 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Reuven Rubinstein,
!    Monte Carlo Optimization, Simulation, and Sensitivity
!    of Queueing Networks,
!    Krieger, 1992,
!    ISBN: 0894647644,
!    LC: QA298.R79.
!
!  Parameters:
!
!    Input, integer DIM_NUM, the dimension of the space.
!
!    Input, integer N, the number of points.
!
!    Output, real ( kind = rk ) X(DIM_NUM,N), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num
  integer n

  real ( kind = rk ) area1
  real ( kind = rk ) area2
  real ( kind = rk ) e(dim_num)
  integer i
  integer i4_uniform_ab
  integer j
  real ( kind = rk ) r
  real ( kind = rk ) x(dim_num,n)
!
!  The construction begins by sampling DIM_NUM points from the
!  exponential distribution with parameter 1.
!
  do j = 1, n

    call random_number ( harvest = e(1:dim_num) )

    e(1:dim_num) = -log ( e(1:dim_num) )
!
!  Based on their relative areas, choose a side of the simplex,
!  or the main face.
!
    x(1:dim_num,j) = e(1:dim_num) / sum ( e(1:dim_num) )

    area1 = sqrt ( real ( dim_num, kind = rk ) )

    area2 = real ( dim_num, kind = rk )

    call random_number ( harvest = r )

    if ( area1 / ( area1 + area2 ) < r ) then
      i = i4_uniform_ab ( 1, dim_num )
      x(i,j) = 0.0D+00
    end if

  end do

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
  call r8mat_normal_01 ( dim_num, n, x )
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
subroutine uniform_on_sphere01_patch_tp ( n, phi1, phi2, theta1, theta2, tp )

!*****************************************************************************80
!
!! uniform_on_sphere01_patch_tp() maps uniform points onto a spherical TP patch.
!
!  Discussion:
!
!    The sphere has center 0 and radius 1.
!
!    A spherical TP patch on the surface of the unit sphere contains those
!    points with radius R = 1 and angles (THETA,PHI) such that
!
!      0.0 <= THETA1 <= THETA <= THETA2 <= 2 * PI
!      0.0 <= PHI1   <= PHI   <= PHI2   <=     PI
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Peter Shirley,
!    Nonuniform Random Point Sets Via Warping,
!    Graphics Gems, Volume III,
!    edited by David Kirk,
!    AP Professional, 1992,
!    ISBN: 0122861663,
!    LC: T385.G6973.
!
!  Parameters:
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) PHI1, PHI2, the latitudinal angle range.
!
!    Input, real ( kind = rk ) THETA1, THETA2, the longitudinal angle range.
!
!    Output, real ( kind = rk ) TP(2,N), the THETA, PHI coordinates of
!    the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) phi1
  real ( kind = rk ) phi2
  real ( kind = rk ) theta1
  real ( kind = rk ) theta2
  real ( kind = rk ) tp(2,n)

  call random_number ( harvest = tp(1:2,1:n) )

  tp(1,1:n) = ( 1.0D+00 - tp(1,1:n) ) * theta1 &
             +            tp(1,1:n)   * theta2

  tp(2,1:n) = acos ( ( 1.0D+00 - tp(2,1:n) ) * cos ( phi1 ) &
                  +              tp(2,1:n)   * cos ( phi2 ) )

  return
end
subroutine uniform_on_sphere01_patch_xyz ( n, phi1, phi2, theta1, theta2, x )

!*****************************************************************************80
!
!! UNIFORM_ON_SPHERE01_PATCH_XYZ maps uniform points to a spherical XYZ patch.
!
!  Discussion:
!
!    The sphere has center 0 and radius 1.
!
!    A sphere XYZ patch on the surface of the unit sphere contains those
!    points with radius R = 1 and angles (THETA,PHI) such that
!
!      0.0 <= THETA1 <= THETA <= THETA2 <= 2 * PI
!      0.0 <= PHI1   <= PHI   <= PHI2   <=     PI
!
!    transformed into Cartesian XYZ coordinates:
!
!      X = cos ( THETA ) * sin ( PHI )
!      Y = sin ( THETA ) * sin ( PHI )
!      Z =                 cos ( PHI )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Peter Shirley,
!    Nonuniform Random Point Sets Via Warping,
!    Graphics Gems, Volume III,
!    edited by David Kirk,
!    AP Professional, 1992,
!    ISBN: 0122861663,
!    LC: T385.G6973.
!
!  Parameters:
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) PHI1, PHI2, the latitudinal angle range.
!
!    Input, real ( kind = rk ) THETA1, THETA2, the longitudinal angle range.
!
!    Output, real ( kind = rk ) X(3,N), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3
  integer n

  real ( kind = rk ) phi(n)
  real ( kind = rk ) phi1
  real ( kind = rk ) phi2
  real ( kind = rk ) theta(n)
  real ( kind = rk ) theta1
  real ( kind = rk ) theta2
  real ( kind = rk ) x(dim_num,n)

  call random_number ( harvest = phi(1:n) )
  phi(1:n) = acos ( ( 1.0D+00 - phi(1:n) ) * cos ( phi1 ) &
                  +             phi(1:n)   * cos ( phi2 ) )

  call random_number ( harvest = theta(1:n) )
  theta(1:n) = ( 1.0D+00 - theta(1:n) ) * theta1 &
             +             theta(1:n)   * theta2

  x(1,1:n) = cos ( theta(1:n) ) * sin ( phi(1:n) )
  x(2,1:n) = sin ( theta(1:n) ) * sin ( phi(1:n) )
  x(3,1:n) = cos ( phi(1:n) )

  return
end
subroutine uniform_on_sphere01_triangle_xyz ( n, v1, v2, v3, x )

!*****************************************************************************80
!
!! UNIFORM_ON_SPHERE01_TRIANGLE_XYZ: sample spherical triangle, XYZ coordinates.
!
!  Discussion:
!
!    The sphere has center 0 and radius 1.
!
!    A spherical triangle on the surface of the unit sphere contains those
!    points with radius R = 1, bounded by the vertices V1, V2, V3.
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
!    John Burkardt
!
!  Reference:
!
!    James Arvo,
!    Stratified sampling of spherical triangles,
!    Computer Graphics Proceedings, Annual Conference Series,
!    ACM SIGGRAPH '95, pages 437-438, 1995.
!
!  Parameters:
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) V1(3), V2(3), V3(3), the XYZ coordinates of
!    the vertices of the spherical triangle.
!
!    Output, real ( kind = rk ) X(3,N), the XYZ coordinates of the
!    sample points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a
  real ( kind = rk ) alpha
  real ( kind = rk ) area
  real ( kind = rk ) area_hat
  real ( kind = rk ) b
  real ( kind = rk ) beta
  real ( kind = rk ) c
  real ( kind = rk ) gamma
  integer j
  real ( kind = rk ) q
  real ( kind = rk ) r
  real ( kind = rk ) r8vec_norm
  real ( kind = rk ) s
  real ( kind = rk ) t
  real ( kind = rk ) u
  real ( kind = rk ) v
  real ( kind = rk ) v1(3)
  real ( kind = rk ) v2(3)
  real ( kind = rk ) v3(3)
  real ( kind = rk ) v31(3)
  real ( kind = rk ) v4(3)
  real ( kind = rk ) v42(3)
  real ( kind = rk ) w
  real ( kind = rk ) x(3,n)
  real ( kind = rk ) xsi1
  real ( kind = rk ) xsi2
  real ( kind = rk ) z
!
!  Compute the sides, angles, and area of the spherical triangle;
!  for now, we assume R = 1.
!
  r = 1.0D+00

  call stri_vertices_to_sides ( r, v1, v2, v3, a, b, c )

  call stri_sides_to_angles ( r, a, b, c, alpha, beta, gamma )

  call stri_angles_to_area ( r, alpha, beta, gamma, area )

  do j = 1, n
!
!  Select the new area.
!
    call random_number ( harvest = xsi1 )
    area_hat = xsi1 * area
!
!  Compute the sine and cosine of the angle phi.
!
    s = sin ( area_hat - alpha )
    t = cos ( area_hat - alpha )
!
!  Compute the pair that determines beta_hat.
!
    u = t - cos ( alpha )
    v = s + sin ( alpha ) * cos ( c )
!
!  Q is the cosine of the new edge length b_hat.
!
    q = ( ( v * t - u * s ) * cos ( alpha ) - v ) &
      / ( ( v * s + u * t ) * sin ( alpha ) )
!
!  V31 = normalized ( V3 - ( V3 dot V1 ) * V1 )
!
    w = dot_product ( v3, v1 )
    v31(1:3) = v3(1:3) - w * v1(1:3)
    v31(1:3) = v31(1:3) / r8vec_norm ( 3, v31(1:3) )
!
!  V4 is the third vertex of the subtriangle V1, V2, V4.
!
    v4(1:3) = q * v1(1:3) + sqrt ( 1.0D+00 - q * q ) * v31(1:3)
!
!  Select cos theta, which will sample along the edge from V2 to V4.
!
    call random_number ( harvest = xsi2 )
    z = 1.0D+00 - xsi2 * ( 1.0D+00 - dot_product ( v4, v2 ) )
!
!  V42 = normalized ( V4 - ( V4 dot V2 ) * V2 )
!
    w = dot_product ( v4, v2 )
    v42(1:3) = v4(1:3) - w * v2(1:3)
    v42(1:3) = v42(1:3) / r8vec_norm ( 3, v42(1:3) )
!
!  Construct the point.
!
    x(1:3,j) = z * v2(1:3) + sqrt ( 1.0D+00 - z * z ) * v42(1:3)

  end do

  return
end
subroutine uniform_on_triangle ( n, v, x )

!*****************************************************************************80
!
!! UNIFORM_ON_TRIANGLE maps uniform points onto the boundary of a triangle.
!
!  Discussion:
!
!    The triangle is defined by the three vertices V1, V2, V3.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) V(2,3), the vertices of the triangle.
!
!    Output, real ( kind = rk ) X(2,N), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 2
  integer n

  integer j
  real ( kind = rk ) l1
  real ( kind = rk ) l2
  real ( kind = rk ) l3
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) t
  real ( kind = rk ) v(2,3)
  real ( kind = rk ) x(m,n)

  l1 = sqrt ( ( v(1,2) - v(1,1) ) ** 2 &
            + ( v(2,2) - v(2,1) ) ** 2 )

  l2 = sqrt ( ( v(1,3) - v(1,2) ) ** 2 &
            + ( v(2,3) - v(2,2) ) ** 2 )

  l3 = sqrt ( ( v(1,1) - v(1,3) ) ** 2 &
            + ( v(2,1) - v(2,3) ) ** 2 )
  
  do j = 1, n
!
!  R can be regarded as the distance of the point on the perimeter,
!  as measured from the origin, along the perimeter.
!
    call random_number ( harvest = r )
    r = ( l1 + l2 + l3 ) * r
!
!  Case 1: between V1 and V2.
!
    if ( r <= l1 ) then

      s = ( l1 - r ) / l1
      t =        r   / l1
      x(1,j) = s * v(1,1) + t * v(1,2)
      x(2,j) = s * v(2,1) + t * v(2,2)
!
!  Case 2: between V2 and V3.
!
    else if ( r <= l1 + l2 ) then

      s = ( l2 - r + l1 ) / l2
      t = (      r - l1 ) / l2
      x(1,j) = s * v(1,2) + t * v(1,3)
      x(2,j) = s * v(2,2) + t * v(2,3)
!
!  Case 3: between V3 and V1.
!
    else

      s = ( l3 - r + l1 + l2 ) / l3
      t = (      r - l1 - l2 ) / l3
      x(1,j) = s * v(1,3) + t * v(1,1)
      x(2,j) = s * v(2,3) + t * v(2,1)

    end if

  end do

  return
end
subroutine uniform_on_triangle01 ( n, x )

!*****************************************************************************80
!
!! UNIFORM_ON_TRIANGLE01: uniform points on the boundary of the unit triangle.
!
!  Discussion:
!
!    The unit triangle is defined by the three vertices (1,0), (0,1) and (0,0).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!
!    Output, real ( kind = rk ) X(2,N), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 2
  integer n

  real ( kind = rk ) a
  real ( kind = rk ) b
  integer j
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) x(m,n)

  s = sqrt ( 2.0D+00 )

  a =   1.0D+00       / ( 2.0D+00 + s )
  b = ( 1.0D+00 + s ) / ( 2.0D+00 + s )
  
  do j = 1, n
!
!  R can be regarded as the distance of the point on the perimeter,
!  as measured from the origin, along the perimeter.
!
    call random_number ( harvest = r )
    r = ( 2.0D+00 + s ) * r
!
!  Case 1: between (0,0) and (1,0).
!
    if ( r <= a ) then

      x(1,j) = 0.0D+00
      x(2,j) = r
!
!  Case 2: between (1,0) and (0,1).
!
    else if ( r <= b ) then

      x(1,j) = 1.0D+00 - ( r - a ) * sqrt ( 2.0D+00 ) / 2.0D+00
      x(2,j) = 0.0D+00 + ( r - a ) * sqrt ( 2.0D+00 ) / 2.0D+00
!
!  Case 3: between (0,1) and (0,0).
!
    else

      x(1,j) = 0.0D+00
      x(2,j) = 1.0D+00 - ( r - b )

    end if

  end do

  return
end
subroutine uniform_walk ( dim_num, n, x )

!*****************************************************************************80
!
!! UNIFORM_WALK generates points on a uniform random walk.
!
!  Discussion:
!
!    The first point is at the origin.  Uniform random numbers are
!    generated to determine the direction of the next step, which
!    is always of length 1, and in a coordinate direction.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM the dimension of the space.
!
!    Input, integer N, the number of points.
!
!    Output, real ( kind = rk ) X(DIM_NUM,N), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim_num
  integer n

  real ( kind = rk ) dir(n-1)
  integer i
  integer j
  real ( kind = rk ) x(dim_num,n)

  x(1:dim_num,1) = 0.0D+00

  call random_number ( harvest = dir(1:n-1) )

  dir(1:n-1) = real ( 2 * dim_num, kind = rk ) * ( dir(1:n-1) - 0.5D+00 )

  do j = 2, n

    x(1:dim_num,j) = x(1:dim_num,j-1)

    i = nint ( abs ( dir(j-1) ) + 0.5D+00 )
    i = min ( i, dim_num )
    i = max ( i, 1 )

    if ( dir(j-1) < 0.0D+00 ) then
      x(i,j) = x(i,j) - 1.0D+00
    else
      x(i,j) = x(i,j) + 1.0D+00
    end if

  end do

  return
end
