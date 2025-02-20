subroutine hypercube_grid ( m, n, ns, a, b, c, x )

!*****************************************************************************80
!
!! HYPERCUBE_GRID: grid points over the interior of a hypercube in M dimensions.
!
!  Discussion:
!
!    In M dimensional space, a logically rectangular grid is to be created.
!    In the I-th dimension, the grid will use S(I) points.
!    The total number of grid points is 
!      N = product ( 1 <= I <= M ) S(I)
!
!    Over the interval [A(i),B(i)], we have 5 choices for grid centering:
!      1: 0,   1/3, 2/3, 1
!      2: 1/5, 2/5, 3/5, 4/5
!      3: 0,   1/4, 2/4, 3/4
!      4: 1/4, 2/4, 3/4, 1
!      5: 1/8, 3/8, 5/8, 7/8
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    29 August 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of points.
!    N = product ( 1 <= I <= M ) NS(I).
!
!    Input, integer NS(M), the number of points along 
!    each dimension.
!
!    Input, real ( kind = rk ) A(M), B(M), the endpoints for each dimension.
!
!    Input, integer C(M), the grid centering for each dimension.
!    1 <= C(*) <= 5.
!
!    Output, real ( kind = rk ) X(M,N) = X(M*S(1),S(2),...,S(M)), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m)
  real ( kind = rk ) b(m)
  integer c(m)
  integer i
  integer j
  integer ns(m)
  integer s
  real ( kind = rk ) x(m,n)
  real ( kind = rk ), allocatable :: xs(:)
!
!  Create the 1D grids in each dimension.
!
  do i = 1, m

    s = ns(i)

    allocate ( xs(1:s) )

    do j = 1, s

      if ( c(i) == 1 ) then

        if ( s == 1 ) then
          xs(j) = 0.5D+00 * ( a(i) + b(i) )
        else
          xs(j) = (   real ( s - j,     kind = rk ) * a(i)   &
                    + real (     j - 1, kind = rk ) * b(i) ) & 
                    / real ( s     - 1, kind = rk )
        end if
      else if ( c(i) == 2 ) then
        xs(j) = (   real ( s - j + 1, kind = rk ) * a(i)   &
                  + real (     j,     kind = rk ) * b(i) ) & 
                  / real ( s     + 1, kind = rk )
      else if ( c(i) == 3 ) then
        xs(j) = (   real ( s - j + 1, kind = rk ) * a(i)   &
                  + real (     j - 1, kind = rk ) * b(i) ) & 
                  / real ( s,         kind = rk )
      else if ( c(i) == 4 ) then
        xs(j) = (   real ( s - j, kind = rk ) * a(i)   &
                  + real (     j, kind = rk ) * b(i) ) & 
                  / real ( s,     kind = rk )
      else if ( c(i) == 5 ) then
        xs(j) = (   real ( 2 * s - 2 * j + 1, kind = rk ) * a(i)   &
                  + real (         2 * j - 1, kind = rk ) * b(i) ) & 
                  / real ( 2 * s,             kind = rk )
      end if

    end do

    call r8vec_direct_product ( i, s, xs, m, n, x )

    deallocate ( xs )

  end do

  return
end
function i4vec_product ( n, a )

!*****************************************************************************80
!
!! I4VEC_PRODUCT returns the product of the entries of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    In FORTRAN90, this facility is offered by the built in
!    PRODUCT function:
!
!      I4VEC_PRODUCT ( N, A ) = PRODUCT ( A(1:N) )
!
!    In MATLAB, this facility is offered by the built in
!    PROD function:
!
!      I4VEC_PRODUCT ( N, A ) = PROD ( A(1:N) )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input, integer A(N), the array.
!
!    Output, integer I4VEC_PRODUCT, the product of the entries.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer a(n)
  integer i4vec_product

  i4vec_product = product ( a(1:n) )

  return
end
subroutine r8mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
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
!    14 June 2004
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
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  character ( len = * ) title

  call r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
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
  integer i2
  integer i2hi
  integer i2lo
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)'
    return
  end if

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(i8,6x)' ) i
    end do

    write ( *, '(''  Row   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc
        i = i2lo - 1 + i2
        write ( ctemp(i2), '(g14.6)' ) a(i,j)
      end do

      write ( *, '(i5,a,5a14)' ) j, ':', ( ctemp(i), i = 1, inc )

    end do

  end do

  return
end
subroutine r8vec_direct_product ( factor_index, factor_order, factor_value, &
  factor_num, point_num, x )

!*****************************************************************************80
!
!! R8VEC_DIRECT_PRODUCT creates a direct product of R8VEC's.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    To explain what is going on here, suppose we had to construct
!    a multidimensional quadrature rule as the product of K rules
!    for 1D quadrature.
!
!    The product rule will be represented as a list of points and weights.
!
!    The J-th item in the product rule will be associated with
!      item J1 of 1D rule 1,
!      item J2 of 1D rule 2,
!      ...,
!      item JK of 1D rule K.
!
!    In particular,
!      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
!    and
!      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)
!
!    So we can construct the quadrature rule if we can properly
!    distribute the information in the 1D quadrature rules.
!
!    This routine carries out that task for the abscissas X.
!
!    Another way to do this would be to compute, one by one, the
!    set of all possible indices (J1,J2,...,JK), and then index
!    the appropriate information.  An advantage of the method shown
!    here is that you can process the K-th set of information and
!    then discard it.
!
!  Example:
!
!    Rule 1:
!      Order = 4
!      X(1:4) = ( 1, 2, 3, 4 )
!
!    Rule 2:
!      Order = 3
!      X(1:3) = ( 10, 20, 30 )
!
!    Rule 3:
!      Order = 2
!      X(1:2) = ( 100, 200 )
!
!    Product Rule:
!      Order = 24
!      X(1:24) =
!        ( 1, 10, 100 )
!        ( 2, 10, 100 )
!        ( 3, 10, 100 )
!        ( 4, 10, 100 )
!        ( 1, 20, 100 )
!        ( 2, 20, 100 )
!        ( 3, 20, 100 )
!        ( 4, 20, 100 )
!        ( 1, 30, 100 )
!        ( 2, 30, 100 )
!        ( 3, 30, 100 )
!        ( 4, 30, 100 )
!        ( 1, 10, 200 )
!        ( 2, 10, 200 )
!        ( 3, 10, 200 )
!        ( 4, 10, 200 )
!        ( 1, 20, 200 )
!        ( 2, 20, 200 )
!        ( 3, 20, 200 )
!        ( 4, 20, 200 )
!        ( 1, 30, 200 )
!        ( 2, 30, 200 )
!        ( 3, 30, 200 )
!        ( 4, 30, 200 )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer FACTOR_INDEX, the index of the factor being
!    processed.  The first factor processed must be factor 1!
!
!    Input, integer FACTOR_ORDER, the order of the factor.
!
!    Input, real ( kind = rk ) FACTOR_VALUE(FACTOR_ORDER), the factor values
!    for factor FACTOR_INDEX.
!
!    Input, integer FACTOR_NUM, the number of factors.
!
!    Input, integer POINT_NUM, the number of elements in the
!    direct product.
!
!    Input/output, real ( kind = rk ) X(FACTOR_NUM,POINT_NUM), the elements of
!    the direct product, which are built up gradually.
!
!  Local Parameters:
!
!    Local, integer START, the first location of a block of 
!    values to set.
!
!    Local, integer CONTIG, the number of consecutive values 
!    to set.
!
!    Local, integer SKIP, the distance from the current value 
!    of START to the next location of a block of values to set.
!
!    Local, integer REP, the number of blocks of values to set.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer factor_num
  integer factor_order
  integer point_num

  integer, save :: contig
  integer factor_index
  real ( kind = rk ) factor_value(factor_order)
  integer j
  integer k
  integer, save :: rep
  integer, save :: skip
  integer start
  real ( kind = rk ) x(factor_num,point_num)

  if ( factor_index == 1 ) then
    contig = 1
    skip = 1
    rep = point_num
    x(1:factor_num,1:point_num) = 0.0D+00
  end if

  rep = rep / factor_order
  skip = skip * factor_order

  do j = 1, factor_order

    start = 1 + ( j - 1 ) * contig

    do k = 1, rep
      x(factor_index,start:start+contig-1) = factor_value(j)
      start = start + skip
    end do

  end do

  contig = contig * factor_order

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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
