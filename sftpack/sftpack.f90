subroutine c4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! C4MAT_PRINT_SOME prints some of a C4MAT.
!
!  Discussion:
!
!    A C4MAT is a matrix of C4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns 
!    in the matrix.
!
!    Input, complex ( kind = 4 ) A(M,N), the matrix.
!
!    Input, integer ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: incx = 4
  integer m
  integer n

  complex ( kind = 4 ) a(m,n)
  character ( len = 20 ) ctemp(incx)
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
  complex ( kind = 4 ) zero

  zero = cmplx ( 0.0E+00, 0.0E+00, kind = 4 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of INCX.
!
  do j2lo = jlo, min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i10,10x)' ) j
    end do

    write ( *, '(a,4a20)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi
!
!  Print out (up to) INCX entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == zero ) then
          ctemp(j2) = '       0.0          '
        else if ( imag ( a(i,j) ) == 0.0E+00 ) then
          write ( ctemp(j2), '(g10.3,10x)' ) real ( a(i,j), kind = 4 )
        else
          write ( ctemp(j2), '(2g10.3)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,a1,4a20)' ) i, ':', ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine c4mat_sftb ( n1, n2, y, x )

!*****************************************************************************80
!
!! C4MAT_SFTB computes a "slow" backward Fourier transform of a C4MAT.
!
!  Discussion:
!
!    SFTF and SFTB are inverses of each other.  If we begin with data
!    X and apply SFTF to get Y, and then apply SFTB to Y,
!    we should get back the original X.
!
!    This routine is provided for illustration and testing.  It is inefficient
!    relative to optimized routines that use fast Fourier techniques.
!
!    For 0 <= I1 <= N1 - 1, 
!        0 <= I2 <= N2 - 1,
!
!      X(I1,I2) = Sum ( 0 <= K2 <= N2 - 1 ) Sum ( 0 <= K1 <= N1 - 1 ) 
!        Y(K1,K2) * exp ( 2 pi i I1 K1 / N1 ) * exp ( 2 pi i I2 K2 / N2 )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, the number of rows and columns of data.
!
!    Input, complex ( kind = 4 ) Y(0:N1-1,0:N2-1), the Fourier coefficients.
!
!    Output, complex ( kind = 4 ) X(0:N1-1,0:N2-1), the data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n1
  integer n2

  complex ( kind = 4 ) cs1
  complex ( kind = 4 ) cs2
  integer i1
  integer i2
  integer j1
  integer j2
  real ( kind = 4 ), parameter :: r4_pi = 3.1415926E+00
  real ( kind = 4 ) theta1
  real ( kind = 4 ) theta2
  complex ( kind = 4 ) x(0:n1-1,0:n2-1)
  complex ( kind = 4 ) y(0:n1-1,0:n2-1)

  x(0:n1-1,0:n2-1) = cmplx ( 0.0E+00, 0.0E+00, kind = 4 )

  do i2 = 0, n2 - 1
    do j2 = 0, n2 - 1
      theta2 = 2.0E+00 * r4_pi * real ( i2 * j2, kind = 4 ) &
        / real ( n2, kind = 4 )
      cs2 = cmplx ( cos ( theta2 ), - sin ( theta2 ), kind = 4 )
      do i1 = 0, n1 - 1
        do j1 = 0, n1 - 1
          theta1 = 2.0E+00 * r4_pi * real ( i1 * j1, kind = 4 ) &
            / real ( n1, kind = 4 )
          cs1 = cmplx ( cos ( theta1 ), - sin ( theta1 ), kind = 4 )
          x(i1,i2) = x(i1,i2) + y(j1,j2) * cs1 * cs2
        end do
      end do
    end do
  end do

  x(0:n1-1,0:n2-1) = x(0:n1-1,0:n2-1) / real ( n1 * n2, kind = 4 )

  return
end
subroutine c4mat_sftf ( n1, n2, x, y )

!*****************************************************************************80
!
!! C4MAT_SFTF computes a "slow" forward Fourier transform of a C4MAT.
!
!  Discussion:
!
!    SFTF and SFTB are inverses of each other.  If we begin with data
!    X and apply SFTF to get Y, and then apply SFTB to Y, 
!    we should get back the original X.
!
!    This routine is provided for illustration and testing.  It is inefficient
!    relative to optimized routines that use fast Fourier techniques.
!
!    For 0 <= I1 <= N1 - 1, 
!        0 <= I2 <= N2 - 1,
!
!      Y(I1,I2) = Sum ( 0 <= K2 <= N2 - 1 ) Sum ( 0 <= K1 <= N1 - 1 ) 
!        X(K1,K2) * exp ( - 2 pi i I1 K1 / N1 ) * exp ( - 2 pi i I2 K2 / N2 )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, the number of rows and columns of data.
!
!    Input, complex ( kind = 4 ) X(0:N1-1,0:N2-1), the data to be transformed.
!
!    Output, complex ( kind = 4 ) Y(0:N1-1,0:N2-1), the Fourier coefficients.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n1
  integer n2

  complex ( kind = 4 ) cs1
  complex ( kind = 4 ) cs2
  integer i1
  integer i2
  integer j1
  integer j2
  real ( kind = 4 ), parameter :: r4_pi = 3.1415926E+00
  real ( kind = 4 ) theta1
  real ( kind = 4 ) theta2
  complex ( kind = 4 ) x(0:n1-1,0:n2-1)
  complex ( kind = 4 ) y(0:n1-1,0:n2-1)

  y(0:n1-1,0:n2-1) = cmplx ( 0.0E+00, 0.0E+00, kind = 4 )

  do i2 = 0, n2 - 1
    do j2 = 0, n2 - 1
      theta2 = - 2.0E+00 * r4_pi * real ( i2 * j2, kind = 4 ) &
        / real ( n2, kind = 4 )
      cs2 = cmplx ( cos ( theta2 ), - sin ( theta2 ), kind = 4 )
      do i1 = 0, n1 - 1
        do j1 = 0, n1 - 1
          theta1 = - 2.0E+00 * r4_pi * real ( i1 * j1, kind = 4 ) &
            / real ( n1, kind = 4 )
          cs1 = cmplx ( cos ( theta1 ), - sin ( theta1 ), kind = 4 )
          y(i1,i2) = y(i1,i2) + x(j1,j2) * cs1 * cs2
        end do
      end do
    end do
  end do

  return
end
subroutine c4mat_uniform_01 ( m, n, seed, c )

!*****************************************************************************80
!
!! C4MAT_UNIFORM_01 returns a unit pseudorandom C4MAT.
!
!  Discussion:
!
!    A C4MAT is a matrix of C4's.
!
!    The angles should be uniformly distributed between 0 and 2 * PI,
!    the square roots of the radius uniformly distributed between 0 and 1.
!
!    This results in a uniform distribution of values in the unit circle.
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
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns 
!    in the matrix.
!
!    Input/output, integer SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, complex ( kind = 4 ) C(M,N), the pseudorandom complex matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  complex ( kind = 4 ) c(m,n)
  integer i
  integer, parameter :: i4_huge = 2147483647
  integer j
  real ( kind = 4 ) r
  integer k
  real ( kind = 4 ), parameter :: r4_pi = 3.1415926E+00
  integer seed
  real ( kind = 4 ) theta

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'C4MAT_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  do j = 1, n
    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge
      end if

      r = sqrt ( real ( seed, kind = 4 ) * 4.656612875E-10 )

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge
      end if

      theta = 2.0E+00 * r4_pi * ( real ( seed, kind = 4 ) * 4.656612875E-10 )

      c(i,j) = r * cmplx ( cos ( theta ), sin ( theta ), kind = 4 )

    end do

  end do

  return
end
subroutine c4vec_print_part ( n, a, max_print, title )

!*****************************************************************************80
!
!! C4VEC_PRINT_PART prints "part" of a C4VEC.
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
!    22 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries of the vector.
!
!    Input, complex ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, integer MAX_PRINT, the maximum number of lines
!    to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  complex ( kind = 4 ) a(n)
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
      write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)
    end do

  else if ( 3 <= max_print ) then

    do i = 1, max_print - 2
      write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)
    end do
    write ( *, '(a)' ) '  ........  ..............  ..............'
    i = n
    write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)

  else

    do i = 1, max_print - 1
      write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)
    end do
    i = max_print
    write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6,2x,a)' ) i, ':', a(i), &
      '...more entries...'

  end if

  return
end
subroutine c4vec_sftb ( n, y, x )

!*****************************************************************************80
!
!! C4VEC_SFTB computes a "slow" backward Fourier transform of a C4VEC.
!
!  Discussion:
!
!    SFTF and SFTB are inverses of each other.  If we begin with data
!    X and apply SFTF to get Y, and then apply SFTB to Y,
!    we should get back the original X.
!
!    This routine is provided for illustration and testing.  It is inefficient
!    relative to optimized routines that use fast Fourier techniques.
!
!    For 0 <= I <= N - 1
!
!      X(I) = 1/N * Sum ( 0 <= J <= N - 1 ) Y(J) * exp ( 2 pi i I J / N )
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
!    Input, integer N, the number of data values.
!
!    Input, complex ( kind = 4 ) Y(0:N-1), the Fourier coefficients.
!
!    Output, complex ( kind = 4 ) X(0:N-1), the data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer j
  real ( kind = 4 ), parameter :: r4_pi = 3.1415926E+00
  real ( kind = 4 ) theta
  complex ( kind = 4 ) x(0:n-1)
  complex ( kind = 4 ) y(0:n-1)

  x(0:n-1) = cmplx ( 0.0E+00, 0.0E+00, kind = 4 )

  do i = 0, n - 1
    do j = 0, n - 1
      theta = - 2.0E+00 * r4_pi * real ( i * j, kind = 4 ) / real ( n, kind = 4 )
      x(i) = x(i) + y(j) * cmplx ( cos ( theta ), sin ( theta ), kind = 4 )
    end do
  end do

  x(0:n-1) = x(0:n-1) / real ( n, kind = 4 )

  return
end
subroutine c4vec_sftf ( n, x, y )

!*****************************************************************************80
!
!! C4VEC_SFTF computes a "slow" forward Fourier transform of a C4VEC.
!
!  Discussion:
!
!    SFTF and SFTB are inverses of each other.  If we begin with data
!    X and apply SFTF to get Y, and then apply SFTB to Y, 
!    we should get back the original X.
!
!    This routine is provided for illustration and testing.  It is inefficient
!    relative to optimized routines that use fast Fourier techniques.
!
!    For 0 <= I <= N - 1
!
!      Y(I) = Sum ( 0 <= J <= N - 1 ) X(J) * exp ( - 2 pi i I J / N )
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
!    Input, integer N, the number of data values.
!
!    Input, complex ( kind = 4 ) X(0:N-1), the data to be transformed.
!
!    Output, complex ( kind = 4 ) Y(0:N-1), the Fourier coefficients.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer j
  real ( kind = 4 ), parameter :: r4_pi = 3.1415926E+00
  real ( kind = 4 ) theta
  complex ( kind = 4 ) x(0:n-1)
  complex ( kind = 4 ) y(0:n-1)

  y(0:n-1) = cmplx ( 0.0E+00, 0.0E+00, kind = 4 )

  do i = 0, n - 1
    do j = 0, n - 1
      theta = - 2.0E+00 * r4_pi * real ( i * j, kind = 4 ) / real ( n, kind = 4 )
      y(i) = y(i) + x(j) * cmplx ( cos ( theta ), - sin ( theta ), kind = 4 )
    end do
  end do

  return
end
subroutine c4vec_uniform_01 ( n, seed, c )

!*****************************************************************************80
!
!! C4VEC_UNIFORM_01 returns a unit pseudorandom C4VEC.
!
!  Discussion:
!
!    A C4VEC is a vector of C4's.
!
!    The angles should be uniformly distributed between 0 and 2 * PI,
!    the square roots of the radius uniformly distributed between 0 and 1.
!
!    This results in a uniform distribution of values in the unit circle.
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
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer N, the number of values to compute.
!
!    Input/output, integer SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, complex ( kind = 4 ) C(N), the pseudorandom complex vector.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  complex ( kind = 4 ) c(n)
  integer i
  integer, parameter :: i4_huge = 2147483647
  integer k
  real ( kind = 4 ), parameter :: r4_pi = 3.1415926E+00
  real ( kind = 4 ) r
  integer seed
  real ( kind = 4 ) theta

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'C4VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    r = sqrt (  real ( seed, kind = 4 ) * 4.656612875E-10 )

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    theta = 2.0E+00 * r4_pi * ( real ( seed, kind = 4 ) * 4.656612875E-10 )

    c(i) = r * cmplx ( cos ( theta ), sin ( theta ), kind = 4 )

  end do

  return
end
subroutine c8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! C8MAT_PRINT_SOME prints some of a C8MAT.
!
!  Discussion:
!
!    A C8MAT is a matrix of C8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns 
!    in the matrix.
!
!    Input, complex ( kind = rk ) A(M,N), the matrix.
!
!    Input, integer ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: incx = 4
  integer m
  integer n

  complex ( kind = rk ) a(m,n)
  character ( len = 20 ) ctemp(incx)
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
  complex ( kind = rk ) zero

  zero = cmplx ( 0.0D+00, 0.0D+00, kind = rk )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of INCX.
!
  do j2lo = jlo, min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i10,10x)' ) j
    end do

    write ( *, '(a,4a20)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi
!
!  Print out (up to) INCX entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == zero ) then
          ctemp(j2) = '       0.0          '
        else if ( imag ( a(i,j) ) == 0.0D+00 ) then
          write ( ctemp(j2), '(g10.3,10x)' ) real ( a(i,j), kind = rk )
        else
          write ( ctemp(j2), '(2g10.3)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,a1,4a20)' ) i, ':', ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine c8mat_sftb ( n1, n2, y, x )

!*****************************************************************************80
!
!! C8MAT_SFTB computes a "slow" backward Fourier transform of a C8MAT.
!
!  Discussion:
!
!    SFTF and SFTB are inverses of each other.  If we begin with data
!    X and apply SFTF to get Y, and then apply SFTB to Y,
!    we should get back the original X.
!
!    This routine is provided for illustration and testing.  It is inefficient
!    relative to optimized routines that use fast Fourier techniques.
!
!    For 0 <= I1 <= N1 - 1, 
!        0 <= I2 <= N2 - 1,
!
!      X(I1,I2) = Sum ( 0 <= K2 <= N2 - 1 ) Sum ( 0 <= K1 <= N1 - 1 ) 
!        Y(K1,K2) * exp ( 2 pi i I1 K1 / N1 ) * exp ( 2 pi i I2 K2 / N2 )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, the number of rows and columns of data.
!
!    Input, complex ( kind = rk ) Y(0:N1-1,0:N2-1), the Fourier coefficients.
!
!    Output, complex ( kind = rk ) X(0:N1-1,0:N2-1), the data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n1
  integer n2

  complex ( kind = rk ) cs1
  complex ( kind = rk ) cs2
  integer i1
  integer i2
  integer j1
  integer j2
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) theta1
  real ( kind = rk ) theta2
  complex ( kind = rk ) x(0:n1-1,0:n2-1)
  complex ( kind = rk ) y(0:n1-1,0:n2-1)

  x(0:n1-1,0:n2-1) = cmplx ( 0.0D+00, 0.0D+00, kind = rk )

  do i2 = 0, n2 - 1
    do j2 = 0, n2 - 1
      theta2 = 2.0D+00 * r8_pi * real ( i2 * j2, kind = rk ) / real ( n2, kind = rk )
      cs2 = cmplx ( cos ( theta2 ), - sin ( theta2 ), kind = rk )
      do i1 = 0, n1 - 1
        do j1 = 0, n1 - 1
          theta1 = 2.0D+00 * r8_pi * real ( i1 * j1, kind = rk ) &
            / real ( n1, kind = rk )
          cs1 = cmplx ( cos ( theta1 ), - sin ( theta1 ), kind = rk )
          x(i1,i2) = x(i1,i2) + y(j1,j2) * cs1 * cs2
        end do
      end do
    end do
  end do

  x(0:n1-1,0:n2-1) = x(0:n1-1,0:n2-1) / real ( n1 * n2, kind = rk )

  return
end
subroutine c8mat_sftf ( n1, n2, x, y )

!*****************************************************************************80
!
!! C8MAT_SFTF computes a "slow" forward Fourier transform of a C8MAT.
!
!  Discussion:
!
!    SFTF and SFTB are inverses of each other.  If we begin with data
!    X and apply SFTF to get Y, and then apply SFTB to Y, 
!    we should get back the original X.
!
!    This routine is provided for illustration and testing.  It is inefficient
!    relative to optimized routines that use fast Fourier techniques.
!
!    For 0 <= I1 <= N1 - 1, 
!        0 <= I2 <= N2 - 1,
!
!      Y(I1,I2) = Sum ( 0 <= K2 <= N2 - 1 ) Sum ( 0 <= K1 <= N1 - 1 ) 
!        X(K1,K2) * exp ( - 2 pi i I1 K1 / N1 ) * exp ( - 2 pi i I2 K2 / N2 )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, the number of rows and columns of data.
!
!    Input, complex ( kind = rk ) X(0:N1-1,0:N2-1), the data to be transformed.
!
!    Output, complex ( kind = rk ) Y(0:N1-1,0:N2-1), the Fourier coefficients.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n1
  integer n2

  complex ( kind = rk ) cs1
  complex ( kind = rk ) cs2
  integer i1
  integer i2
  integer j1
  integer j2
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) theta1
  real ( kind = rk ) theta2
  complex ( kind = rk ) x(0:n1-1,0:n2-1)
  complex ( kind = rk ) y(0:n1-1,0:n2-1)

  y(0:n1-1,0:n2-1) = cmplx ( 0.0D+00, 0.0D+00, kind = rk )

  do i2 = 0, n2 - 1
    do j2 = 0, n2 - 1
      theta2 = - 2.0D+00 * r8_pi * real ( i2 * j2, kind = rk ) &
        / real ( n2, kind = rk )
      cs2 = cmplx ( cos ( theta2 ), - sin ( theta2 ), kind = rk )
      do i1 = 0, n1 - 1
        do j1 = 0, n1 - 1
          theta1 = - 2.0D+00 * r8_pi * real ( i1 * j1, kind = rk ) &
            / real ( n1, kind = rk )
          cs1 = cmplx ( cos ( theta1 ), - sin ( theta1 ), kind = rk )
          y(i1,i2) = y(i1,i2) + x(j1,j2) * cs1 * cs2
        end do
      end do
    end do
  end do

  return
end
subroutine c8mat_uniform_01 ( m, n, seed, c )

!*****************************************************************************80
!
!! C8MAT_UNIFORM_01 returns a unit pseudorandom C8MAT.
!
!  Discussion:
!
!    A C8MAT is a matrix of C8's.
!
!    For now, the input quantity SEED is an integer variable.
!
!    The angles should be uniformly distributed between 0 and 2 * PI,
!    the square roots of the radius uniformly distributed between 0 and 1.
!
!    This results in a uniform distribution of values in the unit circle.
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
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns 
!    in the matrix.
!
!    Input/output, integer SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, complex ( kind = rk ) C(M,N), the pseudorandom complex matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  complex ( kind = rk ) c(m,n)
  integer i
  integer, parameter :: i4_huge = 2147483647
  integer j
  real ( kind = rk ) r
  integer k
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  integer seed
  real ( kind = rk ) theta

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'C8MAT_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  do j = 1, n
    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge
      end if

      r = sqrt ( real ( seed, kind = rk ) * 4.656612875D-10 )

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge
      end if

      theta = 2.0D+00 * r8_pi * ( real ( seed, kind = rk ) * 4.656612875D-10 )

      c(i,j) = r * cmplx ( cos ( theta ), sin ( theta ), kind = rk )

    end do

  end do

  return
end
subroutine c8vec_print_part ( n, a, max_print, title )

!*****************************************************************************80
!
!! C8VEC_PRINT_PART prints "part" of a C8VEC.
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
!    22 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries of the vector.
!
!    Input, complex ( kind = rk ) A(N), the vector to be printed.
!
!    Input, integer MAX_PRINT, the maximum number of lines
!    to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  complex ( kind = rk ) a(n)
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
      write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)
    end do

  else if ( 3 <= max_print ) then

    do i = 1, max_print - 2
      write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)
    end do
    write ( *, '(a)' ) '  ........  ..............  ..............'
    i = n
    write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)

  else

    do i = 1, max_print - 1
      write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)
    end do
    i = max_print
    write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6,2x,a)' ) i, ':', a(i), &
      '...more entries...'

  end if

  return
end
subroutine c8vec_sftb ( n, y, x )

!*****************************************************************************80
!
!! C8VEC_SFTB computes a "slow" backward Fourier transform of a C8VEC.
!
!  Discussion:
!
!    SFTF and SFTB are inverses of each other.  If we begin with data
!    X and apply SFTF to get Y, and then apply SFTB to Y,
!    we should get back the original X.
!
!    This routine is provided for illustration and testing.  It is inefficient
!    relative to optimized routines that use fast Fourier techniques.
!
!    For 0 <= I <= N - 1
!
!      X(I) = 1/N * Sum ( 0 <= J <= N - 1 ) Y(J) * exp ( 2 pi i I J / N )
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
!    Input, integer N, the number of data values.
!
!    Input, complex ( kind = rk ) Y(0:N-1), the Fourier coefficients.
!
!    Output, complex ( kind = rk ) X(0:N-1), the data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer j
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) theta
  complex ( kind = rk ) x(0:n-1)
  complex ( kind = rk ) y(0:n-1)

  x(0:n-1) = cmplx ( 0.0D+00, 0.0D+00, kind = rk )

  do i = 0, n - 1
    do j = 0, n - 1
      theta = - 2.0D+00 * r8_pi * real ( i * j, kind = rk ) / real ( n, kind = rk )
      x(i) = x(i) + y(j) * cmplx ( cos ( theta ), sin ( theta ), kind = rk )
    end do
  end do

  x(0:n-1) = x(0:n-1) / real ( n, kind = rk )

  return
end
subroutine c8vec_sftf ( n, x, y )

!*****************************************************************************80
!
!! C8VEC_SFTF computes a "slow" forward Fourier transform of a C8VEC.
!
!  Discussion:
!
!    SFTF and SFTB are inverses of each other.  If we begin with data
!    X and apply SFTF to get Y, and then apply SFTB to Y, 
!    we should get back the original X.
!
!    This routine is provided for illustration and testing.  It is inefficient
!    relative to optimized routines that use fast Fourier techniques.
!
!    For 0 <= I <= N - 1
!
!      Y(I) = Sum ( 0 <= J <= N - 1 ) X(J) * exp ( - 2 pi i I J / N )
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
!    Input, integer N, the number of data values.
!
!    Input, complex ( kind = rk ) X(0:N-1), the data to be transformed.
!
!    Output, complex ( kind = rk ) Y(0:N-1), the Fourier coefficients.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer j
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) theta
  complex ( kind = rk ) x(0:n-1)
  complex ( kind = rk ) y(0:n-1)

  y(0:n-1) = cmplx ( 0.0D+00, 0.0D+00, kind = rk )

  do i = 0, n - 1
    do j = 0, n - 1
      theta = - 2.0D+00 * r8_pi * real ( i * j, kind = rk ) / real ( n, kind = rk )
      y(i) = y(i) + x(j) * cmplx ( cos ( theta ), - sin ( theta ), kind = rk )
    end do
  end do

  return
end
subroutine c8vec_uniform_01 ( n, seed, c )

!*****************************************************************************80
!
!! C8VEC_UNIFORM_01 returns a unit pseudorandom C8VEC.
!
!  Discussion:
!
!    A C8VEC is a vector of C8's.
!
!    For now, the input quantity SEED is an integer variable.
!
!    The angles should be uniformly distributed between 0 and 2 * PI,
!    the square roots of the radius uniformly distributed between 0 and 1.
!
!    This results in a uniform distribution of values in the unit circle.
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
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer N, the number of values to compute.
!
!    Input/output, integer SEED, the "seed" value,
!    which should NOT be 0.
!    On output, SEED has been updated.
!
!    Output, complex ( kind = rk ) C(N), the pseudorandom complex vector.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  complex ( kind = rk ) c(n)
  integer i
  integer, parameter :: i4_huge = 2147483647
  real ( kind = rk ) r
  integer k
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  integer seed
  real ( kind = rk ) theta

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'C8VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    r = sqrt ( real ( seed, kind = rk ) * 4.656612875D-10 )

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    theta = 2.0D+00 * r8_pi * ( real ( seed, kind = rk ) * 4.656612875D-10 )

    c(i) = r * cmplx ( cos ( theta ), sin ( theta ), kind = rk )

  end do

  return
end
function i4_modp ( i, j )

!*****************************************************************************80
!
!! I4_MODP returns the nonnegative remainder of I4 division.
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
!    An I4 is an integer value.
!
!  Example:
!
!        I     J     MOD I4_MODP    Factorization
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
  integer value

  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MODP - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal divisor J = ', j
    stop 1
  end if

  value = mod ( i, j )

  if ( value < 0 ) then
    value = value + abs ( j )
  end if

  i4_modp = value

  return
end
function i4_wrap ( ival, ilo, ihi )

!*****************************************************************************80
!
!! I4_WRAP forces an I4 to lie between given limits by wrapping.
!
!  Discussion:
!
!    An I4 is an integer value.
!
!    There appears to be a bug in the GFORTRAN compiler which can lead to
!    erroneous results when the first argument of I4_WRAP is an expression.
!    In particular:
!
!    do i = 1, 3
!      if ( test ) then
!        i4 = i4_wrap ( i + 1, 1, 3 )
!      end if
!    end do
!
!    was, when I = 3, returning I4 = 3.  So I had to replace this with
!
!    do i = 1, 3
!      if ( test ) then
!        i4 = i + 1
!        i4 = i4_wrap ( i4, 1, 3 )
!      end if
!    end do
!
!  Example:
!
!    ILO = 4, IHI = 8
!
!    I  Value
!
!    -2     8
!    -1     4
!     0     5
!     1     6
!     2     7
!     3     8
!     4     4
!     5     5
!     6     6
!     7     7
!     8     8
!     9     4
!    10     5
!    11     6
!    12     7
!    13     8
!    14     4
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer IVAL, a value.
!
!    Input, integer ILO, IHI, the desired bounds.
!
!    Output, integer I4_WRAP, a "wrapped" version of the value.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i4_modp
  integer i4_wrap
  integer ihi
  integer ilo
  integer ival
  integer jhi
  integer jlo
  integer value
  integer wide

  jlo = min ( ilo, ihi )
  jhi = max ( ilo, ihi )

  wide = jhi - jlo + 1

  if ( wide == 1 ) then
    value = jlo
  else
    value = jlo + i4_modp ( ival - jlo, wide )
  end if

  i4_wrap = value

  return
end
subroutine r4vec_permute_cyclic ( n, k, a )

!*****************************************************************************80
!
!! R4VEC_PERMUTE_CYCLIC performs a cyclic permutation of an R4VEC.
!
!  Discussion:
!
!    An R4VEC is a vector of R4's.
!
!    For 0 <= K < N, this function cyclically permutes the input vector
!    to have the form
!
!     ( A(K+1), A(K+2), ..., A(N), A(1), ..., A(K) )
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of objects.
!
!    Input, integer K, the increment used.
!
!    Input/output, real ( kind = 4 ) A(N), the array to be permuted.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = 4 ) a(n)
  real ( kind = 4 ) b(n)
  integer i
  integer i4_wrap
  integer ipk
  integer k

  do i = 1, n
    ipk = i4_wrap ( i + k, 1, n )
    b(i) = a(ipk)
  end do

  a(1:n) = b(1:n)

  return
end
subroutine r4vec_print_part ( n, a, max_print, title )

!*****************************************************************************80
!
!! R4VEC_PRINT_PART prints "part" of an R4VEC.
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
!    09 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries of the vector.
!
!    Input, real ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, integer MAX_PRINT, the maximum number of lines
!    to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = 4 ) a(n)
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
    write ( *, '(a)' ) '  ........  ..............'
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
subroutine r4vec_sct ( n, x, y )

!*****************************************************************************80
!
!! R4VEC_SCT computes a "slow" cosine transform of an R4VEC.
!
!  Discussion:
!
!    This routine is provided for illustration and testing.  It is inefficient
!    relative to optimized routines that use fast Fourier techniques.
!
!      Y(1) = Sum ( 1 <= J <= N ) X(J)
!
!      For 2 <= I <= N-1:
!
!        Y(I) = 2 * Sum ( 1 <= J <= N ) X(J)
!          * cos ( PI * ( I - 1 ) * ( J - 1 ) / ( N - 1 ) )
!
!      Y(N) = Sum ( X(1:N:2) ) - Sum ( X(2:N:2) )
!
!    Applying the routine twice in succession should yield the original data,
!    multiplied by 2 * ( N + 1 ).  This is a good check for correctness
!    and accuracy.
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
!    Input, integer N, the number of data values.
!
!    Input, real ( kind = 4 ) X(N), the data sequence.
!
!    Output, real ( kind = 4 ) Y(N), the transformed data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = 4 ) angle
  integer i
  integer j
  real ( kind = 4 ), parameter :: r4_pi = 3.1415926E+00
  real ( kind = 4 ) x(n)
  real ( kind = 4 ) y(n)

  do i = 1, n

    y(i) = x(1) / 2.0E+00

    do j = 2, n - 1
      angle = r4_pi &
        * real ( mod ( ( i - 1 ) * ( j - 1 ), 2 * ( n - 1 ) ), kind = 4 ) &
        / real ( n - 1, kind = 4 )
      y(i) = y(i) + x(j) * cos ( angle )
    end do

    j = n

    angle = r4_pi &
      * real ( mod ( ( i - 1 ) * ( j - 1 ), 2 * ( n - 1 ) ), kind = 4 ) &
      / real ( n - 1, kind = 4 )

    y(i) = y(i) + x(j) * cos ( angle ) / 2.0E+00

  end do

  y(1:n) = 2.0E+00 * y(1:n) &
    * sqrt ( real ( n, kind = 4 ) / real ( n - 1, kind = 4 ) )

  return
end
subroutine r4vec_sftb ( n, azero, a, b, r )

!*****************************************************************************80
!
!! R4VEC_SFTB computes a "slow" backward Fourier transform of an R4VEC.
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
!    10 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data values.
!
!    Input, real ( kind = 4 ) AZERO, the constant Fourier coefficient.
!
!    Input, real ( kind = 4 ) A(N/2), B(N/2), the Fourier coefficients.
!
!    Output, real ( kind = 4 ) R(N), the reconstructed data sequence.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = 4 ) a(n/2)
  real ( kind = 4 ) azero
  real ( kind = 4 ) b(n/2)
  integer i
  integer k
  real ( kind = 4 ) r(n)
  real ( kind = 4 ), parameter :: r4_pi = 3.1415926E+00
  real ( kind = 4 ) theta

  r(1:n) = azero
  do i = 1, n
    do k = 1, n / 2
      theta = real ( k * ( i - 1 ) * 2, kind = 4 ) * r4_pi &
        / real ( n, kind = 4 )
      r(i) = r(i) + a(k) * cos ( theta ) + b(k) * sin ( theta )
    end do
  end do

  return
end
subroutine r4vec_sftf ( n, r, azero, a, b )

!*****************************************************************************80
!
!! R4VEC_SFTF computes a "slow" forward Fourier transform of an R4VEC.
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
!    Input, real ( kind = 4 ) R(N), the data to be transformed.
!
!    Output, real ( kind = 4 ) AZERO, = sum ( 1 <= I <= N ) R(I) / N.
!
!    Output, real ( kind = 4 ) A(N/2), B(N/2), the Fourier coefficients.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = 4 ) a(1:n/2)
  real ( kind = 4 ) azero
  real ( kind = 4 ) b(1:n/2)
  integer i
  integer j
  real ( kind = 4 ), parameter :: r4_pi = 3.1415926E+00
  real ( kind = 4 ) r(n)
  real ( kind = 4 ) theta

  azero = sum ( r(1:n) ) / real ( n, kind = 4 )

  do i = 1, n / 2

    a(i) = 0.0E+00
    b(i) = 0.0E+00

    do j = 1, n
      theta = real ( 2 * i * ( j - 1 ), kind = 4 ) * r4_pi &
        / real ( n, kind = 4 )
      a(i) = a(i) + r(j) * cos ( theta )
      b(i) = b(i) + r(j) * sin ( theta )
    end do

    a(i) = a(i) / real ( n, kind = 4 )
    b(i) = b(i) / real ( n, kind = 4 )

    if ( mod ( n, 2 ) == 1 .or. i /= ( n / 2 ) ) then
      a(i) = 2.0E+00 * a(i)
      b(i) = 2.0E+00 * b(i)
    end if

  end do

  return
end
subroutine r4vec_sht ( n, a, b  )

!*****************************************************************************80
!
!! R4VEC_SHT computes a "slow" Hartley transform of an R4VEC.
!
!  Discussion:
!
!    The discrete Hartley transform B of a set of data A is
!
!      B(I) = 1/sqrt(N) * Sum ( 0 <= J <= N-1 ) A(J) * CAS(2*PI*I*J/N)
!
!    Here, the data and coefficients are indexed from 0 to N-1.
!
!    With the above normalization factor of 1/sqrt(N), the Hartley
!    transform is its own inverse.
!
!    This routine is provided for illustration and testing.  It is inefficient
!    relative to optimized routines.
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
!  Reference:
!
!    Ralph Hartley,
!    A More Symmetrical Fourier Analysis Applied to Transmission Problems,
!    Proceedings of the Institute of Radio Engineers,
!    Volume 30, pages 144-150, 1942.
!
!  Parameters:
!
!    Input, integer N, the number of data values.
!
!    Input, real ( kind = 4 ) A(0:N-1), the data to be transformed.
!
!    Output, real ( kind = 4 ) B(0:N-1), the transformed data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = 4 ) a(0:n-1)
  real ( kind = 4 ) b(0:n-1)
  integer i
  integer j
  real ( kind = 4 ), parameter :: r4_pi = 3.1415926E+00
  real ( kind = 4 ) theta

  b(0:n-1) = 0.0E+00

  do i = 0, n - 1
    do j = 0, n - 1
      theta = 2.0E+00 * r4_pi * real ( mod ( i * j, n ), kind = 4 ) &
        / real ( n, kind = 4 )
      b(i) = b(i) + a(j) * ( cos ( theta ) + sin ( theta ) )
    end do
  end do

  b(0:n-1) = b(0:n-1) / sqrt ( real ( n, kind = 4 ) )

  return
end
subroutine r4vec_sqctb ( n, x, y )

!*****************************************************************************80
!
!! R4VEC_SQCTB computes a "slow" quarter cosine transform backward of an R4VEC.
!
!  Discussion:
!
!    This routine is provided for illustration and testing.  It is inefficient
!    relative to optimized routines that use fast Fourier techniques.
!
!    For 0 <= I <= N-1,
!
!      Y(I) = X(0) + 2 Sum ( 1 <= J <= N-1 ) X(J) * cos ( PI * J * (I+1/2) / N )
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
!  Reference:
!
!    William Briggs, Van Emden Henson,
!    The Discrete Fourier Transform,
!    SIAM, 1995,
!    LC: QA403.5 B75
!
!  Parameters:
!
!    Input, integer N, the number of data values.
!
!    Input, real ( kind = 4 ) X(0:N-1), the data sequence.
!
!    Output, real ( kind = 4 ) Y(0:N-1), the transformed data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer j
  real ( kind = 4 ), parameter :: r4_pi = 3.1415926E+00
  real ( kind = 4 ) theta
  real ( kind = 4 ) x(0:n-1)
  real ( kind = 4 ) y(0:n-1)

  y(0:n-1) = x(0)

  do i = 0, n - 1
    do j = 1, n - 1

      theta = 0.5E+00 * r4_pi * real ( j * ( 2 * i + 1 ), kind = 4 ) &
        / real ( n, kind = 4 )
      y(i) = y(i) + 2.0E+00 * x(j) * cos ( theta  )

    end do

  end do

  return
end
subroutine r4vec_sqctf ( n, x, y )

!*****************************************************************************80
!
!! R4VEC_SQCTF computes a "slow" quarter cosine transform forward of an R4VEC.
!
!  Discussion:
!
!    This routine is provided for illustration and testing.  It is inefficient
!    relative to optimized routines that use fast Fourier techniques.
!
!    For 0 <= I <= N-1,
!
!      Y(I) = (1/N) Sum ( 0 <= J <= N-1 ) X(J) * cos ( PI * I * (J+1/2) / N )
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
!  Reference:
!
!    William Briggs, Van Emden Henson,
!    The Discrete Fourier Transform,
!    SIAM, 1995,
!    QA403.5 B75
!
!  Parameters:
!
!    Input, integer N, the number of data values.
!
!    Input, real ( kind = 4 ) X(0:N-1), the data sequence.
!
!    Output, real ( kind = 4 ) Y(0:N-1), the transformed data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer j
  real ( kind = 4 ), parameter :: r4_pi = 3.1415926E+00
  real ( kind = 4 ) theta
  real ( kind = 4 ) x(0:n-1)
  real ( kind = 4 ) y(0:n-1)

  y(0:n-1) = 0.0E+00

  do i = 0, n - 1
    do j = 0, n - 1
      theta = 0.5E+00 * r4_pi * real ( i * ( 2 * j + 1 ), kind = 4 ) &
        / real ( n, kind = 4 )
      y(i) = y(i) + x(j) * cos ( theta  )
    end do
  end do

  y(0:n-1) = y(0:n-1) / real ( n, kind = 4 )

  return
end
subroutine r4vec_sqstb ( n, x, y )

!*****************************************************************************80
!
!! R4VEC_SQSTB computes a "slow" quarter sine transform backward of an R4VEC.
!
!  Discussion:
!
!    This routine is provided for illustration and testing.  It is inefficient
!    relative to optimized routines that use fast Fourier techniques.
!
!    For 0 <= I <= N-1,
!
!      Y(I) = -2 Sum ( 1 <= J <= N-1 ) X(J) * sin ( PI * J * (I+1/2) / N )
!             - X(N) * cos ( pi * I )
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
!  Reference:
!
!    William Briggs, Van Emden Henson,
!    The Discrete Fourier Transform,
!    SIAM, 1995,
!    QA403.5 B75
!
!  Parameters:
!
!    Input, integer N, the number of data values.
!
!    Input, real ( kind = 4 ) X(N), the data sequence.
!
!    Output, real ( kind = 4 ) Y(0:N-1), the transformed data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer j
  real ( kind = 4 ), parameter :: r4_pi = 3.1415926E+00
  real ( kind = 4 ) theta
  real ( kind = 4 ) x(1:n)
  real ( kind = 4 ) y(0:n-1)

  y(0:n-1) = 0.0E+00

  do i = 0, n - 1
    do j = 1, n - 1

      theta = 0.5E+00 * r4_pi * real ( j * ( 2 * i + 1 ), kind = 4 ) &
        / real ( n, kind = 4 )
      y(i) = y(i) - 2.0E+00 * x(j) * sin ( theta  )

    end do

    theta = r4_pi * real ( i, kind = 4 )
    y(i) = y(i) - x(n) * cos ( theta )

  end do

  return
end
subroutine r4vec_sqstf ( n, x, y )

!*****************************************************************************80
!
!! R4VEC_SQSTF computes a "slow" quarter sine transform forward of an R4VEC.
!
!  Discussion:
!
!    This routine is provided for illustration and testing.  It is inefficient
!    relative to optimized routines that use fast Fourier techniques.
!
!    For 1 <= I <= N,
!
!      Y(I) = -(1/N) Sum ( 0 <= J <= N-1 ) X(J) * sin ( PI * I * (J+1/2) / N )
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
!  Reference:
!
!    William Briggs, Van Emden Henson,
!    The Discrete Fourier Transform,
!    SIAM, 1995,
!    QA403.5 B75
!
!  Parameters:
!
!    Input, integer N, the number of data values.
!
!    Input, real ( kind = 4 ) X(0:N-1), the data sequence.
!
!    Output, real ( kind = 4 ) Y(N), the transformed data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer j
  real ( kind = 4 ), parameter :: r4_pi = 3.1415926E+00
  real ( kind = 4 ) theta
  real ( kind = 4 ) x(0:n-1)
  real ( kind = 4 ) y(n)

  y(1:n) = 0.0E+00

  do i = 1, n
    do j = 0, n - 1
      theta = 0.5E+00 * r4_pi * real ( i * ( 2 * j + 1 ), kind = 4 ) &
        / real ( n, kind = 4 )
      y(i) = y(i) + x(j) * sin ( theta  )
    end do
  end do

  y(1:n) = - y(1:n) / real ( n, kind = 4 )

  return
end
subroutine r4vec_sst ( n, x, y )

!*****************************************************************************80
!
!! R4VEC_SST computes a "slow" sine transform of an R4VEC.
!
!  Discussion:
!
!    This routine is provided for illustration and testing.  It is inefficient
!    relative to optimized routines that use fast Fourier techniques.
!
!    For 1 <= I <= N,
!
!      Y(I) = Sum ( 1 <= J <= N ) X(J) * sin ( PI * I * J / ( N + 1 ) )
!
!    Applying the routine twice in succession should yield the original data,
!    multiplied by N / 2.  This is a good check for correctness and accuracy.
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
!    Input, integer N, the number of data values.
!
!    Input, real ( kind = 4 ) X(N), the data sequence.
!
!    Output, real ( kind = 4 ) Y(N), the transformed data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  real ( kind = 4 ), parameter :: r4_pi = 3.1415926E+00
  real ( kind = 4 ) theta(n)
  real ( kind = 4 ) x(n)
  real ( kind = 4 ) y(n)

  do i = 1, n
    theta(i) = r4_pi * real ( i, kind = 4 ) / real ( n + 1, kind = 4 )
  end do

  y(1:n) = 0.0E+00

  do i = 1, n
    y(1:n) = y(1:n) + 2.0E+00 * x(i) * sin ( real ( i, kind = 4 ) * theta(1:n) )
  end do

  return
end
subroutine r4vec_uniform_ab ( n, a, b, seed, r )

!*****************************************************************************80
!
!! R4VEC_UNIFORM_AB returns a scaled pseudorandom R4VEC.
!
!  Discussion:
!
!    An R4VEC is a vector of R4's.
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
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer M, the number of entries in the vector.
!
!    Input, real ( kind = 4 ) A, B, the lower and upper limits.
!
!    Input/output, integer SEED, the "seed" value, 
!    which should NOT be 0.
!    On output, SEED has been updated.
!
!    Output, real ( kind = 4 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = 4 ) a
  real ( kind = 4 ) b
  integer i
  integer, parameter :: i4_huge = 2147483647
  integer k
  integer seed
  real ( kind = 4 ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4VEC_UNIFORM_AB - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    r(i) = a + ( b - a ) * real ( seed, kind = 4 ) * 4.656612875E-10

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
      write ( *, '(2x,i8,2x,g14.6)' ) i, a(i)
    end do

  else if ( 3 <= max_print ) then

    do i = 1, max_print-2
      write ( *, '(2x,i8,2x,g14.6)' ) i, a(i)
    end do
    write ( *, '(a)' ) '  ........  ..............'
    i = n
    write ( *, '(2x,i8,2x,g14.6)' ) i, a(i)

  else

    do i = 1, max_print - 1
      write ( *, '(2x,i8,2x,g14.6)' ) i, a(i)
    end do
    i = max_print
    write ( *, '(2x,i8,2x,g14.6,2x,a)' ) i, a(i), '...more entries...'

  end if

  return
end
subroutine r8vec_permute_cyclic ( n, k, a )

!*****************************************************************************80
!
!! R8VEC_PERMUTE_CYCLIC performs a cyclic permutation of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    For 0 <= K < N, this function cyclically permutes the input vector
!    to have the form
!
!     ( A(K+1), A(K+2), ..., A(N), A(1), ..., A(K) )
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of objects.
!
!    Input, integer K, the increment used.
!
!    Input/output, real ( kind = rk ) A(N), the array to be permuted.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) b(n)
  integer i
  integer i4_wrap
  integer ipk
  integer k

  do i = 1, n
    ipk = i4_wrap ( i + k, 1, n )
    b(i) = a(ipk)
  end do

  a(1:n) = b(1:n)

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
subroutine r8vec_sct ( n, x, y )

!*****************************************************************************80
!
!! R8VEC_SCT computes a "slow" cosine transform of an R8VEC.
!
!  Discussion:
!
!    This routine is provided for illustration and testing.  It is inefficient
!    relative to optimized routines that use fast Fourier techniques.
!
!      Y(1) = Sum ( 1 <= J <= N ) X(J)
!
!      For 2 <= I <= N-1:
!
!        Y(I) = 2 * Sum ( 1 <= J <= N ) X(J)
!          * cos ( PI * ( I - 1 ) * ( J - 1 ) / ( N - 1 ) )
!
!      Y(N) = Sum ( X(1:N:2) ) - Sum ( X(2:N:2) )
!
!    Applying the routine twice in succession should yield the original data,
!    multiplied by 2 * ( N + 1 ).  This is a good check for correctness
!    and accuracy.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data values.
!
!    Input, real ( kind = rk ) X(N), the data sequence.
!
!    Output, real ( kind = rk ) Y(N), the transformed data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) angle
  integer i
  integer j
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)

  do i = 1, n

    y(i) = x(1) / 2.0D+00

    do j = 2, n - 1
      angle = r8_pi &
        * real ( mod ( ( i - 1 ) * ( j - 1 ), 2 * ( n - 1 ) ), kind = rk ) &
        / real ( n - 1, kind = rk )
      y(i) = y(i) + x(j) * cos ( angle )
    end do

    j = n

    angle = r8_pi &
      * real ( mod ( ( i - 1 ) * ( j - 1 ), 2 * ( n - 1 ) ), kind = rk ) &
      / real ( n - 1, kind = rk )

    y(i) = y(i) + x(j) * cos ( angle ) / 2.0D+00

  end do

  y(1:n) = 2.0D+00 * y(1:n) &
    * sqrt ( real ( n, kind = rk ) / real ( n - 1, kind = rk ) )

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
subroutine r8vec_sht ( n, a, b  )

!*****************************************************************************80
!
!! R8VEC_SHT computes a "slow" Hartley transform of an R8VEC.
!
!  Discussion:
!
!    The discrete Hartley transform B of a set of data A is
!
!      B(I) = 1/sqrt(N) * Sum ( 0 <= J <= N-1 ) A(J) * CAS(2*PI*I*J/N)
!
!    Here, the data and coefficients are indexed from 0 to N-1.
!
!    With the above normalization factor of 1/sqrt(N), the Hartley
!    transform is its own inverse.
!
!    This routine is provided for illustration and testing.  It is inefficient
!    relative to optimized routines.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Ralph Hartley,
!    A More Symmetrical Fourier Analysis Applied to Transmission Problems,
!    Proceedings of the Institute of Radio Engineers,
!    Volume 30, pages 144-150, 1942.
!
!  Parameters:
!
!    Input, integer N, the number of data values.
!
!    Input, real ( kind = rk ) A(0:N-1), the data to be transformed.
!
!    Output, real ( kind = rk ) B(0:N-1), the transformed data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(0:n-1)
  real ( kind = rk ) b(0:n-1)
  integer i
  integer j
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) theta

  b(0:n-1) = 0.0D+00

  do i = 0, n - 1
    do j = 0, n - 1
      theta = 2.0D+00 * r8_pi * real ( mod ( i * j, n ), kind = rk ) &
        / real ( n, kind = rk )
      b(i) = b(i) + a(j) * ( cos ( theta ) + sin ( theta ) )
    end do
  end do

  b(0:n-1) = b(0:n-1) / sqrt ( real ( n, kind = rk ) )

  return
end
subroutine r8vec_sqctb ( n, x, y )

!*****************************************************************************80
!
!! R8VEC_SQCTB computes a "slow" quarter cosine transform backward of an R8VEC.
!
!  Discussion:
!
!    This routine is provided for illustration and testing.  It is inefficient
!    relative to optimized routines that use fast Fourier techniques.
!
!    For 0 <= I <= N-1,
!
!      Y(I) = X(0) + 2 Sum ( 1 <= J <= N-1 ) X(J) * cos ( PI * J * (I+1/2) / N )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Briggs, Van Emden Henson,
!    The Discrete Fourier Transform,
!    SIAM, 1995,
!    LC: QA403.5 B75
!
!  Parameters:
!
!    Input, integer N, the number of data values.
!
!    Input, real ( kind = rk ) X(0:N-1), the data sequence.
!
!    Output, real ( kind = rk ) Y(0:N-1), the transformed data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer j
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) theta
  real ( kind = rk ) x(0:n-1)
  real ( kind = rk ) y(0:n-1)

  y(0:n-1) = x(0)

  do i = 0, n - 1
    do j = 1, n - 1

      theta = 0.5D+00 * r8_pi * real ( j * ( 2 * i + 1 ), kind = rk ) &
        / real ( n, kind = rk )
      y(i) = y(i) + 2.0D+00 * x(j) * cos ( theta  )

    end do

  end do

  return
end
subroutine r8vec_sqctf ( n, x, y )

!*****************************************************************************80
!
!! R8VEC_SQCTF computes a "slow" quarter cosine transform forward of an R8VEC.
!
!  Discussion:
!
!    This routine is provided for illustration and testing.  It is inefficient
!    relative to optimized routines that use fast Fourier techniques.
!
!    For 0 <= I <= N-1,
!
!      Y(I) = (1/N) Sum ( 0 <= J <= N-1 ) X(J) * cos ( PI * I * (J+1/2) / N )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Briggs, Van Emden Henson,
!    The Discrete Fourier Transform,
!    SIAM, 1995,
!    QA403.5 B75
!
!  Parameters:
!
!    Input, integer N, the number of data values.
!
!    Input, real ( kind = rk ) X(0:N-1), the data sequence.
!
!    Output, real ( kind = rk ) Y(0:N-1), the transformed data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer j
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) theta
  real ( kind = rk ) x(0:n-1)
  real ( kind = rk ) y(0:n-1)

  y(0:n-1) = 0.0D+00

  do i = 0, n - 1
    do j = 0, n - 1
      theta = 0.5D+00 * r8_pi * real ( i * ( 2 * j + 1 ), kind = rk ) &
        / real ( n, kind = rk )
      y(i) = y(i) + x(j) * cos ( theta  )
    end do
  end do

  y(0:n-1) = y(0:n-1) / real ( n, kind = rk )

  return
end
subroutine r8vec_sqstb ( n, x, y )

!*****************************************************************************80
!
!! R8VEC_SQSTB computes a "slow" quarter sine transform backward of an R8VEC.
!
!  Discussion:
!
!    This routine is provided for illustration and testing.  It is inefficient
!    relative to optimized routines that use fast Fourier techniques.
!
!    For 0 <= I <= N-1,
!
!      Y(I) = -2 Sum ( 1 <= J <= N-1 ) X(J) * sin ( PI * J * (I+1/2) / N )
!             - X(N) * cos ( pi * I )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Briggs, Van Emden Henson,
!    The Discrete Fourier Transform,
!    SIAM, 1995,
!    QA403.5 B75
!
!  Parameters:
!
!    Input, integer N, the number of data values.
!
!    Input, real ( kind = rk ) X(N), the data sequence.
!
!    Output, real ( kind = rk ) Y(0:N-1), the transformed data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer j
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) theta
  real ( kind = rk ) x(1:n)
  real ( kind = rk ) y(0:n-1)

  y(0:n-1) = 0.0D+00

  do i = 0, n - 1
    do j = 1, n - 1

      theta = 0.5D+00 * r8_pi * real ( j * ( 2 * i + 1 ), kind = rk ) &
        / real ( n, kind = rk )
      y(i) = y(i) - 2.0D+00 * x(j) * sin ( theta  )

    end do

    theta = r8_pi * real ( i, kind = rk )
    y(i) = y(i) - x(n) * cos ( theta )

  end do

  return
end
subroutine r8vec_sqstf ( n, x, y )

!*****************************************************************************80
!
!! R8VEC_SQSTF computes a "slow" quarter sine transform forward of an R8VEC.
!
!  Discussion:
!
!    This routine is provided for illustration and testing.  It is inefficient
!    relative to optimized routines that use fast Fourier techniques.
!
!    For 1 <= I <= N,
!
!      Y(I) = -(1/N) Sum ( 0 <= J <= N-1 ) X(J) * sin ( PI * I * (J+1/2) / N )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Briggs, Van Emden Henson,
!    The Discrete Fourier Transform,
!    SIAM, 1995,
!    QA403.5 B75
!
!  Parameters:
!
!    Input, integer N, the number of data values.
!
!    Input, real ( kind = rk ) X(0:N-1), the data sequence.
!
!    Output, real ( kind = rk ) Y(N), the transformed data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer j
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) theta
  real ( kind = rk ) x(0:n-1)
  real ( kind = rk ) y(n)

  y(1:n) = 0.0D+00

  do i = 1, n
    do j = 0, n - 1
      theta = 0.5D+00 * r8_pi * real ( i * ( 2 * j + 1 ), kind = rk ) &
        / real ( n, kind = rk )
      y(i) = y(i) + x(j) * sin ( theta  )
    end do
  end do

  y(1:n) = - y(1:n) / real ( n, kind = rk )

  return
end
subroutine r8vec_sst ( n, x, y )

!*****************************************************************************80
!
!! R8VEC_SST computes a "slow" sine transform of an R8VEC.
!
!  Discussion:
!
!    This routine is provided for illustration and testing.  It is inefficient
!    relative to optimized routines that use fast Fourier techniques.
!
!    For 1 <= I <= N,
!
!      Y(I) = Sum ( 1 <= J <= N ) X(J) * sin ( PI * I * J / ( N + 1 ) )
!
!    Applying the routine twice in succession should yield the original data,
!    multiplied by N / 2.  This is a good check for correctness and accuracy.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data values.
!
!    Input, real ( kind = rk ) X(N), the data sequence.
!
!    Output, real ( kind = rk ) Y(N), the transformed data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) theta(n)
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)

  do i = 1, n
    theta(i) = r8_pi * real ( i, kind = rk ) / real ( n + 1, kind = rk )
  end do

  y(1:n) = 0.0D+00

  do i = 1, n
    y(1:n) = y(1:n) + 2.0D+00 * x(i) * sin ( real ( i, kind = rk ) * theta(1:n) )
  end do

  return
end
subroutine r8vec_swtb ( n, s, d, x )

!*****************************************************************************80
!
!! R8VEC_SWTB computes a "slow" backward wavelet transform of an R8VEC.
!
!  Discussion:
!
!    This function inverts the D4 Daubechies wavelet.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 June 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data values.
!
!    Input, real ( kind = rk ) S((N+1)/2), D((N+1)/2), the transformed data.
!
!    Output, real ( kind = rk ) X(N), the original data sequence.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) d((n+1)/2)
  real ( kind = rk ) d1((n+1)/2)
  real ( kind = rk ) d2((n+1)/2)
  real ( kind = rk ) s((n+1)/2)
  real ( kind = rk ) s1((n+1)/2)
  real ( kind = rk ) s2((n+1)/2)
  real ( kind = rk ) s3((n+1)/2)
  integer n2
  integer nh
  integer np1h
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(2*((n+1)/2))

  if ( mod ( n, 2 ) == 1 ) then
    n2 = n + 1
  else
    n2 = n
  end if

  np1h = ( n + 1 ) / 2
  nh = n / 2

  d1(1:np1h) = d(1:np1h) / ( ( sqrt ( 3.0D+00 ) + 1.0D+00 ) / sqrt ( 2.0D+00 ) )

  s1(1:np1h) = s(1:np1h) / ( ( sqrt ( 3.0D+00 ) - 1.0D+00 ) / sqrt ( 2.0D+00 ) )

  d2(1:np1h) = d1(1:np1h)
  call r8vec_permute_cyclic ( np1h, +1, d2 )
  s2(1:np1h) = s1(1:np1h) + d2(1:np1h)

  s3(1:np1h) = s2(1:np1h)
  call r8vec_permute_cyclic ( np1h, -1, s3 )

  y(2:n2:2) = d1(1:np1h) &
    + sqrt ( 3.0D+00 ) / 4.0D+00 * s2(1:np1h) &
    + ( sqrt ( 3.0D+00 ) - 2.0D+00 ) / 4.0D+00 * s3(1:np1h)

  y(1:n2-1:2) = s2(1:nh) - sqrt ( 3.0D+00 ) * y(2:n2:2)
 
  x(1:n) = y(1:n)

  return
end
subroutine r8vec_swtf ( n, x, s, d )

!*****************************************************************************80
!
!! R8VEC_SWTF computes a "slow" forward wavelet transform of an R8VEC.
!
!  Discussion:
!
!    This function applies the D4 Daubechies wavelet.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 June 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data values.
!
!    Input, real ( kind = rk ) X(N), the data sequence.
!
!    Output, real ( kind = rk ) S((N+1)/2), D((N+1)/2), the transformed data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) d((n+1)/2)
  real ( kind = rk ) d2((n+1)/2)
  integer n2
  integer np1h
  real ( kind = rk ) s((n+1)/2)
  real ( kind = rk ) s2((n+1)/2)
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(2*((n+1)/2))

  y(1:n) = x(1:n)

  if ( mod ( n, 2 ) == 1 ) then
    n2 = n + 1
    y(n2) = 0.0D+00
  else
    n2 = n
  end if

  np1h = ( n + 1 ) / 2

  s(1:np1h) = y(1:n2-1:2) + sqrt ( 3.0D+00 ) * y(2:n2:2)

  s2(1:np1h) = s(1:np1h)
  call r8vec_permute_cyclic ( np1h, -1, s2 )

  d(1:np1h) = y(2:n2:2) &
    - sqrt ( 3.0D+00 ) / 4.0D+00 * s(1:np1h) &
    - ( sqrt ( 3.0D+00 ) - 2.0D+00 ) / 4.0D+00 * s2(1:np1h)

  d2(1:np1h) = d(1:np1h)
  call r8vec_permute_cyclic ( np1h, +1, d2 )

  s(1:np1h) = s(1:np1h) - d2(1:np1h)

  s = ( sqrt ( 3.0D+00 ) - 1.0D+00 ) / sqrt ( 2.0 ) * s(1:np1h)
  d = ( sqrt ( 3.0D+00 ) + 1.0D+00 ) / sqrt ( 2.0 ) * d(1:np1h)

  return
end
subroutine r8vec_uniform_ab ( n, a, b, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_AB returns a scaled pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    Each dimension ranges from A to B.
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
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real ( kind = rk ) A, B, the lower and upper limits.
!
!    Input/output, integer SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = rk ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a
  real ( kind = rk ) b
  integer i
  integer, parameter :: i4_huge = 2147483647
  integer k
  integer seed
  real ( kind = rk ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_AB - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    r(i) = a + ( b - a ) * real ( seed, kind = rk ) * 4.656612875D-10

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
