function c8_normal_01 ( )

!*****************************************************************************80
!
!! c8_normal_01() returns a unit pseudonormal C8.
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
!  Parameters:
!
!    Output, complex ( kind = ck ) C8_NORMAL_01, a sample of the PDF.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )
  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  complex ( kind = ck ) c8_normal_01
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) v1
  real ( kind = rk ) v2
  real ( kind = rk ) x_c
  real ( kind = rk ) x_r

  call random_number ( harvest = v1 )
  call random_number ( harvest = v2 )

  x_r = sqrt ( - 2.0D+00 * log ( v1 ) ) * cos ( 2.0D+00 * r8_pi * v2 )
  x_c = sqrt ( - 2.0D+00 * log ( v1 ) ) * sin ( 2.0D+00 * r8_pi * v2 )

  c8_normal_01 = cmplx ( x_r, x_c, kind = ck )

  return
end
subroutine c8vec_normal_01 ( n, x )

!*****************************************************************************80
!
!! c8vec_normal_01() returns a unit pseudonormal C8VEC.
!
!  Discussion:
!
!    A C8VEC is an array of double precision complex values.
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
!    07 December 2023
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of values desired.
!
!  Output:
!
!    complex ( kind = ck ) X(N), a sample of the standard normal PDF.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) r(n)
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  complex ( kind = ck ) x(n)

  call random_number ( harvest = r )

  x = sqrt ( - 2.0D+00 * log ( r ) ) * &
    cmplx ( cos ( 2.0D+00 * r8_pi * r ), &
            sin ( 2.0D+00 * r8_pi * r ), kind = ck )

  return
end
subroutine c8vec_print ( n, a, title )

!*****************************************************************************80
!
!! c8vec_print() prints a C8VEC, with an optional title.
!
!  Discussion:
!
!    A C8VEC is a vector of C8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of components of the vector.
!
!    complex ( kind = ck ) A(N), the vector to be printed.
!
!    character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  integer n

  complex ( kind = ck ) a(n)
  integer i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,2x,2g14.6)' ) i, a(i)
  end do

  return
end
function i4_normal_ab ( a, b )

!*****************************************************************************80
!
!! i4_normal_ab() returns a scaled pseudonormal I4.
!
!  Discussion:
!
!    The normal probability distribution function (PDF) is sampled,
!    with mean A and standard deviation B.
!
!    The result is then rounded to the nearest integer.
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
!  Parameters:
!
!    Input, real ( kind = rk ) A, the mean of the PDF.
!
!    Input, real ( kind = rk ) B, the standard deviation of the PDF.
!
!    Output, integer I4_NORMAL_AB, a sample of the normal PDF.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  integer i4_normal_ab
  real ( kind = rk ) r1
  real ( kind = rk ) r2
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) x

  call random_number ( harvest = r1 )
  call random_number ( harvest = r2 )

  x = sqrt ( - 2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * r8_pi * r2 )

  i4_normal_ab = nint ( a + b * x )

  return
end
subroutine i4_normal_ab_test ( )

!*****************************************************************************80
!
!! i4_normal_ab_test() tests I4_NORMAL_AB.
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer i4_normal_ab
  real ( kind = rk ) mu
  integer r
  real ( kind = rk ) sigma

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'I4_NORMAL_AB_TEST'
  write ( *, '(a)' ) '  I4_NORMAL_AB computes integer pseudonormal values '
  write ( *, '(a)' ) '  with mean MU and standard deviation SIGMA.'

  mu = 70.0D+00
  sigma = 10.0D+00

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  MU = ', mu
  write ( *, '(a,g14.6)' ) '  SIGMA = ', sigma
  write ( *, '(a)' ) ''

  do i = 1, 10
    r = i4_normal_ab ( mu, sigma )
    write ( *, '(2x,i8,2x,i8)' ) i, r
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
!  Parameters:
!
!    Output, real ( kind = rk ) R8_NORMAL_01, a normally distributed
!    random value.
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
function r8_normal_ab ( a, b )

!*****************************************************************************80
!
!! r8_normal_ab() returns a scaled pseudonormal R8.
!
!  Discussion:
!
!    The normal probability distribution function (PDF) is sampled,
!    with mean A and standard deviation B.
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
!  Parameters:
!
!    Input, real ( kind = rk ) A, the mean of the PDF.
!
!    Input, real ( kind = rk ) B, the standard deviation of the PDF.
!
!    Output, real ( kind = rk ) R8_NORMAL_AB, a sample of the normal PDF.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) r1
  real ( kind = rk ) r2
  real ( kind = rk ) r8_normal_ab
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) x

  call random_number ( harvest = r1 )
  call random_number ( harvest = r2 )
  x = sqrt ( - 2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * r8_pi * r2 )

  r8_normal_ab = a + b * x

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
!    12 November 2010
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
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns
!    in the array.
!
!    Output, real ( kind = rk ) R(M,N), the array of pseudonormal values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) r(m,n)

  call r8vec_normal_01 ( m * n, r )

  return
end
subroutine r8mat_normal_ab ( m, n, a, b, r )

!*****************************************************************************80
!
!! r8mat_normal_ab() returns a scaled pseudonormal R8MAT.
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
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns
!    in the array.
!
!    Input, real ( kind = rk ) A, B, the mean and standard deviation.
!
!    Output, real ( kind = rk ) R(M,N), the array of pseudonormal values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) r(m,n)

  call r8vec_normal_ab ( m * n, a, b, r )

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! r8mat_print() prints an R8MAT.
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
!! r8mat_print_some() prints some of an R8MAT.
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

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)'
    return
  end if

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
subroutine r8vec_normal_01 ( n, x )

!*****************************************************************************80
!
!! r8vec_normal_01() returns a unit pseudonormal R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
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
!    18 May 2014
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of values desired.
!
!  Output:
!
!    real ( kind = rk ) X(N), a sample of the standard normal PDF.
!
!  Local:
!
!    real ( kind = rk ) R(N+1), is used to store some uniform
!    random values.  Its dimension is N+1, but really it is only needed
!    to be the smallest even number greater than or equal to N.
!
!    integer X_LO_INDEX, X_HI_INDEX, records the range
!    of entries of X that we need to compute
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer m
  real ( kind = rk ) r(n+1)
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) x(n)
  integer x_hi_index
  integer x_lo_index
!
!  Record the range of X we need to fill in.
!
  x_lo_index = 1
  x_hi_index = n
!
!  If we need just one new value, do that here to avoid null arrays.
!
  if ( x_hi_index - x_lo_index + 1 == 1 ) then

    call random_number ( harvest = r(1:2) )

    x(x_hi_index) = &
      sqrt ( - 2.0D+00 * log ( r(1) ) ) * cos ( 2.0D+00 * r8_pi * r(2) )
!
!  If we require an even number of values, that's easy.
!
  else if ( mod ( x_hi_index - x_lo_index, 2 ) == 1 ) then

    m = ( x_hi_index - x_lo_index + 1 ) / 2

    call random_number ( harvest = r(1:2*m) )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( - 2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * cos ( 2.0D+00 * r8_pi * r(2:2*m:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( - 2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * sin ( 2.0D+00 * r8_pi * r(2:2*m:2) )
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
      sqrt ( - 2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * cos ( 2.0D+00 * r8_pi * r(2:2*m-2:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( - 2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * sin ( 2.0D+00 * r8_pi * r(2:2*m-2:2) )

    x(n) = sqrt ( - 2.0D+00 * log ( r(2*m-1) ) ) &
      * cos ( 2.0D+00 * r8_pi * r(2*m) )

  end if

  return
end
subroutine r8vec_normal_ab ( n, a, b, x )

!*****************************************************************************80
!
!! r8vec_normal_ab() returns a scaled pseudonormal R8VEC.
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
!  Parameters:
!
!    Input, integer N, the number of values desired.
!
!    Input, real ( kind = rk ) A, B, the mean and standard deviation.
!
!    Output, real ( kind = rk ) X(N), a sample of the standard normal PDF.
!
!  Local:
!
!    real ( kind = rk ) R(N+1), is used to store some uniform
!    random values.  Its dimension is N+1, but really it is only needed
!    to be the smallest even number greater than or equal to N.
!
!    integer X_LO_INDEX, X_HI_INDEX, records the range
!    of entries of X that we need to compute. 
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a
  real ( kind = rk ) b
  integer m
  real ( kind = rk ) r(n+1)
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) x(n)
  integer x_hi_index
  integer x_lo_index
!
!  Record the range of X we need to fill in.
!
  x_lo_index = 1
  x_hi_index = n
!
!  If we need just one new value, do that here to avoid null arrays.
!
  if ( x_hi_index - x_lo_index + 1 == 1 ) then

    call random_number ( harvest = r(1:2) )

    x(x_hi_index) = &
      sqrt ( - 2.0D+00 * log ( r(1) ) ) * cos ( 2.0D+00 * r8_pi * r(2) )
!
!  If we require an even number of values, that's easy.
!
  else if ( mod ( x_hi_index - x_lo_index, 2 ) == 1 ) then

    m = ( x_hi_index - x_lo_index + 1 ) / 2

    call random_number ( harvest = r(1:2*m) )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( - 2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * cos ( 2.0D+00 * r8_pi * r(2:2*m:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( - 2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * sin ( 2.0D+00 * r8_pi * r(2:2*m:2) )
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
      sqrt ( - 2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * cos ( 2.0D+00 * r8_pi * r(2:2*m-2:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( - 2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * sin ( 2.0D+00 * r8_pi * r(2:2*m-2:2) )

    x(n) = sqrt ( - 2.0D+00 * log ( r(2*m-1) ) ) &
      * cos ( 2.0D+00 * r8_pi * r(2*m) )

  end if

  x(1:n) = a + b * x(1:n)

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! r8vec_print() prints an R8VEC.
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
subroutine r8vec_uniform_01 ( n, r )

!*****************************************************************************80
!
!! r8vec_uniform_01() returns a unit pseudorandom R8VEC.
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
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Output, real ( kind = rk ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) r(n)

  call random_number ( harvest = r(1:n) )

  return
end
subroutine random_seed_initialize ( key )

!*****************************************************************************80
!
!! random_seed_initialize() initializes the FORTRAN90 random number generator.
!
!  Discussion:
!
!    This is the stupidest, most awkward procedure I have seen!
!
!  Modified:
!
!    27 October 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer KEY: an initial seed for the random number generator.
!
  implicit none

  integer key
  integer, allocatable :: seed(:)
  integer seed_size

  call random_seed ( size = seed_size )
  allocate ( seed(seed_size) )
  seed(1:seed_size) = key
  call random_seed ( put = seed )
  deallocate ( seed )

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
