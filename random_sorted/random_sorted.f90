subroutine normal_01_cdf_inv ( p, x )

!*****************************************************************************80
!
!! normal_01_cdf_inv() inverts the standard normal CDF.
!
!  Discussion:
!
!    The result is accurate to about 1 part in 10^16.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 June 2007
!
!  Author:
!
!    Original FORTRAN77 version by Michael Wichura.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Michael Wichura,
!    Algorithm AS241:
!    The Percentage Points of the Normal Distribution,
!    Applied Statistics,
!    Volume 37, Number 3, pages 477-484, 1988.
!
!  Parameters:
!
!    Input, real ( kind = rk ) P, the value of the cumulative probability
!    densitity function.  0 < P < 1.  If P is outside this range, an
!    "infinite" value will be returned.
!
!    Output, real ( kind = rk ) X, the normal deviate value
!    with the property that the probability of a standard normal deviate being
!    less than or equal to the value is P.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), parameter, dimension ( 8 ) :: a = (/ &
    3.3871328727963666080D+00, &
    1.3314166789178437745D+02, &
    1.9715909503065514427D+03, &
    1.3731693765509461125D+04, &
    4.5921953931549871457D+04, &
    6.7265770927008700853D+04, &
    3.3430575583588128105D+04, &
    2.5090809287301226727D+03 /)
  real ( kind = rk ), parameter, dimension ( 8 ) :: b = (/ &
    1.0D+00, &
    4.2313330701600911252D+01, &
    6.8718700749205790830D+02, &
    5.3941960214247511077D+03, &
    2.1213794301586595867D+04, &
    3.9307895800092710610D+04, &
    2.8729085735721942674D+04, &
    5.2264952788528545610D+03 /)
  real ( kind = rk ), parameter, dimension ( 8 ) :: c = (/ &
    1.42343711074968357734D+00, &
    4.63033784615654529590D+00, &
    5.76949722146069140550D+00, &
    3.64784832476320460504D+00, &
    1.27045825245236838258D+00, &
    2.41780725177450611770D-01, &
    2.27238449892691845833D-02, &
    7.74545014278341407640D-04 /)
  real ( kind = rk ), parameter :: const1 = 0.180625D+00
  real ( kind = rk ), parameter :: const2 = 1.6D+00
  real ( kind = rk ), parameter, dimension ( 8 ) :: d = (/ &
    1.0D+00, &
    2.05319162663775882187D+00, &
    1.67638483018380384940D+00, &
    6.89767334985100004550D-01, &
    1.48103976427480074590D-01, &
    1.51986665636164571966D-02, &
    5.47593808499534494600D-04, &
    1.05075007164441684324D-09 /)
  real ( kind = rk ), parameter, dimension ( 8 ) :: e = (/ &
    6.65790464350110377720D+00, &
    5.46378491116411436990D+00, &
    1.78482653991729133580D+00, &
    2.96560571828504891230D-01, &
    2.65321895265761230930D-02, &
    1.24266094738807843860D-03, &
    2.71155556874348757815D-05, &
    2.01033439929228813265D-07 /)
  real ( kind = rk ), parameter, dimension ( 8 ) :: f = (/ &
    1.0D+00, &
    5.99832206555887937690D-01, &
    1.36929880922735805310D-01, &
    1.48753612908506148525D-02, &
    7.86869131145613259100D-04, &
    1.84631831751005468180D-05, &
    1.42151175831644588870D-07, &
    2.04426310338993978564D-15 /)
  real ( kind = rk ) p
  real ( kind = rk ) q
  real ( kind = rk ) r
  real ( kind = rk ) r8poly_value_horner
  real ( kind = rk ), parameter :: split1 = 0.425D+00
  real ( kind = rk ), parameter :: split2 = 5.0D+00
  real ( kind = rk ) x

  if ( p <= 0.0D+00 ) then
    x = - huge ( x )
    return
  end if

  if ( 1.0D+00 <= p ) then
    x = huge ( x )
    return
  end if

  q = p - 0.5D+00

  if ( abs ( q ) <= split1 ) then

    r = const1 - q * q
    x = q * r8poly_value_horner ( 7, a, r ) / r8poly_value_horner ( 7, b, r )

  else

    if ( q < 0.0D+00 ) then
      r = p
    else
      r = 1.0D+00 - p
    end if

    if ( r <= 0.0D+00 ) then

      x = huge ( x )

    else

      r = sqrt ( - log ( r ) )

      if ( r <= split2 ) then

        r = r - const2
        x = r8poly_value_horner ( 7, c, r ) / r8poly_value_horner ( 7, d, r )

      else

        r = r - split2
        x = r8poly_value_horner ( 7, e, r ) / r8poly_value_horner ( 7, f, r )

      end if

    end if

    if ( q < 0.0D+00 ) then
      x = -x
    end if

  end if

  return
end
function r8poly_value_horner ( m, c, x )

!*****************************************************************************80
!
!! R8POLY_VALUE_HORNER evaluates a polynomial using Horner's method.
!
!  Discussion:
!
!    The polynomial 
!
!      p(x) = c0 + c1 * x + c2 * x^2 + ... + cm * x^m
!
!    is to be evaluated at the value X.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 January 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the degree.
!
!    Input, real ( kind = rk ) C(0:M), the polynomial coefficients.  
!    C(I) is the coefficient of X^I.
!
!    Input, real ( kind = rk ) X, the evaluation point.
!
!    Output, real ( kind = rk ) R8POLY_VALUE_HORNER, the polynomial value.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m

  real ( kind = rk ) c(0:m)
  integer i
  real ( kind = rk ) r8poly_value_horner
  real ( kind = rk ) value
  real ( kind = rk ) x

  value = c(m)
  do i = m - 1, 0, -1
    value = value * x + c(i)
  end do

  r8poly_value_horner = value

  return
end
subroutine r8vec_normal_01_sorted ( n, r8vec )

!*****************************************************************************80
!
!! r8vec_normal_01_sorted() returns a sorted normal 01 random vector.
!
!  Discussion:
!
!    The Normal 01 distribution has mean 0 and standard deviation 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 March 2016
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jon Bentley, James Saxe,
!    Generating sorted lists of random numbers,
!    ACM Transactions on Mathematical Software,
!    Volume 6, Number 3, September 1980, pages 359-364.
!
!  Parameters:
!
!    Input, integer N, the number of values to generate.
!
!    Output, real R8VEC(N), a real vector of normal 01 random values
!    in ascending order.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) cdfvec(n)
  integer i
  real ( kind = rk ) r8vec(n)

  call r8vec_uniform_01_sorted1 ( n, cdfvec )
  
  do i = 1, n
    call normal_01_cdf_inv ( cdfvec(i), r8vec(i) )
  end do

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
subroutine r8vec_transpose_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_TRANSPOSE_PRINT prints an R8VEC "transposed".
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Example:
!
!    A = (/ 1.0, 2.1, 3.2, 4.3, 5.4, 6.5, 7.6, 8.7, 9.8, 10.9, 11.0 /)
!    TITLE = 'My vector:  '
!
!    My vector:  1.0    2.1    3.2    4.3    5.4
!                6.5    7.6    8.7    9.8   10.9
!               11.0
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 May 2014
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
  integer ihi
  integer ilo
  character ( len = * ) title
  integer title_length

  title_length = len_trim ( title )

  do ilo = 1, n, 5
    if ( ilo == 1 ) then
      write ( *, '(a)', advance = 'NO' ) trim ( title )
    else
      write ( *, '(a)', advance = 'NO' ) ( ' ', i = 1, title_length )
    end if
    write ( *, '(2x)', advance = 'NO' )
    ihi = min ( ilo + 5 - 1, n )
    write ( *, '(5g14.6)' ) a(ilo:ihi)
  end do

  return
end
subroutine r8vec_uniform_01_sorted1 ( n, r8vec )

!*****************************************************************************80
!
!! r8vec_uniform_01_sorted1() returns a sorted real random vector in [0,1].
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 March 2016
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jon Bentley, James Saxe,
!    Generating sorted lists of random numbers,
!    ACM Transactions on Mathematical Software,
!    Volume 6, Number 3, September 1980, pages 359-364.
!
!  Parameters:
!
!    Input, integer N, the number of values to generate.
!
!    Output, real R8VEC(N), a real vector of random values
!    in ascending order.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  real ( kind = rk ) r
  real ( kind = rk ) r8vec(n)
  real ( kind = rk ) s

  s = 0.0D+00

  do i = 1, n + 1
    call random_number ( harvest = r )
    s = s - log ( r )
    if ( i == n + 1 ) then
      exit
    end if
    r8vec(i) = s
  end do
  
  r8vec(1:n) = r8vec(1:n) / s

  return
end
subroutine r8vec_uniform_01_sorted2 ( n, r8vec )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01_SORTED2() returns a sorted real random vector in [0,1].
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 March 2016
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jon Bentley, James Saxe,
!    Generating sorted lists of random numbers,
!    ACM Transactions on Mathematical Software,
!    Volume 6, Number 3, September 1980, pages 359-364.
!
!  Parameters:
!
!    Input, integer N, the number of values to generate.
!
!    Input, integer SEED, the integer "seed" used to 
!    generate the output random number.  SEED should not be 0.
!    On output, the updated seed.
!
!    Output, real R8VEC(N), a real vector of random values
!    in ascending order.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) curmax
  integer i
  real ( kind = rk ) r(n)
  real ( kind = rk ) r8vec(n)

  call random_number ( harvest = r(1:n) )

  curmax = 1.0D+00

  do i = n, 1, -1
    curmax = curmax * exp ( log ( r(i) ) / real ( i, kind = rk ) )
    r8vec(i) = curmax
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
