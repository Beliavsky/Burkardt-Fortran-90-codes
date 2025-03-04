subroutine cce_order ( l, n )

!*****************************************************************************80
!
!! cce_order(): order of a (Clenshaw Curtis Exponential (CCE) rule from the level.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 June 2012
!
!  Author:
!
!    John Burkardt.
!
!  Parameters:
!
!    Input, integer L, the level of the rule.  
!    1 <= L.
!
!    Output, integer N, the order of the rule.
!
  implicit none

  integer l
  integer n

  if ( l < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CCE_ORDER - Fatal error!'
    write ( *, '(a)' ) '  1 <= L required.'
    write ( *, '(a,i4)' ) '  Input L = ', l
    stop 1
  else if ( l == 1 ) then
    n = 1
  else
    n = 2 ** ( l - 1 ) + 1
  end if

  return
end
subroutine ccl_order ( l, n )

!*****************************************************************************80
!
!! CCL_ORDER: order of a Clenshaw Curtis Linear (CCL) quadrature rule.
!
!  Discussion:
!
!    Our convention is that the abscissas are numbered from left to right.
!
!    The rule is defined on [0,1].
!
!    The integral to approximate:
!
!      Integral ( 0 <= X <= 1 ) F(X) dX
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
!
!     L  2*L-1   N
!    --  -----  --
!     1      1   1
!     2      3   3
!     3      5   5
!     4      7   7
!     5      9   9
!     6     11  11
!     7     13  13
!     8     15  15
!     9     17  17
!    10     19  19
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer L, the level of the rule.
!    1 <= L.
!
!    Output, integer N, the appropriate order.
!
  implicit none

  integer l
  integer n

  if ( l < 1 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'CCL_ORDER - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal value of L = ', l
    stop 1
  end if

  n = 2 * l - 1

  return
end
subroutine ccs_order ( l, n )

!*****************************************************************************80
!
!! CCS_ORDER: order of a Clenshaw Curtis Slow (CCS) quadrature rule.
!
!  Discussion:
!
!    Our convention is that the abscissas are numbered from left to right.
!
!    The rule is defined on [0,1].
!
!    The integral to approximate:
!
!      Integral ( 0 <= X <= 1 ) F(X) dX
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
!
!    The input value L requests a rule of precision at least 2*L-1.
!
!    In order to preserve nestedness, this function returns the order
!    of a rule which is the smallest value of the form 1+2^E which
!    is greater than or equal to 2*L-1.
!
!     L  2*L-1   N
!    --  -----  --
!     1      1   1
!     2      3   3
!     3      5   5
!     4      7   9
!     5      9   9
!     6     11  17
!     7     13  17
!     8     15  17
!     9     17  17
!    10     19  33
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer L, the level of the rule.
!    1 <= L.
!
!    Output, integer N, the appropriate order.
!
  implicit none

  integer l
  integer n

  if ( l < 1 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'CCS_ORDER - Fatal error!'
    write ( *, '(a,i4)' ) '  Illegal value of L = ', l
    stop 1
  end if
!
!  Find the order N that satisfies the precision requirement.
!
  if ( l == 1 ) then
    n = 1
  else
    n = 3
    do while ( n < 2 * l - 1 )
      n = 2 * n - 1
    end do
  end if

  return
end
subroutine cc ( n, x, w )

!*****************************************************************************80
!
!! CC computes a Clenshaw Curtis (CC) quadrature rule based on order.
!
!  Discussion:
!
!    Our convention is that the abscissas are numbered from left to right.
!
!    The rule is defined on [0,1].
!
!    The integral to approximate:
!
!      Integral ( 0 <= X <= 1 ) F(X) dX
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the rule.
!    1 <= N.
!
!    Output, real ( kind = rk ) X(N), the abscissas.
!
!    Output, real ( kind = rk ) W(N), the weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) d
  real ( kind = rk ) e
  integer i
  integer j
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) theta
  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n)

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CC - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    stop 1
  end if

  if ( n == 1 ) then

    x(1) = 0.0D+00
    w(1) = 2.0D+00

  else

    do i = 1, n
      x(i) = cos ( real ( n - i, kind = rk ) * pi &
                 / real ( n - 1, kind = rk ) )
    end do

    x(1) = -1.0D+00
    if ( mod ( n, 2 ) == 1 ) then
      x((n+1)/2) = 0.0D+00
    end if
    x(n) = +1.0D+00

    do i = 1, n

      theta = real ( i - 1, kind = rk ) * pi &
            / real ( n - 1, kind = rk )

      w(i) = 1.0D+00

      do j = 1, ( n - 1 ) / 2

        if ( 2 * j == ( n - 1 ) ) then
          e = 1.0D+00
        else
          e = 2.0D+00
        end if

        w(i) = w(i) - e * cos ( 2.0D+00 * real ( j, kind = rk ) * theta ) &
             / real ( 4 * j * j - 1, kind = rk )

      end do

    end do

    w(1)     =           w(1)     / real ( n - 1, kind = rk )
    w(2:n-1) = 2.0D+00 * w(2:n-1) / real ( n - 1, kind = rk )
    w(n)     =           w(n)     / real ( n - 1, kind = rk )

  end if
!
!  Transform from [-1,+1] to [0,1].
!
  a = -1.0D+00
  b = +1.0D+00
  c =  0.0D+00
  d = +1.0D+00
  call rule_adjust ( a, b, c, d, n, x, w )

  return
end
function fn_integral ( d )

!*****************************************************************************80
!
!! fn_integral evaluates the integral of the Hermite test function.
!
!  Discussion:
!
!    We are working in spatial dimension D.
!
!    We have variables X1, X2, ..., XN.
!
!    Our test function is X1^6.
!
!    Our test integral is 
!      (2 pi)^(-D/2) integral ( -oo, +oo )^D X1^6 exp(-(X1^2+X2^2+..+XD^2)/2) dX1 dX2 ... dXD
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 May 2012
!
!  Author:
!
!    John Burkardt.
!
!  Parameters:
!
!    Input, integer D, the spatial dimension.
!
!    Output, real ( kind = rk ) FN_INTEGRAL, the integral value.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer d
  integer, parameter :: exponent = 6
  real ( kind = rk ) fn_integral
  integer i4_factorial2

  fn_integral = real ( i4_factorial2 ( exponent - 1 ), kind = rk )

  return
end
subroutine fn_value ( d, n, x, fx )

!*****************************************************************************80
!
!! fn_value evaluates a Hermite test function.
!
!  Discussion:
!
!    We are working in spatial dimension D.
!
!    We have variables X1, X2, ..., XN.
!
!    Our test function is X1^6.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 May 2012
!
!  Author:
!
!    John Burkardt.
!
!  Parameters:
!
!    Input, integer D, the spatial dimension.
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) X(D,N), the points.
!
!    Output, real ( kind = rk ) FX(N), the function values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer d
  integer n

  integer, parameter :: exponent = 6
  real ( kind = rk ) fx(n)
  real ( kind = rk ) x(d,n)

  fx(1:n) = x(1,1:n)**exponent

  return
end
function fu_integral ( d )

!*****************************************************************************80
!
!! FU_INTEGRAL is the integral of the test function for the [0,1]^D interval.
!
!  Discussion:
!
!    The same function, integrated over [-1,+1]^D, has an integral that is
!    2^D times larger.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 May 2012
!
!  Author:
!
!    John Burkardt.
!
!  Parameters:
!
!    Input, integer D, the spatial dimension.
!
!    Output, real ( kind = rk ) FU_INTEGRAL, the integral value.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer d
  real ( kind = rk ) fu_integral

  fu_integral = ( 0.5D+00 * erf ( 0.5D+00 / sqrt ( 2.0D+00 ) ) ) ** d

  return
end
subroutine fu_value ( d, n, x, fx )

!*****************************************************************************80
!
!! FU_VALUE is a sample function for the [0,1]^D interval.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 May 2012
!
!  Author:
!
!    John Burkardt.
!
!  Parameters:
!
!    Input, integer D, the spatial dimension.
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = rk ) X(D,N), the points.
!
!    Output, real ( kind = rk ) FX(N), the function values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer d
  integer n

  real ( kind = rk ) fx(n)
  integer i
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) x(d,n)

  fx(1:n) = 1.0D+00

  do i = 1, d
    fx(1:n) = fx(1:n) * exp ( - ( x(i,1:n) / 2.0D+00 )**2 / 2.0D+00 ) &
      / 2.0D+00 / sqrt ( 2.0D+00 * pi )
  end do

  return
end
subroutine get_seq ( d, norm, seq_num, fs )

!*****************************************************************************80
!
!! GET_SEQ generates all positive integer D-vectors that sum to NORM.
!
!  Discussion:
!
!    This function computes a list, in reverse dictionary order, of
!    all D-vectors of positive values that sum to NORM.
!
!    For example, call get_seq ( 3, 5, 6, fs ) returns
!
!      3  1  1
!      2  2  1
!      2  1  2
!      1  3  1
!      1  2  2
!      1  1  3
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 May 2010
!
!  Author:
!
!    Original MATLAB version by Florian Heiss, Viktor Winschel.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Florian Heiss, Viktor Winschel,
!    Likelihood approximation by numerical integration on sparse grids,
!    Journal of Econometrics,
!    Volume 144, 2008, pages 62-80.
!
!  Parameters:
!
!    Input, integer D, the dimension.
!    1 <= D.
!
!    Input, integer NORM, the value that each row must sum to.
!    D <= NORM.
!
!    Input, integer SEQ_NUM, the number of rows of FS.
!
!    Output, integer FS(SEQ_NUM,D).  Each row of FS represents 
!    one vector with all elements positive and summing to NORM.
!
  implicit none

  integer d
  integer seq_num

  integer a
  integer c
  integer fs(seq_num,d)
  integer i
  integer norm
  integer row
  integer seq(d)

  if ( norm < d ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GET_SEQ - Fatal error!'
    write ( *, '(a,i4,a,i4)' ) '  NORM = ', norm, ' < D = ', d
    stop 1
  end if

  seq(1:d) = 0
!
!  The algorithm is written to work with vectors whose minimum value is
!  allowed to be zero.  So we subtract D from NORM at the beginning and
!  then increment the result vectors by 1 at the end!
!
  a = norm - d
  seq(1) = a

  row = 1
  fs(row,1:d) = seq(1:d) + 1
  c = 1

  do while ( seq ( d ) < a )

    if ( c == d ) then
      do i = c - 1, 1, -1
        c = i
        if ( seq(i) /= 0 ) then
          exit
        end if
      end do
    end if

    seq(c) = seq(c) - 1
    c = c + 1
    seq(c) = a - sum ( seq(1:(c-1)) )

    if ( c < d ) then
      seq((c+1):d) = 0
    end if

    row = row + 1
    fs(row,1:d) = seq(1:d) + 1

  end do

  return
end
subroutine gqn ( n, x, w )

!*****************************************************************************80
!
!! GQN provides data for Gauss quadrature with a normal weight.
!
!  Discussion:
!
!    This data assumes integration over the interval (-oo,+oo) with 
!    weight function w(x) = exp(-x*x/2)/sqrt(2*pi).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 December 2012
!
!  Author:
!
!    Original MATLAB version by Florian Heiss, Viktor Winschel.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Florian Heiss, Viktor Winschel,
!    Likelihood approximation by numerical integration on sparse grids,
!    Journal of Econometrics,
!    Volume 144, 2008, pages 62-80.
!
!  Parameters:
!
!    Input, integer N, the number of points and weights.
!    1 <= N <= 25.
!
!    Output, real ( kind = rk ) X(N), the nodes.
!
!    Output, real ( kind = rk ) W(N), the weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n)
  real ( kind = rk ), dimension ( 1 ) :: x01 = (/ &
    0.00000000000000000D+00 /)
  real ( kind = rk ), dimension ( 1 ) :: w01 = (/ &
    1.00000000000000000D+00 /)
  real ( kind = rk ), dimension ( 2 ) :: x02 = (/ &
   -1.00000000000000000D+00,    &
    1.00000000000000000D+00 /)
  real ( kind = rk ), dimension ( 2 ) :: w02 = (/ &
    0.50000000000000000D+00,    &
    0.50000000000000000D+00 /)
  real ( kind = rk ), dimension ( 3 ) :: x03 = (/ &
   -1.73205080756887719D+00,    &
    0.00000000000000000D+00,    &
    1.73205080756887719D+00 /)
  real ( kind = rk ), dimension ( 3 ) :: w03 = (/ &
    0.166666666666666741D+00,  &
    0.66666666666666663D+00,   &
    0.166666666666666741D+00 /)
  real ( kind = rk ), dimension ( 4 ) :: x04 = (/ &
   -2.33441421833897733D+00,   &
   -0.741963784302725915D+00,  &
    0.741963784302725915D+00,  &
    2.33441421833897733D+00 /)
  real ( kind = rk ), dimension ( 4 ) :: w04 = (/ &
    0.0458758547680684983D+00, &
    0.454124145231931453D+00,  &
    0.454124145231931453D+00,  &
    0.0458758547680684983D+00 /)
  real ( kind = rk ), dimension ( 5 ) :: x05 = (/ &
   -2.85697001387280558D+00,   &
   -1.35562617997426593D+00,   &
    0.00000000000000000D+00,   &
    1.35562617997426593D+00,   &
    2.85697001387280558D+00 /)
  real ( kind = rk ), dimension ( 5 ) :: w05 = (/ &
    0.011257411327720691D+00,  &
    0.22207592200561263D+00,   &
    0.533333333333333437D+00,  &
    0.22207592200561263D+00,   &
    0.011257411327720691D+00 /)
  real ( kind = rk ), dimension ( 6 ) :: x06 = (/ &
   -3.32425743355211933D+00,   &
   -1.88917587775371087D+00,   &
   -0.616706590192594217D+00,  &
    0.616706590192594217D+00,  &
    1.88917587775371087D+00,   &
    3.32425743355211933D+00 /)
  real ( kind = rk ), dimension ( 6 ) :: w06 = (/ &
    0.00255578440205624308D+00, &
    0.0886157460419145226D+00,  &
    0.408828469556029195D+00,   &
    0.408828469556029195D+00,   &
    0.0886157460419145226D+00,  &
    0.00255578440205624308D+00 /)
  real ( kind = rk ), dimension ( 7 ) :: x07 = (/ &
   -3.75043971772574247D+00,   &
   -2.36675941073454155D+00,   &
   -1.15440539473996817D+00,   &
    0.00000000000000000D+00,   &
    1.15440539473996817D+00,   &
    2.36675941073454155D+00,   &
    3.75043971772574247D+00 /)
  real ( kind = rk ), dimension ( 7 ) :: w07 = (/ &
    0.000548268855972218754D+00, &
    0.0307571239675864909D+00,   &
    0.240123178605012505D+00,    &
    0.457142857142857573D+00,    &
    0.240123178605012505D+00,    &
    0.0307571239675864909D+00,   &
    0.000548268855972218754D+00 /)
  real ( kind = rk ), dimension ( 8 ) :: x08 = (/ &
   -4.14454718612589446D+00,   &
   -2.80248586128754162D+00,   &
   -1.63651904243510815D+00,   &
   -0.539079811351375171D+00,  &
    0.539079811351375171D+00,  &
    1.63651904243510815D+00,   &
    2.80248586128754162D+00,   &
    4.14454718612589446D+00 /)
  real ( kind = rk ), dimension ( 8 ) :: w08 = (/ &
    0.000112614538375367836D+00, &
    0.00963522012078826297D+00,  &
    0.117239907661758971D+00,    &
    0.373012257679077364D+00,    &
    0.373012257679077364D+00,    &
    0.117239907661758971D+00,    &
    0.00963522012078826297D+00,  &
    0.000112614538375367836D+00 /)
  real ( kind = rk ), dimension ( 9 ) :: x09 = (/ &
   -4.51274586339978256D+00,   &
   -3.20542900285647026D+00,   &
   -2.07684797867783022D+00,   &
   -1.02325566378913257D+00,   &
    0.00000000000000000D+00,   &
    1.02325566378913257D+00,   &
    2.07684797867783022D+00,   &
    3.20542900285647026D+00,   &
    4.51274586339978256D+00 /)
  real ( kind = rk ), dimension ( 9 ) :: w09 = (/ &
    2.23458440077465626D-05,   &
    0.0027891413212317675D+00, &
    0.0499164067652179688D+00, &
    0.244097502894939089D+00,  &
    0.406349206349206848D+00,  &
    0.244097502894939089D+00,  &
    0.0499164067652179688D+00, &
    0.0027891413212317675D+00, &
    2.23458440077465626D-05 /)
  real ( kind = rk ), dimension ( 10 ) :: x10 = (/ &
   -4.85946282833231269D+00,  &
   -3.58182348355192692D+00,  &
   -2.48432584163895465D+00,  &
   -1.46598909439115821D+00,  &
   -0.484935707515497638D+00, &
    0.484935707515497638D+00, &
    1.46598909439115821D+00,  &
    2.48432584163895465D+00,  &
    3.58182348355192692D+00,  &
    4.85946282833231269D+00 /)
  real ( kind = rk ), dimension ( 10 ) :: w10 = (/ &
    4.3106526307183106D-06,      &
    0.000758070934312219725D+00, &
    0.0191115805007703171D+00,   &
    0.135483702980267295D+00,    &
    0.344642334932019401D+00,    &
    0.344642334932019401D+00,    &
    0.135483702980267295D+00,    &
    0.0191115805007703171D+00,   &
    0.000758070934312219725D+00, &
    4.3106526307183106D-06 /)
  real ( kind = rk ), dimension ( 11 ) :: x11 = (/ &
   -5.18800122437487143D+00,  &
   -3.93616660712997746D+00,  &
   -2.8651231606436447D+00,   &
   -1.87603502015484591D+00,  &
   -0.928868997381063877D+00, &
    0.00000000000000000D+00,  &
    0.928868997381063877D+00, &
    1.87603502015484591D+00,  &
    2.8651231606436447D+00,   &
    3.93616660712997746D+00,  &
    5.18800122437487143D+00 /)
  real ( kind = rk ), dimension ( 11 ) :: w11 = (/ &
    8.12184979021490357D-07,     &
    0.000195671930271223244D+00, &
    0.0067202852355372697D+00,   &
    0.0661387460710576441D+00,   &
    0.242240299873970027D+00,    &
    0.369408369408369575D+00,    &
    0.242240299873970027D+00,    &
    0.0661387460710576441D+00,   &
    0.0067202852355372697D+00,   &
    0.000195671930271223244D+00, &
    8.12184979021490357D-07 /)
  real ( kind = rk ), dimension ( 12 ) :: x12 = (/ &
   -5.50090170446774795D+00,   &
   -4.27182584793228148D+00,   &
   -3.22370982877009737D+00,   &
   -2.25946445100079929D+00,   &
   -1.34037519715161668D+00,   &
   -0.444403001944139009D+00,  &
    0.444403001944139009D+00,  &
    1.34037519715161668D+00,   &
    2.25946445100079929D+00,   &
    3.22370982877009737D+00,   &
    4.27182584793228148D+00,   &
    5.50090170446774795D+00 /)
  real ( kind = rk ), dimension ( 12 ) :: w12 = (/ &
    1.49992716763715968D-07,    &
    4.83718492259060763D-05,    &
    0.00220338068753318491D+00, &
    0.0291166879123641378D+00,  &
    0.146967048045329951D+00,   &
    0.321664361512830066D+00,   &
    0.321664361512830066D+00,   &
    0.146967048045329951D+00,   &
    0.0291166879123641378D+00,  &
    0.00220338068753318491D+00, &
    4.83718492259060763D-05,    &
    1.49992716763715968D-07 /)
  real ( kind = rk ), dimension ( 13 ) :: x13 = (/ &
   -5.8001672523865011D+00,     &
   -4.59139844893652072D+00,    &
   -3.56344438028163468D+00,    &
   -2.62068997343221488D+00,    &
   -1.7254183795882394D+00,     &
   -0.856679493519450053D+00,   &
    0.00000000000000000D+00,    &
    0.856679493519450053D+00,   &
    1.7254183795882394D+00,     &
    2.62068997343221488D+00,    &
    3.56344438028163468D+00,    &
    4.59139844893652072D+00,    &
    5.8001672523865011D+00 /)
  real ( kind = rk ), dimension ( 13 ) :: w13 = (/ &
    2.72262764280590389D-08,     &
    1.15265965273338848D-05,     &
    0.000681236350442926191D+00, &
    0.0117705605059965426D+00,   &
    0.0791689558604501409D+00,   &
    0.237871522964135884D+00,    &
    0.340992340992341492D+00,    &
    0.237871522964135884D+00,    &
    0.0791689558604501409D+00,   &
    0.0117705605059965426D+00,   &
    0.000681236350442926191D+00, &
    1.15265965273338848D-05,     &
    2.72262764280590389D-08 /)
  real ( kind = rk ), dimension ( 14 ) :: x14 = (/ &
   -6.08740954690129144D+00,   &
   -4.89693639734556463D+00,   &
   -3.88692457505976963D+00,   &
   -2.96303657983866753D+00,   &
   -2.08834474570194439D+00,   &
   -1.24268895548546432D+00,   &
   -0.412590457954601808D+00,  &
    0.412590457954601808D+00,  &
    1.24268895548546432D+00,   &
    2.08834474570194439D+00,   &
    2.96303657983866753D+00,   &
    3.88692457505976963D+00,   &
    4.89693639734556463D+00,   &
    6.08740954690129144D+00 /)
  real ( kind = rk ), dimension ( 14 ) :: w14 = (/ &
    4.86816125774838718D-09,    &
    2.66099134406763342D-06,    &
    0.00020033955376074381D+00, &
    0.00442891910694740657D+00, &
    0.0386501088242534319D+00,  &
    0.154083339842513656D+00,   &
    0.302634626813019447D+00,   &
    0.302634626813019447D+00,   &
    0.154083339842513656D+00,   &
    0.0386501088242534319D+00,  &
    0.00442891910694740657D+00, &
    0.00020033955376074381D+00, &
    2.66099134406763342D-06,    &
    4.86816125774838718D-09 /)
  real ( kind = rk ), dimension ( 15 ) :: x15 = (/ &
   -6.36394788882983775D+00,   &
   -5.19009359130478209D+00,   &
   -4.19620771126901548D+00,   &
   -3.28908242439876641D+00,   &
   -2.43243682700975805D+00,   &
   -1.60671006902873015D+00,   &
   -0.799129068324548109D+00,  &
    0.00000000000000000D+00,   &
    0.799129068324548109D+00,  &
    1.60671006902873015D+00,   &
    2.43243682700975805D+00,   &
    3.28908242439876641D+00,   &
    4.19620771126901548D+00,   &
    5.19009359130478209D+00,   &
    6.36394788882983775D+00 /)
  real ( kind = rk ), dimension ( 15 ) :: w15 = (/ &
    8.58964989963318053D-10,    &
    5.97541959792059611D-07,    &
    5.64214640518901565D-05,    &
    0.00156735750354995707D+00, &
    0.0173657744921376159D+00,  &
    0.0894177953998444436D+00,  &
    0.232462293609732223D+00,   &
    0.318259518259518204D+00,   &
    0.232462293609732223D+00,   &
    0.0894177953998444436D+00,  &
    0.0173657744921376159D+00,  &
    0.00156735750354995707D+00, &
    5.64214640518901565D-05,    &
    5.97541959792059611D-07,    &
    8.58964989963318053D-10 /)
  real ( kind = rk ), dimension ( 16 ) :: x16 = (/ &
   -6.63087819839312953D+00,    &
   -5.47222570594934332D+00,    &
   -4.49295530252001196D+00,    &
   -3.60087362417154866D+00,    &
   -2.76024504763070189D+00,    &
   -1.95198034571633361D+00,    &
   -1.1638291005549648D+00,     &
   -0.386760604500557381D+00,   &
    0.386760604500557381D+00,   &
    1.1638291005549648D+00,     &
    1.95198034571633361D+00,    &
    2.76024504763070189D+00,    &
    3.60087362417154866D+00,    &
    4.49295530252001196D+00,    &
    5.47222570594934332D+00,    &
    6.63087819839312953D+00 /)
  real ( kind = rk ), dimension ( 16 ) :: w16 = (/ &
    1.49781472316183141D-10,     &
    1.30947321628682029D-07,     &
    1.53000321624872858D-05,     &
    0.000525984926573909786D+00, &
    0.0072669376011847411D+00,   &
    0.0472847523540140674D+00,   &
    0.158338372750949252D+00,    &
    0.286568521238012408D+00,    &
    0.286568521238012408D+00,    &
    0.158338372750949252D+00,    &
    0.0472847523540140674D+00,   &
    0.0072669376011847411D+00,   &
    0.000525984926573909786D+00, &
    1.53000321624872858D-05,     &
    1.30947321628682029D-07,     &
    1.49781472316183141D-10 /)
  real ( kind = rk ), dimension ( 17 ) :: x17 = (/ &
   -6.88912243989533302D+00,  &
   -5.74446007865940711D+00,  &
   -4.77853158962998403D+00,  &
   -3.90006571719801043D+00,  &
   -3.07379717532819408D+00,  &
   -2.28101944025298886D+00,  &
   -1.50988330779674085D+00,  &
   -0.751842600703896302D+00, &
    0.00000000000000000D+00,  &
    0.751842600703896302D+00, &
    1.50988330779674085D+00,  &
    2.28101944025298886D+00,  &
    3.07379717532819408D+00,  &
    3.90006571719801043D+00,  &
    4.77853158962998403D+00,  &
    5.74446007865940711D+00,  &
    6.88912243989533302D+00 /)
  real ( kind = rk ), dimension ( 17 ) :: w17 = (/ &
    2.58431491937491514D-11,     &
    2.80801611793057831D-08,     &
    4.0126794479798725D-06,      &
    0.000168491431551339447D+00, &
    0.00285894606228464989D+00,  &
    0.023086657025711152D+00,    &
    0.0974063711627180806D+00,   &
    0.226706308468978768D+00,    &
    0.299538370126607556D+00,    &
    0.226706308468978768D+00,    &
    0.0974063711627180806D+00,   &
    0.023086657025711152D+00,    &
    0.00285894606228464989D+00,  &
    0.000168491431551339447D+00, &
    4.0126794479798725D-06,      &
    2.80801611793057831D-08,     &
    2.58431491937491514D-11 /)
  real ( kind = rk ), dimension ( 18 ) :: x18 = (/ &
   -7.13946484914647961D+00,   &
   -6.00774591135959746D+00,   &
   -5.05407268544274046D+00,   &
   -4.1880202316294044D+00,    &
   -3.37473653577809074D+00,   &
   -2.59583368891124033D+00,   &
   -1.83977992150864567D+00,   &
   -1.09839551809150127D+00,   &
   -0.365245755507697667D+00,  &
    0.365245755507697667D+00,  &
    1.09839551809150127D+00,   &
    1.83977992150864567D+00,   &
    2.59583368891124033D+00,   &
    3.37473653577809074D+00,   &
    4.1880202316294044D+00,    &
    5.05407268544274046D+00,   &
    6.00774591135959746D+00,   &
    7.13946484914647961D+00 /)
  real ( kind = rk ), dimension ( 18 ) :: w18 = (/ &
    4.41658876935870775D-12,    &
    5.90548847883654844D-09,    &
    1.02155239763698159D-06,    &
    5.17989614411619621D-05,    &
    0.00106548479629164959D+00, &
    0.0105165177519413525D+00,  &
    0.0548966324802226541D+00,  &
    0.160685303893512627D+00,   &
    0.272783234654287887D+00,   &
    0.272783234654287887D+00,   &
    0.160685303893512627D+00,   &
    0.0548966324802226541D+00,  &
    0.0105165177519413525D+00,  &
    0.00106548479629164959D+00, &
    5.17989614411619621D-05,    &
    1.02155239763698159D-06,    &
    5.90548847883654844D-09,    &
    4.41658876935870775D-12 /)
  real ( kind = rk ), dimension ( 19 ) :: x19 = (/ &
   -7.38257902403043165D+00,  &
   -6.2628911565132519D+00,   &
   -5.32053637733603857D+00,  &
   -4.46587262683103159D+00,  &
   -3.66441654745063827D+00,  &
   -2.89805127651575356D+00,  &
   -2.15550276131693508D+00,  &
   -1.4288766760783731D+00,   &
   -0.712085044042379933D+00, &
    0.00000000000000000D+00,  &
    0.712085044042379933D+00, &
    1.4288766760783731D+00,   &
    2.15550276131693508D+00,  &
    2.89805127651575356D+00,  &
    3.66441654745063827D+00,  &
    4.46587262683103159D+00,  &
    5.32053637733603857D+00,  &
    6.2628911565132519D+00,   &
    7.38257902403043165D+00 /)
  real ( kind = rk ), dimension ( 19 ) :: w19 = (/ &
    7.4828300540572308D-13,      &
    1.22037084844747862D-09,     &
    2.53222003209286807D-07,     &
    1.53511459546667444D-05,     &
    0.000378502109414267593D+00, &
    0.00450723542034203555D+00,  &
    0.0286666910301184956D+00,   &
    0.103603657276143998D+00,    &
    0.220941712199143658D+00,    &
    0.283773192751521075D+00,    &
    0.220941712199143658D+00,    &
    0.103603657276143998D+00,    &
    0.0286666910301184956D+00,   &
    0.00450723542034203555D+00,  &
    0.000378502109414267593D+00, &
    1.53511459546667444D-05,     &
    2.53222003209286807D-07,     &
    1.22037084844747862D-09,     &
    7.4828300540572308D-13 /)
  real ( kind = rk ), dimension ( 20 ) :: x20 = (/ &
   -7.61904854167975909D+00,  &
   -6.51059015701365507D+00,  &
   -5.57873880589320148D+00,  &
   -4.73458133404605519D+00,  &
   -3.9439673506573163D+00,   &
   -3.18901481655339003D+00,  &
   -2.45866361117236787D+00,  &
   -1.74524732081412703D+00,  &
   -1.04294534880275092D+00,  &
   -0.346964157081355917D+00, &
    0.346964157081355917D+00, &
    1.04294534880275092D+00,  &
    1.74524732081412703D+00,  &
    2.45866361117236787D+00,  &
    3.18901481655339003D+00,  &
    3.9439673506573163D+00,   &
    4.73458133404605519D+00,  &
    5.57873880589320148D+00,  &
    6.51059015701365507D+00,  &
    7.61904854167975909D+00 /)
  real ( kind = rk ), dimension ( 20 ) :: w20 = (/ &
    1.25780067243793047D-13,     &
    2.4820623623151838D-10,      &
    6.12749025998295973D-08,     &
    4.40212109023086457D-06,     &
    0.000128826279961928981D+00, &
    0.00183010313108049175D+00,  &
    0.0139978374471010428D+00,   &
    0.0615063720639760295D+00,   &
    0.161739333984000255D+00,    &
    0.26079306344955544D+00,     &
    0.26079306344955544D+00,     &
    0.161739333984000255D+00,    &
    0.0615063720639760295D+00,   &
    0.0139978374471010428D+00,   &
    0.00183010313108049175D+00,  &
    0.000128826279961928981D+00, &
    4.40212109023086457D-06,     &
    6.12749025998295973D-08,     &
    2.4820623623151838D-10,      &
    1.25780067243793047D-13 /)
  real ( kind = rk ), dimension ( 21 ) :: x21 = (/ &
   -7.84938289511382248D+00,   &
   -6.75144471871746088D+00,   &
   -5.82938200730447065D+00,   &
   -4.99496394478202532D+00,   &
   -4.21434398168842161D+00,   &
   -3.46984669047537642D+00,   &
   -2.75059298105237326D+00,   &
   -2.0491024682571628D+00,    &
   -1.35976582321123041D+00,   &
   -0.678045692440644054D+00,  &
    0.00000000000000000D+00,   &
    0.678045692440644054D+00,  &
    1.35976582321123041D+00,   &
    2.0491024682571628D+00,    &
    2.75059298105237326D+00,   &
    3.46984669047537642D+00,   &
    4.21434398168842161D+00,   &
    4.99496394478202532D+00,   &
    5.82938200730447065D+00,   &
    6.75144471871746088D+00,   &
    7.84938289511382248D+00 /)
  real ( kind = rk ), dimension ( 21 ) :: w21 = (/ &
    2.09899121956566525D-14,     &
    4.97536860412174643D-11,     &
    1.45066128449307397D-08,     &
    1.22535483614825217D-06,     &
    4.21923474255158655D-05,     &
    0.000708047795481537364D+00, &
    0.00643969705140877684D+00,  &
    0.0339527297865428387D+00,   &
    0.10839228562641938D+00,     &
    0.215333715695059824D+00,    &
    0.270260183572877066D+00,    &
    0.215333715695059824D+00,    &
    0.10839228562641938D+00,     &
    0.0339527297865428387D+00,   &
    0.00643969705140877684D+00,  &
    0.000708047795481537364D+00, &
    4.21923474255158655D-05,     &
    1.22535483614825217D-06,     &
    1.45066128449307397D-08,     &
    4.97536860412174643D-11,     &
    2.09899121956566525D-14 /)
  real ( kind = rk ), dimension ( 22 ) :: x22 = (/ &
   -8.07402998402171157D+00,   &
   -6.98598042401881525D+00,   &
   -6.07307495112289786D+00,   &
   -5.2477244337144251D+00,    &
   -4.47636197731086849D+00,   &
   -3.74149635026651772D+00,   &
   -3.03240422783167629D+00,   &
   -2.34175999628770803D+00,   &
   -1.6641248391179071D+00,    &
   -0.995162422271215541D+00,  &
   -0.331179315715273814D+00,  &
    0.331179315715273814D+00,  &
    0.995162422271215541D+00,  &
    1.6641248391179071D+00,    &
    2.34175999628770803D+00,   &
    3.03240422783167629D+00,   &
    3.74149635026651772D+00,   &
    4.47636197731086849D+00,   &
    5.2477244337144251D+00,    &
    6.07307495112289786D+00,   &
    6.98598042401881525D+00,   &
    8.07402998402171157D+00 /)
  real ( kind = rk ), dimension ( 22 ) :: w22 = (/ &
    3.47946064787714279D-15,     &
    9.84137898234601051D-12,     &
    3.36651415945821088D-09,     &
    3.31985374981400429D-07,     &
    1.33459771268087124D-05,     &
    0.000262283303255964159D+00, &
    0.00280876104757721073D+00,  &
    0.0175690728808057736D+00,   &
    0.0671963114288898905D+00,   &
    0.161906293413675378D+00,    &
    0.250243596586935013D+00,    &
    0.250243596586935013D+00,    &
    0.161906293413675378D+00,    &
    0.0671963114288898905D+00,   &
    0.0175690728808057736D+00,   &
    0.00280876104757721073D+00,  &
    0.000262283303255964159D+00, &
    1.33459771268087124D-05,     &
    3.31985374981400429D-07,     &
    3.36651415945821088D-09,     &
    9.84137898234601051D-12,     &
    3.47946064787714279D-15 /)
  real ( kind = rk ), dimension ( 23 ) :: x23 = (/ &
    -8.29338602741735365D+00,  &
    -7.21465943505186225D+00,  &
    -6.31034985444839958D+00,  &
    -5.49347398647179475D+00,  &
    -4.73072419745147332D+00,  &
    -4.00477532173330442D+00,  &
    -3.30504002175296518D+00,  &
    -2.62432363405918201D+00,  &
    -1.9573275529334242D+00,   &
    -1.29987646830397896D+00,  &
    -0.648471153534495803D+00, &
     0.00000000000000000D+00,  &
     0.648471153534495803D+00, &
     1.29987646830397896D+00,  &
     1.9573275529334242D+00,   &
     2.62432363405918201D+00,  &
     3.30504002175296518D+00,  &
     4.00477532173330442D+00,  &
     4.73072419745147332D+00,  &
     5.49347398647179475D+00,  &
     6.31034985444839958D+00,  &
     7.21465943505186225D+00,  &
     8.29338602741735365D+00 /)
  real ( kind = rk ), dimension ( 23 ) :: w23 = (/ &
    5.73238316780208728D-16,    &
    1.92293531156779128D-12,    &
    7.67088886239990765D-10,    &
    8.77506248386171607D-08,    &
    4.08997724499215494D-06,    &
    9.34081860903129835D-05,    &
    0.00116762863749786134D+00, &
    0.00857967839146566401D+00, &
    0.0388671837034809467D+00,  &
    0.112073382602620911D+00,   &
    0.209959669577542613D+00,   &
    0.258509740808839039D+00,   &
    0.209959669577542613D+00,   &
    0.112073382602620911D+00,   &
    0.0388671837034809467D+00,  &
    0.00857967839146566401D+00, &
    0.00116762863749786134D+00, &
    9.34081860903129835D-05,    &
    4.08997724499215494D-06,    &
    8.77506248386171607D-08,    &
    7.67088886239990765D-10,    &
    1.92293531156779128D-12,    &
    5.73238316780208728D-16 /)
  real ( kind = rk ), dimension ( 24 ) :: x24 = (/ &
   -8.50780351919525835D+00,   &
   -7.43789066602166304D+00,   &
   -6.54167500509863409D+00,   &
   -5.73274717525120092D+00,   &
   -4.97804137463912078D+00,   &
   -4.26038360501990532D+00,   &
   -3.56930676407356096D+00,   &
   -2.89772864322331403D+00,   &
   -2.24046785169175244D+00,   &
   -1.59348042981642024D+00,   &
   -0.953421922932109256D+00,  &
   -0.317370096629452314D+00,  &
    0.317370096629452314D+00,  &
    0.953421922932109256D+00,  &
    1.59348042981642024D+00,   &
    2.24046785169175244D+00,   &
    2.89772864322331403D+00,   &
    3.56930676407356096D+00,   &
    4.26038360501990532D+00,   &
    4.97804137463912078D+00,   &
    5.73274717525120092D+00,   &
    6.54167500509863409D+00,   &
    7.43789066602166304D+00,   &
    8.50780351919525835D+00 /)
  real ( kind = rk ), dimension ( 24 ) :: w24 = (/ &
    9.39019368904192022D-17,     &
    3.71497415276241595D-13,     &
    1.71866492796486901D-10,     &
    2.26746167348046514D-08,     &
    1.21765974544258296D-06,     &
    3.20950056527459886D-05,     &
    0.000464718718779397633D+00, &
    0.00397660892918131129D+00,  &
    0.0211263444089670287D+00,   &
    0.0720693640171784361D+00,   &
    0.161459512867000249D+00,    &
    0.240870115546640562D+00,    &
    0.240870115546640562D+00,    &
    0.161459512867000249D+00,    &
    0.0720693640171784361D+00,   &
    0.0211263444089670287D+00,   &
    0.00397660892918131129D+00,  &
    0.000464718718779397633D+00, &
    3.20950056527459886D-05,     &
    1.21765974544258296D-06,     &
    2.26746167348046514D-08,     &
    1.71866492796486901D-10,     &
    3.71497415276241595D-13,     &
    9.39019368904192022D-17 /)
  real ( kind = rk ), dimension ( 25 ) :: x25 = (/ &
   -8.71759767839958855D+00,  &
   -7.65603795539307619D+00,  &
   -6.76746496380971685D+00,  &
   -5.96601469060670198D+00,  &
   -5.21884809364427937D+00,  &
   -4.50892992296728501D+00,  &
   -3.82590056997249173D+00,  &
   -3.16277567938819271D+00,  &
   -2.51447330395220581D+00,  &
   -1.8770583699478387D+00,   &
   -1.24731197561678919D+00,  &
   -0.622462279186076106D+00, &
    0.00000000000000000D+00,  &
    0.622462279186076106D+00, &
    1.24731197561678919D+00,  &
    1.8770583699478387D+00,   &
    2.51447330395220581D+00,  &
    3.16277567938819271D+00,  &
    3.82590056997249173D+00,  &
    4.50892992296728501D+00,  &
    5.21884809364427937D+00,  &
    5.96601469060670198D+00,  &
    6.76746496380971685D+00,  &
    7.65603795539307619D+00,  &
    8.71759767839958855D+00 /)
  real ( kind = rk ), dimension ( 25 ) :: w25 = (/ &
    1.53003899799868247D-17,    &
    7.10210303700392527D-14,    &
    3.79115000047718706D-11,    &
    5.7380238688993763D-09,     &
    3.53015256024549785D-07,    &
    1.06721949052025363D-05,    &
    0.0001777669069265266D+00,  &
    0.00175785040526379608D+00, &
    0.0108567559914623159D+00,  &
    0.0433799701676449712D+00,  &
    0.114880924303951637D+00,   &
    0.204851025650340413D+00,   &
    0.248169351176485475D+00,   &
    0.204851025650340413D+00,   &
    0.114880924303951637D+00,   &
    0.0433799701676449712D+00,  &
    0.0108567559914623159D+00,  &
    0.00175785040526379608D+00, &
    0.0001777669069265266D+00,  &
    1.06721949052025363D-05,    &
    3.53015256024549785D-07,    &
    5.7380238688993763D-09,     &
    3.79115000047718706D-11,    &
    7.10210303700392527D-14,    &
    1.53003899799868247D-17 /)

  if ( n == 1 ) then
    call r8vec_copy ( n, x01, x )
    call r8vec_copy ( n, w01, w )
  else if ( n == 2 ) then
    call r8vec_copy ( n, x02, x )
    call r8vec_copy ( n, w02, w )
  else if ( n == 3 ) then
    call r8vec_copy ( n, x03, x )
    call r8vec_copy ( n, w03, w )
  else if ( n == 4 ) then
    call r8vec_copy ( n, x04, x )
    call r8vec_copy ( n, w04, w )
  else if ( n == 5 ) then
    call r8vec_copy ( n, x05, x )
    call r8vec_copy ( n, w05, w )
  else if ( n == 6 ) then
    call r8vec_copy ( n, x06, x )
    call r8vec_copy ( n, w06, w )
  else if ( n == 7 ) then
    call r8vec_copy ( n, x07, x )
    call r8vec_copy ( n, w07, w )
  else if ( n == 8 ) then
    call r8vec_copy ( n, x08, x )
    call r8vec_copy ( n, w08, w )
  else if ( n == 9 ) then
    call r8vec_copy ( n, x09, x )
    call r8vec_copy ( n, w09, w )
  else if ( n == 10 ) then
    call r8vec_copy ( n, x10, x )
    call r8vec_copy ( n, w10, w )
  else if ( n == 11 ) then
    call r8vec_copy ( n, x11, x )
    call r8vec_copy ( n, w11, w )
  else if ( n == 12 ) then
    call r8vec_copy ( n, x12, x )
    call r8vec_copy ( n, w12, w )
  else if ( n == 13 ) then
    call r8vec_copy ( n, x13, x )
    call r8vec_copy ( n, w13, w )
  else if ( n == 14 ) then
    call r8vec_copy ( n, x14, x )
    call r8vec_copy ( n, w14, w )
  else if ( n == 15 ) then
    call r8vec_copy ( n, x15, x )
    call r8vec_copy ( n, w15, w )
  else if ( n == 16 ) then
    call r8vec_copy ( n, x16, x )
    call r8vec_copy ( n, w16, w )
  else if ( n == 17 ) then
    call r8vec_copy ( n, x17, x )
    call r8vec_copy ( n, w17, w )
  else if ( n == 18 ) then
    call r8vec_copy ( n, x18, x )
    call r8vec_copy ( n, w18, w )
  else if ( n == 19 ) then
    call r8vec_copy ( n, x19, x )
    call r8vec_copy ( n, w19, w )
  else if ( n == 20 ) then
    call r8vec_copy ( n, x20, x )
    call r8vec_copy ( n, w20, w )
  else if ( n == 21 ) then
    call r8vec_copy ( n, x21, x )
    call r8vec_copy ( n, w21, w )
  else if ( n == 22 ) then
    call r8vec_copy ( n, x22, x )
    call r8vec_copy ( n, w22, w )
  else if ( n == 23 ) then
    call r8vec_copy ( n, x23, x )
    call r8vec_copy ( n, w23, w )
  else if ( n == 24 ) then
    call r8vec_copy ( n, x24, x )
    call r8vec_copy ( n, w24, w )
  else if ( n == 25 ) then
    call r8vec_copy ( n, x25, x )
    call r8vec_copy ( n, w25, w )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GQN - Fatal error!'
    write ( *, '(a)' ) '  Value of N must be between 1 and 25.'
    stop 1
  end if

  return
end
subroutine gqn_order ( l, n )

!*****************************************************************************80
!
!! GQN_ORDER computes the order of a GQN rule from the level.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 June 2012
!
!  Author:
!
!    John Burkardt.
!
!  Parameters:
!
!    Input, integer L, the level of the rule.  
!    1 <= L.
!
!    Output, integer N, the order of the rule.
!
  implicit none

  integer l
  integer n

  if ( l < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GQN_ORDER - Fatal error!'
    write ( *, '(a)' ) '  1 <= L required.'
    write ( *, '(a,i4)' ) '  Input L = ', l
    stop 1
  else if ( l <= 25 ) then
    n = l
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GQN_ORDER - Fatal error!'
    write ( *, '(a)' ) '  L <= 25 required.'
    write ( *, '(a,i4)' ) '  Input L = ', l
    stop 1
  end if

  return
end
subroutine gqn2_order ( l, n )

!*****************************************************************************80
!
!! GQN2_ORDER computes the order of a GQN rule from the level.
!
!  Discussion:
!
!    For this version of the order routine, we have
!
!      n = 2 * l - 1
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 February 2014
!
!  Author:
!
!    John Burkardt.
!
!  Parameters:
!
!    Input, integer L, the level of the rule.  
!    1 <= L.
!
!    Output, integer N, the order of the rule.
!
  implicit none

  integer l
  integer n

  if ( l < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GQN2_ORDER - Fatal error!'
    write ( *, '(a)' ) '  1 <= L required.'
    write ( *, '(a,i4)' ) '  Input L = ', l
    stop 1
  else if ( l <= 13 ) then
    n = 2 * l - 1
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GQN2_ORDER - Fatal error!'
    write ( *, '(a)' ) '  L <= 13 required.'
    write ( *, '(a,i4)' ) '  Input L = ', l
    stop 1
  end if

  return
end
subroutine gqu ( n, x, w )

!*****************************************************************************80
!
!! GQU provides data for Gauss quadrature with a uniform weight.
!
!  Discussion:
!
!    This data assumes integration over the interval [0,1] with 
!    weight function w(x) = 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 December 2012
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Florian Heiss, Viktor Winschel,
!    Likelihood approximation by numerical integration on sparse grids,
!    Journal of Econometrics,
!    Volume 144, 2008, pages 62-80.
!
!  Parameters:
!
!    Input, integer N, the number of points and weights.
!    1 <= N <= 25.
!
!    Output, real ( kind = rk ) X(N), the nodes.
!
!    Output, real ( kind = rk ) W(N), the weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n)
  real ( kind = rk ), dimension ( 1 ) :: x01 = (/ &
  0.500000000000000000D+00 /)
  real ( kind = rk ), dimension ( 1 ) :: w01 = (/ &
  1.000000000000000000D+00 /)
  real ( kind = rk ), dimension ( 2 ) :: x02 = (/ &
  0.211324865405187134D+00,  &
  0.788675134594812866D+00 /)
  real ( kind = rk ), dimension ( 2 ) :: w02 = (/ &
  0.500000000000000000D+00,  &
  0.500000000000000000D+00 /)
  real ( kind = rk ), dimension ( 3 ) :: x03 = (/ &
  0.112701665379258298D+00,  &
  0.500000000000000000D+00,  &
  0.887298334620741702D+00 /)
  real ( kind = rk ), dimension ( 3 ) :: w03 = (/ &
  0.277777777777777124D+00,  &
  0.444444444444445697D+00,  &
  0.277777777777777124D+00 /)
  real ( kind = rk ), dimension ( 4 ) :: x04 = (/ &
  0.0694318442029737692D+00, &
  0.330009478207571871D+00,  &
  0.669990521792428129D+00,  &
  0.930568155797026231D+00 /)
  real ( kind = rk ), dimension ( 4 ) :: w04 = (/ &
  0.173927422568724843D+00,  &
  0.326072577431275157D+00,  &
  0.326072577431275157D+00,  &
  0.173927422568724843D+00 /)
  real ( kind = rk ), dimension ( 5 ) :: x05 = (/ &
  0.0469100770306680737D+00, &
  0.230765344947158502D+00,  &
  0.500000000000000000D+00,  &
  0.769234655052841498D+00,  &
  0.953089922969331926D+00 /)
  real ( kind = rk ), dimension ( 5 ) :: w05 = (/ &
  0.118463442528091739D+00,  &
  0.239314335249685012D+00,  &
  0.284444444444446554D+00,  &
  0.239314335249685012D+00,  &
  0.118463442528091739D+00 /)
  real ( kind = rk ), dimension ( 6 ) :: x06 = (/ &
  0.033765242898423975D+00,  &
  0.169395306766867648D+00,  &
  0.380690406958401506D+00,  &
  0.619309593041598494D+00,  &
  0.830604693233132352D+00,  &
  0.966234757101576025D+00 /)
  real ( kind = rk ), dimension ( 6 ) :: w06 = (/ &
  0.0856622461895818338D+00, &
  0.180380786524070719D+00,  &
  0.233956967286347461D+00,  &
  0.233956967286347461D+00,  &
  0.180380786524070719D+00,  &
  0.0856622461895818338D+00 /)
  real ( kind = rk ), dimension ( 7 ) :: x07 = (/ &
  0.0254460438286208124D+00,  &
  0.12923440720030277D+00,    &
  0.297077424311301463D+00,   &
  0.500000000000000000D+00,   &
  0.702922575688698537D+00,   &
  0.87076559279969723D+00,    &
  0.974553956171379188D+00 /)
  real ( kind = rk ), dimension ( 7 ) :: w07 = (/ &
  0.0647424830844317012D+00,  &
  0.139852695744639349D+00,   &
  0.190915025252560905D+00,   &
  0.2089795918367362D+00,     &
  0.190915025252560905D+00,   &
  0.139852695744639349D+00,   &
  0.0647424830844317012D+00 /)
  real ( kind = rk ), dimension ( 8 ) :: x08 = (/ &
  0.0198550717512319119D+00,  &
  0.101666761293186525D+00,   &
  0.237233795041835505D+00,   &
  0.408282678752175054D+00,   &
  0.591717321247824946D+00,   &
  0.762766204958164495D+00,   &
  0.898333238706813475D+00,   &
  0.980144928248768088D+00 /)
  real ( kind = rk ), dimension ( 8 ) :: w08 = (/ &
  0.0506142681451851803D+00, &
  0.111190517226687935D+00,  &
  0.156853322938944689D+00,  &
  0.181341891689182133D+00,  &
  0.181341891689182133D+00,  &
  0.156853322938944689D+00,  &
  0.111190517226687935D+00,  &
  0.0506142681451851803D+00 /)
  real ( kind = rk ), dimension ( 9 ) :: x09 = (/ &
  0.0159198802461869571D+00, &
  0.0819844463366821152D+00, &
  0.193314283649704821D+00,  &
  0.337873288298095487D+00,  &
  0.500000000000000000D+00,  &
  0.662126711701904513D+00,  &
  0.806685716350295179D+00,  &
  0.918015553663317885D+00,  &
  0.984080119753813043D+00 /)
  real ( kind = rk ), dimension ( 9 ) :: w09 = (/ &
  0.0406371941807845832D+00, &
  0.0903240803474292531D+00, &
  0.130305348201468441D+00,  &
  0.156173538520002264D+00,  &
  0.165119677500630752D+00,  &
  0.156173538520002264D+00,  &
  0.130305348201468441D+00,  &
  0.0903240803474292531D+00, &
  0.0406371941807845832D+00 /)
  real ( kind = rk ), dimension ( 10 ) :: x10 = (/ &
  0.0130467357414141283D+00, &
  0.0674683166555076763D+00, &
  0.160295215850487782D+00,  &
  0.283302302935376393D+00,  &
  0.42556283050918442D+00,   &
  0.57443716949081558D+00,   &
  0.716697697064623607D+00,  &
  0.839704784149512218D+00,  &
  0.932531683344492324D+00,  &
  0.986953264258585872D+00 /)
  real ( kind = rk ), dimension ( 10 ) :: w10 = (/ &
  0.0333356721543420012D+00, &
  0.0747256745752905988D+00, &
  0.109543181257991576D+00,  &
  0.13463335965499873D+00,   &
  0.147762112357377129D+00,  &
  0.147762112357377129D+00,  &
  0.13463335965499873D+00,   &
  0.109543181257991576D+00,  &
  0.0747256745752905988D+00, &
  0.0333356721543420012D+00 /)
  real ( kind = rk ), dimension ( 11 ) :: x11 = (/ &
  0.010885670926971569D+00,  &
  0.0564687001159522861D+00, &
  0.134923997212975322D+00,  &
  0.240451935396594152D+00,  &
  0.365228422023827548D+00,  &
  0.500000000000000000D+00,  &
  0.634771577976172452D+00,  &
  0.759548064603405848D+00,  &
  0.865076002787024678D+00,  &
  0.943531299884047714D+00,  &
  0.989114329073028431D+00 /)
  real ( kind = rk ), dimension ( 11 ) :: w11 = (/ &
  0.0278342835580849164D+00, &
  0.0627901847324526252D+00, &
  0.0931451054638675197D+00, &
  0.11659688229599563D+00,   &
  0.131402272255123881D+00,  &
  0.136462543388950863D+00,  &
  0.131402272255123881D+00,  &
  0.11659688229599563D+00,   &
  0.0931451054638675197D+00, &
  0.0627901847324526252D+00, &
  0.0278342835580849164D+00 /)
  real ( kind = rk ), dimension ( 12 ) :: x12 = (/ &
  0.00921968287664043373D+00,&
  0.0479413718147625456D+00, &
  0.115048662902847654D+00,  &
  0.206341022856691314D+00,  &
  0.316084250500909936D+00,  &
  0.437383295744265488D+00,  &
  0.562616704255734512D+00,  &
  0.683915749499090064D+00,  &
  0.793658977143308686D+00,  &
  0.884951337097152346D+00,  &
  0.952058628185237454D+00,  &
  0.990780317123359566D+00 /)
  real ( kind = rk ), dimension ( 12 ) :: w12 = (/ &
  0.0235876681932543145D+00, &
  0.0534696629976592758D+00, &
  0.0800391642716734436D+00, &
  0.101583713361533282D+00,  &
  0.116746268269177805D+00,  &
  0.124573522906701886D+00,  &
  0.124573522906701886D+00,  &
  0.116746268269177805D+00,  &
  0.101583713361533282D+00,  &
  0.0800391642716734436D+00, &
  0.0534696629976592758D+00, &
  0.0235876681932543145D+00 /)
  real ( kind = rk ), dimension ( 13 ) :: x13 = (/ &
  0.00790847264070593248D+00, &
  0.0412008003885109275D+00,  &
  0.0992109546333450609D+00,  &
  0.178825330279829942D+00,   &
  0.275753624481776649D+00,   &
  0.384770842022432613D+00,   &
  0.500000000000000000D+00,   &
  0.615229157977567387D+00,   &
  0.724246375518223351D+00,   &
  0.821174669720170058D+00,   &
  0.900789045366654939D+00,   &
  0.958799199611489072D+00,   &
  0.992091527359294068D+00 /)
  real ( kind = rk ), dimension ( 13 ) :: w13 = (/ &
  0.0202420023826562281D+00,  &
  0.0460607499188643785D+00,  &
  0.0694367551098938746D+00,  &
  0.0890729903809732021D+00,  &
  0.103908023768444616D+00,   &
  0.113141590131449032D+00,   &
  0.116275776615437407D+00,   &
  0.113141590131449032D+00,   &
  0.103908023768444616D+00,   &
  0.0890729903809732021D+00,  &
  0.0694367551098938746D+00,  &
  0.0460607499188643785D+00,  &
  0.0202420023826562281D+00 /)
  real ( kind = rk ), dimension ( 14 ) :: x14 = (/ &
  0.00685809565159378742D+00, &
  0.0357825581682131855D+00,  &
  0.0863993424651174902D+00,  &
  0.156353547594157316D+00,   &
  0.24237568182092295D+00,    &
  0.340443815536055183D+00,   &
  0.445972525646328166D+00,   &
  0.554027474353671834D+00,   &
  0.659556184463944817D+00,   &
  0.75762431817907705D+00,    &
  0.843646452405842684D+00,   &
  0.91360065753488251D+00,    &
  0.964217441831786815D+00,   &
  0.993141904348406213D+00 /)
  real ( kind = rk ), dimension ( 14 ) :: w14 = (/ &
  0.0175597301658745736D+00,  &
  0.0400790435798802913D+00,  &
  0.0607592853439517105D+00,  &
  0.0786015835790969952D+00,  &
  0.0927691987389691608D+00,  &
  0.102599231860648107D+00,   &
  0.107631926731579161D+00,   &
  0.107631926731579161D+00,   &
  0.102599231860648107D+00,   &
  0.0927691987389691608D+00,  &
  0.0786015835790969952D+00,  &
  0.0607592853439517105D+00,  &
  0.0400790435798802913D+00,  &
  0.0175597301658745736D+00 /)
  real ( kind = rk ), dimension ( 15 ) :: x15 = (/ &
  0.0060037409897573113D+00, &
  0.0313633037996470243D+00, &
  0.0758967082947863414D+00, &
  0.137791134319914965D+00,  &
  0.214513913695730585D+00,  &
  0.302924326461218252D+00,  &
  0.399402953001282701D+00,  &
  0.500000000000000000D+00,  &
  0.600597046998717299D+00,  &
  0.697075673538781748D+00,  &
  0.785486086304269415D+00,  &
  0.862208865680085035D+00,  &
  0.924103291705213659D+00,  &
  0.968636696200352976D+00,  &
  0.993996259010242689D+00 /)
  real ( kind = rk ), dimension ( 15 ) :: w15 = (/ &
  0.0153766209980574341D+00, &
  0.0351830237440541593D+00, &
  0.0535796102335861571D+00, &
  0.0697853389630773147D+00, &
  0.083134602908497196D+00,  &
  0.0930805000077812861D+00, &
  0.0992157426635560391D+00, &
  0.101289120962780907D+00,  &
  0.0992157426635560391D+00, &
  0.0930805000077812861D+00, &
  0.083134602908497196D+00,  &
  0.0697853389630773147D+00, &
  0.0535796102335861571D+00, &
  0.0351830237440541593D+00, &
  0.0153766209980574341D+00 /)
  real ( kind = rk ), dimension ( 16 ) :: x16 = (/ &
  0.00529953250417503074D+00, &
  0.0277124884633836999D+00,  &
  0.0671843988060840669D+00,  &
  0.122297795822498445D+00,   &
  0.191061877798678115D+00,   &
  0.270991611171386371D+00,   &
  0.359198224610370542D+00,   &
  0.452493745081181231D+00,   &
  0.547506254918818769D+00,   &
  0.640801775389629458D+00,   &
  0.729008388828613629D+00,   &
  0.808938122201321885D+00,   &
  0.877702204177501555D+00,   &
  0.932815601193915933D+00,   &
  0.9722875115366163D+00,     &
  0.994700467495824969D+00 /)
  real ( kind = rk ), dimension ( 16 ) :: w16 = (/ &
  0.0135762297058759553D+00, &
  0.0311267619693239538D+00, &
  0.0475792558412465455D+00, &
  0.062314485627767105D+00,  &
  0.0747979944082885623D+00, &
  0.0845782596975014622D+00, &
  0.0913017075224620001D+00, &
  0.0947253052275344315D+00, &
  0.0947253052275344315D+00, &
  0.0913017075224620001D+00, &
  0.0845782596975014622D+00, &
  0.0747979944082885623D+00, &
  0.062314485627767105D+00,  &
  0.0475792558412465455D+00, &
  0.0311267619693239538D+00, &
  0.0135762297058759553D+00 /)
  real ( kind = rk ), dimension ( 17 ) :: x17 = (/ &
  0.00471226234279131795D+00, &
  0.024662239115616158D+00,   &
  0.0598804231365070994D+00,  &
  0.109242998051599316D+00,   &
  0.171164420391654692D+00,   &
  0.243654731456761531D+00,   &
  0.324384118273061794D+00,   &
  0.410757909252076114D+00,   &
  0.500000000000000000D+00,   &
  0.589242090747923886D+00,   &
  0.675615881726938206D+00,   &
  0.756345268543238469D+00,   &
  0.828835579608345308D+00,   &
  0.890757001948400684D+00,   &
  0.940119576863492901D+00,   &
  0.975337760884383842D+00,   &
  0.995287737657208682D+00 /)
  real ( kind = rk ), dimension ( 17 ) :: w17 = (/ &
  0.01207415143427314D+00,   &
  0.0277297646869936118D+00, &
  0.0425180741585896443D+00, &
  0.0559419235967020534D+00, &
  0.0675681842342628902D+00, &
  0.0770228805384053083D+00, &
  0.0840020510782251428D+00, &
  0.0882813526834964474D+00, &
  0.0897232351781034193D+00, &
  0.0882813526834964474D+00, &
  0.0840020510782251428D+00, &
  0.0770228805384053083D+00, &
  0.0675681842342628902D+00, &
  0.0559419235967020534D+00, &
  0.0425180741585896443D+00, &
  0.0277297646869936118D+00, &
  0.01207415143427314D+00 /)
  real ( kind = rk ), dimension ( 18 ) :: x18 = (/ &
  0.00421741578953449547D+00, &
  0.022088025214301199D+00,   &
  0.0536987667512220934D+00,  &
  0.0981475205137384288D+00,  &
  0.154156478469823388D+00,   &
  0.22011458446302623D+00,    &
  0.294124419268578685D+00,   &
  0.374056887154247231D+00,   &
  0.457612493479132354D+00,   &
  0.542387506520867646D+00,   &
  0.625943112845752769D+00,   &
  0.705875580731421315D+00,   &
  0.77988541553697377D+00,    &
  0.845843521530176612D+00,   &
  0.901852479486261571D+00,   &
  0.946301233248777907D+00,   &
  0.977911974785698801D+00,   &
  0.995782584210465505D+00 /)
  real ( kind = rk ), dimension ( 18 ) :: w18 = (/ &
  0.0108080067632407191D+00, &
  0.0248572744474849679D+00, &
  0.0382128651274446646D+00, &
  0.0504710220531437159D+00, &
  0.061277603355739306D+00,  &
  0.0703214573353254518D+00, &
  0.0773423375631328014D+00, &
  0.0821382418729165037D+00, &
  0.084571191481571939D+00,  &
  0.084571191481571939D+00,  &
  0.0821382418729165037D+00, &
  0.0773423375631328014D+00, &
  0.0703214573353254518D+00, &
  0.061277603355739306D+00,  &
  0.0504710220531437159D+00, &
  0.0382128651274446646D+00, &
  0.0248572744474849679D+00, &
  0.0108080067632407191D+00 /)
  real ( kind = rk ), dimension ( 19 ) :: x19 = (/ &
  0.00379657807820787951D+00, &
  0.0198959239325849913D+00,  &
  0.048422048192590994D+00,   &
  0.0886426717314285906D+00,  &
  0.139516911332385307D+00,   &
  0.19972734766915945D+00,    &
  0.267714629312019503D+00,   &
  0.341717950018185057D+00,   &
  0.419820677179887358D+00,   &
  0.500000000000000000D+00,   &
  0.580179322820112642D+00,   &
  0.658282049981814943D+00,   &
  0.732285370687980497D+00,   &
  0.80027265233084055D+00,    &
  0.860483088667614693D+00,   &
  0.911357328268571409D+00,   &
  0.951577951807409006D+00,   &
  0.980104076067415009D+00,   &
  0.99620342192179212D+00 /)
  real ( kind = rk ), dimension ( 19 ) :: w19 = (/ &
  0.00973089411486243415D+00, &
  0.0224071133828498206D+00,  &
  0.0345222713688206687D+00,  &
  0.0457450108112251244D+00,  &
  0.0557833227736671128D+00,  &
  0.064376981269668232D+00,   &
  0.0713033510868034126D+00,  &
  0.0763830210329299597D+00,  &
  0.0794844216969773365D+00,  &
  0.0805272249243919463D+00,  &
  0.0794844216969773365D+00,  &
  0.0763830210329299597D+00,  &
  0.0713033510868034126D+00,  &
  0.064376981269668232D+00,   &
  0.0557833227736671128D+00,  &
  0.0457450108112251244D+00,  &
  0.0345222713688206687D+00,  &
  0.0224071133828498206D+00,  &
  0.00973089411486243415D+00 /)
  real ( kind = rk ), dimension ( 20 ) :: x20 = (/ &
  0.00343570040745255767D+00, &
  0.0180140363610430398D+00,  &
  0.0438827858743370269D+00,  &
  0.0804415140888905533D+00,  &
  0.126834046769924602D+00,   &
  0.181973159636742432D+00,   &
  0.244566499024586381D+00,   &
  0.313146955642290226D+00,   &
  0.386107074429177466D+00,   &
  0.461736739433251331D+00,   &
  0.538263260566748669D+00,   &
  0.613892925570822534D+00,   &
  0.686853044357709774D+00,   &
  0.755433500975413619D+00,   &
  0.818026840363257568D+00,   &
  0.873165953230075398D+00,   &
  0.919558485911109447D+00,   &
  0.956117214125662973D+00,   &
  0.98198596363895696D+00,    &
  0.996564299592547442D+00 /)
  real ( kind = rk ), dimension ( 20 ) :: w20 = (/ &
  0.0088070035695753026D+00, &
  0.0203007149001935561D+00, &
  0.0313360241670545686D+00, &
  0.0416383707883524329D+00, &
  0.0509650599086203179D+00, &
  0.0590972659807592476D+00, &
  0.0658443192245883463D+00, &
  0.0710480546591911871D+00, &
  0.0745864932363019956D+00, &
  0.0763766935653631129D+00, &
  0.0763766935653631129D+00, &
  0.0745864932363019956D+00, &
  0.0710480546591911871D+00, &
  0.0658443192245883463D+00, &
  0.0590972659807592476D+00, &
  0.0509650599086203179D+00, &
  0.0416383707883524329D+00, &
  0.0313360241670545686D+00, &
  0.0203007149001935561D+00, &
  0.0088070035695753026D+00 /)
  real ( kind = rk ), dimension ( 21 ) :: x21 = (/ &
  0.00312391468980521836D+00, &
  0.0163865807168468436D+00,  &
  0.0399503329247996586D+00,  &
  0.0733183177083414073D+00,  &
  0.115780018262161111D+00,   &
  0.166430597901293886D+00,   &
  0.224190582056390086D+00,   &
  0.287828939896280556D+00,   &
  0.355989341598799469D+00,   &
  0.427219072919552412D+00,   &
  0.500000000000000000D+00,   &
  0.572780927080447588D+00,   &
  0.644010658401200531D+00,   &
  0.712171060103719444D+00,   &
  0.775809417943609914D+00,   &
  0.833569402098706114D+00,   &
  0.884219981737838889D+00,   &
  0.926681682291658593D+00,   &
  0.960049667075200341D+00,   &
  0.983613419283153156D+00,   &
  0.996876085310194782D+00 /)
  real ( kind = rk ), dimension ( 21 ) :: w21 = (/ &
  0.00800861412888644909D+00, &
  0.018476894885426285D+00,   &
  0.0285672127134286406D+00,  &
  0.0380500568141897075D+00,  &
  0.0467222117280169935D+00,  &
  0.0543986495835743558D+00,  &
  0.06091570802686435D+00,    &
  0.0661344693166688452D+00,  &
  0.0699436973955366581D+00,  &
  0.0722622019949851341D+00,  &
  0.073040566824845346D+00,   &
  0.0722622019949851341D+00,  &
  0.0699436973955366581D+00,  &
  0.0661344693166688452D+00,  &
  0.06091570802686435D+00,    &
  0.0543986495835743558D+00,  &
  0.0467222117280169935D+00,  &
  0.0380500568141897075D+00,  &
  0.0285672127134286406D+00,  &
  0.018476894885426285D+00,   &
  0.00800861412888644909D+00 /)
  real ( kind = rk ), dimension ( 22 ) :: x22 = (/ &
  0.0028527072588003799D+00, &
  0.0149697510822857094D+00, &
  0.036521613906413064D+00,  &
  0.0670937111398499653D+00, &
  0.106091597010395944D+00,  &
  0.152756368406658627D+00,  &
  0.206179798246544199D+00,  &
  0.265322081006621469D+00,  &
  0.329032089553957907D+00,  &
  0.396069786655889322D+00,  &
  0.465130363340138908D+00,  &
  0.534869636659861092D+00,  &
  0.603930213344110678D+00,  &
  0.670967910446042093D+00,  &
  0.734677918993378531D+00,  &
  0.793820201753455801D+00,  &
  0.847243631593341373D+00,  &
  0.893908402989604056D+00,  &
  0.932906288860150035D+00,  &
  0.963478386093586936D+00,  &
  0.985030248917714291D+00,  &
  0.99714729274119962D+00 /)
  real ( kind = rk ), dimension ( 22 ) :: w22 = (/ &
  0.00731399764913532799D+00, &
  0.0168874507924071104D+00,  &
  0.0261466675763416916D+00,  &
  0.0348982342122602998D+00,  &
  0.0429708031085339753D+00,  &
  0.0502070722214406004D+00,  &
  0.0564661480402697119D+00,  &
  0.0616261884052562506D+00,  &
  0.0655867523935313168D+00,  &
  0.0682707491730076971D+00,  &
  0.0696259364278161291D+00,  &
  0.0696259364278161291D+00,  &
  0.0682707491730076971D+00,  &
  0.0655867523935313168D+00,  &
  0.0616261884052562506D+00,  &
  0.0564661480402697119D+00,  &
  0.0502070722214406004D+00,  &
  0.0429708031085339753D+00,  &
  0.0348982342122602998D+00,  &
  0.0261466675763416916D+00,  &
  0.0168874507924071104D+00,  &
  0.00731399764913532799D+00 /)
  real ( kind = rk ), dimension ( 23 ) :: x23 = (/ &
  0.00261533250122392147D+00, &
  0.013728764390942394D+00,   &
  0.0335144565869919253D+00,  &
  0.061623820864779244D+00,   &
  0.0975557991905799948D+00,  &
  0.140669318434024859D+00,   &
  0.190195062118176939D+00,   &
  0.245249261076996294D+00,   &
  0.304849480984854537D+00,   &
  0.367932159514827495D+00,   &
  0.433371587850766904D+00,   &
  0.500000000000000000D+00,   &
  0.566628412149233096D+00,   &
  0.632067840485172505D+00,   &
  0.695150519015145463D+00,   &
  0.754750738923003706D+00,   &
  0.809804937881823061D+00,   &
  0.859330681565975141D+00,   &
  0.902444200809420005D+00,   &
  0.938376179135220756D+00,   &
  0.966485543413008075D+00,   &
  0.986271235609057606D+00,   &
  0.997384667498776079D+00 /)
  real ( kind = rk ), dimension ( 23 ) :: w23 = (/ &
  0.0067059297435702412D+00,  &
  0.015494002928489686D+00,   &
  0.0240188358655423692D+00,  &
  0.0321162107042629943D+00,  &
  0.0396407058883595509D+00,  &
  0.0464578830300175633D+00,  &
  0.052446045732270824D+00,   &
  0.057498320111205814D+00,   &
  0.0615245421533648154D+00,  &
  0.06445286109404115D+00,    &
  0.0662310197023484037D+00,  &
  0.0668272860930531759D+00,  &
  0.0662310197023484037D+00,  &
  0.06445286109404115D+00,    &
  0.0615245421533648154D+00,  &
  0.057498320111205814D+00,   &
  0.052446045732270824D+00,   &
  0.0464578830300175633D+00,  &
  0.0396407058883595509D+00,  &
  0.0321162107042629943D+00,  &
  0.0240188358655423692D+00,  &
  0.015494002928489686D+00,   &
  0.0067059297435702412D+00 /)
  real ( kind = rk ), dimension ( 24 ) :: x24 = (/ &
  0.00240639000148934468D+00, &
  0.0126357220143452631D+00,  &
  0.0308627239986336566D+00,  &
  0.0567922364977995198D+00,  &
  0.0899990070130485265D+00,  &
  0.129937904210722821D+00,   &
  0.175953174031512227D+00,   &
  0.227289264305580163D+00,   &
  0.283103246186977464D+00,   &
  0.342478660151918302D+00,   &
  0.404440566263191803D+00,   &
  0.467971553568697241D+00,   &
  0.532028446431302759D+00,   &
  0.595559433736808197D+00,   &
  0.657521339848081698D+00,   &
  0.716896753813022536D+00,   &
  0.772710735694419837D+00,   &
  0.824046825968487773D+00,   &
  0.870062095789277179D+00,   &
  0.910000992986951474D+00,   &
  0.94320776350220048D+00,    &
  0.969137276001366343D+00,   &
  0.987364277985654737D+00,   &
  0.997593609998510655D+00 /)
  real ( kind = rk ), dimension ( 24 ) :: w24 = (/ &
  0.00617061489999283508D+00,&
  0.014265694314466934D+00,  &
  0.0221387194087098796D+00, &
  0.0296492924577183847D+00, &
  0.0366732407055402054D+00, &
  0.0430950807659766927D+00, &
  0.0488093260520570393D+00, &
  0.0537221350579829143D+00, &
  0.0577528340268628829D+00, &
  0.0608352364639017928D+00, &
  0.0629187281734143178D+00, &
  0.0639690976733762462D+00, &
  0.0639690976733762462D+00, &
  0.0629187281734143178D+00, &
  0.0608352364639017928D+00, &
  0.0577528340268628829D+00, &
  0.0537221350579829143D+00, &
  0.0488093260520570393D+00, &
  0.0430950807659766927D+00, &
  0.0366732407055402054D+00, &
  0.0296492924577183847D+00, &
  0.0221387194087098796D+00, &
  0.014265694314466934D+00,  &
  0.00617061489999283508D+00 /)
  real ( kind = rk ), dimension ( 25 ) :: x25 = (/ &
  0.00222151510475088187D+00,&
  0.0116680392702412927D+00, &
  0.0285127143855128384D+00, &
  0.052504001060862393D+00,  &
  0.0832786856195830705D+00, &
  0.120370368481321099D+00,  &
  0.163216815763265854D+00,  &
  0.211168534879388581D+00,  &
  0.263498634277142485D+00,  &
  0.319413847095306069D+00,  &
  0.378066558139505737D+00,  &
  0.438567653694644788D+00,  &
  0.500000000000000000D+00,  &
  0.561432346305355212D+00,  &
  0.621933441860494263D+00,  &
  0.680586152904693931D+00,  &
  0.736501365722857515D+00,  &
  0.788831465120611419D+00,  &
  0.836783184236734146D+00,  &
  0.879629631518678901D+00,  &
  0.91672131438041693D+00,   &
  0.947495998939137607D+00,  &
  0.971487285614487162D+00,  &
  0.988331960729758707D+00,  &
  0.997778484895249118D+00 /)
  real ( kind = rk ), dimension ( 25 ) :: w25 = (/ &
  0.00569689925051255347D+00, &
  0.0131774933075161083D+00,  &
  0.0204695783506531476D+00,  &
  0.0274523479879176906D+00,  &
  0.0340191669061785454D+00,  &
  0.0400703501675005319D+00,  &
  0.0455141309914819034D+00,  &
  0.0502679745335253628D+00,  &
  0.0542598122371318672D+00,  &
  0.0574291295728558623D+00,  &
  0.0597278817678924615D+00,  &
  0.0611212214951551217D+00,  &
  0.061588026863357799D+00,   &
  0.0611212214951551217D+00,  &
  0.0597278817678924615D+00,  &
  0.0574291295728558623D+00,  &
  0.0542598122371318672D+00,  &
  0.0502679745335253628D+00,  &
  0.0455141309914819034D+00,  &
  0.0400703501675005319D+00,  &
  0.0340191669061785454D+00,  &
  0.0274523479879176906D+00,  &
  0.0204695783506531476D+00,  &
  0.0131774933075161083D+00,  &
  0.00569689925051255347D+00 /)

  if ( n == 1 ) then
    call r8vec_copy ( n, x01, x )
    call r8vec_copy ( n, w01, w )
  else if ( n == 2 ) then
    call r8vec_copy ( n, x02, x )
    call r8vec_copy ( n, w02, w )
  else if ( n == 3 ) then
    call r8vec_copy ( n, x03, x )
    call r8vec_copy ( n, w03, w )
  else if ( n == 4 ) then
    call r8vec_copy ( n, x04, x )
    call r8vec_copy ( n, w04, w )
  else if ( n == 5 ) then
    call r8vec_copy ( n, x05, x )
    call r8vec_copy ( n, w05, w )
  else if ( n == 6 ) then
    call r8vec_copy ( n, x06, x )
    call r8vec_copy ( n, w06, w )
  else if ( n == 7 ) then
    call r8vec_copy ( n, x07, x )
    call r8vec_copy ( n, w07, w )
  else if ( n == 8 ) then
    call r8vec_copy ( n, x08, x )
    call r8vec_copy ( n, w08, w )
  else if ( n == 9 ) then
    call r8vec_copy ( n, x09, x )
    call r8vec_copy ( n, w09, w )
  else if ( n == 10 ) then
    call r8vec_copy ( n, x10, x )
    call r8vec_copy ( n, w10, w )
  else if ( n == 11 ) then
    call r8vec_copy ( n, x11, x )
    call r8vec_copy ( n, w11, w )
  else if ( n == 12 ) then
    call r8vec_copy ( n, x12, x )
    call r8vec_copy ( n, w12, w )
  else if ( n == 13 ) then
    call r8vec_copy ( n, x13, x )
    call r8vec_copy ( n, w13, w )
  else if ( n == 14 ) then
    call r8vec_copy ( n, x14, x )
    call r8vec_copy ( n, w14, w )
  else if ( n == 15 ) then
    call r8vec_copy ( n, x15, x )
    call r8vec_copy ( n, w15, w )
  else if ( n == 16 ) then
    call r8vec_copy ( n, x16, x )
    call r8vec_copy ( n, w16, w )
  else if ( n == 17 ) then
    call r8vec_copy ( n, x17, x )
    call r8vec_copy ( n, w17, w )
  else if ( n == 18 ) then
    call r8vec_copy ( n, x18, x )
    call r8vec_copy ( n, w18, w )
  else if ( n == 19 ) then
    call r8vec_copy ( n, x19, x )
    call r8vec_copy ( n, w19, w )
  else if ( n == 20 ) then
    call r8vec_copy ( n, x20, x )
    call r8vec_copy ( n, w20, w )
  else if ( n == 21 ) then
    call r8vec_copy ( n, x21, x )
    call r8vec_copy ( n, w21, w )
  else if ( n == 22 ) then
    call r8vec_copy ( n, x22, x )
    call r8vec_copy ( n, w22, w )
  else if ( n == 23 ) then
    call r8vec_copy ( n, x23, x )
    call r8vec_copy ( n, w23, w )
  else if ( n == 24 ) then
    call r8vec_copy ( n, x24, x )
    call r8vec_copy ( n, w24, w )
  else if ( n == 25 ) then
    call r8vec_copy ( n, x25, x )
    call r8vec_copy ( n, w25, w )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GQU - Fatal error!'
    write ( *, '(a)' ) '  Value of N must be between 1 and 25.'
    stop 1
  end if

  return
end
subroutine gqu_order ( l, n )

!*****************************************************************************80
!
!! GQU_ORDER computes the order of a GQU rule from the level.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 June 2012
!
!  Author:
!
!    John Burkardt.
!
!  Parameters:
!
!    Input, integer L, the level of the rule.  
!    1 <= L <= 25.
!
!    Output, integer N, the order of the rule.
!
  implicit none

  integer l
  integer n

  if ( l < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GQU_ORDER - Fatal error!'
    write ( *, '(a)' ) '  1 <= L required.'
    write ( *, '(a,i4)' ) '  Input L = ', l
    stop 1
  else if ( l <= 25 ) then
    n = l
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GQU_ORDER - Fatal error!'
    write ( *, '(a)' ) '  L <= 25 required.'
    write ( *, '(a,i4)' ) '  Input L = ', l
    stop 1
  end if

  return
end
function i4_choose ( n, k )

!*****************************************************************************80
!
!! I4_CHOOSE computes the binomial coefficient C(N,K) as an I4.
!
!  Discussion:
!
!    The value is calculated in such a way as to avoid overflow and
!    roundoff.  The calculation is done in integer arithmetic.
!
!    The formula used is:
!
!      C(N,K) = N! / ( K! * (N-K)! )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    ML Wolfson, HV Wright,
!    Algorithm 160:
!    Combinatorial of M Things Taken N at a Time,
!    Communications of the ACM,
!    Volume 6, Number 4, April 1963, page 161.
!
!  Parameters:
!
!    Input, integer N, K, are the values of N and K.
!
!    Output, integer I4_CHOOSE, the number of combinations of N
!    things taken K at a time.
!
  implicit none

  integer i
  integer i4_choose
  integer k
  integer mn
  integer mx
  integer n
  integer value

  mn = min ( k, n - k )

  if ( mn < 0 ) then

    value = 0

  else if ( mn == 0 ) then

    value = 1

  else

    mx = max ( k, n - k )
    value = mx + 1

    do i = 2, mn
      value = ( value * ( mx + i ) ) / i
    end do

  end if

  i4_choose = value

  return
end
function i4_factorial2 ( n )

!*****************************************************************************80
!
!! I4_FACTORIAL2 computes the double factorial function.
!
!  Discussion:
!
!    The formula is:
!
!      FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
!                      = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
!
!  Example:
!
!     N    Factorial2(N)
!
!     0     1
!     1     1
!     2     2
!     3     3
!     4     8
!     5    15
!     6    48
!     7   105
!     8   384
!     9   945
!    10  3840
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
!    Input, integer N, the argument of the double factorial 
!    function.  If N is less than 1, I4_FACTORIAL2 is returned as 1.
!
!    Output, integer I4_FACTORIAL2, the value of the function..
!
  implicit none

  integer i4_factorial2
  integer n
  integer n_copy

  if ( n < 1 ) then
    i4_factorial2 = 1
    return
  end if

  n_copy = n
  i4_factorial2 = 1

  do while ( 1 < n_copy )
    i4_factorial2 = i4_factorial2 * n_copy
    n_copy = n_copy - 2
  end do

  return
end
function i4_mop ( i )

!*****************************************************************************80
!
!! I4_MOP returns the I-th power of -1 as an I4 value.
!
!  Discussion:
!
!    An I4 is an integer value.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, the power of -1.
!
!    Output, integer I4_MOP, the I-th power of -1.
!
  implicit none

  integer i
  integer i4_mop

  if ( mod ( i, 2 ) == 0 ) then
    i4_mop = 1
  else
    i4_mop = -1
  end if

  return
end
subroutine i4mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! I4MAT_PRINT prints an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4 values.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 June 2003
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
!    Input, integer A(M,N), the matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer m
  integer n

  integer a(m,n)
  integer ihi
  integer ilo
  integer jhi
  integer jlo
  character ( len = * ) title

  ilo = 1
  ihi = m
  jlo = 1
  jhi = n

  call i4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

  return
end
subroutine i4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! I4MAT_PRINT_SOME prints some of an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4 values.
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
!    Input, integer A(M,N), an M by N matrix to be printed.
!
!    Input, integer ILO, JLO, the first row and column to print.
!
!    Input, integer IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: incx = 10
  integer m
  integer n

  integer a(m,n)
  character ( len = 8 ) ctemp(incx)
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
      write ( ctemp(j2), '(i8)' ) j
    end do

    write ( *, '(''  Col '',10a8)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        write ( ctemp(j2), '(i8)' ) a(i,j)

      end do

      write ( *, '(i5,a,10a8)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine i4vec_cum0 ( n, a, a_cum )

!*****************************************************************************80
!
!! I4VEC_CUM0 computes the cumulutive sum of the entries of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    This routine returns a vector of length N+1, with the first value
!    being 0.
!
!  Example:
!
!    Input:
!
!      A = (/ 1, 2, 3, 4 /)
!
!    Output:
!
!      A_CUM = (/ 0, 1, 3, 6, 10 /)
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 December 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, integer A(N), the vector to be summed.
!
!    Output, integer A_CUM(0:N), the cumulative sum of the
!    entries of A.
!
  implicit none

  integer n

  integer a(n)
  integer a_cum(0:n)
  integer i

  a_cum(0) = 0

  do i = 1, n
    a_cum(i) = a_cum(i-1) + a(i)
  end do

  return
end
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 May 2010
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

  integer n

  integer a(n)
  integer i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,a,2x,i12)' ) i, ':', a(i)
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

  integer n

  integer a(n)
  integer i4vec_product

  i4vec_product = product ( a(1:n) )

  return
end
function i4vec_sum ( n, a )

!*****************************************************************************80
!
!! I4VEC_SUM returns the sum of the entries of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    In FORTRAN90, this facility is offered by the built in
!    SUM function:
!
!      I4VEC_SUM ( N, A ) = SUM ( A(1:N) )
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
!    Output, integer I4VEC_SUM, the sum of the entries.
!
  implicit none

  integer n

  integer a(n)
  integer i4vec_sum

  i4vec_sum = sum ( a(1:n) )

  return
end
subroutine i4vec_transpose_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_TRANSPOSE_PRINT prints an I4VEC "transposed".
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Example:
!
!    A = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 /)
!    TITLE = 'My vector:  '
!
!    My vector:
!
!        1    2    3    4    5
!        6    7    8    9   10
!       11
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 April 2005
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

  integer n

  integer a(n)
  integer ihi
  integer ilo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do ilo = 1, n, 5
    ihi = min ( ilo + 5 - 1, n )
    write ( *, '(5i12)' ) a(ilo:ihi)
  end do

  return
end
subroutine kpn ( n, x, w )

!*****************************************************************************80
!
!! KPN provides data for Kronrod-Patterson quadrature with a normal weight.
!
!  Discussion:
!
!    This data assumes integration over the interval (-oo,+oo) with 
!    weight function w(x) = exp(-x*x/2)/sqrt(2*pi).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 December 2012
!
!  Author:
!
!    Original MATLAB version by Florian Heiss, Viktor Winschel.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Florian Heiss, Viktor Winschel,
!    Likelihood approximation by numerical integration on sparse grids,
!    Journal of Econometrics,
!    Volume 144, 2008, pages 62-80.
!
!    Alan Genz, Bradley Keister,
!    Fully symmetric interpolatory rules for multiple integrals
!    over infinite regions with Gaussian weight,
!    Journal of Computational and Applied Mathematics,
!    Volume 71, 1996, pages 299-309.
!
!    Thomas Patterson,
!    The optimal addition of points to quadrature formulae,
!    Mathematics of Computation,
!    Volume 22, Number 104, October 1968, pages 847-856.
!
!  Parameters:
!
!    Input, integer N, the order of the rule.
!
!    Output, real ( kind = rk ) X(N), the nodes.
!
!    Output, real ( kind = rk ) W(N), the weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n)
  real ( kind = rk ), dimension ( 1 ) :: x01 = (/ &
    0.0000000000000000D+00 /)
  real ( kind = rk ), dimension ( 1 ) :: w01 = (/ &
    1.0000000000000000D+00/)
  real ( kind = rk ), dimension ( 3 ) :: x03 = (/ &
   -1.73205080756887719D+00,  &
    0.000000000000000000D+00, &
    1.73205080756887719D+00 /)
  real ( kind = rk ), dimension ( 3 ) :: w03 = (/ &
    0.166666666666666657D+00, &
    0.66666666666666663D+00,  &
    0.166666666666666657D+00/)
  real ( kind = rk ), dimension ( 7 ) :: x07 = (/ &
   -4.18495601767273229D+00,  &
   -1.73205080756887719D+00,  &
   -0.741095349994540853D+00, & 
    0.00000000000000000D+00,  &
    0.741095349994540853D+00, &
    1.73205080756887719D+00,  & 
    4.18495601767273229D+00 /)
  real ( kind = rk ), dimension ( 7 ) :: w07 = (/ &
    0.000695684158369139867D+00, &
    0.138553274729749237D+00,    &
    0.13137860698313561D+00,     & 
    0.458744868257491889D+00,    &
    0.13137860698313561D+00,     &
    0.138553274729749237D+00,    & 
    0.000695684158369139867D+00  /)
  real ( kind = rk ), dimension ( 9 ) :: x09 = (/ &
   -4.18495601767273229D+00,  &
   -2.86127957605705818D+00,  &
   -1.73205080756887719D+00,  &
   -0.741095349994540853D+00, &
    0.00000000000000000D+00,  &
    0.741095349994540853D+00, &
    1.73205080756887719D+00,  &
    2.86127957605705818D+00,  &
    4.18495601767273229D+00 /)
  real ( kind = rk ), dimension ( 9 ) :: w09 = (/ &
    9.42694575565174701D-05,    &
    0.00799632547089352934D+00, &
    0.0948509485094851251D+00,  &
    0.270074329577937755D+00,   &
    0.253968253968254065D+00,   &
    0.270074329577937755D+00,   & 
    0.0948509485094851251D+00,  &
    0.00799632547089352934D+00, &
    9.42694575565174701D-05 /)
  real ( kind = rk ), dimension ( 17 ) :: x17 = (/ &
   -6.36339449433636961D+00,   &
   -5.18701603991365623D+00,   &
   -4.18495601767273229D+00,   & 
   -2.86127957605705818D+00,   &
   -2.59608311504920231D+00,   &
   -1.73205080756887719D+00,   & 
   -1.23042363402730603D+00,   &
   -0.741095349994540853D+00,  &
    0.00000000000000000D+00,   &
    0.741095349994540853D+00,  &
    1.23042363402730603D+00,   &
    1.73205080756887719D+00,   & 
    2.59608311504920231D+00,   &
    2.86127957605705818D+00,   &
    4.18495601767273229D+00,   & 
    5.18701603991365623D+00,   &
    6.36339449433636961D+00 /)
  real ( kind = rk ), dimension ( 17 ) :: w17 = (/ &
    2.11364995054242569D-08,     &
   -8.20492075415092169D-07,     &
    0.000105637836154169414D+00, & 
    0.00703348023782790748D+00,  &
    0.0019656770938777492D+00,   &
    0.0886810021520280101D+00,   & 
    0.0141926548264493645D+00,   &
    0.254561232041712215D+00,    &
    0.266922230335053023D+00,    &
    0.254561232041712215D+00,    &
    0.0141926548264493645D+00,   &
    0.0886810021520280101D+00,   & 
    0.0019656770938777492D+00,   &
    0.00703348023782790748D+00,  &
    0.000105637836154169414D+00, & 
   -8.20492075415092169D-07,     &
    2.11364995054242569D-08 /)
  real ( kind = rk ), dimension ( 19 ) :: x19 = (/ &
   -6.36339449433636961D+00,  &
   -5.18701603991365623D+00,  &
   -4.18495601767273229D+00,  & 
   -3.20533379449919442D+00,  &
   -2.86127957605705818D+00,  &
   -2.59608311504920231D+00,  & 
   -1.73205080756887719D+00,  &
   -1.23042363402730603D+00,  &
   -0.741095349994540853D+00, &
    0.0000000000000000D+00,   &
    0.741095349994540853D+00, &
    1.23042363402730603D+00,  &
    1.73205080756887719D+00,  &
    2.59608311504920231D+00,  &
    2.86127957605705818D+00,  &
    3.20533379449919442D+00,  &
    4.18495601767273229D+00,  &
    5.18701603991365623D+00,  &
    6.36339449433636961D+00 /)
  real ( kind = rk ), dimension ( 19 ) :: w19 = (/ &
    8.62968460222986318D-10,    &
    6.09480873146898402D-07,    &
    6.01233694598479965D-05,    &
    0.00288488043650675591D+00, &
   -0.00633722479337375712D+00, &
    0.0180852342547984622D+00,  &
    0.0640960546868076103D+00,  &
    0.0611517301252477163D+00,  &
    0.208324991649608771D+00,   &
    0.303467199854206227D+00,   &
    0.208324991649608771D+00,   &
    0.0611517301252477163D+00,  &
    0.0640960546868076103D+00,  &
    0.0180852342547984622D+00,  &
   -0.00633722479337375712D+00, &
    0.00288488043650675591D+00, &
    6.01233694598479965D-05,    &
    6.09480873146898402D-07,    &
    8.62968460222986318D-10 /)
  real ( kind = rk ), dimension ( 31 ) :: x31 = (/ &
   -9.0169397898903032D+00,     &
   -7.98077179859056063D+00,    &
   -7.12210670080461661D+00,    &
   -6.36339449433636961D+00,    &
   -5.18701603991365623D+00,    &
   -4.18495601767273229D+00,    &
   -3.63531851903727832D+00,    &
   -3.20533379449919442D+00,    &
   -2.86127957605705818D+00,    &
   -2.59608311504920231D+00,    &
   -2.23362606167694189D+00,    &
   -1.73205080756887719D+00,    &
   -1.23042363402730603D+00,    &
   -0.741095349994540853D+00,   &
   -0.248992297579960609D+00,   &
    0.00000000000000000D+00,    &
    0.248992297579960609D+00,   &
    0.741095349994540853D+00,   &
    1.23042363402730603D+00,    &
    1.73205080756887719D+00,    &
    2.23362606167694189D+00,    &
    2.59608311504920231D+00,    &
    2.86127957605705818D+00,    &
    3.20533379449919442D+00,    &
    3.63531851903727832D+00,    &
    4.18495601767273229D+00,    &
    5.18701603991365623D+00,    &
    6.36339449433636961D+00,    &
    7.12210670080461661D+00,    &
    7.98077179859056063D+00,    &
    9.0169397898903032D+00 /)
  real ( kind = rk ), dimension ( 31 ) :: w31 = (/ &
    1.26184642808151181D-15,    &
   -1.4840835740298868D-13,     &
    5.11580531055042083D-12,    &
    7.92982678648693382D-10,    &
    6.14358432326179133D-07,    &
    5.94749611639316215D-05,    &
    1.50442053909142189D-05,    &
    0.00272984304673340016D+00, &
   -0.00556100630683581572D+00, &
    0.0165924926989360101D+00,  &
    0.00176084755813180017D+00, &
    0.0617185325658671791D+00,  &
    0.0654173928360925611D+00,  &
    0.199688635117345498D+00,   &
    0.0281281015400331666D+00,  &
    0.25890005324151566D+00,    &
    0.0281281015400331666D+00,  &
    0.199688635117345498D+00,   &
    0.0654173928360925611D+00,  &
    0.0617185325658671791D+00,  &
    0.00176084755813180017D+00, &
    0.0165924926989360101D+00,  &
   -0.00556100630683581572D+00, &
    0.00272984304673340016D+00, &
    1.50442053909142189D-05,    &
    5.94749611639316215D-05,    &
    6.14358432326179133D-07,    &
    7.92982678648693382D-10,    &
    5.11580531055042083D-12,    &
   -1.4840835740298868D-13,     &
    1.26184642808151181D-15 /)
  real ( kind = rk ), dimension ( 33 ) :: x33 = (/ &
   -9.0169397898903032D+00,     &
   -7.98077179859056063D+00,    &
   -7.12210670080461661D+00,    &
   -6.36339449433636961D+00,    &
   -5.69817776848810986D+00,    &
   -5.18701603991365623D+00,    &
   -4.18495601767273229D+00,    &
   -3.63531851903727832D+00,    &
   -3.20533379449919442D+00,    &
   -2.86127957605705818D+00,    &
   -2.59608311504920231D+00,    &
   -2.23362606167694189D+00,    &
   -1.73205080756887719D+00,    &
   -1.23042363402730603D+00,    &
   -0.741095349994540853D+00,   &
   -0.248992297579960609D+00,   &
    0.00000000000000000D+00,    &
    0.248992297579960609D+00,   &
    0.741095349994540853D+00,   &
    1.23042363402730603D+00,    &
    1.73205080756887719D+00,    &
    2.23362606167694189D+00,    &
    2.59608311504920231D+00,    &
    2.86127957605705818D+00,    &
    3.20533379449919442D+00,    &
    3.63531851903727832D+00,    &
    4.18495601767273229D+00,    &
    5.18701603991365623D+00,    &
    5.69817776848810986D+00,    &
    6.36339449433636961D+00,    &
    7.12210670080461661D+00,    &
    7.98077179859056063D+00,    &
    9.0169397898903032D+00 /)
  real ( kind = rk ), dimension ( 33 ) :: w33 = (/ &
   -9.93139132868224651D-16,    &
    2.66406251662316506D-13,    &
   -1.93413050008809555D-11,    &
    1.5542195992782658D-09,     &
   -1.34860173485429301D-08,    &
    6.90862611791137378D-07,    &
    5.56911589810814793D-05,    &
    8.32360452957667447D-05,    &
    0.00212022595595963252D+00, &
   -0.00277121890077892431D+00, &
    0.01152924706539879D+00,    &
    0.00735301102049550764D+00, &
    0.0546775561434630422D+00,  &
    0.0774436027462994808D+00,  &
    0.176075987415714591D+00,   &
    0.103876871255742839D+00,   &
    0.139110222363380387D+00,   &
    0.103876871255742839D+00,   &
    0.176075987415714591D+00,   &
    0.0774436027462994808D+00,  &
    0.0546775561434630422D+00,  &
    0.00735301102049550764D+00, &
    0.01152924706539879D+00,    &
   -0.00277121890077892431D+00, &
    0.00212022595595963252D+00, &
    8.32360452957667447D-05,    &
    5.56911589810814793D-05,    &
    6.90862611791137378D-07,    &
   -1.34860173485429301D-08,    &
    1.5542195992782658D-09,     &
   -1.93413050008809555D-11,    &
    2.66406251662316506D-13,    &
   -9.93139132868224651D-16 /)
  real ( kind = rk ), dimension ( 35 ) :: x35 = (/ &
   -9.0169397898903032D+00,     &
   -7.98077179859056063D+00,    &
   -7.12210670080461661D+00,    &
   -6.36339449433636961D+00,    &
   -5.69817776848810986D+00,    &
   -5.18701603991365623D+00,    &
   -4.73643308595229673D+00,    &
   -4.18495601767273229D+00,    &
   -3.63531851903727832D+00,    &
   -3.20533379449919442D+00,    &
   -2.86127957605705818D+00,    &
   -2.59608311504920231D+00,    &
   -2.23362606167694189D+00,    &
   -1.73205080756887719D+00,    &
   -1.23042363402730603D+00,    &
   -0.741095349994540853D+00,   &
   -0.248992297579960609D+00,   &
    0.00000000000000000D+00,    &
    0.248992297579960609D+00,   &
    0.741095349994540853D+00,   &
    1.23042363402730603D+00,    &
    1.73205080756887719D+00,    &
    2.23362606167694189D+00,    &
    2.59608311504920231D+00,    &
    2.86127957605705818D+00,    &
    3.20533379449919442D+00,    &
    3.63531851903727832D+00,    &
    4.18495601767273229D+00,    &
    4.73643308595229673D+00,    &
    5.18701603991365623D+00,    &
    5.69817776848810986D+00,    &
    6.36339449433636961D+00,    &
    7.12210670080461661D+00,    &
    7.98077179859056063D+00,    &
    9.0169397898903032D+00 /)
  real ( kind = rk ), dimension ( 35 ) :: w35 = (/ &
    1.05413265823340136D-18,     &
    5.45004126506381281D-15,     &
    3.09722235760629949D-12,     &
    4.60117603486559168D-10,     &
    2.13941944795610622D-08,     &
    2.46764213457981401D-07,     &
    2.73422068011878881D-06,     &
    3.57293481989753322D-05,     &
    0.000275242141167851312D+00, &
    0.000818953927502267349D+00, &
    0.00231134524035220713D+00,  &
    0.00315544626918755127D+00,  &
    0.015673473751851151D+00,    &
    0.0452736854651503914D+00,   &
    0.0923647267169863534D+00,   &
    0.148070831155215854D+00,    &
    0.191760115888044341D+00,    &
    0.000514894508069213769D+00, &
    0.191760115888044341D+00,    &
    0.148070831155215854D+00,    &
    0.0923647267169863534D+00,   &
    0.0452736854651503914D+00,   &
    0.015673473751851151D+00,    &
    0.00315544626918755127D+00,  &
    0.00231134524035220713D+00,  &
    0.000818953927502267349D+00, &
    0.000275242141167851312D+00, &
    3.57293481989753322D-05,     &
    2.73422068011878881D-06,     &
    2.46764213457981401D-07,     &
    2.13941944795610622D-08,     &
    4.60117603486559168D-10,     &
    3.09722235760629949D-12,     &
    5.45004126506381281D-15,     &
    1.05413265823340136D-18 /)

  if ( n == 1 ) then
    call r8vec_copy ( n, x01, x )
    call r8vec_copy ( n, w01, w )
  else if ( n == 3 ) then
    call r8vec_copy ( n, x03, x )
    call r8vec_copy ( n, w03, w )
  else if ( n == 7 ) then
    call r8vec_copy ( n, x07, x )
    call r8vec_copy ( n, w07, w )
  else if ( n == 9 ) then
    call r8vec_copy ( n, x09, x )
    call r8vec_copy ( n, w09, w )
  else if ( n == 17 ) then
    call r8vec_copy ( n, x17, x )
    call r8vec_copy ( n, w17, w )
  else if ( n == 19 ) then
    call r8vec_copy ( n, x19, x )
    call r8vec_copy ( n, w19, w )
  else if ( n == 31 ) then
    call r8vec_copy ( n, x31, x )
    call r8vec_copy ( n, w31, w )
  else if ( n == 33 ) then
    call r8vec_copy ( n, x33, x )
    call r8vec_copy ( n, w33, w )
  else if ( n == 35 ) then
    call r8vec_copy ( n, x35, x )
    call r8vec_copy ( n, w35, w )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KPN - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of N.'
    stop 1
  end if

  return
end
subroutine kpn_order ( l, n )

!*****************************************************************************80
!
!! KPN_ORDER computes the order of a KPN rule from the level.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 May 2012
!
!  Author:
!
!    John Burkardt.
!
!  Parameters:
!
!    Input, integer L, the level of the rule.  
!    1 <= L <= 25
!
!    Output, integer N, the order of the rule.
!
  implicit none

  integer l
  integer n

  if ( l < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KPN_ORDER - Fatal error!'
    write ( *, '(a)' ) '  1 <= L <= 25 required.'
    write ( *, '(a,i4)' ) '  Input L = ', l
    stop 1
  else if ( l == 1 ) then
    n = 1
  else if ( l <= 3 ) then
    n = 3
  else if ( l == 4 ) then
    n = 7
  else if ( l <= 8 ) then
    n = 9
  else if ( l == 9 ) then
    n = 17
  else if ( l <= 15 ) then
    n = 19
  else if ( l == 16 ) then
    n = 31
  else if ( l == 17 ) then
    n = 33
  else if ( l <= 25 ) then
    n = 35
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KPN_ORDER - Fatal error!'
    write ( *, '(a)' ) '  1 <= L <= 25 required.'
    write ( *, '(a,i4)' ) '  Input L = ', l
    stop 1
  end if

  return
end
subroutine kpu ( n, x, w )

!*****************************************************************************80
!
!! KPU provides data for Kronrod-Patterson quadrature with a uniform weight.
!
!  Discussion:
!
!    This data assumes integration over the interval [0,1] with 
!    weight function w(x) = 1.
!
!    This data was originally provided with 7 digit accuracy, over [0,1].
!    It has been replaced by higher accuracy data over [-1,+1], which
!    is adjusted to [0,1] before return.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 December 2012
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Florian Heiss, Viktor Winschel,
!    Likelihood approximation by numerical integration on sparse grids,
!    Journal of Econometrics,
!    Volume 144, 2008, pages 62-80.
!
!    Alan Genz, Bradley Keister,
!    Fully symmetric interpolatory rules for multiple integrals
!    over infinite regions with Gaussian weight,
!    Journal of Computational and Applied Mathematics,
!    Volume 71, 1996, pages 299-309.
!
!    Thomas Patterson,
!    The optimal addition of points to quadrature formulae,
!    Mathematics of Computation,
!    Volume 22, Number 104, October 1968, pages 847-856.
!
!  Parameters:
!
!    Input, integer N, the order of the rule.
!    Only 1, 3, 7, 15, 31 and 63 are legal input values for N.
!
!    Output, real ( kind = rk ) X(N), the nodes.
!
!    Output, real ( kind = rk ) W(N), the weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) d
  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n)
  real ( kind = rk ), dimension ( 1 ) :: x01 = (/ &
     0.0000000D+00 /)
  real ( kind = rk ), dimension ( 1 ) :: w01 = (/ &
     2.0000000D+00 /)
  real ( kind = rk ), dimension ( 3 ) :: x03 = (/ &
    -0.77459666924148337704D+00, &
     0.0D+00, &
     0.77459666924148337704D+00 /)
  real ( kind = rk ), dimension ( 3 ) :: w03 = (/ &
     0.555555555555555555556D+00, &
     0.888888888888888888889D+00, &
     0.555555555555555555556D+00 /)
  real ( kind = rk ), dimension ( 7 ) :: x07 = (/ &
    -0.96049126870802028342D+00, &
    -0.77459666924148337704D+00, &
    -0.43424374934680255800D+00, &
     0.0D+00, &
     0.43424374934680255800D+00, &
     0.77459666924148337704D+00, &
     0.96049126870802028342D+00 /)
  real ( kind = rk ), dimension ( 7 ) :: w07 = (/ &
     0.104656226026467265194D+00, &
     0.268488089868333440729D+00, &
     0.401397414775962222905D+00, &
     0.450916538658474142345D+00, &
     0.401397414775962222905D+00, &
     0.268488089868333440729D+00, &
     0.104656226026467265194D+00 /)
  real ( kind = rk ), dimension ( 15 ) :: x15 = (/ &
    -0.99383196321275502221D+00, &
    -0.96049126870802028342D+00, &
    -0.88845923287225699889D+00, &
    -0.77459666924148337704D+00, &
    -0.62110294673722640294D+00, &
    -0.43424374934680255800D+00, &
    -0.22338668642896688163D+00, &
     0.0D+00,                    &
     0.22338668642896688163D+00, &
     0.43424374934680255800D+00, &
     0.62110294673722640294D+00, &
     0.77459666924148337704D+00, &
     0.88845923287225699889D+00, &
     0.96049126870802028342D+00, &
     0.99383196321275502221D+00 /)
  real ( kind = rk ), dimension ( 15 ) :: w15 = (/ &
     0.0170017196299402603390D+00, &
     0.0516032829970797396969D+00, &
     0.0929271953151245376859D+00, &
     0.134415255243784220360D+00,  &
     0.171511909136391380787D+00,  &
     0.200628529376989021034D+00,  &
     0.219156858401587496404D+00,  &
     0.225510499798206687386D+00,  &
     0.219156858401587496404D+00,  &
     0.200628529376989021034D+00,  &
     0.171511909136391380787D+00,  &
     0.134415255243784220360D+00,  &
     0.0929271953151245376859D+00, &
     0.0516032829970797396969D+00, &
     0.0170017196299402603390D+00 /)
  real ( kind = rk ), dimension ( 31 ) :: x31 = (/ &
    -0.99909812496766759766D+00, &
    -0.99383196321275502221D+00, &
    -0.98153114955374010687D+00, &
    -0.96049126870802028342D+00, &
    -0.92965485742974005667D+00, &
    -0.88845923287225699889D+00, &
    -0.83672593816886873550D+00, &
    -0.77459666924148337704D+00, &
    -0.70249620649152707861D+00, &
    -0.62110294673722640294D+00, &
    -0.53131974364437562397D+00, &
    -0.43424374934680255800D+00, &
    -0.33113539325797683309D+00, &
    -0.22338668642896688163D+00, &
    -0.11248894313318662575D+00, &
     0.0D+00,                    &
     0.11248894313318662575D+00, &
     0.22338668642896688163D+00, &
     0.33113539325797683309D+00, &
     0.43424374934680255800D+00, &
     0.53131974364437562397D+00, &
     0.62110294673722640294D+00, &
     0.70249620649152707861D+00, &
     0.77459666924148337704D+00, &
     0.83672593816886873550D+00, &
     0.88845923287225699889D+00, &
     0.92965485742974005667D+00, &
     0.96049126870802028342D+00, &
     0.98153114955374010687D+00, &
     0.99383196321275502221D+00, &
     0.99909812496766759766D+00 /)
  real ( kind = rk ), dimension ( 31 ) :: w31 = (/ &
     0.00254478079156187441540D+00, &
     0.00843456573932110624631D+00, &
     0.0164460498543878109338D+00,  &
     0.0258075980961766535646D+00,  &
     0.0359571033071293220968D+00,  &
     0.0464628932617579865414D+00,  &
     0.0569795094941233574122D+00,  &
     0.0672077542959907035404D+00,  &
     0.0768796204990035310427D+00,  &
     0.0857559200499903511542D+00,  &
     0.0936271099812644736167D+00,  &
     0.100314278611795578771D+00,   &
     0.105669893580234809744D+00,   &
     0.109578421055924638237D+00,   &
     0.111956873020953456880D+00,   &
     0.112755256720768691607D+00,   &
     0.111956873020953456880D+00,   &
     0.109578421055924638237D+00,   &
     0.105669893580234809744D+00,   &
     0.100314278611795578771D+00,   &
     0.0936271099812644736167D+00,  &
     0.0857559200499903511542D+00,  &
     0.0768796204990035310427D+00,  &
     0.0672077542959907035404D+00,  &
     0.0569795094941233574122D+00,  &
     0.0464628932617579865414D+00,  &
     0.0359571033071293220968D+00,  &
     0.0258075980961766535646D+00,  &
     0.0164460498543878109338D+00,  &
     0.00843456573932110624631D+00, &
     0.00254478079156187441540D+00 /)
  real ( kind = rk ), dimension ( 63 ) :: x63 = (/ &
    -0.99987288812035761194D+00, &
    -0.99909812496766759766D+00, &
    -0.99720625937222195908D+00, &
    -0.99383196321275502221D+00, &
    -0.98868475754742947994D+00, &
    -0.98153114955374010687D+00, &
    -0.97218287474858179658D+00, &
    -0.96049126870802028342D+00, &
    -0.94634285837340290515D+00, &
    -0.92965485742974005667D+00, &
    -0.91037115695700429250D+00, &
    -0.88845923287225699889D+00, &
    -0.86390793819369047715D+00, &
    -0.83672593816886873550D+00, &
    -0.80694053195021761186D+00, &
    -0.77459666924148337704D+00, &
    -0.73975604435269475868D+00, &
    -0.70249620649152707861D+00, &
    -0.66290966002478059546D+00, &
    -0.62110294673722640294D+00, &
    -0.57719571005204581484D+00, &
    -0.53131974364437562397D+00, &
    -0.48361802694584102756D+00, &
    -0.43424374934680255800D+00, &
    -0.38335932419873034692D+00, &
    -0.33113539325797683309D+00, &
    -0.27774982202182431507D+00, &
    -0.22338668642896688163D+00, &
    -0.16823525155220746498D+00, &
    -0.11248894313318662575D+00, &
    -0.056344313046592789972D+00, &
     0.0D+00, &
     0.056344313046592789972D+00, &
     0.11248894313318662575D+00, &
     0.16823525155220746498D+00, &
     0.22338668642896688163D+00, &
     0.27774982202182431507D+00, &
     0.33113539325797683309D+00, &
     0.38335932419873034692D+00, &
     0.43424374934680255800D+00, &
     0.48361802694584102756D+00, &
     0.53131974364437562397D+00, &
     0.57719571005204581484D+00, &
     0.62110294673722640294D+00, &
     0.66290966002478059546D+00, &
     0.70249620649152707861D+00, &
     0.73975604435269475868D+00, &
     0.77459666924148337704D+00, &
     0.80694053195021761186D+00, &
     0.83672593816886873550D+00, &
     0.86390793819369047715D+00, &
     0.88845923287225699889D+00, &
     0.91037115695700429250D+00, &
     0.92965485742974005667D+00, &
     0.94634285837340290515D+00, &
     0.96049126870802028342D+00, &
     0.97218287474858179658D+00, &
     0.98153114955374010687D+00, &
     0.98868475754742947994D+00, &
     0.99383196321275502221D+00, &
     0.99720625937222195908D+00, &
     0.99909812496766759766D+00, &
     0.99987288812035761194D+00 /)
  real ( kind = rk ), dimension ( 63 ) :: w63 = (/ &
     0.000363221481845530659694D+00, &
     0.00126515655623006801137D+00, &
     0.00257904979468568827243D+00, &
     0.00421763044155885483908D+00, &
     0.00611550682211724633968D+00, &
     0.00822300795723592966926D+00, &
     0.0104982469096213218983D+00, &
     0.0129038001003512656260D+00, &
     0.0154067504665594978021D+00, &
     0.0179785515681282703329D+00, &
     0.0205942339159127111492D+00, &
     0.0232314466399102694433D+00, &
     0.0258696793272147469108D+00, &
     0.0284897547458335486125D+00, &
     0.0310735511116879648799D+00, &
     0.0336038771482077305417D+00, &
     0.0360644327807825726401D+00, &
     0.0384398102494555320386D+00, &
     0.0407155101169443189339D+00, &
     0.0428779600250077344929D+00, &
     0.0449145316536321974143D+00, &
     0.0468135549906280124026D+00, &
     0.0485643304066731987159D+00, &
     0.0501571393058995374137D+00, &
     0.0515832539520484587768D+00, &
     0.0528349467901165198621D+00, &
     0.0539054993352660639269D+00, &
     0.0547892105279628650322D+00, &
     0.0554814043565593639878D+00, &
     0.0559784365104763194076D+00, &
     0.0562776998312543012726D+00, &
     0.0563776283603847173877D+00, &
     0.0562776998312543012726D+00, &
     0.0559784365104763194076D+00, &
     0.0554814043565593639878D+00, &
     0.0547892105279628650322D+00, &
     0.0539054993352660639269D+00, &
     0.0528349467901165198621D+00, &
     0.0515832539520484587768D+00, &
     0.0501571393058995374137D+00, &
     0.0485643304066731987159D+00, &
     0.0468135549906280124026D+00, &
     0.0449145316536321974143D+00, &
     0.0428779600250077344929D+00, &
     0.0407155101169443189339D+00, &
     0.0384398102494555320386D+00, &
     0.0360644327807825726401D+00, &
     0.0336038771482077305417D+00, &
     0.0310735511116879648799D+00, &
     0.0284897547458335486125D+00, &
     0.0258696793272147469108D+00, &
     0.0232314466399102694433D+00, &
     0.0205942339159127111492D+00, &
     0.0179785515681282703329D+00, &
     0.0154067504665594978021D+00, &
     0.0129038001003512656260D+00, &
     0.0104982469096213218983D+00, &
     0.00822300795723592966926D+00, &
     0.00611550682211724633968D+00, &
     0.00421763044155885483908D+00, &
     0.00257904979468568827243D+00, &
     0.00126515655623006801137D+00, &
     0.000363221481845530659694D+00 /)

  if ( n == 1 ) then
    call r8vec_copy ( n, x01, x )
    call r8vec_copy ( n, w01, w )
  else if ( n == 3 ) then
    call r8vec_copy ( n, x03, x )
    call r8vec_copy ( n, w03, w )
  else if ( n == 7 ) then
    call r8vec_copy ( n, x07, x )
    call r8vec_copy ( n, w07, w )
  else if ( n == 15 ) then
    call r8vec_copy ( n, x15, x )
    call r8vec_copy ( n, w15, w )
  else if ( n == 31 ) then
    call r8vec_copy ( n, x31, x )
    call r8vec_copy ( n, w31, w )
  else if ( n == 63 ) then
    call r8vec_copy ( n, x63, x )
    call r8vec_copy ( n, w63, w )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KPU - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of N.'
    stop 1
  end if
!
!  The rule as stored is for the interval [-1,+1].
!  Adjust it to the interval [0,1].
!
  a = -1.0D+00
  b = +1.0D+00
  c =  0.0D+00
  d = +1.0D+00
  call rule_adjust ( a, b, c, d, n, x, w )

  return
end
subroutine kpu_order ( l, n )

!*****************************************************************************80
!
!! KPU_ORDER computes the order of a KPU rule from the level.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 May 2012
!
!  Author:
!
!    John Burkardt.
!
!  Parameters:
!
!    Input, integer L, the level of the rule.  
!    1 <= L <= 25
!
!    Output, integer N, the order of the rule.
!
  implicit none

  integer l
  integer n

  if ( l < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KPU_ORDER - Fatal error!'
    write ( *, '(a)' ) '  1 <= L <= 25 required.'
    write ( *, '(a,i4)' ) '  Input L = ', l
    stop 1
  else if ( l == 1 ) then
    n = 1
  else if ( l <= 3 ) then
    n = 3
  else if ( l <= 6 ) then
    n = 7
  else if ( l <= 12 ) then
    n = 15
  else if ( l <= 24 ) then
    n = 31
  else if ( l <= 25 ) then
    n = 63
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KPU_ORDER - Fatal error!'
    write ( *, '(a)' ) '  1 <= L <= 25 required.'
    write ( *, '(a,i4)' ) '  Input L = ', l
    stop 1
  end if

  return
end
subroutine num_seq ( n, k, seq_num )

!*****************************************************************************80
!
!! NUM_SEQ returns the number of compositions of the integer N into K parts.
!
!  Discussion:
!
!    A composition of the integer N into K parts is an ordered sequence
!    of K nonnegative integers which sum to N.  The compositions (1,2,1)
!    and (1,1,2) are considered to be distinct.
!
!    The 28 compositions of 6 into three parts are:
!
!      6 0 0,  5 1 0,  5 0 1,  4 2 0,  4 1 1,  4 0 2,
!      3 3 0,  3 2 1,  3 1 2,  3 0 3,  2 4 0,  2 3 1,
!      2 2 2,  2 1 3,  2 0 4,  1 5 0,  1 4 1,  1 3 2,
!      1 2 3,  1 1 4,  1 0 5,  0 6 0,  0 5 1,  0 4 2,
!      0 3 3,  0 2 4,  0 1 5,  0 0 6.
!
!    The formula for the number of compositions of N into K parts is
!
!      Number = ( N + K - 1 )! / ( N! * ( K - 1 )! )
!
!    Describe the composition using N '1's and K-1 dividing lines '|'.
!    The number of distinct permutations of these symbols is the number
!    of compositions.  This is equal to the number of permutations of
!    N+K-1 things, with N identical of one kind and K-1 identical of another.
!
!    Thus, for the above example, we have:
!
!      Number = ( 6 + 3 - 1 )! / ( 6! * (3-1)! ) = 28
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Second Edition,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer N, the integer whose compositions are desired.
!
!    Input, integer K, the number of parts in the composition.
!
!    Output, integer SEQ_NUM, the number of compositions of N
!    into K parts.
!
  implicit none

  integer i4_choose
  integer k
  integer n
  integer seq_num

  seq_num = i4_choose ( n + k - 1, n )

  return
end
subroutine nwspgr ( rule, rule_order, dim, k, r_size, s_size, nodes, weights )

!*****************************************************************************80
!
!! NWSPGR generates nodes and weights for sparse grid integration.
!
!  Discussion:
!
!    Thanks to Sergio Pastorello for pointing out an apparent out-of-bounds
!    error that actually indicated a mistake in the dimensions supplied
!    for two calls to R8CVV_RGET, 19 April 2013.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 April 2013
!
!  Author:
!
!    Original MATLAB version by Florian Heiss, Viktor Winschel.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Florian Heiss, Viktor Winschel,
!    Likelihood approximation by numerical integration on sparse grids,
!    Journal of Econometrics,
!    Volume 144, 2008, pages 62-80.
!
!  Parameters:
! 
!    Input, external RULE ( n, x, w ), the name of a subroutine which
!    is given the order N and returns the points X and weights W of the
!    corresponding 1D quadrature rule.
!
!    Input, external RULE_ORDER ( l, n ), the name of a subroutine which
!    is given the level L and returns the order N of the corresponding 1D rule.
!
!    Input, integer DIM, the spatial dimension.
!
!    Input, integer K, the level of the sparse rule.
!
!    Input, integer R_SIZE, the "size" of the sparse rule.
!
!    Output, integer S_SIZE, the size of the sparse rule,
!    after duplicate points have been merged.
!
!    Output, real ( kind = rk ) NODES(DIM,R_SIZE), the nodes of the sparse rule.
!
!    Output, real ( kind = rk ) WEIGHTS(R_SIZE), the weights of the sparse rule.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim
  integer r_size

  integer bq
  integer i
  integer i4_choose
  integer i4_mop
  integer, allocatable :: is(:,:)
  integer j
  integer k
  integer level
  integer, allocatable :: lr(:)
  integer maxq
  integer minq
  integer n
  integer, allocatable :: n1d(:)
  integer n1d_total
  integer nc
  real ( kind = rk ) nodes(dim,r_size)
  integer np
  integer, allocatable :: nr(:)
  integer q
  integer r
  integer, allocatable :: roff(:)
  integer, allocatable :: rq(:)
  integer s_size
  integer seq_num
  real ( kind = rk ), allocatable :: w(:)
  real ( kind = rk ), allocatable :: w1d(:)
  integer, allocatable :: w1d_off(:)
  real ( kind = rk ), allocatable :: wc(:)
  real ( kind = rk ) weights(r_size)
  real ( kind = rk ), allocatable :: wp(:)
  real ( kind = rk ), allocatable :: wr(:)
  real ( kind = rk ), allocatable :: x(:)
  real ( kind = rk ), allocatable :: x1d(:)
  integer, allocatable :: x1d_off(:)
  real ( kind = rk ), allocatable :: xc(:)
  real ( kind = rk ), allocatable :: xp(:,:)
  real ( kind = rk ), allocatable :: xr(:)

  nodes(1:dim,1:r_size) = 0.0D+00
  weights(1:r_size) = 0.0D+00
!
!  Create cell arrays that will contain the points and weights 
!  for levels 1 through K.
!
  allocate ( n1d(k) )
  allocate ( x1d_off(1:k+1) )
  allocate ( w1d_off(1:k+1) )

  x1d_off(1) = 0
  w1d_off(1) = 0

  do level = 1, k

    call rule_order ( level, n )
    n1d(level) = n
    x1d_off(level+1) = x1d_off(level) + n
    w1d_off(level+1) = w1d_off(level) + n

  end do

  n1d_total = x1d_off(k+1)
!
!  Calculate all the 1D rules needed.
!
  allocate ( x1d(n1d_total) )
  allocate ( w1d(n1d_total) )

  do level = 1, k

    n = n1d(level)

    allocate ( x(1:n) )
    allocate ( w(1:n) )

    call rule ( n, x, w )
    call r8cvv_rset ( n1d_total, x1d, k, x1d_off, level, x )
    call r8cvv_rset ( n1d_total, w1d, k, w1d_off, level, w )

    deallocate ( x )
    deallocate ( w )

  end do
!
!  Construct the sparse grid.
!
  minq = max ( 0, k - dim )
  maxq = k - 1
!
!  Q is the level total.
!
  allocate ( lr(1:dim) )
  allocate ( nr(1:dim) )
  allocate ( roff(1:dim+1) )

  r = 0

  do q = minq, maxq
!
!  BQ is the combinatorial coefficient applied to the component
!  product rules which have level Q.
!
    bq = i4_mop ( maxq - q ) * i4_choose ( dim - 1, dim + q - k )
!
!  Compute the D-dimensional row vectors that sum to DIM+Q.
!
    call num_seq ( q, dim, seq_num )

    allocate ( is(1:seq_num,1:dim) )

    call get_seq ( dim, q + dim, seq_num, is )
!
!  Allocate new rows for nodes and weights.
!
    allocate ( rq(1:seq_num) )

    do j = 1, seq_num
      rq(j) = product ( n1d(is(j,1:dim)) )
    end do
!
!  Generate every product rule whose level total is Q.
!
    do j = 1, seq_num

      lr(1:dim) = is(j,1:dim)

      do i = 1, dim
        call rule_order ( lr(i), nr(i) )
      end do

      call r8cvv_offset ( dim, nr, roff )

      nc = sum ( nr(1:dim) )
      allocate ( wc(nc) )
      allocate ( xc(nc) )

      do i = 1, dim
        allocate ( wr(1:nr(i)) )
        allocate ( xr(1:nr(i)) )
!
!  Incorrect dimension "NR(I)" replaced by correct dimension
!  N1D_TOTAL, 19 April 2013, JVB.
!
        call r8cvv_rget ( n1d_total, x1d, k, x1d_off, lr(i), xr )
        call r8cvv_rget ( n1d_total, w1d, k, w1d_off, lr(i), wr )
        call r8cvv_rset ( nc, xc, dim, roff, i, xr )
        call r8cvv_rset ( nc, wc, dim, roff, i, wr )
        deallocate ( wr )
        deallocate ( xr )
      end do

      np = rq(j)
      allocate ( wp(1:np) )
      allocate ( xp(1:dim,1:np) )

      call tensor_product_cell ( nc, xc, wc, dim, nr, roff, np, xp, wp )
!
!  Append the new nodes and weights to the arrays.
!
      nodes(1:dim,r+1:r+np) = xp(1:dim,1:np)
      weights(r+1:r+np) = bq * wp(1:np)
!
!  Update the count.
!
      r = r + rq(j)

      deallocate ( wc )
      deallocate ( wp )
      deallocate ( xc )
      deallocate ( xp )

    end do

    deallocate ( is )
    deallocate ( rq )

  end do

  deallocate ( lr )
  deallocate ( n1d )
  deallocate ( nr )
  deallocate ( roff )
  deallocate ( w1d )
  deallocate ( w1d_off )
  deallocate ( x1d )
  deallocate ( x1d_off )
!
!  Reorder the rule so the points are in ascending lexicographic order.
!
  call rule_sort ( dim, r_size, nodes, weights )
!
!  Suppress duplicate points and merge weights.
!
  r = 1
  do j = 2, r_size
    if ( all ( nodes(1:dim,r) == nodes(1:dim,j) ) ) then
      weights(r) = weights(r) + weights(j)
    else
      r = r + 1
      weights(r) = weights(j)
      nodes(1:dim,r) = nodes(1:dim,j)
    end if
  end do

  s_size = r
!
!  Zero out unneeded entries.
!
  nodes(1:dim,s_size+1:r_size) = 0.0D+00
  weights(s_size+1:r_size) = 0.0D+00
!  
!  Normalize the weights to sum to 1.
!
  weights(1:r) = weights(1:r) / sum ( weights(1:r) )

  return
end
subroutine nwspgr_size ( rule_order, dim, k, r_size )

!*****************************************************************************80
!
!! NWSPGR_SIZE determines the size of a sparse grid rule.
!
!  Discussion:
!
!    This routine does a "raw" count, that is, it does not notice that many
!    points may be repeated, in which case, the size of the rule could be
!    reduced by merging repeated points and combining the corresponding weights.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 December 2012
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Florian Heiss, Viktor Winschel,
!    Likelihood approximation by numerical integration on sparse grids,
!    Journal of Econometrics,
!    Volume 144, 2008, pages 62-80.
!
!  Parameters:
!
!    Input, external RULE_ORDER ( l, n ), the name of a subroutine which
!    is given the level L and returns the order N of the corresponding rule.
!
!    Input, integer DIM, the dimension of the integration problem.
!
!    Input, integer K, the level.  When using the built in 1D 
!    rules, the resulting sparse grid will be exact for polynomials up to total
!    order 2*K-1.  When using the built-in quadrature rules, the maximum value 
!    of K that is available is 25.
!
!    Output, integer R_SIZE, the "size" of the rule, that is,
!    the number of weights and multidimensional quadrature points that will
!    be needed.  The size of the rule will be reduced when duplicate points
!    are merged.
!
  implicit none

  integer dim

  integer i
  integer, allocatable :: is(:,:)
  integer k
  integer level
  integer maxq
  integer minq
  integer n
  integer, allocatable :: n1d(:)
  integer n1d_total
  integer q
  integer r_size
  integer, allocatable :: rq(:)
  external rule_order
  integer seq_num
!
!  Determine the size of each 1D rule.
!
  allocate ( n1d(k) )

  do level = 1, k

    call rule_order ( level, n )
    n1d(level) = n

  end do

  n1d_total = sum ( n1d(1:k) )
!
!  Go through the motions of generating the rules.
!
  minq = max ( 0, k - dim )
  maxq = k - 1
  r_size = 0

  do q = minq, maxq
!
!  Compute the D-dimensional vectors that sum to Q+DIM.
!
    call num_seq ( q, dim, seq_num )

    allocate ( is(1:seq_num,1:dim) )

    call get_seq ( dim, q + dim, seq_num, is )
!
!  Determine the size of each rule.
!
    allocate ( rq(1:seq_num) )

    do i = 1, seq_num
      rq(i) = product ( n1d(is(i,1:dim)) )
    end do
!
!  Add the sizes to the total.
!
    r_size = r_size + sum ( rq )

    deallocate ( is )
    deallocate ( rq )

  end do

  deallocate ( n1d )

  return
end
subroutine quad_rule_print ( m, n, x, w, title )

!*****************************************************************************80
!
!! QUAD_RULE_PRINT prints a multidimensional quadrature rule.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 May 2012
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
!
!    Input, real ( kind = rk ) X(M,N), the abscissas.
!
!    Input, real ( kind = rk ) W(N), the weights.
!
!    Input, character ( len = * ) TITLE, a title for the rule.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer i
  integer j
  character ( len = * ) title
  real ( kind = rk ) w(n)
  real ( kind = rk ) x(m,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do j = 1, n
    write ( *, '(2x,i2,2x,f10.6,a)', advance = 'no' ) j, w(j), ' * f ('
    do i = 1, m
      write ( *, '(f10.6)', advance = 'no' ) x(i,j)
      if ( i < m ) then
        write ( *, '('','')', advance = 'no' )
      else
        write ( *, '('' )'')', advance = 'yes' )
      end if
    end do
  end do

  return
end
subroutine r8cvv_offset ( m, nr, roff )

!*****************************************************************************80
!
!! R8CVV_OFFSET determines the row offsets of an R8CVV.
!
!  Discussion:
!
!    An R8CVV is a "vector of vectors" of R8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows in the array.
!
!    Input, integer NR(M), the row sizes.
!
!    Output, integer ROFF(M+1), the row offsets.
!
  implicit none

  integer m

  integer i
  integer roff(m+1)
  integer nr(m)

  roff(1) = 0
  do i = 1, m
    roff(i+1) = roff(i) + nr(i)
  end do

  return
end
subroutine r8cvv_print ( mn, a, m, roff, title )

!*****************************************************************************80
!
!! R8CVV_PRINT prints an R8CVV.
!
!  Discussion:
!
!    An R8CVV is a "vector of vectors" of R8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MN, the size of the cell array.
!
!    Input, real ( kind = rk ) A(MN), the cell array.
!
!    Input, integer M, the number of rows in the array.
!
!    Input, integer ROFF(M+1), the row offsets.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer mn

  real ( kind = rk ) a(mn)
  integer i
  integer k1
  integer k2
  integer khi
  integer klo
  integer roff(m+1)
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, m

    k1 = roff(i) + 1
    k2 = roff(i+1)

    do klo = k1, k2, 5
      khi = min ( klo + 5 - 1, k2 )
      if ( klo == k1 ) then
        write ( *, '(i5,2x, 5g14.6)' ) i, a(klo:khi)
      else
        write ( *, '(5x,2x, 5g14.6)' )    a(klo:khi)
      end if
    end do

  end do

  return
end
subroutine r8cvv_rget ( mn, a, m, roff, i, ai )

!*****************************************************************************80
!
!! R8CVV_RGET gets row I from an R8CVV.
!
!  Discussion:
!
!    An R8CVV is a "vector of vectors" of R8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MN, the size of the cell array.
!
!    Input, real ( kind = rk ) A(MN), the cell array.
!
!    Input, integer M, the number of rows in the array.
!
!    Input, integer ROFF(M+1), the row offsets.
!
!    Input, integer I, the row.
!    1 <= I <= M.
!
!    Output, real ( kind = rk ) AI(NR(I)), the value of A(I,*).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer mn

  real ( kind = rk ) a(mn)
  real ( kind = rk ) ai(*)
  integer i
  integer k1
  integer k2
  integer nv
  integer roff(m+1)

  k1 = roff(i) + 1
  k2 = roff(i+1)
  nv = k2 + 1 - k1

  ai(1:nv) = a(k1:k2)

  return
end
subroutine r8cvv_rset ( mn, a, m, roff, i, ai )

!*****************************************************************************80
!
!! R8CVV_RSET sets row I from an R8CVV.
!
!  Discussion:
!
!    An R8CVV is a "vector of vectors" of R8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MN, the size of the cell array.
!
!    Input/output, real ( kind = rk ) A(MN), the cell array.
!
!    Input, integer M, the number of rows in the array.
!
!    Input, integer ROFF(M+1), the row offsets.
!
!    Input, integer I, the row.
!    1 <= I <= M.
!
!    Input, real ( kind = rk ) AI(NR(I)), the new value of A(I,*).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer mn

  real ( kind = rk ) a(mn)
  real ( kind = rk ) ai(*)
  integer i
  integer k1
  integer k2
  integer nv
  integer roff(m+1)

  k1 = roff(i) + 1
  k2 = roff(i+1)
  nv = k2 + 1 - k1
  a(k1:k2) = ai(1:nv)

  return
end
subroutine r8mat_normal_01 ( m, n, seed, r )

!*****************************************************************************80
!
!! R8MAT_NORMAL_01 returns a unit pseudonormal R8MAT.
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
!    Input/output, integer SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = rk ) R(M,N), the array of pseudonormal values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer seed
  real ( kind = rk ) r(m,n)

  call r8vec_normal_01 ( m * n, seed, r )

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
!    Output, real ( kind = rk ) R(M,N), the array of pseudorandom values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer i
  integer, parameter :: i4_huge = 2147483647
  integer j
  integer k
  integer seed
  real ( kind = rk ) r(m,n)

  do j = 1, n

    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge
      end if

      r(i,j) = real ( seed, kind = rk ) * 4.656612875D-10

    end do
  end do

  return
end
subroutine r8vec_copy ( n, a1, a2 )

!*****************************************************************************80
!
!! R8VEC_COPY copies an R8VEC.
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
!    17 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the length of the vectors.
!
!    Input, real ( kind = rk ) A1(N), the vector to be copied.
!
!    Output, real ( kind = rk ) A2(N), a copy of A1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a1(n)
  real ( kind = rk ) a2(n)

  a2(1:n) = a1(1:n)

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
!    Local, integer SKIP, the distance from the current value of 
!    START to the next location of a block of values to set.
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
subroutine r8vec_direct_product2 ( factor_index, factor_order, factor_value, &
  factor_num, point_num, w )

!*****************************************************************************80
!
!! R8VEC_DIRECT_PRODUCT2 creates a direct product of R8VEC's.
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
!    This routine carries out the task involving the weights W.
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
!      W(1:4) = ( 2, 3, 5, 7 )
!
!    Rule 2:
!      Order = 3
!      W(1:3) = ( 11, 13, 17 )
!
!    Rule 3:
!      Order = 2
!      W(1:2) = ( 19, 23 )
!
!    Product Rule:
!      Order = 24
!      W(1:24) =
!        ( 2 * 11 * 19 )
!        ( 3 * 11 * 19 )
!        ( 4 * 11 * 19 )
!        ( 7 * 11 * 19 )
!        ( 2 * 13 * 19 )
!        ( 3 * 13 * 19 )
!        ( 5 * 13 * 19 )
!        ( 7 * 13 * 19 )
!        ( 2 * 17 * 19 )
!        ( 3 * 17 * 19 )
!        ( 5 * 17 * 19 )
!        ( 7 * 17 * 19 )
!        ( 2 * 11 * 23 )
!        ( 3 * 11 * 23 )
!        ( 5 * 11 * 23 )
!        ( 7 * 11 * 23 )
!        ( 2 * 13 * 23 )
!        ( 3 * 13 * 23 )
!        ( 5 * 13 * 23 )
!        ( 7 * 13 * 23 )
!        ( 2 * 17 * 23 )
!        ( 3 * 17 * 23 )
!        ( 5 * 17 * 23 )
!        ( 7 * 17 * 23 )
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
!    Input/output, real ( kind = rk ) W(POINT_NUM), the elements of the
!    direct product, which are built up gradually.
!
!  Local Parameters:
!
!    Local, integer START, the first location of a block of values
!    to set.
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
  real ( kind = rk ) w(point_num)

  if ( factor_index == 1 ) then
    contig = 1
    skip = 1
    rep = point_num
    w(1:point_num) = 1.0D+00
  end if

  rep = rep / factor_order
  skip = skip * factor_order

  do j = 1, factor_order

    start = 1 + ( j - 1 ) * contig

    do k = 1, rep
      w(start:start+contig-1) = w(start:start+contig-1) * factor_value(j)
      start = start + skip
    end do

  end do

  contig = contig * factor_order

  return
end
subroutine r8vec_normal_01 ( n, seed, x )

!*****************************************************************************80
!
!! R8VEC_NORMAL_01 returns a unit pseudonormal R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!    This routine can generate a vector of values on one call.  It
!    has the feature that it should provide the same results
!    in the same order no matter how we break up the task.
!
!    Before calling this routine, the user may call RANDOM_SEED
!    in order to set the seed of the random number generator.
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
!  Parameters:
!
!    Input, integer N, the number of values desired.  If N is
!    negative,then the code will flush its internal memory; in particular,
!    if there is a saved value to be used on the next call, it is
!    instead discarded.  This is useful if the user has reset the
!    random number seed, for instance.
!
!    Input/output, integer SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = rk ) X(N), a sample of the standard normal PDF.
!
!  Local parameters:
!
!    Local, integer MADE, records the number of values that have
!    been computed.  On input with negative N, this value overwrites
!    the return value of N, so the user can get an accounting of
!    how much work has been done.
!
!    Local, real ( kind = rk ) R(N+1), is used to store some uniform
!    random values.  Its dimension is N+1, but really it is only needed
!    to be the smallest even number greater than or equal to N.
!
!    Local, integer SAVED, is 0 or 1 depending on whether there is a
!    single saved value left over from the previous call.
!
!    Local, integer X_LO_INDEX, X_HI_INDEX, records the range of entries of
!    X that we need to compute.  This starts off as 1:N, but is adjusted
!    if we have a saved value that can be immediately stored in X(1),
!    and so on.
!
!    Local, real ( kind = rk ) Y, the value saved from the previous call, if
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
  integer seed
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
             sqrt ( -2.0D+00 * log ( r(1) ) ) * cos ( 2.0D+00 * pi * r(2) )
    y =      sqrt ( -2.0D+00 * log ( r(1) ) ) * sin ( 2.0D+00 * pi * r(2) )

    saved = 1

    made = made + 2
!
!  If we require an even number of values, that's easy.
!
  else if ( mod ( x_hi_index - x_lo_index + 1, 2 ) == 0 ) then

    m = ( x_hi_index - x_lo_index + 1 ) / 2

    call r8vec_uniform_01 ( 2*m, seed, r )

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

    call r8vec_uniform_01 ( 2*m, seed, r )

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
!    My vector:
!
!        1.0    2.1    3.2    4.3    5.4
!        6.5    7.6    8.7    9.8   10.9
!       11.0
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
  integer ihi
  integer ilo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do ilo = 1, n, 5
    ihi = min ( ilo + 5 - 1, n )
    write ( *, '(5g14.6)' ) a(ilo:ihi)
  end do

  return
end
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    For now, the input quantity SEED is an integer variable.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 July 2006
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
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input/output, integer SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = rk ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer k
  integer seed
  real ( kind = rk ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + 2147483647
    end if

    r(i) = real ( seed, kind = rk ) * 4.656612875D-10

  end do

  return
end
subroutine rule_adjust ( a, b, c, d, n, x, w )

!*****************************************************************************80
!
!! RULE_ADJUST adjusts a 1D quadrature rule from [A,B] to [C,D].
!
!  Discussion:
!
!    This function is only appropriate for cases involving finite intervals
!    and a uniform weighting function.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, B, the left and right endpoints of the
!    original interval.
!
!    Input, real ( kind = rk ) C, D, the left and right endpoints of the
!    new interval.
!
!    Input, integer N, the order of the rule.
!    1 <= N.
!
!    Input/output, real ( kind = rk ) X(N), the abscissas.
!
!    Input/output, real ( kind = rk ) W(N), the weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) d
  integer i
  real ( kind = rk ) s
  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n)

  do i = 1, n
    x(i) = ( ( b - x(i)     ) * c   &
           + (     x(i) - a ) * d ) &
           / ( b        - a )
  end do

  s = ( d - c ) / ( b - a )

  w(1:n) = s * w(1:n)

  return
end
subroutine rule_sort ( m, n, x, w )

!*****************************************************************************80
!
!! RULE_SORT sorts a multidimensional quadrature rule.
!
!  Discussion:
!
!    This routine simply reindexes the items in the rule so that the points
!    occur in increasing lexicographic order.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input/output, real ( kind = rk ) X(M,N), W(N).
!    The points and weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer i
  integer indx
  integer isgn
  integer j1
  integer j2
  real ( kind = rk ) t(m)
  real ( kind = rk ) w(n)
  real ( kind = rk ) ww
  real ( kind = rk ) x(m,n)

  if ( m <= 0 ) then
    return
  end if

  if ( n <= 1 ) then
    return
  end if
!
!  Initialize.
!
  indx = 0
  isgn = 0
  j1 = 0
  j2 = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( n, indx, j1, j2, isgn )
!
!  Interchange columns J1 and J2.
!
    if ( 0 < indx ) then

      t(1:m)    = x(1:m,j1)
      x(1:m,j1) = x(1:m,j2)
      x(1:m,j2) = t(1:m)

      ww = w(j1)
      w(j1) = w(j2)
      w(j2) = ww      
!
!  Compare columns J1 and J2.
!
    else if ( indx < 0 ) then

      isgn = 0

      do i = 1, m 
 
        if ( x(i,j1) < x(i,j2) ) then
          isgn = -1
          exit
        else if ( x(i,j2) < x(i,j1) ) then
          isgn = +1
          exit
        end if

      end do
!
!  The columns are sorted.
!
    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine sort_heap_external ( n, indx, i, j, isgn )

!*****************************************************************************80
!
!! SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
!
!  Discussion:
!
!    The actual list of data is not passed to the routine.  Hence this
!    routine may be used to sort integers, reals, numbers, names,
!    dates, shoe sizes, and so on.  After each call, the routine asks
!    the user to compare or interchange two items, until a special
!    return value signals that the sorting is completed.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 February 2004
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer N, the number of items to be sorted.
!
!    Input/output, integer INDX, the main communication signal.
!    The user must set INDX to 0 before the first call.
!    Thereafter, the user should not change the value of INDX until
!    the sorting is done.
!    On return, if INDX is
!    * greater than 0,
!      > interchange items I and J;
!      > call again.
!    * less than 0,
!      > compare items I and J;
!      > set ISGN = -1 if I < J, ISGN = +1 if J < I;
!      > call again.
!    * equal to 0, the sorting is done.
!
!    Output, integer I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ISGN, results of comparison of elements
!    I and J. (Used only when the previous call returned INDX less than 0).
!    ISGN <= 0 means I is less than or equal to J;
!    0 <= ISGN means I is greater than or equal to J.
!
  implicit none

  integer i
  integer, save :: i_save = 0
  integer indx
  integer isgn
  integer j
  integer, save :: j_save = 0
  integer, save :: k = 0
  integer, save :: k1 = 0
  integer n
  integer, save :: n1 = 0
!
!  INDX = 0: This is the first call.
!
  if ( indx == 0 ) then

    i_save = 0
    j_save = 0
    k = n / 2
    k1 = k
    n1 = n
!
!  INDX < 0: The user is returning the results of a comparison.
!
  else if ( indx < 0 ) then

    if ( indx == -2 ) then

      if ( isgn < 0 ) then
        i_save = i_save + 1
      end if

      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return

    end if

    if ( 0 < isgn ) then
      indx = 2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then

      if ( n1 == 1 ) then
        i_save = 0
        j_save = 0
        indx = 0
      else
        i_save = n1
        n1 = n1 - 1
        j_save = 1
        indx = 1
      end if

      i = i_save
      j = j_save
      return

    end if

    k = k - 1
    k1 = k
!
!  0 < INDX, the user was asked to make an interchange.
!
  else if ( indx == 1 ) then

    k1 = k

  end if

  do

    i_save = 2 * k1

    if ( i_save == n1 ) then
      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return
    else if ( i_save <= n1 ) then
      j_save = i_save + 1
      indx = -2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then
      exit
    end if

    k = k - 1
    k1 = k

  end do

  if ( n1 == 1 ) then
    i_save = 0
    j_save = 0
    indx = 0
    i = i_save
    j = j_save
  else
    i_save = n1
    n1 = n1 - 1
    j_save = 1
    indx = 1
    i = i_save
    j = j_save
  end if

  return
end
subroutine symmetric_sparse_size ( nr, dim, nodes, x0, nr2 )

!*****************************************************************************80
!
!! SYMMETRIC_SPARSE_SIZE sizes a symmetric sparse rule.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 May 2012
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Florian Heiss, Viktor Winschel,
!    Likelihood approximation by numerical integration on sparse grids,
!    Journal of Econometrics,
!    Volume 144, 2008, pages 62-80.
!
!  Parameters:
! 
!    Input, integer DIM, the dimension.
!
!    Input, integer NR, the dimension of the rule in the 
!    positive orthant.
!
!    Input, real ( kind = rk ) NODES(NR,DIM), the nodes for the positive orthant.
!
!    Input, real ( kind = rk ) X0, the point of symmetry for the 1D rule, 
!    typically 0.
!
!    Output, integer NR2, the dimension of the rule when 
!    "unfolded" to the full space.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim
  integer nr

  integer count
  integer j
  real ( kind = rk ) nodes(nr,dim)
  integer nr2
  integer r
  real ( kind = rk ) x0
!
!  Count the size of the full rule.
!
  nr2 = 0

  do r = 1, nr
    count = 1
    do j = 1, dim
      if ( nodes(r,j) /= x0 ) then
        count = 2 * count
      end if
    end do
    nr2 = nr2 + count
  end do

  return
end
subroutine tensor_product ( d, order1d, n1d, x1d, w1d, n, xnd, wnd )

!*****************************************************************************80
!
!! TENSOR_PRODUCT generates a tensor product quadrature rule.
!
!  Discussion:
!
!    The Kronecker product of an K by L matrix A and an M by N matrix B
!    is the K*M by L*N matrix formed by
!
!      a(1,1) * B,  a(1,2) * B,  ..., a(1,l) * B
!      a(2,1) * B,  a(2,2) * B,  ..., a(2,l) * B
!      ..........   ..........   .... ..........
!      a(k,1) * B,  a(k,2) * B,  ..., a(k,l) * B
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 April 2012
!
!  Author:
!
!    Original MATLAB version by Florian Heiss, Viktor Winschel.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Florian Heiss, Viktor Winschel,
!    Likelihood approximation by numerical integration on sparse grids,
!    Journal of Econometrics,
!    Volume 144, 2008, pages 62-80.
!
!  Parameters:
!
!    Input, integer D, the spatial dimension.
!
!    Input, integer ORDER1D(D), the order of each 1D rule.
!
!    Input, integer N1D, the number of 1D items.
!
!    Input, real ( kind = rk ) X1D(N1D), the 1D nodes.
!
!    Input, real ( kind = rk ) W1D(N1D), the 1D weights.
!
!    Input, integer N, the number of N-dimensional items.
!
!    Output, real ( kind = rk ) XND(D,N), the nodes.
!
!    Output, real ( kind = rk ) WND(N), the weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer d
  integer n
  integer n1d

  integer i
  integer i1
  integer i2
  integer order1d(d)
  real ( kind = rk ) w1d(n1d)
  real ( kind = rk ) wnd(n)
  real ( kind = rk ) x1d(n1d)
  real ( kind = rk ) xnd(d,n)
!
!  Compute the weights.
!
  i2 = 0
  do i = 1, d
    i1 = i2 + 1
    i2 = i2 + order1d(i)
    call r8vec_direct_product2 ( i, order1d(i), w1d(i1:i2), d, n, wnd )
  end do
!
!  Compute the points.
!
  i2 = 0
  do i = 1, d
    i1 = i2 + 1
    i2 = i2 + order1d(i)
    call r8vec_direct_product ( i, order1d(i), x1d(i1:i2), d, n, xnd )
  end do

  return
end
subroutine tensor_product_cell ( nc, xc, wc, dim, nr, roff, np, xp, wp )

!*****************************************************************************80
!
!! TENSOR_PRODUCT_CELL generates a tensor product quadrature rule.
!
!  Discussion:
!
!    The Kronecker product of an K by L matrix A and an M by N matrix B
!    is the K*M by L*N matrix formed by
!
!      a(1,1) * B,  a(1,2) * B,  ..., a(1,l) * B
!      a(2,1) * B,  a(2,2) * B,  ..., a(2,l) * B
!      ..........   ..........   .... ..........
!      a(k,1) * B,  a(k,2) * B,  ..., a(k,l) * B
!
!    The 1D factors are stored in a kind of cell array structure.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 December 2012
!
!  Author:
!
!    John Burkardt.
!
!  Parameters:
!
!    Input, integer NC, the number of items in the cell arrays.
!
!    Input, real ( kind = rk ) XC(NC), a cell array containing points for 
!    1D rules.
!
!    Input, real ( kind = rk ) WC(NC), a cell array containing weights for
!    1D rules.
!
!    Input, integer DIM, the spatial dimension.
!
!    Input, integer NR(DIM), the length of each row of the 
!    cell array.
!
!    Input, integer ROFF(DIM+1), offsets for the cell arrays.
!
!    Input, integer NP, the number of points in the product rule.
!
!    Output, real ( kind = rk ) XP(DIM,NP), the nodes.
!
!    Output, real ( kind = rk ) WP(NP), the weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dim
  integer nc
  integer np

  integer i
  integer n1d
  integer nr(dim)
  integer roff(dim+1)
  real ( kind = rk ), allocatable :: w1d(:)
  real ( kind = rk ) wc(nc)
  real ( kind = rk ) wp(np)
  real ( kind = rk ), allocatable :: x1d(:)
  real ( kind = rk ) xc(nc)
  real ( kind = rk ) xp(dim,np)
!
!  Compute the weights.
!
  do i = 1, dim
    n1d = nr(i)
    allocate ( w1d(1:n1d) )
    call r8cvv_rget ( nc, wc, dim, roff, i, w1d )
    call r8vec_direct_product2 ( i, n1d, w1d, dim, np, wp )
    deallocate ( w1d )
  end do
!
!  Compute the points.
!
  do i = 1, dim
    n1d = nr(i)
    allocate ( x1d(1:n1d) )
    call r8cvv_rget ( nc, xc, dim, roff, i, x1d )
    call r8vec_direct_product ( i, n1d, x1d, dim, np, xp )
    deallocate ( x1d )
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
