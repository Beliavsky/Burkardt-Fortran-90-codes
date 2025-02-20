program main

!*****************************************************************************80
!
!! truncated_normal_rule() computes a truncated normal quadrature rule.
!
!  Discussion:
!
!    The user specifies:
!    * option: 0/1/2/3 for none, lower, upper, double truncation.
!    * N, the number of points in the rule;
!    * MU, the mean of the original normal distribution;
!    * SIGMA, the standard deviation of the original normal distribution,
!    * A, the left endpoint (for options 1 or 3)
!    * B, the right endpoint (for options 2 or 3);
!    * FILENAME, the root name of the output files.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 September 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  integer arg_num
  real ( kind = rk ) b
  character ( len = 255 ) filename
  integer iarg
  integer iargc
  integer ierror
  integer last
  real ( kind = rk ), allocatable, dimension ( : ) :: moment
  real ( kind = rk ) mu
  integer n
  integer option
  real ( kind = rk ) r(2)
  real ( kind = rk ) r8_huge
  real ( kind = rk ) sigma
  character ( len = 255 ) string
  real ( kind = rk ), allocatable, dimension ( : ) :: w
  real ( kind = rk ), allocatable, dimension ( : ) :: x

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'truncated_normal_rule():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For the (truncated) Gaussian probability density function'
  write ( *, '(a)' ) '    pdf(x) = exp(-0.5*((x-MU)/SIGMA)^2) / SIGMA / sqrt ( 2 * pi )'
  write ( *, '(a)' ) '  compute an N-point quadrature rule for approximating'
  write ( *, '(a)' ) '    Integral ( A <= x <= B ) f(x) pdf(x) dx'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The value of OPTION determines the truncation interval [A,B]:'
  write ( *, '(a)' ) '  0: (-oo,+oo)'
  write ( *, '(a)' ) '  1: [A,+oo)'
  write ( *, '(a)' ) '  2: (-oo,B]'
  write ( *, '(a)' ) '  3: [A,B]'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The user specifies OPTION, N, MU, SIGMA, A, B and FILENAME.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FILENAME is used to generate 3 files:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    filename_w.txt - the weight file'
  write ( *, '(a)' ) '    filename_x.txt - the abscissa file.'
  write ( *, '(a)' ) '    filename_r.txt - the region file, listing A and B.'
!
!  Get the number of command line arguments.
!
  arg_num = iargc ( )
  iarg = 0
!
!  Get OPTION.
!
  iarg = iarg + 1
  if ( iarg <= arg_num ) then
    call getarg ( iarg, string )
    call s_to_i4 ( string, option, ierror, last )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter OPTION, 0/1/2/3:'
    read ( *, * ) option
  end if

  if ( option < 0 .or. 3 < option ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRUNCATED_NORMAL_RULE - Fatal error!'
    write ( *, '(a)' ) '  0 <= OPTION <= 3 was required.'
    stop 1
  end if
!
!  Get N.
!
  iarg = iarg + 1
  if ( iarg <= arg_num ) then
    call getarg ( iarg, string )
    call s_to_i4 ( string, n, ierror, last )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter N, the number of quadrature points:'
    read ( *, * ) n
  end if
!
!  Get MU.
!
  iarg = iarg + 1
  if ( iarg <= arg_num ) then
    call getarg ( iarg, string )
    call s_to_r8 ( string, mu, ierror, last )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter MU, the mean value of the normal distribution:'
    read ( *, * ) mu
  end if
!
!  Get SIGMA.
!
  iarg = iarg + 1
  if ( iarg <= arg_num ) then
    call getarg ( iarg, string )
    call s_to_r8 ( string, sigma, ierror, last )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter SIGMA, the standard deviation of the normal distribution:'
    read ( *, * ) sigma
  end if

  sigma = abs ( sigma )
!
!  Get A, perhaps.
!
  if ( option == 1 .or. option == 3 ) then
    iarg = iarg + 1
    if ( iarg <= arg_num ) then
      call getarg ( iarg, string )
      call s_to_r8 ( string, a, ierror, last )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Enter A, the left endpoint:'
      read ( *, * ) a
    end if
  else
    a = - r8_huge ( )
  end if
!
!  Get B.
!
  if ( option == 2 .or. option == 3 ) then
    iarg = iarg + 1
    if ( iarg <= arg_num ) then
      call getarg ( iarg, string )
      call s_to_r8 ( string, b, ierror, last )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Enter B, the right endpoint:'
      read ( *, * ) b
    end if
  else
    b = r8_huge ( )
  end if

  if ( b <= a ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRUNCATED_NORMAL_RULE - Fatal error!'
    write ( *, '(a)' ) '  A < B required.'
    stop 1
  end if
!
!  Get FILENAME.
!
  iarg = iarg + 1
  if ( iarg <= arg_num ) then
    call getarg ( iarg, filename )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) &
      '  Enter FILENAME, the "root name" of the quadrature files).'
    read ( *, '(a)' ) filename
  end if
!
!  Input summary.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' )     '  OPTION = ', option
  write ( *, '(a,i8)' )     '  N = ', n
  write ( *, '(a,g14.6)' ) '  MU = ', mu
  write ( *, '(a,g14.6)' ) '  SIGMA = ', sigma
  if ( option == 1 .or. option == 3 ) then
    write ( *, '(a,g14.6)' ) '  A = ', a
  else
    write ( *, '(a)' ) '  A = -oo'
  end if
  if ( option == 2 .or. option == 3 ) then
    write ( *, '(a,g14.6)' ) '  B = ', b
  else
    write ( *, '(a)' ) '  B = +oo'
  end if
  write ( *, '(a)' )        '  FILENAME = "' // trim ( filename ) // '".'
!
!  Compute the moments.
!
  allocate ( moment(0:2*n) )

  if ( option == 0 ) then
    call moments_normal ( 2 * n + 1, mu, sigma, moment )
  else if ( option == 1 ) then
    call moments_truncated_normal_a ( 2 * n + 1, mu, sigma, a, moment )
  else if ( option == 2 ) then
    call moments_truncated_normal_b ( 2 * n + 1, mu, sigma, b, moment )
  else if ( option == 3 ) then
    call moments_truncated_normal_ab ( 2 * n + 1, mu, sigma, a, b, moment )
  end if
!
!  Construct the rule from the moments.
!
  allocate ( w(n) )
  allocate ( x(n) )

  call moment_method ( n, moment, x, w )
!
!  Write the rule to a file.
!
  r(1) = a
  r(2) = b

  call rule_write ( n, x, w, r, filename )
!
!  Free memory.
!
  deallocate ( moment )
  deallocate ( w )
  deallocate ( x )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'truncated_normal_rule():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine ch_cap ( ch )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
!
!  Discussion:
!
!    Instead of CHAR and ICHAR, we now use the ACHAR and IACHAR functions, 
!    which guarantee the ASCII collating sequence.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character CH, the character to capitalize.
!
  implicit none

  character ch
  integer itemp

  itemp = iachar ( ch )
 
  if ( 97 <= itemp .and. itemp <= 122 ) then
    ch = achar ( itemp - 32 )
  end if
 
  return
end
function ch_eqi ( c1, c2 )

!*****************************************************************************80
!
!! CH_EQI is a case insensitive comparison of two characters for equality.  
!
!  Discussion:
!
!    CH_EQI ( 'A', 'a' ) is TRUE.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    28 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C1, C2, the characters to compare.
!
!    Output, logical CH_EQI, the result of the comparison.
!
  implicit none

  character c1
  character c1_cap
  character c2
  character c2_cap
  logical ch_eqi

  c1_cap = c1
  c2_cap = c2

  call ch_cap ( c1_cap )
  call ch_cap ( c2_cap )

  if ( c1_cap == c2_cap ) then
    ch_eqi = .true.
  else
    ch_eqi = .false.
  end if

  return
end
subroutine ch_to_digit ( ch, digit )

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the value of a base 10 digit.
!
!  Discussion:
!
!    Instead of ICHAR, we now use the IACHAR function, which
!    guarantees the ASCII collating sequence.
!
!  Example:
!
!     CH  DIGIT
!    ---  -----
!    '0'    0
!    '1'    1
!    ...  ...
!    '9'    9
!    ' '    0
!    'X'   -1
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character CH, the decimal digit, '0' through '9' or blank
!    are legal. 
!
!    Output, integer DIGIT, the corresponding value.  
!    If CH was 'illegal', then DIGIT is -1.
!
  implicit none

  character ch
  integer digit

  if ( lle ( '0', ch ) .and. lle ( ch, '9' ) ) then
 
    digit = iachar ( ch ) - 48
 
  else if ( ch == ' ' ) then
 
    digit = 0
 
  else

    digit = -1

  end if
 
  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 January 2008
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
subroutine jacobi_eigenvalue ( n, a, it_max, v, d, it_num, rot_num )

!*****************************************************************************80
!
!! JACOBI_EIGENVALUE carries out the Jacobi eigenvalue iteration.
!
!  Discussion:
!
!    This function computes the eigenvalues and eigenvectors of a
!    real symmetric matrix, using Rutishauser's modfications of the classical
!    Jacobi rotation method with threshold pivoting. 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 September 2013
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real ( kind = rk ) A(N,N), the matrix, which must be square, real,
!    and symmetric.
!
!    Input, integer IT_MAX, the maximum number of iterations.
!
!    Output, real ( kind = rk ) V(N,N), the matrix of eigenvectors.
!
!    Output, real ( kind = rk ) D(N), the eigenvalues, in descending order.
!
!    Output, integer IT_NUM, the total number of iterations.
!
!    Output, integer ROT_NUM, the total number of rotations.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) bw(n)
  real ( kind = rk ) c
  real ( kind = rk ) d(n)
  real ( kind = rk ) g
  real ( kind = rk ) gapq
  real ( kind = rk ) h
  integer i
  integer it_max
  integer it_num
  integer j
  integer k
  integer l
  integer m
  integer p
  integer q
  integer rot_num
  real ( kind = rk ) s
  real ( kind = rk ) t
  real ( kind = rk ) tau
  real ( kind = rk ) term
  real ( kind = rk ) termp
  real ( kind = rk ) termq
  real ( kind = rk ) theta
  real ( kind = rk ) thresh
  real ( kind = rk ) v(n,n)
  real ( kind = rk ) w(n)
  real ( kind = rk ) zw(n)

  do j = 1, n
    do i = 1, n
      v(i,j) = 0.0D+00
    end do
    v(j,j) = 1.0D+00
  end do

  do i = 1, n
    d(i) = a(i,i)
  end do

  bw(1:n) = d(1:n)
  zw(1:n) = 0.0D+00
  it_num = 0
  rot_num = 0

  do while ( it_num < it_max )

    it_num = it_num + 1
!
!  The convergence threshold is based on the size of the elements in
!  the strict upper triangle of the matrix.
!
    thresh = 0.0D+00
    do j = 1, n
      do i = 1, j - 1
        thresh = thresh + a(i,j) ** 2
      end do
    end do

    thresh = sqrt ( thresh ) / real ( 4 * n, kind = rk )

    if ( thresh == 0.0D+00 ) then
      exit 
    end if

    do p = 1, n
      do q = p + 1, n

        gapq = 10.0D+00 * abs ( a(p,q) )
        termp = gapq + abs ( d(p) )
        termq = gapq + abs ( d(q) )
!
!  Annihilate tiny offdiagonal elements.
!
        if ( 4 < it_num .and. &
             termp == abs ( d(p) ) .and. &
             termq == abs ( d(q) ) ) then

          a(p,q) = 0.0D+00
!
!  Otherwise, apply a rotation.
!
        else if ( thresh <= abs ( a(p,q) ) ) then

          h = d(q) - d(p)
          term = abs ( h ) + gapq

          if ( term == abs ( h ) ) then
            t = a(p,q) / h
          else
            theta = 0.5D+00 * h / a(p,q)
            t = 1.0D+00 / ( abs ( theta ) + sqrt ( 1.0D+00 + theta * theta ) )
            if ( theta < 0.0D+00 ) then 
              t = - t
            end if
          end if

          c = 1.0D+00 / sqrt ( 1.0D+00 + t * t )
          s = t * c
          tau = s / ( 1.0D+00 + c )
          h = t * a(p,q)
!
!  Accumulate corrections to diagonal elements.
!
          zw(p) = zw(p) - h                  
          zw(q) = zw(q) + h
          d(p) = d(p) - h
          d(q) = d(q) + h

          a(p,q) = 0.0D+00
!
!  Rotate, using information from the upper triangle of A only.
!
          do j = 1, p - 1
            g = a(j,p)
            h = a(j,q)
            a(j,p) = g - s * ( h + g * tau )
            a(j,q) = h + s * ( g - h * tau )
          end do

          do j = p + 1, q - 1
            g = a(p,j)
            h = a(j,q)
            a(p,j) = g - s * ( h + g * tau )
            a(j,q) = h + s * ( g - h * tau )
          end do

          do j = q + 1, n
            g = a(p,j)
            h = a(q,j)
            a(p,j) = g - s * ( h + g * tau )
            a(q,j) = h + s * ( g - h * tau )
          end do
!
!  Accumulate information in the eigenvector matrix.
!
          do j = 1, n
            g = v(j,p)
            h = v(j,q)
            v(j,p) = g - s * ( h + g * tau )
            v(j,q) = h + s * ( g - h * tau )
          end do

          rot_num = rot_num + 1

        end if

      end do
    end do

    bw(1:n) = bw(1:n) + zw(1:n)
    d(1:n) = bw(1:n)
    zw(1:n) = 0.0D+00

  end do
!
!  Restore upper triangle of input matrix.
!
  do j = 1, n
    do i = 1, j - 1
      a(i,j) = a(j,i)
    end do
  end do
!
!  Ascending sort the eigenvalues and eigenvectors.
!
  do k = 1, n - 1

    m = k

    do l = k + 1, n
      if ( d(l) < d(m) ) then
        m = l
      end if
    end do

    if ( m /= k ) then

      t    = d(m)
      d(m) = d(k)
      d(k) = t

      w(1:n)   = v(1:n,m)
      v(1:n,m) = v(1:n,k)
      v(1:n,k) = w(1:n)

    end if

  end do

  return
end
subroutine moment_method ( n, moment, x, w )

!*****************************************************************************80
!
!! MOMENT_METHOD computes a quadrature rule by the method of moments.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 September 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Gene Golub, John Welsch,
!    Calculation of Gaussian Quadrature Rules,
!    Mathematics of Computation,
!    Volume 23, Number 106, April 1969, pages 221-230.
!
!  Parameters:
!
!    Input, integer N, the order of the quadrature rule.
!
!    Input, real ( kind = rk ) MOMENT(2*N+1), moments 0 through 2*N.
!
!    Output, real ( kind = rk ) X(N), W(N), the points and weights of the 
!    quadrature rule.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ), allocatable :: alpha(:)
  real ( kind = rk ), allocatable :: beta(:)
  logical debug
  real ( kind = rk ) e
  integer flag
  real ( kind = rk ), allocatable :: h(:,:)
  integer i
  integer it_max
  integer it_num
  integer j
  real ( kind = rk ), allocatable :: jacobi(:,:)
  real ( kind = rk ) moment(0:2*n)
  real ( kind = rk ), allocatable :: r(:,:)
  real ( kind = rk ) r8mat_norm_fro
  integer rot_num
  real ( kind = rk ), allocatable :: t(:,:)
  real ( kind = rk ), allocatable :: v(:,:)
  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n)

  debug = .true.

  if ( debug ) then
    call r8vec_print ( 2 * n + 1, moment, '  Moments:' )
  end if
!
!  Define the N+1 by N+1 Hankel matrix H(I,J) = moment(I+J).
!
  allocate ( h(0:n,0:n) )

  do i = 0, n
    do j = 0, n
      h(i,j) = moment(i+j);
    end do
  end do

  if ( debug ) then
    call r8mat_print ( n + 1, n + 1, h, '  Hankel matrix H:' )
  end if
!
!  Compute R, the upper triangular Cholesky factor of H.
!
  allocate ( r(1:n+1,1:n+1) )

  call r8mat_cholesky_factor_upper ( n + 1, h, r, flag )

  if ( flag /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MOMENT_METHOD - Fatal error!'
    write ( *, '(a,i6)' ) '  R8MAT_CHOLESKY_FACTOR_UPPER returned FLAG = ', flag
    stop 1
  end if
  
  allocate ( t(0:n,0:n) )
  t = h - matmul ( transpose ( r ), r )
  e = r8mat_norm_fro ( n + 1, n + 1, t )
  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Froebenius norm H-R''*R = ', e
  deallocate ( t )

  if ( debug ) then
    call r8mat_print ( n + 1, n + 1, r, '  Cholesky factor R:' )
  end if
!
!  Compute ALPHA and BETA from R, using Golub and Welsch's formula.
!
  allocate ( alpha(1:n) )

  alpha(1) = r(1,2) / r(1,1)
  do i = 2, n
    alpha(i) = r(i,i+1) / r(i,i) - r(i-1,i) / r(i-1,i-1)
  end do

  allocate ( beta(1:n-1) )

  do i = 1, n - 1
    beta(i) = r(i+1,i+1) / r(i,i)
  end do
!
!  Compute the points and weights from the moments.
!
  allocate ( jacobi(1:n,1:n) )

  jacobi(1:n,1:n) = 0.0D+00

  do i = 1, n
    jacobi(i,i) = alpha(i)
  end do

  do i = 1, n - 1
    jacobi(i,i+1) = beta(i)
    jacobi(i+1,i) = beta(i)
  end do

  if ( debug ) then
    call r8mat_print ( n, n, jacobi, '  Jacobi matrix J:' )
  end if
!
!  Get the eigendecomposition of the Jacobi matrix.
!
  it_max = 100
  allocate ( v(1:n,1:n) )

  call jacobi_eigenvalue ( n, jacobi, it_max, v, x, it_num, rot_num )

  if ( debug ) then
    call r8mat_print ( n, n, v, '  Eigenvector matrix V:' )
  end if

  w(1:n) = moment(0) * v(1,1:n) ** 2
!
!  Free memory.
!
  deallocate ( alpha )
  deallocate ( beta )
  deallocate ( h )
  deallocate ( jacobi )
  deallocate ( r )
  deallocate ( v )

  return
end
subroutine moments_normal ( m, mu, sigma, w )

!*****************************************************************************80
!
!! MOMENTS_NORMAL returns moments of the standard Normal distribution.
!
!  Discussion:
!
!    pdf(x) = exp ( -((x-mu)/sigma)^2/2 ) / sigma / sqrt ( pi * 2 )
!    mu(k) = integral ( -oo < x < +oo ) x^k pdf(x) dx
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 September 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of moments desired.
!
!    Input, real ( kind = rk ) MU, SIGMA, the mean and standard deviation.
!
!    Output, real ( kind = rk ) W(0:M-1), the weighted integrals of X^0 
!    through X^(M-1).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m

  integer j
  integer j_hi
  integer k
  real ( kind = rk ) mu
  real ( kind = rk ) r8_choose
  real ( kind = rk ) r8_factorial2
  real ( kind = rk ) sigma
  real ( kind = rk ) t
  real ( kind = rk ) w(0:m-1)

  do k = 0, m - 1
    t = 0.0D+00
    j_hi = k / 2
    do j = 0, j_hi
      t = t + r8_choose ( k, 2 * j ) * r8_factorial2 ( 2 * j - 1 ) &
        * sigma ** ( 2 * j ) * mu ** ( k - 2 * j )
    end do
    w(k) = t
  end do

  return
end
subroutine moments_truncated_normal_ab ( m, mu, sigma, a, b, w )

!*****************************************************************************80
!
!! MOMENTS_TRUNCATED_NORMAL_AB: moments of truncated Normal distribution.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 September 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of moments desired.
!
!    Input, real ( kind = rk ) MU, SIGMA, the mean and standard deviation.
!
!    Input, real ( kind = rk ) A, B, the lower and upper truncation limits.
!
!    Output, real ( kind = rk ) W(0:M-1), the weighted integrals of X^0 
!    through X^(M-1).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) mu
  integer order
  real ( kind = rk ) sigma
  real ( kind = rk ) w(0:m-1)

  do order = 0, m - 1
    call truncated_normal_ab_moment ( order, mu, sigma, a, b, w(order) )
  end do

  return
end
subroutine moments_truncated_normal_a ( m, mu, sigma, a, w )

!*****************************************************************************80
!
!! MOMENTS_TRUNCATED_NORMAL_A: moments of lower truncated Normal distribution.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 September 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of moments desired.
!
!    Input, real ( kind = rk ) MU, SIGMA, the mean and standard deviation.
!
!    Input, real ( kind = rk ) A, the lower truncation limit.
!
!    Output, real ( kind = rk ) W(0:M-1), the weighted integrals of X^0 
!    through X^(M-1).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m

  real ( kind = rk ) a
  real ( kind = rk ) mu
  integer order
  real ( kind = rk ) sigma
  real ( kind = rk ) w(0:m-1)

  do order = 0, m - 1
    call truncated_normal_a_moment ( order, mu, sigma, a, w(order) )
  end do

  return
end
subroutine moments_truncated_normal_b ( m, mu, sigma, b, w )

!*****************************************************************************80
!
!! MOMENTS_TRUNCATED_NORMAL_B: moments of upper truncated Normal distribution.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 September 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of moments desired.
!
!    Input, real ( kind = rk ) MU, SIGMA, the mean and standard deviation.
!
!    Input, real ( kind = rk ) B, the upper truncation limit.
!
!    Output, real ( kind = rk ) W(0:M-1), the weighted integrals of X^0 
!    through X^(M-1).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m

  real ( kind = rk ) b
  real ( kind = rk ) mu
  integer order
  real ( kind = rk ) sigma
  real ( kind = rk ) w(0:m-1)

  do order = 0, m - 1
    call truncated_normal_b_moment ( order, mu, sigma, b, w(order) )
  end do

  return
end
subroutine normal_01_cdf ( x, cdf )

!*****************************************************************************80
!
!! NORMAL_01_CDF evaluates the Normal 01 CDF.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    AG Adams,
!    Algorithm 39,
!    Areas Under the Normal Curve,
!    Computer Journal,
!    Volume 12, pages 197-198, 1969.
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, the argument of the CDF.
!
!    Output, real ( kind = rk ) CDF, the value of the CDF.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), parameter :: a1 = 0.398942280444D+00
  real ( kind = rk ), parameter :: a2 = 0.399903438504D+00
  real ( kind = rk ), parameter :: a3 = 5.75885480458D+00
  real ( kind = rk ), parameter :: a4 = 29.8213557808D+00
  real ( kind = rk ), parameter :: a5 = 2.62433121679D+00
  real ( kind = rk ), parameter :: a6 = 48.6959930692D+00
  real ( kind = rk ), parameter :: a7 = 5.92885724438D+00
  real ( kind = rk ), parameter :: b0 = 0.398942280385D+00
  real ( kind = rk ), parameter :: b1 = 3.8052D-08
  real ( kind = rk ), parameter :: b2 = 1.00000615302D+00
  real ( kind = rk ), parameter :: b3 = 3.98064794D-04
  real ( kind = rk ), parameter :: b4 = 1.98615381364D+00
  real ( kind = rk ), parameter :: b5 = 0.151679116635D+00
  real ( kind = rk ), parameter :: b6 = 5.29330324926D+00
  real ( kind = rk ), parameter :: b7 = 4.8385912808D+00
  real ( kind = rk ), parameter :: b8 = 15.1508972451D+00
  real ( kind = rk ), parameter :: b9 = 0.742380924027D+00
  real ( kind = rk ), parameter :: b10 = 30.789933034D+00
  real ( kind = rk ), parameter :: b11 = 3.99019417011D+00
  real ( kind = rk ) cdf
  real ( kind = rk ) q
  real ( kind = rk ) x
  real ( kind = rk ) y
!
!  |X| <= 1.28.
!
  if ( abs ( x ) <= 1.28D+00 ) then

    y = 0.5D+00 * x * x

    q = 0.5D+00 - abs ( x ) * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5 &
      + a6 / ( y + a7 ) ) ) )
!
!  1.28 < |X| <= 12.7
!
  else if ( abs ( x ) <= 12.7D+00 ) then

    y = 0.5D+00 * x * x

    q = exp ( - y ) * b0 / ( abs ( x ) - b1 &
      + b2 / ( abs ( x ) + b3 &
      + b4 / ( abs ( x ) - b5 &
      + b6 / ( abs ( x ) + b7 &
      - b8 / ( abs ( x ) + b9 &
      + b10 / ( abs ( x ) + b11 ) ) ) ) ) )
!
!  12.7 < |X|
!
  else

    q = 0.0D+00

  end if
!
!  Take account of negative X.
!
  if ( x < 0.0D+00 ) then
    cdf = q
  else
    cdf = 1.0D+00 - q
  end if

  return
end
subroutine normal_01_pdf ( x, pdf )

!*****************************************************************************80
!
!! NORMAL_01_PDF evaluates the Normal 01 PDF.
!
!  Discussion:
!
!    The Normal 01 PDF is also called the "Standard Normal" PDF, or
!    the Normal PDF with 0 mean and variance 1.
!
!    PDF(X) = exp ( - 0.5 * X^2 ) / sqrt ( 2 * PI )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, the argument of the PDF.
!
!    Output, real ( kind = rk ) PDF, the value of the PDF.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) pdf
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) x

  pdf = exp ( -0.5D+00 * x * x ) / sqrt ( 2.0D+00 * pi )

  return
end
function r8_choose ( n, k )

!*****************************************************************************80
!
!! R8_CHOOSE computes the binomial coefficient C(N,K) as an R8.
!
!  Discussion:
!
!    The value is calculated in such a way as to avoid overflow and
!    roundoff.  The calculation is done in R8 arithmetic.
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
!    24 March 2008
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
!    Output, real ( kind = rk ) R8_CHOOSE, the number of combinations of N
!    things taken K at a time.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer k
  integer mn
  integer mx
  integer n
  real ( kind = rk ) r8_choose
  real ( kind = rk ) value

  mn = min ( k, n - k )

  if ( mn < 0 ) then

    value = 0.0D+00

  else if ( mn == 0 ) then

    value = 1.0D+00

  else

    mx = max ( k, n - k )
    value = real ( mx + 1, kind = rk )

    do i = 2, mn
      value = ( value * real ( mx + i, kind = rk ) ) / real ( i, kind = rk )
    end do

  end if

  r8_choose = value

  return
end
function r8_factorial ( n )

!*****************************************************************************80
!
!! R8_FACTORIAL computes the factorial of N.
!
!  Discussion:
!
!    factorial ( N ) = product ( 1 <= I <= N ) I
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the argument of the factorial function.
!    If N is less than 1, the function value is returned as 1.
!
!    Output, real ( kind = rk ) R8_FACTORIAL, the factorial of N.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r8_factorial
  integer i
  integer n

  r8_factorial = 1.0D+00

  do i = 1, n
    r8_factorial = r8_factorial * real ( i, kind = rk )
  end do

  return
end
function r8_factorial2 ( n )

!*****************************************************************************80
!
!! R8_FACTORIAL2 computes the double factorial function.
!
!  Discussion:
!
!    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
!                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
!
!  Example:
!
!     N Value
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
!    02 September 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the argument of the double factorial
!    function.  If N is less than 1, the value is returned as 1.0.
!
!    Output, real ( kind = rk ) R8_FACTORIAL2, the value.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  real ( kind = rk ) r8_factorial2
  real ( kind = rk ) r8_n

  if ( n < 1 ) then
    r8_factorial2 = 1.0D+00
    return
  end if

  r8_n = real ( n, kind = rk )
  r8_factorial2 = 1.0D+00

  do while ( 1.0D+00 < r8_n )
    r8_factorial2 = r8_factorial2 * r8_n
    r8_n = r8_n - 2.0D+00
  end do

  return
end
function r8_huge ( )

!*****************************************************************************80
!
!! R8_HUGE returns a very large R8.
!
!  Discussion:
!
!    The value returned by this function is NOT required to be the
!    maximum representable R8.  This value varies from machine to machine,
!    from compiler to compiler, and may cause problems when being printed.
!    We simply want a "very large" but non-infinite number.
!
!    FORTRAN90 provides a built-in routine HUGE ( X ) that
!    can return the maximum representable number of the same datatype
!    as X, if that is what is really desired.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 October 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) R8_HUGE, a "huge" value.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r8_huge

  r8_huge = 1.0D+30

  return
end
function r8_mop ( i )

!*****************************************************************************80
!
!! R8_MOP returns the I-th power of -1 as an R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = rk ) value.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, the power of -1.
!
!    Output, real ( kind = rk ) R8_MOP, the I-th power of -1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  real ( kind = rk ) r8_mop

  if ( mod ( i, 2 ) == 0 ) then
    r8_mop = + 1.0D+00
  else
    r8_mop = - 1.0D+00
  end if

  return
end
subroutine r8mat_cholesky_factor_upper ( n, a, c, flag )

!*****************************************************************************80
!
!! R8MAT_CHOLESKY_FACTOR_UPPER: upper Cholesky factor of a symmetric matrix.
!
!  Discussion:
!
!    The matrix must be symmetric and positive semidefinite.
!
!    For a positive semidefinite symmetric matrix A, the Cholesky factorization
!    is an upper triangular matrix R such that:
!
!      A = R * R'
!
!    The lower Cholesky factor is a lower triangular matrix L such that
!
!      A = L * L'
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of rows and columns of
!    the matrix A.
!
!    Input, real ( kind = rk ) A(N,N), the N by N matrix.
!
!    Output, real ( kind = rk ) C(N,N), the N by N upper triangular
!    Cholesky factor.
!
!    Output, integer FLAG:
!    0, no error occurred.
!    1, the matrix is not positive definite.
!    2, the matrix is not nonnegative definite.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) c(n,n)
  integer flag
  integer i
  integer j
  real ( kind = rk ) sum2

  flag = 0

  c(1:n,1:n) = a(1:n,1:n)

  do j = 1, n

    c(j,1:j-1) = 0.0D+00

    do i = j, n

      sum2 = c(i,j) - dot_product ( c(1:j-1,j), c(1:j-1,i) )

      if ( i == j ) then
        if ( sum2 <= 0.0D+00 ) then
          flag = 1
          return
        else
          c(j,i) = sqrt ( sum2 )
        end if
      else
        if ( c(j,j) /= 0.0D+00 ) then
          c(j,i) = sum2 / c(j,j)
        else
          c(j,i) = 0.0D+00
        end if
      end if

    end do

  end do

  return
end
function r8mat_norm_fro ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_NORM_FRO returns the Frobenius norm of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    The Frobenius norm is defined as
!
!      R8MAT_NORM_FRO = sqrt (
!        sum ( 1 <= I <= M ) sum ( 1 <= j <= N ) A(I,J) * A(I,J) )
!
!    The matrix Frobenius norm is not derived from a vector norm, but
!    is compatible with the vector L2 norm, so that:
!
!      r8vec_norm_l2 ( A * x ) <= r8mat_norm_fro ( A ) * r8vec_norm_l2 ( x ).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 March 2000
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
!    Input, real ( kind = rk ) A(M,N), the matrix whose Frobenius
!    norm is desired.
!
!    Output, real ( kind = rk ) R8MAT_NORM_FRO, the Frobenius norm of A.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) r8mat_norm_fro

  r8mat_norm_fro = sqrt ( sum ( a(1:m,1:n)**2 ) )

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
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
!! R8MAT_PRINT_SOME prints some of an R8MAT.
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
subroutine r8mat_write ( output_filename, m, n, table )

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
!    Input, character ( len = * ) OUTPUT_FILENAME, the output file name.
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
  character ( len = * ) output_filename
  integer output_status
  integer output_unit
  character ( len = 30 ) string
  real ( kind = rk ) table(m,n)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_filename ) // '" on unit ', output_unit
    output_unit = -1
    stop
  end if
!
!  Create a format string.
!
  write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 24, '.', 16, ')'
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
    write ( *, '(2x,i8,a,1x,g24.16)' ) i, ':', a(i)
  end do

  return
end
subroutine rule_write ( n, x, w, r, filename )

!*****************************************************************************80
!
!! RULE_WRITE writes a quadrature rule to a file.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the rule.
!
!    Input, real ( kind = rk ) X(N), the abscissas.
!
!    Input, real ( kind = rk ) W(N), the weights.
!
!    Input, real ( kind = rk ) R(2), defines the region.
!
!    Input, character ( len = * ) FILENAME, specifies the output.
!    'filename_w.txt', 'filename_x.txt', 'filename_r.txt' defining weights,
!    abscissas, and region.
! 
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer  n

  character ( len = * ) filename
  character ( len = 255 ) filename_r
  character ( len = 255 ) filename_w
  character ( len = 255 ) filename_x
  real ( kind = rk ) r(2)
  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n)

  filename_w = trim ( filename ) // '_w.txt'
  filename_x = trim ( filename ) // '_x.txt'
  filename_r = trim ( filename ) // '_r.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Creating quadrature files.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  "Root" file name is   "' // trim ( filename ) // '".'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Weight file will be   "' // trim ( filename_w ) // '".'
  write ( *, '(a)' ) '  Abscissa file will be "' // trim ( filename_x ) // '".'
  write ( *, '(a)' ) '  Region file will be   "' // trim ( filename_r ) // '".'
            
  call r8mat_write ( filename_w, 1, n, w )
  call r8mat_write ( filename_x, 1, n, x )
  call r8mat_write ( filename_r, 1, 2, r )

  return
end
subroutine s_to_i4 ( s, ival, ierror, length )

!*****************************************************************************80
!
!! S_TO_I4 reads an I4 from a string.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 January 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string to be examined.
!
!    Output, integer IVAL, the integer value read from the string.
!    If the string is blank, then IVAL will be returned 0.
!
!    Output, integer IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer LENGTH, the number of characters of S 
!    used to make IVAL.
!
  implicit none

  character c
  integer i
  integer ierror
  integer isgn
  integer istate
  integer ival
  integer length
  character ( len = * ) s

  ierror = 0
  istate = 0
  isgn = 1
  ival = 0

  do i = 1, len_trim ( s )

    c = s(i:i)
!
!  Haven't read anything.
!
    if ( istate == 0 ) then

      if ( c == ' ' ) then

      else if ( c == '-' ) then
        istate = 1
        isgn = -1
      else if ( c == '+' ) then
        istate = 1
        isgn = + 1
      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read the sign, expecting digits.
!
    else if ( istate == 1 ) then

      if ( c == ' ' ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read at least one digit, expecting more.
!
    else if ( istate == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        ival = 10 * ival + ichar ( c ) - ichar ( '0' )
      else
        ival = isgn * ival
        length = i - 1
        return
      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( istate == 2 ) then
    ival = isgn * ival
    length = len_trim ( s )
  else
    ierror = 1
    length = 0
  end if

  return
end
subroutine s_to_r8 ( s, dval, ierror, length )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 value from a string.
!
!  Discussion:
!
!    An "R8" value is simply a real number to be stored as a
!    variable of type "real ( kind = rk )".
!
!    The routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 blanks
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semicolon,
!
!    with most quantities optional.
!
!  Example:
!
!    S                 DVAL
!
!    '1'               1.0
!    '     1   '       1.0
!    '1A'              1.0
!    '12,34,56'        12.0
!    '  34 7'          34.0
!    '-1E2ABCD'        -100.0
!    '-1X2ABCD'        -1.0
!    ' 2E-1'           0.2
!    '23.45'           23.45
!    '-4.2E+2'         -420.0
!    '17d2'            1700.0
!    '-14e-2'         -0.14
!    'e2'              100.0
!    '-12.73e-9.23'   -12.73 * 10.0^(-9.23)
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string containing the
!    data to be read.  Reading will begin at position 1 and
!    terminate at the end of the string, or when no more
!    characters can be read to form a legal real.  Blanks,
!    commas, or other nonnumeric data will, in particular,
!    cause the conversion to halt.
!
!    Output, real ( kind = rk ) DVAL, the value read from the string.
!
!    Output, integer IERROR, error flag.
!    0, no errors occurred.
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer LENGTH, the number of characters read
!    to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character c
  logical ch_eqi
  real ( kind = rk ) dval
  integer ierror
  integer ihave
  integer isgn
  integer iterm
  integer jbot
  integer jsgn
  integer jtop
  integer length
  integer ndig
  real ( kind = rk ) rbot
  real ( kind = rk ) rexp
  real ( kind = rk ) rtop
  character ( len = * ) s
  integer s_length
  character :: TAB = achar ( 9 )

  s_length = len_trim ( s )

  ierror = 0
  dval = 0.0D+00
  length = -1
  isgn = 1
  rtop = 0
  rbot = 1
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

  do

    length = length + 1

    if ( s_length < length+1 ) then
      exit
    end if

    c = s(length+1:length+1)
!
!  Blank character.
!
    if ( c == ' ' .or. c == TAB ) then

      if ( ihave == 2 ) then

      else if ( ihave == 6 .or. ihave == 7 ) then
        iterm = 1
      else if ( 1 < ihave ) then
        ihave = 11
      end if
!
!  Comma.
!
    else if ( c == ',' .or. c == ';' ) then

      if ( ihave /= 1 ) then
        iterm = 1
        ihave = 12
        length = length + 1
      end if
!
!  Minus sign.
!
    else if ( c == '-' ) then

      if ( ihave == 1 ) then
        ihave = 2
        isgn = -1
      else if ( ihave == 6 ) then
        ihave = 7
        jsgn = -1
      else
        iterm = 1
      end if
!
!  Plus sign.
!
    else if ( c == '+' ) then

      if ( ihave == 1 ) then
        ihave = 2
      else if ( ihave == 6 ) then
        ihave = 7
      else
        iterm = 1
      end if
!
!  Decimal point.
!
    else if ( c == '.' ) then

      if ( ihave < 4 ) then
        ihave = 4
      else if ( 6 <= ihave .and. ihave <= 8 ) then
        ihave = 9
      else
        iterm = 1
      end if
!
!  Scientific notation exponent marker.
!
    else if ( ch_eqi ( c, 'E' ) .or. ch_eqi ( c, 'D' ) ) then

      if ( ihave < 6 ) then
        ihave = 6
      else
        iterm = 1
      end if
!
!  Digit.
!
    else if (  ihave < 11 .and. lle ( '0', c ) .and. lle ( c, '9' ) ) then

      if ( ihave <= 2 ) then
        ihave = 3
      else if ( ihave == 4 ) then
        ihave = 5
      else if ( ihave == 6 .or. ihave == 7 ) then
        ihave = 8
      else if ( ihave == 9 ) then
        ihave = 10
      end if

      call ch_to_digit ( c, ndig )

      if ( ihave == 3 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = rk )
      else if ( ihave == 5 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = rk )
        rbot = 10.0D+00 * rbot
      else if ( ihave == 8 ) then
        jtop = 10 * jtop + ndig
      else if ( ihave == 10 ) then
        jtop = 10 * jtop + ndig
        jbot = 10 * jbot
      end if
!
!  Anything else is regarded as a terminator.
!
    else
      iterm = 1
    end if
!
!  If we haven't seen a terminator, and we haven't examined the
!  entire string, go get the next character.
!
    if ( iterm == 1 ) then
      exit
    end if

  end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LENGTH is equal to S_LENGTH.
!
  if ( iterm /= 1 .and. length+1 == s_length ) then
    length = s_length
  end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
  if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then
    ierror = ihave
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'S_TO_R8 - Serious error!'
    write ( *, '(a)' ) '  Illegal or nonnumeric input:'
    write ( *, '(a)' ) '    ' // trim ( s )
    return
  end if
!
!  Number seems OK.  Form it.
!
  if ( jtop == 0 ) then
    rexp = 1.0D+00
  else
    if ( jbot == 1 ) then
      rexp = 10.0D+00 ** ( jsgn * jtop )
    else
      rexp = 10.0D+00 ** ( real ( jsgn * jtop, kind = rk ) &
        / real ( jbot, kind = rk ) )
    end if
  end if

  dval = real ( isgn, kind = rk ) * rexp * rtop / rbot

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
subroutine truncated_normal_ab_moment ( order, mu, s, a, b, moment )

!*****************************************************************************80
!
!! TRUNCATED_NORMAL_AB_MOMENT: moments of the truncated Normal PDF.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 September 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Phoebus Dhrymes,
!    Moments of Truncated Normal Distributions,
!    May 2005.
!
!  Parameters:
!
!    Input, integer ORDER, the order of the moment.
!    0 <= ORDER.
!
!    Input, real ( kind = rk ) MU, S, the mean and standard deviation of the
!    parent Normal distribution.
!    0.0 < S.
!
!    Input, real ( kind = rk ) A, B, the lower and upper truncation limits.
!    A < B.
!
!    Output, real ( kind = rk ) MOMENT, the moment of the PDF.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) a_h
  real ( kind = rk ) a_cdf
  real ( kind = rk ) a_pdf
  real ( kind = rk ) b
  real ( kind = rk ) b_h
  real ( kind = rk ) b_cdf
  real ( kind = rk ) b_pdf
  real ( kind = rk ) ir
  real ( kind = rk ) irm1
  real ( kind = rk ) irm2
  real ( kind = rk ) moment
  real ( kind = rk ) mu
  integer order
  integer r
  real ( kind = rk ) r8_choose
  real ( kind = rk ) s

  if ( order < 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'TRUNCATED_NORMAL_AB_MOMENT - Fatal error!'
    write ( *, '(a)' ) '  ORDER < 0.'
    stop 1
  end if

  if ( s <= 0.0D+00 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'TRUNCATED_NORMAL_AB_MOMENT - Fatal error!'
    write ( *, '(a)' ) '  S <= 0.0.'
    stop 1
  end if

  if ( b <= a ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'TRUNCATED_NORMAL_AB_MOMENT - Fatal error!'
    write ( *, '(a)' ) '  B <= A.'
    stop 1
  end if

  a_h = ( a - mu ) / s
  call normal_01_pdf ( a_h, a_pdf )
  call normal_01_cdf ( a_h, a_cdf )

  if ( a_cdf == 0.0D+00 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'TRUNCATED_NORMAL_AB_MOMENT - Fatal error!'
    write ( *, '(a)' ) '  PDF/CDF ratio fails, because A_CDF is too small.'
    write ( *, '(a,g14.6)' ) '  A_PDF = ', a_pdf
    write ( *, '(a,g14.6)' ) '  A_CDF = ', a_cdf
    stop 1
  end if

  b_h = ( b - mu ) / s
  call normal_01_pdf ( b_h, b_pdf )
  call normal_01_cdf ( b_h, b_cdf )

  if ( b_cdf == 0.0D+00 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'TRUNCATED_NORMAL_AB_MOMENT - Fatal error!'
    write ( *, '(a)' ) '  PDF/CDF ratio fails, because B_CDF is too small.'
    write ( *, '(a,g14.6)' ) '  B_PDF = ', b_pdf
    write ( *, '(a,g14.6)' ) '  B_CDF = ', b_cdf
    stop 1
  end if

  moment = 0.0D+00
  irm2 = 0.0D+00
  irm1 = 0.0D+00

  do r = 0, order

    if ( r == 0 ) then
      ir = 1.0D+00
    else if ( r == 1 ) then
      ir = - ( b_pdf - a_pdf ) / ( b_cdf - a_cdf )
    else
      ir = real ( r - 1, kind = rk ) * irm2 &
        - ( b_h ** ( r - 1 ) * b_pdf - a_h ** ( r - 1 ) * a_pdf ) &
        / ( b_cdf - a_cdf )
    end if

    moment = moment + r8_choose ( order, r ) * mu ** ( order - r ) &
      * ( s ** r ) * ir

    irm2 = irm1
    irm1 = ir

  end do

  return
end
subroutine truncated_normal_a_moment ( order, mu, s, a, moment )

!*****************************************************************************80
!
!! TRUNCATED_NORMAL_A_MOMENT: moments of the lower truncated Normal PDF.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 September 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Phoebus Dhrymes,
!    Moments of Truncated Normal Distributions,
!    May 2005.
!
!  Parameters:
!
!    Input, integer ORDER, the order of the moment.
!    0 <= ORDER.
!
!    Input, real ( kind = rk ) MU, S, the mean and standard deviation of the
!    parent Normal distribution.
!    0.0 < S.
!
!    Input, real ( kind = rk ) A, the lower truncation limit.
!
!    Output, real ( kind = rk ) MOMENT, the moment of the PDF.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) moment
  real ( kind = rk ) mu
  integer order
  real ( kind = rk ) r8_mop
  real ( kind = rk ) s

  call truncated_normal_b_moment ( order, - mu, s, - a, moment )
  moment = r8_mop ( order ) * moment

  return
end
subroutine truncated_normal_b_moment ( order, mu, s, b, moment )

!*****************************************************************************80
!
!! TRUNCATED_NORMAL_B_MOMENT: moments of the upper truncated Normal PDF.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 September 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Phoebus Dhrymes,
!    Moments of Truncated Normal Distributions,
!    May 2005.
!
!  Parameters:
!
!    Input, integer ORDER, the order of the moment.
!    0 <= ORDER.
!
!    Input, real ( kind = rk ) MU, S, the mean and standard deviation of the
!    parent Normal distribution.
!    0.0 < S.
!
!    Input, real ( kind = rk ) B, the upper truncation limit.
!
!    Output, real ( kind = rk ) MOMENT, the moment of the PDF.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) b
  logical debug
  real ( kind = rk ) f
  real ( kind = rk ) h
  real ( kind = rk ) h_cdf
  real ( kind = rk ) h_pdf
  real ( kind = rk ) ir
  real ( kind = rk ) irm1
  real ( kind = rk ) irm2
  real ( kind = rk ) moment
  real ( kind = rk ) mu
  integer order
  integer r
  real ( kind = rk ) r8_choose
  real ( kind = rk ) s

  debug = .true.
  
  if ( order < 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'TRUNCATED_NORMAL_B_MOMENT - Fatal error!'
    write ( *, '(a)' ) '  ORDER < 0.'
    stop 1
  end if

  if ( debug ) then
    write ( *,  '(a,i4,a,g14.6,a,g14.6,a,g14.6)') &
      '  ORDER = ', order, ',  b = ', b, ',  MU = ', mu, ',  S = ', s
  end if
  
  h = ( b - mu ) / s
  call normal_01_pdf ( h, h_pdf )
  call normal_01_cdf ( h, h_cdf )

  if ( debug ) then
    write ( *,  '(a,i4,a,g14.6,a,g14.6,a,g14.6)') &
      '  ORDER = ', order, ',  H = ', h, ',  H_PDF = ', h_pdf, ',  H_CDF = ', h_cdf
  end if
  
  if ( h_cdf == 0.0D+00 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'TRUNCATED_NORMAL_B_MOMENT - Fatal error!'
    write ( *, '(a)' ) '  CDF((B-MU)/S) = 0.'
    stop 1
  end if

  f = h_pdf / h_cdf

  moment = 0.0D+00
  irm2 = 0.0D+00
  irm1 = 0.0D+00

  do r = 0, order

    if ( r == 0 ) then
      ir = 1.0D+00
    else if ( r == 1 ) then
      ir = - f
    else
      ir = - h ** ( r - 1 ) * f + real ( r - 1, kind = rk ) * irm2
    end if

    moment = moment + r8_choose ( order, r ) * mu ** ( order - r ) &
      * ( s ** r ) * ir

    irm2 = irm1
    irm1 = ir

  end do

  if ( debug ) then
    write ( *, '(a,g14.6)' ) '  MOMENT = ', moment
  end if
  
  return
end
