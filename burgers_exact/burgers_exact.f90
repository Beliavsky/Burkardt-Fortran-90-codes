subroutine burgers_viscous_time_exact1 ( nu, vxn, vx, vtn, vt, vu )

!*****************************************************************************80
!
!! burgers_viscous_time_exact1() evaluates solution #1 to the Burgers equation.
!
!  Discussion:
!
!    The form of the Burgers equation considered here is
!
!      du       du        d^2 u
!      -- + u * -- = nu * -----
!      dt       dx        dx^2
!
!    for -1.0 < x < +1.0, and 0 < t.
!
!    Initial conditions are u(x,0) = - sin(pi*x).  Boundary conditions
!    are u(-1,t) = u(+1,t) = 0.  The viscosity parameter nu is taken
!    to be 0.01 / pi, although this is not essential.
!
!    The authors note an integral representation for the solution u(x,t),
!    and present a better version of the formula that is amenable to
!    approximation using Hermite quadrature.  
!
!    This program library does little more than evaluate the exact solution
!    at a user-specified set of points, using the quadrature rule.
!    Internally, the order of this quadrature rule is set to 8, but the
!    user can easily modify this value if greater accuracy is desired.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Claude Basdevant, Michel Deville, Pierre Haldenwang, J Lacroix, 
!    J Ouazzani, Roger Peyret, Paolo Orlandi, Anthony Patera,
!    Spectral and finite difference solutions of the Burgers equation,
!    Computers and Fluids,
!    Volume 14, Number 1, 1986, pages 23-41.
!
!  Input:
!
!    real ( kind = rk ) NU, the viscosity.
!
!    integer VXN, the number of spatial grid points.
!
!    real ( kind = rk ) VX(VXN), the spatial grid points.
!
!    integer VTN, the number of time grid points.
!
!    real ( kind = rk ) VT(VTN), the time grid points.
!
!    real ( kind = rk ) VU(VXN,VTN), the solution of the Burgers
!    equation at each space and time grid point.
!
  implicit none
  
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: qn = 8
  integer vtn
  integer vxn

  real ( kind = rk ) bot
  real ( kind = rk ) c
  real ( kind = rk ) nu
  integer qi
  real ( kind = rk ) qw(qn)
  real ( kind = rk ) qx(qn)
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) vt(vtn)
  integer vti
  real ( kind = rk ) vx(vxn)
  integer vxi
  real ( kind = rk ) vu(vxn,vtn)
  real ( kind = rk ) top
!
!  Compute the rule.
!
  call hermite_ek_compute ( qn, qx, qw )
!
!  Evaluate U(X,T) for later times.
!
  do vti = 1, vtn

    if ( vt(vti) == 0.0D+00 ) then

      do vxi = 1, vxn
        vu(vxi,vti) = - sin ( r8_pi * vx(vxi) )
      end do

    else

      do vxi = 1, vxn

        top = 0.0D+00
        bot = 0.0D+00

        do qi = 1, qn

          c = 2.0D+00 * sqrt ( nu * vt(vti) )

          top = top - qw(qi) * c * sin ( r8_pi * ( vx(vxi) - c * qx(qi) ) ) &
            * exp ( - cos ( r8_pi * ( vx(vxi) - c * qx(qi)  ) ) &
            / ( 2.0D+00 * r8_pi * nu ) )

          bot = bot + qw(qi) * c &
            * exp ( - cos ( r8_pi * ( vx(vxi) - c * qx(qi)  ) ) &
            / ( 2.0D+00 * r8_pi * nu ) )

          vu(vxi,vti) = top / bot

        end do

      end do

    end if

  end do
  
  return
end
subroutine burgers_viscous_time_exact2 ( nu, xn, x, tn, t, u )

!*****************************************************************************80
!
!! burgers_viscous_time_exact2() evaluates solution #2 to the Burgers equation.
!
!  Discussion:
!
!    The form of the Burgers equation considered here is
!
!      du       du        d^2 u
!      -- + u * -- = nu * -----
!      dt       dx        dx^2
!
!    for 0.0 < x < 2 Pi and 0 < t.
!
!    The initial condition is
!
!      u(x,0) = 4 - 2 * nu * dphi(x,0)/dx / phi(x,0)
!
!    where
!
!      phi(x,t) = exp ( - ( x-4*t      ) / ( 4*nu*(t+1) ) )
!               + exp ( - ( x-4*t-2*pi ) / ( 4*nu*(t+1) ) )
!
!    The boundary conditions are periodic:
!
!      u(0,t) = u(2 Pi,t)
!
!    The viscosity parameter nu may be taken to be 0.01, but other values
!    may be chosen.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt.
!
!  Input:
!
!    real ( kind = rk ) NU, the viscosity.
!
!    integer XN, the number of spatial grid points.
!
!    real ( kind = rk ) X(XN), the spatial grid points.
!
!    integer TN, the number of time grid points.
!
!    real ( kind = rk ) T(TN), the time grid points.
!
!    real ( kind = rk ) U(XN,TN), the solution of the Burgers
!    equation at each space and time grid point.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer tn
  integer xn

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) dphi
  integer i
  integer j
  real ( kind = rk ) nu
  real ( kind = rk ) phi
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) t(tn)
  real ( kind = rk ) u(xn,tn)
  real ( kind = rk ) x(xn)

  do j = 1, tn

    do i = 1, xn

      a = ( x(i) - 4.0D+00 * t(j) )
      b = ( x(i) - 4.0D+00 * t(j) - 2.0D+00 * r8_pi )
      c = 4.0D+00 * nu * ( t(j) + 1.0D+00 )
      phi = exp ( - a ** 2 / c ) + exp ( - b ** 2 / c )
      dphi = - 2.0D+00 * a * exp ( - a ** 2 / c ) / c &
             - 2.0D+00 * b * exp ( - b ** 2 / c ) / c
      u(i,j) = 4.0D+00 - 2.0D+00 * nu * dphi / phi

    end do

  end do
  
  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! get_unit() returns a free Fortran unit number.
!
!  Discussion:
!
!    A "free" Fortran unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free Fortran unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free Fortran unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free Fortran unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Output:
!
!    integer IUNIT, the free unit number.
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
subroutine hermite_ek_compute ( n, x, w )

!*****************************************************************************80
!
!! hermite_ek_compute() computes a Gauss-Hermite quadrature rule.
!
!  Discussion:
!
!    The code uses an algorithm by Elhay and Kautsky.
!
!    The abscissas are the zeros of the N-th order Hermite polynomial.
!
!    The integral:
!
!      integral ( -oo < x < +oo ) exp ( - x * x ) * f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    Original Fortran77 version by Sylvan Elhay, Jaroslav Kautsky.
!    Fortran90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, Fortran Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Input:
!
!    integer N, the number of abscissas.
!
!  Output:
!
!    real ( kind = rk ) X(N), the abscissas.
!
!    real ( kind = rk ) W(N), the weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) bj(n)
  integer i
  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n)
  real ( kind = rk ) zemu
!
!  Define the zero-th moment.
!
  zemu = gamma ( 1.0D+00 / 2.0D+00 )
!
!  Define the Jacobi matrix.
!
  do i = 1, n
    bj(i) = real ( i, kind = rk ) / 2.0D+00
  end do
  bj(1:n) = sqrt ( bj(1:n) )

  x(1:n) = 0.0D+00

  w(1) = sqrt ( zemu )
  w(2:n) = 0.0D+00
!
!  Diagonalize the Jacobi matrix.
!
  call imtqlx ( n, x, bj, w )

  w(1:n) = w(1:n)**2

  return
end
subroutine imtqlx ( n, d, e, z )

!*****************************************************************************80
!
!! imtqlx() diagonalizes a symmetric tridiagonal matrix.
!
!  Discussion:
!
!    This routine is a slightly modified version of the EISPACK routine to
!    perform the implicit QL algorithm on a symmetric tridiagonal matrix.
!
!    The authors thank the authors of EISPACK for permission to use this
!    routine.
!
!    It has been modified to produce the product Q' * Z, where Z is an input
!    vector and Q is the orthogonal matrix diagonalizing the input matrix.
!    The changes consist (essentially) of applying the orthogonal
!    transformations directly to Z as they are generated.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    Original Fortran77 version by Sylvan Elhay, Jaroslav Kautsky.
!    Fortran90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, Fortran Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!    Roger Martin, James Wilkinson,
!    The Implicit QL Algorithm,
!    Numerische Mathematik,
!    Volume 12, Number 5, December 1968, pages 377-383.
!
!  Input:
!
!    integer N, the order of the matrix.
!
!    real ( kind = rk ) D(N), the diagonal entries of the matrix.
!
!    real ( kind = rk ) E(N), the subdiagonal entries of the
!    matrix, in entries E(1) through E(N-1).  
!
!    real ( kind = rk ) Z(N), a vector.
!
!  Output:
!
!    real ( kind = rk ) D(N), overwritten.
!
!    real ( kind = rk ) E(N), overwritten.
!
!    real ( kind = rk ) Z(N), the value of Q' * Z, where Q is the matrix
!    that diagonalizes the input symmetric tridiagonal matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) d(n)
  real ( kind = rk ) e(n)
  real ( kind = rk ) f
  real ( kind = rk ) g
  integer i
  integer ii
  integer, parameter :: itn = 30
  integer j
  integer k
  integer l
  integer m
  integer mml
  real ( kind = rk ) p
  real ( kind = rk ) prec
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) z(n)

  prec = epsilon ( prec )

  if ( n == 1 ) then
    return
  end if

  e(n) = 0.0D+00

  do l = 1, n

    j = 0

    do

      do m = l, n

        if ( m == n ) then
          exit
        end if

        if ( abs ( e(m) ) <= prec * ( abs ( d(m) ) + abs ( d(m+1) ) ) ) then
          exit
        end if

      end do

      p = d(l)

      if ( m == l ) then
        exit
      end if

      if ( itn <= j ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IMTQLX - Fatal error!'
        write ( *, '(a)' ) '  Iteration limit exceeded.'
        write ( *, '(a,i8)' ) '  J = ', j
        write ( *, '(a,i8)' ) '  L = ', l
        write ( *, '(a,i8)' ) '  M = ', m
        write ( *, '(a,i8)' ) '  N = ', n
        stop
      end if

      j = j + 1
      g = ( d(l+1) - p ) / ( 2.0D+00 * e(l) )
      r =  sqrt ( g * g + 1.0D+00 )
      g = d(m) - p + e(l) / ( g + sign ( r, g ) )
      s = 1.0D+00
      c = 1.0D+00
      p = 0.0D+00
      mml = m - l

      do ii = 1, mml

        i = m - ii
        f = s * e(i)
        b = c * e(i)

        if ( abs ( g ) <= abs ( f ) ) then
          c = g / f
          r =  sqrt ( c * c + 1.0D+00 )
          e(i+1) = f * r
          s = 1.0D+00 / r
          c = c * s
        else
          s = f / g
          r =  sqrt ( s * s + 1.0D+00 )
          e(i+1) = g * r
          c = 1.0D+00 / r
          s = s * c
        end if

        g = d(i+1) - p
        r = ( d(i) - g ) * s + 2.0D+00 * c * b
        p = s * r
        d(i+1) = g + p
        g = c * r - b
        f = z(i+1)
        z(i+1) = s * z(i) + c * f
        z(i) = c * z(i) - s * f

      end do

      d(l) = d(l) - p
      e(l) = g
      e(m) = 0.0D+00

    end do

  end do
!
!  Sorting.
!
  do ii = 2, n

    i = ii - 1
    k = i
    p = d(i)

    do j = ii, n
      if ( d(j) < p ) then
        k = j
        p = d(j)
      end if
    end do

    if ( k /= i ) then
      d(k) = d(i)
      d(i) = p
      p = z(i)
      z(i) = z(k)
      z(k) = p
    end if

  end do

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! r8mat_print() prints an R8MAT.
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
!    12 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, the number of rows in A.
!
!    integer N, the number of columns in A.
!
!    real ( kind = rk ) A(M,N), the matrix.
!
!    character ( len = * ) TITLE, a title.
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
!  Input:
!
!    integer M, N, the number of rows and columns.
!
!    real ( kind = rk ) A(M,N), an M by N matrix to be printed.
!
!    integer ILO, JLO, the first row and column to print.
!
!    integer IHI, JHI, the last row and column to print.
!
!    character ( len = * ) TITLE, a title.
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
!! r8mat_write() writes an R8MAT file.
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
!    31 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    character ( len = * ) OUTPUT_FILENAME, the output file name.
!
!    integer M, the spatial dimension.
!
!    integer N, the number of points.
!
!    real ( kind = rk ) TABLE(M,N), the data.
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
!  For less precision in the output file, try:
!
!                                            '(', m, 'g', 14, '.', 6, ')'
!
  if ( 0 < m .and. 0 < n ) then

    write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 24, '.', 16, ')'
!
!  Write the data.
!
    do j = 1, n
      write ( output_unit, string ) table(1:m,j)
    end do

  end if
!
!  Close the file.
!
  close ( unit = output_unit )

  return
end
subroutine r8vec_even ( n, alo, ahi, a )

!*****************************************************************************80
!
!! r8vec_even() returns an R8VEC of evenly spaced values.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    If N is 1, then the midpoint is returned.
!
!    Otherwise, the two endpoints are returned, and N-2 evenly
!    spaced points between them.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of values.
!
!    real ( kind = rk ) ALO, AHI, the low and high values.
!
!  Output:
!
!    real ( kind = rk ) A(N), N evenly spaced values.
!    Normally, A(1) = ALO and A(N) = AHI.
!    However, if N = 1, then A(1) = 0.5*(ALO+AHI).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) ahi
  real ( kind = rk ) alo
  integer i

  if ( n == 1 ) then

    a(1) = 0.5D+00 * ( alo + ahi )

  else

    do i = 1, n

      a(i) = ( real ( n - i,     kind = rk ) * alo   &
             + real (     i - 1, kind = rk ) * ahi ) &
             / real ( n     - 1, kind = rk )
    end do

  end if

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
!  Input:
!
!    integer N, the number of components of the vector.
!
!    real ( kind = rk ) A(N), the vector to be printed.
!
!    character ( len = * ) TITLE, a title.
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
!    06 August 2005
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
 
