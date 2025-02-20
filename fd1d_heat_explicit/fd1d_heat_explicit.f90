subroutine fd1d_heat_explicit ( x_num, x, t, dt, cfl, rhs, bc, h, h_new )

!*****************************************************************************80
!
!! fd1d_heat_explicit(): Finite difference solution of 1D heat equation.
!
!  Discussion:
!
!    This program takes one time step to solve the 1D heat equation 
!    with an explicit method.
!
!    This program solves
!
!      dUdT - k * d2UdX2 = F(X,T)
!
!    over the interval [A,B] with boundary conditions
!
!      U(A,T) = UA(T),
!      U(B,T) = UB(T),
!
!    over the time interval [T0,T1] with initial conditions
!
!      U(X,T0) = U0(X)
!
!    The code uses the finite difference method to approximate the
!    second derivative in space, and an explicit forward Euler approximation
!    to the first derivative in time.
!
!    The finite difference form can be written as
!
!      U(X,T+dt) - U(X,T)                  ( U(X-dx,T) - 2 U(X,T) + U(X+dx,T) )
!      ------------------  = F(X,T) + k *  ------------------------------------
!               dt                                   dx * dx
!
!    or, assuming we have solved for all values of U at time T, we have
!
!      U(X,T+dt) = U(X,T) 
!        + cfl * ( U(X-dx,T) - 2 U(X,T) + U(X+dx,T) ) + dt * F(X,T) 
!
!    Here "cfl" is the Courant-Friedrichs-Loewy coefficient:
!
!      cfl = k * dt / dx / dx
!
!    In order for accurate results to be computed by this explicit method,
!    the CFL coefficient must be less than 0.5.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 January 2022
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer X_NUM, the number of points to use in the 
!    spatial dimension.
!
!    real ( kind = rk ) X(X_NUM), the coordinates of the nodes.
!
!    real ( kind = rk ) T, the current time.
!
!    real ( kind = rk ) DT, the size of the time step.
!
!    real ( kind = rk ) CFL, the Courant-Friedrichs-Loewy coefficient,
!    computed by FD1D_HEAT_EXPLICIT_CFL.
!
!    real ( kind = rk ) H(X_NUM), the solution at the current time.
!
!    external RHS, the function which evaluates the right hand side.
!
!    external BC, the function which evaluates the boundary conditions.
!
!  Output:
!
!    real ( kind = rk ) H_NEW(X_NUM), the solution at time T+DT.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer x_num

  real ( kind = rk ) cfl
  real ( kind = rk ) f(x_num)
  real ( kind = rk ) dt
  real ( kind = rk ) h(x_num)
  real ( kind = rk ) h_new(x_num)
  real ( kind = rk ) t
  real ( kind = rk ) x(x_num)

  call rhs ( x_num, x, t, f )

  h_new(1) = 0.0D+00

  h_new(2:x_num-1) = h(2:x_num-1) + dt * f(2:x_num-1) &
    + cfl * (             h(1:x_num-2) &
              - 2.0D+00 * h(2:x_num-1) &
              +           h(3:x_num) )

  h_new(x_num) = 0.0D+00

  call bc ( x_num, x, t + dt, h_new )

  return
end
subroutine fd1d_heat_explicit_cfl ( k, t_num, t_min, t_max, x_num, x_min, &
  x_max, cfl )

!*****************************************************************************80
!
!! fd1d_heat_explicit_cfl() computes the Courant-Friedrichs-Loewy coefficient.
!
!  Discussion:
!
!    The equation to be solved has the form:
!
!      dUdT - k * d2UdX2 = F(X,T)
!
!    over the interval [X_MIN,X_MAX] with boundary conditions
!
!      U(X_MIN,T) = U_X_MIN(T),
!      U(X_MIN,T) = U_X_MAX(T),
!
!    over the time interval [T_MIN,T_MAX] with initial conditions
!
!      U(X,T_MIN) = U_T_MIN(X)
!
!    The code uses the finite difference method to approximate the
!    second derivative in space, and an explicit forward Euler approximation
!    to the first derivative in time.
!
!    The finite difference form can be written as
!
!      U(X,T+dt) - U(X,T)                  ( U(X-dx,T) - 2 U(X,T) + U(X+dx,T) )
!      ------------------  = F(X,T) + k *  ------------------------------------
!               dt                                   dx * dx
!
!    or, assuming we have solved for all values of U at time T, we have
!
!      U(X,T+dt) = U(X,T) 
!        + cfl * ( U(X-dx,T) - 2 U(X,T) + U(X+dx,T) ) + dt * F(X,T) 
!
!    Here "cfl" is the Courant-Friedrichs-Loewy coefficient:
!
!      cfl = k * dt / dx / dx
!
!    In order for accurate results to be computed by this explicit method,
!    the CFL coefficient must be less than 0.5!
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 January 2022
!
!  Author:
! 
!    John Burkardt
!
!  Reference:
!
!    George Lindfield, John Penny,
!    Numerical Methods Using MATLAB,
!    Second Edition,
!    Prentice Hall, 1999,
!    ISBN: 0-13-012641-1,
!    LC: QA297.P45.
!
!  Input:
!
!    real ( kind = rk ) K, the heat conductivity coefficient.
!
!    integer T_NUM, the number of time values, including 
!    the initial value.
!
!    real ( kind = rk ) T_MIN, T_MAX, the minimum and maximum times.
!
!    integer X_NUM, the number of nodes.
!
!    real ( kind = rk ) X_MIN, X_MAX, the minimum and maximum spatial 
!    coordinates.
!
!  Output:
!
!    real ( kind = rk ) CFL, the Courant-Friedrichs-Loewy coefficient.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) cfl
  real ( kind = rk ) dx
  real ( kind = rk ) dt
  real ( kind = rk ) k
  real ( kind = rk ) t_max
  real ( kind = rk ) t_min
  integer t_num
  real ( kind = rk ) x_max
  real ( kind = rk ) x_min
  integer x_num

  dx = ( x_max - x_min ) / real ( x_num - 1, kind = rk )
  dt = ( t_max - t_min ) / real ( t_num - 1, kind = rk )
!
!  Check the CFL condition, print out its value, and quit if it is too large.
!
  cfl = k * dt / dx / dx

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  CFL stability criterion value = ', cfl

  if ( 0.5D+00 <= cfl ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FD1D_HEAT_EXPLICIT_CFL - Fatal error!'
    write ( *, '(a)' ) '  CFL condition failed.'
    write ( *, '(a)' ) '  0.5 <= K * dT / dX / dX = CFL.'
    stop 1
  end if

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! get_unit() returns a free FORTRAN unit number.
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
!    07 January 2022
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
!    07 January 2022
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
!    real ( kind = rk ) TABLE(M,N), the table data.
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
    stop 1
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
subroutine r8vec_linspace ( n, a_first, a_last, a )

!*****************************************************************************80
!
!! r8vec_linspace() creates a vector of linearly spaced values.
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
!    07 January 2022
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of entries in the vector.
!
!    real ( kind = rk ) A_FIRST, A_LAST, the first and last entries.
!
!  Output:
!
!    real ( kind = rk ) A(N), a vector of linearly spaced data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) a_first
  real ( kind = rk ) a_last
  integer i

  if ( n == 1 ) then

    a(1) = ( a_first + a_last ) / 2.0D+00

  else

    do i = 1, n
      a(i) = ( real ( n - i,     kind = rk ) * a_first &
             + real (     i - 1, kind = rk ) * a_last ) &
             / real ( n     - 1, kind = rk )
    end do

  end if

  return
end
subroutine r8vec_write ( output_filename, n, x )

!*****************************************************************************80
!
!! r8vec_write() writes an R8VEC file.
!
!  Discussion:
!
!    An R8VEC is a vector of R8 values.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 January 2022
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    character ( len = * ) OUTPUT_FILENAME, the output file name.
!
!    integer N, the number of points.
!
!    real ( kind = rk ) X(N), the data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer j
  character ( len = * ) output_filename
  integer output_status
  integer output_unit
  real ( kind = rk ) x(n)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_filename ) // '" on unit ', output_unit
    output_unit = -1
    stop 1
  end if

  if ( 0 < n ) then
!
!  Write the data.
!
    do j = 1, n
      write ( output_unit, '(2x,g24.16)' ) x(j)
    end do

  end if
!
!  Close the file.
!
  close ( unit = output_unit )

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
!    07 January 2022
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
 
