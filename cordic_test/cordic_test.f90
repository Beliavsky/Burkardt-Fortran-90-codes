program main

!*****************************************************************************80
!
!! cordic_test() tests cordic().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 June 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'cordic_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test cordic().'

  call arccos_cordic_test ( )
  call arcsin_cordic_test ( )
  call arctan_cordic_test ( )
  call cbrt_cordic_test ( )
  call cossin_cordic_test1 ( )
  call cossin_cordic_test2 ( )
  call exp_cordic_test ( )
  call ln_cordic_test ( )
  call multiply_cordic_test ( )
  call sqrt_cordic_test ( )
  call tan_cordic_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'cordic_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine arccos_cordic_test ( )

!*****************************************************************************80
!
!! arccos_cordic_test() tests arccos_cordic().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a1
  real ( kind = rk ) a2
  real ( kind = rk ) d
  integer n
  integer n_data
  real ( kind = rk ) t

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'arccos_cordic_test():'
  write ( *, '(a)' ) '  arccos_cordic() computes the arccosine of T'
  write ( *, '(a)' ) '  using the CORDIC algorithm.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '          T        N      ArcCos(T)       ArcCos(T)         Difference'
  write ( *, '(a)' ) &
    '                          Tabulated        CORDIC'

  n_data = 0

  do

    call arccos_values ( n_data, t, a1 )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(a)' ) ' '

    do n = 0, 25, 5

      call arccos_cordic ( t, n, a2 )

      d = a1 - a2

      write ( *, '(2x,f12.4,2x,i4,2x,g16.8,2x,g16.8,2x,e12.4)' ) &
        t, n, a1, a2, d

    end do

  end do

  return
end
subroutine arcsin_cordic_test ( )

!*****************************************************************************80
!
!! arcsin_cordic_test() tests arcsin_cordic().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a1
  real ( kind = rk ) a2
  real ( kind = rk ) d
  integer n
  integer n_data
  real ( kind = rk ) t

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'arcsin_cordic_test():'
  write ( *, '(a)' ) '  arcsin_cordic() computes the arcsine of T'
  write ( *, '(a)' ) '  using the CORDIC algorithm.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '          T        N      ArcSin(T)       ArcSin(T)         Difference'
  write ( *, '(a)' ) &
    '                          Tabulated        CORDIC'

  n_data = 0

  do

    call arcsin_values ( n_data, t, a1 )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(a)' ) ' '

    do n = 0, 25, 5

      call arcsin_cordic ( t, n, a2 )

      d = a1 - a2

      write ( *, '(2x,f12.4,2x,i4,2x,g16.8,2x,g16.8,2x,e12.4)' ) &
        t, n, a1, a2, d

    end do

  end do

  return
end
subroutine arctan_cordic_test ( )

!*****************************************************************************80
!
!! arctan_cordic_test() tests arctan_cordic().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a1
  real ( kind = rk ) a2
  real ( kind = rk ) d
  integer n
  integer n_data
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) x
  real ( kind = rk ) y
  real ( kind = rk ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'arctan_cordic_test():'
  write ( *, '(a)' ) '  arctan_cordic() computes the arctangent of Y/X'
  write ( *, '(a)' ) '  using the CORDIC algorithm.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '      X      Y    N       ArcTan(Y/X) ArcTan(Y/X)      Difference'
  write ( *, '(a)' ) &
    '                           Tabulated   CORDIC'
  n_data = 0

  do

    call arctan_values ( n_data, z, a1 )

    if ( n_data == 0 ) then
      exit
    end if

    call random_number ( harvest = r )

    x = r
    y = r * z

    call random_number ( harvest = s )

    if ( s < 0.5D+00 ) then
      x = -x
      y = -y
    end if

    write ( *, '(a)' ) ' '

    do n = 0, 25, 5

      call arctan_cordic ( x, y, n, a2 )

      d = a1 - a2

      write ( *, '(2x,f12.4,2x,f12.4,2x,i4,2x,g16.8,2x,g16.8,2x,e12.4)' ) &
        x, y, n, a1, a2, d

    end do

  end do

  return
end
subroutine cbrt_cordic_test ( )

!*****************************************************************************80
!
!! cbrt_cordic() tests cbrt_cordic().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) d
  real ( kind = rk ) fx1
  real ( kind = rk ) fx2
  integer n
  integer n_data
  real ( kind = rk ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'cbrt_cordic():'
  write ( *, '(a)' ) '  cbrt_cordic() computes the cube root function'
  write ( *, '(a)' ) '  using the CORDIC algorithm.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '          X          N      Cbrt(X)         Cbrt(X)           Difference'
  write ( *, '(a)' ) &
    '                          Tabulated        CORDIC'

  n_data = 0

  do

    call cbrt_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(a)' ) ' '

    do n = 0, 25, 5

      call cbrt_cordic ( x, n, fx2 )

      d = fx1 - fx2

      write ( *, '(2x,g14.6,2x,i4,2x,g16.8,2x,g16.8,2x,e12.4)' ) &
        x, n, fx1, fx2, d

    end do

  end do

  return
end
subroutine cossin_cordic_test1 ( )

!*****************************************************************************80
!
!! cossin_cordic_test1() tests cossin_cordic().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) c1
  real ( kind = rk ) c2
  real ( kind = rk ) d
  integer n
  integer n_data
  real ( kind = rk ) s2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'cossin_cordic_test1():'
  write ( *, '(a)' ) '  cossin_cordic() computes the cosine and sine'
  write ( *, '(a)' ) '  using the CORDIC algorithm.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '          A        N      Cos(A)           Cos(A)           Difference'
  write ( *, '(a)' ) &
    '                          Tabulated        CORDIC'

  n_data = 0

  do

    call cos_values ( n_data, a, c1 )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(a)' ) ' '

    do n = 0, 25, 5

      call cossin_cordic ( a, n, c2, s2 )

      d = c1 - c2

      write ( *, '(2x,f12.4,2x,i4,2x,g16.8,2x,g16.8,2x,e12.4)' ) &
        a, n, c1, c2, d

    end do

  end do

  return
end
subroutine cossin_cordic_test2 ( )

!*****************************************************************************80
!
!! cossin_cordic_test2() tests cossin_cordic().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) c2
  real ( kind = rk ) d
  integer n
  integer n_data
  real ( kind = rk ) s1
  real ( kind = rk ) s2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'cossin_cordic_test2():'
  write ( *, '(a)' ) '  cossin_cordic() computes the cosine and sine'
  write ( *, '(a)' ) '  using the CORDIC algorithm.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '          A        N      Sin(A)           Sin(A)           Difference'
  write ( *, '(a)' ) &
    '                          Tabulated        CORDIC'
  n_data = 0

  do

    call sin_values ( n_data, a, s1 )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(a)' ) ' '

    do n = 0, 25, 5

      call cossin_cordic ( a, n, c2, s2 )

      d = s1 - s2

      write ( *, '(2x,f12.4,2x,i4,2x,g16.8,2x,g16.8,2x,e12.4)' ) &
        a, n, s1, s2, d

    end do

  end do

  return
end
subroutine exp_cordic_test ( )

!*****************************************************************************80
!
!! exp_cordic_test() tests exp_cordic().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) d
  real ( kind = rk ) fx1
  real ( kind = rk ) fx2
  integer n
  integer n_data
  real ( kind = rk ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'exp_cordic_test():'
  write ( *, '(a)' ) '  exp_cordic() computes the exponential function'
  write ( *, '(a)' ) '  using the CORDIC algorithm.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '        X          N       Exp(X)          Exp(X)           Difference'
  write ( *, '(a)' ) &
    '                          Tabulated        CORDIC'

  n_data = 0

  do

    call exp_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(a)' ) ' '

    do n = 0, 25, 5

      call exp_cordic ( x, n, fx2 )

      d = fx1 - fx2

      write ( *, '(2x,f12.4,2x,i4,2x,g16.8,2x,g16.8,2x,e12.4)' ) &
        x, n, fx1, fx2, d

    end do

  end do

  return
end
subroutine ln_cordic_test ( )

!*****************************************************************************80
!
!! ln_cordic_test() tests ln_cordic().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) d
  real ( kind = rk ) fx1
  real ( kind = rk ) fx2
  integer n
  integer n_data
  real ( kind = rk ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ln_cordic_test():'
  write ( *, '(a)' ) '  ln_cordic() computes the natural logarithm function'
  write ( *, '(a)' ) '  using the CORDIC algorithm.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '          X          N        Ln(X)           Ln(X)           Difference'
  write ( *, '(a)' ) &
    '                          Tabulated        CORDIC'

  n_data = 0

  do

    call ln_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(a)' ) ' '

    do n = 0, 25, 5

      call ln_cordic ( x, n, fx2 )

      d = fx1 - fx2

      write ( *, '(2x,g14.6,2x,i4,2x,g16.8,2x,g16.8,2x,e12.4)' ) &
        x, n, fx1, fx2, d

    end do

  end do

  return
end
subroutine multiply_cordic_test ( )

!*****************************************************************************80
!
!! multiply_cordic_test() tests multiply_cordic().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 June 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  real ( kind = rk ) r
  real ( kind = rk ) x
  real ( kind = rk ) y
  real ( kind = rk ) z1
  real ( kind = rk ) z2

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'multiply_cordic_test():'
  write ( *, '(a)' ) '  multiply_cordic() computes Z = X * Y'
  write ( *, '(a)' ) '  using the CORDIC algorithm.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '        X             Y               Z                 Z'
  write ( *, '(a)' ) '                                      (X*Y)             (CORDIC)'
  write ( *, '(a)' ) ''

  do i = 1, 20

    call random_number ( harvest = r )
    x = - 100.0D+00 + 200.0D+00 * r
    call random_number ( harvest = r )
    y = - 100.0D+00 + 200.0D+00 * r
    z1 = x * y
    call multiply_cordic ( x, y, z2 )

    write ( *, '(2x,f12.8,2x,f12.8,2x,f16.8,2x,f16.8)' ) x, y, z1, z2

  end do

  return
end
subroutine sqrt_cordic_test ( )

!*****************************************************************************80
!
!! sqrt_cordic_test() tests sqrt_cordic().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) d
  real ( kind = rk ) fx1
  real ( kind = rk ) fx2
  integer n
  integer n_data
  real ( kind = rk ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'sqrt_cordic_test():'
  write ( *, '(a)' ) '  sqrt_cordic() computes the square root function'
  write ( *, '(a)' ) '  using the CORDIC algorithm.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '          X          N      Sqrt(X)         Sqrt(X)           Difference'
  write ( *, '(a)' ) &
    '                          Tabulated        CORDIC'

  n_data = 0

  do

    call sqrt_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(a)' ) ' '

    do n = 0, 25, 5

      call sqrt_cordic ( x, n, fx2 )

      d = fx1 - fx2

      write ( *, '(2x,g14.6,2x,i4,2x,g16.8,2x,g16.8,2x,e12.4)' ) &
        x, n, fx1, fx2, d

    end do

  end do

  return
end

subroutine tan_cordic_test ( )

!*****************************************************************************80
!
!! tan_cordic_test() tests tan_cordic().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) d
  integer n
  integer n_data
  real ( kind = rk ) t1
  real ( kind = rk ) t2
  real ( kind = rk ) theta

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'tan_cordic_test():'
  write ( *, '(a)' ) '  tan_cordic() computes the tangent of THETA'
  write ( *, '(a)' ) '  using the CORDIC algorithm.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '      THETA        N     Tan(THETA)      Tan(THETA)         Difference'
  write ( *, '(a)' ) &
    '                          Tabulated        CORDIC'

  n_data = 0

  do

    call tan_values ( n_data, theta, t1 )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(a)' ) ' '

    do n = 0, 25, 5

      call tan_cordic ( theta, n, t2 )

      d = t1 - t2

      write ( *, '(2x,f12.4,2x,i4,2x,g16.8,2x,g16.8,2x,e12.4)' ) &
        theta, n, t1, t2, d

    end do

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

