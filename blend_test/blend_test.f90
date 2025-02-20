program main

!*****************************************************************************80
!
!! blend_test() tests blend().
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
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'blend_test()'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test blend().'

  call blend_101_test ( )
  call blend_102_test ( )
  call blend_103_test ( )

  call blend_112_test ( )
  call blend_113_test ( )

  call blend_123_test ( )

  call blend_i_0d1_test ( )

  call blend_ij_0d1_test ( )
  call blend_ij_1d1_test ( )
  call blend_ij_w_1d1_test ( )

  call blend_ijk_0d1_test ( )
  call blend_ijk_1d1_test ( )
  call blend_ijk_2d1_test ( )

  call blend_r_0dn_test ( )
  call blend_r_0dn_identity_test ( )
  call blend_r_0dn_stretch_test ( )

  call blend_rs_0dn_test ( )
  call blend_rs_0dn_identity_test ( )
  call blend_rs_0dn_stretch_test ( )

  call blend_rs_1dn_test ( )
  call blend_rs_1dn_identity_test ( )
  call blend_rs_1dn_stretch_test ( )

  call blend_rst_1dn_identity_test ( )
  call blend_rst_1dn_stretch_test ( )

  call blend_rst_1dn_identity_test ( )
  call blend_rst_1dn_stretch_test ( )

  call blend_rst_2dn_identity_test ( )
  call blend_rst_2dn_stretch_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'blend_test()'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine blend_101_test ( )

!*****************************************************************************80
!
!! blend_101_test() tests blend_101().
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n1 = 10

  real ( kind = rk ) d(n1)
  integer i
  real ( kind = rk ) r

  d(1:n1) = 0.0D+00

  do i = 1, n1
    if ( i == 1 .or. i == n1 ) then
      d(i) = real ( i, kind = rk )
    end if
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'blend_101_test():'
  write ( *, '(a)' ) '  blend_101() blends endpoint values into a list.'

  call r8vec_print ( n1, d, '  Initial data list' )

  do i = 1, n1

    if ( i /= 1 .and. i /= n1 ) then
      r = real (  i - 1, kind = rk ) / real ( n1 - 1, kind = rk )
      call blend_101 ( r, d(1), d(n1), d(i) )
    end if

  end do

  call r8vec_print ( n1, d, '  Interpolated data list' )

  return
end
subroutine blend_102_test ( )

!*****************************************************************************80
!
!! blend_102_test() tests blend_102().
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n1 = 5
  integer, parameter :: n2 = 5

  real ( kind = rk ) d(n1,n2)
  integer i
  integer j
  real ( kind = rk ) r
  real ( kind = rk ) s

  d(1:n1,1:n2) = 0.0D+00

  do i = 1, n1
    do j = 1, n2
      if ( &
        ( i == 1 .or. i == n1 ) .and. &
        ( j == 1 .or. j == n2 ) ) then
        d(i,j) = real ( i + j, kind = rk )
      end if
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'blend_102_test():'
  write ( *, '(a)' ) '  blend_102() blends corner values into a table.'

  call r8mat_print ( n1, n2, d, '  Initial data array' )

  do i = 1, n1

    r = real (  i - 1, kind = rk ) &
      / real ( n1 - 1, kind = rk )

    do j = 1, n2

      s = real (  j - 1, kind = rk ) &
        / real ( n2 - 1, kind = rk )

      if ( &
        ( i == 1 .or. i == n1 ) .and. &
        ( j == 1 .or. j == n2 ) ) then
        cycle
      end if

      call blend_102 ( r, s, d(1,1), d(1,n2), d(n1,1), d(n1,n2), d(i,j) )

    end do

  end do

  call r8mat_print ( n1, n2, d, '  Interpolated data array' )

  return
end
subroutine blend_103_test ( )

!*****************************************************************************80
!
!! blend_103_test() tests blend_103().
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n1 = 3
  integer, parameter :: n2 = 5
  integer, parameter :: n3 = 4

  real ( kind = rk ) d(n1,n2,n3)
  integer i
  integer j
  integer k
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) t

  d(1:n1,1:n2,1:n3) = 0.0D+00

  do i = 1, n1
    do j = 1, n2
      do k = 1, n3
        if ( &
          ( i == 1 .or. i == n1 ) .and. &
          ( j == 1 .or. j == n2 ) .and. &
          ( k == 1 .or. k == n3  ) ) then
          d(i,j,k) = real ( i + j + k, kind = rk )
        end if
      end do
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'blend_103_test():'
  write ( *, '(a)' ) '  blend_103() blends corner values into a table.'

  call r8block_print ( n1, n2, n3, d, '  Initial data array' )

  do i = 1, n1

    r = real (  i - 1, kind = rk ) &
      / real ( n1 - 1, kind = rk )

    do j = 1, n2

      s = real (  j - 1, kind = rk ) &
        / real ( n2 - 1, kind = rk )

      do k = 1, n3

        t = real (  k - 1, kind = rk ) &
          / real ( n3 - 1, kind = rk )

        if ( &
          ( i == 1 .or. i == n1 ) .and. &
          ( j == 1 .or. j == n2 ) .and. &
          ( k == 1 .or. k == n3  ) ) then
          cycle
        end if

        call blend_103 ( r, s, t, &
          d(1,1,1), d(1,1,n3), d(1,n2,1), d(1,n2,n3), &
          d(n1,1,1), d(n1,1,n3), d(n1,n2,1), d(n1,n2,n3), d(i,j,k) )

      end do

    end do

  end do

  call r8block_print ( n1, n2, n3, d, '  Interpolated data array' )

  return
end
subroutine blend_112_test ( )

!*****************************************************************************80
!
!! blend_112_test() tests blend_112().
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n1 = 5
  integer, parameter :: n2 = 5

  real ( kind = rk ) d(n1,n2)
  integer i
  integer j
  real ( kind = rk ) r
  real ( kind = rk ) s

  d(1:n1,1:n2) = 0.0D+00

  do i = 1, n1
    do j = 1, n2
      if ( i == 1 .or. i == n1 .or. j == 1 .or. j == n2 ) then
        d(i,j) = real ( i + j, kind = rk )
      end if
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'blend_112_test():'
  write ( *, '(a)' ) '  blend_112() blends side values into a table.'

  call r8mat_print ( n1, n2, d, '  Initial data array' )

  do i = 1, n1

    r = real (  i - 1, kind = rk ) &
      / real ( n1 - 1, kind = rk )

    do j = 1, n2

      s = real (  j - 1, kind = rk ) &
        / real ( n2 - 1, kind = rk )

      if ( i == 1 .or. i == n1 .or. j == 1 .or. j == n2 ) then
        cycle
      end if

      call blend_112 ( r, s, d(1,1), d(1,n2), d(n1,1), d(n1,n2), &
        d(i,1), d(i,n2), d(1,j), d(n1,j), d(i,j) )

    end do

  end do

  call r8mat_print ( n1, n2, d, '  Interpolated data array' )

  return
end
subroutine blend_113_test ( )

!*****************************************************************************80
!
!! blend_113_test() tests blend_113().
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n1 = 3
  integer, parameter :: n2 = 5
  integer, parameter :: n3 = 4

  real ( kind = rk ) d(n1,n2,n3)
  integer i
  integer j
  integer k
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) t

  d(1:n1,1:n2,1:n3) = 0.0D+00

  do i = 1, n1
    do j = 1, n2
      do k = 1, n3
        if ( &
          ( ( i == 1 .or. i == n1 ) .and. ( j == 1 .or. j == n2 ) ) .or. &
          ( ( i == 1 .or. i == n1 ) .and. ( k == 1 .or. k == n3 ) ) .or. &
          ( ( j == 1 .or. j == n2 ) .and. ( k == 1 .or. k == n3 ) ) ) then
          d(i,j,k) = real ( i + j + k, kind = rk )
        end if
      end do
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'blend_113_test():'
  write ( *, '(a)' ) '  blend_113() blends edge values into a table.'

  call r8block_print ( n1, n2, n3, d, '  Initial data array' )

  do i = 1, n1

    r = real (  i - 1, kind = rk ) &
      / real ( n1 - 1, kind = rk )

    do j = 1, n2

      s = real (  j - 1, kind = rk ) &
        / real ( n2 - 1, kind = rk )

      do k = 1, n3

        t = real (  k - 1, kind = rk ) &
          / real ( n3 - 1, kind = rk )

        if ( &
          ( ( i == 1 .or. i == n1 ) .and. ( j == 1 .or. j == n2 ) ) .or. &
          ( ( i == 1 .or. i == n1 ) .and. ( k == 1 .or. k == n3 ) ) .or. &
          ( ( j == 1 .or. j == n2 ) .and. ( k == 1 .or. k == n3 ) ) ) then
          cycle
        end if

        call blend_113 ( r, s, t, &
          d(1,1,1), d(1,1,n3), d(1,n2,1), d(1,n2,n3), &
          d(n1,1,1), d(n1,1,n3), d(n1,n2,1), d(n1,n2,n3), &
          d(i,1,1), d(i,1,n3), d(i,n2,1), d(i,n2,n3), &
          d(1,j,1), d(1,j,n3), d(n1,j,1), d(n1,j,n3), &
          d(1,1,k), d(1,n2,k), d(n1,1,k), d(n1,n2,k), &
          d(i,j,k) )

      end do

    end do

  end do

  call r8block_print ( n1, n2, n3, d, '  Interpolated data array' )

  return
end
subroutine blend_123_test ( )

!*****************************************************************************80
!
!! blend_123_test() tests blend_123().
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n1 = 3
  integer, parameter :: n2 = 5
  integer, parameter :: n3 = 4

  real ( kind = rk ) d(n1,n2,n3)
  integer i
  integer j
  integer k
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) t

  d(1:n1,1:n2,1:n3) = 0.0D+00

  do i = 1, n1
    do j = 1, n2
      do k = 1, n3
        if ( &
          ( i == 1 .or. i == n1 ) .or. &
          ( j == 1 .or. j == n2 ) .or. &
          ( k == 1 .or. k == n3  ) ) then
          d(i,j,k) = real ( i + j + k, kind = rk )
        end if
      end do
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'blend_123_test():'
  write ( *, '(a)' ) '  blend_123() blends face values into a table.'

  call r8block_print ( n1, n2, n3, d, '  Initial data array' )

  do i = 1, n1

    r = real (  i - 1, kind = rk ) &
      / real ( n1 - 1, kind = rk )

    do j = 1, n2

      s = real (  j - 1, kind = rk ) &
        / real ( n2 - 1, kind = rk )

      do k = 1, n3

        t = real (  k - 1, kind = rk ) &
          / real ( n3 - 1, kind = rk )

        if ( &
          ( i == 1 .or. i == n1 ) .or. &
          ( j == 1 .or. j == n2 ) .or. &
          ( k == 1 .or. k == n3  ) ) then
          cycle
        end if

        call blend_123 ( r, s, t, &
          d(1,1,1), d(1,1,n3), d(1,n2,1), d(1,n2,n3), &
          d(n1,1,1), d(n1,1,n3), d(n1,n2,1), d(n1,n2,n3), &
          d(i,1,1), d(i,1,n3), d(i,n2,1), d(i,n2,n3), &
          d(1,j,1), d(1,j,n3), d(n1,j,1), d(n1,j,n3), &
          d(1,1,k), d(1,n2,k), d(n1,1,k), d(n1,n2,k), &
          d(1,j,k), d(n1,j,k), d(i,1,k), d(i,n2,k), d(i,j,1), d(i,j,n3), &
          d(i,j,k) )

      end do

    end do

  end do

  call r8block_print ( n1, n2, n3, d, '  Interpolated data array' )

  return
end
subroutine blend_i_0d1_test ( )

!*****************************************************************************80
!
!! blend_i_0d1_test() tests blend_i_0d1().
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5

  integer i
  real ( kind = rk ) x(m)

  x(1) = 100.0D+00
  x(m) = 100.0 + real ( ( m - 1 ) * 5, kind = rk )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'blend_i_0d1_test():'
  write ( *, '(a)' ) '  blend_i_0d1() interpolates data in a vector.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  X(1) = ', x(1)
  write ( *, '(a,i2,a,g14.6)' ) '  X(', m, ')= ', x(m)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Interpolated values:'
  write ( *, '(a)' ) ' '

  call blend_i_0d1 ( x, m )

  do i = 1, m
    write ( *, '(2x,i8,2x,g14.6)' ) i, x(i)
  end do

  return
end
subroutine blend_ij_0d1_test ( )

!*****************************************************************************80
!
!! blend_ij_0d1_test() tests blend_ij_0d1().
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m1 = 5
  integer, parameter :: m2 = 4

  external cubic_rs
  integer i
  integer j
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) x(m1,m2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'blend_ij_0d1_test():'
  write ( *, '(a)' ) '  blend_ij_0d1() interpolates data in a table,'
  write ( *, '(a)' ) '  from corner data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a,i8,a)' ) '  The table is ', m1, ' rows by ', &
    m2, ' columns.'
!
!  Load data in the corners only.
!
  i = 1
  j = 1
  r = real (  i - 1, kind = rk ) &
    / real ( m1 - 1, kind = rk )
  s = real (  j - 1, kind = rk ) &
    / real ( m2 - 1, kind = rk )
  call cubic_rs ( r, s, 1, x(i,j) )

  i = m1
  j = 1
  r = real (  i - 1, kind = rk ) &
    / real ( m1 - 1, kind = rk )
  s = real (  j - 1, kind = rk ) &
    / real ( m2 - 1, kind = rk )
  call cubic_rs ( r, s, 1, x(i,j) )

  i = 1
  j = m2
  r = real (  i - 1, kind = rk ) &
    / real ( m1 - 1, kind = rk )
  s = real (  j - 1, kind = rk ) &
    / real ( m2 - 1, kind = rk )
  call cubic_rs ( r, s, 1, x(i,j) )

  i = m1
  j = m2
  r = real (  i - 1, kind = rk ) &
    / real ( m1 - 1, kind = rk )
  s = real (  j - 1, kind = rk ) &
    / real ( m2 - 1, kind = rk )
  call cubic_rs ( r, s, 1, x(i,j) )

  call blend_ij_0d1 ( x, m1, m2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Values interpolated by blend_ij_0d1:'
  write ( *, '(a)' ) ' '

  do i = 1, m1
    write ( *, '(5g14.6)' ) x(i,1:m2)
  end do

  return
end
subroutine blend_ij_1d1_test ( )

!*****************************************************************************80
!
!! blend_ij_1d1_test() tests blend_ij_1d1().
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m1 = 5
  integer, parameter :: m2 = 4

  external cubic_rs
  integer i
  integer j
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) x(m1,m2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'blend_ij_1d1_test():'
  write ( *, '(a)' ) '  blend_ij_1d1() interpolates data in a table,'
  write ( *, '(a)' ) '  from edge data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a,i8,a)' ) '  The table is ', m1, ' rows by ', &
    m2, ' columns.'
!
!  Load data in the edges.
!
  j = 1
  do i = 1, m1
    r = real (  i - 1, kind = rk ) &
      / real ( m1 - 1, kind = rk )
    s = real (  j - 1, kind = rk ) &
      / real ( m2 - 1, kind = rk )
    call cubic_rs ( r, s, 1, x(i,j) )
  end do

  j = m2
  do i = 1, m1
    r = real (  i - 1, kind = rk ) &
      / real ( m1 - 1, kind = rk )
    s = real (  j - 1, kind = rk ) &
      / real ( m2 - 1, kind = rk )
    call cubic_rs ( r, s, 1, x(i,j) )
  end do

  i = 1
  do j = 2, m2 - 1
    r = real (  i - 1, kind = rk ) &
      / real ( m1 - 1, kind = rk )
    s = real (  j - 1, kind = rk ) &
      / real ( m2 - 1, kind = rk )
    call cubic_rs ( r, s, 1, x(i,j) )
  end do

  i = m1
  do j = 2, m2 - 1
    r = real (  i - 1, kind = rk ) &
      / real ( m1 - 1, kind = rk )
    s = real (  j - 1, kind = rk ) &
      / real ( m2 - 1, kind = rk )
    call cubic_rs ( r, s, 1, x(i,j) )
  end do

  call blend_ij_1d1 ( x, m1, m2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Values interpolated by blend_ij_1d1:'
  write ( *, '(a)' ) ' '

  do i = 1, m1
    write ( *, '(5g14.6)' ) x(i,1:m2)
  end do

  return
end
subroutine blend_ij_w_1d1_test ( )

!*****************************************************************************80
!
!! blend_ij_w_1d1_test() tests blend_ij_w_1d1().
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n1 = 5
  integer, parameter :: n2 = 5

  integer i
  integer j
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) r(n1)
  real ( kind = rk ) rad
  real ( kind = rk ) rr
  real ( kind = rk ) s(n2)
  real ( kind = rk ) ss
  real ( kind = rk ) x(n1,n2)
  real ( kind = rk ) y(n1,n2)

  rad = 3.0D+00

  x(1:n1,1:n2) = 0.0D+00
  y(1:n1,1:n2) = 0.0D+00
!
!  Set the boundary values.
!
!  It turns out that our values correspond to the X and Y
!  coordinates of a quarter circle of radius 3, although
!  it is by no means necessary that a formula for the data
!  be known.
!
  do i = 1, n1
    rr = ( real ( i - 1, kind = rk ) / real ( n1 - 1, kind = rk ) )**2
    r(i) = rr
    x(i,1) = 0.0D+00
    y(i,1) = 0.0D+00
    x(i,n2) = rad * cos ( 0.5D+00 * pi * ( 1.0D+00 - rr ) )
    y(i,n2) = rad * sin ( 0.5D+00 * pi * ( 1.0D+00 - rr ) )
  end do

  do j = 1, n2
    ss = ( real ( j - 1, kind = rk ) / real ( n2 - 1, kind = rk ) )**2
    s(j) = ss
    x(1,j) = 0.0D+00
    y(1,j) = rad * ss
    x(n1,j) = rad * ss
    y(n1,j) = 0.0D+00
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'blend_ij_w_1d1_test():'
  write ( *, '(a)' ) '  blend_ij_w_1d1() uses weighted blending to fill in the'
  write ( *, '(a)' ) '  interior of a table.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     R           S           X           Y'
  write ( *, '(a)' ) ' '
  call blend_ij_w_1d1 ( x, r, s, n1, n2 )
  call blend_ij_w_1d1 ( y, r, s, n1, n2 )

  do i = 1, n1
    write ( *, '(a)' ) ' '
    do j = 1, n2
      write ( *, '(4g12.4)' ) r(i), s(j), x(i,j), y(i,j)
    end do
  end do

  return
end
subroutine blend_ijk_0d1_test ( )

!*****************************************************************************80
!
!! blend_ijk_0d1_test() tests blend_ijk_0d1().
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m1 = 4
  integer, parameter :: m2 = 3
  integer, parameter :: m3 = 3

  integer i
  integer j
  integer k
  integer num_extreme
  external quad_rst
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) t
  real ( kind = rk ) x(m1,m2,m3)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'blend_ijk_0d1_test():'
  write ( *, '(a)' ) '  blend_ijk_0d1() interpolates data in a table,'
  write ( *, '(a)' ) '  from corner data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a,i8,a,i8,a)' ) '  The table is ', m1, ' rows by ', &
    m2, ' columns by ', m3, ' layers.'
!
!  Load data on the faces.
!
  do i = 1, m1
    r = real (  i - 1, kind = rk ) &
      / real ( m1 - 1, kind = rk )
    do j = 1, m2
      s = real (  j - 1, kind = rk ) &
        / real ( m2 - 1, kind = rk )
      do k = 1, m3
        t = real (  k - 1, kind = rk ) &
          / real ( m3 - 1, kind = rk )

        num_extreme = 0
        if ( i == 1 .or. i == m1 ) then
          num_extreme = num_extreme + 1
        end if
        if ( j == 1 .or. j == m2 ) then
          num_extreme = num_extreme + 1
        end if
        if ( k == 1 .or. k == m3 ) then
          num_extreme = num_extreme + 1
        end if

        if ( num_extreme == 3 ) then
          call quad_rst ( r, s, t, 1, x(i,j,k) )
        else
          x(i,j,k) = 0.0D+00
        end if

      end do
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data given to blend_ijk_0d1:'
  write ( *, '(a)' ) ' '

  do k = 1, m3
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Layer K = ', k
    write ( *, '(a)' ) ' '
    do i = 1, m1
      write ( *, '(5g14.6)' ) x(i,1:m2,k)
    end do
  end do

  call blend_ijk_0d1 ( x, m1, m2, m3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Values interpolated by blend_ijk_0d1:'
  write ( *, '(a)' ) ' '

  do k = 1, m3
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Layer K = ', k
    write ( *, '(a)' ) ' '
    do i = 1, m1
      write ( *, '(5g14.6)' ) x(i,1:m2,k)
    end do
  end do
!
!  Load all data.
!
  do i = 1, m1
    r = real (  i - 1, kind = rk ) &
      / real ( m1 - 1, kind = rk )
    do j = 1, m2
      s = real (  j - 1, kind = rk ) &
        / real ( m2 - 1, kind = rk )
      do k = 1, m3
        t = real (  k - 1, kind = rk ) &
          / real ( m3 - 1, kind = rk )
        call quad_rst ( r, s, t, 1, x(i,j,k) )
      end do
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Exact data:'

  do k = 1, m3
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Layer K = ', k
    write ( *, '(a)' ) ' '
    do i = 1, m1
      write ( *, '(5g14.6)' ) x(i,1:m2,k)
    end do
  end do

  return
end
subroutine blend_ijk_1d1_test ( )

!*****************************************************************************80
!
!! blend_ijk_1d1_test() tests blend_ijk_1d1().
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m1 = 4
  integer, parameter :: m2 = 3
  integer, parameter :: m3 = 3

  integer i
  integer j
  integer k
  integer num_extreme
  external quad_rst
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) t
  real ( kind = rk ) x(m1,m2,m3)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'blend_ijk_1d1_test():'
  write ( *, '(a)' ) '  blend_ijk_1d1() interpolates data in a table,'
  write ( *, '(a)' ) '  from edge data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a,i8,a,i8,a)' ) '  The table is ', m1, ' rows by ', &
    m2, ' columns by ', m3, ' layers.'
!
!  Load data on the faces.
!
  do i = 1, m1
    r = real (  i - 1, kind = rk ) &
      / real ( m1 - 1, kind = rk )
    do j = 1, m2
      s = real (  j - 1, kind = rk ) &
        / real ( m2 - 1, kind = rk )
      do k = 1, m3
        t = real (  k - 1, kind = rk ) &
          / real ( m3 - 1, kind = rk )

        num_extreme = 0
        if ( i == 1 .or. i == m1 ) then
          num_extreme = num_extreme + 1
        end if
        if ( j == 1 .or. j == m2 ) then
          num_extreme = num_extreme + 1
        end if
        if ( k == 1 .or. k == m3 ) then
          num_extreme = num_extreme + 1
        end if

        if ( 2 <= num_extreme ) then
          call quad_rst ( r, s, t, 1, x(i,j,k) )
        else
          x(i,j,k) = 0.0D+00
        end if

      end do
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data given to blend_ijk_1d1:'
  write ( *, '(a)' ) ' '

  do k = 1, m3
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Layer K = ', k
    write ( *, '(a)' ) ' '
    do i = 1, m1
      write ( *, '(5g14.6)' ) x(i,1:m2,k)
    end do
  end do

  call blend_ijk_1d1 ( x, m1, m2, m3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Values interpolated by blend_ijk_1d1:'
 
  do k = 1, m3
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Layer K = ', k
    write ( *, '(a)' ) ' '
    do i = 1, m1
      write ( *, '(5g14.6)' ) x(i,1:m2,k)
    end do
  end do
!
!  Load all data.
!
  do i = 1, m1
    r = real (  i - 1, kind = rk ) &
      / real ( m1 - 1, kind = rk )
    do j = 1, m2
      s = real (  j - 1, kind = rk ) &
        / real ( m2 - 1, kind = rk )
      do k = 1, m3
        t = real (  k - 1, kind = rk ) &
          / real ( m3 - 1, kind = rk )
        call quad_rst ( r, s, t, 1, x(i,j,k) )
      end do
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Exact data:'

  do k = 1, m3
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Layer K = ', k
    write ( *, '(a)' ) ' '
    do i = 1, m1
      write ( *, '(5g14.6)' ) x(i,1:m2,k)
    end do
  end do

  return
end
subroutine blend_ijk_2d1_test ( )

!*****************************************************************************80
!
!! blend_ijk_2d1_test() tests blend_ijk_2d1().
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m1 = 4
  integer, parameter :: m2 = 3
  integer, parameter :: m3 = 3

  integer i
  integer j
  integer k
  integer num_extreme
  external quad_rst
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) t
  real ( kind = rk ) x(m1,m2,m3)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'blend_ijk_2d1_test():'
  write ( *, '(a)' ) '  blend_ijk_2d1() interpolates data in a table,'
  write ( *, '(a)' ) '  from face data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a,i8,a,i8,a)' ) '  The table is ', m1, ' rows by ', &
    m2, ' columns by ', &
    m3, ' layers.'
!
!  Load data on the faces.
!
  do i = 1, m1
    r = real (  i - 1, kind = rk ) &
      / real ( m1 - 1, kind = rk )
    do j = 1, m2
      s = real (  j - 1, kind = rk ) &
        / real ( m2 - 1, kind = rk )
      do k = 1, m3
        t = real (  k - 1, kind = rk ) &
          / real ( m3 - 1, kind = rk )

        num_extreme = 0
        if ( i == 1 .or. i == m1 ) then
          num_extreme = num_extreme + 1
        end if
        if ( j == 1 .or. j == m2 ) then
          num_extreme = num_extreme + 1
        end if
        if ( k == 1 .or. k == m3 ) then
          num_extreme = num_extreme + 1
        end if

        if ( 1 <= num_extreme ) then
          call quad_rst ( r, s, t, 1, x(i,j,k) )
        else
          x(i,j,k) = 0.0D+00
        end if

      end do
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data given to blend_ijk_2d1:'
  write ( *, '(a)' ) ' '

  do k = 1, m3
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Layer K = ', k
    write ( *, '(a)' ) ' '
    do i = 1, m1
      write ( *, '(5g14.6)' ) x(i,1:m2,k)
    end do
  end do

  call blend_ijk_2d1 ( x, m1, m2, m3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Values interpolated by blend_ijk_2d1:'
 
  do k = 1, m3
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Layer K = ', k
    write ( *, '(a)' ) ' '
    do i = 1, m1
      write ( *, '(5g14.6)' ) x(i,1:m2,k)
    end do
  end do
!
!  Load all data.
!
  do i = 1, m1
    r = real (  i - 1, kind = rk ) &
      / real ( m1 - 1, kind = rk )
    do j = 1, m2
      s = real (  j - 1, kind = rk ) &
        / real ( m2 - 1, kind = rk )
      do k = 1, m3
        t = real (  k - 1, kind = rk ) &
          / real ( m3 - 1, kind = rk )
        call quad_rst ( r, s, t, 1, x(i,j,k) )
      end do
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Exact data:'

  do k = 1, m3
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Layer K = ', k
    write ( *, '(a)' ) ' '
    do i = 1, m1
      write ( *, '(5g14.6)' ) x(i,1:m2,k)
    end do
  end do

  return
end
subroutine blend_r_0dn_test ( )

!*****************************************************************************80
!
!! blend_r_0dn_test() checks out blend_r_0dn().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 October 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 11
  integer, parameter :: n = 2
 
  integer i
  integer j
  external powers_r
  real ( kind = rk ) r
  real ( kind = rk ) x(m,n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'blend_r_0dn_test():'
  write ( *, '(a)' ) '  blend_r_0dn() interpolates endpoint vector data'
  write ( *, '(a)' ) '  into a list.'
  write ( *, '(a)' ) ''
  write ( *, '(a,i3,a)' ) '  The list is ', m, ' rows long.'
  write ( *, '(a,i3)' ) '  Each vector has length', n

  x(1:m,1:n) = 0.0D+00

  do i = 1, m
    r = real ( i - 1, kind = rk ) / real ( m - 1, kind = rk )
    call blend_r_0dn ( r, x(i,1:n), n, powers_r )
  end do

  call r8mat_print ( m, n, x, '  Data blended by blend_r_0dn:' );
!
!  Load all data.
!
  do i = 1, m
    r = real ( i - 1, kind = rk ) / real ( m - 1, kind = rk )
    do j = 1, n
      call powers_r ( r, j, x(i,j) )
    end do
  end do

  call r8mat_print ( m, n, x, '  Exact data:' )

  return
end
subroutine blend_r_0dn_identity_test ( )

!*****************************************************************************80
!
!! blend_r_0dn_identity_test() checks for a gross error in the blend coefficients.
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 1

  external identity_r
  real ( kind = rk ) r
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'blend_r_0dn_identity_test():'
  write ( *, '(a)' ) '  Simple identity test to detect gross errors.'
!
!  Test blend_r_0dn on identity.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Identity test for blend_r_0dn:'
  write ( *, '(a)' ) ' '
  r = 0.0D+00
  call blend_r_0dn ( r, x, n, identity_r )
  write ( *, '(2f8.4)' ) r, x(1:n)

  r = 1.0D+00
  call blend_r_0dn ( r, x, n, identity_r )
  write ( *, '(2f8.4)' ) r, x(1:n)

  r = 0.5D+00
  call blend_r_0dn ( r, x, n, identity_r )
  write ( *, '(2f8.4)' ) r, x(1:n)

  return
end
subroutine blend_r_0dn_stretch_test ( )

!*****************************************************************************80
!
!! blend_r_0dn_stretch_test() checks for simple errors in the blend coefficients.
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 1

  real ( kind = rk ) r
  external stretch_r
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'blend_r_0dn_stretch_test():'
  write ( *, '(a)' ) '  Shift and stretch test to detect simple errors.'
!
!  Test blend_r_0dn on shift by 1, stretch by 2.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Shift and stretch test for blend_r_0dn:'
  write ( *, '(a)' ) ' '
  r = 0.0D+00
  call blend_r_0dn ( r, x, n, stretch_r )
  write ( *, '(2f8.4)' ) r, x(1:n)

  r = 1.0D+00
  call blend_r_0dn ( r, x, n, stretch_r )
  write ( *, '(2f8.4)' ) r, x(1:n)

  r = 0.5D+00
  call blend_r_0dn ( r, x, n, stretch_r )
  write ( *, '(2f8.4)' ) r, x(1:n)

  return
end
subroutine blend_rs_0dn_test ( )

!*****************************************************************************80
!
!! blend_rs_0dn_test() tests blend_rs_0dn().
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m1 = 5
  integer, parameter :: m2 = 4

  external cubic_rs
  integer i
  integer j
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) x(m1,m2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'blend_rs_0dn_test():'
  write ( *, '(a)' ) '  blend_rs_0dn() interpolates data in a table,'
  write ( *, '(a)' ) '  from corner data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a,i8,a)' ) '  The table is ', m1, ' rows by ', &
    m2, ' columns.'

  do i = 1, m1
    r = real (  i - 1, kind = rk ) &
      / real ( m1 - 1, kind = rk )
    do j = 1, m2
      s = real (  j - 1, kind = rk ) &
        / real ( m2 - 1, kind = rk )
      call blend_rs_0dn ( r, s, x(i,j), 1, cubic_rs )
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data blended by blend_rs_0dn:'
  write ( *, '(a)' ) ' '

  do i = 1, m1
    write ( *, '(5g14.6)' ) x(i,1:m2)
  end do
!
!  Load all data.
!
  do i = 1, m1
    do j = 1, m2
      r = real (  i - 1, kind = rk ) &
        / real ( m1 - 1, kind = rk )
      s = real (  j - 1, kind = rk ) &
        / real ( m2 - 1, kind = rk )
      call cubic_rs ( r, s, 1, x(i,j) )
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Exact data:'
  write ( *, '(a)' ) ' '

  do i = 1, m1
    write ( *, '(5g14.6)' ) x(i,1:m2)
  end do

  return
end
subroutine blend_rs_0dn_identity_test ( )

!*****************************************************************************80
!
!! blend_rs_0dn_identity_test() checks for a gross error in the blend coefficients.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 October 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 2

  external identity_rs
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'blend_rs_0dn_identity_test():'
  write ( *, '(a)' ) '  Simple identity test to detect gross errors.'
!
!  Test blend_rs_0dn on identity.
!
  write ( *, '(a)' ) ' '

  write ( *, '(a)' ) '  Identity test for blend_rs_0dn:'
  write ( *, '(a)' ) ' '
  r = 0.0D+00
  s = 0.0D+00
  call blend_rs_0dn ( r, s, x, n, identity_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 1.0D+00
  s = 0.0D+00
  call blend_rs_0dn ( r, s, x, n, identity_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 0.0D+00
  s = 1.0D+00
  call blend_rs_0dn ( r, s, x, n, identity_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 1.0D+00
  s = 1.0D+00
  call blend_rs_0dn ( r, s, x, n, identity_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 0.5D+00
  s = 0.5D+00
  call blend_rs_0dn ( r, s, x, n, identity_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  return
end
subroutine blend_rs_0dn_stretch_test ( )

!*****************************************************************************80
!
!! blend_rs_0dn_stretch_test() checks for simple errors in the blend coefficients.
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 2

  real ( kind = rk ) r
  real ( kind = rk ) s
  external stretch_rs
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'blend_rs_0dn_stretch_test():'
  write ( *, '(a)' ) '  Shift and stretch test to detect simple errors.'
!
!  Test blend_rs_0dn on shift by (1,2), stretch by (3,4).
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Shift and stretch test for blend_rs_0dn:'
  write ( *, '(a)' ) ' '
  r = 0.0D+00
  s = 0.0D+00
  call blend_rs_0dn ( r, s, x, n, stretch_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 1.0D+00
  s = 0.0D+00
  call blend_rs_0dn ( r, s, x, n, stretch_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 0.0D+00
  s = 1.0D+00
  call blend_rs_0dn ( r, s, x, n, stretch_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 1.0D+00
  s = 1.0D+00
  call blend_rs_0dn ( r, s, x, n, stretch_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 0.5D+00
  s = 0.5D+00
  call blend_rs_0dn ( r, s, x, n, stretch_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  return
end
subroutine blend_rs_1dn_test ( )

!*****************************************************************************80
!
!! blend_rs_1dn_test() tests blend_rs_1dn().
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m1 = 5
  integer, parameter :: m2 = 4

  external cubic_rs
  integer i
  integer j
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) x(m1,m2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'blend_rs_1dn_test():'
  write ( *, '(a)' ) '  blend_rs_1dn() interpolates data in a table,'
  write ( *, '(a)' ) '  from edge data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a,i8,a)' ) '  The table is ', m1, ' rows by ', &
    m2, ' columns.'

  do i = 1, m1
    r = real (  i - 1, kind = rk ) &
      / real ( m1 - 1, kind = rk )
    do j = 1, m2
      s = real (  j - 1, kind = rk ) &
        / real ( m2 - 1, kind = rk )
      call blend_rs_1dn ( r, s, x(i,j), 1, cubic_rs )
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data blended by blend_rs_1dn:'
  write ( *, '(a)' ) ' '

  do i = 1, m1
    write ( *, '(5g14.6)' ) x(i,1:m2)
  end do
!
!  Load all data.
!
  do i = 1, m1
    do j = 1, m2
      r = real (  i - 1, kind = rk ) &
        / real ( m1 - 1, kind = rk )
      s = real (  j - 1, kind = rk ) &
        / real ( m2 - 1, kind = rk )
      call cubic_rs ( r, s, 1, x(i,j) )
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Exact data:'
  write ( *, '(a)' ) ' '

  do i = 1, m1
    write ( *, '(5g14.6)' ) x(i,1:m2)
  end do

  return
end
subroutine blend_rs_1dn_identity_test ( )

!*****************************************************************************80
!
!! blend_rs_1dn_identity_test() checks for a gross error in the blend coefficients.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 October 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 2

  external identity_rs
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'blend_rs_1dn_identity_test():'
  write ( *, '(a)' ) '  Simple identity test to detect gross errors.'
!
!  Test blend_rs_1dn on identity.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Identity test for blend_rs_1dn:'
  write ( *, '(a)' ) ' '
  r = 0.0D+00
  s = 0.0D+00
  call blend_rs_1dn ( r, s, x, n, identity_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 1.0D+00
  s = 0.0D+00
  call blend_rs_1dn ( r, s, x, n, identity_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 0.0D+00
  s = 1.0D+00
  call blend_rs_1dn ( r, s, x, n, identity_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 1.0D+00
  s = 1.0D+00
  call blend_rs_1dn ( r, s, x, n, identity_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 0.5D+00
  s = 0.5D+00
  call blend_rs_1dn ( r, s, x, n, identity_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  return
end
subroutine blend_rs_1dn_stretch_test ( )

!*****************************************************************************80
!
!! blend_rs_1dn_stretch_test() checks for simple errors in the blend coefficients.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 October 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 2

  real ( kind = rk ) r
  real ( kind = rk ) s
  external stretch_rs
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'blend_rs_1dn_stretch_test():'
  write ( *, '(a)' ) '  Shift and stretch test to detect simple errors.'
!
!  Test blend_rs_1D on shift by (1,2), stretch by (3,4).
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Shift and stretch test for blend_rs_1dn:'
  write ( *, '(a)' ) ' '
  r = 0.0D+00
  s = 0.0D+00
  call blend_rs_1dn ( r, s, x, n, stretch_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 1.0D+00
  s = 0.0D+00
  call blend_rs_1dn ( r, s, x, n, stretch_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 0.0D+00
  s = 1.0D+00
  call blend_rs_1dn ( r, s, x, n, stretch_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 1.0D+00
  s = 1.0D+00
  call blend_rs_1dn ( r, s, x, n, stretch_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 0.5D+00
  s = 0.5D+00
  call blend_rs_1dn ( r, s, x, n, stretch_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  return
end
subroutine blend_rst_0dn_identity_test ( )

!*****************************************************************************80
!
!! blend_rst_0dn_identity_test() checks for a gross error in the blend coefficients.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 October 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 3

  external identity_rst
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) t
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'blend_rst_0dn_identity_test():'
  write ( *, '(a)' ) '  Simple identity test to detect gross errors.'
!
!  Test blend_rst_0dn on identity.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Identity test for blend_rst_0dn:'
  write ( *, '(a)' ) ' '
  r = 0.0D+00
  s = 0.0D+00
  t = 0.0D+00
  call blend_rst_0dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 1.0D+00
  s = 0.0D+00
  t = 0.0D+00
  call blend_rst_0dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.0D+00
  s = 1.0D+00
  t = 0.0D+00
  call blend_rst_0dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.0D+00
  s = 0.0D+00
  t = 1.0D+00
  call blend_rst_0dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 1.0D+00
  s = 1.0D+00
  t = 1.0D+00
  call blend_rst_0dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.5D+00
  s = 0.5D+00
  t = 0.5D+00
  call blend_rst_0dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  return
end
subroutine blend_rst_0dn_stretch_test ( )

!*****************************************************************************80
!
!! blend_rst_0dn_stretch_test() checks for simple errors in the blend coefficients.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 October 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 3

  real ( kind = rk ) r
  real ( kind = rk ) s
  external stretch_rst
  real ( kind = rk ) t
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'blend_rst_0dn_stretch_test():'
  write ( *, '(a)' ) '  Shift and stretch test to detect simple errors.'
!
!  Test blend_rst_0dn on shift by (1,2,3), stretch by (4,5,6).
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Shift and stretch test for blend_rst_0dn:'
  write ( *, '(a)' ) ' '
  r = 0.0D+00
  s = 0.0D+00
  t = 0.0D+00
  call blend_rst_0dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 1.0D+00
  s = 0.0D+00
  t = 0.0D+00
  call blend_rst_0dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.0D+00
  s = 1.0D+00
  t = 0.0D+00
  call blend_rst_0dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.0D+00
  s = 0.0D+00
  t = 1.0D+00
  call blend_rst_0dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 1.0D+00
  s = 1.0D+00
  t = 1.0D+00
  call blend_rst_0dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.5D+00
  s = 0.5D+00
  t = 0.5D+00
  call blend_rst_0dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  return
end
subroutine blend_rst_1dn_identity_test ( )

!*****************************************************************************80
!
!! blend_rst_1dn_identity_test() checks for a gross error in the blend coefficients.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 October 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 3

  external identity_rst
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) t
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'blend_rst_1dn_identity_test():'
  write ( *, '(a)' ) '  Simple identity test to detect gross errors.'
!
!  Test blend_rst_1dn on identity.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Identity test for blend_rst_1dn:'
  write ( *, '(a)' ) ' '
  r = 0.0D+00
  s = 0.0D+00
  t = 0.0D+00
  call blend_rst_1dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 1.0D+00
  s = 0.0D+00
  t = 0.0D+00
  call blend_rst_1dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.0D+00
  s = 1.0D+00
  t = 0.0D+00
  call blend_rst_1dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.0D+00
  s = 0.0D+00
  t = 1.0D+00
  call blend_rst_1dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 1.0D+00
  s = 1.0D+00
  t = 1.0D+00
  call blend_rst_1dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.5D+00
  s = 0.5D+00
  t = 0.5D+00
  call blend_rst_1dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  return
end
subroutine blend_rst_1dn_stretch_test ( )

!*****************************************************************************80
!
!! blend_rst_1dn_stretch_test() checks for simple errors in the blend coefficients.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 October 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 3

  real ( kind = rk ) r
  real ( kind = rk ) s
  external stretch_rst
  real ( kind = rk ) t
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'blend_rst_1dn_stretch_test():'
  write ( *, '(a)' ) '  Shift and stretch test to detect simple errors.'
!
!  Test blend_rst_1dn on shift by (1,2,3), stretch by (4,5,6).
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Shift and stretch test for blend_rst_1dn:'
  write ( *, '(a)' ) ' '
  r = 0.0D+00
  s = 0.0D+00
  t = 0.0D+00
  call blend_rst_1dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 1.0D+00
  s = 0.0D+00
  t = 0.0D+00
  call blend_rst_1dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.0D+00
  s = 1.0D+00
  t = 0.0D+00
  call blend_rst_1dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.0D+00
  s = 0.0D+00
  t = 1.0D+00
  call blend_rst_1dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 1.0D+00
  s = 1.0D+00
  t = 1.0D+00
  call blend_rst_1dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)
  r = 0.5D+00
  s = 0.5D+00
  t = 0.5D+00
  call blend_rst_1dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  return
end
subroutine blend_rst_2dn_identity_test ( )

!*****************************************************************************80
!
!! blend_rst_2dn_identity_test() checks for a gross error in the blend coefficients.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 October 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 3

  external identity_rst
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) t
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'blend_rst_2dn_identity_test():'
  write ( *, '(a)' ) '  Simple identity test to detect gross errors.'
!
!  Test blend_rst_2dn on identity.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Identity test for blend_rst_2dn:'
  write ( *, '(a)' ) ' '
  r = 0.0D+00
  s = 0.0D+00
  t = 0.0D+00
  call blend_rst_2dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 1.0D+00
  s = 0.0D+00
  t = 0.0D+00
  call blend_rst_2dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.0D+00
  s = 1.0D+00
  t = 0.0D+00
  call blend_rst_2dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.0D+00
  s = 0.0D+00
  t = 1.0D+00
  call blend_rst_2dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 1.0D+00
  s = 1.0D+00
  t = 1.0D+00
  call blend_rst_2dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.5D+00
  s = 0.5D+00
  t = 0.5D+00
  call blend_rst_2dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  return
end
subroutine blend_rst_2dn_stretch_test ( )

!*****************************************************************************80
!
!! blend_rst_2dn_stretch_test() checks for simple errors in the blend coefficients.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 October 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 3

  real ( kind = rk ) r
  real ( kind = rk ) s
  external stretch_rst
  real ( kind = rk ) t
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'blend_rst_2dn_stretch_test():'
  write ( *, '(a)' ) '  Shift and stretch test to detect simple errors.'
!
!  Test blend_rst_2dn on shift by (1,2,3), stretch by (4,5,6).
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Shift and stretch test for blend_rst_2dn:'
  write ( *, '(a)' ) ' '
  r = 0.0D+00
  s = 0.0D+00
  t = 0.0D+00
  call blend_rst_2dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 1.0D+00
  s = 0.0D+00
  t = 0.0D+00
  call blend_rst_2dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.0D+00
  s = 1.0D+00
  t = 0.0D+00
  call blend_rst_2dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.0D+00
  s = 0.0D+00
  t = 1.0D+00
  call blend_rst_2dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 1.0D+00
  s = 1.0D+00
  t = 1.0D+00
  call blend_rst_2dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.5D+00
  s = 0.5D+00
  t = 0.5D+00
  call blend_rst_2dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  return
end
subroutine cubic_rs ( r, s, i, x )

!*****************************************************************************80
!
!! cubic_rs() evaluates a function of R and S used for some tests.
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
!  Input:
!
!    real R, S, the (local) coordinates of a point.
!
!    integer I, the component of X to be returned.
!
!  Output:
!
!    real X, the value of the I-th component of X at the point whose
!    local coordinates are (R,S).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) x

  x = 20.0D+00 * ( r**2 * s**3 )

  return
end
subroutine quad_rst ( r, s, t, i, x )

!*****************************************************************************80
!
!! quad_rst() evaluates a function of (R,S,T) used for some tests.
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
!  Input:
!
!    real R, S, T, the (local) coordinates of a point.
!
!    integer I, the component of X to be returned.
!
!  Output:
!
!    real X, the value of the I-th component of X at the point whose
!    local coordinates are (R,S,T).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) t
  real ( kind = rk ) x

  x = 18.0D+00 * ( r**2 + s + t )

  return
end
subroutine identity_r ( r, i, x )

!*****************************************************************************80
!
!! identity_r() returns a data component given (R).
!
!  Discussion:
!
!    This is a dummy routine, which simply returns (R).
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
!  Input:
!
!    real R, the coordinate of a point that lies on the
!    boundary of the cube.
!
!    integer I, the component of the data which is to be returned.
!
!  Output:
!
!    real X, the I-th component of the data vector X, evaluated
!    at the point (R), which is on an endpoint of the unit line segment.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  real ( kind = rk ) r
  real ( kind = rk ) x

  if ( i == 1 ) then
    x = r
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'identity_r - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal component index I = ', i
    stop 1
  end if

  return
end
subroutine identity_rs ( r, s, i, x )

!*****************************************************************************80
!
!! identity_rs() returns a data component given (R,S).
!
!  Discussion:
!
!    This is a dummy routine, which simply returns (R,S).
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
!  Input:
!
!    real R, S, the coordinates of a point that lies on the
!    boundary of the square.
!
!    integer I, the component of the data which is to be returned.
!
!  Output:
!
!    real X, the I-th component of the data vector X, evaluated
!    at the point (R,S), which is on a corner, or edge, of the unit square.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) x

  if ( i == 1 ) then
    x = r
  else if ( i == 2 ) then
    x = s
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'identity_rs - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal component index I = ', i
    stop 1
  end if

  return
end
subroutine identity_rst ( r, s, t, i, x )

!*****************************************************************************80
!
!! identity_rst() returns a data component given (R,S,T).
!
!  Discussion:
!
!    This is a dummy routine, which simply returns (R,S,T).
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
!  Input:
!
!    real R, S, T, the coordinates of a point that lies on the
!    boundary of the cube.
!
!    integer I, the component of the data which is to be returned.
!
!  Output:
!
!    real X, the I-th component of the data vector X, evaluated
!    at the point (R,S), which is on a corner, edge or face of the unit cube.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) t
  real ( kind = rk ) x

  if ( i == 1 ) then
    x = r
  else if ( i == 2 ) then
    x = s
  else if ( i == 3 ) then
    x = t
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'identity_rst - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal component index I = ', i
    stop 1
  end if

  return
end
subroutine powers_r ( r, j, x )

!*****************************************************************************80
!
!! powers_r() returns a data component given (R).
!
!  Discussion:
!
!    X(R,J) = R^J
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
!  Input:
!
!    real R, the coordinate of a point.
!
!    integer J, the component of the data.
!
!  Output:
!
!    real X, the J-th component of the data vector X, evaluated
!    at the point (R).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer j
  real ( kind = rk ) r
  real ( kind = rk ) x

  x = r ** j

  return
end
subroutine stretch_r ( r, i, x )

!*****************************************************************************80
!
!! stretch_r() returns a data component given (R).
!
!  Discussion:
!
!    This routine shifts by 1 and stretches by 2.
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
!  Input:
!
!    real R, the coordinate of a point that lies on the
!    boundary of the cube.
!
!    integer I, the component of the data which is to be returned.
!
!  Output:
!
!    real X, the I-th component of the data vector X, evaluated
!    at the point (R), which is on an endpoint of the unit line segment.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  real ( kind = rk ) r
  real ( kind = rk ) x

  if ( i == 1 ) then
    x = 2.0D+00 * r + 1.0D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'stretch_r - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal component index I = ', i
    stop 1
  end if

  return
end
subroutine stretch_rs ( r, s, i, x )

!*****************************************************************************80
!
!! stretch_rs() returns a data component given (R,S).
!
!  Discussion:
!
!    This routine shifts by (1,2) and stretches by (3,4).
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
!  Input:
!
!    real R, S, the coordinates of a point that lies on the
!    boundary of the square.
!
!    integer I, the component of the data which is to be returned.
!
!  Output:
!
!    real X, the I-th component of the data vector X, evaluated
!    at the point (R,S), which is on a corner, or edge, of the unit square.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) x

  if ( i == 1 ) then
    x = 3.0D+00 * r + 1.0D+00
  else if ( i == 2 ) then
    x = 4.0D+00 * s + 2.0D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'stretch_rs - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal component index I = ', i
    stop 1
  end if

  return
end
subroutine stretch_rst ( r, s, t, i, x )

!*****************************************************************************80
!
!! stretch_rst() returns a data component given (R,S,T).
!
!  Discussion:
!
!    This routine shifts by (1,2,3) and stretches by (4,5,6)
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
!  Input:
!
!    real R, S, T, the coordinates of a point that lies on the
!    boundary of the cube.
!
!    integer I, the component of the data which is to be returned.
!
!  Output:
!
!    real X, the I-th component of the data vector X, evaluated
!    at the point (R,S), which is on a corner, edge or face of the unit cube.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) t
  real ( kind = rk ) x

  if ( i == 1 ) then
    x = 4.0D+00 * r + 1.0D+00
  else if ( i == 2 ) then
    x = 5.0D+00 * s + 2.0D+00
  else if ( i == 3 ) then
    x = 6.0D+00 * t + 3.0D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'stretch_rst - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal component index I = ', i
    stop 1
  end if

  return
end
subroutine ellipse_rs ( r, s, i, x )

!*****************************************************************************80
!
!! ellipse_rs() maps the boundary of the unit square to an ellipse.
!
!  Discussion:
!
!    The ellipse is ( 3 * cos ( THETA ), 2 * sin ( THETA ) ).
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
!  Input:
!
!    real R, S, the coordinates of a point that lies on the
!    boundary of the square.
!
!    integer I, the component of the data which is to be returned.
!
!  Output:
!
!    real X, the I-th component of the data vector X, evaluated
!    at the point (R,S), which is on a corner, or edge, of the unit square.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) theta
  real ( kind = rk ) x

  if ( r == 0.0D+00 ) then
    theta = 0.25D+00 * pi * ( 5.0D+00 * ( 1.0D+00 - s ) + 3.0D+00 * s  )
  else if ( r == 1.0D+00 ) then
    theta = 0.25D+00 * pi * ( - 1.0D+00 * ( 1.0D+00 - s ) + 1.0D+00 * s )
  else if ( s == 0.0D+00 ) then
    theta = 0.25D+00 * pi * ( 5.0D+00 * ( 1.0D+00 - r ) + 7.0D+00 * r )
  else if ( s == 1.0D+00 ) then
    theta = 0.25D+00 * pi * ( 3.0D+00 * ( 1.0D+00 - r ) + 1.0D+00 * r )
  end if

  if ( i == 1 ) then

    x = 3.0D+00 * cos ( theta )

  else if ( i == 2 ) then

    x = 2.0D+00 * sin ( theta )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ellipse_rs - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal component index I = ', i
    stop 1

  end if

  return
end
subroutine sphere_rst ( r, s, t, i, x )

!*****************************************************************************80
!
!! sphere_rst() maps the boundary of the unit cube to a sphere.
!
!  Discussion:
!
!    The sphere is
!      x = cos ( theta ) * cos ( phi )
!      y = sin ( theta ) * cos ( phi )
!      z = sin ( phi )
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
!  Input:
!
!    real R, S, the coordinates of a point that lies on the
!    boundary of the square.
!
!    integer I, the component of the data which is to be returned.
!
!  Output:
!
!    real X, the I-th component of the data vector X, evaluated
!    at the point (R,S), which is on a corner, or edge, of the unit square.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  real ( kind = rk ) norm
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) t
  real ( kind = rk ) x
!
!  Compute length of vector from ( 0.5, 0.5, 0.5 ) to ( r, s, t )
!
  norm = sqrt ( ( r - 0.5D+00 )**2 + ( s - 0.5D+00 )**2 + ( t - 0.5D+00 )**2 )
!
!  Compute ( x, y, z ) coordinates of a point on the sphere
!  ( x - 0.5 )^2 + ( y - 0.5 )^2 + ( z - 0.5 )^2 = 0.25 that is
!  the projection of the point ( r, s, t ) on the unit cube.
!
  if ( i == 1 ) then

    x = 0.5D+00 + 0.5D+00 * ( r - 0.5D+00 ) / norm

  else if ( i == 2 ) then

    x = 0.5D+00 + 0.5D+00 * ( s - 0.5D+00 ) / norm

  else if ( i == 3 ) then

    x = 0.5D+00 + 0.5D+00 * ( t - 0.5D+00 ) / norm

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'sphere_rst - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal component index I = ', i
    stop 1

  end if

  return
end

