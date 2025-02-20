program main

!*****************************************************************************80
!
!! bvls_test() tests bvls().
!
!  Discussion:
!
!    This program demonstrates the use of BVLS for solving least squares
!    problems with bounds on the variables.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 June 2014
!
!  Author:
!
!    Original FORTRAN90 version by Charles Lawson, Richard Hanson.
!    This version by John Burkardt.
!
!  Reference:
!
!    Charles Lawson, Richard Hanson,
!    Solving Least Squares Problems,
!    SIAM, 1995,
!    ISBN: 0898713560,
!    LC: QA275.L38.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'bvls_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test bvls().'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'bvls_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 runs test case 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 June 2014
!
!  Author:
!
!    Original FORTRAN90 version by Charles Lawson, Richard Hanson.
!    This version by John Burkardt.
!
!  Reference:
!
!    Charles Lawson, Richard Hanson,
!    Solving Least Squares Problems,
!    SIAM, 1995,
!    ISBN: 0898713560,
!    LC: QA275.L38.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: m = 2
  integer, parameter :: n = 2
  integer, parameter :: jstep = 5

  real ( kind = rk8 ) a(m,n)
  real ( kind = rk8 ) a2(m,n)
  real ( kind = rk8 ) b(m)
  real ( kind = rk8 ) b2(m)
  real ( kind = rk8 ) bnd(2,n)
  integer i
  integer ierr
  integer index(n)
  integer j
  integer j1
  integer j2
  integer nsetp
  real ( kind = rk8 ) rnorm
  real ( kind = rk8 ) unbnd
  real ( kind = rk8 ) w(n)
  real ( kind = rk8 ) x(n)

  save bnd
  save unbnd

  data ((bnd(i,j),i=1,2),j=1,2)/&
    1.0D+00, 2.0D+00, &
    3.0D+00, 4.0D+00 /

  data unbnd / 1.0D+06 /

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'

  where ( bnd(1,1:n) == unbnd ) 
    bnd(1,1:n) = -huge(1.0D+00)
  endwhere

  where ( bnd(2,1:n) == unbnd ) 
    bnd(2,1:n) =  huge(1.0D+00)
  endwhere

  write ( *, '(a)' ) ' '
  write ( *, '(a,i5,a,i5,a,g17.5)') &
    '  M =', m,',   N =', n,',   UNBND =', unbnd

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Bounds:'

  do j1 = 1, n, jstep
    j2 = min ( j1 - 1 + jstep, n )
    write ( *, '(a)' ) ' '
    write ( *, '(2x,5g14.6)' ) bnd(1,j1:j2)
    write ( *, '(2x,5g14.6)' ) bnd(2,j1:j2)
  end do

  call random_number ( harvest = b(1:m) )
  call random_number ( harvest = a(1:m,1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Matrix A:'

  do j1 = 1, n, jstep
    j2 = min ( j1 - 1 + jstep, n )
    write ( *, '(a)' ) ' '
    do i = 1,m
      write ( *, '(2x,5g14.6)' ) a(i,j1:j2)
    end do
  end do

  b2(1:m) = b(1:m)
  a2(1:m,1:n) = a(1:m,1:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  RHS B:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,5g14.6)' ) b(1:m)

  call bvls ( m, n, a2, b2, bnd, x, rnorm, nsetp, w, index, ierr )

  call bvls_report ( m, n, a, b, bnd, x, rnorm, nsetp, w, index, ierr )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 runs test case 2.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 June 2014
!
!  Author:
!
!    Original FORTRAN90 version by Charles Lawson, Richard Hanson.
!    This version by John Burkardt.
!
!  Reference:
!
!    Charles Lawson, Richard Hanson,
!    Solving Least Squares Problems,
!    SIAM, 1995,
!    ISBN: 0898713560,
!    LC: QA275.L38.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: m = 2
  integer, parameter :: n = 4
  integer, parameter :: jstep = 5

  real ( kind = rk8 ) a(m,n)
  real ( kind = rk8 ) a2(m,n)
  real ( kind = rk8 ) b(m)
  real ( kind = rk8 ) b2(m)
  real ( kind = rk8 ) bnd(2,n)
  integer i
  integer ierr
  integer index(n)
  integer j
  integer j1
  integer j2
  integer nsetp
  real ( kind = rk8 ) rnorm
  real ( kind = rk8 ) unbnd
  real ( kind = rk8 ) w(n)
  real ( kind = rk8 ) x(n)

  save bnd
  save unbnd

  data ((bnd(i,j),i=1,2),j=1,n)/&
    0.0D+00, 10.0D+00, &
    0.0D+00, 10.0D+00, &
    0.0D+00, 10.0D+00, &
    0.0D+00, 10.0D+00 /

  data unbnd / 1.0D+06 /

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'

  where ( bnd(1,1:n) == unbnd ) 
    bnd(1,1:n) = -huge(1.0D+00)
  endwhere

  where ( bnd(2,1:n) == unbnd ) 
    bnd(2,1:n) =  huge(1.0D+00)
  endwhere

  write ( *, '(a)' ) ' '
  write ( *, '(a,i5,a,i5,a,g17.5)') &
    '  M =', m,',   N =', n,',   UNBND =', unbnd

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Bounds:'

  do j1 = 1, n, jstep
    j2 = min ( j1 - 1 + jstep, n )
    write ( *, '(a)' ) ' '
    write ( *, '(2x,5g14.6)' ) bnd(1,j1:j2)
    write ( *, '(2x,5g14.6)' ) bnd(2,j1:j2)
  end do

  call random_number ( harvest = b(1:m) )
  call random_number ( harvest = a(1:m,1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Matrix A:'

  do j1 = 1, n, jstep
    j2 = min ( j1 - 1 + jstep, n )
    write ( *, '(a)' ) ' '
    do i = 1,m
      write ( *, '(2x,5g14.6)' ) a(i,j1:j2)
    end do
  end do

  b2(1:m) = b(1:m)
  a2(1:m,1:n) = a(1:m,1:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  RHS B:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,5g14.6)' ) b(1:m)

  call bvls ( m, n, a2, b2, bnd, x, rnorm, nsetp, w, index, ierr )

  call bvls_report ( m, n, a, b, bnd, x, rnorm, nsetp, w, index, ierr )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 runs test case 3.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 June 2014
!
!  Author:
!
!    Original FORTRAN90 version by Charles Lawson, Richard Hanson.
!    This version by John Burkardt.
!
!  Reference:
!
!    Charles Lawson, Richard Hanson,
!    Solving Least Squares Problems,
!    SIAM, 1995,
!    ISBN: 0898713560,
!    LC: QA275.L38.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: m = 4
  integer, parameter :: n = 2
  integer, parameter :: jstep = 5

  real ( kind = rk8 ) a(m,n)
  real ( kind = rk8 ) a2(m,n)
  real ( kind = rk8 ) b(m)
  real ( kind = rk8 ) b2(m)
  real ( kind = rk8 ) bnd(2,n)
  integer i
  integer ierr
  integer index(n)
  integer j
  integer j1
  integer j2
  integer nsetp
  real ( kind = rk8 ) rnorm
  real ( kind = rk8 ) unbnd
  real ( kind = rk8 ) w(n)
  real ( kind = rk8 ) x(n)

  save bnd
  save unbnd

  data ((bnd(i,j),i=1,2),j=1,2)/&
       0.0D+00,  100.0D+00, &
    -100.0D+00,  100.0D+00 /

  data unbnd / 1.0D+06 /

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'

  where ( bnd(1,1:n) == unbnd ) 
    bnd(1,1:n) = -huge(1.0D+00)
  endwhere

  where ( bnd(2,1:n) == unbnd ) 
    bnd(2,1:n) =  huge(1.0D+00)
  endwhere

  write ( *, '(a)' ) ' '
  write ( *, '(a,i5,a,i5,a,g17.5)') &
    '  M =', m,',   N =', n,',   UNBND =', unbnd

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Bounds:'

  do j1 = 1, n, jstep
    j2 = min ( j1 - 1 + jstep, n )
    write ( *, '(a)' ) ' '
    write ( *, '(2x,5g14.6)' ) bnd(1,j1:j2)
    write ( *, '(2x,5g14.6)' ) bnd(2,j1:j2)
  end do

  call random_number ( harvest = b(1:m) )
  call random_number ( harvest = a(1:m,1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Matrix A:'

  do j1 = 1, n, jstep
    j2 = min ( j1 - 1 + jstep, n )
    write ( *, '(a)' ) ' '
    do i = 1,m
      write ( *, '(2x,5g14.6)' ) a(i,j1:j2)
    end do
  end do

  b2(1:m) = b(1:m)
  a2(1:m,1:n) = a(1:m,1:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  RHS B:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,5g14.6)' ) b(1:m)

  call bvls ( m, n, a2, b2, bnd, x, rnorm, nsetp, w, index, ierr )

  call bvls_report ( m, n, a, b, bnd, x, rnorm, nsetp, w, index, ierr )

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 runs test case 4.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 June 2014
!
!  Author:
!
!    Original FORTRAN90 version by Charles Lawson, Richard Hanson.
!    This version by John Burkardt.
!
!  Reference:
!
!    Charles Lawson, Richard Hanson,
!    Solving Least Squares Problems,
!    SIAM, 1995,
!    ISBN: 0898713560,
!    LC: QA275.L38.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 10
  integer, parameter :: jstep = 5

  real ( kind = rk8 ) a(m,n)
  real ( kind = rk8 ) a2(m,n)
  real ( kind = rk8 ) b(m)
  real ( kind = rk8 ) b2(m)
  real ( kind = rk8 ) bnd(2,n)
  integer i
  integer ierr
  integer index(n)
  integer j
  integer j1
  integer j2
  integer nsetp
  real ( kind = rk8 ) rnorm
  real ( kind = rk8 ) unbnd
  real ( kind = rk8 ) w(n)
  real ( kind = rk8 ) x(n)

  save bnd
  save unbnd

  data ((bnd(i,j),i=1,2),j=1,10)/&
    0.0D+00,     0.0D+00, &
   -0.3994D+00, -0.3994D+00, &
   -1.0D+00,     1.0D+00, &
   -0.3D+00,    -0.2D+00, &
   21.0D+00,    22.0D+00, &
   -4.0D+00,    -3.0D+00, &
    45.0D+00,   46.0D+00, &
    100.0D+00, 101.0D+00, &
    1.0D+06,     1.0D+06, &
   -1.0D+00,     1.0D+00 /

  data unbnd / 1.0D+06 /

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'

  where ( bnd(1,1:n) == unbnd ) 
    bnd(1,1:n) = -huge(1.0D+00)
  endwhere

  where ( bnd(2,1:n) == unbnd ) 
    bnd(2,1:n) =  huge(1.0D+00)
  endwhere

  write ( *, '(a)' ) ' '
  write ( *, '(a,i5,a,i5,a,g17.5)') &
    '  M =', m,',   N =', n,',   UNBND =', unbnd

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Bounds:'

  do j1 = 1, n, jstep
    j2 = min ( j1 - 1 + jstep, n )
    write ( *, '(a)' ) ' '
    write ( *, '(2x,5g14.6)' ) bnd(1,j1:j2)
    write ( *, '(2x,5g14.6)' ) bnd(2,j1:j2)
  end do

  call random_number ( harvest = b(1:m) )
  call random_number ( harvest = a(1:m,1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Matrix A:'

  do j1 = 1, n, jstep
    j2 = min ( j1 - 1 + jstep, n )
    write ( *, '(a)' ) ' '
    do i = 1,m
      write ( *, '(2x,5g14.6)' ) a(i,j1:j2)
    end do
  end do

  b2(1:m) = b(1:m)
  a2(1:m,1:n) = a(1:m,1:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  RHS B:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,5g14.6)' ) b(1:m)

  call bvls ( m, n, a2, b2, bnd, x, rnorm, nsetp, w, index, ierr )

  call bvls_report ( m, n, a, b, bnd, x, rnorm, nsetp, w, index, ierr )

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 runs test case 5.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 June 2014
!
!  Author:
!
!    Original FORTRAN90 version by Charles Lawson, Richard Hanson.
!    This version by John Burkardt.
!
!  Reference:
!
!    Charles Lawson, Richard Hanson,
!    Solving Least Squares Problems,
!    SIAM, 1995,
!    ISBN: 0898713560,
!    LC: QA275.L38.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: m = 10
  integer, parameter :: n = 5
  integer, parameter :: jstep = 5

  real ( kind = rk8 ) a(m,n)
  real ( kind = rk8 ) a2(m,n)
  real ( kind = rk8 ) b(m)
  real ( kind = rk8 ) b2(m)
  real ( kind = rk8 ) bnd(2,n)
  integer i
  integer ierr
  integer index(n)
  integer j
  integer j1
  integer j2
  integer nsetp
  real ( kind = rk8 ) rnorm
  real ( kind = rk8 ) unbnd
  real ( kind = rk8 ) w(n)
  real ( kind = rk8 ) x(n)

  save bnd
  save unbnd

  data ((bnd(i,j),i=1,2),j=1,5)/&
    0.0D+00,   1.0D+00, &
   -1.0D+00,   0.0D+00, &
    0.0D+00,   1.0D+00, &
    0.3D+00,   0.4D+00, &
    0.048D+00, 0.049D+00 /

  data unbnd / 1.0D+06 /

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'

  where ( bnd(1,1:n) == unbnd ) 
    bnd(1,1:n) = -huge(1.0D+00)
  endwhere

  where ( bnd(2,1:n) == unbnd ) 
    bnd(2,1:n) =  huge(1.0D+00)
  endwhere

  write ( *, '(a)' ) ' '
  write ( *, '(a,i5,a,i5,a,g17.5)') &
    '  M =', m,',   N =', n,',   UNBND =', unbnd

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Bounds:'

  do j1 = 1, n, jstep
    j2 = min ( j1 - 1 + jstep, n )
    write ( *, '(a)' ) ' '
    write ( *, '(2x,5g14.6)' ) bnd(1,j1:j2)
    write ( *, '(2x,5g14.6)' ) bnd(2,j1:j2)
  end do

  call random_number ( harvest = b(1:m) )
  call random_number ( harvest = a(1:m,1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Matrix A:'

  do j1 = 1, n, jstep
    j2 = min ( j1 - 1 + jstep, n )
    write ( *, '(a)' ) ' '
    do i = 1,m
      write ( *, '(2x,5g14.6)' ) a(i,j1:j2)
    end do
  end do

  b2(1:m) = b(1:m)
  a2(1:m,1:n) = a(1:m,1:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  RHS B:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,5g14.6)' ) b(1:m)

  call bvls ( m, n, a2, b2, bnd, x, rnorm, nsetp, w, index, ierr )

  call bvls_report ( m, n, a, b, bnd, x, rnorm, nsetp, w, index, ierr )

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 runs test case 6.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 June 2014
!
!  Author:
!
!    Original FORTRAN90 version by Charles Lawson, Richard Hanson.
!    This version by John Burkardt.
!
!  Reference:
!
!    Charles Lawson, Richard Hanson,
!    Solving Least Squares Problems,
!    SIAM, 1995,
!    ISBN: 0898713560,
!    LC: QA275.L38.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: m = 6
  integer, parameter :: n = 4
  integer, parameter :: jstep = 5

  real ( kind = rk8 ) a(m,n)
  real ( kind = rk8 ) a2(m,n)
  real ( kind = rk8 ) b(m)
  real ( kind = rk8 ) b2(m)
  real ( kind = rk8 ) bnd(2,n)
  integer i
  integer ierr
  integer index(n)
  integer j
  integer j1
  integer j2
  integer nsetp
  real ( kind = rk8 ) rnorm
  real ( kind = rk8 ) unbnd
  real ( kind = rk8 ) w(n)
  real ( kind = rk8 ) x(n)

  save bnd
  save unbnd

  data ((bnd(i,j),i=1,2),j=1,4)/&
    -100.0D+00, 100.0D+00, &
     999.0D+00, 999.0D+00, &
     999.0D+00, 999.0D+00, &
     999.0D+00, 999.0D+00 /

  data unbnd / 999.0D+00 /

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'

  where ( bnd(1,1:n) == unbnd ) 
    bnd(1,1:n) = -huge(1.0D+00)
  endwhere

  where ( bnd(2,1:n) == unbnd ) 
    bnd(2,1:n) =  huge(1.0D+00)
  endwhere

  write ( *, '(a)' ) ' '
  write ( *, '(a,i5,a,i5,a,g17.5)') &
    '  M =', m,',   N =', n,',   UNBND =', unbnd

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Bounds:'

  do j1 = 1, n, jstep
    j2 = min ( j1 - 1 + jstep, n )
    write ( *, '(a)' ) ' '
    write ( *, '(2x,5g14.6)' ) bnd(1,j1:j2)
    write ( *, '(2x,5g14.6)' ) bnd(2,j1:j2)
  end do

  call random_number ( harvest = b(1:m) )
  call random_number ( harvest = a(1:m,1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Matrix A:'

  do j1 = 1, n, jstep
    j2 = min ( j1 - 1 + jstep, n )
    write ( *, '(a)' ) ' '
    do i = 1,m
      write ( *, '(2x,5g14.6)' ) a(i,j1:j2)
    end do
  end do

  b2(1:m) = b(1:m)
  a2(1:m,1:n) = a(1:m,1:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  RHS B:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,5g14.6)' ) b(1:m)

  call bvls ( m, n, a2, b2, bnd, x, rnorm, nsetp, w, index, ierr )

  call bvls_report ( m, n, a, b, bnd, x, rnorm, nsetp, w, index, ierr )

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

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end

