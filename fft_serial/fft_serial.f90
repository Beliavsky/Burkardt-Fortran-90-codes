program main

!*****************************************************************************80
!
!! MAIN is the main program for FFT_SERIAL.
!
!  Discussion:
!
!    The complex data in an N vector is stored as pairs of values in a
!    real vector of length 2*N.
!
!  Modified:
!
!    16 July 2008
!
!  Author:
!
!    Original C version by Wesley Petersen.
!    This version by John Burkardt.
!
!  Reference:
!
!    Wesley Petersen, Peter Arbenz, 
!    Introduction to Parallel Computing - A practical guide with examples in C,
!    Oxford University Press,
!    ISBN: 0-19-851576-6,
!    LC: QA76.58.P47.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) ctime
  real ( kind = rk ) ctime1
  real ( kind = rk ) ctime2
  real ( kind = rk ) error
  logical first
  real ( kind = rk ) flops
  real ( kind = rk ) fnm1
  real ( kind = rk ) ggl
  integer i
  integer icase
  integer it
  integer ln2
  real ( kind = rk ) mflops
  integer n
  integer nits
  real ( kind = rk ) seed
  real ( kind = rk ) sgn
  real ( kind = rk ), allocatable, dimension ( : ) :: w
  real ( kind = rk ), allocatable, dimension ( : ) :: x
  real ( kind = rk ), allocatable, dimension ( : ) :: y
  real ( kind = rk ), allocatable, dimension ( : ) :: z
  real ( kind = rk ) z0
  real ( kind = rk ) z1

  nits = 10000

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FFT_SERIAL():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) &
    '  Demonstrate an implementation of the Fast Fourier Transform'
  write ( *, '(a)' ) '  of a complex data vector.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Accuracy check:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    FFT ( FFT ( X(1:N) ) ) == N * X(1:N)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             N      NITS    Error         Time' // &
    '          Time/Call     MFLOPS'
  write ( *, '(a)' ) ' '

  seed  = 331.0D+00
  n = 1
!
!  LN2 is the log base 2 of N.  Each increase of LN2 doubles N.
!
  do ln2 = 1, 20

    n = 2 * n
!
!  Allocate storage for the complex arrays W, X, Y, Z.  
!
!  We handle the complex arithmetic,
!  and store a complex number as a pair of floats, a complex vector as a doubly
!  dimensioned array whose second dimension is 2. 
!
    allocate ( w(1:n) )
    allocate ( x(1:2*n) )
    allocate ( y(1:2*n) )
    allocate ( z(1:2*n) )

    first = .true.

    do icase = 0, 1

      if ( first ) then

        do i = 1, 2 * n - 1, 2

          z0 = ggl ( seed )
          z1 = ggl ( seed )
          x(i) = z0
          z(i) = z0
          x(i+1) = z1
          z(i+1) = z1
        end do

      else

        do i = 1, 2 * n - 1, 2
          z0 = 0.0D+00
          z1 = 0.0D+00
          x(i) = z0
          z(i) = z0
          x(i+1) = z1
          z(i+1) = z1
        end do

      end if
!
!  Initialize the sine and cosine tables.
!
      call cffti ( n, w )
!
!  Transform forward, then back.
!
      if ( first ) then

        sgn = + 1.0D+00
        call cfft2 ( n, x, y, w, sgn )
        sgn = - 1.0D+00
        call cfft2 ( n, y, x, w, sgn )
! 
!  Results should be same as initial vector multiplied by N. 
!
        fnm1 = 1.0D+00 / real ( n, kind = rk )

        error = 0.0D+00
        do i = 1, 2 * n - 1, 2
          error = error &
          + ( z(i)   - fnm1 * x(i) )**2 &
          + ( z(i+1) - fnm1 * x(i+1) )**2
        end do
        error = sqrt ( fnm1 * error )

        first = .false.

      else

        call cpu_time ( ctime1 )

        do it = 1, nits

          sgn = + 1.0D+00
          call cfft2 ( n, x, y, w, sgn )
          sgn = - 1.0D+00
          call cfft2 ( n, y, x, w, sgn )

        end do

        call cpu_time ( ctime2 )
        ctime = ctime2 - ctime1

        flops = 2.0D+00 * nits * ( 5.0D+00 * n * ln2 )

        mflops = flops / 1.0D+06 / ctime

        write ( *, '(2x,i12,2x,i8,2x,g14.6,2x,g12.4,2x,g12.4,2x,g12.4)' ) &
          n, nits, error, ctime, ctime / real ( 2 * nits, kind = rk ), mflops

      end if

    end do

    if ( mod ( ln2, 4 ) == 0 ) then
      nits = nits / 10
    end if

    if ( nits < 1 ) then
      nits = 1
    end if
!
!  Free memory.
!
    deallocate ( w )
    deallocate ( x )
    deallocate ( y )
    deallocate ( z )

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FFT_SERIAL:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine cfft2 ( n, x, y, w, sgn )

!*****************************************************************************80
!
!! CFFT2 performs a complex Fast Fourier Transform.
!
!  Modified:
!
!    27 April 2008
!
!  Author:
!
!    Original C version by Wesley Petersen.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wesley Petersen, Peter Arbenz, 
!    Introduction to Parallel Computing - A practical guide with examples in C,
!    Oxford University Press,
!    ISBN: 0-19-851576-6,
!    LC: QA76.58.P47.
!
!  Parameters:
!
!    Input, integer N, the size of the array to be transformed.
!
!    Input/output, real ( kind = rk ) X(2*N), the data to be transformed.  
!    On output, the contents of X have been overwritten by work information.
!
!    Output, real ( kind = rk ) Y(2*N), the forward or backward FFT of X.
!
!    Input, real ( kind = rk ) W(N), a table of sines and cosines.
!
!    Input, real ( kind = rk ) SGN, is +1 for a "forward" FFT 
!    and -1 for a "backward" FFT.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer j
  integer m
  integer mj
  real ( kind = rk ) sgn
  logical tgle
  real ( kind = rk ) w(n)
  real ( kind = rk ) x(2*n)
  real ( kind = rk ) y(2*n)

   m = int ( ( log ( real ( n, kind = rk ) ) / log ( 1.99D+00 ) ) )
   mj = 1
!
!  Toggling switch for work array.
!
  tgle = .true.
  call step ( n, mj, x(1), x((n/2)*2+1), y(1), y(mj*2+1), w, sgn )

  if ( n == 2 ) then
    return
  end if

  do j = 1, m - 2  

    mj = mj * 2

    if ( tgle ) then
      call step ( n, mj, y(1), y((n/2)*2+1), x(1), x(mj*2+1), w, sgn )
      tgle = .false.
    else
      call step ( n, mj, x(1), x((n/2)*2+1), y(1), y(mj*2+1), w, sgn )
      tgle = .true.
    end if

  end do
!
!  Last pass through data: move Y to X if needed. 
!
  if ( tgle ) then
    x(1:2*n) = y(1:2*n)
  end if

  mj = n / 2
  call step ( n, mj, x(1), x((n/2)*2+1), y(1), y(mj*2+1), w, sgn )

  return
end
subroutine cffti ( n, w )

!*****************************************************************************80
!
!! CFFTI sets up sine and cosine tables needed for the FFT calculation.
!
!  Modified:
!
!    28 April 2008
!
!  Author:
!
!    Original C version by Wesley Petersen.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wesley Petersen, Peter Arbenz, 
!    Introduction to Parallel Computing - A practical guide with examples in C,
!    Oxford University Press,
!    ISBN: 0-19-851576-6,
!    LC: QA76.58.P47.
!
!  Parameters:
!
!    Input, integer N, the size of the array to be transformed.
!
!    Output, real ( kind = rk ) W(N), a table of sines and cosines.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) arg
  real ( kind = rk ) aw
  integer i
  integer n2
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) w(n)

  n2 = n / 2
  aw = 2.0D+00 * pi / real ( n, kind = rk )

  do i = 1, n2
    arg = aw * real ( i - 1, kind = rk )
    w(2*i-1) = cos ( arg )
    w(2*i)   = sin ( arg )
  end do

  return
end
function ggl ( seed )

!*****************************************************************************80
!
!! GGL generates uniformly distributed pseudorandom numbers. 
!
!  Modified:
!
!    16 July 2008
!
!  Author:
!
!    Original C version by Wesley Petersen, M Troyer, I Vattulainen.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wesley Petersen, Peter Arbenz, 
!    Introduction to Parallel Computing - A practical guide with examples in C,
!    Oxford University Press,
!    ISBN: 0-19-851576-6,
!    LC: QA76.58.P47.
!
!  Parameters:
!
!    Input/output, real ( kind = rk ) SEED, used as a seed for the sequence.
!
!    Output, real ( kind = rk ) GGL, the next pseudorandom value.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), parameter :: d2 = 0.2147483647D+10
  real ( kind = rk ) ggl
  real ( kind = rk ) seed
  real ( kind = rk ) t

  t = mod ( 16807.0D+00 * seed, d2 )
  seed = t
  ggl = ( t - 1.0D+00 ) / ( d2 - 1.0D+00 )

  return
end
subroutine step ( n, mj, a, b, c, d, w, sgn )

!*****************************************************************************80
!
!! STEP carries out one step of the workspace version of CFFT2.
!
!  Modified:
!
!    27 April 2008
!
!  Author:
!
!    Original C version by Wesley Petersen.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Wesley Petersen, Peter Arbenz, 
!    Introduction to Parallel Computing - A practical guide with examples in C,
!    Oxford University Press,
!    ISBN: 0-19-851576-6,
!    LC: QA76.58.P47.
!
!  Parameters:
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) ambr
  real ( kind = rk ) ambu
  real ( kind = rk ) b(n)
  real ( kind = rk ) c(n)
  real ( kind = rk ) d(n)
  integer j
  integer ja
  integer jb
  integer jc
  integer jd
  integer jw
  integer k
  integer lj
  integer mj
  integer mj2
  real ( kind = rk ) sgn
  real ( kind = rk ) w(n)
  real ( kind = rk ) wjw(2)

  mj2 = 2 * mj
  lj  = n / mj2

  do j = 0, lj - 1

    jw = j * mj
    ja  = jw
    jb  = ja
    jc  = j * mj2
    jd  = jc

    wjw(1) = w(jw*2+1) 
    wjw(2) = w(jw*2+2)

    if ( sgn < 0.0D+00 ) then
      wjw(2) = - wjw(2)
    end if

    do k = 0, mj - 1

      c((jc+k)*2+1) = a((ja+k)*2+1) + b((jb+k)*2+1)
      c((jc+k)*2+2) = a((ja+k)*2+2) + b((jb+k)*2+2)

      ambr = a((ja+k)*2+1) - b((jb+k)*2+1)
      ambu = a((ja+k)*2+2) - b((jb+k)*2+2)

      d((jd+k)*2+1) = wjw(1) * ambr - wjw(2) * ambu
      d((jd+k)*2+2) = wjw(2) * ambr + wjw(1) * ambu

    end do
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
