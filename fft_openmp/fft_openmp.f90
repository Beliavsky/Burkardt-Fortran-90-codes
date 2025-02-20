function c8ggl ( seed )

!*****************************************************************************80
!
!! c8ggl() generates uniformly distributed pseudorandom complex values. 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 August 2020
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wesley Petersen, Peter Arbenz, 
!    Introduction to Parallel Computing - A practical guide with examples in C,
!    Oxford University Press,
!    ISBN: 0-19-851576-6,
!    LC: QA76.58.P47.
!
!  Input:
!
!    real ( kind = rk ) SEED, a seed for the sequence.
!
!  Input:
!
!    real ( kind = rk ) SEED, an updated seed value.
!
!    complex ( kind = ck ) C8GGL, a pseudorandom value.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  complex ( kind = ck ) c8ggl
  real ( kind = rk ), parameter :: d2 = 0.2147483647D+10
  real ( kind = rk ) r1
  real ( kind = rk ) r2
  real ( kind = rk ) seed

  seed = mod ( 16807.0D+00 * real ( seed, kind = rk ), d2 )
  r1 = real ( ( seed - 1.0D+00 ) / ( d2 - 1.0D+00 ), kind = rk )
  seed = mod ( 16807.0D+00 * real ( seed, kind = rk ), d2 )
  r2 = real ( ( seed - 1.0D+00 ) / ( d2 - 1.0D+00 ), kind = rk )

  c8ggl = cmplx ( r1, r2, kind = ck )

  return
end
subroutine cfft_1d ( n, x, y, w, sgn )

!*****************************************************************************80
!
!! cfft_1d performs a complex Fast Fourier Transform of 1D data.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
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
!    Input, real ( kind = rk ) SGN:
!    +1 for a "forward" FFT,
!    -1 for a "backward" FFT.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
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
!  Licensing:
!
!    This code is distributed under the MIT license.
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

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
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

!$omp parallel &
!$omp   shared ( aw, n, w ) &
!$omp   private ( arg, i )
!$omp do
  do i = 1, n2
    arg = aw * real ( i - 1, kind = rk )
    w(2*i-1) = cos ( arg )
    w(2*i)   = sin ( arg )
  end do
!$omp end do
!$omp end parallel

  return
end
function r82ggl ( seed )

!*****************************************************************************80
!
!! r82ggl generates uniformly distributed pseudorandom numbers. 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
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
!  Input:
!
!    real ( kind = rk ) SEED, used as a seed for the sequence.
!
!  Parameters:
!
!    real ( kind = rk ) SEED, an updated seed value.
!
!    real ( kind = rk ) r82ggl, the next pseudorandom value.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), parameter :: d2 = 0.2147483647D+10
  real ( kind = rk ) r82ggl
  real ( kind = rk ) seed

  seed = mod ( 16807.0D+00 * real ( seed, kind = rk ), d2 )
  r82ggl = real ( ( seed - 1.0D+00 ) / ( d2 - 1.0D+00 ), kind = rk )

  return
end
subroutine step ( n, mj, a, b, c, d, w, sgn )

!*****************************************************************************80
!
!! STEP carries out one step of the workspace version of cfft_1d.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
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
!  Input:
!
!    integer N, the size of the arrays.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
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

!$omp parallel &
!$omp   shared ( a, b, c, d, lj, mj, mj2, sgn, w ) &
!$omp   private ( ambr, ambu, j, ja, jb, jc, jd, jw, k, wjw )
!$omp do
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
!$omp end do
!$omp end parallel

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

