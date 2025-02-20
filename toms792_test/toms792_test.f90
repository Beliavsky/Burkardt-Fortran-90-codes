program main

!*****************************************************************************80
!
!! toms792_test() tests toms792().
!
!  Discussion:
!
!    This program computes interpolation errors using the
!    scattered data package CSHEP2D for each of ten test
!    functions and a 33 by 33 uniform grid of interpolation
!    points in the unit square.
!
!    This program uses Subroutines TESTDT and TSTFN1 from
!    ACM Algorithm SURVEY to generate a node set and and the
!    test function values.
!
!  Local:
!
!    integer NFUN, the number of test functions.
!
!    integer NSET, the number of node sets.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nmax = 100
  integer, parameter :: nrmax = 10
  integer, parameter :: ni = 33

  real ( kind = rk ) a(9,nmax)
  real ( kind = rk ) ft(ni,ni)
  integer i
  integer ier
  integer j
  integer k
  integer kf
  integer kff
  integer kfl
  integer ks
  integer lcell(nrmax,nrmax)
  integer lnext(nmax)
  integer n
  integer nc
  integer, parameter :: nfun = 10
  integer nfun2
  integer np
  integer nr
  integer, parameter :: nset = 5
  integer nw
  integer nwmax
  real ( kind = rk ) p(ni)
  real ( kind = rk ) rw(nmax)
  real ( kind = rk ) w(nmax)
  real ( kind = rk ) x(nmax)
  real ( kind = rk ) y(nmax)

  double precision dum, dx, dy, ermax, ermean, pw, &
                       rmax, ssa, sse, ssm, sum, xmin, ymin
  double precision cs2val
!
! Input format:
!
  100 format (I2)
!
! Get a user-specified node set number KS.
!
  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS792_TEST():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TOMS792 library.'
  write ( *, '(a)' ) ' '

    1 write (*,110) nset
  110 format (///13X,'CS2TST:  CSHEP2D Test Program'// &
              5X,'Specify a node set number in the range 1', &
                 ' to ',I2,':'/)
      read (*,100,end=999,err=1) ks
      if (ks .lt. 1  .or.  ks .gt. nset) go to 1
!
! Copy N and the nodal coordinates for node set KS.
!
      call testdt ( ks, n, x, y )

      if (n .lt. 10  .or.  n .gt. nmax) then
        write (*,500) n, nmax
        stop
      end if
!
! Allow the user to specify a range of function numbers.
!
      write (*,120) nfun
  120 format (//5X,'Specify the first test function ', &
                   '(1 to ',I2,'):'/)
      read (*,100,err=1) kff
      if (kff .lt. 1  .or.  kff .gt. nfun) go to 1
      write ( *, * ) '  kff = ', kff

      write (*,130) kff, nfun
  130 format (//5X,'Specify the last test function (', &
                   I2,' to ',I2,'):'/)
      read (*,100,err=1) kfl
      if (kfl .lt. kff  .or.  kfl .gt. nfun) go to 1
      write ( *, * ) '  kfl = ', kfl

      nfun2 = kfl-kff+1
!
! Input NC, NW, and NR from the console.
!
      nwmax = min(40,n-1)
    2 write (*,140) n
  140 format (//5x,'n =',i4//5x, &
              'Specify the number of nodes NC for the ', &
              'least squares fits.'/5X,'NC = 17 is ', &
              'recommended.  NC GE 9.'/)
      read (*,100,err=2) nc
      if (nc .lt. 9  .or.  nc .gt. nwmax) go to 2
      write ( *, * ) '  nc = ', nc

    3 write (*,150)
  150 format (///5X,'Specify the number of nodes NW for ', &
              'the weights.  NW = 30 is'/5X,'recommended. ', &
              ' 1  LE  NW  LE  MIN(40,N-1).')
      read (*,100,err=2) nw
      if (1 .gt. nw  .or.  nw .gt. nwmax) go to 3
      write ( *, * ) '  nw = ', nw
!
    4 write (*,160) nrmax
  160 format (///5X,'Specify the number of rows and column', &
              's NR in the uniform grid'/5X,'of cells used', &
              ' to locate nearest neighbors.  NR = Sqrt(N/', &
              '3) is'/5X,'recommended.  1 LE NR LE ',I2)
      read (*,100,err=3) nr
      if (nr .lt. 1  .or.  nr .gt. nrmax) go to 4
      write ( *, * ) '  nr = ', nr
!
!  Set up uniform grid points.
!
      do i = 1,ni
        p(i) = dble ( i - 1 ) / dble ( ni - 1 )
      end do
!
!  Initialize the average SSE/SSM value to zero.
!
      ssa = 0.0D+00
!
!  Print a heading and loop on test functions.
!
      write (*,200) ks, n, ni, nc, nw, nr

      do kf = kff, kfl
!
!  Compute true function values at the nodes.
!
        do k = 1, n
          call tstfn1 (kf,x(k),y(k),0, w(k),dum,dum)
        end do
!
!  Compute true function values FT on the uniform grid, and
!  accumulate the sum of values SUM and sum of squared
!  values SSM.
!
        sum = 0.0D+00
        ssm = 0.0D+00
        do i = 1,ni
          do j = 1,ni
            call tstfn1 (kf,p(i),p(j),0, ft(i,j),dum,dum)
            sum = sum + ft(i,j)
            ssm = ssm + ft(i,j)**2
          end do
        end do
!
!  Compute the sum of squared deviations from the mean SSM.
!
        ssm = ssm - sum * sum / dble ( ni * ni )
!
!  Compute parameters A and RW defining the interpolant.
!
        call cshep2 (n,x,y,w,nc,nw,nr, lcell,lnext,xmin, &
                     ymin,dx,dy,rmax,rw,a,ier)

        if ( ier .ne. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'TOMS792_TEST - Fatal error!'
          write ( *, '(a,i4)' ) '  CSHEP2 returns IER = ', ier
          IF ( ier == 2 ) then
            write ( *, '(a)' ) '  Duplicate nodes encountered.'
          end if
          if ( ier == 3 ) then
            write ( *, '(a)' ) '  All nodes are collinear.'
          end if
          stop
        end if
!
!  Compute interpolation errors.
!
        ermean = 0.0D+00
        ermax = 0.0D+00
        sse = 0.0D+00

        do i = 1,ni
          do j = 1,ni
            pw = cs2val (p(i),p(j),n,x,y,w,nr,lcell,lnext, &
                         xmin,ymin,dx,dy,rmax,rw,a) - &
                 ft(i,j)
            ermean = ermean + abs(pw)
            ermax = max(ermax,abs(pw))
            sse = sse + pw*pw
          end do
        end do

        np = ni*ni
        ermean = ermean/dble(np)
        sse = sse/ssm
        ssa = ssa + sse
        write (*,210) kf, ermax, ermean, sse

      end do
!
!  Print the average SSE/SSM value (averaged over the test
!  functions).
!
      write (*,220) ssa/dble(nfun2)
      go to 1
!
!  Terminate.
!
999   continue
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS792_TEST():'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ''
      call timestamp ( )
      stop 0
!
! Print formats:
!
  200 format (30X,'CS2TST Output:'// &
              1X,16X,'CSHEP2D Interpolation Errors for ', &
              'N nodes and an'/ &
              1X,6X,'NI by NI Uniform Grid of Interpolation', &
              'ion Points in the Unit Square'//1X, &
              6X,'Node set ',I2,4X,'N =',I4,4X,'NI = ',I2, &
              4X,'NC = ',I2,4X,'NW = ',I2,4X,'NR = ',I2/// &
              1X,16X,'Function',4X,'Max Error',4X, &
              'Mean Error',4X,'SSE/SSM'/)
  210 format (1X,19X,I2,9X,F7.4,6X,F8.5,2X,F9.6)
  220 format (//1X,11X,'Average SSE/SSM over the test ', &
              'functions = ',F9.6)
!
! Error message formats:
!
  500 format (///1X,10X,'*** Error in data -- N = ',I4, &
              ', Maximum value =',I4,' ***')
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
!    15 August 2021
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
