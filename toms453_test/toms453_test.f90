program main

!*****************************************************************************80
!
!! toms453_test() tests toms453().
!
!  Modified:
!
!    12 January 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS453_TEST():'
  write ( *, '(a)' ) '  FORTRAN80 version'
  write ( *, '(a)' ) '  Test TOMS453().'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS453_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests BROMIN, ACM TOMS algorithm 453.
!
!  Modified:
!
!    11 July 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nhalf_max
  integer n_num
  integer s_num

  parameter ( nhalf_max = 6 )
  parameter ( n_num = 3 )
  parameter ( s_num = 4 )

  real ( kind = rk ) eps
  integer i
  integer ier
  integer j
  integer k
  integer kk
  integer n
  integer n_half
  integer n_vec(n_num)
  real ( kind = rk ) s
  real ( kind = rk ) s_vec(s_num)
  real ( kind = rk ) tol
  real ( kind = rk ) total
  real ( kind = rk ) wi(nhalf_max)
  real ( kind = rk ) wr(nhalf_max)
  real ( kind = rk ) xi(nhalf_max)
  real ( kind = rk ) xr(nhalf_max)

  save n_vec
  save s_vec

  data n_vec / 6, 9, 12 /
  data s_vec / 0.0D+00, 0.1D+00, 1.0D+00, 4.0D+00 /

  tol = 0.1D-08

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Determine abscissas and weights for'
  write ( *, '(a)' ) '  a variety of values of S and N.'

  do i = 1, n_num

    n = n_vec(i)
    n_half = ( n + 1 ) / 2

    do j = 1, s_num

      s = s_vec(j)

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  N = ', n
      write ( *, '(a,g14.6)' ) '  S = ', s

      call bromin ( n, s, tol, xr, xi, wr, wi, eps, ier )

      if ( 0 < ier ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)' ) 'BROMIN returned IER = ', ier

      else

        if ( ier == -1 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Note that the requested accuracy'
          write ( *, '(a)' ) '  was not achieved.'
        end if

        write ( *, '(a)' ) ' '
        write ( *, '(a,a)' ) '                           ', &
          'XR              XI              WR              WI'
        write ( *, '(a)' ) ' '
        total = 0.0D+00
        do kk = 1, n
          if ( kk .le. ( n - n_half ) ) then
            k = n_half + 1 - kk
            write ( *, '(2x,i8,2x,i8,4(2x,g14.6))' ) &
              kk, k, xr(k), - xi(k), wr(k), - wi(k)
            total = total + wr(k)
          else
            k = kk - ( n - n_half )
            write ( *, '(2x,i8,2x,i8,4(2x,g14.6))' ) &
              kk, k, xr(k),   xi(k), wr(k),   wi(k)
            total = total + wr(k)
          end if
        end do
        write ( *, '(2x,a8,2x,8x,16x,16x,2x,g14.6)' ) 'WR total', total

      end if

    end do

  end do

  return
end
