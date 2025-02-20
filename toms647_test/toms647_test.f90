program main

!*****************************************************************************80
!
!! toms647_test() tests toms647().
!
!  Modified:
!
!    26 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer atmost
  integer dimen
  integer seed

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'toms647_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test TOMS647().'

  atmost = 7000
  dimen = 5
  seed = 12345

  call testf ( dimen, atmost )

  call testh ( dimen, atmost )

  call tests ( dimen, atmost )

  call testu ( dimen, atmost, seed )

  atmost = 7000
  dimen = 10
  seed = 12345

  call testf ( dimen, atmost )

  call testh ( dimen, atmost )

  call tests ( dimen, atmost )

  call testu ( dimen, atmost, seed )

  atmost = 7000
  dimen = 20
  seed = 12345

  call testf ( dimen, atmost )

  call testh ( dimen, atmost )

  call tests ( dimen, atmost )

  call testu ( dimen, atmost, seed )

  atmost = 100
  dimen = 2
  seed = 12345

  call showf ( dimen, atmost )

  call showh ( dimen, atmost )

  call shows ( dimen, atmost )

  call showu ( dimen, atmost, seed )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS647_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine testf ( dimen, atmost )

!*****************************************************************************80
!
!! TESTF uses "GOFAUR" to approximate an integral.
!
!  Discussion:
!
!    This routine tests the accuracy of numerical integration using GOFAUR
!    and integrand #2 of Davis and Rabinowitz, page 406.
!
!  Modified:
!
!    26 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer dimen

  integer atmost
  real correct
  real estimate
  real error
  real f
  logical flag(2)
  integer i
  integer j
  real quasi(dimen)
  real t1
  real t2
  real total

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TESTF'
  write ( *, '(a)' ) '  Test GOFAUR, for quasirandom number generation.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Spatial dimension = ', dimen
  write ( *, '(a,i6)' ) '  ATMOST = ', atmost

  call cpu_time ( t1 )

  call infaur ( flag, dimen, atmost )

  if ( .not. flag(1) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TESTF'
    write ( *, '(a)' ) '  Spatial dimension is not acceptable.'
    return
  end if

  if ( .not. flag(2) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TESTF'
    write ( *, '(a)' ) '  ATMOST is not acceptable.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  INFAUR has initialized the data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Step  Estimate      Error         Time'
  write ( *, '(a)' ) ' '

  correct = 1.0E+00
  total = 0.0E+00

  do i = 1, atmost

    call gofaur ( quasi )

    f = 1.0E+00
    do j = 1, dimen
      f = f * abs ( 4.0E+00 * quasi(j) - 2.0E+00 )
    end do

    total = total + f

    if ( mod ( i, 500 ) == 0 ) then
      call cpu_time ( t2 )
      estimate = total / real ( i )
      error = abs ( estimate - correct )
      write ( *, '(i6,3g14.6)' ) i, estimate, error, t2 - t1
    end if

  end do

  return
end
subroutine testh ( dimen, atmost )

!*****************************************************************************80
!
!! TESTH uses GOHALT to approximate an integral.
!
!  Discussion:
!
!    This routine tests the accuracy of numerical integration using GOHALT
!    and integrand #2 of Davis and Rabinowitz, page 406.
!
!  Modified:
!
!    14 March 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dimen

  integer atmost
  real ( kind = rk ) correct
  real ( kind = rk ) estimate
  real ( kind = rk ) error
  real ( kind = rk ) f
  logical flag(2)
  integer i
  integer j
  real ( kind = rk ) quasi(dimen)
  real ( kind = rk ) t1
  real ( kind = rk ) t2
  real ( kind = rk ) total

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TESTH'
  write ( *, '(a)' ) '  Test GOHALT, for quasirandom number generation.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Spatial dimension = ', dimen
  write ( *, '(a,i6)' ) '  ATMOST = ', atmost

  call cpu_time ( t1 )

  call inhalt ( flag, dimen, atmost, quasi )

  if ( .not. flag(1) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TESTH'
    write ( *, '(a)' ) '  The spatial dimension is not acceptable.'
    return
  end if

  if ( .not. flag(2) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TESTH'
    write ( *, '(a)' ) '  ATMOST is not acceptable.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  INHALT has initialized the data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Step  Estimate      Error         Time'
  write ( *, '(a)' ) ' '

  correct = 1.0D+00
  total = 0.0D+00

  do i = 1, atmost

    call gohalt ( quasi )

    f = 1.0D+00
    do j = 1, dimen
       f = f * abs ( 4.0D+00 * quasi(j) - 2.0D+00 )
    end do

    total = total + f

    if ( mod ( i, 500 ) == 0 ) then
      call cpu_time ( t2 )
      estimate = total / real ( i, kind = rk )
      error = abs ( estimate - correct )
      write ( *, '(i6,3g14.6)' ) i, estimate, error, t2 - t1
    end if

  end do

  return
end
subroutine tests ( dimen, atmost )

!*****************************************************************************80
!
!! TESTS uses "GOSOBL" to approximate an integral.
!
!  Discussion:
!
!    This routine tests the accuracy of numerical integration using GOSOBL
!    and integrand #2 of Davis and Rabinowitz, page 406.
!
!  Modified:
!
!    26 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer dimen

  integer atmost
  real correct
  real estimate
  real error
  real f
  logical flag(2)
  integer i
  integer j
  real quasi(dimen)
  real t1
  real t2
  integer taus
  real total

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TESTS'
  write ( *, '(a)' ) '  Test GOSOBL, for quasirandom number generation.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Spatial dimension = ', dimen
  write ( *, '(a,i6)' ) '  ATMOST = ', atmost

  call cpu_time ( t1 )

  call insobl ( flag, dimen, atmost, taus )

  if ( .not. flag(1) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TESTS'
    write ( *, '(a)' ) '  Spatial dimension is not acceptable.'
    return
  end if

  if ( .not. flag(2) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TESTS'
    write ( *, '(a)' ) '  ATMOST is not acceptable.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  INSOBL has initialized the data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Step  Estimate      Error         Time'
  write ( *, '(a)' ) ' '

  correct = 1.0E+00
  total = 0.0E+00

  do i = 1, atmost

    call gosobl ( quasi )

    f = 1.0E+00
    do j = 1, dimen
      f = f * abs ( 4.0E+00 * quasi(j) - 2.0E+00 )
    end do

    total = total + f

    if ( mod ( i, 500 ) == 0 ) then
      call cpu_time ( t2 )
      estimate = total / real ( i )
      error = abs ( estimate - correct )
      write ( *, '(i6,3g14.6)' ) i, estimate, error, t2 - t1
    end if

  end do

  return
end
subroutine testu ( dimen, atmost, seed )

!*****************************************************************************80
!
!! TESTU uses "UNIF" to approximate an integral.
!
!  Discussion:
!
!    This routine tests the accuracy of numerical integration using UNIF
!    and integrand #2 of Davis and Rabinowitz, page 406.
!
!  Modified:
!
!    26 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer dimen

  integer atmost
  real correct
  real estimate
  real error
  real f
  integer i
  integer k
  real quasi(dimen)
  integer seed
  real t1
  real t2
  real total
  real unif

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TESTU'
  write ( *, '(a)' ) '  Test UNIF, for pseudorandom number generation.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Spatial dimension = ', dimen
  write ( *, '(a,i6)' ) '  ATMOST = ', atmost
  write ( *, '(a,i12)' ) '  SEED = ', seed
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Step  Estimate      Error         Time'
  write ( *, '(a)' ) ' '

  call cpu_time ( t1 )

  correct = 1.0E+00
  total = 0.0E+00

  do i = 1, atmost

    f = 1.0E+00

    do k = 1, dimen
      quasi(k) = unif ( seed )
      f = f * abs ( 4.0E+00 * quasi(k) - 2.0E+00 )
    end do

    total = total + f

    if ( mod ( i, 500 ) == 0 ) then
      call cpu_time ( t2 )
      estimate = total / real ( i )
      error = abs ( estimate - correct )
      write ( *, '(i6,3g14.6)' ) i, estimate, error, t2 - t1
    end if

  end do

  return
end
subroutine showf ( dimen, atmost )

!*****************************************************************************80
!
!! SHOWF uses "GOFAUR" to compute some Faure points and write them out.
!
!  Modified:
!
!    26 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer dimen

  integer atmost
  logical flag(2)
  integer i
  real quasi(dimen)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SHOWF'
  write ( *, '(a)' ) '  Test GOFAUR, for quasirandom number generation.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Spatial dimension = ', dimen
  write ( *, '(a,i6)' ) '  ATMOST = ', atmost

  call infaur ( flag, dimen, atmost )

  if ( .not. flag(1) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SHOWF'
    write ( *, '(a)' ) '  Spatial dimension is not acceptable.'
    return
  end if

  if ( .not. flag(2) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SHOWF'
    write ( *, '(a)' ) '  ATMOST is not acceptable.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  INFAUR has initialized the data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Step  QUASI(1)   QUASI(2)'
  write ( *, '(a)' ) ' '

  do i = 1, atmost

    call gofaur ( quasi )
    write ( *, '(i4,2f7.4)' ) i, quasi(1:2)

  end do

  return
end
subroutine showh ( dimen, atmost )

!*****************************************************************************80
!
!! SHOWH uses "GOHALT" to compute some Halton points and write them out.
!
!  Modified:
!
!    26 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dimen

  integer atmost
  logical flag(2)
  integer i
  real ( kind = rk ) quasi(dimen)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SHOWH'
  write ( *, '(a)' ) '  Test GOHALT, for quasirandom number generation.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Spatial dimension = ', dimen
  write ( *, '(a,i6)' ) '  ATMOST = ', atmost

  call inhalt ( flag, dimen, atmost, quasi )

  if ( .not. flag(1) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SHOWH'
    write ( *, '(a)' ) '  The spatial dimension is not acceptable.'
    return
  end if

  if ( .not. flag(2) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SHOWH'
    write ( *, '(a)' ) '  ATMOST is not acceptable.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  INHALT has initialized the data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Step  QUASI(1)   QUASI(2)'
  write ( *, '(a)' ) ' '

  write ( *, '(i4,2f7.4)' ) 0, quasi(1:2)

  do i = 1, atmost-1

    call gohalt ( quasi )

    write ( *, '(i4,2f7.4)' ) i, quasi(1:2)

  end do

  return
end
subroutine shows ( dimen, atmost )

!*****************************************************************************80
!
!! SHOWS uses "GOSOBL" to compute some Sobol points and write them out.
!
!  Modified:
!
!    26 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer dimen

  integer atmost
  logical flag(2)
  integer i
  real quasi(dimen)
  integer taus

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SHOWS'
  write ( *, '(a)' ) '  Test GOSOBL, for quasirandom number generation.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Spatial dimension = ', dimen
  write ( *, '(a,i6)' ) '  ATMOST = ', atmost

  call insobl ( flag, dimen, atmost, taus )

  if ( .not. flag(1) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TESTS'
    write ( *, '(a)' ) '  Spatial dimension is not acceptable.'
    return
  end if

  if ( .not. flag(2) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TESTH'
    write ( *, '(a)' ) '  ATMOST is not acceptable.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  INSOBL has initialized the data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Step  QUASI(1)   QUASI(2)'
  write ( *, '(a)' ) ' '

  do i = 1, atmost

    call gosobl ( quasi )

    write ( *, '(i4,2f7.4)' ) i, quasi(1:2)

  end do

  return
end
subroutine showu ( dimen, atmost, seed )

!*****************************************************************************80
!
!! SHOWU uses "UNIF" to compute some uniform pseudorandom numbers.
!
!  Modified:
!
!    26 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer dimen

  integer atmost
  integer i
  integer j
  real quasi(dimen)
  integer seed
  real unif

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SHOWU'
  write ( *, '(a)' ) '  Test UNIF, for pseudorandom number generation.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Spatial dimension = ', dimen
  write ( *, '(a,i6)' ) '  ATMOST = ', atmost
  write ( *, '(a,i12)' ) '  SEED = ', seed
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Step  QUASI(1)   QUASI(2)'
  write ( *, '(a)' ) ' '

  do i = 1, atmost

    do j = 1, dimen
      quasi(j) = unif ( seed )
    end do

    write ( *, '(i4,2f7.4)' ) i, quasi(1:2)

  end do

  return
end
