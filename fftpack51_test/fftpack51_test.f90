program main

!*****************************************************************************80
!
!! fftpack51_test() tests fftpack51().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FFTPACK51_TEST():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test FFTPACK51.'

  call cfft1_test ( )
  call cfft2_test ( )
  call cfftm_test ( )
  call cosq1_test ( )
  call cosqm_test ( )
  call cost1_test ( )
  call costm_test ( )
  call rfft1_test ( )
  call rfft2_test ( )
  call rfftm_test ( )
  call sinq1_test ( )
  call sinqm_test ( )
  call sint1_test ( )
  call sintm_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FFTPACK51_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine cfft1_test ( )

!*****************************************************************************80
!
!! cfft1_test tests CFFT1B, CFFT1F and CFFT1I.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 4096

  complex ( kind = ck ) c(n)
  integer ier
  integer inc
  integer lenc
  integer lensav
  integer lenwrk
  integer seed
  real ( kind = rk ), allocatable, dimension ( : ) :: work
  real ( kind = rk ), allocatable, dimension ( : ) :: wsave

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'cfft1_test'
  write ( *, '(a)' ) '  For complex double precision fast Fourier transforms, 1D,'
  write ( *, '(a)' ) '  CFFT1I initializes the transform,'
  write ( *, '(a)' ) '  CFFT1F does a forward transform;'
  write ( *, '(a)' ) '  CFFT1B does a backward transform.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of data items is N = ', n
!
!  Allocate the work arrays.
!
  lenwrk = 2 * n
  lensav = 2 * n + int ( log ( real ( n, kind = rk ) ) / log ( 2.0D+00 ) ) + 4

  write ( *, '(a,i8)' ) '  LENSAV = ', lensav
  write ( *, '(a,i8)' ) '  LENWRK = ', lenwrk

  allocate ( work(1:lenwrk) )
  allocate ( wsave(1:lensav) )

  call cfft1i ( n, wsave, lensav, ier )
!
!  Set the data values.
!
  seed = 1973

  call c8vec_uniform_01 ( n, seed, c )

  call c8vec_print_part ( n, c, 10, '  The original data:' ) 
!
!  Compute the FFT coefficients.
!
  inc = 1
  lenc = n

  call cfft1f ( n, inc, c, lenc, wsave, lensav, work, lenwrk, ier )

  call c8vec_print_part ( n, c, 10, '  The FFT coefficients:' )
!
!  Compute inverse FFT of coefficients.  Should get back the
!  original data.
!
  call cfft1b ( n, inc, c, lenc, wsave, lensav, work, lenwrk, ier )

  call c8vec_print_part ( n, c, 10, '  The retrieved data:' )

  deallocate ( work )
  deallocate ( wsave )

  return
end
subroutine cfft2_test ( )

!*****************************************************************************80
!
!! cfft2_test tests CFFT2B, CFFT2F and CFFT2I.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: l = 32
  integer, parameter :: m = 64

  complex ( kind = ck ) c(l,m)
  integer ier
  integer ldim
  integer lensav
  integer lenwrk
  integer seed
  real ( kind = rk ), allocatable, dimension ( : ) :: work
  real ( kind = rk ), allocatable, dimension ( : ) :: wsave

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'cfft2_test'
  write ( *, '(a)' ) '  For complex double precision fast Fourier transforms, 2D,'
  write ( *, '(a)' ) '  CFFT2I initializes the transform,'
  write ( *, '(a)' ) '  CFFT2F does a forward transform;'
  write ( *, '(a)' ) '  CFFT2B does a backward transform.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The data is stored in an L by M array, with'
  write ( *, '(a,i8)' ) '  L = ', l
  write ( *, '(a,i8)' ) '  M = ', m
!
!  Allocate work arrays.
!
  lenwrk = 2 * l * m

  lensav = 2 * l + int ( log ( real ( l, kind = rk ) ) / log ( 2.0D+00 ) ) &
    + 2 * m + int ( log ( real ( m, kind = rk ) ) / log ( 2.0D+00 ) ) &
    + 8 

  write ( *, '(a,i8)' ) '  LENSAV = ', lensav
  write ( *, '(a,i8)' ) '  LENWRK = ', lenwrk

  allocate ( work(1:lenwrk) )
  allocate ( wsave(1:lensav) )

  call cfft2i ( l, m, wsave, lensav, ier )
!
!  Set the data values.
!
  seed = 1973

  call c8mat_uniform_01 ( l, m, seed, c )

  call c8mat_print_some ( l, m, c, 1, 1, 5, 5, &
    '  Part of the original data:' )
!
!  Compute the FFT coefficients.
!
  ldim = l

  call cfft2f ( ldim, l, m, c, wsave, lensav, work, lenwrk, ier )

  call c8mat_print_some ( l, m, c, 1, 1, 5, 5, &
    '  Part of the FFT coefficients:' )
!
!  Compute inverse FFT of coefficients.  Should get back the
!  original data.
!
  call cfft2b ( ldim, l, m, c, wsave, lensav, work, lenwrk, ier )

  call c8mat_print_some ( l, m, c, 1, 1, 5, 5, '  Part of the retrieved data:' )

  deallocate ( work )
  deallocate ( wsave )

  return
end
subroutine cfftm_test ( )

!*****************************************************************************80
!
!! cfftm_test tests CFFTMB, CFFTMF and CFFTMI.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 32
  integer, parameter :: lot = 6

  complex ( kind = ck ), allocatable, dimension ( : ) :: c
  integer ier
  integer inc
  integer jump
  integer lenc
  integer lensav
  integer lenwrk
  integer seed
  real ( kind = rk ), allocatable, dimension ( : ) :: work
  real ( kind = rk ), allocatable, dimension ( : ) :: wsave

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'cfftm_test'
  write ( *, '(a)' ) '  For complex double precision fast Fourier transforms, 1D, multiple'
  write ( *, '(a)' ) '  CFFTMI initializes the transform,'
  write ( *, '(a)' ) '  CFFTMF does a forward transform;'
  write ( *, '(a)' ) '  CFFTMB does a backward transform.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of sequences is LOT = ', lot
  write ( *, '(a,i8)' ) '  The length of each sequence is N = ', n
!
!  Set work vectors.
!
  lenc = n * lot
  lenwrk = 2 * lot * n
  lensav = 2 * n + int ( log ( real ( n, kind = rk ) ) / log ( 2.0D+00 ) ) + 4

  write ( *, '(a,i8)' ) '  LENC   = ', lenc
  write ( *, '(a,i8)' ) '  LENSAV = ', lensav
  write ( *, '(a,i8)' ) '  LENWRK = ', lenwrk

  allocate ( c(1:lenc) )
  allocate ( wsave(1:lensav) )
  allocate ( work(1:lenwrk) )

  call cfftmi ( n, wsave, lensav, ier )
!
!  Set the data values.
!
  seed = 1973

  call c8mat_uniform_01 ( n, lot, seed, c )

  call c8mat_print_some ( n, lot, c, 1, 1, 5, 5, &
    '  Part of the original data:' )
!
!  Compute the FFT coefficients.
!
  jump = n
  inc = 1

  call cfftmf ( lot, jump, n, inc, c, lenc, wsave, lensav, work, lenwrk, ier )

  call c8mat_print_some ( n, lot, c, 1, 1, 5, 5, &
    '  Part of the FFT coefficients:' )
!
!  Compute inverse FFT of coefficients.  Should get back the
!  original data.
!
  call cfftmb ( lot, jump, n, inc, c, lenc, wsave, lensav, work, lenwrk, ier )

  call c8mat_print_some ( n, lot, c, 1, 1, 5, 5, &
    '  Part of the retrieved data:' )

  deallocate ( c )
  deallocate ( wsave )
  deallocate ( work )

  return
end
subroutine cosq1_test ( )

!*****************************************************************************80
!
!! cosq1_test tests COSQ1B, COSQ1F and COSQ1I.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 4096

  integer ier
  integer inc
  integer lenr
  integer lensav
  integer lenwrk
  real ( kind = rk ) r(n)
  integer seed
  real ( kind = rk ), allocatable, dimension ( : ) :: work
  real ( kind = rk ), allocatable, dimension ( : ) :: wsave

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'cosq1_test'
  write ( *, '(a)' ) '  For real double precision fast cosine transforms, 1D,'
  write ( *, '(a)' ) '  COSQ1I initializes the transform,'
  write ( *, '(a)' ) '  COSQ1F does a forward transform;'
  write ( *, '(a)' ) '  COSQ1B does a backward transform.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of data items is N = ', n
!
!  Set work vectors.
!
  lensav = 2 * n + int ( log ( real ( n, kind = rk ) ) / log ( 2.0D+00 ) ) + 4
  lenwrk = n

  write ( *, '(a,i8)' ) '  LENSAV = ', lensav
  write ( *, '(a,i8)' ) '  LENWRK = ', lenwrk

  allocate ( work(1:lenwrk) )
  allocate ( wsave(1:lensav) )

  call cosq1i ( n, wsave, lensav, ier )
!
!  Set the data values.
!
  seed = 1973

  call random_number ( harvest = r(1:n) )

  call r8vec_print_part ( n, r, 10, '  The original data:' )
!
!  Compute the FFT coefficients.
!
  inc = 1
  lenr = n

  call cosq1f ( n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

  call r8vec_print_part ( n, r, 10, '  The FFT coefficients:' )
!
!  Compute inverse FFT of coefficients.  Should get back the
!  original data.
!
  call cosq1b ( n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

  call r8vec_print_part ( n, r, 10, '  The retrieved data:' )

  deallocate ( work )
  deallocate ( wsave )

  return
end
subroutine cosqm_test ( )

!*****************************************************************************80
!
!! cosqm_test tests COSQMB, COSQMF and COSQMI.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 August 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 32
  integer, parameter :: lot = 6    

  integer ier
  integer inc
  integer jump
  integer lenr
  integer lensav
  integer lenwrk
  real ( kind = rk ), allocatable, dimension ( :, : ) :: r
  integer seed
  real ( kind = rk ), allocatable, dimension ( : ) :: work
  real ( kind = rk ), allocatable, dimension ( : ) :: wsave

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'cosqm_test'
  write ( *, '(a)' ) '  For real double precision fast cosine transform, 1D, multiple'
  write ( *, '(a)' ) '  COSQMI initializes the transform,'
  write ( *, '(a)' ) '  COSQMF does a forward transform;'
  write ( *, '(a)' ) '  COSQMB does a backward transform.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of sequences is LOT =   ', lot
  write ( *, '(a,i8)' ) '  The length of each sequence is N = ', n
!
!  Set work arrays.
!
  lenr = n * lot
  lenwrk = lot * n
  lensav = 2 * n + int ( log ( real ( n, kind = rk ) ) / log ( 2.0D+00 ) ) + 4

  write ( *, '(a,i8)' ) '  LENR   = ', lenr
  write ( *, '(a,i8)' ) '  LENSAV = ', lensav
  write ( *, '(a,i8)' ) '  LENWRK = ', lenwrk

  allocate ( r(1:n,1:lot) )
  allocate ( work(1:lenwrk) )
  allocate ( wsave(1:lensav) )

  call cosqmi ( n, wsave, lensav, ier )
!
!  Set the data values.
!
  seed = 1973

  call random_number ( harvest = r(1:n,1:lot) )

  call r8mat_print_some ( n, lot, r, 1, 1, 5, 5, &
    '  Part of the original data:' )
!
!  Compute the FFT coefficients.
!
  jump = n
  inc = 1

  call cosqmf ( lot, jump, n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

  call r8mat_print_some ( n, lot, r, 1, 1, 5, 5, &
    '  Part of the FFT coefficients:' )
!
!  Compute inverse FFT of coefficients.  Should get back the
!  original data.
!
  call cosqmb ( lot, jump, n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

  call r8mat_print_some ( n, lot, r, 1, 1, 5, 5, &
    '  Part of the retrieved data:' )

  deallocate ( r )
  deallocate ( work )
  deallocate ( wsave )

  return
end
subroutine cost1_test ( )

!*****************************************************************************80
!
!! cost1_test tests COST1B, COST1F and COST1I.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 4096 

  integer ier
  integer inc
  integer lenr
  integer lensav
  integer lenwrk
  real ( kind = rk ) r(n)
  integer seed
  real ( kind = rk ), allocatable, dimension ( : ) :: work
  real ( kind = rk ), allocatable, dimension ( : ) :: wsave

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'cost1_test'
  write ( *, '(a)' ) '  For real double precision fast cosine transforms, 1D,'
  write ( *, '(a)' ) '  COST1I initializes the transform,'
  write ( *, '(a)' ) '  COST1F does a forward transform;'
  write ( *, '(a)' ) '  COST1B does a backward transform.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of data items is N = ', n
!
!  Set work arrays.
!
  lenwrk = n - 1 
  lensav = 2 * n + int ( log ( real ( n, kind = rk ) ) / log ( 2.0D+00 ) ) + 4 

  write ( *, '(a,i8)' ) '  LENSAV = ', lensav
  write ( *, '(a,i8)' ) '  LENWRK = ', lenwrk

  allocate ( work(1:lenwrk) )
  allocate ( wsave(1:lensav) )

  call cost1i ( n, wsave, lensav, ier )
!
!  Set the data values.
!
  seed = 1973

  call random_number ( harvest = r(1:n) )

  call r8vec_print_part ( n, r, 10, '  The original data:' )
!
!  Compute the FFT coefficients.
!
  inc = 1
  lenr = n

  call cost1f ( n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

  call r8vec_print_part ( n, r, 10, '  The FFT coefficients:' )
!
!  Compute inverse FFT of coefficients.
!
  call cost1b ( n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

  call r8vec_print_part ( n, r, 10, '  The retrieved data:' )

  deallocate ( work )
  deallocate ( wsave )

  return
end
subroutine costm_test ( )

!*****************************************************************************80
!
!! costm_test tests COSTMB, COSTMF and COSTMI.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 32
  integer, parameter :: lot = 6

  integer ier
  integer inc
  integer jump
  integer lenr
  integer lensav
  integer lenwrk
  real ( kind = rk ), allocatable, dimension ( :, : ) :: r
  integer seed
  real ( kind = rk ), allocatable, dimension ( : ) :: work
  real ( kind = rk ), allocatable, dimension ( : ) :: wsave

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'costm_test'
  write ( *, '(a)' ) '  For real double precision fast cosine transforms, 1D, multiple'
  write ( *, '(a)' ) '  COSTMI initializes the transform,'
  write ( *, '(a)' ) '  COSTMF does a forward transform;'
  write ( *, '(a)' ) '  COSTMB does a backward transform.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of sequences is LOT =   ', lot
  write ( *, '(a,i8)' ) '  The length of each sequence is N = ', n
!
!  Set work arrays.
!
  lenr = n * lot
  lensav = 2 * n + int ( log ( real ( n, kind = rk ) ) / log ( 2.0D+00 ) ) + 4
  lenwrk = lot * ( n + 1 )

  write ( *, '(a,i8)' ) '  LENR   = ', lenr
  write ( *, '(a,i8)' ) '  LENSAV = ', lensav
  write ( *, '(a,i8)' ) '  LENWRK = ', lenwrk

  allocate ( r(1:n,1:lot) )
  allocate ( wsave(1:lensav) )
  allocate ( work(1:lenwrk) )

  call costmi ( n, wsave, lensav, ier )
!
!  Set the data values.
!
  seed = 1973

  call random_number ( harvest = r(1:n,1:lot) )

  call r8mat_print_some ( n, lot, r, 1, 1, 5, 5, &
    '  Part of the original data:' )
!
!  Compute the FFT coefficients.
!
  jump = n
  inc = 1

  call costmf ( lot, jump, n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

  call r8mat_print_some ( n, lot, r, 1, 1, 5, 5, &
    '  Part of the FFT coefficients:' )
!
!  Compute inverse FFT of coefficients.  Should get back the
!  original data.
!
  call costmb ( lot, jump, n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

  call r8mat_print_some ( n, lot, r, 1, 1, 5, 5, &
    '  Part of the retrieved data:' )

  deallocate ( r )
  deallocate ( work )
  deallocate ( wsave )

  return
end
subroutine rfft1_test ( )

!*****************************************************************************80
!
!! rfft1_test tests RFFT1B, RFFT1F and RFFT1I.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 4096

  integer ier
  integer inc
  integer lenr
  integer lensav
  integer lenwrk
  real ( kind = rk ) r(n)
  integer seed
  real ( kind = rk ), allocatable, dimension ( : ) :: work
  real ( kind = rk ), allocatable, dimension ( : ) :: wsave

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'rfft1_test'
  write ( *, '(a)' ) '  For real double precision fast cosine transforms, 1D,'
  write ( *, '(a)' ) '  RFFT1I initializes the transform,'
  write ( *, '(a)' ) '  RFFT1F does a forward transform;'
  write ( *, '(a)' ) '  RFFT1B does a backward transform.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of data items is N = ', n
!
!  Set work vectors.
!
  lensav = n + int ( log ( real ( n, kind = rk ) ) / log ( 2.0D+00 ) ) + 4
  lenwrk = n

  write ( *, '(a,i8)' ) '  LENSAV = ', lensav
  write ( *, '(a,i8)' ) '  LENWRK = ', lenwrk

  allocate ( work(1:lenwrk) )
  allocate ( wsave(1:lensav) )

  call rfft1i ( n, wsave, lensav, ier )
!
!  Set the data values.
!
  seed = 1973

  call random_number ( harvest = r(1:n) )

  call r8vec_print_part ( n, r, 10, '  The original data:' )
!
!  Compute the FFT coefficients.
!
  inc = 1
  lenr = n

  call rfft1f ( n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

  call r8vec_print_part ( n, r, 10, '  The FFT coefficients:' )
!
!  Compute inverse FFT of coefficients.  Should get back the
!  original data.
!
  call rfft1b ( n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

  call r8vec_print_part ( n, r, 10, '  The retrieved data:' )

  deallocate ( work )
  deallocate ( wsave )

  return
end
subroutine rfft2_test ( )

!*****************************************************************************80
!
!! rfft2_test tests RFFT2B, RFFT2F and RFFT2I.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: l = 32
  integer, parameter :: m = 64

  integer, parameter :: ldim = 2 * ( l / 2 + 1 )

  integer ier
  integer lensav
  integer lenwrk
  real ( kind = rk ) r(ldim,m)
  integer seed
  real ( kind = rk ), allocatable, dimension ( : ) :: work
  real ( kind = rk ), allocatable, dimension ( : ) :: wsave

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'rfft2_test'
  write ( *, '(a)' ) '  For real double precision fast Fourier transform, 2D,'
  write ( *, '(a)' ) '  RFFT2I initializes the transform,'
  write ( *, '(a)' ) '  RFFT2F does a forward transform;'
  write ( *, '(a)' ) '  RFFT2B does a backward transform.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The L by M data is stored in an LDIM by M array, with'
  write ( *, '(a,i8)' ) '  L =    ', l
  write ( *, '(a,i8)' ) '  LDIM = ', ldim
  write ( *, '(a,i8)' ) '  M =    ', m
!
!  Set work arrays.
!
  lenwrk = 2 * ldim * m

  lensav = &
          l + int ( log ( real ( l, kind = rk ) ) / log ( 2.0D+00 ) ) + 4 &
    + 2 * m + int ( log ( real ( m, kind = rk ) ) / log ( 2.0D+00 ) ) + 4 &
    +     m + int ( log ( real ( m, kind = rk ) ) / log ( 2.0D+00 ) ) + 4

  write ( *, '(a,i8)' ) '  LENSAV = ', lensav
  write ( *, '(a,i8)' ) '  LENWRK = ', lenwrk

  allocate ( work(1:lenwrk) )
  allocate ( wsave(1:lensav) )

  call rfft2i ( l, m, wsave, lensav, ier )
!
!  Set the data values.
!
  seed = 1973

  call random_number ( harvest = r(1:ldim,1:m) )

  call r8mat_print_some ( ldim, m, r, 1, 1, 5, 5, &
    '  Part of the original data:' )
!
!  Compute the FFT coefficients.
!
  call rfft2f ( ldim, l, m, r, wsave, lensav, work, lenwrk, ier )

  call r8mat_print_some ( ldim, m, r, 1, 1, 5, 5, &
    '  Part of the FFT coefficients:' )
!
!  Compute inverse FFT of coefficients.  Should get back the
!  original data.
!
  call rfft2b ( ldim, l, m, r, wsave, lensav, work, lenwrk, ier )

  call r8mat_print_some ( ldim, m, r, 1, 1, 5, 5, &
    '  Part of the retrieved data:' )

  deallocate ( work )
  deallocate ( wsave )

  return
end
subroutine rfftm_test ( )

!*****************************************************************************80
!
!! rfftm_test tests RFFTMB, RFFTMF and RFFTMI.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 32
  integer, parameter :: lot = 6

  integer ier
  integer inc
  integer jump
  integer lenr
  integer lensav
  integer lenwrk
  real ( kind = rk ), allocatable, dimension ( :, : ) :: r
  integer seed
  real ( kind = rk ), allocatable, dimension ( : ) :: work
  real ( kind = rk ), allocatable, dimension ( : ) :: wsave

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'rfftm_test'
  write ( *, '(a)' ) '  For real double precision fast Fourier transform, 1D, multiple'
  write ( *, '(a)' ) '  RFFTMI initializes the transform,'
  write ( *, '(a)' ) '  RFFTMF does a forward transform;'
  write ( *, '(a)' ) '  RFFTMB does a backward transform.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of sequences is LOT =   ', lot
  write ( *, '(a,i8)' ) '  The length of each sequence is N = ', n
!
!  Set work arrays.
!
  lenr = n * lot
  lenwrk = lot * n
  lensav = n + int ( log ( real ( n, kind = rk ) ) / log ( 2.0D+00 ) ) + 4

  write ( *, '(a,i8)' ) '  LENR =   ', lenr
  write ( *, '(a,i8)' ) '  LENSAV = ', lensav
  write ( *, '(a,i8)' ) '  LENWRK = ', lenwrk

  allocate ( r(1:n,1:lot) )
  allocate ( wsave(1:lensav) )
  allocate ( work(1:lenwrk) )

  call rfftmi ( n, wsave, lensav, ier )
!
!  Set the data values.
!
  seed = 1973

  call random_number ( harvest = r(1:n,1:lot) )

  call r8mat_print_some ( n, lot, r, 1, 1, 5, 5, '  Part of the original data:' )
!
!  Compute the FFT coefficients.
!
  jump = n
  inc = 1

  call rfftmf ( lot, jump, n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

  call r8mat_print_some ( n, lot, r, 1, 1, 5, 5, &
    '  Part of the FFT coefficients:' )
!
!  Compute inverse FFT of coefficients.  Should get back the
!  original data.
!
  call rfftmb ( lot, jump, n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

  call r8mat_print_some ( n, lot, r, 1, 1, 5, 5, &
    '  Part of the retrieved data:' )

  deallocate ( r )
  deallocate ( work )
  deallocate ( wsave )

  return
end
subroutine sinq1_test ( )

!*****************************************************************************80
!
!! sinq1_test tests SINQ1B, SINQ1F and SINQ1I.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 4096

  integer ier
  integer inc
  integer lenr
  integer lensav
  integer lenwrk
  real ( kind = rk ) r(n)
  integer seed
  real ( kind = rk ), allocatable, dimension ( : ) :: work
  real ( kind = rk ), allocatable, dimension ( : ) :: wsave

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'sinq1_test'
  write ( *, '(a)' ) '  For real double precision fast sine transforms, 1D,'
  write ( *, '(a)' ) '  SINQ1I initializes the transform,'
  write ( *, '(a)' ) '  SINQ1F does a forward transform;'
  write ( *, '(a)' ) '  SINQ1B does a backward transform.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of data items is N = ', n
!
!  Set work arrays.
!
  lenwrk = n
  lensav = 2 * n + int ( log ( real ( n, kind = rk ) ) / log ( 2.0D+00 ) ) + 4

  write ( *, '(a,i8)' ) '  LENSAV = ', lensav
  write ( *, '(a,i8)' ) '  LENWRK = ', lenwrk

  allocate ( work(1:lenwrk) )
  allocate ( wsave(1:lensav) )

  call sinq1i ( n, wsave, lensav, ier )
!
!  Set the data values.
!
  seed = 1973

  call random_number ( harvest = r(1:n) )

  call r8vec_print_part ( n, r, 10, '  The original data:' )
!
!  Compute the FFT coefficients.
!
  inc = 1
  lenr = n

  call sinq1f ( n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

  call r8vec_print_part ( n, r, 10, '  The FFT coefficients:' )
!
!  Compute inverse FFT of coefficients.  Should get back the
!  original data.
!
  call sinq1b ( n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

  call r8vec_print_part ( n, r, 10, '  The retrieved data:' )

  deallocate ( work )
  deallocate ( wsave )

  return
end
subroutine sinqm_test ( )

!*****************************************************************************80
!
!! sinqm_test tests SINQMB, SINQMF and SINQMI.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 32
  integer, parameter :: lot = 6

  integer ier
  integer inc
  integer jump
  integer lenr
  integer lensav
  integer lenwrk
  real ( kind = rk ), allocatable, dimension ( :, : ) :: r
  integer seed
  real ( kind = rk ), allocatable, dimension ( : ) :: work
  real ( kind = rk ), allocatable, dimension ( : ) :: wsave

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'sinqm_test'
  write ( *, '(a)' ) '  For real double precision fast sine transforms, 1D, multiple'
  write ( *, '(a)' ) '  SINQMI initializes the transform,'
  write ( *, '(a)' ) '  SINQMF does a forward transform;'
  write ( *, '(a)' ) '  SINQMB does a backward transform.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of sequences is LOT =   ', lot
  write ( *, '(a,i8)' ) '  The length of each sequence is N = ', n
!
!  Set work arrays.
!
  lenr = n * lot
  lenwrk = lot * n
  lensav = 2 * n + int ( log ( real ( n, kind = rk ) ) / log ( 2.0D+00 ) ) + 4

  write ( *, '(a,i8)' ) '  LENR   = ', lenr
  write ( *, '(a,i8)' ) '  LENSAV = ', lensav
  write ( *, '(a,i8)' ) '  LENWRK = ', lenwrk

  allocate ( r(1:n,1:lot) )
  allocate ( wsave(1:lensav) )
  allocate ( work(1:lenwrk) )

  call sinqmi ( n, wsave, lensav, ier )
!
!  Set the data values.
!
  seed = 1973

  call random_number ( harvest = r(1:n,1:lot) )

  call r8mat_print_some ( n, lot, r, 1, 1, 5, 5, &
    '  Part of the original data:' )
!
!  Compute the FFT coefficients.
!
  jump = n
  inc = 1

  call sinqmf ( lot, jump, n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

  call r8mat_print_some ( n, lot, r, 1, 1, 5, 5, &
    '  Part of the FFT coefficients:' )
!
!  Compute inverse FFT of coefficients.  Should get back the
!  original data.
!
  call sinqmb ( lot, jump, n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

  call r8mat_print_some ( n, lot, r, 1, 1, 5, 5, &
    '  Part of the retrieved data:' )

  deallocate ( r )
  deallocate ( work )
  deallocate ( wsave )

  return
end
subroutine sint1_test ( )

!*****************************************************************************80
!
!! sint1_test tests SINT1B, SINT1F and SINT1I.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 4096

  integer ier
  integer inc
  integer lenr
  integer lensav
  integer lenwrk
  real ( kind = rk ) r(n)
  integer seed
  real ( kind = rk ), allocatable, dimension ( : ) :: work
  real ( kind = rk ), allocatable, dimension ( : ) :: wsave

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'sint1_test'
  write ( *, '(a)' ) '  For real double precision fast sine transforms, 1D,'
  write ( *, '(a)' ) '  SINT1I initializes the transform,'
  write ( *, '(a)' ) '  SINT1F does a forward transform;'
  write ( *, '(a)' ) '  SINT1B does a backward transform.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of data items is N = ', n
!
!  Set work arrays.
!
  lenwrk = 2 * ( n + 1 )
  lensav = n / 2 + n + int ( log ( real ( n, kind = rk ) ) / log ( 2.0D+00 ) ) + 4

  write ( *, '(a,i8)' ) '  LENSAV = ', lensav
  write ( *, '(a,i8)' ) '  LENWRK = ', lenwrk

  allocate ( wsave(1:lensav) )
  allocate ( work(1:lenwrk) )

  call sint1i ( n, wsave, lensav, ier )
!
!  Set the data values.
!
  seed = 1973

  call random_number ( harvest = r(1:n) )

  call r8vec_print_part ( n, r, 10, '  The original data:' )
!
!  Compute the FFT coefficients.
!
  inc = 1
  lenr = n

  call sint1f ( n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

  call r8vec_print_part ( n, r, 10, '  The FFT coefficients:' )
!
!  Compute inverse FFT of coefficients.  Should get back the
!  original data.
!
  call sint1b ( n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

  call r8vec_print_part ( n, r, 10, '  The retrieved data:' )

  deallocate ( work )
  deallocate ( wsave )

  return
end
subroutine sintm_test ( )

!*****************************************************************************80
!
!! sintm_test tests SINTMB, SINTMF and SINTMI.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 32
  integer, parameter :: lot = 6

  integer ier
  integer inc
  integer jump
  integer lenr
  integer lensav
  integer lenwrk
  real ( kind = rk ), allocatable, dimension ( :, : ) :: r
  integer seed
  real ( kind = rk ), allocatable, dimension ( : ) :: work
  real ( kind = rk ), allocatable, dimension ( : ) :: wsave

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'sintm_test'
  write ( *, '(a)' ) '  For real double precision fast sine transforms, 1D, multiple'
  write ( *, '(a)' ) '  SINTMI initializes the transform,'
  write ( *, '(a)' ) '  SINTMF does a forward transform;'
  write ( *, '(a)' ) '  SINTMB does a backward transform.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of sequences is LOT =   ', lot
  write ( *, '(a,i8)' ) '  The length of each sequence is N = ', n
!
!  Set work arrays.
!
  lenr = n * lot
  lenwrk = lot * 2 * ( n + 2 )
  lensav = n / 2 + n + int ( log ( real ( n, kind = rk ) ) / log ( 2.0D+00 ) ) + 4

  write ( *, '(a,i8)' ) '  LENR =   ', lenr
  write ( *, '(a,i8)' ) '  LENSAV = ', lensav
  write ( *, '(a,i8)' ) '  LENWRK = ', lenwrk

  allocate ( r(1:n,1:lot) )
  allocate ( wsave(1:lensav) )
  allocate ( work(1:lenwrk) )

  call sintmi ( n, wsave, lensav, ier )
!
!  Set the data values.
!
  seed = 1973

  call random_number ( harvest = r(1:n,1:lot) )

  call r8mat_print_some ( n, lot, r, 1, 1, 5, 5, &
    '  Part of the original data:' )
!
!  Compute the FFT coefficients.
!
  jump = n
  inc = 1

  call sintmf ( lot, jump, n, inc, r, lenr, wsave, lensav, &
    work, lenwrk, ier )

  call r8mat_print_some ( n, lot, r, 1, 1, 5, 5, &
    '  Part of the FFT coefficients:' )
!
!  Compute inverse FFT of coefficients.  Should get back the
!  original data.
!
  call sintmb ( lot, jump, n, inc, r, lenr, wsave, lensav, &
    work, lenwrk, ier )

  call r8mat_print_some ( n, lot, r, 1, 1, 5, 5, &
    '  Part of the retrieved data:' )

  deallocate ( r )
  deallocate ( work )
  deallocate ( wsave )

  return
end
subroutine c8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! C8MAT_PRINT_SOME prints some of a C8MAT.
!
!  Discussion:
!
!    A C8MAT is a matrix of C8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns 
!    in the matrix.
!
!    Input, complex ( kind = ck ) A(M,N), the matrix.
!
!    Input, integer ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: incx = 4
  integer m
  integer n

  complex ( kind = ck ) a(m,n)
  character ( len = 20 ) ctemp(incx)
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
  character ( len = * )  title
  complex ( kind = ck ) zero

  zero = cmplx ( 0.0D+00, 0.0D+00, kind = rk )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)' 
    return
  end if
!
!  Print the columns of the matrix, in strips of INCX.
!
  do j2lo = jlo, min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i10,10x)' ) j
    end do

    write ( *, '(a,4a20)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi
!
!  Print out (up to) INCX entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == zero ) then
          ctemp(j2) = '       0.0          '
        else if ( imag ( a(i,j) ) == 0.0D+00 ) then
          write ( ctemp(j2), '(g10.3,10x)' ) real ( a(i,j), kind = rk )
        else
          write ( ctemp(j2), '(2g10.3)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,a1,4a20)' ) i, ':', ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine c8mat_uniform_01 ( m, n, seed, c )

!*****************************************************************************80
!
!! C8MAT_UNIFORM_01 returns a unit pseudorandom C8MAT.
!
!  Discussion:
!
!    A C8MAT is a matrix of C8's.
!
!    The angles should be uniformly distributed between 0 and 2 * PI,
!    the square roots of the radius uniformly distributed between 0 and 1.
!
!    This results in a uniform distribution of values in the unit circle.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns
!    in the matrix.
!
!    Input/output, integer SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, complex ( kind = ck ) C(M,N), the pseudorandom complex matrix.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  complex ( kind = ck ) c(m,n)
  integer i
  integer, parameter :: i4_huge = 2147483647
  integer j
  integer k
  real ( kind = rk ) r
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  integer seed
  real ( kind = rk ) theta

  do j = 1, n
    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge
      end if

      r = sqrt ( real ( seed, kind = rk ) * 4.656612875D-10 )

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge
      end if

      theta = 2.0D+00 * r8_pi * ( real ( seed, kind = rk ) * 4.656612875D-10 )

      c(i,j) = r * cmplx ( cos ( theta ), sin ( theta ), kind = rk )

    end do

  end do

  return
end
subroutine c8vec_print_part ( n, a, max_print, title )

!*****************************************************************************80
!
!! C8VEC_PRINT_PART prints "part" of a C8VEC.
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vector, is no more than MAX_PRINT, then
!    the entire vector is printed, one entry per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries of the vector.
!
!    Input, complex ( kind = ck ) A(N), the vector to be printed.
!
!    Input, integer MAX_PRINT, the maximum number of lines
!    to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  complex ( kind = ck ) a(n)
  integer i
  integer max_print
  character ( len = * )  title

  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  if ( n <= max_print ) then

    do i = 1, n
      write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)
    end do

  else if ( 3 <= max_print ) then

    do i = 1, max_print - 2
      write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)
    end do
    write ( *, '(a)' ) '  ........  ..............  ..............'
    i = n
    write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)

  else

    do i = 1, max_print - 1
      write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)
    end do
    i = max_print
    write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6,2x,a)' ) i, ':', a(i), &
      '...more entries...'

  end if

  return
end
subroutine c8vec_uniform_01 ( n, seed, c )

!*****************************************************************************80
!
!! C8VEC_UNIFORM_01 returns a unit pseudorandom C8VEC.
!
!  Discussion:
!
!    A C8VEC is a vector of C8's.
!
!    The angles should be uniformly distributed between 0 and 2 * PI,
!    the square roots of the radius uniformly distributed between 0 and 1.
!
!    This results in a uniform distribution of values in the unit circle.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of values to compute.
!
!    Input/output, integer SEED, the "seed" value, which should 
!    NOT be 0.  On output, SEED has been updated.
!
!    Output, complex ( kind = ck ) C(N), the pseudorandom complex vector.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  complex ( kind = ck ) c(n)
  integer i
  integer, parameter :: i4_huge = 2147483647
  integer k
  real    ( kind = rk ) r
  real    ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  integer seed
  real    ( kind = rk ) theta

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    r = sqrt ( real ( seed, kind = rk ) * 4.656612875D-10 )

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    theta = 2.0D+00 * r8_pi * ( real ( seed, kind = rk ) * 4.656612875D-10 )

    c(i) = r * cmplx ( cos ( theta ), sin ( theta ), kind = rk )

  end do

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
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
!  Parameters:
!
!    Input, integer M, the number of rows in A.
!
!    Input, integer N, the number of columns in A.
!
!    Input, real ( kind = rk ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  character ( len = * )  title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_PRINT_SOME prints some of an R8MAT.
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
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = rk ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ILO, JLO, the first row and column to print.
!
!    Input, integer IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
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
subroutine r8vec_print_part ( n, a, max_print, title )

!*****************************************************************************80
!
!! R8VEC_PRINT_PART prints "part" of an R8VEC.
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vector, is no more than MAX_PRINT, then
!    the entire vector is printed, one entry per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries of the vector.
!
!    Input, real ( kind = rk ) A(N), the vector to be printed.
!
!    Input, integer MAX_PRINT, the maximum number of lines
!    to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  integer i
  integer max_print
  character ( len = * ) title

  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  if ( n <= max_print ) then

    do i = 1, n
      write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)
    end do

  else if ( 3 <= max_print ) then

    do i = 1, max_print - 2
      write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)
    end do
    write ( *, '(a)' ) '  ........  ..............'
    i = n
    write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)

  else

    do i = 1, max_print - 1
      write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)
    end do
    i = max_print
    write ( *, '(2x,i8,a,1x,g14.6,2x,a)' ) i, ':', a(i), '...more entries...'

  end if

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
