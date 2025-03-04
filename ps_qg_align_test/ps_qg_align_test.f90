program main

!*****************************************************************************80
!
!! ps_qg_align_test() tests ps_qg_align().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ps_qg_align_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the PS_QG_ALIGN library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test07 ( )
  call test08 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PS_QG_ALIGN_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests PROFILE_SCORE_READ, PROFILE_SCORE_PRINT.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: acid_num = 23
  integer, parameter :: position_max = 200

  character acid_code(acid_num)
  character conserved(position_max)
  integer entropy(position_max)
  character ( len = 80 ) file_name
  integer gap_extend_percent(0:position_max)
  integer gap_open_percent(0:position_max)
  integer ios
  integer, parameter :: iunit = 1
  integer position_num
  integer score(position_max,acid_num)
  integer score2(acid_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  PROFILE_SCORE_READ reads profile scoring '
  write ( *, '(a)' ) '    information.'
  write ( *, '(a)' ) '  PROFILE_SCORE_PRINT prints the information.'

  file_name = 'profile.txt'
!
!  Read the data from the file.
!
  open ( unit = iunit, file = file_name, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST01 - Fatal error!'
    write ( *, '(a)' ) '  Could not open the data file.'
    return
  end if

  call profile_score_read ( acid_code, acid_num, conserved, entropy, &
    gap_open_percent, gap_extend_percent, iunit, position_max, position_num, &
    score, score2 )

  close ( unit = iunit )
!
!  Print the data.
!
  call profile_score_print ( acid_code, acid_num, conserved, entropy, &
    gap_open_percent, gap_extend_percent, position_max, position_num, &
    score, score2 )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests PS_QG_FOQ, PS_QG_FPQ, PS_QG_FSQ.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 4
  integer, parameter :: n = 8
  integer, parameter :: lds = m

  character b(n)
  real base
  real e(0:lds,0:n)
  real f(0:lds,0:n)
  real, parameter :: gap_extend = -0.5
  integer gap_extend_percent(0:m)
  real, parameter :: gap_open = -2.0
  integer gap_open_percent(0:m)
  real ge(0:m)
  real go(0:m)
  integer i
  integer i2
  integer j
  integer j2
  integer m1
  integer m2
  integer n1
  integer n2
  integer npath
  integer pathi(m+n+1)
  integer pathj(m+n+1)
  real, external :: ps_score_simple
  real s(0:lds,0:n)
  integer t(0:lds,0:n)

  b(1) = 'T'
  b(2) = 'G'
  b(3) = 'A'
  b(4) = 'C'
  b(5) = 'G'
  b(6) = 'C'
  b(7) = 'T'
  b(8) = 'C'

  gap_open_percent(0) =  95
  gap_open_percent(1) =  90
  gap_open_percent(2) =  75
  gap_open_percent(3) =  50
  gap_open_percent(4) =  25

  gap_extend_percent(0) = 100
  gap_extend_percent(1) =  80
  gap_extend_percent(2) =  70
  gap_extend_percent(3) =  40
  gap_extend_percent(4) =  85

  go(0:m) = gap_open * real ( gap_open_percent(0:m) ) / 100.0
  ge(0:m) = gap_extend * real ( gap_extend_percent(0:m) ) / 100.0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02:'
  write ( *, '(a)' ) '  Forward algorithms using quadratic space:'
  write ( *, '(a)' ) '  PS_QG_FSQ determines the alignment scores;'
  write ( *, '(a)' ) '  PS_QG_FOQ finds the optimal score;'
  write ( *, '(a)' ) '  PS_QG_FPQ determines the optimal alignment path.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' I  A  B'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i3,2x,a1)' ) i, b(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Matching Scores:'
  write ( *, '(a)' ) ' '
  do i = 0, m
    if ( i == 0 ) then
      write ( *, '(i3,13f5.1)' ) i, ( 0.0, j = 0, n )
    else
      write ( *, '(i3,13f5.1)' ) i, 0.0, &
        ( ps_score_simple ( i, b(j) ), j = 1, n )
    end if
  end do

  m1 = 0
  m2 = m
  n1 = 0
  n2 = n
  base = 0.0

  call ps_qg_fsq ( b, m, m1, m2, n, n1, n2, ps_score_simple, go, ge, &
    base, lds, s, e, f, t )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  PS_QG_FSQ:'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,13i5)' ) ( j, j = n1, n2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, s(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, e(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, f(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13i5)' ) i, t(i,n1:n2)
  end do

  call ps_qg_foq ( m, m1, m2, n, n1, n2, lds, s, i2, j2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PS_QG_FOQ reports the optimal score is ', s(i2,j2)
  write ( *, '(a,i6)' ) '  with I2 = ', i2
  write ( *, '(a,i6)' ) '  and J2 = ', j2

  call ps_qg_fpq ( i2, j2, m, m1, m2, n, n1, n2, lds, t, npath, pathi, pathj )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Matching path:'
  write ( *, '(a)' ) ' '
  do i = 1, npath
    write ( *, '(3i6)' ) i, pathi(i), pathj(i)
  end do

  call ps_qg_match_print ( b, m, n, npath, pathi, pathj, ps_score_simple, &
    go, ge )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests PS_QG_BOQ, PS_QG_BPQ, PS_QG_BSQ.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 4
  integer, parameter :: n = 8
  integer, parameter :: lds = m

  character b(n)
  real base
  real e(0:lds,0:n)
  real f(0:lds,0:n)
  real, parameter :: gap_extend = -0.5
  integer gap_extend_percent(0:m)
  real, parameter :: gap_open = -2.0
  integer gap_open_percent(0:m)
  real ge(0:m)
  real go(0:m)
  integer i
  integer i1
  integer j
  integer j1
  integer m1
  integer m2
  integer n1
  integer n2
  integer npath
  integer pathi(m+n+1)
  integer pathj(m+n+1)
  real, external :: ps_score_simple
  real s(0:lds,0:n)
  integer t(0:lds,0:n)

  b(1) = 'T'
  b(2) = 'G'
  b(3) = 'A'
  b(4) = 'C'
  b(5) = 'G'
  b(6) = 'C'
  b(7) = 'T'
  b(8) = 'C'

  gap_open_percent(0) =  95
  gap_open_percent(1) =  90
  gap_open_percent(2) =  75
  gap_open_percent(3) =  50
  gap_open_percent(4) =  25

  gap_extend_percent(0) = 100
  gap_extend_percent(1) =  80
  gap_extend_percent(2) =  70
  gap_extend_percent(3) =  40
  gap_extend_percent(4) =  85

  go(0:m) = gap_open * real ( gap_open_percent(0:m) ) / 100.0
  ge(0:m) = gap_extend * real ( gap_extend_percent(0:m) ) / 100.0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03:'
  write ( *, '(a)' ) '  Backward algorithms using quadratic space:'
  write ( *, '(a)' ) '  PS_QG_BSQ determines the alignment scores;'
  write ( *, '(a)' ) '  PS_QG_BOQ finds the optimal score;'
  write ( *, '(a)' ) '  PS_QG_BPQ determines the optimal alignment path.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Verify that these results match those for'
  write ( *, '(a)' ) '  PS_QG_FSQ/PS_QG_FOQ/PS_GG_FPQ on the same data.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' I  A  B'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i3,2x,a1)' ) i, b(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Matching Scores:'
  write ( *, '(a)' ) ' '
  do i = 0, m
    if ( i == 0 ) then
      write ( *, '(i3,13f5.1)' ) i, ( 0.0, j = 0, n )
    else
      write ( *, '(i3,13f5.1)' ) i, 0.0, &
        ( ps_score_simple ( i, b(j) ), j = 1, n )
    end if
  end do

  m1 = 0
  m2 = m
  n1 = 0
  n2 = n
  base = 0.0

  call ps_qg_bsq ( b, m, m1, m2, n, n1, n2, ps_score_simple, go, ge, &
    base, lds, s, e, f, t )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  PS_QG_BSQ:'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,13i5)' ) ( j, j = n1, n2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, s(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, e(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, f(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13i5)' ) i, t(i,n1:n2)
  end do

  call ps_qg_boq ( m, m1, m2, n, n1, n2, lds, s, i1, j1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PS_QG_BOQ reports the optimal score is ', s(i1,j1)
  write ( *, '(a,i6)' ) '  with I1 = ', i1
  write ( *, '(a,i6)' ) '  and J1 = ', j1

  call ps_qg_bpq ( i1, j1, m, m1, m2, n, n1, n2, lds, t, npath, pathi, pathj )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Matching path:'
  write ( *, '(a)' ) ' '
  do i = 1, npath
    write ( *, '(3i6)' ) i, pathi(i), pathj(i)
  end do

  call ps_qg_match_print ( b, m, n, npath, pathi, pathj, ps_score_simple, &
    go, ge )

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests PS_QG_RPL, PS_QG_FSL.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 4
  integer, parameter :: n = 8
  integer, parameter :: lds = m

  character b(n)
  real base
  real e(0:lds,0:n)
  real ev(0:n)
  real f(0:lds,0:n)
  real fv(0:n)
  real, parameter :: gap_extend = -0.5
  integer gap_extend_percent(0:m)
  real, parameter :: gap_open = -2.0
  integer gap_open_percent(0:m)
  real ge(0:m)
  real go(0:m)
  integer i
  integer i2
  integer j
  integer j2
  integer m1
  integer m2
  integer n1
  integer n2
  integer npath
  integer pathi(m+n+1)
  integer pathj(m+n+1)
  real, external :: ps_score_simple
  real s(0:lds,0:n)
  real s_max
  real sv(0:n)
  integer t(0:lds,0:n)
  integer tv(0:n)

  b(1) = 'T'
  b(2) = 'G'
  b(3) = 'A'
  b(4) = 'C'
  b(5) = 'G'
  b(6) = 'C'
  b(7) = 'T'
  b(8) = 'C'

  gap_open_percent(0) =  95
  gap_open_percent(1) =  90
  gap_open_percent(2) =  75
  gap_open_percent(3) =  50
  gap_open_percent(4) =  25

  gap_extend_percent(0) = 100
  gap_extend_percent(1) =  80
  gap_extend_percent(2) =  70
  gap_extend_percent(3) =  40
  gap_extend_percent(4) =  85

  go(0:m) = gap_open * real ( gap_open_percent(0:m) ) / 100.0
  ge(0:m) = gap_extend * real ( gap_extend_percent(0:m) ) / 100.0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04:'
  write ( *, '(a)' ) '  Linear space algorithms:'
  write ( *, '(a)' ) '  PS_QG_FSL determines forward alignment scores;'
  write ( *, '(a)' ) '  PS_QG_RPL determines an optimal alignment path.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Verify that these results match those for'
  write ( *, '(a)' ) '  PS_QG_FSQ/PS_GG_FPQ and PS_GG_BSQ/PS_GG_BPQ'
  write ( *, '(a)' ) '  on the same data.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' I  A  B'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i3,2x,a1)' ) i, b(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Matching Scores:'
  write ( *, '(a)' ) ' '
  do i = 0, m
    if ( i == 0 ) then
      write ( *, '(i3,13f5.1)' ) i, ( 0.0, j = 0, n )
    else
      write ( *, '(i3,13f5.1)' ) i, 0.0, &
        ( ps_score_simple ( i, b(j) ), j = 1, n )
    end if
  end do

  m1 = 0
  n1 = 0
  n2 = n
  base = 0.0

  do m2 = m1, m

    call ps_qg_fsl ( b, m, m1, m2, n, n1, n2, ps_score_simple, go, ge, &
      base, sv, ev, fv, tv, i2, j2, s_max )

    s(m2,0:n) = sv(0:n)
    e(m2,0:n) = ev(0:n)
    f(m2,0:n) = fv(0:n)
    t(m2,0:n) = tv(0:n)

  end do

  m2 = m

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  PS_QG_FSL:'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,13i5)' ) ( j, j = n1, n2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, s(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, e(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, f(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13i5)' ) i, t(i,n1:n2)
  end do

  call ps_qg_rpl ( b, m, m1, m2, n, n1, n2, ps_score_simple, go, ge, &
    npath, pathi, pathj )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Matching path:'
  write ( *, '(a)' ) ' '
  do i = 1, npath
    write ( *, '(3i6)' ) i, pathi(i), pathj(i)
  end do

  call ps_qg_match_print ( b, m, n, npath, pathi, pathj, ps_score_simple, &
    go, ge )

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests PS_QG_RPL, PS_QG_BSL.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 4
  integer, parameter :: n = 8
  integer, parameter :: lds = m

  character b(n)
  real base
  real e(0:lds,0:n)
  real ev(0:n)
  real f(0:lds,0:n)
  real fv(0:n)
  real, parameter :: gap_extend = -0.5
  integer gap_extend_percent(0:m)
  real, parameter :: gap_open = -2.0
  integer gap_open_percent(0:m)
  real ge(0:m)
  real go(0:m)
  integer i
  integer i1
  integer j
  integer j1
  integer m1
  integer m2
  integer n1
  integer n2
  integer npath
  integer pathi(m+n+1)
  integer pathj(m+n+1)
  real, external :: ps_score_simple
  real s(0:lds,0:n)
  real s_max
  real sv(0:n)
  integer t(0:lds,0:n)
  integer tv(0:n)

  b(1) = 'T'
  b(2) = 'G'
  b(3) = 'A'
  b(4) = 'C'
  b(5) = 'G'
  b(6) = 'C'
  b(7) = 'T'
  b(8) = 'C'

  gap_open_percent(0) =  95
  gap_open_percent(1) =  90
  gap_open_percent(2) =  75
  gap_open_percent(3) =  50
  gap_open_percent(4) =  25

  gap_extend_percent(0) = 100
  gap_extend_percent(1) =  80
  gap_extend_percent(2) =  70
  gap_extend_percent(3) =  40
  gap_extend_percent(4) =  85

  go(0:m) = gap_open * real ( gap_open_percent(0:m) ) / 100.0
  ge(0:m) = gap_extend * real ( gap_extend_percent(0:m) ) / 100.0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05:'
  write ( *, '(a)' ) '  Linear space algorithms:'
  write ( *, '(a)' ) '  PS_QG_BSL determines backward alignment scores;'
  write ( *, '(a)' ) '  PS_QG_RPL determines the optimal alignment path.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Verify that these results match those for'
  write ( *, '(a)' ) '  PS_QG_FSQ/PS_GG_FPQ, PS_GG_BSQ/PS_GG_BPQ'
  write ( *, '(a)' ) '  and PS_QG_FSL/PS_GG_RPL on the same data.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' I  A  B'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i3,2x,a1)' ) i, b(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Matching Scores:'
  write ( *, '(a)' ) ' '
  do i = 0, m
    if ( i == 0 ) then
      write ( *, '(i3,13f5.1)' ) i, ( 0.0, j = 0, n )
    else
      write ( *, '(i3,13f5.1)' ) i, 0.0, &
        ( ps_score_simple ( i, b(j) ), j = 1, n )
    end if
  end do

  m1 = 0
  m2 = m
  n1 = 0
  n2 = n
  base = 0.0

  do m1 = 0, m2

    call ps_qg_bsl ( b, m, m1, m2, n, n1, n2, ps_score_simple, go, ge, &
      base, sv, ev, fv, tv, i1, j1, s_max )

    s(m1,0:n) = sv(0:n)
    e(m1,0:n) = ev(0:n)
    f(m1,0:n) = fv(0:n)
    t(m1,0:n) = tv(0:n)

  end do

  m1 = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  PS_QG_BSL:'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,13i5)' ) ( j, j = n1, n2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, s(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, e(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, f(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13i5)' ) i, t(i,n1:n2)
  end do

  call ps_qg_rpl ( b, m, m1, m2, n, n1, n2, ps_score_simple, go, ge, &
    npath, pathi, pathj )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Matching path:'
  write ( *, '(a)' ) ' '
  do i = 1, npath
    write ( *, '(3i6)' ) i, pathi(i), pathj(i)
  end do

  call ps_qg_match_print ( b, m, n, npath, pathi, pathj, ps_score_simple, &
    go, ge )

  return
end
function ps_score_simple ( p1, c2 )

!*****************************************************************************80
!
!! PS_SCORE_SIMPLE computes a single entry profile/sequence matching score.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer P1, the position in the profile.
!
!    Input, character C2, the character in the sequence that is
!    to be matched with profile position P1.
!
!    Output, real PS_SCORE_SIMPLE, the score for matching the given profile
!    position with the sequence character.
!
  implicit none

  character c2
  integer p1
  real ps_score_simple
  real score

  score = 0.0E+00
!
!  Setting P1 = 0 to handle gaps?
!
       if ( p1 == 0 ) then

         if ( c2 == 'A' ) then
      score = 5.0
    else if ( c2 == 'C' ) then
      score = 2.0
    else if ( c2 == 'G' ) then
      score = 4.0
    else if ( c2 == 'T' ) then
      score = 1.0
    else
      score = 0.0
    end if

  else if ( p1 == 1 ) then

         if ( c2 == 'A' ) then
      score = 5.0
    else if ( c2 == 'C' ) then
      score = 2.0
    else if ( c2 == 'G' ) then
      score = 4.0
    else if ( c2 == 'T' ) then
      score = 1.0
    else
      score = 0.0
    end if

  else if ( p1 == 2 ) then

         if ( c2 == 'A' ) then
      score = 2.0
    else if ( c2 == 'C' ) then
      score = 4.0
    else if ( c2 == 'G' ) then
      score = 2.0
    else if ( c2 == 'T' ) then
      score = 1.0
    else
      score = 0.0
    end if

  else if ( p1 == 3 ) then

         if ( c2 == 'A' ) then
      score = 4.0
    else if ( c2 == 'C' ) then
      score = 2.0
    else if ( c2 == 'G' ) then
      score = 4.0
    else if ( c2 == 'T' ) then
      score = 1.0
    else
      score = 0.0
    end if

  else if ( p1 == 4 ) then

         if ( c2 == 'A' ) then
      score = 1.0
    else if ( c2 == 'C' ) then
      score = 1.0
    else if ( c2 == 'G' ) then
      score = 1.0
    else if ( c2 == 'T' ) then
      score = 3.0
    else
      score = 0.0
    end if

  else if ( p1 == 5 ) then

         if ( c2 == 'A' ) then
      score = 0.0
    else if ( c2 == 'C' ) then
      score = 0.0
    else if ( c2 == 'G' ) then
      score = 0.0
    else if ( c2 == 'T' ) then
      score = 0.0
    else
      score = 0.0
    end if

  end if

  ps_score_simple = score

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests PS_QG_FOQ, PS_QG_FPQ, PS_QG_FSQ.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n_max = 200
  integer, parameter :: m_max = 200
  integer, parameter :: lds = m_max

  character b(n_max)
  real base
  real, allocatable :: e(:,:)
  real, allocatable :: f(:,:)
  character ( len = 80 ) file_name
  real, parameter :: gap_extend = -0.5
  integer gap_extend_percent(0:m_max)
  real, parameter :: gap_open = -2.0
  integer gap_open_percent(0:m_max)
  real ge(0:m_max)
  real go(0:m_max)
  integer i
  integer i2
  integer ios
  integer :: iunit = 1
  integer j2
  integer m
  integer m1
  integer m2
  integer n
  integer n1
  integer n2
  integer npath
  integer pathi(m_max+n_max+1)
  integer pathj(m_max+n_max+1)
  real, external :: profile_score
  real, allocatable :: s(:,:)
  character ( len = 100 ) sequence
  character ( len = 20 ) sequence_name
  integer, allocatable :: t(:,:)
  integer test

  allocate ( e(0:lds,0:n_max) )
  allocate ( f(0:lds,0:n_max) )
  allocate ( s(0:lds,0:n_max) )
  allocate ( t(0:lds,0:n_max) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06:'
  write ( *, '(a)' ) '  Forward algorithms using quadratic space:'
  write ( *, '(a)' ) '  PS_QG_FSQ determines the alignment scores;'
  write ( *, '(a)' ) '  PS_QG_FOQ finds the optimal score;'
  write ( *, '(a)' ) '  PS_QG_FPQ determines the optimal alignment path.'

  file_name = 'profile.txt'

  open ( unit = iunit, file = file_name, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST06 - Fatal error!'
    write ( *, '(a)' ) '  Could not open the data file.'
    return
  end if

  call profile_score_read2 ( gap_open_percent, gap_extend_percent, iunit, &
    m_max, m )

  close ( unit = iunit )

  go(0:m) = gap_open * real ( gap_open_percent(0:m) ) / 100.0
  ge(0:m) = gap_extend * real ( gap_extend_percent(0:m) ) / 100.0

  do test = 1, 4
!
!  Set the sequence name and values.
!
    if ( test == 1 ) then

      sequence_name = 's00657 (chunk)'

      sequence = 'NGQSYRGTYSTTVTGRTCQAWSSMTPHSHS'

    else if ( test == 2 ) then

      sequence_name = 'c61545 (chunk)'

      sequence = 'DADKSPWCYTTDPRVRWEFCNLKK'

    else if ( test == 3 ) then

      sequence_name = 'plhu-k (chunk)'

      sequence = 'TPHRHQKTPENYPNAGLTMNYCRNPDADKG'

    else if ( test == 4 ) then

      sequence_name = 'Made-Up (chunk)'

      sequence = 'WSSMTPHRHQCATKTPENYPNAGLTM'

    end if
!
!  Convert the sequence string to a vector of characters.
!
    n = len_trim ( sequence )

    call s_to_chvec ( sequence, n, b )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Sequence name: ' // trim ( sequence_name )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Sequence:'
    write ( *, '(a)' ) ' '
    do i = 1, n
      write ( *, '(i5,5x,a1)' ) i, b(i)
    end do

    m1 = 0
    m2 = m
    n1 = 0
    n2 = n
    base = 0.0

    call ps_qg_fsq ( b, m, m1, m2, n, n1, n2, profile_score, go, ge, &
      base, lds, s, e, f, t )

    call ps_qg_foq ( m, m1, m2, n, n1, n2, lds, s, i2, j2 )

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  PS_QG_FOQ reports optimal score is ', s(i2,j2)
    write ( *, '(a,i6)' ) '  with I2 = ', i2
    write ( *, '(a,i6)' ) '  and J2 = ', j2

    call ps_qg_fpq ( i2, j2, m, m1, m2, n, n1, n2, lds, t, npath, pathi, pathj )

    call ps_qg_match_print ( b, m, n, npath, pathi, pathj, profile_score, &
      go, ge )

  end do

  deallocate ( e )
  deallocate ( f )
  deallocate ( s )
  deallocate ( t )

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests PS_QG_BOQ, PS_QG_BPQ, PS_QG_BSQ.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n_max = 200
  integer, parameter :: m_max = 200
  integer, parameter :: lds = m_max

  character b(n_max)
  real base
  real, allocatable :: e(:,:)
  real, allocatable :: f(:,:)
  character ( len = 80 ) file_name
  real, parameter :: gap_extend = -0.5
  integer gap_extend_percent(0:m_max)
  real, parameter :: gap_open = -2.0
  integer gap_open_percent(0:m_max)
  real ge(0:m_max)
  real go(0:m_max)
  integer i1
  integer ios
  integer :: iunit = 1
  integer j1
  integer m
  integer m1
  integer m2
  integer n
  integer n1
  integer n2
  integer npath
  integer pathi(m_max+n_max+1)
  integer pathj(m_max+n_max+1)
  real, external :: profile_score
  real, allocatable :: s(:,:)
  character ( len = 100 ) sequence
  character ( len = 20 ) sequence_name
  integer, allocatable :: t(:,:)
  integer test

  allocate ( e(0:lds,0:n_max) )
  allocate ( f(0:lds,0:n_max) )
  allocate ( s(0:lds,0:n_max) )
  allocate ( t(0:lds,0:n_max) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07:'
  write ( *, '(a)' ) '  Backward algorithms using quadratic space:'
  write ( *, '(a)' ) '  PS_QG_BSQ determines the alignment scores;'
  write ( *, '(a)' ) '  PS_QG_BOQ finds the optimal score;'
  write ( *, '(a)' ) '  PS_QG_BPQ determines the optimal alignment path.'

  file_name = 'profile.txt'

  open ( unit = iunit, file = file_name, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST07 - Fatal error!'
    write ( *, '(a)' ) '  Could not open the data file.'
    return
  end if

  call profile_score_read2 ( gap_open_percent, gap_extend_percent, iunit, &
    m_max, m )

  close ( unit = iunit )

  go(0:m) = gap_open * real ( gap_open_percent(0:m) ) / 100.0
  ge(0:m) = gap_extend * real ( gap_extend_percent(0:m) ) / 100.0

  do test = 1, 4
!
!  Set the sequence name and values.
!
    if ( test == 1 ) then

      sequence_name = 's00657 (chunk)'

      sequence = 'NGQSYRGTYSTTVTGRTCQAWSSMTPHSHS'

    else if ( test == 2 ) then

      sequence_name = 'c61545 (chunk)'

      sequence = 'DADKSPWCYTTDPRVRWEFCNLKK'

    else if ( test == 3 ) then

      sequence_name = 'plhu-k (chunk)'

      sequence = 'TPHRHQKTPENYPNAGLTMNYCRNPDADKG'

    else if ( test == 4 ) then

      sequence_name = 'Made-Up (chunk)'

      sequence = 'WSSMTPHRHQCATKTPENYPNAGLTM'

    end if
!
!  Convert the sequence string to a vector of characters.
!
    n = len_trim ( sequence )

    call s_to_chvec ( sequence, n, b )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Sequence name: ' // trim ( sequence_name )
    write ( *, '(a)' ) ' '

    m1 = 0
    m2 = m
    n1 = 0
    n2 = n
    base = 0.0

    call ps_qg_bsq ( b, m, m1, m2, n, n1, n2, profile_score, go, ge, &
      base, lds, s, e, f, t )

    call ps_qg_boq ( m, m1, m2, n, n1, n2, lds, s, i1, j1 )

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  PS_QG_BOQ reports optimal score is ', s(i1,j1)
    write ( *, '(a,i6)' ) '  with I1 = ', i1
    write ( *, '(a,i6)' ) '  and J1 = ', j1

    call ps_qg_bpq ( i1, j1, m, m1, m2, n, n1, n2, lds, t, npath, pathi, pathj )

    call ps_qg_match_print ( b, m, n, npath, pathi, pathj, profile_score, &
      go, ge )

  end do

  deallocate ( e )
  deallocate ( f )
  deallocate ( s )
  deallocate ( t )

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests PS_QG_FSL, PS_QG_RPL.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n_max = 200
  integer, parameter :: m_max = 200

  character b(n_max)
  character ( len = 80 ) file_name
  real, parameter :: gap_extend = -0.5
  integer gap_extend_percent(0:m_max)
  real, parameter :: gap_open = -2.0
  integer gap_open_percent(0:m_max)
  real ge(0:m_max)
  real go(0:m_max)
  integer ios
  integer :: iunit = 1
  integer m
  integer m1
  integer m2
  integer n
  integer n1
  integer n2
  integer npath
  integer pathi(m_max+n_max+1)
  integer pathj(m_max+n_max+1)
  real, external :: profile_score
  character ( len = 100 ) sequence
  character ( len = 20 ) sequence_name
  integer test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08:'
  write ( *, '(a)' ) '  Linear space algorithms:'
  write ( *, '(a)' ) '  PS_QG_FSL determines forward alignment scores;'
  write ( *, '(a)' ) '  PS_QG_RPL determines the optimal alignment path.'
  write ( *, '(a)' ) ' '

  file_name = 'profile.txt'

  open ( unit = iunit, file = file_name, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST08 - Fatal error!'
    write ( *, '(a)' ) '  Could not open the data file.'
    return
  end if

  call profile_score_read2 ( gap_open_percent, gap_extend_percent, iunit, &
    m_max, m )

  close ( unit = iunit )

  go(0:m) = gap_open * real ( gap_open_percent(0:m) ) / 100.0
  ge(0:m) = gap_extend * real ( gap_extend_percent(0:m) ) / 100.0

  do test = 1, 4
!
!  Set the sequence name and values.
!
    if ( test == 1 ) then

      sequence_name = 's00657 (chunk)'

      sequence = 'NGQSYRGTYSTTVTGRTCQAWSSMTPHSHS'

    else if ( test == 2 ) then

      sequence_name = 'c61545 (chunk)'

      sequence = 'DADKSPWCYTTDPRVRWEFCNLKK'

    else if ( test == 3 ) then

      sequence_name = 'plhu-k (chunk)'

      sequence = 'TPHRHQKTPENYPNAGLTMNYCRNPDADKG'

    else if ( test == 4 ) then

      sequence_name = 'Made-Up (chunk)'

      sequence = 'WSSMTPHRHQCATKTPENYPNAGLTM'

    end if
!
!  Convert the sequence string to a vector of characters.
!
    n = len_trim ( sequence )

    call s_to_chvec ( sequence, n, b )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Sequence name: ' // trim ( sequence_name )
    write ( *, '(a)' ) ' '

    m1 = 0
    m2 = m
    n1 = 0
    n2 = n

    call ps_qg_rpl ( b, m, m1, m2, n, n1, n2, profile_score, go, ge, &
      npath, pathi, pathj )

    call ps_qg_match_print ( b, m, n, npath, pathi, pathj, profile_score, &
      go, ge )

  end do

  return
end
function profile_score ( p1, c2 )

!*****************************************************************************80
!
!! PROFILE_SCORE computes a single entry profile/sequence matching score.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer P1, the position in the profile.
!
!    Input, character C2, the character in the sequence that is
!    to be matched with profile position P1.
!
!    Output, real PROFILE_SCORE, the score for matching the given profile
!    position with the sequence character.
!
  implicit none

  integer, parameter :: acid_num = 23
  integer, parameter :: m_max = 200

  integer a_to_i4
  character, save, dimension ( acid_num ) :: acid_code
  integer, save, dimension ( 26 ) :: acid_index
  character c2
  character, save, dimension ( m_max ) :: conserved
  logical, parameter :: debug = .false.
  integer, save, dimension ( m_max ) :: entropy
  character ( len = 80 ) file_name
  integer, save, dimension ( m_max ) :: gap_extend_percent
  integer, save, dimension ( m_max ) :: gap_open_percent
  integer i
  character i4_to_a
  integer ic2
  integer index_c2
  logical, save :: initialized = .false.
  integer ios
  integer, parameter :: iunit = 1
  integer, save :: m
  integer p1
  real profile_score
  integer, save, dimension ( m_max, acid_num ) :: score1
  integer, save, dimension ( acid_num ) :: score2

  if ( .not. initialized ) then

    file_name = 'profile.txt'
!
!  Read the data from the file.
!
    open ( unit = iunit, file = file_name, status = 'old', iostat = ios )

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PROFILE_SCORE - Fatal error!'
      write ( *, '(a)' ) '  Could not open the data file.'
      return
    end if

    call profile_score_read ( acid_code, acid_num, conserved, entropy, &
      gap_open_percent, gap_extend_percent, iunit, m_max, m, score1, score2 )

    close ( unit = iunit )

    call a_index ( acid_num, acid_code, acid_index )

    if ( debug ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    I    Acid_code'
      write ( *, '(a)' ) ' '
      do i = 1, acid_num
        write ( *, '(i5,5x,a1)' ) i, acid_code(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    I    Acid_index'
      write ( *, '(a)' ) ' '
      do i = 1, 26
        write ( *, '(i5,5x,a1,4x,i5)' ) i, i4_to_a(i), acid_index(i)
      end do

    end if

    initialized = .true.

  end if

  ic2 = a_to_i4 ( c2 )

  index_c2 = acid_index(ic2)

  if ( index_c2 <= 0 ) then
    profile_score = 0.0
  else if ( p1 <= 0 ) then
    profile_score = 0.0
  else
    profile_score = score1(p1,index_c2)
  end if

  return
end
