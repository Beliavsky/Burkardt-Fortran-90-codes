subroutine a_index ( acid_num, acid_code, acid_index )

!*****************************************************************************80
!
!! A_INDEX sets up a reverse index for the amino acid codes.
!
!  Example:
!
!    Input:
!
!      ACID_CODE =
!        'A', 'R', 'N', 'B', 'D', 'C', 'Q', 'Z', 'E', 'G',
!        'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W',
!        'Y', 'V', 'X'
!
!    Output:
!
!      ACID_INDEX =
!        1,   4,   6,   5,   9,  16,  10,  11,  12,   0,
!       14,  13,  15,   3,   0,  17,   7,   2,  18,  19,
!        0,  22,  20,  23,  21,   8.
!
!  Discussion:
!
!    ACID_CODE allows us to discover the index of item 'R' only
!    by searching through all the entries of ACID_CODE until we
!    encounter an 'R' or reach the end of the array.
!
!    ACID_INDEX allows us to discover the index of item 'R' by
!    converting 'R' to its numeric index of 17 and evaluating
!    ACID_INDEX(17), which is 2.
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
!    Input, integer ACID_NUM, the number of entries in ACID_CODE.
!
!    Input, character ACID_CODE(ACID_NUM), a list of alphabetic
!    characters.  Normally, this list is upper case, and contains a
!    subset of the alphabetic characters 'A' through 'Z'.
!
!    Output, integer ACID_INDEX(26), indicates, for each alphabetic
!    character, the (last) index of ACID_CODE containing that character,
!    or 0 if the character does not occur in ACID_CODE.
!
  implicit none

  integer acid_num

  character a
  integer a_to_i4
  character acid_code(acid_num)
  integer acid_index(26)
  integer i
  integer j

  acid_index(1:26) = 0

  do i = 1, acid_num

    a = acid_code(i)
    j = a_to_i4 ( a )

    if ( j < 1 .or. 26 < j ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'A_INDEX - Fatal error!'
      write ( *, '(a,i6)' ) '  Out-of-bounds acid index J = ', j
      write ( *, '(a)' ) '  Originating acid code character is ' // a
      write ( *, '(a,i6)' ) '  ASCII code = ', ichar ( a )
      stop
    end if

    acid_index(j) = i

  end do

  return
end
function a_to_i4 ( a )

!*****************************************************************************80
!
!! A_TO_I4 returns the index of an alphabetic character.
!
!  Example:
!
!    A  A_TO_I4
!
!    'a'   1
!    'b'   2
!    ...
!    'z'  26
!    'A'   1
!    'B'   2
!    ...
!    'Z'  26
!    '$'  -1
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character A, a character.
!
!    Output, integer A_TO_I4, is the alphabetic index of the 
!    character, between 1 and 26 if the character is alphabetic, and 
!    -1 otherwise.
!
  implicit none

  character a
  integer a_to_i4

  if ( lle ( 'A', a ) .and. lle ( a, 'Z' ) ) then
    a_to_i4 = 1 + ichar ( a ) - ichar ( 'A' )
  else if ( lle ( 'a', a ) .and. lle ( a, 'z' ) ) then
    a_to_i4 = 1 + ichar ( a ) - ichar ( 'a' )
  else
    a_to_i4 = -1
  end if

  return
end
subroutine ch_cap ( c )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, the character to capitalize.
!
  implicit none

  character c
  integer itemp

  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end
subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! I4_SWAP switches two I4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer I, J.  On output, the values of I and
!    J have been interchanged.
!
  implicit none

  integer i
  integer j
  integer k

  k = i
  i = j
  j = k

  return
end
function i4_to_a ( i )

!*****************************************************************************80
!
!! I4_TO_A returns the I-th alphabetic character.
!
!  Example:
!
!    I  I4_TO_A
!
!   -8  ' '
!    0  ' '
!    1  'A'
!    2  'B'
!   ..
!   26  'Z'
!   27  'a'
!   52  'z'
!   53  ' '
!   99  ' '
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, the index of the letter to be returned.
!    0 is a space;
!    1 through 26 requests 'A' through 'Z', (ASCII 65:90);
!    27 through 52 requests 'a' through 'z', (ASCII 97:122);
!
!    Output, character I4_TO_A, the requested alphabetic letter.
!
  implicit none

  integer, parameter :: cap_shift = 64
  integer i
  character i4_to_a
  integer, parameter :: low_shift = 96

  if ( i <= 0 ) then
    i4_to_a = ' '
  else if ( 1 <= i .and. i <= 26 ) then
    i4_to_a = char ( cap_shift + i )
  else if ( 27 <= i .and. i <= 52 ) then
    i4_to_a = char ( low_shift + i - 26 )
  else if ( i >= 53 ) then
    i4_to_a = ' '
  end if

  return
end
subroutine i4vec_reverse ( n, a )

!*****************************************************************************80
!
!! I4VEC_REVERSE reverses the elements of an integer vector.
!
!  Example:
!
!    Input:
!
!      N = 5,
!      A = ( 11, 12, 13, 14, 15 ).
!
!    Output:
!
!      A = ( 15, 14, 13, 12, 11 ).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input/output, integer A(N), the array to be reversed.
!
  implicit none

  integer n

  integer a(n)
  integer i

  do i = 1, n/2
    call i4_swap ( a(i), a(n+1-i) )
  end do

  return
end
subroutine i4vec2_compare ( n, a1, a2, i, j, isgn )

!*****************************************************************************80
!
!! I4VEC2_COMPARE compares pairs of integers stored in two vectors.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data items.
!
!    Input, integer A1(N), A2(N), contain the two components 
!    of each item.
!
!    Input, integer I, J, the items to be compared.
!
!    Output, integer ISGN, the results of the comparison:
!    -1, item I < item J,
!     0, item I = item J,
!    +1, item I > item J.
!
  implicit none

  integer n

  integer a1(n)
  integer a2(n)
  integer i
  integer isgn
  integer j

  isgn = 0

       if ( a1(i) < a1(j) ) then

    isgn = -1

  else if ( a1(i) == a1(j) ) then

         if ( a2(i) < a2(j) ) then
      isgn = -1
    else if ( a2(i) < a2(j) ) then
      isgn = 0
    else if ( a2(i) > a2(j) ) then
      isgn = +1
    end if

  else if ( a1(i) > a1(j) ) then

    isgn = +1

  end if

  return
end
subroutine i4vec2_sort_a ( n, a1, a2 )

!*****************************************************************************80
!
!! I4VEC2_SORT_A ascending sorts a vector of pairs of integers.
!
!  Discussion:
!
!    Each item to be sorted is a pair of integers (I,J), with the I
!    and J values stored in separate vectors A1 and A2.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of items of data.
!
!    Input/output, integer A1(N), A2(N), the data to be sorted.
!
  implicit none

  integer n

  integer a1(n)
  integer a2(n)
  integer i
  integer indx
  integer isgn
  integer j
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( indx > 0 ) then

      call i4_swap ( a1(i), a1(j) )
      call i4_swap ( a2(i), a2(j) )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4vec2_compare ( n, a1, a2, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine profile_score_print ( acid_code, acid_num, conserved, entropy, &
  gap_open_percent, gap_extend_percent, position_max, position_num, score, &
  score2 )

!*****************************************************************************80
!
!! PROFILE_SCORE_PRINT prints profile scoring data.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ACID_CODE(ACID_NUM), the nucleic acid codes.
!
!    Input, integer ACID_NUM, the number of nucleic acids.
!
!    Input, character CONSERVED(POSITION_MAX), contains an indicator
!    of the conserved nucleic acid at each position.
!
!    Input, integer ENTROPY(POSITION_MAX), contains an entropy 
!    measure associated with each position.
!
!    Input, integer GAP_OPEN_PERCENT(0:POSITION_MAX), a percentage
!    multiplier for the gap open penalty for each position.
!
!    Input, integer GAP_EXTEND_PERCENT(0:POSITION_MAX), a 
!    percentage multiplier for the gap extend penalty for each position.
!
!    Input, integer POSITION_MAX, the maximum number of profile 
!    positions for which the calling program has set aside storage.
!
!    Input, integer POSITION_NUM, the number of positions in 
!    the profile.
!
!    Input, integer SCORE(POSITION_MAX,ACID_NUM), the profile 
!    scoring matrix.
!
!    Input, integer SCORE2(ACID_NUM), more scoring data.
!
  implicit none

  integer acid_num
  integer position_max

  character acid_code(acid_num)
  character conserved(position_max)
  integer entropy(position_max)
  integer gap_extend_percent(0:position_max)
  integer gap_open_percent(0:position_max)
  integer i
  integer j
  integer position_num
  integer score(position_max,acid_num)
  integer score2(acid_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Profile Score Print'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of positions in profile = ', position_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Pos Cons  Entropy       Gap_Open    Gap_Extend'
  write ( *, '(a)' ) '                          Percent     Percent'
  write ( *, '(a)' ) ' '

  i = 0
  write ( *, '(i5,4x,1x,10x,2i10)' ) i, &
    gap_open_percent(i), gap_extend_percent(i)
  do i = 1, position_num
    write ( *, '(i5,4x,a1,3i10)' ) i, conserved(i), entropy(i), &
      gap_open_percent(i), gap_extend_percent(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Profile Position Scoring Matrix:'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,a3,2x,23(4x,a1))' ) 'Pos', acid_code(1:acid_num)
  write ( *, '(a)' ) ' '

  do i = 1, position_num
    write ( *, '(i6,2x,23i5)' ) i, score(i,1:acid_num)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(5x,a1,2x,23i5)' ) '*', score2(1:acid_num)

  return
end
subroutine profile_score_read ( acid_code, acid_num, conserved, entropy, &
  gap_open_percent, gap_extend_percent, iunit, position_max, position_num, &
  score, score2 )

!*****************************************************************************80
!
!! PROFILE_SCORE_READ reads profile scoring data from a file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ACID_CODE(ACID_NUM), the nucleic acid codes.
!
!    Input, integer ACID_NUM, the number of nucleic acids.
!
!    Output, character CONSERVED(POSITION_MAX), contains an indicator
!    of the conserved nucleic acid at each position.
!
!    Output, integer ENTROPY(POSITION_MAX), contains an entropy 
!    measure associated with each position.
!
!    Output, integer GAP_OPEN_PERCENT(0:POSITION_MAX), a percentage
!    multiplier for the gap open penalty for each position.
!
!    Output, integer GAP_EXTEND_PERCENT(0:POSITION_MAX), a 
!    percentage multiplier for the gap extend penalty for each position.
!
!    Input, integer IUNIT, the unit number associated with 
!    the file.
!
!    Input, integer POSITION_MAX, the maximum number of profile 
!    positions for which the calling program has set aside storage.
!
!    Output, integer POSITION_NUM, the number of positions in 
!    the profile.
!
!    Output, integer SCORE(POSITION_MAX,ACID_NUM), the profile 
!    scoring matrix.
!
!    Output, integer SCORE2(ACID_NUM), more scoring data.
!
  implicit none

  integer acid_num
  integer position_max

  character acid_code(acid_num)
  integer acid_i
  character conserved(position_max)
  logical, parameter :: debug = .false.
  logical done
  integer entropy(position_max)
  logical found_start
  integer gap_extend_percent(0:position_max)
  integer gap_open_percent(0:position_max)
  integer ierror
  integer ios
  integer iunit
  integer ival
  integer last
  character ( len = 255 ) line
  integer nrec
  integer position_num
  logical s_eqi
  integer score(position_max,acid_num)
  integer score2(acid_num)
  character ( len = 255 ) word
!
!  The file format does not include gap data for entry 0.
!  We want some, so make it up and stick it in early.
!
  gap_open_percent(0) = 100
  gap_extend_percent(0) = 100

  found_start = .false.
  position_num = 0
  nrec = 0

  do

    done = .true.
    read ( iunit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    nrec = nrec + 1

    call word_next_read ( line, word, done )
!
!  If the first word is CONS then we found the starting point of the data.
!
    if ( s_eqi ( word, 'CONS' ) ) then

      if ( found_start ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PROFILE_SCORE_READ - Fatal error!'
        stop
      else
        found_start = .true.
      end if

      do acid_i = 1, acid_num
        call word_next_read ( line, word, done )
        acid_code(acid_i) = word(1:1)
      end do
!
!  If already found start, then...
!
    else if ( found_start ) then
!
!    ...if first character is exclamation, read entropy
!
      if ( word == '!' ) then

        call word_last_read ( line, word )
        call s_to_i4 ( word, ival, ierror, last )
        entropy(position_num) = ival
!
!    ...Else read last line
!
      else if ( word == '*' ) then

        do acid_i = 1, acid_num

          call word_next_read ( line, word, done )
          call s_to_i4 ( word, ival, ierror, last )
          score2(acid_i) = ival

        end do

        exit
!
!    ...Else read scores
!
      else

        position_num = position_num + 1

        if ( position_num > position_max ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'PROFILE_SCORE_READ - Fatal error!'
          write ( *, '(a)' ) '  Not enough space has been provided to'
          write ( *, '(a)' ) '  store data from the file.'
          write ( *, '(a,i6)' ) '  The provided space is POSITION_MAX = ', &
            position_max
          stop
        end if

        conserved(position_num) = word

        do acid_i = 1, acid_num

          call word_next_read ( line, word, done )
          call s_to_i4 ( word, ival, ierror, last )
          score(position_num,acid_i) = ival

        end do

        call word_next_read ( line, word, done )
        call s_to_i4 ( word, ival, ierror, last )
        gap_open_percent(position_num) = ival

        call word_next_read ( line, word, done )
        call s_to_i4 ( word, ival, ierror, last )
        gap_extend_percent(position_num) = ival

      end if
!
!  Else
!
    else

    end if

  end do

  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PROFILE_SCORE_READ - Note:'
    write ( *, '(a,i6)' ) '  Number of records read was ', nrec
  end if

  return
end
subroutine profile_score_read2 ( gap_open_percent, gap_extend_percent, iunit, &
  position_max, position_num )

!*****************************************************************************80
!
!! PROFILE_SCORE_READ2 returns a small amount of information from a profile.
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
!    Output, integer GAP_OPEN_PERCENT(0:POSITION_MAX), a percentage
!    multiplier for the gap open penalty for each position.
!
!    Output, integer GAP_EXTEND_PERCENT(0:POSITION_MAX), a
!    percentage multiplier for the gap extend penalty for each position.
!
!    Input, integer IUNIT, the unit number associated with 
!    the file.
!
!    Input, integer POSITION_MAX, the maximum number of 
!    profile positions for which the calling program has set aside storage.
!
!    Output, integer POSITION_NUM, the number of positions in 
!    the profile.
!
  implicit none

  integer, parameter :: acid_num = 23
  integer position_max

  character acid_code(acid_num)
  character conserved(position_max)
  integer entropy(position_max)
  integer gap_extend_percent(0:position_max)
  integer gap_open_percent(0:position_max)
  integer iunit
  integer position_num
  integer score(position_max,acid_num)
  integer score2(acid_num)

  call profile_score_read ( acid_code, acid_num, conserved, entropy, &
    gap_open_percent, gap_extend_percent, iunit, position_max, position_num, &
    score, score2 )

  return
end
subroutine ps_gg_bsl ( b, m, m1, m2, n, n1, n2, ps_score, go, ge, &
  base, s, e, f, t )

!*****************************************************************************80
!
!! PS_GG_BSL determines a global gap backward alignment score in linear space.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Kun-Mao Chao, Ross Hardison, Webb Miller,
!    Recent Developments in Linear-Space Alignment Methods: A Survey,
!    Journal of Computational Biology,
!    Volume 1, Number 4, 1994, pages 271-291.
!
!    Michael Waterman,
!    Introduction to Computational Biology,
!    Chapman and Hall, 1995.
!
!  Parameters:
!
!    Input, character B(N), the sequence to be aligned.
!
!    Input, integer M, the number of entries in the profile.
!
!    Input, integer M1, M2, the first and last rows of the
!    table to compute.
!
!    Input, integer N, the number of entries in the sequence.
!
!    Input, integer N1, N2, the first and last columns of
!    the table to compute.
!
!    Input, external PS_SCORE, the name of a function of the form
!      function ps_score ( p1, c2 )
!    which returns a real value PS_SCORE for the matching of the position
!    P1 from the profile to character C2 from the sequence B.
!
!    Input, real GO(0:M), the profile position gap open penalty.
!
!    Input, real GE(0:M), the profile position gap extend penalty.
!
!    Input, real BASE, an initial quantity to be added to 
!    certain penalties.
!
!    Output, real S(0:N), the backward optimal score vector.
!    The maximum possible alignment score is in S(N1).
!
!    Output, real E(0:N), the backward final insertion 
!    score vector.
!
!    Output, real F(0:N), the backward final deletion score vector.
!
!    Output, integer T(0:N), the backward pointer vector.
!
  implicit none

  integer, parameter :: DUNNO = 0
  integer, parameter :: MATCH = 1
  integer, parameter :: INSERT = 2
  integer, parameter :: DELETE = 4
  integer, parameter :: TERMINATE = 8

  integer m
  integer n

  character b(n)
  real base
  real diag_new
  real diag_old
  real e(0:n)
  real f(0:n)
  real ge(0:m)
  real go(0:m)
  integer i
  integer j
  integer m1
  integer m2
  integer n1
  integer n2
  real, external :: ps_score
  real s(0:n)
  integer t(0:n)
!
!  The last row, I = M2.
!
  e(n2) = 0.0
  f(n2) = 0.0
  s(n2) = 0.0
  t(n2) = DUNNO

  if ( n2-1 >= n1 ) then
    e(n2-1) = e(n2) + go(m2) + ge(m2)
    f(n2-1) = f(n2) + 2.0 * go(m2) + ge(m2)
    s(n2-1) = s(n2) + go(m2) + ge(m2)
    t(n2-1) = INSERT
  end if

  do j = n2-2, n1, -1
    e(j) = e(j+1) + ge(m2)
    f(j) = f(j+1) + ge(m2)
    s(j) = s(j+1) + ge(m2)
    t(j) = INSERT
  end do
!
!  Upper rectangle.
!
  do i = m2-1, m1, -1

    diag_old = s(n2)

    if ( i == m2-1 ) then
      e(n2) = base + ge(i) + 2.0 * go(i) -       go(i+1)
      f(n2) = base + ge(i) + 2.0 * go(i) - 2.0 * go(i+1)
      s(n2) = base + ge(i)
      t(n2) = DELETE
    else
      e(n2) = e(n2) + ge(i) + 2.0 * go(i) - 2.0 * go(i+1)
      f(n2) = f(n2) + ge(i) + 2.0 * go(i) - 2.0 * go(i+1)
      s(n2) = s(n2) + ge(i)
      t(n2) = DELETE
    end if

    do j = n2-1, n1, -1
!
!  Insertion.
!
      e(j) = e(j+1) + ge(i)

      if ( s(j+1) + go(i) + ge(i) > e(j) ) then
        e(j) = s(j+1) + go(i) + ge(i)
      end if
!
!  Deletion.
!
      f(j) = f(j) + ge(i) + go(i) - go(i+1)

      if ( s(j) + go(i) + ge(i) > f(j) ) then
        f(j) = s(j) + go(i) + ge(i)
      end if
!
!  Best.
!
      diag_new = s(j)

      s(j) = diag_old + ps_score ( i+1, b(j+1) )
      t(j) = MATCH

      if ( s(j) == e(j) ) then
        t(j) = t(j) + INSERT
      else if ( s(j) < e(j) ) then
        s(j) = e(j)
        t(j) = INSERT
      end if

      if ( s(j) == f(j) ) then
        t(j) = t(j) + DELETE
      else if ( s(j) < f(j) ) then
        s(j) = f(j)
        t(j) = DELETE
      end if

      diag_old = diag_new

    end do
  end do

  return
end
subroutine ps_gg_fsl ( b, m, m1, m2, n, n1, n2, ps_score, go, ge, &
  base, s, e, f, t )

!*****************************************************************************80
!
!! PS_GG_FSL determines a global gap forward alignment score in linear space.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Kun-Mao Chao, Ross Hardison, Webb Miller,
!    Recent Developments in Linear-Space Alignment Methods: A Survey,
!    Journal of Computational Biology,
!    Volume 1, Number 4, 1994, pages 271-291.
!
!    Michael Waterman,
!    Introduction to Computational Biology,
!    Chapman and Hall, 1995.
!
!  Parameters:
!
!    Input, character B(N), the sequence to be aligned.
!
!    Input, integer M, the number of entries in the profile.
!
!    Input, integer M1, M2, the first and last rows of the
!    table to compute.
!
!    Input, integer N, the number of entries in the sequence.
!
!    Input, integer N1, N2, the first and last columns of the
!    table to compute.
!
!    Input, external PS_SCORE, the name of a function of the form
!      function ps_score ( p1, c2 )
!    which returns a real value PS_SCORE for the matching of the position
!    P1 from the profile to character C2 from the sequence B.
!
!    Input, real GO(0:M), the profile position gap open penalty.
!
!    Input, real GE(0:M), the profile position gap extend penalty.
!
!    Input, real BASE, an initial quantity to be added to 
!    certain penalties.
!
!    Output, real S(0:N), the forward optimal score vector.
!    The maximum possible alignment score is in S(N2).
!
!    Output, real E(0:N), the forward final insertion score vector.
!
!    Output, real F(0:N), the forward final deletion score vector.
!
!    Output, integer T(0:N), the forward pointer vector.
!
  implicit none

  integer, parameter :: DUNNO = 0
  integer, parameter :: MATCH = 1
  integer, parameter :: INSERT = 2
  integer, parameter :: DELETE = 4
  integer, parameter :: TERMINATE = 8

  integer m
  integer n

  character b(n)
  real base
  real diag_new
  real diag_old
  real e(0:n)
  real f(0:n)
  real ge(0:m)
  real go(0:m)
  real horz
  integer i
  integer j
  integer m1
  integer m2
  integer n1
  integer n2
  real, external :: ps_score
  real s(0:n)
  integer t(0:n)
!
!  The first row, I = M1.
!
  e(n1) = 0.0
  f(n1) = 0.0
  s(n1) = 0.0
  t(n1) = DUNNO

  if ( n1+1 <= n2 ) then
    e(n1+1) = e(n1) + go(m1) + ge(m1)
    f(n1+1) = f(n1) + 2.0 * go(m1) + ge(m1)
    s(n1+1) = s(n1) + go(m1) + ge(m1)
    t(n1+1) = INSERT
  end if

  do j = n1+2, n2
    e(j) = e(j-1) + ge(m1)
    f(j) = f(j-1) + ge(m1)
    s(j) = s(j-1) + ge(m1)
    t(j) = INSERT
  end do
!
!  Subsequent rows.
!
  do i = m1+1, m2

    diag_old = s(n1)

    if ( i == m1+1 ) then
      e(n1) = base + ge(i-1) + go(i-1)
      f(n1) = base + ge(i-1)
      s(n1) = base + ge(i-1)
      t(n1) = DELETE
    else
      e(n1) = e(n1) + ge(i-1)
      f(n1) = f(n1) + ge(i-1)
      s(n1) = s(n1) + ge(i-1)
      t(n1) = DELETE
    end if

    do j = n1+1, n2
!
!  Insertion
!
      e(j) = e(j-1) + ge(i)

      if ( s(j-1) + go(i) + ge(i) > e(j) ) then
        e(j) = s(j-1) + go(i) + ge(i)
      end if
!
!  Deletion.
!
      f(j) = f(j) + ge(i-1)

      if ( s(j) + go(i-1) + ge(i-1) > f(j) ) then
        f(j) = s(j) + go(i-1) + ge(i-1)
      end if
!
!  Best.
!
      diag_new = s(j)

      s(j) = diag_old + ps_score ( i, b(j) )
      t(j) = MATCH

      if ( s(j) == e(j) ) then
        t(j) = t(j) + INSERT
      else if ( s(j) < e(j) ) then
        s(j) = e(j)
        t(j) = INSERT
      end if

      if ( s(j) == f(j) ) then
        t(j) = t(j) + DELETE
      else if ( s(j) < f(j) ) then
        s(j) = f(j)
        t(j) = DELETE
      end if

      diag_old = diag_new

    end do
  end do

  return
end
subroutine ps_qg_boq ( m, m1, m2, n, n1, n2, lds, s, i1, j1 )

!*****************************************************************************80
!
!! PS_QG_BOQ: backward endpoint of a quasiglobal optimal local alignment.
!
!  Discussion:
!
!    The algorithm examines the first column and row of the quadratic score
!    table, and returns I1, J1, the location of the maximum value.
!
!    Actually, for generality, the algorithm examines row M1 and column N1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of entries in sequence A.
!
!    Input, integer M1, M2, the minimum and maximum rows.
!
!    Input, integer N, the number of entries in sequence B.
!
!    Input, integer N1, N2, the minimum and maximum columns.
!
!    Input, integer LDS, the declared upper first dimension of S.
!    LDS must be at least M.
!
!    Output, real S(0:LDS,0:N), the backward optimal score table.
!    Entries in the M1:M2 by N1:N2 block have been computed.
!
!    Output, integer I1, J1, the beginning point of an optimal
!    quasiglobal alignment.
!    I1 is equal to M1, or J1 is equal to N1.
!
  implicit none

  integer lds
  integer n

  integer i
  integer i1
  integer j
  integer j1
  integer m
  integer m1
  integer m2
  integer n1
  integer n2
  integer s(0:lds,0:n)

  i1 = m2
  j1 = n1

  do i = m2-1, m1, -1
    if ( s(i,n1) > s(i1,j1) ) then
      i1 = i
      j1 = n1
    end if
  end do

  do j = n1+1, n2
    if ( s(m1,j) > s(i1,j1) ) then
      i1 = m1
      j1 = j
    end if
  end do

  return
end
subroutine ps_qg_bpq ( i1, j1, m, m1, m2, n, n1, n2, lds, t, npath, pathi, &
  pathj )

!*****************************************************************************80
!
!! PS_QG_BPQ: a quasiglobal gap backward alignment path in quadratic space.
!
!  Discussion:
!
!    The routine must be given a starting point.
!
!    An effort has been made to handle the ambiguous case, where
!    more than one optimal path enters a cell.  In such a case,
!    the code tries to take a delete out if there was a delete in,
!    or an insert out if there was an insert in, since the optimal
!    score calculation includes a penalty for gap opening.
!
!    The routine is called "quadratic" because it uses an M by N array
!    to do the alignment.
!
!    The score table must have been computed before this routine is called.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Kun-Mao Chao, Ross Hardison, Webb Miller,
!    Recent Developments in Linear-Space Alignment Methods: A Survey,
!    Journal of Computational Biology,
!    Volume 1, Number 4, 1994, pages 271-291.
!
!  Parameters:
!
!    Input, integer I1, J1, the starting point of the path.
!    M1 <= I1 <= M2, N1 <= J1 <= N2.
!
!    Input, integer M, the number of entries in sequence A.
!
!    Input, integer M1, M2, the minimum and maximum rows
!    of the computed score table.
!
!    Input, integer N, the number of entries in sequence B.
!
!    Input, integer N1, N2, the minimum and maximum columns
!    of the computed score table.
!
!    Input, integer LDS, the declared upper first dimension of S.
!    LDS must be at least M.
!
!    Input, integer T(0:LDS,0:N), the backward pointer table.
!
!    Output, integer NPATH, the number of points in the matching.
!
!    Input, integer PATHI(M+N+1), PATHJ(M+N+1), contains, in the 
!    first NPATH entries, the indices of the aligned items.
!    The first entries are special marker values:
!      PATHI(1) = 0 and PATHJ(1) = 0;
!    A value of -1 for PATHI or PATHJ indicates a null partner.
!    Otherwise, A(PATHI(I)) is matched to B(PATHJ(I)).
!
  implicit none

  integer lds
  integer m
  integer n

  integer d_new
  integer d_old
  integer i
  integer i_new
  integer i_old
  integer i1
  integer ipath
  integer j
  integer j_new
  integer j_old
  integer j1
  integer m_new
  integer m_old
  integer m1
  integer m2
  integer n1
  integer n2
  integer npath
  integer pathi(m+n+1)
  integer pathj(m+n+1)
  integer t(0:lds,0:n)
  integer t_new
  integer t_old
  integer tij
  integer tij_old

  npath = 0

  if ( i1 < m1 .or. i1 > m2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_QG_BPQ - Fatal error!'
    write ( *, '(a)' ) '  The M1 <= I1 <= M2 requirement was violated.'
    write ( *, '(a,i6)' ) '  M1 = ', m1
    write ( *, '(a,i6)' ) '  I1 = ', i1
    write ( *, '(a,i6)' ) '  M2 = ', m2
    return
  else if ( j1 < n1 .or. j1 > n2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_QG_BPQ - Fatal error!'
    write ( *, '(a)' ) '  The N1 <= J1 <= N2 requirement was violated.'
    write ( *, '(a,i6)' ) '  N1 = ', n1
    write ( *, '(a,i6)' ) '  J1 = ', j1
    write ( *, '(a,i6)' ) '  N2 = ', n2
    return
  end if

  i = i1
  j = j1

  i_old = 0
  d_old = 0
  m_old = 1
  t_old = 0

  do while ( i <= m2 .and. j <= n2 )

    npath = npath + 1
    pathi(npath) = i
    pathj(npath) = j

    if ( i == m2 .and. j == n2 ) then
      exit
    end if

    tij = t(i,j)

    m_new = mod ( tij, 2 )
    tij = tij / 2
    i_new = mod ( tij, 2 )
    tij = tij / 2
    d_new = mod ( tij, 2 )
    tij = tij / 2
    t_new = tij
!
!  Try to handle ambiguous cases.
!
    if ( i_old == 1 ) then

      if ( i_new == 1 ) then
        d_new = 0
        m_new = 0
      end if

    else if ( d_old == 1 ) then

      if ( d_new == 1 ) then
        i_new = 0
        m_new = 0
      end if

    end if
!
!  If we've reached a terminator, accept it.
!
    if ( t_new == 1 ) then

      exit

    else if ( j < n2 .and. i_new == 1 ) then

      j = j + 1
      d_new = 0
      m_new = 0

    else if ( i < m2 .and. d_new == 1 ) then

      i = i + 1
      i_new = 0
      m_new = 0

    else if ( i < m2 .and. j < n2 .and. m_new == 1 ) then

      i = i + 1
      j = j + 1
      i_new = 0
      d_new = 0

    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_QG_BPQ: Unexpected situation!'
      write ( *, '(a,i6)' ) '  I = ', i
      write ( *, '(a,i6)' ) '  J = ', j
      write ( *, '(a,i6)' ) '  T(I,J) = ', t(i,j)
      stop

    end if
!
!  Copy the information.  Only one of these values is now nonzero,
!  recording which direction we actually took.
!
    i_old = i_new
    d_old = d_new
    m_old = m_new
    t_old = t_new

  end do
!
!  Mark DELETEs and INSERTs.
!
  i_new = -1
  j_new = -1

  do ipath = 1, npath

    i_old = i_new
    j_old = j_new

    i_new = pathi(ipath)
    j_new = pathj(ipath)

    if ( i_new == i_old ) then
      pathi(ipath) = -1
    else if ( j_new == j_old ) then
      pathj(ipath) = -1
    end if

  end do

  return
end
subroutine ps_qg_bsl ( b, m, m1, m2, n, n1, n2, ps_score, go, ge, &
  base, s, e, f, t, i1, j1, s_max )

!*****************************************************************************80
!
!! PS_QG_BSL: quasiglobal gap backward alignment score in linear space.
!
!  Discussion:
!
!    An affine gap penalty scheme is used, but there is no penalty for
!    the very first gap (an insertion, or a deletion, but not both), and
!    the very last one.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Kun-Mao Chao, Ross Hardison, Webb Miller,
!    Recent Developments in Linear-Space Alignment Methods: A Survey,
!    Journal of Computational Biology,
!    Volume 1, Number 4, 1994, pages 271-291.
!
!    Michael Waterman,
!    Introduction to Computational Biology,
!    Chapman and Hall, 1995.
!
!  Parameters:
!
!    Input, character B(N), the sequence to be aligned.
!
!    Input, integer M, the number of entries in the profile.
!
!    Input, integer M1, M2, the first and last rows of the
!    table to compute.
!
!    Input, integer N, the number of entries in the sequence.
!
!    Input, integer N1, N2, the first and last columns of the
!    table to compute.
!
!    Input, external PS_SCORE, the name of a function of the form
!      function ps_score ( p1, c2 )
!    which returns a real value PS_SCORE for the matching of the position
!    P1 from the profile to character C2 from the sequence B.
!
!    Input, real GO(0:M), the profile position gap open penalty.
!
!    Input, real GE(0:M), the profile position gap extend penalty.
!
!    Input, real BASE, an initial quantity to be added to certain
!    penalties.
!
!    Output, real S(0:N), the backward optimal score vector.
!    The maximum possible alignment score is in S(N1).
!
!    Output, real E(0:N), the backward final insertion 
!    score vector.
!
!    Output, real F(0:N), the backward final deletion score vector.
!
!    Output, integer T(0:N), the backward pointer vector.
!
!    Output, integer I1, J1, the location of an optimal
!    quasiglobal score.  I1 will equal M1, or J1 will equal N1.
!
!    Output, real S_MAX, the optimal quasiglobal score.
!
  implicit none

  integer, parameter :: DUNNO = 0
  integer, parameter :: MATCH = 1
  integer, parameter :: INSERT = 2
  integer, parameter :: DELETE = 4
  integer, parameter :: TERMINATE = 8

  integer m
  integer n

  character b(n)
  real base
  real diag_new
  real diag_old
  real e(0:n)
  real f(0:n)
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
  real, external :: ps_score
  real s(0:n)
  real s_max
  integer t(0:n)
!
!  The last row, I = M2.
!
  e(n2) = 0.0
  f(n2) = 0.0
  s(n2) = 0.0
  t(n2) = TERMINATE

  if ( n2-1 >= n1 ) then
    e(n2-1) = e(n2)
    f(n2-1) = f(n2) + go(m2)
    s(n2-1) = s(n2)
    t(n2-1) = TERMINATE
  end if

  do j = n2-2, n1, -1
    e(j) = e(j+1)
    f(j) = f(j+1)
    s(j) = s(j+1)
    t(j) = TERMINATE
  end do

  i1 = m2
  j1 = n1
  s_max = s(j1)
!
!  Upper rectangle.
!
  do i = m2-1, m1, -1

    diag_old = s(n2)

    if ( i == m2-1 ) then
      e(n2) = base + 2.0 * go(i) -       go(i+1)
      f(n2) = base
      s(n2) = base
      t(n2) = DELETE
    else
      e(n2) = e(n2)
      f(n2) = f(n2)
      s(n2) = s(n2)
      t(n2) = DELETE
    end if

    do j = n2-1, n1, -1
!
!  Insertion.
!
      e(j) = e(j+1) + ge(i)

      if ( s(j+1) + go(i) + ge(i) > e(j) ) then
        e(j) = s(j+1) + go(i) + ge(i)
      end if
!
!  Deletion.
!
      f(j) = f(j) + ge(i) + go(i) - go(i+1)

      if ( s(j) + go(i) + ge(i) > f(j) ) then
        f(j) = s(j) + go(i) + ge(i)
      end if
!
!  Best.
!
      diag_new = s(j)

      s(j) = diag_old + ps_score ( i+1, b(j+1) )
      t(j) = MATCH

      if ( s(j) == e(j) ) then
        t(j) = t(j) + INSERT
      else if ( s(j) < e(j) ) then
        s(j) = e(j)
        t(j) = INSERT
      end if

      if ( s(j) == f(j) ) then
        t(j) = t(j) + DELETE
      else if ( s(j) < f(j) ) then
        s(j) = f(j)
        t(j) = DELETE
      end if

      diag_old = diag_new

    end do
!
!  Check the first column, or finally, the first row, for a better score.
!  If we have two choices, prefer the location that is nearer (M1,N2).
!
    if ( i > m1 ) then
      if ( s(n1) >= s_max ) then
        i1 = i
        j1 = n1
        s_max = s(n1)
      end if
    else if ( i == m1 ) then
      do j = n1, n2
        if ( s(j) >= s_max ) then
          i1 = m1
          j1 = j
          s_max = s(j)
        end if
      end do
    end if

  end do

  return
end
subroutine ps_qg_bsq ( b, m, m1, m2, n, n1, n2, ps_score, go, ge, &
  base, lds, s, e, f, t )

!*****************************************************************************80
!
!! PS_QG_BSQ: a quasiglobal gap backward alignment score in quadratic space.
!
!  Discussion:
!
!    The routine can compute the full score table, or a sub-block.
!
!    An affine gap penalty scheme is used, but there is no penalty for
!    the very first gap (an insertion, or a deletion, but not both), and
!    the very last one.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Kun-Mao Chao, Ross Hardison, Webb Miller,
!    Recent Developments in Linear-Space Alignment Methods: A Survey,
!    Journal of Computational Biology,
!    Volume 1, Number 4, 1994, pages 271-291.
!
!    Michael Waterman,
!    Introduction to Computational Biology,
!    Chapman and Hall, 1995, page 199.
!
!  Parameters:
!
!    Input, character B(N), the sequence to be aligned.
!
!    Input, integer M, the number of entries in the profile.
!
!    Input, integer M1, M2, the first and last rows of the table
!    to compute.
!
!    Input, integer N, the number of entries in the sequence.
!
!    Input, integer N1, N2, the first and last columns of the
!    table to compute.
!
!    Input, external PS_SCORE, the name of a function of the form
!      function ps_score ( p1, c2 )
!    which returns a real value PS_SCORE for the matching of the position
!    P1 from the profile to character C2 from the sequence B.
!
!    Input, real GO(0:M), the profile position gap open penalty.
!
!    Input, real GE(0:M), the profile position gap extend penalty.
!
!    Input, real BASE, an initial quantity to be added to certain
!    penalties.
!
!    Input, integer LDS, the declared upper first dimension of S.
!    LDS must be at least M.
!
!    Output, real S(0:LDS,0:N), the backward optimal score table.
!    Entries in the M1:M2 by N1:N2 block have been computed.
!
!    Output, real E(0:LDS,0:N), the backward final insertion table.
!    Entries in the M1:M2 by N1:N2 block have been computed.
!
!    Output, real F(0:LDS,0:N), the backward final deletion table.
!    Entries in the M1:M2 by N1:N2 block have been computed.
!
!    Output, integer T(0:LDS,0:N), the backward pointer table.
!    Entries in the M1:M2 by N1:N2 block have been computed.
!
  implicit none

  integer, parameter :: DUNNO = 0
  integer, parameter :: MATCH = 1
  integer, parameter :: INSERT = 2
  integer, parameter :: DELETE = 4
  integer, parameter :: TERMINATE = 8

  integer lds
  integer m
  integer n

  character b(n)
  real base
  real e(0:lds,0:n)
  real f(0:lds,0:n)
  real ge(0:m)
  real go(0:m)
  integer i
  integer j
  integer m1
  integer m2
  integer n1
  integer n2
  real, external :: ps_score
  real s(0:lds,0:n)
  integer t(0:lds,0:n)
!
!  Lower Right corner.
!
  e(m2,n2) = 0.0
  f(m2,n2) = 0.0
  s(m2,n2) = 0.0
  t(m2,n2) = TERMINATE
!
!  Lower Left row.
!
  if ( n2-1 >= n1 ) then
    e(m2,n2-1) = e(m2,n2)
    f(m2,n2-1) = f(m2,n2) + go(m2)
    s(m2,n2-1) = 0.0
    t(m2,n2-1) = TERMINATE
  end if

  do j = n2-2, n1, -1
    e(m2,j) = e(m2,j+1)
    f(m2,j) = f(m2,j+1)
    s(m2,j) = 0.0
    t(m2,j) = TERMINATE
  end do
!
!  Upper rectangle
!
  do i = m2-1, m1, -1

    if ( i == m2-1 ) then
      e(i,n2) = base + 2.0 * go(i) -       go(i+1)
      f(i,n2) = base
      s(i,n2) = 0.0
      t(i,n2) = TERMINATE
    else
      e(i,n2) = e(i+1,n2) + 2.0 * go(i) - 2.0 * go(i+1)
      f(i,n2) = f(i+1,n2)
      s(i,n2) = 0.0
      t(i,n2) = TERMINATE
    end if

    do j = n2-1, n1, -1
!
!  Insertion.
!
      e(i,j) = e(i,j+1) + ge(i)

      if ( s(i,j+1) + go(i) + ge(i) > e(i,j) ) then
        e(i,j) = s(i,j+1) + go(i) + ge(i)
      end if
!
!  Deletion.
!
      f(i,j) = f(i+1,j) + ge(i) + go(i) - go(i+1)

      if ( s(i+1,j) + go(i) + ge(i) > f(i,j) ) then
        f(i,j) = s(i+1,j) + go(i) + ge(i)
      end if
!
!  Best.
!
      s(i,j) = s(i+1,j+1) + ps_score ( i+1, b(j+1) )
      t(i,j) = MATCH

      if ( s(i,j) == e(i,j) ) then
        t(i,j) = t(i,j) + INSERT
      else if ( s(i,j) < e(i,j) ) then
        s(i,j) = e(i,j)
        t(i,j) = INSERT
      end if

       if ( s(i,j) == f(i,j) ) then
        t(i,j) = t(i,j) + DELETE
       else if ( s(i,j) < f(i,j) ) then
        s(i,j) = f(i,j)
        t(i,j) = DELETE
      end if

    end do
  end do

  return
end
subroutine ps_qg_foq ( m, m1, m2, n, n1, n2, lds, s, i2, j2 )

!*****************************************************************************80
!
!! PS_QG_FOQ: the forward endpoint of a quasiglobal optimal local alignment.
!
!  Discussion:
!
!    The algorithm examines the last column and row of the quadratic score
!    table, and returns I2, J2, the location of the maximum value.
!
!    Actually, for generality, the algorithm examines row M2 and column N2.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of entries in sequence A.
!
!    Input, integer M1, M2, the minimum and maximum rows 
!    to consider.
!
!    Input, integer N, the number of entries in sequence B.
!
!    Input, integer N1, N2, the minimum and maximum columns 
!    to consider.
!
!    Input, integer LDS, the declared upper first dimension of S.
!    LDS must be at least M.
!
!    Output, real S(0:LDS,0:N), the forward optimal score table.
!    Entries in the M1:M2 by N1:N2 block have been computed.
!
!    Output, integer I2, J2, the beginning of an optimal
!    quasiglobal alignment.  I2 is equal to M2, or J2 is equal to N2.
!
  implicit none

  integer lds
  integer n

  integer i
  integer i2
  integer j
  integer j2
  integer m
  integer m1
  integer m2
  integer n1
  integer n2
  integer s(0:lds,0:n)

  i2 = m2
  j2 = n1

  do j = n1+1, n2
    if ( s(m2,j) > s(i2,j2) ) then
      i2 = m2
      j2 = j
    end if
  end do

  do i = m2-1, m1, -1
    if ( s(i,n2) > s(i2,j2) ) then
      i2 = i
      j2 = n2
    end if
  end do

  return
end
subroutine ps_qg_fpq ( i2, j2, m, m1, m2, n, n1, n2, lds, t, npath, pathi, &
  pathj )

!*****************************************************************************80
!
!! PS_QG_FPQ: a quasiglobal gap forward alignment path in quadratic space.
!
!  Discussion:
!
!    An effort has been made to handle the ambiguous case, where
!    more than one optimal path enters a cell.  In such a case,
!    the code tries to take a delete out if there was a delete in,
!    or an insert out if there was an insert in, since the optimal
!    score calculation includes a penalty for gap opening.
!
!    The routine is called "quadratic" because it uses an M by N array
!    to do the alignment.
!
!    The score table must have been computed before this routine is called.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Kun-Mao Chao, Ross Hardison, Webb Miller,
!    Recent Developments in Linear-Space Alignment Methods: A Survey,
!    Journal of Computational Biology,
!    Volume 1, Number 4, 1994, pages 271-291.
!
!  Parameters:
!
!    Input, integer I2, J2, the starting point of the path.
!    M1 <= I2 <= M2, N1 <= J2 <= N2.
!
!    Input, integer M, the number of entries in the profile.
!
!    Input, integer M1, M2, the first and last rows of the table
!    to compute.
!
!    Input, integer N, the number of entries in the sequence.
!
!    Input, integer N1, N2, the first and last columns of the
!    table to compute.
!
!    Input, integer LDS, the declared upper first dimension of T.
!    LDS must be at least M.
!
!    Input, integer T(0:LDS,0:N), the forward pointer table.
!
!    Output, integer NPATH, the number of points in the matching.
!
!    Input, integer PATHI(M+N+1), PATHJ(M+N+1), contains, in the 
!    first NPATH entries, the indices of the aligned items.
!    The first entries are special marker values:
!      PATHI(1) = 0 and PATHJ(1) = 0;
!    A value of -1 for PATHI or PATHJ indicates a null partner.
!    Otherwise, A(PATHI(I)) is matched to B(PATHJ(I)).
!
  implicit none

  integer lds
  integer m
  integer n

  integer d_new
  integer d_old
  integer i
  integer i_new
  integer i_old
  integer i2
  integer ipath
  integer j
  integer j_new
  integer j_old
  integer j2
  integer m_new
  integer m_old
  integer m1
  integer m2
  integer n1
  integer n2
  integer npath
  integer pathi(m+n+1)
  integer pathj(m+n+1)
  integer t(0:lds,0:n)
  integer t_new
  integer t_old
  integer tij
  integer tij_old

  npath = 0

  if ( i2 < m1 .or. i2 > m2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_QG_FPQ - Fatal error!'
    write ( *, '(a)' ) '  The M1 <= I2 <= M2 requirement was violated.'
    write ( *, '(a,i6)' ) '  M1 = ', m1
    write ( *, '(a,i6)' ) '  I2 = ', i2
    write ( *, '(a,i6)' ) '  M2 = ', m2
    return
  else if ( j2 < n1 .or. j2 > n2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_QG_FPQ - Fatal error!'
    write ( *, '(a)' ) '  The N1 <= J2 <= N2 requirement was violated.'
    write ( *, '(a,i6)' ) '  N1 = ', n1
    write ( *, '(a,i6)' ) '  J2 = ', j2
    write ( *, '(a,i6)' ) '  N2 = ', n2
    return
  end if

  i = i2
  j = j2

  i_old = 0
  d_old = 0
  m_old = 1
  t_old = 0

  do while ( i >= m1 .and. j >= n1 )

    npath = npath + 1
    pathi(npath) = i
    pathj(npath) = j

    if ( i == m1 .and. j == n1 ) then
      exit
    end if

    tij = t(i,j)

    m_new = mod ( tij, 2 )
    tij = tij / 2
    i_new = mod ( tij, 2 )
    tij = tij / 2
    d_new = mod ( tij, 2 )
    tij = tij / 2
    t_new = tij
!
!  Try to handle ambiguous cases.
!
    if ( i_old == 1 ) then

      if ( i_new == 1 ) then
        d_new = 0
        m_new = 0
      end if

    else if ( d_old == 1 ) then

      if ( d_new == 1 ) then
        i_new = 0
        m_new = 0
      end if

    end if
!
!  If we've reached a terminator, we're happy.  Jump (out) for joy.
!
    if ( t_new == 1 ) then

      exit

    else if ( j > n1 .and. i_new == 1 ) then

      j = j - 1
      d_new = 0
      m_new = 0

    else if ( i > m1 .and. d_new == 1 ) then

      i = i - 1
      i_new = 0
      m_new = 0

    else if ( i > m1 .and. j > n1 .and. m_new == 1 ) then

      i = i - 1
      j = j - 1
      i_new = 0
      d_new = 0

    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_QG_FPQ: Unexpected situation!'
      write ( *, '(a,i6)' ) '  I = ', i
      write ( *, '(a,i6)' ) '  J = ', j
      write ( *, '(a,i6)' ) '  T(I,J) = ', t(i,j)
      stop

    end if
!
!  Copy the information.  Only one of these values is now nonzero,
!  recording which direction we actually took.
!
    i_old = i_new
    d_old = d_new
    m_old = m_new
    t_old = t_new

  end do
!
!  Put the path into proper order.
!
  call i4vec_reverse ( npath, pathi )
  call i4vec_reverse ( npath, pathj )
!
!  Mark DELETEs and INSERTs.
!
  i_new = -1
  j_new = -1

  do ipath = 1, npath

    i_old = i_new
    j_old = j_new

    i_new = pathi(ipath)
    j_new = pathj(ipath)

    if ( i_new == i_old ) then
      pathi(ipath) = -1
    else if ( j_new == j_old ) then
      pathj(ipath) = -1
    end if

  end do

  return
end
subroutine ps_qg_fsl ( b, m, m1, m2, n, n1, n2, ps_score, go, ge, &
  base, s, e, f, t, i2, j2, s_max )

!*****************************************************************************80
!
!! PS_QG_FSL: quasiglobal gap forward alignment score in linear space.
!
!  Discussion:
!
!    An affine gap penalty scheme is used, but there is no penalty for
!    the very first gap (an insertion, or a deletion, but not both), and
!    the very last one.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Kun-Mao Chao, Ross Hardison, Webb Miller,
!    Recent Developments in Linear-Space Alignment Methods: A Survey,
!    Journal of Computational Biology,
!    Volume 1, Number 4, 1994, pages 271-291.
!
!    Michael Waterman,
!    Introduction to Computational Biology,
!    Chapman and Hall, 1995.
!
!  Parameters:
!
!    Input, character B(N), the sequence to be aligned.
!
!    Input, integer M, the number of entries in the profile.
!
!    Input, integer M1, M2, the first and last rows of the
!    table to compute.
!
!    Input, integer N, the number of entries in the sequence.
!
!    Input, integer N1, N2, the first and last columns of the
!    table to compute.
!
!    Input, external PS_SCORE, the name of a function of the form
!      function ps_score ( p1, c2 )
!    which returns a real value PS_SCORE for the matching of the position
!    P1 from the profile to character C2 from the sequence B.
!
!    Input, real GO(0:M), the profile position gap open penalty.
!
!    Input, real GE(0:M), the profile position gap extend penalty.
!
!    Input, real BASE, an initial quantity to be added to certain
!    penalties.
!
!    Output, real S(0:N), the forward optimal score vector.
!    The maximum possible alignment score is in S(N2).
!
!    Output, real E(0:N), the forward final insertion score vector.
!
!    Output, real F(0:N), the forward final deletion score vector.
!
!    Output, integer T(0:N), the forward pointer vector.
!
!    Output, integer I2, J2, the location of an optimal
!    quasiglobal score.  I2 will equal M2, or J2 will equal N2.
!
!    Output, real S_MAX, the optimal quasiglobal score.
!
  implicit none

  integer, parameter :: DUNNO = 0
  integer, parameter :: MATCH = 1
  integer, parameter :: INSERT = 2
  integer, parameter :: DELETE = 4
  integer, parameter :: TERMINATE = 8

  integer m
  integer n

  character b(n)
  real base
  real diag_new
  real diag_old
  real e(0:n)
  real f(0:n)
  real ge(0:m)
  real go(0:m)
  real horz
  integer i
  integer i2
  integer j
  integer j2
  integer m1
  integer m2
  integer n1
  integer n2
  real, external :: ps_score
  real s(0:n)
  real s_max
  integer t(0:n)
!
!  The first row, I = M1.
!
  e(n1) = 0.0
  f(n1) = 0.0
  s(n1) = 0.0
  t(n1) = TERMINATE

  if ( n1+1 <= n2 ) then
    e(n1+1) = e(n1)
    f(n1+1) = f(n1) + go(m1)
    s(n1+1) = 0.0
    t(n1+1) = TERMINATE
  end if

  do j = n1+2, n2
    e(j) = e(j-1)
    f(j) = f(j-1)
    s(j) = 0.0
    t(j) = TERMINATE
  end do

  i2 = m1
  j2 = n2
  s_max = s(j2)
!
!  Subsequent rows.
!
  do i = m1+1, m2

    diag_old = s(n1)

    if ( i == m1+1 ) then
      e(n1) = base + go(i-1)
      f(n1) = base
      s(n1) = 0.0
      t(n1) = TERMINATE
    else
      e(n1) = e(n1)
      f(n1) = f(n1)
      s(n1) = 0.0
      t(n1) = TERMINATE
    end if

    do j = n1+1, n2
!
!  Insertion
!
      e(j) = e(j-1) + ge(i)

      if ( s(j-1) + go(i) + ge(i) > e(j) ) then
        e(j) = s(j-1) + go(i) + ge(i)
      end if
!
!  Deletion.
!
      f(j) = f(j) + ge(i-1)

      if ( s(j) + go(i-1) + ge(i-1) > f(j) ) then
        f(j) = s(j) + go(i-1) + ge(i-1)
      end if
!
!  Best.
!
      diag_new = s(j)

      s(j) = diag_old + ps_score ( i, b(j) )
      t(j) = MATCH

      if ( s(j) == e(j) ) then
        t(j) = t(j) + INSERT
      else if ( s(j) < e(j) ) then
        s(j) = e(j)
        t(j) = INSERT
      end if

      if ( s(j) == f(j) ) then
        t(j) = t(j) + DELETE
      else if ( s(j) < f(j) ) then
        s(j) = f(j)
        t(j) = DELETE
      end if

      diag_old = diag_new

    end do
!
!  Check the last column, or finally, the last row, for a better score.
!  If we have two choices, prefer the location that is nearer (M1,N2).
!
    if ( i < m2 ) then
      if ( s(n2) > s_max ) then
        i2 = i
        j2 = n2
        s_max = s(j2)
      end if
    else if ( i == m2 ) then
      do j = n2, n1, -1
        if ( s(j) > s_max ) then
          i2 = i
          j2 = j
          s_max = s(j2)
        end if
      end do
    end if

  end do

  return
end
subroutine ps_qg_fsq ( b, m, m1, m2, n, n1, n2, ps_score, go, ge, &
  base, lds, s, e, f, t )

!*****************************************************************************80
!
!! PS_QG_FSQ: a quasiglobal gap forward alignment score in quadratic space.
!
!  Discussion:
!
!    An affine gap penalty scheme is used, but there is no penalty for
!    the very first gap (an insertion, or a deletion, but not both), and
!    the very last one.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Kun-Mao Chao, Ross Hardison, Webb Miller,
!    Recent Developments in Linear-Space Alignment Methods: A Survey,
!    Journal of Computational Biology,
!    Volume 1, Number 4, 1994, pages 271-291.
!
!    Michael Waterman,
!    Introduction to Computational Biology,
!    Chapman and Hall, 1995, page 199.
!
!  Parameters:
!
!    Input, character B(N), the sequence to be aligned.
!
!    Input, integer M, the number of entries in the profile.
!
!    Input, integer M1, M2, the first and last rows of the
!    table to compute.
!
!    Input, integer N, the number of entries in the sequence.
!
!    Input, integer N1, N2, the first and last columns of the
!    table to compute.
!
!    Input, external PS_SCORE, the name of a function of the form
!      function ps_score ( p1, c2 )
!    which returns a real value PS_SCORE for the matching of the position
!    P1 from the profile to character C2 from the sequence B.
!
!    Input, real GO(0:M), the profile position gap open penalty.
!
!    Input, real GE(0:M), the profile position gap extend penalty.
!
!    Input, real BASE, an initial quantity to be added to certain
!    penalties.
!
!    Input, integer LDS, the declared upper first dimension of S.
!    LDS must be at least M.
!
!    Output, real S(0:LDS,0:N), the forward optimal score table.
!    Entries in the M1:M2 by N1:N2 block have been computed.
!
!    Output, real E(0:LDS,0:N), the forward final insertion table.
!    Entries in the M1:M2 by N1:N2 block have been computed.
!
!    Output, real F(0:LDS,0:N), the forward final deletion table.
!    Entries in the M1:M2 by N1:N2 block have been computed.
!
!    Output, integer T(0:LDS,0:N), the forward pointer table.
!    Entries in the M1:M2 by N1:N2 block have been computed.
!
  implicit none

  integer, parameter :: DUNNO = 0
  integer, parameter :: MATCH = 1
  integer, parameter :: INSERT = 2
  integer, parameter :: DELETE = 4
  integer, parameter :: TERMINATE = 8

  integer lds
  integer m
  integer n

  character b(n)
  real base
  real e(0:lds,0:n)
  real f(0:lds,0:n)
  real ge(0:m)
  real go(0:m)
  integer i
  integer j
  integer m1
  integer m2
  integer n1
  integer n2
  real, external :: ps_score
  real s(0:lds,0:n)
  integer t(0:lds,0:n)
!
!  Upper Left corner.
!
  e(m1,n1) = 0.0
  f(m1,n1) = 0.0
  s(m1,n1) = 0.0
  t(m1,n1) = TERMINATE
!
!  Upper Right row.
!
  if ( n1+1 <= n2 ) then
    e(m1,n1+1) = e(m1,n1)
    f(m1,n1+1) = f(m1,n1) + go(m1)
    s(m1,n1+1) = 0.0
    t(m1,n1+1) = TERMINATE
  end if

  do j = n1+2, n2
    e(m1,j) = e(m1,j-1)
    f(m1,j) = f(m1,j-1)
    s(m1,j) = 0.0
    t(m1,j) = TERMINATE
  end do
!
!  Lower rectangle.
!
  do i = m1+1, m2

    if ( i == m1+1 ) then
      e(i,n1) = base + go(i-1)
      f(i,n1) = base
      s(i,n1) = 0.0
      t(i,n1) = TERMINATE
    else
      e(i,n1) = e(i-1,n1)
      f(i,n1) = f(i-1,n1)
      s(i,n1) = 0.0
      t(i,n1) = TERMINATE
    end if

    do j = n1+1, n2
!
!  A matching that ends with an insertion could either be
!  *) an extension of an ongoing insertion
!  *) the best matching to the left, followed by a new insertion.
!
      e(i,j) = e(i,j-1) + ge(i)

      if ( s(i,j-1) + go(i) + ge(i) > e(i,j) ) then
        e(i,j) = s(i,j-1) + go(i) + ge(i)
      end if
!
!  A matching that ends with a deletion could either be
!  *) an extension of an ongoing deletion;
!  *) the best matching above, followed by a new deletion.
!
      f(i,j) = f(i-1,j) + ge(i-1)

      if ( s(i-1,j) + go(i-1) + ge(i-1) > f(i,j) ) then
        f(i,j) = s(i-1,j) + go(i-1) + ge(i-1)
      end if
!
!  Best.
!
      s(i,j) = s(i-1,j-1) + ps_score ( i, b(j) )
      t(i,j) = MATCH

      if ( s(i,j) == e(i,j) ) then
        t(i,j) = t(i,j) + INSERT
      else if ( s(i,j) < e(i,j) ) then
        s(i,j) = e(i,j)
        t(i,j) = INSERT
      end if

      if ( s(i,j) == f(i,j) ) then
        t(i,j) = t(i,j) + DELETE
      else if ( s(i,j) < f(i,j) ) then
        s(i,j) = f(i,j)
        t(i,j) = DELETE
      end if

    end do
  end do

  return
end
subroutine ps_qg_match_print ( b, m, n, npath, pathi, pathj, ps_score, &
  go, ge )

!*****************************************************************************80
!
!! PS_QG_MATCH_PRINT prints a quasiglobal gap alignment.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Kun-Mao Chao, Ross Hardison, Webb Miller,
!    Recent Developments in Linear-Space Alignment Methods: A Survey,
!    Journal of Computational Biology,
!    Volume 1, Number 4, 1994, pages 271-291.
!
!  Parameters:
!
!    Input, character B(N), the sequence to be aligned.
!
!    Input, integer M, the number of entries in the profile.
!
!    Input, integer N, the number of entries in sequence B.
!
!    Input, integer NPATH, the number of points in the matching.
!
!    Input, integer PATHI(M+N+1), PATHJ(M+N+1), contains, in the
!    first NPATH entries, the indices of the most recently matched items
!    from sequences A and B.
!
!    Input, external PS_SCORE, the name of a function of the form
!      function ps_score ( p1, c2 )
!    which returns a real value PS_SCORE for the matching of the position
!    P1 from the profile to character C2 from the sequence B.
!
!    Input, real GO(0:M), the profile position gap open penalty.
!
!    Input, real GE(0:M), the profile position gap extend penalty.
!
  implicit none

  integer m
  integer n

  character b(n)
  character ( len = 3 ) c1
  character c3
  character c4
  character ( len = 3 ) c5
  integer col
  real ge(0:m)
  real go(0:m)
  integer i
  integer i_old
  real inc
  integer ipath
  integer j
  integer j_old
  integer npath
  integer p1
  integer pathi(m+n+1)
  integer pathj(m+n+1)
  real, external :: ps_score
  real psum
  integer row

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Profile/sequence matching,'
  write ( *, '(a)' ) 'Affine gap penalty:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' #     B    #      Increm     Score'
  write ( *, '(a)' ) ' '

  psum = 0.0
  ipath = 1

  c1 = ' '
  c3 = ' '
  c4 = ' '
  c5 = ' '

  write ( *, '(a3,2x,a1,2x,a1,2x,a3,2x,10x,f10.2)' ) &
    c1, c3, c4, c5, psum

  i_old = pathi(1)
  j_old = pathj(1)

  row = 0
  col = 0

  do ipath = 2, npath

    i = pathi(ipath)
    j = pathj(ipath)

    if ( i == -1 ) then

      col = j

      c1 = '  |'
      c3 = ' '
      c4 = b(j)
      write ( c5, '(i3)' ) j

      if ( i_old /= -1 ) then
        inc = go(row) + ge(row)
      else
        inc = ge(row)
      end if
!
!  When moving "vertically", the GAP_OPEN penalty and GAP_EXTEND penalty
!  are gotten from the previous (lower index) row.
!
    else if ( j == -1 ) then

      row = i

      write ( c1, '(i3)' ) i
      c3 = ' '
      c4 = '|'
      c5 = ' '

      if ( j_old /= -1 ) then
        inc = go(row-1) + ge(row-1)
      else
        inc = ge(row-1)
      end if

    else

      row = i
      col = j

      write ( c1, '(i3)' ) i
      c3 = '-'
      c4 = b(j)
      write ( c5, '(i3)' ) j

      inc = ps_score ( i, b(j) )

    end if

    i_old = i
    j_old = j

    psum = psum + inc

    write ( *, '(a3,2x,a1,2x,a1,2x,a3,2x,2f10.2)' ) &
      c1, c3, c4, c5, inc, psum

  end do

  return
end
subroutine ps_qg_rpl ( b, m, m1, m2, n, n1, n2, ps_score, go, ge, &
  npath, pathi, pathj )

!*****************************************************************************80
!
!! PS_QG_RPL: quasiglobal gap recursive alignment path in linear space.
!
!  Discussion:
!
!    The routine is called "linear" because it uses a few vectors of
!    dimension N to determine the alignment.
!
!    An affine gap penalty scheme is used, but there is no penalty for
!    the very first gap (an insertion, or a deletion, but not both), and
!    the very last one.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Eugene Myers and Webb Miller,
!    Optimal Alignments in Linear Space,
!    CABIOS, volume 4, number 1, 1988, pages 11-17.
!
!    Kun-Mao Chao, Ross Hardison, Webb Miller,
!    Recent Developments in Linear-Space Alignment Methods: A Survey,
!    Journal of Computational Biology,
!    Volume 1, Number 4, 1994, pages 271-291.
!
!  Parameters:
!
!    Input, character B(N), the sequence to be aligned.
!
!    Input, integer M, the number of entries in the profile.
!
!    Input, integer M1, M2, the minimum and maximum entries of
!    A to align.  For a full alignment, use M1 = 0, M2 = M.
!
!    Input, integer N, the number of entries in the sequence.
!
!    Input, integer N1, N2, the minimum and maximum entries of B
!    to align.  For a full alignment, use N1 = 0, N2 = N.
!
!    Input, external PS_SCORE, the name of a function of the form
!      function ps_score ( p1, c2 )
!    which returns a real value PS_SCORE for the matching of the position
!    P1 from the profile to character C2 from the sequence B.
!
!    Input, real GO(0:M), the profile position gap open penalty.
!
!    Input, real GE(0:M), the profile position gap extend penalty.
!
!    Output, integer NPATH, the number of points in the matching.
!
!    Output, integer PATHI(M+N+1), PATHJ(M+N+1), contains, in the 
!    first NPATH entries, the indices of the matched items from the profile
!    and sequence B.
!
  implicit none

  integer, parameter :: DUNNO = 0
  integer, parameter :: MATCH = 1
  integer, parameter :: INSERT = 2
  integer, parameter :: DELETE = 4
  integer, parameter :: TERMINATE = 8
  integer, parameter :: stack_max = 200

  integer m
  integer m1
  integer m2
  integer n

  character b(n)
  real base
  real eb(0:n)
  real ef(0:n)
  real fb(0:n)
  real ff(0:n)
  real ge(0:m)
  real go(0:m)
  integer i
  integer i1
  integer i1sub
  integer i2
  integer i2sub
  integer i3sub
  integer j
  integer j1
  integer j1sub
  integer j1type
  integer j2
  integer j2sub
  integer j2type
  integer j3sub
  integer j3type
  integer ja
  integer jb
  integer n1
  integer n2
  integer nb1
  integer nb2
  integer npath
  integer pathi(m+n+1)
  integer pathj(m+n+1)
  real, external :: ps_score
  real s_max
  real sb(0:n)
  real sf(0:n)
  integer stack1(6,stack_max)
  real stack2(2,stack_max)
  integer :: stack_num
  real t1
  real t2
  real t3
  integer tb(0:n)
  integer tf(0:n)
  real x
  real y
  real y1
  real y2
!
!  Determine the endpoints of the quasiglobal alignment.
!
  base = 0.0

  call ps_qg_bsl ( b, m, m1, m2, n, n1, n2, ps_score, &
    go, ge, base, sb, eb, fb, tb, i1, j1, s_max )

  base = 0.0

  call ps_qg_fsl ( b, m, m1, m2, n, n1, n2, ps_score, &
    go, ge, base, sf, ef, ff, tf, i2, j2, s_max )

  if ( i1 <= i2 .and. j1 <= j2 ) then
    i1sub = i1
    j1sub = j1
    i2sub = i2
    j2sub = j2
  else if ( i1 >= i2 .and. j1 >= j2 ) then
    i1sub = i2
    j1sub = j2
    i2sub = i1
    j2sub = j1
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_QG_RPL - Fatal error!'
    write ( *, '(a)' ) &
      '  The optimal endpoints don''t have the right orientation!'
    write ( *, '(a,i6)' ) '  I1 = ', i1
    write ( *, '(a,i6)' ) '  J1 = ', j1
    write ( *, '(a,i6)' ) '  I2 = ', i2
    write ( *, '(a,i6)' ) '  J2 = ', j2
    write ( *, '(a,g14.6)' ) '  s_max = ', s_max
    stop
  end if

  j1type = MATCH
  t1 = go(i1sub)

  j2type = MATCH
  t2 = go(i2sub)
!
!  Put the two endpoints on the path.
!
  npath = 1
  pathi(npath) = i1sub
  pathj(npath) = j1sub

  npath = 2
  pathi(npath) = i2sub
  pathj(npath) = j2sub
!
!  Put the initial problem on the stack.
!
  stack_num = 0

  call ps_qg_rpl_push ( i1sub, j1sub, j1type, t1, i2sub, j2sub, j2type, t2, &
    stack1, stack2, stack_num, stack_max )

  do while ( stack_num > 0 )
!
!  Pop the next problem off the stack.
!
    call ps_qg_rpl_pop ( i1sub, j1sub, j1type, t1, i3sub, j3sub, j3type, t3, &
      stack1, stack2, stack_num, stack_max )
!
!  Refuse to handle improperly described subregions.
!
    if ( i1sub > i3sub .or. j1sub > j3sub ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_QG_RPL - Fatal error!'
      write ( *, '(a)' ) '  The indices describing a subregion have the'
      write ( *, '(a)' ) '  wrong sense.'
      write ( *, '(a,i6)' ) '  I1 = ', i1sub
      write ( *, '(a,i6)' ) '  I3 = ', i3sub
      write ( *, '(a,i6)' ) '  J1 = ', j1sub
      write ( *, '(a,i6)' ) '  J3 = ', j3sub
      stop
!
!  Null regions require no processing.
!
    else if ( i1sub == i3sub .and. j1sub == j3sub ) then
!
!  A vertical strip is easy.
!
    else if ( j1sub == j3sub ) then

      do i = i1sub+1, i3sub-1
        npath = npath + 1
        pathi(npath) = i
        pathj(npath) = j1sub
      end do
!
!  A horizontal strip is easy.
!
    else if ( i1sub == i3sub ) then

      do j = j1sub+1, j3sub-1
        npath = npath + 1
        pathi(npath) = i1sub
        pathj(npath) = j
      end do
!
!  For the case where the uncertainty region is two units high, the path
!  can be described as going from J1SUB to JA in row I1SUB, and then
!  from JB in row I3SUB to J3SUB.
!
!  We need to know:
!
!    Does the incoming path to (I1SUB,J1SUB) represent
!      a match,
!      a deletion of A's,
!      an insertion of B's?
!
!    Does the outgoing path from (I3SUB,J3SUB) represent
!      a match,
!      a deletion of A's,
!      an insertion of B's.
!
!    When done, we have to record that the incoming and outgoing paths
!    to (I2SUB,J2SUB) represent:
!      a match,
!      a deletion of A's,
!      an insertion of B's.
!
    else if ( i3sub == i1sub + 1 ) then

      x = - huge ( x )
      ja = 0
      jb = 0
!
!  A: Cost of the path:
!
!    X X X ... X X . ... . . .
!    . . . ... . X X ... X X X
!
!  Note that, when scoring "vertical" gaps (deletions), the GAP_OPEN penalty
!  is computed for the first (lowest index) row, and the GAP_EXTEND penalty
!  going from row I1 to I1+1 is computed using the data associated with
!  row I1.  (This is going to be hell to track down.)
!
      do j = j1sub, j3sub

        nb1 = j - j1sub
        nb2 = j3sub - j

        y = 0.0
!
!  Cost of (possibly null) insertion #1.
!
        if ( nb1 > 0 ) then
          if ( j1type /= INSERT ) then
            y = y + go(i1sub)
          end if
          y = y + ge(i1sub) * nb1
        end if
!
!  Cost of deletion of one position.
!
        y =  y + go(i1sub) + ge(i1sub)
!
!  Cost of (possibly null) insertion #2.
!
        if ( nb2 > 0 ) then
          if ( j3type /= INSERT ) then
            y = y + go(i3sub)
          end if
          y = y + ge(i3sub) * nb2
        end if

        if ( y > x ) then
          x = y
          ja = j
          jb = j
        end if

      end do
!
!  B: Cost of the path:
!
!    X X X ... X . ... . . .
!    . . . ... . X ... X X X
!
      do j = j1sub, j3sub - 1

        nb1 = j - j1sub
        nb2 = j3sub - j - 1

        y = 0.0
!
!  Cost of (possibly null) insertion #1.
!
        if ( nb1 > 0 ) then
          if ( j1type /= INSERT ) then
            y = y + go(i1sub)
          end if
          y = y + ge(i1sub) * nb1
        end if
!
!  Cost of match.
!
        y = y + ps_score ( i3sub, b(j+1) )
!
!  Cost of (possibly null) insertion #2.
!
        if ( nb2 > 0 ) then
          if ( j3type /= INSERT ) then
            y = y + go(i3sub)
          end if
          y = y + ge(i3sub) * nb2
        end if

        if ( y > x ) then
          x = y
          ja = j
          jb = j + 1
        end if

      end do
!
!  Now fill in the path.
!
      do j = j1sub + 1, ja
        npath = npath + 1
        pathi(npath) = i1sub
        pathj(npath) = j
      end do

      do j = jb, j3sub - 1
        npath = npath + 1
        pathi(npath) = i3sub
        pathj(npath) = j
      end do

    else

      i2sub = ( i1sub + i3sub ) / 2

      base = t1

      call ps_gg_fsl ( b, m, i1sub, i2sub, n, j1sub, j3sub, ps_score, &
        go, ge, base, sf, ef, ff, tf )

      base = t3

      call ps_gg_bsl ( b, m, i2sub, i3sub, n, j1sub, j3sub, ps_score, &
        go, ge, base, sb, eb, fb, tb )
!
!  Find J2SUB, the value of J between J1SUB and J3SUB that maximizes
!  SF(J)+SB(J) or FF(J)+FB(J)-GAP_OPEN.
!
      j = j1sub
      y1 = sf(j) + sb(j)

      x = y1
      j2sub = j
      j2type = MATCH
!
!  We subtract GO(I2SUB) here because FB(J) was "paying" this penalty,
!  which is no longer required, since FF(J) has already paid a gap opening fee.
!
      y2 = ff(j) + fb(j) - go(i2sub)

      if ( x < y2 ) then
        x = y2
        j2sub = j
        j2type = DELETE
      end if

      do j = j1sub+1, j3sub

        y1 = sf(j) + sb(j)

        if ( y1 > x ) then
          x = y1
          j2sub = j
          j2type = MATCH
        end if

        y2 = ff(j) + fb(j) - go(i2sub)

        if ( y2 > x ) then
          x = y2
          j2sub = j
          j2type = DELETE
        end if

      end do

      npath = npath + 1
      pathi(npath) = i2sub
      pathj(npath) = j2sub

      if ( j2type == MATCH ) then

        t2 = go(i2sub)

        call ps_qg_rpl_push ( i1sub, j1sub, j1type, t1, i2sub, j2sub, j2type, &
          t2, stack1, stack2, stack_num, stack_max )

        call ps_qg_rpl_push ( i2sub, j2sub, j2type, t2, i3sub, j3sub, j3type, &
          t3, stack1, stack2, stack_num, stack_max )

      else if ( j2type == DELETE ) then

        if ( ( i1sub < i2sub-1 .or. j1sub < j2sub ) .and. &
             ( i2sub+1 < i3sub .or. j2sub < j3sub ) ) then

          npath = npath + 1
          pathi(npath) = i2sub - 1
          pathj(npath) = j2sub

          npath = npath + 1
          pathi(npath) = i2sub + 1
          pathj(npath) = j2sub

          t2 = 0.0

          call ps_qg_rpl_push ( i1sub, j1sub, j1type, t1, i2sub-1, j2sub, &
            j2type, t2, stack1, stack2, stack_num, stack_max )

          call ps_qg_rpl_push ( i2sub+1, j2sub, j2type, t2, i3sub, j3sub, &
            j3type, t3, stack1, stack2, stack_num, stack_max )

        else if ( i2sub+1 < i3sub .or. j2sub < j3sub ) then

          npath = npath + 1
          pathi(npath) = i2sub + 1
          pathj(npath) = j2sub

          t2 = t1

          call ps_qg_rpl_push ( i2sub+1, j2sub, j2type, t2, i3sub, j3sub, &
            j3type, t3, stack1, stack2, stack_num, stack_max )

        else if ( i1sub < i2sub-1 .or. j1sub < j2sub ) then

          npath = npath + 1
          pathi(npath) = i2sub - 1
          pathj(npath) = j2sub

          t2 = t3

          call ps_qg_rpl_push ( i1sub, j1sub, j1type, t1, i2sub-1, j2sub, &
           j2type, t2, stack1, stack2, stack_num, stack_max )

        end if

      end if

    end if

  end do
!
!  When the stack is empty, sort the path indices.
!
  call i4vec2_sort_a ( npath, pathi, pathj )
!
!  Now go through and mark gaps.
!
  do i = npath, 2, -1
    if ( pathi(i) == pathi(i-1) ) then
      pathi(i) = -1
    else if ( pathj(i) == pathj(i-1) ) then
      pathj(i) = -1
    end if
  end do

  return
end
subroutine ps_qg_rpl_pop ( i1sub, j1sub, j1type, t1, i2sub, j2sub, j2type, &
  t2, stack1, stack2, stack_num, stack_max )

!*****************************************************************************80
!
!! PS_QG_RPL_POP pops the data describing a subproblem off of the stack.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer I1SUB, J1SUB, J1TYPE, real T1, 
!    the row and column of the first cell, the type of the match there, and the 
!    appropriate base.
!
!    Output, integer I2SUB, J2SUB, J2TYPE, real T2, 
!    the row and column of the second cell, the type of the match there, and 
!    the appropriate base.
!
!    Input, integer STACK1(6,STACK_MAX), 
!    real STACK2(2,STACK_MAX), two arrays in which stack data 
!    is stored.
!
!    Input/output, integer STACK_NUM, a pointer to the most recent
!    item on the stack.
!
!    Input, integer STACK_MAX, the maximum number of items 
!    in the stack.
!
  implicit none

  integer stack_max

  integer i1sub
  integer i2sub
  integer j1sub
  integer j1type
  integer j2sub
  integer j2type
  integer stack1(6,stack_max)
  real stack2(2,stack_max)
  integer stack_num
  real t1
  real t2

  if ( stack_num < 1 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_QG_RPL_POP - Fatal error!'
    write ( *, '(a)' ) '  No more data to pop.'
    stop

  end if

  i1sub  = stack1(1,stack_num)
  i2sub  = stack1(2,stack_num)
  j1sub  = stack1(3,stack_num)
  j1type = stack1(4,stack_num)
  j2sub  = stack1(5,stack_num)
  j2type = stack1(6,stack_num)

  t1     = stack2(1,stack_num)
  t2     = stack2(2,stack_num)

  stack_num = stack_num - 1

  return
end
subroutine ps_qg_rpl_push ( i1sub, j1sub, j1type, t1, i2sub, j2sub, j2type, &
  t2, stack1, stack2, stack_num, stack_max )

!*****************************************************************************80
!
!! PS_QG_RPL_PUSH pushes the data describing a subproblem onto the stack.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I1SUB, J1SUB, J1TYPE, real T1, 
!    the row and column of the first cell, the type of the match there, and the 
!    appropriate base.
!
!    Input, integer I2SUB, J2SUB, J2TYPE, real T2, 
!    the row and column of the second cell, the type of the match there, and 
!    the appropriate base.
!
!    Input/output, integer STACK1(6,STACK_MAX), 
!    real STACK2(2,STACK_MAX), two arrays in which stack data 
!    is stored.
!
!    Input/output, integer STACK_NUM, a pointer to the most
!    recent item on the stack.
!
!    Input, integer STACK_MAX, the maximum number of items 
!    in the stack.
!
  implicit none

  integer stack_max

  integer i1sub
  integer i2sub
  integer j1sub
  integer j1type
  integer j2sub
  integer j2type
  integer stack1(6,stack_max)
  real stack2(2,stack_max)
  integer stack_num
  real t1
  real t2
!
!  You might be out of stack space.
!
  if ( stack_num >= stack_max ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_QG_RPL_PUSH - Fatal error!'
    write ( *, '(a)' ) '  No more room on the stack.'
    stop
  end if

  stack_num = stack_num + 1

  stack1(1,stack_num) = i1sub
  stack1(2,stack_num) = i2sub
  stack1(3,stack_num) = j1sub
  stack1(4,stack_num) = j1type
  stack1(5,stack_num) = j2sub
  stack1(6,stack_num) = j2type

  stack2(1,stack_num) = t1
  stack2(2,stack_num) = t2

  return
end
subroutine r4vec2_sum_imax ( n, a, b, imax )

!*****************************************************************************80
!
!! R4VEC2_SUM_IMAX returns the index of the maximum sum of two R4VEC's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input, real A(N), B(N), two arrays whose sum is to
!    be examined.
!
!    Output, integer IMAX, the index of the largest entry in A+B.
!
  implicit none

  integer n

  real a(n)
  real b(n)
  integer i
  integer imax
  real sum_max

  if ( n <= 0 ) then

    imax = 0

  else

    imax = 1
    sum_max = a(1) + b(1)

    do i = 2, n
      if ( a(i) + b(i) > sum_max ) then
        sum_max = a(i) + b(i)
        imax = i
      end if
    end do

  end if

  return
end
function s_eqi ( strng1, strng2 )

!*****************************************************************************80
!
!! S_EQI is a case insensitive comparison of two strings for equality.
!
!  Example:
!
!    S_EQI ( 'Anjana', 'ANJANA' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRNG1, STRNG2, the strings to compare.
!
!    Output, logical S_EQI, the result of the comparison.
!
  implicit none

  character c1
  character c2
  integer i
  integer len1
  integer len2
  integer lenc
  logical s_eqi
  character ( len = * ) strng1
  character ( len = * ) strng2

  len1 = len ( strng1 )
  len2 = len ( strng2 )
  lenc = min ( len1, len2 )

  s_eqi = .false.

  do i = 1, lenc

    c1 = strng1(i:i)
    c2 = strng2(i:i)
    call ch_cap ( c1 )
    call ch_cap ( c2 )

    if ( c1 /= c2 ) then
      return
    end if

  end do

  do i = lenc + 1, len1
    if ( strng1(i:i) /= ' ' ) then
      return
    end if
  end do

  do i = lenc + 1, len2
    if ( strng2(i:i) /= ' ' ) then
      return
    end if
  end do

  s_eqi = .true.

  return
end
subroutine s_to_chvec ( s, n, cvec )

!*****************************************************************************80
!
!! S_TO_CHVEC converts a string to a character vector.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string of characters.
!
!    Input/output, integer N.
!    if N is -1, extract characters from 1 to len(S);
!    if N is 0, extract characters up to the last nonblank;
!    if N is positive, extract characters from 1 to N.
!
!    On output, N is the number of characters successfully extracted.
!
!    Output, character CVEC(N), the characters extracted from S.
!
  implicit none

  character cvec(*)
  integer i
  integer n
  character ( len = * ) s

  if ( n <= - 1 ) then
    n = len ( s )
  else if ( n == 0 ) then
    n = len_trim ( s )
  else
    n = min ( n, len ( s ) )
  end if

  do i = 1, n
    cvec(i) = s(i:i)
  end do

  return
end
subroutine s_to_i4 ( s, ival, ierror, last )

!*****************************************************************************80
!
!! S_TO_I4 reads an integer value from a string.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string to be examined.
!
!    Output, integer IVAL, the integer value read from the string.
!    If STRING is blank, then IVAL will be returned 0.
!
!    Output, integer IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer LAST, the last character of S used.
!
  implicit none

  character c
  integer i
  integer ierror
  integer isgn
  integer istate
  integer ival
  integer last
  character ( len = * ) s

  ierror = 0
  istate = 0
  isgn = 1
  ival = 0

  do i = 1, len_trim ( s )

    c = s(i:i)
!
!  Haven't read anything.
!
    if ( istate == 0 ) then

      if ( c == ' ' ) then

      else if ( c == '-' ) then
        istate = 1
        isgn = -1
      else if ( c == '+' ) then
        istate = 1
        isgn = + 1
      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read the sign, expecting digits.
!
    else if ( istate == 1 ) then

      if ( c == ' ' ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read at least one digit, expecting more.
!
    else if ( istate == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        ival = 10 * ival + ichar ( c ) - ichar ( '0' )
      else
        ival = isgn * ival
        last = i - 1
        return
      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( istate == 2 ) then
    ival = isgn * ival
    last = len_trim ( s )
  else
    ierror = 1
    last = 0
  end if

  return
end
subroutine sort_heap_external ( n, indx, i, j, isgn )

!*****************************************************************************80
!
!! SORT_HEAP_EXTERNAL externally sorts a list of items into linear order.
!
!  Discussion:
!
!    The actual list of data is not passed to the routine.  Hence this
!    routine may be used to sort integers, reals, numbers, names,
!    dates, shoe sizes, and so on.  After each call, the routine asks
!    the user to compare or interchange two items, until a special
!    return value signals that the sorting is completed.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 May 1999
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer N, the number of items to be sorted.
!
!    Input/output, integer INDX, the main communication signal.
!
!    The user must set INDX to 0 before the first call.
!    Thereafter, the user should not change the value of INDX until
!    the sorting is done.
!
!    On return, if INDX is
!
!      greater than 0,
!      * interchange items I and J;
!      * call again.
!
!      less than 0,
!      * compare items I and J;
!      * set ISGN = -1 if I < J, ISGN = +1 if I > J;
!      * call again.
!
!      equal to 0, the sorting is done.
!
!    Output, integer I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ISGN, results of comparison of elements I 
!    and J.
!    (Used only when the previous call returned INDX less than 0).
!    ISGN <= 0 means I is less than or equal to J;
!    ISGN => 0 means I is greater than or equal to J.
!
  implicit none

  integer i
  integer indx
  integer isgn
  integer j
  integer, save :: k = 0
  integer, save :: k1 = 0
  integer n
  integer, save :: n1 = 0
!
!  INDX = 0: This is the first call.
!
  if ( indx == 0 ) then

    n1 = n
    k = n / 2
    k1 = k
!
!  INDX < 0: The user is returning the results of a comparison.
!
  else if ( indx < 0 ) then

    if ( indx == -2 ) then

      if ( isgn < 0 ) then
        i = i + 1
      end if

      j = k1
      k1 = i
      indx = - 1
      return

    end if

    if ( isgn > 0 ) then
      indx = 2
      return
    end if

    if ( k <= 1 ) then

      if ( n1 == 1 ) then
        indx = 0
      else
        i = n1
        n1 = n1 - 1
        j = 1
        indx = 1
      end if

      return

    end if

    k = k - 1
    k1 = k
!
!  INDX > 0, the user was asked to make an interchange.
!
  else if ( indx == 1 ) then

    k1 = k

  end if

  do

    i = 2 * k1

    if ( i == n1 ) then
      j = k1
      k1 = i
      indx = - 1
      return
    else if ( i <= n1 ) then
      j = i + 1
      indx = - 2
      return
    end if

    if ( k <= 1 ) then
      exit
    end if

    k = k - 1
    k1 = k

  end do

  if ( n1 == 1 ) then
    indx = 0
  else
    i = n1
    n1 = n1 - 1
    j = 1
    indx = 1
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
subroutine word_last_read ( s, word )

!*****************************************************************************80
!
!! WORD_LAST_READ returns the last word from a string.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string containing words separated
!    by spaces.
!
!    Output, character ( len = * ) WORD, the last word.
!
  implicit none

  integer first
  integer last
  character ( len = * ) s
  character ( len = * ) word

  last = len_trim ( s )

  if ( last <= 0 ) then
    word = ' '
    return
  end if

  first = last

  do

    if ( first <= 1 ) then
      exit
    end if

    if ( s(first-1:first-1) == ' ' ) then
      exit
    end if

    first = first - 1

  end do

  word = s(first:last)

  return
end
subroutine word_next_read ( s, word, done )

!*****************************************************************************80
!
!! WORD_NEXT_READ "reads" words from a string, one at a time.
!
!  Discussion:
!
!    The following characters are considered to be a single word,
!    whether surrounded by spaces or not:
!
!      " ( ) { } [ ]
!
!    Also, if there is a trailing comma on the word, it is stripped off.
!    This is to facilitate the reading of lists.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string, presumably containing words
!    separated by spaces.
!
!    Output, character ( len = * ) WORD.
!
!    If DONE is FALSE, then WORD contains the "next" word read.
!    If DONE is TRUE, then WORD is blank, because there was no more to read.
!
!    Input/output, logical DONE.
!
!    On input with a fresh string, set DONE to TRUE.
!
!    On output, the routine sets DONE:
!      FALSE if another word was read,
!      TRUE if no more words could be read.
!
  implicit none

  logical done
  integer ilo
  integer, save :: lenc = 0
  integer, save :: next = 1
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )
  character ( len = * ) word
!
!  We "remember" LENC and NEXT from the previous call.
!
!  An input value of DONE = TRUE signals a new line of text to examine.
!
  if ( done ) then

    next = 1
    done = .false.
    lenc = len_trim ( s )

    if ( lenc <= 0 ) then
      done = .true.
      word = ' '
      return
    end if

  end if
!
!  Beginning at index NEXT, search the string for the next nonblank,
!  which signals the beginning of a word.
!
  ilo = next
!
!  ...S(NEXT:) is blank.  Return with WORD = ' ' and DONE = TRUE.
!
  do

    if ( ilo > lenc ) then
      word = ' '
      done = .true.
      next = lenc + 1
      return
    end if
!
!  If the current character is blank, skip to the next one.
!
    if ( s(ilo:ilo) /= ' ' .and. s(ilo:ilo) /= TAB ) then
      exit
    end if

    ilo = ilo + 1

  end do
!
!  ILO is the index of the next nonblank character in the string.
!
!  If this initial nonblank is a special character,
!  then that's the whole word as far as we're concerned,
!  so return immediately.
!
  if ( s(ilo:ilo) == '"' .or. &
       s(ilo:ilo) == '(' .or. &
       s(ilo:ilo) == ')' .or. &
       s(ilo:ilo) == '{' .or. &
       s(ilo:ilo) == '}' .or. &
       s(ilo:ilo) == '[' .or. &
       s(ilo:ilo) == ']' ) then

    word = s(ilo:ilo)
    next = ilo + 1
    return

  end if
!
!  Now search for the last contiguous character that is not a
!  blank, TAB, or special character.
!
  next = ilo + 1

  do while ( next <= lenc )

    if ( s(next:next) == ' ' ) then
      exit
    else if ( s(next:next) == TAB ) then
      exit
    else if ( s(next:next) == '"' ) then
      exit
    else if ( s(next:next) == '(' ) then
      exit
    else if ( s(next:next) == ')' ) then
      exit
    else if ( s(next:next) == '{' ) then
      exit
    else if ( s(next:next) == '}' ) then
      exit
    else if ( s(next:next) == '[' ) then
      exit
    else if ( s(next:next) == ']' ) then
      exit
    end if

    next = next + 1

  end do
!
!  Ignore a trailing comma.
!
  if ( s(next-1:next-1) == ',' ) then
    word = s(ilo:next-2)
  else
    word = s(ilo:next-1)
  end if

  return
end
