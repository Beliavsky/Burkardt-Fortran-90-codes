subroutine file_record_count ( file_in_name, record_num )

!*****************************************************************************80
!
!! file_record_count() counts the number of records in a file.
!
!  Discussion:
!
!    It does not count lines that are blank, or that begin with a
!    comment symbol '#'.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_IN_NAME, the name of the input file.
!
!    Output, integer RECORD_NUM, the number of records found.
!
  implicit none

  integer bad_num
  integer comment_num
  character ( len = * ) file_in_name
  integer file_in_unit
  integer ierror
  integer ios
  character ( len = 255 ) line
  integer line_num
  integer record_num

  call get_unit ( file_in_unit )

  open ( unit = file_in_unit, file = file_in_name, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_RECORD_COUNT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file: ' // &
      trim ( file_in_name )
    stop
  end if

  comment_num = 0
  record_num = 0
  line_num = 0
  bad_num = 0

  do

    read ( file_in_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      ierror = line_num
      exit
    end if

    line_num = line_num + 1

    if ( line(1:1) == '#' ) then
      comment_num = comment_num + 1
      cycle
    end if

    if ( len_trim ( line ) == 0 ) then
      comment_num = comment_num + 1
      cycle
    end if

    record_num = record_num + 1

  end do

  close ( unit = file_in_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FILE_RECORD_COUNT:'
  write ( *, '(a,i6)' ) '  Number of lines:           ', line_num
  write ( *, '(a,i6)' ) '  Number of data records:    ', record_num
  write ( *, '(a,i6)' ) '  Number of comment records: ', comment_num

  return
end
subroutine file_record_read ( file_in_name, record_num, word )

!*****************************************************************************80
!
!! file_record_read() reads the records in a file.
!
!  Discussion:
!
!    It does not count lines that are blank, or that begin with a
!    comment symbol '#'.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_IN_NAME, the name of the input file.
!
!    Input, integer RECORD_NUM, the number of records.
!
  implicit none

  integer record_num

  integer comment_num
  character ( len = * ) file_in_name
  integer file_in_unit
  integer ierror
  integer ios
  character ( len = 255 ) line
  integer line_num
  integer record_num2
  character ( len = 255 ) word(record_num)

  call get_unit ( file_in_unit )

  open ( unit = file_in_unit, file = file_in_name, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_RECORD_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file: ' // &
      trim ( file_in_name )
    stop
  end if

  record_num2 = 0

  do 

    read ( file_in_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      ierror = line_num
      exit
    end if

    line_num = line_num + 1

    if ( line(1:1) == '#' ) then
      comment_num = comment_num + 1
      cycle
    end if

    if ( len_trim ( line ) == 0 ) then
      comment_num = comment_num + 1
      cycle
    end if

    record_num2 = record_num2 + 1

    if ( record_num2 <= record_num ) then
      word(record_num2) = line
    end if

  end do

  close ( unit = file_in_unit )

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! get_unit() returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IUNIT.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5 and 6).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
  implicit none

  integer i
  integer ios
  integer iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
function i4_modp ( i, j )

!*****************************************************************************80
!
!! i4_modp() returns the nonnegative remainder of I4 division.
!
!  Discussion:
!
!    If
!      NREM = I4_MODP ( I, J )
!      NMULT = ( I - NREM ) / J
!    then
!      I = J * NMULT + NREM
!    where NREM is always nonnegative.
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
!
!  Example:
!
!        I     J     MOD  I4_MODP    Factorization
!
!      107    50       7       7    107 =  2 *  50 + 7
!      107   -50       7       7    107 = -2 * -50 + 7
!     -107    50      -7      43   -107 = -3 *  50 + 43
!     -107   -50      -7      43   -107 =  3 * -50 + 43
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, the number to be divided.
!
!    Input, integer J, the number that divides I.
!
!    Output, integer I4_MODP, the nonnegative remainder when I is
!    divided by J.
!
  implicit none

  integer i
  integer i4_modp
  integer j

  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MODP - Fatal error!'
    write ( *, '(a,i6)' ) '  I4_MODP ( I, J ) called with J = ', j
    stop
  end if

  i4_modp = mod ( i, j )

  if ( i4_modp < 0 ) then
    i4_modp = i4_modp + abs ( j )
  end if

  return
end
subroutine i4_random ( ilo, ihi, i )

!*****************************************************************************80
!
!! i4_random() returns a random integer in a given range.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ILO, IHI, the minimum and maximum acceptable
!    values.
!
!    Output, integer I, the randomly chosen integer.
!
  implicit none

  integer i
  integer ihi
  integer ilo
  real r
  real rhi
  real rlo
  logical, save :: seed = .false.

  if ( .not. seed ) then
    call random_seed ( )
    seed = .true.
  end if
!
!  Pick a random number in (0,1).
!
  call random_number ( harvest = r )
!
!  Set a real interval [RLO,RHI] which contains the integers [ILO,IHI],
!  each with a "neighborhood" of width 1.
!
  rlo = real ( ilo ) - 0.5E+00
  rhi = real ( ihi ) + 0.5E+00
!
!  Set I to the integer that is nearest the scaled value of R.
!
  i = nint ( ( 1.0E+00 - r ) * rlo + r * rhi )
!
!  In case of oddball events at the boundary, enforce the limits.
!
  i = max ( i, ilo )
  i = min ( i, ihi )

  return
end
subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! i4_swap() switches two I4's.
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
subroutine i4_swap3 ( i, j, k )

!*****************************************************************************80
!
!! i4_swap3() swaps three integer values.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer I, J, K.  On output, the values
!    have been interchanged.
!
  implicit none

  integer i
  integer j
  integer k
  integer l

  l = i
  i = j
  j = k
  k = l

  return
end
subroutine i4_unswap3 ( i, j, k )

!*****************************************************************************80
!
!! i4_unswap3() unswaps three integer values.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer I, J, K.  On output, the values
!    have been interchanged.
!
  implicit none

  integer i
  integer j
  integer k
  integer l

  l = k
  k = j
  j = i
  i = l

  return
end
function i4_wrap ( ival, ilo, ihi )

!*****************************************************************************80
!
!! i4_wrap() forces an I4 to lie between given limits by wrapping.
!
!  Example:
!
!    ILO = 4, IHI = 8
!
!    I  I4_WRAP
!
!    -2     8
!    -1     4
!     0     5
!     1     6
!     2     7
!     3     8
!     4     4
!     5     5
!     6     6
!     7     7
!     8     8
!     9     4
!    10     5
!    11     6
!    12     7
!    13     8
!    14     4
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer IVAL, an integer value.
!
!    Input, integer ILO, IHI, the desired bounds.
!
!    Output, integer I4_WRAP, a "wrapped" version of IVAL.
!
  implicit none

  integer i4_modp
  integer i4_wrap
  integer ihi
  integer ilo
  integer ival
  integer wide

  wide = ihi + 1 - ilo

  if ( wide == 0 ) then
    i4_wrap = ilo
  else
    i4_wrap = ilo + i4_modp ( ival-ilo, wide )
  end if

  return
end
subroutine i4vec_identity ( n, a )

!*****************************************************************************80
!
!! i4vec_identity() sets an I4VEC to the identity vector A(I)=I.
!
!  Modified:
!
!    09 November 2000
!
!  Parameters:
!
!    Input, integer N, the number of elements of A.
!
!    Output, integer A(N), the array to be initialized.
!
  implicit none

  integer n

  integer a(n)
  integer i

  do i = 1, n
    a(i) = i
  end do

  return
end
function lower ( s )

!*****************************************************************************80
!
!! lower() returns a lowercase version of a string.
!
!  Discussion:
!
!    LOWER is a string function of undeclared length.  The length
!    of the argument returned is determined by the declaration of
!    LOWER in the calling routine.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string.
!
!    Output, character ( len = * ) LOWER, a lowercase copy of the string.
!
  implicit none

  integer i
  integer j
  character ( len = * ) lower
  integer n
  character ( len = * ) s

  lower = s

  n = len_trim ( lower )

  do i = 1, n

    j = ichar ( lower(i:i) )

    if ( 65 <= j .and. j <= 90 ) then
      lower(i:i) = char ( j + 32 )
    end if

  end do

  return
end
subroutine overlap_table ( n, word, table )

!*****************************************************************************80
!
!! overlap_table() computes a table of the overlap between pairs of words.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of words.
!
!    Input, character ( len = * ) WORD(N), a list of words.
!
!    Output, integer TABLE(N,N), a table containing, in TABLE(I,J),
!    the number of characters by which the end of WORD(I) and the
!    beginning of word J overlap.
!
  implicit none

  integer n

  integer i
  integer j
  integer table(n,n)
  character ( len = * ) word(n)
!
!  Construct the overlap table.
!
  do i = 1, n
    do j = 1, n

      call s_overlap ( word(i), word(j), table(i,j) )

    end do
  end do

  return
end
subroutine perm_random ( n, p, setup )

!*****************************************************************************80
!
!! perm_random() selects a random permutation of N objects.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of objects to be permuted.
!
!    Input/output, integer P(N), a permutation, in standard 
!    index form.
!
!    If SETUP is .TRUE., then the input value of P contains
!    the "current" labels of the objects.
!
!    Otherwise, P(I) is initialized to I.
!
!    On output, P(I) is the randomly permuted label of the I-th object.
!
!    Input, logical SETUP.
!
!    If SETUP is .TRUE. then the routine assumes the objects
!    are labeled 1, 2, ... N.
!
!    If SETUP is .FALSE., then the input values of P are used
!    as labels; that is, the I-th object is labeled P(I).
!
  implicit none

  integer n

  integer i
  integer j
  integer p(n)
  logical setup

  if ( setup ) then
    call i4vec_identity ( n, p )
  end if

  do i = 1, n
    call i4_random ( i, n, j )
    call i4_swap ( p(i), p(j) )
  end do

  return
end
subroutine s_overlap ( s1, s2, overlap )

!*****************************************************************************80
!
!! s_overlap() determines the overlap between two strings.
!
!  Discussion:
!
!    To determine the overlap, write the first word followed immediately
!    by the second word.  Find the longest substring S which is both
!    a suffix of S1 and a prefix of S2.  The length of this substring
!    is the overlap.
!
!  Example:
!
!    S1              S2        OVERLAP
!
!    'timber'        'beret'   3
!    'timber'        'timber'  6
!    'beret'         'timber'  1
!    'beret'         'berets'  5
!    'beret'         'berth'   0
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S1, S2, the strings to be checked.
!
!    Output, integer OVERLAP, the length of the overlap.
!
  implicit none

  integer i
  integer len1
  integer len2
  integer len3
  integer overlap
  character ( len = * ) s1
  character ( len = * ) s2

  overlap = 0

  len1 = len_trim ( s1 )
  len2 = len_trim ( s2 )
  len3 = min ( len1, len2 )

  do i = 1, len3
    if ( s1(len1+1-i:len1) == s2(1:i) ) then
      overlap = i
    end if
  end do

  return
end
function upper ( s )

!*****************************************************************************80
!
!! upper() returns an uppercase version of a string.
!
!  Discussion:
!
!    UPPER is a string function of undeclared length.  The length
!    of the argument returned is determined by the declaration of
!    UPPER in the calling routine.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string.
!
!    Output, character ( len = * ) UPPER, an uppercase copy of the string.
!
  implicit none

  integer i
  integer j
  integer n
  character ( len = * ) s
  character ( len = * ) upper

  upper = s

  n = len_trim ( upper )

  do i = 1, n

    j = ichar ( upper(i:i) )

    if ( 97 <= j .and. j <= 122 ) then
      upper(i:i) = char ( j - 32 )
    end if

  end do

  return
end
subroutine wordsnake ( n, word, perm )

!*****************************************************************************80
!
!! wordsnake() seeks a high scoring permutation of a set of words.
!
!  Discussion:
!
!    Note that MISHAPEN and TESSILATE, as given in the article,
!    are incorrectly spelled.
!
!    A wordsnake is formed from a list of words.  The words are rearranged
!    so that, where possible, the end of one word matches the beginning of
!    the next.  The highest scoring wordsnake is sought.  This is a hard
!    problem, and is similar, really, to the traveling salesman problem.
!
!  Scoring:
!
!    The score for a wordsnake is the sum of the squares of the number of
!    letters of overlap between pairs of successive words.  (The article
!    does not mention whether to count overlap for the possible wrap-around
!    from the last word to the first again.)
!
!    For instance, consider the wordsnake made up of:
!
!      "writes", "testate", "tater", "remote", "stew"
!
!    and consider the resulting wordsnake:
!
!      writestateremotestew
!
!    which is scored as:
!
!      Word pair     Overlap  Score
!      ------------  -------  -----
!      wri TES tate     3        +9
!      tes TATE r       4       +16
!      tate R emote     1        +1
!      remote stew      0        +0
!      ste W rites      1        +1
!
!      TOTAL SCORE               27
!      TOTAL CHARACTERS          20
!
!    Using this scoring system, and the set of 39 words given below, the 
!    article claims a score of 357 was achieved, with a wordsnake length of 
!    145 characters (that is, printing out the shared letters only once).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    28 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dennis Shasha,
!    Wordsnakes,
!    Dr Dobb's Journal,
!    July, 2000, pages 143-144.
!
!  Parameters:
!
!    Input, integer N, the number of words.
!
!    Input, character ( len = * ) WORD(N), the words to be used.
!
!    Output, integer PERM(N), a permutation of the words to be used
!    in the wordsnake.
!
  implicit none

  integer n

  logical, parameter :: DEBUG = .false.
  integer i
  integer perm(n)
  integer perm_best(n)
  integer score
  integer score_best
  integer score_best_old
  integer table(n,n)
  character ( len = * ), dimension ( n ) :: word
!
!  Compute the overlap table.
!
  call overlap_table ( n, word, table )

  if ( DEBUG ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Overlap table:'
    write ( *, '(a)' ) ' '
    write ( *, '(6x,39i2)' ) ( i, i = 1, n )
    do i = 1, n
      write ( *, '(4x,40i2)' ) i, table(i,1:n)
   end do

  end if
!
!  Wordsnakes will be described by a permutation of our word list.
!  Our initial wordsnake is the identity permutation.
!
  call perm_random ( n, perm, .true. )

  if ( DEBUG ) then
    call wordsnake_print ( n, word, perm )
  end if

  call wordsnake_score ( n, word, perm, score )

  if ( DEBUG ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) 'Initial score is ', score
  end if

  perm_best(1:n) = perm(1:n)
  score_best = score
!
!  Do one naive "greedy" pass, where we follow word(I) by its optimal follower.
!
  call wordsnake_search_greedy ( n, perm, table )

  if ( DEBUG ) then
    call wordsnake_print ( n, word, perm )
  end if

  call wordsnake_score ( n, word, perm, score )

  if ( DEBUG ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) 'Greedy score is ', score
  end if

  if ( score_best < score ) then
    perm_best(1:n) = perm(1:n)
    score_best = score
  end if

  do

    score_best_old = score_best
!
!  Consider all possible inserts.
!
    call wordsnake_search_insert ( n, word, perm_best, score_best )
 
    if ( DEBUG ) then
      call wordsnake_print ( n, word, perm_best )
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) 'INSERT score is ', score_best
    end if
!
!  Consider all possible pair swaps.
!
    call wordsnake_search_swap2 ( n, word, perm_best, score_best )
 
    if ( DEBUG ) then
      call wordsnake_print ( n, word, perm_best )
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) 'SWAP2 score is ', score_best
    end if
!
!  Consider all possible trio swaps.
!
    call wordsnake_search_swap3 ( n, word, perm_best, score_best )
 
    if ( DEBUG ) then
      call wordsnake_print ( n, word, perm_best )
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) 'SWAP3 score is ', score_best
    end if
!
!  Now consider all switches of the form -A-B-C-A- to -A-C-B-A-
!
    call wordsnake_search_transpose ( n, word, perm_best, score_best )

    if ( DEBUG ) then
      call wordsnake_print ( n, word, perm_best )
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) 'TRANSPOSE score is ', score_best
    end if

    if ( score_best == score_best_old ) then
      exit
    end if

  end do

  return
end
subroutine wordsnake_print ( n, word, perm )

!*****************************************************************************80
!
!! wordsnake_print() prints a wordsnake.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of words.
!
!    Input, character ( len = * ) WORD(N), a list of words.
!
!    Input, integer PERM(N), a permutation, the ordering of the
!    words for the wordsnake.
!
  implicit none

  integer n

  integer i
  integer inc
  character ( len = 25 ) lower
  integer na
  integer n_char
  integer n_overlap
  integer overlap
  integer perm(n)
  character ( len = 25 ) s1
  character ( len = 25 ) s2
  character ( len = 25 ) sa
  character ( len = 25 ) sb
  character ( len = 25 ) sc
  integer score
  character, parameter :: space = ' '
  character ( len = 25 ) upper
  character ( len = * ) word(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'wordsnake():'
  write ( *, '(a,i8)' ) '  Number of words = ', n
  write ( *, '(a)' ) ' '

  score = 0
  n_overlap = 0

  do i = 1, n

    s1 = word(perm(i))

    if ( i < n ) then
      s2 = word(perm(i+1))
    else
      s2 = word(perm(1))
    end if

    call s_overlap ( s1, s2, overlap )

    n_overlap = n_overlap + overlap

    inc = overlap**2
    score = score + inc 

    na = len_trim ( s1 )
    sa = lower ( s1(1:na-overlap) )
    sb = upper ( s2(1:overlap) )
    sc = lower ( s2(overlap+1:) )

    write ( *, '(3i8,2x,5a)' ) overlap, inc, score, &
      trim ( sa ), space, trim ( sb ), space, trim ( sc )

  end do

  n_char = 0
  do i = 1, n
    n_char = n_char + len_trim ( word(i) )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Total number of characters = ', n_char
  write ( *, '(a,i8)' ) '  Reduced number of characters = ', n_char - n_overlap

  return
end
subroutine wordsnake_score ( n, word, perm, score )

!*****************************************************************************80
!
!! wordsnake_score() computes the score for a given wordsnake.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of words.
!
!    Input, character ( len = * ) WORD(N), a list of words.
!
!    Input, integer PERM(N), a permutation, the ordering of the
!    words for the wordsnake.
!
!    Output, integer SCORE, the overlap score for this wordsnake.
!
  implicit none

  integer n

  integer i
  integer ip1
  integer j
  integer jp1
  integer overlap
  integer perm(n)
  integer score
  character ( len = * ) word(n)

  score = 0

  do i = 1, n

    j = perm(i)

    if ( i < n ) then
      ip1 = i + 1
    else
      ip1 = 1
    end if

    jp1 = perm(ip1)

    call s_overlap ( word(j), word(jp1), overlap )

    score = score + overlap**2

  end do

  return
end
subroutine wordsnake_search_greedy ( n, perm, table )

!*****************************************************************************80
!
!! wordsnake_search_greedy() constructs a wordsnake using a greedy algorithm.
!
!  Discussion:
!
!    The routine takes the first word, and puts it into the wordsnake.
!    Then it takes the next word that is available, and would produce the
!    best matching score with the previous word, and so on.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    28 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of words.
!
!    Input/output, integer PERM(N), a permutation, the ordering 
!    of the words for the wordsnake.
!
!    Input, TABLE(N,N), a table of the overlap scores for each pair of words.
!
  implicit none

  integer :: n

  integer i
  integer ip
  integer j
  integer jp
  integer k
  integer kp
  integer perm(n)
  integer table(n,n)

  do i = 1, n - 1

    ip = perm(i)

    j = i + 1
    jp = perm(i+1)

    do k = i + 2, n
      kp = perm(k)

      if ( table(ip,kp) > table(ip,jp) ) then
        j = k
        jp = kp
      end if

    end do

    call i4_swap ( perm(i+1), perm(j) )

  end do

  return
end
subroutine wordsnake_search_insert ( n, word, perm_best, score_best )

!*****************************************************************************80
!
!! wordsnake_search_insert() tries to improve the score by inserting a word.
!
!  Discussion:
!
!    The operation considered here modifies the string by moving a single
!    word from one position to another.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of words.
!
!    Input, character ( len = * ) WORD(N), a list of words.
!
!    Input/output, integer PERM_BEST(N), a permutation, the
!    ordering of the words for the wordsnake.
!
!    Input/output, integer SCORE_BEST, the best score so far.
!
  implicit none

  integer n

  logical, parameter :: DEBUG = .false.
  integer i
  integer ii
  integer improve
  integer n_insert
  integer pi_new
  integer pi_old
  integer perm(n)
  integer perm_best(n)
  integer score
  integer score_best
  character ( len = * ) word(n)

  n_insert = 0
  perm(1:n) = perm_best(1:n)
  improve = 0
!
!  For word I, which is now in position PI_OLD = PERM(I), consider inserting
!  it in position PI_NEW.
!
  do i = 1, n

    do pi_new = 1, n

      pi_old = perm(i)

      if ( pi_new /= pi_old ) then

        if ( pi_new < pi_old ) then
          do ii = 1, n
            if ( pi_new <= perm(ii) .and. perm(ii) < pi_old ) then
              perm(ii) = perm(ii) + 1
            else if ( perm(ii) == pi_old ) then
              perm(ii) = pi_new
            end if
          end do
        else if ( pi_new > pi_old ) then
          do ii = 1, n
            if ( pi_old < perm(ii) .and. perm(ii) <= pi_new ) then
              perm(ii) = perm(ii) - 1
            else if ( perm(ii) == pi_old ) then
              perm(ii) = pi_new
            end if
          end do
        end if

        call wordsnake_score ( n, word, perm, score )

        if ( score >= score_best ) then
          improve = improve + score - score_best
          score_best = score
          perm_best(1:n) = perm(1:n)
          n_insert = n_insert + 1
        else
          perm(1:n) = perm_best(1:n)
        end if

      end if
    end do
  end do

  if ( DEBUG ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Number of inserts = ', n_insert
    write ( *, '(a,i6)' ) '  Score improvement = ', improve
  end if

  return
end
subroutine wordsnake_search_swap2 ( n, word, perm_best, score_best )

!*****************************************************************************80
!
!! wordsnake_search_swap2() tries to improve the score by swapping 2 words.
!
!  Discussion:
!
!    The operation essentially checks the possibility that the score would
!    be improved by the following 2-Swap operation:
!
!      Temp  <- Word1
!      Word1 <- Word2
!      Word2 <- Temp.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of words.
!
!    Input, character ( len = * ) WORD(N), a list of words.
!
!    Input/output, integer PERM_BEST(N), a permutation, the 
!    ordering of the words for the wordsnake.
!
!    Input/output, integer SCORE_BEST, the best score so far.
!
  implicit none

  integer n

  logical, parameter :: DEBUG = .false.
  integer i
  integer improve
  integer j
  integer n_swap
  integer perm(n)
  integer perm_best(n)
  integer score
  integer score_best
  character ( len = * ) word(n)

  n_swap = 0
  perm(1:n) = perm_best(1:n)
  improve = 0
!
!  For word I, which is now in position PERM(I), consider swapping
!  it with the word in position PERM(J).
!
  do i = 1, n
    do j = 1, n
      if ( i /= j ) then
        call i4_swap ( perm(i), perm(j) )
        call wordsnake_score ( n, word, perm, score )
        if ( score >= score_best ) then
          improve = improve + score - score_best
          score_best = score
          perm_best(1:n) = perm(1:n)
          n_swap = n_swap + 1
        else
          call i4_swap ( perm(i), perm(j) )
        end if
      end if
    end do
  end do

  if ( DEBUG ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Number of swaps = ', n_swap
    write ( *, '(a,i6)' ) '  Score improvement = ', improve
  end if

  return
end
subroutine wordsnake_search_swap3 ( n, word, perm_best, score_best )

!*****************************************************************************80
!
!! wordsnake_search_swap3() tries to improve the score by swapping 3 words.
!
!  Discussion:
!
!    The operation essentially checks the possibility that the score would
!    be improved by the following 3-Swap operation:
!
!      Temp  <- Word1
!      Word1 <- Word2
!      Word2 <- Word3
!      Word3 <- Temp.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of words.
!
!    Input, character ( len = * ) WORD(N), a list of words.
!
!    Input/output, integer PERM_BEST(N), a permutation, the 
!    ordering of the words for the wordsnake.
!
!    Input/output, integer SCORE_BEST, the best score so far.
!
  implicit none

  integer n

  logical, parameter :: DEBUG = .false.
  integer i
  integer improve
  integer j
  integer k
  integer n_swap
  integer perm(n)
  integer perm_best(n)
  integer score
  integer score_best
  character ( len = * ) word(n)

  n_swap = 0
  perm(1:n) = perm_best(1:n)
  improve = 0
!
!  For word I, which is now in position PERM(I), consider putting it
!  where PERM(J) is.
!
  do i = 1, n
    do j = 1, n
      if ( i /= j ) then
        do k = 1, n
          if ( k /= i .and. k /= j ) then
            call i4_swap3 ( perm(i), perm(j), perm(k) )
            call wordsnake_score ( n, word, perm, score )
            if ( score >= score_best ) then
              improve = improve + score - score_best
              score_best = score
              perm_best(1:n) = perm(1:n)
              n_swap = n_swap + 1
            else
              call i4_unswap3 ( perm(i), perm(j), perm(k) )
            end if
          end if
        end do
      end if
    end do
  end do

  if ( DEBUG ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Number of 3-swaps = ', n_swap
    write ( *, '(a,i6)' ) '  Score improvement = ', improve
  end if

  return
end
subroutine wordsnake_search_transpose ( n, word, perm_best, score_best )

!*****************************************************************************80
!
!! wordsnake_search_transpose() tries to improve the score using transpositions.
!
!  Discussion:
!
!    In a transposition, the wordsnake is essentially cut in three places,
!    and two pieces are swapped.  We can think of the operation as 
!    replacing:
!
!      -A-B-C-A-
!
!    by
!
!      -A-C-B-A-
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    28 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of words.
!
!    Input, character ( len = * ) WORD(N), a list of words.
!
!    Input/output, integer PERM_BEST(N), a permutation, the 
!    ordering of the words for the wordsnake.
!
!    Input/output, integer SCORE_BEST, the best score so far.
!
  implicit none

  integer n

  logical, parameter :: DEBUG = .false.
  integer i1
  integer i2
  integer i3
  integer i4
  integer improve
  integer n_transpose
  integer perm(n)
  integer perm2(n)
  integer perm_best(n)
  integer score
  integer score_best
  character ( len = * ) word(n)

  n_transpose = 0
  perm(1:n) = perm_best(1:n)
  improve = 0
!
!  Consider I0:I1:I2:I3:I4:I5 and transposing to I0:I3:I4:I1:I2:I5.
!
  do i1 = 1, n-1
    do i2 = i1, n-1
      do i4 = i2+1, n

        i3 = i2+1

        perm2(1:n) = perm(1:n)
        perm2(i1:i4+i1-i3) = perm(i3:i4)
        perm2(i4+i1-i2:i4) = perm(i1:i2)

        call wordsnake_score ( n, word, perm2, score )

        if ( score_best <= score ) then
          improve = improve + score - score_best
          score_best = score
          perm_best(1:n) = perm2(1:n)
          n_transpose = n_transpose + 1
        end if

      end do
    end do
  end do

  if ( DEBUG ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Number of transpositions = ', n_transpose
    write ( *, '(a,i6)' ) '  Score improvement = ', improve
  end if

  return
end
