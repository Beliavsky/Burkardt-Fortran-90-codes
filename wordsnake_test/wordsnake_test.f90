program main

!*****************************************************************************80
!
!! wordsnake() computes a good wordsnake from a file of words.
!
!  Usage:
!
!    wordsnake wordsnake.inp
!
!    where "wordsnake.inp" is a list of words, one word per line,
!    to be made into a wordsnake.
!
!  Best so far:
!
!    0    0    0  sea  invent
!    1    1    1  inven T errible
!    3    9   10  terri BLE mish
!    4   16   26  ble MISH apen
!    3    9   35  misha PEN ultimate
!    2    4   39  penultima TE nse
!    2    4   43  ten SE em
!    2    4   47  se EM erge
!    5   25   72  e MERGE r
!    3    9   81  mer GER iatric
!    4   16   97  geria TRIC ky
!    1    1   98  trick Y es
!    2    4  102  y ES sential
!    2    4  106  essenti AL ly
!    1    1  107  all Y et
!    2    4  111  y ET ernal
!    2    4  115  etern AL as
!    3    9  124  a LAS ting
!    5   25  149  la STING er
!    3    9  158  stin GER und
!    3    9  167  ger UND erdevelop
!    4   16  183  underdev ELOP ed
!    3    9  192  elo PED iatric
!    4   16  208  pedia TRIC e
!    3    9  217  tr ICE r
!    3    9  226  i CER tain
!    2    4  230  certa IN credible
!    3    9  239  incredi BLE nd
!    4   16  255  b LEND ing
!    3    9  264  lend ING rate
!    4   16  280  ing RATE s
!    3    9  289  ra TES silate
!    4   16  305  tessi LATE r
!    3    9  314  la TER restrial
!    5   25  339  terres TRIAL s
!    1    1  340  trial S udden
!    3    9  349  sud DEN ude
!    2    4  353  denu DE nse
!    2    4  357  den SE a
!
!    Total number of characters =  254
!    Reduced number of characters =  145
!    Score is  357
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
  implicit none

  integer, parameter :: n = 39

  integer i
  integer iarg
  integer iargc
  character ( len = 255 ) filename  
  integer numarg
  integer perm(n)
  integer record_num
  integer score
  character ( len = 255 ), allocatable, dimension ( : ) :: word
!
!  Say hello.
!
  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'wordsnake_test()'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test wordsnake(), which creates a high scoring'
  write ( *, '(a)' ) '  wordsnake from a given set of words.'
!
!  Get the name of the input file.
!
  numarg = iargc ( )
!
!  Get the input file name.
!
  if ( 1 <= numarg ) then

    iarg = 1
    call getarg ( iarg, filename )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'What is the name of the input file?'
    read ( *, '(a)' ) filename
    if ( filename == ' ' ) then
      stop
    end if

  end if

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Reading word list file "' // trim ( filename ) // '"'
!
!  Count the number of records in the input file.
!
  call file_record_count ( filename, record_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of records in the input file is ', &
    record_num
!
!  Allocate the word array.
!
  allocate ( word(1:record_num) )
!
!  Read the word array.
!
  call file_record_read ( filename, record_num, word )
!
!  Print the word array.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The word list:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i2,2x,a)' ) i, trim ( word(i) )
  end do
!
!  Try to make the word snake.
!
  do i = 1, 3

    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Trial number ', i

    call wordsnake ( n, word, perm )

    call wordsnake_print ( n, word, perm )

    call wordsnake_score ( n, word, perm, score )

    write ( *, '(a,i6)' ) '  The wordsnake score is ', score

  end do
!
!  Free memory.
!
  deallocate ( word )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'wordsnake_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
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
