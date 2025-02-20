program main

!*****************************************************************************80
!
!! matman() is a program for interactive linear algebra demonstrations.
!
!  Discussion:
!
!    Improvements done:
!    * header ended with ****80.
!    * licensing statement inside of each function.
!    * removed lonely "!" ahead of ***80 line.
!    * !! line is followed by just 1 blank line.
!    * replaced "Note:" by "Discussion:".
!    * lowercase function name, followed by ().
!    * removed unused variables.
!    * replaced ** by ^ for exponentiation.
!    * Rephrased ">" inequalities using "<" instead.
!    * Each routine started with implicit none.
!    * Expanded one-line if statements to if / statement / endif.
!    * removed single empty "!" comments.
!    * converted to double precision.
!    Improvements needed:
!    * compare run with F77 version, which is failing.
!    * Rephrase Parameters: as Input: Output:
!
!  Chronicle:
!
!    22 September 2000
!
!    Pitched the key file and authorization.
!
!    25 May 2000
!
!    FORTRAN90 free format.
!
!    15 April 2000
!
!    FORTRAN90: new logical operator names, use LEN_TRIM, 
!    new form of PARAMETER statement.
!
!    An improved random orthogonal matrix routine was added.
!
!    19 February 2000
!
!    Added elementary column operations.
!
!    Experimenting with orthogonal row operations.  
!
!    Deleting baggage: I got rid of the determinant variables and routines, 
!    the transpose operation, and the warnings.
!
!    Version 1.61  25 January 2000
!
!    Making routine names more sensible.
!
!    I noticed that you can't issue the ERO "R4 => R4 * 10".  It wants
!    the multiplier in front.
!
!    Working on the decimal arithmetic.  
!
!    The first decimal arithmetic problem was that I could type "0.01" 
!    and the decimal reading routine would end up reporting "0.009999" 
!    because of internal conversions.  I think I fixed that.
!
!    Found an obscure typo in DEC_MUL that made the exponent wrong.
!
!    Version 1.60  26 September 1998
!
!    * I added eigenvector and one-step Q information to EVJACO.
!      I added ORTRAN to get nicer examples.
!
!    * Renamed "SCADIV" to "ROWDIV", and made division by fractions
!      print as multiplication.
!
!    * I'm getting program failures in DEC_PRINT, which I fixed.
!
!    * The program was simply stopping when encountering overflow
!      in rational calculations.  I fixed RAT_ADD and RAT_MUL.  
!      I replaced calls to RAT_DIV by calls to RAT_MUL.
!
!    Version 1.59   06 July 1998
!
!    * Made the linear algebra sample problems random.
!    * Forced A(I,J) = A(J,I) = 0 exactly in Jacobi method.
!    * Installed newer version of PERM_NEXT routine.
!    * Used new format for argument documentation.
!    * Inserted latest versions of CHRPAK/SUBPAK/SUBSET routines.
!    * Passed maximum integer to DEC and RAT routines, allowed user to set it.
!    * Used a better RANDOM routine, and updated MATKEY.DAT.
!
!    Version 1.58  10 April 1996
!
!    CGC requested the following changes:
!
!    1) The confirmation of the A command should be
!
!       Row 5 <= 3 Row 2 + Row 5, with "Row 5" coming last.  (FIXED)
!
!    2) In LP mode, the feasibility ratios are printed out in
!       both real and rational values, but they disagree.  (FIXED)
!
!    I made the following changes:
!
!    3) I also tried to add some more comments at the beginnings of
!       routines, to define the variables.
!
!    4) I made R_TO_S_RIGHT print out in G14.7 format, to try to get seven
!       digit accuracy where possible, and modified R_PRINT to print out
!       7 decimals as well.
!
!    5) I also altered SETDIG to allow the user to exceed the recommended
!       maximum number of digits, with a warning.
!
!    6) I modified the linear algebra optimization checker to print out
!       the column, as well as the row, where rule 1 is violated.
!
!    7) I changed RELPRN so that if the printed quantities are all integers,
!       the printout is more compact.
!
!    8) I modified AUTO_ERO to eliminate unnecessary row operations,
!       where an entry is already zero.
!
!    9) Added an extra check in ROWADD to skip out immediately if the
!       multiplier is zero.
!
!    10) Added documentation for the DECIMAL, RATIONAL, and real
!       commands, which will eventually replace the "F" command.
!
!    11) Updated my address and EMAIL address.
!
!
!    Version 1.57  15 December 1995
!
!    CGC requested the following changes:
!
!    1) In LP mode, when doing a two phase problem, if you convert
!       from one arithmetic form to another, the objective function
!       data in row NROW+1 is not converted.  So I modified FORM
!       to convert all the rows and columns, not just NROW by NCOL.
!       FIXED.
!
!    2) Move the L command to the short menu, and add to the short
!       menu a description of how to get the long menu.  DONE.
!
!    3) Replace the prompt which follows the "Enter command" prompt
!       with ("H" for short menu, "HELP" for full menu, ? for full help).
!       DONE.
!
!    4) Stumbled across a slight problem.  When in LP mode, using
!       fractional arithmetic, and you enter a problem with artificial
!       variables, the last column of the auxilliary objective function
!       was set to 0/0 instead of 0/1.  FIXED.
!
!    5) Instead of printing the ERO determinant after every operation,
!       I added an EDET command to print it out only on demand.
!       This is a request CGC made earlier, and which we had both
!       forgotten.
!
!    6) Cosmetic change: rational matrices were printed out with two
!       trailing blank lines, which I cut back to one.
!
!    7) Fixed an obscure error in R_READ, which only caused problems
!       on the ALPHA.  CHRCTR was reading numbers all the way to the
!       last blank, setting LCHAR=NCHAR=80, and then R_READ was asking
!       if character LCHAR+1 was a '/'.
!
!    8) I replaced all the "WRITE(*,*)" statements by the more robust
!       "OUTPUT=...", "CALL S_WRITE()" pair.  Some error messages were
!       only going to the screen, and not the permanent output file.
!
!    9) Discovered that when a sample problem is chosen, only NROW
!       by NCOL of the matrix area was set, leaving garbage possibly
!       in other areas.  I rectified this, zeroing out all the
!       rest of the matrix area, via a routine INIMAT.
!
!    Version 1.56  12 October 1995
!
!    I fixed the program, so that you can type "Row 2 <=> Row 3"
!    (spelling out the word "Row") if you want to.
!
!    I found a logic mistake in DEC_MUL which occurs if one of the
!    input quantities is the same as the output, and I fixed it.
!
!    I modified the program so that, if you are working with 4
!    digit decimals, any decimal input is automatically truncated
!    to 4 decimals as well.
!
!    Replaced CHLDEC by a routine which is exact.
!
!    Corrected DEC_TO_RAT, and many other decimal discrepancies.
!
!    Corrected DEC_ADD, so that decimal addition is exact.
!
!    Added the DECIMAL, RATIONAL and real commands, although I
!    did not mention them.
!
!    Corrected the phrase "row reduced echelon form" by replacing
!    it with "reduced row echelon form".
!
!    Added the BASIC command to allow the user to assign a row
!    to a basic variable without using the Change command.
!
!    An unneeded change to CHRINP disabled the "<" command, but
!    I fixed that.
!
!    I modified routine PASS so that its default key corresponds
!    to the value currently stored in MATKEY.DAT.
!
!    Corrected TRANSC so that VMS output files would have correct
!    carriagecontrol.
!
!    Version 1.55  06 October 1995
!
!    Experimenting with allowing longer command names, so that I
!    can have more reasonable names.  Changed COMMAND to 4 characters
!    in length.  Now I can request a determinant with DET, and
!    a transpose with TR.  Note that TR is potentially in conflict
!    with the "Type a Row" command, unless the user puts a space there.
!
!    Added the determinant command.
!
!    Dropped the T C, T R and T E commands.
!
!    Added a square matrix example for the determinant problem.
!
!    Print out the determinant of the ERO's.
!
!    The "<" has a logical flaw.  If you use the "X" command in
!    the input file, you will return to the user input, not the
!    file input, once you're done.  You need to save each input
!    unit number and file name in a stack in order to properly
!    recover.
!
!    If you add NCON, you're going to have to save it in RESTORE
!    and in read/write examples as well...
!
!    CGC requested the ability to add a row or column for linear
!    programming!  In this case, it would correspond to a slack
!    variable.
!
!    Version 1.54  05 September 1994
!
!    EVJACO allows the user to type an integer or a character.
!    But I_READ was allowing an error message to appear if the
!    user typed a "Q", because the IHUSH parameter was reset
!    to 0.  I took out the resetting.
!
!    Because of capability of entering several commands, separated
!    by a semicolon, comments weren't being ignored properly.
!    So now, once a comment is recognized, it's blanked out.
!
!    Version 1.53  25 July 1994
!
!    * Replaced constant "0" by variable "ITERM" in all calls to S_READ.
!    * Added "Y" command, so user can turn autoprinting off or on.
!    * Autoprint after "V" command.
!    * Added capability to enter several commands, separated by a semicolon.
!
!    Version 1.52  12 May 1994
!
!    12 May 1994:
!
!    Corrected SETDIG so that the value the user typed in does
!    not immediately overwrite the current value, until it has
!    been checked.
!
!    11 May 1994:
!
!    Struggling with a "final" problem, in which DEC_PRINT does
!    not work when trying to print out the linear programming
!    solution for the advanced sample problem.  Think it's
!    fixed now.
!
!    10 May 1994:
!
!    Program freezes when going into decimal fraction mode,
!    with simple linear programming problem.
!
!    It would be preferable to delete trailing zeroes from the
!    printouts of real numbers by RELPRN.
!
!    Tracking down a bug in DECREA, I think.
!
!    The linear algebra stuff seems to be working properly with
!    decimal arithmetic.  I still need to check linear
!    programming.
!
!    Added an option to use a default key.
!
!    Inserted new version of CHLDEC which does not try to convert
!    IVAL*10^JVAL into a real first, and so should be able to print
!    out the exact representation, with no trailing blanks or
!    roundoff problems.  Looks a lot better.
!
!    I corrected the sample problem routines, to take account of
!    the 3 different arithmetic modes.
!
!    How about an option to generate a random test problem,
!    similar to the sample?
!
!    09 May 1994:
!
!    Fixed DEC_PRINT.
!    Fixed RELDEC and DEC_TO_R.
!    Now, CHLDEC is printing 0.4 as 0.39999999999, and
!    DECREA is reading too many digits, and overflowing.
!
!    08 May 1994:
!
!    Wrote DBLDEC.  Updated DEC_MUL.
!    Rewrote FORM.  Created new RELDEC, DEC_TO_R, RATDEC, DEC_TO_RAT.
!
!    I checked all the lines where IFORM.EQ.2 occurred,
!    and tried to correct them, but gave up in the LP routines.
!    Check them later!
!
!    I still need to write CHLDEC and fix up DEC_PRINT.
!
!    07 May 1994:
!
!    Rewrote DEC_DIV to account for limited number of decimal digits.
!
!    Added "N" command allowing user to set number of decimal digits.
!
!    06 May 1994:
!
!    Proposal: the decimals should be stored
!    as SMANT (an integer representing the signed mantissa) and
!    IEXP (an integer representing the power of 10).
!
!    All calculations should be carried out by converting to a
!    real value first.  So I need to write just two routines,
!    RELDEC and DEC_TO_R.  Oh, and I forgot to mention that the
!    restriction on the number of decimal digits simply restricts 
!    the size of SMANT.
!
!    05 May 1994:
!
!    Modified AUTO_ERO to use SCADIV rather than
!    SCAMUL.  However, I still have something to complain about.
!    The interim values in the matrix (using decimals)
!    are not themselves decimals.  Once that happens,
!    the whole point is lost.  Is this DEC_DIV's fault?
!
!    03 May 1994:
!
!    What's with this DEC_PRINT routine?
!
!    03 May 1994:
!
!    Fixed error found on 23 April.  It was a tiny mistake in RELDEC.
!
!    23 April 1994:
!
!    Error:  I specified "B" for sample problem,
!    specified "F" and converted to "Decimal", chose "1" digit.
!    4 by 4 matrix got divided by 10, while RHS was correctly rounded.
!
!    Fixed a mistake in FORM which meant that conversions were
!    not being done when going from real to rational.
!
!    Replaced all labeled DO loops by DO/ENDDO pairs.
!
!    Added a "D" command that does division, renaming old
!    D (disk file) command to "K".
!
!    Version 1.51  19 September 1993
!
!    Changed from WRITE(* to WRITE(6 so that output redirection
!    works on DOS machines.
!
!    Also changed READ(* to READ(5.
!
!    I removed all occurrences of real(...), to make it easier
!    to convert the code to double precision if desired.
!    To convert this code to double precision, the only change
!    needed is a global substitute of "double precision" for
!    "real".
!
!    Version 1.50  06 June 1993
!
!    Added a sample linear system solve problem.
!
!    Dropped the "N" command.
!
!    Reversed order of arguments in call to S_WRITE.
!
!    Version 1.49  04 May 1993
!
!    Well guess what, the new version of Language Systems FORTRAN
!    seems to have fixed that bug!
!
!    Renamed SCALE to SCALER, because of a possible conflict
!    with a Macintosh routine (like this is going to help!).
!
!    A continuing bug that occurs on the Macintosh has led me
!    to add all sorts of checks for "NULL" characters in strings.
!    Right now it's just a hunch.
!
!    The "Z" command would claim an error had occurred if you
!    had not yet set up a matrix.
!
!    Changed S_WRITE to print an explicit blank as carriage control
!    on the Mac, since Language Systems "FORTRAN" will otherwise
!    print a null.
!
!    To make life easier for the "<" command, I dropped the
!    initial demand for arithmetic specification.
!
!    Added "<" command to allow user to specify an input file.
!    This is because it's hard to do on VMS via system commands,
!    and impossible on a Macintosh.
!
!    Added autoprint after pivoting.
!
!    Program should no longer fail if using rational arithmetic
!    and an overflow or underflow occurs.  Right now, MATMAN will
!    catch this problem, and halt the computation.  A better
!    solution would allow the user to request that overflows and
!    underflows be "rounded" to decimal and recomputed as ratios.
!
!    Noticed that Macintosh requires FORTRAN carriage control (1X)
!    for output to console, so had to modify S_WRITE.
!
!    Modified advanced LP problem, and corrected label.
!
!    Added S_LEFT, and forced S_READ to flush the string left
!    once it has been read.  This is so that, if I like, I can
!    type "B S" and have it mean the same as "BS".
!
!    Renamed SETSOL to LP_SOL.
!
!    Modified I_READ to have IHUSH parameter, so that I can
!    type a "Q" to quit in the Jacobi iteration.
!
!    Added the # feature, which allows a one line comment
!    beginning with the sharp symbol.
!
!    Added the $ and % commands, which allows me to turn paging off
!    and back on.
!
!    A row or column can be added to the matrix, or deleted
!    from it, in linear algebra mode.  An added row or column
!    can be inserted anywhere in the matrix.
!
!    CGC requested automatic printout of a matrix or table
!    when it is entered, and I have added that.
!
!    In linear algebra mode, the "o" command will now check
!    whether the matrix is in row echelon or reduced row echelon
!    form.
!
!    I added LEQI to this program, to avoid having to capitalize
!    everything.
!
!    Set up new paging routine, which, if it works, will be
!    added to MATALG also.  This also allows me to drop that
!    stupid ICOMP parameter from S_WRITE.
!
!    Cleaned up the main program a bit, so that the commands
!    are all part of one big IF/ELSIF block.
!
!    Ran the CLEANER program, to standardize indentation and
!    statement numbering.
!
!    Version 1.48  16 April 1993
!
!    Ran the STRIPPER program, to lower case all statements.
!
!    Removed last argument of S_READ.
!
!    Replaced old version of CHRCTR by a newer, better one which
!    does not suffer from integer "wrap around" when a large number
!    of decimal places are entered.
!
!    Version 1.47  25 April 1992
!
!    Cleaning up program format:
!    Continuation character is now always "&".
!    Marked the beginning of each routine with a "C****..." line.
!    Routines placed in alphabetical order.
!    Declared all variables.
!    Set SOL, ISLTOP and ISLBOT to be vectors, rather than 2D
!    arrays with a first dimension of 1.
!    Made JACOBI automatically iterative.
!    JACOBI will now accept nonsymmetric matrices.
!    CHRCTI and CHRCTR reset LCHAR = 0 on error.
!
!    Minor modification, 25 September 1991
!
!    Dropped "MAXCOL" from the argument list of JACOBI.
!
!    Version 1.46  01 February 1991
!
!    Tried a modification in LOOP 20 in HLPVMS.
!
!    Version 1.45  29 January 1991
!
!    Corrected error in example file reading.
!
!    Version 1.44  24 December 1990
!
!    Moved initializations to INIT.
!
!    Version 1.43  07 December 1990
!
!    ROWADD will now refuse to add a multiple of a row to itself.
!
!    Version 1.42  03 December 1990
!
!    Added Jacobi example.
!
!    CGC complained that "1" in the P column was missing, for
!    problems with artificial variables.  I tried to restore
!    that.
!
!    CGC complained that in problems with artificial variables, the
!    artificial objective function picks up the constant term of
!    the original objective.
!
!    Also added sample problem with artificial variables.
!
!    Added forced typeout of matrix after each Jacobi step.
!
!    Version 1.41  06 November 1990
!
!    FR and FF commands force real and rational arithmetic.
!
!    Added JACOBI routine to do Jacobi rotations.
!
!    Moved "hello" stuff to routine HELLO.
!
!    Moved transcript stuff to a routine TRANSC.
!
!    Version 1.40  15 October 1990
!
!    Modification made to allow matrix to be entered in entirety,
!    rather than a row at a time.
!
!    Version 1.39  12 October 1990
!
!    S_READ modified so that output need not be capitalized.  This
!    keeps filenames from being forcibly capitalized.
!
!    Version 1.38
!
!    * command allows you to transpose the matrix, LA mode only.
!
!    Version 1.37
!
!    In non linear programming mode, you may enter the entire
!    matrix, or several rows at a time, on one line, if you like.
!
!    Version 1.36
!
!    Updated interface to CHRPAK.
!
!    Inserted obsolete CHRPAK routines into MATMAN source.
!
!    Version 1.35
!
!    Correctly writes out hidden objective row for problems using
!    artificial variables.
!
!    Version 1.34
!
!    Corrected a transposition of variables that meant that, for
!    linear programming problems, constraints and variables were
!    interchanged.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 October 2023
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: maxcol = 30
  integer, parameter :: maxrow = 16

  real ( kind = rk ) a(maxrow,maxcol)
  logical autop
  real ( kind = rk ) b(maxrow,maxcol)
  real ( kind = rk ) c(maxrow,maxcol)
  logical c_is_digit
  character chineq(maxrow)
  character ( len = 20 ) command
  character ( len = 20 ) comold
  character ( len = 255 ) filex
  character ( len = 255 ) file_help
  character ( len = 255 ) filinp
  character ( len = 255 ) file_tran
  integer i_temp
  integer iabot(maxrow,maxcol)
  integer iatop(maxrow,maxcol)
  integer iauthr
  integer iauto
  integer ibase(maxrow)
  integer ibaseb(maxrow)
  integer ibasec(maxrow)
  integer ibbot(maxrow,maxcol)
  integer ibtop(maxrow,maxcol)
  integer icbot(maxrow,maxcol)
  integer icol
  integer icol1
  integer icol2
  integer ictop(maxrow,maxcol)
  integer ierror
  integer iform
  integer ihi
  integer ilo
  integer imat
  integer iopti
  integer iounit(4)
  integer iprint
  integer irow
  integer irow1
  integer irow2
  character isay
  integer isbot
  integer iseed
  integer islbot(maxcol)
  integer isltop(maxcol)
  integer istop
  integer iterm
  integer jhi
  integer jform
  integer jlo
  integer lens
  character ( len = 255 ) line
  character ( len = 255 ) line2
  integer lpmoda
  integer lpmodb
  integer lpmodc
  integer nart
  integer nartb
  integer nartc
  integer ncol
  integer ncolb
  integer ncolc
  integer ncon
  integer nrow
  integer nrowb
  integer nrowc
  integer n_slack
  integer n_slackb
  integer n_slackc
  integer nvar
  integer nvarb
  integer nvarc
  character ( len = 255 ) output
  integer output_page_length
  character ( len = 255 ) prompt
  real ( kind = rk ) rmat1(maxrow,maxcol)
  real ( kind = rk ) rmat2(maxrow,maxcol)
  logical s_eqi
  real ( kind = rk ) sol(maxcol)
  real ( kind = rk ) sval
  character ( len = 255 ) title
!
!  Initializations.
!
  call init ( a, autop, chineq, command, comold, filex, file_help, filinp, &
    file_tran, iabot, iatop, iauthr, ibase, ierror, iform, imat, &
    iounit, iseed, islbot, isltop, line, lpmoda, maxcol, maxrow, nart, ncol, &
    ncon, nrow, n_slack, nvar, sol )
 
  call mat_copy ( a, b, iatop, iabot, ibtop, ibbot, ibase, ibaseb, lpmoda, &
    lpmodb, maxcol, maxrow, nart, nartb, ncol, ncolb, nrow, nrowb, n_slack, &
    n_slackb, nvar, nvarb )
 
  call mat_copy ( a, c, iatop, iabot, ictop, icbot, ibase, ibasec, lpmoda, &
    lpmodc, maxcol, maxrow, nart, nartc, ncol, ncolc, nrow, nrowc, n_slack, &
    n_slackc, nvar, nvarc )
!
!  Say hello.
!
  call hello ( iounit )
!
!  Get the next command from the user.
!
10    continue
 
  iprint = 0
 
  if ( ierror /= 0 ) then

    output = ' '
    call s_write ( iounit, output )
    output = 'MATMAN could not carry out your command:'
    call s_write ( iounit, output )
    lens = len_trim ( command )
    output = '"' // command(1:lens) // '"'

    ierror = 0
!
!  Wipe out the offending command line.
!
    line = ' '
 
    if ( iounit(1) /= 0 ) then
      close ( unit = iounit(1) )
      iounit(1) = 0
      output = ' '
      call s_write ( iounit, output )
      output = 'Because an error occurred, we are closing'
      call s_write ( iounit, output )
      output = 'the input file, and requiring you to respond'
      call s_write ( iounit, output )
      output = 'directly!'
      call s_write ( iounit, output )
    end if
 
  end if

  line = ' '
!
!  Insert a blank line.
!
  if ( command /= '#' ) then
    output = ' '
    call s_write ( iounit, output )
  end if
!
!  Save the name of the previous command as COMOLD, in case we need
!  to undo it.  But only save "interesting" commands.
!
  if ( .not. ( &
    s_eqi ( command, 'H' ) .or. &
    s_eqi ( command, 'HELP' ) .or. & 
    s_eqi ( command, 'N' ) .or. &
    s_eqi ( command, 'O' ) .or. &
    s_eqi ( command, 'S' ) .or. &
    s_eqi ( command, 'T' ) .or. &
    s_eqi ( command, '$' ) .or. &
    s_eqi ( command, '?' ) .or. &
    s_eqi ( command, '%' ) .or. &
    s_eqi ( command, '<' ) .or. &
    s_eqi ( command, '#' ) ) ) then
    comold = command
  end if
 
  if ( command /= '#' ) then
    prompt = 'command? ("H" for help)'
  else
    prompt = ' '
  end if
!
!  No check for terminators.
!
  iterm = 0
  call s_read ( line2, line, prompt, iounit, ierror, iterm )

  line = line2
 
  if ( line2 == ' ' ) then
    go to 10
  end if
 
  if ( ierror /= 0 ) then
    ierror = 0
    command = 'Q'
  end if
!
!  Check to see if the command is an ERO, in a special format.
!
  call s_blank_delete ( line2 )

  if ( s_eqi ( line2, 'R' ) ) then

  else if ( s_eqi ( line2(1:1), 'R' ) .and. &
    ( line2(2:2) == ' ' .or. c_is_digit ( line2(2:2) ) ) ) then

    call row_op_check ( command, ierror, iounit, line2 )

    if ( ierror /= 0 ) then
      go to 10
    end if

    line = line2

  else if ( s_eqi ( line2(1:3), 'ROW' ) .and. &
    ( line2(4:4) == ' ' .or. c_is_digit ( line2(4:4) ) ) ) then

    call row_op_check ( command, ierror, iounit, line2 )
    if ( ierror /= 0 ) then
      go to 10
    end if

    line = line2
!
!  Check to see if the command is an ECO, in a special format.
!
  else if ( s_eqi ( line2, 'C' ) ) then

  else if ( s_eqi ( line2(1:1), 'C' ) .and. &
    ( line2(2:2) == ' ' .or. c_is_digit ( line2(2:2) ) ) ) then

    call col_op_check ( command, ierror, iounit, line2 )

    if ( ierror /= 0 ) then
      go to 10
    end if

    line = line2

  else if ( s_eqi ( line2(1:3), 'COL' ) .and. &
    ( line2(4:4) == ' ' .or. c_is_digit ( line2(4:4) ) ) ) then

    call col_op_check ( command, ierror, iounit, line2 )
    if ( ierror /= 0 ) then
      go to 10
    end if

    line = line2

  else if ( s_eqi ( line2(1:6), 'COLUMN' ) .and. &
    ( line2(7:7) == ' ' .or. c_is_digit ( line2(7:7) ) ) ) then

    call col_op_check ( command, ierror, iounit, line2 )
    if ( ierror /= 0 ) then
      go to 10
    end if

    line = line2
!
!  Maybe this is all I need in order to get my change command.
!
  else if ( s_eqi(line2(1:2), 'A(' ) ) then

    command = 'A('
    line = line2(3:)

  else

    command = ' '

  end if
!
!  If command was not an ERO that had to be translated, read it
!  the regular way.
!
!  Blank, slash, comma, semicolon, equals terminate the command.
!
  if ( command == ' ' ) then

    iterm = 1
    call s_read ( command, line, prompt, iounit, ierror, iterm )

    if ( ierror /= 0 ) then
      ierror = 0
      command = 'Q'
    end if

  end if
!
!  If the "ROW_AUTO" command was issued, the user must give authorization
!  the first time.
!
  if ( ( s_eqi ( command, 'ROW_AUTO' ) .or. s_eqi ( command, 'COL_AUTO' ) ) &
    .and. iauthr == 0 ) then


    if ( iauthr == 0 .or. imat == 0 ) then
      go to 10
    end if

  end if
!
!  Jump here when one command needs to switch to another.
!
20    continue
!
!  Save a copy of the matrix A in B before the operation, but only
!  for certain commands.
!
  if ( s_eqi ( command, 'A(' ) .or. s_eqi ( command, 'B' ) .or. &
       s_eqi ( command, 'BASIC' ) .or. s_eqi ( command, 'COL_ADD' ) .or. &
       s_eqi ( command, 'COL_AUTO' ) .or.s_eqi ( command, 'COL_DIV' ) .or. &
       s_eqi ( command, 'COL_MUL' ) .or. s_eqi ( command, 'COL_SWAP' ) .or. &
       s_eqi ( command, 'E' ) .or. s_eqi ( command, 'G' ) .or. &
       s_eqi ( command, 'J' ) .or. s_eqi ( command, 'L' ) .or. &
       s_eqi ( command, 'P' ) .or. s_eqi ( command, 'RES' ) .or. &
       s_eqi ( command, 'ROW_ADD' ) .or.  s_eqi ( command, 'ROW_AUTO' ) .or. &
       s_eqi ( command, 'ROW_DIV' ) .or.  s_eqi ( command, 'ROW_MUL' ) .or. &
       s_eqi ( command, 'ROW_SWAP' ) .or.  s_eqi ( command, 'TR' ) .or. &
       s_eqi ( command, 'V' ) .or. s_eqi ( command, 'X' )  ) then
  
    call mat_copy ( a, b, iatop, iabot, ibtop, ibbot, ibase, ibaseb, &
      lpmoda, lpmodb, maxcol, maxrow, nart, nartb, ncol, ncolb, &
      nrow, nrowb, n_slack, n_slackb, nvar, nvarb )
 
  end if
!
!  A(I,J)=X
!
  if ( s_eqi ( command(1:2), 'A(' ) ) then
 
    if ( imat /= 1 ) then

      ierror = 1
      output = 'No matrix has been defined yet!'
      call s_write ( iounit, output )

    else

      call change ( a, iatop, iabot, ierror, iform, iounit, line, maxcol, &
        maxrow, ncol, nrow )
 
      iprint = 1

    end if
!
!  B=Set up sample problem.
!
  else if ( s_eqi ( command, 'B' ) ) then
 
    call mat_zero ( a, iabot, iatop, iform, maxcol, maxrow )

    call sample ( a, chineq, iatop, iabot, ibase, ierror, iform, imat, &
      iounit, line, lpmoda, maxcol, maxrow, nart, ncol, nrow, n_slack, &
      nvar, rmat1 )
 
    if ( ierror == 0 ) then
 
      call mat_copy ( a, c, iatop, iabot, ictop, icbot, ibase, ibasec, &
        lpmoda, lpmodc, maxcol, maxrow, nart, nartc, ncol, ncolc, nrow, &
        nrowc, n_slack, n_slackc, nvar, nvarc )
 
      iprint = 1
 
  end if
!
!  BASIC = Assign row I to basic variable J.
!
  else if ( s_eqi ( command, 'BASIC' ) ) then
 
    if ( imat /= 1 ) then

      ierror = 1
      output = 'No matrix has been defined yet!'
      call s_write ( iounit, output )

    else if ( lpmoda /= 1 ) then

    else

      call basic ( ibase, ierror, iounit, line, maxrow, nart, nrow, n_slack, &
        nvar )

    end if
!
!  CHECK = Echelon Form / Optimality check.
!
  else if ( s_eqi ( command, 'CHECK' ) ) then
 
    if ( imat /= 1 ) then

      ierror = 1
      output = 'No matrix has been defined yet!'
      call s_write ( iounit, output )

    else if ( lpmoda == 0 ) then
 
      call la_opt ( a, iabot, iatop, iform, iounit, maxcol, maxrow, ncol, nrow )

    else

      call lp_opt ( a, iatop, iabot, ibase, iform, iopti, iounit, isltop, &
        islbot, lpmoda, maxcol, maxrow, nart, ncol, nrow, n_slack, nvar, sol )

    end if
!
!  COL_ADD=Add a multiple of one column to another.
!
  else if ( s_eqi ( command, 'COL_ADD' ) ) then
 
    if ( imat /= 1 ) then

      ierror = 1
      output = 'No matrix has been defined yet!'
      call s_write ( iounit, output )

    else

      call col_add_param ( ierror, iform, iounit, icol1, icol2, istop, isbot, &
        line, ncol, sval )
 
      if ( ierror == 0 ) then
 
        call col_add ( a, iatop, iabot, ierror, iform, iounit, icol1, icol2, &
          maxcol, maxrow, ncol, sval, istop, isbot )
 
        iprint = 1
 
      end if

    end if
!
!  COL_AUTO=Automatic column reduction.
!
  else if ( s_eqi ( command, 'COL_AUTO' ) ) then
 
    if ( imat /= 1 ) then

      ierror = 1
      output = 'No matrix has been defined yet!'
      call s_write ( iounit, output )

    else if ( lpmoda == 1 ) then
 
      ierror = 1
      output = 'Cannot do this in linear programming mode!'
      call s_write ( iounit, output )

    else

      call col_auto ( a, iatop, iabot, ierror, iform, iounit, maxcol, maxrow, &
        ncol, nrow )
 
    end if
 
    iprint = 1
!
!  COL_DIV = Divide column by scalar.
!
  else if ( s_eqi ( command, 'COL_DIV' ) ) then
 
    if ( imat /= 1 ) then

      ierror = 1
      output = 'No matrix has been defined yet!'
      call s_write ( iounit, output )

    else

      call col_div_param ( icol, ierror, iform, iounit, isbot, istop, line, &
        sval )

      call col_div ( a, iatop, iabot, ierror, iform, iounit, icol, maxcol, &
        maxrow, ncol, nrow, sval, istop, isbot )

      if ( ierror == 0 ) then
        iprint = 1
      end if

    end if
!
!  COL_MUL=Multiply column by scalar.
!
  else if ( s_eqi ( command, 'COL_MUL' ) ) then
 
    if ( imat /= 1 ) then

      ierror = 1
      output = 'No matrix has been defined yet!'
      call s_write ( iounit, output )

    else

      call col_mul_param ( ierror, iform, iounit, icol, istop, isbot, line, &
        sval )
 
      if ( ierror == 0 ) then
 
        call col_mul ( a, iatop, iabot, ierror, iform, iounit, icol, maxcol, &
          maxrow, ncol, nrow, sval, istop, isbot )
 
        iprint = 1

      end if
 
    end if
!
!  COL_SWAP = Swap columns I and J.
!
  else if ( s_eqi ( command, 'COL_SWAP' ) ) then
 
    if ( imat == 0 ) then

      ierror = 1
      output = 'You must set up a matrix first!'
      call s_write ( iounit, output )

    else

      prompt = 'column I, column J.'
      call i_read ( icol1, line, prompt, iounit, ierror )
      if ( ierror /= 0 ) then
        go to 10
      end if
 
      call i_read ( icol2, line, prompt, iounit, ierror )
      if ( ierror /= 0 ) then
        go to 10
      end if
 
      call col_swap ( a, iatop, iabot, ierror, iform, iounit, icol1, icol2, &
        maxcol, maxrow, ncol, nrow )
 
      iprint = 1

    end if
!
!  DEC_DIGIT=Set number of digits.
!
  else if ( s_eqi ( command, 'DEC_DIGIT' ) ) then
 
    call dec_digit_set ( ierror, iounit, line )
!
!  DECimal = use decimal arithmetic
!
  else if ( s_eqi ( command(1:3), 'DEC' ) ) then
 
    jform = 2
 
    call form ( a, b, c, iatop, iabot, ibtop, ibbot, ictop, icbot, iform, &
      imat, iounit, jform, maxcol, maxrow )
 
    iprint = 1
!
!  E=Enter problem definition
!
  else if ( s_eqi ( command, 'E' ) ) then

    call mat_zero ( a, iabot, iatop, iform, maxcol, maxrow )
 
    if ( lpmoda == 0 ) then

      nvar = 0
      call la_inp0 ( ierror, iounit, line, maxcol, maxrow, ncol, nrow )
 
      call la_inp1 ( a, iabot, iatop, ierror, iform, iounit, line, maxcol, &
        maxrow, 1, nrow, 1, ncol )
 
    else
 
      call lp_inp ( a, chineq, iatop, iabot, ibase, ierror, iform, iounit, &
        line, maxcol, maxrow, nart, ncol, ncon, nrow, n_slack, nvar )
 
    end if
 
    if ( iounit(1) == 41 ) then
      close(unit = iounit(1))
      iounit(1) = 0
      output = 'The example has been read.'
      call s_write ( iounit, output )
    end if
 
    if ( ierror /= 0 ) then
      go to 10
    end if
 
    imat = 1
 
    call mat_copy ( a, c, iatop, iabot, ictop, icbot, ibase, ibasec, &
      lpmoda, lpmodc, maxcol, maxrow, nart, nartc, ncol, ncolc, &
      nrow, nrowc, n_slack, n_slackc, nvar, nvarc )
 
    output = 'A copy of this matrix is being saved.'
    call s_write ( iounit, output )
 
    iprint = 1
!
!  G=Add/delete a row or column of the matrix,
!    Add a constraint to the table.
!
  else if ( s_eqi ( command, 'G' ) ) then
 
    if ( imat /= 1 ) then

      ierror = 1
      output = 'No matrix has been defined yet!'
      call s_write ( iounit, output )

    else

      call deladd ( a, iabot, iatop, ibase, ierror, iform, iounit, line, &
        lpmoda, maxcol, maxrow, ncol, ncon, nrow, n_slack, nvar )
 
      iprint = 1

    end if
!
!  H=Help.
!
  else if ( s_eqi ( command, 'H' ) ) then
 
    if ( lpmoda == 0 ) then
      call la_help ( iounit )
    else
      call lp_help ( iounit )
    end if
 
  else if ( s_eqi ( command(1:4), 'HELP' ) ) then
 
    call help ( iounit )
!
!  I_BIG = Set maximum integer.
!
  else if ( s_eqi ( command, 'I_BIG' ) ) then
 
    prompt = 'maximum integer for rational representations.'
    call s_blanks_delete ( prompt )
 
    call i_read ( i_temp, line, prompt, iounit, ierror )

    call i_data ( 'SET', 'I_BIG', i_temp )
!
!  INIT = Initializations.
!
  else if ( s_eqi ( command(1:3), 'INIT' ) ) then

    call init ( a, autop, chineq, command, comold, filex, file_help, filinp, &
      file_tran, iabot, iatop, iauthr, ibase, ierror, iform, imat, &
      iounit, iseed, islbot, isltop, line, lpmoda, maxcol, maxrow, nart, ncol, &
      ncon, nrow, n_slack, nvar, sol )
 
    call mat_copy ( a, b, iatop, iabot, ibtop, ibbot, ibase, ibaseb, &
      lpmoda, lpmodb, maxcol, maxrow, nart, nartb, ncol, ncolb, &
      nrow, nrowb, n_slack, n_slackb, nvar, nvarb )
 
    call mat_copy ( a, c, iatop, iabot, ictop, icbot, ibase, ibasec, &
      lpmoda, lpmodc, maxcol, maxrow, nart, nartc, ncol, ncolc, &
      nrow, nrowc, n_slack, n_slackc, nvar, nvarc )
!
!  J=Jacobi pre and post multiplication by (I,J) plane rotation.
!
  else if ( s_eqi ( command, 'J' ) ) then
 
    if ( imat /= 1 ) then

      ierror = 1
      output = 'No matrix has been defined yet!'
      call s_write ( iounit, output )

    else if ( lpmoda /= 0 ) then

    else if ( iform /= 1 ) then

    else

      call r_jacobi ( a, ibase, ierror, iounit, line, lpmoda, maxcol, maxrow, &
        ncol, nrow, rmat1, rmat2 )
 
      iprint = 0

    end if
!
!  K=Disk file is to be opened or closed.
!
  else if ( s_eqi ( command, 'K' ) ) then
 
    call transc ( file_tran, ierror, iounit, line )
!
!  L=Change between linear algebra and linear programming modes.
!
  else if ( s_eqi ( command, 'L' ) ) then
 
    call lp_set ( ierror, imat, iounit, line, lpmoda, nart, ncol, ncon, nrow, &
      n_slack, nvar )
!
!  O=Orthogonal transformation to zero out A(I,J)..
!
  else if ( s_eqi ( command, 'O' ) ) then
 
    if ( imat == 0 ) then

      ierror = 1
      output = 'You must set up a matrix first!'
      call s_write ( iounit, output )

    else if ( iform /= 1 ) then

      ierror = 1
      output = 'You must use REAL arithmetic!'
      call s_write ( iounit, output )

    else if ( lpmoda == 1 ) then

      ierror = 1
      output = 'You cannot be in linear programming mode!'
      call s_write ( iounit, output )

    else

      prompt = 'row I, column J.'
      call i_read ( irow, line, prompt, iounit, ierror )
      if ( ierror /= 0 ) then
        go to 10
      end if
 
      call i_read ( icol, line, prompt, iounit, ierror )
      if ( ierror /= 0 ) then
        go to 10
      end if
 
      call orth ( a, ierror, iounit, irow, icol, maxcol, maxrow, ncol, nrow )
 
      iprint = 1

    end if
!
!  P=Pivot.
!
  else if ( s_eqi ( command, 'P' ) ) then
 
    if ( imat /= 1 ) then

    else

      iauto = 0
 
      call lp_piv ( a, iatop, iabot, iauto, ibase, ierror, iform, iounit, &
        isltop, islbot, line, lpmoda, maxcol, maxrow, nart, ncol, nrow, &
        n_slack, nvar, sol )
 
      iprint = 1

    end if
!
!  Q = Quit (after confirmation).
!  QUit = QUIT NOW!
!
  else if ( s_eqi ( command(1:1), 'Q' ) ) then
 
    if ( s_eqi ( command(2:2), 'U' ) ) then
      isay = 'Y'
    else
      line = ' '
      prompt = '"Y" to confirm you want to quit.'
      iterm = 0
      call s_read ( isay, line, prompt, iounit, ierror, iterm )
      if ( ierror /= 0) then
        isay = 'Y'
      end if
    end if
 
    if ( s_eqi ( isay, 'Y' ) ) then
      output = ''
      call s_write ( iounit, output )
      output = 'MATMAN():'
      call s_write ( iounit, output )
      output = '  Normal end of execution.'
      call s_write ( iounit, output )
      call timestring ( output )
      call s_write ( iounit, output )

      if ( iounit(3) /= -1 ) then
        call transc ( file_tran, ierror, iounit, line )
      end if
 
      stop
    end if
!
!  RATional = use rational arithmetic
!
  else if ( s_eqi ( command(1:3), 'RAT' ) ) then
 
    jform = 0
 
    call form ( a, b, c, iatop, iabot, ibtop, ibbot, ictop, icbot, iform, &
      imat, iounit, jform, maxcol, maxrow )
 
    iprint = 1
!
!  REAl = use real arithmetic
!
  else if ( s_eqi ( command(1:3), 'REA' ) ) then
 
    jform = 1
 
    call form ( a, b, c, iatop, iabot, ibtop, ibbot, ictop, icbot, iform, &
      imat, iounit, jform, maxcol, maxrow )
 
    iprint = 1
!
!  RES=Restore matrix.
!
  else if ( s_eqi ( command(1:3), 'RES' ) ) then
 
    call mat_restore ( a, c, iabot, iatop, ibase, ibasec, icbot, ictop, &
      ierror, imat, iounit, lpmoda, lpmodc, maxcol, maxrow, nart, nartc, &
      ncol, ncolc, nrow, nrowc, n_slack, n_slackc, nvar, nvarc )
 
    if ( ierror == 0 ) then
      iprint = 1
    end if
!
!  ROW_ADD=Add a multiple of one row to another.
!
  else if ( s_eqi ( command, 'ROW_ADD' ) ) then
 
    if ( imat /= 1 ) then

      ierror = 1
      output = 'No matrix has been defined yet!'
      call s_write ( iounit, output )

    else

      call row_add_param ( ierror, iform, iounit, irow1, irow2, istop, isbot, &
        line, nrow, sval )
 
      if ( ierror == 0 ) then
 
        call row_add ( a, iatop, iabot, ierror, iform, iounit, irow1, irow2, &
          maxcol, maxrow, ncol, sval, istop, isbot )
 
        iprint = 1
 
      end if

    end if
!
!  ROW_AUTO=Automatic row reduction.
!
  else if ( s_eqi ( command, 'ROW_AUTO' ) ) then
 
    if ( imat /= 1 ) then

      ierror = 1
      output = 'No matrix has been defined yet!'
      call s_write ( iounit, output )

    else if ( lpmoda == 0 ) then
 
      call row_auto ( a, iatop, iabot, ibase, ierror, iform, iounit, maxcol, &
        maxrow, ncol, nrow )
 
    else if ( lpmoda == 1 ) then
 
      iauto = 1
 
      call lp_piv ( a, iatop, iabot, iauto, ibase, ierror, iform, iounit, &
        isltop, islbot, line, lpmoda, maxcol, maxrow, nart, ncol, nrow, &
        n_slack, nvar, sol )

    end if
 
    iprint = 1
!
!  ROW_DIV = Divide row by scalar.
!
  else if ( s_eqi ( command, 'ROW_DIV' ) ) then
 
    if ( imat /= 1 ) then

      ierror = 1
      output = 'No matrix has been defined yet!'
      call s_write ( iounit, output )

    else

      call row_div_param ( ierror, iform, iounit, irow, isbot, istop, line, &
        sval )

      call row_div ( a, iatop, iabot, ierror, iform, iounit, irow, maxcol, &
        maxrow, ncol, nrow, sval, istop, isbot )
 
      if ( ierror == 0 ) then
        iprint = 1
      end if

    end if
!
!  ROW_MUL=Multiply row by scalar.
!
  else if ( s_eqi ( command, 'ROW_MUL' ) ) then
 
    if ( imat /= 1 ) then

      ierror = 1
      output = 'No matrix has been defined yet!'
      call s_write ( iounit, output )

    else

      call row_mul_param ( ierror, iform, iounit, irow, istop, isbot, line, &
        sval )
 
      if ( ierror == 0 ) then
 
        call row_mul ( a, iatop, iabot, ierror, iform, iounit, irow, maxcol, &
          maxrow, ncol, nrow, sval, istop, isbot )
 
        iprint = 1

      end if
 
    end if
!
!  ROW_SWAP = Interchange rows I and J.
!
  else if ( s_eqi ( command, 'ROW_SWAP' ) ) then
 
    if ( imat == 0 ) then

      ierror = 1
      output = 'You must set up a matrix first!'
      call s_write ( iounit, output )

    else

      prompt = 'row I, row J.'
      call i_read ( irow1, line, prompt, iounit, ierror )
      if ( ierror /= 0 ) then
        go to 10
      end if
 
      call i_read ( irow2, line, prompt, iounit, ierror )
      if ( ierror /= 0 ) then
        go to 10
      end if
 
      call row_swap ( a, iatop, iabot, ibase, ierror, iform, iounit, &
        irow1, irow2, lpmoda, maxcol, maxrow, ncol, nrow )
 
      iprint = 1

    end if
!
!  S=Store a matrix.
!
  else if ( s_eqi ( command, 'S' ) ) then
 
    if ( imat /= 1 ) then

      ierror = 1
      output = 'No matrix has been defined yet!'
      call s_write ( iounit, output )

    else

      call mat_copy ( a, c, iatop, iabot, ictop, icbot, ibase, &
        ibasec, lpmoda, lpmodc, maxcol, maxrow, nart, nartc, ncol, &
        ncolc, nrow, nrowc, n_slack, n_slackc, nvar, nvarc )
 
      output = 'A copy of the matrix has been stored.'
      call s_write ( iounit, output )

    end if
!
!  T=Type matrix or table.
!
  else if ( s_eqi ( command, 'T' ) ) then
 
    if ( imat /= 1 ) then

      ierror = 1
      output = 'No matrix has been defined yet!'
      call s_write ( iounit, output )

    else

      if ( lpmoda == 0 ) then
        title = 'The current matrix:'
      else if ( lpmoda == 1 ) then
        title = 'The linear programming table:'
      end if

      ilo = 1
      jlo = 1
      jhi = ncol
 
      if ( lpmoda == 1 .and. 0 < nart ) then
        ihi = nrow + 1
      else
        ihi = nrow
      end if

      call mat_print ( a, iabot, iatop, ibase, iform, iounit, ihi, ilo, jhi, &
        jlo, lpmoda, maxcol, maxrow, ncol, nrow, title )

    end if
!
!  TS = Type linear programming solution.
!
  else if ( command(1:2) == 'TS' ) then
 
    if ( imat /= 1 ) then

      ierror = 1
      output = 'No matrix has been defined yet!'
      call s_write ( iounit, output )

    else if ( lpmoda /= 1 ) then

      ierror = 1
      output = 'Error!  There is no linear programming solution'
      call s_write ( iounit, output )
      output = 'to print, because we are not in linear'
      call s_write ( iounit, output )
      output = 'programming mode!'
      call s_write ( iounit, output )

    else

      call lp_sol ( a, iatop, iabot, ibase, iform, isltop, islbot, &
        maxcol, maxrow, ncol, nrow, sol )

      title = 'The linear programming solution'

      call sol_print ( ibase, iform, iounit, islbot, isltop, lpmoda, maxcol, &
        maxrow, nart, ncol, nrow, n_slack, nvar, sol, title )

    end if
!
!  U=Undo last command.
!
  else if ( s_eqi ( command, 'U' ) ) then
 
    if ( s_eqi ( comold, 'K' ) ) then
      comold = 'U'
      command = 'K'
      go to 20
    end if
 
    if ( s_eqi ( comold, 'H' ) ) then
      go to 10
    end if

    if ( s_eqi ( comold, 'HELP' ) ) then
      go to 10
    end if

    if ( s_eqi ( comold, 'N' ) ) then
      go to 10
    end if
 
    if ( s_eqi ( comold, 'L' ) ) then
      comold = 'U'
      command = 'L'
      go to 20
    end if
 
    if ( s_eqi ( comold, 'O' ) ) then
      go to 10
    end if

    if ( s_eqi ( comold, 'T' ) ) then
      go to 10
    end if

    if ( s_eqi ( comold, 'W' ) ) then
      go to 10
    end if
 
    call mat_copy ( b, a, ibtop, ibbot, iatop, iabot, ibaseb, &
      ibase, lpmodb, lpmoda, maxcol, maxrow, nartb, nart, ncolb, &
      ncol, nrowb, nrow, n_slackb, n_slack, nvarb, nvar )
 
    iprint = 1
!
!  V=Remove artificial variables.
!
  else if ( s_eqi ( command, 'V' ) ) then
 
    if ( imat /= 1 ) then

      ierror = 1
      output = 'No matrix has been defined yet!'
      call s_write ( iounit, output )

    else

      call lp_rem ( a, iabot, iatop, ibase, ierror, iform, iounit, lpmoda, &
        maxcol, maxrow, nart, ncol, nrow, n_slack, nvar )
 
      iprint = 1

    end if
!
!  W=Write example to file.
!
  else if ( s_eqi ( command, 'W' ) ) then
 
    if ( imat /= 1 ) then

      ierror = 1
      output = 'No matrix has been defined yet!'
      call s_write ( iounit, output )

    else

      call file_write ( a, chineq, filex, iatop, iabot, ierror, iform, &
        iounit, line, lpmoda, maxcol, maxrow, nart, ncol, nrow, nvar )

    end if
!
!  X=Read example from file.
!
  else if ( s_eqi ( command, 'X' ) ) then
 
    call file_read ( filex, ierror, iform, iounit, line, lpmoda )
 
    if ( ierror /= 0 ) then
      go to 10
    end if
 
    command = 'E'
    go to 20
!
!  Y=Turn autoprint off or on.
!
  else if ( s_eqi ( command, 'Y' ) ) then
 
    autop = .not. autop
 
    if ( autop ) then
      output = 'Autoprinting turned ON.'
    else
      output = 'Autoprinting turned OFF.'
    end if

    call s_write ( iounit, output )
!
!  # = Comment.
!  Blank out the input line so MATMAN doesn't reparse it, looking for commands.
!
  else if ( command == '#' ) then

    line = ' '
!
!  $ sign means no paging.
!
  else if ( command == '$' ) then

    output_page_length = 0
    call i_data ( 'SET', 'OUTPUT_PAGE_LENGTH', output_page_length )
    output = 'Paging turned OFF.'
    call s_write ( iounit, output )
!
!  % means restore paging.
!
  else if ( command == '%' ) then
 
    prompt = 'number of lines to print before pausing.'
    call i_read ( output_page_length, line, prompt, iounit, ierror )
    if ( ierror /= 0 ) then
      go to 10
    end if
 
    call i_data ( 'SET', 'OUTPUT_PAGE_LENGTH', output_page_length )
    output = 'Paging turned ON.'
    call s_write ( iounit, output )
!
!  < means input from a file.
!
  else if ( command == '<' ) then
 
    call infile ( filinp, ierror, iounit, line )
!
! ? Extensive help from file.
!
  else if ( command == '?' ) then
 
    call get_help ( file_help, iounit, line )
!
!  No match!
!
  else if ( command /= ' ' ) then

    output = ' '
    call s_write ( iounit, output )
    output = 'MATMAN did not recogize your command:'
    call s_write ( iounit, output )
    lens = len_trim ( command )
    output = '  "' // command(1:lens) // '"'
    call s_write ( iounit, output )

    ierror = 1

  end if
!
!  After certain operations, print out the matrix.
!
  if ( autop .and. ierror == 0 .and. imat == 1 .and. iprint == 1 ) then
 
    if ( lpmoda == 0 ) then
      title = '  The current matrix:'
    else if ( lpmoda == 1 ) then
      title = '  The linear programming table:'
    end if

    ilo = 1
    jlo = 1
    jhi = ncol
 
    if ( lpmoda == 1 .and. 0 < nart ) then
      ihi = nrow + 1
    else
      ihi = nrow
    end if

    call mat_print ( a, iabot, iatop, ibase, iform, iounit, ihi, ilo, jhi, &
      jlo, lpmoda, maxcol, maxrow, ncol, nrow, title )
 
    iprint = 0
 
  end if

  go to 10
end
subroutine basic ( ibase, ierror, iounit, line, maxrow, nart, nrow, n_slack, &
  nvar )

!*****************************************************************************80
!
!! basic() assigns a row of the table to one of the basic variables.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer IBASE(MAXROW), the basic variables.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Workspace, character ( len = 255 ) LINE, used to hold the user's input.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
!    Input, integer NART, the number of artificial variables.
!
!    Input, integer NROW, the number of rows in the matrix.
!
!    Input, integer N_SLACK, the number of slack variables.
!
!    Input, integer NVAR, the number of basic variables.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxrow

  character ( len = 10 ) chrtmp1
  character ( len = 10 ) chrtmp2
  integer ibase(maxrow)
  integer ierror
  integer iounit(4)
  integer irow
  integer ivar
  character ( len = 255 ) line
  integer nart
  integer nrow
  integer n_slack
  integer nvar
  character ( len = 255 ) output
  character ( len = 255 ) prompt

  prompt = 'row I, basic variable J.'
!
!  Get the row number I.
!
  call i_read ( irow, line, prompt, iounit, ierror )
  if ( ierror /= 0 ) then
    return
  end if
 
  if ( irow < 1 .or. nrow < irow ) then
    output = 'Error!  Illegal row number!'
    call s_write ( iounit, output )
    ierror = 1
    return
  end if
!
!  Get the basic variable index J.
!
  call i_read ( ivar, line, prompt, iounit, ierror )
  if ( ierror /= 0 ) then
    return
  end if
 
  if ( ivar < 1 .or. nvar + n_slack + nart < ivar ) then
    output = 'Error!  Illegal basic variable number!'
    call s_write ( iounit, output )
    ierror = 1
    return
  end if
 
  ibase(irow) = ivar
 
  call i_to_s_left ( irow, chrtmp1 )
  call i_to_s_left ( ivar, chrtmp2 )
  output = 'Assigning row ' // chrtmp1 // ' to basic variable' // chrtmp2
  call s_blanks_delete ( output )
  call s_write ( iounit, output )
 
  return
end
subroutine c_cap ( c )

!*****************************************************************************80
!
!! c_cap() capitalizes a single character.
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

  integer, parameter :: rk = kind ( 1.0D+00 )

  character c
  integer itemp

  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end
function c_eqi ( c1, c2 )

!*****************************************************************************80
!
!! c_eqi() is a case insensitive comparison of two characters for equality.
!
!  Discussion:
!
!    C_EQI ( 'A', 'a' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C1, C2, the characters to compare.
!
!    Output, logical C_EQI, the result of the comparison.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  logical c_eqi
  character c1
  character c2
  character cc1
  character cc2

  cc1 = c1
  cc2 = c2

  call c_cap ( cc1 )
  call c_cap ( cc2 )

  if ( cc1 == cc2 ) then
    c_eqi = .true.
  else
    c_eqi = .false.
  end if

  return
end
function c_is_alpha ( c )

!*****************************************************************************80
!
!! c_is_alpha() returns TRUE if C is an alphabetic character.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, a character to check.
!
!    Output, logical C_IS_ALPHA is TRUE if C is an alphabetic character.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character c
  logical c_is_alpha

  if ( ( lle ( 'a', c ) .and. lle ( c, 'z' ) ) .or. &
       ( lle ( 'A', c ) .and. lle ( c, 'Z' ) ) ) then
    c_is_alpha = .true.
  else
    c_is_alpha = .false.
  end if

  return
end
function c_is_digit ( c )

!*******************************************************************************
!
!! c_is_digit() returns .TRUE. if a character is a decimal digit.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the character to be analyzed.
!
!    Output, logical C_IS_DIGIT, .TRUE. if C is a digit, .FALSE. otherwise.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character c
  logical c_is_digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then
    c_is_digit = .true.
  else
    c_is_digit = .false.
  end if

  return
end
subroutine c_to_digit ( c, digit )

!*****************************************************************************80
!
!! c_to_digit() returns the integer value of a base 10 digit.
!
!  Discussion:
!
!     C   DIGIT
!    ---  -----
!    '0'    0
!    '1'    1
!    ...  ...
!    '9'    9
!    ' '    0
!    'X'   -1
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the decimal digit, '0' through '9' or blank
!    are legal.
!
!    Output, integer DIGIT, the corresponding integer value.  If C was
!    'illegal', then DIGIT is -1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character c
  integer digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

    digit = ichar ( c ) - 48

  else if ( c == ' ' ) then

    digit = 0

  else

    digit = -1

  end if

  return
end
subroutine change ( a, iatop, iabot, ierror, iform, iounit, line, maxcol, &
  maxrow, ncol, nrow )

!*****************************************************************************80
!
!! change() allows the user to change an entry in the array.
!
!  Discussion:
!
!    Expect an input line of the form:
!
!      A(5,12) = 17
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) A(MAXROW,MAXCOL), the current matrix.
!
!    Input/output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
!    IATOP and IABOT represent the rational or decimal matrix
!    whose entry is to be changed.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred.
!
!    Input, integer IFORM, 0/1/2, rational/real/decimal arithmetic.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Workspace, character ( len = 255 ) LINE, used to hold the user's input.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
!    Input, integer NCOL, the number of columns in the matrix.
!
!    Input, integer NROW, the number of rows in the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxcol
  integer maxrow

  real ( kind = rk ) a(maxrow,maxcol)
  character ( len = 10 ) chrtmp1
  character ( len = 10 ) chrtmp2
  character ( len = 22 ) chrtmp3
  integer i_gcd
  integer iabot(maxrow,maxcol)
  integer iatop(maxrow,maxcol)
  integer icol
  integer ierror
  integer iform
  integer iounit(4)
  integer irow
  integer isbot
  integer istop
  integer itemp
  integer len1
  integer len2
  integer len3
  character ( len = 255 ) line
  integer ncol
  integer nrow
  character ( len = 255 ) output
  character ( len = 255 ) prompt
  real ( kind = rk ) rval

  prompt = 'row I, column J, new value S.'
!
!  Get the row number.
!
  call i_read ( irow, line, prompt, iounit, ierror )

  if ( ierror /= 0 ) then
    output = '  I_READ returned error flag!'
    call s_write ( iounit, output )
    return
  end if
 
  if ( irow < 1 .or. nrow < irow ) then
    output = 'Error!  Illegal row value!'
    call s_write ( iounit, output )
    ierror = 1
    return
  end if
!
!  Get the column number.
!
  if ( line(1:1) == ',' ) then
    line(1:1) = ' '
  end if

  call i_read ( icol, line, prompt, iounit, ierror )

  if ( ierror /= 0 ) then
    output = '  I_READ returned error flag!'
    call s_write ( iounit, output )
    return
  end if
 
  if ( icol < 1 .or. ncol < icol ) then
    output = 'Error!  Illegal column value!'
    call s_write ( iounit, output )
    ierror = 1
    return
  end if
!
!  Read the value.
!
  if ( line(1:2) == ')=' ) then
    line(1:2) = ' '
  end if

  if ( iform == 0 ) then
 
    call rat_read ( istop, isbot, line, prompt, iounit, ierror )

    if ( ierror /= 0 ) then
      output = '  RAT_READ returned error flag!'
      call s_write ( iounit, output )
      return
    end if
 
    call rat_to_s_left ( istop, isbot, chrtmp3 )

    itemp = i_gcd ( istop, isbot )
    iatop(irow,icol) = istop / itemp
    iabot(irow,icol) = isbot / itemp
 
  else if ( iform == 1 ) then
 
    call r_read ( rval, line, prompt, iounit, ierror )

    if ( ierror /= 0 ) then
      output = '  R_READ returned error flag!'
      call s_write ( iounit, output )
      return
    end if

    a(irow,icol) = rval

    call r_to_s_left ( rval, chrtmp3 )
 
  else if ( iform == 2 ) then
 
    call dec_read ( istop, isbot, line, prompt, iounit, ierror )

    if ( ierror /= 0 ) then
      output = '  DEC_READ returned error flag!'
      call s_write ( iounit, output )
      return
    end if
 
    call dec_round ( istop, isbot )

    call dec_to_s_left ( istop, isbot, chrtmp3 )

    iatop(irow,icol) = istop
    iabot(irow,icol) = isbot
 
  end if
 
  call i_to_s_left ( irow, chrtmp1 )
  len1 = len_trim ( chrtmp1 )

  call i_to_s_left ( icol, chrtmp2 )
  len2 = len_trim ( chrtmp2 )

  len3 = len_trim ( chrtmp3 )

  output = ' '
  call s_write ( iounit, output )

  output = '  A(' // chrtmp1(1:len1) // ',' // chrtmp2(1:len2) // ') = ' &
    // chrtmp3(1:len3)
  call s_write ( iounit, output )

  return
end
subroutine chrctf ( string, itop, ibot, ierror, lchar )

!*****************************************************************************80
!
!! chrctf() reads an integer or rational fraction from a string.
!
!  Discussion:
!
!    The integer may be in real format, for example '2.25'.  It
!    returns ITOP and IBOT.  If the input number is an integer, ITOP
!    equals that integer, and IBOT is 1.  But in the case of 2.25,
!    the program would return ITOP = 225, IBOT = 100.
!
!    Legal input is
!
!      blanks,
!      initial sign,
!      blanks,
!      integer part,
!      decimal point,
!      fraction part,
!      'E' or 'e' or 'D' or 'd', exponent marker,
!      exponent sign,
!      exponent integer part,
!      blanks,
!      final comma or semicolon,
!
!    with most quantities optional.
!
!  Examples:
!
!    STRING            ITOP      IBOT
!
!    '1'               1         1
!    '     1   '       1         1
!    '1A'              1         1
!    '12,34,56'        12        1
!    '  34 7'          34        1
!    '-1E2ABCD'        -100      1
!    '-1X2ABCD'        -1        1
!    ' 2E-1'           2         10
!    '23.45'           2345      100
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRING, the string containing the
!    data to be read.  Reading will begin at position 1 and
!    terminate when no more characters
!    can be read to form a legal integer.  Blanks, commas,
!    or other nonnumeric data will, in particular, cause
!    the conversion to halt.
!
!    Output, integer ITOP, the integer read from the string,
!    assuming that no negative exponents or fractional parts
!    were used.  Otherwise, the 'integer' is ITOP/IBOT.
!
!    Output, integer IBOT, the integer divisor required to
!    represent numbers which are in real format or have a
!    negative exponent.
!
!    Output, integer IERROR, error flag.
!    0 if no errors,
!    Value of IHAVE when error occurred otherwise.
!
!    Output, integer LCHAR, number of characters read from
!    STRING to form the number.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  logical c_is_digit
  character chrtmp
  integer ibot
  integer ierror
  integer ihave
  integer isgn
  integer iterm
  integer itop
  integer jsgn
  integer jtop
  integer lchar
  integer nchar
  integer ndig
  logical s_eqi
  character ( len = * ) string

  nchar = len ( string )
 
  ierror = 0
  lchar = - 1
  isgn = 1
  itop = 0
  ibot = 1
  jsgn = 1
  jtop = 0
  ihave = 1
  iterm = 0

10    continue

  lchar = lchar + 1
  chrtmp = string(lchar+1:lchar+1)
!
!  Blank.
!
  if ( chrtmp == ' ' ) then
 
    if ( ihave == 2 ) then
 
    else if ( ihave == 6 .or. ihave == 7 ) then
      iterm = 1
    else if ( 1 < ihave ) then
      ihave = 11
    end if
!
!  Comma.
!
  else if ( chrtmp == ',' .or. chrtmp == ';' ) then
 
    if ( ihave /= 1 ) then
      iterm = 1
      ihave = 12
      lchar = lchar + 1
    end if
!
!  Minus sign.
!
  else if ( chrtmp == '-' ) then
 
    if ( ihave == 1 ) then
      ihave = 2
      isgn = - 1
    else if ( ihave == 6 ) then
      ihave = 7
      jsgn = - 1
    else
      iterm = 1
    end if
!
!  Plus sign.
!
  else if ( chrtmp == '+' ) then
 
    if ( ihave == 1 ) then
      ihave = 2
    else if ( ihave == 6 ) then
      ihave = 7
    else
      iterm = 1
    end if
!
!  Decimal point.
!
  else if ( chrtmp == '.' ) then
 
    if ( ihave < 4 ) then
      ihave = 4
    else
      iterm = 1
    end if
!
!  Exponent marker.
!
  else if ( s_eqi ( chrtmp, 'E' ) .or. s_eqi ( chrtmp, 'D' ) ) then
 
    if ( ihave < 6 ) then
      ihave = 6
    else
      iterm = 1
    end if
!
!  Digit.
!
  else if ( c_is_digit ( chrtmp ) .and. ihave < 11 ) then
 
    if ( ihave <= 2 ) then
      ihave = 3
    else if ( ihave == 4 ) then
      ihave = 5
    else if ( ihave == 6 .or. ihave == 7 ) then
      ihave = 8
    end if
 
    call c_to_digit ( chrtmp, ndig )
 
    if ( ihave == 3 ) then
      itop = 10 * itop + ndig
    else if ( ihave == 5 ) then
      itop = 10 * itop + ndig
      ibot = 10*ibot
    else if ( ihave == 8 ) then
      jtop = 10 * jtop + ndig
    end if
!
!  Anything else is regarded as a terminator.
!
  else
    iterm = 1
  end if
 
  if ( iterm /= 1 .and. lchar+1 < nchar ) then
    go to 10
  end if

  if ( iterm /= 1 .and. lchar+1 == nchar ) then
    lchar = nchar
  end if
!
!  Number seems to have terminated.  Have we got a legal number?
!
  if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then
    ierror = ihave
    write ( *, * ) ' '
    write ( *, * ) 'CHRCTF - Serious error!'
    write ( *, * ) '  Illegal or nonnumeric input:'
    write ( *, '(1x,a)' ) string
    return
  end if
!
!  Number seems OK.  Form it.
!
  if ( jsgn == 1 ) then
    itop = itop * 10 ** jtop
  else
    ibot = ibot * 10 ** jtop
  end if
 
  itop = isgn * itop
 
  return
end
subroutine chrctg ( string, itop, ibot, ierror, lchar )

!*****************************************************************************80
!
!! chrctg() reads an integer, decimal fraction or a ratio from a string.
!
!  Discussion:
!
!    CHRCTG returns an equivalent ratio (ITOP/IBOT).
!
!    If the input number is an integer, ITOP equals that integer, and
!    IBOT is 1.   But in the case of 2.25, the program would return
!    ITOP = 225, IBOT = 100.
!
!    A ratio is either
!      a number
!    or
!      a number, "/", a number.
!
!    A "number" is defined as:
!
!      blanks,
!      initial sign,
!      integer part,
!      decimal point,
!      fraction part,
!      E,
!      exponent sign,
!      exponent integer part,
!      blanks,
!      final comma or semicolon,
!
!    with most quantities optional.
!
!  Examples:
!
!    STRING            ITOP      IBOT
!
!    '1'               1         1
!    '     1   '       1         1
!    '1A'              1         1
!    '12,34,56'        12        1
!    '  34 7'          34        1
!    '-1E2ABCD'        -100      1
!    '-1X2ABCD'        -1        1
!    ' 2E-1'           2         10
!    '23.45'           2345      100
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRING, the string containing the
!    data to be read.  Reading will begin at position 1 and
!    terminate when no more characters
!    can be read to form a legal integer.  Blanks, commas,
!    or other nonnumeric data will, in particular, cause
!    the conversion to halt.
!
!    Output, integer ITOP, the integer read from the string,
!    assuming that no negative exponents or fractional parts
!    were used.  Otherwise, the 'integer' is ITOP/IBOT.
!
!    Output, integer IBOT, the integer divisor required to
!    represent numbers which are in decimal format or have a
!    negative exponent.
!
!    Output, integer IERROR, error flag.
!    0 if no errors,
!    Value of IHAVE in CHRCTF when error occurred otherwise.
!
!    Output, integer LCHAR, the number of characters read from
!    STRING to form the number.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer i_gcd
  integer ibot
  integer ibotb
  integer ierror
  integer itemp
  integer itop
  integer itopb
  integer lchar
  integer lchar2
  integer nchar
  character ( len = * ) string

  itop = 0
  ibot = 1
  lchar = 0
 
  call chrctf ( string, itop, ibot, ierror, lchar )

  if ( ierror /= 0) then
    return
  end if
!
!  The number is represented as a fraction.
!  If the next nonblank character is "/", then read another number.
!
  nchar = len_trim ( string )
 
  do i = lchar+1, nchar-1
 
    if ( string(i:i) == '/' ) then
 
      call chrctf ( string(i+1:), itopb, ibotb, ierror, lchar2 )

      if ( ierror /= 0 ) then
        return
      end if
 
      itop = itop * ibotb
      ibot = ibot * itopb
 
      itemp = i_gcd ( itop, ibot )
 
      itop = itop / itemp
      ibot = ibot / itemp
 
      lchar = i + lchar2
 
      return
 
    else if ( string(i:i) /= ' ' ) then
 
      return
 
    end if
 
  end do
 
  return
end
subroutine chrinp ( ierror, iounit, line, prompt )

!*****************************************************************************80
!
!! chrinp() requests new input if the LINE buffer is empty.
!
!  Discussion:
!
!    CHRINP checks to see whether there is any more information in
!    the buffer array LINE.  If so, it simply updates the prompt
!    and returns.  Otherwise, it prints the prompt string out,
!    reads the input from the user, and reprints the prompt and
!    the user input on those I/O units where it is appropriate.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IERROR, error flag.
!    0 if no errors were detected,
!    1 if there was an error in the read,
!    2 if there was an end-of-file in the read.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Input/output, character ( len = 255 ) LINE.
!
!    On input, LINE may contain information that the calling
!    program can use, or LINE may be empty.
!
!    On output, LINE is unchanged if it contained information
!    on input.  But if the input LINE was empty, then the
!    output LINE contains whatever information the user typed.
!
!    Input/output, character ( len = 255 ) PROMPT.
!    On input, the prompt string to be printed.
!    On output, PROMPT has been blanked out, up to the first comma.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer icomma
  integer ierror
  integer iosave
  integer iounit(4)
  integer lchar
  character ( len = 255 ) line
  character ( len = 255 ) output
  integer output_line_count
  character ( len = 255 ) prompt

  ierror = 0
 
10    continue
!
!  If there is nothing in the LINE buffer, then:
!    "turn off" the automatic echo for units between 30 and 39,
!    print the prompt line,
!    "turn on" the automatic echo for units between 30 and 39,
!    read the input line,
!    remove double blanks,
!    don't print a copy of the input on units between 40 and 49.
!
  if ( line == ' ' ) then
 
    do i = 2, 4
      if ( 30 <= iounit(i) .and. iounit(i) <= 39 ) then
        iounit(i) = - iounit(i)
      end if
    end do
 
    lchar = len_trim ( prompt )
    if ( 0 < lchar ) then
      output = 'Enter ' // prompt(1:lchar)
      call s_write ( iounit, output )
    end if
 
    do i = 2, 4
      if ( iounit(i) <= -30 .and. -39 <= iounit(i) ) then
        iounit(i) = - iounit(i)
      end if
    end do
 
    if ( iounit(1) <= 0 ) then
      read ( *, '(a)', end = 50, err = 40 ) line
    else
      read ( iounit(1), '(a)', end = 50, err = 40 ) line
    end if

    call s_blanks_delete ( line )
!
!  Don't echo input to IOUNIT(2).
!
    if ( iounit(1) < 40 .or. 49 < iounit(1) ) then
      iosave = iounit(2)
      if ( iounit(1) <= 0 ) then
        iounit(2) = - 1
      end if
      output = line
      call s_write ( iounit, output )
      iounit(2) = iosave
    end if
 
  end if
!
!  If the user typed something in, reset the line position to 0.
!
  if ( iounit(1) == 0 ) then
    output_line_count = 0
    call i_data ( 'SET', 'OUTPUT_LINE_COUNT', output_line_count )
  end if
!
!  If item was read, remove item from PROMPT list.
!
  if ( line /= ' ' ) then

    icomma = index ( prompt, ',' )

    if ( 0 < icomma .and. icomma < 80 .and. &
         prompt(icomma+1:icomma+1) == ' ' ) then

      icomma = icomma + 1

    end if

    call s_chop ( prompt, 1, icomma )

  end if
 
  return
!
!  Error in input.
!
40    continue

  ierror = 1
  output = ' '
  call s_write ( iounit, output )
  output = 'CHRINP - Error! '
  call s_write ( iounit, output )
  output = '  Error in input line format:'
  call s_write ( iounit, output )
  output = line
  call s_write ( iounit, output )
 
  if ( iounit(1) <= 0 ) then
    line = ' '
    go to 10
  end if
 
  return
!
!  End of input.
!
!  If we are reading from a file, then set IERROR=2 and return.
!  But if we are reading from the user, something is seriously
!  wrong, and we must stop.
!
50    continue

  ierror = 2
  line = ' '
  output = ' '
  call s_write ( iounit, output )
  output = 'CHRINP - Warning!'
  call s_write ( iounit, output )
  output = '  End of input!'
  call s_write ( iounit, output )
 
  if ( iounit(1) == 0 ) then
    output = 'The program is being stopped now!'
    call s_write ( iounit, output )
  else
    close ( unit = iounit(1) )
    iounit(1) = 0
    output = 'Closing current input file!'
    call s_write ( iounit, output )
  end if
 
  return
end
subroutine col_add ( a, iatop, iabot, ierror, iform, iounit, icol1, icol2, &
  maxcol, maxrow, nrow, sval, istop, isbot )

!*****************************************************************************80
!
!! col_add() adds a multiple of one column to another.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = rk ) A(MAXROW,MAXCOL), the current matrix.
!
!    Input/output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
!    IATOP and IABOT represent the current rational or decimal matrix.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred.
!
!    Input, integer IFORM, 0/1/2, rational/real/decimal arithmetic.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Input, integer ICOL1, the column which is to be modified.
!
!    Input, integer ICOL2, the column which is to be multiplied by
!    a given value and added to row ICOL1.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
!    Input, integer NCOL, the number of columns in the matrix.
!
!    Input, real ( kind = rk ) SVAL, the multiplier to use if real 
!    arithmetic is employed.
!
!    Input, integer ISTOP, ISBOT, the fractional or decimal
!    multiplier to use.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxcol
  integer maxrow

  real ( kind = rk ) a(maxrow,maxcol)
  character ( len = 10 ) chrtmp1
  character ( len = 10 ) chrtmp2
  character ( len = 22 ) chrtmp3
  integer i
  integer iabot(maxrow,maxcol)
  integer iatop(maxrow,maxcol)
  integer ibot
  integer icol1
  integer icol2
  integer ierror
  integer iform
  integer iounit(4)
  integer isbot
  integer isbot2
  integer istop
  integer istop2
  integer itop
  integer len1
  integer len2
  integer len3
  integer nrow
  character ( len = 255 ) output
  real ( kind = rk ) sval
!
!  Return immediately if the multiplier is zero.
!
  if ( ( iform == 0 .and. istop == 0 ) .or. &
       ( iform == 1 .and. sval == 0.0D+00 ) .or. &
       ( iform == 2 .and. istop == 0 ) ) then
      return
    end if
!
!  Carry out the operation.
!
  if ( iform == 0 ) then
 
    do i = 1, nrow
 
      call rat_mul ( isbot2, isbot, iabot(i,icol2), istop2, istop, &
        iatop(i,icol2), ierror )
 
      if ( ierror /= 0 ) then
        return
      end if

      call rat_add ( ibot, iabot(i,icol1), isbot2, itop, iatop(i,icol1), &
        istop2 )
 
      iatop(i,icol1) = itop
      iabot(i,icol1) = ibot
 
    end do
 
  else if ( iform == 1 ) then
 
    do i = 1, nrow
      a(i,icol1) = a(i,icol1) + sval * a(i,icol2)
    end do
 
  else if ( iform == 2 ) then
 
    do i = 1, nrow
 
      call dec_mul ( isbot2, isbot, iabot(i,icol2), istop2, istop, &
        iatop(i,icol2) )
 
      call dec_add ( ibot, iabot(i,icol1), isbot2, itop, iatop(i,icol1), &
        istop2 )
 
      iatop(i,icol1) = itop
      iabot(i,icol1) = ibot
 
    end do
 
  end if
!
!  Print out a message.
!
  if ( iform == 0 ) then
 
    if ( istop == isbot ) then
      chrtmp3 = '+'
    else if ( istop == - isbot ) then
      chrtmp3 = '-'
    else
      call rat_to_s_left ( istop, isbot, chrtmp3 )
    end if
 
  else if ( iform == 1 ) then
 
    if ( sval == 1.0D+00 ) then
      chrtmp3 = '+'
    else if ( sval == - 1.0D+00 ) then
      chrtmp3 = '-'
    else
      call r_to_s_left ( sval, chrtmp3 )
    end if
 
  else if ( iform == 2 ) then

    if ( istop == 1 .and. isbot == 0 ) then 
      chrtmp3 = '+'
    else if ( istop == - 1 .and. isbot == 0 ) then
      chrtmp3 = '-'
    else
      call dec_to_s_left ( istop, isbot, chrtmp3 )
    end if
 
  end if
 
  call i_to_s_left ( icol1, chrtmp1 )
  len1 = len_trim ( chrtmp1 )

  call i_to_s_left ( icol2, chrtmp2 )
  len2 = len_trim ( chrtmp2 )

  len3 = len_trim ( chrtmp3 )

  if ( chrtmp3 == '-' .or. chrtmp3 == '+' ) then

    output = '  ECO: Col ' // chrtmp1(1:len1) // ' <=  ' // 'Col ' // &
      chrtmp1(1:len1) // ' ' // chrtmp3(1:1) // ' Col ' // chrtmp2(1:len2)

  else if ( chrtmp3(1:1) == '-' .or. chrtmp3(1:1) == '+' ) then

    output = '  ECO: Col ' // chrtmp1(1:len1) // ' <=  ' // 'Col ' // &
      chrtmp1(1:len1) // ' ' // chrtmp3(1:1) // ' ' // chrtmp3(2:len3) // &
      ' Col ' // chrtmp2(1:len2)

  else

    output = '  ECO: Col ' // chrtmp1(1:len1) // ' <=  ' // 'Col ' // &
      chrtmp1(1:len1) // ' + ' // chrtmp3(1:len3) // ' Col ' // chrtmp2(1:len2)

  end if

  call s_write ( iounit, output )
 
  return
end
subroutine col_add_param ( ierror, iform, iounit, icol1, icol2, istop, isbot, &
  line, ncol, sval )

!*****************************************************************************80
!
!! col_add_param() gets and checks the column add parameters.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred.
!
!    Input, integer IFORM, 0/1/2, rational/real/decimal arithmetic.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Input, integer ICOL1, the column to which the multiple is to be added.
!
!    Input, integer ICOL2, the column which is to be multiplied and
!    added to another row.
!
!    Output, integer ISTOP, ISBOT, the parts of the rational
!    or decimal fraction of the multiplier, if that is the
!    arithmetic being used.
!
!    Workspace, character ( len = 255 ) LINE, used to hold the user's input.
!
!    Input, integer NCOL, the number of columns in the matrix.
!
!    Output, real ( kind = rk ) SVAL, the multiplier, if real 
!    arithmetic is used.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer icol1
  integer icol2
  integer ierror
  integer iform
  integer iounit(4)
  integer isbot
  integer istop
  character ( len = 255 ) line
  integer ncol
  character ( len = 255 ) output
  character ( len = 255 ) prompt
  real ( kind = rk ) sval

  prompt = 'multiplier S, column I to add, target column J.'
!
!  Get the multiplier, SVAL or ISTOP/ISBOT.
!
  if ( iform == 0 ) then
 
    call rat_read ( istop, isbot, line, prompt, iounit, ierror )
    if ( ierror /= 0 ) then
      return
    end if
 
  else if ( iform == 1 ) then
 
    call r_read ( sval, line, prompt, iounit, ierror )
    if ( ierror /= 0 ) then
      return
    end if
 
  else if ( iform == 2 ) then
 
    call dec_read ( istop, isbot, line, prompt, iounit, ierror )

    if ( ierror /= 0 ) then
      output = ' '
      call s_write ( iounit, output )
      output = 'CHECK_ADD(): Error!'
      call s_write ( iounit, output )
      output = '  DEC_READ returns error code.'
      call s_write ( iounit, output )
      return
    end if
 
    call dec_round ( istop, isbot )
 
  end if
!
!  Get the column to add.
!
  call i_read ( icol2, line, prompt, iounit, ierror )

  if ( ierror /= 0 ) then
    ierror = 1
    output = 'Error reading column index from line:'
    call s_write ( iounit, output )
    output = line
    call s_write ( iounit, output )
    return
  end if

  if ( icol2 < 1 .or. ncol < icol2 ) then
    ierror = 1
    output = 'Error!  Column index was not acceptable!'
    call s_write ( iounit, output )
    return
  end if
!
!  Get the column to which we are adding.
!
  call i_read ( icol1, line, prompt, iounit, ierror )
  if ( ierror /= 0 ) then
    return
  end if
 
  if ( icol1 < 1 .or. ncol < icol1 ) then
    ierror = 1
    output = 'Error!  The column index was not acceptable!'
    call s_write ( iounit, output )
    return
  end if
!
!  Make sure the columns are different.
!
  if ( icol1 == icol2 ) then
    output = 'Error!  The columns should not be the same!'
    call s_write ( iounit, output )
    ierror = 1
    return
  end if
 
  ierror = 0
 
  return
end
subroutine col_auto ( a, iatop, iabot, ierror, iform, iounit, maxcol, maxrow, &
  ncol, nrow )

!*****************************************************************************80
!
!! col_auto() automatically column reduces the current matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = rk ) A(MAXROW,MAXCOL), the current matrix.
!
!    Input/output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
!    IATOP and IABOT represent the rational or decimal matrix
!    to which elementary row operations will be applied.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred.
!
!    Input, integer IFORM, 0/1/2, rational/real/decimal arithmetic.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
!    Input, integer NCOL, the number of columns in the matrix.
!
!    Input, integer NROW, the number of rows in the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxcol
  integer maxrow

  real ( kind = rk ) a(maxrow,maxcol)
  real ( kind = rk ) amax
  real ( kind = rk ) atemp
  integer i
  integer iabot(maxrow,maxcol)
  integer iatop(maxrow,maxcol)
  integer ierror
  integer iform
  integer imax
  integer iounit(4)
  integer irow
  integer isbot
  integer istop
  integer j
  integer jcol
  integer kcol
  integer l
  integer lcol
  integer ncol
  integer nrow
  character ( len = 255 ) output
  real ( kind = rk ) sval

  ierror = 0
 
  do j = 1, ncol
 
    jcol = j
 
    do i = 1, nrow
 
      irow = i
!
!  In row IROW, seek the column between ICOL and NCOL with
!  maximum nonzero entry AMAX.
!
      imax = 0
      amax = 0.0D+00
 
      do kcol = jcol, ncol
 
        if ( iform == 0 ) then
          call rat_to_r ( atemp, iatop(irow,kcol), iabot(irow,kcol) )
        else if ( iform == 1 ) then
          atemp = a(irow,kcol)
        else if ( iform == 2 ) then
          call dec_to_r ( atemp, iatop(irow,kcol), iabot(irow,kcol) )
        end if
 
        atemp = abs ( atemp )
 
        if ( amax < atemp ) then
          amax = atemp
          imax = kcol
        end if
 
      end do
 
      if ( imax /= 0 ) then
        kcol = imax
        go to 10
      end if
 
    end do
 
    return
 
10      continue
 
    output = ' '
    call s_write ( iounit, output )
!
!  Interchange the JCOL-th and the pivot columns.
!
    if ( kcol /= jcol ) then
      call col_swap ( a, iatop, iabot, ierror, iform, iounit, kcol, jcol, &
        maxcol, maxrow, ncol, nrow )
    end if
!
!  Divide the pivot column by A(IROW,JCOL) so that A(IROW,JCOL) = 1.
!
    if ( iform == 0 ) then
 
      istop = iatop(irow,jcol)
      isbot = iabot(irow,jcol)
 
    else if ( iform == 1 ) then
 
      sval = a(irow,jcol)
 
    else if ( iform == 2 ) then
 
      istop = iatop(irow,jcol)
      isbot = iabot(irow,jcol)
 
    end if
 
    call col_div ( a, iatop, iabot, ierror, iform, iounit, jcol, &
      maxcol, maxrow, ncol, nrow, sval, istop, isbot )
!
!  Annihilate A(IROW,L) for L not equal to JCOL.
!
    do l = 1, nrow
 
      lcol = l
 
      if ( lcol /= jcol ) then
 
        if ( iform == 0 ) then

          if ( iatop(irow,lcol) /= 0 ) then

            istop = - iatop(irow,lcol)
            isbot = iabot(irow,lcol)

            call col_add ( a, iatop, iabot, ierror, iform, iounit, &
              lcol, jcol, maxcol, maxrow, ncol, sval, istop, isbot )

            iatop(irow,lcol) = 0
            iabot(irow,lcol) = 1

          end if

        else if ( iform == 1 ) then

          if ( a(irow,lcol) /= 0.0D+00 ) then

            sval = - a(irow,lcol)

            call col_add ( a, iatop, iabot, ierror, iform, iounit, &
              lcol, jcol, maxcol, maxrow, ncol, sval, istop, isbot )

            a(irow,lcol) = 0.0D+00

          end if

        else if ( iform == 2 ) then

          if ( iatop(irow,lcol) /= 0 ) then

            istop = - iatop(irow,lcol)
            isbot = iabot(irow,lcol)

            call col_add ( a, iatop, iabot, ierror, iform, iounit, &
              lcol, jcol, maxcol, maxrow, ncol, sval, istop, isbot )

            iatop(irow,lcol) = 0
            iabot(irow,lcol) = 0

          end if

        end if
 
      end if
 
    end do
 
  end do
 
  return
end
subroutine col_del ( a, iabot, iatop, icol, maxcol, maxrow, ncol, nrow )

!*****************************************************************************80
!
!! col_del() deletes a column by shifting other columns to the left.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = rk ) A(MAXROW,MAXCOL), the current matrix.
!
!    Output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
!    IATOP and IABOT represent the rational or decimal matrix
!    to be changed.
!
!    Input, integer ICOL, the column to be deleted.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
!    Input, integer NCOL, the number of columns in the matrix.
!
!    Input, integer NROW, the number of rows in the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxcol
  integer maxrow

  real ( kind = rk ) a(maxrow,maxcol)
  integer i
  integer iabot(maxrow,maxcol)
  integer iatop(maxrow,maxcol)
  integer icol
  integer j
  integer ncol
  integer nrow

  do j = icol, ncol - 1
    do i = 1, nrow
 
      a(i,j) = a(i,j+1)
      iatop(i,j) = iatop(i,j+1)
      iabot(i,j) = iabot(i,j+1)
 
    end do
  end do
 
  return
end
subroutine col_div ( a, iatop, iabot, ierror, iform, iounit, icol, maxcol, &
  maxrow, ncol, nrow, sval, istop, isbot )

!*****************************************************************************80
!
!! col_div() divides a column of the matrix by a scale factor.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = rk ) A(MAXROW,MAXCOL), the current matrix.
!
!    Input/output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
!    IATOP and IABOT represent the current rational or decimal
!    matrix.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred.
!
!    Input, integer IFORM, 0/1/2, rational/real/decimal arithmetic.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Input, integer ICOL, the column to be divided.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
!    Input, integer NCOL, the number of columns in the matrix.
!
!    Input, integer NROW, the number of rows in the matrix.
!
!    Input, real ( kind = rk ) SVAL, the real divisor.
!
!    Input, integer ISTOP, ISBOT, the fractional or decimal divisor.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxcol
  integer maxrow

  real ( kind = rk ) a(maxrow,maxcol)
  character ( len = 10 ) chrtmp1
  character ( len = 24 )  chrtmp2
  integer i
  integer iabot(maxrow,maxcol)
  integer iatop(maxrow,maxcol)
  integer ibot
  integer icol
  integer ierror
  integer iform
  integer iounit(4)
  integer isbot
  integer istop
  integer itop
  integer len1
  integer len2
  integer ncol
  integer nrow
  character ( len = 3 ) op
  character ( len = 255 ) output
  real ( kind = rk ) sval
!
!  Make sure that the column number is legal.
!
  if ( icol < 1 .or. ncol < icol ) then
    output = 'Error!  The column number is out of range!'
    call s_write ( iounit, output )
    ierror = 1
    return
  end if
!
!  Check for an illegal divisor of 0, or a pointless divisor of 1.
!
  if ( iform == 0 ) then
    if ( istop == 0 ) then
      output = 'Error!  It is illegal to divide by 0!'
      call s_write ( iounit, output )
      ierror = 1
      return
    else if ( istop == isbot ) then
      return
    end if
  else if ( iform == 1 ) then
    if ( sval == 0.0D+00 ) then
      output = 'Error!  It is illegal to divide by 0!'
      call s_write ( iounit, output )
      ierror = 1
      return
    else if ( sval == 1.0D+00 ) then
      return
    end if
  else if ( iform == 2 ) then
    if ( istop == 0 ) then
      output = 'Error!  It is illegal to divide by 0!'
      call s_write ( iounit, output )
      ierror = 1
      return
    else if ( istop == 1 .and. isbot==0 ) then
      return
    end if
  end if
!
!  Carry out the division.
!
  if ( iform == 0 ) then
 
    do i = 1, nrow
 
      call rat_mul ( ibot, iabot(i,icol), istop, itop, iatop(i,icol), isbot, &
        ierror )
 
      if ( ierror /= 0 ) then
        return
      end if
 
      iatop(i,icol) = itop
      iabot(i,icol) = ibot
 
    end do
 
  else if ( iform == 1 ) then
 
    do i = 1, nrow
      a(i,icol) = a(i,icol) / sval
    end do
 
  else if ( iform == 2 ) then
 
    do i = 1, nrow
 
      call dec_div ( iatop(i,icol), iabot(i,icol), istop, isbot, &
        iatop(i,icol), iabot(i,icol), ierror )
 
    end do
 
  end if
!
!  Print out a statement about what has been done.
!
  if ( iform == 0 ) then
 
    if ( isbot == 1 ) then

      call i_to_s_left ( istop, chrtmp2 )
      op = ' / '

    else

      call rat_to_s_left ( isbot, istop, chrtmp2 )
      op = ' * '

    end if
 
  else if ( iform == 1 ) then
 
    call r_to_s_left ( sval, chrtmp2 )
    op = ' / '

  else if ( iform == 2 ) then
 
    call dec_to_s_left ( istop, isbot, chrtmp2 )
    op = ' / '

  end if

  call i_to_s_left ( icol, chrtmp1 )
  len1 = len_trim ( chrtmp1 )

  len2 = len_trim ( chrtmp2 )

  output = '  ECO: Col ' // chrtmp1(1:len1) // ' <=  Col ' // chrtmp1(1:len1) &
    // op // chrtmp2(1:len2)

  call s_write ( iounit, output )
 
  return
end
subroutine col_div_param ( icol, ierror, iform, iounit, isbot, istop, line, &
  sval )

!*****************************************************************************80
!
!! col_div_param() gets and checks the column divide parameters.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred.
!
!    Input, integer IFORM, 0/1/2, rational/real/decimal arithmetic.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Workspace, character ( len = 255 ) LINE, used to hold the user's input.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer icol
  integer ierror
  integer iform
  integer iounit(4)
  integer isbot
  integer istop
  character ( len = 255 ) line
  character ( len = 255 ) output
  character ( len = 255 ) prompt
  real ( kind = rk ) sval

  prompt = 'column I, divisor S.'
!
!  Read the row number to be divided.
!
  call i_read ( icol, line, prompt, iounit, ierror )
  if ( ierror /= 0 ) then
    return
  end if
!
!  Read the divisor, either RVAL or ISTOP/ISBOT.
!
  if ( iform == 0 ) then
 
    call rat_read ( istop, isbot, line, prompt, iounit, ierror )
 
  else if ( iform == 1 ) then
 
    call r_read ( sval, line, prompt, iounit, ierror )
 
  else if ( iform == 2 ) then
 
    call dec_read ( istop, isbot, line, prompt, iounit, ierror )

    if ( ierror /= 0 ) then
      output = ' '
      call s_write ( iounit, output )
      output = 'DIVIDE - Fatal error!'
      call s_write ( iounit, output )
      output = '  DEC_READ returned error flag.'
      call s_write ( iounit, output )
      return
    end if

    call dec_round ( istop, isbot )
 
  end if
 
  if ( ierror /= 0 ) then
    return
  end if

  return
end
subroutine col_mul ( a, iatop, iabot, ierror, iform, iounit, icol, maxcol, &
  maxrow, ncol, nrow, sval, istop, isbot )

!*****************************************************************************80
!
!! col_mul() multiplies a column of the matrix by a scale factor.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = rk ) A(MAXROW,MAXCOL), the current matrix.
!
!    Input/output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
!    IATOP and IABOT represent the current rational or decimal matrix.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred.
!
!    Input, integer IFORM, 0/1/2, rational/real/decimal arithmetic.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Input, integer ICOL, the column that is to be multiplied.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
!    Input, integer NCOL, the number of columns in the matrix.
!
!    Input, integer NROW, the number of rows in the matrix.
!
!    Input, real ( kind = rk ) SVAL, the real row multiplier.
!
!    Input, integer ISTOP, ISBOT, the decimal or fractional row multiplier.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxcol
  integer maxrow

  real ( kind = rk ) a(maxrow,maxcol)
  character ( len = 10 ) chrtmp1
  character ( len = 22 ) chrtmp
  integer i
  integer iabot(maxrow,maxcol)
  integer iatop(maxrow,maxcol)
  integer ibot
  integer icol
  integer ierror
  integer iform
  integer iounit(4)
  integer isbot
  integer istop
  integer itop
  integer len1
  integer len2
  integer ncol
  integer nrow
  character ( len = 255 ) output
  real ( kind = rk ) sval
!
!  Make sure column number is OK.
!
  if ( icol < 1 .or. ncol < icol ) then
    output = 'Error!  The column number is out of range!'
    call s_write ( iounit, output )
    ierror = 1
    return
  end if
!
!  For rational arithmetic, make sure bottom of scale factor
!  is not 0.
!
  if ( iform == 0 ) then
    if ( isbot == 0 ) then
      output = 'Error!  Illegal 0 divisor in multiplier!'
      call s_write ( iounit, output )
      ierror = 1
      return
    end if
  end if
!
!  Check for multiplication by 0.
!
  if ( ( iform == 0 .and. istop == 0 ) .or. &
       ( iform == 1 .and. sval == 0.0D+00 ) .or. &
       ( iform == 2 .and. istop == 0 ) ) then
    output = ' '
    call s_write ( iounit, output )
    output = 'Warning - Multiplication by zero is not an ERO.'
    call s_write ( iounit, output )
    ierror = 1
    return
  end if
!
!  Check for multiplication by 1.
!
  if ( iform == 0 ) then
    if ( istop == isbot ) then
      return
    end if
  else if ( iform == 1 ) then
    if ( sval == 1.0D+00 ) then
      return
    end if
  else if ( iform == 2 ) then
    if ( istop == 1 .and. isbot == 0 ) then
      return
    end if
  end if
!
!  Carry out the multiplication.
!
  if ( iform == 0 ) then
 
    do i = 1, nrow
 
      call rat_mul ( ibot, iabot(i,icol), isbot, itop, iatop(i,icol), istop, &
        ierror )

      if ( ierror /= 0 ) then
        return
      end if

      iatop(i,icol) = itop
      iabot(i,icol) = ibot
 
    end do
 
  else if ( iform == 1 ) then
 
    do i = 1, nrow
      a(i,icol) = sval * a(i,icol)
    end do
 
  else if ( iform == 2 ) then
 
    do i = 1, nrow
 
      call dec_mul ( ibot, iabot(i,icol), isbot, itop, iatop(i,icol), istop )

      if ( ierror /= 0 ) then
        return
      end if
 
      iatop(i,icol) = itop
      iabot(i,icol) = ibot
 
    end do
 
  end if
!
!  Confirm the operation.
!
  if ( iform == 0 ) then
 
    if ( istop == - isbot ) then
      chrtmp = '-'
    else
      call rat_to_s_left ( istop, isbot, chrtmp )
    end if
 
  else if ( iform == 1 ) then
 
    if ( sval == - 1.0D+00 ) then
      chrtmp = '-'
    else
      call r_to_s_left ( sval, chrtmp )
    end if
 
  else if ( iform == 2 ) then

    if ( istop == - 1 .and. isbot == 0 ) then
      chrtmp = '-'
    else
      call dec_to_s_left ( istop, isbot, chrtmp )
    end if

  end if
 
  call i_to_s_left ( icol, chrtmp1 )
  len1 = len_trim ( chrtmp1 )

  len2 = len_trim ( chrtmp )
  output = '  ECO: Col ' // chrtmp1(1:len1) // ' <= ' // chrtmp(1:len2) // &
    ' Col ' // chrtmp1(1:len1)
  call s_write ( iounit, output )

  return
end
subroutine col_mul_param ( ierror, iform, iounit, icol, istop, isbot, line, &
  rval )

!*****************************************************************************80
!
!! col_mul_param() gets and checks the column multiply parameters.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred.
!
!    Input, integer IFORM, 0/1/2, rational/real/decimal arithmetic.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Output, integer ICOL, the column to be multiplied.
!
!    Output, integer ISTOP, ISBOT, the multiplier to use for
!    fractional or decimal arithmetic.
!
!    Workspace, character ( len = 255 ) LINE, used to hold the user's input.
!
!    Output, real ( kind = rk ) RVAL, the multiplier to use for real arithmetic.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer icol
  integer ierror
  integer iform
  integer iounit(4)
  integer isbot
  integer istop
  character ( len = 255 ) line
  character ( len = 255 ) prompt
  real ( kind = rk ) rval

  prompt = 'column I, multiplier S.'
!
!  Read the column number to be multiplied.
!
  call i_read ( icol, line, prompt, iounit, ierror )
  if ( ierror /= 0 ) then
    return
  end if
!
!  Read the multiplier, either RVAL or ISTOP/ISBOT.
!
  if ( iform == 0 ) then
 
    call rat_read ( istop, isbot, line, prompt, iounit, ierror )
 
  else if ( iform == 1 ) then
 
    call r_read ( rval, line, prompt, iounit, ierror )
 
  else if ( iform == 2 ) then
 
    call dec_read ( istop, isbot, line, prompt, iounit, ierror )
 
    call dec_round ( istop, isbot )
 
  end if
 
  return
end
subroutine col_op_check ( command, ierror, iounit, line2 )

!*****************************************************************************80
!
!! col_op_check() checks for commands given in the form of ECO's.
!
!  Discussion:
!
!    The form of the elementary column operation commands includes:
!
!    The column interchange command:
!      CI1 <=> CI2
!    This will fail if user types "C I1 <=> C I2"
!
!    The scalar multiply command:
!      CI1 <= S * CI1
!    with or without the "*".
!
!    The scalar divide command:
!      CI1 <= CI1 / S
!
!    The add row command:
!      CI1 <= CI1 + S *CI2
!    or
!      CI1 <= S * CI2 + CI1
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = 4 ) COMMAND.
!    If the routine decides that the user has input an ERO in the
!    natural format, then COMMAND contains the necessary
!    one letter MATMAN command to carry out the ERO.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Workspace, character ( len = 255 ) LINE2, a copy of the user input in LINE.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = 22 ) chrtmp
  character ( len = 10 ) chrtmp1
  character ( len = 10 ) chrtmp2
  character ( len = 20 ) command
  integer icol1
  integer icol2
  integer icol3
  integer idbot2
  integer idbot3
  integer idtop2
  integer idtop3
  integer ierror
  integer iounit(4)
  integer isbot2
  integer isbot3
  integer istop2
  integer istop3
  integer lchar
  logical ldiv
  character ( len = 255 ) line2
  character ( len = 255 ) output
  character ( len = 255 ) string

  command = ' '
!
!  1. Remove all blanks from the line, and capitalize it.
!
  call s_blank_delete ( line2 )
  call s_cap ( line2 )
!
!  2. Is the first character an "C" or "COL" or "COLUMN"?
!
  if ( line2(1:1) /= 'C' ) then
    return
  end if
 
  if ( line2(1:6) == 'COLUMN' ) then
    call s_chop ( line2, 1, 6 )
  else if ( line2(1:3) == 'COL' ) then
    call s_chop ( line2, 1, 3 )
  else
    call s_chop ( line2, 1, 1 )
  end if
!
!  3. The next item should be a column number, ICOL1.
!
  call s_to_i ( line2, icol1, ierror, lchar )
 
  if ( ierror /= 0 ) then
    output = 'Your ECO command could not be understood.'
    call s_write ( iounit, output )
    output = 'The first column number "C1" did not make sense.'
    call s_write ( iounit, output )
    return
  end if
 
  call s_chop ( line2, 1, lchar )
!
!  4. Check for the column interchange string "=", "<>", "<=>" or "<->".
!
  if ( line2(1:2) == '<>' ) then
    string = '<>'
  else if ( line2(1:3) == '<=>' ) then
    string = '<=>'
  else if ( line2(1:3) == '<->' ) then
    string = '<->'
  else if ( line2(1:2) == '<=' ) then
    string = '<='
  else if ( line2(1:2) == '<-' ) then
    string = '<-'
  else if ( line2(1:2) == '=>' ) then
    string = '=>'
  else if ( line2(1:2) == '->' ) then
    string = '->'
  else if ( line2(1:1) == '=' ) then
    string = '='
  else if ( line2(1:2) == ':=' ) then
    string = ':='
  else
    ierror = 1
    output = 'Your ECO command could not be understood.'
    call s_write ( iounit, output )
    output = 'The assignment symbol <=> was missing.'
    call s_write ( iounit, output )
    return
  end if
 
  lchar = len_trim ( string )
 
  call s_chop ( line2, 1, lchar )
!
!  5. The next quantity could be an explicit signed scalar, S2,
!     or an implicit +-1.
!
  if ( line2(1:1) == 'C' ) then

    istop2 = 1.0D+00
    isbot2 = 1.0D+00

  else

    if ( line2(1:2) == '+C' ) then
      istop2 = 1.0D+00
      isbot2 = 1.0D+00
      call s_chop ( line2, 1, 1 )
    else if ( line2(1:2) == '-C' ) then
      istop2 = - 1.0D+00
      isbot2 = 1.0D+00
      call s_chop ( line2, 1, 1 )
    else
      call chrctg ( line2, istop2, isbot2, ierror, lchar )
      call s_chop ( line2, 1, lchar )
 
      if ( ierror /= 0 ) then
        output = 'Your ECO command could not be understood.'
        call s_write ( iounit, output )
        output = 'The multiplier S2 did not make sense.'
        call s_write ( iounit, output )
        ierror = 1
        return
      end if
 
    end if
  end if
!
!  6. Is the next character an optional "*"?
!
  if ( line2(1:1) == '*' ) then
    call s_chop ( line2, 1, 1 )
  end if
!
!  7. Is the next character a "C"?
!
  if ( line2(1:6) == 'COLUMN' ) then
    call s_chop ( line2, 1, 6 )
  else if ( line2(1:3) == 'COL' ) then
    call s_chop ( line2, 1, 3 )
  else if ( line2(1:1) == 'C' ) then
    call s_chop ( line2, 1, 1 )
  else
    ierror = 1
    output = 'Your ECO command could not be understood.'
    call s_write ( iounit, output )
    output = 'Could not find the second column index.'
    call s_write ( iounit, output )
    return
  end if
!
!  8. The next item should be a column number, ICOL2.
!
  call s_to_i ( line2, icol2, ierror, lchar )
 
  if ( ierror /= 0 ) then
    output = 'Your ECO command could not be understood.'
    call s_write ( iounit, output )
    output = 'The second column number "C2" did not make sense.'
    call s_write ( iounit, output )
    ierror = 1
    return
  else
    call s_chop ( line2, 1, lchar )
  end if
!
!  9. If there's nothing more, this must be an interchange
!     or a scaling.  Form the equivalent MATMAN command.
!
  if ( line2 == ' ' ) then
 
    if ( icol1 == icol2 ) then

      command = 'COL_MUL'

      if ( isbot2 < 0 ) then
        isbot2 = - isbot2
        istop2 = - istop2
      end if

      call rat_to_s_left ( istop2, isbot2, chrtmp )
      call i_to_s_left ( icol1, chrtmp1 )
      line2 = chrtmp1 // ' ' // chrtmp
      call s_blanks_delete ( line2 )
      return
    end if
 
    if ( istop2 == 1 .and. isbot2 == 1 ) then
      command = 'COL_SWAP'
      call i_to_s_left ( icol1, chrtmp1 )
      call i_to_s_left ( icol2, chrtmp2 )
      line2 = chrtmp1 // ' ' // chrtmp2
      call s_blanks_delete ( line2 )
      return
    end if
 
    ierror = 1
    output = 'Your ECO command could not be understood.'
    call s_write ( iounit, output )
    output = 'A MULTIPLY command must have C1 and C2 the same.'
    call s_write ( iounit, output )
    output = 'An INTERCHANGE command cannot have a multiplier.'
    call s_write ( iounit, output )
    return
  end if
!
!  10. Is the next quantity a '/', or perhaps a '*'?
!
  ldiv = .false.
 
  if ( line2(1:1) == '/' ) then
 
    ldiv = .true.
    call s_chop ( line2, 1, 1 )
    call chrctg ( line2, idtop2, idbot2, ierror, lchar )
 
    if ( ierror /= 0 ) then
      output = 'Your ECO command could not be understood.'
      call s_write ( iounit, output )
      output = 'The divisor of column 2 did not make sense.'
      call s_write ( iounit, output )
      return
    end if
 
    istop2 = istop2 * idbot2
    isbot2 = isbot2 * idtop2
 
    if ( icol1 == icol2 ) then

      if ( ldiv ) then
        command = 'COL_DIV'
        call i_swap ( istop2, isbot2 )
      else
        command = 'COL_MUL'
      end if

      if ( isbot2 < 0 ) then
        isbot2 = - isbot2
        istop2 = - istop2
      end if

      call i_to_s_left ( icol1, chrtmp1 )
      call rat_to_s_left ( istop2, isbot2, chrtmp )
      line2 = chrtmp1 // ' ' // chrtmp
      call s_blanks_delete ( line2 )
      return

    end if
 
  else if ( line2(1:1) == '*' ) then
 
    call s_chop ( line2, 1, 1 )
    call chrctg ( line2, idtop2, idbot2, ierror, lchar )
 
    if ( ierror /= 0 ) then
      output = 'Your ECO command could not be understood.'
      call s_write ( iounit, output )
      output = 'The multiplier of column 2 did not make sense.'
      call s_write ( iounit, output )
      return
    end if
 
    istop2 = istop2 * idtop2
    isbot2 = isbot2 * idbot2
 
    if ( icol1 == icol2 ) then

      command = 'COL_MUL'

      if ( isbot2 < 0 ) then
        isbot2 = - isbot2
        istop2 = - istop2
      end if

      call i_to_s_left ( icol1, chrtmp1 )
      call rat_to_s_left ( istop2, isbot2, chrtmp )
      line2 = chrtmp1 // ' ' // chrtmp
      call s_blanks_delete ( line2 )
      return

    end if
 
  end if
!
!  11. Is the next quantity a scalar, S3?
!
  if ( line2(1:2) == '+C' ) then
 
    istop3 = 1.0D+00
    isbot3 = 1.0D+00
    call s_chop ( line2, 1, 1) 
 
  else if ( line2(1:2) == '-C' ) then
 
    istop3 = - 1.0D+00
    isbot3 = 1.0D+00
    call s_chop ( line2, 1, 1 )
 
  else
 
    call chrctg ( line2, istop3, isbot3, ierror, lchar )
 
    if ( ierror /= 0 ) then
      output = 'Your ECO command could not be understood.'
      call s_write ( iounit, output )
      output = 'The multiplier S2 did not make sense.'
      call s_write ( iounit, output )
      ierror = 1
      return
    end if
 
    call s_chop ( line2, 1, lchar )
 
  end if
!
!  12. Is the next quantity an optional "*"?
!
  if ( line2(1:1) == '*' ) then
    call s_chop ( line2, 1, 1 )
  end if
!
!  13. Is the next quantity a "C"?
!
  if ( line2(1:6) == 'COLUMN' ) then
    call s_chop ( line2, 1, 6 )
  else if ( line2(1:3) == 'COL' ) then
    call s_chop ( line2, 1, 3 )
  else if ( line2(1:1) == 'C' ) then
    call s_chop ( line2, 1, 1) 
  else
    ierror = 1
    output = 'Your ECO command could not be understood.'
    call s_write ( iounit, output )
    output = 'The "C" marking the third column was misplaced.'
    call s_write ( iounit, output )
    return
  end if
!
!  14. The next item should be a column number, ICOL3.
!
  call s_to_i ( line2, icol3, ierror, lchar )
 
  if ( ierror /= 0 ) then
    output = 'Your ECO command could not be understood.'
    call s_write ( iounit, output )
    output = 'The third column number "C3" did not make sense.'
    call s_write ( iounit, output )
    ierror = 1
    return
  end if
 
  call s_chop ( line2, 1, lchar )
!
!  15. Is the next quantity a '/', or perhaps a '*'?
!
  if ( line2(1:1) == '/' ) then
 
    call s_chop ( line2, 1, 1 )
    call chrctg ( line2, idtop3, idbot3, ierror, lchar )
 
    if ( ierror /= 0 ) then
      output = 'Your ECO command could not be understood.'
      call s_write ( iounit, output )
      output = 'The divisor of column 3 did not make sense.'
      call s_write ( iounit, output )
      return
    end if
 
    istop3 = istop3 * idbot3
    isbot3 = isbot3 * idtop3
 
  else if ( line2(1:1) == '*' ) then
 
    call s_chop ( line2, 1, 1 )
    call chrctg ( line2, idtop3, idbot3, ierror, lchar )
 
    if ( ierror /= 0 ) then
      output = 'Your ECO command could not be understood.'
      call s_write ( iounit, output )
      output = 'The multiplier of column 3 did not make sense.'
      call s_write ( iounit, output )
      return
    end if
 
    istop3 = istop3 * idtop3
    isbot3 = isbot3 * idbot3
 
  end if
!
!  16. Form the equivalent MATMAN command.
!
  if ( icol1 == icol2 ) then

    command = 'COL_ADD'

    if ( isbot3 < 0 ) then
      isbot3 = - isbot3
      istop3 = - istop3
    end if

    call i_to_s_left ( icol3, chrtmp1 )
    call i_to_s_left ( icol1, chrtmp2 )
    call rat_to_s_left ( istop3, isbot3, chrtmp )
    line2 =  chrtmp // ' ' // chrtmp1 // ' ' // chrtmp2
    call s_blanks_delete ( line2 )

  else if ( icol1 == icol3 ) then

    command = 'COL_ADD'

    if ( isbot2 < 0 ) then
      isbot2 = - isbot2
      istop2 = - istop2
    end if

    call rat_to_s_left ( istop2, isbot2, chrtmp )
    call i_to_s_left ( icol2, chrtmp1 )
    call i_to_s_left ( icol1, chrtmp2 )
    line2 = chrtmp // ' ' // chrtmp1 // ' ' // chrtmp2
    call s_blanks_delete ( line2 )

  else

    ierror = 1
    output = 'Your ECO command could not be understood.'
    call s_write ( iounit, output )
    output = 'C2 or C3 must equal C1 in an ECO command.'
    call s_write ( iounit, output )
  end if
 
  return
end
subroutine col_shift ( a, iabot, iatop, icol, maxcol, maxrow, ncol, nrow )

!*****************************************************************************80
!
!! col_shift() allows a new column to be inserted by shifting others right.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = rk ) A(MAXROW,MAXCOL), the current matrix.
!
!    Input/output, integer IABOT(MAXROW,MAXCOL), IATOP(MAXROW,MAXCOL).
!    IATOP and IABOT represent the current rational or decimal
!    matrix.
!
!    Input, integer ICOL, the position of the new column.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
!    Input, integer NCOL, the number of columns in the matrix.
!
!    Input, integer NROW, the number of rows in the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxcol
  integer maxrow

  real ( kind = rk ) a(maxrow,maxcol)
  integer i
  integer iabot(maxrow,maxcol)
  integer iatop(maxrow,maxcol)
  integer icol
  integer j
  integer ncol
  integer nrow

  do j = ncol, icol + 1, -1
    do i = 1, nrow
 
      a(i,j) = a(i,j-1)
      iatop(i,j) = iatop(i,j-1)
      iabot(i,j) = iabot(i,j-1)
 
    end do
  end do
 
  return
end
subroutine col_swap ( a, iatop, iabot, ierror, iform, iounit, icol1, icol2, &
  maxcol, maxrow, ncol, nrow )

!*****************************************************************************80
!
!! col_swap() swaps two columns of a matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = rk ) A(MAXROW,MAXCOL), the current matrix.
!
!    Input/output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
!    IATOP and IABOT represent the current rational or decimal matrix.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred.
!
!    Input, integer IFORM, 0/1/2, rational/real/decimal arithmetic.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Input, integer ICOL1, ICOL2, the rows to swap.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
!    Input, integer NCOL, the number of columns in the matrix.
!
!    Input, integer NROW, the number of rows in the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxcol
  integer maxrow

  real ( kind = rk ) a(maxrow,maxcol)
  character ( len = 10 ) chrtmp1
  character ( len = 10 ) chrtmp2
  integer iabot(maxrow,maxcol)
  integer iatop(maxrow,maxcol)
  integer i
  integer icol1
  integer icol2
  integer ierror
  integer iform
  integer iounit(4)
  integer len1
  integer len2
  integer ncol
  integer nrow
  character ( len = 255 ) output
!
!  Skip out if the two columns are the same.
!
  if ( icol1 == icol2 ) then
    output = 'You have asked to swap a column with itself!'
    call s_write ( iounit, output )
    return
  end if
!
!  Refuse to continue if a row is out of bounds.
!
  if ( ( icol1 < 1 .or. ncol < icol1 ) .or. &
       ( icol2 < 1 .or. ncol < icol2 ) ) then
    ierror = 1
    output = 'One of the columns is illegal!'
    call s_write ( iounit, output )
    return
  end if
!
!  Swap the columns.
!
  do i = 1, nrow
 
    if ( iform == 0 ) then
 
      call i_swap ( iatop(i,icol1), iatop(i,icol2) )
      call i_swap ( iabot(i,icol1), iabot(i,icol2) )
 
    else if ( iform == 1 ) then
 
      call r_swap ( a(i,icol1), a(i,icol2) )
 
    else if ( iform == 2 ) then
 
      call i_swap ( iatop(i,icol1), iatop(i,icol2) )
      call i_swap ( iabot(i,icol1), iabot(i,icol2) )
 
    end if
 
  end do
 
  call i_to_s_left ( icol1, chrtmp1 )
  len1 = len_trim ( chrtmp1 )

  call i_to_s_left ( icol2, chrtmp2 )
  len2 = len_trim ( chrtmp2 )

  output = '  ECO: Col ' // chrtmp1(1:len1) // ' <=> Col ' // chrtmp2(1:len2)

  call s_write ( iounit, output )
 
  return
end
subroutine dec_add ( ibot, ibot1, ibot2, itop, itop1, itop2 )

!*****************************************************************************80
!
!! dec_add() adds two decimal quantities.
!
!  Discussion:
!
!    The routine computes
!  
!      ITOP * 10^IBOT = ITOP1 * 10^IBOT1 + ITOP2 * 10^IBOT2
!
!    while trying to avoid integer overflow.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IBOT, the exponent of the result.
!
!    Input, integer IBOT1, IBOT2, the exponents of the numbers to be added.
!
!    Output, integer ITOP, the coefficient of the result.
!
!    Input, integer ITOP1, ITOP2, the coefficients of the numbers to be added.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ibot
  integer ibot1
  integer ibot2
  integer itop
  integer itop1
  integer itop2
  integer jtop1
  integer jtop2

  if ( itop1 == 0 ) then
    itop = itop2
    ibot = ibot2
    return
  else if ( itop2 == 0 ) then
    itop = itop1
    ibot = ibot1
    return
  else if ( ibot1 == ibot2 ) then
    itop = itop1 + itop2
    ibot = ibot1
    call dec_round ( itop, ibot )
    return
  end if
!
!  Line up the exponents.
!
  jtop1 = itop1
  jtop2 = itop2
 
  if ( ibot1 < ibot2 ) then
    jtop2 = jtop2 * 10 ** ( ibot2 - ibot1 )
  else
    jtop1 = jtop1 * 10 ** ( ibot1 - ibot2 )
  end if
!
!  Add the coefficients.
!
  itop = jtop1 + jtop2
  ibot = min ( ibot1, ibot2 )
!
!  Clean up the result.
!
  call dec_round ( itop, ibot )
 
  return
end
subroutine dec_digit_set ( ierror, iounit, line )

!*****************************************************************************80
!
!! dec_digit_set() allows the user to specify the number of decimal digits.
!
!  Discussion:
!
!    DEC_DIGIT is
!
!    * the number of digits used when converting a real number
!      to a fraction using the "FI" or "FD" command;
!
!    * the maximum number of digits in a decimal.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Workspace, character ( len = 255 ) LINE, used to hold the user's input.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = 10 ) chrtmp1
  integer dec_digit_max
  integer dec_digit
  integer ierror
  integer iounit(4)
  integer itemp
  character ( len = 255 ) line
  character ( len = 255 ) output
  character ( len = 255 ) prompt

  dec_digit_max = 0
  call i_data ( 'GET', 'DEC_DIGIT_MAX', dec_digit_max )

  output = 'How many decimal places should be used in '
  call s_write ( iounit, output )
  output = 'Converting real results to a decimal?'
  call s_write ( iounit, output )
  output = ' '
  call s_write ( iounit, output )
  output = ' 1 means 123.45 becomes 1 * 10^2'
  call s_write ( iounit, output )
  output = ' 2 means 123.45 becomes 12 * 10^1'
  call s_write ( iounit, output )
  output = ' 3 means 123.45 becomes 123'
  call s_write ( iounit, output )
  output = 'and so on.'
  call s_write ( iounit, output )
 
  call i_to_s_left ( dec_digit_max, chrtmp1 )
  prompt = 'number of decimals (1 to ' // chrtmp1 // ').'
  call s_blanks_delete ( prompt )
 
  call i_read ( itemp, line, prompt, iounit, ierror )
 
  if ( ierror /= 0 ) then
    output = 'Your choice was not acceptable!'
    call s_write ( iounit, output )
    return
  end if
!
!  Absolutely do not let DEC_DIGIT be less than 1.
!
  if ( itemp < 1 ) then
    output = 'The number of decimals must be positive!'
    call s_write ( iounit, output )
    ierror = 1
    return
  end if
!
!  Allow user to exceed the maximum, with a warning.
!
  if ( dec_digit_max < itemp ) then
    output = ' '
    call s_write ( iounit, output )
    output = 'Warning!'
    call s_write ( iounit, output )
    output = 'Your choice is larger than the recommended maximum!'
    call s_write ( iounit, output )
    output = 'which is ' // chrtmp1
    call s_blanks_delete ( output )
    call s_write ( iounit, output )
    output = 'It is possible that calculations will break down'
    call s_write ( iounit, output )
    output = 'at any time!  Be careful!'
    call s_write ( iounit, output )
  end if
 
  output = ' '
  call s_write ( iounit, output )
  dec_digit = itemp
  call i_data ( 'SET', 'DEC_DIGIT', dec_digit )

  call i_to_s_left ( dec_digit, chrtmp1 )
  output = 'The number of decimal digits will now be ' // chrtmp1
  call s_blanks_delete ( output )
  call s_write ( iounit, output )
 
  return
end
subroutine dec_div ( itop1, ibot1, itop2, ibot2, itop, ibot, ierror )

!*****************************************************************************80
!
!! dec_div() divides two decimal values.
!
!  Discussion:
!
!    A decimal quantity is stored as
!
!      (ITOP,IBOT) 
!
!    representing the value
!
!      ITOP * 10 ^ IBOT.
!
!    The routine computes 
!
!      ITOP * 10 ^ IBOT = (ITOP1 * 10^IBOT1) / (ITOP2 * 10^IBOT2)
!
!                       = (ITOP1/ITOP2) * 10^(IBOT1-IBOT2)
!
!    while avoiding integer overflow.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ITOP1, IBOT1, the numerator.
!
!    Input, integer ITOP2, IBOT2, the denominator.
!
!    Output, integer ITOP, IBOT, the result.
!
!    Output, integer IERROR.
!    0, no error occurred.
!    1, an error occurred.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) dval
  integer ibot
  integer ibot1
  integer ibot2
  integer ibot3
  integer ierror
  integer itop
  integer itop1
  integer itop2
  integer itop3
!
!  First special case, top fraction is 0.
!
  if ( itop1 == 0 ) then
    itop = 0
    ibot = 0
    return
  end if
!
!  First error, bottom of fraction is 0.
!
  if ( itop2 == 0 ) then
    ierror = 1
    itop = 0
    ibot = 0
    return
  end if
!
!  Second special case, result is 1.
!
  if ( itop1 == itop2 .and. ibot1 == ibot2 ) then
    itop = 1
    ibot = 0
    return
  end if
!
!  Third special case, result is power of 10.
!
  if ( itop1 == itop2 ) then
    itop = 1
    ibot = ibot1 - ibot2
    return
  end if
!
!  Fourth special case: ITOP1/ITOP2 is exact.
!
  if ( ( itop1 / itop2 ) * itop2 == itop1 ) then
    itop = itop1 / itop2
    ibot = ibot1 - ibot2
    return
  end if
!
!  General case.
!
  dval = real ( itop1, kind = rk ) / real ( itop2, kind = rk )
 
  call r_to_dec ( dval, itop3, ibot3 )
 
  itop = itop3
  ibot = ibot3 + ibot1 - ibot2
 
  return
end
subroutine dec_mul ( ibot, ibot1, ibot2, itop, itop1, itop2 )

!*****************************************************************************80
!
!! dec_mul() multiplies two decimals.
!
!  Discussion:
!
!    The routine computes
!
!      ITOP * 10^IBOT = (ITOP1 * 10^IBOT1) * (ITOP2 * 10^IBOT2)
!                     = (ITOP1*ITOP2) * 10^(IBOT1+IBOT2)
!
!    while avoiding integer overflow.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IBOT, the exponent of the result.
!
!    Input, integer IBOT1, the exponent of the first factor.
!
!    Input, integer IBOT2, the exponent of the second factor.
!
!    Output, integer ITOP, the coefficient of the result.
!
!    Input, integer ITOP1, the coefficient of the first factor.
!
!    Input, integer ITOP2, the coefficient of the second factor.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i_big
  real ( kind = rk ) dval
  integer ibot
  integer ibot1
  integer ibot2
  integer ibot3
  integer itop
  integer itop1
  integer itop2
  integer itop3
  real ( kind = rk ) rmax
  real ( kind = rk ) temp

  i_big = 0
  call i_data ( 'GET', 'I_BIG', i_big )
  rmax = real ( i_big, kind = rk )
!
!  The result is zero if either ITOP1 or ITOP2 is zero.
!
  if ( itop1 == 0 .or. itop2 == 0 ) then
    itop = 0
    ibot = 0
    return
  end if
!
!  The result is simple if either ITOP1 or ITOP2 is one.
!
  if ( itop1 == 1 .or. itop2 == 1 ) then
    itop = itop1 * itop2
    ibot = ibot1 + ibot2
    return
  end if
 
  temp = log ( real ( abs ( itop1 ), kind = rk ) ) &
       + log ( real ( abs ( itop2 ), kind = rk ) )
 
  if ( temp < log ( rmax ) ) then
 
    itop = itop1 * itop2
    ibot = ibot1 + ibot2
 
  else
 
    dval = real ( itop1, kind = rk ) * real ( itop2, kind = rk )
 
    call r_to_dec ( dval, itop3, ibot3 )
 
    itop = itop3
    ibot = ibot3 + ( ibot1 + ibot2 )
 
  end if
!
!  Clean up the result.
!
  call dec_round ( itop, ibot )
 
  return
end
subroutine dec_print ( iatop, iabot, ibase, iounit, ihi, ilo, jhi, jlo, &
  lpmoda, maxcol, maxrow, ncol, nrow, title )

!*****************************************************************************80
!
!! dec_print() prints out decimal vectors and matrices.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
!    IATOP and IABOT represent the decimal matrix.
!
!    Input, integer IBASE(MAXROW), the basic variables.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Input, integer IHI, ILO, the last and first rows to print.
!
!    Input, integer JHI, JLO, the last and first columns to print.
!
!    Input, integer LPMODA, 0/1, linear algebra/linear programming mode.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
!    Input, integer NCOL, the number of columns in the matrix.
!
!    Input, integer NROW, the number of rows in the matrix.
!
!    Input, character ( len = * ) TITLE, a label for the object being printed.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: ncolum = 80

  integer maxcol
  integer maxrow

  character ( len = 22 ) chrtmp
  character ( len = 10 ) chrtmp1
  character ( len = 10 ) chrtmp2
  character ( len = 10 ) chrtmp3
  character ( len = 40 ) format2
  integer i
  integer iabot(maxrow,maxcol)
  integer iatop(maxrow,maxcol)
  integer ibase(maxrow)
  integer ichi
  integer iclo
  integer ihi
  integer ilo
  integer imax
  integer imin
  integer iounit(4)
  integer izhi
  integer izlo
  integer j
  integer jhi
  integer jlo
  integer jmax
  integer jmin
  integer khi
  integer klo
  integer kmax
  character ( len = 4 ) lab
  integer lenc
  integer llab
  integer lpmoda
  integer ncol
  integer npline
  integer nrow
  character ( len = 255 ) output
  character ( len = * ) title

  if ( lpmoda == 1 ) then
    llab = 4
  else
    llab = 0
  end if
!
!  Figure out how wide we must make each column.
!
  imax = 0
  jmax = 0
 
  do i = ilo, ihi
    do j = jlo, jhi
 
      call dec_to_s_left ( iatop(i,j), iabot(i,j), chrtmp )
      lenc = len_trim ( chrtmp )
      jmax = max ( jmax, lenc )
 
    end do
  end do
 
  kmax = 2 + imax + 1 + jmax
  npline = ( ncolum - llab ) / kmax
!
!  Set up the format for the heading.
!
  if ( lpmoda == 1 ) then
    call i_to_s_left ( llab, chrtmp1 )
    call i_to_s_left ( npline, chrtmp2 )
    call i_to_s_left ( kmax, chrtmp3 )
    format2 = '(' // chrtmp1 // 'x,' // chrtmp2 // 'i' // chrtmp3 // ')'
  else
    call i_to_s_left ( npline, chrtmp2 )
    call i_to_s_left ( kmax, chrtmp3 )
    format2 = '(' // chrtmp2 // 'i' // chrtmp3 // ')'
  end if
 
  call s_blank_delete ( format2 )
 
  do jmin = jlo, jhi, npline
 
    jmax = min ( jmin+npline-1, jhi )
 
    lab = '    '
!
!  Handle a column vector.
!
    if ( jlo == jhi .and. ilo /= ihi ) then
 
      output = ' '
      call s_write ( iounit, output )
 
      if ( ilo == 1 ) then
        output = title
        call s_write ( iounit, output )
        call i_to_s_left ( jlo, chrtmp1 )
        output = 'Column ' // chrtmp1 // ' transposed.'
        call s_blanks_delete ( output )
        call s_write ( iounit, output )
      end if
 
      do imin = ilo, ihi, npline
 
        imax = min ( imin+npline-1, ihi )
 
        output = ' '
        call s_write ( iounit, output )
 
        do i = imin, imax
          ilo = 4 + ( i - imin ) * kmax + 1
          ihi = 4 + ( i - imin ) * kmax + kmax
          call dec_to_s_left ( iatop(i,jlo), iabot(i,jlo), chrtmp )
          output(ilo:ihi) = adjustr ( chrtmp(1:kmax) )
        end do
 
        call s_write ( iounit, output )
 
      end do
 
      go to 90
    end if
 
    output = ' '
    call s_write ( iounit, output )
 
    if ( jmin == 1 ) then
      output = title
      call s_write ( iounit, output )
      output = ' '
      call s_write ( iounit, output )
    end if
!
!  Print heading for linear programming table.
!
    if ( lpmoda == 1 ) then
 
      write ( output, format2 ) ( j, j = jmin, jmax )
 
      if ( jmin <= ncol-1 .and. ncol-1 <= jmax ) then
        izlo = llab + ((ncol-1)-jmin) * kmax + kmax - 2
        izhi = izlo + 2
        output(izlo:izhi) = '  P'
      end if
 
      if ( jmin <= ncol .and. ncol <= jmax ) then
        iclo = llab + (ncol-jmin) * kmax + kmax - 2
        ichi = iclo+2
        output(iclo:ichi) = '  C'
      end if
 
      call s_write ( iounit, output )
 
      output = ' '
      call s_write ( iounit, output )
!
!  Print heading for linear algebra matrix.
!
    else
 
      if ( 1 < jmin  .or. jmax < ncol .or. 1 < ilo .or. ihi < nrow ) then

        call i_to_s_left ( jmin, chrtmp1 )
        call i_to_s_left ( jmax, chrtmp2 )
        output = 'Columns ' // chrtmp1 // ' to ' // chrtmp2
        call s_blanks_delete ( output )
        call s_write ( iounit, output )
        output = ' '
        call s_write ( iounit, output )
      end if
 
    end if
 
    do i = ilo, ihi
 
      if ( lpmoda == 1 ) then

        if ( i < nrow ) then
          if ( ibase(i) < 10 ) then
            write ( lab, '(a1,i1)' ) 'X', ibase(i)
          else
            write ( lab, '(a1,i2)' ) 'X', ibase(i)
          end if
        else if ( i < ihi ) then
          lab = 'Obj2'
        else
          lab = 'Obj '
        end if

        if ( maxrow == 1 ) then
          lab = '    '
        end if

      end if
 
      if ( lpmoda == 1 ) then
        output(1:4) = lab
      else
        output(1:4) = '    '
      end if
 
      do j = jmin, jmax
        klo = 4 + (j-jmin) * kmax + 1
        khi = 4 + (j-jmin) * kmax + kmax
        call dec_to_s_left ( iatop(i,j), iabot(i,j), chrtmp )
        output(klo:khi) = adjustr ( chrtmp(1:kmax) )
      end do
 
      call s_write ( iounit, output )
 
    end do
 
90      continue
 
  end do
 
  return
end
subroutine dec_read ( itop, ibot, line, prompt, iounit, ierror )

!*****************************************************************************80
!
!! dec_read() reads a decimal, rational or integer, and returns a decimal fraction.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ITOP, IBOT, represents the decimal fraction.
!
!    Workspace, character ( len = 255 ) LINE, used to hold the user's input.
!
!    Input/output, character ( len = 255 ) PROMPT, the prompt string.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Output, integer IERROR.
!    0, no error occurred.
!    1, an error occurred.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ibot
  integer ibot1
  integer ibot2
  integer ierror
  integer iounit(4)
  integer itop
  integer itop1
  integer itop2
  integer lchar
  character ( len = 255 ) line
  character ( len = 255 ) prompt

  ierror = 0
  itop = 0
  ibot = 0
 
10    continue
 
  call chrinp ( ierror, iounit, line, prompt )

  if ( ierror /= 0 ) then
    return
  end if
 
  if ( line == ' ' ) then
    go to 10
  end if
 
  call s_to_dec ( line, itop1, ibot1, lchar )

  if ( len ( line ) <= lchar ) then
    itop = itop1
    ibot = ibot1
  else if ( line(lchar+1:lchar+1) /= '/' ) then
    itop = itop1
    ibot = ibot1
  else
    call s_chop ( line, 1, lchar+1 )
    call s_to_dec ( line, itop2, ibot2, lchar )
    call dec_div ( itop1, ibot1, itop2, ibot2, itop, ibot, ierror )
  end if
 
  call s_chop ( line, 1, lchar )

  return
end
subroutine dec_round ( itop, ibot )

!*****************************************************************************80
!
!! dec_round() rounds a decimal fraction to a given number of digits.
!
!  Discussion:
!
!    The routine takes an arbitrary decimal fraction represented by
!
!      ITOP * 10 ^ IBOT
!
!    and makes sure that ITOP has no more than the allowed number of
!    decimal digits.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ITOP, IBOT, the coefficient and exponent
!    of a decimal fraction.  On return, ITOP has no more than 
!    the allowed number of decimal digits.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dec_digit
  integer ibot
  integer itop

  if ( itop == 0 ) then
    ibot = 0
    return
  end if

  dec_digit = 0
  call i_data ( 'GET', 'DEC_DIGIT', dec_digit )
  
  do while ( 10 ** dec_digit <= abs ( itop ) )
    itop = itop / 10
    ibot = ibot + 1
  end do
  
  do while ( ( itop / 10 ) * 10 == itop )
    itop = itop / 10
    ibot = ibot + 1
  end do
 
  return
end
subroutine dec_to_r ( a, itop, ibot )

!*****************************************************************************80
!
!! dec_to_r() converts a decimal ITOP * 10^IBOT to a real value.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) A, the equivalent real value.
!
!    Input, integer ITOP, IBOT, the coefficient and exponent
!    of the decimal value.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  integer ibot
  integer itop

  a = itop * 10.0D+00 ** ibot
 
  return
end
subroutine dec_to_rat ( iatop, iabot )

!*****************************************************************************80
!
!! dec_to_rat() converts a decimal to a rational representation.
!
!  Discussion:
!
!    On input, a value is represented as IATOP * 10 ^ IABOT.
!
!    On output, approximately the same value is represented as IATOP / IABOT.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer IATOP, IABOT.
!    On input, these quantities represent the value IATOP * 10 ^ IABOT.
!    On output, these quantities represent the value IATOP / IABOT.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i_gcd
  integer iabot
  integer iatop
  integer itmp

  if ( 0 <= iabot ) then
    iatop = iatop * 10 ** iabot
    iabot = 1
  else
    iabot = 10 ** (-iabot)
    itmp = i_gcd ( iatop, iabot )
    iatop = iatop / itmp
    iabot = iabot / itmp
  end if
 
  return
end
subroutine dec_to_s_left ( ival, jval, s )

!*****************************************************************************80
!
!! dec_to_s_left() returns a left-justified representation of IVAL * 10 ^ JVAL.
!
!  Examples:
!
!    IVAL     JVAL       S
!    ----     ----       ------
!       0        0       0
!      21        3       21000
!      -3        0       -3
!     147       -2       14.7
!      16       -5       0.00016
!      34       30       Inf
!     123      -21       0.0000000000000000012
!      34      -30       0.0
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer IVAL, JVAL, integers which represent the decimal.
!
!    Output, character ( len = * ) S, the representation of the value.
!    The string is 'Inf' or '0.0' if the value was too large
!    or small to represent with a fixed point format.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = 22 ) chrrep
  integer i
  integer iget1
  integer iget2
  integer iput1
  integer iput2
  integer ival
  integer jval
  integer maxdigit
  integer ndigit
  integer nleft
  character ( len = * ) s

  s = ' '

  if ( ival == 0 ) then
    s = '0'
    return
  end if

  maxdigit = len ( s )
!
!  Store a representation of IVAL in CHRREP.
!
  write ( chrrep, '(i22)' ) ival
  call s_blank_delete ( chrrep )
  ndigit = len_trim ( chrrep )
!
!  Overflow if JVAL is positive, and MAXDIGIT < NDIGIT + JVAL.
!
  if ( jval > 0 ) then
    if ( maxdigit < ndigit + jval ) then
      s = 'Inf'
      return
    end if
  end if
!
!  Underflow if JVAL is negative, and MAXDIGIT < 3 + NDIGIT - JVAL.
!
  if ( jval < 0 ) then
    if ( ival > 0 ) then
      if ( maxdigit < 3 - ndigit - jval ) then
        s = '0.0'
        return
      end if
    else
      if ( maxdigit < 5 - ndigit - jval ) then
        s = '0.0'
        return
      end if
    end if
  end if
!
!  If JVAL is nonnegative, insert trailing zeros.
!
  if ( 0 <= jval ) then

    s(1:ndigit) = chrrep(1:ndigit)

    do i = ndigit+1, ndigit+jval
      s(i:i) = '0'
    end do

  else if ( jval < 0 ) then

    iput2 = 0
    iget2 = 0
!
!  Sign.
!
    if ( ival < 0 ) then
      iput1 = 1
      iput2 = 1
      iget2 = 1
      s(iput1:iput2) = '-'
      ndigit = ndigit - 1
    end if
!
!  Digits of the integral part.
!
    if ( 0 < ndigit + jval ) then
      iput1 = iput2 + 1
      iput2 = iput1 + ndigit + jval -1
      iget1 = iget2 + 1
      iget2 = iget1 + ndigit+jval - 1
      s(iput1:iput2) = chrrep(iget1:iget2)
    else
      iput1 = iput2 + 1
      iput2 = iput1
      s(iput1:iput2) = '0'
    end if
!
!  Decimal point.
!
    iput1 = iput2 + 1
    iput2 = iput1
    s(iput1:iput2) = '.'
!
!  Leading zeroes.
!
    do i = 1, - jval - ndigit
      iput1 = iput2 + 1
      iput2 = iput1
      s(iput1:iput2) = '0'
    end do

    nleft = min ( -jval, ndigit )
    nleft = min ( nleft, maxdigit - iput2 )
    iput1 = iput2 + 1
    iput2 = iput1 + nleft - 1
    iget1 = iget2 + 1
    iget2 = iget1 + nleft - 1
    s(iput1:iput2) = chrrep(iget1:iget2)

  end if

  return
end
subroutine deladd ( a, iabot, iatop, ibase, ierror, iform, iounit, line, &
  lpmoda, maxcol, maxrow, ncol, ncon, nrow, n_slack, nvar )

!*****************************************************************************80
!
!! deladd() deletes or adds a row or column to the matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = rk ) A(MAXROW,MAXCOL), the current matrix.
!
!    Input/output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
!    IATOP and IABOT represent the rational or decimal matrix
!    to be changed.
!
!    Input/output, integer IBASE(MAXROW), the basic variables.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred.
!
!    Input, integer IFORM, 0/1/2, rational/real/decimal arithmetic.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Workspace, character ( len = 255 ) LINE, used to hold the user's input.
!
!    Input, integer LPMODA, 0/1, linear algebra/linear programming mode.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
!    Input/output, integer NCOL, the number of columns in the matrix.
!
!    Input/output, integer NCON, the number of constraints.
!
!    Input/output, integer NROW, the number of rows in the matrix.
!
!    Input/output, integer N_SLACK, the number of slack variables.
!
!    Input, integer NVAR, the number of basic variables.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxcol
  integer maxrow

  real ( kind = rk ) a(maxrow,maxcol)
  character ( len = 10 ) chrtmp
  integer iabot(maxrow,maxcol)
  integer iatop(maxrow,maxcol)
  integer ibase(maxrow)
  integer icol
  integer ierror
  integer iform
  integer iounit(4)
  integer irow
  character isay
  integer iterm
  character ( len = 255 ) line
  integer lpmoda
  integer ncol
  integer ncon
  integer nrow
  integer n_slack
  integer nvar
  character ( len = 255 ) output
  character ( len = 255 ) prompt
  logical s_eqi

  ierror = 0
 
  if ( lpmoda == 0 ) then
 
    prompt = '"+" to add a row or column, "-" to delete one.'
    iterm = 0
    call s_read ( isay, line, prompt, iounit, ierror, iterm )
 
    if ( isay == '-' ) then
 
      prompt = 'R or C to delete a row or column.'

      iterm = 0
      call s_read ( isay, line, prompt, iounit, ierror, iterm )
!
!  -R: Delete a row in linear algebra mode.
!  --  -----------------------------------
!
      if ( s_eqi ( isay, 'R' ) ) then
!
!  Get row index.
!
        call i_to_s_left ( nrow, chrtmp )
        prompt = 'row to delete, between 1 and ' // chrtmp
        call s_blanks_delete ( prompt )
        call i_read ( irow, line, prompt, iounit, ierror )
        if ( ierror /= 0 ) then
          return
        end if
 
        if ( irow < 1 .or. nrow < irow ) then
          ierror = 1
          output = 'Your row index was not acceptable!'
          return
        end if
!
!  Shift matrix rows.
!
        call row_del ( a, iabot, iatop, irow, maxcol, maxrow, ncol, nrow )
 
        nrow = nrow - 1
 
        output = 'The row has been deleted!'
        call s_write ( iounit, output )
!
!  -C:  Delete a column in linear algebra mode.
!  --   --------------------------------------
!
      else if ( s_eqi ( isay, 'C' ) ) then
!
!  Get column index.
!
        call i_to_s_left ( ncol, chrtmp )
        prompt = 'column to delete, between 1 and ' // chrtmp
        call s_blanks_delete ( prompt )
        call i_read ( icol, line, prompt, iounit, ierror )
        if ( ierror /= 0 ) then
          return
        end if
 
        if ( icol < 1 .or. ncol < icol ) then
          ierror = 1
          output = 'Your column index was not acceptable!'
          return
        end if
!
!  Shift matrix columns.
!
        call col_del ( a, iabot, iatop, icol, maxcol, maxrow, ncol, nrow )
 
        ncol = ncol - 1
 
        output = 'The column has been deleted!'
        call s_write ( iounit, output )
 
      else
 
        ierror = 1
 
      end if
 
    else if ( isay == '+' ) then
 
      prompt = 'R or C to add a row or column.'
      iterm = 0
      call s_read ( isay, line, prompt, iounit, ierror, iterm )
!
!  +R:  Add a row in linear programming mode.
!  --   ------------------------------------
!
      if ( s_eqi ( isay, 'R' ) ) then
 
        if ( nrow < maxrow ) then
 
          nrow = nrow + 1
!
!  Get row index.
!
          call i_to_s_left ( nrow, chrtmp )
          prompt = 'index for new row between 1 and ' // chrtmp
          call s_blanks_delete ( prompt )
          call i_read ( irow, line, prompt, iounit, ierror )
          if ( ierror /= 0 ) then
            return
          end if
 
          if ( irow < 1 .or. nrow < irow ) then
            ierror = 1
            output = 'Your row index was not acceptable!'
            return
          end if
!
!  Shift matrix rows.
!
          call row_shift ( a, iabot, iatop, irow, maxcol, maxrow, ncol, nrow )
!
!  Read in values for new row.
!
          call la_inp1 ( a, iabot, iatop, ierror, iform, iounit, line, maxcol, &
            maxrow, irow, irow, 1, ncol )
 
        else
          ierror = 1
          output = 'There is no space for more rows!'
          call s_write ( iounit, output )
        end if
      end if
!
!  +C: Add a column in linear programming mode.
!  --  ---------------------------------------
! 
      if ( s_eqi ( isay, 'C' ) ) then
 
        if ( ncol < maxcol ) then
 
          ncol=ncol+1
!
!  Get column index.
!
          call i_to_s_left ( ncol, chrtmp )
          prompt = 'index for new column between 1 and ' // chrtmp
          call s_blanks_delete ( prompt )
          call i_read ( icol, line, prompt, iounit, ierror )
          if ( ierror /= 0 ) then
            return
          end if
 
          if ( icol < 1 .or. ncol < icol ) then
            ierror = 1
            output = 'Your column index was not acceptable!'
            return
          end if
!
!  Shift matrix columns.
!
          call col_shift ( a, iabot, iatop, icol, maxcol, maxrow, ncol, nrow )
!
!  Read in values for new column.
!
          irow = 0
 
          call la_inp1 ( a, iabot, iatop, ierror, iform, iounit, line, maxcol, &
            maxrow, 1, nrow, icol, icol )
 
        else
          output = 'Error!  There is no space for more rows!'
          call s_write ( iounit, output )
        end if
 
      end if
    end if
!
!  Add new constraint and slack variable for linear programming.
!
  else
 
    if ( maxrow - 2 <= nrow ) then
      ierror = 1
      output = 'Error!'
      call s_write ( iounit, output )
      output = 'The table cannot be increased in size to'
      call s_write ( iounit, output )
      output = 'make room for the new constraint!'
      call s_write ( iounit, output )
      return
    end if
 
    if ( maxcol <= nvar ) then
      ierror = 1
      output = 'The table cannot be increased in size to'
      call s_write ( iounit, output )
      output = 'make room for the new slack variable!'
      call s_write ( iounit, output )
      return
    end if
 
    output = 'Add a new constraint and slack variable!'
    call s_write ( iounit, output )
!
!  Shift last row down, shift last column to right.
!
!
!  DOES THIS DEPEND ON WHETHER ARTIFICIAL VARIABLES ARE INVOLVED?
!
    ncon = ncon + 1
    irow = ncon
    nrow = nrow + 1
    call row_shift ( a, iabot, iatop, irow, maxcol, maxrow, ncol, nrow )
    n_slack = n_slack + 1
    icol = nvar + n_slack
    ncol = ncol + 1
    call col_shift ( a, iabot, iatop, icol, maxcol, maxrow, ncol, nrow )
    ibase(irow) = nvar + n_slack
!
!  Read in values of constraint.
!
  end if
 
  return
end
subroutine digit_to_c ( digit, c )

!*****************************************************************************80
!
!! digit_to_c() returns the character representation of a decimal digit.
!
!  Example:
!
!    DIGIT   C
!    -----  ---
!      0    '0'
!      1    '1'
!    ...    ...
!      9    '9'
!     17    '*'
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIGIT, the digit value between 0 and 9.
!
!    Output, character C, the corresponding character, or '*' if DIGIT
!    was illegal.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character c
  integer digit

  if ( 0 <= digit .and. digit <= 9 ) then

    c = char ( digit + 48 )

  else

    c = '*'

  end if

  return
end
subroutine file_append ( file_name, ierror, inew, iold, iounit, nrec )

!*****************************************************************************80
!
!! file_append() allows us to append information to a pre-existing file.
!
!  Discussion:
!
!    This routine was created to address the fact that ANSI FORTRAN
!    does not let one easily append information to a sequential
!    access file once it has been closed.  In order to allow a user
!    to append new information, we create a new, writeable copy
!    of the file by means of a temporary copy.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the old file.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred.
!
!    Input, integer INEW, the unit number on which the new copy
!    should be opened.
!
!    Input, integer IOLD, the unit number on which the old file
!    should be opened.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Output, integer NREC, the number of records in the old file.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = 10 ) chrtmp
  character ( len = * ) file_name
  character ( len = 255 ) filtmp
  integer ierror
  integer inew
  integer iold
  integer ios
  integer iounit(4)
  character ( len = 255 ) line
  integer nrec
  character ( len = 255 ) output

  filtmp = 'tmpfil.dat'
  ierror = 0
  nrec = 0
!
!  Open old file as readable.  If it doesn't exist, we can
!  skip ahead.  Otherwise, also open new file as writeable.
!
  open ( unit = iold, file = file_name, status = 'old', err = 50 )

  rewind iold

  open (  unit = inew, file = filtmp,  status = 'new', err = 60 )
!
!  Copy old into temporary, then delete old.
!
  do

    read ( iold, '(a80)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    nrec = nrec + 1
    write ( inew, '(a80)' ) line

  end do
 
  call i_to_s_left ( nrec, chrtmp )
  output = 'The file contains ' // chrtmp // ' lines.'
  call s_blanks_delete ( output )
  call s_write ( iounit, output )
  close ( unit = iold, status = 'delete' )
  close ( unit = inew )
!
!  Reopen old as writeable, write copy of temporary into it.
!
  open (  unit = iold,  file = file_name, status = 'new', err = 60 )

  open ( unit = inew, file = filtmp, status = 'old', err = 60 )

  rewind inew
 
  do

    read ( inew, '(a80)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    write ( iold, '(a80)' ) line

  end do
!
!  Delete temporary file and return.
!
  close ( unit = inew, status = 'delete' )
  return
!
!  The file does not exist.  We may write into it immediately.
!
50    continue

  output = 'Creating a new file.'
  call s_write ( iounit, output )

  open ( unit = iold, file = file_name, status = 'new', err = 60 )

  return
!
!  Delete old copy of FILTMP.
!
60    continue

  ierror = 1
  output = 'The file could not be opened!'
  call s_write ( iounit, output )
 
  return
end
subroutine file_read ( filex, ierror, iform, iounit, line, lpmoda )

!*****************************************************************************80
!
!! file_read() reads an example from a file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) FILEX.
!    On input, the default name of the example file.
!    On output, the chosen name of the example file.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred.
!
!    Input, integer IFORM, 0/1/2, rational/real/decimal arithmetic.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Workspace, character ( len = 255 ) LINE, used to hold the user's input.
!
!    Output, integer LPMODA, 0/1, linear algebra/linear programming mode.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = 10 ) chrtmp
  character ( len = * ) filex
  character ( len = 255 ) file_name
  integer ierror
  integer iform
  integer ilabel
  integer ios
  integer iounit(4)
  integer iterm
  integer iunit
  integer jform
  integer jlabel
  integer jpmoda
  character ( len = 255 ) line
  integer lpmoda
  character ( len = 255 ) output
  character ( len = 255 ) prompt
  logical s_eqi
  character ( len = 255 ) ylabel

  ierror = 0 
  prompt = 'filename to read, default= "' // trim ( filex ) // '".'
  call s_blanks_delete ( prompt )

  iterm = 0
  call s_read ( file_name, line, prompt, iounit, ierror, iterm )

  if ( ierror /= 0 ) then
    return
  end if
 
  if ( file_name(1:1) /= ' ' ) then
    filex = file_name
  end if
 
  call get_unit ( iunit )

  open ( unit = iunit, file = filex, status = 'old', iostat = ios )
 
  if ( ios /= 0 ) then
    ierror = 1
    output = 'Error!  The example file could not be opened!'
    call s_write ( iounit, output )
    return
  end if

  iounit(1) = iunit
!
!  Read and print labels.
!
  ilabel = 0
  ylabel = 'to cancel.'
  call i_to_s_left ( ilabel, chrtmp )
  output = trim ( chrtmp ) // ' ' // trim ( ylabel )
  call s_blanks_delete ( output )
  call s_write ( iounit, output )
  prompt = ' '
 
10    continue

  line = ' '

  iterm = 0
  call s_read ( ylabel, line, prompt, iounit, ierror, iterm )

  if ( ierror == 0 ) then

    call s_blanks_delete ( ylabel )
    call s_cap ( ylabel(1:6) )
 
    if ( s_eqi ( ylabel(1:6), 'LABEL:' ) ) then
      ilabel = ilabel + 1
      ylabel(1:6) = ' '
      call i_to_s_left ( ilabel, chrtmp )
      output = trim ( chrtmp ) // ' ' // trim ( ylabel )
      call s_blanks_delete ( output )
      call s_write ( iounit, output )
    end if
 
    go to 10

  end if
!
!  Close file.
!
  ierror = 0
  close ( unit = iounit(1) )
  iounit(1) = 0
!
!  Get example number from user.
!
30    continue

  line = ' '
  prompt = 'example number.'
  call i_read ( jlabel, line, prompt, iounit, ierror )
  if ( ierror /= 0 ) then
    return
  end if
 
  if ( jlabel <= 0 .or. ilabel < jlabel ) then
    output = 'Your choice was not acceptable.'
    call s_write ( iounit, output )
    go to 30
  end if
!
!  Reopen file, seek that example number.
!
  ilabel = 0

  call get_unit ( iunit )

  open ( unit = iunit, file = filex, status = 'old', iostat = ios )
 
  if ( ios /= 0 ) then
    ierror = 1
    output = 'Error!  The example file could not be opened!'
    call s_write ( iounit, output )
    return
  end if

  iounit(1) = iunit
  prompt = ' '
 
40    continue

  line = ' '

  iterm = 0
  call s_read ( ylabel, line, prompt, iounit, ierror, iterm )

  if ( ierror /= 0 ) then
    go to 50
  end if

  call s_blanks_delete ( ylabel )
  call s_cap ( ylabel(1:6) )
 
  if ( s_eqi ( ylabel(1:6), 'LABEL:' ) ) then
    ilabel = ilabel + 1
    if ( ilabel == jlabel ) then
      go to 60
    end if
  end if
 
  go to 40
 
50    continue

  ierror = 1
  close ( unit = iounit(1) )
  iounit(1) = 0
  output = 'Could not retrieve example.'
  call s_write ( iounit, output )
  return
 
60    continue
!
!  See if the arithmetic mode should be changed.
!
  line = ' '
  call i_read ( jform, line, prompt, iounit, ierror )

  if ( ierror /= 0 ) then
    iounit(1) = 0
    return
  end if
 
  if ( iform /= jform ) then
 
    if ( jform == 0 ) then
      output = 'Arithmetic switched to rational form.'
      iform = jform
    else if ( jform == 1 ) then
      output = 'Arithmetic switched to real form.'
      iform = jform
    else if ( jform == 2 ) then
      output = 'Arithmetic switched to decimal form.'
      iform = jform
    else
      call i_to_s_left ( jform, chrtmp )
      output = 'Illegal value for arithmetic switch:' // chrtmp
      call s_write ( iounit, output )
      ierror = 1
      iounit(1) = 0
      return
    end if
 
    call s_write ( iounit, output )
 
  end if
!
!  See if linear programming mode should be changed.
!
  line = ' '
  call i_read ( jpmoda, line, prompt, iounit, ierror )
  if ( ierror /= 0 ) then
    iounit(1) = 0
    return
  end if
 
  if ( lpmoda /= jpmoda ) then
    if ( jpmoda == 0 ) then
      output = 'Switching to linear algebra mode.'
      lpmoda = jpmoda
    else if ( jpmoda == 1 ) then
      output = 'Switching to linear programming mode.'
      lpmoda = jpmoda
    else
      call i_to_s_left ( jpmoda, chrtmp )
      output = 'Illegal value for linear programming mode:' // chrtmp
      call s_write ( iounit, output )
      iounit(1) = 0
      ierror = 1
      return
    end if
    call s_write ( iounit, output )
  end if
 
  line = ' '

  return
end
subroutine file_write ( a, chineq, filex, iatop, iabot, ierror, iform, &
  iounit, line, lpmoda, maxcol, maxrow, nart, ncol, nrow, nvar )

!*****************************************************************************80
!
!! file_write() writes an example to a file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A(MAXROW,MAXCOL), the current matrix.
!
!    Input, character CHINEQ(MAXROW), the '<', '=', or '>'
!    sign for each linear programming constraint.
!
!    Input/output, character ( len = * ) FILEX.
!    On input, the default name of the example file.
!    On output, the chosen name of the example file.
!
!    Input, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
!    IATOP and IABOT represent the current rational or decimal matrix.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred.
!
!    Input, integer IFORM, 0/1/2, rational/real/decimal arithmetic.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Workspace, character ( len = 255 ) LINE, used to hold the user's input.
!
!    Input, integer LPMODA, 0/1, linear algebra/linear programming mode.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
!    Input, integer NART, the number of artificial variables.
!
!    Input, integer NCOL, the number of columns in the matrix.
!
!    Input, integer NROW, the number of rows in the matrix.
!
!    Input, integer NVAR, the number of basic variables.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxcol
  integer maxrow

  real ( kind = rk ) a(maxrow,maxcol)
  character chineq(maxrow)
  character ( len = 10 ) chrtmp1
  character ( len = 10 ) chrtmp2
  character ( len = * ) filex
  character ( len = 255 ) file_name
  integer i
  integer iabot(maxrow,maxcol)
  integer iatop(maxrow,maxcol)
  integer ierror
  integer iform
  integer ii
  integer inew
  integer iold
  integer iosave
  integer iounit(4)
  character isay
  integer iterm
  integer j
  integer jinc
  integer k
  integer khi
  character ( len = 255 ) line
  integer lpmoda
  integer nart
  integer ncol
  integer nrec
  integer nrow
  integer nvar
  character ( len = 255 ) output
  character ( len = 255 ) prompt
  character ( len = 255 ) xlabel
!
!  Get the filename to use.
!
  prompt = 'file to use, default= "' // trim ( filex ) // '".'
  call s_blanks_delete ( prompt )

  iterm = 0
  call s_read ( file_name, line, prompt, iounit, ierror, iterm )

  if ( ierror /= 0 ) then
    return
  end if

  if ( file_name(1:1) /= ' ' ) then
    filex = file_name
  end if
!
!  Get the label to use.
!
  line = ' '
  prompt = 'label.'

  iterm = 0
  call s_read ( xlabel, line, prompt, iounit, ierror, iterm )
  if ( ierror /= 0 ) then
    return
  end if
!
!  Open the file, whether it is new or old, and prepare to write
!  the new information at the end of the file.
!
  iold = 31
  inew = 32
  call file_append ( filex, ierror, inew, iold, iounit, nrec )
  if ( ierror /= 0 ) then
    return
  end if

  iounit(4) = 31
  output = 'label:    ' // trim ( xlabel )
  call s_write ( iounit, output )
  iounit(2) = -1
  iosave = iounit(3)
  iounit(3) = -1

  call i_to_s_left ( iform, chrtmp1 )
  output = chrtmp1 // ', iform (0 fraction, 1 real, 2 decimal)'
  call s_write ( iounit, output )

  call i_to_s_left ( lpmoda, chrtmp1 )
  output = chrtmp1 // ', lpmode (1 for linear programming)'
  call s_write ( iounit, output )

  if ( lpmoda == 0 ) then
    call i_to_s_left ( nrow, chrtmp1 )
    call i_to_s_left ( ncol, chrtmp2 )
    output = chrtmp1 // ',' // chrtmp2
  else
    call i_to_s_left ( nrow-1, chrtmp1 )
    call i_to_s_left ( nvar, chrtmp2 )
    output = chrtmp1 // ',' // chrtmp2
  end if

  call s_blank_delete ( output )
  call s_write ( iounit, output )
 
  if ( lpmoda == 0 ) then
 
    do i = 1, nrow
 
      isay = ' '
 
      if ( iform == 0 ) then
        jinc = 3
      else if ( iform == 1 ) then
        jinc = 5
      else if ( iform == 2 ) then
        jinc = 3
      end if
 
      do j = 1, ncol, jinc
 
        khi = min ( j+jinc-1, ncol )
 
        if ( iform == 0 ) then
          write ( output, 70 ) isay, ( iatop(i,k), iabot(i,k), k = j, khi )
        else if ( iform == 1 ) then
          write ( output, 80 ) isay, a(i,j:khi)
        else if ( iform == 2 ) then
          write ( output, 70 ) isay, ( iatop(i,k), iabot(i,k), k = j, khi )
        end if
 
        call s_blanks_delete ( output )
        call s_write ( iounit, output )
 
      end do
 
    end do
 
  else
 
    do i = 1, nrow
 
      ii = i

      if ( i == nrow .and. 0 < nart ) then
        ii = nrow + 1
      end if
 
      if ( i == nrow ) then
 
        iatop(ii,1:nvar) = - iatop(ii,1:nvar)
        a(ii,1:nvar) = - a(ii,1:nvar)
 
      end if
 
      isay = chineq(i)
 
      if ( iform == 0 ) then
        jinc = 3
      else if ( iform == 1 ) then
        jinc = 5
      else if ( iform == 2 ) then
        jinc = 3
      end if
 
      do j = 1, nvar, jinc
 
        khi = min ( j+jinc-1, nvar )
 
        if ( iform == 0 ) then
          write ( output, 70 )isay, ( iatop(ii,k), iabot(ii,k), k = j, khi )
        else if ( iform == 1 ) then
          write ( output, 80 ) isay, ( a(ii,k), k = j, khi )
        else if ( iform == 2 ) then
          write ( output, 70 ) isay, ( iatop(ii,k), iabot(ii,k), k = j, khi )
        end if
 
        isay = ' '
        call s_blanks_delete ( output )
        call s_write ( iounit, output )
 
      end do
 
      if ( iform == 0 ) then
        write ( output, 70 ) isay, iatop(ii,ncol), iabot(ii,ncol)
      else if ( iform == 1 ) then
        write ( output, 80 ) isay, a(ii,ncol)
      else if ( iform == 2 ) then
        write ( output, 70 ) isay, iatop(ii,ncol), iabot(ii,ncol)
      end if
 
      call s_blanks_delete ( output )
      call s_write ( iounit, output )
 
      if ( i == nrow ) then
        iatop(ii,1:nvar) = - iatop(ii,1:nvar)
        a(ii,1:nvar) = - a(ii,1:nvar)
      end if
 
    end do
 
  end if
!
!  Write one blank line to avoid end-of-file problems on Macintosh.
!
  output = ' '
  call s_write ( iounit, output )
 
  iounit(2) = 0
  iounit(3) = iosave
  close(unit = iounit(4))
  iounit(4) = - 1
  output = 'The problem has been stored.'
  call s_write ( iounit, output )
 
  return
70    format(a1,3(i12,'/',i12,','))
80    format(a1,5(g14.6,','))
end
subroutine form ( a, b, c, iatop, iabot, ibtop, ibbot, ictop, icbot, iform, &
  imat, iounit, jform, maxcol, maxrow )

!*****************************************************************************80
!
!! form() converts from one arithmetic form to another.
!
!  Discussion:
!
!    On input, IFORM contains a code for the current arithmetic
!    form, and JFORM contains the code for the new form.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = rk ) A(MAXROW,MAXCOL), B(MAXROW,MAXCOL),
!    C(MAXROW,MAXCOL), the current real matrix, and its two
!    backup copies.
!
!    Input/output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL),
!    IBTOP(MAXROW,MAXCOL), IBBOT(MAXROW,MAXCOL), ICTOP(MAXROW,MAXCOL),
!    ICBOT(MAXROW,MAXCOL), the current fractional or decimal matrix
!    and its two backup copies.
!
!    Input/output, integer IFORM, 0/1/2, rational/real/decimal arithmetic.
!
!    Input, integer IMAT, 0/1, a matrix HAS NOT/HAS been defined.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Input, integer JFORM, the arithmetic to be converted to.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxcol
  integer maxrow

  real ( kind = rk ) a(maxrow,maxcol)
  real ( kind = rk ) b(maxrow,maxcol)
  real ( kind = rk ) c(maxrow,maxcol)
  integer i
  integer iabot(maxrow,maxcol)
  integer iatop(maxrow,maxcol)
  integer ibbot(maxrow,maxcol)
  integer ibtop(maxrow,maxcol)
  integer icbot(maxrow,maxcol)
  integer ictop(maxrow,maxcol)
  integer ierror
  integer iform
  integer imat
  integer iounit(4)
  integer j
  integer jform
  character ( len = 255 ) output

  if ( iform == jform ) then
    output = 'You are already using the arithmetic type that'
    call s_write ( iounit, output )
    output = 'you have requested.'
    call s_write ( iounit, output )
    return
  end if
!
!  If there's no matrix, then just set the arithmetic mode
!  and return.
!
  if ( imat == 0 ) then
    iform = jform
    return
  end if
!
!  Convert the matrix data.
!  In a special case, there is data in the matrix beyond row NROW.
!  (Linear programming, with artificial variables).
!  So just convert MAXROW by MAXCOL, rather than the more modest
!  NROW by NCOL.
!
  if ( jform == 0 .and. iform == 1 ) then
 
    do i = 1, maxrow
      do j = 1, maxcol
 
        call r_to_rat ( a(i,j), iatop(i,j), iabot(i,j) )
        call r_to_rat ( b(i,j), ibtop(i,j), ibbot(i,j) )
        call r_to_rat ( c(i,j), ictop(i,j), icbot(i,j) )
 
      end do
    end do
 
  else if ( jform == 0 .and. iform == 2 ) then
 
    do i = 1, maxrow
      do j = 1, maxcol
 
        call dec_to_rat ( iatop(i,j), iabot(i,j) )
        call dec_to_rat ( ibtop(i,j), ibbot(i,j) )
        call dec_to_rat ( ictop(i,j), icbot(i,j) )
  
      end do
    end do
 
  else if ( jform == 1 .and. iform == 0 ) then
 
    do i = 1, maxrow
      do j = 1, maxcol
 
        call rat_to_r ( a(i,j), iatop(i,j), iabot(i,j) )
        call rat_to_r ( b(i,j), ibtop(i,j), ibbot(i,j) )
        call rat_to_r ( c(i,j), ictop(i,j), icbot(i,j) )
 
      end do
    end do
 
  else if ( jform == 1 .and. iform == 2 ) then
 
    do i = 1, maxrow
      do j = 1, maxcol
 
        call dec_to_r ( a(i,j), iatop(i,j), iabot(i,j) )
        call dec_to_r ( b(i,j), ibtop(i,j), ibbot(i,j) )
        call dec_to_r ( c(i,j), ictop(i,j), icbot(i,j) )
 
      end do
    end do
 
  else if ( jform == 2 .and. iform == 0 ) then
 
    do i = 1, maxrow
      do j = 1, maxcol
 
        call rat_to_dec ( iatop(i,j), iabot(i,j), ierror )
        call rat_to_dec ( ibtop(i,j), ibbot(i,j), ierror )
        call rat_to_dec ( ictop(i,j), icbot(i,j), ierror )
 
      end do
    end do
 
  else if ( jform == 2 .and. iform == 1 ) then
 
    do i = 1, maxrow
      do j = 1, maxcol
 
        call r_to_dec ( a(i,j), iatop(i,j), iabot(i,j) )
        call r_to_dec ( b(i,j), ibtop(i,j), ibbot(i,j) )
        call r_to_dec ( c(i,j), ictop(i,j), icbot(i,j) )
 
      end do
    end do
 
  end if
!
!  Update the arithmetic form.
!
  iform = jform
 
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
!    free FORTRAN unit.  GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

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
!
!  No free unit was found.
!
  iunit = 0

  return
end
subroutine hello ( iounit )

!*****************************************************************************80
!
!! hello() greets the user on startup.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 April 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer iounit(4)
  character ( len = 255 ) output

  output = ' '
  call s_write ( iounit, output )
  call timestring ( output )
  call s_write ( iounit, output )
  output = 'MATMAN, version 1.61'
  call s_write ( iounit, output )
  output = 'Last modified on 01 July 2000.'
  call s_write ( iounit, output )
  output = ' '
  call s_write ( iounit, output )
  output = 'An interactive program to perform elementary'
  call s_write ( iounit, output )
  output = 'row and column operations on a matrix, or'
  call s_write ( iounit, output )
  output = 'the simplex method of linear programming.'
  call s_write ( iounit, output )
  output = ' '
  call s_write ( iounit, output )
  output = 'Developed by Charles Cullen and John Burkardt.'
  call s_write ( iounit, output )
  output = 'All rights reserved by the authors.  This program may'
  call s_write ( iounit, output )
  output = 'not be reproduced in any form without written permission.'
  call s_write ( iounit, output )
  output = ' '
  call s_write ( iounit, output )
  output = 'Send comments to burkardt@psc.edu.'
  call s_write ( iounit, output )
 
  return
end
subroutine help ( iounit )

!*****************************************************************************80
!
!! help() prints out a brief list of the available commands.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer iounit(4)
  character ( len = 255 ) output

  output = ' '
  call s_write ( iounit, output )
  output = 'Here is a list of all MATMAN commands:'
  call s_write ( iounit, output )
  output = ' '
  call s_write ( iounit, output )
 
  output = 'A(I,J)=S  Set matrix entry to S.'
  call s_write ( iounit, output )
  output = 'B      Set up sample problem.'
  call s_write ( iounit, output )
  output = 'BASIC I, J changes basic variable I to J.'
  call s_write ( iounit, output )
  output = 'CHECK  Check matrix for reduced row echelon form.'
  call s_write ( iounit, output )
  output = 'CHECK  Check linear program table for optimality.'
  call s_write ( iounit, output )
  output = 'COL_AUTO  Automatic column reduction.'
  call s_write ( iounit, output )
  output = 'DEC    Use decimal arithmetic.'
  call s_write ( iounit, output )
  output = 'DEC_DIGIT Set the number of decimal digits.'
  call s_write ( iounit, output )
  output = 'E      Enter matrix with I rows and J columns.'
  call s_write ( iounit, output )
  output = 'E      Enter a linear program, I constraints, J variables.'
  call s_write ( iounit, output )
  output = 'G      Add/delete a row or column of the matrix.'
  call s_write ( iounit, output )
  output = 'H      for quick help.'
  call s_write ( iounit, output )
  output = 'HELP   for full help (this list).'
  call s_write ( iounit, output )
  output = 'I_BIG  Set size of largest integer for fractions.'
  call s_write ( iounit, output )
  output = 'INIT   Initialize data.'
  call s_write ( iounit, output )
  output = 'J      Jacobi rotation in (I,J) plane.'
  call s_write ( iounit, output )
  output = 'K      Open/close the transcript file.'
  call s_write ( iounit, output )
  output = 'L      To switch between linear algebra and linear programming.'
  call s_write ( iounit, output )
  output = 'P      Pivot linear program, entering I, departing J.'
  call s_write ( iounit, output )
  output = 'Q      Quit.'
  call s_write ( iounit, output )
  output = 'RES    Restore a saved matrix or table'
  call s_write ( iounit, output )
  output = 'RAT    Use rational arithmetic.'
  call s_write ( iounit, output )
  output = 'REAL   Use real arithmetic.'
  call s_write ( iounit, output )
  output = 'ROW_AUTO  Automatic row reduction.'
  call s_write ( iounit, output )
  output = 'S      Store the current matrix or table.'
  call s_write ( iounit, output )
  output = 'T      Type out the matrix'
  call s_write ( iounit, output )
  output = 'TS     Type linear programming solution.'
  call s_write ( iounit, output )
  output = 'U      Undo last operation.'
  call s_write ( iounit, output )
  output = 'V      Remove LP artificial variables.'
  call s_write ( iounit, output )
  output = 'W/X    Write/read example to/from file.'
  call s_write ( iounit, output )
  output = 'Y      Turn automatic printing ON or OFF.'
  call s_write ( iounit, output )
  output = '#      Begins a comment line.'
  call s_write ( iounit, output )
  output = '<      Get input from a file.'
  call s_write ( iounit, output )
  output = '%/$    Turn paging on/off.'
  call s_write ( iounit, output )

  output = ' '
  call s_write ( iounit, output )
  output = 'R1 <=> R2  interchanges two rows.'
  call s_write ( iounit, output )
  output = 'R1 <= S R1 multiplies a row by S.'
  call s_write ( iounit, output )
  output = 'R1 <= R1 + S R2 adds a multiple of another row.'
  call s_write ( iounit, output )
  output = ' '
  call s_write ( iounit, output )
  output = 'C1 <=> C2  interchanges two columns.'
  call s_write ( iounit, output )
  output = 'C1 <= S C1 multiplies a column by S.'
  call s_write ( iounit, output )
  output = 'C1 <= C1 + S C2 adds a multiple of another column.'
  call s_write ( iounit, output )

  return
end
subroutine get_help ( file_help, iounit, line )

!*****************************************************************************80
!
!! get_help() provides extensive help from the MATMAN help file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILHLP, the name of the help file.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Workspace, character ( len = 255 ) LINE, used to hold the user's input.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: maxtop = 40

  logical c_is_digit
  character ( len = 255 ) choice
  character ( len = 255 ) ctemp
  character ( len = 255 ) ctemp2
  character ( len = * ) file_help
  integer i
  integer ierror
  integer iline
  character ( len = 255 ) inline
  integer ios
  integer iounit(4)
  integer iterm
  integer itop
  integer iunit
  integer jerror
  character lab
  integer lchar
  integer lenc
  integer level
  character ( len = 255 ) levelc(maxtop)
  integer levelm(10)
  integer levelo
  integer levelt(maxtop)
  character ( len = * ) line
  integer move
  integer ntop
  integer num
  character ( len = 255 ) output
  integer output_line_count
  character ( len = 255 ) prompt
  logical s_eqi

  ierror = 0
  call get_unit ( iunit )

  output_line_count = 0
  call i_data ( 'SET', 'OUTPUT_LINE_COUNT', output_line_count )
!
!  Open the help file.
!
  open ( unit = iunit, file = file_help, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    output = 'Could not open the help file "' // trim ( file_help ) // '".'
    call s_write ( iounit, output )
    return
  end if

  levelo = 0
  level = 1
  iline = 1
!
!  Move to beginning of current topic by reading MOVE lines from
!  the top of the file.  Record this position, corresponding to
!  the current LEVEL, in LEVELM, in case we later want to back up.
!
!  Print out the heading line of this topic.
!
10    continue

  jerror = 0
  move = iline
  levelm(level) = iline
 
  do i = 1, move-1
    read ( iunit, '(1x)', end = 110, err = 110 )
  end do
 
  output = ' '
  call s_write ( iounit, output )
  read ( iunit, '(a1,a75)', end = 110, err = 110 ) lab, inline
  output = inline
  call s_write ( iounit, output )
!
!  If 'going down' or redisplaying, (as opposed to backing up),
!  display information available under the current topic.
!
!  We stop printing when we hit a numeric label.
!
!  If this label is less than or equal to current level, there are
!  no subtopics.
!
!  Otherwise, we now move ahead to print out the list of subtopics
!  available for this topic.
!
  if ( levelo <= level ) then

    ntop = -1
 
30      continue
 
    read ( iunit, '(a1,a75)', end = 50 ) lab, inline
    move = move + 1
 
    if ( c_is_digit ( lab ) ) then
      read ( lab, '(i1)' ) num
      if ( num <= level ) then
        go to 50
      end if
      ntop = 0
      go to 40
    end if
 
    output = inline
    call s_write ( iounit, output )
    go to 30
  else
    ntop = 0
    inline = ' '
    lab = ' '
  end if
!
!  Locate each subtopic by examining column 1, searching for
!  integer label.
!
!  Assuming we are at level LEVEL, we are searching for labels
!  equal to LEVEL+1.  As we encounter each such label, we want to
!  store the rest of the line as a subtopic.  We ignore labels
!  greater than LEVEL+1 because these are sub-subtopics, and we
!  cease our search when we reach a label less than or equal to
!  LEVEL.
!
40    continue
 
  if ( c_is_digit ( lab ) ) then

    read ( lab, '(i1)' ) num
    if ( num <= level ) then
      go to 50
    end if

    if ( num == level+1 ) then

      ntop = ntop + 1
 
      if ( ntop == 1 ) then
        output = ' '
        call s_write ( iounit, output )
        output = 'Help is available on:'
        call s_write ( iounit, output )
        output = ' '
        call s_write ( iounit, output )
      end if
 
      output = inline
      call s_write ( iounit, output )
      levelt(ntop) = move
      levelc(ntop) = inline

    end if

  end if

  read ( iunit, '(a1,a75)', end = 50, err = 50 ) lab, inline
  move = move + 1
  go to 40
 
50    continue
!
!  Display subtopics.
!
  output = ' '
  call s_write ( iounit, output )
  output = 'Return to back up, ? to redisplay.'
  call s_write ( iounit, output )
!
!  Prompt for user choice of new topic, exit, or back up.
!
60    continue
 
  ierror = 0
  line = ' '
 
  if ( 0 < ntop ) then
    prompt = 'topic you want help on, or RETURN or ?.'
  else
    prompt = 'RETURN or ?.'
  end if
 
  iterm = 0
  call s_read ( choice, line, prompt, iounit, ierror, iterm )

  if ( ierror /= 0 ) then
    ierror = 0
    close ( unit = iunit )
    return
  end if
 
  output_line_count = 0
  call i_data ( 'SET', 'OUTPUT_LINE_COUNT', output_line_count )

  call s_blanks_delete ( choice )
  lenc = len_trim ( choice )
  if ( lenc <= 0 ) then
    choice = '!'
  end if
  ctemp = choice
!
!  Two errors in a row, OK, but three suggests that something is wrong.
!
  if ( ierror /= 0 ) then
    jerror = jerror + 1
    if ( jerror <= 4) then
      go to 60
    end if
    output = 'Too many input errors in a row!'
    call s_write ( iounit, output )
  end if
!
!  Consider ending this help session.
!
  if ( ( ctemp == '!' .and. level == 1 ) .or. 4 < jerror ) then
    close ( unit = iunit )
    return
  end if
!
!  User wants to back up to a supertopic.  We must rewind.
!
  rewind iunit
  levelo = level

  if ( ctemp == '!' ) then

    level = level-1
    iline = levelm(level)
!
!  Redisplay current topic.
!
  else if ( ctemp == '?' ) then

    go to 10
!
!  User wants to go down to a subtopic.
!
  else
 
    do i = 1, ntop

      ctemp2 = levelc(i)
      call s_blanks_delete ( ctemp2 )
      itop = i

      if ( s_eqi ( ctemp(1:lenc), ctemp2(1:lenc) ) ) then
        go to 90
      end if

    end do
 
    lchar = len_trim ( choice )
    output = 'Sorry, no help available on "' // choice(1:lchar) // '".'
    call s_blanks_delete ( output )
    call s_write ( iounit, output )
    jerror = jerror + 1
    go to 60
 
90      continue

    level = level + 1
    iline = levelt(itop)

  end if
 
  go to 10
!
!  Error reading help file.
!
110   continue

  ierror = 1
  output = 'Unexpected error while reading "' // trim ( file_help ) // '".'
  call s_write ( iounit, output )
  close ( unit = iunit )
 
  return
end
subroutine i_data ( op, var, ival )

!*****************************************************************************80
!
!! i_data() stores and retrieves common data items.
!
!  Discussion:
!
!    This routine works like a sort of COMMON block.  It stores or returns
!    the values of certain variables.  Thus, it allows routines
!    to "communicate" without having to have data passed up and
!    down the calling tree in argument lists.
!
!    The variables stored by this version of the routine are:
!
!    'DEC_DIGIT', the number of digits stored for decimals;
!    'DEC_DIGIT_MAX', the maximum number of digits stored for decimals;
!    'I_BIG', the biggest integer to use in calculations.
!    'I_MAX', the maximum integer.
!    'OUTPUT_LINE_COUNT', the number of lines of output since the last pause;
!    'OUTPUT_PAGE_LENGTH', the number of lines per page.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OP, describes the operation to be done.
!    'SET' means set a value.
!    'INC' means increment a value (and return its new value)
!    'GET' means get a value.
!
!    Input, character ( len = * ) VAR, the name of the variable.
!
!    Input/output, integer IVAL.
!    If OP is 'SET', then the variable named in VAR is set to the
!    value IVAL.
!    If OP is 'GET', then the value of IVAL is set to the value of
!    the variable named in VAR.
!    If OP is 'INC', then the value of IVAL is incremented by 1,
!    and its new value is returned in VAR.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, save :: dec_digit = 7
  integer, save :: dec_digit_max = 7
  integer, save :: i_big = 2147483647
  integer, save :: i_max = 2147483647
  integer ival
  character ( len = * ) op
  integer, save :: output_line_count = 0
  integer, save :: output_page_length = 24
  logical s_eqi
  character ( len = * ) var

  if ( s_eqi ( op, 'SET' ) ) then

    if ( s_eqi ( var, 'DEC_DIGIT' ) ) then
      dec_digit = ival
    else if ( s_eqi ( var, 'DEC_DIGIT_MAX' ) ) then
      dec_digit_max = ival
    else if ( s_eqi ( var, 'I_BIG' ) ) then
      i_big = ival
    else if ( s_eqi ( var, 'I_MAX' ) ) then
      i_max = ival
    else if ( s_eqi ( var, 'OUTPUT_LINE_COUNT' ) ) then
      output_line_count = ival
    else if ( s_eqi ( var, 'OUTPUT_PAGE_LENGTH' ) ) then
      output_page_length = ival
    end if

  else if ( s_eqi ( op, 'GET' ) ) then

    if ( s_eqi ( var, 'DEC_DIGIT' ) ) then
      ival = dec_digit
    else if ( s_eqi ( var, 'DEC_DIGIT_MAX' ) ) then
      ival = dec_digit_max
    else if ( s_eqi ( var, 'I_BIG' ) ) then
      ival = i_big
    else if ( s_eqi ( var, 'I_MAX' ) ) then
      ival = i_max
    else if ( s_eqi ( var, 'OUTPUT_LINE_COUNT' ) ) then
      ival = output_line_count
    else if ( s_eqi ( var, 'OUTPUT_PAGE_LENGTH' ) ) then
      ival = output_page_length
    end if

  else if ( s_eqi ( op, 'INC' ) ) then

    if ( s_eqi ( var, 'DEC_DIGIT' ) ) then
      dec_digit = dec_digit + 1
      ival = dec_digit
    else if ( s_eqi ( var, 'DEC_DIGIT_MAX' ) ) then
      dec_digit_max = dec_digit_max + 1
      ival = dec_digit_max
    else if ( s_eqi ( var, 'I_BIG' ) ) then
      i_big = i_big + 1
      ival = i_big
    else if ( s_eqi ( var, 'I_MAX' ) ) then
      i_max = i_max + 1
      ival = i_max
    else if ( s_eqi ( var, 'OUTPUT_LINE_COUNT' ) ) then
      output_line_count = output_line_count + 1
      ival = output_line_count
    else if ( s_eqi ( var, 'OUTPUT_PAGE_LENGTH' ) ) then
      output_page_length = output_page_length + 1
      ival = output_page_length
    end if

  end if
 
  return
end
function i_gcd ( i, j )

!*****************************************************************************80
!
!! i_gcd() finds the greatest common divisor of I and J.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, J, two numbers whose greatest common divisor
!    is desired.
!
!    Output, integer I_GCD, the greatest common divisor of I and J.
!
!    only the absolute values of I and J are
!    considered, so that the result is always nonnegative.
!
!    If I or J is 0, I_GCD is returned as max ( 1, abs ( I ), abs ( J ) ).
!
!    If I and J have no common factor, I_GCD is returned as 1.
!
!    Otherwise, using the Euclidean algorithm, I_GCD is the
!    largest common factor of I and J.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer i_gcd
  integer ip
  integer iq
  integer ir
  integer j

  i_gcd = 1
!
!  Return immediately if either I or J is zero.
!
  if ( i == 0 ) then
    i_gcd = max ( 1, abs ( j ) )
    return
  else if ( j == 0 ) then
    i_gcd = max ( 1, abs ( i ) )
    return
  end if
!
!  Set IP to the larger of I and J, IQ to the smaller.
!  This way, we can alter IP and IQ as we go.
!
  ip = max ( abs ( i ), abs ( j ) )
  iq = min ( abs ( i ), abs ( j ) )
!
!  Carry out the Euclidean algorithm.
!
  do

    ir = mod ( ip, iq )

    if ( ir == 0 ) then
      exit
    end if

    ip = iq
    iq = ir

  end do

  i_gcd = iq

  return
end
subroutine i_random ( ilo, ihi, i )

!*****************************************************************************80
!
!! i_random() returns a random integer in a given range.
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
!    Input, integer ILO, IHI, the minimum and maximum acceptable values.
!
!    Output, integer I, the randomly chosen integer.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer ihi
  integer ilo
  real ( kind = rk ) r
  real ( kind = rk ) rhi
  real ( kind = rk ) rlo
!
!  Pick a random number in (0,1).
!
  call random_number ( harvest = r )
!
!  Set a real interval [RLO,RHI] which contains the integers [ILO,IHI],
!  each with a "neighborhood" of width 1.
!
  rlo = real ( ilo, kind = rk ) - 0.5D+00
  rhi = real ( ihi, kind = rk ) + 0.5D+00
!
!  Set I to the integer that is nearest the scaled value of R.
!
  i = nint ( ( 1.0D+00 - r ) * rlo + r * rhi )
!
!  In case of oddball events at the boundary, enforce the limits.
!
  i = max ( i, ilo )
  i = min ( i, ihi )

  return
end
subroutine i_read ( intval, line, prompt, iounit, ierror )

!*****************************************************************************80
!
!! i_read() reads an integer from the input buffer.
!
!  Discussion:
!
!    The routine accepts LINE which contains input and a PROMPT line.  
!    If LINE is empy, the PROMPT will be printed and LINE read from the 
!    input unit, IOUNIT(1).
! 
!    In either case, the integer INTVAL will be read from LINE,
!    beginning at character 1 and ending at the first comma, slash,
!    blank, or the end of LINE.
!
!    The PROMPT should consist of a string of names of data items,
!    separated by commas, with the current one first.
!
!    The program will print 'ENTER' PROMPT and after reading LINE
!    will strip the characters corresponding to INTVAL from LINE,
!    and the part of PROMPT up to the first comma, leaving LINE and
!    PROMPT ready for another call to I_READ, S_READ, RAT_READ or
!    R_READ.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer INTVAL, the integer that was read.
!
!    Workspace, character ( len = 255 ) LINE, used to hold the user's input.
!
!    Input, character ( len = 255 ) PROMPT, the prompt string.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ierror
  integer intval
  integer iounit(4)
  integer lchar
  character ( len = 255 ) line
  character ( len = 255 ) output
  character ( len = 255 ) prompt

  ierror = 0
  intval = 0
!
!  Retrieve a likely character string from input.
!
10    continue

  call chrinp ( ierror, iounit, line, prompt )

  if ( ierror /= 0 ) then
    output = ' '
    call s_write ( iounit, output )
    output = 'I_READ - Fatal error!'
    call s_write ( iounit, output )
    output = '  CHRINP returned error flag!'
    call s_write ( iounit, output )
    return
  end if

  if ( line == ' ' ) then
    go to 10
  end if
!
!  Convert the character string to an integer.
!
  call s_to_i ( line, intval, ierror, lchar )

  if ( ierror /= 0 ) then
    output = ' '
    call s_write ( iounit, output )
    output = 'I_READ - Fatal error!'
    call s_write ( iounit, output )
    output = '  S_TO_I returned error flag!'
    call s_write ( iounit, output )
    return
  end if
!
!  Remove the character string from the input line.
!
  if ( lchar < len ( line ) ) then
    if ( line(lchar+1:lchar+1) == ',' ) then
      lchar = lchar + 1
    end if
  end if

  call s_chop ( line, 1, lchar )
 
  return
end
subroutine i_swap ( i, j )

!*****************************************************************************80
!
!! i_swap() switches two integer values.
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

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer j
  integer k

  k = i
  i = j
  j = k

  return
end
subroutine i_to_s_left ( intval, s )

!*****************************************************************************80
!
!! i_to_s_left() converts an integer to a left_justified string.
!
!  Examples:
!
!    Assume that S is 6 characters long:
!
!    INTVAL  S
!
!         1  1
!        -1  -1
!         0  0
!      1952  1952
!    123456  123456
!   1234567  ******  <-- Not enough room!
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer INTVAL, an integer to be converted.
!
!    Output, character ( len = * ) S, the representation of the integer.
!    The integer will be left-justified.  If there is not enough space,
!    the string will be filled with stars.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character c
  integer i
  integer idig
  integer ihi
  integer ilo
  integer intval
  integer ipos
  integer ival
  character ( len = * ) s

  s = ' '

  ilo = 1
  ihi = len ( s )

  if ( ihi <= 0 ) then
    return
  end if
!
!  Make a copy of the integer.
!
  ival = intval
!
!  Handle the negative sign.
!
  if ( ival < 0 ) then

    if ( ihi <= 1 ) then
      s(1:1) = '*'
      return
    end if

    ival = - ival
    s(1:1) = '-'
    ilo = 2

  end if
!
!  The absolute value of the integer goes into S(ILO:IHI).
!
  ipos = ihi

10    continue
!
!  Find the last digit of IVAL, strip it off, and stick it into
!  the string.
!
  idig = mod ( ival, 10 )
  ival = ival / 10

  if ( ipos < ilo ) then
    do i = 1, ihi
      s(i:i) = '*'
    end do
    return
  end if

  call digit_to_c ( idig, c )

  s(ipos:ipos) = c
  ipos = ipos - 1

  if ( ival /= 0 ) then
    go to 10
  end if
!
!  Shift the string to the left.
!
  s(ilo:ilo+ihi-ipos-1) = s(ipos+1:ihi)
  s(ilo+ihi-ipos:ihi) = ' '

  return
end
subroutine infile ( filinp, ierror, iounit, line )

!*****************************************************************************80
!
!! infile() handles a new input file.
!
!  Discussion:
!
!    This routine reads the name of an input file, and changes the internal
!    values of IOUNIT(1), and opens the file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) FILINP, the input file name.
!    On input, this is a default value, or the name of a previously
!    used input file.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Workspace, character ( len = 255 ) LINE, used to hold the user's input.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = 255 ) file_name
  character ( len = * ) filinp
  integer ierror
  integer ios
  integer iounit(4)
  integer iterm
  integer iunit
  character ( len = 255 ) line
  character ( len = 255 ) output
  integer output_page_length
  character ( len = 255 ) prompt
!
!  If we were already reading an input file, close it!
!
  if ( iounit(1) /= 0 ) then
    output = 'Closing previous input file "' // trim ( filinp ) // '".'
    call s_write ( iounit, output )
    close ( unit = iounit(1) )
    iounit(1) = 0
  end if
!
!  Get the name of the input file.
!
  prompt = 'file name, default= "' // trim ( filinp ) // '".'
  iterm = 0
  call s_read ( file_name, line, prompt, iounit, ierror, iterm )
  if ( ierror /= 0 ) then
    return
  end if
!
!  If the input file is "*", then the user is typing input.
!
  if ( file_name == '*' ) then
    output = 'MATMAN now expects input directly from the user.'
    call s_write ( iounit, output )
    output_page_length = 24
    call i_data ( 'SET', 'OUTPUT_PAGE_LENGTH', output_page_length )
    output = 'Paging turned ON.'
    call s_write ( iounit, output )
    return
  end if
 
  if ( file_name /= ' ' ) then
    filinp = file_name
  end if
 
  call get_unit ( iunit )

  open ( unit = iunit, file = filinp, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    iounit(1) = 0
    output = 'MATMAN could not open the input file!'
    call s_write ( iounit, output )
    return
  end if

  iounit(1) = iunit
!
!  Turn paging off.
!
  output_page_length = 0
  call i_data ( 'SET', 'OUTPUT_PAGE_LENGTH', output_page_length )
  output = 'Paging turned OFF.'
  call s_write ( iounit, output )
 
  output = 'MATMAN now expects input from "' // trim ( filinp ) // '".'
  call s_write ( iounit, output )

  return
end
subroutine init ( a, autop, chineq, command, comold, filex, file_help, &
  filinp, file_tran, iabot, iatop, iauthr, ibase, ierror, iform, &
  imat, iounit, iseed, islbot, isltop, line, lpmoda, maxcol, maxrow, nart, &
  ncol, ncon, nrow, n_slack, nvar, sol )

!*****************************************************************************80
!
!! init() initializes the data.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) A(MAXROW,MAXCOL), the current matrix.
!
!    Output, logical AUTOP, .TRUE. if the matrix should be
!    automatically printed after most operations, .FALSE. otherwise.
!
!    Output, character CHINEQ(MAXROW), the '<', '=', or '>'
!    sign for each linear programming constraint.
!
!    Output, character ( len = 20 ) COMMAND, the newest command from the user.
!
!    Output, character ( len = 20 ) COMOLD, the previous command from the user.
!
!    Output, character ( len = * ) FILEX, the default examples file.
!
!    Output, character ( len = * ) FILHLP, the default help file.
!
!    Output, character ( len = * ) FILINP, the default input file.
!
!    Output, character ( len = * ) FILE_TRAN, the default transcript file.
!
!    Output, integer IABOT(MAXROW,MAXCOL), IATOP(MAXROW,MAXCOL).
!    IATOP and IABOT represent the current rational or decimal
!    matrix.
!
!    Output, integer IAUTHR,
!    0 if the user has typed the correct password.
!    1 if the user has not typed the correct password.
!
!    Output, integer IBASE(MAXROW), the basic variables.
!
!    Output, integer IERROR.
!    The error flag, which is initialized to zero by this routine.
!
!    Output, integer IFORM, 0/1/2, rational/real/decimal arithmetic.
!
!    Output, integer IMAT, 0/1, a matrix HAS NOT/HAS been defined.
!
!    Output, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Output, integer ISEED, a random number generator seed.
!
!    Output, integer ISLBOT(MAXROW), ISLTOP(MAXROW), the decimal
!    or fractional representation of the linear programming solution.
!
!    Output, character ( len = 255 ) LINE,
!    a buffer used to hold the user's input.
!
!    Output, integer LPMODA, 0/1, linear algebra/linear programming mode.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
!    Output, integer NART, the number of artificial variables.
!
!    Output, integer NCOL, the number of columns in the matrix.
!
!    Output, integer NCON, the number of constraints.
!
!    Output, integer NROW, the number of rows in the matrix.
!
!    Output, integer N_SLACK, the number of slack variables.
!
!    Output, integer NVAR, the number of basic variables.
!
!    Output, real ( kind = rk ) SOL(MAXROW), the current linear 
!    programming solution.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxcol
  integer maxrow

  real ( kind = rk ) a(maxrow,maxcol)
  logical autop
  character chineq(maxrow)
  character ( len = 20 ) command
  character ( len = 20 ) comold
  integer dec_digit
  integer dec_digit_max
  character ( len = * ) filex
  character ( len = * ) file_help
  character ( len = * ) filinp
  character ( len = * ) file_tran
  integer i
  integer i_big
  integer i_max
  integer iabot(maxrow,maxcol)
  integer iatop(maxrow,maxcol)
  integer iauthr
  integer ibase(maxrow)
  integer ierror
  integer iform
  integer imat
  integer iounit(4)
  integer iseed
  integer islbot(maxcol)
  integer isltop(maxcol)
  character ( len = * ) line
  integer lpmoda
  integer nart
  integer ncol
  integer ncon
  integer nrow
  integer n_slack
  integer nvar
  integer output_page_length
  real ( kind = rk ) sol(maxcol)

  call mat_zero ( a, iabot, iatop, iform, maxcol, maxrow )
 
  autop = .true.
  chineq(1:maxrow) = ' '
  command = ' '
  comold = ' '
!
!  The file names, as given here, assume that the MATMAN files
!  are in the directory where the user is working.
!
!  If MATMAN is installed on a computer in a special directory,
!  but the user wishes to run it while working in another
!  directory, then the names of FILHLP and FILE_KEY must be
!  changed to include the directory information.
!
!  Similarly, if a single copy of MATMAN is installed on a
!  multi-user computer, then the file names for FILHLP and FILE_KEY
!  would need to be changed to include the directory information.
!
  filex = 'matman.dat'
  file_help = 'matman.hlp'
  filinp = 'matman.inp'
  file_tran = 'matman.lpt'
!
!  Set IAUTHR to 0 to force students to use the authorization key.
!
  iauthr = 1
 
  do i = 1, maxrow
    ibase(i) = i
  end do

  ierror = 0
  iform = 0
  imat = 0
  iounit(1) = 0
  iounit(2) = 0
  iounit(3) = -1
  iounit(4) = -1
  iseed = 10031952
  isltop(1:maxcol) = 0
  islbot(1:maxcol) = 1
  line = ' '

  output_page_length = 24
  call i_data ( 'SET', 'OUTPUT_PAGE_LENGTH', output_page_length )
 
  lpmoda = 0

  dec_digit_max = 7
  call i_data ( 'SET', 'DEC_DIGIT_MAX', dec_digit_max )

  dec_digit = 4
  call i_data ( 'SET', 'DEC_DIGIT', dec_digit )

  i_max = 2147483647
  call i_data ( 'SET', 'I_MAX', i_max )

  i_big = i_max
  call i_data ( 'SET', 'I_BIG', i_big )

  nart = 0
  ncol = 0
  ncon = 0
  nrow = 0
  n_slack = 0
  nvar = 0
  sol(1:maxcol) = 0.0D+00
 
  return
end
subroutine la_help ( iounit )

!*****************************************************************************80
!
!! la_help() prints out a brief list of useful linear algebra commands.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer iounit(4)
  character ( len = 255 ) output

  output = 'A(I,J)=S  Set matrix entry to S.'
  call s_write ( iounit, output )
  output = 'CHECK checks if the matrix is row reduced.'
  call s_write ( iounit, output )
  output = 'E     enters a matrix to work on.'
  call s_write ( iounit, output )
  output = 'HELP  for full help.'
  call s_write ( iounit, output )
  output = 'L     switches to linear programming.'
  call s_write ( iounit, output )
  output = 'Q     quits.'
  call s_write ( iounit, output )
  output = 'Z     automatic row reduction (requires password).'
  call s_write ( iounit, output )
  output = '?     for interactive help.'
  call s_write ( iounit, output )
  output = ' '
  call s_write ( iounit, output )
  output = 'R1 <=> R2  interchanges two rows'
  call s_write ( iounit, output )
  output = 'R1 <= S R1 multiplies a row by S.'
  call s_write ( iounit, output )
  output = 'R1 <= R1 + S R2 adds a multiple of another row.'
  call s_write ( iounit, output )
  output = ' '
  call s_write ( iounit, output )
 
  return
end
subroutine la_inp0 ( ierror, iounit, line, maxcol, maxrow, ncol, nrow )

!*****************************************************************************80
!
!! la_inp0() begins the process of receiving a matrix from the user.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Workspace, character ( len = 255 ) LINE, used to hold the user's input.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
!    Output, integer NCOL, the number of columns in the matrix.
!
!    Output, integer NROW, the number of rows in the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxcol
  integer maxrow

  character ( len = 10 ) chrtmp
  integer ierror
  integer iounit(4)
  character ( len = 255 ) line
  integer ncol
  integer nrow
  character ( len = 255 ) output
  character ( len = 255 ) prompt

  ierror = 0
  nrow = 0
  ncol = 0
 
  prompt = 'number of rows, number of columns.'
!
!  Get NROW, the number of rows.
!
  call i_read ( nrow, line, prompt, iounit, ierror )

  if ( ierror /= 0 ) then
    output = ' '
    call s_write ( iounit, output )
    output = 'LA_INP0 - Fatal error!'
    call s_write ( iounit, output )
    output = '  I_READ returned error flag.'
    call s_write ( iounit, output )
    return
  end if
 
  if ( nrow < 1 ) then
    output = 'Error!  Negative number of rows not allowed!'
    call s_write ( iounit, output )
    ierror = 1
    return
  else if ( maxrow < nrow ) then
    call i_to_s_left ( maxrow, chrtmp )
    output = 'Number of rows must be less than ' // chrtmp
    call s_blanks_delete ( output )
    ierror = 1
    return
  end if
!
!  Get NCOL, the number of columns.
!
  call i_read ( ncol, line, prompt, iounit, ierror )

  if ( ierror /= 0 ) then
    output = ' '
    call s_write ( iounit, output )
    output = 'LA_INP0 - Fatal error!'
    call s_write ( iounit, output )
    output = '  I_READ returned error flag.'
    call s_write ( iounit, output )
    return
  end if
 
  if ( ncol < 1 ) then
    output = 'Error!  Negative number of columns not allowed!'
    call s_write ( iounit, output )
    ierror = 1
    return
  else if ( maxcol < ncol ) then
    call i_to_s_left ( maxcol, chrtmp )
    output = 'Number of columns must be less than ' // chrtmp
    call s_blanks_delete ( output )
    ierror = 1
    return
  end if

  return
end
subroutine la_inp1 ( a, iabot, iatop, ierror, iform, iounit, line, maxcol, &
  maxrow, row1, row2, col1, col2 )

!*****************************************************************************80
!
!! la_inp1() accepts the values of the entries of a matrix from the user.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) A(MAXROW,MAXCOL), the current matrix.
!
!    Output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
!    IATOP and IABOT represent the current rational or decimal
!    matrix.
!
!    Input, integer ICOL.
!    0, enter a single row of the matrix.
!    1, enter all rows of the matrix.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred.
!
!    Input, integer IFORM, 0/1/2, rational/real/decimal arithmetic.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Input, integer IROW.
!    0, enter a single column of the matrix.
!    1, enter all columns of the matrix.
!
!    Workspace, character ( len = 255 ) LINE, used to hold the user's input.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxcol
  integer maxrow

  real ( kind = rk ) a(maxrow,maxcol)
  character ( len = 10 ) chrtmp1
  character ( len = 10 ) chrtmp2
  character ( len = 10 ) chrtmp3
  integer col1
  integer col2
  integer i
  integer iabot(maxrow,maxcol)
  integer iatop(maxrow,maxcol)
  integer ibot
  integer ierror
  integer iform
  integer iounit(4)
  integer itop
  integer j
  integer len1
  integer len2
  character ( len = 255 ) line
  character ( len = 255 ) prompt
  integer row1
  integer row2
  real ( kind = rk ) rval
!
!  Check that 1 <= ROW1 <= ROW2 <= NROW,
!             1 <= COL1 <= COL2 <= NCOL.
!

!
!  Enter a single value
!
  if ( row1 == row2 .and. col1 == col2 ) then
 
    call i_to_s_left ( row1, chrtmp1 )
    len1 = len_trim ( chrtmp1 )
    call i_to_s_left ( col1, chrtmp2 )
    len2 = len_trim ( chrtmp2 )

    prompt = 'A(' // chrtmp1(1:len1) // ',' // chrtmp2(1:len2) //')'
 
    if ( iform == 0 ) then
 
      call rat_read ( itop, ibot, line, prompt, iounit, ierror )
      if ( ierror /= 0 ) then
        return
      end if
      iatop(row1,col1) = itop
      iabot(row1,col1) = ibot
 
    else if ( iform == 1 ) then
 
      call r_read ( rval, line, prompt, iounit, ierror )
      if ( ierror /= 0 ) then
        return
      end if

      a(row1,col1) = rval
 
    else if ( iform == 2 ) then
 
      call dec_read ( itop, ibot, line, prompt, iounit, ierror )
      if ( ierror /= 0 ) then
        return
      end if
 
      call dec_round ( itop, ibot )
 
      iatop(row1,col1) = itop
      iabot(row1,col1) = ibot
 
    end if
!
!  Enter a single column.
!
  else if ( col1 == col2 ) then
 
    do i = row1, row2
 
      call i_to_s_left ( i, chrtmp1 )
      call i_to_s_left ( row2, chrtmp2 )
      call i_to_s_left ( col1, chrtmp3 )
      prompt = 'entries ' // chrtmp1 // ' to ' // chrtmp2 // ' of column ' // &
        chrtmp3
      call s_blanks_delete ( prompt )
 
      if ( iform == 0 ) then
 
        call rat_read ( itop, ibot, line, prompt, iounit, ierror )
        if ( ierror /= 0 ) then
          return
        end if

        iatop(i,col1) = itop
        iabot(i,col1) = ibot
 
      else if ( iform == 1 ) then
 
        call r_read ( rval, line, prompt, iounit, ierror )
        if ( ierror /= 0 ) then
          return
        end if

        a(i,col1) = rval
 
      else if ( iform == 2 ) then
 
        call dec_read ( itop, ibot, line, prompt, iounit, ierror )
        if ( ierror /= 0 ) then
          return
        end if
 
        call dec_round ( itop, ibot )
 
        iatop(i,col1) = itop
        iabot(i,col1) = ibot
 
      end if
 
    end do
!
!  Enter one or more partial rows of at least 2 entries.
!
  else
 
    do i = row1, row2
      do j = col1, col2
 
        call i_to_s_left ( j, chrtmp1 )
        call i_to_s_left ( col2, chrtmp2 )
        call i_to_s_left ( i, chrtmp3 )
        prompt = 'entries ' // chrtmp1 // ' to ' // chrtmp2 // ' of row ' // &
          chrtmp3
        call s_blanks_delete ( prompt )
 
        if ( iform == 0 ) then
 
          call rat_read ( itop, ibot, line, prompt, iounit, ierror )
          if ( ierror /= 0 ) then
            return
          end if

          iatop(i,j) = itop
          iabot(i,j) = ibot
 
        else if ( iform == 1 ) then
 
          call r_read ( rval, line, prompt, iounit, ierror )
          if ( ierror /= 0 ) then
            return
          end if

          a(i,j) = rval
 
        else if ( iform == 2 ) then
 
          call dec_read ( itop, ibot, line, prompt, iounit, ierror )
          if ( ierror /= 0 ) then
            return
          end if
 
          call dec_round ( itop, ibot )
 
          iatop(i,j) = itop
          iabot(i,j) = ibot
 
        end if
 
      end do
    end do
 
  end if
 
  return
end
subroutine la_opt ( a, iabot, iatop, iform, iounit, maxcol, maxrow, ncol, nrow )

!*****************************************************************************80
!
!! la_opt() checks for row echelon or reduced row echelon form.
!
!  Discussion:
!
!    A matrix is in row echelon form if:
!
!    * The first nonzero entry in each row is 1.
!
!    * The leading 1 in a given row occurs in a column to
!      the right of the leading 1 in the previous row.
!
!    * Rows which are entirely zero must occur last.
!
!    The matrix is in reduced row echelon form if, in addition to
!    the first three conditions, it also satisfies:
!
!    * Each column containing a leading 1 has no other nonzero entries.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A(MAXROW,MAXCOL), the current matrix.
!
!    Input, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
!    IATOP and IABOT represent the current rational or decimal
!    matrix.
!
!    Input, integer IFORM, 0/1/2, rational/real/decimal arithmetic.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
!    Input, integer NCOL, the number of columns in the matrix.
!
!    Input, integer NROW, the number of rows in the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxcol
  integer maxrow

  real ( kind = rk ) a(maxrow,maxcol)
  character ( len = 22 ) chrtmp
  character ( len = 10 ) chrtmp1
  character ( len = 10 ) chrtmp2
  integer i
  integer iabot(maxrow,maxcol)
  integer iatop(maxrow,maxcol)
  integer iform
  integer ii
  integer iounit(4)
  integer izer
  integer j
  integer lead
  integer leadp
  integer len1
  integer len2
  integer len3
  integer ncol
  integer nrow
  character ( len = 255 ) output

  output = ' '
  call s_write ( iounit, output )
  output = 'Checking the matrix for row echelon form...'
  call s_write ( iounit, output )
!
!  Check rule 1.
!
  do i = 1, nrow
    do j = 1, ncol
 
      if ( iform == 0 ) then
 
        if ( iatop(i,j) == 0 ) then
          go to 10
        else if ( iatop(i,j) == iabot(i,j) ) then
          go to 20
        end if
 
      else if ( iform == 1 ) then
 
        if ( a(i,j) == 0 ) then
          go to 10
        else if ( a(i,j) == 1 ) then
          go to 20
        end if
 
      else if ( iform == 2 ) then
 
        if ( iatop(i,j) == 0 ) then
          go to 10
        else if ( iatop(i,j) == 1 .and. iabot(i,j) == 0 ) then
          go to 20
        end if
      end if
 
      output = ' '
      call s_write ( iounit, output )
      output = '  The matrix is NOT in row echelon form.'
      call s_write ( iounit, output )
 
      call i_to_s_left ( i, chrtmp1 )
      output = '  The first nonzero entry in row ' // chrtmp1
      call s_write ( iounit, output )
 
      call i_to_s_left ( j, chrtmp1 )
      output = '  which occurs in column ' // chrtmp1
      call s_write ( iounit, output )
 
      if ( iform == 0 ) then
        call rat_to_s_left ( iatop(i,j), iabot(i,j), chrtmp )
      else if ( iform == 1 ) then
        call r_to_s_left ( a(i,j), chrtmp )
      else if ( iform == 2 ) then
        call dec_to_s_left ( iatop(i,j), iabot(i,j), chrtmp )
      end if
 
      len3 = len_trim ( chrtmp )
      output = '  is ' // chrtmp(1:len3) // ' rather than 1.'
      call s_write ( iounit, output )
      return
 
10        continue
 
      end do
 
20      continue
 
    end do
!
!  Check rule 2.
!
  lead = 0
 
  do i = 1, nrow
    do j = 1, ncol
 
      if ( iform == 0 ) then
 
        if ( iatop(i,j) == 0 ) then
          go to 30
        else if ( iatop(i,j) == iabot(i,j) ) then
          leadp = lead
          lead = j
          if ( leadp < lead ) then
            go to 40
          end if
        end if
 
      else if ( iform == 1 ) then
 
        if ( a(i,j) == 0 ) then
          go to 30
        else if ( a(i,j) == 1 ) then
          leadp = lead
          lead = j
          if ( leadp < lead ) then
            go to 40
          end if
        end if
 
      else if ( iform == 2 ) then
 
        if ( iatop(i,j) == 0 ) then
          go to 30
        else if ( iatop(i,j) == 1 .and. iabot(i,j) == 0 ) then
          leadp = lead
          lead = j
          if ( leadp < lead ) then
            go to 40
          end if
        end if
 
      end if
 
      output = ' '
      call s_write ( iounit, output )
      output = '  The matrix is NOT in row echelon form.'
      call s_write ( iounit, output )
      call i_to_s_left ( i, chrtmp1 )
      len1 = len_trim ( chrtmp1 )
      output = '  The first 1 in row ' // chrtmp1(1:len1) // ' does'
      call s_write ( iounit, output )

      call i_to_s_left ( i-1, chrtmp1 )
      len1 = len_trim ( chrtmp1 )
      output = '  NOT occur to the right of the first 1 in row ' // &
        chrtmp1(1:len1) // '.'
      call s_write ( iounit, output )
      return
 
30      continue
 
    end do
 
40    continue
 
  end do
!
!  Check rule 3.
!
  izer = 0
 
  do i = 1, nrow
 
    if ( izer == 0 ) then
 
      do j = 1, ncol
 
        if ( iform == 0 ) then
          if ( iatop(i,j) /= 0 ) then
            go to 70
          end if
        else if ( iform == 1 ) then
          if ( a(i,j) /= 0 ) then
            go to 70
          end if
        else if ( iform == 2 ) then
          if ( iatop(i,j) /= 0 ) then
            go to 70
          end if
        end if
 
      end do
 
      izer = i
 
    else
 
      do j = 1, ncol
 
        if ( ( iform == 0 .and. iatop(i,j) /= 0 ) .or. &
             ( iform == 1 .and. a(i,j) /= 0 ) .or. &
             ( iform == 2 .and. iatop(i,j) /= 0 ) ) then
 
          output = ' '
          call s_write ( iounit, output )

          output = '  The matrix is NOT in row echelon form.'
          call s_write ( iounit, output )

          call i_to_s_left ( izer, chrtmp1 )
          len1 = len_trim ( chrtmp1 )
          output = '  Row ' // chrtmp1(1:len1) //  ' is entirely zero.'
          call s_write ( iounit, output )

          call i_to_s_left ( i, chrtmp1 )
          len1 = len_trim ( chrtmp1 )
          output = '  Row ' // chrtmp1(1:len1) // ' occurs later, and has'
          call s_write ( iounit, output )

          output = '  nonzero entries in it!'
          call s_write ( iounit, output )

          return

        end if
 
      end do
 
    end if
 
70      continue
 
  end do
 
  output = ' '
  call s_write ( iounit, output )
  output = '  The matrix is in row echelon form.'
  call s_write ( iounit, output )
  output = ' '
  call s_write ( iounit, output )
  output = 'Checking the matrix for row reduced echelon form...'
  call s_write ( iounit, output )
!
!  Check rule 4.
!
  do i = 1, nrow
 
    do j = 1, ncol
!
!  We know first nonzero in this row will be 1.
!
      if ( iform == 0 ) then
        if ( iatop(i,j) == 0 ) then
          go to 90
        end if
      else if ( iform == 1 ) then
        if ( a(i,j) == 0 ) then
          go to 90
        end if
      else if ( iform == 2 ) then
        if ( iatop(i,j) == 0 ) then
          go to 90
        end if
      end if
!
!  The leading 1 of this row is entry (i,j).
!
      do ii = 1, nrow
 
        if ( ii /= i ) then
 
          if ( iform == 0 ) then
            if ( iatop(ii,j) == 0 ) then
              go to 80
            end if
          else if ( iform == 1 ) then
            if ( a(ii,j) == 0 ) then
              go to 80
            end if
          else if ( iform == 2 ) then
            if ( iatop(ii,j) == 0 ) then
              go to 80
            end if
          end if
 
          output = ' '
          call s_write ( iounit, output )

          output = '  The matrix is NOT in reduced row echelon form.'
          call s_write ( iounit, output )

          call i_to_s_left ( i, chrtmp1 )
          len1 = len_trim ( chrtmp1 )

          call i_to_s_left ( j, chrtmp2 )
          len2 = len_trim ( chrtmp2 )

          output = '  Row ' // chrtmp1(1:len1) //  ' has its leading 1 in ' &
            // 'column ' // chrtmp2(1:len2) // '.'
          call s_write ( iounit, output )

          output = '  This means that all other entries of that ' // &
            'column should be zero.'
          call s_write ( iounit, output )
 
          if ( iform == 0 ) then

            call rat_to_s_left ( iatop(ii,j), iabot(ii,j), chrtmp )

          else if ( iform == 1 ) then

            call r_to_s_left ( a(ii,j), chrtmp )

          else if ( iform == 2 ) then

            call dec_to_s_left ( iatop(ii,j), iabot(ii,j), chrtmp )

          end if

          call i_to_s_left ( ii, chrtmp1 )

          output = '  But the entry in row ' // trim ( chrtmp1 ) &
           // ' is ' // trim ( chrtmp )
          call s_write ( iounit, output )
          return

        end if
 
80          continue
 
      end do
 
      go to 100
 
90        continue
 
    end do
 
100     continue

  end do
 
  output = ' '
  call s_write ( iounit, output )
  output = 'The matrix is in reduced row echelon form.'
  call s_write ( iounit, output )
 
  return
end
subroutine lp_help ( iounit )

!*****************************************************************************80
!
!! lp_help() prints out a brief list of useful linear programming commands.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer iounit(4)
  character ( len = 255 ) output

  output = ' '
  call s_write ( iounit, output )
  output = 'A(I,J)=S  Set table entry to S.'
  call s_write ( iounit, output )
  output = 'CHECK checks if the solution is optimal.'
  call s_write ( iounit, output )
  output = 'E     Enters a table to work on.'
  call s_write ( iounit, output )
  output = 'HELP  for full help.'
  call s_write ( iounit, output )
  output = 'L     switches to linear algebra.'
  call s_write ( iounit, output )
  output = 'P     I, J performs a pivot operation.'
  call s_write ( iounit, output )
  output = 'Q     quits.'
  call s_write ( iounit, output )
  output = 'TS    types the linear programming solution.'
  call s_write ( iounit, output )
  output = 'V     removes artificial variables.'
  call s_write ( iounit, output )
  output = 'Z     automatic solution (requires password).'
  call s_write ( iounit, output )
  output = '?     interactive help.'
  call s_write ( iounit, output )
  output = ' '
  call s_write ( iounit, output )
  output = 'R1 <=> R2  interchanges two rows'
  call s_write ( iounit, output )
  output = 'R1 <= S R1 multiplies a row by S.'
  call s_write ( iounit, output )
  output = 'R1 <= R1 + S R2 adds a multiple of another row.'
  call s_write ( iounit, output )
  output = ' '
  call s_write ( iounit, output )
 
  return
end
subroutine lp_inp ( a, chineq, iatop, iabot, ibase, ierror, iform, iounit, &
  line, maxcol, maxrow, nart, ncol, ncon, nrow, n_slack, nvar )

!*****************************************************************************80
!
!! lp_inp() allows the user to enter a linear programming problem.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) A(MAXROW,MAXCOL), the current table.
!
!    Output, character CHINEQ(MAXROW), the '<', '=', or '>'
!    sign for each linear programming constraint.
!
!    Input, integer IABOT(MAXROW,MAXCOL), IATOP(MAXROW,MAXCOL).
!    IATOP and IABOT represent the current rational or decimal
!    matrix.
!
!    Output, integer IBASE(MAXROW), the basic variables.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred.
!
!    Input, integer IFORM, 0/1/2, rational/real/decimal arithmetic.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Workspace, character ( len = 255 ) LINE, used to hold the user's input.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
!    Output, integer NART, the number of artificial variables.
!
!    Output, integer NCOL, the number of columns in the matrix.
!
!    Output, integer NCON, the number of constraints.
!
!    Output, integer NROW, the number of rows in the matrix.
!
!    Output, integer N_SLACK, the number of slack variables.
!
!    Output, integer NVAR, the number of basic variables.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxcol
  integer maxrow

  real ( kind = rk ) a(maxrow,maxcol)
  character chineq(maxrow)
  character ( len = 10 ) chrtmp1
  character ( len = 10 ) chrtmp2
  character ( len = 10 ) chrtmp3
  integer i
  integer iabot(maxrow,maxcol)
  integer iart
  integer iatop(maxrow,maxcol)
  integer ibase(maxrow)
  integer ibot
  integer ierror
  integer iform
  integer iounit(4)
  character isay
  integer islak
  integer iterm
  integer itop
  integer j
  integer j1
  integer j2
  integer jcol
  integer jhi
  integer len1
  character ( len = 255 ) line
  integer nart
  integer ncol
  integer ncon
  integer nrow
  integer n_slack
  integer nvar
  character ( len = 255 ) output
  character ( len = 255 ) prompt
  real ( kind = rk ) rval

  nrow = 0
  ncol = 0
  ncon = 0
  nvar = 0
  n_slack = 0
  nart = 0
!
!  Read two integers defining problem.
!
  prompt = 'number of constraints, number of variables.'
!
!  Get number of constraints.
!
  call i_read ( ncon, line, prompt, iounit, ierror )
  if ( ierror /= 0 ) then
    return
  end if
 
  if ( ncon < 0 ) then
    output = 'Number of constraints must be positive!'
    call s_write ( iounit, output )
    ierror = 1
    return
  else if ( maxrow - 2 < ncon ) then

    call i_to_s_left ( maxrow-2, chrtmp1 )
    len1 = len_trim ( chrtmp1 )
    output = 'Number of constraints must be no greater than ' // chrtmp1(1:len1)
    call s_write ( iounit, output )
    ierror = 1
    return

  end if
 
  nrow = ncon + 1
!
!  Get number of variables.
!
  call i_read ( nvar, line, prompt, iounit, ierror )
  if ( ierror /= 0 ) then
    return
  end if
 
  if ( nvar < 1 ) then
    output = 'A negative number of variables is not allowed!'
    call s_write ( iounit, output )
    ierror = 1
    return
  else if ( maxcol < nvar ) then
    call i_to_s_left ( maxcol, chrtmp1 )
    output = 'Number of variables must be no greater than ' // chrtmp1
    call s_write ( iounit, output )
    ierror = 1
    return
  end if

  line = ' '
 
  do i = 1, nrow
 
    chineq(i) = ' '
 
    if ( i <= ncon ) then

      line = ' '

30        continue

      call i_to_s_left ( i, chrtmp1 )
      prompt = 'sign < > or = and coefficients and RHS of constraint' // chrtmp1
      call s_blanks_delete ( prompt )
      iterm = 0
      call s_read ( isay, line, prompt, iounit, ierror, iterm ) 
      if ( ierror /= 0 ) then
        return
      end if
 
      if ( isay == '<' ) then
        ibase(i) = - 1
        n_slack = n_slack+1
      else if ( isay == '=' ) then
        ibase(i) = 1
        nart=nart+1
      else if ( isay == '>' ) then
        ibase(i) = 0
        n_slack = n_slack+1
        nart = nart+1
      else
        output = 'Huh?  Try again.'
        go to 30
      end if
 
      chineq(i) = isay
    else
      line = ' '
      prompt = 'coefficients and constant of objective function.'
    end if
 
    jhi = nvar+1
 
    do j = 1, jhi
 
      jcol = j
      if ( j == jhi ) then
        jcol = maxcol
      end if

      if ( i <= ncon ) then
 
        if ( j < jhi ) then
          call i_to_s_left ( j, chrtmp1 )
          call i_to_s_left ( nvar, chrtmp2 )
          call i_to_s_left ( i, chrtmp3 )
          prompt = 'entries ' // chrtmp1 // ' to ' // chrtmp2 // &
            ' and RHS of constraint ' // chrtmp3
          call s_blanks_delete ( prompt )
        else
          call i_to_s_left ( i, chrtmp1 )
          prompt = 'right hand side of constraint ' // chrtmp1
        end if
 
      end if
 
      call s_blanks_delete ( prompt )
 
      if ( iform == 0 ) then
 
        call rat_read ( itop, ibot, line, prompt, iounit, ierror )
        if ( ierror /= 0 ) then
          return
        end if

        iatop(i,jcol) = itop
        iabot(i,jcol) = ibot
        if ( (i == nrow) .and. (jcol <= nvar) ) then
          iatop(i,jcol) = - iatop(i,jcol)
        end if
 
      else if ( iform == 1 ) then
 
        call r_read ( rval, line, prompt, iounit, ierror )
        if ( ierror /= 0 ) then
          return
        end if

        a(i,jcol) = rval

        if ( i == nrow .and. jcol <= nvar ) then
          a(i,jcol) = - a(i,jcol)
        end if
 
      else if ( iform == 2 ) then
 
        call dec_read ( itop, ibot, line, prompt, iounit, ierror )
        if ( ierror /= 0 ) then
          return
        end if
 
        call dec_round ( itop, ibot )
 
        iatop(i,jcol) = itop
        iabot(i,jcol) = ibot
        if ( i == nrow .and. jcol <= nvar ) then
          iatop(i,jcol) = - iatop(i,jcol)
        end if
 
      end if
 
    end do
 
  end do
!
!  Move the right hand sides to the proper column.
!
  j2 = nvar + n_slack + nart + 2

  do i = 1, nrow
 
    if ( iform == 0 ) then
      iatop(i,j2) = iatop(i,maxcol)
      iabot(i,j2) = iabot(i,maxcol)
    else if ( iform == 1 ) then
      a(i,j2) = a(i,maxcol)
    else if ( iform == 2 ) then
      iatop(i,j2) = iatop(i,maxcol)
      iabot(i,j2) = iabot(i,maxcol)
    end if
 
  end do
!
!  Place the 1 in the bottom of the P column.
!
  j1 = nvar + n_slack + nart + 1

  if ( iform == 0 ) then
    iatop(nrow,j1) = 1
    iabot(nrow,j1) = 1
  else if ( iform == 1 ) then
    a(nrow,j1) = 1.0D+00
  else if ( iform == 2 ) then
    iatop(nrow,j1) = 1
    iabot(nrow,j1) = 0
  end if
!
!  For artificial variable problems, move the objective row down
!  one row to a "hidden" row, and set up a dummy objective row.
!
  if ( 0 < nart ) then
    output = ' '
    call s_write ( iounit, output )
    output = 'Because we have artificial variables, the objective'
    call s_write ( iounit, output )
    output = 'function is also "artificial".'
    call s_write ( iounit, output )
    output = 'The true objective will be stored away until the'
    call s_write ( iounit, output )
    output = 'artificial variables are gone.'
    call s_write ( iounit, output )
    output = ' '
    call s_write ( iounit, output )
 
    do j = 1, nvar + n_slack
 
      if ( iform == 0 ) then
        iatop(nrow+1,j) = iatop(nrow,j)
        iabot(nrow+1,j) = iabot(nrow,j)
        iatop(nrow,j) = 0
        iabot(nrow,j) = 1
      else if ( iform == 1 ) then
        a(nrow+1,j) = a(nrow,j)
        a(nrow,j) = 0.0D+00
      else if ( iform == 2 ) then
        iatop(nrow+1,j) = iatop(nrow,j)
        iabot(nrow+1,j) = iabot(nrow,j)
        iatop(nrow,j) = 0
        iabot(nrow,j) = 0
      end if
 
    end do
!
!  Move the last two entries of the original row to where they
!  would properly line up in the full problem.
!
    j1 = nvar + n_slack + nart + 1
    j2 = nvar + n_slack + nart + 2

    if ( iform == 0 ) then
 
      iatop(nrow+1,j1) = 1
      iabot(nrow+1,j1) = 1
 
      iatop(nrow+1,j2) = iatop(nrow,j2)
      iabot(nrow+1,j2) = iabot(nrow,j2)
 
      iatop(nrow,j2) = 0
      iabot(nrow,j2) = 1
 
    else if ( iform == 1 ) then
 
      a(nrow+1,j1) = 1.0D+00
      a(nrow+1,j2) = a(nrow,j2)
      a(nrow,j2) = 0.0D+00
 
    else if ( iform == 2 ) then
 
      iatop(nrow+1,j1) = 1
      iabot(nrow+1,j1) = 0
 
      iatop(nrow+1,j2) = iatop(nrow,j2)
      iabot(nrow+1,j2) = iabot(nrow,j2)
 
      iatop(nrow,j2) = 0
      iabot(nrow,j2) = 0
 
    end if
 
  end if
!
!  Set entries corresponding to slack and artificial variables.
!
  islak = 0
  iart = 0
  ncol = nvar + n_slack + nart + 2
 
  do i = 1, ncon
 
    if ( ibase(i) == -1 ) then
      islak = islak + 1
      ibase(i) = nvar + islak
 
      if ( iform == 0 ) then
        iatop(i,nvar+islak) = 1
        iabot(i,nvar+islak) = 1
      else if ( iform == 1 ) then
        a(i,nvar+islak) = 1.0D+00
      else if ( iform == 2 ) then
        iatop(i,nvar+islak) = 1
        iabot(i,nvar+islak) = 0
      end if
 
    else if ( ibase(i) == 0 ) then
 
      islak = islak + 1
      iart = iart + 1
      j = nvar + n_slack + iart
      ibase(i) = j
 
      if ( iform == 0 ) then
        iatop(i,nvar+islak) = - 1
        iabot(i,nvar+islak) = 1
        iatop(i,j) = 1
        iabot(i,j) = 1
        iatop(nrow,j) = 1
        iabot(nrow,j) = 1
      else if ( iform == 1 ) then
        a(i,nvar+islak) = - 1.0D+00
        a(i,j) = 1.0D+00
        a(nrow,j) = 1.0D+00
      else if ( iform == 2 ) then
        iatop(i,nvar+islak) = - 1
        iabot(i,nvar+islak) = 0
        iatop(i,j) = 1
        iabot(i,j) = 0
        iatop(nrow,j) = 1
        iabot(nrow,j) = 0
      end if
 
    else if ( ibase(i) == 1 ) then
 
      iart = iart + 1
      j = nvar + n_slack + iart
      ibase(i) = j
 
      if ( iform == 0 ) then
        iatop(i,j) = 1
        iabot(i,j) = 1
        iatop(nrow,j) = 1
        iabot(nrow,j) = 1
      else if ( iform == 1 ) then
        a(i,j) = 1.0D+00
        a(nrow,j) = 1.0D+00
      else if ( iform == 2 ) then
        iatop(i,j) = 1
        iabot(i,j) = 0
        iatop(nrow,j) = 1
        iabot(nrow,j) = 0
      end if
 
    end if
 
  end do
 
  return
end
subroutine lp_opt ( a, iatop, iabot, ibase, iform, iopti, iounit, isltop, &
  islbot, lpmoda, maxcol, maxrow, nart, ncol, nrow, n_slack, nvar, sol )

!*****************************************************************************80
!
!! lp_opt() checks the current linear programming table for optimality.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A(MAXROW,MAXCOL), the current matrix.
!
!    Input, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
!    IATOP and IABOT represent the current rational or decimal matrix.
!
!    Input, integer IBASE(MAXROW), the basic variables.
!
!    Input, integer IFORM, 0/1/2, rational/real/decimal arithmetic.
!
!    Output, integer IOPTI.
!    0, the current solution is not optimal.
!    1, the current solution is optimal.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Output, integer ISLTOP(MAXROW), ISLBOT(MAXROW), the fractional
!    or decimal representation of the linear programming solution.
!
!    Input, integer LPMODA, 0/1, linear algebra/linear programming mode.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
!    Input, integer NART, the number of artificial variables.
!
!    Input, integer NCOL, the number of columns in the matrix.
!
!    Input, integer NROW, the number of rows in the matrix.
!
!    Input, integer N_SLACK, the number of slack variables.
!
!    Input, integer NVAR, the number of basic variables.
!
!    Output, real ( kind = rk ) SOL(MAXROW), the real representation of the
!    linear programming solution.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxcol
  integer maxrow

  real ( kind = rk ) a(maxrow,maxcol)
  character ( len = 22 ) chrtmp
  character ( len = 10 ) chrtmp1
  integer i
  integer iabot(maxrow,maxcol)
  integer iatop(maxrow,maxcol)
  integer ibase(maxrow)
  integer iform
  integer ihi
  integer ilo
  integer iopti
  integer iounit(4)
  integer islbot(maxcol)
  integer isltop(maxcol)
  integer jhi
  integer jlo
  integer lpmoda
  integer nart
  integer ncol
  integer nrow
  integer n_slack
  integer nvar
  character ( len = 255 ) output
  real ( kind = rk ) sol(maxcol)
  real ( kind = rk ) temp
  character ( len = 255 ) title

  output = ' '
  call s_write ( iounit, output )
  output = 'Optimality test'
  call s_write ( iounit, output )
  output = ' '
  call s_write ( iounit, output )
  iopti = 1
  output = ' '
  call s_write ( iounit, output )
  output = 'Are all objective entries nonnegative?'
  call s_write ( iounit, output )
 
  do i = 1, n_slack + nvar + nart
 
    if ( iform == 0 ) then
      call rat_to_r ( temp, iatop(nrow,i), iabot(nrow,i) )
    else if ( iform == 1 ) then
      temp = a(nrow,i)
    else if ( iform == 2 ) then
      call dec_to_r ( temp, iatop(nrow,i), iabot(nrow,i) )
    end if
 
    if ( temp < 0 ) then
      call i_to_s_left ( i, chrtmp1 )
      output = 'Negative objective coefficient, entry ' // chrtmp1
      call s_blanks_delete ( output )
      call s_write ( iounit, output )
      iopti = 0
    end if
 
  end do
 
  output = ' '
  call s_write ( iounit, output )
 
  if ( iopti == 0 ) then
    output = 'The current solution is NOT optimal.'
    call s_write ( iounit, output )
  else
    output = 'Yes.  The current solution is optimal.'
    call s_write ( iounit, output )
  end if
!
!  Print the current linear programming solution.
!
  call lp_sol ( a, iatop, iabot, ibase, iform, isltop, islbot, maxcol, &
    maxrow, ncol, nrow, sol )
 
  title = 'The linear programming solution:'
 
  jhi = nvar + n_slack + nart
  jlo = 1
  ilo = 1
  ihi = 1
 
  if ( iform == 0 ) then

    call rat_print ( isltop, islbot, ibase, iounit, ihi, ilo, jhi, &
      jlo, lpmoda, jhi, 1, ncol, nrow, title )

  else if ( iform == 1 ) then

    call r_print ( sol, ibase, iounit, ihi, ilo, jhi, jlo, lpmoda, &
      jhi, 1, ncol, nrow, title )

  else if ( iform == 2 ) then

    call dec_print ( isltop, islbot, ibase, iounit, ihi, ilo, &
      jhi, jlo, lpmoda, jhi, 1, ncol, nrow, title )

  end if
 
  output = ' '
  call s_write ( iounit, output )
 
  if ( iform == 0 ) then
    call rat_to_s_left ( iatop(nrow,ncol), iabot(nrow,ncol), chrtmp )
  else if ( iform == 1 ) then
    call r_to_s_left ( a(nrow,ncol), chrtmp )
  else if ( iform == 2 ) then
    call dec_to_s_left ( iatop(nrow,ncol), iabot(nrow,ncol), chrtmp )
  end if

  output = 'Objective = ' // chrtmp
  call s_blanks_delete ( output )
  call s_write ( iounit, output )
!
!  Warn user if artificial variables must be deleted.
!
  if ( 0 < nart ) then
    output = ' '
    call s_write ( iounit, output )
    output = 'This problem has artificial variables.'
    call s_write ( iounit, output )
    output = 'Use the "V" command to remove them.'
    call s_write ( iounit, output )
  end if
 
  return
end
subroutine lp_piv ( a, iatop, iabot, iauto, ibase, ierror, iform, iounit, &
  isltop, islbot, line, lpmoda, maxcol, maxrow, nart, ncol, nrow, n_slack, &
  nvar, sol )

!*****************************************************************************80
!
!! lp_piv() carries out pivoting for a linear programming problem.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = rk ) A(MAXROW,MAXCOL), the current matrix.
!
!    Input/output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
!    IATOP and IABOT represent the current rational or decimal matrix.
!
!    Input, integer IAUTO.
!    0, automatic processing is not being carried out.
!    1, automatic processing is being carried out.
!
!    Input/output, integer IBASE(MAXROW), the basic variables.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred.
!
!    Input, integer IFORM, 0/1/2, rational/real/decimal arithmetic.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Output, integer ISLTOP(MAXROW), ISLBOT(MAXROW), the fractional
!    or decimal representation of the linear programming solution.
!
!    Workspace, character ( len = 255 ) LINE, used to hold the user's input.
!
!    Input, integer LPMODA, 0/1, linear algebra/linear programming mode.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
!    Input, integer NART, the number of artificial variables.
!
!    Input, integer NCOL, the number of columns in the matrix.
!
!    Input, integer NROW, the number of rows in the matrix.
!
!    Input, integer N_SLACK, the number of slack variables.
!
!    Input, integer NVAR, the number of basic variables.
!
!    Output, real ( kind = rk ) SOL(MAXROW), the real representation of the
!    linear programming solution.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxcol
  integer maxrow

  real ( kind = rk ) a(maxrow,maxcol)
  character ( len = 22 ) chrtmp
  character ( len = 22 ) chrtmp2
  integer i
  integer iabot(maxrow,maxcol)
  integer iatop(maxrow,maxcol)
  integer iauto
  integer ibase(maxrow)
  integer ierror
  integer iform
  integer ihi
  integer ilo
  integer iobbot
  integer iobtop
  integer iopti
  integer iounit(4)
  integer ipiv
  integer irow
  integer irow1
  integer irow2
  integer isbot
  integer islbot(maxcol)
  integer isltop(maxcol)
  integer istop
  integer j
  integer jhi
  integer jlo
  integer jpiv
  character ( len = 255 ) line
  integer lpmoda
  integer nart
  integer ncol
  integer nrow
  integer n_slack
  integer nvar
  real ( kind = rk ) objnew
  real ( kind = rk ) objold
  character ( len = 255 ) output
  real ( kind = rk ) sol(maxcol)
  real ( kind = rk ) sval
  real ( kind = rk ) temp
  character ( len = 255 ) title

  if ( lpmoda /= 1 ) then
    ierror = 1
    output = 'This command should only be given during'
    call s_write ( iounit, output )
    output = 'linear programming!'
    call s_write ( iounit, output )
    return
  end if
!
!  For each basic variable, check that the objective row entry is zero.
!  If not, then if notautomatic, complain, else fix it.
!
10    continue
 
  call lp_piv1 ( a, iabot, iatop, iauto, ibase, ierror, iform, iounit, &
    maxcol, maxrow, ncol, nrow )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Save current value of objective function.
!
  if ( iform == 0 ) then

    iobtop = iatop(nrow,ncol)
    iobbot = iabot(nrow,ncol)
    call rat_to_r ( objold, iatop(nrow,ncol), iabot(nrow,ncol) )

  else if ( iform == 1 ) then

    objold = a(nrow,ncol)

  else if ( iform == 2 ) then

    iobtop = iatop(nrow,ncol)
    iobbot = iabot(nrow,ncol)
    call dec_to_r ( objold, iatop(nrow,ncol), iabot(nrow,ncol) )

  end if
!
!  Print out the objective row.
!
  if ( iauto == 0 ) then

    title = 'Objective row'
    ilo = nrow
    ihi = nrow
    jlo = 1
    jhi = ncol
 
    if ( iform == 0 ) then

      call rat_print ( iatop, iabot, ibase, iounit, ihi, ilo, jhi, &
        jlo, lpmoda, maxcol, maxrow, ncol, nrow, title )

    else if ( iform == 1 ) then

      call r_print ( a, ibase, iounit, ihi, ilo, jhi, jlo, lpmoda, &
        maxcol, maxrow, ncol, nrow, title )

    else if ( iform == 2 ) then

      call dec_print ( iatop, iabot, ibase, iounit, ihi, ilo, &
        jhi, jlo, lpmoda, maxcol, maxrow, ncol, nrow, title )

    end if
 
  end if
!
!  Check for optimality.
!
  do j = 1, nvar+n_slack+nart
 
    if ( iform == 0 ) then
      call rat_to_r ( temp, iatop(nrow,j), iabot(nrow,j) )
    else if ( iform == 1 ) then
      temp = a(nrow,j)
    else if ( iform == 2 ) then
      call dec_to_r ( temp, iatop(nrow,j), iabot(nrow,j) )
    end if
 
    if ( temp < 0.0D+00 ) then
      go to 30
    end if
 
  end do
 
  call lp_opt ( a, iatop, iabot, ibase, iform, iopti, iounit, isltop, &
    islbot, lpmoda, maxcol, maxrow, nart, ncol, nrow, n_slack, nvar, sol )
 
  return
!
!  Choose the entering variable.
!
30    continue
 
  call lp_piv2 ( a, iabot, iatop, iauto, ierror, iform, iounit, &
    jpiv, line, maxcol, maxrow, nart, nrow, n_slack, nvar )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Choose the departing variable.
!
  call lp_piv3 ( a, iabot, iatop, iauto, ibase, ierror, iform, &
    iounit, ipiv, jpiv, line, maxcol, maxrow, ncol, nrow )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Pivot on entry (IPIV,JPIV).
!
  irow = ipiv
 
  if ( iform == 0 ) then
    istop = iatop(ipiv,jpiv)
    isbot = iabot(ipiv,jpiv)
  else if ( iform == 1 ) then
    sval = a(ipiv,jpiv)
  else if ( iform == 2 ) then
    istop = iatop(ipiv,jpiv)
    isbot = iabot(ipiv,jpiv)
  end if
 
  call row_div ( a, iatop, iabot, ierror, iform, iounit, irow, &
    maxcol, maxrow, ncol, nrow, sval, istop, isbot )
 
  irow2 = ipiv
 
  do i = 1, nrow
 
    irow1 = i
 
    if ( irow1 /= ipiv ) then
 
      if ( iform == 0 ) then
        istop = - iatop(irow1,jpiv)
        isbot = iabot(irow1,jpiv)
      else if ( iform == 1 ) then
        sval = -a(irow1,jpiv)
      else if ( iform == 2 ) then
        istop = -iatop(irow1,jpiv)
        isbot = iabot(irow1,jpiv)
      end if
 
      call row_add ( a, iatop, iabot, ierror, iform, iounit, &
        irow1, irow2, maxcol, maxrow, ncol, sval, istop, isbot )

    end if
 
  end do
!
!  Print out change in objective.
!
  output = ' '
  call s_write ( iounit, output )
 
  output = 'No change in objective.'
 
  if ( iform == 0 ) then
 
    call rat_to_r ( objnew, iatop(nrow,ncol), iabot(nrow,ncol) )
 
    if ( objold /= objnew ) then
 
      call rat_to_s_left ( iobtop, iobbot, chrtmp )
 
      if ( iobbot /= 1 ) then
        call r_to_s_left ( objold, chrtmp2 )
        output = 'The objective changed from ' // chrtmp // ' = ' // chrtmp2
      else
        output = 'The objective changed from ' // chrtmp
      end if
 
      call s_blanks_delete ( output )
      call s_write ( iounit, output )
 
      call rat_to_s_left ( iatop(nrow,ncol), iabot(nrow,ncol), chrtmp )
 
      if ( iabot(nrow,ncol) /= 1 ) then
        call r_to_s_left ( objnew, chrtmp2 )
        output = 'to ' // chrtmp // ' = ' // chrtmp2
      else
        output = 'to ' // chrtmp
      end if
 
      call s_blanks_delete ( output )
      call s_write ( iounit, output )
 
    end if
 
  else if ( iform == 1 ) then
 
    objnew = a(nrow,ncol)
 
    if ( objold /= objnew ) then
      call r_to_s_left ( objold, chrtmp )
      call r_to_s_left ( objnew, chrtmp2 )
      output = 'The objective changed from ' // chrtmp // ' to ' // chrtmp2
      call s_blanks_delete ( output )
      call s_write ( iounit, output )
    end if
 
  else if ( iform == 2 ) then
 
    call dec_to_r ( objnew, iatop(nrow,ncol), iabot(nrow,ncol) )
 
    if ( objold /= objnew ) then
 
      call dec_to_s_left ( iobtop, iobbot, chrtmp )
      call dec_to_s_left ( iatop(nrow,ncol), iabot(nrow,ncol), chrtmp2 )
 
      output = 'The objective changed from ' // chrtmp // ' to ' // chrtmp2
      call s_blanks_delete ( output )
      call s_write ( iounit, output )
 
    end if
 
  end if
 
  if ( iauto == 1 ) then
    go to 10
  end if
 
  return
end
subroutine lp_piv1 ( a, iabot, iatop, iauto, ibase, ierror, iform, iounit, &
  maxcol, maxrow, ncol, nrow )

!*****************************************************************************80
!
!! lp_piv1() zeroes out objective row entries for basic variables.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = rk ) A(MAXROW,MAXCOL), the current matrix.
!
!    Input/output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
!    IATOP and IABOT represent the current rational or decimal matrix.
!
!    Input, integer IAUTO.
!    0, automatic processing is not being carried out.
!    1, automatic processing is being carried out.
!
!    Input, integer IBASE(MAXROW), the basic variables.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred.
!
!    Input, integer IFORM, 0/1/2, rational/real/decimal arithmetic.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
!    Input, integer NCOL, the number of columns in the matrix.
!
!    Input, integer NROW, the number of rows in the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxcol
  integer maxrow

  real ( kind = rk ) a(maxrow,maxcol)
  real ( kind = rk ) amax
  character ( len = 10 ) chrtmp1
  character ( len = 10 ) chrtmp2
  integer i
  integer iabot(maxrow,maxcol)
  integer iatop(maxrow,maxcol)
  integer iauto
  integer ibase(maxrow)
  integer ierror
  integer iform
  integer imax
  integer iounit(4)
  integer irow
  integer irow1
  integer irow2
  integer isbot
  integer istop
  integer jcol
  integer ncol
  integer nrow
  character ( len = 255 ) output
  real ( kind = rk ) sval
  real ( kind = rk ) temp

  do i = 1, nrow - 1
 
    jcol = ibase(i)
 
    if ( jcol < 1 .or. ncol < jcol ) then
      output = 'Error in the IBASE vector!'
      call s_write ( iounit, output )
      call i_to_s_left ( i, chrtmp1 )
      call i_to_s_left ( jcol, chrtmp2 )
      output = 'Entry ' // chrtmp1 // ' of IBASE = ' // chrtmp2
      call s_write ( iounit, output )
      ierror = 1
      return
    end if
!
!  Check the objective entry in column JCOL.
!
    if ( iform == 0 ) then
      call rat_to_r ( temp, iatop(nrow,jcol), iabot(nrow,jcol) )
    else if ( iform == 1 ) then
      temp = a(nrow,jcol)
    else if ( iform == 2 ) then
      call dec_to_r ( temp, iatop(nrow,jcol), iabot(nrow,jcol) )
    end if
 
    if ( temp /= 0 ) then
 
      irow1 = nrow
      output = ' '
      call s_write ( iounit, output )
      call i_to_s_left ( jcol, chrtmp1 )
      output = 'The objective entry in column ' // chrtmp1 // ' is not zero,'
      call s_blanks_delete ( output )
      call s_write ( iounit, output )
      output = 'but this corresponds to a basic variable.'
      call s_write ( iounit, output )
      output = ' '
      call s_write ( iounit, output )
 
      if ( iauto == 0 ) then
        output = 'Use the "A" command to zero out this entry.'
        call s_write ( iounit, output )
        output = 'THEN you may use the "P" command!'
        call s_write ( iounit, output )
        ierror = 1
        return
      end if
!
!  Search for the IMAX, the row of the maximum entry in column JCOL.
!
      amax = 0.0D+00
      imax = 0
 
      do irow = 1, nrow-1
 
        if ( iform == 0 ) then
          call rat_to_r ( temp, iatop(irow,jcol), iabot(irow,jcol) )
        else if ( iform == 1 ) then
          temp = a(irow,jcol)
        else if ( iform == 2 ) then
          call dec_to_r ( temp, iatop(irow,jcol), iabot(irow,jcol) )
        end if
 
        temp = abs ( temp )
 
        if ( amax < temp ) then
          amax = temp
          imax = irow
        end if
 
      end do
 
      if ( amax == 0.0D+00 ) then
        output = 'The artificial variable cannot be eliminated!'
        call s_write ( iounit, output )
        go to 20
      end if
!
!  Divide row IMAX by entry (IMAX,JCOL) to normalize it.
!
      irow = imax
 
      if ( iform == 0 ) then
        istop = iatop(irow,jcol)
        isbot = iabot(irow,jcol)
      else if ( iform == 1 ) then
        sval = a(irow,jcol)
      else if ( iform == 2 ) then
        istop = iatop(irow,jcol)
        isbot = iabot(irow,jcol)
      end if
 
      call row_div ( a, iatop, iabot, ierror, iform, iounit, irow, &
        maxcol, maxrow, ncol, nrow, sval, istop, isbot )
!
!  Add a multiple of row IMAX to row NROW, to eliminate entry (NROW,JCOL).
!
      irow2 = imax
      irow1 = nrow
 
      if ( iform == 0 ) then
        istop = -iatop(irow1,jcol)
        isbot = iabot(irow1,jcol)
      else if ( iform == 1 ) then
        sval = -a(irow1,jcol)
      else if ( iform == 2 ) then
        istop = -iatop(irow1,jcol)
        isbot = iabot(irow1,jcol)
      end if
 
      call row_add ( a, iatop, iabot, ierror, iform, iounit, &
        irow1, irow2, maxcol, maxrow, ncol, sval, istop, isbot )

      if ( iform == 0 ) then
        iatop(nrow,jcol) = 0
        iabot(nrow,jcol) = 1
      else if ( iform == 1 ) then
        a(nrow,jcol) = 0.0D+00
      else if ( iform == 2 ) then
        iatop(nrow,jcol) = 0
        iabot(nrow,jcol) = 0
      end if
 
    end if
 
20      continue
 
  end do
 
  return
end
subroutine lp_piv2 ( a, iabot, iatop, iauto, ierror, iform, &
  iounit, jpiv, line, maxcol, maxrow, nart, nrow, n_slack, nvar )

!*****************************************************************************80
!
!! lp_piv2() chooses the entering variable for pivoting.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A(MAXROW,MAXCOL), the current matrix.
!
!    Input, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
!    IATOP and IABOT represent the current rational or decimal matrix.
!
!    Input, integer IAUTO.
!    0, automatic processing is not being carried out.
!    1, automatic processing is being carried out.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred.
!
!    Input, integer IFORM, 0/1/2, rational/real/decimal arithmetic.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Output, integer JPIV, the entering variable.
!
!    Workspace, character ( len = 255 ) LINE, used to hold the user's input.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
!    Input, integer NART, the number of artificial variables.
!
!    Input, integer NROW, the number of rows in the matrix.
!
!    Input, integer N_SLACK, the number of slack variables.
!
!    Input, integer NVAR, the number of basic variables.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxcol
  integer maxrow

  real ( kind = rk ) a(maxrow,maxcol)
  character ( len = 10 ) chrtmp
  integer iabot(maxrow,maxcol)
  integer iatop(maxrow,maxcol)
  integer iauto
  integer ierror
  integer iform
  integer iounit(4)
  integer j
  integer jmin
  integer jpiv
  character ( len = 255 ) line
  integer nart
  integer nrow
  integer n_slack
  integer nvar
  character ( len = 255 ) output
  character ( len = 255 ) prompt
  real ( kind = rk ) temp
  real ( kind = rk ) tmin

10    continue
!
!  Set TMIN to the (NROW,1) entry, and JMIN to 1.
!
  if ( iform == 0 ) then
    call rat_to_r ( tmin, iatop(nrow,1), iabot(nrow,1) )
  else if ( iform == 1 ) then
    tmin = a(nrow,1)
  else if ( iform == 2 ) then
    call dec_to_r ( tmin, iatop(nrow,1), iabot(nrow,1) )
  end if
 
  jmin = 1
!
!  Set TMIN to the (NROW,J) entry, if it is smaller.
!
  do j = 2, nvar+n_slack+nart
 
    if ( iform == 0 ) then
      call rat_to_r ( temp, iatop(nrow,j), iabot(nrow,j) )
    else if ( iform == 1 ) then
      temp = a(nrow,j)
    else if ( iform == 2 ) then
      call dec_to_r ( temp, iatop(nrow,j), iabot(nrow,j) )
    end if
 
    if ( temp <= tmin ) then
      tmin = temp
      jmin = j
    end if
 
  end do
!
!  Now get pivot index JPIV.
!
  if ( iauto == 1 ) then
 
    jpiv = jmin
 
  else
 
    output = ' '
    call s_write ( iounit, output )
    output = 'Variable with most negative objective coefficient?'
    call s_write ( iounit, output )
 
    prompt = 'column (=variable number)'
    call i_read ( jpiv, line, prompt, iounit, ierror )
    if ( ierror /= 0 ) then
      return
    end if
 
    if ( jpiv < 1 .or. nvar + n_slack + nart < jpiv ) then
      output = 'Your input was out of bounds.'
      call s_write ( iounit, output )
      go to 10
    end if
 
    if ( iform == 0 ) then
      call rat_to_r ( temp, iatop(nrow,jpiv), iabot(nrow,jpiv) )
    else if ( iform == 1 ) then
      temp = a(nrow,jpiv)
    else if ( iform == 2 ) then
      call dec_to_r ( temp, iatop(nrow,jpiv), iabot(nrow,jpiv) )
    end if
 
    if ( tmin + 0.0001D+00 < temp ) then
      output = 'Not acceptable.'
      call s_write ( iounit, output )
      go to 10
    end if
 
  end if
 
  output = ' '
  call s_write ( iounit, output )
  call i_to_s_left ( jpiv, chrtmp )
  output = 'The entering variable is ' // chrtmp
  call s_blanks_delete ( output )
  call s_write ( iounit, output )
 
  return
end
subroutine lp_piv3 ( a, iabot, iatop, iauto, ibase, ierror, iform, &
  iounit, ipiv, jpiv, line, maxcol, maxrow, ncol, nrow )

!*****************************************************************************80
!
!! lp_piv3() chooses the departing variable for pivoting.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A(MAXROW,MAXCOL), the current matrix.
!
!    Input, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
!    IATOP and IABOT represent the current rational or decimal matrix.
!
!    Input, integer IAUTO.
!    0, automatic processing is not being carried out.
!    1, automatic processing is being carried out.
!
!    Input, integer IBASE(MAXROW), the basic variables.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred.
!
!    Input, integer IFORM, 0/1/2, rational/real/decimal arithmetic.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Output, integer IPIV, the row of the departing variable.
!
!    Input, integer JPIV, the entering variable.
!
!    Workspace, character ( len = 255 ) LINE, used to hold the user's input.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
!    Input, integer NCOL, the number of columns in the matrix.
!
!    Input, integer NROW, the number of rows in the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxcol
  integer maxrow

  real ( kind = rk ) a(maxrow,maxcol)
  real ( kind = rk ) bot
  character ( len = 22 ) chrtmp
  character ( len = 10 ) chrtmp1
  integer i
  integer iabot(maxrow,maxcol)
  integer iatop(maxrow,maxcol)
  integer iauto
  integer ibase(maxrow)
  integer ibot
  integer ierror
  integer iform
  integer imin
  integer iounit(4)
  integer ipiv
  integer itop
  integer jpiv
  character ( len = 255 ) line
  integer ncol
  integer nrow
  character ( len = 255 ) output
  character ( len = 255 ) prompt
  real ( kind = rk ) ratio
  real ( kind = rk ) ratj
  real ( kind = rk ) ratmin
  real ( kind = rk ) temp1
  real ( kind = rk ) temp2
  real ( kind = rk ) top

  if ( iauto == 0 ) then
    output = ' '
    call s_write ( iounit, output )
    output = 'Variable with smallest nonnegative feasibility ratio?'
    call s_write ( iounit, output )
  end if
 
10    continue
 
  imin = 0

  if ( iauto == 0 ) then
    output = ' '
    call s_write ( iounit, output )
    output = 'Nonnegative feasibility ratios:'
    call s_write ( iounit, output )
    output = ' '
    call s_write ( iounit, output )
  end if
 
  ratmin = - 1.0D+00
 
  do i = 1, nrow - 1
 
    if ( iform == 0 ) then
 
      if ( iabot(i,jpiv) < 0 ) then
        iatop(i,jpiv) = - iatop(i,jpiv)
        iabot(i,jpiv) = - iabot(i,jpiv)
      end if
 
      if ( iatop(i,jpiv) <= 0 ) then
        go to 20
      end if
 
      call rat_to_r ( top, iatop(i,ncol), iabot(i,ncol) )
      call rat_to_r ( bot, iatop(i,jpiv), iabot(i,jpiv) )
 
    else if ( iform == 1 ) then
 
      if ( a(i,jpiv) <= 0.0D+00 ) then
        go to 20
      end if
 
      top = a(i,ncol)
      bot = a(i,jpiv)
 
    else if ( iform == 2 ) then
 
      if ( iatop(i,jpiv) <= 0 ) then
        go to 20
      end if
 
      call dec_to_r ( top, iatop(i,ncol), iabot(i,ncol) )
      call dec_to_r ( bot, iatop(i,jpiv), iabot(i,jpiv) )
 
    end if
 
    if ( bot == 0.0D+00 ) then
      go to 20
    end if
 
    ratio = top / bot
 
    if ( iauto == 0 ) then
 
      if ( iform == 0 ) then
 
        call rat_mul ( ibot, iabot(i,ncol), iatop(i,jpiv), &
          itop, iatop(i,ncol), iabot(i,jpiv), ierror )
 
        if ( ibot /= 1 ) then
          call rat_to_s_left ( itop, ibot, chrtmp )
          write ( output, 120 ) i, ibase(i), ratio, chrtmp
        else
          write ( output, 111 ) i, ibase(i), ratio
        end if
 
      else if ( iform == 1 ) then
 
        write ( output, 111 ) i, ibase(i), ratio
 
      else if ( iform == 2 ) then
 
        write ( output, 111 ) i, ibase(i), ratio
 
      end if
 
      call s_blanks_delete ( output )
      call s_write ( iounit, output )
 
    end if
 
    if ( imin == 0 .and. 0.0D+00 <= ratio ) then
      imin = i
      ratmin = ratio
    end if
 
    if ( 0.0D+00 <= ratio .and. ratio < ratmin ) then
      ratmin = ratio
      imin = i
    end if
 
20      continue
 
  end do
 
  if ( imin == 0 ) then
    output = ' '
    call s_write ( iounit, output )
    output = 'Cannot find a departing variable.'
    call s_write ( iounit, output )
    output = 'Presumably, the feasible set is unbounded.'
    call s_write ( iounit, output )
    ierror = 1
    return
  end if
 
30    continue
 
  if ( iauto == 1 ) then

    ipiv = imin

  else

    output = ' '
    call s_write ( iounit, output )
    prompt = 'the row of the departing variable.'
    call i_read ( ipiv, line, prompt, iounit, ierror )
    if ( ierror /= 0 ) then
      return
    end if
 
    if ( ipiv <= 0 .or. nrow -1 < ipiv ) then
      output = 'Illegal row.'
      call s_write ( iounit, output )
      go to 30
    end if
 
    if ( iform == 0 ) then
 
      if ( iatop(ipiv,jpiv) == 0 ) then
        output = 'Illegal zero divisor.'
        call s_write ( iounit, output )
        go to 10
      end if
 
    else if ( iform == 1 ) then
 
      if ( a(ipiv,jpiv) == 0 ) then
        output = 'Illegal zero divisor.'
        call s_write ( iounit, output )
        go to 10
      end if
 
    else if ( iform == 2 ) then
 
      if ( iatop(ipiv,jpiv) == 0 ) then
        output = 'Illegal zero divisor.'
        call s_write ( iounit, output )
        go to 10
      end if
 
    end if
 
    if ( iform == 0 ) then
 
      call rat_to_r ( temp1, iatop(ipiv,ncol), iabot(ipiv,ncol) )
      call rat_to_r ( temp2, iatop(ipiv,jpiv), iabot(ipiv,jpiv) )
 
    else if ( iform == 1 ) then
 
      temp1 = a(ipiv,ncol)
      temp2 = a(ipiv,jpiv)
 
    else if ( iform == 2 ) then
 
      call dec_to_r ( temp1, iatop(ipiv,ncol), iabot(ipiv,ncol) )
      call dec_to_r ( temp2, iatop(ipiv,jpiv), iabot(ipiv,jpiv) )
 
    end if
 
    ratj = temp1 / temp2
 
    if ( ratj < 0.0D+00 ) then
      output = 'The pivot ratio is not acceptable because ' // &
        'it is negative.'
      call s_write ( iounit, output )
      go to 10
    else if ( ratmin < ratj ) then
      output = 'The pivot ratio is not acceptable because ' // &
        'it is not the smallest nonnegative ratio.'
      call s_write ( iounit, output )
      go to 10
    end if
 
  end if
 
  output = ' '
  call s_write ( iounit, output )
  call i_to_s_left ( ibase(ipiv), chrtmp1 )
  call r_to_s_left ( ratmin, chrtmp )
  output = 'The departing variable is ' // chrtmp1 &
    // ' with feasibility ratio ' // chrtmp
  call s_blanks_delete ( output )
  call s_write ( iounit, output )
 
  output = ' '
  call s_write ( iounit, output )
 
  ibase(ipiv) = jpiv
 
  111 format('Row ',i2,', variable ',i2,', ratio = ',g14.6)
  120 format('Row ',i2,', variable ',i2,', ratio = ',g14.6,' = ', a22 )
 
  return
end
subroutine lp_rem ( a, iabot, iatop, ibase, ierror, iform, iounit, lpmoda, &
  maxcol, maxrow, nart, ncol, nrow, n_slack, nvar )

!*****************************************************************************80
!
!! lp_rem() removes the artificial variables.
!
!  Discussion:
!
!    The routine should be called once the artificial objective function has 
!    reached zero.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = rk ) A(MAXROW,MAXCOL), the current matrix.
!
!    Input/output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
!    IATOP and IABOT represent the current rational or decimal matrix.
!
!    Input, integer IBASE(MAXROW), the basic variables.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred.
!
!    Input, integer IFORM, 0/1/2, rational/real/decimal arithmetic.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Input, integer LPMODA, 0/1, linear algebra/linear programming mode.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
!    Input, integer NART, the number of artificial variables.
!
!    Input/output, integer NCOL, the number of columns in the matrix.
!    On output, this number may have changed because of the elimination
!    of artificial variables.
!
!    Input, integer NROW, the number of rows in the matrix.
!
!    Input, integer N_SLACK, the number of slack variables.
!
!    Input, integer NVAR, the number of basic variables.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxcol
  integer maxrow

  real ( kind = rk ) a(maxrow,maxcol)
  character ( len = 10 ) chrtmp1
  integer i
  integer iabot(maxrow,maxcol)
  integer iatop(maxrow,maxcol)
  integer ibase(maxrow)
  integer ierror
  integer iform
  integer inext
  integer iounit(4)
  integer j
  integer jhi
  integer jvar
  integer lpmoda
  integer mart
  integer nart
  integer ncol
  integer nrow
  integer n_slack
  integer nvar
  character ( len = 255 ) output

  if ( lpmoda /= 1 ) then
    ierror = 1
    output = 'This command should only be given during'
    call s_write ( iounit, output )
    output = 'linear programming!'
    call s_write ( iounit, output )
    return
  end if
 
  if ( nart == 0 ) then
    output = 'There aren''t any artificial variables to delete!'
    call s_write ( iounit, output )
    return
  end if
 
  if ( ( iform == 0 .and. iatop(nrow,ncol) /= 0 ) .or. &
       ( iform == 1 .and. a(nrow,ncol) /= 0.0D+00 ) .or. &
       ( iform == 2 .and. iatop(nrow,ncol) /= 0 ) ) then

    output = 'The phase 1 objective function is nonzero.'
    call s_write ( iounit, output )
    output = 'Hence, this problem may have no solution.'
    call s_write ( iounit, output )

  end if
 
  jhi = nvar + n_slack + nart + 2
  mart = nart 
  nart = 0
  inext = nvar + n_slack
 
  do jvar = nvar + n_slack + 1, jhi
 
    if ( jhi - 2 < jvar ) then
      go to 30
    end if
 
      do i = 1, nrow - 1
        if ( ibase(i) == jvar ) then
          nart = nart+1
          go to 30
        end if
      end do
 
      do i = 1, nrow - 1
        if ( jvar < ibase(i) ) then
          ibase(i) = ibase(i) - 1
        end if
      end do
 
    go to 50

30  continue

      inext = inext + 1
 
      if ( iform == 0 ) then
        iatop(1:nrow,inext) = iatop(1:nrow,jvar)
        iabot(1:nrow,inext) = iabot(1:nrow,jvar)
      else if ( iform == 1 ) then
        a(1:nrow,inext) = a(1:nrow,jvar)
      else if ( iform == 2 ) then
        iatop(1:nrow,inext) = iatop(1:nrow,jvar)
        iabot(1:nrow,inext) = iabot(1:nrow,jvar)
      end if
 
50  continue
 
  end do
!
!  If possible, restore the original objective function.
!
  ncol = nvar + n_slack + 2
 
  output = ' '
  call s_write ( iounit, output )
 
  if ( nart /= 0 ) then
    call i_to_s_left ( nart, chrtmp1 )
    output = chrtmp1 // ' artificial variables were not deleted.'
    call s_blanks_delete ( output )
    call s_write ( iounit, output )
    output = 'You must revise the objective row by hand!'
    call s_write ( iounit, output )
  else
    output = 'All the artificial variables were deleted.'
    call s_write ( iounit, output )
    output = 'The original objective function is restored.'
    call s_write ( iounit, output )
 
    do j = 1, nvar+n_slack
 
      if ( iform == 0 ) then
        iatop(nrow,j) = iatop(nrow+1,j)
        iabot(nrow,j) = iabot(nrow+1,j)
      else if ( iform == 1 ) then
        a(nrow,j) = a(nrow+1,j)
      else if ( iform == 2 ) then
        iatop(nrow,j) = iatop(nrow+1,j)
        iabot(nrow,j) = iabot(nrow+1,j)
      end if
 
    end do
 
    do j = nvar+n_slack+1, nvar+n_slack+2
 
      if ( iform == 0 ) then
        iatop(nrow,j) = iatop(nrow+1,j+mart)
        iabot(nrow,j) = iabot(nrow+1,j+mart)
      else if ( iform == 1 ) then
        a(nrow,j) = a(nrow+1,j+mart)
      else if ( iform == 2 ) then
        iatop(nrow,j) = iatop(nrow+1,j+mart)
        iabot(nrow,j) = iabot(nrow+1,j+mart)
      end if
 
    end do
 
  end if
 
  output = ' '
  call s_write ( iounit, output )
  output = 'You must now use the "A" command to zero out'
  call s_write ( iounit, output )
  output = 'objective row entries for all basic variables.'
  call s_write ( iounit, output )
 
  return
end
subroutine lp_sama ( a, chineq, iatop, iabot, ibase, iform, imat, iounit, &
   maxcol, maxrow, nart, ncol, nrow, n_slack, nvar )

!*****************************************************************************80
!
!! lp_sama() sets up an advanced linear programming problem.
!
!  Discussion:
!
!    The problem includes artificial variables.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) A(MAXROW,MAXCOL), the current matrix.
!
!    Output, character CHINEQ(MAXROW), the '<', '=', or '>'
!    sign for each linear programming constraint.
!
!    Output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
!    IATOP and IABOT represent the current rational or decimal matrix.
!
!    Output, integer IBASE(MAXROW), the basic variables.
!
!    Input, integer IFORM, 0/1/2, rational/real/decimal arithmetic.
!
!    Output, integer IMAT, 0/1, a matrix HAS NOT/HAS been defined.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
!    Output, integer NART, the number of artificial variables.
!
!    Output, integer NCOL, the number of columns in the matrix.
!
!    Output, integer NROW, the number of rows in the matrix.
!
!    Output, integer N_SLACK, the number of slack variables.
!
!    Output, integer NVAR, the number of basic variables.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxcol
  integer maxrow

  real ( kind = rk ) a(maxrow,maxcol)
  character chineq(maxrow)
  integer i
  integer iabot(maxrow,maxcol)
  integer iatop(maxrow,maxcol)
  integer ibase(maxrow)
  integer iform
  integer imat
  integer iounit(4)
  integer j
  integer nart
  integer ncol
  integer nrow
  integer n_slack
  integer nvar
  character ( len = 255 ) output

  nvar = 2
  n_slack = 4
  nart = 2
  nrow = n_slack + 1
  ncol = nvar + n_slack + nart + 2
 
  iatop(1,1) = 1
  iatop(1,2) = 2
  iatop(1,3) = -1
  iatop(1,4) = 0
  iatop(1,5) = 0
  iatop(1,6) = 0
  iatop(1,7) = 1
  iatop(1,8) = 0
  iatop(1,9) = 0
  iatop(1,10) = 6
 
  iatop(2,1) = 2
  iatop(2,2) = 1
  iatop(2,3) = 0
  iatop(2,4) = -1
  iatop(2,5) = 0
  iatop(2,6) = 0
  iatop(2,7) = 0
  iatop(2,8) = 1
  iatop(2,9) = 0
  iatop(2,10) = 4
 
  iatop(3,1) = 1
  iatop(3,2) = 1
  iatop(3,3) = 0
  iatop(3,4) = 0
  iatop(3,5) = 1
  iatop(3,6) = 0
  iatop(3,7) = 0
  iatop(3,8) = 0
  iatop(3,9) = 0
  iatop(3,10) = 5
 
  iatop(4,1) = 2
  iatop(4,2) = 1
  iatop(4,3) = 0
  iatop(4,4) = 0
  iatop(4,5) = 0
  iatop(4,6) = 1
  iatop(4,7) = 0
  iatop(4,8) = 0
  iatop(4,9) = 0
  iatop(4,10) = 8
 
  iatop(5,1) = 0
  iatop(5,2) = 0
  iatop(5,3) = 0
  iatop(5,4) = 0
  iatop(5,5) = 0
  iatop(5,6) = 0
  iatop(5,7) = 1
  iatop(5,8) = 1
  iatop(5,9) = 1
  iatop(5,10) = 0
 
  iatop(6,1) = -40
  iatop(6,2) = -30
  iatop(6,3) = 0
  iatop(6,4) = 0
  iatop(6,5) = 0
  iatop(6,6) = 0
  iatop(6,7) = 0
  iatop(6,8) = 0
  iatop(6,9) = 1
  iatop(6,10) = 0
 
  do i = 1, nrow+1
    do j = 1, ncol
      if ( iform == 0 ) then
        iabot(i,j) = 1
      else if ( iform == 2 ) then
        iabot(i,j) = 0
      end if
    end do
  end do
 
  a(1:nrow+1,1:ncol) = real ( iatop(1:nrow+1,1:ncol), kind = rk )
 
  ibase(1) = 7
  ibase(2) = 8
  ibase(3) = 5
  ibase(4) = 6
 
  chineq(1) = '>'
  chineq(2) = '>'
  chineq(3) = '<'
  chineq(4) = '<'
  chineq(5) = ' '
 
  imat = 1
 
  output = ' '
  call s_write ( iounit, output )
  output = 'Advanced linear programming problem:'
  call s_write ( iounit, output )
  output = ' '
  call s_write ( iounit, output )
  output = 'Maximize'
  call s_write ( iounit, output )
  output = '  Z=40 X + 30 Y'
  call s_write ( iounit, output )
  output = 'subject to'
  call s_write ( iounit, output )
  output = '  X + 2 Y > 6'
  call s_write ( iounit, output )
  output = '2 X +   Y > 4'
  call s_write ( iounit, output )
  output = '  X +   Y < 5'
  call s_write ( iounit, output )
  output = '2 X +   Y < 8'
  call s_write ( iounit, output )
  output = ' '
  call s_write ( iounit, output )
 
  return
end
subroutine lp_sams ( a, chineq, iatop, iabot, ibase, iform, imat, iounit, & 
  maxcol, maxrow, nart, ncol, nrow, n_slack, nvar )

!*****************************************************************************80
!
!! lp_sams() sets up a simple linear programming problem.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) A(MAXROW,MAXCOL), the current matrix.
!
!    Output, character CHINEQ(MAXROW), the '<', '=', or '>'
!    sign for each linear programming constraint.
!
!    Output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
!    IATOP and IABOT represent the current rational or decimal matrix.
!
!    Output, integer IBASE(MAXROW), the basic variables.
!
!    Input, integer IFORM, 0/1/2, rational/real/decimal arithmetic.
!
!    Output, integer IMAT, 0/1, a matrix HAS NOT/HAS been defined.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
!    Output, integer NART, the number of artificial variables.
!
!    Output, integer NCOL, the number of columns in the matrix.
!
!    Output, integer NROW, the number of rows in the matrix.
!
!    Output, integer N_SLACK, the number of slack variables.
!
!    Output, integer NVAR, the number of basic variables.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxcol
  integer maxrow

  real ( kind = rk ) a(maxrow,maxcol)
  character chineq(maxrow)
  integer i
  integer iabot(maxrow,maxcol)
  integer iatop(maxrow,maxcol)
  integer ibase(maxrow)
  integer iform
  integer imat
  integer iounit(4)
  integer j
  integer nart
  integer ncol
  integer nrow
  integer n_slack
  integer nvar
  character ( len = 255 ) output

  nvar = 2
  n_slack = 2
  nart = 0
  nrow = n_slack + 1
  ncol = nvar + n_slack + nart + 2
 
  iatop(1,1) = 2
  iatop(1,2) = 2
  iatop(1,3) = 1
  iatop(1,4) = 0
  iatop(1,5) = 0
  iatop(1,6) = 8
 
  iatop(2,1) = 5
  iatop(2,2) = 3
  iatop(2,3) = 0
  iatop(2,4) = 1
  iatop(2,5) = 0
  iatop(2,6) = 15
 
  iatop(3,1) = -120
  iatop(3,2) = -100
  iatop(3,3) = 0
  iatop(3,4) = 0
  iatop(3,5) = 1
  iatop(3,6) = 70
 
  do i = 1, nrow
    do j = 1, ncol
      if ( iform == 0 ) then
        iabot(i,j) = 1
      else if ( iform == 2 ) then
        iabot(i,j) = 0
      end if
    end do
  end do
 
  a(1:nrow,1:ncol) = real ( iatop(1:nrow,1:ncol), kind = rk )
 
  ibase(1) = 3
  ibase(2) = 4
 
  chineq(1) = '<'
  chineq(2) = '<'
  chineq(3) = ' '
 
  imat = 1
 
  output = ' '
  call s_write ( iounit, output )
  output = 'Simple linear programming problem:'
  call s_write ( iounit, output )
  output = ' '
  call s_write ( iounit, output )
  output = 'Maximize:'
  call s_write ( iounit, output )
  output = '  Z = 120 X + 100 Y + 70'
  call s_write ( iounit, output )
  output = 'subject to'
  call s_write ( iounit, output )
  output = '  2 X + 2 Y < 8'
  call s_write ( iounit, output )
  output = '  5 X + 3 Y < 15'
  call s_write ( iounit, output )
  output = ' '
  call s_write ( iounit, output )
 
  return
end
subroutine lp_set ( ierror, imat, iounit, line, lpmoda, nart, ncol, ncon, &
  nrow, n_slack, nvar )

!*****************************************************************************80
!
!! lp_set() switches the linear programming mode.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IERROR.
!    0, no error occurred.
!    1, an error occurred.
!
!    Input, integer IMAT, 0/1, a matrix HAS NOT/HAS been defined.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Workspace, character ( len = 255 ) LINE, used to hold the user's input.
!
!    Input/output, integer LPMODA, 0/1, linear algebra/linear programming mode.
!
!    Input, integer NART, the number of artificial variables.
!
!    Input, integer NCOL, the number of columns in the matrix.
!
!    Output, integer NCON, the number of constraints.
!
!    Input, integer NROW, the number of rows in the matrix.
!
!    Output, integer N_SLACK, the number of slack variables.
!
!    Output, integer NVAR, the number of basic variables.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ierror
  integer imat
  integer iounit(4)
  character isay
  integer iterm
  logical s_eqi
  character ( len = 255 ) line
  integer lpmoda
  integer nart
  integer ncol
  integer ncon
  integer nrow
  integer n_slack
  integer nvar
  character ( len = 255 ) output
  character ( len = 255 ) prompt

  lpmoda = 1 - lpmoda
 
  if ( lpmoda == 0 ) then
    output = 'Switching to linear algebra mode.'
    call s_write ( iounit, output )
    return
  end if
 
  output = 'Switching to linear programming mode.'
  call s_write ( iounit, output )
 
  if ( imat == 0 ) then
    return
  end if
 
  prompt = '"Y" to use current matrix in linear programming.'
  line = ' '
  iterm = 0
  call s_read ( isay, line, prompt, iounit, ierror, iterm )
  if ( ierror /= 0 ) then
    return
  end if
 
  if ( .not. s_eqi ( isay, 'Y' ) ) then
    imat = 0
    nrow = 0
    ncol = 0
    return
  end if
 
  output = ' '
  call s_write ( iounit, output )

  line = ' '
  prompt = '# of slack variables, # of artificial variables.'
  call i_read ( n_slack, line, prompt, iounit, ierror )
  if ( ierror /= 0 ) then
    return
  end if
 
  call i_read ( nart, line, prompt, iounit, ierror )
  if ( ierror /= 0 ) then
    return
  end if
 
  nvar = ncol - 2 - nart - n_slack
  ncon = nrow - 1
 
  if ( nvar <= 0 ) then
    output = 'Values too large or too small!'
    call s_write ( iounit, output )
    return
  end if
 
  output = ' '
  call s_write ( iounit, output )
  output = 'Now please set the row labels (=basic variables)'
  call s_write ( iounit, output )
  output = 'using the "C" command, with I2 = 0.'
  call s_write ( iounit, output )
 
  return
end
subroutine lp_sol ( a, iatop, iabot, ibase, iform, isltop, islbot, maxcol, &
  maxrow, ncol, nrow, sol )

!*****************************************************************************80
!
!! lp_sol() determines the current linear programming solution.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A(MAXROW,MAXCOL), the current matrix.
!
!    Input, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
!    IATOP and IABOT represent the current rational or decimal matrix.
!
!    Input, integer IBASE(MAXROW), the basic variables.
!
!    Input, integer IFORM, 0/1/2, rational/real/decimal arithmetic.
!
!    Output, integer ISLTOP(MAXROW), ISLBOT(MAXROW), the fractional
!    or decimal representation of the linear programming solution.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
!    Input, integer NCOL, the number of columns in the matrix.
!
!    Input, integer NROW, the number of rows in the matrix.
!
!    Output, real ( kind = rk ) SOL(MAXROW), the real representation of the
!    linear programming solution.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxcol
  integer maxrow

  real ( kind = rk ) a(maxrow,maxcol)
  integer i
  integer iabot(maxrow,maxcol)
  integer iatop(maxrow,maxcol)
  integer ibase(maxrow)
  integer iform
  integer islbot(maxcol)
  integer isltop(maxcol)
  integer j
  integer jbase
  integer ncol
  integer nrow
  real ( kind = rk ) sol(maxcol)

  do i = 1, ncol - 2
 
    jbase = 0
 
    do j = 1, nrow-1
      if ( ibase(j) == i ) then
        jbase = j
      end if
    end do
 
    if ( jbase /= 0 ) then
 
      if ( iform == 0 ) then
        isltop(i) = iatop(jbase,ncol)
        islbot(i) = iabot(jbase,ncol)
      else if ( iform == 1 ) then
        sol(i) = a(jbase,ncol)
      else if ( iform == 2 ) then
        isltop(i) = iatop(jbase,ncol)
        islbot(i) = iabot(jbase,ncol)
      end if
 
    else
 
      if ( iform == 0 ) then
        isltop(i) = 0
        islbot(i) = 1
      else if ( iform == 1 ) then
        sol(i) = 0.0D+00
      else if ( iform == 2 ) then
        isltop(i) = 0
        islbot(i) = 0
      end if
 
    end if
 
  end do
 
  return
end
subroutine mat_copy ( a, b, iatop, iabot, ibtop, ibbot, ibase, ibaseb, &
  lpmoda, lpmodb, maxcol, maxrow, nart, nartb, ncol, ncolb, nrow, nrowb, &
  n_slack, n_slackb, nvar, nvarb )

!*****************************************************************************80
!
!! mat_copy() makes a copy of the current problem information.
!
!  Discussion:
!
!    The routine essentially carries out the following copy operations:
!
!    A       --> B
!    IATOP   --> IBTOP
!    IABOT   --> IBBOT
!    IBASE   --> IBASEB
!    LPMODA  --> LPMODB
!    NART    --> NARTB
!    NCOL    --> NCOLB
!    NROW    --> NROWB
!    N_SLACK --> NSLAKB
!    NVAR    --> NVARB
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A(MAXROW,MAXCOL), the current matrix.
!
!    Output, real ( kind = rk ) B(MAXROW,MAXCOL), a copy of the input value of A.
!
!    Input, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
!    IATOP and IABOT represent the rational or decimal matrix.
!
!    Output, integer IBBOT(MAXROW,MAXCOL), IBBOT(MAXROW,MAXCOL),
!    a copy of the input values of IATOP and IABOT.
!
!    Input, integer IBASE(MAXROW), the basic variables.
!
!    Output, integer IBASEB(MAXROW), a copy of the input value of IBASE.
!
!    Input, integer LPMODA, 0/1, linear algebra/linear programming mode.
!
!    Output, integer LPMODB, a copy of the input value of LPMODA.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
!    Input, integer NART, the number of artificial variables.
!
!    Output, integer NARTB, a copy of the input value of NART.
!
!    Input, integer NCOL, the number of columns in the matrix.
!
!    Output, integer NCOLB, a copy of the input value of NCOL.
!
!    Input, integer NROW, the number of rows in the matrix.
!
!    Output, integer NROWB, a copy of the input value of NROW.
!
!    Input, integer N_SLACK, the number of slack variables.
!
!    Output, integer N_SLACKB, a copy of the input value of NSLAK.
!
!    Input, integer NVAR, the number of basic variables.
!
!    Output, integer NVARB, a copy of the input value of NVAR.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxcol
  integer maxrow

  real ( kind = rk ) a(maxrow,maxcol)
  real ( kind = rk ) b(maxrow,maxcol)
  integer iabot(maxrow,maxcol)
  integer iatop(maxrow,maxcol)
  integer ibbot(maxrow,maxcol)
  integer ibtop(maxrow,maxcol)
  integer ibase(maxrow)
  integer ibaseb(maxrow)
  integer lpmoda
  integer lpmodb
  integer nart
  integer nartb
  integer ncol
  integer ncolb
  integer nrow
  integer nrowb
  integer n_slack
  integer n_slackb
  integer nvar
  integer nvarb

  b(1:maxrow,1:maxcol) = a(1:maxrow,1:maxcol)
  ibaseb(1:maxrow) = ibase(1:maxrow)
  ibtop(1:maxrow,1:maxcol) = iatop(1:maxrow,1:maxcol)
  ibbot(1:maxrow,1:maxcol) = iabot(1:maxrow,1:maxcol)
  lpmodb = lpmoda
  nartb = nart
  ncolb = ncol
  nrowb = nrow
  n_slackb = n_slack
  nvarb = nvar
 
  return
end
subroutine mat_print ( a, iabot, iatop, ibase, iform, iounit, ihi, ilo, &
  jhi, jlo, lpmoda, maxcol, maxrow, ncol, nrow, title )

!*****************************************************************************80
!
!! mat_print() prints out the matrix or table and solution.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A(MAXROW,MAXCOL), the current matrix.
!
!    Input, integer IABOT(MAXROW,MAXCOL), IATOP(MAXROW,MAXCOL).
!    IATOP and IABOT represent the current rational or decimal matrix.
!
!    Input, integer IBASE(MAXROW), the basic variables.
!
!    Input, integer IFORM, 0/1/2, rational/real/decimal arithmetic.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Input, integer LPMODA, 0/1, linear algebra/linear programming mode.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
!    Input, integer NCOL, the number of columns in the matrix.
!
!    Input, integer NROW, the number of rows in the matrix.
!
!    Input, character ( len = 255 ) TITLE, the title of the object to be printed.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxcol
  integer maxrow

  real ( kind = rk ) a(maxrow,maxcol)
  integer iabot(maxrow,maxcol)
  integer iatop(maxrow,maxcol)
  integer ibase(maxrow)
  integer iform
  integer ihi
  integer ilo
  integer iounit(4)
  integer jhi
  integer jlo
  integer lpmoda
  integer ncol
  integer nrow
  character ( len = 255 ) title

  if ( iform == 0 ) then
 
    call rat_print ( iatop, iabot, ibase, iounit, ihi, ilo, jhi, &
      jlo, lpmoda, maxcol, maxrow, ncol, nrow, title )
 
  else if ( iform == 1 ) then
 
    call r_print ( a, ibase, iounit, ihi, ilo, jhi, jlo, lpmoda, &
      maxcol, maxrow, ncol, nrow, title )
 
  else if ( iform == 2 ) then
 
    call dec_print ( iatop, iabot, ibase, iounit, ihi, ilo, jhi, &
      jlo, lpmoda, maxcol, maxrow, ncol, nrow, title )
 
  end if
 
  return
end
subroutine mat_restore ( a, c, iabot, iatop, ibase, ibasec, icbot, ictop, &
  ierror, imat, iounit, lpmoda, lpmodc, maxcol, maxrow, nart, nartc, ncol, &
  ncolc, nrow, nrowc, n_slack, n_slackc, nvar, nvarc )

!*****************************************************************************80
!
!! mat_restore() restores a matrix that was saved earlier.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) A(MAXROW,MAXCOL), the input value of C.
!
!    Input, real ( kind = rk ) C(MAXROW,MAXCOL), a saved matrix.
!
!    Output, integer IABOT(MAXROW,MAXCOL), IATOP(MAXROW,MAXCOL).
!    The input value of ICTOP, ICBOT.
!
!    Output, integer IBASE(MAXROW), the input value of IBASEC.
!
!    Input, integer IBASEC(MAXROW), a saved vector to keep track
!    of basic variables.
!
!    Input, integer ICBOT(MAXROW,MAXCOL), ICTOP(MAXROW,MAXCOL).
!    A saved matrix.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred.
!
!    Input, integer IMAT, 0/1, a matrix HAS NOT/HAS been defined.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Output, integer LPMODA, the input value of LPMODC.
!
!    Input, integer LPMODC, a saved linear programming switch.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
!    Output, integer NART, the input value of NARTC.
!
!    Input, integer NARTC, a saved number of artificial variables.
!
!    Output, integer NCOL, the input value of NCOLC.
!
!    Input, integer NCOLC, a saved number of columns.
!
!    Output, integer NROW, the input value of NROWC.
!
!    Input, integer NROWC, a saved number of rows.
!
!    Output, integer N_SLACK, the input value of NSLAKC.
!
!    Input, integer N_SLACKC, a saved number of slack variables.
!
!    Output, integer NVAR, the input value of NVARC.
!
!    Input, integer NVARC, a saved number of basic variables.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxcol
  integer maxrow

  real ( kind = rk ) a(maxrow,maxcol)
  real ( kind = rk ) c(maxrow,maxcol)
  integer iabot(maxrow,maxcol)
  integer iatop(maxrow,maxcol)
  integer ibase(maxrow)
  integer ibasec(maxrow)
  integer icbot(maxrow,maxcol)
  integer ictop(maxrow,maxcol)
  integer ierror
  integer imat
  integer iounit(4)
  integer lpmoda
  integer lpmodc
  integer lpmods
  integer nart
  integer nartc
  integer ncol
  integer ncolc
  integer nrow
  integer nrowc
  integer n_slack
  integer n_slackc
  integer nvar
  integer nvarc
  character ( len = 255 ) output

  if ( imat == 0 ) then
    ierror = 1
    output = 'You must set up a matrix with the "E" command'
    call s_write ( iounit, output )
    output = 'before using the "R" command to restore it!'
    call s_write ( iounit, output )
    return
  end if
!
!  Is there a saved matrix to restore?
!
  if ( ncolc == 0 ) then
    output = 'There is no saved matrix to restore!'
    call s_write ( iounit, output )
    ierror = 1
    return
  end if
!
!  Save a copy of the current linear programming mode.
!
  lpmods = lpmoda
!
!  Overwrite the current information by the old information.
!
  call mat_copy ( c, a, ictop, icbot, iatop, iabot, ibasec, ibase, &
    lpmodc, lpmoda, maxcol, maxrow, nartc, nart, ncolc, ncol, &
    nrowc, nrow, n_slackc, n_slack, nvarc, nvar )
 
  output = 'The saved matrix has been restored.'
  call s_write ( iounit, output )
!
!  Print a warning if linear programming mode has been switched.
!
  if ( lpmods /= lpmoda ) then
    output = 'Note: The linear programming mode has been switched.'
    call s_write ( iounit, output )
  end if
 
  return
end
subroutine mat_zero ( a, iabot, iatop, iform, maxcol, maxrow )

!*****************************************************************************80
!
!! mat_zero() initializes the matrix by zeroing it out.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) A(MAXROW,MAXCOL), the current matrix.
!
!    Output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
!    IATOP and IABOT represent the current rational or decimal
!    matrix.
!
!    Output, integer IFORM, 0/1/2, rational/real/decimal arithmetic.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxcol
  integer maxrow

  real ( kind = rk ) a(maxrow,maxcol)
  integer i
  integer iabot(maxrow,maxcol)
  integer iatop(maxrow,maxcol)
  integer iform
  integer j

  do i = 1, maxrow
    do j = 1, maxcol

      if ( iform == 0 ) then
        iatop(i,j) = 0
        iabot(i,j) = 1
      else if ( iform == 1 ) then
        a(i,j) = 0.0D+00
      else if ( iform == 2 ) then
        iatop(i,j) = 0
        iabot(i,j) = 1
      end if

    end do
  end do
 
  return
end
function npage ( )

!*****************************************************************************80
!
!! npage() determines whether it's time to pause before more printing.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer NPAGE.
!    The current number of pages completed, defined as the
!    number of lines printed, divided by the number of lines per page.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer output_line_count
  integer output_page_length
  integer npage

  output_page_length = 0
  output_line_count = 0
!
!  Get the page length.
!
  call i_data ( 'GET', 'OUTPUT_PAGE_LENGTH', output_page_length )
 
  if ( output_page_length <= 0 ) then
    npage = 0
    return
  end if
!
!  Get the current line number.
!
  call i_data ( 'GET', 'OUTPUT_LINE_COUNT', output_line_count )
 
  npage = output_line_count / output_page_length
  output_line_count = output_line_count - npage * output_page_length

  call i_data ( 'SET', 'OUTPUT_LINE_COUNT', output_line_count )
 
  return
end
subroutine normal_01_sample ( x )

!*****************************************************************************80
!
!! normal_01_sample() samples the standard Normal PDF.
!
!  Discussion:
!
!    The standard normal distribution has mean 0 and standard
!    deviation 1.
!
!  Method:
!
!    The Box-Muller method is used.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) X, a sample of the PDF.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, save :: iset = 0
  real ( kind = rk ), parameter :: PI = 3.14159265358979323846264338327950288419716939937510D+00
  real ( kind = rk ) v1
  real ( kind = rk ) v2
  real ( kind = rk ) x
  real ( kind = rk ), save :: xsave = 0.0D+00

  if ( iset == 0 ) then

    call random_number ( harvest = v1 )

    if ( v1 <= 0.0D+00 ) then
      write ( *, * ) ' '
      write ( *, * ) 'NORMAL_01_SAMPLE - Fatal error!'
      write ( *, * ) '  V1 <= 0.'
      write ( *, * ) '  V1 = ', v1
      stop
    end if

    call random_number ( harvest = v2 )

    if ( v2 <= 0.0D+00 ) then
      write ( *, * ) ' '
      write ( *, * ) 'NORMAL_01_SAMPLE - Fatal error!'
      write ( *, * ) '  V2 <= 0.'
      write ( *, * ) '  V2 = ', v2
      stop
    end if

    x = sqrt ( - 2.0D+00 * log ( v1 ) ) * cos ( 2.0D+00 * PI * v2 )

    xsave = sqrt ( - 2.0D+00 * log ( v1 ) ) * sin ( 2.0D+00 * PI * v2 )

    iset = 1

  else

    x = xsave
    iset = 0

  end if

  return
end
subroutine orth ( a, ierror, iounit, irow, icol, maxcol, maxrow, ncol, nrow )

!*****************************************************************************80
!
!! orth() uses an orthogonal transformation to zero out entry A(IROW,ICOL).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = rk ) A(MAXROW,MAXCOL), the current matrix.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Input, integer IROW, ICOL, the indices of the entry to be zeroed.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
!    Input, integer NCOL, the number of columns in the matrix.
!
!    Input, integer NROW, the number of rows in the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxcol
  integer maxrow

  real ( kind = rk ) a(maxrow,maxcol)
  real ( kind = rk ) bot
  character ( len = 10 ) chrtmp1
  character ( len = 10 ) chrtmp2
  real ( kind = rk ) cj
  integer icol
  integer ierror
  integer iounit(4)
  integer irow
  integer j
  integer len1
  integer len2
  integer ncol
  integer nrow
  character ( len = 255 ) output
  real ( kind = rk ) sj
  real ( kind = rk ) t1
  real ( kind = rk ) t2
!
!  We can't do a diagonal entry.
!
  if ( irow == icol ) then
    output = 'You cannot zero out a diagonal entry!'
    call s_write ( iounit, output )
    return
  end if
!
!  Refuse to continue if a row is out of bounds.
!
  if ( ( irow < 1 .or. nrow < irow ) .or. ( icol < 1 .or. ncol < icol ) ) then
    ierror = 1
    output = 'One of the indices is illegal!'
    call s_write ( iounit, output )
    return
  end if
!
!  A(I,J) should not already be zero.
!
  if ( a(irow,icol) == 0.0D+00 ) then
    output = 'A(I,J) is already zero!'
    call s_write ( iounit, output )
    return
  end if
!
!  Compute CJ and SJ.
!
!  Q = ( C -S )
!      ( S  C )
!
  bot = sqrt ( a(irow,icol) ** 2 + a(icol,icol) ** 2 )
  cj = a(icol,icol) / bot
  sj = a(irow,icol) / bot
  write ( *, * ) 'CJ = ', cj
  write ( *, * ) 'SJ = ', sj
!
!  Mathematics:
!    A = Q' * ( Q * A )
!  Computation:
!    A <- Q * A
!
  do j = 1, ncol
    t1 = a(irow,j)
    t2 = a(icol,j)
    a(irow,j) = cj * t1 - sj * t2
    a(icol,j) = sj * t1 + cj * t2
  end do
 
  call i_to_s_left ( irow, chrtmp1 )
  len1 = len_trim ( chrtmp1 )

  call i_to_s_left ( icol, chrtmp2 )
  len2 = len_trim ( chrtmp2 )

  output = '  ORO: A(' // chrtmp1(1:len1) // ',' // chrtmp2(1:len2) // ') = 0.0'

  call s_write ( iounit, output )
 
  return
end
subroutine orth_random ( lda, n, a )

!*****************************************************************************80
!
!! orth_random() returns a random orthogonal matrix.
!
!  Discussion:
!
!    The inverse of A is equal to A'.
!
!    A * A'  = A' * A = I.
!
!    Columns and rows of A have unit Euclidean norm.
!
!    Distinct pairs of columns of A are orthogonal.
!
!    Distinct pairs of rows of A are orthogonal.
!
!    The L2 vector norm of A*x = the L2 vector norm of x for any vector x.
!
!    The L2 matrix norm of A*B = the L2 matrix norm of B for any matrix B.
!
!    The determinant of A is +1 or -1.
!
!    All the eigenvalues of A have modulus 1.
!
!    All singular values of A are 1.
!
!    All entries of A are between -1 and 1.
!
!    Thanks to Eugene Petrov, B I Stepanov Institute of Physics,
!    National Academy of Sciences of Belarus, for convincingly
!    pointing out the severe deficiencies of an earlier version of
!    this routine.
!
!    Essentially, the computation involves saving the Q factor of the
!    QR factorization of a matrix whose entries are normally distributed.
!    However, it is only necessary to generate this matrix a column at
!    a time, since it can be shown that when it comes time to annihilate
!    the subdiagonal elements of column K, these (transformed) elements of
!    column K are still normally distributed random values.  Hence, there
!    is no need to generate them at the beginning of the process and
!    transform them K-1 times.
!
!    For computational efficiency, the individual Householder transformations
!    could be saved, as recommended in the reference, instead of being
!    accumulated into an explicit matrix format.
!
!  Reference:
!
!    G W Stewart,
!    Efficient Generation of Random Orthogonal Matrices With an Application
!    to Condition Estimators,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 3, June 1980, pages 403-409.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of A.
!
!    Input, integer N, the order of A.
!
!    Output, real ( kind = rk ) A(LDA,N), the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer lda
  integer n

  real ( kind = rk ) a(lda,n)
  integer i
  integer j
  real ( kind = rk ) v(n)
  real ( kind = rk ) x(n)
!
!  Start with A = the identity matrix.
!
  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        a(i,j) = 1.0D+00
      else
        a(i,j) = 0.0D+00
      end if
    end do
  end do
!
!  Now behave as though we were computing the QR factorization of
!  some other random matrix.  Generate the N elements of the first column,
!  compute the Householder matrix H1 that annihilates the subdiagonal elements,
!  and set A := A * H1' = A * H.
!
!  On the second step, generate the lower N-1 elements of the second column,
!  compute the Householder matrix H2 that annihilates them,
!  and set A := A * H2' = A * H2 = H1 * H2.
!
!  On the N-1 step, generate the lower 2 elements of column N-1,
!  compute the Householder matrix HN-1 that annihilates them, and
!  and set A := A * H(N-1)' = A * H(N-1) = H1 * H2 * ... * H(N-1).
!  This is our random orthogonal matrix.
!
  do j = 1, n - 1
!
!  Set the vector that represents the J-th column to be annihilated.
!
    do i = 1, j - 1
      x(i) = 0.0D+00
    end do
    do i = j, n
      call normal_01_sample ( x(i) )
    end do
!
!  Compute the vector V that defines a Householder transformation matrix
!  H(V) that annihilates the subdiagonal elements of X.
!
    call rvec_house_column ( n, x, j, v )
!
!  Postmultiply the matrix A by H'(V) = H(V).
!
    call rmat_house_axh ( lda, n, a, v, a )

  end do

  return
end
subroutine r_jacobi ( a, ibase, ierror, iounit, line, lpmoda, maxcol, maxrow, &
  ncol, nrow, ev, q )

!*****************************************************************************80
!
!! r_jacobi() carries out a Jacobi rotation on a real square matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = rk ) A(MAXROW,MAXCOL), the current matrix.
!
!    Input, integer IBASE(MAXROW), the basic variables.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred.
!
!    Input, integer IFORM, 0/1/2, rational/real/decimal arithmetic.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Workspace, character ( len = 255 ) LINE, used to hold the user's input.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
!    Input, integer NCOL, the number of columns in the matrix.
!
!    Input, integer NROW, the number of rows in the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxcol
  integer maxrow

  real ( kind = rk ) a(maxrow,maxcol)
  real ( kind = rk ) cj
  real ( kind = rk ) ev(maxrow,maxcol)
  integer i
  integer ibase(maxrow)
  integer ierror
  integer ihi
  integer ilo
  integer iounit(4)
  integer irow
  integer j
  integer jcol
  integer jhi
  integer jlo
  integer k
  character ( len = 255 ) line
  integer lpmoda
  integer n
  integer ncol
  integer nrow
  character ( len = 255 ) output
  character ( len = 255 ) prompt
  real ( kind = rk ) q(maxrow,maxcol)
  logical s_eqi
  real ( kind = rk ) sj
  logical sym
  real ( kind = rk ) t1
  real ( kind = rk ) t2
  real ( kind = rk ) temp
  character ( len = 255 ) title
  real ( kind = rk ) tj
  real ( kind = rk ) u
!
!  Return if matrix is not square.
!
  if ( nrow /= ncol ) then
    output = 'Jacobi iteration requires a square matrix!'
    ierror = 1
    call s_write ( iounit, output )
    return
  end if

  n = nrow
!
!  Test for symmetry.
!
  sym = .true.
  do i = 1, n
    do j = 1, i-1
      if ( a(i,j) /= a(j,i) ) then
        sym = .false.
      end if
    end do
  end do
 
  if ( .not. sym ) then
    output = 'Warning!  Because the matrix is not symmetric,'
    call s_write ( iounit, output )
    output = 'Jacobi''s method may not converge!'
    call s_write ( iounit, output )
  end if
!
!  Set up the eigenvector data.
!
  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        ev(i,j) = 1.0D+00
      else
        ev(i,j) = 0.0D+00
      end if
    end do
  end do
!
!  Here is where repetition will begin.
!
30    continue
!
!  Print the current matrix.
!
  ilo = 1
  ihi = n
  jlo = 1
  jhi = n
  title = 'The current matrix'

  call r_print ( a, ibase, iounit, ihi, ilo, jhi, jlo, lpmoda, maxcol, &
    maxrow, ncol, nrow, title )
!
!  Get the row of the entry.
!
40    continue

  output = ' '
  call s_write ( iounit, output )
  prompt = 'row I, column J, or "Q" to quit.'
  call i_read ( irow, line, prompt, iounit, ierror )
 
  if ( ierror /= 0 ) then
    if ( s_eqi ( line(1:1), 'Q' ) ) then
      line = ' '
      ierror = 0
    end if
    return
  end if
 
  if ( irow <= 0 .or. nrow < irow ) then
    output = 'The value of I, the row index, is illegal.'
    call s_write ( iounit, output )
    go to 40
  end if
!
!  Get the column of the entry.
!
50    continue

  prompt = 'column J or "Q" to quit.'
  call i_read ( jcol, line, prompt, iounit, ierror )
 
  if ( ierror /= 0 ) then
    if ( s_eqi ( line(1:1), 'Q' ) ) then
      line = ' '
      ierror = 0
    end if
    return
  end if
 
  if ( jcol <= 0 .or. ncol < jcol ) then
    output = 'The value of J, the column index, is illegal.'
    call s_write ( iounit, output )
    go to 50
  end if
!
!  I and J must not be equal.
!
  if ( irow == jcol ) then
    output = 'Jacobi rotations require I and J to be distinct!'
    call s_write ( iounit, output )
    go to 40
  end if
!
!  A(I,J) should not already be zero.
!
  if ( a(irow,jcol) == 0.0D+00 ) then
    output = 'A(I,J) is already zero!'
    call s_write ( iounit, output )
    go to 40
  end if
!
!  If the matrix is nonsymmetric, we require that A(I,J)+A(J,I)
!  not be zero.
!
  if ( a(irow,jcol) + a(jcol,irow) == 0.0D+00 ) then
    output = ' '
    call s_write ( iounit, output )
    output = 'The algorithm breaks down for this (I,J)'
    call s_write ( iounit, output )
    output = 'because A(I,J) + A(J,I) is zero!'
    call s_write ( iounit, output )
    go to 40
  end if
!
!  Compute CJ and SJ.
!
  u = ( a(jcol,jcol) - a(irow,irow) ) / ( a(irow,jcol) + a(jcol,irow) )
 
  if ( 0.0D+00 <= u ) then
    temp = 1.0D+00
  else
    temp = - 1.0D+00
  end if
 
  tj = temp / ( abs ( u ) + sqrt ( u ** 2 + 1.0D+00 ) )
  cj = 1.0D+00 / sqrt ( tj ** 2 + 1.0D+00 )
  sj = tj * cj
!
!  Set up the Q matrix.
!
  do i = 1, n
    do j = 1, n

      if      ( i == irow .and. j == jcol ) then
        q(i,j) = sj
      else if ( i == irow .and. j == irow ) then
        q(i,j) = cj
      else if ( i == jcol .and. j == irow ) then
        q(i,j) = - sj
      else if ( i == jcol .and. j == jcol ) then
        q(i,j) = cj
      else if ( i == j ) then
        q(i,j) = 1.0D+00
      else
        q(i,j) = 0.0D+00
      end if
      
    end do
  end do
!
!  Print the Q matrix.
!
  ilo = 1
  ihi = n
  jlo = 1
  jhi = n
  title = 'The single step Q transformation matrix'

  call r_print ( q, ibase, iounit, ihi, ilo, jhi, jlo, lpmoda, &
    maxcol, maxrow, n, n, title )
!
!  A -> A * Q
!
  do k = 1, n
    t1 = a(irow,k)
    t2 = a(jcol,k)
    a(irow,k) = t1 * cj - t2 * sj
    a(jcol,k) = t1 * sj + t2 * cj
  end do
!
!  A -> Q' * A
!
  do k = 1, n
    t1 = a(k,irow)
    t2 = a(k,jcol)
    a(k,irow) = cj * t1 - sj * t2
    a(k,jcol) = sj * t1 + cj * t2
  end do
!
!  Premultiply the eigenvectors by Q'.
!    EV -> Q' * EV.
!
  do k = 1, n
    t1 = ev(k,irow)
    t2 = ev(k,jcol)
    ev(k,irow) =    cj * t1 - sj * t2
    ev(k,jcol) =    sj * t1 + cj * t2
  end do
!
!  Print the EV matrix.
!
  ilo = 1
  ihi = n
  jlo = 1
  jhi = n
  title = 'The approximate eigenvector matrix (Q-Transpose)'

  call r_print ( ev, ibase, iounit, ihi, ilo, jhi, jlo, lpmoda, &
    maxcol, maxrow, n, n, title )
 
  go to 30
 
end
subroutine r_print ( a, ibase, iounit, ihi, ilo, jhi, jlo, lpmoda, maxcol, &
  maxrow, ncol, nrow, title )

!*****************************************************************************80
!
!! r_print() prints out real vectors and matrices.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A(MAXROW,MAXCOL), the current matrix.
!
!    Input, integer BASE(MAXROW), keeps track of basic variables.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Input, integer IHI, ILO, the last and first rows to print.
!
!    Input, integer JHI, JLO, the last and first columns to print.
!
!    Input, integer LPMODA, 0/1, linear algebra/linear programming mode.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
!    Input, integer NCOL, the number of columns in the matrix.
!
!    Input, integer NROW, the number of rows in the matrix.
!
!    Input, character ( len = * ) TITLE, a label for the object being printed.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: ncolum = 80

  integer maxcol
  integer maxrow

  real ( kind = rk ) a(maxrow,maxcol)
  logical allint
  character ( len = 10 ) chrtmp1
  character ( len = 10 ) chrtmp2
  character ( len = 10 ) chrtmp3
  character ( len = 40 ) format1
  character ( len = 40 ) format2
  integer i
  integer ibase(maxrow)
  integer ichi
  integer iclo
  integer ihi
  integer ilo
  integer imax
  integer imin
  integer iounit(4)
  integer izhi
  integer izlo
  integer j
  integer jhi
  integer jlo
  integer jmax
  integer jmin
  integer kmax
  integer kmin
  character ( len = 4 ) lab
  integer llab
  integer lpmoda
  integer ncol
  integer npline
  integer nrow
  character ( len = 255 ) output
  character ( len = * ) title

  if ( lpmoda == 1 ) then
    llab = 4
  else
    llab = 0
  end if
!
!  Figure out how many numbers we can fit in (NCOLUM-LLAB) columns.
!
  kmin = 11
  kmax = kmin
 
  do i = ilo, ihi
    do j = jlo, jhi
  
      do while ( 10.0D+00 ** ( kmax - kmin ) <= abs ( a(i,j) ) )
        kmax = kmax + 1
      end do
 
    end do
  end do
 
  npline = ( ncolum - llab ) / kmax
!
!  Check to see if the matrix entries are all integers.
!
  allint = .true.
 
  do i = ilo, ihi
    do j = jlo, jhi
 
      if ( a(i,j) /= real ( int ( a(i,j) ), kind = rk ) ) then
        allint = .false.
      end if

    end do
  end do
!
!  If all integers, cut down KMAX, the width of each number,
!  and update NPLINE, the number of numbers we can print on one line.
!
  if ( allint ) then

    kmax = kmax - 7
    npline = ( ncolum - llab ) / kmax
 
    call i_to_s_left ( llab, chrtmp1 )
    call i_to_s_left ( npline, chrtmp2 )
    call i_to_s_left ( kmax, chrtmp3 )

    if ( lpmoda == 1 ) then
      format1 = '(a' // chrtmp1 // ',' // chrtmp2 // 'f' // chrtmp3 // '.0)'
    else
      format1 = '(' // chrtmp2 // 'f' // chrtmp3 // '.0)'
    end if
!
!  If nonintegral entries, print 7 decimals.
!
  else
 
    call i_to_s_left ( llab, chrtmp1 )
    call i_to_s_left ( npline, chrtmp2 )
    call i_to_s_left ( kmax, chrtmp3 )

    if ( lpmoda == 1 ) then
      format1 = '(a' // chrtmp1 // ',' // chrtmp2 // 'f' // chrtmp3 // '.7)'
    else
      format1 = '(' // chrtmp2 // 'f' // chrtmp3 // '.7)'
    end if
 
  end if
 
  call s_blank_delete ( format1 )
!
!  The second format is for ...
!
  if ( lpmoda == 1 ) then
    format2 = '(' // chrtmp1 // 'x,' // chrtmp2 // 'i' // chrtmp3 // ')'
  else
    format2 = '(' // chrtmp2 // 'i' // chrtmp3 // ')'
  end if

  call s_blank_delete ( format2 )
!
!  Now print the data.
!
  do jmin = jlo, jhi, npline
 
    jmax = min ( jmin + npline - 1, jhi )
    lab = '    '
!
!  Handle a column vector.
!
    if ( jlo == jhi .and. ilo /= ihi ) then
 
      output = ' '
      call s_write ( iounit, output )
 
      if ( ilo == 1 ) then
        output = title
        call s_write ( iounit, output )
        call i_to_s_left ( jlo, chrtmp1 )
        output = 'Column ' // chrtmp1 // ' transposed.'
        call s_blanks_delete ( output )
        call s_write ( iounit, output )
      end if
 
      do imin = ilo, ihi, npline
 
        imax = min ( imin + npline - 1, ihi )
 
        output = ' '
        call s_write ( iounit, output )
 
        if ( lpmoda == 0 ) then
          write ( output, format1 ) ( a(i,jlo), i = imin, imax )
        else
          write ( output, format1 ) lab, ( a(i,jlo), i = imin, imax )
        end if

        call s_write ( iounit, output )

      end do
!
!  Print a matrix.
!
    else
 
      output = ' '
      call s_write ( iounit, output )
 
      if ( jmin == 1 ) then
        output = title
        call s_write ( iounit, output )
        output = ' '
        call s_write ( iounit, output )
      end if
!
!  Print heading for linear programming table.
!
      if ( lpmoda == 1 ) then
 
        write ( output, format2 ) ( j, j = jmin, jmax )
 
        if ( jmin <= ncol-1 .and. ncol-1 <= jmax ) then
          izlo = llab + ( ( ncol - 1 ) - jmin ) * kmax + kmax - 2
          izhi = izlo + 2
          output(izlo:izhi) = '  P'
        end if
 
        if ( jmin <= ncol .and. ncol <= jmax ) then
          iclo = llab + ( ncol - jmin ) * kmax + kmax - 2
          ichi = iclo + 2
          output(iclo:ichi) = '  C'
        end if
 
        call s_write ( iounit, output )
  
        output = ' '
        call s_write ( iounit, output )
!
!  Or print heading for linear algebra matrix.
!
      else

        if ( 1 < jmin .or. jmax < ncol .or. 1 < ilo .or. ihi < nrow ) then
 
          call i_to_s_left ( jmin, chrtmp1 )
          call i_to_s_left ( jmax, chrtmp2 )
          output = 'Columns ' // chrtmp1 // ' to ' // chrtmp2
          call s_blanks_delete ( output )
          call s_write ( iounit, output )
          output = ' '
          call s_write ( iounit, output )
        end if
      end if
 
      do i = ilo, ihi
 
        if ( lpmoda == 1 ) then

          if ( i < nrow ) then
            if ( ibase(i) < 10 ) then
              write ( lab, '(a1,i1)' ) 'X', ibase(i)
            else
              write ( lab, '(a1,i2)' ) 'X', ibase(i)
            end if
          else if ( i < ihi ) then
            lab = 'Obj2'
          else
            lab = 'Obj '
          end if

          if ( maxrow == 1 )then
            lab = '    '
          end if

        end if
 
        if ( lpmoda == 1 ) then
          write ( output, format1 ) lab, ( a(i,j), j = jmin, jmax )
        else
          write ( output, format1 ) ( a(i,j), j = jmin, jmax )
        end if
        call s_write ( iounit, output )
 
      end do
 
    end if
 
  end do
 
  return
end
subroutine r_read ( rval, line, prompt, iounit, ierror )

!*****************************************************************************80
!
!! r_read() "reads" a real value from a line of text.
!
!  Discussion:
!
!    The routine accepts a line of characters which may contain some
!    user input.  If not, it prints out the PROMPT and reads new
!    information into LINE, seeking to find a real number RVAL to
!    return.
!
!    The routine will accept integers, decimals, and ratios of the
!    form R1/R2.  Real numbers may be in scientific notation, as
!    +12.34E-56.78
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) RVAL, the real value found in LINE.
!
!    Workspace, character ( len = 255 ) LINE, used to hold the user's input.
!
!    Workspace, character ( len = 255 ) PROMPT, the prompt string.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Output, integer IERROR.
!    0, no error occurred.
!    nonzero, an error occurred while trying to read the value.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) bot
  integer ierror
  integer iounit(4)
  integer lchar
  character ( len = 255 ) line
  character ( len = 255 ) prompt
  real ( kind = rk ) rval
  real ( kind = rk ) top

  ierror = 0
  rval = 0.0D+00
  top = 0.0D+00
  bot = 1.0D+00
!
!  Read a character string.
!
10    continue
 
  call chrinp ( ierror, iounit, line, prompt )
  if ( ierror /= 0 ) then
    return
  end if

  if ( line == ' ' ) then
    go to 10
  end if
!
!  Convert the character string to a decimal value, TOP.
!
  call s_to_r ( line, top, ierror, lchar )
!
!  If we haven't used up all our characters,
!  and if the next character is '/',
!  then the user means to input the value as a ratio,
!  so prepare to read BOT as well.
!
  if ( lchar+1 < len ( line ) ) then
    if ( line(lchar+1:lchar+1) == '/' ) then
      lchar = lchar + 1
      call s_chop ( line, 1, lchar )
      call s_to_r ( line, bot, ierror, lchar )
      if ( bot == 0.0D+00 ) then
        bot = 1.0D+00
      end if
    end if
  end if
!
!  Set the value of RVAL.
!
  rval = top / bot
!
!  Chop out the characters that were used.
!
  if ( lchar < len ( line ) ) then
    if ( line(lchar+1:lchar+1) == ',' ) then
      lchar = lchar + 1
    end if
  end if

  call s_chop ( line, 1, lchar )
 
  return
end
subroutine r_swap ( x, y )

!*****************************************************************************80
!
!! r_swap() switches two real values.
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
!    Input/output, real ( kind = rk ) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) z
  real ( kind = rk ) x
  real ( kind = rk ) y

  z = x
  x = y
  y = z

  return
end
subroutine r_to_dec ( rval, itop, ibot )

!*****************************************************************************80
!
!! r_to_dec() converts a real value to a decimal fraction form.
!
!  Discussion:
!
!    The routine is given RVAL, and computes ITOP and IBOT, so that 
!    approximately:
!
!      RVAL = ITOP * 10 ^ IBOT
!
!    However, only DEC_DIGIT digits of RVAL are used in constructing the 
!    representation, where DEC_DIGIT is the maximum number of decimal
!    digits used in the decimal representation.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) RVAL, the real number to be converted.
!
!    Output, integer ITOP, IBOT, the approximate decimal representation.
!    ITOP is an integer, strictly between -10 ^ DEC_DIGIT and 10 ^ DEC_DIGIT.
!    IBOT is an integer exponent of 10.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer dec_digit
  integer ibot
  integer itop
  real ( kind = rk ) rtop
  real ( kind = rk ) rval
  real ( kind = rk ) ten1
  real ( kind = rk ) ten2
!
!  Special cases.
!
  if ( rval == 0.0D+00 ) then
    itop = 0
    ibot = 0
    return
  end if

  dec_digit = 0
  call i_data ( 'GET', 'DEC_DIGIT', dec_digit )
!
!  Factor RVAL = RTOP * 10 ^ IBOT
!
  rtop = rval
  ibot = 0
!
!  Now normalize so that 10 ^ (DEC_DIGIT-1) <= ABS(RTOP) < 10 ^ (DEC_DIGIT)
!
  ten1 = 10.0D+00 ** ( dec_digit - 1 )
  ten2 = 10.0D+00 **   dec_digit
  
  do while ( abs ( rtop ) < ten1 )
    rtop = rtop * 10.0D+00
    ibot = ibot - 1
  end do

  do while ( ten2 <= abs ( rtop ) )
    rtop = rtop / 10.0D+00
    ibot = ibot + 1
  end do
!
!  ITOP is the integer part of RTOP, rounded.
!
  itop = nint ( rtop )
!
!  Now divide out any factors of ten from ITOP.
!
  if ( itop /= 0 ) then

    do while ( mod ( itop, 10 ) == 0 )
      itop = itop / 10
      ibot = ibot + 1
    end do

  end if
 
  return
end
subroutine r_to_rat ( a, iatop, iabot )

!*****************************************************************************80
!
!! r_to_rat() converts a real value to a rational value.  
!
!  Discussion:
!
!    The rational value (IATOP/IABOT) is essentially computed by truncating 
!    the decimal representation of the real value after a given number of
!    decimal digits.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, the real value to be converted.
!
!    Output, integer IATOP, IABOT, the numerator and denominator
!    of the rational value that approximates A.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  integer dec_digit
  real ( kind = rk ) factor
  integer i_gcd
  integer iabot
  integer iatop
  integer ibot
  integer ifac
  integer itemp
  integer itop
  integer jfac

  dec_digit = 0
  call i_data ( 'GET', 'DEC_DIGIT', dec_digit )

  factor = 10.0D+00 ** dec_digit
 
  if ( 0 < dec_digit ) then
    ifac = 10 ** dec_digit
    jfac = 1
  else
    ifac = 1
    jfac = 10 ** ( - dec_digit )
  end if
 
  itop = nint ( a * factor ) * jfac
  ibot = ifac
!
!  Factor out the greatest common factor.
!
  itemp = i_gcd ( itop, ibot )
 
  iatop = itop / itemp
  iabot = ibot / itemp
 
  return
end
subroutine r_to_s_left ( rval, s )

!*****************************************************************************80
!
!! r_to_s_left() represents a real using 14 left_justified characters.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) RVAL, a real number.
!
!    Output, character ( len = * ) S, a left-justified character variable 
!    containing the representation of RVAL.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = 14 ) chrtmp
  integer i
  real ( kind = rk ) rval
  character ( len = * ) s
!
!  We can't seem to write directly into the string because of compiler
!  quibbles.
!
  if ( real ( int ( rval ), kind = rk ) == rval .and. abs ( rval ) < 1.0D+13 ) then
 
    write ( chrtmp, '(i14)' ) int ( rval )
 
  else
 
    write ( chrtmp, '(g14.6)' ) rval
 
  end if
 
  do i = 1, len ( chrtmp )
    if ( chrtmp(i:i) /= ' ' ) then
      s = chrtmp(i:)
      return
    end if
  end do

  s = ' '

  return
end
subroutine rat_add ( ibot, ibot1, ibot2, itop, itop1, itop2 )

!*****************************************************************************80
!
!! rat_add() adds two rational values.
!
!  Discussion:
!
!    If numeric overflow would occur, the computation is done in
!    real arithmetic and then converted back to a rational value.
!
!    ITOP / IBOT = ( ITOP1 / IBOT1 ) + ( ITOP2 / IBOT2 )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IBOT, the denominator of the result.
!
!    Input, integer IBOT1, IBOT2, the denominators of the
!    two rational values to be added.
!
!    Output, integer ITOP, the numerator of the result.
!
!    Input, integer ITOP1, ITOP2, the numerators of the
!    two rational values to be added.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i_big
  integer i_gcd
  integer ibot
  integer ibot1
  integer ibot2
  integer ierror
  integer itemp
  integer itop
  integer itop1
  integer itop2
  integer jbot1
  integer jbot2
  integer jbot3
  integer jtop1
  integer jtop2
  real ( kind = rk ) rbot
  real ( kind = rk ) rmax
  real ( kind = rk ) rtop1
  real ( kind = rk ) rtop2
  real ( kind = rk ) rtop3
  real ( kind = rk ) rval
  real ( kind = rk ) rval1
  real ( kind = rk ) rval2

  ierror = 0
 
  i_big = 0
  call i_data ( 'GET', 'I_BIG', i_big )

  rmax = real ( i_big, kind = rk )

  if ( itop1 == 0 ) then
    itop = itop2
    ibot = ibot2
    return
  else if ( itop2 == 0 ) then
    itop = itop1
    ibot = ibot1
    return
  end if
!
!  Make copies of the input arguments, since we will change them.
!
  jbot1 = ibot1
  jbot2 = ibot2
  jtop1 = itop1
  jtop2 = itop2
!
!  Compute the greatest common factor of the two denominators,
!  and factor it out.
!
  jbot3 = i_gcd ( jbot1, jbot2 )
  jbot1 = jbot1 / jbot3
  jbot2 = jbot2 / jbot3
!
!  The fraction may now be formally written as:
!
!    (jtop1*jbot2 + jtop2*jbot1) / (jbot1*jbot2*jbot3)
!
!  Check the tops for overflow.
!
  rtop1 = real ( jtop1, kind = rk ) * real ( jbot2, kind = rk )
 
  if ( rmax < abs ( rtop1 ) ) then
    ierror = 1
    itop = 0
  else
    jtop1 = jtop1 * jbot2
  end if
 
  rtop2 = real ( jtop2, kind = rk ) * real ( jbot1, kind = rk )
 
  if ( rmax < abs ( rtop2 ) ) then
    ierror = 2
    itop = 0
  else
    jtop2 = jtop2 * jbot1
  end if
 
  rtop3 = real ( jtop1, kind = rk ) + real ( jtop2, kind = rk )
 
  if ( rmax < abs ( rtop3 ) ) then
    ierror = 3
    itop = 0
  else
    itop = jtop1 + jtop2
  end if
!
!  Check the bottom for overflow.
!
  rbot = real ( jbot1, kind = rk ) &
       * real ( jbot2, kind = rk ) &
       * real ( jbot3, kind = rk )
 
  if ( rmax < abs ( rbot ) ) then
    ierror = 4
    ibot = 1
  else
    ibot = jbot1 * jbot2 * jbot3
  end if
!
!  If there was potential overflow, then do the computation in
!  real arithmetic and convert back.
!
  if ( ierror /= 0 ) then
    ierror = 0
    rval1 = real ( itop1, kind = rk ) / real ( ibot1, kind = rk )
    rval2 = real ( itop2, kind = rk ) / real ( ibot2, kind = rk )
    rval = rval1 + rval2
    call r_to_rat ( rval, itop, ibot )
  end if
!
!  Put the fraction in lowest terms.
!
  itemp = i_gcd ( itop, ibot )
  itop = itop / itemp
  ibot = ibot / itemp
!
!  Sign of bottom should be positive.
!
  if ( ibot < 0 ) then
    ibot = - ibot
    itop = - itop
  end if
 
  return
end
subroutine rat_mul ( ibot, ibot1, ibot2, itop, itop1, itop2, ierror )

!*****************************************************************************80
!
!! rat_mul() multiplies two fractions.
!
!  Discussion:
!
!    If numeric overflow would occur, the computation is done in
!    real arithmetic and then converted back to a rational value.
!
!    ITOP / IBOT = ( ITOP1 / IBOT1 ) * ( ITOP2 / IBOT2 ).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IBOT, the denominator of the result.
!
!    Input, integer IBOT1, IBOT2, the denominators of the
!    two rational values to be multiplied.
!
!    Output, integer ITOP, the numerator of the result.
!
!    Input, integer ITOP1, ITOP2, the numerators of the
!    two rational values to be multiplied.
!
!    Output, integer IERROR, 0 for no error, 1 if either denominator
!    is zero.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i_big
  integer i_gcd
  integer ibot
  integer ibot1
  integer ibot2
  integer ierror
  integer itemp
  integer itop
  integer itop1
  integer itop2
  integer jbot1
  integer jbot2
  integer jtop1
  integer jtop2
  real ( kind = rk ) rbot
  real ( kind = rk ) rmax
  real ( kind = rk ) rtop
  real ( kind = rk ) temp

  ierror = 0
 
  i_big = 0
  call i_data ( 'GET', 'I_BIG', i_big )

  rmax = real ( i_big, kind = rk )

  if ( itop1 == 0 .or. itop2 == 0 ) then
    itop = 0
    ibot = 1
    return
  end if

  if ( ibot1 == 0 .or. ibot2 == 0 ) then
    ierror = 1
    itop = 0
    ibot = 1
    write ( *, * ) ' '
    write ( *, * ) 'RAT_MUL - Fatal error!'
    write ( *, * ) '  A rational fraction has a zero denominator!'
    return
  end if
!
!  Make copies of the input arguments, since we will change them.
!
  jbot1 = ibot1
  jbot2 = ibot2
  jtop1 = itop1
  jtop2 = itop2
!
!  Get rid of all common factors in top and bottom.
!
  itemp = i_gcd ( jtop1, jbot1 )
  jtop1 = jtop1 / itemp
  jbot1 = jbot1 / itemp
  itemp = i_gcd ( jtop1, jbot2 )
  jtop1 = jtop1 / itemp
  jbot2 = jbot2 / itemp
  itemp = i_gcd ( jtop2, jbot1 )
  jtop2 = jtop2 / itemp
  jbot1 = jbot1 / itemp
  itemp = i_gcd ( jtop2, jbot2 )
  jtop2 = jtop2 / itemp
  jbot2 = jbot2 / itemp
!
!  The fraction (ITOP1*ITOP2)/(IBOT1*IBOT2) is in lowest terms.
!
!  Check the top ITOP1*ITOP2 for overflow.
!
  rtop = real ( jtop1, kind = rk ) * real ( jtop2, kind = rk )
 
  if ( rmax < abs ( rtop ) ) then
    ierror = 1
  else
    itop = jtop1 * jtop2
  end if
!
!  Check the bottom IBOT1*IBOT2 for overflow.
!
  rbot = real ( jbot1, kind = rk ) * real ( jbot2, kind = rk )
 
  if ( rmax < abs ( rbot ) ) then
    ierror = 2
  else
    ibot = jbot1 * jbot2
  end if
!
!  if there was an overflow, then compute RTOP / RBOT and convert
!  back to rational.
!
  if ( ierror == 1 .or. ierror == 2 ) then
    ierror = 0
    temp = rtop / rbot
    call r_to_rat ( temp, itop, ibot )
  end if
!
!  Sign of bottom should be positive.
!
  if ( ibot < 0 ) then
    ibot = - ibot
    itop = - itop
  end if
!
!  The fraction is ITOP/IBOT with no loss of accuracy (unless there
!  was overflow).
!
  return
end
subroutine rat_print ( iatop, iabot, ibase, iounit, ihi, ilo, jhi, &
  jlo, lpmoda, maxcol, maxrow, ncol, nrow, title )

!*****************************************************************************80
!
!! rat_print() prints out rational vectors or matrices.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
!    IATOP and IABOT represent the current rational or decimal
!    matrix.
!
!    Input, integer IBASE(MAXROW), the basic variables.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Input, integer IHI, ILO, the last and first rows to print.
!
!    Input, integer JHI, JLO, the last and first columns to print.
!
!    Input, integer LPMODA, 0/1, linear algebra/linear programming mode.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
!    Input, integer NCOL, the number of columns in the matrix.
!
!    Input, integer NROW, the number of rows in the matrix.
!
!    Input, character ( len = * ) TITLE, a label for the object being printed.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: ncolum = 80

  integer maxcol
  integer maxrow

  character ( len = 10 ) chrtmp1
  character ( len = 10 ) chrtmp2
  character ( len = 10 ) chrtmp3
  character ( len = 40 ) format1
  character ( len = 40 ) format2
  integer i
  integer iabot(maxrow,maxcol)
  integer iatop(maxrow,maxcol)
  integer ibase(maxrow)
  integer ichi
  integer iclo
  integer ihi
  integer ilo
  integer imax
  integer imin
  integer ione
  integer iounit(4)
  integer itemp
  integer izhi
  integer izlo
  integer j
  integer jhi
  integer jlo
  integer jmax
  integer jmin
  integer kmax
  character ( len = 4 ) lab
  integer llab
  integer lpmoda
  integer ncol
  integer none
  integer npline
  integer nrow
  character ( len = 255 ) output
  character ( len = * ) title

  if ( lpmoda == 1 ) then
    llab = 4
  else
    llab = 0
  end if
!
!  Figure out how many rationals we can get in (NCOLUM-LLAB) columns.
!
  lab = '    '
  kmax = 3
 
  do i = ilo, ihi
    do j = jlo, jhi
 
      itemp = abs ( iatop(i,j) )
 
      do while ( 10 ** ( kmax - 2 ) <= itemp )
        kmax = kmax + 1
      end do
 
      itemp = abs ( iabot(i,j) )
 
      do while ( 10 ** ( kmax - 2 ) < itemp )
        kmax = kmax + 1
      end do
 
    end do
  end do
 
  kmax = kmax + 1
  npline = ( ncolum - llab ) / kmax
!
!  Create the formats.
!
  call i_to_s_left ( llab, chrtmp1 )
  call i_to_s_left ( npline, chrtmp2 )
  call i_to_s_left ( kmax, chrtmp3 )

  if ( lpmoda == 1 ) then
    format1 = '(a' // chrtmp1 // ',' // chrtmp2 // 'i' // chrtmp3 // ')'
  else
    format1 = '(' // chrtmp2 // 'i' // chrtmp3 // ')'
  end if
 
  call s_blank_delete ( format1 )
 
  if ( lpmoda == 1 ) then
    format2 = '(' // chrtmp1 // 'x,' // chrtmp2 // 'i' // chrtmp3 // ')'
  else
    format2 = '(' // chrtmp2 // 'i' // chrtmp3 // ')'
  end if
 
  call s_blank_delete ( format2 )
 
  do jmin = jlo, jhi, npline
 
    jmax = min ( jmin+npline-1, jhi )
    lab = '    '
!
!  Handle a column vector.
!
    if ( jlo == jhi .and. ilo /= ihi ) then
 
      output = ' '
      call s_write ( iounit, output )
 
      if ( ilo == 1 ) then
        output = title
        call s_write ( iounit, output )
        output = ' '
        call s_write ( iounit, output )
        call i_to_s_left ( jlo, chrtmp1 )
        output = 'Column ' // chrtmp1 // ' (transposed).'
        call s_blanks_delete ( output )
        call s_write ( iounit, output )
      end if
 
      do imin = ilo, ihi, npline
 
        imax = min ( imin+npline-1, ihi )
 
        output = ' '
        call s_write ( iounit, output )
 
        none = 0
 
        do i = imin, imax
          if ( iabot(i,jlo) == 1 ) then
            ione = 3+(i-imin+1)*kmax
            output(ione:ione) = ' '
          else
            none = 1
          end if
        end do
 
        if ( lpmoda == 1 ) then
          write ( output, format1 ) lab, ( iatop(i,jlo), i = imin, imax )
          call s_write ( iounit, output )
          if ( none == 1 ) then
            write ( output, format1 ) lab, ( iabot(i,jlo), i = imin, imax )
            call s_write ( iounit, output )
          end if
        else
          write ( output, format1 ) ( iatop(i,jlo), i = imin, imax )
          call s_write ( iounit, output )
          write ( output, format1 ) ( iabot(i,jlo), i = imin, imax )
          if ( none == 1 ) then
            write ( output, format1 ) lab, ( iabot(i,jlo), i = imin, imax )
            call s_write ( iounit, output )
          end if
        end if
 
      end do
 
      go to 90
 
    end if
!
!  Handle a 2D array or table.
!
    output = ' '
    call s_write ( iounit, output )
 
    if ( jmin == 1 ) then
      output = title
      call s_write ( iounit, output )
      output = ' '
      call s_write ( iounit, output )
    end if
 
    if ( lpmoda == 1 ) then
 
      write ( output, format2 ) ( j, j = jmin, jmax )
 
      if ( jmin <= ncol-1 .and. ncol-1<=jmax ) then
        izlo = llab + ( ( ncol - 1 ) - jmin ) * kmax + kmax - 2
        izhi = izlo + 2
        output(izlo:izhi)='  P'
      end if
 
      if ( jmin <= ncol .and. ncol <= jmax ) then
        iclo = llab + ( ncol - jmin ) * kmax + kmax - 2
        ichi = iclo + 2
        output(iclo:ichi) = '  C'
      end if
 
      call s_write ( iounit, output )
 
      output = ' '
      call s_write ( iounit, output )
 
    else
 
      if ( 1 < jmin .or. jmax < ncol .or. 1 < ilo .or. ihi < nrow ) then
        call i_to_s_left ( jmin, chrtmp1 )
        call i_to_s_left ( jmax, chrtmp2 )
        output = 'Columns ' // chrtmp1 // ' to ' // chrtmp2
        call s_blanks_delete ( output )
        call s_write ( iounit, output )
        output = ' '
        call s_write ( iounit, output )
      end if
 
    end if
 
    do i = ilo, ihi
 
      if ( lpmoda == 1 ) then
 
        if ( i < nrow ) then
 
          if ( ibase(i) < 10 ) then
            write ( lab, '(''X'',i1)' ) ibase(i)
          else
            write ( lab, '(''X'',i2)' ) ibase(i)
          end if
 
        else if ( i < ihi ) then
          lab = 'Obj2'
        else
          lab = 'Obj '
        end if
 
        if ( maxrow == 1 ) then
          lab = '    '
        end if
 
      end if
 
      if ( lpmoda == 1 ) then
        write ( output, format1 ) lab, ( iatop(i,j), j = jmin, jmax )
        call s_write ( iounit, output )
        lab = '    '
        write ( output, format1 ) lab, ( iabot(i,j), j = jmin, jmax )
      else
        write ( output, format1 ) ( iatop(i,j), j = jmin, jmax )
        call s_write ( iounit, output )
        write ( output, format1 ) ( iabot(i,j), j = jmin, jmax )
      end if
!
!  Delete each denominator that is 1.  If all are 1, don't
!  even print out the line.
!
      none = 0
 
      do j = jmin, jmax
 
        if ( iabot(i,j) == 1 ) then
          ione = llab + (j-jmin+1) * kmax
          output(ione:ione) = ' '
        else
          none = 1
        end if
 
      end do
 
      if ( none == 1 ) then
        call s_write ( iounit, output )
      end if
 
      if ( jmax == jhi .and. i==ihi ) then
      else
        output = ' '
        call s_write ( iounit, output )
      end if
 
    end do
 
90      continue
 
  end do
 
  return
end
subroutine rat_read ( itop, ibot, line, prompt, iounit, ierror )

!*****************************************************************************80
!
!! rat_read() reads a rational value, expressed as integer, decimal or fraction.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ITOP, IBOT, the top and bottom of the
!    fraction that was read.
!
!    Workspace, character ( len = 255 ) LINE, used to hold the user's input.
!
!    Workspace, character ( len = 255 ) PROMPT.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Output, integer IERROR.
!    0, no error occurred.
!    1, an error occurred.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i_gcd
  integer ibot
  integer ibot1
  integer ibot2
  integer ierror
  integer iounit(4)
  integer itemp
  integer itop
  integer itop1
  integer itop2
  integer lchar
  character ( len = 255 ) line
  character ( len = 255 ) prompt

  ierror = 0
  itop = 0
  ibot = 1
 
10    continue
 
  call chrinp ( ierror, iounit, line, prompt )
  if ( ierror /= 0 ) then
    return
  end if
 
  if ( line == ' ' ) then
    go to 10
  end if

  call chrctf ( line, itop1, ibot1, ierror, lchar )
 
  if ( len ( line ) <= lchar ) then
    itop = itop1
    ibot = ibot1
  else if ( line(lchar+1:lchar+1) /= '/' ) then
    itop = itop1
    ibot = ibot1
  else
    lchar = lchar + 1
    call s_chop ( line, 1, lchar )
    call chrctf ( line, itop2, ibot2, ierror, lchar )
    itop = itop1 * ibot2
    ibot = ibot1 * itop2
  end if
 
  call s_chop ( line, 1, lchar )
!
!  Make sure fraction is in lowest terms.
!
  itemp = i_gcd ( itop, ibot )
  itop = itop / itemp
  ibot = ibot / itemp
 
  return
end
subroutine rat_to_dec ( iatop, iabot, ierror )

!*****************************************************************************80
!
!! rat_to_dec() converts a rational value to a decimal value.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer IATOP, IABOT.
!    On input, the rational value (IATOP/IABOT) to be converted.
!    On output, the rational decimal value IATOP * 10^IABOT.
!
!    Output, integer IERROR, 0 for no error, 1 for an error.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r
  integer iabot
  integer iatop
  integer ierror

  if ( iabot /= 0 ) then

    ierror = 0

    r = real ( iatop, kind = rk ) / real ( iabot, kind = rk )
 
    call r_to_dec ( r, iatop, iabot )

  else

    ierror = 1
    iatop = 0
    iabot = 0

  end if
 
  return
end
subroutine rat_to_r ( a, iatop, iabot )

!*****************************************************************************80
!
!! rat_to_r() converts rational values to real values.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) A, the value of the rational quantity.
!
!    Input, integer IATOP, IABOT, the rational quantity
!    (IATOP/IABOT) that is to be converted.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  integer iabot
  integer iatop

  if ( iabot == 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'ratrel(): Warning!'
    write ( *, * ) '  The input fraction had a zero denominator.'
    a = 0.0D+00
  else
    a = real ( iatop, kind = rk ) / real ( iabot, kind = rk )
  end if
 
  return
end
subroutine rat_to_s_left ( ival, jval, string )

!*****************************************************************************80
!
!! rat_to_s_left() returns a left-justified representation of IVAL/JVAL.
!
!  Discussion:
!
!    If the ratio is negative, a minus sign precedes IVAL.
!    A slash separates IVAL and JVAL.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer IVAL, JVAL, the two integers whose
!    ratio IVAL/JVAL is to be represented.
!
!    If IVAL is nonzero and JVAL is 0, STRING will
!    be returned as "Inf" or "-Inf" (Infinity), and if both
!    IVAL and JVAL are zero, STRING will be returned as "NaN"
!    (Not-a-Number).
!
!    Output, character ( len = 22 ) STRING, a left-justified string
!    containing the representation of IVAL/JVAL.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ival
  integer ival2
  integer jval
  integer jval2
  character ( len = * ) string
!
!  Take care of simple cases right away.
!
  if ( ival == 0 ) then
 
    if ( jval /= 0 ) then
      string = '0'
    else
      string = 'NaN'
    end if
 
  else if ( jval == 0 ) then
 
    if ( 0 < ival ) then
      string = 'Inf'
    else
      string = '-Inf'
    end if
!
!  Make copies of IVAL and JVAL.
!
  else
 
    ival2 = ival
    jval2 = jval
 
    if ( jval2 == 1 ) then
      write ( string, '(i11)' ) ival2
    else
      write ( string, '(i11, ''/'', i10)' ) ival2, jval2
    end if
 
    call s_blank_delete ( string )
 
  end if
 
  return
end
subroutine rmat_house_axh ( lda, n, a, v, ah )

!*****************************************************************************80
!
!! rmat_house_axh() computes A*H where H is a compact Householder matrix.
!
!  Discussion:
!
!    The Householder matrix H(V) is defined by
!
!      H(V) = I - 2 * v * v' / ( v' * v )
!
!    This routine is not particularly efficient.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of A.
!
!    Input, integer N, the order of A.
!
!    Output, real ( kind = rk ) A(LDA,N), the matrix.
!
!    Input, real ( kind = rk ) V(N), a vector defining a Householder matrix.
!
!    Output, real ( kind = rk ) AH(LDA,N), the product A*H.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer lda
  integer n

  real ( kind = rk ) a(lda,n)
  real ( kind = rk ) ah(lda,n)
  real ( kind = rk ) ah_temp(n,n)
  integer i
  integer j
  integer k
  real ( kind = rk ) v(n)
  real ( kind = rk ) v_normsq

  v_normsq = 0.0D+00
  do i = 1, n
    v_normsq = v_normsq + v(i) ** 2
  end do
!
!  Compute A*H' = A*H
!
  do i = 1, n
    do j = 1, n
      ah_temp(i,j) = a(i,j)
      do k = 1, n
        ah_temp(i,j) = ah_temp(i,j) - 2.0D+00 * a(i,k) * v(k) * v(j) / v_normsq
      end do
    end do
  end do
!
!  Copy the temporary result into AH.
!  Doing it this way means the user can identify the input arguments A and AH.
!
  ah(1:n,1:n) = ah_temp(1:n,1:n)

  return
end
subroutine row_add ( a, iatop, iabot, ierror, iform, iounit, irow1, irow2, &
  maxcol, maxrow, ncol, sval, istop, isbot )

!*****************************************************************************80
!
!! row_add() adds a multiple of one row to another.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = rk ) A(MAXROW,MAXCOL), the current matrix.
!
!    Input/output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
!    IATOP and IABOT represent the current rational or decimal matrix.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred.
!
!    Input, integer IFORM, 0/1/2, rational/real/decimal arithmetic.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Input, integer IROW1, the row which is to be modified.
!
!    Input, integer IROW2, the row which is to be multiplied by
!    a given value and added to row IROW1.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
!    Input, integer NCOL, the number of columns in the matrix.
!
!    Input, real ( kind = rk ) SVAL, the multiplier to use if real 
!    arithmetic is employed.
!
!    Input, integer ISTOP, ISBOT, the fractional or decimal
!    multiplier to use.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxcol
  integer maxrow

  real ( kind = rk ) a(maxrow,maxcol)
  character ( len = 10 ) chrtmp1
  character ( len = 10 ) chrtmp2
  character ( len = 22 ) chrtmp3
  integer iabot(maxrow,maxcol)
  integer iatop(maxrow,maxcol)
  integer ibot
  integer ierror
  integer iform
  integer iounit(4)
  integer irow1
  integer irow2
  integer isbot
  integer isbot2
  integer istop
  integer istop2
  integer itop
  integer j
  integer len1
  integer len2
  integer len3
  integer ncol
  character ( len = 255 ) output
  real ( kind = rk ) sval
!
!  Return immediately if the multiplier is zero.
!
  if ( ( iform == 0 .and. istop == 0 ) .or. &
       ( iform == 1 .and. sval == 0.0D+00 ) .or. &
       ( iform == 2 .and. istop == 0 ) ) then
      return
    end if
!
!  Carry out the operation.
!
  if ( iform == 0 ) then
 
    do j = 1, ncol
 
      call rat_mul ( isbot2, isbot, iabot(irow2,j), istop2, istop, &
        iatop(irow2,j), ierror )
 
      if ( ierror /= 0 ) then
        return
      end if

      call rat_add ( ibot, iabot(irow1,j), isbot2, itop, iatop(irow1,j), istop2 )
 
      iatop(irow1,j) = itop
      iabot(irow1,j) = ibot
 
    end do
 
  else if ( iform == 1 ) then
 
    do j = 1, ncol
      a(irow1,j) = a(irow1,j) + sval * a(irow2,j)
    end do
 
  else if ( iform == 2 ) then
 
    do j = 1, ncol
 
      call dec_mul ( isbot2, isbot, iabot(irow2,j), istop2, istop, &
        iatop(irow2,j) )
 
      call dec_add ( ibot, iabot(irow1,j), isbot2, itop, &
        iatop(irow1,j), istop2 )
 
      iatop(irow1,j) = itop
      iabot(irow1,j) = ibot
 
    end do
 
  end if
!
!  Print out a message.
!
  if ( iform == 0 ) then
 
    if ( istop == isbot ) then
      chrtmp3 = '+'
    else if ( istop == - isbot ) then
      chrtmp3 = '-'
    else
      call rat_to_s_left ( istop, isbot, chrtmp3 )
    end if
 
  else if ( iform == 1 ) then
 
    if ( sval == 1.0D+00 ) then
      chrtmp3 = '+'
    else if ( sval == - 1.0D+00 ) then
      chrtmp3 = '-'
    else
      call r_to_s_left ( sval, chrtmp3 )
    end if
 
  else if ( iform == 2 ) then

    if ( istop == 1 .and. isbot == 0 ) then 
      chrtmp3 = '+'
    else if ( istop == - 1 .and. isbot == 0 ) then
      chrtmp3 = '-'
    else
      call dec_to_s_left ( istop, isbot, chrtmp3 )
    end if
 
  end if
 
  call i_to_s_left ( irow1, chrtmp1 )
  len1 = len_trim ( chrtmp1 )

  call i_to_s_left ( irow2, chrtmp2 )
  len2 = len_trim ( chrtmp2 )

  len3 = len_trim ( chrtmp3 )

  if ( chrtmp3 == '-' .or. chrtmp3 == '+' ) then

    output = '  ERO: Row ' // chrtmp1(1:len1) // ' <=  ' // 'Row ' // &
      chrtmp1(1:len1) // ' ' // chrtmp3(1:1) // ' Row ' // chrtmp2(1:len2)

  else if ( chrtmp3(1:1) == '-' .or. chrtmp3(1:1) == '+' ) then

    output = '  ERO: Row ' // chrtmp1(1:len1) // ' <=  ' // 'Row ' // &
      chrtmp1(1:len1) // ' ' // chrtmp3(1:1) // ' ' // chrtmp3(2:len3) // &
      ' Row ' // chrtmp2(1:len2)

  else

    output = '  ERO: Row ' // chrtmp1(1:len1) // ' <=  ' // 'Row ' // &
      chrtmp1(1:len1) // ' + ' // chrtmp3(1:len3) // ' Row ' // chrtmp2(1:len2)
  end if

  call s_write ( iounit, output )
 
  return
end
subroutine row_add_param ( ierror, iform, iounit, irow1, irow2, istop, isbot, &
  line, nrow, sval )

!*****************************************************************************80
!
!! row_add_param() gets and checks the row add parameters.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred.
!
!    Input, integer IFORM, 0/1/2, rational/real/decimal arithmetic.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Input, integer IROW1, the row to which the multiple is to be added.
!
!    Input, integer IROW2, the row which is to be multiplied and
!    added to another row.
!
!    Output, integer ISTOP, ISBOT, the parts of the rational
!    or decimal fraction of the multiplier, if that is the
!    arithmetic being used.
!
!    Workspace, character ( len = 255 ) LINE, used to hold the user's input.
!
!    Input, integer NROW, the number of rows in the matrix.
!
!    Output, real ( kind = rk ) SVAL, the multiplier, if real arithmetic is used.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ierror
  integer iform
  integer iounit(4)
  integer irow1
  integer irow2
  integer isbot
  integer istop
  character ( len = 255 ) line
  integer nrow
  character ( len = 255 ) output
  character ( len = 255 ) prompt
  real ( kind = rk ) sval

  prompt = 'multiplier S, row I to add, target row J.'
!
!  Get the multiplier, SVAL or ISTOP/ISBOT.
!
  if ( iform == 0 ) then
 
    call rat_read ( istop, isbot, line, prompt, iounit, ierror )
    if ( ierror /= 0 ) then
      return
    end if
 
  else if ( iform == 1 ) then
 
    call r_read ( sval, line, prompt, iounit, ierror )
    if ( ierror /= 0 ) then
      return
    end if
 
  else if ( iform == 2 ) then
 
    call dec_read ( istop, isbot, line, prompt, iounit, ierror )

    if ( ierror /= 0 ) then
      output = ' '
      call s_write ( iounit, output )
      output = 'CHECK_ADD(): Error!'
      call s_write ( iounit, output )
      output = '  DEC_READ returns error code.'
      call s_write ( iounit, output )
      return
    end if
 
    call dec_round ( istop, isbot )
 
  end if
!
!  Get the row to add, IROW2.
!
  call i_read ( irow2, line, prompt, iounit, ierror )

  if ( ierror /= 0 ) then
    ierror = 1
    output = 'Error reading row index from line:'
    call s_write ( iounit, output )
    output = line
    call s_write ( iounit, output )
    return
  end if

  if ( irow2 < 1 .or. nrow < irow2 ) then
    ierror = 1
    output = 'Error!  Row index was not acceptable!'
    call s_write ( iounit, output )
    return
  end if
!
!  Get the row to which we are adding, IROW1.
!
  call i_read ( irow1, line, prompt, iounit, ierror )
  if ( ierror /= 0 ) then
    return
  end if
 
  if ( irow1 < 1 .or. nrow < irow1 ) then
    ierror = 1
    output = 'Error!  The row index was not acceptable!'
    call s_write ( iounit, output )
    return
  end if
!
!  Make sure the rows are different.
!
  if ( irow1 == irow2 ) then
    output = 'Error!  The rows should not be the same!'
    call s_write ( iounit, output )
    ierror = 1
    return
  end if
 
  ierror = 0
 
  return
end
subroutine row_auto ( a, iatop, iabot, ibase, ierror, iform, iounit, maxcol, &
  maxrow, ncol, nrow )

!*****************************************************************************80
!
!! row_auto() automatically row reduces the current matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = rk ) A(MAXROW,MAXCOL), the current matrix.
!
!    Input/output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
!    IATOP and IABOT represent the rational or decimal matrix
!    to which elementary row operations will be applied.
!
!    Input/output, integer IBASE(MAXROW).  IBASE is information
!    really only used by the linear programming routines.
!    AUTO_ERO only needs it because some lower level routines are
!    shared with the linear programming routines.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred.
!
!    Input, integer IFORM, 0/1/2, rational/real/decimal arithmetic.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
!    Input, integer NCOL, the number of columns in the matrix.
!
!    Input, integer NROW, the number of rows in the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxcol
  integer maxrow

  real ( kind = rk ) a(maxrow,maxcol)
  real ( kind = rk ) amax
  real ( kind = rk ) atemp
  integer i
  integer iabot(maxrow,maxcol)
  integer iatop(maxrow,maxcol)
  integer ibase(maxrow) 
  integer ierror
  integer iform
  integer imax
  integer iounit(4)
  integer irow
  integer isbot
  integer istop
  integer j
  integer jcol
  integer krow
  integer l
  integer lpmoda
  integer lrow
  integer ncol
  integer nrow
  character ( len = 255 ) output
  real ( kind = rk ) sval

  ierror = 0
 
  do i = 1, nrow
 
    irow = i
 
    do j = 1, ncol
 
      jcol = j
!
!  In column JCOL, seek the row between IROW and NROW with
!  maximum nonzero entry AMAX.
!
      imax = 0
      amax = 0.0D+00
 
      do krow = irow, nrow
 
        if ( iform == 0 ) then
          call rat_to_r ( atemp, iatop(krow,jcol), iabot(krow,jcol) )
        else if ( iform == 1 ) then
          atemp = a(krow,jcol)
        else if ( iform == 2 ) then
          call dec_to_r ( atemp, iatop(krow,jcol), iabot(krow,jcol) )
        end if
 
        atemp = abs ( atemp )
 
        if ( amax < atemp ) then
          amax = atemp
          imax = krow
        end if
 
      end do
 
      if ( imax /= 0 ) then
        krow = imax
        go to 10
      end if
 
    end do
 
    return
 
10      continue
 
    output = ' '
    call s_write ( iounit, output )
!
!  Interchange the IROW-th and the pivot rows.
!
    if ( krow /= irow ) then
      lpmoda = 0
      call row_swap ( a, iatop, iabot, ibase, ierror, iform, iounit, &
        krow, irow, lpmoda, maxcol, maxrow, ncol, nrow )
    end if
!
!  Divide the pivot row by A(IROW,JCOL) so that A(IROW,JCOL) = 1.
!
    if ( iform == 0 ) then
 
      istop = iatop(irow,jcol)
      isbot = iabot(irow,jcol)
 
    else if ( iform == 1 ) then
 
      sval = a(irow,jcol)
 
    else if ( iform == 2 ) then
 
      istop = iatop(irow,jcol)
      isbot = iabot(irow,jcol)
 
    end if
 
    call row_div ( a, iatop, iabot, ierror, iform, iounit, irow, &
      maxcol, maxrow, ncol, nrow, sval, istop, isbot )
!
!  Annihilate A(L,JCOL) for L not equal to IROW.
!
    do l = 1, nrow
 
      lrow = l
 
      if ( lrow /= irow ) then
 
        if ( iform == 0 ) then

          if ( iatop(lrow,jcol) /= 0 ) then

            istop = - iatop(lrow,jcol)
            isbot = iabot(lrow,jcol)

            call row_add ( a, iatop, iabot, ierror, iform, iounit, &
              lrow, irow, maxcol, maxrow, ncol, sval, istop, isbot )

            iatop(lrow,jcol) = 0
            iabot(lrow,jcol) = 1

          end if

        else if ( iform == 1 ) then

          if ( a(lrow,jcol) /= 0.0D+00 ) then

            sval = - a(lrow,jcol)

            call row_add ( a, iatop, iabot, ierror, iform, iounit, &
              lrow, irow, maxcol, maxrow, ncol, sval, istop, isbot )

            a(lrow,jcol) = 0.0D+00

          end if

        else if ( iform == 2 ) then

          if ( iatop(lrow,jcol) /= 0 ) then

            istop = - iatop(lrow,jcol)
            isbot = iabot(lrow,jcol)

            call row_add ( a, iatop, iabot, ierror, iform, iounit, &
              lrow, irow, maxcol, maxrow, ncol, sval, istop, isbot )

            iatop(lrow,jcol) = 0
            iabot(lrow,jcol) = 0

          end if

        end if
 
      end if
 
    end do
 
  end do
 
  return
end
subroutine row_del ( a, iabot, iatop, irow, maxcol, maxrow, ncol, nrow )

!*****************************************************************************80
!
!! row_del() deletes a row by shifting other rows up.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = rk ) A(MAXROW,MAXCOL), the current matrix.
!
!    Output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
!    IATOP and IABOT represent the rational or decimal matrix
!    to be changed.
!
!    Input, integer IROW, the row to be deleted.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
!    Input, integer NCOL, the number of columns in the matrix.
!
!    Input, integer NROW, the number of rows in the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxcol
  integer maxrow

  real ( kind = rk ) a(maxrow,maxcol)
  integer i
  integer iabot(maxrow,maxcol)
  integer iatop(maxrow,maxcol)
  integer irow
  integer j
  integer ncol
  integer nrow

  do i = irow, nrow - 1
    do j = 1, ncol
      a(i,j) = a(i+1,j)
      iatop(i,j) = iatop(i+1,j)
      iabot(i,j) = iabot(i+1,j)
    end do
  end do
 
  return
end
subroutine row_div ( a, iatop, iabot, ierror, iform, iounit, irow, &
  maxcol, maxrow, ncol, nrow, sval, istop, isbot )

!*****************************************************************************80
!
!! row_div() divides row IROW of the A matrix by a scale factor.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = rk ) A(MAXROW,MAXCOL), the current matrix.
!
!    Input/output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
!    IATOP and IABOT represent the current rational or decimal
!    matrix.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred.
!
!    Input, integer IFORM, 0/1/2, rational/real/decimal arithmetic.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Input, integer IROW, the row to be divided.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
!    Input, integer NCOL, the number of columns in the matrix.
!
!    Input, integer NROW, the number of rows in the matrix.
!
!    Input, real ( kind = rk ) SVAL, the real divisor.
!
!    Input, integer ISTOP, ISBOT, the fractional or decimal divisor.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxcol
  integer maxrow

  real ( kind = rk ) a(maxrow,maxcol)
  character ( len = 10 ) chrtmp1
  character ( len = 24 ) chrtmp2
  integer iabot(maxrow,maxcol)
  integer iatop(maxrow,maxcol)
  integer ibot
  integer ierror
  integer iform
  integer iounit(4)
  integer irow
  integer isbot
  integer istop
  integer itop
  integer j
  integer len1
  integer len2
  integer ncol
  integer nrow
  character ( len = 3 ) op
  character ( len = 255 ) output
  real ( kind = rk ) sval
!
!  Make sure that the row number is legal.
!
  if ( irow < 1 .or. nrow < irow ) then
    output = 'Error!  The row number is out of range!'
    call s_write ( iounit, output )
    ierror = 1
    return
  end if
!
!  Check for an illegal divisor of 0, or a pointless divisor of 1.
!
  if ( iform == 0 ) then
    if ( istop == 0 ) then
      output = 'Error!  It is illegal to divide by 0!'
      call s_write ( iounit, output )
      ierror = 1
      return
    else if ( istop == isbot ) then
      return
    end if
  else if ( iform == 1 ) then
    if ( sval == 0.0D+00 ) then
      output = 'Error!  It is illegal to divide by 0!'
      call s_write ( iounit, output )
      ierror = 1
      return
    else if ( sval == 1.0D+00 ) then
      return
    end if
  else if ( iform == 2 ) then
    if ( istop == 0 ) then
      output = 'Error!  It is illegal to divide by 0!'
      call s_write ( iounit, output )
      ierror = 1
      return
    else if ( istop == 1 .and. isbot==0 ) then
      return
    end if
  end if
!
!  Carry out the division.
!
  if ( iform == 0 ) then
 
    do j = 1, ncol
 
      call rat_mul ( ibot, iabot(irow,j), istop, itop, iatop(irow,j), isbot, &
        ierror )
 
      if ( ierror /= 0 ) then
        return
      end if
 
      iatop(irow,j) = itop
      iabot(irow,j) = ibot
 
    end do
 
  else if ( iform == 1 ) then
 
    do j = 1, ncol
      a(irow,j) = a(irow,j) / sval
    end do
 
  else if ( iform == 2 ) then
 
    do j = 1, ncol
 
      call dec_div ( iatop(irow,j), iabot(irow,j), istop, isbot, &
        iatop(irow,j), iabot(irow,j), ierror )
 
    end do
 
  end if
!
!  Print out a statement about what has been done.
!

  if ( iform == 0 ) then
 
    if ( isbot == 1 ) then

      call i_to_s_left ( istop, chrtmp2 )
      op = ' / '

    else

      call rat_to_s_left ( isbot, istop, chrtmp2 )
      op = ' * '

    end if
 
  else if ( iform == 1 ) then
 
    call r_to_s_left ( sval, chrtmp2 )
    op = ' / '

  else if ( iform == 2 ) then
 
    call dec_to_s_left ( istop, isbot, chrtmp2 )
    op = ' / '

  end if

  call i_to_s_left ( irow, chrtmp1 )
  len1 = len_trim ( chrtmp1 )

  len2 = len_trim ( chrtmp2 )

  output = '  ERO: Row ' // chrtmp1(1:len1) // ' <=  Row ' // chrtmp1(1:len1) &
    // op // chrtmp2(1:len2)

  call s_write ( iounit, output )
 
  return
end
subroutine row_div_param ( ierror, iform, iounit, irow, isbot, istop, line, &
  sval )

!*****************************************************************************80
!
!! row_div_param() gets and checks the row divide parameters.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred.
!
!    Input, integer IFORM, 0/1/2, rational/real/decimal arithmetic.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Output, integer IROW, the row to be divided.
!
!    Workspace, character ( len = 255 ) LINE, used to hold the user's input.
!
!    Output, real ( kind = rk ) SVAL, the real divisor.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ierror
  integer iform
  integer iounit(4)
  integer irow
  integer isbot
  integer istop
  character ( len = 255 ) line
  character ( len = 255 ) output
  character ( len = 255 ) prompt
  real ( kind = rk ) sval

  prompt = 'row I, divisor S.'
!
!  Read the row number to be divided.
!
  call i_read ( irow, line, prompt, iounit, ierror )
  if ( ierror /= 0 ) then
    return
  end if
!
!  Read the divisor, either RVAL or ISTOP/ISBOT.
!
  if ( iform == 0 ) then
 
    call rat_read ( istop, isbot, line, prompt, iounit, ierror )
 
  else if ( iform == 1 ) then
 
    call r_read ( sval, line, prompt, iounit, ierror )
 
  else if ( iform == 2 ) then
 
    call dec_read ( istop, isbot, line, prompt, iounit, ierror )

    if ( ierror /= 0 ) then
      output = ' '
      call s_write ( iounit, output )
      output = 'DIVIDE - Fatal error!'
      call s_write ( iounit, output )
      output = '  DEC_READ returned error flag.'
      call s_write ( iounit, output )
      return
    end if

    call dec_round ( istop, isbot )
 
  end if
 
  if ( ierror /= 0 ) then
    return
  end if

  return
end
subroutine row_mul ( a, iatop, iabot, ierror, iform, iounit, irow, maxcol, &
  maxrow, ncol, nrow, sval, istop, isbot )

!*****************************************************************************80
!
!! row_mul() multiplies a row of the matrix by a scale factor.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = rk ) A(MAXROW,MAXCOL), the current matrix.
!
!    Input/output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
!    IATOP and IABOT represent the current rational or decimal matrix.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred.
!
!    Input, integer IFORM, 0/1/2, rational/real/decimal arithmetic.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Input, integer IROW, the row that is to be multiplied.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
!    Input, integer NCOL, the number of columns in the matrix.
!
!    Input, integer NROW, the number of rows in the matrix.
!
!    Input, real ( kind = rk ) SVAL, the real row multiplier.
!
!    Input, integer ISTOP, ISBOT, the decimal or fractional row multiplier.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxcol
  integer maxrow

  real ( kind = rk ) a(maxrow,maxcol)
  character ( len = 10 ) chrtmp1
  character ( len = 22 ) chrtmp
  integer iabot(maxrow,maxcol)
  integer iatop(maxrow,maxcol)
  integer ibot
  integer ierror
  integer iform
  integer iounit(4)
  integer irow
  integer isbot
  integer istop
  integer itop
  integer j
  integer len1
  integer len2
  integer ncol
  integer nrow
  character ( len = 255 ) output
  real ( kind = rk ) sval
!
!  Make sure row number is OK.
!
  if ( irow < 1 .or. nrow < irow ) then
    output = 'Error!  The row number is out of range!'
    call s_write ( iounit, output )
    ierror = 1
    return
  end if
!
!  For rational arithmetic, make sure bottom of scale factor
!  is not 0.
!
  if ( iform == 0 ) then
    if ( isbot == 0 ) then
      output = 'Error!  Illegal 0 divisor in multiplier!'
      call s_write ( iounit, output )
      ierror = 1
      return
    end if
  end if
!
!  Check for multiplication by 0.
!
  if (  ( iform == 0 .and. istop == 0 ) .or. &
       ( iform == 1 .and. sval == 0.0D+00 ) .or. &
       ( iform == 2 .and. istop == 0 ) ) then
    output = ' '
    call s_write ( iounit, output )
    output = 'Warning - Multiplication by zero is not an ERO.'
    call s_write ( iounit, output )
    ierror = 1
    return
  end if
!
!  Check for multiplication by 1.
!
  if ( iform == 0 ) then
    if ( istop == isbot ) then
      return
    end if
  else if ( iform == 1 ) then
    if ( sval == 1.0D+00 ) then
      return
    end if
  else if ( iform == 2 ) then
    if ( istop == 1 .and. isbot == 0 ) then
      return
    end if
  end if
!
!  Carry out the multiplication.
!
  if ( iform == 0 ) then
 
    do j = 1, ncol
 
      call rat_mul ( ibot, iabot(irow,j), isbot, itop, iatop(irow,j), istop, &
        ierror )

      if ( ierror /= 0 ) then
        return
      end if

      iatop(irow,j) = itop
      iabot(irow,j) = ibot
 
    end do
 
  else if ( iform == 1 ) then
 
    do j = 1, ncol
      a(irow,j) = sval * a(irow,j)
    end do
 
  else if ( iform == 2 ) then
 
    do j = 1, ncol
 
      call dec_mul ( ibot, iabot(irow,j), isbot, itop, iatop(irow,j), istop )

      if ( ierror /= 0 ) then
        return
      end if
 
      iatop(irow,j) = itop
      iabot(irow,j) = ibot
 
    end do
 
  end if
!
!  Confirm the operation.
!
  if ( iform == 0 ) then
 
    if ( istop == - isbot ) then
      chrtmp = '-'
    else
      call rat_to_s_left ( istop, isbot, chrtmp )
    end if
 
  else if ( iform == 1 ) then
 
    if ( sval == - 1.0D+00 ) then
      chrtmp = '-'
    else
      call r_to_s_left ( sval, chrtmp )
    end if
 
  else if ( iform == 2 ) then

    if ( istop == - 1 .and. isbot == 0 ) then
      chrtmp = '-'
    else
      call dec_to_s_left ( istop, isbot, chrtmp )
    end if

  end if
 
  call i_to_s_left ( irow, chrtmp1 )
  len1 = len_trim ( chrtmp1 )

  len2 = len_trim ( chrtmp )
  output = '  ERO: Row ' // chrtmp1(1:len1) // ' <= ' // chrtmp(1:len2) // &
    ' Row ' // chrtmp1(1:len1)
  call s_write ( iounit, output )

  return
end
subroutine row_mul_param ( ierror, iform, iounit, irow, istop, isbot, line, &
  rval )

!*****************************************************************************80
!
!! row_mul_param() gets and checks the row multiply parameters.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred.
!
!    Input, integer IFORM, 0/1/2, rational/real/decimal arithmetic.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Output, integer IROW, the row to be multiplied.
!
!    Output, integer ISTOP, ISBOT, the multiplier to use for
!    fractional or decimal arithmetic.
!
!    Workspace, character ( len = 255 ) LINE, used to hold the user's input.
!
!    Output, real ( kind = rk ) RVAL, the multiplier to use for real arithmetic.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ierror
  integer iform
  integer iounit(4)
  integer irow
  integer isbot
  integer istop
  character ( len = 255 ) line
  character ( len = 255 ) prompt
  real ( kind = rk ) rval

  prompt = 'row I, multiplier S.'
!
!  Read the row number to be multiplied.
!
  call i_read ( irow, line, prompt, iounit, ierror )
  if ( ierror /= 0 ) then
    return
  end if
!
!  Read the multiplier, either RVAL or ISTOP/ISBOT.
!
  if ( iform == 0 ) then
 
    call rat_read ( istop, isbot, line, prompt, iounit, ierror )
 
  else if ( iform == 1 ) then
 
    call r_read ( rval, line, prompt, iounit, ierror )
 
  else if ( iform == 2 ) then
 
    call dec_read ( istop, isbot, line, prompt, iounit, ierror )
 
    call dec_round ( istop, isbot )
 
  end if
 
  return
end
subroutine row_op_check ( command, ierror, iounit, line2 )

!*****************************************************************************80
!
!! row_op_check() checks for commands given in the form of ERO's.
!
!  Discussion:
!
!    The form of the elementary row operation commands includes:
!
!    The row interchange command:
!      RI1 <=> RI2
!    This will fail if user types "R I1 <=> R I2"
!
!    The scalar multiply command:
!      RI1 <= S * RI1
!    with or without the "*".
!
!    The scalar divide command:
!      RI1 <= RI1 / S
!
!    The add row command:
!      RI1 <= RI1 + S * RI2
!    or
!      RI1 <= S * RI2 + RI1
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = 4 ) COMMAND.
!    If the routine decides that the user has input an ERO in the
!    natural format, then COMMAND contains the necessary
!    one letter MATMAN command to carry out the ERO.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Workspace, character ( len = 255 ) LINE2, a copy of the user input in LINE.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = 22 ) chrtmp
  character ( len = 10 ) chrtmp1
  character ( len = 10 ) chrtmp2
  character ( len = 20 ) command
  integer idbot2
  integer idbot3
  integer idtop2
  integer idtop3
  integer ierror
  integer iounit(4)
  integer irow1
  integer irow2
  integer irow3
  integer isbot2
  integer isbot3
  integer istop2
  integer istop3
  integer lchar
  logical ldiv
  character ( len = 255 ) line2
  character ( len = 255 ) output
  character ( len = 255 ) string

  command = ' '
!
!  1. Remove all blanks from the line, and capitalize it.
!
  call s_blank_delete ( line2 )
  call s_cap ( line2 )
!
!  2. Is the first character an "R" or "ROW"?
!
  if ( line2(1:1) /= 'R' ) then
    return
  end if
 
  if ( line2(1:3) == 'ROW' ) then
    call s_chop ( line2, 1, 3 )
  else
    call s_chop ( line2, 1, 1 )
  end if
!
!  3. The next item should be a row number, IROW1.
!
  call s_to_i ( line2, irow1, ierror, lchar )
 
  if ( ierror /= 0 ) then
    output = 'Your ERO command could not be understood.'
    call s_write ( iounit, output )
    output = 'The first row number "R1" did not make sense.'
    call s_write ( iounit, output )
    return
  end if
 
  call s_chop ( line2, 1, lchar )
!
!  4. Check for the row interchange string "=", "<>", "<=>" or "<->".
!
  if ( line2(1:2) == '<>' ) then
    string = '<>'
  else if ( line2(1:3) == '<=>' ) then
    string = '<=>'
  else if ( line2(1:3) == '<->' ) then
    string = '<->'
  else if ( line2(1:2) == '<=' ) then
    string = '<='
  else if ( line2(1:2) == '<-' ) then
    string = '<-'
  else if ( line2(1:2) == '=>' ) then
    string = '=>'
  else if ( line2(1:2) == '->' ) then
    string = '->'
  else if ( line2(1:1) == '=' ) then
    string = '='
  else if ( line2(1:2) == ':=' ) then
    string = ':='
  else
    ierror = 1
    output = 'Your ERO command could not be understood.'
    call s_write ( iounit, output )
    output = 'The assignment symbol <=> was missing.'
    call s_write ( iounit, output )
    return
  end if
 
  lchar = len_trim ( string )
 
  call s_chop ( line2, 1, lchar )
!
!  5. The next quantity could be an explicit signed scalar, S2,
!     or an implicit +-1.
!
  if ( line2(1:1) == 'R' ) then

    istop2 = 1.0D+00
    isbot2 = 1.0D+00

  else

    if ( line2(1:2) == '+R' ) then
      istop2 = 1.0D+00
      isbot2 = 1.0D+00
      call s_chop ( line2, 1, 1 )
    else if ( line2(1:2) == '-R' ) then
      istop2 = - 1.0D+00
      isbot2 = 1.0D+00
      call s_chop ( line2, 1, 1 )
    else
      call chrctg ( line2, istop2, isbot2, ierror, lchar )
      call s_chop ( line2, 1, lchar )
 
      if ( ierror /= 0 ) then
        output = 'Your ERO command could not be understood.'
        call s_write ( iounit, output )
        output = 'The multiplier S2 did not make sense.'
        call s_write ( iounit, output )
        ierror = 1
        return
      end if
 
    end if
  end if
!
!  6. Is the next character an optional "*"?
!
  if ( line2(1:1) == '*' ) then
    call s_chop ( line2, 1, 1 )
  end if
!
!  7. Is the next character an "R"?
!
  if ( line2(1:3) == 'ROW' ) then
    call s_chop ( line2, 1, 3 )
  else if ( line2(1:1) == 'R' ) then
    call s_chop ( line2, 1, 1 )
  else
    ierror = 1
    output = 'Your ERO command could not be understood.'
    call s_write ( iounit, output )
    output = 'Could not find the second row index.'
    call s_write ( iounit, output )
    return
  end if
!
!  8. The next item should be a row number, IROW2.
!
  call s_to_i ( line2, irow2, ierror, lchar )
 
  if ( ierror /= 0 ) then
    output = 'Your ERO command could not be understood.'
    call s_write ( iounit, output )
    output = 'The second row number "R2" did not make sense.'
    call s_write ( iounit, output )
    ierror = 1
    return
  else
    call s_chop ( line2, 1, lchar )
  end if
!
!  9. If there's nothing more, this must be an interchange
!     or a scaling.  Form the equivalent MATMAN command.
!
  if ( line2 == ' ' ) then
 
    if ( irow1 == irow2 ) then

      command = 'ROW_MUL'

      if ( isbot2 < 0 ) then
        isbot2 = - isbot2
        istop2 = - istop2
      end if

      call rat_to_s_left ( istop2, isbot2, chrtmp )
      call i_to_s_left ( irow1, chrtmp1 )
      line2 = chrtmp1 // ' ' // chrtmp
      call s_blanks_delete ( line2 )
      return
    end if
 
    if ( istop2 == 1 .and. isbot2 == 1 ) then
      command = 'ROW_SWAP'
      call i_to_s_left ( irow1, chrtmp1 )
      call i_to_s_left ( irow2, chrtmp2 )
      line2 = chrtmp1 // ' ' // chrtmp2
      call s_blanks_delete ( line2 )
      return
    end if
 
    ierror = 1
    output = 'Your ERO command could not be understood.'
    call s_write ( iounit, output )
    output = 'A MULTIPLY command must have R1 and R2 the same.'
    call s_write ( iounit, output )
    output = 'An INTERCHANGE command cannot have a multiplier.'
    call s_write ( iounit, output )
    return
  end if
!
!  10. Is the next quantity a '/', or perhaps a '*'?
!
  ldiv = .false.
 
  if ( line2(1:1) == '/' ) then
 
    ldiv = .true.
    call s_chop ( line2, 1, 1 )
    call chrctg ( line2, idtop2, idbot2, ierror, lchar )
 
    if ( ierror /= 0 ) then
      output = 'Your ERO command could not be understood.'
      call s_write ( iounit, output )
      output = 'The divisor of row 2 did not make sense.'
      call s_write ( iounit, output )
      return
    end if
 
    istop2 = istop2 * idbot2
    isbot2 = isbot2 * idtop2
 
    if ( irow1 == irow2 ) then

      if ( ldiv ) then
        command = 'ROW_DIV'
        call i_swap ( istop2, isbot2 )
      else
        command = 'ROW_MUL'
      end if

      if ( isbot2 < 0 ) then
        isbot2 = - isbot2
        istop2 = - istop2
      end if

      call i_to_s_left ( irow1, chrtmp1 )
      call rat_to_s_left ( istop2, isbot2, chrtmp )
      line2 = chrtmp1 // ' ' // chrtmp
      call s_blanks_delete ( line2 )
      return

    end if
 
  else if ( line2(1:1) == '*' ) then
 
    call s_chop ( line2, 1, 1 )
    call chrctg ( line2, idtop2, idbot2, ierror, lchar )
 
    if ( ierror /= 0 ) then
      output = 'Your ERO command could not be understood.'
      call s_write ( iounit, output )
      output = 'The multiplier of row 2 did not make sense.'
      call s_write ( iounit, output )
      return
    end if
 
    istop2 = istop2 * idtop2
    isbot2 = isbot2 * idbot2
 
    if ( irow1 == irow2 ) then

      command = 'ROW_MUL'

      if ( isbot2 < 0 ) then
        isbot2 = - isbot2
        istop2 = - istop2
      end if

      call i_to_s_left ( irow1, chrtmp1 )
      call rat_to_s_left ( istop2, isbot2, chrtmp )
      line2 = chrtmp1 // ' ' // chrtmp
      call s_blanks_delete ( line2 )
      return

    end if
 
  end if
!
!  11. Is the next quantity a scalar, S3?
!
  if ( line2(1:2) == '+R' ) then
 
    istop3 = 1.0D+00
    isbot3 = 1.0D+00
    call s_chop ( line2, 1, 1) 
 
  else if ( line2(1:2) == '-R' ) then
 
    istop3 = - 1.0D+00
    isbot3 = 1.0D+00
    call s_chop ( line2, 1, 1 )
 
  else
 
    call chrctg ( line2, istop3, isbot3, ierror, lchar )
 
    if ( ierror /= 0 ) then
      output = 'Your ERO command could not be understood.'
      call s_write ( iounit, output )
      output = 'The multiplier S2 did not make sense.'
      call s_write ( iounit, output )
      ierror = 1
      return
    end if
 
    call s_chop ( line2, 1, lchar )
 
  end if
!
!  12. Is the next quantity an optional "*"?
!
  if ( line2(1:1) == '*' ) then
    call s_chop ( line2, 1, 1 )
  end if
!
!  13. Is the next quantity an "R" or ROW?
!
  if ( line2(1:3) == 'ROW' ) then
    call s_chop ( line2, 1, 3 )
  else if ( line2(1:1) == 'R' ) then
    call s_chop ( line2, 1, 1) 
  else
    ierror = 1
    output = 'Your ERO command could not be understood.'
    call s_write ( iounit, output )
    output = 'The "R" marking the third row was misplaced.'
    call s_write ( iounit, output )
    return
  end if
!
!  14. The next item should be a row number, IROW3.
!
  call s_to_i ( line2, irow3, ierror, lchar )
 
  if ( ierror /= 0 ) then
    output = 'Your ERO command could not be understood.'
    call s_write ( iounit, output )
    output = 'The third row number "R3" did not make sense.'
    call s_write ( iounit, output )
    ierror = 1
    return
  end if
 
  call s_chop ( line2, 1, lchar )
!
!  15. Is the next quantity a '/', or perhaps a '*'?
!
  if ( line2(1:1) == '/' ) then
 
    call s_chop ( line2, 1, 1 )
    call chrctg ( line2, idtop3, idbot3, ierror, lchar )
 
    if ( ierror /= 0 ) then
      output = 'Your ERO command could not be understood.'
      call s_write ( iounit, output )
      output = 'The divisor of row 3 did not make sense.'
      call s_write ( iounit, output )
      return
    end if
 
    istop3 = istop3 * idbot3
    isbot3 = isbot3 * idtop3
 
  else if ( line2(1:1) == '*' ) then
 
    call s_chop ( line2, 1, 1 )
    call chrctg ( line2, idtop3, idbot3, ierror, lchar )
 
    if ( ierror /= 0 ) then
      output = 'Your ERO command could not be understood.'
      call s_write ( iounit, output )
      output = 'The multiplier of row 3 did not make sense.'
      call s_write ( iounit, output )
      return
    end if
 
    istop3 = istop3 * idtop3
    isbot3 = isbot3 * idbot3
 
  end if
!
!  16. Form the equivalent MATMAN ADD command.
!
  if ( irow1 == irow2 ) then

    command = 'ROW_ADD'

    if ( isbot3 < 0 ) then
      isbot3 = - isbot3
      istop3 = - istop3
    end if

    call i_to_s_left ( irow3, chrtmp1 )
    call i_to_s_left ( irow1, chrtmp2 )
    call rat_to_s_left ( istop3, isbot3, chrtmp )
    line2 =  chrtmp // ' ' // chrtmp1 // ' ' // chrtmp2
    call s_blanks_delete ( line2 )

  else if ( irow1 == irow3 ) then

    command = 'ROW_ADD'

    if ( isbot2 < 0 ) then
      isbot2 = - isbot2
      istop2 = - istop2
    end if

    call rat_to_s_left ( istop2, isbot2, chrtmp )
    call i_to_s_left ( irow2, chrtmp1 )
    call i_to_s_left ( irow1, chrtmp2 )
    line2 = chrtmp // ' ' // chrtmp1 // ' ' // chrtmp2
    call s_blanks_delete ( line2 )

  else

    ierror = 1
    output = 'Your ERO command could not be understood.'
    call s_write ( iounit, output )
    output = 'R2 or R3 must equal R1 in an ERO command.'
    call s_write ( iounit, output )
  end if
 
  return
end
subroutine row_shift ( a, iabot, iatop, irow, maxcol, maxrow, ncol, nrow )

!*****************************************************************************80
!
!! row_shift() allows a new row to be inserted by shifting other rows down.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = rk ) A(MAXROW,MAXCOL), the current matrix.
!
!    Input/output, integer IABOT(MAXROW,MAXCOL), IATOP(MAXROW,MAXCOL).
!    IATOP and IABOT represent the current rational or decimal matrix.
!
!    Input, integer IROW, the position of the new row.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
!    Input, integer NCOL, the number of columns in the matrix.
!
!    Input, integer NROW, the number of rows in the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxcol
  integer maxrow

  real ( kind = rk ) a(maxrow,maxcol)
  integer i
  integer iabot(maxrow,maxcol)
  integer iatop(maxrow,maxcol)
  integer irow
  integer j
  integer ncol
  integer nrow

  do i = nrow, irow + 1, -1
    do j = 1, ncol
 
      a(i,j) = a(i-1,j)
      iatop(i,j) = iatop(i-1,j)
      iabot(i,j) = iabot(i-1,j)
 
    end do
  end do
 
  return
end
subroutine row_swap ( a, iatop, iabot, ibase, ierror, iform, iounit, irow1, &
  irow2, lpmoda, maxcol, maxrow, ncol, nrow )

!*****************************************************************************80
!
!! row_swap() swaps two rows of a matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = rk ) A(MAXROW,MAXCOL), the current matrix.
!
!    Input/output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
!    IATOP and IABOT represent the current rational or decimal matrix.
!
!    Input/output, integer IBASE(MAXROW), the basic variables.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred.
!
!    Input, integer IFORM, 0/1/2, rational/real/decimal arithmetic.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Input, integer IROW1, IROW2, the numbers of the two rows
!    to be swapped.
!
!    Input, integer LPMODA, 0/1, linear algebra/linear programming mode.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
!    Input, integer NCOL, the number of columns in the matrix.
!
!    Input, integer NROW, the number of rows in the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxcol
  integer maxrow

  real ( kind = rk ) a(maxrow,maxcol)
  character ( len = 10 ) chrtmp1
  character ( len = 10 ) chrtmp2
  integer iabot(maxrow,maxcol)
  integer iatop(maxrow,maxcol)
  integer ibase(maxrow)
  integer ierror
  integer iform
  integer iounit(4)
  integer irow1
  integer irow2
  integer j
  integer len1
  integer len2
  integer lpmoda
  integer ncol
  integer nrow
  character ( len = 255 ) output
!
!  Skip out if the two rows are the same.
!
  if ( irow1 == irow2 ) then
    output = 'You have asked to swap a row with itself!'
    call s_write ( iounit, output )
    return
  end if
!
!  Refuse to continue if a row is out of bounds.
!
  if ( ( irow1 < 1 .or. nrow < irow1 ) .or. &
       ( irow2 < 1 .or. nrow < irow2 ) ) then
    ierror = 1
    output = 'One of the rows is illegal!'
    call s_write ( iounit, output )
    return
  end if
!
!  Refuse to swap the last row in linear programming mode.
!
  if ( lpmoda == 1 ) then
    if ( irow1 == nrow .or. irow2 == nrow ) then
      ierror = 1
      output = 'You are in linear programming mode.'
      call s_write ( iounit, output )
      output = 'You may not swap the last row!'
      call s_write ( iounit, output )
      return
    end if
  end if
!
!  Swap the rows.
!
  do j = 1, ncol
 
    if ( iform == 0 ) then
 
      call i_swap ( iatop(irow1,j), iatop(irow2,j) )
      call i_swap ( iabot(irow1,j), iabot(irow2,j) )
 
    else if ( iform == 1 ) then
 
      call r_swap ( a(irow1,j), a(irow2,j) )
 
    else if ( iform == 2 ) then
 
      call i_swap ( iatop(irow1,j), iatop(irow2,j) )
      call i_swap ( iabot(irow1,j), iabot(irow2,j) )
 
    end if
 
  end do
 
  call i_swap ( ibase(irow1), ibase(irow2) )
 
  call i_to_s_left ( irow1, chrtmp1 )
  len1 = len_trim ( chrtmp1 )

  call i_to_s_left ( irow2, chrtmp2 )
  len2 = len_trim ( chrtmp2 )

  output = '  ERO: Row ' // chrtmp1(1:len1) // ' <=> Row ' // chrtmp2(1:len2)

  call s_write ( iounit, output )
 
  return
end
subroutine rvec_house_column ( n, a, k, v )

!*****************************************************************************80
!
!! rvec_house_column() defines a Householder premultiplier that "packs" a column.
!
!  Discussion:
!
!    The routine returns a vector V that defines a Householder
!    premultiplier matrix H(V) that zeros out the subdiagonal entries of
!    column K of the matrix A.
!
!       H(V) = I - 2 * v * v'
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix A.
!
!    Input, real ( kind = rk ) A(N), column K of the matrix A.
!
!    Input, integer K, the column of the matrix to be modified.
!
!    Output, real ( kind = rk ) V(N), a vector of unit L2 norm which defines an
!    orthogonal Householder premultiplier matrix H with the property
!    that the K-th column of H*A is zero below the diagonal.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) alpha
  real ( kind = rk ) column_norm
  integer k
  real ( kind = rk ) rvec_norm2
  real ( kind = rk ) v(n)

  v(1:n) = 0.0D+00

  if ( k < 1 .or. n <= k ) then
    return
  end if

  column_norm = rvec_norm2 ( n + 1 - k, a(k) )

  if ( column_norm == 0.0D+00 ) then
    return
  end if

  alpha = - sign ( 1.0D+00, a(k) ) * column_norm

  v(k) = sqrt ( 0.5D+00 * ( 1.0D+00 - a(k) / alpha ) )
  v(k+1:n) = - 0.5D+00 * a(k+1:n) / ( alpha * v(k) )

  return
end
function rvec_norm2 ( n, a )

!*****************************************************************************80
!
!! rvec_norm2() returns the 2-norm of a vector.
!
!  Discussion:
!
!    The vector 2-norm is defined as:
!
!      RVEC_NORM2 = Sqrt ( Sum ( 1 <= I <= N ) A(I)^2 ).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in A.
!
!    Input, real ( kind = rk ) A(N), the vector.
!
!    Output, real ( kind = rk ) RVEC_NORM2, the 2-norm of A.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  integer i
  real ( kind = rk ) rvec_norm2

  rvec_norm2 = 0.0D+00

  do i = 1, n
    rvec_norm2 = rvec_norm2 + a(i) ** 2
  end do

  rvec_norm2 = sqrt ( rvec_norm2 )

  return
end
subroutine s_blank_delete ( s )

!*****************************************************************************80
!
!! s_blank_delete() removes blanks from a string, left justifying the remainder.
!
!  Discussion:
!
!    All TAB characters are also removed.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
  character c
  integer iget
  integer iput
  character ( len = * ) s
  character TAB
!
  TAB = char ( 9 )
  iput = 0

  do iget = 1, len ( s )

    c = s(iget:iget)

    if ( c /= ' ' .and. c /= TAB ) then
      iput = iput + 1
      s(iput:iput) = c
    end if

  end do

  s(iput+1:) = ' '

  return
end
subroutine s_blanks_delete ( string )

!*****************************************************************************80
!
!! s_blanks_delete() replaces consecutive blanks by one blank.
!
!  Discussion:
!
!    The remaining characters are left justified and right padded with blanks.
!    TAB characters are converted to spaces.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) STRING, the string to be transformed.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer j
  integer nchar
  character newchr
  character oldchr
  character ( len = * ) string
  character TAB

  nchar = len ( string )
  TAB = char ( 9 )
  j = 0
  newchr = ' '

  do i = 1, nchar

    oldchr = newchr
    newchr = string(i:i)

    if ( newchr == TAB ) then
      newchr = ' '
    end if

    string(i:i) = ' '

    if ( oldchr /= ' ' .or. newchr /= ' ' ) then
      j = j + 1
      string(j:j) = newchr
    end if

  end do

  return
end
subroutine s_cap ( string )

!*****************************************************************************80
!
!! s_cap() replaces any lowercase letters by uppercase ones in a string.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) STRING, the string to be transformed.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character c
  integer i
  integer nchar
  character ( len = * ) string

  nchar = len ( string )

  do i = 1, nchar

    c = string(i:i)
    call c_cap ( c )
    string(i:i) = c

  end do

  return
end
subroutine s_chop ( s, ilo, ihi )

!*****************************************************************************80
!
!! s_chop() "chops out" a portion of a string, and closes up the hole.
!
!  Discussion:
!
!    S = 'Fred is not a jerk!'
!
!    call s_chop ( S, 9, 12 )
!
!    S = 'Fred is a jerk!    '
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
!    Input, integer ILO, IHI, the locations of the first and last
!    characters to be removed.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ihi
  integer ihi2
  integer ilo
  integer ilo2
  integer lens
  character ( len = * ) s

  lens = len ( s )

  ilo2 = max ( ilo, 1 )
  ihi2 = min ( ihi, lens )

  if ( ihi2 < ilo2 ) then
    return
  end if

  s(ilo2:lens+ilo2-ihi2-1) = s(ihi2+1:lens)
  s(lens+ilo2-ihi2:lens) = ' '

  return
end
function s_eqi ( strng1, strng2 )

!*****************************************************************************80
!
!! s_eqi() is a case insensitive comparison of two strings for equality.
!
!  Discussion:
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

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer len1
  integer len2
  integer lenc
  logical s_eqi
  character s1
  character s2
  character ( len = * ) strng1
  character ( len = * ) strng2

  len1 = len ( strng1 )
  len2 = len ( strng2 )
  lenc = min ( len1, len2 )

  s_eqi = .false.

  do i = 1, lenc

    s1 = strng1(i:i)
    s2 = strng2(i:i)
    call c_cap ( s1 )
    call c_cap ( s2 )

    if ( s1 /= s2 ) then
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
function s_indexi ( s, sub )

!*****************************************************************************80
!
!! s_indexi() is a case-insensitive INDEX function.
!
!  Discussion:
!
!    The function returns the location in the string at which the
!    substring SUB is first found, or 0 if the substring does not
!    occur at all.
!
!    The routine is also trailing blank insensitive.  This is very
!    important for those cases where you have stored information in
!    larger variables.  If S is of length 80, and SUB is of
!    length 80, then if S = 'FRED' and SUB = 'RED', a match would
!    not be reported by the standard FORTRAN INDEX, because it treats
!    both variables as being 80 characters long!  This routine assumes that
!    trailing blanks represent garbage!
!
!    Because of the suppression of trailing blanks, this routine cannot be
!    used to find, say, the first occurrence of the two-character
!    string 'A '.  However, this routine treats as a special case the
!    occurrence where S or SUB is entirely blank.  Thus you can
!    use this routine to search for occurrences of double or triple blanks
!    in a string, for example, although INDEX itself would be just as
!    suitable for that problem.
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
!    Input, character ( len = * ) S, the string to be searched.
!
!    Input, character ( len = * ) SUB, the substring to search for.
!
!    Output, integer S_INDEXI.  0 if SUB does not occur in
!    the string.  Otherwise S(S_INDEXI:S_INDEXI+LENS-1) = SUB,
!    where LENS is the length of SUB, and is the first place
!    this happens.  S_INDEXI ignores case,
!    unlike the standard FORTRAN INDEX function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer llen1
  integer llen2
  character ( len = * ) s
  logical s_eqi
  integer s_indexi
  character ( len = * ) sub

  s_indexi = 0

  llen1 = len_trim ( s )
  llen2 = len_trim ( sub )
!
!  In case S or SUB is blanks, use LEN.
!
  if ( llen1 == 0 ) then
    llen1 = len ( s )
  end if

  if ( llen2 == 0 ) then
    llen2 = len ( sub )
  end if

  if ( llen1 < llen2 ) then
    return
  end if

  do i = 1, llen1 + 1 - llen2

    if ( s_eqi ( s(i:i+llen2-1), sub ) ) then
      s_indexi = i
      return
    end if

  end do

  return
end
subroutine s_read ( string, line, prompt, iounit, ierror, iterm )

!*****************************************************************************80
!
!! s_read() extracts a character string from the input buffer.
!
!  Discussion:
!
!    The routine accepts an input LINE and a PROMPT line.
!
!    If the LINE is empty, the PROMPT is printed and user input
!    read from IOUNIT(1) into LINE.
!
!    In either case, enough characters are read from LINE to fill
!    STRING and the positions read are removed.
!
!    PROMPT is also updated.  On satisfactory input of STRING,
!    everything in PROMPT up to and including the first comma is removed.
!
!    IOUNIT is assumed to have the following properties, which
!    also apply to routines S_WRITE, CHRINP, R_READ, RELWRT, I_READ
!    and RATREA:
!
!    IOUNIT(1) represents the input unit.  0 is taken to be the user
!    and we READ(*,format) the input.
!
!    IOUNIT(2) is taken to be a standard output unit.  Input is never
!    echoed to IOUNIT(2), but may be to other units.
!
!    Later units:  If their values is between 30 and 39, user input
!    is copied to them, but no output.
!    If between 40 and 49, output is copied to them, but no input.
!    If the unit number is negative, no input is read, nor output
!    written.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) STRING.
!    The user's response to the PROMPT, as read from LINE.
!
!    Input/output, character ( len = 255 ) LINE.
!    A buffer containing the user's input.
!
!    Input/output, character ( len = 255 ) PROMPT.
!    On input, a prompt string that will be printed if necessary.
!
!    On output, if STRING has been read, then PROMPT is cleared out
!    up to, and including, the first comma.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Output, integer IERROR, 
!    0, No error occurred.
!    1, Format error during read.
!    2, End of file during read.
!
!    Input, integer ITERM,
!    0, No check for terminators.
!    1, Blank, slash, comma, semicolon, equals, greater or
!       lesser signs terminate input.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character c
  integer i
  integer ierror
  integer iounit(4)
  integer iterm
  integer lchar
  character ( len = 255 ) line
  integer nchar
  character null
  character ( len = 255 ) prompt
  character ( len = * ) string

  ierror = 0
  null = char(0)
  string = ' '

  call chrinp ( ierror, iounit, line, prompt )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Null input acceptable for character input only.
!
  if ( line == ' ' ) then
    return
  end if

  lchar = 0
  nchar = len ( string )
 
  do i = 1, nchar
 
    if ( lchar /= 0 ) then
      go to 10
    end if
 
    c = line(i:i)
 
    if ( iterm == 1 ) then
 
      if ( c == ' ' .or. c == null .or. c == '/' .or. c == ',' .or. &
           c == ';' .or. c == '=' ) then
        lchar = i
      end if
 
    end if
 
    if ( lchar == 0 ) then
      string(i:i) = c
    end if
 
  end do
 
10    continue
!
!  Chop out the character positions that have been used.
!
  if ( lchar == 0 ) then
    lchar = nchar
  end if

  call s_chop ( line, 1, lchar )
!
!  Force the string to be flush left by removing leading blanks.
!
  line = adjustl ( line )

  return
end
subroutine s_to_dec ( s, itop, ibot, length )

!*****************************************************************************80
!
!! s_to_dec() reads a number from a string, returning a decimal result.
!
!  Discussion:
!
!    The integer may be in real format, for example '2.25'.  It
!    returns ITOP and IBOT.  If the input number is an integer, ITOP
!    equals that integer, and IBOT is 1.  But in the case of 2.25,
!    the program would return ITOP = 225, IBOT = 100.
!
!    Legal input is
!
!          blanks,
!       2  initial sign,
!          blanks,
!       3  whole number,
!       4  decimal point,
!       5  fraction,
!       6  'E' or 'e' or 'D' or 'd', exponent marker,
!       7  exponent sign,
!       8  exponent,
!          blanks
!       9  nonblank
!
!    with most quantities optional.
!
!  Examples:
!
!    S                 ITOP      IBOT     Meaning
!
!    '1'                  1         0        1
!    '     1   '          1         0        1
!    '1A'                 1         0        1
!    '12,34,56'          12         0       12
!    '  34 7'            34         0       34
!    '-1E2ABCD'          -1         2     -100
!    '-1X2ABCD'          -1         0       -1
!    ' 2E-1'              2        -1        0.2
!    '23.45'           2345        -2       23.45
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string containing the
!    data to be read.  Reading will begin at position 1 and
!    terminate when no more characters
!    can be read to form a legal integer.  Blanks, commas,
!    or other nonnumeric data will, in particular, cause
!    the conversion to halt.
!
!    Output, integer ITOP, the integer read from the string,
!    assuming that no negative exponents or fractional parts
!    were used.  Otherwise, the 'integer' is ITOP/IBOT.
!
!    Output, integer IBOT, the integer divisor required to
!    represent numbers which are in real format or have a
!    negative exponent.
!
!    Output, integer LENGTH, number of characters read 
!    to form the number.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  logical c_is_digit
  character c
  integer digit
  integer exponent
  integer exponent_sign
  integer ibot
  integer ihave
  integer iterm
  integer itop
  integer length
  integer mantissa_sign
  character ( len = * ) s
  logical s_eqi

  itop = 0
  ibot = 0

  if ( len ( s ) <= 0 ) then
    length = 0
    return
  end if

  length = - 1
  exponent_sign = 0
  mantissa_sign = 1
  exponent = 0
  ihave = 1
  iterm = 0
!
!  Consider the next character in the string.
!
10    continue

  length = length + 1
  c = s(length+1:length+1)
!
!  Blank.
!
  if ( c == ' ' ) then

    if ( ihave == 1 ) then

    else if ( ihave == 2 ) then

    else 
      iterm = 1
    end if
!
!  Comma or semicolon.
!
  else if ( c == ',' .or. c == ';' ) then

    if ( ihave /= 1 ) then
      iterm = 1
      ihave = 12
      length = length + 1
    end if
!
!  Minus sign.
!
  else if ( c == '-' ) then

    if ( ihave == 1 ) then
      ihave = 2
      mantissa_sign = - 1
    else if ( ihave == 6 ) then
      ihave = 7
      exponent_sign = -1
    else
      iterm = 1
    end if
!
!  Plus sign.
!
  else if ( c == '+' ) then

    if ( ihave == 1 ) then
      ihave = 2
      mantissa_sign = +1
    else if ( ihave == 6 ) then
      ihave = 7
      exponent_sign = +1
    else
      iterm = 1
    end if
!
!  Decimal point.
!
  else if ( c == '.' ) then

    if ( ihave < 4 ) then
      ihave = 4
    else
      iterm = 1
    end if
!
!  Exponent marker.
!
  else if ( s_eqi ( c, 'E' ) .or. s_eqi ( c, 'D' ) ) then

    if ( ihave < 6 ) then
      ihave = 6
    else
      iterm = 1
    end if
!
!  Digit.
!
  else if ( c_is_digit ( c ) ) then

    if ( ihave <= 3 ) then
      ihave = 3
      call c_to_digit ( c, digit )
      itop = 10 * itop + digit
    else if ( ihave <= 5 ) then
      ihave = 5
      call c_to_digit ( c, digit )
      itop = 10 * itop + digit
      ibot = ibot - 1
    else if ( ihave <= 8 ) then
      ihave = 8
      call c_to_digit ( c, digit )
      exponent = 10 * exponent + digit
    else
      write ( *, * ) 'S_TO_DEC: DEBUG: WHAT THE HELL?'
      write ( *, * ) 'IHAVE = ', ihave
      stop
    end if
!
!  Anything else is regarded as a terminator.
!
  else
    iterm = 1
  end if

  if ( iterm == 1 ) then

  else if ( len ( s ) <= length + 1 ) then
    length = len ( s )
  else
    go to 10
  end if
!
!  Number seems to have terminated.
!  Have we got a legal number?
!
  if ( ihave == 1 ) then
    return
  else if ( ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then
    write ( *, * ) ' '
    write ( *, * ) 'S_TO_DEC - Serious error!'
    write ( *, * ) '  Illegal or nonnumeric input:'
    write ( *, '(a)' ) s
    return
  end if
!
!  Normalize.
!
  if ( 0 < itop ) then

    do while ( mod ( itop, 10 ) == 0 )
      itop = itop / 10
      ibot = ibot + 1
    end do

  end if
!
!  Consolidate the number in the form ITOP * 10 ^ IBOT.
!
  ibot = ibot + exponent_sign * exponent
  itop = mantissa_sign * itop

  if ( itop == 0 ) then
   ibot = 0
  end if

  return
end
subroutine s_to_i ( s, ival, ierror, last )

!*****************************************************************************80
!
!! s_to_i() reads an integer value from a string.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 August 1999
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
!    If blank, then IVAL will be returned 0.
!
!    Output, integer IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer LAST, the last character that was
!    part of the representation of IVAL.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character c
  integer i
  integer ierror
  integer isgn
  integer istate
  integer ival
  integer last
  integer lens
  character ( len = * ) s

  ierror = 0
  istate = 0

  isgn = 1
  ival = 0

  lens = len ( s )

  i = 0

10    continue

  i = i + 1

  c = s(i:i)
!
!  ISTATE = 0, haven't read anything except blanks.
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
!  ISTATE = 1, have read a + or minus sign.
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
!  ISTATE = 2, have read at least one digit.
!
  else if ( istate == 2 ) then

    if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
      ival = 10 * ival + ichar ( c ) - ichar ( '0' )
    else
      istate = 3
    end if

  end if
!
!  Continue or exit?
!
  if ( istate == 3 ) then

    ival = isgn * ival
    last = i - 1
    return

  else if ( lens <= i ) then

    if ( istate == 2 ) then
      ival = isgn * ival
      last = lens
    else
      ierror = 1
      last = 0
    end if

    return

  end if

  go to 10
end
subroutine s_to_r ( s, rval, ierror, lchar )

!*****************************************************************************80
!
!! s_to_r() reads a real number from a string.
!
!  Discussion:
!
!    This routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the real number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 spaces
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semicolon.
!
!    with most quantities optional.
!
!  Examples:
!
!    S                RVAL
!
!    '1'               1.0
!    '     1   '       1.0
!    '1A'              1.0
!    '12,34,56'        12.0
!    '  34 7'          34.0
!    '-1E2ABCD'        -100.0
!    '-1X2ABCD'        -1.0
!    ' 2E-1'           0.2
!    '23.45'           23.45
!    '-4.2E+2'         -420.0
!    '17d2'            1700.0
!    '-14e-2'         -0.14
!    'e2'              100.0
!    '-12.73e-9.23'   -12.73 * 10.0 ^ (-9.23)
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
!    Input, character ( len = * ) S, the string containing the
!    data to be read.  Reading will begin at position 1 and
!    terminate at the end of the string, or when no more
!    characters can be read to form a legal real.  Blanks,
!    commas, or other nonnumeric data will, in particular,
!    cause the conversion to halt.
!
!    Output, real ( kind = rk ) RVAL, the real value read from the string.
!
!    Output, integer IERROR, error flag.
!
!    0, no errors occurred.
!
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer LCHAR, the number of characters read from
!    STRING to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  logical c_eqi
  character chrtmp
  integer ierror
  integer ihave
  integer isgn
  integer iterm
  integer jbot
  integer jsgn
  integer jtop
  integer lchar
  integer nchar
  integer ndig
  real ( kind = rk ) rbot
  real ( kind = rk ) rexp
  real ( kind = rk ) rtop
  real ( kind = rk ) rval
  character ( len = * ) s
  character TAB

  nchar = len ( s )
  TAB = char(9)
  ierror = 0
  rval = 0.0D+00
  lchar = - 1
  isgn = 1
  rtop = 0.0D+00
  rbot = 1.0D+00
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

10    continue

  lchar = lchar + 1
  chrtmp = s(lchar+1:lchar+1)
!
!  Blank or TAB character.
!
  if ( chrtmp == ' ' .or. chrtmp == TAB ) then
!
!  20 November 1993
!
!  I would like to allow input like "+ 2", where there is a space
!  between the plus and the number.  So I am going to comment out
!  this line, because I think that's all that's keeping me from
!  doing this.
!
!       if ( ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then
!
    if ( ihave == 2 ) then

    else if ( ihave == 6 .or. ihave == 7 ) then
      iterm = 1
    else if ( 1 < ihave ) then
      ihave = 11
    end if
!
!  Comma.
!
  else if ( chrtmp == ',' .or. chrtmp == ';' ) then

    if ( ihave /= 1 ) then
      iterm = 1
      ihave = 12
      lchar = lchar + 1
    end if
!
!  Minus sign.
!
  else if ( chrtmp == '-' ) then

    if ( ihave == 1 ) then
      ihave = 2
      isgn = - 1
    else if ( ihave == 6 ) then
      ihave = 7
      jsgn = - 1
    else
      iterm = 1
    end if
!
!  Plus sign.
!
  else if ( chrtmp == '+' ) then

    if ( ihave == 1 ) then
      ihave = 2
    else if ( ihave == 6 ) then
      ihave = 7
    else
      iterm = 1
    end if
!
!  Decimal point.
!
  else if ( chrtmp == '.' ) then

    if ( ihave < 4 ) then
      ihave = 4
    else if ( 6 <= ihave .and. ihave <= 8 ) then
      ihave = 9
    else
      iterm = 1
    end if
!
!  Exponent marker.
!
  else if ( c_eqi ( chrtmp, 'E' ) .or. c_eqi ( chrtmp, 'D' ) ) then

    if ( ihave < 6 ) then
      ihave = 6
    else
      iterm = 1
    end if
!
!  Digit.
!
  else if ( ihave < 11 .and. lge ( chrtmp, '0' ) .and. &
    lle ( chrtmp, '9' ) ) then

    if ( ihave <= 2 ) then
      ihave = 3
    else if ( ihave == 4 ) then
      ihave = 5
    else if ( ihave == 6 .or. ihave == 7 ) then
      ihave = 8
    else if ( ihave == 9 ) then
      ihave = 10
    end if

    call c_to_digit ( chrtmp, ndig )

    if ( ihave == 3 ) then
      rtop = 10.0D+00 * rtop + real ( ndig, kind = rk )
    else if ( ihave == 5 ) then
      rtop = 10.0D+00 * rtop + real ( ndig, kind = rk )
      rbot = 10.0D+00 * rbot
    else if ( ihave == 8 ) then
      jtop = 10 * jtop + ndig
    else if ( ihave == 10 ) then
      jtop = 10 * jtop + ndig
      jbot = 10 * jbot
    end if
!
!  Anything else is regarded as a terminator.
!
  else
    iterm = 1
  end if
!
!  If we haven't seen a terminator, and we haven't examined the
!  entire string, go get the next character.
!
  if ( iterm /= 1 .and. lchar+1 < nchar ) then
    go to 10
  end if
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LCHAR is equal to NCHAR.
!
  if ( iterm /= 1 .and. lchar+1 == nchar ) then
    lchar = nchar
  end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
  if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then

    ierror = ihave

    return
  end if
!
!  Number seems OK.  Form it.
!
  if ( jtop == 0 ) then
    rexp = 1.0D+00
  else

    if ( jbot == 1 ) then
      rexp = 10.0D+00 ** ( jsgn * jtop )
    else
      rexp = jsgn * jtop
      rexp = rexp / jbot
      rexp = 10.0D+00 ** rexp
    end if

  end if

  rval = isgn * rexp * rtop / rbot

  return
end
subroutine s_write ( iounit, string )

!*****************************************************************************80
!
!! s_write() writes a character string to one or more output units.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Input, character ( len = * ) STRING, the string to be printed.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer iounit(4)
  integer lchar
  integer npage
  integer output_line_count
  character ( len = * ) string
!
!  If output is to the user, rather than to a file, then
!  see if we need to pause for a new page.
!
  if ( iounit(2) == 0 .and. 0 < npage() ) then
    write ( *, * ) '(more)'
    read ( *, '(a1)', end = 10, err = 10 )
10  continue
  end if
 
  lchar = len_trim ( string )
  if ( lchar <= 0 ) then
    lchar = 1
  end if
 
  do i = 2, 4

    if ( iounit(i) == 0 ) then
!
!  Use the following line for UNIX machines, and the IBM PC:
!
      write ( *, '(a)' ) string(1:lchar)
!
!  Use the following lines for Macintosh and VAX/VMS systems:
!
!         write ( *, '(a,a)' ) ' ', string(1:lchar)
 
    else if ( 0 < iounit(i) ) then
      write ( iounit(i), '(a)' ) string(1:lchar)
    end if
 
  end do
!
!  Update the line count.
!
  output_line_count = 0
  call i_data ( 'INC', 'OUTPUT_LINE_COUNT', output_line_count )
 
  return
end
subroutine sample ( a, chineq, iatop, iabot, ibase, ierror, iform, imat, &
  iounit, line, lpmoda, maxcol, maxrow, nart, ncol, nrow, n_slack, &
  nvar, rmat1 )

!*****************************************************************************80
!
!! sample() allows the user to choose a particular sample problem.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) A(MAXROW,MAXCOL), the current matrix.
!
!    Output, character CHINEQ(MAXROW), the '<', '=', or '>'
!    sign for each linear programming constraint.
!
!    Output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
!    IATOP and IABOT represent the current rational or decimal
!    matrix.
!
!    Output, integer IBASE(MAXROW), the basic variables.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred.
!
!    Output, integer IFORM, 0/1/2, rational/real/decimal arithmetic.
!
!    Input, integer IMAT, 0/1, a matrix HAS NOT/HAS been defined.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Workspace, character ( len = 255 ) LINE, used to hold the user's input.
!
!    Input, integer LPMODA, 0/1, linear algebra/linear programming mode.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
!    Output, integer NART, the number of artificial variables.
!
!    Output, integer NCOL, the number of columns in the matrix.
!
!    Output, integer NROW, the number of rows in the matrix.
!
!    Output, integer N_SLACK, the number of slack variables.
!
!    Output, integer NVAR, the number of basic variables.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxcol
  integer maxrow

  real ( kind = rk ) a(maxrow,maxcol)
  character chineq(maxrow)
  character ( len = 10 ) chrtmp1
  integer i
  integer iabot(maxrow,maxcol)
  integer iatop(maxrow,maxcol)
  integer ibase(maxrow)
  integer ierror
  integer iform
  integer ihi
  integer ilo
  integer imat
  integer iounit(4)
  character isay
  integer iterm
  integer ival
  integer ival2
  integer j
  integer k
  logical s_eqi
  character ( len = 255 ) line
  integer lpmoda
  integer nart
  integer ncol
  integer nrow
  integer n_slack
  integer nvar
  character ( len = 255 ) output
  character ( len = 255 ) prompt
  real ( kind = rk ) rmat1(maxrow,maxcol)

  if ( lpmoda == 0 ) then
 
    output = ' '
    call s_write ( iounit, output )
    output = 'The following examples are available:'
    call s_write ( iounit, output )
    output = '  "E" for eigenvalues;'
    call s_write ( iounit, output )
    output = '  "I" for inverse;'
    call s_write ( iounit, output )
    output = '  "S" for linear solve.'
    call s_write ( iounit, output )
    output = ' '
    call s_write ( iounit, output )
    output = '  "C" to cancel.'
    call s_write ( iounit, output )
    output = ' '
    call s_write ( iounit, output )
 
    prompt = '"D", "E", "I", "S" or "C" to cancel.'
    iterm = 0
    call s_read ( isay, line, prompt, iounit, ierror, iterm )
 
    if ( ierror /= 0 ) then
      return
    end if
!
!  Linear System problem, random square matrix, with RHS vector appended.
!
    if ( s_eqi ( isay, 'S' ) ) then
 
10    continue
 
      prompt = 'number of rows desired.'
 
      call i_read ( nrow, line, prompt, iounit, ierror )
      if ( ierror /= 0 ) then
        return
      end if
 
      ncol = nrow + 1

      if ( nrow < 1 ) then
        output = 'Error!  Negative number of rows not allowed!'
        call s_write ( iounit, output )
        line = ' '
        go to 10
      else if ( maxrow < nrow ) then
        call i_to_s_left ( maxrow, chrtmp1 )
        output = 'Number of rows must be less than ' // chrtmp1
        call s_blanks_delete ( output )
        call s_write ( iounit, output )
        line = ' '
        go to 10
      else if ( maxcol < ncol ) then
        output = 'Please ask for fewer rows NROW, so that '
        call s_write ( iounit, output )
        call i_to_s_left ( maxcol, chrtmp1 )
        output = 'NROW+1 is no more than ' // chrtmp1
        call s_blanks_delete ( output )
        line = ' '
        go to 10
      end if
 
      ilo = -10
      ihi = 10

      do i = 1, nrow
        do j = 1, nrow

          call i_random ( ilo, ihi, ival )

          if ( iform == 0 ) then
            iatop(i,j) = ival
            iabot(i,j) = 1
          else if ( iform == 1 ) then
            a(i,j) = real ( ival, kind = rk )
          else if ( iform == 2 ) then
            iatop(i,j) = ival
            iabot(i,j) = 0
          end if

        end do
      end do

      do i = 1, nrow

        call i_random ( ilo, ihi, ival )
        ibase(i) = ival

      end do

      do i = 1, nrow

        ival = 0
        do j = 1, nrow
          ival = ival + iatop(i,j) * ibase(j)
        end do

        if ( iform == 0 ) then
          iatop(i,nrow+1) = ival
          iabot(i,nrow+1) = 1
        else if ( iform == 1 ) then
          a(i,nrow+1) = real ( ival, kind = rk )
        else if ( iform == 2 ) then
          iatop(i,nrow+1) = ival
          iabot(i,nrow+1) = 0
        end if
 
      end do

      do i = 1, nrow
        ibase(i) = 0
      end do
 
      imat = 1
!
!  Inverse problem, random square matrix, with identity appended.
!
    else if ( s_eqi ( isay, 'I' ) ) then
 
20        continue
 
      nrow = 0
      ncol = 0
 
      prompt = 'number of rows desired.'
 
      call i_read ( nrow, line, prompt, iounit, ierror )
      if ( ierror /= 0 ) then
        return
      end if
 
      if ( nrow < 1 ) then
        output = 'Error!  Negative number of rows not allowed!'
        call s_write ( iounit, output )
        line = ' '
        go to 20
      else if ( maxrow < nrow ) then
        call i_to_s_left ( maxrow, chrtmp1 )
        output = 'Number of rows must be less than ' // chrtmp1
        call s_blanks_delete ( output )
        call s_write ( iounit, output )
        line = ' '
        go to 20
      else if ( maxcol < 2 * nrow ) then
        output = 'Please ask for fewer rows NROW, so that '
        call s_write ( iounit, output )
        call i_to_s_left ( maxcol, chrtmp1 )
        output = '2 * NROW is no more than ' // chrtmp1
        call s_blanks_delete ( output )
        line = ' '
        go to 20
      end if
 
      ncol = 2 * nrow
 
      ilo = -10
      ihi = 10

      do i = 1, nrow
        do j = 1, nrow

          call i_random ( ilo, ihi, ival )

          if ( i == j ) then
            ival2 = 1
          else
            ival2 = 0
          end if

          if ( iform == 0 ) then
            iatop(i,j) = ival
            iabot(i,j) = 1
            iatop(i,j+nrow) = ival2
            iabot(i,j+nrow) = 1
          else if ( iform == 1 ) then
            a(i,j) = real ( ival, kind = rk )
            a(i,j+nrow) = real ( ival2, kind = rk )
          else if ( iform == 2 ) then
            iatop(i,j) = ival
            iabot(i,j) = 0
            iatop(i,j+nrow) = ival2
            iabot(i,j+nrow) = 0
          end if

        end do
      end do
 
      imat = 1
!
!  Eigenvalue sample problem, a random square symmetric matrix.
!
    else if ( s_eqi ( isay, 'E' ) ) then
 
      if ( iform /= 1 ) then
        iform = 1
        output = 'Warning! Automatically switching to REAL arithmetic'
        call s_write ( iounit, output )
      end if

31        continue
 
      nrow = 0
      ncol = 0
 
      prompt = 'number of rows desired.'
 
      call i_read ( nrow, line, prompt, iounit, ierror )
      if ( ierror /= 0 ) then
        return
      end if
 
      if ( nrow < 1 ) then
        output = 'Error!  Negative number of rows not allowed!'
        call s_write ( iounit, output )
        line = ' '
        go to 31
      else if ( maxrow < nrow ) then
        call i_to_s_left ( maxrow, chrtmp1 )
        output = 'Number of rows must be less than ' // chrtmp1
        call s_blanks_delete ( output )
        call s_write ( iounit, output )
        line = ' '
        go to 31
      end if
 
      ncol = nrow
!
!  Get random, integral eigenvalues.
!
      ilo = - 10
      ihi = 10

      do i = 1, nrow
        call i_random ( ilo, ihi, ival )
        iatop(i,1) = ival
      end do
!
!  Get a random orthogonal matrix.
!
      call orth_random ( maxrow, nrow, rmat1 )
!
!  Set A = BT * Lambda * B.
!
      do i = 1, nrow
        do j = 1, nrow
          a(i,j) = 0.0D+00
          do k = 1, nrow
            a(i,j) = a(i,j) + rmat1(k,i) * iatop(k,1) * rmat1(k,j)
          end do
        end do
      end do
 
      imat = 1
 
      output = 'We have a symmetric matrix, whose eigenvalues'
      call s_write ( iounit, output )
      output = 'can be found with the "J" command.'
      call s_write ( iounit, output )

    else
 
      output = 'No problem was selected.'
      call s_write ( iounit, output )
 
    end if
!
!  Linear programming.
!
  else
 
    output = ' '
    call s_write ( iounit, output )
    output = 'The following examples are available:'
    call s_write ( iounit, output )
    output = '  "S" a simple linear programming problem;'
    call s_write ( iounit, output )
    output = '  "A" an advanced linear programming problem.'
    call s_write ( iounit, output )
    output = ' '
    call s_write ( iounit, output )
    output = '  "C" to cancel.'
    call s_write ( iounit, output )
    output = ' '
    call s_write ( iounit, output )
 
    prompt = '"S", "A", or "C" to cancel.'
    iterm = 0
    call s_read ( isay, line, prompt, iounit, ierror, iterm )
 
    if ( ierror /= 0 ) then
      return
    end if
 
    if ( s_eqi ( isay, 'S' ) ) then

      call lp_sams ( a, chineq, iatop, iabot, ibase, iform, imat, &
        iounit, maxcol, maxrow, nart, ncol, nrow, n_slack, nvar )

    else if ( s_eqi ( isay, 'A' ) ) then

      call lp_sama ( a, chineq, iatop, iabot, ibase, iform, imat, &
        iounit, maxcol, maxrow, nart, ncol, nrow, n_slack, nvar ) 

    else
 
      output = 'No problem was selected.'
      call s_write ( iounit, output )
 
    end if
 
  end if
 
  return
end
subroutine timestring ( string )

!*****************************************************************************80
!
!! timestring() writes the current YMDHMS date into a string.
!
!  Example:
!
!    STRING = '31 May 2001   9:45:54.872 AM'
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
!  Output:
!
!    character ( len = * ) STRING, contains the date information.
!    A character length of 40 should always be sufficient.
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
  character ( len = * ) string
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

  write ( string, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end

subroutine transc ( file_tran, ierror, iounit, line )

!*****************************************************************************80
!
!! transc() opens or closes a transcript file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) FILE_TRAN.
!    On input, FILE_TRAN is the current or default transcript file.
!    On output, FILE_TRAN is the file name chosen by the user.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Workspace, character ( len = 255 ) LINE, used to hold the user's input.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = 255 ) file_name
  character ( len = * ) file_tran
  integer ierror
  integer ios
  integer iounit(4)
  integer iosave
  integer iterm
  integer lchar
  character ( len = 255 ) line
  character ( len = 255 ) output
  character ( len = 255 ) prompt
!
!  Get the name of the file.
!
  if ( iounit(3) == -1 ) then

    lchar = len_trim ( file_tran )
    prompt = 'file name, default= "' // file_tran(1:lchar) // '".'
    call s_blanks_delete ( prompt )
    iterm = 0
    call s_read ( file_name, line, prompt, iounit, ierror, iterm )
    if ( ierror /= 0 ) then
      return
    end if

    if ( file_name /= ' ' ) then
      file_tran = file_name
    end if

    call get_unit ( iounit(3) )

    open (  unit = iounit(3), file = file_tran, status = 'replace', &
      form = 'formatted', iostat = ios )

    if ( ios /= 0 ) then
      ierror = 1
      output = 'Unable to open the transcript file!'
      call s_write ( iounit, output )
      return
    end if

    lchar = len_trim ( file_tran )
    output = 'Opening the transcript file "' // file_tran(1:lchar) // '".'
    call s_write ( iounit, output )
 
    iosave = iounit(2)
    iounit(2) = - 1
    call hello ( iounit )
    iounit(2) = iosave
 
  else
 
    lchar = len_trim ( file_tran )
    output = 'Closing the transcript file "' // file_tran(1:lchar) // '".'
    call s_write ( iounit, output )
    close ( unit = iounit(3) )
    iounit(3) = - 1
 
  end if
 
  return
end
subroutine sol_print ( ibase, iform, iounit, islbot, isltop, lpmoda, maxcol, &
  maxrow, nart, ncol, nrow, n_slack, nvar, sol, title )

!*****************************************************************************80
!
!! sol_print() prints out the linear programming solution.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer IBASE(MAXROW), the basic variables.
!
!    Input, integer IFORM, 0/1/2, rational/real/decimal arithmetic.
!
!    Input, integer IOUNIT(4).
!    IOUNIT(1) is the FORTRAN input unit.
!    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
!    IOUNIT(4), if nonzero, are auxilliary output units.
!
!    Output, integer ISLBOT, ISLTOP.
!    Represents the linear programming solution, if fractional
!    or decimal arithmetic is used.
!
!    Input, integer LPMODA, 0/1, linear algebra/linear programming mode.
!
!    Input, integer MAXCOL, the maximum number of matrix columns.
!
!    Input, integer MAXROW, the maximum number of matrix rows.
!
!    Input, integer NART, the number of artificial variables.
!
!    Input, integer NCOL, the number of columns in the matrix.
!
!    Input, integer NROW, the number of rows in the matrix.
!
!    Input, integer N_SLACK, the number of slack variables.
!
!    Input, integer NVAR, the number of basic variables.
!
!    Output, real ( kind = rk ) SOL, the linear programming solution,
!    if real arithmetic is used.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer maxcol
  integer maxrow

  integer ibase(maxrow)
  integer iform
  integer ihi
  integer ilo
  integer iounit(4)
  integer islbot(maxcol)
  integer isltop(maxcol)
  integer jhi
  integer jlo
  integer lpmoda
  integer nart
  integer ncol
  integer nrow
  integer n_slack
  integer nvar
  real ( kind = rk ) sol(maxcol)
  character ( len = 255 ) title

  ilo = 1
  ihi = 1
  jhi = nvar + n_slack + nart
  jlo = 1
 
  if ( iform == 0 ) then

    call rat_print ( isltop, islbot, ibase, iounit, ihi, ilo, jhi, &
      jlo, lpmoda, jhi, 1, ncol, nrow, title )

  else if ( iform == 1 ) then

    call r_print ( sol, ibase, iounit, ihi, ilo, jhi, jlo, lpmoda, &
      jhi, 1, ncol, nrow, title )

  else if ( iform == 2 ) then

    call dec_print ( isltop, islbot, ibase, iounit, ihi, ilo, &
      jhi, jlo, lpmoda, jhi, 1, ncol, nrow, title )

  end if
 
  return
end
