program main

!*****************************************************************************80
!
!! sphere_exactness() investigates the  exactness of a quadrature rule for the sphere.
!
!  Discussion:
!
!    This program investigates the polynomial exactness of a quadrature
!    rule for the unit sphere.
!
!  Usage:
!
!    sphere_exactness files prefix degree_max
!
!    where
!
!    * files explains how the quadrature rule is stored:
!      'XYZW'  for file 'prefix.xyzw' containing (X,Y,Z,Weight);
!      'RTPW'  for file 'prefix.rtpw' containing (Theta, Phi, Weight) (radians);
!      'DTPW'  for file 'prefix.dtpw' containing (Theta, Phi, Weight) (degrees);
!      'XYZ+W' for file 'prefix.xyz' containing (X,Y,Z)
!              and file 'prefix.w' containing Weight;
!      'RTP+W' for file 'prefix.rtp' containing (Theta, Phi ) in radians,
!              and file 'prefix.w' containing Weight;
!      'DTP+W' for file 'prefix.dtp' containing (Theta, Phi ) in degrees,
!              and file 'prefix.w' containing Weight;
!      'XYZ1'  for file 'prefix.xyz' containing (X,Y,Z), 
!              and equal weights, which do not need to be read in.
!      'RTP1'  for file 'prefix.rtp' containing (Theta, Phi ) in radians,
!              and equal weights, which do not need to be read in.
!      'DTP1'  for file 'prefix.dtp' containing (Theta, Phi ) in degrees,'
!              and equal weights, which do not need to be read in.
!    * prefix is the common file prefix;
!    * degree_max is the maximum monomial degree to check.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 September 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer arg_num
  integer degree
  integer degree_max
  integer dim_num
  real ( kind = rk ), allocatable, dimension ( :, : ) :: dtp
  real ( kind = rk ), allocatable, dimension ( :, : ) :: dtpw
  integer, allocatable, dimension ( : ) :: expon
  character ( len = 255 ) filename
  character ( len = 255 ) files
  integer h
  integer iarg
  integer iargc
  integer ierror
  integer last
  logical more
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  integer point_num
  character ( len = 255 ) prefix
  real ( kind = rk ) quad_error
  real ( kind = rk ), allocatable, dimension ( :, : ) :: rtp
  real ( kind = rk ), allocatable, dimension ( :, : ) :: rtpw
  logical s_eqi
  character ( len = 255 ) string
  integer t
  real ( kind = rk ), allocatable, dimension ( : ) :: w
  real ( kind = rk ) w_sum
  real ( kind = rk ), allocatable, dimension ( :, : ) :: xyz
  real ( kind = rk ), allocatable, dimension ( :, : ) :: xyzw

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'sphere_exactness():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Investigate the polynomial exactness of a quadrature'
  write ( *, '(a)' ) '  rule for the unit sphere by integrating all monomials'
  write ( *, '(a)' ) '  of a given degree.'
!
!  Get the number of command line arguments.
!
  arg_num = iargc ( )
!
!  Get the file structure:
!
  if ( 1 <= arg_num ) then

    iarg = 1
    call getarg ( iarg, files )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPHERE_EXACTNESS:'
    write ( *, '(a)' ) '  Describe the files to be read:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  For coordinates and weights in one file:'
    write ( *, '(a)' ) '    XYZW     (X,Y,Z,Weight)'
    write ( *, '(a)' ) '    RTPW     (Theta, Phi, Weight) (radians)'
    write ( *, '(a)' ) '    DTPW     (Theta, Phi, Weight) (degrees)'
    write ( *, '(a)' ) '  For coordinates in one file and weights in another:'
    write ( *, '(a)' ) '    XYZ+W    (X,Y,Z)       + Weight'
    write ( *, '(a)' ) '    RTP+W    (Theta, Phi ) + Weight'
    write ( *, '(a)' ) '    DTP+W    (Theta, Phi ) + Weight'
    write ( *, '(a)' ) '  For coordinates in one file, and equal weights:'
    write ( *, '(a)' ) '    XYZ1     (X,Y,Z)'
    write ( *, '(a)' ) '    RTP1     (Theta, Phi ) (radians)'
    write ( *, '(a)' ) '    DTP1     (Theta, Phi ) (degrees)'

    read ( *, '(a)' ) files

  end if
!
!  Get the file prefix:
!
  if ( 2 <= arg_num ) then

    iarg = 2
    call getarg ( iarg, prefix )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPHERE_EXACTNESS:'
    write ( *, '(a)' ) '  Enter the filename prefix.'

    read ( *, '(a)' ) prefix

  end if
!
!  The third command line argument is the maximum degree.
!
  if ( 3 <= arg_num ) then

    iarg = 3
    call getarg ( iarg, string )
    call s_to_i4 ( string, degree_max, ierror, last )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPHERE_EXACTNESS:'
    write ( *, '(a)' ) '  Please enter the maximum total degree to check.'

    read ( *, * ) degree_max

  end if
!
!  Summarize the input.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPHERE_EXACTNESS: User input:'
  write ( *, '(a)' ) '  File structure = "' // trim ( files ) // '".'
  write ( *, '(a)' ) '  Filename prefix = "' // trim ( prefix ) // '".'
  write ( *, '(a,i8)' ) '  Maximum degree = ', degree_max
!
!  Read data needed to create XYZ and W arrays.
!
  if ( s_eqi ( files, 'xyzw' ) ) then
    
    filename = trim ( prefix ) // '.xyzw'

    call r8mat_header_read ( filename, dim_num, point_num )

    allocate ( xyzw(1:4,1:point_num) )

    call r8mat_data_read ( filename, 4, point_num, xyzw )

    allocate ( xyz(1:3,1:point_num) )
    allocate ( w(1:point_num) )

    xyz(1:3,1:point_num) = xyzw(1:3,1:point_num)
    w(1:point_num) = xyzw(4,1:point_num)

    deallocate ( xyzw )

  else if ( s_eqi ( files, 'xyz+w' ) ) then
    
    filename = trim ( prefix ) // '.xyz'

    call r8mat_header_read ( filename, dim_num, point_num )

    allocate ( xyz(1:3,1:point_num) )

    call r8mat_data_read ( filename, 3, point_num, xyz )

    allocate ( w(1:point_num) )

    filename = trim ( prefix ) // '.w'

    call r8mat_data_read ( filename, 1, point_num, w )

  else if ( s_eqi ( files, 'xyz1' ) ) then
    
    filename = trim ( prefix ) // '.xyz'

    call r8mat_header_read ( filename, dim_num, point_num )

    allocate ( xyz(1:3,1:point_num) )

    call r8mat_data_read ( filename, 3, point_num, xyz )

    allocate ( w(1:point_num) )

    w(1:point_num) = 4.0D+00 * pi / point_num

  else if ( s_eqi ( files, 'rtpw' ) ) then
    
    filename = trim ( prefix ) // '.rtpw'

    call r8mat_header_read ( filename, dim_num, point_num )

    allocate ( rtpw(1:3,1:point_num) )

    call r8mat_data_read ( filename, 3, point_num, rtpw )

    allocate ( xyz(1:3,1:point_num) )

    allocate ( w(1:point_num) )

    xyz(1,1:point_num) = &
      cos ( rtpw(1,1:point_num) ) * sin ( rtpw(2,1:point_num) )
    xyz(2,1:point_num) = &
      sin ( rtpw(1,1:point_num) ) * sin ( rtpw(2,1:point_num) )
    xyz(3,1:point_num) = &
                                    cos ( rtpw(2,1:point_num) )

    w(1:point_num) = rtpw(3,1:point_num)

    deallocate ( rtpw )

  else if ( s_eqi ( files, 'rtp+w' ) ) then
    
    filename = trim ( prefix ) // '.rtp'

    call r8mat_header_read ( filename, dim_num, point_num )

    allocate ( rtp(1:2,1:point_num) )

    call r8mat_data_read ( filename, 2, point_num, rtp )

    filename = trim ( prefix ) // '.w'

    allocate ( w(1:point_num) )

    call r8mat_data_read ( filename, 1, point_num, w )

    allocate ( xyz(1:3,1:point_num) )

    xyz(1,1:point_num) = &
      cos ( rtp(1,1:point_num) ) * sin ( rtp(2,1:point_num) )
    xyz(2,1:point_num) = &
      sin ( rtp(1,1:point_num) ) * sin ( rtp(2,1:point_num) )
    xyz(3,1:point_num) = &
                                   cos ( rtp(2,1:point_num) )

    deallocate ( rtp )

  else if ( s_eqi ( files, 'rtp1' ) ) then
    
    filename = trim ( prefix ) // '.rtp'

    call r8mat_header_read ( filename, dim_num, point_num )

    allocate ( rtp(1:2,1:point_num) )

    call r8mat_data_read ( filename, 2, point_num, rtp )

    allocate ( xyz(1:3,1:point_num) )

    xyz(1,1:point_num) = &
      cos ( rtp(1,1:point_num) ) * sin ( rtp(2,1:point_num) )
    xyz(2,1:point_num) = &
      sin ( rtp(1,1:point_num) ) * sin ( rtp(2,1:point_num) )
    xyz(3,1:point_num) = &
                                   cos ( rtp(2,1:point_num) )

    deallocate ( rtp )

    allocate ( w(1:point_num) )

    w(1:point_num) = 4.0D+00 * pi / point_num

  else if ( s_eqi ( files, 'dtpw' ) ) then
    
    filename = trim ( prefix ) // '.dtpw'

    call r8mat_header_read ( filename, dim_num, point_num )

    allocate ( dtpw(1:3,1:point_num) )

    call r8mat_data_read ( filename, 3, point_num, dtpw )

    allocate ( xyz(1:3,1:point_num) )

    allocate ( w(1:point_num) )

    dtpw(1:2,1:point_num) = dtpw(1:2,1:point_num) * pi / 180.0D+00

    xyz(1,1:point_num) = &
      cos ( dtpw(1,1:point_num) ) * sin ( dtpw(2,1:point_num) )
    xyz(2,1:point_num) = &
      sin ( dtpw(1,1:point_num) ) * sin ( dtpw(2,1:point_num) )
    xyz(3,1:point_num) = &
                                    cos ( dtpw(2,1:point_num) )
    w(1:point_num) = dtpw(3,1:point_num)

    deallocate ( dtpw )

  else if ( s_eqi ( files, 'dtp+w' ) ) then
    
    filename = trim ( prefix ) // '.dtp'

    call r8mat_header_read ( filename, dim_num, point_num )

    allocate ( dtp(1:2,1:point_num) )

    call r8mat_data_read ( filename, 2, point_num, dtp )

    filename = trim ( prefix ) // '.w'

    allocate ( w(1:point_num) )

    call r8mat_data_read ( filename, 1, point_num, w )

    allocate ( xyz(1:3,1:point_num) )

    dtp(1:2,1:point_num) = dtp(1:2,1:point_num) * pi / 180.0D+00

    xyz(1,1:point_num) = &
      cos ( dtp(1,1:point_num) ) * sin ( dtp(2,1:point_num) )
    xyz(2,1:point_num) = &
      sin ( dtp(1,1:point_num) ) * sin ( dtp(2,1:point_num) )
    xyz(3,1:point_num) = &
                                   cos ( dtp(2,1:point_num) )

    deallocate ( dtp )

  else if ( s_eqi ( files, 'dtp1' ) ) then
    
    filename = trim ( prefix ) // '.dtp'

    call r8mat_header_read ( filename, dim_num, point_num )

    allocate ( dtp(1:2,1:point_num) )

    call r8mat_data_read ( filename, 2, point_num, dtp )

    allocate ( xyz(1:3,1:point_num) )

    dtp(1:2,1:point_num) = dtp(1:2,1:point_num) * pi / 180.0D+00

    xyz(1,1:point_num) = &
      cos ( dtp(1,1:point_num) ) * sin ( dtp(2,1:point_num) )
    xyz(2,1:point_num) = &
      sin ( dtp(1,1:point_num) ) * sin ( dtp(2,1:point_num) )
    xyz(3,1:point_num) = &
                                   cos ( dtp(2,1:point_num) )

    deallocate ( dtp )

    allocate ( w(1:point_num) )

    w(1:point_num) = 4.0D+00 * pi / point_num

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPHERE_EXACTNESS - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized file structure choice!'
    stop

  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of points  = ', point_num
!
!  The W's should sum to 4 * PI.
!
  w_sum = sum ( w(1:point_num) )
  
  w(1:point_num) = 4.0D+00 * pi * w(1:point_num) / w_sum
!
!  Explore the monomials.
!
  allocate ( expon(1:3) )
  expon(1:3) = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          Error          Degree  Exponents'

  do degree = 0, degree_max

    write ( *, '(a)' ) ' '

    more = .false.
    h = 0
    t = 0

    do

      call comp_next ( degree, 3, expon, more, h, t )

      call sphere01_monomial_quadrature ( expon, point_num, xyz, w, quad_error )

      write ( *, '(2x,f24.16,3x,i2,4x,10i3)' ) &
        quad_error, degree, expon(1:3)

      if ( .not. more ) then
        exit
      end if

    end do

  end do
!
!  Free memory.
!
  deallocate ( expon )
  deallocate ( w )
  deallocate ( xyz )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPHERE_EXACTNESS:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
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
function ch_eqi ( c1, c2 )

!*****************************************************************************80
!
!! CH_EQI is a case insensitive comparison of two characters for equality.
!
!  Example:
!
!    CH_EQI ( 'A', 'a' ) is .TRUE.
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
!    Input, character C1, C2, the characters to compare.
!
!    Output, logical CH_EQI, the result of the comparison.
!
  implicit none

  logical ch_eqi
  character c1
  character c1_cap
  character c2
  character c2_cap

  c1_cap = c1
  c2_cap = c2

  call ch_cap ( c1_cap )
  call ch_cap ( c2_cap )

  if ( c1_cap == c2_cap ) then
    ch_eqi = .true.
  else
    ch_eqi = .false.
  end if

  return
end
subroutine ch_to_digit ( c, digit )

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the value of a base 10 digit.
!
!  Example:
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
!    Output, integer DIGIT, the corresponding value.  
!    If C was 'illegal', then DIGIT is -1.
!
  implicit none

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
subroutine comp_next ( n, k, a, more, h, t )

!*****************************************************************************80
!
!! COMP_NEXT computes the compositions of the integer N into K parts.
!
!  Discussion:
!
!    A composition of the integer N into K parts is an ordered sequence
!    of K nonnegative integers which sum to N.  The compositions (1,2,1)
!    and (1,1,2) are considered to be distinct.
!
!    The routine computes one composition on each call until there are no more.
!    For instance, one composition of 6 into 3 parts is
!    3+2+1, another would be 6+0+0.
!
!    On the first call to this routine, set MORE = FALSE.  The routine
!    will compute the first element in the sequence of compositions, and
!    return it, as well as setting MORE = TRUE.  If more compositions
!    are desired, call again, and again.  Each time, the routine will
!    return with a new composition.
!
!    However, when the LAST composition in the sequence is computed 
!    and returned, the routine will reset MORE to FALSE, signaling that
!    the end of the sequence has been reached.
!
!    This routine originally used a SAVE statement to maintain the
!    variables H and T.  I have decided that it is safer
!    to pass these variables as arguments, even though the user should
!    never alter them.  This allows this routine to safely shuffle
!    between several ongoing calculations.
!
!
!    There are 28 compositions of 6 into three parts.  This routine will
!    produce those compositions in the following order:
!
!     I         A
!     -     ---------
!     1     6   0   0
!     2     5   1   0
!     3     4   2   0
!     4     3   3   0
!     5     2   4   0
!     6     1   5   0
!     7     0   6   0
!     8     5   0   1
!     9     4   1   1
!    10     3   2   1
!    11     2   3   1
!    12     1   4   1
!    13     0   5   1
!    14     4   0   2
!    15     3   1   2
!    16     2   2   2
!    17     1   3   2
!    18     0   4   2
!    19     3   0   3
!    20     2   1   3
!    21     1   2   3
!    22     0   3   3
!    23     2   0   4
!    24     1   1   4
!    25     0   2   4
!    26     1   0   5
!    27     0   1   5
!    28     0   0   6
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 July 2008
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Second Edition,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer N, the integer whose compositions are desired.
!
!    Input, integer K, the number of parts in the composition.
!
!    Input/output, integer A(K), the parts of the composition.
!
!    Input/output, logical MORE, set by the user to start the computation,
!    and by the routine to terminate it.
!
!    Input/output, integer  H, T, two internal parameters needed 
!    for the computation.  The user should allocate space for these in the 
!    calling program, include them in the calling sequence, but never alter 
!    them!
!
  implicit none

  integer k

  integer a(k)
  integer h
  logical more
  integer n
  integer t
!
!  The first computation.
!
  if ( .not. more ) then

    t = n
    h = 0
    a(1) = n
    a(2:k) = 0
!
!  The next computation.
!
  else

    if ( 1 < t ) then
      h = 0
    end if

    h = h + 1
    t = a(h)
    a(h) = 0
    a(1) = t - 1
    a(h+1) = a(h+1) + 1

  end if
!
!  This is the last element of the sequence if all the
!  items are in the last slot.
!
  more = ( a(k) /= n )

  return
end
subroutine file_column_count ( input_file_name, column_num )

!*****************************************************************************80
!
!! FILE_COLUMN_COUNT counts the number of columns in the first line of a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.
!
!    Most lines of the file is presumed to consist of COLUMN_NUM words,
!    separated by spaces.  There may also be some blank lines, and some 
!    comment lines,
!    which have a "#" in column 1.
!
!    The routine tries to find the first non-comment non-blank line and
!    counts the number of words in that line.
!
!    If all lines are blanks or comments, it goes back and tries to analyze
!    a comment line.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILE_NAME, the name of the file.
!
!    Output, integer COLUMN_NUM, the number of columns in the file.
!
  implicit none

  integer column_num
  logical got_one
  character ( len = * ) input_file_name
  integer input_status
  integer input_unit
  character ( len = 255 ) line
!
!  Open the file.
!
  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_file_name, status = 'old', &
    form = 'formatted', access = 'sequential', iostat = input_status )

  if ( input_status /= 0 ) then
    column_num = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' &
      // trim ( input_file_name ) // '" on unit ', input_unit
    return
  end if
!
!  Read one line, but skip blank lines and comment lines.
!
  got_one = .false.

  do

    read ( input_unit, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
      exit
    end if

    if ( len_trim ( line ) == 0 ) then
      cycle
    end if

    if ( line(1:1) == '#' ) then
      cycle
    end if

    got_one = .true.
    exit

  end do

  if ( .not. got_one ) then

    rewind ( input_unit )

    do

      read ( input_unit, '(a)', iostat = input_status ) line

      if ( input_status /= 0 ) then
        exit
      end if

      if ( len_trim ( line ) == 0 ) then
        cycle
      end if

      got_one = .true.
      exit

    end do

  end if

  close ( unit = input_unit )

  if ( .not. got_one ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Warning!'
    write ( *, '(a)' ) '  The file does not seem to contain any data.'
    column_num = -1
    return
  end if

  call s_word_count ( line, column_num )

  return
end
subroutine file_row_count ( input_file_name, row_num )

!*****************************************************************************80
!
!! FILE_ROW_COUNT counts the number of row records in a file.
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
!    06 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILE_NAME, the name of the input file.
!
!    Output, integer ROW_NUM, the number of rows found.
!
  implicit none

  integer bad_num
  integer comment_num
  integer ierror
  character ( len = * ) input_file_name
  integer input_status
  integer input_unit
  character ( len = 255 ) line
  integer record_num
  integer row_num

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_file_name, status = 'old', &
    iostat = input_status )

  if ( input_status /= 0 ) then
    row_num = -1;
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_ROW_COUNT - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' // &
      trim ( input_file_name ) // '" on unit ', input_unit
    stop
  end if

  comment_num = 0
  row_num = 0
  record_num = 0
  bad_num = 0

  do

    read ( input_unit, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
      ierror = record_num
      exit
    end if

    record_num = record_num + 1

    if ( line(1:1) == '#' ) then
      comment_num = comment_num + 1
      cycle
    end if

    if ( len_trim ( line ) == 0 ) then
      comment_num = comment_num + 1
      cycle
    end if

    row_num = row_num + 1

  end do

  close ( unit = input_unit )

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IUNIT, the free unit number.
!
  implicit none

  integer i
  integer ios
  integer iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

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
subroutine monomial_value ( m, n, e, x, v )

!*****************************************************************************80
!
!! MONOMIAL_VALUE evaluates a monomial.
!
!  Discussion:
!
!    This routine evaluates a monomial of the form
!
!      product ( 1 <= i <= m ) x(i)^e(i)
!
!    The combination 0.0^0 is encountered is treated as 1.0.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 April 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of evaluation points.
!
!    Input, integer E(M), the exponents.
!
!    Input, real ( kind = rk ) X(M,N), the point coordinates.
!
!    Output, real ( kind = rk ) V(N), the monomial values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer e(m)
  integer i
  real ( kind = rk ) v(n)
  real ( kind = rk ) x(m,n)

  v(1:n) = 1.0D+00

  do i = 1, m
    if ( 0 /= e(i) ) then
      v(1:n) = v(1:n) * x(i,1:n) ** e(i)
    end if
  end do

  return
end
subroutine r8mat_data_read ( input_filename, m, n, table )

!*****************************************************************************80
!
!! R8MAT_DATA_READ reads data from an R8MAT file.
!
!  Discussion:
!
!    The file may contain more than N points, but this routine will
!    return after reading N of them.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of points.
!
!    Output, real ( kind = rk ) TABLE(M,N), the table data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer ierror
  character ( len = * ) input_filename
  integer input_status
  integer input_unit
  integer j
  character ( len = 255 ) line
  real ( kind = rk ) table(m,n)
  real ( kind = rk ) x(m)

  ierror = 0

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = input_status )

  if ( input_status /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_DATA_READ - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' // &
      trim ( input_filename ) // '" on unit ', input_unit
    stop
  end if

  j = 0

  do while ( j < n )

    read ( input_unit, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8MAT_DATA_READ - Fatal error!'
      write ( *, '(a)' ) '  Error while reading lines of data.'
      write ( *, '(a,i8)' ) '  Number of values expected per line M = ', m
      write ( *, '(a,i8)' ) '  Number of data lines read, J =         ', j
      write ( *, '(a,i8)' ) '  Number of data lines needed, N =       ', n
      stop
    end if

    if ( line(1:1) == '#' .or. len_trim ( line ) == 0 ) then
      cycle
    end if

    call s_to_r8vec ( line, m, x, ierror )

    if ( ierror /= 0 ) then
      cycle
    end if

    j = j + 1

    table(1:m,j) = x(1:m)

  end do

  close ( unit = input_unit )

  return
end
subroutine r8mat_header_read ( input_filename, m, n )

!*****************************************************************************80
!
!! R8MAT_HEADER_READ reads the header from an R8MAT file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Output, integer M, spatial dimension.
!
!    Output, integer N, the number of points.
!
  implicit none

  character ( len = * ) input_filename
  integer m
  integer n

  call file_column_count ( input_filename, m )

  if ( m <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data columns in'
    write ( *, '(a)' ) '  the file "' // trim ( input_filename ) // '".'
    stop
  end if

  call file_row_count ( input_filename, n )

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data rows in'
    write ( *, '(a)' ) '  the file "' // trim ( input_filename ) // '".'
    stop
  end if

  return
end
function s_eqi ( s1, s2 )

!*****************************************************************************80
!
!! S_EQI is a case insensitive comparison of two strings for equality.  
!
!  Discussion:
!
!    S_EQI ( 'Anjana', 'ANJANA' ) is TRUE.
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
!    Input, character ( len = * ) S1, S2, the strings to compare.
!
!    Output, logical S_EQI, the result of the comparison.
!
  implicit none

  character c1
  character c2
  integer i
  integer lenc
  logical s_eqi
  character ( len = * ) s1
  integer s1_length
  character ( len = * ) s2
  integer s2_length

  s1_length = len ( s1 )
  s2_length = len ( s2 )
  lenc = min ( s1_length, s2_length )
 
  s_eqi = .false.

  do i = 1, lenc

    c1 = s1(i:i)
    c2 = s2(i:i)
    call ch_cap ( c1 )
    call ch_cap ( c2 )

    if ( c1 /= c2 ) then
      return
    end if

  end do
 
  do i = lenc + 1, s1_length
    if ( s1(i:i) /= ' ' ) then
      return
    end if
  end do
 
  do i = lenc + 1, s2_length
    if ( s2(i:i) /= ' ' ) then
      return
    end if
  end do
 
  s_eqi = .true.
 
  return
end
subroutine s_to_i4 ( s, ival, ierror, length )

!*****************************************************************************80
!
!! S_TO_I4 reads an I4 from a string.
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
!    Output, integer IVAL, the value read from the string.
!    If the string is blank, then IVAL will be returned 0.
!
!    Output, integer IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer LENGTH, the number of characters of S 
!    used to make IVAL.
!
  implicit none

  character c
  integer i
  integer ierror
  integer isgn
  integer istate
  integer ival
  integer length
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
        length = i - 1
        return
      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( istate == 2 ) then
    ival = isgn * ival
    length = len_trim ( s )
  else
    ierror = 1
    length = 0
  end if

  return
end
subroutine s_to_r8 ( s, dval, ierror, length )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 from a string.
!
!  Discussion:
!
!    The routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 blanks
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semicolon,
!
!    with most quantities optional.
!
!  Example:
!
!    S                 DVAL
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
!    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 September 2004
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
!    Output, real ( kind = rk ) DVAL, the value read from the string.
!
!    Output, integer IERROR, error flag.
!    0, no errors occurred.
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer LENGTH, the number of characters read
!    to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  logical ch_eqi
  character c
  real ( kind = rk ) dval
  integer ierror
  integer ihave
  integer isgn
  integer iterm
  integer jbot
  integer jsgn
  integer jtop
  integer length
  integer nchar
  integer ndig
  real ( kind = rk ) rbot
  real ( kind = rk ) rexp
  real ( kind = rk ) rtop
  character ( len = * ) s

  nchar = len_trim ( s )

  ierror = 0
  dval = 0.0D+00
  length = -1
  isgn = 1
  rtop = 0
  rbot = 1
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

  do

    length = length + 1

    if ( nchar < length+1 ) then
      exit
    end if

    c = s(length+1:length+1)
!
!  Blank character.
!
    if ( c == ' ' ) then

      if ( ihave == 2 ) then

      else if ( ihave == 6 .or. ihave == 7 ) then
        iterm = 1
      else if ( 1 < ihave ) then
        ihave = 11
      end if
!
!  Comma.
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
        isgn = -1
      else if ( ihave == 6 ) then
        ihave = 7
        jsgn = -1
      else
        iterm = 1
      end if
!
!  Plus sign.
!
    else if ( c == '+' ) then

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
    else if ( c == '.' ) then

      if ( ihave < 4 ) then
        ihave = 4
      else if ( 6 <= ihave .and. ihave <= 8 ) then
        ihave = 9
      else
        iterm = 1
      end if
!
!  Scientific notation exponent marker.
!
    else if ( ch_eqi ( c, 'E' ) .or. ch_eqi ( c, 'D' ) ) then

      if ( ihave < 6 ) then
        ihave = 6
      else
        iterm = 1
      end if
!
!  Digit.
!
    else if (  ihave < 11 .and. lle ( '0', c ) .and. lle ( c, '9' ) ) then

      if ( ihave <= 2 ) then
        ihave = 3
      else if ( ihave == 4 ) then
        ihave = 5
      else if ( ihave == 6 .or. ihave == 7 ) then
        ihave = 8
      else if ( ihave == 9 ) then
        ihave = 10
      end if

      call ch_to_digit ( c, ndig )

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
    if ( iterm == 1 ) then
      exit
    end if

  end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LENGTH is equal to NCHAR.
!
  if ( iterm /= 1 .and. length+1 == nchar ) then
    length = nchar
  end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
  if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then
    ierror = ihave
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'S_TO_R8 - Serious error!'
    write ( *, '(a)' ) '  Illegal or nonnumeric input:'
    write ( *, '(a)' ) '    ' // trim ( s )
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
      rexp = 10.0D+00 ** ( real ( jsgn * jtop, kind = rk ) &
        / real ( jbot, kind = rk ) )
    end if
  end if

  dval = real ( isgn, kind = rk ) * rexp * rtop / rbot

  return
end
subroutine s_to_r8vec ( s, n, rvec, ierror )

!*****************************************************************************80
!
!! S_TO_R8VEC reads an R8VEC from a string.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be read.
!
!    Input, integer N, the number of values expected.
!
!    Output, real ( kind = rk ) RVEC(N), the values read from the string.
!
!    Output, integer IERROR, error flag.
!    0, no errors occurred.
!    -K, could not read data for entries -K through N.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer ierror
  integer ilo
  integer lchar
  real ( kind = rk ) rvec(n)
  character ( len = * ) s

  i = 0
  ierror = 0
  ilo = 1

  do while ( i < n )

    i = i + 1

    call s_to_r8 ( s(ilo:), rvec(i), ierror, lchar )

    if ( ierror /= 0 ) then
      ierror = -i
      exit
    end if

    ilo = ilo + lchar

  end do

  return
end
subroutine s_word_count ( s, nword )

!*****************************************************************************80
!
!! S_WORD_COUNT counts the number of "words" in a string.
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
!    Input, character ( len = * ) S, the string to be examined.
!
!    Output, integer NWORD, the number of "words" in the string.
!    Words are presumed to be separated by one or more blanks.
!
  implicit none

  logical blank
  integer i
  integer lens
  integer nword
  character ( len = * ) s

  nword = 0
  lens = len ( s )

  if ( lens <= 0 ) then
    return
  end if

  blank = .true.

  do i = 1, lens

    if ( s(i:i) == ' ' ) then
      blank = .true.
    else if ( blank ) then
      nword = nword + 1
      blank = .false.
    end if

  end do

  return
end
subroutine sphere01_monomial_integral ( e, integral )

!*****************************************************************************80
!
!! SPHERE01_MONOMIAL_INTEGRAL returns monomial integrals on the unit sphere.
!
!  Discussion:
!
!    The integration region is 
!
!      X^2 + Y^2 + Z^2 = 1.
!
!    The monomial is F(X,Y,Z) = X^E(1) * Y^E(2) * Z^E(3).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    24 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Academic Press, 1984, page 263.
!
!  Parameters:
!
!    Input, integer E(3), the exponents of X, Y and Z in the 
!    monomial.  Each exponent must be nonnegative.
!
!    Output, real ( kind = rk ) INTEGRAL, the integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer e(3)
  integer i
  real ( kind = rk ) integral
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00

  if ( any ( e(1:3) < 0 ) ) then
    integral = - huge ( integral )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPHERE01_MONOMIAL_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  All exponents must be nonnegative.'
    write ( *, '(a,i8)' ) '  E(1) = ', e(1)
    write ( *, '(a,i8)' ) '  E(2) = ', e(2)
    write ( *, '(a,i8)' ) '  E(3) = ', e(3)
    stop
  end if

  if ( all ( e(1:3) == 0 ) ) then

    integral = 2.0D+00 * sqrt ( pi**3 ) / gamma ( 1.5D+00 )

  else if ( any ( mod ( e(1:3), 2 ) == 1 ) ) then

    integral = 0.0D+00

  else

    integral = 2.0D+00

    do i = 1, 3
      integral = integral * gamma ( 0.5D+00 * real ( e(i) + 1, kind = rk ) )
    end do

    integral = integral &
      / gamma ( 0.5D+00 * ( real ( sum ( e(1:3) + 1 ), kind = rk ) ) )

  end if

  return
end
subroutine sphere01_monomial_quadrature ( expon, point_num, xyz, w, quad_error )

!*****************************************************************************80
!
!! SPHERE01_MONOMIAL_QUADRATURE applies quadrature to a monomial in a sphere.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer EXPON(3), the exponents.
!
!    Input, integer POINT_NUM, the number of points in the rule.
!
!    Input, real ( kind = rk ) XYZ(3,POINT_NUM), the quadrature points.
!
!    Input, real ( kind = rk ) W(POINT_NUM), the quadrature weights.
!
!    Output, real ( kind = rk ) QUAD_ERROR, the quadrature error.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) exact
  integer expon(3)
  integer point_num
  real ( kind = rk ) quad
  real ( kind = rk ) quad_error
  real ( kind = rk ) value(point_num)
  real ( kind = rk ) w(point_num)
  real ( kind = rk ) xyz(3,point_num)
!
!  Get the exact value of the integral.
!
  call sphere01_monomial_integral ( expon, exact )
!
!  Evaluate the monomial at the quadrature points.
!
  call monomial_value ( 3, point_num, expon, xyz, value )
!
!  Compute the weighted sum.
!
  quad = dot_product ( w, value )

  quad_error = abs ( quad - exact )

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
