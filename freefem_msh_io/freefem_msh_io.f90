subroutine ch_cap ( ch )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
!
!  Discussion:
!
!    Instead of CHAR and ICHAR, we now use the ACHAR and IACHAR functions,
!    which guarantee the ASCII collating sequence.
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
!    Input/output, character CH, the character to capitalize.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ch
  integer itemp

  itemp = iachar ( ch )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    ch = achar ( itemp - 32 )
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

  integer, parameter :: rk = kind ( 1.0D+00 )

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
!! CH_TO_DIGIT returns the integer value of a base 10 digit.
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
!    Output, integer DIGIT, the corresponding integer value.
!    If C was 'illegal', then DIGIT is -1.
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

  integer, parameter :: rk = kind ( 1.0D+00 )

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
subroutine ffmsh_2d_data_example ( v_num, e_num, t_num, v_xy, v_l, e_v, e_l, &
  t_v, t_l )

!*****************************************************************************80
!
!! FFMSH_2D_DATA_EXAMPLE returns example FFMSH data.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 December 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer V_NUM, the number of vertices.
!
!    Input, integer E_NUM, the number of boundary edges.
!
!    Input, integer T_NUM, the number of triangles.
!
!    Output, real ( kind = rk ) V_XY(2,V_NUM), vertex coordinates.
!
!    Output, integer V_L(V_NUM), vertex labels.
!
!    Output, integer E_V(2,E_NUM), edge vertices.
!
!    Output, integer E_L(E_NUM), vertex labels.
!
!    Output, integer T_V(3,T_NUM), triangle vertices.
!
!    Output, integer T_L(T_NUM), triangle labels.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer e_num
  integer t_num
  integer v_num

  integer e_l(e_num)
  integer, dimension (10) :: e_l_save = (/ &
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1 /)
  integer e_v(2,e_num)
  integer, dimension ( 2, 10 ) :: e_v_save = reshape ( (/ &
  11,  6, &
   6,  4, &
   4,  1, &
   1,  2, &
   2,  5, &
   5,  9, &
   9, 13, &
  13, 15, &
  15, 14, &
  14, 11 /), (/ 2, 10 /) )
  integer t_l(t_num)
  integer, dimension (18) :: t_l_save = (/ &
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
  integer t_v(3,t_num)
  integer, dimension ( 3, 18 ) :: t_v_save = reshape ( (/ &
     1,  3,  4, &
     7,  2,  5, &
     9,  7,  5, &
     8,  6,  4, &
    12,  8,  7, &
    12, 11,  8, &
     3,  1,  2, &
     7,  3,  2, &
     7,  8,  3, &
     4,  3,  8, &
     6,  8, 11, &
    12,  7, 10, &
    11, 12, 14, &
    10,  9, 13, &
    12, 10, 13, &
     7,  9, 10, &
    12, 13, 15, &
    14, 12, 15 /), (/ 3, 18 /) )
  integer v_l(v_num)
  integer, dimension (15) :: v_l_save = (/ &
    1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1 /)
  real ( kind = rk ) v_xy(2,v_num)
  real ( kind = rk ), dimension ( 2, 15 ) :: v_xy_save = reshape ( (/ &
    -0.309016994375D+00,  0.951056516295D+00, &
    -0.809016994375D+00,  0.587785252292D+00, &
    -0.321175165867D+00,  0.475528256720D+00, &
     0.309016994375D+00,  0.951056516295D+00, &
    -1.000000000000D+00,  0.000000000000D+00, &
     0.809016994375D+00,  0.587785252292D+00, &
    -0.333333334358D+00,  0.000000000000D+00, &
     0.237841829972D+00,  0.293892623813D+00, &
    -0.809016994375D+00, -0.587785252292D+00, &
    -0.321175165867D+00, -0.475528259963D+00, &
     1.000000000000D+00,  0.000000000000D+00, &
     0.206011327827D+00, -0.391856835534D+00, &
    -0.309016994375D+00, -0.951056516295D+00, &
     0.809016994375D+00, -0.587785252292D+00, &
     0.309016994375D+00, -0.951056516295D+00 /), (/ 2, 15 /) )

  call i4vec_copy (    v_num, v_l_save,  v_l )
  call r8mat_copy ( 2, v_num, v_xy_save, v_xy )
  call i4vec_copy (    e_num, e_l_save,  e_l )
  call i4mat_copy ( 2, e_num, e_v_save,  e_v )
  call i4vec_copy (    t_num, t_l_save,  t_l )
  call i4mat_copy ( 3, t_num, t_v_save,  t_v )

  return
end
subroutine ffmsh_2d_data_print ( title, v_num, e_num, t_num, v_xy, v_l, e_v, &
  e_l, t_v, t_l )

!*****************************************************************************80
!
!! FFMSH_2D_DATA_PRINT prints FFMSH data.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 December 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) TITLE, a title.
!
!    Input, integer V_NUM, the number of vertices.
!
!    Input, integer E_NUM, the number of boundary edges.
!
!    Input, integer T_NUM, the number of triangles.
!
!    Input, real ( kind = rk ) V_XY(2,V_NUM), vertex coordinates.
!
!    Input, integer V_L(V_NUM), vertex labels.
!
!    Input, integer E_V(2,E_NUM), edge vertices.
!
!    Input, integer E_L(E_NUM), vertex labels.
!
!    Input, integer T_V(3,T_NUM), triangle vertices.
!
!    Input, integer T_L(T_NUM), triangle labels.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer e_num
  integer t_num
  integer v_num

  integer e_l(e_num)
  integer e_v(2,e_num)
  integer t_l(t_num)
  integer t_v(3,t_num)
  character ( len = * ) title
  integer v_l(v_num)
  real ( kind = rk ) v_xy(2,v_num)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) trim ( title )
  
  call i4vec_print (              v_num, v_l,  '  Vertex labels:' )
  call r8mat_transpose_print ( 2, v_num, v_xy, '  Vertex coordinates:' )
  call i4vec_print (              e_num, e_l,  '  Edge labels:' )
  call i4mat_transpose_print ( 2, e_num, e_v,  '  Edge vertices:' )
  call i4vec_print (              t_num, t_l,  '  Triangle labels:' )
  call i4mat_transpose_print ( 3, t_num, t_v,  '  Triangle vertices:' )

  return
end
subroutine ffmsh_2d_data_read ( ffmsh_filename, v_num, e_num, t_num, v_xy, &
  v_l, e_v, e_l, t_v, t_l )

!*****************************************************************************80
!
!! FFMSH_2D_DATA_READ reads data from an FFMSH file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 December 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FFMSH_FILENAME, the FFMSH filename.
!
!    Input, integer V_NUM, the number of vertices.
!
!    Input, integer E_NUM, the number of boundary edges.
!
!    Input, integer T_NUM, the number of triangles.
!
!    Output, real ( kind = rk ) V_XY(2,V_NUM), vertex coordinates.
!
!    Output, integer V_L(V_NUM), vertex labels.
!
!    Output, integer E_V(2,E_NUM), edge vertices.
!
!    Output, integer E_L(E_NUM), vertex labels.
!
!    Output, integer T_V(3,T_NUM), triangle vertices.
!
!    Output, integer T_L(T_NUM), triangle labels.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer e_num
  integer t_num
  integer v_num

  integer e_l(e_num)
  integer e_num2
  integer e_v(2,e_num)
  character ( len = * ) ffmsh_filename
  integer ffmsh_stat
  integer ffmsh_unit
  integer i1
  integer i2
  integer i3
  integer i4
  integer ierror
  integer j
  integer length
  real ( kind = rk ) r1
  real ( kind = rk ) r2
  integer t_l(t_num)
  integer t_num2
  integer t_v(3,t_num)
  character ( len = 255 ) text
  integer v_num2
  integer v_l(v_num)
  real ( kind = rk ) v_xy(2,v_num)

  call get_unit ( ffmsh_unit )

  open ( unit = ffmsh_unit, file = ffmsh_filename, status = 'old', &
    iostat = ffmsh_stat )

  if ( ffmsh_stat /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'FFMSH_2D_DATA_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open input file "' // &
      trim ( ffmsh_filename ) // '"'
    stop 1
  end if
!
!  Read the sizes (again).
!
  read ( ffmsh_unit, '(a)', iostat = ffmsh_stat ) text

  call s_to_i4 ( text, v_num2, ierror, length )
  text = text(length+1:)
  call s_to_i4 ( text, t_num2, ierror, length )
  text = text(length+1:)
  call s_to_i4 ( text, e_num2, ierror, length )
!
!  Read Vertex X, Y, Label
!
  do j = 1, v_num
    read ( ffmsh_unit, '(a)', iostat = ffmsh_stat ) text
    call s_to_r8 ( text, r1, ierror, length )
    text = text(length+1:)
    call s_to_r8 ( text, r2, ierror, length )
    text = text(length+1:)
    call s_to_i4 ( text, i1, ierror, length )
    v_xy(1,j) = r1
    v_xy(2,j) = r2
    v_l(j) = i1
  end do
!
!  Read Triangle V1, V2, V3, Label
!
  do j = 1, t_num
    read ( ffmsh_unit, '(a)', iostat = ffmsh_stat ) text
    call s_to_i4 ( text, i1, ierror, length )
    text = text(length+1:)
    call s_to_i4 ( text, i2, ierror, length )
    text = text(length+1:)
    call s_to_i4 ( text, i3, ierror, length )
    text = text(length+1:)
    call s_to_i4 ( text, i4, ierror, length )
    t_v(1,j) = i1
    t_v(2,j) = i2
    t_v(3,j) = i3
    t_l(j) = i4
  end do
!
!  Read Edge V1, V2, Label
!
  do j = 1, e_num
    read ( ffmsh_unit, '(a)', iostat = ffmsh_stat ) text
    call s_to_i4 ( text, i1, ierror, length )
    text = text(length+1:)
    call s_to_i4 ( text, i2, ierror, length )
    text = text(length+1:)
    call s_to_i4 ( text, i3, ierror, length )
    e_v(1,j) = i1
    e_v(2,j) = i2
    e_l(j) = i3
  end do

  close ( unit = ffmsh_unit )

  return
end
subroutine ffmsh_2d_size_example ( v_num, e_num, t_num )

!*****************************************************************************80
!
!! FFMSH_2D_SIZE_EXAMPLE returns sizes for the 2D example.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 December 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer V_NUM, the number of vertices.
!
!    Output, integer E_NUM, the number of boundary edges.
!
!    Output, integer T_NUM, the number of triangles.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer e_num
  integer t_num
  integer v_num

  e_num = 10
  t_num = 18
  v_num = 15

  return
end
subroutine ffmsh_2d_size_print ( title, v_num, e_num, t_num )

!*****************************************************************************80
!
!! FFMSH_2D_SIZE_PRINT prints the sizes of an FFMSH.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 December 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) TITLE, a title.
!
!    Input, integer V_NUM, the number of vertices.
!
!    Input, integer E_NUM, the number of boundary edges.
!
!    Input, integer T_NUM, the number of triangles.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer e_num
  integer t_num
  character ( len = * ) title
  integer v_num

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ''
  write ( *, '(a,i6)' ) '  Number of vertices = ', v_num
  write ( *, '(a,i6)' ) '  Number of boundary edges = ', e_num
  write ( *, '(a,i6)' ) '  Number of triangles = ', t_num

  return
end
subroutine ffmsh_2d_size_read ( ffmsh_filename, v_num, e_num, t_num )

!*****************************************************************************80
!
!! FFMSH_2D_SIZE_READ reads sizes from a FFMSH file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 December 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FFMSH_FILENAME, the FFMSH filename.
!
!    Output, integer V_NUM, the number of vertices.
!
!    Output, integer E_NUM, the number of boundary edges.
!
!    Output, integer T_NUM, the number of triangles.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer e_num
  character ( len = * ) ffmsh_filename
  integer ffmsh_stat
  integer ffmsh_unit
  integer ierror
  integer length
  integer t_num
  character ( len = 255 ) text
  integer v_num

  call get_unit ( ffmsh_unit )

  open ( unit = ffmsh_unit, file = ffmsh_filename, status = 'old', &
    iostat = ffmsh_stat )

  if ( ffmsh_stat /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'FFMSH_SIZE_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file "' // &
      trim ( ffmsh_filename ) // '"'
    stop 1
  end if

  read ( ffmsh_unit, '(a)', iostat = ffmsh_stat ) text

  call s_to_i4 ( text, v_num, ierror, length )
  text = text(length+1:)
  call s_to_i4 ( text, t_num, ierror, length )
  text = text(length+1:)
  call s_to_i4 ( text, e_num, ierror, length )

  close ( unit = ffmsh_unit )

  return
end
subroutine ffmsh_2d_write ( ffmsh_filename, v_num, e_num, t_num, v_xy, v_l, &
  e_v, e_l, t_v, t_l )

!*****************************************************************************80
!
!! FFMSH_2D_WRITE writes FFMSH data to a file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 December 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FFMSH_FILENAME, the name of the file.
!
!    Input, integer V_NUM, the number of vertices.
!
!    Input, integer E_NUM, the number of boundary edges.
!
!    Input, integer T_NUM, the number of triangles.
!
!    Input, real ( kind = rk ) V_XY(2,V_NUM), vertex coordinates.
!
!    Input, integer V_L(V_NUM), vertex labels.
!
!    Input, integer E_V(2,E_NUM), edge vertices.
!
!    Input, integer E_L(E_NUM), vertex labels.
!
!    Input, integer T_V(3,T_NUM), triangle vertices.
!
!    Input, integer T_L(T_NUM), triangle labels.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer e_num
  integer t_num
  integer v_num

  integer e_l(e_num)
  integer e_v(2,e_num)
  character * ( * ) ffmsh_filename
  integer ffmsh_unit
  integer j
  integer t_l(t_num)
  integer t_v(3,t_num)
  integer v_l(v_num)
  real ( kind = rk ) v_xy(2,v_num)
!
!  Open the file.
!
  call get_unit ( ffmsh_unit )

  open ( unit = ffmsh_unit, file = ffmsh_filename, status = 'replace' )
!
!  Write the data.
!
  write ( ffmsh_unit, '(2x,i6,2x,i6,2x,i6)' ) v_num, t_num, e_num

  do j = 1, v_num
    write ( ffmsh_unit, '(2x,g14.6,2x,g14.6,2x,i6)' ) v_xy(1:2,j), v_l(j)
  end do

  do j = 1, t_num
    write ( ffmsh_unit, '(2x,i6,2x,i6,2x,i6,2x,i6)' ) t_v(1:3,j), t_l(j)
  end do

  do j = 1, e_num
    write ( ffmsh_unit, '(2x,i6,2x,i6,2x,i6)' ) e_v(1:2,j), e_l(j)
  end do

  close ( unit = ffmsh_unit )

  return
end
subroutine i4mat_copy ( m, n, a1, a2 )

!*****************************************************************************80
!
!! I4MAT_COPY copies an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer A1(M,N), the matrix to be copied.
!
!    Output, integer A2(M,N), the copied matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer a1(m,n)
  integer a2(m,n)

  a2(1:m,1:n) = a1(1:m,1:n)

  return
end
subroutine i4mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! I4MAT_TRANSPOSE_PRINT prints an I4MAT, transposed.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer A(M,N), an M by N matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer a(m,n)
  character ( len = * ) title

  call i4mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine i4mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! I4MAT_TRANSPOSE_PRINT_SOME prints some of the transpose of an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4's.
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
!    Input, integer A(M,N), an M by N matrix to be printed.
!
!    Input, integer ILO, JLO, the first row and column to print.
!
!    Input, integer IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: incx = 10
  integer m
  integer n

  integer a(m,n)
  character ( len = 8 ) ctemp(incx)
  integer i
  integer i2
  integer i2hi
  integer i2lo
  integer ihi
  integer ilo
  integer inc
  integer j
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

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(i8)' ) i
    end do

    write ( *, '(''  Row '',10a8)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc

        i = i2lo - 1 + i2

        write ( ctemp(i2), '(i8)' ) a(i,j)

      end do

      write ( *, '(i5,a,10a8)' ) j, ':', ( ctemp(i), i = 1, inc )

    end do

  end do

  return
end
subroutine i4vec_copy ( n, a1, a2 )

!*****************************************************************************80
!
!! I4VEC_COPY copies an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the length of the vectors.
!
!    Input, integer A1(N), the vector to be copied.
!
!    Output, integer A2(N), a copy of A1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer a1(n)
  integer a2(n)

  a2(1:n) = a1(1:n)

  return
end
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 May 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, integer A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer a(n)
  integer i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,a,2x,i12)' ) i, ':', a(i)
  end do

  return
end
subroutine mesh_base_one ( node_num, element_order, element_num, element_node )

!*****************************************************************************80
!
!! MESH_BASE_ONE ensures that the element definition is one-based.
!
!  Discussion:
!
!    The ELEMENT_NODE array contains nodes indices that form elements.
!    The convention for node indexing might start at 0 or at 1.
!    Since a FORTRAN90 program will naturally assume a 1-based indexing, it is
!    necessary to check a given element definition and, if it is actually
!    0-based, to convert it.
!
!    This function attempts to detect 9-based node indexing and correct it.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 October 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, int NODE_NUM, the number of nodes.
!
!    Input, int ELEMENT_ORDER, the order of the elements.
!
!    Input, int ELEMENT_NUM, the number of elements.
!
!    Input/output, int ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM), the element
!    definitions.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer element_num
  integer element_order

  integer element_node(element_order,element_num)
  integer, parameter :: i4_huge = 2147483647
  integer node_max
  integer node_min
  integer node_num

  node_min = + i4_huge
  node_max = - i4_huge

  node_min = minval ( element_node(1:element_order,1:element_num) )
  node_max = maxval ( element_node(1:element_order,1:element_num) )

  if ( node_min == 0 .and. node_max == node_num - 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' )'MESH_BASE_ONE:'
    write ( *, '(a)' )'  The element indexing appears to be 0-based!'
    write ( *, '(a)' )'  This will be converted to 1-based.'
    element_node(1:element_order,1:element_num) = &
      element_node(1:element_order,1:element_num) + 1
  else if ( node_min == 1 .and. node_max == node_num  ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' )'MESH_BASE_ONE:'
    write ( *, '(a)' )'  The element indexing appears to be 1-based!'
    write ( *, '(a)' )'  No conversion is necessary.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MESH_BASE_ONE - Warning!'
    write ( *, '(a)' ) '  The element indexing is not of a recognized type.'
    write ( *, '(a,i8)' ) '  NODE_MIN = ', node_min
    write ( *, '(a,i8)' ) '  NODE_MAX = ', node_max
    write ( *, '(a,i8)' ) '  NODE_NUM = ', node_num
  end if

  return
end
subroutine r8mat_copy ( m, n, a, b )

!*****************************************************************************80
!
!! R8MAT_COPY copies an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> (I+J*M).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 July 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the order of the matrix.
!
!    Input, real ( kind = rk ) A(M,N), the matrix to be copied.
!
!    Output, real ( kind = rk ) B(M,N), a copy of the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) b(m,n)

  b(1:m,1:n) = a(1:m,1:n)

  return
end
subroutine r8mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 June 2004
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
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  character ( len = * ) title

  call r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
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
  integer i2
  integer i2hi
  integer i2lo
  integer ihi
  integer ilo
  integer inc
  integer j
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

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(i8,6x)' ) i
    end do

    write ( *, '(''  Row   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc
        i = i2lo - 1 + i2
        write ( ctemp(i2), '(g14.6)' ) a(i,j)
      end do

      write ( *, '(i5,a,5a14)' ) j, ':', ( ctemp(i), i = 1, inc )

    end do

  end do

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
!    Output, integer IVAL, the integer value read from the string.
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

  integer, parameter :: rk = kind ( 1.0D+00 )

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
!    '-12.73e-9.23'   -12.73 * 10.0^(-9.23)
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

  character c
  logical ch_eqi
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
  if ( iterm /= 1 .and. length + 1 == nchar ) then
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

  integer, parameter :: rk = kind ( 1.0D+00 )

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
