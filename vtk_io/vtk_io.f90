subroutine digit_to_ch ( digit, ch )

!*****************************************************************************80
!
!! digit_to_ch() returns the character representation of a decimal digit.
!
!  Discussion:
!
!    Instead of CHAR, we now use the ACHAR function, which
!    guarantees the ASCII collating sequence.
!
!  Example:
!
!    DIGIT   CH 
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
!    Output, character CH, the corresponding character.
!
  implicit none

  character ch
  integer digit

  if ( 0 <= digit .and. digit <= 9 ) then

    ch = achar ( digit + 48 )

  else

    ch = '*'

  end if
 
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
!    26 October 2008
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
subroutine i4_to_s_left ( i4, s )

!*****************************************************************************80
!
!! i4_to_s_left() converts an I4 to a left-justified string.
!
!  Discussion:
!
!    An I4 is an integer.
!
!  Example:
!
!    Assume that S is 6 characters long:
!
!        I4  S
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
!    28 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I4, an integer to be converted.
!
!    Output, character ( len = * ) S, the representation of the integer.
!    The integer will be left-justified.  If there is not enough space,
!    the string will be filled with stars.
!
  implicit none

  character c
  integer i
  integer i4
  integer idig
  integer ihi
  integer ilo
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
  ival = i4
!
!  Handle the negative sign.
!
  if ( ival < 0 ) then

    if ( ihi <= 1 ) then
      s(1:1) = '*'
      return
    end if

    ival = -ival
    s(1:1) = '-'
    ilo = 2

  end if
!
!  The absolute value of the integer goes into S(ILO:IHI).
!
  ipos = ihi
!
!  Find the last digit of IVAL, strip it off, and stick it into the string.
!
  do

    idig = mod ( ival, 10 )
    ival = ival / 10

    if ( ipos < ilo ) then
      do i = 1, ihi
        s(i:i) = '*'
      end do
      return
    end if

    call digit_to_ch ( idig, c )

    s(ipos:ipos) = c
    ipos = ipos - 1

    if ( ival == 0 ) then
      exit
    end if

  end do
!
!  Shift the string to the left.
!
  s(ilo:ilo+ihi-ipos-1) = s(ipos+1:ihi)
  s(ilo+ihi-ipos:ihi) = ' '
 
  return
end
subroutine vtk_puvw_write ( output_unit, title, node_num, element_num, &
  element_order, xyz, element_node, p, uvw )

!*****************************************************************************80
!
!! vtk_puvw_write() writes pressure and velocity data to a VTK file.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 March 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer OUTPUT_UNIT, the unit number of the file.
!
!    Input, character ( len = * ) TITLE, a title for the data.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer ELEMENT_ORDER, the order of the elements.
!
!    Input, real ( kind = rk8 ) XYZ(3,NODE_NUM), the node coordinates.
!
!    Input, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM), the
!    nodes that make up each element.
!
!    Input, real ( kind = rk8 ) P(NODE_NUM), the pressure at each node.
!
!    Input, real ( kind = rk8 ) UVW(3,NODE_NUM), the velocity at each node.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer element_num
  integer element_order
  integer node_num

  character ( len = 20 ) cell_size_string
  integer element
  integer element_node(element_order,element_num)
  character ( len = 20 ) element_num_string
  integer node
  character ( len = 20 ) node_num_string
  integer output_unit
  real ( kind = rk8 ) p(node_num)
  character ( len = * ) title
  real ( kind = rk8 ) uvw(3,node_num)
  real ( kind = rk8 ) xyz(3,node_num)

  call i4_to_s_left ( node_num, node_num_string )
  call i4_to_s_left ( element_num, element_num_string )
  call i4_to_s_left ( element_num * ( element_order + 1 ), cell_size_string )

  write ( output_unit, '(a)' ) '# vtk DataFile Version 2.0'
  write ( output_unit, '(a)' ) title
  write ( output_unit, '(a)' ) 'ASCII'
  write ( output_unit, '(a)' ) 'DATASET UNSTRUCTURED_GRID'
  write ( output_unit, '(a)' ) 'POINTS ' // trim ( node_num_string ) // ' double'

  do node = 1, node_num
    write ( output_unit, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) xyz(1:3,node)
  end do
!
!  Note that the element node indices must be converted from 1-based to 0-based
!  before being written out.
!
  write ( output_unit, '(a)' ) ' '
  write ( output_unit, '(a)' ) 'CELLS ' // trim ( element_num_string ) &
    // ' ' // trim ( cell_size_string )
  do element = 1, element_num
    write ( output_unit, '(2x,i4,10(2x,i4))' ) &
      element_order, element_node(1:element_order,element) - 1
  end do

  write ( output_unit, '(a)' ) ' '
  write ( output_unit, '(a)' ) 'CELL_TYPES ' // trim ( element_num_string )

  if ( element_order == 4 ) then
    do element = 1, element_num
      write ( output_unit, '(a)' ) '10'
    end do
  else if ( element_order == 10 ) then
    do element = 1, element_num
      write ( output_unit, '(a)' ) '24'
    end do
  end if

  write ( output_unit, '(a)' ) ' '
  write ( output_unit, '(a)' ) 'POINT_DATA ' // trim ( node_num_string )
  write ( output_unit, '(a)' ) 'SCALARS pressure double'
  write ( output_unit, '(a)' ) 'LOOKUP_TABLE default'
  do node = 1, node_num
    write ( output_unit, '(2x,g14.6)' ) p(node)
  end do
  write ( output_unit, '(a)' ) 'VECTORS velocity double'
  do node = 1, node_num
    write ( output_unit, '(3(2x,g14.6))' ) uvw(1:3,node)
  end do

  return
end
