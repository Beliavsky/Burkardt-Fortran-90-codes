program main

!*****************************************************************************80
!
!! FILE_NAME_SEQUENCE_TEST tests the FILE_NAME_SEQUENCE library.
!
!  Discussion:
!
!    There are situations such as animations or parallel processing in which
!    it is necessary to generate a sequence of file names which include
!    an embedded index that increases.  A simple example might be
!
!      "fred0.txt", "fred1.txt", "fred2.txt"
!
!    A side issue arises when the number of files is large enough that the
!    number of digits in the index will vary.  Thus, if we are going to have
!    15 files, do we want to number them as
!
!      "fred00.txt" through "fred14.txt"
!
!    which means, for one thing, that they will alphabetize properly, or
!    will we be satisfied with
!
!      "fred0.txt" through "fred14.txt" ?
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 July 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = 80 ) filename

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FILE_NAME_SEQUENCE_TEST:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) &
    '  Demonstrate ways of generating a numeric sequence of file names.'

  filename = 'frodo_01345_lives.txt'
  call test04 ( filename, 10 )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FILE_NAME_SEQUENCE_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test04 ( filename, filename_num )

!*****************************************************************************80
!
!!  TEST04 uses FILE_NAME_INC.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 July 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = * ) filename
  integer filename_num
  integer i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04:'
  write ( *, '(a)' ) '  FILENAME(I+1) = FILE_NAME_INC ( FILENAME(I) )'
  write ( *, '(a)' ) '  First FILENAME = "' // trim ( filename ) // '".'
  write ( *, '(a,i4)' ) '  Number of filenames = ', filename_num
  write ( *, '(a)' ) '  Numbers may include leading zeros.'
  write ( *, '(a)' ) ' '

  do i = 1, filename_num
    write ( *, '(2x,i4,2x,a)' )  i, '"' // trim ( filename ) // '".'
    call file_name_inc ( filename )
  end do

  return
end
