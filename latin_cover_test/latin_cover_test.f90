program main

!*****************************************************************************80
!
!! latin_cover_test() tests latin_cover().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 June 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'latin_cover_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test latin_cover().'

  call test01 ( )
  call test02 ( )
  call test03 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LATIN_COVER_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests LATIN_COVER.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 June 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, allocatable :: a(:,:)
  integer, parameter :: base = 1
  integer n
  integer, allocatable :: p(:)
  integer seed
  integer test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  LATIN_COVER:'

  do n = 3, 9, 2

    allocate ( p(1:n) )
    allocate ( a(1:n,1:n) )

    seed = 123456789

    do test = 1, 3

      call perm_uniform ( n, base, seed, p )
 
      call perm_print ( n, p, '  Permutation' )

      call latin_cover ( n, p, a )

      call i4mat_print ( n, n, a, '  Latin cover' )

    end do

    deallocate ( a )
    deallocate ( p )

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests LATIN_COVER_2D.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 June 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, allocatable :: a(:,:)
  integer, parameter :: base = 1
  integer n
  integer, allocatable :: p(:,:)
  integer seed
  integer test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  LATIN_COVER_2D:'

  do n = 3, 9, 2

    allocate ( p(2,1:n) )
    allocate ( a(1:n,1:n) )

    seed = 123456789

    do test = 1, 3

      call perm_uniform ( n, base, seed, p(1,1:n) )
 
      call perm_print ( n, p(1,1:n), '  Permutation 1' )

      call perm_uniform ( n, base, seed, p(2,1:n) ) 
 
      call perm_print ( n, p(2,1:n), '  Permutation 2' )
      call latin_cover_2d ( n, p, a )

      call i4mat_print ( n, n, a, '  Latin cover' )

    end do

    deallocate ( a )
    deallocate ( p )

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests LATIN_COVER_3D.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 June 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, allocatable :: a(:,:,:)
  integer, parameter :: base = 1
  integer n
  integer, allocatable :: p(:,:)
  integer seed
  integer test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  LATIN_COVER_3D:'

  do n = 3, 9, 2

    allocate ( p(3,1:n) )
    allocate ( a(1:n,1:n,1:n) )

    seed = 123456789

    do test = 1, 3

      call perm_uniform ( n, base, seed, p(1,1:n) )
 
      call perm_print ( n, p(1,1:n), '  Permutation 1' )

      call perm_uniform ( n, base, seed, p(2,1:n) ) 
 
      call perm_print ( n, p(2,1:n), '  Permutation 2' )

      call perm_uniform ( n, base, seed, p(3,1:n) ) 
 
      call perm_print ( n, p(3,1:n), '  Permutation 3' )

      call latin_cover_3d ( n, p, a )

      call i4block_print ( n, n, n, a, '  Latin cover' )

    end do

    deallocate ( a )
    deallocate ( p )

  end do

  return
end
