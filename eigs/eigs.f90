subroutine eigs ( n, A, lambda )

!*****************************************************************************80
!
!! eigs() computes the eigenvalues of a real matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 June 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer n: the order of the matrix.
!
!    real A(n,n): the matrix.
!
!  Output:
!
!    complex lambda(n): the eigenvalues.
!
  implicit none

  integer, parameter :: ck8 = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  real ( kind = rk8 ) A(n,n)
  real ( kind = rk8 ), allocatable :: A2(:,:)
  integer info
  character jobvl
  character jobvr
  complex ( kind = ck8 ) lambda(n)
  integer lda
  integer ldvl
  integer ldvr
  integer lwork
  real ( kind = rk8 ) vl(0)
  real ( kind = rk8 ) vr(0)
  real ( kind = rk8 ) wi(n)
  real ( kind = rk8 ), allocatable :: work(:)
  real ( kind = rk8 ) wr(n)

  lwork = 3 * n
  allocate ( work(1:lwork) )
!
!  Make a copy of A, since dgeev() will modify the matrix during processing.
!
  allocate ( A2(1:n,1:n) )
  A2(1:n,1:n) = A(1:n,1:n)
!
!  Do not compute left or right eigenvectors!
!  ldvl and ldvr should be 0, but dgeev will reject these as illegal.
!  We set them to 1, which doesn't matter because they aren't used!
!
  jobvl = 'N'
  jobvr = 'N'
  lda = n
  ldvl = 1
  ldvr = 1

  call dgeev ( jobvl, jobvr, n, A2, lda, wr, wi, vl, ldvl, vr, ldvr, work, &
    lwork, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'eigs(): Fatal error!'
    write ( *, '(a,i8)' ) '  dgeev() error flag INFO = ', info
  end if
!
!  Combine wr and wi into the complex values lambda.
!
  lambda = cmplx ( wr, wi, kind = ck8 )
!
!  Release memory.
!
  deallocate ( A2 )
  deallocate ( work )

  return
end

