data = load ( 'cubature_tet_sym_p20_n552_expand_baryc.txt' );
output = fopen ( 'rule20.f90', 'w' );
n = size ( data, 1 );

fprintf ( output, 'subroutine rule20 ( n, a, b, c, d, w )\n' ); 
fprintf ( output, '\n' );
fprintf ( output, '!*****************************************************************************80\n' );
fprintf ( output, '!\n' );
fprintf ( output, '!! rule20() returns the rule of precision 20.\n' );
fprintf ( output, '!\n' );
fprintf ( output, '!  Licensing:\n' );
fprintf ( output, '!\n' );
fprintf ( output, '!    This code is distributed under the GNU LGPL license.\n' );
fprintf ( output, '!\n' );
fprintf ( output, '!  Modified:\n' );
fprintf ( output, '!\n' );
fprintf ( output, '!    13 April 2023\n' );
fprintf ( output, '!\n' );
fprintf ( output, '!  Author:\n' );
fprintf ( output, '!\n' );
fprintf ( output, '!    John Burkardt\n' );
fprintf ( output, '!\n' );
fprintf ( output, '!  Reference:\n' );
fprintf ( output, '!\n' );
fprintf ( output, '!    Jan Jaskowiec, Natarajan Sukumar,\n' );
fprintf ( output, '!    High order cubature rules for tetrahedra and pyramids,\n' );
fprintf ( output, '!    International Journal of Numerical Methods in Engineering,\n' );
fprintf ( output, '!    Volume 121, Number 11, pages 2418-2436, 15 June 2020.\n' );
fprintf ( output, '!\n' );
fprintf ( output, '!  Input:\n' );
fprintf ( output, '!\n' );
fprintf ( output, '!    integer n: the number of quadrature points for this rule.\n' );
fprintf ( output, '!\n' );
fprintf ( output, '!  Output:\n' );
fprintf ( output, '!\n' );
fprintf ( output, '!    real ( kind = rk ) a[n], b[n], c[n], d[n]: the barycentric\n' );
fprintf ( output, '!    coordinates of quadrature points.\n' );
fprintf ( output, '!\n' );
fprintf ( output, '!    real ( kind = rk ) w[n]: the quadrature weights.\n' );
fprintf ( output, '!\n' );
fprintf ( output, '  implicit none\n' );
fprintf ( output, '\n' );
fprintf ( output, '  integer, parameter :: rk = kind ( 1.0D+00)\n' );
fprintf ( output, '\n' );
fprintf ( output, '  integer n\n' );
fprintf ( output, '  integer, parameter :: n_save = 552\n' );
fprintf ( output, '\n' );
fprintf ( output, '  real ( kind = rk ) a(n)\n' );
fprintf ( output, '  real ( kind = rk ), save, dimension ( n_save ) :: a_save = (/  &\n' );
for i = 1 : n - 1
  fprintf ( output, '    %24.20fD+00, &\n', data(i,1) );
end
fprintf ( output, '    %24.20fD+00 /)\n', data(n,1) );

fprintf ( output, '\n' );
fprintf ( output, '  real ( kind = rk ) b(n)\n' );
fprintf ( output, '  real ( kind = rk ), save, dimension ( n_save ) :: b_save = (/ &\n' );
for i = 1 : n - 1
  fprintf ( output, '    %24.20fD+00, &\n', data(i,2) );
end
fprintf ( output, '    %24.20fD+00 /)\n', data(n,2) );

fprintf ( output, '\n' );
fprintf ( output, '  real ( kind = rk ) c(n)\n' );
fprintf ( output, '  real ( kind = rk ), save, dimension ( n_save ) :: c_save = (/ &\n' );
for i = 1 : n - 1
  fprintf ( output, '    %24.20fD+00, &\n', data(i,3) );
end
fprintf ( output, '    %24.20fD+00 /)\n', data(n,3) );

fprintf ( output, '\n' );
fprintf ( output, '  real ( kind = rk ) d(n)\n' );
fprintf ( output, '  real ( kind = rk ), save, dimension ( n_save ) :: d_save = (/ &\n' );
for i = 1 : n - 1
  fprintf ( output, '    %24.20fD+00, &\n', data(i,4) );
end
fprintf ( output, '    %24.20fD+00 /)\n', data(n,4) );

fprintf ( output, '\n' );
fprintf ( output, '  real ( kind = rk ) w(n)\n' );
fprintf ( output, '  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &\n' );
for i = 1 : n - 1
  fprintf ( output, '    %24.20fD+00, &\n', data(i,5) );
end
fprintf ( output, '    %24.20fD+00 /)\n', data(n,5) );
fprintf ( output, '\n' );

fprintf ( output, '  call r8vec_copy ( n, a_save, a )\n' );
fprintf ( output, '  call r8vec_copy ( n, b_save, b )\n' );
fprintf ( output, '  call r8vec_copy ( n, c_save, c )\n' );
fprintf ( output, '  call r8vec_copy ( n, d_save, d )\n' );
fprintf ( output, '  call r8vec_copy ( n, w_save, w )\n' );

fprintf ( output, '\n' );
fprintf ( output, '  return\n' );
fprintf ( output, 'end\n' );

fclose ( output );
