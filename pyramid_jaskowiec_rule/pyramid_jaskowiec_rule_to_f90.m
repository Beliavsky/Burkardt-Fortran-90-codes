data = load ( 'cubature_pyr_sym_p19_n418_expand.txt' );
output = fopen ( 'rule19.f90', 'w' );
n = size ( data, 1 );

fprintf ( output, 'subroutine rule19 ( n, x, y, z, w )\n' ); 
fprintf ( output, '\n' );
fprintf ( output, '!*****************************************************************************80\n' );
fprintf ( output, '!\n' );
fprintf ( output, '!! rule19() returns the pyramid quadrature rule of precision 19.\n' );
fprintf ( output, '!\n' );
fprintf ( output, '!  Licensing:\n' );
fprintf ( output, '!\n' );
fprintf ( output, '!    This code is distributed under the GNU LGPL license.\n' );
fprintf ( output, '!\n' );
fprintf ( output, '!  Modified:\n' );
fprintf ( output, '!\n' );
fprintf ( output, '!    15 April 2023\n' );
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
fprintf ( output, '!    integer n: the number of quadrature points.\n' );
fprintf ( output, '!\n' );
fprintf ( output, '!  Output:\n' );
fprintf ( output, '!\n' );
fprintf ( output, '!    real ( kind = rk ) x[n], y[n], z[n]: the coordinates of\n' );
fprintf ( output, '!    quadrature points.\n' );
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
fprintf ( output, '  real ( kind = rk ) x(n)\n' );
fprintf ( output, '  real ( kind = rk ), save, dimension ( n_save ) :: x_save = (/  &\n' );
for i = 1 : n - 1
  fprintf ( output, '    %24.20fD+00, &\n', data(i,1) );
end
fprintf ( output, '    %24.20fD+00 /)\n', data(n,1) );

fprintf ( output, '\n' );
fprintf ( output, '  real ( kind = rk ) y(n)\n' );
fprintf ( output, '  real ( kind = rk ), save, dimension ( n_save ) :: y_save = (/ &\n' );
for i = 1 : n - 1
  fprintf ( output, '    %24.20fD+00, &\n', data(i,2) );
end
fprintf ( output, '    %24.20fD+00 /)\n', data(n,2) );

fprintf ( output, '\n' );
fprintf ( output, '  real ( kind = rk ) z(n)\n' );
fprintf ( output, '  real ( kind = rk ), save, dimension ( n_save ) :: z_save = (/ &\n' );
for i = 1 : n - 1
  fprintf ( output, '    %24.20fD+00, &\n', data(i,3) );
end
fprintf ( output, '    %24.20fD+00 /)\n', data(n,3) );

fprintf ( output, '\n' );
fprintf ( output, '  real ( kind = rk ) w(n)\n' );
fprintf ( output, '  real ( kind = rk ), save, dimension ( n_save ) :: w_save = (/ &\n' );
for i = 1 : n - 1
  fprintf ( output, '    %24.20fD+00, &\n', data(i,4) );
end
fprintf ( output, '    %24.20fD+00 /)\n', data(n,4) );
fprintf ( output, '\n' );

fprintf ( output, '  x(1:n) = x_save(1:n)\n' );
fprintf ( output, '  y(1:n) = y_save(1:n)\n' );
fprintf ( output, '  z(1:n) = z_save(1:n)\n' );
fprintf ( output, '  w(1:n) = w_save(1:n)\n' );

fprintf ( output, '\n' );
fprintf ( output, '  return\n' );
fprintf ( output, 'end\n' );

fclose ( output );
