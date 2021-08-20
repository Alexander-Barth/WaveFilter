% function test_test()
% Tests filter_vort and filter_var. 

function filter_test()

[X,Y] = ndgrid(linspace(-90,-80,30),linspace(24,31,30));

domain.h = 100 + X+Y;

[m,n] = size(X);

% Earth Radius
Re = 6371e3;
s =  pi*Re/180;

dx = X(2,1) - X(1,1);
dy = Y(1,2) - Y(1,1);

dx = dx * cos(28*pi/180) * s;
dy = dy * s;

domain.pm = ones(m,n)/dx;
domain.pn = ones(m,n)/dy;
domain.mask = ones(m,n);
domain.f = 1e-4;
domain.g = 9.81;

zetai = exp(- (X+86).^2/0.2 - (Y-28).^2/0.2);
ui = zeros(m-1,n);
vi = zeros(m,n-1);

sigma = 50;

tol = 1e-8;

fprintf('Testing filter_vort: ');
[zetaf1,uf1,vf1] = filter_vort(domain,zetai,ui,vi);

if (abs(zetaf1(20,20) - 0.00685792232605)) < tol
  fprintf('OK\n');
else
  fprintf('FAILED\n');
end

fprintf('Testing filter_var:  ');
[zetaf2,uf2,vf2] = filter_var(domain,zetai,ui,vi,sigma);

if (abs(zetaf2(20,20) - 0.00685344988792)) < tol
  fprintf('OK\n');
else
  fprintf('FAILED\n');
end
