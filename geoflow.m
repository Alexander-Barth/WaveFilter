function [uf,vf] = geoflow(domain,zetaf)

g = domain.g;
f = domain.f;
[m,n] = size(zetaf);
Uf = zeros(m-1,n);
Vf = zeros(m,n-1);

pm_v = (domain.pm(:,1:end-1)+domain.pm(:,2:end))/2;
pn_u = (domain.pn(1:end-1,:)+domain.pn(2:end,:))/2;


uf = zeros(m-1,n);
vf = zeros(m,n-1);

%zetaf = x;

i=1:m-1;
j=2:n-1;
uf(i,j) = -g * pn_u(i,j) /(4*f) .* (zetaf(i,j+1) + zetaf(i+1,j+1) - zetaf(i,j-1) - zetaf(i+1,j-1));

i=2:m-1;
j=1:n-1;   
vf(i,j) =  g * pm_v(i,j) /(4*f) .* (zetaf(i+1,j) + zetaf(i+1,j+1) - zetaf(i-1,j) - zetaf(i-1,j+1));

