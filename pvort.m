function qf = pvort(domain,zetaf,uf,vf);

pm = domain.pm;
pn = domain.pn;
h = domain.h;
g = domain.g;
f = domain.f;

pm_v = (pm(:,1:end-1)+pm(:,2:end))/2;
pm_u = (pm(1:end-1,:)+pm(2:end,:))/2;
pn_v = (pn(:,1:end-1)+pn(:,2:end))/2;
pn_u = (pn(1:end-1,:)+pn(2:end,:))/2;

[m,n] = size(zetaf);


uf = uf ./ pm_u;
vf = vf ./ pn_v;

R = sqrt(g * h)/f;

qf = -zetaf ./R.^2;

i=2:m-1;
j=2:n-1;    

qf(i,j) = qf(i,j)  ...
          + 0.25 * f/g* pm(i,j) .* pn(i,j) .* (vf(i+1,j)+vf(i+1,j-1)  - vf(i-1,j)-vf(i-1,j-1)) ...
          - 0.25 * f/g* pm(i,j) .* pn(i,j) .* (uf(i,j+1)+uf(i-1,j+1)  - uf(i,j-1)-uf(i-1,j-1));
