
function [zetaf,Q] = diffequation_var(domain,q,bc,zeta_bc,sigma)

mask = domain.mask;
pm = domain.pm;
pn = domain.pn;
g = domain.g;
f = domain.f;
h = domain.h;
R = sqrt(g*h)/f;

[m,n] = size(q);

zetaf = zeros(m,n);
Q = zeros(m,n);


Op = diffequation_op(m,n,R,pm,pn,mask,bc);


pm_v = (pm(:,1:end-1)+pm(:,2:end))/2;
pm_u = (pm(1:end-1,:)+pm(2:end,:))/2;
pn_v = (pn(:,1:end-1)+pn(:,2:end))/2;
pn_u = (pn(1:end-1,:)+pn(2:end,:))/2;


[hx_u,hx_v,hy_u,hy_v] = gradientc(h,m,n,1,1);

diffx = sparse_diffx(m,n);
diffy = sparse_diffy(m,n);
u2v = sparse_u2v(m,n);
v2u = sparse_v2u(m,n);
diffx_u = sparse_diffx(m-1,n);
diffy_v = sparse_diffy(m,n-1);
trimx = sparse_trimx(m,n-2);
trimy = sparse_trimy(m-2,n);
ext = sparse_ext(m-2,n-2);

mn_u = sparse_diag(pm_u .* pn_u);
mn_v = sparse_diag(pm_v .* pn_v);
mn = sparse_diag(pm .* pn);


P = sigma * (trimy * diffx_u * mn_u * (sparse_diag(hy_u.^2) * diffx  - sparse_diag(hx_u.* hy_u)  * v2u * diffy) + ...
             trimx * diffy_v * mn_v * (sparse_diag(hx_v.^2) * diffy  - sparse_diag(hx_v.* hy_v)  * u2v * diffx));

P = mn * ext * P;

% loop version
% $$$ for j=1:n
% $$$   for i=1:m
% $$$     if ~mask(i,j)
% $$$       l = sub2ind([m n],i,j);    
% $$$       P(l,:) = 0;
% $$$     end      
% $$$   end 
% $$$ end

% or vectorized
P = sparse_diag(mask) * P;


A = [Op -speye(m*n,m*n);  -P  Op ];

% boundary conditions

q([1 m],:) = zeta_bc([1 m],:);
q(:,[1 n]) = zeta_bc(:,[1 n]);


b = zeros(2*m*n,1);
b(1:m*n) = q(:);

x = A \ b;

zetaf = reshape(x(1:m*n),m,n);
Q = reshape(x(m*n+1:end),m,n);

