function [phi,res] = diffequation_vort(domain,q,bc,zeta_bc)

%
% iterative solves the equation
%
% \partial_x \partial_x phi + \partial_y \partial_y phi - phi/R^2= q
%
% spatial vaying R

% R = external Rossby Radius of deformation
% R = c/f = sqrt(g h)/f
% 

pm = domain.pm;
pn = domain.pn;
mask = domain.mask;
h = domain.h;
R = sqrt(domain.g*h)/domain.f;


[m,n] = size(q);
phi = zeros(m,n);

D = diffequation_op(m,n,R,pm,pn,mask,bc);

q(~mask) = zeta_bc(~mask);
q([1 m],:) = zeta_bc([1 m],:);
q(:,[1 n]) = zeta_bc(:,[1 n]);

phi(:) = D \ q(:);

res = D * phi(:) - q(:);
