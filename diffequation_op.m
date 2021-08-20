function D = diffequation_op(m,n,R,pm,pn,mask,bc)

%
% creates the differential operator
%
% \partial_x \partial_x phi + \partial_y \partial_y phi - phi/R^2
%
% spatial varying R

% R = external Rossby Radius of deformation
% R = c/f = sqrt(g h)/f
% 

% Note: ROMS's pm and pn are the inverse of the scale factor h1, h2 
% http://mathworld.wolfram.com/ScaleFactor.html
% Divergence:
% http://mathworld.wolfram.com/Divergence.html

% h1 = 1/pm
% h2 = 1/pn


% compatibility with previous versions
if (isscalar(pm) & isscalar(pn))
  pm = ones(m,n) / pm;
  pn = ones(m,n) / pn;
end

% metric

pm_v = (pm(:,1:end-1)+pm(:,2:end))/2;
pm_u = (pm(1:end-1,:)+pm(2:end,:))/2;
pn_v = (pn(:,1:end-1)+pn(:,2:end))/2;
pn_u = (pn(1:end-1,:)+pn(2:end,:))/2;

% scaled by the grids metric

f_u = pm_u./pn_u;
f_v = pn_v./pm_v;

iS = pm .* pn;
iR2 = R.^(-2);

if 0
  % alternate

  % square of the local external Rossby radius of deformation

  RE = R.^2;

  RE_u = 0.5 * (RE(2:m,:)+ RE(1:m-1,:));
  RE_v = 0.5 * (RE(:,2:n)+ RE(:,1:n-1));
  
  iS = iS.* iR2;
  f_u = f_u .* RE_u;
  f_v = f_v .* RE_v;
end

% loop version

% $$$ Lap = sparse( n*m , n*m );
% $$$ for j=2:n-1
% $$$   for i=2:m-1
% $$$     l = sub2ind([m n],i,j);    
% $$$     if mask(i,j)
% $$$       
% $$$       
% $$$       Lap(l,sub2ind([m n],i+1,j)) = f_u(i,j)   * iS(i,j);
% $$$       Lap(l,sub2ind([m n],i-1,j)) = f_u(i-1,j) * iS(i,j);
% $$$       Lap(l,sub2ind([m n],i,j+1)) = f_v(i,j)   * iS(i,j);
% $$$       Lap(l,sub2ind([m n],i,j-1)) = f_v(i,j-1) * iS(i,j);
% $$$       Lap(l,sub2ind([m n],i,j)) = -iR2(i,j) - (f_u(i,j) + f_u(i-1,j) + f_v(i,j) + f_v(i,j-1)) * iS(i,j);
% $$$     else
% $$$       Lap(l,sub2ind([m n],i,j)) = 1;
% $$$     end      
% $$$   end 
% $$$ end

% vectorized version

mask2 = mask(2:m-1,2:n-1);

% all sea points
[I,J] = ndgrid(2:m-1,2:n-1);
I=I(mask2==1);
J=J(mask2==1);

L   = sub2ind([m n  ],I,J);
L_u = sub2ind([m-1 n],I,J);
L_v = sub2ind([m n-1],I,J);

diag = -iR2(L) - (f_u(L_u) + f_u(L_u-1) + f_v(L_v) + f_v(L_v-m)) .* iS(L);


Lap = sparse( ...
    [L;    L;               L;             L;               L            ],  ...
    [L;    L-m;             L+m;           L-1;             L+1          ],  ...
    [diag; f_v(L_v-m).*iS(L); f_v(L_v).*iS(L); f_u(L_u-1).*iS(L); f_u(L_u).*iS(L)], n*m , n*m );

% all land points

[I,J] = ndgrid(2:m-1,2:n-1);
I=I(mask2==0);
J=J(mask2==0);

if ~isempty(I)
  L = sub2ind([m n],I,J);
  Lap = Lap + sparse(L,L,ones(size(L)), n*m , n*m );
end

BC = sparse( n*m , n*m );

if strcmp(bc,'periodic')

  for j=2:n-1
    % compute the difference between (1,j) and (m-1,j)
    % which must be zero
    
    BC(sub2ind([m n],1,j),sub2ind([m n],1,j)) = 1;
    BC(sub2ind([m n],1,j),sub2ind([m n],m-1,j)) = -1;

    % compute the difference between (m,j) and (2,j)
    % which must be zero
    
    BC(sub2ind([m n],m,j),sub2ind([m n],m,j)) = 1;
    BC(sub2ind([m n],m,j),sub2ind([m n],2,j)) = -1;
  end

  for i=1:m
    % compute the difference between (i,1) and (i,n-1)
    % which must be zero
    
    BC(sub2ind([m n],i,1),sub2ind([m n],i,1)) = 1;
    BC(sub2ind([m n],i,1),sub2ind([m n],i,n-1)) = -1;

    % compute the difference between (i,n) and (i,2)
    % which must be zero
    
    BC(sub2ind([m n],i,n),sub2ind([m n],i,n)) = 1;
    BC(sub2ind([m n],i,n),sub2ind([m n],i,2)) = -1;
  end

elseif strcmp(bc,'dirichlet')
% $$$   for j=2:n-1
% $$$     % values (1,j)
% $$$ 
% $$$     BC(sub2ind([m n],1,j),sub2ind([m n],1,j)) = 1;
% $$$ 
% $$$     % values (m,j)
% $$$     
% $$$     BC(sub2ind([m n],m,j),sub2ind([m n],m,j)) = 1;
% $$$   end

  [I,J] = ndgrid([1 m],2:n-1);
  L = sub2ind([m n],I(:),J(:));  
  BC = BC + sparse(L,L,ones(size(L)),m*n,m*n);
  
% $$$   for i=1:m
% $$$     % values (i,1)
% $$$     
% $$$     BC(sub2ind([m n],i,1),sub2ind([m n],i,1)) = 1;
% $$$ 
% $$$     % values (i,n)
% $$$     
% $$$     BC(sub2ind([m n],i,n),sub2ind([m n],i,n)) = 1;
% $$$   end

  [I,J] = ndgrid(1:m,[1 n]);
  L = sub2ind([m n],I(:),J(:));  
  BC = BC + sparse(L,L,ones(size(L)),m*n,m*n);

else
  error('unknown boundary conditions');
end

D = Lap + BC;

