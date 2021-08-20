function diffy = sparse_diffy(m,n)

% $$$ diffy = sparse( (n-1)*m , n*m );
% $$$ 
% $$$ for j=1:n-1
% $$$   for i=1:m
% $$$     l = sub2ind([m (n-1)],i,j);    
% $$$     
% $$$     diffy(l,sub2ind([m n],i,j+1)) = 1;
% $$$     diffy(l,sub2ind([m n],i,j)) = -1;
% $$$   end 
% $$$ end
% $$$ 


[I,J] = ndgrid(1:m,1:n-1);

L1 = [1:m*(n-1)]';

L2 = sub2ind([m n],I(:),J(:));
one = ones(size(L1));

diffy = sparse( ...
       [L1;     L1;   ],  ...
       [L2;     L2+m; ],  ...
       [-one;  one  ], m*(n-1) , m*n );

%max(abs(diffy(:) - diffy2(:)))