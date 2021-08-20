function trimx = sparse_trimx(m,n)

% $$$ 
% $$$ trimx = sparse( (m-2)*n , n*m );
% $$$ 
% $$$ for j=1:n
% $$$   for i=1:m-2
% $$$     l = sub2ind([(m-2) n],i,j);    
% $$$     trimx(l,sub2ind([m n],i+1,j)) = 1;
% $$$   end 
% $$$ end

[I,J] = ndgrid(1:m-2,1:n);
I=I(:);
J=J(:);

L = sub2ind([m-2 n],I,J);

trimx = sparse(L, sub2ind([m n],I+1,J  ), ones(size(L)), (m-2)*n , n*m);

%max(abs(trimx(:) - S(:)))