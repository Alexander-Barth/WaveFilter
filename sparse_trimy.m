function trimy = sparse_trimy(m,n)


% $$$ trimy = sparse( (n-2)*m , n*m );
% $$$ 
% $$$ for j=1:n-2
% $$$   for i=1:m
% $$$     l = sub2ind([m n-2],i,j);    
% $$$     trimy(l,sub2ind([m n],i,j+1)) = 1;
% $$$   end 
% $$$ end

[I,J] = ndgrid(1:m,1:n-2);
I=I(:);
J=J(:);

L = sub2ind([m n-2],I,J);

trimy = sparse(L, sub2ind([m n],I,J+1), ones(size(L)), (n-2)*m , n*m);

%max(abs(trimy(:) - trimy2(:)))