function ext = sparse_ext(m,n)


% $$$ ext = sparse( (m+2)*(n+2) , n*m );
% $$$ 
% $$$ for j=1:n
% $$$   for i=1:m
% $$$     l = sub2ind([m+2 n+2],i+1,j+1);    
% $$$     ext(l,sub2ind([m n],i,j)) = 1;
% $$$   end 
% $$$ end

[I,J] = ndgrid(1:m,1:n);
I=I(:);
J=J(:);

ext = sparse(sub2ind([m+2 n+2],I+1,J+1), sub2ind([m n],I,J), ones(size(I)), (m+2)*(n+2) , n*m);

%max(abs(ext2(:) - ext(:)))