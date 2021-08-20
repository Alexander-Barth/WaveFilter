function u2v = sparse_u2v(m,n)


% $$$ u2v = sparse(  m*(n-1), (m-1)*n  );
% $$$ 
% $$$ for j=1:n-1
% $$$   for i=2:m-1
% $$$     l = sub2ind([m n-1],i,j);    
% $$$     
% $$$     u2v(l,sub2ind([m-1 n],i-1,j)) = .25;
% $$$     u2v(l,sub2ind([m-1 n],i,j)) = .25;
% $$$     u2v(l,sub2ind([m-1 n],i-1,j+1)) = .25;
% $$$     u2v(l,sub2ind([m-1 n],i,j+1)) = .25;
% $$$   end 
% $$$ 
% $$$    u2v(sub2ind([m n-1],1,j),:) = u2v(sub2ind([m n-1],2,j),:);
% $$$    u2v(sub2ind([m n-1],m,j),:) = u2v(sub2ind([m n-1],m-1,j),:);
% $$$ end


[I,J] = ndgrid(1:m,1:n-1);
I=I(:);
J=J(:);

L = sub2ind([m n-1],I,J);

% special treatement for i=1 and i=m
I = max(I,2);
I = min(I,m-1);

val = 0.25 * ones(size(L));

u2v = sparse(L, sub2ind([m-1 n],I-1,J  ), val, m*(n-1),(m-1)*n) + ...
      sparse(L, sub2ind([m-1 n],I  ,J  ), val, m*(n-1),(m-1)*n) + ...
      sparse(L, sub2ind([m-1 n],I-1,J+1), val, m*(n-1),(m-1)*n) + ...
      sparse(L, sub2ind([m-1 n],I  ,J+1), val, m*(n-1),(m-1)*n);



%max(abs(u2v(:) - S(:)))