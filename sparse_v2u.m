function v2u = sparse_v2u(m,n)

% $$$ 
% $$$ v2u = sparse( (m-1)*n, m*(n-1)  );
% $$$ 
% $$$ for i=1:m-1
% $$$   for j=2:n-1
% $$$     l = sub2ind([m-1 n],i,j);    
% $$$     
% $$$     v2u(l,sub2ind([m n-1],i,j-1)) = .25;
% $$$     v2u(l,sub2ind([m n-1],i,j)) = .25;
% $$$     v2u(l,sub2ind([m n-1],i+1,j-1)) = .25;
% $$$     v2u(l,sub2ind([m n-1],i+1,j)) = .25;
% $$$   end 
% $$$ 
% $$$    v2u(sub2ind([m-1 n],i,1),:) = v2u(sub2ind([m-1 n],i,2),:);
% $$$    v2u(sub2ind([m-1 n],i,n),:) = v2u(sub2ind([m-1 n],i,n-1),:);
% $$$ end


[I,J] = ndgrid(1:m-1,1:n);
I=I(:);
J=J(:);

L = sub2ind([m-1 n],I,J);

% special treatement for j=1 and j=n
J = max(J,2);
J = min(J,n-1);

val = 0.25 * ones(size(L));

v2u = sparse(L, sub2ind([m n-1],I  ,J-1), val, (m-1)*n, m*(n-1)) + ...
      sparse(L, sub2ind([m n-1],I  ,J  ), val, (m-1)*n, m*(n-1)) + ...
      sparse(L, sub2ind([m n-1],I+1,J-1), val, (m-1)*n, m*(n-1)) + ...
      sparse(L, sub2ind([m n-1],I+1,J  ), val, (m-1)*n, m*(n-1));



%max(abs(v2u(:) - S(:)))