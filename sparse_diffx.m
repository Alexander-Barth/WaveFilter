function diffx = sparse_diffx(m,n)

[I,J] = ndgrid(1:m-1,1:n);

L1 = [1:(m-1)*n]';

L2 = sub2ind([m n],I(:),J(:));
one = ones(size(L1));

diffx = sparse( ...
       [L1;     L1;   ],  ...
       [L2;     L2+1; ],  ...
       [-one;  one  ], (m-1)*n , m*n );

