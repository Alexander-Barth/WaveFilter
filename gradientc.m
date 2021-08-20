function [hx_u,hx_v,hy_u,hy_v] = gradientc(h,m,n,dx,dy);

hx_u = (h(2:m,:) - h(1:m-1,:))/dx;
hy_v = (h(:,2:n) - h(:,1:n-1))/dy;

hx_v = zeros(m,n-1);
hx_v(2:m-1,:) = 0.25 * (hx_u(1:m-2,1:n-1) + hx_u(2:m-1,1:n-1) + hx_u(1:m-2,2:n) + hx_u(2:m-1,2:n));
hx_v(1,:) = hx_v(2,:);
hx_v(m,:) = hx_v(m-1,:);

hy_u = zeros(m-1,n);
hy_u(:,2:n-1) = 0.25 * (hy_v(1:m-1,1:n-2) + hy_v(1:m-1,2:n-1) + hy_v(2:m,1:n-2) + hy_v(2:m,2:n-1));
hy_u(:,1) = hy_u(:,2);
hy_u(:,n) = hy_u(:,n-1);

