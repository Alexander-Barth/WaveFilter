% [zetaf,uf,vf] = filter_var(domain,zetai,ui,vi,sigma)
% Variational filter for removing surface gravity waves
%
% Input:
%   domain: structure describing the model domain (see below)
%   zetai:  surface elevation on Arakawa C grid (array m x n)
%   ui:     depth averaged u-velocity on Arakawa C grid (array m-1 x n)
%   vi:     depth averaged v-velocity on Arakawa C grid (array m x n-1)
%   sigma:  scalar indicating the strength of bathymetric constraint.
%           Units of sigma are meters. A typical value of sigma is 50 m.
% Output:
%   zetaf:  filtered surface elevation on Arakawa C grid (array m x n)
%   uf:     filtered depth averaged u-velocity on Arakawa C grid (array m-1 x n)
%   vf:     filtered depth averaged v-velocity on Arakawa C grid (array m x n-1)
%
% The variable domain is a structure with the following fields:
%   domain.h: depth of the domain (array m x n). 
%   domain.pm: inverse of scale factor in 1st dimension (array m x n). For a 
%           regular grid: domain.pm = ones(m,n)/dx
%   domain.pn: inverse of scale factor in 2nd dimension (array m x n). For a 
%           regular grid: domain.pn = ones(m,n)/dy
%   domain.mask: 1 for sea points and 0 land points (array m x n)
%   domain.f: the average Coriolis freqency in s^-1 (scalar)
%   domain.g: the acceleration due to gravity in m s^-2 (scalar)
%
% See:
% Filtering inertia-gravity waves from the initial conditions of the linear shallow water equations
% Alexander Barth, Jean-Marie Beckers, Aida Alvera-Azcarate, Robert H. Weisberg
% Ocean Modelling

function [zetaf,uf,vf] = filter_var(domain,zetai,ui,vi,sigma)

% prescribe the values of the elevation at the boundary to 0
bc = 'dirichlet';

% convert sigma
sigma = 1/sigma^2;

% boundary values
zeta_bc = zeros(size(zetai));

q = pvort(domain,zetai,ui,vi);

[zetaf,res] = diffequation_var(domain,q,bc,zeta_bc,sigma);

[uf,vf] = geoflow(domain,zetaf);