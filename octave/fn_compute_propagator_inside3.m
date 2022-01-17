function [Tj,Ty,Th] = fn_compute_propagator_inside3(varargin)

% global T
D = evalin('caller','b.D');

R = varargin{1};
nR = varargin{2};
area_un = varargin{3};
k_cur = varargin{4};
k_r = varargin{5};
len_I = varargin{6};
dr =  varargin{7};

Ma = bsxfun(@minus,D.ci,R.c.');

% Ar = area_un(ones(len_I,1),:);
Ar = area_un(1);

nx = nR;

Mr = nx(ones(length(D.ci),1),:);
rcos = (real(Ma).*cos(Mr)+imag(Ma).*sin(Mr))./abs(Ma);

% TESTING THESE
z = k_cur*abs(Ma);
% z = k_r*abs(Ma);

bj0 = besselj(0+0,z);
bj1 = 1/2*(besselj(0-1,z)-besselj(0+1,z));

by0 = 1i*bessely(0+0,z);
by1 = 1i*1/2*(bessely(0-1,z)-bessely(0+1,z));

bh0 = besselh(0+0,2,z);
bh1 = 1/2*(besselh(0-1,2,z)-besselh(0+1,2,z));

BUj = 1i/4*Ar.*bj0;
% AUj = -1i/4*Ar.*bj1.*rcos*k_cur;
AUj = 1i/4*Ar.*bj1.*rcos*k_cur/k_r*dr;

BUy = 1i/4*Ar.*by0;
% AUy = -1i/4*Ar.*by1.*rcos*k_cur;
AUy = 1i/4*Ar.*by1.*rcos*k_cur/k_r*dr;

BUh = 1i/4*Ar.*bh0;
% AUy = -1i/4*Ar.*by1.*rcos*k_cur;
AUh = 1i/4*Ar.*bh1.*rcos*k_cur/k_r*dr;

Tj = cat(2,AUj,-BUj);
Ty = -cat(2,AUy,-BUy);
Th = cat(2,AUh,-BUh);