function [Tj,Ty,Th] = fn_compute_propagator_outside3(varargin)

% global T
D = evalin('caller','b.D');

caller = evalin('caller','{g.m,g.hankel_kind}');
[m,kind] = caller{:};

R = varargin{1};
nR = varargin{2};
area_un = varargin{3};
k_cur = varargin{4};
k_r = varargin{5};
len_O = varargin{6};
dr =  varargin{7};

Ma = bsxfun(@minus,D.co,R.c.');

% Ar = area_un(ones(len_O,1),:);
Ar = area_un(1);

ny = nR;

Mr = ny(ones(length(D.co),1),:);
rcos = (real(Ma).*cos(Mr)+imag(Ma).*sin(Mr))./abs(Ma);

% TESTING THESE
z = k_cur*abs(Ma);
% z = k_r*abs(Ma);

% bj0 = besselj(m+0,z);
% bj1 = 1/2*(besselj(m-1,z)-besselj(m+1,z)); % first derivate

% by0 = 1i*bessely(m+0,z);
% by1 = 1i*1/2*(bessely(m-1,z)-bessely(m+1,z)); % first derivate

bh0 = besselh(m+0,kind,z);
bh1 = 1/2*(besselh(m-1,kind,z)-besselh(m+1,kind,z)); % first derivate

% BUj = 1i/4*Ar.*bj0;
% AUj = -1i/4*Ar.*bj1.*rcos*k_cur;
% AUj = 1i/4*Ar.*bj1.*rcos*k_cur/k_r*dr;

% BUy = 1i/4*Ar.*by0;
% AUy = -1i/4*Ar.*by1.*rcos*k_cur;
% AUy = 1i/4*Ar.*by1.*rcos*k_cur/k_r*dr;

BUh = 1i/4*Ar.*bh0;
% AUy = -1i/4*Ar.*by1.*rcos*k_cur;
AUh = 1i/4*Ar.*bh1.*rcos*k_cur/k_r*dr;

% Tj = cat(2,AUj,-BUj);
% Ty = -cat(2,AUy,-BUy);
Tj = [];
Ty = [];
Th = cat(2,AUh,-BUh);