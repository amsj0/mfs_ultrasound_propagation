function [p,v,z] = fn_compute_field_inside_m2(varargin)

caller = evalin('caller','{g.m, g.hankel_kind}');
[m,kind] = caller{:};

R = varargin{1};
nR = varargin{2};
k_cur = varargin{3};
k_r = varargin{4};
dr = varargin{5};

len_R = length(R.c);

Ma = bsxfun(@minus,R.c,R.co.');

ny = nR.';
Mr1 = ny(:,ones(1,len_R));

rcos1 = (real(Ma).*cos(Mr1)+imag(Ma).*sin(Mr1))./abs(Ma);

z = k_cur.*abs(Ma);

bh0 = besselh(m+0,kind,z);
bh1 = 1/2*(besselh(m-1,kind,z)-besselh(m+1,kind,z)); % first derivate

Bh = bh0;
% DTy = by1.*rcos1*k_cur;
% % Dh = bh1.*rcos1*k_cur/k_r*dr;
Dh = bh1.*rcos1*k_cur/(k_r*dr);
% Dh = bh1.*rcos1;

p = Bh;
v = Dh;
end