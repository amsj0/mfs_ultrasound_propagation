function [Tj,Ty,Th,z] = fn_compute_field_boundary3(R,C,nR,k_cur,k_r,factor)

caller = evalin('caller','{g.fac, g.m, g.amp, g.hankel_kind}');
[fac,m,amp,kind] = caller{:};

len_C = length(C.c);
len_R = length(R.c);

Ma = bsxfun(@minus,factor*R.c,C.c);

ny = 0*pi+nR.';

Mr = ny(:,ones(1,len_C));
rcos = (real(Ma).*cos(Mr)+imag(Ma).*sin(Mr))./abs(Ma);
rsin = (real(Ma).*sin(Mr)-imag(Ma).*cos(Mr))./abs(Ma);

z = k_cur*abs(Ma);

ex0 = exp(1i*(m*angle(Ma)-C.p(ones(1,len_R),:)));

bj0 = amp(1)*besselj(m+0,z);
bj1 = amp(1)*1/2*(besselj(m-1,z)-besselj(m+1,z)); % first derivate

by0 = amp(2)*1i*bessely(m+0,z);
by1 = amp(2)*1i*1/2*(bessely(m-1,z)-bessely(m+1,z)); % first derivate

bh0 = amp(2)*besselh(m+0,kind,z);
bh1 = amp(2)*1/2*(besselh(m-1,kind,z)-besselh(m+1,kind,z)); % first derivate

% BVj = bj0.*ex0;
BVj = 1*bj0.*ex0*k_r;
% DVj = (bj1.*rcos*k_cur + 1i*bj0.*m.*rsin).*ex0;
DVj = (1*bj1.*rcos*k_cur + 1*1i*bj0./abs(Ma).*rsin*m).*ex0;
% DVj = (1*bj1.*rcos*k_cur + 1*bj0*1i.*m/RD.*rsin).*ex0;

% BVy = by0.*ex0;
BVy = 1*by0.*ex0*k_r;
% DVy = (by1.*rcos*k_cur + 1i*by0.*m.*rsin).*ex0;
DVy = (1*by1.*rcos*k_cur + 1*1i*by0./abs(Ma).*rsin*m).*ex0;

% BVy = by0.*ex0;
BVh = 1*bh0.*ex0*k_r;
% DVy = (by1.*rcos*k_cur + 1i*by0.*m.*rsin).*ex0;
DVh = (1*bh1.*rcos*k_cur + 1*1i*bh0./abs(Ma).*rsin*m).*ex0;

Tj = fac(1)*cat(1,1*BVj,1*DVj);
Ty = fac(2)*cat(1,1*BVy,1*DVy);
Th = cat(1,1*BVh,1*DVh);