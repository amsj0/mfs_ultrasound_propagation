function p0 = fn_compute_reference_0_2(kc,k_cur,area_un)

% global fac m amp;

caller = evalin('caller','{g.fac0, g.amp}');
[fac,amp] = caller{:};

%         area = f0*RD;
z = k_cur*abs(kc);

%         bh0 = besselj(m+0,z);

% bj0 = amp(1)*besselj(0+0,z);
% by0 = amp(2)*1i*bessely(0+0,z);

%         bh0 = besselh(m+0,going,z);
% p0 = cat(2,fac(1)*bj0,fac(2)*by0);
% p0 = fac(1)*bj0 + fac(2)*by0;
p0 = area_un(1)*besselh(0,2,z);
%         K.p = K.p(K.c);