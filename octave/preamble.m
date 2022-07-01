g.maxminspeed = 1.09645;
% g.maxminspeed = 1.0;

% g.enhance = g.pool_radius(1);
g.enhance = 1.1;
g.n = 3;
% Use this for saved analysis

% g.pool_radius(1) = 0.92;

g.pool_radius(1) = -1;
g.pool_radius(2) = 1;
g.pool_radius(3) = -.015;
g.pool_radius(4) = .015;

g.golden_ratio = 1;

g.grid_ratio = 4/1*0.0625;

g.prop.nfr = str2num(PROPE_NFR);
g.prop.iff = str2num(PROPE_IFF);
g.prop.nfi = str2num(PROPE_NFI);
g.prop.nff = str2num(PROPE_NFF);
g.prop.fff = str2num(PROPE_FFF);
% g.prop.nfr = 1;

g.model_scale = str2num(MODEL_SCL);

% g.prop.fr = 0.5e6;           %frequency
% g.prop.sfr = 0.5e6*linspace(1,1,g.prop.nfr);
% g.prop.rj = 1e3;            %water density
% g.prop.cj = 1481.44134805;   %speed of sound in water
g.prop.att = str2num(ATTEN_RAT);
% g.prop.fr = 2.5e5;             %frequency
g.prop.fr = g.model_scale*1e6;             %frequency
% g.prop.sfr = 2.5e5*linspace(1,1,g.prop.nfr);  
% g.prop.sfr = 1e6*linspace(1,1run,g.prop.nfr);  
g.prop.sfr = g.prop.fr*linspace(g.prop.iff,g.prop.fff,g.prop.nfr);  
g.prop.rj = 999.6150851557516; %water density
g.prop.cj = 1490;              %speed of sound in water

g.tiers = 2;
g.tiers_n = g.tiers - 1;

g.m = 0;
% g.m = 8;
% g.bessz = besselzero(g.m,g.tiers_n,1);
% g.bessz = besselzero(g.m,20-1,1);

% g.bessz = g.bessz([end-1,end]);
% g.bessz = g.bessz([1,end]);
% g.bessz = [0,g.bessz(end)];
% g.bessz = g.bessz(end);

g.prop.zed = g.prop.cj*g.prop.rj;              %acoustic specific impedance
g.prop.wav = g.prop.cj/g.prop.fr;              %acoustic wave number

% STATISTICS
oversqr = @(g,f)@(x,y) g(x,y)./f(x,y);
complex = @(g,f)@(x,y) g(x,y) + 1i*f(x,y);
complex2 = @(g,f)@(x) g(x) + 1i*f(x);

fvar{1} = @(x,y) abs(sum(conj(x).*y)).^2;
fvar{2} = @(x,y) sum(conj(x).*x).*sum(conj(y).*y);
correlat = oversqr(fvar{:});

% nrmse = @(x,y) rms(minus(x,y))./(max(abs(y))-min(abs(y)));
rmse = @(x,y) rms(minus(x,y));

relstd1 = @(x) std(x(1:end/2,:));
relstd2 = @(x) std(x((end/2+1):end,:));

% g.combinedstd = complex2(relstd1,relstd2);

% g.combined = complex(correlat,rmse);
% g.run_analysis = @(func,mode,cf,x,y) cell2mat(cellfun(@(f) (forfun(func,mode,f,x,y)),cf,'uni',0));
% g.run_analysis2 = @(func,mode,cf,x,y,z) cell2mat(cellfun(@(f) (forfun(func,mode,f,x,y,z)),cf,'uni',0));
% g.cn_filter = @(f) cellfun(@(x,y) and(~x,y),f,circshift(f,1,2),'uni',0);

g.hankel_kind = 2;
g.fac0 = [1 -1];
g.fac = [1 -1];
g.amp = [1 1];
% g.fac = [1 -1];
% g.mud = [1 1;1 1];
g.mud = [1 -1;1 -1];
g.tra = [1 1];
% g.tra = [1 -1];

g.TlLIM = 0;
g.BlLIM = -60;

% Recreate FEM solution
matchx = @(p,q) bsxfun(@(x,y) abs(x-y)<eps,p.x(:),q.x.');
matchz = @(p,q) bsxfun(@(x,y) abs(x-y)<eps,p.z(:),q.z.');

% g.recreate = @(p,q) and(matchx(p,q),matchz(p,q));

g.plot_flag = str2num(PLOT_FLAG);

g.extension = STRNG_EXT;
g.convergemod = CONVE_MOD;

g.centre_vector_x = str2num(MODEL_SCL)*(CENTX_VEC).';

g.centre_vector_y = str2num(MODEL_SCL)*(CENTY_VEC).';



g.piston_radius = str2num(MODEL_SCL)*str2num(PISTO_RAD);

g.piston_vector_pitch = str2num(MODEL_SCL)*str2num(PISTO_VCT);

g.piston_vector_catch = str2num(MODEL_SCL)*str2num(PISTO_VCR);

g.piston_distan = str2num(MODEL_SCL)*str2num(PISTO_DST);

g.interface_centre = 0.05;

g.piston_centre = str2num(MODEL_SCL)*str2num(PISTO_CEN);

g.scale = 4;

assignin('base','g',g);
