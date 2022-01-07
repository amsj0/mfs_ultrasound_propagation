%   DEFINES POSTPROCESSING PARAMETERS FUNCTIONS
f = @(x) 20*log10(abs(x));
% f = @(x) (imag(x));

cl = [-20 20];


%   DEFINES PISTON INDEXES 
ppt_per_surface = 1+floor(g.piston_radius*(Neltoverlambda/100));
%ppt_per_surface = 1;

%ndxI_vec = [num2cell(1:length(b.Ti.x),1),cell(1,length(b.To.x))];
%ndxI_mat = cellfun(@(x) circshift(ndxI_vec,x,2),num2cell(0:(ppt_per_surface-1)),'UniformOutput',0);
%ndxI_cat = cat(1,ndxI_mat{:});

ndx.T = fn_integrate_indexing(ppt_per_surface,length(b.Ti.x),length(b.To.x));
ndx.R = fn_integrate_indexing(ppt_per_surface,length(b.Ri.x),length(b.Ro.x));

%ndxO_vec = [cell(1,length(b.Ti.x)),num2cell(1:length(b.To.x),1)];
%ndxO_mat = cellfun(@(x) circshift(ndxO_vec,x,2),num2cell(0:(ppt_per_surface-1)),'UniformOutput',0);
%ndxO_cat = cat(1,ndxO_mat{:});

response.pitch.size = (length(T.a{:})-ppt_per_surface+1);
response.catch.size = (length(R.a{:})-ppt_per_surface+1);

strength.range.ptt = cell(response.pitch.size,response.catch.size);

response.range.pid = zeros(response.pitch.size,response.catch.size);
response.range.prr = zeros(response.pitch.size,response.catch.size);
response.range.prl = zeros(response.pitch.size,response.catch.size);
response.range.ptt = zeros(response.pitch.size,response.catch.size);

response.range.ptx = zeros(response.pitch.size,1);
response.range.prx = zeros(response.catch.size,1);

response.pitch.A = zeros(1,response.pitch.size);
response.catch.A = zeros(1,response.pitch.size);

response.pitch.I = T.z(~b.Ti.ndx)/nRD;
response.pitch.O = T.z(~b.To.ndx)/nRD;

response.catch.I = R.z(~b.Ri.ndx)/nRD;
response.catch.O = R.z(~b.Ro.ndx)/nRD;

response.ndx.RO = find(~b.Ro.ndx);
response.ndx.RI = find(~b.Ri.ndx);

response.ndx.TO = find(~b.To.ndx);
response.ndx.TI = find(~b.Ti.ndx);

draw_flag = 0;
%draw_flag = 1;
surf_flag = 0;
%surf_flag = 1;

if(draw_flag)
  b.A.zeros = zeros(size(b.A.z));

  figure(10),clf,prr = pcolor(b.A.x/nRD,b.A.z/nRD,b.A.zeros);
  % hold on,scatter(S.x/nRD,S.z/nRD,'r'),scatter(real(S.co/nRD),imag(S.co/nRD),'g'),scatter(real(S.ci/nRD),imag(S.ci/nRD),'r'),scatter(T.x(~b.To.ndx)/nRD,T.z(~b.To.ndx)/nRD,'*g'),scatter(T.x(~b.Ti.ndx)/nRD,T.z(~b.Ti.ndx)/nRD,'*r')
  axis equal, axis tight,caxis(cl),shading interp, colormap('gray'),title('Refracted')

  figure(11),clf,prl = pcolor(b.A.x/nRD,b.A.z/nRD,b.A.zeros);
  % hold on,scatter(S.x/nRD,S.z/nRD,'r'),scatter(real(S.co/nRD),imag(S.co/nRD),'g'),scatter(real(S.ci/nRD),imag(S.ci/nRD),'r'),scatter(T.x(~b.To.ndx)/nRD,T.z(~b.To.ndx)/nRD,'*g'),scatter(T.x(~b.Ti.ndx)/nRD,T.z(~b.Ti.ndx)/nRD,'*r')
  axis equal, axis tight,caxis(cl),shading interp, colormap('gray'),title('Reflected')

  figure(13),clf,ptt = pcolor(b.A.x/nRD,b.A.z/nRD,b.A.zeros);
  % hold on,scatter(T.x(~b.To.ndx)/nRD,T.z(~b.To.ndx)/nRD,'*g'),scatter(T.x(~b.Ti.ndx)/nRD,T.z(~b.Ti.ndx)/nRD,'*r'),scatter(R.x(~b.Ro.ndx)/nRD,R.z(~b.Ro.ndx)/nRD,'*g'),scatter(R.x(~b.Ri.ndx)/nRD,R.z(~b.Ri.ndx)/nRD,'*r')
  axis equal, axis tight,caxis(cl),shading interp, colormap('gray'),title('TotalFld')

  figure(14),clf,pid = pcolor(b.A.x/nRD,b.A.z/nRD,b.A.zeros);
  % hold on,scatter(S.x/nRD,S.z/nRD,'r'),scatter(real(S.co/nRD),imag(S.co/nRD),'g'),scatter(real(S.ci/nRD),imag(S.ci/nRD),'r'),scatter(T.x(~b.To.ndx)/nRD,T.z(~b.To.ndx)/nRD,'*g'),scatter(T.x(~b.Ti.ndx)/nRD,T.z(~b.Ti.ndx)/nRD,'*r')
  axis equal, axis tight,caxis(cl),shading interp, colormap('gray'),title('Incident')
  
  F_pid(response.pitch.size) = struct('cdata',[],'colormap',[]);
  F_prr(response.pitch.size) = struct('cdata',[],'colormap',[]);
  F_prl(response.pitch.size) = struct('cdata',[],'colormap',[]);
  F_ptt(response.pitch.size) = struct('cdata',[],'colormap',[]);
end



%%% TODO FIX ME
ndx.tst = reshape(1:(response.pitch.size*response.catch.size),response.pitch.size,response.catch.size);
[n1,n2] = ndgrid(1:response.pitch.size,1:response.catch.size);n0 = (mod(n2+n1+1,response.pitch.size)+1)';
ndx.mat = reshape(ndx.tst((n0(:)-0)+(n1(:)-1)*response.catch.size),response.pitch.size,response.catch.size)';

% pause
for tt = 1:response.pitch.size
        ndxI = [ndx.T.catI{:,ppt_per_surface-1+tt}];
        ndxO = [ndx.T.catO{:,ppt_per_surface-1+tt}];
        response.pitch.A(tt) = mean([response.pitch.I(ndxI),response.pitch.O(ndxO)]);
        ndx.T.i{tt} = ndxI;
        ndx.T.o{tt} = ndxO;
end

for rr = 1:response.catch.size
        ndxI = [ndx.R.catI{:,ppt_per_surface-1+rr}];
        ndxO = [ndx.R.catO{:,ppt_per_surface-1+rr}];
        response.catch.A(rr) = mean([response.catch.I(ndxI),response.catch.O(ndxO)]);
        ndx.R.i{rr} = ndxI;
        ndx.R.o{rr} = ndxO;
end
