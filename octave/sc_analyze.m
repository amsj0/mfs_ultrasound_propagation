resp = [];
doma = [];
sc_prepare

ii = g.prop.nfi;

load(['R',num2str(converge),'_',num2str(ii),'_',num2str(g.prop.nfr),'_',num2str(g.prop.iff*g.model_scale*100),'_',num2str(g.model_scale*100),'.mat'])

sc_integrate

rr0 = diag(response.range.pid);
rr1 = diag(response.range.prl);
rr2 = diag(response.range.prr);
rr = [rr0,rr1,rr2];

dd = field.range.prr+field.range.prl;

doma = zeros(size(dd,1),size(dd,2),g.prop.nfr);
resp = zeros(size(rr,1),size(rr,2),g.prop.nfr);

resp(:,:,ii) = rr;
doma(:,:,ii) = dd;

IS_FILE = true;

while IS_FILE

for ii = (g.prop.nfi+1):g.prop.nfr;
%for unit = Mcmb_l
%  Mcmb = unit{1};
  try
     display(ii)
     load(['R',num2str(converge),'_',num2str(ii),'_',num2str(g.prop.nfr),'_',num2str(g.prop.iff*g.model_scale*100),'_',num2str(g.model_scale*100),'.mat'])
     sc_integrate
     rr0 = diag(response.range.pid);
     rr1 = diag(response.range.prl);
     rr2 = diag(response.range.prr);
     rr = [rr0,rr1,rr2];
%  rr = [rr0,0*rr0+rr1,0*rr0+0*rr1+rr2];
     dd = field.range.prr+field.range.prl;
     resp(:,:,ii) = rr;
     doma(:,:,ii) = dd;
   catch ME
     IS_FILE = false;
     rethrow(ME)
   end
   ii = ii + 1;
end
save(['doma_enh_',num2str(converge),'_',num2str(g.prop.nfr),'_',num2str(g.prop.iff*g.model_scale*100),'_',num2str(g.model_scale*100),'.h5'],'doma','-hdf5')
save(['resp_enh_',num2str(converge),'_',num2str(g.prop.nfr),'_',num2str(g.prop.iff*g.model_scale*100),'_',num2str(g.model_scale*100),'.h5'],'resp','-hdf5')

