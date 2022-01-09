resp = [];
sc_prepare

for ii = g.prop.nfi:g.prop.nfr;
%for unit = Mcmb_l
%  Mcmb = unit{1};
  load(['R',num2str(version),'_',num2str(ii),'_',num2str(g.prop.nfr),'_',num2str(g.prop.iff*g.model_scale*100),'_',num2str(g.model_scale*100),'.h5'])
  sc_integrate
  rr0 = diag(response.range.pid/g.model_scale);
  rr1 = diag(response.range.prl);
  rr2 = diag(response.range.prr);
  rr = [rr0,rr1,rr2];
%  rr = [rr0,0*rr0+rr1,0*rr0+0*rr1+rr2];
  resp = cat(3,resp,rr);
%  pause
end
save(['resp_',num2str(g.prop.nfr),'_',num2str(g.prop.iff*g.model_scale*100),'_',num2str(g.model_scale*100)','.h5'],'resp','-hdf5')
