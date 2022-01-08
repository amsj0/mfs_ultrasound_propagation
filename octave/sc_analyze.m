resp = [];
sc_prepare

for ii = g.prop.nfi:g.prop.nfr;
%for unit = Mcmb_l
%  Mcmb = unit{1};
  load(['RR_',num2str(ii),'_',num2str(g.prop.nfr),'_',num2str(g.prop.iff*g.model_scale*100),'_',num2str(g.model_scale*100),'.mat'])
  sc_integrate
  rr = diag(response.range.pid/g.model_scale+(response.range.prr+response.range.prl));
  resp = [resp,rr];
%  pause
end
save('resp_',num2str(g.prop.nfr),'_',num2str(g.prop.iff*g.model_scale*100),'_',num2str(g.model_scale*100)','.mat','resp','-hdf5')