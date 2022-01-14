resp = [];
doma = [];
sc_prepare

<<<<<<< HEAD
ii = g.prop.nfi;

IS_FILE = true;

while IS_FILE

#for ii = g.prop.nfi:g.prop.nfr;
%for unit = Mcmb_l
%  Mcmb = unit{1};
  try
     load(['R',num2str(converge),'_',num2str(ii),'_',num2str(g.prop.nfr),'_',num2str(g.prop.iff*g.model_scale*100),'_',num2str(g.model_scale*100),'.h5'])
     sc_integrate
     rr0 = diag(response.range.pid/g.model_scale);
     rr1 = diag(response.range.prl);
     rr2 = diag(response.range.prr);
     rr = [rr0,rr1,rr2];
%  rr = [rr0,0*rr0+rr1,0*rr0+0*rr1+rr2];
     dd = field.range.ptt;
     resp = cat(3,resp,rr);
     doma = cat(3,doma,dd);
   catch ME
     IS_FILE = false;
     rethrow(ME)
   end
   ii = ii + 1;
end
save(['doma_enh_',num2str(converge),'_',num2str(g.prop.nfr),'_',num2str(g.prop.iff*g.model_scale*100),'_',num2str(g.model_scale*100),'.h5'],'doma','-hdf5')
save(['resp_enh_',num2str(converge),'_',num2str(g.prop.nfr),'_',num2str(g.prop.iff*g.model_scale*100),'_',num2str(g.model_scale*100),'.h5'],'resp','-hdf5')
=======
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
>>>>>>> 2fd4097502e197651aaa605221aa26e315fd3029
