function fn_matlab_to_h5(version,number,initial,final)

dataset = [num2str(number),'_',num2str(initial),'_',num2str(final)]

load(['R',num2str(version),'_',dataset,'.mat'])

save(['PP_',dataset,'.h5'],'T','R','S','Neltoverlambda','nRD','b','g','-hdf5')

unit = Mcmb_l;

for ii = g.prop.nfi:g.prop.nfr;
  Mcmb = unit{ii};
  save(['R',num2str(version),'_',num2str(ii),'_',num2str(g.prop.nfr),'_',num2str(g.prop.iff*g.model_scale*100),'_',num2str(g.model_scale*100),'.h5'],'Mcmb','-hdf5')
end
