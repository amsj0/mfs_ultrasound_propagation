pkg load parallel


path = 'D:\MATLAB\menisco\';
analysisfile = ['P.mat'];

load([path,analysisfile])

configfile = ['P',num2str(converge),num2str(modifier),'_',num2str(numbr_frequencies),'_',num2str(initi_frequencies),'_',num2str(final_frequencies),'_',num2str(skr),'.mat'];

load([path,configfile])

sc_prepare

datafile = [path,'R',num2str(converge),num2str(modifier),'_',num2str(g.prop.nfi),'_',num2str(g.prop.nfr),'_',num2str(g.prop.iff*g.model_scale*100),'_',num2str(g.model_scale*100),'_',num2str(skr),'.mat'];

input = load(datafile);

Mcmb = input.Mcmb;
sc_integrate

%r0 = diag(response.range.pid);
%r1 = diag(response.range.prl);
%r2 = diag(response.range.prr);

##r0 = response.range.pid(:);
##r1 = response.range.prl(:);
##r2 = response.range.prr(:);

resp0 = zeros(size(response.range.pid,2),g.prop.nfr);
doma0 = zeros(size(field.range.pid,2),g.prop.nfr);


resp = cell(3,size(response.range.pid,1));
doma = cell(3,size(field.range.pid,1));

%resp = cell(1,g.prop.nfr);
%doma = cell(1,g.prop.nfr);

%resp1 = cat(2,r0,r1,r2);
%doma1 = cat(3,field.range.prr,field.range.prl,field.range.pid);

resp(:) = {resp0};
doma(:) = {doma0};

path = 'D:\MATLAB\menisco\';

resp_str = [path,'resp_enh_',num2str(converge),num2str(modifier),'_',num2str(g.prop.nfr),'_',num2str(g.prop.iff*g.model_scale*100),'_',num2str(g.model_scale*100),'_',num2str(skr),'.hdf5']
doma_str = [path,'doma_enh_',num2str(converge),num2str(modifier),'_',num2str(g.prop.nfr),'_',num2str(g.prop.iff*g.model_scale*100),'_',num2str(g.model_scale*100),'_',num2str(skr),'.hdf5']

save(doma_str,'doma','-hdf5')
save(resp_str,'resp','-hdf5')
