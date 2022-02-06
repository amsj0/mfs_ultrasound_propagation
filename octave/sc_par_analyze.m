pkg load parallel

numbr_frequencies = 100
initi_frequencies = 1
final_frequencies = 100
modifier = 's';
converge = 'A';

path = 'D:\MATLAB\menisco\';
configfile = ['PP',num2str(modifier),'_',num2str(numbr_frequencies),'_',num2str(initi_frequencies),'_',num2str(final_frequencies),'.mat'];

load([path,configfile])

sc_prepare

datafile = [path,'R',num2str(converge),num2str(modifier),'_',num2str(g.prop.nfi),'_',num2str(g.prop.nfr),'_',num2str(g.prop.iff*g.model_scale*100),'_',num2str(g.model_scale*100),'.mat'];

input = load(datafile);

Mcmb = input.Mcmb;
sc_integrate

r0 = diag(response.range.pid);
r1 = diag(response.range.prl);
r2 = diag(response.range.prr);

resp0 = zeros(1,g.prop.nfr);
doma0 = zeros(size(field.range.pid,2),g.prop.nfr);


resp = cell(3,size(r0,1));
doma = cell(3,size(field.range.pid,1));

%resp = cell(1,g.prop.nfr);
%doma = cell(1,g.prop.nfr);

%resp1 = cat(2,r0,r1,r2);
%doma1 = cat(3,field.range.prr,field.range.prl,field.range.pid);

resp(:) = {resp0};
doma(:) = {doma0};

%resp = cell2mat(resp);
%doma = cell2mat(doma);

%resp = permute(resp,[3,2,1]);
%doma = permute(doma,[4,2,3,1]);

%resp = mat2cell(resp,ones(1,3),size(rr0,1),g.prop.nfr);
%doma = mat2cell(doma,ones(1,3),ones(1,size(field.range.pid,1)),size(field.range.pid,2),g.prop.nfr);

for ii = 1:size(doma,2)
   resp{1,ii}(1) = r0(ii);
   resp{2,ii}(1) = r1(ii);
   resp{3,ii}(1) = r2(ii);
   doma{1,ii}(:,1) = field.range.prr(ii,:);
   doma{2,ii}(:,1) = field.range.prl(ii,:);
   doma{3,ii}(:,1) = field.range.pid(ii,:);
end


clear input Mcmb r0 r1 r2

memory_factor = 4;

npc = nproc/memory_factor;

proc_tota = g.prop.nfr-g.prop.nfi;

proc_init = rem(proc_tota,npc);
proc_rema = fix(proc_tota/npc);

runs = cell(1,(proc_init!=0)+proc_rema);

runs(:) = {npc};

if (proc_init!=0)
   runs{1} =  proc_init;
end

%while IS_FILE

ii = g.prop.nfi+1;

datafiles = cell(1,npc);

for rr = 1:length(runs)
#(g.prop.nfi+1):npc:(g.prop.nfi+1);
  this_run = runs{rr};
  this_ser = 1:this_run;

  for nsl = this_ser
     datafiles{nsl} = ['R',num2str(converge),'_',num2str(ii+nsl-1),'_',num2str(g.prop.nfr),'_',num2str(g.prop.iff*g.model_scale*100),'_',num2str(g.model_scale*100),'.mat'];
  end

  out = parcellfun(this_run,@(x) fn_analyze(path,x,configfile),datafiles(this_ser));
  for jj = this_ser
    ij = ii+jj-1;
    rr1 = out(jj).rr;
    dd1 = out(jj).dd;
    for kk = 1:size(doma,2)
      resp{1,kk}(ij) = rr1(kk,1);
      resp{2,kk}(ij) = rr1(kk,2);
      resp{3,kk}(ij) = rr1(kk,3);
      doma{1,kk}(:,ij) = dd1(kk,:,1);
      doma{2,kk}(:,ij) = dd1(kk,:,2);
      doma{3,kk}(:,ij) = dd1(kk,:,3);
    end
  end
  ii = ii + this_run;
end
%{
for ii = (g.prop.nfi+2):4:(g.prop.nfr-3);
  display(ii+4)

  input{1} = ['R',num2str(converge),'_',num2str(ii+0),'_',num2str(g.prop.nfr),'_',num2str(g.prop.iff*g.model_scale*100),'_',num2str(g.model_scale*100),'.mat'];
  input{2} = ['R',num2str(converge),'_',num2str(ii+1),'_',num2str(g.prop.nfr),'_',num2str(g.prop.iff*g.model_scale*100),'_',num2str(g.model_scale*100),'.mat'];
  input{3} = ['R',num2str(converge),'_',num2str(ii+2),'_',num2str(g.prop.nfr),'_',num2str(g.prop.iff*g.model_scale*100),'_',num2str(g.model_scale*100),'.mat'];
  input{4} = ['R',num2str(converge),'_',num2str(ii+3),'_',num2str(g.prop.nfr),'_',num2str(g.prop.iff*g.model_scale*100),'_',num2str(g.model_scale*100),'.mat'];

  out = parcellfun(4,@(x) fn_analyze(x),input);
  for jj = 0:3
    ij = ii+jj;
    rr1 = out(jj+1).rr;
    dd1 = out(jj+1).dd;
    for kk = 1:size(doma,2)
      resp{1,kk}(ij) = rr1(kk,1);
      resp{2,kk}(ij) = rr1(kk,2);
      resp{3,kk}(ij) = rr1(kk,3);
      doma{1,kk}(:,ij) = dd1(kk,:,1);
      doma{2,kk}(:,ij) = dd1(kk,:,2);
      doma{3,kk}(:,ij) = dd1(kk,:,3);
   end
  end
end
%}
save([path,'doma_enh_',num2str(converge),'_',num2str(g.prop.nfr),'_',num2str(g.prop.iff*g.model_scale*100),'_',num2str(g.model_scale*100),'.h5'],'doma','-hdf5')
save([path,'resp_enh_',num2str(converge),'_',num2str(g.prop.nfr),'_',num2str(g.prop.iff*g.model_scale*100),'_',num2str(g.model_scale*100),'.h5'],'resp','-hdf5')
