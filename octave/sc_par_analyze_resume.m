for rr = resume:length(runs)
#(g.prop.nfi+1):npc:(g.prop.nfi+1);
 this_run = runs{rr};
  
  this_ser = 1:this_run;

  for nsl = this_ser
     datafiles{nsl} = ['R',num2str(converge),num2str(modifier),'_',num2str(ii+nsl-1),'_',num2str(g.prop.nfr),'_',num2str(g.prop.iff*g.model_scale*100),'_',num2str(g.model_scale*100),'_',num2str(skr),'.mat'];
  end
  
  out = parcellfun(this_run,@(x) fn_analyze(path,x,configfile,analysisfile),datafiles(this_ser));
  for jj = this_ser
    ij = ii+jj-1;
    rr1 = out(jj).rr;
    dd1 = out(jj).dd;
    for kk = 1:size(resp,2)
      %this_kk = ceil(kk/size(doma,2));
      %that_kk = rem(kk,size(doma,2));
      resp{1,kk}(:,ij) = rr1(kk,:,1);
      resp{2,kk}(:,ij) = rr1(kk,:,2);
      resp{3,kk}(:,ij) = rr1(kk,:,3);     
      doma{2,kk}(:,ij) = dd1(kk,:,2);
      doma{3,kk}(:,ij) = dd1(kk,:,3);
      doma{1,kk}(:,ij) = dd1(kk,:,1);
    end
  end
  ii = ii + this_run;
end


save([path,'doma_enh_',num2str(converge),num2str(modifier),'_',num2str(g.prop.nfr),'_',num2str(g.prop.iff*g.model_scale*100),'_',num2str(g.model_scale*100),'_',num2str(skr),'.hdf5'],'doma','-hdf5')
save([path,'resp_enh_',num2str(converge),num2str(modifier),'_',num2str(g.prop.nfr),'_',num2str(g.prop.iff*g.model_scale*100),'_',num2str(g.model_scale*100),'_',num2str(skr),'.hdf5'],'resp','-hdf5')