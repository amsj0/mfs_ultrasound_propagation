%   DEFINES POSTPROCESSING PARAMETERS FUNCTIONS
f = @(x) 20*log10(abs(x));
% f = @(x) (imag(x));

b.A.zeros = zeros(size(b.A.z));

cl = [-2 2];

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

%   DEFINES PISTON INDEXES 
ppt_per_surface = 1+floor(10000*g.piston_radius/(Neltoverlambda));
% ppt_per_surface = 17;

ndxI_vec = [num2cell(1:length(b.Ti.x),1),cell(1,length(b.To.x))];
ndxI_mat = cellfun(@(x) circshift(ndxI_vec,x,2),num2cell(0:(ppt_per_surface-1)),'UniformOutput',0);
ndxI_cat = cat(1,ndxI_mat{:});

ndxO_vec = [cell(1,length(b.Ti.x)),num2cell(1:length(b.To.x),1)];
ndxO_mat = cellfun(@(x) circshift(ndxO_vec,x,2),num2cell(0:(ppt_per_surface-1)),'UniformOutput',0);
ndxO_cat = cat(1,ndxO_mat{:});

response_size = (length(ndxI_vec)-ppt_per_surface+1);

strength_range_ptt = cell(response_size,response_size);

response_range.pid = zeros(response_size,response_size);
response_range.prr = zeros(response_size,response_size);
response_range.prl = zeros(response_size,response_size);
response_range.ptt = zeros(response_size,response_size);

response_range.ptx = zeros(response_size,1);
response_range.prx = zeros(response_size,1);

response_domain = zeros(1,response_size);
response_distaI = T.z(~b.Ti.ndx)/nRD;
response_distaO = T.z(~b.To.ndx)/nRD;

response_ndx.RO = find(~b.Ro.ndx);
response_ndx.RI = find(~b.Ri.ndx);

response_ndx.TO = find(~b.To.ndx);
response_ndx.TI = find(~b.Ti.ndx);

draw_flag = 0;
% draw_flag = 1;
surf_flag = 0;
% surf_flag = 1;

F_pid(response_size) = struct('cdata',[],'colormap',[]);
F_prr(response_size) = struct('cdata',[],'colormap',[]);
F_prl(response_size) = struct('cdata',[],'colormap',[]);
F_ptt(response_size) = struct('cdata',[],'colormap',[]);

ndx.tst = reshape(1:(response_size*response_size),response_size,response_size);
[n1,n2] = ndgrid(1:response_size,1:response_size);n0 = (mod(n2+n1+1,response_size)+1)';
ndx.mat = reshape(ndx.tst((n0(:)-0)+(n1(:)-1)*response_size),response_size,response_size)';

% pause
for tt = 1:response_size
        ndx.i{tt} = [ndxI_cat{:,ppt_per_surface-1+tt}];
        ndx.o{tt} = [ndxO_cat{:,ppt_per_surface-1+tt}];
end
for tt = 1:response_size
    
    %   EXTRACT TRANSMITTER PISTON INDEXES
%     ndx.TI = [ndxI_cat{:,ppt_per_surface-1+tt}];
%     ndx.TO = [ndxO_cat{:,ppt_per_surface-1+tt}];
    ndx.TI = ndx.i{tt};
    ndx.TO = ndx.o{tt};
    
    resTO = response_ndx.TO(ndx.TO);
    resTI = response_ndx.TI(ndx.TI);
    
    %   FILTER FIELD PARAMETERS
    b.D.prr(ndx0,:)  = sum(Mcmb{1,1}(:,ndx.TO),2);  % APPLIED MODEL
    b.D.prr(~ndx0,:) = sum(Mcmb{2,1}(:,ndx.TI),2);  % APPLIED MODEL 
    b.D.prl(ndx0,:)  = sum(Mcmb{3,1}(:,ndx.TI),2);  % APPLIED MODEL 
    b.D.prl(~ndx0,:) = sum(Mcmb{4,1}(:,ndx.TO),2);  % APPLIED MODEL 
    b.D.pid(ndx0,:)  = sum(Mcmb{5,1}(:,ndx.TI),2);  % APPLIED MODEL 
    b.D.pid(~ndx0,:) = sum(Mcmb{6,1}(:,ndx.TO),2);  % APPLIED MODEL 
    b.D.ptt = b.D.pid + b.D.prl + b.D.prr;
    
    response_domain(tt) = mean([response_distaI(ndx.TI),response_distaO(ndx.TO)]);
    
    
    b.T.prr{1} = sum(Mcmb{1,2}(:,ndx.TO),2);  % APPLIED MODEL
    b.T.prr{2} = sum(Mcmb{2,2}(:,ndx.TI),2);  % APPLIED MODEL 
    b.T.prl{1} = sum(Mcmb{3,2}(:,ndx.TI),2);  % APPLIED MODEL 
    b.T.prl{2} = sum(Mcmb{4,2}(:,ndx.TO),2);  % APPLIED MODEL 
    b.T.pid{1} = sum(Mcmb{5,2}(:,ndx.TI),2);  % APPLIED MODEL 
    b.T.pid{2} = sum(Mcmb{6,2}(:,ndx.TO),2);  % APPLIED MODEL    
    b.T.ptt{1} = (b.T.pid{1} + b.T.prl{1} + b.T.prr{1});
    b.T.ptt{2} = (b.T.pid{2} + b.T.prl{2} + b.T.prr{2});
    %   EXTRACT RECEIVER PISTON INDEXES
    for rr = 1:response_size
%     for rr = tt
        
%         ndx.RI = [ndxI_cat{:,ppt_per_surface-1+rr}];
%         ndx.RO = [ndxO_cat{:,ppt_per_surface-1+rr}];
        ndx.RI = ndx.i{rr};
        ndx.RO = ndx.o{rr};
        
%         resRO{jj} = response_ndx.RO(ndx.RO);
%         resRI{jj} = response_ndx.RI(ndx.RI);
        resRO = response_ndx.RO(ndx.RO);
        resRI = response_ndx.RI(ndx.RI);        
        %   EXTRACT PISTON CENTRE
        
%         if((jj+response_size)<ii)
%             strength = strength_range_ptt{jj,ii}.';
%         else
%             strength = ones(length(ndxI_vec),1);
%         end
%         b.R.pid = sparse(response_size,response_size);
%         b.R.prr = sparse(response_size,response_size);
%         b.R.prl = sparse(response_size,response_size);
%         b.R.ptt = sparse(response_size,response_size);
        %   FILTER RESPONSE PARAMETERS    
%         b.R.prr(resRI{jj}) = sum(Mcmb{1,2}(ndx.RI,ndx.TO)*strength(resRO{ii}),2);  % APPLIED MODEL
%         b.R.prr(resRO{jj}) = sum(Mcmb{2,2}(ndx.RO,ndx.TI)*strength(resRI{ii}),2);  % APPLIED MODEL 
%         b.R.prl(resRI{jj}) = sum(Mcmb{3,2}(ndx.RI,ndx.TI)*strength(resRI{ii}),2);  % APPLIED MODEL 
%         b.R.prl(resRO{jj}) = sum(Mcmb{4,2}(ndx.RO,ndx.TO)*strength(resRO{ii}),2);  % APPLIED MODEL 
%         b.R.pid(resRI{jj}) = sum(Mcmb{5,2}(ndx.RI,ndx.TI)*strength(resRI{ii}),2);  % APPLIED MODEL 
%         b.R.pid(resRO{jj}) = sum(Mcmb{6,2}(ndx.RO,ndx.TO)*strength(resRO{ii}),2);  % APPLIED MODEL     
%         b.R.prr(resRI) = sum(Mcmb{1,2}(ndx.RI,ndx.TO),2);  % APPLIED MODEL
%         b.R.prr(resRO) = sum(Mcmb{2,2}(ndx.RO,ndx.TI),2);  % APPLIED MODEL 
%         b.R.prl(resRI) = sum(Mcmb{3,2}(ndx.RI,ndx.TI),2);  % APPLIED MODEL 
%         b.R.prl(resRO) = sum(Mcmb{4,2}(ndx.RO,ndx.TO),2);  % APPLIED MODEL 
%         b.R.pid(resRI) = sum(Mcmb{5,2}(ndx.RI,ndx.TI),2);  % APPLIED MODEL 
%         b.R.pid(resRO) = sum(Mcmb{6,2}(ndx.RO,ndx.TO),2);  % APPLIED MODEL     
        b.R.prr(resRI) = b.T.prr{1}(ndx.RI); % APPLIED MODEL
        b.R.prr(resRO) = b.T.prr{2}(ndx.RO); % APPLIED MODEL
        b.R.prl(resRI) = b.T.prl{1}(ndx.RI); % APPLIED MODEL
        b.R.prl(resRO) = b.T.prl{2}(ndx.RO); % APPLIED MODEL
        b.R.pid(resRI) = b.T.pid{1}(ndx.RI); % APPLIED MODEL
        b.R.pid(resRO) = b.T.pid{2}(ndx.RO); % APPLIED MODEL        
        b.R.ptt = (b.R.pid + b.R.prl + b.R.prr);
%         b.R.prr(:,resTI) = sum(Mcmb{2,2}(:,ndx.TI),2);  % APPLIED MODEL 
%         b.R.prl(:,resTI) = sum(Mcmb{3,2}(:,ndx.TI),2);  % APPLIED MODEL 
%         b.R.prl(:,resTO) = sum(Mcmb{4,2}(:,ndx.TO),2);  % APPLIED MODEL 
%         b.R.pid(:,resTI) = sum(Mcmb{5,2}(:,ndx.TI),2);  % APPLIED MODEL 
%         b.R.pid(:,resTO) = sum(Mcmb{6,2}(:,ndx.TO),2);  % APPLIED MODEL 

%         b.R.ptt = (b.R.pid + b.R.prl + b.R.prr)*strength([resRI{ii},resRO{ii}]);
%         b.R.ptt = (b.R.pid + b.R.prl + b.R.prr);
        
        %   UPDATE FIGURES


    %     b.A.p(~b.D.ndx) = (b.A.p0*(max(abs(b.D.ptt(:)))/max(abs(b.A.p0(:)))));
    %     b.A.p(~b.D.ndx) = (b.A.p0);
    %     b.A.p(isnan(f(b.A.p)))=0;
    %     phf.CData = f(b.A.p);

        %   ANALYSE RESPONSE RANGES
%         strength_range_ptt{tt,rr} = b.R.ptt;
%         if ((tt<(response_size/2) && rr<(response_size/2)) || (tt>(response_size/2) && rr>(response_size/2)))
%             fx = fx1;
%         else
%             fx = fx2;
%         end
        response_range.pid(tt,rr) = (sum(b.R.pid));
        response_range.prr(tt,rr) = (sum(b.R.prr));
        response_range.prl(tt,rr) = (sum(b.R.prl));
        response_range.ptt(tt,rr) = (sum(b.R.ptt)); 
        response_range.prx(rr)    = ~(mod(rr,2*ppt_per_surface))*(sum(b.R.ptt)); 
%         response_range.pid(tt,rr) = fx(sum(b.R.pid));
%         response_range.prr(tt,rr) = fx(sum(b.R.prr));
%         response_range.prl(tt,rr) = fx(sum(b.R.prl));
%         response_range.ptt(tt,rr) = fx(sum(b.R.ptt)); 
%         response_range.prx(rr)    = ~(mod(rr,2*ppt_per_surface))*fx(sum(b.R.ptt)); 
%         response_range_pid(ii,jj) = sum(abs(b.R.pid));
%         response_range_prr(ii,jj) = sum(abs(b.R.prr));
%         response_range_prl(ii,jj) = sum(abs(b.R.prl));
%         response_range_ptt(ii,jj) = sum(abs(b.R.ptt));  
        
        %   CLEAR RESPONSE RANGES
%         b = rmfield(b,'R');

        
%         if(ii==ceil(response_size/2))
%             pause
%         end
    end
%     if(rr==floor(response_size/2))
    response_range.ptx(tt) = (sum(response_range.prx)); 
    if(draw_flag)
            b.A.p(~b.D.ndx) = b.D.prr;
            b.A.p(isnan(f(b.A.p)))=0;
            prr.CData = f(b.A.p);

            b.A.p(~b.D.ndx) = b.D.prl;
            b.A.p(isnan(f(b.A.p)))=0;
            prl.CData = f(b.A.p);

            b.A.p(~b.D.ndx) = b.D.pid;
            b.A.p(isnan(f(b.A.p)))=0;
            pid.CData = f(b.A.p);

            b.A.p(~b.D.ndx) = b.D.ptt;
            b.A.p(isnan(f(b.A.p)))=0;
            ptt.CData = f(b.A.p);
            drawnow

            F_pid(tt) = getframe(pid.Parent);
            F_prr(tt) = getframe(prr.Parent);
            F_prl(tt) = getframe(prl.Parent);
            F_ptt(tt) = getframe(ptt.Parent);
    end
%     end
end

if(surf_flag)
    figure(5),surf(response_domain,response_domain,f(response_range.pid/2))
    title('Incident');axis equal, shading flat,axis tight,view(2), colormap('gray')

    figure(6),surf(response_domain,response_domain,f(response_range.prl/2))
    title('Reflected');axis equal, shading flat,axis tight,view(2), colormap('gray')

    figure(7),surf(response_domain,response_domain,f(response_range.prr/2))
    title('Refracted');axis equal, shading flat,axis tight,view(2), colormap('gray')

    figure(8),surf(response_domain,response_domain,f(response_range.ptt/2))
    title('Total Field');axis equal, shading flat,axis tight,view(2), colormap('gray')

    figure(70),surf(response_domain,response_domain, f((abs(tril(response_range.prr/2,-1)) - abs(triu(response_range.prr/2,1)'))/max(abs(response_range.prr(:)/2))))
    title('Diff Refracted');axis equal, shading flat,axis tight,view(2), colormap('gray')

    figure(71),surf(response_domain,response_domain, f((abs(tril(response_range.prl/2,-1)) - abs(triu(response_range.prl/2,1)'))/max(abs(response_range.prl(:)/2))))
    title('Diff Reflected');axis equal, shading flat,axis tight,view(2), colormap('gray')
    
    figure(72),surf(response_domain,response_domain, f((abs(tril(response_range.pid/2,-1)) - abs(triu(response_range.pid/2,1)'))/max(abs(response_range.pid(:)/2))))
    title('Diff Incident');axis equal, shading flat,axis tight,view(2), colormap('gray')    
    
    figure(73),surf(response_domain,response_domain, f((abs(tril(response_range.ptt/2,-1)) - abs(triu(response_range.ptt/2,1)'))/max(abs(response_range.ptt(:)/2))))
    title('Diff TotalFld');axis equal, shading flat,axis tight,view(2), colormap('gray')
end


% figure(12),clf,pcolor(b.A.x/nRD,b.A.z/nRD,(f(bAp0*(max(abs(b.D.ptt(:)))/max(abs(bAp0(:))))))),axis equal
% figure(12),clf,pcolor(b.A.x/nRD,b.A.z/nRD,(f(b.A.p0))),axis equal
% hold on,scatter(R.x/nRD,R.z/nRD,'r'),scatter(real(S.co/nRD),imag(S.co/nRD),'g'),scatter(real(S.ci/nRD),imag(S.ci/nRD),'r'),scatter(T.x(~b.To.ndx)/nRD,T.z(~b.To.ndx)/nRD,'*g'),scatter(T.x(~b.Ti.ndx)/nRD,T.z(~b.Ti.ndx)/nRD,'*r')
% caxis([-2 2]),shading flat, colormap('gray'),title('HarmFEM')