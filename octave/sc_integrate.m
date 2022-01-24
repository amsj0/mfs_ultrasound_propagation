for tt = 1:response.pitch.size
    
    %   EXTRACT TRANSMITTER PISTON INDEXES
%     ndx.TI = [ndxI_cat{:,ppt_per_surface-1+tt}];
%     ndx.TO = [ndxO_cat{:,ppt_per_surface-1+tt}];
    ndx.TI = ndx.T.i{tt};
    ndx.TO = ndx.T.o{tt};
    
    resTO = response.ndx.TO(ndx.TO);
    resTI = response.ndx.TI(ndx.TI);
    
    %   FILTER FIELD PARAMETERS
    b.D.prr(b.D.ndx0,:)  = sum(Mcmb{1,1}(:,ndx.TO),2);  % APPLIED MODEL
    b.D.prr(~b.D.ndx0,:) = sum(Mcmb{2,1}(:,ndx.TI),2);  % APPLIED MODEL 
    b.D.prl(b.D.ndx0,:)  = sum(Mcmb{3,1}(:,ndx.TI),2);  % APPLIED MODEL 
    b.D.prl(~b.D.ndx0,:) = sum(Mcmb{4,1}(:,ndx.TO),2);  % APPLIED MODEL 
    b.D.pid(b.D.ndx0,:)  = sum(Mcmb{5,1}(:,ndx.TI),2);  % APPLIED MODEL 
    b.D.pid(~b.D.ndx0,:) = sum(Mcmb{6,1}(:,ndx.TO),2);  % APPLIED MODEL 
    b.D.ptt = b.D.pid + b.D.prl + b.D.prr;
    
%    response.pitch.A(tt) = mean([response.pitch.I(ndx.TI),response.pitch.O(ndx.TO)]);
    
    
    b.T.prr{1} = sum(Mcmb{1,2}(:,ndx.TO),2);  % APPLIED MODEL
    b.T.prr{2} = sum(Mcmb{2,2}(:,ndx.TI),2);  % APPLIED MODEL 
    b.T.prl{1} = sum(Mcmb{3,2}(:,ndx.TI),2);  % APPLIED MODEL 
    b.T.prl{2} = sum(Mcmb{4,2}(:,ndx.TO),2);  % APPLIED MODEL 
    b.T.pid{1} = sum(Mcmb{5,2}(:,ndx.TI),2);  % APPLIED MODEL 
    b.T.pid{2} = sum(Mcmb{6,2}(:,ndx.TO),2);  % APPLIED MODEL    
    b.T.ptt{1} = (b.T.pid{1} + b.T.prl{1} + b.T.prr{1});
    b.T.ptt{2} = (b.T.pid{2} + b.T.prl{2} + b.T.prr{2});
    %   EXTRACT RECEIVER PISTON INDEXES
    for rr = 1:response.catch.size
%     for rr = tt
        
%         ndx.RI = [ndxI_cat{:,ppt_per_surface-1+rr}];
%         ndx.RO = [ndxO_cat{:,ppt_per_surface-1+rr}];
        ndx.RI = ndx.R.i{rr};
        ndx.RO = ndx.R.o{rr};
        
%         resRO{jj} = response_ndx.RO(ndx.RO);
%         resRI{jj} = response_ndx.RI(ndx.RI);
        resRO = response.ndx.RO(ndx.RO);
        resRI = response.ndx.RI(ndx.RI);        
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
%        b.R.prr(resRI) = b.T.prr{1}(ndx.RI); % APPLIED MODEL
%        b.R.prr(resRO) = b.T.prr{2}(ndx.RO); % APPLIED MODEL
%        b.R.prl(resRI) = b.T.prl{1}(ndx.RI); % APPLIED MODEL
%        b.R.prl(resRO) = b.T.prl{2}(ndx.RO); % APPLIED MODEL
%        b.R.pid(resRI) = b.T.pid{1}(ndx.RI); % APPLIED MODEL
%        b.R.pid(resRO) = b.T.pid{2}(ndx.RO); % APPLIED MODEL        
%        b.R.ptt = (b.R.pid + b.R.prl + b.R.prr);
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
%        response.range.pid(tt,rr) = (sum(b.R.pid)); b.R.pid = [];
%        response.range.prr(tt,rr) = (sum(b.R.prr)); b.R.prr = [];
%        response.range.prl(tt,rr) = (sum(b.R.prl)); b.R.prl = [];
%        response.range.ptt(tt,rr) = (sum(b.R.ptt)); b.R.ptt = [];
        response.range.pid(tt,rr) = (sum(b.T.pid{1}(ndx.RI)) + sum(b.T.pid{2}(ndx.RO)));
        response.range.prr(tt,rr) = (sum(b.T.prr{1}(ndx.RI)) + sum(b.T.prr{2}(ndx.RO)));
        response.range.prl(tt,rr) = (sum(b.T.prl{1}(ndx.RI)) + sum(b.T.prl{2}(ndx.RO)));
%        response.range.ptt(tt,rr) = (sum(b.T.ptt{1}(ndx.RI)) + sum(b.T.ptt{2}(ndx.RO)));
%        response.range.prx(rr)    = ~(mod(rr,2*ppt_per_surface))*(sum(b.R.ptt)); 
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
    response.range.ptx(tt) = (sum(response.range.prx)); 
    if(draw_flag)
            b.A.p(~b.D.ndx) = b.D.prr;
            b.A.p(isnan(f(b.A.p)))=0;
            set(prr,'CData',f(b.A.p));

            b.A.p(~b.D.ndx) = b.D.prl;
            b.A.p(isnan(f(b.A.p)))=0;
            set(prl,'CData',f(b.A.p));

            b.A.p(~b.D.ndx) = b.D.pid;
            b.A.p(isnan(f(b.A.p)))=0;
            set(pid,'CData',f(b.A.p));

            b.A.p(~b.D.ndx) = b.D.ptt;
            b.A.p(isnan(f(b.A.p)))=0;
            set(ptt,'CData',f(b.A.p));
            drawnow
            if tt == floor(response.pitch.size/2)
                pause
            end

            F_pid(tt) = getframe(get(pid,'Parent'));
            F_prr(tt) = getframe(get(prr,'Parent'));
            F_prl(tt) = getframe(get(prl,'Parent'));
            F_ptt(tt) = getframe(get(ptt,'Parent'));
    end
%     end
end

response.range.ptt = response.range.prl+response.range.prr+response.range.pid;

if(surf_flag)
    figure(5),surf(response.pitch.A,response.catch.A,f(response.range.pid/2))
    title('Incident');axis equal, shading flat,axis tight,view(2), colormap('gray')

    figure(6),surf(response.pitch.A,response.catch.A,f(response.range.prl/2))
    title('Reflected');axis equal, shading flat,axis tight,view(2), colormap('gray')

    figure(7),surf(response.pitch.A,response.catch.A,f(response.range.prr/2))
    title('Refracted');axis equal, shading flat,axis tight,view(2), colormap('gray')

    figure(8),surf(response.pitch.A,response.catch.A,f(response.range.ptt/2))
    title('Total Field');axis equal, shading flat,axis tight,view(2), colormap('gray')

    figure(70),surf(response.pitch.A,response.catch.A, f((abs(tril(response.range.prr/2,-1)) - abs(triu(response.range.prr/2,1)'))/max(abs(response.range.prr(:)/2))))
    title('Diff Refracted');axis equal, shading flat,axis tight,view(2), colormap('gray')

    figure(71),surf(response.pitch.A,response.catch.A, f((abs(tril(response.range.prl/2,-1)) - abs(triu(response.range.prl/2,1)'))/max(abs(response.range.prl(:)/2))))
    title('Diff Reflected');axis equal, shading flat,axis tight,view(2), colormap('gray')
    
    figure(72),surf(response.pitch.A,response.catch.A, f((abs(tril(response.range.pid/2,-1)) - abs(triu(response.range.pid/2,1)'))/max(abs(response.range.pid(:)/2))))
    title('Diff Incident');axis equal, shading flat,axis tight,view(2), colormap('gray')    
    
    figure(73),surf(response.pitch.A,response.catch.A, f((abs(tril(response.range.ptt/2,-1)) - abs(triu(response.range.ptt/2,1)'))/max(abs(response.range.ptt(:)/2))))
    title('Diff TotalFld');axis equal, shading flat,axis tight,view(2), colormap('gray')
end


% figure(12),clf,pcolor(b.A.x/nRD,b.A.z/nRD,(f(bAp0*(max(abs(b.D.ptt(:)))/max(abs(bAp0(:))))))),axis equal
% figure(12),clf,pcolor(b.A.x/nRD,b.A.z/nRD,(f(b.A.p0))),axis equal
% hold on,scatter(R.x/nRD,R.z/nRD,'r'),scatter(real(S.co/nRD),imag(S.co/nRD),'g'),scatter(real(S.ci/nRD),imag(S.ci/nRD),'r'),scatter(T.x(~b.To.ndx)/nRD,T.z(~b.To.ndx)/nRD,'*g'),scatter(T.x(~b.Ti.ndx)/nRD,T.z(~b.Ti.ndx)/nRD,'*r')
% caxis([-2 2]),shading flat, colormap('gray'),title('HarmFEM')