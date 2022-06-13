function rd = fn_analyze(path,datafile,configfile,analysisfile)
    
    load([path,datafile])
    load([path,configfile,'.h5'])
	load([path,analysisfile])
	%   DEFINES PISTON INDEXES
	ppt_per_surface = 1+floor(g.piston_radius*(Neltoverlambda/100));
	%ppt_per_surface = 1;

	ndx.T = fn_integrate_indexing(ppt_per_surface,length(b.Ti.x),length(b.To.x));
	ndx.R = fn_integrate_indexing(ppt_per_surface,length(b.Ri.x),length(b.Ro.x));

	response.pitch.size = (length(T.a)-ppt_per_surface+1);
	response.catch.size = (length(R.a)-ppt_per_surface+1);

	field.range.ptt = zeros(response.pitch.size,length(b.D.c));

	strength.range.ptt = cell(response.pitch.size,response.catch.size);

	response.range.pid = zeros(response.pitch.size,response.catch.size);
	response.range.prr = zeros(response.pitch.size,response.catch.size);
	response.range.prl = zeros(response.pitch.size,response.catch.size);
	response.range.ptt = zeros(response.pitch.size,response.catch.size);

	response.range.ptx = zeros(response.pitch.size,1);
	response.range.prx = zeros(response.catch.size,1);

	response.pitch.A = zeros(1,response.pitch.size);
	response.catch.A = zeros(1,response.pitch.size);

	response.pitch.I = T.z(~b.Ti.ndx)/nRD;
	response.pitch.O = T.z(~b.To.ndx)/nRD;

	response.catch.I = R.z(~b.Ri.ndx)/nRD;
	response.catch.O = R.z(~b.Ro.ndx)/nRD;

	response.ndx.RO = find(~b.Ro.ndx);
	response.ndx.RI = find(~b.Ri.ndx);

	response.ndx.TO = find(~b.To.ndx);
	response.ndx.TI = find(~b.Ti.ndx);

	%%% TODO FIX ME
	ndx.tst = reshape(1:(response.pitch.size*response.catch.size),response.pitch.size,response.catch.size);
	[n1,n2] = ndgrid(1:response.pitch.size,1:response.catch.size);n0 = (mod(n2+n1+1,response.pitch.size)+1)';
	ndx.mat = reshape(ndx.tst((n0(:)-0)+(n1(:)-1)*response.catch.size),response.pitch.size,response.catch.size)';

	% pause
	for tt = 1:response.pitch.size
	        ndxI = [ndx.T.catI{:,ppt_per_surface-1+tt}];
	        ndxO = [ndx.T.catO{:,ppt_per_surface-1+tt}];
	        response.pitch.A(tt) = mean([response.pitch.I(ndxI),response.pitch.O(ndxO)]);
	        ndx.T.i{tt} = ndxI;
	        ndx.T.o{tt} = ndxO;
	end

	for rr = 1:response.catch.size
	        ndxI = [ndx.R.catI{:,ppt_per_surface-1+rr}];
	        ndxO = [ndx.R.catO{:,ppt_per_surface-1+rr}];
	        response.catch.A(rr) = mean([response.catch.I(ndxI),response.catch.O(ndxO)]);
	        ndx.R.i{rr} = ndxI;
	        ndx.R.o{rr} = ndxO;
	end


	for tt = 1:response.pitch.size

	    %   EXTRACT TRANSMITTER PISTON INDEXES
	    ndx.TI = ndx.T.i{tt};
	    ndx.TO = ndx.T.o{tt};

	    resTO = response.ndx.TO(ndx.TO);
	    resTI = response.ndx.TI(ndx.TI);

	    %   FILTER FIELD PARAMETERS
	    b.D.prr(b.D.ndx0,:)  = energy_ratio*sum(Mcmb{1,1}(:,ndx.TO),2);  % APPLIED MODEL
	    b.D.prr(~b.D.ndx0,:) = sum(Mcmb{2,1}(:,ndx.TI),2);  % APPLIED MODEL 
	    b.D.prl(b.D.ndx0,:)  = sum(Mcmb{3,1}(:,ndx.TI),2);  % APPLIED MODEL 
	    b.D.prl(~b.D.ndx0,:) = energy_ratio*sum(Mcmb{4,1}(:,ndx.TO),2);  % APPLIED MODEL 
	    b.D.pid(b.D.ndx0,:)  = sum(Mcmb{5,1}(:,ndx.TI),2);  % APPLIED MODEL 
	    b.D.pid(~b.D.ndx0,:) = energy_ratio*sum(Mcmb{6,1}(:,ndx.TO),2);  % APPLIED MODEL 
	    b.D.ptt = b.D.pid + b.D.prl + b.D.prr;

	    field.range.prr(tt,:) = b.D.prr;
	    field.range.prl(tt,:) = b.D.prl;
	    field.range.pid(tt,:) = b.D.pid;
	    %field.range.ptt(tt,:) = b.D.ptt;

	%    response.pitch.A(tt) = mean([response.pitch.I(ndx.TI),response.pitch.O(ndx.TO)]);

	    b.T.prr{1} = energy_ratio*sum(Mcmb{1,2}(:,ndx.TO),2);  % APPLIED MODEL
	    b.T.prr{2} = sum(Mcmb{2,2}(:,ndx.TI),2);  % APPLIED MODEL
	    b.T.prl{1} = sum(Mcmb{3,2}(:,ndx.TI),2);  % APPLIED MODEL
	    b.T.prl{2} = energy_ratio*sum(Mcmb{4,2}(:,ndx.TO),2);  % APPLIED MODEL
	    b.T.pid{1} = sum(Mcmb{5,2}(:,ndx.TI),2);  % APPLIED MODEL
	    b.T.pid{2} = energy_ratio*sum(Mcmb{6,2}(:,ndx.TO),2);  % APPLIED MODEL
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
	        %   UPDATE FIGURES


	        %   ANALYSE RESPONSE RANGES
	        response.range.pid(tt,rr) = (sum(b.T.pid{1}(ndx.RI)) + sum(b.T.pid{2}(ndx.RO)));
	        response.range.prr(tt,rr) = (sum(b.T.prr{1}(ndx.RI)) + sum(b.T.prr{2}(ndx.RO)));
	        response.range.prl(tt,rr) = (sum(b.T.prl{1}(ndx.RI)) + sum(b.T.prl{2}(ndx.RO)));

	    end
	%     if(rr==floor(response_size/2))
	    response.range.ptx(tt) = (sum(response.range.prx));
	end

	response.range.ptt = response.range.prl+response.range.prr+response.range.pid;

    %rr0 = diag(response.range.pid);
	%rr1 = diag(response.range.prl);
	%rr2 = diag(response.range.prr);

	%rr0 = response.range.pid(:);
	%rr1 = response.range.prl(:);
	%rr2 = response.range.prr(:);	
	
	%rr0 = response.range.pid;
	%rr1 = response.range.prl;
	%rr2 = response.range.prr;	

	%rd.rr = cat(2,rr0,rr1,rr2);
	rd.rr = cat(3,response.range.pid,response.range.prl,response.range.prr);
	rd.dd = cat(3,field.range.prr,field.range.prl,field.range.pid);
end
