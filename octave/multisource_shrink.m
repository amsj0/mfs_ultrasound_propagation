function PF = multisource_shrink(varargin)

g = evalin('base','g');

kr_par = varargin{1}; % Wavenumber parameters
dr_par = varargin{2}; % Density parameters
pr_par = varargin{3}; % Sensitivity parameters
GY =  varargin{4}; % Number of wavelenght in y axis
GX =  varargin{5}; % Number of wavelenght in x axis
nRD =  varargin{6};% Number of wavelenght per radius (NOT USED SHOULD SCALE INTERFACE)
Neltoverlambda =  varargin{7}; % Elements per lambda

% WAVENUMBER VECTOR RANGE
k0 = 2*pi*linspace(nRD,nRD,g.prop.nfr).*(g.prop.sfr/g.prop.fr);
skr = linspace(kr_par(1),kr_par(3),kr_par(2)); 
kr = skr/1000;k0_length = length(k0);kr_length = length(kr);
% DENSITY VECTOR RANGE
d0 = g.prop.rj;
sdr = linspace(dr_par(1),dr_par(3),dr_par(2));
dr = sdr/1000;d0_length = length(d0);dr_length = length(dr);
% SENSITIVITY VECTOR RANGE
p0 = 1;
spr = linspace(pr_par(1),pr_par(3),pr_par(2));
pr = spr/100;p0_length = length(p0);pr_length = length(pr);
% WAVENUMBER VECTOR RANGE
sfr = g.prop.sfr/1000;

RD = nRD*g.prop.wav;
max_k = real(max(k0)*max(kr)); % Maximum wavenumber
max_k0 = real(max(k0)); % Maximum wavenumber
% Np = ceil(max_k*Neltoverlambda/(GY))*GY; % Number of points

% SET DIMENSIONS FOR RECTANGULAR GRID
lambda0 = 2*pi/max_k0;
% gap = 1/(GY*32);
mode = 1;

he_ini = g.interface_centre/g.prop.wav;
he_piston_centre = g.piston_centre;
% we_ini = (0+g.grid_ratio-0.001);
we_ini = -(g.centre_vector_x+10*GX);
%%
% CREATE SURROUNDING SURFACE
% [r,Th,area,norm,Np] = fn_discretize_geometry_3(GY,GX,nRD,Neltoverlambda);
[r,Th,area,norm,Np] = fn_discretize_geometry_plane(g.piston_distan/2,(1+1)*g.piston_distan,Neltoverlambda/100);

gap = 1/(Np);

r_e = cat(2,r{:});
n_e = cat(2,norm{:});
b.t_e = cat(2,Th{:});

normP = cellfun(@(x) x+pi,norm,'UniformOutput',false);

nO = cat(2,normP{:});
nI = cat(2,norm{:});
nR = nO;

S.a = cellfun(@(x) RD*x,area,'uni',0);
S.x = cellfun(@(x,y) RD*x.*cos(y),r,Th,'uni',0);
S.z = cellfun(@(x,y) RD*(x.*sin(y)+he_ini),r,Th,'uni',0);
S.y = cellfun(@(x) zeros(size(x)),S.z,'uni',0);

S.x = permute(cell2mat(S.x),[2 1]);
S.z = permute(cell2mat(S.z),[2 1]);

S.x = S.x(:);
S.z = S.z(:);

area_un = cat(2,S.a{:});
%%
% RADIAL SHIFT BASED ON CIRCULAR PACKING    
f0 = 1.0;
ff  = (f0-2*sin(pi/Np)/(1+sin(pi/Np)));
% RADIAL SHIFT BASED ON CIRCULAR PACKING    
gap = .0015;
% s = 200;
s = 3;      % Distance parameter from surface to point sources
% s = 1;
a = ceil(log((1-sqrt(3)*pi/Np)^s-gap*lambda0)/log((1-sqrt(3)*pi/Np)));
hap = ((1+sqrt(3)*pi/Np)^a-1);    
a = -ceil(log((1-sqrt(3)*pi/Np)^s-gap*lambda0)/log((1-sqrt(3)*pi/Np)));
ham = -(1-(1-sqrt(3)*pi/Np)^a);

f(1) = ff^(s);
f(2) = ff^(-s);

% RECTANGULAR SHIFT BASED ON RECTANGULAR PACKING
f(1) = f0-250/(Neltoverlambda);
f(2) = f0+250/(Neltoverlambda);

S.c = (S.x+1i*S.z)*f0;
S.co = S.c+RD*(f(2)-f0)*1.0*exp(1i*nI.');
S.ci = S.c-RD*(f0-f(1))*1.0*exp(1i*nI.');    

S.ct = S.c*(1-ham);

% eye_T = eye(2*length(area_un));
%%
% CREATE TRANSMITER TRANSDUCER SURFACE
[r,Th,area,norm,Np] = fn_discretize_geometry_plane(0,g.piston_vector_pitch,Neltoverlambda/(100));

T.a = cellfun(@(x) RD*x,area,'uni',0);
T.x = cellfun(@(x,y) x.*sin(y),r,Th,'uni',0);
T.z = cellfun(@(x,y) x.*cos(y),r,Th,'uni',0);
T.y = cellfun(@(x) zeros(size(x)),T.z,'uni',0);

T.a = permute(cell2mat(T.a),[2 1]);
T.x = permute(cell2mat(T.x),[2 1]);
T.y = permute(cell2mat(T.y),[2 1]);
T.z = permute(cell2mat(T.z),[2 1]);

T.a = T.a(:).';
T.x = T.x(:).';
T.y = T.y(:).';
T.z = T.z(:).'+he_piston_centre;

% SLICE EMISSOR SURFACE WITH SURROUNDING SURFACE
ndx0 = fn_enclosure_rectan(lambda0,T.x,T.z,[1,1],[he_ini,we_ini],0,mode);
ndx1 = fn_enclosure_rectan(lambda0,T.x,T.z,[-1,1],[he_ini,we_ini],0,mode);

T.x = RD*T.x;
T.z = RD*T.z;

b.Ti.ndx = ~ndx0;
b.To.ndx = ~ndx1;

% CREATE RECEPTOR TRANSDUCER SURFACE
[r,Th,area,norm,Np] = fn_discretize_geometry_plane(0,g.piston_vector_catch,Neltoverlambda/(100));

R.a = cellfun(@(x) RD*x,area,'uni',0);
R.x = cellfun(@(x,y) x.*sin(y),r,Th,'uni',0);
R.z = cellfun(@(x,y) x.*cos(y),r,Th,'uni',0);
R.y = cellfun(@(x) zeros(size(x)),R.z,'uni',0);

R.a = permute(cell2mat(R.a),[2 1]);
R.x = permute(cell2mat(R.x),[2 1]);
R.y = permute(cell2mat(R.y),[2 1]);
R.z = permute(cell2mat(R.z),[2 1]);

R.a = R.a(:).';
R.x = R.x(:).';
R.y = R.y(:).';
R.z = R.z(:).'+he_piston_centre;

% R = T;

% SLICE RECEPTOR SURFACE WITH SURROUNDING SURFACE
ndx0 = fn_enclosure_rectan(lambda0,R.x,R.z,[1,1],[he_ini,we_ini],0,mode);
ndx1 = fn_enclosure_rectan(lambda0,R.x,R.z,[-1,1],[he_ini,we_ini],0,mode);

R.x = RD*R.x;
R.z = RD*R.z;

R.x = R.x + RD*(g.piston_distan);

b.Ri.ndx = ~ndx0;
b.Ro.ndx = ~ndx1;

[b.Ti,b.To] = copy_filter(T,b.Ti,b.To);
[b.Ri,b.Ro] = copy_filter(R,b.Ri,b.Ro);

% b.Ri.x = b.Ri.x + RD*(g.centre_vector+GX);
% b.Ro.x = b.Ro.x + RD*(g.centre_vector+GX);

r_r(1,:) = g.pool_radius(1)*r_e;
r_r(2,:) = g.pool_radius(2)*r_e;
r_r(3,:) = g.pool_radius(3)*r_e;
r_r(4,:) = g.pool_radius(4)*r_e;

dl = g.grid_ratio;
grid_number_z = floor(GY/dl);
grid_number_x = floor(GX/dl);
% rmax = RD*g.pool_radius(1)*nRD*lambda0/cos(pi/GY);
%%
% CREATE REGION OF INTEREST BASED ON MODEL COMPARISON TAG
if (strcmp(g.extension,'eswsnan'))
    vec.x = g.centre_vector_x+(-(grid_number_x-.5):(grid_number_x-.5))*dl;
    vec.z = g.centre_vector_y+(-(grid_number_z-.5):(grid_number_z-.5))*dl;
    col = grid_number_x*2 + 1*0;
    row = grid_number_z*2 + 1*0;

    M.x = vec.x(ones(row,1),:);
    M.x = M.x(:)';
    M.z = vec.z(ones(col,1),:)';
    M.z = M.z(:)';
    M.y = zeros(1,row*col);
    M.a = zeros(1,row*col);
    
    b.A.x = RD*M.x;
    b.A.z = RD*M.z;
    b.A.y = RD*M.y;

    b.A.x = reshape(b.A.x,[row,col]);
    b.A.y = reshape(b.A.y,[row,col]);
    b.A.z = reshape(b.A.z,[row,col]);
    
    b.A.p0 = zeros(1,row*col);
    b.A.a = zeros(1,row*col);
    
elseif (strcmp(g.extension,'eswsfem'))    
    data = load('pressure-250kHz-y_interf0.035.txt');

    b.A.x = data(1:(g.scale):end, 1:4*(g.scale):end)*nRD;
    [row,col] = size(b.A.x);
    
    b.A.z = data(1:(g.scale):end, 2:4*(g.scale):end)*nRD;
    
    b.A.y = zeros(row,col);
    
    M.x = b.A.x(:)'/RD;
    M.z = b.A.z(:)'/RD;
    M.y = b.A.y(:)'/RD;
    
%     b.A.p = zeros(1,row*col);
%     b.A.p0 = - data(1:(g.scale):end, 3:4*(g.scale):end) + 1j * data(1:(g.scale):end, 4:4*(g.scale):end);
      b.A.p0 = 1j * data(1:(g.scale):end, 3:4*(g.scale):end) + data(1:(g.scale):end, 4:4*(g.scale):end);
%     bAp0 = -1j*conj(b.A.p0);
end
    
b.A.p = zeros(row,col);

ndx1 = fn_enclosure_rectan(lambda0,M.x,M.z,[-1,1],[-1,we_ini],gap,mode);
ndx2 = fn_enclosure_rectan(lambda0,M.x,M.z,[1,1],[1,we_ini],gap,mode);
ndx3 = fn_enclosure_rectan(lambda0,M.x,M.z,[-1,1],[he_ini,we_ini],gap,mode);
ndx4 = fn_enclosure_rectan(lambda0,M.x,M.z,[1,1],[he_ini,we_ini],gap,mode);

b.D.ndx = (ndx1&ndx2)&(ndx3&ndx4);             % Entire domain, for pressure field

[b.D] = copy_filter(M,b.D);

% ndx0 = fn_enclosure_polar(lambda0,b.D.x,b.D.z,[0,0],1,r_r(4,:),n_e,b.t_e,1*gap);
ndx0 = fn_enclosure_rectan(lambda0,b.D.x,b.D.z,[1,1],[he_ini,-we_ini],gap,mode);
% ndx0 = fn_enclosure_rectan(lambda0,M.x,M.z,[1,1],[0,-we_ini],gap,mode);

len_I = sum(ndx0);
len_O = sum(~ndx0);
% len_I = sum(ndx0a);
len_C = len_I + len_O;
clear M

b.D.ndx0 = ndx0;

% len_C = ceil(length(b.C.x));
% len_S = ceil(length(S.x));
len_T = ceil(length(b.D.x));

b.D.c = RD*(b.D.x + 1i*b.D.z).';
b.D.ci = b.D.c(ndx0);
b.D.co = b.D.c(~ndx0);    

b.Ti.c =  (b.Ti.x + 1i*b.Ti.z);
b.Ti.p =  zeros(1,length(b.Ti.x));

b.To.c =  (b.To.x + 1i*b.To.z);
b.To.p =  zeros(1,length(b.To.x));

b.Ri.c =  (b.Ri.x + 1i*b.Ri.z).';
b.Ri.p =  zeros(1,length(b.Ti.x)).';

b.Ro.c =  (b.Ro.x + 1i*b.Ro.z).';
b.Ro.p =  zeros(1,length(b.To.x)).';

% PF = zeros(k0_length,kr_length,pr_length,row,col);
% VF = zeros(k0_length,kr_length,pr_length,row,col);

kci = [];
if size(b.Ti.c)
    kci = bsxfun(@minus,b.D.ci,b.Ti.c);
end
kco = [];
if size(b.To.c)
    kco = bsxfun(@minus,b.D.co,b.To.c);
end
mci = [];
if size(b.Ri.c)
    mci = bsxfun(@minus,b.Ri.c,b.Ti.c);
end
mco = [];
if size(b.Ro.c)
    mco = bsxfun(@minus,b.Ro.c,b.To.c);
end
% % % % % % % % % % % % % TEMPORARLY % % % % % % % % % % %     
%     PF = zero s(len_C,length(tiers_v),kr_length,dr_length);
% % % % % % % % % % % % % TEMPORARLY % % % % % % % % % % %  
d_cur = d0;
% sc_prepare
save([g.path_to_input,'P',g.convergemod,'_',num2str(g.prop.nfr),'_',num2str(g.prop.iff*g.model_scale*100),'_',num2str(g.model_scale*100),'.mat'],'T','R','S','Neltoverlambda','nRD','nR','area_un','b','g')
for ii = g.prop.nfi:g.prop.nff
    k_cur = k0(ii)*lambda0/(RD);
%     k_cur = k0(ii)*lambda0/(g.prop.wav);

    sf_cur = sfr(ii);
    disp('Spectrum Ratio');
    disp(sf_cur)
% COMPUTE FIELD ON BOUNDARY - TEST DOMAIN / BOUNDARY (1 ORDER TEST)
  
    TIh = 0;
    if size(b.Ti.c)
%         [TIj,TIy] = fn_compute_field_boundary_0(S,b.Ti,nR,k_cur,k_cur);
        [~,~,TIh] = fn_compute_field_boundary3(S,b.Ti,nR,area_un,k_cur,k_cur,g.golden_ratio);        
%         [~,~,TIh] = fn_compute_field_boundary3(S,b.Ti,nR,1,k_cur,k_cur,g.golden_ratio);        
    % FIELD ON BOUNDARY - TEST DOMAIN / BOUNDARY
%         TIh = (TIj + TIy);
    end
% COMPUTE MONOPOLE FIELD INSIDE - VIRTUAL DOMAIN / BOUNDARY
    [p2h,v2h,~] = fn_compute_field_inside_m2(S,nR,k_cur,k_cur,1); 

% COMPUTE MONOPOLE FIELD INSIDE - VIRTUAL DOMAIN / DOMAIN
%  [p2i,phI_vd] = fn_compute_field_inside_domain(S,k_cur);

% COMPUTE REFERENCE FIELD INSIDE - INNER DOMAIN (0 ORDER TEST)
    p0kI = fn_compute_reference_0_2(kci,k_cur,area_un);
%     p0kI = fn_compute_reference_0_2(kci,k_cur,1);

% COMPUTE REFERENCE FIELD INSIDE - RECEIVER SURFACE (0 ORDER TEST)    
    p0mI = fn_compute_reference_0_2(mci,k_cur,area_un);

% COMPUTE PROPAGATOR INSIDE - BOUNDARY 2 / INNER DOMAIN
    [~,~,Thi] = fn_compute_propagator_inside3(S,nR,area_un,k_cur,k_cur,len_I,1);         
        
% PROPAGATOR PROPAGATOR INSIDE - BOUNDARY / INNER DOMAIN
%     Thi = TUji + TUyi;  
%     Thi = TUhi;
    
% COMPUTE PROPAGATOR INSIDE - BOUNDARY 2 / INNER DOMAIN
    [~,~,Tri] = fn_compute_propagator_receiverinside_2(S,nR,area_un,k_cur,k_cur,len_I,1);         
        
% PROPAGATOR PROPAGATOR INSIDE - BOUNDARY / INNER DOMAIN
%     Tri = TRji + TRyi;      
    
    k_out = k_cur;
    d_out = d_cur;

    for jj =1:kr_length   

        for pp =1:dr_length   
            d_cur = 1/dr(pp);
            sd_cur = sdr(pp);
    
% % % % % % % % % % % % % TEMPORARLY % % % % % % % % % % %           
%             k_cur = k_out/kr(jj)*(1+1j*g.prop.att)/abs(1+1j*g.prop.att);
            k_curi = k_out/kr(jj);
%             k_cur = k_out/kr(jj)*(1-1j*g.prop.att)/abs(1-1j*g.prop.att);
            k_cur = k_out/kr(jj)*(1-1j*g.prop.att)/abs(1-0*1j*g.prop.att);
%             k_cur = k_out/(kr(jj)*(1-1i*g.prop.att));            
            disp('Wavenumber Ratio');
            disp(kr(jj)/(1-1j*g.prop.att)/abs(1-0*1j*g.prop.att))
%             disp(kr(jj)*(1-1i*g.prop.att))
            disp('Density Ratio');
            disp(dr(pp))
            sk_cur = skr(jj);
        %    sk_cur = skr(ceil(end/2));
        % COMPUTE FIELD ON BOUNDARY - TEST DOMAIN / BOUNDARY (1 ORDER TEST)
        %       [TIj,TIy] = fn_compute_field_boundary(S,b.C,nO,k_cur,k_cur);
            TOh = 0;
            if size(b.To.c)
%                 [TOj,TOy] = fn_compute_field_boundary_0(S,b.To,nR,k_cur,k_cur);
                [~,~,TOh] = fn_compute_field_boundary3(S,b.To,nR,area_un,k_cur,k_cur,g.golden_ratio);        
%                 [~,~,TOh] = fn_compute_field_boundary3(S,b.To,nR,1,k_cur,k_cur,g.golden_ratio);        
            % FIELD ON BOUNDARY - TEST DOMAIN / BOUNDARY
%                 TOh = -(TOj + TOy);
                TOh = -TOh;
            end
            TOhx = 0;
            if numel(TOh)>2
                TOhx = [TOh(1:(end/2),:);-TOh((end/2+1):end,:)];
            end
            TIhx = 0;
            if numel(TIh)>2
                TIhx = [TIh(1:(end/2),:);-TIh((end/2+1):end,:)];
            end
            
        % COMPUTE MONOPOLE FIELD OUTSIDE - VIRTUAL DOMAIN / BOUNDARY            
%             [p1h,v1h,~] = fn_compute_field_outside_m2(S,nR,k_curi,k_out,d_cur);
            [p1h,v1h,~] = fn_compute_field_outside_m2(S,nR,k_curi,k_out,d_cur);
        % COMPUTE MONOPOLE FIELD OUTSIDE - VIRTUAL DOMAIN / DOMAIN
        %   [p2o,phO_vd] = fn_compute_field_outside_domain(S,k_cur);   
        % COMPUTE REFERENCE FIELD OUTSIDE - OUTER DOMAIN (0 ORDER TEST)
             p0kO = fn_compute_reference_0_2(kco,k_cur,area_un);
%             p0kO = fn_compute_reference_0_2(kco,k_cur,1);
        % COMPUTE REFERENCE FIELD OUTSIDE - RECEIVER SURFACE (0 ORDER TEST)             
             p0mO = fn_compute_reference_0_2(mco,k_cur,area_un);
        % COMPUTE PROPAGATOR OUTSIDE - BOUNDARY / INNER DOMAIN
             [~,~,Tho] = fn_compute_propagator_outside3(S,nR,area_un,k_cur,k_out,len_O,d_cur);                
        % PROPAGATOR PROPAGATOR OUTSIDE - BOUNDARY / INNER DOMAIN
%              Tho = TUjo + TUyo;                
        % COMPUTE PROPAGATOR OUTSIDE - BOUNDARY / INNER DOMAIN
             [~,~,Tro] = fn_compute_propagator_receiveroutside_2(S,nR,area_un,k_cur,k_out,len_O,d_cur);                
        % PROPAGATOR PROPAGATOR OUTSIDE - BOUNDARY / INNER DOMAIN
%              Tro = TRjo + TRyo;                 
        %    p2 = [p2o;p2i];
        
        % COMPUTE PROPAGATOR INCIDENT TO REFRACTED - BOUNDARY / BOUNDARY
             [TFI] = fn_propagator_inc_ref(p1h,v1h,-p2h,-v2h);
             [TFO] = fn_propagator_inc_ref(p2h,v2h,-p1h,-v1h);
%              [TFI] = fn_propagator_inc_ref(p1h,v1h,p2h,v2h);
             
        %    TFO = (eye_T - TFI);
        %    TF = fn_propagator_src_beh(p1h,v1h,-p2h,-v2h);
        %    TF = fn_propagator_src_beh(p1h,v1h,p2h,v2h);

        % COMPUTE PROPAGATOR REFRACTED TO INVERTED INCIDENT - BOUNDARY / BOUNDARY
        %    iTFa = fn_invert_propagator(TF);
        %    idTha_Th = iTFa*TIh;  % APPLIED MODEL
        
        % STORE FIELD PARAMETERS
%              PF(ii,jj,pp,:,:) = b.A.p;
        sc_propagate

%         [TFO] = fn_propagator_inc_ref(p2h,-v2h,-p1h,+v1h);
%         sc_propagate2
%        sc_prepare
%         sc_integrate
%        resp_rang(ii) = response_range;
     %     Mcmb_l{ii} = Mcmb;
%         h5write(['RR_',num2str(ii),'_',num2str(g.prop.nfr),'_',num2str(g.prop.iff*g.model_scale*100),'_',num2str(g.model_scale*100),'.mat'],'/Mcmb',Mcmb)
         save([g.path_to_input,'R',g.convergemod,'_',num2str(ii),'_',num2str(g.prop.nfr),'_',num2str(g.prop.iff*g.model_scale*100),'_',num2str(g.model_scale*100),'_',num2str(sd_cur),'_',num2str(sk_cur),'.mat'],'Mcmb')

        end
          %    
    end
end

% save('PF.mat','PF')
% save(['PP_',g.prop.nfr,'_',g.model_scale,'.mat'],'T','R','S','Neltoverlambda','nRD','b','g')
% save(['D:\MATLAB\menisco\RR_',num2str(g.prop.nfr),'_',num2str(g.prop.iff*g.model_scale*100),'_',num2str(g.model_scale*100),'.mat'],'T','R','S','Neltoverlambda','nRD','b','g')
% sc_surf_def
% sc_surf_run