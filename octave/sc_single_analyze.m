pkg load parallel

resp = [];
doma = [];
sc_prepare

ii = g.prop.nfi;

input = load(['R',num2str(converge),'_',num2str(ii),'_',num2str(g.prop.nfr),'_',num2str(g.prop.iff*g.model_scale*100),'_',num2str(g.model_scale*100),'.mat']);

Mcmb = input.Mcmb;
sc_integrate

r0 = diag(response.range.pid);
r1 = diag(response.range.prl);
r2 = diag(response.range.prr);

resp0 = zeros(1,ii);
doma0 = zeros(size(field.range.pid,2),ii);

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
