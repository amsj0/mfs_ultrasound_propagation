function [rC,thC,areaC,normC,M] = fn_discretize_geometry_plane(centre,RN,Nelt,mode)

N = 1;
% x = -nRD/2:1/Nelt:nRD/2;
% y = zeros(len,1);
% r = sqrt(x^2+y^2);
% th = atan2(y,x);
vec = 0;
% number = RN*Nelt;

pts_number = floor(RN / 2 * Nelt) * 2;
RN = (pts_number - 1) / Nelt;

if pts_number > 1
    vec = linspace(-RN/2,RN/2,pts_number);
end
r = centre + vec;
% r = centre + (-RN/2:1/Nelt:RN/2);
th = zeros(1,length(r));

curva = inf*ones(1,length(r));
area = 1/Nelt*ones(1,length(r));
norm = -pi/2*ones(1,length(r));

assignin('base','cm',min(curva));

M = length(r)*N;

rC = cell(N,1);
areaC = cell(N,1);
normC = cell(N,1);
thC = cell(N,1);

for k=1:N
    rC{k} =  r;
    areaC{k} = area;
    normC{k} = norm + (k-1)*(2*pi/N);
%     thC{k} = fn_anglebounds(th + (k-1)*(2*pi/N));
    thC{k} = th + (k-1)*(2*pi/N);
end

