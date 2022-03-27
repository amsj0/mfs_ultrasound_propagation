function varargout = fn_enclosure_rectan(disc,xM,zM,sign,t,gap,~)
 % varargout = fn_enclosure_polar(xM,zM,sign,t)
% th = atan2(zM,xM);
th{1} = zM;
th{2} = xM;
% if (mode)
%  th(th<0) = th(th<0)+2*pi;
% end

loutput = length(sign);

varargout = cell(loutput);

for i = 1:loutput

% if(sign(i)>0)
%     ndx = (th(i)<=(t(1)+disc*gap));
% else
%     ndx = (-(t(i)-disc*gap)>=-th(i));
% end

ndx = (sign(i)*th{i}<=sign(i)*(t(i)+sign(i)*disc*gap));

% if(sign(i)>0)
%     ndx = (th(i)<=(t(1)+sign(i)*disc*gap));
% else
%     ndx = (sign(i)*th(i)<=sign(i)*(t(i)+sign(i)*disc*gap));
% end

varargout{i} = ndx;

end