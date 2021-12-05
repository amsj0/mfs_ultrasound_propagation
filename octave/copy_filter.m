function varargout = copy_filter(varargin)

B = varargin{1};

varargin(1) = [];
varargout = cell(1,nargin-1);

for ii = 1:(nargin-1)
    C = varargin{ii};
    varargout{ii}.x = B.x(~C.ndx);
    varargout{ii}.y = B.y(~C.ndx);
    varargout{ii}.z = B.z(~C.ndx);
    varargout{ii}.ndx = C.ndx;
end