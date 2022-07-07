function ndx = fn_integrate_indexing(points,varargin)

size_I = varargin{1};
size_O = varargin{2};

% MATRIX MAGIC
matrix_ndx = @(v) cellfun(@(x) circshift(v,x,2),num2cell(0:(points-1)),'UniformOutput',0);

vec = [num2cell(1:size_I,1),cell(1,size_O)];
mat = matrix_ndx(vec);
ndx.catI = cat(1,mat{:});

vec = [cell(1,size_I),num2cell(1:size_O,1)];
mat = matrix_ndx(vec);
ndx.catO = cat(1,mat{:});