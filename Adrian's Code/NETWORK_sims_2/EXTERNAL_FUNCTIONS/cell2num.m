function varargout=cell2num(varargin)
%converts cell array of double array to a single array of double concatenating
%cells
if ~iscell(varargin{1})
    return
end
LArgs=length(varargin);
LArr=length(varargin{1});
Vec=cell(1,LArgs);
for i=1:LArgs
    for j=1:LArr
        sz=size(varargin{i}{j});
        if ~isequal(sz(1),1)
            join=varargin{i}{j}';
        else
            join=varargin{i}{j};
        end
        Vec{i}=[Vec{i}  join];
    end
end
varargout=Vec;
end
