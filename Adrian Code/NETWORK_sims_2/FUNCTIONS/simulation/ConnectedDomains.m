function domains=ConnectedDomains(GraphObject)

%% get domains from matlab graph function conncomp
bins=conncomp(GraphObject);


%% calculate entries for every domain and retrieve indices
% an unconnected nanowire is not considered a domain (its sorted out in
% line 11)
[x,~]=histc(bins,unique(bins));
idx=find(x~=1);
domains=cell(length(idx),1);

%%fill the cell with the encountered domains
for i=1:length(idx)
    domains{i}=find(idx(i)==bins);
end
end