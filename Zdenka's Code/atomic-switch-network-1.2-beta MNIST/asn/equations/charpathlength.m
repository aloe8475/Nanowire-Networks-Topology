%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Finds average path length for a given set of nodes, samples taken to
%approximate average for the network

%Author: Miro Astore
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%number of nodes to be tested
function cpl = charPathLength(adjmatrix, trials, number_of_wires)

lengths = zeros(trials, trials);
i = randi(number_of_wires, 1, trials);
j = randi(number_of_wires, 1, trials);

%using already written findPath, calculate shortest distance between each
%trial node
for k=1:trials
%     disp(k)
    for l = 1:trials
        lengths(l,k) = length(findPath(adjmatrix,i(k),j(l)));
    end
end

%calculate average for nodes tested
cpl = sum(sum(lengths))/numel(lengths);
return;