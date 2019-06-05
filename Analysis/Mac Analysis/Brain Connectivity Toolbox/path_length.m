function x=path_length(net_mat)
G=graph(net_mat);
size = numnodes(G);
x = [];
for c = 1:size
      for dest=1:size
          if(c~=dest)
              [P,d] = shortestpath(G,c,dest);
               x = [x, d];
          end
      end
end 