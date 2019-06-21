function x=path_length(net_mat)
G=graph(net_mat);
size = numnodes(G);
x = [];
    time1=tic;
for c = 1:size
      for dest=1:size
          if(c~=dest)
              [P,d] = shortestpath(G,c,dest);
              if mod(c,100)==0 && dest==1
                fprintf(['\n Shortest Path Completed:' num2str(c) '\n']);
                toc(time1)
              end 
               x = [x, d];
          end
      end
end 