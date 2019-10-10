function [D, E]=efficiency_bin(A,local)
%EFFICIENCY_BIN     Global efficiency, local efficiency.
%
%   Eglob = efficiency_bin(A);
%   Eloc = efficiency_bin(A,1);
%
%   The global efficiency is the average of inverse shortest path length,
%   and is inversely related to the characteristic path length.
%
%   The local efficiency is the global efficiency computed on the
%   neighborhood of the node, and is related to the clustering coefficient.
%
%   Inputs:     A,              binary undirected or directed connection matrix
%               local,          optional argument
%                                   local=0 computes global efficiency (default)
%                                   local=1 computes local efficiency
%
%   Output:     Eglob,          global efficiency (scalar)
%               Eloc,           local efficiency (vector)
%
%
%   Algorithm: algebraic path count
%
%   Reference: Latora and Marchiori (2001) Phys Rev Lett 87:198701.
%              Fagiolo (2007) Phys Rev E 76:026107.
%              Rubinov M, Sporns O (2010) NeuroImage 52:1059-69
%
%
%   Mika Rubinov, U Cambridge
%   Jonathan Clayden, UCL
%   2008-2013

% Modification history:
% 2008: Original (MR)
% 2013: Bug fix, enforce zero distance for self-connections (JC)
% 2013: Local efficiency generalized to directed networks
n=length(A);                                %number of nodes
A(1:n+1:end)=0;                             %clear diagonal
A=double(A~=0);                             %enforce double precision

if exist('local','var') && local            %local efficiency
    fprintf('Local Efficiency Started \n'); 
    E=zeros(n,1);   
    for u=1:n
%         fprintf([num2str(u) '\n']);
%         toc
%         fprintf('\n');
        V=find(A(u,:)|A(:,u).');            %neighbors
% %         fprintf('Finished Find \n');        
        sa=A(u,V)+A(V,u).';                 %symmetrized adjacency vector
%         fprintf('Finished Symmetrized Adj Vector \n');        
        [e, D]=distance_inv(A(V,V));             %inverse distance matrix
%         fprintf('Finished Inverse Distance Matrix \n');        
        se=e+e.';                           %symmetrized inverse distance matrix
        numer=sum(sum((sa.'*sa).*se))/2;    %numerator
        if numer~=0
            denom=sum(sa).^2 - sum(sa.^2);  %denominator
            E(u)=numer/denom;               %local efficiency
        end
        tic
    end
else                                        %global efficiency
    fprintf('Global Efficiency Started \n');
    [e, D]=distance_inv(A);
    E=sum(e(:))./(n^2-n); 
end


function [D1,D]=distance_inv(A_)
temp=graph(A_);
D=distances(temp);
n_=length(A_);
D(~D | eye(n_))=inf;                        %assign inf to disconnected nodes and to diagonal
D1=1./D;                                     %invert distance
