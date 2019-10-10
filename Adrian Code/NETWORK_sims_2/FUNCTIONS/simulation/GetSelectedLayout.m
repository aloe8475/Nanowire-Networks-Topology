function SelLayout=GetSelectedLayout(Network,DIDx)
%% dimesionality reduction of netowrk to get only connected nanowires
% for simulations and plot
SelLayout=struct();
AdjMat=adjacency(Network.Graph);
SelLayout.AdjMat=AdjMat(DIDx,DIDx);
SelLayout.SelGraph=graph(SelLayout.AdjMat);
SelLayout.CX=sparse(Network.CrossMat.CX(DIDx,DIDx));
SelLayout.CY=sparse(Network.CrossMat.CY(DIDx,DIDx));
SelLayout.X1=sparse(diag(Network.LayOut.x1(DIDx)));
SelLayout.X2=sparse(diag(Network.LayOut.x2(DIDx)));
SelLayout.Y1=sparse(diag(Network.LayOut.y1(DIDx)));
SelLayout.Y2=sparse(diag(Network.LayOut.y2(DIDx)));

end