function [T,II,cx,cy,yxscale]=ExtractCurrentDataFromSim(ElSel,SimData,type,timeInd)
T=SimData.time;
T=repmat(T,[1 length(ElSel)]);
cx=T(timeInd,:);
II=zeros(length(T),length(ElSel));
yxscale.x='linear';yxscale.y='linear';
for i=1:length(ElSel)
    II(:,i)=SimData.(ElSel{i});
    vel=SimData.(strrep(ElSel{i},'I','V'));
    switch type
        case 'CurrentLog'
            yxscale.y='log';
        case 'Voltage'
            II(:,i)=vel;
        case 'Conductance'
            cond=II(:,i)./vel;
            cond(~isfinite(cond))=NaN;
            II(:,i)=cond;
        case 'ConductanceLog'
            cond=II(:,i)./vel;
            cond(~isfinite(cond))=NaN;
            II(:,i)=cond;
            yxscale.y='log';
        case 'Resistance'
            res=vel./II(:,i);
            res(~isfinite(res))=NaN;
            II(:,i)=res;
        case 'ResistanceLog'
            res=vel./II(:,i);
            res(~isfinite(res))=NaN;
            II(:,i)=res;
            yxscale.y='log';
        case 'AvJuncWidth'
            % not finished
%             LT=height(SimData);
%             mW=zeros(LT,1);
% %             for k=1:LT
% %                 W=SimData.Wmat{k};W=triu(W);
% %                 Wlist=W(W~=0);
% %                 mW(k)=mean(Wlist);
% % 
% %             end
%            
%             t=SimData.time;
%             fast_multiplot(t,vel,II);
% %             cond=II(:,i)./vel;
% %             plot(mW,cond);
%             %II(:,i)=mW;
    end   
end
cy=II(timeInd,:);
end