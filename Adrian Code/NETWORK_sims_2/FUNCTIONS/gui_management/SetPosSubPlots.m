function [OuterPosArr,dim]=SetPosSubPlots(non)
%% subplot function in matlab leaves a lot of blank space between plots
% this is a manual fix to reduce that space, but nonethelss works. 
% for more complex subplot settings, a more general function should be
% implemented
switch non
    case 0
        %2x2
        dim=[2,2];
        OuterPosArr(1,:)=[0 0.50 0.47 0.47];
        OuterPosArr(2,:)=[0.50 0.50 0.47 0.47];
        OuterPosArr(3,:)=[0 0 0.47 0.47];
        OuterPosArr(4,:)=[0.50 0 0.47 0.47];
    case 1
        %3x1
        dim=[3,1];
        OuterPosArr(1,:)=[0 0.66 1 0.33];
        OuterPosArr(2,:)=[0 0.33 1 0.33];
        OuterPosArr(3,:)=[0 0 1 0.33];
    case 2
        %2x1
        dim=[2,1];
        OuterPosArr(1,:)=[0 0.5 1 0.5];
        OuterPosArr(2,:)=[0 0 1 0.5];
    case 3
        %1x1
        dim=[1,1];
        OuterPosArr(1,:)=[0 0 1 1];
        
    otherwise
        OuterPosArray=[];
end

       



end