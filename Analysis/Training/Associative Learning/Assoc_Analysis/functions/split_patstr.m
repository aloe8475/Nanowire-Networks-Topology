function  [digin,digout,target]=split_patstr(StrPat)
pat=StrPat(4:end);
i=1;
BoolText=cell(1,3);remain=pat;
while ~isempty(remain)
    [pat,remain]=strtok(remain,'/');
    numb=pat;remain2=numb;
    while ~isempty(remain2)
        [numb,remain2]=strtok(remain2,'.');
        if isempty(BoolText{i})
            BoolText{i}=numb;
        else
            BoolText{i}=horzcat(BoolText{i},numb);
        end
    end
    i=i+1;
end
digin=BoolText{1};
digout=BoolText{2};
target=BoolText{3};
end