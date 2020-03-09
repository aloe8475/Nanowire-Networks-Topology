function strout=RemoveWhiteSpaces(string,delimiter)


idxcomma=[];
for i=2:length(string)
    if isequal(string(i),' ') && ~isequal(string(i-1),' ')
        idxcomma=[idxcomma i];
    end
end
string(idxcomma)=delimiter;
string(string==' ')=[];
strout=string;
end