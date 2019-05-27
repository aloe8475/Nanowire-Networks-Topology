function [mabObj]=hashtable(wordsRemove)

for i=1:numel(wordsRemove)
    valueSet(i)=i;
end 

mapObj = containers.Map(wordsRemove,valueSet);
end 