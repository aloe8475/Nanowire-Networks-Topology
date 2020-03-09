function barr=NumbToBoolArray(numb)
cArr=dec2bin(numb,9);
cArr=fliplr(cArr);
barr=zeros(1,9);
for i=1:9
    barr(i)=str2double(cArr(i));
end

end