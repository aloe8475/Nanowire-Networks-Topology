function numb=fuzzy_bool(bool_arr,dir)

if isequal(dir,'labview')
    bool_arr=fliplr(bool_arr);
end

lo=length(bool_arr);
numb=0;
for i=1:lo
    
    numb=numb+2^(lo-i)*(bool_arr(i));
    
end


end