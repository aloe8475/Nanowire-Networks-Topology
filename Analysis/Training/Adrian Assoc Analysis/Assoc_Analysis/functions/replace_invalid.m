function corr_field=replace_invalid(fname)

    corr_field=strrep(fname,' ','_');
    corr_field=strrep(corr_field,'?','_');
    corr_field=strrep(corr_field,'(','_');
    corr_field=strrep(corr_field,')','_');
    corr_field=strrep(corr_field,'[','_');
    corr_field=strrep(corr_field,']','_');
    corr_field=strrep(corr_field,':','_');
    corr_field=strrep(corr_field,'\n','_');
    corr_field=strrep(corr_field,'.','_');
    corr_field=strrep(corr_field,'/','_');
    corr_field=strrep(corr_field,'%','_');
    corr_field=strrep(corr_field,'==','_');
    
    corr_field=regexprep(corr_field,'[\n\r]+','');