function [textdata,Vbias]=prepareMetaData(Names,Values)
textdata=struct();
for i=1:length(Names)
    %construc valid field names
    Nam=replace_invalid(Names{i});

    if isempty(Nam)
        continue;
    end
    textdata.(Nam)=Values{i};
end
Vbias=[];


%FOR REFERENCE
% textdata=struct();
% textdata.Filename=Out.propValues{1,1}{1,1};
% textdata.DateHour=Out.propValues{1,3}{1,1};
% textdata.MeasType='Oscillo-TDMS';
% textdata.InputCh='100000000';
% textdata.OutputCh='100000000'; %it does not count for two channel stage yet'
% textdata.Vbias=Out.propValues{1,1}{1,2};
% bias=str2double(textdata.Vbias);
% textdata.Offset='0';
% textdata.SampleRate=Out.propValues{1,1}{1,5};
% textdata.Buffer=Out.propValues{1,1}{1,6};
% textdata.OffArray='NaN';
% textdata.VArray='NaN';
% textdata.Extra=strcat(strjoin(inchan),' , ',strjoin(outchan));


end