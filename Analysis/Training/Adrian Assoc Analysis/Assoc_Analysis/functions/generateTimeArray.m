function [time,data,LChans,ChannelNames]=generateTimeArray(data,textdata,ChannelNames)



LChans=length(ChannelNames);
if any(contains(ChannelNames,'time')) || any(contains(ChannelNames,'Time'))
    idx=or(contains(ChannelNames,'time'),contains(ChannelNames,'Time'));
    time=data{idx};
    data(idx)=[];
    ChannelNames(idx)=[];
    LChans=length(ChannelNames);
else
    type=split(textdata.name,'_');
    type=type{1};
    type=type(3:end);
    switch type
        case 'IVDaq'
            rate=1./textdata.Interpoint_Time;
        case {'OscilloDaq','OscilloDaqTimed'}
            rate=textdata.rate;
        case 'MemoryTask2'
            rate=textdata.Freq_Sampling__Hz_;
        case 'MultiWaveTask'
            rate=textdata.Acq_Rate__s_;
        otherwise
            rate=textdata.rate;
    end
    length_Data=length(data{1});
    timeInt=linspace(0,length_Data-1,length_Data);
    time=(1/str2double(rate)).*timeInt;
end
end
