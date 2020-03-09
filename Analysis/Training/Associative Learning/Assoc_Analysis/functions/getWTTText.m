function [WaitEdit,TrainEdit,TestEdit]=getWTTText(metadata)

WaitTime=metadata.WaitTime_s_;
WaitCond=metadata.WaitConditions;

WaitEdit=strcat('WaitTime: ',WaitTime,' s. WaitCond: ',WaitCond);

TestTime=metadata.TestTime;
TestVolt=metadata.TestVolt;

TestEdit=strcat('TestTime: ',TestTime,' s. TestVolt: ',TestVolt);

TrainTime=metadata.TimeOut;

TrainEdit=strcat('TrainTimeout: ',TrainTime);




end