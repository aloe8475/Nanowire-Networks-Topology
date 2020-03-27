function UpdateSlider(handles)
Sel=handles.PatternList.Value;
Sel=Sel(1);

T=handles.refList{1}.timedata{Sel};
handles.TimeSlider.Value=1;
handles.TimeSlider.Min=1;
handles.TimeSlider.Max=length(T);
handles.TimeSlider.SliderStep=[0.01 0.1];
handles.TimeSlider.UserData=T;

end