close all;

filename = '\article\figures\beta vs cond\workspace - T=3e1, dt=1e-4, V=1.3.mat';
load(filename);
timeStamps = [0, 11.6 ;
              12.3 , 12.9  ;
              13.43, 15    ;
              16.21, 24.41 ;
              24.51, 26.15 ;
              26.3 , 29.9  ];
pinkNoiseAnalysis(t,conductance,timeStamps,true);