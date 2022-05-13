%multi-optimization of GA


clear; close all;
ver=9;%

%% load seed (initial set of parameters)
load('MultiGA_Output66.mat')
para=val(7,:);


startPara = para;
startFval = fitness2(startPara);

[OutputPara1 , OutputFval1]= GA(25, 20,startPara, startFval,ver);
[OutputPara2 , OutputFval2]= GA(30, 25,OutputPara1 , OutputFval1,ver);
[OutputPara3 , OutputFval3]= GA(30, 30,OutputPara2 , OutputFval2,ver);
[OutputPara4 , OutputFval4]= GA(30, 25,OutputPara3 , OutputFval3,ver);