function [Y, time, y0_,TEOUT1,IEOUT1,Y_output1,time1,Y_output2,time2,DNArep]=main1(y0,celltype,ver,mutant)
%main
%11/26/2020

global T_e1 T_term
% T_e1=60;%150;
TDETECT=300;
T_e1=TDETECT;%5/12 300-90 max time of Tini this study can detect
T_Sphase=90;

y0first=y0;
if isempty(celltype)
    celltype='SW';
end

global p;
% parameters(1,ver,'WT');%all parameters
parameters(1,ver,mutant);%all parameters
TS=p.TS; 

 %% simulation
 if strcmp(celltype,'SW')
    % tspan=155;  %60(greatest initiation time)+90(replication time cost)+5(time for Z-ring closed after replication termination)
  tspan=400; %5/12 max sim time
t0  =  0;       % Start time
tf  = tspan;%     % End time  

xoverFcn = @(t,y) podJ_event(t,y,TS);
options  =  odeset('Events',xoverFcn,'RelTol',1e-4,'AbsTol',1e-6);
odefun = @(t, y) ODE_CpdR(t, y);
tout = t0;
y0 = y0.';
yout = y0;

teout  =  [];
yeout  =  [];
ieout  =  [];

while t0<tf
[t,y,te,ye,ie] = ode15s(odefun,[t0 tf],y0,options);
 nt = length(t);
%     
 tout = [tout;t(2:nt)]; 
 yout = [yout;y(2:nt,:)]; 
 teout  =  [teout;te];  
 yeout  =  [yeout;ye]; ieout  =  [ieout;ie];
    y0  =  y(nt,:);
     if isscalar(ie)  ==  0
        ie  =  0;
     end
     if ie  ==  1%time >15min
%check if average [CtrA~P] is lower than TS at t=15min        
CtrAPV  =y(:,81:90)+y(:,321:330)+y(:,331:340);
CtrAPV=CtrAPV';
CtrAPV = sum(CtrAPV)./10;
n=length(CtrAPV);

CtrAPV = CtrAPV(end);

 if mean(CtrAPV)<=TS
     T_e1 = 15;%min(T_e1,20);
%      fprintf('T_e1= %8.5f\n',T_e1)
     T_term = T_e1+T_Sphase;
 end
         elseif ie  ==  2%DNA replication initiates
        if te<T_e1 && te>15
        T_e1 = te;
        elseif te<=15% if replication initiates before 15min
            T_e1=15;
        end
%         fprintf('T_e1= %8.5f\n',T_e1)
        T_term = T_e1+T_Sphase;
%     elseif ie == 2%Fork passes divL
%         y0(395)=1;%SdivL
    elseif ie == 3%Fork passes ctrA
y0(392)=1;%SctrA changed from 0 to 1
    elseif ie==4
y0(393)=1;%SpleC
    elseif ie==5
        y0(394)=1;%SperP
    elseif ie ==6
        y0(391)=1;%SpodJ
         elseif ie==7
       y0(391:395)=0;
    end
    t0 = t(nt);
    if t0 >= tf
        break;
    end
end
 end
 %% simulation results1
 Y_output1=yout';
 time1=tout;
 L=length(yout(:,1));
 M1(:,6)=zeros(1,L);
M1(:,5)=-yout(:,396);
M1(:,4)=-2*yout(:,396);
M1(:,3)=-3*yout(:,396);
M1(:,2)=-4*yout(:,396);
M1(:,1)=-5*yout(:,396);
M1(:,7:11)=-fliplr(M1(:,1:5));
M1=M1';
TEOUT1=teout; IEOUT1=ieout;
%% cut the simulation [0 z-ring-closed]
if T_e1<300
    DNArep=1;
[~,Index] = min(abs(time1-T_term+5));
time1=time1(1:Index);
Y_output1=Y_output1(:,1:Index);
M1=M1(:,1:Index);

T_term - 90 %output DNA initiation time
%% simulation results after z-ring closed

y0 = Y_output1(:,end);
Output2 = main_DIV(y0,ver,mutant);
% Output2 = main_DIV1(T,y0,ver,mutant);
Y_output2 = Output2.yout;
Y_output2 = Y_output2';
M2=Output2.grid;
time2=Output2.time;

Tend1=time1(end);
Y=[Y_output1 Y_output2];
M=[M1 M2];
time=[time1 ;time2+Tend1];
y0_=Y(:,end);
else   %DNA is not initiated before 300 min;
    DNArep=0;
    [~,Index] = min(abs(time1-300));%show plot during 0-180 min
    time1=time1(1:Index);
    Y_output1=Y_output1(:,1:Index);
    M1=M1(:,1:Index);
    Y=Y_output1;
    M=M1;
    time=time1;
    y0_=y0first;
    Y_output2='None';
    time2='None';
end
%% plot
% graphcellcycle2(Y,time)