function [Y, time, y0_,TEOUT1,IEOUT1,Y_output1,time1,Y_output2,time2,DNArep,TDNAini]=main1(y0,celltype,ver,mutant)

global T_e1 T_term
if strcmp(celltype,'SW')
    CHECK = 300;
else 
    CHECK = 270;
end
T_e1=CHECK;%5/12 300-90 max time of Tini this study can detect

T_Sphase=90;%S-phase is arbitrarily set as 90 min

y0first=y0;
if isempty(celltype)
    celltype='SW';
end

global p;
%% load parameters
parameters(1,ver,mutant);%all parameters
TS=p.TS; %threshold of CtrA~P determing the time of chrommosome replication initiation





 %% simulation
 if strcmp(celltype,'SW')
  tspan=400; %max sim time
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
     T_e1 = 15;
 end
         elseif ie  ==  2%DNA replication initiates
        if te<T_e1 && te>15
        T_e1 = te;
        elseif te<=15% if replication initiates before 15min
            T_e1=15;
        end
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
    T_term = T_e1+T_Sphase;
    t0 = t(nt);
    if t0 >= tf
        break;
    end
end

%% 1027 STALKED CELL
 elseif strcmp(celltype,'ST')
     tspan=CHECK;%longer than expected wt st cell cycle
     t0  =  0;       % Start time
tf  = tspan;%     % End time  

xoverFcn = @(t,y) podJ_eventST(t,y,TS);
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
     if ie  ==  1%time >0.0001min
        %check if average [CtrA~P] is lower than TS        
         CtrAPV  =y(:,81:90)+y(:,321:330)+y(:,331:340);
         CtrAPV=CtrAPV';
         CtrAPV = sum(CtrAPV)./10;
         CtrAPV = CtrAPV(end);

         if mean(CtrAPV)<=TS
         T_e1 = 0;
         end
      elseif ie  ==  2%DNA replication initiates
          if T_e1 < 5
          else
            if te>5
            warning('st cell does not initiates DNA replication before 5 min after cell separation')
%             else
              T_e1 = te;
%               mean(CtrAPV) - TS
            end
          end
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
    T_term = T_e1+T_Sphase;
    t0 = t(nt);
    if t0 >= tf
        break;
    end
end
 end
 %% simulation results-1
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
M1=M1'; %grid for spatial plotting
TEOUT1=teout; IEOUT1=ieout;
%% cut the simulation [0 z-ring-closed]
if T_e1<CHECK
    DNArep=1;%DNA replication initiates before 300 min
[~,Index] = min(abs(time1-(T_term+5)));%find the Z-ring closed time in time1
time1=time1(1:Index);
Y_output1=Y_output1(:,1:Index);
M1=M1(:,1:Index);
%% output IE and TE
time1End=time1(end);
TE1end=find(TEOUT1<=time1End);
TE1end=TE1end(end);
TEOUT1=TEOUT1(1:TE1end);
IEOUT1=IEOUT1(1:TE1end);

TDNAini = T_term - 90; %output DNA initiation time
%% simulation results after z-ring closed

y0 = Y_output1(:,end);
Output2 = main_DIV(y0,ver,mutant);
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
    [~,Index] = min(abs(time1-CHECK));%show plot during 0-180 min
    time1=time1(1:Index);
    Y_output1=Y_output1(:,1:Index);
    M1=M1(:,1:Index);
    Y=Y_output1;
    M=M1;
    time=time1;
    y0_=y0first;
    Y_output2='None';
    time2='None';
    TDNAini = 'None';
    warning('DNA does not initiate before CHECK time')
end
