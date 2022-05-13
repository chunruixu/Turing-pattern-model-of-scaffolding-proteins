function [Y, time, y0_,T_e1]=main(y0,para,ver,mutant)
global T_e1 T_term
T_e1=60;%150;
T_Sphase=90;
global TS
TS=para(41);
global p
load_para(ver,mutant);

tspan=170;%125;
t0  =  0;       % Start time
tf  = tspan;%     % End time  Z-ring closed


tout = t0;
y0 = y0.';
yout = y0;

teout  =  [];
yeout  =  [];
ieout  =  [];
xoverFcn = @(t,y) podJ_event(t,y,TS);
options  =  odeset('Events',xoverFcn,'RelTol',1e-4,'AbsTol',1e-6);
% options = odeset('Events', @events5RegCPLX, 'RelTol', 1e-5, 'AbsTol', 1e-7);
odefun = @(t, y) ODE12(t, y, para);

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
        
CtrAPV  =y(:,81:90)+y(:,321:330)+y(:,331:340);
CtrAPV=CtrAPV';
CtrAPV = sum(CtrAPV)./10;
n=length(CtrAPV);
if n>=2
CtrAPV = CtrAPV(n-2:n);
else
    CtrAPV = CtrAPV(end);
end
 if mean(CtrAPV)<=TS
     T_e1 = 15;%min(T_e1,15);
%      fprintf('T_e1= %8.5f\n',T_e1)
     T_term = T_e1+T_Sphase;
 end
         elseif ie  ==  2%DNA replication initiates
        if te<T_e1 && te>15
        T_e1 = te;
        elseif te<=15
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
[~,Index] = min(abs(time1-T_term+5));
time1=time1(1:Index);
Y_output1=Y_output1(:,1:Index);
M1=M1(:,1:Index);

%% simulation results after z-ring closed

y0 = Y_output1(:,end);
Output2 = main_DIV(y0,para,ver,mutant);
% Output2 = main_DIV1(T,y0,ver,mutant);
Y_output2 = Output2.yout;
Y_output2 = Y_output2';
M2=Output2.grid;
time2=Output2.time;


Y=[Y_output1 Y_output2];
M=[M1 M2];
Tend1=time1(end);
time=[time1 ;time2+Tend1];
y0_=Y(:,end);
end