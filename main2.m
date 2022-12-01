function [Y, time,TEOUT,IEOUT]=main(y0,celltype,ver,mutant,tspan)
%main2
%if Z-ring does not closed or closed later?

global T_e1 T_term
T_e1=tspan;%5/12 300-90 max time of Tini this study can detect
T_Sphase=90;%normal is 90

y0first=y0;
if isempty(celltype)
    celltype='SW';
end

global p;

parameters(1,ver,mutant);%all parameters
TS=p.TS; 

 %% simulation
 if strcmp(celltype,'SW')
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
     T_term = T_e1+T_Sphase;
 end
         elseif ie  ==  2%DNA replication initiates
          if te>15  %can take place multiple times
            T_e1 = te;
          elseif te<=15% if replication initiates before 15min
            T_e1=15;
          end
        T_term = T_e1+T_Sphase;
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
%% 1028 stalked cell
 elseif strcmp(celltype, 'ST')
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
              T_e1 = te;
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
%% if no Z-ring closed
    Y=Y_output1;
    time=time1;
TEOUT=TEOUT1;   IEOUT=IEOUT1;


end
