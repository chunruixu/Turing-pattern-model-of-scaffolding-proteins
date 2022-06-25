function output = main_DIV(y0,para,ver,mutant)
% close all
%1/24/2021
% clear all
global p
load_para(ver,mutant);
p.syn_podJ2=para(1);  p.Ji_PodJCtrA=para(2);  p.syn_spmx=para(3);  p.Ja_SpmXCtrA=para(4); p.syn_ctrA1=para(5); 
p.syn_ctrA2=para(6);  p.Ja_CtrACtrA=para(7);  p.deg_ctrA2=para(8); p.Jd_CpdR=para(9); p.syn_pleC=para(10); 
p.fb_PleC=para(11); p.b_PleCDivKP=para(12); p.ub_DivJDivKP=para(13); p.ub_DivJDivK=para(14); p.syn_divK2=para(15); 
p.Ja_DivKCtrA=para(16); p.pho_DivKPleC=para(17);  p.pho_DivKDivJf=para(18); p.pho_DivKDivJb=para(19);  p.depho_DivK=para(20);
p.syn_perP=para(21);  p.Ja_PerPCtrA=para(22); p.syn_divL=para(23);  p.deg_divL=para(24);  p.alpha_DivLPopZ=para(25); 
p.alpha_DivLPodJ=para(26);  p.b_DivLDivKP=para(27); p.syn_cckA=para(28); p.alpha_CckAPopZ=para(29); p.alpha_CckADivL=para(30);
p.pk_cckA1=para(31);  p.kp_cckA2=para(32); p.pk_cckA2=para(33); p.ub_CtrACckAkin=para(34); p.ub_CtrAPCckAph=para(35);
p.syn_cpdR=para(36);  p.deg_cpdR=para(37);  p.JaCpdRCtrA=para(38);  p.fb_cpdR=para(39);  p.phoCpdR =para(40);

p.ph2kin_PleC1=para(42);
p.ph2kin_PleC2=0.1*para(42);

%% Initial values
y0(391:395)=0;%all S ia assumed to be 0 after Z-ring closed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integration parameters
t0  =  0;       % Start time
% tf  = 30;
tf=25;%125+25=150

[tout,yout] = ode15s(@ODE_CpdR_DIV,[t0 tf],y0);
 %% construct grid for plotting
L=length(yout(:,1));
M(:,6)=zeros(1,L);
M(:,5)=-yout(:,396);
M(:,4)=-2*yout(:,396);
M(:,3)=-3*yout(:,396);
M(:,2)=-4*yout(:,396);
M(:,1)=-5*yout(:,396);
M(:,7:11)=-fliplr(M(:,1:5));
M=M';    


output.time = tout;
output.grid = M;
 output.yout = yout;

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%6/8/2021
function dydt = ODE_CpdR_DIV(t,y)

dydt=zeros(396,1);
m=0.1;%0.15;
mplec=0.15;%0.15;
mpodj=0.1;
mctra=0.1;
%% PodJL_m
m_podj=((1-m)*y(391)+m);
% Bin 1
dydt(1) = 0*p.syn_podJ*m_podj+0*p.syn_podJ2*p.Ji_PodJCtrA^p.n_PodJCtrA/(p.Ji_PodJCtrA^p.n_PodJCtrA+y(1+80)^p.n_PodJCtrA)...%CtrAP 81-90
    - (p.deg_podJ1+p.deg_podJm*y(1+210))*y(1)-p.mu*y(1)... % PerP PerP 211-220
    - p.dnv_podJ*y(1)...% denovo polymerization; SpmX (31-40&41-50)
    - p.aut1_podJ/(1+p.alpha_PodJSpmX*(y(1+30)+y(1+40)))*y(1+10)^p.podj*y(1)...% autocatalytic polymerization at poles
+ p.depol_podJ*y(1+10)...%
    + p.D_podJm*(y(2)-y(1))/(y(396)^2); %
% Bin 2 - 9
for i=2:9
    if i==5
                dydt(i) = p.syn_podJ*m_podj+p.syn_podJ2*p.Ji_PodJCtrA^p.n_PodJCtrA/(p.Ji_PodJCtrA^p.n_PodJCtrA+y(i+80)^p.n_PodJCtrA)...
    - (p.deg_podJ1+p.deg_podJm*y(i+210))*y(i)-p.mu*y(i)... % PerP 211-220
    - p.dnv_podJ*y(i)...% denovo polymerization; SpmX 
    - p.aut1_podJ1/(1+p.alpha_PodJSpmX*(y(i+40)+y(i+30)))*y(i+10)^p.podj*y(i)...% autocatalytic polymerization at poles
+ p.depol_podJ*y(i+10)...%
    + p.D_podJm*(y(i-1)-y(i))/(y(396)^2); %
    elseif i==6
                dydt(i) = p.syn_podJ*m_podj+p.syn_podJ2*p.Ji_PodJCtrA^p.n_PodJCtrA/(p.Ji_PodJCtrA^p.n_PodJCtrA+y(i+80)^p.n_PodJCtrA)...
    - (p.deg_podJ1+p.deg_podJm*y(i+210))*y(i)-p.mu*y(i)... % PerP 211-220
    - p.dnv_podJ*y(i)...% denovo polymerization; SpmX 
    - p.aut1_podJ1/(1+p.alpha_PodJSpmX*(y(i+40)+y(i+30)))*y(i+10)^p.podj*y(i)...% autocatalytic polymerization at poles
+ p.depol_podJ*y(i+10)...%
    + p.D_podJm*(y(i+1)-y(i))/(y(396)^2); %
    else
        dydt(i) = p.syn_podJ*m_podj+p.syn_podJ2*p.Ji_PodJCtrA^p.n_PodJCtrA/(p.Ji_PodJCtrA^p.n_PodJCtrA+y(i+80)^p.n_PodJCtrA)...
    - (p.deg_podJ1+p.deg_podJm*y(i+210))*y(i)-p.mu*y(i)... % PerP 211-220
    - p.dnv_podJ*y(i)...% denovo polymerization; SpmX 
    - p.aut1_podJ1/(1+p.alpha_PodJSpmX*(y(i+40)+y(i+30)))*y(i+10)^p.podj*y(i)...% autocatalytic polymerization at poles
+ p.depol_podJ*y(i+10)...%
    + p.D_podJm*((y(i-1)-2*y(i)+y(i+1)))/(y(396)^2); %    
    end

end
% Bin 10
dydt(10) = 0*p.syn_podJ*m_podj+0*p.syn_podJ2*p.Ji_PodJCtrA^p.n_PodJCtrA/(p.Ji_PodJCtrA^p.n_PodJCtrA+y(10+80)^p.n_PodJCtrA)...%CtrAP 81-90
    - (p.deg_podJ1+p.deg_podJm*y(10+210))*y(10)-p.mu*y(10)... % PerP PerP 211-220
    - p.dnv_podJ*y(10)...% denovo polymerization; SpmX (31-40&41-50)
    - p.aut1_podJ/(1+p.alpha_PodJSpmX*(y(10+30)+y(10+40)))*y(10+10)^p.podj*y(10)...% autocatalytic polymerization at poles
+ p.depol_podJ*y(10+10)...%
    + p.D_podJm*(y(9)-y(10))/(y(396)^2); %
 

%% PodJL_p
% Bin 11 
dydt(11) =-(p.deg_podJ1+p.deg_podJp*y(11+200))*y(11)...%PerP:211-220
    + p.dnv_podJ*y(11-10)...%SpmX (31-40&41-50)
    + p.aut1_podJ/(1+p.alpha_PodJSpmX*(y(11+20)+y(11+30)))*y(11)^p.podj*y(11-10)...% autocatalytic polymerization at poles
    - p.depol_podJ*y(11)...%
    + p.D_podJL*(y(12)-y(11))/(y(396)^2)- p.mu*y(11); % polymer diffusion % dilution
% Bin 12-19
for i=12:19
    if i==15
    dydt(i) =-(p.deg_podJ1+p.deg_podJp*y(i+200))*y(i)...%PerP:2i-220
    + p.dnv_podJ*y(i-10)...%SpmX (31-40&41-50)
    + p.aut1_podJ1/(1+p.alpha_PodJSpmX*(y(i+30)+y(i+20)))*y(i)^p.podj*y(i-10)...% autocatalytic polymerization at poles
    - p.depol_podJ*y(i)...%
    + p.D_podJL/(y(396)^2)*(y(i-1)-y(i))- p.mu*y(i); % polymer diffusion % dilution
    elseif i==16
            dydt(i) =-(p.deg_podJ1+p.deg_podJp*y(i+200))*y(i)...%PerP:2i-220
    + p.dnv_podJ*y(i-10)...%SpmX (31-40&41-50)
    + p.aut1_podJ1/(1+p.alpha_PodJSpmX*(y(i+30)+y(i+20)))*y(i)^p.podj*y(i-10)...% autocatalytic polymerization at poles
    - p.depol_podJ*y(i)...%
    + p.D_podJL/(y(396)^2)*(y(i+1)-y(i))- p.mu*y(i); % polymer diffusion % dilution
    else
            dydt(i) =-(p.deg_podJ1+p.deg_podJp*y(i+200))*y(i)...%PerP:2i-220
    + p.dnv_podJ*y(i-10)...%SpmX (31-40&41-50)
    + p.aut1_podJ1/(1+p.alpha_PodJSpmX*(y(i+30)+y(i+20)))*y(i)^p.podj*y(i-10)...% autocatalytic polymerization at poles
    - p.depol_podJ*y(i)...%
    + p.D_podJL/(y(396)^2)*(y(i+1)-2*y(i)+y(i-1))- p.mu*y(i); % polymer diffusion % dilution
    end
end
%Bin 20
dydt(20) =-(p.deg_podJ1+p.deg_podJp*y(20+200))*y(20)...%PerP:211-220
    + p.dnv_podJ*y(20-10)...%SpmX (31-40&41-50)
    + p.aut1_podJ/(1+p.alpha_PodJSpmX*(y(20+20)+y(20+30)))*y(20)^p.podj*y(20-10)...% autocatalytic polymerization at poles
    - p.depol_podJ*y(20)...%
    + p.D_podJL*(y(19)-y(20))/(y(396)^2)- p.mu*y(20); % polymer diffusion % dilution        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PodJS
% no diffussion assumed
for i=21:30
dydt(i)=(p.deg_podJ1+p.deg_podJm*y(i+190))*(y(i-20)+y(i-10))...%PerP 211-220
    -(p.mu+p.deg_s)*y(i);
end

%% SpmXm
%Bin 31
q=2;
dydt(31) = 0*p.syn_spmx*y(31+50)^q/(y(31+50)^q+p.Ja_SpmXCtrA^q)-(p.deg_spmx+p.mu)*y(31)...%CtrAP 81-90
    -p.dnv_spmx*(1+p.alpha_SpmXPopZ*(y(31+20)+y(31+30)))*y(31)+p.depol_spmx*y(31+10)...
    -p.aut_spmx*y(31)*y(31+10)^2 ...%PopZp 61-70  PopZm 51-60
+p.D_spmx*(y(32)-y(31))/(y(396)^2);
%Bin 32-39
for i=32:39
    if i==35
dydt(i) = p.syn_spmx*y(i+50)^q/(y(i+50)^q+p.Ja_SpmXCtrA^q)-(p.deg_spmx+p.mu)*y(i)...%CtrAP 81-90
    -p.dnv_spmx*(1+p.alpha_SpmXPopZ*(y(i+20)+y(i+30)))*y(i)+p.depol_spmx*y(i+10)...
    -p.aut_spmx*y(i)*y(i+10)^2 ...%PopZp 61-70  PopZm 51-60
+p.D_spmx/(y(396)^2)*(-y(i)+y(i-1));
    elseif i==36
        dydt(i) = p.syn_spmx*y(i+50)^q/(y(i+50)^q+p.Ja_SpmXCtrA^q)-(p.deg_spmx+p.mu)*y(i)...%CtrAP 81-90
    -p.dnv_spmx*(1+p.alpha_SpmXPopZ*(y(i+20)+y(i+30)))*y(i)+p.depol_spmx*y(i+10)...
    -p.aut_spmx*y(i)*y(i+10)^2 ...%PopZp 61-70  PopZm 51-60
+p.D_spmx/(y(396)^2)*(y(i+1)-y(i));
    else 
        dydt(i) = p.syn_spmx*y(i+50)^q/(y(i+50)^q+p.Ja_SpmXCtrA^q)-(p.deg_spmx+p.mu)*y(i)...%CtrAP 81-90
    -p.dnv_spmx*(1+p.alpha_SpmXPopZ*(y(i+20)+y(i+30)))*y(i)+p.depol_spmx*y(i+10)...
    -p.aut_spmx*y(i)*y(i+10)^2 ...%PopZp 61-70  PopZm 51-60
+p.D_spmx/(y(396)^2)*(y(i+1)-2*y(i)+y(i-1));
    end
end
%Bin 40
dydt(40) = 0*p.syn_spmx*y(40+50)^q/(y(40+50)^q+p.Ja_SpmXCtrA^q)-(p.deg_spmx+p.mu)*y(40)...%CtrAP 81-90
    -p.dnv_spmx*(1+p.alpha_SpmXPopZ*(y(40+20)+y(40+30)))*y(40)+p.depol_spmx*y(40+10)...
    -p.aut_spmx*y(40)*y(40+10)^2 ...%PopZp 61-70  PopZm 51-60
+p.D_spmx*(y(39)-y(40))/(y(396)^2);

%% SpmXp
%Bin 41
dydt(41)=-(p.deg_spmx+p.mu)*y(41)...
    +p.dnv_spmx*(1+p.alpha_SpmXPopZ*(y(41+10)+y(41+20)))*y(41-10)-p.depol_spmx*y(41)...
      +p.aut_spmx*y(41-10)*y(41)^2 ...%PopZp 61-70  PopZm 51-60
     +p.D_spmxp*(y(42)-y(41))/(y(396)^2);

 %Bin 42-49
 for i=42:49
     if i==45
     dydt(i)=-(p.deg_spmx+p.mu)*y(i)...
    +p.dnv_spmx*(1+p.alpha_SpmXPopZ*(y(i+10)+y(i+20)))*y(i-10)-p.depol_spmx*y(i)...
      +p.aut_spmx*y(i-10)*y(i)^2 ...%
       +p.D_spmxp/(y(396)^2)*(-y(i)+y(i-1));
     elseif i==46
              dydt(i)=-(p.deg_spmx+p.mu)*y(i)...
    +p.dnv_spmx*(1+p.alpha_SpmXPopZ*(y(i+10)+y(i+20)))*y(i-10)-p.depol_spmx*y(i)...
      +p.aut_spmx*y(i-10)*y(i)^2 ...%
       +p.D_spmxp/(y(396)^2)*(y(i+1)-y(i));
     else
              dydt(i)=-(p.deg_spmx+p.mu)*y(i)...
    +p.dnv_spmx*(1+p.alpha_SpmXPopZ*(y(i+10)+y(i+20)))*y(i-10)-p.depol_spmx*y(i)...
      +p.aut_spmx*y(i-10)*y(i)^2 ...%
       +p.D_spmxp/(y(396)^2)*(y(i+1)-2*y(i)+y(i-1));
     end
 end
%Bin 50
dydt(50)=-(p.deg_spmx+p.mu)*y(50)...
    +p.dnv_spmx*(1+p.alpha_SpmXPopZ*(y(50+10)+y(50+20)))*y(50-10)-p.depol_spmx*y(50)...
      +p.aut_spmx*y(50-10)*y(50)^2 ...%PopZp 61-70  PopZm 51-60
     +p.D_spmxp*(y(49)-y(50))/(y(396)^2);

%% PopZm
%Bin 51
dydt(51) =0*p.syn_popz  - (p.mu+p.deg_popzm)*y(51)...
    -p.dnv_popz*(1+p.alpha_PopZPodJ*(y(51-50)+y(51-40)))*y(51)-p.aut_popz*y(51+10)^2*y(51)...
    +p.depol_popz*y(51+10)...%PodJp 11-20& PodJm 1-10
    +p.D_popzm*(y(52)-y(51))/(y(396)^2);

%Bin 52-59
for i=52:59
    if i==55
dydt(i) =p.syn_popz  - (p.mu+p.deg_popzm)*y(i)...
    -p.dnv_popz*(1+p.alpha_PopZPodJ*(y(i-50)+y(i-40)))*y(i)-p.aut_popz1*y(i+10)^2*y(i)...%
    +p.depol_popz*y(i+10)...
       +p.D_popzm/(y(396)^2)*(-y(i)+y(i-1));
    elseif i==56
        dydt(i) =p.syn_popz  - (p.mu+p.deg_popzm)*y(i)...
    -p.dnv_popz*(1+p.alpha_PopZPodJ*(y(i-50)+y(i-40)))*y(i)-p.aut_popz1*y(i+10)^2*y(i)...%
    +p.depol_popz*y(i+10)...
       +p.D_popzm/(y(396)^2)*(y(i+1)-y(i));
    else
        dydt(i) =p.syn_popz  - (p.mu+p.deg_popzm)*y(i)...
    -p.dnv_popz*(1+p.alpha_PopZPodJ*(y(i-50)+y(i-40)))*y(i)-p.aut_popz1*y(i+10)^2*y(i)...%
    +p.depol_popz*y(i+10)...
       +p.D_popzm/(y(396)^2)*(y(i+1)-2*y(i)+y(i-1));
    end
end
%Bin 60
dydt(60) =0*p.syn_popz  - (p.mu+p.deg_popzm)*y(60)...
    -p.dnv_popz*(1+p.alpha_PopZPodJ*(y(60-50)+y(60-40)))*y(60)-p.aut_popz*y(60+10)^2*y(60)...
    +p.depol_popz*y(60+10)...%PodJp 11-20& PodJm 1-10
    +p.D_popzm*(y(59)-y(60))/(y(396)^2);


%% PopZp
%Bin 61
dydt(61) = -(p.mu+p.deg_popzp)*y(61)...
    +p.dnv_popz*(1+p.alpha_PopZPodJ*(y(61-60)+y(61-50)))*y(61-10)...
    +p.aut_popz*y(61)^2*y(61-10)...%PodJp 11-20& PodJm 1-10
    -p.depol_popz*y(61)...
 +p.D_popzp*(y(62)-y(61))/(y(396)^2);

%Bin 62-69
for i=62:69
    if i==65
      dydt(i) = -(p.mu+p.deg_popzp)*y(i)...
    +p.dnv_popz*(1+p.alpha_PopZPodJ*(y(i-60)+y(i-50)))*y(i-10)...
    +p.aut_popz1*y(i)^2*y(i-10)...%PodJp 11-20& PodJm 1-10
    -p.depol_popz*y(i)...
 +p.D_popzp/(y(396)^2)*(-y(i)+y(i-1));
    elseif i==66
              dydt(i) = -(p.mu+p.deg_popzp)*y(i)...
    +p.dnv_popz*(1+p.alpha_PopZPodJ*(y(i-60)+y(i-50)))*y(i-10)...
    +p.aut_popz1*y(i)^2*y(i-10)...%PodJp 11-20& PodJm 1-10
    -p.depol_popz*y(i)...
 +p.D_popzp/(y(396)^2)*(y(i+1)-y(i));
    else
              dydt(i) = -(p.mu+p.deg_popzp)*y(i)...
    +p.dnv_popz*(1+p.alpha_PopZPodJ*(y(i-60)+y(i-50)))*y(i-10)...
    +p.aut_popz1*y(i)^2*y(i-10)...%PodJp 11-20& PodJm 1-10
    -p.depol_popz*y(i)...
 +p.D_popzp/(y(396)^2)*(y(i+1)-2*y(i)+y(i-1));
    end
end

%Bin 70
dydt(70) = -(p.mu+p.deg_popzp)*y(70)...
    +p.dnv_popz*(1+p.alpha_PopZPodJ*(y(70-60)+y(70-50)))*y(70-10)...
    +p.aut_popz*y(70)^2*y(70-10)...%PodJp 11-20& PodJm 1-10
    -p.depol_popz*y(70)...
 +p.D_popzp*(y(69)-y(70))/(y(396)^2);

%% CtrA
m_ctrA = (1-mctra)*y(392)+mctra; %S_ctrA=y(392)
% CckAfkin: 281-290; CckAbkin: 291-300; CckAfph: 301-310; CckAbph: 311-320
%CtrAP:CckAfph 321-330; CtrAP:CckAbph 331-340;
%CtrAu:CckAfkin 341-350; CtrAu:CckAbkin 351-360

%CpdRf 361-370; CpdRb 371-380
%Bin 71
dydt(71)=0*p.syn_ctrA1*(1-y(71+10)/(y(71+10)+y(71)+p.Ji_CtrACtrA))*m_ctrA...
    +0*p.syn_ctrA2*y(71+10)/(y(71+10)+y(71)+p.Ja_CtrACtrA) ...
    -(p.mu+p.deg_ctrA1+p.deg_ctrA2*(y(71+290)+y(71+300))/(p.Jd_CpdR+(y(71+290)+y(71+300))))*y(71) ...%
    -p.b_CtrACckAkin*(y(71+210)+y(71+220))*y(71)+(p.ub_CtrACckAkin+p.deg_cckA)*(y(71+270)+y(71+280))...
    +p.dephoCtrA*(y(71+250)+y(71+260))...
+p.D_CtrA*(y(72)-y(71))/(y(396)^2);

%Bin 72-79
for i=72:79
    if i==75
    dydt(i)=p.syn_ctrA1*(1-y(i+10)/(y(i+10)+y(i)+p.Ji_CtrACtrA))*m_ctrA...
    +p.syn_ctrA2*y(i+10)/(y(i+10)+y(i)+p.Ja_CtrACtrA) ...
    -(p.mu+p.deg_ctrA1+p.deg_ctrA2*(y(i+290)+y(i+300))/(p.Jd_CpdR+(y(i+290)+y(i+300))))*y(i) ...%
    -p.b_CtrACckAkin*(y(i+210)+y(i+220))*y(i)+(p.ub_CtrACckAkin+p.deg_cckA)*(y(i+270)+y(i+280))...
    +p.dephoCtrA*(y(i+250)+y(i+260))...
 +p.D_CtrA/(y(396)^2)*(-y(i)+y(i-1));
    elseif i==76
            dydt(i)=p.syn_ctrA1*(1-y(i+10)/(y(i+10)+y(i)+p.Ji_CtrACtrA))*m_ctrA...
    +p.syn_ctrA2*y(i+10)/(y(i+10)+y(i)+p.Ja_CtrACtrA) ...
    -(p.mu+p.deg_ctrA1+p.deg_ctrA2*(y(i+290)+y(i+300))/(p.Jd_CpdR+(y(i+290)+y(i+300))))*y(i) ...%
    -p.b_CtrACckAkin*(y(i+210)+y(i+220))*y(i)+(p.ub_CtrACckAkin+p.deg_cckA)*(y(i+270)+y(i+280))...
    +p.dephoCtrA*(y(i+250)+y(i+260))...
 +p.D_CtrA/(y(396)^2)*(y(i+1)-y(i));
    else
            dydt(i)=p.syn_ctrA1*(1-y(i+10)/(y(i+10)+y(i)+p.Ji_CtrACtrA))*m_ctrA...
    +p.syn_ctrA2*y(i+10)/(y(i+10)+y(i)+p.Ja_CtrACtrA) ...
    -(p.mu+p.deg_ctrA1+p.deg_ctrA2*(y(i+290)+y(i+300))/(p.Jd_CpdR+(y(i+290)+y(i+300))))*y(i) ...%
    -p.b_CtrACckAkin*(y(i+210)+y(i+220))*y(i)+(p.ub_CtrACckAkin+p.deg_cckA)*(y(i+270)+y(i+280))...
    +p.dephoCtrA*(y(i+250)+y(i+260))...
 +p.D_CtrA/(y(396)^2)*(y(i+1)-2*y(i)+y(i-1));
    end
end
dydt(80)=0*p.syn_ctrA1*(1-y(80+10)/(y(80+10)+y(80)+p.Ji_CtrACtrA))*m_ctrA...
    +0*p.syn_ctrA2*y(80+10)/(y(80+10)+y(80)+p.Ja_CtrACtrA) ...
    -(p.mu+p.deg_ctrA1+p.deg_ctrA2*(y(80+290)+y(80+300))/(p.Jd_CpdR+(y(80+290)+y(80+300))))*y(80) ...%
    -p.b_CtrACckAkin*(y(80+210)+y(80+220))*y(80)+(p.ub_CtrACckAkin+p.deg_cckA)*(y(80+270)+y(80+280))...
    +p.dephoCtrA*(y(80+250)+y(80+260))...
+p.D_CtrA*(y(79)-y(80))/(y(396)^2);


%% CtrAP
% CckAfkin: 281-290; CckAbkin: 291-300; CckAfph: 301-310; CckAbph: 311-320
%CtrAP:CckAfph 321-330; CtrAP:CckAbph 331-340;
%CtrAu:CckAfkin 341-350; CtrAu:CckAbkin 351-360

%CpdRf 361-370; CpdRb 371-380
%Bin 81
dydt(81)=-(p.mu+p.deg_ctrA1+p.deg_ctrA2*(y(81+280)+y(81+290))/(p.Jd_CpdR+(y(81+280)+y(81+290))))*y(81) ...%
      -p.b_CtrAPCckAph*(y(81+220)+y(81+230))*y(81)+(p.ub_CtrAPCckAph+p.deg_cckA)*(y(81+240)+y(81+250))...
    +p.phoCtrA*(y(81+260)+y(81+270))-p.ukdephoCtrA*y(81)...
    + p.D_CtrAP*(y(82)-y(81))/(y(396)^2);

%Bin 82-89
for i=82:89
    if i==85
    dydt(i)=-(p.mu+p.deg_ctrA1+p.deg_ctrA2*(y(i+280)+y(i+290))/(p.Jd_CpdR+(y(i+280)+y(i+290))))*y(i) ...%
      -p.b_CtrAPCckAph*(y(i+220)+y(i+230))*y(i)+(p.ub_CtrAPCckAph+p.deg_cckA)*(y(i+240)+y(i+250))...
    +p.phoCtrA*(y(i+260)+y(i+270))-p.ukdephoCtrA*y(i)...
 +p.D_CtrAP/(y(396)^2)*(-y(i)+y(i-1));
    elseif i==86
            dydt(i)=-(p.mu+p.deg_ctrA1+p.deg_ctrA2*(y(i+280)+y(i+290))/(p.Jd_CpdR+(y(i+280)+y(i+290))))*y(i) ...%
      -p.b_CtrAPCckAph*(y(i+220)+y(i+230))*y(i)+(p.ub_CtrAPCckAph+p.deg_cckA)*(y(i+240)+y(i+250))...
    +p.phoCtrA*(y(i+260)+y(i+270))-p.ukdephoCtrA*y(i)...
 +p.D_CtrAP/(y(396)^2)*(y(i+1)-y(i));
    else
            dydt(i)=-(p.mu+p.deg_ctrA1+p.deg_ctrA2*(y(i+280)+y(i+290))/(p.Jd_CpdR+(y(i+280)+y(i+290))))*y(i) ...%
      -p.b_CtrAPCckAph*(y(i+220)+y(i+230))*y(i)+(p.ub_CtrAPCckAph+p.deg_cckA)*(y(i+240)+y(i+250))...
    +p.phoCtrA*(y(i+260)+y(i+270))-p.ukdephoCtrA*y(i)...
 +p.D_CtrAP/(y(396)^2)*(y(i+1)-2*y(i)+y(i-1));
    end
end
%bin 90
 dydt(90)=-(p.mu+p.deg_ctrA1+p.deg_ctrA2*(y(90+280)+y(90+290))/(p.Jd_CpdR+(y(90+280)+y(90+290))))*y(90) ...%
      -p.b_CtrAPCckAph*(y(90+220)+y(90+230))*y(90)+(p.ub_CtrAPCckAph+p.deg_cckA)*(y(90+240)+y(90+250))...
    +p.phoCtrA*(y(90+260)+y(90+270))-p.ukdephoCtrA*y(90)...
    + p.D_CtrAP*(y(89)-y(90))/(y(396)^2);



%% PleCf
m_pleC = (1-mplec)*y(393)+mplec;%S_pleC=y(393)
% DivKP 231-240
%Bin 91
dydt(91)=0*p.syn_pleC*m_pleC -(p.mu+p.deg_pleC)*y(91) ...
    -p.fb_PleC*y(91-80)*y(91)+p.bf_PleC*y(91+10) ...%PodJp 11-20
     -p.b_PleCDivKP*y(91)*y(91+140)...
     +(p.ub_PleCDivKP+p.deg_divK+p.depho_DivK)*y(91+20)...  %DivKP 231-240; PleCf:DivKP 111-120
    + p.D_PleC*(y(92)-y(91))/(y(396)^2);


%Bin 92-99
for i=92:99
    if i==95
    dydt(i)=p.syn_pleC*m_pleC -(p.mu+p.deg_pleC)*y(i) ...
    -p.fb_PleC*y(i-80)*y(i)+p.bf_PleC*y(i+10) ...%PodJp 11-20
     -p.b_PleCDivKP*y(i)*y(i+140)...
     +(p.ub_PleCDivKP+p.deg_divK+p.depho_DivK)*y(i+20)...  %DivKP 231-240; PleCf:DivKP 111-120
 +p.D_PleC/(y(396)^2)*(-y(i)+y(i-1));
    elseif i==96
            dydt(i)=p.syn_pleC*m_pleC -(p.mu+p.deg_pleC)*y(i) ...
    -p.fb_PleC*y(i-80)*y(i)+p.bf_PleC*y(i+10) ...%PodJp 11-20
     -p.b_PleCDivKP*y(i)*y(i+140)...
     +(p.ub_PleCDivKP+p.deg_divK+p.depho_DivK)*y(i+20)...  %DivKP 231-240; PleCf:DivKP 111-120
 +p.D_PleC/(y(396)^2)*(y(i+1)-y(i));
    else
            dydt(i)=p.syn_pleC*m_pleC -(p.mu+p.deg_pleC)*y(i) ...
    -p.fb_PleC*y(i-80)*y(i)+p.bf_PleC*y(i+10) ...%PodJp 11-20
     -p.b_PleCDivKP*y(i)*y(i+140)...
     +(p.ub_PleCDivKP+p.deg_divK+p.depho_DivK)*y(i+20)...  %DivKP 231-240; PleCf:DivKP 111-120
 +p.D_PleC/(y(396)^2)*(y(i+1)-2*y(i)+y(i-1));
    end
end

%Bin 100
dydt(100)=0*p.syn_pleC*m_pleC -(p.mu+p.deg_pleC)*y(100) ...
    -p.fb_PleC*y(100-80)*y(100)+p.bf_PleC*y(100+10) ...%PodJp 11-20
     -p.b_PleCDivKP*y(100)*y(100+140)...
     +(p.ub_PleCDivKP+p.deg_divK+p.depho_DivK)*y(100+20)...  %DivKP 231-240; PleCf:DivKP 111-120
    + p.D_PleC*(y(99)-y(100))/(y(396)^2);




%% PleCb
for i=101:110
    dydt(i)=-(p.mu+p.deg_pleC)*y(i) ...
        +p.fb_PleC*y(i-90)*y(i-10)-p.bf_PleC*y(i)...%PodJp 11-20
        -p.b_PleCDivKP*y(i)*y(i+130)+(p.ub_PleCDivKP+p.deg_divK+p.depho_DivK)*y(i+20);%DivKP 231-240; PleCb:DivKP 121-130
end


%% PleCf:DivKP
% Bin 111
dydt(111)=-(p.mu+p.deg_pleC+p.deg_divK+p.ub_PleCDivKP+p.depho_DivK+p.ph2kin_PleC1)*y(111)...
    +p.b_PleCDivKP*y(111+120)*y(111-20)...%free DivKP231-240; PleCf 91-100
    -p.fb_PleC*y(111-100)*y(111)+p.bf_PleC*y(111+10)...
    +p.D_PleCDivK*(y(112)-y(111))/(y(396)^2);

% Bin 112-119
for i=112:119
    if i==115
    dydt(i)=-(p.mu+p.deg_pleC+p.deg_divK+p.ub_PleCDivKP+p.depho_DivK+p.ph2kin_PleC1)*y(i)...
    +p.b_PleCDivKP*y(i+120)*y(i-20)...%free DivKP231-240; PleCf 91-100
    -p.fb_PleC*y(i-100)*y(i)+p.bf_PleC*y(i+10)...
    +p.D_PleCDivK*(-y(i)+y(i-1))/(y(396)^2);
    elseif i==116
            dydt(i)=-(p.mu+p.deg_pleC+p.deg_divK+p.ub_PleCDivKP+p.depho_DivK+p.ph2kin_PleC1)*y(i)...
    +p.b_PleCDivKP*y(i+120)*y(i-20)...%free DivKP231-240; PleCf 91-100
    -p.fb_PleC*y(i-100)*y(i)+p.bf_PleC*y(i+10)...
    +p.D_PleCDivK*(y(i+1)-y(i))/(y(396)^2);
    else
            dydt(i)=-(p.mu+p.deg_pleC+p.deg_divK+p.ub_PleCDivKP+p.depho_DivK+p.ph2kin_PleC1)*y(i)...
    +p.b_PleCDivKP*y(i+120)*y(i-20)...%free DivKP231-240; PleCf 91-100
    -p.fb_PleC*y(i-100)*y(i)+p.bf_PleC*y(i+10)...
    +p.D_PleCDivK*(y(i+1)-2*y(i)+y(i-1))/(y(396)^2);
    end
end
% Bin 120
dydt(120)=-(p.mu+p.deg_pleC+p.deg_divK+p.ub_PleCDivKP+p.depho_DivK+p.ph2kin_PleC1)*y(120)...
    +p.b_PleCDivKP*y(120+120)*y(120-20)...%free DivKP231-240; PleCf 91-100
    -p.fb_PleC*y(120-100)*y(120)+p.bf_PleC*y(120+10)...
    +p.D_PleCDivK*(y(119)-y(120))/(y(396)^2);



%% PleCb:DivKP
% Bin 121-130
for i=121:130
    dydt(i)=-(p.mu+p.deg_pleC+p.deg_divK+p.ub_PleCDivKP+p.depho_DivK+p.ph2kin_PleC2)*y(i)...
    +p.b_PleCDivKP*y(i+110)*y(i-20)...%free DivKP231-240; PleCb 101-110
    +p.fb_PleC*y(i-110)*y(i-10)-p.bf_PleC*y(i);
end

%% PleCf,kin
% Bin 131
dydt(131)=-(p.mu+p.deg_pleC)*y(131)-p.fb_PleC*y(131-120)*y(131)+p.bf_PleC*y(131+10)...
    +p.ph2kin_PleC1*y(131-20)...
    +p.D_PleC*(y(132)-y(131))/(y(396)^2);%PleCf:DivKP 111-120

% Bin 132-139
for i=132:139
    if i==135
    dydt(i)=-(p.mu+p.deg_pleC)*y(i)-p.fb_PleC*y(i-120)*y(i)+p.bf_PleC*y(i+10)...
    +p.ph2kin_PleC1*y(i-20)...
    +p.D_PleC*(-y(i)+y(i-1))/(y(396)^2);%PleCf:DivKP 111-120
    elseif i==136
        dydt(i)=-(p.mu+p.deg_pleC)*y(i)-p.fb_PleC*y(i-120)*y(i)+p.bf_PleC*y(i+10)...
    +p.ph2kin_PleC1*y(i-20)...
    +p.D_PleC*(y(i+1)-y(i))/(y(396)^2);%PleCf:DivKP 111-120
    else
        dydt(i)=-(p.mu+p.deg_pleC)*y(i)-p.fb_PleC*y(i-120)*y(i)+p.bf_PleC*y(i+10)...
    +p.ph2kin_PleC1*y(i-20)...
    +p.D_PleC*(y(i+1)-2*y(i)+y(i-1))/(y(396)^2);%PleCf:DivKP 111-120
    end
end
% Bin 140
dydt(140)=-(p.mu+p.deg_pleC)*y(140)-p.fb_PleC*y(140-120)*y(140)+p.bf_PleC*y(140+10)...
    +p.ph2kin_PleC1*y(140-20)...
    +p.D_PleC*(y(139)-y(140))/(y(396)^2);%PleCf:DivKP 111-120


%% PleCb,kin
for i=141:150
    dydt(i)=-(p.mu+p.deg_pleC)*y(i)+p.fb_PleC*y(i-130)*y(i-10)-p.bf_PleC*y(i)...
        +p.ph2kin_PleC2*y(i-20);%PleCb:DivKP 121-130
end

%% DivJf
%Bin 151
dydt(151)=0*p.syn_divJ -(p.mu+p.deg_divJ)*y(151) ...
    -p.fb_DivJ*y(151-110)*y(151)+p.bf_DivJ*y(151+10) ...%SpmXp 41-50
    -p.b_DivJDivKP*y(151)*y(151+80)...%free DivKP231-240;
    +(p.ub_DivJDivKP+p.deg_divK)*y(151+40)...%DivJf:DivKP 191-200
    -p.b_DivJDivK*y(151)*y(151+70)...%free DivK 221-230;
    +(p.ub_DivJDivK+p.deg_divK)*y(151+20)...%DivJf:DivK 171-180
    + p.D_DivJ*(y(152)-y(151))/(y(396)^2);



%Bin 152-159
for i=152:159
    if i==155
    dydt(i)=p.syn_divJ -(p.mu+p.deg_divJ)*y(i) ...
    -p.fb_DivJ*y(i-110)*y(i)+p.bf_DivJ*y(i+10) ...%SpmXp 41-50
    -p.b_DivJDivKP*y(i)*y(i+80)...%free DivKP231-240;
    +(p.ub_DivJDivKP+p.deg_divK)*y(i+40)...%DivJf:DivKP 191-200
    -p.b_DivJDivK*y(i)*y(i+70)...%free DivK 221-230;
    +(p.ub_DivJDivK+p.deg_divK)*y(i+20)...%DivJf:DivK 171-180
    + p.D_DivJ*(-y(i)+y(i-1))/(y(396)^2);
    elseif i==156
            dydt(i)=p.syn_divJ -(p.mu+p.deg_divJ)*y(i) ...
    -p.fb_DivJ*y(i-110)*y(i)+p.bf_DivJ*y(i+10) ...%SpmXp 41-50
    -p.b_DivJDivKP*y(i)*y(i+80)...%free DivKP231-240;
    +(p.ub_DivJDivKP+p.deg_divK)*y(i+40)...%DivJf:DivKP 191-200
    -p.b_DivJDivK*y(i)*y(i+70)...%free DivK 221-230;
    +(p.ub_DivJDivK+p.deg_divK)*y(i+20)...%DivJf:DivK 171-180
    + p.D_DivJ*(y(i+1)-y(i))/(y(396)^2);
    else
            dydt(i)=p.syn_divJ -(p.mu+p.deg_divJ)*y(i) ...
    -p.fb_DivJ*y(i-110)*y(i)+p.bf_DivJ*y(i+10) ...%SpmXp 41-50
    -p.b_DivJDivKP*y(i)*y(i+80)...%free DivKP231-240;
    +(p.ub_DivJDivKP+p.deg_divK)*y(i+40)...%DivJf:DivKP 191-200
    -p.b_DivJDivK*y(i)*y(i+70)...%free DivK 221-230;
    +(p.ub_DivJDivK+p.deg_divK)*y(i+20)...%DivJf:DivK 171-180
    + p.D_DivJ*(y(i+1)-2*y(i)+y(i-1))/(y(396)^2);
    end
end

%Bin 160
dydt(160)=0*p.syn_divJ -(p.mu+p.deg_divJ)*y(160) ...
    -p.fb_DivJ*y(160-110)*y(160)+p.bf_DivJ*y(160+10) ...%SpmXp 41-50
    -p.b_DivJDivKP*y(160)*y(160+80)...%free DivKP231-240;
    +(p.ub_DivJDivKP+p.deg_divK)*y(160+40)...%DivJf:DivKP 191-200
    -p.b_DivJDivK*y(160)*y(160+70)...%free DivK 221-230;
    +(p.ub_DivJDivK+p.deg_divK)*y(160+20)...%DivJf:DivK 171-180
    + p.D_DivJ*(y(159)-y(160))/(y(396)^2);




%% DivJb
for i=161:170
dydt(i)=-(p.mu+p.deg_divJ)*y(i) ...
        +p.fb_DivJ*y(i-120)*y(i-10)-p.bf_DivJ*y(i)...%%SpmXp 41-50 
        -p.b_DivJDivK*y(i+60)*y(i) -p.b_DivJDivKP*y(i+70)*y(i)...%free DivK 221-230;free DivKP231-240;
       +(p.ub_DivJDivK+p.deg_divK)*y(i+20)+(p.ub_DivJDivKP+p.deg_divK)*y(i+40);% DivJb:DivKP 201-210; DivJb:DivK 181-190
end

%% DivJf:DivK
% Bin 171
dydt(171)=-(p.mu+p.deg_divJ+p.deg_divK+p.ub_DivJDivK+p.pho_DivKDivJf)*y(171)...
    +p.b_DivJDivK*y(171-20)*y(171+50)...%free DivK 221-230
    -p.fb_DivJ*y(171-130)*y(171)+p.bf_DivJ*y(171+10)...%SpmXp 41-50 
    + p.D_DivJDivKP*(y(172)-y(171))/(y(396)^2);


% Bin 172-179
for i=172:179
    if i==175
    dydt(i)=-(p.mu+p.deg_divJ+p.deg_divK+p.ub_DivJDivK+p.pho_DivKDivJf)*y(i)...
    +p.b_DivJDivK*y(i-20)*y(i+50)...%free DivK 221-230
    -p.fb_DivJ*y(i-130)*y(i)+p.bf_DivJ*y(i+10)...%SpmXp 41-50 
    + p.D_DivJDivKP*(-y(i)+y(i-1))/(y(396)^2);
    elseif i==176
        dydt(i)=-(p.mu+p.deg_divJ+p.deg_divK+p.ub_DivJDivK+p.pho_DivKDivJf)*y(i)...
    +p.b_DivJDivK*y(i-20)*y(i+50)...%free DivK 221-230
    -p.fb_DivJ*y(i-130)*y(i)+p.bf_DivJ*y(i+10)...%SpmXp 41-50 
    + p.D_DivJDivKP*(y(i+1)-y(i))/(y(396)^2);
    else
        dydt(i)=-(p.mu+p.deg_divJ+p.deg_divK+p.ub_DivJDivK+p.pho_DivKDivJf)*y(i)...
    +p.b_DivJDivK*y(i-20)*y(i+50)...%free DivK 221-230
    -p.fb_DivJ*y(i-130)*y(i)+p.bf_DivJ*y(i+10)...%SpmXp 41-50 
    + p.D_DivJDivKP*(y(i+1)-2*y(i)+y(i-1))/(y(396)^2);
    end
end
%Bin 180
dydt(180)=-(p.mu+p.deg_divJ+p.deg_divK+p.ub_DivJDivK+p.pho_DivKDivJf)*y(180)...
    +p.b_DivJDivK*y(180-20)*y(180+50)...%free DivK 221-230
    -p.fb_DivJ*y(180-130)*y(180)+p.bf_DivJ*y(180+10)...%SpmXp 41-50 
    + p.D_DivJDivKP*(y(179)-y(180))/(y(396)^2);



%% DivJb:DivK
% Bin 181-190
for i=181:190
dydt(i)=-(p.mu+p.deg_divJ+p.deg_divK+p.ub_DivJDivK+p.pho_DivKDivJb)*y(i)...
    +p.b_DivJDivK*y(i-20)*y(i+40)...%free DivK 221-230 DivJb 161-170
    +p.fb_DivJ*y(i-140)*y(i-10)-p.bf_DivJ*y(i);%SpmXp 41-50
end


%% DivJf:DivKP
% Bin 191
dydt(191)=-(p.mu+p.deg_divJ+p.deg_divK+p.ub_DivJDivKP)*y(191)...%DivJf:DivK 171-180
    +p.pho_DivKDivJf*y(191-20)+p.b_DivJDivKP*y(191-40)*y(191+40)...free DivKP 231-240; DivJf 151-160
    -p.fb_DivJ*y(191-150)*y(191)+p.bf_DivJ*y(191+10)...%SpmXp 41-50
    + p.D_DivJDivKP*(y(192)-y(191))/(y(396)^2);


% Bin 192-199
for i=192:199
    if i==195
    dydt(i)=-(p.mu+p.deg_divJ+p.deg_divK+p.ub_DivJDivKP)*y(i)...%DivJf:DivK 171-180
    +p.pho_DivKDivJf*y(i-20)+p.b_DivJDivKP*y(i-40)*y(i+40)...free DivKP 231-240; DivJf 151-160
    -p.fb_DivJ*y(i-150)*y(i)+p.bf_DivJ*y(i+10)...%SpmXp 41-50
    + p.D_DivJDivKP*(-y(i)+y(i-1))/(y(396)^2);
    elseif i==196
        dydt(i)=-(p.mu+p.deg_divJ+p.deg_divK+p.ub_DivJDivKP)*y(i)...%DivJf:DivK 171-180
    +p.pho_DivKDivJf*y(i-20)+p.b_DivJDivKP*y(i-40)*y(i+40)...free DivKP 231-240; DivJf 151-160
    -p.fb_DivJ*y(i-150)*y(i)+p.bf_DivJ*y(i+10)...%SpmXp 41-50
    + p.D_DivJDivKP*(y(i+1)-y(i))/(y(396)^2);
    else
        dydt(i)=-(p.mu+p.deg_divJ+p.deg_divK+p.ub_DivJDivKP)*y(i)...%DivJf:DivK 171-180
    +p.pho_DivKDivJf*y(i-20)+p.b_DivJDivKP*y(i-40)*y(i+40)...free DivKP 231-240; DivJf 151-160
    -p.fb_DivJ*y(i-150)*y(i)+p.bf_DivJ*y(i+10)...%SpmXp 41-50
    + p.D_DivJDivKP*(y(i+1)-2*y(i)+y(i-1))/(y(396)^2);
    end
end

%Bin200
dydt(200)=-(p.mu+p.deg_divJ+p.deg_divK+p.ub_DivJDivKP)*y(200)...%DivJf:DivK 171-180
    +p.pho_DivKDivJf*y(200-20)+p.b_DivJDivKP*y(200-40)*y(200+40)...free DivKP 231-240; DivJf 151-160
    -p.fb_DivJ*y(200-150)*y(200)+p.bf_DivJ*y(200+10)...%SpmXp 41-50
    + p.D_DivJDivKP*(y(199)-y(200))/(y(396)^2);




%% DivJb:DivKP

for i=201:210
    dydt(i)=-(p.mu+p.deg_divJ+p.deg_divK+p.ub_DivJDivKP)*y(i)...%DivJb:DivK 181-190
          +p.pho_DivKDivJb*y(i-20)+p.b_DivJDivKP*y(i-40)*y(i+30)...%DivJb 161-170; DivKP 231-240
            +p.fb_DivJ*y(i-160)*y(i-10)-p.bf_DivJ*y(i);%SpmXp 41-50
end

%% PerP
m_perP=((1-m)*y(394)+m);
n_perp=2;
%Bin 211
dydt(211)=0*p.syn_perP*m_perP*y(211-130)^n_perp/(y(211-130)^n_perp+p.Ja_PerPCtrA^n_perp) ...%CtrAP81-90
    -(p.mu+p.deg_perP)*y(211) ...
    +p.D_PerP*(y(212)-y(211))/(y(396)^2);

%Bin 212-219
for i=212:219
    if i==215
    dydt(i)=p.syn_perP*m_perP*y(i-130)^n_perp/(y(i-130)^n_perp+p.Ja_PerPCtrA^n_perp) ...%CtrAP81-90
    -(p.mu+p.deg_perP)*y(i) ...
    +p.D_PerP*(-y(i)+y(i-1))/(y(396)^2);
    elseif i==216
            dydt(i)=p.syn_perP*m_perP*y(i-130)^n_perp/(y(i-130)^n_perp+p.Ja_PerPCtrA^n_perp) ...%CtrAP81-90
    -(p.mu+p.deg_perP)*y(i) ...
    +p.D_PerP*(y(i+1)-y(i))/(y(396)^2);
    else
            dydt(i)=p.syn_perP*m_perP*y(i-130)^2/(y(i-130)^2+p.Ja_PerPCtrA^2) ...%CtrAP81-90
    -(p.mu+p.deg_perP)*y(i) ...
    +p.D_PerP*(y(i+1)-2*y(i)+y(i-1))/(y(396)^2);
    end
end

%Bin220
dydt(220)=0*p.syn_perP*m_perP*y(220-130)^2/(y(220-130)^2+p.Ja_PerPCtrA^2) ...%CtrAP81-90
    -(p.mu+p.deg_perP)*y(220) ...
    +p.D_PerP*(y(219)-y(220))/(y(396)^2);


%% DivK
%DivJb 161-170 DivJb:DivK 181-190  DivJb:DivKP 201-210
%PleCf 91-100 PleCb 101-110
%DivJf 151-160 DivJf:DivK 171-180
 %Bin 221
dydt(221)=0*p.syn_divK1+0*p.syn_divK2*y(221-140)^2/(y(221-140)^2+p.Ja_DivKCtrA^2) ...%CtrAP 81-90
    -(p.mu+p.deg_divK)*y(221) ...
    -p.pho_DivKPleC*(y(221-90)+y(221-80))*y(221)...%PleCfkin 131-140; PleCbkin 141-150
    +p.depho_DivK*(y(221-110)+y(221-100)) ...%PleCf:DivKP 111-120; PleCb:DivKP 121-130
    -p.b_DivJDivK*(y(221-70)+y(221-60))*y(221)+(p.ub_DivJDivK+p.deg_divJ)*(y(221-50)+y(221-40))...%
+p.D_DivK*(y(222)-y(221))/(y(396)^2);


%Bin 222-229
for i=222:229
    if i==225
    dydt(i)=p.syn_divK1+p.syn_divK2*y(i-140)^2/(y(i-140)^2+p.Ja_DivKCtrA^2) ...%CtrAP 81-90
    -(p.mu+p.deg_divK)*y(i) ...
    -p.pho_DivKPleC*(y(i-90)+y(i-80))*y(i)...%PleCfkin 131-140; PleCbkin 141-150
    +p.depho_DivK*(y(i-110)+y(i-100)) ...%PleCf:DivKP 111-120; PleCb:DivKP 121-130
    -p.b_DivJDivK*(y(i-70)+y(i-60))*y(i)+(p.ub_DivJDivK+p.deg_divJ)*(y(i-50)+y(i-40))...%
+p.D_DivK*(-y(i)+y(i-1))/(y(396)^2);
    elseif i==226
            dydt(i)=p.syn_divK1+p.syn_divK2*y(i-140)^2/(y(i-140)^2+p.Ja_DivKCtrA^2) ...%CtrAP 81-90
    -(p.mu+p.deg_divK)*y(i) ...
    -p.pho_DivKPleC*(y(i-90)+y(i-80))*y(i)...%PleCfkin 131-140; PleCbkin 141-150
    +p.depho_DivK*(y(i-110)+y(i-100)) ...%PleCf:DivKP 111-120; PleCb:DivKP 121-130
    -p.b_DivJDivK*(y(i-70)+y(i-60))*y(i)+(p.ub_DivJDivK+p.deg_divJ)*(y(i-50)+y(i-40))...%
+p.D_DivK*(y(i+1)-y(i))/(y(396)^2);
    else
            dydt(i)=p.syn_divK1+p.syn_divK2*y(i-140)^2/(y(i-140)^2+p.Ja_DivKCtrA^2) ...%CtrAP 81-90
    -(p.mu+p.deg_divK)*y(i) ...
    -p.pho_DivKPleC*(y(i-90)+y(i-80))*y(i)...%PleCfkin 131-140; PleCbkin 141-150
    +p.depho_DivK*(y(i-110)+y(i-100)) ...%PleCf:DivKP 111-120; PleCb:DivKP 121-130
    -p.b_DivJDivK*(y(i-70)+y(i-60))*y(i)+(p.ub_DivJDivK+p.deg_divJ)*(y(i-50)+y(i-40))...%
+p.D_DivK*(y(i+1)-2*y(i)+y(i-1))/(y(396)^2);
    end
end

%Bin230
dydt(230)=0*p.syn_divK1+0*p.syn_divK2*y(230-140)^2/(y(230-140)^2+p.Ja_DivKCtrA^2) ...%CtrAP 81-90
    -(p.mu+p.deg_divK)*y(230) ...
    -p.pho_DivKPleC*(y(230-90)+y(230-80))*y(230)...%PleCfkin 131-140; PleCbkin 141-150
    +p.depho_DivK*(y(230-110)+y(230-100)) ...%PleCf:DivKP 111-120; PleCb:DivKP 121-130
    -p.b_DivJDivK*(y(230-70)+y(230-60))*y(230)+(p.ub_DivJDivK+p.deg_divJ)*(y(230-50)+y(230-40))...%
+p.D_DivK*(y(229)-y(230))/(y(396)^2);


%% DivKP
%DivJb 161-170 DivJb:DivK 181-190  DivJb:DivKP 201-210
%PleCf 91-100 PleCb 101-110
%DivJf 151-160 DivJf:DivK 171-180 DivJf:DivKP 191-200


%DivLf 241-250; DivLb 251-260
%DivLf:DivKP 261-270; DivLb:DivKP 271-280
%Bin 231
dydt(231)= -(p.mu+p.deg_divK)*y(231) ...
     +p.pho_DivKPleC*(y(231-90)+y(231-100))*y(231-10)...%PleCfkin 131-140 PleCbkin 141-150
-p.b_PleCDivKP*y(231)*(y(231-140)+y(231-130))+(p.ub_PleCDivKP+p.deg_pleC)*(y(231-120)+y(231-110))...%PleCfDivKP111-120 PleCbDivKP121-130
    -p.b_DivJDivKP*y(231)*(y(231-70)+y(231-80))+(p.ub_DivJDivKP+p.deg_divJ)*(y(231-40)+y(231-30))...
    -p.b_DivLDivKP*y(231)*(y(231+10)+y(231+20))+(p.ub_DivLDivKP+p.deg_divL)*(y(231+30)+y(231+40))...
         +p.ph2kin_PleC1*y(231-120)+p.ph2kin_PleC2*y(231-110)...
+p.D_DivKP*(y(232)-y(231))/(y(396)^2);


%Bin 232-239
for  i=232:239
    if i==235
    dydt(i)= -(p.mu+p.deg_divK)*y(i) ...
     +p.pho_DivKPleC*(y(i-90)+y(i-100))*y(i-10)...%PleCfkin 131-140 PleCbkin 141-150
-p.b_PleCDivKP*y(i)*(y(i-140)+y(i-130))+(p.ub_PleCDivKP+p.deg_pleC)*(y(i-120)+y(i-110))...%PleCfDivKP111-120 PleCbDivKP121-130
    -p.b_DivJDivKP*y(i)*(y(i-70)+y(i-80))+(p.ub_DivJDivKP+p.deg_divJ)*(y(i-40)+y(i-30))...
    -p.b_DivLDivKP*y(i)*(y(i+10)+y(i+20))+(p.ub_DivLDivKP+p.deg_divL)*(y(i+30)+y(i+40))...
         +p.ph2kin_PleC1*y(i-120)+p.ph2kin_PleC2*y(i-110)...
+p.D_DivKP*(-y(i)+y(i-1))/(y(396)^2);
    elseif i==236
        dydt(i)= -(p.mu+p.deg_divK)*y(i) ...
     +p.pho_DivKPleC*(y(i-90)+y(i-100))*y(i-10)...%PleCfkin 131-140 PleCbkin 141-150
-p.b_PleCDivKP*y(i)*(y(i-140)+y(i-130))+(p.ub_PleCDivKP+p.deg_pleC)*(y(i-120)+y(i-110))...%PleCfDivKP111-120 PleCbDivKP121-130
    -p.b_DivJDivKP*y(i)*(y(i-70)+y(i-80))+(p.ub_DivJDivKP+p.deg_divJ)*(y(i-40)+y(i-30))...
    -p.b_DivLDivKP*y(i)*(y(i+10)+y(i+20))+(p.ub_DivLDivKP+p.deg_divL)*(y(i+30)+y(i+40))...
         +p.ph2kin_PleC1*y(i-120)+p.ph2kin_PleC2*y(i-110)...
+p.D_DivKP*(y(i+1)-y(i))/(y(396)^2);
    else
        dydt(i)= -(p.mu+p.deg_divK)*y(i) ...
     +p.pho_DivKPleC*(y(i-90)+y(i-100))*y(i-10)...%PleCfkin 131-140 PleCbkin 141-150
-p.b_PleCDivKP*y(i)*(y(i-140)+y(i-130))+(p.ub_PleCDivKP+p.deg_pleC)*(y(i-120)+y(i-110))...%PleCfDivKP111-120 PleCbDivKP121-130
    -p.b_DivJDivKP*y(i)*(y(i-70)+y(i-80))+(p.ub_DivJDivKP+p.deg_divJ)*(y(i-40)+y(i-30))...
    -p.b_DivLDivKP*y(i)*(y(i+10)+y(i+20))+(p.ub_DivLDivKP+p.deg_divL)*(y(i+30)+y(i+40))...
         +p.ph2kin_PleC1*y(i-120)+p.ph2kin_PleC2*y(i-110)...
+p.D_DivKP*(y(i+1)-2*y(i)+y(i-1))/(y(396)^2);
    end
end


%Bin 240
dydt(240)= -(p.mu+p.deg_divK)*y(240) ...
     +p.pho_DivKPleC*(y(240-90)+y(240-100))*y(240-10)...%PleCfkin 131-140 PleCbkin 141-150
-p.b_PleCDivKP*y(240)*(y(240-140)+y(240-130))+(p.ub_PleCDivKP+p.deg_pleC)*(y(240-120)+y(240-110))...%PleCfDivKP111-120 PleCbDivKP121-130
    -p.b_DivJDivKP*y(240)*(y(240-70)+y(240-80))+(p.ub_DivJDivKP+p.deg_divJ)*(y(240-40)+y(240-30))...
    -p.b_DivLDivKP*y(240)*(y(240+10)+y(240+20))+(p.ub_DivLDivKP+p.deg_divL)*(y(240+30)+y(240+40))...
         +p.ph2kin_PleC1*y(240-120)+p.ph2kin_PleC2*y(240-110)...
+p.D_DivKP*(y(239)-y(240))/(y(396)^2);





%% DivLf
%DivKP 231-240 %DivLf:DivKP 261-270
%PopZp: 61-70  PodJp: 11-20 
m_divL=1;%((1-m)*y(395)+m);
%Bin 241
dydt(241)=0*p.syn_divL*m_divL-(p.mu+p.deg_divL)*y(241) ...
    -p.fb_divL*(p.alpha_DivLPopZ*y(241-180)+p.alpha_DivLPodJ*y(241-230))*y(241)+p.bf_divL*y(241+10)...
     -p.b_DivLDivKP*y(241)*y(241-10)+(p.ub_DivLDivKP+p.deg_divK)*y(241+20)...
     +p.D_DivL*(y(242)-y(241))/(y(396)^2);


  %Bin 242-249
  for i=242:249
      if i==245
      dydt(i)=p.syn_divL*m_divL-(p.mu+p.deg_divL)*y(i) ...
    -p.fb_divL*(p.alpha_DivLPopZ*y(i-180)+p.alpha_DivLPodJ*y(i-230))*y(i)+p.bf_divL*y(i+10)...
     -p.b_DivLDivKP*y(i)*y(i-10)+(p.ub_DivLDivKP+p.deg_divK)*y(i+20)...
     +p.D_DivL*(-y(i)+y(i-1))/(y(396)^2);
      elseif i==246
                dydt(i)=p.syn_divL*m_divL-(p.mu+p.deg_divL)*y(i) ...
    -p.fb_divL*(p.alpha_DivLPopZ*y(i-180)+p.alpha_DivLPodJ*y(i-230))*y(i)+p.bf_divL*y(i+10)...
     -p.b_DivLDivKP*y(i)*y(i-10)+(p.ub_DivLDivKP+p.deg_divK)*y(i+20)...
     +p.D_DivL*(y(i+1)-y(i))/(y(396)^2);
      else
                dydt(i)=p.syn_divL*m_divL-(p.mu+p.deg_divL)*y(i) ...
    -p.fb_divL*(p.alpha_DivLPopZ*y(i-180)+p.alpha_DivLPodJ*y(i-230))*y(i)+p.bf_divL*y(i+10)...
     -p.b_DivLDivKP*y(i)*y(i-10)+(p.ub_DivLDivKP+p.deg_divK)*y(i+20)...
     +p.D_DivL*(y(i+1)-2*y(i)+y(i-1))/(y(396)^2);
      end
  end

     %Bin 250
dydt(250)=0*p.syn_divL*m_divL-(p.mu+p.deg_divL)*y(250) ...
    -p.fb_divL*(p.alpha_DivLPopZ*y(250-180)+p.alpha_DivLPodJ*y(250-230))*y(250)+p.bf_divL*y(250+10)...
     -p.b_DivLDivKP*y(250)*y(250-10)+(p.ub_DivLDivKP+p.deg_divK)*y(250+20)...
     +p.D_DivL*(y(249)-y(250))/(y(396)^2);

%% DivLb
%DivKP 231-240 %DivLb:DivKP 271-280
%PopZp: 61-70  PodJp: 11-20 

for i=251:260
dydt(i)=-(p.mu+p.deg_divL)*y(i)+p.fb_divL*(p.alpha_DivLPopZ*y(i-190)+p.alpha_DivLPodJ*y(i-240))*y(i-10)-p.bf_divL*y(i)...
        -p.b_DivLDivKP*y(i)*y(i-20)+(p.ub_DivLDivKP+p.deg_divK)*y(i+20);
end

%% DivLf:DivKP
%DivLf 241-250  DivKP 231-240
%PopZp: 61-70  PodJp: 11-20 
%Bin 261
dydt(261)=-(p.mu+p.deg_divL+p.deg_divK+p.ub_DivLDivKP)*y(261)+p.b_DivLDivKP*y(261-20)*y(261-30)...
    -p.fb_divL*(p.alpha_DivLPopZ*y(261-200)+p.alpha_DivLPodJ*y(261-250))*y(261)+p.bf_divL*y(261+10)...
         +p.D_DivLDivK*(y(262)-y(261))/(y(396)^2);
     %Bin 262-269
     for i=262:269
         if i==265
         dydt(i)=-(p.mu+p.deg_divL+p.deg_divK+p.ub_DivLDivKP)*y(i)+p.b_DivLDivKP*y(i-20)*y(i-30)...
    -p.fb_divL*(p.alpha_DivLPopZ*y(i-200)+p.alpha_DivLPodJ*y(i-250))*y(i)+p.bf_divL*y(i+10)...
         +p.D_DivLDivK*(-y(i)+y(i-1))/(y(396)^2);
         elseif i==266
         dydt(i)=-(p.mu+p.deg_divL+p.deg_divK+p.ub_DivLDivKP)*y(i)+p.b_DivLDivKP*y(i-20)*y(i-30)...
    -p.fb_divL*(p.alpha_DivLPopZ*y(i-200)+p.alpha_DivLPodJ*y(i-250))*y(i)+p.bf_divL*y(i+10)...
         +p.D_DivLDivK*(y(i+1)-y(i))/(y(396)^2);
         else
         dydt(i)=-(p.mu+p.deg_divL+p.deg_divK+p.ub_DivLDivKP)*y(i)+p.b_DivLDivKP*y(i-20)*y(i-30)...
    -p.fb_divL*(p.alpha_DivLPopZ*y(i-200)+p.alpha_DivLPodJ*y(i-250))*y(i)+p.bf_divL*y(i+10)...
         +p.D_DivLDivK*(y(i+1)-2*y(i)+y(i-1))/(y(396)^2);
         end
     end

     %Bin 270
dydt(270)=-(p.mu+p.deg_divL+p.deg_divK+p.ub_DivLDivKP)*y(270)+p.b_DivLDivKP*y(270-20)*y(270-30)...
    -p.fb_divL*(p.alpha_DivLPopZ*y(270-200)+p.alpha_DivLPodJ*y(270-250))*y(270)+p.bf_divL*y(270+10)...
         +p.D_DivLDivK*(y(269)-y(270))/(y(396)^2);

%% DivLb:DivKP
for i=271:280
    dydt(i)=-(p.mu+p.deg_divL+p.deg_divK+p.ub_DivLDivKP)*y(i)+p.b_DivLDivKP*y(i-20)*y(i-40)...
    +p.fb_divL*(p.alpha_DivLPopZ*y(i-210)+p.alpha_DivLPodJ*y(i-260))*y(i-10)+p.bf_divL*y(i);
end

%% CckAfkin
%PopZp: 61-70  PodJp: 11-20 
%DivLb251-260; DivLb:DivKP 271-280
%DivLf 241-250
%CckAfph 301-310
%CckAfkin:CtrA 341-350
%CtrA 71-80 CpdRf 361-370 CpdRb 371-380
%Bin281 
dydt(281)=0*p.syn_cckA-(p.mu+p.deg_cckA)*y(281)...
    -p.fb_cckA*(p.alpha_CckAPopZ*y(281-220)+p.alpha_CckADivL*(y(281-30)+y(281-10))*p.alpha_DivLPodJ*y(281-270)/(p.alpha_DivLPopZ*y(281-220)+p.alpha_DivLPodJ*y(281-270)))*y(281)...
    +p.bf_cckA*y(281+10)...
    +(p.pk_cckA1+p.pk_cckA2*(y(281-30)+y(281-40)))*y(281+20)-(p.kp_cckA1+p.kp_cckA2*(y(281-10)+y(281-20)))*y(281)...
    -p.b_CtrACckAkin*y(281)*y(281-210)...
    +(p.deg_ctrA1+p.deg_ctrA2*(y(281+80)+y(281+90))/(p.Jd_CpdR+(y(281+80)+y(281+90)))+p.ub_CtrACckAkin+p.phoCtrA)*y(281+60)...
    +p.D_CckA*(y(282)-y(281))/(y(396)^2);



%Bin282-289
for i=282:289
    if i==285
    dydt(i)=p.syn_cckA-(p.mu+p.deg_cckA)*y(i)...
    -p.fb_cckA*(p.alpha_CckAPopZ*y(i-220)+p.alpha_CckADivL*(y(i-30)+y(i-10))*p.alpha_DivLPodJ*y(i-270)/(p.alpha_DivLPopZ*y(i-220)+p.alpha_DivLPodJ*y(i-270)))*y(i)...
    +p.bf_cckA*y(i+10)...
    +(p.pk_cckA1+p.pk_cckA2*(y(i-30)+y(i-40)))*y(i+20)-(p.kp_cckA1+p.kp_cckA2*(y(i-10)+y(i-20)))*y(i)...
    -p.b_CtrACckAkin*y(i)*y(i-210)...
    +(p.deg_ctrA1+p.deg_ctrA2*(y(i+80)+y(i+90))/(p.Jd_CpdR+(y(i+80)+y(i+90)))+p.ub_CtrACckAkin+p.phoCtrA)*y(i+60)...
    +p.D_CckA*(-y(i)+y(i-1))/(y(396)^2);
    elseif i==286
        dydt(i)=p.syn_cckA-(p.mu+p.deg_cckA)*y(i)...
    -p.fb_cckA*(p.alpha_CckAPopZ*y(i-220)+p.alpha_CckADivL*(y(i-30)+y(i-10))*p.alpha_DivLPodJ*y(i-270)/(p.alpha_DivLPopZ*y(i-220)+p.alpha_DivLPodJ*y(i-270)))*y(i)...
    +p.bf_cckA*y(i+10)...
    +(p.pk_cckA1+p.pk_cckA2*(y(i-30)+y(i-40)))*y(i+20)-(p.kp_cckA1+p.kp_cckA2*(y(i-10)+y(i-20)))*y(i)...
    -p.b_CtrACckAkin*y(i)*y(i-210)...
    +(p.deg_ctrA1+p.deg_ctrA2*(y(i+80)+y(i+90))/(p.Jd_CpdR+(y(i+80)+y(i+90)))+p.ub_CtrACckAkin+p.phoCtrA)*y(i+60)...
    +p.D_CckA*(y(i+1)-y(i))/(y(396)^2);
    else
        dydt(i)=p.syn_cckA-(p.mu+p.deg_cckA)*y(i)...
    -p.fb_cckA*(p.alpha_CckAPopZ*y(i-220)+p.alpha_CckADivL*(y(i-30)+y(i-10))*p.alpha_DivLPodJ*y(i-270)/(p.alpha_DivLPopZ*y(i-220)+p.alpha_DivLPodJ*y(i-270)))*y(i)...
    +p.bf_cckA*y(i+10)...
    +(p.pk_cckA1+p.pk_cckA2*(y(i-30)+y(i-40)))*y(i+20)-(p.kp_cckA1+p.kp_cckA2*(y(i-10)+y(i-20)))*y(i)...
    -p.b_CtrACckAkin*y(i)*y(i-210)...
    +(p.deg_ctrA1+p.deg_ctrA2*(y(i+80)+y(i+90))/(p.Jd_CpdR+(y(i+80)+y(i+90)))+p.ub_CtrACckAkin+p.phoCtrA)*y(i+60)...
    +p.D_CckA*(y(i+1)-2*y(i)+y(i-1))/(y(396)^2);
    end
end


%Bin290
dydt(290)=0*p.syn_cckA-(p.mu+p.deg_cckA)*y(290)...
    -p.fb_cckA*(p.alpha_CckAPopZ*y(290-220)+p.alpha_CckADivL*(y(290-30)+y(290-10))*p.alpha_DivLPodJ*y(290-270)/(p.alpha_DivLPopZ*y(290-220)+p.alpha_DivLPodJ*y(290-270)))*y(290)...
    +p.bf_cckA*y(290+10)...
    +(p.pk_cckA1+p.pk_cckA2*(y(290-30)+y(290-40)))*y(290+20)-(p.kp_cckA1+p.kp_cckA2*(y(290-10)+y(290-20)))*y(290)...
    -p.b_CtrACckAkin*y(290)*y(290-210)...
    +(p.deg_ctrA1+p.deg_ctrA2*(y(290+80)+y(290+90))/(p.Jd_CpdR+(y(290+80)+y(290+90)))+p.ub_CtrACckAkin+p.phoCtrA)*y(290+60)...
    +p.D_CckA*(y(289)-y(290))/(y(396)^2);


%% CckAbkin
%PopZp: 61-70  PodJp: 11-20 
%DivLb251-260; DivLb:DivKP 271-280
%DivLf 241-250
%CckAbph 311-320
%CckAbkin:CtrA 351-360
%CtrA 71-80 CpdRf 361-370 CpdRb 371-380



for i=291:300
    dydt(i)=-(p.mu+p.deg_cckA)*y(i)-p.bf_cckA*y(i)...
        +p.fb_cckA*(p.alpha_CckAPopZ*y(i-230)+p.alpha_CckADivL*(y(i-40)+y(i-20))*p.alpha_DivLPodJ*y(i-280)/(p.alpha_DivLPopZ*y(i-230)+p.alpha_DivLPodJ*y(i-280)))*y(i-10)...
            +(p.pk_cckA1+p.pk_cckA2*(y(i-40)+y(i-50)))*y(i+20)-(p.kp_cckA1+p.kp_cckA2*(y(i-20)+y(i-30)))*y(i)...
            -p.b_CtrACckAkin*y(i)*y(i-220)...
            +(p.deg_ctrA1+p.deg_ctrA2*(y(i+70)+y(i+80))/(p.Jd_CpdR+(y(i+70)+y(i+80)))+p.ub_CtrACckAkin+p.phoCtrA)*y(i+60);
end


%% CckAfph
%PopZp: 61-70  PodJp: 11-20 
%DivLb251-260; DivLb:DivKP 271-280
%DivLf 241-250  CckAfkin281-290
%CtrAP 81-90 CpdRf 361-370 CpdRb 371-380

%CckAfph:CtrAP 301-310
%Bin 301
dydt(301)=-(p.mu+p.deg_cckA)*y(301)+p.bf_cckA*y(301+10)...
    -p.fb_cckA*(p.alpha_CckAPopZ*y(301-240)+p.alpha_CckADivL*(y(301-50)+y(301-30))*p.alpha_DivLPodJ*y(301-290)/(p.alpha_DivLPopZ*y(301-240)+p.alpha_DivLPodJ*y(301-290)))*y(301)...
    -(p.pk_cckA1+p.pk_cckA2*(y(301-50)+y(301-60)))*y(301)+(p.kp_cckA1+p.kp_cckA2*(y(301-30)+y(301-40)))*y(301-20)...
        -p.b_CtrAPCckAph*y(301)*y(301-220)...
        +(p.deg_ctrA1+p.deg_ctrA2*(y(301+60)+y(301+70))/(p.Jd_CpdR+(y(301+60)+y(301+70)))+p.ub_CtrAPCckAph+p.dephoCtrA)*y(301+20)...
        +p.D_CckA*(y(302)-y(301))/(y(396)^2);

%Bin 302-309
for i=302:309
    if i==305
    dydt(i)=-(p.mu+p.deg_cckA)*y(i)+p.bf_cckA*y(i+10)...
    -p.fb_cckA*(p.alpha_CckAPopZ*y(i-240)+p.alpha_CckADivL*(y(i-50)+y(i-30))*p.alpha_DivLPodJ*y(i-290)/(p.alpha_DivLPopZ*y(i-240)+p.alpha_DivLPodJ*y(i-290)))*y(i)...
    -(p.pk_cckA1+p.pk_cckA2*(y(i-50)+y(i-60)))*y(i)+(p.kp_cckA1+p.kp_cckA2*(y(i-30)+y(i-40)))*y(i-20)...
        -p.b_CtrAPCckAph*y(i)*y(i-220)...
        +(p.deg_ctrA1+p.deg_ctrA2*(y(i+60)+y(i+70))/(p.Jd_CpdR+(y(i+60)+y(i+70)))+p.ub_CtrAPCckAph+p.dephoCtrA)*y(i+20)...
        +p.D_CckA*(-y(i)+y(i-1))/(y(396)^2);
    elseif i==306
        dydt(i)=-(p.mu+p.deg_cckA)*y(i)+p.bf_cckA*y(i+10)...
    -p.fb_cckA*(p.alpha_CckAPopZ*y(i-240)+p.alpha_CckADivL*(y(i-50)+y(i-30))*p.alpha_DivLPodJ*y(i-290)/(p.alpha_DivLPopZ*y(i-240)+p.alpha_DivLPodJ*y(i-290)))*y(i)...
    -(p.pk_cckA1+p.pk_cckA2*(y(i-50)+y(i-60)))*y(i)+(p.kp_cckA1+p.kp_cckA2*(y(i-30)+y(i-40)))*y(i-20)...
        -p.b_CtrAPCckAph*y(i)*y(i-220)...
        +(p.deg_ctrA1+p.deg_ctrA2*(y(i+60)+y(i+70))/(p.Jd_CpdR+(y(i+60)+y(i+70)))+p.ub_CtrAPCckAph+p.dephoCtrA)*y(i+20)...
        +p.D_CckA*(y(i+1)-y(i))/(y(396)^2);
    else
        dydt(i)=-(p.mu+p.deg_cckA)*y(i)+p.bf_cckA*y(i+10)...
    -p.fb_cckA*(p.alpha_CckAPopZ*y(i-240)+p.alpha_CckADivL*(y(i-50)+y(i-30))*p.alpha_DivLPodJ*y(i-290)/(p.alpha_DivLPopZ*y(i-240)+p.alpha_DivLPodJ*y(i-290)))*y(i)...
    -(p.pk_cckA1+p.pk_cckA2*(y(i-50)+y(i-60)))*y(i)+(p.kp_cckA1+p.kp_cckA2*(y(i-30)+y(i-40)))*y(i-20)...
        -p.b_CtrAPCckAph*y(i)*y(i-220)...
        +(p.deg_ctrA1+p.deg_ctrA2*(y(i+60)+y(i+70))/(p.Jd_CpdR+(y(i+60)+y(i+70)))+p.ub_CtrAPCckAph+p.dephoCtrA)*y(i+20)...
        +p.D_CckA*(y(i+1)-2*y(i)+y(i-1))/(y(396)^2);
    end
end

%Bin 310
dydt(310)=-(p.mu+p.deg_cckA)*y(310)+p.bf_cckA*y(310+10)...
    -p.fb_cckA*(p.alpha_CckAPopZ*y(310-240)+p.alpha_CckADivL*(y(310-50)+y(310-30))*p.alpha_DivLPodJ*y(310-290)/(p.alpha_DivLPopZ*y(310-240)+p.alpha_DivLPodJ*y(310-290)))*y(310)...
    -(p.pk_cckA1+p.pk_cckA2*(y(310-50)+y(310-60)))*y(310)+(p.kp_cckA1+p.kp_cckA2*(y(310-30)+y(310-40)))*y(310-20)...
        -p.b_CtrAPCckAph*y(310)*y(310-220)...
        +(p.deg_ctrA1+p.deg_ctrA2*(y(310+60)+y(310+70))/(p.Jd_CpdR+(y(310+60)+y(310+70)))+p.ub_CtrAPCckAph+p.dephoCtrA)*y(310+20)...
        +p.D_CckA*(y(309)-y(310))/(y(396)^2);

%% CckAbph
%PopZp: 61-70  PodJp: 11-20 
%DivLb251-260; DivLb:DivKP 271-280
%DivLf 241-250  CckAfkin281-290
%CtrAP 81-90 CpdRf 361-370 CpdRb 371-380

%CckAbph:CtrAP 331-340
for i=311:320
dydt(i)=-(p.mu+p.deg_cckA)*y(i)-p.bf_cckA*y(i)...
    +p.fb_cckA*(p.alpha_CckAPopZ*y(i-250)+p.alpha_CckADivL*(y(i-60)+y(i-40))*p.alpha_DivLPodJ*y(i-300)/(p.alpha_DivLPopZ*y(i-250)+p.alpha_DivLPodJ*y(i-300)))*y(i-10)...
            -(p.pk_cckA1+p.pk_cckA2*(y(i-60)+y(i-70)))*y(i)+(p.kp_cckA1+p.kp_cckA2*(y(i-40)+y(i-50)))*y(i-20)...
            -p.b_CtrAPCckAph*y(i)*y(i-230)...
        +(p.deg_ctrA1+p.deg_ctrA2*(y(i+50)+y(i+60))/(p.Jd_CpdR+(y(i+50)+y(i+60)))+p.ub_CtrAPCckAph+p.dephoCtrA)*y(i+20);
end
%% CtrAP:CckAfph
%CtrAP 81-90 CpdRf 361-370 CpdRb 371-380
%DivLb251-260; DivLb:DivKP 271-280
%Bin 321
dydt(321)=-(p.mu+p.deg_cckA+p.deg_ctrA1+p.deg_ctrA2*(y(321+40)+y(321+50))/(p.Jd_CpdR+(y(321+40)+y(321+50)))+p.ub_CtrAPCckAph+p.dephoCtrA)*y(321)...
    +p.b_CtrAPCckAph*y(321-20)*y(321-240)...%CtrAP 81-90 CckAfph 301-310
    -p.fb_cckA*(p.alpha_CckAPopZ*y(321-260)+p.alpha_CckADivL*(y(321-70)+y(321-50))*p.alpha_DivLPodJ*y(321-310)/(p.alpha_DivLPopZ*y(321-260)+p.alpha_DivLPodJ*y(321-310)))*y(321)...
    +p.bf_cckA*y(321+10)...
     +p.D_CckACtrA*(y(322)-y(321))/(y(396)^2);



 %Bin 322-329
 for i=322:329
     if i==325
     dydt(i)=-(p.mu+p.deg_cckA+p.deg_ctrA1+p.deg_ctrA2*(y(i+40)+y(i+50))/(p.Jd_CpdR+(y(i+40)+y(i+50)))+p.ub_CtrAPCckAph+p.dephoCtrA)*y(i)...
    +p.b_CtrAPCckAph*y(i-20)*y(i-240)...%CtrAP 81-90 CckAfph 301-310
    -p.fb_cckA*(p.alpha_CckAPopZ*y(i-260)+p.alpha_CckADivL*(y(i-70)+y(i-50))*p.alpha_DivLPodJ*y(i-310)/(p.alpha_DivLPopZ*y(i-260)+p.alpha_DivLPodJ*y(i-310)))*y(i)...
    +p.bf_cckA*y(i+10)...
     +p.D_CckACtrA*(-y(i)+y(i-1))/(y(396)^2);
     elseif i==326
         dydt(i)=-(p.mu+p.deg_cckA+p.deg_ctrA1+p.deg_ctrA2*(y(i+40)+y(i+50))/(p.Jd_CpdR+(y(i+40)+y(i+50)))+p.ub_CtrAPCckAph+p.dephoCtrA)*y(i)...
    +p.b_CtrAPCckAph*y(i-20)*y(i-240)...%CtrAP 81-90 CckAfph 301-310
    -p.fb_cckA*(p.alpha_CckAPopZ*y(i-260)+p.alpha_CckADivL*(y(i-70)+y(i-50))*p.alpha_DivLPodJ*y(i-310)/(p.alpha_DivLPopZ*y(i-260)+p.alpha_DivLPodJ*y(i-310)))*y(i)...
    +p.bf_cckA*y(i+10)...
     +p.D_CckACtrA*(y(i+1)-y(i))/(y(396)^2);
     else
         dydt(i)=-(p.mu+p.deg_cckA+p.deg_ctrA1+p.deg_ctrA2*(y(i+40)+y(i+50))/(p.Jd_CpdR+(y(i+40)+y(i+50)))+p.ub_CtrAPCckAph+p.dephoCtrA)*y(i)...
    +p.b_CtrAPCckAph*y(i-20)*y(i-240)...%CtrAP 81-90 CckAfph 301-310
    -p.fb_cckA*(p.alpha_CckAPopZ*y(i-260)+p.alpha_CckADivL*(y(i-70)+y(i-50))*p.alpha_DivLPodJ*y(i-310)/(p.alpha_DivLPopZ*y(i-260)+p.alpha_DivLPodJ*y(i-310)))*y(i)...
    +p.bf_cckA*y(i+10)...
     +p.D_CckACtrA*(y(i+1)-2*y(i)+y(i-1))/(y(396)^2);
     end
 end

%Bin 330
dydt(330)=-(p.mu+p.deg_cckA+p.deg_ctrA1+p.deg_ctrA2*(y(330+40)+y(330+50))/(p.Jd_CpdR+(y(330+40)+y(330+50)))+p.ub_CtrAPCckAph+p.dephoCtrA)*y(330)...
    +p.b_CtrAPCckAph*y(330-20)*y(330-240)...%CtrAP 81-90 CckAfph 301-310
    -p.fb_cckA*(p.alpha_CckAPopZ*y(330-260)+p.alpha_CckADivL*(y(330-70)+y(330-50))*p.alpha_DivLPodJ*y(330-310)/(p.alpha_DivLPopZ*y(330-260)+p.alpha_DivLPodJ*y(330-310)))*y(330)...
    +p.bf_cckA*y(330+10)...
     +p.D_CckACtrA*(y(329)-y(330))/(y(396)^2);
 

%% CtrAP:CckAbph
%CtrAP 81-90 CpdRf 361-370 CpdRb 371-380
%DivLb251-260; DivLb:DivKP 271-280


for i=331:340
    dydt(i)=-(p.mu+p.deg_cckA+p.deg_ctrA1+p.deg_ctrA2*(y(i+30)+y(i+40))/(p.Jd_CpdR+(y(i+30)+y(i+40)))+p.ub_CtrAPCckAph+p.dephoCtrA)*y(i)...
    +p.b_CtrAPCckAph*y(i-10-10)*y(i-250)...% CckAbph 311-320
    +p.fb_cckA*(p.alpha_CckAPopZ*y(i-270)+p.alpha_CckADivL*(y(i-80)+y(i-60))*p.alpha_DivLPodJ*y(i-320)/(p.alpha_DivLPopZ*y(i-270)+p.alpha_DivLPodJ*y(i-320)))*y(i-10)...
    -p.bf_cckA*y(i);
end

%% CtrA:CckAfkin
%CpdRf 361-370 CpdRb 371-380
%DivLb251-260; DivLb:DivKP 271-280

%Bin 341
dydt(341)=-(p.mu+p.deg_cckA+p.deg_ctrA1+p.deg_ctrA2*(y(341+20)+y(341+30))/(p.Jd_CpdR+(y(341+20)+y(341+30)))+p.ub_CtrACckAkin+p.phoCtrA)*y(341)...
    +p.b_CtrACckAkin*y(341-60)*y(341-270)...%CtrA 71-80 CckAfkin 281-290
    -p.fb_cckA*(p.alpha_CckAPopZ*y(341-280)+p.alpha_CckADivL*(y(341-90)+y(341-70))*p.alpha_DivLPodJ*y(341-330)/(p.alpha_DivLPopZ*y(341-280)+p.alpha_DivLPodJ*y(341-330)))*y(341)...
    +p.bf_cckA*y(341+10)...
     +p.D_CckACtrA*(y(342)-y(341))/(y(396)^2);

%Bin 342-349
for i=342:349
    if i==345
    dydt(i)=-(p.mu+p.deg_cckA+p.deg_ctrA1+p.deg_ctrA2*(y(i+20)+y(i+30))/(p.Jd_CpdR+(y(i+20)+y(i+30)))+p.ub_CtrACckAkin+p.phoCtrA)*y(i)...
    +p.b_CtrACckAkin*y(i-60)*y(i-270)...%CtrA 71-80 CckAfkin 281-290
    -p.fb_cckA*(p.alpha_CckAPopZ*y(i-280)+p.alpha_CckADivL*(y(i-90)+y(i-70))*p.alpha_DivLPodJ*y(i-330)/(p.alpha_DivLPopZ*y(i-280)+p.alpha_DivLPodJ*y(i-330)))*y(i)...
    +p.bf_cckA*y(i+10)...
     +p.D_CckACtrA*(-y(i)+y(i-1))/(y(396)^2);
    elseif i==346
        dydt(i)=-(p.mu+p.deg_cckA+p.deg_ctrA1+p.deg_ctrA2*(y(i+20)+y(i+30))/(p.Jd_CpdR+(y(i+20)+y(i+30)))+p.ub_CtrACckAkin+p.phoCtrA)*y(i)...
    +p.b_CtrACckAkin*y(i-60)*y(i-270)...%CtrA 71-80 CckAfkin 281-290
    -p.fb_cckA*(p.alpha_CckAPopZ*y(i-280)+p.alpha_CckADivL*(y(i-90)+y(i-70))*p.alpha_DivLPodJ*y(i-330)/(p.alpha_DivLPopZ*y(i-280)+p.alpha_DivLPodJ*y(i-330)))*y(i)...
    +p.bf_cckA*y(i+10)...
     +p.D_CckACtrA*(y(i+1)-y(i))/(y(396)^2);
    else
        dydt(i)=-(p.mu+p.deg_cckA+p.deg_ctrA1+p.deg_ctrA2*(y(i+20)+y(i+30))/(p.Jd_CpdR+(y(i+20)+y(i+30)))+p.ub_CtrACckAkin+p.phoCtrA)*y(i)...
    +p.b_CtrACckAkin*y(i-60)*y(i-270)...%CtrA 71-80 CckAfkin 281-290
    -p.fb_cckA*(p.alpha_CckAPopZ*y(i-280)+p.alpha_CckADivL*(y(i-90)+y(i-70))*p.alpha_DivLPodJ*y(i-330)/(p.alpha_DivLPopZ*y(i-280)+p.alpha_DivLPodJ*y(i-330)))*y(i)...
    +p.bf_cckA*y(i+10)...
     +p.D_CckACtrA*(y(i+1)-2*y(i)+y(i-1))/(y(396)^2);
    end
end

%Bin 350
dydt(350)=-(p.mu+p.deg_cckA+p.deg_ctrA1+p.deg_ctrA2*(y(350+20)+y(350+30))/(p.Jd_CpdR+(y(350+20)+y(350+30)))+p.ub_CtrACckAkin+p.phoCtrA)*y(350)...
    +p.b_CtrACckAkin*y(350-60)*y(350-270)...%CtrA 71-80 CckAfkin 281-290
    -p.fb_cckA*(p.alpha_CckAPopZ*y(350-280)+p.alpha_CckADivL*(y(350-90)+y(350-70))*p.alpha_DivLPodJ*y(350-330)/(p.alpha_DivLPopZ*y(350-280)+p.alpha_DivLPodJ*y(350-330)))*y(350)...
    +p.bf_cckA*y(350+10)...
     +p.D_CckACtrA*(y(349)-y(350))/(y(396)^2);


%% CtrA:CckAbkin
%CpdRf 361-370 CpdRb 371-380
%DivLb251-260; DivLb:DivKP 271-280

for i=351:360
    dydt(i)=-(p.mu+p.deg_cckA+p.deg_ctrA1+p.deg_ctrA2*(y(i+10)+y(i+20))/(p.Jd_CpdR+(y(i+10)+y(i+20)))+p.ub_CtrACckAkin+p.phoCtrA)*y(i)...
    +p.b_CtrACckAkin*y(i-60)*y(i-280)...%CtrA 71-80 CckAbkin 291-300
    +p.fb_cckA*(p.alpha_CckAPopZ*y(i-290)+p.alpha_CckADivL*(y(i-100)+y(i-80))*p.alpha_DivLPodJ*y(i-340)/(p.alpha_DivLPopZ*y(i-290)+p.alpha_DivLPodJ*y(i-340)))*y(i-10)...
    -p.bf_cckA*y(i);
end

%% CpdRf
% sw=0;
%CtrAP 81-90
%CckAfkin 281-290 CckAbkin 291-300
%CckAfph 301-310 CckAbph 311-320

%Bin 361
dydt(361)=0*p.syn_cpdR*y(361-280)^2/(p.JaCpdRCtrA^2+y(361-280)^2)-(p.mu+p.deg_cpdR)*y(361)...
    -p.fb_cpdR*y(361-300)*y(361)+p.bf_cpdR*y(361+10)...%PopZp 61-70
    -p.phoCpdR*(y(361-70)+y(361-80))*y(361)+p.dephoCpdR*(y(361-60)+y(361-50))*y(361+20)+p.ukdephoCpdR*y(361+20)...%
    +p.D_CpdR*(y(362)-y(361))/(y(396)^2);


%Bin 362-369
for i=362:369
    if i==365
    dydt(i)=p.syn_cpdR*y(i-280)^2/(p.JaCpdRCtrA^2+y(i-280)^2)-(p.mu+p.deg_cpdR)*y(i)...
    -p.fb_cpdR*y(i-300)*y(i)+p.bf_cpdR*y(i+10)...%PopZp 61-70
    -p.phoCpdR*(y(i-70)+y(i-80))*y(i)+p.dephoCpdR*(y(i-60)+y(i-50))*y(i+20)+p.ukdephoCpdR*y(i+20)...%
    +p.D_CpdR*(-y(i)+y(i-1))/(y(396)^2);
    elseif i==366
        dydt(i)=p.syn_cpdR*y(i-280)^2/(p.JaCpdRCtrA^2+y(i-280)^2)-(p.mu+p.deg_cpdR)*y(i)...
    -p.fb_cpdR*y(i-300)*y(i)+p.bf_cpdR*y(i+10)...%PopZp 61-70
    -p.phoCpdR*(y(i-70)+y(i-80))*y(i)+p.dephoCpdR*(y(i-60)+y(i-50))*y(i+20)+p.ukdephoCpdR*y(i+20)...%
    +p.D_CpdR*(y(i+1)-y(i))/(y(396)^2);
    else
        dydt(i)=p.syn_cpdR*y(i-280)^2/(p.JaCpdRCtrA^2+y(i-280)^2)-(p.mu+p.deg_cpdR)*y(i)...
    -p.fb_cpdR*y(i-300)*y(i)+p.bf_cpdR*y(i+10)...%PopZp 61-70
    -p.phoCpdR*(y(i-70)+y(i-80))*y(i)+p.dephoCpdR*(y(i-60)+y(i-50))*y(i+20)+p.ukdephoCpdR*y(i+20)...%
    +p.D_CpdR*(y(i+1)-2*y(i)+y(i-1))/(y(396)^2);
    end
end

%Bin 370
dydt(370)=0*p.syn_cpdR*y(370-280)^2/(p.JaCpdRCtrA^2+y(370-280)^2)-(p.mu+p.deg_cpdR)*y(370)...
    -p.fb_cpdR*y(370-300)*y(370)+p.bf_cpdR*y(370+10)...%PopZp 61-70
    -p.phoCpdR*(y(370-70)+y(370-80))*y(370)+p.dephoCpdR*(y(370-60)+y(370-50))*y(370+20)+p.ukdephoCpdR*y(370+20)...%
    +p.D_CpdR*(y(369)-y(370))/(y(396)^2);


%% CpdRb
for i=371:380
dydt(i)=-(p.mu+p.deg_cpdR)*y(i)...
    +p.fb_cpdR*y(i-310)*y(i-10)-p.bf_cpdR*y(i);
end

%% CpdR~P
%CckAfkin 281-290 CckAbkin 291-300
%CckAfph 301-310 CckAbph 311-320
%Bin 381
dydt(381)=-(p.mu+p.deg_cpdR)*y(381)...
    +p.phoCpdR*(y(381-100)+y(381-90))*y(381-20)-p.dephoCpdR*(y(381-80)+y(381-70))*y(381)-p.ukdephoCpdR*y(381)...%
    +p.D_CpdR*(y(382)-y(381))/(y(396)^2);



%Bin 382-389
for i=382:389
    if i==385
    dydt(i)=-(p.mu+p.deg_cpdR)*y(i)...
    +p.phoCpdR*(y(i-100)+y(i-90))*y(i-20)-p.dephoCpdR*(y(i-80)+y(i-70))*y(i)-p.ukdephoCpdR*y(i)...%
    +p.D_CpdR*(-y(i)+y(i-1))/(y(396)^2);
    elseif i==386
        dydt(i)=-(p.mu+p.deg_cpdR)*y(i)...
    +p.phoCpdR*(y(i-100)+y(i-90))*y(i-20)-p.dephoCpdR*(y(i-80)+y(i-70))*y(i)-p.ukdephoCpdR*y(i)...%
    +p.D_CpdR*(y(i+1)-y(i))/(y(396)^2);
    else
        dydt(i)=-(p.mu+p.deg_cpdR)*y(i)...
    +p.phoCpdR*(y(i-100)+y(i-90))*y(i-20)-p.dephoCpdR*(y(i-80)+y(i-70))*y(i)-p.ukdephoCpdR*y(i)...%
    +p.D_CpdR*(y(i+1)-2*y(i)+y(i-1))/(y(396)^2);
    end
end

%Bin 390
dydt(390)=-(p.mu+p.deg_cpdR)*y(390)...
    +p.phoCpdR*(y(390-100)+y(390-90))*y(390-20)-p.dephoCpdR*(y(390-80)+y(390-70))*y(390)-p.ukdephoCpdR*y(390)...%
    +p.D_CpdR*(y(389)-y(390))/(y(396)^2);

%% S
for i=391:395
% dydt(i)=0;
dydt(i)=0;
end%S - PodJ, CtrA, PleC, PerP, DivL
%% cell growth equation
dydt(396)=p.mu*y(396);
end
end
