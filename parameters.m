function p=parameters(cond,ver,mutant)


global p;
% range_value is synthesis rate of PopZ monomer 
%cond is for growing or ungrowing cell
% Growth conditions
if cond == 0              % cell size is fixed
    p.mu = 0;             % growth rate constant    
elseif cond == 1          % cell is growing 
    p.mu = 0.0053;%0.00526;%0.0055;        % growth rate constant = 0.0055 (units => 1/min)    
end

if ver==1
       H1=0.8; H2=0.8;
%     H1=0.9; H2=0.9;
                     p.D_podJm =100;%for PopZ only & 100 in manuscript;  
p.D_podJL =0.0005;%40;%0.0005;  
 
p.syn_podJ = 0.01;%0.018;%0.2*theta(1);   % monomer synthesis 0.2
% p.syn_podJ2 = 0.012+0.0025;
% p.Ji_PodJCtrA = 0.9;
p.n_PodJCtrA = 4;%2;
p.deg_podJm =0.05;% 0.05;    0.01?- 0.046 paper116
p.deg_podJp =p.deg_podJm;% 0.05;    
p.deg_podJ1 =0.007;%0.007;%0.007;%paper116
 
p.dnv_podJ =2;%0.1*theta(3);% 15              % denovo polymerization 12 in dissertation for table C.2; 1 for table C.1
p.aut1_podJ =75;%+20;%10;               % autocatalytic polymerization (pole)
p.aut1_podJ1=H1*p.aut1_podJ;%central compartment autcatalytic polymerization
p.depol_podJ =0.2;%0.5;             % deplymerization 
 
p.podj=2;
p.deg_s=0.05;%0.05  %0.01 paper116ref
   p.alpha_PodJSpmX =30;%2; %SpmX on PodJ dnv
   
   
%% PopZ
p.syn_popz=0.24;%*theta(8);%1.5;% popZ monomer synthesis
    p.deg_popzm=0.05;%theta(9);%0.05 % popZ monomer degradation
    p.deg_popzp=p.deg_popzm;%theta(9);%0.05;
    p.dnv_popz=2.5;%theta(10);%60;% denovo polymerization
    p.aut_popz=15;%3*theta(11);%3;% autocatalytic polymerization
     p.aut_popz1=H2*p.aut_popz;
       p.alpha_PopZPodJ=100;%1;%PodJ on PopZ
    p.depol_popz=0.15;%0.1; % deplymerization 
    % diffussion rates
    p.D_popzm =835;%750; % PopZ monomer dffusion  
p.D_popzp =0.0005;%40;%0.0005;  % PopZ polymer diffusion

   %% SpmX
%  p.syn_spmx =0.006;%0.04;
   p.deg_spmx = 0.01;%0.1;
     p.dnv_spmx =0.001;%0.1;%
   p.depol_spmx = 0.4*p.dnv_spmx;

   p.aut_spmx=1/10;%50;
%    p.Ja_SpmXCtrA=0.2;
   p.alpha_SpmXPopZ=50;
p.D_spmx =200;%um^2/min
p.D_spmxp=0.0005;

%% CtrA
% p.syn_ctrA1=0.0083*10;%0.0083;%Shenghua2008-exp; 0.026;%Murray-fitted
p.Ji_CtrACtrA=1;
% p.syn_ctrA2=0.073*5;%0.073;%Shenghua2008-exp
% p.Ja_CtrACtrA=8;
p.deg_ctrA1=0.0038;%0.0038;%exp-calculated
% p.deg_ctrA2=0.05;%0.1 exp-estimated
% p.Jd_CpdR=0.3;
p.dephoCtrA=0.14e3;%Murray-estimated~0.14
p.phoCtrA=0.077e3;%7/93*p.dephoCtrA/0.5%~0.021;%Murray-estimated
 p.ukdephoCtrA=0;
p.D_CtrA=427;
p.D_CtrAP=p.D_CtrA;
%% PleC
% p.syn_pleC=0.01;%0.053%Bronson
p.deg_pleC=0.02;%0.028;%Bronson
% p.fb_PleC=100;
p.bf_PleC=0.1;
p.D_PleC=71;
% p.ph2kin_PleC1=4;%2;
% p.ph2kin_PleC2=0.4;
% p.b_PleCDivKP=1;
p.ub_PleCDivKP=2.4/2;
p.D_PleCDivK=63;

%% DivJ
p.syn_divJ=0.008;%0.005;%Bronson
p.deg_divJ=0.035;%0.035;%Bronson
p.fb_DivJ=20;
p.bf_DivJ=0.5;
p.D_DivJ=108;
%%%%%%%%
p.b_DivJDivK=1;
% p.ub_DivJDivK=8/3;
p.b_DivJDivKP=1;
% p.ub_DivJDivKP=16/3;
p.D_DivJDivKP=84.5;


%% DivK
p.syn_divK1=0.001;%0.0004;%Bronson
% p.syn_divK2=0.05;%0.125;%Bronson
% p.Ja_DivKCtrA=2;%1.7/2;%Bronson
p.deg_divK=0.014;%Bronson; 0.002;%Shenghua2008-exp
% p.depho_DivK=4;%PleCph
% p.pho_DivKDivJf=0.1;%0.4;
% p.pho_DivKDivJb=8;
% p.pho_DivKPleC=0.2;
p.D_DivK=1319;
p.D_DivKP=p.D_DivK;%p.D_DivK;

%% DivL
% p.syn_divL =0.05;
% p.deg_divL =0.3;
% p.b_DivLDivKP = 5;
p.ub_DivLDivKP =1.1;% 11;
p.fb_divL=1;
p.bf_divL=1;
% p.alpha_DivLPopZ=0.001;
% p.alpha_DivLPodJ=2+2;
p.D_DivL=76.6;
p.D_DivLDivK=66.4;

%% PerP
% p.syn_perP=0.5;%0.0428;%Bronson
% p.Ja_PerPCtrA=0.2;%2.8/2;%Bronon
p.deg_perP=0.04;%0.04;%Bronson
p.D_PerP=853;

%% CckA
% p.syn_cckA=0.1;
p.deg_cckA=0.2;
p.fb_cckA=1;
% p.alpha_CckAPopZ=0.001;
% p.alpha_CckADivL=100;
p.bf_cckA=0.5;
p.kp_cckA1=0;%0.1;
% p.pk_cckA1=0;%0.1;
% p.kp_cckA2=1*200;%
% p.pk_cckA2=1*20;
p.D_CckA=87.3;
p.b_CtrAPCckAph=1*3;
% p.ub_CtrAPCckAph=0.4;
p.b_CtrACckAkin=3;
% p.ub_CtrACckAkin=0.4;
p.D_CckACtrA=65;

%% CpdR
% p.syn_cpdR=1;
% p.deg_cpdR=0.2;
% p.JaCpdRCtrA=0.5;
% p.fb_cpdR=0.1;
p.bf_cpdR=7+5;
% p.phoCpdR=7;
p.dephoCpdR=1*4;
p.D_CpdR=1.6386e+03;
p.Jd_CpdR1=0;
p.ukdephoCpdR=0.5;  



%%%%%%%%%%%%%%%%%%%%%%%%%
load('MultiGA_Output.mat','val');
para=val(1,:);

% p.syn_podJ = 0.009;
para(19)=145.00;
p.syn_podJ2=para(1);   p.Ji_PodJCtrA=para(2);  p.syn_spmx=para(3);  p.Ja_SpmXCtrA=para(4); p.syn_ctrA1=para(5); 
p.syn_ctrA2=para(6);  p.Ja_CtrACtrA=para(7);  p.deg_ctrA2=para(8); p.Jd_CpdR=para(9); p.syn_pleC=para(10); 
p.fb_PleC=para(11); p.b_PleCDivKP=para(12); p.ub_DivJDivKP=para(13); p.ub_DivJDivK=para(14); p.syn_divK2=para(15); 
p.Ja_DivKCtrA=para(16); p.pho_DivKPleC=para(17);  p.pho_DivKDivJf=para(18); p.pho_DivKDivJb=para(19);  p.depho_DivK=para(20);
p.syn_perP=para(21);  p.Ja_PerPCtrA=para(22); p.syn_divL=para(23);  p.deg_divL=para(24);  p.alpha_DivLPopZ=para(25); 
p.alpha_DivLPodJ=para(26);  p.b_DivLDivKP=para(27); p.syn_cckA=para(28); p.alpha_CckAPopZ=para(29); p.alpha_CckADivL=para(30);
p.pk_cckA1=para(31);  p.kp_cckA2=para(32); p.pk_cckA2=para(33); p.ub_CtrACckAkin=para(34); p.ub_CtrAPCckAph=para(35);
p.syn_cpdR=para(36);  p.deg_cpdR=para(37);  p.JaCpdRCtrA=para(38);  p.fb_cpdR=para(39);  p.phoCpdR =para(40);

p.ph2kin_PleC1=para(42);
p.ph2kin_PleC2=0.1*para(42);

p.TS = para(41);



end
if strcmp(mutant,'deltaPodJ')
    p.syn_podJ = 0;
p.syn_podJ2 = 0;
elseif strcmp(mutant,'PodJ+')
    p.syn_podJ = p.syn_podJ*1.2;
p.syn_podJ2 = p.syn_podJ2*1.2;
elseif strcmp(mutant,'PodJ++')
     p.syn_podJ = p.syn_podJ*12;%3.5;
p.syn_podJ2 = p.syn_podJ2*12;%3.5;
elseif strcmp(mutant,'deltaSpmX')
    p.syn_spmx=0;
elseif strcmp(mutant,'deltaPopZ')
    p.syn_popz=0;
    elseif strcmp(mutant,'PopZ+')
     p.syn_popz =  p.syn_popz*12;
elseif strcmp(mutant,'SpmX+')
    p.syn_spmx=50*p.syn_spmx;
elseif strcmp(mutant,'deltaMmpA')
    p.deg_s=0*p.deg_s;
    elseif strcmp(mutant,'deltaPerP')
    p.syn_perP=0*p.syn_perP;
elseif strcmp(mutant,'deltaDivJ')
   p.syn_divJ=0; %p.alpha_DivLPopZ=0; % p.pho_DivKPleC=0;  %p.b_DivJDivKP=0; p.alpha_CckAPopZ=0;
    elseif strcmp(mutant,'deltaPleC')
    p.syn_pleC=0;
elseif strcmp(mutant,'deltaDivJ&deltaPleC')
    p.syn_divJ=0; p.syn_pleC=0;
        elseif strcmp(mutant,'PleC-H610A')
   p.depho_DivK=0; p.pho_DivKPleC = 0;
elseif strcmp(mutant,'DivJ-H338A')
    p.pho_DivKDivJb=0; p.pho_DivKDivJf=0;
    elseif strcmp(mutant,'DivK-D53A')
    p.pho_DivKDivJb=0; p.pho_DivKDivJf=0;  p.pho_DivKPleC=0; 
 elseif strcmp(mutant,'DivK-D90G')||strcmp(mutant,'DivL_A601')
    p.b_DivLDivKP=0;%.2*p.b_DivLDivKP;
     elseif strcmp(mutant,'PleC-F778L')
    p.pho_DivKPleC=0.01*p.pho_DivKPleC; %0; 
     elseif strcmp(mutant,'DivK-D90G&deltaPleC')
    p.b_DivLDivKP=0;%.3*p.b_DivLDivKP;%p.depho_DivK=0;
    p.syn_pleC=0;
    elseif strcmp(mutant,'p5:overDivLDivKbinding')||strcmp(mutant,'DivL_Y550F')
    p.b_DivLDivKP=1.5*p.b_DivLDivKP;  %p.pk_cckA2
%         p.b_DivLDivKP=1.2*p.b_DivLDivKP;
elseif strcmp(mutant,'p1:deletingPleCphosphatase')
     p.depho_DivK=0.1*p.depho_DivK;
elseif strcmp(mutant,'p2:deletingPleCDivKPbinding')
    p.b_PleCDivKP=0;
elseif strcmp(mutant,'deltaDivL')
    p.syn_divL=0;
% elseif strcmp(mutant,'p3:deletingDivLDivK')
elseif strcmp(mutant,'p3:deletingPodJPleCbinding')
    p.fb_PleC=0;
elseif strcmp(mutant,'p4:deletingPodJDivLbinding')
    p.alpha_DivLPodJ=0;
    elseif strcmp(mutant,'p6:decreaseDivLDivKbinding')
    p.b_DivLDivKP=0.5*p.b_DivLDivKP; %p.kp_cckA2=3*p.kp_cckA2;% p.kp_cckA2=1;%
%     p.b_DivLDivKP=1.5*p.b_DivLDivKP; p.kp_cckA2=0.2*p.kp_cckA2;%p.kp_cckA2=1;%

    elseif strcmp(mutant,'p7:decreasePleCDivKbinding')
    p.b_PleCDivKP=0.1*p.b_DivLDivKP; %p.depho_DivK=5*p.depho_DivK;%1&2
%     p.b_PleCDivKP=1.5*p.b_DivLDivKP; p.depho_DivK=0.1*p.depho_DivK;%1&2
     elseif strcmp(mutant,'p8:DivJDivKbinding')
   p.b_DivJDivK=0.1*p.b_DivJDivK; p.b_DivJDivKP=0.1*p.b_DivJDivKP;
   p.pho_DivKDivJf=100*p.pho_DivKDivJf;  p.pho_DivKDivJb=100*p.pho_DivKDivJb;%1&2
%     p.b_PleCDivKP=1.5*p.b_DivLDivKP; %p.depho_DivK=0.1*p.depho_DivK;%1&2
elseif strcmp(mutant,'p9:p3andp4')
    p.fb_PleC=0;
    p.alpha_DivLPodJ=0;
end





