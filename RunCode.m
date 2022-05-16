
clear all
close all


%% version of params
ver=9;

%% cell type: 'SW' 
celltype='SW';
%% mutant type
mutant='WT';%
 %mutantlist1: 'deltaPodJ';PodJ+;deltaSpmX; SpmX+; deltaPopZ; deltaMmpA; deltaPerP
 %mutantlist2: 'deltaDivJ'  deltaPleC deltaDivJ&deltaPleC PleC-H610A
 %DivJ-H338A DivK-D53A DivK-D90G  PleC-F778L DivL_A601 DivL_Y550F
 %'DivK-D90G&deltaPleC'
 
%   mutant='PleC-F778L';%
 %%%prediction%%%
 %'p1:deletingPleCphosphatase'   'p2:deletingPleCDivKPbinding'  'p5:overDivLDivKbinding'
%  p6:decreaseDivLDivKbinding  p4:deletingPodJDivLbinding  p3:deletingPodJPleCbinding
%p7:decreasePleCDivKbinding  p8:DivJDivKbinding
%  mutant='p4:deletingPodJDivLbinding';


%% initial values
%initial input of initial values:
% y0=zeros(396,1);%SW IC - first cell cycle
 
% y0(11:20)=10e-6;%PodJp
% y0(21:28)=0.001; y0(29:30)=0.1;%PodJS
% y0(61:68)=0.5; y0(69:70)=2;%25;%PopZp
% y0(71:80)=0.2; y0(81:90)=0.5;%CtrA and CtrAP
% y0(91:100)=0.05; %PleCf
% y0(44)=0.1;%PleCb
% % y0(77:80)=0.2;%DivK
% % y0(76)=0.05;%DivJb:DivKP
% y0(97:100)=0.01;%DivLf
% y0(113:128)=0.001;%All single CckA
% y0(145:148)=0.1;%CpdRf
% y0(153:156)=0.1;%CpdRP
% y0(396)=0.02*20; y0(163)=0.02*30;

load('y0_10com_4.mat')
 yori=y0;
   


%%%%%%%%%%%%%%%%%%%%%%%%%%%timer

G=1; %run simulations
if G==1
 CycleNum=3; %number of cell cycles
 for i=1:CycleNum
 TITLE= [num2str(i) 'cellcyle'];

[Y, time, y0_,TE,IE,Y1,time1,Y2,time2,DNArep]=main1(y0,celltype,ver,mutant);%simulation
resultgraph10com1(Y,time,celltype,mutant,TITLE,1) % same WT figures but in different format
result_PlotMutant3(Y,time,celltype,mutant,TITLE,1) %plot mutant figures

if DNArep
y0=IniValue(Y,celltype);%calculate the initial values for next cycle
end
 end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mutant sim
M=0;  %calculate DivK~P in different mutants
if M==1

mutantlist = {'WT','deltaSpmX','deltaDivJ', 'deltaPleC', 'deltaDivJ&deltaPleC', 'deltaPodJ'};


DivKPlist =[];
DivKlist = [];
KP2U = [];KP2T = [];
CtrAPlist = [];
CtrAlist = [];


for i=1:length(mutantlist)

mutant=mutantlist(i);

CycleNum=2;  %output the DivK~P value in the second cycle
 for j=1:CycleNum
[Y, time, y0_,TE,IE,Y1,time1,Y2,time2,DNArep]=main1(y0,celltype,ver,mutant);%simulation
if DNArep
y0=IniValue(Y,celltype);%calculate the initial values for next cycle
end
 end
y0=yori;
[DivKP, DivK,CtrAP, CtrA] = SimVSExp (Y,time);
DivKPlist = [DivKPlist; DivKP];
DivKlist = [DivKlist; DivK];
KP2U = [KP2U; DivKP/DivK];
KP2T = [KP2T; DivKP/(DivK+DivKP)];
CtrAPlist = [CtrAPlist; CtrAP];
CtrAlist = [CtrAlist; CtrA];
end
DivKPlist = DivKPlist./DivKPlist(1);
DivKlist = DivKlist./DivKlist(1);
KP2U = KP2U./KP2U(1);
KP2T = KP2T./KP2T(1);
CtrAPlist = CtrAPlist./CtrAPlist(1);
CtrAlist = CtrAlist./CtrAlist(1);
T1=table(mutantlist.',DivKPlist,DivKlist,KP2T,CtrAPlist, CtrAlist)


end