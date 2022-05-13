 function resultgraph10com1(yout,tout,celltype,mutant,TITLE,concentration)
% close all
%spatial plots
%%simulation
%% Scaffolding spatial and temporal
PodJLp=yout(11:20,:);
PodJLp = flipud(PodJLp);
PodJm=yout(1:10,:);
PodJm = flipud(PodJm);
PodJL=PodJLp+PodJm;
PodJS=yout(21:30,:);
PodJS = flipud(PodJS);
PopZp = yout(61:70,:);
PopZp = flipud(PopZp);
PopZ = yout(51:60,:)+yout(61:70,:);
PopZ = flipud(PopZ);
SpmX=yout(31:40,:)+yout(41:50,:);
SpmX = flipud(SpmX);



DivJbT=yout(161:170,:)+yout(201:210,:)+yout(181:190,:);
DivJT=yout(191:200,:)+yout(151:160,:)+yout(171:180,:)+DivJbT;
DivJT = flipud(DivJT);
DivJbT = flipud(DivJbT);

PleCf=yout(91:100,:);
PleCb=yout(101:110,:);
PleCfDivKP=yout(111:120,:);
PleCbDivKP=yout(121:130,:);
PleCfkin=yout(131:140,:);
PleCbkin=yout(141:150,:);

PleCphT=PleCf+PleCb+PleCfDivKP+PleCbDivKP;%%%%%%%%%%%%%%%%???
PleCkinT=PleCfkin+PleCbkin;
PleCphT=flipud(PleCphT);
PleCkinT=flipud(PleCkinT);
PleC=PleCkinT+PleCphT;

DivK=yout(221:230,:)+yout(171:180,:)+yout(181:190,:);
DivK = flipud(DivK);
DivKP=yout(231:240,:);
DivLDivKPT=yout(261:270,:)+yout(271:280,:);
DivKPT=DivKP+PleCfDivKP+PleCbDivKP+yout(191:200,:)+yout(201:210,:)+DivLDivKPT;
DivKP = flipud(DivKP);
DivKPT = flipud(DivKPT);
DivKtot = DivK+DivKPT;

DivJf = yout(151:160,:) + yout(171:180,:) + yout(191:200,:);
DivJb = yout(161:170,:) + yout(181:190,:) + yout(201:210,:);
DivJf = flipud(DivJf);
DivJb = flipud(DivJb);
DivJtot = DivJf + DivJb;

CtrACckAkin=yout(341:350,:)+yout(351:360,:);
CtrAPCckAph=yout(321:330,:)+yout(331:340,:);
CtrAP=yout(81:90,:)+CtrAPCckAph;
CtrAP = flipud(CtrAP);
CtrA=yout(71:80,:)+CtrACckAkin;
CtrA = flipud(CtrA);
CtrAt = CtrA+CtrAP;
PerP=yout(211:220,:);
PerP=flipud(PerP);

DivLT=yout(241:250,:)+yout(251:260,:);
DivLT=flipud(DivLT);
DivLDivKPT=flipud(DivLDivKPT);
DivLtot = DivLT+DivLDivKPT;
CckAkinT=yout(281:290,:)+yout(291:300,:)+CtrACckAkin;
CckAkinT=flipud(CckAkinT);
CckAphT=yout(301:310,:)+yout(311:320,:)+CtrAPCckAph;
CckAphT=flipud(CckAphT);
CckAT=CckAkinT+CckAphT;

CpdRu=yout(361:370,:)+yout(371:380,:);
CpdRP=yout(381:390,:);
CpdRu=flipud(CpdRu); CpdRP=flipud(CpdRP);

%%Mgrid
Y=yout';
L=length(Y(:,1));
M(:,6)=zeros(1,L);
M(:,5)=-Y(:,396);
M(:,4)=-2*Y(:,396);
M(:,3)=-3*Y(:,396);
M(:,2)=-4*Y(:,396);
M(:,1)=-5*Y(:,396);
M(:,7:11)=-fliplr(M(:,1:5));
M=M';
a=zeros(1,L);
time=tout;

%  if strcmp(celltype,'ST')
% time=time+34;
% end

%% plots
% ha = tight_subplot(2,4,[.06 .03],[.1 .04],[.02 .02]);
% model='flat'; %''
figure()%%PodJ, SpmX, PopZ
% set(gcf,'position',[1000 100 600 600])
set(gcf,'position',[100 100 1800 300])%left, lower, right, upper

set(gcf,'Name',TITLE);
% ha = tight_subplot(5,2,[.1 .03],[.1 .04],[.1 .04]);
ha = tight_subplot(1,4,[.1 .03],[.1 .1],[.05 .04]);

axes(ha(1)); 
PodJL(11,:)=a;
pcolor(time, M, PodJL)
shading flat
% shading interp
colorbar
caxis([0 2]);
% xlabel('time (min)')
title('(a) PodJL')

axes(ha(2));
PodJS(11,:)=a;
pcolor(time, M, PodJS)
shading flat
colorbar
caxis([0 2]);
% xlabel('time (min)')
title('(b) PodJS')


axes(ha(3)); 
PopZ(11,:)=a;
pcolor(time, M, PopZ)
shading flat
colorbar
% caxis([0 4]);
% xlabel('time (min)')
title('(c) PopZ')

axes(ha(4)); 
SpmX(11,:)=a;
pcolor(time, M, SpmX)
shading flat
colorbar
% caxis([0 0.25]);
% xlabel('time (min)')
title('(d) SpmX')

set(findall(gcf,'-property','FontSize'),'FontSize',14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()%%PleC DivJ DivK DivL
% set(gcf,'position',[1000 100 600 600])
set(gcf,'position',[100 100 2500 800])%left, lower, right, upper
set(gcf,'Name',TITLE);
% ha = tight_subplot(5,2,[.1 .03],[.1 .04],[.1 .04]);
ha = tight_subplot(2,4,[.1 .03],[.1 .1],[.05 .04]);

axes(ha(1)); 
PleC(11,:)=a;
pcolor(time, M, PleC)
shading flat
colorbar
% caxis([0 0.2]);
title('(a) PleC')

axes(ha(2));
DivJtot(11,:)=a;
pcolor(time, M, DivJtot)
shading flat
colorbar
% caxis([0.1 0.4]);
% xlabel('time (min)')
title('(b) DivJ')


axes(ha(3)); 
% subplot(1,4,3)
DivKPT(11,:)=a;
pcolor(time, M, DivKPT)
shading flat
colorbar
% caxis([0.1 0.35]);
% xlabel('time (min)')
title('(c) DivKPT')



axes(ha(4)); 
DivLT(11,:)=a;
pcolor(time, M, DivLT)
shading flat
colorbar
% caxis([0.002 0.017]);
% xlabel('time (min)')
title('(d) DivLT')

axes(ha(5));
DivLDivKPT(11,:)=a;
pcolor(time, M, DivLDivKPT)
shading flat
colorbar
% caxis([0.002 0.017]);
% xlabel('time (min)')
title('(e) DivLDivKPT')

axes(ha(6));
CckAT(11,:)=a;
pcolor(time, M, CckAT)
shading flat
colorbar
% caxis([0.05 0.18]);
% xlabel('time (min)')
title('(f) CckA')


axes(ha(7)); 
CpdRu(11,:)=a;
pcolor(time, M, CpdRu)
shading flat
colorbar
% caxis([0.6 1]);
% xlabel('time (min)')
title('(g) CpdRu')

axes(ha(8));
CtrAP(11,:)=a;
pcolor(time, M, CtrAP)
shading flat
colorbar
% caxis([0.05 0.55]);
% xlabel('time (min)')
title('(h) CtrAP')
set(findall(gcf,'-property','FontSize'),'FontSize',14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MUTANT

%%%%%%%%%%%%%%%%%%%
tspan=[tout(1) tout(end)];
%% temperal DATA
SpmX = [0 2340.397; 20 11034.71; 40 13828.5; 60 14676.86; 80 13635.98; 100 9316.033; 120 8849.205; 140 12085.38; 160 13972.28];
SpmX(:,2) = SpmX(:,2)/max(SpmX(:,2));
dpSpmX=SpmX;%paper152
SpmX2 = [0 2857.246; 20 10014.882; 40 13874.083;60 13516.276;80 13973.64; 100 8493.619; 120 8413.962; 140 11470.74;160 11327.459];
SpmX2(:,2) = SpmX2(:,2)/max(SpmX2(:,2));
dpSpmX2=SpmX2;%paper118
SpmX3 =[0 1506.527;20 4945.134;40 8203.305;60 9740.79;80 9703.205;100 10964.376;120 10477.861;140 9446.033;160 4349.246];%paper138
SpmX3(:,2) = SpmX3(:,2)/max(SpmX3(:,2));
dpSpmX3=SpmX3;
PodJL=[0 0; 20 0; 40 1726.012; 60 12068.518; 80 23215.326; 100 22996.104; 120  9466.861; 140  3336.154; 160 5506.548];
PodJL(:,2) = PodJL(:,2)/max(PodJL(:,2));
dpPodJL=PodJL;%paper116
PodJL2=[0 790.719;20 1297.841;40 4612.861;60 9123.397;80 11505.154;100 13792.497;120 15360.276;140 8009.589;160 5775.439];
PodJL2(:,2) = PodJL2(:,2)/max(PodJL2(:,2));
dpPodJL2=PodJL2;%138
PodJS=[0  16869.841; 20  16960.447; 40  4387.841; 60  434.364;  80  0; 100  2627.527; 120  8780.79; 140  8984.033; 160  6333.598];
PodJS(:,2) = PodJS(:,2)/max(PodJS(:,2));
dpPodJS=PodJS;%paper116
PodJS2=[0 7656.953;20 11798.569;40 9196.033;60 7763.104;80 7449.276;100 8469.74; 120 11575.711; 140 11777.761; 160 10382.853];
PodJS2(:,2) = PodJS2(:,2)/max(PodJS2(:,2));
dpPodJS2=PodJS2;
dpDivKTOT =[0 0.67; 20 0.77; 40 0.88; 60 0.87; 80 0.72; 100 0.71; 120 0.89; 140 1];
DivKTOT2=[0 4395.619;20 5939.134;40 6636.012;60 6630.841;80 6025.305;100 5316.184;120 6691.891;140 7454.719;160 5877.569];
DivKTOT2(:,2) = DivKTOT2(:,2)/max(DivKTOT2(:,2));
dpDivKTOT2=DivKTOT2;
DivKTOT3 = [0 17730.68; 24 14609.71; 48 12892.2; 72 15597.2; 96 14817.9; 120 14020.73];%PAPER135
DivKTOT3(:,2) = DivKTOT3(:,2)/max(DivKTOT3(:,2));
dpDivKTOT3=DivKTOT3;
DivKP = [0 19388.08; 30 23675.43; 60 22634; 90 25767.41; 120 26369.41];%PAPER135
DivKP(:,2) = DivKP(:,2)/max(DivKP(:,2));
dpDivKP=DivKP;

dpDivJ = [0 0.364; 20 0.574; 40 0.81; 60 0.997; 80 0.924; 100 0.9; 120 1; 140 0.94]; % Wheeler et al. 1999
dpDivJ2 = [0 0.26; 20 0.56; 40 0.86; 60 0.81; 80 0.97; 100 0.93; 120 0.86; 140 1]; %Sanselicio 2015
DivJ3 = [0 2638.933; 20 5578.497; 40 7750.447;60 10436.912;80 11503.497; 100 9782.79; 120 13502.083; 140 12420.326;160 13399.489];
DivJ3(:,2) = DivJ3(:,2)/max(DivJ3(:,2));%paper118
dpDivJ3=DivJ3;
dpPleCtot = [0 .85/0.85; 20 0.5/0.85; 40 .28/0.85; 60 .38/0.85; 80 .45/0.85; 100 .65/0.85; 120 0.75/0.85; 140 0.8/0.85]; %Viollier 2002 - ommitted last datapoint bc it was at time point 160 and is unlikely behavior
PleCtot2 =[0 7535.004;20 4253.184; 40 3452.841; 60 5912.426; 80 9557.225; 100 11047.447; 120 11043.326; 140 11704.104; 160 6926.953];
PleCtot2(:,2) = PleCtot2(:,2)/max(PleCtot2(:,2));
dpPleCtot2=PleCtot2;
dpCtrAP = [10*1.07 0.4; 80*1.07 0.5; 97*1.07 .8; 130*1.07 1]; %Jacobs 2003 (relative to max = 1) removed second point (t=40, P=0.15) as CtrAT suggests it is incorrect
Mcgrath = [0	12915.912; 20	9480.648; 40	2687.648; 60	515.698; 80 5769.719; 100	9782.426; 120	13602.205;140	1398.439];%Bronson extract
% Mcgrath = [10	0.8; 27	0.6; 45	0.2; 62	0.04; 81	0.4; 98	0.7; 115 0.9; 133	1; 150	0.8];%paper4 I extract
Mcgrath(:,2)=Mcgrath(:,2)/max(Mcgrath(:,2));
dpCtrAT = Mcgrath;%SLOW
dpCtrAT2 = [0 0.8; 20 0; 40 0; 60 0; 80 0.5; 100 1; 120 .85; 140 0.7]; %Collier 2006 - DnaA couples DNA... (Normalized based on max expression)
%QUICK
% CtrAT3=[0 16004.861;20 7879.962;40 1532.690;60 939.447;80 5083.154;100 11343.933; 120 12177.326; 140 13739.397; 160 13562.326];%paper118
CtrAT3 = [0 16183.619; 20 8066.79; 40 1543.861; 60	865.669; 80		5211.74; 100 11540.64; 120	12304.376; 140	14119.569; 160	23971.874];

CtrAT3(:,2) = CtrAT3(:,2)/max(CtrAT3(:,2));
dpCtrAT3=CtrAT3;

CckAkin1 = [0 0.59; 10 0.34; 20 0.29;  40 0.2; 50 0.25; 60 0.38; 80 0.64; 90 0.83; 100 1; 110 1; 120 1; 130 1];%159,fig3A
CckAkin1(:,2) = CckAkin1(:,2)/max(CckAkin1(:,2));
dpCckAkin1=CckAkin1;

CckAkin2 = [10 1145.741; 30 2034.406; 50 519.527; 70 1077.083; 90 14737.3; 110 19101.25; 130 21106.25];%159,fig3B
CckAkin2(:,2) = CckAkin2(:,2)/max(CckAkin2(:,2));
dpCckAkin2=CckAkin2;

dpCpdRtot = [1*0.93 0.97; 25*0.93 1; 50*0.93 0.66; 75*0.93 0.497; 100*0.93 0.502; 125*0.93 0.64; 150*0.93 0.96]; %paper9, Iniesta et al. 2006 Figure 5A by Bronson
dpCpdR = [1*0.93 0.42; 25*0.93 0.938; 50*0.93 0.636; 75*0.93 0.48; 100*0.93 0.07; 125*0.93 0.275; 150*0.93 0]; %Iniesta et al. 2006 Figure 5A estimated by Bronson
dpCpdRtot2 = [0 0.6; 20*0.93 0.41; 40*0.93 0.22; 60*0.93 0.19; 80*0.93 0.19; 100*0.93 0.32; 120*0.93 0.67; 140*0.93 1; 160*0.93 0.83];%paper cpdr1 figS1

%%%%%%%%%%%%%%mutant temporal data
PodJLmutant=[0 0; 20 0; 40 1.81; 60 33.3; 80 99.01; 100 37.09; 120  3.46; 140  13.9; 160 64.48];
PodJLmutant(:,2) = PodJLmutant(:,2)/max(PodJLmutant(:,2));
dpPodJLmutant=PodJLmutant;%paper116 deltaMmpA
PodJSmutant=[0 84.87; 20 95.31; 40 82.08; 60 48.44; 80 62.17; 100 94.98; 120  96.88; 140  99.43; 160 98.02];
PodJSmutant(:,2) = PodJSmutant(:,2)/max(PodJSmutant(:,2));
dpPodJSmutant=PodJSmutant;%paper116 deltaMmpA









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% temporal plots
%%
% if strcmp(mutant,'WT')
if strcmp(mutant,mutant)
LegendSize=8; LabelSize=14;

%  if concentration==1
%     factor1=0.2; factor2=0.3;
% else
%     factor1=yout(162,:); factor2=yout(163,:);
% end


  figure()%DivJ, PleC, DivK, CtrA
set(gcf,'position',[100 400 1400 450])%left, lower, right, upper

ha = tight_subplot(2,3,[.15 .05],[.1 .1],[.05 .05]);
% ha = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
%
%   in:  Nh      number of axes in hight (vertical direction)
%        Nw      number of axes in width (horizontaldirection)
%        gap     gaps between the axes in normalized units (0...1)
%                   or [gap_h gap_w] for different gaps in height and width 
%        marg_h  margins in height in normalized units (0...1)
%                   or [lower upper] for different lower and upper margins 
%        marg_w  margins in width in normalized units (0...1)
%                   or [left right] for different left and right margins 
%
%  out:  ha     array of handles of the axes objects
%                   starting from upper left corner, going row-wise as in
%                   going row-wise as in


% axes(ha(1));
subplot(2,3,1)
PodJL=yout(1:10,:)+yout(11:20,:);   PodJS=yout(21:30,:);
% PodJL_T=sum(PodJL);  PodJS_T=sum(PodJS);
PodJL_T=sum(PodJL);  PodJS_T=sum(PodJS); 
PodJ_T=PodJL_T+PodJS_T;
Fraction_L=PodJL_T./PodJ_T;
Fraction_S=PodJS_T./PodJ_T;
%%data
Time=[0 20 40 60 80 100 120 140 160];
%%Calculate fraction of relative levels by getdata
L=[100 91 19 2 1 12 46 32 53];
S=[0 0 7 50 100 100 41 14 22];
Sum=L+S;

SpmX=sum(yout(31:40,:)+yout(41:50,:));

p5 = line(tout,PodJL_T , 'Color', 'k' , 'LineWidth', 2, 'Linestyle', '-');
p5 = line(tout,PodJS_T , 'Color', 'b' , 'LineWidth', 2, 'Linestyle', '-');

hold on
[score,scalary] = findSquares(tout,PodJL_T, dpPodJL(1:end-1,:));
p6 = scatter(dpPodJL(:,1),dpPodJL(:,2).*scalary, 'k+', 'LineWidth', 1);
hold on
[score,scalary] = findSquares(tout,PodJL_T, dpPodJL2(1:end-1,:));
p6 = scatter(dpPodJL2(:,1),dpPodJL2(:,2).*scalary, 'k^', 'LineWidth', 1);
hold on
[score,scalary] = findSquares(tout,PodJS_T, dpPodJS(1:end-1,:));
p7 = scatter(dpPodJS(:,1),dpPodJS(:,2).*scalary, 'b+', 'LineWidth', 1);
hold on
[score,scalary] = findSquares(tout,PodJS_T, dpPodJS2(1:end-1,:));
p7 = scatter(dpPodJS2(:,1),dpPodJS2(:,2).*scalary, 'b^', 'LineWidth', 1);

hold off
axis([tspan 0 max(PodJL_T+PodJS_T)*1.2]);
% xlabel('time (min)','fontsize',LabelSize)
ylabel('scaled concentration')
xlabel('(a)','fontsize',LabelSize)
% h = legend('SpmX_T','Kaczmarczyk\newlineet al. ''20','Location', 'North');
h = legend('PodJL','PodJS','Chen 05','Guo 04','Chen 05','Guo 04','Location', 'North','fontsize',LegendSize,'NumColumns',3);
  set(h, 'Interpreter', 'tex', 'Box', 'on', 'Orientation', 'vertical');
  
  
% axes(ha(2));
subplot(2,3,2)
  p5 = line(tout,SpmX, 'Color', 'k' , 'LineWidth', 2, 'Linestyle', '--');
hold on
% [score,scalary] = findSquares(tout,sum(SpmX), dpSpmX2(1:end-1,:));
[score,scalary] = findSquares(tout,SpmX, dpSpmX2(1:end-1,:));
p7 = scatter(dpSpmX2(:,1),dpSpmX2(:,2).*scalary, 'k+', 'LineWidth', 1);
hold on
[score,scalary] = findSquares(tout,SpmX, dpSpmX3(1:end-1,:));
p7 = scatter(dpSpmX3(:,1),dpSpmX3(:,2).*scalary, 'k^', 'LineWidth', 1);
hold off
axis([tspan 0 max(SpmX)*1.2]);
% xlabel('time (min)','fontsize',LabelSize)
% ylabel('scaled concentration')
xlabel('(b)','fontsize',LabelSize)
% h = legend('SpmX_T','Kaczmarczyk\newlineet al. ''20','Location', 'North');
h = legend('SpmX','Radhakrishnan 07','Guo 04','Location', 'North','fontsize',LegendSize,'NumColumns',3);
  set(h, 'Interpreter', 'tex', 'Box', 'on', 'Orientation', 'vertical');
  
  
  
% axes(ha(3));
subplot(2,3,3)
DivJfT=sum(yout(151:160,:)+yout(171:180,:)+yout(191:200,:));
   DivJbT=sum(yout(161:170,:)+yout(181:190,:)+yout(201:210,:));
  
DivJT=DivJfT+DivJbT;
 
p5 = line(tout,DivJT , 'Color', 'k' , 'LineWidth', 2, 'Linestyle', '-');
p5 = line(tout,DivJfT , 'Color', 'r' , 'LineWidth', 2, 'Linestyle', '--');
p5 = line(tout,DivJbT , 'Color', 'r' , 'LineWidth', 2, 'Linestyle', ':');
hold on
[score,scalary] = findSquares(tout,DivJT, dpDivJ(1:end-1,:));
p7 = scatter(dpDivJ(:,1),dpDivJ(:,2).*scalary, 'k+', 'LineWidth', 1);
hold on
[score,scalary] = findSquares(tout,DivJT, dpDivJ2(1:end-1,:));
p7 = scatter(dpDivJ2(:,1),dpDivJ2(:,2).*scalary, 'k^', 'LineWidth', 1);
hold on
[score,scalary] = findSquares(tout,DivJT, dpDivJ3(1:end-1,:));
p7 = scatter(dpDivJ3(:,1),dpDivJ3(:,2).*scalary, 'k*', 'LineWidth', 1);
hold off
axis([tspan 0 max(DivJT)*1.8]);
% ylabel('scaled concentration')
xlabel('(c)','fontsize',LabelSize)
h = legend('DivJ_T','DivJ_f','DivJ_b','Wheeler 1999','Sanselicio 2015','Radhakrishnan 07','Location', 'North','fontsize',LegendSize,'NumColumns',3);
  set(h, 'Interpreter', 'tex', 'Box', 'on', 'Orientation', 'vertical');


  
%   axes(ha(4));
subplot(2,3,4)
PleCph=sum(yout(91:100,:)+yout(101:110,:)+yout(111:120,:)+yout(121:130,:));
PleCkin=sum(yout(131:140,:)+yout(141:150,:));

PleCT=PleCph+PleCkin;

p5 = line(tout,PleCT , 'Color', 'k' , 'LineWidth', 2, 'Linestyle', '-');
p5 = line(tout,PleCph , 'Color', 'r' , 'LineWidth', 2, 'Linestyle', '--');
p5 = line(tout,PleCkin , 'Color', 'r' , 'LineWidth', 2, 'Linestyle', ':');
hold on
[score,scalary] = findSquares(tout,PleCT, dpPleCtot(1:end-1,:));
p7 = scatter(dpPleCtot(:,1),dpPleCtot(:,2).*scalary, 'k+', 'LineWidth', 1);
hold on
[score,scalary] = findSquares(tout,PleCT, dpPleCtot2(1:end-1,:));
p7 = scatter(dpPleCtot2(:,1),dpPleCtot2(:,2).*scalary, 'k^', 'LineWidth', 1);
hold off
axis([tspan 0 max(PleCT)*1.6]);
xlabel('(d)','fontsize',LabelSize)
% xlabel('time (min)','fontsize',LabelSize)
ylabel('scaled concentration')
h = legend('PleC_T','PleC_{ph}','PleC_{kin}','Viollier 2002','Guo 04','Location', 'North','fontsize',LegendSize,'NumColumns',3);
 set(h, 'Interpreter', 'tex', 'Box', 'on', 'Orientation', 'vertical');
 
% axes(ha(5));
subplot(2,3,5)
DivK=sum(yout(221:230,:)+yout(171:180,:)+yout(181:190,:));
DivKP=sum(yout(231:240,:)+yout(191:200,:)+yout(201:210,:)+yout(111:120,:)+yout(121:130,:)...
    +yout(261:270,:)+yout(271:280,:));

DivKT=DivK+DivKP;
p5 = line(tout,DivKT , 'Color', 'k' , 'LineWidth', 2, 'Linestyle', '-');
p5 = line(tout,DivK, 'Color', 'b' , 'LineWidth', 2, 'Linestyle', '--');
p5 = line(tout,DivKP , 'Color', 'r' , 'LineWidth', 2, 'Linestyle', ':');
hold on
[score,scalary] = findSquares(tout,DivKT, dpDivKTOT2(1:end-1,:));
p7 = scatter(dpDivKTOT2(:,1),dpDivKTOT2(:,2).*scalary, 'k^', 'LineWidth', 1);
hold on
[score,scalary] = findSquares(tout,DivKT, dpDivKTOT3(1:end-1,:));
p7 = scatter(dpDivKTOT3(:,1),dpDivKTOT3(:,2).*scalary, 'k+', 'LineWidth', 1);
hold on
[score,scalary] = findSquares(tout,DivKP, dpDivKP(1:end-1,:));
p7 = scatter(dpDivKP(:,1),dpDivKP(:,2).*scalary, 'r^', 'LineWidth', 1);
hold off
axis([tspan 0 max(DivKT)*1.6]);
xlabel('(e)','fontsize',LabelSize)
% xlabel('time (min)','fontsize',LabelSize)
% ylabel('scaled concentration')
h = legend('DivK_T','DivK','DivK~P','Gao 04','Jacobs 00','Jacobs 00','Location', 'North','fontsize',LegendSize,'NumColumns',3);
 set(h, 'Interpreter', 'tex', 'Box', 'on', 'Orientation', 'vertical');
 
% axes(ha(6));
subplot(2,3,6)
CtrAP = sum(yout(81:90,:)+yout(321:330,:)+yout(331:340,:)); 
CtrA =sum(yout(71:80,:)+yout(341:350,:)+yout(351:360,:));
CtrAT=CtrA+CtrAP;

p5 = line(tout,CtrAT , 'Color', 'k' , 'LineWidth', 2, 'Linestyle', '-');
p5 = line(tout,CtrAP , 'Color', 'r' , 'LineWidth', 2, 'Linestyle', '-');
hold on
[score,scalary] = findSquares(tout,CtrAT, dpCtrAT(1:end-1,:));
p7 = scatter(dpCtrAT(:,1),dpCtrAT(:,2).*scalary, 'k+', 'LineWidth', 1);
% hold on
% [score,scalary] = findSquares(tout,sum(CtrAT), dpCtrAT2(1:end-1,:));
% p7 = scatter(dpCtrAT2(:,1),dpCtrAT2(:,2).*scalary, 'k^', 'LineWidth', 1);
hold on
[score,scalary] = findSquares(tout,CtrAT, dpCtrAT3(1:end-1,:));
p7 = scatter(dpCtrAT3(:,1),dpCtrAT3(:,2).*scalary, 'k^', 'LineWidth', 1);
hold on
[score,scalary] = findSquares(tout,CtrAP, dpCtrAP(1:end-1,:));
p7 = scatter(dpCtrAP(:,1),dpCtrAP(:,2).*scalary, 'r+', 'LineWidth', 1);
hold off
axis([tspan 0 max(CtrAT)*1.6]);
% xlabel('time (min)','fontsize',LabelSize)
% ylabel('scaled concentration')
xlabel('(f)','fontsize',LabelSize)
% h = legend('SpmX_T','Kaczmarczyk\newlineet al. ''20','Location', 'North');
% h = legend('CtrA_T','CtrA~P','Mcgraph','Collier et al. ''06','Jacobs et al. ''03','Location', 'North');
h = legend('CtrA_T','CtrA~P','Mcgraph','Radhakrishnan 07','Jacobs 03','Location', 'North','fontsize',LegendSize,'NumColumns',3);
set(h, 'Interpreter', 'tex', 'Box', 'on', 'Orientation', 'vertical');
    end



