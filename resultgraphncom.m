function resultgraphncom(yout,tout,TimeIni,TimeZring,celltype,mutant,TITLE,concentration)

% close all
%spatial plots
%%simulation
[PodJm PodJp PodJS SpmXm SpmXp PopZm PopZp CtrA CtrAP PleCf PleCb...
 PleCfDivKP PleCbDivKP PleCfkin PleCbkin DivJf DivJb DivJfDivK DivJbDivK...
 DivJfDivKP DivJbDivKP PerP DivK DivKP DivLf DivLb DivLfDivKP DivLbDivKP...
 CckAfkin CckAbkin CckAfph CckAbph CtrAPCckAfph CtrAPCckAbph CtrACckAfkin CtrACckAbkin...
 CpdRf CpdRb CpdRP N PoleNum n]=INPUT();    
    


FONTSIZE=14;  LegendSize=8;
    %% Scaffolding spatial and temporal
PodJLpS=yout((PodJp-1)*n+1:PodJp*n,:);
PodJLpS = flipud(PodJLpS);
PodJmS=yout((PodJm-1)*n+1:PodJm*n,:);
PodJmS = flipud(PodJmS);
PodJLS=PodJLpS+PodJmS;%%
PodJSS=yout((PodJS-1)*n+1:PodJS*n,:);
PodJSS = flipud(PodJSS);%%
PodJtot = PodJLS + PodJSS;%%


PopZpS = yout((PopZp-1)*n+1:PopZp*n,:);
PopZmS = yout((PopZm-1)*n+1:PopZm*n,:);
PopZmS = flipud(PopZmS);
PopZpS = flipud(PopZpS);
PopZtot = PopZpS+PopZmS;

SpmXtot=yout((SpmXm-1)*n+1:SpmXm*n,:)+yout((SpmXp-1)*n+1:SpmXp*n,:);
SpmXtot = flipud(SpmXtot);



DivJbT=yout((DivJb-1)*n+1:DivJb*n,:)+yout((DivJbDivK-1)*n+1:DivJbDivK*n,:)+yout((DivJbDivKP-1)*n+1:DivJbDivKP*n,:);
DivJfT=yout((DivJf-1)*n+1:DivJf*n,:)+yout((DivJfDivK-1)*n+1:DivJfDivK*n,:)+yout((DivJfDivKP-1)*n+1:DivJfDivKP*n,:);
DivJbT = flipud(DivJbT);
DivJfT = flipud(DivJfT);
DivJtot=DivJfT+DivJbT;



PleCfS=yout((PleCf-1)*n+1:PleCf*n,:);
PleCbS=yout((PleCb-1)*n+1:PleCb*n,:);
PleCfDivKPS=yout((PleCfDivKP-1)*n+1:PleCfDivKP*n,:);
PleCbDivKPS=yout((PleCbDivKP-1)*n+1:PleCbDivKP*n,:);
PleCfkinS=yout((PleCfkin-1)*n+1:PleCfkin*n,:);
PleCbkinS=yout((PleCbkin-1)*n+1:PleCbkin*n,:);

PleCphT=PleCfS+PleCbS+PleCfDivKPS+PleCbDivKPS;%%%%%%%%%%%%%%%%???
PleCkinT=PleCfkinS+PleCbkinS;
PleCphT=flipud(PleCphT);
PleCkinT=flipud(PleCkinT);
PleCtot=PleCkinT+PleCphT;

DivJfDivKS = yout((DivJfDivK-1)*n+1:DivJfDivK*n,:);
DivJbDivKS = yout((DivJbDivK-1)*n+1:DivJbDivK*n,:);
DivJfDivKPS = yout((DivJfDivKP-1)*n+1:DivJfDivKP*n,:);
DivJbDivKPS = yout((DivJbDivKP-1)*n+1:DivJbDivKP*n,:);
DivJfDivKS = flipud(DivJfDivKS);
DivJbDivKS = flipud(DivJbDivKS);
DivJfDivKPS = flipud(DivJfDivKPS);
DivJbDivKPS = flipud(DivJbDivKPS);


DivKS=yout((DivK-1)*n+1:DivK*n,:);
DivKS = flipud(DivKS);
DivKPS=yout((DivKP-1)*n+1:DivKP*n,:);
DivKPS = flipud(DivKPS);

DivJfS = yout((DivJf-1)*n+1:DivJf*n,:);
DivJbS = yout((DivJb-1)*n+1:DivJb*n,:);
DivJfS = flipud(DivJfS);
DivJbS = flipud(DivJbS);

DivLfDivKPS = yout((DivLfDivKP-1)*n+1:DivLfDivKP*n,:);
DivLbDivKPS = yout((DivLbDivKP-1)*n+1:DivLbDivKP*n,:);
DivLfDivKPS = flipud(DivLfDivKPS);
DivLbDivKPS = flipud(DivLbDivKPS);

DivLfS=yout((DivLf-1)*n+1:DivLf*n,:);
DivLbS = yout((DivLb-1)*n+1:DivLb*n,:);
DivLfS = flipud(DivLfS);
DivLbS = flipud(DivLbS);

DivLT=DivLfS+DivLbS;
DivLDivKPT=DivLfDivKPS+DivLbDivKPS;
DivLtot = DivLT+DivLDivKPT;%%%%

DivKT = DivKS + DivJfDivKS + DivJbDivKS;
DivKPT = DivKPS + DivJfDivKPS + DivJbDivKPS + DivLfDivKPS + DivLbDivKPS;
DivKtot = DivKT+DivKPT;%%%%

DivJfT = DivJfS + DivJfDivKS + DivJfDivKPS;
DivJbT = DivJbS + DivJbDivKS + DivJbDivKPS;
DivJtot = DivJfT + DivJbT;%%%%%%%%%%%%

CtrACckAkin=yout((CtrACckAfkin-1)*n+1:CtrACckAfkin*n,:)+yout((CtrACckAbkin-1)*n+1:CtrACckAbkin*n,:);
CtrAPCckAph=yout((CtrAPCckAfph-1)*n+1:CtrAPCckAfph*n,:)+yout((CtrAPCckAbph-1)*n+1:CtrAPCckAbph*n,:);
CtrAPS=yout((CtrAP-1)*n+1:CtrAP*n,:);
CtrAS=yout((CtrA-1)*n+1:CtrA*n,:);
CtrACckAkin = flipud(CtrACckAkin);
CtrAPCckAph = flipud(CtrAPCckAph);
CtrAPS = flipud(CtrAPS);
CtrAS = flipud(CtrAS);
CtrAT = CtrAS + CtrACckAkin;
CtrAPT = CtrAPS + CtrAPCckAph;
CtrAtot = CtrAT+CtrAPT;

PerPS=yout((PerP-1)*n+1:PerP*n,:);
PerPS=flipud(PerPS);


CckAfkinS=yout((CckAfkin-1)*n+1:CckAfkin*n,:);
CckAbkinS=yout((CckAbkin-1)*n+1:CckAbkin*n,:);
CckAfphS=yout((CckAfph-1)*n+1:CckAfph*n,:);
CckAbphS=yout((CckAbph-1)*n+1:CckAbph*n,:);
CckAfkinS=flipud(CckAfkinS);
CckAbkinS=flipud(CckAbkinS);
CckAfphS=flipud(CckAfphS);
CckAbphS=flipud(CckAbphS);

CckAphT=CckAfphS+CckAbphS;
CckAkinT=CckAfkinS+CckAbkinS;
CckAT=CckAkinT+CckAphT;
CckAtot = CckAT + CtrAPCckAph + CtrACckAkin;
CckAphtot = CckAphT + CtrAPCckAph;
CckAkintot = CckAkinT + CtrACckAkin;

CpdRu=yout((CpdRf-1)*n+1:CpdRf*n,:)+yout((CpdRb-1)*n+1:CpdRb*n,:);
CpdRPS=yout((CpdRP-1)*n+1:CpdRP*n,:);
CpdRu=flipud(CpdRu); CpdRPS=flipud(CpdRPS);
CpdRtot = CpdRu + CpdRPS;

% DivKGradient = DivKtot(:,GradientCheckIndex1:GradientCheckIndex2);%%%%%%%%
% DivKPGradient = DivKPT(:,GradientCheckIndex1:GradientCheckIndex2);
% CtrAGradient = CtrAtot(:,GradientCheckIndex1:GradientCheckIndex2);%%%%%%%%
% CtrAPGradient = CtrAPT(:,GradientCheckIndex1:GradientCheckIndex2);%%%%%%%%
% CtrAuGradient = CtrAT(:,GradientCheckIndex1:GradientCheckIndex2);%%%%%%%%
%%Mgrid
Y=yout';
L=length(Y(:,1));
 M(:,n/2+1)=zeros(1,L);%midcell
count=1;
for ii=n/2:-1:1
    M(:,ii)=-count*Y(:,N);
    count=count+1;
end

M(:,(n/2+2):(n+1))=-fliplr(M(:,1:n/2));

M=M';
a=zeros(1,L);
time=tout;



%% plots
% figure()
% subplot(2,4,1)
% CckAkintot(n+1,:)=a;
% pcolor(time, M, CckAkintot)
% shading flat
% colorbar
% % caxis([0 4]);
% title('total CckA kinase')
% subplot(2,4,2)
% CckAphtot(n+1,:)=a;
% pcolor(time, M, CckAphtot)
% shading flat
% colorbar
% % caxis([0 4]);
% title('total CckA phosphatase')
% subplot(2,4,3)
% DivLT(n+1,:)=a;
% pcolor(time, M, DivLT)
% shading flat
% colorbar
% % caxis([0 4]);
% title('DivLT')
% subplot(2,4,4)
% DivLDivKPT(n+1,:)=a;
% pcolor(time, M, DivLDivKPT)
% shading flat
% colorbar
% % caxis([0 4]);
% title('DivLDivKPT')
% subplot(2,4,5)
% DivKPT(n+1,:)=a;
% pcolor(time, M, DivKPT)
% shading flat
% colorbar
% % caxis([0 4]);
% title('DivKPT')
% subplot(2,4,6)
% DivKT(n+1,:)=a;
% pcolor(time, M, DivKT)
% shading flat
% colorbar
% % caxis([0 4]);
% title('DivKuT')
% set(findall(gcf,'-property','FontSize'),'FontSize',14)

%% Spatial Plots WT
if strcmp(mutant,'WT')
    
NN=length(TimeIni);
XX1=TimeIni;
XX2=TimeZring;
YY=zeros(1,NN)-3;
STR1={};
for i=1:NN
    STR1{end+1}='\uparrow Rep Ini';
end
STR2={};
for i=1:NN
    STR2{end+1}='\uparrow Z-ring closed';
end

%%%%%%%%%%%%%%%%%%%%%%%%Spatial Scaffolding Proteins%%%%%%%%%%%%%%%%%%%%%%%
figure()
set(gcf,'position',[100 100 1800 300])%left, lower, right, upper
set(gcf,'Name',TITLE);
ha = tight_subplot(1,4,[.1 .03],[.2 .1],[.05 .04]);
axes(ha(1)); 
PopZtot(n+1,:)=a;
pcolor(time, M, PopZtot)
shading flat
colorbar
% caxis([0 4]);
title('(a) total PopZ','fontweight','bold','fontsize',FONTSIZE)
line([TimeIni,TimeIni],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)
line([TimeZring,TimeZring],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)


axes(ha(2)); 
SpmXtot(n+1,:)=a;
pcolor(time, M, SpmXtot)
shading flat
colorbar
% caxis([0 0.25]);
title('(b) total SpmX','fontweight','bold','fontsize',FONTSIZE)
set(findall(gcf,'-property','FontSize'),'FontSize',14)
line([TimeIni,TimeIni],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)
line([TimeZring,TimeZring],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)

axes(ha(3)); 
PodJLS(n+1,:)=a;
pcolor(time, M, PodJLS)
shading flat
colorbar
caxis([0 2]);
title('(c) total long form PodJ','fontweight','bold','fontsize',FONTSIZE)
set(findall(gcf,'-property','FontSize'),'FontSize',14)
line([TimeIni,TimeIni],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)
line([TimeZring,TimeZring],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)

axes(ha(4));
PodJSS(n+1,:)=a;
pcolor(time, M, PodJSS)
shading flat
colorbar
caxis([0 2]);
% xlabel('time (min)')
title('(d) short form PodJ','fontweight','bold','fontsize',FONTSIZE)
set(findall(gcf,'-property','FontSize'),'FontSize',14)
line([TimeIni,TimeIni],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)
line([TimeZring,TimeZring],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)

%%%%%%%%%%%%%%%%%%%%%%%%Spatial Client Proteins%%%%%%%%%%%%%%%%%%%%%%%
% figure()
% set(gcf,'position',[100 100 2500 800])%left, lower, right, upper
% set(gcf,'Name',TITLE);
% ha = tight_subplot(2,4,[.1 .03],[.1 .1],[.05 .04]);
% 
% axes(ha(1)); 
% PleCtot(n+1,:)=a;
% pcolor(time, M, PleCtot)
% shading flat
% colorbar
% % caxis([0 0.2]);
% title('(a) total PleC')
% set(findall(gcf,'-property','FontSize'),'FontSize',14)
% line([TimeIni,TimeIni],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)
% line([TimeZring,TimeZring],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)
% 
% axes(ha(2));
% DivJtot(n+1,:)=a;
% pcolor(time, M, DivJtot)
% shading flat
% colorbar
% % caxis([0.1 0.4]);
% title('(b) total DivJ')
% set(findall(gcf,'-property','FontSize'),'FontSize',14)
% line([TimeIni,TimeIni],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)
% line([TimeZring,TimeZring],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)
% 
% axes(ha(3)); 
% DivKPT(n+1,:)=a;
% pcolor(time, M, DivKPT)
% shading flat
% colorbar
% % caxis([0.1 0.35]);
% title('(c) total phosphorylated DivK')
% set(findall(gcf,'-property','FontSize'),'FontSize',14)
% line([TimeIni,TimeIni],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)
% line([TimeZring,TimeZring],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)
% 
% 
% axes(ha(4)); 
% DivLT(n+1,:)=a;
% pcolor(time, M, DivLT)
% shading flat
% colorbar
% % caxis([0.002 0.017]);
% title('(d) DivL not bound to DivKP')
% set(findall(gcf,'-property','FontSize'),'FontSize',14)
% line([TimeIni,TimeIni],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)
% line([TimeZring,TimeZring],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)
% 
% 
% axes(ha(5));
% DivLDivKPT(n+1,:)=a;
% pcolor(time, M, DivLDivKPT)
% shading flat
% colorbar
% % caxis([0.002 0.017]);
% title('(e) DivL bound to DivKP')
% set(findall(gcf,'-property','FontSize'),'FontSize',14)
% line([TimeIni,TimeIni],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)
% line([TimeZring,TimeZring],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)
% 
% % axes(ha(6));
% % CckAtot(n+1,:)=a;
% % pcolor(time, M, CckAtot)
% % shading flat
% % colorbar
% % % caxis([0.05 0.18]);
% % title('(f) total CckA')
% % set(findall(gcf,'-property','FontSize'),'FontSize',14)
% % line([TimeIni,TimeIni],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)
% % line([TimeZring,TimeZring],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)
% 
% axes(ha(6));
% CckAkinT(n+1,:)=a;
% pcolor(time, M, CckAkinT)
% shading flat
% colorbar
% % caxis([0.05 0.18]);
% title('(f) CckA kinase')
% set(findall(gcf,'-property','FontSize'),'FontSize',14)
% line([TimeIni,TimeIni],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)
% line([TimeZring,TimeZring],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)
% 
% axes(ha(7)); 
% CpdRu(n+1,:)=a;
% pcolor(time, M, CpdRu)
% shading flat
% colorbar
% % caxis([0.6 1]);
% title('(g) unphosphorylated CpdR')
% set(findall(gcf,'-property','FontSize'),'FontSize',14)
% line([TimeIni,TimeIni],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)
% line([TimeZring,TimeZring],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)
% 
% axes(ha(8));
% CtrAPT(n+1,:)=a;
% pcolor(time, M, CtrAPT)
% shading flat
% colorbar
% % caxis([0.05 0.55]);
% title('(h) total CtrAP')
% line([TimeIni,TimeIni],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)
% line([TimeZring,TimeZring],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)
% set(findall(gcf,'-property','FontSize'),'FontSize',14)

%%%%%%%%%%%%%%%%%%%%%%%%Spatial Client Proteins 2%%%%%%%%%%%%%%%%%%%%%%%
figure()
set(gcf,'position',[100 100 2500 800])%left, lower, right, upper
set(gcf,'Name',TITLE);
ha = tight_subplot(2,4,[.1 .03],[.1 .1],[.05 .04]);

axes(ha(1)); 
PleCtot(n+1,:)=a;
pcolor(time, M, PleCtot)
shading flat
colorbar
% caxis([0 0.2]);
title('(a) total PleC','fontweight','bold','fontsize',FONTSIZE)
set(findall(gcf,'-property','FontSize'),'FontSize',14)
line([TimeIni,TimeIni],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)
line([TimeZring,TimeZring],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)

axes(ha(2));
DivJtot(n+1,:)=a;
pcolor(time, M, DivJtot)
shading flat
colorbar
% caxis([0.1 0.4]);
title('(b) total DivJ','fontweight','bold','fontsize',FONTSIZE)
set(findall(gcf,'-property','FontSize'),'FontSize',14)
line([TimeIni,TimeIni],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)
line([TimeZring,TimeZring],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)

axes(ha(3));
DivKPT(n+1,:)=a;
pcolor(time, M, DivKPT)
shading flat
colorbar
% caxis([0.1 0.35]);
title('(c) total phosphorylated DivK','fontweight','bold','fontsize',FONTSIZE)
set(findall(gcf,'-property','FontSize'),'FontSize',14)
line([TimeIni,TimeIni],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)
line([TimeZring,TimeZring],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)

axes(ha(4)); 
DivLtot(n+1,:)=a;
pcolor(time, M, DivLtot)
shading flat
colorbar
% caxis([0.002 0.017]);
title('(d) total DivL','fontweight','bold','fontsize',FONTSIZE)
set(findall(gcf,'-property','FontSize'),'FontSize',14)
line([TimeIni,TimeIni],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)
line([TimeZring,TimeZring],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)

axes(ha(5));
CckAtot(n+1,:)=a;
pcolor(time, M, CckAtot)
shading flat
colorbar
% caxis([0.05 0.18]);
title('(e) total CckA','fontweight','bold','fontsize',FONTSIZE)
set(findall(gcf,'-property','FontSize'),'FontSize',14)
line([TimeIni,TimeIni],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)
line([TimeZring,TimeZring],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)





axes(ha(6)); 
CpdRu(n+1,:)=a;
pcolor(time, M, CpdRu)
shading flat
colorbar
% caxis([0.6 1]);
title('(f) unphosphorylated CpdR','fontweight','bold','fontsize',FONTSIZE)
set(findall(gcf,'-property','FontSize'),'FontSize',14)
line([TimeIni,TimeIni],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)
line([TimeZring,TimeZring],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)

axes(ha(7));
CtrAtot(n+1,:)=a;
pcolor(time, M, CtrAtot)
shading flat
colorbar
% caxis([0 2]);
title('(g) total CtrA','fontweight','bold','fontsize',FONTSIZE)
set(findall(gcf,'-property','FontSize'),'FontSize',14)
line([TimeIni,TimeIni],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)
line([TimeZring,TimeZring],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)


axes(ha(8));
CtrAPT(n+1,:)=a;
pcolor(time, M, CtrAPT)
shading flat
colorbar
% caxis([0.05 0.55]);
title('(h) total CtrAP','fontweight','bold','fontsize',FONTSIZE)
line([TimeIni,TimeIni],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)
line([TimeZring,TimeZring],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)
set(findall(gcf,'-property','FontSize'),'FontSize',14)



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fig S1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
set(gcf,'position',[100 100 1800 1200])%left, lower, right, upper
set(gcf,'Name',TITLE);
ha = tight_subplot(3,3,[.1 .03],[.1 .05],[.05 .04]);
axes(ha(1)); 
PleCkinT(n+1,:)=a;
pcolor(time, M, PleCkinT)
shading flat
colorbar
% caxis([0 4]);
title('(a) PleC_{kinase}','fontweight','bold','fontsize',FONTSIZE)
line([TimeIni,TimeIni],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)
line([TimeZring,TimeZring],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)

axes(ha(2)); 
PleCphT(n+1,:)=a;
pcolor(time, M, PleCphT)
shading flat
colorbar
% caxis([0 4]);
title('(b) PleC_{phosphatase}','fontweight','bold','fontsize',FONTSIZE)
line([TimeIni,TimeIni],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)
line([TimeZring,TimeZring],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)

axes(ha(3)); 
DivLT(n+1,:)=a;
pcolor(time, M, DivLT)
shading flat
colorbar
% caxis([0.002 0.017]);
title('(c) DivL not bound to DivKP','fontweight','bold','fontsize',FONTSIZE)
set(findall(gcf,'-property','FontSize'),'FontSize',14)
line([TimeIni,TimeIni],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)
line([TimeZring,TimeZring],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)


axes(ha(4));
DivLDivKPT(n+1,:)=a;
pcolor(time, M, DivLDivKPT)
shading flat
colorbar
% caxis([0.002 0.017]);
title('(d) DivL bound to DivKP','fontweight','bold','fontsize',FONTSIZE)
set(findall(gcf,'-property','FontSize'),'FontSize',14)
line([TimeIni,TimeIni],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)
line([TimeZring,TimeZring],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)

axes(ha(5));
CckAkintot(n+1,:)=a;
pcolor(time, M, CckAkintot)
shading flat
colorbar
% caxis([0.05 0.18]);
title('(e) CckA_{kinase}','fontweight','bold','fontsize',FONTSIZE)
set(findall(gcf,'-property','FontSize'),'FontSize',14)
line([TimeIni,TimeIni],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)
line([TimeZring,TimeZring],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)

axes(ha(6));
CckAphtot(n+1,:)=a;
pcolor(time, M, CckAphtot)
shading flat
colorbar
% caxis([0.05 0.18]);
title('(f) CckA_{phosphatase}','fontweight','bold','fontsize',FONTSIZE)
set(findall(gcf,'-property','FontSize'),'FontSize',14)
line([TimeIni,TimeIni],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)
line([TimeZring,TimeZring],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)



axes(ha(7)); 
CpdRtot(n+1,:)=a;
pcolor(time, M, CpdRtot)
shading flat
colorbar
% caxis([0 0.25]);
title('(g) total CpdR','fontweight','bold','fontsize',FONTSIZE)
set(findall(gcf,'-property','FontSize'),'FontSize',14)
line([TimeIni,TimeIni],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)
line([TimeZring,TimeZring],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)

axes(ha(8));
CpdRPS(n+1,:)=a;
pcolor(time, M, CpdRPS)
shading flat
colorbar
% caxis([0 4]);
title('(h) phosphorylated CpdR','fontweight','bold','fontsize',FONTSIZE)
line([TimeIni,TimeIni],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)
line([TimeZring,TimeZring],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)

axes(ha(9)); 
CtrAT(n+1,:)=a;
pcolor(time, M, CtrAT)
shading flat
colorbar
% caxis([0 2]);
title('(i) unphosphorylated CtrA','fontweight','bold','fontsize',FONTSIZE)
set(findall(gcf,'-property','FontSize'),'FontSize',14)
line([TimeIni,TimeIni],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)
line([TimeZring,TimeZring],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)



% %%%%%%%%%%%%%%%%%%%%%%%% background 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
 set(gcf,'position',[100 100 1800 300])%left, lower, right, upper
ha = tight_subplot(1,4,[.1 .03],[.1 .1],[.05 .04]);
axes(ha(1)); 
 PodJLS(n+1,:)=a;
pcolor(time, M, PodJLS)
shading flat
colorbar
caxis([0 4]);
title('PodJL','fontweight','bold','fontsize',FONTSIZE)
axes(ha(2)); 
    PopZtot(n+1,:)=a;
pcolor(time, M, PopZtot)
shading flat
colorbar
caxis([0 30]);
title('PopZ','fontweight','bold','fontsize',FONTSIZE)

axes(ha(3)); 
SpmXtot(n+1,:)=a;
pcolor(time, M, SpmXtot)
shading flat
colorbar
caxis([0 3]);
title('SpmX','fontweight','bold','fontsize',FONTSIZE)

axes(ha(4)); 
DivJtot(n+1,:)=a;
pcolor(time, M, DivJtot)
shading flat
colorbar
% caxis([0.1 0.35]);
title('DivJ','fontweight','bold','fontsize',FONTSIZE)
set(findall(gcf,'-property','FontSize'),'FontSize',14)
% %%%%%%%%%%%%%%%%%%%%%%%% background 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure()
%     set(gcf,'position',[100 100 400 500])%left, lower, right, upper
% subplot(2,1,1)
%  DivKPT(n+1,:) = a;
% pcolor(time, M, DivKPT)
% shading flat
% colorbar
% caxis([0 0.9]);
% title(['(a) WT' ],'fontweight','bold','fontsize',FONTSIZE)
% ylabel('DivKP','fontweight','bold','fontsize',FONTSIZE)
% subplot(2,1,2)
% CtrAPT(n+1,:)=a;
% pcolor(time, M, CtrAPT)
% shading flat
% colorbar
% % caxis([0 0.1]);
% ylabel('CtrAP','fontweight','bold','fontsize',FONTSIZE)
% 
%%%%%%%%%%%%%%%%%%%%%%PodJtot%%%%%%%%%%%%%%%%%%%%%
figure()
PodJtot(n+1,:)=a;
pcolor(time, M, PodJtot)
shading flat
colorbar
caxis([0 4]);
% xlabel('time (min)')
title('total PodJ')
set(findall(gcf,'-property','FontSize'),'FontSize',14)
line([TimeIni,TimeIni],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)
line([TimeZring,TimeZring],[-2.2,2.2],'Color','red','LineStyle','--','LineWidth',2)





%%%%%%%%%%%%%%%%% gradient
% figure()
% ax = gca;
% % set(fig,'defaultAxesColorOrder',['k'; 'r']);
% ComIndex = 1:1:n;
% % Xlabel = categorical({'1','2','3','4','5','6','7','8','9','10'});
% % Xlabel = reordercats(Xlabel,{'1','2','3','4','5','6','7','8','9','10'});
% TRAPZ_X = time(GradientCheckIndex1:GradientCheckIndex2);
% Yleft = [];
% for i=1:n
%     Yleft = [Yleft trapz(TRAPZ_X',DivKPGradient(i,:))/(TRAPZ_X(end)-TRAPZ_X(1))];
% end
% % Yleft=sum(DivKPGradient,2)/timespan;
% Yleft = fliplr (Yleft);
% % Yright=sum(CtrAPGradient,2)/timespan;
% Yright1 = [];
% Yright2 = [];
% for i=1:n
%     Yright1 = [Yright1 trapz(TRAPZ_X',CtrAPGradient(i,:))/(TRAPZ_X(end)-TRAPZ_X(1))];
%     Yright2 = [Yright2 trapz(TRAPZ_X',CtrAuGradient(i,:))/(TRAPZ_X(end)-TRAPZ_X(1))];
% end
% Yright1 = fliplr (Yright1);
% yyaxis left
% plot(ComIndex,Yleft, 'Color', 'k' , 'LineWidth', 2)
% ylabel('Average Scaled Concentration','fontsize',14)
% ax.YColor = 'k';
% yyaxis right
% plot(ComIndex,Yright1, 'Color', 'r' , 'LineWidth', 2)
% % hold on
% % plot(ComIndex,Yright2, 'Color', 'b' , 'LineWidth', 2)
% h = legend('total DivK~P','total CtrA~P');
% ax.YColor = 'r';
% xlabel('Compartment Index (left1: new pole)','fontsize',14)
% xlim([1 n]);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%% mutant case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mutant,'deltaPodJ')   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% deltaPodJ
    %%%%%%%%%%%%%%%%%%%%Fig S2
    figure()
    set(gcf,'position',[100 100 400 500])%left, lower, right, upper
subplot(2,1,1)
 DivKPT(n+1,:) = a;
pcolor(time, M, DivKPT)
shading flat
colorbar
% caxis([0 0.35]);
title('\Delta {\itpodJ}','fontweight','bold','fontsize',FONTSIZE)
ylabel('DivKP','fontweight','bold','fontsize',FONTSIZE)
subplot(2,1,2)
CtrAPT(n+1,:)=a;
pcolor(time, M, CtrAPT)
shading flat
colorbar
% caxis([0 0.1]);
ylabel('CtrAP','fontweight','bold','fontsize',FONTSIZE)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% M1
  
    figure()
 set(gcf,'position',[100 100 1800 300])%left, lower, right, upper
ha = tight_subplot(1,4,[.1 .03],[.1 .1],[.05 .04]);

axes(ha(1)); 
    PopZtot(n+1,:)=a;
pcolor(time, M, PopZtot)
shading flat
colorbar
% caxis([0 4]);
title('total PopZ')

axes(ha(2)); 
SpmXtot(n+1,:)=a;
pcolor(time, M, SpmXtot)
shading flat
colorbar
% caxis([0 0.25]);
title('total SpmX')

axes(ha(3)); 
DivJtot(n+1,:)=a;
pcolor(time, M, DivJtot)
shading flat
colorbar
% caxis([0.1 0.35]);
title('total DivJ')

axes(ha(4)); 
PleCtot(n+1,:)=a;
pcolor(time, M, PleCtot)
shading flat
colorbar
% caxis([0 0.1]);
xlabel('time (min)')
title('total PleC')
set(findall(gcf,'-property','FontSize'),'FontSize',14)

elseif strcmp(mutant,'deltaSpmX')      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% deltaSpmX
   %%%%%%%%%%%%%%%%%%%%Fig S2
    figure()
    set(gcf,'position',[100 100 400 500])%left, lower, right, upper
subplot(2,1,1)
 DivKPT(n+1,:) = a;
pcolor(time, M, DivKPT)
shading flat
colorbar
% caxis([0 0.35]);
title('\Delta {\itspmX}','fontweight','bold','fontsize',FONTSIZE)
ylabel('DivKP','fontweight','bold','fontsize',FONTSIZE)
subplot(2,1,2)
CtrAPT(n+1,:)=a;
pcolor(time, M, CtrAPT)
shading flat
colorbar
% caxis([0 0.1]);
ylabel('CtrAP','fontweight','bold','fontsize',FONTSIZE)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% M1

 figure()
 set(gcf,'position',[100 100 1800 300])%left, lower, right, upper
ha = tight_subplot(1,4,[.1 .03],[.1 .1],[.05 .04]);
        axes(ha(1)); 
        PodJLS(n+1,:)=a;
pcolor(time, M, PodJLS)
shading flat
colorbar
caxis([0 4]);
title('PodJL')

axes(ha(2)); 
PodJSS(n+1,:)=a;
pcolor(time, M, PodJSS)
shading flat
colorbar
% caxis([0 0.5]);
title('PodJS')

axes(ha(3)); 
PopZtot(n+1,:)=a;
pcolor(time, M, PopZtot)
shading flat
colorbar
% caxis([0 4]);
title('total PopZ')

axes(ha(4)); 
DivJtot(n+1,:)=a;
pcolor(time, M, DivJtot)
shading flat
colorbar
% caxis([0.1/2.5 0.35/2.5]);
title('total DivJ')
set(findall(gcf,'-property','FontSize'),'FontSize',14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% deltaPopZ
elseif strcmp(mutant,'deltaPopZ')
    %%%%%%%%%%%%%%%%%%%%Fig S2
    figure()
    set(gcf,'position',[100 100 400 500])%left, lower, right, upper
subplot(2,1,1)
 DivKPT(n+1,:) = a;
pcolor(time, M, DivKPT)
shading flat
colorbar
% caxis([0 0.35]);
title('\Delta {\itpopZ}','fontweight','bold','fontsize',FONTSIZE)
ylabel('DivKP','fontweight','bold','fontsize',FONTSIZE)
subplot(2,1,2)
CtrAPT(n+1,:)=a;
pcolor(time, M, CtrAPT)
shading flat
colorbar
% caxis([0 0.1]);
ylabel('CtrAP','fontweight','bold','fontsize',FONTSIZE)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% M1
    figure()
 set(gcf,'position',[100 100 1800 300])%left, lower, right, upper
ha = tight_subplot(1,4,[.1 .03],[.1 .1],[.05 .04]);
        axes(ha(1)); 
        PodJtot(n+1,:)=a;
pcolor(time, M, PodJtot)
shading flat
colorbar
% caxis([0.1 0.23]);
title('total PodJ')

        axes(ha(2)); 
        SpmXtot(n+1,:)=a;
pcolor(time, M, SpmXtot)
shading flat
colorbar
% caxis([0 0.2]);
title('total SpmX')

        axes(ha(3)); 
DivJtot(n+1,:)=a;
pcolor(time, M, DivJtot)
shading flat
colorbar
caxis([0.08 0.18]);
% xlabel('time (min)')
title('total DivJ')

        axes(ha(4)); 
        CpdRu(n+1,:)=a;
pcolor(time, M, CpdRu)
shading flat
colorbar
caxis([1.5 3.5]);
% xlabel('time (min)')
title('unphosphorylated CpdR')
set(findall(gcf,'-property','FontSize'),'FontSize',14)

elseif strcmp(mutant,'PleC-F778L')%%%%%%%%%%%%%%%%%%%% 
  figure()
    set(gcf,'position',[100 100 400 500])%left, lower, right, upper
subplot(2,1,1)
 DivKPT(n+1,:) = a;
pcolor(time, M, DivKPT)
shading flat
colorbar
caxis([0 0.9]);
title(['(f) {\itpleC-F778L}' ],'fontweight','bold','fontsize',FONTSIZE)
ylabel('DivKP','fontweight','bold','fontsize',FONTSIZE)
subplot(2,1,2)
CtrAPT(n+1,:)=a;
pcolor(time, M, CtrAPT)
shading flat
colorbar
% caxis([0 0.1]);
ylabel('CtrAP','fontweight','bold','fontsize',FONTSIZE)

elseif strcmp(mutant,'p3:deletingPodJPleCbinding')%%%%%%%%%%%%%%%%%%%%
  figure()
    set(gcf,'position',[100 100 400 500])%left, lower, right, upper
subplot(2,1,1)
 DivKPT(n+1,:) = a;
pcolor(time, M, DivKPT)
shading flat
colorbar
caxis([0 0.9]);
title('(g) delocalized PleC','fontweight','bold','fontsize',FONTSIZE)
ylabel('DivKP','fontweight','bold','fontsize',FONTSIZE)
subplot(2,1,2)
CtrAPtot(n+1,:)=a;
pcolor(time, M, CtrAPtot)
shading flat
colorbar
% caxis([0 0.1]);
ylabel('CtrAP','fontweight','bold','fontsize',FONTSIZE)
elseif strcmp(mutant,'p4:deletingPodJDivLbinding')%%%%%%%%%%%%%%%%%%%%
  figure()
    set(gcf,'position',[100 100 400 500])%left, lower, right, upper
subplot(2,1,1)
 DivKPT(n+1,:) = a;
pcolor(time, M, DivKPT)
shading flat
colorbar
caxis([0 0.9]);
title('(h) delocalized DivL','fontweight','bold','fontsize',FONTSIZE)
ylabel('DivKP','fontweight','bold','fontsize',FONTSIZE)
subplot(2,1,2)
CtrAPtot(n+1,:)=a;
pcolor(time, M, CtrAPtot)
shading flat
colorbar
% caxis([0 0.1]);
ylabel('CtrAP','fontweight','bold','fontsize',FONTSIZE)


elseif strcmp(mutant,'DivJ-H338A')
      figure()
    set(gcf,'position',[100 100 400 500])%left, lower, right, upper
subplot(2,1,1)
DivKtot(n+1,:) = a;
        pcolor(time, M, DivKtot)
shading flat
colorbar
% caxis([8 11]);
ylabel('DivKtot','fontweight','bold','fontsize',FONTSIZE)
title(['(c) {\itdivJ-H338A}'],'fontweight','bold','fontsize',FONTSIZE)
subplot(2,1,2)
CtrAPT(n+1,:)=a;
pcolor(time, M, CtrAPT)
shading flat
colorbar
% caxis([6 8]);
ylabel('CtrAP','fontweight','bold','fontsize',FONTSIZE)
    

elseif strcmp(mutant,'deltaPleC')    
    figure()
    set(gcf,'position',[100 100 400 500])%left, lower, right, upper
subplot(2,1,1)
 DivKPT(n+1,:) = a;
pcolor(time, M, DivKPT)
shading flat
colorbar
% caxis([0 1.2]);
title('(d) \Delta {\itpleC}','fontweight','bold','fontsize',FONTSIZE)
ylabel('DivKP','fontweight','bold','fontsize',FONTSIZE)
subplot(2,1,2)
CtrAPT(n+1,:)=a;
pcolor(time, M, CtrAPT)
shading flat
colorbar
% caxis([0 0.1]);
ylabel('CtrAP','fontweight','bold','fontsize',FONTSIZE)


elseif strcmp(mutant,'PleC-H610A')%%%%%%%%%%%%%%%%%%%%
        figure()
    set(gcf,'position',[100 100 400 500])%left, lower, right, upper
subplot(2,1,1)
 DivKPT(n+1,:) = a;
pcolor(time, M, DivKPT)
shading flat
colorbar
% caxis([0 0.35]);
title(['(e) {\itpleC-H610A}'],'fontweight','bold','fontsize',FONTSIZE)
ylabel('DivKP','fontweight','bold','fontsize',FONTSIZE)
subplot(2,1,2)
CtrAPT(n+1,:)=a;
pcolor(time, M, CtrAPT)
shading flat
colorbar
% caxis([0 0.1]);
ylabel('CtrAP','fontweight','bold','fontsize',FONTSIZE)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%
tspan=[tout(1) tout(end)];
%% temperal DATA
SpmXTem = [0 2340.397; 20 11034.71; 40 13828.5; 60 14676.86; 80 13635.98; 100 9316.033; 120 8849.205; 140 12085.38; 160 13972.28];
SpmXTem(:,2) = SpmXTem(:,2)/max(SpmXTem(:,2));
dpSpmX=SpmXTem;%paper152
SpmX2 = [0 2857.246; 20 10014.882; 40 13874.083;60 13516.276;80 13973.64; 100 8493.619; 120 8413.962; 140 11470.74;160 11327.459];
SpmX2(:,2) = SpmX2(:,2)/max(SpmX2(:,2));
dpSpmX2=SpmX2;%paper118
SpmX3 =[0 1506.527;20 4945.134;40 8203.305;60 9740.79;80 9703.205;100 10964.376;120 10477.861;140 9446.033;160 4349.246];%paper138
SpmX3(:,2) = SpmX3(:,2)/max(SpmX3(:,2));
dpSpmX3=SpmX3;
PodJLTem=[0 0; 20 0; 40 1726.012; 60 12068.518; 80 23215.326; 100 22996.104; 120  9466.861; 140  3336.154; 160 5506.548];
PodJLTem(:,2) = PodJLTem(:,2)/max(PodJLTem(:,2));
dpPodJL=PodJLTem;%paper116
PodJL2=[0 790.719;20 1297.841;40 4612.861;60 9123.397;80 11505.154;100 13792.497;120 15360.276;140 8009.589;160 5775.439];
PodJL2(:,2) = PodJL2(:,2)/max(PodJL2(:,2));
dpPodJL2=PodJL2;%138
PodJSTem=[0  16869.841; 20  16960.447; 40  4387.841; 60  434.364;  80  0; 100  2627.527; 120  8780.79; 140  8984.033; 160  6333.598];
PodJSTem(:,2) = PodJSTem(:,2)/max(PodJSTem(:,2));
dpPodJS=PodJSTem;%paper116
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
DivKPTem = [0 19388.08; 30 23675.43; 60 22634; 90 25767.41; 120 26369.41];%PAPER135
DivKPTem(:,2) = DivKPTem(:,2)/max(DivKPTem(:,2));
dpDivKP=DivKPTem;

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
% % if strcmp(mutant,mutant)
% LegendSize=9; LabelSize=14;
% 
% 
% 
  figure()%DivJ, PleC, DivK, CtrA
set(gcf,'position',[100 400 1400 450])%left, lower, right, upper

ha = tight_subplot(2,3,[.15 .05],[.1 .1],[.05 .05]);

axes(ha(1));
% PodJL=yout(1:10,:)+yout(11:20,:);   PodJSS=yout(21:30,:);

PodJL_T=sum(PodJLS);  PodJS_T=sum(PodJSS); 
PodJ_T=PodJL_T+PodJS_T;
Fraction_L=PodJL_T./PodJ_T;
Fraction_S=PodJS_T./PodJ_T;
%%data
Time=[0 20 40 60 80 100 120 140 160];
%%Calculate fraction of relative levels by getdata
L=[100 91 19 2 1 12 46 32 53];
S=[0 0 7 50 100 100 41 14 22];
Sum=L+S;

% SpmX=sum(yout(31:40,:)+yout(41:50,:));

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
% h = legend('SpmX_T','Kaczmarczyk\newlineet al. ''20','Location', 'North');
h = legend('PodJL','PodJS','Chen 05','Guo 04','Chen 05','Guo 04','Location', 'North','fontsize',LegendSize,'NumColumns',3);
  set(h, 'Interpreter', 'tex', 'Box', 'on', 'Orientation', 'vertical');
  
  
axes(ha(2));
SpmX = SpmXtot;
  p5 = line(tout,sum(SpmX), 'Color', 'k' , 'LineWidth', 2, 'Linestyle', '--');
hold on
% [score,scalary] = findSquares(tout,sum(SpmX), dpSpmX2(1:end-1,:));
[score,scalary] = findSquares(tout,sum(SpmX), dpSpmX2(1:end-1,:));
p7 = scatter(dpSpmX2(:,1),dpSpmX2(:,2).*scalary, 'k+', 'LineWidth', 1);
hold on
[score,scalary] = findSquares(tout,sum(SpmX), dpSpmX3(1:end-1,:));
p7 = scatter(dpSpmX3(:,1),dpSpmX3(:,2).*scalary, 'k^', 'LineWidth', 1);
hold off
axis([tspan 0 max(sum(SpmX))*1.2]);
% xlabel('time (min)','fontsize',LabelSize)
% ylabel('scaled concentration')
% h = legend('SpmX_T','Kaczmarczyk\newlineet al. ''20','Location', 'North');
h = legend('SpmX','Radhakrishnan 07','Guo 04','Location', 'North','fontsize',LegendSize,'NumColumns',3);
  set(h, 'Interpreter', 'tex', 'Box', 'on', 'Orientation', 'vertical');
  
  
  
axes(ha(3));
% DivJfT=sum(yout(151:160,:)+yout(171:180,:)+yout(191:200,:));
%    DivJbT=sum(yout(161:170,:)+yout(181:190,:)+yout(201:210,:));
%   
% DivJT=DivJfT+DivJbT;
 


p5 = line(tout,sum(DivJtot) , 'Color', 'k' , 'LineWidth', 2, 'Linestyle', '-');
p5 = line(tout,sum(DivJfT) , 'Color', 'r' , 'LineWidth', 2, 'Linestyle', '--');
p5 = line(tout,sum(DivJbT) , 'Color', 'r' , 'LineWidth', 2, 'Linestyle', ':');
hold on
[score,scalary] = findSquares(tout,sum(DivJtot), dpDivJ(1:end-1,:));
p7 = scatter(dpDivJ(:,1),dpDivJ(:,2).*scalary, 'k+', 'LineWidth', 1);
hold on
[score,scalary] = findSquares(tout,sum(DivJtot), dpDivJ2(1:end-1,:));
p7 = scatter(dpDivJ2(:,1),dpDivJ2(:,2).*scalary, 'k^', 'LineWidth', 1);
hold on
[score,scalary] = findSquares(tout,sum(DivJtot), dpDivJ3(1:end-1,:));
p7 = scatter(dpDivJ3(:,1),dpDivJ3(:,2).*scalary, 'k*', 'LineWidth', 1);
hold off
axis([tspan 0 max(sum(DivJtot))*1.8]);
% ylabel('scaled concentration')
h = legend('DivJ_T','DivJ_f','DivJ_b','Wheeler 1999','Sanselicio 2015','Radhakrishnan 07','Location', 'North','fontsize',LegendSize,'NumColumns',3);
  set(h, 'Interpreter', 'tex', 'Box', 'on', 'Orientation', 'vertical');


  
  axes(ha(4));
% PleCph=sum(yout(91:100,:)+yout(101:110,:)+yout(111:120,:)+yout(121:130,:));
% PleCkin=sum(yout(131:140,:)+yout(141:150,:));

PleCT=PleCphT+PleCkinT;

p5 = line(tout,sum(PleCT) , 'Color', 'k' , 'LineWidth', 2, 'Linestyle', '-');
p5 = line(tout,sum(PleCphT) , 'Color', 'r' , 'LineWidth', 2, 'Linestyle', '--');
p5 = line(tout,sum(PleCkinT) , 'Color', 'r' , 'LineWidth', 2, 'Linestyle', ':');
hold on
[score,scalary] = findSquares(tout,sum(PleCT), dpPleCtot(1:end-1,:));
p7 = scatter(dpPleCtot(:,1),dpPleCtot(:,2).*scalary, 'k+', 'LineWidth', 1);
hold on
[score,scalary] = findSquares(tout,sum(PleCT), dpPleCtot2(1:end-1,:));
p7 = scatter(dpPleCtot2(:,1),dpPleCtot2(:,2).*scalary, 'k^', 'LineWidth', 1);
hold off
axis([tspan 0 max(sum(PleCT))*1.6]);
% xlabel('time (min)','fontsize',LabelSize)
ylabel('scaled concentration')
h = legend('PleC_T','PleC_{ph}','PleC_{kin}','Viollier 2002','Guo 04','Location', 'North','fontsize',LegendSize,'NumColumns',3);
 set(h, 'Interpreter', 'tex', 'Box', 'on', 'Orientation', 'vertical');
 
axes(ha(5));
% DivKS=sum(yout(221:230,:)+yout(171:180,:)+yout(181:190,:));
% DivKPS=sum(yout(231:240,:)+yout(191:200,:)+yout(201:210,:)+yout(111:120,:)+yout(121:130,:)...
%     +yout(261:270,:)+yout(271:280,:));

DivKT=DivKS+DivKPS;
p5 = line(tout,sum(DivKT) , 'Color', 'k' , 'LineWidth', 2, 'Linestyle', '-');
p5 = line(tout,sum(DivKS), 'Color', 'b' , 'LineWidth', 2, 'Linestyle', '--');
p5 = line(tout,sum(DivKPS) , 'Color', 'r' , 'LineWidth', 2, 'Linestyle', ':');
hold on
[score,scalary] = findSquares(tout,sum(DivKT), dpDivKTOT2(1:end-1,:));
p7 = scatter(dpDivKTOT2(:,1),dpDivKTOT2(:,2).*scalary, 'k^', 'LineWidth', 1);
hold on
[score,scalary] = findSquares(tout,sum(DivKT), dpDivKTOT3(1:end-1,:));
p7 = scatter(dpDivKTOT3(:,1),dpDivKTOT3(:,2).*scalary, 'k+', 'LineWidth', 1);
hold on
[score,scalary] = findSquares(tout,sum(DivKPS), dpDivKP(1:end-1,:));
p7 = scatter(dpDivKP(:,1),dpDivKP(:,2).*scalary, 'r^', 'LineWidth', 1);
hold off
axis([tspan 0 max(sum(DivKT))*1.6]);
h = legend('DivK_T','DivK','DivK~P','Gao 04','Jacobs 00','Jacobs 00','Location', 'North','fontsize',LegendSize,'NumColumns',3);
 set(h, 'Interpreter', 'tex', 'Box', 'on', 'Orientation', 'vertical');
 
axes(ha(6));
% CtrAS(n+1,:)=0;
% CtrAT=CtrAS+CtrAPS;

p5 = line(tout,sum(CtrAtot) , 'Color', 'k' , 'LineWidth', 2, 'Linestyle', '-');
p5 = line(tout,sum(CtrAPT) , 'Color', 'r' , 'LineWidth', 2, 'Linestyle', '-');
hold on
[score,scalary] = findSquares(tout,sum(CtrAtot), dpCtrAT(1:end-1,:));
p7 = scatter(dpCtrAT(:,1),dpCtrAT(:,2).*scalary, 'k+', 'LineWidth', 1);
hold on
[score,scalary] = findSquares(tout,sum(CtrAtot), dpCtrAT3(1:end-1,:));
p7 = scatter(dpCtrAT3(:,1),dpCtrAT3(:,2).*scalary, 'k^', 'LineWidth', 1);
hold on
[score,scalary] = findSquares(tout,sum(CtrAPT), dpCtrAP(1:end-1,:));
p7 = scatter(dpCtrAP(:,1),dpCtrAP(:,2).*scalary, 'r+', 'LineWidth', 1);
hold off
axis([tspan 0 max(sum(CtrAtot))*1.6]);
h = legend('CtrA_T','CtrA~P','Mcgraph','Radhakrishnan 07','Jacobs 03','Location', 'North','fontsize',LegendSize,'NumColumns',3);
set(h, 'Interpreter', 'tex', 'Box', 'on', 'Orientation', 'vertical');
    end

