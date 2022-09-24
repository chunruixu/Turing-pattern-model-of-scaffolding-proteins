function result_PlotMutant3(yout,tout,celltype,mutant,TITLE,concentration)
% close all
%spatial plots
%%simulation
%% Scaffolding spatial and temporal
FONTSIZE=14;
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
%%%%%%%%%%%%%Fig 8 WT background
if strcmp(mutant,'WT')
figure()
    set(gcf,'position',[100 100 400 500])%left, lower, right, upper
subplot(2,1,1)
 DivKPT(11,:) = a;
pcolor(time, M, DivKPT)
shading flat
colorbar
caxis([0 0.9]);
title(['(a) WT' ],'fontweight','bold','fontsize',FONTSIZE)
ylabel('DivKP','fontweight','bold','fontsize',FONTSIZE)
subplot(2,1,2)
CtrAP(11,:)=a;
pcolor(time, M, CtrAP)
shading flat
colorbar
% caxis([0 0.1]);
ylabel('CtrAP','fontweight','bold','fontsize',FONTSIZE)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% ha = tight_subplot(2,4,[.06 .03],[.1 .04],[.02 .02]);
if strcmp(mutant,'deltaPodJ')
    figure()
    set(gcf,'position',[100 100 400 500])%left, lower, right, upper
subplot(2,1,1)
 DivKPT(11,:) = a;
pcolor(time, M, DivKPT)
shading flat
colorbar
% caxis([0 0.35]);
title('\Delta {\itpodJ}','fontweight','bold','fontsize',FONTSIZE)
ylabel('DivKP','fontweight','bold','fontsize',FONTSIZE)
subplot(2,1,2)
CtrAP(11,:)=a;
pcolor(time, M, CtrAP)
shading flat
colorbar
% caxis([0 0.1]);
ylabel('CtrAP','fontweight','bold','fontsize',FONTSIZE)



  
    figure()
 set(gcf,'position',[100 100 1800 300])%left, lower, right, upper
ha = tight_subplot(1,4,[.1 .03],[.1 .1],[.05 .04]);

% subplot(1,4,1)
axes(ha(1)); 

    PopZ(11,:)=a;
pcolor(time, M, PopZ)
shading flat
colorbar
% caxis([0 4]);
% xlabel('time (min)')
title('PopZ')
% subplot(1,4,2)
axes(ha(2)); 

SpmX(11,:)=a;
pcolor(time, M, SpmX)
shading flat
colorbar
% caxis([0 0.25]);
% xlabel('time (min)')
title('SpmX')
% subplot(1,4,3)
axes(ha(3)); 

DivJT(11,:)=a;
pcolor(time, M, DivJT)
shading flat
colorbar
% caxis([0.1 0.35]);
% xlabel('time (min)')
title('DivJ')

% subplot(1,4,4)
axes(ha(4)); 

PleC(11,:)=a;
pcolor(time, M, PleC)
shading flat
colorbar
% caxis([0 0.1]);
xlabel('time (min)')
title('PleC')
%%%%%%%%%%%%%%

elseif strcmp(mutant,'deltaSpmX')
    figure()
    set(gcf,'position',[100 100 400 500])%left, lower, right, upper
subplot(2,1,1)
 DivKPT(11,:) = a;
pcolor(time, M, DivKPT)
shading flat
colorbar
% caxis([0 0.35]);
title('\Delta {\itspmX}','fontweight','bold','fontsize',FONTSIZE)
ylabel('DivKP','fontweight','bold','fontsize',FONTSIZE)
subplot(2,1,2)
CtrAP(11,:)=a;
pcolor(time, M, CtrAP)
shading flat
colorbar
% caxis([0 0.1]);
ylabel('CtrAP','fontweight','bold','fontsize',FONTSIZE)
    


 figure()
 set(gcf,'position',[100 100 1800 300])%left, lower, right, upper
ha = tight_subplot(1,4,[.1 .03],[.1 .1],[.05 .04]);
        axes(ha(1)); 
        PodJL(11,:)=a;
pcolor(time, M, PodJL)
shading flat
colorbar
caxis([0 4]);
% xlabel('time (min)')
title('PodJL')
axes(ha(2)); 
PodJS(11,:)=a;
pcolor(time, M, PodJS)
shading flat
colorbar
% caxis([0 0.5]);
% xlabel('time (min)')
title('PodJS')

axes(ha(3)); 
PopZ(11,:)=a;
pcolor(time, M, PopZ)
shading flat
colorbar
% caxis([0 4]);
% xlabel('time (min)')
title('PopZ')

axes(ha(4)); 
DivJT(11,:)=a;
pcolor(time, M, DivJT)
shading flat
colorbar
% caxis([0.1/2.5 0.35/2.5]);
% xlabel('time (min)')
title('DivJ')



elseif strcmp(mutant,'deltaPopZ')
    figure()
    set(gcf,'position',[100 100 400 500])%left, lower, right, upper
subplot(2,1,1)
 DivKPT(11,:) = a;
pcolor(time, M, DivKPT)
shading flat
colorbar
% caxis([0 0.35]);
title('\Delta {\itpopZ}','fontweight','bold','fontsize',FONTSIZE)
ylabel('DivKP','fontweight','bold','fontsize',FONTSIZE)
subplot(2,1,2)
CtrAP(11,:)=a;
pcolor(time, M, CtrAP)
shading flat
colorbar
% caxis([0 0.1]);
ylabel('CtrAP','fontweight','bold','fontsize',FONTSIZE)




    figure()
 set(gcf,'position',[100 100 1800 300])%left, lower, right, upper
ha = tight_subplot(1,4,[.1 .03],[.1 .1],[.05 .04]);
        axes(ha(1)); 
%  subplot(2,1,1)
        PodJL(11,:)=a;
pcolor(time, M, PodJL)
shading flat
colorbar
% caxis([0.1 0.23]);
% xlabel('time (min)')
title('PodJL')
        axes(ha(2)); 
        SpmX(11,:)=a;
pcolor(time, M, SpmX)
shading flat
colorbar
% caxis([0 0.2]);
% xlabel('time (min)')
title('SpmX')
        axes(ha(3)); 
DivJT(11,:)=a;
pcolor(time, M, DivJT)
shading flat
colorbar
caxis([0.08 0.18]);
% xlabel('time (min)')
title('DivJ')
        axes(ha(4)); 
        CpdRu(11,:)=a;
pcolor(time, M, CpdRu)
shading flat
colorbar
caxis([1.5 3.5]);
% xlabel('time (min)')
title('CpdRu')



  

elseif strcmp(mutant,'PleC-F778L')%%%%%%%%%%%%%%%%%%%%
    figure()
 set(gcf,'position',[100 100 1800 300])%left, lower, right, upper
ha = tight_subplot(1,4,[.1 .03],[.1 .1],[.05 .04]);

% subplot(1,4,1)
axes(ha(1)); 
DivKtot(11,:) = a;
        pcolor(time, M, DivKtot)
shading flat
colorbar
% caxis([0 0.9]);
% xlabel('time (min)')
title('DivKtot')

axes(ha(2)); 
DivLDivKPT(11,:)=a;
pcolor(time, M, DivLDivKPT)
shading flat
colorbar
% caxis([0 0.25]);
% xlabel('time (min)')
title('DivLDivKPT')
% subplot(1,4,3)
axes(ha(3)); 
CckAkinT(11,:)=a;
pcolor(time, M, CckAkinT)
shading flat
colorbar
% caxis([0.1 0.35]);
% xlabel('time (min)')
title('CckA kinase')


axes(ha(4)); 
CtrAP(11,:)=a;
pcolor(time, M, CtrAP)
shading flat
colorbar
% caxis([0 0.1]);
% xlabel('time (min)')
title('CtrAP')
    
  figure()
    set(gcf,'position',[100 100 400 500])%left, lower, right, upper
subplot(2,1,1)
 DivKPT(11,:) = a;
pcolor(time, M, DivKPT)
shading flat
colorbar
caxis([0 0.9]);
title(['(f) {\itpleC-F778L}' ],'fontweight','bold','fontsize',FONTSIZE)
ylabel('DivKP','fontweight','bold','fontsize',FONTSIZE)
subplot(2,1,2)
CtrAP(11,:)=a;
pcolor(time, M, CtrAP)
shading flat
colorbar
% caxis([0 0.1]);
ylabel('CtrAP','fontweight','bold','fontsize',FONTSIZE)

elseif strcmp(mutant,'p3:deletingPodJPleCbinding')%%%%%%%%%%%%%%%%%%%%
    figure()
 set(gcf,'position',[100 100 1800 300])%left, lower, right, upper
ha = tight_subplot(1,4,[.1 .03],[.1 .1],[.05 .04]);

% subplot(1,4,1)
axes(ha(1)); 
DivKtot(11,:) = a;
        pcolor(time, M, DivKtot)
shading flat
colorbar
% caxis([0.1 0.35]);
% xlabel('time (min)')
title('DivKtot')

axes(ha(2)); 
DivLDivKPT(11,:)=a;
pcolor(time, M, DivLDivKPT)
shading flat
colorbar
% caxis([0 0.25]);
% xlabel('time (min)')
title('DivLDivKPT')
% subplot(1,4,3)
axes(ha(3)); 
CckAkinT(11,:)=a;
pcolor(time, M, CckAkinT)
shading flat
colorbar
% caxis([0.1 0.35]);
% xlabel('time (min)')
title('CckA kinase')


axes(ha(4)); 
CtrAP(11,:)=a;
pcolor(time, M, CtrAP)
shading flat
colorbar
% caxis([0 0.1]);
% xlabel('time (min)')
title('CtrAP')
    
  figure()
    set(gcf,'position',[100 100 400 500])%left, lower, right, upper
subplot(2,1,1)
 DivKPT(11,:) = a;
pcolor(time, M, DivKPT)
shading flat
colorbar
caxis([0 0.9]);
title('(g) delocalized PleC','fontweight','bold','fontsize',FONTSIZE)
ylabel('DivKP','fontweight','bold','fontsize',FONTSIZE)
subplot(2,1,2)
CtrAP(11,:)=a;
pcolor(time, M, CtrAP)
shading flat
colorbar
% caxis([0 0.1]);
ylabel('CtrAP','fontweight','bold','fontsize',FONTSIZE)
elseif strcmp(mutant,'p4:deletingPodJDivLbinding')%%%%%%%%%%%%%%%%%%%%
    figure()
 set(gcf,'position',[100 100 1800 300])%left, lower, right, upper
ha = tight_subplot(1,4,[.1 .03],[.1 .1],[.05 .04]);

% subplot(1,4,1)
axes(ha(1)); 
% DivKPT(5,:) = a;
% pcolor(time, M, DivKPT)
% shading flat
% colorbar
% caxis([0 0.35]);
% % xlabel('time (min)')
% title('DivKPT')
DivKtot(11,:) = a;
        pcolor(time, M, DivKtot)
shading flat
colorbar
% caxis([0.1 0.35]);
% xlabel('time (min)')
title('DivKtot')

axes(ha(2)); 
DivLDivKPT(11,:)=a;
pcolor(time, M, DivLDivKPT)
shading flat
colorbar
% caxis([0 0.25]);
% xlabel('time (min)')
title('DivLDivKPT')
% subplot(1,4,3)
axes(ha(3)); 
CckAkinT(11,:)=a;
pcolor(time, M, CckAkinT)
shading flat
colorbar
% caxis([0.1 0.35]);
% xlabel('time (min)')
title('CckA kinase')


axes(ha(4)); 
CtrAP(11,:)=a;
pcolor(time, M, CtrAP)
shading flat
colorbar
% caxis([0 0.1]);
% xlabel('time (min)')
title('CtrAP')
    
  figure()
    set(gcf,'position',[100 100 400 500])%left, lower, right, upper
subplot(2,1,1)
 DivKPT(11,:) = a;
pcolor(time, M, DivKPT)
shading flat
colorbar
caxis([0 0.9]);
title('(h) delocalized DivL','fontweight','bold','fontsize',FONTSIZE)
ylabel('DivKP','fontweight','bold','fontsize',FONTSIZE)
subplot(2,1,2)
CtrAP(11,:)=a;
pcolor(time, M, CtrAP)
shading flat
colorbar
% caxis([0 0.1]);
ylabel('CtrAP','fontweight','bold','fontsize',FONTSIZE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mutant,'deltaDivJ')
     figure()
    set(gcf,'position',[100 100 400 500])%left, lower, right, upper
subplot(2,1,1)
DivKtot(11,:) = a;
        pcolor(time, M, DivKtot)
shading flat
colorbar
ylabel('DivKtot','fontweight','bold','fontsize',FONTSIZE)
title('(b) \Delta {\itdivJ}','fontweight','bold','fontsize',FONTSIZE)
subplot(2,1,2)
CtrAP(11,:)=a;
pcolor(time, M, CtrAP)
shading flat
colorbar
% caxis([0 0.1]);
ylabel('CtrAP','fontweight','bold','fontsize',FONTSIZE)
    
         figure()
 set(gcf,'position',[100 100 1800 300])%left, lower, right, upper
ha = tight_subplot(1,4,[.1 .03],[.1 .1],[.05 .04]);

% subplot(1,4,1)
axes(ha(1)); 
DivKtot(11,:) = a;
        pcolor(time, M, DivKtot)
shading flat
colorbar
% caxis([0.1 0.35]);
% xlabel('time (min)')
title('DivKtot')
axes(ha(2)); 
DivKPT(11,:) = a;
pcolor(time, M, DivKPT)
shading flat
colorbar
% caxis([0 4]);
% xlabel('time (min)')
title('DivKPT')

axes(ha(3)); 
DivLDivKPT(11,:)=a;
pcolor(time, M, DivLDivKPT)
shading flat
colorbar
% caxis([0 0.25]);
% xlabel('time (min)')
title('DivLDivKPT')
% subplot(1,4,3)



axes(ha(4)); 
DivLtot(11,:)=a;
pcolor(time, M, DivLtot)
shading flat
colorbar
% caxis([0 0.1]);
% xlabel('time (min)')
title('DivLtot')

figure()
DivLtot(11,:)=a;
pcolor(time, M, DivLtot)
shading flat
colorbar
% caxis([0 0.1]);
ylabel('DivL','fontweight','bold','fontsize',FONTSIZE)
title('\Delta divJ','fontweight','bold','fontsize',FONTSIZE)
elseif strcmp(mutant,'DivJ-H338A')
      figure()
    set(gcf,'position',[100 100 400 500])%left, lower, right, upper
subplot(2,1,1)
%  DivKPT(5,:) = a;
% pcolor(time, M, DivKPT)
% shading flat
% colorbar
% caxis([0 0.35]);
DivKtot(11,:) = a;
        pcolor(time, M, DivKtot)
shading flat
colorbar
% caxis([8 11]);
ylabel('DivKtot','fontweight','bold','fontsize',FONTSIZE)
title(['(c) {\itdivJ-H338A}'],'fontweight','bold','fontsize',FONTSIZE)
subplot(2,1,2)
CtrAP(11,:)=a;
pcolor(time, M, CtrAP)
shading flat
colorbar
% caxis([6 8]);
ylabel('CtrAP','fontweight','bold','fontsize',FONTSIZE)
    
         figure()
 set(gcf,'position',[100 100 1800 300])%left, lower, right, upper
ha = tight_subplot(1,4,[.1 .03],[.1 .1],[.05 .04]);

% subplot(1,4,1)
axes(ha(1)); 
DivKtot(11,:) = a;
        pcolor(time, M, DivKtot)
shading flat
colorbar
title('DivKtot')
axes(ha(2)); 
DivKPT(11,:) = a;
pcolor(time, M, DivKPT)
shading flat
colorbar
% caxis([0 4]);
% xlabel('time (min)')
title('DivKPT')

axes(ha(3)); 
DivLDivKPT(11,:)=a;
pcolor(time, M, DivLDivKPT)
shading flat
colorbar
% caxis([0 0.25]);
% xlabel('time (min)')
title('DivLDivKPT')
% subplot(1,4,3)



axes(ha(4)); 
DivLtot(11,:)=a;
pcolor(time, M, DivLtot)
shading flat
colorbar
% caxis([0 0.1]);
% xlabel('time (min)')
title('DivLtot')

     

elseif strcmp(mutant,'deltaPleC')
    
    figure()
    set(gcf,'position',[100 100 400 500])%left, lower, right, upper
subplot(2,1,1)
 DivKPT(11,:) = a;
pcolor(time, M, DivKPT)
shading flat
colorbar
% caxis([0 1.2]);
title('(d) \Delta {\itpleC}','fontweight','bold','fontsize',FONTSIZE)
ylabel('DivKP','fontweight','bold','fontsize',FONTSIZE)
subplot(2,1,2)
CtrAP(11,:)=a;
pcolor(time, M, CtrAP)
shading flat
colorbar
% caxis([0 0.1]);
ylabel('CtrAP','fontweight','bold','fontsize',FONTSIZE)

     figure()
 set(gcf,'position',[100 100 1800 300])%left, lower, right, upper
ha = tight_subplot(1,4,[.1 .03],[.1 .1],[.05 .04]);


axes(ha(1)); 
DivKtot(11,:) = a;
        pcolor(time, M, DivKtot)
shading flat
colorbar
% caxis([0.1 0.35]);
% xlabel('time (min)')
title('DivKtot')

axes(ha(2)); 
DivLDivKPT(11,:)=a;
pcolor(time, M, DivLDivKPT)
shading flat
colorbar
% caxis([0 0.25]);
% xlabel('time (min)')
title('DivLDivKPT')
% subplot(1,4,3)
axes(ha(3)); 
CckAkinT(11,:)=a;
pcolor(time, M, CckAkinT)
shading flat
colorbar
% caxis([0.1 0.35]);
% xlabel('time (min)')
title('CckA kinase')
axes(ha(4)); 
DivJT(11,:)=a;
pcolor(time, M, DivJT)
shading flat
colorbar
caxis([0 1]);
% xlabel('time (min)')
title('DivJ')
elseif strcmp(mutant,'PleC-H610A')%%%%%%%%%%%%%%%%%%%%
        figure()
    set(gcf,'position',[100 100 400 500])%left, lower, right, upper
subplot(2,1,1)
 DivKPT(11,:) = a;
pcolor(time, M, DivKPT)
shading flat
colorbar
% caxis([0 0.35]);
title(['(e) {\itpleC-H610A}'],'fontweight','bold','fontsize',FONTSIZE)
ylabel('DivKP','fontweight','bold','fontsize',FONTSIZE)
subplot(2,1,2)
CtrAP(11,:)=a;
pcolor(time, M, CtrAP)
shading flat
colorbar
% caxis([0 0.1]);
ylabel('CtrAP','fontweight','bold','fontsize',FONTSIZE)
 figure()
 set(gcf,'position',[100 100 1800 300])%left, lower, right, upper
ha = tight_subplot(1,4,[.1 .03],[.1 .1],[.05 .04]);

% subplot(1,4,1)
axes(ha(1)); 
DivKPT(11,:) = a;
pcolor(time, M, DivKPT)
shading flat
colorbar
% caxis([0 4]);
% xlabel('time (min)')
title('DivKPT')

axes(ha(2)); 
PleC(11,:)=a;
pcolor(time, M, PleC)
shading flat
colorbar
% caxis([0 0.25]);
% xlabel('time (min)')
title('PleC')
axes(ha(3)); 
CckAkinT(11,:)=a;
pcolor(time, M, CckAkinT)
shading flat
colorbar
% caxis([0.1 0.35]);
% xlabel('time (min)')
title('CckA kinase')
axes(ha(4)); 
CtrAP(11,:)=a;
pcolor(time, M, CtrAP)
shading flat
colorbar
% caxis([0 0.1]);
% xlabel('time (min)')
title('CtrAP')
end

