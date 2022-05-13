function [DivKP, DivK,  CtrAP, CtrA] = SimVSExp (yout,time)
%DivKt2mid is the ratio of total and midcell DivK, which is around 3 in WT
%in experiments (paper194)

% DivKpole = yout(111,:)+yout(120,:)+yout(121,:)+yout(130,:)+yout(231,:)+yout(240,:)+yout(261,:)+yout(270,:)...
%     +yout(271,:)+yout(280,:)+yout(191,:)+yout(200,:)+yout(201,:)+yout(210,:);
% 
% DivKpole = trapz(time,DivKpole)/150;
% DivKmid = (yout(90,:)+yout(91,:)+yout(70,:)+yout(71,:)+yout(74,:)+yout(75,:)).*factor2...
%             +(yout(94,:)+yout(95,:)+yout(78,:)+yout(79,:)+yout(82,:)+yout(83,:)+yout(46,:)+yout(47,:)+yout(50,:)+yout(51,:)+yout(106,:)+yout(107,:)+yout(110,:)+yout(111,:)).*factor2;
% DivKmid = trapz(time,DivKmid)/150;
% 
%         % DivKt2mid = sum(DivKpole)/sum(DivKmid)+1;
% % paper126 page164, fig4 and 5 DivKP/DivK AND CtrAP/CtrA in mutants
N=time(end);
DivK=sum(yout(221:230,:)+yout(171:180,:)+yout(181:190,:));

DivK = trapz(time,DivK)/N;

DivKP=sum(yout(111:120,:)+yout(121:130,:)+yout(231:240,:)+yout(261:270,:)+yout(271:280,:)+yout(191:200,:)+yout(201:210,:));
DivKP = trapz(time,DivKP)/N;
    
CtrAP =sum(yout(81:90,:)+yout(321:330,:)+yout(331:340,:));
CtrAP = trapz(time,CtrAP)/N;

CtrA=sum(yout(71:80,:)+yout(341:350,:)+yout(351:360,:));

CtrA = trapz(time,CtrA)/N;






