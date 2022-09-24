function [DivKP, DivK,  CtrAP, CtrA] = SimVSExp (yout,time)

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






