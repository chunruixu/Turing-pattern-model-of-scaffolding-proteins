

function [value,isterminal,direction] = podJ_event(t,y,TS)
global T_e1  T_term

CtrAP =sum(y(81:90,:)+y(321:330,:)+y(331:340,:))*0.1;

T_Sphase=90;

T_term=T_e1+T_Sphase;
% T_e6=(T_term-T_e1)*0.125+T_e1;%S_DivL
T_e2=(T_term-T_e1)*0.37+T_e1;%S_CtrA changed from 0 to 1; the fork passes ctrA
T_e3=(T_term-T_e1)*0.65+T_e1;%S_pleC changed from 0 to 1; the fork passes pleC
T_e4=(T_term-T_e1)*0.74+T_e1;%S_perP changed from 0 to 1; the fork passes perP
T_e5=(T_term-T_e1)*0.87+T_e1;%S_podJ changed from 0 to 1; the fork passes podJ


value = [sign(t-15);sign(TS-CtrAP);sign(t - T_e2); sign(t - T_e3); sign(t - T_e4);sign(t - T_e5); sign(t-T_term)];
isterminal = [1;1; 1;1; 1;1;1];
direction = [+1;+1; +1;+1; +1;+1;+1];

end