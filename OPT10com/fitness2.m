function error=fitness2(para)
warning('error','MATLAB:ode15s:IntegrationTolNotMet');
global p;
ver=9;
%y0



if ver==2||ver==7
load('y0_10com_1.mat')%
elseif ver==8
     load('y0_10com_3.mat')
elseif ver==9
     load('y0_10com_4.mat')
else
    load('y0_10com_2.mat')
end




load_para(ver,'WT');%load fixed parameters


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% PodJS2=[0 7656.953;20 11798.569;40 9196.033;60 7763.104;80 7449.276;100 8469.74; 120 11575.711; 140 11777.761; 160 10382.853];
% PodJS2(:,2) = PodJS2(:,2)/max(PodJS2(:,2));
% dpPodJS2=PodJS2;
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
CtrAT3 = [0 16183.619; 0 16183.619; 0 16183.619; 20 8066.79; 40 1543.861; 60	865.669; 60	865.669;60	865.669; 80		5211.74; 100 11540.64; 120	12304.376; 140	14119.569; 160	23971.874];

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
%%



error(1)=0;
error(2)=0;


YOUT = {};
TOUT = {};
IniT = {};
 CycleNum=3;
for i=1:CycleNum
    try  
        [Y, t, ~,IniTime]=main(y0,para,ver,'WT');%todo: rep ini time
    catch
%         SimPenalty {i}= 1000;
        break
%         if CycleNum==2
%             [DivKP, DivK, DivKpole, DivKmid, CtrAP, CtrA] = SimVSExp (Y,t);
%         end
    end
YOUT{i}=Y;
TOUT{i} = t;
IniT{i} = IniTime;
y0=IniValue(Y,'SW');
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%MOP

    if length(YOUT)< CycleNum
        error(1) = 100;
        error(2) = 100;
    else
        for j=2:CycleNum%sum of errors of sim cycles except for the 1st cycle
            Y=YOUT{j};
            t=TOUT{j};
            inittime = IniT{j};


            
            PodJL=sum(Y(1:10,:)+Y(11:20,:));   PodJS=sum(Y(21:30,:));

            scalePodJ = 20;
            PodJLc = scalePodJ*findSquares(t,PodJL,dpPodJL);%dpPodJL2
            PodJSc = scalePodJ*findSquares(t,PodJS,dpPodJS);%

             CtrAP = sum(Y(81:90,:)+Y(321:330,:)+Y(331:340,:)); 
            CtrA =sum(Y(71:80,:)+Y(341:350,:)+Y(351:360,:));
            CtrAT=CtrA+CtrAP;
            scaleCtrA = 200*2;
            CtrATc =  scaleCtrA*findSquares(t,CtrAT,dpCtrAT);
            CtrAPc = scaleCtrA/3*findSquares(t,CtrAP,dpCtrAP);
            
             
%             concentrationCost = SpmXc + PodJLc + PodJSc + DivJc + CtrATc + CtrAPc + CpdRc+ CckAkinc;
            concentrationCost = PodJLc*4 + PodJSc  + CtrATc ;%+ CckAkinc;
            %%%%%%%%%%%%%%%%%%%%%%
%             SpatialCtrA
% Id = find(t>120); 
% CtrA_sw = sum(Y(81:85,:)+Y(321:325,:)+Y(331:335,:));
% CtrA_st = sum(Y(86:90,:)+Y(326:330,:)+Y(336:340,:));
% DIFF = sum(CtrA_sw(Id(1):Id(end)))-2*sum(CtrA_st(Id(1):Id(end)));
% SpatialCtrAc = mse(min(0,DIFF));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Diff = max(5,abs(inittime-25))-5;
            EventPenalty = mse(max(0,Diff));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%spatial characteristics
%             PodJL_np = (Y(1,:)+Y(5,:));%new pole
            % SpatialPodJc = mse(min(0,PodJL_np*4-PodJL));
            % PopZ_bp = (Y(21,:)+Y(24,:)+Y(25,:)+Y(28,:));%bipolar
            PopZ_center = sum(Y(52:59,:)+Y(62:69,:));
            PopZ_center = sum(PopZ_center)/8;
            PopZ_np = sum(Y(51,:)+Y(61,:));
            PopZ_op = sum(Y(60,:)+Y(70,:));
            PopZtotal = sum(PopZ_center + PopZ_np + PopZ_op)/3;
            SpatialPopZc = mse(min(0,(PopZ_np -4* PopZ_center)/PopZtotal));

            SpatialPopZc2 = mse(min(0,(PopZ_op - 1.5*PopZ_np)/PopZtotal));
            DivKP_np = sum(Y(231,:)+Y(111,:)+Y(121,:)+Y(191,:)+Y(201,:)+Y(261,:)+Y(271,:));
            DivKP_center = sum(Y(232:239,:)+Y(112:119,:)+Y(122:129,:)+Y(192:199,:)+Y(202:209,:)+Y(272:279,:));
            DivKP_center = sum(DivKP_center)/8;
            DivKP_op = sum(Y(240,:)+Y(120,:)+Y(130,:)+Y(200,:)+Y(210,:)+Y(280,:));
            DivKPtotal = sum(DivKP_np + DivKP_op+ DivKP_center)/3;
            SpatialDivKc = mse(min(0,sum(DivKP_np - 2*DivKP_center)/DivKPtotal));
            SpatialCost = SpatialDivKc + SpatialPopZc2*100 + SpatialPopZc*100;
            
             error(1) = error(1)+ concentrationCost+EventPenalty;%+SpatialCtrAc/20;% SpatialCost;%

%             error(2) = error(2)+SpatialPopZc2*10^2+SpatialDivKc + SpatialPopZc*10^2;
            error(2)=error(2)+ SpatialCost;
        end
        
    end
       



end









