function [PodJm PodJp PodJS SpmXm SpmXp PopZm PopZp CtrA CtrAP PleCf PleCb...
    PleCfDivKP PleCbDivKP PleCfkin PleCbkin DivJf DivJb DivJfDivK DivJbDivK...
    DivJfDivKP DivJbDivKP PerP DivK DivKP DivLf DivLb DivLfDivKP DivLbDivKP...
    CckAfkin CckAbkin CckAfph CckAbph CtrAPCckAfph CtrAPCckAbph CtrACckAfkin CtrACckAbkin...
    CpdRf CpdRb CpdRP N PoleNum n]=INPUT()
%%%n should be multiple of ten
n=10;
species=39;%# of species (eg. proteins, mRNA etc)
N=n*species+5+1;%# of equations. 5 is #of S; 1 is # of compartment type
PoleNum=floor(n/10);

PodJm=1; 
PodJp=2; 
PodJS=3; 
SpmXm=4; 
SpmXp=5; 
PopZm=6; 
PopZp=7;
CtrA=8; 
CtrAP=9; 
PleCf=10; 
PleCb=11; 
PleCfDivKP=12; 
PleCbDivKP=13; 
PleCfkin=14;
PleCbkin=15;
DivJf=16; 
DivJb=17;
DivJfDivK=18; 
DivJbDivK=19; 
DivJfDivKP=20; 
DivJbDivKP=21;
PerP=22;
DivK=23; 
DivKP=24;  
DivLf=25; 
DivLb=26; 
DivLfDivKP=27; 
DivLbDivKP=28;
CckAfkin=29; 
CckAbkin=30; 
CckAfph=31; 
CckAbph=32;
CtrAPCckAfph=33;
CtrAPCckAbph=34;
CtrACckAfkin=35;
CtrACckAbkin=36;
CpdRf=37; 
CpdRb=38; 
CpdRP=39;