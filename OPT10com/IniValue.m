function y0=IniValue(yout,celltype)
Y0=yout(:,end);
n=10;%#of compartment
if strcmp(celltype,'SW')
for i=1:39
y0_(i*n-n+1)=Y0(i*n-n+5);
y0_(i*n-n+2)=Y0(i*n-n+5);
y0_(i*n-n+3)=Y0(i*n-n+4);
y0_(i*n-n+4)=Y0(i*n-n+4);
y0_(i*n-n+5)=Y0(i*n-n+3);
y0_(i*n-n+6)=Y0(i*n-n+3);
y0_(i*n-n+7)=Y0(i*n-n+2);
y0_(i*n-n+8)=Y0(i*n-n+2);
y0_(i*n-n+9)=Y0(i*n-n+1);
y0_(i*n-n+10)=Y0(i*n-n+1);
end
y0_(396)=0.02*10; 
elseif strcmp(celltype,'ST')
    for i=1:39
y0_(i*n-n+1)=Y0(i*n-n+6);
y0_(i*n-n+2)=Y0(i*n-n+6);
y0_(i*n-n+3)=Y0(i*n-n+7);
y0_(i*n-n+4)=Y0(i*n-n+7);
y0_(i*n-n+5)=Y0(i*n-n+8);
y0_(i*n-n+6)=Y0(i*n-n+8);
y0_(i*n-n+7)=Y0(i*n-n+9);
y0_(i*n-n+8)=Y0(i*n-n+9);
y0_(i*n-n+9)=Y0(i*n-n+10);
y0_(i*n-n+10)=Y0(i*n-n+10);
    end
y0_(396)=0.024*10; 
end
y0_(391:395)=0;
y0=y0_';