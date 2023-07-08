 
function dat=HmWin(X,a,b)
nd=length(X);
nh=round(a*length(X)/2)*2; 
dat(1:nd,1)=X;
Hwin(1:nh,1)=hamming(nh);
Hwin=Hwin.^b;
% Hwin1=Hwin(1:length(Hwin)/2);
% Hwin2=Hwin(length(Hwin)/2+1:length(Hwin));

Hwin1=sin(0:pi/nh:pi/2);
Hwin2=sin(pi/2:pi/nh:pi);

bdx=1:length(Hwin1);
edx=(length(dat)-length(Hwin2)+1):length(dat);
dat(bdx)=dat(bdx).*Hwin1';
dat(edx)=dat(edx).*Hwin2';