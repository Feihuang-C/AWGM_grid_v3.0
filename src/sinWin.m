 
function dat=sinWin(X,a,b)
nd=length(X);
nh=round(a*length(X)/2)*2;
dat(1:nd,1)=X;
win1=sin(0:pi/nh:pi/2).^b;
win2=sin(pi/2:pi/nh:pi).^b;
bdx=1:length(win1);
edx=(length(dat)-length(win2)+1):length(dat);
dat(bdx)=dat(bdx).*win1';
dat(edx)=dat(edx).*win2';