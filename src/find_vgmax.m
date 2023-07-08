%-------------------------------------------------
%------find the group velocity automatically------
%-------------------------------------------------
function vgmax=find_vgmax(wfc,dis0,tbeg,dt,vg,dvg)
vgmax=[]; 
vgwin=[vg-dvg vg+dvg];
vgwin(vgwin<2)=2;
vgwin(vgwin>5)=5;
tg=dis0./vgwin-tbeg;
twin0=[min(tg) max(tg)];  %time window to be analyzed
kgrp0=ceil((dis0/vg-tbeg)/dt); %--sample # of the estimated velocity

if twin0(1)<0
    twin0(1)=0;
end

nt0=round(twin0/dt)+1;
if nt0(2)>length(wfc)
    %msg=['Bad Wave Length']
    %return
    nt0(2)=length(wfc);
end
% twin0=(nt0-1)*dt;           %time window
nt0=nt0(1):nt0(2); 
nt0(nt0<=0|nt0>length(wfc))=[];
env=abs(hilbert(wfc));
env=env(nt0);  %cut the wave!!
[kmax0]=findmaxima(env,5);  %find the 5 largest points

% figure(1)
% subplot(2,1,1)
% plot(wfc(nt0))
% hold on;
% plot(env,'r')
% plot([kmax0],[env(kmax0)],'ro')

if isempty(kmax0)
    return
end
kmax0=kmax0+nt0(1);
[kmax,imax]=min(abs(kmax0-kgrp0)); %find the one closest to "vg"
kmax=kmax0(imax);                   %the time point
[emax,kmax2]=max(env);              %maximum of env

tmax0=(kmax-1)*dt;
tmax=tmax0+tbeg;
vgmax=dis0/tmax;                    %group velocity of picked pick

%
%-----find the peak value of one envelope----------
%
function [kmax]=findmaxima(dat,mmax)
diff0=diff(dat);
ndf=length(diff0);
nmax=0;
for id=2:ndf-1
    if diff0(id)==0 & diff0(id-1)>0 & diff0(id+1)<0
        nmax=nmax+1;
        pmax(nmax)=id;
        %[blw(nmax),bup(nmax)]=errbars(dat,id,0.95);
    elseif diff0(id)>0 & diff0(id+1)<0
        nmax=nmax+1;
        pmax(nmax)=id+1;
        %[blw(nmax),bup(nmax)]=errbars(dat,id+1,0.95);
    end
end 

if nmax==0
    kmax=[];
    return
end

mmax=min(length(pmax),mmax);
for ii=1:mmax
    [fmax,kk]=max(dat(pmax));
    kmax(ii)=pmax(kk);
    pmax(kk)=[];
end