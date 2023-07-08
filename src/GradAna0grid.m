%% 
%% --------------Wave Gradiometry for a gridpoint-----------------
%%
function [ab,azm0,stc,vgmax]=GradAna0grid(st,ev,par,fwin,kwvmod)         
kredu=par.kredu;
stc=st(1);
dt=stc.dt; %#ok<*NASGU>
azm0=azimuth0(stc.st(1),stc.st(2),stc.ev(1),stc.ev(2))+180; 
if azm0>360
    azm0=azm0-360;
end
if azm0<0
    azm0=azm0+360;
end
nsti=length(st);
vgmax=par.vgmax;
if strcmp(kredu,'vrphs')
    dv=100;
    nloop=0;
    while dv>par.dvmax
        nloop=nloop+1;
        [ab,vphs,stc,uprdt]=GradAna8redu(st,ev,par,fwin,kwvmod);
        if isempty(stc) || isempty(ab) || isnan(vphs) || isempty(vphs) || nloop>25
            ab.pr=[];
            azm0=[];
            stc=[];
            vgmax=[];
            return
        end
        if nloop==1
            vgmax=find_vgmax(stc.dat,stc.dis,stc.tbeg,stc.dt,vgmax,100);
            if isempty(vgmax)
                vgmax=par.vgmax;
            end
        end
        dv=abs(vphs-par.vredu);
        par.vredu=vphs;
    end
end
if isempty(uprdt)
    ab.pr=[];
    azm0=[];
    stc=[];
    vgmax=[];
end


    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WG: par.vredu >0 <0, w/o reducing velocity
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ab,vphs,stc,uprdt]=GradAna8redu(st,ev,par,fwin,kwvmod)         
stc=st(1);
dt=stc.dt;
    vgmax=par.vgmax;
    if par.vredu>0
        [st,sc1]=timeshift(st,ev,par.vredu,dt);
        if isempty(st)||length(st(2).dat)<2*par.wlen
            ab = [];
            vphs = [];
            stc = [];
            uprdt=[];
            return
        end
    end
    [dudx,dudy,stc,uprdt]=WGAgrid_dxdy(st,ev,fwin,par,sc1);
    if isempty(dudx)
            ab = [];
            vphs = [];
            stc = [];
        return
    end
    stc.dst=distance0(stc.st(1),stc.st(2),stc.ev(1),stc.ev(2))*111.1949;
%     clear st;
    [ab,wfc,dt]=WGA_AB(stc,dudx,dudy,par.kazm,kwvmod); %#ok<*ASGLU>
    azm0=azimuth0(stc.st(1),stc.st(2),stc.ev(1),stc.ev(2))+180;
    if azm0>360
        azm0=azm0-360;
    end
    [~,kmax]=find_azm(ab,stc,dt,vgmax);    
    if isempty(kmax)
            ab = [];
            vphs = [];
            stc = [];
        return
    end     
    ab2=WGA_PhyParam(par.vredu,ab,azm0,stc.dst,par.kazm,kmax);
    [azmij,kmax]=find_azm(ab2,stc,dt,vgmax);
    dazm=azmij-azm0;

    if dazm>180
        azmij=azmij-360;
        dazm=azmij-azm0;
    elseif dazm<-180
        azmij=azmij+360;
        dazm=azmij-azm0;
    end

    if abs(dazm)>45
        ab=ab2;
    else
        ab=WGA_PhyParam(par.vredu,ab,azmij,stc.dst,par.kazm,kmax);
    end
    vphs=find_vphase(ab,stc,dt,vgmax,par.wlen);

%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
return
figure(898) %#ok<*UNRCH>
subplot(1,2,1)
stloc=cell2mat({st.st}');
scatter(stloc(2:end,2),stloc(2:end,1),'^','fill')
hold on
scatter(stloc(1,2),stloc(1,1),'or','fill')
plot([stloc(1,2) st(2).ev(2)],[stloc(1,1) st(2).ev(1)],'b')
box on
xlim([min(stloc(:,2)) max(stloc(:,2))])
ylim([min(stloc(:,1)) max(stloc(:,1))])
hold off
subplot(1,2,2)
for j=2:length(st)
    plot((1:length(wfc))-kmax,st(j).dat,'k')
    hold on
end
plot((1:length(wfc))-kmax,stc.dat,'color',[0 153 255]/255,'Linewidth',2)
hold off
title('uobs')
xlim([-200 200])


%%-----------------------------------------------
%% add back the effects of the reducing velocity\
%%-----------------------------------------------
function ab=WGA_PhyParam(vredu,ab,azm0,dst,kazm,kmax)
    if vredu>0
        predu=1/vredu;
        ab.Bx=ab.Bx-predu*sin(azm0/180*pi);
        ab.By=ab.By-predu*cos(azm0/180*pi);
    end
            %-----Azimuth Estimate---------------
        if kazm==0
            ab.theta=atan2(-ab.Bx,-ab.By);
        else
            ab.theta=atan2(ab.dudx,ab.dudy);
        end

            %-----Ray Parameter Estimate-----------
        ab.pr=sqrt(ab.Bx.*ab.Bx+ab.By.*ab.By);
%         am=azm0/180*pi;
        am=ab.theta(kmax);        
        ab.radiation=(ab.Ax.*cos(am)-ab.Ay.*sin(am))*dst;
        ab.gsprd=(ab.Ax.*sin(am)+ab.Ay.*cos(am));
    

%%
%%      -----reduce time-----
%%
function [st,sc1]=timeshift(st,ev,vredu,dt)
nst=length(st);
%shift the waveforms.
nts(1:nst)=nan;
ddis(1:nst-1)=nan;
for is=1:nst
    dis0=distance0(st(is).st(1),st(is).st(2),ev(1),ev(2))*111.1949;
    tsh=dis0/vredu;
    nts(is)=ceil(tsh/dt);
    if is>1
        ddis(is-1) = distance0(st(1).st(1),st(1).st(2),st(is).st(1),st(is).st(2))*111.1949;
    end
end
sc1 = find(ddis==min(ddis))+1;
sc1 = sc1(1);
nts=nts-min(nts);

%make sure all waveforms are in same length
[nmax , ~]=max(nts);
np=length(st(sc1).dat)-nmax;
if np<1
    st = [];
    return
end

for is=2:nst
    st(is).tsh=nts(is)*dt;
    st(is).tbeg=st(is).tbeg+st(is).tsh;
    if nts(is)>0
        st(is).dat(1:nts(is))=[];
    end
    npi=length(st(is).dat);
    if npi>np
        st(is).dat(np+1:npi)=[];
    end
    st(is).dat=sinWin(st(is).dat,0.1,1);
end



%% -----------------------------------
%% ------   -------------
%% -----------------------------------
function vphs=find_vphase(ab,stc,dt,vg,twin)
kmax=ceil((stc.dst/vg-stc.tbeg)/dt); %--point number of estimated velocity
nt=round(twin/dt);                  %time window to be averaged.
nt=(-nt:1:nt)+kmax;
nt(nt<1)=[];
nt(nt>length(ab.pr))=[];
pr=ab.pr(nt);                      %ray-parameter in the time window
vapp=1.0./pr;                 %velocity in the average-window
vphs=mean(vapp);




%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
return
%---test----
subplot(2,1,1)
plot(1.0./ab.pr)
hold on;
plot(nt,vapp,'r')
ylim([0 10])

function [azm,kmax]=find_azm(ab,stc,~,~)
% kmax=ceil((stc.dst/vg-stc.tbeg)/dt); %--point number of estimated velocity
env=abs(hilbert(stc.dat));
[kmax]=findmaxima(env,1); 
azm=ab.theta(kmax)*180/pi;                      %ray-parameter in the time window
azm(azm<0)=azm(azm<0)+360;
azm(azm>360)=azm(azm>360)-360;
return
%---test----
subplot(2,1,1)
plot(1.0./ab.pr)
hold on;
plot(nt,vapp,'r')
ylim([0 10])





%% -----------------------------------
%% ------Output results-------------
%% -----------------------------------
function OutPutVelAzm(fido,ab,azm0,stc,nsti,dt,vgpick,twin,ev) %#ok<*DEFNU>
stn=stc.stn;
st=stc.st;
[envu,~]=envelope(stc.dat);
tpeak=stc.dst/vgpick-stc.tbeg;
nt=round(tpeak/dt);
nts=nt+(-200:200);
nts(nts<1)=[];
nts(nts>length(envu))=[];
kmax=findmaxima(envu(nts),1);
if isempty(kmax)
    kmax=findmaxima(envu,1);
else
    kmax=kmax+nts(1);
end
dnt=round(twin/dt);
nt1=kmax-dnt;
nt2=kmax+dnt;
vgpick=stc.dst/(kmax*dt+stc.tbeg);


if nt2>length(ab.pr); nt2=length(ab.pr); end
if nt1<1; nt1=1; end
        
if nt2<=nt1
    return
end

kmax2=findmaxima(envu,3);
snr_nsts=min(kmax2):max(kmax2);
snr_nsts(snr_nsts<1)=[];
snr_nsts(snr_nsts>length(ab.pr))=[];
NsL1=1:snr_nsts(1);
NsL2=snr_nsts(end):length(ab.pr);
envu0=envu(kmax);
SNR1=envu0/mean(envu([NsL1 NsL2]));

WL=round(5*twin);
snr_nsts=kmax+(-WL:WL);
snr_nsts(snr_nsts<1)=[];
snr_nsts(snr_nsts>length(ab.pr))=[];
NsL1=1:snr_nsts(1);
NsL2=snr_nsts(end):length(ab.pr);
envu0=envu(kmax);
SNR2=envu0/mean(envu([NsL1 NsL2]));
SNR=max([SNR2,SNR1]);

vel=1./ab.pr(nt1:nt2);
azm=ab.theta(nt1:nt2)*180/pi;
azm(azm>360)=azm(azm>360)-360;
azm(azm<0)=azm(azm<0)+360;
dazm=azm-azm0;
dazm(dazm<-180)=dazm(dazm<-180)+360;
dazm(dazm>360)=dazm(dazm>360)-360;

gsp=ab.gsprd(nt1:nt2)*1000;          %Geometrical Spreading * 1000
rad=ab.radiation(nt1:nt2);
ampl0=stc.dat(kmax);
Time=stc.tbeg+kmax/dt;
azmo=ab.theta(kmax)*180/pi;
azmo(azmo>360)=azmo(azmo>360)-360;
azmo(azmo<0)=azmo(azmo<0)+360;
a0=azmo-azm0;
a0(a0>180)=a0(a0>180)-360;
a0(a0<-180)=a0(a0<-180)+360;

v0=1/ab.pr(kmax);
g0=ab.gsprd(kmax)*1000;          %Geometrical Spreading * 1000
r0=ab.radiation(kmax);

dv=sqrt(var(vel));
da=sqrt(var(dazm));
dg=sqrt(var(gsp));
dr=sqrt(var(rad));

Bx = ab.Bx(kmax);
By = ab.By(kmax);
Ax = ab.Ax(kmax);
Ay = ab.Ay(kmax);

fprintf(fido,'%s %7.3f %7.3f %7.2f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %8.3f %2g %10.6f %10.6f %8.3f %7.3f %7.3f %7.3f %11.8f %11.8f %10.6f %10.6f %7.3f\n', ...
    stn,st(1),st(2),st(3),v0,dv,a0,da,r0,dr,g0,dg,vgpick,azmo,nsti,envu0,ampl0,Time,ev(1),ev(2),ev(3),Ax,Ay,Bx,By,SNR);


%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
return
figure(22)
subplot(5,1,1)
set( gca, 'XTick', []);
plot(stc.dat,'k');
hold on
plot(envu,'color',[0.2,0.7,0.2]);
plot(nt1:nt2,stc.dat(nt1:nt2),'b');
scatter(kmax,envu(kmax),'ro','fill');
xlim([1 length(envu)])
set( gca, 'XTick', []);
hold off

subplot(5,1,2)
set( gca, 'XTick', []);
plot(1.0./ab.pr,'k');
hold on
plot(nt1:nt2,1.0./ab.pr(nt1:nt2),'b');
scatter(kmax,1.0./ab.pr(kmax),'ro','fill');
hold off
ylim(1.0./ab.pr(kmax)+[-0.5 0.5])
xlim([1 length(envu)])
set( gca, 'XTick', []);
title(['V_p_h_a_s_e=' num2str(round(v0*100)/100) '¡À' num2str(round(dv*1000)/1000) ' (km/s)'])
hold off

subplot(5,1,3)
set( gca, 'XTick', []);
plot(ab.theta*180/pi,'k');
hold on
scatter(kmax,ab.theta(kmax)*180/pi,'ro','fill');
plot(nt1:nt2,ab.theta(nt1:nt2)*180/pi,'b');
ylim(ab.theta(kmax)*180/pi+[-45 45])
xlim([1 length(envu)])
set( gca, 'XTick', []);
title(['Azimuth=' num2str(round(azmo*100)/100) '¡À' num2str(round(da*100)/100) ' (deg)'])
hold off

subplot(5,1,4)
set( gca, 'XTick', []);
plot(ab.gsprd*1000,'k');
hold on
scatter(kmax,ab.gsprd(kmax)*1000,'ro','fill');
plot(nt1:nt2,ab.gsprd(nt1:nt2)*1000,'b');
ylim(ab.gsprd(kmax)*1000+40*dg*[-1 1])
xlim([1 length(envu)])
set( gca, 'XTick', []);
title(['Geometrical spreading=' num2str(round(g0*100)/100) '¡À' num2str(round(dg*100)/100) ' (1/1000km)'])
hold off

subplot(5,1,5)
plot((1:length(ab.radiation))*dt,ab.radiation,'k');
hold on
scatter(kmax,ab.radiation(kmax),'ro','fill');
plot(nt1:nt2,ab.radiation(nt1:nt2),'b');
ylim(ab.radiation(kmax)+40*dr*[-1 1])
xlim([1 length(envu)])
xlabel('time (s)')
title(['Radiation pattern=' num2str(round(r0*100)/100) '¡À' num2str(round(dr*100)/100) ' (1/deg)'])
hold off



%% ------------------------------------------------
%% -------calculate Ax Bx Ay and By values---------
%% ------------------------------------------------
function [AB,Dc,dt]=WGA_AB(stc,dudx,dudy,kazm,kmodel)
Dc=stc.dat;
dt=stc.dt;
np=length(Dc);
% nh=np/2;
% df=1.0/dt/np;

[envdx,phadx]=envelope0(dudx);
[envdy,phady]=envelope0(dudy);
[envu,phau]=envelope0(Dc);

wt=InsantFrq(Dc,dt);
envudt=diff(envu)/dt;
envudt(np)=envudt(np-1);

if strfind(kmodel,'pla')
    AB.kmodel=0;
    %Bx=envdx./envu./wt.*sin(phadx-phau);
    %By=envdy./envu./wt.*sin(phady-phau);

    dudt=diff(Dc)/dt;
    dudt(np)=dudt(np-1);
    [envdudt,~]=envelope0(dudt);
    Bx=envdx./envdudt;
    By=envdy./envdudt;

    Ax=0.0;
    Ay=0.0;

    AB.Ax=Ax;
    AB.Ay=Ay;
    AB.Bx=Bx;
    AB.By=By;

    AB.xc=stc.st(2);
    AB.yc=stc.st(1);
        %-----Azimuth Estimate---------------
    if kazm~=0
        theta=atan2(dudx,dudy);
    else
        theta=atan2(-AB.Bx,-AB.By);
    end
    AB.theta=theta;
        %-----Ray Parameter Estimate-----------
    %AB.pr=-AB.Bx.*sin((theta))-AB.By.*cos((theta));
    AB.pr=sqrt(AB.Bx.*AB.Bx+AB.By.*AB.By);
    %AB.radiation=(AB.Ax.*cos(theta)-AB.Ay.*sin(theta))*stc.dst;
    %AB.gsprd=(AB.Ax.*sin(theta)+AB.Ay.*cos(theta));
elseif strfind(kmodel,'sph')
    AB.kmodel=1;
    Bx=envdx./envu./wt.*sin(phadx-phau);
    By=envdy./envu./wt.*sin(phady-phau);

    Ax=envdx./envu.*cos(phadx-phau)-envdx./envu./envu./wt.*(envudt).*sin(phadx-phau);
    Ay=envdy./envu.*cos(phady-phau)-envdy./envu./envu./wt.*(envudt).*sin(phady-phau);

    AB.Ax=Ax;
    AB.Ay=Ay;
    AB.Bx=Bx;
    AB.By=By;
    
    AB.xc=stc.st(2);
    AB.yc=stc.st(1);
        %-----Azimuth Estimate---------------
    if kazm~=0
        theta=atan2(dudx,dudy);
    else
        theta=atan2(-AB.Bx,-AB.By);
    end
    AB.theta=theta;
        %-----Ray Parameter Estimate-----------
    %AB.pr=-AB.Bx.*sin((theta))-AB.By.*cos((theta));
    AB.pr=sqrt(AB.Bx.*AB.Bx+AB.By.*AB.By);
    %AB.radiation=(AB.Ax.*cos(theta)-AB.Ay.*sin(theta))*stc.dst;
    %AB.gsprd=(AB.Ax.*sin(theta)+AB.Ay.*cos(theta));
    AB.envu=envu;
    
    
    
    %% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    return
    subplot(4,1,1)
    plot(Dc)
    grid on
    
    subplot(4,1,2)
    plot(dudx)
    grid on

    subplot(4,1,3)
    plot(dudy)
    grid on

    subplot(4,1,4)
    plot(AB.Bx)
    grid on
    
end

function [env,pha]=envelope0(dat0)
    hilb=hilbert(dat0);
    env=abs(hilb);
    pha=atan2(imag(hilb),real(hilb));
    return
    %---------------copy from Chuck's "surface.m" program---------
    hilb=hilbert(dat0);
    ihilb=imag(hilb);
    env=abs(dat0 - i*ihilb);
    pha=atan2(-ihilb,dat0);
    return


%%
%% calculate the instantenous frequency--
%%
function wt=InsantFrq(u,dt)
np=length(u);
ut=diff(u)/dt;
ut(np)=ut(np-1);
uh=-imag(hilbert(u));
uth=-imag(hilbert(ut));
uu=abs(hilbert(u));
uu=uu.*uu;
wt=(ut.*uh-u.*uth)./uu;




%%
%% -----find the peak value of one envelope----------
%%
function [kmax]=findmaxima(dat,mmax)
dat=sinWin(dat,0.1,1);
diff0=diff(dat);
ndf=length(diff0);
nmax=0;
pmax(1:length(2:ndf-1))=nan;
for id=2:ndf-1
    if diff0(id)==0 && diff0(id-1)>0 && diff0(id+1)<0
        nmax=nmax+1;
        pmax(nmax)=id;
        %[blw(nmax),bup(nmax)]=errbars(dat,id,0.95);
    elseif diff0(id)>0 && diff0(id+1)<0
        nmax=nmax+1;
        pmax(nmax)=id+1;
        %[blw(nmax),bup(nmax)]=errbars(dat,id+1,0.95);
    end
end 
pmax(isnan(pmax))=[];
if nmax==0
    kmax=[];
    return
end

mmax=min(length(pmax),mmax);
kmax(1:mmax)=nan;
for ii=1:mmax
    [fmax,kk]=max(dat(pmax));
    kmax(ii)=pmax(kk);
    pmax(kk)=[];
end
    
