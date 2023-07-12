%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Liang Chuntao 2009, v1.0; Cao Feihuang. 2022, v2.0; Cao Feihuang, 2023/07/01,v2.0
% befor runing this program，you need to check the parameters file of aa_WGM_parameters.m firstly.
% before running the program, also check all parameters below
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ac_step2_WGM_grid
clear
close all;
%% setup parameters 
par=aa_WGM_parameters;% load parameters
dataF = dir([par.wvPick_output par.fs '20*.mat']); %#ok<*NODEF>
is_parfor = 0;  % =1,multicore computing; =0,single core computing
is_overwrite= 1;
is_figure = 1;
is_smooth=1; 






%% ————————————————————————————————————
par.is_smooth =is_smooth; 
par.is_overwrite =is_overwrite; 
par.is_figure=is_figure;
par.is_parfor=is_parfor;
par.kphase = 'R'; %%P/R/S/L
par.tsac=[];   %%%T0-9 .or. [];
par.omark='.nz00pct.';
par.knz=0.0; %noise level: 0.0, no noise added to the real data
par.kredu='vrphs'; %% 'vr???/vrphs/vrgrp/vrnan';
par.kwi=1;  %%if kwi<=0, no weight applied
par.kazm=0;  %%%always be 0, other values for test purpose
par.kgauss=[]; %%don't change this parameter in this file
par.kwvmodel='sph';
par.dataF = dataF;
para_WG='para_WG.mat';
save(para_WG)
disp('WGM analysing...')
if is_parfor==1 
    parfor kevi=1:length(dataF) %loops of events
        parforWGcl(kevi,para_WG)
    end
else
    for kevi=1:length(dataF) %loops of events
        parforWGcl(kevi,para_WG)
    end
end
disp('done')
delete(para_WG)






%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   one earthquake, loops in periods
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function parforWGcl(kevi,para_WG)
% try
load(para_WG); %#ok<*LOAD> 
dataF=par.dataF;
fs=par.fs;
is_overwrite=par.is_overwrite;
is_figure=par.is_figure;
WGoutdir = par.WG_output;
periods = par.periods;%center periods
cmp=par.cmp;
kphase=par.kphase;
par.dsmax=par.maxstadist;
par.tsac=[];
omark=par.omark;
kredu=par.kredu;
kwi=par.kwi;
kwvmodel=par.kwvmodel;
temp = load([par.wvPick_output,dataF(kevi).name]);
rsta = temp.rsta;
pickfile=[par.wvPick_output,fs 'wvpick' fs 'wvpick.' dataF(kevi).name];

if ~exist(pickfile,'file')
    disp('no picking info')
    return
end

temp = load(pickfile);
wvpick=temp.wvpick;
if wvpick(1).evflag==0
    return
end
%%
evtstr=strsplit(dataF(kevi).name,'.mat');
par.eloc=cell2mat(evtstr(1));
par.dt=round(median([rsta.dt])*1000)/1000;
dout=[WGoutdir fs par.eloc fs];
gsta=getimeInfo(par.gsta,rsta,par);
pickPRD=cell2mat({wvpick.period});
disp([num2str(kevi) 'th event: ' par.eloc]);
for pbi=1:length(periods) %WGM: loops of periods
    CT=periods(pbi);
    pband = [CT*0.9 CT*1.1];
    fwin=[1/pband(2) 1/pband(1)];
    par.wlen=max(pband)/2;
    pstr=[kphase '.p' num2str2(CT,3,0)];
    par.pband=pband;    
    if kwi==-1
        strwi='.wtn.';
    elseif kwi==0
        strwi='.wt0.';
    elseif kwi==1
        strwi='.wt1.';
    elseif kwi==2
        strwi='.wt2.';
    end
    fout=[dout 'wga.' pstr '.' cmp '.' par.eloc strwi kredu omark 'mat'];
    fout2=[dout 'wga.smr.' pstr '.' cmp '.' par.eloc strwi kredu omark 'mat'];
    if is_overwrite==0 && exist(fout2,'file')
        continue
    end
        picidx=find(pickPRD==CT);
        
    if isempty(picidx)
        [datFlg,vgs,goodstID,kmax0,SNRs] = TSNR2(rsta,fwin,par,is_figure);% find the good waveforms
        goodstID=find(goodstID==1);
        vgs=vgs(goodstID);
        kmax0=kmax0(goodstID);
        SNRs=SNRs(goodstID);
    else
        wvpickP=wvpick(picidx);
        datFlg=wvpickP.datFlag;
        goodstID=find(wvpickP.goodstID==1);
        vgs=wvpickP.vg(goodstID);
        kmax0=wvpickP.kmax0(goodstID);
        SNRs=wvpickP.SNR(goodstID);
    end
    par.vg=median(vgs);
    par.vredu = par.vg+0.2; %%%starting reducing velocity
    mSNR=median(SNRs);
    if datFlg==0||par.vg>par.vglen(2)||par.vg<par.vglen(1)||mSNR<par.mSNR
        continue
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    if is_overwrite==1 
        wga=WGManalysis(kwvmodel,kwi,par,gsta,rsta,goodstID,kmax0,vgs);  %  wave gradiometry anslysis for one reference location
        if wga.datFlg==1 
            if ~exist(dout,'dir')
                mkdir(dout)
            end
            save(fout,'wga')  %  save results of one period
        end             
    else
        if exist(fout,'file') && is_smooth==1% load old results
            load(fout);
        else
            wga=WGManalysis(kwvmodel,kwi,par,gsta,rsta,goodstID,kmax0,vgs);  %   wave gradiometry anslysis for one reference location
            if wga.datFlg==1
                if ~exist(dout,'dir')
                    mkdir(dout)
                end
                save(fout,'wga')  %  save results of one period
            end
        end
    end 
    
    if wga.datFlg==1 && is_smooth==1
        wga1=smoothWG(wga,par,CT,par.eloc); % smooth result
        if ~isempty(wga1)
            if ~exist(dout,'dir')
                mkdir(dout)
            end
            wga=wga1;
            save(fout2,'wga') % save smoothing result
        end
    end
end
% catch
% disp([num2str(kevi) 'th event']);
% end







%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   WGManalysis
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function wga=WGManalysis(kwvmodel,kwi,par,gsta,rsta,goodstID,kmaxs,vgs)
cutWin=par.cutWin; % s
ev = rsta(1).ev;
par.ev=ev;
kdat=par.eloc;
pband=par.pband;
fwin=[1/pband(2) 1/pband(1)];
par.kwi=kwi;
par.ermax=[];
par.vavg=par.vredu;
par.kdat=kdat;
%--find all stations in the given region--
sta = rsta(goodstID);
nst = length(gsta);
k=0;
for ist=1:nst   %station
    [subarr,idx,dsmax]=Findsta(gsta,sta,ist,par);
    kmaxi=kmaxs(idx);
    if length(subarr)<5
        continue
    end
    % % cut waveforms
    subarr=cutwaveforms(subarr,kmaxi,cutWin);
    % %
    par.dt=subarr(2).dt;
    %---match times---
    subarr = timmatch(subarr);
    vgmax = mean(vgs(idx));
    par.vgmax=vgmax;    
    %narrow band filter is applied in stFilter
    subarr=stFilter(subarr,par.tsac,vgmax,pband,200);
    % %
%             figure(2)
%             dt=subarr(2).dt;
%             for kst=2:length(subarr)
%                 wf1=subarr(kst).dat;
%                 plot([(1:length(wf1))*dt+subarr(kst).tbeg],subarr(kst).dat/max(subarr(kst).dat)+2*kst,'r')
%                 hold on
%             end
    % %
    
    if isempty(subarr)||length(subarr)<5
        continue
    end
    if par.knz>0
        subarr=addnoise(subarr,par.knz);
    end
    nst=length(idx);
    par.vgmax=vgmax;
    [ab,azm0,stc,vgmax]=GradAna0grid(subarr,ev,par,fwin,kwvmodel);
    if isempty(ab) || length(ab.pr)<par.wlen
        continue
    end
    k=k+1;
    wgai=Pickwga(ab,azm0,stc,vgmax,par.wlen,ev,nst); 
    wga.stn(k,:)=wgai.stn;
    wga.st(k,:)=wgai.st;
    wga.ev(k,:)=wgai.ev;
    wga.v(k,:)=wgai.v;
    wga.dv(k,:)=wgai.dv;
    wga.a0(k,:)=wgai.a0;
    wga.da(k,:)=wgai.da;
    wga.rd(k,:)=wgai.rd;
    wga.dr(k,:)=wgai.dr;
    wga.gs(k,:)=wgai.gs;
    wga.dg(k,:)=wgai.dg;
    wga.azmo(k,:)=wgai.azmo;
    wga.nst(k,:)=wgai.nst;
    wga.envu0(k,:)=wgai.envu0;
    wga.amp(k,:)=wgai.amp;
    wga.Time(k,:)=wgai.Time;
    wga.Bx (k,:)=wgai.Bx;
    wga.By (k,:)=wgai.By;
    wga.Ax (k,:)=wgai.Ax;
    wga.Ay (k,:)=wgai.Ay;
    wga.SNR(k,:)=wgai.SNR;
    wga.dsmax(k,:)=dsmax;     %#ok<*STRNU>
end
if k>0
    wga.datFlg=1;    
else
    wga.datFlg=0;
end


function wga=Pickwga(ab,azm0,stc,vgpick,twin,ev,nst)
stn=stc.stn;
st=stc.st;
dt=stc.dt;
envu=ab.envu;
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
amp=stc.dat(kmax);
Time=stc.tbeg+kmax/dt;
azmo=ab.theta(kmax)*180/pi;
azmo(azmo>360)=azmo(azmo>360)-360;
azmo(azmo<0)=azmo(azmo<0)+360;
a0=azmo-azm0;
a0(a0>180)=a0(a0>180)-360;
a0(a0<-180)=a0(a0<-180)+360;

v=1/ab.pr(kmax);
gs=ab.gsprd(kmax)*1000;          %Geometrical Spreading * 1000
rd=ab.radiation(kmax);

dv=sqrt(var(vel));
da=sqrt(var(dazm));
dg=sqrt(var(gsp));
dr=sqrt(var(rad));
wga.stn=stn;
wga.st=st;
wga.ev=ev;
wga.v=v;
wga.dv=dv;
wga.a0=a0;
wga.da=da;
wga.rd=rd;
wga.dr=dr;
wga.gs=gs;
wga.dg=dg;
wga.azmo=azmo;
wga.nst=nst;
wga.envu0=envu0;
wga.amp=amp;
wga.Time=Time;
wga.Bx = ab.Bx(kmax);
wga.By = ab.By(kmax);
wga.Ax = ab.Ax(kmax);
wga.Ay = ab.Ay(kmax);
wga.SNR=SNR;

%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
return
figure(22) %#ok<*UNRCH>
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
title(['V_p_h_a_s_e=' num2str(round(v*100)/100) '±' num2str(round(dv*1000)/1000) ' (km/s)'])
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
title(['Azimuth=' num2str(round(azmo*100)/100) '±' num2str(round(da*100)/100) ' (deg)'])
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
title(['Geometrical spreading=' num2str(round(gs*100)/100) '±' num2str(round(dg*100)/100) ' (1/1000km)'])
hold off

subplot(5,1,5)
plot((1:length(ab.radiation))*dt,ab.radiation,'k');
hold on
scatter(kmax,ab.radiation(kmax),'ro','fill');
plot(nt1:nt2,ab.radiation(nt1:nt2),'b');
ylim(ab.radiation(kmax)+40*dr*[-1 1])
xlim([1 length(envu)])
xlabel('time (s)')
title(['Radiation pattern=' num2str(round(rd*100)/100) '±' num2str(round(dr*100)/100) ' (1/deg)'])
hold off


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              cut waveforms
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  stn=cutwaveforms(arr,kmaxs,nT)
stn=arr;
nst=length(stn);
dt=stn(2).dt;
maxkmax=max(kmaxs);
minkmax=min(kmaxs);
timeL=round(nT/dt);
cutTwin=round(minkmax-timeL):round(maxkmax+timeL);
cutTwin(cutTwin<1)=[];
cutTwin(cutTwin>length(stn(2).dat))=[];
%%
% figure(1)
% hold off
%%
for ist=2:nst
    %%
%         wf1=arr(ist).dat;
%         plot([(1:length(wf1))*dt+arr(ist).tbeg],wf1/max(wf1)+2*ist,'k')
%         hold on
%         wf2=arr(ist).dat(cutTwin);
%         plot(cutTwin+arr(ist).tbeg*dt,wf2/max(wf2)+2*ist,'r')
    %%
    stn(ist).dat=stn(ist).dat(cutTwin);
    stn(ist).tbeg=stn(ist).tbeg+cutTwin(1)*dt-dt;
    stn(ist).t0(5)=stn(ist).t0(5)+cutTwin(1)*dt-dt;
    stn(ist).t1(5)=stn(ist).t0(5)+length(stn(ist).dat)*dt;
    stn(ist).npt=length(cutTwin);
    %
%         subplot(1,2,2)
%         wf1=stn(ist).dat;
%         plot([(1:length(wf1))*dt+stn(ist).tbeg],stn(ist).dat/max(stn(ist).dat)+2*ist,'k')
%         hold on
    %%
end



%% -------------------------------------------
%%      add noise                           %
%% -------------------------------------------
function st=addnoise(st,knz)
if isempty(knz)
    return
end

nst=length(st);
for ii=1:nst
    npt=length(st(ii).dat);
    wmax=max(abs(st(ii).dat));
    wfnz=rand(1,npt)*2-1.0;
    wfnz=wfnz*(knz*wmax);
    st(ii).dat=st(ii).dat+wfnz;
end


%
% station selections
%
function subarray=stFilter(stn,tsac,vgmax,pband,twin)
st=stn(2:end);
stc=stn(1);
nst=length(st);
amax(1:nst)=nan;
for ii=1:nst
    dt=st(ii).dt;
    tg=st(ii).dis/vgmax-st(ii).tbeg;
    tw=[tg-twin tg+twin];
    nt=ceil(tw/dt);
    nt=nt(1):nt(2);
    nt(nt<1|nt>length(st(ii).dat))=[];
    if ~isempty(nt)
        %[wf]=gaussfilt(st(ii).dt,0.5,st(ii).dis,mean(pband),st(ii).dat,0);
        wf=bpfilt(st(ii).dat,dt,1/pband(2),1/pband(1));
        st(ii).dat=wf;
        amax(ii)=max(abs(wf(nt)));
    else
        amax(ii)=amax(1)*100;
    end
    
    %%%If no value picked: discard
    if ~isempty(tsac)
        ktim=str2double(tsac(2))+1;
        tpk=st(ii).tphs(ktim);
        if isnan(tpk)
            amax(ii)=amax(1)*100;
        end
    end
end
ist=1:nst;
amax=amax/median(amax);
%amplitude variation should be less than 30%
ist(amax>1.3)=[];
amax=amax(ist);
ist(amax<0.7)=[];
st=st(ist);
subarray(1)=stc;
subarray(2:length(st)+1)=st;

%-------------------------------------------------
%---match the original time and sample rate-------
%-------------------------------------------------
function [subarray]=timmatch(stn)
st=stn(2:end);
stc=stn(1);
% dtm=median([st.dt]);
nst=length(st);
tm2=cell2mat({st.t0}');
t0=tm2(:,2)*24*3600+tm2(:,3)*3600+tm2(:,4)*60+tm2(:,5);
dt=[st.dt]';
npt=[st.npt]';
t1=t0+dt.*npt-dt(1);
dt=mean(round(dt*100000)/100000);

%--match all time series with same the same starting and ending time----
t0max=max(t0);
t1min=min(t1);
tlen=(t1min-t0max);
if tlen<0
    subarray=[];
    return
end

nt0=round((t0max-t0)./dt)+1;
nt1=nt0+round(tlen./dt)-1;

for is=1:nst
    n0=nt0(is);
    n1=nt1(is);
    st(is).dat=st(is).dat(n0:n1);
    st(is).dat=sinWin(st(is).dat,0.1,3);
    st(is).t0(5)=st(is).t0(5)+t0max-t0(is);
    st(is).tbeg=st(is).tbeg+t0max-t0(is);
end
stc.tbeg= st(1).tbeg;
subarray(1)=stc;
subarray(2:nst+1)=st;



%% make  grid begTime
function gsta=getimeInfo(gsta,rsta,par)
stloc=cell2mat({rsta(:).st}');
gstloc=cell2mat({gsta(:).st}');
st_latRang=[min(stloc(:,1))-par.gridsize,max(stloc(:,1))+par.gridsize];
st_lonRang=[min(stloc(:,2))-par.gridsize,max(stloc(:,2))+par.gridsize];
bkID= gstloc(:,1)>=st_latRang(1)&gstloc(:,1)<=st_latRang(2)&...
    gstloc(:,2)>=st_lonRang(1)&gstloc(:,2)<=st_lonRang(2);
gsta=gsta(bkID);


function [stn,idx,dsmax]=Findsta(gsta,rsta,ist,par)
stc=gsta(ist);
dsmax=par.dsmax;
% plot(stc.st(2),stc.st(1),'bo')
% hold on
% gst=cell2mat({gsta(:).st}');
% plot(gst(:,2),gst(:,1),'.k')
% rst=cell2mat({rsta(:).st}');
% plot(rst(:,2),rst(:,1),'^k')
% axis equal
% drawnow

stloc=cell2mat({rsta(:).st}');
yi=stc.st(1);
xi=stc.st(2);
dst=distance(yi,xi,stloc(:,1),stloc(:,2))*111.1949;
idx=find(dst<=dsmax);
while  length(idx)<par.minNST && dsmax<par.maxdsmax
    dsmax=dsmax+10;
    idx=find(dst<=dsmax);
end

if length(idx)<3
    stn=[];
    idx=[];
    %     hold off
    return
end

stID=(find(abs(stloc(:,1)-stc.st(1))<(3*dsmax/111.1949)&abs(stloc(:,2)-stc.st(2))<(3*dsmax/111.1949)));
if isempty(stID)||length(stID)<5
    stn=[];
    idx=[];
    %     hold off
    return
end
stn(2:length(idx)+1)=rsta(idx);
stn(1).ev = rsta(1).ev;
stn(1).st = stc.st;
stn(1).stn = stc.stn;
stn(1).dt = round(rsta(1).dt*1000)/1000;
stn(1).scal = rsta(1).scale;
stn(1).dis = distance(stc.st(1),stc.st(2),rsta(1).ev(1),rsta(1).ev(2))*111.1949;
% stnloc=cell2mat({stn.st}');
% plot(stnloc(:,2),stnloc(:,1),'r^')
% drawnow
% hold off




function  [datFlg,vg,gd1,kmax0,SNRs]= TSNR2(sts,fwin,par,is_figure)
par.vgmin=par.vglen(1);
par.vgmax=par.vglen(2);
kmax0(1:length(sts))=nan;
vg(1:length(sts))=0;
SNRs(1:length(sts))=0;
gd1(1:length(sts))=0;
for i=1:length(sts)
    wf = bpfilt(sts(i).dat,sts(i).dt,fwin(1),fwin(2));
    ei = abs(hilbert(wf));
    wv.ei(i,:)=ei;
    wv.dist(i) = sts(i).dis;
    wv.begtime(i) = sts(i).tbeg;
    wv.pbwfi(i,:)=wf;
    %time windows (par.vbg-par.ved)
    mxt=round((sts(i).dis/par.vgmin-sts(i).tbeg)/sts(i).dt);
    mit=round((sts(i).dis/par.vgmax-sts(i).tbeg)/sts(i).dt);
    int=(mit:mxt);
    int(int>length(sts(i).dat))=[];
    int(int<1)=[];
    if isnan(max(wv.ei(i,int)))
        continue
    end
    kmax0(i) = (find(wv.ei(i,:)==max(wv.ei(i,int))));
    vg(i) = wv.dist(i)/(wv.begtime(i)+kmax0(i)*sts(i).dt);
end
vg(isnan(vg))=[];
vg(abs(vg-median(vg))>par.dvg|vg>=par.vgmax|vg<=par.vgmin)=[];
if isempty(vg)
    datFlg=0;
    return
end

maxVg=min([par.vgmax,median(vg)+1]);
minVg=max([par.vgmin,median(vg)-1]);
mvg=mean(vg);

for i=1:length(sts)
    %time windows (par.vbg-par.ved)
    mxt=round((sts(i).dis/minVg-sts(i).tbeg)/sts(i).dt);
    mit=round((sts(i).dis/maxVg-sts(i).tbeg)/sts(i).dt);
    int=(mit:mxt);
    int(int>length(sts(i).dat))=[];
    int(int<1)=[];
    if isnan(max(wv.ei(i,int)))
        continue
    end
    kmaxs=findmaxima(wv.ei(i,int),5)+mit*sts(i).dt;
    kvgs = wv.dist(i)./(wv.begtime(i)+kmaxs*sts(i).dt);
    dvg=abs(kvgs-median(vg));
    kmxID=dvg==min(dvg);
    kmax0(i)=kmaxs(kmxID);
    vg(i) = kvgs(kmxID);
    %time windows with 20 wavelength
    int2=round(mean(10./fwin)/sts(i).dt);
    int2=(-int2:int2)+kmax0(i);
    int2(int2<1)=[];
    int2(int2>length(wv.ei(i,:)))=[];
    SNRs(i)= max(wv.ei(i,int2))/mean([wv.ei(i,int2(end)+1:length(wv.ei(i,:))),wv.ei(i,1:int2(1))]);
end

mvg= median(vg);
wv.dt=sts(1).dt;
wv. kmax0= kmax0;
wv.mvg = mvg;
wv.SNRs=SNRs;
wv.vg=vg;
wv.goodsta=zeros(1,length(vg));
evstdist=[sts.dis];
gd1 = SNRs>=par.SNR & abs(vg-wv.mvg)<=par.dvg & vg>par.vgmin & vg<par.vgmax & evstdist>=par.ev_minDist & evstdist<=par.ev_maxDist;
wv.goodsta(gd1)=1;
wv.mvg = median(vg(wv.goodsta==1));
wv.ktmax=max(kmax0);
wv.ktmin=min(kmax0);
wv.kmax0=wv.kmax0;
if is_figure
    plotQwave(wv,par);
end
gId=wv.goodsta==1;
wv.dist=wv.dist(gId);
wv.begtime=wv.begtime(gId);
wv.kmax0=wv.kmax0(gId);
wv.ei=wv.ei(gId,:);
wv.pbwfi=wv.pbwfi(gId,:);
wv.SNRs=wv.SNRs(gId);
wv.vg=wv.vg(gId);
wv.goodsta=wv.goodsta(gId);

if length(find(gId==1))<10
    datFlg=0;
else
    raxu=polyfit([sts.dis],kmax0,1);
    if raxu(1)<0
        datFlg=0;
    else
        datFlg=1;
    end
end


function [kmax]=findmaxima(dat,mmax)
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
    [~,kk]=max(dat(pmax));
    kmax(ii)=pmax(kk);
    pmax(kk)=[];
end


function plotQwave(wvQ,par)
fwin=[1/par.pband(2) 1/par.pband(1)];
par.vgmin=par.vglen(1);
par.vgmax=par.vglen(2);
ist=1:length(wvQ.ei(:,1));
mixi(ist)=nan;
maxi(ist)=nan;
ddyi(ist)=nan;
for i=1:length(wvQ.ei(:,1))
    wf2 = wvQ.pbwfi(i,:);
    mxt=round((wvQ.dist(i)/par.vgmin-wvQ.begtime(i))/wvQ.dt);
    mit=round((wvQ.dist(i)/par.vgmax-wvQ.begtime(i))/wvQ.dt);
    int=mit:mxt;
    int(int>length(wf2))=[];
    int(int<1)=[];
    if isempty(int)||isnan(max(wvQ.ei(i,:)))
        continue
    end
    figure(11);
    hold on
    X=(1:length(wf2))*wvQ.dt+wvQ.begtime(i);
    Y=20*wf2/max(wf2)+wvQ.dist(i)-1;
    plot(X,Y,'color',[0.2 0.2 0.2])
    %     X=int*wvQ.dt+wvQ.begtime(i);
    %     Y=40*wf2(int)/max(wf2(int))+wvQ.dist(i)-1;
    %     plot(X,Y,'color',[90 90 90]/255)
    mixi(i)=min(X);
    maxi(i)=max(X);
    ddyi(i)=10+wvQ.dist(i)-1;
end

% plot(min(X),10+wvQ.dist(i)-1,'.r')
% plot(max(X),10+wvQ.dist(i)-1,'.r')
% plot([min(mixi),max(mixi)],[min(ddyi)-100,max(ddyi)+100],'g')
% plot([min(maxi),max(maxi)],[min(ddyi)-100,max(ddyi)+100],'g')
% plot([max(mixi),max(maxi)],[max(ddyi)+100,max(ddyi)+100],'g')
% plot([min(mixi),min(maxi)],[min(ddyi)-100,min(ddyi)-100],'g')

plot(wvQ.kmax0(wvQ.goodsta==1)*wvQ.dt+wvQ.begtime(wvQ.goodsta==1),ddyi(wvQ.goodsta==1),'.g')
plot(wvQ.kmax0(wvQ.goodsta==0)*wvQ.dt+wvQ.begtime(wvQ.goodsta==0),ddyi(wvQ.goodsta==0),'.r')
xlabel('time (s)')
ylabel('epicentral distance (km)')
xlim([min(mixi),max(maxi)])
ylim([min(ddyi)-100,max(ddyi)+100])
box on
hold off

figure(12)
dist=wvQ.dist;
mvg = wvQ.mvg;
hold on
h=plot(dist(wvQ.goodsta==1),wvQ.vg(wvQ.goodsta==1),'o');
set(h,'MarkerFaceColor',[0.4 0.4 0.7],'MarkerEdgeColor',[0.5 0.5 0.5])
scatter(dist(wvQ.goodsta==0),wvQ.vg(wvQ.goodsta==0),'ok')
plot([min(dist) max(dist)],[mvg+par.dvg,mvg+par.dvg],'b')
plot([min(dist) max(dist)],[mvg-par.dvg,mvg-par.dvg],'b')
xlabel('epicentral distance (km)')
ylabel('vg(km/s)')
box on
hold off

if length(find(wvQ.goodsta==1))<1
    %     disp('!!!no good wavefroms')
    figure(15);
    hold on
    for i=1:length(wvQ.ei(:,1))
        figure(15)
        wf = wvQ.pbwfi(i,:);
        mxt=round((wvQ.dist(i)/par.vgmin-wvQ.begtime(i))/wvQ.dt);
        mit=round((wvQ.dist(i)/par.vgmax-wvQ.begtime(i))/wvQ.dt);
        int=mit:mxt;
        int(int>length(wf))=[];
        int(int<1)=[];
        if isempty(int)||isnan(max(wvQ.ei(i,:)))
            continue
        end
        X=(1:length(wf))*wvQ.dt+wvQ.begtime(i);
        Y=40*wf/max(wf)+wvQ.dist(i)-1;
        plot(X,Y,'color',[0.2 0.2 0.2])
        %         X=int*wvQ.dt+wvQ.begtime(i);
        %         Y=40*wf(int)/max(wf(int))+wvQ.dist(i)-1;
        %         plot(X,Y,'color',[90 90 90]/255)
        mixi(i)=min(X);
        maxi(i)=max(X);
        ddyi(i)=10+wvQ.dist(i)-1;
    end
    
    plot(wvQ.kmax0(wvQ.goodsta==1)*wvQ.dt+wvQ.begtime(wvQ.goodsta==1),ddyi(wvQ.goodsta==1),'.g')
    plot(wvQ.kmax0(wvQ.goodsta==0)*wvQ.dt+wvQ.begtime(wvQ.goodsta==0),ddyi(wvQ.goodsta==0),'.r')
    
    xlim([min(mixi),max(maxi)])
    xlabel('epicentral distance (km)')
    ylabel('vg(km/s)')
    box on
    hold off
    return
else
    for i=1:length(wvQ.ei(:,1))
        if wvQ.goodsta(i)==1
            figure(13)
            wf = wvQ.pbwfi(i,:);
            mxt=round((wvQ.dist(i)/par.vgmin-wvQ.begtime(i))/wvQ.dt);
            mit=round((wvQ.dist(i)/par.vgmax-wvQ.begtime(i))/wvQ.dt);
            int=mit:mxt;
            int(int>length(wf))=[];
            int(int<1)=[];
            if isempty(int)||isnan(max(wvQ.ei(i,:)))
                continue
            end
            hold on
            X=(1:length(wf))*wvQ.dt+wvQ.begtime(i);
            Y=40*wf/max(wf)+wvQ.dist(i)-1;
            plot(X,Y,'color',[0.2 0.2 0.2])
            X=int*wvQ.dt+wvQ.begtime(i);
            Y=40*wf(int)/max(wf)+wvQ.dist(i)-1;
            plot(X,Y,'color',[150 150 150]/255)
            
            int2=round(mean(10./fwin)/wvQ.dt);
            int2=[-int2,int2]+wvQ.kmax0(i);
            int2(int2<1)=1;
            int2(int2>length(wvQ.ei(i,:)))=length(wvQ.ei(i,:));
            int22(i,:)=int2; %#ok<AGROW>
        else
            figure(14)
            wf = wvQ.pbwfi(i,:);
            mxt=round((wvQ.dist(i)/par.vglen(1)-wvQ.begtime(i))/wvQ.dt);
            mit=round((wvQ.dist(i)/par.vglen(2)-wvQ.begtime(i))/wvQ.dt);
            int=mit:mxt;
            int(int>length(wf))=[];
            int(int<1)=[];
            if isempty(int)||isnan(max(wvQ.ei(i,:)))
                continue
            end
            hold on
            X=(1:length(wf))*wvQ.dt+wvQ.begtime(i);
            Y=40*wf/max(wf)+wvQ.dist(i)-1;
            plot(X,Y,'color',[0.2 0.2 0.2])
            X=int*wvQ.dt+wvQ.begtime(i);
            Y=40*wf(int)/max(wf)+wvQ.dist(i)-1;
            plot(X,Y,'color',[150 150 150]/255)
            
            %time windows with 10 wavelength
            int2=round(mean(10./fwin)/wvQ.dt);
            int2=[-int2,int2]+wvQ.kmax0(i);
            int2(int2<1)=1;
            int2(int2>length(wvQ.ei(i,:)))=length(wvQ.ei(i,:));
            if isempty(int2)
                plot(X,Y,'r')
                continue
            end
            int22(i,:)=int2;             %#ok<AGROW>
            text(wvQ.kmax0(i)*wvQ.dt,max(Y),num2str(wvQ.SNRs(i)));
            
        end
    end
    figure(13)
    plot(wvQ.kmax0(wvQ.goodsta==1)*wvQ.dt+wvQ.begtime(wvQ.goodsta==1),ddyi(wvQ.goodsta==1),'.g')
    int22s=[int22(wvQ.goodsta==1,:),ddyi(wvQ.goodsta==1)'];
    
	plot(int22s(:,1)*wvQ.dt+(wvQ.begtime(wvQ.goodsta==1))',int22s(:,3),'.cyan')
	plot(int22s(:,2)*wvQ.dt+(wvQ.begtime(wvQ.goodsta==1))',int22s(:,3),'.cyan')
    title('Good waveforms')
    xlabel('time (s)')
    ylabel('epicentral distance (km)')
    xlim([min(mixi),max(maxi)])
    box on
    hold off
    
    figure(14)
    plot(wvQ.kmax0(wvQ.goodsta==0)*wvQ.dt+wvQ.begtime(wvQ.goodsta==0),ddyi(wvQ.goodsta==0),'.r')
    int22s=[int22(wvQ.goodsta==0,:),ddyi(wvQ.goodsta==0)'];    
    plot(int22s(:,1)*wvQ.dt+(wvQ.begtime(wvQ.goodsta==0))',int22s(:,3),'.cyan')
    plot(int22s(:,2)*wvQ.dt+(wvQ.begtime(wvQ.goodsta==0))',int22s(:,3),'.cyan')
    title('Bad waveforms')
    xlabel('time (s)')
    ylabel('epicentral distance (km)')
    xlim([min(mixi),max(maxi)])
    box on
    hold off
end



%%
function  wga1=smoothWG(wga,par,cntrprd,evID)
    periods=par.periods;
    is_figure=par.is_figure;
    vlim=par.vlim;
    vmod=par.vmod;
    smrvs=par.smrvs;
    smras=par.smras;
    smradv=deg2km(smrvs(periods==cntrprd));
    smradaniso=deg2km(smras(periods==cntrprd));    
    cntrV=vmod(vmod(:,2)==cntrprd,1); 
    [wga1,weight,idx]=datselect(wga,cntrV,vlim);

    %%
    if length(wga1.st)<10||weight<0.3
        wga1=[];
        return
    end
    wga1=VcrAmp(wga1,par.gridsize,cntrprd);
    %%
    wga1=smoothavg(wga1,smradv,smradaniso);  
    if is_figure==1&&par.is_parfor==0
        figure(1)
        plot_map2(wga1,cntrprd,evID)
        drawnow
        hold off    
        figure(2)
        plot_map(wga,idx,cntrprd,evID)
        drawnow
        hold off
    end



%%-------------------------------
%% data selection
%%-------------------------------
function [wga,weight,idx]=datselect(wga,cntrV,vlim)
SNR=wga.SNR;
SNRwi=log10(SNR);
v=wga.v;
Ldx=2;
v1=max([2,cntrV-vlim]);
v2=min([5,cntrV+vlim]);
bednb=v>v2 | v<v1;
bednb2=v>5 | v<2;
idx2=~bednb & ~bednb2 ;
weight=1-(length(find(bednb))/length(v)+2*(length(find(bednb2))/length(v)))/3;

V95=wga.Ax;
V1=median(V95)-Ldx*std(V95);
V2=median(V95)+Ldx*std(V95);
idx3=V95>=V1 & V95<=V2;

V95=wga.Ay;
V1=median(V95)-Ldx*std(V95);
V2=median(V95)+Ldx*std(V95);
idx4=V95>=V1 & V95<=V2 ;

V95=wga.v;
V1=median(V95(idx2))-Ldx*std(V95(idx2));
V2=median(V95(idx2))+Ldx*std(V95(idx2));
idx5=V95>=V1 & V95<=V2;

V95=wga.Bx;
V1=median(V95)-Ldx*std(V95);
V2=median(V95)+Ldx*std(V95);
idx6=V95>=V1 & V95<=V2; 

V95=wga.By;
V1=median(V95)-Ldx*std(V95);
V2=median(V95)+Ldx*std(V95);
idx7=V95>=V1 & V95<=V2;

idx=idx2 & idx3 & idx4 & idx5 & idx6 & idx7;
wga.stn=wga.stn(idx,:);
wga.st=wga.st(idx,1:3);
wga.v=wga.v(idx);
wga.dv=wga.dv(idx);
wga.a0=wga.a0(idx);
wga.da=wga.da(idx);
wga.rd=wga.rd(idx);
wga.dr=wga.dr(idx);
wga.gs=wga.gs(idx);
wga.dg=wga.dg(idx);
wga.azmo=wga.azmo(idx);
wga.nst=wga.nst(idx);
wga.envu0=wga.envu0(idx);
wga.Time=wga.Time(idx);
wga.amp=wga.amp(idx);
wga.ev=wga.ev(idx,:);
wga.Ax = wga.Ax(idx);
wga.Ay = wga.Ay(idx);
wga.Bx = wga.Bx(idx);
wga.By = wga.By(idx);
wga.weight = weight*SNRwi(idx);
wga.SNR = wga.SNR(idx);
wga.dsmax = wga.dsmax(idx);




%%
%%average around a station
%%
function bk=smoothavg(wga,dsmax,dsmax2)
st = wga.st;
[nst,~] = size(st);
bk = [];
for is = 1:nst    
    yi = st(is,1);
    xi = st(is,2);
    dst = distance(yi,xi,st(:,1),st(:,2))*111.1949;    
    idx = dst <= dsmax;
    idx2 = dst <= dsmax2;
    if find(idx2)<3
        wga.weight(is,1)=wga.weight(is,1)*0.7;
    end
    bk.stn(is) = {wga.stn(is,:)};
    bk.st(is,1:3) = st(is,1:3);
    bk.azmo(is,1) = wga.azmo(is,1);
    bk.nst(is,1) = wga.nst(is,1);
    bk.envu0(is,1) = wga.envu0(is,1);
    bk.Time(is,1) = wga.Time(is,1);
    bk.amp(is,1) = wga.amp(is,1);
    bk.evla(is,1) = wga.ev(is,1);
    bk.evlo(is,1) = wga.ev(is,2);
    bk.evdpth(is,1) = wga.ev(is,3);
    bk.weight(is,1) = wga.weight(is,1);
    bk.SNR(is,1) = wga.SNR(is,1);
    bk.v(is,1) = median(wga.v(idx));
    bk.va(is,1) = median(wga.v(idx2));
    bk.dv(is,1) = median(wga.dv(idx));
    bk.az(is,1) = median(wga.a0(idx));
    bk.da(is,1) = median(wga.da(idx));
    bk.rd(is,1) = median(wga.rd(idx));
    bk.dr(is,1) = median(wga.dr(idx));
    bk.gs(is,1) = median(wga.gs(idx));
    bk.dg(is,1) = median(wga.dg(idx));
    bk.Ax(is,1) = median(wga.Ax(idx));
    bk.Ay(is,1) = median(wga.Ay(idx));
    bk.Bx(is,1) = median(wga.Bx(idx));
    bk.By(is,1) = median(wga.By(idx));
    bk.dBdx(is,1) = median(wga.dBdx(idx));
    bk.dBdy(is,1) = median(wga.dBdy(idx));
    bk.dAdy(is,1) = median(wga.dAdy(idx));
    bk.dAdx(is,1) = median(wga.dAdx(idx));
    bk.AB2(is,1) = median(wga.AB2(idx));
    bk.vcr(is,1) = median(wga.vcr(idx));
    bk.vcra(is,1) = median(wga.vcr(idx2));    
    bk.dvcr(is,1) = median(wga.dvcr(idx));
    bk.dsmax(is,1) = wga.dsmax(is);
    bk.dsmr(is,:) = [dsmax,dsmax2];
end



function wga=VcrAmp(wga,ds,Ti)
x=wga.st(:,2);
y=wga.st(:,1);
Ax=wga.Ax;
Ay=wga.Ay;
Bx=wga.Bx;
By=wga.By;
[dBxdx,~]=ntrlGradient(x,y,Bx,ds,'natural');
[~,dBydy]=ntrlGradient(x,y,By,ds,'natural');
[dAxdx,~]=ntrlGradient(x,y,Ax,ds,'natural');
[~,dAydy]=ntrlGradient(x,y,Ay,ds,'natural');
AB2=(Ax.*Bx+Ay.*By)*2.0;
wga.dBdx=dBxdx;
wga.dBdy=dBydy;
wga.dAdy=dAydy;
wga.dAdx=dAxdx;
wga.AB2=AB2;
M=Bx.^2+By.^2-(Ax.^2+Ay.^2+dAxdx+dAydy)/(2*pi/Ti)^2;
Vcr=sqrt(1./M);
wga.vcr=Vcr;
wga.dvcr=wga.v-Vcr;

