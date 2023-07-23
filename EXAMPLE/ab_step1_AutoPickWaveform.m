%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cao et al. 2022 v2.0; Cao  Feihuang, 2023/06/01,v3.0
% befor runing this program£¬you need to check the parameters file of aa_WGM_parameters.m firstly.
% before running the program, also check all parameters below
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ab_step1_AutoPickWaveform
clear
close all;
par=aa_WGM_parameters; % load parameters
is_parfor =0; 
par.is_figure = 1;
par.is_overwrite = 1;
par.sacstr=['*' par.cmp '.sac'];
dataF = dir([par.ev_sac_rootpath par.fs '2007*']); 


fs=par.fs;
par.dataF = dataF;
wvPick_output = par.wvPick_output;
dout=[wvPick_output par.fs];
dout2=[wvPick_output par.fs 'wvpick' par.fs];
para_WVPICK='para_wvpick.mat';

if ~exist(dout,'dir')
    mkdir(dout)
end
if ~exist(dout2,'dir')
    mkdir(dout2)
end

save(para_WVPICK)
disp('Picking waveforms...')
if is_parfor==1 
    parfor kevi=1:length(dataF) %loops of events
        autoWVPICK(kevi,para_WVPICK)
    end
else
    for kevi=1:length(dataF) %loops of events
        autoWVPICK(kevi,para_WVPICK)
    end
end
disp('Pick done')
delete(para_WVPICK)


function autoWVPICK(kevi,para_WVPICK)
load(para_WVPICK);
dataF=par.dataF; %#ok<*NODEF>
sacstr=par.sacstr;
fs=par.fs;
is_overwrite=par.is_overwrite;
is_figure=par.is_figure;
periods = par.periods;%center periods
sacPth=[par.ev_sac_rootpath fs dataF(kevi).name fs];
sacFs = dir([sacPth sacstr]);
if length(sacFs)==0 %#ok<ISMT>
    return
end

fout1=[dout2 'wvpick.' dataF(kevi).name '.mat'];
fout2=[dout  dataF(kevi).name '.mat'];

if is_overwrite==0 && exist(fout1,'file')
    return
end
tplSAC=readsac([par.ev_sac_rootpath fs dataF(kevi).name fs sacFs(1).name]);
evdpeth=tplSAC.EVDP;
if length(sacFs)<10 && evdpeth>par.ev_maxDepth
    wvpick(1:length(periods)).evflag=deal(0); %#ok<*STRNU>
    save(fout1,'wvpick')
    return
else
    [wvpick(1:length(periods)).evflag]=deal(0);
end

disp(['Picking...' num2str(kevi) 'th event'])
if ~exist(fout2,'file')||is_overwrite==1
    [rsta]= readevent(sacPth,sacFs,par.stinfo); % read sac file
    save(fout2,'rsta') % save waveform with mat type
else
    load(fout2)
end
dist=[rsta.dis];
for pbi=1:length(periods) %loops of periods
    pband = [periods(pbi)*0.9 periods(pbi)*1.1];
    fwin=[1/pband(2) 1/pband(1)];
    par.wlen=max(pband)/2;
    pband=[periods(pbi)*0.9 periods(pbi)*1.1];
    par.pband=pband;
    [prdflg,vg,goodstID,kmax0,SNR] = TSNR2(rsta,fwin,par,is_figure);% find good waveforms    
    wvpick(pbi).datFlag=prdflg;
    wvpick(pbi).period=periods(pbi);
    wvpick(pbi).vg=vg;
    wvpick(pbi).dist=dist;
    wvpick(pbi).goodstID=goodstID;
    wvpick(pbi).kmax0=kmax0;
    wvpick(pbi).SNR=SNR;    
    wvpick(pbi).evflag=1;
    wvpick(pbi).dist=dist;
end
save(fout1,'wvpick')






function [st] = readevent(datapath,saclist,stinfo)
ist=0;
stnms=stinfo.stnm;
stloc=stinfo.stloc;
[st(1:length(saclist)).staID]=deal(0);
for ist =1:length(saclist)
    sacfilename = [datapath,saclist(ist).name];
    sac = readsac(sacfilename);
    sac.STEL=round(sac.STEL+(rand(1)-0.5)*20);
    sac.STLA=round(sac.STLA*20+(rand(1)-0.5))/20;
    sac.STLO=round(sac.STLO*20+(rand(1)-0.5))/20;
    if ~isempty(sac)
        % initial the event information by the first sac file
        st(ist).stn = sac.KSTNM;
        if isnan(sac.STLA)
            idx=find(ismember(stnms,sac.KSTNM),1);
            if isempty(idx)
                continue
            end            
            st(ist).dis = distance(stloc(idx,1),stloc(idx,2),sac.EVLA,sac.EVLO)*111.1949;
            st(ist).azm = azimuth(sac.EVLA,sac.EVLO,stloc(idx,1),stloc(idx,2)); 
            st(ist).bazm = azimuth(stloc(idx,1),stloc(idx,2),sac.EVLA,sac.EVLO); 
            st(ist).st=stloc(idx,:);
        else
            st(ist).st = [sac.STLA, sac.STLO sac.STEL];
            st(ist).dis = sac.DIST;
            st(ist).azm = sac.AZ;
            st(ist).bazm = sac.BAZ;            
        end
        st(ist).ev = [sac.EVLA  sac.EVLO sac.EVDP];
        st(ist).dt = round(sac.DELTA*1000)/1000;
        st(ist).dat = sac.DATA1;
        st(ist).npt = sac.NPTS;
        st(ist).cmp = sac.KCMPNM;
        st(ist).net=sac.KNETWK;
        st(ist).tphs = [sac.T0 sac.T1 sac.T2 sac.T3 sac.T4...
            sac.T5 sac.T6 sac.T7 sac.T8 sac.T9];
        
        tbeg = sac.B;
        tend = sac.E;
        t0   = [sac.NZYEAR sac.NZJDAY sac.NZHOUR sac.NZMIN sac.NZSEC];
        t0(5)=t0(5)+sac.NZMSEC/1000+tbeg;
        st(ist).t0 = t0;
        t1 = t0;
        t1(5) = t1(5)+tend-tbeg;
        st(ist).t1 = t1;
        st(ist).tbeg = tbeg;
        st(ist).dst=nan;
        st(ist).staID=ist;
        st(ist).scale = 1;

    end
end
badID=[st(ist).scale]==0;
st(badID)=[];
npts=[st.npt];
dnpt=npts-median(npts);

idx1=dnpt==1;
idx2=dnpt==-1;
idx3=dnpt==0;

st1=st(idx1);
st2=st(idx2);
st3=st(idx3);

if ~isempty(find(idx1,1))
    for i=1:length(st1)
        st1(i).dat(end)=[];
        st1(i).t1(5)=st1(i).t1(5)-st1(i).dt;
        st1(i).npt=st1(i).npt-1;
    end
end
if ~isempty(find(idx2,1))
    for i=1:length(st2)
        dat(1:st2(i).npt+1,1)=[st2(i).dat;st2(i).dat(st2(i).npt)];
        st2(i).dat=dat;
        st2(i).t1(5)=st2(i).t1(5)+st2(i).dt;
        st2(i).npt=st2(i).npt+1;
    end
end
st=[st1 st2 st3];


function  [datFlg,vg,gd1,kmax0,SNRs]= TSNR2(sts,fwin,par,is_figure)
par.vgmin=par.vglen(1);
par.vgmax=par.vglen(2);
kmax0(1:length(sts))=nan;
vg(1:length(sts))=0;
SNRs(1:length(sts))=0;
gd1(1:length(sts))=0;
datFlg=0; %#ok<*NASGU>
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
    kmaxs=findmaxima(wv.ei(i,:),5);
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

%% plot picking results
if is_figure
   plotQwave(wv,par);
    figure(11)
    drawnow
    figure(12)
    drawnow
    figure(13)
    drawnow
    figure(14)
    drawnow
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
kmax(isnan(kmax))=[];


function plotQwave(wvQ,par)
fwin=[1/par.pband(2) 1/par.pband(1)];
par.vgmin=par.vglen(1);
par.vgmax=par.vglen(2);
mixi(1:length(wvQ.ei(:,1)))=nan;
maxi(1:length(wvQ.ei(:,1)))=nan;
ddyi(1:length(wvQ.ei(:,1)))=nan;
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
    if i==2
    hold on
    end
    X=(1:length(wf2))*wvQ.dt+wvQ.begtime(i);
    Y=30*wf2/max(wf2)+wvQ.dist(i)-1;
    plot(X,Y,'k')
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
mvg = mean(wvQ.vg);
h=plot(dist(wvQ.goodsta==1),wvQ.vg(wvQ.goodsta==1),'o');
hold on
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
    for i=1:length(wvQ.ei(:,1))
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
        plot(X,Y,'k')
        hold on

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
            X=(1:length(wf))*wvQ.dt+wvQ.begtime(i);
            Y=40*wf/max(wf)+wvQ.dist(i)-1;
            plot(X,Y,'k')
            hold on

            %             X=int*wvQ.dt+wvQ.begtime(i);
            %             Y=40*wf(int)/max(wf(int))+wvQ.dist(i)-1;
            %             plot(X,Y,'color',[90 90 90]/255)
        else
            figure(14)
            wf = wvQ.pbwfi(i,:);
            mxt=round((wvQ.dist(i)/(median(wvQ.vg)-0.5)-wvQ.begtime(i))/wvQ.dt);
            mit=round((wvQ.dist(i)/(median(wvQ.vg)+0.5)-wvQ.begtime(i))/wvQ.dt);
            int=mit:mxt;
            int(int>length(wf))=[];
            int(int<1)=[];
            if isempty(int)||isnan(max(wvQ.ei(i,:)))
                continue
            end
            X=(1:length(wf))*wvQ.dt+wvQ.begtime(i);
            Y=40*wf/max(wf)+wvQ.dist(i)-1;
            plot(X,Y,'k')
            hold on
            X=int*wvQ.dt+wvQ.begtime(i);
            Y=40*wf(int)/max(wf)+wvQ.dist(i)-1;
            plot(X,Y,'color',[0 90 0]/255)
            
            %time windows with 10 wavelength
            int2=round(mean(5./fwin)/wvQ.dt);
            int2=(-int2:int2)+wvQ.kmax0(i);
            int2(int2<1)=[];
            int2(int2>length(wvQ.ei(i,:)))=[];
            if isempty(int2)
                plot(X,Y,'color',[255 90 0]/255)
                continue
            end
            X=int2*wvQ.dt+wvQ.begtime(i);
            Y=40*wvQ.ei(i,int2)/max(wvQ.ei(i,:))+wvQ.dist(i)-1;
            plot(X,Y,'color',[0 0 255]/255)
            scatter(wvQ.kmax0(i)*wvQ.dt+wvQ.begtime(i),max(Y),'red');
            text(wvQ.kmax0(i)*wvQ.dt,max(Y),num2str(wvQ.SNRs(i)));
            
        end
    end
    figure(13)
    plot(wvQ.kmax0(wvQ.goodsta==1)*wvQ.dt+wvQ.begtime(wvQ.goodsta==1),ddyi(wvQ.goodsta==1),'.g')
    title('Good waveforms')
    xlabel('time (s)')
    ylabel('epicentral distance (km)')
    xlim([min(mixi),max(maxi)])
    box on
    hold off
    
    figure(14)
    plot(wvQ.kmax0(wvQ.goodsta==0)*wvQ.dt+wvQ.begtime(wvQ.goodsta==0),ddyi(wvQ.goodsta==0),'.r')
    title('Bad waveforms')
    xlabel('time (s)')
    ylabel('epicentral distance (km)')
    xlim([min(mixi),max(maxi)])
    box on
    hold off
end
